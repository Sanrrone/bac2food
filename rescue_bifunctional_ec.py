#!/usr/bin/env python3
"""rescue_bifunctional_ec.py — restore the second activity of multi-activity CAZymes.

A gene annotator gives a locus at most one EC number. A protein named
"Cellulase/esterase CelE" therefore enters the pipeline as a cellulase and its esterase half
is lost, and a "Xylosidase/arabinosidase" arrives as one of the two. Because the enzyme layer
of this resource is keyed on EC, a lost activity is a lost substrate and so a lost food link.
The activities affected are concentrated in fibre degradation — arabinoxylan backbone and side
chains, ferulate cross-links, starch debranching — which is where the loss matters most, since
those are the substrates the host does not absorb before the colon.

This reads a curated product-name -> EC table (bifunctional_ec.tsv, see its readme for the
curation rules and the exclusions) and adds the missing activities to an existing EC panel.

Two rules keep it conservative:

  ADDITIVE ONLY   an EC already on a locus is never removed or overwritten. Where the orthology
                  annotation and the product name disagree about which member of a family the
                  protein is, the annotation wins; the rescue only supplies what is absent.
  SPECIES-BOUND   a locus binning never assigned to a species is skipped, because there is no
                  bacterium to attribute the activity to.

Every added EC is written to --report with ec_source=product_name, so the share of any
downstream result that rests on name-based inference can be recovered exactly.

Input panel format is the one bac2food_predict.py reads with --mag_tsv: headerless,
<sample> <species> <locus_tag> <comma-separated ECs>. Output is the same shape, so a rescued
panel is a drop-in replacement.

Usage:
    python rescue_bifunctional_ec.py \
        --annot_dir gene_annot --panel_dir gene_annot \
        --out_dir gene_annot_rescued --report bifunctional_added.tsv

    # to also rescue loci that carry no EC at all, supply the binning map:
    python rescue_bifunctional_ec.py ... --locus2species locus_species.tsv

Without --locus2species only loci already present in the panel can be rescued, because the
species of a locus that has no EC row is unknown. That is the larger share of the loss.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

WS = re.compile(r"\s+")
EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def norm_product(s: str) -> str:
    return WS.sub(" ", s.strip()).casefold()


def load_table(path: Path) -> dict[str, list[str]]:
    """product (normalized) -> ECs the name asserts."""
    out: dict[str, list[str]] = {}
    with path.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            prod = norm_product(r["product"])
            ecs = [e.strip() for e in r["ec_numbers"].split(";") if e.strip()]
            bad = [e for e in ecs if not EC_RE.match(e)]
            if bad:
                sys.exit(f"[!] {path}: malformed EC {bad} on product {r['product']!r}")
            if not ecs:
                sys.exit(f"[!] {path}: no EC numbers on product {r['product']!r}")
            out[prod] = ecs
    if not out:
        sys.exit(f"[!] {path} holds no usable rows")
    return out


def load_annot(path: Path, wanted: dict[str, list[str]]) -> dict[str, str]:
    """locus_tag -> normalized product, for products the table covers.

    Accepts any TSV with a header naming a locus column and a product column, which is what
    Prokka and the common alternatives emit.
    """
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (rdr.fieldnames or [])}
        lc = cols.get("locus_tag") or cols.get("locus") or cols.get("gene_id")
        pc = cols.get("product") or cols.get("description") or cols.get("annotation")
        if not lc or not pc:
            sys.exit(f"[!] {path}: need a locus and a product column, got {rdr.fieldnames}")
        out = {}
        for r in rdr:
            p = norm_product(r.get(pc) or "")
            if p in wanted:
                out[(r.get(lc) or "").strip()] = p
        return out


def load_panel(path: Path) -> tuple[list[list[str]], dict[str, int]]:
    """Panel rows, plus locus -> row index."""
    rows, idx = [], {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            idx[p[2]] = len(rows)
            rows.append(p)
    return rows, idx


def load_locus2species(path: Path) -> dict[tuple[str | None, str], str]:
    """(sample|None, locus) -> species. Accepts 2 columns (locus, species), applied to every
    sample, or 3 columns (sample, locus, species)."""
    out: dict[tuple[str | None, str], str] = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                out[(p[0].strip(), p[1].strip())] = p[2].strip()
            elif len(p) == 2:
                out[(None, p[0].strip())] = p[1].strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annot_dir", default="gene_annot",
                    help="directory of <sample>.tsv gene-annotation tables")
    ap.add_argument("--panel_dir", default="gene_annot",
                    help="directory of <sample>_ec.tsv EC panels")
    ap.add_argument("--out_dir", required=True, help="where rescued panels are written")
    ap.add_argument("--table", default="bifunctional_ec.tsv")
    ap.add_argument("--locus2species", default=None,
                    help="TSV mapping loci to species, so that loci carrying no EC at all can "
                         "also be rescued. Without it only loci already in the panel are eligible.")
    ap.add_argument("--report", default=None,
                    help="provenance TSV of every added EC (default <out_dir>/bifunctional_added.tsv)")
    args = ap.parse_args()

    table = load_table(Path(args.table))
    print(f"[*] {len(table)} multi-activity product names, "
          f"{sum(len(v) for v in table.values())} activities", flush=True)

    l2s = load_locus2species(Path(args.locus2species)) if args.locus2species else {}
    if args.locus2species:
        print(f"[*] {len(l2s):,} locus->species assignments", flush=True)
    else:
        print("[*] no --locus2species: loci with no EC row cannot be attributed and are skipped",
              flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = Path(args.report) if args.report else out_dir / "bifunctional_added.tsv"

    annot_dir, panel_dir = Path(args.annot_dir), Path(args.panel_dir)
    samples = sorted(p.stem for p in annot_dir.glob("*.tsv") if not p.name.endswith("_ec.tsv"))
    if not samples:
        sys.exit(f"[!] no annotation tables in {annot_dir}")

    tally = Counter()
    per_product = Counter()
    added_rows: list[list[str]] = []

    for s in samples:
        panel_path = panel_dir / f"{s}_ec.tsv"
        if not panel_path.exists():
            print(f"    {s}: no EC panel, skipped", flush=True)
            continue
        loci = load_annot(annot_dir / f"{s}.tsv", table)
        rows, idx = load_panel(panel_path)
        sample_id = rows[0][0] if rows else s
        new_rows: list[list[str]] = []

        for locus, prod in sorted(loci.items()):
            tally["loci_seen"] += 1
            asserted = table[prod]
            if locus in idx:
                row = rows[idx[locus]]
                have = [e for e in row[3].split(",") if e]
                missing = [e for e in asserted if e not in have]
                if not missing:
                    tally["already_complete"] += 1
                    continue
                row[3] = ",".join(sorted(set(have) | set(missing)))
                tally["partial_rescued"] += 1
                species = row[1]
            else:
                species = l2s.get((sample_id, locus)) or l2s.get((None, locus))
                if not species:
                    tally["unattributable"] += 1
                    continue
                missing = list(asserted)
                new_rows.append([sample_id, species, locus, ",".join(sorted(missing))])
                tally["noec_rescued"] += 1
            tally["ec_added"] += len(missing)
            per_product[prod] += len(missing)
            for ec in missing:
                added_rows.append([sample_id, species, locus, ec, prod, "product_name"])

        rows.extend(new_rows)
        with (out_dir / f"{s}_ec.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerows(rows)

    with report.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "species", "locus_tag", "ec_number", "product", "ec_source"])
        w.writerows(added_rows)

    print(f"\n[*] loci matching the table            {tally['loci_seen']:>8,}")
    print(f"    already carried every activity     {tally['already_complete']:>8,}")
    print(f"    had some EC, missing one added     {tally['partial_rescued']:>8,}")
    print(f"    had no EC, both activities added   {tally['noec_rescued']:>8,}")
    print(f"    no species, skipped                {tally['unattributable']:>8,}")
    print(f"[*] EC assignments added               {tally['ec_added']:>8,}")
    for prod, n in per_product.most_common():
        print(f"      {n:>6}  {prod}")
    print(f"[*] Wrote {out_dir}/ and {report}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
