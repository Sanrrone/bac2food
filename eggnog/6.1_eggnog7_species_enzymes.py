#!/usr/bin/env python3
"""6.1_eggnog7_species_enzymes.py — rebuild bacterium -> EC from eggNOG v7.

Replaces the v6 route (a direct EC annotation, parsed into /data/bac2food/bact_ec.tsv
and flattened by ../5_export/export_resources.py). eggNOG v7's bulk annotation file
carries KEGG Orthology ids instead of EC numbers, so the chain gains one hop:

    OG member protein  ->  tax_id          (e7.og_info_kegg_go.tsv, col 6)
    OG                 ->  KEGG KO         (e7.og_info_kegg_go.tsv, col 7)
    KO                 ->  EC              (KEGG list/ko, via 6.0_kegg_ko_to_ec.py)
    tax_id             ->  organism name   (e7.taxid_info.tsv.gz)

The provenance is therefore bacterium -> KO -> EC, not the v6 bacterium -> EC. The
extra hop is lossy in a way worth stating: only ~39% of KOs are enzymes at all, and
KO -> EC is many-to-many, so the EC set is not a superset of v6's.

Inputs (all on /data; nothing here needs the 12 GB protein FASTA):
    e7.og_info_kegg_go.tsv   2.7 GB, 3,182,553 rows, 9 tab-separated columns
        col 6 = members "taxid.protein,taxid.protein,..."
        col 7 = "K01046|30.00;K11144|8.33"  (KO | % of OG members carrying it)
    e7.taxid_info.tsv.gz     1.2 MB  New_Taxid, Old_Taxid, Sci_Name, Rank, lineages
    kegg_ko_ec.tsv           from 6.0_kegg_ko_to_ec.py

Two details that decide correctness:
  * Member protein prefixes use eggNOG's **Old_Taxid**, so that is the join key; the
    emitted tax_id is the **New_Taxid** (current NCBI), which is the point of moving
    to v7 -- it retires the stale-taxonomy caveat the v6 export carried.
  * A protein sits in several nested OGs, so (tax_id, EC) repeats; pairs are
    deduplicated as integer codes to keep peak memory near 1 GB.

Only 4-level ECs are emitted: partial ECs ("1.1.1.-") cannot join the EC -> substrate
digest downstream, so they would be dead weight.

On --min_consensus: col 7 records the %% of OG members carrying the KO, and an OG's KO
is transferred to all of its members here (eggNOG-mapper's own annotation-transfer
rule). Filtering on that %% looks prudent but measurably is not -- coverage of the ECs
that actually reach an FDC nutrient falls off a cliff, because broad OGs legitimately
mix annotated and unannotated members:

    threshold   ECs    nutrient-reaching EC coverage
    v6 (direct) 4,291            75.9%
    0 (default) 4,819            82.1%     <- best
    >= 20       2,563            51.3%
    >= 50       1,879            39.8%
    >= 90         887            21.1%

So the default is 0; the flag exists to make that trade-off reproducible, not to be
raised casually.

Usage:
    # prototype on the first 300k OG rows, no output file
    python 6.1_eggnog7_species_enzymes.py --sample 300000 --stats_only

    # full build
    python 6.1_eggnog7_species_enzymes.py --out /data/bac2food/exports/species_enzymes.tsv
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA = Path("/data/bac2food")

sys.path.insert(0, str(REPO / "5_export"))
from export_resources import split_organism  # noqa: E402  (same parser as the v6 export)

BACTERIA_TAXID = "2"
FLUSH_EVERY = 40_000_000          # pair codes buffered before a dedup pass


# ==============================================================================
# inputs
# ==============================================================================
def load_ko_ec(path: Path) -> dict[str, list[str]]:
    """KO -> list of full 4-level ECs. Partial ECs are dropped (see module docstring)."""
    ko2ec: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as fh:
        header = fh.readline()
        if not header.startswith("ko\t"):
            raise SystemExit(f"ERROR: {path} is not the 6.0_kegg_ko_to_ec.py output")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[2] == "yes":
                continue
            ko2ec.setdefault(p[0], []).append(p[1])
    if not ko2ec:
        raise SystemExit(f"ERROR: no full-EC pairs in {path}")
    return ko2ec


def load_taxa(path: Path) -> tuple[dict[str, int], list[tuple[str, str]]]:
    """Return (Old_Taxid -> canonical index, [(New_Taxid, Sci_Name), ...]) for bacteria.

    Bacteria are identified by NCBI taxid 2 appearing in Taxid_Lineage, which is
    exact -- unlike matching the word "Bacteria" in the name lineage.

    Member proteins are prefixed with Old_Taxid, but the emitted tax_id is New_Taxid,
    and NCBI has *merged* some taxa: 41 bacterial New_Taxids absorb 49 extra Old_Taxids.
    Indexing on the old id would therefore emit the same (tax_id, EC) fact twice, so the
    canonical index is keyed on New_Taxid and several old ids may point at one entry.
    Where merged records disagree on Sci_Name, the record whose Old_Taxid already equals
    New_Taxid wins, so the name is the current one rather than a retired synonym.
    """
    canon_idx: dict[str, int] = {}
    canon: list[tuple[str, str]] = []
    old2canon: dict[str, int] = {}
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            new, old, sci = p[0], p[1], p[2]
            if BACTERIA_TAXID not in p[5].split(","):
                continue
            i = canon_idx.get(new)
            if i is None:
                i = canon_idx[new] = len(canon)
                canon.append((new, sci))
            elif old == new:                      # authoritative record for a merged taxon
                canon[i] = (new, sci)
            old2canon[old] = i
    if not canon:
        raise SystemExit(f"ERROR: no bacterial taxa parsed from {path}")
    return old2canon, canon


# ==============================================================================
# scan
# ==============================================================================
def scan_og_info(path: Path, ko2ec: dict[str, list[str]], old2canon: dict[str, int],
                 min_consensus: float, sample: int | None) -> tuple[np.ndarray, list[str], dict]:
    """Stream the OG file and return unique (tax, ko) pairs as packed integer codes.

    Non-enzymatic KOs and non-bacterial members are discarded before buffering, which
    is what keeps this tractable: it removes ~61% of KOs and every eukaryote/archaeon
    up front rather than after the join.
    """
    ko_list = sorted(ko2ec)
    ko_idx = {k: i for i, k in enumerate(ko_list)}
    tax_idx = old2canon
    n_ko = len(ko_list)

    buf = np.empty(FLUSH_EVERY, dtype=np.int64)
    fill = 0
    uniq = np.empty(0, dtype=np.int64)
    st = {"rows": 0, "rows_ko": 0, "ko_tokens": 0, "ko_kept": 0,
          "members": 0, "members_bact": 0}

    def flush(u: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.unique(np.concatenate([u, b]) if u.size else b)

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            st["rows"] += 1
            if sample and st["rows"] > sample:
                st["rows"] -= 1
                break
            if st["rows"] % 200_000 == 0:
                print(f"    {st['rows']:,} OG rows | {uniq.size + fill:,} pairs ...",
                      end="\r", flush=True)

            p = line.rstrip("\n").split("\t")
            if len(p) < 7 or not p[6]:
                continue
            st["rows_ko"] += 1

            # ---- KOs of this OG, filtered to enzymes above the consensus cut ----
            kos = []
            for tok in p[6].split(";"):
                if not tok:
                    continue
                st["ko_tokens"] += 1
                ko, _, cons = tok.partition("|")
                if ko not in ko_idx:
                    continue
                if min_consensus > 0.0:
                    try:
                        if float(cons) < min_consensus:
                            continue
                    except ValueError:
                        continue
                kos.append(ko_idx[ko])
            if not kos:
                continue
            st["ko_kept"] += len(kos)

            # ---- bacterial members of this OG ----
            members = p[5]
            if not members:
                continue
            tids = []
            for m in members.split(","):
                st["members"] += 1
                t = m.split(".", 1)[0]
                i = tax_idx.get(t)
                if i is not None:
                    tids.append(i)
            if not tids:
                continue
            st["members_bact"] += len(tids)

            tids = np.unique(np.asarray(tids, dtype=np.int64))
            codes = (tids[:, None] * n_ko + np.asarray(kos, dtype=np.int64)[None, :]).ravel()

            if fill + codes.size > buf.size:
                uniq = flush(uniq, buf[:fill])
                fill = 0
                if codes.size > buf.size:          # single OG larger than the buffer
                    uniq = flush(uniq, np.unique(codes))
                    continue
            buf[fill:fill + codes.size] = codes
            fill += codes.size

    if fill:
        uniq = flush(uniq, buf[:fill])
    print(" " * 70, end="\r")
    return uniq, ko_list, st


# ==============================================================================
# output
# ==============================================================================
def _ec_key(ec: str) -> tuple:
    return tuple(int(x) for x in ec.split("."))


def write_species_enzymes(codes: np.ndarray, ko_list: list[str], canon: list[tuple[str, str]],
                          ko2ec: dict[str, list[str]], out: Path) -> dict:
    """Expand (tax, ko) codes through KO -> EC and write the sorted TSV.

    Sorting is done on integer ranks (organism name, then EC by numeric level) and the
    strings are looked up only while writing, so a ~10M-row table never exists as text
    in memory.
    """
    n_ko = len(ko_list)
    tax_of = (codes // n_ko).astype(np.int32)
    ko_of = (codes % n_ko).astype(np.int32)

    # KO -> EC is many-to-many; expand by repeating each pair once per EC.
    ec_names: list[str] = sorted({e for ecs in ko2ec.values() for e in ecs}, key=_ec_key)
    ec_idx = {e: i for i, e in enumerate(ec_names)}
    ko_ec_idx = [np.asarray([ec_idx[e] for e in ko2ec[k]], dtype=np.int32) for k in ko_list]
    counts = np.asarray([a.size for a in ko_ec_idx], dtype=np.int32)

    rep = counts[ko_of]
    tax_rep = np.repeat(tax_of, rep)
    ec_rep = np.concatenate([ko_ec_idx[i] for i in ko_of]) if ko_of.size else np.empty(0, np.int32)

    pair = np.unique(tax_rep.astype(np.int64) * len(ec_names) + ec_rep.astype(np.int64))
    tax_f = (pair // len(ec_names)).astype(np.int32)
    ec_f = (pair % len(ec_names)).astype(np.int32)

    # per-organism strings, parsed once each (~17k) rather than per row
    new_tax = [c[0] for c in canon]
    org_of = [c[1] for c in canon]
    parsed = [split_organism(o) for o in org_of]
    # rank by (name, tax_id) so organisms sharing a name still sort deterministically
    org_rank = np.empty(len(canon), dtype=np.int64)
    for r, i in enumerate(sorted(range(len(org_of)), key=lambda i: (org_of[i], new_tax[i]))):
        org_rank[i] = r

    order = np.lexsort((ec_f, org_rank[tax_f]))   # ec_names is already in numeric EC order

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        fh.write("tax_id\tgenus\tspecies\tstrain\torganism\tec_number\n")
        rows = []
        for j in order:
            t = tax_f[j]
            g, sp, strn = parsed[t]
            rows.append(f"{new_tax[t]}\t{g}\t{sp}\t{strn}\t{org_of[t]}\t{ec_names[ec_f[j]]}\n")
            if len(rows) >= 500_000:
                fh.writelines(rows)
                rows = []
        fh.writelines(rows)

    return {"rows": int(pair.size),
            "organisms": int(np.unique(tax_f).size),
            "ecs": int(np.unique(ec_f).size),
            "species": len({parsed[t][1] for t in np.unique(tax_f)})}


def main() -> int:
    ap = argparse.ArgumentParser(description="Rebuild bacterium -> EC from eggNOG v7.")
    ap.add_argument("--og_info", default=str(DATA / "e7.og_info_kegg_go.tsv"))
    ap.add_argument("--taxid_info", default=str(DATA / "e7.taxid_info.tsv.gz"))
    ap.add_argument("--ko_ec", default=str(DATA / "kegg_ko_ec.tsv"))
    ap.add_argument("--out", default=str(DATA / "exports/species_enzymes.tsv"))
    ap.add_argument("--min_consensus", type=float, default=0.0,
                    help="Drop KO calls below this %% of OG members. Default 0 (keep all) "
                         "is the measured optimum -- see the docstring; raising it costs "
                         "far more coverage than it buys.")
    ap.add_argument("--sample", type=int, default=None,
                    help="Read only the first N OG rows (prototype mode)")
    ap.add_argument("--stats_only", action="store_true", help="Report, do not write the TSV")
    args = ap.parse_args()

    ko2ec = load_ko_ec(Path(args.ko_ec))
    old2canon, canon = load_taxa(Path(args.taxid_info))
    print(f"[*] bridge: {sum(len(v) for v in ko2ec.values()):,} (KO, EC) pairs "
          f"over {len(ko2ec):,} enzymatic KOs", flush=True)
    print(f"[*] taxa  : {len(canon):,} bacterial tax_ids in eggNOG v7 "
          f"({len(old2canon):,} eggNOG ids after NCBI merges)", flush=True)

    print(f"[*] scanning {args.og_info}"
          + (f" (first {args.sample:,} rows)" if args.sample else "") + " ...", flush=True)
    codes, ko_list, st = scan_og_info(
        Path(args.og_info), ko2ec, old2canon, args.min_consensus, args.sample)

    print(f"    OG rows read        : {st['rows']:,} ({st['rows_ko']:,} with a KO)", flush=True)
    print(f"    KO tokens           : {st['ko_tokens']:,} -> {st['ko_kept']:,} enzymatic & kept", flush=True)
    print(f"    member proteins     : {st['members']:,} -> {st['members_bact']:,} bacterial", flush=True)
    print(f"    unique (tax_id, KO) : {codes.size:,}", flush=True)
    if codes.size == 0:
        raise SystemExit("ERROR: no (tax_id, KO) pairs survived; check inputs/threshold")

    out = Path(args.out)
    if args.stats_only:
        out = Path("/dev/null")
    res = write_species_enzymes(codes, ko_list, canon, ko2ec, out)

    print(f"[*] {'(stats only) ' if args.stats_only else str(out) + ': '}"
          f"{res['rows']:,} rows | {res['organisms']:,} organisms "
          f"({res['species']:,} species) | {res['ecs']:,} EC numbers", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
