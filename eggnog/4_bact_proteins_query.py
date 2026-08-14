#!/usr/bin/env python3
"""
bact_proteins_query.py

Pull the enzyme set of one species / strain out of the bacterial EC table, in the shape
../1_query/ec2food.py expects as --enzyme_tsv.

Input (--tsv): /data/bac2food/exports/species_enzymes.tsv — the eggNOG v6 EC table exported
from /data/bac2food/bact_ec.tsv by ../5_export/export_resources.py, with (at least):
ec_number, genus, species, strain.
It previously read 3_normalized_species.tsv, an incomplete BRENDA SPARQL scrape that has
since been deleted — that table under-reported the enzymes of every organism.

Query behavior (comma-separated list allowed):
- If a query looks like a BINOMIAL species (>=2 tokens), match that species and return ALL its strains.
- If a query looks like a STRAIN query, match ONLY that exact strain (within its species).
  Two supported strain query forms:
    1) "Genus species strain_tokens..."  (>=3 tokens; strain inferred as tokens after binomial)
    2) "strain:NBRC 3283" or "strain=NBRC 3283" (explicit strain query; matches strain column directly)

Note: names follow eggNOG's older taxonomy (Lactobacillus, not the reclassified
Lacticaseibacillus / Levilactobacillus / ...). Query the name as it appears in the table.

Output:
- --rows: full rows
- default: unique enzymes per (species, strain, ec) and includes genus/species/strain columns
Sorting: species, then strain, then ec, then recommended_name (if present)

Examples
  # species -> all strains
  ./4_bact_proteins_query.py --query "Akkermansia muciniphila"

  # strain-only (implicit: 3+ tokens)
  ./4_bact_proteins_query.py --query "Bacteroides thetaiotaomicron VPI-5482" --out btheta.tsv

  # strain-only (explicit)
  ./4_bact_proteins_query.py --query "strain:VPI-5482"

  # multiple queries
  ./4_bact_proteins_query.py --query "Akkermansia muciniphila, Escherichia coli CFT073"
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import unicodedata
from difflib import get_close_matches
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set

GENUS_ABBREV = {
    "e.": "escherichia",
    "b.": "bacillus",
    "p.": "pseudomonas",
    "s.": "staphylococcus",
    "l.": "lactobacillus",
    "c.": "clostridium",
}

_CTRL_RE = re.compile(r"[\x00-\x08\x0B\x0C\x0E-\x1F\x7F]")


def clean_field(v: str) -> str:
    if v is None:
        return ""
    v = unicodedata.normalize("NFKC", str(v))
    v = v.replace("\r\n", "\n").replace("\r", "\n").replace("\n", ". ")
    v = _CTRL_RE.sub(" ", v)
    v = re.sub(r"\s+", " ", v).strip()
    return v


def norm_name(s: str) -> str:
    if s is None:
        return ""
    s = unicodedata.normalize("NFKC", s).strip().casefold()
    s = s.replace("_", " ").replace("-", " ")
    parts = s.split()
    if len(parts) >= 2 and parts[0] in GENUS_ABBREV:
        parts[0] = GENUS_ABBREV[parts[0]]
        s = " ".join(parts)
    s = re.sub(r"[^0-9a-z\s]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def parse_query_list(q: str) -> List[str]:
    return [p.strip() for p in (q or "").split(",") if p.strip()]


def read_header(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8", newline="") as f:
        header = next(csv.reader(f, delimiter="\t"), None)
    if not header:
        raise ValueError("Input TSV has no header.")
    return header


def read_organisms(path: Path, header: List[str], genus_col: str, species_col: str,
                   strain_col: str) -> List[Dict[str, str]]:
    """Stream the table and return only its DISTINCT organisms.

    species_enzymes.tsv has ~9.6M rows but only ~3k organisms. Indexing the distinct
    organisms rather than every row is what keeps this in a few MB of memory instead of
    the gigabytes a row-level index would need.
    """
    gi, si, ti = header.index(genus_col), header.index(species_col), header.index(strain_col)
    seen: Dict[Tuple[str, str], Dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) <= max(gi, si, ti):
                continue
            key = (row[si], row[ti])
            if key not in seen:
                seen[key] = {genus_col: row[gi], species_col: row[si], strain_col: row[ti]}
    return list(seen.values())


def select_rows(path: Path, header: List[str], wanted: Set[Tuple[str, str]],
                species_col: str, strain_col: str) -> List[Dict[str, str]]:
    """Second pass: keep the rows belonging to the matched organisms."""
    si, ti = header.index(species_col), header.index(strain_col)
    out: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) <= max(si, ti):
                continue
            if (row[si], row[ti]) in wanted:
                out.append(dict(zip(header, row)))
    return out


def looks_like_explicit_strain_query(q: str) -> Optional[str]:
    """
    Returns the strain value if query is explicit 'strain:...' or 'strain=...'
    """
    m = re.match(r"^\s*strain\s*[:=]\s*(.+?)\s*$", q, flags=re.IGNORECASE)
    return m.group(1).strip() if m else None


def split_binomial_and_strain(q: str) -> Tuple[str, str]:
    """
    For implicit "Genus species strain..." form:
      returns (binomial, strain_tail)
    """
    qn = norm_name(q)
    toks = qn.split()
    if len(toks) < 3:
        return qn, ""
    binom = " ".join(toks[:2])
    strain_tail = " ".join(toks[2:])
    return binom, strain_tail


def build_indexes(rows: List[Dict[str, str]], genus_col: str, species_col: str, strain_col: str) -> Tuple[Dict[str, List[int]], Dict[Tuple[str, str], List[int]]]:
    """
    Indexes the DISTINCT organisms (not the rows) returned by read_organisms().

    Returns:
      species_idx: normalized species -> organism indices
      strain_idx:  (normalized species, normalized strain) -> organism indices
    """
    species_idx: Dict[str, List[int]] = {}
    strain_idx: Dict[Tuple[str, str], List[int]] = {}

    for i, r in enumerate(rows):
        sp = norm_name(r.get(species_col, ""))
        st = norm_name(r.get(strain_col, ""))
        if sp:
            species_idx.setdefault(sp, []).append(i)
            if st:
                strain_idx.setdefault((sp, st), []).append(i)

    return species_idx, strain_idx


def match_one_query(
    q: str,
    species_keys: List[str],
    species_idx: Dict[str, List[int]],
    strain_idx: Dict[Tuple[str, str], List[int]],
    allow_fuzzy: bool,
) -> Tuple[List[int], List[str]]:
    """
    Returns (row_indices, suggestions)

    Rules:
    - explicit strain: match strain across any species if unique, else require binomial form
    - implicit strain (>=3 tokens): match that strain within that binomial only
    - species (>=2 tokens): match species and return all strains (i.e., all rows for that species)
    - genus-only (1 token): match all species under that genus via startswith()
    """
    q = q.strip()
    if not q:
        return [], []

    explicit_strain = looks_like_explicit_strain_query(q)
    if explicit_strain:
        stn = norm_name(explicit_strain)
        if not stn:
            return [], []
        # find all (species, strain) matches for this strain token
        hits = [idxs for (sp, st), idxs in strain_idx.items() if st == stn]
        if len(hits) == 1:
            return hits[0], []
        elif len(hits) > 1:
            # ambiguous strain across multiple species
            sugg = [f"{sp} {st}" for (sp, st) in strain_idx.keys() if st == stn][:10]
            return [], [f"Ambiguous strain '{explicit_strain}' across species. Use 'Genus species {explicit_strain}'. Examples:"] + sugg
        else:
            return [], get_close_matches(stn, sorted({st for (_, st) in strain_idx.keys()}), n=10, cutoff=0.85)

    qn = norm_name(q)
    toks = qn.split()

    if len(toks) >= 3:
        # implicit strain: "genus species strain..."
        binom, st = split_binomial_and_strain(q)
        if not st:
            return [], []
        key = (binom, st)
        if key in strain_idx:
            return strain_idx[key], []
        # allow exact species but fuzzy strain (optional)
        if allow_fuzzy:
            candidate_strains = [st2 for (sp2, st2) in strain_idx.keys() if sp2 == binom]
            sugg = get_close_matches(st, candidate_strains, n=10, cutoff=0.85)
            if sugg:
                return [], [f"No exact strain match for '{q}'. Close strains under {binom}:"] + sugg
        return [], []

    if len(toks) >= 2:
        # species
        binom = " ".join(toks[:2])
        if binom in species_idx:
            return species_idx[binom], []
        # allow prefix if species column might include extra tokens
        pref = [k for k in species_keys if k == binom or k.startswith(binom + " ")]
        if pref:
            out: List[int] = []
            for k in pref:
                out.extend(species_idx.get(k, []))
            return out, []
        if allow_fuzzy:
            return [], get_close_matches(binom, species_keys, n=10, cutoff=0.85)
        return [], []

    # genus-only
    if len(toks) == 1:
        g = toks[0]
        pref = [k for k in species_keys if k.startswith(g + " ")]
        if pref:
            out: List[int] = []
            for k in pref:
                out.extend(species_idx.get(k, []))
            return out, []
        if allow_fuzzy:
            return [], get_close_matches(g, sorted({k.split()[0] for k in species_keys if k.split()}), n=10, cutoff=0.85)
        return [], []

    return [], []


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="/data/bac2food/exports/species_enzymes.tsv",
                    help="Input TSV: the eggNOG EC table from 5_export (default: %(default)s)")
    ap.add_argument(
        "--query",
        required=True,
        help='Comma-separated queries. Species: "Genus species". Strain: "Genus species STRAIN" or "strain:STRAIN".',
    )
    ap.add_argument("--genus-col", default="genus")
    ap.add_argument("--species-col", default="species")
    ap.add_argument("--strain-col", default="strain")
    ap.add_argument("--ec-col", default="ec_number")
    ap.add_argument("--name-col", default="recommended_name", help="Optional enzyme name column if present")
    ap.add_argument("--allow-fuzzy", action="store_true")
    ap.add_argument("--rows", action="store_true", help="Output original rows")
    ap.add_argument("--out", default="-")
    args = ap.parse_args(argv)

    path = Path(args.tsv)
    fieldnames = read_header(path)
    for col in (args.genus_col, args.species_col, args.strain_col, args.ec_col):
        if col not in fieldnames:
            print(f"ERROR: missing column '{col}'. Columns: {', '.join(fieldnames)}", file=sys.stderr)
            return 2

    organisms = read_organisms(path, fieldnames, args.genus_col, args.species_col, args.strain_col)
    species_idx, strain_idx = build_indexes(organisms, args.genus_col, args.species_col, args.strain_col)
    species_keys = sorted(species_idx.keys())

    queries = parse_query_list(args.query)
    if not queries:
        print("ERROR: --query parsed to empty list", file=sys.stderr)
        return 2

    matched_indices: Set[int] = set()
    had_miss = False

    for q in queries:
        idxs, sugg = match_one_query(q, species_keys, species_idx, strain_idx, allow_fuzzy=args.allow_fuzzy)
        if idxs:
            matched_indices.update(idxs)
        else:
            had_miss = True
            qn = norm_name(q)
            print(f"NOTE: no match for {q!r} (normalized: {qn!r})", file=sys.stderr)
            if sugg:
                for s in sugg[:12]:
                    print(f"  {s}", file=sys.stderr)

    if not matched_indices:
        print("ERROR: no matches.", file=sys.stderr)
        return 1

    wanted = {(organisms[i][args.species_col], organisms[i][args.strain_col])
              for i in matched_indices}
    out_rows = select_rows(path, fieldnames, wanted, args.species_col, args.strain_col)

    # Sort: species, strain, ec, name
    def sort_key(r: Dict[str, str]):
        sp = norm_name(r.get(args.species_col, ""))
        st = norm_name(r.get(args.strain_col, ""))
        ec = clean_field(r.get(args.ec_col, ""))
        nm = clean_field(r.get(args.name_col, "")) if args.name_col in fieldnames else ""
        return (sp, st, ec, nm.casefold())

    out_rows.sort(key=sort_key)

    out_f = sys.stdout if args.out == "-" or not args.out.strip() else open(args.out, "w", encoding="utf-8", newline="")
    close_after = out_f is not sys.stdout

    try:
        if args.rows:
            w = csv.DictWriter(out_f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in out_rows:
                w.writerow({fn: clean_field(r.get(fn, "")) for fn in fieldnames})
        else:
            has_name = args.name_col in fieldnames
            w = csv.writer(out_f, delimiter="\t", lineterminator="\n")
            if has_name:
                w.writerow([args.genus_col, args.species_col, args.strain_col, args.ec_col, args.name_col])
            else:
                w.writerow([args.genus_col, args.species_col, args.strain_col, args.ec_col])

            seen = set()
            for r in out_rows:
                g = clean_field(r.get(args.genus_col, ""))
                sp = clean_field(r.get(args.species_col, ""))
                st = clean_field(r.get(args.strain_col, ""))
                ec = clean_field(r.get(args.ec_col, ""))
                if not sp or not ec:
                    continue
                if has_name:
                    nm = clean_field(r.get(args.name_col, ""))
                    key = (sp, st, ec, nm)
                    if key in seen:
                        continue
                    seen.add(key)
                    w.writerow([g, sp, st, ec, nm])
                else:
                    key = (sp, st, ec)
                    if key in seen:
                        continue
                    seen.add(key)
                    w.writerow([g, sp, st, ec])
    finally:
        if close_after:
            out_f.close()

    if had_miss:
        return 0  # still successful; misses were warned on stderr
    return 0


if __name__ == "__main__":
    raise SystemExit(main())