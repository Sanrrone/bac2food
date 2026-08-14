#!/usr/bin/env python3
"""
1.0_eggnog_ec_substrates_parser.py

Parses ec_species_substrate.tsv (from ec_species_substrate.tar.xz, eggNOG v6) into the
EC -> substrate digest. Only the EC and substrate columns are used downstream; the
bacteria -> EC map comes from /data/bac2food/bact_ec.tsv instead.

Input TSV columns:
  ?ec, ?ec_label, ?species, ?substrates

Output TSV columns:
  ec  species  substrate

- Extracts EC number from URI -> e.g. 1.1.1.1
- Emits ONE ROW PER SUBSTRATE (split on ';')
- Optionally de-duplicates substrates per (ec, species)
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from typing import Iterable, List, Tuple


EC_RE = re.compile(r"(?:/ec/|ec/)(\d+(?:\.\d+){3})\b")  # 1.1.1.1


def extract_ec(s: str) -> str:
    s = (s or "").strip()
    m = EC_RE.search(s)
    if m:
        return m.group(1)
    if re.fullmatch(r"\d+(?:\.\d+){3}", s):
        return s
    return ""


def split_semicolon_list(s: str) -> List[str]:
    parts = [p.strip() for p in (s or "").split(";")]
    return [p for p in parts if p]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True, help="Input TSV")
    ap.add_argument("--out", dest="out", default="-", help="Output TSV (default: stdout)")
    ap.add_argument(
        "--dedup",
        action="store_true",
        help="Deduplicate substrates per (ec,species) across the file",
    )
    args = ap.parse_args()

    fin = open(args.inp, "r", encoding="utf-8", newline="")
    fout = sys.stdout if args.out == "-" else open(args.out, "w", encoding="utf-8", newline="")

    seen: set[Tuple[str, str, str]] = set()

    try:
        r = csv.DictReader(fin, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: Input appears empty or not TSV.")

        ec_col = "?ec" if "?ec" in r.fieldnames else "ec"
        sp_col = "?species" if "?species" in r.fieldnames else "species"
        sub_col = "?substrates" if "?substrates" in r.fieldnames else "substrates"

        w = csv.DictWriter(fout, fieldnames=["ec", "species", "substrate"], delimiter="\t", lineterminator="\n")
        w.writeheader()

        for row in r:
            ec = extract_ec(row.get(ec_col, ""))
            species = (row.get(sp_col, "") or "").strip()
            subs = split_semicolon_list(row.get(sub_col, "") or "")

            for sub in subs:
                key = (ec, species, sub)
                if args.dedup:
                    if key in seen:
                        continue
                    seen.add(key)
                w.writerow({"ec": ec, "species": species, "substrate": sub})

    finally:
        fin.close()
        if fout is not sys.stdout:
            fout.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

