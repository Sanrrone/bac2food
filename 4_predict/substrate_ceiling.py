#!/usr/bin/env python3
"""substrate_ceiling.py — is MRR 0.22 a scoring failure or a data ceiling?

Six nutrient-selection variants all scored 0.12-0.22 on the graded biology gate, including
one with every rarity term switched off. When the answer does not move, the parameter being
swept is not the binding constraint. This asks the prior question instead: for each panel
species, does the resource CONTAIN a nutrient that (a) matches its documented substrate,
(b) is EC-linked to that organism, and (c) is measured in enough foods to lift a ranking?

A species whose substrate fails any of the three cannot be scored correctly by ANY formula.
That fraction is the ceiling, and the gate cannot be read without it.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_biology_panel import PANEL, NEGATIVE  # noqa: E402  same keywords, no fork

HERE = Path(__file__).resolve().parent
PANEL_EC = Path("/data/bac2food/bio_panel_FINAL/panel_ec.tsv")
N2EC = (HERE / ".." / "0_building" / "3_nutrient_to_ec.tsv").resolve()
NUTRIENT_CSV = Path("/data/bac2food/nutrient.csv")
IDX = Path("/data/bac2food/index_modeled")

# A nutrient measured in only a handful of foods cannot move a top-10 ranking unless those
# few foods win outright. Reported, not enforced — the point is to show the distribution.
THIN = 50


def main() -> int:
    nm = {}
    with NUTRIENT_CSV.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            try:
                nm[int(r.get("id") or r.get("nutrient_id"))] = r.get("name") or ""
            except (TypeError, ValueError):
                continue

    ec2n = defaultdict(set)
    with N2EC.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                ec2n[r["ec_number"]].add(int(r["nutrient_id"]))
            except (KeyError, TypeError, ValueError):
                continue

    ec_by_lbl = defaultdict(set)
    with PANEL_EC.open(encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                ec_by_lbl[p[1]].add(p[3])

    ndf = pd.read_parquet(IDX / "nutrient_df.parquet").set_index("nutrient_id")["df_foods"].to_dict()
    modeled = set(int(x) for x in ndf)

    by_species = {}
    for lbl, ecs in ec_by_lbl.items():
        sp = lbl.split("_", 1)[1].replace("_", " ")
        by_species[sp] = {n for e in ecs for n in ec2n.get(e, ()) if n in modeled}

    entries = {**PANEL, **NEGATIVE}
    print("=" * 104)
    print("Substrate ceiling — can the resource express what each organism is documented for?")
    print("=" * 104)
    print(f"  {'species':30} {'match':>5} {'linked':>6} {'best ndf':>9}  nutrient carrying the substrate")
    print("  " + "-" * 100)

    n_none, n_thin, n_ok = 0, 0, 0
    for sp, (activity, kws) in entries.items():
        targs = by_species.get(sp)
        if targs is None:
            continue
        rx = re.compile("|".join(re.escape(k) for k in kws), re.I)
        # every nutrient in the catalogue whose NAME matches the documented substrate
        matching = {n for n, name in nm.items() if rx.search(name)}
        linked = sorted(matching & targs, key=lambda n: -int(ndf.get(n, 0)))
        best = linked[0] if linked else None
        bn = int(ndf.get(best, 0)) if best else 0
        if not linked:
            verdict, n_none = "NO LINK", n_none + 1
        elif bn < THIN:
            verdict, n_thin = f"thin", n_thin + 1
        else:
            verdict, n_ok = "ok", n_ok + 1
        print(f"  {sp:30} {len(matching):5} {len(linked):6} {bn:9}  "
              f"[{verdict:7}] {nm.get(best,'—')[:40] if best else '— nothing EC-linked'}")

    tot = n_none + n_thin + n_ok
    print("  " + "-" * 100)
    print(f"  {n_ok}/{tot} species have a documented substrate that is EC-linked AND measured in "
          f">={THIN} foods")
    print(f"  {n_thin}/{tot} linked but measured in <{THIN} foods (present, but too sparse to rank)")
    print(f"  {n_none}/{tot} have NO EC-linked nutrient matching their substrate at all")
    print("\n  The last two groups are unreachable by any scoring formula. They bound the gate.")
    print("=" * 104)
    return 0


if __name__ == "__main__":
    sys.exit(main())
