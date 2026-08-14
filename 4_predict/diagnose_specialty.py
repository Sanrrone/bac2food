#!/usr/bin/env python3
"""diagnose_specialty.py — why four xylan specialists all rank strawberry first.

The biology panel scores MRR 0.216: thirteen of seventeen well-characterized organisms
never surface their documented substrate, and Bacteroides thetaiotaomicron, ovatus,
uniformis and fragilis return the SAME single food. This reproduces the nutrient
selection that produces that, outside the predictor, so the cause can be read directly.

It recomputes the Phase-9 specialty allowlist exactly as bac2food_predict.py builds it
(spec = log(B_total/bf[n]) + food_w*log(F_total/ff[n]), keep top-K per bacterium) and
prints, per species: what the allowlist contains, and where the substrates the organism
is actually documented for sit in that ranking.
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PANEL_EC = Path("/data/bac2food/bio_panel_FINAL/panel_ec.tsv")
N2EC = (HERE / ".." / "0_building" / "3_nutrient_to_ec.tsv").resolve()
NUTRIENT_CSV = Path("/data/bac2food/nutrient.csv")
IDX = Path("/data/bac2food/index_modeled")
TOPK = 14
FOOD_W = 0.0  # --specialty_food_idf_weight default

# The substrate classes the panel organisms are documented for. Deliberately broad:
# the question is whether ANYTHING in the growth-substrate family survives selection.
SUBSTRATE = re.compile(
    r"\b(?:fib(?:er|re)|starch|amylose|amylopectin|maltodextrin|maltose|oligosaccharide|"
    r"inulin|fructan|pectin|glucan|xylan|arabinan|arabinoxylan|pentosan|hemicellulose|"
    r"cellulose|galactan|mucin|lactose|galactose|glucose|fucose|sialic|"
    r"n-acetylglucosamine|chitin|raffinose|stachyose|galactooligo|fructooligo)\b", re.I)

# What the rankings are ACTUALLY made of, per the panel output.
PHYTO = re.compile(
    r"\b(?:sterol|stanol|campesterol|stigmasterol|brassicasterol|avenasterol|sitosterol|"
    r"lycopene|lutein|carotene|cryptoxanthin|zeaxanthin|carotenoid|"
    r"quercetin|kaempferol|myricetin|morin|catechin|epicatechin|anthocyan|cyanidin|"
    r"delphinidin|malvidin|peonidin|petunidin|pelargonidin|flavon|flavan|isoflav|"
    r"genistein|daidzein|lignan|syringaresinol|secoisolariciresinol|pinoresinol|"
    r"psoralen|isopimpinellin|bergapten|nobiletin|sinensetin|tangeretin|"
    r"phenol|phenolic|coumaric|ferulic|caffeic|chlorogenic|ellagic|gallic|"
    r"resveratrol|oleuropein|hesperidin|naringenin|rutin|apigenin|luteolin)\b", re.I)


def nutrient_names() -> dict[int, str]:
    out = {}
    with NUTRIENT_CSV.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            try:
                out[int(r.get("id") or r.get("nutrient_id"))] = r.get("name") or r.get("nutrient_name") or ""
            except (TypeError, ValueError):
                continue
    return out


def klass(name: str) -> str:
    if SUBSTRATE.search(name):
        return "substrate"
    if PHYTO.search(name):
        return "phytochem"
    return "other"


def main() -> int:
    if not PANEL_EC.exists():
        sys.exit(f"[!] {PANEL_EC} missing — run validate_biology_panel.py first.")

    # EC -> nutrients
    ec2n = defaultdict(set)
    with N2EC.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                ec2n[r["ec_number"]].add(int(r["nutrient_id"]))
            except (KeyError, TypeError, ValueError):
                continue

    # species -> ECs, from the exact input the panel fed the predictor
    ec_by_lbl = defaultdict(set)
    with PANEL_EC.open(encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                ec_by_lbl[p[1]].add(p[3])

    modeled = set(pd.read_parquet(IDX / "nutrient_df.parquet")["nutrient_id"].astype(int))
    ndf = pd.read_parquet(IDX / "nutrient_df.parquet").set_index("nutrient_id")["df_foods"].to_dict()
    F_total = max(1, int(pd.read_parquet(IDX / "modeled_totals.parquet").shape[0]))

    lbl2targ = {lbl: {n for e in ecs for n in ec2n.get(e, ()) if n in modeled}
                for lbl, ecs in ec_by_lbl.items()}
    lbl2targ = {k: v for k, v in lbl2targ.items() if v}

    bf: dict[int, int] = defaultdict(int)
    for targs in lbl2targ.values():
        for n in targs:
            bf[n] += 1
    B_total = max(1, len(lbl2targ))

    def spec(n: int) -> float:
        bn = max(1, bf.get(n, 1))
        fn = int(ndf.get(n, F_total)) or F_total
        return math.log(B_total / bn) + FOOD_W * (math.log(F_total / fn) if fn < F_total else 0.0)

    nm = nutrient_names()
    name = lambda n: nm.get(n, str(n))

    print("=" * 100)
    print(f"Specialty allowlist reconstruction — {B_total} species, top-{TOPK} by "
          f"bacterial-IDF + {FOOD_W}x food-IDF")
    print("=" * 100)

    tot = defaultdict(int)
    for lbl in sorted(lbl2targ):
        targs = lbl2targ[lbl]
        ranked = sorted(targs, key=spec, reverse=True)
        keep = ranked[:TOPK]
        kinds = [klass(name(n)) for n in keep]
        for k in kinds:
            tot[k] += 1
        sp = lbl.split("_", 1)[1].replace("_", " ")
        n_sub_avail = sum(1 for n in targs if klass(name(n)) == "substrate")
        best_sub = next((i for i, n in enumerate(ranked, 1)
                         if klass(name(n)) == "substrate"), None)
        print(f"\n  {sp}   ({len(targs)} targets, {n_sub_avail} of them growth substrates)")
        print(f"    kept: {kinds.count('substrate')} substrate / {kinds.count('phytochem')} phytochemical "
              f"/ {kinds.count('other')} other")
        print(f"    best-ranked growth substrate sits at position "
              f"{best_sub if best_sub else '—'} of {len(ranked)}; cutoff is {TOPK}")
        for n in keep[:6]:
            print(f"       spec {spec(n):5.2f}  bf {bf[n]:3}/{B_total}  "
                  f"[{klass(name(n)):9}] {name(n)[:52]}")

    print("\n" + "-" * 100)
    n_kept = sum(tot.values())
    print(f"  ACROSS ALL SPECIES, of {n_kept} allowlist slots:")
    for k in ("substrate", "phytochem", "other"):
        print(f"    {k:10} {tot[k]:4}  ({100*tot[k]/max(1,n_kept):.1f}%)")

    # How rare is a growth substrate, by construction?
    allt = set().union(*lbl2targ.values())
    sub = [n for n in allt if klass(name(n)) == "substrate"]
    phy = [n for n in allt if klass(name(n)) == "phytochem"]
    def mbf(v):
        return sum(bf[n] for n in v) / max(1, len(v))
    print(f"\n  mean bf (species sharing the nutrient), out of {B_total}:")
    print(f"    growth substrates  {mbf(sub):5.1f}   (n={len(sub)})")
    print(f"    phytochemicals     {mbf(phy):5.1f}   (n={len(phy)})")
    print("  Selection is on 1/bf. A substrate shared by many organisms cannot win it.")
    print("-" * 100)
    return 0


if __name__ == "__main__":
    sys.exit(main())
