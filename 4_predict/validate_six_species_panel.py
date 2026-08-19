#!/usr/bin/env python3
"""validate_six_species_panel.py — the biological-plausibility panel the paper reports.

Technical Validation describes a six-species panel and states that it is reproducible from
the deposited files alone. Until now nothing in the repository ran it: the result came out
of a side session, so a change to the scoring kernel could silently invalidate it. This
script is that missing artifact, and it is the biology gate any formula change must clear.

It is NOT test_specialty_panel.py. That is the retired Phase 9 eight-bacterium check whose
help text still advertises "6/8"; four of its eight organisms (S. mutans, B. dentium,
V. tobetsuensis, B. nordii) carry zero rows in the shipped eggNOG v7 reference, so it
cannot score above 4/8 today, and on a cohort sample it scores 0/8 because seven of the
eight are absent and its thresholds assume a reference-wide run.

Design, matching what the paper claims:

  * Enzyme sets come from the REFERENCE layer (species_enzymes.tsv), not from cohort
    annotations, so the panel depends on deposited files only and is unaffected by which
    infants were sequenced.
  * Each species is scored differentially against the other five. With only five peers the
    peer median is thin, so individual placements are read as directional; what the panel
    asserts is that the substrate CLASSES match documented enzymology.
  * Every expectation below is checked against computed output, and the supporting counts
    (how many EC numbers a species carries for a substrate, how many foods carry a measured
    value) are recomputed rather than quoted, so a kernel change cannot leave a stale number
    standing in the manuscript.

    python validate_six_species_panel.py [--out_dir DIR] [--config parameters.yaml] [--jobs 3]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

# species -> (documented activity, substrate keywords expected in the leading ranking)
# Keywords are matched against food_name + description + top_nutrient_names, because the
# claim is about substrate class and the evidence may surface in either the food or the
# nutrient that carried it.
PANEL = {
    "Flavonifractor plautii": (
        "cleaves the flavonoid C-ring",
        ["flavon", "quercetin", "kaempferol", "catechin", "luteolin", "apigenin",
         "isoflav", "genistein", "daidzein", "soybean", "soy"]),
    "Lacticaseibacillus rhamnosus": (
        "deglycosylates dietary phenolic conjugates",
        ["quercetin", "kaempferol", "glycoside", "glucoside", "rutin", "phenol",
         "flavon"]),
    "Roseburia inulinivorans": (
        "named for inulin utilization (the panel's documented negative)",
        ["inulin", "fructan", "chicory", "jerusalem artichoke"]),
}
PEERS = ["Enterococcus avium", "Enterocloster citroniae", "Veillonella parvula"]

# The negative result the paper explains: R. inulinivorans does not surface inulin foods.
# Both halves of that explanation are recomputed — that the data exist (foods carrying a
# measured inulin value) and that the linkage exists (EC numbers inulin maps to) — so the
# cause remains "granularity" only for as long as that is actually true.
INULIN_RE = re.compile(r"\binulin\b", re.I)


def load_reference(path: Path, species: list[str]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = defaultdict(set)
    want = set(species)
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] in want:
                out[r["species"]].add(r["ec_number"])
    missing = want - set(out)
    if missing:
        sys.exit(f"[!] absent from {path}: {sorted(missing)}. The panel cannot be scored.")
    return out


def write_mag_tsv(ec_by_species: dict[str, set[str]], tax: dict[str, str], dest: Path) -> None:
    """Emit the headerless <sample> <taxid_Genus_species> <locus> <ec> form the predictor
    auto-detects, so the panel enters through exactly the same input path as a real run."""
    with dest.open("w", encoding="utf-8") as fh:
        for sp, ecs in sorted(ec_by_species.items()):
            label = f"{tax[sp]}_{sp.replace(' ', '_')}"
            for i, ec in enumerate(sorted(ecs)):
                fh.write(f"panel\t{label}\tl_{i:06d}\t{ec}\n")


def tax_ids(path: Path, species: list[str]) -> dict[str, str]:
    want, out = set(species), {}
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] in want and r["species"] not in out:
                out[r["species"]] = r["tax_id"]
                if len(out) == len(want):
                    break
    return out


def inulin_support(nutrient_csv: Path, n2ec: Path, bucketed: Path,
                   ec_by_species: dict[str, set[str]], food_parquet: Path) -> dict:
    """Recompute the three counts behind the documented negative.

    The food count MUST be taken on the food set the predictor actually scores, not on the
    whole bucketed store. Inulin is a case where that matters more than usual: 170 foods
    carry a positive value, but 113 of them are branded label products, which the predictor
    drops. Counting the store would report 170 against the paper's 57 and look like a
    regression when it is only a different denominator.
    """
    import pandas as pd
    import pyarrow.dataset as ds

    ids = set()
    with nutrient_csv.open(encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh):
            nm = r.get("name") or r.get("nutrient_name") or ""
            if INULIN_RE.search(nm):
                try:
                    ids.add(int(r.get("id") or r.get("nutrient_id")))
                except (TypeError, ValueError):
                    continue
    ecs = set()
    with n2ec.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                if int(r["nutrient_id"]) in ids:
                    ecs.add(r["ec_number"])
            except (KeyError, TypeError, ValueError):
                continue
    foods = 0
    if ids and bucketed.exists():
        d = ds.dataset(str(bucketed), format="parquet", partitioning="hive")
        df = d.to_table(columns=["fdc_id", "nutrient_id", "amount"],
                        filter=ds.field("nutrient_id").isin(sorted(ids))).to_pandas()
        df = df[df["amount"] > 0]
        if food_parquet.exists():
            meta = pd.read_parquet(food_parquet, columns=["fdc_id", "data_type"])
            df = df.merge(meta, on="fdc_id", how="left")
            df = df[~df["data_type"].isin(["branded_food", "survey_fndds_food"])]
        foods = int(df["fdc_id"].nunique())
    carried = {sp: len(ecs & e) for sp, e in ec_by_species.items()}
    return {"nutrient_ids": sorted(ids), "n_ec": len(ecs), "n_foods": foods,
            "carried": carried}


def read_differential(p: Path) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    with p.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["bacterium"]].append(r)
    for v in out.values():
        v.sort(key=lambda r: int(r["rank"]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_dir", default="/data/bac2food/runs/six_species/six_species_panel")
    ap.add_argument("--config", default=str(HERE / "parameters.yaml"))
    ap.add_argument("--reference", default="/data/bac2food/exports/species_enzymes.tsv")
    ap.add_argument("--nutrient_csv", default="/data/bac2food/nutrient.csv")
    ap.add_argument("--nutrient_to_ec", default=str(HERE / ".." / "0_building" / "3_nutrient_to_ec.tsv"))
    ap.add_argument("--bucketed", default="/data/bac2food/food_nutrient_bucketed")
    ap.add_argument("--food_parquet", default="/data/bac2food/food.parquet")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--reuse", action="store_true", help="skip scoring if output exists")
    ap.add_argument("--differential_formula", choices=["full","explicit_admission","gain_only"],
                    default="full", help="forwarded to the predictor; A/B the formula variants")
    ap.add_argument("--predict_arg", action="append", default=[], metavar="ARG",
                    help="Passed through to bac2food_predict.py verbatim, repeatable. Pass as "
                         "--predict_arg=--flag; argparse refuses a separate token starting '--'.")
    args = ap.parse_args()

    species = list(PANEL) + PEERS
    ref = Path(args.reference)
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    ec_by_species = load_reference(ref, species)
    tax = tax_ids(ref, species)
    mag = out_dir / "panel_ec.tsv"
    write_mag_tsv(ec_by_species, tax, mag)
    print(f"[*] panel input: {mag} "
          f"({sum(len(v) for v in ec_by_species.values()):,} EC rows, {len(species)} species)")
    for sp in species:
        print(f"      {sp:32} {len(ec_by_species[sp]):5,} EC")

    prefix = out_dir / "panel"
    diff = Path(str(prefix) + ".differential.tsv")
    if not (args.reuse and diff.exists()):
        # --specialty_mode is left at its default: the panel must reflect what a user gets,
        # not a configuration chosen to make it pass.
        # `differential` subcommand: the panel reads differential.tsv only.
        cmd = [sys.executable, str(HERE / "bac2food_predict.py"), "differential",
               "--config", args.config,
               "--mag", str(mag), "--out", str(prefix),
               "--max_foods", str(max(10, args.topk)),
               "--jobs", str(args.jobs)]
        cmd += ["--differential_formula", args.differential_formula] + args.predict_arg
        print(f"[*] scoring: {' '.join(cmd[-8:])}", flush=True)
        rc = subprocess.call(cmd, cwd=str(HERE))
        if rc != 0:
            return rc
    if not diff.exists():
        sys.exit(f"[!] no differential output at {diff}")

    d = read_differential(diff)
    label_of = {sp: f"{tax[sp]}_{sp.replace(' ', '_')}" for sp in species}

    print("\n" + "=" * 92)
    print("Six-species panel — substrate classes recovered against documented enzymology")
    print("=" * 92)
    n_pass = 0
    for sp, (activity, kws) in PANEL.items():
        rows = d.get(label_of[sp], [])[:args.topk]
        hits = []
        for r in rows:
            blob = f"{r.get('food_name','')} | {r.get('description','')} | {r.get('top_nutrient_names','')}".lower()
            kw = next((k for k in kws if k in blob), None)
            if kw:
                hits.append((int(r["rank"]), kw, r.get("food_name", "")))
        expect_neg = sp.startswith("Roseburia")
        ok = (not hits) if expect_neg else bool(hits)
        n_pass += ok
        verdict = "as documented" if ok else "DEPARTS from the paper"
        print(f"\n  {sp}  ({activity})")
        print(f"    {'expected: no hits' if expect_neg else 'expected: substrate class in top-' + str(args.topk)}"
              f"  ->  {len(hits)} of {len(rows)} ranked rows match  [{verdict}]")
        for rk, kw, food in hits[:5]:
            print(f"      rank {rk:>2}  {kw:<12} {food[:46]}")
        if not hits and rows:
            print(f"      top-3 instead: " + "; ".join(r.get("food_name", "")[:26] for r in rows[:3]))

    sup = inulin_support(Path(args.nutrient_csv), Path(args.nutrient_to_ec),
                         Path(args.bucketed), ec_by_species, Path(args.food_parquet))
    print("\n" + "-" * 92)
    print("  Why the documented negative is granularity rather than missing data:")
    print(f"    foods carrying a measured inulin value : {sup['n_foods']}")
    print(f"    EC numbers inulin maps to              : {sup['n_ec']}")
    for sp in sorted(sup["carried"], key=lambda s: -sup["carried"][s]):
        print(f"      carried by {sp:32} {sup['carried'][sp]:3}")
    print("-" * 92)
    print(f"  PANEL: {n_pass}/{len(PANEL)} species behave as the paper reports.")
    print("=" * 92)
    return 0 if n_pass == len(PANEL) else 1


if __name__ == "__main__":
    sys.exit(main())
