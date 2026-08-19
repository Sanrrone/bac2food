#!/usr/bin/env python3
"""validate_biology_panel.py — a biology gate strong enough to CHOOSE a scoring formula.

The six-species panel (validate_six_species_panel.py) exists to reproduce what the paper
reports, and it does that job. It cannot do this one. It is binary (a keyword appears in the
top ten or it does not) and it covers three testable species, so every formula variant tried
so far scores 3/3 on it — including one that admits sausages, fast food and alcohol into the
rankings. A gate that passes everything cannot select anything.

This panel differs in three ways:

  * SIZE. Eighteen species with documented substrate specialisms, all present in the shipped
    reference with full EC sets, instead of three.
  * GRADED. Score is the reciprocal rank of the first documented substrate hit (rank 1 -> 1.0,
    rank 4 -> 0.25, absent -> 0), so a formula that puts the right substrate first beats one
    that buries it at rank nine. Binary presence cannot see that difference.
  * PLAUSIBILITY. Reciprocal rank alone can be gamed by a formula that ranks aggressively and
    admits anything, so the whole-plant fraction and the junk/meat fraction are reported
    beside it. A formula that wins on MRR while filling the table with luncheon meat has not
    won.

METHODOLOGICAL COMMITMENT: the keyword sets below are fixed from documented enzymology and
must not be adjusted after seeing how a formula scores. This gate selects the formula; tuning
it to the formula would invert that and make the result meaningless. Each entry carries the
activity it encodes so the choice is auditable.

    python validate_biology_panel.py --differential_formula full|explicit_admission|gain_only
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

# species -> (documented activity, substrate keywords)
# Positives: the organism is characterized for utilizing this substrate class.
PANEL: dict[str, tuple[str, list[str]]] = {
    "Bifidobacterium longum": (
        "human-milk oligosaccharide utilization",
        ["lacto-n-", "fucosyllactose", "sialyllactose", "oligosaccharide", "lactose", "human milk"]),
    "Bifidobacterium bifidum": (
        "HMO and host mucin glycans",
        ["lacto-n-", "fucosyllactose", "sialyllactose", "mucin", "n-acetylglucosamine",
         "oligosaccharide", "lactose", "human milk"]),
    "Bifidobacterium breve": (
        "HMO utilization",
        ["lacto-n-", "fucosyllactose", "sialyllactose", "oligosaccharide", "lactose", "human milk"]),
    "Bifidobacterium adolescentis": (
        "starch and fructan utilization",
        ["starch", "amylopectin", "amylose", "inulin", "fructooligo", "fos", "maltodextrin"]),
    "Akkermansia muciniphila": (
        "mucin degradation",
        ["mucin", "n-acetylglucosamine", "glcnac", "sialic", "fucose", "chitin"]),
    "Faecalibacterium prausnitzii": (
        "inulin, FOS and pectin utilization",
        ["inulin", "fructooligo", "fos", "pectin", "apple", "chicory"]),
    "Roseburia intestinalis": (
        "arabinoxylan and beta-glucan utilization",
        ["arabinoxylan", "xylan", "beta-glucan", "b-glucan", "glucan", "bran", "barley", "oat"]),
    "Bacteroides thetaiotaomicron": (
        "broad plant-polysaccharide utilization",
        ["xylan", "pectin", "arabinan", "arabinoxylan", "polysaccharide", "galactan", "pentosan"]),
    "Bacteroides ovatus": (
        "xylan and hemicellulose utilization",
        ["xylan", "hemicellulose", "arabinoxylan", "pentosan", "polysaccharide"]),
    "Bacteroides uniformis": (
        "xylan utilization",
        ["xylan", "arabinoxylan", "pentosan", "polysaccharide", "hemicellulose"]),
    "Bacteroides fragilis": (
        "host and dietary glycan utilization",
        ["mucin", "sialic", "n-acetylglucosamine", "fucose", "polysaccharide", "galactan"]),
    "Ruminococcus bromii": (
        "resistant-starch degradation",
        ["starch", "amylose", "amylopectin", "maltodextrin", "maltose"]),
    "Clostridium butyricum": (
        "starch fermentation",
        ["starch", "amylopectin", "amylose", "maltodextrin"]),
    "Collinsella aerofaciens": (
        "starch utilization",
        ["starch", "amylopectin", "amylose", "maltose"]),
    "Streptococcus thermophilus": (
        "lactose fermentation (dairy)",
        ["lactose", "galactose", "milk", "dairy", "yogurt", "yoghurt", "cheese", "whey"]),
    "Flavonifractor plautii": (
        "flavonoid C-ring cleavage",
        ["flavon", "quercetin", "kaempferol", "catechin", "luteolin", "apigenin",
         "isoflav", "genistein", "daidzein", "soybean"]),
    "Lacticaseibacillus rhamnosus": (
        "deglycosylation of dietary phenolic conjugates",
        ["quercetin", "kaempferol", "glycoside", "glucoside", "rutin", "phenol", "flavon"]),
}

# Negative control. Named for inulin utilization, but it carries only 3 of the 17 EC numbers
# inulin maps to — fewer than several of its peers — so at EC granularity it holds no
# differential advantage on that substrate. A formula that "finds" inulin here is matching
# noise, so this scores in REVERSE: absence is the correct answer.
NEGATIVE: dict[str, tuple[str, list[str]]] = {
    "Roseburia inulinivorans": (
        "named for inulin, but lacks the EC depth to show it",
        ["inulin", "fructooligo", "chicory", "jerusalem artichoke"]),
}

PEERS = ["Enterococcus avium", "Enterocloster citroniae", "Veillonella parvula",
         "Escherichia coli", "Klebsiella pneumoniae"]

PLANT_CATS = {"Fruits and Fruit Juices", "Vegetables and Vegetable Products",
              "Cereal Grains and Pasta", "Legumes and Legume Products",
              "Nut and Seed Products"}
JUNK_CATS = {"Alcoholic Beverages", "Fast Foods", "Restaurant Foods", "Sweets", "Snacks",
             "Sausages and Luncheon Meats", "Beef Products", "Pork Products",
             "Poultry Products", "Fats and Oils", "Soups, Sauces, and Gravies",
             "Baked Products", "Lamb, Veal, and Game Products"}


def load_reference(path: Path, species: list[str]):
    ec, tax = defaultdict(set), {}
    want = set(species)
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["species"] in want:
                ec[r["species"]].add(r["ec_number"])
                tax.setdefault(r["species"], r["tax_id"])
    return ec, tax


def main() -> int:
    ap = argparse.ArgumentParser()
    # default None, not "full": the gate must test what a user actually gets. Hardcoding a
    # value here silently overrode the predictor's own default and made one sweep row a
    # measurement of the wrong formula.
    ap.add_argument("--differential_formula",
                    choices=["full", "explicit_admission", "gain_only"], default=None)
    ap.add_argument("--out_dir", default=None)
    ap.add_argument("--config", default=str(HERE / "parameters.yaml"))
    ap.add_argument("--reference", default="/data/bac2food/exports/species_enzymes.tsv")
    ap.add_argument("--topk", type=int, default=10,
                    help="GRADING depth: how many ranked rows are searched for the documented "
                         "substrate. Keep this FIXED when sweeping --max_foods, or a longer "
                         "list scores higher purely by having more chances to hit.")
    ap.add_argument("--max_foods", type=int, default=None,
                    help="Predictor scan budget / row cap. Defaults to --topk. Set it "
                         "independently to ask whether a deeper scan improves the TOP of the "
                         "list rather than just lengthening it.")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--label", default=None,
                    help="Names the output directory. Use it to keep sweep variants apart.")
    ap.add_argument("--predict_arg", action="append", default=[], metavar="ARG",
                    help="Passed through to bac2food_predict.py verbatim, repeatable. Lets the "
                         "gate sweep scoring parameters without a second copy of the panel.")
    args = ap.parse_args()

    out_dir = Path(args.out_dir or
                   f"/data/bac2food/runs/bio_panel/bio_panel_{args.label or args.differential_formula or 'default'}")
    out_dir.mkdir(parents=True, exist_ok=True)
    species = list(PANEL) + list(NEGATIVE) + PEERS
    ec, tax = load_reference(Path(args.reference), species)
    present = [s for s in species if ec.get(s)]
    absent = [s for s in species if not ec.get(s)]
    if absent:
        print(f"[!] absent from reference, skipped: {absent}")

    mag = out_dir / "panel_ec.tsv"
    with mag.open("w", encoding="utf-8") as fh:
        for sp in present:
            lbl = f"{tax[sp]}_{sp.replace(' ', '_')}"
            for i, e in enumerate(sorted(ec[sp])):
                fh.write(f"panel\t{lbl}\tl_{i:06d}\t{e}\n")
    print(f"[*] {len(present)} species, {sum(len(ec[s]) for s in present):,} EC rows -> {mag}")

    prefix = out_dir / "panel"
    diff = Path(str(prefix) + ".differential.tsv")
    if not (args.reuse and diff.exists()):
        # `differential` subcommand: this panel grades the differential table only, so the
        # community pass would be pure waste (it is ~40% of a combined run).
        cmd = [sys.executable, str(HERE / "bac2food_predict.py"), "differential",
               "--config", args.config,
               "--mag", str(mag), "--out", str(prefix),
               "--max_foods", str(args.max_foods or args.topk),
               "--jobs", str(args.jobs)]
        if args.differential_formula:
            cmd += ["--differential_formula", args.differential_formula]
        cmd += args.predict_arg
        if subprocess.call(cmd, cwd=str(HERE)) != 0:
            return 1

    rows = defaultdict(list)
    for r in csv.DictReader(diff.open(encoding="utf-8"), delimiter="\t"):
        rows[r["bacterium"]].append(r)
    for v in rows.values():
        v.sort(key=lambda r: int(r["rank"]))

    print("\n" + "=" * 96)
    print(f"Biology panel — {args.label or args.differential_formula or 'default'}: "
          f"differential_formula={args.differential_formula or '(predictor default)'} "
          f"{' '.join(args.predict_arg)}")
    print("=" * 96)
    print(f"  {'species':32} {'rank':>5} {'1/rank':>7}  documented substrate / what surfaced")
    print("  " + "-" * 92)

    rr, n_scored, hits1, hits3, n_empty = [], 0, 0, 0, 0
    for sp, (activity, kws) in PANEL.items():
        lbl = f"{tax.get(sp,'')}_{sp.replace(' ', '_')}"
        v = rows.get(lbl, [])[:args.topk]
        n_scored += 1
        if not v:
            # Producing no rows is a FAILURE, not an absence of evidence. Skipping these
            # would shrink the denominator and let a variant that starves most organisms
            # of output post a higher MRR than one that answers for all of them.
            n_empty += 1
            rr.append(0.0)
            print(f"  {sp:32} {'—':>5} {0.0:7.3f}  NO ROWS ({activity})")
            continue
        hit = None
        for r in v:
            blob = f"{r.get('food_name','')} {r.get('description','')} {r.get('top_nutrient_names','')}".lower()
            kw = next((k for k in kws if k in blob), None)
            if kw:
                hit = (int(r["rank"]), kw, r.get("food_name", ""))
                break
        if hit:
            rr.append(1.0 / hit[0]); hits1 += hit[0] == 1; hits3 += hit[0] <= 3
            print(f"  {sp:32} {hit[0]:5} {1.0/hit[0]:7.3f}  {hit[1]} <- {hit[2][:34]}")
        else:
            rr.append(0.0)
            print(f"  {sp:32} {'—':>5} {0.0:7.3f}  MISS ({activity}); top: {v[0].get('food_name','')[:30]}")

    neg_ok = 0
    for sp, (why, kws) in NEGATIVE.items():
        lbl = f"{tax.get(sp,'')}_{sp.replace(' ', '_')}"
        v = rows.get(lbl, [])[:args.topk]
        found = any(next((k for k in kws if k in
                          f"{r.get('food_name','')} {r.get('description','')} {r.get('top_nutrient_names','')}".lower()),
                         None) for r in v)
        neg_ok += (not found)
        print(f"  {sp:32} {'neg':>5} {'ok' if not found else 'FALSE+':>7}  {why}")

    top = [r for v in rows.values() for r in v[:args.topk]]
    plant = sum(1 for r in top if r.get("food_category") in PLANT_CATS)
    junk = sum(1 for r in top if r.get("food_category") in JUNK_CATS)
    mrr = sum(rr) / len(rr) if rr else 0.0
    print("  " + "-" * 92)
    print(f"  MRR {mrr:.3f} over {n_scored} species   hits@1 {hits1}/{n_scored}   "
          f"hits@3 {hits3}/{n_scored}   negative control {neg_ok}/{len(NEGATIVE)}"
          + (f"   [{n_empty} species produced NO ROWS]" if n_empty else ""))
    print(f"  plausibility: whole-plant {plant}/{len(top)} ({100*plant/max(1,len(top)):.1f}%)   "
          f"junk-or-meat {junk}/{len(top)} ({100*junk/max(1,len(top)):.1f}%)")
    print("=" * 96)
    return 0


if __name__ == "__main__":
    sys.exit(main())
