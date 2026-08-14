#!/usr/bin/env python3
"""rarefy_richness.py — the depth control behind Supplementary Note S3.

Detected richness tracks annotation depth almost perfectly, so the 6-to-12-month rise in
species per sample has to be shown not to be an artifact of the 12-month samples simply
being annotated more deeply. The test: subsample every sample down to a common number of
annotated loci and recompute the ratio. A ratio that stays flat across the whole depth
gradient is not produced by depth.

Until now this analysis existed only as numbers in the manuscript -- no script in the
repository reproduced it, and the numbers were computed on a superseded predictor run. This
file is that missing artifact, and it reads the same inputs as saturation_6_12.py so
"detected species" means the same thing in both.

    python rarefy_richness.py [--pred DIR] [--reps 200] [--seed 0]
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import random
import re
import statistics as st
from collections import defaultdict

ESC = "/data/bac2food/exports/enzyme_substrate_chebi.tsv"
ANNOT = "/data/bac2food/gene_annot_rescued"


def load_ec2nut(path: str) -> dict[str, set[str]]:
    ec2nut: dict[str, set[str]] = defaultdict(set)
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["in_model"] not in ("1", "True", "true", "yes"):
                continue
            for n in re.split(r"[;,]", (row["nutrient_ids"] or "").strip()):
                if n.strip():
                    ec2nut[row["ec_number"]].add(n.strip())
    return ec2nut


def load_samples(pred: str, ec2nut: dict[str, set[str]]):
    """Per sample: the list of (species, has_linked_nutrient) loci the predictor scored.

    Rarefaction has to act on LOCI, not on species: subsampling species would assume the
    answer. Each locus is kept with its species label so that a subsample of loci yields the
    species that survive at that depth.
    """
    out = []
    for path in sorted(glob.glob(os.path.join(ANNOT, "*_ec.tsv"))):
        sample = os.path.basename(path)[:-len("_ec.tsv")]
        m = re.search(r"_(\d+)_months?$", sample)
        if not m:
            continue
        pb = os.path.join(pred, f"{sample}.perBacterium.tsv")
        if not os.path.exists(pb):
            continue
        with open(pb, encoding="utf-8", newline="") as fh:
            scored = {r["bacterium"] for r in csv.DictReader(fh, delimiter="\t")}
        loci = []
        with open(path, encoding="utf-8", newline="") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 4 or f[1] not in scored:
                    continue
                if any(ec.strip() in ec2nut for ec in f[3].split(",")):
                    loci.append(f[1])
        if loci:
            out.append((sample, int(m.group(1)), loci))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", default="/data/bac2food/cohort_cov_phase11")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    ec2nut = load_ec2nut(ESC)
    samples = load_samples(args.pred, ec2nut)
    six = [s for s in samples if s[1] == 6]
    twelve = [s for s in samples if s[1] == 12]
    depths = [len(s[2]) for s in samples]
    print(f"[*] {len(samples)} samples | 6 mo n={len(six)} | 12 mo n={len(twelve)}")
    print(f"[*] linked loci per sample: {min(depths):,} to {max(depths):,} "
          f"(median {st.median(depths):,.0f})")

    # Depth-richness coupling: the reason the control is needed at all.
    rich = [len(set(s[2])) for s in samples]
    n = len(depths)
    rx = {v: i for i, v in enumerate(sorted(depths))}
    ry = {v: i for i, v in enumerate(sorted(rich))}
    dx = [rx[v] for v in depths]
    dy = [ry[v] for v in rich]
    mx, my = st.mean(dx), st.mean(dy)
    cov = sum((a - mx) * (b - my) for a, b in zip(dx, dy))
    rho = cov / ((sum((a - mx) ** 2 for a in dx) * sum((b - my) ** 2 for b in dy)) ** 0.5)
    print(f"[*] Spearman rho(depth, richness) = {rho:.3f}\n")

    full = min(depths)
    grid = [d for d in (250, 500, 1000, 2000, 4000, 8000) if d <= full] + [full]
    grid = sorted(set(grid))
    print(f"{'depth':>8}{'6 mo':>9}{'12 mo':>9}{'ratio':>8}")
    print("-" * 34)
    ratios = []
    for d in grid:
        rng = random.Random(args.seed + d)
        med = {}
        for label, grp in (("6", six), ("12", twelve)):
            per_sample = []
            for _, _, loci in grp:
                vals = [len(set(rng.sample(loci, d))) for _ in range(args.reps)]
                per_sample.append(st.mean(vals))
            med[label] = st.median(per_sample)
        r = med["12"] / med["6"]
        ratios.append(r)
        print(f"{d:>8,}{med['6']:>9.1f}{med['12']:>9.1f}{r:>8.2f}")
    print("-" * 34)
    print(f"ratio across the depth gradient: {min(ratios):.2f} to {max(ratios):.2f}")
    print("A ratio that does not move with depth is not produced by depth.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
