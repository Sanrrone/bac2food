#!/usr/bin/env python3
"""saturation_6_12.py — the union-vs-per-species contrast across the feeding transition.

Usage Notes claims that between 6 and 12 months the UNION of enzyme-linked nutrients a
community can reach does not move, while demand summed PER SPECIES scales with richness.
Both arms have to come from one script or the contrast is not comparable, so everything the
manuscript reports for that paragraph is computed here from the deposited exports plus the
cohort annotations: richness, union, nutrient-species pairs, carriers per nutrient, the
per-species repertoire, and the carrier share.

Only the 598 nutrients that carry a measured value in the composition table count as
"enzyme-linked"; that restriction is what enzyme_substrate_chebi.in_model already encodes.
"""
import csv, glob, os, re, statistics as st, sys
from collections import defaultdict
from scipy.stats import mannwhitneyu

ESC = "/data/bac2food/exports/enzyme_substrate_chebi.tsv"
ANNOT = "/data/bac2food/gene_annot_rescued"
# "Detected species" must mean the same thing here as everywhere else in the paper: the species
# the predictor actually scored, not every label appearing in the annotation. The raw annotation
# carries low-support species the predictor drops, which inflates richness (34/66 rather than
# 30/58) and shifts every per-species quantity derived from it.
PRED = "/data/bac2food/cohort_cov_phase11"

# ---- EC -> {nutrient_id} over the in-model (food-reaching) rows only ------------
ec2nut = defaultdict(set)
with open(ESC, encoding="utf-8", newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["in_model"] not in ("1", "True", "true", "yes"):
            continue
        ids = (row["nutrient_ids"] or "").strip()
        if not ids:
            continue
        for n in re.split(r"[;,]", ids):
            n = n.strip()
            if n:
                ec2nut[row["ec_number"]].add(n)
linked = set().union(*ec2nut.values()) if ec2nut else set()
print(f"[*] {len(ec2nut):,} EC numbers reach {len(linked):,} enzyme-linked nutrients")

# ---- per sample: species -> {nutrient} -----------------------------------------
rows = []
for path in sorted(glob.glob(os.path.join(ANNOT, "*_ec.tsv"))):
    sample = os.path.basename(path)[:-len("_ec.tsv")]
    m = re.search(r"_(\d+)_months?$", sample)
    if not m:
        continue
    month = int(m.group(1))
    pb = os.path.join(PRED, f"{sample}.perBacterium.tsv")
    if not os.path.exists(pb):
        continue
    with open(pb, encoding="utf-8", newline="") as fh:
        scored = {r["bacterium"] for r in csv.DictReader(fh, delimiter="\t")}
    sp2nut = defaultdict(set)
    with open(path, encoding="utf-8", newline="") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[1] not in scored:
                continue
            sp, ecs = f[1], f[3]
            for ec in ecs.split(","):
                ec = ec.strip()
                if ec in ec2nut:
                    sp2nut[sp] |= ec2nut[ec]
    sp2nut = {s: n for s, n in sp2nut.items() if n}
    if not sp2nut:
        continue
    union = set().union(*sp2nut.values())
    carriers = defaultdict(int)
    for n in sp2nut.values():
        for x in n:
            carriers[x] += 1
    rich = len(sp2nut)
    rows.append(dict(
        sample=sample, month=month, richness=rich,
        union=len(union),
        pairs=sum(len(n) for n in sp2nut.values()),
        carriers=st.median(carriers.values()),
        repertoire=st.median(len(n) for n in sp2nut.values()),
        share=st.median(carriers.values()) / rich,
    ))

a = [r for r in rows if r["month"] == 6]
b = [r for r in rows if r["month"] == 12]
print(f"[*] {len(rows)} samples parsed | 6 mo n={len(a)} | 12 mo n={len(b)}\n")
print(f"{'measure':<14}{'6 mo':>10}{'12 mo':>10}{'fold':>8}{'p':>12}")
for k in ("richness", "union", "pairs", "carriers", "repertoire", "share"):
    x = [r[k] for r in a]
    y = [r[k] for r in b]
    mx, my = st.median(x), st.median(y)
    p = mannwhitneyu(x, y, alternative="two-sided").pvalue
    print(f"{k:<14}{mx:>10.4g}{my:>10.4g}{(my/mx if mx else 0):>8.2f}{p:>12.2g}")

allr = [r["richness"] for r in rows]
allu = [r["union"] for r in rows]
print(f"\nacross all {len(rows)} samples: richness {min(allr)}-{max(allr)} | "
      f"union {min(allu)}-{max(allu)} (median {st.median(allu):.0f}, "
      f"{100*st.median(allu)/len(linked):.1f}% of {len(linked)})")
