#!/usr/bin/env python3
"""taxon_level_test.py — build the UNPOOLED taxon profile matrix.

Figure 2b pools each taxon's scores over the infants it was detected in, giving one row
per taxon. That is the right matrix for an ordination (see the pooling note in
prep_figure_data.py) but it makes species untestable as a grouping factor: with one point
per species there are 54 groups of n=1, no within-group variance, and PERMANOVA returns
R^2 = 1 by construction. Family and genus are the only ranks the pooled matrix can test.

This writes the unpooled version instead: one row per (sample, taxon, food), so a species
detected in k samples contributes k replicate profiles and species becomes a factor with
real replication. Every taxon here has >= 5 replicates.

The cost is sparsity. A single (sample, taxon) profile is the predictor's truncated
top-N list — a median of 10 foods out of 737 — where the pooled profile has the union
over samples. The R script quantifies what that does to the distances rather than
assuming it is harmless.

Sample identifiers are anonymized by reusing prep_figure_data.anonymize(); the raw cohort
names must not reach any artifact, and importing rather than reimplementing keeps this
file on the same lineage as the figures.

    BAC2FOOD_PREDICT_DIR=/data/bac2food/cohort_cov_phase11 python3 taxon_level_test.py
"""
from __future__ import annotations

import csv
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from prep_figure_data import PREDICT, anonymize, family_of  # noqa: E402

OUTDIR = Path(os.environ.get("BAC2FOOD_FIGDIR", str(HERE)))

# Genus is the first token of the taxon name. These placeholders are not genera: pooling
# them would fabricate a clade out of taxa whose only shared property is being unresolved.
NOT_A_GENUS = {"uncultured", "unclassified", "candidatus"}


def main() -> int:
    label = anonymize()
    print(f"[*] predictor dir: {PREDICT}")
    print(f"[*] {len(label)} samples")

    rows: list[tuple[str, str, str, str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for stem, lab in label.items():
        p = PREDICT / f"{stem}.differential.tsv"
        if not p.exists():
            continue
        with p.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                taxon, fam = family_of(r["bacterium"])
                k = (lab, taxon, r["food_name"])
                if k in seen:                    # one row per (sample, taxon, food)
                    continue
                seen.add(k)
                genus = taxon.split()[0]
                if genus.lower() in NOT_A_GENUS:
                    genus = "Unclassified"
                rows.append((lab, taxon, genus, fam, r["food_name"]))

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / "taxon_unpooled_long.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "taxon", "genus", "family", "food"])
        w.writerows(rows)

    obs = {(s, t) for s, t, _, _, _ in rows}
    per_obs = Counter((s, t) for s, t, _, _, _ in rows)
    n_taxon = Counter(t for _, t in obs)
    sz = sorted(per_obs.values())
    print(f"  taxon_unpooled_long.csv  {len(rows):,} rows")
    print(f"  [check] {len(obs):,} (sample, taxon) observations, "
          f"{len(n_taxon)} taxa, {len({r[4] for r in rows})} foods")
    print(f"  [check] replicates per taxon: min {min(n_taxon.values())} "
          f"median {sorted(n_taxon.values())[len(n_taxon) // 2]} max {max(n_taxon.values())}")
    print(f"  [check] foods per observation: min {sz[0]} median {sz[len(sz) // 2]} "
          f"max {sz[-1]}  (of {len({r[4] for r in rows})} foods)")
    lvl = defaultdict(set)
    for _, t, g, f, _ in rows:
        lvl["species"].add(t); lvl["genus"].add(g); lvl["family"].add(f)
    print("  [levels] " + ", ".join(f"{k}:{len(v)}" for k, v in lvl.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
