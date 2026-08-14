#!/usr/bin/env python3
"""recompute_match_rate.py — the two cohort "match rates", measured separately.

"Match rate" had been used for two different quantities, which is why the old 34.4 %
could not be reproduced from either one. They are computed side by side here:

  1. FEATURE match rate — share of annotated LOCI reaching an FDC nutrient through
     EC -> substrate -> ChEBI -> nutrient (0_building/3_nutrient_to_ec.tsv). This is what
     the manuscript sentence "N % of annotated features match a reference-map entry"
     describes. It does NOT depend on the bacterium -> EC reference at all, so moving
     eggNOG v6 -> v7 leaves it unchanged; it moves only when the substrate/ChEBI layer or
     the cohort annotation changes. (Exception: --augment_with_reference, which injects
     reference ECs into under-annotated species, would shift it.)

     A locus may carry SEVERAL EC numbers in one cell, and it matches if any of them
     reaches a nutrient. Counting the raw cell as one EC understates the rate by ~9 pp
     (36.3 % vs the correct 45.4 %) -- see `split_ec_cell`.

  2. SPECIES match rate — share of the cohort's species labels that resolve to ANY
     reference EC set via the predictor's three-step ladder (tax_id -> normalized name ->
     Genus+species prefix). This is what --complement_ec prints, and it is exactly what
     the v6 -> v7 switch changes.

The species ladder is imported from the predictor rather than reimplemented, so the
number cannot drift from what the tool actually does.

Usage:
    python recompute_match_rate.py                       # v7 (current) vs v6 (legacy)
    python recompute_match_rate.py --ref /path/to/ref.tsv --label mine
"""
from __future__ import annotations

import argparse
import gc
import re
import importlib.util
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DEFAULT_ANNOT = REPO / "gene_annot"


def load_predictor():
    """Load bac2food_predict.py as a module (its name is not importable as-is)."""
    spec = importlib.util.spec_from_file_location("bac2food_predict", HERE / "bac2food_predict.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bac2food_predict"] = mod
    spec.loader.exec_module(mod)
    return mod


_EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
_EC_PREFIX_RE = re.compile(r"(?i)\bec:?\s*")
_EC_SPLIT_RE = re.compile(r"\s*[,;]\s*")


def split_ec_cell(cell: str) -> list[str]:
    """Split one annotation cell into its EC numbers.

    An annotated locus may carry several EC numbers in one cell ("2.7.2.3,5.3.1.1"):
    178,951 of the cohort's 1,257,938 loci do. `_normalize_ec_frame` in the predictor
    strips any "EC:" prefix, splits on comma/semicolon and explodes, so anything counting
    features here must do the same -- treating the raw cell as a single EC silently
    understates the match rate (it did, by ~9 pp, before this was fixed).
    """
    if not isinstance(cell, str):
        return []
    return [e for e in (x.strip() for x in _EC_SPLIT_RE.split(_EC_PREFIX_RE.sub("", cell)))
            if _EC_RE.match(e)]


def load_cohort(annot_dir: Path) -> pd.DataFrame:
    """Concatenate every gene_annot/*_ec.tsv sample.

    Headerless, 4 columns: sample, species, locus, ec_number. Column 1 is a constant
    placeholder label in every shipped file, so the sample identity is the FILENAME,
    not that column. EC cells are left raw here and split per distinct cell in
    `feature_match_rate`, which keeps the frame small enough to hold in memory.
    """
    files = sorted(annot_dir.glob("*_ec.tsv"))
    if not files:
        raise SystemExit(f"ERROR: no *_ec.tsv under {annot_dir}")
    parts = []
    for f in files:
        df = pd.read_csv(f, sep="\t", header=None, dtype="string",
                         names=["sample", "species", "locus", "ec_number"])
        df["sample"] = f.name
        parts.append(df.dropna(subset=["species", "ec_number"]))
    out = pd.concat(parts, ignore_index=True)
    out.attrs["n_files"] = len(files)
    return out


def feature_match_rate(cohort: pd.DataFrame, nutrient_to_ec: Path) -> dict:
    """Share of annotated features reaching an FDC nutrient (metric 1).

    A FEATURE is one annotated locus. It matches if ANY of its EC numbers reaches a
    nutrient, which is what "annotated features match a reference-map entry" means and
    what the predictor effectively does after exploding multi-EC cells. The exploded and
    distinct-EC views are reported too, since they answer different questions.

    The 1.26M loci carry only a few thousand DISTINCT EC cells, so everything is derived
    from the distinct cells and their frequencies. Materialising a per-row list column and
    exploding it instead costs several GB and gets the process OOM-killed on a 16 GB box.
    """
    n = pd.read_csv(nutrient_to_ec, sep="\t", dtype="string")
    model = set(n["ec_number"].dropna())

    freq = cohort["ec_number"].value_counts()
    cell_ecs = {c: split_ec_cell(c) for c in freq.index}
    cell_hit = {c: any(e in model for e in v) for c, v in cell_ecs.items()}
    cell_n = {c: len(v) for c, v in cell_ecs.items()}

    valid = cohort["ec_number"].map(lambda c: cell_n.get(c, 0) > 0)
    cohort = cohort[valid]
    freq = cohort["ec_number"].value_counts()

    hit = cohort["ec_number"].map(cell_hit)
    per_sample = cohort.assign(_h=hit).groupby("sample")["_h"].mean()

    pairs = pairs_hit = 0
    ec_counts: dict[str, int] = {}
    for cell, k in freq.items():
        for e in cell_ecs[cell]:
            pairs += k
            ec_counts[e] = ec_counts.get(e, 0) + k
            if e in model:
                pairs_hit += k
    ec_set = set(ec_counts)
    multi = int(sum(k for c, k in freq.items() if cell_n[c] > 1))

    return {"rows": int(hit.sum()), "rows_total": len(cohort),
            "rate": 100.0 * hit.mean(),
            "multi_ec_loci": multi,
            "per_sample_mean": 100.0 * per_sample.mean(),
            "per_sample_min": 100.0 * per_sample.min(),
            "per_sample_max": 100.0 * per_sample.max(),
            "pairs": pairs, "pairs_hit": pairs_hit,
            "ec_distinct": len(ec_set), "ec_hit": len(ec_set & model),
            "model_ecs": len(model)}


def main() -> int:
    ap = argparse.ArgumentParser(description="Recompute the cohort reference match rate.")
    ap.add_argument("--annot_dir", default=str(DEFAULT_ANNOT))
    ap.add_argument("--nutrient_to_ec", default=str(REPO / "0_building/3_nutrient_to_ec.tsv"))
    ap.add_argument("--ref", action="append", default=None,
                    help="Reference TSV (repeatable). Default: the v7 export.")
    ap.add_argument("--label", action="append", default=None, help="Label per --ref")
    args = ap.parse_args()

    # v7 only. The v6 layer was retired from the deposit; to reproduce the published
    # v6-vs-v7 comparison, regenerate it from the legacy source and pass it explicitly.
    # --out_dir is required; point it at scratch, NOT at the exports directory. The file is
    # named species_enzymes.tsv either way, so a v6 build written into the deposit replaces
    # the shipped v7 layer.
    #   python 5_export/export_resources.py --only species \
    #          --bact_ec /data/bac2food/bact_ec.tsv --out_dir /data/bac2food/v6_rebuild
    #   python recompute_match_rate.py --ref /data/bac2food/exports/species_enzymes.tsv \
    #          --ref /data/bac2food/v6_rebuild/species_enzymes.tsv \
    #          --label "eggNOG v7" --label "eggNOG v6"
    refs = args.ref or ["/data/bac2food/exports/species_enzymes.tsv"]
    labels = args.label or (["eggNOG v7 (via KEGG KO)"]
                            if not args.ref else [Path(r).name for r in refs])

    cohort = load_cohort(Path(args.annot_dir))
    species = sorted(cohort["species"].unique())
    print(f"[*] cohort: {cohort.attrs['n_files']} samples | {len(cohort):,} annotated EC rows "
          f"| {len(species):,} distinct species labels", flush=True)

    # ---- metric 1: feature match rate (independent of the bacterium -> EC reference) ----
    fm = feature_match_rate(cohort, Path(args.nutrient_to_ec))
    print(f"\n[1] FEATURE match rate (annotated loci reaching an FDC nutrient)")
    print(f"    {fm['rows']:,} / {fm['rows_total']:,} loci = {fm['rate']:.1f} %"
          f"   ({fm['multi_ec_loci']:,} loci carry >1 EC and are split)")
    print(f"    per sample: mean {fm['per_sample_mean']:.1f} % "
          f"(range {fm['per_sample_min']:.1f}-{fm['per_sample_max']:.1f} %)")
    print(f"    exploded (locus, EC): {fm['pairs_hit']:,} / {fm['pairs']:,} "
          f"= {100*fm['pairs_hit']/fm['pairs']:.1f} %")
    print(f"    distinct EC: {fm['ec_hit']:,} / {fm['ec_distinct']:,} "
          f"= {100*fm['ec_hit']/fm['ec_distinct']:.1f} % "
          f"(model covers {fm['model_ecs']:,} EC)")
    print("    -> unchanged by the eggNOG v6 -> v7 switch (see module docstring)")

    # ---- metric 2: species ladder match rate (this is what v7 changes) ----
    # Reduce the cohort to the two things metric 2 needs, then drop the frame: the ladder
    # loads the full 20.5M-row v7 reference and builds three dicts of EC sets from it, so
    # holding the annotation frame alongside it is what pushes this over a 16 GB box.
    weights = cohort["species"].value_counts().to_dict()
    total_rows = len(cohort)
    del cohort
    gc.collect()

    print(f"\n[2] SPECIES match rate (cohort species resolving to a reference EC set)")
    mod = load_predictor()
    results = []
    for ref, lab in zip(refs, labels):
        if not Path(ref).exists():
            print(f"    [!] skipping missing reference: {ref}", flush=True)
            continue
        matched = mod.match_user_to_reference(species, ref)
        hit = [s for s in species if matched.get(s)]
        # Weighted view: species carrying more annotated rows matter more than singletons.
        w = sum(weights[s] for s in hit)
        results.append((lab, len(hit), len(species), w))
        print(f"    {lab}", flush=True)
        print(f"      species matched : {len(hit):,} / {len(species):,} "
              f"({100*len(hit)/len(species):.1f} %)", flush=True)
        print(f"      EC rows covered : {w:,} / {total_rows:,} "
              f"({100*w/total_rows:.1f} %)", flush=True)

    if len(results) == 2:
        (l0, h0, n0, w0), (l1, h1, n1, w1) = results
        print(f"\n[*] delta ({l0} vs {l1}): "
              f"{100*h0/n0 - 100*h1/n1:+.1f} pp species, "
              f"{100*(w0-w1)/total_rows:+.1f} pp EC rows", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
