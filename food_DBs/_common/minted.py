#!/usr/bin/env python3
"""Keep a source's minted nutrient ids stable across ingest re-runs.

A column the FDC catalogue has no id for gets one minted from the source's
reserved block. Every ingester did that POSITIONALLY - first unmapped column
gets START, second START+1 - so adding, removing or renaming one column
renumbered every id after it.

That is not a cosmetic problem. `enzyme_substrate_chebi.tsv` references nutrient
ids by number, and so does anything else built on top of the catalogue. A re-run
that shifted the block by one silently pointed twenty of those references at
nothing.

So the assignment is a REGISTRY, not a counter: the source's own
`*_extra_nutrient_map.tsv` is read first and every label already in it keeps the
id it was given. Only genuinely new labels are minted, above the block's current
high-water mark, in sorted order so a re-run of the same input is deterministic.
This is the same contract `fdc_blocks.assign` gives food accessions.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Tuple

import pandas as pd


def load_registry(map_path, label_col: str = "source_column_raw") -> Dict[str, int]:
    """{label: nutrient_id} from a source's existing extra-nutrient map."""
    p = Path(map_path)
    if not p.exists():
        return {}
    d = pd.read_csv(p, sep="\t", dtype=str)
    if label_col not in d.columns or "nutrient_id" not in d.columns:
        return {}
    out: Dict[str, int] = {}
    for lab, nid in zip(d[label_col], d["nutrient_id"]):
        if pd.notna(lab) and pd.notna(nid):
            out[str(lab)] = int(nid)
    return out


def assign(labels: Iterable[str], registry: Dict[str, int], start: int,
           limit: int | None = None) -> Tuple[Dict[str, int], list[str]]:
    """Return ({label: id}, newly_minted_labels), reusing every known label.

    `registry` is updated in place, so a caller can mint across several sheets
    and still hand out one id per label.
    """
    labels = [str(x) for x in labels]
    known = {lab: registry[lab] for lab in labels if lab in registry}
    nxt = max([*registry.values(), start - 1]) + 1
    fresh = sorted({lab for lab in labels if lab not in registry})
    for lab in fresh:
        if limit is not None and nxt >= limit:
            raise OverflowError(
                f"minted-nutrient block [{start:,}, {limit:,}) is full; widen it "
                f"rather than spilling into the next source's block")
        registry[lab] = nxt
        known[lab] = nxt
        nxt += 1
    return known, fresh
