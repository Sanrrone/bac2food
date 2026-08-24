#!/usr/bin/env python3
"""chain_filter.py - the one definition of "reaches the chain", for every reader.

Three components read the food store and MUST agree about what exists:

    export_resources.py     writes the deposited food_nutrients.tsv
    prune_bucketed_store.py writes the store the predictor scores against
    bac2food_predict.py     scores

If they disagree, the scorer can rank a food the published table does not contain, or
lose a value the table still advertises. That has happened here before, which is why
this lives in one module instead of being spelled out three times.

A nutrient reaches the chain in one of two ways:

  1. DIRECTLY - some EC number acts on it. These are the targets, and they are exactly
     the nutrient_ids in 3_nutrient_to_ec.tsv.

  2. AS A SUBSTITUTE - the scoring kernel's proximity map. When a target is absent from
     a food, the kernel falls back to a chemically related generic and substitutes it at
     reduced efficiency: Pentosan stands in for Xylan, Cellulose for Cellobiose, Starch
     for Pullulan. A substitute carries no EC of its own, so rule 1 does not see it, but
     deleting it silently removes the fallback and the food scores as if the substrate
     were simply not there. Dropping these 69 nutrients cut the candidate pool of one
     cohort sample from 64,387 rows to 35,066 - a 46% loss that looked like a result.

The substitute half is derived from the same two sources the predictor builds `prox`
from - 1_expanded_nutrients.tsv plus the hardcoded form and bacterial aliases - so a
change there reaches the exporters without anyone remembering to mirror it.
"""
from __future__ import annotations

import csv
from pathlib import Path

# Kept in step with the `prox` construction in bac2food_predict.py. Each entry says
# "a food missing <target> may be scored on <generic> instead".
_FORM_ALIASES: list[tuple[list[int], int]] = [
    ([1017, 1021, 1022, 1403, 2058, 1071, 1019], 1079),
    ([1015, 1016, 1020], 1009),
    ([1181, 1182, 1042], 99999),
]

# extra_bacterial_seeds substrates that have no FDC nutrient of their own.
_BACTERIAL_ALIASES: list[tuple[int, int]] = [
    (200001, 96310), (200002, 96310),
    (200007, 1019), (200007, 1021), (200007, 1074),
    (200008, 1019), (200008, 1073),
    (200009, 1069), (200009, 2064),
    (200010, 1022),
    (200011, 1009),
]


def chain_targets(nutrient_to_ec: Path | str) -> set[int]:
    """Nutrients some EC number acts on."""
    with open(nutrient_to_ec, encoding="utf-8") as fh:
        return {int(r["nutrient_id"]) for r in csv.DictReader(fh, delimiter="\t")
                if r.get("nutrient_id")}


def chain_ec(nutrient_to_ec: Path | str) -> set[str]:
    """EC numbers that reach at least one nutrient."""
    with open(nutrient_to_ec, encoding="utf-8") as fh:
        return {r["ec_number"] for r in csv.DictReader(fh, delimiter="\t")
                if r.get("ec_number")}


def chain_nutrients(nutrient_to_ec: Path | str, nutrient_alias: Path | str) -> set[int]:
    """Targets PLUS every generic the kernel may substitute for one of them."""
    targets = chain_targets(nutrient_to_ec)
    keep = set(targets)

    with open(nutrient_alias, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # prox is keyed by the SPECIFIC id and holds generics, so a target that
            # appears as a `specific` can be served by its `generic`.
            try:
                s, g = int(r["specific_nutrient_id"]), int(r["generic_nutrient_id"])
            except (KeyError, TypeError, ValueError):
                continue
            if s in targets:
                keep.add(g)

    for targs, generic in _FORM_ALIASES:
        if any(t in targets for t in targs):
            keep.add(generic)
    for specific, generic in _BACTERIAL_ALIASES:
        if specific in targets:
            keep.add(generic)
    return keep
