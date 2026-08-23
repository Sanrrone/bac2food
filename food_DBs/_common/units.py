#!/usr/bin/env python3
"""Reconcile a source table's units with FDC's canonical unit for the nutrient.

Every ingester maps a source column onto an FDC `nutrient_id`, and the export
then labels the value with FDC's unit for that id - it does not carry the
source's own unit anywhere. So a column whose unit differs from FDC's is
published under the wrong one, silently.

That is what happened to the amino acids. FDC reports them in G; STFCJ, AFCD and
BioFoodComp report them in mg. 43,000 rows were deposited a thousand times too
high - BioFoodComp's quinoa read 873 g of leucine per 100 g - and nothing in the
pipeline could notice, because the number is plausible until you look at the
unit it was measured in.

Use it like this:

    from _common.units import fdc_units, conversion_factor
    UNITS = fdc_units(nutrient_csv)          # {nutrient_id: 'G' | 'MG' | ...}
    f = conversion_factor(src_unit, UNITS[nutrient_id])
    amount = amount * f

`conversion_factor` returns 1.0 when the two units agree, when either is unknown
or when they are not commensurable - it never guesses. A source unit it cannot
parse is reported by `unknown_units()` so the ingester can print it rather than
convert on a hunch.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable

import pandas as pd

# Everything is per 100 g edible portion, so only the SCALE of the unit matters.
# Values are the multiplier that takes the unit to grams.
_TO_GRAM = {
    "kg": 1e3,
    "g": 1.0,
    "mg": 1e-3,
    "ug": 1e-6, "µg": 1e-6, "mcg": 1e-6, "μg": 1e-6,
    "ng": 1e-9,
}
# Units that are not a mass at all. Two of these never convert into each other
# or into a mass, so the factor stays 1.0 and the ingester leaves the value be.
_NON_MASS = {"kcal", "kj", "iu", "re", "rae", "dfe", "ne", "ate", "%", "ph",
             "l", "ml", "cl", "dl", "unit", "score"}

_UNKNOWN: set[str] = set()

# A column can name a basis that is not "per 100 g of food" in its LABEL rather
# than in its unit: McCance reports each fatty acid twice, "/100g food" and
# "/100g fa" (per 100 g of total fatty acids), and its tryptophan column is
# "Tryptophan/60". AFCD writes "(mg/gN)", per gram of nitrogen. None of these
# can be scaled onto the per-100 g basis - they are different measurements, not
# different units - so a column matching this must not be mapped onto an FDC id
# at all. Mint it instead, and the basis travels in the extra nutrient table.
NON_FOOD_BASIS_RE = re.compile(
    r"/\s*100\s*g\s*(?:fa|fatty)\b"      # per 100 g of fatty acids
    r"|/\s*g\s*n\b"                       # per gram of nitrogen
    r"|/\s*100\s*g\s*n\b"
    r"|/\s*60\b"                           # McCance's tryptophan/60
    r"|per\s+100\s*g\s*(?:fa|fatty)\b", re.I)


def is_other_basis(label) -> bool:
    """True when the column is measured against something other than 100 g of food."""
    return bool(label is not None and NON_FOOD_BASIS_RE.search(str(label)))


def _parse(unit: str) -> str | None:
    """Reduce a source's unit string to a bare token: 'mg/100 g' -> 'mg'."""
    if unit is None or (isinstance(unit, float) and pd.isna(unit)):
        return None
    s = str(unit).strip().lower()
    if not s or s in {"nan", "-", "―", "—", "none"}:
        return None
    # "mg/100 g", "mg per 100g", "(mg)", "mg/100 g EP"
    s = s.strip("()[] ")
    parts = re.split(r"\s*(?:/|per)\s*", s)
    if len(parts) > 1:
        # The DENOMINATOR has to be the per-100 g basis everything else uses.
        # AFCD also reports each amino acid "(mg/gN)" - milligrams per gram of
        # NITROGEN - which is a different measurement, not a different unit, and
        # scaling it by 1/1000 would put a number on the wrong basis under the
        # right-looking unit. Refuse, and let the caller mint a separate id.
        den = re.sub(r"[^a-z0-9]", "", parts[1])
        if den not in {"100g", "100gep", "100ml", "100mledible", "100gedible"}:
            _UNKNOWN.add(str(unit))
            return None
    s = re.sub(r"[^a-z%µμ]", "", parts[0])
    if s in _TO_GRAM or s in _NON_MASS:
        return s
    _UNKNOWN.add(str(unit))
    return None


def unknown_units() -> list[str]:
    """Source unit strings this module could not parse, for the ingester to print."""
    return sorted(_UNKNOWN)


def fdc_units(nutrient_csv) -> Dict[int, str]:
    """{nutrient_id: unit_name} from FDC's nutrient.csv."""
    d = pd.read_csv(nutrient_csv, usecols=["id", "unit_name"])
    return {int(i): str(u) for i, u in zip(d["id"], d["unit_name"]) if pd.notna(u)}


def conversion_factor(src_unit, dst_unit) -> float:
    """Multiplier taking a value in `src_unit` to `dst_unit`; 1.0 if unsure.

    Refusing to guess is the point. A wrong factor is worse than no factor,
    because it is just as silent and no longer has an implausible magnitude to
    give it away.
    """
    a, b = _parse(src_unit), _parse(dst_unit)
    if a is None or b is None or a == b:
        return 1.0
    if a in _TO_GRAM and b in _TO_GRAM:
        return _TO_GRAM[a] / _TO_GRAM[b]
    return 1.0


def report(pairs: Iterable[tuple[str, str, str, float]]) -> None:
    """Print the conversions an ingester is about to apply."""
    rows = [p for p in pairs if p[3] != 1.0]
    if not rows:
        print("  -> no unit conversions needed")
        return
    print(f"  -> {len(rows)} columns need a unit conversion:")
    for col, src, dst, f in rows:
        print(f"       {col[:44]:44s} {src:>10s} -> {dst:<4s}  x{f:g}")
