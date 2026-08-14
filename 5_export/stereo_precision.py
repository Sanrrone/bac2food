#!/usr/bin/env python3
"""stereo_precision.py — how often does the name matcher equate two different stereoisomers?

The matcher in ../chebi/dict_to_chebi.py queries several normalized views of a substrate
name, and normalization strips stereochemical descriptors ("D-", "(R)-", "alpha-"). So it can
in principle map a substrate onto a ChEBI entity that differs from it only in configuration.
This quantifies that, because "the matcher might conflate isomers" is the first thing a
reviewer asks and an unquantified answer is worthless.

Three outcomes, and only the third is an error:

  GENERALIZED  the substrate carries a descriptor, the ChEBI entity is the unspecified parent
               ("D-gluconic acid" -> "gluconic acid"). Not an error: composition tables report
               the parent themselves, so this is the resolution the food side actually has.
  SPECIALIZED  the reverse - an unspecified substrate onto a specified entity. Weak, but it
               does not assert the wrong configuration.
  MISMATCH     both are specified and they DISAGREE ("(R)-" onto "(S)-"). The real error.

Reported against the ChEBI-matched rows only, and separately for the rows that reach a food
nutrient, since a mismatch that never reaches a nutrient cannot affect a prediction.

Usage:
    python stereo_precision.py [--export /data/bac2food/exports/enzyme_substrate_chebi.tsv]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

# Stereochemical / configurational descriptors, as they appear in ChEBI and BRENDA names.
# Anchored so "L-" matches the prefix of "L-alanine" but not the "l-" inside "glucosyl-".
DESCRIPTOR = re.compile(
    r"""(?:^|[\s,(\[-])(
        [DL]|d|l                       # D- / L-
        | \(\s*[RS]\s*\)               # (R)- / (S)-
        | \(\s*[+-]\s*\)               # (+)- / (-)-
        | \(\s*\+/-\s*\)               # (+/-)-
        | alpha|beta|α|β               # anomeric / positional
        | cis|trans                    # geometric
        | \(\s*[EZ]\s*\)               # (E)- / (Z)-
        | erythro|threo|rac|meso
        | [RS](?=-)                    # bare R- / S-
    )(?=-|\s|$)""",
    re.VERBOSE,
)

# Strip descriptors and non-alphanumerics: what remains is the configuration-free skeleton.
STRIP = re.compile(r"[^a-z0-9]+")


def descriptors(name: str) -> frozenset[str]:
    return frozenset(m.group(1).lower().replace(" ", "") for m in DESCRIPTOR.finditer(name))


def skeleton(name: str) -> str:
    s = DESCRIPTOR.sub(" ", name.lower())
    return STRIP.sub("", s)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", default="/data/bac2food/exports/enzyme_substrate_chebi.tsv")
    args = ap.parse_args()

    matched = generalized = specialized = mismatch = 0
    mismatch_in_model = 0
    examples: list[tuple[str, str]] = []

    with open(args.export, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cid = (row.get("chebi_id") or "").strip()
            if not cid:
                continue
            matched += 1
            sub = (row.get("substrate") or "").strip()
            lbl = (row.get("chebi_name") or "").strip()
            if not sub or not lbl or sub.lower() == lbl.lower():
                continue
            # Only interested in pairs that are the SAME compound modulo configuration.
            if skeleton(sub) != skeleton(lbl):
                continue
            ds, dl = descriptors(sub), descriptors(lbl)
            if ds == dl:
                continue
            if ds and not dl:
                generalized += 1
            elif dl and not ds:
                specialized += 1
            else:
                mismatch += 1
                if (row.get("in_model") or "").strip().lower() == "yes":
                    mismatch_in_model += 1
                if len(examples) < 8:
                    examples.append((sub, lbl))

    differing = generalized + specialized + mismatch
    pct = lambda n: 100.0 * n / matched if matched else 0.0
    print(f"  ChEBI-matched rows                       {matched:>8,}")
    print(f"  differ ONLY in stereochemistry           {differing:>8,}  ({pct(differing):.1f}%)")
    print(f"    generalized to unspecified parent      {generalized:>8,}   (not an error)")
    print(f"    specialized to a specified isomer      {specialized:>8,}   (weak, not wrong)")
    print(f"    MISMATCH, specified vs specified       {mismatch:>8,}  ({pct(mismatch):.1f}%)")
    print(f"      of those, reaching a food nutrient   {mismatch_in_model:>8,}")
    if examples:
        print("\n  mismatch examples (substrate -> ChEBI entity):")
        for s, l in examples:
            print(f"    {s}  ->  {l}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
