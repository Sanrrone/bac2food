#!/usr/bin/env python3
"""1.0b_brenda_json_to_digest.py — EC -> substrate digest from BRENDA's official JSON download.

Replaces 1.0_eggnog_ec_substrates_parser.py, which parsed a DSMZ SPARQL scrape
(ec_species_substrate.tar.xz). Measured 2026-08-06, that scrape was an incomplete slice of
the very same release: substrates for 5,263 EC where the download carries 6,901, and only
63.5% of the download's (ec, substrate) pairs. The endpoint is also unreachable now, so the
download is the only way to refresh this layer. Same output contract as 1.0_: a two-column
`ec  substrate` TSV, ready for 1.5_reactions_to_digest.py.

Get the input from brenda-enzymes.org/download.php (registration, free). Prefer the JSON to
the flat TXT: they carry identical content -- 8,129 EC blocks each -- but the TXT inlines
protein refs (#10#), citations (<51>), commentaries (...) and continuation lines into the
substrate string, and every one of those is a chance to corrupt a name silently. The JSON
hands the same record over already decomposed, and ijson streams it in ~7 s / 65 MB, so the
709 MB file never lands in memory.

BRENDA data is CC BY 4.0 (verified 2026-08-06 at brenda-enzymes.org/license.php) -- it may be
redistributed, but it must be cited.

Usage:
    python 1.0b_brenda_json_to_digest.py --json /data/bac2food/brenda_2026_1/brenda_2026_1.json \
                                         --out 2_digest_dict.tsv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

import ijson

# BRENDA appends reversibility info to a reaction value in braces: "A + B = C + D {r}".
# The README documents it for SP/NSP. Strip it BEFORE splitting, or it stays glued to the
# last compound and yields phantom substrates -- "nad+ {r}" (252 rows), "h2o {r}" (235),
# "? {r}" (313) -- which then sail past the ubiquitous-cofactor filter in 1.5_ because that
# list quite reasonably contains "nad+", not "nad+ {r}".
REVERSIBILITY = re.compile(r"\s*\{[^{}]*\}\s*$")

# The kinetic fields name their substrate in braces instead: "0.05 {benzyl alcohol}".
BRACE = re.compile(r"\{([^}]*)\}")

# Leading stoichiometric coefficient: "2 NAD+" -> "NAD+".
STOICH = re.compile(r"^\d+\s+")

# BRENDA's JSON values carry a trailing carriage return -- every single one. Left in, it
# becomes part of the substrate name, and every exact-match check downstream silently fails:
# "NAD(P)+\r" is not "NAD(P)+", so the cofactor filter in 1.5_ cannot catch it and
# dict_to_chebi.py cannot resolve it. 1.0_ stripped control characters from the SPARQL
# scrape for the same reason.
CONTROL = re.compile(r"[\x00-\x1f\x7f]")


def clean(part: str) -> str:
    return CONTROL.sub("", part).strip()

# Placeholders and generic acceptors that are not compounds. "more" is BRENDA's own
# placeholder for "further substrates known but not listed".
PLACEHOLDER = {"", "?", "more", "acceptor", "reduced acceptor", "donor", "reduced donor"}

REACTION_FIELDS = ("substrates_products", "natural_substrates_products")
KINETIC_FIELDS = ("km_value", "turnover_number", "kcat_km_value")


def compounds(value: str):
    """Yield every compound in a BRENDA reaction equation, both sides.

    BOTH sides on purpose. BRENDA reactions are largely reversible and the resource asks
    "can this enzyme act on this compound", not "in which direction". Taking only the left
    side loses 36,663 pairs the old scrape had -- e.g. EC 1.1.1.2 acting on
    3-methoxybenzaldehyde, which is written as a product.
    """
    value = REVERSIBILITY.sub("", clean(value))
    for side in value.split("="):
        for part in side.split(" + "):
            part = STOICH.sub("", clean(part)).strip()
            if part.lower() not in PLACEHOLDER and 0 < len(part) < 200:
                yield part


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", required=True, help="brenda_<release>.json (uncompressed)")
    ap.add_argument("--out", default="2_digest_dict.tsv")
    ap.add_argument("--no-kinetic", action="store_true",
                    help="skip the KM/turnover/kcat-KM {substrate} braces; reaction equations only")
    args = ap.parse_args()

    pairs: set[tuple[str, str]] = set()
    n_ec = n_react = n_kin = 0
    release = None

    with open(args.json, "rb") as fh:
        # Two passes would mean reading 709 MB twice; the release string sits before "data".
        for prefix, event, value in ijson.parse(fh):
            if prefix == "release":
                release = value
                break
    with open(args.json, "rb") as fh:
        for ec, payload in ijson.kvitems(fh, "data"):
            if ec == "spontaneous":
                continue
            n_ec += 1
            for field in REACTION_FIELDS:
                for entry in payload.get(field) or []:
                    v = entry.get("value") or ""
                    if "=" not in v:
                        continue
                    n_react += 1
                    for c in compounds(v):
                        pairs.add((ec, c))
            if args.no_kinetic:
                continue
            for field in KINETIC_FIELDS:
                for entry in payload.get(field) or []:
                    for m in BRACE.findall(entry.get("value") or ""):
                        m = clean(m)
                        if m.lower() not in PLACEHOLDER and 0 < len(m) < 200:
                            n_kin += 1
                            pairs.add((ec, m))

    with open(args.out, "w", newline="", encoding="utf-8") as out:
        # lineterminator="\n" is NOT cosmetic. csv.writer defaults to the excel dialect's
        # "\r\n", which parks a carriage return on the last field of every row: the substrate
        # becomes "NAD(P)+\r", so the cofactor filter in 1.5_ never matches it and every
        # exact-name lookup downstream misses. Unix line endings throughout.
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(["ec", "substrate"])
        w.writerows(sorted(pairs))

    print(f"[*] BRENDA release {release}", flush=True)
    print(f"[*] {n_ec:,} EC entries | {n_react:,} reaction records | {n_kin:,} kinetic braces",
          flush=True)
    print(f"[*] wrote {args.out}: {len(pairs):,} unique (ec, substrate) pairs, "
          f"{len({e for e, _ in pairs}):,} EC", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
