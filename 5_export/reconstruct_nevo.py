#!/usr/bin/env python3
"""reconstruct_nevo.py — rebuild the NEVO partition of the composition table locally.

NEVO's terms permit use "only unchanged and stating the source and version number", and state
that while a user "is entitled to make additions" they are "not entitled to make amendment" to
the dataset. Re-keying every food to a shared `nutrient_id` and normalizing amounts to per 100 g
edible portion is an amendment on any reasonable reading, so the deposited `food_nutrients.tsv`
carries no NEVO-derived values. Repackaging them would not help: a separate file, a fork or a
"partition" all redistribute the same amended values.

The deposit ships no NEVO file at all — not even the original release unchanged, which those
terms do appear to permit. Confirmation of that reading has been requested from RIVM and not yet
received, and shipping first and asking afterwards is the wrong order for a source whose terms
are this explicit. What the deposit ships is this script and the mapping it applies. A user who
holds their own copy of NEVO reconstructs the missing rows locally in one command, and the
reconstruction is an addition to their own copy rather than a redistribution of an amended one.

Request the dataset from RIVM — free, and the release this mapping was built against is 2025/9.0:

    https://www.rivm.nl/en/dutch-food-composition-database/use-of-nevo-online/request-dataset

Then:

    python reconstruct_nevo.py --nevo_xlsx NEVO2025_v9.0.xlsx --out food_nutrients_nevo.tsv
    cat food_nutrients_nevo.tsv >> food_nutrients.tsv     # header-skipped; see --append

A later NEVO release will still reconstruct, but the column names it maps may have moved; the
script fails loudly rather than silently dropping components if none of them is recognised.

The output carries the identical 13-column schema, so the reconstructed rows are
indistinguishable from the rest of the table once appended. The mapping is imported from
`food_DBs/NEVO_netherlands/ingest_nevo.py` rather than restated here, so a local reconstruction
cannot drift from what the build pipeline produces.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
INGEST = REPO / "food_DBs/NEVO_netherlands/ingest_nevo.py"

sys.path.insert(0, str(REPO / "food_DBs"))
import fdc_blocks  # noqa: E402
SCHEMA = ["fdc_id", "description", "data_type", "food_category",
          "nutrient_id", "nutrient_name", "unit_name", "amount", "source_db",
          "source_version", "license_id", "redistribution_flag", "attribution_string"]


def load_maps():
    """Import CAT_MAP / NUTRIENT_MAP from the ingest script.

    fdc_ids no longer come from an offset defined here. They are accessions held in
    fdc_id_map.tsv, so a rebuilt NEVO layer lands on exactly the ids the rest of the
    resource already references -- which an offset could not guarantee once NEVO
    republished with renumbered codes.
    """
    if not INGEST.exists():
        sys.exit(f"[!] {INGEST} not found; it defines the mapping this script applies")
    spec = importlib.util.spec_from_file_location("ingest_nevo", INGEST)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.CAT_MAP, mod.NUTRIENT_MAP, mod.DEFAULT_CATEGORY


# NEVO is not a source of the deposited resource, so it has no row in the published licence
# table. Its rights therefore live here, in the one script that touches it. Reconstructed rows
# still carry `redistribution_flag = no`, which is the point: a user who appends them can filter
# them back out before sharing anything derived from the table.
NEVO_RIGHTS = {
    "version": "2025/9.0",
    "licence": "Use permitted only unchanged, stating source and version; "
               "user is not entitled to make amendment",
    "attribution_string": "NEVO-online version 2025/9.0, RIVM, Bilthoven, The Netherlands.",
}


def licence_row(path: Path) -> dict:
    """NEVO's rights. Prefers a `nevo` row in a user-supplied licence table — someone who holds
    a different release should say so there — and falls back to the constants above."""
    if path.exists():
        with path.open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh):
                if r.get("source_db") == "nevo":
                    print(f"[*] Using the 'nevo' row from {path} ({r['version']})")
                    return r
    print(f"[*] No 'nevo' row in {path}; using built-in rights for {NEVO_RIGHTS['version']}")
    return NEVO_RIGHTS


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nevo_xlsx", required=True,
                    help="your own copy of the NEVO release; request it free from RIVM at "
                         "https://www.rivm.nl/en/dutch-food-composition-database/use-of-nevo-online"
                         "/request-dataset (built against 2025/9.0)")
    ap.add_argument("--out", default="food_nutrients_nevo.tsv")
    ap.add_argument("--nutrient", default="/data/bac2food/nutrient.csv",
                    help="FDC nutrient catalogue, for nutrient_name and unit_name")
    ap.add_argument("--licence_table",
                    default=str(Path.home() / "Desktop/paper4/Supplementary_Table_S1.csv"))
    ap.add_argument("--append", default=None,
                    help="append directly to this food_nutrients.tsv instead of writing --out")
    args = ap.parse_args()

    CAT_MAP, NUTRIENT_MAP, DEFAULT_CAT = load_maps()
    lic = licence_row(Path(args.licence_table))

    df = pd.read_excel(args.nevo_xlsx, sheet_name=0)
    df["fdc_id"] = fdc_blocks.assign("nevo", df["NEVO-code"].astype(int))
    df["food_category"] = df["Food group"].map(CAT_MAP).fillna(DEFAULT_CAT)

    cols = [c for c in NUTRIENT_MAP if c in df.columns]
    if not cols:
        sys.exit("[!] none of the expected NEVO component columns are present; "
                 "is this the composition sheet of the NEVO release?")
    for c in cols:
        if df[c].dtype == object:
            df[c] = (df[c].astype(str)
                     .str.replace(r"^[a-zA-Z<>].*", "0.0", regex=True)
                     .str.replace(",", ".", regex=False))
            df[c] = pd.to_numeric(df[c], errors="coerce")

    long = df.melt(id_vars=["fdc_id", "Engelse naam/Food name", "food_category"],
                   value_vars=cols, var_name="nevo_nutrient", value_name="amount")
    long = long.dropna(subset=["amount"])
    long = long[long["amount"] > 0]
    long["nutrient_id"] = long["nevo_nutrient"].map(NUTRIENT_MAP).astype(int)

    nut = pd.read_csv(args.nutrient, usecols=["id", "name", "unit_name"]).rename(
        columns={"id": "nutrient_id", "name": "nutrient_name"})
    out = long.merge(nut, on="nutrient_id", how="left")
    out["description"] = out["Engelse naam/Food name"].astype(str) + " [NEVO]"
    out["data_type"] = "foundation_food"
    out["source_db"] = "nevo"
    out["source_version"] = lic["version"]
    out["license_id"] = lic["licence"]
    out["redistribution_flag"] = "no"
    out["attribution_string"] = lic["attribution_string"]
    out = out[SCHEMA].sort_values(["fdc_id", "nutrient_id"], kind="mergesort")

    if args.append:
        target = Path(args.append)
        with target.open(encoding="utf-8", newline="") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        if header != SCHEMA:
            sys.exit(f"[!] {target} has a different schema; refusing to append\n"
                     f"    got      {header}\n    expected {SCHEMA}")
        with target.open("a", encoding="utf-8", newline="") as fh:
            out.to_csv(fh, sep="\t", index=False, header=False, na_rep="",
                       quoting=csv.QUOTE_MINIMAL)
        print(f"[*] Appended {len(out):,} NEVO rows ({out['fdc_id'].nunique():,} foods) to {target}")
    else:
        out.to_csv(args.out, sep="\t", index=False, na_rep="", quoting=csv.QUOTE_MINIMAL)
        print(f"[*] Wrote {args.out}: {len(out):,} rows, {out['fdc_id'].nunique():,} foods")

    print("[*] These rows carry redistribution_flag = no. They are for local use; filter them "
          "out before sharing any derivative of this table.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
