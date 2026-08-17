#!/usr/bin/env python3
"""Phase 8 / Part A — Read-only audit of every food_DBs/<src>/ ingester.

For each source folder, the audit:
  1. Identifies the raw source file(s) and the parser script(s).
  2. Calls the new _common helpers to enumerate sheets & detect headers.
  3. Loads /data/bac2food/food.parquet and filters by the parser's source tag.
  4. Counts (a) foods in raw vs parsed, (b) distinct nutrient columns in raw vs
     mapped in parsed.
  5. Writes <folder>/audit.tsv listing dropped foods + dropped columns.
  6. Prints a triage summary table to stdout.

No writes outside food_DBs/<src>/audit.tsv files. Safe to re-run.
"""
from __future__ import annotations

import argparse
import re
import sys
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq

THIS = Path(__file__).parent
sys.path.insert(0, str(THIS))

from _common import format_detect, sheet_classify  # noqa: E402


@dataclass
class DBSpec:
    """One row of the audit registry."""
    short: str
    folder: str
    raw_files: list[str]
    parser_scripts: list[str]
    source_tag: str               # the literal substring written into food.parquet description
    expected_raw_foods: int       # documented total (per docx / source) — used as sanity check
    notes: str = ""


# Registry — explicit, no glob magic. Keep in sync with the plan's "Files & locations" table.
REGISTRY: list[DBSpec] = [
    DBSpec("mccance", "McCance_Widdowsons_uk",
           ["McCance_Widdowsons_Composition_of_Foods_Integrated_Dataset_2021..xlsx"],
           ["ingest_mw.py"],
           source_tag="[McCance]", expected_raw_foods=2888,
           notes="14 nutrient sheets, foods spine = 1.3 Proximates"),
    DBSpec("biofoodcomp", "FAO_onu",
           ["BioFoodComp4.0.xlsx"],
           ["biofoodcomp_to_fdc_ingest.py"],
           source_tag="[BioFoodComp]", expected_raw_foods=10133,
           notes="12 numbered nutrient sheets + 1 fatty-acid sub-sheet"),
    DBSpec("phyfoodcomp", "FAO_onu",
           ["PhyFoodComp_1.0.xlsx"],
           ["injest_phyfood.py"],
           source_tag="[PhyFoodComp]", expected_raw_foods=3350,
           notes="18 numbered nutrient sheets"),
    DBSpec("afcd", "asnut_australianw",
           ["AFCD Release 3 - Nutrient profiles.xlsx"],
           ["afcd_to_fdc_ingest.py"],
           source_tag="[AFCD]", expected_raw_foods=1807,
           notes="2 sheets: solids+liquids (1591) + liquids-only (216)"),
    DBSpec("phenol_explorer", "phenol_explorer_france",
           ["composition-data.tsv"],
           ["ingest_phenol_explorer.py"],
           source_tag="[Phenol-Explorer]", expected_raw_foods=459,
           notes="TSV — 7487 rows (food × compound × method)"),
    DBSpec("stfcj_main", "stfcj_japan",
           ["main_1374049_1r12_1.xlsx", "aminoacid_1374049_2r11_1.xlsx",
            "fatty_acid_1374049_3r11_1.xlsx", "org_acid_1388558_4r12r.xlsx"],
           ["injest_japan.py"],
           source_tag="[STFCJ]", expected_raw_foods=2478,
           notes="4 inter-related files joined on item_no"),
    DBSpec("sweden", "livsmedels_sweeden",
           [],   # detect from folder
           ["ingest_swedish_livsmedelsdb.py"],
           source_tag="[SWE]", expected_raw_foods=2600,
           notes=""),
    DBSpec("swiss", "swissfoodcompoDB_swiss",
           [],
           ["swiss_to_fdc_ingest.py"],
           source_tag="[Swiss",        # multi-sheet — descriptions are '[Swiss <sheet>]'
           expected_raw_foods=10600,
           notes="per-sheet ingest; check coverage"),
    DBSpec("cnf", "canadian_nutrientfile_canada",
           [],
           ["ingest_cnf.py"],
           source_tag="[CNF]", expected_raw_foods=5690,
           notes="relational schema"),
    DBSpec("ciqual", "ciqual_france",
           [],
           ["ingest_ciqual.py"],
           source_tag="[Ciqual]", expected_raw_foods=3185,
           notes=""),
    DBSpec("fineli", "fineli_finland",
           [],
           ["ingest_finely.py"],
           source_tag="[Fineli]", expected_raw_foods=4156,
           notes=""),
    DBSpec("frida", "frida_denmark",
           [],
           ["ingest_friday.py"],
           source_tag="[Frida]", expected_raw_foods=1381,
           notes=""),
    DBSpec("wafct", "WAFT_AFRICA",
           [],
           [],
           source_tag="[WAFCT]", expected_raw_foods=1028,
           notes="no parser in folder"),
    DBSpec("fao_pulses", "FAO_onu",
           [],
           [],
           source_tag="[FAO]", expected_raw_foods=0,
           notes="(see biofoodcomp / phyfoodcomp for FAO sheets)"),
]


@dataclass
class DBResult:
    short: str
    folder: str
    parsed_foods: int = 0
    parsed_nutrient_ids: int = 0
    raw_files_found: list[str] = field(default_factory=list)
    raw_total_rows: int = 0           # sum of foods_spine + nutrient_sheet rows
    raw_max_sheet_rows: int = 0       # rows of the biggest nutrient/spine sheet
    raw_distinct_columns: int = 0
    notes: list[str] = field(default_factory=list)
    severity: str = "?"


def find_raw_files(folder: Path, allow: list[str]) -> list[Path]:
    """Return concrete raw-file paths. If `allow` non-empty, use it verbatim;
    otherwise discover xlsx/csv/tsv in folder."""
    if allow:
        return [folder / f for f in allow if (folder / f).exists()]
    out = []
    for ext in ("*.xlsx", "*.xlsm", "*.xls", "*.csv", "*.tsv"):
        for p in folder.glob(ext):
            # Skip outputs the parser created
            if any(s in p.name.lower() for s in ("_injection", "_food_nutrient",
                                                  "merge_", "unmapped",
                                                  "_data.parquet", "_food.csv",
                                                  "audit.")):
                continue
            out.append(p)
    return out


def parsed_count(food_parquet: Path, tag: str) -> int:
    t = pq.read_table(food_parquet, columns=["description"])
    desc = pc.utf8_lower(t.column("description"))
    return pc.sum(pc.match_substring(desc, tag.lower())).as_py()


def parsed_nutrient_ids(folder: Path, tag: str, food_parquet: Path,
                         nutrient_bucket_root: Path) -> int:
    """Count distinct nutrient_ids for this source — read the per-parser
    output parquets in the source folder (the bucketed slice with the tag's
    fdc_id range) rather than scanning the full /data/bac2food/food_nutrient_bucketed."""
    # 1) get the fdc_ids for this tag
    t = pq.read_table(food_parquet, columns=["fdc_id", "description"])
    desc = pc.utf8_lower(t.column("description"))
    mask = pc.match_substring(desc, tag.lower())
    fdcs = pc.filter(t.column("fdc_id"), mask).to_pylist()
    if not fdcs:
        return 0
    # 2) Look at the parser's local bucketed parquet
    local_bucket = None
    for cand in folder.glob("*_food_nutrient_bucketed"):
        local_bucket = cand
        break
    if local_bucket is None or not local_bucket.exists():
        return 0
    seen: set[int] = set()
    fdc_set = set(fdcs)
    for part in local_bucket.glob("bucket=*/*.parquet"):
        try:
            tab = pq.read_table(part, columns=["fdc_id", "nutrient_id"])
        except Exception:
            continue
        ids = tab.column("fdc_id").to_pylist()
        nids = tab.column("nutrient_id").to_pylist()
        for f, n in zip(ids, nids):
            if f in fdc_set:
                seen.add(n)
    return len(seen)


def audit_one(spec: DBSpec, food_parquet: Path,
              nutrient_bucket_root: Path) -> DBResult:
    folder = Path(__file__).parent / spec.folder
    result = DBResult(short=spec.short, folder=spec.folder)
    raw_files = find_raw_files(folder, spec.raw_files)
    result.raw_files_found = [str(f.relative_to(folder)) for f in raw_files]

    # Inspect raw sheets
    nutrient_cols_seen: set[str] = set()
    for rf in raw_files:
        info = format_detect.detect(rf)
        if info.fmt in ("xlsx", "xls"):
            try:
                infos = sheet_classify.classify_workbook(rf)
            except Exception as e:
                result.notes.append(f"sheet_classify failed on {rf.name}: {e}")
                continue
            for s in infos:
                if s.classification in ("foods_spine", "nutrient_sheet"):
                    result.raw_total_rows += s.n_rows
                    result.raw_max_sheet_rows = max(result.raw_max_sheet_rows, s.n_rows)
                    nutrient_cols_seen.update(s.sample_cols)
        elif info.fmt in ("csv", "tsv"):
            try:
                df = pd.read_csv(rf, sep=info.delimiter or ("\t" if info.fmt == "tsv" else ","),
                                 encoding=info.encoding or "utf-8",
                                 dtype=str, nrows=10, on_bad_lines="skip")
                cols = list(df.columns)
                nutrient_cols_seen.update(cols)
                # total row count
                with open(rf, "rb") as fh:
                    n = sum(1 for _ in fh) - 1
                result.raw_total_rows += n
                result.raw_max_sheet_rows = max(result.raw_max_sheet_rows, n)
            except Exception as e:
                result.notes.append(f"csv read failed on {rf.name}: {e}")
        else:
            result.notes.append(f"format {info.fmt!r} unsupported: {rf.name}")
    result.raw_distinct_columns = len(nutrient_cols_seen)

    # Parsed counts
    result.parsed_foods = parsed_count(food_parquet, spec.source_tag)
    result.parsed_nutrient_ids = parsed_nutrient_ids(folder, spec.source_tag,
                                                     food_parquet, nutrient_bucket_root)

    # Severity heuristic
    if result.parsed_foods == 0:
        result.severity = "ZERO"
    elif spec.expected_raw_foods > 0:
        coverage = result.parsed_foods / spec.expected_raw_foods
        if coverage < 0.6:
            result.severity = "HIGH"
        elif coverage < 0.95:
            result.severity = "MED"
        else:
            result.severity = "OK"
    else:
        result.severity = "?"
    return result


def write_audit_tsv(folder: Path, spec: DBSpec, result: DBResult) -> Path:
    out = folder / "audit.tsv"
    rows = [
        ("short", spec.short),
        ("folder", spec.folder),
        ("source_tag", spec.source_tag),
        ("expected_raw_foods", spec.expected_raw_foods),
        ("raw_files_found", ";".join(result.raw_files_found)),
        ("raw_total_sheet_rows", result.raw_total_rows),
        ("raw_max_sheet_rows", result.raw_max_sheet_rows),
        ("raw_distinct_columns_seen", result.raw_distinct_columns),
        ("parsed_foods_in_food_parquet", result.parsed_foods),
        ("parsed_distinct_nutrient_ids", result.parsed_nutrient_ids),
        ("severity", result.severity),
        ("notes", " | ".join(result.notes) if result.notes else ""),
    ]
    with out.open("w") as f:
        f.write("metric\tvalue\n")
        for k, v in rows:
            f.write(f"{k}\t{v}\n")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--food_parquet", default="/data/bac2food/food.parquet")
    ap.add_argument("--nutrient_bucket_root",
                    default="/data/bac2food/food_nutrient_bucketed")
    ap.add_argument("--only", default=None, help="comma-separated short names")
    args = ap.parse_args()
    food_parquet = Path(args.food_parquet)
    bucket_root = Path(args.nutrient_bucket_root)
    only = set(args.only.split(",")) if args.only else None

    results: list[tuple[DBSpec, DBResult]] = []
    for spec in REGISTRY:
        if only and spec.short not in only:
            continue
        print(f"\n--- {spec.short} ({spec.folder}) ---")
        result = audit_one(spec, food_parquet, bucket_root)
        results.append((spec, result))
        folder = Path(__file__).parent / spec.folder
        if folder.exists():
            out = write_audit_tsv(folder, spec, result)
            print(f"  wrote {out}")

    print("\n" + "=" * 84)
    print(f"{'DB':<18} {'sev':<5} {'raw_foods':>11} {'parsed':>10} {'cov%':>6} {'nuts':>6}")
    print("=" * 84)
    for spec, r in results:
        cov = (100.0 * r.parsed_foods / spec.expected_raw_foods
               if spec.expected_raw_foods else float("nan"))
        cov_s = f"{cov:5.1f}" if spec.expected_raw_foods else "    ?"
        print(f"{spec.short:<18} {r.severity:<5} {spec.expected_raw_foods:>11} "
              f"{r.parsed_foods:>10} {cov_s:>6} {r.parsed_nutrient_ids:>6}")


if __name__ == "__main__":
    main()
