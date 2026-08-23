#!/usr/bin/env python3
"""write_source_stubs.py — generate food_DBs/<folder>/SOURCE.md from SOURCES.tsv.

The raw releases are not redistributed here, so each source folder ships a pointer to the
provider instead. The pointers are GENERATED rather than hand-written: fifteen hand-kept
copies of the same facts is fifteen chances to go stale, which is exactly how Supplementary
Table S1 came to describe NEVO under terms the deposit contradicted.

Run after editing SOURCES.tsv:

    python3 food_DBs/write_source_stubs.py

Two folders take more than one source (FAO_onu holds BioFoodComp and PhyFoodComp), so the
stub is written per FOLDER with a section per source_db.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
TSV = HERE / "SOURCES.tsv"


def rows():
    with TSV.open(encoding="utf-8", newline="") as fh:
        body = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(body, delimiter="\t"))


def main() -> None:
    by_folder = defaultdict(list)
    for r in rows():
        if r["folder"] and r["folder"] != "-":
            by_folder[r["folder"]].append(r)

    written = 0
    for folder, rs in sorted(by_folder.items()):
        d = HERE / folder
        if not d.is_dir():
            raise SystemExit(f"SOURCES.tsv names a folder that does not exist: {folder}")
        out = [f"# Source data for `{folder}`", "",
               "The raw release is **not** redistributed with this repository — its terms are the",
               "provider's, not ours. Download it from the page below and drop it in this folder,",
               "then run the ingest script here. Rights per source are in",
               "`5_export/licence_tiers.csv` (= Supplementary Table S1 of the Data Descriptor).", "",
               "Generated from `food_DBs/SOURCES.tsv` by `food_DBs/write_source_stubs.py` — edit",
               "the table, not this file.", ""]
        for r in rs:
            out += [f"## {r['database']} (`{r['source_db']}`)", "",
                    f"* **Provider** — {r['provider']}",
                    f"* **Version this build ingested** — {r['version_used']}",
                    f"* **Download** — <{r['download_page']}>",
                    f"* **Files expected here** — {r['expected_files']}"]
            if r["note"]:
                out += [f"* **Note** — {r['note']}"]
            out += [""]
        (d / "SOURCE.md").write_text("\n".join(out), encoding="utf-8")
        written += 1
    print(f"[*] wrote {written} SOURCE.md stubs from {TSV.name}")


if __name__ == "__main__":
    main()
