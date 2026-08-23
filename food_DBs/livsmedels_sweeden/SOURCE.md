# Source data for `livsmedels_sweeden`

The raw release is **not** redistributed with this repository — its terms are the
provider's, not ours. Download it from the page below and drop it in this folder,
then run the ingest script here. Rights per source are in
`5_export/licence_tiers.csv` (= Supplementary Table S1 of the Data Descriptor).

Generated from `food_DBs/SOURCES.tsv` by `food_DBs/write_source_stubs.py` — edit
the table, not this file.

## Livsmedelsdatabasen (`swedish`)

* **Provider** — Livsmedelsverket (Sweden)
* **Version this build ingested** — 2024
* **Download** — <https://www.livsmedelsverket.se/en/about-us/open-data/food-composition-data/>
* **Files expected here** — LivsmedelsDB_202603061604.xlsx
* **Note** — The export carries the source's editing layer alongside the reference set; the ingest drops the copy records (see food_DBs/_common/non_nutrients.py).
