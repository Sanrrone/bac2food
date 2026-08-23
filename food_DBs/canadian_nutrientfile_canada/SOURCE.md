# Source data for `canadian_nutrientfile_canada`

The raw release is **not** redistributed with this repository — its terms are the
provider's, not ours. Download it from the page below and drop it in this folder,
then run the ingest script here. Rights per source are in
`5_export/licence_tiers.csv` (= Supplementary Table S1 of the Data Descriptor).

Generated from `food_DBs/SOURCES.tsv` by `food_DBs/write_source_stubs.py` — edit
the table, not this file.

## Canadian Nutrient File (`cnf`)

* **Provider** — Health Canada
* **Version this build ingested** — 2015
* **Download** — <https://www.canada.ca/en/health-canada/services/food-nutrition/healthy-eating/nutrient-data/canadian-nutrient-file-2015-download-files.html>
* **Files expected here** — FOOD NAME.csv, FOOD GROUP.csv, FOOD SOURCE.csv, NUTRIENT AMOUNT.csv, NUTRIENT NAME.csv, NUTRIENT SOURCE.csv, MEASURE NAME.csv, CONVERSION FACTOR.csv
* **Note** — The relational CSV set. A 2026 release exists on open.canada.ca; this build ingested 2015.
