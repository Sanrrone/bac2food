# Source data for `fooddatacentral_usa`

The raw release is **not** redistributed with this repository — its terms are the
provider's, not ours. Download it from the page below and drop it in this folder,
then run the ingest script here. Rights per source are in
`5_export/licence_tiers.csv` (= Supplementary Table S1 of the Data Descriptor).

Generated from `food_DBs/SOURCES.tsv` by `food_DBs/write_source_stubs.py` — edit
the table, not this file.

## FoodData Central (`fdc`)

* **Provider** — USDA ARS (USA)
* **Version this build ingested** — 2025
* **Download** — <https://fdc.nal.usda.gov/download-datasets/>
* **Files expected here** — food.csv, food_nutrient.csv, nutrient.csv, food_category.csv
* **Note** — Take the "Full Download" CSV bundle. Only nutrient.csv is kept in this folder; the rest are read from the raw drop at /data/bac2food/fdc_raw.
