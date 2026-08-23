# Source data for `NEVO_netherlands`

The raw release is **not** redistributed with this repository — its terms are the
provider's, not ours. Download it from the page below and drop it in this folder,
then run the ingest script here. Rights per source are in
`5_export/licence_tiers.csv` (= Supplementary Table S1 of the Data Descriptor).

Generated from `food_DBs/SOURCES.tsv` by `food_DBs/write_source_stubs.py` — edit
the table, not this file.

## NEVO-online (`nevo`)

* **Provider** — RIVM (Netherlands)
* **Version this build ingested** — 2025/9.0
* **Download** — <https://www.rivm.nl/en/dutch-food-composition-database/use-of-nevo-online/request-dataset>
* **Files expected here** — NEVO2025_v9.0.xlsx
* **Note** — Request the dataset from RIVM; it is free. The release itself is in neither this repository nor the deposit — RIVM asked that users be sent to their download page. 5_export/reconstruct_nevo.py rebuilds the NEVO partition from it.
