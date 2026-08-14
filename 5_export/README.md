# 5_export — flat reference resources

`export_resources.py` flattens the three reference layers the predictor stands on into plain
TSVs, so they can be read, shared, or deposited without touching the pipeline:

| file | one row is | rows | size |
|---|---|---|---|
| `enzyme_substrate_chebi.tsv` | an enzyme (EC) and one substrate it acts on, with the substrate's ChEBI id | 218,378 | 23 MB |
| `species_enzymes.tsv` | an organism (species + strain) and one EC number it carries | 9,632,315 | 834 MB |
| `food_nutrients.tsv` | a food and one nutrient it contains, with the amount and its source DB | 2,308,817 | 287 MB |

The script only joins and flattens tables the pipeline already built — it recomputes nothing,
so these files are exactly what the predictor sees.

## Run

```bash
python export_resources.py --out_dir /data/bac2food/exports     # all three (~2.1 GB, ~4 min)
python export_resources.py --out_dir /data/bac2food/exports --only enzymes --only species
```

`--out_dir` is **required** and has no default. The output filenames do not encode the input
version — a v6 and a v7 build are both `species_enzymes.tsv` — so a default pointing at the
deposit meant a legacy rebuild silently overwrote the shipped layer. Write rebuilds you do not
intend to ship to a scratch directory instead.

Defaults point at `/data/bac2food`; every input is overridable (`--bact_ec`, `--food`,
`--food_category`, `--nutrient`, `--bucketed_dir`). Peak memory is ~3.7 GB.

## `enzyme_substrate_chebi.tsv`

EC → substrate from the eggNOG v6 digest (`eggnog/2_digest_norm.tsv`), with the substrate resolved
against the ChEBI ontology (`chebi/digest_to_chebi.tsv`).

`ec_number`, `substrate`, `substrate_normalized`, `chebi_id`, `chebi_name`, `chebi_match_type`,
`in_model`, `nutrient_ids`, `nutrient_names`, `model_relation`, `model_score`

* 6,900 EC numbers; 124,538 of the 218,378 rows carry a ChEBI id. The rest are substrates ChEBI
  has no term for — macromolecules, peptides, assay markers, metal clusters. `chebi_match_type`
  says which, and how the match was made (`string_match(100)`, `fuzzy_match(60)`,
  `rhea_confirmed(200)`, `no_match`, …), so rows can be filtered by confidence.
* `in_model` = `yes` on the 5,883 rows that actually reach an FDC nutrient and therefore drive
  food scores (via `0_building/3_nutrient_to_ec.tsv`). `nutrient_ids` / `nutrient_names` list the
  nutrients reached; `model_relation` (`exact`, `is_a`, `conjugate`, …) and `model_score` (0–100)
  describe the best ChEBI link behind that mapping. The other rows are real enzyme–substrate
  facts that no food nutrient measures.

## `species_enzymes.tsv`

`tax_id`, `genus`, `species`, `strain`, `organism`, `ec_number`

* Built from **eggNOG v7** by `eggnog/6.1_eggnog7_species_enzymes.py`. v7 dropped the direct
  EC annotation v6 carried, so the provenance is organism → KEGG KO → EC, bridged through the
  KEGG orthology catalogue (`eggnog/6.0_kegg_ko_to_ec.py`). Only fully specified four-level EC
  numbers are kept — partial ones (`1.1.1.-`) cannot join the substrate digest.
* The legacy v6 route (`/data/bac2food/bact_ec.tsv`, headerless, 58.2M rows / 3.2 GB, exported
  as 9,632,315 deduplicated rows) is **retired** but the source file is retained, so a v6 build
  stays reproducible: `--only species --bact_ec /data/bac2food/bact_ec.tsv`. A `.parquet` with
  the same content is accepted via `--bact_ec` and loads faster.
* The old DSMZ BRENDA SPARQL scrape (`3_normalized_species.tsv`) is **not** used and has been
  deleted: it returned incomplete results, under-reporting the enzymes of every organism. The
  `brenda/` folder is now `eggnog/`, and holds only the EC → substrate digest.
* 20,557,730 rows: 10,751 organisms, 7,317 species, 4,819 EC numbers; 13.27M of the rows
  (64.6%) are strain-resolved.
* `organism` is the name verbatim from the source; `genus` / `species` / `strain` are parsed out
  of it. Infraspecific ranks stay with the taxon (`Acinetobacter calcoaceticus subsp. anitratus`
  is the *species*, not a strain), and `str.` / `cf.` / quoting are handled — see
  `split_organism()`.
* **Nomenclature:** v7 follows current NCBI nomenclature, so reclassified genera
  (`Lactobacillus` → `Lacticaseibacillus`) reconcile directly against modern strain names, and
  taxa merged by NCBI are collapsed onto their surviving identifier. This retired the stale-
  taxonomy caveat the v6 export carried. Joining on `tax_id` is still the more robust practice.

## `food_nutrients.tsv`

`fdc_id`, `source_food_code`, `description`, `data_type`, `food_category`, `nutrient_id`,
`nutrient_name`, `unit_name`, `amount`, `source_db`

* `amount` is **per 100 g edible portion**, in `unit_name` units.
* `fdc_id` is an **accession** from a uniform 3,000,000-wide block per source, so
  `fdc_id // 3_000_000` gives the source's block index (blocks 0–2 are USDA FDC's own ids,
  never reassigned; the 14 other sources run consecutively from 9,000,000, and the resource
  ends below 51,000,000). Allocation is defined once in `food_DBs/fdc_blocks.py` and frozen in
  `fdc_id_map.tsv`; `export_resources.py` asserts containment before writing, because a shared
  block merges two unrelated foods silently instead of failing.
* `source_food_code` is the identifier the **source itself** publishes for that food — CIQUAL's
  `alim_code`, AFCD's `Public Food Key`, McCance's `Food Code`, and so on. It is a column rather
  than something encoded in `fdc_id`: ids used to be `block_base + native_code`, which made the
  block width depend on how sparsely a source numbered its foods instead of on how many it had
  (PhyFoodComp: 3,377 foods spanning 19,020,060 code values, 0.02% density).
* 120,310 foods × 1,779 nutrients. Sorted by `fdc_id`, then `nutrient_id`, so every food's
  nutrients are contiguous.
* **Branded label products are NOT here.** `branded_food` was 1,890,275 of the 2,010,585 foods
  and 25,937,648 of the 28,246,465 values, yet declared only 119 distinct components against
  1,779 for everything else — a nutrition-facts panel, not an analysis. Dropping them makes the
  export match the predictor, which has always run `drop_branded: true`. It also rebalances
  provenance: `fdc` falls from 96% of all values to 50%. Pass `--keep_branded` to get the old
  table back. Nine components exist only on branded labels and go with them (1068 beta-glucans,
  1072 other carbohydrate, 1086 sugar alcohols, 1124/2068 label vitamin E, 1181 inositol,
  1235/1236 added and intrinsic sugars, 1368 EGCG), as do 18 of the 187 `synthetic_bacterial`
  rows that had been injected onto branded foods.
  **After changing this, regenerate `/data/bac2food/live_nutrients.tsv` and rebuild
  `0_building/3_nutrient_to_ec.tsv` and `--only enzymes`** — `in_model` is computed against the
  nutrients that carry a measured value, so it goes stale otherwise (`verify_exports.py` catches
  it as orphaned nutrient_ids).
* `data_type` is the FDC provenance (`foundation_food`, `sr_legacy_food`, `survey_fndds_food`,
  `sub_sample_food`, …). All 14 non-USDA national tables enter as `foundation_food`.
* `food_category` is resolved from whichever of the three places holds it: the merged non-USDA
  foods carry it directly, FDC curated foods carry a numeric id looked up in `food_category.csv`,
  and branded foods carry their own free-text label (only reachable under `--keep_branded`).

### `source_db` — which reference database the value came from

Taken from the file the row sits in inside `food_nutrient_bucketed/`: `merge_phase8_v2.py` writes
each national database to its own `<src>_data.parquet`, and the original FoodData Central rows live
in the `part-*.parquet` files.

| label | database | rows |
|---|---|---|
| `fdc` | FoodData Central (USDA) | 1,156,398 |
| `mccance` | McCance & Widdowson's (UK) | 228,014 |
| `japan` | STFCJ (Japan) | 211,667 |
| `afcd` | AFCD / ASNUT (Australia) | 166,263 |
| `cnf` | Canadian Nutrient File | 130,091 |
| `biofoodcomp` | BioFoodComp (FAO) | 109,417 |
| `fineli` | Fineli (Finland) | 87,880 |
| `ciqual` | CIQUAL (France) | 52,683 |
| `nevo` | NEVO (Netherlands) | 51,522 |
| `swedish` | Livsmedelsdatabasen (Sweden) | 30,916 |
| `swiss` | Swiss Food Composition DB | 25,677 |
| `wafct` | WAFCT (West Africa) | 22,550 |
| `phyfoodcomp` | PhyFoodComp (FAO) | 15,382 |
| `frida` | Frida (Denmark) | 12,904 |
| `phenol_explorer` | Phenol-Explorer (France) | 7,036 |
| `synthetic_bacterial` | curated bacterial substrates injected by `0_building/inject_bacterial_substrates.py` (HMOs, xylan, …) | 169 |
| `mccance;swedish` | reported **identically** by both | 248 |

**When two databases report the same (food, nutrient):** if they give the *same* amount it is one
fact with two witnesses, so the rows collapse into one and `source_db` lists both, `;`-separated
(`mccance;swedish` — the only such pair, 248 rows). If they give *different* amounts, the rows are
**kept separate**, each with its own `source_db`. Merging those would invent a value neither
database reports; the disagreement is real and you should see it.

So a `(fdc_id, nutrient_id)` pair is **not** unique in this file. 42,362 pairs carry more than one
row because the amounts conflict:

* **17,194 conflict across databases** — two national DBs measured the same food differently.
* **25,168 conflict within a single database** — the source itself reports one nutrient twice with
  different values for one food (e.g. McCance gives Niacin 0.6 *and* 1.1 mg for canned ackee).
  Concentrated in `afcd` (37k rows), `japan`, `wafct`, `mccance`. This is a **pre-existing data
  issue in the ingested sources**, not something the export introduces — but it means anything that
  sums or averages this file must decide how to reconcile them first.

Tabs and newlines are stripped from every text field, so each record is exactly one line.

## Verifying the exports

`python verify_exports.py` health-checks the deposit: schema and per-row field counts,
counts against the figures the manuscript quotes, deposit hygiene (only deliverable files
present), and referential integrity across the three files — including that every
`in_model` row reaches a nutrient that actually exists in the composition table. Run it
after any export rebuild; it exits non-zero on failure.

The deposit is exactly three files — `food_nutrients.tsv`, `enzyme_substrate_chebi.tsv` and
`species_enzymes.tsv`. The legacy `species_enzymes.v6.tsv` was retired when the organism→EC
layer moved to eggNOG v7; regenerate it from `/data/bac2food/bact_ec.tsv` (retained) if a v6
comparison is ever needed. Derived caches belong in `/data/bac2food/cache/`, never here.
