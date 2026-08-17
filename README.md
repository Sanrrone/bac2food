# bac2food

Predict which **foods** feed the bacteria in a metagenome, from their enzymes.

The chain is: a bacterium carries **enzymes** (EC numbers) → each enzyme acts on **substrates**
(resolved to ChEBI) → substrates correspond to **nutrients** (FoodData Central) → nutrients are
measured in **foods**. Composing those links gives a bacterium → food score.

```
   eggNOG v7            BRENDA              ChEBI      FoodData Central + 14 national DBs
       │                   │                   │                        │
       ▼                   ▼                   ▼                        ▼
  bacterium → EC     EC → substrate    substrate → ChEBI         food → nutrient
  (eggnog/6.0_+6.1_,  (eggnog/)            (chebi/)               (food.parquet,
   via KEGG KO)             │                   │              food_nutrient_bucketed/)
       │                    └─────────┬─────────┘                        │
       │                              ▼                                  │
       │                0_building: EC ←→ nutrient                       │
       │                   3_nutrient_to_ec.tsv                          │
       │                              └───────────────┬──────────────────┘
       │                                              ▼
       └──────────────────────────────► scoring: bacterium → food
                                          4_predict
```

Scores are predictions of enzymatic **capability**, not of measured growth. A bacterium that
carries the enzyme to degrade a substrate is not thereby shown to do so in vivo.

**`4_predict/bac2food_predict.py` is the entry point for real data.** `0_building` is the
one-time build of the reference maps; `5_export` packages the release.

## Install

```bash
pip install -r requirements.txt        # Python 3.10; pandas + pyarrow and three file readers
```

Or take the environment as a container, which is the reproducible route. The container ships as
a **Singularity/Apptainer image** — a single `.sif` file that runs as you, with no daemon and no
root, so it works on the HPC systems where most of this analysis actually happens:

```bash
apptainer build bac2food.sif bac2food.def          # or: singularity build …

apptainer run --cleanenv \
    --bind /path/to/deposit:/data/bac2food/exports:ro \
    --bind "$PWD/work":/data/bac2food/cache \
    bac2food.sif --mag_tsv sample_ec.tsv --out_prefix run1 --jobs 4
```

`bac2food.def` builds by cloning this repository rather than by copying your working directory.
Apptainer has no ignore-file mechanism, so a wholesale copy would sweep in whatever happens
to be sitting in the tree — the raw source tables, the 248 MB ontology, the cohort annotations.
Cloning makes the image hold exactly what was published, and nothing else.

**Two things are worth knowing before you run it.**

*Use `--cleanenv`.* Apptainer inherits the host environment and auto-binds `$HOME`, `/tmp` and
the current directory. Without it, a host `PYTHONPATH` or a stray package in `~/.local` reaches
into the run, and what you measured is no longer the container.

*Memory is not capped by the container.* Apptainer does not limit memory itself, so the ceiling is
whatever the host or the scheduler gives you — under SLURM, ask for `--mem=8G`. The predictor
peaks near 5.7 GB, driven by loading the reference layer rather than by the community size, and
below roughly 8 GB the Linux OOM killer takes it. That surfaces as an unexplained exit 137, not
as a Python error, so it is worth recognising.

The image holds the code only. Reference tables are bound at run time, never baked in: they are
2.1 GB, they are versioned separately on Zenodo, and not every source's licence permits
redistribution inside an image.

## What this repository does and does not contain

Code and small derived tables are tracked here. Three classes of file are deliberately absent
(see `.gitignore`, which explains each):

| absent | why |
|---|---|
| raw source databases (`food_DBs/**/*.xlsx`, …) | provider-licensed, and the terms differ by source — some carry NonCommercial or ShareAlike conditions, and **at least one does not permit redistribution at all**. Each folder's readme says where to obtain its source. |
| cohort annotations (`gene_annot/`) | third-party metagenome annotations, and their filenames carry sample identifiers |
| `chebi/chebi.obo` (248 MB), `eggnog/*.tar.xz` (117 MB) | third-party downloads that exceed GitHub's 100 MB file limit |
| `brenda/brenda_2026_1.*.tar.gz` (154 MB) | the official BRENDA release. Its own README states the full contents are copyright-protected. The **derived** digest is what the pipeline reads, and that is tracked (`eggnog/2_digest_*.tsv`) — download the release from [brenda-enzymes.org](https://www.brenda-enzymes.org/download.php) only if you need to rebuild it. |
| the exports themselves (`5_export/exports.tar.xz`) | that archive *is* the Zenodo deposit, and `food_nutrients.tsv` inside it carries rows derived from a source whose terms do not permit redistribution. This repository ships the code that builds the deposit, never the deposit. |

What *is* tracked for the restricted sources is `0_building/novel_nutrients/novel_*.csv`: the
component names and units each source declares, with no values. Where a source cannot be
redistributed, that mapping is the whole bring-your-own-source route — hold your own licensed
copy and the pipeline reconstructs those rows locally.

## Layout

| folder | role |
|---|---|
| `eggnog/` | Two layers, two sources. **EC → substrate** digest (`2_digest_norm.tsv`) comes from **BRENDA**, kept, not dropped. **Bacterium → EC** comes from **eggNOG v7** (`6.0_`+`6.1_` → `exports/species_enzymes.tsv`). What was dropped was BRENDA's *species→enzyme* scrape (`3_normalized_species.tsv`), which under-reported every organism's enzymes; eggNOG replaced that role only. The folder name is a leftover of that swap. |
| `chebi/` | substrate name → **ChEBI id** (`digest_to_chebi.tsv`), via `chebi.obo`. |
| `0_building/` | builds the core map **`3_nutrient_to_ec.tsv`** (nutrient ↔ EC, walked over the ChEBI ontology). |
| `food_DBs/` | ingests 14 non-USDA national and regional food databases onto the FDC food tables. |
| `4_predict/` | **the predictor** — a whole metagenome in, three ranked food tables out. |
| `5_export/` | flatten the reference layers into **shareable TSVs** (`export_resources.py`). |
| `analysis/` | **the paper** — every figure and every number reported in the Data Descriptor. Not part of the pipeline; reads the deposit and a directory of predictor outputs. |

### File index — the scripts that matter

| file | what it produces |
|---|---|
| `4_predict/bac2food_predict.py` | community, differential and per-food tables for one metagenome |
| `4_predict/parameters.yaml` | every scoring constant (override with `--config`) |
| `4_predict/chain_coverage.py` | where the enzyme→food chain terminates, per EC |
| `4_predict/recompute_match_rate.py` | the feature and species reference match rates |
| `rescue_bifunctional_ec.py` | restores the second EC of multi-activity CAZymes in an EC panel |
| `bifunctional_ec.tsv` | the curated product-name → EC table it applies (`.readme.txt` explains it) |
| `5_export/export_resources.py` | the three deposited TSVs |
| `5_export/verify_exports.py` | 24 structure, count and join checks over the deposit |
| `0_building/2_nutri2chebi_from_obo.py` | `2_nutrient_to_chebi.tsv` (nutrient → ChEBI) |
| `0_building/3_nutrient_to_ec.py` | `3_nutrient_to_ec.tsv` — **the core map every scorer reads** |
| `food_DBs/bucket_food_nutrient.py` | the Hive-partitioned `food_nutrient_bucketed/` store |
| `food_DBs/merge_phase8_v2.py` | merges every national source into `food.parquet` + the store |
| `food_DBs/audit_parsers.py` | per-source ingestion audit |
| `chebi/dict_to_chebi.py` | `digest_to_chebi.tsv` (substrate → ChEBI) |
| `eggnog/6.0_`, `6.1_` | `species_enzymes.tsv` (organism → EC, via KEGG orthology) |
| `analysis/prep_figure_data.py` | the five figure input CSVs, and `anonymize()` — **run it first** |
| `analysis/make_figures.R` | Figures 1 and 2 |
| `analysis/linkage_walk.py` | the linkage walk over the deposited TSVs alone, importing no part of bac2food |

`analysis/readme.txt` lists the rest and gives the two environment variables that point them
at the deposit and at a predictor run.

Large inputs live in `/data/bac2food`: `food.parquet`, `food_nutrient_bucketed/`, `nutrient.csv`,
`food_category.csv`, `index_modeled/`, and `bact_ec.tsv` (the legacy v6 organism→EC layer, kept
only so a v6 build stays reproducible).

## Run the predictor

```bash
cd 4_predict
python bac2food_predict.py \
    --mag_tsv /path/to/sample_ec.tsv \
    --out_prefix sample \
    --jobs 8
```

`--mag_tsv` auto-detects its input: either a header TSV (`species`, `ec_number`, optional
`strain`) or a headerless annotation TSV whose species and EC columns are found **by content**,
with multi-EC cells (`2.7.2.3,5.3.1.1`) split automatically. One file is one community. Every
constant lives in `4_predict/parameters.yaml`; see `4_predict/README.md` for the four output
tables and the scoring kernel.

Four tables, one per question:

| question | table |
|---|---|
| which foods best feed this whole microbiota | `community.tsv` |
| which foods best suit a given organism | `perBacterium.tsv` |
| which foods a given organism uses better than its peers | `differential.tsv` |
| which organisms a given food favours | `perFood.tsv` |

## Rebuild the reference maps

Only needed if a source database changes. Each folder's `readme.txt` holds the exact commands.
The order matters, and steps 4–6 in particular are not optional — `in_model` is computed against
the nutrients that carry a measured value, so changing the food set without rebuilding the map
leaves the digest citing nutrients no food reports:

```
1. eggnog/     tar -xJf ec_species_substrate.tar.xz          # BRENDA (EC → substrate)
               1.0_eggnog_ec_substrates_parser.py  →  2_digest_dict.tsv
               1.5_reactions_to_digest.py          →  2_digest_norm.tsv

2. chebi/      dict_to_chebi.py --digest ../eggnog/2_digest_norm.tsv
                                                   →  digest_to_chebi.tsv
               (2_digest_NORM, not 2_digest_dict — the norm step is what removes
                ubiquitous cofactors; the unfiltered dict resolves protons and GTP
                as dietary substrates)

3. food_DBs/   merge_phase8_v2.py                  →  food.parquet + food_nutrient_bucketed/

4. 5_export/   export_resources.py --only foods    →  food_nutrients.tsv
               then regenerate /data/bac2food/live_nutrients.tsv from it

5. 0_building/ 0_name_normalizer.py                →  0_nutrient.normalized.tsv
               1_nutrient_expansion.py             →  1_expanded_nutrients.tsv
               2_nutri2chebi_from_obo.py           →  2_nutrient_to_chebi.tsv
               3_nutrient_to_ec.py                 →  3_nutrient_to_ec.tsv     ★ the core map

6. 5_export/   export_resources.py --only enzymes  →  enzyme_substrate_chebi.tsv
               verify_exports.py                   →  expect 24/24
```

After any of this, wipe the predictor's derived caches or it will serve a stale food universe:

```bash
rm /data/bac2food/index_modeled/{static_food_meta.pkl,*.parquet}
```

## Reproducibility check

The derived index is a pure function of the reference tables, so wiping it must change
nothing. That is the cheapest end-to-end check that an environment is sound, and it is the
one to run first inside a fresh container:

```bash
cd 4_predict
md5sum /data/bac2food/index_modeled/* > /tmp/index_before.md5
rm /data/bac2food/index_modeled/{static_food_meta.pkl,*.parquet}
python bac2food_predict.py --mag_tsv <sample>_ec.tsv --out_prefix /tmp/check --jobs 6
md5sum -c /tmp/index_before.md5      # all five files must report OK
diff /tmp/check.community.tsv <previous run>.community.tsv
```

Verified 2026-08-05 on Python 3.10.8: all five index artifacts rebuild with identical
checksums, and all four predictor output tables reproduce byte for byte.

There is a second exact check on the map side. `0_building/3_nutrient_to_ec.py
--no_structural` reproduces the pre-2026-08 map byte for byte (md5
`0148588f0a29e458ffcb3866973737a5`), so any future change to the ontology walk can be
rolled back and confirmed rather than argued about.

## Reference exports

`5_export/export_resources.py` writes the three reference layers as flat TSVs to
`/data/bac2food/exports`: `enzyme_substrate_chebi.tsv`, `species_enzymes.tsv` and
`food_nutrients.tsv`. These are the tables to deposit alongside a paper; see `5_export/README.md`.

**`food_nutrients.tsv` holds analytically characterized foods only.** FDC's branded label
products are excluded by default: they were 1,890,275 of 2,010,585 foods and 92% of all values,
yet declared only 119 distinct components against 1,779 for everything else. The predictor had
always ignored them (`drop_branded: true`), so the deposit now matches the food set the software
scores against, and the USDA share of all values falls from 96% to 50%. Pass `--keep_branded` to
restore them.

## Known caveats

* **Taxonomy.** The shipped organism→EC layer is eggNOG v7: current NCBI nomenclature, 10,751
  organisms, 4,819 EC numbers. Joining on `tax_id` rather than on the name is still the more
  robust practice. The legacy v6 layer (`bact_ec.tsv`, older names, 3,176 organisms) is retained
  only for reproducing v6 builds.
* **Portions.** `/data/bac2food/food_portion.csv` is optional and currently absent; without it
  the predictor assumes 50 g per food and says so at startup.
* **Match rate.** About 45.5% of annotated features in a metagenome reach a food nutrient
  (per-sample range 43.4–46.7% across the 55-sample cohort). The cause of the remainder is the
  nutrient vocabulary, not enzyme–substrate coverage: most bacterial enzymes act on replication,
  signalling and envelope substrates that no composition table reports.
* **Reference augmentation.** `--augment_with_reference` fills a species' EC set from a matched
  reference genome when the metagenome annotated fewer than `--augment_threshold` (default 200)
  ECs for it, on the premise that a shallow annotation understates the organism rather than
  describing it. It is off unless asked for; with it off, the enzyme pool is exactly what the
  supplied annotation contains, and completing it is the caller's responsibility. Measured on
  four cohort samples spanning 14–63 taxa, it fires on 2–5 species per sample. Across the full
  community ranking (77–80 shared foods) the two arms agree at Spearman ρ 0.800–0.992 / Kendall
  τ 0.740–0.947; the top ten is identical in three of four, and the fourth reorders its head
  (Jaccard 0.54) while its body holds (ρ 0.800). It does not homogenize: mean similarity
  *between* different bacteria's top-ten foods moves by −0.030 to +0.003. The threshold is what
  buys that — set it to 0 and every matched species collapses onto the same reference profile.
* **Multi-activity enzymes.** An annotator gives a locus one EC number, so the second activity of
  a bifunctional CAZyme is lost. `rescue_bifunctional_ec.py` restores it from a curated
  product-name table; see `bifunctional_ec.readme.txt` for what it does and does not map.
* **Memory.** The predictor peaks near 5.7 GB, most of it loading the reference layer. Run one
  at a time on a 16 GB machine.
