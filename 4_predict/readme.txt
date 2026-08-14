4_predict/ — Bacteria <-> Food predictor for metagenomes
=========================================================

bac2food_predict.py emits FOUR result files from a single scoring pass, one
per question a user actually asks:

  <prefix>.community.tsv     — which foods best feed this whole microbiota.
                               Each nutrient is weighted by the SHARE of the
                               community able to act on it (--community_weight,
                               default 'membership'). It used to score the plain
                               union of every organism's ECs; that saturates —
                               a gut community reaches ~79% of mappable
                               nutrients regardless of richness — so the union
                               could not tell two communities apart. Pass
                               --community_weight none to reproduce it.
  <prefix>.perBacterium.tsv  — which foods best suit a GIVEN organism, whether
                               or not they also suit its neighbours. Absolute
                               score, not relative. NEW 2026-08.
  <prefix>.differential.tsv  — which foods a given organism exploits BETTER
                               than its peers (comp_score = score -
                               peer_median across bacteria)
  <prefix>.perFood.tsv       — which bacteria a given food favours (top-K)

perBacterium and differential are not interchangeable. A food can be an
organism's single best substrate AND carry its most negative comp_score,
because every co-occurring organism uses it better still. Nor can you recover
perBacterium by inverting perFood: that table keeps only the top
--top_bacteria_per_food organisms per food, so an organism appears for ~8
foods rather than its whole shortlist.

Optional 4th file: <prefix>.complement_ec.tsv (only with --complement_ec)
diagnoses annotator gaps vs the bact_ec.tsv reference.


Input
-----
TSV with columns:
    species       e.g. "1134687_Klebsiella_michiganensis" or "Bacteroides thetaiotaomicron"
    ec_number     e.g. "3.2.1.1"   (or column may be named "ec")
    strain        optional; appended to bacterium label

Extra columns (e.g. a sample id) are ignored.


Canonical invocation
--------------------
python bac2food_predict.py \
    --mag_tsv ../my_cohort_ec.tsv \
    --nutrient_to_ec ../0_building/3_nutrient_to_ec.tsv \
    --food_nutrient /data/bac2food/food_nutrient_bucketed \
    --food          /data/bac2food/food.parquet \
    --nutrient      /data/bac2food/nutrient.csv \
    --food_category /data/bac2food/food_category.csv \
    --food_portion  /data/bac2food/food_portion.csv \
    --index_dir     /data/bac2food/index_modeled \
    --out_prefix    hep \
    --max_foods     10 \
    --spec_alpha    2.0 \
    --augment_with_reference --augment_threshold 200 \
    --complement_ec \
    --jobs          6

Outputs:
    hep.community.tsv      hep.differential.tsv      hep.perFood.tsv
    hep.complement_ec.tsv  (only with --complement_ec)


Just the broader bacterial universe (no user MAGs)
---------------------------------------------------
python bac2food_predict.py \
    --use_reference \
    --out_prefix ref_pool \
    --ref_min_ec 100 \
    --ref_max_species 500 \
    --top_bacteria_per_food 20 \
    --jobs 8


Notes
-----
* The static_food_meta and bucketed parquet at /data/bac2food/index_modeled/
  and /data/bac2food/food_nutrient_bucketed/ are reused as-is. They rebuild
  automatically if missing or if --rebuild-static-meta is set.
* After pulling a new version of bac2food_predict.py, wipe the cache so
  the new canonicalization / category-fallback rules take effect:
      rm /data/bac2food/index_modeled/{static_food_meta.pkl,*.parquet}
* --max_foods is the per-bacterium scan budget AND the output row cap
  for community + differential.
* Default paths point to /data/bac2food/* — override any of them with the
  matching CLI flag.
* See README.md for the full reference (CLI, canonicalization, three-file
  output schemas, validation panel, known limitations).
