# 4_predict — Bacteria ⇄ Food predictor for metagenomes

`bac2food_predict.py` takes a metagenome's enzyme annotation and answers
three distinct questions in three output files:

| file | answers | when to read it |
|---|---|---|
| `<prefix>.community.tsv` | **What foods support the *whole* microbiome?** Treats the union of every input bacterium's ECs as a single pseudo-bacterium; foods are ranked by how broadly they feed the community. | Diet design for a whole microbiome. |
| `<prefix>.differential.tsv` | **What foods does each individual bacterium do *better* at than its peers?** Per-bacterium rows; a food is admitted where `score > peer_median` and the survivors are ranked by absolute `score`. Surfaces specialty foods that absolute scoring buries. | Personalization, finding specialty substrates. |
| `<prefix>.perFood.tsv` | **For each food, which bacteria benefit most?** Inverse of the differential — every scored food paired with its top-K bacteria. | Predicting which species each food upregulates. |

These are the *only* result files. An optional fourth file
`<prefix>.complement_ec.tsv` is emitted when `--complement_ec` is passed —
it's a diagnostic, not a result.


## The two formulas

Both modes share one scoring kernel and then diverge in **three** places:
how nutrients are weighted, whether foods are scored against each other,
and how the final list is chosen. Everything below is the whole model —
there is no fourth hidden term.

### Shared kernel — scoring one food `F` for one organism `B`

**1. Weight each nutrient.** All three factors are properties of the
*nutrient*, never of the organism:

```
w[n] = fibre_boost(n) × nmult[n] × food_IDF(n)

  fibre_boost  = fiber_weight (2.0) if the nutrient is a fibre/starch/
                 oligosaccharide/inulin/pectin/glucan/cellulose, else 1.0
  nmult        = manual override, 1.0 unless set
  food_IDF     = 1 + ln(total_foods / foods_carrying_n)
```

**2. Gain — how much of `B`'s substrate demand `F` actually meets:**

```
gain = Σ_n  w[n] · [ ln(1 + (cur[n] + a_n)/τ_n) − ln(1 + cur[n]/τ_n) ]

  a_n  = amount of nutrient n per 100 g of F
  τ_n  = per-nutrient scale (1 % of the corpus reference level)
  cur  = what the foods already chosen supply — ZERO in differential mode
```

The `ln(1+x)` shape is what makes a second helping of the same substrate
worth less than the first.

**3. Effective gain** — diminishing returns, plus a bonus for meeting
*many* distinct needs rather than one loudly (`ni` = nutrients improved):

```
geff = ln(1 + gain) · √(ni / (ni + 2))
```

**4. Score** — effective gain minus what the food costs you:

```
score = geff · (1 + type_bonus) − food_baseline − amount_cost

  type_bonus    = type_w · (data-quality tier / 600)
  food_baseline = base_cost(F) + art_w · artifact_flag(F)
  amount_cost   = (proc_w · processing(F) + broad_w · (1 − purity)^broad_q)
                  · ln(1 + gain)
```

`food_baseline` is a flat toll per food; `amount_cost` grows with the gain
claimed, so a food cannot buy a high score purely by being nutrient-dense
or by being a broad-spectrum powder.

### Community mode — *"what feeds this whole microbiome?"*

1. Treat the **union** of every organism's targets as one pseudo-organism.
2. Weight each nutrient additionally by **carrier share** — the fraction of
   the community that can actually use it, `× share^coverage_alpha`. This
   is the one term community mode has and differential does not, and it is
   why the table moves with composition instead of mere presence.
3. Choose foods **greedily as a set**: pick the best food, add its nutrient
   amounts into `cur`, repeat. Two extra costs apply while choosing, so the
   result is a *complementary basket* rather than N variations on one food:

```
score_greedy = score − overlap_w · redundancy − 0.20 · category_repeats
```

### Differential mode — *"what does this organism do better than its neighbours?"*

1. Score each organism separately, with `cur = 0` — every food scored from a
   clean slate, independently of every other food.
2. No carrier-share weight, no greedy costs.
3. **Admit** a food on two tests, then **rank** on the score itself:

```
admit  F  if  score > 0  and  score > peer_median[F]
rank   admitted foods by score, descending
cap    at most max(2, 0.4 × max_foods) foods per category or lead substrate
```

The peer median is a **gate, not the ranking key**. Subtracting it would
remove any substrate the peers share — which is exactly what defines a
guild — so four starch degraders would cancel on starch and be ranked on
whatever trace compound happened to differ. See *Ranking rule*.


## Inputs

`--mag_tsv` accepts a bacterium→EC TSV in **either** of two shapes — the format
is auto-detected, so no flag is needed to switch between them.

**1. Header TSV** with these columns (case-insensitive):

| column | required | example |
|---|---|---|
| `species` (or `taxon` / `organism`) | yes | `1134687_Klebsiella_michiganensis` or `Bacteroides thetaiotaomicron` |
| `ec_number` (or `ec`) | yes | `3.2.1.1` |
| `strain` | optional | appended to the displayed bacterium label |

Extra columns (e.g. a `sampleID` column) are ignored.

**2. Headerless annotation TSV** (e.g. `gene_annot/*_ec.tsv`) — no header row;
the **species column and the EC column are recognized by content**, so the
column order is not fixed:

```
sample_A   47715_Lacticaseibacillus_rhamnosus   l_46433   3.4.11.18
sample_A   1522_Clostridium_innocuum            l_07349   2.7.2.3,5.3.1.1
```

- the **species** column is detected as the one holding organism names (a
  `<taxid>_Genus_species` form scores highest); sample-id and gene-locus columns
  are ignored;
- the **EC** column is detected as the one holding EC numbers, and **a single
  cell may list several comma- or semicolon-separated EC numbers** (e.g.
  `2.7.2.3,5.3.1.1`) — each is expanded to its own bacterium→EC edge. An optional
  `EC:` prefix is tolerated.

On load, the script prints which columns it detected, e.g.
`[*] Headerless input (4 cols); auto-detected species=col1, ec=col3.`

> Note: each input file is treated as **one community**. To analyze several
> samples, run the script once per file (it does not loop over a folder).


## Quick start

The two views answer different questions and are produced by different scoring
passes, so each is its own subcommand and exposes only the options that apply
to it.

**Which foods feed the whole microbiome?**

```bash
python bac2food_predict.py community \
    --mag my_metagenome.tsv --out my_run \
    --max_foods 10 --weight membership --augment_with_reference --jobs 6
```
→ `my_run.community.tsv`

**What does each organism use better than its neighbours?**

```bash
python bac2food_predict.py differential \
    --mag my_metagenome.tsv --out my_run \
    --max_foods 10 --rank score --diversity -1 --jobs 6
```
→ `my_run.differential.tsv`, plus `my_run.perFood.tsv` and
`my_run.perBacterium.tsv`, which fall out of the same pass at no extra cost.

Each subcommand runs **only** its own pass. On a 63-organism sample the
community view takes 24 s against 60 s for both, so asking for one view no
longer pays for the other.

Run `bac2food_predict.py community --help` or `... differential --help` to see
just that view's options: `--weight`, `--coverage_alpha` and `--abundance_tsv`
appear only under `community`; `--formula`, `--rank`, `--diversity` and
`--min_peers` only under `differential`.

**Legacy (no subcommand).** Invoked with no mode the script runs both passes
and writes every view, exactly as before:

```bash
python bac2food_predict.py \
    --mag_tsv my_metagenome.tsv --out_prefix my_run \
    --max_foods 10 --augment_with_reference --complement_ec --jobs 6
```

This path is deprecated but supported, and is the faster choice when you want
both views. The deposited run scripts use it, and every long option keeps its
original spelling (`--mag_tsv`, `--out_prefix`, `--community_weight`,
`--differential_rank`, …), so nothing published against the old interface
changes behaviour. Verified byte-identical across all four output files.

Expected runtime: ~1.5 min on 600 species with 6 cores (after one-time
~2 min index build).


## CLI reference

### Required

- `--out_prefix PATH` — output files written as `<PATH>.<suffix>.tsv`
- `--mag_tsv PATH` — input EC annotation. Required unless `--use_reference`
  is set (in which case the bacterial pool comes from `bact_ec.tsv`).

### Mode-specific

Registered only under the subcommand they belong to (and all still present in
the legacy no-subcommand path):

**`community`**
- `--weight` / `--community_weight {membership,abundance,none}` *(default `membership`)*
- `--coverage_alpha FLOAT` *(default 1.0)*
- `--abundance_tsv PATH` — required by `--weight abundance`

**`differential`**
- `--formula` / `--differential_formula {full,explicit_admission,gain_only}` *(default `full`)*
- `--rank` / `--differential_rank {score,comp_score}` *(default `score`)*
- `--diversity` / `--differential_diversity K` *(default `-1` = `max(2, 0.4×max_foods)`)*
- `--min_peers` / `--diff_min_peers N` *(default 20)*
- `--top_bacteria_per_food N` *(default 20)* — perFood cap

### Shared

- `--mag` / `--mag_tsv PATH`, `--out` / `--out_prefix PATH`
- `--max_foods N` *(default 10)* — both the per-bacterium scan budget AND
  the per-bacterium / community row cap in the output. It is **not** a
  quality knob for the differential view, which is fed the clean-slate
  frame where every candidate is scored independently: grading the top 5
  at fixed depth gives an identical MRR 0.296 at 5, 10, 20, 50 and 100.
  Its one biological channel is that it sets the diversity cap, and 10
  (cap 4) is the best point on that curve — MRR 0.308 vs 0.296 uncapped.
- `--allow-spices` / `--no-allow-spices` *(default off)* — include the
  `Spices and Herbs` category.
- `--drop_category CAT` — drop a food category entirely; repeatable.

### Reference (`bact_ec.tsv`)

- `--complement_ec` — write `<prefix>.complement_ec.tsv` diffing the user's
  EC annotation against the reference. Report-only.
- `--augment_with_reference` — for matched species, merge missing
  reference ECs into the user's effective EC set before scoring.
- `--augment_threshold N` *(default 200)* — only augment species with
  fewer than N user ECs (avoids homogenizing well-annotated MAGs).
- `--use_reference` — also score the reference pool from `bact_ec.tsv`
  (in addition to / or instead of `--mag_tsv`).
- `--ref_min_ec N` *(default 20)* / `--ref_max_species N` *(default 0)*
  — filters when `--use_reference` is on.
- `--bact_ec_tsv PATH` *(default `/data/bac2food/bact_ec.tsv`)*

### Data paths

The defaults come from `parameters.yaml` (see below); each can still be
overridden per-run on the CLI:
- `--nutrient_to_ec` *(default `../0_building/3_nutrient_to_ec.tsv`)*
- `--food_nutrient` — bucketed parquet directory
- `--food` / `--nutrient` / `--food_category` / `--food_portion`
- `--nutrient_alias` *(default `../0_building/1_expanded_nutrients.tsv`)*
- `--index_dir` *(default `/data/bac2food/index_modeled`)*

  The derived index is keyed on the **identity** of the store it was built from, recorded in
  `<index_dir>/index_inputs.json`, not only on mtimes. Pointing `--food_nutrient` at a
  different store therefore rebuilds rather than silently reusing the previous store's index.
  It used to reuse it: an mtime check passes whenever the other store's files happen to be
  older than the cached index, which is the normal case for an unmodified canonical store, and
  the run then reports plausible numbers computed against the wrong food universe. If you keep
  several stores, give each its own `--index_dir` anyway — the guard makes the failure loud,
  but a separate directory avoids re-indexing on every switch.

### Runtime

- `--jobs N` *(default = nproc)*
- `--rebuild-static-meta` — force rebuild of the cached static food meta.
- `--config PATH` — use an alternate `parameters.yaml` (see below).

## Configuration — `parameters.yaml`

The constant, non-user-selectable settings live in **`parameters.yaml`** next to
the script, so the model can be tuned without editing the Python source. It holds:

- **`paths`** — default data-file locations (each also overridable by its CLI flag);
- **`flags`** — `drop_branded`, `allow_macro_proxy`, `allow_macro_scan`;
- **`scoring`** — the 13 kernel scalars (`proc_w`, `broad_w`, `broad_q`,
  `type_w`, `art_w`, `overlap_w`, `fiber_weight`, `outlier_ratio`, the
  `k_*` caps, the base weights). Nine former entries (`tau_q`, `min_contrib_frac`, `k_min`,
  `gain_min`, `score_min`, `max_canons`) were removed in Phase 11: the
  first three were dead code — `tau_q` and `min_contrib_frac` sat behind an
  `ndf >= 50000` branch that could never run (the largest EC-linked
  nutrient reaches 15,258 foods) and `k_min` was never referenced. The
  other three are still enforced, as internal constants rather than knobs:
  each was perturbed by a large factor against both a 22-species panel and
  a 63-organism cohort sample without changing one output byte, so they
  guard against pathological input without shaping any result. Three more
  (`lambda_food`, `oligo_floor_mg`, `isoflavone_floor_mg`) were removed
  after a term-removal sweep: both floors sit below every amount observed,
  so the `max()` guarding them never selected them, and `lambda_food` is a
  constant subtracted from every food, which cannot reorder anything.
  Deleting all three left every ranking identical on the panel and on two
  cohort samples, with every score shifted by exactly +0.0800 — the value
  of `lambda_food`. That offset cancels in `comp_score`;
- **`nutrient_ids`** — `oligo_ids`, `isoflavone_ids`, `key_nutrient_ids` (the
  cultivar nutritional-similarity gate);
- **`tables`** — `cat_penalty`, `plant_cats`, `always_drop_cats`, `tp_map`,
  `dangerous_rules`, `safe_rule_pfxs`.

Relative paths in the file resolve against the **`parameters.yaml` directory**
(i.e. `4_predict/`), independent of the working directory. The file is loaded at
import time so multiprocessing workers bind the same values.

Resolution order for the config file:
1. `--config PATH`, else
2. `$BAC2FOOD_PARAMS`, else
3. `parameters.yaml` next to `bac2food_predict.py`.

Per-run knobs that *do* depend on user choice stay as CLI flags and are **not**
in the YAML (`--max_foods`, `--differential_diversity`, `--jobs`,
`--out_prefix`, `--mag_tsv`, …). Logic/structural constants (the EC regex, the
parquet `BUCKETS` count, the text regexes and cultivar/category word-lists) also
stay in the source on purpose.


## How each file is computed

### `<prefix>.community.tsv`

1. Build `community_targs` = union of every input bacterium's effective
   target nutrient set (post-augmentation).
2. Compute a COVERAGE WEIGHT per nutrient: the share of the community able to
   act on it (`--community_weight membership`, the default; `abundance`
   weights organisms by a user-supplied `--abundance_tsv`; `none` restores the
   pre-2026-08 behaviour). This weight multiplies the per-nutrient weight in
   the kernel. Without it the table scores the plain union, which saturates —
   a single infant gut sample reaches 474 of 598 mappable nutrients (median,
   range 438-497) across a 4.5-fold range in taxa — so two communities of very
   different composition scored almost identically.
3. Score a singleton pool `{"[community]": community_targs}` with
   `specificity=False`. Returns the greedy top-`--max_foods` foods.
4. Compute `n_contributing_bacteria` per food: of the food's top
   nutrients, how many user bacteria target them.

Columns: `rank, food_name, representative_fdc_id, description,
food_category, n_variants, score, gain, n_nutrients_improved,
n_targets_total, covered_targets_total, coverage_total_frac,
n_contributing_bacteria, top_nutrient_ids, top_nutrient_names, data_type`.

### `<prefix>.perBacterium.tsv`

Which foods best suit a GIVEN organism, in absolute terms — the question
`differential.tsv` does NOT answer, because that table ranks by advantage over
peers rather than by suitability. A food every organism uses well is ranked
last there while being an excellent food for this one.

It is the same greedy per-bacterium frame `perFood.tsv` inverts, written
un-inverted, so it costs no extra scoring pass. It cannot be recovered from
`perFood.tsv`: that table keeps only the top `--top_bacteria_per_food`
organisms per food, so an organism surfaces for ~8 foods rather than its full
`--max_foods` shortlist, and only where it beat 19 competitors.

Columns: `bacterium, rank, food_name, representative_fdc_id, description,
food_category, n_variants, score, gain, n_nutrients_improved,
top_nutrient_names, data_type`.

### `<prefix>.differential.tsv`

1. Each per-bacterium worker returns two frames: the greedy
   `--max_foods` shortlist AND every candidate scored independently
   (clean slate, no greedy `cur` accumulation).
2. After all workers finish, peer median per food is computed across
   bacteria's clean-slate scores: `comp_score = score − peer_median`.
3. The peer comparison is applied as an **admission test**, not as the
   ranking key: a food is kept for a bacterium only where
   `score > peer_median`, and the survivors are ranked by absolute
   `score` descending, truncated to `--max_foods`.
   Pass `--differential_rank comp_score` for the historical rule (rank on
   the margin itself). See *Ranking rule* below for why the default
   changed.
4. The shortlist is then **diversity-capped**: at most
   `max(2, round(0.4 × --max_foods))` of the rows may share one food
   category or one lead substrate. See *Shortlist diversity*.
4. Per-bacterium diagnostics from the input EC annotation are merged in:
   `n_user_ec`, `n_ec_in_db`, `n_nutrients_targeted`.

Columns: `bacterium, rank, food_name, representative_fdc_id, description,
food_category, n_variants, comp_score, score, peer_median, peer_mean,
peer_n, gain, n_nutrients_improved, top_nutrient_ids, top_nutrient_names,
n_user_ec, n_ec_in_db, n_nutrients_targeted, data_type`.

### `<prefix>.perFood.tsv`

The greedy per-bacterium scores are inverted: group by food, sort
bacteria by score descending, keep top `--top_bacteria_per_food` per
food.

Columns: `representative_fdc_id, food_name, description, food_category,
n_variants, rank, bacterium, score, gain, n_nutrients_improved,
top_nutrient_names, data_type`.


## Food canonicalization

FDC fragments the same food across many rows: (a) per nutrient panel via
`- <panel> - NF<hex>` suffixes, (b) per preparation/state, (c) per
cultivar, (d) per brand, (e) per source vendor (BioFoodComp /
PhyFoodComp / AFCD / Frida / Phenol-Explorer / etc.). For metabolic
modelling almost all of this is noise.

The predictor groups every `fdc_id` sharing a canonical short name and
elects a single representative per group. The remaining variants' nutrient
data folds into the rep via **MEAN aggregation** (mean across actually-
measured variants — the parquet is sparse-by-row so unmeasured nutrients
have no row to drag the average down).

### Three-stage canonicalization

**Stage 1 — strip the noise** (`canonicalize_food_name`):
- `- <panel> - NF<hex>` suffixes
- preparation / state / packaging tokens (`raw`, `cooked`, `frozen`,
  `sliced`, `with/without salt`, `unsweetened`, `bottled`, `whole`,
  `mature`, `ripe`, `overripe`, …)
- trailing source brackets (`[BioFoodComp]`, `[AFCD]`, `[Frida]`,
  `[Phenol-Explorer]`)

**Stage 2 — head-aware chunk trimming** — split on commas:
- If the first chunk is in `_CULTIVAR_STRIP_HEADS` (apples, pears,
  peaches, plums, cherries, apricots, citrus, tropical fruits, grapes,
  melons, all berries, juices, plus BioFoodComp species-level entries
  like `quinoa` / `common bean` / `cassava` / `lettuce` / `tomato`, plus
  Phase 5 additions `nectarine`, `winged bean`, `wild mango`, `baobab
  fruit`, `anchote`, `locust bean`, `mushroom` singular, `potato tuber`,
  etc.), drop everything after the first comma.
- Otherwise keep the first two chunks. Nutritionally meaningful
  varieties survive: `beans, navy` / `rice, brown` / `mushrooms,
  shiitake` (plural) / `wheat, khorasan` / `cheese, cheddar`. Cheese
  family stays distinct: `cheese, parmesan` / `cheese, cottage` /
  `cheese, mozzarella` / `cheese, brie` / `cheese, camembert` / `cheese,
  edam` / `cream, sour` all remain separate canons (confirmed by user).

**Stage 3 — nutritional-similarity gate** (`refine_canons_by_nutrition`)
— for every multi-variant group whose head is **NOT** on the strip list,
fetch the variants' values on a key nutrient subset (protein, fat, carbs,
energy, total sugars, fiber, water, pentosan, starch, cellulose, inulin).
Compute the group median per nutrient. If any variant's amount deviates
by more than 3× from the median on any key nutrient, split it back out
as its own canon with the original FDC descriptive name. Protects
preserve-head groups (cheese / beans / mushrooms / etc.) from accidental
over-merge.

When the head IS on the strip list (`apple`, `nectarine`, `wild mango`,
`wheat flour`, etc.), the gate is **exempted** — the user's intent is
"collapse every variant regardless of nutritional differences". This
ensures ripe vs unripe wild mango both end up under `wild mango/african
mango/bush mango`, and 9 %-protein vs whole-wheat flour both end up
under `wheat flour`.

**Stage 4 — also stripped in Stage 1** (added in Phase 5 amendments 3 & 4):
- `_PAREN_CONTENT_RE` strips ALL parenthetical content (clarifications,
  color descriptors, source notes — `(industrial)`, `(colour of peel:
  olive green)`, `(fat free or skim)`, `(includes foods for USDA's food
  distribution program)`).
- `_QUANT_RE` strips quantitative qualifiers (`9% protein`, `50%
  extraction`, `3.25% milkfat`).
- `_PREP_RE` now also strips cooking methods (`fried`, `pan-fried`,
  `deep-fried`, `breaded`, `battered`, `oven`, `microwaved`), the FDC
  hedge `nfs` ("not further specified"), and dish descriptors `with
  batter`, `with sauce`.

**Stage 5 — smart cultivar-code detection in chunk[1]** (Phase 5
amendment 4) — for preserve-head foods, the second comma-chunk is
examined; if it matches a *cultivar code* pattern, it's dropped and the
next non-code chunk is used. Patterns:
- `[Latin name]` — `rice, [Oryza sativa]` → `rice`
- `'quoted cultivar'` / `"quoted cultivar"` — `lentil, 'CDC Blaze'` → `lentil`
- `var. X` / `cv. X` — `winged bean, var. SLS1` → `winged bean`
- Pure numeric — `rice, 3597` → `rice`
- Alphanumeric — `rice, ADT-21` → `rice`
- Breed cross with `X x Y` — `beef, aberdeen angus x holstein-friesian` → `beef`
- Numeric prefix codes — `pork, 1/2 duroc` → `pork`

Preserves meaningful varieties: `rice, brown` / `beef, sirloin` /
`cheese, cheddar` / `wheat, whole-grain` all stay distinct.

Both Stage 4 and Stage 5 are safe across the board: these patterns are
nearly always non-essential annotations.

### Examples

| raw FDC description | canonical name |
|---|---|
| `Carrots, sliced or crinkle cut, frozen, unprepared - Proximates - NF9913X5` | `carrots` |
| `Carrots, baby, raw` | `carrots, baby` |
| `Beans, navy, mature seeds, raw` | `beans, navy` |
| `Beans, Dry, Navy (0% moisture)` | `beans, navy` |
| `Strawberry, campineiro [BioFoodComp]` | `strawberry` |
| `Apple, Jonagold, fresh, raw [BioFoodComp]` | `apple` |
| `Nectarine, Arctic Star, raw [BioFoodComp]` | `nectarine` |
| `Lowbush blueberry, Blomidon, overripe, entirely blue [BioFoodComp]` | `lowbush blueberry` |
| `Grapefruit juice, white, bottled, unsweetened, ocean spray` | `grapefruit juice` |
| `Winged bean, var. SLS1, seeds, raw [BioFoodComp]` | `winged bean` |
| `Wild mango/African mango/bush mango, unripe, raw (colour of peel: olive green) [biofoodcomp]` | `wild mango/african mango/bush mango` |
| `Wild mango/African mango/bush mango, pulp, unripe, raw (colour of peel: olive green) [biofoodcomp]` | `wild mango/african mango/bush mango` |
| `Wheat flour, white (industrial), 9% protein, bleached, enriched` | `wheat flour` |
| `Milk, fluid, nonfat (fat free or skim)` | `milk, fluid, nonfat` (paren stripped, second chunk kept) |
| `Potato, white (industrial), 50% extraction` | `potato, white` |
| `Mushrooms, shiitake, raw` | `mushrooms, shiitake` (plural preserved as meaningful variety) |
| `Mushroom, 'Agaricus bisporus'` (singular + Latin) | `mushroom` (stripped — Latin codes are cultivar-flood) |
| `Cheese, cheddar, sharp, sliced` | `cheese, cheddar` |
| `Cheese, brie` | `cheese, brie` |
| `Cream, sour` | `cream, sour` |
| `Pasta [Phenol-Explorer]` | `pasta` (category back-filled to `Cereal Grains and Pasta`) |

Net effect on the FDC catalog (v16): **131 930 entries → 25 552 canonical
foods (5.16× merge ratio)**. 1 450 cultivar variants survived the
nutritional gate (only for non-strip heads — strip-list heads collapse
unconditionally). 144 strip-list-head groups exempted from the gate.
Smart cultivar-code detection trimmed `rice` from 321 → 252 canons,
`pork` from 310 → 294, `beef` from 261 → 256 by stripping `[Latin
name]`, quoted cultivars, `var. X`, numeric codes, and breed crosses
without touching meaningful varieties.

### Category back-fill

Several FDC source rows lack `food_category_id` — typically Phenol-Explorer,
AFCD, some BioFoodComp imports. The build pipeline applies a regex-based
keyword fallback (`_DESC_TO_CATEGORY`) covering 25 category-keyword groups
(Baby Foods, Fast Foods, Alcoholic Beverages, Soups+Sauces, Organ Meats,
Snacks, Pasta, Cereal Grains, Dairy + plant-based alternatives, Legumes,
Nuts+Seeds, all meat categories, all fish + shellfish species, Fruits +
juices, Vegetables + tubers + algae, Spices + condiments, Sweets +
sweeteners, Beverages, Fats+Oils). `Pasta [Phenol-Explorer]` → `Cereal
Grains and Pasta`; `Infant formula` → `Baby Foods`; `Pizza with pepperoni,
from restaurant` → `Fast Foods`; `Soup, minestrone` → `Soups, Sauces, and
Gravies`; `Plantains, green` → `Fruits and Fruit Juices`; etc.

**v18 category fill rate: 94.8 %** (89 % in v14, ~60 % before any
back-fill). The build emits `<index_dir>/uncategorized.txt` listing the
remaining ~1340 entries — mostly research-paper titles that leaked into
FDC imports, regional FDC jargon (Korean / Japanese / Scandinavian
abbreviations), and obscure plant compounds. When an output row's
`food_category` is blank, the resolved category from `food_stats` is
used at display time. Recent v18 additions: `quinoa`, `lasagne`,
`falafel`, `hummus`, `frankfurter`, `meatball`, `pike` / `turbot` /
`flounder` / `pangasius`, `celeriac`, `dandelion`, `collard greens`,
`stock cube`, `pesto`, `swiss roll`, `gingerbread`, `candybar`, `sushi`
/ `sashimi`, `horse meat`, `hare`, `meat alternative`.

### Cache invalidation when upgrading

The canonicalization is baked into `static_food_meta.pkl` and the modeled
index parquets. After pulling a new version of `bac2food_predict.py`,
delete the cache:

```bash
rm /data/bac2food/index_modeled/{static_food_meta.pkl,*.parquet}
```

The next predictor invocation rebuilds them (~2 min).


## Nutrient values: duplicate and conflicting rows

`nutrient_id` is the join key between food and enzyme (via
`0_building/3_nutrient_to_ec.tsv`); `fdc_id`, `data_type`, `food_category` and
`source_db` are bookkeeping. So what matters is how the predictor resolves a
`(food, nutrient_id)` pair that carries **more than one amount** — and it often
does. In the food store, 42,362 pairs have conflicting amounts:

- **17,168 conflict across databases** — two national DBs measured the same
  food differently.
- **25,194 conflict within a single database** — the source reports one nutrient
  twice with different values for one food (McCance gives canned ackee Niacin
  0.6 **and** 1.1 mg). Concentrated in AFCD, STFCJ (Japan), WAFCT, McCance. This
  is pre-existing in the ingested sources, not something the pipeline creates.
  `5_export/food_nutrients.tsv` exposes it via its `source_db` column.

**They do not double-count.** `build_modeled_index` aggregates every row for a
`(nutrient_id, canon)` with a **MEAN** — the same mean that folds variant
`fdc_id`s into their canonical food. Ackee's Niacin becomes 0.85 mg. Nothing is
summed, so a duplicated row cannot inflate a nutrient amount.

### The one place row multiplicity does leak: `model_count`

`model_count` counts **rows** per canonical food, and is the denominator of

```
prm = n_nutrients_improved / max(50, model_count)
```

which drives the breadth penalty `BROAD_W * (1 - prm) ** BROAD_Q`. Duplicate rows
inflate the denominator, so `prm` shrinks, the breadth penalty grows, and the
food's score drops slightly.

Two things keep this benign. The bias is **conservative** — it can only lower a
score, never inflate one. And it is second-order: `model_count` already counts
rows across *all* variant `fdc_id`s folded into a canon, so a 545-variant food
like oyster mushrooms already has a denominator in the thousands and `prm` is
near zero whatever the duplicates do.

Worth knowing when reading the code: `prm` is described as a purity ratio, but
`improved nutrients ÷ rows-across-variants` is not literally "the fraction of the
food's nutrients this bacterium improves". Counting distinct `nutrient_id`s per
canon would make it so — a deliberate modelling change that would shift every
score, so it has **not** been made.


## Bacteria → EC reference: cache + matching ladder

The reference is now `exports/species_enzymes.tsv` (**eggNOG v7**, built through KEGG KO
by `eggnog/6.1_eggnog7_species_enzymes.py` — v7 dropped the direct EC annotation that v6
carried). **v7 is now the only layer**: `exports/species_enzymes.v6.tsv` has been retired
from the deposit. The legacy source `/data/bac2food/bact_ec.tsv` is retained and still
readable — `_ref_read_options` detects the 4-column headerless v6 layout as well as the
6-column v7 one, so `--bact_ec_tsv` accepts either if you need to reproduce a v6 run.
The predictor:

1. Converts it to a deduplicated zstd parquet on first use; subsequent reads take
   < 1 s. The cache goes to `/data/bac2food/cache/` (override with `$BAC2FOOD_CACHE`),
   **not** next to the reference TSV — the reference lives in the export/deposit
   directory, and writing `<name>.parquet` beside it shipped a build artifact with the
   published resource. `5_export/verify_exports.py` now fails if anything other than the
   deliverable files appears in that directory.
2. Matches user species → reference species via a three-step ladder
   applied per user species:
   1. Exact `tax_id` match (when the user species has a leading
      `<digits>_` prefix).
   2. Exact normalized species-name match.
   3. Genus + species prefix match (strips strain trailing tokens).

The ladder is applied per user species so collisions don't drop anyone.

### Two different "match rates" — do not conflate them

`recompute_match_rate.py` measures both on the full demonstration cohort (55 samples,
1,257,938 annotated EC rows, 63 distinct species labels):

| metric | what it measures | eggNOG v6 | **eggNOG v7** |
|---|---|---|---|
| **Feature** | annotated loci reaching an FDC nutrient | 45.4 % | **45.4 %** |
| **Species** | cohort species resolving to a reference EC set | 14.3 % | **60.3 %** |
| **Species, row-weighted** | annotated rows belonging to a matched species | 13.5 % | **65.3 %** |

The **feature** rate is the one the paper quotes as "N % of annotated features match a
reference-map entry". It runs through EC → substrate → ChEBI → nutrient and never touches
the Bacteria → EC reference, so v6 → v7 leaves it unchanged (per-sample mean 45.4 %, range
43.4–46.3 %). Only `--augment_with_reference` would move it.

`chain_coverage.py` answers the follow-up question — *why* the other 54.6% never reach a
food. It classifies every EC by the stage at which the chain stops: no substrate (GAP A,
5.0% of cohort EC), no ChEBI id (GAP B, 1.5%), or no FDC nutrient (GAP C, **52.4%**). GAP C
dominates, and its members act on real molecules that simply are not food components, so the
match rate is bounded by the nutrient vocabulary rather than by BRENDA's coverage — a newer
substrate release would move it very little.

A locus may carry several EC numbers in ONE cell — 178,951 of the 1,257,938 do — and it
matches if any of them reaches a nutrient, which is what the predictor does after
`_normalize_ec_frame` explodes those cells. Counting the raw cell as a single EC understates
the rate by ~9 pp (36.3 % instead of 45.4 %); an earlier version of this table did exactly
that. Secondary views: exploded (locus, EC) pairs 43.9 %, distinct EC 41.1 %.

The **species** rate is what `--complement_ec` prints, and it is what v7 improves: 4× more
cohort species now resolve, because v7 carries 10,751 organisms against v6's 3,176 and uses
current NCBI names that match the cohort's MAG labels.

An earlier revision of this file quoted **34.4 %** as "the match rate" without saying
which. It reproduces neither metric against the current model — the model's EC set changed
across the v1→v18 iterations tabulated above — so it is superseded by the table.


## Validation — hep dataset, 8-bacterium panel, v14 differential top picks

| Bacterium | Specialty (literature) | v14 differential top-3 |
|---|---|---|
| *Streptococcus mutans* | sucrose, glucose, fructose (cariogenic) | nectarine / peach / apple — Chlorogenic acid + Epicatechin (the cariogenic phenolic profile) |
| *Bifidobacterium dentium* | HMOs, lactose | apple / peach / nectarine — phenolics + pectin |
| *Akkermansia* | mucin | wild mango / strawberry / wild mango unripe — Pectin/Pentoses + sialic-acid-rich cheese surfacing at #6 |
| *Veillonella tobetsuensis* | lactate, organic acids | **cream sour / cheese brie / cheese camembert** — `Lactose, Lactic acid`. Perfect biological match. |
| *Faecalibacterium* | inulin, FOS, glucose | apple / baobab fruit / wild mango — phenolics + Pectin |
| *Roseburia* | xylan, arabinoxylan | apple / peach / nectarine, then **milk cow / cheese cheddar** with `Galactose/Lactose/Sucrose` |
| *Bacteroides nordii* | broad polysaccharide | **flaxseed / ground flaxseed meal / peanut butter** — Pentosan/Pectin/Inulin |
| *Prevotella timonensis* | xylan, arabinoxylan | apple / cream sour / cheese cream — sugars + dairy |

### Aggregate progression across iterations

| metric | v1 | v11 | v12 | v13 | v14 | v15 | v16 | v17 | **v18** |
|---|---|---|---|---|---|---|---|---|---|
| top-1 max fraction (differential) | 78.7 % | 37 % | 23.2 % | 13.6 % | 26.4 % | 15.8 % | 15.8 % | 15.8 % | **11.9 % milk cow** |
| Shannon entropy top-1 (bits) | 1.54 | 2.79 | 4.29 | 4.70 | 4.23 | 4.62 | 4.62 | 4.62 | **4.73** |
| distinct top-1 foods | 29 | 38 | 70 | 64 | 63 | 71 | 71 | 71 | **70** |
| FDC → canonical foods | none | 47 633 | same | 29 239 | 27 585 | 26 220 | 25 552 | 25 552 | **25 545** |
| Category-fill rate | n/a | ~60 % | ~60 % | ~60 % | 89 % | 89 % | 89 % | 93.2 % | **94.8 %** |
| ChEBI-mapped nutrients (`3_nutrient_to_ec.tsv`) | 434 | 460 | 460 | 460 | 460 | 460 | 460 | **470** | 470 |
| Strip-list-head canons all collapsed to 1 (apple/nectarine/wild mango/wheat flour/tamarillo/terapy bean...) | ✗ | ✗ | ✗ | partial | partial | ✓ | ✓ | ✓ | **✓** |

v14's slight top-1 regression vs v13 reflects more aggressive cultivar
collapse — fewer distinct food rows means more bacteria converge on the
same top picks. This is acceptable because the differential `comp_score`
is recomputed on the new (smaller) food universe; apple's 26.4 % share
means apple is the *truly best differentiator* for that quarter of the
microbiome, not an artifact of redundant cultivar entries.


## Common workflows

### Default run

```bash
python bac2food_predict.py \
    --mag_tsv my_mags.tsv \
    --out_prefix my_run \
    --max_foods 10 \
    --augment_with_reference \
    --complement_ec \
    --jobs 6
```

Three files. Read **community** for diet design (what to feed the
microbiome). Read **differential** for personalization (what each
species is uniquely good at). Read **perFood** to see which species
each food upregulates.

### Broader bacterial universe (no user MAGs)

```bash
python bac2food_predict.py \
    --use_reference \
    --out_prefix ref_pool \
    --ref_min_ec 100 \
    --ref_max_species 500 \
    --top_bacteria_per_food 20 \
    --jobs 8
```

### Conservative augmentation

```bash
python bac2food_predict.py \
    --mag_tsv my_mags.tsv \
    --out_prefix my_run \
    --augment_with_reference --augment_threshold 100 \
    --jobs 6
```


## Substrate DB extensions (companion to the 0_building pipeline)

The substrate ontology lives in `0_building/` and feeds the predictor via
`3_nutrient_to_ec.tsv`. Three companion files extend it:

| file | purpose |
|---|---|
| `0_building/extra_bacterial_seeds.tsv` | 18 ChEBI seeds for substrates absent from FDC (HMOs, GlcNAc, sialic acid, fucose, xylan, arabinoxylan, dextran, cellobiose, pullulan, alginate, agarose) |
| `0_building/extra_nutrient_chebi.tsv` (Phase 6) | 77 hand-curated FDC nutrient → ChEBI mappings (sorbitol, ethanol, resistant starch, all tocopherols/tocotrienols, all carotenoids, all flavonoids/anthocyanins/lignans, organic acids, biogenic amines, etc.) — each row has a `justification` column |
| `0_building/inject_bacterial_substrates.py` | Writes synthetic food_nutrient rows so the novel substrate IDs actually drive food rankings (HMOs in human milk, sialic acid in dairy, fucose+alginate in seaweed, etc.) |

`3_nutrient_to_ec.py` accepts multiple `--extra_seeds` flags (repeatable),
so you can pass both extras files at regenerate time:

```bash
cd 0_building
python 3_nutrient_to_ec.py \
    --nutrient_best 2_nutrient_to_chebi.tsv \
    --extra_seeds extra_bacterial_seeds.tsv \
    --extra_seeds extra_nutrient_chebi.tsv \
    --digest_chebi ../chebi/digest_to_chebi.tsv \
    --chebi_obo ../chebi/chebi.obo \
    --out 3_nutrient_to_ec.tsv \
    --include_simple_sugars --max_cost 1.5
```

**Phase 6 ChEBI impact**: 460 → 470 distinct nutrient_ids with EC links
(+10), 8 870 → 9 234 rows (+364 EC links). The predictor consumes the
extended table without code changes. Soft aliases for GlcNAc / Cellobiose
/ Xylan / Arabinoxylan / Dextran / Pullulan live inside `main()` and
substitute at 25 % efficiency at scoring time.


## Phase 9 differential scoring (current default)

After Phase 8 the substrate graph grew from 470 → 778 nutrient IDs and the
food corpus picked up ~20 K v2 entries (Phenol-Explorer alone added 506
distinct polyphenols with EC mappings). That breadth flooded `score`'s
gain-sum: every bacterium that had any polyphenol-degrading EC scored
polyphenol-rich foods (peach, nectarine, chestnut) highly — `peer_median`
went up in lock-step, `comp_score` lost dynamic range, and specialty
substrate signals (FOS / HMO / mucin / arabinoxylan) got drowned.

**Phase 9 fix — per-bacterium specialty allowlist.** Each bacterium scores
only on its top-K most bacterially-specific nutrients (`spec(B, n) =
log(B_total / bf[n])` ranked desc). Generalist polyphenol signals are
suppressed; specialty substrates dominate.

### Phase 11 — the allowlist did the opposite of what it was for

The allowlist was *measured* against a 17-species graded panel
(`validate_biology_panel.py`, reciprocal rank of the first documented
substrate) and it fails there. Selecting on `log(B_total / bf[n])` ranks
nutrients by how **few** organisms carry them, but a growth substrate is by
definition one a guild **shares**. Across the panel, growth substrates
average `bf` 21.1 of 22 organisms while phytochemicals average 17.5, so:

* **4 of 308 allowlist slots (1.3 %)** held a growth substrate; 41.6 % were
  phytochemicals and 57.1 % other trace analytes.
* the best-ranked growth substrate for a given organism sat at position
  **56–86** of ~450 targets, against a cutoff of **14** — unreachable.
* the visible symptom: *B. thetaiotaomicron*, *ovatus*, *uniformis* and
  *fragilis* all returned the **same single food** (strawberry), ranked on
  the trace flavonol Morin.

`--specialty_mode`, `--spec_alpha`, `--specialty_topk`,
`--specialty_food_idf_weight` and `--specificity` were therefore **removed
outright** in Phase 11, together with the allowlist construction and the
bacterial-IDF factor in the per-nutrient weight. Both had already been
defaulted off, so the removal was verified byte-identical across all four
output files on a 22-species panel and a 63-organism cohort sample.
Reproduce the diagnosis with `diagnose_specialty.py`.

### Ranking rule

`comp_score = score − peer_median` removes any substrate the peers share —
again, exactly what defines a guild. Four starch degraders scored against
each other cancel on starch and get ranked on whatever trace compound
happens to differ. Treating the peer comparison as an admission test and
ranking on absolute `score` fixes that, and is the default.

Measured on the 17-species panel (MRR, and hits within the top 3):

| setting | MRR | hits@3 | whole-plant | junk/meat | neg. control |
|---|---|---|---|---|---|
| Phase 9 defaults | 0.216 | 4/17 | 94.9 % | 0 % | pass |
| **Phase 11 defaults** | **0.312** | **8/17** | **99.1 %** | 0 % | pass |

Three of the four *Bacteroides* now recover their **own** documented
substrate (pectin, hemicellulose, hemicellulose), and all four starch
organisms recover starch or inulin. The cost is convergence: mean pairwise
top-10 Jaccard across cohort organisms rises 0.08 → 0.21, because the old
rule maximised *apparent* distinctness by ranking on whatever differed —
including `Nitrites`, `CU(mg)` and `NA(mg)`, which are not substrates.

**Community output is unaffected**, structurally and not just empirically:
the community pass sets `community_dyn["specificity"] = False` and never
consults the allowlist. Verified byte-identical on three cohort samples.

`--differential_formula` stays at `full`. It is exactly equivalent to
`explicit_admission` while ranking on `comp_score` (the cost stack is a
food property that cancels in the subtraction — both score MRR 0.216), but
under absolute-score ranking nothing cancels and the stack becomes
load-bearing: `full` reaches hits@3 8/17 and passes the negative control,
against 6/17 and a rank-9 false positive.

### Shortlist diversity

A differential table is read as a shortlist to eat from, and the practical
shortlist is 5–10 items — few meals have more ingredients. At that length
redundancy is the dominant failure mode, and differential mode was the one
view with no defence against it: the community/greedy path applies
`OVERLAP_W × redundancy + category-repeat` while building its basket, but
differential is fed the **clean-slate** frame, which by construction has
neither. Measured on a cohort sample, the differential top-5 spanned
**2.08** distinct food categories against the greedy path's **3.71**.

`--differential_diversity` caps how many rows may share one food category
or one lead substrate. It defaults to `max(2, round(0.4 × --max_foods))`
— 2 at 5 foods, 4 at 10. Rows are consumed in score order, so within the
cap the ranking is untouched, and a starved list is backfilled rather than
shortened.

**The cap must scale with list length.** A fixed `K=2` is free at 5 foods
but punitive at 10:

| `--max_foods` | cap | MRR | hits@3 | whole-plant | junk |
|---|---|---|---|---|---|
| 5 | off | 0.296 | 8/17 | 99.1 % | 0 % |
| 5 | 2 | 0.294 | 7/17 | 99.1 % | 0 % |
| 10 | off | 0.312 | 8/17 | 99.1 % | 0 % |
| 10 | 2 *(too tight)* | 0.310 | 7/17 | **85.0 %** | 0.5 % |
| 10 | 4 *(auto)* | 0.315 | 8/17 | 98.2 % | 0 % |

Forcing ≥5 categories out of 10 rows reaches down into categories the
score had correctly declined; at 40 % of the list the cost disappears.

**Known cost.** The cap promotes lower-ranked foods, which is how a false
positive gets in: for *R. inulinivorans* it lifted lettuce, tomato and
broccoli to ranks 7–10, and lettuce carries a measured inulin value, so
the 18-species panel's negative control trips. The six-species panel the
paper reports still passes 3/3, including that negative. Pass
`--differential_diversity 0` for pure score order.

The remaining misses are a **data ceiling, not a scoring failure**. Of the
18 panel organisms, 8 have documented substrates that are EC-linked but
measured in fewer than 50 foods — arabinoxylan in 8, hemicellulose in 9,
sialic acid in 14, inulin in 33 — and a further group collapses onto a
single aggregate (`Starch`, 8 031 foods; `Lactose`, 4 753) that cannot
distinguish the organisms sharing it. Only polyphenols are resolved at
organism-discriminating granularity, which is why they dominated every
ranking. Reproduce with `substrate_ceiling.py`.

`test_specialty_panel.py` and `--validate_panel` are **gone** — four of
that panel's eight organisms carry zero rows in the shipped eggNOG v7
reference, so it could not score above 4/8, and on a cohort sample it
scored 0/8. Use
`validate_six_species_panel.py` (reproduces what the paper reports) and
`validate_biology_panel.py` (selects between formulas).

**Baseline → Phase 9 default → Phase 9 + substrate expansion (hep dataset):**

> **Superseded in Phase 11.** The sweep below tuned `--specialty_topk` and
> `--spec_alpha`, both of which have been removed along with the specialty
> allowlist and bacterial-IDF weighting. It is kept only as a record of how
> the defaults were once chosen. Note what the Phase 9 notes already said of
> the failures: *"Roseburia targets xylan / arabinoxylan / beta-glucan, but
> those substrates are common across many gut bacteria (low bacterial-IDF)
> so they don't land in K=14."* That is the whole defect, recorded at the
> time and treated as out of scope: selecting nutrients by rarity across
> organisms cannot surface a substrate that a guild shares. Phase 11 deleted
> the mechanism instead of retuning it, and Roseburia now recovers
> beta-glucan.



## Known limitations

- **Annotation gaps are the dominant failure mode** for under-annotated
  MAGs whose taxon isn't represented in `bact_ec.tsv` (no reference
  match). Workaround: re-annotate with eggNOG-mapper or dbCAN3, or add
  the missing strain to the reference TSV.
- **Synthetic food rows are literature-based estimates**, not measurements
  for the specific FDC food. Concentrations carry ~2× uncertainty.
- **Branded foods**: FDC `data_type=branded_food` rows are dropped by
  default (`DROP_BRANDED=True`). Some all-caps branded variants slip
  past this filter because they're nominally `survey_fndds_food`; MAX
  aggregation on the canonical foods keeps them harmless.
- **Default paths point at `/data/bac2food/`**. Override every path flag
  if your data lives elsewhere.


## Files in this folder

| file | purpose |
|---|---|
| `bac2food_predict.py` | the predictor |
| `parameters.yaml` | tunable constants — scoring weights, data paths, tables (see **Configuration**) |
| `README.md` | this document |
| `readme.txt` | short CLI quick-reference (legacy) |
