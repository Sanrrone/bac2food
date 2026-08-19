#!/usr/bin/env python3
"""bac2food_predict.py - Bacteria <-> Food predictor for metagenomes.

Two modes:
  --mode bacteria2food : for each bacterium, rank foods it can use
  --mode food2bacteria : for each food, rank bacteria that benefit

Optional --complement_ec : diff user EC set vs /data/bac2food/bact_ec.tsv reference
to flag enzymes likely present in the genome but missed by the annotator.

Scoring kernel, static-food-meta builder, and modeled-index builder are adapted from
the single-organism query prototype this tool supersedes (same math, same parquet
layout). That prototype is not part of the release.
"""
from __future__ import annotations
import argparse, math, os, re, gc, pickle, heapq, sys, zlib, shutil, hashlib
from pathlib import Path
import multiprocessing as mp
import pandas as pd
import numpy as np
import pyarrow.dataset as ds
import pyarrow as pa
import pyarrow.csv as pacsv
try: import resource
except ImportError: resource = None

# CONFIGURATION
# All constant, non-user-selectable parameters live in `parameters.yaml` next to
# this script (override with $BAC2FOOD_PARAMS or --config). They are loaded at
# IMPORT time so that multiprocessing workers - which re-import this module under
# the "forkserver"/"spawn" start methods - bind exactly the same values.
#
# Logic / structural constants that are NOT configuration stay hardcoded here:
#   * EC_RE     - the EC-number validation regex
#   * BUCKETS   - tied to the on-disk bucketed-parquet layout (do not change)
#   * the text regexes and cultivar / category word-lists further down.
import yaml

EC_RE   = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
BUCKETS = 256

# Names bound from parameters.yaml by apply_params(); declared up-front so module
# code and `global` statements resolve them. Real values are set at import below.
PATH_NUTRIENT_TO_EC = PATH_NUTRIENT_CSV = PATH_FOOD_CSV = PATH_FOOD_CATEGORY_CSV = None
PATH_FOOD_PORTION_CSV = PATH_BUCKETED_DIR = PATH_INDEX_DIR = PATH_NUTRIENT_ALIAS = PATH_BACT_EC_REF = None
DROP_BRANDED = DROP_MODELLED = ALLOW_MACRO_PROXY = ALLOW_MACRO_SCAN = None
PROC_W = BROAD_W = BROAD_Q = TYPE_W = ART_W = None
# Guards, not tunables. Each was perturbed by a large factor against both the 22-species
# panel and a 63-organism cohort sample and changed not one output byte, so they are
# constants rather than knobs: they bound pathological input without shaping any result
# we produce. GAIN_MIN drops candidates with no meaningful gain, SCORE_MIN floors the
# score, MAX_CANONS caps the candidate set so a larger food store cannot exhaust memory.
GAIN_MIN   = 0.001
SCORE_MIN  = -2.0
MAX_CANONS = 250_000
K_PLANT = K_OTHER = None
PROC_BASE_W = JUNK_BASE_W = ALC_BASE_W = None
OVERLAP_W = FIBER_WEIGHT = None
OLIGO_IDS = ISOFLAVONE_IDS = None
CAT_PENALTY = PLANT_CATS = TP_MAP = DANGEROUS_RULES = SAFE_RULE_PFXS = ALWAYS_DROP_CATS = None
_KEY_NUTRIENT_IDS = _OUTLIER_RATIO = None


def _resolve_cfg_path(p, base: Path):
    """Resolve a config path: absolute kept as-is; relative resolved against the
    parameters.yaml directory (so paths are independent of the working dir)."""
    if p is None:
        return None
    pp = Path(p)
    return str(pp if pp.is_absolute() else (base / pp))


def apply_params(cfg: dict, base: Path) -> None:
    """Bind every externalized module global from a parsed parameters.yaml dict."""
    g = globals()
    P = cfg["paths"]
    g["PATH_NUTRIENT_TO_EC"]    = _resolve_cfg_path(P["nutrient_to_ec"], base)
    g["PATH_NUTRIENT_CSV"]      = _resolve_cfg_path(P["nutrient_csv"], base)
    g["PATH_FOOD_CSV"]          = _resolve_cfg_path(P["food"], base)
    g["PATH_FOOD_CATEGORY_CSV"] = _resolve_cfg_path(P["food_category"], base)
    g["PATH_FOOD_PORTION_CSV"]  = _resolve_cfg_path(P["food_portion"], base)
    g["PATH_BUCKETED_DIR"]      = _resolve_cfg_path(P["bucketed_dir"], base)
    g["PATH_INDEX_DIR"]         = _resolve_cfg_path(P["index_dir"], base)
    g["PATH_NUTRIENT_ALIAS"]    = _resolve_cfg_path(P["nutrient_alias"], base)
    g["PATH_BACT_EC_REF"]       = _resolve_cfg_path(P["bact_ec_ref"], base)

    F = cfg["flags"]
    g["DROP_BRANDED"]      = bool(F["drop_branded"])
    # Absent from an older parameters.yaml, so default to dropping: a config that predates
    # this flag should get the current food set, not silently score modelled FNDDS rows.
    g["DROP_MODELLED"]     = bool(F.get("drop_modelled", True))
    g["ALLOW_MACRO_PROXY"] = bool(F["allow_macro_proxy"])
    g["ALLOW_MACRO_SCAN"]  = bool(F["allow_macro_scan"])

    S = cfg["scoring"]
    g["PROC_W"]              = float(S["proc_w"])
    g["BROAD_W"]             = float(S["broad_w"])
    g["BROAD_Q"]             = float(S["broad_q"])
    g["TYPE_W"]              = float(S["type_w"])
    g["ART_W"]               = float(S["art_w"])
    g["K_PLANT"]             = int(S["k_plant"])
    g["K_OTHER"]             = int(S["k_other"])
    g["PROC_BASE_W"]         = float(S["proc_base_w"])
    g["JUNK_BASE_W"]         = float(S["junk_base_w"])
    g["ALC_BASE_W"]          = float(S["alc_base_w"])
    g["OVERLAP_W"]           = float(S["overlap_w"])
    g["FIBER_WEIGHT"]        = float(S["fiber_weight"])
    g["_OUTLIER_RATIO"]      = float(S["outlier_ratio"])

    N = cfg["nutrient_ids"]
    g["OLIGO_IDS"]         = set(int(x) for x in N["oligo_ids"])
    g["ISOFLAVONE_IDS"]    = set(int(x) for x in N["isoflavone_ids"])
    g["_KEY_NUTRIENT_IDS"] = tuple(int(x) for x in N["key_nutrient_ids"])

    T = cfg["tables"]
    g["CAT_PENALTY"]      = {str(k): float(v) for k, v in T["cat_penalty"].items()}
    g["PLANT_CATS"]       = set(T["plant_cats"])
    g["ALWAYS_DROP_CATS"] = set(T["always_drop_cats"])
    g["TP_MAP"]           = {str(k): int(v) for k, v in T["tp_map"].items()}
    g["DANGEROUS_RULES"]  = set(T["dangerous_rules"])
    g["SAFE_RULE_PFXS"]   = set(T["safe_rule_pfxs"])


def load_params(path) -> None:
    """Read parameters.yaml at `path` and bind all module globals from it."""
    path = Path(path)
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    apply_params(cfg, path.resolve().parent)


# Default config = parameters.yaml next to this script, or $BAC2FOOD_PARAMS.
# main() re-loads this if the user passes --config.
_PARAMS_PATH = Path(os.environ.get("BAC2FOOD_PARAMS",
                                   str(Path(__file__).with_name("parameters.yaml"))))
load_params(_PARAMS_PATH)

ALCOHOL_RE   = re.compile(r"\b(?:vodka|gin|rum|whiskey|whisky|bourbon|scotch|brandy|cognac|tequila|mezcal|liqueur|beer|lager|ale|cider|wine|champagne|prosecco|vermouth|sake|martini|margarita|mojito|daiquiri|pina colada|long island|bloody mary|cosmopolitan|mimosa|sangria|\balcohol\b|\bethanol\b)\b", re.I)
JUNK_RE      = re.compile(r"\b(?:pastry|fries|frosting|cand(?:y|ies)|syrup|concentrate|dessert|cookie|cake|pie|doughnut|soda|sugar|sweetened|candied|glaze|marshmallow|ice cream|pudding|whip(?:ped)?|topping|chocolate|caramel|gumm(?:y|ies)|crisps?|chips?|puffs?)\b", re.I)
POWDER_RE    = re.compile(r"\b(?:freeze-dried|dehydrated|powder|extract|dried|instant|dry mix|beverage mix)\b", re.I)
PROCESSED_RE = re.compile(r"\b(?:deep[-\s]?fried|pan[-\s]?fried|fried|hash brown|breaded|battered)\b", re.I)
_FIBER       = re.compile(r"\b(?:fiber|starch|oligosaccharide|inulin|pectin|glucan|hemicellulose|cellulose|pentosan)\b", re.I)
NFY_RE       = re.compile(r"(?:^\s*(?:sugars|fatty acids|amino acids|vitamin [a-z0-9 ]+|minerals|organic acids)\s*,)|(?:\s-\s*NFY)")
WHALE_RE     = re.compile(r"whale|seal|muktuk|walrus|blubber|bowhead", re.I)
# ALWAYS_DROP_CATS is configured in parameters.yaml (tables.always_drop_cats).

# Nutrients dropped from consideration regardless of DB content. Designed for
# bacterial-substrate use: simple sugars (sucrose / glucose / fructose / lactose /
# maltose / galactose / mannitol / sorbitol) and free amino acids are
# intentionally KEPT here so that 0_building/3_nutrient_to_ec.tsv built with
# --include_simple_sugars / --include_amino_acids actually reaches the scorer.
# What we still block: vitamins, minerals, energy/proximate aggregates, fatty
# acids and lipid macros, processing artifacts.
_BLOCKED = re.compile(
    r"\b(?:sugar|cholesterol|sfa|mufa|pufa|tfa|trans|fatty acid|fatty acids|sterol|dha|epa|docosahexaenoic|eicosapentaenoic|polyunsaturated|monounsaturated|saturated|triglyceride|lipid|lipids|phospholipids?|glycolipids?|unsaponifiable|trans-monoenoic|trans-dienoic|trans-polyenoic|vitamin\s*[abcdefjkh][-\s\d]*|vitamin\s*b[-\s]?(?:[1-9]|12)|vitamins\b|b[-\s]?2\b|b[-\s]?6\b|b[-\s]?12\b|riboflavin|thiamin|niacin|folate|folic|tetrahydrofolate|5-mthf|pantothenic|biotin|pyridoxin|cobalamin|retinol|retinoic|carotene|tocopherols?|tocotrienols?|total\s+tocopherols?|ascorbic|menaquinone|phylloquinone|calcium|iron|magnesium|phosphorus|potassium|sodium|zinc|copper|selenium|manganese|iodine|choline|betaine|chromium|fluoride|molybdenum|c\d{1,2}:\d+|nacl|salt|salt deklaration|salt labelling|niacin[æa]kvivalent|niacine, pr[ée]form[ée]e|polyols|tannins?|flavans?|flavones?|flavonols?|anthocyanins?|haem|nhaem|hydroxyproline|glycogen|tyramin|putrescin|chitin|polyphenol|nutrient_name)\b|(?:aa \(pernitrogen\)|starch-|pp-|polymers \(>10 mers\))",
    re.I
)

# Food-name canonicalization. Strips two layers of noise:
#  (1) the "- <nutrient panel> - NF<hex>" / ", <panel> -NF<hex>" fragmentation
#      suffix the FDC catalog uses to split the same food into per-panel rows
#      ("Carrots, ... unprepared - Proximates - NF9913X5" + "... - Vitamin E -NF...")
#  (2) preparation / state / cut qualifiers that don't change the metabolic
#      profile from a bacterium's POV (raw, cooked, frozen, sliced, drained,
#      with/without salt, whole, halves, in oil, etc.). "Carrots, sliced or
#      crinkle cut, frozen, unprepared" -> "carrots".
# All fdc_id rows that canonicalize to the same name fold into a single canon,
# their nutrients union'd via MAX in build_modeled_index.
_NF_SUFFIX_RE = re.compile(r"\s*[-,]\s*[A-Z][\w &/,()-]*?\s*-?\s*\bNF\w+\s*$", re.I)
_NF_BARE_RE   = re.compile(r"\s*[-,]\s*\bNF\w+\s*$", re.I)
_PREP_RE      = re.compile(
    r"\b(?:"
    r"raw|cooked|uncooked|boiled|baked|broiled|grilled|roasted|toasted|stewed|steamed|saut[ée]ed|simmered|"
    r"fried|deep[-\s]?fried|pan[-\s]?fried|stir[-\s]?fried|battered|breaded|"
    r"oven|microwaved|"
    r"frozen|fresh|dried|dehydrated|reconstituted|canned|drained|undrained|dry|"
    r"unprepared|prepared|reheated|chilled|shelf\s+stable|nfs|"
    r"cut|sliced|chopped|diced|crinkle\s*cut|julienned|grated|shredded|cubed|halved|quartered|"
    r"halves|quarters|pieces|chunks|wedges|spears|florets|kernels?|with\s+batter|with\s+sauce|"
    r"mature\s+seeds?|seeds?|young|mature|immature|ripe|unripe|overripe|"
    r"with(?:out)?\s+salt|salted|unsalted|"
    r"with(?:out)?\s+skin|peeled|unpeeled|with(?:out)?\s+bones?|boneless|skinless|"
    r"organic|conventional|enriched|fortified|"
    r"ready[-\s]to[-\s]eat|"
    r"in\s+oil|in\s+water|in\s+juice|in\s+syrup|in\s+brine|"
    r"0%\s*moisture"
    r")\b",
    re.I,
)
_EMPTY_PARENS_RE = re.compile(r"\s*\(\s*\)\s*")
_BRACKET_TAG_RE = re.compile(r"\s*\[[^\]]+\]\s*$")  # trailing [BioFoodComp], [AFCD], etc.
# Variant tags anywhere in the description (after source-tag strip): Phenol-
# Explorer and BioFoodComp routinely encode color / cultivar / type variants
# in brackets: "Common cabbage [Purple]", "Orange [Blond] / [Blood]",
# "Tea [Oolong] / [Green] / [Black]", "Common bean [Black] / [Others]",
# "Olive [Black]", "Lettuce [Red] / [Green]", "Swiss chard leaves [Red]",
# "Apple [Cider] / [Dessert]". None of these are identity-defining for the
# bacterial-food predictor; strip them so the canonicalizer collapses cousins.
_BRACKET_VARIANT_RE = re.compile(r"\s*\[[^\]]+\]\s*")
# FDC Foundation Foods sometimes prefix descriptions with their nutrient-panel
# group label ("Proximates, Beef, Eye of Round...", "Beverages, ABBOTT, ENSURE
# PLUS...", "Cereals ready-to-eat, ALPEN"). The real food name lives AFTER
# this prefix; without stripping it, ~360 canon names get mis-headed under
# 'proximates' / 'beverages' / 'cereals' instead of the actual food head.
# Only fires when a second comma chunk follows, so legitimate "Beverages"
# bare entries (no following comma) survive.
_FDC_GROUP_PREFIX_RE = re.compile(
    r"^\s*(?:proximates|beverages|cereals\s+ready[-\s]to[-\s]eat|cereals)\s*,\s*",
    re.I,
)
# Research-paper-title heuristic: descriptions > 120 chars containing any of
# these journal-style keywords. Used in build_static_food_meta to skip the
# ~50 BioFoodComp rows where the source description column leaked the
# literature reference instead of the food name.
_PAPER_TITLE_RE = re.compile(
    r"\b(characterization|determination|extraction|analytical\s+method|"
    r"comprehensive|systematic\s+review|meta-?analysis|isotope\s+analysis|"
    r"hplc-?ms|hplc/ms|gc-?ms|nmr\s+spectroscopy|"
    r"phytochemical\s+screening|nutritional\s+composition\s+study|"
    r"bioactive\s+compounds\s+of|polyphenol\s+content\s+of|"
    r"flavonoid\s+content\s+of|chemical\s+composition\s+of|"
    r"adequate\s+vitamin|status\s+of\s+norwegian|amino\s+acids\s+and\s+other|"
    r"alternative\s+to\s+the\s+iso|solid\s+phase\s+extraction-based)\b",
    re.I,
)
# Strip ALL parenthetical content (clarifications, color descriptors, source
# notes - "(industrial)", "(colour of peel: olive green)", "(fat free or
# skim)", "(includes foods for USDA's food distribution program)"). Parens
# in FDC descriptions are nearly always non-essential annotations, never
# identity-defining.
_PAREN_CONTENT_RE = re.compile(r"\s*\([^)]*\)\s*")
# Strip quantitative qualifiers like "9% protein", "50% extraction",
# "3.25% milkfat". Limited trailing word run so we don't over-eat.
_QUANT_RE = re.compile(r"\b\d+(?:\.\d+)?\s*%\s*[A-Za-z][A-Za-z\s-]{0,25}?(?=,|$)")
# Detect cultivar-code-like chunks for preserve heads. When chunk[1]
# matches one of these patterns, it's a research-level cultivar
# identifier rather than a nutritionally meaningful variety; we drop it
# and look for the next non-code chunk to keep. This catches the rice /
# pork / beef / lentil cultivar-code bloat without forcing those heads
# onto the strip list (which would also kill meaningful varieties like
# rice, brown / cheese, cheddar / mushrooms, shiitake).
_CULTIVAR_CODE_RE = re.compile(
    r"^("
    r"\[[^\]]+\]"                                    # [Oryza sativa]
    r"|'[^']+'"                                      # 'CDC Blaze'
    r"|\"[^\"]+\""                                   # "Gyokuro"
    r"|(?:var|cv)\.?\s+.{1,40}"                      # var. SLS1, cv. Foo
    r"|\d+(?:[-/.]\d+)*"                             # 3597, 80/20, 1.5
    r"|[a-z]{1,5}[-\s]?\d+(?:[-\s]?\d+)*"            # adt-21, rd 6, b-3428
    r"|[a-z\s-]+x\s+[a-z\s-]+"                       # angus x holstein-friesian (hyphen support)
    r"|\d+/\d+\s+\w+"                                # 1/2 duroc, 80/20 trimming
    r"|\*"                                           # bare asterisk markers
    r")$",
    re.I,
)

# Heads where cultivar / variety / color / brand info is metabolically
# irrelevant - all cultivars of strawberry behave the same to a bacterium.
# When a canonicalized name's first comma-separated chunk is in this set,
# everything after the first comma is dropped. For all other heads, we
# preserve up to the first two chunks (so "beans, navy" / "mushrooms,
# shiitake" / "wheat, khorasan" / "rice, brown" stay distinct).
_CULTIVAR_STRIP_HEADS = frozenset({
    # Tree fruits
    "apple", "apples", "pear", "pears", "peach", "peaches",
    "plum", "plums", "cherry", "cherries", "apricot", "apricots",
    "fig", "figs", "date", "dates", "persimmon",
    # Citrus
    "orange", "oranges", "tangerine", "mandarin",
    "lemon", "lemons", "lime", "limes", "grapefruit", "pomelo",
    # Tropical
    "mango", "mangos", "mangoes", "papaya", "pineapple",
    "banana", "bananas", "kiwifruit", "kiwi", "pomegranate",
    "guava", "passion fruit", "tamarind",
    # Grapes & melons
    "grape", "grapes", "watermelon", "cantaloupe", "honeydew",
    "melon", "melons",
    # All berries
    "strawberry", "strawberries",
    "blueberry", "blueberries",
    "lowbush blueberry", "lowbush blueberries",
    "highbush blueberry", "highbush blueberries",
    "southern highbush blueberry",
    "raspberry", "raspberries", "blackberry", "blackberries",
    "cranberry", "cranberries", "boysenberry", "boysenberries",
    "loganberry", "loganberries", "gooseberry", "gooseberries",
    "currant", "currants", "elderberry", "elderberries",
    "mulberry", "mulberries", "bilberry", "lingonberry",
    "bigney berry", "rabiteye", "monkey orange", "miracle fruit",
    # Juices
    "grapefruit juice", "orange juice", "apple juice", "grape juice",
    "pineapple juice", "cranberry juice", "vegetable juice",
    "tomato juice", "carrot juice", "lemonade",
    # BioFoodComp species-level entries with many cultivar rows
    "quinoa", "common bean", "common beans",
    "faba bean", "faba beans", "cowpea", "cowpeas",
    "chickpea", "chickpeas", "cassava", "rice bran",
    # Leafy greens & fruit-vegetables
    "lettuce", "spinach", "kale", "chard",
    "tomato", "tomatoes", "cucumber", "cucumbers",
    "eggplant", "eggplants", "aubergine",
    # Heads added in Phase 5 after a systematic scan of food.parquet:
    # BioFoodComp / PhyFoodComp cultivar codes dominate these and meaningful
    # varieties (if any) survive via the nutritional-similarity gate.
    "nectarine", "nectarines",
    "potato tuber",
    "pea",                            # singular: cultivar floods. Plural "peas" preserved.
    "lentil",                         # singular: BioFoodComp cultivar codes. Plural "lentils" preserved.
    "bean",                           # singular: Latin botanical names. Plural "beans" preserved (navy / kidney / pinto stay).
    "lima bean", "lima beans",
    "mung bean", "mung beans",
    "winged bean", "winged beans",
    "pigeon pea", "pigeon peas",
    "rice bean", "rice beans",
    "lupin", "lupine", "lupins",
    "soybean", "soybeans",
    "pearl millet",
    "date palm",
    "mushroom",                       # singular: Latin species names ("agaricus bisporus"). Plural "mushrooms" preserved.
    "wheat flour",                    # extraction percentages
    "wild mango/african mango/bush mango",
    "baobab fruit /moneky bread",
    "anchote",
    "locust bean", "african locust bean",
    "velvet beans", "velvet bean",
    "mcdonald's",                     # branded item list
    # Phase 5 amendment 4 - surfaced by the build-time under-grouped
    # diagnostic on v15. All have ≥ 5 sub-canons that are almost
    # entirely cultivar codes / breed identifiers / [Latin name] tags.
    "common vetch",                   # IFLVS cultivar codes
    "field pea",                      # quoted cultivar names ('Fidelia', etc.)
    "amaranth globe",                 # ABS / AKS / TVNU codes
    "amaranth",                       # related
    "potato tubers",                  # cultivars (granola, solara, marabel)
    "moth bean",                      # var. codes
    "v. vexillata",                   # TVNU codes (Vigna vexillata)
    "northern highbush blueberry",    # cultivars (bluecrop, jersey, croatan)
    "grasshopper",                    # insect food, sex/state variants
    "wattle",                         # Australian Acacia varieties
    "quince",                         # cultivar codes mixed
    # Phase 6 amendment - surfaced by user spot-check (color/skin/varietal flood,
    # nutritionally equivalent variants):
    "tamarillo", "tamarillo fruit",   # purplish-red vs golden-yellow skin
    "terapy bean",                    # wild / white / brown color variants
    "japanese horse chestnut",        # cultivar floods
    "japanese spikenard",             # cultivation method variants
    # Phase 9 amendment - user spot-check:
    "horse eye bean",                 # testa color variants (dark/light/black brown)
    "rice",                           # user explicitly wants all cultivars collapsed
                                      # (njavara, PTB 39, IR 64, MR 219, etc.).
                                      # Tradeoff: loses brown vs white distinction;
                                      # the user accepts this given the cultivar flood.
})

# Phase 8 amendment - heads that have cultivar / breed floods from
# BioFoodComp + PhyFoodComp + STFCJ (rice cultivar codes like Anjung/Calrose;
# pork/beef/veal breed names like Alentejano/Aberdeen Angus/Iberian-Duroc).
# Unlike the heads in _CULTIVAR_STRIP_HEADS_HARD above, these get the
# nutritional-similarity gate APPLIED so that genuinely distinct variants
# (rice brown vs white vs jasmine; Iberian pork vs commercial pork) split off
# as their own canon. The "trivial" 216 rice cultivars + 402 pig/cattle/veal
# breeds collapse to the bare head; the ~40 legit rice varieties and ~257
# legit meat cuts keep their identity automatically.
_CULTIVAR_STRIP_HEADS_GATED = frozenset({
    # Rice was here in Phase 8 but the gate was splitting 60+ cultivars off
    # because BioFoodComp dehusked/red/whole-grain rice variants legitimately
    # differ on protein/fat from FDC's milled white rice baseline (5-10×).
    # User explicitly wants these collapsed (njavara, PTB 39, IR 64, etc.),
    # so rice moved to HARD strip below.
    "pork",                           # 402 breeds vs 257 cuts
    "beef",                           # similar mix
    "veal",                           # similar mix
    "pasta",                          # Phase 10: verbose McCance recipe descriptors
                                      # ("pasta, egg, fresh, filled with..., boiled in...")
                                      # collapse. Gate splits genuinely distinct pasta types
                                      # (whole-wheat 2-3× fiber vs white) automatically.
    # Phase 8 amendment 5: yoghurt flavor variants. All flavored yoghurts of
    # the same fat content have essentially the same bacterial-substrate
    # profile (lactose + casein + small added-sugar delta). The gate splits
    # off genuinely distinct variants (greek high-protein, nonfat zero-fat,
    # whole milk full-fat) automatically.
    "yoghurt", "yogurt", "yogourt",
})

# Hard strip heads = the gate-exempt ones above. The combined set
# _CULTIVAR_STRIP_HEADS is used by the canonicalize_food_name head-strip
# step; the gate-exemption check uses _CULTIVAR_STRIP_HEADS_HARD only.
_CULTIVAR_STRIP_HEADS_HARD = _CULTIVAR_STRIP_HEADS
_CULTIVAR_STRIP_HEADS = _CULTIVAR_STRIP_HEADS_HARD | _CULTIVAR_STRIP_HEADS_GATED

# Description → category fallback. Some non-FDC source rows (Phenol-Explorer,
# AFCD, certain BioFoodComp imports) come without food_category_id. The
# existing source-specific branches in build_static_food_meta cover only a
# few cases (and the [phenol-explorer] branch defaults to "Fruits and Fruit
# Juices" - which is wrong for `Pasta [Phenol-Explorer]`). This list runs
# AFTER the source-specific branches as a last-resort keyword fallback.
# Ordered: first matching regex wins.
_DESC_TO_CATEGORY = [
    # === Order matters: first match wins. Specific categories go BEFORE
    # generic ones so e.g. "infant formula" hits Baby Foods, not Dairy. ===
    # Baby foods - very specific
    (re.compile(r"\b(infant\s+formula|toddler\s+formula|babyfood|baby\s+food|baby\s+toddler|baby\s+cereal|baby\s+mum\s+mum)\b", re.I), "Baby Foods"),
    # Fast foods - pizzas, burgers, sandwiches, etc.
    (re.compile(r"\b(pizza|pizzas|hamburger|hamburgers|cheeseburger|burger|sandwich|sandwiches|hot\s+dog|corn\s+dog|taco|tacos|burrito|burritos|quesadilla|quesadillas|nachos|nugget|nuggets|fries|french\s+fries|french\s+toast|french\s+dip|meatball\s+sub|wrap|wraps|sub|hoagie|gyro|kebab|kabob|shawarma|slider|sliders|panini|club\s+sandwich|grilled\s+cheese)\b", re.I), "Fast Foods"),
    # Alcoholic beverages - pull out before generic Beverages
    (re.compile(r"\b(wine|beer|ale|lager|stout|porter|spirit|liqueur|vodka|whiskey|whisky|rum|gin|tequila|mezcal|cocktail|black\s+russian|martini|margarita|mojito|cosmopolitan|mimosa|sangria|sake|absinthe|brandy|cognac|champagne|prosecco|vermouth|cider|alcoholic)\b", re.I), "Alcoholic Beverages"),
    # Soups, sauces, gravies - extended with stock + sauces
    (re.compile(r"\b(soup|soups|stew|stews|chowder|bisque|gazpacho|minestrone|broth|consomm[ée]|gumbo|chili|bouillon|miso\s+soup|ramen\s+broth|stock|stock\s+cube|stock\s+from\s+cube|stock\s+gel|pesto|curry\s+sauce|curry\s+paste|hoisin|teriyaki|tomato\s+sauce|tomato\s+pur[ée]e|marinara|alfredo|bolognese|carbonara|tomato\s+coulis|coulis)\b", re.I), "Soups, Sauces, and Gravies"),
    # Organ meats / luncheon meats - before generic beef/pork.
    # Also catches the generic "meat", "meatball", "frankfurter", "wiener",
    # "minced meat", "vegetarian fillet/meat loaf" so they don't fall to
    # the more specific beef/pork/poultry regexes.
    (re.compile(r"\b(liver|kidney|heart|tongue|tripe|intestine|gizzard|brain|sweetbread|offal|organ\s+meat|deli\s+meat|deli-meat|p[âa]t[ée]|pate|haggis|head\s+cheese|meat|meatball|meatballs|meatloaf|meat\s+loaf|wiener|wieners|frankfurter|frankfurters|hot\s+dog|salami|minced\s+meat|vegetarian\s+meat|vegetarian\s+fillet|meat\s+alternative|meatless|jerky|biltong|pastrami|corned\s+beef|brisket|terrine)\b", re.I), "Sausages and Luncheon Meats"),
    # Snacks (popcorn, chips, crisps, etc.) - before cereal generic
    (re.compile(r"\b(popcorn|corn\s+puff|corn\s+chip|chips?|crisps?|pretzel|pretzels|cracker\s+snack|trail\s+mix|granola\s+bar|protein\s+bar|meal\s+replacement\s+bar|cereal\s+bar|rice\s+cake|veggie\s+chip|potato\s+crisps?|potato\s+chips?|tortilla\s+chips?)\b", re.I), "Snacks"),
    # Pasta-specific (priority before generic "Cereal Grains and Pasta")
    (re.compile(r"\b(pasta|spaghetti|macaroni|noodle|noodles|ramen|udon|soba|lasagna|fettuccine|penne|rotini|ziti|orzo|gnocchi|couscous|bulgur|farro|risotto|tortellini|ravioli|vermicelli|linguine|rigatoni|cannelloni|tagliatelle)\b", re.I), "Cereal Grains and Pasta"),
    # Cereal Grains - extended with quinoa, lasagne, more pseudo-cereals
    (re.compile(r"\b(rice|wheat|oat|barley|rye|sorghum|millet|teff|amaranth|spelt|kamut|einkorn|emmer|buckwheat|maize|grain|grains|bran|flour|bread|cracker|crackers|biscuits?|cereal|cereals|tortilla|cornmeal|polenta|grits|pancake|waffle|muffin|bagel|croissant|brioche|baguette|naan|pita|focaccia|crispbread|crisp\s+bread|porridge|muesli|granola|oatmeal|sago|arrowroot|fonio|foniopaddy|hominy|tapioca|matzo|matzah|piroshki|quinoa|lasagne|lasagna|toast|melba|rusk|rusks|breadstick|pretzel|pretzels|popcorn|crouton|pretzel|grits|sourdough|gnocchi|dumplings?|flatbread)\b", re.I), "Cereal Grains and Pasta"),
    # Dairy - plant-based alternatives + traditional dairy
    (re.compile(r"\b(plant-?based\s+alternative\s+to\s+dairy|plant-?based\s+milk|plant-?based\s+cheese|vegan\s+cheese|vegan\s+milk|soymilk|soy\s+milk|oat\s+milk|almond\s+milk|coconut\s+milk|nut\s+milk|rice\s+milk|cashew\s+milk|hemp\s+milk)\b", re.I), "Dairy and Egg Products"),
    (re.compile(r"\b(milk|cheese|yogurt|yoghurt|butter|cream|kefir|buttermilk|whey|curd|dairy|ghee|paneer|labneh|skyr)\b", re.I), "Dairy and Egg Products"),
    (re.compile(r"\b(egg|eggs|omelet|omelette|frittata|quiche)\b", re.I), "Dairy and Egg Products"),
    # Legumes - extended with hummus/falafel/horse gram/sesbania
    (re.compile(r"\b(beans?|lentils?|peas?|chickpeas?|cowpeas?|faba|fava|mung|soy|soybean|edamame|tempeh|tofu|miso|natto|legumes?|lupin|adzuki|cannellini|bambara|groundnuts?|vetch|moth\s+bean|winged\s+bean|pigeon\s+pea|velvet\s+bean|locust\s+bean|hummus|falafel|horse\s+gram|sesbania|black-?eyed\s+pea|garbanzo)\b", re.I), "Legumes and Legume Products"),
    # Nuts & seeds - extended (plurals handled via `s?`)
    (re.compile(r"\b(almonds?|walnuts?|pistachios?|cashews?|pecans?|hazelnuts?|brazil\s+nuts?|macadamias?|peanuts?|chestnuts?|nut|nuts|seed|seeds|chia|flax|flaxseed|sunflower|pumpkin\s+seed|sesame|tahini|poppy|pine\s+nuts?|kola\s+nuts?|cola\s+nuts?)\b", re.I), "Nut and Seed Products"),
    # Beef
    (re.compile(r"\b(beef|steak|brisket|ribeye|sirloin|tenderloin)\b", re.I), "Beef Products"),
    # Pork
    (re.compile(r"\b(pork|ham|bacon|prosciutto|chorizo|salami|sausage)\b", re.I), "Pork Products"),
    # Poultry
    (re.compile(r"\b(chicken|turkey|duck|goose|poultry|hen|quail|partridge|pheasant|guinea\s+fowl|squab|capon)\b", re.I), "Poultry Products"),
    # Lamb / veal / game - extended (horse, hare)
    (re.compile(r"\b(lamb|mutton|veal|deer|venison|bison|elk|rabbit|hare|game|reindeer|caribou|moose|ostrich|kangaroo|antelope|wild\s+boar|emu|alpaca|llama|horse|horse\s+meat)\b", re.I), "Lamb, Veal, and Game Products"),
    # Finfish & shellfish - extended with more species (pike/turbot/flounder/pangasius/whiting/etc.)
    (re.compile(r"\b(fish|salmon|tuna|cod|haddock|tilapia|trout|mackerel|sardine|anchovy|herring|halibut|sole|snapper|sea bass|shrimp|prawn|crab|lobster|crayfish|oyster|mussel|clam|scallop|squid|octopus|seafood|mollusk|mollusks|abalone|albacore|pollock|pollack|sturgeon|hake|saithe|sprat|plaice|eel|swordfish|marlin|monkfish|mullet|kingfish|pompano|escolar|sea\s+urchin|conch|whelk|surimi|caviar|roe|sashimi|fish\s+sauce|amago|hoki|carp|catfish|bass|perch|barramundi|pomfret|smelt|grouper|barracuda|tilefish|john\s+dory|skate|ray|shark|whitefish|crawfish|cuttlefish|fluke|wahoo|mahi[-\s]?mahi|pike|turbot|flounder|pangasius|whiting|hake|burbot|loche|alaska\s+pollock|alaska|red\s+snapper|seabass|tilapia|sturgeon|sushi|sashimi|fish\s+stick|fish\s+fingers|fish\s+ball|fish\s+cake|sea\s+cucumber|krill|prawn|langoustine|crayfish|cockle|periwinkle|whelk|moules|nigiri|maki)\b", re.I), "Finfish and Shellfish Products"),
    # Fruits - extended
    (re.compile(r"\b(apple|apples|pear|pears|peach|peaches|nectarine|plum|plums|cherry|cherries|apricot|apricots|mango|papaya|banana|bananas|pineapple|grape|grapes|berry|berries|strawberry|strawberries|blueberry|blueberries|raspberry|raspberries|blackberry|blackberries|cranberry|cranberries|currant|currants|gooseberry|orange|oranges|lemon|lemons|lime|limes|grapefruit|tangerine|mandarin|kiwi|kiwifruit|melon|melons|watermelon|cantaloupe|honeydew|fig|figs|date|dates|persimmon|pomegranate|guava|jackfruit|durian|lychee|longan|rambutan|soursop|custard apple|baobab|jujube|monkey orange|fruit|fruits|akebia|acerola|prune|prunes|raisin|raisins|olive|olives|wild\s+mango|mulberry|elderberry|gooseberry|bilberry|lingonberry|sea\s+buckthorn|bigney\s+berry|rabiteye)\b", re.I), "Fruits and Fruit Juices"),
    (re.compile(r"\b(juice|nectar|lemonade)\b", re.I), "Fruits and Fruit Juices"),
    # Vegetables - extended with celeriac, swede, beets variations, collards, sauerkraut
    (re.compile(r"\b(carrots?|broccoli|cauliflower|cabbage|brussels\s+sprouts?|kale|spinach|chard|lettuce|arugula|endive|escarole|watercress|celery|cucumbers?|tomatoes?|onions?|shallot|leeks?|garlic|chives?|peppers?|bell\s+pepper|potatoes?|tuber|sweet\s+potato|yams?|cassava|taro|beet|beetroot|beets?|radish|turnip|parsnips?|rutabaga|swede|kohlrabi|jicama|squash|pumpkin|zucchini|eggplant|aubergine|okra|asparagus|artichoke|fennel|celeriac|mushrooms?|seaweed|kelp|nori|wakame|kombu|algae|agar|algaes?|alga|laver|dulse|hijiki|vegetables?|salad|coleslaw|sauerkraut|plantains?|breadfruit|salsify|chayote|bamboo\s+shoots?|dasheen|cocoyam|ulluco|oca|alfalfa\s+sprout|sprouts?|cress|chicory|salsola|samphire|fiddlehead|bok\s+choy|pak\s+choi|napa|tomatillo|nopal|nopales|edamame|abiyuch|amaranth\s+leaves|wattle\s+seed|sea\s+belt|seabelt|spirulina|chlorella|alligator\s+weed|agathi|leaves|moringa|baobab\s+leaves|cowpea\s+leaves|collard|collards|collard\s+greens|dandelion|dandelion\s+greens|spring\s+greens|spring\s+onion|kohlrabi|water\s+chestnut|water\s+spinach|lotus\s+root|kalanchoe|gourd|loofa|loofah|bitter\s+melon|bitter\s+gourd|chinese\s+cabbage|chinese\s+spinach|japanese\s+horse\s+chestnut|japanese\s+spikenard)\b", re.I), "Vegetables and Vegetable Products"),
    # Spices and Herbs - extended with salt + ajowan + relishes
    (re.compile(r"\b(salt|table\s+salt|sodium\s+chloride|sea\s+salt|kosher\s+salt|rock\s+salt|sumac|allspice|nutmeg|mace|star\s+anise|fenugreek|asafoetida|epazote|sansho|wasabi|ajowan|carum\s+copticum|cumin|caraway|cloves?|oregano|basil|thyme|sage|rosemary|parsley|cilantro|mint|dill|fennel\s+seed|ginger|turmeric|cinnamon|saffron|paprika|chili|cardamom|coriander|mustard\s+seed|bay\s+leaf|tarragon|spices?|herbs?|seasoning|relish|ajvar|chutney|condiment|sambal|harissa)\b", re.I), "Spices and Herbs"),
    # Sweets - extended with candybar / icing / gingerbread / ice cream variants
    (re.compile(r"\b(candy|candies|candybar|chocolate|sweet|sweets|syrup|honey|sugar|sugars|marshmallow|caramel|toffee|fudge|gum|gummy|gummies|jelly|jam|preserves|dessert|pudding|cake|cookie|cookies|pie|pies|brownie|brownies|donut|doughnut|pastry|pastries|tart|truffle|gelato|sherbet|sorbet|ice\s+cream|ice\s+lolly|ice\s+pop|popsicle|icing|frosting|gingerbread|agave|stevia|aspartame|saccharin|sucralose|xylitol|erythritol|maple\s+syrup|corn\s+syrup|maltose\s+syrup|sweetener|sweeteners|liquorice|licorice|halva|nougat|baklava|strudel|funnel\s+cake|tres\s+leches|tiramisu|cheesecake|key\s+lime\s+pie|cobbler|crumble|cr[èe]me\s+br[uû]l[ée]e|panna\s+cotta|mousse|flan|cust[a]rd|jellies|frozen\s+yogurt|swiss\s+roll|swiss\s+roll\s+dough|bounty|kitkat|mars\s+bar|snickers|twix|m&m|reese)\b", re.I), "Sweets"),
    # Beverages (non-alcoholic) - extended
    (re.compile(r"\b(coffee|tea|drink|beverage|soda|cola|cocoa|kombucha|smoothie|water|sparkling|chai|matcha|yerba\s+mate|horchata|chocolate\s+milk\s+drink|energy\s+drink|sports\s+drink|protein\s+shake|powdered\s+drink|infusion)\b", re.I), "Beverages"),
    # Fats & Oils - extended
    (re.compile(r"\b(oil|fat|lard|tallow|shortening|margarine|mayonnaise|vinegar|dressing|sauce|gravy|salsa|ketchup|mustard|aioli|olive\s+tapenade|ghee|schmaltz|suet|chicken\s+fat|duck\s+fat|drippings)\b", re.I), "Fats and Oils"),

    # Phase 8 amendment 5: additions surfaced by uncategorized.txt scan
    # Fish/Shellfish - additional species not in the original regex.
    (re.compile(r"\b(anglerfish|anchovy\s+paste|anchovis|ascidian|arctic\s+char|char|bayad|bigeye\s+scad|scad|black\s+seabream|seabream|sea\s+bream|grenadier|grayling|golden\s+redfish|redfish|smelt|dab|gar|kutum|loche|burbot|wolffish|john\s+dory|pilchard|menhaden|sablefish|saithe|orange\s+roughy|sea\s+bream|spot\s+prawn|sea\s+urchin|sea\s+cucumber|sea\s+grape)\b", re.I), "Finfish and Shellfish Products"),
    # Insects (food source - palm weevil larva, ants, grasshoppers, etc.)
    (re.compile(r"\b(ant|ants|bamboo\s+caterpillar|caterpillar|caterpillars|palm\s+weevil|weevil|grasshopper|grasshoppers|cricket|crickets|cricket\s+nymph|nymph|larva|larvae|silkworm|locust|cicada|mealworm|mealworms|edible\s+insect|edible\s+insects|insect|insects)\b", re.I), "Sausages and Luncheon Meats"),
    # Game animals - standalone single-token names that the broader regex misses
    (re.compile(r"\b(armadillo|bear|beaver|kangaroo|wallaby|possum|nutria|capybara|guinea\s+pig|iguana|snake|frog|frog\s+legs|turtle|crocodile|alligator|monkey|deer\s+meat|wild\s+game)\b", re.I), "Lamb, Veal, and Game Products"),
    # Cheeses - specific named cheeses that don't carry the literal "cheese" word
    (re.compile(r"\b(appenzeller|gorgonzola|gruy[eè]re|greyerzer|halloumi|mascarpone|mozzarella|feta|ricotta|manchego|provolone|emmental|emmenthal|camembert|brie|stilton|roquefort|gouda|edam|colby|monterey\s+jack|asiago|gjetost|munster|limburger|fontina|raclette|tilsit|tilsiter|burrata|paneer|labneh|skyr|ymer|ziger|quark|cottage|tomme|reblochon|chevre|chèvre|caciotta|comté|comte|beaufort|cheddar)\b", re.I), "Dairy and Egg Products"),
    # Yeast / fermented bases / nutritional yeast
    (re.compile(r"\b(yeast|baker's\s+yeast|baker.s\s+yeast|brewer's\s+yeast|brewer.s\s+yeast|nutritional\s+yeast|active\s+dry\s+yeast|marmite|vegemite|yeast\s+extract|yeast\s+flakes)\b", re.I), "Spices and Herbs"),
    # Regional grains / legumes / pseudo-cereals not in original regex
    (re.compile(r"\b(bajra|jowar|ragi|makhana|bao\s+bun|basbousa|atole|arepa|bibimbap|bannock|bagel|babaganoush|baba\s+ghanoush|injera|tamale|biryani|paella|chow\s+mein|fried\s+rice|congee|porridge|dahl|dal|daal|dosa|idli|upma|paratha|chapati|roti|naan|pita|focaccia|fougasse|cornpone|hush\s+puppy|hushpuppy|kasha|kheer|payasam|kebab|kibbeh|moussaka)\b", re.I), "Cereal Grains and Pasta"),
    # Indian / African regional legumes - Bengal gram, black gram, green gram, horse gram
    (re.compile(r"\b(bengal\s+gram|black\s+gram|green\s+gram|red\s+gram|horse\s+gram|gram\s+flour|besan|toor|tur|urad|moong|chana|rajma|matki|val|kulthi|benniseed)\b", re.I), "Legumes and Legume Products"),
    # Tropical / niche fruits not in original regex
    (re.compile(r"\b(ackee|atemoya|babaco|bilberry|bilberries|black\s+crowberry|black\s+nightshade|gojiberry|goji\s+berry|wolfberry|greengage|greengages|cloudberry|cloudberries|sea\s+buckthorn|tamarind|tamarillo|cherimoya|sapodilla|sapote|mamey|jicama|loquat|kaffir\s+lime|kumquat|feijoa|carambola|starfruit|santol|noni|breadnut|jackfruit|monstera|miracle\s+fruit|durian|mangosteen|salak|snake\s+fruit|surinam\s+cherry|gandaria|safou|safu|safou\s+fruit|akebia|wood-?sorrel|arecanut)\b", re.I), "Fruits and Fruit Juices"),
    # Niche vegetables / roots / mushrooms
    (re.compile(r"\b(amla|arrowhead|antroewa|bele|yautia|woolly\s+milkcap|morel|porcini|button\s+mushroom|enoki|shimeji|oyster\s+mushroom|maitake|king\s+oyster|chanterelle|trumpet|truffle\s+mushroom|chinese\s+yam|yacon|jicama|callaloo|moringa\s+leaves|nopales|chayote|christophine|fiddlehead\s+ferns?|salsify|water\s+chestnut|achiote|annatto|mushroom)\b", re.I), "Vegetables and Vegetable Products"),
    # Composite dishes / meals
    (re.compile(r"\b(goulash|stew|stews|casserole|lasagne|lasagna|moussaka|paella|risotto|jambalaya|gumbo|sukiyaki|sushi|nigiri|maki|onigiri|bento|kimchi|sauerkraut|coleslaw|salad\s+bowl|grain\s+bowl|buddha\s+bowl|rice\s+bowl|noodle\s+bowl|stir\s+fry|stir-fry|stir-fried|wok|wok\s+mix|tagine|biryani|gallo\s+pinto|huevos\s+rancheros|enchilada|fajita|burrito\s+bowl|frittata\s+bowl|fried\s+rice|bibimbap|bami\s+goreng|bami|nasi\s+goreng|nasi|tabbouleh|hummus\s+bowl|bowl)\b", re.I), "Meals, Entrees, and Side Dishes"),
    # Bitters / regional spirits
    (re.compile(r"\b(bitter|bitters|jägermeister|jagermeister|gammel\s+dansk|aquavit|akvavit|krabask|absinthe|amaretto|grappa|schnapps|raki|arak|ouzo|sambuca|chartreuse|fernet|kvass|cachaça|cachaca|mezcal|baiju|baijiu|shochu|soju|umeshu|amaro|aperitif|digestif|pernod|jagermeister|aperol)\b", re.I), "Alcoholic Beverages"),
    # Pickles / fermented vegetables / condiments missed by base Vegetables regex
    (re.compile(r"\b(pickle|pickles|kimchi|fermented\s+vegetable|miso|fish\s+sauce|fish\s+paste|shrimp\s+paste|patis|garum|colatura|achar|chow\s+chow|relish|chutney|salsa\s+verde)\b", re.I), "Spices and Herbs"),
    # Honey, royal jelly, propolis, bee products
    (re.compile(r"\b(honey|honeydew\s+honey|royal\s+jelly|propolis|bee\s+pollen|beeswax|manuka|honeycomb)\b", re.I), "Sweets"),

    # Phase 8 amendment 5 - round 2: gaps found after running the rebuild
    # Fruits - common terms missed by the base regex due to word boundaries (apple
    # doesn't match applesauce; avocado not in original list at all; plurals)
    (re.compile(r"\b(avocado|avocados|applesauce|apple\s+sauce|apricots?|gooseberries|blackcurrants?|mulberries|cloudberries|loganberries|boysenberries|raspberries|blueberries|cranberries|cherries|strawberries|grapes|pears|peaches|plums|figs|dates|olives|peach\s+halves|fruit\s+salad|fruit\s+cup|fruit\s+cocktail|fruit\s+jelly|fruit\s+leather)\b", re.I), "Fruits and Fruit Juices"),
    # Cereal Grains - bakery + baking aids + bagel plurals + missed cereal types
    (re.compile(r"\b(bagels?|bakery\s+mix|bakery|biscuit\s+dough|cookie\s+dough|cake\s+mix|baking\s+powder|baking\s+soda|baking\s+yeast|raising\s+agent|leavening\s+agent|binding\s+agent|bannocks?|babka|brioche|bun|buns|fritter|fritters|donut|donuts|doughnut|doughnuts|cinnamon\s+roll|cinnamon\s+rolls|empanada|empanadas|knish|pretzel|pretzels|popadom|papad|millet|grits|farina|semolina|injera|chapati|paratha|naan|tortillas?|matzo|matzah|hardtack)\b", re.I), "Cereal Grains and Pasta"),
    # Fish - additional patterns missed (codfish, bigeye, biscayne, plurals)
    (re.compile(r"\b(codfish|bacalao|bacalhau|bacalaitos|biscayne\s+cod|fish\s+sticks|fish\s+balls|fish\s+cakes|crab\s+sticks|imitation\s+crab|herring\s+roe|salmon\s+roe|caviars?)\b", re.I), "Finfish and Shellfish Products"),
    # Vegetables - babycorn / niche missed
    (re.compile(r"\b(babycorn|baby\s+corn|miniature\s+corn|sweetcorn|sweet\s+corn|corn\s+on\s+the\s+cob|corn\s+kernels|antroewa|asiatic\s+dayflower|ambrosia\s+greens|wood-?sorrel|sourgrass|buttercup|cape\s+sorrel|sour\s?sob|bermuda\s+buttercup|sorrel|chinese\s+broccoli|gai\s+lan|brussel\s+sprouts?|broccolini|broccoli\s+rabe|rapini|frisée|frisee|raddichio|radicchio|romaine|iceberg|butter\s+lettuce|tatsoi|mizuna)\b", re.I), "Vegetables and Vegetable Products"),
    # Fats and oils - blended spreads, low-fat spreads, margarine variants
    (re.compile(r"\b(blended\s+spread|spread,\s+blended|low-?fat\s+spread|reduced-?fat\s+spread|table\s+spread|fat\s+spread|spread\s+\d+%|margarine\s+spread)\b", re.I), "Fats and Oils"),
    # Sweets - desserts missed
    (re.compile(r"\b(beignet|beignets|blancmange|baba|babas|baba\s+au\s+rhum|babka|cannoli|cannolis|sufganiyot|profiterole|profiteroles|chocolate\s+egg|chocolate\s+coin|chocolate\s+bar|chocolate\s+bark|chocolate\s+chip|chocolate\s+chunk|chocolate\s+covered|chocolate\s+filled|fruit\s+leather|fruit\s+roll|fruit\s+snack|granola\s+cluster|cereal\s+cluster|popsicle|ice\s+pop|sundae|sorbet|sherbet|halva|halawa|jalebi|gulab\s+jamun|lokum|turkish\s+delight)\b", re.I), "Sweets"),
    # Fast foods - big mac / specific menu items
    (re.compile(r"\b(big\s+mac|whopper|cheeseburger|chicken\s+sandwich|fish\s+sandwich|happy\s+meal|kids\s+meal|fries|onion\s+rings|chicken\s+strips|chicken\s+tenders|chicken\s+wings|wings|buffalo\s+wings|filet[-\s]o[-\s]fish|mcrib|mcmuffin|mc\w+|kfc|burger\s+king)\b", re.I), "Fast Foods"),
    # Soups - common Western dishes missed
    (re.compile(r"\b(borscht|miso\s+soup|tom\s+yum|tom\s+kha|pho|menudo|albondigas|cock-?a-?leekie|scotch\s+broth|cullen\s+skink|mulligatawny|posole|pozole)\b", re.I), "Soups, Sauces, and Gravies"),
    # Spices/other - yeasts, baking aids, leavening, gelatins, agar agar
    (re.compile(r"\b(gelatin|gelatine|agar\s+agar|gum\s+arabic|xanthan\s+gum|guar\s+gum|locust\s+bean\s+gum|carrageenan|pectin|lecithin|citric\s+acid|tartaric\s+acid|malic\s+acid|food\s+colour|food\s+coloring|food\s+color|emulsifier|stabilizer|preservative|food\s+additive)\b", re.I), "Spices and Herbs"),
    # Dairy - egg products missed + ymer + skyr variants
    (re.compile(r"\b(egg\s+white|egg\s+whites|egg\s+yolk|egg\s+yolks|egg\s+powder|powdered\s+eggs?|liquid\s+eggs?|egg\s+substitute|condensed\s+milk|evaporated\s+milk|powdered\s+milk|milk\s+powder|cream\s+cheese|cottage\s+cheese|sour\s+cream|whipped\s+cream|half\s+and\s+half|half-and-half|heavy\s+cream|light\s+cream|skim\s+milk|whole\s+milk|two\s+percent\s+milk|2%\s+milk)\b", re.I), "Dairy and Egg Products"),
    # Final catch-all - single-word foods sometimes have a 'foods' / 'food' / 'meal' descriptor
    (re.compile(r"\b(prepared\s+food|prepared\s+meal|ready\s+meal|ready\s+made\s+meal|frozen\s+meal|microwaveable\s+meal|tv\s+dinner|frozen\s+entree|frozen\s+entr[ée]e|meal\s+kit|home\s+meal)\b", re.I), "Meals, Entrees, and Side Dishes"),
]

# Key nutrients for the similarity gate (parameters.yaml: nutrient_ids.key_nutrient_ids):
# macros, sugars, fiber, water, substrate-relevant fibers. The trace-polyphenol tail is
# left out, measurement noise dominates it. Energy (1008) too: FDC reports kcal and
# BioFoodComp/Ciqual/Frida sometimes kJ under the same nutrient_id, and the 4.184x
# artifact split rice/pork/beef variants as false outliers. It is derivable from
# protein/fat/carb anyway.
# _OUTLIER_RATIO (parameters.yaml: scoring.outlier_ratio) is the deviation from the group
# median, on any one key nutrient, that gives a variant its own canon. 3.0 -> 5.0 in
# Phase 9: 3x split rice cultivars (njavara, IR 64, PTB 39) on small fiber differences,
# 5x still splits rice bran (10x+) and brown vs white rice (>5x).

def canonicalize_food_name(desc: str) -> str:
    if not desc: return ""
    d = str(desc).strip()
    # Strip FDC Foundation Foods nutrient-panel group prefix FIRST
    # ("Proximates, ", "Beverages, ", "Cereals ready-to-eat, ", "Cereals, ")
    # - must run before _NF_SUFFIX_RE because the NF suffix regex is greedy
    # and can eat the entire informative middle of long descriptions.
    d = _FDC_GROUP_PREFIX_RE.sub("", d)
    # Strip FDC NF panel-fragmentation suffix (two flavors).
    d = _NF_SUFFIX_RE.sub("", d)
    d = _NF_BARE_RE.sub("", d)
    # Strip trailing source tag e.g. " [BioFoodComp]" - keep variant info inside
    # the brackets out of the canon name.
    d = _BRACKET_TAG_RE.sub("", d)
    # Then strip any remaining bracketed variant tags anywhere in the description:
    # color / cultivar / sci-name tags that Phenol-Explorer + BioFoodComp use
    # ("Tea [Oolong]", "Orange [Blond]", "Common cabbage [Purple]", etc.).
    d = _BRACKET_VARIANT_RE.sub(" ", d)
    # Strip preparation / state tokens.
    d = _PREP_RE.sub("", d)
    # Strip parenthetical clarifications and quantitative qualifiers - these
    # are non-essential annotations that prevent head-aware grouping when
    # they appear in the second comma-chunk of preserve-head foods. For
    # strip-list heads the second chunk is dropped anyway, so this is a
    # no-op there; the effect is on preserve heads ("potato, white
    # (industrial), 50% extraction" → "potato, white").
    d = _PAREN_CONTENT_RE.sub(" ", d)
    d = _QUANT_RE.sub("", d)
    # Cleanup: drop empty parens left by stripped tokens, collapse repeated
    # commas/spaces, trim punctuation.
    d = _EMPTY_PARENS_RE.sub(" ", d)
    d = re.sub(r"\s+", " ", d)
    d = re.sub(r"\s*,\s*", ", ", d)
    d = re.sub(r"(?:,\s*)+,", ",", d)
    d = re.sub(r"^[\s,-]+|[\s,-]+$", "", d)
    d = d.lower()
    # Head-aware trimming. After all the noise is gone, split on commas and
    # drop cultivar/variety information for known cultivar-flooded heads
    # (apple/strawberry/blueberry/etc.). For everything else, keep at most
    # the first two chunks so genuinely meaningful varieties survive
    # (beans, navy / mushrooms, shiitake / wheat, khorasan) but free-text
    # qualifiers don't accumulate.
    chunks = [c.strip() for c in d.split(",") if c.strip()]
    if not chunks:
        return d
    # Slashed-synonym head: some BioFoodComp / WAFCT entries pack multiple
    # common names of the same food (or several related foods) into the
    # description's first chunk, separated by "/ ". Examples:
    #   "bermuda buttercup/ african wood-sorrel/ ... / soursop, bulb"
    #   "creeping bauhinia/ marama bean/ tamani berry, without testa"
    #   "obscure morning glory/ small white morning glory, leaf"
    # The discriminating signal is **slash followed by space** - the
    # BioFoodComp / WAFCT convention for listing synonyms. Single slashes
    # without a following space are alternation ("Beef, rib eye steak/roast",
    # "Beef, top loin/sirloin/round") which we preserve as-is. Fallback rule
    # ≥3 raw slashes still triggers for densely-slashed entries that may
    # have inconsistent spacing.
    if "/ " in chunks[0] or chunks[0].count("/") >= 3:
        chunks[0] = chunks[0].split("/", 1)[0].strip()
    if chunks[0] in _CULTIVAR_STRIP_HEADS:
        return chunks[0]
    # For preserve heads, walk chunks[1:] and drop any that look like
    # cultivar codes ([Oryza sativa], 'CDC Blaze', var. SLS1, ADT-21,
    # angus x holstein-friesian, etc.); keep the first chunk that doesn't
    # match. This collapses "rice, [oryza sativa]" → "rice" but leaves
    # "rice, brown" → "rice, brown" untouched.
    kept = []
    for c in chunks[1:]:
        if _CULTIVAR_CODE_RE.match(c):
            continue
        kept.append(c)
        break
    return ", ".join([chunks[0]] + kept[:1])


load_smart    = lambda p: pd.read_parquet(str(p)) if str(p).lower().endswith((".parquet",".pq")) else pd.read_csv(str(p), low_memory=False)
classify_rule = lambda r: (False,"macro_proxy") if (r:=(r or "").strip().lower()) in DANGEROUS_RULES else ((True,"form") if any(r.startswith(p) for p in SAFE_RULE_PFXS) else (False,"unknown_proxy"))
get_cov       = lambda x, th: 0.0 if th <= 0 else min(max(x/th, 0.0), 1.0)


def cosine_sim(a, b):
    if not a or not b: return 0.0
    i = set(a) & set(b)
    if not i: return 0.0
    dot = sum(a[k]*b[k] for k in i)
    na = math.sqrt(sum(v*v for v in a.values())); nb = math.sqrt(sum(v*v for v in b.values()))
    return dot/(na*nb) if na > 0 and nb > 0 else 0.0


def round_floats(df, nd=4):
    if df is None or df.empty: return df
    df = df.copy()
    c = df.select_dtypes(include=["float16","float32","float64"]).columns
    if len(c): df[c] = df[c].round(nd)
    return df


class TopKByNutrient:
    """Heap-based per-nutrient top-K food tracker. Source: the query prototype."""
    __slots__ = ("val", "heap")
    def __init__(self):
        self.val = {}; self.heap = {}
    def push(self, nid, canon, amt, K):
        d, h = self.val.setdefault(nid, {}), self.heap.setdefault(nid, [])
        if d.get(canon, -1) >= amt: return
        d[canon] = amt; heapq.heappush(h, (amt, canon))
        if len(d) > K*2:
            newh = [(a, c) for a, c in h if d.get(c) == a]; heapq.heapify(newh); self.heap[nid] = newh
        while len(d) > K and self.heap[nid]:
            a, c = heapq.heappop(self.heap[nid])
            if d.get(c) == a: del d[c]


def refine_canons_by_nutrition(canon_groups: dict, fid_to_desc: dict,
                                food_nutrient_path: str) -> dict:
    """Split multi-variant canon groups when their nutrient profiles disagree.

    Cultivars whose key-nutrient amounts deviate from the group median by
    more than `_OUTLIER_RATIO` on at least one nutrient are pulled out of
    the head-stripped canon and given their own canon - named with the
    variant's original FDC description (lower-cased, light cleanup) so the
    user sees `strawberry, oso grande` instead of `strawberry`.

    The scan is filtered tightly: only the fdc_ids that are in multi-member
    groups and only the ~11 nutrients in `_KEY_NUTRIENT_IDS`. ~10–30 s on
    the full FDC catalog.
    """
    multi = {c: fids for c, fids in canon_groups.items() if len(fids) > 1}
    if not multi:
        return canon_groups
    all_fids = sorted({f for fids in multi.values() for f in fids})
    n_total_variants = sum(len(v) for v in multi.values())
    print(f"[refine] {len(multi)} multi-variant canon groups ({n_total_variants} fdc_ids); "
          f"scanning {len(_KEY_NUTRIENT_IDS)} key nutrients...", flush=True)

    tab = (ds.dataset(food_nutrient_path, format="parquet", partitioning="hive")
             .scanner(
                 columns=["fdc_id", "nutrient_id", "amount"],
                 filter=ds.field("fdc_id").isin(all_fids)
                        & ds.field("nutrient_id").isin(list(_KEY_NUTRIENT_IDS)),
             ).to_table())
    if tab.num_rows == 0:
        print("[refine] no key-nutrient rows found; skipping refinement.", flush=True)
        return canon_groups
    nut_df = tab.to_pandas()
    # Pivot: fid → {nid: amount}. Wide is faster than dict-of-dict here.
    pivot = nut_df.pivot_table(index="fdc_id", columns="nutrient_id",
                                values="amount", aggfunc="mean")

    refined: dict[str, list[int]] = {}
    n_split = 0
    n_exempt = 0
    for cname, fids in canon_groups.items():
        if cname not in multi:
            refined.setdefault(cname, []).extend(fids)
            continue
        # Strip-list head exemption: HARD heads are gate-exempt (user intent
        # is "collapse every variant regardless of nutritional differences"
        # - e.g. ripe vs unripe mango, 9% vs whole-wheat flour, apple
        # cultivars). GATED heads (rice / pork / beef / veal) still go through
        # the outlier check so genuinely distinct varieties (brown vs white
        # rice; Iberian vs commercial pork) split off as their own canon.
        if cname in _CULTIVAR_STRIP_HEADS_HARD:
            refined.setdefault(cname, []).extend(fids)
            n_exempt += 1
            continue
        # Subset pivot to this group's fids
        sub = pivot.reindex(fids).dropna(how="all")
        if sub.empty or sub.shape[0] < 2:
            refined.setdefault(cname, []).extend(fids)
            continue
        median = sub.median(axis=0, skipna=True)  # series indexed by nutrient_id
        # max deviation ratio per fid across measured nutrients
        denom = median.where(median.abs() > 1e-6).abs()
        ratios = (sub.subtract(median, axis=1).abs()).divide(denom, axis=1)
        max_ratio = ratios.max(axis=1, skipna=True)
        outliers = set(max_ratio[max_ratio > _OUTLIER_RATIO].index.tolist())
        keepers = [f for f in fids if f not in outliers]
        if keepers:
            refined.setdefault(cname, []).extend(keepers)
        for f in outliers:
            desc = fid_to_desc.get(f, "").strip()
            out_canon = canonicalize_food_name(desc) or desc.lower() or f"fdc_{f}"
            # If the outlier canonicalizes back to the same name (rare), keep
            # the longer original description so it's visibly distinct.
            if out_canon == cname:
                out_canon = desc.lower().strip(" ,-")
                if not out_canon or out_canon == cname:
                    out_canon = f"{cname} (variant {f})"
            refined.setdefault(out_canon, []).append(f)
            n_split += 1

    print(f"[refine] {n_split} cultivar variants split off as their own canons "
          f"(profile deviation > {_OUTLIER_RATIO}× group median); "
          f"{n_exempt} strip-list-head groups exempted from gate.", flush=True)
    return refined


def build_static_food_meta(args, out_path):
    """Build a dense numpy index over the food catalog. Adapted from the query prototype."""
    print("[*] Building static food meta...", flush=True)
    if args.food_portion and Path(args.food_portion).exists():
        port_df = load_smart(args.food_portion)
        if "gram_weight" not in port_df.columns and "amount" in port_df.columns:
            port_df = port_df.rename(columns={"amount": "gram_weight"})
        port_df = port_df.dropna(subset=["gram_weight"])
        port_df["is_serving"] = port_df.get("modifier", pd.Series([""]*len(port_df))).astype(str).str.lower().str.contains("serving")
        portion_map = port_df.sort_values(["fdc_id","is_serving"], ascending=[True,False]).groupby("fdc_id")["gram_weight"].first().to_dict()
    else:
        print(f"[!] food_portion file missing ({args.food_portion}); using default 50 g per food.", flush=True)
        portion_map = {}
    nut = load_smart(args.nutrient).rename(columns={"id":"nutrient_id","name":"nutrient_name"})
    unit_map = nut.set_index("nutrient_id")["unit_name"].astype("string").str.upper().to_dict()
    modeled = set(nut["nutrient_id"].dropna().astype(int).tolist())
    # Also include every nutrient_id physically present in the bucketed
    # food_nutrient parquet. This covers synthetic substrate IDs (e.g.,
    # 200001+ for HMOs / sialic acid / fucose / alginate / agarose) injected
    # via 0_building/inject_bacterial_substrates.py - without this they'd be
    # silently dropped at score_one_bacterium because the targs filter is
    # `all_t & STATIC_DB["modeled"]`.
    try:
        extra_nids = set(int(x) for x in ds.dataset(args.food_nutrient, format="parquet", partitioning="hive")
                         .to_table(columns=["nutrient_id"])
                         .column("nutrient_id").to_pylist())
        extra_nids -= modeled
        if extra_nids:
            print(f"[*] Modeled: adding {len(extra_nids)} non-FDC nutrient_ids found in food_nutrient parquet "
                  f"(e.g. synthetic bacterial substrates).", flush=True)
            modeled |= extra_nids
    except Exception as e:
        print(f"[!] Could not scan food_nutrient for extra ids: {e}", flush=True)
    id2desc = {}
    for _, row in load_smart(args.food_category).rename(columns={"id":"food_category_id","description":"food_category"}).iterrows():
        raw = str(row["food_category_id"]).strip(); name = str(row["food_category"]); id2desc[raw] = name
        try: id2desc[str(int(float(raw)))] = name
        except: pass
    food_df = load_smart(args.food).drop_duplicates(subset=["fdc_id"], keep="last")
    # Normalize embedded whitespace in descriptions. Some Ciqual / BioFoodComp
    # rows have literal '\n' or '\t' inside the description column (e.g.
    # "Tamarind, mature fruit, flesh without skin,\n with seeds, raw [CIQUAL]"),
    # which breaks any consumer that parses the differential.tsv line-by-line.
    food_df["description"] = food_df["description"].astype(str) \
        .str.replace(r"[\r\n\t]+", " ", regex=True) \
        .str.replace(r"\s{2,}", " ", regex=True) \
        .str.strip()
    meta = food_df.set_index(pd.to_numeric(food_df["fdc_id"], errors="coerce"))[["description","data_type","food_category_id"]].to_dict(orient="index")
    max_fdc = int(max(k for k in meta if pd.notna(k))) + 1000
    c_arr = np.full(max_fdc, -1, dtype=np.int32); p_arr = np.full(max_fdc, 50.0, dtype=np.float32)
    pl_arr = np.zeros(max_fdc, dtype=bool); sp_arr = np.zeros(max_fdc, dtype=bool)
    drop_cats = set(str(x).strip() for x in args.drop_category); food_stats = {}
    _paper_title_skipped = 0
    for fdc_id, m in meta.items():
        if pd.isna(fdc_id): continue
        fid = int(fdc_id); dt = str(m.get("data_type","")); dt_norm = dt.lower().strip().replace(" ","_")
        cat_raw = str(m.get("food_category_id","")).strip(); cat = id2desc.get(cat_raw, "")
        if not cat:
            try: cat = id2desc.get(str(int(float(cat_raw))), "")
            except (ValueError, OverflowError): pass
        desc_lc = str(m.get("description","")).lower().strip()
        # Filter ETL artifacts: 60 fdc_ids in food.parquet have a literal "None"
        # / "nan" / "null" string as their description (stringified Python None
        # from some upstream import). They have no useful identity and were
        # producing the bizarre `food_name="none"` rows in differential output.
        if desc_lc in ("", "none", "nan", "null", "n/a", "na"):
            continue
        # Filter research-paper-title leakage from BioFoodComp: some rows in
        # the source have a literature reference in the description column
        # instead of a food name ("A comprehensive characterization of
        # phenolics, amino acids and other minor bioactives of selected
        # honeys..."). They masquerade as canonical foods, never categorize,
        # and pollute the uncategorized diagnostic. Length > 60 chars AND a
        # journal-style keyword catches them without dropping legitimate
        # long food names (real food entries with journal-style keywords are
        # rare - `comprehensive`, `characterization`, etc. almost never appear
        # in genuine food descriptions).
        if len(desc_lc) > 60 and _PAPER_TITLE_RE.search(desc_lc):
            _paper_title_skipped += 1
            continue
        if (DROP_BRANDED and dt_norm=="branded_food") or (DROP_MODELLED and dt_norm=="survey_fndds_food") or cat in drop_cats or cat in ALWAYS_DROP_CATS or WHALE_RE.search(desc_lc):
            continue
        c_arr[fid] = fid; p_arr[fid] = portion_map.get(fid, 50.0)
        sp_arr[fid] = (cat == "Spices and Herbs"); pl_arr[fid] = (cat in PLANT_CATS)
        if "[phenol-explorer]" in desc_lc:
            # Phenol-Explorer covers pasta, bread, beverages, etc. - use the
            # full _DESC_TO_CATEGORY fallback; only default to "Fruits and
            # Fruit Juices" if literally nothing matches.
            if not cat:
                for pat, fallback_cat in _DESC_TO_CATEGORY:
                    if pat.search(desc_lc):
                        cat = fallback_cat
                        break
                cat = cat or "Fruits and Fruit Juices"
            p = 0.0; bc = 0.0; pl_arr[fid] = True
        elif "[biofoodcomp]" in desc_lc:
            # Use _DESC_TO_CATEGORY directly. The previous version had an
            # inline 4-keyword pre-check that used `substring in desc_lc`,
            # which falsely matched "oat" inside "Goat's-foot" (Bermuda
            # buttercup row → Cereal Grains) and similar partial-word hits.
            # _DESC_TO_CATEGORY uses word-boundary regex so it doesn't
            # have that bug.
            if not cat:
                for pat, fallback_cat in _DESC_TO_CATEGORY:
                    if pat.search(desc_lc):
                        cat = fallback_cat
                        break
                cat = cat or "Vegetables and Vegetable Products"
            p = float(CAT_PENALTY.get(cat, 0.0)); bc = PROC_BASE_W * p + 0.40; pl_arr[fid] = True
        else:
            # Last-resort keyword fallback for any other source missing a
            # food_category_id (AFCD, STFCJ, MERIDA, CNF, miscellaneous).
            if not cat:
                for pat, fallback_cat in _DESC_TO_CATEGORY:
                    if pat.search(desc_lc):
                        cat = fallback_cat
                        break
            p = float(CAT_PENALTY.get(cat, 1.0)); bc = PROC_BASE_W * p
        if NFY_RE.search(desc_lc): p = max(p, 1.5)
        if ALCOHOL_RE.search(desc_lc): p = max(p, 2.0); bc += ALC_BASE_W
        if JUNK_RE.search(desc_lc) or POWDER_RE.search(desc_lc):
            p = max(p, 2.0); bc += (JUNK_BASE_W if JUNK_RE.search(desc_lc) else 0.40)
        if PROCESSED_RE.search(desc_lc): bc += 0.35
        tp = TP_MAP.get(dt_norm, 100)
        food_stats[fid] = {"dh":int(zlib.crc32(desc_lc.encode())),"pl":p,"bc":bc,
                           "art":int(p>=1.5),"tb":TYPE_W*(tp/600.0),"tp":tp,"dt_norm":dt_norm,"cat":cat}
    # Canonical-name grouping
    # FDC fragments the same food across many rows (per nutrient panel via NF
    # suffixes, plus prep/state qualifiers). Group every fdc_id sharing a
    # canonical short name and elect one representative per group; remap
    # c_arr[all variants] -> rep_fid so build_modeled_index folds their
    # nutrients together (MEAN over the variant set).
    canon_groups: dict[str, list[int]] = {}
    fid_to_desc: dict[int, str] = {}
    for fid in food_stats:
        desc_orig = str(meta.get(fid, {}).get("description", ""))
        fid_to_desc[fid] = desc_orig
        cname = canonicalize_food_name(desc_orig) or desc_orig.lower()
        canon_groups.setdefault(cname, []).append(fid)

    # Nutritional-similarity refinement: cultivars with profiles diverging
    # from the group median by > _OUTLIER_RATIO on any key nutrient fall out
    # of the head-stripped canon and keep their original descriptive name.
    canon_groups = refine_canons_by_nutrition(canon_groups, fid_to_desc, args.food_nutrient)

    rep_per_canon: dict[str, int] = {
        cname: max(fids, key=lambda f: (food_stats[f]["tp"], -f))  # highest tp; tiebreak lowest fid
        for cname, fids in canon_groups.items()
    }
    # Remap c_arr: every variant points to its rep.
    for cname, fids in canon_groups.items():
        rep = rep_per_canon[cname]
        for f in fids:
            c_arr[f] = rep
    # Keep food_stats only for reps; attach the canonical display name.
    rep_stats = {}
    for cname, rep in rep_per_canon.items():
        s = dict(food_stats[rep])
        s["canon_name"] = cname
        s["n_variants"] = len(canon_groups[cname])
        rep_stats[rep] = s
    # Validity mask: keep every variant alive so its parquet rows fold into
    # the rep. Only drop branded/spice/whale/etc. that were never added to
    # food_stats in the first place.
    vm = np.zeros(max_fdc, dtype=bool)
    if food_stats: vm[list(food_stats.keys())] = True
    c_arr[~vm] = -1

    n_groups, n_variants = len(rep_per_canon), len(food_stats)
    print(f"[*] Canonicalization: {n_variants} FDC entries → {n_groups} canonical foods "
          f"(merge ratio {n_variants/max(1,n_groups):.2f}x)", flush=True)

    # Under-grouped diagnostic: surface heads (first comma-chunk) that still
    # resolve to multiple canons after refinement, so the user can spot
    # candidates for additional strip-list entries. Sorted by sub-canon
    # count. Top entries should be PRESERVE-head families with meaningful
    # varieties (cheese / mushrooms / beans / soup / cereal / ...); any
    # strip-list head appearing here is a regression.
    heads_to_canons: dict[str, list[str]] = {}
    for cname in rep_per_canon.keys():
        head = cname.split(",")[0].strip() if "," in cname else cname
        heads_to_canons.setdefault(head, []).append(cname)
    flooded = sorted(
        [(h, list(set(cs))) for h, cs in heads_to_canons.items() if len(set(cs)) > 1],
        key=lambda x: -len(x[1]),
    )[:25]
    if flooded:
        print(f"[diag] Under-grouped heads (top 25 by sub-canon count) - "
              f"add to _CULTIVAR_STRIP_HEADS if cultivar-flood, leave if meaningful varieties:",
              flush=True)
        for h, cs in flooded:
            marker = "STRIP-LIST-LEAK" if h in _CULTIVAR_STRIP_HEADS else "preserve?"
            print(f"  [{marker:15s}] {len(cs):>4} canons under '{h}'  e.g. {sorted(cs)[:3]}",
                  flush=True)

    max_nut = max(modeled) + 100; u_arr = np.ones(max_nut, dtype=np.float32)
    for n in modeled:
        u = unit_map.get(n, "").upper()
        u_arr[n] = 1000.0 if u == "G" else (0.001 if u == "UG" else 1.0)
    pickle.dump({"canon_arr":c_arr,"port_arr":p_arr,"plant_arr":pl_arr,"spice_arr":sp_arr,"unit_arr":u_arr,
                 "food_stats":rep_stats,"modeled":modeled,
                 "nut_name":dict(zip(nut["nutrient_id"], nut["nutrient_name"].astype(str))),
                 "rep2can":{f:c_arr[f] for f in rep_stats}}, open(out_path, "wb"))
    print("[*] Static food meta written.", flush=True)
    if _paper_title_skipped:
        print(f"[*] Research-paper-title filter: skipped {_paper_title_skipped} BioFoodComp rows "
              "whose description column held a literature reference, not a food name.", flush=True)

    # Category fill diagnostic - list canonical names that still have no
    # category after the keyword fallback so the user can iterate.
    uncat = sorted({s.get("canon_name","") for s in rep_stats.values() if not s.get("cat","")})
    fill_rate = 1 - len(uncat) / max(1, len(rep_stats))
    print(f"[*] Category fill: {fill_rate*100:.1f}% ({len(rep_stats)-len(uncat)} of {len(rep_stats)} canonical foods)", flush=True)
    if uncat:
        uncat_path = out_path.parent / "uncategorized.txt"
        with open(uncat_path, "w", encoding="utf-8") as f:
            for c in uncat:
                f.write(c + "\n")
        print(f"[*] {len(uncat)} canonical foods still uncategorized; list at {uncat_path}", flush=True)


def build_modeled_index(bdir, static_db, out_dir, batch_rows=750_000):
    """Adapted from the query prototype's build_modeled_index."""
    out_dir.mkdir(parents=True, exist_ok=True)
    nids = sorted(int(x) for x in static_db["modeled"]); buckets = sorted({n % BUCKETS for n in nids})
    scn = ds.dataset(bdir, format="parquet", partitioning="hive").scanner(
        columns=["fdc_id","nutrient_id","amount"],
        filter=ds.field("bucket").isin(buckets) & ds.field("nutrient_id").isin(nids),
        batch_size=batch_rows)
    tot_c, tot_m, df_foods = {}, {}, {}
    tmp = out_dir / "_tmp_scale_parts"
    if tmp.exists(): shutil.rmtree(tmp)
    tmp.mkdir(parents=True, exist_ok=True)
    c_arr, p_arr, sp_arr, u_arr = static_db["canon_arr"], static_db["port_arr"], static_db["spice_arr"], static_db["unit_arr"]
    max_f, max_u = len(c_arr), len(u_arr)
    for bi, b in enumerate(scn.to_batches(), 1):
        f = b["fdc_id"].to_numpy().astype(np.int64, copy=False)
        n = b["nutrient_id"].to_numpy().astype(np.int32, copy=False)
        a = b["amount"].to_numpy().astype(np.float32, copy=False)
        m = (f>=0)&(f<max_f)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f, n, a = f[m], n[m], a[m]; c = c_arr[f]; v = (c != -1) & (~sp_arr[f])
        if not v.any(): continue
        f, n, a, c = f[v], n[v], a[v], c[v]; a = a*(p_arr[f]/50.0)*u_arr[n]
        df = pd.DataFrame({"nutrient_id":n,"canon":c,"amount_norm":a})
        for canon, row in df.groupby("canon", dropna=True).agg(cnt=("nutrient_id","size"), mass=("amount_norm","sum")).iterrows():
            cn = int(canon); tot_c[cn] = tot_c.get(cn, 0) + int(row["cnt"]); tot_m[cn] = tot_m.get(cn, 0.0) + float(row["mass"])
        for nid, canon in df[["nutrient_id","canon"]].drop_duplicates().itertuples(index=False):
            df_foods.setdefault(int(nid), set()).add(int(canon))
        # Aggregate across variant fdc_ids that fold into the same canon.
        # MEAN over actually-measured variants (the parquet is sparse-by-row,
        # so unmeasured nutrients have no row to drag the average down). For
        # cultivar-level data this gives the typical value across cultivars;
        # for NF-fragmented FDC entries each variant carries a disjoint
        # subset of nutrients, so MEAN equals the single measured value.
        sub2 = df.groupby(["nutrient_id","canon"], dropna=True)["amount_norm"].mean().reset_index()
        if not sub2.empty:
            sub2["bucket"] = (sub2["nutrient_id"].astype(int) % BUCKETS).astype(int)
            sub2.to_parquet(tmp/f"part_{bi:06d}.parquet", index=False)
        del df, f, n, a, c
        if bi % 10 == 0: gc.collect()
    pd.DataFrame({"canon":list(tot_c),"model_count":[int(tot_c[c]) for c in tot_c],
                  "model_mass":[float(tot_m.get(c, 0)) for c in tot_c]}).to_parquet(out_dir/"modeled_totals.parquet", index=False)
    pd.DataFrame({"nutrient_id":list(df_foods),"df_foods":[len(df_foods[n]) for n in df_foods]}).to_parquet(out_dir/"nutrient_df.parquet", index=False)
    rows = []; tmp_ds = ds.dataset(str(tmp), format="parquet")
    for nb in range(BUCKETS):
        tab = tmp_ds.to_table(filter=ds.field("bucket") == nb, columns=["nutrient_id","canon","amount_norm"])
        if tab.num_rows == 0: continue
        med = tab.to_pandas().groupby(["nutrient_id","canon"], as_index=False)["amount_norm"].mean() \
                              .groupby("nutrient_id", as_index=False)["amount_norm"].median()
        rows.extend({"nutrient_id":int(r.nutrient_id), "ref_amount_norm":float(r.amount_norm)} for r in med.itertuples(index=False))
    pd.DataFrame(rows).to_parquet(out_dir/"nutrient_scale.parquet", index=False)
    shutil.rmtree(tmp, ignore_errors=True)


# SCORING KERNEL - per-bacterium greedy food selection.
# Adapted from the query prototype's _run_lbl.
STATIC_DB, DYNAMIC_STATE = {}, {}

def _init_w(p, d):
    global STATIC_DB, DYNAMIC_STATE
    STATIC_DB = pickle.load(open(p, "rb")); DYNAMIC_STATE = d


def score_one_bacterium(lbl):
    """Score one bacterium. Returns (greedy_rows_df, candidate_rows_df).

    `greedy_rows_df` is the existing diversified `--max_foods` shortlist
    (absolute output, used by `bacteria2food` mode as-is).

    `candidate_rows_df` is a longer, clean-slate list of every (bacterium,
    food) pair scored, with `score` computed AS IF this food were the first
    pick (no greedy `cur` accumulation, no redundancy penalty). This is what
    the differential-ranking pass needs - peer median per food only makes
    sense if scores are independent of greedy-state.
    """
    all_t = set(int(x) for x in DYNAMIC_STATE["lbl2targ"].get(lbl, set()))
    targs = (all_t & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
    if not targs:
        return pd.DataFrame()
    need = set(targs)
    for t in targs:
        for p in DYNAMIC_STATE["prox"].get(t, []):
            if p[2] or (DYNAMIC_STATE["ams"] and DYNAMIC_STATE["amp"] and DYNAMIC_STATE["ndf"].get(t, 0) < 50000 and p[3] == "macro_proxy"):
                need.add(p[0])

    def _tau(n):
        ref = DYNAMIC_STATE["ref"].get(n, 0.0)
        scale_floor = 0.1 if STATIC_DB["unit_arr"][n] == 1.0 else 0.01
        ref_q = 0.01 * ref
        cap = 500.0 if STATIC_DB["unit_arr"][n] == 1.0 else 50.0
        return min(max(scale_floor, ref_q), cap)
    tau = {n: _tau(n) for n in targs}

    scn = ds.dataset(DYNAMIC_STATE["bdir"], format="parquet", partitioning="hive").scanner(
        columns=["fdc_id","nutrient_id","amount"],
        filter=ds.field("bucket").isin(sorted({n % BUCKETS for n in need})) & ds.field("nutrient_id").isin(sorted(need)),
        batch_size=50_000)
    tk_p, tk_o = TopKByNutrient(), TopKByNutrient()
    c_arr = STATIC_DB["canon_arr"]; p_arr = STATIC_DB["port_arr"]; pl_arr = STATIC_DB["plant_arr"]
    sp_arr = STATIC_DB["spice_arr"]; u_arr = STATIC_DB["unit_arr"]
    max_f, max_u = len(c_arr), len(u_arr)
    for b in scn.to_batches():
        f = b["fdc_id"].to_numpy().astype(np.int64, copy=False)
        n = b["nutrient_id"].to_numpy().astype(np.int32, copy=False)
        a = b["amount"].to_numpy().astype(np.float32, copy=False)
        m = (f>=0)&(f<max_f)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f, n, a = f[m], n[m], a[m]; c = c_arr[f]; v = (c != -1)
        if not DYNAMIC_STATE["asp"]: v &= ~sp_arr[f]
        if not v.any(): continue
        f, n, a, c = f[v], n[v], a[v], c[v]; ports = p_arr[f].copy()
        if DYNAMIC_STATE["asp"]: ports[sp_arr[f]] = 2.0
        a = a*(ports/50.0)*u_arr[n]
        idx = np.lexsort((a, c, n)); n, c, a, f = n[idx], c[idx], a[idx], f[idx]
        m2 = np.empty(len(n), dtype=bool); m2[-1] = True; m2[:-1] = (n[:-1] != n[1:]) | (c[:-1] != c[1:])
        n, c, a, ipl = n[m2], c[m2], a[m2], pl_arr[f[m2]]
        for i in range(len(n)):
            (tk_p if ipl[i] else tk_o).push(int(n[i]), int(c[i]), float(a[i]), K_PLANT if ipl[i] else K_OTHER)

    def get_tk(nid):
        d = dict(tk_p.val.get(nid, {}))
        for k, v in tk_o.val.get(nid, {}).items():
            if k not in d or v > d[k]: d[k] = v
        return d

    c2a = {}
    nd = STATIC_DB["nut_name"]
    for t in targs:
        th = float(tau.get(t, 0.0)); nn = nd.get(t, "")
        mc = 0.0001 if _FIBER.search(nn) else 0.001
        base_cut = th*mc; is_sp = (t in OLIGO_IDS or t in ISOFLAVONE_IDS)
        cut = base_cut
        for c, a in get_tk(int(t)).items():
            if a > cut: c2a.setdefault(c, {})[int(t)] = float(a)
        if not is_sp and DYNAMIC_STATE["ndf"].get(t, 0) < 5000:
            for p in [x for x in DYNAMIC_STATE["prox"].get(t, []) if x[2] or (DYNAMIC_STATE["ams"] and DYNAMIC_STATE["amp"] and x[3] == "macro_proxy")]:
                for c, a in get_tk(int(p[0])).items():
                    v = (float(a)/p[1])*0.25
                    if v > base_cut:
                        d = c2a.setdefault(c, {})
                        d[int(t)] = v if int(t) not in d or v > d[int(t)] else d[int(t)]
    c2a = {c: d for c, d in c2a.items() if d}
    del tk_p, tk_o; gc.collect()
    if not c2a:
        return pd.DataFrame()

    # PER-NUTRIENT WEIGHT - three factors, all of them properties of the NUTRIENT:
    #
    #   w[n] = (fibre boost) × (manual multiplier) × food-IDF
    #
    # A fourth factor, bacterial-IDF ^ spec_alpha, was removed: it demoted nutrients many
    # input bacteria target, which is backwards, since a guild-defining substrate is one
    # every member carries. On a 17-species panel with documented substrates, dropping it
    # raised MRR 0.252 -> 0.325.
    # COMMUNITY COVERAGE WEIGHT (community pass only; `cov_w` is absent otherwise).
    # Scoring the plain union counted a nutrient the same whether 1 organism or 58 could
    # use it, and the union saturates: 471 of 598 mappable nutrients by 6 months, so
    # 30 -> 58 taxa moved it 471 -> 474 and changed no ranking. Weighting each nutrient by
    # the share of members that can use it makes the score follow composition.
    cov_w = DYNAMIC_STATE.get("cov_w")
    cov_alpha = float(DYNAMIC_STATE.get("cov_alpha", 1.0))
    w = {n:(FIBER_WEIGHT if _FIBER.search(nd.get(n, "")) else 1.0)
            * float(DYNAMIC_STATE["nmult"].get(n, 1.0))
            * (1.0 + math.log(DYNAMIC_STATE["N"]/max(1, int(DYNAMIC_STATE["ndf"].get(n, 1)))))
            * ((cov_w.get(n, 0.0) ** cov_alpha) if cov_w else 1.0)
         for n in targs}

    if len(c2a) > MAX_CANONS:
        h = []
        for c, ct in c2a.items():
            s = sum(w.get(n, 1.0)*math.log1p(a/tau.get(n, 1e-9)) for n, a in ct.items() if tau.get(n, 0) > 0)
            if len(h) < MAX_CANONS: heapq.heappush(h, (s, c))
            elif s > h[0][0]: heapq.heapreplace(h, (s, c))
        c2a = {c: c2a[c] for _, c in h}

    cur = {n: 0.0 for n in targs}
    chosen = set(); seen_desc = set(); chosen_imps = []; chosen_cats = []; rows = []
    fs_db = STATIC_DB["food_stats"]
    ni_min = 1 if len(targs) <= 5 else min(5, max(2, int(0.03 * len(targs))))

    for rk in range(1, DYNAMIC_STATE["maxf"] + 1):
        bst = None
        for c, ct in c2a.items():
            fs = fs_db.get(c, {})
            if c in chosen or not fs or fs["dh"] in seen_desc: continue
            gn = 0.0; ni = 0; imp = {}
            for n, a in ct.items():
                th = tau[n]
                mcf = 0.0001 if _FIBER.search(nd.get(n, "")) else 0.001
                if a <= th*mcf: continue
                dg = math.log1p((cur[n] + a)/th) - math.log1p(cur[n]/th)
                if dg > 1e-12:
                    cg = w.get(n, 1.0)*dg; gn += cg; ni += 1; imp[n] = cg
            if ni < ni_min or gn < GAIN_MIN: continue
            imp = {n: v for n, v in imp.items() if v >= 0.05*gn}
            red = sum(sorted([cosine_sim(imp, p) for p in chosen_imps], reverse=True)[:3])/3.0 if chosen_imps else 0.0
            prm = min(1.0, max(0.0, ni/max(50, int(DYNAMIC_STATE["mc"].get(c, 0)))))
            geff = math.log1p(gn)*(math.sqrt(ni/(ni + 2.0)) if ni > 0 else 0.0)
            cat_rep = chosen_cats.count(fs.get("cat", "")); cat_pen = 0.20 * cat_rep
            # Phase 10 simplification - same math, 3 named terms instead of 7.
            #   food_baseline = flat per-food cost (regularizer + base cost + artifact)
            #   amount_cost   = breadth-vs-purity penalty that scales with ln(1+gn)
            #   greedy_extra  = redundancy + category-repeat penalty (greedy path only)
            food_baseline = fs["bc"] + ART_W * fs["art"]
            amount_cost   = (PROC_W * fs["pl"] + BROAD_W * ((1.0 - prm) ** BROAD_Q)) * math.log1p(gn)
            greedy_extra  = OVERLAP_W * red + cat_pen
            sc = geff * (1.0 + fs["tb"]) - food_baseline - amount_cost - greedy_extra
            # Back-compat: greedy-path observability fields expect these two
            # individual penalty names. Recompute them once for the bst tuple.
            proc_cost  = PROC_W * fs["pl"] * math.log1p(gn)
            broad_cost = BROAD_W * math.log1p(gn) * ((1.0 - prm) ** BROAD_Q)
            if sc > SCORE_MIN and (bst is None or sc > bst[0]):
                top_imp = sorted(imp.items(), key=lambda x: x[1], reverse=True)[:3]
                bst = (sc, gn, c, math.log1p(gn), ni, prm, fs["pl"], fs["art"], fs["tb"], geff,
                       proc_cost, broad_cost, red, imp,
                       ",".join(str(n) for n, _ in top_imp),
                       ",".join(nd.get(n, str(n)) for n, _ in top_imp), fs["dh"])
        if not bst: break
        c = bst[2]; chosen.add(c); seen_desc.add(bst[-1]); chosen_imps.append(bst[13])
        chosen_cats.append(fs_db[c].get("cat", ""))
        for n, a in c2a[c].items(): cur[n] += a
        cvs = {t: get_cov(cur[t], tau.get(t, 0.0)) for t in targs}
        cov_tot = sum(1 for v in cvs.values() if v >= 1.0)
        rows.append(dict(rank=rk, bacterium=lbl, representative_fdc_id=c,
                         food_name=fs_db[c].get("canon_name", ""),
                         n_variants=int(fs_db[c].get("n_variants", 1)),
                         score=bst[0], gain=bst[1], n_nutrients_improved=bst[4],
                         n_targets_total=len(targs), covered_targets_total=cov_tot,
                         coverage_total_frac=cov_tot/max(1, len(targs)),
                         top_nutrient_ids=bst[14], top_nutrient_names=bst[15]))

    # Clean-slate per-candidate scoring for differential mode
    # Same formula as the greedy loop but with cur=0 and red=0 - every food
    # is scored independently of selection order, so peer median per food
    # has well-defined semantics.
    if not DYNAMIC_STATE.get("score_pool", True):
        return pd.DataFrame(rows), pd.DataFrame()
    cand_rows = []
    for c, ct in c2a.items():
        fs = fs_db.get(c, {})
        if not fs:
            continue
        gn = 0.0; ni = 0; imp = {}
        for n, a in ct.items():
            th = tau[n]
            mcf = 0.0001 if _FIBER.search(nd.get(n, "")) else 0.001
            if a <= th*mcf:
                continue
            dg = math.log1p(a/th)  # cur is 0 → second term vanishes
            if dg > 1e-12:
                cg = w.get(n, 1.0)*dg; gn += cg; ni += 1; imp[n] = cg
        if ni < ni_min or gn < GAIN_MIN:
            continue
        imp_filt = {n: v for n, v in imp.items() if v >= 0.05*gn}
        prm = min(1.0, max(0.0, ni/max(50, int(DYNAMIC_STATE["mc"].get(c, 0)))))
        geff = math.log1p(gn)*(math.sqrt(ni/(ni + 2.0)) if ni > 0 else 0.0)
        dmode = DYNAMIC_STATE.get("diff_formula", "full")
        if dmode != "full":
            # The cost stack does two jobs.
            #   1. RANKING: dead weight. The stack is a FOOD property, so it shifts `score`
            #      and `peer_median` equally and cancels in score - peer_median. Zeroing
            #      proc_w / fiber_weight / broad_q leaves 85-94% of rankings byte-identical.
            #   2. ADMISSION: load-bearing. cat_penalty 5.0 (Alcoholic Beverages, Fast Foods)
            #      drives those negative, and the pre-median `score > 0` filter then drops
            #      them. Deleting the stack ("gain_only") readmits them: whole-plant foods
            #      91.1% -> 75.0% of differential top-10 rows, 21.4% from categories that
            #      were absent. The six-species panel still scores 3/3 and cannot see it.
            # "explicit_admission" keeps job 2 unchanged (same expression and threshold, so
            # the admitted set is identical) and drops job 1. art_w's 79% effect is then an
            # admission decision, not a scoring one.
            food_cost = (fs["bc"] + ART_W * fs["art"]
                         + (PROC_W * fs["pl"] + BROAD_W * ((1.0 - prm) ** BROAD_Q)) * math.log1p(gn))
            if dmode == "explicit_admission" and geff * (1.0 + fs["tb"]) - food_cost <= 0:
                continue
            sc = geff * (1.0 + fs["tb"])
            proc_cost = broad_cost = 0.0
        else:
            # Phase 10 simplification - same math as the 7-term form, 3 named terms.
            # No greedy-state terms (redundancy, category-repeat) - those are
            # greedy-path only and don't affect differential / perFood / community.
            food_baseline = fs["bc"] + ART_W * fs["art"]
            amount_cost   = (PROC_W * fs["pl"] + BROAD_W * ((1.0 - prm) ** BROAD_Q)) * math.log1p(gn)
            sc = geff * (1.0 + fs["tb"]) - food_baseline - amount_cost
            # Back-compat observability fields for cand_rows.
            proc_cost  = PROC_W * fs["pl"] * math.log1p(gn)
            broad_cost = BROAD_W * math.log1p(gn) * ((1.0 - prm) ** BROAD_Q)
        if sc <= SCORE_MIN:
            continue
        top_imp = sorted(imp_filt.items(), key=lambda x: x[1], reverse=True)[:3]
        cand_rows.append(dict(
            bacterium=lbl, representative_fdc_id=c,
            food_name=fs.get("canon_name", ""),
            n_variants=int(fs.get("n_variants", 1)),
            score=sc, gain=gn, n_nutrients_improved=ni,
            top_nutrient_ids=",".join(str(n) for n, _ in top_imp),
            top_nutrient_names=",".join(nd.get(n, str(n)) for n, _ in top_imp),
        ))
    return pd.DataFrame(rows), pd.DataFrame(cand_rows)


# INPUT LOADERS
# --- input column auto-detection (header-agnostic, multi-EC aware) -----------
# An EC cell: one or more EC numbers (optional "EC:" prefix), comma/semicolon
# separated.  Used both to spot the EC column and to detect a header row.
_EC_CELL_RE = re.compile(
    r"(?i)^(?:ec:?\s*)?\d+\.\d+\.\d+\.\d+"
    r"(?:\s*[,;]\s*(?:ec:?\s*)?\d+\.\d+\.\d+\.\d+)*$")
_TAXID_NAME_RE = re.compile(r"^\d+[_\s]\D")          # "47715_Lacticaseibacillus..."
_LOCUS_RE      = re.compile(r"^[A-Za-z]{1,4}[_-]?\d+$")   # gene/contig id: l_46433, k141_22
_ORG_WORD_RE   = re.compile(r"[A-Za-z]{3,}")         # an alphabetic word
_KNOWN_HEADERS = {"species", "ec_number", "ec", "strain", "bacterium", "taxon", "organism"}


def _col_frac(series: pd.Series, pred) -> float:
    """Fraction of non-empty values in `series` for which pred(value) is true."""
    vals = series.dropna().astype(str).str.strip()
    vals = vals[vals != ""]
    if len(vals) == 0:
        return 0.0
    return float(vals.map(lambda v: bool(pred(v))).mean())


def _species_likeness(series: pd.Series) -> float:
    """Heuristic score that a column holds organism/species names. The taxid-
    prefixed form (`<taxid>_Genus_species`) is unambiguous and scores highest;
    otherwise reward multi-word alphabetic names and penalize gene-locus columns
    and constant columns (e.g. a sample-id repeated on every row)."""
    vals = series.dropna().astype(str).str.strip()
    vals = vals[vals != ""]
    if len(vals) == 0:
        return 0.0
    taxid = float(vals.map(lambda v: bool(_TAXID_NAME_RE.match(v))).mean())
    if taxid >= 0.6:
        return 1.0 + taxid                       # taxid_Species - strongest signal
    org = float(vals.map(
        lambda v: (not _LOCUS_RE.match(v)) and len(_ORG_WORD_RE.findall(v)) >= 2
    ).mean())
    if vals.nunique() == 1:                      # constant column => sample id, not species
        org *= 0.1
    return org


def _detect_columns(df: pd.DataFrame):
    """Pick (species_col, ec_col) from a headerless frame by content."""
    sample = df.head(1000)
    ec_scores = {c: _col_frac(sample[c], _EC_CELL_RE.match) for c in df.columns}
    ec_col = max(ec_scores, key=ec_scores.get)
    if ec_scores[ec_col] < 0.3:
        raise ValueError(
            "Could not find an EC-number column in --mag_tsv. Expected a column of "
            "EC numbers like '3.2.1.1' (optionally comma-separated, e.g. "
            "'2.7.2.3,5.3.1.1'). Got columns: " + ", ".join(map(str, df.columns)))
    sp_candidates = [c for c in df.columns if c != ec_col] or list(df.columns)
    sp_scores = {c: _species_likeness(sample[c]) for c in sp_candidates}
    sp_col = max(sp_scores, key=sp_scores.get)
    return sp_col, ec_col


def _normalize_ec_frame(enz: pd.DataFrame) -> pd.DataFrame:
    """Given a frame with `species`, `ec_number` [, `strain`], split multi-EC
    cells into one EC per row, validate, dedupe, and return the canonical
    [bacterium, species, ec_number] frame consumed by the rest of the pipeline."""
    if "strain" not in enz.columns:
        enz["strain"] = ""
    enz = enz.copy()
    enz["strain"]  = enz["strain"].fillna("").astype("string").str.strip()
    enz["species"] = enz["species"].astype("string").str.strip()
    enz = enz.dropna(subset=["species", "ec_number"])
    # strip optional "EC:" prefixes, then split comma/semicolon EC lists -> rows
    enz["ec_number"] = (enz["ec_number"].astype("string")
                        .str.replace(r"(?i)\bec:?\s*", "", regex=True)
                        .str.split(r"\s*[,;]\s*", regex=True))
    enz = enz.explode("ec_number")
    enz["ec_number"] = enz["ec_number"].astype("string").str.strip()
    enz = enz[enz["ec_number"].str.match(EC_RE, na=False)]
    enz["bacterium"] = enz["species"] + enz["strain"].apply(lambda x: f" | {x}" if x else "")
    enz = enz.drop_duplicates(subset=["bacterium", "ec_number"])
    return enz[["bacterium", "species", "ec_number"]]


def load_user_ec(args) -> pd.DataFrame:
    """Return a normalized frame with columns: bacterium, species, ec_number.

    Accepts two input shapes, auto-detected:
      * Header TSV with `species` and `ec_number`/`ec` columns (+ optional
        `strain`) - the documented format.
      * Headerless annotation TSV (e.g. gene_annot/*_ec.tsv): no header, the
        species and EC columns are recognized by content. EC cells may list
        several comma/semicolon-separated EC numbers (split into one per row).
    The bacterium label keeps the species string verbatim (e.g. the
    `<taxid>_Genus_species` form), matching the rest of the pipeline.
    """
    path = args.mag_tsv
    # Peek at the first row to decide whether there is a header.
    head = pd.read_csv(path, sep="\t", header=None, nrows=1, dtype="string")
    first = [("" if pd.isna(v) else str(v)).strip() for v in head.iloc[0].tolist()]
    last_is_ec = bool(first) and bool(_EC_CELL_RE.match(first[-1]))
    has_header = (not last_is_ec) and any(c.lower() in _KNOWN_HEADERS for c in first)

    if has_header:
        enz = pd.read_csv(path, sep="\t", dtype="string")
        cols = {c.lower(): c for c in enz.columns}
        ec_col = cols.get("ec_number") or cols.get("ec") or enz.columns[-1]
        sp_col = cols.get("species") or cols.get("taxon") or cols.get("organism") or enz.columns[0]
        print(f"[*] Input has a header; species='{sp_col}', ec='{ec_col}'.", flush=True)
        enz = enz.rename(columns={ec_col: "ec_number", sp_col: "species"})
        cols_keep = ["species", "ec_number"] + (["strain"] if "strain" in enz.columns else [])
        enz = enz[cols_keep]
    else:
        enz = pd.read_csv(path, sep="\t", header=None, dtype="string")
        sp_col, ec_col = _detect_columns(enz)
        print(f"[*] Headerless input ({enz.shape[1]} cols); auto-detected "
              f"species=col{sp_col}, ec=col{ec_col}.", flush=True)
        enz = enz[[sp_col, ec_col]].rename(columns={sp_col: "species", ec_col: "ec_number"})

    return _normalize_ec_frame(enz)


# Bact_ec reference: cache the 3.2 GB TSV to a deduplicated parquet on first
# use, then read the parquet for all downstream lookups (complement_ec and
# food2bacteria reference pool).
def _split_taxid(species_str: str):
    """If species starts with `<digits>_<rest>`, return (taxid, rest with underscores->spaces).
    Otherwise (None, species_str).
    """
    s = species_str.strip()
    m = re.match(r"^(\d+)[_\s](.+)$", s)
    if m:
        return m.group(1), m.group(2).replace("_", " ").strip()
    return None, s


def _norm_name(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())


# Derived caches must not land in the export/deposit directory: the reference TSV lives
# there, and writing `<name>.parquet` beside it silently ships a build artifact with the
# published resource. Override with $BAC2FOOD_CACHE.
CACHE_DIR = Path(os.environ.get("BAC2FOOD_CACHE", "/data/bac2food/cache"))


def _newest_mtime(*sources) -> float:
    """Newest mtime among the given files/dirs (dirs scanned one level deep for parquet).

    Missing paths are skipped rather than raising: several inputs are optional
    (food_portion.csv in particular is frequently absent).
    """
    newest = 0.0
    for s in sources:
        if not s:
            continue
        p = Path(s)
        if p.is_dir():
            # Hive-partitioned store: the directory's own mtime does not move when a
            # part file is rewritten in place, so look at the parts themselves.
            for f in p.glob("*/*.parquet"):
                newest = max(newest, f.stat().st_mtime)
        elif p.exists():
            newest = max(newest, p.stat().st_mtime)
    return newest


def _is_stale(target: Path, *sources) -> bool:
    """True if `target` is missing or older than any of its inputs.

    Existence-only checks let a derived index outlive the data it was built from - the
    same silent-staleness failure that put the bact_ec parquet in the deposit directory.
    """
    if not target.exists():
        return True
    return target.stat().st_mtime < _newest_mtime(*sources)


def _index_identity(*paths) -> dict:
    """Fingerprint the inputs a derived index was built from.

    An mtime comparison answers "is the index older than its inputs?" but not "was it built
    from THESE inputs?". Point --food_nutrient at a different store whose files predate the
    cached index and every staleness check passes while the index describes a different food
    universe. That is not hypothetical: it silently corrupted three comparison runs during the
    fdc_id re-key, because the migrated store was written after the index and the canonical
    store was not, so runs against the canonical store reused the migrated store's index.
    """
    out: dict[str, list[int]] = {}
    for p in paths:
        if not p:
            continue
        q = Path(p)
        if q.is_dir():
            fs = sorted(q.rglob("*.parquet"))
            out[str(q.resolve())] = [len(fs), sum(f.stat().st_size for f in fs)]
        elif q.exists():
            out[str(q.resolve())] = [1, q.stat().st_size]
    return out


def _enforce_index_identity(idx_dir: Path, *sources) -> None:
    """Discard a cached index that was not built from `sources`.

    Strict on a missing stamp: an index of unknown provenance is treated as wrong, because the
    failure it guards against is silent and produces plausible numbers. Rebuilding is cheap
    relative to publishing results derived from another store's index.
    """
    import json
    stamp = idx_dir / "index_inputs.json"
    ident = _index_identity(*sources)
    prev = None
    if stamp.exists():
        try:
            prev = json.loads(stamp.read_text())
        except (ValueError, OSError):
            prev = None
    if prev == ident:
        return
    dropped = [p for p in idx_dir.glob("*") if p.name != "index_inputs.json"]
    for p in dropped:
        (shutil.rmtree(p) if p.is_dir() else p.unlink())
    if dropped:
        why = "was built from different inputs" if prev else "has no provenance stamp"
        print(f"[*] index at {idx_dir} {why} - discarded {len(dropped)} cached file(s); "
              f"rebuilding from {sources[0]}", flush=True)
    stamp.write_text(json.dumps(ident, indent=1, sort_keys=True))


def bact_ec_cache_path(bact_ec_tsv: str) -> Path:
    """Return the parquet cache path for a reference TSV, outside the exports directory.

    Keyed by the source file's stem plus a short digest of its directory, so references
    with the same filename in different locations cannot collide.
    """
    p = Path(bact_ec_tsv).resolve()
    key = hashlib.sha1(str(p.parent).encode("utf-8")).hexdigest()[:8]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{p.stem}.{key}.parquet"


def _ref_read_options(bact_ec_tsv: str) -> tuple[pacsv.ReadOptions, dict]:
    """Pick the reader layout for whichever bacteria -> EC reference we were given.

    Two layouts are supported, so the predictor tracks the exported resource instead
    of a private copy of it:

      * species_enzymes.tsv (eggNOG v7, current) - 6 columns WITH a header:
        tax_id, genus, species, strain, organism, ec_number. `organism` is the full
        strain-level name, which is what the legacy file put in its species column, so
        it is the one mapped onto "species" here; the matching ladder's prefix-of-2 step
        collapses it to Genus+species exactly as before.
      * bact_ec.tsv (eggNOG v6, legacy) - 4 columns, NO header:
        tax_id, species, kingdom, ec_number.

    eggNOG v7 dropped its direct EC annotation, so the v7 reference is built through
    KEGG KO (see eggnog/6.1_eggnog7_species_enzymes.py); the schema difference is the
    only thing that reaches the predictor.
    """
    with open(bact_ec_tsv, encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
    cols = first.rstrip("\n").split("\t")
    if cols[:1] == ["tax_id"] and "organism" in cols:          # species_enzymes.tsv
        names = ["tax_id", "genus", "species_binomial", "strain", "species", "ec_number"]
        skip = 1
    else:                                                       # legacy bact_ec.tsv
        names = ["tax_id", "species", "kingdom", "ec_number"]
        skip = 0
    ro = pacsv.ReadOptions(column_names=names, skip_rows=skip, block_size=64 * 1024 * 1024)
    return ro, {c: pa.string() for c in names}


def ensure_bact_ec_cache(bact_ec_tsv: str) -> Path:
    """One-time conversion of the bacteria -> EC reference to a deduplicated parquet
    keyed by (tax_id, species, ec_number). Accepts either reference layout
    (see `_ref_read_options`).

    Subsequent runs read the parquet (< 100 MB, < 1 s).
    """
    cache = bact_ec_cache_path(bact_ec_tsv)
    if cache.exists() and cache.stat().st_mtime >= Path(bact_ec_tsv).stat().st_mtime:
        return cache
    print(f"[*] Building bact_ec cache: {cache} (one-time, streaming {bact_ec_tsv})", flush=True)
    seen = {}  # (taxid, species) -> set(ec_number)
    read_opts, col_types = _ref_read_options(bact_ec_tsv)
    rdr = pacsv.open_csv(
        bact_ec_tsv,
        read_options=read_opts,
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=col_types),
    )
    n_chunks = 0
    try:
        while True:
            try: batch = rdr.read_next_batch()
            except StopIteration: break
            n_chunks += 1
            tids = batch.column("tax_id").to_pylist()
            sps  = batch.column("species").to_pylist()
            ecs  = batch.column("ec_number").to_pylist()
            for tid, sp, ec in zip(tids, sps, ecs):
                if not ec or not EC_RE.match(ec): continue
                if not sp: continue
                seen.setdefault((tid or "", sp), set()).add(ec)
            if n_chunks % 5 == 0:
                print(f"    chunks={n_chunks}  unique_species={len(seen)}", flush=True)
    finally:
        rdr.close()

    rows = [(tid, sp, ec) for (tid, sp), ecs in seen.items() for ec in ecs]
    df = pd.DataFrame(rows, columns=["tax_id","species","ec_number"])
    df["tax_id"] = df["tax_id"].astype("string")
    df.to_parquet(cache, index=False, compression="zstd")
    print(f"[*] Cached {len(df)} unique (tax_id, species, ec_number) rows → {cache}", flush=True)
    return cache


def load_bact_ec_reference(bact_ec_tsv: str,
                           filter_taxids: set | None = None,
                           filter_norm_names: set | None = None) -> pd.DataFrame:
    """Read the cached bact_ec parquet. Optional filters narrow the read via
    pyarrow predicate pushdown so we never materialize rows we don't need.
    Returns a frame with columns: tax_id, species, ec_number.
    """
    cache = ensure_bact_ec_cache(bact_ec_tsv)
    expr = None
    if filter_taxids:
        expr = ds.field("tax_id").isin(sorted(filter_taxids))
    df = ds.dataset(cache, format="parquet").to_table(filter=expr).to_pandas() if expr is not None \
         else pd.read_parquet(cache)
    if filter_norm_names:
        df = df[df["species"].astype(str).map(_norm_name).isin(filter_norm_names)]
    return df


def _prefix2(s: str) -> str:
    """Return the first two whitespace tokens of a normalized name (Genus species).
    Lets `_match_user_to_reference` collapse hep's strain-level names
    ("1134687_Klebsiella_michiganensis_strain") onto bact_ec's species-level
    names ("Klebsiella michiganensis"). Returns the original name unchanged if
    fewer than two tokens.
    """
    parts = s.split()
    return " ".join(parts[:2]) if len(parts) >= 2 else s


def match_user_to_reference(user_species: list[str], bact_ec_tsv: str) -> dict[str, set]:
    """Return {user_species_label -> reference EC set}. Matching ladder
    (tried in order, first hit wins per user species):
      1. Exact tax_id match (user must have leading `<taxid>_` prefix).
      2. Exact normalized species-name match.
      3. Genus + species prefix match (strips strain trailing tokens).

    Crucially, the ladder is applied PER USER SPECIES so collisions don't
    drop anyone: hep can carry both `1150423_Bifidobacterium_dentium` and
    `1689_Bifidobacterium_dentium` and both will resolve to the reference
    "Bifidobacterium dentium Bd1" rows.
    """
    ref = load_bact_ec_reference(bact_ec_tsv, filter_taxids=None, filter_norm_names=None)
    # Normalize per DISTINCT organism name (~10.7k), never per row. Adding `nname`/`pref`
    # columns instead materializes two Python strings for every one of the reference's
    # 20.5M rows, which is several GB and OOM-kills the process on a 16 GB machine - the
    # eggNOG v7 reference is twice the size of the v6 one this was written against.
    ref_by_taxid = ref.groupby("tax_id")["ec_number"].apply(set).to_dict()
    by_species   = ref.groupby("species")["ec_number"].apply(set).to_dict()
    del ref
    ref_by_name: dict[str, set] = {}
    ref_by_pref: dict[str, set] = {}
    for sp_name, ecs in by_species.items():
        nn = _norm_name(str(sp_name))
        ref_by_name.setdefault(nn, set()).update(ecs)
        ref_by_pref.setdefault(_prefix2(nn), set()).update(ecs)
    del by_species

    result: dict[str, set] = {}
    for sp in user_species:
        tid, plain = _split_taxid(sp)
        nname = _norm_name(plain)
        pref  = _prefix2(nname)
        ec_set = (
            (tid and ref_by_taxid.get(tid))
            or ref_by_name.get(nname)
            or ref_by_pref.get(pref)
            or set()
        )
        if ec_set:
            result[sp] = ec_set
    return result


def run_complement_ec(user_ec: pd.DataFrame, bact_ec_tsv: str, out_path: Path):
    """Diff user ECs against the reference per species, using the shared
    taxid → name → prefix-of-2 matching ladder.
    """
    print(f"[*] --complement_ec: querying {bact_ec_tsv}", flush=True)
    user_sp = user_ec.groupby("species")["ec_number"].apply(lambda s: set(s.dropna())).to_dict()
    ref_sets = match_user_to_reference(list(user_sp.keys()), bact_ec_tsv)
    n_matched = sum(1 for sp in user_sp if ref_sets.get(sp))
    print(f"    matched {n_matched} of {len(user_sp)} user species to reference "
          f"({n_matched/max(1,len(user_sp))*100:.1f}%)", flush=True)
    rows = []
    for sp, user_set in user_sp.items():
        rset = ref_sets.get(sp, set())
        missing = sorted(rset - user_set)
        rows.append(dict(species=sp,
                         n_user_ec=len(user_set),
                         n_reference_ec=len(rset),
                         n_missing_ec=len(missing),
                         missing_ec=";".join(missing)))
    pd.DataFrame(rows).sort_values("n_missing_ec", ascending=False).to_csv(out_path, sep="\t", index=False)
    print(f"[*] Wrote {out_path}", flush=True)


def augment_user_ec_from_reference(user_ec: pd.DataFrame, bact_ec_tsv: str,
                                    threshold: int = 0) -> pd.DataFrame:
    """Append rows for reference ECs that the user's annotation is missing.
    Each species's EC set becomes (user ECs) ∪ (matched reference ECs).

    `threshold` gates augmentation by user-EC count: a species is augmented
    only if it has fewer than `threshold` user ECs (genuinely under-annotated).
    `threshold=0` augments every matched species. Recommended >0 because
    blanket augmentation collapses all matched genus-mates onto the same
    ~2500-EC profile, homogenizing the rankings.
    """
    print(f"[*] --augment_with_reference: filling annotation gaps from {bact_ec_tsv}"
          + (f" (threshold={threshold} user ECs)" if threshold else ""), flush=True)
    user_ec = user_ec.copy()
    if "ec_source" not in user_ec.columns:
        user_ec["ec_source"] = "user"
    user_sp = user_ec.groupby("species")["ec_number"].apply(lambda s: set(s.dropna())).to_dict()
    ref_sets = match_user_to_reference(list(user_sp.keys()), bact_ec_tsv)
    added_rows = []
    n_added_total = 0
    n_skipped_well_annotated = 0
    for sp, user_set in user_sp.items():
        ref_set = ref_sets.get(sp, set())
        if not ref_set:
            continue
        if threshold > 0 and len(user_set) >= threshold:
            n_skipped_well_annotated += 1
            continue
        missing = ref_set - user_set
        if not missing:
            continue
        n_added_total += len(missing)
        bact_label = user_ec[user_ec["species"] == sp]["bacterium"].iloc[0]
        for ec in sorted(missing):
            added_rows.append(dict(bacterium=bact_label, species=sp, ec_number=ec, ec_source="reference"))
    if not added_rows:
        print(f"    no missing ECs to add; user annotation already supersets the matched reference.", flush=True)
        return user_ec
    add_df = pd.DataFrame(added_rows)
    n_aug_species = add_df["species"].nunique()
    skip_note = f"  ({n_skipped_well_annotated} matched species skipped - already ≥{threshold} ECs)" if threshold else ""
    print(f"    added {n_added_total} reference ECs across {n_aug_species} species "
          f"(median {add_df.groupby('species').size().median():.0f}/species).{skip_note}", flush=True)
    return pd.concat([user_ec, add_df], ignore_index=True)


# MODE OUTPUTS
_STATIC_FOOD_META_PATH: str | None = None  # set in main() before writers run

def attach_food_meta(df: pd.DataFrame, food_path: str, food_category_path: str,
                     static_food_meta_path: str | None = None) -> pd.DataFrame:
    """Merge FDC food.parquet metadata onto the score frame. Where the FDC
    food_category is empty (Phenol-Explorer / AFCD / some BioFoodComp sources
    don't populate food_category_id), fall back to the resolved category from
    static_food_meta.pkl which has the keyword-fallback applied.
    """
    if df.empty: return df
    fm = load_smart(food_path).drop_duplicates(subset=["fdc_id"], keep="last")
    # Normalize whitespace in description (same fix as build_static_food_meta) -
    # prevents embedded \n in CIQUAL / BioFoodComp descriptions from breaking
    # downstream TSV consumers.
    if "description" in fm.columns:
        fm["description"] = fm["description"].astype(str) \
            .str.replace(r"[\r\n\t]+", " ", regex=True) \
            .str.replace(r"\s{2,}", " ", regex=True) \
            .str.strip()
    if "food_category" in fm.columns: fm = fm.drop(columns=["food_category"])
    cat = load_smart(food_category_path).rename(columns={"id":"food_category_id","description":"food_category"})
    fm = fm.merge(cat.assign(food_category_id=cat["food_category_id"].astype(str)), on="food_category_id", how="left")
    fm["fdc_id"] = pd.to_numeric(fm["fdc_id"], errors="coerce").astype("Int64")
    out = df.merge(fm[["fdc_id","description","data_type","food_category"]],
                   left_on="representative_fdc_id", right_on="fdc_id", how="left")
    out = out.drop(columns=["fdc_id"], errors="ignore")
    # Fallback category from the resolved food_stats (build-time keyword pass).
    sm_path = static_food_meta_path or _STATIC_FOOD_META_PATH
    if sm_path and Path(sm_path).exists():
        sm = pickle.load(open(sm_path, "rb"))
        fid_to_cat = {f: s.get("cat","") for f, s in sm.get("food_stats", {}).items()}
        if fid_to_cat:
            fb = out["representative_fdc_id"].map(fid_to_cat)
            mask = out["food_category"].fillna("").eq("")
            out.loc[mask, "food_category"] = fb[mask].fillna("")
    return out


def write_community(community_scores: pd.DataFrame, n_contributing_per_food: dict,
                    out_prefix: Path, food_path: str, food_category_path: str,
                    keep_negative_scores: bool = False):
    """Write `<prefix>.community.tsv` - top foods for the whole microbiome as
    a single pseudo-bacterium whose EC set is the union of every input
    bacterium's effective ECs.
    """
    if community_scores.empty:
        print("[!] No community scores produced. Are any user ECs mapped to nutrients?", flush=True)
        return
    df = community_scores.copy()
    # Phase 10: drop net-negative-score foods. A negative score means the food's
    # penalty stack exceeds its gain - recommending it is misleading.
    if not keep_negative_scores:
        n0 = len(df)
        df = df[df["score"] > 0].copy()
        if len(df) < n0:
            print(f"[*] Community: filtered {n0-len(df)} negative-score rows "
                  f"({n0} -> {len(df)}); pass --keep_negative_scores to retain.", flush=True)
    df["n_contributing_bacteria"] = df["representative_fdc_id"].map(n_contributing_per_food).fillna(0).astype(int)
    out = attach_food_meta(df, food_path, food_category_path)
    lead = ["rank","food_name","representative_fdc_id","description","food_category",
            "n_variants","score","gain","n_nutrients_improved","n_targets_total",
            "covered_targets_total","coverage_total_frac",
            "n_contributing_bacteria",
            "top_nutrient_ids","top_nutrient_names","data_type"]
    cols = [c for c in lead if c in out.columns]
    round_floats(out[cols], 4).to_csv(out_prefix.with_suffix(".community.tsv"), sep="\t", index=False)
    print(f"[*] Wrote {out_prefix}.community.tsv  ({len(out)} rows)", flush=True)


def _diverse_head(g: pd.DataFrame, max_foods: int, cap: int) -> pd.DataFrame:
    """Take the top `max_foods` rows of one bacterium's ranking, but refuse to spend
    more than `cap` of them on any single food category or any single lead substrate.

    A differential table is read as a shortlist to eat from, and the practical shortlist
    is five to ten items - few meals have more ingredients. At that length redundancy is
    the dominant failure: the clean-slate ranking that feeds this table carries no
    complementarity term (unlike the greedy path, which already applies
    OVERLAP_W * redundancy + category-repeat), so its top five can be five ways of
    eating the same thing. Measured on a cohort sample, the differential top-5 spanned
    2.08 distinct food categories against the greedy path's 3.71.

    A cap rather than a penalty weight, deliberately: "at most `cap` foods from one
    category, at most `cap` sharing a lead substrate" is a rule a reader can apply by
    eye to the output, and it needs no scale calibration against a log-domain score.
    Rows are consumed in score order, so within the cap the ordering is untouched.
    """
    kept, n_cat, n_nut = [], {}, {}
    overflow = []
    for row in g.itertuples(index=False):
        cat = getattr(row, "food_category", None) or ""
        nut = (getattr(row, "top_nutrient_ids", "") or "").split(",")[0]
        if n_cat.get(cat, 0) >= cap or (nut and n_nut.get(nut, 0) >= cap):
            overflow.append(row)
            continue
        kept.append(row)
        n_cat[cat] = n_cat.get(cat, 0) + 1
        if nut:
            n_nut[nut] = n_nut.get(nut, 0) + 1
        if len(kept) >= max_foods:
            break
    # If the cap starved the list - few categories available, or a thin candidate set -
    # backfill in score order rather than returning a short table. The cap shapes the
    # shortlist; it must not silently shrink it.
    if len(kept) < max_foods:
        kept.extend(overflow[:max_foods - len(kept)])
    return pd.DataFrame(kept, columns=g.columns)


def write_differential_bacteria2food(cand_scores: pd.DataFrame, out_prefix: Path,
                                      max_foods: int, food_path: str, food_category_path: str,
                                      min_peers: int = 20,
                                      bact_diag: dict | None = None,
                                      keep_negative_scores: bool = False,
                                      rank_by: str = "comp_score",
                                      diversity: int = 0):
    """Rank each bacterium's foods by how much it beats its peers on that food.

    For each food F, the peer median is computed over all bacteria that
    produced a clean-slate score for F.

    Two ranking rules, selected by `rank_by`:

      comp_score  `comp_score[B,F] = score[B,F] - peer_median[F]`, ranked desc. The
                  original rule: the peer comparison is the score.

      score       the peer comparison is an ADMISSION TEST - keep F only where
                  `score[B,F] > peer_median[F]` - and rank the survivors on absolute
                  `score`. Same question ("where does B beat its peers?"), but the
                  answer is ordered by how good the food is for B rather than by the
                  size of the margin.

    The distinction is not cosmetic. Subtracting the peer median removes any substrate
    the peers share, which is exactly what defines a substrate GUILD: four starch
    degraders scored against each other cancel on starch and are ranked on whatever
    trace phytochemical happens to differ. On the 17-species biology panel the
    subtraction scores MRR 0.216 and the admission test 0.353, recovering starch and
    inulin organisms (B. adolescentis, R. bromii, Collinsella, F. prausnitzii) at or
    near rank 1 with plausibility unchanged (94.1% whole-plant, 0% junk).

    If a food has fewer than `min_peers` peer scores, its peer median is
    statistically thin - we mark `peer_n` so users can filter, but still
    rank it (the rare foods only a handful of bacteria can use are exactly
    the rows differential mode is meant to surface).
    """
    if cand_scores.empty:
        print("[!] No candidate scores produced for differential ranking.", flush=True)
        return
    df = cand_scores.copy()
    # Phase 10: drop net-negative-score (bacterium, food) pairs before peer
    # statistics - they pollute the median for foods where most bacteria score
    # negatively, and they're useless as recommendations regardless.
    if not keep_negative_scores:
        n0 = len(df)
        df = df[df["score"] > 0].copy()
        if len(df) < n0:
            print(f"[*] Differential: filtered {n0-len(df)} negative-score rows "
                  f"({n0} -> {len(df)}); pass --keep_negative_scores to retain.", flush=True)
        if df.empty:
            print("[!] No positive-score rows remain - differential.tsv will be empty.", flush=True)
            return
    peer_stats = df.groupby("representative_fdc_id")["score"].agg(
        peer_median="median", peer_mean="mean", peer_n="size"
    ).reset_index()
    df = df.merge(peer_stats, on="representative_fdc_id", how="left")
    df["comp_score"] = df["score"] - df["peer_median"]
    if rank_by == "score":
        n0 = len(df)
        df = df[df["comp_score"] > 0].copy()
        print(f"[*] Differential: peer test admitted {len(df)} of {n0} rows "
              f"(score > peer median); ranking survivors on absolute score.", flush=True)
        if df.empty:
            print("[!] No row beats its peer median - differential.tsv will be empty.", flush=True)
            return
    df = df.sort_values(["bacterium", rank_by], ascending=[True, False])
    # Food metadata is attached BEFORE truncation because the diversity cap reads
    # food_category. Attaching after would leave nothing to diversify on.
    df = attach_food_meta(df, food_path, food_category_path)
    df = df.sort_values(["bacterium", rank_by], ascending=[True, False])
    if diversity and diversity > 0:
        df = (df.groupby("bacterium", as_index=False, group_keys=False)
                .apply(lambda g: _diverse_head(g, max_foods, diversity))
                .copy())
    else:
        df = (df.groupby("bacterium", as_index=False, group_keys=False)
                .head(max_foods).copy())
    df["rank"] = df.groupby("bacterium").cumcount() + 1
    if bact_diag:
        for col in ("n_user_ec","n_ec_in_db","n_nutrients_targeted"):
            df[col] = df["bacterium"].map({b: d.get(col, 0) for b, d in bact_diag.items()}).fillna(0).astype(int)
    out = df
    lead = ["bacterium","rank","food_name","representative_fdc_id","description","food_category",
            "n_variants","comp_score","score","peer_median","peer_mean","peer_n",
            "gain","n_nutrients_improved","top_nutrient_ids","top_nutrient_names",
            "n_user_ec","n_ec_in_db","n_nutrients_targeted","data_type"]
    cols = [c for c in lead if c in out.columns]
    round_floats(out[cols], 4).to_csv(out_prefix.with_suffix(".differential.tsv"),
                                       sep="\t", index=False)
    print(f"[*] Wrote {out_prefix}.differential.tsv  ({len(out)} rows)", flush=True)


def write_perFood(scores: pd.DataFrame, out_prefix: Path, top_k: int, food_path: str, food_category_path: str,
                  keep_negative_scores: bool = False):
    if scores.empty:
        print("[!] No (bacterium, food) scores were produced.", flush=True)
        return
    # Phase 10: drop net-negative (bacterium, food) pairs - a perFood top-K
    # of bacteria for food F shouldn't include bacteria where F net-hurts them.
    if not keep_negative_scores:
        n0 = len(scores)
        scores = scores[scores["score"] > 0].copy()
        if len(scores) < n0:
            print(f"[*] perFood: filtered {n0-len(scores)} negative-score rows "
                  f"({n0} -> {len(scores)}); pass --keep_negative_scores to retain.", flush=True)
        if scores.empty:
            print("[!] No positive-score rows remain - perFood.tsv will be empty.", flush=True)
            return
    inv = (scores.sort_values(["representative_fdc_id","score"], ascending=[True, False])
                 .groupby("representative_fdc_id", as_index=False, group_keys=False)
                 .head(top_k)
                 .copy())
    inv["rank"] = inv.groupby("representative_fdc_id").cumcount() + 1
    out = attach_food_meta(inv, food_path, food_category_path)
    lead = ["representative_fdc_id","food_name","description","food_category","n_variants",
            "rank","bacterium","score","gain","n_nutrients_improved","top_nutrient_names","data_type"]
    cols = [c for c in lead if c in out.columns]
    round_floats(out[cols], 4).to_csv(out_prefix.with_suffix(".perFood.tsv"), sep="\t", index=False)
    print(f"[*] Wrote {out_prefix}.perFood.tsv  ({len(out)} rows)", flush=True)


def load_abundance(path: str | None, labels: set[str]) -> dict[str, float]:
    """Read relative abundances and match them to the predictor's organism labels.

    Deliberately has no fallback. An earlier version derived abundance from per-organism
    annotated-locus counts, which looks reasonable and is wrong: those counts saturate at
    each genome's gene complement, so they report how completely a genome was assembled,
    not how much of the community it is. On this cohort the locus rank-curve puts rank 10
    at 62% of rank 1 with a top-10 CV of 0.13-0.22 in every sample, where a real
    rank-abundance curve spans orders of magnitude. Silently substituting that for
    abundance would have produced confidently mislabelled output.

    Labels are matched exactly, then by the '<taxid>_Genus_species' form with the taxid
    stripped, then by 'Genus species' with underscores and spaces interchangeable.
    """
    if not path:
        raise SystemExit("[!] --community_weight abundance requires --abundance_tsv "
                         "(see --help; there is no locus-count fallback on purpose).")
    raw: dict[str, float] = {}
    df = load_smart(path)
    cols = {c.lower(): c for c in df.columns}
    sp_col = next((cols[c] for c in ("species", "bacterium", "taxon", "name") if c in cols), None)
    ab_col = next((cols[c] for c in ("abundance", "relative_abundance", "rel_abundance",
                                     "fraction", "coverage") if c in cols), None)
    if sp_col is None or ab_col is None:
        raise SystemExit(f"[!] {path}: need a species/bacterium column and an abundance "
                         f"column; found {list(df.columns)}")
    for sp, ab in zip(df[sp_col].astype(str), pd.to_numeric(df[ab_col], errors="coerce")):
        if pd.notna(ab) and ab > 0:
            raw[sp.strip()] = raw.get(sp.strip(), 0.0) + float(ab)

    def norm(x: str) -> str:
        x = re.sub(r"^\d+_", "", x.strip())
        return re.sub(r"[\s_]+", " ", x).lower()

    by_norm: dict[str, float] = {}
    for k, v in raw.items():
        by_norm[norm(k)] = by_norm.get(norm(k), 0.0) + v
    out: dict[str, float] = {}
    for lbl in labels:
        if lbl in raw:
            out[lbl] = raw[lbl]
        elif norm(lbl) in by_norm:
            out[lbl] = by_norm[norm(lbl)]
    print(f"[*] --abundance_tsv: matched {len(out)} of {len(labels)} organisms "
          f"from {len(raw)} abundance rows.", flush=True)
    if not out:
        raise SystemExit("[!] no organism matched the abundance table; check its species "
                         "naming against the --mag_tsv labels.")
    return out


def write_perBacterium(scores: pd.DataFrame, out_prefix: Path, max_foods: int, food_path: str,
                       food_category_path: str, keep_negative_scores: bool = False):
    """<prefix>.perBacterium.tsv - for each organism, the foods that suit it best.

    Answers "which foods would stimulate THIS microbe", which none of the other three
    tables does. differential.tsv ranks by comp_score, i.e. what an organism exploits
    better than its peers, so a food every organism uses well is ranked LAST there
    despite being an excellent food for this one. perFood.tsv cannot be inverted to
    recover it either: it keeps only the top --top_bacteria_per_food organisms per food,
    so an organism appears for ~8 foods rather than its full shortlist.

    This is the same greedy per-bacterium frame perFood.tsv inverts, written out
    un-inverted - no extra scoring pass.
    """
    if scores.empty:
        print("[!] No (bacterium, food) scores were produced.", flush=True)
        return
    if not keep_negative_scores:
        n0 = len(scores)
        scores = scores[scores["score"] > 0].copy()
        if len(scores) < n0:
            print(f"[*] perBacterium: filtered {n0-len(scores)} negative-score rows "
                  f"({n0} -> {len(scores)}); pass --keep_negative_scores to retain.", flush=True)
        if scores.empty:
            print("[!] No positive-score rows remain - perBacterium.tsv will be empty.", flush=True)
            return
    out = (scores.sort_values(["bacterium", "score"], ascending=[True, False])
                 .groupby("bacterium", as_index=False, group_keys=False)
                 .head(max_foods)
                 .copy())
    out["rank"] = out.groupby("bacterium").cumcount() + 1
    out = attach_food_meta(out, food_path, food_category_path)
    lead = ["bacterium", "rank", "food_name", "representative_fdc_id", "description",
            "food_category", "n_variants", "score", "gain", "n_nutrients_improved",
            "top_nutrient_names", "data_type"]
    cols = [c for c in lead if c in out.columns]
    round_floats(out[cols], 4).to_csv(out_prefix.with_suffix(".perBacterium.tsv"),
                                      sep="\t", index=False)
    print(f"[*] Wrote {out_prefix}.perBacterium.tsv  ({len(out)} rows, "
          f"{out['bacterium'].nunique()} organisms)", flush=True)


# Populated in main() before workers fork. Used by the per-bacterium
# diagnostics merged into `<prefix>.differential.tsv` (n_ec_in_db column).
DB_EC_SET: set = set()


MODES = ("community", "differential")


def _select_mode(argv: list[str]) -> tuple[str, list[str]]:
    """Read the leading subcommand, if there is one.

    The two views answer different questions, need different parameters and - since each
    is produced by its own scoring pass - do not need to be computed together:

        bac2food_predict.py community    --mag ... --out ... [community options]
        bac2food_predict.py differential --mag ... --out ... [differential options]

    Invoked with no subcommand the script keeps its historical behaviour, running both
    passes and writing every view. That path is deprecated but NOT broken: the deposited
    run scripts and the reproducibility harness invoke it, and silently changing what
    they produce would be worse than carrying one branch.
    """
    if argv and argv[0] in MODES:
        return argv[0], argv[1:]
    return "both", argv


def main():
    mode, argv = _select_mode(sys.argv[1:])
    print("==========================================================", flush=True)
    print(f"RUNNING 4_predict/bac2food_predict.py [{mode}]", flush=True)
    print("==========================================================", flush=True)
    if hasattr(mp, "set_start_method"):
        try: mp.set_start_method("forkserver", force=True)
        except RuntimeError: pass

    # Pre-parse --config so a custom parameters.yaml is applied BEFORE the path
    # arguments below read their defaults from the bound globals, and so worker
    # processes (which re-import this module) pick it up via the env var.
    _pre = argparse.ArgumentParser(add_help=False)
    _pre.add_argument("--config", default=None)
    _pre_args, _ = _pre.parse_known_args(argv)
    if _pre_args.config:
        os.environ["BAC2FOOD_PARAMS"] = str(Path(_pre_args.config).resolve())
        load_params(_pre_args.config)
        print(f"[*] Loaded parameters from {_pre_args.config}", flush=True)

    _DESC = {
        "community": "Which foods feed this whole microbiome? Scores the union of every "
                     "input organism's substrate targets as one pseudo-organism, weighting "
                     "each nutrient by the share of the community that can use it, and picks "
                     "foods greedily as a complementary set. Writes <prefix>.community.tsv.",
        "differential": "What does each organism use better than the organisms it co-occurs "
                        "with? Scores every food independently per organism, keeps the foods "
                        "that beat the peer median, and ranks those on the score itself. "
                        "Writes <prefix>.differential.tsv, plus <prefix>.perFood.tsv and "
                        "<prefix>.perBacterium.tsv, which fall out of the same pass.",
        "both": "Bacteria <-> Food predictor (no subcommand: legacy mode). Runs BOTH passes "
                "and writes every view. Prefer `bac2food_predict.py community ...` or "
                "`bac2food_predict.py differential ...`, which expose only the options that "
                "apply to that view and skip the other view's scoring pass.",
    }
    ap = argparse.ArgumentParser(
        prog=f"bac2food_predict.py{'' if mode == 'both' else ' ' + mode}",
        description=_DESC[mode],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("" if mode != "both" else
                "\nSubcommands:\n  community     foods for the whole microbiome\n"
                "  differential  foods each organism uses better than its neighbours\n"))
    ap.add_argument("--config", default=None,
                    help="Path to a parameters.yaml of constants (scoring weights, data "
                         "paths, tables). Default: parameters.yaml next to the script, or "
                         "$BAC2FOOD_PARAMS. CLI flags below still override individual paths.")
    ap.add_argument("--mag", "--mag_tsv", dest="mag_tsv", default=None,
                    help="Bacterium->EC TSV. Either a header TSV with columns "
                         "species, ec_number [, strain], or a headerless annotation "
                         "TSV (e.g. gene_annot/*_ec.tsv) - the species and EC columns "
                         "are auto-detected and EC cells may list several comma- or "
                         "semicolon-separated EC numbers. Required unless --use_reference is set.")
    ap.add_argument("--use_reference", dest="use_reference", action=argparse.BooleanOptionalAction, default=None,
                    help="Use bact_ec.tsv as the bacterial-EC source (in addition to --mag_tsv if both "
                         "given). Default: off (use only the user's MAGs). Set on to score against the "
                         "broader bacterial universe - useful for food→bacteria recommendations.")
    ap.add_argument("--ref_min_ec", type=int, default=20,
                    help="Skip reference species with fewer than this many ECs.")
    ap.add_argument("--ref_max_species", type=int, default=0,
                    help="Cap the number of reference species (0 = no cap).")
    ap.add_argument("--nutrient_to_ec", default=PATH_NUTRIENT_TO_EC)
    ap.add_argument("--food_nutrient", default=PATH_BUCKETED_DIR, help="bucketed food_nutrient parquet dir")
    ap.add_argument("--food", default=PATH_FOOD_CSV)
    ap.add_argument("--nutrient", default=PATH_NUTRIENT_CSV)
    ap.add_argument("--food_category", default=PATH_FOOD_CATEGORY_CSV)
    ap.add_argument("--food_portion", default=PATH_FOOD_PORTION_CSV)
    ap.add_argument("--nutrient_alias", default=PATH_NUTRIENT_ALIAS)
    ap.add_argument("--index_dir", default=PATH_INDEX_DIR)
    ap.add_argument("--out", "--out_prefix", dest="out_prefix", required=True)
    ap.add_argument("--max_foods", type=int, default=10,
                    help="Per-bacterium food scan budget AND output row cap for community + "
                         "differential. Default 10. It is NOT a quality knob for the "
                         "differential view: that view is fed the clean-slate frame, where "
                         "every candidate is scored independently, so there is no exploration "
                         "budget to widen - grading the top 5 at a fixed depth gives an "
                         "identical MRR 0.296 at 5, 10, 20, 50 and 100. Its one biological "
                         "channel is that it sets --differential_diversity, and 10 (cap 4) is "
                         "the best point on that curve: MRR 0.308 against 0.296 uncapped, "
                         "because spreading food categories lifts documented substrates into "
                         "the top 5. Below 10 the cap is too tight; above 20 it stops binding.")
    if mode in ("differential", "both"):
        ap.add_argument("--min_peers", "--diff_min_peers", dest="diff_min_peers", type=int, default=20,
                        help="Differential mode: warn when peer median is computed from fewer than this "
                             "many bacteria for a food (peer statistic is thin). Default 20.")
    if mode in ("differential", "both"):
        ap.add_argument("--top_bacteria_per_food", type=int, default=20)
    if mode in ("community", "both"):
        ap.add_argument("--weight", "--community_weight", dest="community_weight", choices=["membership", "abundance", "none"],
                        default="membership",
                        help="How community.tsv weights a nutrient. 'membership' (default) "
                             "weights it by the SHARE of organisms that can act on it; "
                             "'abundance' does the same but weights organisms by relative "
                             "abundance, and REQUIRES --abundance_tsv; 'none' is the "
                             "pre-2026-08 plain union, in which one organism counts the same "
                             "as all of them. The union saturates - an infant gut community "
                             "reaches 471 of 598 mappable nutrients by 6 months - so it cannot "
                             "distinguish communities and is kept only to reproduce old runs.")
        ap.add_argument("--abundance_tsv", default=None,
                        help="Relative abundances for --community_weight abundance: a TSV with "
                             "columns species/bacterium and abundance (any positive scale; it is "
                             "renormalised). Get it from the profiler you already ran on these "
                             "reads (MetaPhlAn, Bracken, coverM). There is deliberately NO "
                             "fallback to per-organism annotated-locus counts: those saturate at "
                             "each genome's gene complement rather than tracking biomass. "
                             "Measured on this cohort, rank 10 of the locus curve sits at 62%% of "
                             "rank 1 and the top-10 CV is 0.13-0.22 in every sample, where a real "
                             "rank-abundance curve spans orders of magnitude. Using it would "
                             "weight by assembly completeness and genome size while claiming to "
                             "weight by abundance.")
        ap.add_argument("--coverage_alpha", type=float, default=1.0,
                        help="Exponent on the community coverage weight. >1 sharpens the "
                             "preference for foods most of the community can use; <1 flattens "
                             "it toward the old union behaviour. Ignored when "
                             "--community_weight none.")
    ap.add_argument("--drop_category", action="append", default=[])
    ap.add_argument("--allow-spices", dest="allow_spices", action=argparse.BooleanOptionalAction, default=False)
    ap.add_argument("--rebuild-static-meta", action="store_true")
    ap.add_argument("--keep_negative_scores", action=argparse.BooleanOptionalAction, default=False,
                    help="Keep (bacterium, food) pairs with score <= 0 in the three output files. "
                         "Default off - a negative score means the food's penalty stack exceeds its "
                         "gain (the bacterium does not benefit), so the row is misleading as a "
                         "recommendation. Set on for debugging / inspecting the tail.")
    ap.add_argument("--complement_ec", action="store_true",
                    help="Diff user ECs vs the bacterial reference EC set (report only).")
    ap.add_argument("--augment_with_reference", action="store_true",
                    help="Fill annotation gaps: union reference ECs from bact_ec.tsv into the "
                         "user's EC set before scoring. Fixes under-annotated MAGs whose specialty "
                         "enzymes (e.g. fucosidases for HMOs) the annotator missed. The output "
                         "tags augmented rows with ec_source='reference'.")
    ap.add_argument("--augment_threshold", type=int, default=200,
                    help="Only augment species whose user-annotated EC count is below this "
                         "threshold. Above it, the annotation is considered rich enough that "
                         "blanket augmentation would homogenize rankings. Set to 0 to augment "
                         "every matched species; set very high to be conservative. Default 200.")
    ap.add_argument("--bact_ec_tsv", default=PATH_BACT_EC_REF,
                    help="Reference EC TSV for --complement_ec / --augment_with_reference")
    if mode in ("differential", "both"):
        ap.add_argument("--formula", "--differential_formula", dest="differential_formula",
                        choices=["full", "explicit_admission", "gain_only"], default="full",
                        help="How the differential view scores. 'full' (default) subtracts the "
                             "per-food cost stack from the score, so the stack both ranks and "
                             "gates. 'explicit_admission' ranks on gain alone and lets the stack "
                             "gate only - which is exactly equivalent while ranking is on "
                             "comp_score, because the stack is a food property that cancels in "
                             "score - peer_median (both score MRR 0.216 there). It stops being "
                             "equivalent under --differential_rank score, where nothing cancels and "
                             "the stack becomes load-bearing: 'full' then reaches hits@3 8/17 and "
                             "passes the panel's negative control, against 6/17 and a rank-9 false "
                             "positive for 'explicit_admission' (MRR 0.312 vs 0.325 - MRR is the "
                             "closer call, hits@3 the clearer one). "
                             "'gain_only' drops the stack entirely - simplest, but it readmits "
                             "alcohol, fast food and sausages, so it is kept for comparison only. "
                             "Community scoring is unaffected by all three.")
        ap.add_argument("--rank", "--differential_rank", dest="differential_rank", choices=["comp_score", "score"], default="score",
                        help="How the differential view ORDERS the foods it admits. 'score' "
                             "(default) treats the peer comparison as an ADMISSION TEST - keep the "
                             "food only where score > peer median - and ranks survivors on absolute "
                             "score. 'comp_score' is the historical rule: rank on the margin over "
                             "the peer median. Subtracting that median removes any substrate the "
                             "peers SHARE, which is what defines a guild: under 'comp_score' four "
                             "Bacteroides returned one identical food ranked on a trace flavonol, "
                             "and under 'score' three of them return their own documented substrate "
                             "on three different foods. Panel MRR 0.216 -> 0.325. The cost is "
                             "convergence: mean pairwise top-10 Jaccard across cohort organisms "
                             "rises 0.08 -> 0.21, because the old rule maximised apparent "
                             "distinctness by ranking on whatever happened to differ.")
        ap.add_argument("--diversity", "--differential_diversity", dest="differential_diversity", type=int, default=-1, metavar="K",
                        help="Cap how many of a bacterium's ranked foods may share one food "
                             "category or one lead substrate. -1 (default) scales the cap with the "
                             "list: max(2, round(0.4 * --max_foods)) - 2 at 5 foods, 4 at 10. "
                             "0 disables it (pure score order); a positive value sets it directly. "
                             "The table is read as a shortlist to eat from and the practical "
                             "shortlist is 5-10 items, where redundancy dominates: the clean-slate "
                             "scores feeding this table carry no complementarity term, unlike the "
                             "greedy path, so a cohort top-5 spanned 2.08 food categories against "
                             "the greedy path's 3.71 (3.08 with the cap on). The cap MUST scale - "
                             "a fixed K=2 costs nothing at 5 foods (99.1%% whole-plant either way) "
                             "but forces >=5 categories at 10 foods, dropping whole-plant to 85.0%% "
                             "and admitting junk, because spreading a longer list reaches further "
                             "down into worse categories. Rows are taken in score order, so within "
                             "the cap the ranking is unchanged, and the list is backfilled rather "
                             "than shortened if the cap cannot be met.")
    ap.add_argument("--jobs", type=int, default=max(1, os.cpu_count() or 1))
    args = ap.parse_args(argv)
    # Options belonging to the other view are absent from this parser by design; give them
    # their defaults so the shared pipeline below needs no per-mode attribute guards.
    for _k, _v in (("community_weight", "membership"), ("abundance_tsv", None),
                   ("coverage_alpha", 1.0), ("diff_min_peers", 20),
                   ("top_bacteria_per_food", 20), ("differential_formula", "full"),
                   ("differential_rank", "score"), ("differential_diversity", -1)):
        if not hasattr(args, _k):
            setattr(args, _k, _v)

    out_prefix = Path(args.out_prefix)
    idx_dir = Path(args.index_dir); idx_dir.mkdir(parents=True, exist_ok=True)
    # Identity first, mtimes second. The mtime checks below only ask whether the index is older
    # than its inputs; this asks whether it was built from them at all.
    _enforce_index_identity(idx_dir, args.food_nutrient, args.food, args.nutrient,
                            args.food_category, args.food_portion, args.nutrient_alias)
    sm_path = idx_dir / "static_food_meta.pkl"
    global _STATIC_FOOD_META_PATH
    _STATIC_FOOD_META_PATH = str(sm_path)
    if args.rebuild_static_meta and sm_path.exists(): sm_path.unlink()
    # Rebuild when missing OR older than any input it was derived from. An existence-only
    # check silently serves a stale food universe after a source refresh.
    _sm_sources = (args.food, args.nutrient, args.food_category, args.food_portion,
                   args.food_nutrient, args.nutrient_alias)
    if _is_stale(sm_path, *_sm_sources):
        if sm_path.exists():
            print(f"[!] {sm_path.name} is older than its inputs - rebuilding.", flush=True)
        build_static_food_meta(args, sm_path)

    # Decide whether to pull bacteria from the bact_ec reference. Default off
    # for all modes - uses only the user's MAGs. Pass --use_reference to also
    # include the broader bact_ec universe.
    use_ref = bool(args.use_reference) if args.use_reference is not None else False
    if args.mag_tsv is None and not use_ref:
        ap.error("Either --mag_tsv or --use_reference must be provided.")

    # --- assemble the input EC frame ---
    parts = []
    if args.mag_tsv is not None:
        u = load_user_ec(args)
        if not u.empty:
            parts.append(u)
            print(f"[*] Loaded {len(u)} EC rows for {u['bacterium'].nunique()} bacteria from --mag_tsv.", flush=True)
    if use_ref:
        ref = load_bact_ec_reference(args.bact_ec_tsv)
        # Filter low-coverage / capped species.
        counts = ref.groupby("species")["ec_number"].nunique()
        keep_species = counts[counts >= args.ref_min_ec].index
        if args.ref_max_species > 0:
            keep_species = counts.loc[keep_species].sort_values(ascending=False).head(args.ref_max_species).index
        ref = ref[ref["species"].isin(keep_species)].copy()
        ref["strain"] = ""
        ref["bacterium"] = "[ref] " + ref["species"].astype(str)
        ref = ref[["bacterium","species","ec_number"]]
        parts.append(ref)
        print(f"[*] Loaded {len(ref)} EC rows for {ref['bacterium'].nunique()} reference bacteria from bact_ec.", flush=True)

    if not parts:
        print("[!] No valid EC rows in input. Exiting.", flush=True); return
    user_ec = pd.concat(parts, ignore_index=True).drop_duplicates(subset=["bacterium","ec_number"])
    print(f"[*] Total: {len(user_ec)} EC rows for {user_ec['bacterium'].nunique()} bacteria.", flush=True)

    if args.augment_with_reference:
        user_ec = augment_user_ec_from_reference(user_ec, args.bact_ec_tsv,
                                                  threshold=args.augment_threshold)
        user_ec = user_ec.drop_duplicates(subset=["bacterium","ec_number"], keep="first")
        print(f"[*] After --augment_with_reference: {len(user_ec)} EC rows.", flush=True)

    # --- nutrient targets per bacterium ---
    nut_df = load_smart(args.nutrient)
    blocked = set(nut_df[nut_df["name"].astype(str).str.contains(_BLOCKED, regex=True, na=False)]["id"].dropna().astype(int).tolist())
    n2ec = pd.read_csv(args.nutrient_to_ec, sep="\t")
    global DB_EC_SET
    DB_EC_SET = set(n2ec["ec_number"].dropna().astype(str).unique())
    merged = user_ec.merge(n2ec, on="ec_number")
    lbl2targ = {lbl: set(sub["nutrient_id"].astype(int).unique()) for lbl, sub in merged.groupby("bacterium")}
    print(f"[*] {len(lbl2targ)} bacteria have at least one EC mapped to a nutrient.", flush=True)

    # --- proximity map (same logic as the query prototype) ---
    prox = {}
    if args.nutrient_alias and Path(args.nutrient_alias).exists():
        a = pd.read_csv(args.nutrient_alias, sep="\t")
        div_map = a[a["rule"].str.strip().str.lower().isin(DANGEROUS_RULES)].groupby(["generic_nutrient_id","rule"])["specific_nutrient_id"].nunique().to_dict()
        for g, s, r in a[["generic_nutrient_id","specific_nutrient_id","rule"]].itertuples(index=False):
            so, k = classify_rule(str(r).strip().lower())
            prox.setdefault(int(s), []).append((int(g), float(div_map.get((int(g), r), 1.0)), so, k))
    for t in [1017,1021,1022,1403,2058,1071,1019]: prox.setdefault(t, []).append((1079, 1.0, True, "form"))
    for t in [1015,1016,1020]:                    prox.setdefault(t, []).append((1009, 1.0, True, "form"))
    for t in [1181,1182,1042]:                    prox.setdefault(t, []).append((99999, 1.0, True, "form"))

    # Aliases for bacterial substrates from extra_bacterial_seeds.tsv that don't
    # have their own FDC nutrient_id but ARE chemically equivalent to one (or
    # are a depolymerization product of one). Each tuple (specific_id, generic_id)
    # tells the scoring kernel: "when looking for specific_id in a food, also
    # check generic_id and substitute at 25% efficiency". This is the same
    # mechanism the original pipeline uses for vitamin / acid-base / form
    # aliases - inherited from the query prototype.
    #
    #   200001 GlcNAc            ← FDC 96310 CHITIN          (chitin → GlcNAc monomer)
    #   200002 GalNAc            ← FDC 96310 CHITIN          (loose proxy: chitin family)
    #   200007 Xylan             ← FDC 1019  Pentosan,
    #                             FDC 1021  Hemicellulose,
    #                             FDC 1074  Xylose
    #   200008 Arabinoxylan      ← FDC 1019  Pentosan,
    #                             FDC 1073  Arabinose
    #   200009 Dextran           ← FDC 1069  Oligosaccharides,
    #                             FDC 2064  Oligosaccharides (alt FDC id)
    #   200010 Cellobiose        ← FDC 1022  Cellulose
    #   200011 Pullulan          ← FDC 1009  Starch
    bacterial_aliases = [
        (200001, 96310), (200002, 96310),
        (200007, 1019),  (200007, 1021), (200007, 1074),
        (200008, 1019),  (200008, 1073),
        (200009, 1069),  (200009, 2064),
        (200010, 1022),
        (200011, 1009),
    ]
    for specific, generic in bacterial_aliases:
        prox.setdefault(specific, []).append((generic, 1.0, True, "form"))

    # --- build modeled index if missing ---
    # Derived from the bucketed parquet and the static meta, so it is stale if either is
    # newer than it - not merely if it is absent.
    if any(_is_stale(idx_dir/f, args.food_nutrient, sm_path)
           for f in ["modeled_totals.parquet", "nutrient_df.parquet", "nutrient_scale.parquet"]):
        build_modeled_index(args.food_nutrient, pickle.load(open(sm_path, "rb")), idx_dir)
    tot_pq = pd.read_parquet(idx_dir/"modeled_totals.parquet")


    dyn = {
        "bdir": args.food_nutrient, "nmult": {}, "lbl2targ": lbl2targ,
        "prox": {k: [(p[0], float(p[1]) if float(p[1]) > 0 else 1.0, p[2], p[3]) for p in v] for k, v in prox.items()},
        "mc": tot_pq.set_index("canon")["model_count"].to_dict(),
        "mm": tot_pq.set_index("canon")["model_mass"].to_dict(),
        "ndf": pd.read_parquet(idx_dir/"nutrient_df.parquet").set_index("nutrient_id")["df_foods"].to_dict(),
        "N":   len(tot_pq),
        "ref": pd.read_parquet(idx_dir/"nutrient_scale.parquet").set_index("nutrient_id")["ref_amount_norm"].to_dict(),
        "maxf": args.max_foods, "ams": ALLOW_MACRO_SCAN, "amp": ALLOW_MACRO_PROXY,
        "asp":  args.allow_spices, "blocked": blocked,
        "score_pool": True,   # required for the differential output
        # Differential formula selector. Off = the historical food-cost stack; on = the
        # cancellation-aware form. Kept as a switch rather than a replacement so the two
        # can be A/B-ed against the six-species biology gate on identical inputs.
        "diff_formula": args.differential_formula,
    }

    # --- per-bacterium scoring (drives differential.tsv + perFood.tsv) ---
    # Skipped entirely in community mode: the community view is produced by its own single
    # pass over a pseudo-organism, so scoring every organism separately would be wasted work.
    # Hoisted out of the per-bacterium branch: community mode skips that branch entirely and
    # still needs a context for its own single-worker pool.
    ctx = mp.get_context("forkserver") if hasattr(mp, "get_context") else mp
    if mode == "community" or not lbl2targ:
        all_scores = pd.DataFrame()
        cand_scores = pd.DataFrame()
    else:
        ordered = sorted(lbl2targ.keys(), key=lambda k: len(lbl2targ[k]), reverse=True)
        with ctx.Pool(args.jobs, initializer=_init_w, initargs=(sm_path, dyn), maxtasksperchild=1) as pool:
            results = pool.map(score_one_bacterium, ordered)
        greedy = [r[0] for r in results if isinstance(r, tuple) and not r[0].empty]
        cands  = [r[1] for r in results if isinstance(r, tuple) and not r[1].empty]
        all_scores  = pd.concat(greedy, ignore_index=True) if greedy else pd.DataFrame()
        cand_scores = pd.concat(cands,  ignore_index=True) if cands  else pd.DataFrame()

    # --- community scoring (drives community.tsv) --- skipped in differential mode.
    # Treat the union of every input bacterium's targets as one pseudo-
    # bacterium "[community]".
    community_scores = pd.DataFrame()
    n_contributing_per_food: dict = {}
    if mode != "differential" and lbl2targ:
        community_targs = set().union(*lbl2targ.values())
        # Per-target count: how many of the user's bacteria contribute it.
        target_to_n_bact: dict[int, int] = {}
        for targs in lbl2targ.values():
            for t in targs:
                target_to_n_bact[t] = target_to_n_bact.get(t, 0) + 1
        # Coverage weight per nutrient: what SHARE of the community can act on it.
        #   membership  each organism counts once. Answers "how many members benefit".
        #   abundance   organisms weighted by their share of annotated loci, the only
        #               abundance proxy the annotation input carries. Answers "how much of
        #               the community, by mass, benefits" - closer to the biology, but it
        #               inherits any assembly/annotation depth bias per organism.
        #   none        the pre-2026-08 plain union. Kept so old runs reproduce.
        cov_w = None
        if args.community_weight != "none":
            if args.community_weight == "abundance":
                ab = load_abundance(args.abundance_tsv, set(lbl2targ))
                tot = sum(ab.values()) or 1.0
                wt = {l: ab.get(l, 0.0) / tot for l in lbl2targ}
                miss = [l for l in lbl2targ if l not in ab]
                if miss:
                    print(f"[!] --abundance_tsv: {len(miss)} of {len(lbl2targ)} organisms have "
                          f"no abundance and contribute nothing to the community weight "
                          f"(e.g. {miss[0]}).", flush=True)
            else:
                wt = {l: 1.0 / max(1, len(lbl2targ)) for l in lbl2targ}
            cov_w = {}
            for lbl, targs in lbl2targ.items():
                for t in targs:
                    cov_w[int(t)] = cov_w.get(int(t), 0.0) + wt[lbl]
            share = sorted(cov_w.values())
            print(f"[*] community_weight={args.community_weight}: {len(cov_w)} nutrients, "
                  f"carrier share median {share[len(share)//2]:.3f}, max {share[-1]:.3f} "
                  f"(alpha {args.coverage_alpha})", flush=True)

        community_dyn = dict(dyn)
        community_dyn["lbl2targ"] = {"[community]": community_targs}
        community_dyn["score_pool"] = False
        community_dyn["cov_w"] = cov_w
        community_dyn["cov_alpha"] = args.coverage_alpha
        with ctx.Pool(1, initializer=_init_w, initargs=(sm_path, community_dyn), maxtasksperchild=1) as pool:
            cres = pool.map(score_one_bacterium, ["[community]"])
        if cres and isinstance(cres[0], tuple) and not cres[0][0].empty:
            community_scores = cres[0][0].drop(columns=["bacterium"], errors="ignore")
        # n_contributing_bacteria per food: how many input bacteria have at
        # least one of the nutrients that this food covers in our top-N.
        if not community_scores.empty:
            # Per-food top-N nutrient ids string → set of ints.
            for _, row in community_scores.iterrows():
                nids = {int(x) for x in str(row.get("top_nutrient_ids","")).split(",") if x.strip().isdigit()}
                n_contributing_per_food[int(row["representative_fdc_id"])] = max(
                    (target_to_n_bact.get(n, 0) for n in nids), default=0)

    # --- per-bacterium diagnostics merged into differential output ---
    ec_per_bact = user_ec.groupby("bacterium")["ec_number"].apply(lambda s: set(s)).to_dict()
    bact_diag = {
        b: {"n_user_ec": len(ec_set),
            "n_ec_in_db": len({e for e in ec_set if e in DB_EC_SET}),
            "n_nutrients_targeted": len(lbl2targ.get(b, set()))}
        for b, ec_set in ec_per_bact.items()
    }

    # --- write the result files for this mode ---
    if mode != "differential":
        write_community(community_scores, n_contributing_per_food, out_prefix, args.food,
                        args.food_category, keep_negative_scores=args.keep_negative_scores)
    if mode == "community":
        if args.complement_ec:
            run_complement_ec(user_ec, args.bact_ec_tsv, out_prefix.with_suffix(".complement_ec.tsv"))
        print("[*] Done.", flush=True)
        return
    write_differential_bacteria2food(
        cand_scores, out_prefix, args.max_foods,
        args.food, args.food_category,
        min_peers=args.diff_min_peers,
        bact_diag=bact_diag,
        keep_negative_scores=args.keep_negative_scores,
        rank_by=args.differential_rank,
        diversity=(max(2, round(0.4 * args.max_foods))
                   if args.differential_diversity < 0 else args.differential_diversity),
    )
    write_perFood(all_scores, out_prefix, args.top_bacteria_per_food, args.food, args.food_category,
                  keep_negative_scores=args.keep_negative_scores)

    write_perBacterium(all_scores, out_prefix, args.max_foods, args.food, args.food_category,
                       keep_negative_scores=args.keep_negative_scores)

    # --- optional: complement_ec diagnostic ---
    if args.complement_ec:
        run_complement_ec(user_ec, args.bact_ec_tsv, out_prefix.with_suffix(".complement_ec.tsv"))

    print("[*] Done.", flush=True)


if __name__ == "__main__":
    main()
