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
import pyarrow.parquet as pq
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

# The food set and the nutrient set this script scores against must be the same
# ones 5_export/export_resources.py publishes, or the predictor recommends foods
# and nutrients the released table does not contain. Both read the same bucketed
# store, so the shared exclusion policy has to be applied by both readers rather
# than by whichever one happens to run first.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "food_DBs"))
from _common.non_nutrients import (COPY_RECORD_RE, NON_NUTRIENT_IDS,  # noqa: E402
                                   find_exact_relistings, relisting_candidates,
                                   source_of_bucket_file)

EC_RE   = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
# A spice is eaten by the gram against the 100 g everything is measured on.
# Only read under --allow-spices; the category is excluded by default.
SPICE_SERVING_WEIGHT = 0.02   # 2 g serving / 100 g basis

BUCKETS = 256

# Names bound from parameters.yaml by apply_params(); declared up-front so module
# code and `global` statements resolve them. Real values are set at import below.
PATH_NUTRIENT_TO_EC = PATH_NUTRIENT_CSV = PATH_FOOD_CSV = PATH_FOOD_CATEGORY_CSV = None
PATH_BUCKETED_DIR = PATH_INDEX_DIR = PATH_NUTRIENT_ALIAS = PATH_BACT_EC_REF = None
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
# and the two readers of that canon aggregate it differently. score_one_bacterium
# takes the MAXIMUM over the members for each nutrient independently, so a canon
# is scored at its richest member; build_modeled_index takes the MEAN over the
# members that actually measured each one, for the modeled tables. An over-merged
# canon is therefore wrong twice over - the score inherits the most extreme
# member, the modeled value blends foods that are not the same food - which is
# what every composition axis below exists to prevent.
# The body must not cross a comma. With "," inside the character class the lazy
# quantifier still found the LEFTMOST start that could complete the match, so
# "Cheese, cheddar, natural shredded sharp, store brand, GREAT VALUE (CA1,NE) -
# NFY120WVO" lost everything from ", cheddar" onwards and 252 rows of branded
# cheddar landed in the bare 'cheese' canon. The sampling parenthetical is the
# one place a comma is allowed, and it is matched explicitly.
_NF_SUFFIX_RE = re.compile(
    # [A-Z0-9], not [A-Z]: the brand is stripped BEFORE this runs, so what is
    # left in front of the code can start with a digit ("Parmesan cheese,
    # grated, KRAFT 100% (AL,CA1) - NFY120DQP" becomes ", 100% (AL,CA1) - NF...")
    r"\s*[-,]\s*[A-Z0-9][\w &/%.'-]*(?:\s*\([^)]*\))?\s*-?\s*\bNF\w+\s*$", re.I)
# FDC sample codes come in two prefixes, not one: 6,096 rows end in an NF code and
# 446 in a CY code ("Broccoli, raw (IN1,NY1) - CY0906E"). Only NF was stripped, so
# the CY rows kept the code and produced canons like "rice - cy120co". CY requires
# a digit so the pattern cannot swallow an ordinary word that happens to start CY.
_NF_BARE_RE   = re.compile(r"\s*[-,]\s*\b(?:NF\w+|CY\w*\d\w*)\s*$", re.I)
# On an FDC lab row the ALL-CAPS chunks behind the food are the brand and the
# product line ("Salsa, PACE CHUNKY, MEDIUM", "Cheddar cheese, sliced, store
# brand, CRYSTAL FARMS & SHULLSBURG WISCONSIN"). Applied ONLY where an NF or CY
# sample code was actually present, because outside those rows an ALL-CAPS chunk
# can be the food itself - "BURGER KING - HAM" is ham, and "OIL, OLIVE, EXTRA
# LIGHT" is olive oil. The first chunk is never taken.
_CAPS_CHUNK_RE = re.compile(r",\s*[A-Z][A-Z0-9 &/'.%-]*[A-Z0-9][A-Z0-9 &/'.%-]*\s*$")


# "Starch," heads 60 FDC lab rows the same way "Minerals," does - but AFCD files
# a real food as "Starch, potato", and _PANEL_HEADS deliberately leaves it out
# for exactly that reason. The sample code is the discriminator, as it already
# is for "Sugars," and "Sweets,".
_CODED_PANEL_RE = re.compile(r"^\s*(?:starch|sugars?|sweets?)\s*,\s*", re.I)


def _strip_caps_brand(d: str) -> str:
    for _ in range(4):
        nd = _CAPS_CHUNK_RE.sub("", d)
        if nd == d or "," not in nd:
            return nd if nd.strip(" ,") else d
        d = nd
    return d
_PREP_RE      = re.compile(
    r"\b(?:"
    # the method must be consumed WHOLE: stripping only "roasted" from
    # "oil-roasted" left a bare "oil" chunk, which reads as peanut OIL.
    r"(?:oil|dry|honey|dark|light)[-\s]roasted|"
    r"raw|cooked|uncooked|boiled|baked|broiled|grilled|roasted|toasted|stewed|steamed|saut[ée]ed|simmered|"
    # braised, poached and barbecued are cooking methods like the rest; without
    # them FDC's lab rows kept "braised" as a chunk of the food's name
    r"braised|poached|barbecued|"
    r"fried|deep[-\s]?fried|pan[-\s]?fried|pan[-\s]?broiled|stir[-\s]?fried|battered|breaded|"
    r"oven|microwaved|"
    # how it is HEATED, not what it is - and FDC spells it both ways, which
    # split 'macaroni and cheese, microwavable' from the "microwaveable" one
    r"microwave(?:a)?ble|"
    r"(?:sun|freeze|spray|air|oven|vacuum)[-\s]?dried|"
    # "tinned" is the UK/NEVO spelling of "canned" and was never stripped, so
    # 50 canons carried it ("apricots tinned, sweetened") and stood apart from
    # the "canned" twin they name the same food as. It also blocked the plural
    # head fold, which is why those canons still read "apricots" and "cherries".
    r"frozen|fresh|dried|dehydrated|reconstituted|canned|tinned|drained|undrained|rinsed|dry|"
    r"unprepared|prepared|reheated|chilled|shelf\s+stable|nfs|"
    r"cut|sliced|chopped|diced|crinkle\s*cut|julienned|grated|shredded|cubed|halved|quartered|"
    r"halves|quarters|pieces|chunks|wedges|spears|florets|"
    # same guard as "seeds": NEVO's "Seeds and kernels" is one aggregate food
    r"(?<!and )(?<!or )kernels?(?!\s+(?:and|or)\b)"
    # ...nor when the kernel is what the product is PRESSED from. Palm KERNEL
    # oil is lauric where palm oil is palmitic - a different fat entirely - and
    # WAFCT's fortified palm kernel oil was landing in the palm oil canon.
    r"(?!\s+(?:oil|flour|meal|cake|butter|milk)\b)|"
    r"with\s+batter|with\s+sauce|"
    # "seeds?" is stripped as FDC's legume convention ("Beans, black, mature
    # seeds"), but not when it is half of a compound name: CIQUAL's "Seeds and
    # peanuts, dried" canonicalised to 'and peanut, dried', and a bread "with
    # pumpkin seeds" lost the seeds and kept "with pumpkin".
    r"mature\s+seeds?|"
    r"(?<!sesame )(?<!sunflower )(?<!pumpkin )(?<!poppy )(?<!hemp )(?<!chia )"
    r"(?<!flax )(?<!linseed )(?<!mustard )(?<!caraway )(?<!fennel )(?<!nigella )"
    r"(?<!squash )(?<!melon )(?<!cotton )(?<!niger )(?<!fenugreek )"
    r"seeds?(?!\s+(?:and|or)\b)"
    # ...nor when the seed is what the product is MADE of: Fineli's "Seed
    # Bread" came out as bare 'bread', which is a different food.
    r"(?!\s+(?:bread|loaf|roll|bun|mix|blend|oil|butter|paste|cracker|cake)\b)|"
    r"young|mature|immature|ripe|unripe|overripe|"
    r"half[-\s]?salted|demi[-\s]?sel|with(?:out)?\s+salt\s+added|salt\s+added|"
    r"no\s+salt\s+added|with(?:out)?\s+salt|salted|unsalted|"
    r"(?:meat\s+)?with(?:out)?\s+skin\s+and\s+(?:stone|pit|seeds?|core)s?|"
    r"(?:and\s+)?(?:stone|pit|core|peel|skin)s?\s+removed|"
    # STFCJ welds the claim to the word it qualifies - "thigh, meat with skin",
    # "breast, meat without skin" - and taking only the "with skin" half left
    # "meat" standing as a chunk of the name, on nine canons. "Meat with skin"
    # is the whole bird, which is what a bare "chicken" is; "meat without skin"
    # is the flesh and the part axis has already read it off the probe.
    r"(?:meat\s+)?with(?:out)?\s+skin|peeled|unpeeled|with(?:out)?\s+bones?|boneless|skinless|"
    # FDC's poultry phrasing for "without skin". The part axis reads it and
    # re-appends the "flesh" label; left in the name it blocked that label,
    # because _append_state declines a label whose token is still standing.
    r"meat\s+only|"
    r"organic|conventional|enriched|fortified|"
    r"ready[-\s]to[-\s]eat|"
    r"in\s+oil|in\s+water|in\s+juice|in\s+syrup|in\s+brine|"
    r"0%\s*moisture"
    r")\b",
    re.I,
)
# --- Preparation states that CHANGE composition -------------------------------
# _PREP_RE above treats every preparation token as nutritionally neutral and
# deletes it, so "Apples, dried" and "Apples, raw" both canonicalize to "apple".
# For form tokens (sliced / peeled / halved) that is correct. For tokens that
# remove water, remove fibre or add sugar or fat it is not: a canon is scored at
# its highest member for each nutrient, so one dehydrated
# variant hands the whole canon its per-100 g values. Measured on the shipped
# index: 114 canons carried fibre inflated >1.5x this way (median 2.9x, mushroom
# 16x from dried cloud ears, apple 3.9x from dehydrated rings), and 89 carried
# fat inflated >1.5x from a fried variant. Fibre is doubled by
# scoring.fiber_weight, so the error lands squarely on the substrate signal.
#
# These tokens are still stripped by _PREP_RE (so chunking is unchanged); they
# are additionally DETECTED here and re-attached as a canon suffix, which keeps
# "apple" and "apple, dried" as separate canons with separate nutrient vectors.
# Ordered: the first match wins, so "cranberry, dried, sweetened" reads as dried
# (water removal dominates the per-100 g change).
_STATE_PATTERNS = [
    ("dried", re.compile(
        r"\b(?:dried|dehydrated|desiccated|freeze[-\s]?dried|sun[-\s]?dried|sundried|"
        r"powdered|powder)\b|0%\s*moisture"
        # Bare "dry" is the dry-mix / dry-cereal / dry-legume sense ("Oat Bran,
        # dry", "Soybeans, dry, raw", "dry pasta, uncooked"), worth ~43 more
        # inflated canons. The lookahead drops the three collocations where it
        # is not dehydration at all: "dry roasted" (a roasting method), "dry
        # heat" (a cooking method) and "dry type" (vermouth, i.e. not sweet).
        r"|\bdry\b(?![-\s]+(?:roast|roasted|heat|type|sausage|cured))", re.I)),
    ("flour", re.compile(r"\bflour\b", re.I)),
    # "in juice" is a canning medium, not a juice product - the food is still
    # the fruit, so those phrases are masked out before this runs.
    # "juice pack" is the canning MEDIUM, not the form of the food: canned
    # pineapple packed in juice was read as pineapple juice.
    ("juice", re.compile(r"\b(?:juice(?!\s+pack)|nectar)\b", re.I)),
    ("paste", re.compile(r"\b(?:paste|concentrate)\b", re.I)),
    # "crisps" (plural noun) is the snack; bare "crisp" is an adjective and was
    # turning "Pear, crisp pear, Suli, ripe, raw" into 'pear, fried'.
    ("fried", re.compile(r"\b(?:fried|crisps|fritters?)\b", re.I)),
    # (?<!un) keeps "unsweetened" out. Bare "chips" is deliberately absent:
    # it would mislabel chocolate-chip products as fried.
    ("sweetened", re.compile(
        r"\b(?:jam|jelly|marmalade|candied|glac[ée]|in\s+(?:heavy\s+|light\s+)?syrup)\b"
        r"|(?<!un)\bsweetened\b", re.I)),
    # Sodium is the point. _PREP_RE strips salted/unsalted as a preparation
    # word, so "Almonds w skin salted" and "Almonds w skin unsalted" both
    # landed on 'almond', and the canon's sodium became the salted almond's -
    # the unsalted one cannot be scored below it. Two labels rather than one "salt" state, because a
    # single label cannot tell the two apart. Last in the list: a food that is
    # both salted and dried keeps the bigger compositional change.
]

# Salt is an INDEPENDENT axis, not one of the mutually-exclusive states above.
# Folding it into _STATE_PATTERNS was not enough: only the first matching state
# is kept, so "Cod, dried, salted" and "Cod, dried, unsalted" both resolve to
# the "dried" state and merge anyway. Detected separately and appended on top.
# The canning medium changes the food's composition and must survive: oil-packed
# and water-packed tuna differ by roughly an order of magnitude in fat, and both
# were resolving to the single canon "tuna tinned". Detected on its own axis for
# the same reason salt is - only one _STATE_PATTERNS label is ever kept.
# FDC writes the medium two ways - "canned in juice" and "canned, juice pack" -
# and only the first was read. The second left "juice pack" as a chunk of its
# own, so canned pineapple packed in juice canonicalised to 'pineapple juice',
# which is a different food entirely.
# Livsmedelsverket writes the medium with "w/" rather than "in" - "Green peas
# canned w/ brine", "Pineapple canned w/ juice" - which left 13 canons spelling
# the same medium as "X with brine" beside 105 spelling it "X, in brine". The
# "with" forms are only read inside _PACKED_CONTEXT_RE, so "porridge, made with
# water" and "macaroni, boiled with oil" are not read as packed.
_PACK_PATTERNS = [
    # the medium is often named - "canned in pear juice", "packed in olive oil" -
    # and the qualifier has to be inside the match or it is left behind as a
    # chunk: "Peach, canned in pear juice" canonicalised to 'peach juice'.
    ("in oil",   re.compile(r"\b(?:in|with)\s+(?:\w+\s+){0,2}oils?\b|\boils?\s+pack(?:ed)?\b", re.I)),
    ("in brine", re.compile(r"\b(?:in|with)\s+(?:\w+\s+){0,2}brine\b|\bbrine\s+pack(?:ed)?\b", re.I)),
    ("in juice", re.compile(r"\b(?:in|with)\s+(?:\w+\s+){0,2}juices?\b|\bjuices?\s+pack(?:ed)?\b", re.I)),
    ("in water", re.compile(r"\b(?:in|with)\s+(?:\w+\s+){0,2}water\b|\bwater\s+pack(?:ed)?\b", re.I)),
]
# "in syrup" is deliberately absent: _STATE_PATTERNS already reads it as the
# "sweetened" state, and listing it twice would print both labels.
#
# The medium only counts when the food is actually PACKED in it. Without this
# guard "Lobster, boiled/cooked in water" reads as water-packed lobster, and
# "Tuna, canned in water, drained solids" keeps a medium that has been poured
# away - the drained solids are what was measured.
_PACKED_CONTEXT_RE = re.compile(r"\b(?:canned|tinned|preserved|jarred|packed|bottled)\b", re.I)
_DRAINED_RE = re.compile(r"\bdrained\b", re.I)
# "prepared with water" is reconstitution, not a packing medium, and the two
# phrasings collide: FDC writes "Soup, vegetable chicken, canned, prepared with
# water", which has both a packed context and a "with water". The wording is
# blanked out of the copy the medium is read from.
_RECONSTITUTED_RE = re.compile(
    r"\b(?:prepared|made|mixed|diluted|reconstituted|heated|cooked|boiled|"
    r"simmered|thinned)\s+(?:up\s+)?with\b", re.I)

# Trim level, an axis of its own for the same reason salt is. "Veal" is a hard
# strip head, so every chunk after it was dropped and ONE canon absorbed 281
# descriptions - "separable fat" (~70% fat) sitting alongside "separable lean
# only" (~2%), and the canon reports the fat one. Beef and lamb carry the same defect at larger scale.
#
# "separable lean and fat" is the WHOLE cut and stays unlabelled, which is why
# the lean pattern refuses to match when "and fat" follows it.
_TRIM_PATTERNS = [
    # CNF writes "composite cuts, fat" where FDC writes "separable fat", so a
    # chunk that is nothing BUT the word fat counts too.
    ("separable fat",  re.compile(r"\bseparable\s+fat\b|\bfat\s+only\b|,\s*fat\s*(?=,|$)", re.I)),
    # a chunk that is nothing but "lean" says the same thing FDC's "separable
    # lean only" does, and it was taking a qualifier slot instead of a label
    ("lean",           re.compile(
        r"\bseparable\s+lean\b(?!\s+and\s+fat)|\blean\s+only\b|"
        r",\s*lean\s*(?=,|$)", re.I)),
]

# Wild or farmed is an axis, and it was being read by accident. The word
# usually survives as a chunk (613 canons carry ", wild"), but where the
# two-chunk strip head pushes it out it is lost silently, and the canon then
# holds both: 'salmon, atlantic' carried 8 wild members and 10 farmed ones,
# 'trout, rainbow' 4 and 7, 'catfish, channel' 4 and 4. Farmed Atlantic salmon
# runs to roughly twice the fat of wild and a very different omega-3 to omega-6
# ratio, so a canon holding both reports the farmed fat and describes neither fish.
#
# Only the FARMED wording is struck from the name (see the note where this is
# read). "Wild" is left standing, because it is half of a species name far more
# often than it is a provenance - wild rice, wild boar, wild mango, wild garlic,
# wild apricot - and _append_state already declines a label whose word is still
# in the canon, so those keep their names untouched and only the rows that
# actually lost the chunk get the label back.
_PROVENANCE_PATTERNS = [
    # "organic farming" is Frida's organic claim, not aquaculture; "cultured"
    # and "cultivated" are deliberately absent - the first is dairy
    # fermentation (cultured buttermilk, sour cream) and the second a growing
    # method for mushrooms and seaweed, not the fish axis.
    ("farmed", re.compile(
        r"\bfarmed\b|\bfarm[-\s]raised\b|\baquacultured?\b", re.I)),
    ("wild", re.compile(r"\bwild\b", re.I)),
]
# "Wild" only says something where the food is also FARMED or grown; where the
# species is only ever taken from the wild it is the unmarked case and the label
# just splits the canon. Measured on the corpus: 460 heads carry a "wild" row
# and only 16 carry a farmed one as well, the rest being capture fisheries and
# game - haddock, cod, mackerel, sardine, seal, moose, caribou, beaver,
# ptarmigan and some three hundred more fish. Labelling those took 16 of the 18
# 'seal, ringed' rows out of their own canon for nothing.
#
# Curated, and read off the whole description rather than the head, because the
# species name can sit anywhere in it. The list is the farmed and cultivated
# foods this corpus actually contains; a species that is only wild is absent by
# design, and so is rice - "wild rice" is a different GRAIN (Zizania, not
# Oryza), not a provenance, and it keeps the chunk it always had.
_FARMED_FOOD_RE = re.compile(
    r"\b(?:salmon|trout|char|carp|catfish|tilapia|perch|bass|seabass|"
    r"sea\s*bream|seabream|gilthead|turbot|halibut|cod|sturgeon|eel|"
    r"barramundi|kingfish|milkfish|pangasius|whitefish|"
    r"prawn|shrimp|crayfish|crab|lobster|oyster|mussel|clam|scallop|abalone|"
    r"rabbit|duck|goose|"
    r"blueberr(?:y|ies)|blackberr(?:y|ies)|raspberr(?:y|ies)|cranberr(?:y|ies)|"
    r"strawberr(?:y|ies)|gooseberr(?:y|ies)|cloudberr(?:y|ies))\b", re.I)
_FARMED_PHRASE_RE = re.compile(
    r"\s*,?\s*\b(?:farmed|farm[-\s]raised|aquacultured?)\b", re.I)


def _detect_provenance(d: str) -> tuple[str, str]:
    """Return (wild/farmed label, matched token), or ("",""). See above."""
    for label, rx in _PROVENANCE_PATTERNS:
        hit = rx.search(d)
        if hit:
            if label == "wild" and not _FARMED_FOOD_RE.search(d):
                return "", ""
            return label, hit.group(0)
    return "", ""


_SALT_PATTERNS = [
    # Order is load-bearing: first match wins, so every negative form has to be
    # listed before the positive one it contains. "no sodium added" must not be
    # read as "sodium added".
    ("unsalted", re.compile(
        r"\b(?:un|non)-?salted\b|\bwithout\s+salt\b|\bno\s+salt\s+added\b|"
        r"\bno\s+sodium\s+added\b|\bsodium[-\s]?free\b|"
        # CIQUAL and Fineli put the "added" in front - "Peanut, no added salt" -
        # and only the trailing form was listed, so the claim was not read at
        # all and 37 canons carried the residue "no" as a chunk of their name.
        r"\bno\s+added\s+(?:salt|sodium)\b|"
        # Fineli coordinates the two claims - "Fried Without Fat And Salt" -
        # and the adjacent form below cannot reach the second half
        r"\bwithout\s+(?:\w+\s+){0,2}and\s+salt\b|"
        r"\bsalt[-\s]?free\b", re.I)),
    ("half-salted", re.compile(r"\bhalf[-\s]?salted\b|\bdemi[-\s]?sel\b", re.I)),
    # A reduced-sodium canned bean is not an unsalted one and not a fully salted
    # one either. Without a label of its own the range collapses onto its
    # salted end:
    # 'bean, great northern' held "canned, sodium added", "canned, no salt added"
    # and "canned, reduced sodium" as one food, and 'beef' and 'pork' did the
    # same across 1,292 and 894 members.
    ("low-salt", re.compile(
        # "less salt" is Fineli's wording and was not read at all, so the
        # reduced-salt rye bread kept the phrase as a chunk of its name
        r"\b(?:reduced|low|lower|less)[-\s](?:sodium|salt)"
        r"(?:\s*(?:and|&|/)\s*sugars?)?\b|"
        r"\blightly\s+salted\b", re.I)),
    # FDC writes the positive form as "sodium added" rather than "salted", which
    # is why 'blackeye pea, sodium added' ended up carrying the phrase in its
    # name instead of resolving to the salt axis.
    ("salted", re.compile(
        r"(?<!un)(?<!half.)(?<!half)\bsalted\b|\bwith\s+salt\b|"
        # ...and the front-loaded spelling, which nine rows use and none of the
        # three trailing forms below reach: "Pasta, cooked, with added salt" and
        # "Oat Macaroni, Boiled With Added Salt" were reading as unmarked and
        # merging into the unsalted canon. Safe here because the unsalted
        # pattern is tried first and takes "no added salt" out of the running.
        r"\b(?:with\s+)?added\s+(?:salt|sodium)\b|"
        r"\bsalt\s+added\b|\bsodium\s+added\b|\badded\s+sodium\b", re.I)),
]

# Fat level is an axis too, and the least marked case is the full-fat one. None
# of skim / nonfat / low-fat is a preparation word, so they were surviving in
# some names and being dropped from others, and 43 canons ended up holding both
# ends: yoghurt (281 members), cheese, cheddar (163) and cheese, mozzarella
# (162) each pulled the skimmed version's numbers toward the full-fat ones.
# "light" is deliberately absent - light syrup, light beer and light cream are
# three different claims and only one of them is about fat.
_FAT_PATTERNS = [
    # The noun the marker qualifies has to be INSIDE the match, or removing the
    # wording leaves it standing on its own: "part skim milk" left 'milk' as a
    # chunk ('cheese, mozzarella milk') and "milk-fat free" left 'milk'.
    ("fat-free", re.compile(
        r"(?<!part[-\s])(?<!semi[-\s])(?<!half[-\s])(?<!partly[-\s])"
        r"(?<!partially[-\s])\bskim(?:med)?(?:\s+milk)?\b|"
        r"\b(?:milk[-\s]?)?non[-\s]?fat(?:\s+milk)?\b|"
        r"\b(?:milk[-\s]?)?fat[-\s]free(?:\s+milk)?\b|\b0\s*%\s*fat\b|"
        r"\bdefatted\b", re.I)),
    ("low-fat", re.compile(
        r"\blow[-\s]?fat(?:\s+milk)?\b|\breduced[-\s]?fat(?:\s+milk)?\b|"
        r"\bsemi[-\s]?skim(?:med)?(?:\s+milk)?\b|"
        r"\bpart(?:ly|ially)?[-\s]?skim(?:med)?(?:\s+milk)?\b|"
        r"\bhalf[-\s]?skim(?:med)?(?:\s+milk)?\b|\bhalf[-\s]?fat\b", re.I)),
]

# On a dairy food "whole" written as a CHUNK is the fat level, and it was the
# third spelling of one claim: 'greek yoghurt' held the unmarked rows, 'greek
# yoghurt, whole' the ones spelling it out and 'greek yoghurt, full-fat' the
# ones using the house word. Seven foods carried both "X" and "X, whole" - cow,
# goat, sheep and buffalo milk, condensed milk, milk and greek yoghurt.
#
# Only as a chunk, and only on a dairy head. Adding "whole milk" and "full-fat"
# to _FAT_PATTERNS above looks equivalent and is not: the phrase is usually an
# INGREDIENT ("sweet yeast dough for pie, whole milk", "multigrain bread ...
# full-fat"), which labelled breads and doughs by the milk in them, and striking
# the wording took the food's own noun with it - "condensed whole milk,
# sweetened" came out as 'condensed, sweetened, full-fat'.
# Everywhere else "whole" is a FORM - a whole almond, a whole carrot, a whole
# chicken - and per 100 g that is the same food as the sliced one.
_DAIRY_HEAD_RE = re.compile(
    r"^\s*(?:[\w'\u2019-]+\s+){0,2}"
    r"(?:milks?|creams?|yogh?o?urts?|quarks?|kefir|buttermilks?|curds?|skyr)\b",
    re.I)
# the source tag has not been stripped from `d` yet where this fires, so the
# chunk can be the last one and still be followed by " [CIQUAL]"
_WHOLE_CHUNK_RE = re.compile(r"(?:^|,)\s*whole\s*(?=\s*(?:,|\[|$))", re.I)

# Enrichment is stripped as a preparation word, so the distinction was always
# lost: 52 canons hold both an enriched and an unenriched member, among them
# rice (634 members), wheat flour (98) and pasta (95). Enriched white rice
# carries roughly five times the iron and thirty times the folate of the
# unenriched grain, which for a nutrient resource is the whole point.
_FORTIFY_PATTERNS = [
    # The negative forms have to be excluded explicitly: "salt, not fortified
    # with iodine" and "dextrose tablets, non-fortified" both contain the word.
    # The nutrient name has to be INSIDE the match. Without it "Milk, skimmed,
    # with added vitamin D" lost the phrase but kept the "D", which then welded
    # itself to the food as 'milk d'.
    ("enriched", re.compile(
        r"(?<!un)(?<!non[-\s])(?<!not\s)\b(?:enriched|fortified)\b|"
        # "ascorbic acid" is how FDC and CNF spell added vitamin C on a juice -
        # "Apple juice, canned or bottled, unsweetened, with added ascorbic acid" -
        # and it was absent from every alternation below, so 17 rows read as
        # unfortified. 'apple juice' held 9 of them beside 23 plain members, and
        # under the max the plain juice reported the fortified one's vitamin C.
        r"\bwith\s+added\s+(?:vitamins?(?:\s+[a-k]\d*\b)?|minerals?|calcium|iron|"
        r"folate|folic\s+acid|iodine|zinc|fibre|fiber|ascorbic\s+acid)"
        r"(?:\s*(?:,|and/or|and|or|&)\s*(?!without\b)(?:added\s+)?"
        r"(?:vitamins?(?:\s+[a-k]\d*\b)?|minerals?|"
        r"calcium|iron|folate|folic\s+acid|iodine|zinc|fibre|fiber|"
        r"ascorbic\s+acid))*|"
        # CIQUAL and NEVO write the claim WITHOUT the word "added" - "Cereal bar
        # with fruit with vitamins and minerals", "fortified w Ca and Vit B12" -
        # so the label printed but the wording stayed, and 120 canons carried the
        # nutrient list beside their own "enriched". "fibre" is deliberately not
        # in this half: fibre is a property of the food as often as an addition.
        # "calcium sulphate" is tofu's coagulant, not a fortificant, and eating
        # the "calcium" half of it left 'tofu, sulfate' calling itself enriched
        # Fineli and FDC drop the "with" as well as the "added" - "Orange Juice,
        # Unsweetened, Added Calcium", "Mineral Water, Added Zinc And Vitamin E",
        # "reduced calorie, with added vitamin E" - and the phrase was reached by
        # neither half above, so 68 rows read as unfortified. The negatives have
        # to be shut out by hand here: this branch has no "with" to hang a
        # lookbehind on, so "no added iron" and "without added vitamin A" would
        # otherwise match on their tail. _UNFORTIFIED_PHRASE_RE still owns the
        # wording of those, and reads them before this can.
        r"(?<!no\s)(?<!not\s)(?<!without\s)(?<!and\s)\badded\s+"
        r"(?:vitamins?(?:\s+[a-k]\d*\b)?|minerals?|calcium|iron|folate|"
        r"folic\s+acid|iodine|zinc|ascorbic\s+acid|thiamin[e]?|niacin|"
        r"riboflavin|selenium|fibre|fiber)"
        r"(?:\s*(?:,|and/or|and|or|&)\s*(?!without\b)(?:added\s+)?"
        r"(?:vitamins?e?(?:\s+[a-k]\d*\b)?|minerals?|calcium|iron|folate|"
        r"folic\s+acid|iodine|zinc|ascorbic\s+acid|thiamin[e]?|niacin|"
        r"riboflavin|selenium|fibre|fiber|[a-k]\d{0,2}))*|"
        r"\bwith\s+(?:vitamins?e?(?:\s+[a-k]\d*)?|vit\.?\s*[a-k]\d*|minerals?|"
        r"(?:calcium|ca)(?!\s*(?:sulph|sulf|chlor|carbon|lact|phosph|citr|salt))|"
        r"iron|folate|folic\s+acid|iodine|zinc|ascorbic\s+acid)\b"
        # the continuation has to accept a BARE letter, or "with vitamin E and
        # C" consumed only its first half and left " and c" as the food's name
        r"(?:\s*(?:,|and/or|and|or|&)\s*(?!without\b)"
        r"(?:vitamins?e?(?:\s+[a-k]\d*)?|vit\.?\s*[a-k]\d*|minerals?|"
        r"calcium|ca|iron|folate|folic\s+acid|iodine|zinc|fibre|fiber|"
        r"[a-k]\d{0,2})\b)*"
        # ...and a trailing "added", or "with vit C added" left "added" behind
        r"(?:\s+added)?", re.I)),
]

# FDC deposits a family of analytical rows whose values are reported on a DRY
# MATTER basis - "Beans, Dry, Black (0% moisture)" - rather than per 100 g of
# food as eaten, which is what every other row in the export means and what the
# table header promises. 17 canons hold them and in several they are the
# majority: 253 of the 261 members of 'bean, pinto, dried'. They are real
# measurements, so they are kept and labelled rather than dropped; the label is
# what stops a dry-matter figure standing in for the as-eaten values of the same bean.
# Six claims below survive nowhere else, because every food that carries them
# is a HARD cultivar strip head - lima bean, cowpea, broad bean, banana, rice
# and wheat flour all drop every chunk after the head - so a label is the only
# form in which they can reach the canon at all. Measured on the round-16
# export; each count is canons holding BOTH ends of the claim as one food.
#
# FIBRE, as a figure and as a claim. The figure was the only nutrient
# percentage _QUANT_RE ate: "Cheese, 13% fat" keeps its number because
# _detect_fat reads it into a label and _append_states puts it back, while
# "Crispbread, Rye, 17% Fibre" came out bare. Read the same way here, for the
# same reason. "high fibre" is the wording without a number and goes with it.
_FIBRE_PCT_RE = re.compile(r"\d+(?:[.,]\d+)?\s*%\s*fib(?:re|er)\b", re.I)
_HIGH_FIBRE_RE = re.compile(
    r"\bhigh[-\s]?fib(?:re|er)\b|\bfib(?:re|er)[-\s]enriched\b|"
    r"\bfib(?:re|er)[-\s]rich\b", re.I)
# "app." is Swedish FoodComp's approximation word - "Bread white unsweetened
# app. 2.5% fibre" - and it sits in front of the figure exactly as "around" and
# "approx." do in front of a fat one, where _FAT_PHRASE_RE already removes it.
_FIBRE_PHRASE_RE = re.compile(
    r"\s*,?\s*(?:\b(?:app|approx|about|around|env)\.?\s+)?"
    + _FIBRE_PCT_RE.pattern + r"|\s*,?\s*(?:" + _HIGH_FIBRE_RE.pattern + r")", re.I)

# HARVEST MATURITY. FDC files the green legume and the dried one as "immature
# seeds" and "mature seeds" of the same food, and both chunks die on the strip
# head: 'lima bean' holds 8 immature-seed rows beside 7 mature-seed ones,
# 'cowpea' 8 beside 7, 'broad bean' 4 beside 4, 'sweet corn' 5 immature-kernel
# rows. Raw immature lima beans carry about 4.9 g fibre and 6.8 g protein per
# 100 g against 19 g and 21.5 g for the mature seed - a four-fold gap inside one
# canon, and the score reports the larger of the two. Mature is left unmarked:
# it is the unqualified bean, and _MATURITY_SKIP_RE already keeps "mature seeds"
# out of the cheese-ageing axis for the same reason.
_IMMATURE_RE = re.compile(
    r"\bimmature\s+(?:seeds?|kernels?|beans?|pods?|grains?|fruits?)\b|"
    r"(?:^|,)\s*immature\s*(?=\s*(?:,|\[|$))|\bimmature\s+green\b", re.I)

# RIPENESS, which is starch against sugar and therefore the substrate question
# this resource exists to answer. FDC measured the stages separately and named
# them: 'banana' holds 224 rows saying ripe, 164 saying slightly ripe, 164
# saying overripe and 15 saying green or unripe, as ONE 594-member canon, so
# under the max a banana reports the green fruit's starch and the overripe
# fruit's sugars at the same time. Only the off-default ends are labelled - a
# plain "Bananas, raw" IS the ripe fruit, so printing "ripe" would split the
# canon against itself, the same argument that makes unsalted and unsweetened
# silent. No context gate: unlike "green", these words are never a colour, a
# cultivar or a cut anywhere in the corpus.
_RIPENESS_PATTERNS = [
    ("overripe", re.compile(r"\bover[-\s]?ripe\b", re.I)),
    ("slightly ripe", re.compile(
        r"\bslightly\s+ripe\b|\bpartially\s+ripe\b|\bhalf[-\s]?ripe\b", re.I)),
    ("unripe", re.compile(r"\bunripe\b|\bgreen[-\s]?mature\b", re.I)),
]

# DECAFFEINATED. Caffeine and theobromine are both nutrients in this export, so
# a decaffeinated coffee folded into the regular one is handed the regular one's
# caffeine outright - 19 canons, 'coffee, instant, dried' holding 6 decaf rows
# beside 11 regular. The word is unambiguous wherever it appears.
_DECAF_RE = re.compile(
    r"\bdecaff?eina(?:ted|ed)?\b|\bdecaf\b|\bcaffeine[-\s]free\b", re.I)

# GRAIN REFINEMENT, which is the fibre and mineral axis of every cereal and had
# no axis at all. 'wheat flour' holds 21 whole-grain rows beside 5 refined ones
# (Frida's "extraction rate 75%" among them) and 'rice' 12 beside 11 of
# Phenol-Explorer's "Rice, refined"; WAFCT's three maize porridges each hold the
# wholegrain meal beside the degermed grit, and degerming is what removes the
# thiamin and the magnesium. Read only on a CEREAL and never on a fat, so
# _OIL_GRADE_PATTERNS keeps "refined" for the pressing grade of an oil, which is
# a different claim about a different food. The wholegrain end stays unmarked:
# the sources that bother to say it usually say it in the head, where
# _append_state already declines to print a label whose word is standing.
_CEREAL_CONTEXT_RE = re.compile(
    r"\b(?:flours?|breads?|rice|wheat|pastas?|noodles?|spaghetti|macaroni|"
    r"oats?|barley|rye|maize|corn|millet|sorghum|teff|spelt|semolina|couscous|"
    r"bulgur|crackers?|crispbread|biscuits?|cereals?|porridge|tortillas?|"
    r"chapati|grits?|meal|bran|groats?)\b", re.I)
_REFINED_GRAIN_RE = re.compile(
    r"(?<!un)\brefined\b|\bdegermed\b|\bsifted\b|\bpolished\b|"
    r"\bextraction\s+rate\b|\b\d+\s*%\s*extraction\b", re.I)
# BOTH ends are labelled here, unlike salt and sugar where one end is silent,
# because neither is the unmarked case: most rows state no refinement at all.
# Labelling only "refined" would have left 'wheat flour' holding its 21
# whole-grain rows beside 73 unstated white ones - the 3x fibre gap the axis
# exists to separate - so the whole-grain end has to be named too. What is left
# unlabelled is then honestly what the sources did not say.
_WHOLEGRAIN_RE = re.compile(
    r"\bwholemeal\b|\bwholegrain\b|\bwhole[-\s]grain\b|\bwhole[-\s]wheat\b|"
    r"\bwholewheat\b|\bgraham\s+flour\b", re.I)

# GLUTEN-FREE, which on a pasta or a bread is not a label on the same food but a
# different grain: gluten-free pasta is rice or maize, and 'pasta' held 9 of
# them beside 74 durum-wheat rows, 'pasta, dried' 6 beside 36. 33 canons.
_GLUTEN_FREE_RE = re.compile(
    r"\b(?:naturally\s+)?gluten[-\s]free\b", re.I)

# The CUT of a grain is not its composition: rolled, steel-cut, quick and
# old-fashioned oats are the same groat per 100 g. They were invisible until
# the wholegrain wording above was lifted off the chunk in front of them, and
# then they split one honest 'oat, wholegrain' into three - one of them named
# 'oat, steel', because _PREP_RE eats the "cut" and leaves the adjective.
# Deliberately short: "cracked" and "flaked" are NOT here, because cracked
# wheat is bulgur and flaked rice is poha, and both are their own food.
_GRAIN_FORM_RE = re.compile(
    r"^(?:steel(?:[-\s]?cut)?|rolled|old[-\s]?fashioned|"
    r"quick(?:[-\s]?cooking)?|pinhead|jumbo)$", re.I)

_MOISTURE_PATTERNS = [
    ("0% moisture basis", re.compile(
        r"\b0\s*%\s*moisture\b|\bmoisture[-\s]free\s+basis\b|"
        r"\bdry\s+matter\s+basis\b|\bmoisture\s+basis\b", re.I)),
]

# The peel measured on its own is not the fruit. BioFoodComp separates them and
# files each as its own row, so 70 of the 309 members of 'apple' were apple PEEL
# - several times the fibre and polyphenol content of the flesh - folded
# into whole apple, whose fibre then read as the peel's. 38 canons carry the same defect, including nectarine (20 of
# 175), peach (20 of 181) and, for skin rather than peel, 'chicken fat' (15 of
# 18) and 'chicken, broiler' (8 of 116), where the skin is about 30 g fat/100 g.
#
# The guard is the shape rather than a lookahead: the part must be a WHOLE comma
# chunk. "Pear, Bosc, peel and flesh, green, raw" is the entire fruit and says so
# in one chunk, so it does not match and is left alone.
# Peel and skin were the first two; the rest are the same defect measured again
# on the round-9 index. 247 rows named FLESH (the peeled food) in a canon that
# did not, 60 named BRAN - 54 of them in 'rice', where rice bran carries about
# 20 g of fat against milled rice's 0.7 - 48 named a SPROUT and 38 a LEAF.
#
# "root" and "tuber" are deliberately absent: for a carrot the root IS the food
# and for a potato the tuber is, so labelling them would split the plain rows
# from the ones whose source spelled the default out. _PART_PRESERVE already
# carries the head/part pairs where a root or tuber is the odd organ.
_PART_PATTERNS = [
    ("peel", re.compile(r"(?:^|,)\s*(?:peel|rind|zest)s?\s*(?:only)?\s*(?=,|$)", re.I)),
    ("skin", re.compile(r"(?:^|,)\s*skins?\s*(?:only)?\s*(?=,|$)", re.I)),
    # "without skin" says the same thing "flesh" does - 98 rows spelled it that
    # way inside a canon that also held the with-skin rows, and the difference
    # is real: an apple's skin carries much of its fibre.
    # FDC spells the same statement "meat only" on poultry where Fineli and
    # NEVO write "without skin", and the two forms sat in different canons:
    # a skinless chicken breast carries about a sixth of the fat of one with
    # the skin on, and the canon reported the skin-on figure for both.
    ("flesh", re.compile(
        r"(?:^|,)\s*(?:flesh|pulp)\s*(?:only)?\s*(?=,|$)"
        # ...and this one need not be a whole chunk. WAFCT welds it to the cut -
        # "Chicken, light meat without skin, raw" - and the whole-chunk test put
        # four skinless rows in the with-skin canon. "Without skin" can only
        # ever mean skinless, wherever it sits.
        # "without skin AND stone" is how a peach is eaten and _PREP_RE already
        # owns that phrase; only the bare claim is read here.
        r"|\bwithout\s+skins?\b(?!\s+and\s+(?:stone|pit|core|seed))"
        # the lookahead has to allow a trailing parenthetical: FDC writes
        # "Bear, black, meat (Alaska Native)" and the phrase is the last thing
        # before it
        r"|(?:^|,)\s*meat\s+only\s*(?=\s*(?:,|\(|\[|$))|\bskinless\b", re.I)),
    ("bran", re.compile(r"(?:^|,)\s*brans?\s*(?=,|$)", re.I)),
    ("germ", re.compile(r"(?:^|,)\s*germs?\s*(?=,|$)", re.I)),
    ("husk", re.compile(r"(?:^|,)\s*(?:husks?|hulls?)\s*(?=,|$)", re.I)),
    ("sprouted", re.compile(r"(?:^|,)\s*sprout(?:ed|s)?\s*(?=,|$)", re.I)),
    ("leaf", re.compile(r"(?:^|,)\s*(?:leaf|leaves)\s*(?=,|$)", re.I)),
    # (the label a match produces is the matched word itself, singularised, so
    #  "stalks" stays 'stalk' rather than being renamed to the pattern's key)
    ("stem", re.compile(r"(?:^|,)\s*(?:stems?|stalks?)\s*(?=,|$)", re.I)),
    ("pod", re.compile(r"(?:^|,)\s*pods?\s*(?=,|$)", re.I)),
    ("shell", re.compile(r"(?:^|,)\s*shells?\s*(?=,|$)", re.I)),
]

# Sugar is an independent axis for the same reason salt is. "Pears, average,
# stewed with sugar" and "Pears, average, stewed without sugar" were one canon
# with 222 members, and because a canon is scored at its highest member the
# plain stewed pear was being handed the sweetened one's sugar. 21 canons held both
# forms, among them apple (312 members), yoghurt (288) and fig, dried (88).
# Separate from _STATE_PATTERNS, which is mutually exclusive: a food can be
# dried AND sweetened, and only one state slot exists there.
# Pressing grade is composition, not marketing. Extra virgin olive oil carries
# an order of magnitude more polyphenols than refined oil pressed from the same
# fruit, and Phenol-Explorer files the three grades as separate entries for
# exactly that reason - but FDC and CIQUAL write the grade as a trailing chunk,
# where the two-chunk rule dropped it and refined oil then inherited the extra
# virgin polyphenol maxima. Read only in a fat context, because "refined"
# also grades flour, sugar and salt, where it means something else.
_OIL_CONTEXT_RE = re.compile(
    r"\b(?:oils?|fats?|ghee|shortening|margarine|tallow|lard|dripping)\b", re.I)
_OIL_GRADE_PATTERNS = [
    ("extra virgin", re.compile(r"\bextra[-\s]?virgin\b", re.I)),
    ("virgin", re.compile(r"(?<!extra\s)(?<!extra-)\bvirgin\b", re.I)),
    ("refined", re.compile(r"(?<!un)\brefined\b", re.I)),
    ("cold-pressed", re.compile(r"\bcold[-\s]?press(?:ed)?\b", re.I)),
]


# Colour is composition wherever a pigment is the nutrient. Black rice carries
# anthocyanins that white rice has none of; the same holds for black against
# white beans (Phenol-Explorer files them separately for exactly that reason),
# red against white grapefruit juice (lycopene), red against yellow tomato, and
# green against black olive. 110 canons held BOTH a dark- and a light-coloured
# member, 'rice' (636 members) and 'common bean' (223) among them, and under
# the max each handed the pale food the pigmented one's values outright.
#
# Only a WHOLE chunk counts, which is what keeps a cultivar name out: "Apples,
# raw, red delicious" and "Potato tuber, Red LaSoda" name a variety, not a
# colour, and stripping to 'apple, red' would claim something the row does not.
# Phenol-Explorer writes the colour as a bracket variant instead - "Common bean
# [Black]" - so that spelling is read too.
#
# Bare "light" and "dark" are NOT colours here: on a dairy or a spread they mean
# reduced fat, which the fat axis already owns. They are read only as
# "<dark|light> meat", where FDC and CNF use them for the poultry cut - and
# there the difference is real, chicken dark meat carrying about twice the fat
# of light.
_COLOUR_WORDS = (r"black|red|white|green|yellow|purple|pink|brown|golden|blue|"
                 r"ivory|cream[-\s]?coloured")
_COLOUR_PATTERNS = [
    ("dark meat", re.compile(r"\bdark\s+meat\b", re.I)),
    ("light meat", re.compile(r"\blight\s+meat\b", re.I)),
    ("colour", re.compile(
        rf"(?:^|,)\s*(?:{_COLOUR_WORDS})"
        # "or" and "/" join two names for ONE colour ("pink or red"); "and"
        # joins two different ones ("tea, black and green" is a blend), and a
        # blend is not either of them
        rf"(?:\s*(?:or|/)\s*(?:{_COLOUR_WORDS}))?"
        rf"(?:\s+(?:flesh|skinned|fleshed))?\s*(?=,|$)"
        rf"|\[\s*(?:{_COLOUR_WORDS})\s*\]", re.I)),
]


# "Egg, white" is the albumen, not a white egg, and "Eggs, yolk" its pair. The
# colour axis took the chunk and left 'egg, dried' holding dried egg white.
_COLOUR_HEAD_SKIP_RE = re.compile(r"^\s*eggs?\s*,", re.I)
# On a refined cereal or a sugar, "white" is not a variety - it is the refining
# state, and it is the unmarked one. Printing it would split the plain rows from
# the ones whose source happened to spell the default out, so the label is
# suppressed while the DETECTION stands: that is what still keeps white rice
# apart from the black, brown and red rices that carry the pigment. Everywhere
# else white IS a variety - white beans, white cabbage, white grapefruit - and
# the label is printed.
_WHITE_IS_DEFAULT_RE = re.compile(
    r"^\s*(?:rice|flour|sugars?|bread|maize|corn|corn\s*meal|wheat|barley|"
    r"sorghum|millet|teff|fonio|pasta|noodles?|semolina|quinoa|couscous|"
    r"polenta|grits)\b", re.I)
_WHITE_RE = re.compile(r"^(?:white|ivory|cream[-\s]?coloured)$", re.I)


def _detect_colour(d: str) -> tuple[str, str]:
    """Return (colour label, matched token), or ("",""). See _COLOUR_PATTERNS."""
    if _COLOUR_HEAD_SKIP_RE.match(d):
        return "", ""
    for label, rx in _COLOUR_PATTERNS:
        hit = rx.search(d)
        if hit:
            if label != "colour":
                return label, hit.group(0)
            word = hit.group(0).strip(" ,[]").lower()
            word = re.sub(r"\s+(?:flesh|skinned|fleshed)$", "", word).strip()
            # CIQUAL writes a disjunction for one colour - "pink or red" - and
            # the second name is the common one, so it is what the label uses.
            word = re.split(r"\s*(?:or|/)\s*", word)[-1]
            if _WHITE_RE.match(word) and _WHITE_IS_DEFAULT_RE.match(d):
                # detected, so the split against black and brown holds, but not
                # printed: on a refined cereal white is the unmarked form
                return "", hit.group(0)
            return word, hit.group(0)
    return "", ""


_SUGAR_PATTERNS = [
    # Negative first, as above: "no sugar added" contains "sugar added".
    ("unsweetened", re.compile(
        r"\b(?:un|non)-?sweetened\b|\bno\s+(?:added\s+)?sugars?\s+added\b|"
        r"\bno\s+added\s+sugars?\b|\bwithout\s+(?:added\s+)?sugars?\b|"
        r"\bsugar[-\s]?free\b|\bwithout\s+sweetener\b", re.I)),
    # A reduced-sugar claim is a level on this axis, not a chunk of the name:
    # "Milk, chocolate, lowfat, reduced sugar" kept it as one.
    ("low-sugar", re.compile(
        r"\b(?:reduced|low|lower|less)\s+sugars?\b|\blow[-\s]sugar\b|"
        r"\breduced[-\s]sugar\b", re.I)),
    # Light syrup carries roughly half the added sugar of heavy syrup, and the
    # two were separate canons before the wording was stripped; folding them
    # into one "sweetened" would hand light-syrup fruit the heavy-syrup sugar.
    ("lightly sweetened", re.compile(r"\blight\s+syrup\b", re.I)),
    ("sweetened", re.compile(
        r"(?<!un)(?<!non)\bsweetened\b|\bwith\s+(?:added\s+)?sugars?\b|"
        # Fineli coordinates the two claims - "With Added Salt And Sugar" - and
        # the sugar half sits behind the conjunction where no form below reaches
        # it; the canon came out 'peanut butter with and sugar'.
        r"\badded\s+(?:salt|sodium)\s+and\s+sugars?\b|"
        r"\bsugars?\s+added\b|\badded\s+sugars?\b|\bwith\s+sweetener\b|"
        r"\b(?:extra\s+)?(?:heavy|medium)\s+syrup\b|\bsyrup\s+pack\b", re.I)),
]
# Plant parts that must survive a strip head. The strip-head rule drops
# everything after the first comma, which is right for a cultivar and wrong for
# a different organ of the same plant: cassava leaves are a leafy vegetable and
# cassava tuber is a starch (leaves carry 5.6x the tuber's protein), tamarind
# leaves carry 28x the fibre of tamarind pulp. Curated rather than a blanket
# part list on purpose - for lettuce, spinach and kale the leaf IS the food, so
# protecting "lettuce, leaf" would fragment 354 records for nothing. Only
# head/part pairs where the part is a DIFFERENT organ from the food the head
# names belong here.
# The CUT is composition too, and it was going the same way the organ went.
# "Beef" is a cultivar strip-head - it has to be, or every breed name ("Japanese
# beef cattle", "dairy fattened steer", "Belgian Blue") becomes a canon - and
# the head drops everything behind it, so chuck, brisket, tenderloin and eye of
# round all landed in one 1,291-member 'beef'. Beef chuck carries roughly four
# times the fat of eye of round, and a single 'beef' canon reports the chuck
# figure for all of them.
#
# Read only in a MEAT context and only as a whole chunk. Outside one, "round",
# "breast", "leg" and "plate" are ordinary words, and "rib" is a vegetable part.
_MEAT_CONTEXT_RE = re.compile(
    r"^\s*(?:beef|pork|lamb|veal|mutton|goat|venison|bison|buffalo|horse|"
    r"rabbit|hare|chicken|turkey|duck|goose|quail|pheasant|partridge|guinea|"
    r"ostrich|emu|game\s+meat|deer|elk|moose|caribou|reindeer|antelope|boar|"
    r"kangaroo|camel|alpaca|llama|yak|seal|whale|walrus)\b", re.I)
# Longest first: the alternation is first-match-wins, so "spare ribs" has to be
# offered before "ribs" or a spare rib is filed as a rib.
_CUT_WORDS = (
    r"short\s*ribs?|spare\s*ribs?|rib\s*eye|ribeye|chuck(?:\s+eye)?|ribs?|"
    r"short\s+loin|tenderloin|striploin|longissimus(?:\s+dorsi)?|loin|"
    r"top\s+sirloin|sirloin|"
    r"eye\s+of\s+round|top\s+round|bottom\s+round|tip\s+round|round|"
    r"brisket|flank|plate|fore\s*shank|hind\s*shank|shank|skirt|"
    r"flat\s+iron|hanger|porterhouse|t[-\s]?bone|chump|"
    r"blade|neck|rump|topside|silverside|thick\s+flank|knuckle|"
    r"leg|shoulder|breast|belly|hock|hip|picnic|ham|wing|thigh|drumstick|"
    r"back|saddle|rack|escalope|cutlet|medallion")
# FDC names the retail cut, not the primal: "Beef, shoulder top blade steak,
# boneless, ...". The primal word is inside the chunk rather than the whole of
# it, so a whole-chunk test missed 182 of the 579 members of 'beef'. Up to three
# words may sit in front of it and a retail noun behind it; the LABEL is the
# primal, which is the granularity that changes composition.
_CUT_RETAIL = r"steaks?|roasts?|chops?|joints?|cutlets?|fillets?|medallions?|muscles?"
_CUT_PATTERNS = [
    # Mince is cut from trim, so it carries more fat than the joints it is made
    # from, and it presents orders of magnitude more surface to a gut bacterium.
    # Three sources state it and the corpus handled it three different ways:
    # "Turkey, ground" became 'ground turkey', "Lamb, ground" kept the chunk,
    # and "Beef, ground" / "Pork, ground" - both gated strip heads - lost it
    # outright and sat in the bare canon beside a strip steak. Read FIRST, ahead
    # of the primal cuts: on a mince the grind is the larger statement. The
    # animal gate in _detect_cut is what keeps ground SPICES out of it.
    ("ground", re.compile(r",\s*(grounds?|minced?)(?:\s+meat)?\s*(?=,|$)", re.I)),
    # same guard as the organ axis: the cut can only sit BEHIND the animal
    ("cut", re.compile(
        # the prefix is LAZY: greedy, it swallowed "spare " and then matched the
        # bare "ribs", filing a spare rib as a rib
        rf",\s*(?:[a-z]+\s+){{0,3}}?({_CUT_WORDS})(?:\s+(?:{_CUT_RETAIL}))?\s*(?=,|$)",
        re.I)),
]


def _detect_cut(d: str) -> tuple[str, str]:
    """Return (primal-cut label, matched token) for a meat, or ("","")."""
    if not _MEAT_CONTEXT_RE.match(d):
        return "", ""
    for _, rx in _CUT_PATTERNS:
        hit = rx.search(d)
        if hit:
            _cutw = _fold_plural(re.sub(r"\s+", " ", hit.group(1).strip().lower()))
            # one label for the two spellings: FDC grinds and McCance minces
            if _cutw in ("ground", "minced", "mince"):
                _cutw = "ground"
            # the longissimus dorsi IS the loin muscle; FDC and CNF name it
            # anatomically on 115 pork rows and by the cut everywhere else
            if _cutw.startswith("longissimus"):
                _cutw = "loin"
            return _cutw, hit.group(0)
    return "", ""


# Maturity is composition on a cheese: an aged cheddar has lost water, so it
# carries more fat, protein and sodium per 100 g than a mild one. 116 of the 400
# members of 'cheese, cheddar' said which they were and the canon did not.
# "mature seeds" is NOT this axis - it is FDC's phrase for the dry legume seed,
# and _PREP_RE already owns it.
# Matched as a WORD rather than a whole chunk, because the grade is rarely one:
# FDC writes "cheddar, natural shredded sharp" and puts it inside the brand on
# "Cheddar cheese, sliced, SARGENTO SHARP". 124 of the 284 members of
# 'cheese, cheddar' said sharp and the canon did not. The cheese gate is what
# makes a bare word safe here.
_MATURITY_PATTERNS = [
    ("mature", re.compile(
        r"\b(?:extra\s+)?(?:sharp|matured|aged|vintage|ripened)\b"
        r"|(?:^|,)\s*(?:extra\s+)?mature\s*(?=,|$)", re.I)),
    ("mild", re.compile(r"\b(?:mild|young|baby)\b", re.I)),
]
_MATURITY_SKIP_RE = re.compile(r"\bmature\s+seeds?\b", re.I)
# Read only on a CHEESE. Everywhere else "mature" is a ripeness stage rather
# than an ageing one - BioFoodComp writes it on beans, peaches and vetch - and
# _PREP_RE already owns that word, so reading it here split 150 rows off their
# canons for nothing.
_MATURITY_CONTEXT_RE = re.compile(r"\bchees[e]?\b|\bcheddar\b|\bgouda\b", re.I)


def _detect_maturity(d: str) -> tuple[str, str]:
    """Return (maturity label, matched token), or ("",""). See _MATURITY_PATTERNS."""
    if _MATURITY_SKIP_RE.search(d) or not _MATURITY_CONTEXT_RE.search(d):
        return "", ""
    for label, rx in _MATURITY_PATTERNS:
        hit = rx.search(d)
        if hit:
            # the young end keeps the source's own word - a baby spinach and a
            # young carrot are not "mild" - while sharp, aged and vintage all
            # converge on "mature", which is what they mean on a cheese
            if label == "mild":
                return hit.group(0).strip(" ,").lower(), hit.group(0)
            return label, hit.group(0)
    return "", ""


# Flavour is sugar. A strawberry Greek yoghurt carries roughly twice the sugar
# of the plain one, and 106 of the 238 members of 'greek yoghurt, fat-free' were
# flavoured - FDC writes the flavour INSIDE the brand chunk ("CHOBANI STRAWBERRY
# NON-FAT"), where the brand strip took it with the brand.
#
# Read only on a product whose plain form is the unmarked one. Outside that set
# these words name the food rather than a flavour added to it: "Tomatoes,
# orange" is a colour, "Melon, banana" a variety and "Plantain banana" the fruit.
_FLAVOURED_HEAD_RE = re.compile(
    r"^\s*(?:gree[kc]\s+)?(?:yogh?o?urts?|kefir|quark|skyr|"
    r"ice\s+creams?|frozen\s+yogh?o?urts?|desserts?|puddings?|custards?|"
    r"mousses?|milkshakes?|smoothies?|cream\s+cheese)\b", re.I)
_FLAVOUR_WORDS = (
    r"strawberr(?:y|ies)|blueberr(?:y|ies)|raspberr(?:y|ies)|blackberr(?:y|ies)|"
    r"cherr(?:y|ies)|peach|apricot|mango|banana|pineapple|passion\s*fruit|"
    r"vanilla|chocolate|coffee|caramel|toffee|fudge|honey|maple|cinnamon|mint|"
    r"lemon|lime|orange|coconut|hazelnut|pistachio|almond|mixed\s+berry|"
    r"forest\s+fruit|fruit")
_FLAVOUR_PATTERNS = [
    ("flavour", re.compile(rf"\b(?:{_FLAVOUR_WORDS})\b", re.I)),
]


def _detect_flavour(d: str) -> tuple[str, str]:
    """Return (flavour label, matched token) for a flavoured product, else ("","")."""
    if not _FLAVOURED_HEAD_RE.match(d):
        return "", ""
    for _, rx in _FLAVOUR_PATTERNS:
        hit = rx.search(d)
        if hit:
            return _fold_plural(hit.group(0).lower()), hit.group(0)
    return "", ""


# An organ is not the muscle it was cut out of. FDC files them behind the head
# as "Beef, variety meats and by-products, liver", so the beef strip-head
# dropped the chunk and 60 organ rows landed in the 1,291-member 'beef' canon -
# 40 more in 'veal', 39 in 'pork', 296 rows across 44 canons. Beef liver carries
# roughly 9,000 ug RAE of vitamin A per 100 g against about none in muscle, and
# every cut of beef inside that canon inherited the liver's vitamin A.
#
# Only a WHOLE chunk counts, and a legume head is excluded outright: "Common
# bean, Kidney" is a kidney BEAN, and "Heart of palm" is a palm.
_ORGAN_WORDS = (r"livers?|kidneys?|hearts?|brains?|tongues?|lungs?|spleen|"
                r"pancreas|tripe|sweetbreads?|thymus|chitterlings?|marrow|"
                r"oxtails?|gizzards?|giblets?|offal|blood")
_ORGAN_PATTERNS = [
    # a leading comma is required: where the organ is the FIRST chunk it IS the
    # food ("Kidney, boiled, salted"), and taking it left the canon as 'salt'
    ("organ", re.compile(rf",\s*(?:{_ORGAN_WORDS})\s*(?=,|$)", re.I)),
]
_ORGAN_HEAD_SKIP_RE = re.compile(
    r"^\s*(?:beans?|peas?|lentils?|grams?|pulses?|cowpeas?|soy\w*|chickpeas?|"
    r"common\s+bean|palms?|hearts?\s+of\s+palm|artichokes?|cabbages?|"
    r"lettuces?|celery|celeriac)\b", re.I)


def _detect_organ(d: str) -> tuple[str, str]:
    """Return (organ label, matched token), or ("",""). See _ORGAN_PATTERNS."""
    if _ORGAN_HEAD_SKIP_RE.match(d):
        return "", ""
    for _, rx in _ORGAN_PATTERNS:
        hit = rx.search(d)
        if hit:
            return _singularize(hit.group(0).strip(" ,").lower()), hit.group(0)
    return "", ""


_PART_PRESERVE = {
    "cassava":  {"leaf", "leaves"},     # tuber is the food
    "amaranth": {"leaf", "leaves"},     # grain is the food
    "cowpea":   {"leaf", "leaves"},     # seed is the food
    "tamarind": {"leaf", "leaves"},     # pulp is the food
    "eggplant": {"leaf", "leaves"},     # fruit is the food
    "rice":     {"bran"},               # endosperm is the food
    "soybean":  {"sprout", "sprouts"},  # seed is the food
    "soybeans": {"sprout", "sprouts"},
}


# Canning / packing media, removed before state detection so "canned in juice"
# and "packed in water" do not register as juice products.
_PACK_MEDIUM_RE = re.compile(r"\bin\s+(?:juice|water|brine|oil)\b", re.I)
# "Soup, onion, dry, mix, prepared with water" is the RECONSTITUTED soup, not
# the sachet: the dry-mix wording describes the ingredient it was made from.
# When "prepared with" is present the dry-mix tokens are dropped, so the
# prepared food stays with its own kind and only the unprepared sachet splits.
_PREPARED_WITH_RE = re.compile(r"\bprepared\s+with\b", re.I)
_DRY_MIX_RE = re.compile(r"\bdry\b[\s,]*(?:mix)?", re.I)


# On a drink, "dry" is a style (not sweet), not a preparation state. Two
# independent audits found the same corruption: "Sherry, dry" -> 'sherry, dried',
# "Alcoholic beverage, wine, dessert, dry" -> 'wine, dried', "dry ginger ale" ->
# 'ginger ale, dried'. None of these drinks is dehydrated.
_DRY_STYLE_CONTEXT_RE = re.compile(
    r"\b(?:wine|sherry|vermouth|port|madeira|marsala|champagne|prosecco|cava|"
    r"riesling|cider|perry|sake|gin|martini|ale|beer|stout|lager|"
    r"ginger\s+ale|ginger\s+beer|tonic|soft\s+drink|alcoholic\s+beverage)\b", re.I)
_BARE_DRY_RE = re.compile(r"\bdry\b", re.I)
# "with" introduces what is IN the dish, not what was done to it.
# "without" is deliberately not matched - it is a claim about the dish.
_INGREDIENT_LIST_RE = re.compile(r"\s+with\s+.*$", re.I | re.S)


def _detect_state(d: str) -> tuple[str, str]:
    """Return (state label, the token that matched) for a description, or ("","")."""
    probe = _PACK_MEDIUM_RE.sub(" ", d)
    if _DRY_STYLE_CONTEXT_RE.search(probe):
        probe = _BARE_DRY_RE.sub(" ", probe)
    if _PREPARED_WITH_RE.search(probe):
        probe = _DRY_MIX_RE.sub(" ", probe)
    # Everything after " with " is an INGREDIENT list, and an ingredient's state
    # is not the dish's state: "Ciabatta sandwich w/ mozzarella cheese sundried
    # tomato lettuce" is not a dried food, it is a sandwich containing a dried
    # tomato. ~25 canons carried a false ", dried" from exactly this. Applied
    # after the dry-mix handling above, which needs to see "prepared with".
    probe = _INGREDIENT_LIST_RE.sub("", probe)
    for label, rx in _STATE_PATTERNS:
        hit = rx.search(probe)
        if hit:
            return label, hit.group(0)
    return "", ""


def _detect_trim(d: str) -> tuple[str, str]:
    """Return (trim label, matched token), or ("",""). See _TRIM_PATTERNS."""
    for label, rx in _TRIM_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


# Which oil matters. Sardines in olive oil, in sunflower oil and in peanut oil
# carry different fatty acids, and folding all three into "in oil" lets whichever
# is richest stand in for the other two. Water, brine and juice are not kept apart the same
# way: there the medium's own composition is close to nil.
_NAMED_OIL_RE = re.compile(
    r"\b(?:in|with)\s+((?:\w+\s+){1,2})oils?\b", re.I)


def _pack_wording(d: str) -> str:
    """The medium's wording, whether or not the label survives _detect_pack.

    Draining pours the medium away, so _detect_pack rightly refuses to LABEL a
    drained food by it - but the phrase is not part of the food's name either
    way. Left standing it renamed the food: "Apricot, canned in pear juice,
    drained" lost "canned" and "drained" to _PREP_RE and came out as 'apricot
    juice', which is a different food with a different sugar profile.
    """
    if not _PACKED_CONTEXT_RE.search(d):
        return ""
    d = _RECONSTITUTED_RE.sub(" ", d)
    for _, rx in _PACK_PATTERNS:
        hit = rx.search(d)
        if hit:
            return hit.group(0)
    return ""


# Oil and brine are media in their own right, and several sources name one
# without ever saying "canned": Fineli files "Cheese, Feta Cheese In Oil" and
# "Tomato, Sun-Dried, In Oil", CIQUAL "Olive, green, in brine", and Swiss
# Generic Foods "Tuna in oil, drained". All four failed the packed-context test
# and merged into the plain food - the same defect the comment below already
# argues against for canned tuna, where drained tuna in oil carries about eight
# times the fat of drained tuna in water. 8 canons.
#
# Water and juice are deliberately NOT admitted this way: without the packed
# context "boiled in water" is a cooking method, which is exactly what the
# guard was written to catch. Nor is a food that was COOKED in the oil - "Fast
# foods, potato, french fried in vegetable oil", "Pork, Strips, Fried In Oil" -
# where the oil is the frying medium and _detect_cook_fat owns the claim.
# "in", never "with". Without the packed context a "with <oil>" is an
# INGREDIENT, not a medium, and admitting it labelled four foods by something
# they are made of: "Puddings, banana, dry mix, with added oil" came out as
# 'pudding, banana, dried, in added oil', filled milk "with lauric acid oil" as
# 'milk, filled, in lauric acid oil', and agutuk - which is whipped WITH seal
# oil - as 'agutuk, in seal oil'. An explicit "added" is excluded for the same
# reason even after "in".
_MEDIUM_NO_CONTEXT_RE = re.compile(
    r"\bin\s+(?!added\b)(?:\w+\s+){0,2}(?:oils?|brine)\b", re.I)
_COOKED_IN_RE = re.compile(
    r"\b(?:deep[-\s]?fried|pan[-\s]?fried|stir[-\s]?fried|fried|frying|"
    r"saut(?:e|\u00e9)ed|browned|boiled|simmered|poached|blanched|heated|"
    r"grilled|griddled|roasted|baked|cooked|stewed|braised|steamed|"
    r"casseroled|toasted)\b", re.I)


def _detect_pack(d: str) -> tuple[str, str]:
    """Return (packing-medium label, matched token), or ("",""). See _PACK_PATTERNS."""
    if not _PACKED_CONTEXT_RE.search(d):
        if not (_MEDIUM_NO_CONTEXT_RE.search(d) and not _COOKED_IN_RE.search(d)):
            return "", ""
    d = _RECONSTITUTED_RE.sub(" ", d)
    if _DRAINED_RE.search(d) and not _PACK_PATTERNS[0][1].search(d):
        # Draining pours off water and brine, so the medium stops counting. Oil
        # is different: the flesh absorbs it and drained tuna in oil carries
        # roughly eight times the fat of drained tuna in water. Dropping the
        # medium for both put 107 rows of each into one 'tuna' canon.
        return "", ""
    for label, rx in _PACK_PATTERNS:
        hit = rx.search(d)
        if hit:
            if label == "in oil":
                named = _NAMED_OIL_RE.search(hit.group(0))
                if named:
                    label = f"in {named.group(1).strip().lower()} oil"
            return label, hit.group(0)
    return "", ""


def _detect_fat_whole(d: str) -> tuple[str, str]:
    """"whole" as a chunk on a dairy head is the full-fat level. See above."""
    if _DAIRY_HEAD_RE.match(d):
        hit = _WHOLE_CHUNK_RE.search(d)
        if hit:
            return "full-fat", hit.group(0).strip(" ,")
    return "", ""


def _detect_salt(d: str) -> tuple[str, str]:
    """Return (salt label, matched token), or ("",""). See _SALT_PATTERNS."""
    for label, rx in _SALT_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


# For a leafy vegetable or a herb the leaf IS the food, so labelling it would
# split the plain rows from the ones whose source spelled the default out -
# 354 lettuce records for nothing.
_LEAF_IS_DEFAULT_RE = re.compile(
    r"^\s*(?:lettuces?|spinach|kales?|chard|swiss\s+chard|cabbages?|rocket|"
    r"arugula|watercress|endives?|chicor(?:y|ies)|basil|parsle(?:y|ies)|mint|"
    r"coriander|cilantro|dill|sages?|thyme|oregano|rosemary|bay|tarragon|"
    r"chervil|sorrel|nettles?|purslane|collards?|pak\s*choi|bok\s*choy|"
    r"cress|radicchio|escarole|mizuna|shiso|curry\s+leaf|vine\s+leaf|teas?)\b",
    re.I)


def _detect_part(d: str) -> tuple[str, str]:
    """Return (anatomical-part label, matched token), or ("",""). See _PART_PATTERNS."""
    for label, rx in _PART_PATTERNS:
        hit = rx.search(d)
        if hit:
            if label == "leaf" and _LEAF_IS_DEFAULT_RE.match(d):
                return "", ""
            # the matched word IS the label: "Broccoli, stalks" reads better as
            # 'broccoli, stalk' than as the pattern's key 'stem'
            word = re.sub(r"\s+only$", "", hit.group(0).strip(" ,").lower())
            # BioFoodComp writes "pulp" where CIQUAL and FDC write "flesh"; the
            # two named the same part of the same fruit under two canons
            # "meat" here is only ever what is left of "meat only" - the "only"
            # is stripped a line above - and it says what "without skin" says.
            if (word in ("pulp", "meat", "skinless")
                    or word.startswith("without skin")):
                word = "flesh"
            return _fold_plural(word), hit.group(0)
    return "", ""


def _detect_moisture(d: str) -> tuple[str, str]:
    """Return ("0% moisture basis", token), or ("",""). See _MOISTURE_PATTERNS."""
    for label, rx in _MOISTURE_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


_INGREDIENT_CTX_RE = re.compile(r"\b(?:with|from|containing)\s+(?:\w+\s+){0,1}$", re.I)


# The fat level stated as a NUMBER. 1,163 rows carry one and _QUANT_RE deleted
# every one of them, so "Cheese, Edam, 17% Fat" and a 30% edam were one canon,
# and FDC's ground beef - sold at 5, 7, 10, 15, 20 and 30% fat - landed whole
# inside the bare 'beef'. Fat is the axis that exists to stop exactly that.
#
# The number is kept rather than banded, because no single threshold is right
# for two foods at once: 17% fat is a LOW-fat cheese and a high-fat yoghurt.
# "97% fat-free" is a claim about the other 97% and is excluded by name; 0% is
# left to the "fat-free" pattern, which already reads it and reads it better.
# A stated RANGE is kept whole - "Cottage Cheese, 2-5% Fat" is labelled
# '2-5% fat'. Reading only its top would be the very error the axis exists to
# prevent, and reading only its bottom understates the food.
_FAT_NUM = r"\d{1,2}(?:[.,]\d)?(?:\s*[-\u2013]\s*\d{1,2}(?:[.,]\d)?)?"
# The approximation word belongs to the figure, not to the food. CIQUAL writes
# "Tomme cheese, reduced fat, around 13% fat" and ">= 35% fat"; taking only the
# figure left "around" standing as a chunk of the name.
_FAT_PCT_RE = re.compile(
    rf"(?:\b(?:around|about|approx\.?|approximately|env\.?)\s+|[<>]=?\s*)?"
    rf"(?<![\d.,])({_FAT_NUM})\s*%\s*fat\b(?![-\s]?free)"
    rf"|\bfat\s*({_FAT_NUM})\s*%", re.I)


def _detect_fat(d: str) -> tuple[str, str]:
    """Return (fat-level label, matched token), or ("",""). See _FAT_PATTERNS.

    The EARLIEST marker wins, not the first pattern in the list. FDC writes
    "Milk, lowfat, fluid, 1% milkfat, with added nonfat milk solids": taking the
    patterns in list order found the "nonfat" of the solids and labelled 1% milk
    fat-free. Where two markers start in the same place the list order decides,
    which is what still keeps "part-skim" from being read as "skim".
    """
    hit = _FAT_PCT_RE.search(d)
    if hit and not _INGREDIENT_CTX_RE.search(d[:hit.start()]):
        num = re.sub(r"\s*[-\u2013]\s*", "-", (hit.group(1) or hit.group(2)).replace(",", "."))
        num = "-".join(n[:-2] if n.endswith(".0") else n for n in num.split("-"))
        if any(float(n) > 0 for n in num.split("-")):
            return f"{num}% fat", hit.group(0)
    best = None
    for i, (label, rx) in enumerate(_FAT_PATTERNS):
        for hit in rx.finditer(d):
            # A marker inside an INGREDIENT clause grades the ingredient, not
            # the food: "Babyfood, banana juice with low fat yogurt" is not a
            # low-fat babyfood, and reading it as one both mislabels the food
            # and takes the wording, leaving 'with yoghurt' behind.
            if _INGREDIENT_CTX_RE.search(d[:hit.start()]):
                continue
            if best is None or hit.start() < best[0]:
                best = (hit.start(), i, label, hit.group(0))
            break
    return (best[2], best[3]) if best else ("", "")


def _detect_oil_grade(d: str) -> tuple[str, str]:
    """Return (pressing-grade label, matched token) for a fat, or ("","")."""
    if not _OIL_CONTEXT_RE.search(d):
        return "", ""
    for label, rx in _OIL_GRADE_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


# The claim, denied. The word-level lookbehinds inside _FORTIFY_PATTERNS only
# guard the bare "fortified"; "Salt, not fortified with iodine" then matched on
# its "with iodine" half instead and came out as 'salt, enriched'.
_NOT_FORTIFIED_RE = re.compile(
    r"\b(?:not|un|non[-\s]?)\s*fortified\b|\bunenriched\b|\bnot\s+enriched\b", re.I)


def _detect_fortify(d: str) -> tuple[str, str]:
    """Return ("enriched", matched token), or ("",""). See _FORTIFY_PATTERNS."""
    if _NOT_FORTIFIED_RE.search(d):
        return "", ""
    for label, rx in _FORTIFY_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


def _detect_fibre(d: str) -> tuple[str, str]:
    """Return (fibre figure or high-fibre claim, matched token), or ("","").

    Mirrors the fat-percentage half of _detect_fat: the figure becomes the
    label, so it survives a strip head the way "13% fat" does.
    """
    hit = _FIBRE_PCT_RE.search(d)
    if hit and not _INGREDIENT_CTX_RE.search(d[:hit.start()]):
        num = re.sub(r"\s*[-–]\s*", "-", hit.group(0).split("%")[0].strip().replace(",", "."))
        num = "-".join(n[:-2] if n.endswith(".0") else n for n in num.split("-"))
        try:
            if any(float(n) > 0 for n in num.split("-")):
                return f"{num}% fibre", hit.group(0)
        except ValueError:
            pass
    hit = _HIGH_FIBRE_RE.search(d)
    if hit and not _INGREDIENT_CTX_RE.search(d[:hit.start()]):
        return "high-fibre", hit.group(0)
    return "", ""


def _detect_immature(d: str) -> tuple[str, str]:
    """Return ("immature", matched token) for a pre-maturity harvest, or ("","")."""
    hit = _IMMATURE_RE.search(d)
    return ("immature", hit.group(0)) if hit else ("", "")


def _detect_ripeness(d: str) -> tuple[str, str]:
    """Return (ripeness label, matched token), or ("",""). See _RIPENESS_PATTERNS."""
    for label, rx in _RIPENESS_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


def _detect_decaf(d: str) -> tuple[str, str]:
    """Return ("decaffeinated", matched token), or ("","")."""
    hit = _DECAF_RE.search(d)
    return ("decaffeinated", hit.group(0)) if hit else ("", "")


def _detect_refined_grain(d: str) -> tuple[str, str]:
    """Return ("refined", matched token) for a milled cereal, or ("","").

    Gated to a cereal and away from a fat: "refined" on an oil is a pressing
    grade, which _detect_oil_grade owns.
    """
    if not _CEREAL_CONTEXT_RE.search(d) or _OIL_CONTEXT_RE.search(d):
        return "", ""
    # No ingredient guard here, unlike the fat and fibre figures: on a cereal the
    # "from <grain>" clause names what the food IS, not something added to it -
    # "Porridge, soft, from degermed yellow maize grit" is a degermed-maize
    # porridge - and the guard read it as an ingredient and dropped the claim.
    hit = _REFINED_GRAIN_RE.search(d)
    if hit:
        return "refined", hit.group(0)
    hit = _WHOLEGRAIN_RE.search(d)
    if hit:
        return "wholegrain", hit.group(0)
    return "", ""


def _detect_gluten_free(d: str) -> tuple[str, str]:
    """Return ("gluten-free", matched token), or ("","")."""
    hit = _GLUTEN_FREE_RE.search(d)
    return ("gluten-free", hit.group(0)) if hit else ("", "")


def _detect_sugar(d: str) -> tuple[str, str]:
    """Return (sugar label, matched token), or ("",""). See _SUGAR_PATTERNS."""
    for label, rx in _SUGAR_PATTERNS:
        hit = rx.search(d)
        if hit:
            return label, hit.group(0)
    return "", ""


# The labels _append_states actually PRINTED on the canon being built. Read by
# canonicalize_food_name to look the curated tables up a second time with the
# labels lifted off; see the note there. Module-level because the labels are
# produced deep inside _canonicalize_core and consumed after it returns.
_PRINTED_LABELS: list[str] = []
# A curated name that already states the fortification in words. See the
# label-aware lookup in canonicalize_food_name.
_FORT_SAID_RE = re.compile(r"\b(?:fortified|enriched|with\s+added)\b", re.I)


def _is_rename(key: str, value: str) -> bool:
    """Is this table entry a RE-SPELLING of the same name, or a substitution?

    Only a re-spelling may be carried onto a labelled sibling. An entry that
    swaps the food for another one is written for the exact row it names:
    'wheat' -> 'cereal, ready-to-eat, shredded wheat, lightly frosted' is FDC's
    breakfast cereal, and letting the label-aware lookup reach it turned
    "wheat, bran" into that cereal with bran on the end.
    """
    w = lambda x: {_singularize(t) for t in re.findall(r"[a-z0-9]+", x.lower())}
    if w(key) == w(value):
        return True
    # ...and a respelling that only moves a space is still a respelling. The
    # word-set test cannot see it - {'pigeonpea'} against {'pigeon','pea'} - so
    # 'pigeonpea' -> 'pigeon pea' and 'broadbean' -> 'broad bean' stopped firing
    # the moment a label was appended, and the canon shipped as 'pigeonpea,
    # immature' beside the curated 'pigeon pea'. Compared on letters alone,
    # which is narrow enough that 'wheat' cannot reach FDC's frosted cereal.
    flat = lambda x: re.sub(r"[^a-z0-9]", "", x.lower())
    return flat(key) == flat(value)


def _strip_whole_chunk(d: str, token: str) -> str:
    """Remove a token only where it is an ENTIRE comma chunk.

    The unconditional strip the fat and fortification axes use is safe because
    their wordings are phrases in their own right. A refinement or gluten claim
    is usually an ADJECTIVE on the food's own noun, and cutting it out of the
    middle of a chunk destroys the name the downstream head rules read:
    "Flour, whole grain oat" became 'flour oat' instead of 'oat flour,
    wholegrain' (217 rows), "Bread, whole wheat-based" became 'bread-based',
    and "Rye Bread, Wholegrain Rye, Dark Wheat Flour" became 'rye bread rye'.
    Where the wording is NOT its own chunk it is simply left alone: either it
    survives into the canon, and _append_state then declines to print a label
    the name already carries, or the two-chunk rule drops it and the label
    prints. Both outcomes are right; only cutting it out mid-chunk is wrong.
    """
    if not token:
        return d
    want = re.sub(r"[^a-z0-9]", "", token.lower())
    if not want:
        return d
    out, hit = [], False
    for c in d.split(","):
        # the source tag has not been stripped yet where this fires, so the
        # last chunk can be "Wholegrain [Fineli]" and still be the whole chunk
        bare = re.sub(r"[^a-z0-9]", "", _BRACKET_TAG_RE.sub("", c).lower())
        if not hit and bare == want:
            hit = True
            continue
        out.append(c)
    return ",".join(out) if hit and any(c.strip() for c in out) else d


def _append_states(canon: str, *states) -> str:
    """Apply several independent state labels, each guarded by _append_state."""
    for st in states:
        out = _append_state(canon, st)
        if out != canon and st and st[0]:
            _PRINTED_LABELS.append(st[0].strip().lower())
        canon = out
    return canon


# Labels that are detected but never printed. Unsalted is the unmarked case:
# a raw carrot and an "unsalted" carrot are the same food, so "carrot, unsalted"
# and "carrot" were two canons for one thing. Detection still MATTERS, because
# _SALT_PATTERNS checks unsalted first and that is what stops "unsalted" being
# read as "salted" - so the split against the genuinely salted version survives.
# Unsweetened joins it for the same reason: an unsweetened apple sauce and a
# plain one are the same food, so printing the label would split a canon that
# detection has already protected from the sweetened version.
_SILENT_STATE_LABELS = frozenset({"unsalted", "unsweetened"})


def _append_state(canon: str, state: tuple[str, str]) -> str:
    """Re-attach a detected state to a canon name, unless the canon already
    carries it. Both spellings have to be checked: the label ("chickpea, flour"
    must not become "chickpea, flour, flour") and the token that actually
    matched ("chili powder" must not become "chili powder, dried", since the
    word surviving in the canon is "powder", not the label "dried")."""
    label, token = state
    if not label or not canon or label in _SILENT_STATE_LABELS:
        return canon
    for w in (label, token.strip()):
        if w and re.search(rf"\b{re.escape(w)}\b", canon, re.I):
            return canon
    return f"{canon}, {label}"


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
# The alcohol arm of the same class prefix. "Beverages, almond milk" already
# loses its class noun above, but "Alcoholic beverage, beer" kept one, so 14
# canons were headed by the category rather than the drink and 'beer' and
# 'alcoholic beverage, beer' were two foods.
#
# The prefix goes only when what follows NAMES a drink. Two of the fourteen do
# not: "Alcoholic beverage, rice (sake)" is rice wine and stripping the class
# would file sake under the grain, and "Alcoholic beverage, distilled, all (gin,
# rum, vodka, whiskey) 86 proof" is the unnamed spirits average. Both are
# renamed by hand in _CANON_OVERRIDES instead.
# Wild rice is Zizania, not Oryza: a different grain carrying roughly twice the
# protein and four times the fibre of white rice. Five sources write it "Rice,
# wild" and rice is a strip head, so every one of those 21 rows was landing
# inside the 361-member 'rice' canon, and the sources that write it the other
# way round had a canon of their own. This is a NAME, not the provenance axis
# below - which is why rice is absent from the farmed-food gate there.
_WILD_RICE_RE = re.compile(r"\brices?\s*,\s*wild\b", re.I)

_ALCOHOL_CLASS_RE = re.compile(
    r"^\s*alcoholic\s+beverages?\s*,\s*(?=(?:"
    r"beer|ale|lager|stout|porter|wine|champagne|prosecco|cider|perry|sake|"
    r"mead|vermouth|sherry|port|liqueur|liquor|spirits?|whisk(?:e)?y|vodka|"
    r"rum|gin|tequila|mezcal|brandy|cognac|armagnac|absinthe|schnapps|aquavit|"
    r"akvavit|ouzo|sambuca|grappa|pisco|soju|shochu|baijiu|cocktail|daiquiri|"
    r"margarita|martini|mojito|pina\s+colada|creme\s+de\s+menthe|"
    r"black\s+russian|cosmopolitan|mimosa|sangria|screwdriver|bloody\s+mary|"
    r"tom\s+collins|malt\s+beer)\b)", re.I)
# Research-paper-title heuristic: descriptions > 120 chars containing any of
# these journal-style keywords. Used in build_static_food_meta to skip the
# ~50 BioFoodComp rows where the source description column leaked the
# literature reference instead of the food name.
# Matched against the lower-cased description, so re.I is load-bearing here even
# though the source spelling is "Kopia av".
_COPY_RECORD_RE = re.compile(COPY_RECORD_RE, re.I)

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
# FDC's experimental_food type is a BIBLIOGRAPHY, not a food list. Of its 114
# rows, 113 are literature references ("Effects of Prune Supplementation on
# Cardiometabolic Health in Postmenopausal Women: ...", "USDA, FDA, and ODS-NIH
# Database for the Iodine Content of Common Foods (Release 4.0)"). Exactly one is
# a food: "Wild Turkey, White and Dark Meat, Raw and Roasted", at 49 characters.
# The shortest title is 52, so length separates them cleanly on this type.
#
# This replaced a keyword sieve. Chasing titles by vocabulary does not converge -
# a list broad enough to catch "diversity in", "as affected by" and "correlates
# with" also caught "Fruits dessert, all types (sugar content is less than fruits
# compote...)", which is a real food. The data type is the honest signal.
_EXPERIMENTAL_TITLE_LEN = 50

# "95% extraction" is a flour MILLING rate, not solvent extraction, and the
# keyword sieve above has always mistaken it for one: two real PhyFoodComp maize
# foods were being dropped from the index before this exemption existed. Seven
# rows carry a milling rate; the two long enough to trip the 60-char gate are the
# ones that were lost.
_MILLING_EXTRACTION_RE = re.compile(r"\d+\s*%\s*extraction", re.I)

# Strip ALL parenthetical content (clarifications, color descriptors, source
# notes - "(industrial)", "(colour of peel: olive green)", "(fat free or
# skim)", "(includes foods for USDA's food distribution program)"). Parens
# in FDC descriptions are nearly always non-essential annotations, never
# identity-defining.
_PAREN_CONTENT_RE = re.compile(r"\s*\([^)]*\)\s*")
# Strip quantitative qualifiers like "9% protein", "50% extraction",
# "3.25% milkfat". Limited trailing word run so we don't over-eat.
# Quantities and their property. The trailing run is capped at 25 characters,
# which is right for "37% fat" and "0% moisture basis". It over-reaches on an
# INGREDIENT stated as a percentage - Livsmedelsverket's "Crepe with 37% shrimp
# stuffing heated" loses the filling - but only 121 rows across the corpus are
# affected and most of those are marketing text ("100% natural italian"), so
# the narrower forms tried against it left worse names behind than they fixed.
# A FIBRE percentage is composition, not a qualifier, and it was the only
# nutrient figure this rule ate: "Cheese, 13% fat" keeps its figure because
# _FAT_PCT_RE reads it first and puts it back, while "Crispbread, Rye, 17%
# Fibre" came out as bare 'crispbread, rye'. 57 rows lost the figure, and the
# curated table had already ruled the other way on the same wording - the
# override for 'crisp bread wholegrain rye 15% fibre flatbröd' spells the canon
# 'crisp bread, wholegrain rye, 15% fibre' - so the rule was overwriting a
# decision a curator had made by hand. Fibre is the double-weighted nutrient in
# the score and the score takes the max, so folding a 17%-fibre crispbread into
# the plain one hands the plain one the figure outright.
# _FIBRE_PCT_RE is defined with the other composition axes above.
_QUANT_RE = re.compile(r"\b\d+(?:\.\d+)?(?:\s*[-–]\s*\d+(?:\.\d+)?)?\s*%\s*(?!fib(?:re|er)\b)[A-Za-z][A-Za-z\s-]{0,25}?(?=,|$)")
# Orphan connectives left behind by the prep strip. "Durian, raw or frozen"
# loses "raw" and "frozen" but keeps the "or"; "Bread, with seeds" loses
# "seeds" but keeps the "with". Either way a chunk survives that is nothing
# but a joining word. 39 canons across the corpus, none a real food name.
_CONNECTIVE     = r"(?:or|and|with|without|in|of|the|a|an|to|from|on|for|by)"
_CONN_ONLY_RE   = re.compile(rf"^{_CONNECTIVE}(?:\s+{_CONNECTIVE})*$", re.I)
_CONN_TRAIL_RE  = re.compile(rf"(?:\s+{_CONNECTIVE})+$", re.I)

# Plural folding, applied to the last word of each comma chunk. Without it the
# same food arrives under two spellings and scores twice: apples/apple splits
# 1,857 foods, tomatoes/tomato 1,534, strawberries/strawberry 603. 223 such
# pairs split 17,963 foods.
# The -es rule fires only after a sibilant (peaches -> peach, radishes ->
# radish) so that olives falls through to the plain -s rule and gives olive
# rather than "oliv". -ves is a table, not a pattern, for the same reason:
# the pattern would take olives -> "olif".
_VES_SINGULAR = {"leaves": "leaf", "halves": "half", "loaves": "loaf",
                 "calves": "calf", "hooves": "hoof", "knives": "knife",
                 "wolves": "wolf"}
# Singular words ending in s. The ss/us/is endings are handled by rule; these
# are the ones that are not, plus mass nouns whose singular is not a food
# ("turnip greens" -> "turnip green").
_NOT_PLURAL = frozenset({"molasses", "christmas", "shoes", "grits", "greens",
                         "species", "brussels", "swiss",
                         # cheeses and drinks whose names end in -s in the
                         # singular; taking the s off gives 'causs', 'maroille'
                         # and 'bitter', none of which is the food
                         "causses", "maroilles", "bitters", "speculoos",
                         "pommes", "frijoles", "sports", "childrens",
                         "dolichos", "foothills", "selects"})
# The -is guard below protects Latin epithets (ensiformis, aestivalis), French
# (frais, coulis, liegeois) and Finnish (haerkis, taikaruis) singulars, 28 words
# in all. These four are the only genuine plurals it catches by mistake.
_IS_PLURAL = frozenset({"litchis", "rostis", "chapatis", "blinis"})
_SIBILANT_ES_RE = re.compile(r"(?:s|x|z|ch|sh)es$")

# CSV quoting that leaked into the stored description. 1,720 descriptions carry a
# stray double quote, and it survives into the canon: "Beef, flank, steak, trimmed
# to 0"" fat" gives '"beef, flank', which sorts and groups separately from 'beef'.
_QUOTE_RE = re.compile(r'["\u201c\u201d]')

# Heads that name a place or a drinks aisle rather than a food. FDC files restaurant
# dishes under "Restaurant, <cuisine>, <dish>" and CNF files spirits under
# "Alcohol, <spirit>", so the canon came out as 'restaurant, italian' and
# 'alcohol, whisky'. The cuisine list is closed: FDC uses exactly these five.
# Restaurant chains and meal programmes that FDC files as a PREFIX to the dish:
# "Carrabba's Italian Grill, spaghetti with meat sauce", "School Lunch, chicken
# nuggets", "Domino's 14\" cheese pizza". The prefix names who served the food, not
# what it is, so it splits one dish across as many canons as there are vendors.
#
# Curated rather than derived from a possessive: "cat's whisker" and "jew's mallow"
# are vegetables, "shepherd's pie" and "Lindström's steak" are dishes, and
# "chipotle dip" is an ingredient, not the chain. Every name below was confirmed to
# head at least three different dishes in food.parquet. Bare "chipotle" is
# deliberately absent for that reason.
_VENUE_BRANDS = (
    r"mcdonald's|wendy's|denny's|applebee's|arby's|popeyes?|burger\s+king|"
    r"taco\s+bell|pizza\s+hut|little\s+caesars|papa\s+john's|domino's|kfc|"
    r"subway|chick-fil-a|olive\s+garden|cracker\s+barrel|on\s+the\s+border|"
    r"t\.g\.i\.?\s*friday's|carrabba's\s+italian\s+grill|campbell's(?:\s+chunky)?"
)
# The prefix must be followed by a comma or by a pack size ("14 cheese pizza", the
# quote already removed by _QUOTE_RE). Requiring one or the other keeps a bare brand
# word from eating the head of a food that merely starts with the same letters.
_VENUE_PREFIX_RE = re.compile(
    r"^\s*(?:"
    r"(?:" + _VENUE_BRANDS + r")"
    r")"
    # Separator: a comma ("KFC, biscuit"), a hyphen ("BURGER KING - HAM"), a pack
    # size ("Domino's 14 cheese pizza", the quote already gone), or plain
    # whitespace ("McDONALD'S Bacon Ranch Salad"). Whitespace is safe only
    # because every name in _VENUE_BRANDS is a brand rather than a word that
    # could open a food; the generic heads below it still require punctuation.
    r"(?:\s*[,-]\s*|\s+\d+\s*['\u2019]?\s+|\s+)",
    re.I)

# The generic heads cannot take a bare space: "restaurant" and "school lunch"
# are ordinary words and "alcohol" is a nutrient, so they need punctuation.
_VENUE_GENERIC_RE = re.compile(
    r"^\s*(?:"
    r"restaurant\s*,\s*(?:chinese|latino|family\s+style|mexican|italian)"
    r"|chinese\s+restaurant|restaurant"
    r"|fast\s+foods?|school\s+lunch|alcohol"
    r")\s*[,-]\s*",
    re.I)

# Label copy that is not composition. "Pillsbury Golden Layer Buttermilk Biscuit,
# Artificial Flavor" is one biscuit, not a variety of one. Added sugar, added
# calcium and added vitamins are deliberately NOT here: those do change what a
# gram of the food contains, so they keep their own canon.
_LABEL_ONLY_RE = re.compile(
    r"^(?:artificial|natural)(?:ly)?\s+flavou?r(?:s|ed|ing)?$", re.I)

# Chunks left as pure punctuation. The prep strip cuts inside slash alternations
# ("Anglerfish, raw/ without added fat, fried" loses "raw" and keeps "/"), which
# leaves 68 canons carrying a chunk that is only a separator.
_PUNCT_ONLY_RE = re.compile(r"^[\s/\\|;:.\u2013\u2014<>\[\]()*&+-]+$")

# Characters that carry no meaning at the edge of a chunk. Kept off this list:
# "%" (real: "acetic acid 12%") and any interior character, since only the ends
# are trimmed ("filet-o-fish" and "m&m's candy" survive intact).
_POLISH_EDGE = " \t,;:/\\|.\u2013\u2014<>=[](){}*&+-"

# Chunks that classify rather than name. "Vermouth, dry type" keeps "type" once
# the dry-type collocation is excluded from the dried state; "Beef, all cuts"
# and "Peanuts, all types" read as 'beef, all cut' and 'peanut, all type'. The
# match must cover the WHOLE chunk, so a real qualifier survives - "butter,
# unknown fat content" is a genuine aggregate entry and is left alone.
_FILLER_CHUNK_RE = re.compile(
    r"^(?:all\s+|assorted\s+|mixed\s+)?(?:type|kind|variet(?:y|ies)|form|sort|class|grade|cut)s?$"
    # "general" and "assorted" are NOT filler: "Restaurant food, general" is an
    # aggregate entry distinct from the named restaurant dishes, and dropping
    # the word would merge it into them.
    r"|^(?:average(?:\s+value)?|unspecified|miscellaneous)$"
    # FDC grades some oils by what they are FOR rather than what they are:
    # "Oil, soybean, salad or cooking" is the same oil as "Oil, soybean".
    r"|^(?:salad\s+or\s+cooking|cooking\s+or\s+salad|salad\s+and\s+cooking)$"
    # grain size is not composition; it appears only on sugar and salt
    r"|^granulated$"
    # FDC marks brine-injected poultry "enhanced" and the plain kind
    # "non-enhanced"; the negative is the unmarked case and names nothing
    r"|^non[-\s]?enhanced$"
    # texture and retail format, not composition: FDC grades salsa "chunky" and
    # sells cheese as a "chunk" or a "block"
    r"|^(?:chunky|chunk|block|wedge|creamy|smooth|crunchy|solids?)$"
    # "meat and skin" is FDC's whole bird, which is what a bare "chicken" is;
    # spelling the default out took the qualifier slot the CUT needed, so the
    # breast and the thigh both came out as plain 'chicken, meat and skin'.
    r"|^meat\s+and\s+skin$"
    # McCance splits the same statement over two chunks - "Duck, raw, meat, fat
    # and skin" - so the half left standing has to go too
    r"|^(?:fat\s+and\s+skins?|skins?\s+and\s+fat)$"
    # FDC files its bean breeding accessions as a chunk of their own -
    # "Beans, Dry, Pink, 11F-8082 (0% moisture)" - and 44 of them became canons.
    # Letters BETWEEN digits are what makes it a code: "0-50 mg calcium per
    # litre" and "20-30 g fat" are quantities and are left alone.
    r"|^\d{1,3}[a-z]{1,2}[-\u2013]\d{2,}$"
    # "store brand" is the ABSENCE of a brand, and "mixed composite" is how the
    # sample was drawn. Both took a qualifier slot: 91 rows of white bread came
    # out as 'white bread, store brand'.
    r"|^(?:store|name|national|other|own)[-\s]?brands?"
    r"(?:\s*/\s*(?:store|name|national|other|own)?[-\s]?brands?)?$"
    r"|^storebrands?$|^(?:mixed\s+)?composite$"
    # Fragments left behind when a preparation word is stripped out of a
    # compound: "sun-dried" -> "sun", "seeds removed" -> "removed",
    # "not fortified" -> "not", "pre-cooked" -> "pre". None names a food.
    r"|^(?:not|only|just|fully|removed|purchased|weighed|added|pre|post|semi|"
    r"deep|shallow|home|freeze|spray|sun|air|vacuum|single|regular)$"
    # A bare number with no unit: NEVO grades cheese "45+" and the "+" is
    # stripped, leaving canons like "cheese 45" / "cheese brie 50".
    r"|^\d+\+?$", re.I)

# WAFCT / Burkina Faso entries pair a local name with an English gloss:
# "babenda-1 *: sauce from green leaf, groundnut powder". The " *: " is a
# separator between two names, not noise, so it becomes a comma and the normal
# two-chunk rule does the rest. Any asterisk left over is a footnote marker.
# Only the gloss's FIRST chunk is kept, and it is trimmed here rather than by
# the later two-chunk rule, because state detection runs in between: the
# discarded tail of "babenda-1 *: sauce from green leaf, groundnut powder"
# contains "powder", which the dried-state probe would read as the state of
# a sauce that is not dried.
_LOCAL_GLOSS_RE = re.compile(r"\s*\*\s*:\s*([^,]*).*$")
_FOOTNOTE_STAR_RE = re.compile(r"\*+")
# "Bread, multigrain and/or with seeds" loses "with seeds" to _PREP_RE and keeps
# the compound conjunction, which _CONN_TRAIL_RE cannot see: it expects
# whitespace-separated words and "and/or" is a single token.
# NEVO grades cheese "45+"; the "+" is stripped upstream and the bare number
# is left welded to the name ("cheese 45", "cheese brie 50").
_TRAIL_GRADE_RE = re.compile(r"\s+\d{1,3}\+?$")
_LINE_NUMBER_RE = re.compile(
    r"(?<=\s)\d{1,2}\s+(?!%|g\b|mg\b|ug\b|ml\b|cl\b|oz\b|lb\b|kg\b|cm\b|vol\b|grain\b|month\b|year\b)(?=[a-z])", re.I)
# A serving diameter is not composition: FDC files pizzas as '14" pizza',
# which became the canon and pushed the topping out of the name entirely.
_DIAMETER_RE = re.compile(r"\b\d{1,2}\s*(?:\"\"|\"|inch|in\.)?\s*(?=pizza\b|sub\b|sandwich\b)", re.I)
_ANDOR_TRAIL_RE = re.compile(r"\s+(?:and|or)\s*/\s*(?:and|or)\s*$", re.I)

# STFCJ writes synonym tags with full-width CJK brackets, which no ASCII
# bracket rule can match: "fish, cod, walleye pollock* \uff3b*syn. alaska pollock\uff3d"
# kept the tag and then shed only half of it, leaving 27 canons with a dangling
# "\uff3d" and one canon that was the single character "\uff0d".
_FULLWIDTH_MAP = str.maketrans({
    "\uff3b": "[", "\uff3d": "]", "\uff08": "(", "\uff09": ")",
    "\uff0c": ",", "\uff1b": ";", "\uff1a": ":", "\uff0f": "/",
    "\uff0d": "-", "\u3000": " ",
    # STFCJ also marks its synonym footnotes with the FULL-WIDTH asterisk, which
    # _FOOTNOTE_STAR_RE cannot see: 'horse mackerel, japanese jack mackerel\uff0a'.
    "\uff0a": "*",
    # Fineli and Frida write the possessive with an acute accent or a curly
    # quote rather than an apostrophe. The brand list matches "kellogg" and
    # strips it, but the "\u00b4S" that follows survives as a word, so
    # "Kellogg\u00b4S Frost Rice Krispies" canonicalised to '\u00b4s frost rice krispy'.
    "\u00b4": "'", "\u2019": "'", "\u2018": "'", "\u0060": "'",
})

# NEVO and CIQUAL abbreviate "with" as w and "without" as wo ("Apple w skin av",
# "Bechamel sauce, w butter", "Apple pie Dutch w shortbread wo butter"). Left
# untranslated they produced 487 canons ending in a dangling "w". Only a token
# with whitespace on BOTH sides is the abbreviation - "S&W PREMIUM",
# "DMR-ESR-W", "CA9890234W" and "CO92059-8W" all carry a W that is part of a
# code, and each is preceded by "&", "-" or a digit rather than a space.
# Livsmedelsverket (SWE) writes the same abbreviation with a slash - "Apple w/
# skin", "Asparagus green boiled w/ salt" - 741 rows, and "w/o" for without.
# Handled before the bare forms so "w/o" is not read as "w" + stray slash.
_WO_ABBREV_RE = re.compile(r"(?<= )(?:wo|w/o)(?=[\s,]|$)", re.I)
_W_ABBREV_RE  = re.compile(r"(?<= )w/?(?=[\s,]|$)", re.I)
# NEVO closes aggregate entries with "av" (average): "Apple w skin av". It names
# no food, and left in place it produced canons like 'apple av'. Trailing only -
# anchored to a comma or end of string so no interior word is touched.
_AVERAGE_TOKEN_RE = re.compile(r"\s*\b(?:av|avg)\b\.?(?=\s*(?:,|$))", re.I)

# A parenthetical opening with with/without/no/plain names an INGREDIENT
# difference, not an annotation. _PAREN_CONTENT_RE deleted it wholesale, which
# silently merged "BIG MAC" with "BIG MAC (without Big Mac Sauce)",
# "FILET-O-FISH" with "(without tartar sauce)", "Hotcakes (plain)" with
# "Hotcakes (with 2 pats margarine & syrup)" and two more - pairs whose fat
# differs by exactly the omitted component. Promoting it to a comma chunk keeps
# the distinction and reads correctly: "big mac, without big mac sauce".
_PAREN_QUALIFIER_RE = re.compile(r"\s*\(\s*((?:with|without|no|plain)\b[^)]*)\)", re.I)

# FDC sample/analysis provenance that survived _NF_BARE_RE because it is neither
# an NF nor a CY code: "fat, chicken, skin, braised ... - 14b-03-04-totalfat".
_SAMPLE_CODE_RE = re.compile(r"\s*[-,]\s*\b\d+[a-z]?(?:-\d+)+-[a-z]+\w*\s*$", re.I)
# "Mild chicken strip, analyzed 2006" - the year is provenance, not identity.
_ANALYZED_YEAR_RE = re.compile(r"\s*[-,]?\s*\banaly[sz]ed\s+\d{4}\b", re.I)

# A venue is not always a prefix: "Sweet And Sour Pork, Chinese restaurant" puts
# the dish first and the place second, where no ^-anchored rule can reach it. A
# chunk that is nothing but a venue name carries no food information wherever it
# sits, so it is dropped outright.
_VENUE_CHUNK_RE = re.compile(
    r"^(?:" + _VENUE_BRANDS + r"|chinese\s+restaurant|restaurant|fast\s+foods?|"
    r"pizza\s+chain|school\s+lunch)$", re.I)

# FDC subsample records repeat one food once per nutrient panel and put the PANEL
# first: "Minerals, Salsa, TOSTITOS CHUNKY, MEDIUM - NFY090KVS". Unstripped, the
# canon becomes the nutrient ('mineral', 'niacin', 'fatty acid'), which is not a
# food at all. 2,400 foods across 24 panel heads plus the group names below.
# Guarded two ways, because several of these are also real foods:
#   - only for FDC rows (no source tag). Fineli and McCance title-case their own
#     food names, so "Sugar, Demerara [McCance]" must not trip the rule.
#   - only when the NEXT chunk is capitalized or is itself a panel head. That is
#     what separates "Sugars, Salsa, TOSTITOS" from "Sugars, brown", "Salt, table",
#     "Water, non-carbonated" and "Pectin, liquid", which are foods.

# FDC files some rows under the analyte that was measured, joined to the food by
# " - " or a comma: "FA - Beef, porterhouse steak", "B12, B6, B3, B2 - Beef, Eye
# of Round", "Choline - Beef, top loin steak", "Total lipid, Sugar, Granulated".
# The canon then names a nutrient, not a food. The vocabulary is closed and
# deliberately excludes every analyte word that is also a food (salt, sugar,
# water, starch), so a real food name can never match.
_ANALYTE = (r"(?:f\.?a\.?|tdf|chole|cholesterol(?:[-\s]\w+)?|choline|betaine|"
            # element symbols, only ever matched with the " - " separator below
            r"se|fe|zn|ca|mg|na|cu|mn|cr|mo|se\s*&\s*\w+|"
            # ... and spelled out, which "se" above does not cover
            r"selenium|iodine|iron|zinc|calcium|magnesium|phosphorus|"
            r"potassium|copper|manganese|chromium|molybdenum|"
            r"carotenoids?|lycopene|lutein|isoflavones?|phytosterols?|"
            r"total\s+(?:lipid|fat|dietary\s+fibre?|sugars?)|"
            r"vit\.?\s*[a-ek]\d*|vitamins?\s*[a-ek]?\d*|b\d{1,2}|"
            r"niacin|thiamin[e]?|riboflavin|folates?|retinol|tocopherols?|"
            r"proximates?|minerals?|fatty\s+acids?|amino\s+acids?|"
            r"moisture|ash|energy|carbohydrates?|pantothenic\s+acid|"
            # FDC's own spelling slips ("panthothenic") and the panels the NF
            # rows use that were missing: starch heads 60 rows, and the
            # vitamin D panel writes both forms of the analyte in one label
            r"panthothenic\s+acid|"
            r"vitamins?\s*d3?\s*(?:and|&)\s*25[-\s]?oh[-\s]?d3?)")
# "Sugars," and "Sweets," are FDC/CNF panel labels on lab rows ("Sugars, Cheese,
# swiss, slices (CA2, CO) - 18c-17-03-Sug") but ALSO real foods ("Sugars, brown").
# The lab rows are the ones carrying a sample code, so that is the discriminator.
# Applied once, at the very start, before the analyte prefixes are stripped -
# otherwise "Minerals, Sugar, Granulated, White - NFY040XEG" would lose BOTH
# labels and canonicalize to the bare adjective "granulated".
# Country of origin does not change composition per 100 g. Restricted to meat
# heads on purpose: "swiss" in "cheese, swiss" and "danish" in "pastry, danish"
# name the food, not where it came from.
# CIQUAL and BioFoodComp write the skinless meat as a bare "meat" chunk where
# FDC and McCance write "meat only", so ten animals carried both spellings:
# 'chicken, meat' beside 'chicken, flesh', and the same for duck, goose, hare,
# pheasant, pigeon, poultry, quail, turkey and the chicken leg. Rewritten to the
# phrase the part axis already reads.
#
# Two guards. The head must be the ANIMAL and not a dish made of it - "Babyfood,
# meat, beef", "Sauce, meat", "Pie, meat" and "Frankfurter, meat" are the 30-odd
# rows where the word names an ingredient. And the source must not have said the
# skin or the fat is in as well: "Duck, raw, meat, fat and skin" is the whole
# bird, which carries about five times the fat of the flesh, and it was landing
# in the same canon as CIQUAL's skinless duck.
# CNF files 171 traditional-foods rows with "native" as a chunk of the name -
# "Fish, burbot (loche), native, raw", "Game meat, native, bear, raw" - and the
# word marks the dataset, not the food: a native arctic char and an arctic char
# are the same fish. Twelve foods carried both spellings, and in the game rows
# it took the head outright, giving 40 canons like 'native, caribou, liver'.
#
# Two foods keep it, because there "native" is half of a species name: FDC's
# native persimmon is Diospyros virginiana against the Japanese kaki, and AFCD's
# native oyster is the flat oyster, not the Pacific one.
_NATIVE_RE = re.compile(r"\s*,?\s*\bnative\b", re.I)
_NATIVE_IS_SPECIES_RE = re.compile(r"^\s*(?:persimmons?|oysters?)\b", re.I)

_MEAT_DISH_HEAD_RE = re.compile(
    r"^\s*(?:baby\s*food|bouillon|broth|stock|soups?|sauces?|gravy|gravies|"
    r"p[a\u00e2]t[e\u00e9]s?|pies?|pasty|pasties|frankfurters?|sausages?|cakes?|"
    r"sandwiche?s?|rolls?|spring\s+rolls?|seasoning|stews?|casseroles?|curr(?:y|ies)|"
    r"balls?|loaf|loaves|burgers?|patt(?:y|ies)|nuggets?|dumplings?|fillings?|"
    r"spreads?|salads?|confit|"
    # ...nor on a plant, where "meat" is the SUBSTITUTE and not a part of an
    # animal: "Soy, meat" is textured soy protein and came out as 'soy, flesh'
    r"soya?|tofu|tempeh|seitan|gluten|mycoprotein|quorn|vegan|vegetarian|"
    r"plant[-\s]based)\b", re.I)
# the source tag and FDC's "(Alaska Native)" have not been stripped from `d`
# where this fires, so the chunk can be the last one and still be followed by
# something: "Bear, black, meat (Alaska Native)"
_MEAT_CHUNK_RE = re.compile(r"(?<=,)\s*meat\s*(?=\s*(?:,|\(|\[|$))", re.I)
_SKIN_INCLUDED_RE = re.compile(
    r"\b(?:and|with|\+|,)\s*(?:fat\s+and\s+)?skins?\b|\bskin\s+(?:and|on)\b|"
    r"\bfat\s+and\s+skins?\b", re.I)

_MEAT_HEAD_RE = re.compile(
    r"^(?:beef|veal|lamb|mutton|pork|chicken|turkey|duck|goat|game)\b", re.I)
_ORIGIN_CHUNK_RE = re.compile(
    r"^(?:new\s+zealand|australian|american|canadian|imported|domestic|"
    r"u\.?s\.?a?\.?)$", re.I)

# A comma chunk that names a legal corporate entity is a manufacturer, not a
# food. FDC files beverages as "Beverages, <COMPANY>, <PRODUCT>, <food>", and
# when the brand list strips the trading name the entity form is left standing
# on its own: "The COCA-COLA company, DASANI, water, bottled" had become the
# canon 'the company, dasani'. Matched on the whole chunk so a food can never
# be hit - no food name ends in "gmbh" or "inc".
_CORP_ENTITY_RE = re.compile(
    r"^(?:the\s+)?[\w\s.&'\u2019/-]*?\b(?:compan(?:y|ies)|gmbh|co\.\s*kg|"
    r"inc\.?|corp\.?|corporation|ltd\.?|limited|llc|l\.l\.c\.?|"
    r"b\.\s*v\.?|n\.\s*v\.?|plc|pty|s\.\s*a\.\s*s\.?)\b[\w\s.&'\u2019/-]*$",
    re.I)

# NOTE: a bare "brand" marker ("Soybean flour, Arrowhead Mills brand",
# "Mars bar and own brand equivalents", "Farina, assorted brands including
# CREAM OF WHEAT") is deliberately NOT stripped by rule. Round 4 tried it and
# it made things worse: _CANON_OVERRIDES already carries fourteen curated
# entries keyed on the marker being present, and removing the word before the
# override lookup invalidated every one of them - "arrowhead mills brand"
# became the meaningless "arrowhead mill" instead of resolving to "soybean
# flour". Stripping the marker cannot generalise anyway, because the trading
# name it leaves behind still needs a curated target. Keep the marker so the
# override key matches, and add an override for each new brand that appears.

# Country of origin stated as provenance rather than as part of the food name.
# "Sweet cherries, imported from the U.S.A., raw" is the same cherry as any
# other sweet cherry per 100 g. Distinct from _ORIGIN_CHUNK_RE, which handles
# the bare adjective form ("lamb, new zealand") and is gated to meat heads;
# an explicit "imported from" is unambiguous, so it needs no gate.
_PROVENANCE_RE = re.compile(
    r"\s*,?\s*\b(?:imported|exported|shipped|sourced|grown|produced|harvested)"
    r"\s+(?:from|in)\b[^,]*", re.I)

# Database bookkeeping that leaked into the food name. CIQUAL marks superseded
# rows by appending "-> ARCHIVE" to the description itself.
_DB_MARKER_RE = re.compile(
    r"\s*-+>\s*archive\b.*$|\s*\b(?:see\s+(?:also|entry)|formerly|obsolete|"
    r"do\s+not\s+use|deleted)\b[^,]*", re.I)

# BioFoodComp ships rows whose species field is literally "N/A". They are real
# measurements of an unidentified plant, so they are kept rather than dropped,
# but a canon named 'n/a, leaf' reads as a food code in a published index.
# "ready to feed" and "ready-to-feed" are the same claim written two ways, and
# the databases use both: 91 canons carry the hyphens, 6 do not, so the split is
# pure spelling. Normalising to the majority form merges each stray back into
# its twin - "Infant formula, MEAD JOHNSON, Enfamil 24, ready to feed" stopped
# reaching 'infant formula, ready-to-feed' once the code chunk before it was
# dropped, because nothing else in the pipeline hyphenates the phrase.
# Applied to the finished canon, not to the description (see the note at
# the call site).
_READY_TO_RE = re.compile(
    r"\bready\s+to\s+(eat|drink|feed|serve|use|cook|bake|heat)\b", re.I)

_NA_HEAD_RE = re.compile(r"^n\s*/\s*a\b", re.I)
# The same placeholder away from the head is just noise and must go, not be
# renamed: FDC files "Mission Figs, Dried, Region 1, n/a, No, Carotenoids",
# and once the region code is dropped as a cultivar code the "n/a" is next
# in line to become the qualifier ("mission fig, n/a, dried").
# The "<medium> pack" wording has to be removed once it has been read, exactly
# as the salt and sugar phrasings are. _detect_pack declines to label a DRAINED
# food - the medium has been poured off - but the words stayed behind, so
# "Pineapple, canned, juice pack, drained" kept 'juice' as a chunk and the
# material rule then folded it into 'pineapple juice'.
_PACK_PHRASE_RE = re.compile(
    r"\s*,?\s*\b(?:extra\s+)?(?:heavy|light|medium)?\s*"
    r"(?:syrup|juice|water|brine|oil)\s+pack(?:ed)?\b", re.I)
# The fat, fortification and pressing-grade wordings have to go the same way
# the salt and sugar phrasings do. _append_state declines to add a label when
# the token it matched is still in the name - a guard that stops "chili powder"
# becoming "chili powder, dried" - so leaving the wording in place printed the
# SOURCE spelling instead of the house label, and the same food ended up under
# both: 'cheese, low fat' beside 'cheese, low-fat' (143 members), 'milk,
# skimmed' beside 'milk, fat-free', 'soy flour, defatted' beside
# 'soy flour, fat-free'.
def _marker_phrase(patterns) -> re.Pattern:
    return re.compile("|".join(rf"\s*,?\s*(?:{rx.pattern})"
                               for _, rx in patterns), re.I)


# "baked without fat" repeats a low-fat claim that has already been read, but
# where NO fat level was read it is the only thing saying the food was fried
# dry: "Atlantic herring, without added fat, fried" is not fried herring.
# The trim wording goes the same way the fat wording does: FDC writes "lean
# only" and "separable fat", and leaving the phrase in place made _append_state
# decline the label, so 372 rows kept 'chicken, lean only' instead of the house
# 'lean'.
_TRIM_PHRASE_RE = _marker_phrase(_TRIM_PATTERNS)
_FAT_PHRASE_RE = re.compile(
    _marker_phrase(_FAT_PATTERNS).pattern
    + r"|\s*,?\s*\bwithout\s+(?:added\s+)?fats?\b"
    # the approximation word goes WITH the figure it qualifies, or it is left
    # behind as a chunk: 'tomme cheese, around, 13% fat'
    + r"|\s*,?\s*\b(?:around|about|approx\.?|approximately|env\.?)(?=\s*,?\s*[\d<>])", re.I)

# Frida splits the figure from its unit across two chunks - "Milk, whole, 3.5,
# (UHT), % fat" - so the fat axis saw no number and "% fat" survived as a chunk
# of the name. Put the figure back on its unit before anything reads the string.
_SPLIT_PCT_RE = re.compile(
    r",\s*(\d+(?:[.,]\d+)?)\s*,\s*((?:\([^)]*\)\s*,\s*)?)%\s*fat\b", re.I)
_FORT_PHRASE_RE  = _marker_phrase(_FORTIFY_PATTERNS)
_GRADE_PHRASE_RE = _marker_phrase(_OIL_GRADE_PATTERNS)
_POWDER_CHUNK_RE = re.compile(r"\s*,\s*powdere?d?\s*(?=,|$)", re.I)
# The UNfortified form is the unmarked one, so its wording is dropped outright.
# FDC pairs the two claims - "evaporated, with added vitamin D and without added
# vitamin A" - and taking only the positive half left the negative half standing
# as 'milk, evaporated and without added vitamin'.
_UNFORTIFIED_PHRASE_RE = re.compile(
    # the denial and whatever it denies: "Salt, not fortified with iodine" kept
    # "not ... with iodine" once the bare word was stripped out from under it
    r"\s*,?\s*\b(?:not|un|non[-\s]?)\s*fortified"
    r"(?:\s+with\s+(?:vitamins?e?(?:\s+[a-k]\d*)?|minerals?|calcium|ca|iron|"
    r"folate|folic\s+acid|iodine|zinc|fibre|fiber)\b"
    r"(?:\s*(?:,|and/or|and|or|&)\s*(?:vitamins?e?(?:\s+[a-k]\d*)?|minerals?|calcium|ca|"
    r"iron|folate|folic\s+acid|iodine|zinc|fibre|fiber)\b)*)?|"
    r"\s*,?\s*\b(?:and\s+)?(?:without|no)\s+added\s+"
    r"(?:vitamins?(?:\s+[a-k]\d*\b)?|minerals?|calcium|iron|folate|"
    r"folic\s+acid|iodine|zinc|fibre|fiber)"
    r"(?:\s*(?:,|and/or|and|or|&)\s*(?:added\s+)?"
    r"(?:vitamins?(?:\s+[a-k]\d*\b)?|minerals?|calcium|iron|folate|"
    r"folic\s+acid|iodine|zinc|fibre|fiber))*", re.I)
_NA_CHUNK_RE = re.compile(r"^n\s*/\s*a\.?$", re.I)
# Where a chilled product sits in the shop is not what it is. CIQUAL and FDC
# both file the storage form as a chunk of its own ("Soymilk, refrigerated",
# "Chocolate mousse, refrigerated"), which split the chilled product from its
# ambient twin though per 100 g they are the same food. Matched as a WHOLE
# chunk only, so "buttermilk biscuit, refrigerated dough" - where the raw
# dough IS the food and is not a baked biscuit - keeps its name.
# FDC writes some compound names with the qualifier AFTER the comma ("Yogurt,
# Greek, strawberry, non-fat, CHOBANI"), while every other database writes it in
# front ("Greek yogurt, plain, non-fat"). Only the head chunk and one more reach
# the canon, so the comma form pushed the flavour out of the name and 108 rows
# of Greek strawberry yoghurt landed on 'yoghurt, fat-free' beside the plain
# ones, handing plain fat-free yoghurt the strawberry sugar. Moving the
# qualifier to the front puts both databases on the same shape. Curated, not a
# blanket rule: it is only correct where the qualifier names a STYLE of the
# food, not a variety of it.
_FRONTED_QUALIFIER_RE = re.compile(
    r"^(yogh?o?urt)\s*,\s*(greek|icelandic|skyr|bulgarian|swiss)\b", re.I)

# A taxonomic or menu CATEGORY in front of a food is not the food. FDC and CNF
# file whole chapters this way - "Crustaceans, crab, alaska king", "Fish,
# salmon, atlantic", "Game meat, bison, chuck" - and the head then occupies the
# name while the species sits behind it. Six heads qualify, on the evidence that
# each fronts dozens of DIFFERENT foods and the food behind it usually stands as
# a canon on its own: fish (234 canons), mollusk (57), spice (57), grain (47),
# game meat (30) and crustacean (29). It also settles what "game" means, which
# FDC is loose about: its game chapter holds farmed bison, beefalo, goat, horse
# and "rabbit, domesticated".
#
# "deli-meat" is deliberately NOT here, though it has the same shape: it carries
# a processing claim the product name does not always repeat, and dropping it
# would fold "deli-meat, chicken" into the fresh-chicken canon.
#
# The lookahead is the guard. Where the chunk after the head is a state rather
# than a species the head IS the food, and stripping it would leave a name with
# no noun in it: "fish, dried", "fish, roe", "game meat, dried, salted".
_CATEGORY_HEAD_RE = re.compile(
    r"^(?:fish(?:es)?|spices?|grains?|mollus[ck]s?|crustaceans?|game\s+meats?|nuts?)\s*,\s*"
    # \s* inside the lookahead as well as before it: without it the engine
    # backtracks the separator to zero width, the guard then sees a leading
    # space and fails to match, and the head is stripped after all.
    r"(?!\s*(?:dried|dry|salted|unsalted|smoked|canned|cooked|raw|boiled|fried|frozen|"
    r"fresh|whole|lean|fillets?|flesh|roe|liver|average|mixed\s+species|steaks?|"
    r"in\s+oil|in\s+water|in\s+brine|"
    # ... and where the chunk names an INGREDIENT rather than a species. Fineli
    # writes one dish as "Fish, Cooking Cream (19% Fat), Oven-Baked"; stripping
    # the head there left it colliding with the three real cooking-cream rows.
    r"cooking\s+cream|"
    # "nut" joins the list on the same evidence (31 canons: 'nut, acorn',
    # 'nut, chestnut', 'nut, walnut, dried', and 'nut, pecan, salted' already
    # colliding with 'pecan nut, salted'). Its own three exceptions are
    # aggregate or dish entries where "nut" is an ingredient, not a chapter:
    # "Nuts, formulated, wheat-based", "Nuts, simulated product" and McCance's
    # "Nut, mushroom and rice roast".
    r"formulated|simulated|mushroom\s+and\s+rice|"
    # ... and where the chunk is only HALF a compound nut name. Frida and AFCD
    # write "Nut, brazil", "Nut, pine", "Nut, pea", "Nut, coco"; taking the head
    # off leaves 'brazil', 'pine', 'pea' and 'coco', none of which is the food.
    r"brazil|pine|peas?|coco|hazel|beech|chest|ground|monkey|kola|cola|betel|"
    r"areca|tiger|water|butter)\b)", re.I)

# "ready to bake / fry / cook" says the food still has to be cooked, which is the
# unmarked raw state - and in FDC's phrasing the claim takes the single qualifier
# slot a canon has, so "Tortillas, ready-to-bake or -fry, corn" became
# 'tortilla, ready-to-bake or -fry' and lost the corn-versus-flour distinction
# that actually changes the composition. Deliberately NOT applied to
# ready-to-eat / -drink / -feed / -serve / -use, where the claim marks the food
# as being in its final form rather than a concentrate or a dry mix.
_READY_RAW_RE = re.compile(
    r"\s*,?\s*\bready[-\s]?to[-\s]?(?:bake|fry|cook|heat|roast|grill)\b"
    r"(?:\s*(?:or|/|,)\s*-?\s*(?:bake|fry|cook|heat|roast|grill)\b)*", re.I)

# BioFoodComp writes one part three ways - "flesh" (512 rows), "pulp" (147) and
# "fruit flesh" (117). _detect_part already maps all three to the "flesh" label,
# but where the CHUNK survives the label is declined and the source's spelling
# stands: 'coconut, pulp' sat beside 'coconut, flesh', and seven canons carried
# both words at once ('baobab fruit/monkey bread, pulp, flesh'). Whole chunk
# only - no canon in the corpus uses "pulp" in the purée sense, and a bare word
# inside a longer chunk is not the part.
_PULP_CHUNK_RE = re.compile(r"^(?:fruit\s+flesh|pulp)$", re.I)

_STORAGE_CHUNK_RE = re.compile(
    r"^(?:refrigerated|chilled|shelf[-\s]?stable|ambient)$", re.I)

# A qualifier that states the UNMARKED case says nothing, and every one of them
# split a food from itself. Measured on the shipped index: 36 foods carried both
# "X" and "X, plain" ('almond milk', 'bagel', 'butter', 'tofu', 'yoghurt'), and
# 20 carried both "X" and "X, mixed species" - FDC's way of writing that the
# species behind a genus value was not recorded ('shrimp', 'squid', 'trout').
# Whole chunk only: "beans, baked, plain or vegetarian" and "biscuit,
# plain/buttermilk" are alternations, and Phenol-Explorer's "Common cabbage" and
# "Common octopus" are species names, which is why "common" is not on this list.
_NULL_QUALIFIER_RE = re.compile(
    # "regular" and "average" are already struck by _FILLER_CHUNK_RE; "all" is
    # not, and FDC's "Alcoholic beverage, wine, table, all" carried it into the
    # name.
    r"^(?:plain|mixed\s+species|ordinary|all)$"
    # A chunk that says the axis was not recorded says nothing at all. CIQUAL
    # writes it out - "Grenadier, from any fishing spot", "Ice cream or sorbet
    # or ice pop, any flavour", "Fruit compote and similar, any fruit" - and the
    # phrase became part of the name. No food name begins with "any".
    r"|^(?:from\s+any|any)\b", re.I)
# ...except on chocolate, where British and Nordic sources use "plain" for DARK:
# Fineli spells it out ("Chocolate, Plain, Dark Chocolate"), McCance and Frida
# leave it implied. Dropping the word there would pour dark chocolate into the
# milk-chocolate canon.
_PLAIN_IS_DARK_RE = re.compile(r"^chocolates?\b", re.I)
# NOTE on "average": the rules already drop it (_FILLER_CHUNK_RE), so every
# 'beef, average' / 'lettuce, average' / 'potato, average' in the index came in
# through a curated table VALUE, which is never re-canonicalised. Those are
# normalised in the table itself rather than by a rule that could not fire.

# A frying fat named only as "vegetable oil" names nothing - it is whatever the
# processor happened to use - and the phrase survived as a chunk, splitting the
# food from its own canon: FDC's "potato, french fried in vegetable oil" sat at
# 'french fry, in vegetable oil, fried' beside Fineli's 'french fry, fried', and
# 15 more did the same. A NAMED fat stays: what is fried in rapeseed and what is
# fried in olive oil absorb different fatty acids, and 49 canons state one.
_GENERIC_OIL_RE = re.compile(
    r"\s*,?\s*\bin\s+(?:vegetable|blended|cooking|frying|unspecified)\s+"
    r"(?:oils?|fats?|shortenings?)\b", re.I)

# The fat a food was cooked IN, where the source names it. McCance and Fineli
# state it on 176 rows and the corpus was reading it by accident: "Doughnut,
# fried in rapeseed oil" kept the phrase while "Aubergine, fried in rapeseed
# oil" lost it to the two-chunk strip head, so 19 canons held two or three
# different fats at once - 'potato chip, homemade, fried' held corn, rapeseed and
# sunflower at once. What a deep-fried food absorbs is the frying fat, 5-15 g of
# it per 100 g, and corn and sunflower oil carry three times rapeseed's linoleic
# acid while rapeseed carries the alpha-linolenic.
#
# Read as an axis so the label is printed in one spelling wherever the source
# states it. The unnamed fats never reach here - _GENERIC_OIL_RE has already
# struck "in vegetable oil" and its kind, which name nothing.
_COOK_FAT_RE = re.compile(
    r"\b(?:deep[-\s]|pan[-\s]|stir[-\s]|shallow[-\s]|oven[-\s])?"
    r"(?:fried|roasted|saut[e\u00e9]ed|sauteed|cooked|baked|browned)\s+"
    r"(?:in|with)\s+([a-z]+)\s+(?:oils?|fats?)\b", re.I)


def _detect_cook_fat(d: str) -> tuple[str, str]:
    """Return ("in <fat> oil", matched phrase), or ("",""). See above."""
    hit = _COOK_FAT_RE.search(d)
    if not hit:
        return "", ""
    return f"in {hit.group(1).lower()} oil", hit.group(0)


_TRADEMARK_RE = re.compile(r"[\u00ae\u2122\u2117]")

# CIQUAL writes the plural as an option: "with sugar(s)", "fortified with
# vitamin(s) and/or mineral(s)". Normalised before any probe is taken.
_OPTIONAL_PLURAL_RE = re.compile(r"(?<=[a-z])\(s\)", re.I)

# NEVO names its dairy substitutes by what they REPLACE and what they are made
# of, in one long sentence: "Plant-based alternative to Gouda cheese based on
# coconut oil fortified w Ca and Vit B12". The whole sentence became the canon.
# Rewritten to the house shape - the food first, the base in the qualifier slot -
# which also keeps the base: "plant-based gouda cheese, coconut oil". The two
# halves stay in ONE name rather than two chunks, because the two-chunk rule
# would then drop the base and merge the oat cream into the soy one.
_PLANT_ALT_BASE_RE = re.compile(
    r"\bplant[-\s]?based\s+alternatives?\s+to\s+(.+?)\s+based\s+on\s+", re.I)
_PLANT_ALT_RE = re.compile(r"\bplant[-\s]?based\s+alternatives?\s+to\s+", re.I)
# "Soy milk, non-dairy alternative to milk" says twice what "soy milk" says once.
_DAIRY_ALT_RE = re.compile(
    r"\s*,?\s*\bnon[-\s]?dairy\s+alternatives?\s+to\s+milk\b", re.I)

# ---------------------------------------------------------------------------
# ALTERNATIONS
#
# 488 canons carry an "A or B" in their name. They are not one class. Where the
# two sides are DIFFERENT foods the "or" is the source telling you it does not
# know which one was measured - "Sausage, beef or pork meat", "Milk, sheep or
# goat" - and that honesty is worth keeping: collapsing to the first side would
# file a pork sausage under beef and pull its iron toward beef's.
#
# Where the two sides are the SAME food the "or" is a source's house style, and
# it split that food across as many canons as there are spellings:
# 'yoghurt or fermented milk, plain' sat beside 'yoghurt, plain', and FDC's
# "Macaroni or noodles with cheese" beside plain macaroni cheese. Only those are
# collapsed, and only from this hand-checked list - there is no rule that can
# tell "chicken or turkey" from "grissini or bread stick" without knowing what
# the words mean.
#
# The replacement is the side more canons already spell, counted rather than
# preferred, exactly as _SYNONYM_PAIRS is built.
_ALTERNATION_PHRASES = [
    # dairy: NEVO/CIQUAL class labels. 24 + 12 canons.
    (r"yogh?urts?\s*or\s+fermented\s+milks?",            "yoghurt"),
    (r"fermented\s+milks?\s+or\s+yogh?urts?",            "yoghurt"),
    (r"quarks?\s+or\s+dairy\s+special(?:i)?t(?:y|ies)",  "quark"),
    # FDC's survey wording; the user's spot-check that started this round
    (r"macaronis?\s+or\s+noodles?",                      "macaroni"),
    # soft drinks: three spellings of one category
    (r"carbonated\s+beverages?\s+or\s+fruit\s+soft\s+drinks?", "soft drink"),
    (r"carbonated\s+beverages?\s+or\s+soft\s+drinks?",   "soft drink"),
    # a crepe IS a thin pancake; CIQUAL spells the buckwheat one both ways
    (r"crepes?\s+or\s+buckwheat\s+(?:pan)?cakes?",       "crepe"),
    (r"crepes?\s+or\s+buckwheat\s+crepes?",              "crepe"),
    (r"pikelets?\s+or\s+pancakes?",                      "pancake"),
    # Sardina pilchardus under three names
    (r"(?:european\s+)?pilchards?\s+or\s+sardines?",     "sardine"),
    (r"sardines?\s+or\s+(?:european\s+)?pilchards?",     "sardine"),
    # Corylus avellana under two
    (r"filberts?\s+or\s+hazelnuts?",                     "hazelnut"),
    (r"hazelnuts?\s+or\s+filberts?",                     "hazelnut"),
    # Spondias dulcis under three, two of them French
    (r"prune\s+de\s+cyth[eè]re\s+or\s+pomme\s+cyth[eè]re?\s+or\s+golden\s+apple", "golden apple"),
    (r"pomme\s+cyth[eè]re?\s+or\s+golden\s+apple",       "golden apple"),
    # same food, different national name
    (r"chapatis?\s+or\s+rotis?",                         "chapati"),
    (r"grissinis?\s+or\s+bread\s+sticks?",               "bread stick"),
    (r"marzipans?\s+or\s+almond\s+pastes?",              "marzipan"),
    (r"phyllos?\s+or\s+filo\s+doughs?",                  "filo dough"),
    (r"nuoc\s+mam\s+sauces?\s+or\s+fish\s+sauces?",      "fish sauce"),
    (r"sakes?\s+or\s+rice\s+wines?",                     "sake"),
    (r"water\s+kefirs?\s+or\s+tibicos?",                 "water kefir"),
    (r"kombus?\s+or\s+japanese\s+kelps?",                "kombu"),
    (r"head[-\s]?cheese\s+p[aâ]t[ée]s?\s+or\s+brawns?",  "brawn"),
    (r"egg\s+rolls?\s+or\s+nems?",                       "egg roll"),
    (r"samoo?sas?\s+or\s+samoo?sas?\s+or\s+filled\s+filo\s+pastry", "samosa"),
    (r"sponge\s+fingers?\s+or\s+lady\s+fingers?",        "sponge finger"),
    (r"vine\s+leaf\s+stuffed\s+wi?th?\s+rice\s+or\s+dolmas?", "dolma"),
    (r"meatballs?\s+or\s+rissoles?",                     "meatball"),
    (r"veggie\s+burgers?\s+or\s+soyburgers?",            "veggie burger"),
    (r"black\s+puddings?\s+or\s+sausages?",              "black pudding"),
    (r"baking\s+powders?\s+or\s+raising\s+agents?",      "baking powder"),
    (r"masetholes?\s+or\s+white\s+milkwoods?",           "white milkwood"),
    (r"rosebay\s+willow\s+herbs?\s+or\s+fireweeds?",     "fireweed"),
    (r"star\s+apples?\s+or\s+milk\s+fruits?",            "star apple"),
    (r"norway\s+lobsters?\s+or\s+scampis?",              "norway lobster"),
    (r"chinese\s+or\s+japanese\s+artichokes?",           "chinese artichoke"),
    # Stachys, Agaricus, Boletus, Brassica: one organism, two common names
    (r"caesar'?s?\s+mushrooms?\s+or\s+royal\s+agarics?", "caesar's mushroom"),
    (r"ceps?\s+or\s+boletus\s+mushrooms?",               "cep"),
    (r"button\s+mushrooms?\s+or\s+cultivated\s+mushrooms?", "button mushroom"),
    (r"romanesco\s+cauliflowers?\s+or\s+romanesco\s+broccolis?", "romanesco"),
    (r"fruit\s+cocktails?\s+or\s+salads?",               "fruit cocktail"),
    (r"pop-?\s?corn\s+or\s+(?:air-)?popped\s+corn",      "popcorn"),
    (r"shrimps?\s+or\s+prawns?",                         "shrimp"),
    # the spelling of ONE word differs, not the food
    (r"vegetarian\s+meat\s?loaf\s+or\s+patty",           "vegetarian meatloaf"),
    # form, not food: "Barley flour or meal" is one milled product
    (r"flours?\s+or\s+meals?",                           "flour"),
    (r"puffed\s+or\s+extruded",                          "puffed"),
    # FDC's own category label, not two foods: "pasteurized process cheese food
    # or product" is one processed cheese.
    (r"cheese\s+foods?\s+or\s+products?",              "cheese food"),
]
_ALTERNATION_RE = [(re.compile(r"\b" + pat + r"\b", re.I), rep)
                   for pat, rep in _ALTERNATION_PHRASES]

# Retail formats offered as a choice. Neither side is a food and the pair took
# the qualifier slot: "Butter, stick or tub" and "Cheese, square or brick".
_FORM_ALTERNATION_RE = re.compile(
    r"^(?:ice\s+cream\s+)?(?:bar|stick|tub|square|brick|block|dice|strip|slice|"
    r"shredded|diced|whole|crushed|granulated|lump)s?"
    r"\s+or\s+"
    r"(?:bar|stick|tub|square|brick|block|dice|strip|slice|"
    r"shredded|diced|whole|crushed|granulated|lump)s?(?:\s+grain)?$", re.I)

# Several sources restate the head in the qualifier slot - Fineli's "Pea Stew,
# Dried Peas", STFCJ's "Bun with filling, fried bun", NEVO's "Spinach Soup,
# Spinach". The chunk adds no word the head does not already carry, and it took
# the ONE qualifier slot a canon has, so the real qualifier behind it was lost.
#
# Applied to the SOURCE chunks only, never to a label an axis appended. That is
# the whole guard: "Water shield, young leaves, bottled in water" ends up as
# 'water shield, leaf, in water', where both trailing chunks look redundant
# against a head that happens to contain the word "water" - and both are real.
_RESTATES_STOP = frozenset({"in", "of", "with", "from", "the", "a", "an",
                            "and", "or", "made", "type", "style"})
# A chunk that OPENS with a preposition is a phrase about the food, not another
# name for it: "Vanilla ice cream, without cream" must not become "without".
_CHUNK_PREP_RE = re.compile(r"^(?:with|without|in|from|for|on|to|containing)\b", re.I)


def _drops_restatement(chunks: list[str]) -> list[str]:
    """Drop chunks whose informative words all appear in the head already."""
    # compared on the singular: the plural fold runs at the very end, so at this
    # point Fineli's "Pea Stew, Dried Peas" still has "peas" against a "pea" head
    def _words(x: str) -> set:
        return {_singularize(w) for w in
                re.findall(r"[a-z0-9\u00c0-\u024f]+", x.lower())} - _RESTATES_STOP
    head = _words(chunks[0])
    if not head:
        return chunks
    # Only for a head of two words or more. Where the head is the class noun on
    # its own - "Cake, cheese cake", "Sausage, liver sausage" - it is the HEAD
    # that is redundant, not the chunk, and the fold to 'cheese cake' belongs to
    # _SELF_NAMING_TAIL_RE; stripping here would give 'cake, cheese' instead.
    last = re.findall(r"[a-z0-9\u00c0-\u024f]+", chunks[0].lower())
    last = last[-1] if len(last) > 1 else ""
    out = [chunks[0]]
    for c in chunks[1:]:
        w = _words(c)
        if w and w <= head:
            continue
        # The chunk names the head's own class again: "Salad dressing, french
        # dressing", "Breakfast cereal, corn cereal", "Multigrain bread, oat
        # bread". Saying it twice is what split 'salad dressing, french' from
        # 'salad dressing, french dressing'. The repeat goes, the qualifier
        # stays. Run here, on the SOURCE chunks, so no axis label is at risk:
        # "apricot, in juice" and "spreadable fat, 41% fat" are appended later
        # and would each lose their noun to this.
        cw = c.split()
        if (len(cw) >= 2 and last and cw[-1].lower().rstrip(",") == last
                and not _CHUNK_PREP_RE.match(c)):
            c = " ".join(cw[:-1])
        out.append(c)
    return out

# USDA files every ordinary chicken under its retail class - "Chicken, broilers
# or fryers, breast, meat only, raw" - and CNF writes the same bird as
# "Chicken, broiler, back". The term names a market size, not a food, and it
# split one bird across two families: 38 canons spelling it "broilers or fryer"
# stood beside 41 spelling it "broiler". Only the DEFAULT class is dropped. A
# stewing hen is a spent layer carrying several times the fat of a young fryer,
# and "roasting", "stewing" and "capon" stay on the name for that reason - as
# does turkey's "fryer-roaster", which is the young bird rather than the norm.
_POULTRY_CLASS_RE = re.compile(
    r"\b(chickens?|turkeys?)\s*,?\s*"
    r"(?:broilers?\s+or\s+fryers?|broilers?|fryers?(?![-\s]roaster)|"
    r"all\s+classes)\b", re.I)

# Livsmedelsverket records how the item is SOLD, not what it is: "Vegetarian
# sausage w/ soy protein chilled or frozen product". _PREP_RE takes "chilled"
# and "frozen" and leaves "or product" standing, which became a chunk of the
# name on 25 canons - 'plant-based nugget soy protein or product'. Stripped as
# a whole phrase, and only behind chilled/frozen, so NEVO's bare "product"
# (which its curated overrides read as the frozen ready meal) is untouched.
# CNF writes "light meat only" where FDC writes "light meat, meat only"; the
# whole-chunk test could not see the claim and four skinless turkey rows sat in
# the with-skin canon.
_LIGHT_MEAT_ONLY_RE = re.compile(r"\b(light|dark)\s+meat\s+only\b", re.I)

_POULTRY_ALONE_RE = re.compile(
    r",\s*(?:broilers?\s+or\s+fryers?|broilers?|fryers?(?![-\s]roaster))\s*(?=,|$)", re.I)

_SUPPLY_FORM_RE = re.compile(
    r"\s*,?\s*\b(?:chilled|frozen|refrigerated)\s*"
    r"(?:or\s+(?:chilled|frozen|refrigerated)\s*)?products?\b", re.I)

# NEVO and Livsmedelsverket write no punctuation at all: "Beans broad raw",
# "Nuts macadamia unsalted", "Crackers cream", "Gherkins pickled wo sugar". The
# canonicalizer splits on commas, so the whole string stayed one chunk and 154
# canons kept a plural head welded to their modifier - 'beans broad' sitting
# beside the 103-member 'broad bean'. Restoring the comma the source omitted is
# enough; every rule downstream then works as it does on any other source.
#
# The head has to be the FIRST word and genuinely plural, and what follows it
# must be a word rather than a connective or punctuation - "Nuts and seeds",
# "Beans, snap" and "Lentils /glass" are all left alone. _NOT_PLURAL carries the
# names that only look plural (Brussels, Maroilles, Causses, bitters, sports).
_PLURAL_HEAD_RE = re.compile(r"^([A-Za-z]{4,})\s+(?=[a-z])", re.I)
_PLURAL_HEAD_SKIP_RE = re.compile(
    r"^(?:with|without|and|or|in|of|for|from|on|to|the|a|an|av|w|wo|"
    # NEVO closes a processed entry with "product" ("Strawberries product" is
    # the frozen one), and several curated entries are keyed on that wording
    r"products?)\b", re.I)


def _split_plural_head(d: str) -> str:
    m = _PLURAL_HEAD_RE.match(d)
    if not m:
        return d
    head = m.group(1)
    single = _singularize(head.lower())
    if single == head.lower() or len(single) < 3:
        return d
    rest = d[m.end():]
    if _PLURAL_HEAD_SKIP_RE.match(rest):
        return d
    # A DISH is not a modifier of the fruit in front of it. CIQUAL writes
    # "Strawberries tart, from bakery"; inserting the comma exposed "tart" to
    # the strawberry strip-head, which drops everything behind it, and the tart
    # landed in the fruit canon.
    if _MATERIAL_TAIL_BLOCK_RE.match(rest):
        return d
    return f"{head}, {rest}"

# One spelling per compound. The sources hyphenate these differently and the
# index carried both: 'bread, gluten free' beside 'bread, gluten-free',
# 'soy flour, full fat' beside 'full-fat', 'barley flour, whole grain' beside
# 'barley flour, wholegrain' (107 canons spell it closed against 49 open).
_COMPOUND_SPELLINGS = [
    (re.compile(r"\bwhole[-\s]?grains?\b", re.I), "wholegrain"),
    (re.compile(r"\bwhole[-\s]?meal\b", re.I), "wholemeal"),
    (re.compile(r"\bgluten[-\s]free\b", re.I), "gluten-free"),
    (re.compile(r"\bdairy[-\s]free\b", re.I), "dairy-free"),
    (re.compile(r"\blactose[-\s]free\b", re.I), "lactose-free"),
    (re.compile(r"\balcohol[-\s]free\b", re.I), "alcohol-free"),
    (re.compile(r"\bcaffeine[-\s]free\b", re.I), "caffeine-free"),
    (re.compile(r"\bfull[-\s]fat\b", re.I), "full-fat"),
    (re.compile(r"\bhalf[-\s]fat\b", re.I), "half-fat"),
    (re.compile(r"\bwhole[-\s]wheat\b", re.I), "whole wheat"),
    (re.compile(r"\bmulti[-\s]grain\b", re.I), "multigrain"),
    # 264 canons spell it closed and 54 open, for the same claim
    (re.compile(r"\bhome[-\s]made\b", re.I), "homemade"),
    # "flesh only" is the part axis's own label with the word "only" still on
    # it: _detect_part trims it, but where the chunk SURVIVES the label is
    # declined and the source spelling stands. 16 foods carried both spellings -
    # 'chicken, flesh' beside 'chicken, flesh only', and the same for duck,
    # turkey, rabbit, pheasant, coconut, durian, lychee and eight more.
    (re.compile(r"\bflesh\s+only\b", re.I), "flesh"),
    # A table wine is an ordinary wine - the word separates it from fortified
    # and dessert wines, which have names of their own here - so 'wine, table,
    # red' and 'table wine' both belong on 'wine'.
    # three spellings for one biscuit: 'breadstick' (4 canons), 'bread stick'
    # (2) and 'bread, stick' (1)
    (re.compile(r"\bbread\s*,\s*stick\b|\bbread\s+stick\b", re.I), "breadstick"),
    (re.compile(r"\btable\s+wine\b", re.I), "wine"),
    (re.compile(r"\bwine\s*,\s*table\b", re.I), "wine"),
    # ...and the provenance axis below is spelled three ways in the sources
    (re.compile(r"\bwild[-\s](?:caught|harvested)\b", re.I), "wild"),
    (re.compile(r"\bready[-\s]to[-\s]", re.I), "ready-to-"),
]


# The cut and organ axes append their label behind the animal, so every source
# that writes the same thing WITHOUT a comma - "Beef round", "Pork shoulder",
# "Beef liver" - or with the parts the other way round - "Liver, pork",
# "Tongue, beef" - ended up under a second canon. 20 pairs, among them
# 'beef, round' (122) against 'beef round', and 'beef liver' (12) against
# 'beef, liver' (7).
_MEAT_NAMES = (r"beef|pork|lamb|veal|mutton|goat|venison|bison|buffalo|horse|"
               r"rabbit|chicken|turkey|duck|goose|ox|pig|calf|sheep")
# McCance writes the animal, the others the meat: an ox liver is a beef liver
# and a calf's is a veal one.
_ANIMAL_TO_MEAT = {"ox": "beef", "pig": "pork", "calf": "veal", "sheep": "lamb"}

# One grind, two words. FDC and CNF write "Beef, ground"; McCance, NEVO and
# Livsmedelsverket write "Beef mince". 96 canons carry the label as "ground"
# against 19 spelling it "mince", so the mince forms fold in. Applied to a
# WHOLE chunk only and gated on the animal, because a mince PIE and mincemeat
# are fruit, and "minced beef ball" is the name of a dish.
_MINCE_ANIMAL = (r"beef|pork|ham|bacon|lamb|veal|mutton|chicken|turkey|duck|"
                 r"goose|venison|elk|reindeer|moose|game|rabbit|horse")
_MINCE_RE = re.compile(rf"\b({_MINCE_ANIMAL})\s*,?\s*minced?\s*(?=,|$)", re.I)
# ...and the fronted spelling, which curated values still carry:
# "ground turkey, patty, 7% fat" beside "turkey, patty, ground, 7% fat".
_GROUND_HEAD_RE = re.compile(rf"^ground\s+({_MINCE_ANIMAL})\b", re.I)

# Bread reads type-first in this corpus and by a wide margin - 35 canons spell
# it "rye bread" against 5 for "bread, rye", 74 "multigrain bread" against 2,
# 46 "wheat bread" against 12 - but a curated entry pushed the other way and
# left both forms standing. Only the GRAIN types flip: "bread, bagel" and
# "bread, naan" name a bread that is not "<grain> bread", and "bread, wheat
# bun" is a bun, so the modifier has to be the whole chunk and on this list.
_BREAD_TYPE_RE = re.compile(
    r"^bread,\s*(rye|wheat|white|brown|wholegrain|wholemeal|whole\s+wheat|"
    r"multigrain|spelt|sourdough|oat|barley|corn|potato|rice|soda|"
    r"pumpernickel|seeded|unleavened|raisin|flat)\s*(?=,|$)", re.I)


def _meat_name(w: str) -> str:
    w = w.lower()
    return _ANIMAL_TO_MEAT.get(w, w)


_MEAT_PART_JOIN_RE = re.compile(
    rf"^({_MEAT_NAMES})\s+({_CUT_WORDS}|{_ORGAN_WORDS})(?=,|$)", re.I)
_MEAT_PART_SWAP_RE = re.compile(
    rf"^({_ORGAN_WORDS}|{_CUT_WORDS})\s*,\s*({_MEAT_NAMES})(?=,|$)", re.I)


# Twenty cheese types are written BOTH ways in the index, and the head-first
# form wins every one of them: 'cheese, mozzarella' against 'mozzarella cheese,
# from cow's milk', 'cheese, cottage' (146 members) against 'cottage cheese,
# full fat', 'cheese, camembert' (12) against 'camembert cheese, from cow's
# milk'. CIQUAL and the Swiss table write the modifier first, FDC writes the
# head first, and the two never met.
#
# A LIST rather than a pattern, because "X cheese" is not always a cheese type:
# "cauliflower cheese" and "macaroni cheese" are dishes, and turning them into
# 'cheese, cauliflower' would file a gratin as a cheese.
_CHEESE_TYPES = (
    r"brie|camembert|cottage|cream|edam|emmental|feta|firm|fontina|goat|gouda|"
    r"gruy[eè]re|hard|mascarpone|mozzarella|processed|provolone|ricotta|"
    r"roquefort|semi-hard|semi-soft|soft|cheddar|parmesan|halloumi|paneer|"
    r"manchego|gorgonzola|pecorino|romano|asiago|havarti|jarlsberg|raclette|"
    r"tilsiter|appenzeller|sbrinz|comt[eé]|beaufort|abondance|cantal|chaource|"
    r"reblochon|morbier|taleggio|burrata|bocconcini|cotija|panela|sheep|ewe|"
    r"buffalo")
_CHEESE_ORDER_RE = re.compile(rf"^({_CHEESE_TYPES})\s+cheese(?=,|$)", re.I)


# Two OPPOSITE conventions, each following the majority the index has already
# settled on, so that one food does not end up under two spellings.
#
# 1. MATERIAL heads read modifier-first. 'X oil' has 135 canons against 30 for
#    'oil, X'; juice 125/11, flour 190/23, milk 398/88, vinegar 14/4. The head
#    moves behind its modifier - "Oil, coconut" joins the 22-member
#    'coconut oil', "Juice, prune, dried" joins 'prune juice'.
#
#    Only the FIRST chunk moves, and only when it is a plain modifier. Anything
#    that names a state, a composition axis or a grade stays where it is: a
#    "milk, whole" must not become "whole milk" while "milk, skimmed" becomes
#    "skimmed milk", because the pair is then split across two spellings again,
#    and "oil, industrial" / "juice, cocktail" / "flour, all purpose" name no
#    food at all once the head is taken off the front.
_MATERIAL_HEAD_BLOCK = (
    r"whole|skim(?:med)?|semi[-\s]?skimmed|half\s+skimmed|nonfat|non[-\s]?fat|"
    r"lowfat|low[-\s]?fat|reduced\s+fat|fat[-\s]?free|full[-\s]?fat|light|"
    r"dried|dry|raw|fresh|frozen|canned|bottled|cooked|boiled|powdered|powder|"
    r"sweetened|unsweetened|salted|unsalted|enriched|fortified|flavoured|"
    r"flavored|pasteuri[sz]ed|uht|sterili[sz]ed|homogeni[sz]ed|evaporated|"
    r"condensed|cultured|fermented|filled|imitation|plain|bulk|producer|"
    r"lactose[-\s]?free|low[-\s]?lactose|industrial|commercial|cocktail|blend|"
    r"mixed|assorted|average|formulated|simulated|all[-\s]purpose|cooking|"
    r"vegetable|animal|unknown|fat\s+content|with|without|from|in\b|"
    # a FORM, not a modifier: "Milk, fluid, nonfat" is milk, not "fluid milk"
    r"fluid|liquid|solid|semi[-\s]?solid|"
    # a nutrient claim is not a modifier either: Fineli writes
    # "Juice, Unsweetened, Vitamin C, Vitamin D, Probiotics"
    r"vitamins?\s*[a-e]?\d*|minerals?|probiotics?"
)
_MATERIAL_HEAD_RE = re.compile(
    rf"^(oil|juice|flour|vinegar|milk)\s*,\s*(?!\s*(?:{_MATERIAL_HEAD_BLOCK})\b)"
    # the modifier itself: letters, spaces, hyphens and ampersands only, so a
    # percentage ("milk, 1.5% fat"), a code or a grade can never be moved
    r"([a-z][a-z'\-]*(?:[ &][a-z'\-]+){0,2})\s*(?=,|$)", re.I)


# The same head also turns up one chunk LATE, which is Phenol-Explorer's shape:
# "Olive, oil, extra virgin" put the grade in the third chunk, where the
# two-chunk rule dropped it, and the food read as 'olive, oil'. Folding the head
# onto the food frees the slot the grade needs. Only a head immediately after
# the food is folded, so "Salad dressing, french, cottonseed, oil" is untouched.
_MATERIAL_TAIL_RE = re.compile(
    r"^([a-z][a-z'\- ]{2,28}?)\s*,\s*(oils?|juices?|flours?|milk|butter)\s*(?=,|$)",
    re.I)
# Behind a DISH the same chunk is an ingredient or a cooking medium, not the
# material the food is made of. Fineli writes its recipes as ingredient lists -
# "Rice Porridge, Milk, Salt", "Mashed Potato Casserole, Milk, Eggs" - and
# Livsmedelsverket writes "Macaroni boiled w/ oil"; folding the head on gives
# 'rice porridge milk' and 'macaroni oil'. "butter" is not a tail head at all:
# every one of its 11 rows is either butter-flavoured ("Cookies, butter",
# "Icing, butter") or a butter BEAN, and none is a dairy compound.
_MATERIAL_TAIL_BLOCK_RE = re.compile(
    r"\b(?:porridges?|gruels?|casseroles?|stews?|gravy|gravies|purees?|rusks?|"
    r"crackers?|cookies?|biscuits?|croissants?|icings?|sweets?|cakes?|breads?|"
    r"macaroni|pastas?|noodles?|spaghetti|pop\s?corns?|soups?|sauces?|salads?|"
    r"puddings?|pies?|desserts?|beverages?|drinks?|mixe?s?|dishe?s?|"
    r"sandwiches?|pizzas?|omelettes?|risotto|paella|loaf|loaves|bars?|rolls?|"
    r"buns?|pastry|pastries|tarts?|waffles?|pancakes?|crepes?|dumplings?|"
    r"patty|patties|burgers?|nuggets?|cutlets?|schnitzels?|gratins?|quiches?|"
    r"lasagne|moussaka|curry|curries|stir[-\s]?fry|tortillas?|tacos?|wraps?|"
    r"pitta?s?|naans?|chapatis?|rotis?|tostadas?|blinis?)\b", re.I)
# ... and "Beans, butter" is the butter BEAN, a lima, not a dairy compound.
_MATERIAL_TAIL_BEAN_RE = re.compile(r"\bbeans?$", re.I)


# "Chocolate, milk" is a milk chocolate BAR - four sources write it that way -
# and folding it gave 'chocolate milk', which merged 30 g fat/100 g of
# confectionery into NEVO's 3 g/100 g drink. Phenol-Explorer files the drink as
# "Chocolate, milk, beverage", so the beverage word is what tells them apart.
_MILK_BEVERAGE_RE = re.compile(r"\s*,\s*(?:beverage|drink|shake)", re.I)


def _fold_material_tail(d: str) -> str:
    m = _MATERIAL_TAIL_RE.match(d)
    if not m:
        return d
    food, head = m.group(1).strip().lower(), m.group(2).lower().rstrip("s")
    rest = d[m.end():]
    if _MATERIAL_TAIL_BLOCK_RE.search(food):
        return d
    # On an animal the material is what the meat was cooked IN or coated WITH,
    # not what it is made of: FDC's "Chicken, ... cooked, fried, flour" folded
    # to 'chicken flour', and there is no such food.
    if _MEAT_CONTEXT_RE.match(food):
        return d
    if head == "butter" and _MATERIAL_TAIL_BEAN_RE.search(food):
        return d
    if head == "milk" and food == "chocolate" and not _MILK_BEVERAGE_RE.match(rest):
        return d
    # the food keeps its own number: STFCJ writes "Lemons, juice" and McCance
    # "Milk, goats", which folded to 'lemons juice' and 'goats milk'
    food = _fold_plural(food)
    if head in food.split():
        return food + rest
    return f"{food} {head}{rest}"


# Fineli repeats the material in the modifier: "Flour, Cornstarch, Cornflour",
# "Flour, Potato Starch, Potato Flour", "Flour, Corn Meal, Coarsely Ground".
# Moving the head in front of those gives 'cornstarch flour' and 'corn meal
# flour', which name the material twice.
_SELF_NAMING_TAIL_RE = re.compile(
    r"(?:flour|starch|meal|crumbs?|bran|semolina|oil|juice|milk|butter)$", re.I)


def _material_head_swap(m: re.Match) -> str:
    # the modifier keeps its own number: McCance writes "Milk, goats"
    head, mod = m.group(1).lower(), _fold_plural(m.group(2).strip().lower())
    # "Flour, cornflour" and "Flour, wheat flour" name the material twice; the
    # head is redundant rather than misplaced, so it is dropped, not moved.
    if mod.endswith(head) or _SELF_NAMING_TAIL_RE.search(mod):
        return mod
    return f"{mod} {head}"


# 2. SPECIES read genus-first, which is the opposite order and the opposite
#    majority: 'salmon, coho' has 11 members against 1 for 'coho salmon',
#    'herring, atlantic' 8 against 1, 'pike, northern' 6 against 2. Genus-first
#    also keeps the species of one fish together in an alphabetical index,
#    which is what a composition table is read from. Only these heads, and only
#    when the whole canon is "<modifier> <head>" with nothing else in it.
_SPECIES_GENUS_RE = re.compile(
    # the modifier: one or two plain words
    r"^([a-z]+(?:\s+[a-z]+)?)\s+"
    r"(salmon|cod|herring|mackerel|tuna|mullet|pike|crab|octopus|pollock|"
    r"sprat|roughy|sole|flounder|snapper|perch|bass|trout|whiting|shark|"
    r"squid|shrimp|prawn|lobster|clam|oyster|mussel|scallop|abalone)(?=,|$)", re.I)
# ... except where the two words are the accepted name of the fish rather than
# a modifier plus a genus. "Horse mackerel" and "jack mackerel" are species
# names in their own right, and "rock salmon" is a trade name for dogfish.
# "Norway lobster" is Nephrops and "pike-perch" is Sander: filing either behind
# the genus in front of it would put it with animals it is not related to.
# The modifier is never a connective. Without this "pizza with tuna" reads as
# modifier "pizza with" plus genus "tuna" and becomes 'tuna, pizza with', and
# "shrimp or prawn" becomes 'prawn, shrimp or'. Tested against the MODIFIER
# only: an earlier version tested the whole canon, so a "without added fat"
# further along the name silently switched the rule off and left
# 'atlantic herring, without added fat, fried' beside 'herring, atlantic'.
_SPECIES_CONNECTIVE_RE = re.compile(
    r"\b(?:with|without|or|and|in|of|from|for|on)$", re.I)
_SPECIES_GENUS_KEEP = re.compile(
    r"^(?:horse|jack|blue\s+jack|rock|spanish\s+jack|round\s+scad|norway|pike|"
    # a DISH is not a modifier: "sauce oyster" and "salad tuna" name the dish,
    # and swapping them puts the fish in front of its own recipe
    r"salad|sandwich|soup|stew|pizza|sauce|pate|paste|cake|burger|salt|"
    # nor is a PREPARATION - those belong to the state machinery, which appends
    # them behind the food anyway, so moving them here would race it
    r"breaded|battered|smoked|dried|cooked|fried|grilled|canned|salted|"
    r"pickled|raw|fresh|frozen|boiled|steamed|minced|potted|imitation)$", re.I)

# The salt and sugar axes are read off the description before it is stripped, so
# the wording that carried them has to be removed from the name afterwards.
# Without this the FDC phrasing survives as a chunk of its own: "Blackeye pea,
# canned, sodium added" produced the canon 'blackeye pea, sodium added', and
# _append_state then declined to add ", salted" because the token it matched was
# already present in the name. Only the axis PHRASES go - a bare "sugar" or
# "salt" is left alone, because there the word can be the food.
_AXIS_PHRASE_RE = re.compile(
    r"\s*,?\s*\b(?:no\s+)?(?:sodium|sugars?|salt)\s+added\b(?:\s+in\s+processing)?"
    # The negative has to be swallowed WITH the phrase. Matching only "added
    # salt" inside "no added salt" left the bare word "no" standing as a chunk,
    # so 'peanut, no' and 'muesli, no' sat beside 'peanut' and 'muesli' - the
    # very split the axis exists to prevent, since unsalted and unsweetened are
    # silent labels and both rows belong on the bare name.
    r"|\s*,?\s*\b(?:with|no|without)?\s*added\s+(?:sodium|sugars?|salt)"
    r"(?:\s*(?:and|&|/)\s*(?:sodium|sugars?|salt))?\b"
    # FDC says WHERE it was added - "Potatoes, french fried, salt added in
    # processing" - and the trailing half was left standing as a chunk of the
    # name ('potato, in processing, salted')
    r"(?:\s+in\s+processing)?"
    # ...and the same claim spelled with "salt". Fineli writes "Less Salt" and
    # NEVO "low salt"; the label was read but the wording stayed, so
    # _append_state saw its own token still standing and printed nothing.
    r"|\s*,?\s*\b(?:reduced|low|lower|less)[-\s](?:sodium|salt)"
    r"(?:\s*(?:and|&|/)\s*sugars?)?\b"
    r"|\s*,?\s*\bwith(?:out)?\s+(?:added\s+)?sugars?\b"
    r"|\s*,?\s*\b(?:un|non)-?sweetened\b"
    r"|\s*,?\s*\bsugar[-\s]?free\b|\s*,?\s*\bsodium[-\s]?free\b"
    # "lightly salted" is read as the low-salt level, then _PREP_RE takes the
    # "salted" and leaves "lightly" standing as a chunk of its own
    r"|\s*,?\s*\blightly\s+salted\b"
    # "canned, in heavy syrup" is read as sweetened, but only the comma-less
    # spelling reached _PREP_RE, so "canned in heavy syrup" kept the wording and
    # _append_state then declined to print the label
    r"|\s*,?\s*\bin\s+(?:extra\s+)?(?:heavy|light|medium)?\s*syrups?\b"
    # FDC couples the two claims - "Popcorn, microwave, low fat and sodium" -
    # and taking the fat half left "and sodium" standing as a chunk
    r"|\s*\band\s+sodium\b"
    # ...and Fineli's "Fried Without Fat And Salt", where the fat half is read
    # by the fat axis and the rest was left standing as a chunk
    r"|\s*\band\s+salt\b(?=\s*(?:,|\[|$))"
    # "baked without fat" is the same claim as the "low fat" already read
    r"|\s*,?\s*\b(?:reduced|low|lower|less)\s+sugars?\b", re.I)

_LAB_ROW_RE = re.compile(
    r"[-,]\s*(?:NF\w+|CY\w*\d\w*|\d+[a-z]?(?:-\d+)+-[a-z]+\w*)\s*$", re.I)
_PANEL_FOOD_PREFIX_RE = re.compile(r"^\s*(?:sugars|sweets|total\s+sugars)\s*,\s*", re.I)

_ANALYTE_PREFIX_RE = re.compile(
    r"^\s*" + _ANALYTE + r"(?:\s*(?:,|&|and)\s*" + _ANALYTE + r")*\s*(?:[-\u2013]\s+|,\s*)",
    re.I)

_PANEL_HEADS = frozenset({
    "proximates", "minerals", "vitamins", "fatty acids", "amino acids", "ash",
    "moisture", "nitrogen", "energy", "carbohydrates", "niacin", "riboflavin",
    "thiamin", "cholesterol", "pantothenic acid", "retinol", "folate", "choline",
    "selenium", "tocopherol",
    # NOT panel heads, though FDC uses them as such: each is also a food in
    # its own right, and stripping it deletes the only name the row has.
    #   "starch", "sugars", "sugar", "salt", "water",
    "pectin", "protein", "vitamin a", "vitamin b6", "vitamin b12", "vitamin c",
    "vitamin d", "vitamin e", "vitamin k", "calcium", "iron", "zinc", "magnesium",
    "phosphorus", "potassium", "sodium", "copper", "manganese",
})


def _is_panel_head(x: str) -> bool:
    # Singular form too: this runs before plural folding, so the source spelling is
    # whatever the panel used ("Tocopherols, Salsa", "Folates, Lettuce, Romaine").
    x = x.lower().lstrip('& ').strip()
    return x in _PANEL_HEADS or _singularize(x) in _PANEL_HEADS


def _strip_panel_prefix(d: str) -> str:
    if _BRACKET_TAG_RE.search(d):
        return d
    parts = [x.strip() for x in d.split(",")]
    i = 0
    while i + 1 < len(parts):
        nxt = parts[i + 1]
        if not _is_panel_head(parts[i]):
            break
        if not (nxt[:1].isupper() or _is_panel_head(nxt)):
            break
        i += 1
    return ", ".join(parts[i:]) if i else d


# "cookies" -> "cooky" is the single most visible naming bug in the corpus: one
# slice alone carried 40 canons, and it also split the food in two, because CNF
# spells it "cookie" and folds correctly. English cannot distinguish these from
# berry/berries by pattern, so the -ie words are listed.
_IE_PLURALS = {
    "cookies": "cookie", "brownies": "brownie", "calories": "calorie",
    "smoothies": "smoothie", "hoagies": "hoagie", "veggies": "veggie",
    "kedgerees": "kedgeree", "brie": "brie",
}
# Foreign and dish names that merely END in s. Stripping it invents a non-word
# and, worse, a second canon for a food that already has one.
_FOREIGN_NOT_PLURAL = frozenset({
    "speculaas", "kalops", "sprits", "couscous", "hummus", "houmous", "tapas",
    "nachos", "tacos", "churros", "gyros", "hoummos", "falafels", "biryanis",
    "foie gras", "gras", "tostitos", "doritos", "cheerios", "pringles",
    "bitterleaves", "deurmekaarbos", "moringaleaves", "calyces", "matzos",
    "sos", "mars", "twix", "ritz", "oreos", "lays",
})
_IRREGULAR_SINGULAR = {
    "calyces": "calyx", "bitterleaves": "bitterleaf",
    "moringaleaves": "moringaleaf", "deurmekaarbos": "deurmekaarbos",
}


def _singularize(w: str) -> str:
    # Possessives are not plurals. Without this McDONALD'S -> "mcdonald'", and
    # 93 canons across the corpus end in a bare apostrophe.
    if w.endswith("'s") or w.endswith("\u2019s"):   return w
    if w in _IRREGULAR_SINGULAR:          return _IRREGULAR_SINGULAR[w]
    if w in _IE_PLURALS:                  return _IE_PLURALS[w]
    if w in _FOREIGN_NOT_PLURAL:          return w
    if w in _VES_SINGULAR:                return _VES_SINGULAR[w]
    if w in _NOT_PLURAL:                  return w
    if len(w) < 4 or not w.endswith("s"): return w
    if w in _IS_PLURAL:                   return w[:-1]
    if re.search(r"(?:ss|us|is)$", w):    return w
    if w.endswith("ies") and len(w) > 4:  return w[:-3] + "y"
    if w.endswith("oes") and len(w) > 4:  return w[:-2]
    if _SIBILANT_ES_RE.search(w):         return w[:-2]
    return w[:-1]


def _fold_plural(chunk: str) -> str:
    """Singularize the final word only. "green beans" -> "green bean", but
    "solids and liquids" keeps its leading plural, which is enough to merge
    the pair since both spellings share it."""
    if " " not in chunk:
        return _singularize(chunk)
    head, _, last = chunk.rpartition(" ")
    return f"{head} {_singularize(last)}"


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
    # Accession / breeding codes that the shapes above miss. Each alternative
    # is anchored tightly enough that no food qualifier can match it:
    r"|[a-z]\d[a-z]\d"                               # T1A1, T1B2 (maize hybrids)
    r"|[a-z]{1,12}[-\s]?\d+[a-z]?"                   # bosworth-3, 17w, adt-21
    r"|\w+\s+\d+\s*[a-z]\b"                        # texas 17 w
    # DMR-ESR-W. The two lookaheads are what keep real hyphenated food
    # words out: a breeding acronym always carries at least one vowel-less
    # segment, which "ready-to-drink", "ready-to-eat" and "oil-in-water"
    # never do.
    r"|(?=[a-z-]{3,30}$)(?=(?:[a-z]*-)*[bcdfghjklmnpqrstvwxz]{3,6}(?:-|$))"
    r"[a-z]{1,6}(?:-[a-z]{1,6}){1,4}"
    r"|[^,]*\bn\s*\u00b0\s*\d+[^,]*"               # IRNAS n° 11
    r"|[a-z]{2,}(?:/[a-z0-9]{1,6})+/(?:\s*[a-z0-9-]{1,6})?"  # KARI/BN/ BK- (slash
                                                     #  required: never matches
                                                     #  the peach/nectarine form)
    r"|[^,]*\b(?:h[\u00ed i]brido|hybrid|cultivar|landrace|accession|clone)\b[^,]*"
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
    "faba bean", "faba beans", "broad bean", "broad beans", "cowpea", "cowpeas",
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

def _strip_venue(d: str) -> str:
    """Remove a leading restaurant / brand / programme prefix, if any."""
    d = _VENUE_PREFIX_RE.sub("", d, count=1)
    return _VENUE_GENERIC_RE.sub("", d, count=1)




# ---------------------------------------------------------------------------
# MANUAL CURATION
#
# Some names cannot be reached by rule. A brand is only recognisable as a brand
# if you know the brand; "all flavors except chocolate" is filler only because a
# person reads it as filler. Rather than bend the general rules until they
# overfit, the irreducible cases are listed here. Two layers, cheapest first.
# ---------------------------------------------------------------------------

# Commercial product lines that survived the venue/brand prefix rules because
# they sit MID-name rather than at the front. Stripped as whole words so a real
# food word is never eaten.
# Curated by hand from ~1,200 BRAND findings across ten parallel audits of the
# full canon list. Token-frequency mining alone was NOT safe to accept: its top
# hits included "cheese", "milk", "candy", "crisp", "flake" and "snack", which
# are food words, and deleting them would have destroyed real names. Every entry
# below is a company or product line, matched as a whole word.
_BRAND_NAMES = [
    # spreads / margarine
    "smart beat", "smart balance", "benecol", "shedd's spread", "country crock",
    "i can't believe it's not butter", "brummel and brown", "fleischmann's",
    "blue bonnet", "parkay", "promise", "take control", "olivio", "becel",
    "blue band", "flora", "keiju", "vita d'or", "l\u00e4tta", "latta",
    # cereals / snacks
    "quaker", "malt-o-meal", "kellogg's", "kellogg", "nabisco", "ralston",
    "health valley", "nature's path", "post", "general foods",
    "cap'n crunch", "krispy kreme", "keebler", "pillsbury", "general mills",
    "oreo", "ryvita", "finn crisp", "kexchoklad", "pringles", "tostitos",
    "doritos", "m&m's", "reese's", "mars", "twix", "spam",
    # dairy / soy
    "silk", "vitasoy", "alpro", "nasoya", "valio", "gefilus", "chobani",
    "breakstone", "daisy", "oumph", "isostar",
    # brands / manufacturers
    "kraft", "nestl\u00e9", "nestle", "heinz", "hormel", "oscar mayer",
    "swanson", "snapple", "coca-cola", "guinness", "weight watcher",
    "healthy choice", "stove top", "udi's", "van's", "atria", "snellman",
    "saarioinen", "kokkikartano", "pirkka", "liga",
    # infant / toddler formula
    "abbott", "gerber", "mead johnson", "nutrilon", "hipp", "hero baby",
    "bonbebe", "aptamil", "nutricia",
    # dutch retailers (NEVO files store brands as the product name)
    "albert heijn", "kruidvat", "etos", "jumbo", "bonus",
    # Round 4: trading names left standing once the corporate entity chunk
    # around them was removed, or filed by FDC as the product name itself.
    # A name belongs here only when it is the WHOLE informative chunk. Where the
    # brand prefixes a product line ("ENFAMIL, ENFACARE", "Flax Plus Maple Pecan
    # Crunch", "Nissin, Cup Noodle") removing it does not clean the name, it just
    # promotes the sub-brand: "infant formula, ready-to-feed" shattered into
    # fifteen canons named after formula lines. Those are already handled by
    # curated overrides and are deliberately absent from this list, as are the
    # names that only leave a bare "brand" behind ("soymilk, westsoy brand" ->
    # "soymilk, brand", and "gensoy"/"good sense" then collide on "soynut,
    # brand", which would fold two different products into one).
    "dasani", "minute maid", "nos zero", "slimfast", "unilever", "digiorno",
    "hot pocket", "croissant pockets", "fuze",
    # Round 14: Fineli files the manufacturer and the product line as chunks of
    # the name, so a food's canon was its trademark - twelve of them for rye
    # bread alone ('rye bread, oululainen jalkiuunipala', 'rye bread, reissumy',
    # 'rye bread, taikaruis ruisviipale'), and NO bare 'rye bread' canon existed
    # at all. Mined by taking the tokens that appear in Fineli rows and in no
    # other source, then reading every one: the mine also returned "talkkuna"
    # (roasted grain flour), "powan" (a whitefish), "kassler" (smoked loin),
    # "rypsi" (turnip rape), "kissel" (a fruit soup) and "milkcap" (a Lactarius
    # mushroom), which are foods, and none of those is in this list.
    "fazer", "vaasan", "oululainen", "leip\u00e4aitta", "linkosuo",
    "p\u00e5gen", "pagen", "perheleipurit", "jyv\u00e4shyv\u00e4",
    "lauantai", "ruispalat", "reissumies", "reissumy", "taikaruis",
    "maalaisviipaleet", "moniviljaviipaleet", "t\u00e4yshyv\u00e4t",
    "aittaleip\u00e4", "aamiaisleip\u00e4", "oivallus", "ruiskakko",
    "j\u00e4lkiuunipala", "j\u00e4lkiuunileip\u00e4",
    "j\u00e4lkiuuniviipaleet", "ruislimppu", "tosi rukiinen",
    # retailers and own-labels
    # ("rainbow", an S-market own-label, is deliberately absent: rainbow trout)
    "lidl", "x-tra",
    # dairy, drinks and desserts
    "oatly", "yosa", "risifrutti", "fr\u00f6dinge", "mehukatti", "jacky",
    "novelle", "oivariini", "ingmariini", "luonto",
    # meals, meat and catering
    "felix", "antell", "hesburger", "kartanon", "ullan pakari", "pakari",
    "verso food",
    # meal-replacement lines
    "allevo", "nutrilett", "naturdiet", "profeel",
]
_BRAND_RE = re.compile(r"\b(?:" + "|".join(re.escape(b) for b in _BRAND_NAMES) + r")\b", re.I)

# "Margarine-like spread" is FDC's category label, not a food name. Ten canons
# carried it; a person calls all of them margarine.
_SPREAD_LIKE_RE = re.compile(r"\bmargarine[-\s]like(?:\s+spread)?\b", re.I)

# Non-food fragments that survived into canon names across many sources:
# NEVO's serving basis ("... standaard 2 p 100 ml"), the Swedish "e.g." marker
# that introduces a dish name, and Dutch/Finnish grading words.
_SERVING_BASIS_RE = re.compile(r"\bp\s*\d+\s*(?:ml|g)\b", re.I)
_EG_MARKER_RE = re.compile(r"\be\.\s*g\.\s*", re.I)
# Heat treatment and homogenisation are process, not composition, and CIQUAL
# writes them as a chunk of their own: "Milk, semi-skimmed, pasteurised" took
# the qualifier slot the fat level needed.
_GRADE_WORD_RE = re.compile(
    r"\b(?:standaard|biologisch|naturel|pasteuri[sz]ed|pasteuri[sz]ee|"
    r"homogeni[sz]ed|sterili[sz]ed|ultra[-\s]?pasteuri[sz]ed)\b", re.I)

# Marketing intensifiers. They grade a product against its own range, they do
# not describe the food: "SMART BEAT Super Light", "Light Buttery Spread".
# "extra virgin" is a legal pressing grade, not an intensifier, and it is the
# grade that carries the polyphenols: dropping the "extra" folded FDC's extra
# virgin olive oil into the plain 'olive oil' canon beside refined oil.
_MARKETING_RE = re.compile(
    r"\b(?:super|extra(?![-\s]virgin\b)|ultra|premium|deluxe|classic|gourmet|buttery|"
    r"original\s+recipe|new\s+improved)\b", re.I)
# Once "margarine-like spread" has become "margarine", a trailing "spread" is
# the same word twice: "margarine, light spread" -> "margarine, light".
# Applied only when the name already says margarine, so "cheese spread" and
# "sandwich spread" - where the word IS the food - are untouched.
_REDUNDANT_SPREAD_RE = re.compile(r"\s*\bspreads?\b", re.I)
# Hedges left stranded when the quantity they qualify is removed:
# "Margarine-like spread, approximately 60% fat, tub" -> "margarine, approximately".
_HEDGE_RE = re.compile(r"\b(?:approximately|approx\.?|app\.|about|circa|ca\.|at\s+least|at\s+most|min\.|max\.)\s*", re.I)

# Exclusionary flavour phrases classify a group, they do not name a food:
# "Puddings, all flavors except chocolate, low calorie, instant, dry mix".
_FLAVOUR_FILLER_RE = re.compile(
    r"\b(?:all|any)\s+(?:flavou?rs?|varieties|kinds?)\s+"
    r"(?:except|other\s+than|but)\s+[a-z]+\b", re.I)

# Final say. Keyed on the canon the rules produce, applied last, so a rule
# change cannot silently strand an entry: if the left-hand side stops being
# produced the override simply stops firing (audited by test_canonicalize.py).
# ---------------------------------------------------------------------------
# CROSS-DATABASE SYNONYMS
#
# Sixteen national databases name the same food differently, so one food arrives
# as several canons and its evidence is split. Direction of each merge is chosen
# by COUNTING the canons that already use each spelling, not by preference, so
# the table changes as few names as possible.
#
# Deliberately NOT in this table, though they look like obvious pairs:
#   maize -> corn      "corn" is also sweetcorn, a fresh vegetable; merging the
#                      dry grain into it would cross two very different foods.
#   beet <-> beetroot  "beet" is also sugar beet.
#   prawn -> shrimp    strictly different taxa; awaiting the seafood audit.
#   rocket/arugula     only 14 canons, and "rocket" is ambiguous.
# ---------------------------------------------------------------------------
_SYNONYM_PAIRS = {
    # variant             : canonical (the spelling more canons already use)
    "yogurt":               "yoghurt",     # 60 canons -> 99
    "yogourt":              "yoghurt",     # CNF's spelling; 13 canons stood apart
    "donut":                "doughnut",    # 1 -> 34
    "garbanzo beans":       "chickpea",    # 2 -> 18
    "garbanzo bean":        "chickpea",
    "garbanzo":             "chickpea",
    "groundnut":            "peanut",      # 50 -> 119
    "courgette":            "zucchini",    # 10 -> 22
    "aubergine":            "eggplant",    # 15 -> 19
    # Round 12: the same food under two regional names, both present in the
    # index. Only pairs that are the SAME SPECIES are here.
    "rocket":               "arugula",     # Eruca vesicaria; 5 canons -> 2
    "swede":                "rutabaga",    # Brassica napobrassica; 16 -> 4
    "capsicum":             "bell pepper",  # Capsicum annuum; the merges then
                                            # take it to 'pepper, sweet'
    "faba bean":            "broad bean",  # Vicia faba; 4 -> 25
    "fava bean":            "broad bean",
    "haricot bean":         "navy bean",   # Phaseolus vulgaris white; 9 rows -> 118
    "linseed":              "flaxseed",    # Linum usitatissimum
    "cornflour":            "cornstarch",  # UK cornflour IS US cornstarch
    "bicarbonate of soda":  "baking soda",
    "pak choi":             "bok choy",    # Brassica rapa chinensis
    "sweetcorn":            "sweet corn",
    "cooking banana":       "plantain",
    "soya beans":           "soybean",     # ... and "soy bean" is "soybean"
    "soya bean":            "soybean",
    "soy beans":            "soybean",
    "soy bean":             "soybean",
    "soya":                 "soy",         # 76 canons -> 196
    "maize":                "corn",        # Zea mays; 92 canons -> 248
    # NOT here, and each for a reason: "biscuit"/"cookie" (a US biscuit is a
    # bread), "marrow"/"squash" (marrow is also the bone), "saithe"/"pollock"
    # and "hake"/"whiting" (different species), "chicory"/"endive" (Cichorium
    # intybus against C. endivia), "treacle tart" (a dish, not molasses), and
    # "prawn"/"shrimp", where the two are different families.
}
_SYNONYM_RE = [
    (re.compile(r"\b" + re.escape(k) + r"\b", re.I), v)
    # longest first, so "garbanzo bean" is matched before "garbanzo"
    for k, v in sorted(_SYNONYM_PAIRS.items(), key=lambda kv: -len(kv[0]))
]

_NAME_REWRITES = [
    # Botanically NOT spinach: New Zealand spinach is Tetragonia tetragonioides,
    # not Spinacia oleracea, and its composition differs. Merged at the explicit
    # request of the project owner (2026-08-21); delete this entry to separate
    # them again.
    (re.compile(r"\bnew\s+zealand\s+spinach\b", re.I), "spinach"),
    # CNF files this sausage as "Sausage, Knackwurst (Knockwurst)" while FDC
    # files it as "Knackwurst, knockwurst" - one food, two canons.
    (re.compile(r"\bsausage,\s*knackwurst\b", re.I), "knackwurst"),
    # Fast-food product LINES. Once the chain name is stripped, what is left is
    # the marketing size, not the food: "WENDY'S, Jr. Hamburger, with cheese"
    # became 'jr. hamburger, with cheese' and "WENDY'S, CLASSIC DOUBLE, with
    # cheese" became 'double, with cheese' - a canon with no food noun at all.
    # Values are per 100 g, so a Jr. or a Single is compositionally the same
    # burger and the size word goes. "Double" stays: two patties against one bun
    # is a real change in the meat-to-bread ratio, hence in the composition.
    (re.compile(r"\bclassic\s+double\b|\bdouble\s+stack\b", re.I), "double hamburger"),
    (re.compile(r"\bclassic\s+single\s+hamburger\b|\bclassic\s+single\b", re.I), "hamburger"),
    (re.compile(r"\bjr\.?\s+hamburger\b", re.I), "hamburger"),
    (re.compile(r"\bdave's\s+hot\s*'?n'?\s*juicy\s*[\d/]*\s*lb,?\s*", re.I), "hamburger, "),
    (re.compile(r"\bwhopper\s+jr\.?\b", re.I), "whopper"),
    # Same claim, two spellings: "no cheese" / "without cheese" split one food.
    (re.compile(r"\bno\s+(?=cheese\b|mayo\b|mayonnaise\b|sauce\b|salt\b)", re.I), "without "),
]

_CANON_OVERRIDES = {
    "'s free sour cream, fat-free": 'sour cream, fat-free',
    "'s halloween crunch": 'cereal, ready-to-eat, sweetened corn and oat',
    "'s oops! all berries cereal": 'cereal, ready-to-eat, sweetened corn and oat with berry',
    "'s peanut butter crunch": 'cereal, ready-to-eat, peanut butter flavoured corn and oat',
    "'s sour cream, low-fat": 'sour cream, low-fat',
    '4-grain flake, riihikosken vehnämylly': '4-grain flake',
    '9 oz house sirloin steak': 'sirloin steak',
    'abalone or ormer or sea ear': 'abalone',
    'abricot, and soaked, dried': 'apricot, soaked, dried',
    'achicoria azul': 'blue chicory',
    'african starapple, flesh': 'african star apple, flesh',
    'african starapple, peel': 'african star apple, peel',
    'after eight chocolate mint': 'chocolate mint thins',
    "afzelia 's oil": 'afzelia seed oil',
    'afzelia africana, shelled': 'afzelia bean, shelled',
    # Ackee is named like any other fruit - apple, banana - not "ackee fruit".
    # Three sources, three spellings: McCance says "Ackee", PhyFoodComp says
    # "Ackee fruit", WAFCT says "Akee, fruit". Both of the latter are pointed
    # straight at 'ackee' rather than chained through one another, so neither
    # depends on the single hop _consult_tables allows.
    'ackee fruit': 'ackee',
    'akee, fruit': 'ackee',
    'almond cake product': 'almond cake, frozen',
    'almond cookie bittersweet mandelkubb': 'almond cookie, bittersweet',
    "almond, ann's house of nut": 'almond, blanched, oil-roasted',
    'almond, blue diamond - ky, salted': 'almond, dry roasted, salted',
    'almond, emerald - ga, salted': 'almond, dry roasted, salted',
    'almond, nonpareil marketing': 'almond, nonpareil variety',
    'almond, unblanched, salted': 'almond, salted',
    'alpen': 'muesli cereal, ready-to-eat, with dried fruit',
    'alpha-bit': 'cereal, ready-to-eat, oat and corn alphabet shape',
    'always tender, center chop': 'pork chop, center cut',
    'always tender, pork loin': 'pork, loin, boneless',
    'always tender, pork loin filet': 'pork loin fillet, lemon garlic-flavoured',
    'always tender, pork tenderloin': 'pork, tenderloin, peppercorn-flavoured',
    'amaranthus, leaf': 'amaranth leaf',
    'amber, hard cider': 'hard cider',
    'american process cheese product, individually wrapped': 'american process cheese product',
    'american yam bean root, irnas n°': 'american yam bean root, irnas n° 11',
    'american yam bean root, local': 'american yam bean root',
    'and peanut, dried, unsalted': 'seed and peanut, dried',
    "andrea's, gluten-free soft dinner roll": 'dinner roll, soft, gluten-free',
    'animal fat, native, dried': 'moose fat, dried',
    'appel pie dutch with shortbread with butter': 'apple pie dutch with shortbread with butter',
    'apple kissel, apple soup, dried': 'apple kissel (apple soup), made from dried apple',
    'apple pie dutch with shortbread average': 'apple pie, dutch, with shortbread',
    'apple sauce without sugar tinned': 'applesauce, unsweetened',
    'apple zing': 'cereal, ready-to-eat, apple cinnamon flavoured',
    'apples with oil or liquid margarine sugar almond paste walnut': 'baked apple with oil or liquid margarine sugar almond paste walnut',
    'apricot oil': 'apricot kernel oil',
    'archway home style cookie, dutch cocoa': 'cookie, dutch cocoa',
    'arizona, tea': 'tea, black, ready-to-drink, lemon',
    'aronia berry': 'chokeberry',
    'arracacia, tuber': 'arracacha, tuber',
    "arrow leaf's elephants ear": 'arrowleaf elephant ear',
    'arroz con frijoles negro': 'arroz con frijoles negros (rice with black bean)',
    'arroz con grandule': 'arroz con gandules (rice with pigeon pea)',
    'arroz con habichuelas colorado': 'arroz con habichuelas coloradas (rice with red bean)',
    'artemisinin, leaf': 'sweet wormwood, leaf',
    'arugula arugula salad': 'arugula',
    'arugula, rocket': 'arugula',
    'avocado fruit': 'avocado',
    'avocado pear': 'avocado',
    'avocado, all commercial variety': 'avocado',
    'b.l.t. sub on white bread with bacon, lettuce and tomato': 'sub sandwich on white bread with bacon, lettuce and tomato',
    'baby food, commercially produced, flour': 'baby food, commercial, with carrot, rice flour, beef and pea',
    'babyfood, baby mum mum rice biscuit': 'babyfood, rice biscuit',
    'babyfood, graduates lil biscuits vanilla wheat': 'babyfood, vanilla wheat biscuit',
    'babyfood, product': 'babyfood, baked cereal finger snack, fortified',
    'babyfood, product, enriched': 'babyfood, baked cereal finger snack, fortified',
    'back velvet tamarind': 'black velvet tamarind',
    'bacon rasher, back': 'bacon, back',
    'bacon, egg & cheese mcgriddle': 'breakfast sandwich on pancakes with bacon, egg and cheese',
    'bacon-wrapped frankfurt sausage, freid': 'bacon-wrapped frankfurt sausage, fried',
    "baker's yeast granualted, dried": "baker's yeast, dried",
    'baking chocolate, snackfood us': 'baking chocolate, milk chocolate mini baking bits',
    'bambara groundnut, dried': 'bambara peanut, white, dried',
    'bambara groundnut, dried, red': 'bambara peanut, red, dried',
    'bambara peanut, combined variety, dried': 'bambara peanut, dried',
    'bambara peanut, for 120 min, dried': 'bambara peanut, boiled for 120 min, dried',
    'bambara peanut, for 30 min, dried': 'bambara peanut, boiled for 30 min, dried',
    'bambara peanut, for 60 min, dried': 'bambara peanut, boiled for 60 min, dried',
    'bambara peanut, for 90 min, dried': 'bambara peanut, boiled for 90 min, dried',
    'bamboo shoot, nonsalted': 'bamboo shoot',
    'banakou né barâand kanss saagbo, foutou of cassava and plantain': 'banakou né barâand kanss saagbo, foutou of cassava and unripe plantain',
    'banquet, salisbury steak with gravy': 'salisbury steak with gravy, frozen, unprepared',
    'baobab fruit, flesh, dried': 'baobab, flesh, dried',
    'baobab, dead-rat tree/monkey-bread tree/upside-down tree/cream of tartar tree': 'baobab fruit/monkey bread, flesh',
    'bar and own brand equivalent': 'chocolate bar with caramel and nougat',
    "barbara's puffin, original": 'cereal, ready-to-eat, puffed corn and oat, original',
    'barcelona nut': 'hazelnut',
    'barley bran, extruded at 115ºc and 20%h20': 'barley bran, extruded at 115ºc and 20% h2o',
    'barley flour wholegrain': 'barley flour, wholegrain',
    'barley, plpa': 'barley, wholegrain',
    'barley, plpb': 'barley, wholegrain',
    'barley, sunnita': 'barley, irradiated, water-soaked and incubated',
    'bean curd': 'tofu',
    'bean, black eyed, dried': 'bean, blackeye, dried',
    'bean, blackeye pea, salted': 'blackeye pea, salted',
    'bean, chick pea': 'chickpea',
    'bean, chick pea, dried': 'chickpea, dried',
    'bean, chick peas/garbanzo': 'chickpea',
    'bean, chick peas/garbanzo, salted': 'chickpea, salted',
    'bean, chickpeas/garbanzo, dried': 'chickpea, dried',
    'bean, edamame': 'edamame',
    'bean, fava': 'broad bean',
    'bean, french tinned': 'french bean, tinned',
    'bean, french, salted': 'bean, french, mature seed, boiled, salted',
    'bean, in barbecue sauce': 'baked bean, canned in barbecue sauce',
    'bean, in tomato sauce': 'baked bean, canned in tomato sauce',
    'bean, kidney red': 'bean, kidney, red',
    'bean, kidney red, dried': 'bean, kidney, red, dried',
    'bean, legume': 'mung bean',
    'bean, legume, salted': 'mung bean, salted, sprouted',
    'bean, long yard kousenband': 'yardlong bean',
    'bean, mung': 'mung bean',
    'bean, pigeon pea, dried': 'pigeon pea, dried',
    'bean, white in tomato sauce': 'baked bean, in tomato sauce',
    'bean, with pork and sweet sauce': 'baked bean, canned with pork and sweet sauce',
    'bean, with pork and tomato sauce': 'baked bean, canned with pork and tomato sauce',
    'beansprout, soy': 'soybean sprout',
    "bearnaise sauce mix, campbell's/blå band": 'bearnaise sauce mix',
    'bearnaise sauce mix, knorr': 'bearnaise sauce mix, dried',
    'bearnaise sauce rte heated or instant powder': 'bearnaise sauce, from frozen or instant powder, heated',
    'bearnaise sauce rte heated product or instant powder': 'bearnaise sauce, from frozen or instant powder, heated',
    'beef av': 'beef',
    'beef bonless beef outside': 'beef, boneless outside round',
    'beef breakfast sausage, banquet brown n serve sausage link': 'beef breakfast sausage link, brown and serve',
    'beef breakfast sausage, eckrich sausage link': 'beef breakfast sausage link',
    'beef meat, ground': 'beef mince',
    'beef ox tongue': 'beef, tongue',
    'beef seasoned pastirma turkish, dried, salted': 'beef pastirma, dried, salted',
    'beef smoke, dried': 'beef, smoke-dried',
    'beef smoke- lightly, dried, salted': 'beef, smoke-dried, lightly salted',
    'beef steak without gravy produkt': 'beef steak, roasted without gravy, frozen',
    'beef stew homemade kalops': 'beef stew, homemade',
    'beef stew kalops': 'beef stew, frozen',
    'beef stew product kalops': 'beef stew, frozen',
    'beef stew with potatoes onion beer produkt sjömansbiff': 'beef stew with potato, onion and beer, frozen',
    'beef stock paste or powder with large-scale, low-salt': 'beef stock paste or powder with reduced salt',
    'beef stock paste or powder with reduced salt large-scale': 'beef stock paste or powder with reduced salt',
    'beef stroganoff product': 'beef stroganoff, frozen',
    'beef tender loin': 'beef, tenderloin',
    'beer alcohol-free <0, 1 vol%': 'beer, alcohol-free, <0.1 vol%',
    'beer alcohol-free vol % 0.5': 'beer, alcohol-free, 0.5 vol%',
    'beer low alcohol, 2 vol%': 'beer, low alcohol, 0.1-1.2 vol%',
    'beer vol. % 5.4': 'beer, 5.4 vol%',
    'beer, non-alcoholic': 'beer, alcohol-free',
    'beesting, colostrum': 'bovine colostrum',
    'beet': 'beetroot',
    'beet root': 'beetroot',
    'beet, pickled': 'beetroot, pickled',
    'beetroot, red beet': 'beetroot',
    'beetroot, salt': 'beetroot, salted',
    'beetrot, leaf': 'beetroot, leaf',
    'belgium endive': 'chicory, witloof',
    'bell pepper green': 'pepper, sweet, green',
    'bell pepper green red product': 'bell pepper, green and red, frozen',
    'bell pepper red': 'pepper, sweet, red',
    'bell pepper yellow': 'pepper, sweet, yellow',
    'bell pepper, green': 'pepper, sweet, green',
    'bell pepper, red': 'pepper, sweet, red',
    'bell pepper, yellow': 'pepper, sweet, yellow',
    'bengal dayflower, leaf': 'benghal dayflower, leaf',
    'bengal gram, whole': 'chickpea',
    'benghal, leaf': 'benghal dayflower, leaf',
    'berberis, fruit': 'barberry, fruit',
    'berry colossal crunch': 'cereal, ready-to-eat, sweetened corn with berry flavour',
    'berry drink, proviva': 'berry drink with probiotics',
    'beverage, coconut water': 'coconut water',
    'big breakfast': 'breakfast platter with egg, sausage, hash brown and biscuit',
    "big daddy's ls 16 51% wholegrain rolled edge cheese pizza": 'pizza, cheese topping, wholegrain crust',
    "big daddy's ls 16 51% wholegrain rolled edge turkey pepperoni pizza": 'pizza, turkey pepperoni topping, wholegrain crust',
    'big mac': 'double beef hamburger with cheese, lettuce and sauce',
    'big mac, without big mac sauce': 'double beef hamburger with cheese and lettuce, without sauce',
    'bilberry or blueberry, bilberry powder, dried': 'bilberry, dried',
    'biscuit 12-36 month': 'biscuit baby 12-36 month',
    'biscuit belvita ontbijtbiscuit': 'breakfast biscuit',
    'biscuit dutch shortbread sprits with choc': 'shortbread biscuit, sprits, with chocolate',
    'biscuit for children, carneval': 'biscuit for children',
    'biscuit for children, carneval prinsessa': 'biscuit for children',
    'biscuit for children, moomin': 'biscuit for children',
    'biscuit for children, tivoli': 'biscuit for children',
    'biscuit haverkick': 'biscuit, oat',
    "biscuit jaffa cakes/cake pim's": 'jaffa cake',
    'biscuit lu time out': 'biscuit, chocolate coated wafer',
    'biscuit milkbiscuit ah': 'milk biscuit',
    'biscuit milkbreak milkbiscuit': 'milk biscuit',
    'biscuit shortbread bastogne': 'biscuit, spiced caramelised shortbread',
    'biscuit spiced small kruidnoten with chocolate': 'biscuit, spiced, kruidnoten, with chocolate',
    'biscuit spiced small kruidnoten with dark choc': 'biscuit spiced small kruidnoten with dark chocolate',
    'biscuit with chocolate layer scholiertje': 'biscuit with chocolate layer',
    'biscuit with currants evergreen': 'biscuit with currant',
    'biscuit, bastogne': 'biscuit, spiced caramelised shortbread',
    'biscuit, biscuits with jam topping': 'biscuit with jam topping',
    'biscuit, cracotte': 'crispbread biscuit with sesame seed',
    'biscuit, digestive type': 'biscuit, digestive',
    'biscuit, elovena': 'oat snack biscuit',
    'biscuit, gingernut': 'biscuit, ginger nut, homemade',
    'biscuit, jyväshyvä cracotte': 'crispbread biscuit with sesame seed',
    'biscuit, savoury, dried': 'savoury rice cracker with vegetable powder',
    'biscuit, savoury, flour': 'biscuit, savoury',
    'biscuit, savoury, salted': 'savoury corn cake, salted',
    'biscuit, sweet, sweetened': 'sweet sandwich biscuit with cream and jam filling',
    'biscuit, wafers with cream': 'wafer biscuit with cream filling',
    'biscuits & snacks cheesy': 'cheese-flavoured savoury biscuit',
    'bisquit, wholegrain gluten free digestive': 'biscuit, wholegrain gluten-free, digestive type',
    'bitter': 'bitters liqueur',
    'bitter, gammel dansk bitter dram': 'bitter, herbal',
    'bitterleaf, leaf': 'bitter leaf',
    'bitterleaf, leaf, dried': 'bitter leaf, dried',
    'black currant': 'blackcurrant',
    'black currant, black currant powder, dried': 'black currant, dried',
    'blackberries with sugar product': 'blackberry, sweetened, frozen',
    'blackeye pea, dried': 'bean, blackeye, dried',
    'blancmange powder to cook pudding': 'blancmange pudding powder, dried',
    'blended spread 43% fat e.g bregott mindre': 'blended spread, 43% fat',
    'blended spread 43% fat e.g bregott mindre, enriched': 'blended spread, 43% fat, enriched',
    'blended spread 60%, voi & rypsi': 'blended spread 60%, butter and rapeseed oil',
    'blended spread 70%, creme bonjour voi & rypsiöljy, salted': 'blended spread 70%, butter and rapeseed oil, salted',
    'blended spread 70%, voi & kasviöljyt': 'blended spread 70%, butter and vegetable oils',
    'blended spread 70%, voi, salted': 'blended spread 70%, butter, salted',
    'blood pancake, kartanon': 'blood pancake, industrial',
    'blueberry mini spooner': 'cereal, ready-to-eat, frosted shredded wheat with blueberry',
    'blueberry muffin tops cereal': 'cereal, ready-to-eat, blueberry muffin flavoured',
    'boerhavia, root': 'red spiderling, root',
    'bok choy, pak-choi': 'bok choy',
    'bolognese gravy mix, knorr': 'bolognese gravy mix, dried',
    'boost plus, nutritional drink': 'nutritional supplement drink, ready-to-drink',
    'bounty bar and own brand equivalent': 'coconut bar with chocolate coating',
    "brain, calf's": 'calf, brain',
    'brazilnut': 'brazil nut',
    'brazilnut, unblanched, dried': 'brazil nut, dried',
    'bread white milk 3% fibre hönökaka': 'white milk bread, 3% fibre',
    'bread wholegrain rye rallarhalvor': 'rye bread, wholegrain',
    'bread with blood blodbröd, salted': 'blood bread, boiled with salt',
    'bread with blood paltbröd': 'blood bread, frozen',
    'bread with blood product paltbröd': 'blood bread, frozen',
    'bread with sifted rye 4% fibre rågkaka': 'sifted rye bread, 4% fibre',
    'bread, flat, flour': 'flatbread, from white wheat and whole barley flour',
    'bread, from white flour, dried': 'bread, from white flour, with dried fruit',
    'bread, from white wheat flour, dried': 'bread, from white wheat flour, with dried fruit, toasted',
    'bread, wheat and whole amarantht': 'bread, wheat and whole amaranth, with phytase',
    'bread, wholemeal': 'whole wheat bread',
    'breadcrumb, grounded crisp bread wholegrain wheat rye sugar': 'breadcrumb, ground wholegrain wheat and rye crispbread with sugar',
    'breakfast cereal cornflakes plus/1 de beste': 'breakfast cereal, cornflake',
    'breakfast cereal musli wheat oat rye barley wholegrain with fruit nut': 'breakfast cereal muesli wheat oat rye barley wholegrain with fruit nut',
    'breakfast cereal porridge olvarit fijne granen 6+ month': 'breakfast cereal porridge, fine grain, 6+ month',
    'breakfast cereal porridge olvarit tarwe en rogge 8+ month': 'breakfast cereal porridge, wheat and rye, 8+ month',
    'breakfast cereal puffed rice with sugar cocoa powder coco pop': 'breakfast cereal, puffed rice with sugar and cocoa powder',
    'breakfast cereal puffed rice with sugar rice krispy': 'breakfast cereal, puffed rice with sugar, fortified',
    'breakfast cereal rice gluten free special flake': 'breakfast cereal, rice flake, gluten-free, fortified',
    'breakfast cereal rice gluten-free special flake, enriched': 'breakfast cereal, rice flake, gluten-free, fortified',
    'breakfast cereal wheat bran wholegrain all-bran': 'breakfast cereal, wholegrain wheat bran, fortified',
    'breakfast cereal wheat bran wholegrain all-bran, enriched': 'breakfast cereal, wholegrain wheat bran, fortified',
    'breakfast cereal wheat wholegrain weetabix': 'breakfast cereal, wholegrain wheat biscuit',
    'breakfast cereal wholegrain special flake': 'breakfast cereal, wholegrain flake',
    'breakfast cereal wholegrain with sugar cheerios': 'breakfast cereal, wholegrain oat ring with sugar, fortified',
    'breakfast cereal wholegrain with sugar raisin coconut crunchy': 'breakfast cereal, wholegrain with sugar, raisin and coconut',
    'breakfast cereal, wheat, dried': 'breakfast cereal, wheat and rice, with dried fruit',
    'breakfast cereal, whole wheat, dried': 'breakfast cereal, whole wheat flake, with dried fruit and nut',
    'breakfast drink coolbest fruitontbijt': 'fruit breakfast drink',
    'breakfast drink goede morgen vifit': 'fermented milk breakfast drink',
    'breakfast, english muffin with butter': 'english muffin with butter',
    'breakfast, english muffin with cheese and sausage': 'english muffin, with cheese and sausage',
    'breakfast, english muffin with egg': 'english muffin, with egg',
    'breakfast, french toast stick': 'french toast stick',
    'breakfast, french toast with butter': 'french toast with butter',
    'breakfast, with syrup and margarine': 'breakfast platter, with syrup and margarine',
    'breast milk': 'human milk',
    'breastmilk': 'human milk',
    'broad bean, bean': 'broad bean, immature',
    'broadbean': 'broad bean',
    'broadbean, dried': 'broad bean, dried',
    'broadbean, dried, salted': 'broad bean, dried, salted',
    'broadbean, in pod': 'broad bean, in pod',
    'broadbean, salted': 'broad bean, salted',
    'broadbean, solids and liquid': 'broad bean, solids and liquid',
    'broadbean, unsalted': 'broad bean',
    'broccoli product': 'broccoli, frozen',
    'broccoli, salt': 'broccoli, salted',
    'brown french bread, with flour type or': 'brown french bread',
    'brussel sprout': 'brussels sprout',
    'brussel sprout, without addition of salt': 'brussels sprout, without addition of salt',
    'brussels sprout, automn': 'brussels sprout, autumn',
    'brussels sprout, salt': 'brussels sprout, salted',
    'brussels sprouts product': 'brussels sprout, frozen',
    'build-up powder, soup': 'nutritional supplement powder, savoury soup flavour',
    'burrito supreme with beef': 'burrito with beef, bean, cheese and sour cream',
    'burrito supreme with chicken': 'burrito with chicken, bean, cheese and sour cream',
    'burrito supreme with steak': 'burrito with steak, bean, cheese and sour cream',
    'butter biscuit, sondey butterkeks classic': 'butter biscuit',
    'butter biscuit, wholewheat': 'butter biscuit, whole wheat',
    'butter bisquit': 'butter biscuit',
    'butter oil or concentrated butter': 'butter oil',
    'butterhead lettuce, cabbage lettuce': 'butterhead lettuce',
    'button mushroom product': 'button mushroom, sliced, frozen',
    'cabbage pak-choi': 'bok choy',
    'cabbage, bok choy': 'bok choy',
    'cabbage, chinesse, salted': 'cabbage, chinese, salted',
    'cabbage, for 30 sec': 'cabbage, microwaved for 30 sec',
    'cabbage, for 45 sec': 'cabbage, microwaved for 45 sec',
    'cabbage, for 60 sec': 'cabbage, microwaved for 60 sec',
    "caesar's salad, prepacked": 'caesar salad, prepacked',
    'cake and pastry, danish pastry': 'danish pastry',
    'cake frosting with almond butter sugar toscaglasyr': 'cake frosting with almond, butter and sugar',
    'cake, cheese cake almondy': 'cheese cake, frozen',
    'cake, coffeecake': 'coffee cake with cheese',
    'cake, fruitcake': 'fruit cake, commercial',
    'cake, pudding-type, dried': 'cake mix, pudding-type, dried',
    'cake, swiss roll': 'swiss roll',
    'cake, yellow, dried': 'cake mix, yellow, dried',
    'calamus, rhizom': 'sweet flag, rhizome',
    'canavalia ensiforme, mature, dried': 'jack bean, dried',
    'canavalia ensiformis': 'jack bean',
    'canavalia ensiformis, dried': 'jack bean, dried',
    'canavalia gladiata dc, dried': 'sword bean, dried',
    'canavalia gladiata, dried': 'sword bean, dried',
    'canavalia gladiata, dried, brown': 'sword bean, brown, dried',
    'canavalia gladiata, dried, red': 'sword bean, red, dried',
    'candy': 'candy, peanut butter pieces, candy-coated',
    'candy, after eight': 'chocolate mint thins',
    'candy, after eight mint': 'candy, mint fondant, dark chocolate-coated',
    'candy, almond joy candy bar': 'milk chocolate coconut and almond bar',
    'candy, baby ruth bar': 'chocolate bar with caramel, nougat and peanuts',
    "candy, bit-o'-honey candy chew": 'candy, honey chew',
    'candy, bite': 'candy, peanut butter cup bites',
    'candy, butterfinger bar': 'chocolate-coated peanut butter crisp bar',
    'candy, caramello candy bar': 'milk chocolate bar with caramel',
    'candy, chunky bar': 'milk chocolate bar with raisins and peanuts',
    "candy, confectioner's coating": "confectioner's coating chips, butterscotch",
    'candy, crunch bar and dessert topping': 'candy, chocolate bar with crisped rice',
    'candy, fast break': 'candy, milk chocolate with peanut butter and nougat',
    'candy, golden almond solitaire': 'milk chocolate with almonds',
    'candy, goober': 'chocolate covered peanuts',
    'candy, goobers chocolate covered peanut': 'candy, chocolate covered peanut',
    'candy, gum drop': 'candy, gumdrops, sorbitol sweetened',
    'candy, hershey': 'chocolate-coated wafer bar',
    "candy, hershey's": 'chocolate-coated coconut and almond bites',
    "candy, hershey's golden almond solitaire": 'chocolate bar, milk chocolate with almond',
    "candy, hershey's milk chocolate with almond bite": 'milk chocolate with almond bites',
    "candy, hershey's pot of gold almond bar": 'milk chocolate almond bar',
    "candy, hershey's skor toffee bar": 'toffee bar with chocolate coating',
    'candy, jelly bean': 'candy, jellybean',
    'candy, kit kat wafer bar': 'chocolate-coated wafer bar',
    'candy, krackel chocolate bar': 'milk chocolate bar with crisped rice',
    'candy, mounds candy bar': 'dark chocolate coconut bar',
    'candy, mr. goodbar chocolate bar': 'milk chocolate bar with peanuts',
    'candy, oh henry! bar': 'candy, chocolate bar with fudge, caramel and peanuts',
    'candy, raisinets chocolate covered raisin': 'chocolate covered raisins',
    'candy, rolo caramels in milk chocolate': 'caramels in milk chocolate',
    'candy, skittle': 'fruit-flavoured chewy candy',
    'candy, snackfood us': 'candy, chocolate-covered whipped nougat bar',
    'candy, snickers almond bar': 'chocolate bar with nougat, caramel and almonds',
    'candy, special dark chocolate bar': 'dark chocolate, bar',
    'candy, starburst fruit chew': 'fruit chews',
    'candy, symphony milk chocolate bar': 'milk chocolate, bar',
    'candy, toblerone': 'milk chocolate with honey and almond nougat',
    'candy, tootsie roll': 'chocolate-flavour chewy roll',
    'candy, twizzler': 'liquorice-type strawberry twist',
    'candy, whatchamacallit candy bar': 'candy bar, chocolate-coated peanut butter crisp with caramel',
    'candy, york peppermint pattie': 'candy, peppermint pattie, chocolate-coated',
    'candybar': 'chocolate candy bar, average',
    'candybar bounty': 'coconut bar with chocolate coating',
    'candybar kitkat': 'chocolate-coated wafer bar',
    'candybar lion': 'wafer bar with caramel, puffed rice and chocolate',
    'candybar milky way': 'chocolate bar with whipped nougat',
    'candybar nut': 'chocolate bar with nuts and caramel',
    'candybar snicker': 'chocolate bar with caramel, nougat and peanuts',
    'cane molasses': 'molasses',
    'canola oil': 'rapeseed oil',
    'carbonated beverage, cola': 'cola',
    'carbonated drink, club soda': 'soda water',
    'carbonated drink, cola': 'cola',
    'carbonated drink, cream soda': 'cream soda',
    'carbonated drink, ginger ale': 'ginger ale',
    'carbonated drink, grape soda': 'grape soda',
    'carbonated drink, lemon-lime soda': 'lemon-lime soda',
    'carbonated drink, orange soda': 'orange soda',
    'carbonated drink, pepper type': 'pepper-type soda',
    'carbonated drink, root beer': 'root beer',
    'carbonated drink, tonic water': 'tonic water',
    'carbonated, sprite': 'carbonated lemon-lime soda, without caffeine',
    'carrot, danish': 'carrot',
    'carrot, salt': 'carrot, salted',
    'carrot, solids ans liquid, unsalted': 'carrot, solids and liquid',
    'cashew': 'cashew nut',
    "cashew, planter's, dried": 'cashew nut, roasted, dried',
    'cassava or manioc, root': 'cassava, root',
    "cat's ear flour": "cat's ear, seed flour, fat-free",
    "cat's ear flour, fat-free": "cat's ear, seed flour, fat-free",
    'catbrier, greenbrier': 'catbrier, shoot',
    'catfish b, wild': 'catfish, wild',
    'catsup': 'ketchup',
    'catsup, low sodium': 'ketchup, low sodium',
    'cauliflower product': 'cauliflower, frozen',
    'cauliflower, danish': 'cauliflower',
    "causses blue cheese, from cow's milk": 'causses blue cheese',
    'caviar, danish': 'lumpfish roe, salted',
    'celeriac, celery root': 'celeriac',
    'celeriac, turnip-rooted celery': 'celeriac',
    'celery, root': 'celeriac',
    'cereal bar with chocolate, special k bar chocolate, enriched': 'cereal bar with chocolate, fortified',
    'cereal bar with chocolate, ´s special k bar chocolate': 'cereal bar with chocolate, fortified',
    'cereal flake, plain': 'cereal flake, with sugar',
    'cereal, alpen': 'muesli cereal, ready-to-eat',
    "cereal, heritage o's": 'cereal, ready-to-eat, multigrain',
    "cereal, honest o's": 'cereal, ready-to-eat, oat ring',
    'cereal, kashi: golean crunch': 'cereal, ready-to-eat, high-protein wholegrain cluster',
    'cereal, kashi: high fibre flakes and granola': 'cereal, ready-to-eat, high fibre flake and granola',
    'cereal, kashi: honey puffed': 'cereal, ready-to-eat, honey puffed wholegrain',
    'cereal, kashi: puffed': 'cereal, ready-to-eat, puffed wholegrain',
    "cereal, o's": 'cereal, ready-to-eat, toasted oat',
    "cereal, oat o's": 'cereal, ready-to-eat, toasted oat',
    "cereal, whole o's": 'cereal, ready-to-eat, wholegrain oat ring',
    "cereals ready-to-eat, quaker, cap'n crunch": 'cereal, ready-to-eat, sweetened corn and oat',
    'champignon': 'button mushroom',
    'champignon, fried': 'button mushroom, fried',
    'chantarelle': 'chanterelle',
    'chanterelle or girolle mushroom': 'chanterelle',
    'chanterelle, kokt': 'chanterelle, boiled',
    'chard': 'swiss chard',
    'cheerios': 'cereal, ready-to-eat, toasted wholegrain oat ring',
    'cheese 20+ leidse with cumin/fries clove': 'cheese, leidse with cumin or frisian with clove, 20+',
    'cheese 30+ age 4-7 mth': 'cheese, 30+, aged 4-7 month',
    'cheese 40+ leiden with cumin/fries clove': 'cheese 40+ leiden with cumin/frisian clove',
    'cheese bluefort': 'blue cheese',
    'cheese cake with cottage cheese ostkaka': 'cheesecake with cottage cheese',
    'cheese cream soft boursin': 'cheese, cream, soft, with garlic and herbs',
    'cheese cream soft mon chou': 'cheese, cream, soft',
    'cheese cream soft paturain': 'cheese, cream, soft, with herbs',
    'cheese dutch in swiss-style': 'cheese, swiss-style 45+',
    'cheese emmenthaler': 'cheese, emmental',
    "cheese goat's milk chèvre": "chèvre cheese, from goat's milk, 25% fat",
    'cheese hard parmesan': 'cheese, parmesan',
    "cheese mozzarella made from cow's milk": "cheese, mozzarella, from cow's milk",
    'cheese old amsterdam': 'cheese gouda, extra mature',
    'cheese pasty, vaasan kotiuunin': 'cheese pasty',
    "cheese pizza, thin 'n crispy crust": 'cheese pizza, thin crispy crust',
    'cheese pizza, ultimate deep dish crust': 'cheese pizza, deep dish crust',
    'cheese soufflé': 'cheese souffle',
    'cheese soup, vegetable bouillon, flour': 'cheese soup with vegetable bouillon, thickened with wheat flour',
    'cheese spread 15+ balans eru': 'cheese spread 15+',
    'cheese spread 60+ kiri': 'cheese spread 60+',
    'cheese spread kids eru': 'cheese spread, for children',
    'cheese spread, american or cheddar cheese base': 'cheese spread, american or cheddar cheese base, low-fat',
    'cheese, american cheddar': 'imitation cheese, american cheddar',
    "cheese, cheddar, from cow's milk": 'cheese, cheddar',
    'cheese, cottage, small curd': 'cheese, cottage, small curd, 4% milkfat',
    'cheese, dried': 'cheese, queso seco',
    "cheese, emmental, from cow's milk": 'cheese, emmental',
    'cheese, fetta': 'cheese, feta',
    'cheese, goat, 20-25%': 'cheese, goat, in brine, 20-25%',
    "cheese, gorgonzola, from cow's milk": 'cheese, gorgonzola',
    "cheese, mozzarella, from buffala's milk": "cheese, mozzarella, from buffalo's milk",
    "cheese, mozzarella, from cow's milk": 'cheese, mozzarella',
    'cheese, oltermanni rypsi': 'cheese with rapeseed oil, 24% fat',
    "cheese, parmesan, from cow's milk": 'cheese, parmesan',
    "cheese, pecorino, from ewe's milk": 'cheese, pecorino',
    'cheese, polar': 'cheese, semi-hard, 20% fat',
    'cheese, polar täyteläinen': 'cheese, semi-hard, full-flavoured, 28% fat',
    'cheese, port de salut': 'cheese, port salut',
    'cheese, raclette': "cheese, raclette, from cow's milk",
    'cheese, reblochon': "cheese, reblochon, from cow's milk",
    "cheese, reblochon, from cow's milk": 'cheese, reblochon',
    'cheese, soft, around': 'cheese, soft, around 6% fat',
    'cheeseprod with veg fat kees oud': 'cheese product with vegetable fat, mature',
    'cheeseproduct with veg fat kees jong belegen': 'cheese product with vegetable fat, medium-matured',
    'cheez whiz light process cheese product': 'process cheese product, low-fat',
    'cheez whiz process cheese sauce': 'process cheese sauce',
    'chewing gum without sugar': 'chewing gum, sugar-free',
    'chewing gum, sugarless': 'chewing gum, sugar-free',
    'chick pea': 'chickpea',
    'chick pea flour': 'chickpea flour',
    'chick pea patty, falafel, fried': 'falafel, fried',
    'chick pea, dried': 'chickpea, dried',
    'chick peas leblebi turkish': 'chickpea, roasted',
    'chick-n-strip': 'chicken strip, breaded, fried',
    'chicken filet sandwich, with lettuce': 'chicken fillet sandwich, grilled, with lettuce and tomato',
    "chicken finger, from kid's menu": 'chicken finger, breaded, fried',
    "chicken finger, from kids' menu": 'chicken finger, breaded, fried',
    'chicken mcnugget': 'chicken nugget, breaded, fried',
    'chicken meatball, kariniemen': 'chicken meatball',
    'chicken nugget, chicken mcnugget, fried': 'chicken nugget, deep fried',
    'chicken nugget, star shaped': 'chicken nugget, breaded, fried',
    'chicken parmigiana without pasta': 'chicken parmesan without pasta',
    'chicken rice, with pea, fried': 'chicken fried rice with pea, onion and egg',
    'chicken schnitzel, kariniemi': 'chicken schnitzel, breaded, oven-baked',
    "chicken tender, from kids' menu": 'chicken tender, breaded, fried',
    "chicken tenderloin platter, from kid's menu, fried": 'chicken tenderloin platter, fried',
    'chicken, brest, fried': 'chicken, breast, fried',
    'chicken, no broth': 'chicken, canned, without broth',
    'chicken, product': 'chicken, canned roast meat with seasoning',
    'chicken, tsukune': 'chicken meatball, tsukune',
    'chicken, white race': 'chicken, white breed, meat and skin',
    'chickory, raddichio': 'radicchio',
    'child formula, nutrition': 'child formula, ready-to-feed',
    'childrens food lasagna veg': "children's food, vegetable lasagna, canned",
    'childrens food spaghetti bolognese': "children's food, spaghetti bolognese, canned",
    'chinese cabbage, pak choï ou pak choy': 'bok choy',
    'chinese cabbageor bok choï brede, rods and leaf': 'bok choy, stems and leaf',
    'chinese dish, chicken and vegetable': 'chicken with vegetable, chinese-style',
    'chinese dish, lo mein': 'vegetable lo mein, without meat',
    'chinese dish, orange chicken': 'orange chicken',
    'chinese dish, rice without meat, fried': 'fried rice without meat, restaurant prepared',
    'chinese dish, shrimp and vegetable': 'shrimp with vegetable, restaurant-prepared',
    'chinese dish, sweet and sour chicken': 'sweet and sour chicken',
    'chinese dish, sweet and sour pork': 'sweet and sour pork',
    'chinese pear, namphung': 'chinese pear, peeled',
    'chinese privet flour': 'chinese privet, seed flour, fat-free',
    'chinese privet flour, fat-free': 'chinese privet, seed flour, fat-free',
    'chip': 'oven chips, baked from frozen',
    'chip, fried': 'chips, pre-fried, frozen',
    'chip, vegetable pre, fried': 'vegetable chips, pre-fried, frozen',
    'chip, wheat with bacontaste': 'wheat snack with bacon flavour',
    'cho cho': 'chayote',
    'cho-cho-marrow': 'chayote',
    'chocolate and cappucino filled biscuit, catago medaillon cappuccino': 'biscuit filled with chocolate and cappuccino',
    "chocolate and milk bar, hershey's": 'chocolate and milk bar',
    'chocolate bar filled kinder': 'chocolate bar with milk filling',
    'chocolate biscuit with vanilla filling cookie': 'chocolate biscuit, with vanilla filling',
    'chocolate chip cookie, jyväshyvä suklaapisara': 'chocolate chip cookie, industrial',
    'chocolate chip cookie, suklaapisara': 'chocolate chip cookie, industrial',
    'chocolate filled biscuit with meringue, catago chocobit': 'biscuit with meringue filled with chocolate',
    'chocolate filled with caramel rolo': 'chocolate filled with caramel',
    'chocolate marshmallow matey': 'cereal, ready-to-eat, chocolate oat with marshmallow',
    'chocolate swiss roll, naturally gluten-free': 'chocolate swiss roll, gluten-free',
    'chocolate, dried': 'hot wheat cereal, chocolate flavour, dried',
    'chocolate, hot': 'hot chocolate mix, instant powder',
    'chocolate, with water, unsalted': 'hot wheat cereal, chocolate flavour, prepared with water',
    'chokeberry, chokeberry powder, dried': 'chokeberry powder, dried',
    'choko': 'chayote',
    'chole, fa - beef, lean': 'beef, eye of round roast/steak, lean',
    'chop': 'chop, average of veal, pork and lamb',
    'chorizo pork sausage, johnsonville': 'chorizo pork sausage',
    'chorizo pork sausage, mixed brand': 'chorizo pork sausage',
    'christmas crunch': 'cereal, ready-to-eat, sweetened corn and oat, seasonal',
    'cider vol. %': 'cider, 1 vol%',
    'cider, half': 'cider, half-dry',
    'cider, traditionnal': 'cider, traditional',
    'cinnamon toaster': 'cereal, ready-to-eat, cinnamon wheat and rice square',
    'clam, pippies': 'clam, pipi',
    'coco-roo': 'cereal, ready-to-eat, cocoa puffed corn',
    'cocoa dyno-bite': 'cereal, ready-to-eat, cocoa crisp rice',
    'cocoa pebble': 'cereal, ready-to-eat, cocoa flavoured crisp rice',
    'cocoa powder sweetened nesquik': 'cocoa powder, sweetened',
    'cocoa product powder ovomaltine': 'malted cocoa drink powder, fortified',
    'cocoa, instant, dried': 'cocoa, instant, with milk, powder',
    'coconut fat or oil': 'coconut oil',
    'coconut, coconut water': 'coconut water',
    'cocyam, leaf': 'cocoyam, leaf',
    'cocyam, leaf, dried': 'cocoyam, leaf, dried',
    'cod, filet': 'cod, fillet',
    'coffee whitener, powdered': 'coffee whitener powder',
    'coffeemate, whitener powder': 'coffee whitener powder',
    'cognac or brandy vol. %': 'cognac or brandy, 40 vol%',
    'cold sub on white bread with lettuce and tomato': 'cold cut sub on white bread with lettuce and tomato',
    'coley, flesh': 'saithe, flesh',
    'colombo papaya, flesh, mature': 'papaya, flesh',
    'colossal crunch': 'cereal, ready-to-eat, sweetened corn and oat puff',
    'commom bean, dried': 'common bean, dried',
    'common wheat, bread': 'bagel, from common wheat',
    'common wheat, instant, dried': 'instant chinese noodle, non-fried, cup, dried',
    'common wheat, outer wheat jiaozi dough': 'jiaozi dumpling wrapper, steamed wheat dough',
    'common wheat, outer wheat shumai dough': 'shumai dumpling wrapper, steamed wheat dough',
    'common wheat, premixed flour for food': 'premixed flour for fried food',
    'common wheat, somen and hiyamugi, dried': 'somen and hiyamugi noodle, dried, boiled',
    'common wheat, udon': 'udon noodle, boiled',
    'common wheat, udon, dried': 'udon noodle, dried, boiled',
    'common wheat, wholegrain flour': 'whole wheat flour',
    'complan powder, sweet': 'nutritional supplement drink powder, sweet',
    'condensed milk, no added sugar': 'evaporated milk, whole',
    'condensed whole milk, sweetened': 'milk, condensed, sweetened',
    'continental mill, krusteaz almond poppyseed muffin mix, dried': 'almond poppyseed muffin mix, dried',
    'cookie, different variety': 'cookie, assorted',
    'cookie, vanilla sandwich with cream filling': 'cookie, vanilla sandwich with creme filling',
    'cooking fat liq vlees&jus': 'cooking fat, liquid, for meat and gravy',
    'cooking fat solid 97% fat >17 g sat fa, salted': 'cooking fat, solid, 97% fat, >17 g saturated fatty acid, salted',
    'coquille, scallop shell st. jacque': 'scallop, coquille saint-jacques',
    'cordial base juice': 'cordial base, 25% citrus fruit juice',
    'cordial juice': 'cordial, 25% citrus fruit juice, diluted',
    'corn burst': 'cereal, ready-to-eat, frosted corn flake',
    'corn cobs product': 'corn cob, frozen',
    'corn flour, wholegrain': 'corn flour, wholegrain, blue corn',
    'corn oil, maize oil': 'corn oil',
    'corn product': 'corn kernel, frozen',
    "corn salad, lamb's lettuce": "lamb's lettuce",
    'corn, combined variety': 'corn, whole kernel, boiled',
    'corn, combined variety, dried': 'corn, dried',
    'corn, corn flour': 'corn flour (ground corn), white kernel',
    'corn, corn grits': 'corn grits, white kernel',
    'corn, corn meal': 'cornmeal, white kernel',
    'corn, cultivar: cuzco, salted': 'corn kernel, cuzco, oil-roasted and salted',
    'corn, dmr-esr-w variety, dried': 'corn, dried',
    'corn, gbaévè variety, dried': 'corn, dried',
    'corn, gnonli variety, dried': 'corn, dried',
    'corn, gougba variety, dried': 'corn, dried',
    'corn, salted': 'corn kernel, oil-roasted and salted',
    'corn, tzpb-sr variety, dried': 'corn, dried',
    'cos or romaine lettuce': 'romaine lettuce',
    'cotton oil': 'cottonseed oil',
    'country-style bread, french bread': 'french bread, country-style, multigrain and/or seeds',
    'courgette': 'zucchini',
    'cowberry': 'lingonberry',
    'crab blue': 'crab, blue',
    'cracker vitalu': 'cracker, fortified',
    'cracker vitalu with calcium': 'cracker, fortified with calcium',
    'cracker vitalu with calcium, enriched': 'cracker, fortified with calcium',
    'cracker, fortified, enriched': 'cracker, fortified with calcium',
    'cracker, matz': 'cracker, matzo',
    'cracker, wheat thin': 'cracker, baked wheat, thin',
    'crayfish freshwater': 'freshwater crayfish',
    'cream and similar, fat content unknownthick or semi-thick': 'cream and similar, fat content unknown, thick or semi-thick',
    'cream cheese, block': 'cheese, cream, full-fat',
    'cream of wheat, 1 minute cook time, dried': 'farina, 1 minute cook time, dried',
    'cream of wheat, 2 1/2 minute cook time': 'farina, quick cooking, cooked with water',
    'cream of wheat, 2 1/2 minute cook time, dried': 'farina, dried',
    'cream of wheat, dried': 'farina, dried',
    'cream of wheat, instant, dried': 'farina, instant, dried',
    'cream of wheat, salted': 'farina, cooked with water, salted',
    'cream type product finesse voor koken': 'cream-type product for cooking',
    'cream, substitute, dried': 'cream substitute, flavoured, dried',
    'creamy dressing, made with sour cream and/or buttermilk and oil': 'creamy dressing, made with sour cream and/or buttermilk and oil, reduced calorie',
    'crepe filled chocolate and cerels ball, prepacked': 'crepe filled with chocolate and cereal balls',
    'cripsbakes dutch low sodium': 'crispbakes dutch low sodium',
    'crisp bread sandwich wholegrain rye diff. filling': 'crisp bread sandwich, wholegrain rye, various fillings',
    'crisp bread wheat sandwich diff. filling': 'crisp bread sandwich, wheat, various fillings',
    'crisp bread wheat with poppy seed': 'crisp bread, wheat with poppy seed, 6% fibre',
    'crisp bread wholegrain rye 15% fibre flatbröd': 'crisp bread, wholegrain rye, 15% fibre',
    'crisp bread wholegrain rye 18% fibre mörkt': 'crisp bread, wholegrain rye, dark, 18% fibre',
    'crisp bread wholegrain rye 20% fiber e.g sport': 'crisp bread, wholegrain rye, 20% fibre',
    'crisp bread wholegrain rye wheat corn with sour dough 15% fibre spisbröd': 'sourdough crisp bread, wholegrain rye, wheat and corn, 15% fibre',
    'crisp bread wholegrain wheat with poppy seeds fibre 5.5%': 'crisp bread, wholegrain wheat with poppy seed, 5.5% fibre',
    'crisp rice': 'crispy rice',
    'crispbread gluten-free fette croccanti schar': 'crispbread, gluten-free',
    'crispbread sandwich wasa': 'crispbread sandwich',
    'crispbread, whole rye bread': 'rye crispbread, wholegrain',
    'crispy chicken filet sandwich, with lettuce and mayonnaise': 'crispy chicken fillet sandwich, with lettuce and mayonnaise',
    'crispy hexagon': 'cereal, ready-to-eat, corn and oat hexagon',
    "croissan'wich with egg and cheese": 'croissant breakfast sandwich with egg and cheese',
    "croissan'wich with sausage and cheese": 'croissant breakfast sandwich with sausage and cheese',
    "croissan'wich with sausage, egg and cheese": 'croissant breakfast sandwich with sausage, egg and cheese',
    'crown, sprout': 'crown daisy, sprout',
    'crunchy bran': 'cereal, ready-to-eat, crunchy wheat bran',
    'crustecean, crayfish, farmed': 'crustacean, crayfish',
    'crustecean, spiny lobster': 'crustacean, spiny lobster',
    'crêpe with button mushroom filling heated product': 'crêpe with button mushroom filling, frozen, heated',
    'cumquat': 'kumquat',
    'curly dock, yellow dock, leaf': 'curly dock, leaf',
    "curry sauce, uncle ben's": 'curry sauce, jarred',
    'curry tree, leaf': 'curry leaf',
    'cyca, tuber': 'cycad, tuber',
    'cytosport, muscle milk': 'milk protein drink, ready-to-drink',
    'daim cake': 'almond cake with chocolate and toffee',
    'dairy drink campina fruitmelk': 'dairy drink with fruit',
    'dandelion greens': 'dandelion, leaf',
    'danish pastry, square': 'danish pastry',
    'dannon, water': 'bottled water, non-carbonated, with fluoride',
    'dark chocolate, more than': 'dark chocolate, more than 40% cocoa, for cooking',
    'desert horse purslane, just leaf': 'desert horse purslane, mature leaf',
    'desert horsepurslane, leaf': 'desert horse purslane, leaf',
    'dessert date, flesh, dried': 'desert date, flesh, dried',
    'dessert date, ndofane1, dried': 'desert date, ndofane1, dried',
    'dessert, fudgesicle bar': 'dessert, frozen, fudge bar',
    'dessert, popsicle pop, orange': 'ice lolly, sugar-free',
    'dextrose tablets non': 'dextrose tablets, non-fortified',
    'dfr modifast intensive milkshake': 'meal replacement milkshake',
    'dip mix powder diff. flavour': 'dip mix powder, various flavours',
    "dip, frito's": 'bean dip, original flavour',
    'dip, hummus': 'hummus, commercial',
    'distilled alcoholic beverage, whisky': 'whisky',
    'dolichos bean, cofee-brown': 'dolichos bean, coffee-brown',
    'dolichos lablab bean': 'hyacinth bean',
    'double quarter pounder with cheese': 'cheeseburger, large double beef patty',
    'double whopper, with cheese': 'double hamburger, flame-grilled, with cheese',
    'double whopper, without cheese': 'double hamburger, flame-grilled, without cheese',
    'doughnut roasting fat 100%, sunnuntai': 'doughnut frying fat, 100% fat',
    'drink mix, oat, dried': 'sports drink mix, orange flavour, dried',
    'drink soy groeidrink 1-3': 'soy growing-up drink, 1-3 years',
    'drinking chocolate, powder': 'drinking chocolate powder',
    'drumstick or moringa pod': 'drumstick pod',
    'drumstick or moringaleaf': 'drumstick leaf',
    'drumstick, just leaf': 'drumstick, mature leaf',
    'duchess potato': 'duchesse potato',
    'dulce de leche or confiture de lait': 'dulce de leche',
    'durain': 'durian',
    'durain, chanee': 'durian, chanee',
    'durum wheat pre, wholegrain': 'durum wheat, pre-cooked, wholegrain',
    'dutch brand loaf, chicken': 'dutch brand loaf, chicken, pork and beef',
    'eas soy protein powder': 'soy protein powder',
    'eas whey protein powder': 'whey protein powder',
    'edamame bean, shelled': 'edamame, shelled',
    'edamame beans parbolied without pod': 'edamame bean, parboiled, shelled',
    'edamame, hearty brand': 'edamame',
    'edible burdock, root': 'burdock root',
    'effect drink': 'probiotic milk drink, fat-free',
    'effect drink, fat-free': 'probiotic milk drink, fat-free',
    'egg mcmuffin': 'breakfast muffin sandwich with egg, cheese and ham',
    'egg powder chicken': 'egg, dried',
    'egg white chicken': 'egg white',
    'egg whole chicken av': 'egg, chicken, whole, boiled',
    'egg yolk chicken': 'egg yolk',
    'egg, hard': 'egg, hard-boiled',
    'egg, soft': 'egg, soft-boiled',
    'elongate glass-perchlet, large size': 'elongate glass-perchlet',
    'elovena muru, broad bean product': 'broad bean crumb product',
    'emblic': 'amla',
    'emmentaler cheese': "cheese, emmental, from cow's milk",
    'energy drink golden power/bullit/freeway': 'energy drink',
    'energy drink red bull': 'energy drink, with sugar',
    'energy drink red bull sugarfree': 'energy drink, sugar free',
    'energy drink, amp': 'energy drink, with sugar',
    'energy drink, full throttle': 'energy drink, with sugar',
    'energy drink, monster': 'energy drink, with sugar',
    'energy drink, monster with vitamins c, enriched': 'energy drink, with sugar, enriched',
    'energy drink, rockstar': 'energy drink, with sugar',
    'energy drink, vault': 'energy drink, with sugar',
    'energy drink, vault zero': 'energy drink, sugar free',
    'english muffin, from white flour': 'muffin, english',
    'english muffin, from white flour, with dried fruit, toasted': 'muffin, english, with dried fruit, toasted',
    'ensure plus, ready-to-drink': 'nutritional supplement drink, ready-to-drink',
    'ensure, nutritional shake': 'nutritional shake, ready-to-drink',
    'escarole or endive': 'escarole',
    'european pilchard or sardine': 'european pilchard',
    'european rowan': 'rowanberry',
    'exotic fruit, flesh': 'mango/papaya, without skin',
    'extravaganzza feast pizza, hand-tossed crust': 'pizza with meat and vegetable topping, hand-tossed crust',
    'falafel chickpea croquette, fried': 'falafel, fried',
    'familia': 'swiss-style muesli, ready-to-eat',
    'farina, assorted brands including cream of wheat, dried': 'farina, enriched, quick cooking, dried',
    'farina, assorted brands including cream of wheat, dried, enriched': 'farina, enriched, quick cooking, dried',
    'farina, assorted brands including cream of wheat, unsalted': 'farina, enriched, quick cooking, cooked',
    'fat animal': 'animal fat',
    'fat blend liquid 80% fat arla smör- och rapsolja': 'fat blend, liquid, 80% fat, butter and rapeseed oil',
    'fat blend liquid 80% fat ica raps- och smörolja': 'fat blend, liquid, 80% fat, rapeseed and butter oil',
    'fat blend liquid 80% fat ica raps- och smörolja, enriched': 'fat blend, liquid, 80% fat, rapeseed and butter oil, enriched',
    'fat for gravy smeltjus': 'fat for gravy',
    'fenugreek seed, 10 minute': 'fenugreek seed, boiled 10 minutes',
    'fenugreek seed, 15 minute': 'fenugreek seed, boiled 15 minutes',
    'fenugreek seed, 5 minute': 'fenugreek seed, boiled 5 minutes',
    "feta-type cheese from ewe's milk, and spice": "feta-type cheese from ewe's milk, in oil and spice",
    'filberts or hazelnut, dried': 'hazelnut, blanched, dried',
    'filberts or hazelnut, unblanched': 'hazelnut, dry roasted',
    'filet-o-fish': 'fish fillet sandwich with cheese and tartar sauce',
    'filet-o-fish, without tartar sauce': 'fish fillet sandwich with cheese, without tartar sauce',
    'fine rye bread flour': 'fine rye bread, from fine rye and dark wheat flour',
    'finnish cheese': 'finnish oven cheese',
    'fish fat >5 g fat': 'fish, fatty, >5 g fat, raw',
    'fish finger, and pre, fried': 'fish finger, breaded and pre-fried',
    'fish or seafood au gratin, intended to be': 'fish or seafood au gratin, intended to be cooked',
    'fish soup, fish': 'fish soup, fresh fish, milk-based',
    'fish stew brasilian style with couconut milk tomato pulp bell pepper': 'fish stew, brazilian style, with coconut milk, tomato pulp and bell pepper',
    'fish stock paste or powder large-scale': 'fish stock, paste or powder',
    'fish stock paste or powder with large-scale, low-salt': 'fish stock paste or powder, reduced salt',
    'fish stock paste or powder with reduced salt large-scale': 'fish stock paste or powder, reduced salt',
    'fish, cooking cream': 'fish, oven-baked in cooking cream',
    'fish, filet': 'fish, fillet',
    'fish, gefiltefish': 'gefilte fish',
    'fish, hand': 'fish fillet, hand-battered, fried',
    'fish, japanese smelt': 'fish, japanese smelt, simmered in soy and syrup',
    'flax seed': 'flaxseed, roasted',
    'flaxseed or flaxseed': 'flaxseed',
    'flesh, rhus macowanii': 'karee, flesh',
    'flounder, filet': 'flounder, fillet',
    'flour oat wholemeal': 'oat flour, wholegrain',
    'flour omelette with nutella': 'flour omelette with chocolate hazelnut spread',
    'flour rice with vit b1': 'rice flour, fortified with vitamin b1',
    'flour rice with vit b1, enriched': 'rice flour, fortified with vitamin b1',
    'flour rye wholemeal': 'rye flour, wholegrain',
    'flour spelt wholemeal': 'spelt flour, wholegrain',
    'flour wheat wholemeal': 'whole wheat flour',
    'flour, wholegrain oat': 'oat flour, wholegrain',
    'fluted guord, dried': 'fluted gourd, dried',
    'formulated bar, luna bar': 'formulated nutrition bar, chocolate nut',
    'formulated bar, power bar': 'formulated nutrition bar, chocolate',
    'formulated bar, slim-fast optima meal bar': 'meal replacement bar, milk chocolate peanut',
    'formulated bar, snackfood us': 'formulated bar, chocolate almond snack',
    'formulated bar, south beach protein bar': 'formulated protein bar',
    'formulated bar, zone perfect crunch bar': 'formulated nutrition bar, mixed flavours',
    'fragrant manjack, wihout': 'fragrant manjack, without seeds',
    'frankfurt sausage': 'frankfurter',
    'frankfurt sausage mix, with frensh fry, fried': 'frankfurt sausage mix, with french fries, fried',
    'frankfurter sausage': 'frankfurter',
    'free singles american process cheese product, fat-free': 'american process cheese product, fat-free',
    'french bread flour': 'french bread',
    'french fries and onion': 'french fries with sauces and onion',
    "french's yellow mustard": 'mustard, yellow',
    "friday's shrimp": 'shrimp, breaded',
    'frosted mini spooner': 'cereal, ready-to-eat, frosted shredded wheat',
    'frosty dairy dessert': 'frozen dairy dessert',
    "fruit 'mimusops obovata', large": 'mimusops obovata fruit, large',
    'fruit av excl citrus': 'fruit, excluding citrus',
    'fruit av including citrus': 'fruit, including citrus',
    'fruit drink conc diluted 1 to': 'fruit drink concentrate, diluted 1 to 7',
    'fruit drink conc with 45-50 mg vit c': 'fruit drink concentrate, with 45-50 mg vitamin c',
    'fruit drink conc with sugar and sweetener 10-15 g cho': 'fruit drink concentrate with sugar and sweetener, 10-15 g carbohydrate',
    'fruit drink conc with sugar and sweetener 30-35 g cho': 'fruit drink concentrate, with sugar and sweetener, 30-35 g carbohydrate',
    'fruit drink concentrate karvan cevitam': 'fruit drink concentrate, undiluted',
    'fruit drink concentrate roosvicee fruitkracht ferro': 'fruit drink concentrate with iron',
    'fruit drink concentrate roosvicee fruitkracht pruimen': 'fruit drink concentrate, prune',
    'fruit drink concentrate with sweetener karvan cevitam': 'fruit drink concentrate with sweetener',
    'fruit drink with dairy taksi with sugar': 'fruit drink with dairy, with sugar',
    'fruit drink with dairy taksi with sweetener': 'fruit drink with dairy, with sweetener',
    'fruit drink, hevi shot': 'fruit drink, apple-carrot-strawberry',
    'fruit drink, hyvää päivää, sweetened': 'fruit drink, citrus, partially artificially sweetened',
    'fruit drink, mehukatti fruit snack': 'fruit drink',
    'fruit drink, proviva': 'fruit drink with probiotics',
    'fruit flavored drink containing less than, with high vitamin c, juice': 'fruit flavored drink, less than 3% fruit juice, with vitamin C',
    'fruit flavored drink, less than, juice': 'fruit flavored drink, less than 3% juice',
    'fruit juice concentrated diluted 1 to': 'fruit juice concentrate, diluted 1 to 13',
    'fruit juice drink dubbeldrank': 'juice drink, fruit',
    'fruit juice drink roosvicee multivit': 'fruit juice drink, multivitamin',
    'fruit juice drink with sweetener 5-<8 g cho': 'fruit juice drink with sweetener, 5-<8 g carbohydrate',
    'fruit juice drink, greater than': 'fruit juice drink, greater than 3% fruit juice',
    'fruit juice smoothie, bolthouse farm': 'fruit juice smoothie, berry',
    'fruit juice smoothie, naked juice': 'fruit juice smoothie, mixed fruit',
    'fruit juice smoothie, odwalla': 'fruit juice smoothie, mixed fruit',
    "fruit pod 'sarcostemma viminale', large": 'sarcostemma viminale fruit pod, large',
    'fruit soup pasteurized rte mixed fruit, dried': 'fruit soup, ready-to-eat, with mixed dried fruit',
    'fruit soup rte mixed fruit, dried, enriched': 'fruit soup, ready-to-eat, with mixed dried fruit, enriched',
    'fruity dyno-bite': 'cereal, ready-to-eat, fruit flavoured crisp rice',
    'fruity pebble': 'cereal, ready-to-eat, fruit flavoured crisp rice',
    'frying fat horeca': 'frying fat, catering',
    'game meat, native': 'agutuk, meat-caribou',
    'game meat, native, dried': 'bearded seal (oogruk), air-dried',
    'gammon deboned': 'gammon ham, boiled and deboned',
    'george weston bakery, brownberry sage and onion stuffing mix, dried': 'stuffing mix, sage and onion, dried',
    'george weston bakery, thomas english muffin': 'muffin, english',
    'gerolsteiner brunnen gmbh & co. kg, gerolsteiner naturally sparkling mineral water': 'sparkling mineral water',
    'gerolsteiner naturally sparkling mineral water': 'sparkling mineral water',
    'gex blue cheese, or jura blue cheese or septmoncel blue cheese': 'gex blue cheese',
    'ghee butter': 'clarified butter',
    'ghee, clarified butter': 'clarified butter',
    'gin vol. %': 'gin, 40 vol%',
    'glaceau vitamin water': 'vitamin water, fruit punch flavoured, fortified',
    'glaceau vitamin water, enriched': 'vitamin water, fruit punch flavoured, fortified',
    'glass perchlet, large size': 'glass perchlet, large',
    'glass perchlet, medium size': 'glass perchlet',
    'glass perchlet, small size': 'glass perchlet, small',
    'glucose liquid, bp': 'liquid glucose, pharmaceutical grade',
    'gluten-free, french dinner roll': 'dinner roll, french, gluten-free',
    'gluten-free, soft & delicious white sandwich bread': 'white sandwich bread, gluten-free, soft',
    'gluten-free, soft & hearty wholegrain bread': 'wholegrain bread, gluten-free, soft',
    'glutino, gluten-free cookie': 'cookie, gluten-free, chocolate vanilla creme',
    'glutino, gluten-free wafer': 'wafer, gluten-free, lemon flavoured',
    'goblet of ice cream, coffe or chocolate ice cream topped with whipped cream': 'goblet of ice cream, coffee or chocolate, topped with whipped cream',
    'golden crisp': 'cereal, ready-to-eat, sugar-coated puffed wheat',
    'golden layer buttermilk biscuit, refrigerated dough': 'layered buttermilk biscuit, refrigerated dough',
    'golden puff': 'cereal, ready-to-eat, sugar-coated puffed wheat',
    'gorgonzola': 'cheese, gorgonzola',
    'goutweed, ground elder, leaf': 'ground elder, leaf',
    'grain, corn flour': 'corn flour',
    'grand, buttermilk biscuit': 'buttermilk biscuit, large, refrigerated dough',
    'granola with nuts and/or only': 'granola with nuts and/or seeds only',
    'grape oil': 'grape seed oil',
    'grape-nuts cereal': 'cereal, ready-to-eat, baked wheat and barley nugget',
    'grape-nuts flake': 'cereal, ready-to-eat, baked wheat and barley flake',
    'gratin with cabbage kalpudding': 'gratin with cabbage',
    'gratin with herring potatoes egg milk sillpudding, salted': 'gratin with salted herring, potato, egg and milk',
    'gravy, from powder with water, dried': 'gravy, prepared from powder with water',
    'great grain': 'cereal, ready-to-eat, wholegrain with raisin, date and pecan',
    'great grains crunchy pecan cereal': 'cereal, ready-to-eat, wholegrain with crunchy pecan',
    'greek sallad with feta cheese': 'greek salad with feta cheese',
    'greek-style marinated mushroom, prepacked': 'greek-style marinated mushroom',
    'green beans product': 'green bean, frozen',
    'green chilisauce basbaas cagaar': 'green chilli sauce',
    'green gram': 'mung bean',
    'green peas product': 'green pea, frozen',
    'grooseberry, wild, green': 'gooseberry, green',
    'grooseberry, wild, purple': 'gooseberry, purple',
    'ground nut': 'peanut',
    'ground turkey, pan- crumble': 'ground turkey, 93% lean, pan-broiled crumble',
    'groundcherry': 'cape gooseberry',
    'guava, common': 'guava',
    'guinea, flesh': 'guinea hen, flesh',
    'guinea, meat and skin': 'guinea hen, meat and skin',
    'guinea, meat only': 'guinea hen, meat only',
    'gum': 'seed gum (locust bean, guar)',
    'halibut wild': 'halibut, wild',
    'halibut, weighed and skin': 'halibut, weighed with bones and skin',
    'ham shoulder medium fat': 'ham shoulder, medium fat, boiled',
    'ham stock paste or powder large-scale': 'ham stock, paste or powder',
    'hamburger double meat with bun cheese pickled cucumber from restaurant': 'hamburger, double meat, with bun, cheese and pickled cucumber',
    'hamburger, double whopper': 'double hamburger, flame-grilled, takeaway',
    'hamburger, large': 'hamburger, single patty, with condiments',
    'hamburger, mcfeast': 'hamburger with cheese, lettuce and tomato, takeaway',
    'hamburger, whopper': 'hamburger, flame-grilled, takeaway',
    'hand': 'fish fillet, hand-battered, fried',
    'hard biscuit': 'biscuit, hard, with almond',
    'hard margarine 80%, sunnuntai': 'hard margarine, 80% fat, for baking',
    'hard wheat, semolina': 'durum wheat semolina',
    'hare, filet': 'hare, fillet',
    'haricot bambara, paste of cowpea and okra powder with egg and onion': 'haricot bambara, fried paste of cowpea and okra powder with egg and onion',
    'havregurt, oat product, enriched': 'oat-based yoghurt alternative, flavoured',
    'hazelnut, ham': 'hazelnut, ham variety',
    'hazelnuts or filbert': 'hazelnut',
    'hemp seed, without hull': 'hemp, hulled',
    'hemp seeds with husk': 'hemp, whole',
    'hemp seeds without husk': 'hemp, hulled',
    'hempseed oil': 'hemp oil',
    'herb paste boemboe': 'indonesian spice paste',
    'herb, cream cheese': 'herb spread with cream cheese and garlic',
    "hi-c flashin' fruit punch": 'drink, fruit punch',
    'homestyle chicken fillet sandwich': 'breaded chicken fillet sandwich, takeaway',
    'honey bunches of oat': 'cereal, ready-to-eat, honey roasted oat and corn flake',
    'honey bunches of oats with vanilla bunch': 'cereal, ready-to-eat, oat and corn flake with vanilla cluster',
    'honey buzzer': 'cereal, ready-to-eat, honey-sweetened corn and oat',
    'honey graham life cereal': 'cereal, ready-to-eat, honey graham oat and corn square',
    'honey graham oh!': 'cereal, ready-to-eat, honey graham corn and oat',
    'honey graham square': 'cereal, ready-to-eat, honey graham square',
    'honey nut scooter': 'cereal, ready-to-eat, honey nut oat ring',
    'honey nut wheat': 'cereal, ready-to-eat, shredded wheat with honey and nut',
    'honeycomb sphere with chocolate couverture malteser': 'honeycomb sphere with chocolate couverture',
    'honeydew': 'honeydew melon',
    'horlick, made up with semi-skimmed milk, dried': 'malted milk powder, made up with semi-skimmed milk',
    'horlick, made up with skimmed milk, dried': 'malted milk powder, made up with skimmed milk',
    'horlick, made up with whole milk, dried': 'malted milk powder, made up with whole milk',
    'horse pursland, leaf': 'horse purslane, leaf',
    'horse, filet': 'horse, fillet',
    'horse-radish': 'horseradish',
    'hot dogs beef': 'hot dog, beef',
    'hot pocket, meatballs & mozzarella stuffed sandwich': 'stuffed sandwich with meatball and mozzarella, frozen',
    "hot pockets ham 'n cheese stuffed sandwich": 'stuffed sandwich with ham and cheese, frozen',
    'hotdog': 'hot dog',
    'houmous': 'hummus',
    'hourse eye bean, black testa type': 'horse eye bean, black testa',
    'hourse eye bean, dark brown testa type': 'horse eye bean, dark brown testa',
    'hourse eye bean, light brown testa type': 'horse eye bean, light brown testa',
    'house foods firm tofu': 'tofu, firm',
    'house foods soft tofu': 'tofu, soft',
    'human milk, colostrum': 'breastmilk, colostrum',
    'hummus mashed chickpea': 'hummus',
    'hummus natural': 'hummus',
    'hummus, other - az': 'hummus',
    'hummus, sabra - ga': 'hummus, classic',
    'hummus, sabra - ky': 'hummus, classic',
    'hummus, tribe - ga': 'hummus, classic',
    'hummus, tribe - ky': 'hummus, classic',
    'hungry man, salisbury steak with gravy': 'salisbury steak with gravy',
    'härkis bolognese': 'broad bean mince bolognese',
    'härkis lasagne': 'broad bean mince lasagne',
    'härkis, broad bean product': 'broad bean product',
    'ice lolly e.g saftis': 'ice lolly, fruit juice',
    'ice lolly festini': 'ice lolly',
    'ice tea with sugar 4-<5 g cho': 'iced tea with sugar, 4-<5 g carbohydrate',
    'icing, chocolate': 'frosting, chocolate',
    'icing, chocolate, dried': 'chocolate frosting, dry mix',
    'indian penny worth, leaf': 'indian pennywort, leaf',
    'indian, papad': 'papad',
    'infant formula combiotik': 'infant formula, organic',
    'infant formula little steps': 'infant formula',
    'infant formula nutrasense': 'infant formula',
    'infant formula, enfamil': 'infant formula, ready-to-feed',
    'infant formula, enfamil enspire powder': 'infant formula, dried',
    'infant formula, enfamil for supplementing': 'infant formula for supplementing, ready-to-feed',
    'infant formula, enfamil for supplementing, dried': 'infant formula for supplementing, dried',
    'infant formula, enfamil gentlease': 'infant formula, reduced-lactose, dried',
    'infant formula, enfamil lipil': 'infant formula with iron, ready-to-feed',
    'infant formula, enfamil premature 30 calorie': 'preterm infant formula, 30 kcal per fl oz',
    'infant formula, enfamil premature high protein 24 calorie': 'preterm infant formula, high protein, 24 kcal per fl oz',
    'infant formula, enfamil reguline': 'infant formula with prebiotics, ready-to-feed',
    'infant formula, enfamil reguline powder': 'infant formula with prebiotics, powder',
    'infant formula, enfamil, dried': 'infant formula, dried',
    'infant formula, enfamil, paste': 'infant formula, liquid concentrate',
    'infant formula, gentlease': 'infant formula, reduced-lactose, ready-to-feed',
    'infant formula, good start': 'infant formula, ready-to-feed',
    'infant formula, good start essentials soy': 'soy infant formula, dried',
    'infant formula, good start soy': 'soy infant formula, ready-to-feed',
    'infant formula, good start supreme': 'infant formula with iron, ready-to-feed',
    'infant formula, good start, dried': 'infant formula, dried',
    'infant formula, next step': 'follow-on soy formula, ready-to-feed',
    'infant formula, next step prosobee, dried': 'follow-on soy formula, dried',
    'infant formula, next step, dried': 'follow-on soy formula, dried',
    'infant formula, nutrition': 'infant formula with iron, ready-to-feed',
    'infant formula, nutrition, dried': 'hydrolysed infant formula, dried',
    'infant formula, pbm product': 'infant formula, store brand, ready-to-feed',
    'infant formula, pbm product, dried': 'infant formula, store brand, dried',
    'infant formula, pregestimil': 'hydrolysed infant formula, dried',
    'infant formula, pregestimil 20 calorie': 'hydrolysed infant formula, 20 kcal per fl oz',
    'infant formula, pregestimil 24 calorie': 'hydrolysed infant formula, 24 kcal per fl oz',
    'infant formula, prosobee': 'soy infant formula with iron, ready-to-feed',
    'inobuta, lean and fat, flesh': 'pig-wild boar crossbred (inobuta), flesh',
    'instant noodle soup with spice mix diff. flavour': 'instant noodle soup with spice mix, various flavours',
    'interstate brands corp, wonder hamburger roll': 'hamburger roll',
    'irish potato': 'potato',
    'ironweed, just leaf': 'ironweed, mature leaf',
    'italian meringue with biscuit base milk chocolate couverture couconut kokostopp': 'italian meringue on biscuit base with milk chocolate and shredded coconut',
    'italian meringue with biscuit base milk chocolate couverture mums-mum': 'italian meringue with biscuit base and milk chocolate couverture',
    'italian pork sausage': 'italian pork sausage, cooked',
    'italian pork sausage, great value': 'italian pork sausage',
    'italian pork sausage, johnsonville hot': 'italian pork sausage, hot',
    'italian pork sausage, johnsonville mild': 'italian pork sausage, mild',
    'italian pork sausage, kroger': 'italian pork sausage',
    'italian pork sausage, store/other': 'italian pork sausage',
    'italian sausage salsiccia': 'italian sausage salsiccia, raw',
    'jack bean, for 1 h, dried': 'jack bean, boiled for 1 h, dried',
    'jack bean, for 2 h, dried': 'jack bean, boiled for 2 h, dried',
    'jack bean, for 3 h, dried': 'jack bean, boiled for 3 h, dried',
    'jack bean, for 4 h, dried': 'jack bean, boiled for 4 h, dried',
    'jack bean, until brownish colour': 'jack bean, roasted until brownish colour',
    'jacque, flesh': 'jackfruit, flesh',
    'jagermeister herb liquor': 'herb liqueur',
    'jakfruit': 'jackfruit, seed',
    "jansson's casserole, jansson's temptation": "jansson's temptation",
    'japanese radish, daikon, paste': 'japanese radish, daikon, pickled in rice bran',
    'japanese worcester sauce, sweet thick type for okonomiyaki': 'japanese worcester sauce, sweet thick type',
    'java-plum': 'jambolan',
    'jelly sweet': 'jelly candy',
    'jerusalem-artichoke, tuber': 'jerusalem artichoke, tuber',
    'jimmy dean, sausage': 'sausage, egg and cheese breakfast biscuit',
    'juice drink capri-sun multivitamin': 'juice drink, multivitamin',
    'juice drink dubbelfrisss': 'juice drink with sugar',
    'juice drink ocean spray cranberry': 'juice drink, cranberry',
    'juice drink tintelfruit light': 'juice drink light',
    'juice drink tintelfruit with sugar': 'juice drink with sugar',
    'juice drink with sugar and sweetener 4-6 g cho': 'juice drink with sugar and sweetener, 4-6 g carbohydrate',
    'juice fruit healthy people cranberry': 'cranberry juice',
    'juice tomato appelsientje zontomaat': 'tomato juice',
    'juice tomato/vegetable appelsientje tomatientje': 'juice tomato/vegetable',
    'jujube, wihout': 'jujube, without seeds',
    'kaki': 'persimmon',
    'kale product': 'kale, frozen',
    'kidney pork': 'pork, kidney',
    'kidney, salt': 'kidney, boiled, salted',
    'kielbasa, polish': 'kielbasa, turkey and beef, smoked',
    'kinder chocolate egg': 'milk chocolate egg with milk filling',
    'king vitaman': 'cereal, ready-to-eat, vitamin-fortified corn and oat',
    'kiwano': 'horned melon',
    'kiwi': 'kiwi fruit',
    'kiwifruit': 'kiwi fruit',
    'knol-khol': 'kohlrabi',
    'kohl rabi': 'kohlrabi',
    'kovai, small': 'ivy gourd',
    'krabask bitter': 'bitter liqueur',
    'lablab bean': 'hyacinth bean',
    'lablab purpureus, dried': 'hyacinth bean, dried',
    'lager, alcohol-free': 'beer, alcohol-free',
    'lager, strong': 'lager, extra strong',
    'lamb flank of lamb': 'lamb, flank',
    "lamb's, leaf": "lamb's quarters, leaf",
    'lamb, filet': 'lamb, fillet',
    'lamb, pancrea': 'lamb, pancreas',
    'langoustine': 'norway lobster',
    'langoustine, wild': 'norway lobster, wild',
    'lard or pork fat': 'lard',
    'large round sandwich with lettuce, tuna': 'large round sandwich with lettuce, tuna, anchovy and black olive',
    'large white beans without brine': 'bean, large white, canned, without brine',
    'lasagna classico': 'lasagna with meat and cheese, restaurant',
    'lasagna heated product': 'lasagna, frozen, heated',
    'lasagna veg. with soy protein': 'vegetable lasagna with soy protein',
    'lasagna w. zucchini eggplant beef mince meat': 'lasagna with zucchini, eggplant and minced beef',
    'lasagne bolognese, ready meal': 'lasagna bolognese, ready meal',
    'lemon juice, or bottled': 'lemon juice, canned or bottled',
    'lemon speciality to be diluted, for beverages or culinary use': 'lemon speciality to be diluted, for beverages or culinary use, no added sugars',
    'lemonada, limeade': 'limeade',
    'lemonade with syrupwith sugar': 'lemonade with syrup, with sugar',
    'lettuce av': 'lettuce',
    'lichee': 'lychee',
    'lime juice, or bottled': 'lime juice, canned or bottled',
    'limea bean': 'lima bean',
    'limea bean, water-soaked': 'lima bean, water-soaked',
    'lingonberry, cowberry': 'lingonberry',
    'linseed oil, flaxseed oil': 'flaxseed oil',
    'linseed or flaxseed oil': 'flaxseed oil',
    'lipton brisk, tea, black': 'tea, black, ready-to-drink, lemon',
    'liqueur purified or vodka vol. %': 'plain spirit or vodka, 40 vol%',
    'liqueur reduced sweetness vol. %': 'liqueur, reduced sweetness, 24 vol%',
    'liqueur sweet vol. %': 'liqueur, sweet, 24 vol%',
    'liqueur vol. % 38 kaptenlöjtnant': 'liqueur, 38 vol%',
    'liquid margarine non-diary': 'liquid margarine, non-dairy',
    'liquid margarine non-diary, enriched': 'liquid margarine, non-dairy, enriched',
    'liquid milk': 'milk, raw, whole',
    'liquid milk, holstein': 'milk, raw, whole',
    'liquid milk, jersey': 'milk, raw, whole',
    'liquor flavoured vol. % 40 brännvin': 'flavoured spirit, 40 vol%',
    'liquorice dutch type double, salted': 'liquorice, dutch double-salted',
    'liquorice stophoest': 'liquorice cough pastille',
    'litchi': 'lychee',
    'lithchi, flesh': 'lychee, flesh',
    'little hogweed, just leaf': 'little hogweed, mature leaf',
    "liver calf's": 'calf, liver',
    'liver casserole, kartanon': 'liver casserole, ready meal',
    'liver haddock tinned': 'haddock liver, tinned',
    'liver sausage, liver pate': 'liver sausage',
    'liver sausage, liverwurst': 'liver sausage',
    'liver, salt': 'liver, boiled, salted',
    "lmab's, leaves and shoot": "lamb's quarters, leaves and shoot",
    'lolly, hard variety': 'lolly, hard',
    'long bean': 'yardlong bean',
    'loquat, nécyar de cristal': 'loquat, néctar de cristal',
    'lotus tuber': 'lotus root',
    'low-fat margarine 40% fat latt & lagom': 'margarine, 40% fat, low-fat',
    'lucozade': 'glucose energy drink, carbonated',
    'lump fish': 'lumpfish',
    'lump roe, semi-preserved': 'lumpfish roe, semi-preserved',
    'm&m pretzel chocolate candy': 'pretzel coated with chocolate',
    "macadamia nut, ann's house of nut, dried": 'macadamia nut, dried',
    'macadamia, salted': 'macadamia nut, salted',
    "macaroni & cheese, from kid's menu": 'macaroni cheese, restaurant',
    "macaroni & cheese, from kids' menu": 'macaroni cheese, restaurant',
    'macaroni and cheese, mix': 'macaroni and cheese, prepared from dry mix',
    'macaroni and härkis casserole, low-fat': 'macaroni and broad bean mince casserole, low-fat',
    'macaroni from the alp': 'macaroni from the alps (älplermagronen)',
    "macaroni n' cheese plate, from kid's menu": 'macaroni cheese plate, restaurant',
    'macaroni or noodles with cheese, made from reduced fat packaged mix': 'macaroni or noodle with cheese dry mix, unprepared, low-fat',
    'macaroni or noodles with cheese, microwaveable': 'macaroni or noodle with cheese, microwaveable, unprepared',
    'macaroni, dark, flour': 'macaroni, dark, from whole wheat and dark wheat flour',
    'macaroni, pasta': 'macaroni, egg pasta',
    'macaroni, pasta, unsalted': 'macaroni, egg pasta, boiled',
    'macaroni, salt': 'macaroni, boiled, salted',
    'maccheronccini pasta, durum wheat semolina': 'maccheroncini pasta, durum wheat semolina',
    'maccheronccini pasta, durum wheat semolina and broad bean flour': 'maccheroncini pasta, durum wheat semolina and broad bean flour',
    'mackerel, filet': 'mackerel, fillet',
    'madeira wine vol. %': 'madeira wine, 18 vol%',
    'maize oil': 'corn oil',
    'malted milk drink, made with whole milk': 'malted milk powder, made up with whole milk',
    'malted milk powder, made up with whole milk': 'malted milk drink, made with whole milk',
    'maltesers and similar product': 'honeycomb spheres with chocolate coating',
    'mandarin orange': 'mandarin',
    'manderin': 'mandarin',
    'mange-tout': 'snow pea',
    'mangetout pea': 'snow pea',
    'mango bassignac, flesh': 'mango, flesh',
    'mango josé, flesh': 'mango, flesh',
    'mango moussache, flesh': 'mango, flesh',
    'mangue julie, flesh': 'mango julie, flesh',
    'manioc': 'cassava',
    'maple brown sugar life cereal': 'cereal, ready-to-eat, maple and brown sugar oat square',
    'margarine 40% fat <17 g sat fa, salted': 'margarine, 40% fat, <17 g saturated fatty acid, salted, low-fat',
    'margarine 40% fat <17 g sat fa, salted, low-fat': 'margarine, 40% fat, <17 g saturated fatty acid, salted, low-fat',
    'margarine 40% fat latt & lagom, enriched': 'margarine, 40% fat, enriched, low-fat',
    'margarine 40% fat latt & lagom, low-fat, enriched': 'margarine, 40% fat, enriched, low-fat',
    'margarine 40%, rainbow': 'margarine, 40% fat',
    'margarine 60%, alentaja': 'margarine, 60% fat, with plant stanols',
    'margarine 60%, kultarypsi': 'margarine, 60% fat, rapeseed oil',
    'margarine 60%, rainbow': 'margarine, 60% fat, lactose-free',
    'margarine 80% fat >24 g sat fatty acids for nevo recipe, salted': 'margarine 80% fat >24 g sat fa, salted',
    'margarine 80%, kulta': 'margarine, 80% fat',
    'margarine 80%, rypsi': 'margarine, 80% fat, rapeseed oil, fluid',
    'margarine 80%, sunnuntai': 'margarine, 80% fat, rapeseed oil, fluid',
    'margarine for baking 80%, rainbow': 'margarine for baking, 80% fat',
    'margarine for cooking and baking non-diary': 'margarine for cooking and baking, non-dairy',
    'margarine for cooking and baking non-diary, enriched': 'margarine for cooking and baking, non-dairy, enriched',
    'margarine goed begin, low-fat': 'margarine, low-fat',
    'margarine liq 80% fat <17 g sat fa, salted': 'margarine, liquid, 80% fat, <17 g saturated fatty acid, salted',
    'margarine non-diary': 'margarine spread, non-dairy',
    'margarine non-diary, enriched': 'margarine spread, non-dairy, enriched',
    'margarine or blended for bread': 'margarine or blended spread for bread',
    'margarine prod ah chol reducing, low-fat': 'margarine, cholesterol-lowering, low-fat',
    'margarine prod goede start, low-fat': 'margarine, low-fat',
    'margarine prod omega3 plus, low-fat': 'margarine with omega-3, low-fat',
    'margarine product 60% fat >17 g sat fa, salted': 'margarine 60% fat >17 g sat fa, salted',
    'margarine product balan': 'margarine product, low-fat',
    'margarine product bewust light, low-fat': 'margarine, light, low-fat',
    'margarine product goed begin, low-fat': 'margarine, low-fat',
    'margarine product romig': 'margarine product, creamy, 60% fat',
    'margarine used in cakes raisa': 'margarine used in cakes',
    'margarine, light without saturated fat': 'margarine, light',
    'margarine, margarine vegetable oil': 'margarine-like vegetable oil spread, 67-70% fat',
    'margarine, margarine-type vegetable oil': 'margarine-type vegetable oil spread, 70% fat',
    'margarine, omega plus': 'margarine-like spread with plant sterols and fish oil',
    'margarine, smart squeeze': 'margarine',
    'marli smoothie': 'berry smoothie with added vitamins',
    'marmalade diff. flavour': 'marmalade, various flavours',
    "maroilles laitier cheese, from cow's milk": 'maroilles laitier cheese',
    'marshmallow matey': 'cereal, ready-to-eat, oat with marshmallow',
    "martha white food, martha white's buttermilk biscuit mix, dried": 'buttermilk biscuit mix, dried',
    "martha white food, martha white's chewy fudge brownie mix, dried": 'fudge brownie mix, dried',
    "mary's gone cracker, original cracker": 'cracker, gluten-free',
    'marzipanproduct, chocolatecover': 'marzipan product with chocolate coating',
    'mashed potato powder, eldorado': 'mashed potato powder',
    'mashed potato powder, felix': 'mashed potato powder',
    'mashed potato powder, rainbow': 'mashed potato powder',
    'mcchicken sandwich': 'breaded chicken patty sandwich',
    'mcchicken sandwich, without mayonnaise': 'breaded chicken patty sandwich, without mayonnaise',
    'mckee baking, little debbie nutty bar': 'wafer bar with peanut butter, chocolate covered',
    'meat av excl liver': 'meat, excluding liver',
    'meat pastie, portti, fried': 'meat pasty, deep fried',
    'meatball, enebacken': 'meatball, industrial',
    'meatball, hk': 'meatball, industrial',
    'meatball, kartanon': 'meatball, industrial',
    'meatball, mummon lihapulla': 'meatball, industrial',
    'melon cantaloupe': 'cantaloupe',
    'melon water': 'watermelon',
    'mentha, leaf': 'field mint, leaf',
    'mento, cheewy dragee': 'chewy mint dragee',
    'mesclun or salad': 'mesclun',
    'mexican, burrito with bean': 'burrito with bean',
    'mexican, burrito with beans and beef': 'burrito, with beans and beef',
    'mexican, burrito with beans and cheese': 'burrito, with beans and cheese',
    'mexican, burrito with beans and chili pepper': 'burrito with beans and chili pepper',
    'mexican, burrito with beef': 'burrito with beef',
    'mexican, burrito with beef and chili pepper': 'burrito with beef and chili pepper',
    'mexican, chimichanga with beef': 'chimichanga with beef',
    'mexican, chimichanga with beef and cheese': 'chimichanga with beef and cheese',
    'mexican, chimichanga with beef and red chili': 'chimichanga with beef and red chili',
    'mexican, enchilada with cheese': 'cheese enchilada',
    'mexican, enchilada with cheese and beef': 'enchilada with cheese and beef',
    'mexican, enchirito with cheese': 'enchirito with cheese, beef and bean',
    'mexican, frijoles with cheese': 'frijoles with cheese',
    'mexican, nachos with cheese': 'nachos, with cheese',
    'mexican, taco salad': 'taco salad',
    'mexican, taco salad with chili con carne': 'taco salad with chili con carne',
    'mexican, taco with beef': 'taco with beef, cheese and lettuce',
    'mexican, taco with chicken': 'taco with chicken, lettuce and cheese',
    'mexican, tostada with bean': 'tostada with beans, beef and cheese',
    'mexican, tostada with beans and cheese': 'tostada with beans and cheese',
    'mexican, tostada with beef and cheese': 'tostada with beef and cheese',
    'mexican, tostada with guacamole': 'tostada with guacamole',
    'mifu, lactose-free milk protein product': 'milk protein product, lactose-free',
    'milk based drink yakult balance': 'fermented milk drink with probiotics, reduced sugar',
    'milk based drink yakult original': 'fermented milk drink with probiotics',
    'milk chocolate candy with peppermint flavour, marianne': 'milk chocolate candy with peppermint filling',
    'milk chocolate covered biscuits with filling, amanie': 'milk chocolate covered biscuit with filling',
    'milk chocolate with chrunchy almond caramel centre': 'milk chocolate with crunchy almond caramel centre',
    'milk chocolate-flavoured with s-sk milk and sweetened cocoa powder': 'milk chocolate-flavoured with semi-skimmed milk and sweetened cocoa powder',
    'milk condensed with sugar': 'milk, condensed, sweetened',
    'milk goats- full-fat': 'goat milk, whole',
    'milk of mare, iitalian saddle': 'mare milk, italian saddle breed',
    'milk of toggenburg': 'goat milk, toggenburg breed',
    'milk pudding with milk powder egg kalvdan': 'milk pudding with milk powder and egg',
    'milk, 1 ug added vitamin d': 'milk, 1.5% fat, with added vitamin D',
    'milk, i ug added vitamin d': 'milk, with added vitamin D, fat-free',
    'milk, i ug added vitamin d, fat-free': 'milk, with added vitamin D, fat-free',
    'milk, uht, low-fat': 'milk, low-fat',
    'milkbased drink yakult plus': 'fermented milk drink with probiotics',
    'milkshake, hesburger, vanilla': 'milkshake, from fast food',
    'milkshake, power cow, chocolate': 'milkshake, chocolate flavour',
    'milky way and own brand equivalent': 'chocolate bar with nougat',
    'minced beef': 'beef mince',
    'minced liver steak, potato': 'minced liver steak with potato and onion',
    'minced meat beefburger, horeca': 'minced meat beefburger, food service',
    'minced meat hamburger patty, grillitassu, fried': 'minced meat hamburger patty, fried without fat',
    'minced meat hamburger patty, hk': 'minced meat hamburger patty',
    'minced meat patties with pickled beetroot cucumber product biff à la lindström': 'minced meat patties with pickled beetroot and cucumber, biff à la lindström, fried',
    'minced meat vegetarian based on soy with iron and vitamin b12 de vegetarische slager, enriched': 'minced meat vegetarian based on soy with iron and vitamin b12, enriched',
    'minced meat vegetarian based on soya with iron and vitamin b12 de vegetarische slager': 'minced meat vegetarian based on soy with iron and vitamin b12',
    'mineral water, novelle plus': 'mineral water with added calcium',
    'mineral water, novelle plus, enriched': 'mineral water with added calcium, enriched',
    'mineral water, novelle plus, sweetened': 'mineral water, sweetened with sugar',
    'mineral water, plus, enriched': 'mineral water with added calcium, enriched',
    'mineral water, plus, sweetened': 'mineral water, sweetened with sugar',
    'mini pretzel rods salta pinnar': 'mini pretzel rod',
    'minibread, retail': 'mini bread, toasted, retail',
    'mission food, mission flour tortilla': 'flour tortilla, soft taco size',
    'mix seasoning mexican': 'mexican seasoning mix',
    'mixed vegetable, salt': 'vegetable, mixed, salted',
    'mixed, salted': 'mixed nut, salted',
    'mollusc, blue mussel': 'mussel, blue',
    "mom's best, honey nut toasty o's": 'cereal, ready-to-eat, honey nut oat ring',
    "mom's best, sweetened wheat-ful": 'cereal, ready-to-eat, sweetened shredded wheat',
    'monosodium glutamate ve-tsin': 'monosodium glutamate',
    'monster energy drink, low carb': 'energy drink, low carbohydrate',
    'moose product': 'moose, sautéed, frozen product',
    'mori-nu, tofu': 'tofu, silken, extra firm',
    'mothbean': 'moth bean',
    'mothbean, salted': 'moth bean, salted',
    'mothbean, unsalted': 'moth bean',
    "mother's cinnamon oat crunch": 'cereal, ready-to-eat, cinnamon oat crunch',
    "mother's cocoa bumper": 'cereal, ready-to-eat, cocoa puff',
    "mother's graham bumper": 'cereal, ready-to-eat, graham flavoured puff',
    "mother's oat bran": 'oat bran cereal, prepared with water and salt',
    "mother's oat bran cereal": 'cereal, ready-to-eat, toasted oat bran',
    "mother's oat bran, dried": 'oat bran cereal, dried',
    "mother's oat bran, unsalted": 'oat bran cereal, prepared with water',
    "mother's peanut butter bumpers cereal": 'cereal, ready-to-eat, peanut butter flavoured puff',
    'mott, apple juice light': 'apple juice, light, fortified with vitamin C',
    'mott, apple juice light with vitamin c, enriched': 'apple juice, light, fortified with vitamin C',
    'moui béinré, rice porridge with water': 'moui béinré, rice porridge with water, milk and sugar',
    'mozzarella': 'cheese, mozzarella',
    'mucuna monosperma dc ex wight, dried': 'mucuna monosperma, dried',
    'mucuna prurien, dried': 'velvet bean, dried',
    'muesli, added fruit & nut, dried': 'muesli, toasted, with added dried fruit and nut',
    'muesli, kilomysli, dried': 'muesli with raisin, corn, oat and dried fruit',
    'muesli, untoasted or natural style, dried': 'muesli, untoasted, with added dried fruit',
    'muffin, english style': 'muffin, english',
    'muffin, english style, dried': 'english muffin, from white flour, with dried fruit, toasted',
    'muffin, english style, flour': 'english muffin, from white flour',
    'muksu soy drink': 'soy drink for children',
    'mulled wine, non-alcololic': 'mulled wine, non-alcoholic',
    'multigrain bread av white and brown': 'multigrain bread with seed, average of white and brown',
    'multigrain bread av white and brown with with flaxseed': 'multigrain bread, average of white and brown, with seed and extra flaxseed',
    'multigrain bread brown with vikorn vitaminebrood, enriched': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread brown with vikorn volvezel, enriched': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread brown with with iron and vitamins vikorn vitaminebrood': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread brown with with iron and vitamins vikorn vitaminebrood, enriched': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread brown with with iron and vitamins vikorn volvezel': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread brown with with iron and vitamins vikorn volvezel, enriched': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread, buttermilk bread, flour': 'multigrain buttermilk bread, from wheat and dark wheat flour',
    'multigrain bread, carrot roll, flour': 'multigrain bread, carrot roll',
    'multigrain bread, faltbread': 'multigrain flatbread with pearl barley and sour milk',
    'multigrain bread, sweet-and-sour bread, flour': 'multigrain sweet-and-sour bread, from wheat and rye flour',
    'multigrain bread, water, flour': 'multigrain bread, from graham, rye and wheat flour',
    'mungo bean': 'black gram',
    'mushroom, chanterelle': 'chanterelle',
    'mushroom, hon-shimiji, white': 'mushroom, honshimeji, white',
    'mushroom, lions mane': "mushroom, lion's mane",
    'musk melon': 'muskmelon',
    'mustard sauce with oil vinegar dill gravlaxså': 'mustard sauce with oil, vinegar and dill, chilled',
    'mustard sauce with oil vinegar dill product gravlaxså': 'mustard sauce with oil, vinegar and dill, chilled',
    'nachos supreme': 'nachos with beef, bean and sour cream',
    'naturally sparkling': 'sparkling mineral water',
    "nature's path": 'cereal, ready-to-eat, flaxseed flake',
    "nature's path, pumpkin granola": 'granola, ready-to-eat, with pumpkin seed and flaxseed',
    'navy bean': 'bean, navy',
    'nestea, tea, black': 'tea, black, ready-to-drink, lemon',
    'non-heading chinese cabbage, osaka-shirona, salted, leaf': 'non-heading chinese cabbage, salted',
    'nos energy drink': 'energy drink, with sugar',
    'nos energy drink, enriched': 'energy drink, with sugar, enriched',
    'novelty, klondike': 'novelty, ice cream type',
    'novelty, no sugar added creamsicle pop': 'ice pop, orange cream, no sugar added',
    'nut and fruit mix, nut, dried, salted': 'nut and fruit mix, dried, salted',
    'nut and roast, mixed nuts and sunflower seed': 'nut and seed roast, mixed nuts and sunflower seeds',
    'nut and roast, mixed nuts and sunflower seeds with egg': 'nut and seed roast, mixed nuts and sunflower seeds, with egg',
    'nut mix, salted': 'mixed nut, salted',
    'nut mix, unsalted': 'mixed nut',
    'nut, brazil': 'brazil nut',
    'nut, coco': 'coconut',
    'nut, formulated, salted': 'wheat-based nut substitute, salted',
    'nut, formulated, unsalted': 'wheat-based nut substitute',
    'nut, pea, salted': 'peanut, oil-roasted, salted',
    'nut, pine': 'pine nut',
    'nut, simulated product, salted': 'wheat-based nut substitute, salted',
    'nut, simulated product, unsalted': 'wheat-based nut substitute',
    'nutritional supplement for people with diabete, liquid': 'nutritional supplement for people with diabetes, liquid',
    'nuts macadamia, unsalted': 'macadamia nut',
    'nuts mixed, unsalted': 'mixed nut',
    'oat and vanillla biscuit, jyväshyvä': 'oat and vanilla biscuit',
    'oat and vanillla biscuit, lu': 'oat and vanilla biscuit',
    'oat biscuit wholegrain hobnob': 'oat biscuit, wholegrain',
    'oat biscuit, jyväshyvä luonnonhyvät': 'oat biscuit, rapeseed oil-based',
    'oat biscuit, luonnonhyvät': 'oat biscuit, rapeseed oil-based',
    'oat blenders with honey': 'cereal, ready-to-eat, oat cluster with honey',
    'oat blenders with honey & almond': 'cereal, ready-to-eat, oat cluster with honey and almond',
    'oat bran flake, health valley': 'cereal, ready-to-eat, oat bran flake',
    'oat bran, danish': 'oat bran',
    'oat bran, extruded at 115ºc and 20%h20': 'oat bran, extruded at 115ºc and 20% h2o',
    'oat cinnamon life': 'cereal, ready-to-eat, oat and corn square with cinnamon',
    'oat crispbread, oululainen': 'oat crispbread',
    'oat drink without milk, ikaffe': 'oat drink without milk, for coffee',
    'oat drink without milk, nalle chocolate': 'oat drink without milk, chocolate, for children',
    'oat drink without milk, nalle vanilla': 'oat drink without milk, vanilla, for children',
    'oat drink without milk, nordic': 'oat drink, with calcium and vitamin D',
    'oat drink without milk, oatly': 'oat drink without milk, with calcium and vitamins',
    'oat drink without milk, oatly chocolate': 'oat drink without milk, chocolate',
    'oat drink without milk, oatly ikaffe': 'oat drink without milk, for coffee',
    'oat drink without milk, rainbow': 'oat drink without milk, with calcium and vitamins',
    'oat groats whole bolied, salted': 'oat groat, whole, boiled, salted',
    'oat ice cream non-diary flavoured': 'oat ice cream, non-dairy, flavoured',
    'oat ice cream with vanilla non-diary': 'oat ice cream with vanilla, non-dairy',
    'oat life': 'cereal, ready-to-eat, oat and corn square',
    'oat rice bolied, salted': 'oat rice, boiled, salted',
    'oat spread, oatly påmackan': 'oat spread, 20% fat',
    'oat, meal': 'oat meal, ground',
    'oat, wheat and honey': 'granola with oat, wheat and honey',
    'oatly havregurt, oat product': 'oat-based yoghurt alternative, flavoured',
    'oatmeal cookies with raisin': 'oatmeal cookies with raisin, packaged',
    'oatmeal cookies with raisin, pepperidge farm soft': 'oatmeal cookie with raisin, soft baked',
    'oatmeal cookies with raisin, store/other brand': 'oatmeal cookies with raisin',
    'oatmeal square': 'cereal, ready-to-eat, oat square',
    'oatmeal square, cinnamon': 'cereal, ready-to-eat, oat square with cinnamon',
    'oatmeal square, golden maple': 'cereal, ready-to-eat, oat square with maple',
    'oatmeal, real medley, dried': 'oatmeal with apple and walnut, dried',
    'oats toast flour': 'oat toast bread, from oat flour and oat flake',
    'ocean spray, cran cherry': 'cranberry-cherry juice drink',
    'ocean spray, cran grape': 'juice drink, cranberry-grape',
    'ocean spray, cran lemonade': 'cranberry lemonade drink',
    'ocean spray, cran pomegranate': 'cranberry-pomegranate juice drink',
    'ocean spray, cran raspberry juice drink': 'cranberry-raspberry juice drink',
    'ocean spray, cran-energy, juice': 'cranberry energy juice drink',
    'ocean spray, cranberry-apple juice drink': 'cranberry-apple juice drink, bottled',
    'ocean spray, diet cran cherry': 'cranberry-cherry juice drink, diet',
    'ocean spray, diet cranberry juice': 'cranberry juice drink, diet',
    'ocean spray, light cranberry': 'cranberry juice drink, light',
    'ocean spray, light cranberry and raspberry flavored juice': 'cranberry-raspberry juice drink, light',
    'ocean spray, ruby red cranberry': 'cranberry-grapefruit juice drink',
    'ocean spray, white cranberry peach': 'white cranberry-peach juice drink',
    'ocean spray, white cranberry strawberry flavored juice drink': 'white cranberry-strawberry juice drink',
    'oil blend': 'vegetable oil, blend',
    'oil for frying and deep fat frying, kultasula paistoöljy': 'oil for frying and deep-fat frying',
    'oil palmkernel': 'palm kernel oil',
    'oil soy': 'soybean oil',
    'oil wok': 'wok oil',
    "oil, from cow's milk butter": 'ghee',
    'olive oil virgine': 'olive oil, extra virgin',
    'olive, green manzanilla - pa': 'olive, manzanilla, green',
    'olive, small to large': 'olive, ripe, canned',
    'onion yellow product': 'onion, yellow, frozen',
    'onion, spring or scallion': 'onion, spring',
    'or bottled juice': 'juice, canned or bottled, sweetened',
    'orange mango with vitamin, c, enriched': 'orange mango juice drink, fortified',
    'orange mango, enriched': 'orange mango juice drink, fortified',
    'orange mango, with vitamin': 'orange mango juice drink, fortified',
    'orange-grapefruit juice, or bottled': 'orange-grapefruit juice, canned or bottled, unsweetened',
    'oriental sauce with vegetables uncle bens orientalisk sås': 'oriental sauce with vegetable, canned',
    'original chicken sandwich': 'breaded chicken sandwich, takeaway',
    'original round cheese pizza, regular crust': 'cheese pizza, regular crust',
    'original round meat and vegetable pizza, regular crust': 'meat and vegetable pizza, regular crust',
    'original round pepperoni pizza, regular crust': 'pepperoni pizza, regular crust',
    'original taco with beef, cheese and lettuce': 'hard shell taco with beef, cheese and lettuce',
    'other milk, goat milk': 'goat milk',
    'other milk, human milk': 'human milk',
    'ovaltine beverage with reduced fat milk, without sugar': 'malted chocolate beverage with reduced fat milk, without sugar',
    'ovaltine beverage with skimmed milk, without sugar': 'malted chocolate beverage with skimmed milk, without sugar',
    'ovaltine beverage with whole milk, without sugar': 'malted chocolate beverage with whole milk, without sugar',
    'ovaltine powder': 'malted chocolate drink powder, fortified',
    'ovaltine powder, made up with whole milk': 'malted chocolate drink powder, made up with whole milk',
    'ovaltine, chocolate malt powder': 'chocolate malt drink powder, fortified',
    'ovaltine, malt powder': 'malt drink powder',
    'pak choy': 'bok choy',
    'palm fat': 'palm kernel fat',
    'palm oil with vitamin, 1100–2400 mcg/100g, enriched': 'palm kernel oil, with vitamin A, enriched',
    'palm oil, palm fat': 'palm oil',
    'palm oil, unfortified': 'palm kernel oil, unfortified',
    'palm oil, with vitamin': 'palm kernel oil, with vitamin A',
    'palm, tenera': 'palm kernel, tenera variety',
    'pam cooking spray oil': 'cooking spray oil',
    'pancake, aunt jemima': 'pancake, commercial',
    'paneer': 'cheese, paneer',
    'pangasius, filet': 'pangasius, fillet',
    'panto - beef, top loin steak, lean': 'beef, top loin steak, boneless, lean',
    'panto - beef, top round roast/steak, lean': 'beef, top round roast/steak, lean',
    'parmesan cheese, 100%': 'cheese, parmesan, grated',
    'parsnip, danish': 'parsnip',
    'partly skimmed milk': 'milk, low-fat',
    'passionfruit': 'passion fruit',
    'pasta salat with italian dressing': 'pasta salad with italian dressing',
    'pasta sauce diff. brand': 'pasta sauce, various brands',
    'pasta white av': 'pasta, white, boiled',
    'pasta-roe-shrimpsalad, avocado': 'pasta, roe and shrimp salad with avocado',
    'pastry with arrak oatmeal crushed bisquits arraksboll': 'pastry with arrack, oatmeal and crushed biscuit',
    'pastry with soft almond macaroon buttercream chocolate chokladbiskvi': 'pastry with almond macaroon, buttercream and chocolate',
    'pastry, greek': 'baklava',
    'pate, chicken liver': 'pâté, chicken liver',
    'pate, liver': 'pâté, liver',
    'pawpaw': 'papaya',
    'pea protein mince or product': 'pea protein mince',
    'pea sprout': 'pea, sprouted',
    'pea, chick': 'chickpea',
    'pea, mange-tout': 'snow pea',
    "peach pie, marie callender's": 'peach pie, frozen',
    'peanut butter, skippy': 'peanut butter, extra crunchy',
    'peanut oil with vitamin, 600–1000 mcg/100g, enriched': 'peanut oil, fortified with vitamin a',
    'peanut oil, with vitamin': 'peanut oil, fortified with vitamin a',
    'peanut, bean': 'peanut, immature, boiled',
    'peanut, chinese, dried': 'peanut, dried',
    'peanut, coarsely': 'peanut, coarsely ground',
    'peanut, combined variety, dried': 'peanut, dried',
    'peanut, for 30 sec': 'peanut, microwaved 30 seconds',
    'peanut, for 45 sec': 'peanut, microwaved 45 seconds',
    'peanut, for 60 sec': 'peanut, microwaved 60 seconds',
    'peanut, manipintar variety, dried': 'peanut, dried',
    'peanut, sinkarzie variety, dried': 'peanut, dried',
    'peanut, spanish type': 'peanut, spanish',
    'peanut, spanish, salted': 'peanut, spanish, oil-roasted, salted',
    'peanut, with oil, salted': 'peanut, oil-roasted, salted',
    'peanutsauce homemade with s-sk milk': 'peanut sauce, homemade with semi-skimmed milk, without added fat',
    'peanutsauce homemade with water': 'peanut sauce, homemade with water, without added fat',
    'pear cider': 'pear cider, alcohol-free',
    'peasoup with beef steak mince': 'pea soup with minced beef steak',
    'pecannuts oil, salted': 'pecan nut, oil-roasted, salted',
    'peiti beurre with chocolate, cookie': 'petit beurre biscuit with chocolate',
    'pencil yam, maloga bean': 'pencil yam, root',
    'pepperidge farm, goldfish': 'baked cheddar snack cracker',
    "pepperoni pizza, thin 'n crispy crust": 'pepperoni pizza, thin crispy crust',
    'pepperoni pizza, ultimate deep dish crust': 'pepperoni pizza, deep dish crust',
    'pepsico, gatorade': 'sports drink, ready-to-drink',
    'physalis': 'cape gooseberry',
    'pie, apple, flour': 'pie, apple',
    'pigeonpea': 'pigeon pea',
    'pigeonpea, salted': 'pigeon pea, salted',
    'pigeonpea, unsalted': 'pigeon pea',
    'pillow pak turkey pepperoni': 'turkey pepperoni, sliced',
    'pine nut, pine': 'pine nut',
    'pirogue filled with meat heated produkt': 'pirog filled with meat, frozen, heated',
    'pistachio': 'pistachio nut',
    "pistachio, ann's house of nut": 'pistachio nut',
    'pistachio, dried': 'pistachio nut, dried',
    'pistachio, paramount farm': 'pistachio nut, raw',
    'pizza capricciosa with smoked ham button mushrooms from restaurant': 'pizza capricciosa with smoked ham and button mushroom',
    "pizza, tony's breakfast pizza sausage": 'pizza, breakfast, with sausage',
    "pizza, tony's smartpizza wholegrain 4x6 cheese pizza 50/50 cheese": 'pizza, cheese, wholegrain',
    "pizza, tony's smartpizza wholegrain 4x6 pepperoni pizza 50/50 cheese": 'pizza, pepperoni, wholegrain crust',
    'plain cake': 'cake',
    'plain cake, chocolate cake': 'cake, chocolate',
    'plain cake, frangipani cake': 'frangipani cake',
    'plain cake, lemon cake': 'lemon cake',
    'plain chocolate candy': 'chocolate candy, plain',
    'plain, soymilk': 'soymilk',
    'plant-based alternative to fruit/vanilla yoghurt based on soy sweetened with calcium and vitamin, enriched': 'plant-based alternative to fruit/vanilla yoghurt based on soy, sweetened, fortified with calcium and vitamin',
    'plant-based alternative to fruit/vanilla yoghurt based on soya sweetened with calcium and vitamin': 'plant-based alternative to fruit/vanilla yoghurt based on soy, sweetened, fortified with calcium and vitamin',
    'plant-based alternative to yoghurt based on soy sweetened with calcium and vitamin, enriched': 'soy yoghurt alternative, sweetened, fortified with calcium and vitamin',
    'plant-based alternative to yoghurt based on soya sweetened with calcium and vitamin': 'soy yoghurt alternative, sweetened, fortified with calcium and vitamin',
    'plant-based bacon imitation soy- and wheat protein or product': 'plant-based bacon imitation soy- and wheat protein',
    'plant-based bits soy protein or product tzay': 'plant-based soy protein bits, chilled or frozen',
    'plant-based bits soy protein tzay': 'plant-based soy protein bits',
    'plant-based bits soy protein tzay, fried': 'plant-based soy protein bits, fried',
    'plant-based bits soy protein with thyme garlic product ®': 'plant-based soy protein bits with thyme and garlic',
    'plant-based kebab soy protein product ®': 'plant-based soy protein kebab',
    'plantain banana': 'plantain',
    'plantain cooking banana': 'plantain',
    'plus milk': 'milk, with added vitamin D, fat-free',
    'plus omega-3 dha, soymilk': 'soymilk, with added omega-3 dha',
    'pod, microloma saggitum': 'microloma sagittatum, pod',
    'pollock, saithe': 'saithe',
    'pollock, walley, wild': 'pollock, walleye, wild',
    "pomme d'eau or malaca apple, flesh": 'water apple/malay apple, flesh',
    "pont l'evêque cheese, from cow's milk": "pont l'evêque cheese",
    'poppy': 'poppy seed',
    'poppy seed, mccormick': 'poppy seed',
    'pork av': 'pork',
    'pork berkshire boars x, longissimus dorsi': 'pork, loin, berkshire cross',
    'pork chop, stekt': 'pork chop, breaded, fried',
    'pork ham aspic-gelly': 'pork ham in aspic jelly, canned',
    'pork picninc shoulder': 'pork picnic shoulder',
    'pork schnitzel not': 'pork schnitzel, not breaded',
    'pork tenderloin medaillon': 'pork tenderloin medallion',
    'pork, loin, longissimus dorsi': 'pork, loin, trimmed',
    'pork, sparerib': 'pork, spare rib',
    'porrdige buttermilk groat': 'porridge, buttermilk, groats',
    'porridge graham flour grahamsgröt': 'porridge, graham flour',
    'porridge milk with wheat flour white lammetjespap': 'porridge, milk with white wheat flour',
    'porridge non-dairy barley with raisins bessola': 'porridge, non-dairy barley with raisin',
    'porridge oat, unfortified': 'oat porridge, unfortified, made with semi-skimmed milk',
    'porridge pyjamapapje': 'porridge for children, ready-to-eat',
    'porridge, tô, flour': 'porridge, tô, from corn and cassava flour',
    'port wine white red vol. %': 'port wine, white or red, 20 vol%',
    'post bran flake': 'cereal, ready-to-eat, bran flake',
    'post great grains banana nut crunch': 'cereal, ready-to-eat, wholegrain with banana and nut',
    'post great grains cranberry almond crunch': 'cereal, ready-to-eat, wholegrain with cranberry and almond',
    'post honey bunches of oats with cinnamon bunch': 'cereal, ready-to-eat, oat and corn flake with cinnamon cluster',
    'post raisin bran cereal': 'cereal, ready-to-eat, raisin bran',
    'post selects blueberry morning': 'cereal, ready-to-eat, wholegrain flake with blueberry',
    'post selects maple pecan crunch': 'cereal, ready-to-eat, wholegrain flake with maple and pecan',
    'potat tuber, toyoshiro': 'potato tuber, toyoshiro',
    'potato biled commercial kitchen, salted': 'potato, boiled, commercial kitchen, salted',
    'potato pancake homemade raggmunk': 'potato pancake, homemade',
    'potato product': 'potato, diced, frozen',
    'potato, automn': 'potato, autumn',
    'potato, avarage': 'potato',
    'potato, floury av': 'potato, floury',
    'potato, french fry': 'potato, french, fried',
    'potato, french, fried': 'french fry, fried',
    'potato, french, fried, salted': 'french fry, salted',
    'potato, german butterball': 'potato, unpeeled, raw',
    'potato, waxy av': 'potato, waxy',
    'potatoe': 'potato',
    'potatoe, without addition of fat and salt': 'potato, without addition of fat and salt',
    'potatoe, without addition of salt': 'potato, without addition of salt',
    'powan, whitefish': 'whitefish, farmed',
    'powerade, lemon-lime flavored': 'sports drink, lemon-lime flavored',
    'powerade, zero': 'sports drink, calorie-free',
    'prawn': 'shrimp or prawn',
    'processed meat prod 20-30 g fat excl liver': 'processed meat product, 20-30 g fat, excluding liver',
    'processed meat prod >30 g fat excl liver': 'processed meat product, >30 g fat, excluding liver',
    'processed meat prod excl liver': 'processed meat product, excluding liver',
    'profeel protein drink, sugar-free': 'protein drink, sugar-free',
    'profeel voimamaitojuoma milk': 'milk protein drink, lactose-free, fat-free',
    'profeel voimamaitojuoma milk, fat-free': 'milk protein drink, lactose-free, fat-free',
    'prune de cythère, or pomme cythere or golden apple, mature': 'cythere apple',
    'pudding, low calorie, dried': 'pudding, low calorie',
    'puff pastry av': 'puff pastry, baked',
    'pumpkin and squash seed, salted': 'pumpkin and squash seed kernel, roasted, salted',
    'pumpkin and squash seed, whole': 'pumpkin and squash seed, whole, roasted',
    'pumpkin and squash seed, whole, salted': 'pumpkin and squash seed, whole, roasted, salted',
    'pumpkin and squash, dried': 'pumpkin and squash seed, dried',
    'pumpkin and squash, japanese squash': 'japanese squash',
    'pumpkin and squash, whole, unsalted': 'pumpkin and squash seed, whole, roasted',
    'pumpkin, dehulled': 'pumpkin seed, dehulled, fat-free',
    'pumpkin, dehulled, fat-free': 'pumpkin seed, dehulled, fat-free',
    'pumpkin, hulled, dried': 'pumpkin seed, hulled, dried',
    'pupusa, con frijoles': 'pupusa with bean',
    'pupusa, con queso': 'pupusa with cheese',
    'pupusa, del cerdo': 'pupusa with pork',
    'quark flavoured artificial sweetener': 'quark, flavoured, with artificial sweetener',
    'quark with fruit danoontje': 'quark with fruit, for children',
    'quark, quark yoghurt': 'quark yoghurt',
    'quarter pounder': 'hamburger, large single beef patty',
    'quarter pounder with cheese': 'cheeseburger, large single beef patty',
    'queen vicoria anana, flesh': 'pineapple, queen victoria, flesh',
    'quenelle, veal': 'veal quenelle',
    'quorn bolognese': 'mycoprotein bolognese',
    'raddiccio': 'radicchio',
    'raisin, golden seedless': 'raisin, golden',
    'rajmah, black': 'kidney bean, black',
    'rajmah, brown': 'kidney bean, brown',
    'rajmah, red': 'bean, kidney, red',
    'ranch snack wrap': 'grilled chicken wrap with ranch sauce',
    'ranch snack wrap, crispy': 'crispy chicken wrap with ranch sauce',
    'rape oil': 'rapeseed oil',
    'rapeseed oil holl': 'rapeseed oil, high oleic low linolenic',
    'raspberries and blueberries product': 'raspberry and blueberry, frozen',
    'raspberries with sugar product': 'raspberry, sweetened, frozen',
    'red beet': 'beetroot',
    'red beetroot': 'beetroot',
    'red currant': 'redcurrant',
    'red currant, red dutch': 'redcurrant',
    'red fish, wild': 'redfish, wild',
    'red gram': 'pigeon pea',
    'reddi wip whipped topping, fat-free': 'whipped topping, fat-free',
    'redfish, ocean perch': 'redfish',
    'redfish, rosefish': 'redfish',
    'reindeer product': 'reindeer, sautéed, frozen product',
    'reindeer, salt': 'reindeer, boiled, salted',
    'rennin, chocolate': 'rennet dessert, chocolate, prepared with 2% milk',
    'rice cakes puffed with fruit flavour organix': 'rice cake, puffed, with fruit flavour',
    'rice drink, nordic': 'rice drink with calcium and vitamin D',
    'rice drink, rice dream chocolate': 'rice drink, chocolate, with calcium',
    'rice drink, rice dream hazelnut and almond': 'rice drink with hazelnut and almond',
    'rice flour, enriched': 'rice flour, fortified with vitamin b1',
    'rice-a-roni, chicken flavor': 'rice and pasta mix, chicken flavour, dried',
    'rich chocolate, dried': 'hot chocolate beverage powder, rich',
    'ricy': 'rice',
    'ritz cracker': 'cracker, round buttery',
    'rockmelon': 'cantaloupe',
    'roe, ragout': 'roe deer, ragout',
    'roe, schnitzel': 'roe deer, schnitzel',
    'roll, hamburger': 'hamburger roll, wholegrain white, calcium-fortified',
    'roll, hamburger, enriched': 'hamburger roll, wholegrain white, calcium-fortified',
    'rolled oats wholegrain': 'oat, rolled',
    'root, hemidesmus indicus var. indicus': 'indian sarsaparilla, root',
    'rose hip, rose hip powder, dried': 'rose hip, dried',
    'rosehip soup rte instant powder with o sugar': 'rosehip soup from instant powder, without sugar, fortified',
    'rosehip soup rte instant powder with o sugar, enriched': 'rosehip soup from instant powder, without sugar, fortified',
    'rowanberry, rowanberry powder, dried': 'rowanberry, dried',
    'rowanberry, sorbus': 'rowanberry',
    'ruby red grapefruit juice blend, ocean spray': 'grapefruit juice blend, with added vitamin C',
    'ruby red grapefruit juice blend, ocean spray, enriched': 'grapefruit juice blend, with added vitamin C, enriched',
    "rudi's, gluten-free bakery": 'sandwich bread, gluten-free',
    'rum vol. %': 'rum, 40 vol%',
    'rutabaga or rutabaga': 'rutabaga',
    'rye bread, dark with high fat seeds and whole': 'rye bread, dark, with oil seed and whole kernel',
    'rye bread, water, flour': 'rye bread, from rye and dark wheat flour',
    'rye bread, wholegrain rye, flour': 'rye bread, wholegrain rye with dark wheat flour, industrial',
    'rye macaroni, salt': 'rye macaroni, boiled, salted',
    'rye, wholegrain flour': 'rye flour, wholegrain',
    'ryebread, water, flour': 'rye bread, from wholegrain rye flour',
    "saint-felicien cheese, from cow's milk": 'saint-felicien cheese',
    "saint-nectaire cheese, from cow's milk": 'saint-nectaire cheese',
    'salad crab': 'crab salad',
    'salad dressing with mayonnaise tomato homemade rhode islandsa': 'rhode island dressing with mayonnaise and tomato, homemade',
    'salad dressing, blue or roquefort cheese dressing': 'salad dressing, blue or roquefort cheese',
    'salad dressing, french dressing, juice': 'salad dressing, french, homemade',
    'salad dressing, mayo mayonnaise dressing, fat-free': 'salad dressing mayonnaise type, fat-free',
    'salad dressing, miracle whip free dressing, fat-free': 'salad dressing, fat-free',
    'salad dressing, russian dressing': 'salad dressing, russian',
    'salad fish': 'fish salad',
    'salad salmon': 'salmon salad',
    'salad shrimp': 'shrimp salad',
    'salad tuna': 'tuna salad',
    'sallad with tuna shellfish mayonnaise dressing lettuce tomato': 'salad with tuna, shellfish, mayonnaise dressing, lettuce and tomato',
    'salmon norwegian fjord, farmed': 'salmon, farmed',
    'salmon norwegian, farmed': 'salmon, farmed, raw',
    'salmon nowegian fjord, farmed': 'salmon farmed norwegian fjord',
    'salmon, red': 'salmon, sockeye',
    'salt bakers with added iodin': "baker's salt with added iodine",
    'salt low sodium losalt': 'salt, low sodium',
    'salt not with iodine': 'salt, not fortified with iodine',
    'salty liqourice': 'salty liquorice',
    'salty liqourice pastille': 'salty liquorice pastille',
    'salty liqourice pastille, unsweetened': 'salty liquorice pastille, unsweetened',
    'sandwich white bread with mince meat patty pickled cucumber beetroot parisersmörgå': 'open sandwich on white bread with minced meat patty, pickled cucumber and beetroot',
    'sandwiches and burger, fish sandwich with tartar sauce and cheese': 'fish sandwich, with tartar sauce and cheese',
    'sandwiches and burger, submarine sandwich on white bread with coldcut': 'submarine sandwich, coldcut on white bread with lettuce and tomato',
    'sandwiches and burger, submarine sandwich on white bread with roast beef': 'submarine sandwich, roast beef on white bread with lettuce and tomato',
    'sandwiches and burger, submarine sandwich on white bread with turkey': 'submarine sandwich, turkey on white bread',
    'santen creamed coconut block': 'creamed coconut, block',
    'sapodilla, skin and removed': 'sapodilla, skin and seeds removed',
    'sapota': 'sapodilla',
    'sardine, in tomatosauce': 'sardine, in tomato sauce',
    'sauce': 'sauce, pasta',
    'sauce bolognese homemade': 'bolognese sauce, homemade',
    'sauce bolognese jar': 'bolognese sauce, jarred, ready-to-eat',
    'sauce butter': 'butter sauce',
    'sauce cheese- based on roux': 'cheese sauce, based on roux',
    'sauce with oil mustard dill for spiced salmon homemade gravlaxsas hovmastarsa': 'mustard and dill sauce for cured salmon, homemade',
    'sauce, chinese cook': 'sauce, chinese cook-in, sweet and sour',
    'sauce, indian cook': 'indian cook-in sauce',
    'sauce, tabasco': 'hot pepper sauce',
    'sausage app 15% fat': 'sausage, about 15% fat',
    'sausage dutch frikandel and onion': 'sausage, dutch frikandel, with sauce and onion',
    'sausage dutch frikandel unprep': 'sausage dutch frikandel, unprepared',
    'sausage excl liver product': 'sausage, excluding liver products',
    'sausage hash, rign bologna': 'sausage hash, ring bologna',
    'sausage incl liver product': 'sausage, including liver products',
    'sausage mcgriddle': 'breakfast sandwich on pancakes with sausage',
    'sausage mcmuffin': 'breakfast muffin sandwich with sausage and cheese',
    'sausage mcmuffin with egg': 'breakfast muffin sandwich with sausage, cheese and egg',
    "sausage pizza, thin 'n crispy crust": 'sausage pizza, thin crispy crust',
    'sausage pizza, ultimate deep dish crust': 'sausage pizza, deep dish crust',
    'sausage roast- vegetarian based on soy/wheat with iron and vit b12, enriched': 'vegetarian roast sausage based on soy/wheat, fortified with iron and vit b12',
    'sausage roast- vegetarian based on soy/wheat, enriched': 'vegetarian roast sausage based on soy/wheat, fortified with iron and vit b12',
    'sausage roast- vegetarian based on soya/wheat with iron and vit b12': 'vegetarian roast sausage based on soy/wheat, fortified with iron and vit b12',
    'sausage salam turkish': 'turkish salami sausage',
    'sausage soup, sausage': 'sausage soup, fresh sausage',
    'sausage soup, siskonmakkara soup': 'sausage soup, fresh sausage',
    'sausage spicy 73-75% meat e.g bratwurst': 'sausage, spicy, 73-75% meat, cooked',
    'sausage varmkrov, fried': 'sausage varmkorv, fried',
    'sausage, egg & cheese mcgriddle': 'breakfast sandwich on pancakes with sausage, egg and cheese',
    'sausage. swisswurst, pork and beef': 'swisswurst, pork and beef',
    'sawo, peel and removed': 'sawo, peel and seeds removed',
    'schar, gluten-free': 'roll, gluten-free, white',
    'schnitzel vegetarian based on milk filled several flavour, enriched': 'schnitzel, vegetarian, milk-based, filled, unprepared, fortified with iron',
    'schnitzel vegetarian based on milk filled several flavours with iron': 'schnitzel, vegetarian, milk-based, filled, unprepared, fortified with iron',
    'schnitzel vegetarian based on milk filled several flavours with iron, enriched': 'schnitzel, vegetarian, milk-based, filled, unprepared, fortified with iron',
    'sea buckthorn, sea buckthorn powder, dried': 'sea buckthorn, dried',
    'sea-buckthornberry': 'sea buckthorn berry',
    'seafood sauce rte instant powder heated with milk water butter smögen': 'seafood sauce from instant powder, prepared with milk, water and butter',
    'seafood, cuttlefish, fried': 'breaded cuttlefish, fried',
    'seafood, prawn, fried': 'breaded prawn, fried',
    'seafood, white fish, fried': 'breaded white fish, fried',
    'seaweed, canadian cultivated emi-tsunomata': 'seaweed, tsunomata',
    'seaweed, canadian cultivated emi-tsunomata, dried': 'seaweed, tsunomata, dried',
    'seaweed, irishmoss': 'seaweed, irish moss',
    'semolina cookies filled with dates nuts mamoul': 'semolina cookie filled with date and nut',
    'semolina porridge with sugar mamonia': 'semolina porridge with sugar',
    'sesame butter, from whole, paste': 'sesame butter, from whole seed, paste',
    'sesame butter, tahini': 'tahini',
    'sesame drink water extract of with husk': 'sesame drink, water extract of seed with husk',
    'sesame drink water extract of without husk': 'sesame drink, water extract of seed without husk',
    'sesame nuggets with soy, hälsans kök': 'sesame nuggets with soy',
    'sesame seed, and hulled, dried': 'sesame seed, hulled, dried',
    'sesame seed, grown': 'sesame seed, organically grown',
    'sesame seed, grown inorganic': 'sesame seed, conventionally grown',
    'sesame seeds with husk': 'sesame, with hull',
    'sesame seeds without husk': 'sesame without hull',
    'sharon': 'persimmon',
    'sharon, kaki': 'persimmon',
    'sharp-snout sea bream, wild': 'sharpsnout sea bream, wild',
    'shea butter or shea oil': 'shea butter',
    'shellfish salad with mussles shrimps button mushrooms dressing homemade vastkustsallad': 'shellfish salad with mussels, shrimps, button mushrooms and dressing, homemade',
    'shellfish salad with mussles shrimps button mushrooms homemade västkustsallad': 'shellfish salad with mussels, shrimps and button mushrooms, homemade',
    'sherry medium vol. %': 'sherry, medium dry, 17 vol%',
    'short crust pastryt, pure butter': 'short crust pastry, pure butter',
    'shortening confectionery, coconut and or palm': 'shortening, confectionery, hydrogenated coconut and palm kernel',
    'shrimp, dutch': 'shrimp, boiled',
    'silverbeet': 'swiss chard',
    'six-spotted cockroach, large': 'six-spotted cockroach',
    'small beer vol. % 2.3': 'small beer, 2.3% alcohol by volume',
    'small citrus fruit clementine mandarin satsuma': 'clementine/mandarin/satsuma',
    'small fruit tomato, blended with disrupted': 'tomato, small fruit, blended with seed, skin and pulp',
    'small locust, inedible parts removed, salted': 'small locust, inedible parts removed, roasted, salted',
    'smart soup, french lentil': 'soup, french lentil',
    'smart soup, greek minestrone': 'soup, greek minestrone',
    'smart soup, indian bean masala': 'soup, indian bean masala',
    'smart soup, moroccan chick pea': 'soup, moroccan chick pea',
    'smart soup, santa fe corn chowder': 'soup, corn chowder',
    'smart soup, thai coconut curry': 'soup, thai coconut curry',
    'smart soup, vietnamese carrot lemongrass': 'soup, vietnamese carrot lemongrass',
    'smokies sausage little cheese': 'little smokies cheese sausage, pork and turkey',
    'snack, balance': 'energy bar',
    'snack, betty crocker fruit roll ups': 'snack, fruit roll, berry flavoured',
    'snack, clif bar': 'energy bar',
    'snack, farley candy': 'snack, fruit snack, with added vitamins',
    'snack, fritolay': 'multigrain snack chip, french onion flavour',
    'snack, m&m': 'pretzel snack with cheddar cheese filling',
    'snack, nutri-grain cereal bar': 'cereal bar with fruit',
    'snack, nutri-grain fruit and nut bar': 'snack bar, fruit and nut',
    'snack, potato chip': 'potato, chip',
    'snack, sunkist': 'fruit roll, strawberry',
    "snackbar, granny's": 'snack bar, fruit and nut',
    "snackwell's devil's food cookie cake, fat-free": "devil's food cookie cake, fat-free",
    'snap pea': 'sugar snap pea',
    'soft drink m sugar without caffeine': 'soft drink with sugar, without caffeine',
    'soft drink with sugar and sweetener 2-<5 cho with caffeine': 'soft drink with sugar and sweetener, 2-<5 g carbohydrate, with caffeine',
    'soft drink with sugar and sweetener 5-<8 g cho with caffeine': 'soft drink with sugar and sweetener, 5-<8 g carbohydrate, with caffeine',
    'soft toffee covered with milk chocolate, dumle': 'soft toffee covered with milk chocolate',
    'softdrink with milk serum and sweetener rivella': 'soft drink with milk serum and sweetener',
    'soup legume based packet, dried': 'soup, legume based, prepared from dried packet',
    'soup meat based packet, dried': 'soup, meat based, prepared from dried packet',
    'soup with chicken pasta corn mexicanasoppa': 'soup with chicken, pasta and corn',
    'soup with pork root vegetables fransk bondsoppa': 'soup with pork and root vegetable',
    'soup with potato mince meat cowboysoppa': 'soup with potato and minced meat',
    'soup with white cabbage minced beef carrot hedvigsoppa': 'soup with white cabbage, minced beef and carrot',
    'soup with white cabbage minced beef nikkaluoktasoppa': 'soup with white cabbage and minced beef',
    'soup, nissin, dried': 'soup, ramen noodle, dried',
    'sour milk, butter milk': 'buttermilk',
    'sour sop, fruit pulp': 'soursop, fruit pulp',
    'sour, natal or money plum': 'sour plum, fruit flesh',
    'soured milk 0.1 fat, fat-free': 'soured milk, 0.1% fat, fat-free',
    'sourplum, flesh': 'sour plum, flesh',
    'southern style chicken biscuit': 'breaded chicken biscuit sandwich',
    'soy bacon bit, giant': 'soy bacon bits',
    'soy bacon, manischewitz': 'soy bacon',
    'soy dessert, flavoured': 'soy dessert, flavoured, with sugar',
    'soy dessert, flavoured, sweetened': 'soy dessert, flavoured, with sugar',
    'soy drink calcium, vitamin': 'soy drink, with added calcium and vitamins',
    'soy drink, aplro soy vanilla': 'soy drink, vanilla',
    'soy drink, calcium': 'soy drink, with added calcium',
    'soy drink, flavoured': 'soy drink, flavoured, with sugar',
    'soy drink, flavoured, sweetened': 'soy drink, flavoured, with sugar',
    'soy drink, rainbow': 'soy drink, organic',
    'soy drink, sugar': 'soy drink, sweetened',
    'soy drink, vitamin': 'soy drink, with added vitamins',
    'soy ice cream non-diary': 'soy ice cream, non-dairy',
    'soy oil': 'soybean oil',
    'soy product, berrysoya': 'soy product with berries',
    'soy product, go': 'soured soy product, flavoured',
    'soy protein bits with thyme garlic product ®': 'soy protein bits with thyme and garlic',
    'soy protein concentrate, produced by acid wash': 'soy protein concentrate, with acid and water wash',
    'soy protein concentrate, produced by alcohol extraction': 'soy protein concentrate, with alcohol',
    'soy protein isolate, potassium type': 'soy protein isolate, with potassium',
    'soy protein kebab product ®, fried': 'soy protein kebab, fried',
    'soy protein pulled soy product ®': 'pulled soy protein, frozen',
    'soy protein pulled soy product ®, fried': 'pulled soy protein, fried',
    'soy turkey, ive': 'soy turkey slices',
    'soy, bacon bit': 'soy bacon bits',
    'soy, burger': 'soy burger',
    'soy, cheese': 'soy cheese',
    'soy, commercially': 'soy drink, unfortified',
    'soy, flake': 'soy flakes',
    'soy, grits': 'soy grits',
    'soy, meat': 'soy meat substitute',
    'soy, milk powder': 'soy milk powder',
    'soy, sausage': 'soy sausage',
    'soy, tempe': 'tempeh',
    'soy, tofu': 'tofu, fermented',
    'soy, yoghurt': 'soy yoghurt',
    'soya bean, anidaso variety, dried': 'soybean, dried',
    'soya bean, combined variety, dried': 'soybean, dried',
    'soya bean, dried': 'soybean, dried',
    'soya bean, jenguma variety, dried': 'soybean, dried',
    'soya bean, quarshie variety, dried': 'soybean, dried',
    'soya flour, fat-free': 'soy flour, fat-free',
    'soybean curd': 'tofu',
    'soybean flour, arrowhead mills brand': 'soybean flour',
    'soybean flour, bulk brand': 'soybean flour',
    'soybean protein, bulk brand': 'soy protein',
    'soybean protein, textured vegetable': 'textured vegetable soy protein',
    'soyghurt plain ca vitd vitb12': 'soyghurt, fortified with calcium, vitamin d and vitamin b12',
    'soyghurt plain ca vitd vitb12, enriched': 'soyghurt, fortified with calcium, vitamin d and vitamin b12',
    'soymilk, pacific soy brand': 'soymilk',
    'soymilk, westsoy brand': 'soymilk',
    'soynut, gensoy brand': 'soynut, raw',
    'soynut, good sense brand': 'soynut, roasted',
    'spagetti bolognese, hk': 'spaghetti bolognese, ready meal',
    'spam, canned': 'canned luncheon meat, pork',
    'spanish fish, large': 'spanish fish',
    'spanish rice mix, mix, dried': 'spanish rice mix, dried',
    'spead pro-activ calorie light, low-fat': 'spread with plant sterols, reduced calorie, low-fat',
    'spelt flour, wholemeal': 'spelt flour, wholegrain',
    'spelt wheat flour, whole': 'spelt flour, wholegrain',
    'spice biscuit sprinkles bolletje': 'spiced biscuit sprinkles',
    'spice, poppy': 'poppy seed',
    'spicemix taco': 'spice mix, taco',
    'spider flower, just leaf': 'spider flower, mature leaf',
    'spinach product': 'spinach, frozen',
    'spiny amaranth, just leaf': 'spiny amaranth, mature leaf',
    'spit- pork and vegetable': 'spit-roasted pork and vegetable',
    'sports drink aa isotone': 'sports drink, isotonic',
    'sports drink aquarius': 'sports drink, isotonic',
    'sports drink extran energy': 'sports drink with carbohydrate',
    'sports drink extran hydro': 'sports drink, hypotonic',
    'spread ah omega-3, low-fat': 'spread with omega-3, low-fat',
    'spread pro-activ, low-fat': 'spread with plant sterols, low-fat',
    'spread, 20% butter': 'spread, 20% butter and 80% canola oil',
    'spreadable fat, with plant strol, unsalted': 'spreadable fat 41% fat or less, with plant sterols',
    'squash, melonnette jaspée from vendée': 'squash, melonnette jaspée',
    'squid or calamari': 'squid',
    'squid or calamari, no added fat, fried': 'squid, no added fat, fried',
    "st. john's bread powder, carob powder": 'carob powder',
    'st. paulin cheese': 'saint-paulin, firm cheese',
    'star apple or milk truit, flesh': 'star apple or milk fruit, flesh',
    'star fruit': 'carambola',
    'steak & cheese sub on white bread with american cheese, lettuce and tomato': 'submarine sandwich, steak and cheese on white bread with cheese, lettuce and tomato',
    'still soft drink with tea extract, sugar and sweetener contents unknown': 'still soft drink with tea extract, sugar and sweetener',
    'stock, liquid, dried': 'stock, liquid, prepared from powder or cube',
    'strawberries with sugar product': 'strawberry, sweetened, frozen',
    'submarine sandwich, steak and cheese on white bread with cheese': 'submarine sandwich, steak and cheese on white bread with cheese, lettuce and tomato',
    'suffeli chocolate bar, waffle': 'chocolate wafer bar with toffee filling',
    'suffeli puffi snack, puffed corn and chocolate flavored coating': 'puffed corn snack with chocolate flavoured coating',
    'sugar castor brown': 'caster sugar, brown',
    'sugar castor white': 'caster sugar, white',
    'sugar granulated': 'sugar',
    'sugar powdered': 'sugar, icing',
    'sugar, sucrose': 'sugar',
    'sultana': 'raisin, golden',
    'sun country, kretschmer honey crunch wheat germ': 'wheat germ, honey crunch',
    'sun country, kretschmer wheat bran': 'wheat bran, toasted',
    'sun country, kretschmer wheat germ': 'wheat germ',
    'sunflower': 'sunflower seed',
    'sunflower oil ho, refined': 'sunflower oil, high oleic, refined',
    'supreme pizza, hand-tossed crust': 'pizza with meat and vegetable topping, hand-tossed crust',
    'surmullet or red mullet': 'mullet, red',
    'swede, unsalted': 'rutabaga',
    'swedish brown beans rte': 'swedish brown bean, ready-to-eat',
    'swedish brown beans rte product': 'swedish brown bean, ready-to-eat',
    'swedish cheese cake': 'swedish cheese cake, 7% fat',
    'swedish flatbread with pastrami potatoes lettuce mayonnaise tunnbrödsrulle': 'swedish flatbread roll with pastrami, potatoes, lettuce and mayonnaise',
    'swedish flatbread with sausage mashed potatoes shrimp- cucumber mayonnaise tunnbrödrulle': 'swedish flatbread roll with sausage, mashed potato, shrimp, cucumber and mayonnaise',
    'swedish hash product': 'swedish hash, cooked, frozen',
    'swedish punch arrack vol. %': 'swedish arrack punch, 26 vol%',
    "sweet & sour sauce, uncle ben's": 'sweet and sour sauce, jarred',
    'sweet butter': 'maple butter',
    'sweet crunch/quisp': 'cereal, ready-to-eat, sweetened corn puff',
    'sweet pepper av': 'pepper, sweet',
    'sweet potato, after baking': 'sweet potato, baked, skin removed',
    'sweet potato, for 30 sec': 'sweet potato, microwaved for 30 sec',
    'sweet potato, for 45 sec': 'sweet potato, microwaved for 45 sec',
    'sweet potato, for 60 sec': 'sweet potato, microwaved for 60 sec',
    "sweet potato, or bbq'd, fried": 'sweet potato, white flesh, cooked, no added fat',
    'sweet potato, purple flesh type': 'sweet potato, purple flesh',
    'sweet potato, salt, flesh': 'sweet potato, salted',
    'sweet potatoe': 'sweet potato',
    'sweet yeast sough for pie, whole milk': 'sweet yeast dough for pie, whole milk',
    'sweet, baking chocolate': 'baking chocolate, milk chocolate mini baking bits',
    'sweet, candy': 'boiled sweets',
    'sweet, chocolate': 'milk chocolate',
    'sweet, cocoa, dried': 'cocoa powder, unsweetened',
    "sweet, confectioner's coating or chip": "confectioner's coating chips, peanut butter",
    'sweet, fancy molasses': 'molasses, fancy',
    'sweet, fruit butter': 'apple butter',
    'sweet, fruit pectinbased': 'fruit sweets, pectin-based',
    'sweet, gelatin, dried': 'gelatin, unsweetened, dry powder',
    'sweet, honey': 'honey',
    'sweet, jams and preserve': 'jam, apricot',
    'sweet, marshmallow type spekkie': 'marshmallow sweets',
    'sweetener, aspartam': 'sweetener, aspartame and acesulfame K',
    'sweetener, aspartam, dried': 'sweetener, aspartame, dried',
    'sweetener, canderel, dried': 'sweetener, aspartame and acesulfame-K, dried',
    'sweetener, hermesetas gold': 'sweetener, saccharin, aspartame and cyclamate blend',
    'sweetener, hermesetas liquid': 'sweetener, saccharin and cyclamate blend, liquid',
    'swiss chard, varcicla': 'swiss chard',
    'swiss potato fritter rösti heated product': 'swiss potato fritter rösti, heated, frozen',
    'swiss roll, napakymppi': 'chocolate swiss roll',
    'syrup apple rinse': 'apple syrup',
    'syrup keukenstroop': 'syrup, kitchen',
    'syrup, golden': 'golden syrup',
    'taco flavoured minced meat unspec': 'taco-flavoured minced meat, unspecified',
    'taco-flavoured pork mince mvu': 'taco-flavoured pork mince',
    'tahini, sesame seed pulp': 'tahini',
    'tamarillo fruit': 'tamarillo',
    'tasteeo': 'cereal, ready-to-eat, toasted oat ring',
    'tea herbal instant sweetend': 'herbal tea, instant, sweetened, prepared',
    'tea herbal instant sweetend powder': 'herbal tea, instant, sweetened, dried',
    'tea, black and green': 'tea, black and green, ready-to-drink, lemon, diet',
    'tea, chinese': 'tea, chinese, infusion',
    'tempeh, turtle islan brand': 'tempeh',
    'tempeh, white wave brand': 'tempeh',
    'terapy bean': 'tepary bean, black',
    'the company, glaceau vitamin water': 'vitamin water, fruit punch flavoured, fortified',
    "the company, hi-c flashin' fruit punch": 'drink, fruit punch',
    'the company, nos energy drink': 'energy drink, with sugar',
    'the horned melon, flesh, sweetened': 'horned melon, flesh, sweetened',
    'the perfect, crispy six wholegrain + four cracker': 'cracker, six wholegrain and four seed',
    'the works pizza, original crust': 'pizza with meat and vegetable topping, thin crust',
    'thin pancakes with nutella': 'thin pancakes with chocolate hazelnut spread',
    'thinly- or vegetable': 'vegetables, shredded or diced, frozen',
    "tikka masala, uncle ben's": 'tikka masala sauce, jarred',
    'tilsiter milk': 'cheese, tilsiter, from raw milk',
    'toast melba other variety': 'toast melba',
    'toddler drink, puramino toddler powder': 'toddler formula, amino acid based, dried',
    'toddler formula combiotik': 'toddler formula, organic',
    'toddler formula groeimelk': 'toddler growing-up milk',
    'toddler formula little steps': 'toddler formula',
    'toddler formula nutrasense': 'toddler formula',
    'toddler formula, enfagrow': 'toddler formula',
    'toddler formula, nutramigen toddler with lgg powder': 'hydrolysed toddler formula with probiotics, powder',
    'tofu bean curd': 'tofu',
    'tofu ferm': 'tofu, firm, with calcium salt',
    'tofu, and fermented, salted': 'tofu, salted and fermented',
    'tofu, ferm': 'tofu, firm, with nigari',
    'tofu, from germinated': 'tofu, from germinated soybean',
    'tofu, silken style': 'tofu, silken',
    'tofu, silky': 'tofu, silken',
    'tofu, soya bean': 'tofu, steamed',
    'tofu, soybean curd': 'tofu',
    'tofu, with calcium sulfate, fried': 'tofu, with calcium sulphate, fried',
    'tom collin, gin cocktail': 'tom collins, gin cocktail',
    'tomato av': 'tomato',
    'tomato flesh': 'tomato, flesh, canned',
    'tomato red cherry tomato': 'tomato, cherry, red',
    'tootie fruity': 'cereal, ready-to-eat, fruit flavoured corn and oat ring',
    'topic/snickers and own brand equivalent': 'chocolate bar with nuts and caramel',
    "topping, smucker's magic shell": 'topping, chocolate',
    'traditional confectionery, isobe-senbei, flour': 'isobe-senbei, soft wheat flour cracker',
    'traditional confectionery, kawara senbei, flour': 'kawara senbei, hard wheat flour cracker',
    'traditional confectionery, nanbu-senbei, flour': 'nanbu-senbei, round wheat flour cracker with peanut',
    'tradtion french bread': 'french bread, traditional',
    'treacle or molasses': 'molasses',
    'tuber, dioscorea pentaphylla l. var. pentaphylla': 'five-leaf yam, tuber',
    'tuber, momordica diocia': 'spine gourd, tuber',
    'tuber, nymphaea pubescen': 'pink water lily, tuber',
    'tuber, nymphaea rubra': 'red water lily, tuber',
    'tulumba tatlisi eclair turkish, fried': 'tulumba tatlisi, fried',
    'tumbleweed, just leaf': 'tumbleweed, mature leaf',
    'tuna, ahi or yellowfin': 'tuna, yellowfin',
    'turkey breakfast sausage, butterball': 'turkey breakfast sausage',
    'turkey breakfast sausage, hillshire smoked': 'turkey breakfast sausage, smoked',
    'turkey breakfast sausage, honeysuckle': 'turkey breakfast sausage, fresh',
    'turkey breakfast sausage, honeysuckle white': 'turkey breakfast sausage, white meat',
    'turkey breakfast sausage, jennie o - mild': 'turkey breakfast sausage, mild',
    'turkey breakfast sausage, jimmy dean': 'turkey breakfast sausage',
    'turkey, brest, fried': 'turkey, breast, fried',
    'turkey, brest, ground, fried': 'turkey, breast, ground, fried',
    'turkey, rotisserie': 'turkey, white meat, rotisserie deli cut',
    'turnip top': 'turnip greens',
    'turnip, for 30 sec': 'turnip, microwaved for 30 sec',
    'turnip, for 45 sec': 'turnip, microwaved for 45 sec',
    'turnip, for 60 sec': 'turnip, microwaved for 60 sec',
    'twizzlers strawberry twists candy': 'strawberry twist candy',
    'ultimate chicken grill sandwich': 'grilled chicken sandwich, takeaway',
    'uncle sam cereal': 'cereal, ready-to-eat, wheat berry with flaxseed',
    'usa azumaya, firm tofu': 'tofu, extra firm',
    'usa azumaya, silken tofu': 'tofu, silken',
    'usa sprouted, tofu plus firm': 'tofu, super firm, sprouted',
    'usa, firm tofu': 'tofu, extra firm, organic',
    'usa, lite firm tofu': 'tofu, firm, low-fat',
    'usa, lite silken tofu': 'tofu, silken, low-fat',
    'usa, silken tofu': 'tofu, silken, organic',
    'usa, soft tofu': 'tofu, soft, organic',
    'usa, tofu plus firm': 'tofu, extra firm',
    'v. aconitifolia': 'moth bean',
    'v. ambacensis, cream, dried': 'vigna ambacensis seed, cream, dried',
    'v. luteola, brown mottled, dried': 'hairypod cowpea, brown mottled, dried',
    'v. luteola, dried, brown': 'hairypod cowpea, brown, dried',
    'v. oblongifolia, brown-black, dried': 'vigna oblongifolia, brown-black, dried',
    'v. oblongifolia, dried, brown': 'vigna oblongifolia, brown, dried',
    'v. racemosa, brown mottled, dried': 'vigna racemosa, brown mottled, dried',
    'v. racemosa, light brown, dried': 'vigna racemosa, light brown, dried',
    'v. reticulata, brown mottled, dried': 'vigna reticulata, brown mottled, dried',
    'v. reticulata, light brown, dried': 'vigna reticulata, light brown, dried',
    'v. unguiculata dekindtiana, cream, dried': 'wild cowpea seed, cream, dried',
    'v. unguiculata dekindtiana, dried, brown': 'wild cowpea, brown, dried',
    'v. unguiculata dekindtiana, light brown, dried': 'wild cowpea, light brown, dried',
    'v. unguiculata, dried, black': 'cowpea, black, dried',
    'v. unguiculata, maroon, mature, dried': 'cowpea, maroon, dried',
    'v. vexillata': 'zombi pea',
    'v. vexillata macrosperma, green-brown, dried': 'zombi pea, green-brown, dried',
    'v. vexillata, dried': 'zombi pea, dried',
    'vanilla sauce, wheat flour, low-fat': 'vanilla sauce with low-fat milk, thickened with potato flour',
    'vanilla, light': 'soft-serve ice cream cone, vanilla, light',
    'vanillacrème': 'vanilla cream',
    'vanille cream, homemade': 'vanilla cream, homemade',
    'veal av': 'veal',
    'veal schnitzel not': 'veal schnitzel, not breaded',
    'veal stock paste or powder large-scale': 'veal stock, paste or powder',
    'veal stock paste or powder with large-scale, low-salt': 'veal stock paste or powder with reduced salt',
    'veal stock paste or powder with reduced salt large-scale': 'veal stock paste or powder with reduced salt',
    'veal tenderloin medaillon': 'veal tenderloin medallion',
    'vegatable oil product 82%, culinesse': 'vegetable oil product 82%, fluid, for cooking',
    'vegatable oil product 82%, fluid': 'vegetable oil product 82%, fluid',
    'vegetable fat spread 28%, mini': 'vegetable fat spread, 28% fat',
    'vegetable fat spread 30%, kevyt': 'vegetable fat spread, 30% fat, light',
    'vegetable fat spread 35%, pro-activ': 'vegetable fat spread, 35% fat, with plant sterols',
    'vegetable fat spread 39 %, original': 'vegetable fat spread, 39% fat',
    'vegetable fat spread 55 %, rypsi': 'vegetable fat spread, 55% fat, rapeseed oil',
    'vegetable fat spread 70%, gold': 'vegetable fat spread, 70% fat',
    'vegetable fat spread 75 %, rypsi': 'vegetable fat spread, 75% fat, rapeseed oil',
    'vegetable mix with peas beans carrot cauliflower product, salted': 'vegetable mix with peas, beans, carrot and cauliflower, boiled, salted',
    'vegetable oil, almond': 'almond oil',
    'vegetable oil, apricot': 'apricot oil',
    'vegetable oil, avocado': 'avocado oil',
    'vegetable oil, babassu': 'babassu oil',
    'vegetable oil, canola': 'rapeseed oil',
    'vegetable oil, cocoa butter': 'cocoa butter',
    'vegetable oil, coconut': 'coconut oil',
    'vegetable oil, corn': 'corn oil',
    'vegetable oil, corn and canola': 'corn and canola oil',
    'vegetable oil, cottonseed': 'cottonseed oil',
    'vegetable oil, grapeseed': 'grapeseed oil',
    'vegetable oil, hazelnut': 'hazelnut oil',
    'vegetable oil, mustard': 'mustard oil',
    'vegetable oil, oat': 'oat oil',
    'vegetable oil, olive': 'olive oil',
    'vegetable oil, palm': 'palm oil',
    'vegetable oil, peanut': 'peanut oil',
    'vegetable oil, poppyseed': 'poppyseed oil',
    'vegetable oil, rice bran': 'rice bran oil',
    'vegetable oil, safflower': 'safflower oil',
    'vegetable oil, sesame': 'sesame oil',
    'vegetable oil, soybean': 'soybean oil',
    'vegetable oil, soybean lecithin': 'soybean lecithin oil',
    'vegetable oil, sunflower': 'sunflower oil, mid-oleic',
    'vegetable oil, walnut': 'walnut oil',
    'vegetable oil, wheat germ': 'wheat germ oil',
    'vegetable patties caribbean with carrot brown rice chickpeas broccoli': 'vegetable patties, caribbean, with carrot brown rice chickpeas broccoli, fried',
    'vegetable stock paste or powder large-scale': 'vegetable stock paste or powder',
    'vegetable, pickled turkish tursu': 'pickled vegetable, turkish tursu',
    'vegetables av': 'vegetable',
    'vegetables for stir-frying dutch': 'vegetables for stir-frying',
    'vegetarian based on soy/wheat with iron and vit b12, enriched': 'vegetarian chunk based on soy/wheat, fortified with iron and vit b12',
    'vegetarian based on soy/wheat, enriched': 'vegetarian chunk based on soy/wheat, fortified with iron and vit b12',
    'vegetarian based on soya/wheat with iron and vit b12': 'vegetarian chunk based on soy/wheat, fortified with iron and vit b12',
    'vegetarian sausage with soy and wheat protein heated e.g middagskorv': 'vegetarian sausage with soy and wheat protein, heated',
    'vegetarian sausage with sunflower protein and pea protein heated or product': 'vegetarian sausage with sunflower protein and pea protein, heated',
    'vegetarian sausage with sunflower protein and pea protein or product': 'vegetarian sausage with sunflower protein and pea protein',
    'velveeta light process cheese product, low-fat': 'process cheese product, low-fat',
    'velveeta process cheese spread': 'cheese spread, process',
    'verdolaga, leaf': 'purslane, leaf',
    'vernonia, leaf': 'bitter leaf',
    'vigna aconitifolia, mature, dried': 'moth bean, dried',
    'vigna radiata, mature, dried': 'mung bean, dried',
    'vigna umbellata, mature, dried': 'rice bean, dried',
    'vigna unguiculata, mature, dried': 'cowpea, dried',
    'vigna vexillata, mature, dried': 'zombi pea, dried',
    'vine spinach': 'malabar spinach',
    'vinespinach': 'malabar spinach',
    'vitamin water, vitamin well, sweetened': 'vitamin water, sweetened with fructose',
    'voimamaitojuoma milk, fat-free': 'milk protein drink, lactose-free, fat-free',
    'wafer with caramel milk puffed rice chocolate couverture lion': 'wafer with caramel, puffed rice and chocolate couverture',
    'waffle crisp': 'cereal, ready-to-eat, maple waffle flavoured corn',
    'waffle luikse with chocolate': 'liege waffle with chocolate',
    'walnut, blue diamond grower': 'walnut, raw',
    'walnut, english': 'walnut',
    'water 0-50 mg calcium p litre': 'water, 0-50 mg calcium per litre',
    'water 50-100 mg calcium p litre': 'water, 50-100 mg calcium per litre',
    'water >100 mg calcium p litre': 'water, >100 mg calcium per litre',
    'water convolvulus': 'water spinach',
    'water cress': 'watercress',
    'water cress, leaf': 'watercress, leaf',
    'water melon, dried': 'watermelon, dried',
    'waterapple, skin and removed': 'water apple, skin and seeds removed',
    'waterleaf leaf, wild': 'waterleaf, leaf, wild',
    'wax beans product': 'wax bean, frozen',
    'weetabix wholegrain cereal': 'cereal, ready-to-eat, wholegrain wheat biscuit',
    'weighted average': 'fish',
    'weighted average, fried': 'fish, fried',
    "wend'y, crispy chicken sandwich": 'crispy chicken sandwich',
    'wheat': 'cereal, ready-to-eat, shredded wheat, lightly frosted',
    'wheat bran, extruded at 115ºc and 20%h20': 'wheat bran, extruded at 115ºc and 20% h2o',
    'wheat bread av brown and wholemeal': 'wheat bread, average of brown and wholemeal',
    'wheat bread wholemeal av fine and coarse': 'wholemeal wheat bread, average of fine and coarse',
    'wheat bread wholemeal vollerkoren, enriched': 'wheat bread, wholemeal, fortified with vitamin D and fibre',
    'wheat bread wholemeal with vit d and fibre vollerkoren': 'wheat bread, wholemeal, fortified with vitamin D and fibre',
    'wheat bread wholemeal with vit d and fibre vollerkoren, enriched': 'wheat bread, wholemeal, fortified with vitamin D and fibre',
    'wheat bread, large': 'wheat bread, coarse meal',
    'wheat flour wholegrain graham flour': 'wheat flour, wholegrain graham',
    'wheat, bagged cereal': 'cereal, ready-to-eat, shredded wheat',
    'wheat, honey': 'granola with oat, wheat, honey and raisin',
    'wheatena, dried': 'hot wheat cereal, dried',
    'wheatena, with water': 'toasted wheat cereal, with water',
    'wheatena, with water, salted': 'hot wheat cereal, cooked with water, salted',
    'wheatgerm oil': 'wheat germ oil',
    'whey powder. for food industry': 'whey powder, for food industry',
    'whey protein powder isolate': 'whey protein isolate powder',
    'whiskey sour mix, bottled': 'whisky sour mix, bottled',
    'whiskey sour mix, dried': 'whisky sour mix, dried',
    'whisky vol. %': 'whisky, 40 vol%',
    'white bread, french bread, flour': 'bread, french, white',
    'white bread, french bread, flour, low-fat': 'bread, french, white, low-fat',
    'white currant': 'whitecurrant',
    'white currant, white dutch': 'white currant',
    'white shrimp or tropical shrimp': 'shrimp, white',
    'whitecurrant': 'white currant',
    'whiting rå': 'whiting, raw',
    'whole hearts oat cereal': 'cereal, ready-to-eat, wholegrain oat heart',
    'wholewheat bread': 'whole wheat bread',
    'whopper, with cheese': 'hamburger, flame-grilled, with cheese',
    'whopper, without cheese': 'hamburger, flame-grilled, without cheese',
    'wiener': 'frankfurter',
    'wiener, beef': 'frankfurter, beef',
    'wiener, chicken': 'frankfurter, chicken',
    'wiener, turkey': 'frankfurter, turkey',
    'wild raddish, shoots and leaf': 'wild radish, shoots and leaf',
    'wine medicinal pleegzusterbloedwijn': 'medicinal tonic wine',
    'wine red or rosé vol. %': 'wine, red or rosé, 1 vol%',
    'wine red vol. %': 'wine, red, 14 vol%',
    'wine white or rhine wine vol. %': 'wine, white or rhine, 10 vol%',
    'wine white vol. %': 'wine, white, 12 vol%',
    'wine, rose': 'wine, rosé',
    'winged beans leaf': 'winged bean leaf',
    'winged beans tuber': 'winged bean tuber',
    'with crunchberry': 'cereal, ready-to-eat, sweetened corn and oat with berry piece',
    'with mycoprotein': 'plant-based piece, with mycoprotein',
    'without testa, bauhinia petersiana': 'kalahari coffee bean, without testa',
    'wok-vegetable, green bean-carrot-sweet pepper- sweet corn-parsnip': 'wok-vegetable, green bean, carrot, sweet pepper, sweet corn and parsnip',
    'wonder hamburger roll': 'hamburger roll',
    'yandlong bean, pod': 'yardlong bean, pod',
    'yeast extract, marmite': 'yeast extract spread',
    'yello mealworm, large': 'yellow mealworm, large',
    'yello mealworm, medium': 'yellow mealworm',
    'yello mealworm, small': 'yellow mealworm, small',
    'yellow banana, flesh': 'banana, flesh, steamed',
    'yoghurt bulgarian whole milk': 'yoghurt, bulgarian, whole milk',
    'yoghurt cream- with fruit': 'cream yoghurt with fruit',
    'yoghurt drink 7-9g cho': 'yoghurt drink, 7-9 g carbohydrate',
    'yoghurt drink actimel fruit': 'yoghurt drink with probiotics, fruit',
    'yoghurt drink actimel plain': 'yoghurt drink with probiotics',
    'yoghurt drink activia': 'yoghurt drink with probiotics',
    'yoghurt drink ayran': 'ayran, salted yoghurt drink',
    'yoghurt drink fristi with sweetener': 'yoghurt fruit drink, with sweetener',
    'yoghurt drink pro-activ': 'yoghurt drink with plant sterols',
    'yoghurt drink vifit fruit': 'yoghurt drink with probiotics, fruit',
    'yoghurt full-fat plain activia': 'yoghurt, full-fat, with probiotics',
    'yoghurt greek full-fat': 'greek yoghurt, full-fat',
    'yoghurt greek, fat-free': 'greek yoghurt, fat-free',
    'yoghurt natural': 'yoghurt',
    'yoghurt or fermented milk, with chocolate flake': 'yoghurt or fermented milk, with chocolate flake, with cream, with sugar',
    'yoghurt or fermented milk, with chocolate flake, sweetened': 'yoghurt or fermented milk, with chocolate flake, with cream, with sugar',
    'yoghurt snack breaker': 'yoghurt snack pot',
    'yoghurt substitute, soy-based. with fruit or flavour': 'yoghurt substitute, soy-based, with fruit or flavour',
    'yoghurtor fermented milk, plain or with fruit': 'yoghurt or fermented milk, plain or with fruit',
    'yogurt, chocolate': 'frozen yoghurt, chocolate',
    'yogurt, chocolate, sweetened': 'frozen yoghurt, chocolate, sweetened',
    'yogurt, chocolate, sweetened, fat-free': 'frozen yoghurt, chocolate, sweetened, fat-free',
    'yogurt, flavors other than chocolate': 'frozen yoghurt, flavours other than chocolate',
    'yogurt, vanilla': 'frozen yoghurt, vanilla, soft-serve',
    'zander, pike-perch': 'pike-perch',
    'zevia, cola': 'cola, without added sugars and with sweetener',
    'zizyphus': 'jujube',
    'zucchini or zucchini, flesh and skin': 'zucchini, flesh and skin',
    'zucchini or zucchini, puree': 'zucchini, puree',
    'zucchini, in rapeseed oil': 'zucchini, roasted in rapeseed oil',
}



# ---------------------------------------------------------------------------
# CROSS-DATABASE CANON MERGES  (generated, then reviewed; see the note below)
#
# Sixteen national databases write the same food in different orders, so one
# food arrived as several canons and its evidence was split: 'olive oil'
# (fineli, frida, japan, swedish, swiss, wafct), 'oil, olive' (afcd, fdc,
# mccance) and 'oil olive' (nevo) were three separate foods. This table folds
# 982 such canons into 848 groups.
#
# Built by grouping canons whose words match once punctuation and word order
# are ignored, then picking as target the form the MOST source databases
# already use - so the direction is data-driven and churn is minimised.
#
# Three guards, each added after a real false positive was caught in review:
#   * packing medium is kept apart from the food: 'pineapple juice' (the
#     juice) must not meet 'pineapple, in juice' (fruit packed in juice);
#   * comparison operators are kept apart: 'dark chocolate cocoa < 70%' had
#     merged into '> 70%', its opposite;
#   * two groups are denied outright, because word order carries the meaning:
#     chocolate milk (a drink) vs milk chocolate (a confectionery), and
#     beef roast (a cut) vs roast beef (a cooked dish);
#   * a slashed dual common name never merges with the spaced form:
#     'peach/nectarine' is the African-DB synonym convention and the
#     slash is meaningful, not punctuation to be normalised away.
# ---------------------------------------------------------------------------
_CANON_MERGES = {
    'all-purpose flour': 'flour, all purpose',
    'all-purpose flour, enriched': 'flour, all purpose, enriched',
    'amaranth, leaf': 'amaranth leaf',
    'amaranth, leaf, unsalted': 'amaranth leaf, unsalted',
    'apple pie': 'pie, apple',
    'apple strudel': 'strudel, apple',
    'archway home style cookies oatmeal': 'archway home style cookie, oatmeal',
    'artichoke, heart': 'artichoke heart',
    'artichoke, jerusalem': 'jerusalem artichoke',
    'arugula, arugula': 'rocket, arugula',
    'asparagus green': 'asparagus, green',
    'asparagus green, salted': 'asparagus, green, salted',
    'asparagus soup': 'soup, asparagus',
    'asparagus white': 'asparagus, white',
    'atlantic halibut': 'halibut, atlantic',
    'avocado fruit': 'avocado, fruit',
    'baby carrot': 'carrot, baby',
    'babyfood, fruit and vegetable': 'babyfood, vegetables and fruit',
    'bacon rasher, back, low-salt': 'bacon, back, low-salt',
    'bacon rasher, streaky': 'bacon rasher streaky',
    'bacon, middle rasher': 'bacon rasher, middle',
    'bacon, middle rasher, fried': 'bacon rasher, middle, fried',
    'bacon, turkey': 'turkey bacon',
    'bacon, turkey, low-salt': 'turkey bacon, low-salt',
    'bagel, with egg': 'bagel, egg',
    'banana bread': 'bread, banana',
    'baobab, fruit': 'baobab fruit',
    'bar, milk and cereal': 'milk and cereal bar',
    'bar, muesli': 'muesli bar',
    'barbecue sauce': 'sauce, barbecue',
    'barbeque sausage, cheese': 'barbeque sausage with cheese',
    'barley flake': 'flake, barley',
    'barley malt flour': 'barley, malt flour',
    'barley pearl': 'barley, pearl',
    'barley wholegrain': 'barley, wholegrain',
    'barley, flake': 'flake, barley',
    'barley, pearled': 'pearled barley',
    'barley, wholegrain flour': 'barley flour, wholegrain',
    'bay leaf': 'bay, leaf',
    'bay leaf, dried': 'bay, leaf, dried',
    'bean sprout, mung': 'mung bean sprout',
    'bean, broad': 'broad bean',
    'bean, broad, dried': 'broad bean, dried',
    'bean, broad, unsalted': 'broad bean, unsalted',
    'bean, brown, dried': 'brown bean, dried',
    'bean, brown, dried, 0% moisture basis': 'brown bean, dried, 0% moisture basis',
    'bean, chili': 'chili with bean',
    'bean, green': 'green bean',
    'bean, green, dried': 'green bean, dried',
    'bean, haricot, dried': 'bean, navy, dried',
    'bean, hyacinth': 'hyacinth bean',
    'bean, hyacinth, salted': 'hyacinth bean, salted',
    'bean, in tomato sauce, low-salt': 'baked bean, canned in tomato sauce, low-salt',
    'bean, lima': 'lima bean',
    'bean, lima, dried': 'lima bean, dried',
    'bean, lima, dried, salted': 'lima bean, dried, salted',
    'bean, lima, salted': 'lima bean, salted',
    'bean, mung': 'mung bean',
    'bean, mung, dried': 'mung bean, dried',
    'bean, mung, dried, unsalted': 'mung bean, dried, unsalted',
    'bean, mung, salted, sprouted': 'mung bean, salted, sprouted',
    'bean, mung, sprouted': 'mung bean, sprouted',
    'bean, red kidney': 'bean, kidney, red',
    'bean, small white': 'bean, small, white',
    'bean, soy, dried': 'soybean, dried',
    'bean, sword': 'sword bean',
    'bean, wax': 'wax bean',
    'bean, winged, dried': 'winged bean, dried',
    'bean, winged, dried, salted': 'winged bean, dried, salted',
    'bean, yardlong': 'yardlong bean',
    'bean, yardlong, salted': 'yardlong bean, salted',
    'beef gravy': 'gravy, beef',
    'beef stew': 'stew, beef',
    'beef stew with potatoes and vegetable': 'vegetable and beef stew with potato',
    'beer, white': 'beer white',
    'beet, salted': 'beetroot, salted',
    'bell pepper green': 'bell pepper, green',
    'bell pepper red': 'bell pepper, red',
    'berry pie, double-crust': 'berry pie, double crust',
    'berry pie, shortbread crust with sour milk': 'berry pie with sour milk, shortbread crust',
    'betel, leaf': 'betel leaf',
    'biscuit digestive': 'biscuit, digestive',
    'biscuit sponge finger': 'biscuit, sponge finger',
    'biscuit with nuts and chocolate': 'biscuit, with chocolate and nut',
    'biscuit, wafer': 'wafer biscuit',
    'biscuit, with chocolate': 'biscuit with chocolate',
    'biscuit, with chocolate and with sweetener, enriched': 'biscuit with chocolate, enriched',
    'biscuit, with egg and sausage': 'sausage biscuit with egg',
    'biscuit, with sausage': 'sausage biscuit',
    'biscuit, with sesame': 'sesame biscuit',
    'bitterleaf, leaf, salted': 'bitter leaf, salted',
    'black bean': 'bean, black',
    'black bean soup': 'soup, black bean',
    'black gram, salted': 'bean, mungo, salted',
    'black pepper': 'pepper, black',
    'black tea, brewed': 'tea, brewed, black',
    'blackberries product, sweetened': 'blackberry, frozen, sweetened',
    'blueberry pie': 'pie, blueberry',
    'bouillon, meat': 'meat bouillon',
    'bran, wheat': 'wheat bran',
    'bread crumb, dried': 'bread, crumb, dried',
    'bread roll': 'bread, roll',
    'bread white': 'white bread',
    'bread white wheat': 'white bread, wheat',
    'bread wholegrain rye': 'rye bread, wholegrain',
    'bread with sourdough gluten-free': 'gluten-free bread sourdough',
    'bread, brown': 'brown bread',
    'bread, garlic and herb': 'bread, garlic or herb',
    'bread, rye': 'rye bread',
    'bread, sprouted wheat': 'bread, wheat, sprouted',
    'bread, wheat': 'wheat bread',
    'bread, white': 'white bread',
    'bread, white wheat': 'white bread, wheat',
    'bread, whole wheat': 'whole wheat bread',
    'bread-crumb': 'bread crumb',
    'bread/roll': 'bread, roll',
    'breakfast cereal corn flake': 'breakfast cereal, corn flake',
    'breakfast cereal cornflake': 'breakfast cereal, cornflake',
    'breakfast cereal puffed rice': 'breakfast cereal, puffed rice',
    'breakfast cereal rice krispy': 'breakfast cereal, rice krispy',
    'breakfast cereal wholegrain special flake, enriched': 'breakfast cereal, wholegrain flake, enriched',
    'breakfast cereal, with bran': 'breakfast cereal, bran',
    'breakfast tart, low-fat': 'tart, breakfast, low-fat',
    'broccoli product, salted': 'broccoli, frozen, salted',
    'brown sugar': 'sugar, brown',
    'brownie, with nut': 'brownie with nut',
    'buckwheat, groat': 'buckwheat groat',
    'bulgur wheat': 'wheat, bulgur',
    'burger, chicken': 'chicken burger',
    'burrito, with bean': 'bean burrito',
    'burritos with beef': 'burrito with beef',
    'butter, clarified': 'clarified butter',
    'butternut squash': 'squash, butternut',
    'cabbage chinese': 'cabbage, chinese',
    'cabbage green': 'cabbage, green',
    'cabbage red': 'red cabbage',
    'cabbage red with apple': 'cabbage, with apple, red',
    'cabbage white': 'cabbage, white',
    'cabbage, mustard': 'mustard cabbage',
    'cabbage, mustard, salted': 'mustard cabbage, salted',
    'cabbage, napa': 'napa cabbage',
    'cabbage, red': 'red cabbage',
    'cabbage, salted, red': 'red cabbage, salted',
    'cake chocolate': 'cake, chocolate',
    'cake mix, pudding-type, white, dried, enriched': 'cake, pudding-type, dried, enriched, white',
    'cake mix, yellow, dried, enriched': 'cake, mix, dried, enriched, yellow',
    'cake plain': 'cake',
    'cake with carrot': 'carrot cake',
    'cake, carrot': 'carrot cake',
    'cake, coffeecake, enriched': 'coffee cake with cheese, enriched',
    'cake, fruit': 'fruit cake',
    'cake, mud': 'mud cake',
    'cake, pudding-type, white, dried, enriched': 'cake mix, pudding-type, dried, enriched',
    'caramel candy': 'candy, caramel',
    'carbonated beverage or fruit soft drink, fruit content < 10%': 'carbonated beverage or fruit soft drink, fruit content <10%',
    'carbonated beverage or fruit soft drink, fruit content < 10%, sweetened': 'carbonated beverage or fruit soft drink, fruit content <10%, sweetened',
    'cashew nut, unsalted': 'nut, cashew, unsalted',
    'cashew, salted': 'cashew nut, salted',
    'catsup, low-salt': 'ketchup, low-salt',
    'cereal bar with milk': 'milk and cereal bar',
    'cereal flake, plain with vitamins and mineral, sweetened, enriched': 'cereal flake, with sugar, enriched',
    'cereal, harvest crunch: raisins with almond': 'cereal, harvest crunch: raisins and almond',
    'cereal, harvest crunch: raisins with almond, low-fat': 'cereal, harvest crunch: raisins and almond, low-fat',
    'cereal, honeycomb': 'honeycomb cereal',
    'cereal, raisin bran': 'raisin bran cereal',
    'channel catfish, wild': 'catfish, channel, wild',
    'chard, salted': 'swiss chard, salted',
    'chayote fruit': 'chayote, fruit',
    'cheese brie': 'cheese, brie',
    'cheese camembert': 'cheese, camembert',
    'cheese cheddar': 'cheese, cheddar',
    'cheese cottage': 'cheese, cottage',
    'cheese cracker': 'cracker, cheese',
    'cheese danish blue': 'cheese, danish blue',
    'cheese edam': 'cheese, edam',
    'cheese feta': 'cheese, feta',
    'cheese goat': 'cheese, goat',
    'cheese gorgonzola': 'cheese, gorgonzola',
    'cheese gouda': 'cheese, gouda',
    'cheese gruyere': 'cheese, gruyere',
    'cheese hard': 'cheese, hard',
    'cheese limburger': 'cheese, limburger',
    'cheese parmesan': 'cheese, parmesan',
    'cheese processed': 'cheese, processed',
    'cheese ricotta': 'cheese, ricotta',
    'cheese roquefort': 'cheese, roquefort',
    'cheese sauce': 'sauce, cheese',
    'cheese stilton': 'cheese, stilton',
    'cheese, american': 'american cheese',
    'cheese, american, fat-free': 'american cheese, fat-free',
    'cheese, blue': 'blue cheese',
    'cheese, goat, hard': 'cheese goat hard',
    'cheese, goats milk': 'cheese, goat milk',
    'cheese, parmesan, store brand/other brand': 'cheese, parmesan',
    'cheese, port-salut': 'cheese, port salut',
    'cheese, processed, in slice': 'cheese, processed, slice',
    'cheese, ricotta, whole milk': 'ricotta cheese whole milk',
    'cheese, whey': 'whey cheese',
    'cherries sweet': 'sweet cherry',
    'chewing gum sugar free': 'chewing gum, sugar-free',
    'chicken and rice soup': 'soup, chicken rice',
    'chicken egg, whole': 'egg whole chicken',
    'chicken gravy': 'gravy, chicken',
    'chicken meat': 'chicken, meat',
    'chicken mince': 'chicken, mince',
    'chicken mince, fried': 'chicken, mince, fried',
    'chicken noodle soup': 'soup, chicken noodle',
    'chicken pizza': 'pizza, chicken',
    'chicken soup': 'soup, chicken',
    'chicken tender': 'chicken, tender',
    'chicken, breast, mesquite flavor, fat-free': 'chicken, breast, fat-free',
    'chicken, broilers or fryer, breast': 'chicken, broiler or fryer, breast',
    'chicken, marinated wing': 'chicken, wing, marinated',
    'chicken, nugget': 'chicken nugget',
    'chicken, sausage': 'chicken sausage',
    'chicken, schnitzel': 'chicken schnitzel',
    'chili bean': 'chili with bean',
    'chinese cabbage': 'cabbage, chinese',
    'chinese chestnut': 'chestnut, chinese',
    'chinese chestnut, dried': 'chestnut, chinese, dried',
    'chocolate bar with nougat, caramel and peanuts': 'chocolate bar with caramel, nougat and peanuts',
    'chocolate biscuit with vanilla filling': 'chocolate biscuit, with vanilla filling',
    'chocolate cake': 'cake, chocolate',
    'chocolate chip cookie': 'cookie, chocolate chip',
    'chocolate cookie': 'cookie, chocolate',
    'chocolate dark': 'chocolate, dark',
    'chocolate frosting': 'frosting, chocolate',
    'chocolate malt powder': 'chocolate malt, powder',
    'chocolate malt, with fat free milk, dried': 'chocolate malt powder, fat-free',
    'chocolate pudding': 'pudding, chocolate',
    'chocolate sandwich cookie': 'cookie, chocolate sandwich',
    'chocolate sponge cake': 'sponge cake, chocolate',
    'chocolate spread with hazelnut': 'chocolate hazelnut spread',
    'chocolate white': 'white chocolate',
    'chocolate, bar': 'chocolate bar',
    'chocolate, milk': 'milk chocolate',
    'chocolate, white': 'white chocolate',
    'chow mein, chicken': 'chicken chow mein',
    'chutney, mango': 'mango chutney',
    'cider, apple': 'cider apple',
    'cocoa powder sweetened': 'cocoa powder, sweetened',
    'cocoa, powder': 'cocoa powder',
    'coconut cream': 'coconut, cream',
    'cod, fillet, fried': 'cod fillet, fried',
    'cod, liver': 'cod liver',
    'coffee brewed': 'coffee, brewed',
    'coffee instant powder': 'instant coffee, powder',
    'coffee whitener, milk fat, dried': 'coffee whitener powder',
    'coffee, beverage': 'coffee beverage',
    'coffee, irish': 'irish coffee',
    'common cabbage, green': 'cabbage, common, green',
    'common cabbage, purple': 'cabbage, common, purple',
    'common cabbage, red': 'cabbage, common, red',
    'common cabbage, white': 'cabbage, common, white',
    'complan powder, sweet': 'nutritional supplement drink powder, sweet, fat-free',
    'condensed milk, sweetened': 'milk, condensed, sweetened',
    'coriander, leaf': 'coriander leaf',
    'coriander, leaf, dried': 'coriander leaf, dried',
    'corn grain': 'corn, grain',
    'corn tortilla': 'tortilla, corn',
    'corn with red and green pepper, solids and liquid': 'corn with red or green pepper, solids and liquid',
    'corn, flake': 'corn flake',
    'corn, on the cob': 'corn, on cob',
    'corn, refined flour': 'corn flour, white',
    'corn, sweet': 'sweet corn',
    'corn, sweet, salted': 'sweet corn, salted',
    'corned-beef': 'corned beef',
    'cottage cheese, low-lactose, fat-free': 'cheese, cottage, fat-free',
    'cowpea leaf, dried': 'cowpea, leaf, dried',
    'crab, brown meat': 'crab brown meat',
    'crab, white meat': 'crab white meat',
    'cracker, cream': 'cream cracker',
    'cracker, meal': 'cracker meal',
    'cracker, standard snack type': 'cracker, standard snack-type',
    'cracker, standard snack type, low-salt': 'cracker, standard snack-type, low-salt',
    'cranberry-grape juice drink': 'juice drink, cranberry-grape',
    'cranberry-orange relish': 'cranberry orange relish',
    'cream cheese, plain, low-fat': 'cheese, cream, low-fat',
    'cream cooking': 'cooking cream',
    'cream sour': 'sour cream',
    'cream whipped, sweetened': 'whipped cream, sweetened',
    'cream whipping': 'cream, whipping',
    'cream, double': 'double cream',
    'cream, heavy': 'heavy cream',
    'cream, sour': 'sour cream',
    'cream, sour, fat-free': 'sour cream, fat-free',
    'cream, sour, low-fat': 'sour cream, low-fat',
    'cream, substitute': 'cream substitute',
    'creamy dressing, made with sour cream and/or buttermilk and oil, fat-free': 'creamy dressing, made with sour cream and/or buttermilk and oil, reduced calorie, fat-free',
    'cress garden': 'garden cress',
    'cress, garden': 'garden cress',
    'cress, garden, salted': 'garden cress, salted',
    'crispbread wholemeal': 'crispbread, wholemeal',
    'croissant butter': 'croissant, butter',
    'croissant cheese': 'croissant, cheese',
    'croissant with butter': 'croissant, butter',
    'crumble, apple': 'apple crumble',
    'curry, cauliflower and potato': 'potato cauliflower curry',
    'curry, lamb': 'lamb curry',
    'curry, potato and pea': 'curry, pea and potato',
    'curry, powder': 'curry powder',
    'custard with egg, refrigerated': 'egg custard, refrigerated',
    'custard-apple': 'custard apple',
    'dandelion greens, salted': 'dandelion, leaf, salted',
    'dandelion leaf': 'dandelion, leaf',
    'dark chocolate': 'chocolate, dark',
    'dark chocolate bar': 'dark chocolate, bar',
    'dessert soy, enriched': 'soy dessert, enriched',
    'dessert wine': 'wine, dessert',
    'dessert, fudgesicle bar, low-fat': 'dessert, frozen, fudge bar, low-fat',
    'dessert, ice cream': 'ice cream dessert',
    'double cheeseburger': 'cheeseburger, double',
    'double hamburger': 'hamburger, double',
    'doughnut, iced': 'doughnut iced',
    'doughnut, plain': 'doughnut plain',
    'dressing, french': 'french dressing',
    'dressing, french, fat-free': 'french dressing, fat-free',
    'drink almond': 'almond drink',
    'drink almond unsweetened': 'almond drink, unsweetened',
    'drink almond, enriched': 'almond drink, enriched',
    'drink coconut': 'coconut drink',
    'drink coconut without sugar': 'coconut drink without sugar',
    'drink oat': 'oat drink',
    'drink oat, enriched': 'oat drink, enriched',
    'drink soy': 'soy drink',
    'drink, grape': 'grape drink',
    'drinking chocolate, dried, low-fat': 'drinking chocolate powder, low-fat',
    'drumstick, leaf': 'drumstick leaf',
    'drumstick, leaf, salted': 'drumstick leaf, salted',
    'drumstick, pod': 'drumstick pod',
    'drumstick, pod, salted': 'drumstick pod, salted',
    'duck egg': 'egg, duck',
    'eel in homemade': 'eel homemade',
    'eel, smoked': 'eel smoked',
    'egg custard': 'custard, egg',
    'egg noodle, dried': 'noodle, egg, dried',
    'egg noodle, salted': 'noodle, egg, salted',
    'egg, white': 'egg white',
    'egg, white, dried': 'egg white, dried',
    'egg, without shell': 'egg without shell',
    'egg, yolk': 'egg yolk',
    'egg, yolk, dried': 'egg yolk, dried',
    'egg, yolk, salted': 'egg yolk, salted',
    'energy drink red bull': 'energy drink, red bull',
    'english muffin': 'muffin, english',
    'english muffin, enriched': 'muffin, english, enriched',
    'european chestnut': 'chestnut, european',
    'european chestnut, dried': 'chestnut, european, dried',
    'fat chicken': 'chicken fat',
    'fat, chicken': 'chicken fat',
    'fat, chicken, skin': 'chicken fat, skin',
    'fat, goose': 'goose fat',
    'fat, turkey': 'turkey fat',
    'fennel, florence': 'florence fennel',
    'fenugreek, leaf': 'fenugreek leaf',
    'fermented milk drink': 'milk drink, fermented',
    'filo pastry': 'pastry, filo',
    'fish fillet': 'fish, fillet',
    'fish sandwich with tartar sauce and cheese': 'fish sandwich, with tartar sauce and cheese',
    'flour almond': 'almond flour',
    'flour buckwheat': 'buckwheat flour',
    'flour cassava': 'cassava flour',
    'flour coconut': 'coconut flour',
    'flour corn': 'corn flour',
    'flour rice': 'rice flour',
    'flour rye': 'rye flour',
    'flour soy': 'soy flour',
    'flour spelt': 'spelt flour',
    'flour spelt wholemeal': 'spelt, wholemeal flour',
    'flour, all-purpose': 'flour, all purpose',
    'flour, all-purpose, enriched': 'flour, all purpose, enriched',
    'flour, whole wheat': 'whole wheat flour',
    'flour, wholegrain oat': 'oat flour, wholegrain',
    'foxtail millet flour': 'millet flour, foxtail',
    'french bean': 'bean, french',
    'french bread, country-style': 'country-style bread, french',
    'french fries, fried in vegetable oil': 'french fry, fried',
    'frog, leg': 'frog leg',
    'frog, leg, fried': 'frog leg, fried',
    'fruit citrus': 'citrus fruit',
    'fruit drink concentrate undiluted': 'fruit drink concentrate, undiluted',
    'fruit drink, hyvää päivää, sweetened, enriched': 'fruit drink, citrus, partially artificially sweetened, enriched',
    'fruit juice drink': 'juice drink, fruit',
    'fruit juice, concentrated': 'fruit juice concentrated',
    'fruit juice, concentrated, enriched': 'fruit juice concentrated, enriched',
    'fruit juice, mixed': 'mixed fruit juice',
    'fruit mixed, dried': 'mixed fruit, dried',
    'fruit punch drink': 'drink, fruit punch',
    'fruit yoghurt': 'yoghurt, fruit',
    'fruit yoghurt, enriched': 'yoghurt, fruit, enriched',
    'fufu of and plantain with red palm oil': 'fufu of plantain with red palm oil',
    'garlic bread': 'bread, garlic',
    'gelatin dessert, mix': 'dessert, gelatin mix',
    'ginger, root': 'ginger root',
    'gingerbread cake': 'cake, gingerbread',
    'globe artichoke': 'artichoke, globe',
    'goat meat': 'goat, meat',
    'goat meat, lean': 'goat, meat, lean',
    'goose egg': 'egg, goose',
    'gourd, ash': 'ash gourd',
    'gourd, bottle': 'bottle gourd',
    'gourd, ridge': 'ridge gourd',
    'gourd, white-flowered': 'gourd, white flowered',
    'gourd, white-flowered, salted': 'gourd, white flowered, salted',
    'grahams cracker': 'cracker, graham',
    'grain, amaranth': 'amaranth grain',
    'grain, barley': 'barley, grain',
    'grain, barley flour': 'grain, barley, flour',
    'grain, buckwheat': 'buckwheat, grain',
    'grain, rice flour': 'grain, rice, flour',
    'grain, rye': 'rye, grain',
    'grain, sorghum': 'sorghum, grain',
    'grain, triticale': 'triticale, grain',
    'grain, wheat': 'wheat grain',
    'grapefruit juice white': 'grapefruit juice, white',
    'grapefruit red': 'grapefruit, red',
    'great northern bean': 'bean, great northern',
    'greek yoghurt-style': 'greek yoghurt style',
    'green banana': 'banana, green',
    'green cabbage': 'cabbage, green',
    'green pea': 'pea, green',
    'green peas product, salted': 'green pea, frozen, salted',
    'green tea': 'tea, green',
    'green tea, brewed': 'tea, brewed, green',
    'ground turkey': 'turkey, ground',
    'guinea fowl meat': 'guinea fowl, meat',
    'halibut, greenland': 'greenland halibut',
    'halibut, smoked': 'halibut smoked',
    'ham lean': 'ham, lean',
    'ham smoked': 'ham, smoked',
    'hamburger with cheese': 'hamburger, with cheese',
    'hard biscuit': 'biscuit, hard',
    'hare meat': 'hare, meat',
    'hazelnut-chocolate spread': 'chocolate hazelnut spread',
    'herbal tea': 'tea, herbal',
    'hollandaise sauce, homemade': 'hollandaise sauce homemade',
    'horse meat': 'horse, meat',
    'hot chocolate, cocoa': 'cocoa, hot chocolate',
    'hp-sauce, brown sauce': 'hp sauce brown sauce',
    'hulled barley': 'barley, hulled',
    'hyacinth-bean': 'hyacinth bean',
    'hyacinth-bean, salted': 'hyacinth bean, salted',
    'hyacinth-bean, unsalted': 'hyacinth bean, unsalted',
    'ice cream chocolate': 'chocolate ice cream',
    'ice cream dairy': 'ice cream, dairy',
    'ice cream, chocolate': 'chocolate ice cream',
    'ice cream, cone': 'ice cream cone',
    'ice cream, soft-serve': 'ice cream, soft serve',
    'ice cream, stick, chocolate': 'ice cream stick',
    'ice cream, strawberry': 'strawberry ice cream',
    'ice cream, vanilla': 'vanilla ice cream',
    'ice cream, vanilla flavour': 'ice cream vanilla flavour',
    'ice cream, vanilla, fat-free': 'vanilla ice cream, fat-free',
    'icing sugar': 'sugar, icing',
    'instant coffee': 'coffee, instant',
    'instant coffee powder': 'instant coffee, powder',
    'jam, blackberry': 'blackberry jam',
    'jam, blueberry': 'blueberry jam',
    'jam, cherry': 'cherry jam',
    'jam, fruit': 'fruit jam',
    'jam, raspberry': 'raspberry jam',
    'jam, raspberry, low-sugar': 'raspberry jam, low-sugar',
    'jam, strawberry': 'strawberry jam',
    'jam, strawberry, low-sugar': 'strawberry jam, low-sugar',
    'japanese chestnut': 'chestnut, japanese',
    'japanese chestnut, dried': 'chestnut, japanese, dried',
    'jelly sweet': 'sweet, jelly',
    'jerusalem-artichoke': 'jerusalem artichoke',
    'juice apple': 'apple juice',
    'juice carrot': 'carrot juice',
    'juice drink, cranberry-apple low calorie': 'cranberry-apple juice drink, low calorie',
    'juice drink, grape': 'grape juice drink',
    'juice drink, orange': 'orange juice drink',
    'juice drink, orange and apricot': 'orange and apricot juice drink',
    'juice drink, pineapple and grapefruit': 'pineapple and grapefruit juice drink',
    'juice drink, pineapple and orange': 'pineapple and orange juice drink',
    'juice grape': 'grape juice',
    'juice grapefruit': 'grapefruit juice',
    'juice lemon': 'lemon juice',
    'juice orange': 'orange juice',
    'juice orange freshly squeezed': 'orange juice, freshly squeezed',
    'juice pear': 'pear juice',
    'juice pineapple': 'pineapple juice',
    'kale curly': 'curly kale',
    'kenaf leaf, ife-ken': 'kenaf, ife-ken, leaf',
    'ketchup tomato': 'tomato ketchup',
    'kidney bean': 'bean, kidney',
    'kidney bean, red': 'bean, kidney, red',
    'kidney lamb': 'lamb, kidney',
    'kidney ox': 'ox, kidney',
    'kidney pork': 'pork, kidney',
    'kiwi fruit green': 'kiwi fruit, green',
    'kiwi fruit yellow': 'kiwi fruit, yellow',
    'kola, bitter': 'bitter kola',
    'ladyfinger cookie': 'cookie, ladyfinger',
    'lamb meat': 'lamb, meat',
    'lamb mince': 'lamb, mince',
    'lamb mince, fried': 'lamb, mince, fried',
    'lamb, chop': 'lamb chop',
    'lamb, chop, fried': 'lamb chop, fried',
    'lamb, lean, breast': 'lamb, breast, lean',
    'lamb, lean, leg': 'lamb, leg, lean',
    'lamb, lean, loin': 'lamb, loin, lean',
    'lamb, lean, neck': 'lamb, neck, lean',
    'lamb, lean, shoulder': 'lamb, shoulder, lean',
    'lamb, ragout': 'lamb ragout',
    'lamb, steak': 'lamb steak',
    'lard, pork': 'pork lard',
    'lasagna with meat and sauce': 'lasagna with meat sauce',
    'lasagna with meat and sauce, low-fat': 'lasagna with meat & sauce, low-fat',
    'lasagna with vegetable': 'lasagna, vegetable',
    'lasagne, pork-beef mince': 'lasagne, pork beef mince',
    'lasagne, vegetable': 'vegetable lasagne',
    'lecithin, soy': 'soy lecithin',
    'lemon juice from concentrate, or bottled': 'lemon juice from concentrate, bottled',
    'lentil soup': 'soup, lentil',
    'lettuce butterhead': 'butterhead lettuce',
    'lettuce iceberg': 'iceberg lettuce',
    'lettuce red': 'lettuce, red',
    'lettuce romaine': 'romaine lettuce',
    'limburger cheese': 'cheese, limburger',
    'liqueur, coffee with cream': 'liqueur, coffee and cream',
    'liqueur, cream': 'cream liqueur',
    'liquorice, allsort': 'liquorice allsort',
    'liver chicken': 'chicken, liver',
    'liver ox': 'ox, liver',
    'liver pate': 'pate, liver',
    'liver pork': 'pork, liver',
    'liver pâté, pork': 'pork liver pâté',
    'liver, reindeer': 'reindeer liver',
    'low-fat margarine': 'margarine, low-fat',
    'macadamia': 'macadamia nut',
    'macadamia, salted': 'macadamia nut, salted',
    'macaroni & cheese, homemade': 'macaroni cheese, homemade',
    'maize, refined flour with vitamin, enriched': 'corn flour, white, enriched',
    'malabar spinach, leaf': 'malabar spinach leaf',
    'marble cake': 'cake marble',
    'margarine light, low-fat': 'margarine, light, low-fat',
    'margarine vegetable-oil, 60% fat, enriched': 'margarine, vegetable oil, 60% fat, enriched',
    'margarine vegetable-oil, enriched': 'margarine, vegetable oil, enriched',
    'margarine, low-fat, enriched': 'margarine, enriched, low-fat',
    'margarine, margarine-butter blend': 'margarine, butter-margarine blend',
    'mashed potato': 'potato, mashed',
    'mashed potato, dried': 'potato, mashed, dried',
    'meat loaf in homemade': 'meat loaf, homemade',
    'meat, with herb, red': 'red meat',
    'meatless, bacon': 'bacon, meatless',
    'meatless, bacon bit': 'bacon bit, meatless',
    'meatless, chicken': 'chicken, meatless',
    'meatless, chicken, fried': 'chicken, meatless, fried',
    'meatless, luncheon slice': 'luncheon slice, meatless',
    'meatless, meatball': 'meatball, meatless',
    'melon honeydew': 'honeydew melon',
    'meringue with chocolate': 'chocolate meringue',
    'milk chocolate bar': 'milk chocolate, bar',
    'milk chocolate with hazelnut': 'chocolate milk with hazelnut',
    'milk condensed with sugar': 'condensed milk, with sugar',
    'milk condensed, sweetened': 'milk, condensed, sweetened',
    'milk goat, full-fat': 'goat milk, full-fat',
    'milk human': 'human milk',
    'milk powder, fat-free': 'milk, dried, fat-free',
    'milk whole': 'milk, whole',
    'milk whole, dried': 'milk, whole, dried',
    'milk, enriched, low-fat': 'milk, low-fat, enriched',
    'millet, grain': 'grain, millet',
    'mince chicken': 'chicken, mince',
    'minestrone soup': 'soup, minestrone',
    'mint leaf': 'mint, leaf',
    'miso soy paste': 'miso, soy paste',
    'mix of nut, salted': 'nut mix, salted',
    'mixed': 'mixed nut',
    'mixed berries': 'mixed berry',
    'mixed vegetable': 'vegetable, mixed',
    'mixed vegetable, salted': 'vegetable, mixed, salted',
    'mixed vegetable, unsalted': 'vegetable, mixed, unsalted',
    'monterey jack cheese': 'cheese, monterey jack',
    'moussaka, pork-beef mince': 'moussaka, pork beef mince',
    'moussaka, pork-beef mince, low-fat': 'moussaka, pork beef mince, low-fat',
    'mousse chocolate': 'chocolate mousse',
    'mousse, chocolate': 'chocolate mousse',
    'mousse, chocolate, low-fat': 'chocolate mousse, low-fat',
    'muesli biscuit': 'biscuit muesli',
    'muesli, cereal': 'cereal, muesli',
    'muesli, fruit': 'muesli with fruit',
    'muesli, with fruit': 'muesli with fruit',
    'muesli, with fruit, enriched': 'muesli with fruit, enriched',
    'muffin, american': 'american muffin',
    'multigrain bread brown with seed, fortified with iron and vitamins': 'multigrain bread, brown, with seeds, fortified with iron and vitamins',
    'multigrain bread, carrot roll, flour, fat-free': 'multigrain bread, carrot roll, fat-free',
    'multigrain bread, graham flour': 'multigrain bread, graham, flour',
    'mung bean, salted': 'mung bean, salted, sprouted',
    'mungo bean': 'bean, mungo',
    'mungo bean, salted': 'bean, mungo, salted',
    'mushroom soup': 'soup, mushroom',
    'mushroom, king oyster': 'king oyster mushroom',
    'mushroom, oyster': 'oyster mushroom',
    'mushroom, shiitake': 'shiitake mushroom',
    'mushroom, shiitake, dried': 'shiitake mushroom, dried',
    'mushroom, shiitake, salted': 'shiitake mushroom, salted',
    'mushroom-barley soup': 'soup, mushroom barley',
    'mustard, leaf': 'mustard leaf',
    'mustard, powder': 'mustard powder',
    'navy bean': 'bean, navy',
    'navy bean, dried': 'bean, navy, dried',
    'navy bean, dried, 0% moisture basis': 'bean, navy, dried, 0% moisture basis',
    'nectar, mango': 'mango nectar',
    'new potato': 'potato, new',
    'new potato, salted': 'potato, new, salted',
    'new potato, unsalted': 'potato, new, unsalted',
    'non-alcoholic, wine': 'wine, non-alcoholic',
    'novelties creamsicle pop': 'novelty, creamsicle pop',
    'novelty, klondike, fat-free': 'novelty, ice cream type, fat-free',
    'nut and dried fruit mix, salted': 'nut and fruit mix, dried, salted',
    'nut, brazil, unsalted': 'brazil nut, unsalted',
    'nutmeg ground': 'nutmeg, ground',
    'oat biscuit': 'biscuit, oat',
    'oat flour wholegrain': 'oat flour, wholegrain',
    'oat, porridge': 'oat porridge',
    'oat, wholegrain flour': 'oat flour, wholegrain',
    'oil coconut': 'coconut oil',
    'oil corn': 'corn oil',
    'oil flaxseed': 'flaxseed oil',
    'oil olive': 'olive oil',
    'oil palm': 'palm oil',
    'oil peanut': 'peanut oil',
    'oil rapeseed': 'rapeseed oil',
    'oil rice bran': 'rice bran oil',
    'oil safflower': 'safflower oil',
    'oil sesame': 'sesame oil',
    'oil sunflower seed': 'sunflower oil',
    'oil vegetable': 'oil, vegetable',
    'olive, green manzanilla': 'olive, manzanilla, green',
    'olive, green stuffed with pimento': 'olive, stuffed with pimento, green',
    'omelette with ham and cheese': 'ham and cheese omelette',
    'omelette, cheese': 'cheese omelette',
    'onion red': 'red onion',
    'onion soup': 'soup, onion',
    'onion soup, french': 'soup, french onion',
    'onion welsh': 'onion, welsh',
    'onion yellow': 'onion, yellow',
    'onion, pickled': 'pickled onion',
    'onion, red': 'red onion',
    'onion, without added fat, fried, red': 'red onion, fried',
    'or bottled juice, sweetened': 'juice, canned or bottled, sweetened',
    'orange juice freshly-squeezed': 'orange juice, freshly squeezed',
    'orange marmalade': 'marmalade, orange',
    'ovaltine powder': 'ovaltine, powder',
    'ovaltine powder, enriched': 'malted chocolate drink powder, fortified',
    'oxheart cabbage': 'cabbage oxheart',
    'pancake plain': 'pancake',
    'paprika, powder': 'paprika powder',
    'parsley, leaf': 'leaf parsley',
    'passion-fruit juice, purple': 'passion fruit juice, purple',
    'passion-fruit juice, yellow': 'passion fruit juice, yellow',
    'passion-fruit, purple': 'passion fruit, purple',
    'pastry cream': 'cream pastry',
    'pastry, danish': 'danish pastry',
    'pea soup': 'soup, pea',
    'pea sprout': 'pea, sprout',
    'pea, carrot, salted': 'peas and carrot, salted',
    'pea, pigeon': 'pigeon pea',
    'pea, pigeon, salted': 'pigeon pea, salted',
    'pea, split': 'split pea',
    'pea, split green': 'split pea, green',
    'pea, split, dried': 'split pea, dried',
    'pea, split, green': 'split pea, green',
    'pea, split, salted': 'split pea, salted',
    'pea, sugar-snap': 'sugar snap pea',
    'peach nectarine': 'peach/nectarine',
    'peanut, chocolate coated': 'chocolate-coated peanut',
    'peanut, peanut butter': 'peanut butter with peanut',
    'pearl barley': 'barley, pearl',
    'pearl onions pickled': 'pearl onion, pickled',
    'pecan': 'pecan nut',
    'pecan, salted': 'pecan nut, salted',
    'pennywort, indian': 'indian pennywort',
    'pepper sauce': 'sauce, pepper',
    'pepper spice, black': 'pepper, black',
    'pepper spice, green': 'pepper, green',
    'pepper spice, white': 'pepper, white',
    'pepper, bell': 'bell pepper',
    'pepper, bell, green': 'bell pepper, green',
    'pepper, bell, red': 'bell pepper, red',
    'pepper, chili, green': 'chili pepper',
    'pepper, chilli': 'chilli pepper',
    'persimmon, japanese, dried': 'japanese persimmon, dried',
    'pesto, green': 'pesto green',
    'pesto, red': 'pesto red',
    'pie crust, standard-type': 'pie crust, standard type',
    'pie crust, standard-type, enriched': 'pie crust, standard type, enriched',
    'pie, apple, flour, enriched': 'pie, apple, enriched',
    'pie, cottage': 'cottage pie',
    'pie, fruit': 'fruit pie',
    'pie, steak and kidney': 'pie, steak & kidney',
    'pie, vegetable': 'vegetable pie',
    'pigeon-pea': 'pigeon pea',
    'pike-perch, wild': 'pike perch, wild',
    'pinto bean, dried': 'bean, pinto, dried',
    'pistachio nut, unsalted': 'nut, pistachio, unsalted',
    'pizza margherita': 'margherita pizza',
    'pizza salami': 'salami pizza',
    'pizza with ham': 'ham pizza',
    'pizza with ham and mushroom': 'ham and mushroom pizza',
    'pizza with ham and pineapple': 'ham and pineapple pizza',
    'pizza with tomato and cheese': 'pizza, cheese and tomato',
    'pizza with vegetable': 'pizza, vegetable',
    'pizza, cheese & tomato': 'pizza, cheese and tomato',
    'pizza, cheese and vegetable': 'pizza with cheese and vegetable',
    'pizza, ham & pineapple': 'ham and pineapple pizza',
    'pizza, ham and cheese': 'ham and cheese pizza',
    'pizza, ham and pineapple': 'ham and pineapple pizza',
    'pizza, kebab': 'kebab pizza',
    'pizza, seafood': 'seafood pizza',
    'pizza, tuna': 'pizza with tuna',
    'pizza, vegetarian': 'vegetarian pizza',
    'plain cake': 'cake',
    'plant-based bits soy protein': 'plant-based soy protein bits',
    'plant-based bits soy protein, fried': 'plant-based soy protein bits, fried',
    'plantain cooking banana': 'plantain, plantain',
    'plus milk, fat-free': 'milk, with added vitamin D, fat-free',
    'polar bear, flesh': 'bear, polar, flesh',
    'popcorn, microwave': 'microwave popcorn',
    'popcorn, microwave, low-fat': 'microwave popcorn, low-fat',
    'pork sausage': 'sausage, pork',
    'pork sausage, low-salt': 'sausage, pork, low-salt',
    'porridge oat': 'oat porridge',
    'porridge, tô, flour, enriched': 'porridge, tô, from corn and cassava flour, enriched',
    'potato chip': 'potato, chip',
    'potato chip, low-fat': 'potato, chip, low-fat',
    'potato crisp, fried': 'crisp, potato, fried',
    'potato crisps flavoured': 'crisps potato flavoured',
    'potato old': 'potato, old',
    'potato old, salted': 'potato, old, salted',
    'potato old, unsalted': 'potato, old, unsalted',
    'potato sweet': 'sweet potato',
    'potato, asterix': 'potato asterix',
    'potato, crisp, fried': 'crisp, potato, fried',
    'potato, croquette': 'potato croquette',
    'potato, duchesse': 'duchesse potato',
    'potato, hashed-brown': 'potato, hashed brown',
    'potato, in skin, salted': 'potato, skin, salted',
    'potato, in skin, unsalted': 'potato, skin, unsalted',
    'potato, old, salted, flesh only': 'potato, old, salted',
    'potato, tuber': 'potato tuber',
    'potato, vacuum-packed': 'potato, vacuum packed',
    'potato, white flesh and skin': 'potato, flesh and skin, white',
    'potato, white skin and flesh': 'potato, flesh and skin, white',
    'pound cake': 'cake, pound',
    'prawn, shrimp, wild': 'shrimp or prawn',
    'process cheese spread': 'cheese spread, process',
    'protein powder whey based': 'protein powder, whey based',
    'pudding, caramel': 'caramel pudding',
    'pudding, rice': 'rice pudding',
    'puff pastry': 'pastry, puff',
    'pumpkin and squash, unsalted': 'pumpkin, squash, unsalted',
    'pumpkin leaf': 'pumpkin, leaf',
    'pumpkin leaf, salted': 'pumpkin, leaf, salted',
    'pumpkin leaf, unsalted': 'pumpkin, leaf, unsalted',
    'pumpkin, flower, salted': 'pumpkin flower, salted',
    'pumpkin, squash': 'pumpkin and squash',
    'quail egg': 'egg, quail',
    'quark flavoured, sweetened': 'quark, flavoured, sweetened',
    'quark plain': 'quark',
    'quark with fruit': 'fruit quark',
    'quark, plain, sweetened': 'quark plain, sweetened',
    'quesadilla with cheese': 'cheese quesadilla',
    'quesadilla, with chicken': 'quesadilla, chicken',
    'quiche, lorraine': 'quiche lorraine',
    'rabbit domesticated': 'rabbit, domesticated',
    'rabbit wild': 'rabbit, wild',
    'radish black': 'radish, black',
    'radish, leaf': 'radish leaf',
    'rapeseed oil, cold-pressed': 'rapeseed oil cold pressed',
    'raspberries product, sweetened': 'raspberry, frozen, sweetened',
    'ravioli, cheese with tomato sauce': 'ravioli, cheese and tomato sauce',
    'red beet': 'beet, red',
    'red kidney bean': 'bean, kidney, red',
    'red kidney bean, dried': 'bean, kidney, red, dried',
    'red raspberry': 'raspberry, red',
    'red wine': 'wine, red',
    'rice brown': 'rice, brown',
    'rice flour white': 'rice flour, white',
    'rice noodle, dried': 'noodle, rice, dried',
    'rice, bran': 'rice bran',
    'rice, brown, flour': 'rice flour, brown',
    'ringed seal, brain': 'seal, ringed, brain',
    'ringed seal, flesh': 'seal, ringed, flesh',
    'ringed seal, heart': 'seal, ringed, heart',
    'ringed seal, liver': 'seal, ringed, liver',
    'risotto, chicken': 'chicken risotto',
    'risotto, vegetable': 'vegetable risotto',
    'risotto, with vegetable': 'vegetable risotto',
    'ritz cracker': 'cracker, ritz',
    'roe, cod': 'cod, roe',
    'root parsley': 'parsley root',
    'rusk, wholemeal': 'rusk wholemeal',
    'rutabaga or rutabaga': 'rutabaga, rutabaga',
    'rye and wheat bread': 'bread wheat rye',
    'rye crispbread': 'crispbread, rye',
    'rye flakes wholegrain': 'rye, wholegrain flake',
    'rye flour wholegrain': 'rye flour, wholegrain',
    'rye grain': 'rye, grain',
    'rye, wholegrain flour': 'rye flour, wholegrain',
    'rye-flour, bolted': 'rye flour, bolted',
    'saccharin sweetener': 'sweetener, saccharin',
    'saithe, pollock': 'pollock, saithe',
    'salad coleslaw': 'coleslaw salad',
    'salad dressing honey/mustard': 'salad dressing, honey mustard',
    'salad dressing with mayonnaise': 'salad dressing, mayonnaise',
    'salad dressing, blue or roquefort cheese dressing, fat-free': 'salad dressing, blue or roquefort cheese, fat-free',
    'salad dressing, honey/mustard': 'salad dressing, honey mustard',
    'salad dressing, mayonnaise and mayonnaise type': 'salad dressing, mayonnaise and mayonnaise-type',
    'salad dressing, mayonnaise and mayonnaise type, fat-free': 'salad dressing, mayonnaise and mayonnaise-type, fat-free',
    'salad dressing, mayonnaise-type': 'salad dressing, mayonnaise type',
    'salad egg with potato': 'potato salad with egg',
    'salad potato': 'potato salad',
    'salad russian': 'russian salad',
    'salad tuna': 'tuna salad',
    'salad, beetroot': 'beetroot salad',
    'salad, caesar': 'caesar salad',
    'salad, caesar with chicken': 'caesar salad with chicken',
    'salad, coleslaw': 'coleslaw salad',
    'salad, potato': 'potato salad',
    'salad, red cabbage and apple': 'red cabbage salad with apple',
    'salad, rice': 'rice salad',
    'salad, taco': 'taco salad',
    'salad, waldorf': 'waldorf salad',
    'salmon smoked': 'salmon, smoked',
    'salsify, black': 'black salsify',
    'salt sea': 'sea salt',
    'salty snack, cracker': 'cracker, salty snack',
    'salty snack, cracker, low-fat': 'cracker, salty snack, low-fat',
    'sandwich baguette': 'sandwich, baguette',
    'sandwich spread, meatless': 'meatless, sandwich spread',
    'sandwich with ham': 'sandwich, ham',
    'sandwich, ham and cheese': 'sandwich with ham and cheese',
    'sardines in tomato sauce': 'sardine, in tomato sauce',
    'sauce barbecue': 'sauce, barbecue',
    'sauce chilli': 'chilli sauce',
    'sauce garlic': 'garlic sauce',
    'sauce oyster': 'sauce, oyster',
    'sauce, barbeque': 'barbeque sauce',
    'sauce, cranberry, sweetened': 'cranberry sauce, sweetened',
    'sauce, hollandaise': 'hollandaise sauce',
    'sauce, horseradish': 'horseradish sauce',
    'sauce, peanut': 'peanut sauce',
    'sauce, soy': 'soy sauce',
    'sauce, soy, low-salt': 'soy sauce, low-salt',
    'sauce, tomato': 'tomato sauce',
    'sauce, tomato, low-salt': 'tomato sauce, low-salt',
    'sausage beef': 'sausage, beef',
    'sausage chorizo': 'sausage, chorizo',
    'sausage smoked': 'smoked sausage',
    'sausage, bratwurst': 'sausage bratwurst',
    'sausage, bratwurst, low-fat': 'sausage bratwurst, low-fat',
    'sausage, chicken': 'chicken sausage',
    'sausage, meatless': 'meatless, sausage',
    'sausage, merguez': 'merguez sausage',
    'sausage, salami': 'sausage salami',
    'sausage, smoked': 'smoked sausage',
    'sausage, smoked, low-fat': 'smoked sausage, low-fat',
    'sausage, turkey and pork': 'sausage, pork and turkey',
    'savoury biscuit': 'biscuit, savoury',
    'savoy cabbage': 'cabbage, savoy',
    'scrambled egg': 'egg, scrambled',
    'seaweed kelp': 'seaweed, kelp',
    'seaweed nori, dried': 'seaweed, nori, dried',
    'seeds and peanut, dried, salted': 'peanut, dried, salted',
    'semi-skimmed milk': 'milk, low-fat',
    'semi-skimmed milk 1.5% fat fortified': 'milk, low-fat, enriched',
    'semolina porridge': 'porridge semolina',
    'semolina, wheat': 'semolina wheat',
    'sesame seed paste tahini, salted': 'tahini sesame paste',
    'sesame seeds whole': 'sesame, whole',
    'shake, vanilla': 'vanilla shake',
    "shepherd's-purse, leaf": "shepherd's purse, leaf",
    'skimmed milk': 'milk, fat-free',
    'smoothie fruit': 'fruit smoothie',
    'snack, potato chip, salted': 'potato, chip, salted',
    'soft caramel candy': 'caramel soft candy',
    'sorghum grain': 'sorghum, grain',
    'souffle, cheese': 'cheese souffle',
    'soup clear with vegetable': 'vegetable soup clear',
    'soup vegetable': 'soup, vegetable',
    'soup, bean & ham': 'soup, bean with ham',
    'soup, bean & ham, low-salt': 'soup, bean with ham, low-salt',
    'soup, beef and vegetable': 'soup, vegetable beef',
    'soup, beef and vegetable, low-salt': 'soup, vegetable beef, low-salt',
    'soup, beef mushroom': 'soup, beef and mushroom',
    'soup, carrot': 'carrot soup',
    'soup, chicken & noodle, dried': 'soup, chicken noodle, dried',
    'soup, chicken and vegetable': 'soup, chicken vegetable',
    'soup, chicken and vegetable, dried': 'soup, chicken vegetable, dried',
    'soup, chicken with rice': 'soup, chicken rice',
    'soup, cream of chicken': 'chicken soup, cream',
    'soup, cream of chicken, low-salt': 'chicken soup, cream, low-salt',
    'soup, cream of vegetable': 'vegetable soup with cream',
    'soup, goulash': 'goulash soup',
    'soup, potato and leek': 'soup, leek and potato',
    'soup, pumpkin': 'pumpkin soup',
    'soup, spinach': 'spinach soup',
    'soup, tomato noodle': 'soup tomato with noodle',
    'soup, vegetable chicken': 'soup, chicken vegetable',
    'soup, vegetable chicken, low-salt': 'soup, chicken vegetable, low-salt',
    'soup, vegetable with beef': 'soup, vegetable beef',
    'soy dessert, flavoured with calcium and vitamin d, sweetened, enriched': 'soy dessert, flavoured, with sugar, enriched',
    'soy drink, flavoured in calcium, sweetened, enriched': 'soy drink, flavoured, with sugar, enriched',
    'soy paste, miso': 'miso, soy paste',
    'soy sauce shoyu': 'soy sauce, shoyu',
    'soy, pudding': 'soy pudding',
    'soy-textured vegetable protein': 'textured vegetable soy protein',
    'soya bean': 'bean, soy',
    'soybean, sprout': 'soybean sprout',
    'soymilk, chocolate': 'chocolate, soymilk',
    'soymilk, chocolate and d, enriched': 'chocolate, soymilk, enriched',
    'spelt flour wholegrain': 'spelt flour, wholegrain',
    'spelt flour, wholemeal': 'spelt, wholemeal flour',
    'sponge cake': 'cake, sponge',
    'spread chocolate hazelnut': 'chocolate hazelnut spread',
    'sprout, alfalfa': 'alfalfa sprout',
    'sprout, bean': 'bean sprout',
    'starch sweetener, glucose fructose syrup': 'starch sweetener, fructose glucose syrup',
    'starch, potato': 'potato starch',
    'strawberries product, sweetened': 'strawberry, frozen, sweetened',
    'sugar powdered': 'sugar, dried',
    'sundae, hot fudge': 'hot fudge sundae',
    'sundae, strawberry': 'strawberry sundae',
    'super-skimmed milk <0.1% fat fortified': 'milk, enriched, fat-free',
    'sushi, salmon nigiri': 'sushi nigiri with salmon',
    'sushi, tuna nigiri': 'sushi nigiri with tuna',
    'sweet and sour sauce': 'sauce, sweet and sour',
    'sweet pepper': 'pepper, sweet',
    'sweet pepper green': 'pepper, sweet, green',
    'sweet pepper red': 'pepper, sweet, red',
    'sweet pepper yellow': 'pepper, sweet, yellow',
    'sweet pepper, green': 'pepper, sweet, green',
    'sweet pepper, red': 'pepper, sweet, red',
    'sweet pepper, unsalted': 'pepper, sweet, unsalted',
    'sweet pepper, yellow': 'pepper, sweet, yellow',
    'sweet potato, after baking, salted': 'sweet potato, baked, skin removed, salted',
    'sweet potato, leaf': 'sweet potato leaf',
    'sweet potato, leaf, unsalted': 'sweet potato leaf, unsalted',
    'sweet whey': 'whey, sweet',
    'sweet, marmalade': 'marmalade sweet',
    'sweetener, stevia': 'stevia sweetener',
    'swiss chard, leaf': 'swiss chard leaf',
    'syrup corn': 'syrup, corn',
    'syrup maple': 'maple syrup',
    'syrup, chocolate': 'chocolate syrup',
    'syrup, golden': 'golden syrup',
    'syrup, maple': 'maple syrup',
    'syrup, rice malt': 'syrup rice malt',
    'syrup, with sugar': 'syrup sugar',
    'tahini, paste': 'tahini paste',
    'tallow, beef fat': 'fat, beef tallow',
    'tamale, cheese': 'cheese tamale',
    'tamarind leaf': 'tamarind, leaf',
    'tap water': 'water, tap',
    'tapioca, pearl, dried': 'tapioca pearl, dried',
    'taro, leaf': 'taro leaf',
    'taro, leaf, salted': 'taro leaf, salted',
    'taro, leaf, unsalted': 'taro leaf, unsalted',
    'taro, shoot, salted': 'taro shoot, salted',
    'tart cherry juice': 'cherry juice, tart',
    'tartar sauce': 'sauce, tartar',
    'tea brewed': 'tea, brewed',
    'tea, leaf': 'tea leaf',
    'teriyaki sauce': 'sauce, teriyaki',
    'thousand island dressing': 'dressing, thousand island',
    'tigernut, tuber- flour': 'tigernut, tuber-flour',
    'tilapia fish': 'fish, tilapia',
    'toaster pastry, brown sugar and cinnamon': 'toaster pastry, brown-sugar-cinnamon',
    'tomato and mozzarella pizza': 'mozzarella and tomato pizza',
    'tomato sauce, with mushroom': 'sauce, tomato and mushroom',
    'tomato soup': 'soup, tomato',
    'tomato, paste': 'tomato paste',
    'tomatoes with oil, dried': 'tomato with oil, dried',
    'tongue calf': 'calf, tongue',
    'tongue sausage': 'sausage tongue',
    'tortilla, wheat': 'wheat tortilla',
    'trout, rainbow, fillet': 'rainbow trout fillet',
    'trout, steelhead, flesh, dried': 'trout, steelhead, dried, flesh',
    'tuna pizza': 'pizza with tuna',
    'tuna yellowfin': 'tuna, yellowfin',
    'turkey and gravy': 'gravy, turkey',
    'turkey egg': 'egg, turkey',
    'turkey vegetable gravy with cooking cream': 'turkey vegetable gravy, cooking cream',
    'turkey, crumble, ground, 7% fat': 'turkey, ground, crumble, 7% fat',
    'turkey, ground, fat-free': 'ground turkey, fat-free',
    'turkey, light, dark meat': 'turkey, light or dark meat',
    'turkey, patty, ground, 7% fat': 'turkey, ground, patty, 7% fat',
    'turnip, prairie': 'prairie turnip',
    'turnip, root': 'turnip root',
    'unsweetened, soymilk': 'soymilk, unsweetened',
    'vanilla pudding': 'pudding, vanilla',
    'vegetable au gratin': 'vegetable, au gratin',
    'vegetable oil': 'oil, vegetable',
    'vegetable oil blend': 'vegetable oil, blend',
    'vegetable soup': 'soup, vegetable',
    'vegetable soup, low-salt': 'soup, vegetable, low-salt',
    'vegetables': 'vegetable',
    'vegetarian sausage with soy and wheat protein heated': 'vegetarian sausage with soy and wheat protein, heated',
    'vienna sausage': 'sausage, vienna',
    'vinegar balsamic': 'balsamic vinegar',
    'waffle, gluten-free': 'waffle gluten-free',
    'wasabi paste': 'wasabi, paste',
    'wasabi root': 'wasabi, root',
    'water convolvulus, salted': 'water spinach, salted',
    'water, coconut': 'coconut water',
    'water, mineral': 'mineral water',
    'water, soda': 'soda water',
    'waterleaf leaf': 'waterleaf, leaf',
    'welsh onion': 'onion, welsh',
    'whale, ventral groove meat, flesh': 'whale, ventral groove, flesh',
    'wheat bread, with whole': 'whole wheat bread',
    'wheat bulgur': 'wheat, bulgur',
    'wheat cracker': 'cracker, wheat',
    'wheat flour durum': 'durum wheat flour',
    'wheat, bran': 'wheat bran',
    'wheat, durum': 'durum wheat',
    'wheat, germ': 'wheat germ',
    'wheat, khorasan': 'khorasan wheat',
    'wheat, puffed': 'puffed wheat',
    'wheat, puffed, enriched': 'puffed wheat, enriched',
    'wheat, semolina': 'semolina wheat',
    'wheat, with gluten': 'wheat, gluten',
    'whipped cream': 'cream, whipped',
    'whipping cream': 'cream, whipping',
    'white bean, dried': 'bean, white, dried',
    'white beans with tomato sauce': 'beans white in tomato sauce',
    'white cabbage': 'cabbage, white',
    'white pepper': 'pepper, white',
    'white radish': 'radish, white',
    'white roll, gluten-free': 'roll, gluten-free, white',
    'white sauce, savoury': 'sauce, savoury, white',
    'white wine': 'wine, white',
    'white wine, sweet': 'wine, sweet, white',
    'white yam': 'yam, white',
    'whole egg': 'egg, whole',
    'whole egg, dried': 'egg, whole, dried',
    'whole milk': 'milk, whole',
    'wholegrain barley flour': 'barley flour, wholegrain',
    'wholegrain mustard': 'mustard, wholegrain',
    'wholegrain rye bread': 'rye bread, wholegrain',
    'wholegrain rye flour': 'rye flour, wholegrain',
    'wiener, beef, fat-free': 'frankfurter, beef, fat-free',
    'wild boar': 'boar, wild',
    'wild boar meat': 'wild boar, meat',
    'wine red': 'wine, red',
    'wine rose': 'wine, rose',
    'wine white': 'wine, white',
    'wine white sweet': 'wine, sweet, white',
    'wine, mulled': 'mulled wine',
    'wine, port': 'port wine',
    'wine, port, enriched': 'port wine, enriched',
    'with soy protein': 'soy protein',
    'worcestershire sauce': 'sauce, worcestershire',
    'yeast extract marmite': 'yeast extract, marmite',
    "yeast, baker's": 'yeast, baker',
    "yeast, bakers'": 'yeast, baker',
    'yellow mealworm, larva': 'mealworm, larva, yellow',
    'yellow mustard': 'mustard, yellow',
    'yellow onion': 'onion, yellow',
    'yellow plum': 'plum, yellow',
    'yoghurt cream- with fruit': 'cream yoghurt, with fruit',
    'yoghurt drinkable flavoured 1% fat': 'yoghurt drinkable flavoured, 1% fat',
    'yoghurt greek, full-fat': 'greek yoghurt, full-fat',
    'yoghurt mild vanilla flavour': 'mild yoghurt vanilla flavour',
    'yoghurt mild vanilla flavour, enriched': 'mild yoghurt vanilla flavour, enriched',
    'yoghurt plain': 'yoghurt',
    'yoghurt vanilla, low-fat': 'yoghurt, vanilla, low-fat',
}


def _edit_distance_le(a: str, b: str, cap: int = 2) -> bool:
    """True if Levenshtein(a, b) <= cap. Small strings only; no numpy."""
    if abs(len(a) - len(b)) > cap:
        return False
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        if min(cur) > cap:
            return False
        prev = cur
    return prev[-1] <= cap


def _is_spelling_alias(kept: str, cand: str) -> bool:
    """True if `cand` is the same word as `kept` under a different spelling.

    FDC records both spellings of a food in one description - "Knackwurst,
    knockwurst, pork, beef" - and the second is not a qualifier, it is a repeat.
    Guarded hard: single words only, both at least 6 characters, and within two
    edits, so genuinely different short foods ("pea" / "tea") can never collide.
    """
    if " " in cand or len(cand) < 6 or len(kept) < 6:
        return False
    head = kept.rsplit(" ", 1)[-1] if " " in kept else kept
    if len(head) < 6:
        return False
    return _edit_distance_le(head, cand, 2)


def _final_polish(name: str) -> str:
    """Tidy an assembled canon name. Runs on every return path.

    Each rule in the core edits `d` at the point where it belongs, so a token
    removed late can leave punctuation that an earlier filter has already walked
    past: "Lobster, boiled/cooked in water" loses both prep words and comes out
    as 'lobster, /'. Re-ordering the core is not the fix - each step sits where
    it does for its own documented reason - so the finished name is tidied once,
    here, which is also the only place that can see a chunk repeated
    ('butter, butter', 'milk, milk', 'beef tenderloin, beef tenderloin').

    Cheap and idempotent: polishing an already-clean name is a no-op.
    """
    if not name:
        return name
    chunks: list[str] = []
    seen: set[str] = set()
    for c in name.split(","):
        c = re.sub(r"\s+", " ", c.strip(_POLISH_EDGE).strip())
        # A slash joins two names for one food and the sources space it three
        # ways: "baobab fruit /monkey bread" was a canon of its own beside
        # "baobab fruit/monkey bread", and 20 canons carried a spaced slash.
        c = re.sub(r"\s*/\s*", "/", c)
        c = _ANDOR_TRAIL_RE.sub("", c)
        c = _TRAIL_GRADE_RE.sub("", c)
        if len(c.split()) > 2:          # only when a real name remains
            c = _LINE_NUMBER_RE.sub("", c).strip()
        if (not c or _PUNCT_ONLY_RE.match(c) or _CONN_ONLY_RE.match(c)
                or _FILLER_CHUNK_RE.match(c) or _CORP_ENTITY_RE.match(c)
                or _NA_CHUNK_RE.match(c) or _STORAGE_CHUNK_RE.match(c)):
            continue
        # Never the head - a name has to keep an identity, and "plain" alone is
        # not one. The chocolate gate reads the head that is already standing.
        if _PULP_CHUNK_RE.match(c):
            c = "flesh"
        if (chunks and _NULL_QUALIFIER_RE.match(c)
                and not _PLAIN_IS_DARK_RE.match(chunks[0])):
            continue
        if c.lower() in seen:
            continue
        if any(_is_spelling_alias(k, c.lower()) for k in chunks):
            continue
        seen.add(c.lower())
        chunks.append(c)
    # A name made ENTIRELY of filler still has to keep an identity; falling
    # through to "" would silently drop the food from the index.
    if not chunks:
        return re.sub(r"\s+", " ", name.strip(_POLISH_EDGE).strip())
    # Doubled head noun. Fineli, CNF and STFCJ file foods as "class, item", and
    # the item usually repeats the class: "cake, cheese cake", "pudding, rice
    # pudding", "salad, potato salad", "sausage, liver sausage", "starch, corn
    # starch". The class adds nothing a reader needs, so the item stands alone.
    # Guarded on a whole-word suffix match, so "bean, green bean" collapses but
    # "bean, black" (a real variety, not a repeat) does not.
    if len(chunks) > 1 and chunks[1].lower().endswith(" " + chunks[0].lower()):
        chunks = chunks[1:]
    if len(chunks) > 1 and _MEAT_HEAD_RE.match(chunks[0]):
        chunks = [chunks[0]] + [c for c in chunks[1:] if not _ORIGIN_CHUNK_RE.match(c)]
    return ", ".join(chunks)


def canonicalize_food_name(desc: str) -> str:
    canon = _final_polish(_canonicalize_core(desc))
    # Degenerate result: every informative token was stripped and only a joining
    # word survived. "Seeds and kernels av" loses both nouns to _PREP_RE (which
    # carries "seeds" and "kernels" as preparation words) and canonicalized to
    # the bare conjunction "and". Fall back to the tidied raw description - a
    # long name beats a name that is not a food.
    if not canon or _CONN_ONLY_RE.match(canon) or _PUNCT_ONLY_RE.match(canon):
        # translate first: _clean_variant_name does not normalise full-width
        # punctuation, so the placeholder survives the emptiness test as "－".
        alt = _clean_variant_name(desc).translate(_FULLWIDTH_MAP)
        if not alt.strip(_POLISH_EDGE).strip():
            # The name field itself can be a placeholder. STFCJ files one alga
            # as "– [Saccharina japonica [Syn. Laminaria jaonica]] [STFCJ]":
            # strip the source tag and the only identity left is the binomial
            # inside the brackets. Read up to the next bracket of either kind so
            # a nested "[Syn. ...]" cannot be mistaken for the primary name.
            hit = re.search(r"\[\s*([^\[\]]+)", _BRACKET_TAG_RE.sub("", str(desc)))
            if hit:
                alt = hit.group(1)
        canon = _final_polish(_AVERAGE_TOKEN_RE.sub("", alt).lower())
    # Spelling only, and last, so that no upstream rule sees the hyphenated form:
    # applied earlier it turns "Cereal, ready to eat, ..." into FDC's category
    # prefix and the prefix strip then promotes the product line behind it.
    canon = _READY_TO_RE.sub(lambda m: "ready-to-" + m.group(1).lower(), canon)
    # ... and again on the finished canon. In the description the head can sit
    # too deep for either rule to reach it - CNF writes "Grains, rice, white,
    # glutinous, flour" and Phenol-Explorer "Orange [Blond], juice from
    # concentrate" - but by here the two-chunk rule has brought it to the front.
    canon = _MATERIAL_HEAD_RE.sub(_material_head_swap, canon)
    canon = _fold_material_tail(canon)
    for _rx, _rep in _COMPOUND_SPELLINGS:
        canon = _rx.sub(_rep, canon)
    canon = _CHEESE_ORDER_RE.sub(lambda m: f"cheese, {m.group(1).lower()}", canon)
    canon = _MEAT_PART_JOIN_RE.sub(
        lambda m: f"{_meat_name(m.group(1))}, {_fold_plural(m.group(2).lower())}", canon)
    canon = _MEAT_PART_SWAP_RE.sub(
        lambda m: f"{_meat_name(m.group(2))}, {_fold_plural(m.group(1).lower())}", canon)
    _sg = _SPECIES_GENUS_RE.match(canon)
    if (_sg and not _SPECIES_GENUS_KEEP.match(_sg.group(1))
            and not _SPECIES_CONNECTIVE_RE.search(_sg.group(1))):
        canon = f"{_sg.group(2)}, {_sg.group(1)}{canon[_sg.end():]}"
    # merges first (mechanical, generated), then the hand-curated overrides,
    # which therefore always have the final say.
    canon = _consult_tables(canon)
    # last, and after the tables, so a curated value spelling it "mince" folds
    # in exactly as a rule-produced one does
    canon = _GROUND_HEAD_RE.sub(lambda m: f"{m.group(1).lower()}, ground", canon)
    canon = _MINCE_RE.sub(lambda m: f"{m.group(1).lower()}, ground", canon)
    canon = _BREAD_TYPE_RE.sub(lambda m: f"{m.group(1).lower()} bread", canon)
    canon = re.sub(r"(?:\s*,\s*)+", ", ", canon).strip(" ,")
    # the fold can print the label the canon already carried
    _seen, _out = set(), []
    for _c in canon.split(", "):
        if _c in _seen and _c in _PRINTED_LABELS:
            continue
        _seen.add(_c); _out.append(_c)
    return ", ".join(_out)


def _consult_tables(canon: str) -> str:
    """_CANON_MERGES first, then the hand-curated _CANON_OVERRIDES."""
    hit = _CANON_MERGES.get(canon)
    if hit is None:
        hit = _CANON_OVERRIDES.get(canon)
    if hit is not None:
        return _CANON_OVERRIDES.get(hit, hit)
    # Both tables are keyed on the canon the rules make, LABELS AND ALL, so a
    # new axis kills every curated entry for the foods it labels - silently,
    # because the entry is still there. Adding the fat percentage turned
    # 'cheese hard' into 'cheese hard, 24% fat' and the curated rename to
    # 'cheese, hard' stopped firing on forty canons at once.
    #
    # So the tables are consulted a second time with the labels lifted off and
    # then put back. Only the labels _append_states actually printed on THIS
    # canon are peeled, which is what stops the fallback reaching a shorter
    # name that means something else: "potato" is not "potato, diced" with a
    # label removed, and no label was printed there to peel.
    if _PRINTED_LABELS:
        parts = [c.strip() for c in canon.split(",")]
        tail = []
        while len(parts) > 1 and parts[-1].lower() in _PRINTED_LABELS:
            tail.insert(0, parts.pop())
            base = ", ".join(parts)
            hit = _CANON_MERGES.get(base)
            if hit is None:
                hit = _CANON_OVERRIDES.get(base)
            if hit is not None and _is_rename(base, hit):
                hit = _CANON_OVERRIDES.get(hit, hit)
                # ...unless the curated name already says it: "cream, whipping"
                # -> "whipping cream" must not become "whipping cream, 38% fat,
                # 38% fat" when the entry itself carries the label.
                _said = {c.strip().lower() for c in hit.split(",")}
                keep = [t for t in tail if t.lower() not in _said]
                # ...and the curated name can say it in words rather than as a
                # chunk: "apple juice, light, fortified with vitamin C" already
                # carries the claim, so re-appending "enriched" said it twice.
                if _FORT_SAID_RE.search(hit):
                    keep = [t for t in keep if t.lower() != "enriched"]
                out = ", ".join([hit] + keep) if keep else hit
                # Consulted once more, exactly: the name the fallback builds can
                # itself be a curated key ("red kidney bean, dried"), and without
                # this it stood as its own canon beside the form it merges into.
                # Once only - the tables must not be walked recursively.
                hit2 = _CANON_MERGES.get(out) or _CANON_OVERRIDES.get(out)
                return _CANON_OVERRIDES.get(hit2, hit2) if hit2 else out
    return canon


def _canonicalize_core(desc: str) -> str:
    _PRINTED_LABELS.clear()
    if not desc: return ""
    d = str(desc).strip().translate(_FULLWIDTH_MAP)
    _desc0 = d
    # "Starch," / "Sugars," / "Sweets," head an FDC lab row exactly as "Minerals,"
    # does, but each is also a real food, so the sample code is the
    # discriminator - and this has to run while the head is still at the FRONT.
    # Applied later it would strip the food out of "Minerals, Sugar, Granulated,
    # White - NFY040XEG", where the panel is "Minerals" and the food IS sugar.
    if _CODED_PANEL_RE.match(d) and (_NF_SUFFIX_RE.search(d) or _NF_BARE_RE.search(d)):
        _nd = _CODED_PANEL_RE.sub("", d)
        if _nd.strip(" ,"):
            d = _nd
    d = _OPTIONAL_PLURAL_RE.sub("s", d)
    d = _QUOTE_RE.sub("", d)
    # Semicolons separate chunks exactly as commas do, but only commas were
    # split on, so "Tuna; chunk light; canned in water; drained solids" was
    # never chunked and came out as 'tuna; chunk light; ; solid' - the hollow
    # ";  ;" is the packing-medium rule emptying a chunk it could not remove.
    d = d.replace(";", ",")
    if _LAB_ROW_RE.search(d):
        d = _PANEL_FOOD_PREFIX_RE.sub("", d)
    d = _LOCAL_GLOSS_RE.sub(r", \1", d)
    d = _FOOTNOTE_STAR_RE.sub("", d)
    # A trademark mark is not part of the food's name; six canons ended in one
    # ("soy protein kebab \u00ae, fried").
    d = _TRADEMARK_RE.sub(" ", d)
    # Bookkeeping and provenance run before the synonym/rewrite loop: each can
    # sit at the head of the string and would otherwise become the canon.
    d = _DB_MARKER_RE.sub("", d)
    d = _PROVENANCE_RE.sub("", d)
    if _NA_HEAD_RE.match(d):
        d = _NA_HEAD_RE.sub("unidentified plant", d, count=1)
    for _rx, _rep in _SYNONYM_RE:
        d = _rx.sub(_rep, d)
    for _rx, _rep in _NAME_REWRITES:
        d = _rx.sub(_rep, d)
    for _rx, _rep in _ALTERNATION_RE:
        d = _rx.sub(_rep, d)
    d = _DAIRY_ALT_RE.sub("", d)
    d = _PLANT_ALT_BASE_RE.sub(r"plant-based \1, ", d)
    d = _PLANT_ALT_RE.sub("plant-based ", d)
    d = _SPREAD_LIKE_RE.sub("margarine", d)
    d = _BRAND_RE.sub(" ", d)
    d = _MARKETING_RE.sub(" ", d)
    d = _HEDGE_RE.sub("", d)
    d = _SERVING_BASIS_RE.sub(" ", d)
    d = _DIAMETER_RE.sub("", d)
    d = _EG_MARKER_RE.sub("", d)
    d = _GRADE_WORD_RE.sub(" ", d)
    if re.search(r"\bmargarine\b", d, re.I):
        d = _REDUNDANT_SPREAD_RE.sub("", d)
    d = _FLAVOUR_FILLER_RE.sub("", d)
    d = _WO_ABBREV_RE.sub("without", d)
    d = _W_ABBREV_RE.sub("with", d)
    d = _LIGHT_MEAT_ONLY_RE.sub(r"\1 meat, meat only", d)
    d = _ANALYZED_YEAR_RE.sub("", d)
    d = _strip_venue(d)
    d = _strip_panel_prefix(d)
    # Strip FDC Foundation Foods nutrient-panel group prefix FIRST
    # ("Proximates, ", "Beverages, ", "Cereals ready-to-eat, ", "Cereals, ")
    # - must run before _NF_SUFFIX_RE because the NF suffix regex is greedy
    # and can eat the entire informative middle of long descriptions.
    d = _FDC_GROUP_PREFIX_RE.sub("", d)
    d = _ALCOHOL_CLASS_RE.sub("", d)
    d = _WILD_RICE_RE.sub("wild rice", d)
    d = _SPLIT_PCT_RE.sub(lambda m: f", {m.group(2)}{m.group(1)}% fat", d)
    if not _NATIVE_IS_SPECIES_RE.match(d):
        d = _NATIVE_RE.sub("", d)
    if not _MEAT_DISH_HEAD_RE.match(d):
        # the whole animal says nothing a bare "duck" does not; the flesh does
        d = _MEAT_CHUNK_RE.sub(" " if _SKIN_INCLUDED_RE.search(d) else " meat only",
                               d, count=1)
    for _ in range(3):        # "Chole, FA - Beef" stacks two analytes
        nd = _ANALYTE_PREFIX_RE.sub("", d)
        if nd == d or not nd.strip():
            break
        d = nd
    # after the analyte prefix, so the qualifier is actually at the head
    d = _split_plural_head(d)
    d = _FRONTED_QUALIFIER_RE.sub(lambda mm: f"{mm.group(2)} {mm.group(1)}", d)
    d = _CATEGORY_HEAD_RE.sub("", d)
    d = _MATERIAL_HEAD_RE.sub(_material_head_swap, d)
    d = _fold_material_tail(d)
    d = _READY_RAW_RE.sub("", d)
    d = _SUPPLY_FORM_RE.sub("", d)
    d = _POULTRY_CLASS_RE.sub(r"\1", d)
    # Frida files the class term behind the ORGAN rather than the bird, where
    # it IS the animal: "Liver, broiler or fryer, raw" is chicken liver, and
    # dropping the term outright put it in the generic 'liver' canon beside
    # beef liver.
    d = _POULTRY_ALONE_RE.sub(", chicken", d)
    # Again, now that any panel/group prefix is gone: 58 rows are filed as
    # "Beverages, WENDY'S, tea" or "Minerals, Chinese Restaurant, Fried Rice",
    # where the venue sits behind the prefix and the first pass could not see it.
    d = _strip_venue(d)
    # Detect composition-changing state on a SEPARATE probe taken here, before
    # _NF_SUFFIX_RE and _PREP_RE run; `d` itself is untouched, so no existing
    # canon name shifts. The order matters twice over: _PREP_RE deletes the
    # state token outright, and _NF_SUFFIX_RE is greedy back to the first
    # capitalised chunk, so "Figs, Dried, Pass 2, Region 3, ... - NFY010CFL"
    # collapses to "Figs" and takes the state with it. That path alone hid the
    # state of 238 catalog rows, freeze-dried romaine and freeze-dried Bartlett
    # pears among them. The bracket strips are applied to the probe so a source
    # tag can never be read as a state.
    _probe = _BRACKET_TAG_RE.sub("", _BRACKET_VARIANT_RE.sub(" ", d))
    # The medium is read FIRST and struck out of the probe the state is read
    # from. "Peach, canned in pear juice" is canned peach, not peach juice, and
    # the juice state fired on the medium: the canon came out 'peach juice'.
    _pack = _detect_pack(_probe)
    # The wording is blanked whether or not the LABEL survived. _detect_pack
    # rightly refuses to label a DRAINED food by a medium that has been poured
    # away, but the state probe then read the medium as the food's own state:
    # "Apricot, canned in pear juice, drained" came out as 'apricot juice'.
    _pw = _pack[1] or _pack_wording(_probe)
    _sprobe = (re.sub(rf"\b{re.escape(_pw)}\b", " ", _probe, flags=re.I)
               if _pw else _probe)
    _state = _detect_state(_sprobe)
    _salt = _detect_salt(_probe)
    _sugar = _detect_sugar(_probe)
    _fat = _detect_fat(_probe) or ("", "")
    _whole = _detect_fat_whole(_probe)
    if _whole[0]:
        # On a dairy head the word IS the fat statement, so it goes whether or
        # not the source also gave a figure. Frida gives both - "Milk, whole,
        # 3.5, (UHT), % fat" - and leaving the word behind printed the claim
        # twice, as 'milk, whole, 3.5% fat' beside 'milk, full-fat'. The figure
        # wins where there is one; the label stands in where there is not.
        d = _WHOLE_CHUNK_RE.sub(", ", d, count=1)
        if not _fat[0]:
            _fat = _whole
    _fort = _detect_fortify(_probe)
    _grade = _detect_oil_grade(_probe)
    # read off d, not the probe: Phenol-Explorer writes the colour as a bracket
    # variant ("Common bean [Black]") and the probe has already lost it
    _colour = _detect_colour(d)
    _organ = _detect_organ(d)
    _cut = _detect_cut(d) if not _organ[0] else ("", "")
    _flav = _detect_flavour(d)
    # read off the ORIGINAL description: the brand chunk that carries the grade
    # ("SARGENTO SHARP") has already been stripped off d by this point
    _matur = _detect_maturity(_desc0)
    _moist = _detect_moisture(_probe)
    # read off the probe, like the fat figure: the source tag must not be in
    # reach of the percentage, and the bracket variants are already gone
    _fibre = _detect_fibre(_probe)
    _immat = _detect_immature(_probe)
    # ...but ripeness off the ORIGINAL description: FDC writes it in the caps
    # brand chunk ("BANANAS, SLIGHTLY RIPE, MEDIUM SIZE"), which _strip_caps_brand
    # has already taken off d by the time the probe is built
    _ripe = _detect_ripeness(_desc0)
    _decaf = _detect_decaf(_probe)
    _rgrain = _detect_refined_grain(_probe)
    _gfree = _detect_gluten_free(_probe)
    _part = _detect_part(_probe)
    _trim = _detect_trim(_probe)
    # read off the ORIGINAL description: BioFoodComp writes "Cod, wild, dorsal
    # muscle" and the chunk can be gone by the time _probe is built
    _prov = _detect_provenance(_desc0)
    # read AFTER _GENERIC_OIL_RE has struck the unnamed fats, so only a fat the
    # source actually named can reach the label
    _cfat = _detect_cook_fat(_GENERIC_OIL_RE.sub(" ", _probe))
    # read off _probe above; drop the wording so it cannot become a chunk
    d = _AXIS_PHRASE_RE.sub("", d)
    d = _PACK_PHRASE_RE.sub("", d)
    d = _UNFORTIFIED_PHRASE_RE.sub("", d)
    d = _GENERIC_OIL_RE.sub("", d)
    # only the axes that actually FIRED lose their wording: "Sugar, refined" is
    # not a pressing grade, and deleting the word there leaves bare 'sugar'
    if _matur[1]:
        d = d.replace(_matur[1], ", " if _matur[1].startswith(",") else " ", 1)
    if _cut[1]:
        d = d.replace(_cut[1], ", " if _cut[1].startswith(",") else " ", 1)
    if _organ[1]:
        d = d.replace(_organ[1], ", " if _organ[1].startswith(",") else " ", 1)
    if _colour[1]:
        d = d.replace(_colour[1], ", " if _colour[1].startswith(",") else " ", 1)
    if _state[0] == "dried":
        # "powder" IS the dried state here, so carrying both spellings split one
        # food across up to five canons - turmeric had 'turmeric, dried',
        # 'turmeric, powder', 'turmeric, powdered, dried', 'turmeric, ground'
        # and 'turmeric, ground, dried'. Only a WHOLE chunk goes: "chili powder"
        # and "cocoa powder" keep the word, because there it is the name.
        d = _POWDER_CHUNK_RE.sub("", d)
    if _trim[0]:
        d = _TRIM_PHRASE_RE.sub("", d)
    if _fat[0]:
        d = _FAT_PHRASE_RE.sub("", d)
    if _fort[0]:
        d = _FORT_PHRASE_RE.sub("", d)
    if _grade[0]:
        d = _GRADE_PHRASE_RE.sub("", d)
    # Each new axis loses its wording once it has been READ, exactly as the fat
    # and fortification phrasings do: _append_state declines a label whose own
    # token is still standing in the name, so leaving the words in place would
    # print the source's spelling instead of the house label and split the canon
    # in two again ('crispbread, 17% fibre' beside 'crispbread, rye').
    if _fibre[0]:
        d = _FIBRE_PHRASE_RE.sub("", d)
    if _immat[1]:
        d = d.replace(_immat[1], ", " if _immat[1].startswith(",") else " ", 1)
    if _ripe[1]:
        d = re.sub(rf"\s*,?\s*{re.escape(_ripe[1])}\b", "", d, count=1, flags=re.I)
    if _decaf[1]:
        d = re.sub(rf"\s*,?\s*{re.escape(_decaf[1])}\b", "", d, count=1, flags=re.I)
    # whole chunk only, for the reason set out on _strip_whole_chunk
    if _rgrain[1]:
        d = _strip_whole_chunk(d, _rgrain[1])
    if _gfree[1]:
        d = _strip_whole_chunk(d, _gfree[1])
    # Only the FARMED half loses its wording. Striking "wild" would rename the
    # species that are called wild - 'wild rice' would become 'rice, wild' and
    # sit beside the grain - so it is left in place and _append_state declines
    # the label wherever the word is still standing.
    if _prov[0] == "farmed":
        d = _FARMED_PHRASE_RE.sub("", d)
    if _cfat[1]:
        # the cooking VERB has to survive - it is the state - so only the
        # "in <fat> oil" tail of the phrase is struck
        d = re.sub(r"\s*\b(?:in|with)\s+[a-z]+\s+(?:oils?|fats?)\b", "", d,
                   count=1, flags=re.I)
    if _pack[1]:
        # only the wording that was actually READ as a medium goes; a "with
        # water" that failed the packed-context test is still a real ingredient.
        # The medium is often NAMED - "canned in pear juice", "canned in olive
        # oil" - and the qualifier has to go with it, or _append_state sees its
        # own token still in the name and prints no label at all, leaving
        # 'peach juice' and 'tuna, in olive oil' beside 'tuna, in oil'.
        d = re.sub(rf"\s*,?\s*{re.escape(_pw)}\b", "", d, flags=re.I)
    elif _pw:
        # ...and a DECLINED medium loses its wording too. Only the packed
        # phrasings are touched, so a genuine ingredient ("porridge made with
        # water") is never at risk.
        d = re.sub(rf"\s*,?\s*{re.escape(_pw)}\b", "", d, flags=re.I)
    # Strip FDC NF panel-fragmentation suffix (two flavors).
    _had_code = bool(_NF_SUFFIX_RE.search(d) or _NF_BARE_RE.search(d))
    d = _NF_SUFFIX_RE.sub("", d)
    d = _NF_BARE_RE.sub("", d)
    if _had_code:
        d = _strip_caps_brand(d)
    # Strip trailing source tag e.g. " [BioFoodComp]" - keep variant info inside
    # the brackets out of the canon name.
    d = _BRACKET_TAG_RE.sub("", d)
    # Then strip any remaining bracketed variant tags anywhere in the description:
    # color / cultivar / sci-name tags that Phenol-Explorer + BioFoodComp use
    # ("Tea [Oolong]", "Orange [Blond]", "Common cabbage [Purple]", etc.).
    d = _BRACKET_VARIANT_RE.sub(" ", d)
    d = _AVERAGE_TOKEN_RE.sub("", d)
    # Strip preparation / state tokens.
    d = _PREP_RE.sub("", d)
    # Strip parenthetical clarifications and quantitative qualifiers - these
    # are non-essential annotations that prevent head-aware grouping when
    # they appear in the second comma-chunk of preserve-head foods. For
    # strip-list heads the second chunk is dropped anyway, so this is a
    # no-op there; the effect is on preserve heads ("potato, white
    # (industrial), 50% extraction" → "potato, white").
    # Promote ingredient-difference parentheticals to a comma chunk BEFORE the
    # blanket paren strip below, and after _PREP_RE so the promoted text keeps
    # words _PREP_RE would otherwise eat ("with sauce" is one of its tokens).
    d = _PAREN_QUALIFIER_RE.sub(r", \1", d)
    d = _PAREN_CONTENT_RE.sub(" ", d)
    d = _QUANT_RE.sub("", d)
    d = _SAMPLE_CODE_RE.sub("", d)
    # Cleanup: drop empty parens left by stripped tokens, collapse repeated
    # commas/spaces, trim punctuation.
    d = _EMPTY_PARENS_RE.sub(" ", d)
    # Unclosed parenthetical. _NF_SUFFIX_RE can cut INSIDE a parenthesis, because a
    # comma is one of its anchors: "KRAFT 100% (AL,CA1) - NFY120DQP" loses
    # ",CA1) - NFY120DQP" and leaves "(AL", which _PAREN_CONTENT_RE cannot match
    # because there is no closing bracket left to match it against.
    if d.count("(") > d.count(")"):
        d = re.sub(r"\s*\([^()]*$", "", d)
    if d.count(")") > d.count("("):
        d = d.replace(")", " ")
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
    chunks = [c for c in chunks
              if not _PUNCT_ONLY_RE.match(c) and not _LABEL_ONLY_RE.match(c)
              and not _VENUE_CHUNK_RE.match(c)]
    # ...and the grain form, on a cereal head only. See _GRAIN_FORM_RE.
    if chunks and _CEREAL_CONTEXT_RE.search(chunks[0]):
        chunks = chunks[:1] + [c for c in chunks[1:] if not _GRAIN_FORM_RE.match(c)]
    # ...but never the FIRST chunk: "Ice cream bar or stick, chocolate coated"
    # is an ice cream, and dropping the head left the canon 'chocolate coated'.
    chunks = chunks[:1] + [c for c in chunks[1:] if not _FORM_ALTERNATION_RE.match(c)]
    chunks = [c.strip(" /") for c in chunks]
    chunks = [_CONN_TRAIL_RE.sub("", c) for c in chunks
              if c and not _CONN_ONLY_RE.match(c)]
    chunks = [c for c in chunks if c]
    if not chunks:
        return d
    chunks = _drops_restatement(chunks)
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
        # First NON-EMPTY part: the prep strip can leave a leading slash
        # ("Pieces/chunks vegetarian based on soya" loses both "pieces" and
        # "chunks"), and splitting that took the empty left side.
        parts = [x.strip() for x in chunks[0].split("/") if x.strip()]
        if parts:
            chunks[0] = parts[0]
    # A strip head still yields to a protected plant part (see _PART_PRESERVE):
    # dropping it would merge two different organs of the same plant.
    _protected = (len(chunks) > 1
                  and chunks[1] in _PART_PRESERVE.get(chunks[0], frozenset()))
    if chunks[0] in _CULTIVAR_STRIP_HEADS and not _protected:
        return _append_states(_fold_plural(chunks[0]), _flav, _matur, _cut, _organ, _colour, _state, _trim, _pack, _salt, _sugar, _fat, _fort, _moist, _part, _grade, _prov, _cfat, _fibre, _immat, _ripe, _decaf, _rgrain, _gfree)
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
    return _append_states(", ".join(_fold_plural(c) for c in [chunks[0]] + kept[:1]), _flav, _matur, _cut, _organ, _colour, _state, _trim, _pack, _salt, _sugar, _fat, _fort, _moist, _part, _grade, _prov, _cfat, _fibre, _immat, _ripe, _decaf, _rgrain, _gfree)


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


def _clean_variant_name(desc: str) -> str:
    """Tidy a raw description for use as a canon name.

    refine_canons_by_nutrition falls back to the description when a nutritional
    outlier canonicalizes straight back to its parent - the detail is the point,
    since it is what makes the split visible. Taking it verbatim also took the
    source tags with it, so 750 canon names carried "[stfcj]", "[wafct]" or a
    binomial. Only the tag noise is removed here; the descriptive middle stays.

    The FDC sample code goes too, for the same reason the tags do: it is provenance,
    not identity, and it left canons like "haddock, triad 2 - cy060r4".

    The bracket strip loops because a few STFCJ names nest them
    ("[Sargassum fusiforme [Syn. Hizikia fusiformis]]"), which a single
    non-recursive pass leaves half-eaten, and any orphaned bracket is dropped
    afterwards.
    """
    d = _NF_BARE_RE.sub("", _NF_SUFFIX_RE.sub("", _QUOTE_RE.sub("", str(desc))))
    for _ in range(4):
        nd = _BRACKET_VARIANT_RE.sub(" ", _BRACKET_TAG_RE.sub("", d))
        if nd == d:
            break
        d = nd
    d = d.replace("[", " ").replace("]", " ")
    d = re.sub(r"\s+", " ", d)
    d = re.sub(r"\s*,\s*", ", ", d)
    d = re.sub(r"(?:,\s*)+,", ",", d)
    return d.strip(" ,-").lower()


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
                out_canon = _clean_variant_name(desc)
                if not out_canon or out_canon == cname:
                    out_canon = f"{cname} (variant {f})"
            refined.setdefault(out_canon, []).append(f)
            n_split += 1

    print(f"[refine] {n_split} cultivar variants split off as their own canons "
          f"(profile deviation > {_OUTLIER_RATIO}× group median); "
          f"{n_exempt} strip-list-head groups exempted from gate.", flush=True)
    return refined


class _FoodLookup:
    """Compact stand-in for the dense fdc_id-indexed arrays.

    canon_arr, plant_arr and spice_arr used to be dense over the fdc_id space. That
    space runs to 48,002,027 because fdc_ids are block-allocated per source, while
    only ~128,000 slots are ever populated - 0.27%. The three arrays cost 288 MB to
    carry 1.3 MB of information, and every worker process paid it.

    This keeps the live ids sorted with one parallel array per attribute and resolves
    a whole batch with searchsorted - O(log n) vectorised rather than the O(1) gather
    it replaces, which on a 50,000-row batch is not measurable next to the parquet read.

    An fdc_id the index does not know resolves to canon -1, exactly as a dense slot
    that was never filled did, so the `(c != -1)` guard downstream is unchanged. The
    old `f < max_f` bound is no longer needed either: an id past the end misses instead
    of indexing out of range.
    """
    __slots__ = ("ids", "canon_of", "plant_of", "spice_of")

    def __init__(self, ids, canon_of, plant_of, spice_of):
        self.ids = ids; self.canon_of = canon_of
        self.plant_of = plant_of; self.spice_of = spice_of

    def resolve(self, f):
        """(canon, is_spice, is_plant) for a batch of fdc_ids. canon is -1 if unknown."""
        if len(self.ids) == 0:
            z = np.zeros(len(f), dtype=bool)
            return np.full(len(f), -1, dtype=np.int64), z, z
        pos = np.searchsorted(self.ids, f)
        np.clip(pos, 0, len(self.ids) - 1, out=pos)
        hit = self.ids[pos] == f
        return (np.where(hit, self.canon_of[pos], -1),
                hit & self.spice_of[pos],
                hit & self.plant_of[pos])


def _compact_food_arrays(c_arr, pl_arr, sp_arr):
    """Dense fdc-indexed arrays -> the sorted-id form _FoodLookup reads."""
    ids = np.nonzero(c_arr != -1)[0].astype(np.int64)
    return ids, c_arr[ids].astype(np.int64), pl_arr[ids].copy(), sp_arr[ids].copy()


def _lookup_from_static(db) -> "_FoodLookup":
    """Build the lookup from either pickle form.

    An index written before the compaction still carries the dense arrays; rather
    than force a rebuild for a format change alone, convert it on load. The identity
    stamp still governs whether the index is VALID - this governs only its shape.
    """
    if "food_ids" in db:
        return _FoodLookup(db["food_ids"], db["food_canon"], db["food_plant"], db["food_spice"])
    return _FoodLookup(*_compact_food_arrays(db["canon_arr"], db["plant_arr"], db["spice_arr"]))


def build_static_food_meta(args, out_path):
    """Build a dense numpy index over the food catalog. Adapted from the query prototype."""
    print("[*] Building static food meta...", flush=True)
    nut = load_smart(args.nutrient).rename(columns={"id":"nutrient_id","name":"nutrient_name"})
    # A nutrient CAN have no unit, and 19 do. They are the ids minted for a
    # different measurement BASIS by the ingest unit reconciliation - AFCD's
    # amino acids reported "(mg/gN)", McCance's "Tryptophan/60", one "/100g fa"
    # fatty acid - which cannot be scaled onto per-100 g and so carry no
    # canonical unit at all. Left as NA they reached .upper() as None and killed
    # the whole index build. No unit means no scaling, which is the same 1.0 the
    # else-branch below already applies to everything that is neither G nor UG.
    unit_map = (nut.set_index("nutrient_id")["unit_name"]
                .astype("string").str.upper().fillna("").to_dict())
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
    # Drop the source columns that are not composition values. They have to come out
    # AFTER the store scan above, which re-adds every id it finds and would otherwise
    # put them straight back. Each one is modeled by construction - it sits in
    # nutrient.csv like any other - and each would then be counted as a nutrient the
    # food supplies: `Latest Revision in Version = 2.1` is carried by all 10,133
    # BioFoodComp foods and lands in model_count and model_mass as 2.1 mg of
    # something. See food_DBs/_common/non_nutrients.py.
    _bogus = modeled & set(NON_NUTRIENT_IDS)
    if _bogus:
        modeled -= _bogus
        print(f"[*] Modeled: excluded {len(_bogus)} non-composition source columns "
              f"(identifiers, version stamps, conversion factors, as-purchased yields).",
              flush=True)
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
    # Exact re-listings: the same food written twice with the same values. Dropped here
    # for the same reason the export drops them - a food the released table does not
    # contain must not be scorable either - and by the same shared rule, so the two
    # cannot reach different verdicts. Only the candidates are read out of the store
    # (a duplicated description is necessary for a match), which is ~19% of eligible
    # foods rather than the whole 28M-row store.
    _dt_all = food_df["data_type"].astype("string").fillna("").str.lower().str.strip().str.replace(" ", "_")
    _keep = np.ones(len(food_df), dtype=bool)
    if DROP_BRANDED:  _keep &= (_dt_all != "branded_food").to_numpy(dtype=bool)
    if DROP_MODELLED: _keep &= (_dt_all != "survey_fndds_food").to_numpy(dtype=bool)
    _elig = food_df.loc[_keep, ["fdc_id", "description", "data_type"]].copy()
    _elig["fdc_id"] = pd.to_numeric(_elig["fdc_id"], errors="coerce")
    _elig = _elig.dropna(subset=["fdc_id"]).astype({"fdc_id": "int64"})
    _relisted, _cand = set(), relisting_candidates(_elig)
    if _cand:
        _rows = []
        for _f in sorted(Path(args.food_nutrient).glob("*/*.parquet")):
            _t = pq.read_table(_f, columns=["fdc_id", "nutrient_id", "amount"],
                               filters=[("fdc_id", "in", sorted(_cand))]).to_pandas()
            # The non-composition columns have to go before the value vectors are
            # compared, not merely before scoring. 572 BioFoodComp foods carry the
            # version stamp and nothing else, so on the raw store their whole vector
            # is "96926:2.1" and they all look like re-listings of each other: the
            # count comes out at 175 instead of 74, and 101 real foods are discarded.
            _t = _t[~_t["nutrient_id"].isin(NON_NUTRIENT_IDS)]
            if len(_t):
                _t["source_code"] = source_of_bucket_file(_f)
                _rows.append(_t)
        if _rows:
            _relisted = find_exact_relistings(_elig, pd.concat(_rows, ignore_index=True))
        del _rows
    del _elig, _cand, _dt_all, _keep

    meta = food_df.set_index(pd.to_numeric(food_df["fdc_id"], errors="coerce"))[["description","data_type","food_category_id"]].to_dict(orient="index")
    max_fdc = int(max(k for k in meta if pd.notna(k))) + 1000
    c_arr = np.full(max_fdc, -1, dtype=np.int32)
    pl_arr = np.zeros(max_fdc, dtype=bool); sp_arr = np.zeros(max_fdc, dtype=bool)
    drop_cats = set(str(x).strip() for x in args.drop_category); food_stats = {}
    _paper_title_skipped = _copy_skipped = _relist_skipped = 0
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
        if not _MILLING_EXTRACTION_RE.search(desc_lc) and (
                (len(desc_lc) > 60 and _PAPER_TITLE_RE.search(desc_lc)) or
                (dt_norm == "experimental_food"
                 and len(desc_lc) > _EXPERIMENTAL_TITLE_LEN)):
            _paper_title_skipped += 1
            continue
        # Source editing-layer copies, dropped for the same reason the export drops
        # them: the name is the original's and the values are not.
        if _COPY_RECORD_RE.search(desc_lc):
            _copy_skipped += 1
            continue
        if fid in _relisted:
            _relist_skipped += 1
            continue
        if (DROP_BRANDED and dt_norm=="branded_food") or (DROP_MODELLED and dt_norm=="survey_fndds_food") or cat in drop_cats or cat in ALWAYS_DROP_CATS or WHALE_RE.search(desc_lc):
            continue
        c_arr[fid] = fid
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
        # Category flags are read HERE, below the fallback chain, not above it.
        # Only FDC and a few sources ship a food_category_id; for everything
        # else - Phenol-Explorer, BioFoodComp, AFCD, STFCJ, MERIDA, CNF - `cat`
        # is empty until one of the branches above fills it. Reading the flags
        # any earlier tests them against "" and silently answers False.
        # Measured on the shipped index before this moved: 317 of the 335 foods
        # that end up categorised "Spices and Herbs" were left unflagged, so
        # --allow-spices held back 18 of them and 378 Phenol-Explorer spice rows
        # (clove 265, rosemary 76, cumin 16, coriander 13) ranked in the cohort
        # differential tables with the flag OFF. Clove at 100 g is not a food.
        # The OR keeps the explicit pl_arr=True the Phenol-Explorer and
        # BioFoodComp branches set: those two are plant sources whatever their
        # category resolves to, and that patch is what kept this bug hidden on
        # the plant side while the spice side had no equivalent.
        sp_arr[fid] = (cat == "Spices and Herbs")
        pl_arr[fid] = bool(pl_arr[fid]) or (cat in PLANT_CATS)
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
    # c_arr[all variants] -> rep_fid so both readers fold their nutrients
    # together - build_modeled_index by MEAN over the variant set,
    # score_one_bacterium by MAX over it.
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
        u = str(unit_map.get(n) or "").upper()
        u_arr[n] = 1000.0 if u == "G" else (0.001 if u == "UG" else 1.0)
    _ids, _can, _pl, _sp = _compact_food_arrays(c_arr, pl_arr, sp_arr)
    print(f"[*] Food lookup: {len(_ids):,} live ids of an {len(c_arr):,}-wide fdc space "
          f"({(c_arr.nbytes+pl_arr.nbytes+sp_arr.nbytes)/1e6:.0f} MB dense -> "
          f"{(_ids.nbytes+_can.nbytes+_pl.nbytes+_sp.nbytes)/1e6:.1f} MB).", flush=True)
    pickle.dump({"food_ids":_ids,"food_canon":_can,"food_plant":_pl,"food_spice":_sp,
                 "unit_arr":u_arr,
                 "food_stats":rep_stats,"modeled":modeled,
                 "nut_name":dict(zip(nut["nutrient_id"], nut["nutrient_name"].astype(str))),
                 "rep2can":{f:c_arr[f] for f in rep_stats}}, open(out_path, "wb"))
    print("[*] Static food meta written.", flush=True)
    if _paper_title_skipped:
        print(f"[*] Literature filter: skipped {_paper_title_skipped} rows whose description "
              "column held a citation rather than a food name (almost all of FDC's "
              "experimental_food type, which is a bibliography).", flush=True)
    if _copy_skipped:
        print(f"[*] Editing-layer filter: skipped {_copy_skipped} source copy records "
              "(name of the original, values of something else).", flush=True)
    if _relist_skipped:
        print(f"[*] Re-listing filter: skipped {_relist_skipped} foods repeated verbatim "
              "(same source, same name, same complete value vector).", flush=True)

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
    fl = _lookup_from_static(static_db); u_arr = static_db["unit_arr"]
    max_u = len(u_arr)
    for bi, b in enumerate(scn.to_batches(), 1):
        f = b["fdc_id"].to_numpy().astype(np.int64, copy=False)
        n = b["nutrient_id"].to_numpy().astype(np.int32, copy=False)
        a = b["amount"].to_numpy().astype(np.float32, copy=False)
        m = (f>=0)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f, n, a = f[m], n[m], a[m]
        c, _is_sp, _ = fl.resolve(f); v = (c != -1) & (~_is_sp)
        if not v.any(): continue
        f, n, a, c = f[v], n[v], a[v], c[v]; a = a*u_arr[n]
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
        # NOTE the asymmetry: score_one_bacterium reduces the SAME canon by MAX,
        # not by mean. Deliberate. This table answers "what does this food
        # typically contain", where an outlying dried variant should not set the
        # number; the score answers "can this canon supply this nutrient at all",
        # where the richest member is the honest answer. Read one number as the
        # other and every worked example in this file stops making sense.
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
    STATIC_DB["food_lookup"] = _lookup_from_static(STATIC_DB)


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
    fl = STATIC_DB["food_lookup"]; u_arr = STATIC_DB["unit_arr"]
    max_u = len(u_arr)
    for b in scn.to_batches():
        f = b["fdc_id"].to_numpy().astype(np.int64, copy=False)
        n = b["nutrient_id"].to_numpy().astype(np.int32, copy=False)
        a = b["amount"].to_numpy().astype(np.float32, copy=False)
        m = (f>=0)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f, n, a = f[m], n[m], a[m]
        c, is_sp, is_pl = fl.resolve(f); v = (c != -1)
        if not DYNAMIC_STATE["asp"]: v &= ~is_sp
        if not v.any(): continue
        n, a, c = n[v], a[v], c[v]; is_sp, is_pl = is_sp[v], is_pl[v]; a = a*u_arr[n]
        # Spices only. Everything else is compared on the 100 g basis it is
        # stored and published on, with no rescaling at all.
        if DYNAMIC_STATE["asp"]: a = np.where(is_sp, a*SPICE_SERVING_WEIGHT, a)
        # Collapse each (nutrient, canon) to ONE value, and that value is the
        # MAXIMUM: the lexsort orders by amount within the run and the mask
        # below keeps the LAST element of it. A canon is therefore scored at
        # its richest member for every nutrient independently - which is why
        # an over-merged canon is a composition error and not merely a naming
        # one, and why the axes above exist. (build_modeled_index takes the
        # MEAN over the same members for the modeled tables; the two readers
        # disagree on purpose - see the note there.)
        idx = np.lexsort((a, c, n)); n, c, a, is_pl = n[idx], c[idx], a[idx], is_pl[idx]
        m2 = np.empty(len(n), dtype=bool); m2[-1] = True; m2[:-1] = (n[:-1] != n[1:]) | (c[:-1] != c[1:])
        n, c, a, ipl = n[m2], c[m2], a[m2], is_pl[m2]
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

    Missing paths are skipped rather than raising: not every input is required
    on every run.
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


def _canonicalizer_fingerprint() -> str:
    """Hash the rules that decide which fdc_ids share a canon.

    The input fingerprint below answers "was this index built from these files?".
    It cannot answer "was it built by this canonicalizer?" - and that failure is
    just as silent. Editing a strip-head list, a preparation regex or the state
    patterns changes which foods fold together while every input file stays byte
    identical, so an index cached from the previous rules is served unchanged and
    the run reports plausible numbers for a food universe the code no longer
    builds. Hashing the rule surface makes that edit invalidate the cache.
    """
    import hashlib, inspect
    parts = []
    for fn in (canonicalize_food_name, _canonicalize_core, _final_polish,
               _material_head_swap, _fold_material_tail, _split_plural_head,
               _strip_caps_brand,
               _fold_plural, _singularize,
               _detect_state, _detect_salt, _detect_sugar, _detect_fat, _detect_fortify,
               _detect_moisture, _detect_part, _detect_oil_grade, _detect_colour, _detect_organ, _detect_cut, _detect_flavour, _detect_maturity,
               _detect_pack, _detect_trim, _append_state, _append_states,
               _strip_whole_chunk, _is_rename,
               _detect_fibre, _detect_immature, _detect_ripeness, _detect_decaf,
               _detect_refined_grain, _detect_gluten_free,
               _is_spelling_alias, _edit_distance_le, _clean_variant_name,
               _strip_panel_prefix, _strip_venue):
        try:
            parts.append(inspect.getsource(fn))
        except (OSError, TypeError):
            parts.append(getattr(fn, "__name__", "?"))
    parts.append(repr(sorted(_CULTIVAR_STRIP_HEADS)))
    parts.append(repr(sorted(_CULTIVAR_STRIP_HEADS_HARD)))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _STATE_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _SALT_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _PROVENANCE_PATTERNS]))
    parts.append(_FARMED_FOOD_RE.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _SUGAR_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _FAT_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _FORTIFY_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _MOISTURE_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _RIPENESS_PATTERNS]))
    for _rx in (_FIBRE_PCT_RE, _HIGH_FIBRE_RE, _IMMATURE_RE, _DECAF_RE,
                _CEREAL_CONTEXT_RE, _REFINED_GRAIN_RE, _WHOLEGRAIN_RE,
                _GLUTEN_FREE_RE, _GRAIN_FORM_RE):
        parts.append(_rx.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _PART_PATTERNS]))
    parts.append(_LEAF_IS_DEFAULT_RE.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _OIL_GRADE_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _COLOUR_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _ORGAN_PATTERNS]))
    parts.append(_ORGAN_HEAD_SKIP_RE.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _CUT_PATTERNS]))
    parts.append(_MEAT_CONTEXT_RE.pattern)
    parts.append(_CUT_RETAIL)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _FLAVOUR_PATTERNS]))
    parts.append(_FLAVOURED_HEAD_RE.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _MATURITY_PATTERNS]))
    parts.append(_MATURITY_SKIP_RE.pattern)
    parts.append(_MATURITY_CONTEXT_RE.pattern)
    parts.append(_COLOUR_HEAD_SKIP_RE.pattern)
    parts.append(_WHITE_IS_DEFAULT_RE.pattern)
    parts.append(_OIL_CONTEXT_RE.pattern)
    parts.append(repr([(lab, rx.pattern) for lab, rx in _PACK_PATTERNS]))
    parts.append(repr([(lab, rx.pattern) for lab, rx in _TRIM_PATTERNS]))
    parts.append(_PACKED_CONTEXT_RE.pattern)
    parts.append(_MEDIUM_NO_CONTEXT_RE.pattern)
    parts.append(_COOKED_IN_RE.pattern)
    for rx in (_PREP_RE, _CULTIVAR_CODE_RE, _PACK_MEDIUM_RE, _PREPARED_WITH_RE,
               _DRY_MIX_RE, _BRACKET_VARIANT_RE, _BRACKET_TAG_RE,
               _VENUE_PREFIX_RE, _VENUE_GENERIC_RE, _VENUE_CHUNK_RE, _LABEL_ONLY_RE, _PUNCT_ONLY_RE,
               _PAPER_TITLE_RE, _MILLING_EXTRACTION_RE,
               _NF_SUFFIX_RE, _NF_BARE_RE, _FDC_GROUP_PREFIX_RE, _ALCOHOL_CLASS_RE,
               _WILD_RICE_RE, _SPLIT_PCT_RE, _NATIVE_RE, _NATIVE_IS_SPECIES_RE, _MEAT_DISH_HEAD_RE, _MEAT_CHUNK_RE, _SKIN_INCLUDED_RE,
               _DAIRY_HEAD_RE, _WHOLE_CHUNK_RE,
               _PAREN_QUALIFIER_RE, _SAMPLE_CODE_RE, _ANALYZED_YEAR_RE,
               _WO_ABBREV_RE, _W_ABBREV_RE, _AVERAGE_TOKEN_RE,
               # round 4
               _CORP_ENTITY_RE, _PROVENANCE_RE, _DB_MARKER_RE,
               _NA_HEAD_RE, _NA_CHUNK_RE, _READY_TO_RE, _STORAGE_CHUNK_RE, _AXIS_PHRASE_RE,
               _NULL_QUALIFIER_RE, _PLAIN_IS_DARK_RE, _PULP_CHUNK_RE,
               _FARMED_PHRASE_RE, _GENERIC_OIL_RE, _COOK_FAT_RE,
               _PACK_PHRASE_RE, _NAMED_OIL_RE, _FAT_PCT_RE, _NOT_FORTIFIED_RE, _POWDER_CHUNK_RE, _TRIM_PHRASE_RE, _INGREDIENT_CTX_RE, _UNFORTIFIED_PHRASE_RE, _FAT_PHRASE_RE, _FORT_PHRASE_RE, _GRADE_PHRASE_RE,
               _FRONTED_QUALIFIER_RE, _CATEGORY_HEAD_RE, _READY_RAW_RE,
               # round 5
               _PLURAL_HEAD_RE, _PLURAL_HEAD_SKIP_RE,
               _MATERIAL_HEAD_RE, _MATERIAL_TAIL_RE, _MATERIAL_TAIL_BLOCK_RE,
               _MATERIAL_TAIL_BEAN_RE, _SELF_NAMING_TAIL_RE, _RECONSTITUTED_RE,
               _MILK_BEVERAGE_RE,
               _SPECIES_GENUS_RE, _SPECIES_GENUS_KEEP, _SPECIES_CONNECTIVE_RE,
               _CHEESE_ORDER_RE, _MEAT_PART_JOIN_RE, _MEAT_PART_SWAP_RE):
        parts.append(rx.pattern)
    parts.append(repr(_POLISH_EDGE))
    parts.append(_FILLER_CHUNK_RE.pattern)
    parts.append(repr(sorted(_BRAND_NAMES)))
    parts.append(repr(sorted(_CANON_OVERRIDES.items())))
    parts.append(repr([(rx.pattern, rep) for rx, rep in _NAME_REWRITES]))
    parts.append(repr(sorted(_SYNONYM_PAIRS.items())))
    parts.append(repr(sorted(_ANIMAL_TO_MEAT.items())))
    parts.append(repr([(rx.pattern, rep) for rx, rep in _COMPOUND_SPELLINGS]))
    parts.append(repr(sorted(_CANON_MERGES.items())))
    # round 14 rule surfaces: each of these decides which fdc_ids share a canon
    # and none is reachable through the sources hashed above
    parts.append(repr(_ALTERNATION_PHRASES))
    parts.append(_FORM_ALTERNATION_RE.pattern)
    parts.append(_POULTRY_CLASS_RE.pattern)
    parts.append(_POULTRY_ALONE_RE.pattern)
    parts.append(_SUPPLY_FORM_RE.pattern)
    parts.append(_PLANT_ALT_BASE_RE.pattern)
    parts.append(_PLANT_ALT_RE.pattern)
    parts.append(_DAIRY_ALT_RE.pattern)
    parts.append(_MINCE_RE.pattern)
    parts.append(_GROUND_HEAD_RE.pattern)
    parts.append(_BREAD_TYPE_RE.pattern)
    parts.append(_LIGHT_MEAT_ONLY_RE.pattern)
    parts.append(_TRADEMARK_RE.pattern)
    parts.append(_CHUNK_PREP_RE.pattern)
    parts.append(_FORT_SAID_RE.pattern)
    parts.append(repr(sorted(_RESTATES_STOP)))
    parts.append(repr(_BRAND_NAMES))
    parts.append(_SPREAD_LIKE_RE.pattern)
    parts.append(_FLAVOUR_FILLER_RE.pattern)
    parts.append(_MARKETING_RE.pattern)
    parts.append(_REDUNDANT_SPREAD_RE.pattern)
    parts.append(_HEDGE_RE.pattern)
    parts.append(_DRY_STYLE_CONTEXT_RE.pattern)
    parts.append(_INGREDIENT_LIST_RE.pattern)
    parts.append(_ANALYTE_PREFIX_RE.pattern)
    parts.append(_LAB_ROW_RE.pattern)
    parts.append(_PANEL_FOOD_PREFIX_RE.pattern)
    parts.append(_ORIGIN_CHUNK_RE.pattern)
    parts.append(_MEAT_HEAD_RE.pattern)
    parts.append(repr(sorted(_SILENT_STATE_LABELS)))
    parts.append(repr(sorted(_IE_PLURALS.items())))
    parts.append(repr(sorted(_FOREIGN_NOT_PLURAL)))
    parts.append(repr(sorted(_IRREGULAR_SINGULAR.items())))
    parts.append(_ANDOR_TRAIL_RE.pattern)
    parts.append(_TRAIL_GRADE_RE.pattern)
    parts.append(_SERVING_BASIS_RE.pattern)
    parts.append(_DIAMETER_RE.pattern)
    parts.append(_LINE_NUMBER_RE.pattern)
    parts.append(_EG_MARKER_RE.pattern)
    parts.append(_GRADE_WORD_RE.pattern)
    parts.append(_FOOTNOTE_STAR_RE.pattern)
    parts.append(_LOCAL_GLOSS_RE.pattern)
    parts.append(repr(sorted(_FULLWIDTH_MAP.items())))
    parts.append(repr(sorted((k, sorted(v)) for k, v in _PART_PRESERVE.items())))
    parts.append(repr(sorted(_PANEL_HEADS)))
    parts.append(_CODED_PANEL_RE.pattern)
    parts.append(_CAPS_CHUNK_RE.pattern)
    parts.append(repr(_EXPERIMENTAL_TITLE_LEN))
    parts.append(repr(_OUTLIER_RATIO))
    parts.append(repr(sorted(_KEY_NUTRIENT_IDS)))
    return hashlib.sha256("\x00".join(parts).encode("utf-8")).hexdigest()[:16]


def _index_identity(*paths) -> dict:
    """Fingerprint the inputs a derived index was built from.

    An mtime comparison answers "is the index older than its inputs?" but not "was it built
    from THESE inputs?". Point --food_nutrient at a different store whose files predate the
    cached index and every staleness check passes while the index describes a different food
    universe. That is not hypothetical: it silently corrupted three comparison runs during the
    fdc_id re-key, because the migrated store was written after the index and the canonical
    store was not, so runs against the canonical store reused the migrated store's index.
    """
    out: dict[str, list] = {"__canonicalizer__": [_canonicalizer_fingerprint()]}
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
    # A spice is eaten by the gram, not by the 100 g the table reports it on.
    # Cinnamon at 250 mg calcium/100 g outranks milk on a per-100 g comparison
    # and nobody eats 100 g of cinnamon, so under --allow-spices a spice's
    # values are scaled to a 2 g serving against the same 100 g basis. The
    # whole category is excluded by default and no cohort run in the paper
    # enables it; this only decides what happens when someone asks for it.
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
                            args.food_category, args.nutrient_alias)
    sm_path = idx_dir / "static_food_meta.pkl"
    global _STATIC_FOOD_META_PATH
    _STATIC_FOOD_META_PATH = str(sm_path)
    if args.rebuild_static_meta and sm_path.exists(): sm_path.unlink()
    # Rebuild when missing OR older than any input it was derived from. An existence-only
    # check silently serves a stale food universe after a source refresh.
    _sm_sources = (args.food, args.nutrient, args.food_category,
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
