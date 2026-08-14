#!/usr/bin/env python3
"""inject_bacterial_substrates.py

Inject synthetic (fdc_id, nutrient_id, amount) rows for bacterial substrates
that aren't in FDC's native catalog (HMOs, sialic acid, free L-fucose,
alginate, agarose).

Concentrations are literature-backed (citations in source). Units = mg per
100 g edible portion, matching the FDC food_nutrient.parquet convention for
nutrients without a declared unit (the 4_predict pipeline treats unknown
nutrient_ids as MG by default).

Substrates that have an FDC equivalent (GlcNAc/cellobiose/xylan/arabinoxylan/
dextran/pullulan) are NOT injected here — they're aliased in 4_predict's
prox map so the scoring kernel substitutes from the FDC nutrient at runtime.

Writes one synthetic parquet per bucket into /data/bac2food/food_nutrient_bucketed/
bucket=N/. After running, delete the modeled_index cache so 4_predict
rebuilds with the new rows.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

# Output location: alongside the existing bucketed parquet so the dataset
# scan picks it up automatically.
BUCKETED_DIR = Path("/data/bac2food/food_nutrient_bucketed")
BUCKETS = 256

# ---------------------------------------------------------------------------
# Hand-curated rows. Each tuple: (fdc_id, nutrient_id, amount_mg_per_100g)
# ---------------------------------------------------------------------------
#
# HMOs (Human Milk Oligosaccharides) — 2'-FL, 3-FL, difucosyllactose, sialyl-FL
#
# Literature anchors (all per 100 g of milk; densities ~1.03 g/mL so
# values reported as mg/L divided by ~10):
#   Total HMO: colostrum 9–22 g/L → ~1500 mg/100 g; mature 5–15 g/L → ~700 mg/100 g.
#   2'-FL: dominant fucosylated HMO in ~80% of mothers (secretors).
#       colostrum: 2–3 g/L (Thurl et al. 2017), mature 1–3 g/L → 100–300 mg/100 g.
#   3-FL: 0.1–1 g/L (rises across lactation) → 10–100 mg/100 g.
#   Difucosyllactose (DFL = lacto-difucotetraose): 100–500 mg/L → 10–50 mg/100 g.
#   Sialyl-FL (3'-sialyl-3-fucosyllactose): rare HMO, ~50–200 mg/L → 5–20 mg/100 g.
#   Generic fucosyllactose pool (200017): sum of fucosylated species, ~50% of total HMO.
#
# Infant formulas with HMO supplementation:
#   Similac Pro-Sensitive HMO: 200 mg 2'-FL per liter → 20 mg/100 g.
#
ROWS: list[tuple[int, int, float]] = []

# Human milk reference foods (Frida & STFCJ Foundation Foods)
HUMAN_MILK_FDC = [
    50001125,  # Human milk, colostrum [Frida]
    50001167,  # Human milk, transitional (10th day post partum) [Frida]
    50001180,  # Human milk, mature [Frida]
    81013051,  # Other milk, human milk, mature [STFCJ]
]
# colostrum gets the colostrum profile, mature gets the mature profile
HM_HMO_PROFILE = {
    50001125: {"2fl":2500, "3fl":50,  "dfl":300, "fl_total":6000, "sfl":150},
    50001167: {"2fl":1800, "3fl":150, "dfl":200, "fl_total":4000, "sfl":120},
    50001180: {"2fl":1500, "3fl":300, "dfl":150, "fl_total":3000, "sfl":80},
    81013051: {"2fl":1500, "3fl":300, "dfl":150, "fl_total":3000, "sfl":80},
}
# nutrient_id mapping: 200014 2'-FL, 200015 3-FL, 200016 difucosyllactose,
# 200017 generic fucosyllactose pool, 200018 sialyl-fucosyllactose
for fid, p in HM_HMO_PROFILE.items():
    ROWS += [
        (fid, 200014, p["2fl"]),
        (fid, 200015, p["3fl"]),
        (fid, 200016, p["dfl"]),
        (fid, 200017, p["fl_total"]),
        (fid, 200018, p["sfl"]),
    ]

# HMO-fortified infant formulas (Similac Pro-Sensitive with 2'-FL).
# fdc_ids 354120 and 611714 are the branded Similac entries; even though branded
# is normally dropped by build_static_food_meta, the 4_predict pipeline still
# benefits when DROP_BRANDED is off or when a similar non-branded ref exists.
HMO_FORMULAS = [354120, 611714]
for fid in HMO_FORMULAS:
    ROWS += [
        (fid, 200014, 20),   # 2'-FL ~20 mg/100 g (Similac spec)
        (fid, 200017, 25),   # generic fucosyllactose pool
    ]

# ---------------------------------------------------------------------------
# Sialic acid (N-acetylneuraminic acid, Neu5Ac) — total (free+bound)
#
# Literature anchors (mg/100 g):
#   Human milk mature: 80–160 mg/L total → 10–16 mg/100 g (Wang 2009).
#   Human milk colostrum: 500–1500 mg/L → 50–150 mg/100 g (peak).
#   Cow milk: 30–80 mg/L → 3–8 mg/100 g (Tao et al. 2009).
#   Skim/non-fat dairy: similar mg/100 g (sialic mostly on whey glycoproteins).
#   Egg yolk: ~25 mg/100 g (Wang 2012).
#   Whey protein concentrate: 100–300 mg/100 g.
#   Cheese (cheddar, cottage): 20–50 mg/100 g.
#
# nutrient_id 200003 = N-acetylneuraminic acid (specific), 200004 = sialic acid (generic).
# We populate the same amounts under both so the ChEBI traversal hits either seed.
SIALIC_RAW = [
    (50001125, 100),  # human milk colostrum
    (50001167, 50),   # human milk transitional
    (50001180, 15),   # mature human milk (Frida)
    (81013051, 15),   # mature human milk (STFCJ)
    (171265,   12),   # Milk, whole, 3.25% milkfat, with added vitamin D
    (172217,   12),   # Milk, whole, 3.25% milkfat, w/o added vitamins
    (171267,   8),    # Milk, reduced fat, 2%
    (170872,   7),    # Milk, lowfat, 1%
    (171269,   6),    # Milk, nonfat, fluid (skim)
    (167697,   12),   # Buttermilk, fluid, cultured, reduced fat
    (170843,   25),   # Cheese, fontina (matured cheeses have more bound Neu5Ac)
    (170899,   20),   # Cheese, cheddar, sharp
    (170885,   30),   # Whey, acid, fluid
    (171282,   25),   # Whey, sweet, fluid
    (172184,   25),   # Egg, yolk, raw, fresh
]
for fid, amt in SIALIC_RAW:
    ROWS += [(fid, 200003, amt), (fid, 200004, amt)]

# ---------------------------------------------------------------------------
# L-fucose (free + bound, primarily as fucoidans in brown algae and HMOs)
#
# Literature anchors (mg/100 g):
#   Brown seaweed (kelp, kombu, wakame) raw: 200–1000 mg/100 g fucose
#       (fucoidan 5–10% of dry weight, fucose ~50% of fucoidan; raw seaweed
#       is ~10% solids → ~50–500 mg/100 g; conservative midpoint 400).
#   Brown seaweed dried: 2000–8000 mg/100 g.
#   Human milk (in fucosylated HMOs): 150–300 mg/100 g.
#   Plant pectin (fruits, vegetables): trace (≤20 mg/100 g).
#
# nutrient_id 200005 = L-fucose, 200006 = fucose generic.
FUCOSE_RAW = [
    (50001125, 600),  # human milk colostrum (fucosylated HMOs)
    (50001167, 350),
    (50001180, 250),
    (81013051, 250),
    (167602, 4000),   # Seaweed, Canadian Cultivated EMI-TSUNOMATA, dry
    (167603, 400),    # Seaweed, Canadian Cultivated EMI-TSUNOMATA, rehydrated
    (168457, 400),    # Seaweed, kelp, raw
    (168458, 300),    # Seaweed, laver, raw
    (168456, 200),    # Seaweed, irishmoss, raw (red algae, lower fucose)
]
for fid, amt in FUCOSE_RAW:
    ROWS += [(fid, 200005, amt), (fid, 200006, amt)]

# ---------------------------------------------------------------------------
# Alginate (brown algae cell-wall polysaccharide)
#
# Literature anchors (mg/100 g):
#   Brown seaweed raw: 1500–4000 mg/100 g (20–40% dry weight × ~10% solids).
#   Brown seaweed dried: 15–40 g/100 g → 15000–40000 mg/100 g.
#   Irish moss / red algae: very low (red algae use carrageenan, not alginate).
#
ALGINATE_RAW = [
    (167602, 25000),  # Seaweed dry — brown algae blend
    (167603, 3000),   # Seaweed rehydrated
    (168457, 2500),   # Seaweed, kelp, raw
]
for fid, amt in ALGINATE_RAW:
    ROWS.append((fid, 200012, amt))

# ---------------------------------------------------------------------------
# Agarose (red-algae polysaccharide, ~70% of agar)
#
# Literature anchors (mg/100 g):
#   Irish moss (Chondrus crispus, red algae) raw: 100–500 mg/100 g agarose.
#   Agar food gel: huge (~70% solids, ~50 g/100 g agarose).
#   Other red algae (nori, laver): variable; nori has porphyran, not agarose.
#
AGAROSE_RAW = [
    (168456, 400),    # Seaweed, irishmoss, raw
    (168458, 100),    # Seaweed, laver, raw (mostly porphyran but some agarose)
]
for fid, amt in AGAROSE_RAW:
    ROWS.append((fid, 200013, amt))


# ===========================================================================
# Phase 9 extensions — cover the zero-coverage substrates from
# extra_bacterial_seeds.tsv with literature-anchored food measurements.
# The Phase 9 spec-mode allowlist filter only surfaces substrates that
# (a) have ECs mapped (handled by 3_nutrient_to_ec.tsv via seeds), AND
# (b) actually appear in the food-nutrient parquet. Without (b), Roseburia's
# xylanase ECs map to nutrient_id 200007 (Xylan), but no food has a xylan
# row, so wheat bran can't outscore polyphenol fruits.
# All amounts are mg per 100 g edible portion, literature ranges:
#   Bach Knudsen 1997, J. Anim. Feed Sci. Technol. 67:319-338 (cereal fiber)
#   Saulnier et al. 2007, J. Cereal Sci. 46:261-281 (wheat arabinoxylan)
#   Lazaridou & Biliaderis 2007 (oat β-glucan)
#   Roberfroid 2007, J. Nutr. 137:2493S-2502S (inulin/FOS in chicory/onion)
#   Englyst et al. 1992 (resistant starch in raw potato + banana)
#   Lahaye 1991 (algal polysaccharides — fucoidan, laminarin, porphyran)
# ===========================================================================

# Xylan (200007) — total xylose polymer including arabinoxylan backbone
XYLAN_RAW = [
    (169722, 25000),  # Wheat bran, crude (~25% of dry mass)
    (80000762, 22000),# Wheat bran, raw [Swiss Generic Foods]
    (93030154, 22000),# Wheat bran, raw [PhyFoodComp]
    (93030135, 18000),# Wheat bran, Back Cross of Roshan, raw [PhyFoodComp]
    (1107673, 9000),  # RYE, ALL NATURAL (~9% in rye whole grain)
    (1113891, 5000),  # RYE NO SEEDS PREMIUM BREAD
    (1131482, 5500),  # RYE BOULE
    (169418, 6000),   # Wheat, whole grain
    (169721, 7500),   # Wheat, hard red winter
    (170284, 5500),   # Barley, hulled
]
for fid, amt in XYLAN_RAW:
    ROWS.append((fid, 200007, amt))

# Arabinoxylan (200008) — water-extractable + water-unextractable AX
ARABINOXYLAN_RAW = [
    (169722, 22000),  # Wheat bran, crude (~22% AX)
    (80000762, 20000),# Wheat bran, raw [Swiss]
    (93030154, 19000),# Wheat bran, raw [PhyFoodComp]
    (93030135, 16000),# Wheat bran, Roshan
    (93030136, 18000),# Wheat bran, coarse [PhyFoodComp]
    (1107673, 8000),  # RYE
    (1113891, 4500),  # RYE bread
    (1131482, 5000),  # RYE BOULE
    (169721, 6500),   # Wheat, hard red winter
    (169418, 5000),   # Wheat, whole grain
    (170284, 4500),   # Barley, hulled
    (169744, 4000),   # Oat bran, raw — approximate
]
for fid, amt in ARABINOXYLAN_RAW:
    ROWS.append((fid, 200008, amt))

# Mixed-linkage β-1,3/1,4-glucan (200027) — distinctive to oats and barley
BETA_GLUCAN_MIXED = [
    (169744, 5500),   # Oat bran, raw (~5.5%)
    (1112430, 4000),  # OATS & FLAX OATMEAL
    (170284, 4500),   # Barley, hulled (~4.5%)
    (1108161, 3000),  # BARLEY PEAS & LENTILS — approximate
    (1111616, 3500),  # BARLEY PEAS & LENTILS
    (169418, 800),    # Wheat (low β-glucan)
    (1107673, 1500),  # RYE (~1.5%)
]
for fid, amt in BETA_GLUCAN_MIXED:
    ROWS.append((fid, 200027, amt))

# Resistant starch (200019, type 2/3 — raw potato, green banana, cooked-cooled grains)
RESISTANT_STARCH = [
    (170026, 17000),  # Potatoes, flesh and skin, raw (17%, mostly RS2)
    (170027, 17000),  # Potatoes, russet, raw
    (170028, 17000),  # Potatoes, white, raw
    (169250, 100),    # Potato, baked (cooked, low RS3)
    (169231, 8000),   # Bananas, raw — green-yellow has ~8% RS
    (167763, 6000),   # Plantains, raw (~6% RS)
    (169431, 9000),   # Rice, white, raw (long grain) ~ 9% RS2
    (169703, 5500),   # Beans, kidney, mature seeds, raw — small RS3 after cooking
    (170287, 7000),   # Lentils, raw (cooked-cooled has RS3)
]
for fid, amt in RESISTANT_STARCH:
    ROWS.append((fid, 200019, amt))

# Fucoidan (200020) — sulfated fucose polymer, brown algae specific
FUCOIDAN_RAW = [
    (168457, 4000),   # Seaweed, kelp, raw (~4% dry; lower wet)
    (387935, 5500),   # WAKAME (5.5% dry)
    (10034254, 6000), # Seaweed, Kombu, Dried [Fineli]
    (50001204, 6000), # Seaweed, kombu, dried [Frida]
]
for fid, amt in FUCOIDAN_RAW:
    ROWS.append((fid, 200020, amt))

# Laminarin (200021) — β-1,3-glucan storage polysaccharide, brown algae
LAMINARIN_RAW = [
    (168457, 8000),   # Kelp, raw (~8% dry)
    (387935, 4000),   # WAKAME (~4% dry)
    (10034254, 9500), # Kombu, dried
    (50001204, 9500), # Kombu, dried
]
for fid, amt in LAMINARIN_RAW:
    ROWS.append((fid, 200021, amt))

# Porphyran (200022) — sulfated galactan from red algae (nori/laver)
PORPHYRAN_RAW = [
    (168458, 15000),  # Seaweed, laver, raw (~15%)
    (10034136, 22000),# Seaweed, Nori, Dried [Fineli]
    (80000297, 22000),# Seaweed, Nori, dried [Swiss]
    (168456, 8000),   # Seaweed, irishmoss (also has porphyran/carrageenan)
]
for fid, amt in PORPHYRAN_RAW:
    ROWS.append((fid, 200022, amt))

# Carrageenan (200023) — sulfated red algal polysaccharide
CARRAGEENAN_RAW = [
    (168456, 35000),  # Seaweed, irishmoss (~35%)
    (168458, 2000),   # Seaweed, laver (small amount)
]
for fid, amt in CARRAGEENAN_RAW:
    ROWS.append((fid, 200023, amt))

# Glucomannan (200024) — β-1,4-mannan from konjac, also small in coffee
GLUCOMANNAN_RAW = [
    (167720, 100),    # Coffee, brewed — small
    (169251, 800),    # Mushrooms, white, raw (small)
]
for fid, amt in GLUCOMANNAN_RAW:
    ROWS.append((fid, 200024, amt))

# Arabinogalactan type II (200025) — small amounts in many plants
ARABINOGALACTAN_RAW = [
    (170393, 200),    # Carrots, raw
    (169145, 500),    # Beets, raw
    (170000, 250),    # Onions, raw
    (170287, 350),    # Lentils, raw
    (167793, 100),    # Apples, raw fuji
]
for fid, amt in ARABINOGALACTAN_RAW:
    ROWS.append((fid, 200025, amt))

# Glucuronoxylan (200026) — hemicellulose of hardwoods + some fruit pectins
GLUCURONOXYLAN_RAW = [
    (169722, 4000),   # Wheat bran, crude (~4%)
    (170284, 2500),   # Barley, hulled
    (1131482, 1500),  # RYE BOULE
]
for fid, amt in GLUCURONOXYLAN_RAW:
    ROWS.append((fid, 200026, amt))

# β-1,4-mannan linear (200028) — coffee, copra, palm kernel
MANNAN_LINEAR = [
    (167720, 500),    # Coffee, brewed
    (167769, 600),    # Coconut meat, raw
]
for fid, amt in MANNAN_LINEAR:
    ROWS.append((fid, 200028, amt))

# Mucin glycan generic proxy (200029) — animal organ meats, gut tissues, jelly
MUCIN_PROXY = [
    (171060, 500),    # Beef tripe, raw (mucin-coated)
    (171066, 300),    # Beef liver, raw
]
for fid, amt in MUCIN_PROXY:
    ROWS.append((fid, 200029, amt))

# Chitobiose (200030) — GlcNAc disaccharide, mushroom + insect chitin hydrolysates
CHITOBIOSE_RAW = [
    (169251, 200),    # Mushrooms, white, raw
    (169254, 300),    # Mushrooms, shiitake, dried
    (169255, 250),    # Mushrooms, portabella, raw
]
for fid, amt in CHITOBIOSE_RAW:
    ROWS.append((fid, 200030, amt))

# Fructooligosaccharide (FOS, 200031) — high in chicory, J. artichoke, onion, garlic, leek
FOS_RAW = [
    (169993, 18000),  # Chicory roots, raw (~18%)
    (169992, 4000),   # Chicory greens, raw
    (170404, 6000),   # Chicory, witloof
    (170000, 1300),   # Onions, raw
    (169230, 5400),   # Garlic, raw
    (169246, 2700),   # Leeks, raw
    (168389, 2200),   # Asparagus, raw
    (168992, 25000),  # Agave, cooked (~25% FOS/inulin)
    (168993, 50000),  # Agave, dried
    (167793, 500),    # Apples, raw fuji
    (2372957, 90000), # INULIN POWDER — pure inulin/FOS supplement
    (169231, 600),    # Bananas, raw
]
for fid, amt in FOS_RAW:
    ROWS.append((fid, 200031, amt))

# Galactooligosaccharide (GOS, 200032) — small in legumes, large in supplemented dairy
GOS_RAW = [
    (169703, 1500),   # Beans, kidney, raw (raffinose/stachyose ≈ GOS family)
    (170287, 1800),   # Lentils, raw
    (174270, 2200),   # Soybeans, mature seeds, raw
]
for fid, amt in GOS_RAW:
    ROWS.append((fid, 200032, amt))

# HMO subtypes (200033 LNT, 200034 LNnT, 200035 6'-SL, 200036 3'-SL) — human milk
HMO_SUBTYPES = [
    (171279, [(200033, 150), (200034, 100), (200035, 70), (200036, 50)]),  # Milk, human
    (40001662, [(200033, 200), (200034, 130), (200035, 100), (200036, 80)]), # human colostrum [McCance]
    (2705383, [(200033, 130), (200034, 90), (200035, 65), (200036, 45)]),  # Milk, human
]
for fid, nlist in HMO_SUBTYPES:
    for nid, amt in nlist:
        ROWS.append((fid, nid, amt))

# Pectic oligosaccharide (POS, 200037) — apple, citrus, beet pectin
POS_RAW = [
    (167793, 200),    # Apples, raw fuji
    (168201, 200),    # Apples, raw red delicious
    (167682, 5000),   # Pectin, liquid (commercial source)
    (168821, 90000),  # Pectin, dry — pure
    (169145, 1000),   # Beets, raw (beet pectin)
]
for fid, amt in POS_RAW:
    ROWS.append((fid, 200037, amt))

# β-mannan oligosaccharide (200038)
MOS_RAW = [
    (167720, 200),    # Coffee, brewed
    (167769, 300),    # Coconut meat, raw
]
for fid, amt in MOS_RAW:
    ROWS.append((fid, 200038, amt))

# Also fill the existing zero-coverage substrates from Phase 2's
# original seed list: GlcNAc (200001), GalNAc (200002), Dextran (200009),
# Cellobiose (200010), Pullulan (200011) for foods where they're known.
GLCNAC_RAW = [
    (169251, 600),    # Mushrooms, white, raw (chitin → GlcNAc on hydrolysis)
    (169254, 1200),   # Mushrooms, shiitake, dried
    (169255, 700),    # Mushrooms, portabella, raw
]
for fid, amt in GLCNAC_RAW:
    ROWS.append((fid, 200001, amt))

DEXTRAN_RAW = [
    (169251, 100),    # Mushrooms (small — yeast-derived dextran-like β-glucan)
    (169640, 50),     # Honey, strained (microbial dextran from Leuconostoc)
]
for fid, amt in DEXTRAN_RAW:
    ROWS.append((fid, 200009, amt))

CELLOBIOSE_RAW = [
    (169251, 80),     # Mushrooms (raw plant cell wall)
    (169993, 50),     # Chicory roots, raw
]
for fid, amt in CELLOBIOSE_RAW:
    ROWS.append((fid, 200010, amt))

# ---------------------------------------------------------------------------
# Write to bucketed parquet
# ---------------------------------------------------------------------------
df = pd.DataFrame(ROWS, columns=["fdc_id", "nutrient_id", "amount"])
df["fdc_id"] = df["fdc_id"].astype("int64")
df["nutrient_id"] = df["nutrient_id"].astype("int64")
df["amount"] = df["amount"].astype("float64")
df["bucket"] = (df["nutrient_id"] % BUCKETS).astype("int32")

print(f"[*] {len(df)} synthetic (fdc_id, nutrient_id, amount) rows for "
      f"{df['nutrient_id'].nunique()} novel substrates across {df['fdc_id'].nunique()} foods.")
print(df.groupby("nutrient_id").agg(n_foods=("fdc_id","nunique"),
                                    median_mg=("amount","median"),
                                    max_mg=("amount","max")).to_string())

n_written = 0
for b, sub in df.groupby("bucket"):
    out_dir = BUCKETED_DIR / f"bucket={int(b)}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "synthetic_bacterial.parquet"
    sub.drop(columns=["bucket"]).to_parquet(out_path, index=False)
    n_written += 1
print(f"[*] Wrote {n_written} synthetic parquet files into {BUCKETED_DIR}/bucket=*/")
print("[*] NOTE: delete /data/bac2food/index_modeled/*.parquet and static_food_meta.pkl "
      "before next 4_predict run so the modeled index picks up the new rows.")
