#!/usr/bin/env python3
"""
Reproducibility probe: can a third party rebuild a bac2food-style predictor
using ONLY the three exported TSVs? No repo code, no private tables.

Chain attempted:  organism -> ec_number -> nutrient_ids -> fdc_id (ranked foods)
"""
import csv, sys, collections

EXP = "/data/bac2food/exports"
ORG_SUBSTR = sys.argv[1] if len(sys.argv) > 1 else "Bacteroides thetaiotaomicron"

csv.field_size_limit(10**9)

# ---- STEP 1: organism -> ECs  (species_enzymes.tsv only) -------------------
ecs = set()
orgs_hit = set()
with open(f"{EXP}/species_enzymes.tsv", newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if ORG_SUBSTR.lower() in (row["organism"] or "").lower():
            ecs.add(row["ec_number"])
            orgs_hit.add(row["organism"])
print(f"STEP 1  organism->EC : {len(ecs):,} distinct EC across {len(orgs_hit):,} matching organism strings")

# ---- STEP 2: EC -> nutrient_ids  (enzyme_substrate_chebi.tsv only) ---------
nut_ids = set()
nut_names = {}
ec_with_nut = set()
edges = 0
with open(f"{EXP}/enzyme_substrate_chebi.tsv", newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        if row["ec_number"] not in ecs:
            continue
        if row["in_model"] != "yes" or not row["nutrient_ids"]:
            continue
        ids = [x.strip() for x in row["nutrient_ids"].split(",") if x.strip()]
        names = [x.strip() for x in (row["nutrient_names"] or "").split(";")]
        for i, nid in enumerate(ids):
            nut_ids.add(nid)
            if i < len(names) and names[i]:
                nut_names.setdefault(nid, names[i])
            ec_with_nut.add(row["ec_number"])
            edges += 1
print(f"STEP 2  EC->nutrient : {len(nut_ids):,} distinct nutrient_id reachable "
      f"({len(ec_with_nut):,}/{len(ecs):,} of the organism's ECs reach a food nutrient; {edges:,} edges)")

# ---- STEP 3: nutrient -> foods  (food_nutrients.tsv only, streamed) --------
# Naive score: sum over matched nutrients of the MEAN amount (mean = the export's
# own advice for reconciling the retained conflicts). Analytical foods only.
KEEP = {"foundation_food", "sr_legacy_food", "survey_fndds_food"}
acc = collections.defaultdict(lambda: collections.defaultdict(list))  # fdc -> nut -> [amounts]
desc = {}
seen_rows = 0
with open(f"{EXP}/food_nutrients.tsv", newline="") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        seen_rows += 1
        if row["nutrient_id"] not in nut_ids:
            continue
        if row["data_type"] not in KEEP:
            continue
        try:
            amt = float(row["amount"])
        except (TypeError, ValueError):
            continue
        if amt <= 0:
            continue
        fid = row["fdc_id"]
        acc[fid][row["nutrient_id"]].append(amt)
        desc[fid] = row["description"]
print(f"STEP 3  nutrient->food: streamed {seen_rows:,} rows -> {len(acc):,} analytical foods touched")

# rank: breadth (how many of the organism's nutrients the food supplies) then mass
ranked = []
for fid, nuts in acc.items():
    breadth = len(nuts)
    mass = sum(sum(v) / len(v) for v in nuts.values())   # MEAN per nutrient, then sum
    ranked.append((breadth, mass, fid))
ranked.sort(reverse=True)

print(f"\n=== naive top-15 foods for '{ORG_SUBSTR}' (built from 3 TSVs alone) ===")
for breadth, mass, fid in ranked[:15]:
    print(f"  {breadth:2d} nutrients | {mass:12.2f} | {desc[fid][:62]}")

print(f"\n--- conflict evidence is visible to a third party too ---")
multi = sum(1 for fid in acc for n in acc[fid] if len(acc[fid][n]) > 1)
print(f"  (food,nutrient) cells with >1 retained value in this slice: {multi:,}")
