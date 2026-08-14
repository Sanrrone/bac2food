#!/usr/bin/env python3
"""ec2food_differential.py — Ecosystem-Aware food recommendations"""
from __future__ import annotations
import argparse, math, os, re, gc, pickle, heapq, zlib, shutil
from pathlib import Path
import pandas as pd
import numpy as np
import pyarrow.dataset as ds
try: import resource
except ImportError: resource = None

#usage: python ec2food_differential.py   --enzyme_tsv ../sample_A.tsv   --out_prefix my_probiotic_blend   --promote_species "1134687_Klebsiella_michiganensis,137838_Clostridium_neonatale"   --inhibit_species "1311_Streptococcus_agalactiae,1309_Streptococcus_mutans"   --inhibit_weight 2.0   --auto_tune
# ==============================================================================
# PIPELINE CONFIGURATION
# ==============================================================================
PATH_NUTRIENT_TO_EC     = "../0_building/3_nutrient_to_ec.tsv"
PATH_NUTRIENT_CSV       = "/data/bac2food/nutrient.csv"
PATH_FOOD_CSV           = "/data/bac2food/food.parquet"
PATH_FOOD_CATEGORY_CSV  = "/data/bac2food/food_category.csv"
PATH_BUCKETED_DIR       = "/data/bac2food/food_nutrient_bucketed"
PATH_INDEX_DIR          = "/data/bac2food/index_modeled"

DROP_BRANDED = True
ALLOW_MACRO_PROXY = False
ALLOW_MACRO_SCAN = False

EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
_FIBER = re.compile(r"\b(?:fiber|starch|oligosaccharide|inulin|pectin|glucan|hemicellulose|cellulose|pentosan)\b", re.I)

BUCKETS = 256
TAU_Q, MIN_CONTRIB_FRAC, GAIN_MIN, SCORE_MIN = 0.05, 0.01, 0.001, -2.0
PROC_W, BROAD_W, BROAD_Q, TYPE_W, ART_W = 0.20, 0.20, 2.0, 0.25, 8.0
K_PLANT, K_OTHER, MAX_CANONS = 800, 200, 250000
PROC_BASE_W, JUNK_BASE_W, ALC_BASE_W = 0.15, 0.20, 0.30
LAMBDA_FOOD, OVERLAP_W, FIBER_WEIGHT = 0.08, 0.60, 2.0

CAT_PENALTY = {
    "Fast Foods":5.0,"Restaurant Foods":5.0,"Sweets":5.0,"Snacks":4.0,"Alcoholic Beverages":5.0,
    "Baked Products":3.0,"Beverages":4.0,"Spices and Herbs":5.0,"Fats and Oils":4.0,
    "Fruits and Fruit Juices":0.05,"Cereal Grains and Pasta":0.30,"Nut and Seed Products":0.10,
    "Vegetables and Vegetable Products":0.0,"Legumes and Legume Products":0.0,
}
PLANT_CATS = {"Vegetables and Vegetable Products","Legumes and Legume Products",
              "Fruits and Fruit Juices","Cereal Grains and Pasta","Nut and Seed Products"}
TP_MAP = {"foundation_food":600,"sr_legacy_food":500,"survey_fndds_food":450,
          "experimental_food":300,"branded_food":100}

# The Master Blacklist
_BLOCKED = re.compile(
    r"\b(?:sucrose|glucose|fructose|lactose|maltose|galactose|mannitol|xylitol|sorbitol|sugar|alanine|arginine|asparagine|aspartic acid|cysteine|glutamic acid|glutamine|glycine|histidine|isoleucine|leucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|cystine|cholesterol|sfa|mufa|pufa|tfa|trans|fatty acid|fatty acids|sterol|dha|epa|docosahexaenoic|eicosapentaenoic|polyunsaturated|monounsaturated|saturated|triglyceride|lipid|lipids|phospholipids?|glycolipids?|unsaponifiable|trans-monoenoic|trans-dienoic|trans-polyenoic|vitamin\s*[abcdefjkh][-\s\d]*|vitamin\s*b[-\s]?(?:[1-9]|12)|vitamins\b|b[-\s]?2\b|b[-\s]?6\b|b[-\s]?12\b|riboflavin|thiamin|niacin|folate|folic|tetrahydrofolate|5-mthf|pantothenic|biotin|pyridoxin|cobalamin|retinol|retinoic|carotene|tocopherols?|tocotrienols?|total\s+tocopherols?|ascorbic|menaquinone|phylloquinone|calcium|iron|magnesium|phosphorus|potassium|sodium|zinc|copper|selenium|manganese|iodine|choline|betaine|chromium|fluoride|molybdenum|c\d{1,2}:\d+|nacl|salt|salt deklaration|salt labelling|niacin[æa]kvivalent|niacine, pr[ée]form[ée]e|polyols|tannins?|flavans?|flavones?|flavonols?|anthocyanins?|haem|nhaem|val|thr|tyr|glu|gly|his|pro|ser|arg|asp|phe|ala|met|cys|ile|leu|trp|valin|threonin|tyrosin|glycin|prolin|serin|arginin|phenylalanin|alanin|methionin|cystein|isoleucin|leucin|hydroxyproline|glycogen|tyramin|putrescin|chitin|polyphenol|nutrient_name)\b|(?:aa \(pernitrogen\)|starch-|pp-|polymers \(>10 mers\))", 
    re.I
)

load_smart = lambda p: pd.read_parquet(str(p)) if str(p).lower().endswith((".parquet",".pq")) else pd.read_csv(str(p),low_memory=False)
def cosine_sim(a,b):
    if not a or not b: return 0.0
    i=set(a)&set(b)
    if not i: return 0.0
    dot=sum(a[k]*b[k] for k in i); na=math.sqrt(sum(v*v for v in a.values())); nb=math.sqrt(sum(v*v for v in b.values()))
    return dot/(na*nb) if na>0 and nb>0 else 0.0

def round_floats(df, nd=4):
    if df is None or df.empty: return df
    c = df.select_dtypes(include=["float16", "float32", "float64"]).columns
    if len(c): df[c] = df[c].round(nd)
    return df

class TopKByNutrient:
    __slots__=("val","heap")
    def __init__(self): self.val={}; self.heap={}
    def push(self,nid,canon,amt,K):
        d,h=self.val.setdefault(nid,{}),self.heap.setdefault(nid,[])
        if d.get(canon,-1)>=amt: return
        d[canon]=amt; heapq.heappush(h,(amt,canon))
        if len(d)>K*2: newh=[(a,c) for a,c in h if d.get(c)==a]; heapq.heapify(newh); self.heap[nid]=newh
        while len(d)>K and self.heap[nid]:
            a,c=heapq.heappop(self.heap[nid])
            if d.get(c)==a: del d[c]

def _run_differential_diet(STATIC_DB, DYNAMIC_STATE, promote_species, inhibit_species, inh_weight, coll_weight, maxf):
    if inhibit_species: print(f"    - Inhibiting: {', '.join(inhibit_species)} (Penalty Weight: {inh_weight}x)")

    targs_prom, targs_inh, targs_all = set(), set(), set()
    bystander_counts = {}
    
    total_bystanders = len(DYNAMIC_STATE["lbl2targ"]) - len(promote_species) - len(inhibit_species)
    total_bystanders = max(1, total_bystanders)
    
    # 1. Target aggregation and Bystander Frequency mapping
    for lbl, targets in DYNAMIC_STATE["lbl2targ"].items():
        targs_all |= targets
        if lbl in promote_species: 
            targs_prom |= targets
        elif lbl in inhibit_species: 
            targs_inh |= targets
        else:
            for n in targets:
                bystander_counts[n] = bystander_counts.get(n, 0) + 1

    # 2. Filter against modeled nutrients and blocklist
    targs_prom = (targs_prom & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
    targs_inh = (targs_inh & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
    need = (targs_all & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
    
    # 3. Calculate EXACTLY what fraction of the untargeted gut eats each nutrient
    bystander_frac = {n: bystander_counts.get(n, 0) / total_bystanders for n in need}
    
    species_targets = {sp: (DYNAMIC_STATE["lbl2targ"].get(sp, set()) & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"] for sp in promote_species}
    
    if not targs_prom: return pd.DataFrame()

    def _tau(n): return min(max(0.1 if STATIC_DB["unit_arr"][n]==1.0 else 0.01,(0.01 if DYNAMIC_STATE["ndf"].get(n,1e9)<50000 else TAU_Q)*DYNAMIC_STATE["ref"].get(n,0.0)),500.0 if STATIC_DB["unit_arr"][n]==1.0 else 50.0)
    tau = {n: _tau(n) for n in need}
    
    scn = ds.dataset(DYNAMIC_STATE["bdir"], format="parquet", partitioning="hive").scanner(
        columns=["fdc_id", "nutrient_id", "amount"],
        filter=ds.field("bucket").isin(sorted({n % BUCKETS for n in need})) & ds.field("nutrient_id").isin(sorted(need)),
        batch_size=100_000
    )

    tk_p, tk_o = TopKByNutrient(), TopKByNutrient()
    c_arr, p_arr, pl_arr, sp_arr, u_arr = STATIC_DB["canon_arr"], STATIC_DB["port_arr"], STATIC_DB["plant_arr"], STATIC_DB["spice_arr"], STATIC_DB["unit_arr"]
    max_f, max_u = len(c_arr), len(u_arr)
    
    for b in scn.to_batches():
        f, n, a = b["fdc_id"].to_numpy().astype(np.int64), b["nutrient_id"].to_numpy().astype(np.int32), b["amount"].to_numpy().astype(np.float32)
        m = (f >= 0) & (f < max_f) & (n >= 0) & (n < max_u) & np.isfinite(a)
        if not m.any(): continue
        f, n, a = f[m], n[m], a[m]; c = c_arr[f]; v = (c != -1) & (~sp_arr[f])
        if not v.any(): continue
        f, n, a, c = f[m], n[m], a[m], c[m]
        a = a * (p_arr[f] / 100.0) * u_arr[n]
        for i in range(len(n)): (tk_p if pl_arr[f[i]] else tk_o).push(int(n[i]), int(c[i]), float(a[i]), K_PLANT if pl_arr[f[i]] else K_OTHER)

    def get_tk(nid):
        d = dict(tk_p.val.get(nid, {}))
        for k, v in tk_o.val.get(nid, {}).items():
            if k not in d or v > d[k]: d[k] = v
        return d

    c2a = {}
    for t in need:
        th, nn = float(tau.get(t, 0.0)), STATIC_DB["nut_name"].get(t, "")
        mc = 0.0001 if _FIBER.search(nn) else (0.001 if DYNAMIC_STATE["ndf"].get(t, 1e9) < 50000 else MIN_CONTRIB_FRAC)
        base_cut = th * mc
        for c, a in get_tk(int(t)).items():
            if a > base_cut: c2a.setdefault(c, {})[int(t)] = float(a)
            
    c2a = {c: ct for c, ct in c2a.items() if any(n in targs_prom for n in ct)}

    w = {n: (FIBER_WEIGHT if _FIBER.search(STATIC_DB["nut_name"].get(n, "")) else 1.0) * (1.0 + math.log(DYNAMIC_STATE["N"] / max(1, int(DYNAMIC_STATE["ndf"].get(n, 1))))) for n in need}
    
    ST = type('S', (), {'cur': {n: 0.0 for n in need}, 'chosen': set(), 'seen_desc': set(), 'chosen_imps': [], 'rows': [], 'chosen_cats': []})()
    fs_db, nd = STATIC_DB["food_stats"], STATIC_DB["nut_name"]

    for rk in range(1, maxf + 1):
        bst = None
        for c, ct in c2a.items():
            fs = fs_db.get(c, {})
            if c in ST.chosen or not fs or fs["dh"] in ST.seen_desc: continue
            
            # THE ROBUST CATEGORY CAP: Only 1 item per category!
            # THE FLEXIBLE CATEGORY CAP: Only 1 item per standard category!
            cat_str = fs.get("cat", "")
            # If it has a real category, enforce the cap of 1. If it's blank, let it pass!
            if cat_str and not pd.isna(cat_str) and str(cat_str).strip() != "":
                if ST.chosen_cats.count(cat_str) >= 1: 
                    continue
            
            gn_prom, gn_inh, gn_col = 0.0, 0.0, 0.0
            ni_prom = 0
            imp = {}
            gn_prom_per_sp = {sp: 0.0 for sp in promote_species}
            
            for n, a in ct.items():
                th = tau[n]
                dg = math.log1p((ST.cur[n] + a) / th) - math.log1p(ST.cur[n] / th)
                if dg > 1e-12:
                    cg = w.get(n, 1.0) * dg
                    
                    if n in targs_prom:
                        gn_prom += cg
                        ni_prom += 1
                        imp[n] = cg
                        for sp in promote_species:
                            if n in species_targets[sp]: gn_prom_per_sp[sp] += cg
                            
                    if n in targs_inh: 
                        gn_inh += cg
                        
                    # THE TRUE FIX: Always calculate collateral penalty!
                    # bystander_frac ALREADY perfectly tracks how many non-target bugs eat this nutrient.
                    gn_col += cg * bystander_frac.get(n, 0.0)
                        
                    # THE NEW DIFFERENTIAL FORMULA
                        
                    # NO MORE COMMENSAL PASS. All collateral damage is tracked.
                    if n in (need - targs_prom - targs_inh):
                        gn_col += cg * bystander_frac.get(n, 0.0)
                        
            # THE NEW DIFFERENTIAL FORMULA
            net_gain = gn_prom - (inh_weight * gn_inh) - (coll_weight * gn_col)
            
            if ni_prom < 2 or net_gain < GAIN_MIN: continue
            
            imp = {n: v for n, v in imp.items() if v >= 0.05 * gn_prom}
            red = sum(sorted([cosine_sim(imp, p) for p in ST.chosen_imps], reverse=True)[:3]) / 3.0 if ST.chosen_imps else 0.0
            prm = min(1.0, max(0.0, ni_prom / max(50, int(DYNAMIC_STATE["mc"].get(c, 0)))))
            geff = math.log1p(net_gain) * (math.sqrt(ni_prom / (ni_prom + 2.0)) if ni_prom > 0 else 0.0)
            
            cat_rep = ST.chosen_cats.count(fs.get("cat", ""))
            cat_pen = 0.20 * cat_rep
            proc_cost = PROC_W * fs["pl"] * math.log1p(gn_prom)
            broad_cost = BROAD_W * math.log1p(gn_prom) * ((1.0 - prm)**BROAD_Q)
            
            sc = geff - LAMBDA_FOOD - proc_cost - broad_cost - (ART_W * fs["art"]) - (OVERLAP_W * red) - fs["bc"] + (fs["tb"] * geff) - cat_pen
            
            if sc > SCORE_MIN and (bst is None or sc > bst[0]):
                top = sorted(imp.items(), key=lambda x: x[1], reverse=True)[:3]
                bst = (sc, net_gain, gn_prom, gn_inh, gn_col, c, ni_prom, fs["pl"], geff, red, imp,
                       ",".join(str(n) for n, _ in top), ",".join(nd.get(n, str(n)) for n, _ in top), fs["dh"], gn_prom_per_sp, cat_str)
                       
        if not bst: break
        c = bst[5]; ST.chosen.add(c); ST.seen_desc.add(bst[13]); ST.chosen_imps.append(bst[10])
        # Track the category string (including our fix for blank categories)
        ST.chosen_cats.append(bst[15]) 
        
        for n, a in c2a[c].items(): ST.cur[n] += a
        
        indiv_gains_str = " | ".join([f"{sp.split('_', 1)[-1] if '_' in sp else sp}: {round(bst[14][sp], 4)}" for sp in promote_species])
        
        ST.rows.append(dict(
            rank=rk,
            intervention=", ".join(promote_species),
            description="",
            top_beneficial_nutrients=bst[12],
            differential_score=bst[0],
            net_gain=bst[1],
            gross_promote_gain=bst[2],
            individual_gains=indiv_gains_str,
            gross_inhibit_penalty=bst[3],
            gross_collateral_penalty=bst[4],
            canonical_fdc_id=c,
            representative_fdc_id=c,
            n_improved_promote=bst[6],
            redundancy=bst[9]
        ))
        
    return pd.DataFrame(ST.rows)

def _generate_collateral_report(STATIC_DB, DYNAMIC_STATE, chosen_fdcs, out_path):
    print("\n[*] Running Post-Diet Collateral Evaluation for entire community...")
    
    scn = ds.dataset(DYNAMIC_STATE["bdir"], format="parquet", partitioning="hive").scanner(
        columns=["fdc_id", "nutrient_id", "amount"],
        filter=ds.field("fdc_id").isin(chosen_fdcs)
    )
    
    diet_nutrients = {}
    p_arr, u_arr = STATIC_DB["port_arr"], STATIC_DB["unit_arr"]
    
    for b in scn.to_batches():
        df = b.to_pandas()
        for _, row in df.iterrows():
            f, n, a = int(row['fdc_id']), int(row['nutrient_id']), float(row['amount'])
            if n >= len(u_arr) or f >= len(p_arr): continue
            scaled_amount = a * (p_arr[f] / 100.0) * u_arr[n]
            if scaled_amount > 0: diet_nutrients[n] = diet_nutrients.get(n, 0.0) + scaled_amount

    results = []
    dyn_N, dyn_ndf, dyn_ref = DYNAMIC_STATE["N"], DYNAMIC_STATE["ndf"], DYNAMIC_STATE["ref"]
    
    for species, targets in DYNAMIC_STATE["lbl2targ"].items():
        clean_targs = (targets & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
        if not clean_targs: continue
        
        species_gain, fed_targets = 0.0, 0
        for n in clean_targs:
            a = diet_nutrients.get(n, 0.0)
            if a > 0:
                th = min(max(0.1 if u_arr[n]==1.0 else 0.01,(0.01 if dyn_ndf.get(n,1e9)<50000 else TAU_Q)*dyn_ref.get(n,0.0)),500.0 if u_arr[n]==1.0 else 50.0)
                nn = STATIC_DB["nut_name"].get(n, "")
                w = (FIBER_WEIGHT if _FIBER.search(nn) else 1.0) * (1.0 + math.log(dyn_N / max(1, int(dyn_ndf.get(n, 1)))))
                species_gain += (w * math.log1p(a / th))
                fed_targets += 1
                
        results.append({"species": species, "total_targets": len(clean_targs), "fed_targets": fed_targets, "collateral_gain": round(species_gain, 4)})
        
    res_df = pd.DataFrame(results).sort_values(by="collateral_gain", ascending=False)
    res_df.to_csv(out_path, sep="\t", index=False)
    print(f"[+] Output saved: {out_path}")
    print("\n--- TOP 5 FED SPECIES (Community Check) ---")
    print(res_df.head(5).to_string(index=False))

def main():
    print("==========================================================")
    print("RUNNING ec2food_differential.py (Ecosystem-Aware Therapeutics)")
    print("==========================================================")
    ap = argparse.ArgumentParser()
    ap.add_argument("--enzyme_tsv", required=True, help="Full TSV of species and EC numbers")
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--promote_species", required=True, help="Comma-separated list of species to feed")
    ap.add_argument("--inhibit_species", required=False, default="", help="Comma-separated list of species to starve")
    ap.add_argument("--inhibit_weight", type=float, default=1.0, help="Penalty multiplier for feeding pathogens")
    ap.add_argument("--collateral_weight", type=float, default=0.2, help="Manual penalty multiplier for bystander feeding")
    ap.add_argument("--auto_tune", action="store_true", help="Dynamically loop to find the optimal collateral weight")
    ap.add_argument("--max_foods", type=int, default=5)
    args = ap.parse_args()
    
    sm_path = Path(PATH_INDEX_DIR) / "static_food_meta.pkl"
    nut_df = load_smart(PATH_NUTRIENT_CSV)
    
    # Load Auto-Aliaser Map
    name_to_id = {str(r.name).strip().lower(): int(r.id) for r in nut_df.itertuples() if pd.notna(r.name) and pd.notna(r.id)}
    prox = {}
    auto_alias_map = {
        "tyramin": "tyramine", "putrescin": "putrescine", "cholecalciferol (vit d3)": "cholecalciferol",
        "ferrulic acid": "ferulic acid", "zink": "zinc", "sulfur, s": "sulfur", "mercury, hg": "mercury",
        "rubidium, rb": "rubidium", "titanium, ti": "titanium", "5-formyltetrahydrofolic acid (5-hcoh4": "5-formyltetrahydrofolic acid",
        "tetrahydrofolic acid (thf)": "tetrahydrofolic acid", "(+)-catechin": "catechin", "(-)-epicatechin": "epicatechin",
        "(-)-epigallocatechin": "epigallocatechin", "(-)-epicatechin 3-o-gallate": "epicatechin gallate",
        "epicatechin-3-gallate": "epicatechin gallate", "p-hydroxy benzoic acid": "4-hydroxybenzoic acid",
        "mk4": "menaquinone-4", "mk5": "menaquinone-5", "mk6": "menaquinone-6", "mk7": "menaquinone-7",
        "mk8": "menaquinone-8", "mk9": "menaquinone-9", "mk10": "menaquinone-10", "ip6": "phytic acid",
        "aconitate": "aconitic acid", "benzoate": "benzoic acid", "cinnamate": "cinnamic acid", "citrate": "citric acid",
        "formate": "formic acid", "fumarate": "fumaric acid", "gentisate": "gentisic acid", "gluconate": "gluconic acid",
        "glycolate": "glycolic acid", "isocitrate": "isocitric acid", "lactate": "lactic acid", "malonate": "malonic acid",
        "orotate": "orotic acid", "oxalate": "oxalic acid", "oxaloacetate": "oxaloacetic acid", "pantothenate": "pantothenic acid",
        "phytate": "phytic acid", "pyruvate": "pyruvic acid", "quinate": "quinic acid", "salicylate": "salicylic acid",
        "succinate": "succinic acid", "tartrate": "tartaric acid", "vanillate": "vanillic acid", "alpha-ketoglutarate": "alpha-ketoglutaric acid"
    }
    for alias, canon_name in auto_alias_map.items():
        if alias in name_to_id and canon_name in name_to_id:
            prox.setdefault(name_to_id[alias], []).append((name_to_id[canon_name], 1.0, True, "form"))

    blocked = set(nut_df[nut_df["name"].astype(str).str.contains(_BLOCKED, regex=True, na=False)]["id"].dropna().astype(int).tolist())

    enz = pd.read_csv(args.enzyme_tsv, sep="\t", dtype="string")
    cols = {c.lower():c for c in enz.columns}
    if "ec_number" not in enz.columns: enz = enz.rename(columns={cols.get("ec_number", cols.get("ec", enz.columns[-1])): "ec_number"})
    if "species" not in enz.columns: enz = enz.rename(columns={enz.columns[0]: "species"})
    enz["strain"] = enz.get("strain", pd.Series([""]*len(enz))).fillna("").astype("string").str.strip()
    enz = enz[enz["ec_number"].str.match(EC_RE, na=False)].dropna(subset=["species", "ec_number"])
    enz_nut = enz.merge(pd.read_csv(PATH_NUTRIENT_TO_EC, sep="\t"), on="ec_number")
    
    enz_nut["label"] = enz_nut["species"] + enz_nut["strain"].apply(lambda x: f" | {x}" if x else "")
    lbl2targ = {lbl: set(sub["nutrient_id"].astype(int).unique()) for lbl, sub in enz_nut.groupby("label")}
    
    tot_pq = pd.read_parquet(Path(PATH_INDEX_DIR) / "modeled_totals.parquet")
    dyn = {
        "bdir": PATH_BUCKETED_DIR, "lbl2targ": lbl2targ, "prox": prox, "blocked": blocked,
        "mc": tot_pq.set_index("canon")["model_count"].to_dict(),
        "ndf": pd.read_parquet(Path(PATH_INDEX_DIR) / "nutrient_df.parquet").set_index("nutrient_id")["df_foods"].to_dict(),
        "N": len(tot_pq),
        "ref": pd.read_parquet(Path(PATH_INDEX_DIR) / "nutrient_scale.parquet").set_index("nutrient_id")["ref_amount_norm"].to_dict(),
    }
    
    promote_list = [s.strip() for s in args.promote_species.split(",") if s.strip()]
    inhibit_list = [s.strip() for s in args.inhibit_species.split(",") if s.strip()] if args.inhibit_species else []
    
    STATIC_DB = pickle.load(open(sm_path, "rb"))
    
    # ---------------------------------------------------------
    # THE DYNAMIC AUTO-TUNER (Hill Climbing with Patience)
    # ---------------------------------------------------------
    best_df = None
    best_score = -1.0
    best_w = 0.0
    
    if args.auto_tune:
        print(f"\n[*] Starting Dynamic Auto-Tuner (Hill Climbing Optimization)...")
        current_w = 0
        step_size = 0.2
        patience = 50       # Reverted to a reasonable hill-climb timeout to prevent endless loops
        stagnant_steps = 0
        
        while True:
            test_df = _run_differential_diet(STATIC_DB, dyn, promote_list, inhibit_list, args.inhibit_weight, current_w, args.max_foods)
            if test_df.empty: break
            
            tot_prom = test_df["gross_promote_gain"].sum()
            tot_col = test_df["gross_collateral_penalty"].sum()
            
            # THE DIET COMPLETENESS MULTIPLIER
            # If it only finds 1 food out of 5, its score is divided by 5!
            completeness = len(test_df) / args.max_foods
            
            # The Specificity-Weighted Index: (Ratio) * sqrt(Good) * Completeness
            specificity_ratio = tot_prom / max(1.0, tot_col)
            thera_index = specificity_ratio * math.sqrt(max(0.0, tot_prom)) * completeness
            
            print(f"[*] TESTED WEIGHT {current_w:5.2f}x | Foods: {len(test_df)}/{args.max_foods} | Promote: {tot_prom:<6.1f} | Collateral: {tot_col:<6.1f} | Index: {thera_index:.2f}")
            
            # Check for improvement (must improve by a margin to count as non-stagnant)
            if thera_index > best_score + 0.01:
                best_score = thera_index
                best_df = test_df
                best_w = current_w
                stagnant_steps = 0  # Reset patience!
            else:
                stagnant_steps += 1
                
            if stagnant_steps >= patience:
                print(f"    [STOP] No improvement for {patience} consecutive steps. True peak found.")
                break
                
            # Safety abort: If the weight is so high the diet is starving the target bug entirely
            if tot_prom < 5.0:
                print(f"    [STOP] Diet collapsed (Target Gain too low).")
                break
                
            current_w += step_size
            
    else:
        # Run standard single weight
        print(f"\n[*] Running single evaluation at Collateral Weight: {args.collateral_weight}x")
        best_w = args.collateral_weight
        best_df = _run_differential_diet(STATIC_DB, dyn, promote_list, inhibit_list, args.inhibit_weight, best_w, args.max_foods)
        if not best_df.empty:
            tot_prom = best_df["gross_promote_gain"].sum()
            tot_col = best_df["gross_collateral_penalty"].sum()
            best_score = (tot_prom ** 2) / max(1.0, (tot_prom + tot_col))

    if best_df is not None and not best_df.empty:
        print(f"\n[+] OPTIMAL DIET EXTRACTED at Collateral Weight: {best_w:.2f}x (Therapeutic Index: {best_score:.1f})")
        
        fm = load_smart(PATH_FOOD_CSV).drop_duplicates(subset=["fdc_id"], keep="last")
        if "food_category" in fm.columns: fm = fm.drop(columns=["food_category"])
        cat = load_smart(PATH_FOOD_CATEGORY_CSV).rename(columns={"id":"food_category_id", "description":"food_category"})
        fm["food_category_id"] = fm.get("food_category_id", pd.Series([""]*len(fm))).astype(str)
        cat["food_category_id"] = cat["food_category_id"].astype(str)
        fm = fm.merge(cat[["food_category_id", "food_category"]], on="food_category_id", how="left")
        
        out_df = best_df.drop(columns=["description", "food_category"], errors="ignore")
        out_df = round_floats(out_df.merge(fm[["fdc_id", "description", "food_category"]], left_on="representative_fdc_id", right_on=pd.to_numeric(fm["fdc_id"], errors="coerce"), how="left"), 4)
        
        final_cols = ["intervention", "rank", "description", "food_category", "top_beneficial_nutrients", "differential_score", "net_gain", "gross_promote_gain", "individual_gains", "gross_inhibit_penalty", "gross_collateral_penalty", "canonical_fdc_id", "representative_fdc_id", "n_improved_promote", "redundancy"]
        
        diet_path = Path(args.out_prefix).with_suffix(".differential_diet.tsv")
        out_df[final_cols].to_csv(diet_path, sep="\t", index=False)
        print(f"[*] Wrote successful optimal diet to {diet_path}")
        
        col_path = Path(args.out_prefix).with_suffix(".collateral_impact.tsv")
        _generate_collateral_report(STATIC_DB, dyn, out_df["canonical_fdc_id"].dropna().astype(int).tolist(), col_path)
        
    else:
        print("[!] No diet could be found.")

if __name__ == "__main__": main()