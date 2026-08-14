#!/usr/bin/env python3
"""ec2food.py — EC -> foods via greedy weighted cover (v45 - Final Container Edition)"""
from __future__ import annotations
import argparse, math, os, re, gc, pickle, heapq, zlib, shutil
from pathlib import Path
import multiprocessing as mp
import pandas as pd
import numpy as np
import pyarrow.dataset as ds
try: import resource
except ImportError: resource = None

# ==============================================================================
# PIPELINE CONFIGURATION (HARDCODED FOR CONTAINER ISOLATION)
# ==============================================================================
PATH_NUTRIENT_TO_EC     = "../0_building/3_nutrient_to_ec.tsv"
PATH_NUTRIENT_CSV       = "/data/bac2food/nutrient.csv"
PATH_FOOD_CSV           = "/data/bac2food/food.parquet"
PATH_FOOD_CATEGORY_CSV  = "/data/bac2food/food_category.csv"
PATH_FOOD_PORTION_CSV   = "/data/bac2food/food_portion.csv"
PATH_BUCKETED_DIR       = "/data/bac2food/food_nutrient_bucketed"
PATH_INDEX_DIR          = "/data/bac2food/index_modeled"
PATH_NUTRIENT_ALIAS     = "../0_building/1_expanded_nutrients.tsv"
PATH_FOODB_DIR          = "/data/bac2food/foodb_ca"

DROP_BRANDED = True
ALLOW_MACRO_PROXY = False  # Mathematically strict: Requires explicit measurements
ALLOW_MACRO_SCAN = False   # Disabled by default for precision
# ==============================================================================

###############################################################################
# MATH PARAMS - TOUCH ONLY IN CASE OF MATH EMERGENCY - 

EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")
BUCKETS = 256
TAU_Q, MIN_CONTRIB_FRAC, K_MIN, GAIN_MIN, SCORE_MIN = 0.05, 0.01, 1, 0.001, -2.0
PROC_W, BROAD_W, BROAD_Q, TYPE_W, ART_W = 0.20, 0.20, 2.0, 0.25, 8.0
K_PLANT, K_OTHER, MAX_CANONS = 800, 200, 250000
PROC_BASE_W, JUNK_BASE_W, ALC_BASE_W = 0.15, 0.20, 0.30

# --- DIVERSITY TUNING ---
LAMBDA_FOOD, OVERLAP_W = 0.08, 0.60
FIBER_WEIGHT = 2.0
# ------------------------
################################################################################
OLIGO_FLOOR_MG, ISOFLAVONE_FLOOR_MG = 10.0, 1.0
OLIGO_IDS, ISOFLAVONE_IDS = {1076,1077,2063}, {2049,2050,2051}

CAT_PENALTY = {
    "Fast Foods":5.0,"Restaurant Foods":5.0,"Sweets":5.0,"Snacks":4.0,"Alcoholic Beverages":5.0,
    "Baked Products":3.0,"Beverages":4.0,"Meals, Entrees, and Side Dishes":1.2,
    "Soups, Sauces, and Gravies":3.0,"Baby Foods":4.0,"Spices and Herbs":5.0,"Fats and Oils":4.0,
    "Beef Products":3.0,"Pork Products":3.0,"Poultry Products":3.0,"Finfish and Shellfish Products":3.0,
    "Lamb, Veal, and Game Products":3.0,"Sausages and Luncheon Meats":3.0,"American Indian/Alaska Native Foods":3.0,
    "Fruits and Fruit Juices":0.05,"Cereal Grains and Pasta":0.30,"Nut and Seed Products":0.10,
    "Vegetables and Vegetable Products":0.0,"Legumes and Legume Products":0.0,
}
PLANT_CATS = {"Vegetables and Vegetable Products","Legumes and Legume Products",
              "Fruits and Fruit Juices","Cereal Grains and Pasta","Nut and Seed Products"}
TP_MAP = {"foundation_food":600,"sr_legacy_food":500,"survey_fndds_food":450,
          "experimental_food":300,"branded_food":100}
DANGEROUS_RULES = {"protein_to_amino_acids","fat_to_fatty_acids","carb_to_sugars","minerals_total_to_specific","amino_acids_total_to_specific"}
SAFE_RULE_PFXS  = {"vitamin_","folate_","niacin_","tocopherol_","tocotrienol_","salt_","isomer_","stereo_","hydrate_","acid_base_","fiber_","carb_to_complex"}

ALCOHOL_RE  = re.compile(r"\b(?:vodka|gin|rum|whiskey|whisky|bourbon|scotch|brandy|cognac|tequila|mezcal|liqueur|beer|lager|ale|cider|wine|champagne|prosecco|vermouth|sake|martini|margarita|mojito|daiquiri|pina colada|long island|bloody mary|cosmopolitan|mimosa|sangria|\balcohol\b|\bethanol\b)\b",re.I)
JUNK_RE     = re.compile(r"\b(?:pastry|fries|frosting|cand(?:y|ies)|syrup|concentrate|dessert|cookie|cake|pie|doughnut|soda|sugar|sweetened|candied|glaze|marshmallow|ice cream|pudding|whip(?:ped)?|topping|chocolate|caramel|gumm(?:y|ies)|crisps?|chips?|puffs?)\b",re.I)
POWDER_RE   = re.compile(r"\b(?:freeze-dried|dehydrated|powder|extract|dried|instant|dry mix|beverage mix)\b",re.I)
PROCESSED_RE= re.compile(r"\b(?:deep[-\s]?fried|pan[-\s]?fried|fried|hash brown|breaded|battered)\b",re.I)
_FIBER      = re.compile(r"\b(?:fiber|starch|oligosaccharide|inulin|pectin|glucan|hemicellulose|cellulose|pentosan)\b",re.I)
NFY_RE      = re.compile(r"(?:^\s*(?:sugars|fatty acids|amino acids|vitamin [a-z0-9 ]+|minerals|organic acids)\s*,)|(?:\s-\s*NFY)")
WHALE_RE    = re.compile(r"whale|seal|muktuk|walrus|blubber|bowhead",re.I)
ALWAYS_DROP_CATS = {"American Indian/Alaska Native Foods"}
_BLOCKED = re.compile(
    r"\b(?:sucrose|glucose|fructose|lactose|maltose|galactose|mannitol|xylitol|sorbitol|sugar|alanine|arginine|asparagine|aspartic acid|cysteine|glutamic acid|glutamine|glycine|histidine|isoleucine|leucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|cystine|cholesterol|sfa|mufa|pufa|tfa|trans|fatty acid|fatty acids|sterol|dha|epa|docosahexaenoic|eicosapentaenoic|polyunsaturated|monounsaturated|saturated|triglyceride|lipid|lipids|phospholipids?|glycolipids?|unsaponifiable|trans-monoenoic|trans-dienoic|trans-polyenoic|vitamin\s*[abcdefjkh][-\s\d]*|vitamin\s*b[-\s]?(?:[1-9]|12)|vitamins\b|b[-\s]?2\b|b[-\s]?6\b|b[-\s]?12\b|riboflavin|thiamin|niacin|folate|folic|tetrahydrofolate|5-mthf|pantothenic|biotin|pyridoxin|cobalamin|retinol|retinoic|carotene|tocopherols?|tocotrienols?|total\s+tocopherols?|ascorbic|menaquinone|phylloquinone|calcium|iron|magnesium|phosphorus|potassium|sodium|zinc|copper|selenium|manganese|iodine|choline|betaine|chromium|fluoride|molybdenum|c\d{1,2}:\d+|nacl|salt|salt deklaration|salt labelling|niacin[æa]kvivalent|niacine, pr[ée]form[ée]e|polyols|tannins?|flavans?|flavones?|flavonols?|anthocyanins?|haem|nhaem|val|thr|tyr|glu|gly|his|pro|ser|arg|asp|phe|ala|met|cys|ile|leu|trp|valin|threonin|tyrosin|glycin|prolin|serin|arginin|phenylalanin|alanin|methionin|cystein|isoleucin|leucin|hydroxyproline|glycogen|tyramin|putrescin|chitin|polyphenol|nutrient_name)\b|(?:aa \(pernitrogen\)|starch-|pp-|polymers \(>10 mers\))", 
    re.I
)
load_smart    = lambda p: pd.read_parquet(str(p)) if str(p).lower().endswith((".parquet",".pq")) else pd.read_csv(str(p),low_memory=False)
classify_rule = lambda r: (False,"macro_proxy") if (r:=(r or "").strip().lower()) in DANGEROUS_RULES else ((True,"form") if any(r.startswith(p) for p in SAFE_RULE_PFXS) else (False,"unknown_proxy"))
get_cov       = lambda x,th: 0.0 if th<=0 else min(max(x/th,0.0),1.0)
fmt_nlist     = lambda ids,nd,lim=5: ",".join([nd.get(t,str(t)) for t in ids][:lim]) + (f" (+{len(ids)-lim} more)" if len(ids)>lim else "")

def cosine_sim(a,b):
    if not a or not b: return 0.0
    i=set(a)&set(b)
    if not i: return 0.0
    dot=sum(a[k]*b[k] for k in i); na=math.sqrt(sum(v*v for v in a.values())); nb=math.sqrt(sum(v*v for v in b.values()))
    return dot/(na*nb) if na>0 and nb>0 else 0.0

def round_floats(df,nd=4):
    if df is None or df.empty: return df
    c=df.select_dtypes(include=["float16","float32","float64"]).columns
    if len(c): df[c]=df[c].round(nd)
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

def build_static_food_meta(args,out_path):
    print("[*] Building Dense Numpy Array Index...",flush=True)
    port_df=load_smart(PATH_FOOD_PORTION_CSV)
    if "gram_weight" not in port_df.columns and "amount" in port_df.columns: port_df=port_df.rename(columns={"amount":"gram_weight"})
    port_df=port_df.dropna(subset=["gram_weight"])
    port_df["is_serving"]=port_df.get("modifier",pd.Series([""]*len(port_df))).astype(str).str.lower().str.contains("serving")
    portion_map=port_df.sort_values(["fdc_id","is_serving"],ascending=[True,False]).groupby("fdc_id")["gram_weight"].first().to_dict()
    nut=load_smart(PATH_NUTRIENT_CSV).rename(columns={"id":"nutrient_id","name":"nutrient_name"})
    unit_map=nut.set_index("nutrient_id")["unit_name"].astype("string").str.upper().to_dict()
    modeled=set(nut["nutrient_id"].dropna().astype(int).tolist())
    id2desc={}
    for _,row in load_smart(PATH_FOOD_CATEGORY_CSV).rename(columns={"id":"food_category_id","description":"food_category"}).iterrows():
        raw=str(row["food_category_id"]).strip(); name=str(row["food_category"]); id2desc[raw]=name
        try: id2desc[str(int(float(raw)))]=name
        except: pass
    food_df=load_smart(PATH_FOOD_CSV).drop_duplicates(subset=["fdc_id"],keep="last")
    meta=food_df.set_index(pd.to_numeric(food_df["fdc_id"],errors="coerce"))[["description","data_type","food_category_id"]].to_dict(orient="index")
    max_fdc=int(max(k for k in meta if pd.notna(k)))+1000
    c_arr=np.full(max_fdc,-1,dtype=np.int32); p_arr=np.full(max_fdc,50.0,dtype=np.float32)
    pl_arr=np.zeros(max_fdc,dtype=bool); sp_arr=np.zeros(max_fdc,dtype=bool)
    drop_cats=set(str(x).strip() for x in args.drop_category); food_stats={}
    for fdc_id,m in meta.items():
        if pd.isna(fdc_id): continue
        fid=int(fdc_id); dt=str(m.get("data_type","")); dt_norm=dt.lower().strip().replace(" ","_")
        cat_raw=str(m.get("food_category_id","")).strip(); cat=id2desc.get(cat_raw,"")
        if not cat:
            try: cat=id2desc.get(str(int(float(cat_raw))),"")
            except (ValueError,OverflowError): pass
        desc_lc=str(m.get("description","")).lower().strip()
        if (DROP_BRANDED and dt_norm=="branded_food") or cat in drop_cats or cat in ALWAYS_DROP_CATS or WHALE_RE.search(desc_lc): continue
        c_arr[fid]=fid; p_arr[fid]=portion_map.get(fid,50.0)
        sp_arr[fid]=(cat=="Spices and Herbs"); pl_arr[fid]=(cat in PLANT_CATS)

        if "[phenol-explorer]" in desc_lc:
            cat = "Fruits and Fruit Juices" if not cat else cat
            p = 0.0
            bc = 0.0
            pl_arr[fid] = True
        elif "[biofoodcomp]" in desc_lc:
            if not cat:
                if any(k in desc_lc for k in ["bean", "pea", "lentil", "chickpea", "soy"]):
                    cat = "Legumes and Legume Products"
                elif any(k in desc_lc for k in ["berry", "fruit", "apple", "grape", "apricot", "peach", "plum"]):
                    cat = "Fruits and Fruit Juices"
                elif any(k in desc_lc for k in ["grain", "bran", "rice", "barley", "wheat", "rye", "oat", "quinoa"]):
                    cat = "Cereal Grains and Pasta"
                else:
                    cat = "Vegetables and Vegetable Products"
            p = float(CAT_PENALTY.get(cat, 0.0))
            bc = PROC_BASE_W * p + 0.40
            pl_arr[fid] = True
        else:
            p = float(CAT_PENALTY.get(cat, 1.0))
            bc = PROC_BASE_W * p
            
        if NFY_RE.search(desc_lc): p=max(p,1.5)
        if ALCOHOL_RE.search(desc_lc): p=max(p,2.0); bc+=ALC_BASE_W
        if JUNK_RE.search(desc_lc) or POWDER_RE.search(desc_lc): p=max(p,2.0); bc+=(JUNK_BASE_W if JUNK_RE.search(desc_lc) else 0.40)
        if PROCESSED_RE.search(desc_lc): bc+=0.35
        tp=TP_MAP.get(dt_norm,100)
        food_stats[fid]={"dh":int(zlib.crc32(desc_lc.encode())),"pl":p,"bc":bc,"art":int(p>=1.5),"tb":TYPE_W*(tp/600.0),"tp":tp,"dt_norm":dt_norm,"cat":cat}
    if DROP_BRANDED:
        assert all(fs["dt_norm"]!="branded_food" for fs in food_stats.values()),"Branded foods leaked into the index!"
    canon_best={}
    for f,fs in food_stats.items():
        c=c_arr[f]
        if canon_best.get(c,-1)<fs["tp"]: canon_best[c]=fs["tp"]
    valid_set={f for f,fs in food_stats.items() if canon_best.get(c_arr[f],-1)==fs["tp"]}
    vm=np.zeros(max_fdc,dtype=bool)
    if valid_set: vm[list(valid_set)]=True
    c_arr[~vm]=-1
    max_nut=max(modeled)+100; u_arr=np.ones(max_nut,dtype=np.float32)
    for n in modeled:
        u=unit_map.get(n,"").upper(); u_arr[n]=1000.0 if u=="G" else (0.001 if u=="UG" else 1.0)
    pickle.dump({"canon_arr":c_arr,"port_arr":p_arr,"plant_arr":pl_arr,"spice_arr":sp_arr,"unit_arr":u_arr,
                 "food_stats":{f:food_stats[f] for f in valid_set},"modeled":modeled,
                 "nut_name":dict(zip(nut["nutrient_id"],nut["nutrient_name"].astype(str))),
                 "rep2can":{f:c_arr[f] for f in valid_set}},open(out_path,"wb"))
    print("[*] Static Food Meta successfully created.",flush=True)

def build_modeled_index(bdir,static_db,out_dir,batch_rows=750_000):
    out_dir.mkdir(parents=True,exist_ok=True)
    nids=sorted(int(x) for x in static_db["modeled"]); buckets=sorted({n%BUCKETS for n in nids})
    scn=ds.dataset(bdir,format="parquet",partitioning="hive").scanner(
        columns=["fdc_id","nutrient_id","amount"],
        filter=ds.field("bucket").isin(buckets)&ds.field("nutrient_id").isin(nids),batch_size=batch_rows)
    tot_c,tot_m,df_foods={},{},{}
    tmp=out_dir/"_tmp_scale_parts"
    if tmp.exists(): shutil.rmtree(tmp)
    tmp.mkdir(parents=True,exist_ok=True)
    c_arr,p_arr,sp_arr,u_arr=static_db["canon_arr"],static_db["port_arr"],static_db["spice_arr"],static_db["unit_arr"]
    max_f,max_u=len(c_arr),len(u_arr)
    for bi,b in enumerate(scn.to_batches(),1):
        f,n,a=b["fdc_id"].to_numpy().astype(np.int64,copy=False),b["nutrient_id"].to_numpy().astype(np.int32,copy=False),b["amount"].to_numpy().astype(np.float32,copy=False)
        m=(f>=0)&(f<max_f)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f,n,a=f[m],n[m],a[m]; c=c_arr[f]; v=(c!=-1)&(~sp_arr[f])
        if not v.any(): continue
        f,n,a,c=f[v],n[v],a[v],c[v]; a=a*(p_arr[f]/50.0)*u_arr[n]
        df=pd.DataFrame({"nutrient_id":n,"canon":c,"amount_norm":a})
        for canon,row in df.groupby("canon",dropna=True).agg(cnt=("nutrient_id","size"),mass=("amount_norm","sum")).iterrows():
            cn=int(canon); tot_c[cn]=tot_c.get(cn,0)+int(row["cnt"]); tot_m[cn]=tot_m.get(cn,0.0)+float(row["mass"])
        for nid,canon in df[["nutrient_id","canon"]].drop_duplicates().itertuples(index=False): df_foods.setdefault(int(nid),set()).add(int(canon))
        sub2=df.groupby(["nutrient_id","canon"],dropna=True)["amount_norm"].max().reset_index()
        if not sub2.empty: sub2["bucket"]=(sub2["nutrient_id"].astype(int)%BUCKETS).astype(int); sub2.to_parquet(tmp/f"part_{bi:06d}.parquet",index=False)
        del df,f,n,a,c
        if bi%10==0: gc.collect()
    pd.DataFrame({"canon":list(tot_c),"model_count":[int(tot_c[c]) for c in tot_c],"model_mass":[float(tot_m.get(c,0)) for c in tot_c]}).to_parquet(out_dir/"modeled_totals.parquet",index=False)
    pd.DataFrame({"nutrient_id":list(df_foods),"df_foods":[len(df_foods[n]) for n in df_foods]}).to_parquet(out_dir/"nutrient_df.parquet",index=False)
    rows=[]; tmp_ds=ds.dataset(str(tmp),format="parquet")
    for nb in range(BUCKETS):
        tab=tmp_ds.to_table(filter=ds.field("bucket")==nb,columns=["nutrient_id","canon","amount_norm"])
        if tab.num_rows==0: continue
        med=tab.to_pandas().groupby(["nutrient_id","canon"],as_index=False)["amount_norm"].max().groupby("nutrient_id",as_index=False)["amount_norm"].median()
        rows.extend({"nutrient_id":int(r.nutrient_id),"ref_amount_norm":float(r.amount_norm)} for r in med.itertuples(index=False))
    pd.DataFrame(rows).to_parquet(out_dir/"nutrient_scale.parquet",index=False); shutil.rmtree(tmp,ignore_errors=True)

STATIC_DB,DYNAMIC_STATE={},{}
def _init_w(p,d): global STATIC_DB,DYNAMIC_STATE; STATIC_DB=pickle.load(open(p,"rb")); DYNAMIC_STATE=d

def _run_lbl(lbl):
    all_t=set(int(x) for x in DYNAMIC_STATE["lbl2targ"].get(lbl,set()))
    targs=(all_t & STATIC_DB["modeled"]) - DYNAMIC_STATE["blocked"]
    if not targs: return pd.DataFrame(),pd.DataFrame()
    need=set(targs)
    for t in targs:
        for p in DYNAMIC_STATE["prox"].get(t,[]):
            if p[2] or (DYNAMIC_STATE["ams"] and DYNAMIC_STATE["amp"] and DYNAMIC_STATE["ndf"].get(t,0)<50000 and p[3]=="macro_proxy"): need.add(p[0])
    def _tau(n): return min(max(0.1 if STATIC_DB["unit_arr"][n]==1.0 else 0.01,(0.01 if DYNAMIC_STATE["ndf"].get(n,1e9)<50000 else TAU_Q)*DYNAMIC_STATE["ref"].get(n,0.0)),500.0 if STATIC_DB["unit_arr"][n]==1.0 else 50.0)
    tau={n:_tau(n) for n in targs}
    scn=ds.dataset(DYNAMIC_STATE["bdir"],format="parquet",partitioning="hive").scanner(
        columns=["fdc_id","nutrient_id","amount"],
        filter=ds.field("bucket").isin(sorted({n%BUCKETS for n in need}))&ds.field("nutrient_id").isin(sorted(need)),batch_size=50_000)
    tk_p,tk_o=TopKByNutrient(),TopKByNutrient()
    c_arr,p_arr,pl_arr,sp_arr,u_arr=STATIC_DB["canon_arr"],STATIC_DB["port_arr"],STATIC_DB["plant_arr"],STATIC_DB["spice_arr"],STATIC_DB["unit_arr"]
    max_f,max_u=len(c_arr),len(u_arr)
    for b in scn.to_batches():
        f,n,a=b["fdc_id"].to_numpy().astype(np.int64,copy=False),b["nutrient_id"].to_numpy().astype(np.int32,copy=False),b["amount"].to_numpy().astype(np.float32,copy=False)
        m=(f>=0)&(f<max_f)&(n>=0)&(n<max_u)&np.isfinite(a)
        if not m.any(): continue
        f,n,a=f[m],n[m],a[m]; c=c_arr[f]; v=(c!=-1)
        if not DYNAMIC_STATE["asp"]: v&=~sp_arr[f]
        if not v.any(): continue
        f,n,a,c=f[v],n[v],a[v],c[v]; ports=p_arr[f].copy()
        if DYNAMIC_STATE["asp"]: ports[sp_arr[f]]=2.0
        a=a*(ports/50.0)*u_arr[n]
        idx=np.lexsort((a,c,n)); n,c,a,f=n[idx],c[idx],a[idx],f[idx]
        m2=np.empty(len(n),dtype=bool); m2[-1]=True; m2[:-1]=(n[:-1]!=n[1:])|(c[:-1]!=c[1:])
        n,c,a,ipl=n[m2],c[m2],a[m2],pl_arr[f[m2]]
        for i in range(len(n)): (tk_p if ipl[i] else tk_o).push(int(n[i]),int(c[i]),float(a[i]),K_PLANT if ipl[i] else K_OTHER)
    def get_tk(nid):
        d=dict(tk_p.val.get(nid,{}))
        for k,v in tk_o.val.get(nid,{}).items():
            if k not in d or v>d[k]: d[k]=v
        return d
    c2a={}
    for t in targs:
        th=float(tau.get(t,0.0)); nn=STATIC_DB["nut_name"].get(t,"")
        mc=0.0001 if _FIBER.search(nn) else (0.001 if DYNAMIC_STATE["ndf"].get(t,1e9)<50000 else MIN_CONTRIB_FRAC)
        base_cut=th*mc; is_sp=(t in OLIGO_IDS or t in ISOFLAVONE_IDS)
        cut=max(base_cut,OLIGO_FLOOR_MG if t in OLIGO_IDS else ISOFLAVONE_FLOOR_MG if t in ISOFLAVONE_IDS else 0)
        for c,a in get_tk(int(t)).items():
            if a>cut: c2a.setdefault(c,{})[int(t)]=float(a)
        if not is_sp and DYNAMIC_STATE["ndf"].get(t, 0) < 5000:
            for p in [x for x in DYNAMIC_STATE["prox"].get(t,[]) if x[2] or (DYNAMIC_STATE["ams"] and DYNAMIC_STATE["amp"] and x[3]=="macro_proxy")]:
                for c,a in get_tk(int(p[0])).items():
                    v=(float(a)/p[1])*0.25
                    if v>base_cut: d=c2a.setdefault(c,{}); d[int(t)]=v if int(t) not in d or v>d[int(t)] else d[int(t)]
    c2a={c:d for c,d in c2a.items() if d}; del tk_p,tk_o; gc.collect()
    m_ram=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0 if resource else 0.0
    if not c2a:
        print(f"[DEBUG] {lbl}: 0 candidates. (RAM:{m_ram:.0f}MB)",flush=True)
        return pd.DataFrame(),pd.DataFrame([dict(species=lbl,nutrient_id=n,nutrient_name=STATIC_DB["nut_name"].get(n,""),coverage=0,cur_amount=0,tau=tau.get(n,0),reason="not_covered") for n in targs])
    print(f"[DEBUG] {lbl}: {len(c2a)} candidates. RAM:{m_ram:.0f}MB",flush=True)
    
    w={n:(FIBER_WEIGHT if _FIBER.search(STATIC_DB["nut_name"].get(n,"")) else 1.0)*float(DYNAMIC_STATE["nmult"].get(n,1.0))*(1.0+math.log(DYNAMIC_STATE["N"]/max(1,int(DYNAMIC_STATE["ndf"].get(n,1))))) for n in targs}
    
    if len(c2a)>MAX_CANONS:
        h=[]
        for c,ct in c2a.items():
            s=sum(w.get(n,1.0)*math.log1p(a/tau.get(n,1e-9)) for n,a in ct.items() if tau.get(n,0)>0)
            if len(h)<MAX_CANONS: heapq.heappush(h,(s,c))
            elif s>h[0][0]: heapq.heapreplace(h,(s,c))
        c2a={c:c2a[c] for _,c in h}
        
    ST=type('S',(),{'cur':{n:0.0 for n in targs},'chosen':set(),'seen_desc':set(),'chosen_imps':[],'rows':[],'chosen_cats':[]})()
    fs_db=STATIC_DB["food_stats"]
    
    ni_min = 1 if len(targs) <= 5 else min(5, max(2, int(0.03 * len(targs))))
    
    nd=STATIC_DB["nut_name"]
    measured_ids=[t for t in targs if DYNAMIC_STATE["ndf"].get(t,0)>0]
    missing_ids=[t for t in targs if DYNAMIC_STATE["ndf"].get(t,0)==0]
    n_meas,n_miss=len(measured_ids),len(missing_ids)
    miss_str=fmt_nlist(missing_ids,nd); meas_str=fmt_nlist(measured_ids,nd)
    for rk in range(1,DYNAMIC_STATE["maxf"]+1):
        bst=None
        for c,ct in c2a.items():
            fs=fs_db.get(c,{})
            if c in ST.chosen or not fs or fs["dh"] in ST.seen_desc: continue
            gn,ni,imp=0.0,0,{}
            for n,a in ct.items():
                th=tau[n]; mcf=0.0001 if _FIBER.search(nd.get(n,"")) else (0.001 if DYNAMIC_STATE["ndf"].get(n,1e9)<50000 else MIN_CONTRIB_FRAC)
                floor=OLIGO_FLOOR_MG if n in OLIGO_IDS else (ISOFLAVONE_FLOOR_MG if n in ISOFLAVONE_IDS else 0)
                if a<=max(th*mcf,floor): continue
                dg=math.log1p((ST.cur[n]+a)/th)-math.log1p(ST.cur[n]/th)
                if dg>1e-12: cg=w.get(n,1.0)*dg; gn+=cg; ni+=1; imp[n]=cg
            if ni<ni_min or gn<GAIN_MIN: continue
            imp={n:v for n,v in imp.items() if v>=0.05*gn}
            red=sum(sorted([cosine_sim(imp,p) for p in ST.chosen_imps],reverse=True)[:3])/3.0 if ST.chosen_imps else 0.0
            prm=min(1.0,max(0.0,ni/max(50,int(DYNAMIC_STATE["mc"].get(c,0)))))
            geff=math.log1p(gn)*(math.sqrt(ni/(ni+2.0)) if ni>0 else 0.0)
            
            cat_rep = ST.chosen_cats.count(fs.get("cat", ""))
            cat_pen = 0.20 * cat_rep
            
            proc_cost=PROC_W*fs["pl"]*math.log1p(gn); broad_cost=BROAD_W*math.log1p(gn)*((1.0-prm)**BROAD_Q)
            
            sc=geff-LAMBDA_FOOD-proc_cost-broad_cost-(ART_W*fs["art"])-(OVERLAP_W*red)-fs["bc"]+(fs["tb"]*geff)-cat_pen
            
            if sc>SCORE_MIN and (bst is None or sc>bst[0]):
                top=sorted(imp.items(),key=lambda x:x[1],reverse=True)[:3]
                bst=(sc,gn,c,math.log1p(gn),ni,prm,fs["pl"],fs["art"],fs["tb"],geff,proc_cost,broad_cost,red,imp,
                     ",".join(str(n) for n,_ in top),",".join(nd.get(n,str(n)) for n,_ in top),fs["dh"])
        if not bst: break
        c=bst[2]; ST.chosen.add(c); ST.seen_desc.add(bst[-1]); ST.chosen_imps.append(bst[13])
        ST.chosen_cats.append(fs_db[c].get("cat", "")) 
        
        for n,a in c2a[c].items(): ST.cur[n]+=a
        cvs={t:get_cov(ST.cur[t],tau.get(t,0.0)) for t in targs}
        cov_tot=sum(1 for v in cvs.values() if v>=1.0)
        cov_meas=sum(1 for t,v in cvs.items() if DYNAMIC_STATE["ndf"].get(t,0)>0 and v>=1.0)
        mcvs=[v for t,v in cvs.items() if DYNAMIC_STATE["ndf"].get(t,0)>0]
        ST.rows.append(dict(rank=rk,species=lbl,description="",food_category="",
            top_nutrient_names=bst[15],score=bst[0],gain=bst[1],
            n_targets=len(targs),n_targets_measured=n_meas,n_targets_missing=n_miss,
            covered_targets_measured=cov_meas,coverage_measured_frac=cov_meas/max(1,n_meas),
            top_nids=bst[14],canonical_fdc_id=c,representative_fdc_id=c,data_type="",
            norm_gain=bst[3],n_improved=bst[4],purity=bst[5],proc_level=bst[6],
            proc_cost=bst[10],broad_cost=bst[11],artifact=bst[7],artifact_cost=ART_W*bst[7],type_bonus=bst[8],
            breadth_factor=0,gain_eff=bst[9],redundancy=bst[12],redundancy_cost=OVERLAP_W*bst[12],fallback=0,fallback_reason="",
            covered_targets_total=cov_tot,coverage_total_frac=cov_tot/max(1,len(targs)),
            avg_coverage_measured=sum(mcvs)/max(1,len(mcvs)),min_coverage_measured=min(mcvs) if mcvs else 0.0,
            measured_nutrient_names=meas_str,missing_nutrient_names=miss_str))
    unc=[dict(species=lbl,nutrient_id=n,nutrient_name=nd.get(n,""),coverage=get_cov(ST.cur.get(n,0),tau.get(n,0)),
              cur_amount=ST.cur.get(n,0),tau=tau.get(n,0),reason="not_covered")
         for n in targs if get_cov(ST.cur.get(n,0),tau.get(n,0))<1.0]
    return pd.DataFrame(ST.rows),pd.DataFrame(unc)

def rescue_missing_targets(unc_df, foodb_dir):
    """Scan FooDB for missing nutrients and append results to the DataFrame."""
    if not foodb_dir or not os.path.exists(foodb_dir):
        print(f"[*] Skipping FooDB Rescue. Directory not found: {foodb_dir}", flush=True)
        return unc_df
        
    print(f"[*] Initiating FooDB Rescue for missing targets...", flush=True)
    missing_targets = unc_df['nutrient_name'].dropna().unique()
    
    custom_map = {
        "inositol": "myo-inositol", 
        "phytic acid": "phytic acid",
        "vitamin b-12": "vitamin b12",
        "vitamin b-6": "vitamin b6",
        "vitamin c": "ascorbic acid",
        "folate": "folic acid",
        "retinol": "vitamin a",
        "beta carotene": "beta-carotene"
    }
    
    clean_name = lambda x: str(x).strip().lower() if pd.notna(x) else ""
    
    comp_df = pd.read_csv(os.path.join(foodb_dir, "Compound.csv"), usecols=['id', 'name'], low_memory=False)
    comp_df['clean'] = comp_df['name'].apply(clean_name)
    nut_df = pd.read_csv(os.path.join(foodb_dir, "Nutrient.csv"), usecols=['id', 'name'], low_memory=False)
    nut_df['clean'] = nut_df['name'].apply(clean_name)
    syn_df = pd.read_csv(os.path.join(foodb_dir, "CompoundSynonym.csv"), usecols=['synonym', 'source_id', 'source_type'], low_memory=False)
    syn_df['clean'] = syn_df['synonym'].apply(clean_name)

    nutrient_to_foodb = {}
    for orig_name in missing_targets:
        cl = clean_name(orig_name)
        lookups = [cl]
        if cl in custom_map: lookups.append(custom_map[cl])
        matches = set()
        for lk in lookups:
            for i in comp_df[comp_df['clean'] == lk]['id'].tolist(): matches.add((i, 'Compound'))
            for i in nut_df[nut_df['clean'] == lk]['id'].tolist(): matches.add((i, 'Nutrient'))
            for i in syn_df[(syn_df['clean'] == lk) & (syn_df['source_type'] == 'Compound')]['source_id'].unique()[:5]: matches.add((i, 'Compound'))
        if matches: nutrient_to_foodb[orig_name] = list(matches)

    if not nutrient_to_foodb:
        unc_df['FooDB_Rescue_Candidates'] = "No match in FooDB"
        return unc_df

    content_path = os.path.join(foodb_dir, "Content.csv")
    needed_compounds = set(m[0] for matches in nutrient_to_foodb.values() for m in matches if m[1] == 'Compound')
    needed_nutrients = set(m[0] for matches in nutrient_to_foodb.values() for m in matches if m[1] == 'Nutrient')
    
    food_hits = []
    cols_to_use = ['source_id', 'source_type', 'orig_food_common_name', 'standard_content']
    
    for chunk in pd.read_csv(content_path, usecols=cols_to_use, chunksize=250000, low_memory=False):
        c_mask = (chunk['source_type'] == 'Compound') & (chunk['source_id'].isin(needed_compounds))
        n_mask = (chunk['source_type'] == 'Nutrient') & (chunk['source_id'].isin(needed_nutrients))
        valid = chunk[c_mask | n_mask].copy()
        valid['standard_content'] = pd.to_numeric(valid['standard_content'], errors='coerce')
        valid = valid.dropna(subset=['standard_content'])
        valid = valid[valid['standard_content'] > 0]
        if not valid.empty: food_hits.append(valid)
            
    if not food_hits:
        unc_df['FooDB_Rescue_Candidates'] = "Found in ontology, but no quantified foods"
        return unc_df
        
    hits_df = pd.concat(food_hits, ignore_index=True)
    
    def get_foods(orig_name):
        if orig_name not in nutrient_to_foodb: return "No match in FooDB"
        targets = nutrient_to_foodb[orig_name]
        mask = pd.Series(False, index=hits_df.index)
        for sid, stype in targets:
            mask |= (hits_df['source_id'] == sid) & (hits_df['source_type'] == stype)
        sub = hits_df[mask]
        if sub.empty: return "Found in ontology, but no quantified foods"
        top_foods = sub.groupby('orig_food_common_name')['standard_content'].max().sort_values(ascending=False).head(5)
        return " | ".join([f"{f} ({amt:.2f} mg/100g)" for f, amt in top_foods.items()])

    unc_df['FooDB_Rescue_Candidates'] = unc_df['nutrient_name'].apply(get_foods)
    print("[*] FooDB Rescue complete!", flush=True)
    return unc_df

def main():
    print("==========================================================\nRUNNING ec2food_v45.py (Clean Container Edition)\n==========================================================",flush=True)
    if hasattr(mp,"set_start_method"):
        try: mp.set_start_method("forkserver",force=True)
        except RuntimeError: pass
        
    ap=argparse.ArgumentParser()
    ap.add_argument("--enzyme_tsv",required=True)
    ap.add_argument("--out_prefix",required=True)
    ap.add_argument("--jobs",type=int,default=max(1,os.cpu_count() or 1))
    ap.add_argument("--max_foods",type=int,default=10)
    ap.add_argument("--force_pick",action="store_true")
    ap.add_argument("--drop_category",action="append",default=[])
    ap.add_argument("--allow-spices",dest="allow_spices",action=argparse.BooleanOptionalAction,default=False)
    ap.add_argument("--rebuild-static-meta",action="store_true")
    args=ap.parse_args()
    sm_path=Path(PATH_INDEX_DIR)/"static_food_meta.pkl"
    nut_df=load_smart(PATH_NUTRIENT_CSV)
    blocked=set(nut_df[nut_df["name"].astype(str).str.contains(_BLOCKED,regex=True,na=False)]["id"].dropna().astype(int).tolist())
    if args.rebuild_static_meta and sm_path.exists(): sm_path.unlink()
    if not sm_path.exists(): build_static_food_meta(args,sm_path)
        
    enz=pd.read_csv(args.enzyme_tsv,sep="\t",dtype="string")
    cols={c.lower():c for c in enz.columns}
    if "ec_number" not in enz.columns: enz=enz.rename(columns={cols.get("ec_number",cols.get("ec",enz.columns[-1])):"ec_number"})
    if "species" not in enz.columns: enz=enz.rename(columns={enz.columns[0]:"species"})
    enz["strain"]=enz.get("strain",pd.Series([""]*len(enz))).fillna("").astype("string").str.strip()
    enz=enz[enz["ec_number"].str.match(EC_RE,na=False)].dropna(subset=["species","ec_number"])
    enz_nut=enz.merge(pd.read_csv(PATH_NUTRIENT_TO_EC,sep="\t"),on="ec_number")
    enz_nut["label"]=enz_nut["species"]+enz_nut["strain"].apply(lambda x:f" | {x}" if x else "")
    lbl2targ={lbl:set(sub["nutrient_id"].astype(int).unique()) for lbl,sub in enz_nut.groupby("label")}
    prox={}
    
    if PATH_NUTRIENT_ALIAS:
        a=pd.read_csv(PATH_NUTRIENT_ALIAS,sep="\t")
        div_map=a[a["rule"].str.strip().str.lower().isin(DANGEROUS_RULES)].groupby(["generic_nutrient_id","rule"])["specific_nutrient_id"].nunique().to_dict()
        for g,s,r in a[["generic_nutrient_id","specific_nutrient_id","rule"]].itertuples(index=False):
            so,k=classify_rule(str(r).strip().lower()); prox.setdefault(int(s),[]).append((int(g),float(div_map.get((int(g),r),1.0)),so,k))
    for t in [1017,1021,1022,1403,2058,1071,1019]: prox.setdefault(t,[]).append((1079,1.0,True,"form"))
    for t in [1015,1016,1020]: prox.setdefault(t,[]).append((1009,1.0,True,"form"))
    for t in [1181,1182,1042]: prox.setdefault(t,[]).append((99999,1.0,True,"form"))

    # ==============================================================================
    # THE AUTO-ALIASER: Instantly rescue specific synonyms (NO MACROS/AMINO ACIDS)
    # ==============================================================================
    name_to_id = {str(r.name).strip().lower(): int(r.id) for r in nut_df.itertuples() if pd.notna(r.name) and pd.notna(r.id)}
    
    auto_alias_map = {
            # --- Typos, Formatting, and Trace Elements ---
            "tyramin": "tyramine",
            "putrescin": "putrescine",
            "cholecalciferol (vit d3)": "cholecalciferol",
            "ferrulic acid": "ferulic acid",
            "zink": "zinc",
            "sulfur, s": "sulfur",
            "mercury, hg": "mercury",
            "rubidium, rb": "rubidium",
            "titanium, ti": "titanium",
            "5-formyltetrahydrofolic acid (5-hcoh4": "5-formyltetrahydrofolic acid",
            "tetrahydrofolic acid (thf)": "tetrahydrofolic acid",

            # --- Stereoisomers & Chemical Prefixes (Flattening to standard names) ---
            "(+)-catechin": "catechin",
            "(-)-epicatechin": "epicatechin",
            "(-)-epigallocatechin": "epigallocatechin",
            "(-)-epicatechin 3-o-gallate": "epicatechin gallate",
            "epicatechin-3-gallate": "epicatechin gallate",
            "p-hydroxy benzoic acid": "4-hydroxybenzoic acid",

            # --- The Vitamin K2 Series (Menaquinones) ---
            "mk4": "menaquinone-4",
            "mk5": "menaquinone-5",
            "mk6": "menaquinone-6",
            "mk7": "menaquinone-7",
            "mk8": "menaquinone-8",
            "mk9": "menaquinone-9",
            "mk10": "menaquinone-10",

            # --- Specific Molecule Shorthand ---
            "ip6": "phytic acid", # Inositol hexaphosphate is Phytic acid

            # --- The Organic Acid / Conjugate Base Pairs (-ate to -ic acid) ---
            "aconitate": "aconitic acid",
            "benzoate": "benzoic acid",
            "cinnamate": "cinnamic acid",
            "citrate": "citric acid",
            "formate": "formic acid",
            "fumarate": "fumaric acid",
            "gentisate": "gentisic acid",
            "gluconate": "gluconic acid",
            "glycolate": "glycolic acid",
            "isocitrate": "isocitric acid",
            "lactate": "lactic acid",
            "malonate": "malonic acid",
            "orotate": "orotic acid",
            "oxalate": "oxalic acid",
            "oxaloacetate": "oxaloacetic acid",
            "pantothenate": "pantothenic acid",
            "phytate": "phytic acid",
            "pyruvate": "pyruvic acid",
            "quinate": "quinic acid",
            "salicylate": "salicylic acid",
            "succinate": "succinic acid",
            "tartrate": "tartaric acid",
            "vanillate": "vanillic acid",
            "alpha-ketoglutarate": "alpha-ketoglutaric acid"
        }
    
    for alias, canon_name in auto_alias_map.items():
        if alias in name_to_id and canon_name in name_to_id:
            alias_id = name_to_id[alias]
            canon_id = name_to_id[canon_name]
            prox.setdefault(alias_id, []).append((canon_id, 1.0, True, "form"))
    # ==============================================================================

    idx_dir=Path(PATH_INDEX_DIR)
    if not all((idx_dir/f).exists() for f in ["modeled_totals.parquet","nutrient_df.parquet","nutrient_scale.parquet"]):
        build_modeled_index(PATH_BUCKETED_DIR,pickle.load(open(sm_path,"rb")),idx_dir)
    tot_pq=pd.read_parquet(idx_dir/"modeled_totals.parquet")
    dyn={"bdir":PATH_BUCKETED_DIR,"nmult":{},"lbl2targ":lbl2targ,
         "prox":{k:[(p[0],float(p[1]) if float(p[1])>0 else 1.0,p[2],p[3]) for p in v] for k,v in prox.items()},
         "mc":tot_pq.set_index("canon")["model_count"].to_dict(),"mm":tot_pq.set_index("canon")["model_mass"].to_dict(),
         "ndf":pd.read_parquet(idx_dir/"nutrient_df.parquet").set_index("nutrient_id")["df_foods"].to_dict(),
         "N":len(tot_pq),"ref":pd.read_parquet(idx_dir/"nutrient_scale.parquet").set_index("nutrient_id")["ref_amount_norm"].to_dict(),
         "maxf":args.max_foods,"force_pick":args.force_pick,
         "ams":ALLOW_MACRO_SCAN,
         "amp":ALLOW_MACRO_PROXY,
         "asp":args.allow_spices,"blocked":blocked}
         
    with mp.get_context("forkserver" if hasattr(mp,"get_context") else None).Pool(args.jobs,initializer=_init_w,initargs=(sm_path,dyn),maxtasksperchild=1) as pool:
        res=pool.map(_run_lbl,sorted(lbl2targ.keys(),key=lambda k:len(lbl2targ[k]),reverse=True))
    g_all=[g for g,_ in res if not g.empty]; u_all=[u for _,u in res if not u.empty]
    if g_all:
            out=pd.concat(g_all,ignore_index=True)
            fm=load_smart(PATH_FOOD_CSV).drop_duplicates(subset=["fdc_id"],keep="last")
            if "food_category" in fm.columns: fm=fm.drop(columns=["food_category"])
            cat=load_smart(PATH_FOOD_CATEGORY_CSV).rename(columns={"id":"food_category_id","description":"food_category"})
            fm=fm.merge(cat.assign(food_category_id=cat["food_category_id"].astype(str)),on="food_category_id",how="left")
            out=out.drop(columns=["description","data_type","food_category"],errors="ignore")
            out=round_floats(out.merge(fm[["fdc_id","description","data_type","food_category"]],left_on="representative_fdc_id",right_on=pd.to_numeric(fm["fdc_id"],errors="coerce"),how="left"),4)
            
            lead=["species","rank","description","food_category","top_nutrient_names","score","gain",
                  "n_targets","n_targets_measured","n_targets_missing","covered_targets_measured","coverage_measured_frac",
                  "top_nids","canonical_fdc_id","representative_fdc_id","data_type","norm_gain","n_improved",
                  "purity","proc_level","proc_cost","broad_cost","artifact","artifact_cost","type_bonus","breadth_factor",
                  "gain_eff","redundancy","redundancy_cost","fallback","fallback_reason","covered_targets_total",
                  "coverage_total_frac","avg_coverage_measured","min_coverage_measured","measured_nutrient_names","missing_nutrient_names"]
            out_cols=[c for c in lead if c in out]+[c for c in out if c not in lead+["_imp","fdc_id","_desc_hash"]]
            out[out_cols].to_csv(Path(args.out_prefix).with_suffix(".detailed_metrics.tsv"),sep="\t",index=False)
            
            human_map = {
                "species": "Bacterium",
                "rank": "Rank",
                "description": "Recommended_Food",
                "score": "Score",
                "top_nutrient_names": "Key_Nutrients_Provided",
                "coverage_measured_frac": "Cumulative_Diet_Coverage",
                "n_targets": "Total_Enzymatic_Targets",
                "n_targets_measured": "Databases_Quantified_Targets",
                "n_targets_missing": "Dark_Matter_Targets"
            }
            human_out = out[list(human_map.keys())].rename(columns=human_map)
            human_out.to_csv(Path(args.out_prefix).with_suffix(".clean_recommendations.tsv"),sep="\t",index=False)
            
            print(f"[*] Wrote detailed metrics to {args.out_prefix}.detailed_metrics.tsv",flush=True)
            print(f"[*] Wrote human-readable table to {args.out_prefix}.clean_recommendations.tsv",flush=True)
    if u_all: 
        unc_df = pd.concat(u_all,ignore_index=True)
        unc_df = round_floats(unc_df, 4)
        if PATH_FOODB_DIR:
            unc_df = rescue_missing_targets(unc_df, PATH_FOODB_DIR)
        unc_df.to_csv(Path(args.out_prefix).with_suffix(".uncovered_nutrients.tsv"),sep="\t",index=False)
        print(f"[*] Wrote {args.out_prefix}.uncovered_nutrients.tsv",flush=True)

if __name__=="__main__": main()