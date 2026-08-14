#!/usr/bin/env python3
"""
expand_nutrients_to_standard_units.py

Key Upgrades:
- REAL FDC ID BINDING: Maps specific fibers directly to "Fiber, total dietary" (1079) instead of a virtual "Fiber" node.
- REAL FDC CARB BINDING: Maps specific carbs directly to "Carbohydrate, by difference" (1005).
"""

from __future__ import annotations
import argparse
import csv
import re
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple

amino_acid_names = [
    "Alanine", "Arginine", "Asparagine", "Aspartic acid", "Cysteine", "Cystine",
    "Glutamic acid", "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine",
    "Lysine", "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine",
    "Tryptophan", "Tyrosine", "Valine",
]

_SIMPLE_SUGARS = re.compile(r"\b(sucrose|glucose|fructose|lactose|maltose|galactose|mannitol|xylitol|sorbitol|sugar)\b", re.I)

def norm_compact(s: str) -> str:
    s = (s or "").strip().lower()
    s = re.sub(r"[\(\)\[\]\{\}]", " ", s)
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

def norm_norm(s: str) -> str:
    s = (s or "").strip().lower()
    s = re.sub(r"\s+", " ", s).strip()
    return s

def strip_total_suffix(name: str) -> str:
    s = (name or "").strip()
    s = re.sub(r"\s*,\s*total\s*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s+total\s*$", "", s, flags=re.IGNORECASE)
    return s.strip()

def apply_generic_alias(name: str) -> str:
    s = (name or "").strip()
    low = s.lower().strip()
    if re.search(r"\btotal\b", low):
        cleaned = re.sub(r"\s*,\s*total\s*$", "", s, flags=re.IGNORECASE)
        cleaned = re.sub(r"\s+total\s*$", "", cleaned, flags=re.IGNORECASE)
        if cleaned != s: return cleaned.strip()

    if re.search(r"\bvitamin\s*a\b", low) and re.search(r"\b(rae|re|iu|international units)\b", low): return "Vitamin A"
    if low.startswith("vitamin c") and "," in low: return "Vitamin C"
    if "wax esters" in low: return "Wax esters"
    if "vitamin b-6" in low and ("," in low or "pyridox" in low or re.search(r"\bn\d+\b", low) or "+" in low): return "Vitamin B-6"
    
    # CRITICAL FIX: Ensure fiber generic name exactly matches FDC 1079
    if low.startswith("fiber") or "dietary fiber" in low: return "Fiber, total dietary"
    if low in {"sugars", "total sugars"} or low.startswith("sugars, total"): return "Carbohydrate, by difference"
    if "reducing sugars" in low: return "Carbohydrate, by difference"
    if low.strip() == "oligosaccharides": return "Carbohydrate, by difference"
    if "total sugar alcohols" in low or low.strip() == "sugar alcohols": return "Carbohydrate, by difference"
    
    if low.startswith(("cis-", "trans-")):
        if "beta-carotene" in low: return "Carotene, beta"
        if "lycopene" in low: return "Lycopene"
        if "lutein/zeaxanthin" in low or ("lutein" in low and "zeaxanthin" in low): return "Lutein + zeaxanthin"
    if re.search(r"\bvitamin\s*c\b", low): return "Vitamin C"
    if re.search(r"\bvitamin\s*a\b", low) and re.search(r"\b(rae|re)\b", low): return "Vitamin A"
    if re.search(r"\bvitamin\s*a\b", low) and (re.search(r"\biu\b", low) or "international units" in low): return "Vitamin A"
    if re.search(r"\b,\s*(added|intrinsic)\s*$", low): return re.sub(r"\s*,\s*(added|intrinsic)\s*$", "", s, flags=re.IGNORECASE).strip()
    if low.startswith("adjusted "): return s[len("adjusted "):].strip()
    if re.search(r"\blow\b.*\bmolecular\b.*\bfiber\b", low) or "low molecular weight dietary fiber" in low: return "Fiber, total dietary"
    if re.search(r"\bfiber\b", low) or re.search(r"\bfibre\b", low):
        if any(x in low for x in ["total dietary", "dietary fiber", "neutral detergent", "soluble", "insoluble", "aoac", "idf", "sdf", "sdfp", "sdfs", "hmwdf", "lmwdf", "high molecular weight", "low molecular weight", "crude"]):
            return "Fiber, total dietary"
    if re.fullmatch(r"(dietary\s+fiber|dietary\s+fibre|fiber|fibre)", low): return "Fiber, total dietary"
    if "nitrosamine" in low: return strip_total_suffix("Nitrosamine")
    if low.startswith("vitamin e") and "label entry" in low: return "Vitamin E"
    if re.search(r"\bvitamin\s*d\b", low) and re.search(r"\b(d2|d3|d2\s*\+\s*d3)\b", low): return "Vitamin D (D2 + D3)"
    if "25-hydroxycholecalciferol" in low: return "Vitamin D3 (cholecalciferol)"
    if "25-hydroxyergocalciferol" in low: return "Vitamin D2 (ergocalciferol)"
    if "lutein" in low and "zeaxanthin" in low: return "Lutein + Zeaxanthin"
    if low.strip() == "beta-glucans": return "Beta-glucan"
    if low.startswith("cis-"): return s[4:].strip()
    return s

@dataclass
class Nutr:
    nid: int
    name: str
    compact: str
    row: Dict[str, str]

def ensure_norm_cols(row: Dict[str, str], name: str) -> Dict[str, str]:
    out = dict(row)
    nn = (out.get("normalized_norm") or "").strip()
    nnc = (out.get("normalized_norm_no_commas") or "").strip()
    nnp = (out.get("normalized_norm_no_paren") or "").strip()
    nnx = (out.get("normalized_norm_compact") or "").strip()
    if not nn: nn = norm_norm(name)
    if not nnc: nnc = norm_norm(re.sub(r",", " ", name))
    if not nnp: nnp = norm_norm(re.sub(r"\([^)]*\)", " ", name))
    if not nnx: nnx = norm_compact(name)
    out["normalized_norm"] = nn
    out["normalized_norm_no_commas"] = nnc
    out["normalized_norm_no_paren"] = nnp
    out["normalized_norm_compact"] = nnx
    return out

def read_nutrients(fn: str) -> Dict[int, Nutr]:
    nutrients: Dict[int, Nutr] = {}
    with open(fn, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            nid = int(row["id"])
            name = row.get("name", "")
            row = ensure_norm_cols(row, name)
            comp = (row.get("normalized_norm_compact") or "").strip() or norm_compact(name)
            nutrients[nid] = Nutr(nid=nid, name=name, compact=comp, row=row)
    return nutrients

def index_by_compact(nutrients: Dict[int, Nutr]) -> Dict[str, List[int]]:
    idx: Dict[str, List[int]] = {}
    for nid, n in nutrients.items():
        idx.setdefault(n.compact, []).append(nid)
    return idx

def find_ids_by_name_regex(nutrients: Dict[int, Nutr], pattern: re.Pattern) -> List[int]:
    out: List[int] = []
    for nid, n in nutrients.items():
        if pattern.search(n.name.lower()):
            out.append(nid)
    return out

def choose_existing_or_create(
    nutrients: Dict[int, Nutr],
    compact_index: Dict[str, List[int]],
    next_id: int,
    desired_name: str,
) -> Tuple[int, str, int]:
    comp = norm_compact(desired_name)
    if comp in compact_index and compact_index[comp]:
        return compact_index[comp][0], "existing", next_id
    nid = next_id
    next_id += 1
    row = {"id": str(nid), "name": desired_name}
    row = ensure_norm_cols(row, desired_name)
    nutrients[nid] = Nutr(nid=nid, name=desired_name, compact=row["normalized_norm_compact"], row=row)
    compact_index.setdefault(comp, []).append(nid)
    return nid, "new", next_id

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nutrients", required=True, help="nutrient.normalized.tsv")
    ap.add_argument("--out_expansions", default="expansions.tsv")
    ap.add_argument("--out_units", default="nutrients.units.tsv")
    args = ap.parse_args()

    nutrients = read_nutrients(args.nutrients)
    compact_index = index_by_compact(nutrients)

    original_ids: Set[int] = set(nutrients.keys())
    next_id = (max(original_ids) if original_ids else 0) + 1

    generic_roots: Set[str] = {
        "protein", "carbohydrate", "carbohydrates", "carbohydrate by difference", "lipid", "lipids", "fat", "fiber", "fiber total dietary", "nitrosamine",
        "vitamin a", "vitamin e", "vitamin d d2 d3", "organic acids", "organic acid",
        "fatty acids total trans", "fatty acids total saturated", "fatty acids total monounsaturated",
        "fatty acids total polyunsaturated", "total lipid fat", "total fat nlea", "lutein zeaxanthin",
        "cysteine methionine", "phenylalanine tyrosine", "sugars", "total sugars", "sugars total",
        "reducing sugars", "oligosaccharides", "sugar alcohols", "total sugar alcohols",
        "other saccharides", "carbohydrate other", "phenolic acids", "polyphenols", "phenolics",
        "beta glucans", "amino acids", "minerals", "proximates", "proximate", "vitamins and other components",
        "tocotrienols", "other carotenoids", "tocopherols", "total tocopherols", "total tocotrienols",
        "folate", "folate total", "wax esters total wax", "wax esters", "tocopherols and tocotrienols",
        "glycerides", "phospholipids", "glycolipids",
    }

    discard_compacts: Set[str] = {"energy", "water", "solids", "solids soluble"}

    carb_name_pattern = re.compile(
        r"\b(starch|sucrose|glucose|fructose|lactose|maltose|amylose|amylopectin|pectin|cellulose|"
        r"hemicellulose|pentose|pentoses|pentosan|glycogen|arabinose|xylose|galactose|raffinose|"
        r"stachyose|ribose|mannose|inulin|triose|tetrose|xylitol|mannitol|sorbitol|verbascose)\b"
    )

    lipid_name_pattern = re.compile(
        r"\b(cholesterol|glyceride|triglyceride|phospholipid|glycolipid|sterol|phytosterol|ergosterol|"
        r"stigmasterol|campesterol|fatty acid|omega|dha|epa|retinol|carotene|lycopene|lutein|zeaxanthin|"
        r"phytoene|phytofluene)\b"
    )

    minerals_pattern = re.compile(
        r"\b(calcium|iron|magnesium|phosphorus|potassium|sodium|zinc|copper|selenium|manganese|iodine|"
        r"chlorine|chromium|cobalt|fluoride|molybdenum|sulfur|aluminum|antimony|arsenic|barium|beryllium|"
        r"boron|bromine|cadmium|gold|lead|lithium|mercury|nickel|rubidium|silicon|silver|strontium|tin|"
        r"titanium|vanadium)\b", re.IGNORECASE
    )

    organic_acid_name_pattern = re.compile(r"\b(acid)\b")
    fatty_acids_total_pattern = re.compile(r"\bfatty acids\b")
    fa_unit_name_pattern = re.compile(r"\b(sfa|mufa|pufa|tfa)\b|\b\d{1,2}:\d\b", re.IGNORECASE)
    vit_e_family_pattern = re.compile(r"\b(tocopherol|tocotrienol)\b", re.IGNORECASE)
    vit_a_family_pattern = re.compile(r"\b(retinol|retinal|retinoic|carotene|cryptoxanthin)\b", re.IGNORECASE)

    sugar_parent_name_pattern = re.compile(
        r"\b(sugars?\b|total sugars?\b|reducing sugars?\b|oligosaccharides?\b|sugar alcohols?\b|other saccharides?\b)\b",
        re.IGNORECASE,
    )

    expansions: List[Tuple[int, str, str, int, str, str, str, str]] = []

    def add_edge(gid: int, sid: int, source: str, rule: str):
        g = nutrients[gid]
        s = nutrients[sid]
        expansions.append((g.nid, g.name, g.compact, s.nid, s.name, s.compact, source, rule))

    def is_sugar_parent(name: str, comp: str) -> bool:
        low = (name or "").lower()
        if comp in {"sugars", "total sugars", "reducing sugars", "oligosaccharides", "sugar alcohols", "total sugar alcohols", "other saccharides", "carbohydrate other"}: return True
        if sugar_parent_name_pattern.search(low):
            if "glucose" in low or "fructose" in low or "galactose" in low or "sucrose" in low or "lactose" in low or "maltose" in low: return False
            return True
        return False

    def carb_units_from_table(exclude_ids: Set[int]) -> List[int]:
        cand = find_ids_by_name_regex(nutrients, carb_name_pattern)
        out: List[int] = []
        for cid in cand:
            if cid in exclude_ids: continue
            n = nutrients[cid]
            if n.compact in discard_compacts: continue
            if is_sugar_parent(n.name, n.compact): continue
            aliased_comp = norm_compact(apply_generic_alias(n.name))
            if aliased_comp in generic_roots: continue
            out.append(cid)
        return sorted(set(out))

    generic_ids: Set[int] = set()
    for nid, n in nutrients.items():
        if n.compact in discard_compacts: continue
        aliased = apply_generic_alias(n.name)
        aliased_comp = norm_compact(aliased)
        aliased_low = aliased.lower().strip()

        if "tocopherols and tocotrienols" in aliased_low: generic_ids.add(nid); continue
        if aliased_comp in {"amino acids", "minerals", "proximates", "proximate", "vitamins and other components", "total tocotrienols", "tocotrienols"}: generic_ids.add(nid); continue
        if aliased_comp == "other carotenoids": generic_ids.add(nid); continue
        if aliased_comp in {"tocopherols", "tocotrienols", "total tocopherols", "total tocotrienols"}: generic_ids.add(nid); continue
        if aliased_comp in {"folate", "folate total"} or aliased_low.startswith("folate, total"): generic_ids.add(nid); continue
        if aliased_low.strip() in {"glycerides", "phospholipids", "glycolipids"}: generic_ids.add(nid); continue
        
        nm = n.name.lower()
        if "cysteine" in nm and "methionine" in nm: generic_ids.add(nid); continue
        if "phenylalanine" in nm and "tyrosine" in nm: generic_ids.add(nid); continue
        if is_sugar_parent(aliased, aliased_comp): generic_ids.add(nid); continue
        if aliased_comp in generic_roots: generic_ids.add(nid); continue
        if "carbohydrate" in aliased_low or aliased_low == "protein": generic_ids.add(nid); continue
        if aliased_low in {"total lipid (fat)", "total fat (nlea)", "fat", "lipid"}: generic_ids.add(nid); continue
        if "organic acid" in aliased_low or "organic acids" in aliased_low: generic_ids.add(nid); continue
        if "fatty acids" in aliased_low and "total" in aliased_low: generic_ids.add(nid); continue
        if aliased_low.startswith("vitamin d") and ("d2 + d3" in aliased_low or ("d2" in aliased_low and "d3" in aliased_low)): generic_ids.add(nid); continue
        if "lutein" in n.name.lower() and "zeaxanthin" in n.name.lower() and "+" in n.name: generic_ids.add(nid); continue
        if aliased.strip() and aliased.strip() != n.name.strip(): generic_ids.add(nid); continue
        if "nonstarch polysaccharides" in aliased_low: generic_ids.add(nid); continue
        if aliased_low.startswith("phenolic acids") or aliased_low.startswith("polyphenols"): generic_ids.add(nid); continue

    for gid in sorted(generic_ids):
        g = nutrients[gid]
        aliased_name = apply_generic_alias(g.name)
        gnm = aliased_name.lower().strip()
        gcomp = norm_compact(aliased_name)

        if aliased_name.strip() and aliased_name.strip() != g.name.strip():
            sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, aliased_name)
            if sid != gid: add_edge(gid, sid, src, "alias_to_canonical")
            continue

        orig_low = (g.name or "").lower().strip()
        aliased_low = (aliased_name or "").lower().strip()

        if gcomp == "other carotenoids" or gnm == "other carotenoids":
            carotenoid_pattern = re.compile(r"\b(carotene|lycopene|lutein|zeaxanthin|cryptoxanthin|phytoene|phytofluene)\b", re.IGNORECASE)
            cand = find_ids_by_name_regex(nutrients, carotenoid_pattern)
            cand2 = [cid for cid in cand if cid != gid and nutrients[cid].compact not in discard_compacts and norm_compact(apply_generic_alias(nutrients[cid].name)) not in generic_roots]
            for sid in sorted(set(cand2)): add_edge(gid, sid, "existing", "other_carotenoids_to_components_existing")
            continue

        if gnm.startswith("folate") and ("total" in gnm or gcomp in {"folate", "folate total"}):
            targets = ["Folic acid", "Folate, food", "5-methyl tetrahydrofolate (5-MTHF)", "Tetrahydrofolic acid (THF)", "10-Formyl folic acid (10HCOFA)", "5-Formyltetrahydrofolic acid (5-HCOH4", "Folate, free"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "folate_total_to_forms")
            continue

        if orig_low.startswith("phenolic acids") or aliased_low.startswith("phenolic acids"):
            phenolic_acids = ["Chlorogenic acid", "Caffeic acid", "p-Coumaric acid", "Ferulic acid", "Gallic acid", "Gentisic acid", "Vanillic acid", "Cinnamic acid", "Ellagic acid", "p-Hydroxy benzoic acid"]
            for t in phenolic_acids:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "phenolic_acids_total_to_components")
            continue

        if gcomp in {"tocopherols", "total tocopherols"} or gnm in {"tocopherols", "total tocopherols"}:
            for t in ["Tocopherol, alpha", "Tocopherol, beta", "Tocopherol, gamma", "Tocopherol, delta"]:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "tocopherols_total_to_isoforms")
            continue

        if gcomp in {"tocotrienols", "total tocotrienols"} or gnm in {"tocotrienols", "total tocotrienols"}:
            for t in ["Tocotrienol, alpha", "Tocotrienol, beta", "Tocotrienol, gamma", "Tocotrienol, delta"]:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "tocotrienols_total_to_isoforms")
            continue

        if gnm == "tocopherols and tocotrienols" or gcomp == "tocopherols and tocotrienols":
            targets = ["Tocopherol, alpha", "Tocopherol, beta", "Tocopherol, gamma", "Tocopherol, delta", "Tocotrienol, alpha", "Tocotrienol, beta", "Tocotrienol, gamma", "Tocotrienol, delta"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "tocopherols_tocotrienols_to_isoforms")
            continue

        if gcomp in {"total tocotrienols", "tocotrienols"} or "tocotrienol" in gnm and "total" in gnm:
            targets = ["Tocotrienol, alpha", "Tocotrienol, beta", "Tocotrienol, gamma", "Tocotrienol, delta"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "tocotrienols_total_to_components")
            continue

        if orig_low.startswith("polyphenols") or aliased_low.startswith("polyphenols"):
            polyphenol_pattern = re.compile(r"\b(phenolic|flavonoid|flavon|flavan|anthocyan|catechin|procyanidin|isoflavone|daidzein|genistein|glycitein|daidzin|genistin|glycitin|naringenin|quercetin|kaempferol|myricetin|apigenin|luteolin|theaflav|thearubigin|theogallin|tann|coumestrol|biochanin|formononetin)\b", re.IGNORECASE)
            cand = find_ids_by_name_regex(nutrients, polyphenol_pattern)
            cand2 = [cid for cid in cand if cid != gid and not (re.search(r"\btotal\b", nutrients[cid].name.lower()) and "," in nutrients[cid].name.lower())]
            for sid in sorted(set(cand2)): add_edge(gid, sid, "existing", "polyphenols_total_to_components_existing")
            continue

        if gnm.startswith("vitamin d") and ("d2 + d3" in gnm or ("d2" in gnm and "d3" in gnm)):
            d2_id, d2_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Vitamin D2 (ergocalciferol)")
            d3_id, d3_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Vitamin D3 (cholecalciferol)")
            add_edge(gid, d2_id, d2_src, "vitamin_d_split")
            add_edge(gid, d3_id, d3_src, "vitamin_d_split")
            continue

        if gcomp == "protein" or gnm == "protein":
            for aa in amino_acid_names:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, aa)
                add_edge(gid, sid, src, "protein_to_amino_acids")
            continue

        # --- FIXED: PREBIOTIC FIBER EXPANSION TO REAL FDC IDs ---
        if gcomp == "fiber" or "fiber" in gnm or "fiber total dietary" in gcomp:
            fid, fsrc, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Fiber, total dietary")
            if fid != gid: add_edge(gid, fid, fsrc, "fiber_alias_to_fiber")
            
            targets = ["Cellulose", "Hemicellulose", "Pectin", "Beta-glucan", "Inulin", "Resistant starch"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                # CRITICAL: Link the specific fibers to the REAL FDC id (fid), not the generic 'gid'!
                add_edge(fid, sid, src, "fiber_to_specific_fibers")
            continue

        if gcomp == "nitrosamine" or "nitrosamine" in gnm:
            nid2, src2, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Nitrosamine")
            if nid2 != gid: add_edge(gid, nid2, src2, "nitrosamine_total_to_nitrosamine")
            continue

        # --- FIXED: CARBOHYDRATE SPLITTING TO REAL FDC IDs ---
        if (gcomp == "carbohydrate" or "carbohydrate" in gnm or gcomp == "carbohydrates" or is_sugar_parent(g.name, g.compact) or "carbohydrate by difference" in gcomp):
            cid, csrc, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Carbohydrate, by difference")
            if cid != gid: add_edge(gid, cid, csrc, "carb_alias")

            cand = carb_units_from_table(exclude_ids={gid, cid})
            for sid in cand:
                s_name = nutrients[sid].name
                if _SIMPLE_SUGARS.search(s_name):
                    add_edge(cid, sid, "existing", "carb_to_sugars") 
                else:
                    add_edge(cid, sid, "existing", "carb_to_complex_carbs") 
            continue

        if gcomp in {"lipid", "lipids", "fat", "total lipid fat", "total fat nlea"} or "lipid" in gnm or "total fat" in gnm:
            cand = find_ids_by_name_regex(nutrients, lipid_name_pattern)
            cand = [cid for cid in cand if cid != gid and nutrients[cid].compact not in discard_compacts]
            cand2 = [cid for cid in cand if norm_compact(apply_generic_alias(nutrients[cid].name)) not in generic_roots]
            for sid in sorted(set(cand2)):
                add_edge(gid, sid, "existing", "fat_to_fatty_acids")
            continue

        if gcomp in {"organic acids", "organic acid"} or "organic acid" in gnm:
            cand = [nid2 for nid2, n2 in nutrients.items() if nid2 != gid and n2.compact not in discard_compacts and not fatty_acids_total_pattern.search(n2.name.lower()) and organic_acid_name_pattern.search(n2.name.lower()) and "acid" in n2.name.lower()]
            for sid in sorted(set(cand)): add_edge(gid, sid, "existing", "organic_acids_to_specific_existing")
            continue

        if "fatty acids" in gnm and "total" in gnm:
            cand = [nid2 for nid2, n2 in nutrients.items() if nid2 != gid and n2.compact not in discard_compacts and fa_unit_name_pattern.search(n2.name.lower()) and "total" not in n2.name.lower()]
            for sid in sorted(set(cand)): add_edge(gid, sid, "existing", "fatty_acids_total_to_individual_existing")
            continue

        if gcomp == "vitamin e" or gnm == "vitamin e":
            cand = [cid for cid in find_ids_by_name_regex(nutrients, vit_e_family_pattern) if cid != gid and nutrients[cid].compact not in discard_compacts]
            for sid in sorted(set(cand)): add_edge(gid, sid, "existing", "vitamin_e_to_vit_e_forms_existing")
            continue

        if gcomp == "vitamin a" or gnm == "vitamin a":
            cand = [cid for cid in find_ids_by_name_regex(nutrients, vit_a_family_pattern) if cid != gid and nutrients[cid].compact not in discard_compacts]
            for sid in sorted(set(cand)): add_edge(gid, sid, "existing", "vitamin_a_to_vit_a_forms_existing")
            continue

        if "lutein" in g.name.lower() and "zeaxanthin" in g.name.lower() and "+" in g.name:
            for nid2, n2 in nutrients.items():
                if n2.name.lower().strip() == "lutein": add_edge(gid, nid2, "existing", "lutein_zeaxanthin_split")
                if n2.name.lower().strip() == "zeaxanthin": add_edge(gid, nid2, "existing", "lutein_zeaxanthin_split")
            continue

        if gcomp == "amino acids" or gnm == "amino acids":
            for aa in amino_acid_names:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, aa)
                add_edge(gid, sid, src, "amino_acids_to_components")
            continue

        if "cysteine" in g.name.lower() and "methionine" in g.name.lower():
            cys_id, cys_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Cysteine")
            met_id, met_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Methionine")
            add_edge(gid, cys_id, cys_src, "aa_combo_split")
            add_edge(gid, met_id, met_src, "aa_combo_split")
            continue

        if "phenylalanine" in g.name.lower() and "tyrosine" in g.name.lower():
            phe_id, phe_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Phenylalanine")
            tyr_id, tyr_src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, "Tyrosine")
            add_edge(gid, phe_id, phe_src, "aa_combo_split")
            add_edge(gid, tyr_id, tyr_src, "aa_combo_split")
            continue

        if "nonstarch polysaccharides" in gnm:
            targets = ["Cellulose", "Hemicellulose", "Pectin", "Beta-glucan", "Inulin"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "nonstarch_polysaccharides_to_components")
            continue

        if gnm.startswith("phenolic acids"):
            phenolic_acids = ["Chlorogenic acid", "Caffeic acid", "p-Coumaric acid", "Ferulic acid", "Gallic acid", "Gentisic acid", "Vanillic acid", "Cinnamic acid", "Ellagic acid", "p-Hydroxy benzoic acid"]
            for t in phenolic_acids:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "phenolic_acids_total_to_components")
            continue

        if gcomp == "minerals" or gnm == "minerals":
            cand = find_ids_by_name_regex(nutrients, minerals_pattern)
            cand2 = [cid for cid in cand if cid != gid and nutrients[cid].compact not in discard_compacts and norm_compact(apply_generic_alias(nutrients[cid].name)) not in generic_roots]
            for sid in sorted(set(cand2)): add_edge(gid, sid, "existing", "minerals_total_to_specific")
            continue

        if gcomp == "vitamins and other components" or gnm == "vitamins and other components":
            targets = ["Vitamin A", "Vitamin C", "Vitamin D (D2 + D3)", "Vitamin E", "Vitamin K", "Vitamin B-6", "Vitamin B-12"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "vitamins_and_other_components_to_vitamin_parents")
            continue

        if gcomp in {"proximates", "proximate"} or gnm in {"proximates", "proximate"}:
            targets = ["Protein", "Total lipid (fat)", "Carbohydrate, by difference", "Ash", "Fiber, total dietary"]
            for t in targets:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "proximates_to_macros")
            continue

        if gnm == "glycerides":
            for t in ["Triglyceride", "Triglycerides"]:
                sid, src, next_id = choose_existing_or_create(nutrients, compact_index, next_id, t)
                add_edge(gid, sid, src, "glycerides_to_triglycerides")
            continue

    added_ids: Set[int] = set(nutrients.keys()) - original_ids
    parent_of: Dict[int, Tuple[int, str, str]] = {}
    out_degree: Dict[int, int] = {}

    for (gid, _gname, gcomp, sid, _sname, _scomp, _src, rule) in expansions:
        out_degree[gid] = out_degree.get(gid, 0) + 1
        if sid not in parent_of:
            parent_of[sid] = (gid, gcomp, rule)

    generic_ids_to_drop: Set[int] = set()
    for nid, n in nutrients.items():
        if n.compact in discard_compacts: continue
        aliased_comp = norm_compact(apply_generic_alias(n.name))
        if aliased_comp in generic_roots or n.compact in generic_roots:
            generic_ids_to_drop.add(nid)

    unit_ids = sorted(set([nid for nid, n in nutrients.items() if n.compact not in discard_compacts and nid not in generic_ids_to_drop and out_degree.get(nid, 0) == 0]))

    with open(args.out_expansions, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["generic_nutrient_id", "generic_name", "generic_compact", "specific_nutrient_id", "specific_name", "specific_compact", "specific_source", "rule"])
        for row in expansions: w.writerow(list(row))

    with open(args.out_units, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["nutrient_id", "name", "normalized_norm", "normalized_norm_no_commas", "normalized_norm_no_paren", "normalized_norm_compact", "is_added_by_expansion", "parent_generic_id", "parent_generic_compact", "rule"])
        for nid in unit_ids:
            n = nutrients[nid]
            row = ensure_norm_cols(n.row, n.name)
            is_added = 1 if nid in added_ids else 0
            pgid, pgcomp, rule = parent_of[nid] if nid in parent_of else ("", "", "")
            w.writerow([nid, n.name, row["normalized_norm"], row["normalized_norm_no_commas"], row["normalized_norm_no_paren"], row["normalized_norm_compact"], is_added, pgid, pgcomp, rule])

if __name__ == "__main__":
    main()