#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import difflib
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Set, Iterable

# ============================================================
# Scoring / thresholds
# ============================================================
SCORE_OVERRIDE = 200
SCORE_RHEA_CONFIRMED = 150
SCORE_EXACT_NAME = 100
SCORE_EXACT_SYNONYM = 90
SCORE_FUZZY_MATCH = 60

FUZZY_STRICT = 0.90
FUZZY_RELAXED_RHEA = 0.80  # only used in rhea-gated pool

EC_EXTRACT_RE = re.compile(r"(\d+\.\d+\.\d+\.\d+)")

# ============================================================
# Ignore ubiquitous molecules (drop entirely)
# ============================================================
IGNORE_CHEBI_IDS: Set[str] = {
    "CHEBI:57540",   # NAD
    "CHEBI:57945",   # NADH
    "CHEBI:58349",   # NADP
    "CHEBI:57783",   # NADPH
    "CHEBI:30616",   # ATP
    "CHEBI:456216",  # ADP
    "CHEBI:15377",   # water
}

IGNORE_QUERY_RE = re.compile(
    r"""(?ix)
    \b(
        nadp?h?                    # nad, nadh, nadp, nadph
      | nadp?\+                    # nad+, nadp+
      | nad\( ?\+ ?\)              # nad(+)
      | nadp\( ?\+ ?\)             # nadp(+)
      | thio[\-\s]?nadp?\+?        # thio-nad+, thio nadp+
      | atp|adp
      | h2o|water
      | abts
      | tmpd
      | tetramethyl.*phenylenediamine
    )\b
    """
)

# ============================================================
# Non-CHEBI classifiers (ONLY applied when matching fails)
# ============================================================
ASSAY_MARKERS_RE = re.compile(
    r"""(?ix)
    (4\-nitroanilide|p\-nitroanilide|nitroanilide|
     4\-methylcoumaryl|7\-amido\-4\-methylcoumarin|amc|
     azocasein|azocoll|
     \bp\-nitrophenyl\b|\b4\-nitrophenyl\b|
     benzyloxycarbonyl|z\-|boc|fmoc|t\-butyloxycarbonyl|
     dinitrophenyl|methylumbelliferyl|4-methylumbelliferyl)
    """
)

ASSAY_REAGENT_RE = re.compile(
    r"""(?ix)\b(
      abts|azino-bis|benzthiazol|benzothiazol|
      tmpd|tetramethyl.*phenylenediamine|
      tetramethylbenzidine|
      dichlorophenolindophenol|dcpip|
      remazol|brilliant\ blue|
      sulfhydryl\s+reagents
    )\b"""
)

MACROMOL_RE = re.compile(
    r"""(?ix)\b(
      protein|albumin|casein|insulin|collagen|laminin|fibrinogen|plasminogen|
      actin|tubulin|vimentin|receptor|kininogen|mucin|decorin|
      ferredoxin|flavodoxin|cytochrome|amyloid|
      ovalbumin|fetuin|asialofetuin|vitronectin|
      cadherin|catenin|synuclein|histone|
      trypsinogen|interleukin|complement|
      ezrin|protamine|phosvitin|glycoprotein|
      lysozyme|ubiquitin|caspase|procaspase|kallikrein|
      myosin|tropomyosin|
      factor\s+(?:viii|ix|x|xi|xii)|igg|thyroglobulin|
      aggrecan|osteopontin|galectin|syndecan|
      # enzyme/protein-ish words
      synthase|kinase|dehydrogenase|reductase|transferase|ligase|polymerase|
      phosphatase|isomerase|mutase|protease|peptidase|transporter|channel|
      subunit|component|complex
    )\b"""
)

PEPTIDE_3LETTER_RE = re.compile(
    r"""(?ix)
    \b(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)
    (?:-(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)){2,}\b
    """
)

FE_S_RE = re.compile(r"(?ix)\[\s*\d+\s*Fe\s*-\s*\d+\s*S\s*\]\s*[-\w]*\b|iron[-\s]?sulfur")

def classify_non_chebi_substrate(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return ""
    s_norm = norm_basic(s)

    if ASSAY_REAGENT_RE.search(s_norm):
        return "assay_reagent_no_chebi"
    if FE_S_RE.search(raw):
        return "metal_cluster_no_chebi"
    if PEPTIDE_3LETTER_RE.search(raw):
        return "peptide_no_chebi"
    if MACROMOL_RE.search(s_norm):
        return "macromolecule_no_chebi"
    if ASSAY_MARKERS_RE.search(s_norm):
        return "assay_marker_no_chebi"
    if re.search(r"(?i)\bprotein\b.*\bphospho\b", s_norm) or re.search(r"(?i)\bN\-phospho\-L\-histidine\b", s_norm):
        return "protein_modification_no_chebi"
    return ""

# ============================================================
# Patterns for moiety extraction
# ============================================================
COA_RE = re.compile(r"(?i)\b(coa|coenzyme a)\b")
TRNA_RE = re.compile(r"(?i)\btrna\b")

# ============================================================
# Normalization & helpers
# ============================================================
_CTRL_RE = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
_WS_RE = re.compile(r"\s+")
_PARENS_RE = re.compile(r"\([^)]*\)")
_BRACKETS_RE = re.compile(r"\[[^\]]*\]")

GREEK_MAP = {
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "κ": "kappa", "ω": "omega",
}

def asciiize_biochem(s: str) -> str:
    x = s or ""
    for k, v in GREEK_MAP.items():
        x = x.replace(k, v)
    x = x.replace("′", "'").replace("’", "'")
    return x

def norm_basic(s: str) -> str:
    if not s:
        return ""
    x = str(s).replace("\u00ad", "").replace("’", "'").replace("−", "-")
    x = asciiize_biochem(x)
    x = x.strip().lower()
    x = _WS_RE.sub(" ", x).strip()
    return x

def norm_no_commas(s: str) -> str:
    x = norm_basic(s).replace(",", " ")
    return _WS_RE.sub(" ", x).strip()

def norm_no_paren(s: str) -> str:
    x = norm_basic(s)
    x = _PARENS_RE.sub(" ", x)
    x = _BRACKETS_RE.sub(" ", x)
    return _WS_RE.sub(" ", x).strip()

def norm_compact(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", norm_basic(s))

def norm_ec(s: str) -> str:
    m = EC_EXTRACT_RE.search(s or "")
    return m.group(1) if m else ""

def safe_chebi_int(cid: str) -> int:
    m = re.search(r"CHEBI:(\d+)", cid or "")
    return int(m.group(1)) if m else 0

def len_bin(n: int) -> int:
    return min((n // 5) * 5, 25)

def first_alnum(s: str) -> str:
    for ch in s:
        if ch.isalnum():
            return ch.lower()
    return ""

def tsv_clean_cell(s: object) -> str:
    if s is None:
        return ""
    x = str(s)
    x = x.replace("\r", "")
    x = x.replace("\t", " ").replace("\n", " ")
    x = _CTRL_RE.sub("", x)
    x = _WS_RE.sub(" ", x).strip()
    return x

def tsv_clean_row(row: dict) -> dict:
    return {k: tsv_clean_cell(v) for k, v in (row or {}).items()}

SALT_WORDS_RE = re.compile(
    r"""(?ix)
    \b(
      sodium|potassium|lithium|ammonium|calcium|magnesium|zinc|iron|ferric|ferrous|
      chloride|bromide|iodide|fluoride|nitrate|sulfate|sulphate|phosphate|carbonate|acetate|
      salt|disodium|dipotassium|monosodium|monopotassium
    )\b
    """
)

HYDRATE_RE = re.compile(
    r"""(?ix)
    \b(
      monohydrate|dihydrate|trihydrate|tetrahydrate|pentahydrate|hexahydrate|heptahydrate|octahydrate|
      hydrate
    )\b
    """
)

def norm_salt_stripped(s: str) -> str:
    x = norm_no_paren(s)
    x = SALT_WORDS_RE.sub(" ", x)
    x = HYDRATE_RE.sub(" ", x)
    x = re.sub(r"\b\d+\s*(?:mM|uM|µM|nM|M)\b", " ", x, flags=re.I)
    return _WS_RE.sub(" ", x).strip()

# ============================================================
# Data structures
# ============================================================
@dataclass(frozen=True)
class ChebiTerm:
    chebi_id: str
    name: str

# ============================================================
# Globals (per worker)
# ============================================================
_worker_terms: Dict[str, ChebiTerm] = {}
_worker_term_texts: Dict[str, List[str]] = {}

_worker_idx_basic: Dict[str, List[str]] = {}
_worker_idx_no_commas: Dict[str, List[str]] = {}
_worker_idx_no_paren: Dict[str, List[str]] = {}
_worker_idx_compact: Dict[str, List[str]] = {}
_worker_idx_salt: Dict[str, List[str]] = {}

_worker_keys_basic: List[str] = []
_worker_keys_no_commas: List[str] = []
_worker_keys_no_paren: List[str] = []
_worker_keys_compact: List[str] = []
_worker_keys_salt: List[str] = []

_worker_fuzzy_buckets_basic: Dict[Tuple[str, int], List[str]] = {}
_worker_fuzzy_buckets_no_commas: Dict[Tuple[str, int], List[str]] = {}
_worker_fuzzy_buckets_no_paren: Dict[Tuple[str, int], List[str]] = {}
_worker_fuzzy_buckets_compact: Dict[Tuple[str, int], List[str]] = {}
_worker_fuzzy_buckets_salt: Dict[Tuple[str, int], List[str]] = {}

_worker_rhea: Dict[str, Set[str]] = {}
_worker_overrides: Dict[str, str] = {}
_worker_fuzzy_cache: Dict[Tuple[str, str, str, float], str] = {}

# ============================================================
# Override normalization / resolution (no "loose")
# ============================================================
def normalize_overrides(
    overrides_raw: Dict[str, str],
    terms: Dict[str, ChebiTerm],
    idx_basic: Dict[str, List[str]],
    idx_no_commas: Dict[str, List[str]],
    idx_no_paren: Dict[str, List[str]],
    idx_compact: Dict[str, List[str]],
    idx_salt: Dict[str, List[str]],
) -> Dict[str, str]:
    def strip_ctl(s: str) -> str:
        s = (s or "").replace("\r", "").replace("\u00ad", "")
        s = "".join(ch for ch in s if (ch == "\t" or ch == "\n" or ord(ch) >= 32))
        return s.strip()

    def looks_like_chebi_id(s: str) -> bool:
        return bool(re.match(r"^CHEBI:\d+$", s or ""))

    def pick_best(cids: List[str]) -> str:
        cids2 = [c for c in cids if c in terms]
        if not cids2:
            return ""
        return min(cids2, key=safe_chebi_int)

    def resolve_to_chebi_id(target: str) -> Tuple[str, str]:
        t = strip_ctl(target)
        if not t:
            return "", "empty"
        if looks_like_chebi_id(t):
            return (t, "ok_id") if t in terms else ("", "missing_from_terms")

        keys = [
            ("basic", idx_basic, norm_basic(t)),
            ("no_commas", idx_no_commas, norm_no_commas(t)),
            ("no_paren", idx_no_paren, norm_no_paren(t)),
            ("compact", idx_compact, norm_compact(t)),
            ("salt", idx_salt, norm_salt_stripped(t)),
        ]
        for view_name, idx, key in keys:
            if key and key in idx:
                cid = pick_best(idx[key])
                if cid:
                    return cid, f"ok_label_{view_name}"
        return "", "could_not_resolve"

    resolved: Dict[str, str] = {}
    dropped: List[Tuple[str, str, str]] = []

    for k_raw, v_raw in (overrides_raw or {}).items():
        k = strip_ctl(k_raw)
        if not k:
            continue
        cid, reason = resolve_to_chebi_id(v_raw)
        if cid:
            resolved[norm_basic(k)] = cid
        else:
            dropped.append((k_raw, v_raw, reason))

    if dropped:
        print(f"WARNING: dropped {len(dropped)} overrides that could not be resolved", file=sys.stderr)
        for k_raw, v_raw, reason in dropped[:50]:
            print(f"  override {k_raw!r} -> {v_raw!r} ({reason})", file=sys.stderr)
    return resolved

# ============================================================
# Query variants (keep conservative; avoid exploding variants)
# ============================================================
def substrate_variants(raw: str) -> List[str]:
    s = (raw or "").replace("\r", " ").replace("\n", " ").replace("\t", " ")
    s = re.sub(r"[\x00-\x1f\x7f]", " ", s)
    s = _WS_RE.sub(" ", s).strip()
    if not s:
        return []

    v: List[str] = [s]
    v.append(re.sub(r"\s*\[[^\]]+\]\s*", " ", s).strip())
    v.append(re.sub(r"\s*[\(\[]\s*[a-z0-9_+\-\s]+\s*[\)\]]\s*$", "", s, flags=re.I).strip())
    v.append(re.sub(r"^\(\s*\d+\s*\)\s*", "", s).strip())
    v.append(re.sub(r"^\d+\s+", "", s).strip())
    v.append(s.replace("(+)", "+").replace("(-)", "-").replace("NADP(+)", "NADP+").replace("NAD(+)", "NAD+"))
    v.append(re.sub(r"^\([rs]\)-", "", s, flags=re.I).strip())
    v.append(re.sub(r"^[dlrse]-", "", s, flags=re.I).strip())
    v.append(re.sub(r"(?i)\bp[-\s]?nitrophenyl\b", "p-nitrophenyl", s).strip())
    v.append(re.sub(r"(?i)\balpha[-\s]?naphthyl\b", "alpha-naphthyl", s).strip())
    v.append(norm_no_paren(s))
    
    # Halobenzaldehydes: 2-/3-/4- <-> o/m/p- and ortho/meta/para
    m = re.search(r"(?i)\b([234])[-\s]?(fluoro|chloro|bromo|iodo)\s*benzaldehyde\b", s)
    if m:
        pos = m.group(1)
        hal = m.group(2).lower()
        omp = {"2": "o", "3": "m", "4": "p"}[pos]
        v.append(f"{omp}-{hal}benzaldehyde")
        v.append(f"{omp} {hal}benzaldehyde")
        full = {"o": "ortho", "m": "meta", "p": "para"}[omp]
        v.append(f"{full}-{hal}benzaldehyde")
        v.append(f"{full} {hal}benzaldehyde")

    # methyl-1,4-benzoquinone synonyms
    if re.search(r"(?i)\bmethyl[-\s]?1,4[-\s]?benzoquinone\b", s):
        v.append("2-methyl-1,4-benzoquinone")
        v.append("toluquinone")
        v.append("tolquinone")

    # 2,5-dimethyl benzoquinone normalization
    if re.search(r"(?i)\b2,5[-\s]?dimethyl[-\s]?(?:1,4[-\s]?)?benzoquinone\b", s):
        v.append("2,5-dimethyl-1,4-benzoquinone")
        v.append("2,5-dimethyl-4-benzoquinone")


    out: List[str] = []
    seen: Set[str] = set()
    for x in v:
        x = (x or "").replace("\r", " ").replace("\n", " ").replace("\t", " ")
        x = re.sub(r"[\x00-\x1f\x7f]", " ", x)
        x = _WS_RE.sub(" ", x).strip()
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out

# ============================================================
# Worker init / fuzzy buckets
# ============================================================
def _bucketize(keys: List[str]) -> Dict[Tuple[str, int], List[str]]:
    buckets: Dict[Tuple[str, int], List[str]] = {}
    for k in keys:
        if not k:
            continue
        fc = first_alnum(k)
        if not fc:
            continue
        buckets.setdefault((fc, len_bin(len(k))), []).append(k)
    return buckets

def init_worker(
    term_texts,
    terms,
    idx_basic,
    idx_no_commas,
    idx_no_paren,
    idx_compact,
    idx_salt,
    keys_basic,
    keys_no_commas,
    keys_no_paren,
    keys_compact,
    keys_salt,
    rhea,
    overrides,
):
    global _worker_terms, _worker_term_texts
    global _worker_idx_basic, _worker_idx_no_commas, _worker_idx_no_paren, _worker_idx_compact, _worker_idx_salt
    global _worker_keys_basic, _worker_keys_no_commas, _worker_keys_no_paren, _worker_keys_compact, _worker_keys_salt
    global _worker_fuzzy_buckets_basic, _worker_fuzzy_buckets_no_commas, _worker_fuzzy_buckets_no_paren, _worker_fuzzy_buckets_compact, _worker_fuzzy_buckets_salt
    global _worker_rhea, _worker_overrides, _worker_fuzzy_cache

    _worker_term_texts = term_texts
    _worker_terms = terms

    _worker_idx_basic = idx_basic
    _worker_idx_no_commas = idx_no_commas
    _worker_idx_no_paren = idx_no_paren
    _worker_idx_compact = idx_compact
    _worker_idx_salt = idx_salt

    _worker_keys_basic = keys_basic
    _worker_keys_no_commas = keys_no_commas
    _worker_keys_no_paren = keys_no_paren
    _worker_keys_compact = keys_compact
    _worker_keys_salt = keys_salt

    _worker_fuzzy_buckets_basic = _bucketize(_worker_keys_basic)
    _worker_fuzzy_buckets_no_commas = _bucketize(_worker_keys_no_commas)
    _worker_fuzzy_buckets_no_paren = _bucketize(_worker_keys_no_paren)
    _worker_fuzzy_buckets_compact = _bucketize(_worker_keys_compact)
    _worker_fuzzy_buckets_salt = _bucketize(_worker_keys_salt)

    _worker_rhea = rhea
    _worker_overrides = overrides
    _worker_fuzzy_cache = {}

def _token_sort(s: str) -> str:
    toks = [t for t in re.split(r"\s+", s.strip()) if t]
    toks.sort()
    return " ".join(toks)

def _bucket_pool(buckets: Dict[Tuple[str, int], List[str]], q: str) -> List[str]:
    fc = first_alnum(q)
    b = len_bin(len(q))
    pool: List[str] = []
    for bb in [b, max(b - 5, 0), min(b + 5, 25)]:
        pool.extend(buckets.get((fc, bb), []))
    return pool

def _get_close_best(mode: str, view: str, q: str, pool: List[str], cutoff: float) -> str:
    ck = (mode, view, q, cutoff)
    if ck in _worker_fuzzy_cache:
        return _worker_fuzzy_cache[ck]
    m = difflib.get_close_matches(q, pool, n=1, cutoff=cutoff)
    best = m[0] if m else ""
    _worker_fuzzy_cache[ck] = best
    return best

# ============================================================
# Matching worker
# ============================================================
def match_substrate_worker(row: dict) -> dict:
    ec = norm_ec(tsv_clean_cell(row.get("ec") or row.get("ec_number") or ""))
    raw_sub = tsv_clean_cell(row.get("substrate") or row.get("sub") or "")
    norm_sub_in = tsv_clean_cell(row.get("substrate_normalized") or "")

    if not raw_sub and not norm_sub_in:
        return {"ec": ec, "sub": raw_sub, "cid": "", "cname": "", "mtype": "no_match"}

    if IGNORE_QUERY_RE.search(raw_sub) or (norm_sub_in and IGNORE_QUERY_RE.search(norm_sub_in)):
        return {"ec": ec, "sub": raw_sub, "cid": "", "cname": "", "mtype": "ignored_ubiquitous"}

    # Build candidates (conservative)
    def extract_moieties(s: str) -> List[str]:
        s0 = (s or "").strip()
        if not s0:
            return []
        out: List[str] = []

        # X-CoA / coenzyme A -> keep X
        if COA_RE.search(s0):
            x = re.sub(r"(?i)\b(coa|coenzyme a)\b", " ", s0)
            x = re.sub(r"[-–—]", " ", x)
            x = _BRACKETS_RE.sub(" ", x)
            x = _PARENS_RE.sub(" ", x)
            x = _WS_RE.sub(" ", x).strip()
            if x:
                out.append(x)

        # tRNA adducts -> remove bracket and keep remainder
        if TRNA_RE.search(s0):
            x = _BRACKETS_RE.sub(" ", s0)
            x = _PARENS_RE.sub(" ", x)
            x = _WS_RE.sub(" ", x).strip()
            if x and not re.fullmatch(r"(?i)(2'?[-\s]?phospho|phospho|ligated|trna)+", x):
                out.append(x)

        # ferrocenium salts: try stripping counterion and common spelling
        if re.search(r"(?i)\bhexafluorophosphate\b", s0):
            core = re.sub(r"(?i)\bhexafluorophosphate\b", " ", s0)
            core = _WS_RE.sub(" ", core).strip()
            if core and core != s0:
                out.append(core)
            out.append("ferrocenium")
            out.append("ferricenium")  # common misspelling

        # remove oxidation-state prefixes as candidate
        out.append(re.sub(r"^(reduced|oxidized)\s+", "", s0, flags=re.I).strip())

        uniq: List[str] = []
        seen: Set[str] = set()
        for x in out:
            x = x.strip()
            if x and x != s0 and x not in seen:
                seen.add(x)
                uniq.append(x)
        return uniq

    candidates: List[str] = []
    if raw_sub:
        candidates.append(raw_sub)
    if norm_sub_in:
        candidates.append(norm_sub_in)
    candidates.extend(extract_moieties(raw_sub))
    candidates.extend(substrate_variants(raw_sub))

    # de-dup candidates
    uniq2: List[str] = []
    seen2: Set[str] = set()
    for q in candidates:
        q = (q or "").strip()
        if q and q not in seen2:
            seen2.add(q)
            uniq2.append(q)
    candidates = uniq2

    results: Dict[str, Tuple[int, int, int]] = {}
    q_pen = norm_basic(norm_sub_in if norm_sub_in else raw_sub)

    def add_hit(cid: str, score: int, is_primary: bool = False):
        if not cid or cid in IGNORE_CHEBI_IDS:
            return
        term = _worker_terms.get(cid)
        if not term:
            return

        # rhea confirmation bumps score
        if ec and ec in _worker_rhea and cid in _worker_rhea[ec]:
            score = max(score, SCORE_RHEA_CONFIRMED)

        # mild penalties (keep from your working versions)
        name_l = term.name.lower()
        if "zwitterion" in name_l and "zwitterion" not in q_pen:
            score -= 10
        if ("pyranose" in name_l or "furanose" in name_l) and not ("pyranose" in q_pen or "furanose" in q_pen):
            score -= 5

        key = (score, 1 if is_primary else 0, -safe_chebi_int(cid))
        if cid not in results or key > results[cid]:
            results[cid] = key

    # overrides (multiple normalizations)
    for q in candidates:
        kb = norm_basic(q)
        knc = norm_no_commas(q)
        knp = norm_no_paren(q)
        kc = norm_compact(q)
        kss = norm_salt_stripped(q)
        for k in (kb, knc, knp, kc, kss):
            if k and k in _worker_overrides:
                add_hit(_worker_overrides[k], SCORE_OVERRIDE, True)

    # exact matches across views
    for q in candidates:
        kb = norm_basic(q)
        knc = norm_no_commas(q)
        knp = norm_no_paren(q)
        kc = norm_compact(q)
        kss = norm_salt_stripped(q)

        view_checks = [
            (_worker_idx_basic, kb, lambda cid: norm_basic(_worker_terms[cid].name) == kb),
            (_worker_idx_no_commas, knc, lambda cid: norm_no_commas(_worker_terms[cid].name) == knc),
            (_worker_idx_no_paren, knp, lambda cid: norm_no_paren(_worker_terms[cid].name) == knp),
            (_worker_idx_compact, kc, lambda cid: norm_compact(_worker_terms[cid].name) == kc),
            (_worker_idx_salt, kss, lambda cid: norm_salt_stripped(_worker_terms[cid].name) == kss),
        ]

        for idx, key, is_primary_fn in view_checks:
            if key and key in idx:
                for cid in idx[key]:
                    primary = is_primary_fn(cid)
                    add_hit(cid, SCORE_EXACT_NAME if primary else SCORE_EXACT_SYNONYM, primary)

    # fuzzy fallback (strict, then rhea-gated relaxed)
    if not results:
        q_src = norm_sub_in if (norm_sub_in and len(norm_sub_in) > 5) else raw_sub

        fuzzy_views = [
            ("salt", norm_salt_stripped(q_src), _worker_fuzzy_buckets_salt, _worker_idx_salt),
            ("no_paren", norm_no_paren(q_src), _worker_fuzzy_buckets_no_paren, _worker_idx_no_paren),
            ("basic", norm_basic(q_src), _worker_fuzzy_buckets_basic, _worker_idx_basic),
            ("compact", norm_compact(q_src), _worker_fuzzy_buckets_compact, _worker_idx_compact),
        ]

        for view_name, qv, buckets, idx in fuzzy_views:
            if results:
                break
            if not qv or len(qv) <= 4:
                continue
            pool = _bucket_pool(buckets, qv)
            if not pool:
                continue

            qs = [qv]
            if " " in qv and view_name in ("basic", "no_paren", "salt"):
                qs.append(_token_sort(qv))

            best_key = ""
            for qq in qs:
                bk = _get_close_best("strict", view_name, qq, pool, FUZZY_STRICT)
                if bk:
                    best_key = bk
                    break

            if best_key:
                for cid in idx.get(best_key, []):
                    add_hit(cid, SCORE_FUZZY_MATCH)

        # rhea-gated relaxed fuzzy (participant-only)
        if not results and ec and ec in _worker_rhea:
            qv = norm_no_paren(q_src)
            if qv and len(qv) > 4:
                pool_keys: List[str] = []
                key_to_cids: Dict[str, List[str]] = {}

                for cid in _worker_rhea[ec]:
                    texts = _worker_term_texts.get(cid)
                    if not texts:
                        term = _worker_terms.get(cid)
                        texts = [term.name] if term else []
                    for t in texts:
                        k = norm_no_paren(t)
                        if not k:
                            continue
                        pool_keys.append(k)
                        key_to_cids.setdefault(k, []).append(cid)

                if pool_keys:
                    best = difflib.get_close_matches(qv, pool_keys, n=1, cutoff=FUZZY_RELAXED_RHEA)
                    if best:
                        for cid in key_to_cids.get(best[0], []):
                            add_hit(cid, SCORE_EXACT_SYNONYM)

    # Only NOW apply "no_chebi" classification (does not hurt mapping)
    if not results:
        cls = classify_non_chebi_substrate(raw_sub)
        if cls:
            return {"ec": ec, "sub": raw_sub, "cid": "", "cname": "", "mtype": cls}
        return {"ec": ec, "sub": raw_sub, "cid": "", "cname": "", "mtype": "no_match"}

    best_id = max(results.keys(), key=lambda k: results[k])
    score = results[best_id][0]
    mtype = "rhea_confirmed" if score >= SCORE_RHEA_CONFIRMED else ("string_match" if score >= SCORE_EXACT_SYNONYM else "fuzzy_match")

    return {
        "ec": ec,
        "sub": raw_sub,
        "cid": best_id,
        "cname": _worker_terms[best_id].name,
        "mtype": f"{mtype}({score})",
    }

# ============================================================
# Parse ChEBI OBO
# ============================================================
def parse_chebi_obo(path: str):
    terms: Dict[str, ChebiTerm] = {}
    term_texts: Dict[str, List[str]] = {}

    idx_basic: Dict[str, List[str]] = {}
    idx_no_commas: Dict[str, List[str]] = {}
    idx_no_paren: Dict[str, List[str]] = {}
    idx_compact: Dict[str, List[str]] = {}
    idx_salt: Dict[str, List[str]] = {}

    _id_re = re.compile(r"^id:\s*(CHEBI:\d+)")
    _nm_re = re.compile(r"^name:\s*(.+)")
    _sy_re = re.compile(r'^synonym:\s*"([^"]+)"')
    _al_re = re.compile(r"^alt_id:\s*(CHEBI:\d+)")

    def add_to_indexes(cid: str, txt: str):
        b = norm_basic(txt)
        nc = norm_no_commas(txt)
        np = norm_no_paren(txt)
        c = norm_compact(txt)
        ss = norm_salt_stripped(txt)

        if b:
            idx_basic.setdefault(b, []).append(cid)
        if nc:
            idx_no_commas.setdefault(nc, []).append(cid)
        if np:
            idx_no_paren.setdefault(np, []).append(cid)
        if c:
            idx_compact.setdefault(c, []).append(cid)
        if ss:
            idx_salt.setdefault(ss, []).append(cid)

    def flush(cid: Optional[str], cnm: Optional[str], syns: List[str], alts: List[str]):
        if not cid or not cnm:
            return
        terms[cid] = ChebiTerm(cid, cnm)
        term_texts[cid] = [cnm] + syns

        add_to_indexes(cid, cnm)
        for s in syns:
            add_to_indexes(cid, s)
        for a in alts:
            add_to_indexes(cid, a)

    cur_id: Optional[str] = None
    cur_nm: Optional[str] = None
    cur_syns: List[str] = []
    cur_alts: List[str] = []
    in_term = False

    print("Parsing ChEBI OBO...", file=sys.stderr)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                if in_term:
                    flush(cur_id, cur_nm, cur_syns, cur_alts)
                cur_id, cur_nm, cur_syns, cur_alts = None, None, [], []
                in_term = True
                continue
            if line == "[Typedef]":
                if in_term:
                    flush(cur_id, cur_nm, cur_syns, cur_alts)
                in_term = False
                continue
            if not in_term:
                continue

            m = _id_re.match(line)
            if m:
                cur_id = m.group(1)
                continue
            m = _nm_re.match(line)
            if m:
                cur_nm = m.group(1)
                continue
            m = _sy_re.match(line)
            if m:
                cur_syns.append(m.group(1))
                continue
            m = _al_re.match(line)
            if m:
                cur_alts.append(m.group(1))
                continue

    if in_term:
        flush(cur_id, cur_nm, cur_syns, cur_alts)

    keys_basic = list(idx_basic.keys())
    keys_no_commas = list(idx_no_commas.keys())
    keys_no_paren = list(idx_no_paren.keys())
    keys_compact = list(idx_compact.keys())
    keys_salt = list(idx_salt.keys())

    print(f"OBO loaded terms: {len(terms):,}", file=sys.stderr)
    print(f"Index keys basic:     {len(idx_basic):,}", file=sys.stderr)
    print(f"Index keys no_commas: {len(idx_no_commas):,}", file=sys.stderr)
    print(f"Index keys no_paren:  {len(idx_no_paren):,}", file=sys.stderr)
    print(f"Index keys compact:   {len(idx_compact):,}", file=sys.stderr)
    print(f"Index keys salt:      {len(idx_salt):,}", file=sys.stderr)

    return (
        term_texts, terms,
        idx_basic, idx_no_commas, idx_no_paren, idx_compact, idx_salt,
        keys_basic, keys_no_commas, keys_no_paren, keys_compact, keys_salt
    )

# ============================================================
# Rhea map
# ============================================================
def load_rhea_metabolite_map(rhea_dir: str) -> Dict[str, Set[str]]:
    ec_to_rhea: Dict[str, Set[str]] = {}
    rhea_to_smiles: Dict[str, Set[str]] = {}
    smiles_to_chebi: Dict[str, Set[str]] = {}
    ec_to_chebi: Dict[str, Set[str]] = {}

    f_ec = os.path.join(rhea_dir, "rhea2ec.tsv")
    if os.path.exists(f_ec):
        with open(f_ec, "r", encoding="utf-8", errors="replace") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                rid = row.get("RHEA_ID") or row.get("MASTER_ID") or ""
                ec_raw = row.get("ID") or row.get("EC") or ""
                ec = norm_ec(ec_raw)
                if rid and ec:
                    ec_to_rhea.setdefault(ec, set()).add(rid)

    f_rxn = os.path.join(rhea_dir, "rhea-reaction-smiles.tsv")
    if os.path.exists(f_rxn):
        with open(f_rxn, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
                rid = parts[0]
                rhea_to_smiles.setdefault(rid, set()).update(re.split(r"\.|\s*>>\s*", parts[1]))

    f_chebi = os.path.join(rhea_dir, "rhea-chebi-smiles.tsv")
    if os.path.exists(f_chebi):
        with open(f_chebi, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    smiles_to_chebi.setdefault(parts[1], set()).add(parts[0])

    print("Building Rhea biochemical knowledge map...", file=sys.stderr)
    for ec, rids in ec_to_rhea.items():
        participants: Set[str] = set()
        for rid in rids:
            for sm in rhea_to_smiles.get(rid, []):
                if sm in smiles_to_chebi:
                    participants.update(smiles_to_chebi[sm])
        if participants:
            ec_to_chebi[ec] = participants

    return ec_to_chebi

# ============================================================
# Main (fast + low RAM)
# ============================================================
def read_rows(path: str) -> Tuple[int, List[dict]]:
    # simplest: read once; main memory use is dominated by OBO indexes anyway
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = [tsv_clean_row(r) for r in reader]
    return len(rows), rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--digest", required=True)
    ap.add_argument("--obo", required=True)
    ap.add_argument("--rhea_dir", default="rheaDB/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunksize", type=int, default=500)  # key perf knob
    args = ap.parse_args()

    (
        term_texts, terms,
        idx_basic, idx_no_commas, idx_no_paren, idx_compact, idx_salt,
        keys_basic, keys_no_commas, keys_no_paren, keys_compact, keys_salt
    ) = parse_chebi_obo(args.obo)

    rhea_map = load_rhea_metabolite_map(args.rhea_dir)

    overrides_raw = {
        # pinned / polymers
        "alpha-naphthyl acetate": "CHEBI:156230",
        "1-naphthylacetate": "CHEBI:156230",
        "hemicellulose": "CHEBI:61266",
        "xylan": "CHEBI:37166",
        "arabinoxylan": "CHEBI:22603",
        "arabinoxylans": "CHEBI:22603",
        "arabinan": "CHEBI:22590",
        "chitin": "CHEBI:17029",
        "colloidal chitin": "CHEBI:17029",
        "alpha-chitin": "CHEBI:17029",
        "glycol chitin": "CHEBI:17029",
        "microcrystalline cellulose": "CHEBI:62968",
        "phenyl alpha-glucoside": "CHEBI:91122",
        "starch": "CHEBI:28017",

        # resolve by label (will be resolved to a CHEBI id via normalize_overrides)
        "filter paper": "cellulose",
        "phosphoric acid swollen cellulose": "cellulose",
        "phosphoric acid-swollen cellulose": "cellulose",
        "cellooligomers": "cellulose",
        "sweet potato starch": "starch",
        "maize starch": "starch",
        "cassava starch": "starch",
    }

    overrides_raw.update({
        # starch variants
        "soluble starch": "starch",
        "raw starch": "starch",
        "potato starch": "starch",
        "corn starch": "starch",
        "wheat starch": "starch",
        "rice starch": "starch",
        "waxy corn starch": "starch",

        # xylan variants
        "beechwood xylan": "xylan",
        "birchwood xylan": "xylan",
        "oat spelt xylan": "xylan",
        "beechwood xylan (soluble)": "xylan",

        # arabinoxylan variants
        "wheat arabinoxylan": "arabinoxylan",
        "rye arabinoxylan": "arabinoxylan",

        # cellulose/materials
        "avicel": "microcrystalline cellulose",   # if present; else fallback resolves to cellulose if you prefer
        "avisel": "microcrystalline cellulose",
        "filter paper": "cellulose",

        # arabinan / pectin
        "1,5-alpha-l-arabinan": "arabinan",
        "sugar beet arabinan": "arabinan",
        "apple pectin": "pectin",
    })

    overrides = normalize_overrides(
        overrides_raw, terms,
        idx_basic, idx_no_commas, idx_no_paren, idx_compact, idx_salt
    )

    print(f"Reading input file: {args.digest}...", file=sys.stderr)
    total, rows = read_rows(args.digest)

    print(f"Starting {args.workers} workers for {total} rows...", file=sys.stderr)
    start = time.time()
    n_hit = 0

    with open(args.out, "w", encoding="utf-8", newline="") as out_f:
        writer = csv.DictWriter(
            out_f,
            delimiter="\t",
            fieldnames=["ec_number", "substrate", "chebi_id", "chebi_name", "match_type"],
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(
                term_texts, terms,
                idx_basic, idx_no_commas, idx_no_paren, idx_compact, idx_salt,
                keys_basic, keys_no_commas, keys_no_paren, keys_compact, keys_salt,
                rhea_map,
                overrides,
            ),
        ) as ex:
            # IMPORTANT: map + chunksize avoids 95k futures overhead
            for i, res in enumerate(ex.map(match_substrate_worker, rows, chunksize=args.chunksize), 1):
                if res.get("cid"):
                    n_hit += 1

                writer.writerow({
                    "ec_number": tsv_clean_cell(res.get("ec", "")),
                    "substrate": tsv_clean_cell(res.get("sub", "")),
                    "chebi_id": tsv_clean_cell(res.get("cid", "")),
                    "chebi_name": tsv_clean_cell(res.get("cname", "")),
                    "match_type": tsv_clean_cell(res.get("mtype", "")),
                })

                if i % 500 == 0 or i == total:
                    elapsed = time.time() - start
                    rate = i / elapsed if elapsed > 0 else 0.0
                    pct = (i / total) * 100
                    eta = (total - i) / rate if rate > 0 else 0.0
                    bar = "=" * int(pct / 2)
                    sys.stderr.write(f"\r|{bar:<50}| {pct:5.1f}% Hits: {n_hit} ETA: {int(eta)}s  ")
                    sys.stderr.flush()

    sys.stderr.write("\n")
    print(f"Final: Matched {n_hit}/{total} ({n_hit/total:.1%})", file=sys.stderr)

if __name__ == "__main__":
    main()