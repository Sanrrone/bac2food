#!/usr/bin/env python3
"""
nutri2chebi_from_obo.py

Map nutrients -> ChEBI IDs by matching nutrient queries against ChEBI OBO
names + synonyms.

Important behavior
- We ONLY use ChEBI-side synonyms (from chebi.obo) as match targets.
- On the nutrient side, we DO allow controlled query expansions (vitamin D2/D3, salts, etc.)
  to increase recall, but your base query comes from the nutrient TSV columns.

Matching
- ChEBI is indexed in 4 "views": exact_norm / exact_no_commas / exact_no_paren / exact_compact
- For each nutrient, we take up to 4 base queries from your TSV:
    normalized_norm
    normalized_norm_no_commas
    normalized_norm_no_paren
    normalized_norm_compact
  Each base query is matched ONLY in its corresponding view (view-aligned).
- For each base query, we optionally generate query variants; each variant is normalized and
  matched in the SAME view as the base query (still view-aligned).

Post-filter
- Deduplicate per CHEBI ID (keep best row per CHEBI).
- Winner logic:
    if best - second >= --clear-gap : keep only best (ties allowed)
    else keep all hits within --keep-within of best.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple


# -------------------- Utilities --------------------

def _collapse_ws(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip()


# -------------------- Normalizer --------------------

_ADMIN_PATTERNS = [
    r"\bdo\s*not\s*use\b",
    r"\barchiv(?:ed|e)\b",
    r"\bdeprecated\b",
    r"\bobsolete\b",
    r"\bnot\s+for\s+use\b",
]

_META_PATTERNS = [
    r"\batwater\b",
    r"\bgeneral\s+factors?\b",
    r"\bspecific\s+factors?\b",
    r"\bby\s+difference\b",
    r"\bintrinsic\b",
]

_TRAILING_DROP = [
    r"\bequivalents?\b",
    r"\beq\.?\b",
    r"\btotal\b",
]

_LEADING_DROP = [
    r"^\s*total\s+",
    r"^\s*sum\s+of\s+",
]

_GREEK_MAP = str.maketrans({
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "κ": "kappa", "μ": "mu",
})
_GREEK_FORMS = ("alpha", "beta", "gamma", "delta")

_ARROW_MAP = {"→": "->", "⟶": "->", "⟹": "->", "⇒": "->"}

_CODELIKE_RE = re.compile(r"^[A-Z]{1,3}\d{1,3}$")  # G12, B12, etc.


def greek_form_in_text(s: str) -> Optional[str]:
    low = (s or "").lower()
    for g in _GREEK_FORMS:
        if re.search(rf"\b{g}\b", low) or re.search(rf"\b{g}-", low):
            return g
    return None


def is_codelike_label(lbl: str) -> bool:
    t = (lbl or "").strip()
    return bool(_CODELIKE_RE.fullmatch(t)) or (len(t) <= 3 and any(c.isdigit() for c in t))


@dataclass(frozen=True)
class NormalizedName:
    norm: str
    no_commas: str
    no_paren: str
    compact: str


def _remove_admin_phrases_in_parens(inner: str) -> bool:
    low = inner.lower()
    return any(re.search(p, low) for p in _ADMIN_PATTERNS)


def _paren_rewriter(match: re.Match) -> str:
    inner = match.group(1).strip()
    if not inner:
        return " "
    if _remove_admin_phrases_in_parens(inner):
        return " "
    # keep only short stereochem tokens (D/L/DL); drop other metadata in parens
    if len(inner) <= 12 and re.fullmatch(r"[A-Za-z \-]+", inner):
        if re.search(r"\b(D|L|DL)\b", inner):
            return f" ({inner}) "
        return " "
    return f" ({inner}) "


def _strip_parentheticals_selective(s: str) -> str:
    return re.sub(r"\(([^)]*)\)", _paren_rewriter, s)


def _remove_all_parentheticals(s: str) -> str:
    return _collapse_ws(re.sub(r"\([^)]*\)", " ", s))


def _basic_cleanup(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return ""

    # greek letters, arrows, dashes
    s = s.translate(_GREEK_MAP)
    for k, v in _ARROW_MAP.items():
        s = s.replace(k, v)
    s = s.replace("–", "-").replace("—", "-").replace("−", "-")

    # keep only selected parens content
    s = _strip_parentheticals_selective(s)

    # element patterns like "Phosphorus, P" -> "Phosphorus"
    s = re.sub(r"\b([A-Za-z]+),\s*[A-Z][a-z]?\b$", r"\1", s)

    # formula patterns like "Salt, NaCl" -> "Salt"
    s = re.sub(r"^(.+?),\s*([A-Z][A-Za-z0-9]{0,10})\s*$", r"\1", s)

    # leading/trailing junk
    for pat in _LEADING_DROP:
        s = re.sub(pat, "", s, flags=re.IGNORECASE)
    for pat in _META_PATTERNS:
        s = re.sub(pat, " ", s, flags=re.IGNORECASE)
    for pat in _TRAILING_DROP:
        s = re.sub(rf"\s+{pat}\s*$", " ", s, flags=re.IGNORECASE)

    # spacing normalization
    s = s.replace("/", " / ")
    s = s.replace(",", ", ")
    s = s.replace(";", " ; ")
    s = _collapse_ws(s)
    s = re.sub(r"\(\s*\)", " ", s)
    s = _collapse_ws(s)
    return s


def _normalize_core(s: str) -> str:
    s = (s or "").lower().strip()
    s = re.sub(r"\s*-\s*", "-", s)
    s = re.sub(r"\s*/\s*", "/", s)
    s = re.sub(r"\s*;\s*", "; ", s)
    s = re.sub(r"\s*,\s*", ", ", s)
    s = _collapse_ws(s)
    s = re.sub(r",\s*$", "", s)
    s = re.sub(r"^(.+?)\s*\1$", r"\1", s)  # collapse exact duplication artifacts
    return s


def _compact_last_resort(s: str) -> str:
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return _collapse_ws(s)


def normalize_text(raw: str) -> NormalizedName:
    cleaned = _basic_cleanup(raw)
    norm = _normalize_core(cleaned)
    no_commas = _collapse_ws(norm.replace(",", " "))
    no_paren = _normalize_core(_remove_all_parentheticals(cleaned))
    compact = _compact_last_resort(cleaned)
    return NormalizedName(norm=norm, no_commas=no_commas, no_paren=no_paren, compact=compact)


# -------------------- Query expansion (nutrient-side) --------------------

_GREEK_WORDS = {"alpha", "beta", "gamma", "delta"}

_QUALIFIER_DROP = [
    r"\bdetermined\b",
    r"\bequivalen(ts?|ce)\b",
    r"\bdfe\b",
    r"\bfood\b",
    r"\btotal\b",
]

ALIASES = {
    "thiamin": ["thiamine", "vitamin b1", "aneurin"],
    "thiamine": ["thiamin", "vitamin b1", "aneurin"],

    "niacin": ["nicotinic acid", "nicotinamide", "vitamin b3", "vitamin pp"],
    "vitamin b3": ["niacin", "nicotinic acid", "nicotinamide", "vitamin pp"],

    "tetrahydrofolic acid": ["tetrahydrofolate"],  # do NOT inject "thf" globally
    # remove the "thf": ... key entirely


    "5-methyl tetrahydrofolate": ["5-methyltetrahydrofolate", "5-methyl-THF", "5-MTHF", "levomefolate"],
    "folate": ["folic acid", "pteroylglutamic acid", "tetrahydrofolate", "THF"],

    "choline": ["2-hydroxyethyl(trimethyl)azanium"],
}


def alias_variants(q: str) -> List[str]:
    low = (q or "").lower()
    out: List[str] = []
    for k, vals in ALIASES.items():
        if k in low:
            out.extend(vals)
    return out


def strip_qualifiers(q: str) -> str:
    s = q or ""
    for pat in _QUALIFIER_DROP:
        s = re.sub(pat, " ", s, flags=re.IGNORECASE)
    return _collapse_ws(s)


def expand_query_variants(query: str) -> List[str]:
    """
    Returns nutrient-side query variants (NOT ChEBI synonyms; ChEBI synonyms are indexed
    from chebi.obo as match targets).

    Goal: generate a compact, non-duplicative set of plausible spellings / aliases that
    reflect common nutrient naming, stereochemistry markers, primes, fatty-acid shorthand,
    and sterol "-dienol" conventions.
    """
    # --- normalize input ---
    q = _collapse_ws((query or "").strip())
    if not q:
        return []
    q = q.replace("，", ",")
    q = re.sub(r"[,\.;:]+$", "", q).strip()
    if not q:
        return []

    low = q.lower()

    out: Set[str] = set()

    # --- always include base + basic hyphen/space variants ---
    out.add(q)
    out.add(q.replace("-", " "))
    out.add(q.replace(" ", "-"))

    # --- remove leading stereochem markers: "(+)-", "(-)-", "(+)", "(−)" etc. ---
    out.add(re.sub(r"^\(\s*[\+\-]\s*\)\s*[-]?\s*", "", q).strip())
    out.add(re.sub(r"^[\+\-]\s*[-]?\s*", "", q).strip())

    # --- qualifier-stripped variant (your drop list handles nutrition qualifiers) ---
    out.add(strip_qualifiers(q))

    # --- digits cleanup for artifacts like "niacin0" (only if present) ---
    q_nodigits = _collapse_ws(re.sub(r"\d+\b", "", q))
    if q_nodigits and q_nodigits != q:
        out.add(q_nodigits)

    # --- curated aliases (your ALIASES dict) ---
    for av in alias_variants(q):
        out.add(av)

    # --- normalize prime/apostrophe variants (3′ / 4′) ---
    if any(ch in q for ch in ("'", "’", "′")):
        q1 = q.replace("’", "'").replace("′", "'")
        out.add(q1)
        out.add(q1.replace("'", "′"))              # unicode prime
        out.add(q1.replace("'", "’"))              # curly apostrophe
        out.add(_collapse_ws(q1.replace("'", " prime ")))
        out.add(_collapse_ws(q1.replace("'", "")))  # drop primes entirely (sometimes helps)

    # --- THF disambiguation: only folate context ---
    if re.search(r"\b(thf|tetrahydrofolate|tetrahydrofolic)\b", low) and re.search(r"\bfolate|folic\b", low):
        out |= {"tetrahydrofolate", "tetrahydrofolic acid", "thf"}

    # --- targeted compounds / nutrition shorthands ---
    # --- Sterol double-bond comma variants (7,22 vs 7 22 vs 7-22) ---
    if re.search(r"\b(ergosta|ergost)[-\s].*\d,\d", low):
        out.add(re.sub(r"(\d),(\d)", r"\1 \2", low))   # 7,22 -> 7 22
        out.add(re.sub(r"(\d),(\d)", r"\1-\2", low))   # 7,22 -> 7-22
        out.add(re.sub(r",", "", low))                 # remove commas entirely

    # --- Cellulose ---
    # In many ChEBI OBOs, "cellulose" may not exist as a primary name;
    # instead you see "Cellulose, microcrystalline" or cellulose derivatives.
    # We expand to common label forms that exist as *names*.
    if low == "cellulose":
        out |= {
            "cellulose microcrystalline",
            "cellulose, microcrystalline",
            "cellulose powder",
        }

    
    # tyrosol
    if low == "tyrosol":
        out |= {
            "tyrosol",
            "4-hydroxyphenethyl alcohol",
            "4-hydroxyphenethylalcohol",
            "2-(4-hydroxyphenyl)ethanol",
            "4-hydroxyphenylethanol",
            "p-hydroxyphenethyl alcohol",
        }

    # hydroxybenzoic / ferulic
    if "p-hydroxy" in low and "benzoic" in low:
        out |= {"p-hydroxybenzoic acid", "4-hydroxybenzoic acid", "para-hydroxybenzoic acid"}
    if "ferrulic acid" in low:
        out |= {"ferulic acid", "4-hydroxy-3-methoxycinnamic acid"}

    # p-coumaric
    if "p-coumaric acid" in low or "p-coumaric" in low:
        out |= {"p-coumaric acid", "4-hydroxycinnamic acid", "para-coumaric acid"}

    # phosphotidyl misspelling
    if "phosphotidyl" in low:
        out.add(re.sub(r"phosphotidyl", "phosphatidyl", low))

    # "X from Y" constructs: keep base and source
    m_from = re.search(r"\bfrom\s+(.+)$", low)
    if m_from:
        src = m_from.group(1).strip()
        base = low.split(" from ")[0].strip()
        if src:
            out.add(src)
        if base:
            out.add(base)

    # vitamin D split + explicit
    if low.startswith("vitamin "):
        toks = q.split()
        if len(toks) >= 3:
            out.add(" ".join(toks[:2]))
            out.add(" ".join(toks[2:]))
    if "vitamin d2" in low or "ergocalciferol" in low:
        out |= {"vitamin d2", "ergocalciferol"}
    if "vitamin d3" in low or "cholecalciferol" in low:
        out |= {"vitamin d3", "cholecalciferol"}

    # heme / non-heme / salt
    if "heme" in low or "haem" in low:
        out |= {
            "heme", "haem", "heme iron", "haem iron", "iron heme", "iron haem",
            "heme b", "protoheme", "protoporphyrin ix iron", "ferriprotoporphyrin ix",
        }
    if "non heme" in low or "non-heme" in low or "nonheme" in low:
        out |= {"iron", "dietary iron", "non-heme iron", "nonheme iron"}
    if "salt" in low or "nacl" in low or "sodium chloride" in low:
        out |= {"sodium chloride", "nacl", "table salt", "salt"}

    # tocopherol/tocotrienol -> vitamin E fallback (specific alpha/beta/etc handled elsewhere via scoring)
    if re.search(r"\b(tocopherol|tocotrienol)\b", low):
        out.add("vitamin e")

    # trailing greek word flip: "tocopherol beta" -> "beta tocopherol" etc.
    toks = q.split()
    if len(toks) >= 2 and toks[-1].lower() in _GREEK_WORDS:
        greek = toks[-1]
        head = " ".join(toks[:-1])
        out.add(f"{greek} {head}")
        out.add(f"{greek}-{head.replace(' ', '-')}")
        out.add(f"{greek} {head.replace('-', ' ')}")

    # catechin / gallate patterns
    if "catechin" in low and ("gallate" in low or "galloyl" in low):
        out |= {
            "catechin 3-gallate",
            "catechin-3-gallate",
            "(+)-catechin 3-gallate",
            "(+)-catechin-3-gallate",
            "catechin 3-o-gallate",
            "catechin-3-o-gallate",
            "(+)-catechin-3-o-gallate",
            "catechin gallate",
        }

    # --- Pentosan (often represented as pentosan polysaccharide / xylan-like) ---
    if low == "pentosan":
        out |= {
            "pentosan",
            "pentosan polysaccharide",
            "pentosan polysaccharides",
            "pentose polymer",
            "xylan",          # common “pentosan” class in nutrition contexts
            "arabinoxylan",   # common cereal pentosan
        }

    # gallocatechin (handles "gallo catechin" phrasing)
    if re.search(r"\bgallo\s*catechin\b", low) or "gallocatechin" in low:
        out |= {"gallocatechin", "(+)-gallocatechin", "gallo catechin", "gallo-catechin"}

    # epicatechin family
    if low == "epigallocatechin":
        out |= {"epigallocatechin", "(-)-epigallocatechin", "egc"}
    if low == "epicatechin":
        out |= {"epicatechin", "(-)-epicatechin", "ec"}
    if low in {"epicatechin-3-gallate", "epicatechin 3 gallate"}:
        out |= {
            "epicatechin-3-gallate",
            "epicatechin 3-o-gallate",
            "epicatechin 3 o gallate",
            "(-)-epicatechin-3-o-gallate",
            "epicatechin gallate",
            "ecg",
        }

    # isosakuranetin / chrysoeriol (include systematic variants + prime variants)
    if low == "isosakuranetin":
        out |= {
            "isosakuranetin",
            "4′-methoxynaringenin",
            "4'-methoxynaringenin",
            "5,7-dihydroxy-4′-methoxyflavanone",
            "5,7-dihydroxy-4'-methoxyflavanone",
            "5,7-dihydroxy-4-methoxyflavanone",
            "4'-methoxyflavanone",
        }
    if low == "chrysoeriol":
        out |= {
            "chrysoeriol",
            "3′-methoxyluteolin",
            "3'-methoxyluteolin",
            "3 prime methoxyluteolin",
            "3-methoxyluteolin",
            "3',4',5,7-tetrahydroxy-3-methoxyflavone",
            "3′,4′,5,7-tetrahydroxy-3-methoxyflavone",
            "3'-methoxy-4',5,7-trihydroxyflavone",
        }

    # daidzin / genistin
    if low == "daidzin":
        out |= {
            "daidzin",
            "daidzein 7-o-glucoside",
            "daidzein-7-o-glucoside",
            "daidzein 7-glucoside",
            "daidzein-7-glucoside",
        }
    if low == "genistin":
        out |= {
            "genistin",
            "genistein 7-o-glucoside",
            "genistein-7-o-glucoside",
            "genistein 7-glucoside",
            "genistein-7-glucoside",
        }

    # stigmastadiene: add common canonical stigmasta-* forms + class fallback
    if low in {"stigmastadiene", "stigmasta diene"}:
        out |= {
            "stigmastadiene",
            "stigmasta-diene",
            "stigmasta  diene",
            "stigmasta-5,22-diene",
            "stigmasta-7,22-diene",
            "stigmasta-5,24-diene",
            "stigmasta-7,24-diene",
            # if you accept class-level mapping when underspecified:
            "stigmastane",
            "stigmastane sterol",
        }

    # delta-7 stigmastenol patterns
    if low in {"delta-7-stigmastenol", "delta 7 stigmastenol", "delta7-stigmastenol"}:
        out |= {
            "delta-7-stigmastenol",
            "Δ7-stigmastenol",
            "delta7-stigmastenol",
            "stigmast-7-enol",
            "stigmasta-7-enol",
            "stigmast-7-en-3-ol",
            "stigmast-7-en-3beta-ol",
            "stigmast-7-en-3β-ol",
        }
        # delta symbol normalization
        out.add(re.sub(r"\bdelta\b", "Δ", q, flags=re.IGNORECASE))
        out.add(q.replace("Δ", "delta-").replace("Δ", "delta "))

    # ergosta-*-dienol -> ergosta-*-dien-3β-ol (ChEBI style)
    if re.search(r"\b(ergosta|ergost)[-\s].*dienol\b", low):
        out.add(re.sub(r"dienol\b", "dien-3-ol", low))
        out.add(re.sub(r"dienol\b", "dien-3beta-ol", low))
        out.add(re.sub(r"dienol\b", "dien-3β-ol", low))
    if low == "ergosta-7,22-dienol":
        out |= {"ergosta-7,22-dien-3-ol", "ergosta-7,22-dien-3beta-ol", "ergosta-7,22-dien-3β-ol"}
    if low == "ergosta-5,7-dienol":
        out |= {"ergosta-5,7-dien-3-ol", "ergosta-5,7-dien-3beta-ol", "ergosta-5,7-dien-3β-ol"}

    # ergosta-7-enol style
    if low in {"ergosta-7-enol", "ergosta 7 enol", "ergosta7-enol"}:
        out |= {
            "ergosta-7-enol",
            "ergost-7-enol",
            "ergosta-7-en-3-ol",
            "ergost-7-en-3-ol",
            "ergost-7-en-3beta-ol",
            "ergost-7-en-3β-ol",
        }

    # thearubigins: known absent in your chebi.obo, but keep harmless variants
    if low == "thearubigins":
        out |= {"thearubigin", "thearubigins", "tea pigments", "black tea pigments"}

    # --- fatty-acid shorthand (supports c/t flag) ---
    m_fa = re.match(
        r"^(SFA|MUFA|PUFA|TFA)\s+(\d{1,2}):(\d)(?:\s+n-(\d))?(?:\s*\(([^)]+)\))?(?:\s+([ct]))?$",
        q.strip(),
        flags=re.IGNORECASE,
    )
    if m_fa:
        cls = m_fa.group(1).upper()
        c = int(m_fa.group(2))
        db = int(m_fa.group(3))
        paren = (m_fa.group(5) or "").strip().upper()
        ct = (m_fa.group(6) or "").lower()

        if "DHA" in paren:
            out |= {"docosahexaenoic acid", "DHA"}
        if "EPA" in paren:
            out |= {"eicosapentaenoic acid", "EPA"}
        if "DPA" in paren:
            out |= {"docosapentaenoic acid", "DPA"}

        if cls == "SFA" and db == 0:
            sfa = {
                4: ["butyric acid", "butanoic acid"],
                5: ["valeric acid", "pentanoic acid"],
                6: ["caproic acid", "hexanoic acid"],
                7: ["enanthic acid", "heptanoic acid"],
                8: ["caprylic acid", "octanoic acid"],
                9: ["pelargonic acid", "nonanoic acid"],
                10: ["capric acid", "decanoic acid"],
                11: ["undecylic acid", "undecanoic acid"],
                12: ["lauric acid", "dodecanoic acid"],
                13: ["tridecylic acid", "tridecanoic acid"],
                14: ["myristic acid", "tetradecanoic acid"],
                15: ["pentadecylic acid", "pentadecanoic acid"],
                16: ["palmitic acid", "hexadecanoic acid"],
                17: ["margaric acid", "heptadecanoic acid"],
                18: ["stearic acid", "octadecanoic acid"],
                19: ["nonadecylic acid", "nonadecanoic acid"],
                20: ["arachidic acid", "eicosanoic acid"],
                21: ["heneicosanoic acid"],
                22: ["behenic acid", "docosanoic acid"],
                23: ["tricosanoic acid"],
                24: ["lignoceric acid", "tetracosanoic acid"],
            }
            for nm in sfa.get(c, []):
                out.add(nm)

        # generic systematic class-like names for unsaturated
        chain = {
            10: "dec",
            11: "undec",
            12: "dodec",
            13: "tridec",
            14: "tetradec",
            15: "pentadec",
            16: "hexadec",
            17: "heptadec",
            18: "octadec",
            19: "nonadec",
            20: "eicos",
            21: "heneicos",
            22: "docos",
            23: "tricos",
            24: "tetracos",
        }.get(c)

        unsat = {
            1: "enoic",
            2: "adienoic",
            3: "atrienoic",
            4: "atetraenoic",
            5: "apentaenoic",
            6: "ahexaenoic",
        }.get(db)

        if chain and unsat and db > 0:
            base = f"{chain}{unsat} acid"
            out.add(base)
            out.add(f"{c}:{db}")
            if ct == "c":
                out |= {"cis", f"cis-{base}"}
            elif ct == "t":
                out |= {"trans", f"trans-{base}"}

        if cls == "PUFA" and c == 20 and db == 4:
            out.add("arachidonic acid")

    # plural -> singular fallback
    if low.endswith("s") and len(low) > 4:
        out.add(low[:-1])

    # p- prefix variants (general)
    if re.match(r"^p[-\s]", low):
        out.add(re.sub(r"^p[-\s]+", "para ", low))
        out.add(re.sub(r"^p[-\s]+", "4-", low))
        out.add(re.sub(r"^p[-\s]+", "4 ", low))

    # return stable-ish order (set -> list)
    return [x for x in out if x]


# -------------------- TSV helpers --------------------

def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        return [{k: (v if v is not None else "") for k, v in row.items()} for row in r]


def write_tsv(path: str, fieldnames: List[str], rows: Iterable[Dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


# -------------------- OBO parsing (ChEBI) --------------------

_OBO_TERM_START = re.compile(r"^\[Term\]\s*$")
_OBO_KEYVAL = re.compile(r"^([a-zA-Z_]+):\s*(.*)$")


def _parse_synonym_line(v: str) -> Optional[Tuple[str, str]]:
    # synonym: "Cellulose" RELATED [kegg.compound]
    m = re.match(r'^"(.+?)"\s+(EXACT|BROAD|NARROW|RELATED)\b', v)
    if not m:
        return None
    return (m.group(1), m.group(2))



def parse_chebi_obo(path: str) -> Dict[str, Dict[str, List]]:
    terms: Dict[str, Dict[str, List]] = {}
    cur_id: Optional[str] = None
    cur_name: Optional[str] = None
    cur_syns: List[Tuple[str, str]] = []
    cur_subsets: List[str] = []
    in_term = False

    def flush() -> None:
        nonlocal cur_id, cur_name, cur_syns, cur_subsets
        if cur_id:
            terms[cur_id] = {
                "name": ([cur_name] if cur_name else []),
                "syn": list(cur_syns),  # list of (text, scope)
                "subset": list(cur_subsets),  # e.g. ["3:STAR"] - ChEBI's curation status
            }
        cur_id = None
        cur_name = None
        cur_syns = []
        cur_subsets = []

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if _OBO_TERM_START.match(line):
                if in_term:
                    flush()
                in_term = True
                continue
            if not in_term or not line.strip():
                continue

            m = _OBO_KEYVAL.match(line)
            if not m:
                continue
            k, v = m.group(1), m.group(2)

            if k == "id" and v.startswith("CHEBI:"):
                cur_id = v.strip()
            elif k == "name":
                cur_name = v.strip()
            elif k == "synonym":
                parsed = _parse_synonym_line(v.strip())
                if parsed:
                    cur_syns.append(parsed)
            elif k == "subset":
                cur_subsets.append(v.strip())

    if in_term:
        flush()
    return terms



# -------------------- Matching + Scoring --------------------

MatchKey = Tuple[str, str]  # (match_view, key)
SymbolKey = str             # exact case-sensitive symbol, e.g. "Co"


def unique_preserve(seq: Iterable[Tuple[str, str, str]]) -> List[Tuple[str, str, str]]:
    seen: Set[Tuple[str, str]] = set()
    out: List[Tuple[str, str, str]] = []
    for cid, lbl, term in seq:
        key = (cid, term)
        if key in seen:
            continue
        seen.add(key)
        out.append((cid, lbl, term))
    return out


VIEW_WEIGHT = {
    "exact_norm": 40,
    "exact_no_commas": 35,
    "exact_no_paren": 30,
    "exact_compact": 20,
    "symbol_exact": 22,
    "element_symbol": 50,      # deterministic: outranks every string-similarity view
}

_ELEMENT_SYMBOL_RE = re.compile(r",\s*([A-Z][a-z]?)\s*$")

# --- INFOODS MINERAL TAGNAMES ----------------------------------------------------------
# Several sources name their mineral columns with the INFOODS tagname plus a unit, e.g.
# "CS(mcg)", "MG(mg)", "K(mg)". Those tagnames are element symbols UPPERCASED, which makes
# them collide head-on with the IUPAC one-letter amino-acid codes that ChEBI uses for
# peptide synonyms. Left to the ordinary name matcher this produced, in the shipped map:
#     K(mg)   -> lysine / Lys-Lys      P(mg)  -> L-proline / poly(butyl acrylate)
#     CS(mcg) -> Cys-Ser               AS(mcg)-> Ala-Ser / artesunate
#     TI(mcg) -> Thr-Ile               MG(mg) -> Met-Gly / monoacylglycerols
#     S(mg)   -> Ser-Ser / STRANGE QUARK        B(mcg) -> BOTTOM QUARK
# i.e. minerals measured in micrograms were being handed the enzyme repertoires of
# peptides, polymers, drugs and subatomic particles. Recognising the tagname and forcing
# the element reading is the fix; doing it here keeps the correction in the script that
# creates the link rather than in the emitted TSV.
_ELEMENT_SYMBOLS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu",
}
# INFOODS codes that are NOT the bare element symbol. "I" would be read as a roman numeral,
# so iodine is tagged ID; that is the only such alias in the sources ingested here.
_INFOODS_ALIAS = {"ID": "I"}
_TAGNAME_UNIT_RE = re.compile(r"^\s*([A-Za-z]{1,3})\s*\([^)]{1,6}\)\s*$")

# The general "<code>(<unit>)" column name, of which the element columns are one case.
# Sources also use it for fatty-acid chain lengths ("C10 (%T)" = C10 as % of total fatty
# acids), vitamins ("PP-(mg)" = niacin) and classes ("WAX(mg)", "GABA(mg)").
_TAGNAME_ANY_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9]{0,5})\s*[-–]?\s*\([^)]{1,8}\)\s*$")


def uninformative_tagname(raw_name: str) -> str:
    """The code of a `<code>(<unit>)` column too short to be matched by name, else "".

    A one- or two-letter code, or any code carrying a digit, conveys nothing a string
    matcher can use, and ChEBI is full of one-letter synonyms for it to collide with. Left
    to the ordinary matcher these produced, in the shipped map:
        C10 (%T) .. C24 (%T)  ->  cysteine / cytosine / calciol   (17 columns, ~150 EC each)
        PP-(mg)   [= niacin]  ->  proton / Pro-Pro / L-proline
        OA(g)     [= oxalic]  ->  oxolinic acid
    i.e. fatty-acid chain-length columns were handed the enzyme repertoire of an amino acid
    and a nucleobase. There is no evidence in a two-character code to decide between those
    readings, so the honest output is NO ChEBI id: an unmatched nutrient is a known gap,
    whereas a wrongly matched one silently propagates into every food score.

    Longer codes (WAX, GABA, VITA) are left to the normal matcher, which resolves them on
    their own merits. Element columns are handled earlier and never reach this test.
    """
    m = _TAGNAME_ANY_RE.match((raw_name or "").strip())
    if not m:
        return ""
    code = m.group(1)
    return code if (len(code) <= 2 or any(c.isdigit() for c in code)) else ""


def extract_element_symbol(raw_name: str) -> str:
    """Element symbol from the FDC style "Calcium, Ca", or "" if absent."""
    m = _ELEMENT_SYMBOL_RE.search((raw_name or "").strip())
    return m.group(1) if m else ""


def extract_tagname_element(raw_name: str) -> str:
    """Element symbol from the INFOODS style "CA(mg)" / "CS(mcg)", or "" if not one.

    Deliberately separate from extract_element_symbol. The FDC style spells the element
    out ("Calcium, Ca"), so the ordinary name matcher already resolves it correctly and
    with better specificity - "Fluoride, F" reaches the fluoride ION, which is the species
    actually measured, where an atom-only rule would coarsen it to "fluorine atom". Only
    the bare tagname form is ambiguous enough to need overriding, so only it is locked.
    """
    m = _TAGNAME_UNIT_RE.match((raw_name or "").strip())
    if not m:
        return ""
    tag = m.group(1).upper()
    sym = _INFOODS_ALIAS.get(tag, tag.capitalize())
    return sym if sym in _ELEMENT_SYMBOLS else ""


def element_atom_hits(sym: str, sym_idx, label_by_id) -> List[Tuple[str, str, str]]:
    """The ChEBI neutral-atom entity for an element symbol, as [(chebi_id, label, term)].

    The symbol index alone is not enough to identify an element: it is keyed on 1-2 letter
    ChEBI synonyms, and "K" is a synonym of BOTH potassium and lysine, "S" of sulfur and
    serine, "V" of vanadium and valine. ChEBI does, however, label every neutral atom
    "<name> atom" ("potassium atom", "caesium atom"), and no peptide or polymer carries
    that suffix - so the suffix is what disambiguates. Returns [] when the symbol has no
    atom entity, in which case the caller falls back to ordinary name matching.
    """
    out: List[Tuple[str, str, str]] = []
    for cid, _lbl, term in sym_idx.get(sym, []):
        lbl = label_by_id.get(cid, _lbl) or _lbl
        if (lbl or "").strip().lower().endswith(" atom"):
            out.append((cid, lbl, term))
    return out


def looks_like_element_query(q_norm: str) -> bool:
    q = (q_norm or "").strip().lower()
    return bool(re.fullmatch(r"[a-z]{3,20}", q))


def complexity_penalty(term_text: str) -> int:
    t = (term_text or "")
    p = 0
    if "(" in t or ")" in t:
        p += 20
    if re.search(r"\(\s*\d+\s*[\+\-]\s*\)", t):
        p += 25
    if re.search(r"\b[IVX]{1,4}\b", t):
        p += 10
    if re.search(r"\b\d+\b", t):
        p += 10
    if "+" in t or "-" in t:
        p += 10
    if "," in t or ";" in t:
        p += 4
    if len(t) >= 25:
        p += min(15, ((len(t) - 25) // 10) * 3)
    return p


def parse_source_from_term(chebi_term: str) -> str:
    if chebi_term.endswith("[name]"):
        return "name"
    m = re.search(r"\[(synonym:(EXACT|BROAD|NARROW|RELATED))\]$", chebi_term)
    if m:
        return m.group(1)  # e.g. synonym:RELATED
    return "obo"



def score_hit(nutrient_query: str, view: str, chebi_label: str, chebi_term: str) -> int:
    base = VIEW_WEIGHT.get(view, 0)

    src_bonus = 0
    if chebi_term.endswith("[name]"):
        src_bonus = 10
    else:
        m = re.search(r"\[synonym:(EXACT|BROAD|NARROW|RELATED)\]$", chebi_term)
        if m:
            scope = m.group(1)
            src_bonus = {
                "EXACT": 7,
                "BROAD": 4,
                "NARROW": 4,
                "RELATED": 1,
            }[scope]


    s = base + src_bonus

    n_q = normalize_text(nutrient_query).norm
    n_lbl = normalize_text(chebi_label).norm if chebi_label else ""

    qlow = n_q.lower()
    lbl_low = (chebi_label or "").lower()
    term_low = (chebi_term or "").lower()

    # Greek form agreement (prevents beta queries matching alpha candidates)
    q_form = greek_form_in_text(nutrient_query)
    if q_form:
        cand_form = greek_form_in_text(chebi_label) or greek_form_in_text(chebi_term)
        if cand_form and cand_form != q_form:
            s -= 80
        if not cand_form:
            s -= 20

    # Substituted folate specificity: if query indicates 5-methyl folate, prefer matching candidates
    if re.search(r"\b5[-\s]?methyl\b", qlow) or "5-mthf" in qlow or "methyltetrahydrofolate" in qlow:
        # boost candidates that also indicate 5-methyl tetrahydrofolate
        if re.search(r"\b5[-\s]?methyl\b", lbl_low) or "methyltetrahydrofolate" in lbl_low:
            s += 40
        else:
            # penalize generic folate/THF hits
            if any(t in lbl_low for t in ["folic acid", "folate", "tetrahydrofolate", "tetrahydrofolic"]):
                s -= 25


    # Vitamin query heuristic: avoid code-like entries and prefer "vitamin" in label/term
    if "vitamin" in qlow:
        if "vitamin" in lbl_low or "vitamin" in term_low:
            s += 25
        else:
            s -= 25
        if is_codelike_label(chebi_label):
            s -= 35

    # exact normalized label match boost
    if n_lbl and n_lbl == n_q:
        s += 25

    # element-only logic
    if looks_like_element_query(n_q):
        if " atom" in lbl_low or lbl_low.endswith(" atom"):
            s += 30

        # ChEBI names every curated element entity "<element> atom", and the ion where one
        # is measured - that naming convention is the same one element_atom_hits() relies
        # on to tell potassium from lysine. For an element query such a label IS the exact
        # match, so it must collect what a bare-name candidate collects from the
        # exact-label (+25) and single-token (+15) rules; otherwise a record simply named
        # "copper" outscores "copper atom" 90 to 77 on those two bonuses alone.
        # ChEBI release 253 imported exactly such records from ChEMBL. CHEBI:748257
        # "copper" is_a organic molecular entity - wrong for a metal, no element or metal
        # parentage, and its only other edges are has_role, which the nutrient->EC walk
        # excludes on purpose. Letting those win severed copper, chromium, gold, arsenic
        # and antimony from their enzymes while the map still looked fully populated.
        # "Fluoride, F" keeps reaching the fluoride ION, which is the species actually
        # measured, because " ion" is credited the same way.
        if n_lbl in (f"{n_q} atom", f"{n_q} ion"):
            s += 40

        if n_lbl and n_lbl != n_q and not n_lbl.endswith(" atom") and not n_lbl.endswith(" ion"):
            if n_q not in n_lbl:
                s -= 60

        if n_lbl and n_q not in n_lbl:
            s -= 25

        s -= complexity_penalty(chebi_label or chebi_term)

        if re.fullmatch(r"[A-Za-z]{3,20}", (chebi_label or "").strip()):
            s += 15

    if n_q == "cellulose":
        bad = {"acetate", "propionate", "phthalate", "inhibitor", "synthesis"}
        cand_text = f"{chebi_label} {chebi_term}".lower()
        if any(w in cand_text for w in bad):
            s -= 80


    return s


def build_chebi_indices(
    chebi_terms: Dict[str, Dict[str, List[str]]]
) -> Tuple[
    Dict[MatchKey, List[Tuple[str, str, str]]],
    Dict[SymbolKey, List[Tuple[str, str, str]]],
    Dict[str, str],
]:
    idx: Dict[MatchKey, List[Tuple[str, str, str]]] = defaultdict(list)
    sym_idx: Dict[SymbolKey, List[Tuple[str, str, str]]] = defaultdict(list)
    label_by_id: Dict[str, str] = {}

    def add_term(cid: str, term_text: str, source: str) -> None:
        t0 = (term_text or "").strip()
        if not t0:
            return

        # drop long ALLCAPS noise tokens
        if re.fullmatch(r"[A-Z0-9]+", t0) and len(t0) > 3:
            return

        tagged = f"{t0} [{source}]"

        # case-sensitive symbol index only for 1–2 letters
        if re.fullmatch(r"[A-Za-z]{1,2}", t0):
            sym_idx[t0].append((cid, label_by_id.get(cid, ""), tagged))

        n = normalize_text(t0)
        for view, key in [
            ("exact_norm", n.norm),
            ("exact_no_commas", n.no_commas),
            ("exact_no_paren", n.no_paren),
            ("exact_compact", n.compact),
        ]:
            if key:
                idx[(view, key)].append((cid, label_by_id.get(cid, ""), tagged))

    for cid, payload in chebi_terms.items():
        name_list = payload.get("name", [])
        syn_list = payload.get("syn", [])

        if name_list:
            label_by_id[cid] = name_list[0]

        for nm in name_list:
            add_term(cid, nm, "name")
        for sy, scope in syn_list:
            add_term(cid, sy, f"synonym:{scope}")



    return idx, sym_idx, label_by_id


# -------------------- Main --------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nutrients", required=True, help="nutrients.standard.tsv (tab-separated)")
    ap.add_argument("--chebi-obo", required=True, help="chebi.obo")
    ap.add_argument("--chebi-table", default=None, help="optional chebi_table.tsv for label backfill")
    ap.add_argument("--out", required=True, help="output nutrient_to_chebi.tsv (long format)")
    ap.add_argument("--min-score", type=int, default=30, help="minimum score to keep a hit (default: 30)")
    ap.add_argument("--clear-gap", type=int, default=20,
                    help="if best score exceeds 2nd best by >= this, keep only best (default: 20)")
    ap.add_argument("--keep-within", type=int, default=8,
                    help="otherwise keep hits within this many points of best (default: 8)")
    ap.add_argument("--fallback-bar", type=int, default=40,
                    help="run the fallback rescue pass unless the primary pass already "
                         "produced a hit scoring at least this (default: 40)")
    args = ap.parse_args()

    nutrients = read_tsv(args.nutrients)
    chebi_terms = parse_chebi_obo(args.chebi_obo)
    idx, sym_idx, label_by_id = build_chebi_indices(chebi_terms)

    # ChEBI curation status. "3:STAR" is manually curated by the ChEBI team; "2:STAR" is
    # an unvetted third-party submission. See the curation-status rescue below.
    curated_ids = {cid for cid, payload in chebi_terms.items()
                   if "3:STAR" in (payload.get("subset") or [])}

    def is_curated(cid: str) -> bool:
        return cid in curated_ids

    # optional label backfill
    if args.chebi_table:
        chebi_table = read_tsv(args.chebi_table)
        for row in chebi_table:
            cid = (row.get("chebi_id", "") or "").strip()
            lbl = (row.get("lbl", "") or "").strip()
            if cid and lbl and (cid not in label_by_id or not label_by_id[cid] or label_by_id[cid] == "NA"):
                label_by_id[cid] = lbl

    fieldnames = [
        "nutrient_id", "nutrient_name", "query",
        "match_view", "match_source",
        "chebi_id", "chebi_label", "chebi_term", "score",
    ]
    all_rows: List[Dict[str, str]] = []
    n_element_locked = 0        # element columns settled by the symbol pass (see 1b)
    n_code_skipped = 0          # <code>(unit) columns too short to match by name
    n_star_rescued = 0          # nutrients where a curated 3:STAR term was added back

    VIEW_TO_COL = {
        "exact_norm": "normalized_norm",
        "exact_no_commas": "normalized_norm_no_commas",
        "exact_no_paren": "normalized_norm_no_paren",
        "exact_compact": "normalized_norm_compact",
    }
    VIEWS = ["exact_norm", "exact_no_commas", "exact_no_paren", "exact_compact"]

    def _keys_for_view(view: str, query_variant: str) -> str:
        nn = normalize_text(query_variant)
        if view == "exact_norm":
            return nn.norm
        if view == "exact_no_commas":
            return nn.no_commas
        if view == "exact_no_paren":
            return nn.no_paren
        if view == "exact_compact":
            return nn.compact
        return nn.norm

    for r in nutrients:
        nid = (r.get("nutrient_id") or r.get("id") or "").strip()
        nname = (r.get("name") or "").strip()

        # stable output query label
        out_query = (r.get("normalized_norm") or "").strip()
        if not out_query:
            out_query = normalize_text(nname).norm

        hits: List[Tuple[str, str, str, str, int]] = []  # (view, chebi_id, chebi_label, chebi_term, score)

        # 1) ELEMENT COLUMNS. Resolved by the atom-suffix rule, and once resolved the
        # reading is SETTLED: the name passes are skipped entirely. Running them as well
        # is what let "CS(mcg)" pick up Cys-Ser alongside caesium - both score, both
        # survive --keep-within, and the peptide then inherits caesium's enzyme set. There
        # is no reading of a mineral column under which a dipeptide is a candidate.
        # Scored 100 rather than through score_hit(): this is a deterministic
        # identification from the column name, not a fuzzy string similarity, and
        # score_hit's symbol_exact view returns 23, below the default --min-score of 30,
        # which is why the symbol pass had never once fired.
        tag_sym = extract_tagname_element(nname)
        el_hits = element_atom_hits(tag_sym, sym_idx, label_by_id) if tag_sym else []
        bad_code = "" if tag_sym else uninformative_tagname(nname)
        if el_hits:
            for cid, lbl, term in unique_preserve(el_hits):
                hits.append(("element_symbol", cid, lbl, term, 100))
            n_element_locked += 1
        elif bad_code:
            n_code_skipped += 1          # leave `hits` empty -> emitted as unmatched
        else:
            # 2) PRIMARY PASS: view-aligned matching using your nutrient columns
            for view, col in VIEW_TO_COL.items():
                base_q = (r.get(col) or "").strip()
                if not base_q:
                    continue

                # Scored against base_q, NOT against the variant that produced the key.
                # That is deliberate: expand_query_variants() broadens as well as
                # rewrites, so a variant judged on its own terms lets the broadening win -
                # "alpha-tocopherol" expands to "vitamin e" and would take an exact-label
                # 100, burying the specific alpha-tocopherol at 75 and collapsing all
                # eight tocopherol/tocotrienol isomers onto one class. Scoring the alias
                # against the ORIGINAL name is what keeps the specific reading. The
                # fallback pass below scores against the variant precisely because it is a
                # last resort; see the confidence bar there.
                variants = expand_query_variants(base_q)
                keys: Set[str] = set()
                for qv in variants:
                    k = _keys_for_view(view, qv)
                    if k:
                        keys.add(k)

                for key in keys:
                    for cid, _lbl, term in unique_preserve(idx.get((view, key), [])):
                        lbl = label_by_id.get(cid, _lbl) or _lbl
                        sc = score_hit(base_q, view, lbl, term)  # score against base query for this column
                        if sc >= args.min_score:
                            hits.append((view, cid, lbl, term, sc))

        # 3) FALLBACK PASS: search expanded variants of out_query/name across ALL views.
        # Skipped for uninformative codes -- the fallback is broader than the primary pass, so
        # running it there would reinstate exactly the collisions the skip exists to prevent.
        #
        # Gated on a CONFIDENCE BAR, not on `hits` being empty. Emptiness made the rescue
        # depend on the primary pass failing completely, so any hit one point over
        # --min-score cancelled it. ChEBI release 253 supplied exactly that: it added the
        # synonym "Vitamin b1 (as thiamine nitrate)" to thiamine mononitrate, whose
        # no-paren form scores 31 against "thiamine". Under the old gate that 31 counted as
        # success, the fallback never ran, and Thiamine/Thiamin dropped from a curated
        # score-100 vitamin B1 to a score-31 fortification salt -- with no error and no
        # change in row count. The bar is set just above --min-score so only barely-over
        # hits reopen the rescue; anything the primary pass resolved with real confidence
        # (a clean exact-label match scores 66-100) still suppresses it, which is what
        # keeps "alpha-tocopherol" from being rescued into the broader "vitamin E".
        primary_best = max((h[4] for h in hits), default=0)
        if primary_best < args.fallback_bar and not bad_code:
            fallback_seed = out_query or nname
            variants = expand_query_variants(fallback_seed)
            for qv in variants:
                for view in VIEWS:
                    key = _keys_for_view(view, qv)
                    if not key:
                        continue
                    for cid, _lbl, term in unique_preserve(idx.get((view, key), [])):
                        lbl = label_by_id.get(cid, _lbl) or _lbl
                        sc = score_hit(qv, view, lbl, term)  # score against the exact variant used
                        if sc >= args.min_score:
                            hits.append((view, cid, lbl, term, sc))

        # no hits at all
        if not hits:
            all_rows.append({
                "nutrient_id": nid,
                "nutrient_name": nname,
                "query": out_query,
                "match_view": "",
                "match_source": "",
                "chebi_id": "",
                "chebi_label": "",
                "chebi_term": "",
                "score": "",
            })
            continue

        # per-CHEBI dedup: keep best row per chebi_id
        best_per_chebi: Dict[str, Tuple[str, str, str, str, int]] = {}
        for view, cid, lbl, term, sc in hits:
            prev = best_per_chebi.get(cid)
            if prev is None:
                best_per_chebi[cid] = (view, cid, lbl, term, sc)
                continue
            if (
                sc > prev[4] or
                (sc == prev[4] and VIEW_WEIGHT.get(view, 0) > VIEW_WEIGHT.get(prev[0], 0)) or
                (sc == prev[4] and VIEW_WEIGHT.get(view, 0) == VIEW_WEIGHT.get(prev[0], 0)
                 and len(lbl or "") < len(prev[2] or ""))
            ):
                best_per_chebi[cid] = (view, cid, lbl, term, sc)

        hits2 = list(best_per_chebi.values())
        hits2.sort(key=lambda x: (-x[4], -VIEW_WEIGHT.get(x[0], 0), len(x[2] or "")))

        best = hits2[0][4]
        second = hits2[1][4] if len(hits2) > 1 else None
        if second is not None and (best - second) >= args.clear_gap:
            kept = [h for h in hits2 if h[4] == best]
        else:
            kept = [h for h in hits2 if h[4] >= (best - args.keep_within)]

        # CURATION-STATUS RESCUE. ChEBI marks manually curated entries 3:STAR and
        # unvetted third-party submissions 2:STAR. The two are not interchangeable for an
        # ontology walk: a 2:STAR record typically carries a name, a formula and has_role
        # edges, but none of the structural relationships the nutrient->EC walk traverses.
        # Release 253 imported ~13,600 such records, many bearing the plain common name of
        # a compound whose curated entry is named more precisely - so the submission
        # outscores the curated term on exact-label alone and evicts it here.
        #
        # Lactic acid is the case that exposed it. CHEBI:749277 "lactic acid" (2:STAR) has
        # ONE is_a edge, to racemate, and no structural relationships at all; CHEBI:28358
        # "rac-lactic acid" (3:STAR) carries has_part (R)- and (S)-lactic acid plus
        # is_conjugate_acid_of lactate, and those edges are the entire reason lactic acid
        # reaches any enzyme. Scored 75 to 41, the dead end won and the nutrient silently
        # dropped to zero ECs - which would have broken the Veillonella-to-dairy case.
        #
        # So: never let an all-2:STAR result evict the curated alternative. This ADDS the
        # best 3:STAR candidate back rather than demoting the 2:STAR one, because the
        # submission is often the better NAME match while the curated term is the better
        # CHEMISTRY - the walk wants both, and a nutrient carrying several seeds is
        # already normal here. Nothing is removed, so no existing mapping can be lost.
        #
        # Restore them as a GROUP, under the same keep_within rule that picks among equals
        # everywhere else. Rescuing only the single best curated term silently halves any
        # nutrient whose chemistry is split across several curated entries: "lactic acid"
        # has three (rac-lactic acid, (R)-lactic acid, 2-hydroxypropanoic acid), all 3:STAR
        # and all tied at 41, and they carry 177 EC between them - keeping one leaves 66.
        if kept and not any(is_curated(h[1]) for h in kept):
            # hits2 is already sorted by descending score, so the first is the best curated.
            cur_hits = [h for h in hits2 if is_curated(h[1]) and h not in kept]
            if cur_hits:
                best_cur = cur_hits[0][4]
                kept.extend(h for h in cur_hits if h[4] >= (best_cur - args.keep_within))
                n_star_rescued += 1
        hits2 = kept

        for view, cid, lbl, term, sc in hits2:
            all_rows.append({
                "nutrient_id": nid,
                "nutrient_name": nname,
                "query": out_query,
                "match_view": view,
                "match_source": parse_source_from_term(term),
                "chebi_id": cid,
                "chebi_label": lbl,
                "chebi_term": term,
                "score": str(sc),
            })

    print(f"[*] {n_element_locked} nutrient columns resolved as elements and locked; "
          f"{n_code_skipped} left unmatched as uninformative <code>(unit) names", flush=True)
    print(f"[*] {n_star_rescued} nutrients kept a curated 3:STAR term that an unvetted "
          f"2:STAR submission would otherwise have evicted", flush=True)
    write_tsv(args.out, fieldnames, all_rows)


if __name__ == "__main__":
    main()
