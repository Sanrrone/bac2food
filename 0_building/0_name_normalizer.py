#!/usr/bin/env python3
"""
name_normalizer.py

Deterministic, regex-driven normalizer for nutrient / enzyme-substrate names.
Goal: maximize exact/synonym matching against ontologies (e.g., ChEBI OBO),
while avoiding meaning-changing transformations.

It outputs multiple normalized "views" so you can try matches in a safe order:
  1) norm            (conservative)
  2) norm_no_commas  (helps comma-phrases)
  3) norm_no_paren   (removes ALL parentheses, if needed)
  4) norm_compact    (last resort for fuzzy matching; removes most punctuation)

USAGE
-----
1) Import as a module:

    from name_normalizer import normalize_name
    n = normalize_name("Fiber, crude (DO NOT USE - Archived)")
    print(n.norm)

2) CLI: normalize a column from TSV/CSV and write a new TSV.

    # Example for your FoodDataCentral nutrient.csv (quoted CSV)
    python 0_name_normalizer.py --input nutrient.csv --format csv --column name --output nutrient.normalized.tsv


The output TSV always contains:
  raw, cleaned, norm, norm_no_commas, norm_no_paren, norm_compact
plus any original columns from the input file.
"""

from __future__ import annotations
import argparse
import csv
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

# ---------------------------------------------------------------------------
# Regex rules (conservative)
# ---------------------------------------------------------------------------

# Admin / curation notes that should never be used for matching.
_ADMIN_PATTERNS = [
    r"\bdo\s*not\s*use\b",
    r"\barchiv(?:ed|e)\b",
    r"\bdeprecated\b",
    r"\bobsolete\b",
    r"\bnot\s+for\s+use\b",
    r"\bplaceholder\b",
    r"\bunknown\b",
]

# Nutrition metadata phrases that usually aren't part of the chemical name.
_META_PATTERNS = [
    r"\batwater\b",
    r"\bgeneral\s+factors?\b",
    r"\bspecific\s+factors?\b",
    r"\bby\s+difference\b",
    r"\bas\s+(?:n|as)\b",  # "as N", "as Na", etc. (handled more carefully below)
]

# Common trailing units/qualifiers that block matches (keep conservative).
_TRAILING_DROP = [
    r"\bequivalents?\b",
    r"\beq\.?\b",
    r"\btotal\b",
]

# Leading quantifier words frequently present in nutrients.
_LEADING_DROP = [
    r"^\s*total\s+",
    r"^\s*sum\s+of\s+",
]

# Normalize Greek letters and some symbols.
_GREEK_MAP = str.maketrans({
    "α": "alpha",
    "β": "beta",
    "γ": "gamma",
    "δ": "delta",
    "κ": "kappa",
    "μ": "mu",
})

# Known shorthand expansions (high precision).
# Add to this gradually as you find recurring enzyme assay substrates.
_SHORTHANDS = {
    "pnp": "p-nitrophenyl",
    "4np": "4-nitrophenyl",
    "onp": "o-nitrophenyl",
    "xos": "xylooligosaccharide",
    "gos": "galactooligosaccharide",
    "fos": "fructooligosaccharide",
}

# Preserve chemically meaningful tokens like D-/L-, alpha/beta, numbers, linkage arrows.
# We'll normalize arrows and dashes but avoid stripping these.
_ARROW_MAP = {
    "→": "->",
    "⟶": "->",
    "⟹": "->",
    "⇒": "->",
}

# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class NormalizedName:
    raw: str
    cleaned: str
    norm: str
    norm_no_commas: str
    norm_no_paren: str
    norm_compact: str


def _collapse_ws(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def _remove_admin_phrases_in_parens(inner: str) -> bool:
    low = inner.lower()
    return any(re.search(p, low) for p in _ADMIN_PATTERNS)


def _paren_rewriter(match: re.Match) -> str:
    """
    Decide what to do with parentheticals.
    - Remove parentheses that are clearly administrative (DO NOT USE, Archived, etc.)
    - Otherwise keep, but space-pad them so downstream normalization can handle.
    """
    inner = match.group(1).strip()
    if not inner:
        return " "
    if _remove_admin_phrases_in_parens(inner):
        return " "
    # If very short plain word like "(fat)" or "(total)" we drop it (often redundant)
    if len(inner) <= 12 and re.fullmatch(r"[A-Za-z \-]+", inner):
        # keep if it's stereochem (rare in this pattern) - simple heuristic:
        if re.search(r"\b(D|L|DL)\b", inner):
            return f" ({inner}) "
        return " "
    return f" ({inner}) "


def _strip_parentheticals_selective(s: str) -> str:
    # Replace "(...)" using the policy above
    return re.sub(r"\(([^)]*)\)", _paren_rewriter, s)


def _basic_cleanup(raw: str) -> str:
    s = (raw or "").strip()
    if not s:
        return ""

    # normalize unicode symbols
    s = s.translate(_GREEK_MAP)
    for k, v in _ARROW_MAP.items():
        s = s.replace(k, v)
    s = s.replace("–", "-").replace("—", "-").replace("−", "-")

    # selectively strip parentheses
    s = _strip_parentheticals_selective(s)

    # remove leading phrases
    for pat in _LEADING_DROP:
        s = re.sub(pat, "", s, flags=re.IGNORECASE)

    # drop obvious metadata phrases anywhere (conservative)
    for pat in _META_PATTERNS:
        s = re.sub(pat, " ", s, flags=re.IGNORECASE)

    # drop trailing qualifiers (as whole words)
    for pat in _TRAILING_DROP:
        s = re.sub(rf"\s+{pat}\s*$", " ", s, flags=re.IGNORECASE)

    # expand some shorthands as separate tokens (avoid changing chemical tokens inside longer words)
    def expand_token(m: re.Match) -> str:
        tok = m.group(0)
        return _SHORTHANDS.get(tok.lower(), tok)

    s = re.sub(r"\b(pnp|4np|onp|xos|gos|fos)\b", expand_token, s, flags=re.IGNORECASE)

    # normalize separators but keep chemistry tokens
    s = s.replace("/", " / ")
    s = s.replace(",", ", ")
    s = s.replace(";", " ; ")
    s = _collapse_ws(s)
    
    # remove empty parentheses created by stripping metadata
    s = re.sub(r"\(\s*\)", " ", s)
    s = re.sub(r"\b([A-Za-z]+),\s*[A-Z][a-z]?\b$", r"\1", s)

    s = _collapse_ws(s)

    return s


def _normalize_core(s: str) -> str:
    s = s.lower().strip()

    # normalize spacing around punctuation while preserving chemistry tokens
    s = re.sub(r"\s*-\s*", "-", s)       # keep hyphenated tokens
    s = re.sub(r"\s*/\s*", "/", s)
    s = re.sub(r"\s*;\s*", "; ", s)
    s = re.sub(r"\s*,\s*", ", ", s)
    s = _collapse_ws(s)

    # remove duplicated commas/semicolons spacing issues
    s = re.sub(r",\s*,+", ", ", s)
    s = re.sub(r";\s*;+", "; ", s)
    s = re.sub(r"^(.+?)\s*\1$", r"\1", s)
    s = re.sub(r",\s*$", "", s)

    return s


def _remove_all_parentheticals(s: str) -> str:
    return _collapse_ws(re.sub(r"\([^)]*\)", " ", s))


def _compact_last_resort(s: str) -> str:
    # Remove almost everything except alphanumerics; intended for fuzzy candidates only.
    s = s.lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return _collapse_ws(s)


def normalize_name(raw: str) -> NormalizedName:
    cleaned = _basic_cleanup(raw)
    norm = _normalize_core(cleaned)
    norm_no_commas = _collapse_ws(norm.replace(",", " "))
    norm_no_paren = _normalize_core(_remove_all_parentheticals(cleaned))
    norm_compact = _compact_last_resort(cleaned)

    return NormalizedName(
        raw=str(raw or ""),
        cleaned=cleaned,
        norm=norm,
        norm_no_commas=norm_no_commas,
        norm_no_paren=norm_no_paren,
        norm_compact=norm_compact,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _read_rows(path: str, fmt: str) -> Iterable[Dict[str, str]]:
    if fmt == "csv":
        with open(path, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f)
            for row in r:
                yield {k: (v if v is not None else "") for k, v in row.items()}
    elif fmt == "tsv":
        with open(path, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                yield {k: (v if v is not None else "") for k, v in row.items()}
    else:
        raise ValueError("format must be csv or tsv")


def _write_tsv(path: str, fieldnames: List[str], rows: Iterable[Dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Input file (CSV or TSV).")
    ap.add_argument("--format", required=True, choices=["csv", "tsv"], help="Input format.")
    ap.add_argument("--column", required=True, help="Column name to normalize.")
    ap.add_argument("--output", required=True, help="Output TSV path.")
    args = ap.parse_args()

    out_rows = []
    n_total = 0
    n_empty = 0

    for row in _read_rows(args.input, args.format):
        n_total += 1
        raw = row.get(args.column, "")
        nn = normalize_name(raw)

        if not nn.cleaned:
            n_empty += 1

        row_out = dict(row)
        row_out["normalized_raw"] = nn.raw
        row_out["normalized_cleaned"] = nn.cleaned
        row_out["normalized_norm"] = nn.norm
        row_out["normalized_norm_no_commas"] = nn.norm_no_commas
        row_out["normalized_norm_no_paren"] = nn.norm_no_paren
        row_out["normalized_norm_compact"] = nn.norm_compact
        out_rows.append(row_out)

    # output columns: original + added (preserve original order if possible)
    if out_rows:
        orig_fields = list(out_rows[0].keys())
        # ensure normalized fields are last
        norm_fields = [
            "normalized_raw",
            "normalized_cleaned",
            "normalized_norm",
            "normalized_norm_no_commas",
            "normalized_norm_no_paren",
            "normalized_norm_compact",
        ]
        # dedupe while preserving order
        seen = set()
        fieldnames = []
        for f in orig_fields:
            if f not in seen:
                seen.add(f)
                fieldnames.append(f)
        # make sure normalized fields exist and are at end
        fieldnames = [f for f in fieldnames if f not in norm_fields] + norm_fields

        _write_tsv(args.output, fieldnames, out_rows)

    print(f"Wrote {args.output} (rows={n_total}, empty_cleaned={n_empty})")


if __name__ == "__main__":
    main()

