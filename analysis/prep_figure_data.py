#!/usr/bin/env python3
"""prep_figure_data.py — build the small CSVs the two manuscript figures plot.

Everything here is derived from the deposited exports and the predictor outputs; no
figure number is typed by hand. Cohort sample identifiers are ANONYMIZED at this
boundary: the source cohort name must not appear in any published artifact, so the
mapping to "Infant N" is applied here and the raw identifiers never leave this script.

Figure 1 keeps two panels (source contributions, ChEBI resolution); the conflict counts
and the linkage funnel are reported in the text instead.
Figure 2: (a) cohort-level food ranking, (b) one pooled taxon x food NMDS, (c) differential view.

Outputs land beside this file as fig1*.csv / fig2*.csv.
"""
from __future__ import annotations

import csv
import os
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Export directory and output directory are both overridable so that the NEVO-free
# deposit and the NEVO-inclusive local build are produced by THIS script rather than by
# a forked copy of it: two lineages that must stay in step cannot be allowed to drift.
EXPORTS = Path(os.environ.get("BAC2FOOD_EXPORTS", "/data/bac2food/exports"))
OUTDIR = Path(os.environ.get("BAC2FOOD_FIGDIR", str(HERE)))
# Predictor output directory. Overridable so the 4-sample figures in earlier drafts stay
# reproducible after the full 55-sample run moved outputs elsewhere.
PREDICT = Path(os.environ.get("BAC2FOOD_PREDICT_DIR",
                              "/home/svalenzuela/Desktop/bac2food/4_predict"))


def write(name: str, header: list[str], rows: list) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    p = OUTDIR / name
    with p.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  {name:26} {len(rows):5} rows")


# --------------------------------------------------------------------------- fig 1
# Row counts are COUNTED from the shipped export rather than listed here. They used to be
# a hand-kept literal, which silently went stale the moment branded foods left the export
# (FDC alone moved 27,094,028 -> 1,156,398). Only the display names are hard-coded.
SOURCE_LABEL = {
    "fdc": "FoodData Central (USA)", "mccance": "McCance & Widdowson (UK)",
    "japan": "STFCJ (Japan)", "afcd": "AFCD (Australia)",
    "cnf": "Canadian Nutrient File", "biofoodcomp": "BioFoodComp (FAO)",
    "fineli": "Fineli (Finland)", "ciqual": "CIQUAL (France)",
    "nevo": "NEVO (Netherlands)", "swedish": "Livsmedelsdatabasen (Sweden)",
    "swiss": "Swiss FCDB", "wafct": "WAFCT (West Africa)",
    "phyfoodcomp": "PhyFoodComp (FAO)", "frida": "Frida (Denmark)",
    "phenol_explorer": "Phenol-Explorer (France)",
    "synthetic_bacterial": "Synthetic substrates",
    "mccance;swedish": "McCance;Swedish (agreeing)",
}
TIER = {
    "NEVO (Netherlands)": "Restricted (not redistributed)",
    "AFCD (Australia)": "Copyleft / NonCommercial",
    "BioFoodComp (FAO)": "Copyleft / NonCommercial",
    "PhyFoodComp (FAO)": "Copyleft / NonCommercial",
    "WAFCT (West Africa)": "Copyleft / NonCommercial",
}


def fig1a() -> None:
    c = Counter()
    with (EXPORTS / "food_nutrients.tsv").open(encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            c[line.rstrip("\n").rsplit("\t", 1)[-1]] += 1
    unknown = set(c) - set(SOURCE_LABEL)
    if unknown:
        raise SystemExit(f"fig1a: no display name for source_db {sorted(unknown)}")
    print(f"  [check] source rows sum to {sum(c.values()):,} across {len(c)} labels")
    write("fig1a_sources.csv", ["source", "rows", "tier"],
          [(SOURCE_LABEL[k], n, TIER.get(SOURCE_LABEL[k], "Open (redistributable)"))
           for k, n in c.most_common()])


def fig1b_chebi() -> None:
    """ChEBI match-confidence tiers, counted from the shipped digest."""
    c = Counter()
    with (EXPORTS / "enzyme_substrate_chebi.tsv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c[r["chebi_match_type"] or "no_match"] += 1

    def bucket(k: str) -> str:
        if k.startswith("rhea"):         return "Rhea-confirmed"
        if k.startswith("string_match"): return "Exact name match"
        if k.startswith("fuzzy_match"):  return "Fuzzy name match"
        if k.endswith("_no_chebi"):      return "No ChEBI entity exists*"
        if k == "ignored_ubiquitous":    return "Ubiquitous, ignored"
        return "Unmatched"

    g = Counter()
    for k, v in c.items():
        g[bucket(k)] += v
    order = ["Rhea-confirmed", "Exact name match", "Fuzzy name match",
             "No ChEBI entity exists*", "Ubiquitous, ignored", "Unmatched"]
    write("fig1b_chebi.csv", ["tier", "rows"], [(t, g.get(t, 0)) for t in order])


# --------------------------------------------------------------------------- taxonomy
# No NCBI taxonomy dump is available locally, and only 26 genera occur in the cohort, so
# family assignment is a curated genus->family map. Clostridium is polyphyletic and several
# of its former members have been formally reassigned, so those are handled per species
# BEFORE the genus lookup rather than lumped under Clostridiaceae.
GENUS_FAMILY = {
    "Anaerostipes": "Lachnospiraceae", "Atopobium": "Atopobiaceae",
    "Bifidobacterium": "Bifidobacteriaceae", "Citrobacter": "Enterobacteriaceae",
    "Clostridium": "Clostridiaceae", "Enterocloster": "Lachnospiraceae",
    "Enterococcus": "Enterococcaceae", "Erysipelatoclostridium": "Erysipelotrichaceae",
    "Eubacterium": "Eubacteriaceae", "Faecalimonas": "Lachnospiraceae",
    "Flavonifractor": "Oscillospiraceae", "Haemophilus": "Pasteurellaceae",
    "Intestinibacter": "Peptostreptococcaceae", "Klebsiella": "Enterobacteriaceae",
    "Lacticaseibacillus": "Lactobacillaceae", "Ligilactobacillus": "Lactobacillaceae",
    "Limosilactobacillus": "Lactobacillaceae", "Megasphaera": "Veillonellaceae",
    "Micrococcaceae": "Micrococcaceae", "Roseburia": "Lachnospiraceae",
    "Ruminococcus": "Oscillospiraceae", "Streptococcus": "Streptococcaceae",
    "Subdoligranulum": "Oscillospiraceae", "Tyzzerella": "Lachnospiraceae",
    "Veillonella": "Veillonellaceae", "uncultured": "Unclassified",
}
SPECIES_FAMILY = {
    "Clostridium_innocuum": "Erysipelotrichaceae",     # reassigned; Erysipelatoclostridium
    "Clostridium_clostridioforme": "Lachnospiraceae",  # reassigned; Enterocloster
}


def family_of(bacterium: str) -> tuple[str, str]:
    """(display name, family) from a '<taxid>_Genus_species' predictor label."""
    name = re.sub(r"^\d+_", "", bacterium)
    parts = name.split("_")
    genus = parts[0]
    for k, v in SPECIES_FAMILY.items():
        if name.startswith(k):
            return name.replace("_", " "), v
    return name.replace("_", " "), GENUS_FAMILY.get(genus, "Other")


# --------------------------------------------------------------------------- fig 2
def anonymize() -> dict[str, str]:
    """Map raw predictor output prefixes onto neutral infant labels.

    Sorted by (timepoint, name) so the assignment is deterministic across reruns and
    carries no information about the original identifiers.
    """
    keys = []
    for f in sorted(PREDICT.glob("*.community.tsv")):
        stem = f.name.replace(".community.tsv", "")
        m = re.search(r"_(\d+)_month", stem)
        keys.append((int(m.group(1)) if m else 0, stem))
    keys.sort()
    return {stem: f"Infant {i} ({mo} mo)" for i, (mo, stem) in enumerate(keys, 1)}


def fig2a(label: dict[str, str]) -> None:
    """Cohort-level view: foods ranked across all sampled infants, not one infant."""
    per_food: dict[str, list[float]] = defaultdict(list)
    cat: dict[str, str] = {}
    for stem in label:
        with (PREDICT / f"{stem}.community.tsv").open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                per_food[r["food_name"]].append(float(r["score"]))
                cat.setdefault(r["food_name"], r["food_category"] or "Uncategorized")
    n_cohort = len(label)
    rows = [(f, round(sum(v) / len(v), 4), len(v), n_cohort,
             round(min(v), 4), round(max(v), 4), cat[f])
            for f, v in per_food.items()]
    rows.sort(key=lambda x: -x[1])                   # by score; n_infants is shown, not sorted on
    rows = rows[:12]
    write("fig2a_cohort_foods.csv",
          ["food", "mean_score", "n_infants", "n_cohort", "min_score", "max_score",
           "category"], rows)
    # The community view turns out to be near-invariant between infants; quantify it rather
    # than let a range bar imply a spread that is not there. The MEDIAN is the figure the
    # manuscript quotes, so print it: this line used to report only the max, which meant a
    # rebuild could never surface drift in the number actually published.
    spread = sorted((x[5] - x[4]) / x[1] for x in rows if x[1])
    med = statistics.median(spread)
    print(f"  [check] within-food score spread across infants: median {100*med:.2f} % "
          f"(quoted in text), max {100*spread[-1]:.2f} %")


def fig2b_nmds(label: dict[str, str]) -> None:
    """One long-format taxon x food profile, pooled over the cohort, for a single NMDS.

    Each taxon contributes ONE row per food: the mean of the scores it received in the
    infants where it was detected. Pooling is not cosmetic. The same species scores
    differently between infants (max s.d. 3.6) because per-sample gene recovery yields a
    slightly different EC set, and stratifying leaves each taxon with only the ~20 foods
    of one infant's truncated list -- 95 % sparsity, which drives metaMDS into a
    degenerate near-zero-stress solution. The union over infants gives every taxon a
    fuller profile and a genuine ordination (stress ~0.23, converged).

    NOTE this is the predictor's per-bacterium TOP-ranked food list, i.e. a truncated
    profile, not an exhaustive bacterium x food matrix. The legend says so, and the
    ordination is therefore run on presence/absence (Sorensen), because rank depth is
    comparable between taxa while raw score magnitude within a truncated list is not.
    """
    acc: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    fam_of: dict[str, str] = {}
    for stem, lab in label.items():
        p = PREDICT / f"{stem}.differential.tsv"
        if not p.exists():
            continue
        seen = set()
        with p.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                taxon, fam = family_of(r["bacterium"])
                k = (lab, taxon, r["food_name"])
                if k in seen:            # one score per (infant, taxon, food)
                    continue
                seen.add(k)
                fam_of[taxon] = fam
                acc[(taxon, r["food_name"], r["food_category"] or "Uncategorized")].append(
                    float(r["score"]))

    rows = [(t, fam_of[t], food, cat, round(sum(v) / len(v), 4), len(v))
            for (t, food, cat), v in sorted(acc.items())]
    write("fig2b_nmds_long.csv",
          ["taxon", "family", "food", "category", "mean_score", "n_infants"], rows)

    taxa = {r[0] for r in rows}
    foods = {r[2] for r in rows}
    print(f"  [check] pooled profile: {len(taxa)} taxa x {len(foods)} foods, "
          f"sparsity {1 - len(rows) / (len(taxa) * len(foods)):.3f}")
    fams = Counter(fam_of[t] for t in taxa)
    print("  [families] " + ", ".join(f"{k}:{v}" for k, v in fams.most_common()))


def fig2c(label: dict[str, str]) -> None:
    want = ["Bifidobacterium breve", "Bifidobacterium longum", "Bifidobacterium adolescentis",
            "Veillonella parvula", "Streptococcus thermophilus", "Limosilactobacillus fermentum"]
    # Foods are taken in the order the PREDICTOR ranks them (the `rank` column), not
    # re-sorted here. The differential table now admits a food where score > peer_median and
    # then ranks the survivors on absolute score; re-sorting on comp_score would make the
    # figure disagree with the table a reader reproduces from the deposited files.
    # Both quantities are carried through so the panel can be drawn on either.
    best: dict[str, list] = defaultdict(list)
    for stem in label:
        p = PREDICT / f"{stem}.differential.tsv"
        if not p.exists():
            continue
        with p.open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                taxon, _ = family_of(r["bacterium"])
                if taxon in want:
                    best[taxon].append((int(r["rank"]), r["food_name"],
                                        float(r["score"]), float(r["comp_score"])))
    rows = []
    for nm, items in best.items():
        seen: dict[str, tuple[float, float]] = {}
        for _rk, food, sc, comp in sorted(items, key=lambda x: x[0]):
            seen.setdefault(food, (sc, comp))
            if len(seen) == 5:
                break
        g, sp = nm.split()[0], " ".join(nm.split()[1:])
        rows.extend((f"{g[0]}. {sp}", f, v[0], v[1]) for f, v in seen.items())
    write("fig2c_differential.csv", ["taxon", "food", "score", "comp_score"], rows)
    missing = [t for t in want if t not in best]
    if missing:
        print(f"  [check] fig2c: {len(best)} of {len(want)} taxa have admitted foods; "
              f"absent: {missing}")


if __name__ == "__main__":
    print("figure data:")
    fig1a(); fig1b_chebi()
    lab = anonymize()
    fig2a(lab); fig2b_nmds(lab); fig2c(lab)
    print(f"  [anon] {len(lab)} samples -> {sorted(lab.values())}")
    print("done.")
