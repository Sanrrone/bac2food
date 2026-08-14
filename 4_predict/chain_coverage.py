#!/usr/bin/env python3
"""chain_coverage.py — where does bacterium -> EC -> substrate -> ChEBI -> nutrient break?

The feature match rate (see recompute_match_rate.py) says how much of a cohort reaches a
food nutrient. This says WHY the rest does not, which is what decides whether a newer
source release would help:

  GAP A  EC has no substrate at all           -> BRENDA coverage gap; a refresh COULD help
  GAP B  substrate present, no ChEBI id       -> name-matching gap in ../chebi/
  GAP C  ChEBI id present, no FDC nutrient    -> the nutrient vocabulary stops here;
                                                 no substrate-source refresh can close it

Measured 2026-08: GAP C dominates (52% of cohort EC, 49% of resource EC), and its members
act on real ChEBI-identified molecules that simply are not food components (tRNA ligases,
cell-wall enzymes, phosphorylated intermediates). The 45.4% feature match rate is therefore
close to a structural ceiling set by what counts as a nutrient, not a fixable coverage
deficiency -- BRENDA already supplies substrates for ~95% of the cohort's EC numbers.

Usage:
    python chain_coverage.py
    python chain_coverage.py --top 20
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import importlib.util
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent


def _load_split_ec_cell():
    """Reuse recompute_match_rate.split_ec_cell so both tools split multi-EC cells alike."""
    spec = importlib.util.spec_from_file_location("_rmr", HERE / "recompute_match_rate.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_rmr"] = mod
    spec.loader.exec_module(mod)
    return mod.split_ec_cell


split_ec_cell = _load_split_ec_cell()

# Substrates that are not food: nucleic acid, protein/tRNA, cell envelope, or pure
# energy/redox cofactors. Used only to characterise GAP C, never to filter data.
NONDIET_RE = re.compile(
    r"(?i)\b(dna|rna|trna|rrna|mrna|oligonucleotide|"
    r"peptidoglycan|murein|lipopolysaccharide|lipid a|teichoic|"
    r"protein|polypeptide|peptidyl|ribosom|thioredoxin|glutaredoxin|ferredoxin|"
    r"ubiquinone|menaquinone|quinone|cytochrome|heme|haem|"
    r"atp|adp|amp|gtp|gdp|nad|nadp|fad|fmn|coenzyme a|acyl-carrier|s-adenosyl)\b")


def load_chain(digest_chebi: Path, nutrient_to_ec: Path):
    ec_sub: dict[str, set] = collections.defaultdict(set)
    ec_sub_chebi: dict[str, set] = collections.defaultdict(set)
    with digest_chebi.open(encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ec, su = r["ec_number"], r["substrate"]
            ec_sub[ec].add(su)
            if r.get("chebi_id"):
                ec_sub_chebi[ec].add(su)
    ec_nut: set[str] = set()
    with nutrient_to_ec.open(encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("ec_number"):
                ec_nut.add(r["ec_number"])
    return ec_sub, ec_sub_chebi, ec_nut


def load_cohort_ecs(annot_dir: Path) -> collections.Counter:
    """EC -> annotated-locus count, splitting multi-EC cells as the predictor does."""
    c: collections.Counter = collections.Counter()
    for fp in sorted(glob.glob(str(annot_dir / "*_ec.tsv"))):
        with open(fp, encoding="utf-8", errors="replace") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) == 4 and p[3]:
                    for ec in split_ec_cell(p[3]):
                        c[ec] += 1
    return c


def load_resource_ecs(path: Path) -> set:
    ecs = set()
    with path.open(encoding="utf-8") as f:
        next(f)
        for line in f:
            ecs.add(line.rsplit("\t", 1)[-1].strip())
    return ecs


def report(name, ecs, ec_sub, ec_sub_chebi, ec_nut, weights=None):
    tot = len(ecs)
    wsum = (lambda s: sum(weights[e] for e in s)) if weights else None
    wt = wsum(ecs) if weights else 0
    has_sub = {e for e in ecs if ec_sub.get(e)}
    has_chebi = {e for e in has_sub if ec_sub_chebi.get(e)}
    reach = {e for e in ecs if e in ec_nut}

    print(f"\n=== {name} ===")

    def line(tag, s):
        extra = f"   |  {wsum(s):>9,} loci ({100*wsum(s)/wt:5.1f}%)" if weights else ""
        print(f"  {tag:<46}{len(s):>5,} ({100*len(s)/tot:5.1f}%){extra}")

    line("total EC", ecs)
    line("1. has >=1 substrate", has_sub)
    line("2. ...>=1 substrate with a ChEBI id", has_chebi)
    line("3. ...reaches an FDC nutrient", reach)
    gapA, gapB, gapC = ecs - has_sub, has_sub - has_chebi, has_chebi - reach
    print("  " + "-" * 68)
    line("GAP A no substrate      (refresh COULD help)", gapA)
    line("GAP B no ChEBI id       (name-matching gap)", gapB)
    line("GAP C no nutrient       (refresh WON'T help)", gapC)
    return gapA, gapB, gapC


def main() -> int:
    ap = argparse.ArgumentParser(description="Diagnose where the linkage chain breaks.")
    ap.add_argument("--annot_dir", default=str(REPO / "gene_annot"))
    ap.add_argument("--digest_chebi", default=str(REPO / "chebi/digest_to_chebi.tsv"))
    ap.add_argument("--nutrient_to_ec", default=str(REPO / "0_building/3_nutrient_to_ec.tsv"))
    ap.add_argument("--species_enzymes", default="/data/bac2food/exports/species_enzymes.tsv")
    ap.add_argument("--top", type=int, default=12, help="How many gap EC to list")
    args = ap.parse_args()

    ec_sub, ec_sub_chebi, ec_nut = load_chain(Path(args.digest_chebi), Path(args.nutrient_to_ec))
    print(f"chain: {len(ec_sub):,} EC with substrates | {len(ec_sub_chebi):,} with a ChEBI id "
          f"| {len(ec_nut):,} reaching a nutrient")

    cohort = load_cohort_ecs(Path(args.annot_dir))
    gapA, _, gapC = report("COHORT EC (what a cohort run sees)", set(cohort),
                           ec_sub, ec_sub_chebi, ec_nut, weights=cohort)
    report("RESOURCE EC (exported species_enzymes)", load_resource_ecs(Path(args.species_enzymes)),
           ec_sub, ec_sub_chebi, ec_nut)

    # Is GAP C a data deficiency, or enzymes with no dietary substrate to find?
    nd = sum(1 for e in gapC if all(NONDIET_RE.search(s) for s in ec_sub[e]))
    print(f"\nGAP C character: {nd:,} of {len(gapC):,} act only on nucleic acid / protein / "
          f"cofactor substrates;")
    print(f"  the other {len(gapC)-nd:,} act on real small molecules that are simply not FDC "
          f"nutrients (intracellular\n  intermediates, cell-wall precursors) -- a nutrient-vocabulary "
          f"boundary, not a source gap.")

    print(f"\ntop GAP A EC by cohort abundance (no substrate recorded):")
    for ec, n in sorted(((e, cohort[e]) for e in gapA), key=lambda x: -x[1])[:args.top]:
        print(f"  {ec:<12} {n:>8,} loci")
    print(f"\ntop GAP C EC by cohort abundance (ChEBI but no nutrient):")
    for ec, n in sorted(((e, cohort[e]) for e in gapC), key=lambda x: -x[1])[:args.top]:
        print(f"  {ec:<12} {n:>8,} loci")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
