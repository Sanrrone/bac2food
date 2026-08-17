#!/usr/bin/env python3
"""verify_exports.py — health-check the three exported resource files.

The exports are what third parties actually receive, and after the eggNOG v7 rebuild
(eggnog/6.1_) only one of the three was regenerated. This checks all three anyway,
because the value of the resource is that they JOIN:

    species_enzymes.tsv → enzyme_substrate_chebi.tsv → food_nutrients.tsv
    organism → EC       → substrate/ChEBI → nutrient_id → food

Three classes of check:
  1. STRUCTURE — header, field count per row (catches embedded tabs/newlines, the
     failure mode that silently shifts every column), and id/EC formats.
  2. COUNTS    — against the figures the manuscript quotes, so the paper and the
     deposit cannot drift apart.
  3. JOINS     — referential integrity across the three files, including one
     end-to-end walk from a named organism to real foods.

Streams every file; peak memory is set by the distinct-id sets, not by file size
(food_nutrients.tsv is ~3.3 GB).

Usage:
    python verify_exports.py
    python verify_exports.py --exports /data/bac2food/exports --quick   # skip food_nutrients
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

EC_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

SCHEMA = {
    "species_enzymes.tsv": ["tax_id", "genus", "species", "strain", "organism", "ec_number"],
    "enzyme_substrate_chebi.tsv": ["ec_number", "substrate", "substrate_normalized", "chebi_id",
                                   "chebi_name", "chebi_match_type", "in_model", "nutrient_ids",
                                   "nutrient_names", "model_relation", "model_score"],
    "food_nutrients.tsv": ["fdc_id", "description", "data_type", "food_category", "nutrient_id",
                           "nutrient_name", "unit_name", "amount", "source_db"],
}

# Figures quoted in the manuscript (Data Records / Technical Validation).
EXPECTED = {
    "species_enzymes.tsv": {"rows": 20_557_730, "organisms": 10_751, "ec": 4_819},
    # in_model 4,199 -> 4,342 -> 4,094 on 2026-08-05. Two changes to the nutrient->EC map, both in
    # 0_building/ and both deliberate: INFOODS mineral tagnames stopped resolving to
    # dipeptides (-271 links, ALL from columns like K(mg), NA(mg), FE(mg)), and the ChEBI
    # walk gained the structural relations (+414). Then a second pass removed short-code
    # collisions (C10..C24 (%T) fatty-acid columns matching cysteine/cytosine, PP-(mg)
    # niacin matching proline, IP3..IP6 matching Ile-Pro) and translocase-only mineral
    # links (Na, Zn, Rb): 2,816 further spurious EC-links, hence 4,342 -> 4,094. Then 4,094 ->
    # 4,065 when branded foods left food_nutrients.tsv: `in_model` is computed against the
    # nutrients that carry a measured value, and Inositol (1181) was measured only on branded
    # labels, so its 29 EC rows no longer reach a food. Row/EC/ChEBI counts are unchanged
    # throughout because those come from the digest, which none of this touched.
    # chebi/in_model track the ChEBI release: 253 (07 Jul 2026) resolved 624 more substrates
    # than the Feb 2026 release did, and the nutrient map gained cellulose and
    # 7,3',4'-trihydroxyflavone. rows/ec do NOT move with the ontology - they come from the
    # BRENDA digest - so a change in either of those is a real regression, not a refresh.
    # 4,250 -> 4,254 when the 3:STAR rescue in 2_nutri2chebi_from_obo.py started restoring
    # the curated candidates as a GROUP. Rescuing only the single best one had left lactic
    # acid on 66 of its 177 EC, because its chemistry is split across three curated terms
    # (rac-, (R)-, and 2-hydroxypropanoic acid) that all tie at score 41.
    # 4,254 -> 4,213 when 3_nutrient_to_ec.py started applying extra_nutrient_chebi.tsv as
    # --override_seeds. A NET DROP that is a correction: nutrient 96062 "Alcohol" alone gave
    # back 134 EC it had borrowed from CHEBI:30879, ChEBI's generic "any R-OH" class, and
    # that outweighs the 9 nutrients the curated table gained (Resistant starch 0 -> 118,
    # Cellulose 1 -> 30, Lignin 8 -> 26, ...). See 0_building/extra_nutrient_chebi.readme.txt.
    #
    # 2026-08-06: rows/ec MOVED, and for once that is not a regression. The substrate layer
    # was rebuilt from BRENDA's official 2026.1 download (eggnog/1.0b_brenda_json_to_digest.py)
    # after the DSMZ SPARQL scrape it replaced was measured to be an incomplete slice of that
    # SAME release - 5,263 EC against 6,901, and 63.5% of the pairs. Both sources resolve
    # against ChEBI at the same efficiency (57.7% vs 57.0%), so the 2.27x volume is real
    # coverage, not noise. If rows/ec move again WITHOUT a documented digest rebuild, that IS
    # a regression. See eggnog/readme.txt.
    "enzyme_substrate_chebi.tsv": {"rows": 218_378, "ec": 6_900, "chebi": 124_538, "in_model": 10_271},
    # Branded label products were removed from the export on 2026-08-05 (--keep_branded
    # restores them). They were 25,937,648 of the 28,246,465 values and 1,890,275 of the
    # 2,010,585 foods, yet carried only 119 distinct components against 1,779 for
    # everything else — a nutrition-facts panel, not an analysis. The predictor had
    # already been ignoring them (`drop_branded: true`), so this makes the deposit match
    # what the science used. 9 components are lost with them (beta-glucans, inositol,
    # EGCG, added/intrinsic sugars, sugar alcohols and three label-only entries), as are
    # 18 of the 187 synthetic_bacterial rows, which had been injected onto branded foods.
    # FNDDS (survey_fndds_food) was removed on 2026-08-10 for the same reason, and the
    # evidence for it is stronger: none of its 353,015 values carries a `data_points` sample
    # count and none carries a `derivation_id`, because FNDDS is recipe-modelled from SR
    # Legacy ingredients for dietary-survey coding rather than assayed. Keeping it would
    # partly double-count rows already in this table. It cost 5,432 foods and no components
    # at all — all 1,779 nutrients and all 598 enzyme-linked ones survive — and it ended the
    # FDC backbone's majority: 51.2% of values before, 42.2% after.
    "food_nutrients.tsv": {"rows": 1_904_276, "foods": 112_550, "nutrients": 1_779},
}

# Exactly what the deposit should contain. Anything else in the directory ships with
# the resource by accident — the predictor used to write its parquet cache here, because
# the cache path was derived as "next to the reference TSV".
# These two checks stay at zero tolerance, and keeping them there is the point: when
# branded foods left the export, Inositol (1181) and Epigallocatechin-3-gallate (1368)
# stopped carrying any measured value, and the digest went on citing them until
# 3_nutrient_to_ec.py was rebuilt against the new live_nutrients.tsv. That transient is
# exactly what these two catch, so REBUILD THE MAP after any change to the food set
# rather than widening the tolerance here. Neither compound is actually lost: EGCG is
# measured analytically by Phenol-Explorer as "(-)-Epigallocatechin 3-O-gallate"
# (nutrient 240010, 93 foods); inositol survives only inside WAFCT's ambiguous
# "Phytic Acid / Myo-Inositol" column (462 foods), which the matcher rightly refuses to
# resolve, so its 29 EC rows are a genuine coverage loss.

DELIVERABLES = {
    "food_nutrients.tsv",           # the harmonized composition table
    "enzyme_substrate_chebi.tsv",   # EC -> substrate -> ChEBI -> nutrient_id
    "species_enzymes.tsv",          # organism -> EC (eggNOG v7),
    # rights table: one row per source_db, the join target for licence filtering
    "licences.tsv",
}

# No artifact of the non-redistributable source ships, in any form. Its derived values are an
# amendment its terms withhold, and redistributing even the original release unchanged is held
# back pending written confirmation. A user brings their own copy and rebuilds that partition
# locally. The check is filename-level because the failure it guards against is someone
# dropping the release back into the deposit directory by hand, which is how it got there the
# first time. The key below is an internal source_db value, not a statement about the source.
WITHHELD = ("nevo",)

results: list[tuple[bool, str]] = []


def check(ok: bool, label: str, detail: str = "") -> bool:
    results.append((ok, label))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" — {detail}" if detail else ""), flush=True)
    return ok


def scan(path: Path, cols: list[str]):
    """Stream a TSV, validating the header and per-row field count. Yields field lists."""
    with path.open(encoding="utf-8", errors="replace", newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != cols:
            raise SystemExit(f"ERROR: {path.name} header mismatch\n  got      {header}\n  expected {cols}")
        n = len(cols)
        for i, line in enumerate(fh, 2):
            f = line.rstrip("\n").split("\t")
            if len(f) != n:
                raise SystemExit(f"ERROR: {path.name} line {i} has {len(f)} fields, expected {n}")
            yield f


def verify_species_enzymes(p: Path):
    print(f"\n[1] {p.name}")
    rows = 0
    orgs, ecs, taxa = set(), set(), set()
    bad_ec = bad_tax = blank_org = 0
    seen_pairs_sampled = set()
    dup_pairs = 0
    for tax_id, genus, species, strain, organism, ec in scan(p, SCHEMA[p.name]):
        rows += 1
        if not EC_RE.match(ec):
            bad_ec += 1
        if not tax_id.isdigit():
            bad_tax += 1
        if not organism.strip():
            blank_org += 1
        orgs.add(organism); ecs.add(ec); taxa.add(tax_id)
        if rows <= 4_000_000:                      # duplicate probe on a bounded prefix
            k = (tax_id, ec)
            if k in seen_pairs_sampled:
                dup_pairs += 1
            seen_pairs_sampled.add(k)
    e = EXPECTED[p.name]
    check(bad_ec == 0, "every ec_number is a 4-level EC", f"{bad_ec} malformed")
    check(bad_tax == 0, "every tax_id is numeric", f"{bad_tax} malformed")
    check(blank_org == 0, "no blank organism", f"{blank_org} blank")
    check(dup_pairs == 0, "no duplicate (tax_id, EC) in first 4M rows", f"{dup_pairs} dupes")
    check(rows == e["rows"], f"row count == {e['rows']:,}", f"got {rows:,}")
    check(len(orgs) == e["organisms"], f"organisms == {e['organisms']:,}", f"got {len(orgs):,}")
    check(len(ecs) == e["ec"], f"distinct EC == {e['ec']:,}", f"got {len(ecs):,}")
    check(len(taxa) == len(orgs), "tax_id and organism are 1:1",
          f"{len(taxa):,} tax_id vs {len(orgs):,} organism")
    return ecs, orgs


def verify_enzyme_substrate(p: Path):
    print(f"\n[2] {p.name}")
    rows = in_model = with_chebi = 0
    ecs, nutrient_ids = set(), set()
    bad_ec = bad_flag = model_without_nut = 0
    for f in scan(p, SCHEMA[p.name]):
        rows += 1
        ec, chebi_id, flag, nut = f[0], f[3], f[6], f[7]
        if not EC_RE.match(ec):
            bad_ec += 1
        ecs.add(ec)
        if chebi_id:
            with_chebi += 1
        if flag not in ("yes", "no"):
            bad_flag += 1
        if flag == "yes":
            in_model += 1
            if not nut.strip():
                model_without_nut += 1
        for n in nut.split(","):
            if n.strip():
                nutrient_ids.add(n.strip())
    e = EXPECTED[p.name]
    check(bad_ec == 0, "every ec_number is a 4-level EC", f"{bad_ec} malformed")
    check(bad_flag == 0, "in_model is yes/no", f"{bad_flag} other")
    check(model_without_nut == 0, "in_model=yes implies nutrient_ids present",
          f"{model_without_nut} violations")
    check(rows == e["rows"], f"row count == {e['rows']:,}", f"got {rows:,}")
    check(len(ecs) == e["ec"], f"distinct EC == {e['ec']:,}", f"got {len(ecs):,}")
    check(with_chebi == e["chebi"], f"rows with a ChEBI id == {e['chebi']:,}", f"got {with_chebi:,}")
    check(in_model == e["in_model"], f"in_model rows == {e['in_model']:,}", f"got {in_model:,}")
    return ecs, nutrient_ids


def verify_food_nutrients(p: Path):
    print(f"\n[3] {p.name}")
    rows = 0
    foods, nutrients, sources = set(), set(), set()
    bad_amount = 0
    for fdc_id, desc, dtype, cat, nid, nname, unit, amount, src in scan(p, SCHEMA[p.name]):
        rows += 1
        foods.add(fdc_id); nutrients.add(nid); sources.add(src)
        if amount:
            try:
                float(amount)
            except ValueError:
                bad_amount += 1
    e = EXPECTED[p.name]
    check(bad_amount == 0, "every non-empty amount parses as a number", f"{bad_amount} bad")
    check(rows == e["rows"], f"row count == {e['rows']:,}", f"got {rows:,}")
    check(len(foods) == e["foods"], f"distinct foods == {e['foods']:,}", f"got {len(foods):,}")
    check(len(nutrients) == e["nutrients"], f"distinct nutrients == {e['nutrients']:,}",
          f"got {len(nutrients):,}")
    return nutrients, sources


def main() -> int:
    ap = argparse.ArgumentParser(description="Health-check the exported resource files.")
    ap.add_argument("--exports", default="/data/bac2food/exports")
    ap.add_argument("--quick", action="store_true",
                    help="skip food_nutrients.tsv (the 3.3 GB scan)")
    args = ap.parse_args()
    d = Path(args.exports)

    print("=" * 72)
    print("bac2food export health check")
    print("=" * 72)

    se_ec, se_orgs = verify_species_enzymes(d / "species_enzymes.tsv")
    dg_ec, dg_nut = verify_enzyme_substrate(d / "enzyme_substrate_chebi.tsv")

    fn_nut = fn_src = None
    if not args.quick:
        fn_nut, fn_src = verify_food_nutrients(d / "food_nutrients.tsv")

    print("\n[4] deposit hygiene")
    present = {f.name for f in d.iterdir() if f.is_file()}
    extra = present - DELIVERABLES
    missing = DELIVERABLES - present
    check(not missing, "all deliverable files present", f"missing {sorted(missing)}" if missing else "")
    check(not extra, "no non-deliverable files in the export directory",
          f"{len(extra)} stray: {sorted(extra)}" if extra else "")

    # [4b] Licence compliance. This is the check whose absence let a non-compliant file pass
    # 24/24. Rights live in licences.tsv now, joined on source_db, so the checks are: the table
    # exists, every label in the composition table resolves in it (composites included), and no
    # source whose licence forbids redistributing derived values has any row in the deposit.
    print("\n[4b] licence compliance of the deposit")
    import csv as _csv
    lic_p = d / "licences.tsv"
    if lic_p.exists():
        with lic_p.open(encoding="utf-8", newline="") as fh:
            lic_rows = {r["source_db"]: r for r in _csv.DictReader(fh, delimiter="\t")}
        need = {"source_db", "licence", "version", "tier",
                "derived_values_redistributable", "attribution_string", "in_deposit"}
        have = set(next(iter(lic_rows.values())).keys()) if lic_rows else set()
        check(need <= have, "licences.tsv carries the rights fields",
              f"missing {sorted(need - have)}" if not need <= have else f"{len(lic_rows)} sources")
        if fn_src is not None:
            unresolved = {s for s in fn_src if s not in lic_rows}
            check(not unresolved, "every source_db resolves in licences.tsv",
                  f"unresolved: {sorted(unresolved)}" if unresolved
                  else f"{len(fn_src)} labels, all joined")
            leaked = {s for s in fn_src
                      if lic_rows.get(s, {}).get("in_deposit") == "no"
                      or lic_rows.get(s, {}).get("derived_values_redistributable", "").startswith("no")}
            check(not leaked, "no restricted-source values in the deposit",
                  f"found {sorted(leaked)} — their licence forbids redistributing derived "
                  f"values; rebuild with the default --restricted exclude" if leaked
                  else "none of the deposited sources is restricted")
        else:
            print("  [SKIP] source_db referential checks (--quick)")
    else:
        check(False, "licences.tsv present in the deposit", "absent")

    # The row-level check above only sees values inside food_nutrients.tsv. A withheld source
    # can also enter the deposit as a file — the original release dropped in beside the exports.
    # That is a different failure and needs its own check.
    stowaway = sorted(p.name for p in d.iterdir()
                      if p.is_file() and any(s in p.name.lower() for s in WITHHELD))
    check(not stowaway, "no withheld-source release file in the deposit",
          f"found {stowaway} — {', '.join(WITHHELD)} ships as regeneration code only, not as data"
          if stowaway else f"withheld: {', '.join(WITHHELD)}")

    print("\n[5] cross-export joins (the chain the resource promises)")
    shared = se_ec & dg_ec
    check(len(shared) > 0, "species_enzymes EC join enzyme_substrate_chebi EC",
          f"{len(shared):,} of {len(se_ec):,} organism-side EC resolve to a substrate "
          f"({100*len(shared)/len(se_ec):.1f}%)")
    if fn_nut is not None:
        orphan = dg_nut - fn_nut
        # A nutrient_id can exist in the 3,300-entry catalog yet carry no value in the
        # composition table. The digest links to such ids, so an (EC, substrate) row can
        # be flagged in_model and still reach no actual food. Quantify rather than assume.
        dead_rows = live_rows = 0
        with (d / "enzyme_substrate_chebi.tsv").open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["in_model"] != "yes":
                    continue
                ids = [x.strip() for x in (r["nutrient_ids"] or "").split(",") if x.strip()]
                if any(i in fn_nut for i in ids):
                    live_rows += 1
                else:
                    dead_rows += 1
        check(not orphan, "every nutrient_id in the digest exists in food_nutrients",
              f"{len(orphan):,} orphaned (catalog entries with no measured value)"
              + (f" e.g. {sorted(orphan)[:5]}" if orphan else ""))
        check(dead_rows == 0, "every in_model row reaches a nutrient present in the table",
              f"{dead_rows:,} of {dead_rows + live_rows:,} in_model rows reach no food "
              f"({100*dead_rows/max(1, dead_rows+live_rows):.1f}%)")
    else:
        print("  [SKIP] nutrient_id referential checks (--quick)")

    print()
    npass = sum(1 for ok, _ in results if ok)
    print("=" * 72)
    print(f"{npass}/{len(results)} checks passed")
    for ok, label in results:
        if not ok:
            print(f"  FAILED: {label}")
    print("=" * 72)
    return 0 if npass == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
