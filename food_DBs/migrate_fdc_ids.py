#!/usr/bin/env python3
"""migrate_fdc_ids.py — one-time migration from offset-encoded fdc_ids to block accessions.

Builds `fdc_id_map.tsv`, the frozen accession registry every ingest reads from after this
runs, and (with --apply) rewrites the store to use it.

    python migrate_fdc_ids.py --build            # write the map + validation report, touch nothing
    python migrate_fdc_ids.py --apply --out DIR  # rewrite food.parquet + bucketed store into DIR

WHY THE NEW IDS ARE ORDER-PRESERVING
------------------------------------
Within a source, the accession is assigned by ascending native code, so the *order* of a
source's ids is unchanged. Across sources, blocks keep their old relative sequence. Both
matter because the predictor elects a representative food per canonical group with

    max(fids, key=lambda f: (food_stats[f]["tp"], -f))     # highest type-priority, lowest fid

and nearly every non-USDA food is `foundation_food`, so `tp` ties and the tie-break resolves
to block order. An order-preserving remap therefore elects the same representative for every
group, which is what keeps scores, rankings and food names identical to the released run.
`--build` asserts this rather than assuming it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fdc_blocks as B  # noqa: E402

DEFAULT_FOOD = "/data/bac2food/food.parquet"
DEFAULT_BUCKETED = "/data/bac2food/food_nutrient_bucketed"


def decode_source(fdc_ids: pd.Series) -> pd.Series:
    """Assign each id to a source using the LEGACY boundaries (pre-migration layout)."""
    bounds = [(0, B.FDC_RESERVED_BLOCKS * B.BLOCK_WIDTH, B.FDC_KEY)]
    legacy = sorted(((B.legacy_base(k), k) for _, k, _, _ in B.SOURCES))
    for n, (lo, key) in enumerate(legacy):
        hi = legacy[n + 1][0] if n + 1 < len(legacy) else 1 << 62
        bounds.append((lo, hi, key))
    out = pd.Series(pd.NA, index=fdc_ids.index, dtype="object")
    for lo, hi, key in bounds:
        out[(fdc_ids >= lo) & (fdc_ids < hi)] = key
    return out


def build_map(food_path: Path) -> pd.DataFrame:
    f = pd.read_parquet(food_path, columns=["fdc_id", "description", "data_type"])
    f["source_key"] = decode_source(f["fdc_id"])
    if f["source_key"].isna().any():
        raise SystemExit(f"[!] {int(f['source_key'].isna().sum())} ids match no legacy block")

    parts = []
    for key in [B.FDC_KEY] + [k for _, k, _, _ in B.SOURCES]:
        s = f[f["source_key"] == key].copy()
        if s.empty:
            continue
        s["source_food_code"] = s["fdc_id"] - B.legacy_base(key)
        if key == B.FDC_KEY:
            # USDA ids are the upstream authority: never reassigned.
            s["fdc_id_new"] = s["fdc_id"]
        else:
            s = s.sort_values("source_food_code", kind="mergesort")
            if s["source_food_code"].duplicated().any():
                raise SystemExit(f"[!] {key}: duplicate native codes; cannot key on them")
            s["fdc_id_new"] = B.base(key) + pd.RangeIndex(len(s))
            over = int((s["fdc_id_new"] >= B.limit(key)).sum())
            if over:
                raise SystemExit(f"[!] {key}: {over} accessions overflow a "
                                 f"{B.BLOCK_WIDTH:,}-wide block")
        parts.append(s[["source_key", "source_food_code", "fdc_id", "fdc_id_new",
                        "data_type", "description"]])
    m = pd.concat(parts, ignore_index=True)
    return m.rename(columns={"fdc_id": "fdc_id_legacy"})


def validate(m: pd.DataFrame) -> list[str]:
    """Every property the migration has to preserve, checked rather than assumed."""
    errs: list[str] = []
    ext = m[m["source_key"] != B.FDC_KEY]

    if ext["fdc_id_new"].duplicated().any():
        errs.append(f"{int(ext['fdc_id_new'].duplicated().sum())} duplicate new accessions")
    if m["fdc_id_new"].duplicated().any():
        errs.append("new accession collides with a USDA id")

    # Containment: no id outside its own block.
    ids_by_source = {k: set(g["fdc_id_new"]) for k, g in m.groupby("source_key")}
    try:
        B.assert_disjoint(ids_by_source)
    except AssertionError as e:
        errs.append(str(e))

    # Order preservation within each source -- the property representative election needs.
    for key, g in ext.groupby("source_key"):
        g = g.sort_values("fdc_id_legacy")
        if not g["fdc_id_new"].is_monotonic_increasing:
            errs.append(f"{key}: remap is not order-preserving within the source")

    # Order preservation across sources: legacy block sequence must equal new block sequence.
    legacy_order = [k for _, k in sorted((B.legacy_base(k), k) for _, k, _, _ in B.SOURCES)]
    new_order = [k for _, k in sorted((B.base(k), k) for _, k, _, _ in B.SOURCES)]
    if legacy_order != new_order:
        errs.append(f"block sequence changed: {legacy_order} -> {new_order}")

    # source_of() must round-trip every new id back to the source that owns it.
    bad = [k for k, g in m.groupby("source_key")
           if any(B.source_of(i) != k for i in g["fdc_id_new"].head(2000))]
    if bad:
        errs.append(f"source_of() misroutes: {bad}")
    return errs


def report(m: pd.DataFrame) -> None:
    ext = m[m["source_key"] != B.FDC_KEY]
    print(f"\n{'source':<18}{'foods':>9}{'block base':>13}{'last used':>13}"
          f"{'occupancy':>11}  {'ids kept'}")
    print("-" * 78)
    for _, key, lab, _ in [(0, B.FDC_KEY, B.FDC_LABEL, 0)] + list(B.SOURCES):
        g = m[m["source_key"] == key]
        if g.empty:
            continue
        kept = int((g["fdc_id_legacy"] == g["fdc_id_new"]).sum())
        occ = "" if key == B.FDC_KEY else f"{len(g) / B.BLOCK_WIDTH:>10.2%}"
        print(f"{lab[:17]:<18}{len(g):>9,}{B.base(key):>13,}"
              f"{int(g['fdc_id_new'].max()):>13,}{occ:>11}  "
              f"{kept:,}{' (all)' if kept == len(g) else ''}")
    print("-" * 78)
    print(f"{'total':<18}{len(m):>9,}   changed: {int((m['fdc_id_legacy'] != m['fdc_id_new']).sum()):,}"
          f"   unchanged: {int((m['fdc_id_legacy'] == m['fdc_id_new']).sum()):,}")
    print(f"\nold span 0..{int(m['fdc_id_legacy'].max()):,}   "
          f"new span 0..{int(m['fdc_id_new'].max()):,}   "
          f"({1 - m['fdc_id_new'].max() / m['fdc_id_legacy'].max():.1%} narrower)")


def apply_store(m: pd.DataFrame, food_path: Path, bucketed: Path, out: Path) -> None:
    """Rewrite food.parquet and the bucketed values into `out`. Never writes in place.

    Only the 2% of ids that actually move are carried in the remap dict, so the value store
    is rewritten bucket by bucket with a partial map rather than a 28M-row join. Bucketing is
    on `nutrient_id % 256`, so no row changes partition and each file is independent.
    """
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq

    out.mkdir(parents=True, exist_ok=True)
    moved = m[m["fdc_id_legacy"] != m["fdc_id_new"]]
    remap = dict(zip(moved["fdc_id_legacy"].tolist(), moved["fdc_id_new"].tolist()))
    print(f"\n[*] remapping {len(remap):,} moved ids "
          f"({len(m) - len(remap):,} keep their accession)")

    # --- food.parquet: remap the key, and record provenance the old scheme only implied ---
    f = pd.read_parquet(food_path)
    key = m.set_index("fdc_id_legacy")
    f["source_db"] = f["fdc_id"].map(key["source_key"])
    f["source_food_code"] = f["fdc_id"].map(key["source_food_code"])
    f["fdc_id"] = f["fdc_id"].map(lambda i: remap.get(i, i))
    f = f.sort_values("fdc_id", kind="mergesort").reset_index(drop=True)
    f.to_parquet(out / "food.parquet", index=False)
    print(f"[*] wrote {out / 'food.parquet'}  ({len(f):,} foods, "
          f"+source_db +source_food_code)")

    # --- value store: one pass per bucket ---
    dest = out / "food_nutrient_bucketed"
    dest.mkdir(exist_ok=True)
    dataset = ds.dataset(str(bucketed), format="parquet", partitioning="hive")
    n_rows = n_moved = 0
    for frag in dataset.get_fragments():
        t = frag.to_table()
        # The partition column is encoded in the path, not in the fragment's own schema.
        b = int(re.search(r"bucket=(\d+)", frag.path).group(1))
        ids = t.column("fdc_id").to_pandas()
        hit = ids.isin(remap)
        n_rows += len(ids)
        n_moved += int(hit.sum())
        if hit.any():
            ids = ids.mask(hit, ids[hit].map(remap))
            t = t.set_column(t.schema.get_field_index("fdc_id"), "fdc_id",
                             pa.array(ids, type=pa.int64()))
        # One file PER FRAGMENT, keeping its original basename. Each bucket holds several
        # files (one per source DB, plus compaction output); writing them all to a single
        # part-0.parquet silently keeps only the last and drops the rest.
        d = dest / f"bucket={b}"
        d.mkdir(exist_ok=True)
        keep = [n for n in t.schema.names if n != "bucket"]
        pq.write_table(t.select(keep), d / Path(frag.path).name)
    # A migration that drops rows still produces a valid-looking store, so count both sides.
    written = ds.dataset(str(dest), format="parquet", partitioning="hive").count_rows()
    if written != n_rows:
        raise AssertionError(
            f"value store lost rows: read {n_rows:,}, wrote {written:,} "
            f"({n_rows - written:,} missing)")
    print(f"[*] wrote {dest}  ({n_rows:,} values, {n_moved:,} remapped, count verified)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--food", default=DEFAULT_FOOD)
    ap.add_argument("--bucketed", default=DEFAULT_BUCKETED)
    ap.add_argument("--map", default="/data/bac2food/rekey/fdc_id_map.tsv")
    ap.add_argument("--build", action="store_true", help="write the map, change nothing else")
    ap.add_argument("--apply", action="store_true", help="rewrite the store into --out")
    ap.add_argument("--out", default="/data/bac2food/rekey/store")
    args = ap.parse_args()

    m = build_map(Path(args.food))
    errs = validate(m)
    report(m)

    print("\nvalidation")
    print("-" * 78)
    if errs:
        for e in errs:
            print(f"  FAIL  {e}")
        return 1
    print("  ok  accessions unique and inside their blocks")
    print("  ok  order preserved within every source")
    print("  ok  block sequence unchanged -> representative election unchanged")
    print("  ok  source_of() round-trips")

    Path(args.map).parent.mkdir(parents=True, exist_ok=True)
    m[["source_key", "source_food_code", "fdc_id_new", "fdc_id_legacy"]] \
        .rename(columns={"fdc_id_new": "fdc_id"}) \
        .to_csv(args.map, sep="\t", index=False)
    print(f"\n[*] wrote {args.map}  ({len(m):,} accessions)")

    if args.apply:
        apply_store(m, Path(args.food), Path(args.bucketed), Path(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
