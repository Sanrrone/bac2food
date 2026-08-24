#!/usr/bin/env python3
"""prune_bucketed_store.py - drop from the bucketed store every nutrient no enzyme acts on.

The store is what the PREDICTOR reads; food_nutrients.tsv is what the deposit ships. They
have to agree about what exists, or the scorer can rank a food the published table does not
contain - the failure mode that bit this project once already. export_resources.py applies
the chain filter with --span_chain; this applies the same filter, derived from the same
3_nutrient_to_ec.tsv, to the store behind it.

The keep-set is NOT just the EC targets: it carries the generics the scoring kernel
substitutes for a missing target. Dropping those 69 cut one sample's candidate pool by
46%. See _common/chain_filter.py, which export_resources.py shares.

Filtering on NUTRIENT alone is deliberate and sufficient. A food whose every value is an
unlinked nutrient ends up with zero rows and disappears on its own, so no food-level
exclusion policy (branded, modelled, re-listings) has to be restated here - restating it
is how the two readers drift apart.

Layout is preserved exactly, including per-source filenames: source_of_bucket_file() reads
the SOURCE out of the filename, so renaming or merging these files silently destroys the
source_db column downstream.

    python3 prune_bucketed_store.py --in_dir  /data/bac2food/food_nutrient_bucketed \
                                    --out_dir /data/bac2food/food_nutrient_bucketed_pruned
"""
from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common.chain_filter import chain_nutrients, chain_targets  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
DEFAULT_MAP = REPO / "0_building" / "3_nutrient_to_ec.tsv"
DEFAULT_ALIAS = REPO / "0_building" / "1_expanded_nutrients.tsv"




def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--nutrient_to_ec", default=str(DEFAULT_MAP))
    ap.add_argument("--nutrient_alias", default=str(DEFAULT_ALIAS))
    ap.add_argument("--force", action="store_true", help="overwrite out_dir if it exists")
    args = ap.parse_args()

    src, dst = Path(args.in_dir), Path(args.out_dir)
    if dst.exists():
        if not args.force:
            raise SystemExit(f"ERROR: {dst} exists; pass --force to replace it.")
        shutil.rmtree(dst)

    keep = chain_nutrients(Path(args.nutrient_to_ec), Path(args.nutrient_alias))
    _vs: dict = {}

    def keep_arr(typ):
        """The value_set, built once per arrow type rather than once per file."""
        if typ not in _vs:
            _vs[typ] = pa.array(sorted(keep), type=typ)
        return _vs[typ]

    n_t = len(chain_targets(Path(args.nutrient_to_ec)))
    print(f"[*] {n_t:,} EC targets + {len(keep)-n_t} kernel substitutes = "
          f"{len(keep):,} nutrients kept")

    files = sorted(src.glob("bucket=*/*.parquet"))
    if not files:
        raise SystemExit(f"ERROR: no bucket=*/*.parquet under {src}")

    tot_in = tot_out = 0
    kept_files = dropped_files = 0
    for i, f in enumerate(files, 1):
        t = pq.read_table(f)
        tot_in += t.num_rows
        t = t.filter(pc.is_in(t["nutrient_id"], value_set=keep_arr(t["nutrient_id"].type)))
        tot_out += t.num_rows
        if t.num_rows == 0:
            dropped_files += 1                     # an all-unlinked source in this bucket
        else:
            out = dst / f.parent.name / f.name     # bucket=N/<source>_data.parquet
            out.parent.mkdir(parents=True, exist_ok=True)
            pq.write_table(t, out, compression="zstd")
            kept_files += 1
        if i % 200 == 0 or i == len(files):
            print(f"    {i}/{len(files)} files | {tot_out:,}/{tot_in:,} rows kept", flush=True)

    print(f"[*] {src} -> {dst}")
    print(f"[*] rows  {tot_in:,} -> {tot_out:,}  ({100*tot_out/tot_in:.1f}% kept)")
    print(f"[*] files {len(files):,} -> {kept_files:,}  ({dropped_files:,} became empty)")
    a = sum(p.stat().st_size for p in src.rglob("*.parquet"))
    b = sum(p.stat().st_size for p in dst.rglob("*.parquet"))
    print(f"[*] size  {a/1e6:.1f} MB -> {b/1e6:.1f} MB")


if __name__ == "__main__":
    main()
