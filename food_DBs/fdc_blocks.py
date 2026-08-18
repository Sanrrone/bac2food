#!/usr/bin/env python3
"""fdc_blocks.py — the single definition of how `fdc_id` is allocated.

Every ingest script imports from here. Nothing mints an `fdc_id` by writing a literal
offset, because that is how the two collisions in this repository's history happened:

  * Swedish foods were minted at 40,000,000, the block McCance already fills positionally.
    fdc_id 40000001 was "Agar, dried" (McCance) and "Beef tallow" (Swedish) at the same time.
    1,788 Swedish records were overwritten, their nutrient values left on the McCance food.
    (food_DBs/livsmedels_sweeden/ingest_swedish_livsmedelsdb.py)
  * WAFCT was minted at 94,000,000, which sits *inside* PhyFoodComp's range: PhyFoodComp adds
    a native code reaching 19,020,060 to a 92M base, so it sprawls from 93M to 111M. 27 WAFCT
    foods landed on PhyFoodComp ones. (food_DBs/WAFT_AFRICA/injest_wafct.py)

Both were silent. A shared id merges two unrelated foods; nothing raises.

WHY BLOCKS ARE UNIFORM AND WHY IDS ARE ACCESSIONS
-------------------------------------------------
The old scheme set `fdc_id = block_base + native_source_code`. That made the id reversible —
subtract the base, recover the code printed in the source's own table — but it tied the block
WIDTH to how sparsely a source numbers its foods rather than to how many foods it has:

    source           foods    max native code    density
    PhyFoodComp      3,377         19,020,060      0.02%   <- needs a 20M block
    BioFoodComp     10,133          1,200,005      0.84%
    CNF              5,690            503,380      1.13%
    WAFCT            1,028             14,034      7.32%

So blocks had to be 10M+ and were still overrun. Instead, `fdc_id` is now an **accession**:
opaque, dense, assigned once, and never recomputed. The source's own identifier is no longer
encoded in the arithmetic — it is carried explicitly in the `source_food_code` column, which
is strictly more honest, because a reader no longer has to know the base to recover it.

Density then follows the food count, not the code space. The largest source needs 10,133
slots, so a 3,000,000-wide block runs at 0.34% occupancy: ~296x headroom.

STABILITY
---------
Accessions are frozen in `fdc_id_map.tsv` (source_key, source_food_code, fdc_id). Ingests look
the code up there; only codes absent from the map are assigned, and they are appended at the
high-water mark of their block. A source republishing with new or renumbered foods therefore
cannot disturb the ids of foods already released — which a bare positional `range()` could not
promise, and three of the old ingests (Phenol-Explorer, McCance, AFCD, Swiss) quietly did not.

LAYOUT
------
    block index  =  fdc_id // 3_000_000

    blocks 0-2   0          .. 8,999,999    USDA FoodData Central, native ids, NOT reassigned
    block  3     9,000,000  .. 11,999,999   Phenol-Explorer
    block  4     12,000,000 ..              Fineli
    ...          one block per source, in ingestion order

USDA keeps blocks 0-2 because its ids are the upstream authority and must stay untouched; it
currently reaches 2,751,503, so the reservation carries 3.3x headroom. `assert_disjoint()`
enforces the boundary at export, so an FDC release crossing 9M fails loudly instead of
colliding with Phenol-Explorer.

Source order is the order the sources were ingested, and is deliberately NOT re-sorted. The
predictor elects one representative fdc_id per canonical food group by highest type-priority,
tie-broken on lowest fdc_id (bac2food_predict.py, `rep_per_canon`). Nearly every non-USDA food
carries data_type `foundation_food`, so that tie-break resolves to *block order*. Preserving
the sequence keeps representative election — and therefore every score, ranking and food name
in the outputs — identical to the released version.
"""
from __future__ import annotations

import os
from pathlib import Path

BLOCK_WIDTH = 3_000_000

# Blocks reserved for USDA FoodData Central's own identifiers. FDC ids are upstream and are
# never reassigned; this range only has to be wide enough that FDC cannot grow into block 3.
FDC_RESERVED_BLOCKS = 3

# (block_index, key, label, legacy_base) in ingestion order.
# `legacy_base` is the offset the pre-migration scheme used; it is kept ONLY so the migration
# can decode an old id back to its native source code, and is not used to mint anything.
SOURCES: tuple[tuple[int, str, str, int], ...] = (
    (3,  "phenol_explorer", "Phenol-Explorer",                  9_000_000),
    (4,  "fineli",          "Fineli",                          10_000_000),
    (5,  "ciqual",          "CIQUAL",                          20_000_000),
    (6,  "nevo",            "NEVO-online",                     30_000_000),
    (7,  "mccance",         "McCance & Widdowson's",           40_000_000),
    (8,  "swedish",         "Livsmedelsdatabasen",             41_000_000),
    (9,  "frida",           "Frida",                           50_000_000),
    (10, "cnf",             "Canadian Nutrient File",          60_000_000),
    (11, "afcd",            "AFCD",                            70_000_000),
    (12, "swiss",           "Swiss Food Composition Database",  80_000_000),
    (13, "stfcj",           "STFCJ",                           81_000_000),
    (14, "biofoodcomp",     "FAO BioFoodComp",                 82_000_000),
    (15, "phyfoodcomp",     "FAO PhyFoodComp",                 92_000_000),
    (16, "wafct",           "FAO WAFCT",                      120_000_000),
)

FDC_KEY = "fdc"
FDC_LABEL = "USDA FoodData Central"

_BY_KEY = {k: (i, lab, legacy) for i, k, lab, legacy in SOURCES}
_BY_BLOCK = {i: (k, lab, legacy) for i, k, lab, legacy in SOURCES}


def base(source_key: str) -> int:
    """First fdc_id of a source's block."""
    if source_key == FDC_KEY:
        return 0
    try:
        return _BY_KEY[source_key][0] * BLOCK_WIDTH
    except KeyError:
        raise KeyError(f"unknown source key {source_key!r}; known: {sorted(_BY_KEY)}") from None


def limit(source_key: str) -> int:
    """One past the last fdc_id of a source's block."""
    if source_key == FDC_KEY:
        return FDC_RESERVED_BLOCKS * BLOCK_WIDTH
    return base(source_key) + BLOCK_WIDTH


def legacy_base(source_key: str) -> int:
    """The offset the pre-migration scheme used. Migration only."""
    return 0 if source_key == FDC_KEY else _BY_KEY[source_key][2]


def label(source_key: str) -> str:
    return FDC_LABEL if source_key == FDC_KEY else _BY_KEY[source_key][1]


def source_of(fdc_id: int) -> str:
    """Which source an fdc_id belongs to. The whole point of a uniform block width."""
    b = int(fdc_id) // BLOCK_WIDTH
    if b < FDC_RESERVED_BLOCKS:
        return FDC_KEY
    if b not in _BY_BLOCK:
        raise ValueError(f"fdc_id {fdc_id} falls in unallocated block {b} "
                         f"({b * BLOCK_WIDTH:,}..{(b + 1) * BLOCK_WIDTH - 1:,})")
    return _BY_BLOCK[b][0]


def assert_disjoint(ids_by_source: dict[str, set[int]]) -> None:
    """Fail loudly if any source's ids leave its block or overlap another's.

    Called at export. The historical collisions were silent overwrites, so the check is on
    containment, not on whether a duplicate happens to be visible in the merged table.
    """
    problems: list[str] = []
    for key, ids in ids_by_source.items():
        if not ids:
            continue
        lo, hi = base(key), limit(key)
        out = [i for i in ids if not (lo <= i < hi)]
        if out:
            problems.append(
                f"{key}: {len(out)} id(s) outside block [{lo:,}, {hi:,}) "
                f"e.g. {sorted(out)[:3]}")
        used = hi - lo
        if len(ids) > used:
            problems.append(f"{key}: {len(ids):,} ids exceed block capacity {used:,}")
    seen: dict[int, str] = {}
    for key, ids in ids_by_source.items():
        for i in ids:
            if i in seen and seen[i] != key:
                problems.append(f"id {i} claimed by both {seen[i]} and {key}")
                break
            seen[i] = key
    if problems:
        raise AssertionError("fdc_id block allocation is invalid:\n  " + "\n  ".join(problems))


ACCESSION_MAP = Path(os.environ.get("BAC2FOOD_FDC_MAP", "/data/bac2food/fdc_id_map.tsv"))


def load_map(path: Path | None = None) -> dict[tuple[str, str], int]:
    """Read the frozen accession registry as {(source_key, source_food_code): fdc_id}.

    Codes are compared as strings. Sources key their foods with everything from plain
    plain integers to zero-padded composites (WAFCT '01_172'), and int() would collapse
    '007' and '7' onto one accession.
    """
    p = Path(path) if path else ACCESSION_MAP
    if not p.exists():
        return {}
    out: dict[tuple[str, str], int] = {}
    with p.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_src, i_code, i_id = (header.index(c) for c in
                               ("source_key", "source_food_code", "fdc_id"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[(f[i_src], f[i_code])] = int(f[i_id])
    return out


def assign(source_key: str, codes, path: Path | None = None,
           persist: bool = True) -> list[int]:
    """Map a source's own food codes to accessions, minting only what is genuinely new.

    Existing codes keep the accession they were first given, which is what makes ids stable
    across source releases -- the property a bare positional `range()` cannot offer, and that
    four of the old ingests silently lacked. New codes are appended above the block's current
    high-water mark in ascending code order, so a re-run of the same input is deterministic.

    Minting is recorded: new accessions are appended to the registry unless persist=False.
    """
    p = Path(path) if path else ACCESSION_MAP
    reg = load_map(p)
    codes = [str(c).strip() for c in codes]

    lo, hi = base(source_key), limit(source_key)
    used = [v for (s, _), v in reg.items() if s == source_key]
    nxt = (max(used) + 1) if used else lo

    fresh = sorted({c for c in codes if (source_key, c) not in reg})
    minted: list[tuple[str, str, int]] = []
    for c in fresh:
        if nxt >= hi:
            raise OverflowError(
                f"{source_key}: block [{lo:,}, {hi:,}) is full at {len(used) + len(minted):,} "
                f"foods; widen BLOCK_WIDTH and re-key, do not spill into the next source")
        reg[(source_key, c)] = nxt
        minted.append((source_key, c, nxt))
        nxt += 1

    if minted and persist:
        new_file = not p.exists()
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("a", encoding="utf-8") as fh:
            if new_file:
                fh.write("source_key\tsource_food_code\tfdc_id\tfdc_id_legacy\n")
            for s, c, v in minted:
                fh.write(f"{s}\t{c}\t{v}\t\n")
        print(f"[fdc_blocks] {source_key}: minted {len(minted):,} new accession(s) "
              f"-> {p} (block now {len(used) + len(minted):,}/{hi - lo:,})")

    return [reg[(source_key, c)] for c in codes]


def describe() -> str:
    rows = [f"{'block':>5}  {'range':>27}  {'source':<32}",
            f"{'-' * 5}  {'-' * 27}  {'-' * 32}"]
    rows.append(f"{'0-' + str(FDC_RESERVED_BLOCKS - 1):>5}  "
                f"{0:>12,} .. {FDC_RESERVED_BLOCKS * BLOCK_WIDTH - 1:>11,}  {FDC_LABEL:<32}")
    for i, k, lab, _ in SOURCES:
        rows.append(f"{i:>5}  {i * BLOCK_WIDTH:>12,} .. "
                    f"{(i + 1) * BLOCK_WIDTH - 1:>11,}  {lab:<32}")
    return "\n".join(rows)


if __name__ == "__main__":
    print(describe())
    print(f"\nblock width {BLOCK_WIDTH:,}   allocated through "
          f"{(SOURCES[-1][0] + 1) * BLOCK_WIDTH - 1:,}")
