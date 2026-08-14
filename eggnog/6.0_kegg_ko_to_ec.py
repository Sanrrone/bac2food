#!/usr/bin/env python3
"""6.0_kegg_ko_to_ec.py — build the KEGG KO -> EC bridge.

eggNOG v7 dropped the direct EC annotation that v6 carried: e7.og_info_kegg_go.tsv
annotates orthologous groups with KEGG Orthology (KO) ids, not EC numbers. To keep
bac2food's bacterium -> EC layer we need a KO -> EC map, and KEGG publishes one for
free inside the KO list itself:

    K00001  ADH; alcohol dehydrogenase [EC:1.1.1.1]
    K00121  frmA, ...; S-(hydroxymethyl)glutathione dehydrogenase [EC:1.1.1.284 1.1.1.1]

so a single ~2 MB request (https://rest.kegg.jp/list/ko) yields the whole bridge.
The per-KO endpoint (/link/ec/K00001) returns the same thing but would need ~28,000
requests; the bulk /link/ec/ko endpoint is not served (HTTP 400).

Only ~39% of KOs are enzymes; the rest (transporters, structural proteins,
regulators) carry no EC and are dropped here -- that loss is inherent to the bridge,
not a parsing bug, so it is reported explicitly.

Partial ECs ("1.1.1.-", an incompletely characterised activity) cannot join the
EC -> substrate digest, which is keyed on 4-level ECs. They are written to the output
with a `partial` flag so the loss stays auditable; --full-only drops them instead.

Usage:
    python 6.0_kegg_ko_to_ec.py --out /data/bac2food/kegg_ko_ec.tsv
    python 6.0_kegg_ko_to_ec.py --ko_list /data/bac2food/kegg_ko_list.tsv --out ...
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import urllib.request
from pathlib import Path

KEGG_KO_LIST = "https://rest.kegg.jp/list/ko"

# "[EC:1.1.1.1 2.7.1.-]" -> one bracket, one-or-more space-separated ECs.
_EC_BLOCK_RE = re.compile(r"\[EC:([^\]]+)\]")
# A single EC token: 4 dot-separated levels, any of which may be "-" (partial).
_EC_TOKEN_RE = re.compile(r"\d+(?:\.(?:\d+|-)){3}")


def fetch_ko_list(path: Path | None, timeout: int = 180) -> list[str]:
    """Return the KEGG KO list as lines, from a cached file or the REST API."""
    if path and path.exists():
        print(f"[*] Reading cached KO list: {path}", flush=True)
        return path.read_text(encoding="utf-8").splitlines()

    print(f"[*] Downloading {KEGG_KO_LIST} ...", flush=True)
    with urllib.request.urlopen(KEGG_KO_LIST, timeout=timeout) as r:
        text = r.read().decode("utf-8")
    if path:
        path.write_text(text, encoding="utf-8")
        print(f"    cached to {path}", flush=True)
    return text.splitlines()


def parse_ko_ec(lines: list[str]) -> tuple[list[tuple[str, str, str]], dict]:
    """Parse KO list lines into (ko, ec_number, is_partial) rows + stats.

    A KO may carry several ECs (multifunctional / ambiguous assignments), so this is
    a genuine many-to-many bridge: one KO can fan out to N ECs, and one EC can be
    reached from several KOs.
    """
    rows: list[tuple[str, str, str]] = []
    seen: set[tuple[str, str]] = set()
    n_ko = n_enzyme = 0

    for line in lines:
        line = line.rstrip("\n")
        if not line:
            continue
        ko, _, definition = line.partition("\t")
        ko = ko.strip()
        if not ko.startswith("K"):
            continue
        n_ko += 1

        ecs: list[str] = []
        for block in _EC_BLOCK_RE.findall(definition):
            ecs += _EC_TOKEN_RE.findall(block)
        if not ecs:
            continue
        n_enzyme += 1

        for ec in ecs:
            key = (ko, ec)
            if key in seen:
                continue
            seen.add(key)
            rows.append((ko, ec, "yes" if "-" in ec else "no"))

    stats = {
        "kos_total": n_ko,
        "kos_with_ec": n_enzyme,
        "pairs": len(rows),
        "pairs_full": sum(1 for r in rows if r[2] == "no"),
        "ecs_full": len({r[1] for r in rows if r[2] == "no"}),
    }
    return rows, stats


def main() -> int:
    ap = argparse.ArgumentParser(description="Build the KEGG KO -> EC bridge.")
    ap.add_argument("--out", default="/data/bac2food/kegg_ko_ec.tsv",
                    help="Output TSV (ko, ec_number, partial)")
    ap.add_argument("--ko_list", default="/data/bac2food/kegg_ko_list.tsv",
                    help="Cache path for the raw KEGG KO list; downloaded if absent")
    ap.add_argument("--full-only", action="store_true",
                    help="Drop partial ECs (1.1.1.-) instead of flagging them")
    args = ap.parse_args()

    lines = fetch_ko_list(Path(args.ko_list) if args.ko_list else None)
    rows, st = parse_ko_ec(lines)
    if args.full_only:
        rows = [r for r in rows if r[2] == "no"]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["ko", "ec_number", "partial"])
        w.writerows(sorted(rows))

    pct = 100.0 * st["kos_with_ec"] / st["kos_total"] if st["kos_total"] else 0.0
    print(f"[*] {out}: {len(rows):,} (KO, EC) pairs", flush=True)
    print(f"    KOs in KEGG        : {st['kos_total']:,}", flush=True)
    print(f"    KOs with an EC     : {st['kos_with_ec']:,} ({pct:.1f}%) "
          f"-- the rest are non-enzymatic and cannot reach an EC", flush=True)
    print(f"    pairs w/ full EC   : {st['pairs_full']:,} "
          f"-> {st['ecs_full']:,} distinct 4-level ECs", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
