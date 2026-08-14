#!/usr/bin/env python3
"""
nutrient2ec.py

Map FoodData Central nutrients to EC numbers using ChEBI ontology.

Key Upgrades:
- STRICT PREBIOTIC FILTERING: Prevents host-absorbed nutrients (simple sugars, free aminos) 
  and ubiquitous solvents (water) from generating false-positive microbiome targets.
- Preserves the optimized max_cost=1.5 graph traversal.
"""

from __future__ import annotations
import argparse
import csv
import heapq
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional, Set

CHEBI_RE = re.compile(r"^CHEBI:\d+$")
# ChEBI labels a neutral element "<name> atom" and a monatomic ion "<name>(2+)".
_ELEMENT_ENTITY_RE = re.compile(r"^[a-z]+(?:\s+atom|\(\d*[+-]\))$", re.I)

# --- ChEBI RELATION VOCABULARY ---------------------------------------------------------
# chebi.obo declares nine relationship typedefs. Traversing only is_a + the two conjugate
# relations leaves ~40k structural edges on the table, which is the single largest source
# of nutrients that carry a ChEBI id but reach no enzyme. Each entry below is
#   RO/BFO id -> (forward edge label, reverse edge label)
# where "forward" is the direction asserted in the file (the term that DECLARES the
# relationship -> its target).
#
# Direction matters and is not symmetric in cost. Walking a derivative UP to its parent
# scaffold generalizes and is safe: one derivative has one functional parent. Walking a
# parent DOWN to its derivatives specializes and fans out: quercetin alone has dozens.
# The two directions therefore get separate weights (see --w_deriv_up / --w_deriv_down).
#
# has_role (RO:0000087, 45,252 edges - by far the most numerous) is DELIBERATELY EXCLUDED.
# It is a semantic classification ("is an antioxidant", "is a food component"), not a
# structural one; traversing it would connect every polyphenol to every other polyphenol
# through a shared role node and manufacture enzyme links that no chemistry supports.
REL_SYMMETRIC = {
    "RO:0018033": ("conjugate_base_of", "conjugate_acid_of"),   # same species, pH shift
    "RO:0018034": ("conjugate_acid_of", "conjugate_base_of"),
    "RO:0018036": ("tautomer_of",       "tautomer_of"),         # same compound, proton shift
    "RO:0018039": ("enantiomer_of",     "enantiomer_of"),       # same constitution, chirality
}
REL_DIRECTED = {
    # declaring term -> target is the "up"/generalizing direction in every case here.
    "RO:0018038": ("has_functional_parent", "functional_parent_of"),
    "RO:0018040": ("has_parent_hydride",    "parent_hydride_of"),
    "RO:0018037": ("substituent_group_from", "substituent_group_to"),
    "BFO:0000051": ("has_part",             "part_of"),
}
REL_ALL = set(REL_SYMMETRIC) | set(REL_DIRECTED)

# --- BIOLOGICAL FILTERS ---
# These regexes catch nutrients that will never reach the colon or act as noise.
_SIMPLE_SUGARS = re.compile(r"\b(sucrose|glucose|fructose|lactose|maltose|galactose|mannitol|xylitol|sorbitol|sugar)\b", re.I)
_AMINO_ACIDS = re.compile(r"\b(alanine|arginine|asparagine|aspartic acid|cysteine|glutamic acid|glutamine|glycine|histidine|isoleucine|leucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|cystine)\b", re.I)
_UBIQUITOUS_NOISE = re.compile(r"^(water|energy|energy \(atwater.*|ash|solids|nitrogen|proximates)$", re.I)

@dataclass(frozen=True)
class DigestRow:
    ec_number: str
    substrate: str
    chebi_id: str          
    chebi_name: str

@dataclass
class Edge:
    dst: str
    etype: str             
    direction: str         
    weight: float

def parse_chebi_obo(obo_path: str,
                    w_is_a: float = 1.0,
                    w_conj: float = 0.5,
                    w_taut: float = 0.25,
                    w_enant: float = 0.5,
                    w_deriv_up: float = 0.75,
                    w_deriv_down: float = 1.25,
                    w_part: float = 1.0,
                    structural: bool = True,
                    ) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, List[Edge]]]:
    alt2primary: Dict[str, str] = {}
    name_of: Dict[str, str] = {}
    graph: Dict[str, List[Edge]] = defaultdict(list)

    current_id: Optional[str] = None
    current_name: Optional[str] = None
    current_alt_ids: List[str] = []
    current_is_a_parents: List[str] = []
    current_relations: List[Tuple[str, str]] = []  

    def flush_term():
        nonlocal current_id, current_name, current_alt_ids, current_is_a_parents, current_relations
        if not current_id:
            return

        if current_name:
            name_of[current_id] = current_name

        for a in current_alt_ids:
            alt2primary[a] = current_id

        for parent in current_is_a_parents:
            if CHEBI_RE.match(parent):
                graph[current_id].append(Edge(parent, "is_a", "up", w_is_a))
                graph[parent].append(Edge(current_id, "is_a", "down", w_is_a))

        for ro, target in current_relations:
            if not CHEBI_RE.match(target):
                continue
            is_conj = ro in ("RO:0018033", "RO:0018034")
            if ro in REL_SYMMETRIC and (is_conj or structural):
                fwd, rev = REL_SYMMETRIC[ro]
                w = w_conj if is_conj else (w_taut if ro == "RO:0018036" else w_enant)
                graph[current_id].append(Edge(target, fwd, "fwd", w))
                graph[target].append(Edge(current_id, rev, "rev", w))
            elif structural and ro in REL_DIRECTED:
                fwd, rev = REL_DIRECTED[ro]
                # "up" = declaring term -> target: derivative to its scaffold, whole to
                # its part. Cheap. The reverse fans out, so it costs more.
                w_up = w_part if ro == "BFO:0000051" else w_deriv_up
                w_dn = w_part if ro == "BFO:0000051" else w_deriv_down
                graph[current_id].append(Edge(target, fwd, "up", w_up))
                graph[target].append(Edge(current_id, rev, "down", w_dn))

        current_id = None
        current_name = None
        current_alt_ids = []
        current_is_a_parents = []
        current_relations = []

    with open(obo_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            if line == "[Term]":
                flush_term()
                continue
            if not line:
                continue

            if line.startswith("id: "):
                current_id = line.split("id: ", 1)[1].strip()
                continue
            if line.startswith("name: "):
                current_name = line.split("name: ", 1)[1].strip()
                continue
            if line.startswith("alt_id: "):
                current_alt_ids.append(line.split("alt_id: ", 1)[1].strip())
                continue
            if line.startswith("is_a: "):
                parent = line.split("is_a: ", 1)[1].split(" ! ", 1)[0].strip()
                current_is_a_parents.append(parent)
                continue
            if line.startswith("relationship: "):
                rel = line.split("relationship: ", 1)[1].strip()
                parts = rel.split()
                if len(parts) >= 2:
                    ro = parts[0].strip()
                    target = parts[1].strip()
                    if ro in REL_ALL:            # has_role is not in REL_ALL, by design
                        current_relations.append((ro, target))
                continue
            if line in ("[Term]", "[Typedef]"):
                flush_term()
                continue

    flush_term()

    all_nodes = set(graph.keys()) | set(name_of.keys()) | set(alt2primary.values())
    for pid in all_nodes:
        alt2primary.setdefault(pid, pid)

    return alt2primary, name_of, graph

CHEBI_ANY_RE = re.compile(r"CHEBI:\d+", re.IGNORECASE)

def canonicalize_chebi(chebi_id: str, alt2primary: Dict[str, str]) -> Optional[str]:
    s = (chebi_id or "").strip()
    if not s or s.upper() == "NA":
        return None
    m = CHEBI_ANY_RE.search(s)
    if not m:
        return None
    cid = m.group(0).upper()
    return alt2primary.get(cid, cid)

def load_digest_chebi(digest_path: str, alt2primary: Dict[str, str]) -> Dict[str, List[DigestRow]]:
    idx: Dict[str, List[DigestRow]] = defaultdict(list)
    with open(digest_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"ec_number", "substrate", "chebi_id", "chebi_name"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"{digest_path}: missing columns: {sorted(missing)}")

        for row in r:
            c = canonicalize_chebi(row["chebi_id"], alt2primary)
            if not c:
                continue
            idx[c].append(
                DigestRow(
                    ec_number=row["ec_number"].strip(),
                    substrate=row["substrate"].strip(),
                    chebi_id=c,
                    chebi_name=row["chebi_name"].strip(),
                )
            )
    return idx

def load_nutrient_seeds(best_path: str, alt2primary: Dict[str, str]) -> Dict[str, Dict[str, str]]:
    seeds: Dict[str, Dict[str, str]] = defaultdict(dict)
    with open(best_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        required = {"nutrient_id", "nutrient_name", "chebi_id"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise ValueError(f"{best_path}: missing columns: {sorted(missing)}")

        for row in r:
            nid = row["nutrient_id"].strip()
            nname = row["nutrient_name"].strip()
            c = canonicalize_chebi(row["chebi_id"], alt2primary)
            if not c:
                continue
            seeds[nid][c] = nname
    return seeds

def reconstruct_path(prev: Dict[str, Tuple[str, Edge]], src: str, dst: str) -> List[Edge]:
    if dst == src:
        return []
    edges: List[Edge] = []
    cur = dst
    while cur != src and cur in prev:
        p, e = prev[cur]
        edges.append(e)
        cur = p
    edges.reverse()
    return edges

_CONJ = {"conjugate_base_of", "conjugate_acid_of"}
# Same molecule under a different proton/charge/tautomer state: an enzyme acting on one
# acts on the other, so these are reported as "identity-preserving" rather than as an
# inference step. Enantiomers are NOT in this set - enzymes are routinely stereospecific.
_IDENTITY = _CONJ | {"tautomer_of"}
_SCAFFOLD = {"has_functional_parent", "functional_parent_of",
             "has_parent_hydride", "parent_hydride_of",
             "substituent_group_from", "substituent_group_to"}


def classify_relation(path_edges: List[Edge]) -> str:
    """Coarse label for the path taken, so downstream can filter by inference strength.

    Order matters: the most specific single-type labels are tested first, then the broad
    families, then "mixed". Callers that want only high-confidence links should keep
    exact / conjugate / tautomer and treat scaffold and part_of as weaker evidence.
    """
    if not path_edges:
        return "exact"
    etypes = {e.etype for e in path_edges}
    if etypes == {"is_a"}:
        return "is_a"
    if etypes <= _CONJ:
        return "conjugate" if len(path_edges) > 1 else path_edges[0].etype
    if etypes == {"tautomer_of"}:
        return "tautomer"
    if etypes == {"enantiomer_of"}:
        return "enantiomer"
    if etypes <= _SCAFFOLD:
        return "scaffold"
    if etypes <= {"has_part", "part_of"}:
        return "part_of"
    if etypes <= _IDENTITY:
        return "identity"
    return "mixed"

def dijkstra_limited(graph: Dict[str, List[Edge]],
                     src: str,
                     max_cost: float
                     ) -> Tuple[Dict[str, float], Dict[str, Tuple[str, Edge]]]:
    dist: Dict[str, float] = {src: 0.0}
    prev: Dict[str, Tuple[str, Edge]] = {}
    heap: List[Tuple[float, str]] = [(0.0, src)]

    while heap:
        d, u = heapq.heappop(heap)
        if d != dist.get(u, float("inf")):
            continue
        if d > max_cost:
            continue

        for e in graph.get(u, []):
            nd = d + e.weight
            if nd > max_cost:
                continue
            if nd < dist.get(e.dst, float("inf")):
                dist[e.dst] = nd
                prev[e.dst] = (u, e)
                heapq.heappush(heap, (nd, e.dst))

    return dist, prev

def score_from_cost(cost: float, max_cost: float, relation: str) -> int:
    if max_cost <= 0:
        return 100
    base = max(0.0, 100.0 * (1.0 - cost / max_cost))
    if relation == "exact":
        bonus = 10.0
    elif relation.startswith("conjugate") or relation in ("tautomer", "identity"):
        bonus = 5.0                       # same molecule, different protonation state
    elif relation in ("scaffold", "part_of"):
        bonus = -5.0                      # a structural inference, not an identity
    else:
        bonus = 0.0
    return int(round(max(0.0, min(100.0, base + bonus))))
    
def load_live_nutrients(path: str) -> Set[str]:
    """nutrient_ids that actually carry a measured value in the composition table.

    The nutrient catalog holds ~3,300 entries but only ~1,788 of them carry any value in
    the harmonized table. Mapping an EC to a catalog-only nutrient produces an edge that
    can never reach a food: it was counted as "in the model" while being a dead end.
    Filtering here, at the source of the map, is what makes the downstream `in_model`
    flag mean what it says.

    Accepts either the exported food_nutrients.tsv (any TSV carrying a `nutrient_id`
    column) or a plain newline-delimited id list.
    """
    ids: Set[str] = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        first = fh.readline().rstrip("\n").split("\t")
        col = first.index("nutrient_id") if "nutrient_id" in first else None
        if col is None and first and first[0].strip():
            ids.add(first[0].strip())          # headerless list: keep line 1
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = f[col] if col is not None and col < len(f) else (f[0] if f else "")
            if v.strip():
                ids.add(v.strip())
    if not ids:
        raise SystemExit(f"ERROR: no nutrient_ids read from {path}")
    return ids


def main():
    ap = argparse.ArgumentParser(
        description="Map nutrients to EC numbers via ChEBI graph (is_a + conjugate), alt_id-aware."
    )
    ap.add_argument("--nutrient_best", required=True, help="nutrient_to_chebi.best.tsv")
    ap.add_argument("--digest_chebi", required=True, help="digest_chebi.tsv")
    ap.add_argument("--chebi_obo", required=True, help="chebi.obo")
    ap.add_argument("--out", required=True, help="Output TSV: nutrient_to_ec.tsv")
    ap.add_argument("--max_cost", type=float, default=1.5, help="Max traversal cost (default 1.5)")
    ap.add_argument("--w_is_a", type=float, default=1.0, help="Edge weight for is_a (default 1.0)")
    ap.add_argument("--w_conj", type=float, default=0.5, help="Edge weight for conjugate edges (default 0.5)")

    # --- structural relations (see REL_SYMMETRIC / REL_DIRECTED at the top of the file) ---
    ap.add_argument("--no_structural", action="store_true",
                    help="Traverse only is_a + the two conjugate relations, i.e. reproduce the "
                         "pre-2026 graph exactly. Use to regenerate the legacy map or to measure "
                         "how much the structural edges contribute.")
    ap.add_argument("--w_taut", type=float, default=0.25,
                    help="is_tautomer_of. Same compound, shifted proton: near-free (default 0.25).")
    ap.add_argument("--w_enant", type=float, default=0.5,
                    help="is_enantiomer_of. Same constitution, opposite chirality. NOT free: many "
                         "enzymes are stereospecific, so the hop is priced like a conjugate "
                         "(default 0.5) and is labelled 'enantiomer' in the relation column so it "
                         "can be filtered downstream.")
    ap.add_argument("--w_deriv_up", type=float, default=0.75,
                    help="Derivative -> scaffold (has_functional_parent / has_parent_hydride / "
                         "is_substituent_group_from). Generalizing and one-to-one, so cheap "
                         "(default 0.75, i.e. two hops fit inside the default max_cost).")
    ap.add_argument("--w_deriv_down", type=float, default=2.0,
                    help="Scaffold -> derivative, the reverse of --w_deriv_up. One scaffold has "
                         "many derivatives, so this fans out badly. The default 2.0 exceeds the "
                         "default --max_cost of 1.5 and so switches the direction OFF. That is a "
                         "measured choice, not caution: enabling it (1.25) added exactly ONE "
                         "further nutrient while adding 5,748 rows, i.e. it produced links like "
                         "stigmastane -> alpha-spinasterol that share a skeleton and nothing "
                         "else. Lower it below --max_cost only if you want that recall.")
    ap.add_argument("--w_part", type=float, default=1.0,
                    help="has_part / part_of, e.g. a glycoside to its sugar moiety (default 1.0).")
    ap.add_argument("--keep_mineral_transport", action="store_true",
                    help="Keep EC class 7 (translocase) links for nutrients that ARE a chemical "
                         "element. Off by default. Minerals are, with few exceptions, not "
                         "metabolised: a bacterium either moves the ion across its membrane "
                         "(transport) or uses it as a cofactor that lets some OTHER enzyme act, "
                         "and in neither case is the mineral a dietary substrate being broken "
                         "down. Scoring a food because it is high in sodium and the organism has "
                         "a sodium pump describes ion homeostasis, not digestion. The exceptions "
                         "are real and are KEPT, because they are genuine chemistry on the "
                         "element and are not class 7: sulfur -> sulfur/sulfide oxidoreductases "
                         "(EC 1.12.x, 1.13.11.18) and mercury -> mercuric reductase (EC 1.16.1.1). "
                         "Non-element nutrients keep their transporters too -- an ABC importer "
                         "for arabinose or galacturonic acid IS evidence the sugar is used.")
    ap.add_argument("--max_hub_ec", type=int, default=40,
                    help="Promiscuity guard. A digest node carrying more than this many distinct "
                         "EC numbers is a ubiquitous metabolite (glucose 125, L-glutamate 134, "
                         "dioxygen 617), and inheriting its whole enzyme set through a STRUCTURAL "
                         "inference is how punicalagin acquires 125 enzymes for containing a "
                         "glucose. Such nodes are therefore only accepted when the nutrient IS "
                         "the compound (exact / conjugate / tautomer), never via is_a, scaffold "
                         "or part_of. Default 40 = the 99th percentile of the digest (1.0% of "
                         "nodes). Set 0 to disable.")
    ap.add_argument("--keep_unmatched", action="store_true", help="Output nutrients with no enzyme matches.")
    ap.add_argument("--live_nutrients", default=None,
                    help="Path to food_nutrients.tsv (or a nutrient_id list). Nutrients absent "
                         "from it carry no measured value, so edges to them can never reach a "
                         "food; they are dropped. Strongly recommended — without it the map "
                         "contains dead-end edges that inflate every downstream coverage figure.")
    
    # NEW ARGUMENT: Purge host-absorbed targets
    ap.add_argument("--strict_prebiotic", action=argparse.BooleanOptionalAction, default=True,
                    help="If True, prevents Simple Sugars, Amino Acids, and Water/Energy from mapping to bacterial enzymes.")
    ap.add_argument("--include_simple_sugars", action="store_true",
                    help="Override --strict_prebiotic for simple sugars (glucose, fructose, lactose, sucrose, etc.) — "
                         "needed when modelling colonic bacteria that ferment them or when host-absorption is irrelevant.")
    ap.add_argument("--include_amino_acids", action="store_true",
                    help="Override --strict_prebiotic for free amino acids — needed for proteolytic bacteria (e.g. Clostridium).")
    ap.add_argument("--extra_seeds", action="append", default=[],
                    help="Optional TSV(s) with nutrient_id/nutrient_name/chebi_id rows to append to the seed map. "
                         "Use this to inject novel bacterial substrates (HMOs, mucin, sialic acid, GlcNAc, etc.) "
                         "without modifying the upstream 0/1/2 pipeline. ADDITIVE: a nutrient keeps whatever the "
                         "name matcher already gave it. Repeatable: pass once per file.")
    ap.add_argument("--override_seeds", action="append", default=[],
                    help="Optional TSV(s) in the same format whose rows REPLACE the auto-derived seeds for the "
                         "nutrients they name, instead of adding to them. This is the only way a hand-curated "
                         "mapping can beat the name matcher: --extra_seeds cannot, because a nutrient the matcher "
                         "already resolved keeps that resolution and the curated id merely joins it. Nutrient "
                         "96062 'Alcohol' is the case in point - it matches CHEBI:30879, ChEBI's generic "
                         "'any R-OH' CLASS, and inherits every alcohol-acting enzyme (197 EC under ChEBI 253, up "
                         "from 22) when the food-table column means ethanol. Applied after --extra_seeds, so an "
                         "override wins. Repeatable: pass once per file.")
    args = ap.parse_args()

    live_nutrients = load_live_nutrients(args.live_nutrients) if args.live_nutrients else None
    if live_nutrients is not None:
        print(f"[*] --live_nutrients: {len(live_nutrients):,} nutrients carry a measured value", flush=True)
    n_dropped_dead = 0

    alt2primary, name_of, graph = parse_chebi_obo(
        args.chebi_obo, w_is_a=args.w_is_a, w_conj=args.w_conj,
        w_taut=args.w_taut, w_enant=args.w_enant,
        w_deriv_up=args.w_deriv_up, w_deriv_down=args.w_deriv_down,
        w_part=args.w_part, structural=not args.no_structural)
    n_edges = Counter(e.etype for edges in graph.values() for e in edges)
    print(f"[*] ChEBI graph: {len(graph):,} nodes, {sum(n_edges.values()):,} edges "
          f"({'structural relations ON' if not args.no_structural else 'is_a + conjugate only'})",
          flush=True)
    for k, v in n_edges.most_common():
        print(f"      {k:24} {v:>9,}")
    digest_idx = load_digest_chebi(args.digest_chebi, alt2primary)
    nutrient_seeds = load_nutrient_seeds(args.nutrient_best, alt2primary)
    for seeds_path in (args.extra_seeds or []):
        extra = load_nutrient_seeds(seeds_path, alt2primary)
        n_added = 0
        for nid, smap in extra.items():
            for cid, nname in smap.items():
                if cid not in nutrient_seeds.get(nid, {}):
                    nutrient_seeds[nid][cid] = nname
                    n_added += 1
        print(f"[*] --extra_seeds: merged {n_added} (nutrient, ChEBI) pairs "
              f"({len(extra)} extra nutrients) from {seeds_path}", flush=True)

    # Applied AFTER --extra_seeds so a curated override wins over both the name matcher and
    # any additive seed file. Replaces the whole seed set for a nutrient rather than adding
    # to it, because the failure being corrected is the matcher picking a broad ChEBI CLASS:
    # leaving that class in place alongside the curated compound changes nothing.
    for seeds_path in (args.override_seeds or []):
        override = load_nutrient_seeds(seeds_path, alt2primary)
        n_repl = n_new = 0
        for nid, smap in override.items():
            if nid in nutrient_seeds:
                n_repl += 1
            else:
                n_new += 1
            nutrient_seeds[nid] = dict(smap)
        print(f"[*] --override_seeds: replaced the seeds of {n_repl} nutrients and seeded "
              f"{n_new} new ones from {seeds_path}", flush=True)

    out_fields = [
        "nutrient_id",
        "nutrient_name",
        "seed_chebi_id",
        "seed_chebi_name",
        "matched_chebi_id",
        "matched_chebi_name",
        "relation",
        "cost",
        "steps",
        "score",
        "ec_number",
        "substrate",
        "digest_chebi_name",
        "path_edge_types",
    ]

    wrote_any_for_nutrient: Set[str] = set()
    n_considered = 0                      # nutrients that survived the live/biology filters
    n_hub_blocked = 0
    n_transport_dropped = 0
    hub_ec = {c: len({d.ec_number for d in rows}) for c, rows in digest_idx.items()}
    rel_rows = Counter()                  # rows written, by relation label
    rel_nutrients = defaultdict(set)      # nutrients reached, by relation label

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        seen = set()            

        for nutrient_id, seed_map in nutrient_seeds.items():
            nutrient_name = next(iter(seed_map.values())) if seed_map else ""
            nname_lc = nutrient_name.lower().strip()

            # Drop nutrients with no measured value: an edge to one is unreachable.
            if live_nutrients is not None and nutrient_id not in live_nutrients:
                n_dropped_dead += 1
                continue

            # BIOLOGICAL FILTER: Drop Host-Absorbed and Ubiquitous Noise if flag is active.
            # --include_simple_sugars / --include_amino_acids override each sub-filter
            # independently of --strict_prebiotic so the noise filter stays on.
            if args.strict_prebiotic:
                if _UBIQUITOUS_NOISE.match(nname_lc):
                    continue
                if not args.include_simple_sugars and _SIMPLE_SUGARS.search(nname_lc):
                    continue
                if not args.include_amino_acids and _AMINO_ACIDS.search(nname_lc):
                    continue

            matched_any = False
            n_considered += 1

            for seed_chebi_id in seed_map.keys():
                seed_name = name_of.get(seed_chebi_id, "")
                # ChEBI names elements "<name> atom" and their monatomic ions "<name>(2+)".
                seed_is_element = bool(_ELEMENT_ENTITY_RE.match(seed_name or ""))

                dist, prev = dijkstra_limited(graph, seed_chebi_id, args.max_cost)

                for node, cost in dist.items():
                    if node not in digest_idx:
                        continue

                    path_edges = reconstruct_path(prev, seed_chebi_id, node)
                    relation = classify_relation(path_edges)
                    steps = len(path_edges)

                    # Promiscuity guard: a structural inference may not route through a
                    # ubiquitous metabolite and collect its whole enzyme set.
                    #
                    # EXEMPT: the identity relations (the nutrient IS the compound) and
                    # is_a. Excluding is_a is not a softening - guarding it silently
                    # deleted 437 CORRECT links that the legacy map already shipped,
                    # because an is_a step is usually a specialization to the very same
                    # compound: Maltose -> alpha-maltose, Cellobiose -> beta-cellobiose,
                    # Ribose -> D-ribofuranose, Ornithine -> L-ornithine. Those substrates
                    # have many enzymes because they genuinely have many enzymes. The
                    # guard exists to stop the relations ADDED in 2026 (scaffold, part_of)
                    # from reaching a hub they merely contain, so it applies to those only.
                    if (args.max_hub_ec > 0
                            and relation in ("scaffold", "part_of", "enantiomer", "mixed")
                            and hub_ec.get(node, 0) > args.max_hub_ec):
                        n_hub_blocked += 1
                        continue
                    score = score_from_cost(cost, args.max_cost, relation)

                    path_edge_types = ",".join([e.etype for e in path_edges]) if path_edges else ""
                    matched_name = name_of.get(node, "")

                    for dr in digest_idx[node]:
                        # Translocases move an ion; they do not degrade it. See
                        # --keep_mineral_transport for why this is off for elements only.
                        if (seed_is_element and not args.keep_mineral_transport
                                and dr.ec_number.startswith("7.")):
                            n_transport_dropped += 1
                            continue
                        key = (nutrient_id, seed_chebi_id, node, dr.ec_number)
                        if key in seen:
                            continue
                        seen.add(key)
                        rel_rows[relation] += 1
                        rel_nutrients[relation].add(nutrient_id)
                        w.writerow({
                            "nutrient_id": nutrient_id,
                            "nutrient_name": nutrient_name,
                            "seed_chebi_id": seed_chebi_id,
                            "seed_chebi_name": seed_name,
                            "matched_chebi_id": node,
                            "matched_chebi_name": matched_name,
                            "relation": relation,
                            "cost": f"{cost:.3f}",
                            "steps": str(steps),
                            "score": str(score),
                            "ec_number": dr.ec_number,
                            "substrate": dr.substrate,
                            "digest_chebi_name": dr.chebi_name,
                            "path_edge_types": path_edge_types,
                        })

                    matched_any = True
                    wrote_any_for_nutrient.add(nutrient_id)

            if args.keep_unmatched and (not matched_any):
                w.writerow({
                    "nutrient_id": nutrient_id,
                    "nutrient_name": nutrient_name,
                    "seed_chebi_id": "",
                    "seed_chebi_name": "",
                    "matched_chebi_id": "",
                    "matched_chebi_name": "",
                    "relation": "",
                    "cost": "",
                    "steps": "",
                    "score": "",
                    "ec_number": "",
                    "substrate": "",
                    "digest_chebi_name": "",
                    "path_edge_types": "",
                })

    if live_nutrients is not None:
        print(f"[*] dropped {n_dropped_dead:,} nutrients with no measured value "
              f"(catalog-only entries whose edges could never reach a food)", flush=True)

    # Coverage summary. The point of the structural relations is to shrink the last line;
    # print it every run so a change to the weights is measured, not assumed.
    if not args.keep_mineral_transport:
        print(f"[*] mineral transport filter: dropped {n_transport_dropped:,} translocase "
              f"(EC 7.x) links from element nutrients", flush=True)
    if args.max_hub_ec > 0:
        print(f"[*] promiscuity guard (--max_hub_ec {args.max_hub_ec}): blocked "
              f"{n_hub_blocked:,} structural matches onto ubiquitous metabolites", flush=True)
    n_hit = len(wrote_any_for_nutrient)
    print(f"[*] coverage: {n_hit:,}/{n_considered:,} nutrients reached >=1 EC "
          f"({100.0 * n_hit / n_considered if n_considered else 0:.1f} %); "
          f"{n_considered - n_hit:,} still unmatched", flush=True)
    print("[*] by relation (rows / distinct nutrients):", flush=True)
    for rel, n in rel_rows.most_common():
        excl = len(rel_nutrients[rel] - set().union(
            *(v for k, v in rel_nutrients.items() if k != rel)) ) if len(rel_rows) > 1 else 0
        print(f"      {rel:14} {n:>8,} / {len(rel_nutrients[rel]):>5}"
              f"   ({excl} reachable ONLY this way)")
    print(f"[*] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()