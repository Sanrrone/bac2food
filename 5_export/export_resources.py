#!/usr/bin/env python3
"""export_resources.py — flatten the three bac2food reference resources into TSVs.

  1. enzyme_substrate_chebi.tsv : EC number -> substrate it acts on -> ChEBI id
  2. species_enzymes.tsv        : organism (species + strain) -> EC numbers it carries
  3. food_nutrients.tsv         : food -> nutrient -> amount (per 100 g), + source DB

The script only joins and flattens tables the pipeline already built (the enzyme
substrate digest, the ChEBI mapping, the bact_ec reference, and the FDC-derived food
tables). It recomputes nothing, so the exports are exactly what the predictor sees.

Provenance — the two microbial layers have DIFFERENT sources by design:
  * Bacteria -> EC   : eggNOG v7, built by eggnog/6.1_eggnog7_species_enzymes.py via KEGG KO.
                       This script's own `--only species` route is the LEGACY v6 one
                       (/data/bac2food/bact_ec.tsv, direct EC annotation), retained only so a
                       v6 build stays reproducible; it is no longer what ships.
  * EC -> substrate  : BRENDA (eggnog/2_digest_norm.tsv, from ec_species_substrate.tar.xz).
                       eggNOG publishes no substrate data at any version, so this layer cannot
                       move to v7 — that is correct, not an oversight.
The old DSMZ BRENDA SPARQL scrape supplied species -> enzyme only; it returned incomplete
results, under-reporting every organism, and has been removed (eggNOG replaced that link).

food_nutrients.tsv holds the ANALYTICAL food set: branded label products are excluded by
default (--keep_branded restores them). They were 92% of the values but only ~119 of the
1,788 components, and the predictor has always ignored them.

The outputs are large (species_enzymes.tsv is ~1.8 GB); write them to a volume with
room, e.g.:

    python export_resources.py --out_dir /data/bac2food/exports
"""
from __future__ import annotations
import argparse, csv, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parent.parent
DATA = Path("/data/bac2food")

sys.path.insert(0, str(REPO / "food_DBs"))
import fdc_blocks  # noqa: E402
# The canon name is the predictor's own grouping, imported rather than reimplemented so the
# published column and the rankings can never describe two different groupings.
sys.path.insert(0, str(REPO / "4_predict"))
from bac2food_predict import canonicalize_food_name  # noqa: E402
from _common.non_nutrients import (COPY_RECORD_RE, NON_NUTRIENT_IDS,  # noqa: E402
                                   find_exact_relistings, relisting_candidates,
                                   source_of_bucket_file)


# ------------------------------------------------------------------------------
# The chain filter: only what reaches food -> nutrient -> enzyme -> organism
# ------------------------------------------------------------------------------
# 3_nutrient_to_ec.tsv is the SAME file the predictor loads (bac2food_predict.
# PATH_NUTRIENT_TO_EC), so the released tables and the scorer cannot disagree about
# what "reaches" means. Deriving the set here from a second source would recreate
# exactly the drift the module docstring warns about, where the two readers of the
# bucketed store apply different exclusion policies.
#
# Both directions of the filter are the same statement read from opposite ends:
#   * a nutrient is kept when at least one EC acts on it;
#   * an EC is kept when it acts on at least one nutrient the food table measures.
# A row failing either could never take part in a prediction - the scorer filters to
# the EC-linked nutrient set before it reads a single amount - so pruning them removes
# unreachable records, not usable ones.
#
# What this DOES cost, and it is a real cost: the aggregate composition fields go.
# "Dietary fibre, total", "Energy", "Protein" and "Carbohydrate by difference" are
# mixtures or computed values with no single ChEBI referent, so none of them survives.
# The specific fibres do - inulin, pectin, beta-glucan, cellulose, hemicellulose,
# resistant starch - which is the distinction the resource exists to make. The
# unpruned table remains reproducible with --span_chain keep.
PATH_NUTRIENT_TO_EC = REPO / "0_building" / "3_nutrient_to_ec.tsv"


def chain_linked_ids(path=None) -> tuple[set[int], set[str]]:
    """(nutrient_ids reached by >=1 EC, EC numbers reaching >=1 nutrient)."""
    src = Path(path or PATH_NUTRIENT_TO_EC)
    m = pd.read_csv(src, sep="\t", usecols=["nutrient_id", "ec_number"],
                    dtype={"nutrient_id": "int64", "ec_number": "string"})
    nut = set(m["nutrient_id"].dropna().astype(int).tolist())
    ec = set(m["ec_number"].dropna().astype(str).tolist())
    print(f"[*] chain map {src.name}: {len(nut):,} nutrients reach an EC; "
          f"{len(ec):,} EC reach a nutrient.", flush=True)
    return nut, ec


def assert_fdc_blocks(fdc_ids: pd.Series) -> None:
    """Fail the export if any fdc_id sits outside its source's block.

    Both historical id collisions were silent: a shared id merges two unrelated foods and
    nothing raises, so the damage only surfaced later as nutrient values attached to the
    wrong food. Containment is therefore checked at the point the resource is published,
    not left to whoever writes the next ingest.
    """
    ids = pd.unique(fdc_ids)
    by_source: dict[str, set[int]] = {}
    for i in ids:
        by_source.setdefault(fdc_blocks.source_of(int(i)), set()).add(int(i))
    fdc_blocks.assert_disjoint(by_source)
    print(f"[*] fdc_id blocks OK: {len(ids):,} ids across {len(by_source)} sources, "
          f"all inside their {fdc_blocks.BLOCK_WIDTH:,}-wide block", flush=True)

# Infraspecific ranks belong to the taxon name, not to the strain label; "str."
# is the marker that the strain designation starts. "cf."/"aff." flag an uncertain
# identification and sit between the genus and the epithet.
_INFRA_RANKS = {"subsp.", "subsp", "pv.", "bv.", "var.", "serovar", "biovar", "f.sp."}
_STRAIN_MARK = {"str.", "strain"}
_UNCERTAIN = {"cf.", "aff."}
_TAB_RE = re.compile(r"[\t\r\n]+")
_QUOTE_RE = re.compile(r"^['\"]+|['\"]+$")


def _sanitize(s: pd.Series) -> pd.Series:
    """Collapse embedded tabs/newlines so every record stays on one TSV line."""
    return s.astype("string").fillna("").str.replace(_TAB_RE, " ", regex=True).str.strip()


def _ec_sort_key(ec: pd.Series) -> pd.DataFrame:
    """Numeric sort key so 1.1.1.2 < 1.1.1.10 (lexicographic would invert them)."""
    parts = ec.str.split(".", expand=True)
    return parts.apply(pd.to_numeric, errors="coerce")


def split_organism(name: str) -> tuple[str, str, str]:
    """Split a UniProt-style organism name into (genus, species, strain).

    'Acetobacter aceti NBRC 3283'                  -> ('Acetobacter', 'Acetobacter aceti', 'NBRC 3283')
    'Propionibacterium freudenreichii subsp. shermanii JS'
        -> ('Propionibacterium', 'Propionibacterium freudenreichii subsp. shermanii', 'JS')
    'Bacteriovorax sp. BAL6_X'                     -> ('Bacteriovorax', 'Bacteriovorax sp.', 'BAL6_X')
    'Francisella cf. tularensis subsp. novicida 3523'
        -> ('Francisella', 'Francisella cf. tularensis subsp. novicida', '3523')
    """
    toks = str(name).split()
    if not toks:
        return "", "", ""
    # Quotes are stripped from genus/species only ("'Nostoc azollae' 0708"); a quoted
    # strain designation ("str. 'Deep ecotype'") is kept verbatim.
    genus = _QUOTE_RE.sub("", toks[0])
    if len(toks) == 1:
        return genus, genus, ""
    sp_toks, i = [genus], 1
    if toks[i].lower() in _UNCERTAIN:      # 'Francisella cf. tularensis'
        sp_toks.append(toks[i])
        i += 1
    if i < len(toks):                      # the specific epithet (or 'sp.')
        sp_toks.append(toks[i])
        i += 1
    # absorb infraspecific ranks (subsp. shermanii, pv. tomato, ...) into the taxon name
    while i + 1 < len(toks) and toks[i].lower() in _INFRA_RANKS:
        sp_toks += [toks[i], toks[i + 1]]
        i += 2
    if i < len(toks) and toks[i].lower() in _STRAIN_MARK:
        i += 1
    species = _QUOTE_RE.sub("", " ".join(sp_toks))
    return genus, species, " ".join(toks[i:])


# ==============================================================================
# 1. enzyme -> substrate -> ChEBI
# ==============================================================================
def export_enzyme_substrate(args) -> None:
    out = Path(args.out_dir) / "enzyme_substrate_chebi.tsv"
    print(f"[*] Building {out.name} ...", flush=True)

    d = pd.read_csv(REPO / "chebi/digest_to_chebi.tsv", sep="\t", dtype="string")
    norm = pd.read_csv(REPO / "eggnog/2_digest_norm.tsv", sep="\t", dtype="string")
    norm = norm.rename(columns={"ec": "ec_number"})[["ec_number", "substrate", "substrate_normalized"]]
    d = d.merge(norm.drop_duplicates(["ec_number", "substrate"]), on=["ec_number", "substrate"], how="left")

    # Which (EC, substrate) pairs actually reach an FDC nutrient in the model, and via which one.
    n = pd.read_csv(REPO / "0_building/3_nutrient_to_ec.tsv", sep="\t", dtype="string")
    n["score_num"] = pd.to_numeric(n["score"], errors="coerce")
    g = n.sort_values("score_num", ascending=False).groupby(["ec_number", "substrate"], sort=False)
    link = pd.DataFrame({
        "nutrient_ids":   g["nutrient_id"].agg(lambda s: ",".join(dict.fromkeys(s))),
        "nutrient_names": g["nutrient_name"].agg(lambda s: "; ".join(dict.fromkeys(s))),
        "model_relation": g["relation"].first(),
        "model_score":    g["score"].first(),
    }).reset_index()

    d = d.merge(link, on=["ec_number", "substrate"], how="left")
    d["in_model"] = d["nutrient_ids"].notna().map({True: "yes", False: "no"})
    d = d.rename(columns={"match_type": "chebi_match_type"})

    for c in ("substrate", "substrate_normalized", "chebi_name", "nutrient_names"):
        d[c] = _sanitize(d[c])

    cols = ["ec_number", "substrate", "substrate_normalized", "chebi_id", "chebi_name",
            "chebi_match_type", "in_model", "nutrient_ids", "nutrient_names",
            "model_relation", "model_score"]
    key = _ec_sort_key(d["ec_number"])
    d = d.assign(**{f"_k{i}": key[i] for i in key.columns}).sort_values(
        [f"_k{i}" for i in key.columns] + ["substrate"], kind="mergesort")

    d[cols].to_csv(out, sep="\t", index=False, na_rep="")
    print(f"[*] {out}: {len(d):,} rows | {d.ec_number.nunique():,} EC numbers | "
          f"{d.chebi_id.notna().sum():,} rows with a ChEBI id | "
          f"{(d.in_model == 'yes').sum():,} rows used by the model", flush=True)


# ==============================================================================
# 2. species / strain -> enzymes
# ==============================================================================
_BACT_EC_COLS = ["tax_id", "organism", "kingdom", "ec_number"]


def read_bact_ec(path) -> pd.DataFrame:
    """Load the bacterial EC reference as unique (tax_id, organism, EC) rows.

    bact_ec.tsv is headerless and carries one row per annotated protein (58M rows,
    3.2 GB), so the same (organism, EC) fact repeats thousands of times; the distinct
    set is only ~9.6M rows. Read it in chunks, encode the three low-cardinality text
    columns as integer codes, and deduplicate as we go, so peak memory stays flat.
    A .parquet holding the same table is accepted too (already deduplicated).
    """
    path = Path(path)
    if path.suffix == ".parquet":
        u = pd.read_parquet(path).rename(columns={"species": "organism"})
        return u[["tax_id", "organism", "ec_number"]].drop_duplicates()

    cols = ["tax_id", "organism", "ec_number"]          # column 3 ("Bacteria") is constant
    maps: dict[str, dict] = {c: {} for c in cols}
    parts = []
    reader = pd.read_csv(path, sep="\t", header=None, names=_BACT_EC_COLS, usecols=cols,
                         dtype={c: "category" for c in cols}, chunksize=8_000_000)
    for ch in reader:
        ch = ch.dropna()
        for c in cols:
            m = maps[c]
            codes = np.array([m.setdefault(v, len(m)) for v in ch[c].cat.categories],
                             dtype="int32")
            ch[c] = codes[ch[c].cat.codes.to_numpy()]
        parts.append(ch.drop_duplicates())
        print(f"    read {sum(len(p) for p in parts):,} unique rows so far ...", end="\r", flush=True)

    df = pd.concat(parts, ignore_index=True).drop_duplicates()
    for c in cols:                                       # integer codes -> back to text
        inv = np.empty(len(maps[c]), dtype=object)
        for v, i in maps[c].items():
            inv[i] = v
        df[c] = inv[df[c].to_numpy()]
    print(f"    bact_ec: {len(df):,} unique (tax_id, organism, EC) rows" + " " * 20, flush=True)
    return df


def export_species_enzymes(args) -> None:
    out = Path(args.out_dir) / "species_enzymes.tsv"
    print(f"[*] Building {out.name} ...", flush=True)

    # Which organism->EC layer is the source of record.
    #
    # There are two, they carry the same filename, and picking the wrong one silently
    # DOWNGRADES the resource. The published layer is the eggNOG v7 build (4,819 EC,
    # routed through KEGG KO by eggnog/6.1_eggnog7_species_enzymes.py). The legacy v6
    # route parsed a direct EC annotation into bact_ec.tsv and reaches only 4,291 EC,
    # covering 75.9% of nutrient-reaching EC against v7's 82.1%.
    #
    # This function used to read the v6 file unconditionally, so `--only species` quietly
    # replaced a v7 export with a v6 one - the hazard the module docstring names but did
    # not defend against. It now reads the v7 layer by default and the v6 route has to be
    # asked for. The v7 PIPELINE cannot currently be re-run (e7.taxid_info.tsv.gz is not
    # on disk), which is exactly why its output is kept as a file of record rather than
    # regenerated on demand.
    if args.species_from_bact_ec:
        print(f"[!] Reading the LEGACY v6 layer {args.bact_ec} - 4,291 EC against v7's "
              f"4,819. This is not what the deposit ships.", flush=True)
        df = read_bact_ec(args.bact_ec)
    else:
        src = Path(args.species_source)
        if not src.exists():
            raise SystemExit(f"ERROR: v7 organism->EC layer not found at {src}. Pass "
                             f"--species_source, or --species_from_bact_ec to accept the "
                             f"v6 downgrade deliberately.")
        df = pd.read_csv(src, sep="\t", dtype="string")
        print(f"[*] v7 organism->EC layer: {len(df):,} rows, "
              f"{df['ec_number'].nunique():,} EC from {src}", flush=True)
    df["organism"] = _sanitize(df["organism"])

    # The v7 layer already carries genus/species/strain - it was written with them - so
    # re-deriving them here would collide on the merge (genus_x/genus_y) and, worse,
    # would overwrite the pipeline's own taxonomy with a re-parse of the display name.
    # Only the v6 route arrives as bare (organism, EC) and needs the split.
    if not {"genus", "species", "strain"} <= set(df.columns):
        # Parse each distinct organism once (~3k), then broadcast back over the 9.6M rows.
        orgs = pd.DataFrame({"organism": df["organism"].unique()})
        parsed = orgs["organism"].map(split_organism)
        orgs["genus"], orgs["species"], orgs["strain"] = (parsed.str[i] for i in range(3))
        df = df.merge(orgs, on="organism", how="left")
    for c in ("genus", "species", "strain"):
        df[c] = _sanitize(df[c])

    if args.span_chain == "require":
        _, keep_ec = chain_linked_ids()
        n0, ec0 = len(df), df["ec_number"].nunique()
        df = df[df["ec_number"].isin(keep_ec)]
        print(f"[*] Excluded {n0 - len(df):,} rows carrying an EC that reaches no nutrient "
              f"({ec0:,} -> {df['ec_number'].nunique():,} EC); an organism-to-enzyme edge "
              f"with no enzyme-to-food edge cannot take part in a prediction. "
              f"Pass --span_chain keep to retain them.", flush=True)

    key = _ec_sort_key(df["ec_number"])
    df = df.assign(**{f"_k{i}": key[i] for i in key.columns}).sort_values(
        ["organism"] + [f"_k{i}" for i in key.columns], kind="mergesort")

    cols = ["tax_id", "genus", "species", "strain", "organism", "ec_number"]
    df[cols].to_csv(out, sep="\t", index=False, na_rep="")
    print(f"[*] {out}: {len(df):,} rows | {df.organism.nunique():,} organisms "
          f"({df.species.nunique():,} species, {(df.strain != '').sum():,} strain-resolved rows) | "
          f"{df.ec_number.nunique():,} EC numbers", flush=True)


# ==============================================================================
# 3. food -> nutrient -> amount
# ==============================================================================
def read_food_nutrient_with_source(bucketed_dir) -> tuple[pd.DataFrame, list[str]]:
    """Read the bucketed food_nutrient store, tagging every row with its source DB.

    Returns (frame with an int source_code column, list mapping code -> label).

    A (food, nutrient) measured by several DBs appears once per DB. Where they report the
    SAME amount it is one fact with several witnesses, so the rows collapse into one and the
    label becomes "ciqual;fdc". Where they report DIFFERENT amounts they are kept as separate
    rows: the disagreement is real data, and merging it would invent a value neither DB gives.
    """
    files = sorted(Path(bucketed_dir).glob("*/*.parquet"))
    if not files:
        raise SystemExit(f"No parquet files under {bucketed_dir}")
    labels = sorted({source_of_bucket_file(f) for f in files})
    code = {s: i for i, s in enumerate(labels)}

    frames = []
    for f in files:
        t = pq.read_table(f, columns=["fdc_id", "nutrient_id", "amount"]).to_pandas()
        t["source_code"] = np.int16(code[source_of_bucket_file(f)])
        frames.append(t)
    df = pd.concat(frames, ignore_index=True)
    del frames

    dup = df.duplicated(["fdc_id", "nutrient_id"], keep=False)
    uniq, d = df[~dup], df[dup]
    del df

    merged = (d.groupby(["fdc_id", "nutrient_id", "amount"], dropna=False, sort=False)["source_code"]
              .agg(lambda s: ";".join(sorted({labels[c] for c in s})))
              .reset_index())
    shared = int((merged["source_code"].str.contains(";")).sum())

    # Give every distinct combined label ("ciqual;fdc") its own code, so the whole column
    # stays an int and only ~a dozen label strings ever exist.
    for lab in merged["source_code"].unique():
        if lab not in code:
            code[lab] = len(labels)
            labels.append(lab)
    merged["source_code"] = merged["source_code"].map(code).astype("int16")

    out = pd.concat([uniq, merged[["fdc_id", "nutrient_id", "amount", "source_code"]]],
                    ignore_index=True)
    print(f"[*] Sources: {', '.join(l for l in labels if ';' not in l)}", flush=True)
    print(f"[*] Collapsed {len(d) - len(merged):,} duplicate rows; {shared:,} values are "
          f"reported identically by >1 DB (joined with ';'). Rows where DBs disagree on the "
          f"amount are kept separate.", flush=True)
    return out, labels


def export_food_nutrients(args) -> None:
    out = Path(args.out_dir) / "food_nutrients.tsv"
    print(f"[*] Building {out.name} ...", flush=True)

    food = pd.read_parquet(args.food, columns=["fdc_id", "data_type", "description",
                                               "food_category_id", "food_category",
                                               "source_food_code"])
    # fdc_id is an opaque accession, so the identifier a user needs in order to look a food up
    # in the table it came from has to travel with it. Under the old offset scheme this was
    # implicit -- recoverable only by knowing the block base -- which meant it was effectively
    # not published at all.
    assert_fdc_blocks(food["fdc_id"])
    # The category lives in three different places depending on where the food came from:
    # the merged non-USDA foods already carry a food_category label; FDC curated foods carry a
    # numeric food_category_id that has to be looked up; branded foods carry their own
    # free-text category in that same id column. Resolve in that order.
    cat = pd.read_csv(args.food_category, dtype="string").set_index("id")["description"]
    fcid = food["food_category_id"].astype("string")
    is_num = fcid.str.fullmatch(r"\d+", na=False)
    food["food_category"] = (food["food_category"].astype("string")
                             .fillna(fcid.where(is_num).map(cat))
                             .fillna(fcid.where(~is_num)))
    # `canon` folds preparation FORM and spelling variants of one food onto a single name
    # ("Carrots, sliced, frozen, unprepared" and "carrot, raw" both -> "carrot"), so a user
    # can group rows the way the predictor does without reimplementing the rules. It does
    # NOT fold preparation STATE: drying, juicing, frying and sweetening change per-100 g
    # composition, so "Carrot, dried" keeps a canon of its own instead of handing the
    # carrot group a dehydrated food's values. It is a grouping key, NOT an identity
    # claim: it is derived from `description` alone and never crosses sources
    # deliberately. 116,053 foods fold to 22,094 canon names
    # (the 8 blank-description foods carry no canon and are not among them).
    # Computed on the distinct descriptions, not per row: 2.1M rows, 52k distinct names.
    _canon_of = {d: canonicalize_food_name(d) for d in food["description"].dropna().unique()}
    food["canon"] = _sanitize(food["description"].map(_canon_of).astype("string"))
    food["description"] = _sanitize(food["description"])
    food["food_category"] = _sanitize(food["food_category"])
    food["source_food_code"] = _sanitize(food["source_food_code"].astype("string"))
    food = food[["fdc_id", "source_food_code", "description", "canon", "data_type",
                 "food_category"]]

    # Per-record licence and provenance. The 16 sources fall under incompatible terms, so
    # the table cannot carry one blanket licence: rights travel with the value. The tier
    # table is the single source of truth and is shipped as Supplementary Table S1 — do not
    # duplicate its contents here, or the paper and the data will drift apart.
    lic = load_licence_table(args.licence_table)

    # Branded label products are excluded by default. They dominate the count (1.89M of
    # 2.01M foods, 92% of all values) but carry only the ~119 components a nutrition-facts
    # panel declares, so they add scale without analytical depth and would misrepresent
    # what the resource actually measures. This matches the predictor, which has shipped
    # with `drop_branded: true` in 4_predict/parameters.yaml from the start — so the
    # default export is now exactly the food set the software scores against.
    n_before = len(food)
    if not args.keep_branded:
        keep = (food["data_type"].astype("string").fillna("") != "branded_food")
        food = food[keep.to_numpy(dtype=bool)]
        print(f"[*] Excluded {n_before - len(food):,} branded_food entries "
              f"({len(food):,} foods retained); pass --keep_branded to include them.",
              flush=True)

    # FNDDS is excluded on the same grounds, and the evidence for it is stronger than for
    # branded. Not one of its 353,015 values carries a `data_points` sample count and not one
    # carries a `derivation_id`, because FNDDS values are recipe-modelled from SR Legacy
    # ingredients for dietary-survey coding rather than assayed. A branded label is at least an
    # assertion about a specific product; an FNDDS value is computed from rows already in this
    # table, so keeping it would partly double-count. Its nutrient set is 65 wide against
    # Foundation's 227. Removing it also ends the FDC backbone's majority: 51.2% -> 42.2% of
    # all values. No nutrient and no enzyme link is lost — all 1,749 nutrients and all 598
    # enzyme-linked ones survive, because every one of them is reported by some other source.
    n_before = len(food)
    if not args.keep_modelled:
        keep = (food["data_type"].astype("string").fillna("") != "survey_fndds_food")
        food = food[keep.to_numpy(dtype=bool)]
        print(f"[*] Excluded {n_before - len(food):,} survey_fndds_food entries "
              f"({len(food):,} foods retained); modelled from recipes, not measured. "
              f"Pass --keep_modelled to include them.", flush=True)

    n_before = len(food)
    keep = ~food["description"].astype("string").fillna("").str.contains(COPY_RECORD_RE,
                                                                        regex=True, na=False)
    food = food[keep.to_numpy(dtype=bool)]
    if n_before != len(food):
        print(f"[*] Excluded {n_before - len(food):,} source editing-layer copies "
              f"({len(food):,} foods retained); a copied recipe keeps the original's name "
              f"plus a timestamp and is then edited away from it.", flush=True)

    nut = pd.read_csv(args.nutrient, usecols=["id", "name", "unit_name"])
    nut = nut.rename(columns={"id": "nutrient_id", "name": "nutrient_name"})
    nut["nutrient_name"] = _sanitize(nut["nutrient_name"])

    df, labels = read_food_nutrient_with_source(args.bucketed_dir)
    write_licence_table(lic, Path(out).parent / "licences.tsv", args.restricted, labels)

    tbl = pa.Table.from_pandas(df, preserve_index=False).sort_by(
        [("fdc_id", "ascending"), ("nutrient_id", "ascending")])
    del df
    n_before = tbl.num_rows
    tbl = tbl.filter(pc.invert(pc.is_in(tbl["nutrient_id"],
                                        value_set=pa.array(sorted(NON_NUTRIENT_IDS)))))
    if n_before != tbl.num_rows:
        print(f"[*] Excluded {n_before - tbl.num_rows:,} rows carrying one of "
              f"{len(NON_NUTRIENT_IDS)} source columns that are not composition values "
              f"(identifiers, version stamps, conversion factors, as-purchased yields); "
              f"see food_DBs/_common/non_nutrients.py.", flush=True)

    # Exact re-listings: same source, same name, same complete value vector. The rule
    # and its rationale live in food_DBs/_common/non_nutrients.py because the predictor
    # has to reach the same verdict - it scores against this same store, so a food this
    # table does not contain must not be scorable either.
    _cand = relisting_candidates(food)
    _sub = tbl.filter(pc.is_in(tbl["fdc_id"], value_set=pa.array(sorted(_cand)))) \
              .select(["fdc_id", "nutrient_id", "amount", "source_code"]).to_pandas()
    _redundant = find_exact_relistings(food, _sub)
    del _sub, _cand
    if _redundant:
        food = food[~food["fdc_id"].isin(_redundant)]
        print(f"[*] Excluded {len(_redundant):,} exact re-listings ({len(food):,} foods "
              f"retained); same source, same name, same complete value vector.", flush=True)
        tbl = tbl.filter(pc.invert(pc.is_in(tbl["fdc_id"],
                                            value_set=pa.array(sorted(_redundant)))))

    if args.span_chain == "require":
        keep_nut, _ = chain_linked_ids()
        n0 = tbl.num_rows
        tbl = tbl.filter(pc.is_in(tbl["nutrient_id"],
                                  value_set=pa.array(sorted(keep_nut), type=tbl["nutrient_id"].type)))
        print(f"[*] Excluded {n0 - tbl.num_rows:,} rows whose nutrient reaches no EC "
              f"(aggregate and computed fields: total dietary fibre, energy, protein, "
              f"carbohydrate by difference). The specific fibres are retained.", flush=True)
        # A food whose every value has just been dropped is no longer a food this
        # resource says anything about, so it leaves the food table too - otherwise
        # the header count and the row count describe different sets.
        live = set(tbl["fdc_id"].to_pylist())
        n_food = len(food)
        food = food[food["fdc_id"].isin(live)]
        print(f"[*] Excluded {n_food - len(food):,} foods left with no value that reaches "
              f"an EC ({len(food):,} foods retained).", flush=True)

    label_arr = np.array(labels, dtype=object)
    total = tbl.num_rows
    print(f"[*] {total:,} food_nutrient rows; joining metadata and writing ...", flush=True)

    written = orphan = 0
    dropped_restricted = [0]
    with open(out, "w", newline="") as fh:
        w = None
        for start in range(0, total, args.chunk_rows):
            ch = tbl.slice(start, args.chunk_rows).to_pandas()
            ch["source_db"] = label_arr[ch.pop("source_code").to_numpy()]
            ch = ch.merge(food, on="fdc_id", how="inner").merge(nut, on="nutrient_id", how="left")
            orphan += len(tbl.slice(start, args.chunk_rows)) - len(ch)
            # The Restricted tier is empty as of 2026-08-18: RIVM's permission moved NEVO to
            # `Provider permission`, and it was the only member. The mechanism stays because
            # the next source with terms like NEVO's old ones needs it, and because a tier
            # named in the licence table is how that decision gets recorded rather than coded.
            if args.restricted == "exclude":
                keep = ~ch["source_db"].isin(lic["restricted"])
                dropped_restricted[0] += int((~keep).sum())
                ch = ch[keep]
            ch = ch[["fdc_id", "source_food_code", "description", "canon", "data_type",
                     "food_category", "nutrient_id", "nutrient_name", "unit_name", "amount",
                     "source_db"]]
            ch.to_csv(fh, sep="\t", index=False, header=(w is None), na_rep="",
                      quoting=csv.QUOTE_MINIMAL)
            w = True
            written += len(ch)
            print(f"    {written:,}/{total:,} rows", end="\r", flush=True)

    # `orphan` counts every row the inner join discarded, which is now two different
    # things: foods carrying no entry in food.parquet at all, and branded foods removed
    # above. Keep them separable or the next reader will read one as the other.
    reason = "fdc_id absent from food table" if args.keep_branded else \
             "fdc_id absent from food table, or branded"
    print(f"\n[*] {out}: {written:,} rows | {orphan:,} rows dropped ({reason})", flush=True)
    if args.restricted == "exclude":
        print(f"[*] Excluded {dropped_restricted[0]:,} rows from restricted sources "
              f"({', '.join(sorted(lic['restricted'])) or 'none'}); their terms forbid "
              f"redistribution of amended values. Pass --restricted include to build the "
              f"LOCAL, non-redistributable copy.", flush=True)
    else:
        print(f"[!] --restricted include: this file contains restricted-source values and "
              f"MUST NOT be deposited or shared. Local use only.", flush=True)
    print("[*] amount is per 100 g edible portion, in unit_name units, except\n    where the canon ends in '0% moisture basis': those are FDC analytical\n    rows reported on a dry matter basis.", flush=True)


# Most restrictive tier wins when sources agree on a value and are `;`-joined: a composite
# row inherits the strictest terms of its constituents, never the loosest.
# "Provider permission" ranks above the public licences and below "Restricted": those
# rows ship, but on a permission granted to this deposit rather than on terms a reader
# inherits, so a composite touching one must resolve to it rather than to CC BY.
_TIER_RANK = {"Open": 0, "Copyleft": 1, "NonCommercial + ShareAlike": 2,
              "Provider permission": 3, "Restricted": 4}


def load_licence_table(path: str) -> dict:
    """source_db -> (version, licence, redistribution_flag, attribution). Composite
    `a;b` labels are resolved on demand to the strictest tier of their parts."""
    import csv as _csv
    rows = {}
    with open(path, encoding="utf-8", newline="") as fh:
        for r in _csv.DictReader(fh):
            rows[r["source_db"]] = r
    if not rows:
        raise SystemExit(f"[!] {path} holds no licence rows")
    def resolve(label):
        parts = [p for p in label.split(";") if p in rows]
        if not parts:
            return ("", "unknown", "no", "")
        worst = max(parts, key=lambda p: _TIER_RANK.get(rows[p]["tier"], max(_TIER_RANK.values())))
        r = rows[worst]
        flag = "yes" if r["derived_values_redistributable"].startswith("yes") else (
               "no" if r["derived_values_redistributable"].startswith("no") else "conditional")
        return ("; ".join(rows[p]["version"] for p in parts),
                r["licence"],
                flag,
                " ".join(rows[p]["attribution_string"] for p in parts))
    return {"rows": rows, "resolve": resolve,
            "restricted": {k for k, r in rows.items() if r["tier"] == "Restricted"}}


def write_licence_table(lic: dict, out, restricted_policy: str, observed=()) -> None:
    """Write `licences.tsv`: one row per source_db, joined to food_nutrients on `source_db`.

    The rights belong here rather than on every value row. Repeating a licence and an
    attribution string across 2.26M rows expressed sixteen facts in 265 MB, and a normalized
    table is the form a user can actually read, cite and diff. `source_db` is already the join
    key, so nothing is lost: filtering the composition table by licence is a single join.
    """
    import csv as _csv
    cols = ["source_db", "database", "provider_country", "licence", "version", "tier",
            "derived_values_redistributable", "attribution_string", "in_deposit"]
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        # Every label that can appear in source_db gets a row, composites included: a value
        # reported identically by two sources carries the `;`-joined label, and a user joining
        # naively on source_db must not hit a missing key. Composites inherit the strictest
        # terms of their parts.
        # A restricted source contributes no rows to an `exclude` build, so it must not appear
        # in the published rights table either — a row for a source nobody can join to is just
        # a puzzle for the reader. It stays in the INTERNAL tier table, which is how the build
        # knows to drop it; dropping it from there instead silently readmits its values.
        labels = sorted(set(lic["rows"]) | {l for l in observed if l})
        if restricted_policy == "exclude":
            labels = [l for l in labels
                      if not any(p in lic["restricted"] for p in l.split(";"))]
        for src in labels:
            r = lic["rows"].get(src)
            if r is None:
                ver, licence, flag, attr = lic["resolve"](src)
                parts = [p for p in src.split(";") if p in lic["rows"]]
                r = {"database": " + ".join(lic["rows"][p]["database"] for p in parts),
                     "provider_country": " + ".join(lic["rows"][p]["provider_country"] for p in parts),
                     "licence": licence, "version": ver, "tier": "composite (strictest of parts)",
                     "derived_values_redistributable": flag, "attribution_string": attr}
            row = {c: r.get(c, "") for c in cols}
            row["source_db"] = src
            # Says plainly whether this source's DERIVED values are in the deposited table,
            # so a reader never has to infer it from the licence prose.
            row["in_deposit"] = ("no" if (restricted_policy == "exclude"
                                          and any(p in lic["restricted"]
                                                  for p in src.split(";"))) else "yes")
            w.writerow(row)
    print(f"[*] Wrote {out} ({len(lic['rows'])} sources; rights join on source_db)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # REQUIRED, deliberately no default. The output filenames do not encode which input
    # version produced them (a v6 build and a v7 build are both "species_enzymes.tsv"),
    # so a default of <DATA>/exports meant that regenerating the legacy layer silently
    # overwrote the shipped v7 deposit. Make the target an explicit choice every time.
    ap.add_argument("--out_dir", required=True,
                    help="where the TSVs are written (REQUIRED). Pass "
                         f"{DATA / 'exports'} to write the deposit itself; use a scratch "
                         "directory for rebuilds you do not intend to ship.")
    ap.add_argument("--only", choices=["enzymes", "species", "foods"], action="append",
                    help="export only this resource (repeatable; default: all three)")
    ap.add_argument("--bact_ec", default=str(DATA / "bact_ec.tsv"),
                    help="bacterial EC reference; .tsv or .parquet (default: %(default)s)")
    ap.add_argument("--food", default=str(DATA / "food.parquet"))
    ap.add_argument("--food_category", default=str(DATA / "food_category.csv"))
    ap.add_argument("--nutrient", default=str(DATA / "nutrient.csv"))
    ap.add_argument("--bucketed_dir", default=str(DATA / "food_nutrient_bucketed"))
    ap.add_argument("--licence_table",
                    default=str(Path(__file__).resolve().parent / "licence_tiers.csv"),
                    help="per-source licence tiers; the single source of truth for the four "
                         "rights columns and for which sources are restricted")
    ap.add_argument("--restricted", choices=["exclude", "include"], default="exclude",
                    help="'exclude' (default) drops sources whose licence forbids "
                         "redistributing derived values — this is the DEPOSIT build. "
                         "'include' builds the full local table, which must not be shared.")
    ap.add_argument("--keep_branded", action="store_true",
                    help="keep FDC branded_food entries in food_nutrients.tsv. They are "
                         "excluded by default: nutrition-facts label declarations, ~119 "
                         "components, no analytical depth (see export_food_nutrients)")
    ap.add_argument("--keep_modelled", action="store_true",
                    help="keep FDC survey_fndds_food entries. Excluded by default: recipe-"
                         "modelled for dietary-survey coding, 0%% carry a sample count or a "
                         "derivation record, 65 components (see export_food_nutrients)")
    ap.add_argument("--species_source", default=str(DATA / "species_enzymes_v7_full.tsv"),
                    help="the unpruned eggNOG v7 organism->EC layer, the source of record "
                         "for species_enzymes.tsv")
    ap.add_argument("--species_from_bact_ec", action="store_true",
                    help="build the organism->EC layer from the LEGACY v6 bact_ec.tsv "
                         "instead (4,291 EC vs v7's 4,819). Not what the deposit ships.")
    ap.add_argument("--span_chain", choices=["require", "keep"], default="require",
                    help="require (default): publish only records that span "
                         "food->nutrient->enzyme->organism. Nutrients no EC acts on, foods "
                         "left with no such value, and EC that reach no nutrient are "
                         "excluded - none of them can take part in a prediction. "
                         "keep: the unfiltered resource, including the aggregate "
                         "composition fields (energy, protein, total dietary fibre).")
    ap.add_argument("--chunk_rows", type=int, default=2_000_000,
                    help="rows per write chunk for the food table (default: %(default)s)")
    args = ap.parse_args()

    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    wanted = set(args.only or ["enzymes", "species", "foods"])
    if "enzymes" in wanted:
        export_enzyme_substrate(args)
    if "species" in wanted:
        export_species_enzymes(args)
    if "foods" in wanted:
        export_food_nutrients(args)
    print("[*] Done.", flush=True)


if __name__ == "__main__":
    sys.exit(main())
