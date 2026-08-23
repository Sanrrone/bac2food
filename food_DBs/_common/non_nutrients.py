"""Shared data policy for the bucketed food_nutrient store.

Both readers of that store apply it: 5_export/export_resources.py, which
publishes the table, and 4_predict/bac2food_predict.py, which scores against
it. A food or a column excluded by one and kept by the other means the
predictor recommends something the released table does not contain, so the
rules live here rather than in either caller.

SOURCE COLUMNS THAT ARE NOT COMPOSITION VALUES

Every ingester mints a new nutrient_id for any source column it cannot map onto
FDC's catalogue. That is right for a real measurement the catalogue lacks
(dextrin, erythritol, C4-C22 fatty acids all arrived that way) and wrong for the
metadata columns that sit in the same sheets. Ten of them reached the export and
were published as nutrients, 26,350 rows in all, including a database version
stamp carried on every one of BioFoodComp's 10,133 foods.

None of the ten reaches an EC number, so no prediction ever used one; the damage
is to the table's meaning, not to the biology.

Two kinds are listed here. The first are not measurements at all: identifiers,
version stamps and the conversion factors a compiler needs to derive protein from
nitrogen. The second are real measurements of the wrong thing: `Edible proportion`
and `Refuse` describe the whole food as purchased, and `Dry matter` restates
`Water`, so none of them is an amount per 100 g of edible portion, which is what a
row of file 1 is defined to be.

Keyed by nutrient_id because that is what survives into the export; the names are
what the sources call them.
"""
from __future__ import annotations

NON_NUTRIENT_IDS: dict[int, str] = {
    # identifiers, version stamps, conversion factors
    96926:  "Latest Revision in Version (biofoodcomp database version stamp)",
    220001: "Code (biofoodcomp, 'as used in version 1.0 and 1.1')",
    220002: "XN (biofoodcomp nitrogen -> protein conversion factor)",
    96894:  "ISSCAAP (biofoodcomp FAO fish species classification code)",
    260202: "Nitrogen conversion factor (mccance)",
    260193: "Glycerol conversion factor (mccance)",
    220051: "XFA (biofoodcomp, labelled 'internal use' by the source)",
    # measurements, but not per 100 g of edible portion
    260188: "Edible proportion (mccance, yield of the food as purchased)",
    250002: "Refuse (stfcj, inedible fraction of the food as purchased)",
    96421:  "Dry matter (biofoodcomp, the complement of Water)",
}

# Livsmedelsverket ships its editing layer alongside the reference set. Copying a
# recipe in their system produces a row named after the original plus a timestamp,
# and the copy is then edited away from it: 24002560 is filed as mashed potato but
# carries 6.1 g water against the original's 78.8, i.e. it was rebased onto the dry
# mix. The name describes neither food.
COPY_RECORD_RE = r"^Kopia av |\(Copy of "


# ---------------------------------------------------------------------------
# Exact re-listings.
#
# A compilation can carry the same food more than once. BioFoodComp repeats whole
# cultivar blocks across its per-nutrient sheets, so codes 301975-301984 (ten
# chickpea cultivars) reappear verbatim at 301985-301994. Where the repeat also
# reports the same values it is the same measurement written twice; where it
# reports DIFFERENT values it is a different measurement and is kept, which is why
# the test is the complete value vector rather than the name.
#
# Two callers need the same answer or they disagree about what exists: the export
# writes the table, and the predictor scores against the store the table is built
# from. It matters more than the count suggests. The 74 re-listings carry 591 rows,
# 0.03% of the corpus, but they are not spread evenly - they land on 35 canonical
# foods, and for four of them (crispbread, barley bran, cracker whole wheat-based,
# bean light red kidney) the duplicate IS half the group, so model_count doubles.
# The per-nutrient MEAN is unaffected either way, since averaging a value with
# itself returns it.
# ---------------------------------------------------------------------------

# Distinct physical samples of one product, not one record written many times: 22
# of the "MILK, 2%" subsamples each report a single value of 8.0, and agreement on
# one rounded number says nothing at all. 15,610 of the 15,817 matches inside this
# data type carry exactly one nutrient. Collapsing them would throw away the
# measurement-level detail the data type exists to preserve.
DEDUP_EXEMPT_TYPES = frozenset({"sub_sample_food"})


def relisting_candidates(food, *, id_col="fdc_id", desc_col="description",
                         type_col="data_type"):
    """fdc_ids that could possibly be re-listings: a duplicated description.

    A necessary condition, not a sufficient one - the value vectors still have to
    match. Exposed separately so a caller can narrow what it reads out of the
    bucketed store before paying for the value comparison: 12,027 of 62,900
    eligible foods qualify, so it removes 81% of the work.
    """
    ded = food.loc[~food[type_col].astype("string").isin(DEDUP_EXEMPT_TYPES)]
    dup = ded.loc[ded.duplicated(subset=[desc_col], keep=False), id_col]
    return set(dup.tolist())


def find_exact_relistings(food, values, *, id_col="fdc_id", desc_col="description",
                          type_col="data_type", src_col="source_code",
                          nutrient_col="nutrient_id", amount_col="amount"):
    """fdc_ids to drop: same source, same name, same complete value vector.

    Keeps the lowest fdc_id of each matching set. Keyed on source as well as name,
    because two sources reporting an identical value for the same food is the
    borrowed-value case (CNF derives from USDA SR) and is two reports, not one.

    `values` may be pre-narrowed to relisting_candidates(food); it is subset again
    here, so passing the whole store is also correct.
    """
    cand = relisting_candidates(food, id_col=id_col, desc_col=desc_col, type_col=type_col)
    if not cand:
        return set()
    v = values.loc[values[id_col].isin(cand), [id_col, nutrient_col, amount_col, src_col]]
    if v.empty:
        return set()
    v = v.sort_values([id_col, nutrient_col], kind="stable")
    kv = v[nutrient_col].astype(str) + ":" + v[amount_col].astype(str)
    sig = kv.groupby(v[id_col], sort=True).agg("|".join)
    src = v.groupby(id_col, sort=True)[src_col].first()
    desc = food.set_index(id_col)[desc_col]
    key = src.astype(str) + "\x00" + desc.loc[sig.index].astype(str) + "\x00" + sig
    return set(int(x) for x in key.index[key.duplicated(keep="first")])


def source_of_bucket_file(path):
    """Which reference DB a bucketed parquet came from.

    merge_phase8_v2.py writes each non-USDA source to its own `<src>_data.parquet`
    inside every bucket; the original FoodData Central rows sit in the `part-*.parquet`
    files, several per bucket. The filename is therefore the provenance - but only if
    every FDC part collapses onto one label, or two identical FDC foods landing in
    different parts read as two sources and stop matching each other.
    """
    b = path.name if hasattr(path, "name") else str(path).rsplit("/", 1)[-1]
    if b.startswith("part-"):
        return "fdc"
    return b.replace("_data.parquet", "").replace(".parquet", "")
