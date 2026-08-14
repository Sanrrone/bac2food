eggnog/ — EC -> substrate digest (BRENDA) + bacterium -> EC rebuild (eggNOG v7)
==============================================================================

TWO LAYERS, TWO SOURCES -- read this before "updating" anything (clarified 2026-08-03).

  * EC -> substrate   from BRENDA. KEPT, never deprecated.
  * bacterium -> EC   from eggNOG (v6 originally, v7 now).

BRENDA was dropped for ONE role only: its DSMZ SPARQL scrape 3_normalized_species.tsv
supplied species -> enzyme and under-reported every organism's enzymes, so eggNOG replaced
*that* link (see ../5_export/README.md). BRENDA's EC -> substrate digest survived the swap
and is still what the shipped enzyme_substrate_chebi.tsv is built from. The folder name is a
leftover of that role swap; the top-level README briefly mislabelled the surviving digest
"eggNOG v6", which is what this note corrects.

Why the substrate layer CANNOT be moved to eggNOG v7: eggNOG is an orthology resource and
publishes no substrate data at any version -- there is nothing to move it to. Evidence that
the digest is BRENDA's, verified against ec_species_substrate.tar.xz: its most frequent
"substrate" value is BRENDA's placeholder "more" (1,027,040 of 35,355,079 rows); the
substrates are literature assay compounds (e.g. "ethyl (2R)-(4-fluorophenyl)(hydroxy)acetate"
for EC 1.1.1.1); and only 4.7% of them match a Rhea/ChEBI name, ruling out a Rhea origin.
BRENDA is CC BY 4.0 -- redistributable, but it must be cited.

The extract carries three columns (ec, species, substrate). Both digest builders keep only
(ec, substrate); the species column -- the deprecated link -- is dropped before anything
downstream sees it.

Extract taken 2026-02-19. The note that used to sit here -- "one release behind BRENDA 2026.1
(March 2026)" -- was wrong twice over, corrected 2026-08-06 against the official download:
BRENDA 2026.1 is dated 2026-02-10, i.e. it PRECEDES the scrape, and it is still the current
release as of 2026-08-06. So the scrape is not a release behind; there is no newer release.

What it IS, measured against brenda_2026_1.json, is INCOMPLETE. The DSMZ SPARQL endpoint
yielded substrates for 5,263 EC where the download carries 6,901, and it covers only 63.5% of
the download's (ec, substrate) pairs verbatim. That is the same lossiness that got the
endpoint's species -> enzyme link replaced by eggNOG; it affects the substrate link too.
The endpoint is unreachable as of 2026-08-06 (no route to sparql.brenda-enzymes.org), so the
official download is now the only way to refresh this layer.

Note the two sources do NOT share a vocabulary. The endpoint served normalized ligand names
("NADP(+)", "water"); the download serves raw literature strings ("NADP+", "H2O"). A refresh
therefore changes how ../chebi/dict_to_chebi.py matches, and the ubiquitous-cofactor list in
1.5_reactions_to_digest.py must be re-checked against the new spellings before trusting any
row count. BRENDA data is CC BY 4.0 (verified 2026-08-06 at brenda-enzymes.org/license.php).

Build the substrate digest (EC -> substrate) from the BRENDA extract; feeds ../chebi/ and
../0_building/. Audited 2026-08-03 and found correct: recomputing unique (ec, substrate)
pairs from the 35,355,079-row source reproduces all 105,933 shipped pairs (one cosmetic
difference, a stripped control character). 453,474 source rows carry an EMPTY ec field and
are silently dropped by the `^\d` filter -- correct, but silent. The ubiquitous/"more"
filter in 1.5_ works: zero such rows reach the exports. Note 1.5_'s normalization is
ADDITIVE, not destructive -- ../chebi/dict_to_chebi.py queries the raw name and the
normalized name (and further variants), so stripping "D-"/"L-" cannot lose a match:

CURRENT build, from the official download (2026-08-06 onward). The release file is NOT in
this repository -- BRENDA's own README states the full contents are copyright-protected, so
what ships here is the derived digest (2_digest_dict.tsv, 2_digest_norm.tsv) and not its
source. Fetch brenda_2026_1.json.tar.gz from https://www.brenda-enzymes.org/download.php
first if you actually need to rebuild:

  tar -xzf brenda_2026_1.json.tar.gz -C /data/bac2food/brenda_2026_1
  python 1.0b_brenda_json_to_digest.py --json /data/bac2food/brenda_2026_1/brenda_2026_1.json \
                                       --out 2_digest_dict.tsv
  python 1.5_reactions_to_digest.py          # 2_digest_dict.tsv -> 2_digest_norm.tsv

RETIRED build, from the SPARQL scrape (kept so an old map can be reproduced):

  tar -xJf ec_species_substrate.tar.xz
  python 1.0_eggnog_ec_substrates_parser.py --in ec_species_substrate.tsv --out 2_digest_dict.tsv --dedup

  (2_gen_digestdict.py is a chunked alternative to 1.0_ for the same input.)

Two traps cost a full afternoon when 1.0b_ was written; both are silent, and both are now
guarded in that script:
  * BRENDA appends reversibility to a reaction value in braces -- "A + NADH = B + NAD+ {r}".
    Split the equation without stripping it and the marker stays glued to the last compound,
    giving phantom substrates ("nad+ {r}", "h2o {r}", "? {r}") that the cofactor filter below
    cannot catch, because it quite reasonably lists "nad+" and not "nad+ {r}".
  * csv.writer defaults to the excel dialect's "\r\n", so every substrate ends up as
    "NAD(P)+\r" and NO exact-name lookup downstream can ever match it. Use
    lineterminator="\n".
Take BOTH sides of the "=" as well. BRENDA reactions are largely reversible and the question
here is "can this enzyme act on this compound", not "in which direction"; left-side-only
loses 36,663 pairs, e.g. EC 1.1.1.2 on 3-methoxybenzaldehyde, which is written as a product.

Downstream:
  ../chebi/dict_to_chebi.py --digest 2_digest_norm.tsv --obo chebi.obo --out digest_to_chebi.tsv
      ^^^^^^^^^^^^^^^^ 2_digest_NORM, never 2_digest_dict. This line said `dict` until
      2026-08-05 and it is the one mistake here that silently corrupts the resource:
      1.5_ is what strips the ubiquitous cofactors, so feeding the unfiltered dict adds
      10,607 rows and resolves 792 `hydron`, 247 GTP, 172 UTP, 134 dATP, 102 coenzyme A,
      83 ammonia, 63 FAD and 57 CO2 as though they were dietary substrates.
      dict_to_chebi.py's own IGNORE_CHEBI_IDS covers only 7 compounds (NAD(P)(H), ATP,
      ADP, water) and will NOT save you. Sanity check: the output must have 218,378 rows,
      matching 2_digest_norm.tsv, not 105,933.
  ../0_building/3_nutrient_to_ec.py ... --digest_chebi ../chebi/digest_to_chebi.tsv

Look up the enzymes of one organism (feeds ../1_query/ec2food.py --enzyme_tsv):

  python 4_bact_proteins_query.py --query "Bacteroides thetaiotaomicron VPI-5482" --out btheta.tsv


--------------------------------------------------------------------
bacterium -> EC : two routes (v6 direct, v7 via KEGG KO)
--------------------------------------------------------------------

ROUTE A (v6, original). /data/bac2food/bact_ec.tsv — eggNOG v6 carried an EC annotation
directly. ../5_export/export_resources.py flattens it into exports/species_enzymes.tsv.
Preserved as exports/species_enzymes.v6.tsv. Names follow eggNOG's older taxonomy
(Lactobacillus, not Lacticaseibacillus); tax_id is the stable join key.

ROUTE B (v7, current export). eggNOG v7 REMOVED the direct EC annotation — its bulk file
e7.og_info_kegg_go.tsv carries KEGG Orthology (KO) ids instead. So the chain gains a hop:

    6.0_kegg_ko_to_ec.py          KEGG list/ko -> kegg_ko_ec.tsv   (the KO -> EC bridge)
    6.1_eggnog7_species_enzymes.py  OG members + KOs + bridge -> exports/species_enzymes.tsv

  python 6.0_kegg_ko_to_ec.py --out /data/bac2food/kegg_ko_ec.tsv
  python 6.1_eggnog7_species_enzymes.py --sample 300000 --stats_only     # prototype
  python 6.1_eggnog7_species_enzymes.py --out /data/bac2food/exports/species_enzymes.tsv

Inputs, all on /data (none of the 12 GB FASTA / 27 GB tree downloads are needed):
    e7.og_info_kegg_go.tsv  2.7 GB   col6 = "taxid.protein,..."  col7 = "K01046|30.00;..."
    e7.taxid_info.tsv.gz    1.2 MB   New_Taxid, Old_Taxid, Sci_Name, Rank, lineages
    kegg_ko_ec.tsv          from 6.0 (one ~2 MB KEGG request; /link/ec/ko is not served)
Full build ~1 min, ~2 GB RAM.

Two things that decide correctness (do not "simplify" them away):
  * Member proteins are prefixed with Old_Taxid, so that is the JOIN key, but the emitted
    tax_id is New_Taxid — this is the point of v7, it retires the stale-taxonomy caveat.
    NCBI merged some taxa (41 New_Taxids absorb 49 Old_Taxids), so the canonical index is
    keyed on New_Taxid; keying on the old id emits duplicate (tax_id, EC) rows.
  * --min_consensus defaults to 0 and should stay there. Filtering on col7's consensus %
    looks prudent but measurably wrecks coverage (>=20 already falls below v6).


--------------------------------------------------------------------
v6 vs v7 — measured, so the trade-off is a choice and not a surprise
--------------------------------------------------------------------

                                    v6 (direct EC)     v7 (via KEGG KO)
  rows                                 9,632,315          20,557,730
  organisms                                3,176              10,751
  distinct EC numbers                      4,291               4,819
  EC coverage of the digest                72.6%               79.7%
  EC coverage of nutrient-reaching ECs     75.9%               82.1%
  median ECs per organism                  3,155               2,010

v7 is BROADER: 3.4x the organisms (22 v6 taxa lost, 7,597 new), +528 distinct ECs, and it
covers more of the ECs that actually reach an FDC nutrient — the metric that matters, since
an EC with no nutrient link is inert in bac2food.

v7 is SHALLOWER PER ORGANISM: all 3,154 shared taxa lose EC depth (100%), median 3,155 ->
2,010. That is the bridge tax — only ~39% of KEGG KOs are enzymes at all, and KO -> EC is
many-to-many, so the v7 EC set is NOT a superset of v6's (40 digest ECs, 20 nutrient-reaching
ECs are lost outright). Predictions for well-characterised gut species rest on fewer enzymes
than before, while the long tail of species gains coverage it never had.

To roll back to v6:  regenerate from the retained legacy source —
  python ../5_export/export_resources.py --only species \
         --bact_ec /data/bac2food/bact_ec.tsv --out_dir /data/bac2food/v6_rebuild
  then point 4_predict/parameters.yaml:bact_ec_ref at the result.
--out_dir is REQUIRED (argparse rejects the run without it, before any I/O). Point it at a
scratch directory, NOT at /data/bac2food/exports: the output is named species_enzymes.tsv
regardless of the input version, so writing a v6 build into the deposit replaces the shipped
v7 layer. It used to default to the deposit, which made that overwrite silent.
(exports/species_enzymes.v6.tsv was retired from the deposit once v7 became the only
shipped layer; bact_ec.tsv is kept precisely so this stays reproducible.)

The EC -> substrate -> ChEBI layer is untouched by any of this: it is keyed on EC, so
exports/enzyme_substrate_chebi.tsv is unchanged and both routes plug into it identically.
