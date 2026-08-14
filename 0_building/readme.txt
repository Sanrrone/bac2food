in eggnog folder (see eggnog/readme.txt)
python 1.0_eggnog_ec_substrates_parser.py --in ec_species_substrate.tsv --out 2_digest_dict.tsv --dedup
python 1.5_reactions_to_digest.py

in chebi folder
python dict_to_chebi.py --digest ../eggnog/2_digest_norm.tsv --obo chebi.obo --out digest_to_chebi.tsv
# NORM, not dict: 1.5_reactions_to_digest.py is where the ubiquitous cofactors are removed.
# See eggnog/readme.txt. Output must be 218,378 rows (BRENDA 2026.1 download; it was
# 95,326 under the SPARQL scrape, and 94,886 once the cofactor filter learned NAD(P)+/
# NAD(P)H and bracketed macromolecules).

in 0_building folder
python 0_name_normalizer.py --input /data/bac2food/nutrient.csv --format csv --column name --output 0_nutrient.normalized.tsv
python 1_nutrient_expansion.py --nutrients 0_nutrient.normalized.tsv --out_expansions 1_expanded_nutrients.tsv --out_units 1_nutrients_units.tsv
python 2_nutri2chebi_from_obo.py --nutrients 1_nutrients_units.tsv --chebi-obo ../chebi/chebi.obo --out 2_nutrient_to_chebi.tsv --min-score 30
python 3_nutrient_to_ec.py --nutrient_best 2_nutrient_to_chebi.tsv --digest_chebi  ../chebi/digest_to_chebi.tsv --chebi_obo ../chebi/chebi.obo --out 3_nutrient_to_ec.tsv --max_cost 1.5 --w_is_a 1.0 --w_conj 0.5 --include_simple_sugars --extra_seeds extra_bacterial_seeds.tsv --override_seeds extra_nutrient_chebi.tsv --live_nutrients /data/bac2food/live_nutrients.tsv
# --override_seeds is NOT interchangeable with --extra_seeds: the first REPLACES a nutrient's
# seeds, the second only adds to them. This line carried neither until 2026-08-06, so the 77
# curated rows in extra_nutrient_chebi.tsv were absent from every shipped map even though
# 4_predict/README.md documented passing them. See extra_nutrient_chebi.readme.txt.

# 2026-08 CHANGES TO THE NUTRIENT->ENZYME LINK (both are script fixes, not TSV edits)
#
# 2_nutri2chebi_from_obo.py — INFOODS mineral tagnames are now locked to their element.
#   Sources that name mineral columns "CS(mcg)", "K(mg)", "MG(mg)" collide with the IUPAC
#   one-letter amino-acid codes ChEBI uses for peptide synonyms, so the shipped map had
#   K(mg)->lysine, P(mg)->L-proline, CS(mcg)->Cys-Ser, TI(mcg)->Thr-Ile, MG(mg)->Met-Gly,
#   S(mg)->strange quark, B(mcg)->bottom quark. 28 columns are affected and 16 of them had
#   acquired enzymes (potassium alone carried 92 EC). The element reading now wins and the
#   name passes are skipped. FDC-style "Calcium, Ca" is deliberately NOT locked: it already
#   resolves correctly and with better specificity (Fluoride, F -> the fluoride ion).
#
# 3_nutrient_to_ec.py — the walk now uses the whole structural ChEBI vocabulary.
#   chebi.obo declares nine relationship types; the walk used three. Adding
#   has_functional_parent, has_parent_hydride, is_substituent_group_from, is_tautomer_of,
#   is_enantiomer_of and has_part reaches +33 nutrients, chiefly the flavonoid glycosides
#   and hydroxycinnamate esters whose deglycosylation IS the gut-bacterial reaction of
#   interest (myricetin 3-O-glucoside -> myricetin, 5-feruloylquinic acid -> quinic acid).
#   has_role is excluded on purpose: it is semantic, not structural, and would link every
#   polyphenol to every other through a shared role node.
#   Net effect on the map: 548 -> 565 nutrients with an enzyme, of which 16 removed were
#   the mineral artefacts above, so the real gain is +33 and the map is smaller in error.
#
# Regenerating both is ~70 s and ~1.2 GB peak. The outputs are inputs to the enzymes
# export and to every cohort prediction, so a rebuild invalidates the reported match rate.

--live_nutrients is REQUIRED for a correct map. The nutrient catalog has ~3,300 entries but only
1,788 carry a measured value in the composition table; without the filter the map emits edges to
the other ~243, which can never reach a food but were still counted as "in the model". That
inflated every downstream coverage figure (the cohort match rate read 45.4% instead of 36.9%).
Regenerate the id list after any rebuild of the food export:

  python -c "import csv;ids={r[4] for r in csv.reader(open('/data/bac2food/exports/food_nutrients.tsv'),delimiter='\t')};\
ids.discard('nutrient_id');open('/data/bac2food/live_nutrients.tsv','w').write('nutrient_id\n'+'\n'.join(sorted(ids,key=int)))"

Then rebuild the dependent export:
  python ../5_export/export_resources.py --only enzymes --out_dir /data/bac2food/exports
and check it:
  python ../5_export/verify_exports.py
(--out_dir is required and has no default; here we DO mean the deposit itself.)

# Flag notes:
#   --include_simple_sugars       Allow glucose / fructose / lactose / sucrose / maltose /
#                                 galactose / mannitol / sorbitol through the strict_prebiotic
#                                 filter (needed for bacterial-substrate modelling).
#   --include_amino_acids         Same for free amino acids (proteolytic bacteria).
#   --max_hub_ec <n>              Promiscuity guard (default 40 = 99th percentile of the
#                                 digest). A ChEBI node carrying more EC than this is a
#                                 ubiquitous metabolite (D-glucose 125, L-glutamate 134,
#                                 dioxygen 617). Structural inferences may not route
#                                 through one; identity relations (exact/conjugate/
#                                 tautomer) still may. Without it punicalagin inherits
#                                 125 enzymes for containing a glucose. 0 disables.
#   --no_structural               Reproduce the pre-2026 graph (is_a + conjugate only).
#   --w_taut/--w_enant/           Weights for the structural relations added in 2026.
#     --w_deriv_up/--w_deriv_down/  --w_deriv_down defaults ABOVE --max_cost, i.e. the
#     --w_part                      scaffold->derivative direction is off; measured at
#                                   +1 nutrient for +5,748 rows, so it is not worth it.
#   --extra_seeds <tsv>           Append curated (nutrient_id, name, chebi_id) rows for
#                                 substrates absent from FDC (HMOs, GlcNAc, GalNAc, sialic
#                                 acid, fucose, xylan, arabinoxylan, dextran, pullulan,
#                                 cellobiose, alginate, agarose). See extra_bacterial_seeds.tsv.
#                                 NOTE: only seeds whose nutrient_id ALSO exists in FDC
#                                 food_nutrient data will affect food rankings; novel
#                                 substrate IDs (200001+) currently propagate EC mappings
#                                 but cannot drive food scores until matching
#                                 food_nutrient rows exist.




# 2026-08-05: branded foods left the composition export
#
# 5_export/export_resources.py now drops FDC branded_food entries by default
# (--keep_branded restores them). That changes which nutrients carry a measured value,
# so the nutrient->EC map is NOT valid across the change. Rebuild in this order:
#
#   1. python ../5_export/export_resources.py --only foods --out_dir /data/bac2food/exports
#   2. regenerate /data/bac2food/live_nutrients.tsv (recipe above)   1,788 -> 1,779 ids
#   3. rerun the 3_nutrient_to_ec.py command above                   533 -> 531 nutrients
#   4. python ../5_export/export_resources.py --only enzymes ...     in_model 4,094 -> 4,065
#   5. python ../5_export/verify_exports.py                          expect 24/24
#
# Skipping 2-4 leaves the digest citing nutrients no food measures. Only Inositol (1181)
# and EGCG (1368) were affected, both branded-label-only; EGCG loses nothing because
# Phenol-Explorer measures it analytically as nutrient 240010, but inositol's 29 EC rows
# are a real loss, since the only analytical source for it is WAFCT's ambiguous
# "Phytic Acid / Myo-Inositol" column, which the matcher correctly refuses to resolve.
#
# Every cohort prediction must also be re-run: wipe the index cache first
#   rm /data/bac2food/index_modeled/{static_food_meta.pkl,*.parquet}
# The predictor already dropped branded, so the FOOD side is unchanged, but the map is not:
# modeled nutrients per run went 311 -> 309 and food rankings shifted.
