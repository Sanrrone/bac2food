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
    # `source_food_code` (the identifier the source itself publishes) has been in the export
    # since the re-key, and is documented in Data Records; this list had not been updated, so
    # the check failed the correct file.
    "food_nutrients.tsv": ["fdc_id", "source_food_code", "description", "canon", "data_type",
                           "food_category", "nutrient_id", "nutrient_name", "unit_name",
                           "amount", "source_db"],
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
    # 2026-08-18: NEVO rejoined the deposit under written RIVM permission (+51,522 values,
    # +2,328 foods, no new nutrients), and this baseline was also carrying pre-re-key
    # numbers from before migrate_fdc_ids.py. Both corrected together against the
    # measured export, so a drift from here is a real regression again.
    # 1,956,046 -> 1,929,050 rows and 116,693 -> 116,043 foods on 2026-08-19. Three
    # exclusions, all in export_food_nutrients: 10 source columns that are not composition
    # values (food_DBs/_common/non_nutrients.py), 4 Livsmedelsverket editing-layer copies,
    # and 74 exact re-listings. 572 BioFoodComp foods went with them: their ONLY value was
    # the "Latest Revision in Version" database stamp, so once it went they carried no
    # measurement and the export's inner join dropped them. nutrients 1,779 -> 1,769 is
    # exactly the 10 removed columns.
    # canon_blank: 8 FDC foods whose description is empty in food.parquet itself. They
    # cannot be canonicalized because there is nothing to canonicalize; the count is
    # pinned so a future ingest that starts dropping descriptions is caught here.
    # canon 20,329 -> 21,567 on 2026-08-20, with rows/foods/nutrients all unchanged.
    # canonicalize_food_name stopped folding composition-changing preparation states into
    # the base food: drying, juicing, frying and sweetening alter per-100 g values, and the
    # predictor unions a canon's nutrients by MAX, so one dehydrated variant was handing its
    # whole group the concentrated figures (apple read 12.4 g fibre against a true 3.2).
    # Only the grouping moved, which is why this is the single count that changed.
    # canon 21,567 -> 20.480 on 2026-08-21, again with rows/foods/nutrients unchanged.
    # Two classes of silent over-merge were undone. (1) A venue prefix could swallow the
    # dish: 44 distinct McDonald's items, from french fries to a side salad, shared the
    # single canon "mcdonald's". (2) Identity-bearing detail was being deleted rather than
    # kept - an ingredient parenthetical ("BIG MAC (without Big Mac Sauce)"), the w/wo and
    # w/ abbreviations for with/without that NEVO, CIQUAL and Livsmedelsverket use, and
    # salted vs unsalted, which is now an axis of its own because a single state slot let
    # "Cod, dried, salted" and "Cod, dried, unsalted" both resolve to "dried" and merge.
    # Then to 20,096, after unsalted became the unmarked default (a raw carrot is
    # unsalted, so 'carrot, unsalted' and 'carrot' were one food under two names)
    # and trim became an axis - 'veal' alone had absorbed 281 descriptions, pure
    # separable fat MAX-unioned with separable lean only.
    # A further fall to 20,480 on 2026-08-21, when 1,764 curated overrides
    # from five parallel domain audits corrected misnamed canons and folded the
    # duplicates they exposed (brands, typos, and species mislabels such as
    # WAFCT palm KERNEL oil filed as palm oil).
    # The count first ROSE to 22,090 as those merges were undone, then FELL to
    # 21,093 when _CANON_MERGES folded 968 cross-database duplicates: sixteen
    # national databases write the same food in different orders, so 'olive oil',
    # 'oil, olive' and 'oil olive' were three separate foods with their evidence
    # split three ways. Both moves are grouping only; no food was added or removed,
    # which is why rows/foods/nutrients are unchanged throughout.
    # 20,095 -> 19,986 on 2026-08-22, rows/foods/nutrients unchanged. The
    # five naming defects reported by hand turned out to be classes, so every
    # canon in the index was scored against them and the 2,216 that scored were
    # read. Six rules absorbed the mechanical part - a corporate entity chunk is
    # a manufacturer ("The COCA-COLA company, DASANI, water" was the canon "the
    # company, dasani"), "imported from the U.S.A." is provenance not identity,
    # CIQUAL's "-> ARCHIVE" is bookkeeping, "n/a" is a placeholder, the cultivar
    # regex was widened to reach accession codes (IRNAS n° 11, KARI/BN/, Texas
    # 17W, DMR-ESR-W), and "ready to feed" is normalised to the hyphenated form
    # 91 other canons already use. The remaining 378 are curated overrides, one
    # name each, from reading the 2,216 flagged canons twelve slices at a time.
    # The whole fall is duplicate names collapsing onto the correct one: 96
    # targets absorbed a defective twin ("ricy" into rice, "cauliflower, danish"
    # into cauliflower, a 97-row canon called "sauce" that is entirely bottled
    # spaghetti sauce into "sauce, pasta"), and none crosses from one food to
    # another. Requests to SPLIT an over-merged canon are a different and larger
    # job; they are filed in AUDIT_REMAINING_ROUND4.tsv, not applied here.
    #
    # Then 19,986 -> 20,558 the same day, and this one RISES because the pass
    # recovered distinctions rather than merging names. Five composition axes
    # were missing or incomplete, so the values either side of each were being
    # MAX-unioned into one canon:
    #   sugar        "Pears, stewed with sugar" and "stewed without sugar" were
    #                one canon of 222; 21 canons held both, among them apple
    #                (312) and yoghurt (288). Unsweetened joins unsalted as an
    #                unmarked label - detected, so the split holds, never printed.
    #   sodium       FDC writes the positive form as "sodium added", which was
    #                landing in the name ('blackeye pea, sodium added') instead
    #                of on the axis, and "reduced sodium" had no label at all.
    #   fat level    43 canons held both ends: skim yoghurt was being handed the
    #                full-fat figures.
    #   enrichment   stripped as a preparation word, so enriched and unenriched
    #                rice sat together in one 634-member canon.
    #   oil pack     draining pours off water and brine but not absorbed oil, so
    #                the medium now survives a drain for oil only; 107 rows of
    #                each had been merged into 'tuna'.
    # 673 rows moved onto a correct salt or sugar label alone. Storage form
    # ("refrigerated") went the other way and is no longer part of a name.
    # Three follow-on corrections are in the same number: bare "skim" was
    # matching inside "part-skim" and calling 288 rows of part-skim mozzarella
    # fat-free, a cheese that is roughly 16 % fat; "Selenium" was missing from
    # the analyte prefixes, so 217 lab rows kept it as their food name; and
    # FDC's "Yogurt, Greek, strawberry" form was pushing the flavour out of a
    # two-chunk canon, which is now fronted to "Greek yogurt, strawberry" the
    # way every other database writes it.
    #
    # Finally 20,601 -> 20,592, two opposite moves that nearly cancel. A sixth
    # axis labels the ~700 FDC rows reported on a DRY MATTER basis ("Beans, Dry,
    # Black (0% moisture)"), which are not per 100 g as eaten and were the
    # majority of several canons - 253 of the 261 members of 'bean, pinto,
    # dried'. They are kept and labelled rather than dropped, so they can no
    # longer MAX-union with the as-eaten values of the same bean. Against that,
    # 14 canons under a generic head ("nut, pistachio", "bean, mung") folded
    # into the family that already held their food. Three candidates of the same
    # mechanical shape were rejected on reading: a cranberry BEAN is borlotti
    # and not the fruit, a butter BEAN is a lima and not dairy, and lentil flour
    # is not generic flour.
    #
    # And 20,592 -> 20,744, again mostly recovery. Three more classes, all found
    # by asking the same question - what does a member say that the name does
    # not?
    #   category head   A taxonomic or menu chapter in front of the food is not
    #                   the food. Six heads qualified on the evidence that each
    #                   fronts dozens of DIFFERENT foods: fish (234 canons),
    #                   mollusk (57), spice (57), grain (47), game meat (30),
    #                   crustacean (29). Freeing the qualifier slot also let the
    #                   species through, so "Fish, salmon, atlantic" is now
    #                   'salmon, atlantic' rather than 'fish, salmon'. It settles
    #                   what "game" means too, which FDC is loose about: its game
    #                   chapter holds farmed bison, goat, horse and "rabbit,
    #                   domesticated". "deli-meat" was rejected - it carries a
    #                   curing claim the product name does not always repeat.
    #   anatomical part 70 of the 309 members of 'apple' were apple PEEL, which
    #                   carries several times the fibre of the flesh; 38 canons
    #                   held the same defect, and 15 of the 18 members of
    #                   'chicken fat' were skin.
    #   raw readiness   "ready to bake / fry" says the food is still raw, and in
    #                   FDC's phrasing it took the one qualifier slot a canon
    #                   has: 'tortilla, ready-to-bake or -fry' had lost the corn
    #                   versus flour distinction. Not applied to ready-to-eat /
    #                   -drink / -feed, where the claim marks the final form.
    # A note for the next pass: _CANON_MERGES is keyed on the canon the rules
    # produce, exactly as _CANON_OVERRIDES is, and 113 keys across the two
    # tables silently stopped firing when the new labels appeared. Check both.
    #
    # And 20,744 -> 20,992. This round asked the opposite question of the last
    # one: not "what does a member say that the name does not?" but "how many
    # names does one food have?" - the whole index was swept for pairs that name
    # the same thing twice.
    #   one spelling per   The composition markers were being printed in the
    #   marker             SOURCE's spelling, because _append_state declines to
    #                      add a label when the token it matched is still in the
    #                      name. 'cheese, low fat' stood beside 'cheese, low-fat'
    #                      (143 members), 'milk, skimmed' beside 'milk, nonfat'
    #                      and 'milk, fat-free', 'soy flour, defatted' beside
    #                      'soy flour, fat-free'. The wording is now dropped once
    #                      read, for the fat, fortification and pressing-grade
    #                      axes as it already was for salt and sugar - but only
    #                      when that axis actually fired, or "Sugar, refined"
    #                      loses its own name.
    #   head order         Two OPPOSITE house conventions, each following the
    #                      majority the index had already settled on. MATERIALS
    #                      read modifier-first ('X oil' 135 canons against 30 for
    #                      'oil, X'; juice 125/11, flour 190/23, milk 398/88),
    #                      SPECIES read genus-first ('salmon, coho' 11 against 1
    #                      for 'coho salmon'; herring 8/1, pike 6/2). 52 pairs
    #                      named one food twice.
    #   nut chapter        "nut" joins the six category heads on the same
    #                      evidence (31 canons), with its own guards: "Nuts,
    #                      formulated" is an aggregate, and Frida's "Nut, brazil"
    #                      / "Nut, pine" / "Nut, pea" are half a compound name.
    #   packing medium     Read BEFORE the state and struck out of the probe the
    #                      state is read from. "Peach, canned in pear juice" was
    #                      coming out 'peach juice'. FDC's "<medium> pack"
    #                      phrasing and Livsmedelsverket's "canned w/ brine" now
    #                      reach the same label as "canned in brine" - but
    #                      "prepared with water" is reconstitution, not a medium,
    #                      and is excluded. Which OIL is kept: sardines in olive,
    #                      sunflower and peanut oil carry different fatty acids.
    #   powder is dried    'turmeric, dried', 'turmeric, powder', 'turmeric,
    #                      powdered, dried', 'turmeric, ground' - five canons for
    #                      one spice. A whole chunk that is only "powder" is the
    #                      dried label spelled out; inside a name ("chili
    #                      powder") the word stays.
    #   ingredient guard   A marker inside a "with ..." clause grades the
    #                      INGREDIENT: "Babyfood, banana juice with low fat
    #                      yogurt" is not a low-fat babyfood. And the EARLIEST
    #                      fat marker wins, not the first pattern in the list -
    #                      "Milk, lowfat, fluid, 1% milkfat, with added nonfat
    #                      milk solids" was being called fat-free.
    #   character noise    The full-width asterisk STFCJ footnotes with, and the
    #                      acute accent Fineli writes possessives with - the
    #                      brand strip took "Kellogg" and left "\u00b4S" standing as
    #                      a word of its own.
    #
    # And 20,992 -> 21,425. Two classes, one a split and one a convergence.
    #   colour            Colour is composition wherever a pigment is the
    #                     nutrient. 77 canons held both a dark- and a
    #                     light-coloured member: 'rice' (636 members) carried
    #                     black, brown and red rice beside white, 'common bean'
    #                     (223) black beside white, and 'tea, infusion' green
    #                     beside black. Only a WHOLE chunk counts, which is what
    #                     keeps a cultivar name out - "Apples, raw, red
    #                     delicious" and "Potato tuber, Red LaSoda" name a
    #                     variety, not a colour. On a refined cereal or a sugar
    #                     "white" is the refining state and the unmarked one, so
    #                     it is detected but not printed; on a bean or a cabbage
    #                     it is a variety and it is. "Egg, white" is the albumen
    #                     and is excluded outright, and bare "light"/"dark" are
    #                     read only as poultry cuts, where FDC uses them.
    #   missing commas    NEVO and Livsmedelsverket write no punctuation at all
    #                     ("Beans broad raw", "Nuts macadamia unsalted",
    #                     "Crackers cream"), so 154 canons kept a plural head
    #                     welded to its modifier - 'beans broad' beside the
    #                     103-member 'broad bean'. Restoring the comma the
    #                     source omitted is enough for every rule downstream.
    #                     Guarded three ways: the head must be genuinely plural
    #                     (_NOT_PLURAL carries Brussels, Maroilles, Causses,
    #                     bitters, sports), what follows must not be a
    #                     connective, and it must not be a DISH - "Strawberries
    #                     tart" is not a strawberry.
    # Two traps this round, both from the same mistake in opposite directions:
    # "Chocolate, milk" is a milk chocolate BAR (four sources write it that way)
    # and folding it gave 'chocolate milk', merging 30 g fat/100 g of
    # confectionery into a 3 g/100 g drink; and re-keying an orphaned override
    # must never GENERALISE the key - "blackberries product" is NEVO's frozen
    # one, and moving the entry to bare 'blackberry' froze every blackberry.
    #
    # And 21,425 -> 21,803, from the largest MAX-union left in the index. 'beef'
    # held 1,291 members: brain, heart, liver, lungs and tongue beside chuck,
    # brisket, tenderloin and eye of round. Beef liver carries roughly 9,000 ug
    # RAE of vitamin A per 100 g against about none in muscle, and beef chuck
    # about four times the fat of eye of round - so the union handed every cut
    # of beef the liver's vitamin A and the chuck's fat.
    #   organ            296 organ rows across 44 canons, 60 of them in 'beef',
    #                    40 in 'veal', 39 in 'pork'. FDC files them behind the
    #                    head as "variety meats and by-products, liver", and the
    #                    breed strip-head dropped the chunk.
    #   cut              The head has to strip breeds ("Japanese beef cattle",
    #                    "dairy fattened steer", "Belgian Blue") or every one
    #                    becomes a canon - but it was taking the cut with them.
    # Both are read only as a whole chunk BEHIND the animal, and only in a meat
    # context: at the front the organ IS the food ("Kidney, boiled, salted"
    # canonicalised to 'salt'), and outside a meat context "round", "breast",
    # "plate" and "rib" are ordinary words - celery has ribs and broccoli has a
    # plate of stalks. A legume head is excluded outright, because "Common bean,
    # Kidney" is a kidney BEAN.
    #
    # And 21,803 -> 21,857, from FDC's lab rows, where a greedy regex had been
    # eating the food.
    #   NF suffix        _NF_SUFFIX_RE allowed a comma inside its body, so the
    #                    lazy quantifier still took the LEFTMOST start that could
    #                    complete the match: "Cheese, cheddar, natural shredded
    #                    sharp, store brand, GREAT VALUE (CA1,NE) - NFY120WVO"
    #                    lost everything from ", cheddar" on, and 241 rows of
    #                    branded cheddar sat in the bare 'cheese' canon.
    #   brand chunks     With the middle back, the ALL-CAPS brand chunks behind
    #                    the food had to go instead - but only on a row that
    #                    carried a sample code, because off one an ALL-CAPS
    #                    chunk can BE the food ("BURGER KING - HAM").
    #                    "store brand" is the absence of a brand, not one.
    #   flavour          A strawberry Greek yoghurt carries about twice the
    #                    sugar of the plain one, and FDC writes the flavour
    #                    inside the brand chunk ("CHOBANI STRAWBERRY NON-FAT"):
    #                    106 of the 238 members of 'greek yoghurt, fat-free'
    #                    were flavoured. Read only where the plain form is the
    #                    unmarked one, so "Tomatoes, orange" and "Melon, banana"
    #                    keep their own names.
    #   cheese order     Twenty cheese types were written both ways and the
    #                    head-first form won every one ('cheese, cottage' 146
    #                    members against 'cottage cheese, full fat'). A LIST,
    #                    not a pattern: "cauliflower cheese" is a dish.
    # One trap: "Starch," heads an FDC lab row the way "Minerals," does, but
    # AFCD files a real food as "Starch, potato" - _PANEL_HEADS leaves it out
    # for that reason. The sample code is the discriminator, and the strip has
    # to run while the head is still at the FRONT, or it takes the food out of
    # "Minerals, Sugar, Granulated, White - NFY040XEG".
    #
    # And 21,857 -> 22,031, from the plant-part axis and the retail cut.
    #   plant part       247 rows named the FLESH of a food whose canon did not,
    #                    60 named BRAN - 54 of them inside 'rice', where rice
    #                    bran carries about 20 g of fat against milled rice's
    #                    0.7 - 48 named a SPROUT and 38 a LEAF. "root" and
    #                    "tuber" are deliberately absent: for a carrot the root
    #                    IS the food. So is the leaf of a lettuce, which is why
    #                    the leaf label is suppressed on leafy heads and herbs.
    #                    BioFoodComp writes "pulp" where the others write
    #                    "flesh"; the two named one part under two canons.
    #   retail cut       FDC names the retail cut, not the primal - "Beef,
    #                    shoulder top blade steak" - so a whole-chunk test
    #                    missed 182 of the 579 members of 'beef'. The primal is
    #                    the label, because that is the granularity that changes
    #                    composition. The prefix has to be LAZY, or the pattern
    #                    swallows "spare " and files a spare rib as a rib.
    #
    # And 22,031 -> 21,831, the first round that took the count DOWN: it asked
    # only "how many names does one food have?" and answered it 205 times.
    #   part order       The cut and organ axes append their label behind the
    #                    animal, so every source writing the same thing without
    #                    a comma ("Beef round", "Beef liver") or the other way
    #                    round ("Liver, pork") had a canon of its own. An ox
    #                    liver is a beef liver and a calf's a veal one.
    #   compounds        One spelling each: 'bread, gluten free' beside
    #                    'gluten-free', 'barley flour, whole grain' beside
    #                    'wholegrain' (107 canons closed against 49 open).
    #   merge hygiene    _CANON_MERGES is consulted ONCE, so a chain a->b->c
    #                    only ever moved a to b, and a 2-cycle a->b->a swapped
    #                    two spellings without converging them - 27 of those had
    #                    accumulated. Every entry now points straight at its
    #                    group's representative, and 37 override values that
    #                    were themselves merge keys were followed through.
    # Three same-key groups are deliberately NOT merged: "milk chocolate" and
    # "chocolate milk" (a bar and a drink), "roast beef" and "beef roast" (a
    # cooked dish and a cut), and the low-lactose/lactose-free milk pair.
    #
    # And 21,831 -> 21,857. Three more classes, all measured the same way:
    #   maturity         An aged cheddar has lost water, so it carries more fat,
    #                    protein and sodium per 100 g than a mild one; 116 of
    #                    the 400 members of 'cheese, cheddar' said which they
    #                    were. Read ONLY on a cheese - everywhere else "mature"
    #                    is a ripeness stage that _PREP_RE already owns, and
    #                    reading it there split 150 rows off for nothing.
    #   without skin     Says the same thing "flesh" does, and 98 rows spelled
    #                    it that way inside a canon that also held the with-skin
    #                    rows. Much of an apple's fibre is in the skin.
    #   longissimus      The longissimus dorsi IS the loin muscle: FDC and CNF
    #                    name it anatomically on 115 pork rows and by the cut
    #                    everywhere else.
    #
    # And 21,857 -> 21,836. Two strips that were taking half a name:
    #   seed / kernel    _PREP_RE drops "seeds" as FDC's legume convention
    #                    ("Beans, black, mature seeds"), but it was dropping the
    #                    word wherever it stood: CIQUAL's "Seeds and peanuts,
    #                    dried" canonicalised to 'and peanut, dried', and seven
    #                    rows of PUMPKIN SEED (49 g fat, 19 g protein) sat in
    #                    the 'pumpkin' canon beside the flesh (0.1 g fat). Now
    #                    guarded on both sides - the seed-crop names in front of
    #                    it, and a following "and"/"or".
    #   accession code   FDC files its bean breeding accessions as a chunk of
    #                    their own ("Beans, Dry, Pink, 11F-8082 (0% moisture)")
    #                    and 44 became canons of one row each. Letters BETWEEN
    #                    digits are what makes it a code, so "0-50 mg calcium
    #                    per litre" and "20-30 g fat" are left alone.
    #
    # And 21,836 -> 21,801: sixteen regional-synonym pairs, both halves already
    # in the index and never meeting. maize/corn (92 canons against 248),
    # soya/soy (76/196), swede/rutabaga, capsicum/bell pepper, faba and fava and
    # broad bean, haricot/navy bean, linseed/flaxseed, cornflour/cornstarch,
    # rocket/arugula, pak choi/bok choy, bicarbonate of soda/baking soda.
    # Only pairs that are the SAME SPECIES went in. Deliberately left apart:
    # biscuit/cookie (a US biscuit is a bread), marrow/squash (marrow is also
    # the bone), saithe/pollock and hake/whiting (different species),
    # chicory/endive (Cichorium intybus against C. endivia), treacle tart (a
    # dish, not molasses) and prawn/shrimp (different families).
    #
    # 2026-08-22, UNITS. Three sources report amino acids in mg/100 g and the
    # export labels every value with FDC's unit for the nutrient id - it does
    # not carry the source's own unit anywhere. 43,224 amino-acid rows were
    # deposited a thousand times too high (BioFoodComp's quinoa read 873 g of
    # leucine per 100 g), and STFCJ's "Amino acids, total" column, also in mg,
    # was mapped onto nutrient id 1003, so 929 canons carried an impossible
    # protein maximum - Parmesan read 48,000 g/100 g beside its real 44.
    # Now 637 amino-acid rows above 5 g (all high-protein foods) and none
    # above 100. See food_DBs/_common/units.py; the reconciliation runs in each
    # ingester, not here, so the bucketed store the predictor reads is fixed
    # too. It runs in both directions: AFCD's organic acids and beta-carotene
    # were in g and µg against FDC ids whose unit is MG, so those were a
    # thousand times too LOW.
    #
    # The counts moved with it, and every part of the move is accounted for:
    #   rows    1,929,050 -> 1,929,627. AFCD reports each amino acid twice, per
    #           100 g and "(mg/gN)" - per gram of NITROGEN - and both landed on
    #           the same FDC id. A different BASIS is a different measurement,
    #           not a different unit, so the second is minted separately now.
    #           McCance's "Tryptophan/60" and its "/100g fa" fatty acids are the
    #           same shape and go the same way.
    #   foods   116,043 -> 116,053. Ten BioFoodComp foods the live store's build
    #           did not have; re-running the committed ingester finds them.
    #   nutrients 1,769 -> 1,749. Net of 21 minted ids added (18 AFCD per-gram-N
    #           columns, STFCJ's "Amino acids, total" now on its own id instead
    #           of overwriting protein, 2 McCance) and 41 gone, 40 of them
    #           BioFoodComp ids the re-run renumbered.
    #   canon   21,801 -> 21,803, from the ten new foods.
    #
    # Minted ids are now assigned from a REGISTRY rather than a counter
    # (food_DBs/_common/minted.py). They were positional, so adding or dropping
    # one column renumbered the whole block - and enzyme_substrate_chebi.tsv
    # references them BY NUMBER, so one re-run silently pointed 20 of those at
    # nothing. Two consecutive runs of an ingester now mint zero new ids.
    #
    # 21,803 -> 21,799 on a review pass. The maturity grade is rarely a whole
    # chunk - FDC writes "cheddar, natural shredded sharp" and puts it inside
    # the brand on "Cheddar cheese, sliced, SARGENTO SHARP" - so it is matched
    # as a WORD now, and read off the description as it ARRIVED, before the
    # brand chunk carrying it is stripped. 187 of the 284 members of
    # 'cheese, cheddar' said which they were; they are now 124 mature and 63
    # mild. The cheese gate is what makes a bare word safe.
    #
    # 21,799 -> 22,152 in round 14, and the row, food and nutrient counts do not
    # move with it: only NAMES changed, no value did. 3,721 of the 53,983
    # descriptions landed on a different canon. The count went UP because two of
    # the round's readings SPLIT canons that were mixing foods, and both are the
    # defect the axes exist for - a canon's members are averaged, so a
    # distinction the name does not carry is a distinction the numbers lose:
    #   * the fat level stated as a NUMBER. 1,163 rows carry one and _QUANT_RE
    #     deleted every one, so a 17%-fat edam and a 30% one were one canon and
    #     FDC's ground beef - sold at 5, 7, 10, 15, 20 and 30% fat - sat whole
    #     inside the bare 'beef'. The number is kept rather than banded: no one
    #     threshold is right for two foods at once, since 17% fat is a LOW-fat
    #     cheese and a high-fat yoghurt. A stated range is kept whole.
    #   * skin. FDC's "meat only" is Fineli's "without skin" and CNF's "light
    #     meat only"; "meat and skin" is the whole bird and now says nothing.
    #     A skinless chicken breast carries about a sixth of the fat of one with
    #     the skin on. Canons mixing marked and unmarked members: 89 -> 3.
    # Against those, the round CONVERGED as much as it split: the USDA poultry
    # class term ("Chicken, broilers or fryers") that split one bird across 79
    # canons, "tinned" (50 canons that stood apart from their "canned" twins),
    # 46 same-food alternations ("Yoghurt or fermented milk", "Macaroni or
    # noodles with cheese"), Fineli's trademarks (twelve canons of rye bread,
    # and no bare 'rye bread' at all), and one grind spelled four ways.
    # Duplicate-name groups are back to 3, all of them deliberately blocked.
    #
    # 22,152 -> 22,093 in round 15. Rows, foods and nutrients again do not move:
    # 1,790 of the 53,983 descriptions changed canon and not one value did. The
    # net is DOWN because this round was mostly convergence:
    #   * "no added salt" / "no added sugars" was matched only from "added"
    #     onwards, so the negative was left standing and 37 canons ended in the
    #     bare word "no" - 'peanut, no' beside 'peanut'. The front-loaded
    #     POSITIVE ("Pasta, cooked, with added salt") was not read at all, so
    #     nine salted rows were merging INTO their unsalted canon.
    #   * a qualifier that states the unmarked case: 36 foods carried both "X"
    #     and "X, plain" ('almond milk', 'bagel', 'butter', 'tofu'), and 20 both
    #     "X" and "X, mixed species" ('shrimp', 'squid', 'trout').
    #   * one part, four spellings. 'chicken, flesh', 'chicken, flesh only',
    #     'chicken, meat' and (on fruit) 'coconut, pulp' / 'jujube, fruit flesh'
    #     were the same part; 16 foods carried two of them and seven canons
    #     carried two at once ('baobab fruit/monkey bread, pulp, flesh').
    #   * FDC's alcohol class noun ('alcoholic beverage, beer' beside 'beer'),
    #     CNF's traditional-foods marker "native" (171 rows, and in the game
    #     rows it took the head outright - 'native, caribou, liver'), and
    #     "home-made" against "homemade" (54 canons against 264).
    # Two readings SPLIT, both for the usual reason - the members are averaged,
    # so a distinction the name does not carry is one the numbers lose:
    #   * wild or farmed. 'salmon, atlantic' held 8 wild members and 10 farmed,
    #     'trout, rainbow' 4 and 7. Farmed Atlantic salmon runs to about twice
    #     the fat of wild. Printed only for foods that come both ways: 460 heads
    #     carry a wild row and only 16 a farmed one too, the rest being capture
    #     fisheries and game where wild is the unmarked case.
    #   * the fat a food was COOKED in, where the source names it. 19 canons
    #     averaged two or three different fats ('potato chip, homemade, fried'
    #     held corn, rapeseed and sunflower); the unnamed ones ("vegetable",
    #     "blended") name nothing and are struck. Clash count 19 -> 0.
    # 22,093 -> 22,094 on the re-run: three residues the index build's own
    # diagnostic surfaced. CIQUAL's approximation word was being left behind by
    # the fat figure it qualifies ('tomme cheese, around, 13% fat'), Frida
    # splits a figure from its unit across two chunks ("Milk, whole, 3.5, (UHT),
    # % fat", which had produced the meaningless 'milk, % fat'), and the word
    # "whole" was printing the fat claim a second time beside the figure.
    # Ten curated merges were also reversed, which renames without regrouping so
    # the count does not move: each had folded the spelling MORE rules produce
    # into the one fewer do, and each of those was the form that reads as a
    # run-on or puts the material before its food - 'flour soy' for soy flour,
    # 'sausage chorizo', 'nutmeg ground', 'cheese stilton', 'radish black'.
    "food_nutrients.tsv": {"rows": 1_929_627, "foods": 116_053, "nutrients": 1_749,
                           "canon": 22_094, "canon_blank": 8},
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

# RIVM's permission (18 Aug 2026) covers NEVO's DERIVED values, which now ship inside
# food_nutrients.tsv, and not the NEVO release itself: they asked that users be pointed at their
# download page instead. So this guard narrows rather than lifts — no NEVO *file* may appear in
# the deposit directory. It stays filename-level because the failure it guards against is
# someone dropping the release back in by hand, which is how it got there the first time.
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
    foods, nutrients, sources, canons = set(), set(), set(), set()
    bad_amount = 0
    blank_canon = set()
    for (fdc_id, code, desc, canon, dtype, cat, nid, nname, unit, amount,
         src) in scan(p, SCHEMA[p.name]):
        rows += 1
        foods.add(fdc_id); nutrients.add(nid); sources.add(src)
        canons.add(canon) if canon else blank_canon.add(fdc_id)
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
    check(len(canons) == e["canon"], f"distinct canon names == {e['canon']:,}",
          f"got {len(canons):,}")
    check(len(blank_canon) == e["canon_blank"],
          f"foods with no canon == {e['canon_blank']} (blank description upstream)",
          f"got {len(blank_canon)}")
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
          f"found {stowaway} — the derived values ship inside food_nutrients.tsv; the release "
          f"itself does not, and users are pointed at the provider's download page"
          if stowaway else f"no release file for: {', '.join(WITHHELD)}")

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
