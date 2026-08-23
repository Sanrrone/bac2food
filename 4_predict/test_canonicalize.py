#!/usr/bin/env python3
"""Regression cases for canonicalize_food_name.

Run: python3 4_predict/test_canonicalize.py

Every case here is one that was wrong at some point and was fixed. The canon name
decides which fdc_ids share a nutrient vector, and build_modeled_index AVERAGES that
vector across members, so a canon that merges two different foods pulls both
toward a number that is right for neither. That failure is silent - it produces plausible
numbers - which is why these are pinned rather than left to inspection.
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("b2f", HERE / "bac2food_predict.py")
b2f = importlib.util.module_from_spec(spec)
sys.modules["b2f"] = b2f
try:
    spec.loader.exec_module(b2f)
except SystemExit:
    pass
canon = b2f.canonicalize_food_name

CASES = [
    # --- venue / brand / programme prefixes name who served the food, not what it is
    ("Carrabba's Italian Grill, spaghetti with meat sauce", "spaghetti with meat sauce"),
    ("Olive Garden, spaghetti with meat sauce",             "spaghetti with meat sauce"),
    ("School Lunch, chicken nuggets, frozen",               "chicken nugget"),
    ("McDonald's, filet-o-fish",       "fish fillet sandwich with cheese and tartar sauce"),
    ("KFC, biscuit",                                        "biscuit"),
    ('Domino\'s 14" cheese pizza, original crust',          "cheese pizza, original crust"),
    ("Chinese Restaurant, egg rolls, vegetarian",           "egg roll, vegetarian"),
    # --- the venue can sit BEHIND an FDC panel or group prefix
    ("Beverages, WENDY'S, tea, ready-to-drink, unsweetened", "tea, ready-to-drink"),
    ("Minerals, Chinese Restaurant, Fried Rice (No Meat) (NY) - NFY1209PL", "rice, fried"),
    # --- a venue is not always a prefix; a chunk that is only a venue goes wherever it sits
    ("Sweet And Sour Pork, Chinese restaurant (IN1,NY1) - CY09041", "sweet and sour pork"),
    ("Gravy, CAMPBELL'S, chicken",                                  "gravy, chicken"),
    # --- separators FDC actually uses: hyphen, and none at all
    ("BURGER KING - HAM",                                "ham"),
    ("CHICK-FIL-A - AMERICAN CHEESE",                    "american cheese"),
    ("McDONALD'S Bacon Ranch Salad with Crispy Chicken", "bacon ranch salad with crispy chicken"),
    # --- generic heads need punctuation, so an ordinary word is safe
    ("Restaurant food, general",                         "restaurant food, general"),
    ("Alcohol, beer",                                    "beer"),
    # --- sample codes come in two prefixes
    ("Broccoli, raw (IN1,NY1) - CY0906E",                "broccoli"),
    ("Salsa, TOSTITOS CHUNKY, MEDIUM - NFY090KVS",       "salsa"),
    ("Cypress nut",                                      "cypress nut"),   # not a CY code
    # --- panel heads that are themselves foods must survive _strip_panel_prefix
    ("Sugars, brown",         "sugar, brown"),
    ("Salt, table",           "salt, table"),
    ("Water, non-carbonated", "water, non-carbonated"),
    ("Pectin, liquid",        "pectin, liquid"),
    # --- label copy is not composition
    # CONFLICT, flagged to the project owner 2026-08-21. On 2026-08-20 they asked
    # for this exact string to keep "pillsbury" and lose only the trailing
    # "artificial flavor". On 2026-08-21 they asked for brands to be stripped
    # ("SMART BEAT ..." -> "margarine"). The later, more general rule wins here.
    # Drop "pillsbury" from _BRAND_NAMES to restore the earlier behaviour.
    ("Pillsbury Golden Layer Buttermilk Biscuit, Artificial Flavor",
     "golden layer buttermilk biscuit"),
    ("Yogurt, natural flavors",                             "yoghurt"),   # spelling converged
    # --- separators left behind by the prep strip
    ("Asparagus, /",                                        "asparagus"),
    ("Broccoli, /",                                         "broccoli"),
    # --- NOT venues: real foods that merely look possessive or brand-like
    ("Chipotle dip, yogurt based",   "chipotle dip, yoghurt based"),  # the pepper, not the chain
    ("Shepherd's pie, beef",         "shepherd's pie, beef"),         # a dish
    ("Cat's whisker, leaves, raw",   "cat's whisker, leaf"),          # a vegetable
    ("Jew's mallow, raw",            "jew's mallow"),                 # a vegetable
    ("Cafe au lait",                 "cafe au lait"),                 # a drink
    ("Babyfood, apple",              "babyfood, apple"),              # genuinely not raw apple
    # --- composition-changing state must survive; form must not
    ("Apples, raw, with skin",                 "apple"),
    ("Carrots, sliced, frozen, unprepared",    "carrot"),
    ("carrot, raw",                            "carrot"),
    ("Apples, dried, sulfured, uncooked",      "apple, dried"),
    ("Apple [Dessert], jam [Phenol-Explorer]", "apple, sweetened"),
    ("Carrot juice, canned",                   "carrot juice"),
    ("Beef, ground, 80/20, pan-fried",         "beef, ground, fried"),
    ("Milk, dry, whole",                       "milk, dried, full-fat"),
    ("Chili powder",                           "chili powder"),        # already dry; no suffix
    ("Soup, onion, dry, mix, prepared with water", "soup, onion"),     # the reconstituted soup
    # --- a leaf is not a tuber
    ("Cassava, leaves, raw [WAFCT]",  "cassava, leaf"),
    ("Cassava, tuber, raw [WAFCT]",   "cassava"),
    ("Rice, bran, crude",             "rice bran"),   # _CANON_MERGES: majority form
    # --- source tags and cultivar codes never reach the canon
    ("Almonds, raw [Prunus dulcis] [STFCJ]", "almond"),
    # --- "95% extraction" is a milling rate, not solvent extraction
    # ("white" is silent on the colour axis - it is the unmarked form of nearly
    #  every food that carries it, so printing it would split the plain rows
    #  from the ones whose source spelled the default out)
    ("Maize flour, white, 95% extraction [Zea mays L.] [PhyFoodComp]", "corn flour"),
    ("Sugar, granulated, white", "sugar"),
    # ... but on a BEAN or a CABBAGE white is a variety, not a refining state
    ("Cabbage, white, raw [Frida]", "cabbage, white"),
    ("Beans, white, mature seeds, raw", "bean, white"),
    # --- venue prefix must not swallow the dish (user-reported)
    ("McDONALD'S, Bacon Ranch Salad without chicken", "bacon ranch salad without chicken"),
    ("McDONALD'S, french fries", "french fry"),
    # --- an ingredient-difference parenthetical is identity, not annotation:
    #     these five pairs merged, and the omitted component carries the fat
    # NOTE: "Big Mac" and "Filet-O-Fish" are McDonald's TRADEMARKS, so the brand
    # policy replaces them with the generic dish. The pair still has to stay
    # SPLIT on the sauce - that is what these cases pin. Delete the matching
    # _CANON_OVERRIDES rows to keep the trademark names instead.
    ("McDONALD'S, BIG MAC", "double beef hamburger with cheese, lettuce and sauce"),
    ("McDONALD'S, BIG MAC (without Big Mac Sauce)",
     "double beef hamburger with cheese and lettuce, without sauce"),
    ("McDONALD'S, FILET-O-FISH (without tartar sauce)",
     "fish fillet sandwich with cheese, without tartar sauce"),
    ("McDONALD'S, Hotcakes (plain)", "hotcake"),
    # --- semicolons separate chunks exactly as commas do (user-reported)
    ("Tuna; chunk light; canned in water; drained solids", "tuna, chunk light"),
    # --- a token stripped late must not leave orphan punctuation (user-reported)
    ("Lobster, boiled/cooked in water [CIQUAL]", "lobster"),
    ("Quark, <", "quark"),
    ("Bamboo shoot, ]", "bamboo shoot"),
    ("Butter, butter", "butter"),
    # --- w = with, wo/w/o = without (NEVO, CIQUAL, SWE); 741 rows use "w/"
    ("Apple w/ skin [SWE]", "apple"),
    ("Apple w skin av [NEVO]", "apple"),
    # The w/o abbreviation must still expand to "without" rather than be lost
    # as a slash - that is what this case pins. Where the expansion LANDS
    # changed in the sugar-axis pass: "without sugar" is now the unsweetened
    # label, and unsweetened is unmarked, so the name is the bare food. The
    # split against the sweetened yoghurt survives on the axis instead.
    ("Yoghurt w/o sugar [SWE]", "yoghurt"),
    ("Bechamel sauce, w butter, homemade [CIQUAL]", "bechamel sauce, with butter"),
    # --- a W that is part of a code is NOT the abbreviation. Round 4 goes
    #     further and drops the breeding code itself, so the maize converges
    #     with every other maize instead of standing alone; the point the case
    #     pins is unchanged - the trailing W must never become "with".
    ("Maize, DMR-ESR-W, dried, raw [BioFoodComp]", "corn, dried"),
    # --- salt is an independent axis: these merged and sodium was max-unioned.
    #     UNSALTED is the unmarked default (a raw carrot is unsalted), so only
    #     "salted" prints - but the two must still land on DIFFERENT canons.
    ("Almonds w skin salted [NEVO]", "almond, salted"),
    ("Almonds w skin unsalted [NEVO]", "almond"),          # unsalted is unmarked
    ("Cod, dried, salted", "cod, dried, salted"),
    ("Cod, dried, unsalted", "cod, dried"),
    # --- "dry-roasted" is a roasting method, not dehydration (hyphen slipped past)
    ("Peanuts, all types, dry-roasted, with salt", "peanut, salted"),
    ("Peanuts, all types, dry-roasted, without salt", "peanut"),
    ("Oat Bran, dry", "oat bran, dried"),
    # --- full-width CJK brackets and the placeholder name field
    # The taxonomic head is stripped in the round-5 pass, so the species now
    # leads. What this case pins is unchanged: the asterisk and its bracketed
    # synonym must not survive into the name.
    ("fish, cod, walleye pollock*, raw \uff3b*syn. alaska pollock\uff3d",
     "cod, walleye pollock"),
    ("\uff0d [Saccharina japonica [Syn. Laminaria jaonica]] [STFCJ]", "saccharina japonica"),
    # --- percentage RANGES; the single-number rule left the orphan "2-"
    ("Cottage Cheese, 2-5% Fat, Lactose-Free [Fineli]", "cheese, cottage, lactose-free, 2-5% fat"),
    # --- WAFCT "local name *: english gloss"; the tail must not set the state
    ("babenda-1 *: sauce from green leaf, groundnut powder", "babenda-1, sauce from green leaf"),
    # --- compound conjunction _CONN_TRAIL_RE cannot see
    ("Bread, multigrain and/or with seeds [CIQUAL]", "multigrain bread"),
    # --- provenance, not identity
    ("Mild chicken strip, analyzed 2006", "mild chicken strip"),
    ("Fat, chicken, skin, braised, from drumsticks and thighs, non-enhanced (national) - 14b-03-01-totalfat",
     # the skin rows split off onto the part axis in round 5; chicken skin is
     # about 30 g fat/100 g and was being averaged into rendered chicken fat
     "chicken fat, skin"),
    # --- classifier chunks that name nothing
    ("Vermouth, dry type", "vermouth"),
    # --- ...but "general" IS an aggregate marker, and must survive
    ("Restaurant food, general", "restaurant food, general"),
    ("Butter, unknown fat content", "butter, unknown fat content"),
    # --- every informative token stripped: fall back, never emit a conjunction
    # the plural fold takes the LAST word, which is what keeps the two spellings
    # of this aggregate together
    ("Seeds and kernels av [NEVO]", "seeds and kernel"),
    # ... and neither "seeds" nor "kernels" is stripped when it is half of a
    # compound name: this row used to canonicalise to 'and peanut, dried'
    ("Seeds and peanuts, dried, unsalted (average) [CIQUAL]", "seeds and peanut, dried"),
    # FDC files its bean breeding accessions as a chunk of their own; 44 became
    # canons of one row each
    ("Beans, Dry, Pink, 11F-8082 (0% moisture)", "bean, pink, dried, 0% moisture basis"),
    ("Water, 0-50 mg calcium per litre", "water, 0-50 mg calcium per litre"),
    # --- MANUAL CURATION, reported by the project owner 2026-08-21
    ("Puddings, all flavors except chocolate, low calorie, instant, dry mix", "pudding, low calorie"),
    ("Puddings, all flavors except chocolate, low calorie, regular, dry mix", "pudding, low calorie"),
    ("Margarine-like spread, SMART BEAT Super Light without saturated fat", "margarine, light"),
    ("Margarine-like spread, SMART BEAT Smart Squeeze", "margarine"),
    ("Margarine-like spread, BENECOL Light Spread", "margarine, light"),
    ("Margarine-like spread, SMART BALANCE Light Buttery Spread", "margarine, light"),
    # ...but "spread" is the food itself here, and must survive
    ("Cheese spread, cream cheese", "cheese spread, cream cheese"),
    ("Sandwich spread, pork, beef", "sandwich spread, pork"),
    # --- a repeated spelling is not a qualifier
    ("Knackwurst, knockwurst, pork, beef", "knackwurst"),
    # --- merged at the owner's explicit request; botanically a different plant
    ("New Zealand spinach, raw", "spinach"),
    # ...and every variant of it, which a canon-keyed override could not reach
    ("New Zealand spinach, boiled, drained, with salt [CNF]", "spinach, salted"),
    ("Sausage, Knackwurst (Knockwurst), cooked [CNF]", "knackwurst"),
    # --- FDC writes both "Margarine-like spread," and bare "Margarine-like,"
    # Fat level became an axis of its own in the round-5 pass, so the claim
    # now survives into the name. It has to: a fat-free spread and a
    # margarine are not the same food, and merging them handed the fat-free
    # version 80 g of fat per 100 g.
    ("Margarine-like, vegetable oil spread, fat-free, tub",
     "margarine, vegetable oil, fat-free"),
    # --- analyte prefixes name the measurement, not the food
    ("FA - Beef, porterhouse steak, lean, raw", "beef, porterhouse, lean"),
    ("Carotenoids, American cheese - NFY090HK9", "american cheese"),
    ("Salt, Iodized", "salt, iodized"),
    # --- "dry" is a style on a drink, not dehydration
    ("Sherry, dry", "sherry"),
    # ...and the class noun goes, which also stops a dessert wine - a sweet,
    # high-sugar drink - being poured into the same canon as a table wine.
    ("Alcoholic beverage, wine, dessert, dry", "wine, dessert"),
    # --- the singulariser invented non-words and split foods in two
    ("Cookies, oatmeal", "cookie, oatmeal"),
    ("Chives, freeze-dried", "chive, dried"),
    # --- the canning medium is composition; oil- and water-packed tuna differ ~10x in fat
    ("Tuna in oil tinned [NEVO]", "tuna, in oil"),
    ("Tuna in water tinned [NEVO]", "tuna, in water"),
    # --- doubled head noun from the "class, item" convention
    ("Cake, cheese cake", "cheese cake"),
    ("Sausage, liver sausage", "liver sausage"),
    # a variety, not a repeat, so the doubled-head rule leaves it alone - and the
    # colour axis now keeps it apart from the pale beans it was folded in with
    ("Bean, black", "bean, black"),
    # --- colour is composition wherever a pigment is the nutrient
    ("Rice, black, unenriched, raw", "rice, black"),
    ("Rice, brown, long-grain, raw", "rice, brown"),
    ("Rice, white, long-grain, regular, enriched, cooked", "rice, enriched"),
    ("Common bean [Black], dehulled, raw [Phenol-Explorer]", "common bean, black"),
    ("Common bean [White], dehulled, raw [Phenol-Explorer]", "common bean, white"),
    ("Tea, black, brewed", "tea, brewed, black"),
    ("Tea, green, brewed", "tea, brewed, green"),
    ("Wine, table, red", "wine, red"),
    ("Wine, table, white", "wine, white"),
    ("Chicken, broilers or fryers, dark meat, meat only, raw",
     "chicken, dark meat, flesh"),
    # ... but only as a WHOLE chunk, so a cultivar name is never read as one
    ("Apples, raw, red delicious, with skin", "apple"),
    ("Potato tuber, Red LaSoda, flesh, raw [BioFoodComp]", "potato tuber, flesh"),
    # --- a plant part is composition: rice bran carries about 20 g of fat
    #     against milled rice's 0.7, and 54 bran rows sat inside 'rice'
    ("Rice, BETA, bran, raw [BioFoodComp]", "rice bran"),
    ("Bean, adzuki, sprouted, raw [Frida]", "bean, sprouted"),
    ("Cowpea, leaves, raw [BioFoodComp]", "cowpea, leaf"),
    # faba, fava and broad bean are Vicia faba under three names
    ("Faba bean, Vesuvio, hulls, raw [BioFoodComp]", "broad bean, hull"),
    # ... but for a leafy vegetable the leaf IS the food, and for a potato the
    # tuber is
    ("Lettuce, leaves, raw", "lettuce"),
    ("Potatoes, tuber, raw [Solanum tuberosum] [STFCJ]", "potato tuber"),
    ("Carrot, regular (European type), root, frozen [STFCJ]", "carrot"),
    # --- "without skin" IS the flesh, and it is a real difference: much of an
    #     apple's fibre is in the skin
    ("Apples, raw, without skin", "apple, flesh"),
    ("Apples, raw, with skin", "apple"),
    # --- the longissimus dorsi IS the loin muscle
    ("Pork, longissimus dorsi muscle, raw", "pork, loin"),
    # --- an aged cheddar has lost water, so it carries more of everything per
    #     100 g than a mild one
    ("Cheese, cheddar, sharp, sliced", "cheese, cheddar, mature"),
    ("Cheese, cheddar, mild, block/chunk", "cheese, cheddar, mild"),
    # ... read only on a cheese: elsewhere "mature" is a ripeness stage, and
    # "mature seeds" is FDC's phrase for the dry legume seed
    ("Carrots, mature, raw", "carrot"),
    ("Peach, mature, raw [BioFoodComp]", "peach"),
    ("Beans, black, mature seeds, raw", "bean, black"),
    ("Pork, spare ribs, in black bean sauce, meat only", "pork, spare rib, flesh"),
    # ... and bare "light" on a spread is the fat axis, not a colour
    ("Butter, light", "butter, light"),
    # --- NEVO and Livsmedelsverket write no punctuation; restoring the comma is
    #     enough for every rule downstream to work
    ("Beans broad raw [NEVO]", "broad bean"),
    ("Nuts macadamia unsalted [NEVO]", "macadamia nut"),
    ("Crackers cream [NEVO]", "cream cracker"),
    ("Nuts and seeds, mixed", "nut, mixed"),          # a connective is not a modifier
    ("Causses blue cheese [SWE]", "causses blue cheese"),   # not a plural
    ("Sports drink [NEVO]", "sports drink"),
    # --- "Chocolate, milk" is a milk chocolate BAR, not a drink
    ("Chocolate, milk [Phenol-Explorer]", "milk chocolate"),
    ("Chocolate, milk, beverage [Phenol-Explorer]", "chocolate milk, beverage"),
    ("Chocolate milk [NEVO]", "chocolate milk"),
    # --- the food keeps its own number when the material head folds onto it
    ("Lemons, juice, fresh [STFCJ]", "lemon juice"),
    ("Milk, goats, pasteurised [McCance]", "goat milk"),
    # --- "crisp" is an adjective; only the plural noun means the snack
    ("Pear, crisp pear, Suli, ripe, raw", "pear"),
    # --- dropping "without" reversed the meaning
    ("Peach/nectarine, Without Skin And Stone [Fineli]", "peach/nectarine"),
    # --- NEVO grades cheese "45+"; the bare number is not a name
    ("Cheese 45+ [NEVO]", "cheese"),
    # --- CROSS-DATABASE CONVERGENCE: same food, different national word order
    # Pressing grade is composition: extra virgin olive oil carries an order of
    # magnitude more polyphenols than refined oil from the same fruit, and
    # Phenol-Explorer files the three grades separately for that reason. Both
    # spellings of the head have to land on the same canon.
    ("Oil, olive, extra virgin", "olive oil, extra virgin"),
    ("Olive, oil, extra virgin [Phenol-Explorer]", "olive oil, extra virgin"),
    ("Olive, oil, refined [Phenol-Explorer]", "olive oil, refined"),
    ("Olive Oil [Fineli]", "olive oil"),
    # "refined" outside a fat context is a different word and must not be read
    # as a pressing grade
    ("Sugar, refined", "sugar, refined"),
    # material heads read modifier-first, and both source spellings converge
    ("Oil, coconut", "coconut oil"),
    ("Oil, palm", "palm oil"),
    ("Flour, rye", "rye flour"),
    ("Rye, flour", "rye flour"),
    ("Vinegar, balsamic", "balsamic vinegar"),
    ("Milk, cow", "cow milk"),
    # ... but a state or a grade in the tail slot never moves in front
    # "whole" on a dairy food is the fat level, and the house word for it
    ("Milk, whole", "milk, full-fat"),
    # "skimmed", "nonfat" and "fat-free" are three spellings of one food; the
    # label is what the canon carries, so the wording is dropped once read
    ("Milk, skimmed", "milk, fat-free"),
    ("Milk, nonfat", "milk, fat-free"),
    ("Cheese, low fat", "cheese, low-fat"),
    ("Flour, soy, defatted", "soy flour, fat-free"),
    ("Milk, skimmed, with added vitamin D", "milk, fat-free, enriched"),
    # the noun the marker qualifies goes WITH it, or it is left standing alone
    ("Cheese, mozzarella, part skim milk", "cheese, mozzarella, low-fat"),
    ("Milk dessert, frozen, milk-fat free, chocolate",
     "milk dessert, chocolate, fat-free"),
    # the earliest marker wins: the "nonfat" here belongs to the added solids
    ("Milk, lowfat, fluid, 1% milkfat, with added nonfat milk solids, "
     "vitamin A and vitamin D", "milk, fluid, low-fat"),
    # ... and part-skim still beats skim, because they start in the same place
    ("Mozzarella cheese, low-moisture, part-skim",
     "cheese, mozzarella, low-moisture, low-fat"),
    # unfortified is the unmarked form; FDC pairs the two claims in one string
    ("Milk, canned, evaporated, with added vitamin D and without added vitamin A",
     "milk, evaporated, enriched"),
    ("Milk, fluid, 1% fat, without added vitamin A and vitamin D", "milk, fluid, 1% fat"),
    # process is not composition
    ("Milk, semi-skimmed, pasteurised [CIQUAL]", "milk, low-fat"),
    # "powder" IS the dried state, and carrying both spellings split turmeric
    # across five canons
    ("Turmeric, powder", "turmeric, dried"),
    ("Turmeric, powdered, dried", "turmeric, dried"),
    ("Chili powder", "chili powder"),          # there the word is the name
    # the packing medium is often named, and the qualifier has to go with it
    ("Peach, canned in pear juice [AFCD]", "peach, in juice"),
    # which oil is kept: sardines in olive oil and in sunflower oil carry
    # different fatty acids, so folding them into "in oil" would hand each the
    # average of the other
    ("Tuna, canned in olive oil, drained", "tuna, in olive oil"),
    ("European pilchard or sardine, canned in sunflower oil [CIQUAL]",
     "sardine, in sunflower oil"),
    # light syrup is about half the added sugar of heavy syrup
    ("Fruit cocktail, canned, light syrup pack", "fruit cocktail, lightly sweetened"),
    ("Fruit cocktail, canned, heavy syrup pack", "fruit cocktail, sweetened"),
    # "without added fat" is only redundant once a fat level HAS been read
    ("Atlantic herring, without added fat, fried [BioFoodComp]",
     "herring, atlantic, without added fat, fried"),
    ("Peaches, canned, juice pack, solids and liquids", "peach, in juice"),
    ("Fruit cocktail, canned in heavy syrup", "fruit cocktail, sweetened"),
    # ... but "prepared with water" is reconstitution, not a medium
    ("Soup, vegetable chicken, canned, prepared with water, low sodium",
     "soup, chicken vegetable, low-salt"),
    ("Porridge, made with water [McCance]", "porridge, made with water"),
    # a reduced-sugar claim is a level on the sugar axis
    ("Milk, chocolate, lowfat, reduced sugar", "chocolate milk, low-sugar, low-fat"),
    # FDC couples the two claims in one phrase
    ("Popcorn, microwave, low fat and sodium", "microwave popcorn, low-fat"),
    ("Tortilla chips, low fat, baked without fat", "tortilla chip, low-fat"),
    # CNF spells it yogourt
    ("Yogourt, plain [CNF]", "yoghurt"),
    # a material head folds onto the food, but not behind a DISH, where the
    # chunk is an ingredient or a cooking medium
    ("Olive, oil, extra virgin [Phenol-Explorer]", "olive oil, extra virgin"),
    ("Rice Porridge, Milk, Salt [Fineli]", "rice porridge, milk"),
    ("Beans, butter, canned, re-heated, drained [McCance]", "bean, butter"),
    ("Cookies, butter, commercially prepared, unenriched", "cookie, butter"),
    ("Peanut, butter [Phenol-Explorer]", "peanut butter"),
    ("Tortillas, ready-to-bake or -fry, flour, shelf stable", "tortilla, flour"),
    # Fineli repeats the material in the modifier
    # UK cornflour IS US cornstarch, so the two chunks are one word twice
    ("Flour, Cornstarch, Cornflour [Fineli]", "cornstarch"),
    ("Flour, all purpose", "flour, all purpose"),
    ("Juice, cocktail, cranberry", "juice, cocktail"),
    ("Oil, industrial, canola, high oleic", "oil, industrial"),
    # "Flour, cornflour" names the material twice - drop the head, do not move it
    ("Flour, cornflour", "cornstarch"),
    # species read genus-first, which is the opposite majority
    ("Coho salmon", "salmon, coho"),
    ("Atlantic herring", "herring, atlantic"),
    ("Rainbow trout, farmed", "trout, rainbow, farmed"),
    # ... but a dish, a connective or a preparation is not a species modifier
    ("Pizza with tuna", "pizza with tuna"),
    ("Shrimp or prawn", "shrimp"),
    ("Tuna salad", "tuna salad"),
    ("Norway lobster", "norway lobster"),
    ("Horse mackerel [STFCJ]", "horse mackerel"),
    # "nut" is a category chapter like fish and spice
    ("Nuts, almonds", "almond"),
    # the curated table spells the pecan with its noun ('pecan nut'); the head
    # strip and the table agree on ONE canon, which is what matters here
    ("Nuts, pecans, dry roasted, with salt added", "pecan nut, salted"),
    # ... except where it is an ingredient in a dish or an aggregate entry
    ("Nuts, formulated, wheat-based, all flavors except macadamia, without salt",
     "nut, formulated"),
    ("Nut, mushroom and rice roast, homemade [McCance]",
     "nut, mushroom and rice roast"),
    # Fineli writes the possessive with an acute accent, so the brand strip took
    # "Kellogg" and left the "\u00b4S" behind as a word of its own
    ("Kellogg\u00b4S Frost Rice Krispies, Sugar-Coated, Vitamins [Fineli]",
     "frost rice krispy, sugar-coated"),
    # STFCJ marks its synonym footnotes with the FULL-WIDTH asterisk
    ("Fish, cod, walleye pollock\uff0a [STFCJ]", "cod, walleye pollock"),
    ("Olive oil [Fineli]", "olive oil"),
    ("Oil olive [NEVO]", "olive oil"),
    ("Chicken, liver, raw", "chicken, liver"),
    # the organ axis files the part behind the animal, whichever way the source
    # wrote it - and an ox liver is a beef liver, a calf's a veal one
    ("Liver chicken [NEVO]", "chicken, liver"),
    ("Liver, ox, raw [McCance]", "beef, liver"),
    ("Kidney, pig, raw [McCance]", "pork, kidney"),
    ("Beef round, raw", "beef, round"),
    ("Pork shoulder, raw", "pork, shoulder"),
    # ...but word order carries the MEANING in these two, so they stay apart
    ("Chocolate milk [NEVO]", "chocolate milk"),
    ("Milk chocolate [Fineli]", "milk chocolate"),
    # ...and a slashed dual common name is a convention, not stray punctuation
    ("Peach/nectarine, raw", "peach/nectarine"),
    # ...and opposite comparison operators must never converge
    ("Dark chocolate cocoa < 70% [SWE]", "dark chocolate cocoa < 70%"),
    ("Dark chocolate cocoa >70% [SWE]", "dark chocolate cocoa >70%"),
    # --- fast-food chains name product LINES, not foods (user-reported)
    ("WENDY'S, Jr. Hamburger, with cheese", "hamburger, with cheese"),
    ("WENDY'S, CLASSIC DOUBLE, with cheese", "double hamburger, with cheese"),
    ("WENDY'S, CLASSIC SINGLE Hamburger, no cheese", "hamburger, without cheese"),
    ("WENDY'S, Double Stack, with cheese", "double hamburger, with cheese"),
    # --- spelling converged on the variant more canons already used
    ("Aubergine, raw [McCance]", "eggplant"),
    ("Courgette, boiled", "zucchini"),
    ("Groundnut, raw [WAFCT]", "peanut"),
    ("Garbanzo beans, canned", "chickpea"),
    # --- "Sugars," / "Sweets," are panel labels on LAB rows but real foods
    #     otherwise. Removing sugar from _PANEL_HEADS (to save "Salt, Iodized")
    #     regressed these; the sample code is the discriminator.
    ("Sugars, Cheese, swiss, slices (CA2, CO) - 18c-17-03-Sug", "cheese, swiss"),
    ("Sugars, Fried rice, Chinese restaurant (NY1) - CY0906E", "rice, fried"),
    ("Sugars, brown", "sugar, brown"),
    ("Minerals, Sugar, Granulated, White, Name brand (AL1, CA1, MI1) - NFY040XEG", "sugar"),
    # --- an INGREDIENT's state is not the dish's state: "sundried tomato" in a
    #     sandwich was making the sandwich a dried food (~25 canons)
    ("Ciabatta sandwich w/ mozzarella cheese sundried tomato lettuce [SWE]",
     "ciabatta sandwich with mozzarella cheese tomato lettuce"),
    ("Muesli bar with dried fruit [SWE]", "muesli bar with fruit"),
    ("Sandwich, baguette, dry sausage and butter, homemade [CIQUAL]", "sandwich, baguette"),
    # --- TRIM is an axis: "veal" is a hard strip head and had absorbed 281
    #     descriptions, pure fat sitting with lean only. "lean and fat" is the
    #     whole cut and stays unlabelled.
    ("Veal, all cuts, separable fat, cooked [AFCD]", "veal, separable fat"),
    # the cut is composition: beef chuck carries about four times the fat of eye
    # of round, and the breed strip-head was dropping it along with the breed
    ("Veal, Australian, rib, rib roast, separable lean only, raw", "veal, rib, lean"),
    ("Veal, Australian, rib, rib roast, separable lean and fat, raw", "veal, rib"),
    ("Beef, chuck, arm pot roast, separable lean and fat, choice, raw", "beef, chuck"),
    ("Beef, round, eye of round, separable lean only, raw", "beef, round, lean"),
    ("Beef, brisket, whole, separable lean and fat, raw", "beef, brisket"),
    ("Pork, fresh, loin, whole, separable lean and fat, raw", "pork, loin"),
    ("Chicken, broiler, breast, meat only, raw", "chicken, breast, flesh"),
    # ... but only in a MEAT context: outside one these are ordinary words
    ("Broccoli, stalks, raw", "broccoli, stalk"),
    ("Celery, ribs, raw", "celery, rib"),
    # --- an organ is not the muscle it was cut out of. Beef liver carries about
    #     9,000 ug RAE of vitamin A per 100 g against roughly none in muscle.
    ("Beef, variety meats and by-products, liver, cooked, braised", "beef, liver"),
    ("Beef, variety meats and by-products, brain, raw", "beef, brain"),
    ("Turkey, all classes, gizzard, cooked, simmered", "turkey, gizzard"),
    # ... but "Common bean, Kidney" is a kidney BEAN
    ("Common bean, Kidney, 10kGy irradiated, raw [Phaseolus vulgaris L]",
     "common bean"),
    ("Lamb, New Zealand, composite cuts, fat, cooked [CNF]", "lamb, separable fat"),
    # --- country of origin is not composition, but only on a MEAT head
    ("Lamb, New Zealand, composite cuts, lean and fat, cooked [CNF]", "lamb"),
    # --- a serving diameter is not composition, and the TOPPING is
    ('Fast Food, Pizza Chain, 14" pizza, pepperoni topping, regular crust',
     "pizza, pepperoni topping"),
    ('Fast Food, Pizza Chain, 14" pizza, cheese topping, thin crust',
     "pizza, cheese topping"),
    # --- a product-line number is not composition; a unit-bearing one is kept
    ("Cereals ready-to-eat, HEALTH VALLEY, FIBER 7 Flakes", "fiber flake"),
    ("9 oz house sirloin steak", "sirloin steak"),
    ("Aquavit, 40 % vol.", "aquavit, 40 % vol"),
    # --- unsalted is the unmarked default (user-reported)
    ("Carrots, boiled, drained, without salt", "carrot"),
    ("Carrots, raw", "carrot"),
    ("Pistachio nuts, dry roasted, without salt", "pistachio nut"),

    # --- round 4: the class sweep ------------------------------------------
    # A corporate entity chunk is a manufacturer, never a food.
    ("The COCA-COLA company, DASANI, water, bottled", "water"),
    # Provenance is not identity; the same cherry per 100 g either way.
    ("Sweet cherries, imported from the U.S.A., raw", "sweet cherry"),
    # CIQUAL marks superseded rows in the description itself.
    ("Doughnut filled with fruit filling (e.g. jam), prepacked-> ARCHIVE [CIQUAL]",
     "doughnut filled with fruit filling, prepacked"),
    # "n/a" at the head is an unidentified plant and is named as one; away from
    # the head it is bookkeeping and goes, or it inherits the qualifier slot
    # when the region code in front of it is dropped as an accession code.
    ("N/A, leaf, raw [BioFoodComp]", "unidentified plant, leaf"),
    ("Mission Figs, Dried, Region 1, n/a, No, Carotenoids  - NF",
     "mission fig, dried"),
    # Accession and breeding codes, which the round-3 regex did not reach.
    ("Corn, Texas 17W, raw [BioFoodComp]", "corn"),
    ("American yam bean root, IRNAS n° 11 [BioFoodComp]", "american yam bean root"),
    ("Bambara nut, KARI/BN/ BK- [BioFoodComp]", "bambara nut"),
    ("Corn, cultivar: Cuzco, oil-roasted and salted [STFCJ]",
     "corn kernel, oil-roasted and salted"),
    # Spelling only, and applied to the finished canon: earlier in the pipeline
    # it turns the phrase into FDC's category prefix and the prefix strip then
    # promotes the product line behind it ("cereal, ready to eat, On Track").
    ("Coffee, ready to drink, sweetened", "coffee, ready-to-drink, sweetened"),
    ("Cereal, ready to  eat, On Track, President's Choice [CNF]",
     "cereal, ready-to-eat"),

    # --- round 4 regression guards -----------------------------------------
    # A bare "brand" marker is deliberately NOT stripped by rule. Removing it
    # before the override lookup invalidates fourteen curated keys: these four
    # each resolved to a junk name ("soybean flour, arrowhead mill") until the
    # rule that ate the marker was taken back out.
    ("Soybean flour, Arrowhead Mills brand [PhyFoodComp]", "soybean flour"),
    ("Tempeh, Turtle Islan brand [PhyFoodComp]", "tempeh"),
    ("Bread, wholemeal, store/other brands", "whole wheat bread"),
    ("Cereals, farina, enriched, assorted brands including CREAM OF WHEAT, "
     "quick (1-3 minutes), dry", "farina, enriched, quick cooking, dried"),
    # A brand name belongs in _BRAND_NAMES only when it is the whole chunk.
    # Adding "enfamil" promoted the formula line behind it and shattered this
    # canon into fifteen; adding "flax plus" did the same to "cereal".
    ("Infant formula, MEAD JOHNSON, ENFAMIL, ENFACARE, ready-to-feed",
     "infant formula, ready-to-feed"),
    ("Cereal, ready to eat, Flax Plus Maple Pecan Crunch, Nature's Path [CNF]",
     "cereal"),
    # _CORP_ENTITY_RE strips the entity chunk, so the five override keys that
    # carried one need their post-rule form too; this is the one that reaches
    # a curated target rather than a rule-built name.
    ("Beverages, THE COCA-COLA COMPANY, NOS energy drink, original",
     "energy drink, with sugar"),

    # --- round 5: the axes the sweep was missing ----------------------------
    # FDC writes the salt axis as "sodium added"; the phrase was ending up in
    # the name instead of on the axis.
    ("Blackeye pea, canned, sodium added, drained and rinsed", "blackeye pea, salted"),
    ("beans, great northern, canned, sodium added", "bean, great northern, salted"),
    # Negative forms have to beat the positive ones they contain.
    ("Beans, great northern, canned, solids and liquid, no salt added [CNF]",
     "bean, great northern"),
    ("Cheese, cottage, lowfat, 1% milkfat, no sodium added", "cheese, cottage, low-fat"),
    # Reduced sodium is neither end of the range and needs its own label.
    ("Beans, great northern, canned, solids and liquid, reduced sodium [CNF]",
     "bean, great northern, low-salt"),
    # Sugar is an axis: these two were one canon with 222 members.
    ("Pears, average, stewed with sugar [McCance]", "pear, sweetened"),
    ("Pears, average, stewed without sugar [McCance]", "pear"),
    # FDC writes the canning medium as "heavy syrup pack", with no "in".
    ("Peaches, canned, heavy syrup pack, solids and liquids", "peach, sweetened"),
    # Water pack is a MEDIUM, and ", in water" is the label 104 other canons
    # already carry ("tuna, in water", "bamboo shoot, in water"). Canned peaches
    # in water are not raw peaches - they are blanched and leached - so the two
    # stay apart. This case expected bare 'peach' while only the syrup form was
    # read; the water form fell through and the wording was left as a chunk.
    ("Peaches, canned, water pack, solids and liquids", "peach, in water"),
    # ... and once the medium has been poured off there is no medium to record,
    # but the WORDING still has to go or it becomes the food
    ("Pineapple, canned, juice pack, drained", "pineapple"),
    # ... but a syrup that IS the food keeps its name.
    ("Maple syrup", "maple syrup"),
    # Storage form is not identity: a chilled product and its ambient twin are
    # the same food per 100 g, so a whole chunk of "refrigerated" is dropped.
    ("Soymilk, refrigerated, unsweetened", "soymilk"),
    ("Soymilk, refrigerated, sweetened", "soymilk, sweetened"),
    # _STORAGE_CHUNK_RE deliberately does NOT reach inside a larger chunk, so it
    # leaves "refrigerated dough" - raw dough is not a baked biscuit - alone.
    # The wording is still lost here, but to the preparation strip rather than
    # to this rule; both canons hold only the dough row, so nothing is merged.
    # Filed in AUDIT_REMAINING_ROUND4.tsv rather than fixed.
    ("Biscuit, multigrain, refrigerated dough [CNF]", "biscuit, multigrain"),
    # Fat level is an axis; full fat is the unmarked case.
    ("Yogurt, plain, skim milk [CNF]", "yoghurt, fat-free"),
    # Draining pours off water, not absorbed oil.
    ("Fish, tuna, light, canned in oil, drained solids", "tuna, light, in oil"),
    ("Fish, tuna, light, canned in water, drained solids", "tuna, light"),
    # Enrichment is stripped as a preparation word, so it needs its own axis -
    # and the negated form must not fire it.
    ("Rice, white, long-grain, regular, enriched, cooked", "rice, enriched"),
    ("Almond drink, plain, no added sugars, not fortified, prepacked [CIQUAL]",
     "almond drink"),
    # A slogan is not a food: "KRAFT 100%" left the figure standing alone.
    # grated is a FORM, so it goes; what mattered here was that the brand and
    # the "100%" that survived its removal both go too
    ("Parmesan cheese, grated, KRAFT 100% (AL,CA1) - NFY120DQP", "cheese, parmesan"),
    # Garbanzo is chickpea; filing it under "bean" stranded it from the family.
    ("beans, chickpeas/garbanzos, dry", "chickpea, dried"),
    ("beans, chick peas/garbanzos, canned, sodium added", "chickpea, salted"),

    # Dry-matter rows are kept but labelled: their values are not per 100 g as
    # eaten, and 253 of the 261 members of 'bean, pinto, dried' are of this kind.
    ("Beans, Dry, Black (0% moisture)", "bean, black, dried, 0% moisture basis"),
    ("Beans, black, mature seeds, raw", "bean, black"),
    # A generic head in front of a legume that has its own family strands it.
    ("Beans, mung, mature seeds, raw", "mung bean"),
    # ... but a cranberry BEAN is borlotti, not the fruit, and a butter BEAN is
    # a lima, not dairy. Both are deliberately left under the generic head.
    ("Beans, cranberry (roman), mature seeds, raw", "bean, cranberry"),
    # part-skim is low-fat, not fat-free: bare "skim" must not match inside it.
    # The NF-suffix strip used to eat everything from the second chunk on, so
    # this row lost "low-moisture" - and 252 rows of branded cheddar lost the
    # word "cheddar" and landed in the bare 'cheese' canon.
    ("Cholesterol, Mozzarella cheese, low-moisture, part-skim, KRAFT (FL,MO) - NFY120DJH",
     "cheese, mozzarella, low-moisture, low-fat"),
    ("Minerals, Cheese, cheddar, natural shredded sharp, store brand, "
     "GREAT VALUE (CA1,NE) - NFY120WVO", "cheese, cheddar, mature"),
    # the grade is rarely a whole chunk - FDC puts it inside the brand here
    ("Minerals, Cheddar cheese, sliced, SARGENTO SHARP  (NY) - NFY120DBN",
     "cheese, cheddar, mature"),
    # twenty cheese types are written both ways; the head-first form wins
    # cow's milk is the default mozzarella; the buffalo one stays apart
    ("Mozzarella cheese, from cow's milk [CIQUAL]", "cheese, mozzarella"),
    ("Mozzarella cheese, from buffalo's milk [CIQUAL]",
     "cheese, mozzarella, from buffalo's milk"),
    ("Cottage cheese, full fat [Fineli]", "cheese, cottage, full-fat"),
    # one spelling per compound, whichever way the source hyphenates it
    ("Bread, gluten free [Frida]", "bread, gluten-free"),
    ("Barley flour, whole grain", "barley flour, wholegrain"),
    # ... but "cauliflower cheese" is a dish, not a cheese type
    ("Cauliflower cheese, retail [McCance]", "cauliflower cheese, retail"),
    # --- "Starch," heads an FDC lab row the way "Minerals," does, but AFCD
    #     files a real food as "Starch, potato". The sample code tells them
    #     apart, and the strip has to run while the head is still at the FRONT.
    ("Starch, potato [AFCD]", "potato starch"),
    ("Starch, Salsa, PACE CHUNKY, MEDIUM (AL,CA1) - NFY090KXF", "salsa"),
    ("Minerals, Sugar, Granulated, White, Name brand (AL1, CA1) - NFY040XEG", "sugar"),
    # --- on a coded row the ALL-CAPS chunks behind the food are the brand
    ("Minerals, Cheddar cheese, sliced, store brand, "
     "CRYSTAL FARMS & SHULLSBURG WISCONSIN (MO,NY) - NFY120DE1", "cheese, cheddar"),
    ("Minerals, Peanut butter, creamy, SKIPPY (MI) - NFY120C7W", "peanut butter"),
    # ... but off one an ALL-CAPS chunk can be the food itself
    ("BURGER KING - HAM", "ham"),
    ("OIL, OLIVE, EXTRA LIGHT", "olive oil, light"),
    # --- flavour is sugar: a strawberry Greek yoghurt carries about twice the
    #     sugar of the plain one, and FDC writes the flavour inside the brand
    ("Proximates, Greek yogurt, CHOBANI STRAWBERRY NON-FAT (NY1) - NFY0910ZZ",
     "greek yoghurt, strawberry, fat-free"),
    ("Minerals, Greek yogurt, CHOBANI PLAIN NON-FAT (NY1) - NFY0910ZW",
     "greek yoghurt, fat-free"),
    # ... read only where the plain form is the unmarked one, so these three
    # keep their own names
    ("Tomatoes, orange, raw", "tomato"),
    ("Melon, banana (Navajo)", "melon"),
    ("Plantain banana, raw [CIQUAL]", "plantain"),
    # "Selenium" is an FDC analyte panel head, not a food.
    ("Selenium, Yogurt, Greek, strawberry, non-fat, CHOBANI (CA2,NC) - NFY120OPI",
     "greek yoghurt, strawberry, fat-free"),
    # FDC puts the style after the comma where every other database fronts it;
    # unfronted, the flavour was pushed out of the two-chunk canon.
    ("Greek Yogurt, strawberry, non-fat", "greek yoghurt, strawberry, fat-free"),


    # ---- round 14 --------------------------------------------------------
    # "tinned" is "canned"; it also blocked the plural head fold
    ("Apricots in syrup tinned [NEVO]", "apricot, sweetened"),
    ("Mushrooms tinned [NEVO]", "mushroom"),
    # the USDA poultry class term names a market size, not a food - but only the
    # DEFAULT class goes, because a stewing hen is a different bird
    ("Chicken, broilers or fryers, wing, meat and skin, roasted", "chicken, wing"),
    ("Chicken, stewing, meat only, raw", "chicken, stewing, flesh"),
    ("Turkey, fryer-roasters, breast, meat only, raw", "turkey, fryer-roaster, breast, flesh"),
    # FDC's "meat only" is Fineli's "without skin"; "meat and skin" is the whole bird
    ("Chicken, Breast, Without Skin, Oven-Baked [Fineli]", "chicken, breast, flesh"),
    ("Chicken, broilers or fryers, breast, meat and skin, roasted", "chicken, breast"),
    # the fat level stated as a number. A RANGE is kept whole; 0% stays fat-free
    ("Cheese, Edam, 17% Fat [Fineli]", "cheese, edam, 17% fat"),
    ("Beef, ground, 70% lean meat / 30% fat, raw", "beef, ground, 30% fat"),
    ("Milk, 0 % Fat, Boiled [Fineli]", "milk, fat-free"),
    ("Turkey, breast, smoked, lemon pepper flavor, 97% fat-free", "turkey, smoked, breast, fat-free"),
    # one grind, four spellings, and ground SPICES are not it
    ("Beef, minced, raw", "beef, ground"),
    ("Pork, ground, raw", "pork, ground"),
    ("Chicken, mince [McCance]", "chicken, ground"),
    ("Ginger, ground", "ginger, ground"),
    ("Minced beef ball with egg and breadcrumb [SWE]", "minced beef ball with egg and breadcrumb"),
    # a source chunk that restates the head carries nothing - but an axis label
    # that happens to repeat a word of the head does
    ("Pea Stew, Dried Peas [Fineli]", "pea stew, dried"),
    ("Frozen Vegetable Mix, Grilled Vegetables, Salt, Oil [Fineli]", "vegetable mix, salt"),
    ("Water shield, young leaves, bottled in water [Brasenia schreberi] [STFCJ]",
     "water shield, leaf, in water"),
    # alternations: the same food under two names collapses, two foods do not
    ("Macaroni or noodles with cheese, microwaveable, unprepared", "macaroni with cheese"),
    ("Yoghurt or fermented milk, plain [CIQUAL]", "yoghurt"),
    ("Sausage, beef or pork meat, fried [SWE]", "sausage, beef or pork meat, fried"),
    ("Butter, stick or tub", "butter"),
    ("Ice cream bar or stick, chocolate coated [CIQUAL]", "ice cream bar or stick, chocolate coated"),
    # NEVO names its substitutes by what they replace and what they are made of
    ("Plant-based alternative to Gouda cheese based on coconut oil [NEVO]",
     "plant-based gouda cheese, coconut oil"),
    ("Soy milk, non-dairy alternative to milk, enriched [NEVO]", "soy milk, enriched"),
    # a DRAINED food loses the medium's label, and its wording with it
    ("Apricot, canned in pear juice, drained [AFCD]", "apricot"),
    ("Apricot, canned in pear juice [AFCD]", "apricot, in juice"),
    ("Tuna, canned in oil, drained solids", "tuna, in oil"),
    # the seed and the kernel are the food when something is made OF them
    ("Seed Bread [Fineli]", "seed bread"),
    ("Beans, black, mature seeds, raw", "bean, black"),
    ("Palm kernel oil, fortified with vitamin A, 600\u20131000 mcg/100g [WAFCT]",
     "palm kernel oil, 600\u20131000 mcg/100g, enriched"),
    # the claim, and the claim denied
    ("Cereal bar with fruit with vitamins and minerals [CIQUAL]", "cereal bar with fruit, enriched"),
    ("Salt, not fortified with iodine [CIQUAL]", "salt"),
    ("Tofu, with calcium sulfate, fried", "tofu, with calcium sulfate"),
    ("Cabbage with carrot [SWE]", "cabbage with carrot"),
    ("Lemonade fruit juice drink light with vitamin E and C [SWE]",
     "lemonade fruit juice drink light, enriched"),
    # "less salt" is the same claim as "reduced sodium"
    ("Rye Bread, Less Salt [Fineli]", "rye bread, low-salt"),
    ("Brown sauce, reduced salt/sugar [McCance]", "brown sauce, low-salt"),
    # Fineli files the trademark as the food's name
    ("Rye Bread, Reissumy, Vaasan [Fineli]", "rye bread"),
    ("Rufous Milkcap [Fineli]", "rufous milkcap"),
    ("Rainbow trout, farmed", "trout, rainbow, farmed"),
    # a curated entry has to survive a new axis label: the tables are consulted
    # again with the labels lifted off
    ("Cheese, Hard, 24% Fat [Fineli]", "cheese, hard, 24% fat"),
    ("Cream, Whipping, 38% Fat [Fineli]", "cream, whipping, 38% fat"),
    # Livsmedelsverket records how the item is SOLD
    ("Vegetarian sausage w/ soy protein chilled or frozen product [SWE]",
     "vegetarian sausage with soy protein"),
    # the material behind an animal is what it was cooked in, not what it is
    ("Chicken, broilers or fryers, meat and skin, cooked, fried, flour", "chicken, flour"),
    # CIQUAL's optional plural and its "and/or" nutrient list
    ("Chocolate flavoured milk beverage, with sugar(s), partially skimmed, fortified with vitamin(s) and/or mineral(s) [CIQUAL]",
     "chocolate flavoured milk beverage, sweetened, low-fat, enriched"),
    ("Milk, semi-skimmed [CIQUAL]", "milk, low-fat"),
    # CNF welds the claim to the cut; BioFoodComp welds it to the fillet
    ("Turkey, all classes, light meat only, raw [CNF]", "turkey, light meat, flesh"),
    ("Cod, wild, skinless fillet, raw [BioFoodComp]", "cod, wild, flesh"),
    ("Peach/nectarine, Without Skin And Stone [Fineli]", "peach/nectarine"),
    # a lone class term IS the bird
    ("Liver, broiler or fryer, raw [Frida]", "chicken, liver"),
    ("Beef, ground meat, raw [STFCJ]", "beef, ground"),
    # the class noun said twice
    ("Salad dressing, french dressing", "salad dressing, french"),
    ("Cake, cheese cake", "cheese cake"),
    ("Vanilla ice cream, without cream [SWE]", "vanilla ice cream, without cream"),
    # bread reads type-first, but only for the grain types
    ("Bread, rye", "rye bread"),
    ("Bread, bagel", "bread, bagel"),

    # --- round 15 -----------------------------------------------------------
    # "no added salt/sugar" was matched only from "added" onwards, so the
    # negative was left standing as a chunk: 37 canons ended in ", no"
    ("Peanut, no added salt [CIQUAL]", "peanut"),
    ("Muesli, No Added Sugar [Fineli]", "muesli"),
    ("Oat Macaroni, Boiled, No Added Salt [Fineli]", "oat macaroni"),
    # ...and the front-loaded POSITIVE was not read at all, so the salted rows
    # were merging into the unsalted canon
    ("Pasta, cooked, unenriched, with added salt", "pasta, salted"),
    ("Peanut Butter With Added Salt And Sugar [Fineli]",
     "peanut butter, salted, sweetened"),
    ("Potatoes, french fried, all types, salt added in processing, frozen, unprepared",
     "french fry, salted"),
    # the class noun goes when what follows names the drink
    ("Alcoholic beverage, beer, regular, all", "beer"),
    ("Alcoholic beverage, beer, light", "beer, light"),
    ("Alcoholic beverage, pina colada, canned", "pina colada"),
    # ...and stays when it does not: sake is not rice
    ("Alcoholic beverage, rice (sake)", "alcoholic beverage, rice"),
    # a qualifier that states the unmarked case says nothing
    ("Almond milk, unsweetened, plain, shelf stable", "almond milk"),
    ("Bagels, plain", "bagel"),
    ("Abalone, mixed species, raw", "abalone"),
    ("Wheat germ, regular [CNF]", "wheat germ"),
    ("Grenadier, from any fishing spot, raw [CIQUAL]", "grenadier"),
    ("Ice cream or sorbet or ice pop, any flavour [CIQUAL]",
     "ice cream or sorbet or ice pop"),
    # ...except on chocolate, where "plain" is British and Nordic for dark
    ("Chocolate, plain [McCance]", "chocolate, plain"),
    # ...and inside an alternation, where it is one of the two options
    ("Beans, baked, canned, plain or vegetarian", "bean, plain or vegetarian"),
    # an unnamed frying fat names nothing; a named one changes the fatty acids
    ("Fast foods, potato, french fried in vegetable oil", "french fry, fried"),
    ("Cod, in batter, fried in rapeseed oil [McCance]",
     "cod, in batter, fried, in rapeseed oil"),
    # wild or farmed, on the foods that come both ways
    ("Salmon, Atlantic, wild, raw", "salmon, atlantic, wild"),
    ("Salmon, Atlantic, farmed, raw", "salmon, atlantic, farmed"),
    ("Salmon, atlantic, aquaculture, raw [Frida]", "salmon, atlantic, farmed"),
    ("Fish, halibut, Greenland, wild caught, raw", "greenland halibut, wild"),
    # ...and not on the ones that are only ever taken from the wild
    ("Haddock, wild, dorsal muscle, raw [BioFoodComp]", "haddock, wild"),
    ("Cereal, hot, oats, instant: Wild Berry Medley, dry, Quaker [CNF]",
     "cereal, hot, dried"),
    # "wild" next to rice is the GRAIN, not a provenance
    ("Rice, wild, raw [CIQUAL]", "wild rice"),
    ("Rice, wild, WR-1, brown, raw [BioFoodComp]", "wild rice, brown"),
    # one part, four spellings
    ("Chicken, flesh only, raw [McCance]", "chicken, flesh"),
    ("Coconut, pulp, raw [BioFoodComp]", "coconut, flesh"),
    ("Jujube, fruit flesh, raw [BioFoodComp]", "jujube, flesh"),
    ("Chicken, meat, raw [CIQUAL]", "chicken, flesh"),
    # ...but a bare "meat" is an INGREDIENT on a dish, and the whole bird when
    # the source has said the skin and fat are in as well
    ("Babyfood, meat, beef, junior", "babyfood, meat"),
    ("Duck, raw, meat, fat and skin [McCance]", "duck"),
    # "whole" is the fat level on a dairy food and a form everywhere else
    ("Cow milk, whole [CIQUAL]", "cow milk, full-fat"),
    ("Greek yoghurt, whole [McCance]", "greek yoghurt, full-fat"),
    ("Almonds, whole", "almond, whole"),
    # ...and "whole milk" in a recipe is the milk, not the food's own fat level
    ("Sweet yeast dough for pie, whole milk [Fineli]",
     "sweet yeast dough for pie, whole milk"),
    ("Condensed whole milk, sweetened [NEVO]", "milk, condensed, sweetened"),
    # CNF's traditional-foods marker is not part of the food's name
    ("Game meat, native, caribou (reindeer), liver, raw [CNF]", "caribou, liver"),
    ("Fish, burbot (loche), native, raw [CNF]", "burbot"),
    # ...except where it is half of a species name
    ("Persimmons, native, raw", "persimmon, native"),
    # one claim, two spellings
    ("Mayonnaise, home-made [CIQUAL]", "mayonnaise, homemade"),
    ("Alcoholic beverage, wine, table, red, Merlot", "wine, red"),
    # STFCJ welds the skin claim to the word it qualifies, and taking only the
    # "with skin" half left "meat" standing on nine canons. The two must not
    # meet: a skin-on chicken breast carries several times the fat.
    ("Chicken, broiler, breast, meat with skin, raw [STFCJ]", "chicken, breast"),
    ("Chicken, broiler, breast, meat without skin [STFCJ]", "chicken, breast, flesh"),
    ("Bear, black, meat (Alaska Native)", "bear, black, flesh"),
    # ...and on a PLANT head "meat" is the substitute, not a part of an animal
    ("Soy, meat [SWE]", "soy meat substitute"),
    # Fineli coordinates the two claims and only the fat half was being read
    ("Beef Mince, 17% Fat, Fried Without Fat And Salt [Fineli]",
     "beef, ground, fried, 17% fat"),
    # the approximation word belongs to the figure, not to the food
    ("Tomme cheese, reduced fat, around 13% fat [CIQUAL]", "tomme cheese, 13% fat"),
    ("Cream, from Isigny, >= 35% fat [CIQUAL]", "cream, from isigny, 35% fat"),
    # Frida splits the figure from its unit across two chunks
    ("Milk, whole, 3.5, (UHT), % fat [Frida]", "milk, 3.5% fat"),
]
def main() -> int:
    failed = []
    for desc, want in CASES:
        got = canon(desc)
        if got != want:
            failed.append((desc, want, got))
    print(f"{len(CASES) - len(failed)}/{len(CASES)} canonicalize_food_name cases pass")
    for desc, want, got in failed:
        print(f"  FAIL {desc!r}\n       expected {want!r}\n       got      {got!r}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
