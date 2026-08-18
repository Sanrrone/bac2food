# Data licence

The deposited resource is four files: `food_nutrients.tsv`, `species_enzymes.tsv`,
`enzyme_substrate_chebi.tsv` and `licences.tsv`. They are not all under one term,
because the 16 source databases they are assembled from are not.

Rights are machine-readable. `licences.tsv` carries one row per `source_db` and
joins to the composition table on that column, so any partition below can be
selected with a single join.

## The four partitions

**1. Open — CC BY 4.0.**
Values derived from USDA FoodData Central, CIQUAL, McCance & Widdowson's CoFID,
the Canadian Nutrient File, Fineli, Livsmedelsdatabasen, the Swiss FCDB, STFCJ
(Japan), Frida, Phenol-Explorer 3.0 and the curated synthetic bacterial-substrate
set — together with **all curation original to this work**: the nutrient→ChEBI→EC
map, the harmonization itself, the licence table, and the bifunctional-EC
curation — are released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Attribute using the `attribution_string` column of `licences.tsv`.

**2. ShareAlike — CC BY-SA 3.0 AU.**
Values with `source_db = afcd` (Australian Food Composition Database, FSANZ)
retain [CC BY-SA 3.0 AU](https://creativecommons.org/licenses/by-sa/3.0/au/).
Derivatives that include these rows must themselves be ShareAlike. Commercial
use is permitted.

**3. NonCommercial + ShareAlike — CC BY-NC-SA 3.0 IGO.**
Values with `source_db` in `wafct`, `biofoodcomp` or `phyfoodcomp` (FAO/INFOODS)
retain [CC BY-NC-SA 3.0 IGO](https://creativecommons.org/licenses/by-nc-sa/3.0/igo/).
These rows may **not** be used commercially, by us or by anyone else.

**4. Provider permission — RIVM terms, not a public licence.**
Values with `source_db = nevo` are derived from NEVO-online and are here by written
permission from RIVM, granted 18 August 2026 for this resource. They are not
relicensed by us, and that permission does not travel with the file: **redistributing
these rows onward needs RIVM's own consent**, which is what separates this partition
from the three above. Any use must carry the reference

> NEVO-online version 2025/9.0, RIVM, Bilthoven, The Netherlands

and RIVM's conditions also forbid charging end users for the use of NEVO data. The
NEVO dataset itself is not redistributed here; request it free from RIVM at
<https://www.rivm.nl/en/dutch-food-composition-database/use-of-nevo-online/request-dataset>.
Three changes were made to the source values and no others — identifiers mapped onto
the shared keys with NEVO's own code kept in `source_food_code`, amounts expressed per
100 g edible portion, and notations made machine-readable — and no NEVO value is
combined with any other database's.

A value reported identically by two sources carries a `;`-joined `source_db`
label, has its own row in `licences.tsv`, and inherits the **strictest** terms of
its constituents.

## If you need uniformly permissive data

Filter to `tier = Open` in `licences.tsv` and join. What remains is CC BY 4.0
throughout and is the great majority of the table.

## The microbial reference layers

`species_enzymes.tsv` derives from eggNOG v7 and `enzyme_substrate_chebi.tsv`
from BRENDA (CC BY 4.0). Both sit outside the food-source tiering and are
redistributed with the attribution their licences require.

## Why there is no NonCommercial option for the whole resource

A ShareAlike source (AFCD) cannot be relicensed as NonCommercial, and a
NonCommercial source (FAO) cannot be relicensed as permissive or sold. NEVO is not
ours to relicense in either direction. No single term covers the table, so
partitioning is not a stylistic choice here; it is the only lawful arrangement.
