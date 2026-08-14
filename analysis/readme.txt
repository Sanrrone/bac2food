analysis/ — every analysis reported in the Data Descriptor.

These scripts are not part of the pipeline. They read the deposited exports and a directory
of predictor outputs, and they produce the figures and the numbers quoted in the article.
They are here because the paper states that no number it reports rests on an unscripted
analysis, and this folder is what makes that true.

Two inputs, both overridable by environment variable:

  BAC2FOOD_EXPORTS      the deposited TSVs        default /data/bac2food/exports
  BAC2FOOD_PREDICT_DIR  predictor outputs, one    default /data/bac2food/cohort_cov_phase11
                        set of tables per sample

The published run is cohort_cov_phase11. Earlier directories in the same style exist and
give different numbers; they are superseded.

Order matters only in one place: prep_figure_data.py writes the five CSVs that
make_figures.R reads, and it also holds anonymize(), which taxon_level_test.py imports.
Run it first.

  python3 prep_figure_data.py          # -> fig1a_sources.csv, fig1b_chebi.csv,
                                       #    fig2a_cohort_foods.csv, fig2b_nmds_long.csv,
                                       #    fig2c_differential.csv
  Rscript make_figures.R               # -> Figure1.{png,pdf}, Figure2.{png,pdf}

Sample labels never appear. anonymize() assigns "Infant N (M mo)" by sorting the predictor
output filenames on (timepoint, name), so the mapping is deterministic across reruns and
carries nothing from the original identifiers. Every script that touches per-sample data
goes through it rather than reimplementing it.

  file                       what it answers
  -------------------------  ------------------------------------------------------------
  prep_figure_data.py        figure inputs; anonymize(); family_of()
  make_figures.R             Figures 1 and 2 (needs R + ggplot2, patchwork, vegan)
  linkage_walk.py            the software-free linkage walk over the exports alone.
                             Imports no part of bac2food -- that is the point of it:
                             a reader can confirm the resource without the predictor.
  saturation_6_12.py         Usage Notes: why community capability saturates between
                             6 and 12 months while per-species demand does not
  rarefy_richness.py         the depth control on that contrast -- rarefies on LOCI, not
                             on species, since subsampling species would assume the answer
  core_stripping.R           Note S3: does a shared core conceal a taxonomic signal in the
                             ordination (progressive rare-food thresholds)
  taxon_level_test.py        Note S4: builds the two replicated rank designs
  taxon_level_test.R         Note S4: PERMANOVA / ANOSIM / dispersion over those designs
  taxon_level_effectsize.R   Note S4: R2 floors and group-count-free statistics

saturation_6_12.py needs scipy; nothing else here does. The pipeline's requirements.txt
does not pin scipy, because the pipeline itself never imports it.
