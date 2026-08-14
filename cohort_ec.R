# cohort_ec.R — build the bacterium -> EC panel for a metagenomic cohort.
#
# In : gene_annot/<sample>.tsv     per-sample gene annotation (locus_tag ... EC_number ...)
#      gene_coords.tsv             locus_tag -> species map, from the HGT project
# Out: cohort_bac_ec.tsv              species / EC_number panel, ready for
#                                  4_predict/bac2food_predict.py --mag_tsv
#
# Note: gene_annot/ also holds <sample>_ec.tsv files (headerless: sample, taxid_species,
# locus, EC). Those are a DIFFERENT format that bac2food_predict.py reads directly, and
# they must not be globbed here — hence the _ec exclusion below.

library(tibble)
library(dplyr)
library(readr)
library(stringr)

ANNOT_DIR   <- "gene_annot"
GENE_COORDS <- "~/Desktop/horizontal_gene_transfer/gene_coords.tsv"
OUT_TSV     <- "cohort_bac_ec.tsv"

l <- list.files(path = ANNOT_DIR, pattern = "\\.tsv$", full.names = TRUE)
l <- l[!str_detect(basename(l), "_ec\\.tsv$")]

sample_ecdf <- lapply(l, function(tsv) {
  sname <- sub("\\.tsv$", "", basename(tsv))
  read_tsv(tsv, comment = "", quote = "", show_col_types = FALSE) %>%
    select(locus_tag, EC_number) %>%
    filter(!is.na(EC_number)) %>%
    filter(str_detect(EC_number, "^\\d+\\.\\d+\\.\\d+\\.\\d+$")) %>%
    mutate(sampleID = sname) %>%
    unique()
}) %>% bind_rows()

gcoords <- read_tsv(GENE_COORDS, show_col_types = FALSE)
gcoords$locus_tag <- sapply(gcoords$gID, function(g) strsplit(g, "__")[[1]][2])

cohort_ec_df <- merge(sample_ecdf, gcoords %>% select(sampleID, species, locus_tag),
                   by = c("sampleID", "locus_tag"), all.x = TRUE) %>% tibble()
cohort_ec_df <- cohort_ec_df %>% select(sampleID, species, EC_number) %>% unique()

write_tsv(cohort_ec_df, OUT_TSV)
