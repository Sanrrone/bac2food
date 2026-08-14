#!/usr/bin/env Rscript
# taxon_level_test.R — is family/genus the wrong rank for the Figure 2b tests?
#
# The manuscript reports PERMANOVA and ANOSIM at family and genus on the POOLED taxon x food
# matrix (one row per taxon, scores averaged over the infants it was detected in). Species
# cannot be tested there: one point per species means 54 groups of n=1, zero within-group
# variance, R^2 = 1 by construction. So this script tests all three ranks on the UNPOOLED
# matrix built by taxon_level_test.py, where a species detected in k samples has k replicate
# profiles.
#
# Two things make the unpooled design a different experiment, not a drop-in refinement:
#
#  * SPARSITY. An unpooled profile is one sample's truncated top-N list — a median of 10
#    foods out of 737 — where a pooled profile is the union over samples. Sorensen distance
#    between two ~10-food binary vectors is near 1 unless they overlap, so the distance
#    matrix is compressed against its ceiling. The mean within/between contrasts printed
#    below are there to show how much room is left, rather than let a p-value stand alone.
#  * NON-INDEPENDENCE. The same species recurs across samples and every taxon in one sample
#    is scored from the same gene set. Permutations are therefore BLOCKED by sample
#    (how(blocks = )), so labels are shuffled only within a sample and the test asks whether
#    a species profile is consistent ACROSS samples rather than whether samples differ.
#
# The nested model genus/species is the direct answer to "is genus mixing species with
# different enzyme repertoires": it partitions variance into a between-genus term and a
# species-within-genus term, sequentially. A large second term means genus is lumping.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(vegan); library(tibble)
})

here <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
figdir <- Sys.getenv("BAC2FOOD_FIGDIR", unset = here)
setwd(figdir)
set.seed(42)

NPERM <- as.integer(Sys.getenv("NPERM", unset = "999"))
cat(sprintf("vegan %s | %d permutations\n\n", packageVersion("vegan"), NPERM))

# ------------------------------------------------------------------ unpooled matrix
d <- read_csv("taxon_unpooled_long.csv", show_col_types = FALSE)
w <- d |> mutate(present = 1L) |>
  distinct(sample, taxon, genus, family, food, present) |>
  pivot_wider(names_from = food, values_from = present, values_fill = 0L)

meta <- w |> select(sample, taxon, genus, family)
mat  <- as.matrix(w[, -(1:4)])
cat(sprintf("unpooled: %d observations x %d foods; density %.3f\n",
            nrow(mat), ncol(mat), mean(mat > 0)))

bd <- vegdist(mat, "bray", binary = TRUE)
cat(sprintf("Sorensen distances: mean %.3f, median %.3f, %% at 1.0 (no shared food) = %.1f%%\n\n",
            mean(bd), median(bd), 100 * mean(bd == 1)))

# Mean within- vs between-group distance. A rank that carves the data usefully should show a
# gap here; if within == between the p-value is measuring something other than cohesion.
contrast <- function(dm, g) {
  m <- as.matrix(dm); same <- outer(g, g, "==")
  diag(same) <- NA
  ok <- !is.na(same) & upper.tri(same)
  c(within = mean(m[ok & same]), between = mean(m[ok & !same]))
}

# A rank is testable only where it has replication; singleton groups carry no within-group
# variance and would inflate the between-group term with points that were never at risk of
# being placed elsewhere.
test_rank <- function(rank_name, g, min_n = 3) {
  keep <- g %in% names(table(g))[table(g) >= min_n] & g != "Unclassified"
  if (sum(keep) < 6 || length(unique(g[keep])) < 2) {
    cat(sprintf("%-8s: too few groups with n>=%d to test\n", rank_name, min_n)); return(NULL)
  }
  sub <- as.dist(as.matrix(bd)[keep, keep])
  gg  <- factor(g[keep]); blk <- factor(meta$sample[keep])
  ctl <- how(nperm = NPERM, blocks = blk)
  set.seed(42); pm <- adonis2(sub ~ gg, permutations = ctl)
  set.seed(42); an <- anosim(sub, gg, permutations = ctl)
  set.seed(42); bp <- permutest(betadisper(sub, gg), permutations = NPERM)
  ct <- contrast(sub, g[keep])
  cat(sprintf(paste0("%-8s: R2=%.3f p=%.4f | ANOSIM R=%.3f p=%.4f | betadisper p=%.3f",
                     " | within %.3f vs between %.3f (gap %.3f) | n=%d obs, %d groups\n"),
              rank_name, pm$R2[1], pm$`Pr(>F)`[1], an$statistic, an$signif,
              bp$tab$`Pr(>F)`[1], ct[["within"]], ct[["between"]],
              ct[["between"]] - ct[["within"]], sum(keep), nlevels(gg)))
  invisible(list(R2 = pm$R2[1], p = pm$`Pr(>F)`[1], aR = an$statistic, ap = an$signif))
}

cat("--- unpooled, one rank at a time (blocked by sample) ---\n")
test_rank("family",  meta$family)
test_rank("genus",   meta$genus)
test_rank("species", meta$taxon)

# ------------------------------------------------------------------ nested: does species
# add anything beyond genus? Sequential terms, genus first, so the species term is the
# variance left over WITHIN genera.
cat("\n--- nested genus/species: what species explains that genus does not ---\n")
gk <- meta$genus %in% names(table(meta$genus))[table(meta$genus) >= 3] &
      meta$genus != "Unclassified"
# A genus represented by a single species contributes nothing to the nested term and would
# only dilute it, so the nested model is fitted on multi-species genera.
multi <- meta |> filter(gk) |> distinct(genus, taxon) |> count(genus) |> filter(n >= 2)
nk <- gk & meta$genus %in% multi$genus
sub <- as.dist(as.matrix(bd)[nk, nk])
gg <- factor(meta$genus[nk]); ss <- factor(meta$taxon[nk])
cat(sprintf("multi-species genera: %d (%s); %d observations, %d species\n",
            nrow(multi), paste(multi$genus, collapse = ", "), sum(nk), nlevels(ss)))
set.seed(42)
nest <- adonis2(sub ~ gg / ss, permutations = how(nperm = NPERM,
                                                  blocks = factor(meta$sample[nk])), by = "terms")
print(nest)

# The same comparison expressed as distances, which does not depend on a permutation model:
# how far apart are two profiles of the SAME species, two species of the same genus, and two
# different genera?
m <- as.matrix(sub); iu <- upper.tri(m)
same_sp  <- outer(ss, ss, "==")[iu]
same_gen <- outer(gg, gg, "==")[iu]
v <- m[iu]
cat(sprintf("\nmean Sorensen distance:\n  same species          %.3f (n=%s)\n",
            mean(v[same_sp]), format(sum(same_sp), big.mark = ",")))
cat(sprintf("  same genus, diff spp  %.3f (n=%s)\n",
            mean(v[same_gen & !same_sp]), format(sum(same_gen & !same_sp), big.mark = ",")))
cat(sprintf("  different genus       %.3f (n=%s)\n",
            mean(v[!same_gen]), format(sum(!same_gen), big.mark = ",")))

# ------------------------------------------------------------------ timepoint-pooled
# The middle design, and the one that answers the sparsity objection to the block above.
# Pooling every sample of a species gives one dense profile and no replication; pooling
# within TIMEPOINT gives three profiles per species, each a union over the infants sampled
# at 1, 6 and 12 months. Species stays testable and the distances come off the ceiling.
cat("\n--- timepoint-pooled: 3 replicates per species, denser profiles ---\n")
tp <- d |> mutate(month = sub(".*\\((\\d+) mo\\)", "\\1", sample)) |>
  distinct(month, taxon, genus, family, food) |> mutate(present = 1L) |>
  pivot_wider(names_from = food, values_from = present, values_fill = 0L)
tmeta <- tp |> select(month, taxon, genus, family)
tmat  <- as.matrix(tp[, -(1:4)])
tbd   <- vegdist(tmat, "bray", binary = TRUE)
cat(sprintf("timepoint-pooled: %d observations (%d taxa x %d timepoints) x %d foods; ",
            nrow(tmat), length(unique(tmeta$taxon)), length(unique(tmeta$month)), ncol(tmat)))
cat(sprintf("density %.3f; mean Sorensen %.3f\n", mean(tmat > 0), mean(tbd)))

test_tp <- function(rank_name, g, min_n = 3) {
  keep <- g %in% names(table(g))[table(g) >= min_n] & g != "Unclassified"
  if (sum(keep) < 6 || length(unique(g[keep])) < 2) {
    cat(sprintf("%-8s: too few groups with n>=%d to test\n", rank_name, min_n)); return(invisible())
  }
  sub <- as.dist(as.matrix(tbd)[keep, keep]); gg <- factor(g[keep])
  ctl <- how(nperm = NPERM, blocks = factor(tmeta$month[keep]))
  set.seed(42); pm <- adonis2(sub ~ gg, permutations = ctl)
  set.seed(42); an <- anosim(sub, gg, permutations = ctl)
  set.seed(42); bp <- permutest(betadisper(sub, gg), permutations = NPERM)
  ct <- contrast(sub, g[keep])
  cat(sprintf(paste0("%-8s: R2=%.3f p=%.4f | ANOSIM R=%.3f p=%.4f | betadisper p=%.3f",
                     " | within %.3f vs between %.3f (gap %.3f) | n=%d obs, %d groups\n"),
              rank_name, pm$R2[1], pm$`Pr(>F)`[1], an$statistic, an$signif,
              bp$tab$`Pr(>F)`[1], ct[["within"]], ct[["between"]],
              ct[["between"]] - ct[["within"]], sum(keep), nlevels(gg)))
}
test_tp("family",  tmeta$family)
test_tp("genus",   tmeta$genus)
test_tp("species", tmeta$taxon)

# ------------------------------------------------------------------ pooled matrix check
# The matrix the manuscript actually uses. Species is untestable here (n=1 per group); what
# CAN be asked is whether genus lumps dissimilar species, using the same distance contrast.
cat("\n--- pooled matrix (the one Figure 2b and the manuscript use) ---\n")
pl_long <- read_csv("fig2b_nmds_long.csv", show_col_types = FALSE)

# LINEAGE GUARD. These two CSVs are written by different scripts and nothing but the working
# directory ties them together, so a leftover fig2b_nmds_long.csv from an earlier predictor
# run sits happily beside a fresh taxon_unpooled_long.csv and gets compared to it. That
# already happened once: a pre-collision-fix pooled matrix (737 foods) was contrasted against
# post-fix unpooled data (736), and the mismatch showed up only as statistics that would not
# reconcile. Compare the food sets and stop rather than mix lineages.
f_pool <- sort(unique(pl_long$food)); f_unpool <- sort(unique(d$food))
if (!identical(f_pool, f_unpool)) {
  stop(sprintf(paste0("lineage mismatch: fig2b_nmds_long.csv has %d foods, ",
                      "taxon_unpooled_long.csv has %d (%d shared). They come from different ",
                      "predictor runs. Set BAC2FOOD_FIGDIR to the directory holding the ",
                      "current run and regenerate both."),
               length(f_pool), length(f_unpool), length(intersect(f_pool, f_unpool))))
}
pl <- pl_long |> select(taxon, family, food, mean_score) |>
  pivot_wider(names_from = food, values_from = mean_score, values_fill = 0)
pm_mat <- as.matrix(pl[, -(1:2)]); rownames(pm_mat) <- pl$taxon
keep <- rowSums(pm_mat > 0) >= 50           # same inclusion rule as make_figures.R
pm_mat <- pm_mat[keep, , drop = FALSE]; pl <- pl[keep, ]
pg <- sub(" .*$", "", pl$taxon)
pbd <- vegdist(pm_mat, "bray", binary = TRUE)
cat(sprintf("pooled: %d taxa x %d foods; density %.3f; mean Sorensen %.3f\n",
            nrow(pm_mat), ncol(pm_mat), mean(pm_mat > 0), mean(pbd)))
okg <- !grepl("^(uncultured|unclassified|candidatus)$", tolower(pg))
mg <- names(table(pg[okg]))[table(pg[okg]) >= 2]
sel <- okg & pg %in% mg
pc <- contrast(as.dist(as.matrix(pbd)[sel, sel]), pg[sel])
cat(sprintf("within-genus %.3f vs between-genus %.3f (gap %.3f) over %d taxa in %d multi-species genera\n",
            pc[["within"]], pc[["between"]], pc[["between"]] - pc[["within"]],
            sum(sel), length(mg)))
