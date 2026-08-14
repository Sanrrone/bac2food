# core_stripping.R — does a shared core of common foods conceal a taxonomic signal?
#
# Supplementary Note S5 asserts that the family overlap in Figure 2b is not an artifact of
# every taxon sharing a common core: restricting the ordination to progressively rarer foods
# leaves the family effect flat while stress rises. Until now that check had NO script behind
# it — the numbers came from a side session, so a change to the predictor could invalidate
# them silently, which is exactly what happened when the differential ranking rule changed.
# This is the missing artifact.
#
# Design, matching the claim: keep foods carried by at most X of the taxa, for X in
# 90/75/50/25 percent, and at each threshold refit the same test Figure 2b uses — PERMANOVA
# and ANOSIM on Sorensen distances over families with at least three members — plus the NMDS
# stress, since the claim is that stress rises while the effect does not.
#
#   Rscript core_stripping.R          # reads fig2b_nmds_long.csv in the working directory
suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(vegan)})

set.seed(1)
NPERM <- 9999
THRESH <- c(0.90, 0.75, 0.50, 0.25)

long <- read_csv("fig2b_nmds_long.csv", show_col_types = FALSE)
fam_of <- long |> select(taxon, family) |> distinct()

wide <- long |> select(taxon, food, mean_score) |>
  pivot_wider(names_from = food, values_from = mean_score, values_fill = 0) |>
  as.data.frame()
rownames(wide) <- wide$taxon; wide$taxon <- NULL
n_taxa <- nrow(wide)

# prevalence = fraction of taxa carrying the food at all (presence, matching the Sorensen
# distance the ordination uses)
prev <- colSums(wide > 0) / n_taxa

cat(sprintf("matrix: %d taxa x %d foods; prevalence median %.3f, max %.3f\n\n",
            n_taxa, ncol(wide), median(prev), max(prev)))
cat(sprintf("%-10s %7s %7s %8s %8s %9s %8s\n",
            "keep<=", "foods", "taxa", "R2", "p", "ANOSIM_R", "stress"))
cat(strrep("-", 62), "\n")

for (th in THRESH) {
  keep <- names(prev)[prev <= th]
  sub <- wide[, keep, drop = FALSE]
  # a taxon left with no food cannot enter a Sorensen distance
  sub <- sub[rowSums(sub > 0) > 0, , drop = FALSE]
  fam <- fam_of$family[match(rownames(sub), fam_of$taxon)]
  # same rule as Figure 2b: real families with >= 3 members; singletons carry no
  # within-group variance and would inflate the between-group term
  ok <- !is.na(fam) & fam != "Unclassified"
  fn <- table(fam[ok])
  test <- ok & fam %in% names(fn)[fn >= 3]
  if (sum(test) < 5 || length(unique(fam[test])) < 2) {
    cat(sprintf("%-10.2f %7d %7d   too few taxa/families to test\n", th, ncol(sub), sum(test)))
    next
  }
  d <- vegdist(sub[test, , drop = FALSE], method = "bray", binary = TRUE)
  gg <- factor(fam[test])
  pm <- adonis2(d ~ gg, permutations = NPERM)
  an <- anosim(d, gg, permutations = NPERM)
  st <- tryCatch(metaMDS(d, k = 2, trymax = 50, trace = 0)$stress, error = function(e) NA_real_)
  cat(sprintf("%-10.2f %7d %7d %8.3f %8.4f %9.3f %8.3f\n",
              th, ncol(sub), sum(test), pm$R2[1], pm$`Pr(>F)`[1], an$statistic, st))
}
cat("\nClaim under test: the family effect stays flat while stress rises, i.e. no subset of\n")
cat("the foods separates the families better than the whole.\n")
