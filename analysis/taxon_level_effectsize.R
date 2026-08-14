#!/usr/bin/env Rscript
# taxon_level_effectsize.R — compare ranks on statistics that are not inflated by group count.
#
# Raw PERMANOVA R^2 cannot be compared across taxonomic ranks: it rises with the number of
# groups whether or not the grouping means anything. With g groups and n observations, a
# random factor already absorbs about (g-1)/(n-1) of the variance, so species (54 groups)
# starts from a much higher floor than family (16). Every rank comparison below is therefore
# reported as F, as R^2 against its own null floor, and as ANOSIM R and the within/between
# distance gap, both of which are free of the group-count effect.
#
# Runs only the cheap designs (timepoint-pooled, 147 rows; pooled, 53 rows). The unpooled
# design costs ~25 min at 999 permutations and is handled by taxon_level_test.R.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(vegan)
})

here <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
setwd(Sys.getenv("BAC2FOOD_FIGDIR", unset = here))
set.seed(42)
NPERM <- as.integer(Sys.getenv("NPERM", unset = "9999"))

d <- read_csv("taxon_unpooled_long.csv", show_col_types = FALSE)

tp <- d |> mutate(month = sub(".*\\((\\d+) mo\\)", "\\1", sample)) |>
  distinct(month, taxon, genus, family, food) |> mutate(present = 1L) |>
  pivot_wider(names_from = food, values_from = present, values_fill = 0L)
meta <- tp |> select(month, taxon, genus, family)
mat  <- as.matrix(tp[, -(1:4)])
bd   <- vegdist(mat, "bray", binary = TRUE)

cat(sprintf("timepoint-pooled: %d obs x %d foods, density %.3f, mean Sorensen %.3f\n\n",
            nrow(mat), ncol(mat), mean(mat > 0), mean(bd)))
cat(sprintf("%-8s %6s %7s %8s %8s %8s %9s %9s\n",
            "rank", "groups", "F", "R2", "R2null", "R2/null", "ANOSIM_R", "gap"))

for (r in c("family", "genus", "species")) {
  g <- if (r == "species") meta$taxon else if (r == "genus") meta$genus else meta$family
  keep <- g %in% names(table(g))[table(g) >= 3] & g != "Unclassified"
  sub <- as.dist(as.matrix(bd)[keep, keep]); gg <- factor(g[keep])
  ctl <- how(nperm = NPERM, blocks = factor(meta$month[keep]))
  set.seed(42); pm <- adonis2(sub ~ gg, permutations = ctl)
  set.seed(42); an <- anosim(sub, gg, permutations = ctl)
  m <- as.matrix(sub); same <- outer(gg, gg, "=="); iu <- upper.tri(m)
  gap <- mean(m[iu & !same]) - mean(m[iu & same])
  # The variance a grouping of this many levels absorbs by construction, given the sample
  # size. R2 below this line is not evidence of anything.
  r2null <- (nlevels(gg) - 1) / (sum(keep) - 1)
  cat(sprintf("%-8s %6d %7.2f %8.3f %8.3f %8.2f %9.3f %9.3f\n",
              r, nlevels(gg), pm$F[1], pm$R2[1], r2null, pm$R2[1] / r2null,
              an$statistic, gap))
}
