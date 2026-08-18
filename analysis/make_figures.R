#!/usr/bin/env Rscript
# make_figures.R — the two manuscript figures.
#
# Fig 1  the resource: source contributions and substrate-to-ChEBI mapping confidence.
#        (Conflict counts and the linkage funnel are reported in the text, not plotted.)
# Fig 2  bac2food applied to an infant gut-metagenome cohort. Sample labels are already
#        anonymized upstream in prep_figure_data.py; no cohort identifier reaches this file.
#
# Multi-panel figures are emitted as ONE combined file with embedded a/b/c tags, which is
# what Scientific Data requires (panels may not be uploaded separately).
#
# TYPEFACE: deliberately none. These use ggplot2's default device font rather than a
# Helvetica/Arial clone -- the plain default is what is wanted here. Do not reintroduce a
# `base_family=` / `family=` argument.
#
# ORDINATION NOTE: the NMDS runs on ONE pooled taxon x food matrix, presence/absence,
# Sorensen (= binary Bray-Curtis). Two choices are load-bearing:
#  * POOLED, not stratified by infant. Per infant a taxon carries only ~20 foods from a
#    truncated top-N list (95% sparse) and metaMDS collapses to a degenerate near-zero
#    stress solution with an "insufficient data" warning. The union over infants gives
#    54 taxa x 257 foods and a genuine fit (stress 0.23 in 2-D, 0.17 in 3-D, converged).
#  * BINARY. Rank depth is comparable between taxa; raw score magnitude inside a truncated
#    list is not. The abundance-weighted version also fails to converge here.
# Do not "restore" per-infant faceting or abundance weighting: both look finer-grained and
# are actually meaningless.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(forcats); library(scales); library(patchwork); library(vegan); library(tibble)
})

here <- dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])))
# Read the CSVs from, and write the figures into, BAC2FOOD_FIGDIR when it is set. Matches
# the same override in prep_figure_data.py so the public and fuller builds are
# drawn by one script rather than by two copies that can drift apart.
figdir <- Sys.getenv("BAC2FOOD_FIGDIR", unset = here)
setwd(figdir)
set.seed(42)

# Okabe-Ito derived: colourblind-safe and legible in greyscale.
OK <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73", red = "#D55E00",
        purple = "#CC79A7", sky = "#56B4E9", yellow = "#F0E442", grey = "#8C8C8C",
        brown = "#8B5A2B", teal = "#117777")

base <- theme_minimal(base_size = 7.4) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = .22, colour = "grey88"),
        axis.title   = element_text(size = 7),
        axis.text    = element_text(size = 6.6, colour = "grey20"),
        # Panel titles NAME what is plotted; they do not state the finding. Every
        # qualifier that used to sit in a subtitle (axis transform, distance metric,
        # what the point labels mean, free facet scales) is in the manuscript caption,
        # which is where a reader of the printed figure looks for it. Do not
        # reintroduce plot.subtitle: it duplicates the caption and eats panel height.
        plot.title    = element_text(size = 7.8, face = "plain", colour = "black",
                                     margin = margin(b = 3)),
        plot.title.position = "plot",
        legend.key.size = unit(2.9, "mm"),
        legend.text  = element_text(size = 6.2),
        legend.title = element_blank(),
        plot.margin  = margin(3, 5, 3, 3))

lbl <- function(...) geom_text(..., size = 1.95, colour = "grey30", hjust = -0.12)

# ---------------------------------------------------------------- Figure 1
s <- read_csv("fig1a_sources.csv", show_col_types = FALSE) |>
  mutate(source = fct_reorder(source, rows),
         tier = factor(recode(tier,
                  "Open (redistributable)"   = "Open",
                  "Copyleft / NonCommercial" = "Copyleft / NC",
                  "Provider permission"      = "By permission"),
                levels = c("Open", "Copyleft / NC", "By permission")))

p1a <- ggplot(s, aes(rows, source, fill = tier)) +
  geom_col(width = .70) +
  lbl(aes(label = comma(rows))) +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                breaks = c(1e2, 1e4, 1e6), expand = expansion(mult = c(0, .40))) +
  scale_fill_manual(values = unname(OK[c("blue", "orange", "purple")])) +
  labs(title = "Food–nutrient values by source database",
       x = "food–nutrient values", y = NULL) +
  base + theme(legend.position = "bottom", legend.margin = margin(t = -6, b = -2))

ch <- read_csv("fig1b_chebi.csv", show_col_types = FALSE) |>
  mutate(tier = factor(tier, levels = rev(tier)))

p1b <- ggplot(ch, aes(rows, tier, fill = tier)) +
  geom_col(width = .68, show.legend = FALSE) +
  lbl(aes(label = paste0(comma(rows), " (", percent(rows / sum(rows), accuracy = .1), ")"))) +
  # .52 was not enough headroom for the longest label: "101,617 (46.5%)" on the widest bar
  # had its closing bracket clipped at the panel edge. The label is bar length + text, so
  # the widest bar needs room for the widest string, not an average one.
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale()),
                     expand = expansion(mult = c(0, .62))) +
  scale_fill_manual(values = rev(unname(OK[c("green", "blue", "sky", "grey", "yellow", "orange")]))) +
  # The asterisk in the "No ChEBI entity exists*" tier label is expanded in the Figure 1
  # caption ("macromolecules, peptides, assay markers and metal clusters"). It is not
  # expanded here: the full wording overruns the 90 mm half-panel and is clipped.
  labs(title = "Substrate-to-ChEBI mapping confidence",
       x = "enzyme–substrate rows", y = NULL) +
  base

fig1 <- (p1a | p1b) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave("Figure1.pdf", fig1, width = 180, height = 78, units = "mm", device = cairo_pdf)
ggsave("Figure1.png", fig1, width = 180, height = 78, units = "mm", dpi = 400)

# ---------------------------------------------------------------- Figure 2
co <- read_csv("fig2a_cohort_foods.csv", show_col_types = FALSE) |>
  mutate(food = fct_reorder(food, mean_score),
         category = sub(" and .*", "", category))

p2a <- ggplot(co, aes(mean_score, food)) +
  geom_segment(aes(x = 0, xend = mean_score, y = food, yend = food),
               colour = "grey85", linewidth = .35) +
  geom_point(aes(colour = category), size = 1.6) +
  geom_text(aes(label = paste0(n_infants, "/", n_cohort)), hjust = -0.7, size = 1.85,
            colour = "grey40") +
  # .16 was sized for a 4-infant cohort ("4/4"); "55/55" is twice as wide and the top
  # row's label ran off the panel. Keep headroom for the widest denominator.
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .24))) +
  # Named, not positional. A positional palette silently remaps colours whenever the
  # ranked set changes category composition, and errors outright when a new category
  # appears -- which is what coverage weighting did by promoting Spices into the top 12.
  scale_colour_manual(values = c("Beverages" = OK[["orange"]], "Cereal Grains" = OK[["sky"]],
                                 "Fruits" = OK[["green"]], "Legumes" = OK[["purple"]],
                                 "Vegetables" = OK[["blue"]], "Spices" = OK[["red"]],
                                 "Nut and Seed Products" = OK[["brown"]],
                                 "Dairy" = OK[["teal"]], "Baked Products" = OK[["grey"]]),
                      na.value = OK[["grey"]], drop = TRUE) +
  labs(title = "Community-level food ranking",
       x = "community suitability score", y = NULL) +
  base + theme(legend.position = "bottom", legend.margin = margin(t = -6, b = -2)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.8)))

# --- NMDS: ONE ordination, every taxon a point, coloured by taxonomic family
d <- read_csv("fig2b_nmds_long.csv", show_col_types = FALSE)
w <- d |> select(taxon, family, food, mean_score) |>
  pivot_wider(names_from = food, values_from = mean_score, values_fill = 0)
mat <- as.matrix(w[, -(1:2)]); rownames(mat) <- w$taxon

# Inclusion criterion, applied as a rule rather than to a named point: a taxon needs a
# profile dense enough to place. Sorensen distance from a near-empty presence/absence vector
# is close to 1 against everything, so an under-observed taxon lands far from all others for
# want of data rather than for biology.
#
# The cut is a FRACTION of the food set, not a count. It was written as a literal 50 when the
# matrix held 737 foods (6.8% of them), which was correct then and silently wrong later: the
# differential ranking rule reduced the matrix to 434 foods, where a literal 50 is a 11.5%
# bar. It then dropped 7 taxa instead of 1 — including three of the four Enterococcaceae,
# taking that family below the n>=3 rule and out of the family test altogether. A criterion
# meant to remove unplaceable points had begun removing a clade.
#
# Holding the fraction fixed reproduces the original intent on any matrix, and the assertion
# below enforces the property the original comment claimed but never checked: that the cut
# sits in a gap in the foods-per-taxon distribution rather than slicing through a dense part
# of it, so the excluded points are genuinely isolated rather than merely the tail.
KEEP_FRAC <- 0.068
n_food_per_taxon <- rowSums(mat > 0)
cut_at <- max(2, floor(KEEP_FRAC * ncol(mat)))
keep <- n_food_per_taxon >= cut_at
srt <- sort(n_food_per_taxon)
below <- srt[srt < cut_at]; above <- srt[srt >= cut_at]
if (length(below) && length(above)) {
  gap <- min(above) - max(below)
  cat(sprintf("  [nmds] cut at %d foods (%.1f%% of %d); gap below/above = %d/%d (width %d)\n",
              cut_at, 100 * KEEP_FRAC, ncol(mat), max(below), min(above), gap))
  if (gap < 5)
    warning(sprintf(paste0("sparse-taxon cut at %d slices a dense part of the distribution ",
                           "(gap only %d); the excluded taxa are a tail, not isolated points. ",
                           "Re-derive the threshold before trusting the family test."),
                    cut_at, gap))
}
if (any(!keep)) cat(sprintf("  [nmds] excluded %d sparse taxa: %s\n", sum(!keep),
                            paste(sprintf("%s (%d foods)", rownames(mat)[!keep],
                                          rowSums(mat > 0)[!keep]), collapse = ", ")))
mat <- mat[keep, , drop = FALSE]; w <- w[keep, ]

# k=3: in 2-D the retained taxa alone fit at stress 0.235, above the conventional 0.2 ceiling,
# because dropping the far point leaves only the tight core to satisfy. Three dimensions bring
# it back to 0.160. Panels show axes 1 and 2; the distance matrix, and therefore every
# statistic below, is unaffected by the choice of k.
nm  <- metaMDS(mat, distance = "bray", binary = TRUE, k = 3, trymax = 200,
               trace = 0, autotransform = FALSE)

# Does family explain the ordination? Test it rather than asserting it from the eye.
# Singleton families carry no within-group variance, so the test uses families with n>=3.
# "Unclassified" is excluded: it is a bin of taxa whose genus is unresolved, not a clade,
# and treating it as one group would fabricate within-group cohesion that does not exist.
fam_n  <- table(w$family[w$family != "Unclassified"])
testab <- w$family %in% names(fam_n)[fam_n >= 3]
bd     <- vegdist(mat, "bray", binary = TRUE)
bd_sub <- as.dist(as.matrix(bd)[testab, testab])
set.seed(42)
pm  <- adonis2(bd_sub ~ w$family[testab], permutations = 9999)
r2  <- pm$R2[1]; pv <- pm$`Pr(>F)`[1]; fstat <- pm$F[1]
# R2 rises with group count whether or not the grouping means anything, so report
# the null floor (g-1)/(n-1) beside it; the ratio is the comparable quantity.
r2null <- (sum(fam_n >= 3) - 1) / (sum(testab) - 1)
# A significant PERMANOVA can come from unequal WITHIN-group spread rather than separated
# centroids, so report the dispersion test alongside it. Here it is far from significant
# (F = 0.42, P = 0.91), which is what licenses reading the result as real separation.
set.seed(42)
bdp <- permutest(betadisper(bd_sub, factor(w$family[testab])), permutations = 9999)
# ANOSIM is quoted alongside PERMANOVA in the text; compute it here rather than in a
# side session, so every statistic the manuscript reports comes out of this one script.
set.seed(42)
an <- anosim(bd_sub, factor(w$family[testab]), permutations = 9999)
cat(sprintf("PERMANOVA family: F=%.2f R2=%.3f (null %.3f, ratio %.2f) p=%.4f (n=%d taxa, %d families) | betadisper p=%.3f | ANOSIM R=%.3f p=%.3f\n",
            fstat, r2, r2null, r2/r2null, pv, sum(testab), sum(fam_n >= 3), bdp$tab$`Pr(>F)`[1],
            an$statistic, an$signif))

# The same test one rank finer. The manuscript quotes these alongside the family result, and
# they previously came from a side session with no script behind them — the reason they are
# computed here now. Genus is the first token of the taxon name; "uncultured"/"unclassified"
# placeholders are not genera and are dropped rather than pooled into a fake clade.
gen      <- sub(" .*$", "", rownames(mat))
gen_ok   <- !grepl("^(uncultured|unclassified|candidatus)$", tolower(gen))
gen_n    <- table(gen[gen_ok])
gtest    <- gen_ok & gen %in% names(gen_n)[gen_n >= 3]
if (sum(gtest) >= 6 && length(unique(gen[gtest])) >= 2) {
  bd_g <- as.dist(as.matrix(bd)[gtest, gtest])
  set.seed(42)
  pmg <- adonis2(bd_g ~ gen[gtest], permutations = 9999)
  set.seed(42)
  ang <- anosim(bd_g, factor(gen[gtest]), permutations = 9999)
  cat(sprintf("PERMANOVA genus : F=%.2f R2=%.3f (null %.3f, ratio %.2f) p=%.4f (n=%d taxa, %d genera) | ANOSIM R=%.3f p=%.3f\n",
              pmg$F[1], pmg$R2[1], (sum(gen_n >= 3)-1)/(sum(gtest)-1), pmg$R2[1]/((sum(gen_n >= 3)-1)/(sum(gtest)-1)), pmg$`Pr(>F)`[1], sum(gtest), sum(gen_n >= 3),
              ang$statistic, ang$signif))
} else {
  cat("PERMANOVA genus : too few genera with n>=3 to test\n")
}

TOPF <- sort(names(sort(fam_n, decreasing = TRUE))[seq_len(min(8, length(fam_n)))])
ord <- tibble(taxon = w$taxon,
              family = factor(ifelse(w$family %in% TOPF, w$family,
                                     ifelse(w$family == "Unclassified",
                                            "Unclassified", "Other")),
                              levels = c(TOPF, "Unclassified", "Other")),
              NMDS1 = vegan::scores(nm, display = "sites")[, 1],
              NMDS2 = vegan::scores(nm, display = "sites")[, 2])

# Hulls only for real families with >=3 taxa; with fewer, a polygon is a line or a dot,
# and neither "Unclassified" nor "Other" is a clade whose spread would mean anything.
hull <- ord |> filter(!family %in% c("Other", "Unclassified")) |>
  group_by(family) |> filter(n() >= 3) |> slice(chull(NMDS1, NMDS2)) |> ungroup()

p2b <- ggplot(ord, aes(NMDS1, NMDS2, colour = family)) +
  geom_polygon(data = hull, aes(fill = family), alpha = .07, linewidth = .18) +
  geom_point(size = 1.5, alpha = .95) +
  # Headroom for the stats box. Without it the second line runs into the top of the taxon
  # cloud, which sits high in this ordination.
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.25, size = 1.85,
           colour = "grey35", lineheight = 1.15,
           label = sprintf("stress = %.2f\nPERMANOVA family: R² = %.2f, P = %.3f",
                           nm$stress, r2, pv)) +
  scale_colour_manual(values = c(unname(OK[c("blue", "orange", "green", "purple",
                                             "sky", "red", "brown", "teal")]),
                                 "#5A5A5A", "#BEBEBE")) +
  scale_fill_manual(values = c(unname(OK[c("blue", "orange", "green", "purple",
                                           "sky", "red", "brown", "teal")]),
                               "#5A5A5A", "#BEBEBE"), guide = "none") +
  labs(title = "Ordination of predicted food-utilisation profiles",
       x = "NMDS1", y = "NMDS2") +
  base + theme(panel.grid.major.y = element_line(linewidth = .22, colour = "grey92"),
               panel.border = element_rect(fill = NA, colour = "grey85", linewidth = .3),
               legend.position = "right",
               legend.text = element_text(size = 5.9, face = "italic")) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.7)))

df <- read_csv("fig2c_differential.csv", show_col_types = FALSE) |>
  group_by(taxon) |> slice_max(comp_score, n = 4, with_ties = FALSE) |> ungroup() |>
  # "Baobab fruit /Moneky bread" is a misspelling in the FAO BioFoodComp source (fdc_id
  # 42006685, BioFoodComp code 601346). The resource preserves source strings verbatim and
  # should keep doing so; here we only TRUNCATE the display label rather than reproduce the
  # typo in print.
  mutate(food = sub(" ?/moneky bread", "", food, ignore.case = TRUE),
         # Wrap long source strings onto two lines. facet_wrap sizes each panel to its
         # widest y-axis label, so ONE long name (FAO BioFoodComp catalogues edible
         # insects, e.g. "madagascar hissing cockroach, adult") narrows every other panel
         # until ggplot clips the facet strips and the species names read ". adolescent".
         # Wrapping the label is a display fix only; the data keeps the source string.
         food = stringr::str_wrap(food, width = 22),
         key = fct_reorder(paste0(food, "​", taxon), comp_score))

p2c <- ggplot(df, aes(comp_score, key, fill = taxon)) +
  geom_col(width = .66, show.legend = FALSE) +
  scale_y_discrete(labels = function(x) sub("​.*", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0, .10)), n.breaks = 4) +
  scale_fill_manual(values = unname(OK[c("blue", "orange", "green", "purple", "sky", "red")])) +
  facet_wrap(~taxon, scales = "free", ncol = 3) +
  labs(title = "Differential food scores by taxon",
       x = "differential score", y = NULL) +
  base + theme(strip.text = element_text(size = 6.8, face = "italic",
                                         margin = margin(t = 1, b = 1)),
               panel.spacing.x = unit(4, "mm"), panel.spacing.y = unit(2, "mm"))

fig2 <- (p2a | p2b) / p2c +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

ggsave("Figure2.pdf", fig2, width = 180, height = 176, units = "mm", device = cairo_pdf)
ggsave("Figure2.png", fig2, width = 180, height = 176, units = "mm", dpi = 400)

cat(sprintf("wrote Figure1 (2 panels) and Figure2 (3 panels); NMDS stress %.3f, %d taxa\n",
            nm$stress, nrow(ord)))
