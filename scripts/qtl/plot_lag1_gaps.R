#!/usr/bin/env Rscript
# plot_lag1_gaps.R -- distribution of lag-1 (adjacent-marker) genetic distances in the
# refined map. This is the statistic that decides whether the map is usable: the bulk
# should sit near the expected spacing, and any excess length shows up as a tail.
#
# Current map (error.prob = 2.5e-3): 1,713 markers, 3,388.6 cM total, median gap
# 0.816 cM against 0.869 expected for ~1,700 markers over 1,500 cM. So the BULK is
# right and the ~1,800 cM of excess lives in ~22 gaps over 10 cM -- which the
# error.prob sweep showed are immune to that parameter (max_gap moved only 59.6 ->
# 51.9 cM across a 400x change).
#
# Panels:
#   A  histogram, log10 x -- the whole distribution, so the tail is visible at all
#   B  ECDF with the expected spacing marked
#   C  per chromosome (boxplot, log10) -- is the tail localised?
#   D  cumulative cM contributed, ordered by gap size -- HOW MUCH LENGTH the tail owns
#
# Panel D is the one that matters: it answers "if I capped gaps at X cM, what would the
# map length be", without picking X.
#
# READ-ONLY apart from the plot.
# Usage: Rscript scripts/qtl/plot_lag1_gaps.R > agent/plot_lag1_gaps.log 2>&1

suppressMessages({library(qtl); library(data.table); library(ggplot2); library(patchwork)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")
D <- "data/qtl/derived"
dir.create("output/qtl", showWarnings = FALSE, recursive = TRUE)

cr <- readRDS(file.path(D, "rqtl_cross_map_teonamqc.rds"))
m  <- pull.map(cr)
info <- fread(file.path(D, "abh_all_marker_info.tsv"),
              colClasses = list(character = "chr_v5"))

g <- rbindlist(lapply(names(m), function(ch) {
  v <- m[[ch]][order(m[[ch]])]
  if (length(v) < 2) return(NULL)
  nm <- names(v)
  data.table(chr = ch, from = head(nm, -1), to = tail(nm, -1),
             cM = as.numeric(diff(v)))
}))
g <- merge(g, info[, .(from = marker, pos_from = pos_v5)], by = "from")
g <- merge(g, info[, .(to   = marker, pos_to   = pos_v5)], by = "to")
g[, Mb := (pos_to - pos_from) / 1e6]
g[, chr := factor(chr, levels = as.character(1:10))]
setorder(g, chr, pos_from)

n_mk  <- totmar(cr)
tot   <- sum(sapply(m, function(v) max(v) - min(v)))
exp_sp <- 1500 / n_mk

rule("1. LAG-1 GENETIC DISTANCE")
cat("intervals:", nrow(g), " markers:", n_mk, " total:", round(tot,1), "cM\n")
cat(sprintf("expected spacing at 1500 cM / %d markers: %.3f cM\n", n_mk, exp_sp))
print(round(quantile(g$cM, c(0,.25,.5,.75,.90,.95,.99,1)), 3))
cat("\nzero-distance intervals (co-located markers):", g[cM < 1e-6, .N], "\n")
for (t in c(1, 2, 5, 10, 25, 50))
  cat(sprintf("  gaps > %2d cM: %4d  contributing %7.1f cM (%.1f%% of total)\n",
              t, g[cM > t, .N], g[cM > t, sum(cM)], 100*g[cM > t, sum(cM)]/tot))

rule("2. HOW MUCH LENGTH DOES THE TAIL OWN?")
gs <- g[order(-cM)]
gs[, cum_cM := cumsum(cM)]
gs[, rank := .I]
for (k in c(5, 10, 22, 50, 100))
  cat(sprintf("  the %3d largest intervals hold %7.1f cM (%.1f%% of the map)\n",
              k, gs$cum_cM[k], 100*gs$cum_cM[k]/tot))
cat(sprintf("\nmap length if every gap were capped at 10 cM: %.1f cM\n",
            sum(pmin(g$cM, 10))))
cat(sprintf("map length if every gap were capped at  5 cM: %.1f cM\n",
            sum(pmin(g$cM, 5))))

rule("3. PER CHROMOSOME")
print(g[, .(n = .N, median = round(median(cM),3), q95 = round(quantile(cM,.95),2),
            max = round(max(cM),1), n_gt10 = sum(cM > 10),
            cM_in_gt10 = round(sum(cM[cM > 10]),1)), by = chr][order(as.integer(chr))])

rule("4. PLOT")
pA <- ggplot(g, aes(cM)) +
  geom_histogram(bins = 60, fill = "grey55", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = exp_sp, linetype = 2, colour = "grey20") +
  annotate("text", x = exp_sp, y = Inf, label = sprintf(" expected %.2f cM", exp_sp),
           hjust = 0, vjust = 1.6, size = 3.1) +
  scale_x_log10() +
  labs(title = "A  lag-1 distance (log scale)", x = "cM to next marker", y = "intervals") +
  theme_classic(base_size = 11)

pB <- ggplot(g, aes(cM)) + stat_ecdf(linewidth = 0.5) +
  geom_vline(xintercept = exp_sp, linetype = 2, colour = "grey20") +
  scale_x_log10() +
  labs(title = "B  cumulative distribution", x = "cM to next marker", y = "proportion") +
  theme_classic(base_size = 11)

pC <- ggplot(g, aes(chr, cM)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey90", linewidth = 0.3) +
  geom_hline(yintercept = exp_sp, linetype = 2, colour = "grey20") +
  scale_y_log10() +
  labs(title = "C  per chromosome", x = "chromosome (v5)", y = "cM to next marker") +
  theme_classic(base_size = 11)

pD <- ggplot(gs, aes(rank, 100*cum_cM/tot)) +
  geom_line(linewidth = 0.5) +
  geom_vline(xintercept = 22, linetype = 3, colour = "grey40") +
  annotate("text", x = 22, y = 5, label = " 22 gaps >10 cM", hjust = 0, size = 3.1) +
  labs(title = "D  length owned by the largest intervals",
       x = "intervals, largest first", y = "% of total map length") +
  theme_classic(base_size = 11)

p <- (pA | pB) / (pC | pD) +
  patchwork::plot_annotation(
    title = sprintf("Lag-1 genetic distance -- %d markers, %.0f cM, error.prob = 0.01",
                    n_mk, tot))
ggsave("output/qtl/f2_lag1_gaps.png", p, width = 11, height = 8, dpi = 150, bg = "white")
fwrite(g[, .(chr, from, to, cM, Mb)], file.path(D, "lag1_gaps.tsv"), sep = "\t")
cat("wrote output/qtl/f2_lag1_gaps.png and", file.path(D, "lag1_gaps.tsv"), "\n")
