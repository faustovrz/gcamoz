#!/usr/bin/env Rscript
# qtl_plot_map.R -- plot the refined F2 map, same style as zealhmm's
# analysis/teonam-genetic-map.qmd (chromosome ideograms with marker ticks).
# Usage: Rscript scripts/qtl/qtl_plot_map.R > agent/qtl_plot_map.log 2>&1

suppressMessages({library(qtl); library(data.table); library(ggplot2)})
D <- "data/qtl/derived"
dir.create("output/qtl", showWarnings = FALSE, recursive = TRUE)

cr  <- readRDS(file.path(D, "rqtl_cross_map.rds"))
map <- pull.map(cr)
chrs <- names(map)

# zero each chromosome. est.map PRESERVES the input chromosome offset, and the input
# was read.cross's dummy 5 cM spacing (chr minima came out as 25, 5, 60, 5, 50, 130,
# 75, 10, 135, 55). Plotting absolute cM put chr6 at 130-170 with a bar from 0.
mk <- rbindlist(lapply(chrs, function(ch)
  data.table(chr = ch, cm = map[[ch]] - min(map[[ch]]))))
mk[, chr := factor(chr, levels = as.character(1:10))]

tab <- mk[, .(markers = .N, length_cM = round(max(cm) - min(cm), 1)), by = chr][order(chr)]
print(rbind(tab[, .(chr = as.character(chr), markers, length_cM)],
            data.table(chr = "total", markers = sum(tab$markers),
                       length_cM = round(sum(tab$length_cM), 1))))

p <- ggplot(mk, aes(x = chr, y = cm)) +
  geom_segment(aes(xend = chr, y = 0, yend = cm), linewidth = 6, colour = "grey85") +
  geom_segment(aes(xend = chr, y = cm - 0.15, yend = cm + 0.15),
               linewidth = 6, colour = "grey30", alpha = 0.20) +
  scale_y_reverse(expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "chromosome (v5)", y = "position (cM)",
       title = sprintf("F2 refined map -- %d markers, %.1f cM total",
                       sum(tab$markers), sum(tab$length_cM))) +
  theme_classic(base_size = 14)
ggsave("output/qtl/f2_map_ideogram.png", p, width = 9, height = 6, dpi = 150, bg = "white")

# Marey map: physical vs genetic. Non-monotonic = order problem; flat = no recombination.
info <- fread(file.path(D, "abh_marker_info.tsv"), colClasses = list(character = "chr_v5"))
mk2 <- rbindlist(lapply(chrs, function(ch)
  data.table(chr = ch, marker = names(map[[ch]]),
             cm = as.numeric(map[[ch]]) - min(map[[ch]]))))
mk2 <- merge(mk2, info[, .(marker, pos_v5)], by = "marker")
mk2[, `:=`(mb = pos_v5 / 1e6, chr = factor(chr, levels = as.character(1:10)))]
q <- ggplot(mk2[order(chr, mb)], aes(mb, cm)) +
  geom_line(colour = "grey60", linewidth = 0.4) +
  geom_point(size = 0.5, colour = "grey20") +
  facet_wrap(~chr, scales = "free", ncol = 5) +
  labs(x = "physical position (Mb, B73 v5)", y = "genetic position (cM)",
       title = "Marey map -- monotonic and sigmoid is healthy") +
  theme_classic(base_size = 11)
ggsave("output/qtl/f2_map_marey.png", q, width = 12, height = 5.5, dpi = 150, bg = "white")
cat("\nwrote output/qtl/f2_map_ideogram.png and f2_map_marey.png\n")
