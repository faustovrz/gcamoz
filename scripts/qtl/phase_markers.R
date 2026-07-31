#!/usr/bin/env Rscript
# qtl_phase_markers.R -- make marker polarity consistent in LINKAGE PHASE.
#
# WHY: rf between two markers is invariant only if BOTH are flipped. Flipping ONE
# turns rf into 1-rf, so an out-of-phase marker looks unlinked to its neighbours.
# With ~47% of junctions anti-phase, est.map returned 161,170 cM against a maize
# expectation of ~1,400 (agent/qtl_build_map.log). Reading CML530's allele directly
# from the assembly cut that to ~6.6% (93.3% in phase). This closes the rest.
#
# METHOD: start from the CML530 polarity already in place -- it is tied to a real
# parental haplotype, so it is globally consistent, not merely internally
# consistent -- then flip only markers that disagree with their LOCAL NEIGHBOURHOOD.
#
# Deliberately NOT a naive chain (flip i+1 if cor(i, i+1) < 0): one bad marker in
# such a chain propagates a flip to every marker downstream of it. Instead each
# marker is compared against a WINDOW of up to K neighbours on each side, and
# flipped only if its mean correlation with them is negative. Single bad markers
# get corrected; they cannot corrupt the rest of the chromosome. Iterated to
# convergence.
#
# The CML530 allele still does the job it is actually good for: NAMING which
# haplotype is the parent's, i.e. the SIGN of allele effects. Phase comes from
# linkage. I had this backwards originally.
#
# Reads   data/qtl/derived/rqtl_gen_abh.csv, rqtl_phe.csv, abh_marker_info.tsv
# Writes  data/qtl/derived/rqtl_gen_abh_phased.csv
#         data/qtl/derived/phase_flips.tsv
# Usage: Rscript scripts/qtl/qtl_phase_markers.R > agent/qtl_phase_markers.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

D       <- "data/qtl/derived"
GENFILE <- file.path(D, "rqtl_gen_abh.csv")
PHEFILE <- file.path(D, "rqtl_phe.csv")
MKINFO  <- file.path(D, "abh_marker_info.tsv")
OUTGEN  <- file.path(D, "rqtl_gen_abh_phased.csv")
OUTFLIP <- file.path(D, "phase_flips.tsv")
K       <- 5L      # neighbours each side
MAXIT   <- 20L

rule("1. LOAD")
cr <- read.cross(format = "csvsr", dir = "", genfile = GENFILE, phefile = PHEFILE,
                 genotypes = c("A","H","B"), na.strings = c("-","NA"),
                 crosstype = "f2", estimate.map = FALSE)
mkinfo <- fread(MKINFO, colClasses = list(character = "chr_v5"))
G <- pull.geno(cr)                                  # ind x marker, 1/2/3
mk_chr <- rep(names(cr$geno), nmar(cr))
mk_nm  <- colnames(G)
dose <- G - 1L                                      # 0,1,2
cat("individuals:", nrow(dose), " markers:", ncol(dose), "\n")
cat("marker order is PHYSICAL (chr_v5, pos_v5), fixed by qtl_make_abh.R\n")

phase_consistency <- function(d, chrv) {
  r <- unlist(lapply(unique(chrv), function(ch) {
    i <- which(chrv == ch); if (length(i) < 2) return(NULL)
    sapply(seq_len(length(i) - 1), function(k) {
      a <- d[, i[k]]; b <- d[, i[k + 1]]
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 30) return(NA_real_)
      cor(a[ok], b[ok])
    })
  }))
  r[is.finite(r)]
}
r0 <- phase_consistency(dose, mk_chr)
cat(sprintf("\nBEFORE: %.1f%% of adjacent junctions in phase (median cor %.3f)\n",
            100 * mean(r0 > 0), median(r0)))
cat(sprintf("        anti-phase junctions: %d of %d\n", sum(r0 < 0), length(r0)))

rule("1b. DROP BADLY-LINKED MARKERS *BEFORE* PHASING")
# ORDER MATTERS AND I HAD IT WRONG. Phasing first, filtering second, meant the
# windowed correction was applied to markers whose correlations are pure noise: of
# the 53 markers the first version flipped, 39.6% were weakly linked, against 2.4%
# of unflipped markers -- a 16x enrichment (agent/qtl_diagnose_chr4_chr9.log).
# Flipping a noisy marker relabels noise; it does not fix it, and it corrupted the
# genotypes the downstream filter then had to judge.
#
# THRESHOLD IS DATA-DRIVEN, NOT FITTED TO THE MAIZE EXPECTATION. The per-marker
# max |cor| distribution has a cliff: 2nd percentile 0.141, 5th percentile 0.870.
# A cut at 0.80 removes the broken tail (~4%); 0.90 would remove 8.5% and start
# taking healthy markers.
#
# Using the MAX over a +/-K window (not the adjacent junction) matters: a marker
# sitting beside a 35 Mb physical gap is legitimately weakly correlated with THAT
# neighbour, and weak-junction rate does rise with gap size (6.6% under 0.5 Mb ->
# 33% over 30 Mb). Taking the max over several neighbours avoids penalising a
# marker for one large gap.
LINK_MIN <- 0.80
LINK_K   <- 5L
maxcor <- function(d, chrv) {
  out <- rep(NA_real_, ncol(d))
  for (ch in unique(chrv)) {
    idx <- which(chrv == ch); n <- length(idx); dd <- d[, idx, drop = FALSE]
    for (k in seq_len(n)) {
      nb <- setdiff(max(1, k - LINK_K):min(n, k + LINK_K), k)
      cs <- sapply(nb, function(j) {
        a <- dd[, k]; b <- dd[, j]; ok <- is.finite(a) & is.finite(b)
        if (sum(ok) < 30) return(NA_real_)
        abs(cor(a[ok], b[ok]))
      })
      out[idx[k]] <- suppressWarnings(max(cs, na.rm = TRUE))
    }
  }
  out
}
best <- maxcor(dose, mk_chr)
cat("per-marker max |cor| with any neighbour:\n")
print(round(quantile(best, c(0, .01, .02, .05, .10, .25, .50, 1), na.rm = TRUE), 3))
drop_link <- is.finite(best) & best < LINK_MIN
cat(sprintf("\ndropping %d markers (%.1f%%) below |cor| %.2f\n",
            sum(drop_link), 100 * mean(drop_link), LINK_MIN))
cat("per chromosome:\n"); print(table(mk_chr[drop_link]))
keep_mk <- mk_nm[!drop_link]
dose   <- dose[, !drop_link, drop = FALSE]
mk_chr <- mk_chr[!drop_link]
mk_nm  <- mk_nm[!drop_link]
cat("markers retained:", length(mk_nm), "\n")
r0b <- phase_consistency(dose, mk_chr)
cat(sprintf("phase consistency after the drop, before phasing: %.1f%% (was %.1f%%)\n",
            100 * mean(r0b > 0), 100 * mean(r0 > 0)))

rule("2. ITERATIVE WINDOWED PHASE CORRECTION")
flipped <- rep(FALSE, ncol(dose))
for (it in seq_len(MAXIT)) {
  nflip <- 0L
  for (ch in unique(mk_chr)) {
    idx <- which(mk_chr == ch); n <- length(idx)
    if (n < 3) next
    d <- dose[, idx, drop = FALSE]
    # mean correlation of each marker with up to K neighbours each side
    score <- sapply(seq_len(n), function(k) {
      nb <- setdiff(max(1, k - K):min(n, k + K), k)
      cs <- sapply(nb, function(j) {
        a <- d[, k]; b <- d[, j]; ok <- is.finite(a) & is.finite(b)
        if (sum(ok) < 30) return(NA_real_)
        cor(a[ok], b[ok])
      })
      mean(cs, na.rm = TRUE)
    })
    bad <- which(is.finite(score) & score < 0)
    if (length(bad)) {
      gi <- idx[bad]
      dose[, gi] <- 2L - dose[, gi]          # flip: 0<->2, 1 unchanged
      flipped[gi] <- !flipped[gi]
      nflip <- nflip + length(bad)
    }
  }
  r <- phase_consistency(dose, mk_chr)
  cat(sprintf("iter %2d: flipped %4d  ->  %.1f%% in phase, median cor %.3f\n",
              it, nflip, 100 * mean(r > 0), median(r)))
  if (nflip == 0L) break
}
r1 <- phase_consistency(dose, mk_chr)
cat(sprintf("\nAFTER: %.1f%% in phase (median cor %.3f)\n",
            100 * mean(r1 > 0), median(r1)))
cat(sprintf("       anti-phase junctions: %d of %d  (was %d)\n",
            sum(r1 < 0), length(r1), sum(r0 < 0)))
cat("total markers flipped:", sum(flipped), sprintf("(%.1f%%)\n", 100*mean(flipped)))

rule("3. WHICH MARKERS WERE FLIPPED?")
fl <- data.table(marker = mk_nm, chr = mk_chr, flipped = flipped)  # mk_nm already filtered
fl <- merge(fl, mkinfo[, .(marker, pos_v5, polarity_source, cml_class)],
            by = "marker", all.x = TRUE)
print(fl[, .(n = .N, n_flipped = sum(flipped),
             pct = round(100*mean(flipped), 1)), by = polarity_source])
cat("\nby CML530 class:\n")
print(fl[, .(n = .N, n_flipped = sum(flipped),
             pct = round(100*mean(flipped), 1)), by = cml_class][order(-n)])
cat("\nper chromosome:\n")
print(fl[, .(n = .N, n_flipped = sum(flipped),
             pct = round(100*mean(flipped), 1)), by = chr][order(as.integer(chr))])
cat("\nINTERPRETATION: markers on ARBITRARY polarity should be flipped ~50% of the\n")
cat("time (they were coin flips). A LOW flip rate among CML530-polarized markers\n")
cat("confirms the assembly read was right; a high rate would indict it.\n")

rule("4. WRITE PHASED GENOTYPES")
# back to A/H/B; dose 0=A, 1=H, 2=B
ABH <- matrix(c("A","H","B")[dose + 1L], nrow = nrow(dose),
              dimnames = dimnames(dose))
ABH[is.na(dose)] <- NA
gen <- fread(GENFILE, sep = ",", header = TRUE)
ids <- names(gen)[-(1:2)]
# Do NOT assume row order matches. read.cross sorts chromosomes NUMERICALLY
# (1,2,...,10) whereas the file is written in lexicographic order (1,10,2,...),
# because chr_v5 is a character column. So pull.geno's column order differs from
# the file's row order. Match by NAME, and assert full coverage both ways.
# MARKERS: match by name -- read.cross sorts chromosomes numerically (1,2,...,10)
#   while the file is lexicographic (1,10,2,...), so orders differ.
# INDIVIDUALS: positional. read.cross does NOT reorder individuals, and
#   pull.geno() does not attach their IDs as rownames, so name matching is not
#   available. Assert the count instead.
cat("pull.geno rownames present:", !is.null(rownames(dose)), "\n")
gen <- gen[id %in% keep_mk]                     # drop filtered markers from output
stopifnot(setequal(gen$id, mk_nm), length(ids) == nrow(dose))
ABHt <- t(ABH)                                   # markers x individuals
ABHt <- ABHt[match(gen$id, rownames(ABHt)), , drop = FALSE]
stopifnot(identical(rownames(ABHt), gen$id), ncol(ABHt) == length(ids))
out <- data.table(id = gen$id, chr = gen$chr)
out <- cbind(out, as.data.table(ABHt))
setnames(out, c("id", "chr", ids))
fwrite(out, OUTGEN, na = "-", quote = FALSE)
fwrite(fl, OUTFLIP, sep = "\t")
cat("wrote", OUTGEN, " dim:", paste(dim(out), collapse = " x "), "\n")
cat("wrote", OUTFLIP, "\n")

rule("5. VERIFY ROUND-TRIP")
cr2 <- read.cross(format = "csvsr", dir = "", genfile = OUTGEN, phefile = PHEFILE,
                  genotypes = c("A","H","B"), na.strings = c("-","NA"),
                  crosstype = "f2", estimate.map = FALSE)
cat("nind:", nind(cr2), " markers:", totmar(cr2), " nchr:", nchr(cr2), "\n")
g2 <- pull.geno(cr2)
r2 <- phase_consistency(g2 - 1L, rep(names(cr2$geno), nmar(cr2)))
cat(sprintf("phase consistency after round-trip: %.1f%%\n", 100*mean(r2 > 0)))
# marker count is INTENTIONALLY lower than the input: section 1b drops badly-linked
# markers. Assert against the retained set, not the original.
stopifnot(nind(cr2) == nind(cr), totmar(cr2) == length(keep_mk))
tb <- table(factor(g2, levels = 1:3))
cat(sprintf("genotype composition AA %.1f%% AB %.1f%% BB %.1f%%\n",
            100*tb[1]/sum(tb), 100*tb[2]/sum(tb), 100*tb[3]/sum(tb)))
cat("\nNOTE: A/B are now LINKAGE-PHASED haplotypes. Which one is CML530's is\n")
cat("recorded per marker in phase_flips.tsv (polarity_source x flipped) and is\n")
cat("needed only for the SIGN of allele effects, not for the map or QTL positions.\n")
