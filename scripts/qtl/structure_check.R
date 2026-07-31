#!/usr/bin/env Rscript
# qtl_structure_check.R -- two questions the file inventory alone cannot answer:
#
#   Q-A) Are either of the two parental inbreds present among the 186 samples?
#        An inbred parent must show near-zero heterozygosity genome-wide.
#        An F2 individual shows ~0.5 het at informative markers.
#
#   Q-B) Are the 186 samples ONE biparental F2, or F2s pooled from several of
#        the 49 crosses? Stelio's email says "F2 populations ... from all
#        combinations" (plural), which would break the R/qtl `f2` cross type.
#        Three independent diagnostics below.
#
# READ-ONLY. stdout -> agent/qtl_structure_check.log
# Run:  Rscript agent/qtl_structure_check.R > agent/qtl_structure_check.log 2>&1

suppressMessages({library(data.table)})
SNP1ROW <- "data/qtl/Report_DMz26-3123_SNP_mapping_2.csv"
DART_SKIP <- 6
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep="")

d <- fread(SNP1ROW, skip = DART_SKIP, header = TRUE, na.strings = c("", "NA"))
mc <- 1:which(names(d) == "RepAvg")
G  <- as.matrix(d[, -mc, with = FALSE])
samples <- colnames(G)
plate <- unlist(fread(SNP1ROW, nrows = DART_SKIP, header = FALSE,
                      colClasses = "character")[3, (length(mc) + 1):(length(mc) + ncol(G)),
                                                with = FALSE])

chr <- sub(" .*", "", d$Chrom_Maize_B73_V4.0.assembly)
pos <- suppressWarnings(as.numeric(d$ChromPosSnp_Maize_B73_V4.0.assembly))
n0 <- rowSums(G == "0"); n1 <- rowSums(G == "1"); n2 <- rowSums(G == "2")
tot <- n0 + n1 + n2
# p = frequency of the reference allele, counting the het class as one of each
pref <- (2 * n0 + n2) / (2 * tot)
maf  <- pmin(pref, 1 - pref)
cs   <- (n0 - tot/4)^2/(tot/4) + (n1 - tot/4)^2/(tot/4) + (n2 - tot/2)^2/(tot/2)
p121 <- pchisq(cs, 2, lower.tail = FALSE)

keep <- chr %in% as.character(1:10) & d$CallRate >= 0.90 & d$RepAvg >= 0.95 &
        maf >= 0.15 & p121 > 0.01
keep[is.na(keep)] <- FALSE
dt <- data.table(i = which(keep), clone = d$CloneID[keep], cr = d$CallRate[keep],
                 chr = chr[keep], pos = pos[keep])
fw <- dt[order(-cr), .SD[1], by = clone][order(as.integer(chr), pos)]
Gf <- G[fw$i, , drop = FALSE]
cat("framework markers:", nrow(Gf), " samples:", ncol(Gf), "\n")

## ============================================================ Q-A: any parent?
rule("Q-A. IS EITHER PARENTAL INBRED AMONG THE 186 SAMPLES?")
# genome-wide het over ALL markers (the relevant scale for spotting an inbred)
het_all <- colMeans(G == "2") / colMeans(G != "-")
# het over informative (framework) markers only: F2 expectation ~0.5,
# inbred expectation ~0.0
het_inf <- colMeans(Gf == "2") / colMeans(Gf != "-")
cat("heterozygosity over ALL 5,208 markers:\n"); print(round(summary(het_all), 4))
cat("\nheterozygosity over the", nrow(Gf), "INFORMATIVE framework markers:\n")
print(round(summary(het_inf), 4))
cat("\nF2 expectation at informative markers = 0.50",
    "(less, here, because low read depth undercalls hets)\n")
cat("inbred-parent expectation              = 0.00 - 0.05\n\n")
o <- order(het_inf)[1:10]
cat("TEN LOWEST-heterozygosity samples (the only possible parent candidates):\n")
print(data.frame(sample = samples[o], het_informative = round(het_inf[o], 4),
                 het_all = round(het_all[o], 4), callrate = round(colMeans(G != "-")[o], 3),
                 plate = plate[o]), row.names = FALSE)
cat("\nsamples with informative-marker het < 0.10 :", sum(het_inf < 0.10), "\n")
cat("samples with informative-marker het < 0.20 :", sum(het_inf < 0.20), "\n")
cat("=> if this count is 0, no inbred parent was genotyped in this order.\n")

## =============================================== Q-B1: allele-frequency test
rule("Q-B1. ONE CROSS OR MANY? -- allele frequency at informative markers")
# In a single F2 (F1 selfed), every informative marker must have p(ref) ~ 0.5.
# Pooling F2s from crosses with different parents gives a SPREAD of p, because a
# marker informative in cross A may be fixed for either allele in cross B.
cat("p(ref allele) across the", nrow(Gf), "framework markers:\n")
print(round(summary(pref[fw$i]), 4))
cat("\nSD of p(ref) =", round(sd(pref[fw$i]), 4), "\n")
cat("expected SD for a single F2, n=186, purely binomial sampling =",
    round(sqrt(0.25 / (2 * 186)), 4), "\n")
cat("\nfraction of framework markers with p(ref) in 0.45-0.55 :",
    round(mean(abs(pref[fw$i] - 0.5) < 0.05), 3), "\n")
cat("fraction in 0.40-0.60                                  :",
    round(mean(abs(pref[fw$i] - 0.5) < 0.10), 3), "\n")
cat("\nhistogram of p(ref):\n")
print(table(cut(pref[fw$i], breaks = seq(0, 1, 0.05))))

## ==================================================== Q-B2: PCA for structure
rule("Q-B2. ONE CROSS OR MANY? -- PCA on the framework genotypes")
# dosage: 0 = hom ref, 1 = het, 2 = hom alt; missing -> marker mean
X <- matrix(NA_real_, nrow = ncol(Gf), ncol = nrow(Gf),
            dimnames = list(samples, fw$clone))
X[] <- t(ifelse(Gf == "0", 0, ifelse(Gf == "2", 1, ifelse(Gf == "1", 2, NA))))
for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
pc <- prcomp(X, center = TRUE, scale. = FALSE)
ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat("variance explained by PC1-10 (%):\n"); print(round(ve[1:10], 2))
cat("\nA single F2 gives a smooth decay with a small PC1 (no discrete clusters).\n")
cat("k pooled F2 families give k-1 large, well-separated leading PCs.\n\n")
cat("PC1 range:", round(range(pc$x[,1]), 1), " PC2 range:", round(range(pc$x[,2]), 1), "\n")
# a crude cluster test: is 2-means separation on PC1-2 better than chance?
km <- kmeans(pc$x[, 1:2], centers = 2, nstart = 25)
cat("\n2-means on PC1-2: between/total SS =", round(km$betweenss / km$totss, 3),
    "(cluster sizes", paste(km$size, collapse = "/"), ")\n")
cat("gap statistic proxy -- ratio of PC1 sd to PC2 sd:",
    round(pc$sdev[1] / pc$sdev[2], 3), "(~1 means no dominant axis)\n")
cat("\nPC1 vs PLATE (checks for a batch effect, not population structure):\n")
print(round(tapply(pc$x[, 1], plate, mean), 2))
print(t.test(pc$x[, 1] ~ plate)$p.value)

## ============================================== Q-B3: duplicate / clone check
rule("Q-B3. DUPLICATE SAMPLES")
# identity-by-state on shared non-missing calls
n <- ncol(Gf); best <- data.table()
ibs <- matrix(NA_real_, n, n, dimnames = list(samples, samples))
for (a in 1:(n-1)) for (b in (a+1):n) {
  ok <- Gf[, a] != "-" & Gf[, b] != "-"
  ibs[a, b] <- mean(Gf[ok, a] == Gf[ok, b])
}
v <- ibs[upper.tri(ibs)]
cat("pairwise IBS across", length(v), "pairs:\n"); print(round(summary(v), 4))
w <- which(ibs > 0.90, arr.ind = TRUE)
cat("\npairs with IBS > 0.90 (would indicate duplicated samples):", nrow(w), "\n")
if (nrow(w)) print(data.frame(a = samples[w[,1]], b = samples[w[,2]],
                              ibs = round(ibs[w], 4)), row.names = FALSE)
cat("\nexpected IBS between two unrelated F2 sibs at informative markers ~0.375-0.40\n")

cat("\nstructure check complete.\n")
