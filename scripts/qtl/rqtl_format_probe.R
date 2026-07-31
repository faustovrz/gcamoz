#!/usr/bin/env Rscript
# qtl_rqtl_format_probe.R -- settle the R/qtl csvsr file layout EMPIRICALLY.
#
# WHY: agent/qtl_make_abh.R wrote a csvsr genotype file as
#        marker, chr, pos, <one column per individual>
# and read.cross() reported 167 individuals instead of 166, warning
#   "1 individuals with genotypes but no phenotypes: pos"
# and treating the Mb positions as unknown genotype codes. So read.cross consumed
# only TWO leading columns from the genfile, not three. Rather than guess from the
# docs (which describe a `mapfile` for non-csv formats), build a TINY cross with a
# KNOWN answer, write it in each candidate layout, and see which one round-trips.
#
# Ground truth used below: 4 markers on 2 chromosomes, 3 individuals, 2 phenotypes.
# A layout is CORRECT only if read.cross returns nind=3, totmar=4, nchr=2 AND the
# genotypes come back matching what was written.
#
# READ-ONLY with respect to project data; writes only into a temp dir.
# Usage: Rscript agent/qtl_rqtl_format_probe.R > agent/qtl_rqtl_format_probe.log 2>&1

suppressMessages({library(qtl); library(data.table)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

TD <- tempfile("rqtlfmt"); dir.create(TD)
IND    <- c("ind1", "ind2", "ind3")
MARK   <- c("m1", "m2", "m3", "m4")
CHR    <- c("1", "1", "2", "2")
POS    <- c(0.0, 10.5, 0.0, 7.25)
# genotypes: markers x individuals
GT <- matrix(c("A","H","B",
               "H","B","A",
               "B","A","H",
               "A","A","B"), nrow = 4, byrow = TRUE,
             dimnames = list(MARK, IND))
PHE <- data.frame(trait1 = c(1.5, 2.5, 3.5), trait2 = c(10, 20, 30), row.names = IND)

rule("0. GROUND TRUTH")
cat("individuals:", length(IND), " markers:", length(MARK), " chromosomes:",
    length(unique(CHR)), " phenotypes:", ncol(PHE), "\n\n")
print(GT); cat("\n"); print(PHE)

# csvsr phenotype file: rotated -> one ROW per phenotype, one column per individual
phefile <- file.path(TD, "phe.csv")
phe_dt <- data.table(id = colnames(PHE))
phe_dt <- cbind(phe_dt, as.data.table(t(as.matrix(PHE))))
setnames(phe_dt, c("id", IND))
fwrite(phe_dt, phefile, quote = FALSE)
cat("\nphefile written:\n"); cat(readLines(phefile), sep = "\n"); cat("\n")

try_layout <- function(label, gen_dt) {
  genfile <- file.path(TD, paste0("gen_", gsub("[^a-z0-9]", "", tolower(label)), ".csv"))
  fwrite(gen_dt, genfile, quote = FALSE, na = "-")
  cat("\n", strrep("-", 68), "\n", label, "\n", strrep("-", 68), "\n", sep = "")
  cat("genfile:\n"); cat(readLines(genfile), sep = "\n"); cat("\n\n")
  res <- tryCatch({
    cr <- suppressWarnings(read.cross(format = "csvsr", dir = "", genfile = genfile,
              phefile = phefile, genotypes = c("A","H","B"),
              na.strings = c("-","NA"), crosstype = "f2", estimate.map = FALSE))
    ok_dim <- nind(cr) == length(IND) && totmar(cr) == length(MARK) && nchr(cr) == length(unique(CHR))
    # do the genotypes round-trip? pull.geno returns ind x marker, coded 1/2/3
    pg <- pull.geno(cr)
    back <- matrix(c("A","H","B")[pg], nrow = nrow(pg), dimnames = dimnames(pg))
    ok_gt <- identical(unname(back), unname(t(GT)[, colnames(back), drop = FALSE]))
    mp <- unlist(pull.map(cr))
    cat(sprintf("nind=%d  totmar=%d  nchr=%d  dims_ok=%s  genotypes_ok=%s\n",
                nind(cr), totmar(cr), nchr(cr), ok_dim, ok_gt))
    cat("map positions read back:", paste(round(mp, 3), collapse = " "), "\n")
    cat("expected positions     :", paste(POS, collapse = " "), "\n")
    cat("positions_ok:", isTRUE(all.equal(unname(mp), POS)), "\n")
    ok_dim && ok_gt
  }, error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); FALSE })
  invisible(res)
}

rule("1. CANDIDATE LAYOUTS")

# (a) three leading columns: marker, chr, pos  -- what qtl_make_abh.R did
a <- data.table(id = MARK, chr = CHR, pos = POS); a <- cbind(a, as.data.table(GT))
r_a <- try_layout("A) genfile = marker, chr, pos, <individuals>", a)

# (b) two leading columns: marker, chr  -- what the failure implied
b <- data.table(id = MARK, chr = CHR); b <- cbind(b, as.data.table(GT))
r_b <- try_layout("B) genfile = marker, chr, <individuals>", b)

rule("2. VERDICT")
cat("layout A (marker, chr, pos, inds):", ifelse(r_a, "CORRECT", "WRONG"), "\n")
cat("layout B (marker, chr, inds)     :", ifelse(r_b, "CORRECT", "WRONG"), "\n\n")
if (r_a && !r_b) {
  cat("=> keep three leading columns; the 167-individual bug in qtl_make_abh.R\n")
  cat("   must have another cause (check for a stray column or name mismatch).\n")
} else if (r_b && !r_a) {
  cat("=> csvsr genfile takes TWO leading columns (marker, chr). Positions are\n")
  cat("   NOT carried in the genfile for this format. Fix qtl_make_abh.R by\n")
  cat("   dropping the pos column and supplying the map separately (or rely on\n")
  cat("   est.map, which is the plan anyway since cM comes from recombination).\n")
} else if (r_a && r_b) {
  cat("=> both parsed; prefer A since it preserves positions.\n")
} else {
  cat("=> NEITHER layout round-tripped. Inspect the output above before trusting\n")
  cat("   any csvsr file. Consider format='csvs' (unrotated) instead.\n")
}

rule("3. NOTE ON POSITIONS")
cat("cM positions are a placeholder regardless: the agreed approach is to fix\n")
cat("marker ORDER from the B73 physical map and derive distances with est.map.\n")
cat("So a layout that cannot carry positions is not a real loss -- but it must be\n")
cat("known, not assumed, because a silently mis-parsed column becomes a fake\n")
cat("genotype column (exactly what produced the 167th 'individual' named 'pos').\n")
unlink(TD, recursive = TRUE)
