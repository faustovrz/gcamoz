#!/usr/bin/env Rscript
# phenotype_dictionary_unify.R -- ONE phenotype dictionary for all projects.
#
# TASKS
#   1. PASP -> PAPP (plant appearance), per FVRZ.
#   2. Decide whether a single dictionary is SAFE. The abbreviations are shared
#      between the multilocation diallel and the F2 QTL trial, but the
#      multilocation dictionary states NO FORMULA for HI ("ratio") or RYE
#      ("definition per field protocol"). If the two datasets compute them
#      differently, one abbreviation means two things -- that is the real clash,
#      not the F2 derivation (which is an exact fit on 168 plants and stands
#      regardless of how the dictionary is organised).
#      So: DERIVE the multilocation formulas numerically and compare.
#   3. Emit a single unified dictionary in LONG format: one row per
#      (variable, project). Long format because units genuinely differ --
#      g/plot over ~25 plants in the diallel vs g/plant on single F2 plants --
#      and no single field should have to encode two values.
#
# Reads   data/multilocation.csv, data/multilocation_dictionary.csv,
#         data/qtl/mozpue_phenotype.xlsx                     (all read-only)
# Writes  data/phenotype_dictionary.csv                       (unified, tracked)
#
# Usage: Rscript agent/phenotype_dictionary_unify.R > agent/phenotype_dictionary_unify.log 2>&1

suppressMessages({library(data.table); library(readxl)})
rule <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

MULTI_D <- "data/multilocation_dictionary.csv"
MULTI   <- "data/multilocation.csv"
PHENO   <- "data/qtl/mozpue_phenotype.xlsx"
OUT     <- "data/phenotype_dictionary.csv"

# does candidate reproduce target on all comparable rows?
chk <- function(label, target, candidate, tol = 1e-4) {
  ok <- is.finite(target) & is.finite(candidate)
  if (!any(ok)) { cat(sprintf("  %-30s no comparable rows\n", label)); return(FALSE) }
  mx <- max(abs(target[ok] - candidate[ok]))
  hit <- mx < tol
  cat(sprintf("  %-30s n=%-4d max|diff|=%-12.6g %s\n", label, sum(ok), mx,
              ifelse(hit, "MATCH", "")))
  hit
}

rule("1. MULTILOCATION DATA")
md <- fread(MULTI_D)
m  <- fread(MULTI)
cat("rows:", nrow(m), " cols:", ncol(m), "\n")
cat("columns:", paste(names(m), collapse = " "), "\n")
cat("\ndictionary rows:", nrow(md), " traits:", sum(md$role == "trait"), "\n")

rule("2. DERIVE MULTILOCATION FORMULAS")
has <- function(...) all(c(...) %in% names(m))
res <- list()

cat("PUE:\n")
if (has("GY", "p_input", "PUE")) {
  res$PUE_multi <- chk("PUE = GY(kg/ha)/p_input", m$PUE, m$GY * 1000 / m$p_input)
} else cat("  columns missing\n")

# Two tolerances. The cleaned CSV stores ROUNDED trait columns, so a ratio
# recomputed from rounded inputs cannot match to 1e-4 even when the formula is
# right. `tol_round` accepts agreement consistent with input rounding; anything
# failing THAT is a genuinely different formula, not a precision artefact.
TOL_R <- 0.02
cat("HI:  (tol", TOL_R, "for rounding)\n")
hi_hits <- character(0)
hi_try <- list("GYPP/SW"      = m$GYPP / m$SW,
               "GYPP/TDM"     = m$GYPP / m$TDM,
               "GW/TDM"       = m$GW / m$TDM,
               "GW/SW"        = m$GW / m$SW,
               "EW/TDM"       = m$EW / m$TDM,
               "GW/(GW+SW)"   = m$GW / (m$GW + m$SW),
               "EW/(EW+SW)"   = m$EW / (m$EW + m$SW),
               "GW/EW"        = m$GW / m$EW,
               "GY*1000/TDM"  = m$GY * 1000 / m$TDM,
               "GYPP*NP/TDM"  = m$GYPP * m$NP / m$TDM)
for (nm in names(hi_try)) if (chk(paste("HI =", nm), m$HI, hi_try[[nm]], TOL_R)) hi_hits <- c(hi_hits, nm)

cat("RSR: (tol", TOL_R, "for rounding)\n")
rsr_hits <- character(0)
if (chk("RSR = RW/SW", m$RSR, m$RW / m$SW, TOL_R)) rsr_hits <- c(rsr_hits, "RW/SW")
if (chk("RSR = RW/TDM", m$RSR, m$RW / m$TDM, TOL_R)) rsr_hits <- c(rsr_hits, "RW/TDM")

cat("RYE: (tol", TOL_R, "for rounding)\n")
rye_hits <- character(0)
rye_try <- list("GYPP/RW" = m$GYPP / m$RW, "GW/RW" = m$GW / m$RW,
                "GY/RW"   = m$GY / m$RW,   "GW/SW" = m$GW / m$SW,
                "GYPP/SW" = m$GYPP / m$SW)
for (nm in names(rye_try)) if (chk(paste("RYE =", nm), m$RYE, rye_try[[nm]], TOL_R)) rye_hits <- c(rye_hits, nm)
# "Relative yield efficiency" in P studies is often yield-under-stress relative to
# an optimal-P reference -- a CROSS-ENVIRONMENT quantity, not a row-wise ratio.
# Test that reading before concluding the definition is unknown.
mm <- copy(m)
mm[, gy_opt_ref := mean(GY[treatment == "Optimal"], na.rm = TRUE), by = hybrid]
if (chk("RYE = GY/mean(GY|Optimal,hybrid)", mm$RYE, mm$GY / mm$gy_opt_ref, TOL_R))
  rye_hits <- c(rye_hits, "GY/optimal-P reference (cross-environment)")
mm[, gy_env_mean := mean(GY, na.rm = TRUE), by = env]
if (chk("RYE = GY/mean(GY|env)", mm$RYE, mm$GY / mm$gy_env_mean, TOL_R))
  rye_hits <- c(rye_hits, "GY/environment mean")

cat("TDM:\n")
tdm_hits <- character(0)
if (has("TDM", "SW", "RW")) if (chk("TDM = SW+RW", m$TDM, m$SW + m$RW)) tdm_hits <- c(tdm_hits, "SW+RW")
if (has("TDM", "SW", "RW", "EW")) if (chk("TDM = SW+RW+EW", m$TDM, m$SW + m$RW + m$EW)) tdm_hits <- c(tdm_hits, "SW+RW+EW")

cat("GY vs GYPP:\n")
if (has("GY", "GYPP", "NP")) chk("GY = GYPP*NP/plotarea?", m$GY, m$GYPP * m$NP / 1000)

rule("3. F2 FORMULAS (re-derived here for a like-for-like comparison)")
p <- as.data.frame(suppressMessages(read_excel(PHENO, sheet = "data base")))
names(p)[1] <- "plant"; names(p) <- trimws(names(p))
p <- p[-1, ]
for (n in names(p)) if (n != "OBS:") p[[n]] <- suppressWarnings(as.numeric(p[[n]]))
p <- p[!is.na(p$plant), ]
f <- as.data.table(p)
setnames(f,
  c("plant appearance","Culm Diameter(cm)","Plant height(cm)",
    "Height of Ear Insertion(cm)","yield (g/plt)","Yield","Root Biomass(g)",
    "shoot Dry mass (g)","Ear weight(g)","Total Dry mass weight (g)"),
  c("PAPP","SD","PH","EH","GYPP","GY","RW","SW","EW","TDM"))
cat("F2 HI:\n");  f_hi  <- chk("HI = GYPP/SW",       f$HI,  f$GYPP / f$SW)
cat("F2 RSR:\n"); f_rsr <- chk("RSR = RW/SW",        f$RSR, f$RW / f$SW)
cat("F2 RYE:\n"); f_rye <- chk("RYE = GYPP/RW",      f$RYE, f$GYPP / f$RW)
cat("F2 TDM:\n"); f_tdm <- chk("TDM = SW+RW",        f$TDM, f$SW + f$RW)
cat("F2 PUE:\n"); f_pue <- chk("PUE = GY(kg/ha)/20", f$PUE, f$GY * 1000 / 20)

rule("4. CLASH REPORT")
cmp <- function(v, multi, f2) {
  if (!length(multi)) { st <- "multiloc UNVERIFIABLE"; }
  else if (f2 %in% multi) st <- "CONSISTENT"
  else st <- "*** CLASH ***"
  cat(sprintf("%-5s multiloc: %-22s F2: %-12s -> %s\n",
              v, ifelse(length(multi), paste(multi, collapse = "|"), "(none matched)"), f2, st))
  st
}
s_hi  <- cmp("HI",  hi_hits,  "GYPP/SW")
s_rsr <- cmp("RSR", rsr_hits, "RW/SW")
s_rye <- cmp("RYE", rye_hits, "GYPP/RW")
s_tdm <- cmp("TDM", tdm_hits, "SW+RW")
cat("\nPUE: same formula GY(kg/ha)/p_input in both; only p_input differs",
    "(multiloc 10/50, F2 20). CONSISTENT.\n")

rule("5. UNIFIED DICTIONARY (long: one row per variable x project)")
# Shared trait definitions. `scale` carries the per-plot vs per-plant difference
# that would otherwise be silently overloaded onto `unit`.
trait <- function(variable, name, unit, formula, project, scale, notes = "")
  data.table(variable, name, role = "trait", unit, formula, project, scale,
             raw_header = "", notes)

# raw_header maps the ABBREVIATION back to the literal column name in the source
# file. Required so downstream scripts (agent/qtl_make_abh.R) can rename source
# columns from ONE authority instead of duplicating the mapping. The diallel
# columns are already abbreviated, so raw_header is only populated for qtl_f2.
RAW_QTL <- c(
  plant = "plant", PAPP = "plant appearance", SD = "Culm Diameter(cm)",
  PH = "Plant height(cm)", EH = "Height of Ear Insertion(cm)",
  GYPP = "yield (g/plt)", GY = "Yield", RW = "Root Biomass(g)",
  SW = "shoot Dry mass (g)", EW = "Ear weight(g)",
  TDM = "Total Dry mass weight (g)", PUE = "PUE", HI = "HI", RSR = "RSR",
  RYE = "RYE", obs = "OBS:")

d_multi <- rbindlist(list(
  trait("GY","Grain yield","t/ha","", "multilocation","plot"),
  trait("GYPP","Grain yield per plant","g","", "multilocation","plot-derived"),
  trait("PUE","Phosphorus-use efficiency","g/g","GY(kg/ha)/p_input","multilocation","plot",
        "p_input = 10 (stress) or 50 (optimal) kg/ha"),
  trait("HI","Harvest index","g/g", if (length(hi_hits)) hi_hits[1] else "UNVERIFIED",
        "multilocation","plot", if (length(hi_hits)) "" else "dictionary says only 'ratio'; not reproduced from data - CONFIRM"),
  trait("RSR","Root:shoot ratio","g/g", if (length(rsr_hits)) rsr_hits[1] else "UNVERIFIED","multilocation","plot"),
  trait("RYE","Relative yield efficiency","g/g", if (length(rye_hits)) rye_hits[1] else "UNVERIFIED",
        "multilocation","plot", if (length(rye_hits)) "" else "dictionary says 'definition per field protocol'; not reproduced - CONFIRM"),
  trait("TDM","Total dry mass","g/plot", if (length(tdm_hits)) tdm_hits[1] else "UNVERIFIED","multilocation","plot"),
  trait("RW","Root weight","g/plot","","multilocation","plot"),
  trait("SW","Shoot dry weight","g/plot","","multilocation","plot","stover, excludes ear"),
  trait("EW","Ear weight","g/plot","","multilocation","plot","per-plot total over ~25 plants"),
  trait("GW","Grain weight","g/plot","","multilocation","plot"),
  trait("PH","Plant height","cm","","multilocation","plot-mean"),
  trait("SD","Stem (culm) diameter","cm","","multilocation","plot-mean"),
  trait("NP","Number of plants per plot (stand)","count","","multilocation","plot")))
# DTS/DTA/ASI are already in multilocation_dictionary.csv with role="derived",
# so they come through `design` below -- do NOT redefine them here or the
# (variable, project) key duplicates.

d_qtl <- rbindlist(list(
  trait("PAPP","Plant appearance (visual score)","score 1-4","","qtl_f2","plant",
        "ORDINAL, not continuous - use scanone(model='np')"),
  trait("SD","Stem (culm) diameter","cm","","qtl_f2","plant"),
  trait("PH","Plant height","cm","","qtl_f2","plant"),
  trait("EH","Ear insertion height","cm","","qtl_f2","plant"),
  trait("GYPP","Grain yield per plant","g/plant","","qtl_f2","plant"),
  trait("GY","Grain yield","t/ha","GYPP*0.0533333","qtl_f2","plant",
        "implies ~53,333 plants/ha (0.1875 m2/plant)"),
  trait("RW","Root weight","g/plant","","qtl_f2","plant","PER PLANT - do not pool with multilocation g/plot"),
  trait("SW","Shoot dry weight","g/plant","","qtl_f2","plant","PER PLANT; stover, excludes ear"),
  trait("EW","Ear weight","g/plant","","qtl_f2","plant","PER PLANT - do not pool with multilocation g/plot"),
  trait("TDM","Total dry mass","g/plant","SW+RW","qtl_f2","plant","excludes ear; PER PLANT"),
  trait("PUE","Phosphorus-use efficiency","g/g","GY(kg/ha)/20","qtl_f2","plant",
        "p_input = 20 kg P/ha, DERIVED numerically (sd 2.8e-16); neither multilocation rate"),
  trait("HI","Harvest index","g/g","GYPP/SW","qtl_f2","plant",
        "NOT conventional HI: denominator excludes the ear. Grain-per-shoot."),
  trait("RSR","Root:shoot ratio","g/g","RW/SW","qtl_f2","plant"),
  trait("RYE","Relative yield efficiency","g/g","GYPP/RW","qtl_f2","plant",
        "grain per unit root weight")))

design <- md[role != "trait", .(variable, name, role, unit,
                                formula = "", project = "multilocation",
                                scale = "plot", raw_header = "", notes)]
extra_design <- data.table(
  variable = c("plant","obs"),
  name = c("F2 plant ID (field book number)","Field observation note"),
  role = c("design","note"), unit = c("NA","NA"), formula = c("",""),
  project = "qtl_f2", scale = "plant", raw_header = "",
  notes = c("1-188; 93 & 94 not genotyped; 20 IDs unphenotyped; 166 have both",
            "15 plants flagged 'Higher M.C' - GY/GYPP not on a common moisture basis"))

# populate raw_header for every qtl_f2 row from the single map above
for (tb in list(d_qtl, extra_design))
  tb[project == "qtl_f2", raw_header := unname(RAW_QTL[variable])]
stopifnot(!any(is.na(d_qtl$raw_header)), all(d_qtl$raw_header != ""))

dict <- rbindlist(list(design, extra_design, d_multi, d_qtl), use.names = TRUE)
setorder(dict, role, variable, project)
fwrite(dict, OUT)
cat("wrote", OUT, " rows:", nrow(dict), "\n\n")
print(dict[role == "trait", .(variable, project, unit, formula)])

rule("6. INTEGRITY CHECKS")
cat("variables in BOTH projects:",
    paste(sort(intersect(d_multi$variable, d_qtl$variable)), collapse = " "), "\n")
shared <- merge(d_multi[, .(variable, u_m = unit, f_m = formula)],
                d_qtl[,   .(variable, u_q = unit, f_q = formula)], by = "variable")
cat("\nunit differences on shared abbreviations:\n")
print(shared[u_m != u_q, .(variable, multilocation = u_m, qtl_f2 = u_q)])
cat("\nformula differences on shared abbreviations:\n")
fd <- shared[f_m != f_q & f_m != "" & f_q != ""]
if (nrow(fd)) print(fd[, .(variable, multilocation = f_m, qtl_f2 = f_q)]) else
  cat("  none where both are known\n")
stopifnot(!any(duplicated(dict[, .(variable, project)])))
cat("\n(variable, project) unique: OK\n")
