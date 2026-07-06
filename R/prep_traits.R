# =============================================================================
# prep_traits.R  —  shared data preparation for the PUE hybrid-analysis paper
# -----------------------------------------------------------------------------
# Single source of truth sourced by every analysis notebook. It:
#   * reads data/multilocation.csv (cleaned; see data/multilocation_dictionary.csv),
#   * builds the design factors and the serpentine field row/column,
#   * derives the flowering-time traits (DTA, DTS, ASI) from the recorded dates,
#   * defines the 11 analysed phenotypes (Joe Gage's list) with groups/units,
#   * exposes two analysis frames: `plots` (all, incl. check) and `hybrids`.
# Objects created: raw_multiloc, field_multiloc, plots, hybrids, trait_meta,
#                  analysis_traits, iqr_outliers().
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# --- read the cleaned data of record ----------------------------------------
raw_multiloc <- readr::read_csv("data/multilocation.csv", show_col_types = FALSE) %>%
  mutate(
    file_order = row_number(),                                  # as-entered plot order
    env        = factor(env),                                   # 5 environments (loc x P-treatment)
    loc        = factor(loc),                                   # 4 physical sites
    treatment  = factor(treatment, levels = c("Optimal", "Stress")),
    rep        = factor(rep, levels = c("I", "II", "III")),
    gen        = factor(gen, levels = as.character(0:49)),
    female     = factor(female,
                        levels = c("Check", "CML364", "CML366", "CML434",
                                   "CML435", "CML439", "CML530", "CML532")),
    male       = factor(male,                                   # IIAM donor IDs (legacy MOZL)
                        levels = c("Check", "EN17", "EN20", "EN21",
                                   "EN25", "EN31", "EN32", "EN64")),
    across(c(silking_date, anthesis_date, planting_date), as.Date)
  )

# --- serpentine field coordinates (5 cols x 10 rows per env x rep) -----------
# ENTRY is the planting order and increases down the file within each env x rep,
# so file order recovers the layout even where the `entry` column has gaps.
field_multiloc <- raw_multiloc %>%
  group_by(env, rep) %>%
  arrange(file_order, .by_group = TRUE) %>%
  mutate(
    plot_order = row_number(),
    field_row  = ceiling(plot_order / 5),
    pos_in_row = ((plot_order - 1) %% 5) + 1,
    field_col  = if_else(field_row %% 2 == 1, pos_in_row, 6L - pos_in_row)
  ) %>%
  ungroup() %>%
  # derived flowering traits: days from that environment's planting date
  mutate(
    DTA = as.integer(anthesis_date - planting_date),   # days to anthesis (male flowering)
    DTS = as.integer(silking_date  - planting_date),   # days to silking (female flowering)
    ASI = DTS - DTA                                     # anthesis-silking interval
  )

# --- the 11 analysed phenotypes (Joe Gage's list), grouped -------------------
# group A = hybrid selection, B = physiological explanation, C = auxiliary.
trait_meta <- tribble(
  ~trait, ~name,                        ~unit,              ~group,
  "GY",   "Grain yield",                "t ha⁻¹", "A – Selection",
  "PUE",  "Phosphorus-use efficiency",  "kg kg⁻¹ P", "A – Selection",
  "HI",   "Harvest index",              "g g⁻¹",  "B – Physiological",
  "RSR",  "Root:shoot ratio",           "g g⁻¹",  "B – Physiological",
  "RW",   "Root weight",                "g plot⁻¹", "B – Physiological",
  "SW",   "Shoot dry weight",           "g plot⁻¹", "B – Physiological",
  "TDM",  "Total dry mass",             "g plot⁻¹", "B – Physiological",
  "PH",   "Plant height",               "cm",               "C – Auxiliary",
  "DTA",  "Days to anthesis",           "days",             "C – Auxiliary",
  "DTS",  "Days to silking",            "days",             "C – Auxiliary",
  "ASI",  "Anthesis–silking interval", "days",         "C – Auxiliary"
) %>%
  mutate(group = factor(group, levels = c("A – Selection",
                                          "B – Physiological",
                                          "C – Auxiliary")),
         trait = factor(trait, levels = trait))     # fixes A->B->C display order

analysis_traits <- as.character(trait_meta$trait)

# --- analysis frames --------------------------------------------------------
design_cols <- c("entry", "loc", "env", "treatment", "p_input",
                 "rep", "gen", "female", "male", "hybrid",
                 "field_row", "field_col")

plots <- field_multiloc %>%                          # every plot, incl. the check
  select(all_of(design_cols), all_of(analysis_traits))

hybrids <- plots %>%                                 # 49 line x tester hybrids only
  filter(gen != "0") %>%
  mutate(env = droplevels(env), female = droplevels(female),
         male = droplevels(male), gen = droplevels(gen),
         row = factor(field_row, levels = 1:10),
         col = factor(field_col, levels = 1:5),
         cross = factor(hybrid))

# --- Tukey IQR outlier flag (Joe's data-quality rule: boxplot / IQR) ---------
# fence = "inner" (1.5*IQR, boxplot whiskers) or "outer" (3*IQR, extreme).
iqr_outliers <- function(x, k = 1.5) {
  q <- quantile(x, c(.25, .75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  !is.na(x) & (x < q[1] - k * iqr | x > q[2] + k * iqr)
}
