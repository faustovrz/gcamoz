suppressPackageStartupMessages({
  library(FielDHub); library(readr); library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})

## ---- STEP 1: build the design exactly as the notebook does ----
alpha_design <- alpha_lattice(t = 50, k = 5, r = 3, l = 1,
                              plotNumber = 101, locationNames = "ALPHA-DEMO", seed = 1234)
field_book <- as_tibble(alpha_design$fieldBook)

cat("=== STEP 1: FielDHub fieldBook columns ===\n")
print(names(field_book))
cat("\nN plots:", nrow(field_book), "\n")
cat("N reps:", length(unique(field_book$REP)), "\n")
cat("N entries:", length(unique(field_book$ENTRY)), " range:",
    min(field_book$ENTRY), "-", max(field_book$ENTRY), "\n")
cat("IBLOCK per rep:\n"); print(field_book %>% distinct(REP, IBLOCK) %>% count(REP, name = "n_iblock"))
cat("Plots per (REP,IBLOCK):\n"); print(field_book %>% count(REP, IBLOCK) %>% count(n, name = "n_blocks"))

## ---- shared co-occurrence machinery ----
all_keys_for <- function(gens) {
  cb <- combn(sort(gens), 2)
  paste(cb[1, ], cb[2, ], sep = "_")
}
cooc_counts <- function(data, block_cols, gen_col, all_keys) {
  data <- data %>% mutate(.blk = interaction(across(all_of(block_cols)), drop = TRUE),
                          .g = .data[[gen_col]])
  cnt <- data %>%
    group_by(.blk) %>%
    summarise(pk = list({
      g <- sort(.g)
      if (length(g) < 2) character(0) else {
        cb <- combn(g, 2); paste(cb[1, ], cb[2, ], sep = "_")
      }
    }), .groups = "drop") %>%
    tidyr::unnest(pk) %>%
    count(pk, name = "cooc")
  full <- tibble(pk = all_keys) %>%
    left_join(cnt, by = "pk") %>%
    mutate(cooc = replace_na(cooc, 0L))
  full$cooc
}
dist_table <- function(counts, maxc) as.integer(table(factor(counts, levels = 0:maxc)))

## ---- STEP 2: co-occurrences in FielDHub TRUE incomplete blocks ----
fh_gens <- sort(unique(field_book$ENTRY))
stopifnot(length(fh_gens) == 50)
all_keys <- all_keys_for(fh_gens)
stopifnot(length(all_keys) == 1225)

fh_cooc <- cooc_counts(field_book, c("REP", "IBLOCK"), "ENTRY", all_keys)
fh_max <- max(fh_cooc)
cat("\n=== STEP 2: FielDHub alpha lattice co-occurrence distribution ===\n")
fh_full_tab <- dist_table(fh_cooc, fh_max)
names(fh_full_tab) <- paste0("count=", 0:fh_max)
print(fh_full_tab)
cat("Max co-occurrence:", fh_max, "\n")
cat("# pairs exceeding 1:", sum(fh_cooc > 1), "\n")

## ---- STEP 3: Stelio CHOKWE STS ----
df <- read_csv("data/multilocation.csv", show_col_types = FALSE) %>%
  rename(ENV = loc, REP = rep, GEN = gen) %>%   # new schema -> internal names
  mutate(file_order = row_number())             # entry is damaged; use file order
cat("\nENV values present:\n"); print(sort(unique(df$ENV)))

chk <- df %>%
  filter(ENV == "CHOKWE STS") %>%
  group_by(REP) %>%
  arrange(file_order, .by_group = TRUE) %>%
  mutate(plot_order = row_number(),
         field_row  = ceiling(plot_order / 5)) %>%
  ungroup() %>%
  mutate(GEN = as.integer(GEN))

chk_gens <- sort(unique(chk$GEN))
cat("\nCHOKWE STS rows:", nrow(chk), " n genotypes:", length(chk_gens), "\n")
chk_keys <- all_keys_for(chk_gens)
chk_cooc <- cooc_counts(chk, c("REP", "field_row"), "GEN", chk_keys)
chk_max <- max(chk_cooc)

## common axis
maxc <- max(fh_max, chk_max, 9)

cat("\n=== STEP 3: CHOKWE STS distribution ===\n")
chk_tab <- dist_table(chk_cooc, maxc); names(chk_tab) <- paste0("count=", 0:maxc)
print(chk_tab)
cat("Max:", chk_max, " | # pairs > 1:", sum(chk_cooc > 1), "\n")

## ---- STEP 4: 3-facet bar chart ----
fh_df <- tibble(count = 0:maxc, pairs = dist_table(fh_cooc, maxc),
                panel = "FielDHub alpha lattice (constructed)", grp = "fh")
ideal <- numeric(maxc + 1); ideal[1] <- 925; ideal[2] <- 300
ideal_df <- tibble(count = 0:maxc, pairs = ideal,
                   panel = "Alpha(0,1) ideal", grp = "ideal")
chk_df <- tibble(count = 0:maxc, pairs = dist_table(chk_cooc, maxc),
                 panel = "Stelio CHOKWE STS (observed)", grp = "chk")

plot_all <- bind_rows(fh_df, ideal_df, chk_df)
plot_all$panel <- factor(plot_all$panel,
                         levels = c("FielDHub alpha lattice (constructed)",
                                    "Alpha(0,1) ideal",
                                    "Stelio CHOKWE STS (observed)"))
plot_all$grp <- factor(plot_all$grp, levels = c("fh", "ideal", "chk"))

p <- ggplot(plot_all, aes(x = factor(count), y = pairs, fill = grp)) +
  geom_col() +
  geom_text(aes(label = round(pairs)), vjust = -0.3, size = 2.6) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c(fh = "#31a354", ideal = "#31a354", chk = "#2c7fb8"),
                    guide = "none") +
  labs(title = "Co-occurrence test on a TRUE alpha lattice vs. Stelio's data",
       x = "Co-occurrence count (# blocks a pair shares)",
       y = "Number of genotype pairs (of 1225)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/fieldhub_alpha_cooccurrence.png", p, width = 11, height = 4.5, dpi = 150)
cat("\nPNG saved to output/fieldhub_alpha_cooccurrence.png\n")
