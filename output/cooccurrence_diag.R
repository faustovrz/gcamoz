suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

set.seed(20260616)

df <- read_csv("data/multilocation.csv", show_col_types = FALSE) %>%
  rename(ENV = env, REP = rep, GEN = gen) %>%   # env = location x treatment
  mutate(file_order = row_number())             # entry is damaged; use file order

field <- df %>%
  group_by(ENV, REP) %>%
  arrange(file_order, .by_group = TRUE) %>%
  mutate(
    plot_order = row_number(),
    field_row  = ceiling(plot_order / 5),
    pos_in_row = ((plot_order - 1) %% 5) + 1,
    field_col  = if_else(field_row %% 2 == 1, pos_in_row, 6L - pos_in_row)
  ) %>%
  ungroup() %>%
  mutate(GEN = as.integer(GEN))

stopifnot(nrow(field) == 750)   # 5 environments x 3 reps x 50 entries
gens <- sort(unique(field$GEN))
stopifnot(length(gens) == 50)
cat("N genotypes:", length(gens), " range:", min(gens), "-", max(gens), "\n")

# all 1225 unordered pairs
pairs <- t(combn(gens, 2))
pair_key <- function(i, j) paste(pmin(i, j), pmax(i, j), sep = "_")
all_keys <- pair_key(pairs[, 1], pairs[, 2])
stopifnot(length(all_keys) == 1225)

# count co-occurrences within blocks for a given set of (block-id, GEN) rows
cooc_counts <- function(data, block_cols) {
  data <- data %>% mutate(.blk = interaction(across(all_of(block_cols)), drop = TRUE))
  # for each block, all pairs of GEN within it
  cnt <- data %>%
    group_by(.blk) %>%
    summarise(pk = list({
      g <- sort(GEN)
      if (length(g) < 2) character(0) else {
        cb <- combn(g, 2)
        paste(cb[1, ], cb[2, ], sep = "_")
      }
    }), .groups = "drop") %>%
    tidyr::unnest(pk) %>%
    count(pk, name = "cooc")
  # join to full pair list, fill 0
  full <- tibble(pk = all_keys) %>%
    left_join(cnt, by = "pk") %>%
    mutate(cooc = replace_na(cooc, 0L))
  full$cooc
}

dist_table <- function(counts, maxc) {
  tab <- table(factor(counts, levels = 0:maxc))
  as.integer(tab)
}

# (a) per environment, blocks = ENV,REP,field_row, counting across the 3 reps of that ENV
envs <- sort(unique(field$ENV))
maxc_overall <- 0
per_env <- lapply(envs, function(e) {
  sub <- field %>% filter(ENV == e)
  cooc_counts(sub, c("REP", "field_row"))
})
names(per_env) <- envs
maxc_env <- max(sapply(per_env, max))

# (b) pooled across all 15 replicates (5 env x 3 reps)
pooled <- cooc_counts(field, c("ENV", "REP", "field_row"))
maxc_pool <- max(pooled)

maxc <- max(maxc_env, maxc_pool)

cat("\n=== (a) PER-ENVIRONMENT distribution (3 reps each) ===\n")
env_tab <- sapply(per_env, dist_table, maxc = maxc)
rownames(env_tab) <- paste0("count=", 0:maxc)
print(env_tab)

cat("\n=== (b) POOLED (15 reps) distribution ===\n")
pool_tab <- dist_table(pooled, maxc = maxc)
names(pool_tab) <- paste0("count=", 0:maxc)
print(pool_tab)

cat("\nMax co-occurrence (per-env):", maxc_env,
    " | (pooled):", maxc_pool, "\n")
cat("# pairs exceeding 1 (per-env, any ENV):",
    sum(sapply(per_env, function(x) sum(x > 1))), "\n")
cat("# pairs exceeding 1 (pooled):", sum(pooled > 1), "\n")

# (c) simulated random RCBD: per-environment (3 reps), permute 50 gens into 10 rows of 5
n_sim <- 1000
sim_maxc <- maxc
sim_accum <- matrix(0, nrow = sim_maxc + 1, ncol = 1)
# template: one ENV with 3 reps x 10 rows x 5 = 150 plots
sim_block <- expand.grid(REP = 1:3, field_row = 1:10, slot = 1:5)
for (s in seq_len(n_sim)) {
  d <- sim_block
  d$GEN <- unlist(lapply(1:3, function(r) sample(gens, 50)))
  cc <- cooc_counts(d %>% mutate(across(c(REP, field_row), as.integer)),
                    c("REP", "field_row"))
  t <- table(factor(cc, levels = 0:sim_maxc))
  sim_accum <- sim_accum + as.integer(t)
}
sim_mean <- sim_accum / n_sim
rownames(sim_mean) <- paste0("count=", 0:sim_maxc)
cat("\n=== (c) SIMULATED RANDOM RCBD per-env expectation (mean # pairs, n_sim=", n_sim, ") ===\n", sep = "")
print(round(sim_mean[, 1], 2))

# ---- PLOT ----
plot_df <- bind_rows(lapply(envs, function(e) {
  tibble(count = 0:maxc, pairs = dist_table(per_env[[e]], maxc),
         panel = e, type = "Observed")
}))
sim_df <- tibble(count = 0:maxc, pairs = sim_mean[, 1],
                 panel = "Random RCBD (sim mean)", type = "RCBD")

# (d) ideal alpha(0,1)-lattice expectation (deterministic)
#   v=50, k=5, s=10, r=3 reps/env -> within-block pairs = r*s*choose(k,2) = 300
#   resolvable alpha(0,1): no pair concurs >1, so count=1 -> 300, count=0 -> 925
alpha_pairs <- numeric(maxc + 1)
alpha_pairs[1] <- 925  # count = 0
alpha_pairs[2] <- 300  # count = 1  (counts >=2 stay 0)
alpha_df <- tibble(count = 0:maxc, pairs = alpha_pairs,
                   panel = "Alpha(0,1) lattice (ideal)", type = "Alpha")

plot_all <- bind_rows(plot_df, sim_df, alpha_df)
# 5 env panels + ideal alpha + random-RCBD expectation, in a 3-column grid
plot_all$panel <- factor(plot_all$panel,
                         levels = c(envs, "Alpha(0,1) lattice (ideal)",
                                    "Random RCBD (sim mean)"))
plot_all$type <- factor(plot_all$type, levels = c("Observed", "Alpha", "RCBD"))

p <- ggplot(plot_all, aes(x = factor(count), y = pairs, fill = type)) +
  geom_col() +
  geom_text(aes(label = round(pairs)), vjust = -0.3, size = 2.6) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c(Observed = "#2c7fb8", Alpha = "#31a354",
                               RCBD = "#d95f0e")) +
  labs(
    title = "Genotype pair co-occurrence within field rows (putative k=5 incomplete blocks)",
    subtitle = "Per-environment over 3 reps; bottom-right = random-RCBD expectation",
    x = "Co-occurrence count (# blocks a pair shares)",
    y = "Number of genotype pairs (of 1225)",
    fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/cooccurrence_distribution.png", p,
       width = 10, height = 6, dpi = 150)
cat("\nPNG saved.\n")