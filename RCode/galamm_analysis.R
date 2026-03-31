# ==============================================================================
# galamm Factor Models for Critical Appraisal Quality Assessment
# Systematic Map: AI use in Ecology
#
# Model: 18 binary/ordinal critical appraisal items (Yes/Partially/No/NA)
# scored across 72 secondary reviews are modelled as indicators of latent
# quality constructs via Generalised Additive Latent and Mixed Models (galamm).
#
# Factor structure (5 theoretically motivated factors):
#   F1 Transparency  : protocol, PRISMA diagram, included/excluded lists
#   F2 Search quality: search string, comprehensiveness, languages reported
#   F3 Methods       : screening methodology, full-text counts, data extraction
#   F4 Open science  : raw data, metadata, supplementary files, code
#   F5 Inclusivity   : non-English consideration, competing interests, intended
#                      non-English inclusion
#
# Strategy:
#   1. Fit a single-factor model (general quality index) as a sanity check
#   2. Fit 5 separate single-factor models, one per domain (robust)
#   3. Attempt a joint 5-factor confirmatory model
#   4. Extract factor scores and visualise by study characteristics
# ==============================================================================

library(galamm)
library(tidyverse)
library(readxl)
library(ggplot2)

# ==============================================================================
# 1. Load data
# ==============================================================================

crit_raw <- read_excel(
  "DataExtracted/Critical appraisal data extraction form (Responses) March 12.xlsx"
)

mapping_raw <- read_excel(
  "DataExtracted/Systematic mapping data extraction form (Responses).xlsx"
)

# Standardise the study ID column name
names(crit_raw)[3]    <- "study_id"
names(mapping_raw)[3] <- "study_id"

crit_raw$study_id    <- trimws(crit_raw$study_id)
mapping_raw$study_id <- trimws(mapping_raw$study_id)

cat("Critical appraisal data: ", nrow(crit_raw), "studies,",
    ncol(crit_raw), "columns\n")
cat("Systematic mapping data: ", nrow(mapping_raw), "studies,",
    ncol(mapping_raw), "columns\n")

# ==============================================================================
# 2. Recode and reshape critical appraisal items
#
# Column 17 (Non-English search terms) is excluded: all responses are "No" or
# "NA" (not applicable), giving essentially zero variance — this item cannot
# contribute to any latent construct.
#
# Remaining 17 items (Excel columns 5–16, 18–22):
#   item01  col 5  : Protocol / citation provided          → F1
#   item02  col 6  : Boolean search string                 → F2
#   item03  col 7  : Non-English consideration             → F5
#   item04  col 8  : Search comprehensiveness assessed     → F2
#   item05  col 9  : Screening methodology described       → F3
#   item06  col 10 : Full-text inclusion counts provided   → F3
#   item07  col 11 : Raw data provided                     → F4
#   item08  col 12 : Metadata provided                     → F4
#   item09  col 13 : Competing interests reported          → F5
#   item10  col 14 : Supplementary files provided          → F4
#   item11  col 15 : Languages searched reported           → F2
#   item12  col 16 : Intended non-English inclusion        → F5
#   item13  col 18 : Data extraction process described     → F3
#   item14  col 19 : PRISMA / ROSES diagram                → F1
#   item15  col 20 : Included reviews listed               → F1
#   item16  col 21 : Excluded reviews with reasons         → F1
#   item17  col 22 : Analysis code / scripts provided      → F4
# ==============================================================================

# Column indices to keep (0-based from crit_raw; col 17 = index 17 is skipped)
item_col_idx <- c(5:16, 18:22)

item_labels <- c(
  "item01_protocol",        # col 5
  "item02_search_string",   # col 6
  "item03_nonenglish",      # col 7
  "item04_search_comp",     # col 8
  "item05_screening",       # col 9
  "item06_fulltext_n",      # col 10
  "item07_raw_data",        # col 11
  "item08_metadata",        # col 12
  "item09_competing",       # col 13
  "item10_suppl_files",     # col 14
  "item11_languages",       # col 15
  "item12_nonenglish_incl", # col 16
  "item13_extraction",      # col 18 (col 17 skipped)
  "item14_prisma",          # col 19
  "item15_incl_list",       # col 20
  "item16_excl_list",       # col 21
  "item17_code"             # col 22
)

# Which factor each item belongs to (same order as item_labels)
item_factor <- c(
  "F1", "F2", "F5", "F2", "F3", "F3",  # items 01–06
  "F4", "F4", "F5", "F4", "F2", "F5",  # items 07–12
  "F3", "F1", "F1", "F1", "F4"          # items 13–17
)

# Extract and rename
crit_items <- crit_raw[, c(3, item_col_idx)]
names(crit_items) <- c("study_id", item_labels)

# Recode: Yes/Partially = 1 (lenient), No = 0, "NA" string / true NA = missing
recode_response <- function(x) {
  case_when(
    x == "Yes"       ~ 1L,
    x == "Partially" ~ 1L,
    x == "No"        ~ 0L,
    TRUE             ~ NA_integer_   # "NA" string and true NA both → missing
  )
}

crit_coded <- crit_items |>
  mutate(across(all_of(item_labels), recode_response))

# Reshape to long format
long_data <- crit_coded |>
  pivot_longer(
    cols      = all_of(item_labels),
    names_to  = "item",
    values_to = "response"
  ) |>
  filter(!is.na(response)) |>
  # item must be an ordered factor — levels define row order in lambda
  mutate(item = factor(item, levels = item_labels))

cat("\nLong data:", nrow(long_data), "observations from",
    n_distinct(long_data$study_id), "studies\n")
cat("Mean pass rate per item:\n")
long_data |>
  group_by(item) |>
  summarise(pass_rate = mean(response), n = n()) |>
  print(n = 20)

# ==============================================================================
# 3. Lambda (loading) matrices
# ==============================================================================

n_items <- length(item_labels)  # 17

# ── 3a. Single-factor lambda (general quality) ──────────────────────────────
# Reference item = item05_screening (49% pass rate) — most moderate, best
# discriminating item. All other loadings are free (NA).
lambda_1f <- matrix(NA, nrow = n_items, ncol = 1)
lambda_1f[5, 1] <- 1   # item05_screening ← reference

# ── 3b. Five-factor lambda (confirmatory) ────────────────────────────────────
# Rows = items (in item_labels order), columns = factors (F1–F5).
# 0 = item does not load on this factor (constrained).
# 1 = fixed reference loading (first item in each factor's block).
# NA = free loading to estimate.
lambda_5f <- matrix(0, nrow = n_items, ncol = 5,
                    dimnames = list(item_labels,
                                   c("F1_transparency", "F2_search",
                                     "F3_methods", "F4_openscience",
                                     "F5_inclusivity")))

# Reference items are chosen to have the most moderate pass rates within each
# domain — extreme pass rates (< 5% or > 85%) lead to Hauck-Donner separation
# in binary models and must be avoided as references.
#
# Pass rates (approximate):
#   item01 1%, item02 96%, item03 3%, item04 7%,  item05 49%, item06 79%
#   item07 21%, item08 7%, item09 82%, item10 39%, item11 32%, item12 3%
#   item13 17%, item14 35%, item15 39%, item16 4%, item17 3%

# F1: Transparency (reference = item15_incl_list, 39% pass rate)
lambda_5f[1,  1] <- NA   # item01_protocol          free (1%)
lambda_5f[14, 1] <- NA   # item14_prisma             free (35%)
lambda_5f[15, 1] <- 1    # item15_incl_list          ← reference (39%)
lambda_5f[16, 1] <- NA   # item16_excl_list          free (4%)

# F2: Search quality (reference = item11_languages, 32% pass rate)
lambda_5f[2,  2] <- NA   # item02_search_string      free (96%)
lambda_5f[4,  2] <- NA   # item04_search_comp        free (7%)
lambda_5f[11, 2] <- 1    # item11_languages          ← reference (32%)

# F3: Methods reporting (reference = item05_screening, 49% pass rate)
lambda_5f[5,  3] <- 1    # item05_screening          ← reference (49%)
lambda_5f[6,  3] <- NA   # item06_fulltext_n         free (79%)
lambda_5f[13, 3] <- NA   # item13_extraction         free (17%)

# F4: Open science (reference = item10_suppl_files, 39% pass rate)
lambda_5f[7,  4] <- NA   # item07_raw_data           free (21%)
lambda_5f[8,  4] <- NA   # item08_metadata           free (7%)
lambda_5f[10, 4] <- 1    # item10_suppl_files        ← reference (39%)
lambda_5f[17, 4] <- NA   # item17_code               free (3%)

# F5: Inclusivity (reference = item09_competing, 82% — least extreme option)
# Note: item03 (3%) and item12 (3%) have extreme pass rates; model estimates
# for these items should be interpreted cautiously.
lambda_5f[3,  5] <- NA   # item03_nonenglish         free (3%)
lambda_5f[9,  5] <- 1    # item09_competing          ← reference (82%)
lambda_5f[12, 5] <- NA   # item12_nonenglish_incl    free (3%)

cat("\n=== Lambda matrix (17 items × 5 factors) ===\n")
print(lambda_5f)

# ── 3c. Domain-specific lambda matrices (one per factor) ────────────────────
domains <- list(
  F1_transparency = which(item_factor == "F1"),
  F2_search       = which(item_factor == "F2"),
  F3_methods      = which(item_factor == "F3"),
  F4_openscience  = which(item_factor == "F4"),
  F5_inclusivity  = which(item_factor == "F5")
)

# ==============================================================================
# 4. Single-factor model — general quality index
# ==============================================================================

cat("\n=== MODEL 1: Single-factor (general quality) ===\n")

mod_1f <- galamm(
  formula  = response ~ 0 + item + (0 + eta | study_id),
  data     = long_data,
  family   = binomial,
  load_var = "item",
  lambda   = lambda_1f,
  factor   = "eta"
)

cat("Converged:", mod_1f$model$opt$convergence == 0, "\n")
print(summary(mod_1f))

# Extract factor loadings from the 1-factor model summary
# summary()$Lambda is a data.frame with columns eta and SE, rows = lambda1..N
lam_summary_1f  <- summary(mod_1f)$Lambda
loadings_1f <- data.frame(
  item       = item_labels,
  factor     = item_factor,
  pass_rate  = sapply(item_labels,
                 function(i) mean(long_data$response[long_data$item == i])),
  loading    = lam_summary_1f[, 1],  # estimated lambda values
  loading_se = lam_summary_1f[, 2]   # SE (NaN for items with separation)
)
cat("\nEstimated loadings (1-factor model):\n")
print(loadings_1f[order(loadings_1f$factor, -loadings_1f$loading), ],
      row.names = FALSE)

# ==============================================================================
# 5. Five separate domain models (one factor per domain)
# ==============================================================================

cat("\n=== MODEL SET: Domain-specific single-factor models ===\n")

domain_models  <- list()
domain_scores  <- list()

for (dom_name in names(domains)) {

  dom_items <- item_labels[domains[[dom_name]]]
  dat_dom   <- long_data |> filter(item %in% dom_items) |>
    mutate(item = factor(item, levels = dom_items))

  n_dom <- length(dom_items)

  # Choose reference item: the item with pass rate closest to 0.5
  pass_rates  <- sapply(dom_items,
    function(i) mean(dat_dom$response[dat_dom$item == i]))
  ref_idx     <- which.min(abs(pass_rates - 0.5))  # most moderate

  lam_dom <- matrix(NA, nrow = n_dom, ncol = 1)
  lam_dom[ref_idx, 1] <- 1   # fix reference

  cat("\n-- Domain:", dom_name, "(", n_dom, "items ) --\n")
  cat("  Reference item:", dom_items[ref_idx],
      sprintf("(pass rate = %.1f%%)\n", 100 * pass_rates[ref_idx]))
  cat("  Pass rates:", paste(sprintf("%s=%.0f%%", dom_items,
                                     100 * pass_rates), collapse = ", "), "\n")

  # Try binomial first; fall back to gaussian if optimizer fails
  fit_domain <- function(fam) {
    tryCatch(
      galamm(
        formula  = response ~ 0 + item + (0 + eta | study_id),
        data     = dat_dom,
        family   = fam,
        load_var = "item",
        lambda   = lam_dom,
        factor   = "eta"
      ),
      error = function(e) {
        cat("  [", deparse(substitute(fam)), "] ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )
  }

  mod_dom <- fit_domain(binomial)
  if (is.null(mod_dom)) {
    cat("  Retrying with Gaussian family...\n")
    mod_dom <- fit_domain(gaussian)
  }

  if (!is.null(mod_dom)) {
    cat("  Converged:", mod_dom$model$opt$convergence == 0, "\n")
    cat("  Family:", mod_dom$model$lmod$resp$family$family, "\n")
    print(summary(mod_dom))

    scores_dom <- ranef(mod_dom)$study_id |>
      rownames_to_column("study_id") |>
      rename(!!dom_name := eta)

    domain_models[[dom_name]] <- mod_dom
    domain_scores[[dom_name]] <- scores_dom
  } else {
    cat("  Both families failed — domain omitted.\n")
  }
}

# ==============================================================================
# 6. Joint 5-factor confirmatory model
# ==============================================================================

cat("\n=== MODEL 2: Joint 5-factor confirmatory model ===\n")

mod_5f <- tryCatch(
  galamm(
    formula  = response ~ 0 + item +
      (0 + eta1 + eta2 + eta3 + eta4 + eta5 | study_id),
    data     = long_data,
    family   = binomial,
    load_var = "item",
    lambda   = lambda_5f,
    factor   = c("eta1", "eta2", "eta3", "eta4", "eta5")
  ),
  error = function(e) {
    cat("5-factor model failed:", conditionMessage(e), "\n")
    cat("Using domain-specific factor scores instead.\n")
    NULL
  }
)

if (!is.null(mod_5f)) {
  cat("Converged:", mod_5f$model$opt$convergence == 0, "\n")
  print(summary(mod_5f))
}

# ==============================================================================
# 7. Extract factor scores and join with mapping data
# ==============================================================================

# Use 5-factor model scores if available, otherwise use domain model scores
if (!is.null(mod_5f)) {
  scores_wide <- ranef(mod_5f)$study_id |>
    rownames_to_column("study_id") |>
    rename(
      F1_transparency = eta1,
      F2_search       = eta2,
      F3_methods      = eta3,
      F4_openscience  = eta4,
      F5_inclusivity  = eta5
    )
  score_source <- "joint 5-factor model"
} else {
  # Merge scores from the 5 separate domain models
  scores_wide <- Reduce(
    function(a, b) left_join(a, b, by = "study_id"),
    domain_scores
  )
  score_source <- "5 separate domain models"
}

cat("\nFactor scores extracted from:", score_source, "\n")

# Also add general quality index from 1-factor model
scores_1f_df <- ranef(mod_1f)$study_id |>
  rownames_to_column("study_id") |>
  rename(general_quality = eta)

scores_wide <- left_join(scores_wide, scores_1f_df, by = "study_id")

# Subset of mapping variables as study-level predictors
mapping_sub <- mapping_raw |>
  select(
    study_id,
    year          = `Year of publication`,
    review_type   = `The main type of secondary review`,
    ai_category   = `Broad category of AI models discussed`,
    data_modality = `Categories of the data used for the AI models`,
    n_articles    = `Number of articles synthesized`
  ) |>
  mutate(
    study_id   = trimws(study_id),
    year       = as.integer(year),
    n_articles = as.numeric(n_articles)
  )

scores_full <- left_join(scores_wide, mapping_sub, by = "study_id")

cat("\n=== Factor scores (first 10 studies) ===\n")
print(head(scores_full, 10))

# ==============================================================================
# 8. Visualise
# ==============================================================================

factor_nice <- c(
  F1_transparency = "F1: Transparency",
  F2_search       = "F2: Search Quality",
  F3_methods      = "F3: Methods",
  F4_openscience  = "F4: Open Science",
  F5_inclusivity  = "F5: Inclusivity",
  general_quality = "General Quality"
)

scores_long <- scores_full |>
  pivot_longer(
    cols      = c(starts_with("F"), general_quality),
    names_to  = "factor",
    values_to = "score"
  ) |>
  mutate(factor = recode(factor, !!!factor_nice))

# ── Plot 1: factor scores by publication year ─────────────────────────────────
p1 <- scores_long |>
  filter(!is.na(year), factor != "General Quality") |>
  ggplot(aes(x = factor(year), y = score, fill = factor)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  facet_wrap(~factor, scales = "free_y", ncol = 2) +
  labs(
    title   = "Latent quality factor scores by publication year",
    x       = "Year",
    y       = "Factor score (latent scale)",
    caption = paste("Factor scores from:", score_source)
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"))

ggsave("RCode/factor_scores_by_year.png", p1,
       width = 10, height = 10, dpi = 150)
cat("Saved: RCode/factor_scores_by_year.png\n")

# ── Plot 2: factor loading heatmap (from 1-factor model) ─────────────────────
# Cap very large loadings (separation artefacts) for display purposes
loading_long <- loadings_1f |>
  mutate(
    item         = factor(item, levels = rev(item_labels)),
    loading_disp = pmin(pmax(loading, -3), 3),   # cap at ±3 for colour scale
    label        = ifelse(abs(loading) > 3,
                          paste0(round(loading, 1), "*"),
                          as.character(round(loading, 2)))
  )

p2 <- ggplot(loading_long, aes(x = "General\nQuality", y = item,
                                fill = loading_disp)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 0, limits = c(-3, 3), name = "Loading\n(capped)") +
  labs(
    title    = "Factor loadings — single-factor (general quality) model",
    subtitle = "* = loading capped for display; large values indicate near-separation",
    x        = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0))

ggsave("RCode/factor_loadings_1f.png", p2,
       width = 5, height = 8, dpi = 150)
cat("Saved: RCode/factor_loadings_1f.png\n")

# ── Plot 3: general quality score vs number of articles synthesised ───────────
p3 <- scores_full |>
  filter(!is.na(n_articles)) |>
  ggplot(aes(x = n_articles, y = general_quality,
             colour = review_type, shape = review_type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  scale_x_log10() +
  labs(
    title  = "General quality score vs number of articles synthesised",
    x      = "Number of articles synthesised (log scale)",
    y      = "General quality factor score",
    colour = "Review type", shape = "Review type"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("RCode/quality_vs_n_articles.png", p3,
       width = 8, height = 6, dpi = 150)
cat("Saved: RCode/quality_vs_n_articles.png\n")

# ==============================================================================
# 9. Restructured 5-factor model
#
# Changes from original:
#   • item09_competing (ethics/transparency) moved from F5 → F1
#   • F5 (Language inclusivity) retains only items 03 and 12
#     But both have ~3% pass rate — variance is near zero, so F5 is dropped
#     from the joint model and replaced by a single "inclusivity pass" indicator
#   • item02_search_string (96% pass) removed from F2 (no discrimination)
#   → 4-factor restructured model: F1*, F2*, F3, F4
# ==============================================================================

cat("\n=== MODEL 3: Restructured 4-factor model ===\n")
cat("Changes: item09→F1 (ethics/transparency); item02 dropped (96% pass rate);\n")
cat("         F5 dropped (items 03+12 both 3% pass, no between-study variance)\n\n")

# Items retained and their factor assignments
item_labels_r <- c(
  "item01_protocol",        # F1
  "item04_search_comp",     # F2*
  "item05_screening",       # F3
  "item06_fulltext_n",      # F3
  "item07_raw_data",        # F4
  "item08_metadata",        # F4
  "item09_competing",       # F1* (moved from F5)
  "item10_suppl_files",     # F4
  "item11_languages",       # F2*
  "item13_extraction",      # F3
  "item14_prisma",          # F1
  "item15_incl_list",       # F1
  "item16_excl_list",       # F1
  "item17_code"             # F4
)
item_factor_r <- c("F1","F2","F3","F3","F4","F4","F1","F4","F2","F3","F1","F1","F1","F4")

long_r <- long_data |>
  filter(item %in% item_labels_r) |>
  mutate(item = factor(item, levels = item_labels_r))

n_r <- length(item_labels_r)  # 14

# Lambda matrix (14 × 4); references chosen at moderate pass rates
lambda_r <- matrix(0, nrow = n_r, ncol = 4,
                   dimnames = list(item_labels_r,
                                   c("F1_transp", "F2_search",
                                     "F3_methods", "F4_openscience")))

# F1* Transparency+Ethics  (ref = item15_incl_list, 39%)
lambda_r[1,  1] <- NA  # item01_protocol       (1%)
lambda_r[7,  1] <- NA  # item09_competing      (82%)
lambda_r[11, 1] <- NA  # item14_prisma         (35%)
lambda_r[12, 1] <- 1   # item15_incl_list      ← reference (39%)
lambda_r[13, 1] <- NA  # item16_excl_list      (4%)

# F2* Search  (ref = item11_languages, 32%)
lambda_r[2,  2] <- NA  # item04_search_comp    (7%)
lambda_r[9,  2] <- 1   # item11_languages      ← reference (32%)

# F3 Methods  (ref = item05_screening, 49%)
lambda_r[3,  3] <- 1   # item05_screening      ← reference (49%)
lambda_r[4,  3] <- NA  # item06_fulltext_n     (79%)
lambda_r[10, 3] <- NA  # item13_extraction     (17%)

# F4 Open science  (ref = item10_suppl_files, 39%)
lambda_r[5,  4] <- NA  # item07_raw_data       (21%)
lambda_r[6,  4] <- NA  # item08_metadata       (7%)
lambda_r[8,  4] <- 1   # item10_suppl_files    ← reference (39%)
lambda_r[14, 4] <- NA  # item17_code           (3%)

cat("Restructured lambda (14 items × 4 factors):\n")
print(lambda_r)

mod_r <- tryCatch(
  galamm(
    formula  = response ~ 0 + item +
      (0 + eta1 + eta2 + eta3 + eta4 | study_id),
    data     = long_r,
    family   = binomial,
    load_var = "item",
    lambda   = lambda_r,
    factor   = c("eta1", "eta2", "eta3", "eta4")
  ),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(mod_r)) {
  cat("Converged:", mod_r$model$opt$convergence == 0, "\n")
  print(summary(mod_r))

  scores_r <- ranef(mod_r)$study_id |>
    rownames_to_column("study_id") |>
    rename(F1_transp_eth = eta1, F2_search2 = eta2,
           F3_methods2   = eta3, F4_openscience2 = eta4)
  scores_full <- left_join(scores_full, scores_r, by = "study_id")

  # ── Plot 4: restructured 4-factor loading heatmap ──────────────────────────
  lam_r_summ <- summary(mod_r)$Lambda
  lam_r_df   <- data.frame(
    item    = item_labels_r[rep(seq_len(n_r), 4)],
    factor  = rep(c("F1*","F2*","F3","F4"), each = n_r),
    loading = as.numeric(lambda_r)
  )
  # Replace NAs with estimated values from summary
  est_vals <- lam_r_summ[, 1]
  k <- 1
  for (j in seq_len(4)) {
    for (i in seq_len(n_r)) {
      if (is.na(lambda_r[i, j])) {
        lam_r_df$loading[lam_r_df$item == item_labels_r[i] &
                         lam_r_df$factor == c("F1*","F2*","F3","F4")[j]] <- est_vals[k]
        k <- k + 1
      }
    }
  }
  lam_r_df <- lam_r_df |>
    filter(loading != 0) |>
    mutate(
      item    = factor(item, levels = rev(item_labels_r)),
      loading_disp = pmin(pmax(loading, -5), 5),
      label   = ifelse(abs(loading) > 5,
                       paste0(round(loading, 1), "*"),
                       as.character(round(loading, 2)))
    )

  p4 <- ggplot(lam_r_df, aes(x = factor, y = item, fill = loading_disp)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 3.2) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                         midpoint = 0, limits = c(-5, 5),
                         name = "Loading\n(capped ±5)") +
    labs(
      title    = "Factor loadings — restructured 4-factor model",
      subtitle = "F1* = Transparency+Ethics; F2* = Search (item02 removed)",
      x = "Factor", y = NULL
    ) +
    theme_bw(base_size = 11)

  ggsave("RCode/factor_loadings_4f.png", p4,
         width = 7, height = 7, dpi = 150)
  cat("Saved: RCode/factor_loadings_4f.png\n")
}

# ==============================================================================
# 10. Model with study-level predictors
#
# Add year and review type as fixed effects to the 1-factor model.
# These are study-level covariates; we join them to the long data so that
# each item observation inherits the study's predictors.
# ==============================================================================

cat("\n=== MODEL 4: General quality ~ year + review type ===\n")

long_pred <- long_data |>
  left_join(mapping_sub, by = "study_id") |>
  mutate(
    year_c       = scale(year)[, 1],          # centre and scale year
    log_n        = scale(log1p(n_articles))[, 1],
    is_sys_map   = as.integer(review_type == "Systematic map")
  )

mod_pred <- tryCatch(
  galamm(
    formula  = response ~ 0 + item + year_c + is_sys_map +
                          (0 + eta | study_id),
    data     = long_pred,
    family   = binomial,
    load_var = "item",
    lambda   = lambda_1f,
    factor   = "eta"
  ),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(mod_pred)) {
  cat("Converged:", mod_pred$model$opt$convergence == 0, "\n")
  print(summary(mod_pred))
  cat("\nInterpretation:\n")
  cat("  year_c    : effect of publication year on item pass rate (log-odds)\n")
  cat("  is_sys_map: difference between systematic maps and systematic reviews\n")
}

# ==============================================================================
# 11. Summary loading table (print-ready)
# ==============================================================================

cat("\n=== SUMMARY: Factor loadings from 1-factor model ===\n")
loadings_1f |>
  arrange(factor, desc(loading)) |>
  mutate(
    loading    = round(loading, 3),
    loading_se = round(loading_se, 3),
    pass_pct   = sprintf("%.1f%%", 100 * pass_rate),
    flag       = case_when(
      pass_rate < 0.05 | pass_rate > 0.90 ~ "extreme",
      abs(loading) > 5                    ~ "separation",
      TRUE                                ~ ""
    )
  ) |>
  select(factor, item, pass_pct, loading, loading_se, flag) |>
  print(row.names = FALSE)

# ==============================================================================
# 12. Save all results
# ==============================================================================

write_csv(scores_full, "RCode/factor_scores.csv")
cat("\nFactor scores saved to: RCode/factor_scores.csv\n")

write_csv(loadings_1f, "RCode/loadings_1factor.csv")
cat("Loadings (1-factor) saved to: RCode/loadings_1factor.csv\n")

cat("\n=== Done ===\n")
