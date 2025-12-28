# ============================================================
# PANEL DATA ECONOMETRICS - Code File Part 1 (R translation)
# Source: /mnt/data/PDE_codefile_part1.ipynb
# ============================================================

# ---------------------------
# 0) REQUIREMENTS SET-UP
# ---------------------------

# PY (original):
# # Requirements.txt file installation
# # !pip install -r requirements.txt

# R equivalent: install required packages once (uncomment if needed)
# install.packages(c(
#   "readxl","writexl","dplyr","tidyr","lubridate","stringr","ggplot2",
#   "purrr","broom","plm","psych","car","lmtest","sandwich","patchwork"
# ))

suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(broom)
  library(plm)
  library(psych)
  library(car)
  library(lmtest)
  library(sandwich)
})

# ---------------------------
# 1) LIBRARIES / HELPERS
# ---------------------------

# PY (original imports):
# import warnings, pandas as pd, numpy as np, random
# import matplotlib.pyplot as plt, plotly.graph_objects as go, seaborn as sns
# import scipy.stats as stats
# from scipy.stats import norm, levene, ks_2samp, kstest, pearsonr, spearmanr
# import statsmodels.api as sm
# import statsmodels.formula.api as smf
# from statsmodels.nonparametric.smoothers_lowess import lowess
# from linearmodels.panel import PanelOLS, RandomEffects, PooledOLS
#
# R equivalents:
# - warnings: base options / suppressWarnings
# - pandas/numpy: dplyr + base vectors
# - scipy/stats: stats, car, lmtest
# - seaborn/matplotlib: ggplot2
# - statsmodels/linearmodels: lm(), plm(), fixest (optional)

options(warn = -1) # PY: warnings.filterwarnings("ignore") (approx)

# ------------------------------------------------------------
# Statistical Significance labelling
# ------------------------------------------------------------
# PY:
# def significance_stars(p):
#   if p < 0.001: return "***"
#   elif p < 0.01: return "**"
#   elif p < 0.05: return "*"
#   else: return ""

significance_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# ------------------------------------------------------------
# Plot helper: histogram + KDE + normal curve overlay
# ------------------------------------------------------------
# Translation note:
# - Python used seaborn histplot + kdeplot and then plotted a normal curve with same mean/sd.
# - R uses ggplot2: geom_histogram(aes(y=after_stat(density))) + geom_density() + stat_function(dnorm).

overlay_one_gg <- function(data_vec, label = "", bins = 30) {
  data_vec <- data_vec[!is.na(data_vec)]
  n <- length(data_vec)
  if (n < 5) {
    return(ggplot() + ggtitle(paste(label, "(n<5)")) + theme_minimal())
  }
  mu <- mean(data_vec)
  sd_ <- sd(data_vec)
  ggplot(data.frame(x = data_vec), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.25) +
    geom_density(linewidth = 0.8) +
    stat_function(fun = dnorm, args = list(mean = mu, sd = sd_), linewidth = 0.8, linetype = "dashed") +
    labs(title = label, y = "Density", x = NULL) +
    theme_minimal()
}

# ------------------------------------------------------------
# Correlation matrix with t-stats, p-values and significance stars
# ------------------------------------------------------------
# Python computed correlation matrices, and in some places built annotation matrices with t-stats/stars.
# In R we implement a utility that returns:
# - corr matrix
# - p-value matrix
# - annotation matrix (rounded r + stars)

corr_with_p <- function(df_num, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  df_num <- as.data.frame(df_num)
  cols <- names(df_num)
  k <- length(cols)

  r_mat <- matrix(NA_real_, nrow = k, ncol = k, dimnames = list(cols, cols))
  p_mat <- matrix(NA_real_, nrow = k, ncol = k, dimnames = list(cols, cols))

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      xi <- df_num[[i]]
      xj <- df_num[[j]]
      ok <- complete.cases(xi, xj)
      if (sum(ok) >= 3) {
        ct <- suppressWarnings(cor.test(xi[ok], xj[ok], method = method))
        r_mat[i, j] <- unname(ct$estimate)
        p_mat[i, j] <- ct$p.value
      }
    }
  }

  annot <- matrix("", nrow = k, ncol = k, dimnames = list(cols, cols))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (!is.na(r_mat[i, j])) {
        annot[i, j] <- sprintf("%.2f%s", r_mat[i, j], significance_stars(p_mat[i, j]))
      }
    }
  }

  list(r = r_mat, p = p_mat, annot = annot)
}

# ============================================================
# 1) PART 1 - DATASET, UNIVARIATE & BIVARIATE DESCRIPTIVE STATS
# ============================================================

# ---------------------------
# Data loading
# ---------------------------

# PY:
# raw_df = pd.read_excel("Data/Dataset_MP_Impact_functional_Distribution.xlsx")
# raw_df["year"] = pd.to_datetime(raw_df["year"], format="%Y")

raw_df <- read_excel("Data/Dataset_MP_Impact_functional_Distribution.xlsx")

# Ensure expected columns exist
stopifnot(all(c("year", "country") %in% names(raw_df)))

# Convert year to Date (Jan 1 of that year) so it is not treated as numeric in descriptive stats
# Translation note: Python used datetime with format %Y; in R we build a Date.
raw_df <- raw_df %>%
  mutate(year = as.Date(paste0(as.integer(year), "-01-01")))

# ---------------------------
# Df Symbols Legend (labels)
# ---------------------------

# PY: labels_mapper = {...}
labels_mapper <- c(
  year    = "Year",
  country = "Country",
  i       = "Short-Term Interest Rate",
  P       = "GDP Deflator",
  W       = "Nominal Compensation Per Employee",
  WR      = "Real Compensation Per Employee",
  GDP     = "Real Gross Domestic Product",
  LS      = "Adjusted Labour Share",
  PCOM    = "Energy Commodities Price Index",
  UN      = "Unemployment Rate",
  SHORTUN = "Short-Term Unemployment",
  LONGUN  = "Long-Term Unemployment",
  LF      = "Labor Force",
  REER    = "Real Effective Exchange Rate",
  SH      = "Shadow Interest Rate"
)

# ---------------------------
# 1.1 / 1.2: panel shape basics
# ---------------------------

# PY pattern:
# raw_df["country"].nunique(), raw_df["year"].nunique(), min/max years, etc.
n_countries <- n_distinct(raw_df$country)
n_years <- n_distinct(raw_df$year)

year_min <- min(raw_df$year, na.rm = TRUE)
year_max <- max(raw_df$year, na.rm = TRUE)

message("Number of countries: ", n_countries)
message("Number of years: ", n_years)
message("Time span: ", format(year_min), " to ", format(year_max))

# ---------------------------
# 1.3: quick comments / correlation motivation
# ---------------------------
# Translation note: notebook question prompt; no executable code needed.

# ---------------------------
# 1.4 VARIABLE TRANSFORMATIONS: between / within
# ---------------------------

# PY:
# trans_df = raw_df.copy()
# numeric_cols = raw_df.select_dtypes(include="number").columns
# for c in numeric_cols:
#   trans_df[f"{c}_between"] = raw_df.groupby("country")[c].transform("mean")
# for c in numeric_cols:
#   trans_df[f"{c}_within"] = raw_df[c] - raw_df.groupby("country")[c].transform("mean")

trans_df <- raw_df

numeric_cols <- names(raw_df)[vapply(raw_df, is.numeric, logical(1))]
numeric_cols <- setdiff(numeric_cols, c()) # keep as-is (year is Date, not numeric)

# Between transformation (country means repeated per row)
for (c in numeric_cols) {
  trans_df[[paste0(c, "_between")]] <- ave(raw_df[[c]], raw_df$country, FUN = function(x) mean(x, na.rm = TRUE))
}

# Within (one-way FE demeaning by country)
for (c in numeric_cols) {
  mu_i <- ave(raw_df[[c]], raw_df$country, FUN = function(x) mean(x, na.rm = TRUE))
  trans_df[[paste0(c, "_within")]] <- raw_df[[c]] - mu_i
}

# ---------------------------
# 1.5: Plot distributions - within vs between
# ---------------------------

# Python created a grid of subplots and overlaid hist/kde/normal for within and between.
# R approach: create a list of ggplots and (optionally) arrange with patchwork.

target_variables <- c("GDP", "UN", "P", "WR", "LS", "i")

plots_within_between <- map(target_variables, function(var) {
  within_col  <- paste0(var, "_within")
  between_col <- paste0(var, "_between")

  p_between <- overlay_one_gg(trans_df[[between_col]], label = paste0(var, "_between"))
  p_within  <- overlay_one_gg(trans_df[[within_col]],  label = paste0(var, "_within"))

  # Combine using facet-like layout (two plots stacked) without extra packages:
  # Difficulty note: Python overlaid within & between on same axes; ggplot overlay is possible but more fiddly.
  # Here we present them as a pair.
  list(between = p_between, within = p_within)
})

# If you want a single-page arrangement, install patchwork and use:
# library(patchwork)
# wrap_plots(map(plots_within_between, ~ .x$between), ncol = 3) / wrap_plots(map(plots_within_between, ~ .x$within), ncol = 3)

# ---------------------------
# 1.6 FD. FIRST DIFFERENCE TRANSFORMATION
# ---------------------------

# PY pattern:
# trans_diff_df = raw_df.sort_values(["country","year"]).groupby("country")[numeric].diff()
# keep year/country, drop first obs per country, etc.

trans_diff_df <- raw_df %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(across(all_of(numeric_cols), ~ .x - lag(.x), .names = "{.col}_1diff")) %>%
  ungroup()

# ---------------------------
# 1.7 / 1.8: First difference scatter + correlations
# ---------------------------

# Python plotted GDPG vs EDA/GDP after transformations.
# In this dataset, variable names differ; adapt by changing ycol/xcol below.
# Difficulty note: the original notebook likely defined specific variables later (GDPG, EDA_GDP etc.).
# Here we keep a generic example using GDP and i if those exist.

if (all(c("GDP_1diff", "i_1diff") %in% names(trans_diff_df))) {
  df_scatter <- trans_diff_df %>% select(country, year, GDP_1diff, i_1diff)
  ct <- cor.test(df_scatter$GDP_1diff, df_scatter$i_1diff, use = "complete.obs")

  p_fd_scatter <- ggplot(df_scatter, aes(x = i_1diff, y = GDP_1diff)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = sprintf("First Differences: corr=%.3f%s", unname(ct$estimate), significance_stars(ct$p.value)),
      x = "i_1diff", y = "GDP_1diff"
    ) +
    theme_minimal()
}

# ---------------------------
# 1.9 - 1.11: TWFE balanced panel transformations and variance ordering
# ---------------------------

# Python logic (typical):
# - Build balanced panel subset
# - Compute TWFE: x_it - x_i. - x_.t + x_..
# - Compute variances of transformed vars and order

# Balanced panel detection
# A panel is balanced if each country has same number of time observations (all years in sample).
panel_counts <- raw_df %>% count(country, name = "nT")
T_max <- max(panel_counts$nT, na.rm = TRUE)
balanced_countries <- panel_counts %>% filter(nT == T_max) %>% pull(country)

balanced_df <- raw_df %>% filter(country %in% balanced_countries)

# Compute TWFE for each numeric variable: x_it - x_i. - x_.t + x_..
# Translation note: year is Date; we treat it as time index for group means.
TWFE_balanced_df <- balanced_df %>%
  group_by(country) %>%
  mutate(across(all_of(numeric_cols), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_c_demean")) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(across(ends_with("_c_demean"), ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}_t_demean")) %>%
  ungroup()

# The above yields: (x_it - x_i.) - mean_t(x_it - x_i.) which equals x_it - x_i. - x_.t + x_..
# We expose final names like VAR_TWFE.
for (c in numeric_cols) {
  TWFE_balanced_df[[paste0(c, "_TWFE")]] <- TWFE_balanced_df[[paste0(c, "_c_demean_t_demean")]]
}
TWFE_balanced_df <- TWFE_balanced_df %>%
  select(year, country, ends_with("_TWFE"), everything()) %>%
  select(-ends_with("_c_demean"), -ends_with("_t_demean"))

# Variance ranking (TWFE)
twfe_vars <- paste0(numeric_cols, "_TWFE")
twfe_var_tbl <- TWFE_balanced_df %>%
  summarise(across(all_of(twfe_vars), ~ var(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "variance") %>%
  arrange(variance)

# ---------------------------
# 1.12 - 1.14: TWFE unbalanced + boxplots per country
# ---------------------------

# For unbalanced TWFE, the same formula applies using available observations.
# (x_it - x_i.) - (x_.t - x_..) equivalently x_it - x_i. - x_.t + x_..
# We compute it directly.

overall_means <- raw_df %>% summarise(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)))
# helper to pull overall mean for a column:
overall_mean <- function(col) as.numeric(overall_means[[col]])

TWFE_unbalanced_df <- raw_df %>%
  group_by(country) %>%
  mutate(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_i_mean")) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_t_mean")) %>%
  ungroup()

for (c in numeric_cols) {
  TWFE_unbalanced_df[[paste0(c, "_TWFE")]] <-
    TWFE_unbalanced_df[[c]] -
    TWFE_unbalanced_df[[paste0(c, "_i_mean")]] -
    TWFE_unbalanced_df[[paste0(c, "_t_mean")]] +
    overall_mean(c)
}

# Clean helper mean columns
TWFE_unbalanced_df <- TWFE_unbalanced_df %>%
  select(-ends_with("_i_mean"), -ends_with("_t_mean"))

# Example boxplot: distribution of one TWFE variable by country (customize variable)
if ("GDP_TWFE" %in% names(TWFE_unbalanced_df)) {
  p_box_gdp <- ggplot(TWFE_unbalanced_df, aes(x = country, y = GDP_TWFE)) +
    geom_boxplot(outlier.alpha = 0.3) +
    coord_flip() +
    theme_minimal() +
    labs(title = "GDP_TWFE by country", x = NULL, y = NULL)
}

# ---------------------------
# 1.15 / 1.16 / 1.17: Univariate + bivariate descriptive stats and correlations
# ---------------------------

# Univariate descriptive stats (similar to pandas describe)
desc_tbl <- raw_df %>%
  summarise(across(all_of(numeric_cols), list(
    n = ~ sum(!is.na(.x)),
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    p25 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    p75 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}__{.fn}"))

# Long format for easier reading/export
desc_long <- desc_tbl %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "stat"), sep = "__", remove = FALSE) %>%
  select(variable, stat, value) %>%
  arrange(variable, stat)

# Correlation matrices
corr_raw <- corr_with_p(raw_df %>% select(all_of(numeric_cols)), method = "pearson")
corr_within <- corr_with_p(trans_df %>% select(ends_with("_within")) %>% select(where(is.numeric)), method = "pearson")
corr_between <- corr_with_p(trans_df %>% select(ends_with("_between")) %>% select(where(is.numeric)), method = "pearson")

# Difficulty note: Python plotted annotated heatmaps with seaborn.
# In R, you can use ggplot2 + geom_tile, or packages like corrplot/ggcorrplot.
# Here we export matrices; plotting is optional.

# Example heatmap (optional, simple)
plot_corr_heatmap <- function(r_mat, title = "Correlation heatmap") {
  df <- as.data.frame(as.table(r_mat))
  names(df) <- c("Var1", "Var2", "r")
  ggplot(df, aes(x = Var1, y = Var2, fill = r)) +
    geom_tile() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL)
}

p_corr_raw <- plot_corr_heatmap(corr_raw$r, "Raw correlations")

# ============================================================
# 2) DETRENDING + LAG EXTRACTION (as in later notebook cells)
# ============================================================

# Python later constructed "within_detrended" columns and then added lags up to max_lag.
# We replicate with generic utilities.

# ---------------------------
# Country-specific detrending (within variables)
# ---------------------------

# Translation note:
# - Python likely used statsmodels OLS per country for each within variable on time trend, then took residuals.
# - In R we do: residuals(lm(var_within ~ t_numeric)) by country.
#
# Choose within variables to detrend. Adjust list as needed.
within_vars_to_detrend <- paste0(target_variables, "_within")
within_vars_to_detrend <- within_vars_to_detrend[within_vars_to_detrend %in% names(trans_df)]

trans_detrended_df <- trans_df %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(t_numeric = as.numeric(year)) %>%
  group_modify(function(.x, .g) {
    for (v in within_vars_to_detrend) {
      fit <- lm(.x[[v]] ~ .x$t_numeric)
      .x[[paste0(v, "_detrended")]] <- resid(fit)
    }
    .x
  }) %>%
  ungroup() %>%
  select(-t_numeric)

# ---------------------------
# Lag extraction
# ---------------------------

# Python used max_lag (variable); if not defined, choose a default.
max_lag <- 4

detrended_cols <- paste0(within_vars_to_detrend, "_detrended")
detrended_cols <- detrended_cols[detrended_cols %in% names(trans_detrended_df)]

trans_detrended_lag_df <- trans_detrended_df %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(across(all_of(detrended_cols), list(
    lag1 = ~ lag(.x, 1),
    lag2 = ~ lag(.x, 2),
    lag3 = ~ lag(.x, 3),
    lag4 = ~ lag(.x, 4)
  ), .names = "{.col}_L{.fn}")) %>%
  ungroup()

# If max_lag != 4, generalize:
add_lags <- function(df, cols, max_lag = 4) {
  df %>%
    arrange(country, year) %>%
    group_by(country) %>%
    mutate(across(all_of(cols),
                  .fns = set_names(as.list(seq_len(max_lag)), paste0("L", seq_len(max_lag))) |>
                    lapply(function(L) function(x) lag(x, L)),
                  .names = "{.col}_{.fn}")) %>%
    ungroup()
}

trans_detrended_lag_df <- add_lags(trans_detrended_df, detrended_cols, max_lag = max_lag)

# Correlation matrix (within + lags), as in Python cell:
within_lag_corr_matrix <- cor(trans_detrended_lag_df %>% select(where(is.numeric)), use = "pairwise.complete.obs")

# ============================================================
# 3) EXPORTS (final transformed df)
# ============================================================

# Python final merge:
# final_trans_df = trans_df.merge(trans_diff_df, on=["year","country"], how="outer")...
# .merge(TWFE_unbalanced_df, ...)
# drop duplicated columns
# export csv and xlsx

final_trans_df <- trans_df %>%
  full_join(trans_diff_df, by = c("year", "country")) %>%
  full_join(TWFE_unbalanced_df %>% select(year, country, ends_with("_TWFE")), by = c("year", "country")) %>%
  arrange(country, year)

# Export
dir.create("Data", showWarnings = FALSE)
write.csv(final_trans_df, "Data/final_trans_df.csv", row.names = FALSE)
write_xlsx(final_trans_df, "Data/final_trans_df.xlsx")