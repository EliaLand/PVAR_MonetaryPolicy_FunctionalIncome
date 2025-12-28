# ============================================================
# PANEL DATA ECONOMETRICS - Code File Part 2 (R translation)
# Source: /mnt/data/PDE_codefile_part2.ipynb
# ============================================================

# ---------------------------
# 0) REQUIREMENTS SET-UP
# ---------------------------

# PY (original):
# # Requirements.txt file installation
# # !pip install -r requirements.txt

# R: install packages once (uncomment if needed)
# install.packages(c(
#   "readr","dplyr","tidyr","lubridate","stringr",
#   "plm","lmtest","sandwich","broom",
#   "AER","tseries","ggplot2"
# ))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(plm)
  library(lmtest)
  library(sandwich)
  library(broom)
  library(AER)
  library(tseries)
  library(ggplot2)
})

options(warn = -1) # PY: warnings.filterwarnings("ignore") (approx)

# ---------------------------
# 1) HELPERS
# ---------------------------

# PY: significance stars helper
significance_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Robust SE helper (cluster by group, like common panel robust covariance)
robust_se <- function(model, cluster = c("group", "time")) {
  cluster <- match.arg(cluster)
  vc <- vcovHC(model, type = "HC1", cluster = cluster)
  sqrt(diag(vc))
}

# ---------------------------
# 2) DATA LOADING
# ---------------------------

# PY (original):
# final_trans_df = pd.read_csv("Data/final_trans_df.csv")
# final_trans_df["year"] = pd.to_datetime(final_trans_df["year"])
# final_trans_df = final_trans_df.set_index(["country", "year"])

final_trans_df <- read_csv("Data/final_trans_df.csv", show_col_types = FALSE) %>%
  mutate(year = as.Date(year))

# plm needs a pdata.frame (panel structure)
pdf <- pdata.frame(final_trans_df, index = c("country", "year"), drop.index = FALSE, row.names = TRUE)

# ---------------------------
# 3) PANEL MODELS: BETWEEN / WITHIN / RE (Mundlak) / TWFE / FD
# ---------------------------

# Python created Models 1..5 by swapping one labour variable.
# We implement a reusable function that runs the same estimator set.

run_panel_suite <- function(pdf, y_var, X_vars) {

  # BETWEEN estimator
  between_mod <- plm(
    as.formula(paste0(y_var, " ~ ", paste(X_vars, collapse = " + "))),
    data = pdf, model = "between", effect = "individual"
  )

  # OWFE-Within (one-way FE, entity effects)
  within_mod <- plm(
    as.formula(paste0(y_var, " ~ ", paste(X_vars, collapse = " + "))),
    data = pdf, model = "within", effect = "individual"
  )

  # RE-Mundlak:
  # Python used RandomEffects with Mundlak augmentation (adding country means of regressors).
  # In R we create mean(X|i) and include them as extra regressors in a random effects model.
  df_tmp <- as.data.frame(pdf)
  for (x in X_vars) {
    df_tmp[[paste0(x, "_mean_i")]] <- ave(df_tmp[[x]], df_tmp$country, FUN = function(z) mean(z, na.rm = TRUE))
  }
  mundlak_vars <- c(X_vars, paste0(X_vars, "_mean_i"))
  mundlak_pdf <- pdata.frame(df_tmp, index = c("country", "year"), drop.index = FALSE, row.names = TRUE)

  mundlak_mod <- plm(
    as.formula(paste0(y_var, " ~ ", paste(mundlak_vars, collapse = " + "))),
    data = mundlak_pdf, model = "random", effect = "individual"
  )

  # TWFE-Within: two-way fixed effects
  twfe_mod <- plm(
    as.formula(paste0(y_var, " ~ ", paste(X_vars, collapse = " + "))),
    data = pdf, model = "within", effect = "twoways"
  )

  # First Differences
  fd_mod <- plm(
    as.formula(paste0(y_var, " ~ ", paste(X_vars, collapse = " + "))),
    data = pdf, model = "fd"
  )

  list(
    Between = between_mod,
    OWFE_Within = within_mod,
    RE_Mundlak = mundlak_mod,
    TWFE_Within = twfe_mod,
    First_Differences = fd_mod
  )
}

print_panel_suite <- function(mods) {
  # Print coefficient tables with HC1 robust SE (cluster=group)
  for (nm in names(mods)) {
    cat("\n====================\n", nm, "\n====================\n", sep = "")
    m <- mods[[nm]]
    se <- robust_se(m, cluster = "group")
    print(coeftest(m, vcov. = vcovHC(m, type = "HC1", cluster = "group")))
  }
}

# ---------------------------
# Model 1
# ---------------------------
# PY (original intent):
# X_vars = ["UN", "P", "WR", "LS", "i"]; y_var="GDP"
mods_m1 <- run_panel_suite(pdf, y_var = "GDP", X_vars = c("UN", "P", "WR", "LS", "i"))
print_panel_suite(mods_m1)

# ---------------------------
# Model 2
# ---------------------------
mods_m2 <- run_panel_suite(pdf, y_var = "GDP", X_vars = c("UN", "P", "W", "LS", "i"))
print_panel_suite(mods_m2)

# ---------------------------
# Model 3
# ---------------------------
mods_m3 <- run_panel_suite(pdf, y_var = "GDP", X_vars = c("SHORTUN", "P", "WR", "LS", "i"))
print_panel_suite(mods_m3)

# ---------------------------
# Model 4
# ---------------------------
mods_m4 <- run_panel_suite(pdf, y_var = "GDP", X_vars = c("LONGUN", "P", "WR", "LS", "i"))
print_panel_suite(mods_m4)

# ---------------------------
# Model 5
# ---------------------------
mods_m5 <- run_panel_suite(pdf, y_var = "GDP", X_vars = c("LF", "P", "WR", "LS", "i"))
print_panel_suite(mods_m5)

# DIFFICULTY NOTE:
# The Python notebook used linearmodels.compare(results) to create a neat side-by-side table.
# In R, you can get a similar table with packages like modelsummary, stargazer, or texreg.
# This translation prints robust coeftest tables for each estimator to remain dependency-light.

# ---------------------------
# 4) AUTOCORRELATION CHECK IN FIRST DIFFERENCES (MODEL 1)
# ---------------------------

# PY idea:
# - compute first differences by country
# - estimate FD model, extract residuals
# - regress resid_t on resid_{t-1} (AR(1))

df_fd <- as.data.frame(pdf) %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(across(c("GDP","UN","P","WR","LS","i"), ~ .x - lag(.x), .names = "{.col}_1diff")) %>%
  ungroup() %>%
  filter(!is.na(GDP_1diff))

# FD regression (PanelOLS on first-differenced data in Python; here plain OLS on diff rows is fine)
fd_fit <- lm(GDP_1diff ~ UN_1diff + P_1diff + WR_1diff + LS_1diff + i_1diff, data = df_fd)
fd_resid <- resid(fd_fit)

# Lag residuals within country
df_res <- df_fd %>%
  mutate(resid = fd_resid) %>%
  group_by(country) %>%
  mutate(resid_lag = lag(resid, 1)) %>%
  ungroup() %>%
  filter(!is.na(resid_lag))

ar1_fit <- lm(resid ~ resid_lag, data = df_res)
cat("\n====================\nAR(1) regression of FD residuals\n====================\n")
print(summary(ar1_fit))

# ---------------------------
# 5) ANDERSON-HSIAO ARDL SET-UP (DIFFERENCES + LAGS + LEVEL IVs)
# ---------------------------

# PY:
# - use *_1diff columns already available in final_trans_df
# - create lag1 of 1diff vars
# - create lag2 of LEVEL vars as instruments

df2 <- as.data.frame(pdf) %>% arrange(country, year)

# Ensure the first-difference columns exist; if not, construct them.
need_1diff <- c("GDP","UN","P","WR","LS","i")
for (v in need_1diff) {
  nm <- paste0(v, "_1diff")
  if (!nm %in% names(df2)) {
    df2 <- df2 %>%
      group_by(country) %>%
      mutate("{nm}" := .data[[v]] - lag(.data[[v]], 1)) %>%
      ungroup()
  }
}

# Lagged first differences (t-1)
df2 <- df2 %>%
  group_by(country) %>%
  mutate(
    GDP_1diff_lag1 = lag(GDP_1diff, 1),
    UN_1diff_lag1  = lag(UN_1diff, 1),
    P_1diff_lag1   = lag(P_1diff, 1),
    WR_1diff_lag1  = lag(WR_1diff, 1),
    LS_1diff_lag1  = lag(LS_1diff, 1),
    i_1diff_lag1   = lag(i_1diff, 1)
  ) %>%
  ungroup()

# Lagged levels (t-2) as IVs
df2 <- df2 %>%
  group_by(country) %>%
  mutate(
    GDP_lag2 = lag(GDP, 2),
    P_lag2   = lag(P, 2),
    UN_lag2  = lag(UN, 2),
    WR_lag2  = lag(WR, 2),
    LS_lag2  = lag(LS, 2),
    i_lag2   = lag(i, 2)
  ) %>%
  ungroup()

diff_ah_df <- df2 %>%
  filter(complete.cases(GDP_1diff, GDP_1diff_lag1, GDP_lag2, P_1diff, P_1diff_lag1, P_lag2,
                        UN_1diff, WR_1diff, LS_1diff, i_1diff))

# ---------------------------
# 6) IV RELEVANCE CHECK (first-stage on GDP_1diff_lag1 ~ GDP_lag2 + P_lag2)
# ---------------------------

fs_y <- lm(GDP_1diff_lag1 ~ GDP_lag2 + P_lag2, data = diff_ah_df)
cat("\n====================\nFirst stage relevance check\n====================\n")
print(summary(fs_y))

# ---------------------------
# 7) ANDERSON-HSIAO ARDL MODEL (IV2SLS)
# ---------------------------

# PY (with truncated "..." inside regressor list):
# X included GDP_1diff_lag1, P_1diff_lag1, P_1diff, UN_1diff, WR_1diff, LS_1diff, ... (likely i_1diff)
# Z included exogenous controls + P terms + GDP_lag2 + P_lag2
#
# DIFFICULTY NOTE:
# The notebook cell contains a literal "..." so one regressor was omitted.
# Based on later cells (OLS/IV specs), we include i_1diff in the structural equation.

# Structural (second stage) equation: treat GDP_1diff_lag1 as endogenous
# Exogenous: P_1diff + P_1diff_lag1 + UN_1diff + WR_1diff + LS_1diff + i_1diff
# Instruments: GDP_lag2 + P_lag2 (excluded) plus included exogenous
iv_formula <- GDP_1diff ~ P_1diff + P_1diff_lag1 + UN_1diff + WR_1diff + LS_1diff + i_1diff +
  GDP_1diff_lag1 | P_1diff + P_1diff_lag1 + UN_1diff + WR_1diff + LS_1diff + i_1diff + GDP_lag2 + P_lag2

iv_fit <- ivreg(iv_formula, data = diff_ah_df)
cat("\n====================\nIV (AER::ivreg) - Anderson-Hsiao style\n====================\n")
print(summary(iv_fit, diagnostics = TRUE))

# ---------------------------
# 8) UNIVARIATE DESCRIPTIVE STATS (AH variables)
# ---------------------------

target_variables <- c(
  "GDP_1diff", "GDP_1diff_lag1", "GDP_lag2",
  "P_1diff", "P_1diff_lag1", "P_lag2",
  "UN_1diff", "WR_1diff", "LS_1diff", "i_1diff"
)

desc_df <- diff_ah_df %>%
  select(all_of(target_variables)) %>%
  summarise(across(everything(), list(
    n = ~ sum(!is.na(.x)),
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    min = ~ min(.x, na.rm = TRUE),
    p25 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    p75 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE),
    skewness = ~ {
      x <- .x[!is.na(.x)]
      if (length(x) < 3) NA_real_ else mean((x - mean(x))^3) / (sd(x)^3)
    },
    kurtosis = ~ {
      x <- .x[!is.na(.x)]
      if (length(x) < 4) NA_real_ else mean((x - mean(x))^4) / (sd(x)^4)
    }
  ), .names = "{.col}__{.fn}"))

diff_ah_desc_df <- desc_df %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("variable", "stat"), sep = "__") %>%
  arrange(variable, stat)

print(diff_ah_desc_df)

# ---------------------------
# 9) CORRELATION HEATMAP (with t-stats) - AH variables
# ---------------------------

# PY built an annotation matrix with r and t-stats and then used seaborn heatmap.
# R implementation using ggplot2 tiles + text labels.

df_corr <- diff_ah_df %>% select(all_of(target_variables)) %>% drop_na()
corr_mat <- cor(df_corr, use = "pairwise.complete.obs")
n <- nrow(df_corr)

t_from_r <- function(r, n) r * sqrt((n - 2) / pmax(1e-12, 1 - r^2))

corr_long <- as.data.frame(as.table(corr_mat)) %>%
  rename(Var1 = Var1, Var2 = Var2, r = Freq) %>%
  mutate(t = t_from_r(r, n),
         label = sprintf("%.2f\n(%.2f)", r, t))

p_corr <- ggplot(corr_long, aes(x = Var1, y = Var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = label), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Anderson-Hsiao ARDL Model Instrumental Variables - Correlation Heatmap\n(r-value with t-statistics in parentheses)",
    x = NULL, y = NULL
  )

print(p_corr)

# ---------------------------
# 10) UNIT-ROOT TESTING (Fisher-ADF panel test, manual aggregation)
# ---------------------------

# PY:
# - run ADF per country series (maxlag=1, autolag=None)
# - collect p-values and combine with Fisher: stat = -2 sum(log p_i) ~ Chi-square(2k)

fisher_adf <- function(df, var, country_col = "country", min_len = 4, k_lag = 1) {
  pvals <- c()
  for (cty in unique(df[[country_col]])) {
    s <- df %>% filter(.data[[country_col]] == cty) %>% pull(.data[[var]])
    s <- s[!is.na(s)]
    if (length(s) >= min_len) {
      # tseries::adf.test uses k = number of lags of the differenced series
      # DIFFICULTY NOTE: Python's adfuller(maxlag=1, autolag=None) is closest to adf.test(k=1),
      # but the exact regression specification differs slightly across implementations.
      p <- tryCatch(adf.test(s, k = k_lag)$p.value, error = function(e) NA_real_)
      if (!is.na(p)) pvals <- c(pvals, p)
    }
  }
  stat <- -2 * sum(log(pvals))
  df_chi <- 2 * length(pvals)
  pval <- 1 - pchisq(stat, df = df_chi)
  list(stat = stat, df = df_chi, p.value = pval, k = length(pvals))
}

gdp_fisher <- fisher_adf(diff_ah_df, "GDP_1diff", k_lag = 1)
cat(sprintf("\nFisher-ADF (panel) for GDP_1diff: stat=%.3f, df=%d, p-value=%.6f (countries used=%d)\n",
            gdp_fisher$stat, gdp_fisher$df, gdp_fisher$p.value, gdp_fisher$k))

p_fisher <- fisher_adf(diff_ah_df, "P_1diff", k_lag = 1)
cat(sprintf("Fisher-ADF (panel) for P_1diff: stat=%.3f, df=%d, p-value=%.6f (countries used=%d)\n",
            p_fisher$stat, p_fisher$df, p_fisher$p.value, p_fisher$k))

# ---------------------------
# 11) OLS with lagged 1diff of GDP and explanatory vars
# ---------------------------

ols_dyn_df <- diff_ah_df %>% drop_na(
  GDP_1diff, GDP_1diff_lag1, P_1diff_lag1, UN_1diff_lag1, WR_1diff_lag1, LS_1diff_lag1, i_1diff_lag1
)

ols_dyn <- lm(
  GDP_1diff ~ GDP_1diff_lag1 + P_1diff_lag1 + UN_1diff_lag1 + WR_1diff_lag1 + LS_1diff_lag1 + i_1diff_lag1,
  data = ols_dyn_df
)

cat("\n====================\nOLS dynamic (HC1 robust)\n====================\n")
print(coeftest(ols_dyn, vcov. = vcovHC(ols_dyn, type = "HC1")))

# ---------------------------
# 12) OLS benchmark + IV2SLS (levels instruments)  (Python cell contained "..." truncation)
# ---------------------------

# Python intended an OLS benchmark:
# GDP_1diff ~ GDP_1diff_lag1 + P_1diff + P_1diff_lag1 + UN_1diff + WR_1diff + LS_1diff + i_1diff

ols_bench <- lm(
  GDP_1diff ~ GDP_1diff_lag1 + P_1diff + P_1diff_lag1 + UN_1diff + WR_1diff + LS_1diff + i_1diff,
  data = diff_ah_df
)

cat("\n====================\nOLS benchmark (HC1 robust)\n====================\n")
print(coeftest(ols_bench, vcov. = vcovHC(ols_bench, type = "HC1")))

# IV already estimated above as iv_fit (AER::ivreg). For convenience:
iv_model <- iv_fit

# ---------------------------
# 13) PARAMETER DEVIATION: OLS vs IV
# ---------------------------

ols_params <- coef(ols_bench)
iv_params  <- coef(iv_model)

common <- intersect(names(ols_params), names(iv_params))
comp <- tibble(
  term = common,
  OLS = unname(ols_params[common]),
  IV  = unname(iv_params[common])
) %>%
  mutate(
    abs_diff = abs(IV - OLS),
    pct_diff = 100 * (IV - OLS) / OLS
  )

print(comp)

# ---------------------------
# 14) IV QUALITY TESTS
# ---------------------------

# PY:
# print(fs_y.f_test("GDP_lag2 = 0, P_lag2 = 0"))
# print(iv.first_stage)
# print(iv.wu_hausman())

# R equivalents:
# - joint relevance: linearHypothesis on first stage
# - diagnostics: summary(ivreg, diagnostics=TRUE) includes weak instruments and Wu-Hausman.

cat("\n====================\nJoint F-test on excluded instruments (first stage)\n====================\n")
print(car::linearHypothesis(fs_y, c("GDP_lag2 = 0", "P_lag2 = 0")))

cat("\n====================\nIV diagnostics (weak instruments / Wu-Hausman)\n====================\n")
print(summary(iv_model, diagnostics = TRUE))

# ---------------------------
# 15) IMPULSE RESPONSE FUNCTION (IRF) for IV and OLS benchmark
# ---------------------------

# PY cells had "..." where IRF computation was; we reconstruct from the ARDL form:
# GDP_1diff_t = beta_y * GDP_1diff_{t-1} + beta1 * P_1diff_t + beta2 * P_1diff_{t-1} + ...
#
# For a 1-unit increase in P_1diff at horizon 1 (t=1) with P_1diff_{0}=0:
# h=1: beta1
# h=2: beta_y*beta1 + beta2
# h=3: beta_y^2*beta1 + beta_y*beta2
# h=4: beta_y^3*beta1 + beta_y^2*beta2

compute_irf <- function(beta_y, beta1, beta2, H = 4) {
  h <- 1:H
  irf <- numeric(H)
  for (k in 1:H) {
    if (k == 1) irf[k] <- beta1
    if (k >= 2) irf[k] <- (beta_y^(k-1)) * beta1 + (beta_y^(k-2)) * beta2
  }
  tibble(h = h, irf = irf)
}

# IV coefficients
b_iv <- coef(iv_model)
beta_y_iv <- unname(b_iv["GDP_1diff_lag1"])
beta1_iv  <- unname(b_iv["P_1diff"])
beta2_iv  <- unname(b_iv["P_1diff_lag1"])
irf_iv <- compute_irf(beta_y_iv, beta1_iv, beta2_iv, H = 4)

p_irf_iv <- ggplot(irf_iv, aes(x = h, y = irf)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = 1:4) +
  theme_minimal() +
  labs(
    title = "IV Model - IRF\n1-unit increase in P on GDP_1diff (horizon range 1-4)",
    x = "Horizon (t)",
    y = "Impulse response of GDP_1diff to 1-unit increase in P"
  )
print(p_irf_iv)

# OLS benchmark coefficients
b_ols <- coef(ols_bench)
beta_y_ols <- unname(b_ols["GDP_1diff_lag1"])
beta1_ols  <- unname(b_ols["P_1diff"])
beta2_ols  <- unname(b_ols["P_1diff_lag1"])
irf_ols <- compute_irf(beta_y_ols, beta1_ols, beta2_ols, H = 4)

p_irf_ols <- ggplot(irf_ols, aes(x = h, y = irf)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(breaks = 1:4) +
  theme_minimal() +
  labs(
    title = "OLS Benchmark - IRF\n1-unit increase in P on GDP_1diff (horizon range 1-4)",
    x = "Horizon (t)",
    y = "Impulse response of GDP_1diff to 1-unit increase in P"
  )
print(p_irf_ols)

# ---------------------------
# 16) LONG-RUN COEFFICIENTS (IV and OLS)
# ---------------------------

# PY:
# beta_LT = (beta1 + beta2) / (1 - beta_y)
beta_LT_iv  <- (beta1_iv + beta2_iv) / (1 - beta_y_iv)
beta_LT_ols <- (beta1_ols + beta2_ols) / (1 - beta_y_ols)

cat(sprintf("\nLong-run coefficient (IV):  %.4f\n", beta_LT_iv))
cat(sprintf("Long-run coefficient (OLS): %.4f\n", beta_LT_ols))
