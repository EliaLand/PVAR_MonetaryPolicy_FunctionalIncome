# ============================================================
# PANEL DATA ECONOMETRICS - Code File Part 3 (R translation)
# Source: /mnt/data/PDE_codefile_part3.ipynb
# ============================================================

# ---------------------------
# Requirements / packages
# ---------------------------

# PY (original):
# # Requirements.txt file installation
# # !pip install -r requirements.txt

# R: install packages once (uncomment if needed)
# install.packages(c("readr","dplyr","tidyr","lubridate","stringr","plm","lmtest","sandwich","car","ggplot2"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(plm)
  library(lmtest)
  library(sandwich)
  library(car)
  library(ggplot2)
})

# ---------------------------
# Statistical Significance labelling  (PY: significance_stars)
# ---------------------------
significance_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# We supress potential warnings with this command (PY: warnings.filterwarnings("ignore"))
options(warn = -1)

# ---------------------------
# Data loading from CodeFile Part 1 for raw data and transformed variables
# ---------------------------

# PY:
# final_trans_df = pd.read_csv("Data/final_trans_df.csv")
# final_trans_df["year"] = pd.to_datetime(final_trans_df["year"])
# final_trans_df = final_trans_df.set_index(["country", "year"])

final_trans_df <- read_csv("Data/final_trans_df.csv", show_col_types = FALSE) %>%
  mutate(year = as.Date(year))

# (!!!) Ensure a panel index (i="country", t="year")
pdf <- pdata.frame(final_trans_df, index = c("country", "year"), drop.index = FALSE, row.names = TRUE)

# ============================================================
# Panel Local Projections (LP) IRFs
# ============================================================
# We build a "monetary policy shock" as a policy-rule residual
# (!!!) We want to simulate an “unexpected tightening” measure. The idea is:
#       policy rate i_it = systematic part (function of macro controls + FE) + residual
#       residual = “innovation/shock”.
# (!!!) We also include: country fixed effects (entity_effects) and year fixed effects (time_effects)
#       to absorb global shocks common to all countries.
# We cluster SEs by country because macro panels have within-country serial dependence.

# ---------------------------
# Helper: robust clustered SE (cluster by country)
# ---------------------------
vcov_cluster_country <- function(plm_model) {
  vcovHC(plm_model, type = "HC1", cluster = "group")
}

# ---------------------------
# Helper: run a two-way FE LP for a given outcome and horizon
# ---------------------------
# PY loop had an ellipsis "..." inside; reconstructed here:
# - create y_{t+h} by shifting within country
# - regress y_{t+h} on mp_shock_t + controls (GDP, P) with entity and time FE
# - cluster SE by country

run_lp_twfe <- function(df, outcome, shock, controls = c("GDP","P"), horizons = 0:10) {
  rows <- list()

  for (h in horizons) {
    df_h <- df %>%
      arrange(country, year) %>%
      group_by(country) %>%
      mutate(y_lead = dplyr::lead(.data[[outcome]], n = h)) %>%
      ungroup()

    df_h <- df_h %>%
      filter(!is.na(y_lead), !is.na(.data[[shock]])) %>%
      filter(if_all(all_of(controls), ~ !is.na(.x)))

    if (nrow(df_h) == 0) {
      rows[[length(rows)+1]] <- data.frame(
        Horizon = h, beta = NA_real_, se = NA_real_, p = NA_real_, nobs = 0
      )
      next
    }

    pdf_h <- pdata.frame(df_h, index = c("country","year"), drop.index = FALSE)

    fml <- as.formula(
      paste0("y_lead ~ ", shock, " + ", paste(controls, collapse = " + "))
    )

    mod <- plm(fml, data = pdf_h, model = "within", effect = "twoways")

    ct <- coeftest(mod, vcov. = vcov_cluster_country(mod))
    beta <- unname(ct[shock, "Estimate"])
    se   <- unname(ct[shock, "Std. Error"])
    pval <- unname(ct[shock, "Pr(>|t|)"])

    rows[[length(rows)+1]] <- data.frame(
      Horizon = h, beta = beta, se = se, p = pval, nobs = nobs(mod)
    )
  }

  bind_rows(rows)
}

# ---------------------------
# 1) Monetary policy shock from policy rule (two-way FE)
# ---------------------------

df <- final_trans_df %>% as.data.frame()

# Controls
rhs <- c("P", "GDP")

# ///////////
# Policy Rule
# ///////////
# (!!!) i_it = a + b*P_it + c*GDP_it + FE_i + FE_t + u_it
# R/plm equivalent: two-way FE ("twoways") within estimator
rule_mod <- plm(
  i ~ P + GDP,
  data = pdata.frame(df, index = c("country","year"), drop.index = FALSE),
  model = "within",
  effect = "twoways"
)

cat("\n====================\nPolicy rule (two-way FE) with clustered SEs\n====================\n")
print(coeftest(rule_mod, vcov. = vcov_cluster_country(rule_mod)))

# (!!!) The residuals are the "monetary policy shocks" (unexpected component)
df$mp_shock <- residuals(rule_mod)

mon_pol_df <- df

# ============================================================
# Local Projections for WR (horizon from 0 to 10)
# ============================================================
# (!!!) We deploy local projections, instead of imposing VAR dynamics, as LP estimates each horizon directly.
# (!!!) WR_i,t+h = a_h + b_h*shock_i,t + controls + FE_i + FE_t + u_i,t+h
# (!!!) GDP and P as controls to avoid attributing systematic macro movements

lp_wr <- run_lp_twfe(mon_pol_df, outcome = "WR", shock = "mp_shock", controls = c("GDP","P"), horizons = 0:10) %>%
  transmute(
    Horizon = Horizon,
    `Beta_shock (WR)` = beta,
    `Standard Error` = se,
    `p-value` = p,
    `N-obs` = nobs
  )

lp_wr

# ============================================================
# Local Projections for LS (horizon from 0 to 10)
# ============================================================
# (!!!) LS_i,t+h = a_h + b_h*shock_i,t + controls + FE_i + FE_t + u_i,t+h

lp_ls <- run_lp_twfe(mon_pol_df, outcome = "LS", shock = "mp_shock", controls = c("GDP","P"), horizons = 0:10) %>%
  transmute(
    Horizon = Horizon,
    `Beta_shock (LS)` = beta,
    `Standard Error` = se,
    `p-value` = p,
    `N-obs` = nobs
  )

lp_ls

# ============================================================
# IRF-curves Plotting  (PY: matplotlib; R: ggplot2)
# ============================================================
# (!!!) Code from Part 2

# Real Compensation (WR)
p_wr <- ggplot(lp_wr, aes(x = Horizon, y = `Beta_shock (WR)`)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Panel Local Projection IRF (WR) to 1-unit mp_shock",
    x = "Horizon h",
    y = "LP coefficient on mp_shock"
  )
print(p_wr)

# Labour Share (LS)
p_ls <- ggplot(lp_ls, aes(x = Horizon, y = `Beta_shock (LS)`)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Panel Local Projection IRF (LS) to 1-unit mp_shock",
    x = "Horizon h",
    y = "LP coefficient on mp_shock"
  )
print(p_ls)

# ============================================================
# WR - Diagnostic Test - Placebo/pre-trend check (shock lead)
# ============================================================
# (!!!) If the shock is truly unanticipated, a future shock should not predict today’s outcomes.
# (!!!) So when we regress WR_t and LS_t on mp_shock_{t+1}, we want the coefficient to be close to zero.

df_placebo <- mon_pol_df %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(shock_lead1 = dplyr::lead(mp_shock, 1)) %>%
  ungroup()

# Place test - WR
placebo_wr <- plm(
  WR ~ shock_lead1 + GDP + P,
  data = pdata.frame(df_placebo, index = c("country","year"), drop.index = FALSE),
  model = "within",
  effect = "twoways"
)

cat("\n====================\nPlacebo / pre-trend check (WR)\n====================\n")
print(coeftest(placebo_wr, vcov. = vcov_cluster_country(placebo_wr)))

# ============================================================
# LS - Diagnostic Test - Placebo/pre-trend check (shock lead)
# ============================================================
placebo_ls <- plm(
  LS ~ shock_lead1 + GDP + P,
  data = pdata.frame(df_placebo, index = c("country","year"), drop.index = FALSE),
  model = "within",
  effect = "twoways"
)

cat("\n====================\nPlacebo / pre-trend check (LS)\n====================\n")
print(coeftest(placebo_ls, vcov. = vcov_cluster_country(placebo_ls)))

# ============================================================
# State-dependent LP (High vs Low inflation)
# ============================================================
# Regime Definition = high inflation vs low inflation
# (!!!) State dependence: responses may differ depending on macro conditions.
# We proxy inflation with P_1diff (Price Deflator)
# (in PY: infl = groupby(country)["P"].diff())

df_sd <- final_trans_df %>%
  as.data.frame() %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(infl = P - dplyr::lag(P, 1)) %>%
  ungroup()

# (!!!) Split at median
med <- median(df_sd$infl, na.rm = TRUE)
df_sd$high_infl <- ifelse(df_sd$infl > med, 1L, 0L)

# ============================================================
# mp_shock Computation (recompute for new df)
# ============================================================
rule_mod_sd <- plm(
  i ~ P + GDP,
  data = pdata.frame(df_sd, index = c("country","year"), drop.index = FALSE),
  model = "within",
  effect = "twoways"
)

df_sd$mp_shock <- residuals(rule_mod_sd)

# ============================================================
# Regime-specific Shock Variables
# ============================================================
# (!!!) From the coefficients on shock_high and shock_low, we get the regime-specific effects.
df_sd$shock_high <- df_sd$mp_shock * df_sd$high_infl
df_sd$shock_low  <- df_sd$mp_shock * (1 - df_sd$high_infl)

# ============================================================
# Helper: state-dependent LP (two shocks in one regression)
# ============================================================

run_sd_lp_twfe <- function(df, outcome, controls = c("GDP","P"), horizons = 0:10) {
  rows <- list()

  for (h in horizons) {
    df_h <- df %>%
      arrange(country, year) %>%
      group_by(country) %>%
      mutate(y_lead = dplyr::lead(.data[[outcome]], n = h)) %>%
      ungroup()

    df_h <- df_h %>%
      filter(!is.na(y_lead), !is.na(shock_high), !is.na(shock_low)) %>%
      filter(if_all(all_of(controls), ~ !is.na(.x)))

    if (nrow(df_h) == 0) {
      rows[[length(rows)+1]] <- data.frame(
        Horizon = h,
        Beta_high = NA_real_, SE_high = NA_real_, p_high = NA_real_,
        Beta_low  = NA_real_, SE_low  = NA_real_, p_low  = NA_real_,
        N_obs = 0
      )
      next
    }

    pdf_h <- pdata.frame(df_h, index = c("country","year"), drop.index = FALSE)

    fml <- as.formula(
      paste0("y_lead ~ shock_high + shock_low + ", paste(controls, collapse = " + "))
    )

    mod <- plm(fml, data = pdf_h, model = "within", effect = "twoways")
    ct <- coeftest(mod, vcov. = vcov_cluster_country(mod))

    rows[[length(rows)+1]] <- data.frame(
      Horizon = h,
      Beta_high = unname(ct["shock_high", "Estimate"]),
      SE_high   = unname(ct["shock_high", "Std. Error"]),
      p_high    = unname(ct["shock_high", "Pr(>|t|)"]),
      Beta_low  = unname(ct["shock_low", "Estimate"]),
      SE_low    = unname(ct["shock_low", "Std. Error"]),
      p_low     = unname(ct["shock_low", "Pr(>|t|)"]),
      N_obs     = nobs(mod)
    )
  }

  bind_rows(rows)
}

# ============================================================
# State-dependent LP for WR (horizon from 0 to 10)
# ============================================================
sd_wr <- run_sd_lp_twfe(df_sd, outcome = "WR", controls = c("GDP","P"), horizons = 0:10)
sd_wr

# ============================================================
# State-dependent LP for LS (horizon from 0 to 10)
# ============================================================
sd_ls <- run_sd_lp_twfe(df_sd, outcome = "LS", controls = c("GDP","P"), horizons = 0:10)
sd_ls

# ============================================================
# WR - Wald Test (statistical test)
# ============================================================
# (!!!) Test H0: beta_high(h) = beta_low(h)

run_wald_by_horizon <- function(df, outcome, controls = c("GDP","P"), horizons = 0:10) {
  out <- list()

  for (h in horizons) {
    df_h <- df %>%
      arrange(country, year) %>%
      group_by(country) %>%
      mutate(y_lead = dplyr::lead(.data[[outcome]], n = h)) %>%
      ungroup() %>%
      filter(!is.na(y_lead), !is.na(shock_high), !is.na(shock_low)) %>%
      filter(if_all(all_of(controls), ~ !is.na(.x)))

    if (nrow(df_h) == 0) {
      out[[length(out)+1]] <- data.frame(Horizon = h, `Wald-Stat` = NA_real_, `p-value` = NA_real_, note = "no data")
      next
    }

    mod <- plm(
      as.formula(paste0("y_lead ~ shock_high + shock_low + ", paste(controls, collapse = " + "))),
      data = pdata.frame(df_h, index = c("country","year"), drop.index = FALSE),
      model = "within",
      effect = "twoways"
    )

    # H0: shock_high - shock_low = 0
    # Use robust clustered vcov like in the main LP regressions
    lh <- car::linearHypothesis(mod, "shock_high - shock_low = 0", vcov. = vcov_cluster_country(mod), test = "Chisq")

    stat <- as.numeric(lh[2, "Chisq"])
    pval <- as.numeric(lh[2, "Pr(>Chisq)"])

    out[[length(out)+1]] <- data.frame(
      Horizon = h,
      `Wald-Stat` = stat,
      `p-value` = pval,
      note = "dropped/absorbed"
    )
  }

  bind_rows(out) %>% arrange(Horizon) %>% mutate(Horizon = as.integer(Horizon))
}

wt_wr <- run_wald_by_horizon(df_sd, outcome = "WR", controls = c("GDP","P"), horizons = 0:10)
wt_wr

# ============================================================
# LS - Wald Test (statistical test)
# ============================================================
wt_ls <- run_wald_by_horizon(df_sd, outcome = "LS", controls = c("GDP","P"), horizons = 0:10)
wt_ls

# ============================================================
# IRF-curves Plotting (state-dependent)
# ============================================================

# Real Compensation (WR)
p_sd_wr <- ggplot(sd_wr, aes(x = Horizon)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = Beta_high), linewidth = 1) +
  geom_point(aes(y = Beta_high), size = 2) +
  geom_line(aes(y = Beta_low), linewidth = 1) +
  geom_point(aes(y = Beta_low), size = 2) +
  theme_minimal() +
  labs(
    title = "State-dependent LP IRF (WR) at Different Regimes\nHigh against Low Inflation",
    x = "Horizon h",
    y = "LP coefficient on mp_shock"
  )
print(p_sd_wr)

# Labour Share (LS)
p_sd_ls <- ggplot(sd_ls, aes(x = Horizon)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = Beta_high), linewidth = 1) +
  geom_point(aes(y = Beta_high), size = 2) +
  geom_line(aes(y = Beta_low), linewidth = 1) +
  geom_point(aes(y = Beta_low), size = 2) +
  theme_minimal() +
  labs(
    title = "State-dependent LP IRF (LS) at Different Regimes\nHigh against Low Inflation",
    x = "Horizon h",
    y = "LP coefficient on mp_shock"
  )
print(p_sd_ls)
