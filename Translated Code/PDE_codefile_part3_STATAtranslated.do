
********************************************************************************
* PDE_codefile_part3_translated.do
********************************************************************************

version 16.0
clear all
set more off
set linesize 255

capture mkdir "Graphs"
capture mkdir "Outputs"

********************************************************************************
* 0) Input
********************************************************************************
global IN_CSV "Data/final_trans_df.csv"
capture confirm file "$IN_CSV"
if _rc {
    di as error "Input file not found: $IN_CSV"
    di as error "Please run Part 1 translation first, or edit global IN_CSV."
    exit 601
}

import delimited using "$IN_CSV", clear varnames(1)

* Create country_id and year_int robustly
capture confirm string variable country
if _rc {
    tostring country, replace
}
encode country, gen(country_id)

capture confirm variable year_int
if _rc {
    capture confirm numeric variable year
    if !_rc {
        gen int year_int = year
    }
    else {
        capture confirm string variable year
        if _rc {
            di as error "Could not locate year/year_int in the CSV."
            exit 498
        }
        gen double __tmpdate = date(year,"YMD")
        replace __tmpdate = date(year,"DMY") if missing(__tmpdate)
        replace __tmpdate = date(year,"MDY") if missing(__tmpdate)
        gen int year_int = year(__tmpdate)
        drop __tmpdate
    }
}
label var year_int "Year (numeric)"
xtset country_id year_int

********************************************************************************
* 1) Monetary policy shock (policy-rule residual)  (Python cell 11)
********************************************************************************
* Python:
*   rule = PanelOLS(i ~ const + P + GDP + FE_i + FE_t, clustered by entity)
*   mp_shock = residuals
capture confirm variable i
capture confirm variable P
capture confirm variable GDP
if _rc {
    di as error "Missing required variables i, P, GDP. Please verify input dataset."
    exit 498
}

quietly regress i P GDP i.country_id i.year_int, vce(cluster country_id)
predict double mp_shock, resid
label var mp_shock "Monetary policy shock: residual from i on P,GDP + FE_i + FE_t"

********************************************************************************
* 2) Local Projections (LP) for WR, horizons 0..10 (Python cell 12)
********************************************************************************
* Specification:
*   WR_{i,t+h} = a_h + b_h*mp_shock_{i,t} + controls(GDP,P) + FE_i + FE_t + u
* Cluster SEs by country.
capture confirm variable WR
if _rc {
    di as error "Variable WR not found."
    exit 498
}

tempfile lp_wr
tempname post_wr
postfile `post_wr' int Horizon double Beta_shock_SE Standard_Error p_value N_obs using `lp_wr', replace

forvalues h=0/10 {
    * Lead outcome within country
    by country_id (year_int): gen double WR_lead`h' = F`h'.WR

    * Restrict to valid obs (drop missing in y or regressors)
    quietly regress WR_lead`h' mp_shock GDP P i.country_id i.year_int ///
        if !missing(WR_lead`h', mp_shock, GDP, P), vce(cluster country_id)

    * Store coefficient on mp_shock
    local b  = _b[mp_shock]
    local se = _se[mp_shock]
    local p  = 2*ttail(e(df_r), abs(`b'/`se'))
    local N  = e(N)

    post `post_wr' (`h') (`b') (`se') (`p') (`N')

    drop WR_lead`h'
}
postclose `post_wr'

use `lp_wr', clear
rename Beta_shock_SE Beta_shock_WR
label var Beta_shock_WR "LP coefficient on mp_shock (WR)"
label var Standard_Error "Std. Error (clustered by country)"
label var p_value "p-value"
label var N_obs "N observations"
sort Horizon

list, abbrev(24)
export delimited using "Outputs/lp_wr.csv", replace

* Plot IRF curve (WR)
twoway (line Beta_shock_WR Horizon, msymbol(o)), ///
    yline(0, lpattern(dash)) ///
    title("LP IRF: WR response to mp_shock") ///
    xtitle("Horizon h") ytitle("Coefficient on mp_shock")
graph export "Graphs/lp_irf_wr.png", replace

********************************************************************************
* 3) Local Projections (LP) for LS, horizons 0..10 (Python cell 13)
********************************************************************************
capture confirm variable LS
if _rc {
    di as error "Variable LS not found."
    exit 498
}

tempfile lp_ls
tempname post_ls
postfile `post_ls' int Horizon double Beta_shock_LS Standard_Error p_value N_obs using `lp_ls', replace

forvalues h=0/10 {
    by country_id (year_int): gen double LS_lead`h' = F`h'.LS

    quietly regress LS_lead`h' mp_shock GDP P i.country_id i.year_int ///
        if !missing(LS_lead`h', mp_shock, GDP, P), vce(cluster country_id)

    local b  = _b[mp_shock]
    local se = _se[mp_shock]
    local p  = 2*ttail(e(df_r), abs(`b'/`se'))
    local N  = e(N)

    post `post_ls' (`h') (`b') (`se') (`p') (`N')

    drop LS_lead`h'
}
postclose `post_ls'

use `lp_ls', clear
label var Beta_shock_LS "LP coefficient on mp_shock (LS)"
label var Standard_Error "Std. Error (clustered by country)"
label var p_value "p-value"
label var N_obs "N observations"
sort Horizon

list, abbrev(24)
export delimited using "Outputs/lp_ls.csv", replace

* Plot IRF curve (LS)
twoway (line Beta_shock_LS Horizon, msymbol(o)), ///
    yline(0, lpattern(dash)) ///
    title("LP IRF: LS response to mp_shock") ///
    xtitle("Horizon h") ytitle("Coefficient on mp_shock")
graph export "Graphs/lp_irf_ls.png", replace

********************************************************************************
* 4) Diagnostic placebo/pre-trend check: future shock should not predict today
*    (Python cells 15-16)
********************************************************************************
by country_id (year_int): gen double shock_lead1 = F1.mp_shock

* Placebo test - WR_t on mp_shock_{t+1}
quietly regress WR shock_lead1 GDP P i.country_id i.year_int ///
    if !missing(WR, shock_lead1, GDP, P), vce(cluster country_id)
estimates store placebo_wr

* Placebo test - LS_t on mp_shock_{t+1}
quietly regress LS shock_lead1 GDP P i.country_id i.year_int ///
    if !missing(LS, shock_lead1, GDP, P), vce(cluster country_id)
estimates store placebo_ls

estimates table placebo_wr placebo_ls, ///
    b(%10.4f) se(%10.4f) stats(N r2, fmt(%9.0g %9.3f) labels("N" "R2")) ///
    title("Placebo test: outcome_t on future shock_{t+1}") star

drop shock_lead1

********************************************************************************
* 5) Regime definition: High vs Low inflation (Python cells 14/17)
********************************************************************************
* Python proxy: inflation = diff(P) within country; define high_infl by median split.
by country_id (year_int): gen double infl = D.P
summarize infl, detail
scalar infl_med = r(p50)
gen byte high_infl = (infl > infl_med) if !missing(infl)
label var infl "Inflation proxy: first difference of P"
label var high_infl "High inflation regime (infl > sample median)"

* Regime-specific shocks (Python cell 19)
gen double shock_high = mp_shock * high_infl
gen double shock_low  = mp_shock * (1 - high_infl)
label var shock_high "Regime shock: mp_shock * high_infl"
label var shock_low  "Regime shock: mp_shock * (1-high_infl)"

********************************************************************************
* 6) Recompute mp_shock on this df (Python cell 18) — already computed above
********************************************************************************
* The notebook recomputed mp_shock after copying df; in Stata we already have mp_shock.

********************************************************************************
* 7) State-dependent LP for WR (Python cell 20)
********************************************************************************
tempfile sd_wr
tempname post_sdwr
postfile `post_sdwr' int Horizon ///
    double Beta_high SE_high p_high ///
    double Beta_low  SE_low  p_low ///
    double N_obs using `sd_wr', replace

forvalues h=0/10 {
    by country_id (year_int): gen double WR_lead`h' = F`h'.WR

    quietly regress WR_lead`h' shock_high shock_low GDP P i.country_id i.year_int ///
        if !missing(WR_lead`h', shock_high, shock_low, GDP, P), vce(cluster country_id)

    local bh  = _b[shock_high]
    local seh = _se[shock_high]
    local ph  = 2*ttail(e(df_r), abs(`bh'/`seh'))

    local bl  = _b[shock_low]
    local sel = _se[shock_low]
    local pl  = 2*ttail(e(df_r), abs(`bl'/`sel'))

    local N   = e(N)

    post `post_sdwr' (`h') (`bh') (`seh') (`ph') (`bl') (`sel') (`pl') (`N')

    drop WR_lead`h'
}
postclose `post_sdwr'

use `sd_wr', clear
label var Beta_high "Coef on shock_high (WR)"
label var Beta_low  "Coef on shock_low (WR)"
sort Horizon
export delimited using "Outputs/sd_wr.csv", replace

twoway ///
 (line Beta_high Horizon, msymbol(o)) ///
 (line Beta_low  Horizon, msymbol(o)), ///
 yline(0, lpattern(dash)) ///
 title("State-dependent LP IRF (WR): High vs Low inflation") ///
 xtitle("Horizon h") ytitle("LP coefficient") ///
 legend(order(1 "High inflation" 2 "Low inflation"))
graph export "Graphs/sd_lp_irf_wr.png", replace

********************************************************************************
* 8) State-dependent LP for LS (Python cell 21)
********************************************************************************
tempfile sd_ls
tempname post_sdls
postfile `post_sdls' int Horizon ///
    double Beta_high SE_high p_high ///
    double Beta_low  SE_low  p_low ///
    double N_obs using `sd_ls', replace

forvalues h=0/10 {
    by country_id (year_int): gen double LS_lead`h' = F`h'.LS

    quietly regress LS_lead`h' shock_high shock_low GDP P i.country_id i.year_int ///
        if !missing(LS_lead`h', shock_high, shock_low, GDP, P), vce(cluster country_id)

    local bh  = _b[shock_high]
    local seh = _se[shock_high]
    local ph  = 2*ttail(e(df_r), abs(`bh'/`seh'))

    local bl  = _b[shock_low]
    local sel = _se[shock_low]
    local pl  = 2*ttail(e(df_r), abs(`bl'/`sel'))

    local N   = e(N)

    post `post_sdls' (`h') (`bh') (`seh') (`ph') (`bl') (`sel') (`pl') (`N')

    drop LS_lead`h'
}
postclose `post_sdls'

use `sd_ls', clear
label var Beta_high "Coef on shock_high (LS)"
label var Beta_low  "Coef on shock_low (LS)"
sort Horizon
export delimited using "Outputs/sd_ls.csv", replace

twoway ///
 (line Beta_high Horizon, msymbol(o)) ///
 (line Beta_low  Horizon, msymbol(o)), ///
 yline(0, lpattern(dash)) ///
 title("State-dependent LP IRF (LS): High vs Low inflation") ///
 xtitle("Horizon h") ytitle("LP coefficient") ///
 legend(order(1 "High inflation" 2 "Low inflation"))
graph export "Graphs/sd_lp_irf_ls.png", replace

********************************************************************************
* 9) Wald tests: H0 beta_high(h) = beta_low(h) for each horizon
*    (Python cells 22-23)
********************************************************************************

* WR Wald tests
tempfile wt_wr
tempname post_wtwr
postfile `post_wtwr' int Horizon double WaldStat p_value using `wt_wr', replace

forvalues h=0/10 {
    by country_id (year_int): gen double WR_lead`h' = F`h'.WR

    quietly regress WR_lead`h' shock_high shock_low GDP P i.country_id i.year_int ///
        if !missing(WR_lead`h', shock_high, shock_low, GDP, P), vce(cluster country_id)

    * H0: shock_high - shock_low = 0  <=>  shock_high = shock_low
    test shock_high = shock_low
    local F = r(F)
    local p = r(p)

    post `post_wtwr' (`h') (`F') (`p')

    drop WR_lead`h'
}
postclose `post_wtwr'
use `wt_wr', clear
sort Horizon
export delimited using "Outputs/wald_test_wr.csv", replace

* LS Wald tests
tempfile wt_ls
tempname post_wtls
postfile `post_wtls' int Horizon double WaldStat p_value using `wt_ls', replace

forvalues h=0/10 {
    by country_id (year_int): gen double LS_lead`h' = F`h'.LS

    quietly regress LS_lead`h' shock_high shock_low GDP P i.country_id i.year_int ///
        if !missing(LS_lead`h', shock_high, shock_low, GDP, P), vce(cluster country_id)

    test shock_high = shock_low
    local F = r(F)
    local p = r(p)

    post `post_wtls' (`h') (`F') (`p')

    drop LS_lead`h'
}
postclose `post_wtls'
use `wt_ls', clear
sort Horizon
export delimited using "Outputs/wald_test_ls.csv", replace

