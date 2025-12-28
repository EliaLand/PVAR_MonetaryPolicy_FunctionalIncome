
********************************************************************************
* PDE_codefile_part2_translated.do
********************************************************************************

version 16.0
clear all
set more off
set linesize 255

********************************************************************************
* 0) Paths / inputs
********************************************************************************
global IN_CSV   "Data/final_trans_df.csv"
capture confirm file "$IN_CSV"
if _rc {
    di as error "Input file not found: $IN_CSV"
    di as error "Please run Part 1 translation first, or edit global IN_CSV."
    exit 601
}

capture mkdir "Graphs"
capture mkdir "Outputs"

********************************************************************************
* 1) "Requirements.txt installation" (Python cell 4)
********************************************************************************
* In Stata there is no pip-install step.
* If you need extra Stata packages, use:
*   ssc install <packagename>, replace

********************************************************************************
* 2) "Libraries import" (Python cell 5)
********************************************************************************
* Stata relies on built-in commands + optional SSC packages.

********************************************************************************
* 3) Statistical Significance labelling (Python cell 6)
********************************************************************************
capture program drop significance_stars
program define significance_stars, rclass
    * Usage: significance_stars, p(<pvalue>)
    syntax , p(real)
    local stars ""
    if (`p' < 0.001) local stars "***"
    else if (`p' < 0.01) local stars "**"
    else if (`p' < 0.05) local stars "*"
    return local stars "`stars'"
end

********************************************************************************
* 4) Suppress warnings (Python cell 7)
********************************************************************************
* Already handled by "set more off". No direct equivalent of warnings.filterwarnings.

********************************************************************************
* 5) Data loading (Python cell 8)
********************************************************************************
import delimited using "$IN_CSV", clear varnames(1)

* The notebook converts "year" to datetime and uses it as time index.
* In Part 1, "year_date" was created as Jan-1 of year; exports may contain:
* - year (string like 1985-01-01 or 1985)
* - year_int (numeric)
* We'll robustly create a numeric time variable year_int.

capture confirm variable year_int
if _rc {
    capture confirm numeric variable year
    if !_rc {
        gen int year_int = year
    }
    else {
        * Try to parse common date formats from "year" string.
        capture confirm string variable year
        if !_rc {
            gen double __tmpdate = date(year,"YMD")
            replace __tmpdate = date(year,"DMY") if missing(__tmpdate)
            replace __tmpdate = date(year,"MDY") if missing(__tmpdate)
            * If still missing, try destring
            capture destring year, gen(__tmpyearnum) ignore(" -/") force
            if _rc==0 {
                replace __tmpdate = mdy(1,1,__tmpyearnum) if missing(__tmpdate) & !missing(__tmpyearnum)
                drop __tmpyearnum
            }
            gen int year_int = year(__tmpdate)
            drop __tmpdate
        }
        else {
            di as error "Could not find year or year_int variable in CSV."
            exit 498
        }
    }
}
label var year_int "Year (numeric for panel time)"

* Country id
capture confirm string variable country
if _rc {
    tostring country, replace
}
encode country, gen(country_id)
label var country_id "Country id"

* Panel settings
xtset country_id year_int

********************************************************************************
* 6) PART 2A - Classic benchmark multivariate panel estimators
*    (Python cells 11-15 show 5 model variants; many estimation lines were
*     elided by "..." in the notebook. We implement the standard Stata analogs.)
********************************************************************************

* Helper: create Mundlak (RE with means) variables for a set of regressors
capture program drop mundlak_means
program define mundlak_means
    syntax varlist, prefix(name)
    foreach x of local varlist {
        capture drop `prefix'`x'
        by country_id: egen double `prefix'`x' = mean(`x')
        label var `prefix'`x' "Mundlak mean of `x' (by country)"
    }
end

* A small helper to run and store the 5 estimators for a given model
capture program drop run_benchmark_estimators
program define run_benchmark_estimators
    * args: depvar indepvars modelname
    syntax , depvar(name) indepvars(string asis) modelname(string)

    local y "`depvar'"
    local X "`indepvars'"

    di as txt "============================================================"
    di as txt "Model: `modelname'   Dependent: `y'   Regressors: `X'"
    di as txt "============================================================"

    * Between estimator
    quietly xtreg `y' `X', be vce(robust)
    estimates store `modelname'_between

    * One-way FE (within)
    quietly xtreg `y' `X', fe vce(robust)
    estimates store `modelname'_owfe

    * Random effects with Mundlak correction
    mundlak_means `X', prefix(m_)
    quietly xtreg `y' `X' m_*, re vce(robust)
    estimates store `modelname'_re_mundlak

    * Two-way FE: entity FE + year FE dummies
    quietly xtreg `y' `X' i.year_int, fe vce(robust)
    estimates store `modelname'_twfe

    * First differences estimator
    quietly xtreg `y' `X', fd vce(robust)
    estimates store `modelname'_fd

    * Comparison table (base Stata)
    estimates table `modelname'_between `modelname'_owfe `modelname'_re_mundlak `modelname'_twfe `modelname'_fd, ///
        b(%10.4f) se(%10.4f) stats(N r2, fmt(%9.0g %9.3f) labels("N" "R2")) ///
        title("Benchmark estimators - `modelname'") star

    * Clean up Mundlak means to avoid clashes across models
    drop m_*
end

********************************************************************************
* MODEL 1 (Python cell 11)
* y = GDP ; X = UN P WR LS i
* NOTE: Python comments mention PCOM excluded, and other variables excluded to
* avoid multicollinearity. The actual listed X-vars are below.
********************************************************************************
run_benchmark_estimators, depvar(GDP) indepvars("UN P WR LS i") modelname("M1")

********************************************************************************
* MODEL 2 (Python cell 12)
* y = GDP ; X = UN P W LS i   (nominal comp W instead of WR)
* NOTE (difficulty): The notebook later loops over ["UN","P","WR","LS","i"] when
* building lag2 instruments even though X_vars uses W. We keep the Model 2 X_vars
* as stated in the comments.
********************************************************************************
run_benchmark_estimators, depvar(GDP) indepvars("UN P W LS i") modelname("M2")

********************************************************************************
* MODEL 3 (Python cell 13)
* y = GDP ; X = SHORTUN P WR LS i
********************************************************************************
run_benchmark_estimators, depvar(GDP) indepvars("SHORTUN P WR LS i") modelname("M3")

********************************************************************************
* MODEL 4 (Python cell 14)
* y = GDP ; X = LONGUN P WR LS i
********************************************************************************
run_benchmark_estimators, depvar(GDP) indepvars("LONGUN P WR LS i") modelname("M4")

********************************************************************************
* MODEL 5 (Python cell 15)
* y = GDP ; X = LF P WR LS i
********************************************************************************
run_benchmark_estimators, depvar(GDP) indepvars("LF P WR LS i") modelname("M5")

********************************************************************************
* 7) Autocorrelation check in first-difference residuals (Python cell 17)
********************************************************************************
* Python: estimate FD model for Model 1, get residuals, regress resid on lag(resid).
preserve
    keep if !missing(GDP, UN, P, WR, LS, i)
    xtset country_id year_int
    quietly xtreg GDP UN P WR LS i, fd vce(robust)
    predict double resid_fd, e
    gen double resid_fd_L1 = L1.resid_fd
    drop if missing(resid_fd, resid_fd_L1)
    regress resid_fd resid_fd_L1
    estimates store ar1_resid_fd
restore

********************************************************************************
* 8) PART 2B - Dynamic model (Anderson-Hsiao style) on Model 1
*    (Python cells 18-24, 26, 28, 30-34, 36-40)
********************************************************************************

* Ensure required first-difference variables exist.
* Part 1 export likely already contains *_1diff. If not, create from levels.
foreach v in GDP UN P WR LS i W SHORTUN LONGUN LF {
    capture confirm variable `v'
    if _rc continue
    capture confirm variable `v'_1diff
    if _rc {
        gen double `v'_1diff = D.`v'
        label var `v'_1diff "First difference of `v' (created in Part 2)"
    }
}

* Lagged first differences (lag 1)
foreach v in GDP UN P WR LS i W SHORTUN LONGUN LF {
    capture confirm variable `v'_1diff
    if _rc continue
    capture drop `v'_1diff_lag1
    gen double `v'_1diff_lag1 = L1.`v'_1diff
    label var `v'_1diff_lag1 "Lag 1 of `v'_1diff"
}

* Lag 1 of GDP_1diff specifically (Python: GDP_1diff_lag1)
capture confirm variable GDP_1diff_lag1
if _rc {
    gen double GDP_1diff_lag1 = L1.GDP_1diff
    label var GDP_1diff_lag1 "Lag 1 of GDP_1diff"
}

* Instruments in levels lagged twice (t-2): GDP_lag2, P_lag2, etc.
foreach v in GDP UN P WR LS i W SHORTUN LONGUN LF {
    capture confirm variable `v'
    if _rc continue
    capture drop `v'_lag2
    gen double `v'_lag2 = L2.`v'
    label var `v'_lag2 "Lag 2 of `v' (level)"
}

* Create analysis sample as Python: diff_ah_df = df.dropna()
* Here: keep obs with all needed vars non-missing.
preserve
    keep country country_id year_int ///
         GDP_1diff GDP_1diff_lag1 GDP_lag2 ///
         P_1diff P_1diff_lag1 P_lag2 ///
         UN_1diff WR_1diff LS_1diff i_1diff ///
         UN_1diff_lag1 WR_1diff_lag1 LS_1diff_lag1 i_1diff_lag1
    drop if missing(GDP_1diff, GDP_1diff_lag1, GDP_lag2, ///
                    P_1diff, P_1diff_lag1, P_lag2, ///
                    UN_1diff, WR_1diff, LS_1diff, i_1diff)

    tempfile diff_ah
    save `diff_ah', replace
restore

use `diff_ah', clear
xtset country_id year_int

********************************************************************************
* 8.1) Univariate descriptive stats for target variables (Python cell 22)
********************************************************************************
summarize GDP_1diff GDP_1diff_lag1 GDP_lag2 ///
          P_1diff P_1diff_lag1 P_lag2 ///
          UN_1diff WR_1diff LS_1diff i_1diff, detail

* Skewness / kurtosis (Stata: summarize, detail provides skewness/kurtosis)
* If you need a separate table, you can loop and post results.

********************************************************************************
* 8.2) Correlations (Python cell 24 heatmap)
********************************************************************************
corr GDP_1diff GDP_1diff_lag1 GDP_lag2 ///
     P_1diff P_1diff_lag1 P_lag2 ///
     UN_1diff WR_1diff LS_1diff i_1diff
pwcorr GDP_1diff GDP_1diff_lag1 GDP_lag2 ///
       P_1diff P_1diff_lag1 P_lag2 ///
       UN_1diff WR_1diff LS_1diff i_1diff, sig obs

********************************************************************************
* 8.3) Unit-root testing (Python cell 26: manual Fisher-ADF)
********************************************************************************
* Stata has built-in Fisher-type panel unit root tests.
* This is the direct Stata equivalent of the intended Fisher-ADF logic.
xtunitroot fisher GDP_1diff, dfuller lags(1)
xtunitroot fisher P_1diff,   dfuller lags(1)

********************************************************************************
* 8.4) OLS dynamic regression using lagged first differences (Python cell 28)
********************************************************************************
regress GDP_1diff GDP_1diff_lag1 P_1diff_lag1 UN_1diff_lag1 WR_1diff_lag1 LS_1diff_lag1 i_1diff_lag1, vce(robust)
estimates store ols_dyn

********************************************************************************
* 8.5) OLS benchmark used for OLS vs IV comparison (Python cell 30)
********************************************************************************
regress GDP_1diff GDP_1diff_lag1 P_1diff P_1diff_lag1 UN_1diff WR_1diff LS_1diff i_1diff, vce(robust)
estimates store ols

********************************************************************************
* 8.6) IVREG (2SLS) with instruments in levels (Python cell 30)
* Endogenous: GDP_1diff_lag1
* Instruments: GDP_lag2 P_lag2
********************************************************************************
ivregress 2sls GDP_1diff (GDP_1diff_lag1 = GDP_lag2 P_lag2) ///
    P_1diff P_1diff_lag1 UN_1diff WR_1diff LS_1diff i_1diff, vce(robust)
estimates store iv

********************************************************************************
* 8.7) Parameter deviation OLS vs IV (Python cell 31)
********************************************************************************
* We build a small dataset comparing common coefficients.
tempname b_ols b_iv
matrix `b_ols' = e(b)
* Save IV first, then reload OLS/IV matrices safely:
estimates restore ols
matrix `b_ols' = e(b)
estimates restore iv
matrix `b_iv'  = e(b)

* Collect coefficient names intersection
local cols_ols : colnames `b_ols'
local cols_iv  : colnames `b_iv'

* Create dataset of comparisons
preserve
    clear
    set obs 0
    gen str64 term = ""
    gen double b_ols = .
    gen double b_iv  = .
    gen double abs_diff = .
    gen double pct_diff = .

    foreach t of local cols_ols {
        * only keep if also in IV
        local iniv : list t in cols_iv
        if "`iniv'"=="0" continue

        set obs `=_N+1'
        replace term = "`t'" in L
        replace b_ols = `b_ols'[1,"`t'"] in L
        replace b_iv  = `b_iv'[1,"`t'"]  in L
        replace abs_diff = abs(b_iv - b_ols) in L
        replace pct_diff = 100*(b_iv - b_ols)/b_ols in L
    }

    order term b_ols b_iv abs_diff pct_diff
    list, abbrev(32)

    export delimited using "Outputs/ols_vs_iv_coeff_comparison.csv", replace
restore

********************************************************************************
* 8.8) IV relevance check (Python cell 33) + IV quality tests (Python cell 34)
********************************************************************************
* First stage regression (simple):
regress GDP_1diff_lag1 GDP_lag2 P_lag2, vce(robust)
test GDP_lag2 P_lag2

* Display first-stage diagnostics from Stata 2SLS:
* "first" shows first-stage results; "estat firststage" provides diagnostics.
ivregress 2sls GDP_1diff (GDP_1diff_lag1 = GDP_lag2 P_lag2) ///
    P_1diff P_1diff_lag1 UN_1diff WR_1diff LS_1diff i_1diff, vce(robust) first
estat firststage
estat endogenous

********************************************************************************
* 8.9) Long-run coefficient (Python cells 39-40)
* beta_LT = (beta_1 + beta_2) / (1 - beta_y)
********************************************************************************
* IV long-run coefficient
estimates restore iv
scalar beta_y_iv = _b[GDP_1diff_lag1]
scalar beta1_iv  = _b[P_1diff]
scalar beta2_iv  = _b[P_1diff_lag1]
scalar beta_LT_iv = (beta1_iv + beta2_iv) / (1 - beta_y_iv)
di as result "Long-run coefficient (IV): " %9.4f beta_LT_iv

* OLS long-run coefficient
estimates restore ols
scalar beta_y_ols = _b[GDP_1diff_lag1]
scalar beta1_ols  = _b[P_1diff]
scalar beta2_ols  = _b[P_1diff_lag1]
scalar beta_LT_ols = (beta1_ols + beta2_ols) / (1 - beta_y_ols)
di as result "Long-run coefficient (OLS): " %9.4f beta_LT_ols

********************************************************************************
* 8.10) Impulse Response Function (IRF) for horizons 1..4 (Python cells 36-37)
*
* NOTE (difficulty / reconstruction):
* The notebook's IRF computation lines were replaced by "...". Given the dynamic
* specification:
*   GDP_1diff_t = beta_y*GDP_1diff_{t-1} + beta1*P_1diff_t + beta2*P_1diff_{t-1} + ...
* a 1-unit shock to P_1diff at t=1 implies responses:
*   h=1: beta1
*   h=2: beta2 + beta_y*beta1
*   h=3: beta_y*(beta2 + beta_y*beta1) = beta_y^2*beta1 + beta_y*beta2
*   h=4: beta_y^3*beta1 + beta_y^2*beta2
* This is what we implement below for both IV and OLS.
********************************************************************************

* IV IRF
estimates restore iv
scalar by = _b[GDP_1diff_lag1]
scalar b1 = _b[P_1diff]
scalar b2 = _b[P_1diff_lag1]

preserve
    clear
    set obs 4
    gen int horizon = _n
    gen double irf_iv = .
    replace irf_iv = b1                         if horizon==1
    replace irf_iv = b2 + by*b1                 if horizon==2
    replace irf_iv = (by^2)*b1 + by*b2          if horizon==3
    replace irf_iv = (by^3)*b1 + (by^2)*b2      if horizon==4

    twoway (line irf_iv horizon), ///
        yline(0, lpattern(dash)) ///
        title("IV Model - IRF: 1-unit increase in P_1diff") ///
        xtitle("Horizon (t)") ytitle("Impulse response of GDP_1diff")
    graph export "Graphs/irf_iv_Pshock.png", replace
restore

* OLS IRF
estimates restore ols
scalar by = _b[GDP_1diff_lag1]
scalar b1 = _b[P_1diff]
scalar b2 = _b[P_1diff_lag1]

preserve
    clear
    set obs 4
    gen int horizon = _n
    gen double irf_ols = .
    replace irf_ols = b1                         if horizon==1
    replace irf_ols = b2 + by*b1                 if horizon==2
    replace irf_ols = (by^2)*b1 + by*b2          if horizon==3
    replace irf_ols = (by^3)*b1 + (by^2)*b2      if horizon==4

    twoway (line irf_ols horizon), ///
        yline(0, lpattern(dash)) ///
        title("OLS Model - IRF: 1-unit increase in P_1diff") ///
        xtitle("Horizon (t)") ytitle("Impulse response of GDP_1diff")
    graph export "Graphs/irf_ols_Pshock.png", replace
restore

