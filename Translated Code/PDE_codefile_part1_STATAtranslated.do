
********************************************************************************
* PDE_codefile_part1_translated.do
********************************************************************************

version 16.0
clear all
set more off
set linesize 255

********************************************************************************
* 0) Paths / inputs
********************************************************************************
* NOTE: Adjust these paths if your folder structure differs.
* Python used: "Data/Dataset_MP_Impact_functional_Distribution.xlsx"
global DATA_XLSX  "Data/Dataset_MP_Impact_functional_Distribution.xlsx"
global OUT_CSV    "Data/final_trans_df.csv"
global OUT_XLSX   "Data/final_trans_df.xlsx"

********************************************************************************
* 1) "Requirements.txt installation" (Python cell 4)
********************************************************************************
* In Stata there is no pip-install step. If you need extra Stata packages:
* ssc install <packagename>, replace

********************************************************************************
* 2) "Libraries import" (Python cell 5)
********************************************************************************
* In Stata, functionality is provided by built-in commands + optional SSC pkgs.
* We'll rely on base Stata commands (import excel, egen, xtset, corr, regress, etc.)

********************************************************************************
* 3) Statistical Significance labelling (Python cell 6)
********************************************************************************
capture program drop significance_stars
program define significance_stars, rclass
    * Usage:
    *   significance_stars, p(<pvalue>)
    * Returns r(stars) = "***", "**", "*", or ""
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
* Stata: optionally reduce output. We already use "set more off".

********************************************************************************
* 5) Data loading (Python cell 10)
********************************************************************************
capture confirm file "$DATA_XLSX"
if _rc {
    di as error "Input file not found: $DATA_XLSX"
    di as error "Please place the Excel file in the expected folder or edit global DATA_XLSX."
    exit 601
}

import excel using "$DATA_XLSX", firstrow clear

* Python converted year to datetime to avoid including it as numeric in describe().
* In Stata we keep a numeric year for panel/time ops, and also create a daily date.
capture confirm numeric variable year
if _rc {
    destring year, replace ignore(" ")
}
gen int year_int = year
label var year_int "Year (numeric)"
gen double year_date = mdy(1,1,year_int)
format year_date %td
label var year_date "Year (Stata date Jan-1)"

* Ensure country is string (as in Python). If country is numeric, decode if labeled.
capture confirm string variable country
if _rc {
    tostring country, replace
}
label var country "Country"

********************************************************************************
* 6) Variable labels (Python cell 11: labels_mapper)
********************************************************************************
* Python renamed columns for plots; Stata cannot use spaces in variable names, so
* we apply VARIABLE LABELS instead (these show in graphs and tables).
capture label var i        "Short-Term Interest Rate"
capture label var P        "GDP Deflator"
capture label var W        "Nominal Compensation Per Employee"
capture label var WR       "Real Compensation Per Employee"
capture label var GDP      "Real Gross Domestic Product"
capture label var LS       "Adjusted Labour Share"
capture label var PCOM     "Energy Commodities Price Index"
capture label var UN       "Unemployment Rate"
capture label var SHORTUN  "Short-Term Unemployment"
capture label var LONGUN   "Long-Term Unemployment"
capture label var LF       "Labor Force"
capture label var REER     "Real Effective Exchange Rate"
capture label var SH       "Shadow Interest Rate"

********************************************************************************
* 7) Descriptive statistics (Python cells 12-13)
********************************************************************************
* Total dataset statistics (similar to pandas .describe())
ds year country year_int year_date, not
local numvars `r(varlist)'
di as txt "===================="
di as txt "Overall descriptive statistics"
di as txt "===================="
summarize `numvars', detail

* Single country statistics
di as txt "===================="
di as txt "Descriptive statistics by country"
di as txt "===================="
levelsof country, local(ctrys)
foreach c of local ctrys {
    di as txt "=== `c' ==="
    summarize `numvars' if country=="`c'", detail
}

********************************************************************************
* 8) Pre-plotting adjustments (Python cell 14)
********************************************************************************
* Python created target_countries and a color map for Matplotlib.
* In Stata we will store the target country list. Color mapping is not ported
* 1:1 (Stata schemes differ); you may customize graph options later.
local target_countries "Australia Germany Italy Japan UK USA"

********************************************************************************
* 9) Over-time plotting (Python cell 15)
********************************************************************************
* Python looped each variable and plotted time series per target country.
* Stata equivalent: twoway line per country. We generate graphs only if desired.
* NOTE (difficulty): A generic loop creating one multi-line graph per variable
* can be heavy; you can comment out this block if you only need data outputs.

capture mkdir "Graphs"

foreach v of local numvars {
    * Skip if variable does not exist
    capture confirm variable `v'
    if _rc continue

    * Keep only target countries for the graph
    preserve
        keep if inlist(country, "Australia","Germany","Italy","Japan","UK","USA")
        sort country year_int
        * Create separate lines per country using by(country) or overlay lines.
        * Overlay approach:
        twoway ///
            (line `v' year_int if country=="Australia") ///
            (line `v' year_int if country=="Germany") ///
            (line `v' year_int if country=="Italy") ///
            (line `v' year_int if country=="Japan") ///
            (line `v' year_int if country=="UK") ///
            (line `v' year_int if country=="USA"), ///
            legend(order(1 "Australia" 2 "Germany" 3 "Italy" 4 "Japan" 5 "UK" 6 "USA")) ///
            title("`v' evolution over time") xtitle("Year") ytitle("`v'")
        graph export "Graphs/ts_`v'.png", replace
    restore
}

********************************************************************************
* 10) Boxplots (Python cell 16)
********************************************************************************
* Python created per-variable boxplots by country.
* Stata: graph box var, over(country)
foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue
    preserve
        keep if inlist(country, "Australia","Germany","Italy","Japan","UK","USA")
        graph box `v', over(country) title("`v' - distribution by country")
        graph export "Graphs/box_`v'.png", replace
    restore
}

********************************************************************************
* 11) Deviation from cross-sectional mean plotting & boxplots (Python cells 17-18)
********************************************************************************
* NOTE (difficulty): Python constructed time-varying cross-sectional means by year
* and computed each country's deviation from the contemporaneous mean, plus a
* time-varying dispersion measure. We'll compute deviations; plots are optional.

* Ensure panel identifiers
encode country, gen(country_id)
xtset country_id year_int

foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue
    by year_int: egen double `v'_csmean = mean(`v')
    gen double `v'_dev_csmean = `v' - `v'_csmean
    label var `v'_csmean "Cross-sectional mean of `v' (by year)"
    label var `v'_dev_csmean "Deviation from cross-sectional mean of `v' (by year)"
}

********************************************************************************
* 12) Between transformation (Python cell 24)
********************************************************************************
foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue
    by country_id: egen double `v'_between = mean(`v')
    label var `v'_between "Between transformation: mean over country of `v'"
}

********************************************************************************
* 13) One-way within (demeaning) transformation (Python cell 25)
********************************************************************************
foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue
    gen double `v'_within = `v' - `v'_between
    label var `v'_within "Within transformation: `v' - country mean(`v')"
}

********************************************************************************
* 14) Variance analysis (Python cell 26)
********************************************************************************
* Python intended: var_pooled = Var(x_it), var_between = Var(x_i.), var_within = Var(x_it-x_i.)
* We'll compute these variances and store them in a Stata dataset.
tempfile vartable
preserve
    clear
    set obs 0
    gen str64 Variable = ""
    gen double var_pooled  = .
    gen double var_between = .
    gen double var_within  = .
    gen double share_between = .
    gen double share_within  = .

    tempfile tmpstats
    save `vartable', replace
restore

foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue

    quietly summarize `v'
    * If all missing, skip
    if r(N)==0 continue

    quietly summarize `v', detail
    local vp = r(Var)

    quietly summarize `v'_between, detail
    local vb = r(Var)

    quietly summarize `v'_within, detail
    local vw = r(Var)

    local sb = .
    local sw = .
    if (`vp'!=0 & `vp'!=.) {
        local sb = `vb'/`vp'
        local sw = `vw'/`vp'
    }

    preserve
        use `vartable', clear
        set obs `=_N+1'
        replace Variable      = "`v'" in L
        replace var_pooled    = `vp' in L
        replace var_between   = `vb' in L
        replace var_within    = `vw' in L
        replace share_between = `sb' in L
        replace share_within  = `sw' in L
        save `vartable', replace
    restore
}

use `vartable', clear
order Variable var_pooled var_between var_within share_between share_within
list, abbrev(20)

* 3D variance surface plotting (Python cell 27)
* NOTE (difficulty): Plotly 3D surface has no direct Stata base equivalent.
* You can export vartable and plot externally, or use Stata's graph3d (user-written)
* if installed. We omit the 3D surface.

********************************************************************************
* 15) First difference transformation (Python cell 36)
********************************************************************************
* Requires panel ts settings (already xtset).
foreach v of local numvars {
    capture confirm variable `v'
    if _rc continue
    gen double `v'_1diff = D.`v'
    label var `v'_1diff "First difference of `v'"
}

********************************************************************************
* 16) Balanced panel set-up (Python cell 45) + TWFE (Python cell 46)
********************************************************************************
* Python defined balanced_raw_df = raw_df[year >= 1985].
preserve
    keep if year_int >= 1985

    * Two-way fixed effects (TWFE) transformation:
    * x_it - x_i. - x_.t + x_..
    foreach v of local numvars {
        capture confirm variable `v'
        if _rc continue

        by country_id: egen double `v'_i_mean = mean(`v')
        by year_int:   egen double `v'_t_mean = mean(`v')
        egen double `v'_all_mean = mean(`v')

        gen double `v'_TWFE = `v' - `v'_i_mean - `v'_t_mean + `v'_all_mean
        label var `v'_TWFE "TWFE transformed `v' (balanced sample year>=1985)"
    }

    tempfile TWFE_balanced
    save `TWFE_balanced', replace
restore

********************************************************************************
* 17) Time component plotting and TWFE plotting (Python cells 47-48)
********************************************************************************
* NOTE (difficulty): Python created an explicit (-x_.t + x_..) component and plotted.
* In Stata, that component is (`v'_all_mean - `v'_t_mean). We compute it if needed.
preserve
    use `TWFE_balanced', clear
    foreach v of local numvars {
        capture confirm variable `v'_all_mean
        if _rc continue
        gen double `v'_timecomp = -`v'_t_mean + `v'_all_mean
        label var `v'_timecomp "Time component for `v': -x_.t + x_.."
    }
    tempfile TWFE_balanced2
    save `TWFE_balanced2', replace
restore

********************************************************************************
* 18) TWFE balanced panel descriptive stats / boxplots / correlations (Python 50-55)
********************************************************************************
preserve
    use `TWFE_balanced2', clear
    * Descriptive for TWFE variables
    local twfevars ""
    foreach v of local numvars {
        capture confirm variable `v'_TWFE
        if _rc continue
        local twfevars `twfevars' `v'_TWFE
    }
    summarize `twfevars', detail

    * By country
    levelsof country, local(ctrys2)
    foreach c of local ctrys2 {
        di as txt "=== `c' (TWFE) ==="
        summarize `twfevars' if country=="`c'", detail
    }

    * Boxplots (may be many)
    foreach v of local twfevars {
        graph box `v', over(country) title("`v' - TWFE distribution by country (year>=1985)")
        graph export "Graphs/box_TWFE_`v'.png", replace
    }

    * Correlation matrix (heatmap not replicated; use corr/pwcorr)
    corr `twfevars'
    * If you want significance:
    pwcorr `twfevars', sig obs
restore

********************************************************************************
* 19) Unbalanced panel TWFE (Python cells 57-58: Wansbeek–Kapteyn)
********************************************************************************
* NOTE (difficulty): Wansbeek–Kapteyn (1989) provides a specific transformation
* for TWFE in unbalanced panels. In Stata, an exact reproduction is non-trivial
* without custom code. A practical equivalent for obtaining the TWFE "purged"
* variable is: regress v on i.country i.year and use residuals.
* This yields a two-way FE residualized series that matches the idea of TWFE
* demeaning under unbalancedness.
preserve
    * Start from full sample (unbalanced)
    foreach v of local numvars {
        capture confirm variable `v'
        if _rc continue
        quietly regress `v' i.country_id i.year_int if !missing(`v')
        predict double `v'_TWFEu, resid
        label var `v'_TWFEu "TWFE residual of `v' from i.country + i.year (unbalanced approx.)"
        drop _hat
    }
    tempfile TWFE_unbalanced
    save `TWFE_unbalanced', replace
restore

********************************************************************************
* 20) TWFE unbalanced descriptive stats (Python cells 60, 62-63)
********************************************************************************
preserve
    use `TWFE_unbalanced', clear
    local twfeuvars ""
    foreach v of local numvars {
        capture confirm variable `v'_TWFEu
        if _rc continue
        local twfeuvars `twfeuvars' `v'_TWFEu
    }
    summarize `twfeuvars', detail

    levelsof country, local(ctrys3)
    foreach c of local ctrys3 {
        di as txt "=== `c' (TWFE unbalanced residuals) ==="
        summarize `twfeuvars' if country=="`c'", detail
    }
restore

********************************************************************************
* 21) Correlation heatmaps for Between / Within / First-diff etc. (Python 68, 70, 71, 79, 82, 83)
********************************************************************************
* NOTE: Stata does not natively plot Seaborn-style correlation heatmaps.
* We compute correlation matrices using corr/pwcorr.

* Between transformed (Python cell 68)
local betweenvars ""
foreach v of local numvars {
    capture confirm variable `v'_between
    if _rc continue
    local betweenvars `betweenvars' `v'_between
}
corr `betweenvars'
pwcorr `betweenvars', sig obs

********************************************************************************
* 22) Within-transformed lag extraction (Python cells 69, 73)
********************************************************************************
* Python created lag1, lag2, lag5, lag10 (and sometimes all lags to 10).
* Stata: after xtset, use L#.var.
* We create explicit lag variables matching the Python naming.
foreach v of local numvars {
    capture confirm variable `v'_within
    if _rc continue
    gen double `v'_within_Lag1  = L1.`v'_within
    gen double `v'_within_Lag2  = L2.`v'_within
    gen double `v'_within_Lag5  = L5.`v'_within
    gen double `v'_within_Lag10 = L10.`v'_within
}

* Create trend within country (Python used cumcount + 1)
by country_id (year_int): gen long trend = _n
label var trend "Within-country time trend (1..T)"

* "All lags" up to 10 for within transformed variables
forvalues L=1/10 {
    foreach v of local numvars {
        capture confirm variable `v'_within
        if _rc continue
        capture drop `v'_within_Lag`L'
        gen double `v'_within_Lag`L' = L`L'.`v'_within
    }
}

********************************************************************************
* 23) Country-specific detrending for Model 1 variables (Python cells 74-75)
********************************************************************************
* Python detrended within variables via country-specific linear trend (polyfit degree=1).
* Stata: run country-specific regressions on trend and store residuals.

* Define "Model 1" target vars (from notebook):
local model1 "i UN P GDP WR LS"

foreach base of local model1 {
    capture confirm variable `base'_within
    if _rc continue

    gen double `base'_within_detrended = .
    levelsof country_id, local(ids)
    foreach id of local ids {
        quietly regress `base'_within trend if country_id==`id' & !missing(`base'_within, trend)
        * If not enough obs, skip
        if e(N)<2 continue
        tempvar r
        predict double `r' if country_id==`id', resid
        replace `base'_within_detrended = `r' if country_id==`id'
        drop `r'
    }
    label var `base'_within_detrended "`base'_within detrended by country-specific linear trend"
}

* Lags (all lags up to 10) for detrended variables (Python cell 75)
forvalues L=1/10 {
    foreach base of local model1 {
        capture confirm variable `base'_within_detrended
        if _rc continue
        gen double `base'_within_detrended_Lag`L' = L`L'.`base'_within_detrended
    }
}

********************************************************************************
* 24) First differences lag extraction (Python cells 80-81)
********************************************************************************
* Create lags for first differences for the same Model 1 set (and more if needed)
foreach base of local model1 {
    capture confirm variable `base'_1diff
    if _rc continue
    forvalues L=1/10 {
        gen double `base'_1diff_Lag`L' = L`L'.`base'_1diff
    }
}

********************************************************************************
* 25) Bivariate graphs with linear, quadratic and lowess fit (Python cell 85)
********************************************************************************
* Python produced many bivariate plots. We'll implement one representative set:
* GDP_1diff vs i_1diff (as in Python cell 41) with lfit, qfit, and lowess.
capture confirm variable GDP_1diff
capture confirm variable i_1diff
if !_rc {
    twoway ///
        (scatter GDP_1diff i_1diff, msymbol(o)) ///
        (lfit    GDP_1diff i_1diff) ///
        (qfit    GDP_1diff i_1diff) ///
        (lowess  GDP_1diff i_1diff), ///
        title("GDP_1diff vs i_1diff (linear, quadratic, lowess)") ///
        xtitle("i_1diff") ytitle("GDP_1diff")
    graph export "Graphs/bivar_GDP1diff_i1diff.png", replace
}

********************************************************************************
* 26) Final transformed dataset export (Python cell 86)
********************************************************************************
* Python merged: trans_df (within/between), trans_diff_df (1diff), TWFE_unbalanced_df.
* In Stata we already have within/between/1diff in the current dataset, and we saved
* TWFE_unbalanced to a tempfile. We merge on (year_int, country).

tempfile base_full
save `base_full', replace

* Merge unbalanced TWFE residuals
preserve
    use `TWFE_unbalanced', clear
    keep year_int country country_id
    foreach v of local numvars {
        capture confirm variable `v'_TWFEu
        if _rc continue
        keep `v'_TWFEu
    }
restore

use `base_full', clear
merge 1:1 country_id year_int using `TWFE_unbalanced', nogen

* Sort and export
sort country_id year_int

* Export CSV
export delimited using "$OUT_CSV", replace

* Export Excel
export excel using "$OUT_XLSX", firstrow(variables) replace

di as result "Exported final transformed dataset to:"
di as result " - $OUT_CSV"
di as result " - $OUT_XLSX"

