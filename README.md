# simplydata_tests

> simplydata: **Learn basic Python, and simplydata takes you from your first table to publication-ready statistics, in plain words, with every number checked against the field's standard tools.**

Jupyter notebooks that test the Python package [simplydata](https://github.com/bioforkgithub/simplydata)
on published reference datasets, side by side with the packages it calls underneath (pandas, NumPy,
SciPy, statsmodels) and with independent programs (R and its field-standard packages, lifelines). Every number is compared with an
explicit tolerance, and each comparison says what it found, so the results can be cited:

- **agrees**: within the tolerance;
- **simplydata error – being fixed**: a simplydata bug (FINDINGS.md class A), until the fix is released;
- **explained difference – reference**: the reference rounds, stops early, uses another definition or has its own
  known bug (class C; the explanation names which);
- **explained difference – machine**: the last digits differ between CPUs or systems;
- **explained difference – precision limit**: the limit of double precision or of the package simplydata calls, the
  same on every system (class B; and C29, a reference that stopped at R's default convergence);
- **difference not yet explained**: outside the tolerance, and FINDINGS.md does not explain it yet (listed there);
- **waiting for word** and **skipped**: no comparison could be made yet (the reason is given).

Every difference keeps its finding id, a one-line explanation and a link to [FINDINGS.md](FINDINGS.md); the results
files keep the machine column "status" (PASS, FAIL, SKIP) beside the reader's column "result", so nothing is hidden.

The notebooks are saved with their outputs, so they can be read on GitHub without running anything. The outputs
are those of the GitHub Actions run of the pull request (Ubuntu); the stored results files name in their column
"ran on" the system that made them.

| notebook | what is checked | data |
|---|---|---|
| [00_summary](notebooks/00_summary.ipynb) | every comparison in one table: by notebook, by kind of reference, every difference with its explanation, every check waiting for a word or skipped | results of 01 to 32 |
| [01_describing_data](notebooks/01_describing_data.ipynb) | `find_average_of`, `find_standard_deviation_of`, `find_autocorrelation_of`, `find_median_of` | NIST StRD univariate (9 datasets, certified) |
| [02_comparing_groups](notebooks/02_comparing_groups.ipynb) | `compare_averages_of` (Welch, ANOVA, Welch ANOVA, Kruskal-Wallis), `compare_before_and_after`, `compare_spreads_of` | Student's sleep data; NIST StRD ANOVA (11, certified); iris; warpbreaks |
| [03_regression_NIST](notebooks/03_regression_NIST.ipynb) | `fit_line`, `fit_curve`, `explain` (least squares), `fit_formula` through 0 | NIST StRD linear regression (11, certified); Anscombe's quartet; mtcars; stackloss |
| [04_counts_and_yes_no](notebooks/04_counts_and_yes_no.ipynb) | `explain` (logistic, Poisson, negative binomial), `compare_counts_of`, `compare_risks` | birthwt; NMES 1988 |
| [05_survival](notebooks/05_survival.ipynb) | `find_survival_over_time` (Kaplan-Meier), `compare_survival` (log-rank), `explain_survival` (Cox) | NCCTG lung cancer |
| [06_curves_and_dose_response](notebooks/06_curves_and_dose_response.ipynb) | `fit_formula` (nonlinear least squares, both NIST starts), `explain(out_of=)`, `fit_dose_response` (LD50) | NIST StRD nonlinear (27, certified); Puromycin; Bliss's beetles |
| [07_distributions](notebooks/07_distributions.ipynb) | `distribution(...).find_chance_that_value`, `find_cut_off` | NIST/SEMATECH e-Handbook tables |
| [08_combining_studies](notebooks/08_combining_studies.ipynb) | `combine_studies` (REML, DL, Paule-Mandel, common effect, Mantel-Haenszel, Peto, Hartung-Knapp, subgroups, meta-regression, Egger, trim and fill) | metadat: BCG trials (Viechtbauer 2010, printed), passive smoking (Hackshaw 1997, printed), Normand, Curtis, Hart, McDaniel; against metafor and meta |
| [09_survival_and_person_time](notebooks/09_survival_and_person_time.ipynb) | survival added after 0.2.0: time points, restricted mean, life table, adjusted curves, late entry, competing risks, Gray's test, Fine-Gray, stratified and start-stop Cox, Weibull | lung, veteran, mgus2, heart, udca (survival); bfeed, channing, larynx (KMsurv); against R survival, cmprsk, survRM2, KMsurv, adjustedCurves, statsmodels, lifelines |
| [10_power_and_sample_size](notebooks/10_power_and_sample_size.ipynb) | `find_how_many_needed`, `find_chance_of_finding` for averages, percents, ANOVA, chi-square, correlation, regression, logistic, equivalence, survival, clusters | Cohen (1988) tables, G*Power 3.1 manual, R help pages, pwr vignette (printed); against statsmodels, SciPy, R pwr, TOSTER, presize, gsDesign, powerSurvEpi, pwrss |
| [11_communities](notebooks/11_communities.ipynb) | diversity, rarefaction, Chao1, Hill numbers, Bray-Curtis, PERMANOVA, spread check, ANOSIM, PCoA, NMDS, CA, db-RDA, CCA, Mantel, clusters, clr | vegan BCI, dune, varespec, mite (vegan tutorial, printed); against vegan, zCompositions, compositions, cluster, SciPy/NumPy formulas |
| [12_factor_analysis](notebooks/12_factor_analysis.ipynb) | confirmatory and structural models, FIML, exploratory factor analysis, parallel analysis, rotations, Cronbach's alpha | Holzinger & Swineford (1939) (Rosseel 2012, printed fit); against lavaan, psych, GPArotation, factanal, statsmodels |
| [13_circular_spatial_extremes](notebooks/13_circular_spatial_extremes.ipynb) | average direction, Rayleigh test, Moran's I, variograms, GEV and generalised Pareto fits, return levels | circular datasets; Columbus (Anselin 1988); Meuse; Port Pirie, Fremantle and daily rain (Coles 2001, printed); against SciPy, circular, ape, spdep, gstat, extRemes |
| [14_small_additions](notebooks/14_small_additions.ipynb) | Storey q-values, automatic ARIMA, Nemenyi and Conover, pairwise percents, dose-response, Box-Cox, line-fitting methods, covariance, MANOVA | Hedenfalk p-values (pi0 0.67), AirPassengers (Box & Jenkins), rounding times, Fleiss smokers, ryegrass, poisons (Box & Cox 1964), leaflife, creatinine, iris, mtcars; against qvalue, forecast, PMCMRplus, drc, car, smatr, mcr, deming, statsmodels |
| [15_surveys_and_designs](notebooks/15_surveys_and_designs.ipynb) | `use_survey_design` (means, totals, quantiles, percents, groups, models, tests), `make_design_table`, `make_random_assignments`, `keep_random_rows` | api stratified and cluster samples (Lumley 2010, printed); against R survey, textbook formulas, agricolae, blockrand, pandas (same seeds) |
| [16_cause_and_effect](notebooks/16_cause_and_effect.ipynb) | balanced groups by weighting and matching, mediation, attributable fractions, direct and indirect standardisation, crossover bioequivalence | lalonde; JOBS II; dry cat food and FLUTD; Down syndrome by birth order (Selvin 2004); Schütz, Labes & Fuglsang (2014) reference datasets; against WeightIt, MatchIt, cobalt, statsmodels, mediation, epiR, epitools, R lm |
| [17_bayesian](notebooks/17_bayesian.ipynb) | `bayesian=True` on tests and correlations, `explain(bayesian=True)` | sleep, iris, mtcars, race dolls (BayesFactor manual, printed), sleepstudy; against BayesFactor, Rouder et al. (2009) and Gunel & Dickey (1974) formulas in SciPy, rstanarm (same seed, and long runs within 4 Monte Carlo standard errors) |
| [18_summaries_real_data](notebooks/18_summaries_real_data.ipynb) | number words (`nums.`), table summaries by group and `every=`, intervals of averages, medians, percentiles and standard deviations, R's nine percentile rules, frequency and analytic weights, `find_summary_of` (every measure), bootstrap intervals, `count_each`, `find_percent_of_rows_where`, `find_correlation_between` (every method, one-sided, partial, shuffling, bootstrap) | NIST StRD univariate (certified); airquality; API schools (survey); R `cor.test` help page (Hollander & Wolfe); mtcars; against base R, DescTools, e1071, psych, WRS2, Hmisc, survey, boot, ppcor, SciPy (same seed) |
| [19_tests_averages_spreads_shapes](notebooks/19_tests_averages_spreads_shapes.ipynb) | `compare_average_to`, `compare_averages_of` (every `test_method=`, one-sided, shuffling, equivalence, repeated measures, Bayes factor with `prior_width=`), `compare_before_and_after`, `sdata.compare_averages`, `compare_spreads_of`, `compare_shapes_of`, `check_if_random_order`, `compare_shape_to` | SciPy docstring examples (McDonald's mussels, Darwin's maize, Brunner & Munzel, Mood, `alexandergovern`, `levene`, `mannwhitneyu`); R help-page examples (`ansari.test`, `mood.test`, `wilcox.test`, `mcnemar.test`, randtests, TOSTER, kSamples); morley, sleep, PlantGrowth, InsectSprays, warpbreaks, CO2, circular wind; against R, WRS2, onewaytests, PMCMRplus, perm, car, goftest, Rouder et al.'s JZS integral |
| [20_tests_counts_pairs_matched](notebooks/20_tests_counts_pairs_matched.ipynb) | `compare_counts_of` (Fisher, G, Barnard, Boschloo, trend, cells, exact and shuffled p-values), `sdata.compare_counts`, `compare_counts_to`, `compare_percent_to`, `compare_each_pair_of` (every follow-up test, Dunnett, letters), `compare_matched_columns` (Friedman, Quade, Page, Cochran's Q), `add_adjusted_p_value_column`, `compare_correlations` | R `fisher.test`, `prop.trend.test`, `binom.test`, `friedman.test`, DescTools `CochranQTest` help pages; SciPy `tukey_hsd`, `dunnett`, `page_trend_test` docstrings; HairEyeColor, PlantGrowth, InsectSprays, Hedenfalk p-values, mtcars; against R, DescTools, multcomp, PMCMRplus, qvalue, crank, and formulas in NumPy |
| [21_agreement_and_diagnostic_tests](notebooks/21_agreement_and_diagnostic_tests.ipynb) | `agreement_between` (Cohen's, weighted and Fleiss' kappa, every ICC form, Krippendorff, Lin), `limits_of_agreement`, `reliability_of`, `compare_to_truth` (yes/no tests and ROC), `sdata.compare_to_truth`, `compare_accuracy` | R irr (`diagnoses`, `anxiety`, Krippendorff's published example), pROC (`aSAH`) and epiR (`epi.tests`, Scott et al. 2008, its stated answer) help pages; creatinine (mcr); Holzinger & Swineford; against irr, psych, pROC, `binom.test` and the published formulas in NumPy |
| [22_risks_and_rates](notebooks/22_risks_and_rates.ipynb) | `compare_rates` and `sdata.compare_rates` (exact test and interval, one-sided, equivalence, non-inferiority, `per=`, `interval=`), `rate_of`, `compare_risks` and `sdata.compare_risks` (Newcombe, Wald, Fisher, one-sided, equivalence, Mantel-Haenszel), `compare_balanced_groups_of` (`against=`, `interval=`) | R `poisson.test` (Andersen 1977) and `mantelhaen.test` (UCB admissions) help pages; NCCTG lung (survival); lalonde; against R `poisson.test`, `prop.test`, `fisher.test`, `mantelhaen.test`, DescTools, and formulas in SciPy |
| [23_distributions](notebooks/23_distributions.ipynb) | every `dist.` word on 15 distributions (density, percentiles, cut-offs, middle range, moments, summary, parameters, random numbers, chart), fits with `fixed=`, `compare_distributions_for` | R `d`/`q` functions and circular; textbook moment formulas; airquality ozone and `discoveries`; against the likelihood equations in SciPy, MASS `fitdistr`, fitdistrplus (tight) |
| [24_charts_and_directions](notebooks/24_charts_and_directions.ipynb) | the numbers drawn by the chart words (density, cumulative share, Q-Q dots, every smooth line, error bars, curves, fitted lines, curves and distributions, dots, bars, pairs, `difference_dots`, `trade_off_curve`, richness curves, threshold check, variogram, dose response, cluster tree), `find_variogram_of`, `find_average_direction_of` by group | mtcars, airquality, Lake Huron, creatinine, aSAH, vegan BCI, Meuse, Coles's rain, ryegrass, USArrests, circular pigeons; against R `lowess`, `ecdf`, `lm`, `nls`, `smooth.spline`, `hclust`, `stats::filter`, circular, and formulas in NumPy |
| [25_models_outcomes](notebooks/25_models_outcomes.ipynb) | `explain` for counts (Poisson, quasi, negative binomial, zero-inflated, hurdle), yes/no with links, shares (beta), positive numbers (Gamma), ordered and unordered groups, percentiles, robust errors, weights, profile intervals, tables of tests, standardised effects, `add_predicted_column` | NMES 1988 (Zeileis, Kleiber & Jackman 2008, printed), Bliss's beetles (Dobson & Barnett 2008, printed), GasolineYield (Ferrari & Cribari-Neto 2004, printed), clotting (McCullagh & Nelder 1989, printed), engel (quantreg vignette, printed), stackloss (statsmodels docstring, printed), Moore, housing, bioChemists; against R stats, MASS, pscl, betareg, quantreg, sandwich, car, nnet, effectsize |
| [26_models_mixed_and_more](notebooks/26_models_mixed_and_more.ipynb) | `explain(grouped_by=)` (linear and yes/no mixed models, GEE, correlated errors), splines and smooth terms, multiple imputation, bootstrap, `find_model_averages_of` (contrasts, within groups, over a range), `compare_models`, cross-validation, penalised models | sleepstudy (Bates et al. 2015, printed), pigs and mtcars (emmeans vignette, printed), cbpp, dietox, Ovary, warpbreaks, housing, nhanes; against lme4, lmerTest, geepack, nlme, mgcv, mice, emmeans, AICcmodavg, boot, glmnet, SciPy (same seed) |
| [27_fitting](notebooks/27_fitting.ipynb) | `fit_formula` with every `test_method=` and weights, `fit_line` (Passing-Bablok, Siegel), power transforms, segmented lines, `fit_dose_response` (3, 4, 5 parameters, Weibull, binomial) | NIST StRD nonlinear (26 problems x 3 methods, certified) and Pontius; creatinine; ryegrass, selenium; against R nls, mcr, segmented, drc, SciPy |
| [28_patterns_across_columns](notebooks/28_patterns_across_columns.ipynb) | principal components, k-means and hierarchical groups, silhouettes, distances (Bray-Curtis, Jaccard, Aitchison), NMDS, PCoA, CCA, discriminant analysis, log-ratio columns, db-RDA, RDA, PERMANOVA, SIMPER, Mantel, Procrustes, exploratory factor analysis, richness | USArrests, iris, dune, varespec, Hydrochem, Holzinger & Swineford (Rosseel 2012, printed); against R stats, cluster, vegan, MASS, compositions, zCompositions, psych |
| [29_forecasts_smoothing_trend](notebooks/29_forecasts_smoothing_trend.ipynb) | `forecast` and `compare_forecasts_for` (exponential smoothing, seasonal naive, ARIMA), every smoothing method, trend and season, trend removal, autocorrelation, trend tests, unit-root tests, change points, Moran's I on the globe | maxau and nottem (trend vignette, printed), Cobb (1978) Nile change point, co2, USAccDeaths, lh, Nile, mdeaths, quakes; against R stats, forecast, signal, trend, urca, tseries, changepoint, spdep |
| [30_survival_more](notebooks/30_survival_more.ipynb) | time and event from dates, competing risks (Aalen-Johansen, Gray's test, cause-specific Cox, Fine-Gray), Weibull and log-normal models, restricted mean time, adjusted curves, Cox with multiple imputation | mgus2, lung, ovarian (survival); against R survival, cmprsk, survRM2, adjustedCurves, mice |
| [31_combining_planning_surveys](notebooks/31_combining_planning_surveys.ipynb) | `combine_studies` (more methods and effect sizes), power and sample size for more tests (same-seed simulation too), group-sequential and survival designs, survey designs and models, `make_design_table` (split plot), random assignments | BCG trials (Viechtbauer 2010, printed), pwr vignette (printed), api (Lumley 2010, printed); against metafor, stats, gsDesign, powerSurvEpi, pwr, pwrss, TOSTER, survey, agricolae, blockrand, SciPy (same seed) |
| [32_models_planning_more](notebooks/32_models_planning_more.ipynb) | the rows the earlier notebooks left open: `compare_models` for groups, ordered groups and shares, model averages for zero-inflated and Gamma models, crossed random effects, Aitchison PCoA, power for matched pairs and counts, power by simulation from distributions, CFA with empty cells kept, design tables (completely random, blocks, Latin square), Bayesian Poisson model | housing, GasolineYield, bioChemists, clotting, carrots, Hydrochem, bfi, warpbreaks, pwr vignette (printed); against MASS, nnet, betareg, lmtest, emmeans, lme4, lavaan, agricolae, rstanarm, SciPy (same seed) |

## What the comparisons found

The results files as stored here (the pull request's run 45443130321 on GitHub Actions, every notebook with simplydata main
[3f0192b](https://github.com/bioforkgithub/simplydata/commit/3f0192b207f63836c505e4da99e0c2bc8b7cbd9e), 2026-11-24, fifth pass):

| notebooks | comparisons | agrees | simplydata error – being fixed | explained difference – reference | explained difference – machine | explained difference – precision limit | difference not yet explained | waiting for word | skipped |
|---|---|---|---|---|---|---|---|---|---|
| 01 to 08 | 1884 | 1743 | 0 | 5 | 10 | 125 | 0 | 0 | 1 |
| 09 to 17 | 2732 | 2689 | 0 | 20 | 1 | 22 | 0 | 0 | 0 |
| 18 to 24 | 2368 | 2345 | 0 | 21 | 0 | 1 | 0 | 0 | 1 |
| 25 to 32 | 19734 | 19487 | 0 | 18 | 11 | 218 | 0 | 0 | 0 |
| all | 26718 | 26264 | 0 | 64 | 22 | 366 | 0 | 0 | 2 |

**Fourth pass (simplydata main a5a8ff0).** The 29 rows left as "simplydata error – being fixed" turned out not to be
simplydata errors, and none is left: the textbook Wald intervals of the attributable fractions against R epiR's, which
simplydata gives and names (C31, 4 rows, was A5); two exact percentile intervals that pick different order statistics
(C32, 10, was A8); R `smooth.spline`, which lands slightly off the penalised optimum that SciPy and simplydata reach (C33,
1, was A17; not because x repeats, as A17 said); and NIST Lanczos1 with trust region and dogbox, at the double-precision
limit of a 9e-14 residual standard deviation (B1, 14, were A21). 22 older A21 and A26 rows of notebook 27, kept from other
systems, pass on both Ubuntu runs of this pull request and on a Windows run on GitHub (the new workflow `run-windows`),
and were removed from `results/known_fails.csv`. See [FINDINGS.md](FINDINGS.md), "Fourth pass".

**Third pass (simplydata main a5a8ff0, PR #42).** The last 5 checks that waited for a word compare and agree: the
correspondence-analysis eigenvalues and total inertia (notebook 11), the partial correlation's p-value (18) and the best
cut-off of the ROC chart (24). The newer simplydata (PRs #37 to #42) fixed 71 rows that were worth fixing (A7 37, A26 16,
A1 6, A6 5, A21 5, A25 2). Its fits now run to convergence (A7), so the references run with R's default stopping rule
differ from it: each of them got a converged reference next to it, and the 77 where simplydata agrees with the converged
one to 1e-8 are explained (C29, precision limit). 20 drc comparisons (notebook 14 and the Weibull fit of 27) also differ
from drc at relTol 1e-12; an independent optimum (R optim and SciPy least_squares, 1e-15) shows that drc stops short of
the least-squares optimum in its parameters even at relTol 1e-15, and simplydata gives exactly what drc gives: by
Manish's decision of 2026-11-23 these 44, and 7 rows of notebook 14 that were A7, are explained difference – precision
limit (C30). The previous run (simplydata 3fd3b85, PR #39) is described below and
in [FINDINGS.md](FINDINGS.md), "Second pass" and "Third pass".

**Second pass (simplydata main 3fd3b85, PR #39).** The 153 checks that waited for a simplydata word (the numbers a chart
draws, the median survival time, the note of `fit_formula`, the two-column correlation table, the rows that pull a
model, the numbers of `check_model_of`, the NMDS stress and the correspondence-analysis shares) are switched on: 143
agree, 4 are differences already explained (A7, A17, C10, C28), 5 still wait for a word (the correspondence-analysis
eigenvalues and total inertia, the partial correlation's p-value, the cut-off written on the ROC chart) and 1 reference
row is no longer compared (see FINDINGS.md, A9). Many of them now compare a whole drawn layer point by point, so there
are more comparisons than before (23,438 on simplydata 47594a0). PR #39 also fixed A9 to A16, A20, A22, A23 and most
of A21; it brought two new points worth fixing (A25, A26), and R irr, kept as a second reference, now differs (C26, C27).
See [FINDINGS.md](FINDINGS.md), "Second pass".

No comparison is outside its tolerance without an explanation in [FINDINGS.md](FINDINGS.md). The "simplydata
error – being fixed" rows become "agrees" when the fixes are released and the notebooks run again.

**Results (simplydata main, commit [37afb81](https://github.com/bioforkgithub/simplydata/commit/37afb8135bee95a1e307d563f420fc6fc45ce194), after release 0.2.0; 2026-11-21, Windows 11): 5692 comparisons, 5383 agree, 295 outside the tolerance (each
explained), 14 skipped.**
On Ubuntu (the GitHub Actions run, OpenBLAS pinned as below) the same notebooks give 5695 comparisons, 5383 agree,
298 outside the tolerance, 14 skipped (5388 and 293 on the other kind of runner); NIST nonlinear problems near the
6-digit line and one Monte Carlo comparison differ between the systems. Notebooks 01 to 07 test the methods of 0.2.0
(on 0.2.0 they gave 2052 comparisons, 1842 agree, 194 outside the tolerance, 16 skipped); notebooks 08 to 17 test the
methods added after it. Every difference is explained in [FINDINGS.md](FINDINGS.md).
In short: three of the four points found in 0.2.0 are fixed (NIST Filip standard errors, the misleading "converged" on
BoxBOD's far start, `arctan` in formulas); five points are worth fixing in simplydata (Longley standard errors from
`explain`, two near-exact polynomial fits, the name of the attributable-fraction interval, the power for two percents,
and R fits called with their default stopping rules); the rest are limits shared with the package underneath or
references that are themselves less accurate. Every published value in notebooks 08 to 17 passes.

**All notebooks, 01 to 32, fifth pass (GitHub Actions run 45443130321, 2026-11-24, simplydata main 3f0192b): 26,718
comparisons, 26,264 agree, 0 simplydata error – being fixed, 64 explained difference – reference, 22 explained
difference – machine, 366 explained difference – precision limit, 0 difference not yet explained, 0 waiting for word, 2
skipped.** The same as the fourth pass but for 4 B25 rows (notebook 26) that move from agrees to machine on this runner.

**All notebooks, 01 to 32, fourth pass (GitHub Actions run 45184089331, 2026-11-23, simplydata main a5a8ff0): 26,718
comparisons, 26,268 agree, 0 simplydata error – being fixed, 64 explained difference – reference, 18 explained
difference – machine, 366 explained difference – precision limit, 0 difference not yet explained, 0 waiting for word, 2
skipped.**

**All notebooks, 01 to 32 (GitHub Actions, 2026-11-23, simplydata main a5a8ff0): 26,718 comparisons, 26,265 agree, 29
simplydata error – being fixed, 48 explained difference – reference, 27 explained difference – machine, 347 explained
difference – precision limit, 0 difference not yet explained, 0 waiting for word, 2 skipped.**

**All notebooks, 01 to 32 (GitHub Actions, 2026-11-22, simplydata main 3fd3b85): 26,478 comparisons, 26,067 agree, 110 simplydata error – being fixed, 48 explained difference – reference, 23 explained difference – machine, 223 explained difference – precision limit, 5 waiting for word, 2 skipped; no difference not yet explained** (before, on simplydata 47594a0: 23,438 comparisons, 22,727 agree, 153 waiting for word; the machine differences move between runner machines). `00_summary` reads every results file with simplydata words only and counts them.

**Notebooks 01 to 08 in the simplydata-only layout (simplydata main 47594a0; 2026-11-22, GitHub Actions): 1,785 comparisons, 1,650 agree, 28 simplydata error – being fixed (A1, A4, A7), 5 explained difference – reference (C2, C5), 10 explained difference – machine and 84 explained difference – precision limit (B1 to B5; B1 moves between runner machines), 7 waiting for word, 1 skipped** (before the new layout: 2,987 comparisons, 239 of them outside the tolerance).
Like notebooks 18 to 24 they now contain only simplydata, with their references in `references/<notebook>.csv`: printed
and certified values first (NIST StRD and the e-Handbook tables, Student, Anscombe, Bates & Watts, Dobson & Barnett,
Viechtbauer, Hackshaw et al., and the NMES fits printed by Zeileis, Kleiber & Jackman 2008), and only where nothing is
printed one program run once (mostly R). A number that a printed or certified value covers is no longer also compared
with the parent package and R, which is why there are fewer comparisons; every difference is one already explained, and
the dropped comparisons are listed in [FINDINGS.md](FINDINGS.md). The 7 waiting for word need the median survival time
as a table and the note of `fit_formula`; the 1 skipped is NIST NumAcc1's autocorrelation (3 values, 4 needed).

**Notebooks 18 to 24 (simplydata main, commit [47594a0](https://github.com/bioforkgithub/simplydata/commit/47594a0153ee1f8b026dd5c94f6dfc9b8bdcbc74); 2026-11-21, GitHub Actions, the same on Windows 11): 2366 comparisons, 2181 agree, 69 outside the tolerance (each explained),
116 skipped (115 of them waiting for a word).**
These notebooks contain only simplydata: the reference values are computed by `scripts/references/<notebook>.py`
and stored in `references/<notebook>.csv`, each with its source, its kind (printed, certified, or a run of the
reference package), a link to the manual page and the package versions; the notebook reads them with `sdata.table`
and prints the source of each check. The SKIPs are comparisons whose numbers simplydata does not give back yet (the
numbers a chart draws, and the p-value that `find_correlation_between` only prints for two columns): 115 checks marked
"waiting for word: ..." in their reason, with the simplydata word they wait for, to be switched on when the word exists.
They check, on real data, the words that the coverage audit of 2026-11-21 found checked only on made-up data or only
against the package simplydata calls. Each word is run on the data of its parent package's own manual example where
one exists (R help pages and SciPy docstrings; the reference scripts read them from the installed packages and write
them to `references/data/`, which is not stored in this repository) and on real datasets shipped with R, against independent references (R functions
simplydata does not call for that word, formulas written out, published or certified values). Eleven new points are
worth fixing in simplydata (A8 to A18 in [FINDINGS.md](FINDINGS.md); among them Dunn's test attaching p-values to
the wrong pairs with more than three groups, Barnard's and Boschloo's tests conditioning on the wrong margins in
`compare_counts_of`, and standardised residuals computed on a table with 0.5 added when one cell is empty).
Each results file of these notebooks has a column "ran on" (GitHub runner or local, and the CPU).

**Notebooks 09 to 17, now also simplydata-only (simplydata main 47594a0; 2026-11-22, GitHub Actions): 2705 comparisons, 2636 agree, 33 simplydata error – being fixed (A5, A6, A7), 16 explained difference – reference (C7, C8, C9), 7 explained difference – precision limit (B6), 1 explained difference – machine (B7), 12 waiting for word** (before the new layout: 2705 comparisons, 56 of them outside the tolerance; B7 is a Monte Carlo comparison that lands just outside 4 Monte Carlo standard errors on some Ubuntu runners).
Their reference values are in `references/<notebook>.csv` (98 printed in manuals, papers and books, 2607 from one run
of the original library), each with its source and a link. The 12 waiting for word (notebook 11) need:
`add_nmds(show="numbers")`, `add_correspondence_axes(show="numbers")` and `chart.get_drawn_numbers()`. See FINDINGS.md.
**Notebooks 25 to 32 (simplydata main, commit [47594a0](https://github.com/bioforkgithub/simplydata/commit/47594a0153ee1f8b026dd5c94f6dfc9b8bdcbc74); 2026-11-22, GitHub Actions): 16,582 comparisons, 16,260 agree, 148 simplydata error – being fixed, 18 explained difference – reference, 3 explained difference – machine, 133 explained difference – precision limit, 19 waiting for word, 1 skipped.**
They check, on real data, the groups Models, Fitting, Patterns across columns, Forecasts, Smoothing and trend,
Survival, Combining studies and Planning, in the same layout as notebooks 18 to 24 (only simplydata in the notebook;
reference values from `scripts/references/<notebook>.py`). The printed values come first (papers, vignettes, help pages
and docstrings, compared to half a unit of the last printed digit) and NIST certified values next; only where nothing
is printed does the reference script run the help-page or vignette example of the function simplydata calls, once.
Five points are worth fixing in simplydata (A20 to A24 in FINDINGS.md): the p-values of `explain(table_of_tests="in
order")`, the accuracy of `fit_formula`'s "trust region" and "dogbox" methods, the MASE scale of
`compare_forecasts_for`, the `of_degree=` option of Savitzky-Golay smoothing, and one-column CSV files (the same as
A18). On Windows 11 the same notebooks give 16,245 agree and 317 outside the tolerance, each explained (the NIST problems
near the digit line land on either side of it depending on the CPU).

## How each comparison is judged

- **Certified values (NIST StRD)**: by the log relative error (LRE), the number of significant digits that
  agree. Linear problems, ANOVA and summary statistics need at least 9 digits; nonlinear problems at least
  6, the levels used in simplydata's own reliability study (after McCullough 1999).
- **Printed values** (Student 1908, Anscombe 1973, Bates & Watts 1988, Dobson & Barnett 2008, the
  e-Handbook tables, and in notebooks 08 to 17 Viechtbauer 2010, Cohen 1988, Rosseel 2012, Coles 2001, Lumley 2010,
  Schütz et al. 2014, the BayesFactor manual and the R help pages): within half a unit of the last printed digit,
  or as stated in the notebook.
- **Parent packages and R**: a relative difference of 1e-12 (simple summaries), 1e-10 (tests from closed
  formulas), 1e-9 (least squares) or 1e-6 (models fitted by iteration). Results that come from random draws
  (Bayesian posteriors, mediation simulations) are compared exactly where the same seed reproduces the draws, and
  otherwise within 4 Monte Carlo standard errors computed in the notebook. The tolerance is set per notebook,
  is the same for every dataset in it, and is printed beside each comparison.

The tolerances are not loosened when a comparison falls outside them; the difference is reported and explained.

**CI fails on new or unexplained differences; explained ones are listed in FINDINGS.md.** Every difference stays
visible in the notebooks and in `00_summary`, with what it is (simplydata error, reference, machine, precision limit),
its FINDINGS.md section (A1 ... C25) and a one-line explanation with a link. `results/known_fails.csv` lists each explained difference (machine
status FAIL) with its finding id, the reader's words, the explanation, the link and the system it occurs on (`any`, `windows` or `ubuntu`; some nonlinear results differ between
systems). The check fails when a notebook stops with an error or when a difference is not in that list (it is shown as
"difference not yet explained"). When a listed difference of kind A (worth fixing in simplydata) agrees again, it reports "fixed - remove from
known_fails" as a warning. To update the list after changing FINDINGS.md:
`python scripts/make_known_fails.py <windows all_comparisons.csv> <ubuntu all_comparisons.csv> [...]`; after
changing only FINDINGS.md's wording, `python scripts/make_known_fails.py --labels`, then
`python scripts/label_results.py` refreshes the reader's columns of the stored results.

GitHub runs the workflow on different machines, and results that lie close to a tolerance (NIST nonlinear
and high-degree polynomial problems) can land on either side of it depending on the machine. The workflow
sets `OPENBLAS_CORETYPE=Haswell` and one thread (`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`,
`MKL_NUM_THREADS` = 1), which removed one source of this, and its "Machine" step prints the CPU and the
OpenBLAS kernels in use. Two kinds of runner still give slightly different results: for commit 37afb81, AMD EPYC
7763 runners gave 5383 PASS / 298 FAIL, and Intel Xeon Platinum 8573C and AMD EPYC 9V74 runners 5388 PASS /
293 FAIL, out of 5695.
`known_fails.csv` holds the FAILs of Windows and of both kinds; the column "same on every runner" is "no" for the
19 (all NIST nonlinear, B1) that appear on one kind only. Such a row is never reported as fixed, and a FAIL that appears on no system listed still
fails the check.

## Open in VS Code (recommended)

1. Install [VS Code](https://code.visualstudio.com/) with the **Python** and **Jupyter** extensions
   (both from Microsoft).
2. Install the packages (below) into a virtual environment or conda environment.
3. Open this folder in VS Code, open a notebook in `notebooks/`, click **Select Kernel** (top right) and
   choose the environment you installed into. Then **Run All**.

The notebooks import `notebooks/checks.py`, so open them from the `notebooks` folder (VS Code and
JupyterLab do this by default).

## Open in JupyterLab

```
jupyter lab notebooks/
```

## Install

With pip (Python 3.11 to 3.13):

```
python -m venv .venv
# Windows: .venv\Scripts\activate    macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
```

With conda (also installs R and the R packages the notebooks use, except four that are on CRAN only):

```
conda env create -f environment.yml
conda activate simplydata_tests
```

`requirements.txt` pins simplydata to commit 3fd3b85 of its main branch (installed from GitHub, as it is not yet
released on PyPI; notebooks 18 to 24 were run with it, notebooks 01 to 17 with commit 37afb81) and the versions the
saved results were made with.

### R

The R comparisons need R 4.x (https://cran.r-project.org/) with the packages listed in the workflow
(`.github/workflows/notebooks.yml`, step "Install R packages": simplydata's own R packages plus mediation and
spData, and for notebooks 18 to 24 e1071, WRS2, Hmisc, fitdistrplus, onewaytests, crank, goftest, kSamples,
randtests, perm, fBasics, multcomp and VGAM; rstanarm, rstan and posterior; qvalue from Bioconductor); rpy2 must be able to find R (set `R_HOME` if it
cannot; on Windows, run from a normal command prompt or PowerShell rather than Git Bash). Without R, the
R comparisons are skipped and the rest still run. Set `SIMPLYDATA_TESTS_REQUIRE_R=1` to make a missing R
an error, as the automatic checks do.

## Run everything

```
python scripts/get_data.py             # write the datasets not stored here to data/fetched/, check every checksum
python scripts/run_notebooks.py        # run all notebooks, save them with outputs; exit 1 on a notebook
                                       # error or a FAIL not listed in results/known_fails.csv
```

`run_notebooks.py` runs `get_data.py` itself first, so the first line is only needed to run a single notebook
in Jupyter or VS Code.

The GitHub Actions workflow (`.github/workflows/notebooks.yml`) does the same on Ubuntu with Python 3.12
and R on every push and weekly.

## Data

The datasets that are free to share are stored in `data/`: the NIST Statistical Reference Datasets (US Government
works, public domain), data published long enough ago to be in the public domain or plain published facts (Student's
sleep data, iris, Anscombe's quartet, warpbreaks, mtcars, stackloss, Puromycin, AirPassengers, Bliss 1935, and the
small tables we typed from Hollander & Wolfe, Fleiss, Willeberg and Selvin), the CC0 columbus data of spData and the
poisons data of boot (licence "Unlimited"). [data/SOURCES.md](data/SOURCES.md) gives each one's origin, citation,
licence, URL and sha256 checksum, and says for the copies taken from R's `datasets` package why they may be stored.

### Data we tested on but do not host

The other 40 datasets are not stored in this repository, because their copies are under GPL or LGPL, under
copyright, or carry no licence: they come from the R packages survival, metadat, KMsurv, vegan, lavaan, circular,
sp, ismev, qvalue (Bioconductor), drc, smatr, mcr, MASS, AER, survey, MatchIt, mediation, lme4 and BayesFactor, and
from the supplementary files of Schütz, Labes & Fuglsang (2014, The AAPS Journal 16, 1292-1297).
`scripts/get_data.py` writes them to `data/fetched/` (ignored by git) when the tests run:

- from the installed R package that ships each one (with Rscript, on your own computer; nothing is downloaded);
- the bioequivalence datasets of Schütz et al. (2014) from the publisher's supplementary files.

Each file is written the same way on every system and checked against its sha256 checksum in `data/SHA256SUMS`; if
a package or the publisher changes the data, `get_data.py` stops and says which file and source changed, so the
tests never run on other values without notice. To run it locally:

```
python scripts/get_data.py             # needs R with the packages above (see "Install"); writes data/fetched/
python scripts/get_data.py --refresh   # write every fetched file again (and download the NIST files again)
python scripts/get_data.py --check     # only check the checksums (no R, no network)
```

If R is not on your PATH, set `R_HOME`. The script is part of these tests only: simplydata itself has no network code.

## Versions in the saved results

simplydata main at commit 37afb81 (its `__version__` still reads 0.2.0), Python 3.12.10, NumPy 2.5.5, pandas 2.3.3,
SciPy 1.18.1, statsmodels 0.15.0, lifelines 0.30.3, rpy2 3.6.9, R 4.6.2 (survival 3.8-12, metafor 5.2-1,
vegan 2.7-6, lavaan 0.7-3, rstanarm 2.32.2 and the others each notebook prints). Each notebook prints its versions in its first
output.

## Citation

If you use these notebooks or their results, please cite them (see [CITATION.cff](CITATION.cff)) and cite
simplydata itself:

Victor, M. P. (2026). simplydata: simple words for data science in Python (version 0.2.0).
https://doi.org/10.5281/zenodo.24210045

Cite each dataset as given in [data/SOURCES.md](data/SOURCES.md).

## Licence

MIT (see [LICENSE](LICENSE)). The datasets keep their own licences (see `data/SOURCES.md`).
