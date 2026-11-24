# Findings

Every comparison marked FAIL in the notebooks, what the numbers are, and what we think causes it. No tolerance was loosened to make a comparison pass. Where a reference program was run to tighter convergence than its default (statsmodels, R `glm`, lifelines, R `coxph`, R `nls`, metafor, drc, extRemes), the notebook says so; that only makes the comparison stricter. The full list of every FAIL, with all numbers, is at the end of [00_summary.ipynb](notebooks/00_summary.ipynb) and in `results/all_comparisons.csv`.

simplydata main at commit 37afb81 (after 0.2.0; `simplydata.__version__` still reads 0.2.0), Python 3.12, NumPy 2.5.5, pandas 2.3.3, SciPy 1.18.1, statsmodels 0.15.0, lifelines 0.30.3, R 4.6.2 (each notebook prints its package versions), run on 2026-11-21.

**Totals on Windows 11: 5692 comparisons, 5383 PASS, 295 FAIL, 14 SKIP** (simplydata 0.2.0, notebooks 01 to 07 only: 2052, 1842 PASS, 194 FAIL, 16 SKIP). On Ubuntu (GitHub Actions): 5695 comparisons, 5383 PASS, 298 FAIL, 14 SKIP on AMD EPYC runners, and 5388 PASS, 293 FAIL on Intel Xeon 8573C and AMD EPYC 9V74 runners; the FAILs that occur on one system only are NIST nonlinear problems near the 6-digit line (B1), one R `lm` comparison (C4) and one Monte Carlo comparison (B7). `results/known_fails.csv` lists the FAILs of both systems.

Notebooks 01 to 07 test the methods of 0.2.0, notebooks 08 to 17 the methods added after it, and notebooks 18 to 24 (run on simplydata main 47594a0; their own section is at the end, before the SKIPs) the words that the coverage audit of 2026-11-21 found checked only on made-up data or only against the package simplydata calls. The FAILs are of three kinds:

- **A. Worth fixing in simplydata** (simplydata could do better than it does now): A1 (Longley only), A4, A5, A6, A7; from notebooks 18 to 24, A8 to A17. (Fourth pass, at the end: A5, A8, A17 and the last A21 rows turned out not to be simplydata errors and are now C31, C32, C33 and B1.)
- **B. Limits of the package underneath, of double-precision numbers or of Monte Carlo error**, shared by simplydata and its parent: B1 to B7.
- **C. The reference program, not simplydata**, is the one that differs: C1 to C9; from notebooks 18 to 24, C10 and C11.

From notebooks 25 to 32: A20 to A24, B20 to B25, C20 to C25 (their own section is near the end). From the second pass
(simplydata 3fd3b85, the section at the end): A25, A26, C26 to C28; A9 to A16, A20, A22 and A23 are fixed there.

**What the notebooks print for each kind.** The notebooks, `00_summary` and the column "result" of the results files
use words a reader can take at face value, and keep the machine column "status" (PASS, FAIL, SKIP) and the finding
id beside them (`results/known_fails.csv` has the words, a one-line explanation and a link to the section here):

| finding | words shown |
|---|---|
| A | simplydata error – being fixed (becomes "agrees" when the fix is released) |
| C | explained difference – reference (the reference rounds, stops early, uses another definition or has its own bug) |
| B, on every system | explained difference – precision limit (double precision, or the package simplydata calls) |
| B, on one system or one kind of runner only, and B7 | explained difference – machine (the last digits differ between CPUs or systems) |
| not in known_fails.csv | difference not yet explained (CI fails) |

## Fixed in simplydata after 0.2.0

These three were section A in the results for 0.2.0. Rerun on commit 37afb81:

| finding (0.2.0) | 0.2.0 | main 37afb81 |
|---|---|---|
| A1, NIST Filip with `fit_curve(of_degree=10)` | estimates 7.84 to 7.86 digits; standard errors wrong (intercept 45.68 for 298.08, x^4 4.277 for 71.65) or NaN (x^5 to x^10) | all 11 estimates 13.35 to 13.57 digits, all 11 standard errors 13.59 to 13.65 digits (intercept 298.0845, x^4 71.6479); the 22 FAILs against NIST now pass |
| A2, NIST BoxBOD from NIST's far start 1 | ended at b1 = 172.5, b2 = 110.9 with R squared 0.00 and still said `converged = "yes"` | still ends there (as SciPy does; R `nls` stops with an error), but now says `converged = "no"` and prints that the search may have stopped in a flat region, naming other start values or `search="whole range"`; notebook 06 checks both (2 comparisons pass). The wrong estimates are now B5 |
| A3, `arctan` in formulas (NIST Roszman1) | not in the formula reader: 2 SKIPs | Roszman1 fitted from both starts: estimates 6.22 to 7.32 digits, standard errors 6.40 to 7.27, residual SD 11.3; equal to SciPy; R `nls` within 1.2e-7 to 6.5e-7 (34 comparisons pass) |

How A1 was fixed: `fit_curve_core` now maps x onto [-1, 1] (NumPy's polynomial domain), fits statsmodels WLS there, and turns the estimates and covariance back into those of the raw powers with statsmodels `t_test`. This also lifted NIST Wampler4 and the Wampler5 x^4 term above 9 digits (B4). On near-exact fits it costs a little (A4).

## A. Worth fixing in simplydata

### A1. `table.explain` standard errors on NIST Longley still agree to 7.97 to 8.79 digits (needs 9)

The Filip part of A1 is fixed (above); the Longley part is not. Estimates pass (10.9 to 13.2 digits); standard errors:

| term | simplydata standard error | NIST certified | LRE |
|---|---|---|---|
| intercept | 890420.393 | 890420.383607373 | 7.97 |
| x2 | 0.0334910 | 0.0334910077722432 | 8.16 |
| x3 | 0.488400 | 0.488399681651699 | 8.13 |
| x4 | 0.214274 | 0.214274163161675 | 8.48 |
| x5 | 0.226073 | 0.226073200069370 | 8.79 |
| x6 | 455.478504 | 455.478499142212 | 7.97 |

Cause (unchanged from 0.2.0): `explain` fits statsmodels with `method="qr"`, and statsmodels builds the covariance as `inv(R'R)`, which squares the condition number. R `lm` and statsmodels' default (pinv) agree with NIST to about 13 digits. Suggestion: take the covariance from `inv(R) inv(R)'`, or from the pinv fit, or centre and scale the columns as `fit_curve` now does.

### A4. `table.fit_curve` on near-exact fits: Wampler1 standard errors and the Wampler3 intercept just under 9 digits

The mapping back to raw powers (the A1 fix) amplifies rounding where the fit is (almost) exact:

| quantity | simplydata (main) | simplydata 0.2.0 | NIST certified | LRE main / 0.2.0 |
|---|---|---|---|---|
| Wampler1 intercept standard error | 1.557e-09 | 1.913e-10 | 0 | 8.81 / 9.72 |
| Wampler1 x standard error | 1.710e-09 | 2.101e-10 | 0 | 8.77 / 9.68 |
| Wampler3 intercept estimate | 1.0000000015 | 1.0000000007 | 1 | 8.82 / 9.15 |

**Resolution (simplydata PR #37, 2026-11-22): explained difference, precision limit.** No existing least-squares routine was found that brings Wampler1, 3 and 4 past 9 digits without losing digits on Filip or Wampler5 (tested on the raw powers: LAPACK gelsd, gelsy and gelss through SciPy lstsq, NumPy's QR as R `lm` fits them; on the mapped powers: NumPy `Polynomial.fit`). An own iterative-refinement step reached 9.1 to 10.3 digits but was not kept, because simplydata does not write its own numerical methods. So `fit_curve` is unchanged, and simplydata's reliability study lists these numbers as explained differences with the reason: turning the mapped fit back into raw powers adds rounding of about 1e-9. Wampler5 stays at 8.0 digits, as before.

Small, but a step back from 0.2.0 on these two problems. Wampler1's points lie exactly on the curve and NIST certifies standard errors of 0, so the LRE measures the rounding noise. Suggestion: use the mapped domain only when the raw design is badly conditioned (for example a condition number above 1e10), or follow the mapped fit with one step of iterative refinement on the raw powers: residuals from the raw-power estimates, with the correction solved on the mapped design.

### A5. `compare_risks(show="attributable fractions")`: the intervals are not the "Wald intervals" the printout names

On Willeberg's dry cat food and FLUTD table (notebook 16), simplydata passes on epiR's numbers correctly: the 6 comparisons with R `epi.2by2` and the 4 printed values of its help page pass. But the printout says "Wald intervals from the risk ratio (Rothman, Greenland & Lash 2008)":

| quantity | simplydata (= epiR) | Wald, textbook formula | relative difference |
|---|---|---|---|
| attributable fraction in the exposed, 95 % CI low | 0.329129 | 0.301058 ((RR - 1)/RR on the Wald RR interval 1.4307 to 11.225) | 0.093 |
| same, high | 0.907214 | 0.910915 | 0.0041 |
| population attributable fraction, 95 % CI low | 0.323759 | 0.036138 (Jewell's Wald interval on log(1 - PAF)) | 7.96 |
| same, high | 0.749138 | 0.782377 | 0.042 |

Cause: with a single 2 x 2 table, epiR fills `AFRisk.strata.wald` from the crude **score** interval of the risk ratio (1.4906 to 10.777; (RR - 1)/RR gives exactly 0.3291289 and 0.9072135). The population fraction's interval comes from the extreme combinations of the two risks' intervals. Suggestion: either say so in the printout ("score interval of the risk ratio, carried through; epiR, one table") and name the PAF method, or take the Wald interval from `RR.strata.wald` and Jewell's formula for the PAF.

**Resolution (fourth pass): not a simplydata error, now C31.** simplydata PR #37 (commit f9e2150, in a5a8ff0) changed the
printout: it no longer says "Wald" but names the intervals epiR gives for one table (the score interval of the risk ratio
among the exposed, the exact intervals of the two risks for the population fraction). The numbers are epiR's (the 6
comparisons with R `epi.2by2` agree to 1e-16); the 4 rows against the textbook Wald formula compare another
definition and are listed as C31.

### A6. Power for two percents counts both tails, while the word says it is R `power.prop.test`

`find_chance_of_finding(percents=...)` and the "power with that many" of `find_how_many_needed(percents=...)`, in notebook 10:

| case | simplydata | R `power.prop.test` | relative difference |
|---|---|---|---|
| 50 per group, 50 % vs 75 % (R help page example, printed 0.7401659) | 0.7401672 | 0.7401659 | 1.7e-6 |
| power with that many, 50 % vs 75 % | | | 5.0e-8 |
| power with that many, 30 % vs 15 % | | | 8.4e-7 |
| power with that many, 10 % vs 5 % | | | 1.1e-6 |

Cause: the power comes from statsmodels `power_proportions_2indep(alternative="two-sided")`, which adds the far tail. That is R's `strict = TRUE` (checked: 0.7401672). R's default, and simplydata's own docstring, count only the direction of the difference. The number of people needed is not affected, because it comes from R `power.prop.test` itself. Suggestion: compute the power one-sided at half the level, as simplydata already does for the t-tests, or take `power.prop.test(n = ...)$power` from R.

### A7. R fits are called with their default stopping rules and stop 1e-6 to 1e-4 short of the optimum

simplydata gives exactly what the R function gives with its default settings. Run tighter, the same R function, and an independent Python solution, agree on a different answer. So simplydata could do better by passing a tighter stopping rule:

| notebook, word | R function and its default | simplydata (= R default) | tight reference (R, and Python) | relative difference |
|---|---|---|---|---|
| 08, `combine_studies` Paule-Mandel, BCG trials | `metafor::rma(method = "PM")`, tol 1e-5 | tau^2 0.3180937 | 0.3180685 (metafor with tol 1e-12; statsmodels `combine_effects` iterated to 1e-12; the PM equation solved with brentq) | tau^2 7.9e-5, standard error 3.2e-5, prediction interval up to 9.6e-5, p-value 4.8e-4 (23 FAILs) |
| 08, `combine_studies` REML, BCG trials | `metafor::rma(method = "REML")`, tol 1e-5 | p-value | metafor with tol 1e-12 | 1.4e-6 (1 FAIL; estimates, tau^2 and every published value pass) |
| 14, `fit_dose_response`, ryegrass | `drc::drm`, relTol 1e-7 | b 2.982219, ED10 1.463706 | drc with relTol 1e-12: b 2.982232, ED10 1.463712; SciPy `curve_fit` agrees with it | 1e-6 to 1.4e-5 (13 FAILs) |
| 13, `find_return_levels_of(values_above=30)`, Coles's daily rain | `extRemes::fevd(type = "GP")`, optim's default reltol | 100-year level 106.342311 | fevd with reltol 1e-14: 106.342764 | 1e-6 to 9e-6 (11 FAILs) |

All of these lie far below the printed precision (every published value passes), but they are avoidable. Suggestion: pass `control = list(tol = 1e-10, threshold = 1e-10, maxiter = 1000)` to `metafor::rma`, `control = drmc(relTol = 1e-12)` to `drc::drm`, and `optim.args = list(control = list(reltol = 1e-12))` to `extRemes::fevd`.

## B. Limits shared with the package underneath, or of double-precision numbers

### B1. NIST nonlinear problems: some estimates and standard errors agree to 3 to 6 digits

`fit_formula` is SciPy `curve_fit(method="lm")` (MINPACK, with the Jacobian by finite differences). On every problem simplydata gives exactly SciPy's numbers, and every SciPy comparison passes, Roszman1 included. Against NIST (which needs 6 digits, as in simplydata's own reliability study) most problems pass fully from both starts. The rest, unchanged from 0.2.0 (88 FAILs on Windows):

| problem (NIST difficulty) | start | lowest LRE: estimates | lowest LRE: standard errors |
|---|---|---|---|
| Lanczos1 (average) | 1, 2 | passes | 3.0 (and residual SD 3.0) |
| Lanczos2 (average) | 1, 2 | passes | 4.0 |
| Lanczos3 (lower) | 1, 2 | 4.9 / 4.9 | 5.0 / 4.3 |
| Bennett5 (higher) | 1, 2 | 5.6 / 5.4 | 5.2 / 5.1 |
| ENSO (average) | 1, 2 | 5.3 / 5.1 (b8) | passes |
| Nelson (average) | 1, 2 | 5.8 / 5.5 (b2) | 5.8 / 5.0 |
| Thurber (higher) | 1, 2 | passes | 5.8 / 5.7 |
| MGH10 (higher) | 1, 2 | passes | 5.9 / 5.3 |
| MGH17 (average) | 2 | passes | 5.5 |
| MGH09 (higher) | 1 | passes | 6.0 (5.97) |
| Misra1b (lower) | 2 | passes | 6.0 (5.97) |

These results lie close to the 6-digit line and depend on the floating-point library. On Ubuntu (GitHub Actions) the same notebook gives LREs up to about 1 digit higher or lower, so the set of FAILs in notebook 06 differs between the two systems; `results/known_fails.csv` lists both. Lanczos1 has a residual standard deviation of 9e-14 (an almost exact fit), so finite-difference derivatives cannot give its standard errors to more than 3 digits. A tighter stopping rule or exact derivatives (the formula reader could differentiate the formula) would likely gain digits.

### B2. One-way ANOVA on NIST SmLs07, SmLs08, SmLs09: F agrees to 4.2 to 4.4 digits

| dataset | simplydata F | NIST certified | LRE | SciPy `f_oneway` | R `anova(lm())` |
|---|---|---|---|---|---|
| SmLs07 | 21.000811887818774 | 21 | 4.41 | same as simplydata | 21.000510702936396 |
| SmLs08 | 201.01300409594847 | 201 | 4.19 | same as simplydata | 200.60120778530217 |
| SmLs09 | 2001.134926220951 | 2001 | 4.17 | same as simplydata | 687.06013841679112 |

These datasets are values near 1e12 that differ by 0.1, which cannot be stored exactly as doubles. Shifting the data before the test does not help (tested), and R is further off. This is the limit of reading these numbers into double precision.

### B3. Standard deviation of NIST NumAcc4 agrees to 8.25 digits (needs 9)

simplydata gives 0.10000000055879354; certified 0.1. pandas `std()` and R `sd()` give exactly the same number. NumAcc4 holds 1001 values near 1e7 that differ by 0.1 or 0.2, which double precision cannot store exactly.

### B4. Polynomials and a line through 0 at 7.9 to 8.9 digits

- NIST Wampler5 (the hardest of the five): the `fit_curve` estimates of the intercept to x^3 agree to 7.86 to 8.90 digits (0.2.0: intercept to x^4 at 5.68 to 8.28). Wampler4 now passes (0.2.0: 7.67 to 8.13). Here simplydata is closer to NIST than statsmodels' QR on the raw powers (C6) and R `lm` (C4).
- NIST NoInt1, fitted with `fit_formula("y = b1 * x")` because `fit_line` has no line through 0: the slope passes, and its standard error agrees to 8.70 digits (finite-difference Jacobian; R `lm(y ~ 0 + x)` matches NIST). A through-zero option for `fit_line` would be exact.

### B5. NIST BoxBOD from the far start: the estimates miss, as in SciPy and R

From NIST's start 1, `fit_formula` ends at b1 = 172.5, b2 = 110.9 (certified 213.81 and 0.5472), with residual SD 49.4 (certified 17.09) and infinite standard errors: 5 FAILs against NIST. SciPy `curve_fit` ends at the same point and R `nls` stops with an error. Since the A2 fix, simplydata says `converged = "no"` and tells the user to try other start values. From start 2 everything passes.

### B6. Standard errors from finite-difference second derivatives (GEV fits, the Box-Cox interval)

- Notebook 13, GEV fits to Port Pirie and Fremantle: the standard errors agree with R `extRemes::fevd` to 4e-6 to 1.3e-4 relative (5 FAILs; Fremantle scale 0.0114983 against 0.0114969). Both programs take them from a Hessian found by finite differences: simplydata through SciPy, extRemes through `optim`. An exact reference with `numDeriv::hessian` moved in the fourth digit with its step size, so none is used. The estimates, and the values Coles (2001) printed (0.028, 0.020, 0.098), pass.
- Notebook 14, Box-Cox power on the poisons data: simplydata's power and interval equal R `car::powerTransform` exactly. car takes the interval's standard error from `optim`'s finite-difference Hessian, while the Python reference uses a Richardson-extrapolated difference; the interval ends differ by 6e-6 and 2e-5 (2 FAILs). The printed power -0.75 and the suggested power -1 pass.

### B7. One Monte Carlo comparison of the Bayesian mixed model lands just outside 4 Monte Carlo standard errors (Ubuntu only)

Notebook 17 fits `explain("Reaction", using=["Days"], grouped_by="Subject", bayesian=True)` on sleepstudy (rstanarm `stan_glmer`, 4 chains of 2000 iterations, seed as simplydata sets it). simplydata reproduces rstanarm with the same seed exactly on both systems (every median, interval, R-hat and effective draw passes at 1e-9). Against a long run (4 chains of 10,000 iterations, another seed), every quantity must lie within 4 combined Monte Carlo standard errors (MCSEs from `posterior::mcse_median` and `mcse_quantile`). On Ubuntu, where the same seed gives different draws than on Windows, the intercept median is 250.176 against 251.865, a difference of 1.689 with an allowed 1.669 (about 4.05 MCSEs; 1 FAIL). The intercept has the fewest effective draws in that model (663, above simplydata's warning level of 400; R-hat 1.007), and the MCSE estimated from so few draws is itself uncertain. The same comparison passes on Windows, and all 26 other Monte Carlo comparisons pass on both systems. This is the Monte Carlo error of rstanarm's default chain length, not a wrong number; more iterations (simplydata could offer a longer-run option) would shrink it.

## C. The reference program differs, not simplydata

| # | comparison | what happens |
|---|---|---|
| C1 | NIST NumAcc3, NumAcc4: lag-1 autocorrelation against R `acf` | simplydata equals statsmodels `acf` exactly and passes NIST (11.9 and 10.7 digits). R `acf` is closer to NIST, so the two differ by 1.2e-12 and 1.9e-11 relative (allowed 1e-12). |
| C2 | NIST AtmWtAg: F and p-value against R `anova(lm())` | simplydata (= SciPy) agrees with NIST to 10.2 digits, R to 9.8; they differ by 1.5e-10 relative (allowed 1e-10). |
| C3 | NIST Wampler1, Wampler2: standard errors against R `lm` and against statsmodels' QR on the raw powers | the points lie exactly on the curve and NIST certifies standard errors of 0. Every program returns rounding noise (1e-9 to 1e-19), which differs between programs. Against NIST all pass except the two in A4. |
| C4 | NIST Filip, Wampler3, Wampler4, Wampler5 and Longley against R `lm` | R `lm` drops two Filip terms (SKIP) and is less accurate than simplydata on the other Filip terms and on Wampler3 to 5 (it differs from NIST by 1e-9 to 1e-6). |
| C5 | NMES 1988, negative binomial: 4 p-values against R `MASS::glm.nb` (1e-6 to 1.5e-5 relative) | estimates and standard errors agree with R, and every number agrees with statsmodels. R's `theta.ml` stops at a tolerance of about 1e-4, so p-values near 1e-65 move in the sixth digit. |
| C6 | NIST Filip, Wampler4 and Wampler5 against statsmodels `OLS(...).fit(method="qr")` on the raw powers (30 FAILs, new) | this was simplydata's own call in 0.2.0. Since the A1 fix, simplydata agrees with NIST to 13 digits on Filip. statsmodels' QR on the raw powers gets standard errors of 45.68 for 298.08 and NaN for x^5 to x^10, and its Wampler4/5 estimates differ from NIST by up to 2e-6. |
| C7 | notebook 10: Cohen (1988) Table 2.4.1 at d = 0.2, and the power of R `pwr.anova.test` and `pwr.f2.test` (5 FAILs) | Cohen prints 393 per group. The exact n is 393.4, which simplydata, R `power.t.test` and statsmodels all round up to 394; Cohen's tables use an approximation. pwr's power differs from simplydata by 1.3e-10 to 6.2e-10 relative (allowed 1e-10). simplydata equals statsmodels exactly, and pwr's own power at simplydata's n gives the power asked for. R's non-central F (`pnbeta`) is accurate to about 1e-9. |
| C8 | notebook 14: the AirPassengers airline model against statsmodels `SARIMAX` with the order fixed (5 FAILs) | simplydata equals R `forecast::auto.arima` exactly and R `arima` (ML) to within 4e-7. statsmodels' state-space likelihood, already set as R does it (`simple_differencing=True`), ends 4e-6 to 1.1e-4 away (log-likelihood 244.6965 against R's 244.6995). |
| C9 | notebook 14: MANOVA on iris (2 and 4 columns), the Hotelling-Lawley F against R `summary.manova` (6 FAILs) | simplydata equals statsmodels `MANOVA.mv_test()` in every number, and equals R in Pillai's trace (simplydata's main test), Wilks' lambda, Roy's root and the Hotelling-Lawley statistic itself (23.36465, 32.47732). Only the F approximation of the Hotelling-Lawley trace differs. statsmodels uses McKeon's (1974) approximation (F 850.90 on 4 and 174.17 degrees of freedom for 2 columns; 582.20 on 8 and 203.40 for 4); R uses another (846.97 on 4 and 290; 580.53 on 8 and 286). Both are documented approximations, and every p-value is below 1e-110. simplydata could say which one its table uses. |

## Notebooks 18 to 24: the words of the coverage audit (stream A), on real data, simplydata main 47594a0

Notebooks 18 to 24 check the summary, test, agreement, diagnostic, risk, distribution and chart words that the coverage audit of 2026-11-21 found checked only on made-up data, or only against the package simplydata calls. Each word is run on the data of its parent package's own manual example where one exists (R help pages and SciPy docstrings) and on real datasets.

**Layout.** The notebooks contain only simplydata: they load the data with `sdata.table`, run the word, read the numbers with simplydata's own words (`get_row`, `get_column`, `get_number_at`, through `checks.read_row`/`read_rows`/`read_numbers`) and compare. Every reference value is computed by `scripts/references/<notebook>.py` (R functions simplydata does not call for that word, SciPy and statsmodels, formulas written out, and the numbers printed in docstrings, help pages, NIST and papers) and stored in `references/<notebook>.csv`, with its source, its kind (`printed`, `certified` or `run reference`), a link to the manual page and the package versions; each check in a notebook prints its "Source:" line. The data the notebooks read from R packages are written by the same scripts to `references/data/` (not stored in this repository). Where simplydata has no word that gives a number back (the numbers a chart draws; the p-value and interval that `find_correlation_between` only prints for two columns), the reference is kept and the comparison is a SKIP that says so.

**Totals (notebooks 18 to 24): 2366 comparisons, 2181 PASS, 69 FAIL, 116 SKIP** (the same on Windows 11 and on GitHub Actions). The FAILs: A8 (10), A10 (3), A11 (16), A12 (16), A13 (13), A14 (2), A15 (4), all worth fixing in simplydata; B3 (1, NumAcc4 again); C10 (3) and C11 (1). The SKIPs: 63 chart comparisons of notebook 24 and 22 of `dist.draw_chart` (no simplydata word gives back the numbers a chart draws), 30 of notebook 18 (the p-values and intervals that `find_correlation_between` only prints for two columns, and their Monte Carlo check), and `compare_matched_columns(based_on="averages")`, which simplydata says is not here yet. A9, A16 and A17 were found in the first run of these notebooks, which read the printed text and the drawn figures; they can be checked again in the notebooks once simplydata gives those numbers back.

### A8. `find_percentile_of(..., interval=95)`: the quartiles' intervals are not DescTools' exact intervals, as the help says

simplydata's help says its exact interval (SciPy `quantile_test`) is "for the median and the quartiles the same as R DescTools::QuantileCI(method = "exact")". For the median it is (notebook 18: every month of airquality passes); for the quartiles it is not (10 FAILs of 20):

```python
t = sdata.table(airquality)   # R datasets::airquality
t.find_percentile_of("Ozone", 25, for_each_group_in="Month", interval=95)
# May: high end 16; R DescTools::QuantileCI(x, probs = 0.25, method = "exact") gives 14
# (75th percentile: May low end 19 against 23; September high end 78 against 73)
```

simplydata equals SciPy (the 90th percentile's interval equals `quantile_test` to 1e-12). Both intervals are exact (order statistics with at least 95 % coverage) but choose different pairs of order statistics away from the median. Suggestion: say "for the median" only, or name the rule.

**Resolution (fourth pass): not a simplydata error, now C32.** simplydata commit 54dbbf9 (in 3fd3b85 and a5a8ff0)
changed the help and the method basis: the interval is SciPy `quantile_test`'s (Thompson 1936), and for percents other
than 50 DescTools picks another pair of order statistics, so the ends can differ; both are exact. The 10 rows are listed
as C32.

### A9. Spearman's p-value for few rows is the t approximation; R gives the exact p-value

`find_correlation_between(test_method="spearman")` takes the p-value from SciPy (t approximation). R `cor.test` gives the exact p-value (AS 89) below 1290 rows without ties. On the example of R's `cor.test` help page (Hollander & Wolfe, 9 rows) simplydata prints p = 0.088 and R gives 0.0968 (1 FAIL; against R's `exact = FALSE` every Spearman p-value passes). Suggestion: the exact p-value for small samples, or say that it is the t approximation.

### A10. `table.compare_counts_of(test_method="barnard" / "boschloo")` fixes the margins of the wrong variable

Barnard's and Boschloo's tests are not symmetric: the two groups are the two samples. On the twins table of R's `fisher.test` help page (Fisher 1962: convicted or not, by twin type) `sdata.compare_counts(counts)` gives Barnard's p = 0.000348, as R Exact `exact.test` with the twin types as the samples and as Barnard's test written out in NumPy (a grid of 20 001 values). `table.compare_counts_of` on the same data as rows gives 0.000280, which is Exact with the outcome as the samples (3 FAILs: Barnard against both references, Boschloo against Exact):

```python
rows = ...  # 30 rows: twins (Dizygotic/Monozygotic) and convicted (yes/no), R's fisher.test example
sdata.table(rows).compare_counts_of("convicted", for_each_group_in="twins", test_method="barnard")   # 0.000280
sdata.compare_counts([[2, 15], [10, 3]], test_method="barnard")                                      # 0.000348
```

On a symmetric table (tea tasting) both agree. Suggestion: give `exact.test` the table with the `for_each_group_in=` groups as its rows.

### A11. `compare_counts_of(show="cells")` with an empty cell: the residuals of a table with 0.5 added to every cell

statsmodels `Table` adds 0.5 to every cell when one cell is 0 (`shift_zeros=True`, its default). On the job satisfaction table of R's `fisher.test` help page (Agresti 2002; one cell is 0) all 16 standardised residuals differ from R `chisq.test()$stdres` (16 FAILs; "< 15k, VeryD" 0.0802 against 0.2096), and the expected counts shown are those of the shifted table. Suggestion: `Table(counts, shift_zeros=False)`.

### A12. `compare_each_pair_of(test_method="dunn")` with more than three groups: p-values attached to the wrong pairs

On InsectSprays (six sprays) the Holm-adjusted p-values are right numbers attached to other pairs (16 FAILs: 8 against R PMCMRplus `kwAllPairsDunnTest`, 8 against Dunn's test written out in NumPy; the two references agree):

```python
sdata.table(insectsprays).compare_each_pair_of("count", for_each_group_in="spray", test_method="dunn")
# C - B: 0.586 (right: 5.1e-06; 0.586 is D - C's); E - D: 3.3e-06 (right: 1); D - A: 5.1e-06 (right: 0.0128)
```

With three groups (PlantGrowth) every pair is right. Probably FSA `dunnTest`'s order of comparisons is assumed to be the order of simplydata's pairs; matching by the comparison names would fix it.

### A13. Pairwise Mann-Whitney with ties: the normal approximation, where R (since 4.4.0) is exact

`compare_each_pair_of(based_on="ranks")` gives a pair the exact p-value only without ties, R's rule before 4.4.0. R 4.4.0 made `wilcox.test` exact with ties too (below 50 values), and simplydata's help names R `pairwise.wilcox.test`. Pairs without ties agree to 1e-15; pairs with ties differ (13 FAILs; InsectSprays C - A after Holm 5.1e-4 against 1.6e-5). `compare_averages_of(based_on="ranks")` likely uses the same rule. Suggestion: the exact conditional p-value with ties, or say which rule is used.

### A14. Krippendorff's alpha with three or more raters and no empty cells: irr's weighting, inherited

`agreement_between(test_method="krippendorff")` calls R irr `kripp.alpha`. With complete data irr weights each pair of values in a unit by 1 instead of 1/(m - 1) (its code: `mc <- rep(1, ...)` when there is no NA), right for two raters only. Krippendorff's coincidence-matrix formula written out (checked first on Krippendorff's own published example, 0.743 nominal and 0.849 interval, which irr also gets right because that example has gaps) gives 0.1701 for irr's `anxiety` (3 raters, interval; simplydata 0.1631) and 0.4334 for `diagnoses` (6 raters, nominal; simplydata 0.4309): 2 FAILs. Suggestion: compute alpha in Python (a few lines) or work round irr.

### A15. ICC(A,k), the average of raters with absolute agreement: irr's interval, inherited

`agreement_between(..., checking_for="agreement", of="average of raters")` gives irr's interval (-0.1373, 0.7463) for `anxiety`; McGraw & Wong's (1996) formula written out and R psych `ICC(lmer = FALSE)` agree with each other on (-0.1266, 0.7451) (4 FAILs). irr puts the average-rater ICC, times the number of raters, into the degrees of freedom meant for the single-rater ICC. The estimate, F, p-value and every other form's interval agree with all three.

### A16. `draw_smooth_line()` (LOWESS) is statsmodels with `delta = 0`, not R `lowess`'s default

simplydata's help says the line uses the nearest 2/3 of the dots "as R lowess()". R `lowess` by default interpolates between points closer than 1 % of the x range (`delta`); statsmodels fits every point. On mtcars (weight, mpg) the lines differ by up to 5e-4 mpg (2.9e-5 relative; 1 FAIL); with R's `delta = 0` they agree to 1e-10. Suggestion: `delta=0.01 * np.ptp(x)`, or say `delta = 0`.

### A18. `sdata.table` on a CSV file with one column guesses a letter as the separator (the same as A24 of stream B)

Found while writing the new layout: `sdata.table("darwin.csv")` on a file whose only column is `difference` gives the columns `['diff', 'r', 'nc', 'Unnamed: 3']` (it takes `e` as the separator). The notebooks pass `separator=","`. Suggestion: never guess a letter or digit as a separator; with one column, there is none.

### A17. The smoothing spline with a number for `smoothness=` differs from R `smooth.spline(lambda = )` (not because x repeats; see C33)

On the Lake Huron years (no repeated x) simplydata agrees with R to better than 1e-6; on mtcars weight (some weights repeat) the curves differ by up to 1.2e-5 relative (1 FAIL). The cause was not found; R collapses repeated x values into one point with a weight, which simplydata may do differently.

**Resolution (fourth pass): not a simplydata error, now C33; the cause above was wrong.** The repeats are not the
cause: R collapses repeated x exactly as simplydata does (average y, weight = how many; default weights are not
rescaled), and with every weight made distinct (+ 1e-4 per row) R and SciPy still differ by 1.17e-5. R's curve gives a
larger penalised sum than SciPy's, so R lands slightly off the optimum; see C33.

### A7 and B3, again

- A7: `fit_dose_response` (the chart) on ryegrass is R drc with its default stopping rule and differs from a tight SciPy fit by 7.8e-6 relative (1 FAIL), as `table.fit_dose_response` in notebook 14.
- B3: `nums.find_standard_deviation` on NIST NumAcc4 agrees to 8.25 digits, as `table.find_standard_deviation_of` (1 FAIL).

### C10 and C11 (the reference differs)

| # | comparison | what happens |
|---|---|---|
| C10 | notebook 23: `compare_distributions_for` on the airquality ozone against R fitdistrplus `gofstat` (3 FAILs); notebook 24: the gamma curve of `fit_distribution` against `dgamma` at the fitdistrplus fit (1 FAIL) | fitdistrplus fits by `optim`, which even with `reltol = 1e-14` stops about 1e-7 short of the maximum; the goodness-of-fit statistics and the curve move by 1.2e-6 to 2e-6. simplydata agrees to 1e-10 with the same quantities at the exact maximum likelihood fit (the likelihood equations solved in SciPy), and with fitdistrplus's estimates and log-likelihood at 1e-6. |
| C11 | notebook 21: Cohen's kappa p-value on irr's `diagnoses` against R irr `kappa2` (1 FAIL) | the z statistic agrees; irr computes the p-value as 2 * (1 - pnorm(z)), which loses digits at p = 2.6e-12 (relative difference 4e-5). simplydata equals the normal tail computed directly to 1e-15. |

The one SKIP of notebooks 18 to 24: `compare_matched_columns(based_on="averages")`, which simplydata says is not here yet (it points to `compare_averages_of(matched_by=)`, checked in notebook 19).

## Notebooks 09 to 17 in the simplydata-only layout (simplydata main 47594a0, 2026-11-22)

**Layout.** Notebooks 09 to 17 now hold only simplydata, as notebooks 18 to 24 do: the reference values are made by `scripts/references/<notebook>.py` and stored in `references/<notebook>.csv`; every check prints its source with a link (a version-pinned manual page, or the reference script itself for a formula written out; `scripts/references/links.py` names the links `common.url_of` cannot). The tables a notebook needs with recoding (yes/no for died, labels for groups) are written by its reference script to `references/data/`. Of the 2705 reference values, 98 are printed in a manual, paper or book (Cohen 1988, the G*Power manual, R help pages, the pwr and vegan vignettes, Rosseel 2012, Coles 2001, Box & Jenkins, Box & Cox, Storey & Tibshirani, Lumley 2010, Schütz et al. 2014, the BayesFactor manual and others) and 2607 come from running the original library once.

**Totals: 2705 comparisons before and after.** Before the new layout: 2649 within the tolerance, 56 outside it. Now (GitHub Actions run 44963165536): 2636 agrees, 33 simplydata error – being fixed (A6 in notebook 10, A7 in 13 and 14, A5 in 16), 16 explained difference – reference (C7 in 10, C8 and C9 in 14), 7 explained difference – precision limit (B6 in 13 and 14), 1 explained difference – machine (B7 in 17: a Monte Carlo comparison that lands just outside 4 Monte Carlo standard errors on some Ubuntu runners), 12 waiting for word, 0 skipped, no difference not yet explained. Every difference keeps its key and its row in `results/known_fails.csv` (make_known_fails.py gives the same 57 rows for 09 to 17 as before). No new disagreement was found.

What changed in the comparisons, and why:

- **Twelve checks of notebook 11 now wait for a simplydata word** ("waiting for word"). The old notebook read these numbers from simplydata's printout or from the drawn chart, which a simplydata-only notebook cannot do: the NMDS stress of dune and varespec (`add_nmds(show="numbers")`, 2), the eigenvalues, total inertia and percent per axis of the correspondence analysis (`add_correspondence_axes(show="numbers")`, 7, three of which compared NumPy's SVD with the printed tutorial and never involved simplydata; the printed values are kept as references), and the heights of the cluster tree (`chart.get_drawn_numbers()`, 3).
- **Notebook 10:** four checks ran R pwr at simplydata's own exact sample size, which needs a simplydata number on the reference side. They now compare simplydata's exact size with the size solved from pwr's own power function (`uniroot`, tol 1e-12); "power with that many" is worked out at the reference's own rounded-up size, and simplydata's rounded-up size is checked against it.
- **Whole-column checks renamed:** where the old notebook computed one largest difference or one count of mismatches between two tables (UDCA times in 09, Bray-Curtis distances, clr values and cluster labels in 11, q-values and local FDRs in 14, design layouts, randomisation lists and sampled rows in 15), the reference script now stores the whole column and the notebook compares it value by value; the comparison is the same, its name says "largest relative difference over N values".
- **Notebook 15** reads the survey degrees of freedom from the column `find_average_of` gives back, instead of a private attribute.
- **Notebook 16:** the limit a × b of the mediation check is now the product of simplydata's own two coefficients, against R `lm` (before, statsmodels against R only).

## Notebooks 01 to 08 in the simplydata-only layout (2026-11-22, simplydata main 47594a0)

Notebooks 01 to 08 now have the layout of notebooks 18 to 24: the cells contain only simplydata (data read with
`sdata.table`, numbers read back with `get_row` through `checks.read_row`/`read_rows`/`read_rows_by`, log and Fisher z
scales made with `add_column`), and every reference value comes from `scripts/references/<notebook>.py`, stored in
`references/<notebook>.csv` with its source, kind and link. 00_summary reads the results with simplydata words only.

Each reference is now, in this order: a printed or certified value (NIST StRD, the NIST/SEMATECH e-Handbook tables,
Student 1908, Anscombe 1973, Bates & Watts 1988, Dobson & Barnett 2008, Viechtbauer 2010, Hackshaw et al. 1997, and,
new, the NMES 1988 Poisson and negative binomial fits printed by Zeileis, Kleiber & Jackman 2008, JSS 27(8)); only
where nothing is printed, one program run once by the reference script (R where the notebook used R: base R, MASS,
survival, metafor, meta; SciPy or statsmodels only where R gives no such number), with the tolerances as before. A
quantity that a printed or certified value covers is no longer also compared with the parent package and with R.

**Totals (notebooks 01 to 08, GitHub Actions run of 2026-11-22): 1785 comparisons, 1650 agree, 28 simplydata error –
being fixed, 5 explained difference – reference, 10 explained difference – machine, 84 explained difference – precision
limit, 7 waiting for word, 1 skipped** (PASS/FAIL/SKIP 1650/130/8; before: 2987 comparisons, 2734 PASS, 239 FAIL, 14 SKIP). Every FAIL is one already explained: A1 (6),
A4 (3), A7 (19), B1 (78 to 81, by system and runner machine), B2 (5), B3 (1), B4 (5), B5 (5), C2 (1), C5 (4). No new finding.

| notebook | before | after | what was dropped, and why |
|---|---|---|---|
| 01 | 98 (94 PASS, 3 FAIL, 1 SKIP) | 37 (35, 1, 1) | the pandas, NumPy, statsmodels and R comparisons of the certified mean, standard deviation and autocorrelation (61, with the two C1 FAILs); the median against R only |
| 02 | 118 (108, 10, 0) | 66 (60, 6, 0) | the SciPy and statsmodels comparisons (41), and R's F where NIST certifies it (11, with the C2 F and two B2 FAILs) |
| 03 | 436 (329, 105, 2) | 168 (154, 14, 0) | every SciPy, statsmodels (QR) and R comparison of a NIST certified number (268, with the C3, C4 and C6 FAILs and the 2 R SKIPs); mtcars and stackloss against R only |
| 04 | 159 (155, 4, 0) | 87 (83, 4, 0) | the statsmodels and SciPy comparisons (75); R's NMES estimates and standard errors, now compared with the printed fits instead; the "theta" check, which compared statsmodels with R and not simplydata |
| 05 | 197 (197, 0, 0) | 36 (30, 0, 6) | lifelines (163; R survival is the original); the event-time comparisons as whole columns; the median time waits for a word (6) |
| 06 | 994 (890, 93, 11) | 542 (455, 86, 1) | SciPy `curve_fit` and R `nls` on the NIST problems (434 comparisons and 11 SKIPs), Puromycin and the beetles against SciPy, statsmodels and R where printed values exist; the LD50 is now read from `fit_dose_response` (the same model) instead of being worked out in the notebook; the flat-region note waits for a word |
| 07 | 84 (84, 0, 0) | 44 (44, 0, 0) | SciPy on the grid (40); R kept |
| 08 | 901 (877, 24, 0) | 805 (786, 19, 0) | statsmodels and the Python formula (51, with the 5 statsmodels A7 FAILs), R metafor with its default stopping rule for the BCG iterative fits, R where the published values cover the quantity, and R meta for Egger's test (metafor `regtest` kept) |

The findings whose comparisons are no longer made stay in this file, as they were found: C1 (R `acf`), C3 and C4
(R `lm`), C6 (statsmodels QR), the R side of C2's F, and the statsmodels side of A7. The simplydata numbers they
concern are still compared with the certified or tight reference (A7's R metafor comparisons remain).

Two checks wait for a simplydata word: `find_survival_over_time(..., show="median")`, the median time and its 95 % CI
per group as a table (6 comparisons; simplydata prints them); and a "note" column in the `fit_formula` result, so
that the note "the search may have stopped in a flat region" (FINDINGS A2) can be read (1 comparison).

## SKIPs (14, before the simplydata-only layout of notebooks 01 to 08)

- NIST NumAcc1 lag-1 autocorrelation: NumAcc1 has 3 values and `find_autocorrelation_of` asks for at least 4 (a stated limit, not a wrong number). Still a SKIP.
- NIST Filip: R `lm` gives NA for two terms (it drops them as collinear). No longer compared (NIST certifies them).
- NIST nonlinear, R `nls` (11): R stopped with an error, or did not itself reach the certified values from that start, so it could not serve as a reference there. No longer compared (NIST certifies them).
- Roszman1 is no longer skipped (A3 fixed).

## Notebooks 25 to 32: real-data checks of the groups Models, Fitting, Patterns across columns, Forecasts, Smoothing and trend, Survival, Combining studies and Planning

simplydata main at commit 47594a0 (method choice by `test_method=`), run on GitHub Actions (one job per notebook,
`.github/workflows/run-gaps-b.yml`; the column "ran on" of each results file names the runner) and on Windows 11
while the notebooks were written. The notebooks hold only simplydata: data read with `sdata.table`, the words, and the
numbers read back with simplydata's own words (`get_column_names`, `get_column`, `get_row`, `get_cell_at`, through
`notebooks/reference_checks.py`). Their reference values come from `scripts/references/<notebook>.py`, saved to
`references/<notebook>.csv` in the columns of `scripts/references/common.py`: first the numbers printed in a paper,
vignette, help page or docstring (compared to half a unit of the last printed digit) or certified by NIST; only where
nothing is printed, the help-page or vignette example of the function simplydata calls, run once (Zeileis, Kleiber &
Jackman 2008; Dobson & Barnett 2008; Ferrari & Cribari-Neto 2004; McCullagh & Nelder 1989; Bates et al. 2015; the
emmeans, quantreg, betareg, trend and pwr vignettes; Rosseel 2012; Viechtbauer 2010; Lumley 2010; Cobb 1978; NIST
StRD). Every comparison prints its "Source:" with a link. Checks that need a word simplydata does not have yet are
written as they will read and recorded as SKIP with the reason "waiting for word: ...". The findings of these
notebooks are numbered from 20, so that they cannot clash with findings added for other notebooks.

Totals (GitHub Actions, ubuntu-latest): 16,582 comparisons, 16,260 agree, 302 outside the tolerance (148 simplydata
error – being fixed, 18 explained difference – reference, 3 machine, 133 precision limit), 19 waiting for word, 1
skipped (the A23 refusal); on Windows 11: 16,245 agree, 317 outside the tolerance. Every FAIL is explained below and listed in `results/known_fails.csv`.

### A20. `explain(table_of_tests="in order")`: the p-values use the wrong denominator degrees of freedom

Moore's conformity data (the help-page example of `car::Anova`), `conformity ~ fcategory * partner.status` (45 rows,
39 residual degrees of freedom). The sums of squares and F statistics equal R `anova()` and statsmodels
`anova_lm(typ=1)` exactly, but the p-values do not:

| term | F (= R) | simplydata p | R `anova()` p (F on 39 df) | p from F on 41 df |
|---|---|---|---|---|
| fcategory | 0.089023 on 2 | 0.914996 | 0.915010 | 0.914996 |
| partner.status | 10.1207 on 1 | 0.0027934 | 0.0028742 | 0.0027934 |

The "in order" p-values are taken on 41 degrees of freedom (the residual df of the model without the interaction),
although the F statistics divide by the full model's residual mean square. Types 2 and 3 are right. Minimal
reproduction: `sdata.table(Moore).explain("conformity", using=["fcategory * partner.status"], table_of_tests="in
order")` against `anova(lm(conformity ~ fcategory * partner.status, Moore))`.

### A21. `fit_formula(test_method="trust region" / "dogbox")` reach only 2 to 5 digits where the default reaches 7

From NIST's second start, the Levenberg-Marquardt default of `fit_formula` gets every lower-difficulty problem to 7 or
more digits (Misra1a: standard errors 7.4 and 8.3 digits). With `test_method="trust region"` or `"dogbox"` the same
problems stop earlier: Misra1a's standard errors agree to 4.6 and 4.7 digits (2.70694 for 2.70701), Kirby2's estimates
to 5.0 to 5.6 and standard errors to 3.5 to 5.3 digits, Hahn1's to 1.6 to 2.7 digits (147 FAILs on Windows, 139 on
the GitHub runner; about 380 of these comparisons pass). The cause is SciPy `least_squares`' default stopping rules
(`ftol`, `xtol`, `gtol` 1e-8) and the finite-difference Jacobian behind the standard errors. Suggestion: pass
`ftol = xtol = gtol = 1e-15` (as for the default method) and take the covariance from a more accurate Jacobian at the
solution (`jac="3-point"` or the formula's derivatives). Minimal reproduction: NIST Misra1a,
`fit_formula("y = b1*(1-exp(-b2*x))", start={"b1": 250, "b2": 0.0005}, test_method="trust region")`.

**Fourth pass.** After PRs #39 and #41 the 14 rows left on the pull request's run of a5a8ff0 are all NIST Lanczos1
(standard errors and residual standard deviation, 3.3 to 3.6 digits). They are the double-precision limit that B1
already describes for the default method, and are now listed as B1 (see "Fourth pass" at the end). A21 has no rows
left.

### A22. `compare_forecasts_for` scales MASE by the last value, not by the same season last time it promises

For the seasonal series co2 and USAccDeaths, `compare_forecasts_for(..., repeats_every=12)` says MASE is scaled by
the in-sample error of the same season last time. R `forecast::accuracy(f, test)` with the held-out readings as a
seasonal `ts` does that; simplydata passes them as plain numbers, so accuracy() scales by the one-step naive error:

| series, method | simplydata MASE | R accuracy, test as ts (seasonal scale) | R accuracy, test as numbers |
|---|---|---|---|
| co2, exponential smoothing | 0.277792 | 0.234368 | 0.277792 |
| USAccDeaths, ARIMA(0,1,1)(0,1,1)12 | 0.362187 | 0.480829 | |

MAE, RMSE and MAPE are right; the ranking can change when MASE is near 1. Suggestion: pass the held-out readings as
`ts(test, frequency = m, start = ...)`. The same code serves `forecast(check_on_last=)`.

### A23. Savitzky-Golay smoothing refuses `of_degree=`, the option its docstrings name

`table.add_smoothed_column(..., test_method="savitzky-golay", of_degree=2)` and `nums.smooth(test_method="savitzky-golay",
of_degree=2)` stop with "does not take of_degree= (it takes degree=)", and `degree=` is renamed to `of_degree=` by the
words' spelling rule and refused the same way, so the degree cannot be set at all. The default degree 2 is right
(SciPy `savgol_filter` and R `signal::sgolayfilt` pass). The cause is the option table in `_smoothing.py`
(`OPTIONS["savitzky-golay"] = {"degree"}`), not updated when the option was renamed. Notebook 29 records the two
refusals as SKIPs.

### A24. A CSV file with one column is read with "s" as its separator

Same as A18 (stream A): `sdata.table()` guesses the separator of a one-column file wrongly. Notebook 28 avoids it by
writing the dune species table with a second column.

### Found with the first layout of notebooks 25 to 31 (results at commit 9c5a918)

The first layout also compared each word with tight reruns of the reference (R with `tol`/`reltol` 1e-12, formulas in
NumPy, statsmodels, lifelines). The printed-first layout no longer runs these, so the findings A7 (more cases), B20,
B21, B6 (more cases), C5 (more cases), C20, C22, C23 and C25 below have no FAIL in the present results; they stay
here as what those comparisons showed.

### A7, more cases: simplydata's R calls with default stopping rules

| notebook | word | R call | differences from the tight reference |
|---|---|---|---|
| 25 | `explain(outcome="counts", robust_errors=True)` | GLM refitted by `glm()` with its default epsilon 1e-8, then `sandwich::vcovHC(type = "HC3")` | 1.9e-6 to 2.4e-6 in 3 standard errors (against the same call on a tight fit and the HC3 formula in NumPy); equal to the default call |
| 27 | `fit_dose_response(of_kind="5 parameter log-logistic" / "weibull")` | `drc::drm`, relTol 1e-7 | 2e-6 to 1e-3 (78 FAILs; the 3-parameter and binomial curves pass) |
| 31 | `combine_studies` (REML) on Hackshaw, Curtis, McDaniel, BCG | `metafor::rma`, tol 1e-5 | tau^2 2.7e-6 to 2.9e-5, estimates and intervals 2e-6 to 1e-5, one p-value 1.3e-6 (15 FAILs against REML solved by Fisher scoring to 1e-14; every metafor comparison and every published value passes) |

### B20. p-values far in the tail magnify differences far below the tolerance

Where estimates and standard errors agree with the reference within 1e-6, a p-value near 1e-40 or 1e-146 moves by about
z^2 times the relative difference in z: the probit and complementary log-log slopes of Bliss's beetles (p 3.7e-40 and
1.7e-34; 2.2e-5 and 4.2e-6 relative), the intercept of the sleepstudy random-intercept model (p 1.0e-146; 8.4e-6) and
one hurdle-model p-value on a GitHub run. The z statistics agree to about 1e-7.

### B21. R `confint()` of a GLM interpolates its profile on a grid

`explain(..., interval_by="profile")` gives exactly R `confint()`. The exact profile-likelihood interval, solved here
with SciPy `brentq`, differs by 3.9e-6 to 4.2e-6, because R evaluates the profile at a few points and interpolates.

### B22. A chi-square p-value below the smallest normal double is 0

`compare_models(outcome="counts")` on NMES 1988 (likelihood-ratio statistic 1466.8 on 1 df): simplydata and SciPy
`chi2.sf` give 0, R gives 6.3e-321 (a subnormal number). The statistic agrees exactly.

### B23. The smoothing spline's generalised cross-validation optimum, SciPy against R

`nums.smooth(test_method="spline", smoothness="auto")` (SciPy `make_smoothing_spline`, GCV) and R `smooth.spline()`
choose slightly different penalties; the smoothed Nile flow differs by up to 9.7e-6 relative. With a fixed penalty
(`smoothness=1000`, R's lambda) every comparison passes (to 2e-8).

### B24. Oblimin rotation converges to GPArotation's tolerance

`find_hidden_factors_behind(test_method="minres", rotation="oblimin")` on Holzinger & Swineford: loadings within 4.8e-6
of `psych::fa` run directly, at GPArotation's convergence tolerance (1e-5). Uniquenesses pass.

### B6 and B1, more cases

B6 (first layout): the standard error of one ordinal-model term (housing, Type Atrium) differs by 1.1e-6 from both tight references
(numerical Hessian). B1: `fit_formula` with the default method from NIST start 2 on Lanczos1 to 3, Thurber, MGH17, ENSO,
Bennett5, MGH10 and Misra1b, as in notebook 06.

### C5, more cases

The negative binomial p-values on NMES 1988 against `MASS::glm.nb` (1.4e-6 to 2e-5); the estimates and standard errors
pass and match Zeileis et al. (2008).

### B25. Two glmer fits of the same model agree to about 1e-5

`find_model_averages_of(outcome="yes or no", grouped_by="herd")` on cbpp (lme4): the model averages and their standard
errors differ from `emmeans(glmer(...), type = "response")` run directly by up to 1.3e-5 relative (4 standard errors on every system; the 4
model averages, by 2e-6 to 1e-5, on Windows and on some GitHub runners only). simplydata calls the same `glmer`; the Laplace fit stops at lme4's convergence tolerance, and the herd
coded as text instead of numbers starts the optimiser from another point. The estimates on the link scale pass at
1e-6.

### C7, more cases

`find_chance_of_finding(for_test="explain", effect_size=0.15, using=3, values_per_group=77)` against R
`pwr.f2.test(u = 3, v = 73, f2 = 0.15)`: 1.5e-10 apart (pwr integrates the noncentral F; simplydata uses SciPy's), as
in notebook 10. In notebook 32 the vignette prints `w = 0.2886751` rounded but passes it unrounded, so the printed power
(0.6739834) is compared within half a unit plus the effect of w's rounding (4.5e-7).

### C20. statsmodels references less accurate than simplydata

statsmodels `MixedLM` (REML) standard errors 1.5e-6 to 1.5e-5 from simplydata and R lme4 (lme4 with a tight optimiser
agrees with simplydata to 1e-8); statsmodels `GEE` with an exchangeable working correlation 1e-6 to 1.2e-5, and with
its AR-1 structure (another moment estimator of the correlation) up to 7.7e-2, while simplydata equals R
`geepack::geeglm` in every number.

### C21. R and lifelines references at their default stopping rules

R `MASS::polr` with its defaults (up to 7e-5; with `reltol = 1e-12` it agrees with simplydata and statsmodels), R
`nnet::multinom` (1.3e-6 to 1.6e-6; tight: pass), R `factanal` (up to 1e-5; `psych::fa(fm = "ml")` passes), lifelines
`CoxPHFitter`, `WeibullAFTFitter` and `LogNormalAFTFitter` (1.5e-6 to 7.8e-6; R `coxph` and `survreg` run tight pass).
In the present layout: `polr` (21 FAILs, housing), `multinom` (4) and the help-page `glm(family = Gamma(link = "log"))`
on McCullagh & Nelder's clotting data (up to 2.2e-5; the printed coefficients -0.01655 and 0.01534 pass).

### C22. R `MASS::rlm` and statsmodels `RLM` converge to different points

Huber's M-estimate of the stack loss data: simplydata equals statsmodels `RLM` and the output printed in its docstring
(8 decimals); R `rlm` (acc 1e-12) differs by 1.7e-6 to 7.1e-6.

### C23. `mcr`'s Passing-Bablok line is not exactly the 1983 algorithm

simplydata equals `mcr::mcreg(method.reg = "PaBa", method.ci = "analytical")` in every number. The algorithm of Passing
& Bablok (1983) written out in NumPy gives a slope 4.3e-9 and an intercept 5.8e-8 away, and interval ends up to 3.9e-4
away (1.17308 = 61/52, a pairwise slope, where mcr gives 1.17300).

### C24. spdep's `knearneigh(longlat = TRUE)` picks other neighbours than great-circle distances on a sphere

Moran's I of the Fiji earthquake magnitudes with 8 nearest neighbours: simplydata (sf points, crs 4326) equals the
statistic written out in NumPy with haversine distances (1e-10); spdep with `longlat = TRUE` on the plain coordinates
picks other neighbours for a few places and differs by 2.9e-3.

### C25. Fine-Gray by `cmprsk::crr` and by `survival::finegray` + `coxph`

simplydata equals `cmprsk::crr`; the survival package's route (its competing-risks vignette) gives coefficients 1.1e-3
to 2.1e-3 apart on mgus2, whose times in months have many ties, which the two handle differently.

### Notes (not FAILs)

- `explain()` under `use_survey_design()` names the intercept "(Intercept)"; elsewhere it is "(intercept)".
- Two R references follow what simplydata documents rather than the help-page default: `arima(method = "ML")` (the
  help page uses "CSS-ML", which differs in the sixth digit) and the cross-correlation with the sign of the lag
  reversed.
- Gray's test (`compare_survival(competing_with=)`) has no independent implementation here beyond `cmprsk::cuminc`, the
  call simplydata makes.

## Second pass: simplydata main 3fd3b85 (PR #39), 2026-11-22

simplydata PR #39 (merge 3fd3b85) added the words that 153 checks of notebooks 01 to 32 waited for
(`chart.get_numbers_drawn()`, the table that `show_rows_that_pull_model` now gives back, `check_model_of(show="numbers")`,
`find_correlation_between(a, b, show="result table")`, `find_survival_over_time(show="median")`, the "note" column of
`fit_formula`, `add_nmds(show="numbers")`, `add_correspondence_axes(show="numbers")`; also `repeat_each_row_by` and
`sdata.get_example(from_r_package=)`) and fixes for A20 to A24. The notebooks now install simplydata at 3fd3b85
(`requirements.txt`) and the waiting checks are switched on. Every notebook was first run unchanged on 3fd3b85 (GitHub
Actions run 45069126571, `.github/workflows/run-second-pass.yml`) to see what the fixes alone change.

**What PR #39 fixed** (every comparison of the finding passes on the GitHub runner): A10 (3), A11 (16), A12 (16), A13
(13), A14 (2), A15 (4), A20 (2), A22 (7); A23 (Savitzky-Golay takes `of_degree=` again: the skipped refusal of notebook 29 is
now a check against SciPy `savgol_filter(x, 11, 2)`, 100 values, all agree); A16 (the LOWESS line of `draw_smooth_line`, read with
`chart.get_numbers_drawn()`, now agrees with R `lowess` at its default delta to 1.5e-10); and most of A21 (120 of its FAILs pass; 24 still fail, plus A26 below). These rows are taken
out of `results/known_fails.csv`; rows that depend on the runner or were seen on Windows only are kept.

**R package icr.** simplydata 3fd3b85 computes Krippendorff's alpha with R `icr::krippalpha`; the workflows and
`environment.yml` now install icr (notebook 21 stopped without it).

**Reference names with a package version.** Some reference sources name the R package version (such as "R spdep
moran.mc ... (spdep 1.4.3)"); a new version on CRAN changes the name, so `results/known_fails.csv` also lists the
explained FAILs under the new name (C24: 1.4.2 and 1.4.3).

**A9 fixed too.** `find_correlation_between(test_method="spearman")` now takes R's exact p-value (AS 89) for few rows
without ties, as its help says, and agrees with R `cor.test` on the Hollander & Wolfe example and mtcars. The two-sided
reference with `exact = FALSE` (the t approximation, another definition) is no longer compared; the one-sided rows name
R's default exact p-value as their source.

**The checks that waited for a word.** Of the 153 "waiting for word" results of 47594a0, the ones whose word now gives
the number back are switched on; their notebooks read the number with simplydata words only and compare it with the
reference file. Two references were corrected while doing so, before any comparison was counted: the shuffled
correlation p-value of notebook 18 (the reference shuffled the first column, simplydata shuffles the second; with the
same seed the two give 0.1256 and 0.126, and the reference now shuffles the same column), and the correspondence-analysis
shares of notebook 11 (compared unrounded to 1e-10, instead of rounded to one decimal with tolerance 0). Still waiting:
the eigenvalues and total inertia of `add_correspondence_axes` (its `show="numbers"` gives the shares only), and the
partial correlation's p-value (`find_correlation_between(allowing_for=)` refuses `show="result table"`).

### A25. Dunn's test after the A12 fix: one Holm-adjusted p-value below 1 where it must be 1

On InsectSprays (six sprays) `compare_each_pair_of("count", for_each_group_in="spray", test_method="dunn")` now
attaches every p-value to the right pair (A12's 16 FAILs pass), but one pair that passed before now fails:

| pair | simplydata | R PMCMRplus `kwAllPairsDunnTest` (Holm) | Dunn (1964) written out, NumPy + Holm |
|---|---|---|---|
| F - B (and B - F) | 0.926028 | 1 | 1 |

Holm's step-down takes, for each pair, the largest adjusted value of the pairs before it in the order of the raw
p-values (capped at 1), so the pair with the largest raw p-value cannot get less than the pairs before it; both
references give 1. simplydata's 0.926 looks like that pair's value before the step (not yet checked against its raw
p-value). 2 FAILs (notebook 20). Suggestion: take the cumulative maximum over the sorted pairs (statsmodels `multipletests(method="holm")`
does) before attaching the values to the pairs.

### A26. The A21 fix moved 19 comparisons of NIST's second start just under 6 digits

PR #39 gave `fit_formula(test_method="trust region" / "dogbox")` a central-difference Jacobian with a relative step
(`jac="3-point"`, `diff_step = eps^(1/3)`). That brought 120 of A21's FAILs to 6 digits or more, but 19 comparisons
that passed on 47594a0 now land just under the line (LRE against the NIST certified values):

| problem | method | quantities | LRE now |
|---|---|---|---|
| Thurber | trust region and dogbox | standard errors of b2 to b7 (12) | 5.77 to 5.85 |
| MGH09 | trust region | b4 estimate, b4 standard error | 5.98, 5.99 |
| MGH09 | dogbox | b2 standard error, b4 estimate, b4 standard error | 5.90, 5.98, 5.98 |
| ENSO | trust region and dogbox | b6 estimate | 5.94 |

(Thurber b2 standard error: 39.5712219 for the certified 39.5711561.) The estimates of MGH09 and ENSO move too (about
1e-6 relative), not only standard errors; the cause is not yet pinned down. Notebook 27; the Levenberg-Marquardt
default is not affected. Suggestion: as in A21, tighten `ftol`, `xtol` and `gtol`
for these two methods too.

### C26. R irr `kripp.alpha` with complete data (the weighting of A14), now that simplydata no longer calls it

simplydata 3fd3b85 computes Krippendorff's alpha with R `icr::krippalpha` and agrees with Krippendorff's
coincidence-matrix formula written out (A14 fixed): `anxiety` 0.170099, `diagnoses` 0.433410. R irr `kripp.alpha`,
kept in notebook 21 as a second reference, gives 0.163066 and 0.430878, for the reason given in A14 (each pair of
values weighted by 1 instead of 1/(m - 1) when no cell is empty). 2 FAILs, explained by the reference.

### C27. R irr `icc` interval for the average of raters with absolute agreement (the interval of A15)

simplydata 3fd3b85 gives McGraw & Wong's (1996) interval for ICC(A,k) on `anxiety`, (-0.126583, 0.745149), equal to
the formula written out and to R psych `ICC` (A15 fixed). R irr `icc`, kept as a reference, gives (-0.137308,
0.746326), for the reason given in A15. 2 FAILs, explained by the reference; the estimate, F and p-value agree.

### C28. R `lowess(delta = 0)` is another line than simplydata's, now that simplydata uses R's default delta

A16 is fixed: simplydata 3fd3b85 passes R's default `delta` (1 % of the x range) to statsmodels, and its LOWESS line on
mtcars (weight, mpg), read with `chart.get_numbers_drawn()`, agrees with R `lowess()` to 1.5e-10. The second reference
of notebook 24, R `lowess(delta = 0)`, fits every point and differs by up to 2.9e-5 relative (the same 2.9e-5 as A16,
the other way round). 1 FAIL, explained by the reference.

**Notebook 24 with `chart.get_numbers_drawn()`.** The other three chart FAILs repeat known findings with the same
numbers: A17 (the smoothing spline on mtcars weight, where x repeats, 1.17e-5 against R `smooth.spline(lambda = 0.01)`;
Lake Huron, without repeats, agrees to 2.0e-8), C10 (the gamma curve against R `dgamma` at the fitdistrplus fit,
1.43e-6; against simplydata's own maximum likelihood fit it agrees to 1.0e-15) and A7 (the ryegrass dose-response curve
against a tight SciPy fit, 7.8e-6). Still waiting: the best cut-off of `trade_off_curve`, written on the chart as text
(the marked point itself is read and agrees with Youden's index).

**Totals of the second pass** (simplydata 3fd3b85, GitHub Actions run 45072494056 of the pull request; notebooks 01 to 32; on the per-notebook runs the machine differences were 16, with 7 more agrees): 26,478 comparisons, 26,067
agree, 110 simplydata error – being fixed, 48 explained difference – reference, 23 explained difference – machine, 223
explained difference – precision limit, 0 difference not yet explained, 5 waiting for word, 2 skipped. Of the 153
"waiting for word" results of 47594a0: 143 agree, 2 are simplydata errors already known (A7, A17), 2 are explained by
the reference (C10, C28), 5 still wait (3 in notebook 11, 1 in 18, 1 in 24, above) and 1 reference row is no longer
compared (Spearman with `exact = FALSE`, under A9). Most of them now compare a whole drawn layer or table point by
point, so there are more comparisons than before (23,438 on 47594a0).


## Third pass: simplydata main a5a8ff0 (PR #42), 2026-11-23

**The last five checks that waited for a word now compare.** simplydata PR #42 gives back the numbers they waited
for: the correspondence-analysis eigenvalues and total inertia of varespec (`add_correspondence_axes(show="numbers")`,
notebook 11: 0.524932, 0.356798 and 2.083198, against Oksanen's printed 0.5249, 0.3568 and 2.083), the partial
correlation's p-value (`find_correlation_between("mpg", "wt", allowing_for="hp", show="result table")`, notebook 18:
1.1196471362000158e-06, equal to R ppcor `pcor.test` and within 2.7e-14 of the t test of the R `lm` residuals on n - 3)
and the best cut-off of `trade_off_curve` (the `value` column of the marked dot in `chart.get_numbers_drawn()`,
notebook 24: 0.205, against Youden's index over the midpoints and R pROC `coords(..., "best")`). All five agree.
"waiting for word" is now 0.

**What the newer simplydata changed elsewhere.** The pin moved from 3fd3b85 to a5a8ff0, which brings simplydata PRs
#37 (the A1, A4 to A7 fixes), #40, #41 (A25, A26) and #42. On both Ubuntu runs of the pull request, 71 FAILs that were
listed as worth fixing now agree and are removed from `results/known_fails.csv`: A7 37, A26 16, A1 6, A6 5, A21 5, A25 2.
Thirteen listed rows of B1 (NIST nonlinear standard errors near the 6-digit line) and C21 pass on one of the two runners
and fail on the other; they depend on the runner machine and are kept, as is every row seen on another system.

The eight `find_model_averages_of(outcome="yes or no", grouped_by="herd")` rows of notebook 26 are B25 as before:
simplydata's values and R's are the same to the last digit at 3fd3b85 and at a5a8ff0; only the reference's label changed
(emmeans 2.0.4 to 2.0.5 on CRAN), and the label is part of a comparison's key.

### C29. References run with R's default stopping rule; simplydata now converges tighter (its A7 fix)

simplydata's A7 fix (PR #37, commit f9e2150) runs R's iterative fits to convergence, with stopping rules only:
metafor `rma` with REML or Paule-Mandel gets `control = list(tol = 1e-10, threshold = 1e-10, maxiter = 1000)` (metafor's
default: tol and threshold 1e-5, maxiter 100), drc `drm` gets `control = drmc(relTol = 1e-15, maxIt = 10000)` (default
relTol 1e-7, maxIt 500) and extRemes `fevd` for the generalised Pareto gets `optim.args = list(method = "BFGS", control =
list(reltol = 1e-12, maxit = 10000))` (optim's default reltol is about 1.5e-8). The references that call these functions
with their defaults now differ from simplydata by 1e-6 to 1e-3 relative, more than their tolerance of 1e-6.

To tell whether the reference or simplydata stopped short, each of these comparisons got a second reference: the same R
call run to convergence (metafor `control = list(tol = 1e-12, threshold = 1e-12, maxiter = 1000)`, drc
`drmc(relTol = 1e-12, maxIt = 10000)`; notebook 13 already had extRemes at `reltol = 1e-14`). The default-rule rows are
kept. For 77 of the 97 new FAILs simplydata agrees with the converged reference to 1e-8 or better (largest: 2.0e-9 in
notebook 08, 3.7e-10 in 31, 0 in 13, and under 1e-8 for drc's 5-parameter log-logistic fit in 27): these are explained
by the reference stopping early, and listed as **explained difference – precision limit**. `make_known_fails.py` names
exactly these rows (C29_ROWS); any other FAIL against a default-rule reference stays unexplained.

| notebook | comparisons | reference with R's default rule | converged reference | simplydata against the converged one |
|---|---|---|---|---|
| 08 combining studies | 38 | R metafor rma | R metafor rma (tol 1e-12) | at most 2.0e-9 |
| 31 combining, planning, surveys | 14 | R metafor escalc + rma (REML, test = adhoc) | the same, tol 1e-12 | at most 3.7e-10 |
| 13 circular, spatial, extremes | 2 | R extRemes return.level (default stopping rule) | R extRemes return.level (reltol 1e-14) | 0 |
| 27 fitting | 23 | R drc drm + ED (help-page ryegrass example) | R drc drm + ED (relTol 1e-12) | at most 1e-8 |

**20 drc comparisons are not C29** (14: 4, against drc's default rule; 27: 16 of the Weibull fit): simplydata also
differs from drc run to relTol 1e-12, by 1.3e-8 to 2.0e-6 relative. They are explained by C30 below:

| notebook | dataset | method | quantity | simplydata against drc relTol 1e-12 |
|---|---|---|---|---|
| 14_small_additions | ryegrass (drc) | fit_dose_response | b estimate | 1.18e-06 |
| 14_small_additions | ryegrass (drc) | fit_dose_response | c estimate | 1.98e-06 |
| 14_small_additions | ryegrass (drc) | fit_dose_response | ED10 estimate | 1.35e-06 |
| 14_small_additions | ryegrass (drc) | fit_dose_response | ED90 estimate | 3.93e-07 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | b standard error | 1.12e-07 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | c standard error | 1.29e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | d standard error | 1.21e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | e estimate | 1.90e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | e standard error | 8.05e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED10 | 1.87e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED10 standard error | 8.72e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED10 95 % low | 2.27e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED10 95 % high | 3.74e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED50 | 1.89e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED50 standard error | 1.36e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED50 95 % low | 1.96e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED50 95 % high | 1.84e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED90 | 1.91e-08 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED90 standard error | 1.62e-07 |
| 27_fitting | ryegrass (drc) | fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90]) | ED90 95 % high | 4.41e-08 |

**An independent optimum (not drc).** To find out whether drc's reference or simplydata stopped short, the same
least-squares problems were solved to convergence by two other optimisers, both started from drc's relTol 1e-12
estimates (`scripts/references/drc_optimum.py`): R `optim` (BFGS, reltol 1e-15, with the analytic gradient; its Hessian
gives the standard errors as drc defines them, sigma^2 times the inverse of half the Hessian of the residual sum of
squares) and SciPy `least_squares` (xtol = ftol = gtol = 1e-15, analytic Jacobian). The effective doses are drc's
formulas written out, with the delta method for their standard errors and drc's t interval; at drc's own estimates they
give drc's `ED()` to 3e-16. The two optimisers agree with each other to 2.3e-9 (LL.4) and 1.0e-8 (Weibull) in the
parameters and to 2.8e-9 and 7.6e-9 in the doses, so the optimum is determined to about 1e-8. The residual sum of
squares at each set of estimates (in R; the reference scripts write it to
`references/data/14_small_additions/ryegrass_ll4_rss.csv` and `references/data/27_fitting/ryegrass_weibull_rss.csv` on
every run, a folder that is not stored here):

| estimates | LL.4 (notebook 14) | Weibull W1.4 (notebook 27) |
|---|---|---|
| drc, default stopping rule (relTol 1e-7) | 5.400214628430577 | 6.024150858627328 |
| drc, relTol 1e-12 (the reference) | 5.400214627833775 | 6.024150853791260 |
| drc, relTol 1e-15 = simplydata (the same estimates to the last digit) | 5.400214627775889 | 6.024150853781385 |
| SciPy curve_fit (LL.4) | 5.400214627771515 | |
| R optim BFGS (reltol 1e-15) | 5.400214627771515 | 6.024150853086786 |
| SciPy least_squares (1e-15) | 5.400214627771515 | 6.024150853086787 |

simplydata's estimates have a lower residual sum of squares than drc's relTol 1e-12 reference in both fits, so the
reference has not converged in its parameters. But simplydata is not at the optimum either. On LL.4 it is 3.0e-8 to
6.8e-7 from the independent optimum (drc at relTol 1e-12: 2.1e-7 to 1.3e-6); on the Weibull fit 2.4e-7 to 5.7e-6, about
as far as drc at relTol 1e-12. Both are far outside the optimisers' own spread of about 1e-8, so these 20 rows do not
meet the rule for C29 (simplydata at the independent optimum within the optimisers' spread) and stay unexplained; nor is
simplydata farther from the optimum than the reference it fails against. The cause is in drc: even at relTol 1e-15, drc's
optimiser stops 1e-7 to 6e-6 short of the least-squares optimum in the parameters, and simplydata calls drc. The
independent optimum is now compared in notebooks 14 and 27 at the tolerance for iterative fits (1e-6, README). On LL.4
simplydata passes for every estimate and dose; two standard errors miss (b 1.4e-6, ED10 1.5e-6; at the same point drc's
numerical Hessian and optim's differ by about 1e-6). On the Weibull fit 22 of 27 comparisons miss, by up to 5.7e-6 for the
estimates and 2.2e-5 for the standard errors. Manish decided how to label them on 2026-11-23: C30.

### C30. R drc stops short of the least-squares optimum even at relTol 1e-15; simplydata gives exactly what drc gives

R drc stops short of the least-squares optimum even at relTol 1e-15; simplydata gives exactly what drc gives (no
refinement code, by Manish's decision of 2026-11-23, the same rule as for Wampler), and is nearer the optimum than drc's
default and relTol 1e-12 answers. The evidence is above: two independent optimisers run to 1e-15, R `optim` (BFGS,
analytic gradient) and SciPy `least_squares` (analytic Jacobian), agree with each other to 2.3e-9 (LL.4, notebook 14)
and 1.0e-8 (Weibull, notebook 27) in the parameters, so the optimum is known to about 1e-8; the residual sums of squares
order drc's default rule above drc at relTol 1e-12, above simplydata (= drc at relTol 1e-15), above the optimum, in both
fits (the table above). simplydata's distance from the optimum (up to 6.8e-7 on LL.4 and 5.7e-6 on the Weibull fit, 2.2e-5
for a Weibull standard error) is drc's own.

Labelled **explained difference – precision limit**, by the rows `make_known_fails.py` lists (C30_ROWS), 51 in all: the 20
comparisons with drc's default rule and with drc at relTol 1e-12 above, the 24 with the independent optimum (2 standard
errors of LL.4, 22 rows of the Weibull fit), and the 7 rows of notebook 14 against drc at relTol 1e-12 (b and c estimates,
b standard error, ED10 estimate, standard error and interval) that were listed as A7. They have the same cause: simplydata
is nearer the optimum than that reference, and its own distance from the optimum is drc's. Of the 13 A7 rows notebook 14
had, the other 6 now agree and were removed.

**Totals of the third pass** (simplydata a5a8ff0, GitHub Actions run 45169408948 of the pull request; notebooks 01 to
32): 26,718 comparisons, 26,265 agree, 29 simplydata error – being fixed, 48 explained difference – reference, 27
explained difference – machine, 347 explained difference – precision limit (77 of them C29, 51 C30), 0 difference not
yet explained, 0 waiting for word, 2 skipped. Before (simplydata 3fd3b85, run 45072494056):
26,478 comparisons, 26,067 agree, 110 being fixed, 48 reference, 23 machine, 223 precision limit, 0 not yet explained, 5
waiting for word, 2 skipped. There are 240 more comparisons: the 5 checks that waited, the converged references added
next to the default-rule ones and the independent optimum of the two drc fits. The push run of the same commit gave the
same counts; on other runner machines the B1 rows near the 6-digit line move between agree, machine and precision limit
(run 45152693369: 26,268 agree, 18 machine).


## Fourth pass: the last "being fixed" rows looked into, simplydata a5a8ff0, 2026-11-23

The 29 rows still labelled simplydata error – being fixed on the third pass (A5 4, A8 10, A17 1, A21 14) were traced in
simplydata's code (main 050d199; the code of these words is the same as at a5a8ff0). None of them needs a change in
simplydata: simplydata gives the numbers of the package it calls, and says which method it uses. They are relabelled
below. simplydata PR #46 adds tests that pin A17 and Lanczos1 to the package underneath (SciPy `make_smoothing_spline` to
1e-12; SciPy `curve_fit` with simplydata's options).

### C31. Textbook Wald intervals of the attributable fractions; simplydata gives R epiR's intervals and says so

Notebook 16, Willeberg's dry cat food and FLUTD table (was A5). The reference is the textbook Wald interval; simplydata
gives the intervals R epiR `epi.2by2` gives for one table and names them in its printout (simplydata PR #37):

| quantity | simplydata (= R epiR) | Wald, textbook formula |
|---|---|---|
| AF 95 % CI low | 0.329129 | 0.301058 |
| AF 95 % CI high | 0.907214 | 0.910915 |
| PAF 95 % CI low | 0.323759 | 0.036138 |
| PAF 95 % CI high | 0.749138 | 0.782377 |

The same rows against R `epi.2by2` agree to 1e-16. Explained difference – reference (another definition), 4 rows.

### C32. Two exact percentile intervals pick different order statistics (SciPy Thompson against DescTools)

Notebook 18, airquality Ozone by month (was A8). simplydata gives SciPy `quantile_test`'s interval (Thompson 1936: each
end leaves at most 2.5 % beyond it) and its help says that R DescTools `QuantileCI(method = "exact")` picks, for percents
other than 50, the nearby pair whose coverage is closest above the level (simplydata 54dbbf9). Both are exact; for the
median they are the same (every month agrees).

| quantity | simplydata (= SciPy) | DescTools |
|---|---|---|
| p25 95 % CI high, months 5, 6, 7, 9 | 16, 29, 52, 21 | 14, 23, 50, 23 |
| p75 95 % CI high, month 9 | 78 | 73 |
| p75 95 % CI low, months 5, 6, 7, 8, 9 | 19, 21, 63, 65, 24 | 23, 23, 64, 66, 23 |

Explained difference – reference (another rule), 10 rows.

### C33. R smooth.spline lands slightly off the penalised optimum; SciPy and simplydata reach a lower penalised sum

Notebook 24, mtcars weight and mpg, `draw_smooth_line(test_method="spline", smoothness=0.01)` against R
`smooth.spline(lambda = 0.01)` (was A17): largest relative difference 1.17e-5 over 200 values. Both minimise
sum w (y - g(x))^2 + lambda * integral g''^2 with x mapped onto 0 to 1. Checked with R 4.6.2 and SciPy 1.18.1:

- **Repeated x are not the cause.** R collapses repeated x as simplydata does (the average y, with weight = how many;
  default weights are not rescaled). With every weight made distinct (+ 1e-4 per row) R and SciPy still differ by
  1.171e-5 at lambda = 0.01 and 9.6e-6 at lambda = 0.001.
- **R is the one off the optimum.** The penalised sum of SciPy's curve is 209.6387742 against R's 209.6387766 at
  lambda = 0.01, and 194.7932780 against 194.7932853 at lambda = 0.001 (each computed from the curve's own fitted values
  and the exact integral of its second derivative squared).
- simplydata's curve is SciPy `make_smoothing_spline`'s at the same lambda, to 1e-12 (simplydata PR #46).

Lake Huron (no repeated x) agrees to 2.0e-8 at its smoothness. Explained difference – reference, 1 row.

### NIST Lanczos1 with trust region and dogbox: B1 (the double-precision limit), not A21

Notebook 27, NIST Lanczos1 from start 2 (14 rows, were A21): the standard errors and the residual standard deviation
agree with NIST to 3.34 digits (trust region) and 3.57 (dogbox); the estimates agree to 10.6. The Levenberg-Marquardt
default gives 3.5 on the same quantities, listed as B1 since 0.2.0.

- The certified residual standard deviation is 8.9e-14: the data were generated to 14 digits, so the residuals are the
  rounding of y to 13 digits. Evaluating the formula in double precision moves each residual by about 1e-16, about
  0.1 % of its size. At the same fitted values, adding the three exponential terms in another order gives a residual
  standard deviation of 8.9075e-14, 8.9050e-14 or 8.9139e-14: it moves in its 4th digit, which is 3.3 to 3.6 digits.
  The standard errors scale with it.
- simplydata's values are SciPy `curve_fit`'s with the options simplydata passes (the same on Windows to the digits
  shown); SciPy's own `x_scale="jac"` gives 3.0 to 3.1. Exact derivatives would not change the residual standard
  deviation, and simplydata does not compute its own (PLAN.md, section 7).

Listed as B1, explained difference – precision limit (`make_known_fails.py`, LANCZOS1_B1_ROWS).

### The 22 older A21 and A26 rows of notebook 27: removed

`results/known_fails.csv` still listed 17 A21 rows (Bennett5, Lanczos2, Lanczos3, MGH10, MGH17) and 5 A26 rows (MGH09
trust region b4 estimate and standard error, Thurber dogbox b3, b4 and b6 standard errors; listed as A21 before the A26
rule, which now gives them A26), kept because they had been seen on another system or runner. Each was removed only after
it passed on all three runs of this pull request's commits with simplydata a5a8ff0: the pull request's run 45184089331
(Intel Xeon 8370C), the push run 45184049817 (AMD EPYC 9V45) and run-windows 45232851547 (windows-latest, AMD; notebook 27
only). Their lowest LRE on the three is 6.52 (Lanczos3 dogbox b1, Ubuntu); the Thurber standard errors are 6.57 to 6.80.
On Windows notebook 27 has no FAIL that `results/known_fails.csv` does not explain.

`run-windows` is a new workflow: it sets up Python, R and the R packages as `notebooks` does and runs the notebooks named
in `.github/windows-notebooks.txt` on windows-latest, so the Windows rows of `results/known_fails.csv` can come from
GitHub. `scripts/run_notebooks.py` now checks only the notebooks it was asked to run, so the stored Ubuntu results of
the others are not read as Windows FAILs.

**Totals of the fourth pass** (simplydata a5a8ff0, GitHub Actions run 45184089331 of the pull request; notebooks 01 to
32): 26,718 comparisons, 26,268 agree, 0 simplydata error – being fixed, 64 explained difference – reference, 18
explained difference – machine, 366 explained difference – precision limit, 0 difference not yet explained, 0 waiting
for word, 2 skipped. Before (third pass, run 45169408948): 26,265 agree, 29 being fixed, 48 reference, 27 machine, 347
precision limit. The 15 rows of C31 to C33 move to reference and the 14 Lanczos1 rows to precision limit; the B1 rows
near the 6-digit line move between agree, machine and precision limit with the runner machine.

## Fifth pass: the confirming rerun with simplydata main 3f0192b (PR #56), 2026-11-24

The pin moved from a5a8ff0 to 3f0192b, the combined statistics pull request of simplydata (lesson fixes, a note about
repeated rows after some statistics results, other names for some words, the renames of `turn_rows_into_columns` to
`spread_into_columns` and `turn_columns_into_rows` to `gather_into_rows`, `in_sheet=`). `00_summary` now uses
`spread_into_columns`; no other notebook used a renamed word. Every notebook ran again on GitHub Actions (the pull
request's run 45443130321, AMD EPYC 9V74), and notebook 27 on windows-latest (run-windows 45443124789): no new FAIL.

Every comparison keeps its label but 4: the cbpp model averages of `find_model_averages_of (outcome="yes or no",
grouped_by="herd")` in notebook 26 (B25) move from agrees to machine on this runner, as B25 says they do; the
simplydata values move by less than 1e-5 relative. The other changes of simplydata's values are in the last digits
(at most 4e-12 relative). The run references were made again with the R packages of the day (survival 3.8.13,
segmented 2.3.0 among them); their values moved in the last digits only, and no label changed.

**Totals of the fifth pass** (simplydata 3f0192b, GitHub Actions run 45443130321 of the pull request; notebooks 01 to
32): 26,718 comparisons, 26,264 agree, 0 simplydata error – being fixed, 64 explained difference – reference, 22
explained difference – machine, 366 explained difference – precision limit, 0 difference not yet explained, 0 waiting
for word, 2 skipped.
