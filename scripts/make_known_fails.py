"""Write results/known_fails.csv: every FAIL that FINDINGS.md explains, with its finding id and platform.

    python scripts/make_known_fails.py WINDOWS_all_comparisons.csv UBUNTU_all_comparisons.csv [MORE_UBUNTU ...]

Each argument is a results/all_comparisons.csv from a full run on that system (the Ubuntu ones from the
GitHub Actions artifact "executed-notebooks"; give several when GitHub's runner machines give different
results). A FAIL found on Windows and on Ubuntu gets platform "any"; on one only, "windows" or "ubuntu".
"same on every runner" is "no" for a FAIL that turned up on some Ubuntu runs but not others; such a row
is never reported as fixed when it passes. Each FAIL is given its finding id (A1 ... C5 in FINDINGS.md) by the rules
in finding_id(); a FAIL that no rule explains stops the script, so nothing is listed unexplained.
For maintainers, after FINDINGS.md has been updated.

Each row also gets the words a reader sees (column "result"), a one-line explanation and a link to its FINDINGS.md
section (see result_of()); the machine column "finding" keeps the class (A, B, C) and number.

    python scripts/make_known_fails.py --labels    rewrites only "result", "explanation" and "link" of the rows there
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
KEY = ["notebook", "dataset", "method", "quantity", "reference source"]
FINDINGS_URL = "https://github.com/bioforkgithub/simplydata_tests/blob/main/FINDINGS.md"

# The words a reader sees for each comparison are defined in notebooks/checks.py (which prints them and writes them
# to the column "result"); the machine column "status" keeps PASS / FAIL / SKIP and "finding" the FINDINGS.md id.
sys.path.insert(0, str(ROOT / "notebooks"))
from checks import MACHINE, PRECISION, REFERENCE, SIMPLYDATA_ERROR  # noqa: E402

MACHINE_FINDINGS = {"B7"}   # Monte Carlo draws that differ between systems with the same seed
# A findings that were looked into and kept as they are, with their "Resolution" paragraph in FINDINGS.md: the limit of
# double precision, not a bug left to fix (A4: simplydata PR #37)
RESOLVED_AS_PRECISION = {"A4"}
# C findings where the reference, not simplydata, stopped short of the optimum: the difference is a limit of the
# reference's precision (C29: R's default stopping rules; simplydata converges tighter since its A7 fix), or of the
# package simplydata calls (C30: drc stops short of the least-squares optimum; Manish's decision of 2026-11-23)
REFERENCE_PRECISION = {"C29", "C30"}


def result_of(finding: str, platform: str = "any", same_on_every_runner: str = "yes") -> str:
    """The reader's words for an explained FAIL: A -> simplydata error (being fixed); C -> the reference; B ->
    machine when it shows on one system only or not on every runner (or is B7), otherwise the precision limit. An A
    finding resolved as a precision limit (RESOLVED_AS_PRECISION) is the precision limit."""
    if finding in RESOLVED_AS_PRECISION or finding in REFERENCE_PRECISION:
        return PRECISION
    if finding.startswith("A"):
        return SIMPLYDATA_ERROR
    if finding.startswith("C"):
        return REFERENCE
    if finding in MACHINE_FINDINGS or platform != "any" or same_on_every_runner == "no":
        return MACHINE
    return PRECISION


def _slug(heading: str) -> str:
    import re

    text = re.sub(r"[^\w\- ]", "", heading.strip().lstrip("#").strip().lower())
    return text.replace(" ", "-")


def explanations() -> dict:
    """{finding id: (one-line explanation, link)} read from FINDINGS.md: the section heading of the finding, or for a
    finding in a table (C1 ... C11) its row's "comparison" cell, linked to the section that holds the table."""
    import re

    found: dict = {}
    heading = ""
    for line in (ROOT / "FINDINGS.md").read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            heading = line
            m = re.match(r"###\s+([ABC]\d+)\.\s+(.*)", line)
            if m and m.group(1) not in found:
                found[m.group(1)] = (m.group(2).replace("`", ""), f"{FINDINGS_URL}#{_slug(line)}")
        m = re.match(r"\|\s*([ABC]\d+)\s*\|\s*([^|]+)\|", line)
        if m and m.group(1) not in found:
            found[m.group(1)] = (m.group(2).strip().replace("`", ""), f"{FINDINGS_URL}#{_slug(heading)}")
    return found


def add_labels(table: pd.DataFrame) -> pd.DataFrame:
    words = explanations()
    table = table.copy()
    table["result"] = [result_of(f, p, s) for f, p, s in
                       zip(table["finding"], table["platform"], table["same on every runner"])]
    table["explanation"] = [words.get(f, ("", ""))[0] for f in table["finding"]]
    table["link"] = [words.get(f, ("", ""))[1] for f in table["finding"]]
    missing = sorted(set(table.loc[table["explanation"] == "", "finding"]))
    if missing:
        raise SystemExit(f"no FINDINGS.md section found for {missing}")
    return table


def _without_version(ref: str) -> str:
    """The reference source without the package version the reference script adds at its end ("(drc 4.0.0)")."""
    import re

    return re.sub(r" \([A-Za-z][\w.]* [\d.-]+\)$", "", ref)


# Third pass (simplydata a5a8ff0): comparisons with a reference run with R's default stopping rule where simplydata
# agrees to 1e-8 or better with the same R function run to convergence (the tight reference next to it), checked on
# the pull request's run. Only these rows; any other FAIL against a default-rule reference stays unexplained.
C29_ROWS = {  # (notebook, dataset, method, default-rule reference): quantities within 1e-8 of the tight reference
    ('08_combining_studies', 'BCG trials', 'combine_studies (REML, risk differences)', 'R metafor rma'):
        ['95 % CI high', '95 % CI low', 'estimate', 'p-value', 'prediction interval high', 'prediction interval low', 'standard error', 'tau squared'],
    ('08_combining_studies', 'BCG trials', 'combine_studies (for_each_group_in alloc)', 'R metafor rma'):
        ['random: tau squared', 'systematic: tau squared'],
    ('08_combining_studies', 'BCG trials', 'combine_studies (using ablat)', 'R metafor rma'):
        ['intercept: 95 % CI high', 'intercept: 95 % CI low', 'intercept: estimate', 'intercept: p-value', 'per 1 more "ablat": 95 % CI high', 'per 1 more "ablat": 95 % CI low', 'per 1 more "ablat": estimate', 'per 1 more "ablat": p-value'],
    ('08_combining_studies', 'passive smoking', 'combine_studies (REML, Hartung-Knapp)', 'R metafor rma'):
        ['95 % CI high', 'I squared %', 'estimate', 'p-value', 'prediction interval high', 'prediction interval low', 'standard error', 'tau squared'],
    ('08_combining_studies', 'Curtis 1998', 'combine_studies (ratios of averages)', 'R metafor rma'):
        ['95 % CI high', '95 % CI low', 'I squared %', 'estimate', 'p-value', 'prediction interval high', 'prediction interval low', 'standard error', 'tau squared'],
    ('08_combining_studies', 'McDaniel 1994', 'combine_studies (correlations)', 'R metafor rma'):
        ['p-value', 'prediction interval low', 'tau squared'],
    ('13_circular_spatial_extremes', 'Coles rain (ismev)', 'find_return_levels_of (values_above=30)', 'R extRemes return.level (default stopping rule)'):
        ['10-year level', '100-year level'],
    ('27_fitting', 'ryegrass (drc)', 'fit_dose_response (of_kind="5 parameter log-logistic", doses_for_effect=[10, 50, 90])', 'R drc drm + ED (help-page ryegrass example)'):
        ['ED10', 'ED10 95 % high', 'ED10 95 % low', 'ED10 standard error', 'ED50', 'ED50 95 % high', 'ED50 95 % low', 'ED50 standard error', 'ED90', 'ED90 95 % high', 'ED90 95 % low', 'ED90 standard error', 'b estimate', 'b standard error', 'c standard error', 'd standard error', 'e estimate', 'e standard error', 'f estimate', 'f standard error'],
    ('27_fitting', 'ryegrass (drc)', 'fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90])', 'R drc drm + ED (help-page ryegrass example)'):
        ['b estimate', 'c estimate', 'd estimate'],
    ('31_combining_planning_surveys', 'Hackshaw (1997), passive smoking', 'combine_studies (estimates= with interval columns, a ratio)', 'R metafor escalc + rma (REML, test = adhoc)'):
        ['estimate', 'high', 'standard error', 'tau squared'],
    ('31_combining_planning_surveys', 'Hackshaw (1997), passive smoking', 'combine_studies (standard_errors=, estimates_are_ratios=False)', 'R metafor escalc + rma (REML, test = adhoc)'):
        ['estimate', 'high', 'standard error', 'tau squared'],
    ('31_combining_planning_surveys', 'Curtis et al. (1998), CO2', 'combine_studies (averages=, standard_deviations=, show="ratios of averages")', 'R metafor escalc + rma (REML, test = adhoc)'):
        ['estimate', 'high', 'low', 'standard error', 'tau squared'],
    ('31_combining_planning_surveys', 'McDaniel et al. (1994), interviews', 'combine_studies (correlations=, counts=)', 'R metafor escalc + rma (REML, test = adhoc)'):
        ['tau squared'],
}

# Third pass, Manish's decision of 2026-11-23: the drc fits of notebooks 14 (LL.4) and 27 (Weibull). R drc stops short of
# the least-squares optimum even at relTol 1e-15 (R optim and SciPy least_squares, both at 1e-15, agree with each other
# to 2.3e-9 and 1.0e-8); simplydata gives exactly what drc gives and is nearer the optimum than drc's default and
# relTol 1e-12 answers. Only these rows.
C30_ROWS = {  # (notebook, dataset, method, reference without its version): quantities (the drc fits; see C30)
    ('14_small_additions', 'ryegrass (drc)', 'fit_dose_response', 'R drc drm (default stopping rule)'):
        ['b estimate', 'c estimate'],
    ('14_small_additions', 'ryegrass (drc)', 'fit_dose_response', 'R optim BFGS (reltol 1e-15, analytic gradient)'):
        ['ED10 standard error', 'b standard error'],
    ('14_small_additions', 'ryegrass (drc)', 'fit_dose_response', 'R drc ED (default stopping rule)'):
        ['ED10 estimate', 'ED90 estimate'],
    ('27_fitting', 'ryegrass (drc)', 'fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90])', 'R optim BFGS (reltol 1e-15, analytic gradient)'):
        ['ED10 95 % high', 'ED10 95 % low', 'ED10 standard error', 'ED50', 'ED50 95 % high', 'ED50 95 % low', 'ED50 standard error', 'ED90', 'ED90 95 % high', 'ED90 95 % low', 'ED90 standard error', 'b estimate', 'b standard error', 'c estimate', 'c standard error', 'e estimate', 'e standard error'],
    ('27_fitting', 'ryegrass (drc)', 'fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90])', 'SciPy least_squares (1e-15)'):
        ['ED50', 'ED90', 'b estimate', 'c estimate', 'e estimate'],
    ('27_fitting', 'ryegrass (drc)', 'fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90])', 'R drc drm + ED (help-page ryegrass example)'):
        ['ED10', 'ED10 95 % high', 'ED10 95 % low', 'ED10 standard error', 'ED50', 'ED50 95 % high', 'ED50 95 % low', 'ED50 standard error', 'ED90', 'ED90 95 % high', 'ED90 standard error', 'b standard error', 'c standard error', 'd standard error', 'e estimate', 'e standard error'],
    ('14_small_additions', 'ryegrass (drc)', 'fit_dose_response', 'R drc ED (delta, relTol 1e-12)'):
        ['ED10 95 % CI high', 'ED10 95 % CI low', 'ED10 estimate', 'ED10 standard error'],
    ('14_small_additions', 'ryegrass (drc)', 'fit_dose_response', 'R drc drm (relTol 1e-12)'):
        ['b estimate', 'b standard error', 'c estimate'],
}


def finding_id(row: pd.Series) -> str | None:
    """The FINDINGS.md section that explains this FAIL, or None."""
    nb, ds, method, ref = row["notebook"], row["dataset"], row["method"], row["reference source"]
    quantity = row["quantity"]
    if quantity in C29_ROWS.get((nb, ds, method, _without_version(ref)), ()):
        return "C29"
    if quantity in C30_ROWS.get((nb, ds, method, _without_version(ref)), ()):
        return "C30"
    is_r = ref.startswith("R")
    if nb == "01_describing_data":
        if ds == "NIST NumAcc4" and method == "find_standard_deviation_of" and ref == "NIST certified":
            return "B3"
        if ds in ("NIST NumAcc3", "NIST NumAcc4") and method == "find_autocorrelation_of" and is_r:
            return "C1"
    if nb == "02_comparing_groups" and method == "compare_averages_of (ANOVA)":
        if ds == "NIST AtmWtAg" and is_r:
            return "C2"
        if ds in ("NIST SmLs07", "NIST SmLs08", "NIST SmLs09"):
            return "B2"
    is_qr = ref == "statsmodels OLS (QR)"
    if nb == "03_regression_NIST":
        if ds == "NIST Longley":
            return "A1" if ref == "NIST certified" else ("C4" if is_r else None)
        if ds == "NIST Filip":
            return "C4" if is_r else ("C6" if is_qr else None)
        if ds == "NIST NoInt1":
            return "B4"
        if ds in ("NIST Wampler1", "NIST Wampler2"):
            if is_r or is_qr:
                return "C3"
            if ds == "NIST Wampler1" and ref == "NIST certified" and quantity.endswith("standard error"):
                return "A4"
        if ds == "NIST Wampler3":
            if is_r:
                return "C4"
            if ref == "NIST certified" and quantity.endswith("estimate"):
                return "A4"
        if ds in ("NIST Wampler4", "NIST Wampler5"):
            return "B4" if ref == "NIST certified" else ("C4" if is_r else ("C6" if is_qr else None))
    if nb == "04_counts_and_yes_no" and ds == "NMES 1988" and method == "explain (counts that vary a lot)" and is_r:
        return "C5"
    if nb == "06_curves_and_dose_response" and ds.startswith("NIST") and method.startswith("fit_formula"):
        if ds == "NIST BoxBOD" and method == "fit_formula (start 1)" and ref == "NIST certified":
            return "B5"
        if ref.startswith("expected"):
            return None
        return "B1"
    return new_notebooks(nb, ds, method, quantity, ref)


def new_notebooks(nb: str, ds: str, method: str, quantity: str, ref: str) -> str | None:
    """Notebooks 08 onward (methods added after simplydata 0.2.0)."""
    if nb == "08_combining_studies" and ds == "BCG trials":
        if method.startswith("combine_studies (Paule-Mandel") and ref in (
                "R metafor rma (tol 1e-12)", "statsmodels combine_effects (iterated)"):
            return "A7"
        if method == "combine_studies (REML, normal)" and ref == "R metafor rma (tol 1e-12)":
            return "A7"
    if nb == "13_circular_spatial_extremes":
        if (ds in ("Port Pirie (ismev)", "Fremantle (ismev)") and ref == "R extRemes fevd"
                and quantity.startswith("standard error of")):
            return "B6"
        if ds == "Coles rain (ismev)" and method.startswith("find_return_levels_of") and ref in (
                "the formula on the R fit", "R extRemes return.level", "R extRemes ci.fevd (same grid)"):
            return "A7"
    if nb == "14_small_additions":
        if ds == "AirPassengers (log)" and ref == "statsmodels SARIMAX (fixed order)":
            return "C8"
        if ds.startswith("poisons") and ref == "statsmodels OLS + SciPy, the profile (Wald interval)":
            return "B6"
        if ds == "ryegrass (drc)" and ref in ("R drc drm (relTol 1e-12)", "R drc ED (delta, relTol 1e-12)",
                                              "SciPy curve_fit (LL.4)"):
            return "A7"
        if (ds == "iris" and method.startswith("compare_averages_of (MANOVA") and ref == "R summary.manova"
                and quantity.startswith("Hotelling-Lawley trace")):
            return "C9"
    if nb == "10_power_and_sample_size":
        if ds.startswith("d = 0.2") and ref == "Cohen (1988), Table 2.4.1":
            return "C7"
        if ref in ("R pwr.anova.test", "R pwr.f2.test") and quantity in ("power with that many", "power"):
            return "C7"
        if (method in ("find_how_many_needed (two percents)", "find_chance_of_finding (two percents)")
                and quantity in ("power with that many", "power")
                and ref in ("R power.prop.test", "R power.prop.test help page, printed")):
            return "A6"
    if nb == "17_bayesian" and ds == "sleepstudy" and ref.startswith("R rstanarm (4 x 10,000 iterations"):
        return "B7"
    if (nb == "16_cause_and_effect" and ds == "dry cat food and FLUTD"
            and method == "compare_risks (attributable fractions)" and ref == "NumPy, textbook formula"
            and quantity in C31_QUANTITIES):
        return "C31"   # fourth pass: was A5; simplydata gives epiR's one-table intervals and names them (PR #37)
    return gaps_a_notebooks(nb, ds, method, quantity, ref) or gaps_b_notebooks(nb, ds, method, quantity, ref)


def gaps_a_notebooks(nb: str, ds: str, method: str, quantity: str, ref: str) -> str | None:
    """Notebooks 18 to 24 (real-data checks of the words the 2026-11-21 coverage audit listed, stream A)."""
    if nb == "18_summaries_real_data":
        if ds == "NIST NumAcc4" and method == "nums.find_standard_deviation" and ref == "NIST certified":
            return "B3"
        if (method == "find_percentile_of (interval=)" and ref == "R DescTools QuantileCI(exact)"
                and (ds, quantity) in C32_ROWS):
            return "C32"   # fourth pass: was A8; SciPy's Thompson interval, named as such (simplydata 54dbbf9)
        if method == 'find_correlation_between (test_method="spearman")' and ref.startswith("R cor.test (default: exact)"):
            return "A9"
    if nb == "20_tests_counts_pairs_matched":
        if ds == "Convictions (R fisher.test example)" and method in ('compare_counts_of (test_method="barnard")',
                                                                       'compare_counts_of (test_method="boschloo")'):
            return "A10"
        if method == 'compare_counts_of (show="cells")' and ref == "R chisq.test()$stdres":
            return "A11"
        if (method == 'compare_each_pair_of (test_method="dunn")' and ds.startswith("InsectSprays")
                and quantity in ("adjusted p-value F - B", "adjusted p-value B - F")):
            return "A25"   # second pass (simplydata 3fd3b85): the A12 fix leaves Holm's maximum off this pair
        if method == 'compare_each_pair_of (test_method="dunn")' and ds.startswith("InsectSprays"):
            return "A12"
        if method == 'compare_each_pair_of (based_on="ranks")' and ref == "R pairwise.wilcox.test":
            return "A13"
    if nb == "21_agreement_and_diagnostic_tests":
        if method == 'agreement_between (test_method="krippendorff")' and ref.startswith("Krippendorff (2011)"):
            return "A14"
        if (method == 'agreement_between (raters="a sample", checking_for="agreement", of="average of raters")'
                and "CI" in quantity and ref in ("McGraw & Wong (1996), NumPy", "R psych ICC")):
            return "A15"
        if ref == "R irr kappa2 (help page example)" and quantity == "p-value":
            return "C11"
        # second pass (simplydata 3fd3b85): simplydata no longer inherits irr's errors; irr kept as a reference
        if method == 'agreement_between (test_method="krippendorff")' and ref == "R irr kripp.alpha":
            return "C26"
        if (method == 'agreement_between (raters="a sample", checking_for="agreement", of="average of raters")'
                and "CI" in quantity and ref == "R irr icc (help page example)"):
            return "C27"
    if nb == "23_distributions" and method == "table.compare_distributions_for" and ref.startswith("R fitdistrplus"):
        return "C10"
    if nb == "24_charts_and_directions":
        if method == "draw_smooth_line" and ref.startswith("R lowess (its default delta"):
            return "A16"
        if method == "draw_smooth_line" and ref == "R lowess(delta = 0)":
            return "C28"   # second pass: simplydata now uses R's default delta (A16 fixed); delta = 0 is another line
        if (method == 'draw_smooth_line (test_method="spline", smoothness=0.01)' and ds == "mtcars weight and mpg"
                and ref == "R smooth.spline(lambda = 0.01)"
                and quantity == "spline: largest relative difference over 200 values"):
            return "C33"   # fourth pass: was A17; R lands slightly off the penalised optimum
        if method == 'fit_distribution (distribution="gamma")' and ref.startswith("R dgamma at fitdistrplus"):
            return "C10"
        if method == "fit_dose_response" and ref == "SciPy curve_fit (tight)":
            return "A7"
    return None


# Second pass (simplydata 3fd3b85): the comparisons of NIST's second start that passed on 47594a0 and land just under
# 6 digits after the A21 fix (FINDINGS.md A26)
A26_ROWS = {(f"NIST {ds}", f'fit_formula (start 2, test_method="{m}")', q)
            for ds, m, q in
            [("Thurber", m, f"b{i} standard error") for m in ("trust region", "dogbox") for i in range(2, 8)]
            + [("MGH09", "trust region", "b4 estimate"), ("MGH09", "trust region", "b4 standard error"),
               ("MGH09", "dogbox", "b2 standard error"), ("MGH09", "dogbox", "b4 estimate"),
               ("MGH09", "dogbox", "b4 standard error"),
               ("ENSO", "trust region", "b6 estimate"), ("ENSO", "dogbox", "b6 estimate")]}

# Fourth pass (simplydata a5a8ff0, FINDINGS "Fourth pass"): the A findings of the earlier passes that were looked into
# and are not simplydata errors. Only these rows; any other FAIL of the same method stays unexplained.
# C31 (was A5): notebook 16, the textbook Wald intervals of the attributable fractions
C31_QUANTITIES = {"AF 95 % CI low", "AF 95 % CI high", "PAF 95 % CI low", "PAF 95 % CI high"}
# C32 (was A8): notebook 18, the quartiles' exact intervals of airquality Ozone against DescTools
C32_ROWS = {("airquality", f"Ozone {q}, month {m}") for q, m in
            [("p25 95 % CI high", 5), ("p25 95 % CI high", 6), ("p25 95 % CI high", 7), ("p25 95 % CI high", 9),
             ("p75 95 % CI high", 9), ("p75 95 % CI low", 5), ("p75 95 % CI low", 6), ("p75 95 % CI low", 7),
             ("p75 95 % CI low", 8), ("p75 95 % CI low", 9)]}
# B1 (were A21): NIST Lanczos1 from start 2 with trust region and dogbox, the standard errors and the residual standard
# deviation (3.3 to 3.6 digits; the Levenberg-Marquardt default gives 3.5): the double-precision limit of a 9e-14
# residual standard deviation
LANCZOS1_B1_ROWS = {("NIST Lanczos1", f'fit_formula (start 2, test_method="{m}")', q)
                    for m in ("trust region", "dogbox")
                    for q in [f"b{i} standard error" for i in range(1, 7)] + ["residual standard deviation"]}


def gaps_b_notebooks(nb: str, ds: str, method: str, quantity: str, ref: str) -> str | None:
    """Notebooks 25 to 32 (real-data checks of the groups Models, Fitting, Patterns across columns, Forecasts,
    Smoothing and trend, Survival, Combining studies, Planning), in the simplydata-only layout: their reference
    sources name the R function and package version, such as "R MASS polr (help page example) (MASS 7.3.66)".
    Their findings are numbered from 20 (A20, B20, C20 ...) so that they never clash with other notebooks'."""
    if nb == "25_models_outcomes":
        if method == 'explain (table_of_tests="in order")':
            return "A20"
        if method == 'explain (outcome="ordered groups", in_order=)' and ref.startswith("R MASS polr"):
            return "C21"
        if method == 'explain (outcome="groups")' and ref.startswith("R nnet multinom"):
            return "C21"
        if method == 'explain (outcome="positive numbers")' and ref.startswith("R stats glm"):
            return "C21"
    if nb == "26_models_mixed_and_more":
        if method.startswith("compare_models") and quantity.endswith("p-value"):
            return "B22"
        if method.startswith("find_model_averages_of") and 'grouped_by="herd"' in method and ref.startswith("R emmeans"):
            return "B25"
    if nb == "27_fitting" and ds.startswith("NIST") and method.startswith("fit_formula (start 2"):
        if (ds, method, quantity) in A26_ROWS:
            return "A26"
        if "levenberg-marquardt" in method or (ds, method, quantity) in LANCZOS1_B1_ROWS:
            return "B1"
        return "A21"
    if nb == "28_patterns_across_columns" and method.startswith("find_hidden_factors_behind") and ref.startswith("R psych fa"):
        return "B24"
    if nb == "29_forecasts_smoothing_trend":
        if method.startswith("compare_forecasts_for"):
            return "A22"
        if method.startswith("find_autocorrelation_of (on_globe=True") and ref.startswith("R spdep"):
            return "C24"
        if method == 'nums.smooth (test_method="spline", smoothness="auto")':
            return "B23"
    if nb == "31_combining_planning_surveys" and ref.startswith("R pwr pwr.f2.test") and quantity == "power":
        return "C7"
    return None


def main(windows: str, *ubuntu: str) -> int:
    win = pd.read_csv(windows)
    win_fails = win.loc[win["status"] == "FAIL", KEY].drop_duplicates()
    ubuntu_runs = [pd.read_csv(path) for path in ubuntu]
    # a FAIL on any of the Ubuntu runs counts; the runs can come from different runner machines
    counted = pd.concat([t.loc[t["status"] == "FAIL", KEY].drop_duplicates() for t in ubuntu_runs])
    counted = counted.value_counts().rename("ubuntu runs").reset_index()
    both = win_fails.merge(counted, on=KEY, how="outer", indicator=True)
    both["platform"] = both["_merge"].map({"both": "any", "left_only": "windows", "right_only": "ubuntu"})
    # "no" when the FAIL turned up on some Ubuntu runs but not on others (it depends on the machine)
    seen_on_ubuntu = both["ubuntu runs"].fillna(0)
    both["same on every runner"] = ((seen_on_ubuntu == 0) | (seen_on_ubuntu == len(ubuntu_runs))).map(
        {True: "yes", False: "no"})
    both["finding"] = both.apply(finding_id, axis=1)
    unexplained = both[both["finding"].isna()]
    if len(unexplained):
        print("FAILs that no finding explains (add them to FINDINGS.md and to finding_id first):")
        print(unexplained[KEY].to_string(index=False))
        return 1
    out = add_labels(both[KEY + ["platform", "same on every runner", "finding"]].sort_values(KEY))
    out.to_csv(ROOT / "results" / "known_fails.csv", index=False)
    print(f"wrote {len(out)} rows to results/known_fails.csv: "
          + ", ".join(f"{k} {v}" for k, v in out["platform"].value_counts().items()))
    print(out["finding"].value_counts().sort_index().to_string())
    return 0


def labels_only() -> int:
    path = ROOT / "results" / "known_fails.csv"
    table = pd.read_csv(path)
    table = add_labels(table[KEY + ["platform", "same on every runner", "finding"]])
    table.to_csv(path, index=False)
    print(table["result"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(labels_only() if sys.argv[1:] == ["--labels"] else main(*sys.argv[1:]))
