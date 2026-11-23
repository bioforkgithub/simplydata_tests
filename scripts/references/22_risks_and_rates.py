"""Reference values for notebooks/22_risks_and_rates.ipynb, saved to references/22_risks_and_rates.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/22_risks_and_rates/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/22_risks_and_rates.py      (from the repository root; needs R with the packages listed)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks and realdata)

import re

import numpy as np
import pandas as pd
import simplydata as sdata
from scipy import stats

import realdata as rd
from checks import Checks, quiet, r_available, r_run, read_row, read_rows, versions

versions("pandas", "numpy", "scipy", "statsmodels", "rpy2")
checks = Recorder("22_risks_and_rates")
HAS_R = r_available()
CLOSED = 1e-10
ITERATIVE = 1e-6   # values from models fitted by iteration (the propensity scores of compare_balanced_groups_of)
if HAS_R:
    print(rd.r_version_of("survival", "DescTools", "epiR"))


def table_of(word, *args, **options):
    """Every row of a simplydata result table, each as {column: value} (read with table.get_row)."""
    return read_rows(quiet(word, *args, **options))


def keyed(rows, column="measure"):
    """The rows read by table_of, keyed by one column (the measure)."""
    return {row[column]: row for row in rows}


def columns_of(found):
    """The column names of rows read by table_of (a list of rows, or rows keyed by keyed())."""
    first = found[0] if isinstance(found, list) else next(iter(found.values()))
    return list(first)


def garwood(k, t, level=0.95):
    """Garwood's exact interval of a Poisson rate k / t (chi-square percentiles)."""
    a = 1 - level
    low = stats.chi2.ppf(a / 2, 2 * k) / 2 / t if k > 0 else 0.0
    return low, stats.chi2.ppf(1 - a / 2, 2 * k + 2) / 2 / t


def ratio_interval(k1, t1, k2, t2, level=0.95):
    """The exact interval of the rate ratio (k1/t1)/(k2/t2), from the Clopper-Pearson interval of k1/(k1 + k2)."""
    a = 1 - level
    n = k1 + k2
    lo = stats.beta.ppf(a / 2, k1, n - k1 + 1) if k1 > 0 else 0.0
    hi = stats.beta.ppf(1 - a / 2, k1 + 1, n - k1) if k1 < n else 1.0
    return lo / (1 - lo) * t2 / t1, hi / (1 - hi) * t2 / t1


example = rd.r_example_code("stats", "poisson.test")
m = re.search(r"poisson\.test\(c\((\d+),\s*([\d+]+)\),\s*c\((\d+),\s*([\d+]+)\)\)", example)
k1, k2, t1, t2 = (int(eval(g)) for g in m.groups())   # the help page writes sums such as 6+8+7
print("from the help page:", k1, k2, t1, t2)
cities = pd.DataFrame({"city": ["Fredericia", "other cities"], "cases": [k1, k2], "person_years": [t1, t2]})
towns = sdata.table(cities)
towns.compare_rates("cases", time="person_years", for_each_group_in="city", against="other cities")

ds = "lung cancer, Fredericia and other cities (R poisson.test example; Andersen 1977)"


def rate_rows(found):
    rates = {row["group"]: row for row in found if row["measure"] == "rate"}
    ratio = next(row for row in found if row["measure"].startswith("rate ratio"))
    return rates, ratio


for word_name, call in [("table.compare_rates", lambda **o: table_of(towns.compare_rates, "cases", time="person_years",
                                                                     for_each_group_in="city", against="other cities", **o)),
                        ("sdata.compare_rates", lambda **o: table_of(sdata.compare_rates, events=[k1, k2], time=[t1, t2],
                                                                     names=["Fredericia", "other cities"],
                                                                     against="other cities", **o))]:
    for per, interval in [(1000, 95), (100000, 90)]:
        rates, ratio = rate_rows(call(per=per, interval=interval))
        lo_c = [c for c in columns_of(rates) if c.startswith("rate ") and "CI low" in c][0]
        hi_c = [c for c in columns_of(rates) if c.startswith("rate ") and "CI high" in c][0]
        label = f"{word_name} (per={per}, interval={interval})"
        for city, k, t in [("Fredericia", k1, t1), ("other cities", k2, t2)]:
            low, high = garwood(k, t, interval / 100)
            checks.compare(ds, label, f"rate, {city}", rates[city]["rate"], per * k / t, "events / time", CLOSED, "relative")
            checks.compare(ds, label, f"rate {lo_c[5:]}, {city}", rates[city][lo_c], per * low, "Garwood exact, SciPy",
                           CLOSED, "relative")
            checks.compare(ds, label, f"rate {hi_c[5:]}, {city}", rates[city][hi_c], per * high, "Garwood exact, SciPy",
                           CLOSED, "relative")
        lo_r = [c for c in ratio if c.startswith("estimate") and "CI low" in c][0]
        hi_r = [c for c in ratio if c.startswith("estimate") and "CI high" in c][0]
        low, high = ratio_interval(k1, t1, k2, t2, interval / 100)
        checks.compare(ds, label, "rate ratio", ratio["estimate"], (k1 / t1) / (k2 / t2), "ratio of the rates", CLOSED, "relative")
        checks.compare(ds, label, f"rate ratio {lo_r[9:]}", ratio[lo_r], low, "exact conditional interval, SciPy", CLOSED, "relative")
        checks.compare(ds, label, f"rate ratio {hi_r[9:]}", ratio[hi_r], high, "exact conditional interval, SciPy", CLOSED, "relative")
        if HAS_R:
            r = r_run(f"p <- poisson.test(c({k1}, {k2}), c({t1}, {t2}), conf.level = {interval / 100}); list(p = p$p.value)")
            checks.compare(ds, label, "p-value", ratio["p-value"], r["p"], "R poisson.test", CLOSED, "relative")
    for opts, alternative in [({"sides": 1, "higher": "Fredericia"}, "greater"), ({"sides": 1, "lower": "Fredericia"}, "less"),
                              ({"alternative": "greater"}, "greater")]:
        _, ratio = rate_rows(call(**opts))
        p_cond = (stats.binom.sf(k1 - 1, k1 + k2, t1 / (t1 + t2)) if alternative == "greater"
                  else stats.binom.cdf(k1, k1 + k2, t1 / (t1 + t2)))
        label = f"{word_name} (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")"
        checks.compare(ds, label, "p-value", ratio["p-value"], p_cond, "conditional binomial test, SciPy", CLOSED, "relative")
        if HAS_R:
            r = r_run(f'list(p = poisson.test(c({k1}, {k2}), c({t1}, {t2}), alternative = "{alternative}")$p.value)')
            checks.compare(ds, label, "p-value", ratio["p-value"], r["p"], "R poisson.test", CLOSED, "relative")


def tost_rate(low, high):
    def chance(r):
        return r * t1 / (r * t1 + t2)
    return stats.binom.sf(k1 - 1, k1 + k2, chance(low)), stats.binom.cdf(k1, k1 + k2, chance(high))


p_low, p_high = tost_rate(0.5, 2)
for word_name, call in [("table.compare_rates", lambda **o: table_of(towns.compare_rates, "cases", time="person_years",
                                                                     for_each_group_in="city", against="other cities", **o)),
                        ("sdata.compare_rates", lambda **o: table_of(sdata.compare_rates, events=[k1, k2], time=[t1, t2],
                                                                     names=["Fredericia", "other cities"],
                                                                     against="other cities", **o))]:
    _, ratio = rate_rows(call(same_within=[0.5, 2]))
    checks.compare(ds, f"{word_name} (same_within=)", "p-value against the low margin", ratio["p-value against the low margin"],
                   p_low, "exact conditional test written out, SciPy", CLOSED, "relative")
    checks.compare(ds, f"{word_name} (same_within=)", "p-value against the high margin",
                   ratio["p-value against the high margin"], p_high, "exact conditional test written out, SciPy", CLOSED,
                   "relative")
    _, ratio = rate_rows(call(same_within=[0.5, 2], not_lower="Fredericia"))
    checks.compare(ds, f'{word_name} (same_within=, not_lower="Fredericia")', "p-value", ratio["p-value"], p_low,
                   "exact conditional test written out, SciPy", CLOSED, "relative")
    _, ratio = rate_rows(call(same_within=[0.5, 2], not_higher="Fredericia"))
    checks.compare(ds, f'{word_name} (same_within=, not_higher="Fredericia")', "p-value", ratio["p-value"], p_high,
                   "exact conditional test written out, SciPy", CLOSED, "relative")


lung = rd.r_dataset("lung", "survival")
lung["died"] = (lung["status"] == 2).astype(int)
lung["sex"] = np.where(lung["sex"] == 1, "male", "female")
trial = sdata.table(lung)
ds = "NCCTG lung cancer trial (R survival lung)"
totals = lung.groupby("sex")[["died", "time"]].sum()
for per, interval in [(1000, 95), (365.25, 99)]:
    found = keyed(table_of(trial.rate_of, "died", time="time", for_each_group_in="sex", per=per, interval=interval), "group")
    lo_c = [c for c in columns_of(found) if "CI low" in c][0]
    hi_c = [c for c in columns_of(found) if "CI high" in c][0]
    for sex in ("female", "male"):
        k, t = totals.loc[sex, "died"], totals.loc[sex, "time"]
        low, high = garwood(k, t, interval / 100)
        label = f"rate_of (per={per}, interval={interval})"
        checks.compare(ds, label, f"rate, {sex}", found[sex]["rate"], per * k / t, "events / time", CLOSED, "relative")
        checks.compare(ds, label, f"{lo_c}, {sex}", found[sex][lo_c], per * low, "Garwood exact, SciPy", CLOSED, "relative")
        checks.compare(ds, label, f"{hi_c}, {sex}", found[sex][hi_c], per * high, "Garwood exact, SciPy", CLOSED, "relative")
        if HAS_R:
            r = r_run(f"p <- poisson.test({k}, {t}, conf.level = {interval / 100}); list(low = p$conf.int[1], high = p$conf.int[2])")
            checks.compare(ds, label, f"{lo_c}, {sex}", found[sex][lo_c], per * r["low"], "R poisson.test", CLOSED, "relative")
            checks.compare(ds, label, f"{hi_c}, {sex}", found[sex][hi_c], per * r["high"], "R poisson.test", CLOSED, "relative")
found = table_of(trial.compare_rates, "died", time="time", for_each_group_in="sex", against="female")
_, ratio = rate_rows(found)
if HAS_R:
    r = r_run(f"""p <- poisson.test(c({totals.loc['male', 'died']}, {totals.loc['female', 'died']}),
        c({totals.loc['male', 'time']}, {totals.loc['female', 'time']})); list(p = p$p.value, low = p$conf.int[1], high = p$conf.int[2], e = unname(p$estimate))""")
    for quantity, key in [("estimate", "e"), ("estimate 95 % CI low", "low"), ("estimate 95 % CI high", "high"), ("p-value", "p")]:
        checks.compare(ds, 'compare_rates (against="female")', f"rate ratio {quantity}", ratio[quantity], r[key],
                       "R poisson.test", CLOSED, "relative")


ucb = rd.r_dataset("UCBAdmissions")
applicants = pd.DataFrame([(row.Admit, row.Gender, row.Dept) for row in ucb.itertuples() for _ in range(int(row.Freq))],
                          columns=["admit", "sex", "department"])
berkeley = sdata.table(applicants)
berkeley.compare_risks("admit", for_each_group_in="sex", against="Male", value="Admitted")

ds = "UCB admissions 1973 (Bickel et al. 1975)"
counts = applicants.groupby("sex")["admit"].apply(lambda s: (s == "Admitted").sum())
sizes = applicants.groupby("sex").size()
x1, n1, x2, n2 = counts["Female"], sizes["Female"], counts["Male"], sizes["Male"]


def wilson(x, n, z):
    p = x / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def newcombe(x1, n1, x2, n2, level=0.95):
    z = stats.norm.ppf(1 - (1 - level) / 2)
    p1, p2 = x1 / n1, x2 / n2
    l1, u1 = wilson(x1, n1, z)
    l2, u2 = wilson(x2, n2, z)
    d = p1 - p2
    return d - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), d + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2)


def risk_rows(found):
    return keyed(found)


for word_name, call in [("table.compare_risks", lambda **o: table_of(berkeley.compare_risks, "admit", for_each_group_in="sex",
                                                                     against="Male", value="Admitted", **o)),
                        ("sdata.compare_risks", lambda **o: table_of(sdata.compare_risks, events=[int(x1), int(x2)],
                                                                     out_of=[int(n1), int(n2)], names=["Female", "Male"],
                                                                     against="Male", **o))]:
    for interval in (95, 90):
        found = call(interval=interval)
        lo_c = [c for c in columns_of(found) if "CI low" in c and c.startswith("estimate")][0]
        hi_c = [c for c in columns_of(found) if "CI high" in c and c.startswith("estimate")][0]
        level = interval / 100
        z = stats.norm.ppf(1 - (1 - level) / 2)
        label = f"{word_name} (interval={interval})"
        risks = {row["group"]: row for row in found if row["measure"] == "risk"}
        if HAS_R:
            for sex, x, n in [("Female", x1, n1), ("Male", x2, n2)]:
                r = r_run(f"b <- binom.test({x}, {n}, conf.level = {level}); list(low = b$conf.int[1], high = b$conf.int[2])")
                checks.compare(ds, label, f"risk {lo_c[9:]}, {sex}", risks[sex][lo_c], r["low"], "R binom.test", CLOSED, "relative")
                checks.compare(ds, label, f"risk {hi_c[9:]}, {sex}", risks[sex][hi_c], r["high"], "R binom.test", CLOSED, "relative")
        rows = risk_rows(found)
        low, high = newcombe(x1, n1, x2, n2, level)
        diff = rows["risk difference"]
        checks.compare(ds, label, "risk difference", diff["estimate"], x1 / n1 - x2 / n2, "difference of the risks", CLOSED, "relative")
        checks.compare(ds, label, f"risk difference {lo_c[9:]}", diff[lo_c], low, "Newcombe (1998) hybrid score, NumPy", CLOSED, "relative")
        checks.compare(ds, label, f"risk difference {hi_c[9:]}", diff[hi_c], high, "Newcombe (1998) hybrid score, NumPy", CLOSED, "relative")
        rr = (x1 / n1) / (x2 / n2)
        se_rr = np.sqrt(1 / x1 - 1 / n1 + 1 / x2 - 1 / n2)
        odds = (x1 / (n1 - x1)) / (x2 / (n2 - x2))
        se_or = np.sqrt(1 / x1 + 1 / (n1 - x1) + 1 / x2 + 1 / (n2 - x2))
        for measure, value, se in [("risk ratio", rr, se_rr), ("odds ratio", odds, se_or)]:
            row = rows[measure]
            checks.compare(ds, label, measure, row["estimate"], value, "formula, NumPy", CLOSED, "relative")
            checks.compare(ds, label, f"{measure} {lo_c[9:]}", row[lo_c], value * np.exp(-z * se), "Wald on the log scale, NumPy",
                           CLOSED, "relative")
            checks.compare(ds, label, f"{measure} {hi_c[9:]}", row[hi_c], value * np.exp(z * se), "Wald on the log scale, NumPy",
                           CLOSED, "relative")
        if HAS_R and interval == 95:
            r = r_run(f"list(p = prop.test(c({x1}, {x2}), c({n1}, {n2}))$p.value)")
            checks.compare(ds, word_name, "p-value", diff["p-value"], r["p"], "R prop.test (Yates)", CLOSED, "relative")
    if HAS_R:
        found = risk_rows(call(p_value_by="exact"))
        r = r_run(f"""f <- fisher.test(matrix(c({x1}, {n1 - x1}, {x2}, {n2 - x2}), 2, byrow = TRUE)); list(p = f$p.value, or = unname(f$estimate))""")
        p_row = next(row for row in found.values() if row["p-value"] is not None and row["p-value"] == row["p-value"])
        checks.compare(ds, f'{word_name} (p_value_by="exact")', "p-value", p_row["p-value"], r["p"], "R fisher.test", CLOSED, "relative")
        cond = [m for m in found if "conditional" in m]
        if cond:
            checks.compare(ds, f'{word_name} (p_value_by="exact")', "conditional odds ratio", found[cond[0]]["estimate"],
                           r["or"], "R fisher.test", 1e-6, "relative")
        for opts, alternative in [({"sides": 1, "higher": "Female"}, "greater"), ({"sides": 1, "lower": "Female"}, "less"),
                                  ({"alternative": "less"}, "less")]:
            found = risk_rows(call(**opts))
            r = r_run(f'list(p = prop.test(c({x1}, {x2}), c({n1}, {n2}), alternative = "{alternative}")$p.value)')
            checks.compare(ds, f"{word_name} (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")", "p-value",
                           found["risk difference"]["p-value"], r["p"], "R prop.test (Yates)", CLOSED, "relative")
    # equivalence and non-inferiority, in percentage points (Wald, unpooled)
    p1, p2 = x1 / n1, x2 / n2
    se = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    d = p1 - p2
    p_lo, p_hi = stats.norm.sf((d + 0.05) / se), stats.norm.cdf((d - 0.05) / se)
    found = call(same_within=5)
    row = next(r_ for r_ in found if "equivalence" in r_["measure"])
    checks.compare(ds, f"{word_name} (same_within=)", "p-value against the low margin", row["p-value against the low margin"],
                   p_lo, "two one-sided Wald z tests, SciPy", CLOSED, "relative")
    checks.compare(ds, f"{word_name} (same_within=)", "p-value against the high margin", row["p-value against the high margin"],
                   p_hi, "two one-sided Wald z tests, SciPy", CLOSED, "relative")
    for opt, ref in [("not_lower", p_lo), ("not_higher", p_hi)]:
        found = call(same_within=5, **{opt: "Female"})
        row = [r_ for r_ in found if r_["p-value"] is not None and r_["p-value"] == r_["p-value"]][-1]
        checks.compare(ds, f'{word_name} (same_within=, {opt}="Female")', "p-value", row["p-value"], ref,
                       "one-sided Wald z test, SciPy", CLOSED, "relative")


found = keyed(table_of(berkeley.compare_risks, "admit", for_each_group_in="sex", against="Male", value="Admitted",
                 allowing_for="department"))
print(found)
strata = []
for dept, d in applicants.groupby("department"):
    a = ((d["sex"] == "Female") & (d["admit"] == "Admitted")).sum(); b = ((d["sex"] == "Female") & (d["admit"] != "Admitted")).sum()
    c = ((d["sex"] == "Male") & (d["admit"] == "Admitted")).sum(); dd = ((d["sex"] == "Male") & (d["admit"] != "Admitted")).sum()
    strata.append((a, b, c, dd))
S = np.array(strata, dtype=float)
a, b, c, dd = S.T
n = S.sum(axis=1)
n1s, n0s = a + b, c + dd
rr_mh = (a * n0s / n).sum() / (c * n1s / n).sum()
var_log = (((n1s * n0s * (a + c) - a * c * n) / n ** 2).sum()) / ((a * n0s / n).sum() * (c * n1s / n).sum())
rr_rows = [m for m in found if m.startswith("risk ratio allowing")]
or_rows = [m for m in found if m.startswith("odds ratio allowing")]
z = stats.norm.ppf(0.975)
ds_mh = "UCB admissions by department (R mantelhaen.test example)"
checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel risk ratio", found[rr_rows[0]]["estimate"],
               rr_mh, "Mantel-Haenszel formula, NumPy", CLOSED, "relative")
checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel risk ratio 95 % CI low",
               found[rr_rows[0]]["estimate 95 % CI low"], rr_mh * np.exp(-z * np.sqrt(var_log)),
               "Greenland & Robins (1985) variance, NumPy", CLOSED, "relative")
checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel risk ratio 95 % CI high",
               found[rr_rows[0]]["estimate 95 % CI high"], rr_mh * np.exp(z * np.sqrt(var_log)),
               "Greenland & Robins (1985) variance, NumPy", CLOSED, "relative")
if HAS_R:
    r = r_run("""t <- UCBAdmissions[c("Admitted", "Rejected"), c("Female", "Male"), ]
        t <- aperm(t, c(2, 1, 3))
        m <- mantelhaen.test(t, correct = FALSE); mc <- mantelhaen.test(t)
        bd <- DescTools::BreslowDayTest(t)
        list(or = unname(m$estimate), low = m$conf.int[1], high = m$conf.int[2], p = m$p.value, pc = mc$p.value,
             bd = unname(bd$statistic), bdp = bd$p.value)""")
    row = found[or_rows[0]]
    checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel odds ratio", row["estimate"], r["or"],
                   "R mantelhaen.test", CLOSED, "relative")
    checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel odds ratio 95 % CI low",
                   row["estimate 95 % CI low"], r["low"], "R mantelhaen.test", CLOSED, "relative")
    checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Mantel-Haenszel odds ratio 95 % CI high",
                   row["estimate 95 % CI high"], r["high"], "R mantelhaen.test", CLOSED, "relative")
    checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Cochran-Mantel-Haenszel p-value", row["p-value"],
                   r["pc"], "R mantelhaen.test (continuity correction, its default)", CLOSED, "relative")
    bd_row = [m for m in found if "Breslow-Day" in m][0]
    checks.compare(ds_mh, 'compare_risks (allowing_for="department")', "Breslow-Day p-value", found[bd_row]["p-value"],
                   r["bdp"], "R DescTools BreslowDayTest", CLOSED, "relative")


from checks import load_csv

lalonde = load_csv("lalonde.csv")
COLUMNS = ["age", "educ", "race", "married", "nodegree", "re74", "re75"]
jobs = sdata.table(lalonde)
ds = "lalonde (Dehejia & Wahba 1999)"
for method in ("weighting", "matching"):
    extra = {"effect_for": 1} if method == "matching" else {}
    base = table_of(jobs.compare_balanced_groups_of, "re78", for_each_group_in="treat", allowing_for=COLUMNS,
                    test_method=method, **extra)[0]
    print(base)
    if method == "weighting":
        other = table_of(jobs.compare_balanced_groups_of, "re78", for_each_group_in="treat", allowing_for=COLUMNS,
                         test_method=method, against=1)[0]
        checks.compare(ds, f'compare_balanced_groups_of (test_method="{method}", against=1)', "effect, the other way round",
                       other["estimate"], -base["estimate"], "minus the effect against 0", ITERATIVE, "relative")
        checks.compare(ds, f'compare_balanced_groups_of (test_method="{method}", against=1)', "standard error",
                       other["standard error"], base["standard error"], "the same standard error", ITERATIVE, "relative")
    ninety = table_of(jobs.compare_balanced_groups_of, "re78", for_each_group_in="treat", allowing_for=COLUMNS,
                      test_method=method, interval=90, **extra)[0]
    lo_c = [c for c in ninety if "CI low" in c][0]
    hi_c = [c for c in ninety if "CI high" in c][0]
    lo95 = [c for c in base if "CI low" in c][0]
    q = (base["estimate"] - base[lo95]) / base["standard error"]
    factor = stats.norm.ppf(0.95) / stats.norm.ppf(0.975) if abs(q - stats.norm.ppf(0.975)) < 1e-6 else None
    if factor is None:   # a t interval: find its degrees of freedom from the 95 % interval
        from scipy.optimize import brentq
        df_ = brentq(lambda v: stats.t.ppf(0.975, v) - q, 1, 1e7)
        half = stats.t.ppf(0.95, df_) * base["standard error"]
        ref_name = f"estimate ± t(0.95, {df_:.0f} df) SE"
    else:
        half = stats.norm.ppf(0.95) * base["standard error"]
        ref_name = "estimate ± z(0.95) SE"
    checks.compare(ds, f'compare_balanced_groups_of (test_method="{method}", interval=90)', lo_c, ninety[lo_c],
                   base["estimate"] - half, ref_name, 1e-9, "relative")
    checks.compare(ds, f'compare_balanced_groups_of (test_method="{method}", interval=90)', hi_c, ninety[hi_c],
                   base["estimate"] + half, ref_name, 1e-9, "relative")


checks.data("cities", cities)
checks.data("lung", lung)
checks.data("applicants", applicants)
checks.save()
