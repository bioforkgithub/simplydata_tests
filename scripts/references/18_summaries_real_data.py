"""Reference values for notebooks/18_summaries_real_data.ipynb, saved to references/18_summaries_real_data.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/18_summaries_real_data/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/18_summaries_real_data.py      (from the repository root; needs R with the packages listed)
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
from statsmodels.stats.proportion import multinomial_proportions_confint, proportion_confint
from statsmodels.stats.weightstats import DescrStatsW

import realdata as rd
from checks import Checks, quiet, r_available, r_run, read_nist, read_numbers, read_rows, read_rows_by, versions

versions("pandas", "numpy", "scipy", "statsmodels", "rpy2")
checks = Recorder("18_summaries_real_data")
HAS_R = r_available()
EXACT = 1e-12      # sums, counts, extremes, order statistics
CLOSED = 1e-10     # closed formulas
if HAS_R:
    print(rd.r_version_of("DescTools", "e1071", "psych", "WRS2", "Hmisc", "survey", "boot"))


# The two-column numbers of find_correlation_between (the p-value, the bootstrap interval) come back in the one-row table
# of show="result table" (simplydata 3fd3b85 and later), so their references are kept at full precision. With
# allowing_for= the same table gives the partial correlation's p-value (simplydata a5a8ff0 and later).


def record(ds, method, quantity, reference, name, tolerance=None, how="relative"):
    """A reference kept at full precision, for a number the notebook reads from show="result table"."""
    checks.compare(ds, method, quantity, None, reference, name, CLOSED if tolerance is None else tolerance, how)
NIST_DATA = {}


def rows(result):
    """A simplydata result table as {column: [values]}, read row by row with simplydata's table.get_row
    (checks.read_rows); a number stays a number. Only the reference side uses pandas."""
    if not hasattr(result, "get_row"):
        return result
    found = read_rows(result)
    return {c: [row[c] for row in found] for c in found[0]}


def col(found, k):
    """The k-th column (from 0; -1 the last) of a result read by rows()."""
    return found[list(found)[k]]


NIST_UNIVARIATE = ["PiDigits", "Lottery", "Lew", "Mavro", "Michelso", "NumAcc1", "NumAcc2", "NumAcc3", "NumAcc4"]
for name in NIST_UNIVARIATE:
    problem = read_nist(name)
    values = problem.data["y"].to_numpy()
    NIST_DATA[name] = problem.data[["y"]]
    nums = sdata.numbers(values.tolist())
    ds = f"NIST {name}"
    checks.compare(ds, "nums.find_average", "mean", quiet(nums.find_average), problem.certified["mean"],
                   "NIST certified", 9, "lre")
    checks.compare(ds, "nums.find_standard_deviation", "standard deviation", quiet(nums.find_standard_deviation),
                   problem.certified["standard deviation"], "NIST certified", 9, "lre")
    if HAS_R:
        r = r_run("""q <- quantile(y, c(0.1, 0.5, 0.9), type = 7)
            tb <- table(y); top <- as.numeric(names(tb)[which.max(tb)])
            list(total = sum(y), biggest = max(y), smallest = min(y), p10 = q[[1]], p50 = q[[2]], p90 = q[[3]],
                 kinds = length(tb), top = top, topcount = max(tb))""", {"y": values})
        checks.compare(ds, "nums.find_total", "total", quiet(nums.find_total), r["total"], "R sum", EXACT, "relative")
        checks.compare(ds, "nums.find_biggest", "biggest", quiet(nums.find_biggest), r["biggest"], "R max", EXACT,
                       "relative")
        checks.compare(ds, "nums.find_smallest", "smallest", quiet(nums.find_smallest), r["smallest"], "R min",
                       EXACT, "relative")
        for percent, key in [(10, "p10"), (50, "p50"), (90, "p90")]:
            checks.compare(ds, "nums.find_percentile (percent=)", f"percentile {percent}",
                           float(quiet(nums.find_percentile, percent)), r[key], "R quantile(type = 7)", EXACT,
                           "relative")
        counted = rows(quiet(nums.count_each))
        checks.compare(ds, "nums.count_each", "different values", len(col(counted, 0)), r["kinds"], "R table", 0)
        checks.compare(ds, "nums.count_each", "biggest count", col(counted, 1)[0], r["topcount"], "R table", 0)

for _n, _f in NIST_DATA.items():
    checks.data(f"nist_{_n}", _f)


airquality = rd.r_dataset("airquality")
print(airquality.head())
ozone = sdata.numbers(airquality["Ozone"].tolist())
ozone.count_each()

ds = "airquality Ozone"
if HAS_R:
    r = r_run("""tb <- table(Ozone, useNA = "ifany")
        list(mean = mean(Ozone, na.rm = TRUE), sd = sd(Ozone, na.rm = TRUE), total = sum(Ozone, na.rm = TRUE),
             biggest = max(Ozone, na.rm = TRUE), smallest = min(Ozone, na.rm = TRUE),
             p25 = quantile(Ozone, 0.25, na.rm = TRUE)[[1]], kinds = length(tb), empty = sum(is.na(Ozone)),
             empty_percent = 100 * mean(is.na(Ozone)))""", airquality[["Ozone"]])
    checks.compare(ds, "nums.find_average", "mean", quiet(ozone.find_average), r["mean"], "R mean", EXACT, "relative")
    checks.compare(ds, "nums.find_standard_deviation", "standard deviation", quiet(ozone.find_standard_deviation),
                   r["sd"], "R sd", CLOSED, "relative")
    checks.compare(ds, "nums.find_total", "total", quiet(ozone.find_total), r["total"], "R sum", EXACT, "relative")
    checks.compare(ds, "nums.find_biggest", "biggest", quiet(ozone.find_biggest), r["biggest"], "R max", 0)
    checks.compare(ds, "nums.find_smallest", "smallest", quiet(ozone.find_smallest), r["smallest"], "R min", 0)
    checks.compare(ds, "nums.find_percentile (percent=)", "percentile 25", float(quiet(ozone.find_percentile, 25)),
                   r["p25"], "R quantile(type = 7)", EXACT, "relative")
    counted = rows(quiet(ozone.count_each))
    checks.compare(ds, "nums.count_each", "different values (empty counted)", len(col(counted, 0)), r["kinds"],
                   'R table(useNA = "ifany")', 0)
    empty_at = [i for i, v in enumerate(col(counted, 0)) if v is None or v != v][0]
    checks.compare(ds, "nums.count_each", "empty cells", col(counted, 1)[empty_at], r["empty"], "R sum(is.na)", 0)
    checks.compare(ds, "nums.count_each", "empty cells percent (1 decimal)", col(counted, 2)[empty_at],
                   r["empty_percent"], "R, rounded to 1 decimal as printed", 0.05)


air = airquality.copy()
air["date"] = pd.to_datetime(dict(year=1973, month=air["Month"].astype(int), day=air["Day"].astype(int)))
table = sdata.table(air)
table.find_biggest_of("Ozone", for_each_group_in="Month")

R_BY_MONTH = """g <- Month
f <- function(v, fun) as.numeric(tapply(v, g, fun))
mode1 <- function(v) { v <- v[!is.na(v)]; tb <- table(v); as.numeric(names(tb)[which.max(tb)]) }
modecount <- function(v) { v <- v[!is.na(v)]; max(table(v)) }
modeties <- function(v) { v <- v[!is.na(v)]; tb <- table(v); sum(tb == max(tb)) }
list(biggest = f(Ozone, function(v) max(v, na.rm = TRUE)), smallest = f(Ozone, function(v) min(v, na.rm = TRUE)),
     total = f(Ozone, function(v) sum(v, na.rm = TRUE)), rows = as.numeric(table(g)),
     kinds = f(Temp, function(v) length(unique(v[!is.na(v)]))),
     mode = f(Temp, mode1), modecount = f(Temp, modecount), modeties = f(Temp, modeties),
     modepercent = f(Temp, function(v) 100 * modecount(v) / sum(!is.na(v))),
     mean = f(Ozone, function(v) mean(v, na.rm = TRUE)), sd = f(Ozone, function(v) sd(v, na.rm = TRUE)),
     median = f(Ozone, function(v) median(v, na.rm = TRUE)),
     allbiggest = max(Ozone, na.rm = TRUE), allsmallest = min(Ozone, na.rm = TRUE), alltotal = sum(Ozone, na.rm = TRUE),
     allrows = length(Ozone), allkinds = length(unique(Temp)), allmode = mode1(Temp))"""
if HAS_R:
    r = r_run(R_BY_MONTH, air[["Month", "Ozone", "Temp"]])
    months = [5, 6, 7, 8, 9]
    for every in (None, "month"):
        group = "date" if every else "Month"
        how = f'every="{every}"' if every else 'for_each_group_in="Month"'
        ref = "R tapply by cut(date, \"month\")" if every else "R tapply"

        def by_group(word, *args, column):
            found = rows(quiet(word, *args, for_each_group_in=group, every=every))
            return np.asarray(found[column], dtype=float)

        pairs = [("find_biggest_of", table.find_biggest_of, ("Ozone",), "biggest"),
                 ("find_smallest_of", table.find_smallest_of, ("Ozone",), "smallest"),
                 ("find_total_of", table.find_total_of, ("Ozone",), "total"),
                 ("find_average_of", table.find_average_of, ("Ozone",), "mean"),
                 ("find_median_of", table.find_median_of, ("Ozone",), "median"),
                 ("find_standard_deviation_of", table.find_standard_deviation_of, ("Ozone",), "sd")]
        for label, word, args, key in pairs:
            found = rows(quiet(word, *args, for_each_group_in=group, every=every))
            values = np.asarray(col(found, -1), dtype=float)
            tol = CLOSED if key == "sd" else EXACT
            for month, ours, ref_value in zip(months, values, r[key]):
                checks.compare(ds := "airquality", f"{label} ({how})", f"Ozone {key}, month {month}", ours, ref_value,
                               ref, tol, "relative")
        found = rows(quiet(table.count_rows, for_each_group_in=group, every=every))
        for month, ours, ref_value in zip(months, col(found, -1), r["rows"]):
            checks.compare("airquality", f"count_rows ({how})", f"rows, month {month}", ours, ref_value, "R table", 0)
        found = rows(quiet(table.count_different_values_in, "Temp", for_each_group_in=group, every=every))
        for month, ours, ref_value in zip(months, col(found, -1), r["kinds"]):
            checks.compare("airquality", f"count_different_values_in ({how})", f"Temp values, month {month}", ours,
                           ref_value, "R length(unique())", 0)
        found = rows(quiet(table.find_most_common_of, "Temp", for_each_group_in=group, every=every))
        for i, month in enumerate(months):
            line = [col(found, k)[i] for k in range(5)]
            checks.compare("airquality", f"find_most_common_of ({how})", f"Temp mode, month {month}", line[1],
                           r["mode"][i], "R table, smallest of the tied", 0)
            checks.compare("airquality", f"find_most_common_of ({how})", f"Temp mode count, month {month}",
                           line[2], r["modecount"][i], "R table", 0)
            checks.compare("airquality", f"find_most_common_of ({how})", f"Temp mode percent, month {month}",
                           line[3], r["modepercent"][i], "R table", EXACT, "relative")
            checks.compare("airquality", f"find_most_common_of ({how})", f"Temp values tied at the top, month {month}",
                           line[4], r["modeties"][i], "R table", 0)
    # without groups
    checks.compare("airquality", "find_biggest_of", "Ozone biggest", quiet(table.find_biggest_of, "Ozone"),
                   r["allbiggest"], "R max", 0)
    checks.compare("airquality", "find_smallest_of", "Ozone smallest", quiet(table.find_smallest_of, "Ozone"),
                   r["allsmallest"], "R min", 0)
    checks.compare("airquality", "find_total_of", "Ozone total", quiet(table.find_total_of, "Ozone"), r["alltotal"],
                   "R sum", EXACT, "relative")
    checks.compare("airquality", "count_rows", "rows", quiet(table.count_rows), r["allrows"], "R nrow", 0)
    checks.compare("airquality", "count_different_values_in", "Temp values",
                   quiet(table.count_different_values_in, "Temp"), r["allkinds"], "R length(unique())", 0)
    checks.compare("airquality", "find_most_common_of", "Temp mode", quiet(table.find_most_common_of, "Temp"),
                   r["allmode"], "R table, smallest of the tied", 0)


table.find_median_of("Ozone", for_each_group_in="Month", interval=95)

if HAS_R:
    r = r_run("""g <- Month
        tt <- lapply(split(Ozone, g), function(v) t.test(v[!is.na(v)])$conf.int)
        med <- lapply(split(Ozone, g), function(v) DescTools::MedianCI(v, method = "exact", na.rm = TRUE))
        q <- function(p) lapply(split(Ozone, g), function(v) DescTools::QuantileCI(v, probs = p, method = "exact", na.rm = TRUE))
        sdci <- lapply(split(Ozone, g), function(v) { v <- v[!is.na(v)]; n <- length(v)
             sqrt((n - 1) * var(v) / qchisq(c(0.975, 0.025), n - 1)) })
        list(tlow = sapply(tt, `[`, 1), thigh = sapply(tt, `[`, 2), mlow = sapply(med, `[`, 2), mhigh = sapply(med, `[`, 3),
             q25low = sapply(q(0.25), `[`, 2), q25high = sapply(q(0.25), `[`, 3),
             q75low = sapply(q(0.75), `[`, 2), q75high = sapply(q(0.75), `[`, 3),
             sdlow = sapply(sdci, `[`, 1), sdhigh = sapply(sdci, `[`, 2))""", air[["Month", "Ozone"]])
    found = rows(quiet(table.find_average_of, "Ozone", for_each_group_in="Month", interval=95))
    lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
    for i, month in enumerate(months):
        checks.compare("airquality", "find_average_of (interval=)", f"Ozone 95 % CI low, month {month}",
                       found[lo][i], r["tlow"][i], "R t.test", CLOSED, "relative")
        checks.compare("airquality", "find_average_of (interval=)", f"Ozone 95 % CI high, month {month}",
                       found[hi][i], r["thigh"][i], "R t.test", CLOSED, "relative")
    for word, label, keys in [(table.find_median_of, "find_median_of (interval=)", ("mlow", "mhigh")),
                              (table.find_standard_deviation_of, "find_standard_deviation_of (interval=)",
                               ("sdlow", "sdhigh"))]:
        found = rows(quiet(word, "Ozone", for_each_group_in="Month", interval=95))
        lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
        ref = "R DescTools MedianCI(exact)" if keys[0] == "mlow" else "chi-square formula in R"
        for i, month in enumerate(months):
            checks.compare("airquality", label, f"Ozone 95 % CI low, month {month}", found[lo][i],
                           r[keys[0]][i], ref, CLOSED, "relative")
            checks.compare("airquality", label, f"Ozone 95 % CI high, month {month}", found[hi][i],
                           r[keys[1]][i], ref, CLOSED, "relative")
    for percent, keys in [(25, ("q25low", "q25high")), (75, ("q75low", "q75high"))]:
        found = rows(quiet(table.find_percentile_of, "Ozone", percent, for_each_group_in="Month", interval=95))
        lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
        for i, month in enumerate(months):
            for col, key, end in [(lo, keys[0], "low"), (hi, keys[1], "high")]:
                checks.compare("airquality", "find_percentile_of (interval=)", f"Ozone p{percent} 95 % CI {end}, month {month}",
                               found[col][i], r[key][i], "R DescTools QuantileCI(exact)", EXACT, "relative")
    # the parent, run directly: SciPy quantile_test (for the whole column)
    x = air["Ozone"].dropna().to_numpy()
    found = rows(quiet(table.find_percentile_of, "Ozone", 90, interval=95))
    ci = stats.quantile_test(x, q=np.quantile(x, 0.9), p=0.9).confidence_interval(0.95)
    lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
    checks.compare("airquality", "find_percentile_of (interval=)", "Ozone p90 95 % CI low", found[lo][0], ci.low,
                   "SciPy quantile_test", EXACT, "relative")
    checks.compare("airquality", "find_percentile_of (interval=)", "Ozone p90 95 % CI high", found[hi][0], ci.high,
                   "SciPy quantile_test", EXACT, "relative")


RULES = {1: "inverted cdf", 2: "averaged inverted cdf", 3: "closest observation", 4: "interpolated inverted cdf",
         5: "hazen", 6: "weibull", 7: "linear", 8: "median unbiased", 9: "normal unbiased"}
PERCENTS = [5, 25, 50, 80, 99]
if HAS_R:
    for column in ("Ozone", "Wind"):
        r = r_run("x <- v[!is.na(v)]; list(" + ", ".join(
            f"t{t} = as.numeric(quantile(x, c({', '.join(str(p / 100) for p in PERCENTS)}), type = {t}))"
            for t in RULES) + ")", {"v": air[column]})
        for t, name in RULES.items():
            for p, ref in zip(PERCENTS, r[f"t{t}"]):
                ours = quiet(table.find_percentile_of, column, p, test_method=name)
                checks.compare(f"airquality {column}", f'find_percentile_of (test_method="{name}")', f"percentile {p}",
                               ours, ref, f"R quantile(type = {t})", EXACT, "relative")
        ours = quiet(table.find_percentile_of, column, 90, test_method="inverted_cdf")
        r1 = r_run("list(q = quantile(v[!is.na(v)], 0.9, type = 1)[[1]])", {"v": air[column]})
        checks.compare(f"airquality {column}", 'find_percentile_of (test_method="inverted_cdf")', "percentile 90", ours,
                       r1["q"], "R quantile(type = 1)", EXACT, "relative")


counts = air["Temp"].value_counts().sort_index()
frequency_rows = (pd.DataFrame({"temperature": counts.index.to_numpy(dtype=float),
                               "days": counts.to_numpy(dtype=float)}))
frequency = sdata.table(pd.DataFrame({"temperature": counts.index.to_numpy(dtype=float),
                                      "days": counts.to_numpy(dtype=float)}))
if HAS_R:
    r = r_run("""list(p10 = quantile(Temp, 0.1, type = 7)[[1]], p90 = quantile(Temp, 0.9, type = 7)[[1]],
        med = median(Temp), mean = mean(Temp), sd = sd(Temp))""", air[["Temp"]])
    ds = "airquality Temp as a frequency table"
    for p, key in [(10, "p10"), (90, "p90")]:
        checks.compare(ds, 'find_percentile_of (weights=, weight_kind="frequency")', f"percentile {p}",
                       quiet(frequency.find_percentile_of, "temperature", p, weights="days", weight_kind="frequency"),
                       r[key], "R quantile of the 153 days", EXACT, "relative")
    checks.compare(ds, 'find_median_of (weights=, weight_kind="frequency")', "median",
                   quiet(frequency.find_median_of, "temperature", weights="days", weight_kind="frequency"), r["med"],
                   "R median of the 153 days", EXACT, "relative")
    checks.compare(ds, "find_average_of (weights=)", "average",
                   quiet(frequency.find_average_of, "temperature", weights="days"), r["mean"],
                   "R mean of the 153 days", EXACT, "relative")
    checks.compare(ds, 'find_standard_deviation_of (weights=, weight_kind="frequency")', "standard deviation",
                   quiet(frequency.find_standard_deviation_of, "temperature", weights="days", weight_kind="frequency"),
                   r["sd"], "R sd of the 153 days", CLOSED, "relative")

apistrat = rd.r_dataset("apistrat", "survey", data_set="api")
schools_rows = apistrat[["api00", "pw", "stype"]]
schools = sdata.table(schools_rows)
if HAS_R:
    r = r_run("""library(survey)
        d <- svydesign(id = ~1, weights = ~pw, data = data.frame(api00 = api00, pw = pw))
        q <- svyquantile(~api00, d, c(0.25, 0.5, 0.9), qrule = "math", ci = FALSE)$api00
        list(q = as.numeric(q), sd = sqrt(Hmisc::wtd.var(api00, pw, normwt = TRUE)), wmean = weighted.mean(api00, pw))""",
              apistrat[["api00", "pw"]])
    ds = "API stratified sample (apistrat), api00"
    for p, ref in zip((25, 50, 90), r["q"]):
        ours = quiet(schools.find_percentile_of, "api00", p, weights="pw", weight_kind="analytic")
        checks.compare(ds, 'find_percentile_of (weights=, weight_kind="analytic")', f"percentile {p}", ours, ref,
                       'R survey svyquantile(qrule = "math")', EXACT, "relative")
        parent = np.percentile(apistrat["api00"], p, weights=apistrat["pw"], method="inverted_cdf")
        checks.compare(ds, 'find_percentile_of (weights=, weight_kind="analytic")', f"percentile {p}", ours, parent,
                       "NumPy percentile(weights=)", EXACT, "relative")
    checks.compare(ds, 'find_median_of (weights=, weight_kind="analytic")', "median",
                   quiet(schools.find_median_of, "api00", weights="pw", weight_kind="analytic"), r["q"][1],
                   'R survey svyquantile(qrule = "math")', EXACT, "relative")
    checks.compare(ds, 'find_standard_deviation_of (weights=, weight_kind="analytic")', "standard deviation",
                   quiet(schools.find_standard_deviation_of, "api00", weights="pw", weight_kind="analytic"), r["sd"],
                   "R Hmisc wtd.var(normwt = TRUE)", CLOSED, "relative")
    checks.compare(ds, "find_average_of (weights=)", "average", quiet(schools.find_average_of, "api00", weights="pw"),
                   r["wmean"], "R weighted.mean", EXACT, "relative")
    checks.compare(ds, "find_average_of (weights=)", "average", quiet(schools.find_average_of, "api00", weights="pw"),
                   DescrStatsW(apistrat["api00"], weights=apistrat["pw"]).mean, "statsmodels DescrStatsW", EXACT,
                   "relative")


SHOW = ["count", "average", "median", "total", "smallest", "biggest", "standard deviation", "variance",
        "standard error", "confidence interval", "range", "interquartile range", "empty cells", "geometric mean",
        "skewness", "kurtosis", "coefficient of variation", "median absolute deviation", "trimmed average",
        "trimmed average 20", "most common", "percentile 90", "harmonic mean", "winsorised average",
        "winsorised average 10"]
summary = rows(quiet(table.find_summary_of, "Ozone", for_each_group_in="Month", show=SHOW))
print(list(summary))
summary

R_SUMMARY = """one <- function(v) { e <- sum(is.na(v)); v <- v[!is.na(v)]; n <- length(v); tb <- table(v)
    ci <- t.test(v)$conf.int
    c(count = n, average = mean(v), median = median(v), total = sum(v), smallest = min(v), biggest = max(v),
      sd = sd(v), var = var(v), se = sd(v) / sqrt(n), cilow = ci[1], cihigh = ci[2], range = max(v) - min(v),
      iqr = IQR(v, type = 7), empty = e, gmean = psych::geometric.mean(v), skew = e1071::skewness(v, type = 2),
      kurt = e1071::kurtosis(v, type = 2), cv = sd(v) / mean(v), mad = mad(v), trim10 = mean(v, trim = 0.1),
      trim20 = mean(v, trim = 0.2), mode = as.numeric(names(tb)[which.max(tb)]), p90 = quantile(v, 0.9, type = 7)[[1]],
      hmean = psych::harmonic.mean(v), win20 = WRS2::winmean(v, tr = 0.2), win10 = WRS2::winmean(v, tr = 0.1)) }
m <- sapply(split(Ozone, Month), one)
setNames(lapply(rownames(m), function(k) as.numeric(m[k, ])), rownames(m))"""
KEYS = {"count": "count", "average": "average", "median": "median", "total": "total", "smallest": "smallest",
        "biggest": "biggest", "standard deviation": "sd", "variance": "var", "standard error": "se",
        "range": "range", "interquartile range": "iqr", "empty cells": "empty", "geometric mean": "gmean",
        "skewness": "skew", "kurtosis": "kurt", "coefficient of variation": "cv", "median absolute deviation": "mad",
        "trimmed average 10 %": "trim10", "trimmed average 20 %": "trim20", "most common": "mode",
        "percentile 90": "p90", "harmonic mean": "hmean", "winsorised average 20 %": "win20",
        "winsorised average 10 %": "win10"}
if HAS_R:
    r = r_run(R_SUMMARY, air[["Month", "Ozone"]])
    for every in (None, "month"):
        found = rows(quiet(table.find_summary_of, "Ozone", for_each_group_in="date" if every else "Month",
                           every=every, show=SHOW))
        how = 'every="month"' if every else 'for_each_group_in="Month"'
        for i, month in enumerate(months):
            for measure, key in KEYS.items():
                column = [c for c in list(found) if c == measure]
                if not column:
                    checks.skip("airquality", f"find_summary_of ({how})", f"Ozone {measure}, month {month}", "R",
                                f"no column {measure!r} in the result")
                    continue
                ref = r[key][i]
                tol = EXACT if key in ("count", "total", "smallest", "biggest", "median", "range", "iqr", "empty",
                                       "mode", "p90", "average") else CLOSED
                checks.compare("airquality", f"find_summary_of ({how}, show=)", f"Ozone {measure}, month {month}",
                               found[measure][i], ref, "R (e1071, psych, WRS2, base)", tol, "relative")
            lo = [c for c in list(found) if "95 % CI low" in c]
            hi = [c for c in list(found) if "95 % CI high" in c]
            if lo and hi:
                checks.compare("airquality", f"find_summary_of ({how}, show=)", f"Ozone 95 % CI low, month {month}",
                               found[lo[0]][i], r["cilow"][i], "R t.test", CLOSED, "relative")
                checks.compare("airquality", f"find_summary_of ({how}, show=)", f"Ozone 95 % CI high, month {month}",
                               found[hi[0]][i], r["cihigh"][i], "R t.test", CLOSED, "relative")


table.find_interval_of("Ozone", "median")

x = air["Ozone"].dropna().to_numpy()
METHODS = {"bootstrap": "BCa", "percentile bootstrap": "percentile", "reverse percentile bootstrap": "basic"}
SUMMARIES = {"average": np.mean, "median": np.median, "standard deviation": lambda v: np.std(v, ddof=1)}
for summary_name, function in SUMMARIES.items():
    for method, scipy_method in METHODS.items():
        found = rows(quiet(table.find_interval_of, "Ozone", summary_name, test_method=method, resamples=9999,
                           with_seed=1))
        lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
        res = stats.bootstrap((x,), lambda v, axis=-1: np.apply_along_axis(function, axis, v), n_resamples=9999,
                              method=scipy_method, confidence_level=0.95, rng=np.random.default_rng(1))
        label = f'find_interval_of ("{summary_name}", test_method="{method}")'
        checks.compare("airquality Ozone", label, "95 % CI low", found[lo][0], res.confidence_interval.low,
                       "SciPy bootstrap, same seed", EXACT, "relative")
        checks.compare("airquality Ozone", label, "95 % CI high", found[hi][0], res.confidence_interval.high,
                       "SciPy bootstrap, same seed", EXACT, "relative")

if HAS_R:
    R_TYPES = {"bootstrap": ("bca", 4, 5), "percentile bootstrap": ("percent", 4, 5),
               "reverse percentile bootstrap": ("basic", 4, 5)}
    r = r_run("""library(boot); set.seed(1)
        b <- boot(x, function(v, i) mean(v[i]), R = 200000)
        out <- list()
        for (t in c("bca", "perc", "basic")) { ci <- boot.ci(b, type = t)[[c(bca = "bca", perc = "percent", basic = "basic")[[t]]]]
            out[[paste0(t, "low")]] <- ci[4]; out[[paste0(t, "high")]] <- ci[5] }
        out""", {"x": x})
    for method, (scipy_method, key) in {"bootstrap": ("BCa", "bca"), "percentile bootstrap": ("percentile", "perc"),
                                        "reverse percentile bootstrap": ("basic", "basic")}.items():
        runs = np.array([stats.bootstrap((x,), np.mean, n_resamples=9999, method=scipy_method,
                                         rng=np.random.default_rng(100 + k)).confidence_interval for k in range(20)])
        mcse = runs.std(axis=0, ddof=1)
        found = rows(quiet(table.find_interval_of, "Ozone", "average", test_method=method, resamples=9999, with_seed=1))
        lo, hi = [c for c in list(found) if "CI low" in c][0], [c for c in list(found) if "CI high" in c][0]
        for end, col, k in (("low", lo, 0), ("high", hi, 1)):
            allowed = 4 * np.sqrt(mcse[k] ** 2 * (1 + 9999 / 200000))
            checks.compare("airquality Ozone", f'find_interval_of ("average", test_method="{method}")',
                           f"95 % CI {end}", found[col][0], r[f"{key}{end}"],
                           "R boot.ci, 200 000 resamples (4 MCSE)", allowed)


measured_rows = air.dropna(subset=["Ozone"]).assign(month=lambda d: "month " + d["Month"].astype(int).astype(str))
measured = sdata.table(measured_rows)
measured.count_each("month", interval=95)

if HAS_R:
    found = rows(quiet(measured.count_each, "month", interval=95))
    order = sorted(range(len(found["month"])), key=lambda i: found["month"][i])   # by month, as R's table
    found = {c: [v[i] for i in order] for c, v in found.items()}
    lo, hi = [c for c in list(found) if "low" in c][0], [c for c in list(found) if "high" in c][0]
    month_counts = air.dropna(subset=["Ozone"])["Month"].astype(int).value_counts().sort_index()
    r = r_run("""m <- DescTools::MultinomCI(n, conf.level = 0.95, method = "sisonglaz")
        list(low = 100 * m[, 2], high = 100 * m[, 3])""", {"n": month_counts.to_numpy(dtype=float)})
    parent = 100 * multinomial_proportions_confint(month_counts.to_numpy(), alpha=0.05, method="sison-glaz")
    for i, month in enumerate(month_counts.index):
        checks.compare("airquality, days with ozone", "count_each (interval=)", f"month {month} 95 % low",
                       found[lo][i], r["low"][i], 'R DescTools MultinomCI("sisonglaz")', 1e-6, "relative")
        checks.compare("airquality, days with ozone", "count_each (interval=)", f"month {month} 95 % high",
                       found[hi][i], r["high"][i], 'R DescTools MultinomCI("sisonglaz")', 1e-6, "relative")
        checks.compare("airquality, days with ozone", "count_each (interval=)", f"month {month} 95 % low",
                       found[lo][i], parent[i, 0], "statsmodels multinomial_proportions_confint", EXACT, "relative")


table.find_percent_of_rows_where("Ozone", is_above=80, for_each_group_in="Month")

if HAS_R:
    r = r_run("""res <- lapply(split(Ozone, Month), function(v) { v <- v[!is.na(v)]; k <- sum(v > 80); n <- length(v)
            b <- binom.test(k, n)$conf.int
            w <- DescTools::BinomCI(k, n, method = "wilson"); a <- DescTools::BinomCI(k, n, method = "agresti-coull")
            c(100 * k / n, 100 * b, 100 * w[2:3], 100 * a[2:3]) })
        m <- do.call(rbind, res)
        list(p = m[, 1], elow = m[, 2], ehigh = m[, 3], wlow = m[, 4], whigh = m[, 5], alow = m[, 6], ahigh = m[, 7])""",
              air[["Month", "Ozone"]])
    for every in (None, "month"):
        how = 'every="month"' if every else 'for_each_group_in="Month"'
        for interval_by, keys, ref in [("exact", ("elow", "ehigh"), "R binom.test"),
                                       ("wilson", ("wlow", "whigh"), 'R DescTools BinomCI("wilson")'),
                                       ("agresti-coull", ("alow", "ahigh"), 'R DescTools BinomCI("agresti-coull")')]:
            found = rows(quiet(table.find_percent_of_rows_where, "Ozone", is_above=80,
                               for_each_group_in="date" if every else "Month", every=every, interval_by=interval_by))
            lo, hi = [c for c in list(found) if "low" in c][0], [c for c in list(found) if "high" in c][0]
            pc = [c for c in list(found) if c not in (lo, hi) and "percent" in c][0]
            for i, month in enumerate(months):
                label = f'find_percent_of_rows_where ({how}, interval_by="{interval_by}")'
                if interval_by == "exact":
                    checks.compare("airquality", label, f"percent of days above 80, month {month}", found[pc][i],
                                   r["p"][i], "R", EXACT, "relative")
                checks.compare("airquality", label, f"95 % CI low, month {month}", found[lo][i], r[keys[0]][i],
                               ref, CLOSED, "absolute")
                checks.compare("airquality", label, f"95 % CI high, month {month}", found[hi][i], r[keys[1]][i],
                               ref, CLOSED, "absolute")


made = rd.r_example_objects("stats", "cor.test")
print("objects made by the help page's examples:", made)
print(rd.r_example_code("stats", "cor.test")[:600])
hw = r_run("list(x = x, y = y)")
hollander_rows = pd.DataFrame({"x": hw["x"], "y": hw["y"]})
hollander = sdata.table(hollander_rows)
hollander.find_correlation_between("x", "y")

TEST = ["pearson", "spearman", "kendall"]
SIDES = [("two.sided", {}), ("greater", {"sides": 1, "higher": "positive"}), ("less", {"sides": 1, "higher": "negative"}),
         ("greater", {"alternative": "greater"})]


def check_correlations(name, t, data, a, b, exact_pvalue=True):
    """r (given back, full precision) and the p-value (printed in the sentence: 3 or 4 decimals) against R."""
    for method in TEST:
        for alternative, opts in SIDES:
            r = r_run(f"""ct <- suppressWarnings(cor.test(a, b, method = "{method}", alternative = "{alternative}",
                          exact = {'NULL' if exact_pvalue and method != "spearman" else 'FALSE'}))
                list(r = unname(ct$estimate), p = ct$p.value)""", {"a": data[a], "b": data[b]})
            found, text = rd.run_printed(t.find_correlation_between, a, b, test_method=method, **opts)
            option_text = ", ".join(f"{k}={v!r}" for k, v in opts.items())
            label = f'find_correlation_between (test_method="{method}"{", " + option_text if option_text else ""})'
            if not opts:
                checks.compare(name, label, "r", found, r["r"], "R cor.test", CLOSED, "relative")
            if method == "spearman" and exact_pvalue:
                # without ties and below 1290 rows R cor.test's default is the exact p-value (AS 89), which
                # simplydata takes from R; exact = FALSE (the t approximation) is another definition, not compared
                exact = r_run(f'list(p = cor.test(a, b, method = "spearman", alternative = "{alternative}")$p.value)',
                              {"a": data[a], "b": data[b]})
                record(name, label, "p-value", exact["p"], "R cor.test (default: exact)")
            else:
                record(name, label, "p-value", r["p"],
                       "R cor.test" + (" (exact = FALSE)" if method == "spearman" else ""))


if HAS_R:
    hw_frame = pd.DataFrame({"x": hw["x"], "y": hw["y"]})
    check_correlations("Hollander & Wolfe (R cor.test example)", hollander, hw_frame, "x", "y")
    mtcars = rd.r_dataset("mtcars")
    cars = sdata.table(mtcars)
    check_correlations("mtcars mpg, wt", cars, mtcars, "mpg", "wt", exact_pvalue=False)


if HAS_R:
    for method in TEST:
        pairs = rows(quiet(cars.find_correlation_between, "mpg", "wt", "hp", "qsec", test_method=method,
                           show="each pair"))
        if method == TEST[0]:
            print(pairs)
        for line in read_rows(quiet(cars.find_correlation_between, "mpg", "wt", "hp", "qsec", test_method=method,
                                    show="each pair")):
            a, b = list(line.values())[:2]
            r = r_run(f"""ct <- suppressWarnings(cor.test(a, b, method = "{method}", exact = FALSE))
                list(r = unname(ct$estimate), p = ct$p.value)""", {"a": mtcars[a], "b": mtcars[b]})
            label = f'find_correlation_between (test_method="{method}", show="each pair")'
            checks.compare("mtcars", label, f"r {a}, {b}", line["r"], r["r"], "R cor.test", CLOSED, "relative")
            checks.compare("mtcars", label, f"p-value {a}, {b}", line["p-value"], r["p"], "R cor.test", CLOSED,
                           "relative")


if HAS_R:
    for name, data, a, b in [("mtcars cyl, gear", mtcars, "cyl", "gear"), ("mtcars carb, am", mtcars, "carb", "am"),
                             ("airquality Month, Temp", air, "Month", "Temp")]:
        t = sdata.table(data)
        r = r_run("""C <- 0; D <- 0; Tx <- 0; n <- length(a)
            for (i in 1:(n - 1)) for (j in (i + 1):n) { s <- sign(a[i] - a[j]) * sign(b[i] - b[j])
                if (s > 0) C <- C + 1 else if (s < 0) D <- D + 1
                if (a[i] == a[j]) Tx <- Tx + 1 }
            pairs <- n * (n - 1) / 2
            list(gamma = (C - D) / (C + D), somers = (C - D) / (pairs - Tx),
                 dt_gamma = DescTools::GoodmanKruskalGamma(a, b), dt_somers = DescTools::SomersDelta(a, b, direction = "column"))""",
                  {"a": data[a], "b": data[b]})
        for method, key in [("gamma", "gamma"), ("somers", "somers")]:
            found = quiet(t.find_correlation_between, a, b, test_method=method)
            checks.compare(name, f'find_correlation_between (test_method="{method}")', "estimate", found, r[key],
                           "pairs counted in R", CLOSED, "relative")
            checks.compare(name, f'find_correlation_between (test_method="{method}")', "estimate", found,
                           r[f"dt_{key}"], "R DescTools", CLOSED, "relative")


if HAS_R:
    found, text = rd.run_printed(cars.find_correlation_between, "mpg", "wt", allowing_for="hp")
    r = r_run("""e1 <- resid(lm(mpg ~ hp)); e2 <- resid(lm(wt ~ hp)); r <- cor(e1, e2); n <- length(mpg)
        t <- r * sqrt((n - 3) / (1 - r^2))
        pc <- ppcor::pcor.test(mpg, wt, hp)
        list(r = r, p = 2 * pt(-abs(t), n - 3), ppcor_r = pc$estimate, ppcor_p = pc$p.value)""", mtcars)
    checks.compare("mtcars mpg, wt | hp", "find_correlation_between (allowing_for=)", "partial r", found, r["r"],
                   "residuals of R lm", CLOSED, "relative")
    checks.compare("mtcars mpg, wt | hp", "find_correlation_between (allowing_for=)", "p-value", None, r["p"],
                   "residuals of R lm, t on n - 3", CLOSED, "relative")
    checks.compare("mtcars mpg, wt | hp", "find_correlation_between (allowing_for=)", "p-value", None,
                   r["ppcor_p"], "R ppcor pcor.test", CLOSED, "relative")
    checks.compare("mtcars mpg, wt | hp", "find_correlation_between (allowing_for=)", "partial r", found,
                   r["ppcor_r"], "R ppcor pcor.test", CLOSED, "relative")

    p_table = quiet(cars.find_correlation_between, "mpg", "wt", "hp", show="p-value")
    first = p_table.get_column_names()[0]
    pvalues = read_rows_by(p_table, first)
    print(pvalues)
    for a, b in [("mpg", "wt"), ("mpg", "hp"), ("wt", "hp")]:
        r = r_run("list(p = cor.test(a, b)$p.value)", {"a": mtcars[a], "b": mtcars[b]})
        checks.compare("mtcars", 'find_correlation_between (show="p-value")', f"p-value {a}, {b}",
                       float(pvalues[a][b]), r["p"], "R cor.test", CLOSED, "relative")

def pearson_statistic(a, axis=-1):
    return stats.pearsonr(a, data_b, axis=axis).statistic


for name, frame, a, b in [("Hollander & Wolfe (R cor.test example)", pd.DataFrame(hw) if HAS_R else None, "x", "y")]:
    if frame is None:
        break
    t = sdata.table(frame)
    found, text = rd.run_printed(t.find_correlation_between, a, b, p_value_by="shuffling", shuffles=9999,
                                 with_seed=1)
    print(text)
    data_b = frame[b].to_numpy()
    # the same seed: simplydata shuffles the second column against the first ("the values of one column were shuffled
    # against the other"), so the reference shuffles b with a held fixed; shuffling a instead uses the random numbers
    # differently (0.1256 against 0.1260 here), and the million shuffles below are the check that does not depend on it
    data_a = frame[a].to_numpy()
    parent = stats.permutation_test((data_b,), lambda v, axis=-1: stats.pearsonr(np.broadcast_to(data_a, np.shape(v)), v,
                                                                                axis=axis).statistic,
                                    permutation_type="pairings", vectorized=True, n_resamples=9999,
                                    rng=np.random.default_rng(1))
    record(name, 'find_correlation_between (p_value_by="shuffling", shuffles=)', "p-value", parent.pvalue,
           "SciPy permutation_test, same seed", 1e-12)
    long = stats.permutation_test((frame[a].to_numpy(),), pearson_statistic, permutation_type="pairings",
                                  n_resamples=1_000_000, rng=np.random.default_rng(20261121))
    p = long.pvalue
    allowed = 4 * np.sqrt(p * (1 - p) / 9999 + p * (1 - p) / 1_000_000)
    record(name, 'find_correlation_between (p_value_by="shuffling", shuffles=)', "p-value", p,
           "1 000 000 shuffles, other seed (4 MCSE)", allowed, "absolute")
    found, text = rd.run_printed(t.find_correlation_between, a, b, interval_by="bootstrap", resamples=9999,
                                 with_seed=1)
    print(text)
    parent = stats.bootstrap((frame[a].to_numpy(), frame[b].to_numpy()),
                             lambda u, v, axis=-1: stats.pearsonr(u, v, axis=axis).statistic, paired=True,
                             vectorized=True, n_resamples=9999, method="BCa", rng=np.random.default_rng(1))
    record(name, 'find_correlation_between (interval_by="bootstrap", resamples=)', "95 % CI low",
           parent.confidence_interval.low, "SciPy bootstrap BCa, same seed")
    record(name, 'find_correlation_between (interval_by="bootstrap", resamples=)', "95 % CI high",
           parent.confidence_interval.high, "SciPy bootstrap BCa, same seed")


checks.data("airquality", airquality)
checks.data("air", air)
checks.data("frequency_rows", frequency_rows)
checks.data("schools_rows", schools_rows)
checks.data("measured_rows", measured_rows)
checks.data("hollander_rows", hollander_rows)
checks.data("mtcars", mtcars)
checks.save()
