"""Reference values for notebooks/19_tests_averages_spreads_shapes.ipynb, saved to references/19_tests_averages_spreads_shapes.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/19_tests_averages_spreads_shapes/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/19_tests_averages_spreads_shapes.py      (from the repository root; needs R with the packages listed)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks and realdata)

import numpy as np
import pandas as pd
import simplydata as sdata
from scipy import integrate, stats

import realdata as rd
from checks import Checks, quiet, r_available, r_run, read_row, versions

versions("pandas", "numpy", "scipy", "statsmodels", "rpy2")
checks = Recorder("19_tests_averages_spreads_shapes")
HAS_R = r_available()
CLOSED = 1e-10
EXACT = 1e-12
ITERATIVE = 1e-6
DOCSTRING = CLOSED   # a number a SciPy docstring prints (16-17 digits, from one SciPy build): the closed-formula class
if HAS_R:
    print(rd.r_version_of("TOSTER", "BSDA", "WRS2", "onewaytests", "PMCMRplus", "perm", "car", "randtests",
                          "nortest", "goftest", "fBasics", "tseries", "kSamples", "circular", "BayesFactor", "afex"))


def one(word, *args, **options):
    """The one-row answer of a simplydata test word, without its printout, read with table.get_row(1)."""
    return read_row(quiet(word, *args, **options))


def compare_many(dataset, method, found, reference, reference_name, tolerance=CLOSED, how="relative"):
    for quantity, value in reference.items():
        checks.compare(dataset, method, quantity, found[quantity], value, reference_name, tolerance, how)


def docstring_value(shown, source_part, index=0):
    """The number(s) a docstring prints for the example whose source contains source_part."""
    want = next(w for s, w in shown if source_part in s)
    numbers = [float(v) for v in __import__("re").findall(r"-?\d+\.\d+(?:e-?\d+)?|-?\d+(?:e-?\d+)", want)]
    return numbers[index]


def jzs_bayes_factor(t, n_effective, df, r):
    """Rouder et al. (2009), equation 1 with the Cauchy scale r: the JZS Bayes factor from t, by SciPy quad."""
    def integrand(g):
        a = 1 + n_effective * g * r ** 2
        return (a ** -0.5 * (1 + t ** 2 / (a * df)) ** (-(df + 1) / 2)
                * (2 * np.pi) ** -0.5 * g ** -1.5 * np.exp(-1 / (2 * g)))
    top = integrate.quad(integrand, 0, np.inf, epsabs=0, epsrel=1e-12, limit=500)[0]
    return top / (1 + t ** 2 / df) ** (-(df + 1) / 2)


morley = rd.r_dataset("morley")
light = sdata.table(morley)
light.compare_average_to("Speed", 792.458)

DS = "morley (Michelson 1879)"
if HAS_R:
    for opts, alternative in [({}, "two.sided"), ({"sides": 1, "higher": "average"}, "greater"),
                              ({"sides": 1, "lower": "average"}, "less"), ({"alternative": "less"}, "less")]:
        found = one(light.compare_average_to, "Speed", 792.458, **opts)
        r = r_run(f"""tt <- t.test(Speed, mu = 792.458, alternative = "{alternative}")
            list(statistic = unname(tt$statistic), `p-value` = tt$p.value, `degrees of freedom` = unname(tt$parameter),
                 `average 95 % CI low` = tt$conf.int[1], `average 95 % CI high` = tt$conf.int[2],
                 `effect size` = (mean(Speed) - 792.458) / sd(Speed))""", morley[["Speed"]])
        label = "compare_average_to" + (" (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")" if opts else "")
        if alternative != "two.sided":
            r = {k: v for k, v in r.items() if k in ("statistic", "p-value") or (alternative == "greater" and k.endswith("low"))
                 or (alternative == "less" and k.endswith("high"))}
        compare_many(DS, label, found, r, "R t.test")
    found = one(light.compare_average_to, "Speed", 792.458, test_method="z", standard_deviation=80)
    r = r_run("""n <- length(Speed); z <- (mean(Speed) - 792.458) / (80 / sqrt(n))
        b <- BSDA::z.test(Speed, mu = 792.458, sigma.x = 80)
        list(statistic = z, `p-value` = 2 * pnorm(-abs(z)), `average 95 % CI low` = mean(Speed) - qnorm(0.975) * 80 / sqrt(n),
             `average 95 % CI high` = mean(Speed) + qnorm(0.975) * 80 / sqrt(n), bsda_p = b$p.value)""", morley[["Speed"]])
    checks.compare(DS, 'compare_average_to (test_method="z", standard_deviation=)', "p-value", found["p-value"],
                   r.pop("bsda_p"), "R BSDA z.test", CLOSED, "relative")
    compare_many(DS, 'compare_average_to (test_method="z", standard_deviation=)', found, r, "z formula in R")


if HAS_R:
    r = r_run("""n <- length(Speed); d <- mean(Speed) - 792.458; se <- sd(Speed) / sqrt(n)
        plow <- pt((d + 50) / se, n - 1, lower.tail = FALSE); phigh <- pt((d - 50) / se, n - 1)
        ci <- t.test(Speed, mu = 792.458, conf.level = 0.9)$conf.int - 792.458
        list(low = plow, high = phigh, p = max(plow, phigh), cilow = ci[1], cihigh = ci[2])""", morley[["Speed"]])
    found = one(light.compare_average_to, "Speed", 792.458, same_within=50)
    compare_many(DS, "compare_average_to (same_within=)", found,
                 {"p-value against the low margin": r["low"], "p-value against the high margin": r["high"],
                  "p-value": r["p"], "difference 90 % CI low": r["cilow"], "difference 90 % CI high": r["cihigh"]},
                 "TOST written out in R")
    found = one(light.compare_average_to, "Speed", 792.458, same_within=50, not_lower="average")
    checks.compare(DS, 'compare_average_to (same_within=, not_lower="average")', "p-value", found["p-value"], r["low"],
                   "one-sided t written out in R", CLOSED, "relative")
    found = one(light.compare_average_to, "Speed", 792.458, same_within=50, not_higher="average")
    checks.compare(DS, 'compare_average_to (same_within=, not_higher="average")', "p-value", found["p-value"],
                   r["high"], "one-sided t written out in R", CLOSED, "relative")


darwin = pd.DataFrame({"difference": np.asarray(rd.docstring_data(stats.wilcoxon, "d"), dtype=float)})
names, shown = rd.docstring_example(stats.wilcoxon)
maize = sdata.table(darwin)
printed = [w for s, w in shown if s.startswith("res.statistic, res.pvalue")]
print(printed)
two_sided = [float(v) for v in printed[0].strip("()").split(",")]
greater = [float(v) for v in printed[1].strip("()").split(",")]
ds = "Darwin's maize (SciPy wilcoxon docstring)"
found = one(maize.compare_average_to, "difference", 0, based_on="ranks")
checks.compare(ds, 'compare_average_to (based_on="ranks")', "p-value", found["p-value"],
               two_sided[1], "SciPy wilcoxon docstring, printed", DOCSTRING, "relative")
found_greater = one(maize.compare_average_to, "difference", 0, based_on="ranks", sides=1, higher="average")
checks.compare(ds, 'compare_average_to (based_on="ranks", sides=1, higher="average")', "p-value",
               found_greater["p-value"], greater[1], "SciPy wilcoxon docstring, printed", DOCSTRING, "relative")
checks.compare(ds, 'compare_average_to (based_on="ranks", sides=1, higher="average")', "statistic (V)",
               found_greater["statistic"], greater[0], "SciPy wilcoxon docstring, printed", 0)
if HAS_R:
    for opts, alternative in [({}, "two.sided"), ({"sides": 1, "higher": "average"}, "greater"),
                              ({"alternative": "less"}, "less")]:
        r = r_run(f"""w <- wilcox.test(d, mu = 0, alternative = "{alternative}")
            list(statistic = unname(w$statistic), `p-value` = w$p.value)""", {"d": darwin["difference"]})
        found = one(maize.compare_average_to, "difference", 0, based_on="ranks", **opts)
        label = 'compare_average_to (based_on="ranks"' + "".join(f", {k}={v!r}" for k, v in opts.items()) + ")"
        compare_many(ds, label, found, r, "R wilcox.test (exact)")


def bayes_checks(dataset, word, word_args, t, n_eff, df, r_call, label, value=None):
    for width_name, width in [("wide", 1.0), (0.5, 0.5)]:
        found = one(word, *word_args, test_method="bayesian", prior_width=width_name, with_seed=7)
        checks.compare(dataset, f"{label} (test_method=\"bayesian\", prior_width={width_name!r})", "Bayes factor",
                       found["Bayes factor"], jzs_bayes_factor(t, n_eff, df, width),
                       "Rouder et al. (2009) integral, SciPy quad", ITERATIVE, "relative")
        if HAS_R:
            r = r_run(f"""suppressPackageStartupMessages(library(BayesFactor))
                b <- {r_call.format(width=width)}
                set.seed(7, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
                p <- posterior(b, iterations = 10000, progress = FALSE); q <- summary(p)$quantiles
                list(bf = exp(b@bayesFactor$bf), low = q[1, 1], mid = q[1, 3], high = q[1, 5])""")
            checks.compare(dataset, f"{label} (test_method=\"bayesian\", prior_width={width_name!r})", "Bayes factor",
                           found["Bayes factor"], r["bf"], "R BayesFactor ttestBF", ITERATIVE, "relative")
            estimate = [c for c in found if c.endswith("95 % credible interval low")][0].replace(
                " 95 % credible interval low", "")
            for key, column in [("low", f"{estimate} 95 % credible interval low"), ("mid", estimate),
                                ("high", f"{estimate} 95 % credible interval high")]:
                ref = r[key] + (value or 0)
                checks.compare(dataset, f"{label} (test_method=\"bayesian\", prior_width={width_name!r}, with_seed=7)", column, found[column], ref,
                               "R BayesFactor posterior, same seed", CLOSED, "relative")


d = darwin["difference"].to_numpy()
t_darwin = stats.ttest_1samp(d, 0).statistic
bayes_checks(ds, maize.compare_average_to, ("difference", 0), t_darwin, len(d), len(d) - 1,
             "ttestBF(x = c(" + ", ".join(map(str, d)) + "), rscale = {width})", "compare_average_to")


names, shown = rd.docstring_example(stats.f_oneway)
places = ["tillamook", "newport", "petersburg", "magadan", "tvarminne"]
mussels = pd.DataFrame([(p, v) for p in places for v in rd.docstring_data(stats.f_oneway, p)],
                       columns=["place", "length"])
shells = sdata.table(mussels)
shells.compare_averages_of("length", for_each_group_in="place", test_method="anova")

ds = "mussel shells (SciPy f_oneway docstring; McDonald et al. 1991)"
want = next(w for s, w in shown if s.startswith("f_oneway("))
F_printed, p_printed = rd.printed_number(want, "statistic"), rd.printed_number(want, "pvalue")
found = one(shells.compare_averages_of, "length", for_each_group_in="place", test_method="anova")
checks.compare(ds, 'compare_averages_of (test_method="anova")', "statistic", found["statistic"], F_printed,
               "SciPy f_oneway docstring, printed", DOCSTRING, "relative")
checks.compare(ds, 'compare_averages_of (test_method="anova")', "p-value", found["p-value"], p_printed,
               "SciPy f_oneway docstring, printed", DOCSTRING, "relative")
R_GROUPS = """g <- factor(group); m <- median(y); tab <- table(g, y > m)
    a <- anova(lm(y ~ g)); w <- oneway.test(y ~ g); k <- kruskal.test(y ~ g); md <- suppressWarnings(chisq.test(tab, correct = FALSE))
    list(aF = a$`F value`[1], ap = a$`Pr(>F)`[1], wF = unname(w$statistic), wp = w$p.value, wdf2 = unname(w$parameter[2]),
         kH = unname(k$statistic), kp = k$p.value, mchi = unname(md$statistic), mp = md$p.value)"""
if HAS_R:
    r = r_run(R_GROUPS, {"group": mussels["place"], "y": mussels["length"]})
    for test, keys, ref in [("anova", {"statistic": "aF", "p-value": "ap"}, "R anova(lm())"),
                            ("welch anova", {"statistic": "wF", "p-value": "wp", "degrees of freedom 2": "wdf2"},
                             "R oneway.test"),
                            ("kruskal-wallis", {"statistic": "kH", "p-value": "kp"}, "R kruskal.test"),
                            ("median", {"statistic": "mchi", "p-value": "mp"}, "R chisq.test on above/not above")]:
        found = one(shells.compare_averages_of, "length", for_each_group_in="place", test_method=test)
        compare_many(ds, f'compare_averages_of (test_method="{test}")', found, {k: r[v] for k, v in keys.items()}, ref)


names, shown = rd.docstring_example(stats.median_test)
groups3 = pd.DataFrame([(g, v) for g in ("g1", "g2", "g3") for v in rd.docstring_data(stats.median_test, g)],
                       columns=["group", "y"])
ds = "SciPy median_test docstring"
found = one(sdata.table(groups3).compare_averages_of, "y", for_each_group_in="group", test_method="median")
checks.compare(ds, 'compare_averages_of (test_method="median")', "p-value", found["p-value"], 0.12609082774093244,
               "SciPy median_test docstring, printed", DOCSTRING, "relative")
if HAS_R:
    r = r_run(R_GROUPS, {"group": groups3["group"], "y": groups3["y"]})
    checks.compare(ds, 'compare_averages_of (test_method="median")', "p-value", found["p-value"], r["mp"],
                   "R chisq.test on above/not above", CLOSED, "relative")

names, shown = rd.docstring_example(stats.alexandergovern)
cities = ["atlanta", "chicago", "houston", "memphis"]
houses = pd.DataFrame([(c, v) for c in cities for v in rd.docstring_data(stats.alexandergovern, c)],
                      columns=["city", "price"])
ds = "house prices (SciPy alexandergovern docstring)"
found = one(sdata.table(houses).compare_averages_of, "price", for_each_group_in="city", test_method="alexander-govern")
checks.compare(ds, 'compare_averages_of (test_method="alexander-govern")', "statistic", found["statistic"],
               4.65087071883494, "SciPy alexandergovern docstring, printed", DOCSTRING, "relative")
checks.compare(ds, 'compare_averages_of (test_method="alexander-govern")', "p-value", found["p-value"],
               0.19922132490385214, "SciPy alexandergovern docstring, printed", DOCSTRING, "relative")
if HAS_R:
    r = r_run("""d <- data.frame(y = y, g = factor(group)); a <- onewaytests::ag.test(y ~ g, data = d, verbose = FALSE)
        list(statistic = unname(a$statistic), `p-value` = a$p.value)""", {"group": houses["city"], "y": houses["price"]})
    compare_many(ds, 'compare_averages_of (test_method="alexander-govern")', found, r, "R onewaytests ag.test")


names, shown = rd.docstring_example(stats.mannwhitneyu)
ages = pd.DataFrame([("men", v) for v in rd.docstring_data(stats.mannwhitneyu, "males")]
                    + [("women", v) for v in rd.docstring_data(stats.mannwhitneyu, "females")],
                    columns=["sex", "age"])
ds = "ages (SciPy mannwhitneyu docstring)"
found = one(sdata.table(ages).compare_averages_of, "age", for_each_group_in="sex", test_method="mann-whitney")
checks.compare(ds, 'compare_averages_of (test_method="mann-whitney")', "statistic (U of men)", found["statistic"],
               17.0, "SciPy mannwhitneyu docstring, printed", 0)
if HAS_R:
    r = r_run("""w <- wilcox.test(y[g == "men"], y[g == "women"], exact = TRUE)
        list(statistic = unname(w$statistic), `p-value` = w$p.value)""", {"g": ages["sex"], "y": ages["age"]})
    compare_many(ds, 'compare_averages_of (test_method="mann-whitney")', found, r, "R wilcox.test (exact)")

names, shown = rd.docstring_example(stats.brunnermunzel)
pain = pd.DataFrame([("x1", v) for v in rd.docstring_data(stats.brunnermunzel, "x1")]
                    + [("x2", v) for v in rd.docstring_data(stats.brunnermunzel, "x2")], columns=["group", "score"])
ds = "pain scores (SciPy brunnermunzel docstring; Brunner & Munzel 2000)"
found = one(sdata.table(pain).compare_averages_of, "score", for_each_group_in="group", test_method="brunner-munzel")
checks.compare(ds, 'compare_averages_of (test_method="brunner-munzel")', "statistic", found["statistic"],
               3.1374674823029505, "SciPy brunnermunzel docstring, printed", DOCSTRING, "relative")
checks.compare(ds, 'compare_averages_of (test_method="brunner-munzel")', "p-value", found["p-value"],
               0.0057862086661515377, "SciPy brunnermunzel docstring, printed", DOCSTRING, "relative")
if HAS_R:
    r = r_run("""x <- y[g == "x1"]; z <- y[g == "x2"]; n1 <- length(x); n2 <- length(z)
        rk <- rank(c(x, z)); r1 <- rk[1:n1]; r2 <- rk[-(1:n1)]; i1 <- rank(x); i2 <- rank(z)
        s1 <- sum((r1 - i1 - mean(r1) + (n1 + 1) / 2)^2) / (n1 - 1); s2 <- sum((r2 - i2 - mean(r2) + (n2 + 1) / 2)^2) / (n2 - 1)
        w <- n1 * n2 * (mean(r2) - mean(r1)) / ((n1 + n2) * sqrt(n1 * s1 + n2 * s2))
        df <- (n1 * s1 + n2 * s2)^2 / ((n1 * s1)^2 / (n1 - 1) + (n2 * s2)^2 / (n2 - 1))
        list(statistic = w, `p-value` = 2 * pt(-abs(w), df))""", {"g": pain["group"], "y": pain["score"]})
    compare_many(ds, 'compare_averages_of (test_method="brunner-munzel")', found, r, "Brunner-Munzel formula in R")


sleep = rd.r_dataset("sleep")
sleep["drug"] = "drug " + sleep["group"].astype(str)
sleepers = sdata.table(sleep)
ds = "sleep (Cushny & Peebles 1905), as two groups"
if HAS_R:
    for test, var_equal in [("welch", "FALSE"), ("student", "TRUE")]:
        for opts, alternative in [({}, "two.sided"), ({"sides": 1, "higher": "drug 2"}, "less"),
                                  ({"sides": 1, "lower": "drug 2"}, "greater"), ({"alternative": "less"}, "less"),
                                  ({"alternative": "two-sided"}, "two.sided")]:
            found = one(sleepers.compare_averages_of, "extra", for_each_group_in="drug", test_method=test, **opts)
            r = r_run(f"""tt <- t.test(extra[drug == "drug 1"], extra[drug == "drug 2"], var.equal = {var_equal},
                              alternative = "{alternative}")
                list(statistic = unname(tt$statistic), `p-value` = tt$p.value, `degrees of freedom` = unname(tt$parameter),
                     `difference 95 % CI low` = tt$conf.int[1], `difference 95 % CI high` = tt$conf.int[2])""",
                      sleep[["extra", "drug"]])
            if alternative == "less":
                r.pop("difference 95 % CI low")
            if alternative == "greater":
                r.pop("difference 95 % CI high")
            label = f'compare_averages_of (test_method="{test}"' + "".join(f", {k}={v!r}" for k, v in opts.items()) + ")"
            compare_many(ds, label, found, r, "R t.test")
    found = one(sleepers.compare_averages_of, "extra", for_each_group_in="drug", test_method="yuen")
    r = r_run("""y <- WRS2::yuen(extra ~ drug, data = data.frame(extra = extra, drug = drug), tr = 0.2)
        list(statistic = unname(y$test), `p-value` = y$p.value, `degrees of freedom` = unname(y$df))""", sleep[["extra", "drug"]])
    # WRS2 reports the size of the test statistic (no sign); simplydata keeps the sign of group 1 minus group 2
    checks.compare(ds, 'compare_averages_of (test_method="yuen")', "statistic (size)", abs(found["statistic"]),
                   r["statistic"], "R WRS2 yuen", CLOSED, "relative")
    r.pop("statistic")
    compare_many(ds, 'compare_averages_of (test_method="yuen")', found, r, "R WRS2 yuen")

a = sleep.loc[sleep["drug"] == "drug 1", "extra"].to_numpy()
b = sleep.loc[sleep["drug"] == "drug 2", "extra"].to_numpy()
found = one(sleepers.compare_averages_of, "extra", for_each_group_in="drug", p_value_by="shuffling", shuffles=9999,
            with_seed=1)
parent = stats.permutation_test((a, b), lambda x, y, axis=-1: np.mean(x, axis=axis) - np.mean(y, axis=axis),
                                n_resamples=9999, vectorized=True, rng=np.random.default_rng(1))
checks.compare(ds, 'compare_averages_of (p_value_by="shuffling", shuffles=9999, with_seed=1)', "p-value",
               found["p-value"], parent.pvalue, "SciPy permutation_test, same seed", EXACT, "relative")
found = one(sleepers.compare_averages_of, "extra", for_each_group_in="drug", p_value_by="shuffling", shuffles=200000)
if HAS_R:
    r = r_run("""suppressPackageStartupMessages(library(perm))
        list(p = permTS(extra[drug == "drug 1"], extra[drug == "drug 2"], method = "exact.network")$p.value)""",
              sleep[["extra", "drug"]])
    checks.compare(ds, 'compare_averages_of (p_value_by="shuffling", shuffles=200000: every split)', "p-value",
                   found["p-value"], r["p"], "R perm permTS (exact)", EXACT, "relative")


mtcars = rd.r_dataset("mtcars")
mtcars["transmission"] = np.where(mtcars["am"] == 1, "manual", "automatic")
cars = sdata.table(mtcars)
ds = "mtcars mpg by transmission (TOSTER t_TOST example 1)"
if HAS_R:
    r = r_run("""x <- mpg[am == 0]; y <- mpg[am == 1]; d <- mean(x) - mean(y)
        se <- sqrt(var(x) / length(x) + var(y) / length(y))
        df <- se^4 / ((var(x) / length(x))^2 / (length(x) - 1) + (var(y) / length(y))^2 / (length(y) - 1))
        plow <- pt((d + 3) / se, df, lower.tail = FALSE); phigh <- pt((d - 3) / se, df)
        res <- TOSTER::t_TOST(mpg ~ am, data = data.frame(mpg = mpg, am = am), eqb = 3)
        list(low = plow, high = phigh, tlow = res$TOST$p.value[2], thigh = res$TOST$p.value[3])""", mtcars[["mpg", "am"]])
    found = one(cars.compare_averages_of, "mpg", for_each_group_in="transmission", same_within=3)
    print(found)
    compare_many(ds, "compare_averages_of (same_within=)", found,
                 {"p-value against the low margin": r["low"], "p-value against the high margin": r["high"]},
                 "Welch TOST written out in R")
    compare_many(ds, "compare_averages_of (same_within=)", found,
                 {"p-value against the low margin": r["tlow"], "p-value against the high margin": r["thigh"]},
                 "R TOSTER t_TOST (help page example)")
    found = one(cars.compare_averages_of, "mpg", for_each_group_in="transmission", same_within=3, not_lower="manual")
    print(found["what it compares"], found["p-value"])
    checks.compare(ds, 'compare_averages_of (same_within=, not_lower="manual")', "p-value", found["p-value"],
                   r["high"], "one-sided Welch t written out in R", CLOSED, "relative")
    found = one(cars.compare_averages_of, "mpg", for_each_group_in="transmission", same_within=3, not_higher="manual")
    checks.compare(ds, 'compare_averages_of (same_within=, not_higher="manual")', "p-value", found["p-value"],
                   r["low"], "one-sided Welch t written out in R", CLOSED, "relative")

t_sleep = stats.ttest_ind(a, b).statistic
for width_name, width in [("wide", 1.0), (0.5, 0.5)]:
    found = one(sleepers.compare_averages_of, "extra", for_each_group_in="drug", test_method="bayesian",
                prior_width=width_name, with_seed=7)
    label = f'compare_averages_of (test_method="bayesian", prior_width={width_name!r})'
    checks.compare(ds := "sleep (Cushny & Peebles 1905), as two groups", label, "Bayes factor", found["Bayes factor"],
                   jzs_bayes_factor(t_sleep, len(a) * len(b) / (len(a) + len(b)), len(a) + len(b) - 2, width),
                   "Rouder et al. (2009) integral, SciPy quad", ITERATIVE, "relative")
    if HAS_R:
        r = r_run(f"""suppressPackageStartupMessages(library(BayesFactor))
            b <- ttestBF(x = x, y = y, paired = FALSE, rscale = {width})
            set.seed(7, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
            p <- posterior(b, iterations = 10000, progress = FALSE); q <- summary(p)$quantiles
            list(bf = exp(b@bayesFactor$bf), low = q["beta (x - y)", 1], mid = q["beta (x - y)", 3], high = q["beta (x - y)", 5])""",
                  {"x": a, "y": b})
        checks.compare(ds, label, "Bayes factor", found["Bayes factor"], r["bf"], "R BayesFactor ttestBF", ITERATIVE,
                       "relative")
        for key, column in [("low", "difference 95 % credible interval low"), ("mid", "difference"),
                            ("high", "difference 95 % credible interval high")]:
            checks.compare(ds, f'compare_averages_of (test_method="bayesian", prior_width={width_name!r}, with_seed=7)', column, found[column],
                           r[key], "R BayesFactor posterior, same seed", CLOSED, "relative")


warpbreaks = rd.r_dataset("warpbreaks")
ds = "warpbreaks (Tippett 1950)"
found = one(sdata.table(warpbreaks).compare_averages_of, "breaks", for_each_group_in="tension", test_method="jonckheere",
            in_order=["L", "M", "H"])
print(found)
if HAS_R:
    r = r_run("""j <- suppressWarnings(PMCMRplus::jonckheereTest(breaks ~ factor(tension, levels = c("L", "M", "H")),
                    data = data.frame(breaks = breaks, tension = tension), alternative = "two.sided"))
        list(p = j$p.value, z = unname(j$statistic))""", warpbreaks)
    checks.compare(ds, 'compare_averages_of (test_method="jonckheere", in_order=)', "p-value", found["p-value"],
                   r["p"], "R PMCMRplus jonckheereTest", CLOSED, "relative")

co2 = rd.r_dataset("CO2")
co2["conc"] = co2["conc"].astype(int).astype(str)
ds = "CO2 uptake (Potvin et al. 1990)"
found = one(sdata.table(co2).compare_averages_of, "uptake", for_each_group_in="conc", matched_by="Plant")
print(found)
if HAS_R:
    r = r_run("""d <- data.frame(uptake = uptake, conc = factor(conc), Plant = factor(Plant))
        s <- summary(aov(uptake ~ conc + Error(Plant/conc), data = d))[["Error: Plant:conc"]][[1]]
        F <- s$`F value`[1]; df1 <- s$Df[1]; df2 <- s$Df[2]
        w <- reshape(d, idvar = "Plant", timevar = "conc", direction = "wide")[, -1]
        k <- ncol(w); C <- t(contr.poly(k)); S <- C %*% cov(w) %*% t(C)
        eps <- sum(diag(S))^2 / ((k - 1) * sum(S^2))
        list(F = F, p = pf(F, df1, df2, lower.tail = FALSE), eps = eps, pgg = pf(F, eps * df1, eps * df2, lower.tail = FALSE))""",
              co2[["uptake", "conc", "Plant"]])
    checks.compare(ds, 'compare_averages_of (matched_by="Plant")', "statistic (F)", found["statistic"], r["F"],
                   "R aov with Error()", CLOSED, "relative")
    pcols = [c for c in found if "p-value" in c]
    gg = [c for c in pcols if "Greenhouse" in c or "corrected" in c]
    checks.compare(ds, 'compare_averages_of (matched_by="Plant")', "p-value (Greenhouse-Geisser)",
                   found[gg[0] if gg else "p-value"], r["pgg"], "epsilon from the covariance matrix, R pf", CLOSED,
                   "relative")
    eps = [c for c in found if "epsilon" in c.lower()]
    if eps:
        checks.compare(ds, 'compare_averages_of (matched_by="Plant")', "Greenhouse-Geisser epsilon", found[eps[0]],
                       r["eps"], "epsilon from the covariance matrix", CLOSED, "relative")


made = rd.r_example_objects("stats", "wilcox.test")
hw = r_run("list(x = x, y = y)")
depression = pd.DataFrame({"first": hw["x"], "second": hw["y"]})
visits = sdata.table(depression)
ds = "depression scores (R wilcox.test example; Hollander & Wolfe 1973)"
if HAS_R:
    for opts, alternative in [({"sides": 1, "higher": "before"}, "greater"), ({}, "two.sided"),
                              ({"sides": 1, "lower": "before"}, "less"), ({"alternative": "less"}, "greater")]:
        r = r_run(f"""w <- wilcox.test(x, y, paired = TRUE, alternative = "{alternative}")
            list(`p-value` = w$p.value, statistic = unname(w$statistic))""", {"x": depression["first"], "y": depression["second"]})
        found = one(visits.compare_before_and_after, "first", "second", based_on="ranks", **opts)
        label = 'compare_before_and_after (based_on="ranks"' + "".join(f", {k}={v!r}" for k, v in opts.items()) + ")"
        checks.compare(ds, label, "p-value", found["p-value"], r["p-value"], "R wilcox.test (exact)", CLOSED, "relative")


wide = sleep.pivot(index="ID", columns="group", values="extra")
paired = pd.DataFrame({"drug1": wide["1"].to_numpy(), "drug2": wide["2"].to_numpy()})
patients = sdata.table(paired)
ds = "sleep (Cushny & Peebles 1905), paired"
if HAS_R:
    for opts, alternative in [({}, "two.sided"), ({"sides": 1, "higher": "after"}, "greater"),
                              ({"sides": 1, "lower": "after"}, "less"), ({"alternative": "greater"}, "greater")]:
        r = r_run(f"""tt <- t.test(drug2, drug1, paired = TRUE, alternative = "{alternative}")
            list(statistic = unname(tt$statistic), `p-value` = tt$p.value)""", paired)
        found = one(patients.compare_before_and_after, "drug1", "drug2", **opts)
        label = "compare_before_and_after" + (" (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")" if opts else "")
        compare_many(ds, label, found, r, "R t.test(paired = TRUE)")
    found = one(patients.compare_before_and_after, "drug1", "drug2", test_method="sign")
    r = r_run("""d <- drug2 - drug1; list(p = binom.test(sum(d > 0), sum(d != 0))$p.value, up = sum(d > 0))""", paired)
    checks.compare(ds, 'compare_before_and_after (test_method="sign")', "p-value", found["p-value"], r["p"],
                   "R binom.test on the signs", CLOSED, "relative")
    checks.compare(ds, 'compare_before_and_after (test_method="sign")', "rows that went up", found["rows that went up"],
                   r["up"], "R", 0)
change = (paired["drug2"] - paired["drug1"]).to_numpy()
signs = np.array(np.meshgrid(*[[-1, 1]] * len(change))).reshape(len(change), -1).T
means = (signs * change).mean(axis=1)
exact_p = np.mean(np.abs(means) >= abs(change.mean()) * (1 - 1e-12))
found = one(patients.compare_before_and_after, "drug1", "drug2", p_value_by="shuffling", shuffles=9999, with_seed=1)
checks.compare(ds, 'compare_before_and_after (p_value_by="shuffling", shuffles=, with_seed=)', "p-value",
               found["p-value"], exact_p, "all 1024 sign patterns, NumPy", EXACT, "relative")
if HAS_R:
    r = r_run("""res <- TOSTER::t_TOST(x = drug1, y = drug2, paired = TRUE, eqb = c(-0.5, 2))
        d <- drug1 - drug2; n <- length(d); se <- sd(d) / sqrt(n)
        list(tlow = res$TOST$p.value[2], thigh = res$TOST$p.value[3],
             low = pt((mean(d) + 0.5) / se, n - 1, lower.tail = FALSE), high = pt((mean(d) - 2) / se, n - 1))""", paired)
    # TOSTER tests drug1 - drug2; simplydata's change is after - before, so before = drug2 and after = drug1
    found = one(patients.compare_before_and_after, "drug2", "drug1", same_within=[-0.5, 2])
    compare_many("sleep, paired (TOSTER t_TOST example 2)", "compare_before_and_after (same_within=)", found,
                 {"p-value against the low margin": r["low"], "p-value against the high margin": r["high"]},
                 "paired TOST written out in R")
    compare_many("sleep, paired (TOSTER t_TOST example 2)", "compare_before_and_after (same_within=)", found,
                 {"p-value against the low margin": r["tlow"], "p-value against the high margin": r["thigh"]},
                 "R TOSTER t_TOST (help page example)")
    found = one(patients.compare_before_and_after, "drug2", "drug1", same_within=[-0.5, 2], not_lower="after")
    checks.compare("sleep, paired (TOSTER t_TOST example 2)", 'compare_before_and_after (same_within=, not_lower="after")',
                   "p-value", found["p-value"], r["low"], "one-sided paired t written out in R", CLOSED, "relative")
    found = one(patients.compare_before_and_after, "drug2", "drug1", same_within=[-0.5, 2], not_higher="after")
    checks.compare("sleep, paired (TOSTER t_TOST example 2)", 'compare_before_and_after (same_within=, not_higher="after")',
                   "p-value", found["p-value"], r["high"], "one-sided paired t written out in R", CLOSED, "relative")

t_paired = stats.ttest_rel(paired["drug2"], paired["drug1"]).statistic
for width_name, width in [("wide", 1.0), (0.5, 0.5)]:
    found = one(patients.compare_before_and_after, "drug1", "drug2", test_method="bayesian", prior_width=width_name,
                with_seed=7)
    checks.compare(ds, f'compare_before_and_after (test_method="bayesian", prior_width={width_name!r})', "Bayes factor",
                   found["Bayes factor"], jzs_bayes_factor(t_paired, len(change), len(change) - 1, width),
                   "Rouder et al. (2009) integral, SciPy quad", ITERATIVE, "relative")


made = rd.r_example_objects("stats", "mcnemar.test")
perf = r_run("list(n = as.vector(Performance))")["n"]  # column-major: (A,A), (D,A), (A,D), (D,D)
cells = [("approve", "approve", perf[0]), ("disapprove", "approve", perf[1]), ("approve", "disapprove", perf[2]),
         ("disapprove", "disapprove", perf[3])]
approval = pd.DataFrame([(a, b) for a, b, n in cells for _ in range(int(n))], columns=["first survey", "second survey"])
surveys = sdata.table(approval)
ds = "presidential approval (R mcnemar.test example; Agresti 1990)"
found = one(surveys.compare_before_and_after, "first survey", "second survey")
print(found)
if HAS_R:
    r = r_run("""m <- mcnemar.test(Performance); list(statistic = unname(m$statistic), `p-value` = m$p.value,
        exact = binom.test(Performance[1, 2], Performance[1, 2] + Performance[2, 1])$p.value)""")
    compare_many(ds, "compare_before_and_after (yes/no columns: McNemar)", found,
                 {"statistic": r["statistic"], "p-value": r["p-value"]}, "R mcnemar.test")
    found = one(surveys.compare_before_and_after, "first survey", "second survey", p_value_by="exact")
    checks.compare(ds, 'compare_before_and_after (p_value_by="exact")', "p-value", found["p-value"], r["exact"],
                   "R binom.test on the changes", CLOSED, "relative")


plants = rd.r_dataset("PlantGrowth")
summary = plants.groupby("group")["weight"].agg(["mean", "std", "count"])
ds = "PlantGrowth (Dobson 1983), from summary numbers"
if HAS_R:
    two = summary.loc[["ctrl", "trt1"]]
    for test, var_equal in [(None, "FALSE"), ("student", "TRUE")]:
        found = one(sdata.compare_averages, averages=two["mean"].tolist(), standard_deviations=two["std"].tolist(),
                    counts=two["count"].tolist(), test_method=test)
        r = r_run(f"""tt <- t.test(weight[group == "ctrl"], weight[group == "trt1"], var.equal = {var_equal})
            list(statistic = unname(tt$statistic), `p-value` = tt$p.value, `degrees of freedom` = unname(tt$parameter),
                 `difference 95 % CI low` = tt$conf.int[1], `difference 95 % CI high` = tt$conf.int[2])""", plants)
        compare_many(ds, "sdata.compare_averages" + (f' (test_method="{test}")' if test else ""), found, r,
                     "R t.test on the rows")
    for test, var_equal in [(None, "FALSE"), ("anova", "TRUE")]:
        found = one(sdata.compare_averages, averages=summary["mean"].tolist(),
                    standard_deviations=summary["std"].tolist(), counts=summary["count"].tolist(), test_method=test)
        r = r_run(f"""o <- oneway.test(weight ~ group, data = data.frame(weight = weight, group = group), var.equal = {var_equal})
            list(statistic = unname(o$statistic), `p-value` = o$p.value)""", plants)
        compare_many(ds, "sdata.compare_averages" + (f' (test_method="{test}")' if test else " (3 groups)"), found, r,
                     "R oneway.test on the rows")


made = rd.r_example_objects("stats", "ansari.test")
iron = r_run("list(a = ramsay, b = jung.parekh)")
serum = pd.DataFrame([("ramsay", v) for v in iron["a"]] + [("jung.parekh", v) for v in iron["b"]],
                     columns=["method", "iron"])
ds = "serum iron (R ansari.test / mood.test example; Hollander & Wolfe 1973)"
serum_t = sdata.table(serum)
if HAS_R:
    r = r_run("""a <- ansari.test(iron[method == "jung.parekh"], iron[method == "ramsay"])
        m <- mood.test(iron[method == "jung.parekh"], iron[method == "ramsay"])
        v <- var.test(iron[method == "jung.parekh"], iron[method == "ramsay"])
        list(ap = a$p.value, astat = unname(a$statistic), mp = m$p.value, mz = unname(m$statistic), vp = v$p.value, vF = unname(v$statistic))""",
              serum)
    found = one(serum_t.compare_spreads_of, "iron", for_each_group_in="method", test_method="ansari")
    print(found)
    checks.compare(ds, 'compare_spreads_of (test_method="ansari")', "p-value", found["p-value"], r["ap"],
                   "R ansari.test", CLOSED, "relative")
    found = one(serum_t.compare_spreads_of, "iron", for_each_group_in="method", test_method="mood")
    checks.compare(ds, 'compare_spreads_of (test_method="mood")', "p-value", found["p-value"], r["mp"],
                   "R mood.test", CLOSED, "relative")
    found = one(serum_t.compare_spreads_of, "iron", for_each_group_in="method", test_method="f")
    first, second = sorted(serum["method"].unique())
    va, vb = (serum.loc[serum["method"] == g, "iron"].var(ddof=1) for g in (first, second))
    na, nb = (int((serum["method"] == g).sum()) for g in (first, second))
    F = va / vb
    p_f = 2 * min(stats.f.cdf(F, na - 1, nb - 1), stats.f.sf(F, na - 1, nb - 1))
    checks.compare(ds, 'compare_spreads_of (test_method="f")', "statistic (F)", found["statistic"], F,
                   "variance ratio, NumPy", CLOSED, "relative")
    checks.compare(ds, 'compare_spreads_of (test_method="f")', "p-value", found["p-value"], p_f, "F formula, SciPy",
                   CLOSED, "relative")
    checks.compare(ds, 'compare_spreads_of (test_method="f")', "p-value", found["p-value"], r["vp"], "R var.test",
                   CLOSED, "relative")

names, shown = rd.docstring_example(stats.levene)
lev = pd.DataFrame([(g, v) for g in ("a", "b", "c") for v in rd.docstring_data(stats.levene, g)],
                   columns=["group", "y"])
for name, frame, column, group, printed in [("SciPy levene docstring", lev, "y", "group", 0.002431505967249681),
                                            ("PlantGrowth (Dobson 1983)", plants, "weight", "group", None)]:
    for test in ("levene", "brown-forsythe"):
        found = one(sdata.table(frame).compare_spreads_of, column, for_each_group_in=group, test_method=test)
        if printed is not None and test == "levene":
            checks.compare(name, f'compare_spreads_of (test_method="{test}")', "p-value", found["p-value"], printed,
                           "SciPy levene docstring, printed", DOCSTRING, "relative")
        if HAS_R:
            r = r_run("""l <- car::leveneTest(y ~ factor(g), data = data.frame(y = y, g = g), center = median)
                list(statistic = l$`F value`[1], `p-value` = l$`Pr(>F)`[1])""", {"y": frame[column], "g": frame[group]})
            compare_many(name, f'compare_spreads_of (test_method="{test}")', found, r, "R car leveneTest(center = median)")

insects = rd.r_dataset("InsectSprays")
if HAS_R:
    r = r_run("""b <- bartlett.test(count ~ spray, data = data.frame(count = count, spray = spray))
        f <- fligner.test(count ~ spray, data = data.frame(count = count, spray = spray))
        list(bs = unname(b$statistic), bp = b$p.value, fs = unname(f$statistic), fp = f$p.value)""", insects)
    for test, keys, ref in [("bartlett", ("bs", "bp"), "R bartlett.test (help page example)"),
                            ("fligner", ("fs", "fp"), "R fligner.test (help page example)")]:
        found = one(sdata.table(insects).compare_spreads_of, "count", for_each_group_in="spray", test_method=test)
        compare_many("InsectSprays (Beall 1942)", f'compare_spreads_of (test_method="{test}")', found,
                     {"statistic": r[keys[0]], "p-value": r[keys[1]]}, ref)


if HAS_R:
    two_groups = plants[plants["group"].isin(["ctrl", "trt1"])]
    for name, frame, column, group in [("PlantGrowth ctrl, trt1", two_groups, "weight", "group"),
                                       ("serum iron (Hollander & Wolfe)", serum, "iron", "method")]:
        found = one(sdata.table(frame).compare_shapes_of, column, for_each_group_in=group)
        g1, g2 = sorted(frame[group].unique())
        r = r_run(f"""k <- suppressWarnings(ks.test(y[g == "{g1}"], y[g == "{g2}"]))
            list(statistic = unname(k$statistic), `p-value` = k$p.value)""", {"y": frame[column], "g": frame[group]})
        compare_many(name, "compare_shapes_of", found, r, "R ks.test")

made = rd.r_example_objects("kSamples", "ad.test")
u = r_run("list(u1 = u1, u2 = u2, u3 = u3)")
three = pd.DataFrame([(k, v) for k in ("u1", "u2", "u3") for v in u[k]], columns=["sample", "value"])
ds = "R kSamples ad.test example"
found = one(sdata.table(three).compare_shapes_of, "value", for_each_group_in="sample", shuffles=9999, with_seed=1)
print(found)
parent = stats.anderson_ksamp([np.asarray(u[k]) for k in ("u1", "u2", "u3")], variant="midrank",
                              method=stats.PermutationMethod(n_resamples=9999, rng=np.random.default_rng(1)))
checks.compare(ds, "compare_shapes_of (3 groups, shuffles=9999, with_seed=1)", "p-value", found["p-value"],
               parent.pvalue, "SciPy anderson_ksamp, same seed", EXACT, "relative")
if HAS_R:
    r = r_run("""set.seed(2627); a <- kSamples::ad.test(u1, u2, u3, method = "exact", dist = FALSE, Nsim = 1000)
        set.seed(1); long <- kSamples::ad.test(u1, u2, u3, method = "simulated", dist = FALSE, Nsim = 200000)
        list(stat = a$ad[2, 2], p = a$ad[2, 4], longp = long$ad[2, 4])""")
    checks.compare(ds, "compare_shapes_of (3 groups)", "statistic (standardised, version 2)", found["statistic"],
                   r["stat"], "R kSamples ad.test (help page example)", 1e-4, "relative")
    p = r["p"]
    checks.compare(ds, "compare_shapes_of (3 groups, shuffles=9999, with_seed=1)", "p-value", found["p-value"], p,
                   "R kSamples ad.test, help page run (1000 draws; 4 MCSE)",
                   4 * np.sqrt(p * (1 - p) / 9999 + p * (1 - p) / 1000))
    p = r["longp"]
    checks.compare(ds, "compare_shapes_of (3 groups, shuffles=9999, with_seed=1)", "p-value", found["p-value"], p,
                   "R kSamples ad.test, 200 000 draws (4 MCSE)", 4 * np.sqrt(p * (1 - p) / 9999 + p * (1 - p) / 200000))

made = rd.r_example_objects("randtests", "runs.test")
earth = pd.DataFrame({"density": r_run("list(x = earthden)")["x"]})
potatoes = rd.r_dataset("sweetpotato", "randtests")
for name, frame, column in [("density of the earth (R randtests runs.test example; Brownlee 1965)", earth, "density"),
                            ("sweet potato yield 1868-1937 (R randtests)", potatoes, "yield")]:
    found = one(sdata.table(frame).check_if_random_order, column)
    if HAS_R:
        r = r_run("""rt <- randtests::runs.test(y); list(statistic = unname(rt$statistic), `p-value` = rt$p.value,
            runs = unname(rt$runs))""", {"y": frame[column]})
        compare_many(name, "check_if_random_order", found, r, "R randtests runs.test")


for name, frame, column in [("morley (Michelson 1879)", morley, "Speed"),
                            ("density of the earth (Brownlee 1965)", earth, "density")]:
    t = sdata.table(frame)
    x = frame[column].to_numpy()
    if not HAS_R:
        break
    r = r_run("""list(sw = shapiro.test(y)$statistic, swp = shapiro.test(y)$p.value,
        lil = nortest::lillie.test(y)$statistic, lilp = nortest::lillie.test(y)$p.value,
        dag = fBasics::dagoTest(y)@test$statistic[1], dagp = fBasics::dagoTest(y)@test$p.value[1],
        jb = unname(tseries::jarque.bera.test(y)$statistic), jbp = tseries::jarque.bera.test(y)$p.value,
        x = 0)""", {"y": x})
    found = one(t.compare_shape_to, column)
    compare_many(name, "compare_shape_to (Shapiro-Wilk)", found, {"statistic": r["sw"], "p-value": r["swp"]},
                 "R shapiro.test")
    found = one(t.compare_shape_to, column, test_method="lilliefors")
    z = np.sort((x - x.mean()) / x.std(ddof=1))
    n = len(x)
    cdf = stats.norm.cdf(z)
    D = max(np.max(np.arange(1, n + 1) / n - cdf), np.max(cdf - np.arange(n) / n))
    checks.compare(name, 'compare_shape_to (test_method="lilliefors")', "statistic", found["statistic"], D,
                   "largest gap, NumPy", CLOSED, "relative")
    compare_many(name, 'compare_shape_to (test_method="lilliefors")', found, {"statistic": r["lil"], "p-value": r["lilp"]},
                 "R nortest lillie.test")
    found = one(t.compare_shape_to, column, test_method="dagostino")
    compare_many(name, 'compare_shape_to (test_method="dagostino")', found, {"statistic": r["dag"], "p-value": r["dagp"]},
                 "R fBasics dagoTest")
    found = one(t.compare_shape_to, column, test_method="jarque-bera")
    compare_many(name, 'compare_shape_to (test_method="jarque-bera")', found, {"statistic": r["jb"], "p-value": r["jbp"]},
                 "R tseries jarque.bera.test")

# against a fixed normal distribution, on values without ties (Darwin's maize differences: the exact
# Kolmogorov-Smirnov p-value assumes no ties; morley and the earth densities have ties)
x = darwin["difference"].to_numpy()
fixed = sdata.distribution("normal", average=float(x.mean()), standard_deviation=float(1.1 * x.std(ddof=1)))
if HAS_R:
    r = r_run("""m <- mean(y); s <- 1.1 * sd(y); k <- ks.test(y, "pnorm", m, s, exact = TRUE)
        cv <- goftest::cvm.test(y, "pnorm", mean = m, sd = s)
        list(ks = unname(k$statistic), ksp = k$p.value, cvm = unname(cv$statistic), cvmp = cv$p.value)""", {"y": x})
    found = one(maize.compare_shape_to, "difference", fixed)
    compare_many("Darwin's maize (SciPy wilcoxon docstring)", "compare_shape_to (shape=a fixed distribution: Kolmogorov-Smirnov)",
                 found, {"statistic": r["ks"], "p-value": r["ksp"]}, "R ks.test(exact = TRUE)")
    found = one(maize.compare_shape_to, "difference", fixed, test_method="cramer-von-mises")
    compare_many("Darwin's maize (SciPy wilcoxon docstring)", 'compare_shape_to (test_method="cramer-von-mises")', found,
                 {"statistic": r["cvm"], "p-value": r["cvmp"]}, "R goftest cvm.test", ITERATIVE)


t = sdata.table(potatoes)
fit = quiet(t.fit_distribution_to, "yield", "gamma").show_parameters()
print(fit)
ds = "sweet potato yield (R randtests)"
found = one(t.compare_shape_to, "yield", "gamma", samples=9999, with_seed=1)
again = one(t.compare_shape_to, "yield", "gamma", samples=9999, with_seed=1)
fewer = one(t.compare_shape_to, "yield", "gamma", samples=999, with_seed=2)
print(found)
checks.compare(ds, 'compare_shape_to (shape="gamma", samples=9999, with_seed=1)', "p-value, run again",
               found["p-value"], again["p-value"], "simplydata, same seed", 0)
p = found["p-value"]
checks.compare(ds, 'compare_shape_to (shape="gamma", samples=)', "p-value, 999 against 9999 samples", fewer["p-value"],
               p, "simplydata 9999 samples (4 MCSE)", 4 * np.sqrt(p * (1 - p) / 999 + p * (1 - p) / 9999))
if HAS_R:
    params = list(fit.values())
    shape_value = fit.get("shape")
    scale_value = fit.get("scale", None)
    rate_value = fit.get("rate", 1 / scale_value if scale_value else None)
    r = r_run(f"""list(ad = unname(goftest::ad.test(y, "pgamma", shape = {shape_value!r}, rate = {rate_value!r},
                   estimated = FALSE)$statistic))""", {"y": potatoes["yield"]})
    checks.compare(ds, 'compare_shape_to (shape="gamma")', "statistic (Anderson-Darling)", found["statistic"], r["ad"],
                   "R goftest ad.test at simplydata's fit", CLOSED, "relative")


wind = rd.r_dataset("wind", "circular")
wind["degrees"] = np.degrees(wind["d"])
ds = "wind directions (R circular)"
if HAS_R:
    r = r_run("""n <- length(d); C <- sum(cos(d)); S <- sum(sin(d)); Rbar <- sqrt(C^2 + S^2) / n; Rn <- n * Rbar
        Z <- n * Rbar^2
        series <- if (n < 50) 1 + (2 * Z - Z^2) / (4 * n) - (24 * Z - 132 * Z^2 + 76 * Z^3 - 9 * Z^4) / (288 * n^2) else 1
        list(Rbar = Rbar, p = exp(-Z) * series,
             circ = circular::rayleigh.test(circular::circular(d))$p.value)""", wind[["d"]])
    for column, turn, opts in [("d", "radians", {"shape": "even around the circle"}),
                               ("degrees", 360, {"shape": "even around the circle"}),
                               ("degrees", 360, {"shape": "uniform"}),
                               ("degrees", 360, {"shape": "even around the circle", "test_method": "rayleigh"})]:
        found = one(sdata.table(wind).compare_shape_to, column, full_turn=turn, **opts)
        label = "compare_shape_to (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + f", full_turn={turn!r})"
        checks.compare(ds, label, "p-value", found["p-value"], r["p"], "Rayleigh test written out in R", CLOSED,
                       "relative")
        checks.compare(ds, label, "p-value", found["p-value"], r["circ"], "R circular rayleigh.test", CLOSED,
                       "relative")


checks.data("morley", morley)
checks.data("darwin", darwin)
checks.data("mussels", mussels)
checks.data("groups3", groups3)
checks.data("houses", houses)
checks.data("ages", ages)
checks.data("pain", pain)
checks.data("sleep", sleep)
checks.data("mtcars", mtcars)
checks.data("warpbreaks", warpbreaks)
checks.data("co2", co2)
checks.data("depression", depression)
checks.data("paired", paired)
checks.data("approval", approval)
checks.data("plants", plants)
checks.data("serum", serum)
checks.data("lev", lev)
checks.data("insects", insects)
checks.data("three", three)
checks.data("earth", earth)
checks.data("potatoes", potatoes)
checks.data("wind", wind)
checks.save()
