"""Reference values for notebooks/23_distributions.ipynb, saved to references/23_distributions.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/23_distributions/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/23_distributions.py      (from the repository root; needs R with the packages listed)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks and realdata)

import numpy as np
import pandas as pd
import simplydata as sdata
from scipy import optimize, special, stats

import realdata as rd
from checks import Checks, quiet, r_available, r_run, read_numbers, read_row, read_rows_by, versions

versions("pandas", "numpy", "scipy", "rpy2")
checks = Recorder("23_distributions")
HAS_R = r_available()
EXACT = 1e-12
QUANTILE = 1e-10
ITERATIVE = 1e-6
if HAS_R:
    print(rd.r_version_of("MASS", "fitdistrplus", "circular"))


G = special.gamma
w_k, w_s = 1.5, 2.0
wm = w_s * G(1 + 1 / w_k)
wv = w_s ** 2 * (G(1 + 2 / w_k) - G(1 + 1 / w_k) ** 2)
w_skew = (w_s ** 3 * G(1 + 3 / w_k) - 3 * wm * wv - wm ** 3) / wv ** 1.5
w_kurt = (w_s ** 4 * G(1 + 4 / w_k) - 4 * w_skew * wv ** 1.5 * wm - 6 * wm ** 2 * wv - wm ** 4) / wv ** 2 - 3
ln_m, ln_s = 1.0, 0.5
e2 = np.exp(ln_s ** 2)
d1, d2 = 5, 12
f_var = 2 * d2 ** 2 * (d1 + d2 - 2) / (d1 * (d2 - 2) ** 2 * (d2 - 4))
f_skew = (2 * d1 + d2 - 2) * np.sqrt(8 * (d2 - 4)) / ((d2 - 6) * np.sqrt(d1 * (d1 + d2 - 2)))
f_kurt = 12 * (d1 * (5 * d2 - 22) * (d1 + d2 - 2) + (d2 - 4) * (d2 - 2) ** 2) / (d1 * (d2 - 6) * (d2 - 8) * (d1 + d2 - 2))
a_, b_ = 2.0, 5.0
beta_var = a_ * b_ / ((a_ + b_) ** 2 * (a_ + b_ + 1))
beta_skew = 2 * (b_ - a_) * np.sqrt(a_ + b_ + 1) / ((a_ + b_ + 2) * np.sqrt(a_ * b_))
beta_kurt = 6 * ((a_ - b_) ** 2 * (a_ + b_ + 1) - a_ * b_ * (a_ + b_ + 2)) / (a_ * b_ * (a_ + b_ + 2) * (a_ + b_ + 3))
nb_r, nb_p = 3, 0.4
CASES = [
    # name, parameters, R d, R q, values, mean, variance, skewness, excess kurtosis, counts?
    ("normal", dict(average=50, standard_deviation=10), "dnorm({x}, 50, 10)", "qnorm({p}, 50, 10)", [35, 50, 61.5],
     50, 100, 0, 0, False),
    ("lognormal", dict(meanlog=ln_m, sdlog=ln_s), "dlnorm({x}, 1, 0.5)", "qlnorm({p}, 1, 0.5)", [0.5, 2.7, 6],
     np.exp(ln_m + ln_s ** 2 / 2), (e2 - 1) * np.exp(2 * ln_m + ln_s ** 2), (e2 + 2) * np.sqrt(e2 - 1),
     e2 ** 4 + 2 * e2 ** 3 + 3 * e2 ** 2 - 6, False),
    ("exponential", dict(rate=0.5), "dexp({x}, 0.5)", "qexp({p}, 0.5)", [0.1, 2, 7], 2, 4, 2, 6, False),
    ("gamma", dict(shape=2, rate=0.5), "dgamma({x}, shape = 2, rate = 0.5)", "qgamma({p}, shape = 2, rate = 0.5)",
     [0.3, 4, 11], 4, 8, 2 / np.sqrt(2), 3, False),
    ("beta", dict(shape1=a_, shape2=b_), "dbeta({x}, 2, 5)", "qbeta({p}, 2, 5)", [0.05, 0.3, 0.8], a_ / (a_ + b_),
     beta_var, beta_skew, beta_kurt, False),
    ("uniform", dict(min=2, max=5), "dunif({x}, 2, 5)", "qunif({p}, 2, 5)", [2.5, 3.5, 4.9], 3.5, 9 / 12, 0, -1.2, False),
    ("t", dict(degrees_of_freedom=7), "dt({x}, 7)", "qt({p}, 7)", [-2.5, 0, 1.3], 0, 7 / 5, 0, 6 / 3, False),
    ("chi-square", dict(degrees_of_freedom=4), "dchisq({x}, 4)", "qchisq({p}, 4)", [0.5, 3, 12], 4, 8, np.sqrt(2), 3,
     False),
    ("F", dict(df1=d1, df2=d2), "df({x}, 5, 12)", "qf({p}, 5, 12)", [0.2, 1, 4], d2 / (d2 - 2), f_var, f_skew, f_kurt,
     False),
    ("weibull", dict(shape=w_k, scale=w_s), "dweibull({x}, 1.5, 2)", "qweibull({p}, 1.5, 2)", [0.2, 1.7, 5], wm, wv,
     w_skew, w_kurt, False),
    ("logistic", dict(location=1, scale=2), "dlogis({x}, 1, 2)", "qlogis({p}, 1, 2)", [-4, 1, 6], 1, 4 * np.pi ** 2 / 3,
     0, 1.2, False),
    ("binomial", dict(tries=10, chance=0.3), "dbinom({x}, 10, 0.3)", "qbinom({p}, 10, 0.3)", [0, 3, 7], 3, 2.1,
     (1 - 0.6) / np.sqrt(2.1), (1 - 6 * 0.21) / 2.1, True),
    ("poisson", dict(average=3), "dpois({x}, 3)", "qpois({p}, 3)", [0, 3, 8], 3, 3, 1 / np.sqrt(3), 1 / 3, True),
    ("negative binomial", dict(size=nb_r, prob=nb_p), "dnbinom({x}, 3, 0.4)", "qnbinom({p}, 3, 0.4)", [0, 4, 12],
     nb_r * (1 - nb_p) / nb_p, nb_r * (1 - nb_p) / nb_p ** 2, (2 - nb_p) / np.sqrt(nb_r * (1 - nb_p)),
     6 / nb_r + nb_p ** 2 / (nb_r * (1 - nb_p)), True),
    ("geometric", dict(chance=0.25), "dgeom({x}, 0.25)", "qgeom({p}, 0.25)", [0, 2, 9], 0.75 / 0.25, 0.75 / 0.0625,
     (2 - 0.25) / np.sqrt(0.75), 6 + 0.0625 / 0.75, True),
]
PERCENTS = [5, 50, 95]

for name, params, r_d, r_q, values, mean, var, skew, kurt, counts in CASES:
    dist = sdata.distribution(name, **params)
    ds = f"{name} " + ", ".join(f"{k}={v:g}" for k, v in params.items())
    shown = dist.show_parameters()
    print(ds, shown)
    summary = read_row(quiet(dist.show_summary))
    for quantity, ours, ref in [("average", quiet(dist.find_average), mean),
                                ("variance", quiet(dist.find_variance), var),
                                ("standard deviation", quiet(dist.find_standard_deviation), np.sqrt(var)),
                                ("show_summary average", summary["average"], mean),
                                ("show_summary variance", summary["variance"], var),
                                ("show_summary skewness", summary["skewness"], skew),
                                ("show_summary excess kurtosis", summary["excess kurtosis"], kurt)]:
        word = quantity.split(" ")[0] if quantity.startswith("show_summary") else {
            "average": "find_average", "variance": "find_variance", "standard deviation": "find_standard_deviation"}[quantity]
        how = "absolute" if ref == 0 else "relative"
        checks.compare(ds, f"dist.{word}", quantity, ours, ref, "textbook formula (Johnson, Kotz & Balakrishnan)", EXACT, how)
    if not HAS_R:
        continue
    r = r_run("list(" + ", ".join(f"d{i} = {r_d.format(x=x)}" for i, x in enumerate(values)) + ", "
              + ", ".join(f"q{p} = {r_q.format(p=p / 100)}" for p in PERCENTS) + ", "
              + f"c20 = {r_q.format(p=0.2)}, c80 = {r_q.format(p=0.8)}, m05 = {r_q.format(p=0.05)}, m95 = {r_q.format(p=0.95)}, "
              + f"m025 = {r_q.format(p=0.025)}, m975 = {r_q.format(p=0.975)})")
    for i, x in enumerate(values):
        ours = quiet(dist.find_density_at, x) if not counts else quiet(dist.find_chance_that_value, x)
        word = "dist.find_density_at" if not counts else "dist.find_chance_that_value (a count)"
        checks.compare(ds, word, f"at {x:g}", ours, r[f"d{i}"], "R " + r_d.split("(")[0], EXACT, "relative")
    found = quiet(dist.find_percentile, PERCENTS)
    for p, ours in zip(PERCENTS, read_numbers(found)):
        checks.compare(ds, "dist.find_percentile (percent=)", f"percentile {p}", ours, r[f"q{p}"], "R " + r_q.split("(")[0],
                       QUANTILE, "relative")
    checks.compare(ds, "dist.find_median", "median", quiet(dist.find_median), r["q50"], "R " + r_q.split("(")[0], QUANTILE,
                   "relative")
    checks.compare(ds, "dist.find_cut_off (with_chance_below=)", "cut-off with 20 % below",
                   quiet(dist.find_cut_off, with_chance_below=0.2), r["c20"], "R " + r_q.split("(")[0], QUANTILE, "relative")
    if not counts:
        checks.compare(ds, "dist.find_cut_off (with_chance_above=)", "cut-off with 20 % above",
                       quiet(dist.find_cut_off, with_chance_above=0.2), r["c80"], "R " + r_q.split("(")[0], QUANTILE,
                       "relative")
    for percent, keys in [(90, ("m05", "m95")), (95, ("m025", "m975"))]:
        low, high = read_numbers(quiet(dist.find_middle_range, percent))
        label = "dist.find_middle_range" + (" (percent=90)" if percent == 90 else "")
        checks.compare(ds, label, f"middle {percent} % low", low, r[keys[0]], "R " + r_q.split("(")[0], QUANTILE, "relative")
        checks.compare(ds, label, f"middle {percent} % high", high, r[keys[1]], "R " + r_q.split("(")[0], QUANTILE, "relative")


vm = sdata.distribution("von mises", average=1, concentration=2)
ds = "von mises average=1, concentration=2"
print(vm.show_parameters())
if HAS_R:
    for x in (0.0, 1.0, 2.5, 5.0):
        r = r_run(f"""suppressPackageStartupMessages(library(circular))
            list(d = dvonmises(circular({x}), mu = circular(1), kappa = 2))""")
        checks.compare(ds, "dist.find_density_at", f"at {x:g} radians", quiet(vm.find_density_at, x), r["d"],
                       "R circular dvonmises", EXACT, "relative")
checks.compare(ds, "dist.find_density_at", "at 1 radian", quiet(vm.find_density_at, 1.0),
               np.exp(2) / (2 * np.pi * special.i0(2)), "exp(k) / (2 pi I0(k))", EXACT, "relative")


PLAIN = {"meanlog": "log_average", "sdlog": "log_standard_deviation", "shape1": "shape_1", "shape2": "shape_2",
         "df1": "degrees_of_freedom_1", "df2": "degrees_of_freedom_2", "size": "size", "prob": "chance",
         "min": "lowest", "max": "highest"}
for name, params, r_d, r_q, values, mean, var, skew, kurt, counts in CASES:
    dist = sdata.distribution(name, **params)
    ds = f"{name} " + ", ".join(f"{k}={v:g}" for k, v in params.items())
    shown = dist.show_parameters()
    given = list(params.values())
    print(ds, "->", shown)
    for key, value in params.items():
        checks.compare(ds, "dist.show_parameters", f"{key} is shown (as {PLAIN.get(key, key)} or its plain name)",
                       float(any(np.isclose(v, value, rtol=0, atol=0) for v in shown.values())), 1.0, "the value given", 0)
    draws = np.asarray(read_numbers(quiet(dist.make_random_numbers, how_many=1000, with_seed=3)), dtype=float)
    frozen = dist._frozen if hasattr(dist, "_frozen") else None
    if frozen is not None:
        parent = frozen.rvs(size=1000, random_state=np.random.default_rng(3))
        checks.vector(ds, "dist.make_random_numbers (how_many=, with_seed=)", "1000 draws", parent,
                      "SciPy rvs, same seed", 1e-15)
    big = np.asarray(read_numbers(quiet(dist.make_random_numbers, how_many=400_000, with_seed=4)), dtype=float)
    checks.compare(ds, "dist.make_random_numbers", "average of 400 000 draws", big.mean(), mean,
                   "true mean (4 MCSE)", 4 * np.sqrt(var / 400_000))

# dist.draw_chart: the curve simplydata draws, read in the notebook with chart.get_numbers_drawn(). The reference is
# R's density (d...) or distribution function (p...) at the x of each drawn point (the x are read from the drawn
# matplotlib line; only the y is compared).
R_PAGES = {"normal": "Normal", "lognormal": "Lognormal", "exponential": "Exponential", "gamma": "GammaDist",
           "beta": "Beta", "uniform": "Uniform", "t": "TDist", "chi-square": "Chisquare", "F": "Fdist",
           "weibull": "Weibull", "logistic": "Logistic"}
R_VERSION = r_run("list(v = R.version.string)")["v"][0] if HAS_R else ""
for name, params, r_d, r_q, values, mean, var, skew, kurt, counts in CASES:
    if counts:
        continue
    ds = f"{name} " + ", ".join(f"{k}={v:g}" for k, v in params.items())
    dist = sdata.distribution(name, **params)
    for show in ("density", "cumulative"):
        method = "dist.draw_chart" + ("" if show == "density" else ' (show="cumulative")')
        r_f = r_d if show == "density" else r_d.replace("d", "p", 1)
        source = "R " + r_f.split("(")[0]
        if not HAS_R:
            checks.skip(ds, method, "the drawn curve", source, "R is not available")
            continue
        chart = quiet(dist.draw_chart) if show == "density" else quiet(dist.draw_chart, show="cumulative")
        fig = quiet(lambda: chart)._render()
        line = max(fig.axes[0].get_lines(), key=lambda l: len(l.get_xdata()))
        x = np.asarray(line.get_xdata(), dtype=float)
        reference = np.asarray(r_run("list(y = " + r_f.format(x="x") + ")", data={"x": x})["y"], dtype=float)
        checks.vector(ds, method, "the drawn curve", reference, source, QUANTILE)
        checks.rows[-1]["url"] = (f"https://search.r-project.org/R/refmans/stats/html/{R_PAGES[name]}.html"
                                  f" ({R_VERSION.strip()})")


airquality = rd.r_dataset("airquality")
ozone = airquality["Ozone"].dropna().to_numpy()
air = sdata.table(airquality)
n = len(ozone)
lx = np.log(ozone)


def loglik(dist_scipy):
    return dist_scipy.logpdf(ozone).sum()


k_hat = optimize.brentq(lambda k: np.log(k) - special.digamma(k) - (np.log(ozone.mean()) - lx.mean()), 1e-3, 1e3,
                        xtol=1e-15, rtol=1e-15)
b_hat = k_hat / ozone.mean()
info = n * np.array([[special.polygamma(1, k_hat), -1 / b_hat], [-1 / b_hat, k_hat / b_hat ** 2]])
se_gamma = np.sqrt(np.diag(np.linalg.inv(info)))


def weibull_shape_equation(k):
    xk = ozone ** k
    return (xk * np.log(ozone)).sum() / xk.sum() - 1 / k - lx.mean()


wk = optimize.brentq(weibull_shape_equation, 0.1, 20, xtol=1e-15, rtol=1e-15)
wsc = (np.mean(ozone ** wk)) ** (1 / wk)
mu, sd = ozone.mean(), ozone.std(ddof=0)
lmu, lsd = lx.mean(), lx.std(ddof=0)
FITS = {
    "normal": ({"average": mu, "standard_deviation": sd}, {"average": sd / np.sqrt(n), "standard_deviation": sd / np.sqrt(2 * n)},
               stats.norm(mu, sd), 2, "normal"),
    "lognormal": ({"log_average": lmu, "log_standard_deviation": lsd},
                  {"log_average": lsd / np.sqrt(n), "log_standard_deviation": lsd / np.sqrt(2 * n)},
                  stats.lognorm(lsd, scale=np.exp(lmu)), 2, "lognormal"),
    "exponential": ({"rate": 1 / mu}, {"rate": (1 / mu) / np.sqrt(n)}, stats.expon(scale=mu), 1, "exponential"),
    "gamma": ({"shape": k_hat, "rate": b_hat}, {"shape": se_gamma[0], "rate": se_gamma[1]},
              stats.gamma(k_hat, scale=1 / b_hat), 2, "gamma"),
    "weibull": ({"shape": wk, "scale": wsc}, {}, stats.weibull_min(wk, scale=wsc), 2, "weibull"),
}
ds = "airquality ozone (Chambers et al. 1983)"
for name, (estimates, ses, frozen, k_params, r_name) in FITS.items():
    fit = quiet(air.fit_distribution_to, "Ozone", name)
    row = read_row(quiet(fit.turn_into_table))
    for p_name, value in estimates.items():
        column = p_name if p_name in row else p_name.replace("_", " ")
        checks.compare(ds, f'table.fit_distribution_to ("{name}")', p_name, row[column], value,
                       "likelihood equations solved in SciPy" if name in ("gamma", "weibull") else "closed-form MLE",
                       ITERATIVE if name in ("gamma", "weibull") else EXACT, "relative")
    for p_name, value in ses.items():
        column = f"standard error of {p_name if p_name in row else p_name.replace('_', ' ')}"
        checks.compare(ds, f'table.fit_distribution_to ("{name}")', f"standard error of {p_name}", row[column], value,
                       "observed information, formula", ITERATIVE, "relative")
    ll = loglik(frozen)
    checks.compare(ds, f'table.fit_distribution_to ("{name}")', "log-likelihood", row["log-likelihood"], ll,
                   "at the reference estimates", ITERATIVE, "relative")
    checks.compare(ds, f'table.fit_distribution_to ("{name}")', "AIC", row["AIC"], 2 * k_params - 2 * ll,
                   "2k - 2 log L", ITERATIVE, "relative")
    checks.compare(ds, f'table.fit_distribution_to ("{name}")', "BIC", row["BIC"], k_params * np.log(n) - 2 * ll,
                   "k log n - 2 log L", ITERATIVE, "relative")
    via_numbers = quiet(sdata.numbers(ozone.tolist()).fit_distribution, distribution=name).show_parameters()
    for p_name, value in estimates.items():
        key = p_name if p_name in via_numbers else p_name.replace("_", " ")
        checks.compare(ds, f'nums.fit_distribution (distribution="{name}")', p_name, via_numbers[key], value,
                       "likelihood equations solved in SciPy" if name in ("gamma", "weibull") else "closed-form MLE",
                       ITERATIVE if name in ("gamma", "weibull") else EXACT, "relative")
    if HAS_R:
        r = r_run(f"""f <- fitdistrplus::fitdist(x, "{ {'normal': 'norm', 'lognormal': 'lnorm', 'exponential': 'exp', 'gamma': 'gamma', 'weibull': 'weibull'}[name] }",
                      control = list(reltol = 1e-14)); list(e = unname(f$estimate), ll = f$loglik)""", {"x": ozone})
        checks.compare(ds, f'table.fit_distribution_to ("{name}")', "log-likelihood", row["log-likelihood"], r["ll"],
                       "R fitdistrplus fitdist (reltol 1e-14)", ITERATIVE, "relative")
        for p_name, value in zip(estimates, np.atleast_1d(r["e"])):
            column = p_name if p_name in row else p_name.replace("_", " ")
            checks.compare(ds, f'table.fit_distribution_to ("{name}")', p_name, row[column], value,
                           "R fitdistrplus fitdist (reltol 1e-14)", ITERATIVE, "relative")

for word, call in [("table.fit_distribution_to", lambda: quiet(air.fit_distribution_to, "Ozone", "gamma", fixed={"shape": 2})),
                   ("nums.fit_distribution", lambda: quiet(sdata.numbers(ozone.tolist()).fit_distribution, "gamma",
                                                           fixed={"shape": 2}))]:
    row = read_row(quiet(call().turn_into_table))
    rate = 2 / ozone.mean()
    checks.compare(ds, f'{word} ("gamma", fixed={{"shape": 2}})', "rate", row["rate"], rate, "2 / mean (closed form)",
                   ITERATIVE, "relative")
    checks.compare(ds, f'{word} ("gamma", fixed={{"shape": 2}})', "standard error of rate", row["standard error of rate"],
                   rate / np.sqrt(2 * n), "rate / sqrt(2 n)", ITERATIVE, "relative")
    checks.compare(ds, f'{word} ("gamma", fixed={{"shape": 2}})', "log-likelihood", row["log-likelihood"],
                   stats.gamma(2, scale=1 / rate).logpdf(ozone).sum(), "at the closed-form rate", ITERATIVE, "relative")

discoveries = rd.r_dataset("discoveries")
discoveries_counts = discoveries.rename(columns={discoveries.columns[-1]: "count"})[["count"]]
counts = discoveries.iloc[:, -1].to_numpy(dtype=float)
ds = "discoveries 1860-1959 (R datasets)"
inventions = sdata.numbers(counts.tolist())
row = read_row(quiet(quiet(inventions.fit_distribution, distribution="poisson").turn_into_table))
print(row)
pcol = [c for c in row if c in ("average", "lambda", "rate")][0]
checks.compare(ds, 'nums.fit_distribution (distribution="poisson")', "average", row[pcol], counts.mean(), "the mean (closed form)",
               EXACT, "relative")
checks.compare(ds, 'nums.fit_distribution (distribution="poisson")', "standard error", row[f"standard error of {pcol}"],
               np.sqrt(counts.mean() / len(counts)), "sqrt(mean / n)", ITERATIVE, "relative")
if HAS_R:
    found = read_row(quiet(quiet(inventions.fit_distribution, distribution="negative binomial").turn_into_table))
    print(found)
    r = r_run("""f <- suppressWarnings(MASS::fitdistr(x, "negative binomial", control = list(reltol = 1e-14)))
        list(size = f$estimate[["size"]], mu = f$estimate[["mu"]], ll = f$loglik)""", {"x": counts})
    checks.compare(ds, 'nums.fit_distribution (distribution="negative binomial")', "log-likelihood", found["log-likelihood"],
                   r["ll"], "R MASS fitdistr (reltol 1e-14)", ITERATIVE, "relative")
    checks.compare(ds, 'nums.fit_distribution (distribution="negative binomial")', "size", found["size"], r["size"],
                   "R MASS fitdistr (reltol 1e-14)", ITERATIVE, "relative")


found = read_rows_by(quiet(air.compare_distributions_for, "Ozone"), "distribution")
print(found)
x_sorted = np.sort(ozone)
i = np.arange(1, n + 1)
for name in found:
    frozen = FITS[name][2]
    F = frozen.cdf(x_sorted)
    ks = max(np.max(i / n - F), np.max(F - (i - 1) / n))
    cvm = 1 / (12 * n) + np.sum((F - (2 * i - 1) / (2 * n)) ** 2)
    ad = -n - np.mean((2 * i - 1) * (np.log(F) + np.log(1 - F[::-1])))
    ll = loglik(frozen)
    k_params = FITS[name][3]
    for quantity, value in [("log-likelihood", ll), ("AIC", 2 * k_params - 2 * ll), ("BIC", k_params * np.log(n) - 2 * ll),
                            ("Kolmogorov-Smirnov statistic", ks), ("Cramer-von Mises statistic", cvm),
                            ("Anderson-Darling statistic", ad)]:
        checks.compare(ds := "airquality ozone (Chambers et al. 1983)", "table.compare_distributions_for", f"{name}: {quantity}",
                       found[name][quantity], value, "formula at the reference fit, NumPy", ITERATIVE, "relative")
if HAS_R:
    r = r_run("""fits <- lapply(c("gamma", "weibull", "lnorm", "norm"), function(d) fitdistrplus::fitdist(x, d, control = list(reltol = 1e-14)))
        g <- fitdistrplus::gofstat(fits); list(ad = unname(g$ad), ks = unname(g$ks), cvm = unname(g$cvm))""", {"x": ozone})
    for j, name in enumerate(["gamma", "weibull", "lognormal", "normal"]):
        for quantity, key in [("Anderson-Darling statistic", "ad"), ("Kolmogorov-Smirnov statistic", "ks"),
                              ("Cramer-von Mises statistic", "cvm")]:
            checks.compare(ds, "table.compare_distributions_for", f"{name}: {quantity}", found[name][quantity], r[key][j],
                           "R fitdistrplus gofstat (reltol 1e-14)", ITERATIVE, "relative")


checks.data("airquality", airquality)
checks.data("discoveries_counts", discoveries_counts)
checks.save()
