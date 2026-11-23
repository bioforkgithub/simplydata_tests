"""Reference values for notebooks/17_bayesian.ipynb, saved to references/17_bayesian.csv.

- Printed: the Bayes factors of the BayesFactor manual (Morey & Rouder, vignette "Using the 'BayesFactor' package",
  version 0.9.12-4.8): sleep, paired, 17.259; iris, Sepal.Length with Sepal.Width, 0.509; Hraba & Grant's dolls,
  1.81. Half a unit of the last printed digit.
- Run here: R BayesFactor with the same call (a Bayes factor is worked out without random draws: 1e-9 relative);
  the published formulas in Python, independent of R (the JZS t-test integral of Rouder et al. 2009, eq. 1; the
  exact Jeffreys correlation Bayes factor of Ly, Verhagen & Wagenmakers 2016, both by SciPy quad, within 4 times
  BayesFactor's own reported integration error, at least 1e-8; the Gunel & Dickey 1974 closed form, 1e-9).
- Credible intervals come from random draws. (1) Same draws: R with the same seed (BayesFactor posterior(iterations
  = 10000) after set.seed(1); rstanarm seed = 1, 4 chains of 2000 iterations), 1e-9 relative. (2) Monte Carlo: an
  independent long run with another seed (BayesFactor 200,000 draws; rstanarm 4 chains of 10,000 iterations); the
  allowed difference is 4 * sqrt(MCSE_same^2 + MCSE_long^2), the MCSEs from posterior::mcse_median() and
  posterior::mcse_quantile() on each run's own draws (the same-seed draws are simplydata's draws). R-hat and the
  effective draws (rstanarm) against the same-seed run only (1e-9 relative).

This script runs rstanarm (six long fits): run it on GitHub Actions, not on a laptop. The tables the notebook reads
(the sleep data with its groups named, the two drugs side by side, drug 1 alone, mtcars with gearbox and engine
named) are written to references/data/17_bayesian/.

    python scripts/references/17_bayesian.py      (from the repository root; needs R with BayesFactor, rstanarm, posterior)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import integrate, special, stats  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "17_bayesian"
checks = Recorder(STEM, ("BayesFactor", "rstanarm", "posterior"))
HAS_R = r_available()
SAME_DRAWS = 1e-9
CLOSED_FORM = 1e-9
N_MCSE = 4
SEED = 'set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")'
LONG_SEED = 'set.seed(20261121, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")'
MANUAL = "BayesFactor manual (0.9.12-4.8), printed"
MANUAL_TITLE = "Morey RD, Rouder JN. Using the 'BayesFactor' package, version 0.9.12-4.8 (vignette)"
MANUAL_URL = "https://cran.r-project.org/web/packages/BayesFactor/vignettes/manual.html"

sleep = load_csv("sleep.csv")
sleep["group"] = sleep["group"].map({1: "drug 1", 2: "drug 2"})
checks.data("sleep", sleep)
drug_1 = sleep.loc[sleep["group"] == "drug 1", "extra"].to_numpy(dtype=float)
drug_2 = sleep.loc[sleep["group"] == "drug 2", "extra"].to_numpy(dtype=float)
mtcars = load_csv("mtcars.csv")
iris = load_csv("iris.csv")
dolls = load_csv("race_dolls.csv")
sleepstudy = load_csv("sleepstudy.csv")


def ref(dataset, method, quantity, value, source, tolerance, how="relative"):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, how)


def manual(dataset, method, shown):
    checks.printed(dataset, method, "Bayes factor", shown, MANUAL, title=MANUAL_TITLE, url=MANUAL_URL)


def jzs(t, n_eff, df, r):
    # Rouder et al. (2009), equation 1, with g ~ inverse gamma(1/2, r^2 / 2)
    def f(g):
        return ((1 + n_eff * g) ** -0.5 * (1 + t ** 2 / ((1 + n_eff * g) * df)) ** (-(df + 1) / 2)
                * r / np.sqrt(2 * np.pi) * g ** -1.5 * np.exp(-r ** 2 / (2 * g)))
    return integrate.quad(f, 0, np.inf, limit=500, epsabs=0, epsrel=1e-12)[0] / (1 + t ** 2 / df) ** (-(df + 1) / 2)


def jeffreys_correlation(r, n, width):
    # Ly, Verhagen & Wagenmakers (2016): exact likelihood ratio h(rho), stretched beta prior of width w on rho
    a = 1 / width
    log_norm = -special.betaln(a, a) - (2 * a - 1) * np.log(2)
    def f(rho):
        h = (1 - rho ** 2) ** ((n - 1) / 2) * special.hyp2f1((n - 1) / 2, (n - 1) / 2, 0.5, (r * rho) ** 2)
        return h * np.exp(log_norm + (a - 1) * np.log1p(-rho ** 2))
    return 2 * integrate.quad(f, 0, 1, limit=500, epsabs=0, epsrel=1e-12)[0]


def gunel_dickey(table):
    # rows: the groups (their totals fixed); columns: the values; a = 1
    table = np.asarray(table, dtype=float)
    def log_b(alpha):
        return special.gammaln(alpha).sum() - special.gammaln(alpha.sum())
    a = np.ones(table.shape[1])
    log_bf = sum(log_b(row + a) - log_b(a) for row in table) - (log_b(table.sum(axis=0) + a) - log_b(a))
    return float(np.exp(log_bf))


def integration_tolerance(error):
    return max(4 * float(error), 1e-8)


def interval_references(dataset, method, bf_code, parameter, data, transform=""):
    """The 2.5 %, 50 % and 97.5 % quantiles of the posterior: R with simplydata's seed and draws, and a long run."""
    pick = transform or f'x <- as.numeric(p[, "{parameter}"])'
    run = ("q3 <- function(x) unname(quantile(x, c(0.025, 0.5, 0.975)))\n"
           "se3 <- function(x) c(posterior::mcse_quantile(x, probs = 0.025), posterior::mcse_median(x), "
           "posterior::mcse_quantile(x, probs = 0.975))\n"
           f"{bf_code}\n{SEED}\np <- BayesFactor::posterior(b, iterations = 10000, progress = FALSE)\n"
           "q <- summary(p)$quantiles\n"
           f"{pick}\nsame <- q3(x); same_se <- se3(x)\n"
           f'same_coda <- unname(q["{parameter}", c("2.5%", "50%", "97.5%")])\n'
           f"{LONG_SEED}\np <- BayesFactor::posterior(b, iterations = 200000, progress = FALSE)\n"
           f"{pick}\nlong <- q3(x); long_se <- se3(x)\n"
           "list(same = same_coda, same_se = unname(same_se), long = long, long_se = unname(long_se))")
    r = r_run(run, data)
    for k, name in enumerate(("2.5 % quantile", "median", "97.5 % quantile")):
        ref(dataset, method, name, r["same"][k], "R BayesFactor posterior (same seed)", SAME_DRAWS)
        allowed = N_MCSE * np.hypot(r["same_se"][k], r["long_se"][k])
        ref(dataset, method, name, r["long"][k], "R BayesFactor posterior (200,000 draws, other seed)", allowed, "absolute")


# ---- two groups
t_two = stats.ttest_ind(drug_1, drug_2).statistic
n_eff = len(drug_1) * len(drug_2) / (len(drug_1) + len(drug_2))
for width, r_width, scale in (("medium", '"medium"', np.sqrt(2) / 2), ("wide", '"wide"', 1.0),
                              ("ultrawide", '"ultrawide"', np.sqrt(2))):
    dataset, method = "Student's sleep data", f"compare_averages_of (bayesian, prior {width})"
    if HAS_R:
        r = r_run(f'b <- BayesFactor::ttestBF(x = extra[group == "drug 1"], y = extra[group == "drug 2"], '
                  f'rscale = {r_width}); e <- BayesFactor::extractBF(b)\nlist(bf = e$bf, error = e$error)', sleep)
        ref(dataset, method, "Bayes factor", r["bf"], "R BayesFactor ttestBF", SAME_DRAWS)
        ref(dataset, method, "Bayes factor", jzs(t_two, n_eff, len(drug_1) + len(drug_2) - 2, scale),
            "SciPy quad of the JZS integral (Rouder et al. 2009, eq. 1)", integration_tolerance(r["error"]))
        if width == "medium":
            interval_references(dataset, method + ": difference",
                                'b <- BayesFactor::ttestBF(x = extra[group == "drug 1"], y = extra[group == "drug 2"], '
                                'rscale = "medium")', "beta (x - y)", sleep)

# ---- before and after
wide = pd.DataFrame({"drug_1": drug_1, "drug_2": drug_2})
checks.data("sleep_wide", wide)
dataset, method = "Student's sleep data", "compare_before_and_after (bayesian)"
change = drug_2 - drug_1
t_paired = change.mean() / (change.std(ddof=1) / np.sqrt(len(change)))
manual(dataset, method, "17.259")
if HAS_R:
    r = r_run('b <- BayesFactor::ttestBF(x = drug_2, y = drug_1, paired = TRUE); e <- BayesFactor::extractBF(b)\n'
              'list(bf = e$bf, error = e$error)', wide)
    ref(dataset, method, "Bayes factor", r["bf"], "R BayesFactor ttestBF", SAME_DRAWS)
    ref(dataset, method, "Bayes factor", jzs(t_paired, len(change), len(change) - 1, np.sqrt(2) / 2),
        "SciPy quad of the JZS integral (Rouder et al. 2009, eq. 1)", integration_tolerance(r["error"]))
    interval_references(dataset, method + ": average change",
                        "b <- BayesFactor::ttestBF(x = drug_2, y = drug_1, paired = TRUE)", "mu", wide)

# ---- one column against a number
one = pd.DataFrame({"extra": drug_1})
checks.data("drug_1", one)
dataset, method = "Student's sleep data (drug 1)", "compare_average_to 1 (bayesian, prior wide)"
t_one = (drug_1.mean() - 1) / (drug_1.std(ddof=1) / np.sqrt(len(drug_1)))
if HAS_R:
    r = r_run('b <- BayesFactor::ttestBF(x = extra, mu = 1, rscale = "wide"); e <- BayesFactor::extractBF(b)\n'
              'list(bf = e$bf, error = e$error)', one)
    ref(dataset, method, "Bayes factor", r["bf"], "R BayesFactor ttestBF", SAME_DRAWS)
    ref(dataset, method, "Bayes factor", jzs(t_one, len(drug_1), len(drug_1) - 1, 1.0),
        "SciPy quad of the JZS integral (Rouder et al. 2009, eq. 1)", integration_tolerance(r["error"]))
    interval_references(dataset, method + ": difference", 'b <- BayesFactor::ttestBF(x = extra - 1, rscale = "wide")',
                        "mu", one)

# ---- correlation
column = {c.lower().replace(".", "_"): c for c in iris.columns}
pairs = [("Fisher's iris", iris, column["sepal_length"], column["sepal_width"], "0.509"), ("mtcars", mtcars, "mpg", "wt", None)]
for dataset, data, y, x, printed in pairs:
    method = f"find_correlation_between {y}, {x} (bayesian)"
    if printed is not None:
        manual(dataset, method, printed)
    if HAS_R:
        pair = pd.DataFrame({"a": data[y], "b": data[x]})
        r = r_run('b <- BayesFactor::correlationBF(y = a, x = b, rscale = "medium"); e <- BayesFactor::extractBF(b)\n'
                  'list(bf = e$bf, error = e$error)', pair)
        ref(dataset, method, "Bayes factor", r["bf"], "R BayesFactor correlationBF", SAME_DRAWS)
        sample_r = np.corrcoef(data[y], data[x])[0, 1]
        ref(dataset, method, "Bayes factor", jeffreys_correlation(sample_r, len(data), 1 / 3),
            "SciPy quad of the exact Jeffreys Bayes factor (Ly et al. 2016)", integration_tolerance(r["error"]))
        interval_references(dataset, method + ": r", 'b <- BayesFactor::correlationBF(y = a, x = b, rscale = "medium")',
                            "rho", pair)

# ---- counts
cars = mtcars.assign(gearbox=mtcars["am"].map({0: "automatic", 1: "manual"}),
                     engine=mtcars["vs"].map({0: "V-shaped", 1: "straight"}))
checks.data("mtcars", cars)
tables = [("Hraba & Grant dolls", dolls, "doll", "child", "1.81"), ("mtcars", cars, "gearbox", "engine", None)]
for dataset, data, value, group, printed in tables:
    method = f"compare_counts_of {value} by {group} (bayesian)"
    groups, values = sorted(data[group].unique()), sorted(data[value].unique())
    counts = pd.crosstab(data[group], data[value]).loc[groups, values]
    if printed is not None:
        manual(dataset, method, printed)
    ref(dataset, method, "Bayes factor", gunel_dickey(counts.to_numpy()),
        "Gunel & Dickey (1974) closed form (SciPy gammaln)", CLOSED_FORM)
    if HAS_R:
        frame = pd.DataFrame({"g": data[group].astype(str), "v": data[value].astype(str)})
        levels = (f'tab <- table(factor(g, levels = c({", ".join(repr(str(x)) for x in groups)})), '
                  f'factor(v, levels = c({", ".join(repr(str(x)) for x in values)})))\n'
                  'b <- BayesFactor::contingencyTableBF(tab, sampleType = "indepMulti", fixedMargin = "rows")')
        r = r_run(levels + "\nlist(bf = BayesFactor::extractBF(b)$bf)", frame)
        ref(dataset, method, "Bayes factor", r["bf"], "R BayesFactor contingencyTableBF", SAME_DRAWS)
        for i, g in enumerate(groups, start=1):
            for j, v in enumerate(values, start=1):
                interval_references(dataset, f"{method}: share {v} in {g}", levels, f"omega[{i},{j}]", frame)

# ---- Bayesian models (rstanarm)
MODELS = [
    ("mtcars", mtcars, "mpg ~ wt + hp", ["wt", "hp"], "gaussian()", "stan_glm", None),
    ("mtcars", mtcars, "am ~ wt", ["wt"], "binomial()", "stan_glm", "odds ratio"),
    ("sleepstudy", sleepstudy, "Reaction ~ Days + (1 | Subject)", ["Days"], "gaussian()", "stan_glmer", None),
]
checks.data("sleepstudy", sleepstudy)
STAN = ("draws_of <- function(m) posterior::as_draws_array(as.array(m))\n"
        "q3 <- function(x) unname(quantile(x, c(0.025, 0.5, 0.975)))\n"
        "se3 <- function(x) c(posterior::mcse_quantile(x, probs = 0.025), posterior::mcse_median(x), "
        "posterior::mcse_quantile(x, probs = 0.975))\n")
for dataset, data, formula, using, family, fit, ratio in MODELS:
    method = f"explain {formula} (bayesian, rstanarm {fit})"
    if not HAS_R:
        continue
    columns = [c for c in data.columns if c in formula]
    frame = data[columns].copy()
    if "Subject" in frame:
        frame["Subject"] = frame["Subject"].astype(str)
    fixed_names = ["(Intercept)"] + using
    extra = ["sigma"] if family == "gaussian()" else []
    code_r = (STAN + f"d <- data.frame({', '.join(columns)})\n"
              f"same <- rstanarm::{fit}({formula}, data = d, family = {family}, chains = 4, iter = 2000, seed = 1, "
              "refresh = 0, cores = 1)\n"
              f"long <- rstanarm::{fit}({formula}, data = d, family = {family}, chains = 4, iter = 10000, "
              "seed = 20261121, refresh = 0, cores = 1)\n"
              "ds <- draws_of(same); dl <- draws_of(long)\n"
              "s <- posterior::summarise_draws(ds, 'rhat', 'ess_bulk', 'ess_tail')\n"
              f"pars <- c({', '.join(repr(p) for p in fixed_names + extra)})\n"
              "out <- list(pars = pars)\n"
              "for (k in pars) { xs <- as.numeric(posterior::extract_variable(ds, k)); "
              "xl <- as.numeric(posterior::extract_variable(dl, k)); "
              "out[[paste0(k, '|same')]] <- q3(xs); out[[paste0(k, '|same_se')]] <- se3(xs); "
              "out[[paste0(k, '|long')]] <- q3(xl); out[[paste0(k, '|long_se')]] <- se3(xl); "
              "out[[paste0(k, '|diag')]] <- c(s$rhat[s$variable == k], s$ess_bulk[s$variable == k], "
              "s$ess_tail[s$variable == k]) }\n"
              "out")
    r = r_run(code_r, frame)
    pars = r["pars"] if isinstance(r["pars"], list) else [r["pars"]]
    for par in pars:
        term = {"(Intercept)": "(intercept)", "sigma": "spread around the model (sigma)"}.get(par, par)
        convert = np.exp if (ratio and par != "sigma") else (lambda v: v)
        for k, name in enumerate(("2.5 % quantile", "median", "97.5 % quantile")):
            quantity = f"{term} {name}" + (f" ({ratio})" if ratio else "")
            ref(dataset, method, quantity, convert(r[f"{par}|same"][k]), "R rstanarm (same seed)", SAME_DRAWS)
            allowed = N_MCSE * np.hypot(r[f"{par}|same_se"][k], r[f"{par}|long_se"][k])
            if ratio:  # MCSE on the log scale, carried to the ratio scale by the delta method
                allowed *= convert(r[f"{par}|long"][k])
            ref(dataset, method, quantity, convert(r[f"{par}|long"][k]), "R rstanarm (4 x 10,000 iterations, other seed)",
                allowed, "absolute")
        for k, name in enumerate(("R-hat", "effective draws (middle)", "effective draws (tails)")):
            ref(dataset, method, f"{term} {name}", r[f"{par}|diag"][k], "R rstanarm (same seed)", SAME_DRAWS)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
