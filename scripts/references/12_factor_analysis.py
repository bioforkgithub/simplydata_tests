"""Reference values for notebooks/12_factor_analysis.ipynb, saved to references/12_factor_analysis.csv.

Printed first: the fit of the three-factor model of Holzinger & Swineford's nine tests as Rosseel (2012, Journal of
Statistical Software 48(2), lavaan paper) prints it (chi-square 85.306 on 24 df, CFI 0.931, TLI 0.896, RMSEA 0.092
with 90 % CI 0.071 to 0.114, SRMR 0.065), within half a unit of the last printed digit. Everything else is the
original library run once: R lavaan (cfa, sem, fitMeasures, standardizedSolution), psych (fa, fa.parallel, alpha,
alpha.ci) with GPArotation, base R factanal, statsmodels Factor(method="ml"), and the closed formulas for Cronbach's
alpha (Cronbach 1951) and Feldt's (1965) interval (NumPy and SciPy).

Tolerances: 1e-6 relative for an iterative fit against the same program run the same way (lavaan, psych); 1e-12 for
alpha against its formula; 1e-10 for Feldt's interval; 1e-3 for maximum likelihood fitted by a different optimiser
(factanal, statsmodels). psych's parallel analysis and oblimin use random numbers: set.seed(1) first, as simplydata.

The loadings are stored as psych prints them (its factor order and signs): simplydata calls psych::fa itself and
gives the factors in the same order and with the same signs, so no alignment is needed.

Data written for the notebook (references/data/12_factor_analysis/): the nine tests, and the same with the first 30
scores of x1 emptied.

    python scripts/references/12_factor_analysis.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.multivariate.factor import Factor  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "12_factor_analysis"
checks = Recorder(STEM, ("lavaan", "psych", "GPArotation"))
HAS_R = r_available()
ITERATIVE = 1e-6
FORMULA = 1e-12
QUANTILES = 1e-10
OTHER_OPTIMISER = 1e-3
ROSSEEL = "https://doi.org/10.18637/jss.v048.i02"

hs = load_csv("holzinger_swineford_1939.csv")
TESTS = [f"x{i}" for i in range(1, 10)]
MEASURED_BY = {"visual": ["x1", "x2", "x3"], "textual": ["x4", "x5", "x6"], "speed": ["x7", "x8", "x9"]}
FRAME = "d <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9)\n"
MODEL = "model <- 'visual =~ x1 + x2 + x3; textual =~ x4 + x5 + x6; speed =~ x7 + x8 + x9'\n"
DATASET = "Holzinger & Swineford 1939"
checks.data("holzinger_swineford", hs)
emptied = hs.copy()
emptied.loc[emptied.index[:30], "x1"] = np.nan
checks.data("holzinger_swineford_x1_emptied", emptied)

# ---- the three-factor model
method = "check_hidden_factors (measured_by=)"
printed = {"chi-square": "85.306", "degrees of freedom": "24", "CFI": "0.931", "TLI": "0.896", "RMSEA": "0.092",
           "RMSEA 90 % CI low": "0.071", "RMSEA 90 % CI high": "0.114", "SRMR": "0.065"}
for quantity, shown in printed.items():
    checks.compare(DATASET, method, quantity, None, float(shown), "Rosseel (2012), printed",
                   0 if quantity == "degrees of freedom" else 5e-4)
    checks.rows[-1].update({"digits as printed": shown, "url": ROSSEEL,
                            "source title": "Rosseel (2012), lavaan: An R package for structural equation modeling, "
                                            "J. Stat. Softw. 48(2)"})
if HAS_R:
    r = r_run(FRAME + MODEL + "f <- lavaan::cfa(model, data = d, estimator = 'ML'); m <- lavaan::fitMeasures(f)\n"
              "s <- lavaan::standardizedSolution(f)\n"
              "list(chisq = m[['chisq']], df = m[['df']], cfi = m[['cfi']], tli = m[['tli']], rmsea = m[['rmsea']], "
              "low = m[['rmsea.ci.lower']], high = m[['rmsea.ci.upper']], srmr = m[['srmr']], p = m[['pvalue']], "
              "lhs = s$lhs, op = s$op, rhs = s$rhs, est = s$est.std, se = s$se)", hs[TESTS])
    for quantity, key in (("chi-square", "chisq"), ("degrees of freedom", "df"), ("CFI", "cfi"), ("TLI", "tli"),
                          ("RMSEA", "rmsea"), ("RMSEA 90 % CI low", "low"), ("RMSEA 90 % CI high", "high"),
                          ("SRMR", "srmr")):
        checks.compare(DATASET, method, quantity, None, r[key], "R lavaan::fitMeasures", ITERATIVE, "relative")
    checks.compare(DATASET, method, "chi-square p-value", None, r["p"], "R lavaan::fitMeasures", ITERATIVE, "relative")
    solution = pd.DataFrame({k: r[k] for k in ("lhs", "op", "rhs", "est", "se")})
    for _, s in solution.iterrows():
        if s["op"] == "=~":
            term = f"{s['lhs']} measured by {s['rhs']}"
        elif s["op"] == "~~" and s["lhs"] != s["rhs"] and s["lhs"] in MEASURED_BY:
            term = f"{s['lhs']} with {s['rhs']}"
        else:
            continue
        checks.compare(DATASET, method, f"{term} (standardised)", None, s["est"], "R lavaan::standardizedSolution",
                       ITERATIVE, "relative")
        checks.compare(DATASET, method, f"{term} standard error", None, s["se"], "R lavaan::standardizedSolution",
                       ITERATIVE, "relative")

    # ---- a structural path, and empty cells
    method = "check_hidden_factors (formula=)"
    r = r_run(FRAME + "f <- lavaan::sem('visual =~ x1 + x2 + x3; speed =~ x7 + x8 + x9; speed ~ visual', data = d, "
              "estimator = 'ML')\n"
              "s <- lavaan::standardizedSolution(f); s <- s[s$op == '~', ]; m <- lavaan::fitMeasures(f)\n"
              "list(path = s$est.std, se = s$se, chisq = m[['chisq']], cfi = m[['cfi']], rmsea = m[['rmsea']])",
              hs[TESTS])
    checks.compare(DATASET, method, "speed explained by visual", None, r["path"], "R lavaan::sem", ITERATIVE, "relative")
    checks.compare(DATASET, method, "its standard error", None, r["se"], "R lavaan::sem", ITERATIVE, "relative")
    for quantity, key in (("chi-square", "chisq"), ("CFI", "cfi"), ("RMSEA", "rmsea")):
        checks.compare(DATASET, method, quantity, None, r[key], "R lavaan::sem", ITERATIVE, "relative")

    method = "check_hidden_factors (empty cells, FIML)"
    ds = DATASET + ", 30 cells of x1 emptied"
    r = r_run(FRAME + MODEL + "d$x1[1:30] <- NA\n"
              "f <- lavaan::cfa(model, data = d, estimator = 'ML', missing = 'ml')\n"
              "s <- lavaan::standardizedSolution(f); m <- lavaan::fitMeasures(f)\n"
              "list(chisq = m[['chisq']], cfi = m[['cfi']], rmsea = m[['rmsea']], x1 = s$est.std[1], x9 = s$est.std[9])",
              hs[TESTS])
    for quantity, key in (("chi-square", "chisq"), ("CFI", "cfi"), ("RMSEA", "rmsea")):
        checks.compare(ds, method, quantity, None, r[key], "R lavaan::cfa(missing = 'ml')", ITERATIVE, "relative")
    checks.compare(ds, method, "visual measured by x1 (standardised)", None, r["x1"], "R lavaan::cfa(missing = 'ml')",
                   ITERATIVE, "relative")
    checks.compare(ds, method, "speed measured by x9 (standardised)", None, r["x9"], "R lavaan::cfa(missing = 'ml')",
                   ITERATIVE, "relative")


# ---- exploratory factor analysis
def efa(label, r_code):
    r = r_run(FRAME + r_code + "\nL <- unclass(f$loadings)\n"
              "list(k = ncol(L), L = as.numeric(L), h2 = as.numeric(f$communality), u2 = as.numeric(f$uniquenesses))",
              hs[TESTS])
    method = f"find_hidden_factors_behind ({label})"
    k = int(r["k"])
    checks.compare(DATASET, method, "number of factors", None, k, "R psych", 0)
    loadings = np.asarray(r["L"]).reshape(len(TESTS), k, order="F")
    for j in range(k):
        for i, test in enumerate(TESTS):
            checks.compare(DATASET, method, f"{test} on factor {j + 1}", None, loadings[i, j], "R psych::fa", ITERATIVE,
                           "relative")
    for i, test in enumerate(TESTS):
        checks.compare(DATASET, method, f"shared part of {test}", None, r["h2"][i], "R psych::fa", ITERATIVE, "relative")
        checks.compare(DATASET, method, f"own part of {test}", None, r["u2"][i], "R psych::fa", ITERATIVE, "relative")


if HAS_R:
    efa("defaults: parallel analysis, minres, oblimin",
        "set.seed(1); pa <- suppressWarnings(psych::fa.parallel(d, fm = 'minres', fa = 'fa', plot = FALSE))\n"
        "f <- psych::fa(d, nfactors = pa$nfact, fm = 'minres', rotate = 'oblimin')")
    efa("3 factors, principal axis, promax", "set.seed(1); f <- psych::fa(d, nfactors = 3, fm = 'pa', rotate = 'promax')")
    efa("3 factors, maximum likelihood, varimax",
        "set.seed(1); f <- psych::fa(d, nfactors = 3, fm = 'ml', rotate = 'varimax')")

method = "find_hidden_factors_behind (3 factors, maximum likelihood, varimax)"
statsmodels_fit = Factor(hs[TESTS].to_numpy(dtype=float), n_factor=3, method="ml").fit()
for i, test in enumerate(TESTS):
    checks.compare(DATASET, method, f"own part of {test}", None, statsmodels_fit.uniqueness[i], "statsmodels Factor (ml)",
                   OTHER_OPTIMISER, "relative")
if HAS_R:
    r = r_run(FRAME + "f <- stats::factanal(d, factors = 3, rotation = 'varimax'); list(u = as.numeric(f$uniquenesses))",
              hs[TESTS])
    for i, test in enumerate(TESTS):
        checks.compare(DATASET, method, f"own part of {test}", None, r["u"][i], "R stats::factanal", OTHER_OPTIMISER,
                       "relative")


# ---- Cronbach's alpha with Feldt's interval
def cronbach(values):
    k = values.shape[1]
    return k / (k - 1) * (1 - values.var(axis=0, ddof=1).sum() / values.sum(axis=1).var(ddof=1))


for factor, items in MEASURED_BY.items():
    method = "reliability_of"
    dataset = f"{DATASET}, {factor} tests"
    values = hs[items].to_numpy(dtype=float)
    n, k = values.shape
    alpha = cronbach(values)
    checks.compare(dataset, method, "Cronbach's alpha", None, alpha, "Cronbach (1951) formula", FORMULA, "relative")
    low = 1 - (1 - alpha) * stats.f.ppf(0.975, n - 1, (n - 1) * (k - 1))
    high = 1 - (1 - alpha) * stats.f.ppf(0.025, n - 1, (n - 1) * (k - 1))
    checks.compare(dataset, method, "95 % CI low", None, low, "Feldt (1965) formula (SciPy F)", QUANTILES, "relative")
    checks.compare(dataset, method, "95 % CI high", None, high, "Feldt (1965) formula (SciPy F)", QUANTILES, "relative")
    for j, item in enumerate(items):
        checks.compare(dataset, method, f"alpha without {item}", None, cronbach(np.delete(values, j, axis=1)),
                       "Cronbach (1951) formula", FORMULA, "relative")
    if HAS_R:
        r = r_run("x <- data.frame(" + ", ".join(items) + ")\n"
                  "a <- psych::alpha(x, check.keys = FALSE)\n"
                  "f <- psych::alpha.ci(a$total$raw_alpha, n.obs = nrow(x), n.var = ncol(x), p.val = 0.05)\n"
                  "list(alpha = a$total$raw_alpha, low = f$lower.ci, high = f$upper.ci, drop = a$alpha.drop$raw_alpha)",
                  hs[items])
        checks.compare(dataset, method, "Cronbach's alpha", None, r["alpha"], "R psych::alpha", FORMULA, "relative")
        checks.compare(dataset, method, "95 % CI low", None, r["low"], "R psych::alpha.ci", QUANTILES, "relative")
        checks.compare(dataset, method, "95 % CI high", None, r["high"], "R psych::alpha.ci", QUANTILES, "relative")
        for j, item in enumerate(items):
            checks.compare(dataset, method, f"alpha without {item}", None, r["drop"][j], "R psych::alpha", FORMULA,
                           "relative")

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
