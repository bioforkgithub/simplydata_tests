"""Reference values for notebooks/14_small_additions.ipynb, saved to references/14_small_additions.csv.

In this order: the published values (Storey & Tibshirani 2003: pi0 = 0.67; Box & Cox 1964: the power -0.75 and
the suggested -1; Box & Jenkins 1976: the airline model's order), half a unit of the last printed digit; the Python
parent where there is one (statsmodels, SciPy, pandas); the formula written out in NumPy/SciPy where one exists; and
R, run once here (qvalue, forecast, stats, PMCMRplus, drc, car, smatr, mcr, deming).

Tolerances: relative 1e-10 for closed formulas, relative 1e-6 for anything fitted by iteration (ARIMA, the
dose-response curve, the Box-Cox power and its interval).

The tables the notebook reads (the data files of this repository with the columns the notebook needs: the log of
the passengers, the rounding times one row per player and route, one row per patient for the smokers, the log10
leaf measures, the poison and treatment labels as text) are written to references/data/14_small_additions/.

    python scripts/references/14_small_additions.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy import optimize, stats  # noqa: E402
from statsmodels.multivariate.manova import MANOVA  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "14_small_additions"
R_PACKAGES = ("qvalue", "forecast", "PMCMRplus", "drc", "car", "smatr", "mcr", "deming")
checks = Recorder(STEM, R_PACKAGES)
HAS_R = r_available()
CLOSED = 1e-10
ITERATIVE = 1e-6


def ref(dataset, method, quantity, value, source, tolerance, how="relative"):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, how)


def published(dataset, method, quantity, value, source, tolerance, url="", printed=""):
    """A published value (kind "printed"), with the link to the paper."""
    checks.compare(dataset, method, quantity, None, value, source, tolerance, "absolute")
    checks.rows[-1].update({"kind": "printed", "url": url, "digits as printed": printed})


# ---- Storey's q-values and the local FDR: Hedenfalk genes
hedenfalk = load_csv("hedenfalk_p.csv")
checks.data("hedenfalk", hedenfalk)
dataset = "Hedenfalk genes (qvalue)"
p = hedenfalk["p"].to_numpy()
largest = int(np.argmax(p))
published(dataset, 'add_adjusted_p_value_column (method="storey")', "pi0", 0.67, "Storey & Tibshirani (2003), printed",
          0.005, "https://doi.org/10.1073/pnas.1530509100", "0.67")


def storey_q(p, pi0):
    order = np.argsort(p)
    m = len(p)
    q = pi0 * m * p[order] / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1)
    return out


ROWS = [0, 1, 2, int(np.argmin(p)), largest]
if HAS_R:
    r = r_run("q <- qvalue::qvalue(p); list(pi0 = q$pi0, q = q$qvalues, lfdr = q$lfdr)", {"p": p}, packages=("qvalue",))
    method = 'add_adjusted_p_value_column (method="storey")'
    ref(dataset, method, "pi0", r["pi0"], "R qvalue", CLOSED)
    mine = storey_q(p, r["pi0"])
    r_q, r_lfdr = np.asarray(r["q"]), np.asarray(r["lfdr"])
    ref(dataset, method, "q-values below 0.05", float((r_q < 0.05).sum()), "R qvalue", 0, "absolute")
    for i in ROWS:
        ref(dataset, method, f"q-value of row {i + 1}", r_q[i], "R qvalue", CLOSED)
        ref(dataset, method, f"q-value of row {i + 1}", mine[i], "NumPy, the formula (R's pi0)", CLOSED)
    checks.vector(dataset, method, "all 3170 q-values", r_q, "R qvalue", CLOSED)
    method = 'add_adjusted_p_value_column (method="local fdr")'
    ref(dataset, method, "local FDR below 0.05", float((r_lfdr < 0.05).sum()), "R qvalue", 0, "absolute")
    for i in ROWS:
        ref(dataset, method, f"local FDR of row {i + 1}", r_lfdr[i], "R qvalue", CLOSED)
    checks.vector(dataset, method, "all 3170 local FDRs", r_lfdr, "R qvalue", CLOSED)

# ---- the airline model
air = load_csv("air_passengers.csv")
air["log_passengers"] = np.log(air["passengers"])
checks.data("air", air)
dataset = "AirPassengers (log)"
method = 'forecast (order="automatic")'
ORDER = ("p", "d", "q", "seasonal p", "seasonal d", "seasonal q", "period")
for name, expected in zip(ORDER, (0, 1, 1, 0, 1, 1, 12)):
    published(dataset, method, f"order: {name}", expected, "Box & Jenkins (1976), the airline model", 0,
              printed=str(expected))
sar = sm.tsa.SARIMAX(air["log_passengers"].to_numpy(), order=(0, 1, 1), seasonal_order=(0, 1, 1, 12),
                     simple_differencing=True).fit(disp=False, method="lbfgs", maxiter=1000, cov_type="approx",
                                                   optim_kwds={"pgtol": 1e-12, "factr": 10})
parent = {"ma1": sar.params[0], "sma1": sar.params[1], "ma1 standard error": sar.bse[0],
          "sma1 standard error": sar.bse[1], "log-likelihood": sar.llf}
for quantity, value in parent.items():
    ref(dataset, method, quantity, value, "statsmodels SARIMAX (fixed order)", ITERATIVE)
if HAS_R:
    r = r_run('''y <- stats::ts(log(passengers), frequency = 12)
                 a <- forecast::auto.arima(y, stepwise = FALSE, approximation = FALSE)
                 f <- forecast::forecast(a, h = 12, level = 95)
                 b <- stats::arima(y, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12),
                                   method = "ML", optim.control = list(reltol = 1e-14, maxit = 1000))
                 list(order = a$arma[c(1, 6, 2, 3, 7, 4, 5)], auto = unname(a$coef), auto_se = unname(sqrt(diag(a$var.coef))),
                      auto_ll = a$loglik, mean = as.numeric(f$mean), low = as.numeric(f$lower), high = as.numeric(f$upper),
                      fixed = unname(b$coef), fixed_se = unname(sqrt(diag(b$var.coef))), fixed_ll = b$loglik)''',
              {"passengers": air["passengers"]}, packages=("forecast",))
    for name, expected in zip(ORDER, r["order"]):
        ref(dataset, method, f"order: {name}", expected, "R forecast auto.arima", 0, "absolute")
    auto = {"ma1": r["auto"][0], "sma1": r["auto"][1], "ma1 standard error": r["auto_se"][0],
            "sma1 standard error": r["auto_se"][1], "log-likelihood": r["auto_ll"],
            "forecast 1961-01": r["mean"][0], "forecast 1961-12": r["mean"][11], "1961-01 95 % low": r["low"][0],
            "1961-12 95 % high": r["high"][11]}
    for quantity, value in auto.items():
        ref(dataset, method, quantity, value, "R forecast auto.arima, forecast", ITERATIVE)
    fixed = {"ma1": r["fixed"][0], "sma1": r["fixed"][1], "ma1 standard error": r["fixed_se"][0],
             "sma1 standard error": r["fixed_se"][1], "log-likelihood": r["fixed_ll"]}
    for quantity, value in fixed.items():
        ref(dataset, method, quantity, value, "R arima (fixed order, ML)", ITERATIVE)

# ---- after Friedman's test: Nemenyi and Conover
rounding = load_csv("rounding_times.csv")
ROUTES = ["round_out", "narrow_angle", "wide_angle"]
long = rounding.melt(id_vars="player", value_vars=ROUTES, var_name="route", value_name="time")
checks.data("rounding_long", long)
if HAS_R:
    print("Friedman chi-squared of these data (R friedman.test):",
          round(r_run("list(s = unname(friedman.test(as.matrix(cbind(a, b, c)))$statistic))",
                      {"a": rounding["round_out"], "b": rounding["narrow_angle"], "c": rounding["wide_angle"]})["s"], 3),
          "(printed in R's friedman.test help: 11.143)")
dataset = "rounding first base (Hollander & Wolfe)"
mean_rank = rounding[ROUTES].rank(axis=1).mean()
n, k = len(rounding), 3
PAIRS = [(a, b) for i, a in enumerate(sorted(ROUTES)) for b in sorted(ROUTES)[i + 1:]]   # simplydata's order
STATISTIC = {"nemenyi": "q", "conover": "t"}
for test in ("nemenyi", "conover"):
    method = f'compare_each_pair_of (matched_by=, test="{test}")'
    for a, b in PAIRS:
        pair = f"{a} vs {b}"
        for g, name in (("1", a), ("2", b)):
            ref(dataset, method, f"{pair}: mean rank in group {g}", mean_rank[name], "pandas rank (within player)", CLOSED)
        if test == "nemenyi":
            q = abs(mean_rank[a] - mean_rank[b]) / np.sqrt(k * (k + 1) / (12 * n))
            ref(dataset, method, f"{pair}: q", q, "SciPy, the formula", CLOSED)
            ref(dataset, method, f"{pair}: p-value", stats.studentized_range.sf(q, k, np.inf), "SciPy, the formula", CLOSED)
    if HAS_R:
        function = {"nemenyi": "frdAllPairsNemenyiTest(y, g, b)",
                    "conover": 'frdAllPairsConoverTest(y, g, b, p.adjust.method = "holm")'}[test]
        r = r_run(f'''t <- PMCMRplus::{function}
                      list(names = c(rownames(t$p.value), colnames(t$p.value)), p = as.numeric(t$p.value),
                           s = as.numeric(t$statistic))''',
                  {"y": long["time"], "g": long["route"], "b": long["player"].astype(str)}, packages=("PMCMRplus",))
        rows, columns = r["names"][:2], r["names"][2:]
        p_table = pd.DataFrame(np.reshape(r["p"], (2, 2), order="F"), index=rows, columns=columns)
        s_table = pd.DataFrame(np.reshape(r["s"], (2, 2), order="F"), index=rows, columns=columns)
        source = f"R PMCMRplus {function.split('(')[0]}"
        for a, b in PAIRS:
            where = (b, a) if b in rows and a in columns else (a, b)
            ref(dataset, method, f"{a} vs {b}: p-value", p_table.loc[where], source, CLOSED)
            ref(dataset, method, f"{a} vs {b}: {STATISTIC[test]} (size)", abs(s_table.loc[where]), source, CLOSED)

# ---- pairs of percents: Fleiss's smokers
smokers = load_csv("smokers_fleiss.csv")
rows = [{"group": str(g), "smoker": "yes"} for g, s in zip(smokers["group"], smokers["smokers"]) for _ in range(s)]
rows += [{"group": str(g), "smoker": "no"} for g, s, n in zip(smokers["group"], smokers["smokers"], smokers["patients"])
         for _ in range(n - s)]
checks.data("patients", pd.DataFrame(rows))
dataset = "smokers (Fleiss 1981)"
counts = smokers.set_index(smokers["group"].astype(str))
GROUPS = sorted(counts.index)
SMOKER_PAIRS = [(a, b) for i, a in enumerate(GROUPS) for b in GROUPS[i + 1:]]
for test in ("chi-square", "fisher"):
    method = f'compare_each_pair_of (based_on="percents"{", test=" + chr(34) + "fisher" + chr(34) if test == "fisher" else ""})'
    raw = []
    for a, b in SMOKER_PAIRS:
        table = [[counts.loc[g, "smokers"], counts.loc[g, "patients"] - counts.loc[g, "smokers"]] for g in (a, b)]
        raw.append(stats.chi2_contingency(table, correction=True)[1] if test == "chi-square"
                   else stats.fisher_exact(table)[1])
    holm = multipletests(raw, method="holm")[1]
    for i, (a, b) in enumerate(SMOKER_PAIRS):
        ref(dataset, method, f"{a} vs {b}: p-value", raw[i], "SciPy", CLOSED)
        ref(dataset, method, f"{a} vs {b}: p-value (Holm)", holm[i], "SciPy + statsmodels multipletests (Holm)", CLOSED)
    if HAS_R:
        if test == "chi-square":
            r = r_run("list(p = as.numeric(pairwise.prop.test(s, n)$p.value))",
                      {"s": smokers["smokers"], "n": smokers["patients"]})
            name = "R pairwise.prop.test"
        else:
            r = r_run('''raw <- c(); for (i in 1:3) for (j in (i + 1):4)
                           raw <- c(raw, fisher.test(rbind(c(s[i], n[i] - s[i]), c(s[j], n[j] - s[j])))$p.value)
                         list(p = p.adjust(raw, "holm"))''', {"s": smokers["smokers"], "n": smokers["patients"]})
            name = "R fisher.test + p.adjust (Holm)"
        # pairwise.prop.test gives a lower triangle (rows 2-4, columns 1-3), filled column by column
        values = np.asarray(r["p"])
        if test == "chi-square":
            matrix = values.reshape(3, 3, order="F")
            values = [matrix[j - 2, i - 1] for i in range(1, 4) for j in range(i + 1, 5)]
        for i, (a, b) in enumerate(SMOKER_PAIRS):
            ref(dataset, method, f"{a} vs {b}: p-value (Holm)", values[i], name, CLOSED)

# ---- the dose-response curve: ryegrass
ryegrass = load_csv("ryegrass.csv")
checks.data("ryegrass", ryegrass)
dataset, method = "ryegrass (drc)", "fit_dose_response"


def ll4(x, b, c, d, e):
    with np.errstate(divide="ignore"):
        return c + (d - c) / (1 + np.exp(b * (np.log(x) - np.log(e))))


if HAS_R:
    default = r_run('''m <- drc::drm(y ~ x, fct = drc::LL.4())
                       e <- drc::ED(m, c(10, 50, 90), interval = "delta", display = FALSE)
                       list(est = unname(coef(m)), ed = unname(e[, 1]))''',
                    {"x": ryegrass["conc"], "y": ryegrass["root_length"]}, packages=("drc",))
    for k, short in enumerate(("b", "c", "d", "e")):
        ref(dataset, method, f"{short} estimate", default["est"][k], "R drc drm (default stopping rule)", ITERATIVE)
    for k, ed in enumerate(("ED10", "ED50", "ED90")):
        ref(dataset, method, f"{ed} estimate", default["ed"][k], "R drc ED (default stopping rule)", ITERATIVE)
    r = r_run('''m <- drc::drm(y ~ x, fct = drc::LL.4(), control = drc::drmc(relTol = 1e-12))
                 e <- drc::ED(m, c(10, 50, 90), interval = "delta", display = FALSE)
                 list(est = unname(coef(m)), se = unname(sqrt(diag(vcov(m)))), ed = unname(e[, 1]), ed_se = unname(e[, 2]),
                      ed_low = unname(e[, 3]), ed_high = unname(e[, 4]))''',
              {"x": ryegrass["conc"], "y": ryegrass["root_length"]}, packages=("drc",))
    for k, short in enumerate(("b", "c", "d", "e")):
        ref(dataset, method, f"{short} estimate", r["est"][k], "R drc drm (relTol 1e-12)", ITERATIVE)
        ref(dataset, method, f"{short} standard error", r["se"][k], "R drc drm (relTol 1e-12)", ITERATIVE)
    for k, ed in enumerate(("ED10", "ED50", "ED90")):
        for quantity, key in (("estimate", "ed"), ("standard error", "ed_se"), ("95 % CI low", "ed_low"),
                              ("95 % CI high", "ed_high")):
            ref(dataset, method, f"{ed} {quantity}", r[key][k], "R drc ED (delta, relTol 1e-12)", ITERATIVE)
    p_fit, _ = optimize.curve_fit(ll4, ryegrass["conc"], ryegrass["root_length"], p0=r["est"], method="lm",
                                  xtol=1e-14, ftol=1e-14, maxfev=100000)
    for k, short in enumerate(("b", "c", "d", "e")):
        ref(dataset, method, f"{short} estimate", p_fit[k], "SciPy curve_fit (LL.4)", ITERATIVE)
    # An independent converged optimum (not drc), started from drc's estimates (drc_optimum.py): R optim and SciPy
    # least_squares; and the residual sum of squares at each set of estimates, kept as evidence in references/data.
    import drc_optimum

    x, y = ryegrass["conc"], ryegrass["root_length"]
    opt = drc_optimum.optim_r(x, y, "LL.4", r["est"])
    lsq = drc_optimum.least_squares_scipy(x, y, "LL.4", r["est"])
    OPTIM = "R optim BFGS (reltol 1e-15, analytic gradient)"
    for k, short in enumerate(("b", "c", "d", "e")):
        ref(dataset, method, f"{short} estimate", opt["est"][k], OPTIM, ITERATIVE)
        ref(dataset, method, f"{short} standard error", opt["se"][k], OPTIM, ITERATIVE)
        ref(dataset, method, f"{short} estimate", lsq["est"][k], "SciPy least_squares (1e-15)", ITERATIVE)
    for k, ed in enumerate(("ED10", "ED50", "ED90")):
        for quantity, key in (("estimate", "ed"), ("standard error", "ed_se"), ("95 % CI low", "ed_low"),
                              ("95 % CI high", "ed_high")):
            ref(dataset, method, f"{ed} {quantity}", opt[key][k], OPTIM, ITERATIVE)
        ref(dataset, method, f"{ed} estimate", lsq["ed"][k], "SciPy least_squares (1e-15)", ITERATIVE)
    tight15 = r_run('''m <- drc::drm(y ~ x, fct = drc::LL.4(), control = drc::drmc(relTol = 1e-15, maxIt = 10000))
                       list(est = unname(coef(m)))''', {"x": x, "y": y}, packages=("drc",))
    sets = {"drc, default stopping rule (relTol 1e-7)": default["est"], "drc, relTol 1e-12": r["est"],
            "drc, relTol 1e-15 (simplydata's setting)": tight15["est"], "SciPy curve_fit (LL.4)": p_fit,
            OPTIM: opt["est"], "SciPy least_squares (1e-15)": lsq["est"]}
    found = drc_optimum.rss(x, y, "LL.4", sets)
    checks.data("ryegrass_ll4_rss", pd.DataFrame({"estimates": list(sets), "b": [v[0] for v in sets.values()],
                                                  "c": [v[1] for v in sets.values()], "d": [v[2] for v in sets.values()],
                                                  "e": [v[3] for v in sets.values()],
                                                  "residual sum of squares": [found[n] for n in sets]}))

# ---- the Box-Cox power: poisons
poisons = load_csv("poisons.csv").astype({"poison": str, "treatment": str})
checks.data("poisons", poisons)
dataset = "poisons (Box & Cox 1964)"
y = poisons["time"].to_numpy()
X = pd.get_dummies(poisons[["poison", "treatment"]], drop_first=True).astype(float)
X.insert(0, "const", 1.0)


def profile(lam):
    z = np.log(y) if abs(lam) < 1e-12 else (y ** lam - 1) / lam
    rss = sm.OLS(z, X.to_numpy()).fit().ssr
    return -len(y) / 2 * np.log(rss / len(y)) + (lam - 1) * np.log(y).sum()


best = optimize.minimize_scalar(lambda l: -profile(l), bounds=(-3, 3), method="bounded", options={"xatol": 1e-12}).x
second = lambda h: (profile(best + h) - 2 * profile(best) + profile(best - h)) / h ** 2  # noqa: E731
curvature = (4 * second(5e-4) - second(1e-3)) / 3  # Richardson extrapolation of the central difference
se = 1 / np.sqrt(-curvature)
z = stats.norm.ppf(0.975)
mine = {"power": best, "95 % CI low": best - z * se, "95 % CI high": best + z * se}
method = "find_best_power_of (using=)"
BOX_COX = "https://doi.org/10.1111/j.2517-6161.1964.tb00553.x"
published(dataset, method, "power", -0.75, "Box & Cox (1964), printed", 0.005, BOX_COX, "-0.75")
published(dataset, method, "simple power", -1.0, "Box & Cox (1964), suggested", 0, BOX_COX, "-1")
for quantity, value in mine.items():
    ref(dataset, method, quantity, value, "statsmodels OLS + SciPy, the profile (Wald interval)", ITERATIVE)
if HAS_R:
    r = r_run('''s <- summary(car::powerTransform(lm(time ~ poison + treatment)))$result
                 list(power = s[1, 1], low = s[1, 3], high = s[1, 4], simple = s[1, 2])''',
              {"time": poisons["time"], "poison": poisons["poison"], "treatment": poisons["treatment"]}, packages=("car",))
    for quantity, key in (("power", "power"), ("95 % CI low", "low"), ("95 % CI high", "high"), ("simple power", "simple")):
        ref(dataset, method, quantity, r[key], "R car powerTransform", ITERATIVE)
method = "find_best_power_of"
_, power, (low, high) = stats.boxcox(y, alpha=0.05)
for quantity, value in (("power", power), ("95 % CI low", low), ("95 % CI high", high)):
    ref(dataset, method, quantity, value, "SciPy stats.boxcox", ITERATIVE)
if HAS_R:
    r = r_run("list(power = unname(car::powerTransform(time ~ 1)$lambda))", {"time": poisons["time"]}, packages=("car",))
    ref(dataset, method, "power", r["power"], "R car powerTransform", ITERATIVE)

# ---- lines when both columns have errors, or with outliers
leaflife = load_csv("leaflife.csv")
leaflife["log_lma"] = np.log10(leaflife["lma"])
leaflife["log_longev"] = np.log10(leaflife["longev"])
checks.data("leaflife", leaflife)
creatinine = load_csv("creatinine.csv")
checks.data("creatinine", creatinine)
anscombe = load_csv("anscombe.csv")
checks.data("anscombe", anscombe)

dataset, method = "leaflife (smatr)", 'fit_line (method="standardised major axis")'
x, y = leaflife["log_lma"], leaflife["log_longev"]
slope = np.sign(np.corrcoef(x, y)[0, 1]) * y.std() / x.std()
ref(dataset, method, "slope", slope, "the formula", CLOSED)
ref(dataset, method, "intercept", y.mean() - slope * x.mean(), "the formula", CLOSED)
if HAS_R:
    r = r_run("k <- smatr::sma(y ~ x)$coef[[1]]; list(k = as.numeric(as.matrix(k)))", {"x": x, "y": y},
              packages=("smatr",))
    kk = np.reshape(r["k"], (2, 3), order="F")  # rows: intercept, slope; columns: estimate, lower, upper
    for quantity, value in (("intercept", kk[0, 0]), ("intercept 95 % CI low", kk[0, 1]),
                            ("intercept 95 % CI high", kk[0, 2]), ("slope", kk[1, 0]), ("slope 95 % CI low", kk[1, 1]),
                            ("slope 95 % CI high", kk[1, 2])):
        ref(dataset, method, quantity, value, "R smatr sma", CLOSED)

dataset = "creatinine (mcr)"
x, y = creatinine["serum"], creatinine["plasma"]
sxx, syy, sxy = x.var(), y.var(), np.cov(x, y)[0, 1]
for lam in (4, 1):
    method = f'fit_line (method="deming", y_error_over_x_error={lam})'
    slope = (syy - lam * sxx + np.sqrt((syy - lam * sxx) ** 2 + 4 * lam * sxy ** 2)) / (2 * sxy)
    ref(dataset, method, "slope", slope, "the formula", CLOSED)
    ref(dataset, method, "intercept", y.mean() - slope * x.mean(), "the formula", CLOSED)
    if HAS_R:
        r = r_run(f'''m <- mcr::mcreg(x, y, error.ratio = 1 / {lam}, method.reg = "Deming", method.ci = "analytical")
                      list(est = unname(m@para[, "EST"]), se = unname(m@para[, "SE"]), low = unname(m@para[, "LCI"]),
                           high = unname(m@para[, "UCI"]))''', {"x": x, "y": y}, packages=("mcr",))
        for k, term in enumerate(("intercept", "slope")):
            for quantity, key in (("", "est"), (" standard error", "se"), (" 95 % CI low", "low"),
                                  (" 95 % CI high", "high")):
                ref(dataset, method, f"{term}{quantity}", r[key][k], "R mcr mcreg (Deming)", CLOSED)
method = 'fit_line (method="passing-bablok")'
if HAS_R:
    r = r_run('''m <- mcr::mcreg(x, y, method.reg = "PaBa", method.ci = "analytical")
                 list(est = unname(m@para[, "EST"]), low = unname(m@para[, "LCI"]), high = unname(m@para[, "UCI"]))''',
              {"x": x, "y": y}, packages=("mcr",))
    for k, term in enumerate(("intercept", "slope")):
        for quantity, key in (("", "est"), (" 95 % CI low", "low"), (" 95 % CI high", "high")):
            ref(dataset, method, f"{term}{quantity}", r[key][k], "R mcr mcreg (PaBa)", CLOSED)

method = 'fit_line (method="theil-sen")'
for dataset, data, xc, yc in (("Anscombe set 3", anscombe, "x3", "y3"), ("creatinine (mcr)", creatinine, "serum", "plasma")):
    f = stats.theilslopes(data[yc], data[xc], alpha=0.95, method="separate")
    for quantity, value in (("slope", f.slope), ("intercept", f.intercept), ("slope 95 % CI low", f.low_slope),
                            ("slope 95 % CI high", f.high_slope)):
        ref(dataset, method, quantity, value, "SciPy theilslopes", CLOSED)
    if HAS_R:
        r = r_run("list(slope = unname(coef(deming::theilsen(y ~ x))[2]))", {"x": data[xc], "y": data[yc]},
                  packages=("deming",))
        ref(dataset, method, "slope", r["slope"], "R deming theilsen", CLOSED)

# ---- the covariance table and MANOVA
mtcars = load_csv("mtcars.csv")
checks.data("mtcars", mtcars)
iris = load_csv("iris.csv")
checks.data("iris", iris)
dataset, method = "mtcars", 'find_correlation_between (show="covariance")'
parent = mtcars[["mpg", "wt", "hp"]].cov()
r = r_run("list(c = as.numeric(cov(cbind(mpg, wt, hp))))", mtcars[["mpg", "wt", "hp"]]) if HAS_R else None
for i, a in enumerate(("mpg", "wt", "hp")):
    for j, b in enumerate(("mpg", "wt", "hp")):
        if j < i:
            continue
        ref(dataset, method, f"{a} with {b}", parent.loc[a, b], "pandas DataFrame.cov", CLOSED)
        if r is not None:
            ref(dataset, method, f"{a} with {b}", r["c"][i + 3 * j], "R cov", CLOSED)

TESTS = {"Pillai's trace": "Pillai", "Wilks' lambda": "Wilks", "Hotelling-Lawley trace": "Hotelling-Lawley",
         "Roy's greatest root": "Roy"}
QUANTITIES = ("statistic", "F", "degrees of freedom 1", "degrees of freedom 2", "p-value")
dataset = "iris"
for columns in (["Sepal.Length", "Petal.Length"], ["Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"]):
    method = f"compare_averages_of (MANOVA of {len(columns)} columns)"
    names = {c: f"v{i}" for i, c in enumerate(columns)}
    fit = MANOVA.from_formula(" + ".join(names.values()) + " ~ Species", data=iris.rename(columns=names))
    s = fit.mv_test().results["Species"]["stat"]
    for test in TESTS:
        for quantity, key in zip(QUANTITIES, ("Value", "F Value", "Num DF", "Den DF", "Pr > F")):
            ref(dataset, method, f"{test}: {quantity}", s.loc[test, key], "statsmodels MANOVA", CLOSED)
    if HAS_R:
        for test, r_name in TESTS.items():
            r = r_run(f'''m <- manova(cbind({", ".join(names.values())}) ~ Species)
                          s <- summary(m, test = "{r_name}")$stats
                          list(v = unname(s[1, 2:6]))''', iris.rename(columns=names))
            for quantity, value in zip(QUANTITIES, r["v"]):
                ref(dataset, method, f"{test}: {quantity}", value, "R summary.manova", CLOSED)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
