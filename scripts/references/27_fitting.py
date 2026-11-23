"""Reference values for notebooks/27_fitting.ipynb.

Printed or certified first: the NIST StRD certified values. Where nothing is printed, the help-page example of the
library simplydata calls is run once (R lm, nls, mcr, segmented, drc; SciPy siegelslopes, boxcox, yeojohnson) and
its value saved with the link to that page. Also saves the data the notebook reads (references/data/27_fitting/).
"""

import pandas as pd
from scipy import stats

from _tools import References, NIST_LLS, NIST_NLS, load_csv, r_available, r_dataset, r_run, read_nist, rman

refs = References("27_fitting")
HAS_R = r_available()
CLOSED, ITERATIVE = 1e-10, 1e-6

# 1. NIST nonlinear problems: certified values, data and NIST's second start
PROBLEMS = ["Misra1a", "Chwirut2", "Chwirut1", "Lanczos3", "Gauss1", "Gauss2", "DanWood", "Misra1b", "Kirby2", "Hahn1",
            "MGH17", "Lanczos1", "Lanczos2", "Gauss3", "Misra1c", "Misra1d", "Roszman1", "ENSO", "MGH09", "Thurber",
            "BoxBOD", "Rat42", "MGH10", "Eckerle4", "Rat43", "Bennett5"]
starts = []
for name in PROBLEMS:
    nist = read_nist(name)
    refs.data(f"nist_{name}", nist.data)
    url = NIST_NLS.format(name=name.lower())
    for p, (s1, s2) in nist.starts.items():
        starts.append({"problem": name, "parameter": p, "start 1": s1, "start 2": s2})
    for test_method in ("trust region", "dogbox", "levenberg-marquardt"):
        check = f'fit_formula (start 2, test_method="{test_method}")'
        for p in nist.starts:
            refs.certified(f"NIST {name}", check, f"{p} estimate", nist.certified[p], 6, "NIST StRD certified values", url)
            refs.certified(f"NIST {name}", check, f"{p} standard error", nist.certified[f"{p} standard error"], 6,
                           "NIST StRD certified values", url)
        refs.certified(f"NIST {name}", check, "residual standard deviation", nist.certified["residual standard deviation"], 6,
                       "NIST StRD certified values", url)
refs.data("nist_starts", pd.DataFrame(starts))

# 2. Pontius (certified); cars and Puromycin (run: R lm and nls help pages)
pontius = read_nist("Pontius")
refs.data("nist_Pontius", pontius.data)
for i, term in enumerate(["(intercept)", "x", "x**2"]):
    for suffix, key in (("estimate", f"B{i}"), ("standard error", f"B{i} standard error")):
        refs.certified("NIST Pontius", "fit_curve (of_degree=2)", f"{term} {suffix}", pontius.certified[key], 9,
                       "NIST StRD certified values", NIST_LLS.format(name="Pontius"))

cars = r_dataset("cars")
refs.data("cars", cars)
treated = load_csv("puromycin.csv")
if HAS_R:
    r = r_run('''s <- summary(lm(dist ~ poly(speed, 2, raw = TRUE), weights = speed))$coefficients
                 list(b = s[, 1], se = s[, 2], p = s[, 4])''', cars)
    refs.run_all("cars", 'fit_curve (of_degree=2, weights="speed")',
                 {f"{t} {c}": r[k][i] for i, t in enumerate(["(intercept)", "speed", "speed**2"])
                  for c, k in (("estimate", "b"), ("standard error", "se"), ("p-value", "p"))},
                 "R stats lm (weights =)", rman("stats", "lm"), "stats", CLOSED)
    r = r_run('''m <- nls(rate ~ Vm * conc / (K + conc), start = list(Vm = 200, K = 0.1), weights = conc)
                 s <- summary(m)
                 list(Vm = coef(m)[[1]], K = coef(m)[[2]], `Vm standard error` = s$coefficients[1, 2],
                      `K standard error` = s$coefficients[2, 2], `residual standard error` = s$sigma)''', treated)
    refs.run_all("Puromycin (treated)", 'fit_formula (weights="conc")', r, "R stats nls (weights =), the help page's Puromycin fit",
                 rman("stats", "nls"), "stats", ITERATIVE)

# 3. Passing-Bablok (run: R mcr) and Siegel (run: SciPy)
creat = load_csv("creatinine.csv")
if HAS_R:
    r = r_run('''m <- mcr::mcreg(serum, plasma, method.reg = "PaBa", method.ci = "analytical", alpha = 0.05)
                 p <- m@para
                 list(slope = p[2, 1], `slope 95 % low` = p[2, 3], `slope 95 % high` = p[2, 4], intercept = p[1, 1],
                      `intercept 95 % low` = p[1, 3], `intercept 95 % high` = p[1, 4])''', creat, packages=("mcr",))
    refs.run_all("creatinine (mcr)", 'fit_line (test_method="passing-bablok")', r, "R mcr mcreg (PaBa, analytical), help-page data",
                 rman("mcr", "mcreg"), "mcr", CLOSED)
# The confidence band chart.fit_line(with_confidence_band=True) draws (notebook section 7): R predict.lm(interval =
# "confidence") at the 200 evenly spaced x from the smallest to the largest serum value where the band is drawn; the
# notebook compares the drawn low ends, then the high ends, by their largest relative difference.
if HAS_R:
    import json

    import numpy as np

    grid = np.linspace(creat["serum"].min(), creat["serum"].max(), 200)
    r = r_run(f"""m <- lm(plasma ~ serum)
                  p <- predict(m, data.frame(serum = c({", ".join(repr(float(v)) for v in grid)})), interval = "confidence")
                  list(low = unname(p[, "lwr"]), high = unname(p[, "upr"]))""", creat)
    refs.run("creatinine (mcr)", "chart fit_line (band=True)", "the band drawn (R predict.lm, interval = confidence)", 0.0,
             "R stats predict.lm (interval = confidence)", rman("stats", "predict.lm"), "stats", CLOSED)
    refs.rows[-1].update({"value": float("nan"), "values": json.dumps([float(v) for v in r["low"] + r["high"]]),
                          "tolerance type": "largest relative difference"})
res = stats.siegelslopes(creat["plasma"], creat["serum"], method="hierarchical")
refs.run_all("creatinine (mcr)", 'fit_line (test_method="siegel")', {"slope": res.slope, "intercept": res.intercept},
             "SciPy stats.siegelslopes (hierarchical)",
             "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.siegelslopes.html", "scipy", CLOSED)

# 4. segmented (run: the R help page's call on its stagnant data)
stag = r_dataset("stagnant", "segmented")
refs.data("stagnant", stag)
if HAS_R:
    r = r_run('''m <- lm(y ~ x); s <- segmented::segmented(m, seg.Z = ~ x)
                 ci <- confint(s, digits = 15); sl <- segmented::slope(s, digits = 15)$x; d <- segmented::davies.test(m, seg.Z = ~ x)
                 list(`bend 1` = s$psi[1, 2], `bend 1 95 % CI low` = ci[1, 2], `bend 1 95 % CI high` = ci[1, 3],
                      `slope 1` = sl[1, 1], `slope 1 standard error` = sl[1, 2], `slope 2` = sl[2, 1],
                      `slope 2 standard error` = sl[2, 2], intercept = coef(s)[[1]], `Davies p-value` = d$p.value,
                      `residual standard error` = summary(s)$sigma)''', stag, packages=("segmented",))
    refs.run_all("stagnant (segmented)", "fit_line (with_bends=1)", r, "R segmented segmented + davies.test",
                 rman("segmented", "segmented"), "segmented", ITERATIVE)

# 5. Box-Cox and Yeo-Johnson (run: SciPy, which simplydata calls)
poisons = load_csv("poisons.csv")
t = poisons["time"].to_numpy()
_, lam, ci = stats.boxcox(t, alpha=0.05)
refs.run_all("poisons (Box & Cox 1964)", 'find_best_power_of (test_method="box-cox")', {"power": lam, "95 % low": ci[0], "95 % high": ci[1]},
             "SciPy stats.boxcox", "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.boxcox.html", "scipy", ITERATIVE)
_, lam_yj = stats.yeojohnson(t)
refs.run("poisons (Box & Cox 1964)", 'find_best_power_of (test_method="yeo-johnson")', "power", lam_yj, "SciPy stats.yeojohnson",
         "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.yeojohnson.html", "scipy", ITERATIVE)

# 6. drc (run: the help-page examples of drm and ED on ryegrass; selenium binomial)
rye = load_csv("ryegrass.csv")
DRC = {"3 parameter log-logistic": ("LL.3()", ["b", "d", "e"]), "5 parameter log-logistic": ("LL.5()", ["b", "c", "d", "e", "f"]),
       "weibull": ("W1.4()", ["b", "c", "d", "e"])}
selenium = r_dataset("selenium", "drc")
sel = selenium[selenium["type"] == 1].reset_index(drop=True)
refs.data("selenium_type1", sel)
if HAS_R:
    # drc's default stopping rule (relTol 1e-7), and run to convergence (relTol 1e-12)
    for (kind, (fct, params)), (control, source) in [(k, c) for k in DRC.items() for c in (
            ("", "R drc drm + ED (help-page ryegrass example)"),
            (", control = drc::drmc(relTol = 1e-12, maxIt = 10000)", "R drc drm + ED (relTol 1e-12)"))]:
        r = r_run(f'''m <- drc::drm(root_length ~ conc, fct = drc::{fct}{control})
                      s <- summary(m)$coefficients
                      e <- drc::ED(m, c(10, 50, 90), interval = "delta", display = FALSE)
                      list(b = unname(s[, 1]), se = unname(s[, 2]), ed = unname(e[, 1]), edse = unname(e[, 2]),
                           lo = unname(e[, 3]), hi = unname(e[, 4]))''', rye, packages=("drc",))
        values = {f"{p} estimate": r["b"][i] for i, p in enumerate(params)}
        values.update({f"{p} standard error": r["se"][i] for i, p in enumerate(params)})
        for i, ed in enumerate(("ED10", "ED50", "ED90")):
            values.update({ed: r["ed"][i], f"{ed} standard error": r["edse"][i], f"{ed} 95 % low": r["lo"][i], f"{ed} 95 % high": r["hi"][i]})
        refs.run_all("ryegrass (drc)", f'fit_dose_response (of_kind="{kind}", doses_for_effect=[10, 50, 90])', values,
                     source, rman("drc", "ED.drc"), "drc", ITERATIVE)
    # The Weibull fit: an independent converged optimum (not drc), started from drc's relTol 1e-12 estimates
    # (drc_optimum.py: R optim and SciPy least_squares), and the residual sum of squares at each set of estimates,
    # kept as evidence in references/data.
    import drc_optimum

    WEIBULL = 'fit_dose_response (of_kind="weibull", doses_for_effect=[10, 50, 90])'
    x, y = rye["conc"], rye["root_length"]
    drc_sets = {}
    for label, control in (("drc, default stopping rule (relTol 1e-7)", ""),
                           ("drc, relTol 1e-12", ", control = drc::drmc(relTol = 1e-12, maxIt = 10000)"),
                           ("drc, relTol 1e-15 (simplydata's setting)", ", control = drc::drmc(relTol = 1e-15, maxIt = 10000)")):
        drc_sets[label] = r_run(f"m <- drc::drm(y ~ x, fct = drc::W1.4(){control}); list(est = unname(coef(m)))",
                                {"x": x, "y": y}, packages=("drc",))["est"]
    opt = drc_optimum.optim_r(x, y, "W1.4", drc_sets["drc, relTol 1e-12"])
    lsq = drc_optimum.least_squares_scipy(x, y, "W1.4", drc_sets["drc, relTol 1e-12"])
    OPTIM = "R optim BFGS (reltol 1e-15, analytic gradient)"
    values = {f"{p} estimate": opt["est"][i] for i, p in enumerate("bcde")}
    values.update({f"{p} standard error": opt["se"][i] for i, p in enumerate("bcde")})
    for i, ed in enumerate(("ED10", "ED50", "ED90")):
        values.update({ed: opt["ed"][i], f"{ed} standard error": opt["ed_se"][i], f"{ed} 95 % low": opt["ed_low"][i],
                       f"{ed} 95 % high": opt["ed_high"][i]})
    refs.run_all("ryegrass (drc)", WEIBULL, values, OPTIM, rman("stats", "optim"), "stats", ITERATIVE)
    values = {f"{p} estimate": lsq["est"][i] for i, p in enumerate("bcde")}
    values.update({ed: lsq["ed"][i] for i, ed in enumerate(("ED10", "ED50", "ED90"))})
    refs.run_all("ryegrass (drc)", WEIBULL, values, "SciPy least_squares (1e-15)",
                 "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.optimize.least_squares.html",
                 "scipy", ITERATIVE)
    sets = dict(drc_sets, **{OPTIM: opt["est"], "SciPy least_squares (1e-15)": lsq["est"]})
    found = drc_optimum.rss(x, y, "W1.4", sets)
    refs.data("ryegrass_weibull_rss", pd.DataFrame({"estimates": list(sets), **{p: [v[i] for v in sets.values()]
                                                                               for i, p in enumerate("bcde")},
                                                    "residual sum of squares": [found[n] for n in sets]}))
    r = r_run('''m <- drc::drm(dead / total ~ conc, weights = total, fct = drc::LL.2(), type = "binomial")
                 s <- summary(m)$coefficients; e <- drc::ED(m, 50, interval = "delta", display = FALSE)
                 list(`b estimate` = s[1, 1], `e estimate` = s[2, 1], `b standard error` = s[1, 2],
                      `e standard error` = s[2, 2], `ED50 95 % low` = e[1, 3], `ED50 95 % high` = e[1, 4])''', sel, packages=("drc",))
    refs.run_all("selenium, type 1 (drc)", 'fit_dose_response (outcome="yes or no", out_of="total")', r,
                 "R drc drm (binomial), help-page selenium data", rman("drc", "selenium"), "drc", ITERATIVE)

refs.save()
