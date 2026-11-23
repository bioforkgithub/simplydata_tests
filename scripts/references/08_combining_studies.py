"""Reference values for notebooks/08_combining_studies.ipynb, saved to references/08_combining_studies.csv.

simplydata's combine_studies calls R metafor, so metafor is the original. In order of preference:
- printed: the BCG worked example of Viechtbauer (2010), JSS 36(3), and the pooled odds ratio of passive smoking
  published by Hackshaw et al. (1997), BMJ 315, to half a unit of the last printed digit;
- run once here, only for what nothing prints: R metafor (escalc, rma, rma.mh, rma.peto, regtest, trimfill) and, for
  Peters' test, which metafor does not have, R meta (metabias), with the settings simplydata states it uses (REML,
  Hartung-Knapp with metafor's ad hoc rule test = "adhoc", the prediction interval from predict()). Relative
  difference 1e-10 for results from closed formulas (each study's estimate, a common effect, DerSimonian-Laird,
  Mantel-Haenszel, Peto, Q, the regression tests) and 1e-6 for results found by iteration (REML, Paule-Mandel, trim
  and fill). The BCG models fitted by iteration are compared with metafor run to convergence, control = list(tol =
  1e-12, threshold = 1e-12): metafor's defaults stop at a tolerance of 1e-5 (FINDINGS A7 is about simplydata
  stopping early there); the other iterative fits with metafor's defaults, as before.
The values on a log scale (risk, odds and rate ratios) or on Fisher's z scale (correlations) are compared on that
scale; the notebook puts simplydata's answers on it with table.add_column.
The data the notebook reads are written to references/data/08_combining_studies/.

    python scripts/references/08_combining_studies.py      (from the repository root; needs R with metafor and meta)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402

from checks import load_csv, r_available, r_eval_text, r_run  # noqa: E402

NOTEBOOK = "08_combining_studies"
CLOSED = 1e-10
ITERATIVE = 1e-6
VIECHTBAUER = "Viechtbauer (2010), printed"
V_TITLE = "Viechtbauer (2010), Conducting meta-analyses in R with the metafor package, J. Stat. Softw. 36(3)"
V_URL = "https://doi.org/10.18637/jss.v036.i03"
HACKSHAW = "Hackshaw et al. (1997), printed"
H_TITLE = ("Hackshaw, Law & Wald (1997), The accumulated evidence on lung cancer and environmental tobacco smoke, "
           "BMJ 315, 980-988")
H_URL = "https://doi.org/10.1136/bmj.315.7114.980"
TIGHT = ", control = list(tol = 1e-12, threshold = 1e-12, maxiter = 1000))"
TIGHT_SOURCE = "R metafor rma (tol 1e-12)"   # metafor run to convergence (its default stops at tol and threshold 1e-5)

ref = Recorder(NOTEBOOK, r_packages=("metafor", "meta"))
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")
VERSIONS = {p: r_eval_text(f'as.character(packageVersion("{p}"))') for p in ("metafor", "meta")}


def rman(package, topic):
    return f"https://search.r-project.org/CRAN/refmans/{package}/html/{topic}.html (version {VERSIONS[package]})"


def add(dataset, method, quantity, value, source, tolerance, url):
    ref.compare(dataset, method, quantity, None, value, source, tolerance, "relative" if tolerance else "absolute")
    ref.rows[-1]["url"] = url


def add_all(dataset, method, values, source, tolerance, url, leave_out=()):
    for quantity, value in values.items():
        if quantity not in leave_out:
            add(dataset, method, quantity, value, source, tolerance, url)


def r_rma(fit_call, data, escalc_call):
    """R metafor: escalc, the model, predict() (on the scale it is combined on)."""
    code = f"""e <- {escalc_call}
f <- {fit_call}
p <- predict(f)
out <- list(estimate = f$b[1], `95 % CI low` = f$ci.lb[1], `95 % CI high` = f$ci.ub[1], `standard error` = f$se[1],
            `p-value` = f$pval[1], Q = f$QE)
if (f$method != "EE") {{
  out$`tau squared` <- f$tau2; out$`I squared %` <- f$I2
  out$`prediction interval low` <- p$pi.lb; out$`prediction interval high` <- p$pi.ub }}
out"""
    return r_run(code, data, packages=("metafor",))


def printed(dataset, method, quantity, shown, source, title, url):
    ref.printed(dataset, method, quantity, shown, source, title, url)


# ---------------------------------------------------------------- BCG vaccine trials
bcg = load_csv("dat_bcg.csv")
bcg["name"] = bcg["author"] + " " + bcg["year"].astype(str)
bcg["vaccinated"] = bcg["tpos"] + bcg["tneg"]
bcg["unvaccinated"] = bcg["cpos"] + bcg["cneg"]
ref.data("bcg", bcg)

# the published values (REML, normal interval)
for quantity, shown in (("log risk ratio", "-0.7145"), ("standard error", "0.1798"), ("95 % CI low (log)", "-1.0669"),
                        ("95 % CI high (log)", "-0.3622"), ("tau squared", "0.3132"), ("I squared %", "92.22"),
                        ("Q", "152.2330"), ("risk ratio", "0.4894"), ("risk ratio 95 % CI low", "0.3441"),
                        ("risk ratio 95 % CI high", "0.6962"), ("prediction interval low", "0.1546"),
                        ("prediction interval high", "1.5490")):
    printed("BCG trials", "combine_studies (REML, normal)", quantity, shown, VIECHTBAUER, V_TITLE, V_URL)
printed("BCG trials", "combine_studies (DerSimonian-Laird, normal)", "tau squared", "0.3088", VIECHTBAUER, V_TITLE,
        V_URL)
printed("BCG trials", "combine_studies (using ablat)", "intercept", "0.2515", VIECHTBAUER, V_TITLE, V_URL)
printed("BCG trials", "combine_studies (using ablat)", "slope", "-0.0291", VIECHTBAUER, V_TITLE, V_URL)

# every study's log risk ratio and standard error
e = r_run("e <- metafor::escalc(measure = 'RR', ai = tpos, n1i = vaccinated, ci = cpos, n2i = unvaccinated)\n"
          "list(yi = as.numeric(e$yi), vi = as.numeric(e$vi))", bcg, packages=("metafor",))
for i, name in enumerate(bcg["name"]):
    add("BCG trials", "combine_studies (each study)", f"{name}: log risk ratio", e["yi"][i], "R metafor escalc", CLOSED,
        rman("metafor", "escalc"))
    add("BCG trials", "combine_studies (each study)", f"{name}: standard error", np.sqrt(e["vi"][i]),
        "R metafor escalc", CLOSED, rman("metafor", "escalc"))

# every method= and interval_by=
ESCALC_RR = "metafor::escalc(measure = 'RR', ai = tpos, n1i = vaccinated, ci = cpos, n2i = unvaccinated)"
PRINTED_REML = ("estimate", "standard error", "95 % CI low", "95 % CI high", "tau squared", "I squared %", "Q",
                "prediction interval low", "prediction interval high")
RUNS = [  # label, R model, iterative, quantities a printed value already covers
    ("REML, Hartung-Knapp", 'metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc")', True, ()),
    ("REML, normal", 'metafor::rma(yi, vi, data = e, method = "REML", test = "z")', True, PRINTED_REML),
    ("DerSimonian-Laird, normal", 'metafor::rma(yi, vi, data = e, method = "DL", test = "z")', False,
     ("tau squared",)),
    ("DerSimonian-Laird, Hartung-Knapp", 'metafor::rma(yi, vi, data = e, method = "DL", test = "adhoc")', False, ()),
    ("Paule-Mandel, normal", 'metafor::rma(yi, vi, data = e, method = "PM", test = "z")', True, ()),
    ("Paule-Mandel, Hartung-Knapp", 'metafor::rma(yi, vi, data = e, method = "PM", test = "adhoc")', True, ()),
    ("common effect", 'metafor::rma(yi, vi, data = e, method = "EE")', False, ()),
]
for label, model, iterative, printed_already in RUNS:
    if iterative:
        values = r_rma(model[:-1] + TIGHT, bcg, ESCALC_RR)
        source, tolerance = "R metafor rma (tol 1e-12)", ITERATIVE
    else:
        values = r_rma(model, bcg, ESCALC_RR)
        source, tolerance = "R metafor rma", CLOSED
    add_all("BCG trials", f"combine_studies ({label})", values, source, tolerance, rman("metafor", "rma.uni"),
            printed_already)

# Mantel-Haenszel and Peto; odds ratios and risk differences
for show, measure in (("risk ratios", "RR"), ("odds ratios", "OR")):
    r = r_run(f"f <- metafor::rma.mh(ai = tpos, bi = tneg, ci = cpos, di = cneg, measure = '{measure}')\n"
              "list(estimate = f$b[1], `95 % CI low` = f$ci.lb, `95 % CI high` = f$ci.ub, `standard error` = f$se, "
              "`p-value` = f$pval, Q = f$QE)", bcg, packages=("metafor",))
    add_all("BCG trials", f"combine_studies (Mantel-Haenszel, {show})", r, "R metafor rma.mh", CLOSED,
            rman("metafor", "rma.mh"))
r = r_run("f <- metafor::rma.peto(ai = tpos, bi = tneg, ci = cpos, di = cneg)\n"
          "list(estimate = f$b[1], `95 % CI low` = f$ci.lb, `95 % CI high` = f$ci.ub, `standard error` = f$se, "
          "`p-value` = f$pval, Q = f$QE)", bcg, packages=("metafor",))
add_all("BCG trials", "combine_studies (Peto)", r, "R metafor rma.peto", CLOSED, rman("metafor", "rma.peto"))
for show, measure in (("odds ratios", "OR"), ("risk differences", "RD")):
    r = r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc")', bcg,
              f"metafor::escalc(measure = '{measure}', ai = tpos, bi = tneg, ci = cpos, di = cneg)")
    add_all("BCG trials", f"combine_studies (REML, {show})", r, "R metafor rma", ITERATIVE, rman("metafor", "rma.uni"))
    r = r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc"' + TIGHT, bcg,
              f"metafor::escalc(measure = '{measure}', ai = tpos, bi = tneg, ci = cpos, di = cneg)")
    add_all("BCG trials", f"combine_studies (REML, {show})", r, TIGHT_SOURCE, ITERATIVE, rman("metafor", "rma.uni"))

# subgroups (metafor's default stopping rule, and run to convergence with TIGHT)
METHOD = "combine_studies (for_each_group_in alloc)"
for control, source in (("", "R metafor rma"), (TIGHT[:-1], TIGHT_SOURCE)):
    r = r_run("e <- " + ESCALC_RR + f"""
g <- sort(unique(alloc))
z <- lapply(g, function(x) metafor::rma(yi, vi, data = e, subset = alloc == x, test = "z"{control}))
k <- lapply(g, function(x) metafor::rma(yi, vi, data = e, subset = alloc == x, test = "adhoc"{control}))
q <- metafor::rma(sapply(z, function(m) m$b[1]), sei = sapply(z, function(m) m$se), mods = ~ factor(g), method = "EE")
list(`Q between` = q$QM, `p between` = q$QMp, est = sapply(k, function(m) m$b[1]), tau2 = sapply(k, function(m) m$tau2),
     low = sapply(k, function(m) m$ci.lb), groups = g)""", bcg, packages=("metafor",))
    add("BCG trials", METHOD, "Q between", r["Q between"], source, ITERATIVE, rman("metafor", "rma.uni"))
    add("BCG trials", METHOD, "p between", r["p between"], source, ITERATIVE, rman("metafor", "rma.uni"))
    for i, group in enumerate(r["groups"]):
        for quantity, key in (("log risk ratio", "est"), ("tau squared", "tau2"), ("95 % CI low (log)", "low")):
            add("BCG trials", METHOD, f"{group}: {quantity}", r[key][i], source, ITERATIVE,
                rman("metafor", "rma.uni"))

# meta-regression on latitude (Hartung-Knapp; the default stopping rule, and TIGHT)
for control, source in (("", "R metafor rma"), (TIGHT[:-1], TIGHT_SOURCE)):
    r = r_run("e <- " + ESCALC_RR + f"""
f <- metafor::rma(yi, vi, mods = ~ ablat, data = e, method = "REML", test = "adhoc"{control})
list(b = as.numeric(f$b), low = f$ci.lb, high = f$ci.ub, p = f$pval)""", bcg, packages=("metafor",))
    for j, term in enumerate(("intercept", 'per 1 more "ablat"')):
        for quantity, key in (("estimate", "b"), ("95 % CI low", "low"), ("95 % CI high", "high"), ("p-value", "p")):
            add("BCG trials", "combine_studies (using ablat)", f"{term}: {quantity}", r[key][j], source, ITERATIVE,
                rman("metafor", "rma.uni"))

# Peters' test (R meta) and trim and fill (R metafor)
r = r_run("""m <- meta::metabin(tpos, tpos + tneg, cpos, cpos + cneg, sm = "RR")
b <- meta::metabias(m, method.bias = "Peters", k.min = 10)
list(t = unname(b$statistic), p = b$pval)""", bcg, packages=("meta",))
add("BCG trials", "combine_studies (Peters' test)", "t", r["t"], "R meta metabias", CLOSED, rman("meta", "metabias"))
add("BCG trials", "combine_studies (Peters' test)", "p-value", r["p"], "R meta metabias", CLOSED,
    rman("meta", "metabias"))
r = r_run("e <- " + ESCALC_RR + """
f <- metafor::trimfill(metafor::rma(yi, vi, data = e, test = "adhoc"))
list(k0 = f$k0, est = f$b[1], low = f$ci.lb)""", bcg, packages=("metafor",))
add("BCG trials", "combine_studies (trim and fill)", "studies added", r["k0"], "R metafor trimfill", 0,
    rman("metafor", "trimfill"))
add("BCG trials", "combine_studies (trim and fill)", "log risk ratio", r["est"], "R metafor trimfill", ITERATIVE,
    rman("metafor", "trimfill"))
add("BCG trials", "combine_studies (trim and fill)", "95 % CI low (log)", r["low"], "R metafor trimfill", ITERATIVE,
    rman("metafor", "trimfill"))

# ---------------------------------------------------------------- passive smoking
smoke = load_csv("dat_hackshaw1998.csv")
smoke["name"] = smoke["author"] + " " + smoke["year"].astype(str) + " (" + smoke["study"].astype(str) + ")"
smoke["odds ratio"] = smoke["or"]
smoke["odds ratio 95 % CI low"] = smoke["or.lb"]
smoke["odds ratio 95 % CI high"] = smoke["or.ub"]
ref.data("smoke", smoke)
for quantity, shown in (("odds ratio", "1.24"), ("95 % CI low", "1.13"), ("95 % CI high", "1.36")):
    printed("passive smoking", "combine_studies (DerSimonian-Laird, normal)", quantity, shown, HACKSHAW, H_TITLE, H_URL)
CONV = ("metafor::conv.wald(out = or, ci.lb = or.lb, ci.ub = or.ub, transf = log, "
        "data = data.frame(or = or, or.lb = or.lb, or.ub = or.ub))")
add_all("passive smoking", "combine_studies (DerSimonian-Laird, normal)",
        r_rma('metafor::rma(yi, vi, data = e, method = "DL", test = "z")', smoke, CONV), "R metafor rma", CLOSED,
        rman("metafor", "rma.uni"), leave_out=("estimate", "95 % CI low", "95 % CI high"))
add_all("passive smoking", "combine_studies (REML, Hartung-Knapp)",
        r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc")', smoke, CONV), "R metafor rma",
        ITERATIVE, rman("metafor", "rma.uni"))
add_all("passive smoking", "combine_studies (REML, Hartung-Knapp)",
        r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc"' + TIGHT, smoke, CONV), TIGHT_SOURCE,
        ITERATIVE, rman("metafor", "rma.uni"))
r = r_run(f"""e <- {CONV}
g <- metafor::regtest(e$yi, e$vi, model = "lm")
f <- metafor::trimfill(metafor::rma(yi, vi, data = e, test = "adhoc"))
list(t = g$zval, p = g$pval, k0 = f$k0, est = f$b[1])""", smoke, packages=("metafor",))
add("passive smoking", "combine_studies (Egger's test)", "t", r["t"], "R metafor regtest", CLOSED,
    rman("metafor", "regtest"))
add("passive smoking", "combine_studies (Egger's test)", "p-value", r["p"], "R metafor regtest", CLOSED,
    rman("metafor", "regtest"))
add("passive smoking", "combine_studies (trim and fill)", "studies added", r["k0"], "R metafor trimfill", 0,
    rman("metafor", "trimfill"))
add("passive smoking", "combine_studies (trim and fill)", "log odds ratio", r["est"], "R metafor trimfill", ITERATIVE,
    rman("metafor", "trimfill"))

# ---------------------------------------------------------------- averages, events in time, correlations
normand = load_csv("dat_normand1999.csv")
curtis = load_csv("dat_curtis1998.csv")
curtis["name"] = "experiment " + curtis["id"].astype(str)
hart = load_csv("dat_hart1999.csv")
mcdaniel = load_csv("dat_mcdaniel1994.csv")
mcdaniel["name"] = "study " + mcdaniel["study"].astype(str)
for name, frame in (("normand", normand), ("curtis", curtis), ("hart", hart), ("mcdaniel", mcdaniel)):
    ref.data(name, frame)
MEANS = "m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i"
SETS = [  # dataset label, data, method label, escalc call, study names
    ("Normand 1999", normand, "default", f"metafor::escalc(measure = 'SMD', {MEANS})", normand["source"]),
    ("Normand 1999", normand, "differences in averages", f"metafor::escalc(measure = 'MD', {MEANS})", normand["source"]),
    ("Curtis 1998", curtis, "ratios of averages", f"metafor::escalc(measure = 'ROM', {MEANS})", curtis["name"]),
    ("Hart 1999", hart, "default", "metafor::escalc(measure = 'IRR', x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i)",
     hart["study"]),
    ("McDaniel 1994", mcdaniel, "correlations", "metafor::escalc(measure = 'ZCOR', ri = ri, ni = ni)", mcdaniel["name"]),
]
for label, data, show, escalc_call, names in SETS:
    method = f"combine_studies ({show})"
    values = r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc")', data, escalc_call)
    add_all(label, method, values, "R metafor rma", ITERATIVE, rman("metafor", "rma.uni"))
    values = r_rma('metafor::rma(yi, vi, data = e, method = "REML", test = "adhoc"' + TIGHT, data, escalc_call)
    add_all(label, method, values, TIGHT_SOURCE, ITERATIVE, rman("metafor", "rma.uni"))
    e = r_run(f"e <- {escalc_call}\nlist(yi = as.numeric(e$yi), vi = as.numeric(e$vi))", data, packages=("metafor",))
    for i, name in enumerate(names):
        add(label, f"{method} (each study)", f"{name}: estimate", e["yi"][i], "R metafor escalc", CLOSED,
            rman("metafor", "escalc"))
        add(label, f"{method} (each study)", f"{name}: standard error", np.sqrt(e["vi"][i]), "R metafor escalc", CLOSED,
            rman("metafor", "escalc"))
ref.save()
