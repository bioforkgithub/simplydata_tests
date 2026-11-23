"""Reference values for notebooks/16_cause_and_effect.ipynb, saved to references/16_cause_and_effect.csv.

- Printed: the attributable fractions of the first example of R epiR's epi.2by2 help page (75.05 %, 32.91 to 90.72;
  54.20 %), the directly standardised rate of the R epitools ageadjust.direct help page (92.3, 80.4 to 105.8 per
  100,000) and the ratios and 90 % intervals of Schutz, Labes & Fuglsang (2014, Table II): half a unit of the last
  printed digit.
- Run here: R WeightIt, cobalt, MatchIt (+ lm + sandwich vcovCL), mediation, epiR, epitools, poisson.test and lm;
  statsmodels (GLM, OLS, Mediation with the same seed); the textbook formulas in NumPy and SciPy.
- Monte Carlo (mediation): the allowed difference is 4 Monte Carlo standard errors of the simulated average,
  worked out from statsmodels' own simulation (times sqrt(2) against R's independent simulation).

The tables the notebook reads that are not data files of this repository as they are (the cats one row per cat;
each bioequivalence dataset alone) are written to references/data/16_cause_and_effect/.

    python scripts/references/16_cause_and_effect.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.mediation import Mediation  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "16_cause_and_effect"
checks = Recorder(STEM, ("WeightIt", "MatchIt", "cobalt", "mediation", "epiR", "epitools"))
HAS_R = r_available()
CLOSED = 1e-10
LEAST_SQUARES = 1e-9
ITERATIVE = 1e-6


def ref(dataset, method, quantity, value, source, tolerance, how="relative"):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, how)


def refs(dataset, method, values, source, tolerance, how="relative"):
    for quantity, value in values.items():
        ref(dataset, method, quantity, value, source, tolerance, how)


# ---- lalonde
lalonde = load_csv("lalonde.csv")
COLUMNS = ["age", "educ", "race", "married", "nodegree", "re74", "re75"]
LALONDE_R = ("d <- data.frame(t = treat, y = re78, age = age, educ = educ, race = race, married = married, "
             "nodegree = nodegree, re74 = re74, re75 = re75)\n"
             "f <- t ~ age + educ + race + married + nodegree + re74 + re75\n")
DS = "lalonde"
M = "compare_balanced_groups_of (weighting, everyone)"
X = pd.get_dummies(lalonde[COLUMNS], columns=["race"], drop_first=True, dtype=float)
score = sm.GLM(lalonde["treat"], sm.add_constant(X), family=sm.families.Binomial()).fit(tol=1e-12).predict()
w = np.where(lalonde["treat"] == 1, 1 / score, 1 / (1 - score))
treated = lalonde["treat"] == 1
effect = (np.average(lalonde["re78"][treated], weights=w[treated])
          - np.average(lalonde["re78"][~treated], weights=w[~treated]))
ref(DS, M, "effect", effect, "statsmodels GLM + weighted means", ITERATIVE)
if HAS_R:
    r = r_run(LALONDE_R + """W <- WeightIt::weightit(f, data = d, method = "glm", estimand = "ATE")
        co <- summary(WeightIt::lm_weightit(y ~ t, data = d, weightit = W))$coefficients
        list(effect = co["t", 1], `standard error` = co["t", 2], `p-value` = co["t", 4])""", lalonde)
    refs(DS, M, r, "R WeightIt lm_weightit", ITERATIVE)

    r = r_run(LALONDE_R + """W <- WeightIt::weightit(f, data = d, method = "glm", estimand = "ATE")
        b <- cobalt::bal.tab(W, un = TRUE, binary = "std", continuous = "std")$Balance
        b <- b[rownames(b) != "prop.score", ]
        list(name = rownames(b), before = b$Diff.Un, after = b$Diff.Adj)""", lalonde)
    names = {"race_black": "race: black", "race_hispan": "race: hispan", "race_white": "race: white"}
    for name, before, after in zip(r["name"], r["before"], r["after"]):
        name = names.get(name, name)
        ref(DS, "compare_balanced_groups_of (show=balance)", f"{name} difference before", before, "R cobalt bal.tab",
            ITERATIVE)
        ref(DS, "compare_balanced_groups_of (show=balance)", f"{name} difference after", after, "R cobalt bal.tab",
            ITERATIVE)

    for caliper in (None, 0.2):
        M = "compare_balanced_groups_of (matching, trained)" + ("" if caliper is None else ", match_within=0.2")
        how = "" if caliper is None else ', link = "linear.logit", caliper = 0.2'
        r = r_run(LALONDE_R + f"""m <- MatchIt::matchit(f, data = d, method = "nearest", distance = "glm"{how})
            md <- MatchIt::match.data(m)
            fit <- stats::lm(y ~ t, data = md, weights = weights)
            co <- lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = ~subclass))
            list(effect = co["t", 1], `standard error` = co["t", 2])""", lalonde)
        refs(DS, M, r, "R MatchIt + lm + vcovCL", ITERATIVE)

# ---- JOBS II mediation
jobs = load_csv("jobs_ii.csv")
DS, M = "JOBS II", "explain (through=)"
mediator = smf.ols("job_seek ~ treat + econ_hard + sex + age", data=jobs)
outcome = smf.ols("depress2 ~ treat + job_seek + econ_hard + sex + age", data=jobs)
fit = Mediation(outcome, mediator, "treat", "job_seek").fit(method="parametric", n_rep=1000, rng=1)
summary = fit.summary(alpha=0.05)
labels = {"indirect effect": "ACME (average)", "direct effect": "ADE (average)", "total effect": "Total effect"}
parent = {f"{k} {c}": summary.loc[v, col] for k, v in labels.items()
          for c, col in (("estimate", "Estimate"), ("95 % CI low", "Lower CI bound"), ("95 % CI high", "Upper CI bound"))}
refs(DS, M, parent, "statsmodels Mediation (same seed)", CLOSED)
draws = {"indirect effect": (fit.ACME_ctrl + fit.ACME_tx) / 2, "direct effect": (fit.ADE_ctrl + fit.ADE_tx) / 2,
         "total effect": fit.total_effect}
mcse = {k: float(np.std(v, ddof=1) / np.sqrt(len(v))) for k, v in draws.items()}
print("Monte Carlo standard errors:", {k: round(v, 5) for k, v in mcse.items()})
a = mediator.fit().params["treat"]
fit_out = outcome.fit()
limit = {"indirect effect": a * fit_out.params["job_seek"], "direct effect": fit_out.params["treat"]}
limit["total effect"] = limit["indirect effect"] + limit["direct effect"]
for k, value in limit.items():
    ref(DS, M, f"{k} estimate", value, "exact limit (product of coefficients, statsmodels OLS)", 4 * mcse[k], "absolute")
if HAS_R:
    r = r_run("""set.seed(1)
        med <- lm(job_seek ~ treat + econ_hard + sex + age)
        out <- lm(depress2 ~ treat + job_seek + econ_hard + sex + age)
        m <- mediation::mediate(med, out, treat = "treat", mediator = "job_seek", sims = 1000)
        list(`indirect effect estimate` = m$d.avg, `direct effect estimate` = m$z.avg, `total effect estimate` = m$tau.coef,
             a = coef(med)[["treat"]] * coef(out)[["job_seek"]])""", jobs)
    for k in ("indirect effect", "direct effect", "total effect"):
        ref(DS, M, f"{k} estimate", r[f"{k} estimate"], "R mediation::mediate", 4 * np.sqrt(2) * mcse[k], "absolute")
    # the notebook multiplies simplydata's own two linear-model coefficients (explain), as the limit does
    ref(DS, "R lm (check of the limit)", "a * b", r["a"], "R lm", LEAST_SQUARES)

# ---- attributable fractions: dry cat food and FLUTD (Willeberg 1977)
flutd = load_csv("flutd.csv")
rows = [{"food": f, "flutd": "yes"} for f, e in zip(flutd["food"], flutd["flutd"]) for _ in range(e)]
rows += [{"food": f, "flutd": "no"} for f, e, n in zip(flutd["food"], flutd["flutd"], flutd["cats"]) for _ in range(n - e)]
checks.data("cats", pd.DataFrame(rows))
DS, M = "dry cat food and FLUTD", "compare_risks (attributable fractions)"
EPIR_URL = "https://search.r-project.org/CRAN/refmans/epiR/html/epi.2by2.html"
for quantity, shown in {"AF among the exposed": "75.05", "AF 95 % CI low": "32.91", "AF 95 % CI high": "90.72",
                        "PAF": "54.20"}.items():
    checks.printed(DS, M, quantity + " (percent)", shown, "epiR help, printed", title="R epiR epi.2by2, help page example 1",
                   url=EPIR_URL)
a1, n1 = flutd["flutd"][0], flutd["cats"][0]
a0, n0 = flutd["flutd"][1], flutd["cats"][1]
rr = (a1 / n1) / (a0 / n0)
se = np.sqrt(1 / a1 - 1 / n1 + 1 / a0 - 1 / n0)
low, high = rr * np.exp(-stats.norm.ppf(0.975) * se), rr * np.exp(stats.norm.ppf(0.975) * se)
risk_all = (a1 + a0) / (n1 + n0)
paf = (risk_all - a0 / n0) / risk_all
b1, b0 = n1 - a1, n0 - a0
paf_var = (b1 + paf * (a1 + b0)) / ((n1 + n0) * a0)  # Jewell (2004): a = a1, b = b1, c = a0, d = b0
formula = {"AF among the exposed": (rr - 1) / rr, "AF 95 % CI low": (low - 1) / low, "AF 95 % CI high": (high - 1) / high,
           "PAF": paf, "PAF 95 % CI low": 1 - np.exp(np.log(1 - paf) + stats.norm.ppf(0.975) * np.sqrt(paf_var)),
           "PAF 95 % CI high": 1 - np.exp(np.log(1 - paf) - stats.norm.ppf(0.975) * np.sqrt(paf_var))}
refs(DS, M, formula, "NumPy, textbook formula", CLOSED)
if HAS_R:
    r = r_run("""r <- epiR::epi.2by2(as.table(matrix(c(flutd[1], cats[1] - flutd[1], flutd[2], cats[2] - flutd[2]),
                                                 nrow = 2, byrow = TRUE)), method = "cohort.count")$massoc.detail
        list(`AF among the exposed` = r$AFRisk.strata.wald$est, `AF 95 % CI low` = r$AFRisk.strata.wald$lower,
             `AF 95 % CI high` = r$AFRisk.strata.wald$upper, PAF = r$PAFRisk.strata.wald$est,
             `PAF 95 % CI low` = r$PAFRisk.strata.wald$lower, `PAF 95 % CI high` = r$PAFRisk.strata.wald$upper)""", flutd)
    refs(DS, M, r, "R epiR epi.2by2", ITERATIVE)

# ---- standardised rates: Down syndrome (Selvin 2004)
down = load_csv("down_syndrome.csv")
down["birth_order"] = down["birth_order"].astype(str)
DS, M = "Down syndrome by birth order", "rate_of (standard_population=)"
standard = down.groupby("mother_age")["births"].sum()
ORDERS = list(dict.fromkeys(down["birth_order"]))
for order in ORDERS:
    part = down[down["birth_order"] == order].set_index("mother_age")
    w = standard[part.index] / standard[part.index].sum()
    rate = part["cases"] / part["births"]
    dsr = (w * rate).sum()
    var = (w ** 2 * part["cases"] / part["births"] ** 2).sum()
    wm = (w / part["births"]).max()
    formula = {"crude rate": part["cases"].sum() / part["births"].sum() * 1e5, "adjusted rate": dsr * 1e5,
               "adjusted rate 95 % CI low": stats.gamma.ppf(0.025, dsr ** 2 / var, scale=var / dsr) * 1e5,
               "adjusted rate 95 % CI high": stats.gamma.ppf(0.975, (dsr + wm) ** 2 / (var + wm ** 2),
                                                             scale=(var + wm ** 2) / (dsr + wm)) * 1e5}
    refs(DS, f"{M}, birth order {order}", formula, "SciPy, Fay & Feuer formula", CLOSED)
    if HAS_R:
        r = r_run(f"""s <- tapply(births, mother_age, sum); k <- birth_order == "{order}"
            a <- epitools::ageadjust.direct(cases[k], births[k], stdpop = s[mother_age[k]]) * 1e5
            list(`crude rate` = a[["crude.rate"]], `adjusted rate` = a[["adj.rate"]],
                 `adjusted rate 95 % CI low` = a[["lci"]], `adjusted rate 95 % CI high` = a[["uci"]])""", down)
        refs(DS, f"{M}, birth order {order}", r, "R epitools ageadjust.direct", CLOSED)
for quantity, shown in {"adjusted rate": "92.3", "adjusted rate 95 % CI low": "80.4",
                        "adjusted rate 95 % CI high": "105.8"}.items():
    checks.printed(DS, f"{M}, birth order 1", quantity, shown, "epitools help, printed",
                   title="R epitools ageadjust.direct, help page example",
                   url="https://search.r-project.org/CRAN/refmans/epitools/html/ageadjust.direct.html")

M = "compare_rates (standard_rates=)"
pooled = down.groupby("mother_age")[["cases", "births"]].sum()
rates = (pooled["cases"] / pooled["births"] * 100000).to_dict()
for order in ORDERS:
    part = down[down["birth_order"] == order]
    observed = part["cases"].sum()
    expected = (part["births"] * part["mother_age"].map(rates) / 100000).sum()
    formula = {"observed / expected": observed / expected, "expected": expected,
               "observed / expected 95 % CI low": stats.chi2.ppf(0.025, 2 * observed) / 2 / expected,
               "observed / expected 95 % CI high": stats.chi2.ppf(0.975, 2 * (observed + 1)) / 2 / expected}
    refs(DS, f"{M}, birth order {order}", formula, "SciPy, exact Poisson formula", CLOSED)
    if HAS_R:
        r = r_run(f"""s <- tapply(cases, mother_age, sum); p <- tapply(births, mother_age, sum); k <- birth_order == "{order}"
            a <- epitools::ageadjust.indirect(cases[k], births[k], stdcount = s[mother_age[k]], stdpop = p[mother_age[k]])$sir
            pt <- stats::poisson.test(a[["observed"]], a[["exp"]])
            list(`observed / expected` = a[["sir"]], expected = a[["exp"]], `observed / expected 95 % CI low` = pt$conf.int[1],
                 `observed / expected 95 % CI high` = pt$conf.int[2], `p-value` = pt$p.value)""", down)
        refs(DS, f"{M}, birth order {order}", r, "R epitools + poisson.test", CLOSED)

# ---- 2 x 2 crossover: the bioequivalence reference datasets (Schutz, Labes & Fuglsang 2014)
be = load_csv("bioequivalence_reference.csv")
PUBLISHED = {"A": ("0.9509", "0.9076", "0.9962"), "B": ("0.7110", "0.5145", "0.9826"), "C": ("0.5856", "0.3941", "0.8703")}
SCHUTZ = "Schütz et al. (2014), Table II"
SCHUTZ_TITLE = ("Schütz H, Labes D, Fuglsang A (2014). Reference datasets for 2-treatment, 2-sequence, 2-period "
                "bioequivalence studies. AAPS J 16(6):1292-1297, doi:10.1208/s12248-014-9661-0")
for name, printed in PUBLISHED.items():
    DS, M = f"bioequivalence reference dataset {name}", "compare_averages_of (period=, crossover)"
    part = be[be["dataset"] == name].copy()
    checks.data(f"bioequivalence_{name}", part)
    for quantity, shown in zip(("ratio", "90 % CI low", "90 % CI high"), printed):
        checks.printed(DS, M, quantity, shown, SCHUTZ, title=SCHUTZ_TITLE, url="https://doi.org/10.1208/s12248-014-9661-0")
    part["log_auc"] = np.log(part["auc"])
    ols = smf.ols("log_auc ~ C(subject) + C(period) + C(formulation, Treatment('reference'))", data=part).fit()
    term = "C(formulation, Treatment('reference'))[T.test]"
    ci = ols.conf_int(alpha=0.10).loc[term]
    refs(DS, M, {"ratio": np.exp(ols.params[term]), "90 % CI low": np.exp(ci[0]), "90 % CI high": np.exp(ci[1]),
                 "within-subject CV": np.sqrt(np.exp(ols.scale) - 1)}, "statsmodels OLS", LEAST_SQUARES)
    if HAS_R:
        r = r_run("""s <- factor(subject); p <- factor(period); f <- relevel(factor(formulation), "reference")
            first <- tapply(as.character(f), list(s, p), function(v) v)[, 1]
            q <- factor(first[as.character(s)])
            m <- stats::lm(log(auc) ~ q + s + p + f)
            ci <- stats::confint(m, level = 0.90)["ftest", ]
            list(ratio = exp(coef(m)[["ftest"]]), `90 % CI low` = exp(ci[[1]]), `90 % CI high` = exp(ci[[2]]),
                 `within-subject CV` = sqrt(exp(summary(m)$sigma^2) - 1))""", part)
        refs(DS, M, r, "R lm (EMA model)", LEAST_SQUARES)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
