"""Reference values for notebooks/04_counts_and_yes_no.ipynb, saved to references/04_counts_and_yes_no.csv.

In order of preference:
- printed: the Poisson and negative binomial fits of the NMES 1988 office visits that Zeileis, Kleiber & Jackman
  (2008) printed (the Poisson estimates and standard errors with 6 decimals in the R summary on page 13; the negative
  binomial estimates and standard errors with 3 decimals in their Table 2), to half a unit of the last printed digit;
- exact by construction: the two risks of the 2 x 2 table, 30 of 74 smokers and 29 of 115 non-smokers;
- run once here, only where nothing is printed: R glm (logistic; the p-values of the count models) and MASS::glm.nb,
  run to tight convergence (glm.control(epsilon = 1e-15): R takes its standard errors from the weights of the step
  before the last, so with its default 1e-8 it can stop one step early and p-values move by 1e-6), relative
  difference 1e-6; R chisq.test, fisher.test, binom.test and prop.test (closed formulas, 1e-10); SciPy for the G-test
  and statsmodels for Newcombe's interval of the risk difference and the Wald intervals of the risk and odds ratios,
  which base R does not give (1e-10).
The data the notebook reads are written to references/data/04_counts_and_yes_no/.

    python scripts/references/04_counts_and_yes_no.py      (from the repository root; needs R with MASS)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.contingency_tables import Table2x2  # noqa: E402
from statsmodels.stats.proportion import confint_proportions_2indep  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

NOTEBOOK = "04_counts_and_yes_no"
ITERATIVE = 1e-6
CLOSED_FORM = 1e-10
R_TIGHT = "control = glm.control(epsilon = 1e-15, maxit = 100)"
ZKJ = "Zeileis, Kleiber & Jackman (2008), printed"
ZKJ_TITLE = "Zeileis, Kleiber & Jackman (2008), Regression models for count data in R, J. Stat. Softw. 27(8)"
ZKJ_URL = "https://doi.org/10.18637/jss.v027.i08"
RMAN = "https://search.r-project.org/R/refmans/stats/html/{topic}.html"
SCIPY = "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.{name}.html"
SM = "https://www.statsmodels.org/v0.15.0/generated/statsmodels.stats.{name}.html"

ref = Recorder(NOTEBOOK, r_packages=("MASS",))
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")


def r_name(term: str) -> str:
    """simplydata's 'race (white)' is R's 'racewhite'; '(intercept)' is '(Intercept)'."""
    if term == "(intercept)":
        return "(Intercept)"
    if term.endswith(")") and " (" in term:
        column, level = term[:-1].split(" (", 1)
        return f"{column}{level}"
    return term


R_COEFFICIENTS = """s <- coef(summary(f))
list(estimate = s[, 1], `standard error` = s[, 2], `p-value` = s[, 4], names = rownames(s))"""


def r_model(dataset, method, code, data, terms, quantities, url, packages=()):
    r = r_run(code + "\n" + R_COEFFICIENTS, data, packages)
    for term in terms:
        k = r["names"].index(r_name(term))
        for quantity in quantities:
            ref.compare(dataset, method, f"{term} {quantity}", None, r[quantity][k], "R", ITERATIVE, "relative")
            ref.rows[-1]["url"] = url
    return r


# ---- logistic regression on birthwt
birthwt = load_csv("birthwt.csv")
birthwt["race"] = birthwt["race"].map({1: "white", 2: "black", 3: "other"})
ref.data("birthwt", birthwt)
r_model("birthwt", "explain (yes or no)", f"f <- glm(low ~ age + lwt + race + smoke, family = binomial, {R_TIGHT})",
        birthwt, ["(intercept)", "age", "lwt", "race (other)", "race (white)", "smoke"],
        ["estimate", "standard error", "p-value"], RMAN.format(topic="glm"))

# ---- the count models on NMES 1988 (the printed fits use the names of the paper's data)
nmes = load_csv("nmes1988.csv")
formula = "visits ~ hospital + health + chronic + gender + school + insurance"
TERMS = {"(intercept)": "(Intercept)", "hospital": "hosp", "health (poor)": "healthpoor",
         "health (excellent)": "healthexcellent", "chronic": "numchron", "gender (male)": "gendermale",
         "school": "school", "insurance (yes)": "privinsyes"}
POISSON = {"(Intercept)": ("1.028874", "0.023785"), "hosp": ("0.164797", "0.005997"),
           "healthpoor": ("0.248307", "0.017845"), "healthexcellent": ("-0.361993", "0.030304"),
           "numchron": ("0.146639", "0.004580"), "gendermale": ("-0.112320", "0.012945"),
           "school": ("0.026143", "0.001843"), "privinsyes": ("0.201687", "0.016860")}
NEGBIN = {"(Intercept)": ("0.929", "0.055"), "hosp": ("0.218", "0.020"), "healthpoor": ("0.305", "0.049"),
          "healthexcellent": ("-0.342", "0.061"), "numchron": ("0.175", "0.012"), "gendermale": ("-0.126", "0.031"),
          "school": ("0.027", "0.004"), "privinsyes": ("0.224", "0.039")}
for method, printed, where in (("explain (counts)", POISSON, "summary(fm_pois), page 13"),
                               ("explain (counts that vary a lot)", NEGBIN, "Table 2, column fm_nbin, page 17")):
    for term, name in TERMS.items():
        estimate, se = printed[name]
        for quantity, shown in (("estimate", estimate), ("standard error", se)):
            ref.printed("NMES 1988", method, f"{term} {quantity}", shown, ZKJ, f"{ZKJ_TITLE}, {where}", ZKJ_URL)
r_model("NMES 1988", "explain (counts)", f"f <- glm({formula}, family = poisson, {R_TIGHT})", nmes, list(TERMS),
        ["p-value"], RMAN.format(topic="glm"))
r_model("NMES 1988", "explain (counts that vary a lot)", f"f <- MASS::glm.nb({formula}, {R_TIGHT})", nmes,
        list(TERMS), ["p-value"], "https://search.r-project.org/CRAN/refmans/MASS/html/glm.nb.html", ("MASS",))

# ---- the 2 x 2 table of smoking against low birth weight
counts = birthwt.assign(low=birthwt["low"].map({1: "yes", 0: "no"}),
                        smoke=birthwt["smoke"].map({1: "smoker", 0: "non-smoker"}))
ref.data("counts", counts)
for method, code, source, topic in (
        ("compare_counts_of (chi-square)",
         "t <- chisq.test(table(smoke, low)); list(statistic = t$statistic, `p-value` = t$p.value)", "R chisq.test",
         "chisq.test"),
        ("compare_counts_of (fisher)", "t <- fisher.test(table(smoke, low)); list(`p-value` = t$p.value)",
         "R fisher.test", "fisher.test")):
    for quantity, value in r_run(code, counts).items():
        ref.compare("birthwt", method, quantity, None, value, source, CLOSED_FORM, "relative")
        ref.rows[-1]["url"] = RMAN.format(topic=topic)
table = pd.crosstab(counts["smoke"], counts["low"]).to_numpy()
g = stats.chi2_contingency(table, correction=False, lambda_="log-likelihood")
for quantity, value in (("statistic", g.statistic), ("p-value", g.pvalue)):
    ref.compare("birthwt", "compare_counts_of (g)", quantity, None, value,
                "SciPy chi2_contingency(correction=False, lambda_='log-likelihood')", CLOSED_FORM, "relative")
    ref.rows[-1]["url"] = SCIPY.format(name="chi2_contingency")

# ---- risks: 30 of 74 smokers, 29 of 115 non-smokers
low, high = confint_proportions_2indep(30, 74, 29, 115, method="newcomb", compare="diff")
t2 = Table2x2(np.array([[30, 74 - 30], [29, 115 - 29]]))
for measure, (estimate, lo, hi, name, url) in {
        "risk difference": (30 / 74 - 29 / 115, low, high, "statsmodels confint_proportions_2indep(newcomb)",
                            SM.format(name="proportion.confint_proportions_2indep")),
        "risk ratio": (t2.riskratio, *t2.riskratio_confint(), "statsmodels Table2x2",
                       SM.format(name="contingency_tables.Table2x2")),
        "odds ratio": (t2.oddsratio, *t2.oddsratio_confint(), "statsmodels Table2x2",
                       SM.format(name="contingency_tables.Table2x2"))}.items():
    for quantity, value in (("estimate", estimate), ("95 % CI low", lo), ("95 % CI high", hi)):
        ref.compare("birthwt", "compare_risks", f"{measure} {quantity}", None, value, name, CLOSED_FORM, "relative")
        ref.rows[-1]["url"] = url
for label, events, total in (("smoker", 30, 74), ("non-smoker", 29, 115)):
    ref.compare("birthwt", "compare_risks", f"{label} risk", None, events / total,
                f"exact by construction ({events} of {total})", CLOSED_FORM, "relative")
    ref.rows[-1]["kind"] = "printed"
    r = r_run(f"b <- binom.test({events}, {total})$conf.int; list(low = b[1], high = b[2])")
    for quantity, value in (("risk 95 % CI low", r["low"]), ("risk 95 % CI high", r["high"])):
        ref.compare("birthwt", "compare_risks", f"{label} {quantity}", None, value, "R binom.test", CLOSED_FORM,
                    "relative")
        ref.rows[-1]["url"] = RMAN.format(topic="binom.test")
ref.compare("birthwt", "compare_risks", "p-value", None, r_run("list(p = prop.test(c(30, 29), c(74, 115))$p.value)")["p"],
            "R prop.test", CLOSED_FORM, "relative")
ref.rows[-1]["url"] = RMAN.format(topic="prop.test")
ref.save()
