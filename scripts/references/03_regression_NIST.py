"""Reference values for notebooks/03_regression_NIST.ipynb, saved to references/03_regression_NIST.csv.

In order of preference: the NIST StRD certified values of the eleven linear regression datasets (estimates, their
standard errors, R squared and the residual standard deviation; LRE of at least 9 digits); the fit Anscombe (1973)
printed for his four sets (within about one unit of the last printed digit, as the sets themselves agree with it);
and, only where nothing is printed or certified (mtcars and stackloss), R lm() run once here (relative difference
1e-9). The data the notebook reads are written to references/data/03_regression_NIST/.

    python scripts/references/03_regression_NIST.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

from checks import load_csv, r_available, r_run, read_nist  # noqa: E402

NOTEBOOK = "03_regression_NIST"
MIN_LRE = 9
RELATIVE = 1e-9
NIST_URL = "https://www.itl.nist.gov/div898/strd/lls/data/{name}.shtml"
ANSCOMBE = "Anscombe (1973), printed"
ANSCOMBE_TITLE = "Anscombe (1973), Graphs in statistical analysis, The American Statistician 27(1), 17-21"
ANSCOMBE_URL = "https://doi.org/10.1080/00031305.1973.10478966"
R_LM = "https://search.r-project.org/R/refmans/stats/html/lm.html"

ref = Recorder(NOTEBOOK)
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")


def certified(dataset, method, quantity, value, name):
    ref.compare(dataset, method, quantity, None, value, "NIST certified", MIN_LRE, "lre")
    ref.rows[-1]["url"] = NIST_URL.format(name=name)
    ref.rows[-1]["source title"] = f"NIST StRD linear regression: {name}"


# ---- a straight line: Norris
norris = read_nist("Norris")
ref.data("nist_Norris", norris.data)
for ours, key in {"intercept": "B0", "slope": "B1", "intercept standard error": "B0 standard error",
                  "slope standard error": "B1 standard error", "R squared": "R squared"}.items():
    certified("NIST Norris", "fit_line", ours, norris.certified[key], "Norris")

# ---- polynomials
POLYNOMIALS = {"Pontius": 2, "Wampler1": 5, "Wampler2": 5, "Wampler3": 5, "Wampler4": 5, "Wampler5": 5, "Filip": 10}
for name, degree in POLYNOMIALS.items():
    d = read_nist(name)
    ref.data(f"nist_{name}", d.data)
    terms = ["(intercept)", "x"] + [f"x**{k}" for k in range(2, degree + 1)]
    for k, term in enumerate(terms):
        certified(f"NIST {name}", "fit_curve", f"{term} estimate", d.certified[f"B{k}"], name)
        certified(f"NIST {name}", "fit_curve", f"{term} standard error", d.certified[f"B{k} standard error"], name)
    certified(f"NIST {name}", "fit_curve", "R squared", d.certified["R squared"], name)

# ---- several predictors: Longley
longley = read_nist("Longley")
ref.data("nist_Longley", longley.data)
for k, term in enumerate(["(intercept)"] + [f"x{i}" for i in range(1, 7)]):
    certified("NIST Longley", "explain", f"{term} effect", longley.certified[f"B{k}"], "Longley")
    certified("NIST Longley", "explain", f"{term} standard error", longley.certified[f"B{k} standard error"],
              "Longley")

# ---- a line through 0
for name in ("NoInt1", "NoInt2"):
    d = read_nist(name)
    ref.data(f"nist_{name}", d.data)
    for ours, key in (("b1", "B1"), ("standard error of b1", "B1 standard error"),
                      ("residual standard error", "residual standard deviation")):
        certified(f"NIST {name}", "fit_formula (y = b1 x)", ours, d.certified[key], name)

# ---- Anscombe's quartet: one printed fit for all four sets
for i in (1, 2, 3, 4):
    for quantity, value, tolerance in (("slope", 0.5, 0.001), ("intercept", 3.0, 0.01),
                                       ("slope standard error", 0.118, 0.001), ("R squared", 0.667, 0.001)):
        ref.compare(f"Anscombe set {i}", "fit_line", quantity, None, value, ANSCOMBE, tolerance)
        ref.rows[-1].update(kind="printed", url=ANSCOMBE_URL, **{"source title": ANSCOMBE_TITLE,
                                                                  "digits as printed": f"{value:g}"})

# ---- mtcars and stackloss: R lm
for dataset, (table, outcome, using) in {
        "mtcars": (load_csv("mtcars.csv"), "mpg", ["wt", "hp", "am"]),
        "stackloss": (load_csv("stackloss.csv"), "stack_loss", ["air_flow", "water_temp", "acid_conc"])}.items():
    r = r_run(f"f <- summary(lm({outcome} ~ {' + '.join(using)})); "
              "list(effect = coef(f)[, 1], se = coef(f)[, 2], t = coef(f)[, 3], p = coef(f)[, 4])",
              table[[outcome] + using])
    for k, term in enumerate(["(intercept)"] + using):
        for column, key in (("effect", "effect"), ("standard error", "se"), ("statistic", "t"), ("p-value", "p")):
            ref.compare(dataset, "explain", f"{term} {column}", None, r[key][k], "R lm", RELATIVE, "relative")
            ref.rows[-1]["url"] = R_LM
ref.save()
