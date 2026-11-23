"""Reference values for notebooks/06_curves_and_dose_response.ipynb, saved to references/06_curves_and_dose_response.csv.

In order of preference:
- certified: the NIST StRD nonlinear regression problems (27 of them, each from NIST's two starting points): the
  parameters, their standard errors and the residual standard deviation, to an LRE of at least 6 digits (the level of
  simplydata's own reliability study, after McCullough 1999);
- printed: Puromycin's Michaelis-Menten fit printed by Bates & Watts (1988), and the logistic fit of Bliss's beetles
  printed by Dobson & Barnett (2008), with the LD50 that follows from it, to half a unit of the last printed digit
  (0.01 for the LD50, which carries the rounding of both printed numbers);
- expected: what simplydata says of its own fit (BoxBOD from NIST's far start must be reported as not converged);
- run once here, only for the LD50 to more digits: R glm() on (killed, survived), 10^(-intercept / slope), run to
  tight convergence (relative difference 1e-6).
The data the notebook reads are written to references/data/06_curves_and_dose_response/.

    python scripts/references/06_curves_and_dose_response.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402

from checks import load_csv, r_available, r_run, read_nist  # noqa: E402

NOTEBOOK = "06_curves_and_dose_response"
MIN_LRE = 6
ITERATIVE = 1e-6
NIST_URL = "https://www.itl.nist.gov/div898/strd/nls/data/{name}.shtml"
NIST_PROBLEMS = ["Misra1a", "Chwirut2", "Chwirut1", "Lanczos3", "Gauss1", "Gauss2", "DanWood", "Misra1b", "Kirby2",
                 "Hahn1", "Nelson", "MGH17", "Lanczos1", "Lanczos2", "Gauss3", "Misra1c", "Misra1d", "Roszman1", "ENSO",
                 "MGH09", "Thurber", "BoxBOD", "Rat42", "MGH10", "Eckerle4", "Rat43", "Bennett5"]
BATES = "Bates & Watts (1988), printed"
BATES_TITLE = "Bates & Watts (1988), Nonlinear Regression Analysis and Its Applications, Wiley (the Puromycin example)"
BATES_URL = "https://doi.org/10.1002/9780470316757"
DOBSON = "Dobson & Barnett (2008), printed"
DOBSON_TITLE = "Dobson & Barnett (2008), An Introduction to Generalized Linear Models, 3rd ed., Chapman & Hall/CRC, ISBN 978-1-58488-950-2 (the beetle mortality example)"
DOBSON_URL = ""

ref = Recorder(NOTEBOOK)
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")

# ---- the NIST problems
for name in NIST_PROBLEMS:
    nist = read_nist(name)
    data = nist.data.copy()
    if name == "Nelson":
        data["log_y"] = np.log(data["y"])
    ref.data(f"nist_{name}", data)
    for start_number in (1, 2):
        method = f"fit_formula (start {start_number})"
        for b in nist.starts:
            for quantity, key in ((b, b), (f"{b} standard error", f"{b} standard error")):
                ref.compare(f"NIST {name}", method, quantity, None, nist.certified[key], "NIST certified", MIN_LRE,
                            "lre")
                ref.rows[-1]["url"] = NIST_URL.format(name=name.lower())
        ref.compare(f"NIST {name}", method, "residual standard deviation", None,
                    nist.certified["residual standard deviation"], "NIST certified", MIN_LRE, "lre")
        ref.rows[-1]["url"] = NIST_URL.format(name=name.lower())
    # the starting values, for the notebook (NIST gives two per parameter)
    ref.data(f"nist_{name}_starts", {"parameter": list(nist.starts), "start 1": [s[0] for s in nist.starts.values()],
                                     "start 2": [s[1] for s in nist.starts.values()]})

# ---- BoxBOD from the far start: what simplydata must say
ref.compare("NIST BoxBOD", "fit_formula (start 1)", "says not converged when the fit is wrong (1 = yes)", None, 1.0,
            "expected (A2)", 0)
ref.rows[-1].update(kind="expected", **{"source title": "FINDINGS A2: fit_formula must report converged = no when the "
                                                        "search stops in a flat region"})
ref.compare("NIST BoxBOD", "fit_formula (start 1)", "note says the search stopped in a flat region (1 = yes)", None,
            1.0, "expected (A2)", 0)
ref.rows[-1].update(kind="expected", **{"source title": 'the "note" column of the fit_formula result must say "the search '
                                                        'may have stopped in a flat region"'})

# ---- Puromycin and Bliss's beetles
ref.data("puromycin", load_csv("puromycin.csv"))
ref.printed("Puromycin (treated)", "fit_formula", "Vm", "212.7", BATES, BATES_TITLE, BATES_URL)
ref.printed("Puromycin (treated)", "fit_formula", "K", "0.06412", BATES, BATES_TITLE, BATES_URL)
bliss = load_csv("bliss.csv")
bliss["dose"] = 10 ** bliss["log10_dose"]
ref.data("bliss", bliss)
ref.printed("Bliss beetles", "explain (yes or no, out_of)", "intercept", "-60.717", DOBSON, DOBSON_TITLE, DOBSON_URL)
ref.printed("Bliss beetles", "explain (yes or no, out_of)", "slope", "34.270", DOBSON, DOBSON_TITLE, DOBSON_URL)
ref.compare("Bliss beetles", "fit_dose_response (yes or no, out_of)", "LD50", None, 10 ** (60.717 / 34.270),
            "Dobson & Barnett (2008), printed fit", 0.01)
ref.rows[-1].update(kind="printed", url=DOBSON_URL, **{"source title": DOBSON_TITLE + "; LD50 = 10^(60.717 / 34.270)"})
r = r_run("f <- glm(cbind(killed, beetles - killed) ~ log10_dose, family = binomial, "
          "control = glm.control(epsilon = 1e-15, maxit = 100))\n"
          "list(LD50 = 10^(-coef(f)[[1]] / coef(f)[[2]]))", bliss[["killed", "beetles", "log10_dose"]])
ref.compare("Bliss beetles", "fit_dose_response (yes or no, out_of)", "LD50", None, r["LD50"],
            "R glm, 10^(-intercept / slope)", ITERATIVE, "relative")
ref.rows[-1]["url"] = "https://search.r-project.org/R/refmans/stats/html/glm.html"
ref.save()
