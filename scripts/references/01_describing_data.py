"""Reference values for notebooks/01_describing_data.ipynb, saved to references/01_describing_data.csv.

The nine NIST StRD univariate datasets: NIST certifies the mean, the standard deviation and the lag-1
autocorrelation of each (to 15 significant digits), so those are the references (certified, LRE of at least 9
digits). NIST certifies no median: the median is the one value run in a program here (R median(), relative
difference 1e-12), with the median of NumAcc4, which is exact by construction (10000000.2).
The data the notebook reads (one column y per dataset) are written to references/data/01_describing_data/.

    python scripts/references/01_describing_data.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

from checks import r_available, r_run, read_nist  # noqa: E402

NOTEBOOK = "01_describing_data"
NIST_UNIVARIATE = ["PiDigits", "Lottery", "Lew", "Mavro", "Michelso", "NumAcc1", "NumAcc2", "NumAcc3", "NumAcc4"]
MIN_LRE = 9
RELATIVE = 1e-12
NIST_URL = "https://www.itl.nist.gov/div898/strd/univ/data/{name}.dat"
R_MEDIAN = "https://search.r-project.org/R/refmans/stats/html/median.html"

ref = Recorder(NOTEBOOK)
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")


def certified(dataset, method, quantity, value, name):
    ref.compare(dataset, method, quantity, None, value, "NIST certified", MIN_LRE, "lre")
    ref.rows[-1]["url"] = NIST_URL.format(name=name)
    ref.rows[-1]["source title"] = f"NIST StRD univariate summary statistics: {name}"


for name in NIST_UNIVARIATE:
    d = read_nist(name)
    ref.data(f"nist_{name}", d.data[["y"]])
    ds = f"NIST {name}"
    certified(ds, "find_average_of", "mean", d.certified["mean"], name)
    certified(ds, "find_standard_deviation_of", "standard deviation", d.certified["standard deviation"], name)
    if len(d.data) < 4:
        ref.skip(ds, "find_autocorrelation_of", "lag-1 autocorrelation", "NIST certified",
                 f"simplydata gives no number: find_autocorrelation_of asks for at least 4 values and {name} has "
                 f"{len(d.data)}")
        ref.rows[-1]["url"] = NIST_URL.format(name=name)
    else:
        certified(ds, "find_autocorrelation_of", "lag-1 autocorrelation", d.certified["autocorrelation"], name)
    ref.compare(ds, "find_median_of", "median", None, r_run("list(m = median(y))", d.data)["m"], "R median",
                RELATIVE, "relative")
    ref.rows[-1]["url"] = R_MEDIAN

ref.compare("NIST NumAcc4", "find_median_of", "median", None, 10000000.2, "exact by construction (NIST NumAcc4)",
            MIN_LRE, "lre")
ref.rows[-1].update(kind="certified", url=NIST_URL.format(name="NumAcc4"),
                    **{"source title": "NIST NumAcc4: 500 values 0.1 below 10000000.2, 500 values 0.1 above, one "
                                       "equal to it"})
ref.save()
