"""Reference values for notebooks/07_distributions.ipynb, saved to references/07_distributions.csv.

In order of preference: the printed tables of the NIST/SEMATECH e-Handbook of Statistical Methods, section 1.3.6.7
(the standard normal, t, chi-square and F tables), to half a unit of their last printed digit; and, for the grid of
values the tables do not cover, R (pnorm, pt, pchisq, pf and their quantile functions) run once here, to a relative
difference of 1e-12.

    python scripts/references/07_distributions.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

from checks import r_available, r_run  # noqa: E402

NOTEBOOK = "07_distributions"
RELATIVE = 1e-12
HANDBOOK = "NIST/SEMATECH e-Handbook 1.3.6.7, printed"
PAGE = "https://www.itl.nist.gov/div898/handbook/eda/section3/eda367{n}.htm"
RMAN = "https://search.r-project.org/R/refmans/stats/html/{topic}.html"

ref = Recorder(NOTEBOOK)
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")

# The normal table gives the area from 0 to 1.96 as 0.47500, so the chance above 1.96 is 0.02500.
for dataset, method, quantity, shown, n, title in (
        ("standard normal table", "find_chance_that_value", "chance above 1.96", "0.02500", 1,
         "Table 1.3.6.7.1, cumulative distribution function of the standard normal: 0.47500 from 0 to 1.96"),
        ("t table", "find_cut_off", "upper 2.5 % point, 10 df", "2.228", 2,
         "Table 1.3.6.7.2, upper critical values of Student's t distribution"),
        ("chi-square table", "find_cut_off", "upper 5 % point, 10 df", "18.307", 4,
         "Table 1.3.6.7.4, critical values of the chi-square distribution"),
        ("F table", "find_cut_off", "upper 5 % point, 5 and 10 df", "3.326", 3,
         "Table 1.3.6.7.3, upper critical values of the F distribution")):
    ref.printed(dataset, method, quantity, shown, HANDBOOK, f"NIST/SEMATECH e-Handbook of Statistical Methods, {title}",
                PAGE.format(n=n))

CASES = {"normal(10, 2)": ("pnorm({x}, 10, 2)", "qnorm({p}, 10, 2)", "Normal", [4.0, 7.5, 10.0, 12.3, 17.0]),
         "t(3)": ("pt({x}, 3)", "qt({p}, 3)", "TDist", [-6.0, -1.0, 0.0, 0.7, 4.5]),
         "chi-square(4)": ("pchisq({x}, 4)", "qchisq({p}, 4)", "Chisquare", [0.2, 1.0, 3.36, 9.49, 20.0]),
         "F(3, 12)": ("pf({x}, 3, 12)", "qf({p}, 3, 12)", "Fdist", [0.1, 0.5, 1.0, 3.49, 8.0])}
for label, (r_p, r_q, topic, values) in CASES.items():
    for x in values:
        ref.compare(label, "find_chance_that_value", f"chance below {x:g}", None,
                    r_run(f"list(p = {r_p.format(x=x)})")["p"], "R", RELATIVE, "relative")
        ref.rows[-1]["url"] = RMAN.format(topic=topic)
    for p in [0.001, 0.025, 0.5, 0.9, 0.999]:
        ref.compare(label, "find_cut_off", f"cut-off with {p:g} below", None,
                    r_run(f"list(q = {r_q.format(p=p)})")["q"], "R", RELATIVE, "relative")
        ref.rows[-1]["url"] = RMAN.format(topic=topic)
ref.save()
