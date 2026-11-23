"""Reference values for notebooks/05_survival.ipynb, saved to references/05_survival.csv.

The lung data (Loprinzi et al. 1994) ship with R's survival package, and simplydata computes these words as R's
survival package does, so R survival (Therneau) is the original. Its manual pages print no numbers for these models
on lung, so R is run once here: summary(survfit()) at every event time (closed formula, relative difference 1e-10;
whole columns compared by their largest relative difference), the median time with its interval, survdiff
(1e-10), and coxph with Efron's ties run to tight convergence (eps 1e-12) with cox.zph (iterative, 1e-6).
The data the notebook reads are written to references/data/05_survival/.

    python scripts/references/05_survival.py      (from the repository root; needs R with survival)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402

from checks import load_csv, r_available, r_eval_text, r_run  # noqa: E402

NOTEBOOK = "05_survival"
CLOSED_FORM = 1e-10
ITERATIVE = 1e-6
SURV = "https://search.r-project.org/CRAN/refmans/survival/html/{topic}.html"

ref = Recorder(NOTEBOOK, r_packages=("survival",))
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")
version = r_eval_text('as.character(packageVersion("survival"))')


def url(topic):
    return f"{SURV.format(topic=topic)} (version {version})"


lung = load_csv("lung.csv")
lung["died"] = (lung["status"] == 2).astype(int)        # status: 1 censored, 2 dead
lung["sex"] = lung["sex"].map({1: "male", 2: "female"})
ref.data("lung", lung)

# ---- Kaplan-Meier, each sex
for sex in ("female", "male"):
    group = lung[lung["sex"] == sex]
    r = r_run("s <- summary(survival::survfit(survival::Surv(time, died) ~ 1))\n"
              "list(time = s$time, survival = s$surv, low = s$lower, high = s$upper, risk = s$n.risk, "
              "events = s$n.event)", group[["time", "died"]].astype(float))
    ds = f"lung ({sex})"
    ref.compare(ds, "find_survival_over_time", "event times", None, len(np.atleast_1d(r["time"])),
                "R summary(survfit)", 0)
    ref.rows[-1]["url"] = url("summary.survfit")
    for quantity, key in (("time", "time"), ("survival", "survival"), ("95 % CI low", "low"),
                          ("95 % CI high", "high"), ("at risk", "risk"), ("events", "events")):
        ref.vector(ds, "find_survival_over_time", f"{quantity} at every event time", r[key], "R summary(survfit)",
                   CLOSED_FORM)
        ref.rows[-1]["url"] = url("summary.survfit")
    r = r_run("t <- summary(survival::survfit(survival::Surv(time, died) ~ 1))$table\n"
              "list(median = t[['median']], low = t[['0.95LCL']], high = t[['0.95UCL']])",
              group[["time", "died"]].astype(float))
    for quantity, key in (("median", "median"), ("median 95 % CI low", "low"), ("median 95 % CI high", "high")):
        ref.compare(ds, "find_survival_over_time (median)", quantity, None, r[key], "R survfit", CLOSED_FORM,
                    "relative")
        ref.rows[-1]["url"] = url("print.survfit")

# ---- the log-rank test
r = r_run("t <- survival::survdiff(survival::Surv(time, died) ~ sex)\n"
          "list(statistic = t$chisq, `p-value` = pchisq(t$chisq, length(t$n) - 1, lower.tail = FALSE), "
          "`expected events in female` = t$exp[1], `expected events in male` = t$exp[2])",
          lung[["time", "sex", "died"]].astype({"died": float}))
for quantity, value in r.items():
    ref.compare("lung", "compare_survival", quantity, None, value, "R survdiff", CLOSED_FORM, "relative")
    ref.rows[-1]["url"] = url("survdiff")

# ---- the Cox model
complete = lung[["time", "died", "age", "sex", "ph.ecog"]].dropna()
design = complete.assign(male=(complete["sex"] == "male").astype(float), died=complete["died"].astype(float))
r = r_run("f <- survival::coxph(survival::Surv(time, died) ~ age + male + ph.ecog, ties = 'efron', "
          "control = survival::coxph.control(eps = 1e-12, iter.max = 100))\n"
          "s <- summary(f)$coefficients; z <- survival::cox.zph(f)$table\n"
          "list(coefficient = s[, 1], `standard error` = s[, 3], `p-value` = s[, 5], zph = z[1:3, 3])",
          design[["time", "died", "age", "male", "ph.ecog"]])
for k, term in enumerate(("age", "sex: male vs female", "ph.ecog")):
    for quantity in ("coefficient", "standard error", "p-value"):
        ref.compare("lung", "explain_survival", f"{term} {quantity}", None, r[quantity][k], "R coxph", ITERATIVE,
                    "relative")
        ref.rows[-1]["url"] = url("coxph")
    ref.compare("lung", "explain_survival", f"{term} proportional hazards check p-value", None, r["zph"][k],
                "R cox.zph", ITERATIVE, "relative")
    ref.rows[-1]["url"] = url("cox.zph")
ref.save()
