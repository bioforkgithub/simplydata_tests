"""Reference values for notebooks/09_survival_and_person_time.ipynb, saved to references/09_survival_and_person_time.csv.

Nothing is printed for these datasets in a form the notebook can use (the textbook tables they come from print
results rounded to 3 or 4 digits and were themselves made with these R packages), so every reference here is the
original library run once: R survival (survfit, survdiff, coxph, survreg), cmprsk (cuminc, crr), survRM2 (rmst2),
KMsurv (lifetab), adjustedCurves (adjustedsurv); statsmodels (SurvfuncRight, CumIncidenceRight, survdiff, PHReg);
lifelines (KaplanMeierFitter, CoxPHFitter, CoxTimeVaryingFitter, WeibullAFTFitter); and the closed formulas
written out (person-time, the actuarial life table).

Tolerances (relative difference): 1e-10 for results from closed formulas, 1e-6 for models fitted by iteration.
lifelines is run with tighter convergence than its defaults (CoxPHFitter precision=1e-12; WeibullAFTFitter L-BFGS-B
with ftol=1e-15, gtol=1e-12), because its defaults stop 1e-6 to 1e-5 short of the solution.

The tables the notebook reads (the data files of this repository with the recoding the notebook needs: yes/no for
died, male/female for sex ...) are written to references/data/09_survival_and_person_time/.

    python scripts/references/09_survival_and_person_time.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import CoxPHFitter, CoxTimeVaryingFitter, KaplanMeierFitter, WeibullAFTFitter  # noqa: E402
from lifelines.utils import restricted_mean_survival_time  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402
from statsmodels.duration.survfunc import CumIncidenceRight, SurvfuncRight, survdiff  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "09_survival_and_person_time"
R_PACKAGES = ("survival", "cmprsk", "survRM2", "KMsurv", "adjustedCurves")
checks = Recorder(STEM, R_PACKAGES)
HAS_R = r_available()
CLOSED = 1e-10
ITERATIVE = 1e-6
TIGHT_COX = {"precision": 1e-12, "max_steps": 500}  # lifelines' default stops at 1e-7


def at(times, values, t):
    """The value of a step function (given at its jump times) at time t."""
    return np.asarray(values)[np.searchsorted(np.asarray(times), t, side="right") - 1]


def ref(dataset, method, quantity, value, source, tolerance):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, "relative")


lung = load_csv("lung.csv")
lung["died"] = np.where(lung["status"] == 2, "yes", "no")
lung["sex"] = lung["sex"].map({1: "male", 2: "female"})
checks.data("lung", lung)

# ---- person-time from dates: UDCA trial
udca = load_csv("udca.csv")
checks.data("udca", udca)
entry, last, death = (pd.to_datetime(udca[c]) for c in ("entry.dt", "last.dt", "death.dt"))
end = pd.concat([death, last], axis=1).min(axis=1)
years = (end - entry).dt.days / 365.25
python = {"person-years": years.sum(), "deaths": float(((death <= last) & death.notna()).sum()), "rows": float(len(udca))}
M = "add_time_and_event_columns"
for quantity, value in python.items():
    ref("UDCA trial", M, quantity, value, "pandas dates (written out)", CLOSED)
checks.vector("UDCA trial", M, "time of each row (years)", years, "pandas dates (written out)", 1e-12)
if HAS_R:
    r = r_run("""s <- as.Date(entry.dt); l <- as.Date(last.dt); d <- as.Date(ifelse(death.dt == "nan", NA, death.dt))
e <- pmin(d, l, na.rm = TRUE)
list(`person-years` = sum(as.numeric(difftime(e, s, units = "days"))) / 365.25, deaths = sum(!is.na(d) & d <= l),
     rows = length(s))""", udca[["entry.dt", "last.dt", "death.dt"]].fillna("nan"))
    for quantity, value in r.items():
        ref("UDCA trial", M, quantity, value, "R difftime", CLOSED)

# ---- survival at chosen times: NCCTG lung
M = "find_survival_over_time (at_time_points)"
if HAS_R:
    r = r_run("""s <- summary(survival::survfit(survival::Surv(time, died == "yes") ~ sex), times = c(180, 365))
list(surv = s$surv, lower = s$lower, upper = s$upper, risk = s$n.risk, events = ave(s$n.event, s$strata, FUN = cumsum))""",
              lung)
for i, (sex, t) in enumerate((("female", 180), ("female", 365), ("male", 180), ("male", 365))):
    group = lung[lung["sex"] == sex]
    sm_fit = SurvfuncRight(group["time"], group["died"] == "yes")
    ref("NCCTG lung", M, f"{sex} at {t}: survival", at(sm_fit.surv_times, sm_fit.surv_prob, t),
        "statsmodels SurvfuncRight", CLOSED)
    km = KaplanMeierFitter().fit(group["time"], group["died"] == "yes")
    ref("NCCTG lung", M, f"{sex} at {t}: survival", float(km.survival_function_at_times(t).iloc[0]),
        "lifelines KaplanMeierFitter", CLOSED)
    if HAS_R:
        for quantity, key in (("survival", "surv"), ("95 % CI low", "lower"), ("95 % CI high", "upper"),
                              ("at risk", "risk"), ("events by then", "events")):
            ref("NCCTG lung", M, f"{sex} at {t}: {quantity}", r[key][i], "R survfit", CLOSED)

# ---- restricted mean survival time
M = "find_survival_over_time (average_time_up_to)"
for sex in ("female", "male"):
    group = lung[lung["sex"] == sex]
    km = KaplanMeierFitter().fit(group["time"], group["died"] == "yes")
    ref("NCCTG lung", M, f"{sex}: average days up to 365", restricted_mean_survival_time(km, t=365),
        "lifelines restricted_mean_survival_time", CLOSED)
if HAS_R:
    r = r_run("""t <- summary(survival::survfit(survival::Surv(time, died == "yes") ~ sex), rmean = 365)$table
a <- survRM2::rmst2(time, as.numeric(died == "yes"), as.numeric(sex == "male"), tau = 365)
list(rmean = unname(t[, "rmean"]), se = unname(t[, "se(rmean)"]),
     rm2 = c(a$RMST.arm0$rmst[1], a$RMST.arm1$rmst[1]), rm2se = c(a$RMST.arm0$rmst[2], a$RMST.arm1$rmst[2]),
     diff = unname(a$unadjusted.result[1, ]), ratio = unname(a$unadjusted.result[2, ]))""", lung)
    for i, sex in enumerate(("female", "male")):
        for quantity, key, name in (("average days up to 365", "rmean", "R survfit (rmean)"),
                                    ("standard error", "se", "R survfit (rmean)"),
                                    ("average days up to 365", "rm2", "R survRM2 rmst2"),
                                    ("standard error", "rm2se", "R survRM2 rmst2")):
            ref("NCCTG lung", M, f"{sex}: {quantity}", r[key][i], name, CLOSED)
    M = "compare_survival (average_time_up_to)"
    for quantity, value in (("difference", r["diff"][0]), ("difference 95 % CI low", r["diff"][1]),
                            ("difference 95 % CI high", r["diff"][2]), ("difference p-value", r["diff"][3]),
                            ("ratio", r["ratio"][0]), ("ratio 95 % CI low", r["ratio"][1]),
                            ("ratio 95 % CI high", r["ratio"][2])):
        ref("NCCTG lung", M, quantity, value, "R survRM2 rmst2", CLOSED)

# ---- actuarial life table: bfeed
bfeed = load_csv("bfeed.csv")
checks.data("bfeed", bfeed)
CUTS = [0, 2, 3, 5, 7, 11, 17, 25, 37, 53]
k = np.searchsorted(CUTS, bfeed["duration"], side="right") - 1
events = np.bincount(k[bfeed["delta"] == 1], minlength=10).astype(float)
lost = np.bincount(k[bfeed["delta"] == 0], minlength=10).astype(float)
at_start = len(bfeed) - np.concatenate([[0], np.cumsum(events + lost)[:-1]])
effective = at_start - lost / 2
q = events / effective
survival = np.concatenate([[1.0], np.cumprod(1 - q)[:-1]])
with np.errstate(divide="ignore", invalid="ignore"):
    greenwood = survival * np.sqrt(np.concatenate([[0.0], np.cumsum(q / (effective - events))[:-1]]))
    width = np.diff(CUTS + [np.nan])
    hazard = events / (width * (effective - events / 2))
    hazard_se = np.sqrt((hazard ** 2 / events) * (1 - (hazard * width / 2) ** 2))
python = {"at start": at_start, "lost": lost, "events": events, "effective at risk": effective,
          "survival at the start": survival, "standard error of survival": greenwood, "hazard": hazard,
          "standard error of hazard": hazard_se}
if HAS_R:
    r = r_run("""tis <- c(0, 2, 3, 5, 7, 11, 17, 25, 37, 53, NA)
k <- findInterval(duration, tis[1:10])
t <- KMsurv::lifetab(tis, length(duration), tabulate(k[delta == 0], 10), tabulate(k[delta == 1], 10))
list(`at start` = t$nsubs, lost = t$nlost, `effective at risk` = t$nrisk, events = t$nevent,
     `survival at the start` = t$surv, `standard error of survival` = t$se.surv, hazard = t$hazard,
     `standard error of hazard` = t$se.hazard)""", bfeed)
M = "find_survival_over_time (life table)"
for i, start in enumerate(CUTS):
    for quantity, values in python.items():
        if np.isfinite(values[i]):
            ref("bfeed (KMsurv)", M, f"from week {start}: {quantity}", values[i], "actuarial formulas (Python)", CLOSED)
        if HAS_R and np.isfinite(r[quantity][i]):
            ref("bfeed (KMsurv)", M, f"from week {start}: {quantity}", r[quantity][i], "R KMsurv lifetab", CLOSED)

# ---- adjusted curves
if HAS_R:
    r = r_run("""d <- na.omit(data.frame(time = time, event = as.numeric(died == "yes"), group = factor(sex), age = age,
                         ecog = ph.ecog))
m <- survival::coxph(survival::Surv(time, event) ~ group + age + ecog, data = d, x = TRUE)
a <- adjustedCurves::adjustedsurv(data = d, variable = "group", ev_time = "time", event = "event", method = "direct",
                                  outcome_model = m, conf_int = TRUE, times = c(180, 365))
list(group = as.character(a$adj$group), time = a$adj$time, surv = a$adj$surv, low = a$adj$ci_lower, high = a$adj$ci_upper)""",
              lung[["time", "died", "sex", "age", "ph.ecog"]], packages=("survival",))
    for i, (sex, t) in enumerate(zip(r["group"], r["time"])):
        for quantity, key in (("survival", "surv"), ("95 % CI low", "low"), ("95 % CI high", "high")):
            ref("NCCTG lung", "find_survival_over_time (allowing_for)", f"{sex} at {t:g}: {quantity}", r[key][i],
                "R adjustedCurves adjustedsurv", ITERATIVE)

# ---- delayed entry: Channing House
channing = load_csv("channing.csv")
channing = channing[channing["age"] > channing["ageentry"]].copy()
channing["gender"] = np.where(channing["gender"] == 1, "male", "female")
checks.data("channing", channing)
M = "find_survival_over_time (followed_from)"
sm_fit = SurvfuncRight(channing["age"], channing["death"], entry=channing["ageentry"])
ref("Channing House", M, "survival at 960 months", at(sm_fit.surv_times, sm_fit.surv_prob, 960),
    "statsmodels SurvfuncRight", CLOSED)
km = KaplanMeierFitter().fit(channing["age"], channing["death"], entry=channing["ageentry"])
ref("Channing House", M, "survival at 960 months", float(km.survival_function_at_times(960).iloc[0]),
    "lifelines KaplanMeierFitter", CLOSED)
chi, p = survdiff(channing["age"], channing["death"], channing["gender"], entry=channing["ageentry"])
ref("Channing House", "compare_survival (followed_from)", "chi-square", chi, "statsmodels survdiff", CLOSED)
ref("Channing House", "compare_survival (followed_from)", "p-value", p, "statsmodels survdiff", CLOSED)
if HAS_R:
    r = r_run("""s <- summary(survival::survfit(survival::Surv(ageentry, age, death) ~ 1), times = 960)
m <- survival::coxph(survival::Surv(ageentry, age, death) ~ gender, ties = "exact")
list(surv = s$surv, low = s$lower, high = s$upper, risk = s$n.risk, score = unname(m$score))""", channing)
    for quantity, key in (("survival at 960 months", "surv"), ("95 % CI low", "low"), ("95 % CI high", "high"),
                          ("at risk", "risk")):
        ref("Channing House", M, quantity, r[key], "R survfit", CLOSED)
    ref("Channing House", "compare_survival (followed_from)", "chi-square", r["score"], "R coxph score test", CLOSED)

# ---- competing events: mgus2
mgus = load_csv("mgus2.csv")
mgus["etime"] = np.where(mgus["pstat"] == 1, mgus["ptime"], mgus["futime"])
mgus["outcome"] = np.where(mgus["pstat"] == 1, "pcm", np.where(mgus["death"] == 1, "death", "censored"))
mgus = mgus[["age", "sex", "etime", "outcome"]]
checks.data("mgus", mgus)
M = "find_risk_over_time (competing_with)"
codes = {"censored": 0, "pcm": 1, "death": 2}
for sex in ("F", "M"):
    group = mgus[mgus["sex"] == sex]
    cif = CumIncidenceRight(group["etime"], group["outcome"].map(codes))
    for j, event in enumerate(("pcm", "death")):
        ref("mgus2", M, f"{sex}, {event} by 120 months: risk", at(cif.times, cif.cinc[j], 120),
            "statsmodels CumIncidenceRight", CLOSED)
if HAS_R:
    r = r_run("""st <- factor(outcome, levels = c("censored", "pcm", "death"))
s <- summary(survival::survfit(survival::Surv(etime, st) ~ sex), times = 120)
list(risk = c(s$pstate[, 2], s$pstate[, 3]), low = c(s$lower[, 2], s$lower[, 3]), high = c(s$upper[, 2], s$upper[, 3]))""",
              mgus)
    for i, (event, sex) in enumerate((("pcm", "F"), ("pcm", "M"), ("death", "F"), ("death", "M"))):
        for quantity, key in (("risk", "risk"), ("95 % CI low", "low"), ("95 % CI high", "high")):
            ref("mgus2", M, f"{sex}, {event} by 120 months: {quantity}", r[key][i], "R survfit (Aalen-Johansen)", CLOSED)

TERMS = ["age", "sex: M vs F"]
x = pd.DataFrame({"age": mgus["age"], "male": (mgus["sex"] == "M").astype(float)})
ph = PHReg(mgus["etime"], x, status=(mgus["outcome"] == "pcm").astype(float), ties="efron").fit()
cph = CoxPHFitter().fit(x.assign(T=mgus["etime"], E=(mgus["outcome"] == "pcm").astype(int)), "T", "E", fit_options=TIGHT_COX)
M = "explain_survival (competing_with, cause-specific)"
for j, term in enumerate(TERMS):
    for name, coef, se in (("statsmodels PHReg", ph.params[j], ph.bse[j]),
                           ("lifelines CoxPHFitter", cph.params_.iloc[j], cph.standard_errors_.iloc[j])):
        ref("mgus2", M, f"{term}: coefficient", coef, name, ITERATIVE)
        ref("mgus2", M, f"{term}: standard error", se, name, ITERATIVE)
if HAS_R:
    r = r_run("""s <- match(outcome, c("pcm", "death")); s[is.na(s)] <- 0
t <- cmprsk::cuminc(etime, s, sex, cencode = 0)$Tests
m <- survival::coxph(survival::Surv(etime, outcome == "pcm") ~ age + sex)
f <- cmprsk::crr(etime, s, cbind(age, as.numeric(sex == "M")), failcode = 1, cencode = 0)
list(chi = t["1", "stat"], p = t["1", "pv"], cs = unname(coef(m)), cs_se = unname(sqrt(diag(vcov(m)))),
     fg = unname(f$coef), fg_se = unname(sqrt(diag(f$var))))""", mgus)
    ref("mgus2", "compare_survival (Gray's test)", "chi-square", r["chi"], "R cmprsk cuminc", CLOSED)
    ref("mgus2", "compare_survival (Gray's test)", "p-value", r["p"], "R cmprsk cuminc", CLOSED)
    for j, term in enumerate(TERMS):
        ref("mgus2", M, f"{term}: coefficient", r["cs"][j], "R coxph", ITERATIVE)
        ref("mgus2", M, f"{term}: standard error", r["cs_se"][j], "R coxph", ITERATIVE)
        ref("mgus2", "explain_survival (fine-gray)", f"{term}: coefficient", r["fg"][j], "R cmprsk crr", ITERATIVE)
        ref("mgus2", "explain_survival (fine-gray)", f"{term}: standard error", r["fg_se"][j], "R cmprsk crr", ITERATIVE)

# ---- stratified log-rank and Cox: veteran
veteran = load_csv("veteran.csv")
veteran["trt"] = veteran["trt"].map({1: "standard", 2: "test"})
checks.data("veteran", veteran)
chi, p = survdiff(veteran["time"], veteran["status"], veteran["trt"], strata=veteran["celltype"])
ref("veteran", "compare_survival (allowing_for)", "chi-square", chi, "statsmodels survdiff", CLOSED)
ref("veteran", "compare_survival (allowing_for)", "p-value", p, "statsmodels survdiff", CLOSED)
TERMS = ["trt: test vs standard", "karno"]
x = pd.DataFrame({"test": (veteran["trt"] == "test").astype(float), "karno": veteran["karno"].astype(float)})
ph = PHReg(veteran["time"], x, status=veteran["status"], ties="efron", strata=veteran["celltype"]).fit()
cph = CoxPHFitter().fit(x.assign(T=veteran["time"], E=veteran["status"], cell=veteran["celltype"]), "T", "E",
                        strata=["cell"], fit_options=TIGHT_COX)
M = "explain_survival (allowing_for)"
for j, term in enumerate(TERMS):
    for name, coef, se in (("statsmodels PHReg", ph.params[j], ph.bse[j]),
                           ("lifelines CoxPHFitter", cph.params_.iloc[j], cph.standard_errors_.iloc[j])):
        ref("veteran", M, f"{term}: coefficient", coef, name, ITERATIVE)
        ref("veteran", M, f"{term}: standard error", se, name, ITERATIVE)
if HAS_R:
    r = r_run("""s <- survival::survdiff(survival::Surv(time, status) ~ trt + survival::strata(celltype))
m <- survival::coxph(survival::Surv(time, status) ~ trt + karno + survival::strata(celltype))
list(chi = s$chisq, p = stats::pchisq(s$chisq, 1, lower.tail = FALSE), coef = unname(coef(m)),
     se = unname(sqrt(diag(vcov(m)))))""", veteran)
    ref("veteran", "compare_survival (allowing_for)", "chi-square", r["chi"], "R survdiff", CLOSED)
    ref("veteran", "compare_survival (allowing_for)", "p-value", r["p"], "R survdiff", CLOSED)
    for j, term in enumerate(TERMS):
        ref("veteran", M, f"{term}: coefficient", r["coef"][j], "R coxph", ITERATIVE)
        ref("veteran", M, f"{term}: standard error", r["se"][j], "R coxph", ITERATIVE)

# ---- start-stop rows (Stanford heart) and a Weibull model (larynx)
heart = load_csv("heart.csv")
checks.data("heart", heart)
TERMS = ["age", "surgery", "transplant"]
M = "explain_survival (followed_from)"
ctv = CoxTimeVaryingFitter().fit(heart[["id", "start", "stop", "event", "age", "surgery", "transplant"]].astype(float),
                                 id_col="id", event_col="event", start_col="start", stop_col="stop")
for term in TERMS:
    ref("Stanford heart", M, f"{term}: coefficient", ctv.params_[term], "lifelines CoxTimeVaryingFitter", ITERATIVE)
    ref("Stanford heart", M, f"{term}: standard error", ctv.standard_errors_[term], "lifelines CoxTimeVaryingFitter",
        ITERATIVE)
larynx = load_csv("larynx.csv")
larynx["stage"] = "stage " + larynx["stage"].astype(str)
checks.data("larynx", larynx)
WTERMS = ["stage: stage 2 vs stage 1", "stage: stage 3 vs stage 1", "stage: stage 4 vs stage 1", "age"]
W = "explain_survival (distribution weibull)"
design = pd.get_dummies(larynx["stage"], drop_first=True, dtype=float).assign(age=larynx["age"], T=larynx["time"],
                                                                             E=larynx["delta"])
aft = WeibullAFTFitter()
aft._scipy_fit_method = "L-BFGS-B"
aft.fit(design, "T", "E", fit_options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 10000})
for term, column in zip(WTERMS, ["stage 2", "stage 3", "stage 4", "age"]):
    ref("larynx (KMsurv)", W, f"{term}: coefficient", aft.params_[("lambda_", column)], "lifelines WeibullAFTFitter",
        ITERATIVE)
if HAS_R:
    r = r_run("""m <- survival::coxph(survival::Surv(start, stop, event) ~ age + surgery + transplant)
list(coef = unname(coef(m)), se = unname(sqrt(diag(vcov(m)))))""", heart)
    for j, term in enumerate(TERMS):
        ref("Stanford heart", M, f"{term}: coefficient", r["coef"][j], "R coxph", ITERATIVE)
        ref("Stanford heart", M, f"{term}: standard error", r["se"][j], "R coxph", ITERATIVE)
    r = r_run("""m <- survival::survreg(survival::Surv(time, delta) ~ factor(stage) + age, dist = "weibull")
list(coef = unname(coef(m))[-1], se = unname(sqrt(diag(vcov(m))))[2:5])""", larynx)
    for j, term in enumerate(WTERMS):
        ref("larynx (KMsurv)", W, f"{term}: coefficient", r["coef"][j], "R survreg", ITERATIVE)
        ref("larynx (KMsurv)", W, f"{term}: standard error", r["se"][j], "R survreg", ITERATIVE)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
