"""Reference values for notebooks/30_survival_more.ipynb.

No printed numbers match these checks, so the help-page or vignette example of the R function simplydata calls is run
once (survival survfit, coxph, survreg; cmprsk cuminc and crr; survRM2 rmst2; adjustedCurves; mice) and saved with the
link to its page. Also saves the data the notebook reads: mgus2 with the competing-risks set-up of the survival
package's vignette "compete", the dates made from its times, and ovarian.
"""

import numpy as np
import pandas as pd

from _tools import References, load_csv, r_available, r_dataset, r_run, rman

refs = References("30_survival_more")
HAS_R = r_available()
CLOSED, ITERATIVE = 1e-10, 1e-6

mgus = r_dataset("mgus2", "survival")
mgus["etime"] = np.where(mgus["pstat"] == 0, mgus["futime"], mgus["ptime"])
mgus["event"] = np.where(mgus["pstat"] == 1, "pcm", np.where(mgus["death"] == 1, "death", "censor"))
refs.data("mgus2", mgus)
R_MGUS = 'ev <- factor(event, levels = c("censor", "pcm", "death"))\n'

# 1. dates made from the times (follow-up from 2000-01-01, months of 30.4375 days), and the times they give back
d = mgus[["id", "etime", "event"]].copy()
start = pd.Timestamp("2000-01-01")
days = np.round(d["etime"] * 30.4375)
d["start"] = start.strftime("%Y-%m-%d")
d["end"] = (start + pd.to_timedelta(days, unit="D")).dt.strftime("%Y-%m-%d")
d["pcm_date"] = d["end"].where(d["event"] == "pcm")
d["death_date"] = d["end"].where(d["event"] == "death")
refs.data("mgus2_dates", d[["start", "end", "pcm_date", "death_date"]])
refs.run_vector("mgus2 (survival), dates made from the times", 'add_time_and_event_columns (competing_with=["death_date"], unit="months")',
                "time (months)", days / 30.4375, "the vignette's etime, in whole days", rman("survival", "mgus2"), "survival", CLOSED)
refs.run_vector("mgus2 (survival), dates made from the times", 'add_time_and_event_columns (competing_with=["death_date"], unit="months")',
                "event", d["event"].map({"pcm": "pcm_date", "death": "death_date", "censor": "still followed"}), "the vignette's event",
                rman("survival", "mgus2"), "survival", 0, "text")

ovarian = r_dataset("ovarian", "survival")
refs.data("ovarian", ovarian)
lung = load_csv("lung.csv")
lung_cc = lung.dropna(subset=["age", "sex", "ph.ecog"]).reset_index(drop=True)
refs.data("lung_complete", lung_cc)
late = lung[lung["time"] > 30].assign(entry=30.0)
refs.data("lung_after_day_30", late)
lung_a = lung.dropna(subset=["age"]).assign(sex_group=lambda t: t["sex"].map({1: "male", 2: "female"}))
refs.data("lung_sex_group", lung_a)

if HAS_R:
    POINTS = [60, 120, 240]
    r = r_run(R_MGUS + '''f <- survival::survfit(survival::Surv(etime, ev) ~ 1); s <- summary(f, times = c(60, 120, 240)); k <- match(c("pcm", "death"), f$states)
                 list(pcm = s$pstate[, k[1]], death = s$pstate[, k[2]], lo = s$lower[, k[1]], hi = s$upper[, k[1]])''', mgus, packages=("survival",))
    check = 'find_risk_over_time (competing_with="death", at_time_points=[60, 120, 240])'
    for i, p in enumerate(POINTS):
        refs.run_all("mgus2 (survival)", check, {f"risk of pcm by {p} months": r["pcm"][i], f"risk of death by {p} months": r["death"][i],
                                                 f"pcm by {p} months: 95 % low": r["lo"][i], f"pcm by {p} months: 95 % high": r["hi"][i]},
                     "R survival survfit (Aalen-Johansen; vignette compete)", rman("survival", "survfit.formula"), "survival", CLOSED)
    r = r_run(R_MGUS + "c1 <- cmprsk::cuminc(etime, as.integer(ev) - 1, group = sex, cencode = 0); list(s = c1$Tests[1, 1], p = c1$Tests[1, 2])", mgus,
              packages=("cmprsk",))
    refs.run_all("mgus2 (survival)", 'compare_survival (competing_with="death", for_each_group_in="sex")', {"Gray's chi-square": r["s"], "p-value": r["p"]},
                 "R cmprsk cuminc", rman("cmprsk", "cuminc"), "cmprsk", CLOSED)
    for test_method, code, label, url in (
            ("cox", 'm <- survival::coxph(survival::Surv(etime, ev == "pcm") ~ age + sex); b <- coef(m); s <- sqrt(diag(vcov(m)))',
             "R survival coxph (death as censoring)", rman("survival", "coxph")),
            ("fine-gray", 'm <- cmprsk::crr(etime, as.integer(ev) - 1, cbind(age, sex == "M"), failcode = 1, cencode = 0); b <- m$coef; s <- sqrt(diag(m$var))',
             "R cmprsk crr", rman("cmprsk", "crr"))):
        r = r_run(R_MGUS + code + "\nlist(a = b[[1]], b = b[[2]], c = s[[1]], d = s[[2]])", mgus, packages=("survival", "cmprsk"))
        refs.run_all("mgus2 (survival)", f'explain_survival (competing_with="death", test_method="{test_method}")',
                     {"age coefficient": r["a"], "sex coefficient": r["b"], "age standard error": r["c"], "sex standard error": r["d"]},
                     label, url, "survival" if test_method == "cox" else "cmprsk", ITERATIVE)
    r = r_run('''m <- survival::coxph(survival::Surv(time, status == 2) ~ age + sex + ph.ecog); list(b = unname(coef(m)), s = unname(sqrt(diag(vcov(m)))))''',
              lung_cc, packages=("survival",))
    refs.run_all("lung (survival)", 'explain_survival (test_method="cox")',
                 {**{f"{t} coefficient": v for t, v in zip(["age", "sex", "ph.ecog"], r["b"])},
                  **{f"{t} standard error": v for t, v in zip(["age", "sex", "ph.ecog"], r["s"])}},
                 "R survival coxph (help page lung)", rman("survival", "coxph"), "survival", ITERATIVE)
    for distribution, dist in (("weibull", "weibull"), ("log-normal", "lognormal")):
        r = r_run(f'''s <- summary(survival::survreg(survival::Surv(futime, fustat) ~ ecog.ps + rx, dist = "{dist}"))$table
                      list(a = s[2, 1], b = s[3, 1], c = s[2, 2], d = s[3, 2])''', ovarian, packages=("survival",))
        refs.run_all("ovarian (survival)", f'explain_survival (distribution="{distribution}")',
                     {"ecog.ps coefficient": r["a"], "rx coefficient": r["b"], "ecog.ps standard error": r["c"], "rx standard error": r["d"]},
                     "R survival survreg (help page ovarian)", rman("survival", "survreg"), "survival", ITERATIVE)
    r = r_run('''r <- survRM2::rmst2(time, as.integer(status == 2), arm = as.integer(sex == 2), tau = 365)
                 list(a0 = r$RMST.arm0$rmst[1], a1 = r$RMST.arm1$rmst[1], diff = r$unadjusted.result[1, 1], ratio = r$unadjusted.result[2, 1])''', lung,
              packages=("survRM2",))
    refs.run_all("lung (survival)", "compare_survival (average_time_up_to=365)",
                 {"average time up to 365 days, sex 1": r["a0"], "average time up to 365 days, sex 2": r["a1"], "difference (2 minus 1)": r["diff"],
                  "ratio (2 over 1)": r["ratio"]}, "R survRM2 rmst2", rman("survRM2", "rmst2"), "survRM2", CLOSED)
    r = r_run('''d <- data.frame(time = time, status = as.integer(status == 2), age = age, sex_group = factor(sex_group))
                 m <- survival::coxph(survival::Surv(time, status) ~ sex_group + age, data = d, x = TRUE)
                 a <- adjustedCurves::adjustedsurv(data = d, variable = "sex_group", ev_time = "time", event = "status", method = "direct",
                                                   outcome_model = m, times = c(180, 365))$adj
                 list(group = as.character(a$group), time = a$time, surv = a$surv)''', lung_a, packages=("survival", "adjustedCurves"))
    refs.run_all("lung (survival)", 'find_survival_over_time (allowing_for=["age"], at_time_points=[180, 365])',
                 {f"{g}, day {t:g}": s for g, t, s in zip(r["group"], r["time"], r["surv"])}, "R adjustedCurves adjustedsurv (direct)",
                 rman("adjustedCurves", "adjustedsurv"), "adjustedCurves", ITERATIVE)
    r = r_run('''d <- data.frame(time, status = as.integer(status == 2), age, ph.ecog, wt.loss); d$H <- mice::nelsonaalen(d, time, status)
                 imp <- mice::mice(d, m = 400, method = "pmm", seed = 1, printFlag = FALSE)
                 p <- mice::pool(with(imp, survival::coxph(survival::Surv(time, status) ~ age + ph.ecog + wt.loss)))$pooled
                 list(b = p$estimate, B = p$b)''', lung, packages=("mice", "survival"))
    for i, c in enumerate(["age", "ph.ecog", "wt.loss"]):
        refs.run("lung (survival)", 'explain_survival (empty_cells="multiple imputation", with_seed=)', f"{c} pooled coefficient", r["b"][i],
                 "R mice (m = 400) + pool(coxph), within 4 MCSE", rman("mice", "pool"), "mice", 4 * np.sqrt(r["B"][i] / 20 + r["B"][i] / 400), "absolute")

# 7. the numbers the survival charts draw (chart.get_numbers_drawn()): each curve starts at time 0 (survival 1, risk 0),
# steps at every time with an event and runs on to the last time followed, as R plot(survfit) draws it
STEPS = """k <- f$n.event > 0; if (is.matrix(k)) k <- rowSums(f$n.event) > 0
           list(time = c(0, f$time[k], max(f$time)), y = c(Y0, Y[k], Y[length(Y)]))"""
if HAS_R:
    for sex in (1, 2):
        r = r_run(f"f <- survival::survfit(survival::Surv(time[sex == {sex}], status[sex == {sex}] == 2) ~ 1); Y <- f$surv; Y0 <- 1\n" + STEPS,
                  lung, packages=("survival",))
        for key, name in (("time", "time"), ("y", "survival")):
            refs.run_vector("lung (survival)", 'chart draw_survival_curves (colour="sex")', f"curve of sex {sex}: {name} at each step", r[key],
                            "R survival survfit (Kaplan-Meier)", rman("survival", "survfit.formula"), "survival", CLOSED)
    r = r_run("f <- survival::survfit(survival::Surv(entry, time, status == 2) ~ 1); Y <- f$surv; Y0 <- 1\n" + STEPS, late, packages=("survival",))
    for key, name in (("time", "time"), ("y", "survival")):
        refs.run_vector("lung after day 30 (survival)", 'chart draw_survival_curves (followed_from="entry")', f"curve: {name} at each step", r[key],
                        "R survival survfit (Surv(entry, time, event), late entry)", rman("survival", "survfit.formula"), "survival", CLOSED)
    r = r_run('''d <- data.frame(time = time, status = as.integer(status == 2), age = age, sex_group = factor(sex_group))
                 m <- survival::coxph(survival::Surv(time, status) ~ sex_group + age, data = d, x = TRUE)
                 a <- adjustedCurves::adjustedsurv(data = d, variable = "sex_group", ev_time = "time", event = "status", method = "direct",
                                                   outcome_model = m)$adj
                 list(group = as.character(a$group), time = a$time, surv = a$surv)''', lung_a, packages=("survival", "adjustedCurves"))
    groups = np.asarray(r["group"])
    for g in ("female", "male"):
        for key, name in (("time", "time"), ("surv", "survival")):
            refs.run_vector("lung (survival)", 'chart draw_survival_curves (allowing_for=["age"])', f"adjusted curve of {g}: {name} at each step",
                            np.asarray(r[key])[groups == g], "R adjustedCurves adjustedsurv (direct, Cox model)", rman("adjustedCurves", "adjustedsurv"),
                            "adjustedCurves", ITERATIVE)
    r = r_run(R_MGUS + "f <- survival::survfit(survival::Surv(etime, ev) ~ 1); Y <- f$pstate[, match('pcm', f$states)]; Y0 <- 0\n" + STEPS, mgus,
              packages=("survival",))
    for key, name in (("time", "time"), ("y", "risk of pcm")):
        refs.run_vector("mgus2 (survival)", 'chart draw_risk_curves (competing_with="death")', f"curve: {name} at each step", r[key],
                        "R survival survfit (Aalen-Johansen; vignette compete)", rman("survival", "survfit.formula"), "survival", CLOSED)

refs.save()
