"""Reference values for notebooks/26_models_mixed_and_more.ipynb.

Printed first: Bates, Maechler, Bolker & Walker (2015) for sleepstudy (the lme4 vignette "lmer"); the emmeans vignette
"basics" for the pigs and mtcars marginal means. Where nothing is printed, the help-page example of the library
simplydata calls is run once (lme4, lmerTest, geepack, nlme, mgcv, splines, mice, SciPy bootstrap, emmeans, anova,
AICcmodavg, boot, glmnet) and saved with the link to its page. Also saves the data the notebook reads.
"""

import numpy as np
import pandas as pd
from scipy import stats

from _tools import References, VIGNETTE, load_csv, r_available, r_dataset, r_run, r_version, rman

refs = References("26_models_mixed_and_more")
HAS_R = r_available()
CLOSED, ITERATIVE, SAME_SEED = 1e-10, 1e-6, 1e-12

# 1-2. sleepstudy
sleep = r_dataset("sleepstudy", "lme4")
sleep["Days2"] = sleep["Days"] ** 2
refs.data("sleepstudy", sleep)
B15 = "Bates, Maechler, Bolker & Walker (2015), J. Stat. Softw. 67(1), section 3.6 and 5.2.4"
B15_URL = "https://doi.org/10.18637/jss.v067.i01 ; " + VIGNETTE.format(package="lme4", file="lmer.pdf")
lme4 = ("lme4", r_version("lme4"))
refs.printed_all("sleepstudy (lme4)", 'explain (grouped_by=, own_slope_for="Days")',
                 {"(intercept) estimate": "251.41", "(intercept) standard error": "6.82", "Days estimate": "10.47", "Days standard error": "1.55"},
                 B15, B15_URL, "published", *lme4)
refs.printed_all("sleepstudy (lme4)", "compare_models (grouped_by=)",
                 {"AIC of Reaction ~ Days + (1 | Subject), ML": "1802", "BIC of Reaction ~ Days + (1 | Subject), ML": "1815",
                  "log-likelihood of Reaction ~ Days + (1 | Subject), ML": "-897"}, B15 + " (anova of fm3)", B15_URL, "published", *lme4)
if HAS_R:
    r = r_run('''d <- data.frame(Reaction = Reaction, Days = Days, Days2 = Days2, Subject = factor(Subject))
                 m1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = d)
                 s <- summary(lmerTest::as_lmerModLmerTest(m1))$coefficients
                 m0 <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = d); m2 <- lme4::lmer(Reaction ~ Days + Days2 + (1 | Subject), data = d)
                 a <- suppressMessages(anova(m0, m2))
                 mm <- lme4::lmer(Reaction ~ Days + (Days | Subject), data = d)
                 re <- lme4::ranef(mm, condVar = TRUE)$Subject; pv <- attr(re, "postVar"); cf <- coef(mm)$Subject
                 list(b1 = unname(lme4::fixef(m1)), se1 = unname(sqrt(diag(as.matrix(vcov(m1))))), df = unname(s[, "df"]),
                      ps = unname(s[, "Pr(>|t|)"]), aic = a$AIC, bic = a$BIC, ll = a$logLik, chi = a$Chisq[2], pchi = a[["Pr(>Chisq)"]][2],
                      b0 = re[, 1], bs = re[, 2], sd0 = sqrt(pv[1, 1, ]), sd1 = sqrt(pv[2, 2, ]), own0 = cf[, 1], own1 = cf[, 2],
                      fit1 = unname(fitted(m1)))''', sleep, packages=("lme4", "lmerTest"))
    url = rman("lme4", "lmer")
    for i, term in enumerate(["(intercept)", "Days"]):
        refs.run_all("sleepstudy (lme4)", "explain (grouped_by=)", {f"{term} estimate": r["b1"][i], f"{term} standard error": r["se1"][i]},
                     "R lme4 lmer (help page sleepstudy)", url, "lme4", ITERATIVE)
        for label in ('explain (grouped_by=, p_value_by="satterthwaite")', 'explain (grouped_by=, p_value_by="few groups")'):
            refs.run_all("sleepstudy (lme4)", label, {f"{term} degrees of freedom": r["df"][i], f"{term} p-value": r["ps"][i]},
                         "R lmerTest (Satterthwaite)", rman("lmerTest", "summary.lmerModLmerTest"), "lmerTest", ITERATIVE)
    refs.run_all("sleepstudy (lme4)", "compare_models (grouped_by=)",
                 {"AIC bigger": r["aic"][1], "BIC bigger": r["bic"][1], "log-likelihood bigger": r["ll"][1], "chi-square": r["chi"], "p-value": r["pchi"]},
                 "R anova of two lmer fits (refitted by ML)", rman("lme4", "anova.merMod"), "lme4", ITERATIVE)
    check = 'find_group_levels_of (own_slope_for="Days")'
    for name, key in (("starting level", "b0"), ("slope", "bs"), ("starting level conditional SD", "sd0"), ("slope conditional SD", "sd1"),
                      ("own starting level", "own0"), ("own slope", "own1")):
        refs.run_vector("sleepstudy (lme4)", check, name, r[key], "R lme4 ranef(condVar = TRUE) and coef", rman("lme4", "ranef"), "lme4", ITERATIVE)
    refs.run_vector("sleepstudy (lme4)", "add_predicted_column (grouped_by=)", "predicted Reaction", r["fit1"], "R fitted(lmer)",
                    rman("lme4", "lmer"), "lme4", ITERATIVE)

# 3. cbpp (help page of glmer)
cbpp = r_dataset("cbpp", "lme4")
cbpp["period"] = cbpp["period"].astype(str)
refs.data("cbpp", cbpp)
refs.data("cbpp_one_row_per_animal", pd.DataFrame([(h, p, int(k < i)) for h, p, i, n in cbpp[["herd", "period", "incidence", "size"]].itertuples(index=False)
                                                   for k in range(int(n))], columns=["herd", "period", "sick"]))
if HAS_R:
    r = r_run('''d <- data.frame(incidence = incidence, size = size, period = factor(period), herd = factor(herd))
                 s <- summary(lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd), family = binomial, data = d))$coefficients
                 list(b = unname(s[, 1]), se = unname(s[, 2]), p = unname(s[, 4]))''', cbpp, packages=("lme4",))
    for i, term in enumerate(["(intercept)", "period (2)", "period (3)", "period (4)"]):
        refs.run_all("cbpp (lme4)", 'explain (outcome="yes or no", out_of=, grouped_by=)',
                     {f"{term} log odds ratio": r["b"][i], f"{term} standard error": r["se"][i], f"{term} p-value": r["p"][i]},
                     "R lme4 glmer (help page example)", rman("lme4", "glmer"), "lme4", ITERATIVE)

# 4. dietox (geepack)
diet = r_dataset("dietox", "geepack").sort_values(["Pig", "Time"]).reset_index(drop=True)
diet["Cu"] = diet["Cu"].astype(str)
refs.data("dietox", diet)
if HAS_R:
    for within, corstr in (("alike", "exchangeable"), ("in time order", "ar1"), ("independent", "independence")):
        r = r_run(f'''d <- data.frame(Weight = Weight, Time = Time, Cu = factor(Cu), Pig = Pig)
                      s <- coef(summary(geepack::geeglm(Weight ~ Time + Cu, id = Pig, data = d, corstr = "{corstr}")))
                      list(b = unname(s[, 1]), se = unname(s[, 2]))''', diet, packages=("geepack",))
        terms = ["(intercept)", "Time", "Cu (Cu035)", "Cu (Cu175)"]
        refs.run_all("dietox (geepack)", f'explain (effects="for the whole population", within_groups="{within}")',
                     {**{f"{t} estimate": v for t, v in zip(terms, r["b"])}, **{f"{t} standard error": v for t, v in zip(terms, r["se"])}},
                     "R geepack geeglm", rman("geepack", "geeglm"), "geepack", ITERATIVE)

# 5. Ovary (nlme gls)
ovary = r_dataset("Ovary", "nlme")
ovary["s"] = np.sin(2 * np.pi * ovary["Time"])
ovary["c"] = np.cos(2 * np.pi * ovary["Time"])
ovary["order"] = np.arange(1, len(ovary) + 1)
ovary["Mare"] = "mare " + ovary["Mare"].astype(str)
refs.data("Ovary", ovary)
if HAS_R:
    for option, value, part in (("neighbours_alike_in_order_of", "order", "correlation = nlme::corAR1(form = ~ order)"),
                                ("spread_changes_with", "Mare", "weights = nlme::varIdent(form = ~ 1 | Mare)")):
        r = r_run(f'''d <- data.frame(follicles = follicles, s = s, c = c, order = order, Mare = factor(Mare))
                      m <- nlme::gls(follicles ~ s + c, data = d, {part}); t <- summary(m)$tTable
                      list(b = t[, 1], se = t[, 2], p = t[, 4], rho = if ("{option}" == "neighbours_alike_in_order_of")
                           coef(m$modelStruct$corStruct, unconstrained = FALSE)[[1]] else NA)''', ovary, packages=("nlme",))
        values = {f"{t} {c}": r[k][i] for i, t in enumerate(["(intercept)", "s", "c"]) for c, k in (("estimate", "b"), ("standard error", "se"), ("p-value", "p"))}
        if option == "neighbours_alike_in_order_of":
            values["correlation of neighbours"] = r["rho"]
        refs.run_all("Ovary (nlme)", f'explain ({option}="{value}")', values, "R nlme gls", rman("nlme", "gls"), "nlme", ITERATIVE)

# 6. GAM (mcycle) and natural spline (cars)
mcycle = r_dataset("mcycle", "MASS")
refs.data("mcycle", mcycle)
cars = r_dataset("cars")
refs.data("cars", cars)
if HAS_R:
    r = r_run('''s <- summary(mgcv::gam(accel ~ s(times), method = "REML"))
                 list(intercept = s$p.table[1, 1], `intercept standard error` = s$p.table[1, 2],
                      `effective degrees of freedom of s(times)` = s$s.table[1, 1], `F of s(times)` = s$s.table[1, 3])''', mcycle, packages=("mgcv",))
    refs.run_all("mcycle (MASS)", 'explain (curve_for="times")', r, "R mgcv gam (REML)", rman("mgcv", "gam"), "mgcv", ITERATIVE)
    r = r_run('''a <- anova(lm(dist ~ 1), lm(dist ~ splines::ns(speed, df = 4)))
                 list(`F of the whole curve` = a$F[2], `p-value of the whole curve` = a[["Pr(>F)"]][2])''', cars, packages=("splines",))
    refs.run_all("cars", 'explain (curve_for="speed", with_knots=3)', r, "R lm with splines::ns(df = 4)", rman("splines", "ns"), "splines", CLOSED)

# 7. multiple imputation (nhanes): R mice, m = 400; tolerance 4 Monte Carlo standard errors (R's between-imputation
#    variance for simplydata's 20 imputations and R's 400)
nhanes = r_dataset("nhanes", "mice")
refs.data("nhanes", nhanes)
if HAS_R:
    r = r_run('''imp <- mice::mice(data.frame(age = age, bmi = bmi, hyp = hyp, chl = chl), m = 400, method = "pmm", seed = 1, printFlag = FALSE)
                 p <- mice::pool(with(imp, lm(chl ~ age + bmi)))$pooled
                 list(b = p$estimate, B = p$b)''', nhanes, packages=("mice",))
    for i, t in enumerate(["(intercept)", "age", "bmi"]):
        refs.run("nhanes (mice)", 'explain (empty_cells="multiple imputation", with_seed=)', f"{t} pooled estimate", r["b"][i],
                 "R mice (m = 400) + pool, within 4 MCSE", rman("mice", "mice"), "mice", 4 * np.sqrt(r["B"][i] / 20 + r["B"][i] / 400), "absolute")

# 8. bootstrap: SciPy stats.bootstrap (BCa), the call simplydata makes, with the same seed
stackloss = load_csv("stackloss.csv")
x, y = stackloss["air_flow"].to_numpy(), stackloss["stack_loss"].to_numpy()


def slope_and_intercept(xs, ys, axis=-1):
    xm, ym = xs.mean(axis=axis, keepdims=True), ys.mean(axis=axis, keepdims=True)
    b = ((xs - xm) * (ys - ym)).sum(axis=axis) / ((xs - xm) ** 2).sum(axis=axis)
    return np.stack([np.squeeze(ym, axis) - b * np.squeeze(xm, axis), b])


res = stats.bootstrap((x, y), slope_and_intercept, n_resamples=2000, paired=True, vectorized=True, method="BCa",
                      random_state=np.random.default_rng(7))
for i, t in enumerate(["(intercept)", "air_flow"]):
    refs.run_all("stackloss", 'explain (interval_by="bootstrap", resamples=, with_seed=)',
                 {f"{t} 95 % low": res.confidence_interval.low[i], f"{t} 95 % high": res.confidence_interval.high[i]},
                 "SciPy stats.bootstrap (BCa, same seed)", "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.bootstrap.html",
                 "scipy", SAME_SEED)

# 9. model averages: the emmeans vignette "basics" (printed); the other cases run emmeans once
EMB = "emmeans vignette 'basics' (Lenth), printed output"
EMB_URL = VIGNETTE.format(package="emmeans", file="basics.html")
emv = ("emmeans", r_version("emmeans"))
for g, (m, s) in {"fish": ("0.0337", "0.000926"), "soy": ("0.0257", "0.000945"), "skim": ("0.0229", "0.000994")}.items():
    refs.printed_all("pigs (emmeans)", "find_model_averages_of",
                     {f"{g} model average": m, f"{g} standard error": s, f"{g} degrees of freedom": "23"}, EMB, EMB_URL, "printed", *emv)
# R prints this column with 3 significant digits padded to 6 decimals ("0.001030" is 0.00103)
for g, (m, s) in {"09": ("0.0322", "0.00103"), "12": ("0.0270", "0.000969"), "15": ("0.0263", "0.00110"), "18": ("0.0241", "0.00134")}.items():
    refs.printed_all("pigs (emmeans)", "find_model_averages_of", {f"percent {g} model average": m, f"percent {g} standard error": s},
                     EMB, EMB_URL, "printed", *emv)
for g, (m, s) in {"4": ("20.8", "2.05"), "6": ("18.7", "1.19"), "8": ("20.2", "1.77")}.items():
    refs.printed_all("mtcars", "find_model_averages_of", {f"{g} cylinders model average": m, f"{g} cylinders standard error": s},
                     EMB, EMB_URL, "printed", *emv)
pigs = r_dataset("pigs", "emmeans")
pigs["inverse_conc"] = 1 / pigs["conc"]
pigs["percent_group"] = "p" + pigs["percent"].astype(int).astype(str).str.zfill(2)
refs.data("pigs", pigs)
mt = load_csv("mtcars.csv")
mt["cylinders"] = "c" + mt["cyl"].astype(int).astype(str)
mt["dispsq"] = mt["disp"] ** 2
refs.data("mtcars", mt)
warp = load_csv("warpbreaks.csv")
R_WARP = "d <- data.frame(breaks = breaks, wool = factor(wool), tension = factor(tension, levels = c('L', 'M', 'H')))\n"
EMM = rman("emmeans", "emmeans")
if HAS_R:
    r = r_run(R_WARP + '''m <- lm(breaks ~ wool * tension, data = d); e <- emmeans::emmeans(m, ~ tension); s <- summary(e)
                 pp <- summary(pairs(e, adjust = "tukey"))
                 list(mean = s$emmean, se = s$SE, est = pp$estimate, p = pp$p.value, contrast = as.character(pp$contrast))''', warp, packages=("emmeans",))
    values = {f"{t} model average": v for t, v in zip(["L", "M", "H"], r["mean"])}
    values.update({f"{t} standard error": v for t, v in zip(["L", "M", "H"], r["se"])})
    for label, est, p in zip(r["contrast"], r["est"], r["p"]):
        a, b = [s.strip() for s in label.split(" - ")]
        key, sign = (f"{a} - {b}", 1) if a < b else (f"{b} - {a}", -1)   # simplydata names a pair in sorted order
        values[f"{key} difference"] = sign * est
        values[f"{key} p-value (Tukey)"] = p
    refs.run_all("warpbreaks", "find_model_averages_of (wool * tension)", values, "R emmeans + pairs (help page example)", EMM, "emmeans", CLOSED)
    for comparing, how in {"each with the first": 'method = "trt.vs.ctrl"', "each with the next": 'method = "consec"',
                           "each with the average": 'method = "eff"', "trend": 'method = "poly"',
                           "each with the average of the ones before": 'method = list("M vs L" = c(-1, 1, 0), "H vs L and M" = c(-1/2, -1/2, 1))'}.items():
        r = r_run(R_WARP + f'''c <- summary(emmeans::contrast(emmeans::emmeans(lm(breaks ~ wool + tension, data = d), ~ tension), {how}))
                     list(est = c$estimate, se = c$SE, p = c$p.value)''', warp, packages=("emmeans",))
        values = {}
        for i in range(len(r["est"])):
            values[f"contrast {i + 1} estimate"] = r["est"][i]
            values[f"contrast {i + 1} standard error"] = r["se"][i]
            if comparing in ("trend", "each with the average of the ones before"):
                values[f"contrast {i + 1} p-value"] = r["p"][i]
        refs.run_all("warpbreaks", f'find_model_averages_of (comparing="{comparing}")', values, "R emmeans contrast",
                     rman("emmeans", "contrast"), "emmeans", CLOSED)
    r = r_run(R_WARP + '''s <- summary(emmeans::emmeans(lm(breaks ~ wool * tension, data = d), ~ tension | wool))
                 list(mean = s$emmean, se = s$SE, wool = as.character(s$wool), tension = as.character(s$tension))''', warp, packages=("emmeans",))
    values = {}
    for m, s, w, t in zip(r["mean"], r["se"], r["wool"], r["tension"]):
        values[f"{t} (within {w}) model average"] = m
        values[f"{t} (within {w}) standard error"] = s
    refs.run_all("warpbreaks", 'find_model_averages_of (within_each_group_in="wool")', values, "R emmeans (~ tension | wool)", EMM, "emmeans", CLOSED)
    r = r_run('''d <- data.frame(mpg = mpg, wt = wt, cylinders = factor(cylinders)); m <- lm(mpg ~ cylinders + wt, data = d)
                 s <- summary(emmeans::emmeans(m, ~ wt | cylinders, at = list(wt = seq(min(d$wt), max(d$wt), length.out = 50))))
                 list(mean = s$emmean, se = s$SE)''', mt, packages=("emmeans",))
    refs.run_vector("mtcars", 'find_model_averages_of (over_range_of="wt")', "model average", r["mean"], "R emmeans (at = 50 values)", EMM, "emmeans", CLOSED)
    refs.run_vector("mtcars", 'find_model_averages_of (over_range_of="wt")', "standard error", r["se"], "R emmeans (at = 50 values)", EMM, "emmeans", CLOSED)

birthwt = load_csv("birthwt.csv")
birthwt["race_group"] = birthwt["race"].map({1: "white", 2: "black", 3: "other"})
refs.data("birthwt", birthwt)
nmes = load_csv("nmes1988.csv")
gas = r_dataset("GasolineYield", "betareg")
gas["batch"] = "b" + gas["batch"].astype(str).str.zfill(2)
refs.data("GasolineYield", gas)
housing = r_dataset("housing", "MASS")
expanded = housing.loc[housing.index.repeat(housing["Freq"].astype(int))].reset_index(drop=True)
refs.data("housing_one_row_per_tenant", expanded)
cbpp_rows = pd.read_csv(refs.data_dir / "cbpp_one_row_per_animal.csv")
CASES = [
    ("birthwt", birthwt, 'outcome="yes or no"', "d <- data.frame(low = low, race_group = factor(race_group), lwt = lwt)\nm <- glm(low ~ race_group + lwt, family = binomial, data = d)",
     "race_group", ["black", "other", "white"]),
    ("NMES 1988", nmes, 'outcome="counts"', "d <- data.frame(visits = visits, health = factor(health), chronic = chronic)\nm <- glm(visits ~ health + chronic, family = poisson, data = d)",
     "health", ["average", "excellent", "poor"]),
    ("NMES 1988", nmes, 'outcome="counts (overdispersed)"', "d <- data.frame(visits = visits, health = factor(health), chronic = chronic)\nm <- glm(visits ~ health + chronic, family = quasipoisson, data = d)",
     "health", ["average", "excellent", "poor"]),
    ("NMES 1988", nmes, 'outcome="counts that vary a lot"', "d <- data.frame(visits = visits, health = factor(health), chronic = chronic)\nm <- MASS::glm.nb(visits ~ health + chronic, data = d)",
     "health", ["average", "excellent", "poor"]),
    ("GasolineYield (betareg)", gas, 'outcome="shares"', "d <- data.frame(yield = yield, batch = factor(batch), temp = temp)\nm <- betareg::betareg(yield ~ batch + temp, data = d)",
     "batch", sorted(gas["batch"].unique())),
    ("cbpp (lme4), one row per animal", cbpp_rows, 'outcome="yes or no", grouped_by="herd"',
     "d <- data.frame(sick = sick, period = factor(period), herd = factor(herd))\nm <- lme4::glmer(sick ~ period + (1 | herd), family = binomial, data = d)",
     "period", ["1", "2", "3", "4"]),
]
if HAS_R:
    for dataset, data, options, model, groups, levels in CASES:
        r = r_run(model + f'''\ns <- summary(emmeans::emmeans(m, ~ {groups}, type = "response")); list(mean = s[[2]], se = s$SE)''', data,
                  packages=("emmeans", "MASS", "betareg", "lme4"))
        values = {}
        for g, m, s in zip(levels, r["mean"], r["se"]):
            values[f"{g} model average"] = m
            values[f"{g} standard error"] = s
        refs.run_all(dataset, f"find_model_averages_of ({options})", values, "R emmeans (type = response)", EMM, "emmeans", ITERATIVE)
    for outcome, model in (("ordered groups", "m <- MASS::polr(Sat ~ Infl + Type + Cont, data = d, Hess = TRUE)"),
                           ("groups", "m <- nnet::multinom(Sat ~ Infl + Type + Cont, data = d, trace = FALSE)")):
        levels = 'levels = c("Low", "Medium", "High")' + (", ordered = TRUE" if outcome == "ordered groups" else "")
        r = r_run(f'''d <- data.frame(Sat = factor(Sat, {levels}), Infl = factor(Infl), Type = factor(Type), Cont = factor(Cont))
                      {model}
                      s <- summary(emmeans::emmeans(m, ~ Sat | Infl, mode = "prob"))
                      list(p = s$prob, se = s$SE, sat = as.character(s$Sat), infl = as.character(s$Infl))''', expanded, packages=("emmeans", "MASS", "nnet"))
        values = {}
        for p, se, sat, infl in zip(r["p"], r["se"], r["sat"], r["infl"]):
            values[f"chance of {sat} for Infl {infl}"] = p
            values[f"standard error of the chance of {sat} for Infl {infl}"] = se
        refs.run_all("housing (MASS)", f'find_model_averages_of (outcome="{outcome}")', values, "R emmeans (mode = prob)", EMM, "emmeans", ITERATIVE)

# 10. comparing models
if HAS_R:
    r = r_run('''m0 <- lm(stack_loss ~ air_flow); m1 <- lm(stack_loss ~ air_flow + water_temp + acid_conc); a <- anova(m0, m1)
                 list(`AIC smaller` = AIC(m0), `AIC bigger` = AIC(m1), statistic = a$F[2], `p-value` = a[["Pr(>F)"]][2])''', stackloss)
    refs.run_all("stackloss", 'compare_models (outcome="numbers")', r, "R anova(m0, m1) and AIC", rman("stats", "anova.lm"), "stats", CLOSED)
bio = r_dataset("bioChemists", "pscl")
refs.data("bioChemists", bio)
birthwt["bwt_kg"] = birthwt["bwt"] / 1000
refs.data("birthwt", birthwt)
GLMS = [("NMES 1988", nmes, "visits", ["hospital"], ["hospital", "chronic"], "counts", "glm(visits ~ {f}, family = poisson, control = glm.control(epsilon = 1e-12, maxit = 100))", "Chisq", "glm"),
        ("NMES 1988", nmes, "visits", ["hospital"], ["hospital", "chronic"], "counts (overdispersed)", "glm(visits ~ {f}, family = quasipoisson, control = glm.control(epsilon = 1e-12, maxit = 100))", "F", "glm"),
        ("NMES 1988", nmes, "visits", ["hospital"], ["hospital", "chronic"], "counts that vary a lot", "MASS::glm.nb(visits ~ {f})", "LRT", "anova.negbin"),
        ("birthwt", birthwt, "low", ["lwt"], ["lwt", "smoke"], "yes or no", "glm(low ~ {f}, family = binomial, control = glm.control(epsilon = 1e-12, maxit = 100))", "Chisq", "glm"),
        ("birthwt weight", birthwt, "bwt_kg", ["lwt"], ["lwt", "smoke"], "positive numbers", 'glm(bwt_kg ~ {f}, family = Gamma(link = "log"), control = glm.control(epsilon = 1e-12, maxit = 100))', "F", "glm"),
        ("bioChemists (pscl)", bio, "art", ["ment"], ["ment", "kid5"], "counts with many zeros", "pscl::zeroinfl(art ~ {f})", "lrtest", "zeroinfl")]
if HAS_R:
    for dataset, data, column, using, against, outcome, model, test, topic in GLMS:
        tail = {"lrtest": "a <- lmtest::lrtest(m0, m1); stat <- a$Chisq[2]; p <- a[['Pr(>Chisq)']][2]",
                "LRT": "a <- anova(m0, m1); stat <- a[['LR stat.']][2]; p <- a[['Pr(Chi)']][2]",
                "F": "a <- anova(m0, m1, test = 'F'); stat <- a$F[2]; p <- a[['Pr(>F)']][2]",
                "Chisq": "a <- anova(m0, m1, test = 'Chisq'); stat <- a$Deviance[2]; p <- a[['Pr(>Chi)']][2]"}[test]
        r = r_run(f"m0 <- {model.format(f=' + '.join(using))}; m1 <- {model.format(f=' + '.join(against))}; {tail}\nlist(statistic = stat, `p-value` = p)",
                  data, packages=("MASS", "pscl", "lmtest"))
        refs.run_all(dataset, f'compare_models (outcome="{outcome}")', r, f"R {test} test of the two fits",
                     rman("pscl" if topic == "zeroinfl" else ("MASS" if topic == "anova.negbin" else "stats"), topic), "stats", ITERATIVE)
mammals = r_dataset("mammals", "MASS")
mammals["log_body"], mammals["log_brain"] = np.log(mammals["body"]), np.log(mammals["brain"])
mammals["log_body2"] = mammals["log_body"] ** 2
refs.data("mammals", mammals)
if HAS_R:
    r = r_run('''m1 <- glm(log_brain ~ log_body); m2 <- glm(log_brain ~ log_body + log_body2)
                 t <- AICcmodavg::aictab(list(m1, m2), modnames = c("a", "b"), sort = FALSE)
                 list(aicc = t$AICc, delta = t$Delta_AICc, w = t$AICcWt, ll = t$LL, aic = c(AIC(m1), AIC(m2)))''', mammals, packages=("AICcmodavg",))
    values = {}
    for i, name in enumerate(["log_brain ~ log_body", "log_brain ~ log_body + log_body2"]):
        for column, key in (("AIC", "aic"), ("AICc", "aicc"), ("difference in AICc", "delta"), ("Akaike weight", "w"), ("log-likelihood", "ll")):
            values[f"{name}: {column}"] = r[key][i]
    refs.run_all("mammals (MASS)", "compare_models (using= a list of models: AICc)", {k: v for k, v in values.items() if v != 0},
                 "R AICcmodavg aictab", rman("AICcmodavg", "aictab"), "AICcmodavg", CLOSED)
    refs.run_all("mammals (MASS)", "compare_models (using= a list of models: AICc)", {k: v for k, v in values.items() if v == 0},
                 "R AICcmodavg aictab", rman("AICcmodavg", "aictab"), "AICcmodavg", CLOSED, "absolute")
    for seed in (1, 7):
        r = r_run(f'''d <- data.frame(log_brain = log_brain, log_body = log_body, log_body2 = log_body2)
                      m1 <- glm(log_brain ~ log_body, data = d); m2 <- glm(log_brain ~ log_body + log_body2, data = d)
                      set.seed({seed}); a <- boot::cv.glm(d, m1, K = 10)$delta; set.seed({seed}); b <- boot::cv.glm(d, m2, K = 10)$delta
                      list(a = a, b = b)''', mammals, packages=("boot",))
        values = {}
        for name, key in (("log_brain ~ log_body", "a"), ("log_brain ~ log_body + log_body2", "b")):
            values[f"{name}: error on left-out rows"] = r[key][0]
            values[f"{name}: error on left-out rows (bias-corrected)"] = r[key][1]
        refs.run_all("mammals (MASS)", f"compare_models (check_on_left_out_rows=True, with_seed={seed})", values,
                     "R boot cv.glm (K = 10, same seed; help-page mammals data)", rman("boot", "cv.glm"), "boot", SAME_SEED)

# 11. lasso and elastic net (glmnet vignette data)
if HAS_R:
    import checks

    for name, family in (("QuickStartExample", "gaussian"), ("BinomialExample", "binomial")):
        checks._R.r(f'data("{name}", package = "glmnet")')
        x = np.asarray(checks._R.r(f"{name}$x"))
        n = int(checks._R.r(f"nrow({name}$x)")[0])
        x = x.reshape(-1, order="F").reshape(n, -1, order="F")
        frame = pd.DataFrame(x, columns=[f"V{i + 1}" for i in range(x.shape[1])])
        frame["y"] = np.asarray(checks._R.r(f"as.numeric({name}$y)"))
        refs.data(name, frame)
        columns = [c for c in frame.columns if c != "y"]
        for test_method, alpha, seed in (("lasso", 1, 1), ("elastic net", 0.5, 1), ("lasso", 1, 3)):
            if family == "binomial" and test_method == "elastic net":
                continue
            outcome = "numbers" if family == "gaussian" else "yes or no"
            r = r_run(f'''x <- as.matrix(data.frame({", ".join(columns)})); set.seed({seed})
                          f <- glmnet::cv.glmnet(x, y, family = "{family}", alpha = {alpha}, nfolds = 10)
                          list(b = as.vector(coef(f, s = "lambda.1se")), lambda = f$lambda.1se)''', frame, packages=("glmnet",))
            values = {"penalty chosen (lambda.1se)": r["lambda"], "columns kept (number)": sum(1 for b in r["b"][1:] if b != 0)}
            values.update({f"{c} effect": b for c, b in zip(["(intercept)"] + columns, r["b"]) if b != 0})
            refs.run_all(f"{name} (glmnet)", f'find_columns_that_predict (test_method="{test_method}", outcome="{outcome}", with_seed={seed})',
                         values, "R glmnet cv.glmnet (same seed)", rman("glmnet", "cv.glmnet"), "glmnet", SAME_SEED)

# 12. the groups that pull a mixed model: nothing printed, so R lme4 is run once: influence(m, groups =) refits the model
# without each subject and cooks.distance() says how far the fixed effects move. show_rows_that_pull_model(grouped_by=)
# gives back the groups with a Cook's distance above 4 / (number of groups), the biggest first.
if HAS_R:
    r = r_run("""d <- data.frame(y = Reaction, Days = Days, g = factor(Subject))
                 m <- lme4::lmer(y ~ Days + (1 | g), data = d, REML = TRUE)
                 inf <- stats::influence(m, groups = "g"); cd <- stats::cooks.distance(inf)
                 g <- as.character(inf$deleted); cd <- unname(as.numeric(cd))
                 keep <- which(cd > 4 / length(cd)); keep <- keep[order(-cd[keep])]
                 list(group = g[keep], cooks = cd[keep])""", sleep, packages=("lme4",))
    check = "show_rows_that_pull_model (grouped_by=)"
    refs.run_vector("sleepstudy (lme4)", check, "Subject", list(np.atleast_1d(r["group"])), "R lme4 influence (groups =) and cooks.distance",
                    rman("lme4", "influence.merMod"), "lme4", 0, "text")
    refs.run_vector("sleepstudy (lme4)", check, "Cook's distance", np.atleast_1d(r["cooks"]),
                    "R lme4 influence (groups =) and cooks.distance", rman("lme4", "influence.merMod"), "lme4", ITERATIVE)

refs.save()
