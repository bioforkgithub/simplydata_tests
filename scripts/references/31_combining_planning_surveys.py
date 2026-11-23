"""Reference values for notebooks/31_combining_planning_surveys.ipynb.

Printed first: Viechtbauer (2010) for the BCG random-effects model (the metafor JSS article); the pwr vignette for two
sample sizes; Lumley (2010) for the stratified api mean. Elsewhere the help-page example of the library simplydata calls
is run once (metafor escalc/rma/rma.peto/conv.wald, R power.t.test, gsDesign, powerSurvEpi, pwr, pwrss, TOSTER, SciPy
stats.power with the same seed, survey, agricolae, blockrand) and saved with the link to its page.
"""

import numpy as np
from scipy import stats

from _tools import References, VIGNETTE, load_csv, r_available, r_run, r_version, rman

refs = References("31_combining_planning_surveys")
HAS_R = r_available()
CLOSED, ITERATIVE, SAME = 1e-10, 1e-6, 1e-12
MF = rman("metafor", "rma.uni")
TIGHT = ", control = list(tol = 1e-12, threshold = 1e-12, maxiter = 1000)"   # metafor run to convergence

# 1. combining studies
V10 = "Viechtbauer (2010), Conducting meta-analyses in R with the metafor package, J. Stat. Softw. 36(3), section 4"
refs.printed_all("BCG trials (metadat)", 'combine_studies (show="risk ratios", test_method="reml", interval_by="normal")',
                 {"estimate": "-0.7145", "standard error": "0.1798", "low": "-1.0669", "high": "-0.3622", "tau squared": "0.3132",
                  "I squared %": "92.22", "Q": "152.2330"}, V10, "https://doi.org/10.18637/jss.v036.i03", "published", "metafor", r_version("metafor"))
bcg = load_csv("dat_bcg.csv")
hack = load_csv("dat_hackshaw1998.csv")
if HAS_R:
    r = r_run('''y <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg)
                 m <- metafor::rma(yi, vi, data = y, method = "REML", test = "adhoc")
                 list(estimate = m$b[[1]], `standard error` = m$se, low = m$ci.lb, high = m$ci.ub, `p-value` = m$pval, `tau squared` = m$tau2)''', bcg,
              packages=("metafor",))
    refs.run_all("BCG trials (metadat)", 'combine_studies (show="risk ratios", test_method="reml", interval_by="Hartung-Knapp")', r,
                 "R metafor rma (REML, Hartung-Knapp with the ad hoc rule, test = adhoc), help-page BCG data", MF, "metafor", ITERATIVE)
    r = r_run('''m <- metafor::rma.peto(ai = tpos, bi = tneg, ci = cpos, di = cneg)
                 list(estimate = m$b[[1]], `standard error` = m$se, low = m$ci.lb, high = m$ci.ub)''', bcg, packages=("metafor",))
    refs.run_all("BCG trials (metadat)", 'combine_studies (show="odds ratios", test_method="peto")', r, "R metafor rma.peto (help page example)",
                 rman("metafor", "rma.peto"), "metafor", CLOSED)
    CASES = [("Hackshaw (1997), passive smoking", "dat_hackshaw1998.csv",
              "y <- metafor::conv.wald(out = or, ci.lb = or.lb, ci.ub = or.ub, data = data.frame(or, or.lb, or.ub), transf = log)",
              "combine_studies (estimates= with interval columns, a ratio)", rman("metafor", "conv.wald")),
             ("Hackshaw (1997), passive smoking", "dat_hackshaw1998.csv", "y <- data.frame(yi = yi, vi = vi)",
              "combine_studies (standard_errors=, estimates_are_ratios=False)", MF),
             ("Curtis et al. (1998), CO2", "dat_curtis1998.csv",
              "y <- metafor::escalc(measure = 'ROM', m1i = m1i, sd1i = sd1i, n1i = n1i, m2i = m2i, sd2i = sd2i, n2i = n2i)",
              'combine_studies (averages=, standard_deviations=, show="ratios of averages")', rman("metafor", "escalc")),
             ("McDaniel et al. (1994), interviews", "dat_mcdaniel1994.csv", "y <- metafor::escalc(measure = 'ZCOR', ri = ri, ni = ni)",
              "combine_studies (correlations=, counts=)", rman("metafor", "escalc")),
             ("Hart et al. (1999), warfarin", "dat_hart1999.csv", "y <- metafor::escalc(measure = 'IRR', x1i = x1i, t1i = t1i, x2i = x2i, t2i = t2i)",
              "combine_studies (events=, time=)", rman("metafor", "escalc"))]
    for dataset, file, escalc, check, url in CASES:
        data = load_csv(file)
        if "hackshaw" in file:
            data = data.rename(columns={"or.lb": "or.lb", "or.ub": "or.ub"})
        back = "metafor::transf.ztor" if "ZCOR" in escalc else "identity"   # correlations: the estimate and interval turned back
        # metafor's default stopping rule (tol and threshold 1e-5), and run to convergence (TIGHT)
        for control, source in (("", "R metafor escalc + rma (REML, test = adhoc)"),
                                (TIGHT, "R metafor escalc + rma (REML, test = adhoc, tol 1e-12)")):
            r = r_run(escalc + f'''
                  m <- metafor::rma(yi, vi, data = y, method = "REML", test = "adhoc"{control}); p <- predict(m, transf = {back})
                  list(estimate = p$pred, `standard error` = m$se, low = p$ci.lb, high = p$ci.ub, `tau squared` = m$tau2)''', data, packages=("metafor",))
            refs.run_all(dataset, check, r, source, url, "metafor", ITERATIVE)

# 3. planning: the pwr vignette (printed) and the R functions simplydata calls (run)
PV = "pwr vignette (Champely et al.), printed output"
PV_URL = VIGNETTE.format(package="pwr", file="pwr-vignette.html")
refs.printed("planning settings", "find_how_many_needed (difference=0.75, spread=2.25)", "values per group (exact)", "142.2466",
             PV + ": power.t.test(delta = 0.75, sd = 2.25, power = 0.8)", PV_URL, "printed", "pwr", r_version("pwr"))
refs.printed("planning settings", 'find_how_many_needed (for_test="compare_average_to", sides=1, power=0.9)', "values (exact)", "174.7796",
             PV + ": pwr.t.test(d = 0.5/2.25, power = 0.9, alternative = 'greater', type = 'one.sample')", PV_URL, "printed", "pwr", r_version("pwr"))
if HAS_R:
    r = r_run('''list(c = power.t.test(delta = 0.5, power = 0.8, sig.level = 0.01)$n, d = power.t.test(n = 30, delta = 0.5, sig.level = 0.01)$power,
                      cl = power.t.test(n = 100 / 1.95, delta = 0.5)$power)''', None)
    refs.run("planning settings", "find_how_many_needed (effect_size=0.5, significance=0.01)", "values per group (exact)", r["c"], "R stats power.t.test",
             rman("stats", "power.t.test"), "stats", ITERATIVE)
    refs.run("planning settings", "find_chance_of_finding (effect_size=0.5, values_per_group=30, significance=0.01)", "power", r["d"], "R stats power.t.test",
             rman("stats", "power.t.test"), "stats", CLOSED)
    refs.run("planning settings", "find_chance_of_finding (values_per_group=100, in_clusters_of=20, alike_within_clusters=0.05)", "power", r["cl"],
             "R stats power.t.test at n / (1 + (m - 1) ICC)", rman("stats", "power.t.test"), "stats", CLOSED)
    r = r_run('''e <- gsDesign::nEvents(hr = 0.7, alpha = 0.025, beta = 0.2)
                 f <- powerSurvEpi::ssizeCT.default(power = 0.8, k = 1, pE = 0.30, pC = 0.40, RR = 0.7, alpha = 0.05)
                 list(e = e, n = max(f[[1]], f[[2]]))''', None, packages=("gsDesign", "powerSurvEpi"))
    refs.run("planning settings", 'find_how_many_needed (for_test="compare_survival", effect_size=0.7, percents=[30, 40])', "events (exact)", r["e"],
             "R gsDesign nEvents (Schoenfeld)", rman("gsDesign", "nEvents"), "gsDesign", CLOSED)
    refs.run("planning settings", 'find_how_many_needed (for_test="compare_survival", test_method="freedman")', "people per group", r["n"],
             "R powerSurvEpi ssizeCT.default (Freedman)", rman("powerSurvEpi", "ssizeCT.default"), "powerSurvEpi", 0, "absolute")
    for test, u, total in (("explain", 3, 3), ("compare_models", 2, 3)):
        r = r_run(f"p <- pwr::pwr.f2.test(u = {u}, f2 = 0.15, power = 0.8); list(n = p$v + {total} + 1)", None, packages=("pwr",))
        label = f'find_how_many_needed (for_test="{test}", effect_size=0.15' + (", against=)" if test == "compare_models" else ")")
        refs.run("planning settings", label, "rows (exact)", r["n"], "R pwr pwr.f2.test", rman("pwr", "pwr.f2.test"), "pwr", ITERATIVE)
    r = r_run("p <- pwr::pwr.f2.test(u = 3, v = 77 - 3 - 1, f2 = 0.15); list(p = p$power)", None, packages=("pwr",))
    refs.run("planning settings", 'find_chance_of_finding (for_test="explain", effect_size=0.15, using=3, values_per_group=77)', "power", r["p"],
             "R pwr pwr.f2.test", rman("pwr", "pwr.f2.test"), "pwr", CLOSED)
    r = r_run('''a <- pwrss::power.z.logistic(base.prob = 0.2, odds.ratio = 1.5, power = 0.8, alpha = 0.05, distribution = "normal", verbose = FALSE)
                 b <- pwrss::power.z.logistic(base.prob = 0.2, odds.ratio = 1.5, n = 200, alpha = 0.05, distribution = "normal", verbose = FALSE)
                 list(n = a$n, p = b$power)''', None, packages=("pwrss",))
    refs.run("planning settings", 'find_how_many_needed (for_test="explain", outcome="yes or no", percents=[20], effect_size=1.5)', "values (rounded up)",
             r["n"], "R pwrss power.z.logistic (n rounded up)", rman("pwrss", "power.z.logistic"), "pwrss", 0, "absolute")
    refs.run("planning settings", 'find_chance_of_finding (for_test="explain", outcome="yes or no", percents=[20], effect_size=1.5, values_per_group=200)',
             "power", r["p"], "R pwrss power.z.logistic", rman("pwrss", "power.z.logistic"), "pwrss", ITERATIVE)
    r = r_run("p <- TOSTER::power_t_TOST(delta = 0, sd = 10, eqb = 5, power = 0.8, alpha = 0.05, type = 'two.sample'); list(n = p$n)", None, packages=("TOSTER",))
    refs.run("planning settings", "find_how_many_needed (same_within=5, spread=10)", "values per group (exact)", r["n"], "R TOSTER power_t_TOST",
             rman("TOSTER", "power_t_TOST"), "TOSTER", ITERATIVE)

# rank power by simulation: SciPy stats.power (the call simplydata makes) with a fresh generator default_rng(1) for each size
from simplydata.power_more import BATCH, EXACT_LIMIT  # noqa: E402  simplydata's settings for the call (not results)


def mann_whitney(x, y, axis=-1):
    return stats.mannwhitneyu(x, y, axis=axis, method="exact" if x.shape[-1] < EXACT_LIMIT else "asymptotic", use_continuity=True).pvalue


for n in (20, 30):
    rng = np.random.default_rng(1)
    power = stats.power(mann_whitney, (lambda size: rng.standard_normal(size=size), lambda size: rng.normal(0.5, 1.0, size=size)),
                        (np.array([n]), np.array([n])), significance=0.05, n_resamples=2000, batch=BATCH).power[0]
    refs.run("planning settings", 'find_chance_of_finding (based_on="ranks", values_per_group=[20, 30], tries=2000, with_seed=1)', f"power with {n} per group",
             power, "SciPy stats.power (Mann-Whitney, same seed)", "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.power.html",
             "scipy", SAME)

# 4. surveys: Lumley (2010) printed; R survey svymean and svyglm (help page api designs)
LUM = "Lumley (2010), Complex Surveys: A Guide to Analysis Using R, section 2.3 (svymean of api00 in the stratified sample)"
refs.printed_all("apistrat (survey)", 'use_survey_design (weights="pw", strata="stype", population_sizes="fpc"), find_average_of',
                 {"average of api00": "662.29", "standard error": "9.4089"}, LUM, "https://doi.org/10.1002/9780470580066", "published", "survey", r_version("survey"))
DESIGNS = [("apistrat (survey)", "apistrat.csv", 'weights="pw"', "ids = ~1, weights = ~pw"),
           ("apistrat (survey)", "apistrat.csv", 'weights="pw", strata="stype"', "ids = ~1, strata = ~stype, weights = ~pw"),
           ("apistrat (survey)", "apistrat.csv", 'weights="pw", strata="stype", population_sizes="fpc"', "ids = ~1, strata = ~stype, weights = ~pw, fpc = ~fpc"),
           ("apiclus1 (survey)", "apiclus1.csv", 'weights="pw", clusters="dnum"', "ids = ~dnum, weights = ~pw"),
           ("apiclus1 (survey)", "apiclus1.csv", 'weights="pw", clusters="dnum", population_sizes="fpc"', "ids = ~dnum, weights = ~pw, fpc = ~fpc")]
if HAS_R:
    for dataset, file, options, design in DESIGNS:
        r = r_run(f'''d <- survey::svydesign({design}, data = data.frame(api00, pw, stype, dnum, fpc)); m <- survey::svymean(~api00, d)
                      list(mean = coef(m)[[1]], se = survey::SE(m)[[1]])''', load_csv(file), packages=("survey",))
        refs.run_all(dataset, f"use_survey_design ({options}), find_average_of", {"average of api00": r["mean"], "standard error": r["se"]},
                     "R survey svymean (help page designs)", rman("survey", "svydesign"), "survey", CLOSED)
    r = r_run('''d <- survey::svydesign(ids = ~1, strata = ~stype, weights = ~pw, fpc = ~fpc, data = data.frame(api00, pw, stype, fpc, ell, meals))
                 s <- summary(survey::svyglm(api00 ~ ell + meals, design = d))$coefficients
                 list(b = unname(s[, 1]), se = unname(s[, 2]), p = unname(s[, 4]))''', load_csv("apistrat.csv"), packages=("survey",))
    refs.run_all("apistrat (survey)", "use_survey_design (weights=, strata=, population_sizes=), explain",
                 {f"{t} {q}": r[k][i] for i, t in enumerate(["(Intercept)", "ell", "meals"]) for q, k in (("estimate", "b"), ("standard error", "se"), ("p-value", "p"))},
                 "R survey svyglm", rman("survey", "svyglm"), "survey", CLOSED)
    # 5. designs and randomisation lists with the same seed
    for seed in (4, 11):
        r = r_run(f'''d <- agricolae::design.split(c("1", "2"), c("x", "y", "z"), r = 3, design = "rcbd", serie = 0, seed = {seed}, kinds = "Mersenne-Twister")$book
                      list(a = as.numeric(as.character(d[[4]])), b = as.character(d[[5]]), block = as.numeric(as.character(d$block)))''', None, packages=("agricolae",))
        check = f'make_design_table (of_kind="split plot", whole_plots="a", repeats=3, with_seed={seed})'
        refs.run_vector("designs (same seed)", check, "a", r["a"], "R agricolae design.split (same seed)", rman("agricolae", "design.split"), "agricolae", 0, "absolute")
        refs.run_vector("designs (same seed)", check, "b", r["b"], "R agricolae design.split (same seed)", rman("agricolae", "design.split"), "agricolae", 0, "text")
        refs.run_vector("designs (same seed)", check, "block", r["block"], "R agricolae design.split (same seed)", rman("agricolae", "design.split"), "agricolae",
                        0, "absolute")
    for how_many, seed in ((20, 3), (30, 8)):
        r = r_run(f'''set.seed({seed}, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
                      b <- blockrand::blockrand(n = {how_many}, num.levels = 2, levels = c("A", "B"), block.sizes = c(2, 3))
                      list(group = as.character(b$treatment), size = as.numeric(b$block.size))''', None, packages=("blockrand",))
        check = f"make_random_assignments (how_many={how_many}, in_blocks_of=[4, 6], with_seed={seed})"
        refs.run_vector("designs (same seed)", check, "group", r["group"], "R blockrand (same seed)", rman("blockrand", "blockrand"), "blockrand", 0, "text")
        refs.run_vector("designs (same seed)", check, "block size", r["size"], "R blockrand (same seed)", rman("blockrand", "blockrand"), "blockrand", 0, "absolute")

# 6. the numbers the forest and funnel plots draw (chart.get_numbers_drawn()) on the BCG trials: each study's risk ratio
# and its 95 % interval (metafor escalc, normal interval), the line at 1, the combined risk ratio with its interval (the
# diamond; REML with Hartung-Knapp, test = adhoc, as simplydata's default) and the prediction interval; the funnel: each
# study at (risk ratio, standard error), the line at the combined effect and the pseudo 95 % limits, as metafor funnel()
if HAS_R:
    r = r_run('''y <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg); s <- summary(y)
                 m <- metafor::rma(yi, vi, data = y, method = "REML", test = "adhoc"); p <- predict(m, transf = exp); z <- qnorm(0.975)
                 se <- sqrt(y$vi); top <- max(se)
                 list(dot = exp(as.numeric(y$yi)), low = exp(s$ci.lb), high = exp(s$ci.ub), se = se, est = p$pred, lb = p$ci.lb, ub = p$ci.ub,
                      pilb = p$pi.lb, piub = p$pi.ub, centre = exp(m$b[[1]]),
                      fx = exp(c(m$b[[1]] - z * top, m$b[[1]], m$b[[1]] + z * top)), fy = c(top, 0, top))''', bcg, packages=("metafor",))
    D, F = "chart forest_plot (line_at=1, log_scale=)", "chart draw_funnel_plot (log_scale=)"
    ES, RMA, FUN = rman("metafor", "escalc"), MF, rman("metafor", "funnel")
    refs.run_vector("BCG trials (metadat)", D, "study dot (risk ratio)", r["dot"], "R metafor escalc (RR)", ES, "metafor", CLOSED)
    refs.run_vector("BCG trials (metadat)", D, "study interval low", r["low"], "R metafor escalc + summary (95 % normal interval)", ES, "metafor", CLOSED)
    refs.run_vector("BCG trials (metadat)", D, "study interval high", r["high"], "R metafor escalc + summary (95 % normal interval)", ES, "metafor", CLOSED)
    refs.run_vector("BCG trials (metadat)", D, "line at", [1.0, 1.0], "line_at=1", RMA, "metafor", CLOSED)
    refs.run_all("BCG trials (metadat)", D, {"diamond: combined risk ratio": r["est"], "diamond: left end (95 % low)": r["lb"],
                                              "diamond: right end (95 % high)": r["ub"], "prediction interval low": r["pilb"],
                                              "prediction interval high": r["piub"]},
                 "R metafor rma (REML, Hartung-Knapp, test = adhoc) + predict(transf = exp)", rman("metafor", "predict.rma"), "metafor", ITERATIVE)
    refs.run_vector("BCG trials (metadat)", F, "study dot across (risk ratio)", r["dot"], "R metafor escalc (RR)", ES, "metafor", CLOSED)
    refs.run_vector("BCG trials (metadat)", F, "study dot up (standard error)", r["se"], "R metafor escalc, sqrt(vi)", ES, "metafor", CLOSED)
    refs.run_vector("BCG trials (metadat)", F, "line at the combined effect", [r["centre"], r["centre"]], "R metafor rma (REML) + funnel()", FUN,
                    "metafor", ITERATIVE)
    refs.run_vector("BCG trials (metadat)", F, "funnel across", r["fx"], "R metafor funnel(): combined effect +/- qnorm(0.975) x standard error",
                    FUN, "metafor", ITERATIVE)
    refs.run_vector("BCG trials (metadat)", F, "funnel up", r["fy"], "R metafor funnel(): from the largest standard error to 0", FUN, "metafor", CLOSED)

refs.save()
