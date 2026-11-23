"""Reference values for notebooks/32_models_planning_more.ipynb (rows left open by notebooks 25 to 31).

Printed first: the pwr vignette (paired t-test and chi-square power). Elsewhere the help-page example of the library
simplydata calls is run once (MASS polr, nnet multinom, betareg with lmtest, emmeans on pscl and glm fits, vegan and
zCompositions, SciPy stats.power with the same seed, R power.prop.test, lavaan with missing values, agricolae designs,
lme4 with crossed groups on lmerTest's carrots, rstanarm with the same seed) and saved with the link to its page.
"""

import numpy as np
import pandas as pd
from scipy import stats

from _tools import References, VIGNETTE, load_csv, r_available, r_dataset, r_run, r_version, rman

refs = References("32_models_planning_more")
HAS_R = r_available()
CLOSED, ITERATIVE, SAME = 1e-10, 1e-6, 1e-12

# 1. compare_models for "ordered groups", "groups" and "shares"; model averages for "counts with many zeros" and "positive numbers"
housing = r_dataset("housing", "MASS")
expanded = housing.loc[housing.index.repeat(housing["Freq"].astype(int))].reset_index(drop=True)
refs.data("housing_one_row_per_tenant", expanded)
gas = r_dataset("GasolineYield", "betareg")
gas["batch"] = "b" + gas["batch"].astype(str).str.zfill(2)
refs.data("GasolineYield", gas)
bio = r_dataset("bioChemists", "pscl")
refs.data("bioChemists", bio)
clot = pd.DataFrame({"u": [5, 10, 15, 20, 30, 40, 60, 80, 100] * 2, "time": [118, 58, 42, 35, 27, 25, 21, 19, 18, 69, 35, 26, 21, 18, 16, 13, 12, 12],
                     "lot": ["lot 1"] * 9 + ["lot 2"] * 9})
refs.data("clotting_two_lots", clot)
if HAS_R:
    HOUSE = '''d <- data.frame(Sat = factor(Sat, levels = c("Low", "Medium", "High"){o}), Infl = factor(Infl), Type = factor(Type), Cont = factor(Cont))\n'''
    for outcome, fit, url in (("ordered groups", "MASS::polr(Sat ~ {f}, data = d)", rman("MASS", "polr")),
                              ("groups", "nnet::multinom(Sat ~ {f}, data = d, trace = FALSE)", rman("nnet", "multinom"))):
        r = r_run(HOUSE.format(o=", ordered = TRUE" if outcome == "ordered groups" else "") +
                  f"m0 <- {fit.format(f='Infl + Type')}; m1 <- {fit.format(f='Infl + Type + Cont')}; a <- anova(m0, m1)\n"
                  "list(statistic = a[['LR stat.']][2], p = a[['Pr(Chi)']][2])", expanded, packages=("MASS", "nnet"))
        refs.run_all("housing (MASS)", f'compare_models (outcome="{outcome}")', {"statistic": r["statistic"], "p-value": r["p"]},
                     "R anova (likelihood-ratio test) of the two fits", url, "MASS" if outcome == "ordered groups" else "nnet", ITERATIVE)
    r = r_run('''d <- data.frame(yield = yield, batch = factor(batch), temp = temp)
                 a <- lmtest::lrtest(betareg::betareg(yield ~ temp, data = d), betareg::betareg(yield ~ batch + temp, data = d))
                 list(statistic = a$Chisq[2], p = a[["Pr(>Chisq)"]][2])''', gas, packages=("betareg", "lmtest"))
    refs.run_all("GasolineYield (betareg)", 'compare_models (outcome="shares")', {"statistic": r["statistic"], "p-value": r["p"]},
                 "R lmtest lrtest of two betareg fits (the betareg vignette compares models this way)", rman("lmtest", "lrtest"), "lmtest", ITERATIVE)
    r = r_run('''d <- data.frame(art = art, mar = factor(mar), ment = ment); m <- pscl::zeroinfl(art ~ mar + ment, data = d)
                 s <- summary(emmeans::emmeans(m, ~ mar)); list(mean = s$emmean, se = s$SE)''', bio, packages=("pscl", "emmeans"))
    refs.run_all("bioChemists (pscl)", 'find_model_averages_of (outcome="counts with many zeros")',
                 {"Married model average": r["mean"][0], "Single model average": r["mean"][1], "Married standard error": r["se"][0],
                  "Single standard error": r["se"][1]}, "R emmeans on pscl zeroinfl", rman("emmeans", "emmeans"), "emmeans", ITERATIVE)
    r = r_run('''d <- data.frame(time = time, lot = factor(lot), log_u = log(u)); m <- glm(time ~ lot + log_u, family = Gamma(link = "log"), data = d,
                                                                          control = glm.control(epsilon = 1e-12, maxit = 100))
                 s <- summary(emmeans::emmeans(m, ~ lot, type = "response")); list(mean = s$response, se = s$SE)''', clot, packages=("emmeans",))
    refs.run_all("clotting times, two lots (R glm help page)", 'find_model_averages_of (outcome="positive numbers")',
                 {"lot 1 model average": r["mean"][0], "lot 2 model average": r["mean"][1], "lot 1 standard error": r["se"][0],
                  "lot 2 standard error": r["se"][1]}, "R emmeans on glm(Gamma, log link)", rman("emmeans", "emmeans"), "emmeans", ITERATIVE)

# 2. crossed groups: lmerTest's carrots (Consumer and Product crossed), the lmerTest help-page data
carrots = r_dataset("carrots", "lmerTest")
carrots = carrots.dropna(subset=["Preference", "sens2", "Homesize"]).reset_index(drop=True)
refs.data("carrots", carrots)
if HAS_R:
    r = r_run('''d <- data.frame(Preference = Preference, sens2 = sens2, Homesize = factor(Homesize), Consumer = factor(Consumer), Product = factor(Product))
                 m <- lme4::lmer(Preference ~ sens2 + Homesize + (1 | Consumer) + (1 | Product), data = d)
                 s <- summary(m)$coefficients; list(b = unname(s[, 1]), se = unname(s[, 2]))''', carrots, packages=("lme4",))
    terms = ["(intercept)", "sens2", "Homesize (3)"]
    refs.run_all("carrots (lmerTest)", 'explain (grouped_by="Consumer", crossed_with="Product")',
                 {**{f"{t} estimate": v for t, v in zip(terms, r["b"])}, **{f"{t} standard error": v for t, v in zip(terms, r["se"])}},
                 "R lme4 lmer with (1 | Consumer) + (1 | Product)", rman("lmerTest", "carrots"), "lme4", ITERATIVE)

# 3. Aitchison principal coordinates (dune, 27 species)
dune = r_dataset("dune", "vegan").drop(columns=["row"], errors="ignore")
DUNE = [c for c in dune.columns if (dune[c] > 0).sum() >= 2]
refs.data("dune", dune)
refs.data("DUNE_27_species", pd.DataFrame({"number": range(1, len(DUNE) + 1), "species": DUNE}))
if HAS_R:
    r = r_run(f'''x <- data.frame({", ".join(DUNE)})
                  z <- zCompositions::cmultRepl(x, method = "GBM", output = "p-counts", z.warning = 1, suppress.print = TRUE)
                  list(p = as.vector(cmdscale(vegan::vegdist(z, method = "aitchison"), k = 2)))''', dune, packages=("vegan", "zCompositions"))
    refs.run_matrix("dune (vegan), 27 species", 'add_principal_coordinates (measure="aitchison", how_many=2)', "coordinates",
                    np.reshape(r["p"], (len(dune), 2), order="F"), "R cmdscale on zCompositions + vegan vegdist (aitchison)", rman("stats", "cmdscale"),
                    "stats", CLOSED, "signs")

# 4. power: the pwr vignette (printed); R power.prop.test; simulated power from your own distributions (SciPy, same seed)
PV = "pwr vignette (Champely et al.), printed output"
PV_URL = VIGNETTE.format(package="pwr", file="pwr-vignette.html")
pwr = ("pwr", r_version("pwr"))
refs.printed("planning settings", 'find_chance_of_finding (for_test="compare_before_after", effect_size=0.32, values_per_group=24, sides=1)', "power",
             "0.4508691", PV + ": pwr.t.test(n = 24, d = 0.08/0.25, type = 'paired', alternative = 'greater')", PV_URL, "printed", *pwr)
refs.printed("planning settings", 'find_how_many_needed (for_test="compare_before_after", effect_size=0.32, sides=1)', "pairs (exact)",
             "61.75209", PV + ": pwr.t.test(d = 0.08/0.25, power = 0.8, type = 'paired', alternative = 'greater')", PV_URL, "printed", *pwr)
refs.printed("planning settings", 'find_chance_of_finding (for_test="compare_counts", effect_size=0.2886751, values_per_group=25, groups=(4, 2))',
             "power", "0.6739834", PV + ": pwr.chisq.test(w = 0.2886751, N = 100, df = 3)", PV_URL, "printed", *pwr)
# The vignette passes w unrounded (an ES.w2 result) and prints it to 7 decimals, so the power it prints carries the
# rounding of w as well: d(power)/dw is about 9 here, so half a unit of w's last digit (5e-8) moves the power by up to
# 4.5e-7. The tolerance is half a unit of the printed power plus that.
refs.rows[-1]["tolerance"] = 5e-8 + 9 * 5e-8
refs.rows[-1]["reason"] = ""
refs.printed("planning settings", 'find_how_many_needed (for_test="compare_counts", effect_size=0.2886751, groups=(4, 2))', "values in all (exact)",
             "130.8308", PV + ": pwr.chisq.test(w = 0.2886751, df = 3, power = 0.8)", PV_URL, "printed", *pwr)
refs.printed("planning settings", 'find_chance_of_finding (for_test="compare_counts", effect_size=0.2182179, values_per_group=50, groups=(2, 2), significance=0.01)',
             "power", "0.3469206", PV + ": pwr.chisq.test(w = 0.2182179, N = 100, df = 1, sig.level = 0.01)", PV_URL, "printed", *pwr)
if HAS_R:
    r = r_run('''a <- power.prop.test(p1 = 0.30, p2 = 0.15, power = 0.8)$n; b <- power.prop.test(n = 120, p1 = 0.30, p2 = 0.15)$power
                 list(a = a, b = b)''', None)
    refs.run("planning settings", 'find_how_many_needed (for_test="compare_counts", percents=[30, 15])', "values per group (exact)", r["a"],
             "R stats power.prop.test", rman("stats", "power.prop.test"), "stats", ITERATIVE)
    refs.run("planning settings", 'find_chance_of_finding (for_test="compare_counts", percents=[30, 15], values_per_group=120)', "power", r["b"],
             "R stats power.prop.test (one tail, R's default)", rman("stats", "power.prop.test"), "stats", ITERATIVE)
from simplydata.power_more import BATCH, EXACT_LIMIT  # noqa: E402  simplydata's settings for the call (not results)


def welch(x, y, axis=-1):
    return stats.ttest_ind(x, y, equal_var=False, axis=axis).pvalue


for n in (20, 40):
    rng = np.random.default_rng(1)
    a, b = stats.lognorm(s=0.5, scale=np.exp(0.0)), stats.lognorm(s=0.5, scale=np.exp(0.4))
    power = stats.power(welch, (lambda size: a.rvs(size=size, random_state=rng), lambda size: b.rvs(size=size, random_state=rng)),
                        (np.array([n]), np.array([n])), significance=0.05, n_resamples=2000, batch=BATCH).power[0]
    refs.run("planning settings", "find_chance_of_finding (drawing_from=[two lognormal distributions], values_per_group=[20, 40], tries=2000, with_seed=1)",
             f"power with {n} per group", power, "SciPy stats.power (Welch's t-test, the same distributions and seed)",
             "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.stats.power.html", "scipy", SAME)

# 5. confirmatory factor analysis with empty cells (psych::bfi, which has empty cells): lavaan, missing = "ml"
bfi = r_dataset("bfi", "psych")[["A1", "A2", "A3", "A4", "A5", "C1", "C2", "C3", "C4", "C5"]]
refs.data("bfi", bfi)
if HAS_R:
    r = r_run('''f <- lavaan::cfa("agree =~ A1 + A2 + A3 + A4 + A5; consc =~ C1 + C2 + C3 + C4 + C5",
                                  data = data.frame(A1, A2, A3, A4, A5, C1, C2, C3, C4, C5), missing = "ml")
                 m <- lavaan::fitMeasures(f, c("chisq", "df", "cfi", "rmsea", "srmr"))
                 s <- lavaan::standardizedSolution(f); k <- s$op == "=~"
                 list(chisq = m[["chisq"]], df = m[["df"]], cfi = m[["cfi"]], rmsea = m[["rmsea"]], std = s$est.std[k], rhs = s$rhs[k])''', bfi, packages=("lavaan",))
    values = {"chi-square": r["chisq"], "degrees of freedom": r["df"], "CFI": r["cfi"], "RMSEA": r["rmsea"]}
    values.update({f"standardised loading of {v}": s for v, s in zip(r["rhs"], r["std"])})
    refs.run_all("bfi (psych)", "check_hidden_factors (measured_by=), empty cells kept (FIML)", values, "R lavaan cfa (missing = 'ml')",
                 rman("lavaan", "cfa"), "lavaan", ITERATIVE)

# 6. other kinds of design table (R agricolae, same seed)
if HAS_R:
    for kind, code, columns in (
            ("completely random", 'agricolae::design.crd(c("x", "y", "z"), r = 4, serie = 0, seed = SEED, kinds = "Mersenne-Twister")$book', ["f1"]),
            ("randomised blocks", 'agricolae::design.rcbd(c("x", "y", "z"), r = 4, serie = 0, seed = SEED, kinds = "Mersenne-Twister")$book', ["block", "f1"]),
            ("latin square", 'agricolae::design.lsd(c("x", "y", "z", "w"), serie = 0, seed = SEED, kinds = "Mersenne-Twister")$book', ["row", "col", "f1"])):
        for seed in (2, 9):
            r = r_run(f'''b <- {code.replace("SEED", str(seed))}
                          list(f1 = as.character(b[[ncol(b)]]), block = if ("block" %in% names(b)) as.character(b$block) else "",
                               row = if ("row" %in% names(b)) as.character(b$row) else "", col = if ("col" %in% names(b)) as.character(b$col) else "")''',
                      None, packages=("agricolae",))
            check = f'make_design_table (of_kind="{kind}", repeats=4, with_seed={seed})' if kind != "latin square" else f'make_design_table (of_kind="latin square", with_seed={seed})'
            for c in columns:
                refs.run_vector("designs (same seed)", check, c, r[c], f"R agricolae (same seed)", rman("agricolae", code.split("::")[1].split("(")[0]),
                                "agricolae", 0, "text")

# 7. the Bayesian model (rstanarm, the same seed): the posterior medians
warp = load_csv("warpbreaks.csv")
if HAS_R:
    r = r_run('''d <- data.frame(y = breaks, wool = factor(wool), tension = factor(tension))
                 m <- rstanarm::stan_glm(y ~ wool + tension, data = d, family = stats::poisson(), chains = 4L, iter = 2000L, seed = 1, refresh = 0, cores = 1L)
                 list(b = unname(apply(as.matrix(m), 2, median))[1:4])''', warp, packages=("rstanarm",))
    refs.run_all("warpbreaks", 'explain (outcome="counts", test_method="bayesian", with_seed=1)',
                 {f"{t} posterior median (log scale)": v for t, v in zip(["(intercept)", "wool (B)", "tension (L)", "tension (M)"], r["b"])},
                 "R rstanarm stan_glm (the same seed)", rman("rstanarm", "stan_glm"), "rstanarm", 1e-10)

# 8. check_model_of(test_method="bayesian"): show="numbers" (R-hat, effective draws, divergences and the posterior
# predictive check) and show="prior check" (the prior predictive picture, read with chart.get_numbers_drawn()).
# Nothing is printed for these, so R is run once with the calls simplydata makes, with the same seed: stan_glm (4 chains
# of 2000, seed 1), posterior_predict(draws = 50, seed = 1), posterior::summarise_draws, and the count of each value
# in the data against the median and 90 % interval of its count in the 50 made-up sets (as bayesplot::ppc_bars).
if HAS_R:
    for prior_only in (False, True):
        code = f'''d <- data.frame(y = breaks, wool = factor(wool), tension = factor(tension))
                     m <- rstanarm::stan_glm(y ~ wool + tension, data = d, family = stats::poisson(), chains = 4L, iter = 2000L,
                                             seed = 1, refresh = 0, cores = 1L{", prior_PD = TRUE" if prior_only else ""})
                     yrep <- rstanarm::posterior_predict(m, draws = 50L, seed = 1)
                     s <- posterior::summarise_draws(posterior::as_draws_array(as.array(m)), "rhat", "ess_bulk", "ess_tail")
                     lv <- sort(unique(c(d$y, as.vector(yrep))))
                     observed <- as.vector(table(factor(d$y, levels = lv)))
                     k <- t(apply(yrep, 1, function(v) as.vector(table(factor(v, levels = lv)))))
                     made <- t(apply(k, 2, stats::quantile, probs = c(0.05, 0.5, 0.95), names = FALSE))
                     list(rhat = s$rhat, bulk = s$ess_bulk, tail = s$ess_tail, divergent = rstan::get_num_divergent(m$stanfit),
                          grid = lv, observed = observed, low = made[, 1], middle = made[, 2], high = made[, 3])'''
        # run inside local() and withCallingHandlers(), as simplydata runs its R code: rstanarm's posterior_predict looks
        # at the function that called it (deparse(sys.function())), which fails when it is called from the top level here
        r = r_run("""local({ out <- NULL; utils::capture.output(out <- withCallingHandlers(eval(parse(text = sd_code)),
                       warning = function(w) invokeRestart("muffleWarning"), message = function(m) invokeRestart("muffleMessage")))
                     out })""", {**{c: warp[c] for c in warp.columns}, "sd_code": [code]}, packages=("rstanarm",))
        source, url = "R rstanarm posterior_predict and posterior summarise_draws (the same seed)", rman("rstanarm", "posterior_predict")
        if not prior_only:
            check = 'check_model_of (test_method="bayesian", show="numbers")'
            for key, name in (("rhat", "R-hat"), ("bulk", "effective draws in the middle (bulk)"), ("tail", "effective draws in the tails")):
                refs.run_vector("warpbreaks", check, name, np.atleast_1d(r[key]), source, url, "rstanarm", 1e-10)
            refs.run("warpbreaks", check, "divergent transitions", r["divergent"], source, url, "rstanarm", 0, "absolute")
            every = 1
            names = {"grid": "posterior predictive check: value", "observed": "posterior predictive check: count in the data",
                     "middle": "posterior predictive check: median count in the made-up sets",
                     "low": "posterior predictive check: 5 % of the made-up sets", "high": "posterior predictive check: 95 % of the made-up sets"}
        else:
            check = 'check_model_of (test_method="bayesian", show="prior check")'
            every = max(1, len(np.atleast_1d(r["grid"])) // 60)
            names = {"grid": "drawn bars: value", "observed": "drawn bars: count in the data",
                     "middle": "drawn dots: median count in the made-up sets",
                     "low": "drawn error bars: 5 % of the made-up sets", "high": "drawn error bars: 95 % of the made-up sets"}
            refs.run("warpbreaks", check, "drawn bars (number)", len(np.atleast_1d(r["grid"])), source, url, "rstanarm", 0, "absolute")
        for key, name in names.items():
            refs.run_vector("warpbreaks", check, name, np.atleast_1d(r[key]), source, url, "rstanarm", 1e-9, "absolute", every=every)

refs.save()
