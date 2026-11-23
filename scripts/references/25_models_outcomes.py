"""Reference values for notebooks/25_models_outcomes.ipynb.

Printed first: Zeileis, Kleiber & Jackman (2008) for the NMES 1988 count models; Dobson & Barnett (2008) for the
probit and complementary log-log beetle fits; Ferrari & Cribari-Neto (2004) Table 1 as printed in the betareg vignette;
McCullagh & Nelder (1989) for the clotting Gamma fit; the quantreg vignette for Engel's median regression; the
statsmodels RLM docstring for the stack loss data. Where nothing is printed, the help-page example of the library
simplydata calls is run once (R glm/sandwich, pscl, confint, glm Gamma, car Anova and anova, effectsize, MASS polr,
nnet multinom, quantreg rq, effectsize, lm, predict) and saved with the link to its page.
"""

import numpy as np

from _tools import References, VIGNETTE, load_csv, py_version, r_available, r_dataset, r_run, r_version, rman

refs = References("25_models_outcomes")
HAS_R = r_available()
CLOSED, ITERATIVE = 1e-10, 1e-6

# 1. NMES 1988: Zeileis, Kleiber & Jackman (2008), J. Stat. Softw. 27(8) (the pscl vignette "countreg")
ZKJ = "Zeileis, Kleiber & Jackman (2008), J. Stat. Softw. 27(8)"
ZKJ_URL = "https://doi.org/10.18637/jss.v027.i08 ; " + VIGNETTE.format(package="pscl", file="countreg.pdf")
pscl = ("pscl", r_version("pscl"))
POISSON = {"(intercept)": ("1.028874", "0.023785"), "hospital": ("0.164797", "0.005997"), "health (poor)": ("0.248307", "0.017845"),
           "health (excellent)": ("-0.361993", "0.030304"), "chronic": ("0.146639", "0.004580"),
           "gender (male)": ("-0.112320", "0.012945"), "school": ("0.026143", "0.001843"), "insurance (yes)": ("0.201687", "0.016860")}
for term, (b, s) in POISSON.items():
    refs.printed("NMES 1988", 'explain (outcome="counts")', f"{term} estimate", b, ZKJ + ", summary(fm_pois)", ZKJ_URL, "published", *pscl)
    refs.printed("NMES 1988", 'explain (outcome="counts")', f"{term} standard error", s, ZKJ + ", summary(fm_pois)", ZKJ_URL, "published", *pscl)
refs.printed("NMES 1988", 'explain (outcome="counts (overdispersed)")', "dispersion", "6.706", ZKJ + ", section 3.3", ZKJ_URL, "published", *pscl)
for term, s in {"(intercept)": "0.062", "hospital": "0.016", "health (poor)": "0.046", "health (excellent)": "0.078", "chronic": "0.012",
                "gender (male)": "0.034", "school": "0.005", "insurance (yes)": "0.044"}.items():
    refs.printed("NMES 1988", 'explain (outcome="counts (overdispersed)")', f"{term} standard error", s, ZKJ + ", Table 2", ZKJ_URL,
                 "published", *pscl)
NB = {"(intercept)": ("0.929", "0.055"), "hospital": ("0.218", "0.020"), "health (poor)": ("0.305", "0.049"),
      "health (excellent)": ("-0.342", "0.061"), "chronic": ("0.175", "0.012"), "gender (male)": ("-0.126", "0.031"),
      "school": ("0.027", "0.004"), "insurance (yes)": ("0.224", "0.039")}
for term, (b, s) in NB.items():
    refs.printed("NMES 1988", 'explain (outcome="counts that vary a lot")', f"{term} estimate", b, ZKJ + ", Table 2", ZKJ_URL, "published", *pscl)
    refs.printed("NMES 1988", 'explain (outcome="counts that vary a lot")', f"{term} standard error", s, ZKJ + ", Table 2", ZKJ_URL,
                 "published", *pscl)
HURDLE = {
    "count": {"(intercept)": ("1.197699", "0.058973"), "hospital": ("0.211898", "0.021396"), "health (poor)": ("0.315958", "0.048056"),
              "health (excellent)": ("-0.331861", "0.066093"), "chronic": ("0.126421", "0.012452"),
              "gender (male)": ("-0.068317", "0.032416"), "school": ("0.020693", "0.004535"), "insurance (yes)": ("0.100172", "0.042619")},
    "more than 0": {"(intercept)": ("0.043147", "0.139852"), "hospital": ("0.312449", "0.091437"),
                    "health (poor)": ("-0.008716", "0.161024"), "health (excellent)": ("-0.289570", "0.142682"),
                    "chronic": ("0.535213", "0.045378"), "gender (male)": ("-0.415658", "0.087608"),
                    "school": ("0.058541", "0.011989"), "insurance (yes)": ("0.747120", "0.100880")}}
for part, values in HURDLE.items():
    for term, (b, s) in values.items():
        refs.printed("NMES 1988", 'explain (outcome="counts in two steps that vary a lot")', f"{part}: {term} estimate", b,
                     ZKJ + ", summary(fm_hurdle0)", ZKJ_URL, "published", *pscl)
        refs.printed("NMES 1988", 'explain (outcome="counts in two steps that vary a lot")', f"{part}: {term} standard error", s,
                     ZKJ + ", summary(fm_hurdle0)", ZKJ_URL, "published", *pscl)
ZINB = {
    "count": {"(intercept)": ("1.194", "0.057"), "hospital": ("0.201", "0.020"), "health (poor)": ("0.285", "0.045"),
              "health (excellent)": ("-0.319", "0.060"), "chronic": ("0.129", "0.012"), "gender (male)": ("-0.080", "0.031"),
              "school": ("0.021", "0.004"), "insurance (yes)": ("0.126", "0.042")},
    "extra zeros": {"(intercept)": ("-0.047", "0.269"), "hospital": ("-0.800", "0.421"), "chronic": ("-1.248", "0.178"),
                    "insurance (yes)": ("-1.176", "0.220"), "school": ("-0.084", "0.026"), "gender (male)": ("0.648", "0.200")}}
for part, values in ZINB.items():
    for term, (b, s) in values.items():
        check = 'explain (outcome="counts with many zeros that vary a lot", zero_part_using=)'
        refs.printed("NMES 1988", check, f"{part}: {term} estimate", b, ZKJ + ", Table 2 (fm_zinb)", ZKJ_URL, "published", *pscl)
        refs.printed("NMES 1988", check, f"{part}: {term} standard error", s, ZKJ + ", Table 2 (fm_zinb)", ZKJ_URL, "published", *pscl)

nmes = load_csv("nmes1988.csv")
ORDER = 'c("(Intercept)", "hosp", "healthpoor", "healthexcellent", "numchron", "gendermale", "school", "privinsyes")'
TERMS = ["(intercept)", "hospital", "health (poor)", "health (excellent)", "chronic", "gender (male)", "school", "insurance (yes)"]
if HAS_R:
    # robust (HC3) standard errors: nothing printed (the article prints HC0); run sandwich::vcovHC on the article's glm
    r = r_run(f'''d <- data.frame(ofp = visits, hosp = hospital, health = factor(health, levels = c("average", "poor", "excellent")),
                       numchron = chronic, gender = factor(gender), school = school, privins = factor(insurance))
                  m <- glm(ofp ~ ., data = d, family = poisson)
                  list(hc3 = unname(sqrt(diag(sandwich::vcovHC(m, type = "HC3")))[{ORDER}]))''', nmes, packages=("sandwich",))
    refs.run_all("NMES 1988", 'explain (outcome="counts", robust_errors=True)', {f"{t} standard error (HC3)": v for t, v in zip(TERMS, r["hc3"])},
                 "R sandwich vcovHC (type HC3) of the article's Poisson glm", rman("sandwich", "vcovHC"), "sandwich", ITERATIVE)

# 2. bioChemists: run the help-page example of pscl::zeroinfl
bio = r_dataset("bioChemists", "pscl")
refs.data("bioChemists", bio)
BIO_TERMS = ["(intercept)", "fem (Women)", "mar (Single)", "kid5", "phd", "ment"]
if HAS_R:
    r = r_run('''d <- data.frame(art = art, fem = factor(fem, levels = c("Men", "Women")), mar = factor(mar, levels = c("Married", "Single")),
                                 kid5 = kid5, phd = phd, ment = ment)
                 s <- summary(pscl::zeroinfl(art ~ ., data = d))$coefficients
                 list(cb = unname(s$count[, 1]), cs = unname(s$count[, 2]), cp = unname(s$count[, 4]),
                      zb = unname(s$zero[, 1]), zs = unname(s$zero[, 2]), zp = unname(s$zero[, 4]))''', bio, packages=("pscl",))
    for part, key in (("count", "c"), ("extra zeros", "z")):
        for i, term in enumerate(BIO_TERMS):
            for column, k in (("log ratio", "b"), ("standard error", "s"), ("p-value", "p")):
                refs.run("bioChemists (pscl)", 'explain (outcome="counts with many zeros")', f"{part}: {term} {column}", r[key + k][i],
                         "R pscl zeroinfl (help page example)", rman("pscl", "zeroinfl"), "pscl", ITERATIVE)

# 3. Bliss's beetles: Dobson & Barnett (2008), section 7.3 (printed); profile intervals: run R confint (stats)
DB = "Dobson & Barnett (2008), An Introduction to Generalized Linear Models, 3rd ed., Table 7.3"
DB_URL = "https://doi.org/10.1201/9780367807849"
for link, ((a, sa), (b, sb)) in {"probit": (("-34.935", "2.648"), ("19.728", "1.487")),
                                 "cloglog": (("-39.572", "3.240"), ("22.041", "1.799"))}.items():
    check = f'explain (outcome="yes or no", out_of=, link="{link}")'
    refs.printed_all("Bliss beetles", check, {"intercept": a, "intercept standard error": sa, "slope": b, "slope standard error": sb},
                     DB, DB_URL, "published")
bliss = load_csv("bliss.csv")
if HAS_R:
    r = r_run('''m <- glm(cbind(killed, beetles - killed) ~ log10_dose, family = binomial); ci <- suppressMessages(confint(m))
                 list(a = ci[1, 1], b = ci[1, 2], c = ci[2, 1], d = ci[2, 2])''', bliss)
    refs.run_all("Bliss beetles", 'explain (outcome="yes or no", out_of=, interval_by="profile")',
                 {"intercept 95 % low": r["a"], "intercept 95 % high": r["b"], "slope 95 % low": r["c"], "slope 95 % high": r["d"]},
                 "R confint of the glm (profile likelihood)", rman("MASS", "confint-MASS"), "MASS", ITERATIVE)

# 4. Gasoline yields: Ferrari & Cribari-Neto (2004) Table 1, printed in the betareg vignette (batch 10 is R's baseline)
FC = "Ferrari & Cribari-Neto (2004) Table 1, printed in the betareg vignette (Cribari-Neto & Zeileis 2010)"
FC_URL = "https://doi.org/10.18637/jss.v034.i02 ; " + VIGNETTE.format(package="betareg", file="betareg.html")
betareg = ("betareg", r_version("betareg"))
values = {"intercept (batch 10 coding)": "-6.159571", "batch 1 (batch 10 coding)": "1.727729", "batch 1 standard error": "0.101229",
          "temp": "0.010967", "temp standard error": "0.000413"}
values.update({f"batch {k} (batch 10 coding)": p for k, p in zip(range(2, 10), ("1.322597", "1.572310", "1.059714", "1.133752", "1.040162",
                                                                                 "0.543692", "0.495901", "0.385793"))})
refs.printed_all("GasolineYield (betareg)", 'explain (outcome="shares")', values, FC, FC_URL, "published", *betareg)
refs.data("GasolineYield", r_dataset("GasolineYield", "betareg"))

# 5. Clotting times: McCullagh & Nelder (1989), section 8.4.2 (printed, inverse link); log link: run the R glm help page
MN = "McCullagh & Nelder (1989), Generalized Linear Models, 2nd ed., section 8.4.2"
refs.printed_all("clotting times (R glm help page)", 'explain (outcome="positive numbers", link="inverse")',
                 {"intercept": "-0.01655", "slope": "0.01534"}, MN, "https://doi.org/10.1201/9780203753736", "published")
if HAS_R:
    r = r_run('''u <- c(5, 10, 15, 20, 30, 40, 60, 80, 100); lot1 <- c(118, 58, 42, 35, 27, 25, 21, 19, 18)
                 s <- summary(glm(lot1 ~ log(u), family = Gamma(link = "log")))$coefficients
                 list(intercept = s[1, 1], slope = s[2, 1], `intercept standard error` = s[1, 2], `slope standard error` = s[2, 2],
                      `slope p-value` = s[2, 4])''', None)
    refs.run_all("clotting times (R glm help page)", 'explain (outcome="positive numbers")', r,
                 "R stats glm (help page clotting example, Gamma with log link)", rman("stats", "glm"), "stats", ITERATIVE)

# 6. Moore: run the help-page example of car::Anova (types 2 and 3) and anova() (in order), with effectsize
moore = r_dataset("Moore", "carData")
refs.data("Moore", moore)
EFFECTS = ["fcategory", "partner.status", "fcategory : partner.status"]
if HAS_R:
    for kind, call in (("type 2", "car::Anova(m, type = 2)"), ("type 3", "car::Anova(m, type = 3)"), ("in order", "anova(m)")):
        r = r_run(f'''d <- data.frame(conformity = conformity, fcategory = factor(fcategory), partner.status = factor(partner.status))
                      m <- lm(conformity ~ fcategory * partner.status, data = d,
                              contrasts = list(fcategory = contr.sum, partner.status = contr.sum))
                      a <- {call}; k <- if (nrow(a) == 5) 2:4 else 1:3
                      e <- effectsize::eta_squared(a, partial = TRUE)
                      list(ss = a[k, "Sum Sq"], f = a[k, "F value"], df = a[k, "Df"], p = a[k, "Pr(>F)"],
                           eta = e$Eta2_partial[e$Parameter != "(Intercept)"])''', moore, packages=("car",))
        values = {}
        for i, e in enumerate(EFFECTS):
            values.update({f"{e}: sum of squares": r["ss"][i], f"{e}: F": r["f"][i], f"{e}: degrees of freedom": r["df"][i],
                           f"{e}: p-value": r["p"][i]})
        url = rman("car", "Anova") if kind != "in order" else rman("stats", "anova")
        refs.run_all("Moore (carData)", f'explain (table_of_tests="{kind}")', values, f"R {call.split('(')[0]} (help page example)",
                     url, "car" if kind != "in order" else "stats", CLOSED)
        refs.run_all("Moore (carData)", f'explain (table_of_tests="{kind}")',
                     {f"{e}: partial eta-squared": r["eta"][i] for i, e in enumerate(EFFECTS)}, "R effectsize eta_squared (partial)",
                     rman("effectsize", "eta_squared"), "effectsize", CLOSED)

# 7. housing: run the help-page example of MASS::polr (weights = Freq) and nnet::multinom
housing = r_dataset("housing", "MASS")
refs.data("housing_one_row_per_tenant", housing.loc[housing.index.repeat(housing["Freq"].astype(int))].reset_index(drop=True))
PO = ["Infl (Low)", "Infl (Medium)", "Type (Atrium)", "Type (Terrace)", "Type (Tower)", "Cont (Low)"]
R_HOUSE = '''d <- data.frame(Sat = factor(Sat, levels = c("Low", "Medium", "High"), ordered = TRUE), Infl = factor(Infl, levels = c("High", "Low", "Medium")),
                             Type = factor(Type, levels = c("Apartment", "Atrium", "Terrace", "Tower")), Cont = factor(Cont, levels = c("High", "Low")), Freq = Freq)\n'''
if HAS_R:
    r = r_run(R_HOUSE + '''m <- MASS::polr(Sat ~ Infl + Type + Cont, weights = Freq, data = d, Hess = TRUE)
                 list(b = unname(coef(m)), se = unname(sqrt(diag(vcov(m)))[1:6]), z = unname(m$zeta))''', housing, packages=("MASS",))
    values = {f"{t} estimate": v for t, v in zip(PO, r["b"])}
    values.update({f"{t} standard error": v for t, v in zip(PO, r["se"])})
    values.update({"cut-point Low | Medium": r["z"][0], "cut-point Medium | High": r["z"][1]})
    refs.run_all("housing (MASS)", 'explain (outcome="ordered groups", in_order=)', values, "R MASS polr (help page example)",
                 rman("MASS", "polr"), "MASS", ITERATIVE)
    r = r_run(R_HOUSE.replace('levels = c("Low", "Medium", "High"), ordered = TRUE', 'levels = c("High", "Low", "Medium")') + '''
                 m <- nnet::multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = d, trace = FALSE, Hess = TRUE)
                 list(b = as.vector(t(coef(m))), se = as.vector(t(summary(m)$standard.errors)))''', housing, packages=("nnet",))
    values = {}
    for j, g in enumerate(("Low", "Medium")):
        values.update({f"{g}: {t} estimate": r["b"][7 * j + i] for i, t in enumerate(["(intercept)"] + PO)})
        values.update({f"{g}: {t} standard error": r["se"][7 * j + i] for i, t in enumerate(["(intercept)"] + PO)})
    refs.run_all("housing (MASS)", 'explain (outcome="groups")', values, "R nnet multinom", rman("nnet", "multinom"), "nnet", ITERATIVE)

# 8. Engel: the quantreg vignette prints the median regression (printed); other percentiles: run rq
QR = "Koenker, Quantile Regression in R: a vignette (quantreg), sections 3 and 5"
QR_URL = VIGNETTE.format(package="quantreg", file="rq.pdf")
quantreg = ("quantreg", r_version("quantreg"))
refs.printed_all("engel (quantreg)", "explain (for_percentile=)",
                 {"percentile 50: (intercept)": "81.4822474", "percentile 50: income": "0.5601806",
                  "percentile 50: (intercept) standard error": "19.25066", "percentile 50: income standard error": "0.02828"},
                 QR, QR_URL, "printed", *quantreg)
engel = r_dataset("engel", "quantreg")
refs.data("engel", engel)
if HAS_R:
    r = r_run('''m <- quantreg::rq(foodexp ~ income, tau = c(0.05, 0.25, 0.75, 0.95))
                 s <- lapply(summary(m, se = "nid"), function(x) coef(x))
                 list(b = as.vector(coef(m)), se = unlist(lapply(s, function(x) x[, 2])))''', engel, packages=("quantreg",))
    for i, p in enumerate((5, 25, 75, 95)):
        for j, term in enumerate(["(intercept)", "income"]):
            refs.run("engel (quantreg)", "explain (for_percentile=)", f"percentile {p}: {term}", r["b"][2 * i + j],
                     "R quantreg rq (help page example)", rman("quantreg", "rq"), "quantreg", ITERATIVE)
            refs.run("engel (quantreg)", "explain (for_percentile=)", f"percentile {p}: {term} standard error", r["se"][2 * i + j],
                     "R quantreg summary.rq (se = nid)", rman("quantreg", "summary.rq"), "quantreg", ITERATIVE)

# 9. stack loss: the statsmodels RLM docstring (printed)
RLM_URL = f"https://www.statsmodels.org/v{py_version('statsmodels')}/generated/statsmodels.robust.robust_linear_model.RLM.html"
for term, (b, s) in {"air_flow": ("0.82938433", "0.11100521"), "water_temp": ("0.92606597", "0.30293016"),
                     "acid_conc": ("-0.12784672", "0.12864961"), "(intercept)": ("-41.02649835", "9.79189854")}.items():
    refs.printed("stackloss", "explain (giving_outliers_less_weight=True)", f"{term} estimate", b, "statsmodels RLM docstring (HuberT)",
                 RLM_URL, "printed", "statsmodels", py_version("statsmodels"))
    refs.printed("stackloss", "explain (giving_outliers_less_weight=True)", f"{term} standard error", s, "statsmodels RLM docstring (HuberT)",
                 RLM_URL, "printed", "statsmodels", py_version("statsmodels"))

# 10. standardised effects (run effectsize) and precision weights (run lm)
stackloss = load_csv("stackloss.csv")
if HAS_R:
    r = r_run('''s <- effectsize::standardize_parameters(lm(stack_loss ~ air_flow + water_temp + acid_conc), method = "refit")
                 list(b = s$Std_Coefficient[-1], lo = s$CI_low[-1])''', stackloss, packages=("effectsize",))
    COLS = ["air_flow", "water_temp", "acid_conc"]
    refs.run_all("stackloss", 'explain (show="standardised effects")',
                 {**{f"{t} estimate": v for t, v in zip(COLS, r["b"])}, **{f"{t} 95 % low": v for t, v in zip(COLS, r["lo"])}},
                 "R effectsize standardize_parameters (refit)", rman("effectsize", "standardize_parameters"), "effectsize", CLOSED)
    r = r_run('''s <- summary(lm(I(killed / beetles) ~ log10_dose, weights = beetles))$coefficients
                 list(intercept = s[1, 1], slope = s[2, 1], `slope standard error` = s[2, 2], `slope p-value` = s[2, 4])''', bliss)
    refs.run_all("Bliss beetles", "explain (weights=)", r, "R stats lm (weights =)", rman("stats", "lm"), "stats", CLOSED)

# 11. predictions: run R predict
birthwt = load_csv("birthwt.csv")
if HAS_R:
    for new_rows in (False, True):
        r = r_run(f'''p <- predict(lm(stack_loss ~ air_flow + water_temp), interval = "{'prediction' if new_rows else 'confidence'}")
                      list(predicted = p[, 1], low = p[, 2], high = p[, 3])''', stackloss)
        check = f"add_predicted_column (interval=95{', for_new_rows=True' if new_rows else ''})"
        for key in ("predicted", "low", "high"):
            refs.run_vector("stackloss", check, key, r[key], "R stats predict.lm", rman("stats", "predict.lm"), "stats", CLOSED)
    r = r_run('''m <- glm(low ~ age + lwt + smoke, family = binomial); p <- predict(m, se.fit = TRUE)
                 list(predicted = plogis(p$fit), low = plogis(p$fit - qnorm(0.975) * p$se.fit), high = plogis(p$fit + qnorm(0.975) * p$se.fit))''',
              birthwt)
    for key in ("predicted", "low", "high"):
        refs.run_vector("birthwt", 'add_predicted_column (outcome="yes or no", interval=95)', key, r[key],
                        "R stats predict.glm (link scale, turned back)", rman("stats", "predict.glm"), "stats", ITERATIVE)
    r = r_run('''m <- glm(visits ~ hospital + chronic, family = poisson); p <- predict(m, se.fit = TRUE)
                 list(predicted = exp(p$fit), low = exp(p$fit - qnorm(0.975) * p$se.fit), high = exp(p$fit + qnorm(0.975) * p$se.fit))''', nmes)
    for key in ("predicted", "low", "high"):
        refs.run_vector("NMES 1988", 'add_predicted_column (outcome="counts", interval=95)', key, r[key],
                        "R stats predict.glm (link scale, turned back)", rman("stats", "predict.glm"), "stats", ITERATIVE, every=100)
    r = r_run('''s <- killed / beetles; m <- nls(s ~ a + b * log10_dose, start = list(a = 0, b = 1), weights = beetles)
                 list(predicted = unname(predict(m)))''', bliss)
    refs.run_vector("Bliss beetles", "add_predicted_column (formula, weights=)", "predicted", r["predicted"], "R stats predict(nls(weights =))",
                    rman("stats", "nls"), "stats", ITERATIVE)

# 12. rows that pull a model and the numbers of check_model_of: nothing printed, so R is run once. The model of
# LifeCycleSavings is the example of R's influence.measures help page; the logistic model of birthwt is section 11's.
# show_rows_that_pull_model gives back the rows with a Cook's distance above 4 / n, the biggest first.
lcs = r_dataset("LifeCycleSavings")
refs.data("LifeCycleSavings", lcs)
if HAS_R:
    PULL = """keep <- which(cd > 4 / length(cd)); keep <- keep[order(-cd[keep])]
              list(row = keep, cooks = unname(cd[keep]), hat = unname(h[keep]), residual = unname(res[keep]))"""
    r = r_run("m <- lm(sr ~ pop15 + pop75 + dpi + ddpi); cd <- cooks.distance(m); h <- hatvalues(m); res <- rstudent(m); " + PULL, lcs)
    for key, name in (("row", "row"), ("cooks", "Cook's distance"), ("hat", "leverage"), ("residual", "studentised residual")):
        refs.run_vector("LifeCycleSavings", 'show_rows_that_pull_model (outcome="numbers")', name, np.atleast_1d(r[key]),
                        "R stats cooks.distance, hatvalues, rstudent", rman("stats", "influence.measures"), "stats",
                        0 if key == "row" else CLOSED, "absolute" if key == "row" else "relative")
    r = r_run('m <- glm(low ~ age + lwt + smoke, family = binomial); cd <- cooks.distance(m); h <- hatvalues(m); '
              'res <- rstandard(m, type = "pearson"); ' + PULL, birthwt)
    for key, name in (("row", "row"), ("cooks", "Cook's distance"), ("hat", "leverage"), ("residual", "standardised residual")):
        refs.run_vector("birthwt", 'show_rows_that_pull_model (outcome="yes or no")', name, np.atleast_1d(r[key]),
                        'R stats cooks.distance, hatvalues, rstandard(type = "pearson") of a glm', rman("stats", "influence.measures"),
                        "stats", 0 if key == "row" else ITERATIVE, "absolute" if key == "row" else "relative")
    # check_model_of(show="numbers"): car::vif, lmtest::bptest (Koenker's studentised form, its default) and White's test
    # as bptest on the columns, their squares and products; and the rows with a Cook's distance above 0.5
    r = r_run("""m <- lm(sr ~ pop15 + pop75 + dpi + ddpi); v <- car::vif(m); bp <- lmtest::bptest(m)
                 w <- lmtest::bptest(m, ~ (pop15 + pop75 + dpi + ddpi)^2 + I(pop15^2) + I(pop75^2) + I(dpi^2) + I(ddpi^2))
                 list(vif = unname(v), bp = unname(bp$statistic), bp_df = unname(bp$parameter), bp_p = bp$p.value,
                      w = unname(w$statistic), w_df = unname(w$parameter), w_p = w$p.value, big = sum(cooks.distance(m) > 0.5))""",
              lcs, packages=("car", "lmtest"))
    check = 'check_model_of (show="numbers")'
    refs.run_all("LifeCycleSavings", check, {f"variance inflation factor, {c}": v for c, v in zip(["pop15", "pop75", "dpi", "ddpi"], r["vif"])},
                 "R car vif", rman("car", "vif"), "car", CLOSED)
    for name, key in (("Breusch-Pagan test (chi-square)", "bp"), ("White's test (chi-square)", "w")):
        refs.run_all("LifeCycleSavings", check, {name: r[key], f"{name} df": r[f"{key}_df"], f"{name} p-value": r[f"{key}_p"]},
                     "R lmtest bptest (studentised)", rman("lmtest", "bptest"), "lmtest", CLOSED)
    refs.run("LifeCycleSavings", check, "rows with a Cook's distance above 0.5", r["big"], "R stats cooks.distance",
             rman("stats", "influence.measures"), "stats", 0, "absolute")

refs.save()
