"""Reference values for notebooks/15_surveys_and_designs.ipynb, saved to references/15_surveys_and_designs.csv.

- Printed: Lumley (2010, Complex Surveys, chapter 2) prints the stratified mean API 2000 662.29 (standard error
  9.4089), the total enrolment 3,687,178 (114,642) and the cluster mean 644.17 (23.542): half a unit of the last
  printed digit (the design degrees of freedom, 197 and 14, exactly).
- Run here: R survey (svydesign, svymean, svytotal, svyquantile, svyciprop, svyby, svyglm, svyttest, svychisq) with
  the same design (1e-9 relative); the textbook formulas of Lohr (2010, sections 3.2 and 9.1) in NumPy/pandas
  (1e-10 relative); R agricolae and R blockrand with the same seed and generator, and pandas groupby().sample()
  with the same random_state (exact).
- Definitions (latin square, split plot, every combination as often as the design says, permuted blocks, Yates'
  standard order): the count of plots or blocks that break the definition must be 0; the notebook counts them on
  simplydata's table with simplydata words.

A layout is stored as one list of numbers, column after column (a level given as text is stored as its position
in the list of levels the notebook passes, counting from 1), and compared by its largest relative difference
(tolerance 0: every cell the same).

    python scripts/references/15_surveys_and_designs.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "15_surveys_and_designs"
checks = Recorder(STEM, ("survey", "agricolae", "blockrand"))
HAS_R = r_available()
LUMLEY = "Lumley (2010), printed"
LUMLEY_TITLE = "Lumley T (2010). Complex Surveys: A Guide to Analysis Using R. Wiley, chapter 2"
FORMULA = "NumPy, textbook formula (Lohr 2010)"
CLOSED = 1e-10
R_SURVEY = 1e-9

apistrat = load_csv("apistrat.csv")
apistrat["sch"] = (apistrat["sch.wide"] == "Yes").astype(float)
STRAT_R = ('d <- survey::svydesign(ids = ~1, strata = ~stype, weights = ~pw, fpc = ~fpc, '
           'data = data.frame(stype, pw, fpc, api00, enroll, meals, ell, sch, awards))\n')
DS = "api stratified sample"
M = "after use_survey_design (strata), find_average_of / find_total_of"
for quantity, shown in {"mean of api00": "662.29", "standard error of the mean": "9.4089", "total of enroll": "3687178",
                        "standard error of the total": "114642"}.items():
    checks.printed(DS, M, quantity, shown, LUMLEY, title=LUMLEY_TITLE)
checks.compare(DS, M, "degrees of freedom", None, 197, LUMLEY, 0)
checks.rows[-1]["source title"] = LUMLEY_TITLE


def linearised(y: str, data: pd.DataFrame, strata, psu, ratio: bool):
    """Lohr (2010), sections 3.2 and 9.1: the weighted total (or ratio mean) and its linearised standard error."""
    w, values = data["pw"], data[y]
    estimate = (w * values).sum() / (w.sum() if ratio else 1)
    z = w * (values - estimate) / w.sum() if ratio else w * values
    frame = pd.DataFrame({"z": z, "h": strata, "psu": psu, "N": data["fpc"]})
    units = frame.groupby(["h", "psu"]).agg(z=("z", "sum"), N=("N", "first")).reset_index()
    variance = 0.0
    for _, g in units.groupby("h"):
        n_h, N_h = len(g), g["N"].iloc[0]
        variance += (1 - n_h / N_h) * n_h / (n_h - 1) * ((g["z"] - g["z"].mean()) ** 2).sum()
    return estimate, float(np.sqrt(variance))


m, se_m = linearised("api00", apistrat, apistrat["stype"], np.arange(len(apistrat)), ratio=True)
t, se_t = linearised("enroll", apistrat, apistrat["stype"], np.arange(len(apistrat)), ratio=False)
for quantity, value in {"mean of api00": m, "standard error of the mean": se_m, "total of enroll": t,
                        "standard error of the total": se_t}.items():
    checks.compare(DS, M, quantity, None, value, FORMULA, CLOSED, "relative")
if HAS_R:
    r = r_run(STRAT_R + """m <- survey::svymean(~api00, d); cm <- confint(m, df = survey::degf(d))
        t <- survey::svytotal(~enroll, d); ct <- confint(t, df = survey::degf(d))
        list(`mean of api00` = coef(m)[[1]], `standard error of the mean` = survey::SE(m)[[1]],
             `95 % CI low of the mean` = cm[1], `total of enroll` = coef(t)[[1]],
             `standard error of the total` = survey::SE(t)[[1]], `95 % CI high of the total` = ct[2],
             `degrees of freedom` = survey::degf(d))""", apistrat)
    checks.compare_all(DS, M, None, r, "R survey svymean / svytotal", R_SURVEY, "relative")

    M = "after use_survey_design (strata), other words"
    r = r_run(STRAT_R + """q <- survey::svyquantile(~api00, d, quantiles = 0.5, qrule = "math", ci = TRUE, df = survey::degf(d))$api00
        q9 <- survey::svyquantile(~api00, d, quantiles = 0.9, qrule = "math", ci = TRUE, df = survey::degf(d))$api00
        p <- survey::svyciprop(~sch, d, method = "logit"); pci <- attr(p, "ci")
        by <- survey::svyby(~one, ~awards, stats::update(d, one = 1), survey::svytotal)
        sh <- survey::svymean(~awards, d)
        g <- summary(survey::svyglm(api00 ~ ell + meals, design = d))$coefficients
        tt <- survey::svyttest(api00 ~ awards, d)
        ch <- survey::svychisq(~awards + stype, d)
        list(median = q[1, 1], `median 95 % CI low` = q[1, 2], `median 95 % CI high` = q[1, 3],
             `median standard error` = q[1, 4], `percentile 90` = q9[1, 1], `percentile 90 95 % CI low` = q9[1, 2],
             `percent sch.wide Yes` = 100 * coef(p)[[1]], `percent 95 % CI low` = 100 * pci[[1]],
             `percent 95 % CI high` = 100 * pci[[2]], `schools with an award` = coef(by)[[2]],
             `schools with an award, standard error` = survey::SE(by)[[2]], `percent with an award` = 100 * coef(sh)[[2]],
             intercept = g[1, 1], ell = g[2, 1], meals = g[3, 1], `meals standard error` = g[3, 2], `ell p-value` = g[2, 4],
             `difference (awards)` = tt$estimate[[1]], `t (awards)` = tt$statistic[[1]], `p-value (awards)` = tt$p.value,
             `Rao-Scott F` = ch$statistic[[1]], `Rao-Scott p-value` = ch$p.value,
             `Rao-Scott design degrees of freedom` = ch$parameter[[2]])""", apistrat)
    checks.compare_all(DS, M, None, r, "R survey", R_SURVEY, "relative")

# ---- the one-stage cluster sample
apiclus1 = load_csv("apiclus1.csv")
DS, M = "api cluster sample", "after use_survey_design (clusters), find_average_of"
for quantity, shown in {"mean of api00": "644.17", "standard error of the mean": "23.542"}.items():
    checks.printed(DS, M, quantity, shown, LUMLEY, title=LUMLEY_TITLE)
checks.compare(DS, M, "degrees of freedom", None, 14, LUMLEY, 0)
checks.rows[-1]["source title"] = LUMLEY_TITLE
ybar, se = linearised("api00", apiclus1, np.zeros(len(apiclus1)), apiclus1["dnum"], ratio=True)
checks.compare(DS, M, "mean of api00", None, ybar, FORMULA, CLOSED, "relative")
checks.compare(DS, M, "standard error of the mean", None, se, FORMULA, CLOSED, "relative")
if HAS_R:
    r = r_run('d <- survey::svydesign(ids = ~dnum, weights = ~pw, fpc = ~fpc, data = data.frame(dnum, pw, fpc, api00))\n'
              'm <- survey::svymean(~api00, d)\n'
              'list(`mean of api00` = coef(m)[[1]], `standard error of the mean` = survey::SE(m)[[1]], '
              '`degrees of freedom` = survey::degf(d))', apiclus1)
    checks.compare_all(DS, M, None, r, "R survey svymean", R_SURVEY, "relative")


# ---- design tables: R agricolae with the same seed; definitions
def codes(values, levels=None):
    """A column as numbers: a text level becomes its position in `levels` (counting from 1)."""
    return [float(levels.index(v) + 1) if levels else float(v) for v in values]


DESIGNS = [  # (title, kind, factors, R code giving one vector per compared column, {column: levels})
    ("2 x 2 x 2 factorial, completely random, seed 1", "completely random", ["a", "b", "c"],
     'b <- agricolae::design.ab(c(2, 2, 2), r = 2, design = "crd", serie = 0, seed = 1, kinds = "Mersenne-Twister")$book\n'
     'list(a = as.numeric(as.character(b$A)) - 1, b = as.numeric(as.character(b$B)) - 1, c = as.numeric(as.character(b$C)) - 1)',
     {}),
    ("2 x 2 factorial in 4 randomised blocks, seed 3", "randomised blocks", ["a", "b"],
     'b <- agricolae::design.ab(c(2, 2), r = 4, design = "rcbd", serie = 0, seed = 3, kinds = "Mersenne-Twister")$book\n'
     'list(block = as.numeric(as.character(b$block)), a = as.numeric(as.character(b$A)) - 1, b = as.numeric(as.character(b$B)) - 1)',
     {}),
    ("one factor (3 varieties) in 4 randomised blocks, seed 5", "randomised blocks", ["variety"],
     'b <- agricolae::design.rcbd(c("A", "B", "C"), r = 4, serie = 0, seed = 5, kinds = "Mersenne-Twister")$book\n'
     'list(block = as.numeric(as.character(b$block)), variety = as.character(b[[3]]))', {"variety": ["A", "B", "C"]}),
    ("one factor (4 levels), completely random, 3 repeats, seed 6", "completely random", ["v"],
     'b <- agricolae::design.crd(c("p", "q", "r", "s"), r = 3, serie = 0, seed = 6, kinds = "Mersenne-Twister")$book\n'
     'list(v = as.character(b[[3]]))', {"v": ["p", "q", "r", "s"]}),
    ("latin square of 4 levels, seed 2", "latin square", ["v"],
     'b <- agricolae::design.lsd(c("A", "B", "C", "D"), serie = 0, seed = 2, kinds = "Mersenne-Twister")$book\n'
     'list(row = as.numeric(as.character(b$row)), column = as.numeric(as.character(b$col)), v = as.character(b[[4]]))',
     {"v": ["A", "B", "C", "D"]}),
    ("split plot, irrigation on whole plots, 3 blocks, seed 4", "split plot", ["irrigation", "variety"],
     'b <- agricolae::design.split(c("dry", "wet"), c("A", "B", "C"), r = 3, design = "rcbd", serie = 0, seed = 4, '
     'kinds = "Mersenne-Twister")$book\n'
     'list(block = as.numeric(as.character(b$block)), irrigation = as.character(b[[4]]), variety = as.character(b[[5]]))',
     {"irrigation": ["dry", "wet"], "variety": ["A", "B", "C"]}),
]
for title, kind, factors, r_code, levels in DESIGNS:
    if HAS_R:
        r = r_run(r_code)
        layout = [v for column, values in r.items() for v in codes(values if isinstance(values, list) else [values],
                                                                      levels.get(column))]
        checks.vector(title, "make_design_table", "agricolae's layout, cell by cell", layout,
                      "R agricolae (same seed, Mersenne-Twister)", 0.0)
    if kind == "latin square":
        checks.compare(title, "make_design_table", "rows missing a level", None, 0.0, "definition of a latin square", 0)
        checks.compare(title, "make_design_table", "columns missing a level", None, 0.0, "definition of a latin square", 0)
    elif kind == "split plot":
        checks.compare(title, "make_design_table", "whole plots with more than one irrigation", None, 0.0,
                       "definition of a split plot", 0)
        checks.compare(title, "make_design_table", "whole plots missing a variety", None, 0.0, "definition of a split plot", 0)
    else:
        checks.compare(title, "make_design_table", "combinations not used as often as the design says", None, 0.0,
                       "definition of the design", 0)
        if len(factors) > 1:
            checks.compare(title, "make_design_table", "plots with a wrong standard order", None, 0.0,
                           "Yates' standard order (Montgomery 2017)", 0)

# ---- randomisation lists: R blockrand with the same seed
for how_many, sizes, seed, groups in ((60, [4, 6], 1, ["A", "B"]), (45, [3, 6], 7, ["A", "B", "C"])):
    title = f"{how_many} places, {len(groups)} groups, blocks of {sizes}, two centres, seed {seed}"
    checks.compare(title, "make_random_assignments", "blocks where the groups are not equally often", None, 0.0,
                   "definition of permuted blocks", 0)
    if HAS_R:
        k = len(groups)
        r = r_run(f'''set.seed({seed}, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
            out <- NULL
            for (s in c("Bergen", "Oslo")) {{
              b <- blockrand::blockrand(n = {how_many}, num.levels = {k}, levels = c({", ".join(repr(g) for g in groups)}),
                                        block.sizes = c({", ".join(str(s // k) for s in sizes)}), stratum = s)
              out <- rbind(out, data.frame(list = s, block = as.numeric(as.character(b$block.id)),
                                           size = as.numeric(b$block.size), group = as.character(b$treatment)))
            }}
            list(list = out$list, block = out$block, size = out$size, group = out$group)''')
        layout = codes(r["list"], ["Bergen", "Oslo"]) + codes(r["block"]) + codes(r["size"]) + codes(r["group"], groups)
        checks.vector(title, "make_random_assignments", "blockrand's list, cell by cell", layout,
                      "R blockrand (same seed)", 0.0)

# ---- random rows within groups: pandas groupby().sample() with the same seed
for call, parent in ((dict(how_many=10, seed=2), dict(n=10, random_state=2)),
                     (dict(percent=10, seed=2), dict(frac=0.10, random_state=2)),
                     (dict(how_many=5, seed=11), dict(n=5, random_state=11))):
    if "percent" in call:
        title = f"api stratified sample, {call['percent']} % of each school type, seed {call['seed']}"
    else:
        title = f"api stratified sample, {call['how_many']} of each school type, seed {call['seed']}"
    reference = sorted(apistrat.groupby("stype").sample(**parent)["snum"])
    checks.vector(title, "keep_random_rows", "school numbers kept, in order", reference, "pandas groupby().sample()", 0.0)
    checks.compare(title, "keep_random_rows", "rows kept", None, float(len(reference)), "pandas groupby().sample()", 0)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
