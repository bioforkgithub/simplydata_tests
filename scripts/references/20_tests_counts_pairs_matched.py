"""Reference values for notebooks/20_tests_counts_pairs_matched.ipynb, saved to references/20_tests_counts_pairs_matched.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/20_tests_counts_pairs_matched/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/20_tests_counts_pairs_matched.py      (from the repository root; needs R with the packages listed)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks and realdata)

import re

import numpy as np
import pandas as pd
import simplydata as sdata
from scipy import stats

import realdata as rd
from checks import Checks, load_csv, quiet, r_available, r_run, read_numbers, read_row, read_rows, versions

versions("pandas", "numpy", "scipy", "statsmodels", "rpy2")
checks = Recorder("20_tests_counts_pairs_matched")
HAS_R = r_available()
CLOSED = 1e-10
EXACT = 1e-12
ITERATIVE = 1e-6
DOCSTRING = CLOSED
if HAS_R:
    print(rd.r_version_of("DescTools", "Exact", "PMCMRplus", "FSA", "rstatix", "multcompView", "effectsize",
                          "cocor", "qvalue", "crank", "multcomp"))


def one(word, *args, **options):
    return read_row(quiet(word, *args, **options))


def compare_many(dataset, method, found, reference, reference_name, tolerance=CLOSED, how="relative"):
    for quantity, value in reference.items():
        checks.compare(dataset, method, quantity, found[quantity], value, reference_name, tolerance, how)


def pair_row(rows, a, b):
    """The row of a pairs table (read with checks.read_rows) for groups a and b, in either order."""
    return next(r for r in rows if {r["group 1"], r["group 2"]} == {a, b})


def barnard_by_grid(table, points=20001):
    """Barnard's unconditional test (pooled z, two-sided) for a 2 x 2 table whose ROWS are the two samples,
    written out in NumPy: every table with the same row totals, and the largest p-value over a grid of the
    common chance of the event (an independent implementation, not SciPy's or R Exact's search)."""
    (a, b), (c, d) = np.asarray(table, dtype=int)
    n1, n2 = a + b, c + d
    x1, x2 = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
    p1, p2 = x1 / n1, x2 / n2
    pooled = (x1 + x2) / (n1 + n2)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (p1 - p2) / np.sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
    z = np.nan_to_num(z, nan=0.0)
    extreme = np.abs(z) >= abs(z[a, c]) - 1e-12
    best = 0.0
    for pi in np.linspace(0, 1, points)[1:-1]:
        prob = stats.binom.pmf(x1, n1, pi) * stats.binom.pmf(x2, n2, pi)
        best = max(best, prob[extreme].sum())
    return best


def rows_from_table(counts, row_name, column_name, row_labels, column_labels):
    """A table of counts (rows x columns) turned into one row per person."""
    out = [(r, c) for i, r in enumerate(row_labels) for j, c in enumerate(column_labels)
           for _ in range(int(counts[i][j]))]
    return pd.DataFrame(out, columns=[row_name, column_name])


made = rd.r_example_objects("stats", "fisher.test")
print("made by the help page:", made)
tables = {}
for name in ("TeaTasting", "Convictions", "Job"):
    t = r_run(f"list(n = as.vector({name}), r = rownames({name}), c = colnames({name}), nr = nrow({name}))")
    nr = int(t["nr"])
    counts = np.array(t["n"]).reshape(-1, nr).T   # R stores column by column
    tables[name] = (counts, t["r"], t["c"])
    print(name, "\n", pd.DataFrame(counts, index=t["r"], columns=t["c"]))

for name in ("TeaTasting", "Convictions"):
    counts, row_labels, column_labels = tables[name]
    frame = rows_from_table(counts, "row", "column", row_labels, column_labels)
    checks.data(f"{name}_rows", frame)
    t = sdata.table(frame)
    ds = f"{name} (R fisher.test example)"
    found = one(t.compare_counts_of, "column", for_each_group_in="row", test_method="fisher")
    if HAS_R:
        r = r_run(f"""f <- fisher.test({name}); list(`p-value` = f$p.value, or = unname(f$estimate))""")
        checks.compare(ds, 'compare_counts_of (test_method="fisher")', "p-value", found["p-value"], r["p-value"],
                       "R fisher.test", CLOSED, "relative")
    via_counts = one(sdata.compare_counts, counts.tolist(), test_method="fisher")
    checks.compare(ds, 'sdata.compare_counts (test_method="fisher")', "p-value", via_counts["p-value"],
                   found["p-value"], "table.compare_counts_of on the rows", EXACT, "relative")
    # Barnard: the two rows (the groups of for_each_group_in=) are the two samples
    reference = barnard_by_grid(counts)
    for word, call in [("compare_counts_of", lambda: one(t.compare_counts_of, "column", for_each_group_in="row",
                                                         test_method="barnard")),
                       ("sdata.compare_counts", lambda: one(sdata.compare_counts, counts.tolist(), test_method="barnard"))]:
        checks.compare(ds, f'{word} (test_method="barnard")', "p-value", call()["p-value"], reference,
                       "Barnard written out, rows as samples (NumPy grid of 20 001)", 1e-4, "relative")
    for test, method in [("barnard", "z-pooled"), ("boschloo", "boschloo")]:
        if not HAS_R:
            break
        r = r_run(f"""suppressPackageStartupMessages(library(Exact))
            list(rows = exact.test({name}, method = "{method}", to.plot = FALSE)$p.value,
                 columns = exact.test(t({name}), method = "{method}", to.plot = FALSE)$p.value)""")
        found = one(t.compare_counts_of, "column", for_each_group_in="row", test_method=test)
        via_counts = one(sdata.compare_counts, counts.tolist(), test_method=test)
        checks.compare(ds, f'compare_counts_of (test_method="{test}")', "p-value", found["p-value"], r["rows"],
                       "R Exact exact.test, groups as the samples", ITERATIVE, "relative")
        checks.compare(ds, f'sdata.compare_counts (test_method="{test}")', "p-value", via_counts["p-value"], r["rows"],
                       "R Exact exact.test, rows as the samples", ITERATIVE, "relative")
    via_counts = one(sdata.compare_counts, counts.tolist())
    if HAS_R:
        r = r_run(f"""x <- suppressWarnings(chisq.test({name})); list(statistic = unname(x$statistic), `p-value` = x$p.value)""")
        compare_many(ds, "sdata.compare_counts", via_counts, r, "R chisq.test")


counts, row_labels, column_labels = tables["Job"]
job = rows_from_table(counts, "income", "satisfaction", row_labels, column_labels)
jobs = sdata.table(job)
ds = "job satisfaction (R fisher.test example; Agresti 2002)"
found = one(jobs.compare_counts_of, "satisfaction", for_each_group_in="income", test_method="fisher", with_seed=1)
print(found)
for orientation, tab in [("groups as rows", pd.crosstab(job["income"], job["satisfaction"]).to_numpy()),
                         ("groups as columns", pd.crosstab(job["satisfaction"], job["income"]).to_numpy())]:
    parent = stats.fisher_exact(tab, method=stats.MonteCarloMethod(n_resamples=100000, rng=np.random.default_rng(1)))
    print(orientation, parent.pvalue)
    if abs(parent.pvalue - found["p-value"]) < 1e-12:
        break
checks.compare(ds, 'compare_counts_of (test_method="fisher", 4 x 4, with_seed=1)', "p-value", found["p-value"],
               parent.pvalue, f"SciPy fisher_exact (Monte Carlo), same seed ({orientation}, sorted labels)", EXACT,
               "relative")
if HAS_R:
    r = r_run("""f <- fisher.test(Job); g <- DescTools::GTest(Job); x <- suppressWarnings(chisq.test(Job))
        set.seed(1); s <- chisq.test(Job, simulate.p.value = TRUE, B = 200000)
        list(fp = f$p.value, gs = unname(g$statistic), gp = g$p.value, stdres = as.vector(t(x$stdres)), sp = s$p.value)""")
    p = r["fp"]
    checks.compare(ds, 'compare_counts_of (test_method="fisher", 4 x 4, with_seed=1)', "p-value", found["p-value"], p,
                   "R fisher.test exact (4 MCSE)", 4 * np.sqrt(p * (1 - p) / 100000))
    found = one(jobs.compare_counts_of, "satisfaction", for_each_group_in="income", test_method="fisher",
                p_value_by="exact")
    checks.compare(ds, 'compare_counts_of (test_method="fisher", p_value_by="exact")', "p-value", found["p-value"], p,
                   "R fisher.test", CLOSED, "relative")
    found = one(jobs.compare_counts_of, "satisfaction", for_each_group_in="income", test_method="g")
    compare_many(ds, 'compare_counts_of (test_method="g")', found, {"statistic": r["gs"], "p-value": r["gp"]},
                 "R DescTools GTest")
    cells = {(row["income"], row["satisfaction"]): row
             for row in read_rows(quiet(jobs.compare_counts_of, "satisfaction", for_each_group_in="income", show="cells"))}
    for k, (i, j) in enumerate([(i, j) for i in range(4) for j in range(4)]):
        checks.compare(ds, 'compare_counts_of (show="cells")', f"standardised residual {row_labels[i]}, {column_labels[j]}",
                       cells[(row_labels[i], column_labels[j])]["standardised residual"], r["stdres"][k],
                       "R chisq.test()$stdres", CLOSED, "relative")
    found = one(jobs.compare_counts_of, "satisfaction", for_each_group_in="income", p_value_by="shuffling",
                shuffles=9999, with_seed=1)
    p = r["sp"]
    checks.compare(ds, 'compare_counts_of (p_value_by="shuffling", shuffles=9999, with_seed=1)', "p-value",
                   found["p-value"], p, "R chisq.test(simulate.p.value, B = 200000) (4 MCSE)",
                   4 * np.sqrt(p * (1 - p) / 9999 + p * (1 - p) / 200000))


made = rd.r_example_objects("stats", "prop.trend.test")
fl = r_run("list(s = smokers, n = patients)")
groups = ["g1", "g2", "g3", "g4"]
smoke = pd.DataFrame([(g, "smoker") for g, s in zip(groups, fl["s"]) for _ in range(int(s))]
                     + [(g, "non-smoker") for g, s, n in zip(groups, fl["s"], fl["n"]) for _ in range(int(n - s))],
                     columns=["age group", "smoking"])
found = one(sdata.table(smoke).compare_counts_of, "smoking", for_each_group_in="age group", test_method="trend",
            in_order=groups)
x, n = np.array(fl["s"]), np.array(fl["n"])
scores = np.arange(1, 5)
pbar = x.sum() / n.sum()
tbar = (n * scores).sum() / n.sum()
chi = ((x * (scores - tbar)).sum()) ** 2 / (pbar * (1 - pbar) * (n * (scores - tbar) ** 2).sum())
compare_many("smokers (R prop.trend.test example; Fleiss 1981)", 'compare_counts_of (test_method="trend", in_order=)',
             found, {"statistic": chi, "p-value": stats.chi2.sf(chi, 1)}, "Cochran-Armitage written out, NumPy")
if HAS_R:
    r = r_run("p <- prop.trend.test(smokers, patients); list(statistic = unname(p$statistic), `p-value` = p$p.value)")
    compare_many("smokers (R prop.trend.test example; Fleiss 1981)", 'compare_counts_of (test_method="trend", in_order=)',
                 found, r, "R prop.trend.test")


hec = rd.r_dataset("HairEyeColor")
eyes = pd.DataFrame({"eye": np.repeat(hec["Eye"].to_numpy(), hec["Freq"].astype(int).to_numpy())})
students = sdata.table(eyes)
ds = "HairEyeColor eye colours (Snee 1974)"
order = ["Brown", "Blue", "Hazel", "Green"]
observed = eyes["eye"].value_counts().reindex(order).to_numpy(dtype=float)
for expected in (None, {"Brown": 0.4, "Blue": 0.3, "Hazel": 0.15, "Green": 0.15}):
    found = one(students.compare_counts_to, "eye", expected=expected)
    label = "compare_counts_to" + (" (expected=)" if expected else "")
    if HAS_R:
        p = "rep(0.25, 4)" if expected is None else "c(0.4, 0.3, 0.15, 0.15)"
        r = r_run(f"""x <- chisq.test(o, p = {p}); list(statistic = unname(x$statistic), `p-value` = x$p.value,
            `degrees of freedom` = unname(x$parameter), `effect size` = sqrt(unname(x$statistic) / sum(o)))""",
                  {"o": observed})
        compare_many(ds, label, found, r, "R chisq.test")


example = rd.r_example_code("stats", "binom.test")
giant, dwarf = map(int, re.search(r"binom\.test\(c\((\d+),\s*(\d+)\),\s*p\s*=\s*3/4\)", example).groups())
print("from the help page:", giant, "giant,", dwarf, "dwarf")
pea_rows = pd.DataFrame({"plant": ["giant"] * giant + ["dwarf"] * dwarf})
peas = sdata.table(pea_rows)
ds = "Mendel's peas (R binom.test example; Conover 1971)"
if HAS_R:
    for opts, alternative in [({}, "two.sided"), ({"sides": 1, "higher": "percent"}, "greater"),
                              ({"sides": 1, "lower": "percent"}, "less"), ({"alternative": "less"}, "less")]:
        found = one(peas.compare_percent_to, "plant", value="giant", percent=75, **opts)
        r = r_run(f"""b <- binom.test(c({giant}, {dwarf}), p = 3/4, alternative = "{alternative}")
            list(`p-value` = b$p.value, low = 100 * b$conf.int[1], high = 100 * b$conf.int[2])""")
        label = "compare_percent_to" + (" (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")" if opts else "")
        checks.compare(ds, label, "p-value", found["p-value"], r["p-value"], "R binom.test", CLOSED, "relative")
        if not opts:
            checks.compare(ds, label, "percent 95 % CI low", found["percent 95 % CI low"], r["low"], "R binom.test",
                           CLOSED, "relative")
            checks.compare(ds, label, "percent 95 % CI high", found["percent 95 % CI high"], r["high"], "R binom.test",
                           CLOSED, "relative")
            h = 2 * np.arcsin(np.sqrt(giant / (giant + dwarf))) - 2 * np.arcsin(np.sqrt(0.75))
            checks.compare(ds, label, "effect size (Cohen's h)", found["effect size"], h, "formula, NumPy", CLOSED,
                           "relative")
    for method in ("wilson", "agresti-coull"):
        found = one(peas.compare_percent_to, "plant", value="giant", percent=75, interval_by=method)
        r = r_run(f"""b <- DescTools::BinomCI({giant}, {giant + dwarf}, method = "{method}")
            list(low = 100 * b[2], high = 100 * b[3])""")
        checks.compare(ds, f'compare_percent_to (interval_by="{method}")', "percent 95 % CI low",
                       found["percent 95 % CI low"], r["low"], "R DescTools BinomCI", CLOSED, "relative")
        checks.compare(ds, f'compare_percent_to (interval_by="{method}")', "percent 95 % CI high",
                       found["percent 95 % CI high"], r["high"], "R DescTools BinomCI", CLOSED, "relative")


tk = pd.DataFrame([(f"group{i}", v) for i in range(3) for v in rd.docstring_data(stats.tukey_hsd, f"group{i}")],
                  columns=["group", "y"])
want = rd.docstring_printed(stats.tukey_hsd, "print(res)")
print(want)
printed = {m.group(1): m.groups()[1:] for m in re.finditer(r"\((\d - \d)\)\s+(-?[\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", want)}
found = read_rows(quiet(sdata.table(tk).compare_each_pair_of, "y", for_each_group_in="group", test_method="tukey"))
print(found)
ds = "SciPy tukey_hsd docstring"
pcol_tukey = [c for c in found[0] if c.startswith("p-value")][0]
for line in found:
    i, j = line["group 1"][-1], line["group 2"][-1]
    stat, p, low, high = printed[f"{i} - {j}"]
    for quantity, value, shown in [("difference", line["difference"], stat), ("p-value", line[pcol_tukey], p),
                                   ("difference 95 % CI low", line["difference 95 % CI low"], low),
                                   ("difference 95 % CI high", line["difference 95 % CI high"], high)]:
        checks.printed(ds, 'compare_each_pair_of (test_method="tukey")', f"{quantity} {i} - {j}", shown,
                       "SciPy tukey_hsd docstring, printed")
if HAS_R:
    r = r_run("""t <- TukeyHSD(aov(y ~ group, data = data.frame(y = y, group = factor(group))))$group
        list(name = rownames(t), diff = t[, 1], low = t[, 2], high = t[, 3], p = t[, 4])""", tk)
    for name, diff, low, high, p in zip(r["name"], r["diff"], r["low"], r["high"], r["p"]):
        g2, g1 = name.split("-")      # R gives "group1-group0": group1 minus group0
        line = pair_row(found, g1, g2)
        sign = 1 if line["group 1"] == g2 else -1   # simplydata gives group 1 minus group 2
        checks.compare(ds, 'compare_each_pair_of (test_method="tukey")', f"p-value {name}", line[pcol_tukey], p,
                       "R TukeyHSD", ITERATIVE, "relative")
        checks.compare(ds, 'compare_each_pair_of (test_method="tukey")', f"difference {name}", sign * line["difference"],
                       diff, "R TukeyHSD", CLOSED, "relative")


samples = rd.docstring_data(stats.dunnett, "samples")
control = rd.docstring_data(stats.dunnett, "control")
dn = pd.DataFrame([("control", v) for v in control] + [(f"treatment {k + 1}", v) for k, s in enumerate(samples)
                                                      for v in s], columns=["group", "y"])
ds = "SciPy dunnett docstring (Dunnett 1955)"
printed_stat = [float(v) for v in re.findall(r"-?\d+\.\d+", rd.docstring_printed(stats.dunnett, "res.statistic", 0))]
printed_p = [float(v) for v in re.findall(r"-?\d+\.\d+", rd.docstring_printed(stats.dunnett, "res.pvalue", 0))]
found = quiet(sdata.table(dn).compare_each_pair_of, "y", for_each_group_in="group", against="control",
              with_seed=1)
found = {row["group 1"]: row for row in read_rows(found)}
print(found)
for k in range(2):
    g = f"treatment {k + 1}"
    checks.printed(ds, "compare_each_pair_of (against=, with_seed=1)", f"statistic {g}", f"{printed_stat[k]:.8f}",
                   "SciPy dunnett docstring, printed")
    checks.compare(ds, "compare_each_pair_of (against=, with_seed=1)", f"p-value {g}", found[g]["p-value (adjusted)"],
                   printed_p[k], "SciPy dunnett docstring, printed (8 digits; Monte Carlo)", 1e-3)
    parent = stats.dunnett(*[np.asarray(s) for s in samples], control=np.asarray(control),
                           rng=np.random.default_rng(1))
    checks.compare(ds, "compare_each_pair_of (against=, with_seed=1)", f"p-value {g}", found[g]["p-value (adjusted)"],
                   parent.pvalue[k], "SciPy dunnett, same seed", EXACT, "relative")
if HAS_R:
    r = r_run("""suppressPackageStartupMessages(library(multcomp))
        d <- data.frame(y = y, group = relevel(factor(group), ref = "control"))
        s <- summary(glht(aov(y ~ group, data = d), linfct = mcp(group = "Dunnett")), test = adjusted("single-step"))
        list(t = unname(s$test$tstat), p = unname(s$test$pvalues))""", dn)
    for k in range(2):
        g = f"treatment {k + 1}"
        checks.compare(ds, "compare_each_pair_of (against=)", f"statistic {g}", found[g]["statistic"], r["t"][k],
                       "R multcomp glht Dunnett", CLOSED, "relative")
        checks.compare(ds, "compare_each_pair_of (against=)", f"p-value {g}", found[g]["p-value (adjusted)"],
                       r["p"][k], "R multcomp glht Dunnett (both by integration; 1e-3)", 1e-3)


plants = rd.r_dataset("PlantGrowth")
insects = rd.r_dataset("InsectSprays")
R_PAIRS = {
    "welch": ("pairwise.t.test(y, g, pool.sd = FALSE, p.adjust.method = 'holm')$p.value", "R pairwise.t.test(pool.sd = FALSE)"),
    "mann-whitney": ("suppressWarnings(pairwise.wilcox.test(y, g, p.adjust.method = 'holm'))$p.value",
                     "R pairwise.wilcox.test"),
    "games-howell": ("{ m <- PMCMRplus::gamesHowellTest(y, factor(g)); m$p.value }", "R PMCMRplus gamesHowellTest"),
    "lsd": ("pairwise.t.test(y, g, pool.sd = TRUE, p.adjust.method = 'none')$p.value",
            "R pairwise.t.test(pool.sd = TRUE, no adjustment)"),
    "dunn": ("suppressWarnings(PMCMRplus::kwAllPairsDunnTest(y, factor(g), p.adjust.method = 'holm'))$p.value",
             "R PMCMRplus kwAllPairsDunnTest"),
    "conover": ("suppressWarnings(PMCMRplus::kwAllPairsConoverTest(y, factor(g), p.adjust.method = 'holm'))$p.value",
                "R PMCMRplus kwAllPairsConoverTest"),
    "scheffe": ("""{ a <- anova(lm(y ~ g)); mse <- a$`Mean Sq`[2]; dfe <- a$Df[2]; k <- length(unique(g))
                     m <- tapply(y, g, mean); n <- tapply(y, g, length); lv <- names(m)
                     out <- matrix(NA, k - 1, k - 1, dimnames = list(lv[-1], lv[-k]))
                     for (i in 2:k) for (j in 1:(i - 1)) { F <- (m[i] - m[j])^2 / (mse * (1 / n[i] + 1 / n[j])) / (k - 1)
                         out[lv[i], lv[j]] <- pf(F, k - 1, dfe, lower.tail = FALSE) }
                     out }""", "Scheffe formula written out in R"),
}
for name, frame, column, group in [("PlantGrowth (Dobson 1983)", plants, "weight", "group"),
                                   ("InsectSprays (Beall 1942)", insects, "count", "spray")]:
    t = sdata.table(frame)
    for test, (code, ref) in R_PAIRS.items():
        opts = {"based_on": "ranks"} if test == "mann-whitney" else {"test_method": test}
        found = read_rows(quiet(t.compare_each_pair_of, column, for_each_group_in=group, **opts))
        pcol = [c for c in found[0] if c.startswith("p-value (") or c == "p-value (adjusted)"]
        pcol = pcol[0] if pcol else "p-value"
        if not HAS_R:
            break
        r = r_run(f"""m <- {code}; list(p = as.vector(m), row = rep(rownames(m), ncol(m)), col = rep(colnames(m), each = nrow(m)))""",
                  {"y": frame[column], "g": frame[group]})
        option = 'based_on="ranks"' if test == "mann-whitney" else f'test_method="{test}"'
        label = f"compare_each_pair_of ({option})"
        for p, a, b in zip(r["p"], r["row"], r["col"]):
            if np.isnan(p):
                continue
            line = pair_row(found, a, b)
            tol = ITERATIVE if test == "games-howell" else CLOSED
            checks.compare(name, label, f"adjusted p-value {a} - {b}", line[pcol], p, ref, tol, "relative")


from statsmodels.stats.multitest import multipletests

for name, frame, column, group in [("PlantGrowth (Dobson 1983)", plants, "weight", "group"),
                                   ("InsectSprays (Beall 1942)", insects, "count", "spray")]:
    y = frame[column].to_numpy(dtype=float)
    g = frame[group].to_numpy()
    ranks = stats.rankdata(y)
    N = len(y)
    _, ties = np.unique(y, return_counts=True)
    variance = N * (N + 1) / 12 - (ties ** 3 - ties).sum() / (12 * (N - 1))
    levels = sorted(set(g))
    pairs = [(a, b) for i, a in enumerate(levels) for b in levels[i + 1:]]
    raw = []
    for a, b in pairs:
        ra, rb = ranks[g == a], ranks[g == b]
        z = (ra.mean() - rb.mean()) / np.sqrt(variance * (1 / len(ra) + 1 / len(rb)))
        raw.append(2 * stats.norm.sf(abs(z)))
    holm = multipletests(raw, method="holm")[1]
    found = read_rows(quiet(sdata.table(frame).compare_each_pair_of, column, for_each_group_in=group, test_method="dunn"))
    pcol = [c for c in found[0] if c.startswith("p-value (")][0]
    for (a, b), p in zip(pairs, holm):
        line = pair_row(found, a, b)
        checks.compare(name, 'compare_each_pair_of (test_method="dunn")', f"adjusted p-value {a} - {b}", line[pcol],
                       p, "Dunn (1964) written out, NumPy + Holm", CLOSED, "relative")


for name, frame, column, group in [("PlantGrowth (Dobson 1983)", plants, "weight", "group"),
                                   ("InsectSprays (Beall 1942)", insects, "count", "spray")]:
    t = sdata.table(frame)
    letters = read_rows(quiet(t.compare_each_pair_of, column, for_each_group_in=group, show="letters"))
    pairs = read_rows(quiet(t.compare_each_pair_of, column, for_each_group_in=group))
    print(letters)
    share = {(la["group"], lb["group"]): bool(set(la["letters"]) & set(lb["letters"])) for la in letters for lb in letters}
    for line in pairs:
        a, b = line["group 1"], line["group 2"]
        checks.compare(name, 'compare_each_pair_of (show="letters")', f"{a} and {b} share a letter",
                       float(share[(a, b)]), float(line["p-value (Holm)"] >= 0.05),
                       "adjusted p-value of 0.05 or more (Python)", 0)


made = rd.r_example_objects("stats", "friedman.test")
rt = r_run("list(a = RoundingTimes[, 1], b = RoundingTimes[, 2], c = RoundingTimes[, 3], n = colnames(RoundingTimes))")
rounding = pd.DataFrame({rt["n"][0]: rt["a"], rt["n"][1]: rt["b"], rt["n"][2]: rt["c"]})
methods = list(rounding.columns)
ds = "rounding first base (R friedman.test example; Hollander & Wolfe 1973)"
found = one(sdata.table(rounding).compare_matched_columns, *methods)
n, k = rounding.shape
if HAS_R:
    r = r_run("""f <- friedman.test(RoundingTimes); w <- effectsize::kendalls_w(RoundingTimes, verbose = FALSE)
        list(statistic = unname(f$statistic), `p-value` = f$p.value, w = w$Kendalls_W[1])""")
    compare_many(ds, "compare_matched_columns", found, {"statistic": r["statistic"], "p-value": r["p-value"]},
                 "R friedman.test")
    checks.compare(ds, "compare_matched_columns", "effect size (Kendall's W)", found["effect size"],
                   r["statistic"] / (n * (k - 1)), "chi-square / (n (k - 1))", CLOSED, "relative")
    checks.compare(ds, "compare_matched_columns", "effect size (Kendall's W)", found["effect size"], r["w"],
                   "R effectsize kendalls_w", CLOSED, "relative")
    found = one(sdata.table(rounding).compare_matched_columns, *methods, test_method="quade")
    r = r_run("""y <- RoundingTimes; b <- nrow(y); k <- ncol(y)
        Q <- rank(apply(y, 1, function(v) max(v) - min(v)))
        R <- t(apply(y, 1, rank)); S <- Q * (R - (k + 1) / 2)
        A <- sum(S^2); B <- sum(colSums(S)^2) / b; F <- (b - 1) * B / (A - B)
        q <- quade.test(y)
        list(F = F, p = pf(F, k - 1, (b - 1) * (k - 1), lower.tail = FALSE), qF = unname(q$statistic), qp = q$p.value)""")
    compare_many(ds, 'compare_matched_columns (test_method="quade")', found, {"statistic": r["F"], "p-value": r["p"]},
                 "Quade's statistic written out in R")
    compare_many(ds, 'compare_matched_columns (test_method="quade")', found, {"statistic": r["qF"], "p-value": r["qp"]},
                 "R quade.test")

table = np.asarray(rd.docstring_data(stats.page_trend_test, "table"))
page = pd.DataFrame(table, columns=["first", "second", "third"])
names_, shown_ = rd.docstring_example(stats.page_trend_test)
want = next(w for s, w in shown_ if w.startswith("PageTrendTestResult"))
ds = "SciPy page_trend_test docstring (Page 1963)"
found = one(sdata.table(page).compare_matched_columns, "first", "second", "third", test_method="page")
checks.compare(ds, 'compare_matched_columns (test_method="page")', "statistic (L)", found["statistic"],
               rd.printed_number(want, "statistic"), "SciPy page_trend_test docstring, printed", DOCSTRING, "relative")
checks.compare(ds, 'compare_matched_columns (test_method="page")', "p-value", found["p-value"],
               rd.printed_number(want, "pvalue"), "SciPy page_trend_test docstring, printed", DOCSTRING, "relative")
if HAS_R:
    r = r_run("""p <- crank::page.trend.test(cbind(a, b, c), ranks = FALSE); list(L = unname(p$L))""",
              {"a": page["first"], "b": page["second"], "c": page["third"]})
    checks.compare(ds, 'compare_matched_columns (test_method="page")', "statistic (L)", found["statistic"], r["L"],
                   "R crank page.trend.test", CLOSED, "relative")

if HAS_R:
    r_run("list(ok = 1)", packages=("DescTools",))
    made = rd.r_example_objects("DescTools", "CochranQTest")
    cq = r_run("list(A = as.character(d.frm$A), B = as.character(d.frm$B), C = as.character(d.frm$C))")
    answers = pd.DataFrame(cq)
    found = one(sdata.table(answers).compare_matched_columns, "A", "B", "C")
    r = r_run("""q <- CochranQTest(resp ~ time | id, data = d.long); list(statistic = unname(q$statistic), `p-value` = q$p.value)""")
    compare_many("SAS FREQ example (R DescTools CochranQTest example)", "compare_matched_columns (yes/no columns: Cochran's Q)",
                 found, r, "R DescTools CochranQTest")
checks.skip("rounding first base", 'compare_matched_columns (based_on="averages")', "repeated-measures ANOVA",
            "R", 'simplydata says based_on="averages" is not here yet (it points to compare_averages_of(matched_by=))')


hedenfalk = load_csv("hedenfalk_p.csv")
hedenfalk_300 = hedenfalk.iloc[:300].rename(columns={hedenfalk.columns[0]: "p-value"})
pcol = hedenfalk.columns[0]
ds = "Hedenfalk et al. (2001) p-values"
for method, r_method, size in [("holm", "holm", None), ("bonferroni", "bonferroni", None),
                               ("benjamini-hochberg", "BH", None), ("bh", "BH", None), ("fdr", "BH", None),
                               ("benjamini-yekutieli", "BY", None), ("hochberg", "hochberg", None),
                               ("hommel", "hommel", 300)]:
    data = hedenfalk.iloc[:size] if size else hedenfalk
    t = quiet(sdata.table(data.rename(columns={pcol: "p-value"})).add_adjusted_p_value_column, "adjusted",
              test_method=method)
    adjusted = np.asarray(read_numbers(quiet(t.get_column, "adjusted")), dtype=float)
    if not HAS_R:
        break
    r = r_run(f'list(q = p.adjust(p, method = "{r_method}"))', {"p": data[pcol]})["q"]
    diffs = np.abs(adjusted - np.asarray(r)) / np.maximum(np.asarray(r), 1e-300)
    worst = int(np.argmax(diffs))
    checks.vector(ds, f'add_adjusted_p_value_column (test_method="{method}")', f"adjusted p-values ({len(data)} rows)",
                  r, 'R p.adjust', CLOSED)
if HAS_R:
    for method, column in [("storey", "q"), ("q-values", "q"), ("local fdr", "lfdr"), ("lfdr", "lfdr")]:
        try:
            t = quiet(sdata.table(hedenfalk.rename(columns={pcol: "p-value"})).add_adjusted_p_value_column,
                      test_method=method)
        except Exception as problem:  # noqa: BLE001
            checks.skip(ds, f'add_adjusted_p_value_column (test_method="{method}")', column, "R qvalue", str(problem)[:150])
            continue
        new = t.get_column_names()[-1]
        values = np.asarray(read_numbers(quiet(t.get_column, new)), dtype=float)
        r = r_run(f"q <- qvalue::qvalue(p); list(v = q${'qvalues' if column == 'q' else 'lfdr'})", {"p": hedenfalk[pcol]})["v"]
        diffs = np.abs(values - np.asarray(r)) / np.maximum(np.asarray(r), 1e-300)
        worst = int(np.argmax(diffs))
        checks.vector(ds, f'add_adjusted_p_value_column (test_method="{method}")', f"{column} ({len(values)} rows)",
                      r, "R qvalue (same package)", CLOSED)


mtcars = rd.r_dataset("mtcars")
mtcars["transmission"] = np.where(mtcars["am"] == 1, "manual", "automatic")
cars = sdata.table(mtcars)
found = one(cars.compare_correlations, "mpg", "wt", for_each_group_in="transmission")
print(found)
groups = [mtcars[mtcars["transmission"] == g] for g in ("automatic", "manual")]
rs = [stats.pearsonr(g["mpg"], g["wt"]).statistic for g in groups]
ns = [len(g) for g in groups]
z = (np.arctanh(rs[0]) - np.arctanh(rs[1])) / np.sqrt(1 / (ns[0] - 3) + 1 / (ns[1] - 3))
compare_many("mtcars mpg, wt by transmission", "compare_correlations (for_each_group_in=)", found,
             {"statistic": z, "p-value": 2 * stats.norm.sf(abs(z))}, "Fisher's z written out, NumPy")
found = one(cars.compare_correlations, ["mpg", "wt"], ["mpg", "hp"])
r12, r13 = stats.pearsonr(mtcars["mpg"], mtcars["wt"]).statistic, stats.pearsonr(mtcars["mpg"], mtcars["hp"]).statistic
r23 = stats.pearsonr(mtcars["wt"], mtcars["hp"]).statistic
n = len(mtcars)
rbar = (r12 + r13) / 2
psi = r23 * (1 - 2 * rbar ** 2) - 0.5 * rbar ** 2 * (1 - 2 * rbar ** 2 - r23 ** 2)
c = psi / (1 - rbar ** 2) ** 2
z = (np.arctanh(r12) - np.arctanh(r13)) * np.sqrt(n - 3) / np.sqrt(2 - 2 * c)
compare_many("mtcars mpg with wt and with hp", "compare_correlations (two overlapping pairs)", found,
             {"statistic": z, "p-value": 2 * stats.norm.sf(abs(z))}, "Steiger (1980) eq. 15 written out, NumPy")


plants_t = sdata.table(plants)
cases = [("compare_each_pair_of", lambda: plants_t.compare_each_pair_of("weight", for_each_group_in="group",
                                                                       test_method="turkey"), "tukey"),
         ("add_adjusted_p_value_column", lambda: sdata.table(hedenfalk.rename(columns={pcol: "p-value"}))
          .add_adjusted_p_value_column(test_method="hollm"), "holm"),
         ("find_percent_of_rows_where", lambda: plants_t.find_percent_of_rows_where("weight", is_above=5,
                                                                                   interval_by="wilsn"), "wilson")]
for word, call, right in cases:
    try:
        quiet(call)
        said = ""
    except Exception as problem:  # noqa: BLE001
        said = str(problem)
    checks.compare("PlantGrowth / Hedenfalk", f"{word} (misspelt value)", f'stops and suggests "{right}"',
                   float(f'"{right}"' in said and "Did you mean" in said), 1.0, "simplydata's message", 0)


checks.data("job", job)
checks.data("smoke", smoke)
checks.data("eyes", eyes)
checks.data("pea_rows", pea_rows)
checks.data("tk", tk)
checks.data("dn", dn)
checks.data("plants", plants)
checks.data("insects", insects)
checks.data("rounding", rounding)
checks.data("page", page)
checks.data("answers", answers)
checks.data("mtcars", mtcars)
checks.data("hedenfalk_300", hedenfalk_300)
checks.save()
