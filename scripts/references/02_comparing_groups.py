"""Reference values for notebooks/02_comparing_groups.ipynb, saved to references/02_comparing_groups.csv.

In order of preference: the averages Student (1908) printed for the sleep data (half a unit of the last printed
digit); the NIST StRD certified F statistics and degrees of freedom of the eleven one-way ANOVA datasets (LRE of at
least 9 digits); and, only where nothing is printed or certified, R run once here (relative difference 1e-10):
t.test, anova(lm()), oneway.test, kruskal.test, bartlett.test, fligner.test, and Levene's test written out in base R
(an ANOVA of the distances from each group's median, as car::leveneTest(center = median)).
The data the notebook reads are written to references/data/02_comparing_groups/.

    python scripts/references/02_comparing_groups.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

from checks import load_csv, r_available, r_run, read_nist  # noqa: E402

NOTEBOOK = "02_comparing_groups"
MIN_LRE = 9
RELATIVE = 1e-10
STUDENT = "Student (1908), 2 decimals, printed"
STUDENT_TITLE = "Student (1908), The probable error of a mean, Biometrika 6(1), 1-25, doi:10.1093/biomet/6.1.1"
STUDENT_URL = "https://doi.org/10.1093/biomet/6.1.1"
NIST_ANOVA_URL = "https://www.itl.nist.gov/div898/strd/anova/{name}.html"
RMAN = "https://search.r-project.org/R/refmans/stats/html/{topic}.html"

ref = Recorder(NOTEBOOK)
if not r_available():
    raise SystemExit("R is needed for the reference values of this notebook")


def run(dataset, method, values: dict, source, topic):
    for quantity, value in values.items():
        ref.compare(dataset, method, quantity, None, value, source, RELATIVE, "relative")
        ref.rows[-1]["url"] = RMAN.format(topic=topic)


# ---- Student's sleep data
sleep = load_csv("sleep.csv")
sleep["group"] = "drug " + sleep["group"].astype(str)
wide = sleep.pivot(index="ID", columns="group", values="extra").reset_index()
ref.data("sleep", sleep)
ref.data("sleep_wide", wide)
ds = "Student's sleep data"
for quantity, shown in {"average in drug 1": "0.75", "average in drug 2": "2.33", "difference": "-1.58"}.items():
    ref.printed(ds, "compare_averages_of (Welch)", quantity, shown, STUDENT, STUDENT_TITLE, STUDENT_URL)
r = r_run("""t <- t.test(extra[group == "drug 1"], extra[group == "drug 2"])
             list(statistic = t$statistic, `degrees of freedom` = t$parameter, `p-value` = t$p.value)""", sleep)
run(ds, "compare_averages_of (Welch)", r, "R t.test", "t.test")
ref.printed(ds, "compare_before_and_after", "average change", "1.58", STUDENT, STUDENT_TITLE, STUDENT_URL)
r = r_run("""t <- t.test(`drug 2`, `drug 1`, paired = TRUE)
             list(statistic = t$statistic, `degrees of freedom` = t$parameter, `p-value` = t$p.value)""",
          wide[["drug 1", "drug 2"]])
run(ds, "compare_before_and_after", r, "R t.test(paired = TRUE)", "t.test")

# ---- NIST one-way ANOVA
for name in ["SiRstv", "SmLs01", "SmLs02", "SmLs03", "AtmWtAg", "SmLs04", "SmLs05", "SmLs06", "SmLs07", "SmLs08",
             "SmLs09"]:
    d = read_nist(name)
    ref.data(f"nist_{name}", d.data)
    ds = f"NIST {name}"
    for quantity, key, tolerance, how in (("F", "F", MIN_LRE, "lre"), ("between df", "between df", 0, "absolute"),
                                          ("within df", "within df", 0, "absolute")):
        ref.compare(ds, "compare_averages_of (ANOVA)", quantity, None, d.certified[key], "NIST certified", tolerance,
                    how)
        ref.rows[-1]["url"] = NIST_ANOVA_URL.format(name=name)
        ref.rows[-1]["source title"] = f"NIST StRD analysis of variance: {name}"
    r = r_run("a <- anova(lm(y ~ factor(group))); list(`p-value` = a$`Pr(>F)`[1])", d.data)
    run(ds, "compare_averages_of (ANOVA)", r, "R anova(lm())", "anova")

# ---- Welch's ANOVA on iris, Kruskal-Wallis and the spreads on warpbreaks
iris = load_csv("iris.csv")
r = r_run("""t <- oneway.test(Sepal.Length ~ Species)
             list(statistic = t$statistic, `degrees of freedom 1` = t$parameter[1],
                  `degrees of freedom 2` = t$parameter[2], `p-value` = t$p.value)""", iris)
run("iris", "compare_averages_of (Welch ANOVA)", r, "R oneway.test", "oneway.test")
warp = load_csv("warpbreaks.csv")
r = r_run("t <- kruskal.test(breaks ~ tension); list(statistic = t$statistic, `p-value` = t$p.value)", warp)
run("warpbreaks", "compare_averages_of (Kruskal-Wallis)", r, "R kruskal.test", "kruskal.test")
for test, code, source, topic in (
        ("levene", "a <- anova(lm(abs(breaks - ave(breaks, tension, FUN = median)) ~ tension)); "
                   "list(statistic = a$`F value`[1], `p-value` = a$`Pr(>F)`[1])", "R anova of distances from medians",
         "anova"),
        ("bartlett", "t <- bartlett.test(breaks ~ tension); list(statistic = t$statistic, `p-value` = t$p.value)",
         "R bartlett.test", "bartlett.test"),
        ("fligner", "t <- fligner.test(breaks ~ tension); list(statistic = t$statistic, `p-value` = t$p.value)",
         "R fligner.test", "fligner.test")):
    run("warpbreaks", f"compare_spreads_of ({test})", r_run(code, warp), source, topic)
ref.save()
