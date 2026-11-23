"""Reference values for notebooks/10_power_and_sample_size.ipynb, saved to references/10_power_and_sample_size.csv.

Printed values first: Cohen (1988) Tables 2.4.1 and 8.4.4, Faul et al. (2007) and the G*Power 3.1 manual, the help
pages of R power.t.test and power.prop.test, the pwr vignette, the Stata power manuals ([PSS-2] power logrank,
power twomeans cluster), the presize paper (Haynes et al. 2021) and Bland's note on cluster trials. Then the
original libraries run once: statsmodels' power classes (the Python parent), R power.t.test and power.prop.test
(tol = 1e-12), pwr, TOSTER, presize, gsDesign, powerSurvEpi and pwrss, and SciPy stats.power with the same seed.

Tolerances: a whole number of values exact; a power at a given number of values, or events from a closed formula,
relative 1e-10; an exact number of values from a root search, relative 1e-6; a printed value as before (half a
unit of its last printed digit).

pwr finds n with uniroot at its default tolerance (about 1e-4 in n). The first layout checked pwr's power at
simplydata's exact n; here (where simplydata's numbers are not used) the exact n is solved from pwr's own power
function with uniroot(tol = 1e-12), and compared with simplydata's exact n (relative 1e-6).

The power "with that many" is computed at the whole number of values the reference itself finds (the ceiling of
its exact n); simplydata's whole number is compared with that same ceiling.

    python scripts/references/10_power_and_sample_size.py      (from the repository root; needs R)
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.power import FTestAnovaPower, FTestPowerF2, GofChisquarePower, TTestIndPower, TTestPower  # noqa: E402

from checks import r_available, r_run  # noqa: E402

STEM = "10_power_and_sample_size"
R_PACKAGES = ("pwr", "TOSTER", "presize", "gsDesign", "powerSurvEpi", "pwrss")
checks = Recorder(STEM, R_PACKAGES)
HAS_R = r_available()
CLOSED = 1e-10
SOLVED = 1e-6
COHEN = "https://doi.org/10.4324/9780203771587"
FAUL = "https://doi.org/10.3758/BF03193146"
GPOWER = "https://www.psychologie.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Psychologie/AAP/gpower/GPowerManual.pdf"
STATA_LOGRANK = "https://www.stata.com/manuals/psspowerlogrank.pdf"
STATA_CLUSTER = "https://www.stata.com/manuals/psspowertwomeanscluster.pdf"
HAYNES = "https://doi.org/10.21105/joss.03118"
BLAND = "https://www-users.york.ac.uk/~mb55/msc/clinbio/week9/samp_size.pdf"
R_STATS = "https://search.r-project.org/R/refmans/stats/html/{}.html"
PWR_VIGNETTE = "https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html"


def rel(dataset, method, quantity, value, source, tolerance):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, "relative")


def exact(dataset, method, quantity, value, source):
    checks.compare(dataset, method, quantity, None, value, source, 0)


def printed(dataset, method, quantity, value, source, tolerance, url):
    """A printed value with the tolerance the first layout used (absolute)."""
    checks.compare(dataset, method, quantity, None, value, source, tolerance)
    checks.rows[-1]["url"] = url
    checks.rows[-1]["kind"] = "printed"
    checks.rows[-1]["digits as printed"] = repr(value)


def pwr_root(call_with_n, low):
    """The n at which an R pwr call (with `n` left as the variable) gives power 0.8, by uniroot(tol = 1e-12)."""
    return r_run(f"f <- function(n) {call_with_n}$power - 0.8\n"
                 f"list(n = uniroot(f, c({low}, 1e6), tol = 1e-12, extendInt = 'upX')$root)")["n"]


# ---- two averages
def t_two(d, power=0.8, sides=2, label=None):
    dataset = label or f"d = {d}, power {power}, {sides}-sided"
    method = "find_how_many_needed (two averages)"
    alpha = 0.05 / 2 if sides == 2 else 0.05
    n = TTestIndPower().solve_power(effect_size=d, power=power, alpha=alpha, alternative="larger")
    rel(dataset, method, "values per group (exact)", n, "statsmodels TTestIndPower", SOLVED)
    rel(dataset, method, "power with that many",
        TTestIndPower().power(effect_size=d, nobs1=math.ceil(n), alpha=alpha, alternative="larger"),
        "statsmodels TTestIndPower", CLOSED)
    if HAS_R:
        alt = "two.sided" if sides == 2 else "one.sided"
        r = r_run(f"a <- power.t.test(delta = {d}, sd = 1, power = {power}, alternative = '{alt}', tol = 1e-12)\n"
                  f"b <- power.t.test(n = ceiling(a$n), delta = {d}, sd = 1, alternative = '{alt}')\n"
                  "list(n = a$n, whole = ceiling(a$n), power = b$power)")
        rel(dataset, method, "values per group (exact)", r["n"], "R power.t.test", SOLVED)
        exact(dataset, method, "values per group", r["whole"], "R power.t.test")
        rel(dataset, method, "power with that many", r["power"], "R power.t.test", CLOSED)


for d, value in ((0.2, 393), (0.5, 64), (0.8, 26)):
    t_two(d)
    printed(f"d = {d}, power 0.8, 2-sided", "find_how_many_needed (two averages)", "values per group", value,
            "Cohen (1988), Table 2.4.1", 0, COHEN)
t_two(0.5, power=0.95)
printed("d = 0.5, power 0.95, 2-sided", "find_how_many_needed (two averages)", "values per group", 105,
        "Faul et al. (2007), G*Power 3", 0, FAUL)
t_two(1.0, power=0.9)
printed("d = 1.0, power 0.9, 2-sided", "find_how_many_needed (two averages)", "values per group (exact)", 22.02110,
        "R power.t.test help page, printed", 5e-6, R_STATS.format("power.t.test"))
t_two(0.5, sides=1)
t_two(0.5, power=0.9, label="difference 5, spread 10, power 0.9")

for d, n in ((1.0, 20), (0.5, 30), (0.3, 100)):
    dataset, method = f"d = {d}, {n} per group", "find_chance_of_finding (two averages)"
    rel(dataset, method, "power", TTestIndPower().power(effect_size=d, nobs1=n, alpha=0.025, alternative="larger"),
        "statsmodels TTestIndPower", CLOSED)
    if HAS_R:
        r = r_run(f"list(power = power.t.test(n = {n}, delta = {d}, sd = 1)$power)")
        rel(dataset, method, "power", r["power"], "R power.t.test", CLOSED)
    if (d, n) == (1.0, 20):
        printed(dataset, method, "power", 0.8689528, "R power.t.test help page, printed", 5e-8,
                R_STATS.format("power.t.test"))

# ---- one group and before and after
for test, kind in (("compare_average_to", "one.sample"), ("compare_before_after", "paired")):
    for d in (0.3, 0.5):
        dataset, method = f"d = {d}, power 0.8", f"find_how_many_needed ({test})"
        n = TTestPower().solve_power(effect_size=d, power=0.8, alpha=0.025, alternative="larger")
        rel(dataset, method, "values (exact)", n, "statsmodels TTestPower", SOLVED)
        rel(dataset, method, "power with that many",
            TTestPower().power(effect_size=d, nobs=math.ceil(n), alpha=0.025, alternative="larger"),
            "statsmodels TTestPower", CLOSED)
        if HAS_R:
            r = r_run(f"a <- power.t.test(delta = {d}, sd = 1, power = 0.8, type = '{kind}', tol = 1e-12)\n"
                      f"b <- power.t.test(n = ceiling(a$n), delta = {d}, sd = 1, type = '{kind}')\n"
                      "list(n = a$n, whole = ceiling(a$n), power = b$power)")
            rel(dataset, method, "values (exact)", r["n"], f"R power.t.test ({kind})", SOLVED)
            exact(dataset, method, "values", r["whole"], f"R power.t.test ({kind})")
            rel(dataset, method, "power with that many", r["power"], f"R power.t.test ({kind})", CLOSED)

# ---- ANOVA
for f, k in ((0.25, 4), (0.4, 3), (0.1, 5)):
    dataset, method = f"f = {f}, {k} groups, power 0.8", "find_how_many_needed (ANOVA)"
    n = FTestAnovaPower().solve_power(effect_size=f, power=0.8, alpha=0.05, k_groups=k) / k
    rel(dataset, method, "values per group (exact)", n, "statsmodels FTestAnovaPower", SOLVED)
    rel(dataset, method, "power with that many",
        FTestAnovaPower().power(effect_size=f, nobs=k * math.ceil(n), alpha=0.05, k_groups=k),
        "statsmodels FTestAnovaPower", CLOSED)
    if HAS_R:
        rel(dataset, method, "values per group (exact)", pwr_root(f"pwr::pwr.anova.test(k = {k}, n = n, f = {f})", 2),
            "R pwr.anova.test power, solved by uniroot (tol 1e-12)", SOLVED)
        r = r_run(f"n <- ceiling(pwr::pwr.anova.test(k = {k}, f = {f}, power = 0.8)$n)\n"
                  f"list(n = n, power = pwr::pwr.anova.test(k = {k}, n = n, f = {f})$power)")
        exact(dataset, method, "values per group", r["n"], "R pwr.anova.test")
        rel(dataset, method, "power with that many", r["power"], "R pwr.anova.test", CLOSED)
    if (f, k) == (0.25, 4):
        printed(dataset, method, "values per group", 45, "Cohen (1988), Table 8.4.4", 0, COHEN)

# ---- chi-square
for w, groups in ((0.3, (2, 2)), (0.3, (3, 4)), (0.1, (2, 3))):
    df = (groups[0] - 1) * (groups[1] - 1)
    dataset, method = f"w = {w}, {groups[0]} x {groups[1]} table, power 0.8", "find_how_many_needed (chi-square)"
    n = GofChisquarePower().solve_power(effect_size=w, power=0.8, alpha=0.05, n_bins=df + 1) / groups[0]
    rel(dataset, method, "values per group (exact)", n, "statsmodels GofChisquarePower", SOLVED)
    rel(dataset, method, "power with that many",
        GofChisquarePower().power(effect_size=w, nobs=groups[0] * math.ceil(n), alpha=0.05, n_bins=df + 1),
        "statsmodels GofChisquarePower", CLOSED)
    if HAS_R:
        N = pwr_root(f"pwr::pwr.chisq.test(w = {w}, N = n, df = {df})", df + 1)
        rel(dataset, method, "values per group (exact)", N / groups[0],
            "R pwr.chisq.test power, solved by uniroot (tol 1e-12)", SOLVED)
        r = r_run(f"list(power = pwr::pwr.chisq.test(w = {w}, N = {groups[0] * math.ceil(n)}, df = {df})$power)")
        rel(dataset, method, "power with that many", r["power"], "R pwr.chisq.test", CLOSED)

# ---- two percents
for p1, p2, power in ((50, 75, 0.9), (30, 15, 0.8), (10, 5, 0.8)):
    dataset, method = f"{p1} % against {p2} %, power {power}", "find_how_many_needed (two percents)"
    if HAS_R:
        r = r_run(f"a <- power.prop.test(p1 = {p1 / 100}, p2 = {p2 / 100}, power = {power}, tol = 1e-12)\n"
                  f"b <- power.prop.test(n = ceiling(a$n), p1 = {p1 / 100}, p2 = {p2 / 100})\n"
                  "list(n = a$n, whole = ceiling(a$n), power = b$power)")
        rel(dataset, method, "values per group (exact)", r["n"], "R power.prop.test", SOLVED)
        exact(dataset, method, "values per group", r["whole"], "R power.prop.test")
        rel(dataset, method, "power with that many", r["power"], "R power.prop.test", CLOSED)
    if (p1, p2) == (50, 75):
        printed(dataset, method, "values per group (exact)", 76.70693, "R power.prop.test help page, printed", 5e-6,
                R_STATS.format("power.prop.test"))
dataset, method = "50 % against 75 %, 50 per group", "find_chance_of_finding (two percents)"
printed(dataset, method, "power", 0.7401659, "R power.prop.test help page, printed", 5e-8, R_STATS.format("power.prop.test"))
if HAS_R:
    r = r_run("list(power = power.prop.test(n = 50, p1 = .50, p2 = .75)$power)")
    rel(dataset, method, "power", r["power"], "R power.prop.test", CLOSED)

# ---- correlation
if HAS_R:
    for r_size in (0.1, 0.3, 0.5):
        dataset, method = f"r = {r_size}, power 0.8", "find_how_many_needed (correlation)"
        rel(dataset, method, "pairs (exact)", pwr_root(f"pwr::pwr.r.test(n = n, r = {r_size})", 4),
            "R pwr.r.test power, solved by uniroot (tol 1e-12)", SOLVED)
        r = r_run(f"n <- ceiling(pwr::pwr.r.test(r = {r_size}, power = 0.8)$n)\n"
                  f"list(n = n, power = pwr::pwr.r.test(n = n, r = {r_size})$power)")
        exact(dataset, method, "pairs", r["n"], "R pwr.r.test")
        rel(dataset, method, "power with that many", r["power"], "R pwr.r.test", CLOSED)
        r = r_run(f"list(power = pwr::pwr.r.test(n = 50, r = {r_size})$power)")
        rel(f"r = {r_size}, 50 pairs", "find_chance_of_finding (correlation)", "power", r["power"], "R pwr.r.test", CLOSED)

# ---- regression
dataset, method = "f² = 0.1/0.9, 5 columns, 95 rows", "find_chance_of_finding (regression)"
printed(dataset, method, "power", 0.673586, "G*Power 3.1 manual 13.3.1, printed", 5e-7, GPOWER)
rel(dataset, method, "power", FTestPowerF2().power(effect_size=0.1 / 0.9, df_num=5, df_denom=89, alpha=0.05),
    "statsmodels FTestPowerF2", CLOSED)
if HAS_R:
    rel(dataset, method, "power", r_run("list(p = pwr::pwr.f2.test(u = 5, v = 89, f2 = 0.1 / 0.9)$power)")["p"],
        "R pwr.f2.test", CLOSED)

dataset, method = "f² = 0.15, 3 columns, power 0.8", "find_how_many_needed (regression)"
rel(dataset, method, "rows (exact)", FTestPowerF2().solve_power(effect_size=0.15, df_num=3, alpha=0.05, power=0.8) + 3 + 1,
    "statsmodels FTestPowerF2", SOLVED)
if HAS_R:
    rel(dataset, method, "rows (exact)", pwr_root("pwr::pwr.f2.test(u = 3, v = n, f2 = 0.15)", 1) + 4,
        "R pwr.f2.test power, solved by uniroot (tol 1e-12)", SOLVED)

dataset, method = "2 columns added to 1, f² = 0.3/0.7, significance 0.001", "find_how_many_needed (compare_models)"
printed(dataset, method, "rows (exact)", 49.88971 + 4, "pwr vignette, printed v + 4", 5e-6, PWR_VIGNETTE)
rel(dataset, method, "rows (exact)", FTestPowerF2().solve_power(effect_size=0.3 / 0.7, df_num=2, alpha=0.001, power=0.8) + 4,
    "statsmodels FTestPowerF2", SOLVED)

# ---- logistic regression, equivalence, interval width
dataset, method = "odds ratio 1.5, 50 %, power 0.95", "find_how_many_needed (logistic regression)"
printed(dataset, method, "values", 337, "G*Power 3.1 manual 29.4, printed", 0, GPOWER)
if HAS_R:
    r = r_run("r <- pwrss::power.z.logistic(base.prob = 0.5, odds.ratio = 1.5, power = 0.95, alpha = 0.05, "
              "distribution = 'normal', method = 'demidenko(vc)', ceil.n = FALSE, verbose = 0)\n"
              "list(n = r$n)")
    rel(dataset, method, "values (exact)", r["n"], "R pwrss::power.z.logistic", SOLVED)

dataset, method = "within 5, spread 10, power 0.8", "find_how_many_needed (equivalence, TOST)"
if HAS_R:
    r = r_run("a <- TOSTER::power_t_TOST(delta = 0, sd = 10, low_eqbound = -5, high_eqbound = 5, alpha = 0.05, "
              "power = 0.8, type = 'two.sample')\n"
              "b <- TOSTER::power_t_TOST(n = ceiling(a$n), delta = 0, sd = 10, low_eqbound = -5, "
              "high_eqbound = 5, alpha = 0.05, type = 'two.sample')\n"
              "list(n = a$n, whole = ceiling(a$n), power = b$power)")
    rel(dataset, method, "values per group (exact)", r["n"], "R TOSTER::power_t_TOST", SOLVED)
    exact(dataset, method, "values per group", r["whole"], "R TOSTER::power_t_TOST")
    rel(dataset, method, "power with that many", r["power"], "R TOSTER::power_t_TOST", CLOSED)
    r = r_run("list(n = presize::prec_mean(mean = 0, sd = 10, conf.width = 4)$n)")
    rel("interval 4 wide, spread 10", "find_how_many_needed (interval width)", "values (exact)", r["n"],
        "R presize::prec_mean", SOLVED)
dataset, method = "Wilson interval 5 points wide near 15 %", "find_how_many_needed (interval width)"
printed(dataset, method, "values (exact)", 783.4897, "Haynes et al. (2021), printed", 5e-5, HAYNES)
printed(dataset, method, "values", 784, "Haynes et al. (2021), printed", 0, HAYNES)
if HAS_R:
    r = r_run("list(n = presize::prec_prop(p = 0.15, conf.width = 0.05, method = 'wilson')$n)")
    rel(dataset, method, "values (exact)", r["n"], "R presize::prec_prop", SOLVED)


# ---- survival
def events(hr, sides, value=None):
    dataset, method = f"hazard ratio {hr:.4g}, {sides}-sided", "find_how_many_needed (survival, events)"
    if value is not None:
        printed(dataset, method, "events", value, "Stata [PSS-2] power logrank, printed", 0, STATA_LOGRANK)
    if HAS_R:
        r = r_run(f"d <- gsDesign::nEvents(hr = {hr!r}, alpha = 0.05, beta = 0.2, ratio = 1, sided = {sides})\n"
                  "list(d = d, whole = ceiling(d))")
        rel(dataset, method, "events (exact)", r["d"], "R gsDesign::nEvents", CLOSED)
        exact(dataset, method, "events", r["whole"], "R gsDesign::nEvents")


events(0.5, 2, 66)
events(math.log(0.5) / math.log(0.6), 1, 266)
events(0.7, 2)
if HAS_R:
    r = r_run("list(n = powerSurvEpi::ssizeEpi.default(power = 0.8, theta = 0.7, p = 0.5, psi = 0.35, rho2 = 0))")
    exact("hazard ratio 0.7, 30 % and 40 % with the event", "find_how_many_needed (survival, people)", "people in all",
          r["n"], "R powerSurvEpi::ssizeEpi.default")
dataset, method = "hazard ratio 0.5, everyone followed to the event", "find_how_many_needed (survival, Freedman)"
printed(dataset, method, "people per group", 36, "Stata [PSS-2] power logrank, printed", 0, STATA_LOGRANK)
if HAS_R:
    r = r_run("n <- powerSurvEpi::ssizeCT.default(power = 0.8, k = 1, pE = 1, pC = 1, RR = 0.5, alpha = 0.05)\n"
              "list(n = n[[1]])")
    exact(dataset, method, "people per group", r["n"], "R powerSurvEpi::ssizeCT.default")

# ---- clusters and the rank test
dataset, method = "clusters of 20, ICC 0.025", "find_how_many_needed (clusters)"
printed(dataset, method, "clusters per group", 13, "Stata [PSS-2] power twomeans cluster, printed", 0, STATA_CLUSTER)
printed(dataset, method, "values in whole clusters", 260, "Stata [PSS-2] power twomeans cluster, printed", 0, STATA_CLUSTER)
rel(dataset, method, "design effect", 1 + 19 * 0.025, "1 + (m - 1) ICC", CLOSED)
printed("clusters of 30, ICC 0.05", method, "design effect", 2.45, "Bland, sample size in cluster trials, printed", 5e-3,
        BLAND)

rng = np.random.default_rng(1)


def mann_whitney(x, y, axis=-1):
    return stats.mannwhitneyu(x, y, axis=axis, method="asymptotic", use_continuity=True).pvalue


simulated = stats.power(mann_whitney, (rng.standard_normal, lambda size: rng.normal(0.5, 1, size)), (64, 64),
                        significance=0.05, n_resamples=10000, batch=1000)
rel("d = 0.5, 64 per group", "find_chance_of_finding (ranks, simulated)", "power", float(simulated.power),
    "SciPy stats.power (seed 1)", CLOSED)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
