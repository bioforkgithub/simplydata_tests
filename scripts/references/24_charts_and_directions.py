"""Reference values for notebooks/24_charts_and_directions.ipynb, saved to references/24_charts_and_directions.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/24_charts_and_directions/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/24_charts_and_directions.py      (from the repository root; needs R with the packages listed)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks and realdata)

import numpy as np
import pandas as pd
import simplydata as sdata
from matplotlib.collections import LineCollection, PathCollection
from scipy import optimize, signal, special, stats

import realdata as rd
from checks import Checks, load_csv, quiet, r_available, r_run, read_rows, versions

versions("pandas", "numpy", "scipy", "statsmodels", "matplotlib", "seaborn", "rpy2")
checks = Recorder("24_charts_and_directions")
HAS_R = r_available()
CLOSED = 1e-10
ITERATIVE = 1e-6
if HAS_R:
    print(rd.r_version_of("vegan", "gstat", "extRemes", "drc", "circular", "pROC", "fitdistrplus"))


def drawn(chart):
    """The matplotlib figure of a simplydata chart (its numbers are read from it; nothing is changed)."""
    return quiet(lambda: chart)._render()


def longest_line(ax):
    return max(ax.get_lines(), key=lambda l: len(l.get_xdata()))


def xy(line):
    return np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)


def dots(ax):
    found = [c for c in ax.collections if isinstance(c, PathCollection)]
    return np.vstack([c.get_offsets() for c in found]) if found else np.empty((0, 2))


def band(ax):
    """The lower and upper edge of a shaded band, at each x of its outline."""
    poly = [c for c in ax.collections if type(c).__name__ in ("FillBetweenPolyCollection", "PolyCollection")][0]
    v = poly.get_paths()[0].vertices
    frame = pd.DataFrame(v, columns=["x", "y"]).groupby("x")["y"].agg(["min", "max"])
    return frame.index.to_numpy(), frame["min"].to_numpy(), frame["max"].to_numpy()


def largest_relative(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))


def compare_curve(dataset, method, quantity, ours, reference, reference_name, tolerance=CLOSED):
    """In the reference script: the whole reference curve (the notebook compares the largest relative difference)."""
    checks.vector(dataset, method, quantity, reference, reference_name, tolerance)


# The notebook reads every chart's numbers with chart.get_numbers_drawn() (simplydata main 3fd3b85 and later); the
# data it draws them from are written to references/data/24_charts_and_directions/ here.
mtcars = rd.r_dataset("mtcars")
airquality = rd.r_dataset("airquality")
checks.data("mtcars", mtcars)
checks.data("airquality", airquality)


ds = "mtcars mpg (Henderson & Velleman 1981)"
x = mtcars["mpg"].to_numpy()
n = len(x)
fig = drawn(sdata.chart(mtcars).draw_density(x="mpg"))
gx, gy = xy(longest_line(fig.axes[0]))
h = x.std(ddof=1) * n ** (-1 / 5)
kde = np.array([np.mean(stats.norm.pdf((g - x) / h)) / h for g in gx])
compare_curve(ds, "draw_density", "density", gy, kde, "Gaussian kernel, Scott's bandwidth, NumPy")

fig = drawn(sdata.chart(mtcars).draw_cumulative_share(x="mpg"))
sx, sy = xy(longest_line(fig.axes[0]))
step = pd.DataFrame({"x": sx, "y": sy}).groupby("x")["y"].max()
step = step[np.isin(step.index, x)]
if HAS_R:
    r = r_run("e <- ecdf(x); list(v = e(u))", {"x": x, "u": step.index.to_numpy()})
    compare_curve(ds, "draw_cumulative_share", "share at or below each value", step.to_numpy(), r["v"], "R ecdf")


def filliben(n):
    m = (np.arange(1, n + 1) - 0.3175) / (n + 0.365)
    m[-1] = 0.5 ** (1 / n)
    m[0] = 1 - m[-1]
    return m


for distribution in ("normal", "gamma"):
    options = {} if distribution == "normal" else {"distribution": "gamma"}
    fig = drawn(sdata.chart(mtcars).draw_qq_dots("mpg", **options))
    d = dots(fig.axes[0])
    d = d[np.argsort(d[:, 1])]
    if distribution == "normal":
        theory = stats.norm.ppf(filliben(n))
    else:
        lx = np.log(x)
        k = optimize.brentq(lambda k: np.log(k) - special.digamma(k) - (np.log(x.mean()) - lx.mean()), 1e-3, 1e4,
                            xtol=1e-15, rtol=1e-15)
        theory = stats.gamma.ppf(filliben(n), k, scale=x.mean() / k)
    label = "draw_qq_dots" + ("" if distribution == "normal" else ' (distribution="gamma")')
    compare_curve(ds, label, "dots along (theoretical quantiles)", np.sort(d[:, 0]), theory,
                  "Filliben (1975) positions" + (", gamma fitted by ML in SciPy" if distribution == "gamma" else ""),
                  CLOSED if distribution == "normal" else ITERATIVE)
    compare_curve(ds, label, "dots up (the sorted values)", np.sort(d[:, 1]), np.sort(x), "the data, sorted")
    lx_, ly_ = xy(longest_line(fig.axes[0]))
    slope, intercept = np.polyfit(theory, np.sort(x), 1)
    compare_curve(ds, label, "fitted line", ly_, intercept + slope * lx_, "least squares through the dots, NumPy",
                  CLOSED if distribution == "normal" else ITERATIVE)


ds = "mtcars weight and mpg"
base = sdata.chart(mtcars).draw_dots(x="wt", y="mpg")
fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").draw_smooth_line())
lx_, ly_ = xy(longest_line(fig.axes[0]))
if HAS_R:
    ours = pd.DataFrame({"x": lx_, "y": ly_}).groupby("x")["y"].mean()
    for delta, name in [("0", "R lowess(delta = 0)"), ("0.01 * diff(range(wt))", "R lowess (its default delta, 1 % of the range)")]:
        r = r_run(f"l <- lowess(wt, mpg, delta = {delta}); list(x = l$x, y = l$y)", mtcars[["wt", "mpg"]])
        ref = pd.DataFrame({"x": r["x"], "y": r["y"]}).groupby("x")["y"].mean()
        compare_curve(ds, "draw_smooth_line", "LOWESS line", ours.reindex(ref.index).to_numpy(), ref.to_numpy(), name, ITERATIVE)
    for name_, frame_, xcol, ycol in [("mtcars weight and mpg", mtcars, "wt", "mpg"),
                                      ("Lake Huron levels (no ties)", None, "year", "level")]:
        if frame_ is None:
            frame_ = pd.DataFrame({"year": np.arange(1875, 1973, dtype=float),
                                   "level": rd.r_dataset("LakeHuron").iloc[:, -1].to_numpy()})
        fig = drawn(sdata.chart(frame_).draw_dots(x=xcol, y=ycol).draw_smooth_line(test_method="spline", smoothness=0.01))
        sx, sy = xy(longest_line(fig.axes[0]))
        r = r_run("s <- smooth.spline(a, b, lambda = 0.01); list(y = predict(s, u)$y)", {"a": frame_[xcol], "b": frame_[ycol], "u": sx})
        compare_curve(name_, 'draw_smooth_line (test_method="spline", smoothness=0.01)', "spline", sy, r["y"],
                      "R smooth.spline(lambda = 0.01)", ITERATIVE)
    fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").draw_smooth_line(test_method="loess", with_confidence_band=True))
    sx, sy = xy(longest_line(fig.axes[0]))
    r = r_run("""m <- loess(mpg ~ wt, data = data.frame(wt = wt, mpg = mpg)); p <- predict(m, data.frame(wt = u), se = TRUE)
        q <- qt(0.975, p$df); list(y = p$fit, low = p$fit - q * p$se.fit, high = p$fit + q * p$se.fit)""",
              {"wt": mtcars["wt"], "mpg": mtcars["mpg"], "u": sx})
    compare_curve(ds, 'draw_smooth_line (test_method="loess")', "loess line", sy, r["y"], "R loess (same package)")
    bx, blow, bhigh = band(fig.axes[0])
    r = r_run("""m <- loess(mpg ~ wt, data = data.frame(wt = wt, mpg = mpg)); p <- predict(m, data.frame(wt = u), se = TRUE)
        q <- qt(0.975, p$df); list(low = p$fit - q * p$se.fit, high = p$fit + q * p$se.fit)""",
              {"wt": mtcars["wt"], "mpg": mtcars["mpg"], "u": bx})
    compare_curve(ds, 'draw_smooth_line (with_confidence_band=True)', "band low", blow, r["low"], "R predict(loess, se = TRUE) (same package)")
    compare_curve(ds, 'draw_smooth_line (with_confidence_band=True)', "band high", bhigh, r["high"], "R predict(loess, se = TRUE) (same package)")

huron = rd.r_dataset("LakeHuron")
huron = pd.DataFrame({"year": np.arange(1875, 1875 + len(huron)), "level": huron.iloc[:, -1].to_numpy()})
checks.data("huron", huron)
ds = "Lake Huron levels 1875-1972 (Brockwell & Davis)"
lv = huron["level"].to_numpy()


def smooth_line(**options):
    fig = drawn(sdata.chart(huron).draw_lines(x="year", y="level").draw_smooth_line(**options))
    lines = sorted(fig.axes[0].get_lines(), key=lambda l: len(l.get_xdata()))
    smooth = [l for l in fig.axes[0].get_lines() if not np.allclose(np.asarray(l.get_ydata(), dtype=float), lv, equal_nan=True)
              or len(l.get_ydata()) != len(lv)]
    return xy(smooth[-1] if smooth else lines[-1])


if HAS_R:
    for opts, r_filter in [({"test_method": "moving average", "window": 5}, "rep(1/5, 5)"),
                           ({"test_method": "moving average", "weights": [1, 2, 3, 2, 1]}, "c(1, 2, 3, 2, 1) / 9")]:
        sx, sy = smooth_line(**opts)
        r = np.asarray(r_run(f"list(v = as.numeric(stats::filter(x, {r_filter}, sides = 2)))", {"x": lv})["v"], dtype=float)
        ref = pd.Series(r, index=huron["year"])
        ours = pd.Series(sy, index=sx)
        inner = ref.dropna().index.intersection(ours.dropna().index)
        label = "draw_smooth_line (" + ", ".join(f"{k}={v!r}" for k, v in opts.items()) + ")"
        compare_curve(ds, label, "moving average (inside the ends)", ours.loc[inner].to_numpy(), ref.loc[inner].to_numpy(),
                      "R stats::filter(sides = 2)")
sx, sy = smooth_line(test_method="savitzky-golay", of_degree=2)
parent = signal.savgol_filter(lv, 5, 2)
compare_curve(ds, 'draw_smooth_line (test_method="savitzky-golay", of_degree=2)', "smoothed line", sy, parent,
              "SciPy savgol_filter (window 5, the parent)")
sx, sy = smooth_line(test_method="exponential", newest_weight=0.3)
s = np.empty_like(lv)
s[0] = lv[0]
for t in range(1, len(lv)):
    s[t] = 0.3 * lv[t] + 0.7 * s[t - 1]
compare_curve(ds, 'draw_smooth_line (test_method="exponential", newest_weight=0.3)', "smoothed line", sy, s,
              "the recursion written out, NumPy")


ds = "mtcars mpg by cylinders"
if HAS_R:
    r = r_run("""g <- factor(cyl); m <- tapply(mpg, g, mean); s <- tapply(mpg, g, sd); n <- tapply(mpg, g, length)
        ci <- sapply(split(mpg, g), function(v) t.test(v)$conf.int)
        list(m = as.numeric(m), sd = as.numeric(s), se = as.numeric(s / sqrt(n)), lo = ci[1, ], hi = ci[2, ])""",
              mtcars[["cyl", "mpg"]])
    for kind, low, high in [("standard deviation", np.array(r["m"]) - r["sd"], np.array(r["m"]) + r["sd"]),
                            ("standard error", np.array(r["m"]) - r["se"], np.array(r["m"]) + r["se"]),
                            ("confidence interval", r["lo"], r["hi"])]:
        fig = drawn(sdata.chart(mtcars).draw_dots(x="cyl", y="mpg").draw_error_bars(kind))
        segs = [seg for c in fig.axes[0].collections if isinstance(c, LineCollection) for seg in c.get_segments()]
        ends = sorted((float(s_[:, 0].mean()), float(s_[:, 1].min()), float(s_[:, 1].max())) for s_ in segs if len(s_) == 2
                      and np.isclose(s_[0, 0], s_[1, 0]))
        label = f'draw_error_bars (of_kind="{kind}")'
        compare_curve(ds, label, "bar bottoms", [e[1] for e in ends], low, "R tapply / t.test")
        compare_curve(ds, label, "bar tops", [e[2] for e in ends], high, "R tapply / t.test")

ds = "mtcars weight and mpg"
fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").draw_curve("y = 40 - 5 * x", between=[1, 6], steps=50))
cx, cy = xy(longest_line(fig.axes[0]))
checks.compare(ds, "draw_curve (steps=50)", "points drawn", len(cx), 50, "steps=", 0)
compare_curve(ds, "draw_curve (formula_or_distribution=a formula)", "curve", cy, 40 - 5 * cx, "the formula, NumPy")
if HAS_R:
    for show, f in [("density", "dnorm"), ("cumulative", "pnorm")]:
        dist = sdata.distribution("normal", average=20, standard_deviation=6)
        fig = drawn(sdata.chart(mtcars).draw_curve(dist, between=[5, 35], show=show))
        cx, cy = xy(longest_line(fig.axes[0]))
        r = r_run(f"list(v = {f}(u, 20, 6))", {"u": cx})
        compare_curve(ds, f'draw_curve (a distribution, show="{show}")', show, cy, r["v"], f"R {f}")
    fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").fit_line(with_confidence_band=True))
    fx, fy = xy(longest_line(fig.axes[0]))
    bx, blow, bhigh = band(fig.axes[0])
    r = r_run("""m <- lm(mpg ~ wt); p <- predict(m, data.frame(wt = u)); b <- predict(m, data.frame(wt = v), interval = "confidence")
        list(y = p, low = b[, 2], high = b[, 3])""", {"wt": mtcars["wt"], "mpg": mtcars["mpg"], "u": fx, "v": bx})
    compare_curve(ds, "fit_line", "line", fy, r["y"], "R lm")
    compare_curve(ds, "fit_line (with_confidence_band=True)", "band low", blow, r["low"], 'R predict(interval = "confidence")')
    compare_curve(ds, "fit_line (with_confidence_band=True)", "band high", bhigh, r["high"], 'R predict(interval = "confidence")')
    for degree in (2, 3):
        fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").fit_curve(of_degree=degree))
        fx, fy = xy(longest_line(fig.axes[0]))
        r = r_run(f"list(y = predict(lm(mpg ~ poly(wt, {degree}, raw = TRUE)), data.frame(wt = u)))",
                  {"wt": mtcars["wt"], "mpg": mtcars["mpg"], "u": fx})
        compare_curve(ds, f"fit_curve (of_degree={degree})", "curve", fy, r["y"], "R lm(poly(raw = TRUE))")
    fig = drawn(sdata.chart(mtcars).draw_dots(x="wt", y="mpg").fit_formula("mpg = a * exp(b * wt)", start={"a": 40, "b": -0.3}))
    fx, fy = xy(longest_line(fig.axes[0]))
    r = r_run("""m <- nls(mpg ~ a * exp(b * wt), start = list(a = 40, b = -0.3))
        list(y = predict(m, data.frame(wt = u)))""", {"wt": mtcars["wt"], "mpg": mtcars["mpg"], "u": fx})
    compare_curve(ds, "fit_formula (formula=, start=)", "curve", fy, r["y"], "R nls (its default stopping rule)", ITERATIVE)
    p_, _ = optimize.curve_fit(lambda w, a, b: a * np.exp(b * w), mtcars["wt"], mtcars["mpg"], p0=[40, -0.3], xtol=1e-14,
                               ftol=1e-14)
    compare_curve(ds, "fit_formula (formula=, start=)", "curve", fy, p_[0] * np.exp(p_[1] * fx), "SciPy curve_fit (tight)",
                  ITERATIVE)
    ozone = airquality["Ozone"].dropna().to_numpy()
    fig = drawn(sdata.chart(airquality).draw_histogram(x="Ozone", show="density").fit_distribution("gamma"))
    fx, fy = xy(longest_line(fig.axes[0]))
    r = r_run("""f <- fitdistrplus::fitdist(x, "gamma", control = list(reltol = 1e-14)); list(y = dgamma(u, f$estimate[1], f$estimate[2]))""",
              {"x": ozone, "u": fx})
    compare_curve("airquality ozone", 'fit_distribution (distribution="gamma")', "curve", fy, r["y"],
                  "R dgamma at fitdistrplus fit (reltol 1e-14)", ITERATIVE)
    lo_ = np.log(ozone)
    k_ = optimize.brentq(lambda k: np.log(k) - special.digamma(k) - (np.log(ozone.mean()) - lo_.mean()), 1e-3, 1e4,
                         xtol=1e-15, rtol=1e-15)
    compare_curve("airquality ozone", 'fit_distribution (distribution="gamma")', "curve", fy,
                  stats.gamma.pdf(fx, k_, scale=ozone.mean() / k_), "gamma density at the ML fit (likelihood equation, SciPy)",
                  ITERATIVE)


ds = "mtcars"
fig = drawn(sdata.chart(mtcars).draw_all_dots(x="cyl", y="mpg"))
d = dots(fig.axes[0])
compare_curve(ds, "draw_all_dots", "heights of the dots, sorted", np.sort(d[:, 1]), np.sort(mtcars["mpg"].to_numpy()),
              "the values, sorted")
fig = drawn(sdata.chart(mtcars).draw_count_bars(x="cyl"))
heights = sorted(p.get_height() for p in fig.axes[0].patches)
if HAS_R:
    r = r_run("list(n = as.numeric(table(cyl)))", mtcars[["cyl"]])
    compare_curve(ds, "draw_count_bars", "bar heights, sorted", heights, sorted(r["n"]), "R table")
fig = quiet(sdata.draw_all_pairs, mtcars, columns=["mpg", "wt"])._render()
panels = [ax for ax in fig.axes if len(dots(ax))]
pairs = [dots(ax) for ax in panels]
ok = [float(np.allclose(np.sort(p[:, 0]), np.sort(mtcars[c].to_numpy())) and np.allclose(np.sort(p[:, 1]), np.sort(mtcars[r_].to_numpy())))
      for p, (c, r_) in zip(pairs, [("wt", "mpg"), ("mpg", "wt")])]
for k_, value in enumerate(ok):
    checks.compare(ds, "sdata.draw_all_pairs", f"panel {k_ + 1}: the dots are the two columns' values", value, 1.0,
                   "the data", 0)


creatinine = load_csv("creatinine.csv")
ds = "creatinine, serum and plasma (R mcr)"
fig = drawn(sdata.chart(creatinine).difference_dots(first="serum", second="plasma"))
d = dots(fig.axes[0])
mean_pair = (creatinine["serum"] + creatinine["plasma"]).to_numpy() / 2
diff = (creatinine["serum"] - creatinine["plasma"]).to_numpy()
order = np.lexsort((d[:, 1], d[:, 0]))
ref_order = np.lexsort((diff, mean_pair))
compare_curve(ds, "difference_dots", "dots along (averages)", d[order, 0], mean_pair[ref_order], "(first + second) / 2")
compare_curve(ds, "difference_dots", "dots up (differences)", d[order, 1], diff[ref_order], "first - second")
levels = sorted({round(float(l.get_ydata()[0]), 12) for l in fig.axes[0].get_lines() if len(l.get_ydata()) == 2
                 and np.isclose(l.get_ydata()[0], l.get_ydata()[1])})
ref_levels = sorted([diff.mean() - 1.96 * diff.std(ddof=1), diff.mean(), diff.mean() + 1.96 * diff.std(ddof=1)])
compare_curve(ds, "difference_dots", "the three lines (lower limit, average, upper limit)", levels, ref_levels,
              "Bland & Altman (1986), NumPy")

asah = rd.r_dataset("aSAH", "pROC")
checks.data("asah", asah)
ds = "aSAH (R pROC example; Turck et al. 2010)"
fig = drawn(sdata.chart(asah).trade_off_curve(score="s100b", truth="outcome"))
roc_line = max(fig.axes[0].get_lines(), key=lambda l: len(l.get_xdata()))
rx, ry = xy(roc_line)
score = asah["s100b"].to_numpy()
case = (asah["outcome"] == "Poor").to_numpy()
values = np.unique(score)
cuts = np.concatenate([[-np.inf], (values[:-1] + values[1:]) / 2, [np.inf]])
fpr = np.array([(score[~case] > c).mean() for c in cuts])
tpr = np.array([(score[case] > c).mean() for c in cuts])
# The drawn points and the thresholds' (false alarms, found) pairs, each set without repeats and in order: the same
# two lists mean that every drawn point is a threshold's pair and every threshold's pair is drawn.
reference_points = sorted({(round(a, 12), round(b, 12)) for a, b in zip(fpr, tpr)})
checks.vector(ds, "trade_off_curve", "the drawn points without repeats, in order: false alarms",
              [a for a, _ in reference_points], "the shares counted at each threshold, NumPy", CLOSED)
checks.vector(ds, "trade_off_curve", "the drawn points without repeats, in order: found",
              [b for _, b in reference_points], "the shares counted at each threshold, NumPy", CLOSED)
auc = np.mean((score[case][:, None] > score[~case][None, :]) + 0.5 * (score[case][:, None] == score[~case][None, :]))
legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()][0]
rd.compare_printed(checks, ds, "trade_off_curve", "AUC in the legend", legend.split("AUC ")[1].split(" ")[0], auc,
                   "Mann-Whitney AUC, NumPy")
best = 1 + int(np.argmax(tpr[1:-1] - fpr[1:-1]))
youden = cuts[best]
checks.vector(ds, "trade_off_curve", "best cut-off marked: its point (false alarms, found)", [fpr[best], tpr[best]],
              "Youden's index over midpoints, NumPy", CLOSED)
# The cut-off value written next to the marked point: chart.get_numbers_drawn() gives it in the value column of the
# marked dot (simplydata a5a8ff0 and later).
checks.compare(ds, "trade_off_curve", "best cut-off marked", None, youden, "Youden's index over midpoints, NumPy",
               CLOSED, "relative")
if HAS_R:
    r = r_run('r <- pROC::roc(as.numeric(outcome == "Poor"), s100b, levels = c(0, 1), direction = "<", quiet = TRUE); '
              'b <- pROC::coords(r, "best", best.method = "youden", ret = "threshold", transpose = FALSE); '
              "list(cut = b$threshold[1])", asah[["outcome", "s100b"]])
    checks.compare(ds, "trade_off_curve", "best cut-off marked", None, r["cut"], "R pROC coords (best, Youden)",
                   CLOSED, "relative")


bci = load_csv("bci.csv")
species = list(bci.columns)
plots = bci.iloc[:5].reset_index(drop=True)
plots.insert(0, "plot", [f"plot {i + 1}" for i in range(5)])
checks.data("bci_plots_1_to_5", plots)
ds = "Barro Colorado Island, plots 1-5 (vegan BCI)"
try:
    fig = drawn(sdata.chart(plots).draw_richness_curves(columns=species, names="plot"))
    curves = sorted(fig.axes[0].get_lines(), key=lambda l: -len(l.get_xdata()))[:5]
    for line in curves:
        cx_, cy_ = xy(line)
        total = int(round(cx_.max()))
        row = plots[species].to_numpy(dtype=float)
        which = int(np.argmin(np.abs(row.sum(axis=1) - total)))
        counts = row[which][row[which] > 0]
        N = counts.sum()
        def lchoose(a, b):
            return special.gammaln(a + 1) - special.gammaln(b + 1) - special.gammaln(a - b + 1)
        expected = [np.sum(1 - np.exp(np.where(N - counts >= m, lchoose(N - counts, m), -np.inf) - lchoose(N, m)))
                    for m in cx_]
        compare_curve(ds, "draw_richness_curves", f"curve of plot {which + 1}", cy_, expected, "Hurlbert (1971), NumPy")
except Exception as problem:  # noqa: BLE001
    checks.skip(ds, "draw_richness_curves", "curves", "Hurlbert (1971)", f"the chart could not be drawn: {str(problem)[:120]}")

rain = load_csv("coles_rain.csv")
ds = "daily rainfall, south-west England (Coles 2001)"
try:
    fig = drawn(sdata.chart(rain).draw_threshold_check(x="rain", between=[10, 30]))
    top = fig.axes[0].child_axes[0]   # the two panels are drawn as panels inside the chart
    lines = [xy(l) for l in top.get_lines() if len(l.get_xdata()) > 2]
    middle = sorted(lines, key=lambda p: np.mean(p[1]))[1]
    low_l, high_l = sorted(lines, key=lambda p: np.mean(p[1]))[0], sorted(lines, key=lambda p: np.mean(p[1]))[2]
    tx, ty = middle
    r_ = rain["rain"].to_numpy()
    excess = [r_[r_ > u] - u for u in tx]
    me = np.array([e.mean() for e in excess])
    half = np.array([stats.norm.ppf(0.975) * e.std(ddof=1) / np.sqrt(len(e)) for e in excess])
    compare_curve(ds, "draw_threshold_check", "mean excess", ty, me, "mean of x - u above u, NumPy")
    compare_curve(ds, "draw_threshold_check", "mean excess, 95 % interval low", low_l[1], me - half,
                  "mean - 1.96 sd / sqrt(n), NumPy")
    compare_curve(ds, "draw_threshold_check", "mean excess, 95 % interval high", high_l[1], me + half,
                  "mean + 1.96 sd / sqrt(n), NumPy")
except Exception as problem:  # noqa: BLE001
    checks.skip(ds, "draw_threshold_check", "mean excess", "NumPy", f"the chart could not be drawn: {str(problem)[:120]}")

meuse = load_csv("meuse.csv")
ds = "Meuse zinc (sp / gstat)"
xy_ = meuse[["x", "y"]].to_numpy(dtype=float)
z = meuse["zinc"].to_numpy(dtype=float)
cutoff = np.hypot(np.ptp(xy_[:, 0]), np.ptp(xy_[:, 1])) / 3
width = cutoff / 15
i, j = np.triu_indices(len(z), 1)
dist = np.hypot(*(xy_[i] - xy_[j]).T)
sq = (z[i] - z[j]) ** 2
keep = dist <= cutoff
bins = np.minimum((dist[keep] / width).astype(int), 14)
emp = pd.DataFrame({"bin": bins, "d": dist[keep], "g": sq[keep] / 2}).groupby("bin").agg(d=("d", "mean"), g=("g", "mean"),
                                                                                        n=("g", "size"))
for kind in ("spherical", "exponential"):
    found_rows = read_rows(quiet(sdata.table(meuse).find_variogram_of, "zinc", x="x", y="y", of_kind=kind))
    found = {c: [row[c] for row in found_rows] for c in found_rows[0]}   # each column, read with table.get_row
    print(found_rows[:3])
    label = f'find_variogram_of (of_kind="{kind}")'
    compare_curve(ds, label, "distance of each bin", np.asarray(found["distance"], dtype=float), emp["d"].to_numpy(),
                  "gstat's bins written out, NumPy")
    compare_curve(ds, label, "semivariance of each bin", np.asarray(found["semivariance"], dtype=float), emp["g"].to_numpy(),
                  "gstat's bins written out, NumPy")
    pcol = [c for c in found if c == "pairs"]
    if pcol:
        compare_curve(ds, label, "pairs in each bin", np.asarray(found["pairs"], dtype=float), emp["n"].to_numpy(), "counted, NumPy")
    if HAS_R:
        model = {"spherical": "Sph", "exponential": "Exp"}[kind]
        r = r_run(f"""suppressPackageStartupMessages(library(gstat)); library(sp)
            d <- data.frame(x = x, y = y, zinc = zinc); coordinates(d) <- ~x + y
            v <- variogram(zinc ~ 1, d); f <- fit.variogram(v, vgm("{model}"))
            list(nugget = f$psill[1], psill = f$psill[2], range = f$range[2])""", meuse[["x", "y", "zinc"]])
        for quantity, column in [("nugget", [c for c in found if "nugget" in c]),
                                 ("psill", [c for c in found if "partial sill" in c]),
                                 ("range", [c for c in found if c.startswith("fitted range") or c == "range"])]:
            if column:
                checks.compare(ds, label, f"fitted {quantity}", found[column[0]][0], r[quantity],
                               "R gstat fit.variogram (same package)", ITERATIVE, "relative")
    fig = drawn(sdata.chart(meuse).draw_variogram("zinc", x="x", y="y", of_kind=kind))
    d = dots(fig.axes[0])
    d = d[np.argsort(d[:, 0])]
    compare_curve(ds, f'draw_variogram (of_kind="{kind}")', "dots up (semivariances)", d[:, 1], emp["g"].to_numpy(),
                  "gstat's bins written out, NumPy")

rye = load_csv("ryegrass.csv")
ds = "ryegrass (drc; Inderjit et al. 2002)"


def ll4(x, b, c, d, e):
    with np.errstate(divide="ignore"):
        return c + (d - c) / (1 + np.exp(b * (np.log(x) - np.log(e))))


try:
    fig = drawn(sdata.chart(rye).draw_dots(x="conc", y="root_length").fit_dose_response())
    fx, fy = xy(longest_line(fig.axes[0]))
    params, _ = optimize.curve_fit(ll4, rye["conc"].replace(0, 1e-12), rye["root_length"], p0=[2, 0.5, 8, 3],
                                   xtol=1e-14, ftol=1e-14, maxfev=100000)
    compare_curve(ds, "fit_dose_response", "curve", fy[fx > 0], ll4(fx[fx > 0], *params), "SciPy curve_fit (tight)",
                  ITERATIVE)
except Exception as problem:  # noqa: BLE001
    checks.skip(ds, "fit_dose_response", "curve", "SciPy", f"the chart could not be drawn: {str(problem)[:120]}")

usarrests = rd.r_dataset("USArrests")
checks.data("usarrests", usarrests)
ds = "USArrests (McNeil 1977)"


def join_heights(fig):
    heights = []
    for line in fig.axes[0].get_lines():
        y = np.asarray(line.get_ydata(), dtype=float)
        if len(y) == 4:   # each join is drawn as a bracket: up, across, down
            heights.append(y[1])
    for c in fig.axes[0].collections:
        if isinstance(c, LineCollection):
            for seg in c.get_segments():
                if len(seg) == 4:
                    heights.append(seg[1, 1])
    return np.sort(heights)


columns = ["Murder", "Assault", "UrbanPop", "Rape"]
fig = drawn(sdata.chart(usarrests).draw_cluster_tree(columns=columns, names="row name"))
if HAS_R:
    r = r_run("list(h = sort(hclust(dist(scale(cbind(a, b, c, d))), method = 'ward.D2')$height))",
              {"a": usarrests["Murder"], "b": usarrests["Assault"], "c": usarrests["UrbanPop"], "d": usarrests["Rape"]})
    compare_curve(ds, "draw_cluster_tree", "join heights (Ward, scaled columns)", join_heights(fig), r["h"],
                  'R hclust(method = "ward.D2") on scaled columns')
    small = bci.iloc[:12]
    sites = small.copy()
    sites.insert(0, "plot", [f"p{i + 1}" for i in range(len(small))])
    checks.data("bci_plots_1_to_12", sites)
    cols = {f"s{k}": small[c] for k, c in enumerate(species)}
    for linkage, r_method in [(None, "average"), ("complete", "complete")]:
        opts = {"linkage": linkage} if linkage else {}
        fig = drawn(sdata.chart(sites).draw_cluster_tree(columns=species, measure="bray-curtis", names="plot", **opts))
        r = r_run(f"""m <- as.matrix(data.frame({", ".join(cols)})); list(h = sort(hclust(vegan::vegdist(m, "bray"), method = "{r_method}")$height))""",
                  cols)
        label = 'draw_cluster_tree (measure="bray-curtis"' + (f', linkage="{linkage}")' if linkage else ")")
        compare_curve("Barro Colorado Island, plots 1-12 (vegan BCI)", label, f"join heights ({r_method})",
                      join_heights(fig), r["h"], f'R hclust("{r_method}") on vegan vegdist("bray")')


pigeons = rd.r_dataset("pigeons", "circular")
print(pigeons.head())
direction = [c for c in pigeons.columns if pigeons[c].dtype.kind == "f"][0]
group = [c for c in pigeons.columns if pigeons[c].dtype.kind != "f" and c != "row name"][0]
ds = "homing pigeons (R circular pigeons)"
found = read_rows(quiet(sdata.table(pigeons).find_average_direction_of, direction, full_turn=360, for_each_group_in=group))
print(found)
if HAS_R:
    for line in found:
        g = line[group]
        values = pigeons.loc[pigeons[group] == g, direction].to_numpy()
        r = r_run("""suppressPackageStartupMessages(library(circular))
            x <- circular(v, units = "degrees", template = "none", modulo = "2pi")
            m <- as.numeric(mean(x)) %% 360
            list(mean = m, rho = rho.circular(x), sd = as.numeric(sd.circular(x)) * 180 / pi)""", {"v": values})
        avg_col = [c for c in line if c.startswith("average direction")][0]
        rho_col = [c for c in line if "resultant length" in c][0]
        sd_col = [c for c in line if "circular standard deviation" in c][0]
        label = "find_average_direction_of (for_each_group_in=)"
        checks.compare(ds, label, f"average direction, {g}", line[avg_col], r["mean"], "R circular mean.circular", CLOSED, "relative")
        checks.compare(ds, label, f"resultant length, {g}", line[rho_col], r["rho"], "R circular rho.circular", CLOSED, "relative")
        checks.compare(ds, label, f"circular standard deviation, {g}", line[sd_col], r["sd"], "R circular sd.circular", CLOSED,
                       "relative")


checks.data("pigeons", pigeons)
checks.save()
