"""Reference values for notebooks/13_circular_spatial_extremes.ipynb, saved to references/13_circular_spatial_extremes.csv.

Printed first: Coles (2001, An Introduction to Statistical Modeling of Extreme Values) for the GEV fits of Port Pirie
(3.4.1) and Fremantle (6.3.1), the generalised Pareto fit of the daily rain (4.4.1) and their return levels; the
spherical variogram of log zinc printed in the gstat vignette. Each within half a unit of the last printed digit.
Elsewhere the original library is run once: SciPy (circmean, circstd, directional_stats; genextreme and genpareto
fitted to a tight stopping rule) where SciPy is the parent; NumPy written out from the definition for Moran's I and
the empirical variogram; R circular, ape, spdep, gstat and extRemes.

Tolerances: relative 1e-10 for closed formulas, relative 1e-6 for fits by iteration; the exact ends of the profile
likelihood within one step of extRemes's grid.

Data written for the notebook (references/data/13_circular_spatial_extremes/): the repository's data files, with the
three samples of Fisher B.10 as separate tables, the Meuse metals on the log scale, the rain excesses above 30 mm and
the daily rain with a date column.

    python scripts/references/13_circular_spatial_extremes.py      (from the repository root; needs R)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402  (also puts notebooks/ on the path for checks)

import math  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import optimize, stats  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "13_circular_spatial_extremes"
checks = Recorder(STEM, ("circular", "spdep", "ape", "gstat", "extRemes"))
HAS_R = r_available()
CLOSED = 1e-10
ITERATIVE = 1e-6
COLES = ("Coles (2001), An Introduction to Statistical Modeling of Extreme Values, Springer",
         "https://doi.org/10.1007/978-1-4471-3675-0")
GSTAT = ("Pebesma, gstat: spatial and spatio-temporal geostatistical modelling (gstat vignette)",
         "https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf")


def ref(dataset, method, quantity, value, source, tolerance=CLOSED, how="relative"):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, how)


def printed(dataset, method, quantity, shown, source, half_unit, where):
    """A printed number with the half-unit tolerance the notebook has always used."""
    checks.compare(dataset, method, quantity, None, float(shown), source, half_unit, "absolute")
    checks.rows[-1].update({"digits as printed": shown, "source title": where[0], "url": where[1], "kind": "printed"})


# ---- directions
turtles = load_csv("turtles.csv")
icu = load_csv("icu_arrivals.csv")
wind = load_csv("wind_directions.csv")
fisher_b10 = load_csv("fisher_b10.csv")
for name, frame in (("turtles", turtles), ("icu_arrivals", icu), ("wind_directions", wind), ("fisher_b10", fisher_b10)):
    checks.data(name, frame)
TWO_PI = 2 * math.pi
R_UNITS = {360: "degrees", 24: "hours", TWO_PI: "radians"}
DIRECTION_SETS = [("turtles (circular)", turtles, "arrival", 360),
                  ("ICU arrivals (circular fisherB1)", icu, "hour", 24),
                  ("wind (circular)", wind, "direction_radians", "radians")]
for s in ("set1", "set2", "set3"):
    part = fisher_b10[fisher_b10["set"] == s].reset_index(drop=True)
    checks.data(f"fisher_b10_{s}", part)
    DIRECTION_SETS.append((f"Fisher B.10 {s} (circular fisherB10)", part, "direction", 360))


def scipy_direction(values, turn):
    t = np.asarray(values, dtype=float) / turn * TWO_PI
    return {"average direction": stats.circmean(values, high=turn, low=0) % turn,
            "resultant length": stats.directional_stats(np.column_stack([np.cos(t), np.sin(t)])).mean_resultant_length,
            "circular standard deviation": stats.circstd(values, high=turn, low=0)}


def r_direction(values, turn):
    code = f'''a <- circular::circular(v, units = "{R_UNITS[turn]}", template = "none")
               list(m = (as.numeric(circular::mean.circular(a)) + {turn!r}) %% {turn!r},
                    rho = circular::rho.circular(a), sd = as.numeric(circular::sd.circular(a)) * {turn!r} / (2 * pi))'''
    r = r_run(code, {"v": values}, packages=("circular",))
    return {"average direction": r["m"], "resultant length": r["rho"], "circular standard deviation": r["sd"]}


for dataset, data, column, full_turn in DIRECTION_SETS:
    turn = TWO_PI if full_turn == "radians" else full_turn
    references = [("SciPy circmean, directional_stats, circstd", scipy_direction(data[column], turn))]
    if HAS_R:
        references.append(("R circular", r_direction(data[column], turn)))
    for reference_name, reference in references:
        for quantity, value in reference.items():
            ref(dataset, "find_average_direction_of", quantity, value, reference_name)
for s in ("set1", "set2", "set3"):
    reference = scipy_direction(fisher_b10.loc[fisher_b10["set"] == s, "direction"], 360)
    for quantity, value in reference.items():
        ref(f"Fisher B.10 {s} (circular fisherB10)", "find_average_direction_of (for_each_group_in=)", quantity, value,
            "SciPy circmean, directional_stats, circstd")

# ---- Rayleigh test
method = "compare_shape_to (Rayleigh)"
for dataset, data, column, full_turn in DIRECTION_SETS:
    turn = TWO_PI if full_turn == "radians" else full_turn
    ref(dataset, method, "resultant length", scipy_direction(data[column], turn)["resultant length"],
        "SciPy directional_stats")
    if HAS_R:
        r = r_run(f'''t <- circular::rayleigh.test(circular::circular(v, units = "{R_UNITS[turn]}", template = "none"))
                      list(statistic = unname(t$statistic), p = t$p.value)''', {"v": data[column]}, packages=("circular",))
        ref(dataset, method, "resultant length", r["statistic"], "R circular rayleigh.test")
        ref(dataset, method, "p-value", r["p"], "R circular rayleigh.test")

# ---- Moran's I
columbus = load_csv("columbus.csv")
meuse = load_csv("meuse.csv")
for metal in ("zinc", "copper", "cadmium"):
    meuse[f"log_{metal}"] = np.log(meuse[metal])
checks.data("columbus", columbus)
checks.data("meuse", meuse)


def numpy_moran(z, x, y, nearest=None):
    xy = np.column_stack([x, y]).astype(float)
    d = np.sqrt(((xy[:, None, :] - xy[None, :, :]) ** 2).sum(axis=2))
    if nearest is None:
        with np.errstate(divide="ignore"):
            w = np.where(d > 0, 1 / d, 0.0)
    else:
        w = np.zeros_like(d)
        for i in range(len(d)):
            order = [j for j in np.argsort(d[i], kind="stable") if j != i][:nearest]
            w[i, order] = 1.0
    w = w / w.sum(axis=1, keepdims=True)
    z = np.asarray(z, dtype=float) - np.mean(z)
    return len(z) / w.sum() * (z @ w @ z) / (z @ z)


MORAN_SETS = [("Columbus (spData)", columbus, c) for c in ("crime", "house_value", "income")]
MORAN_SETS += [("Meuse (sp)", meuse, "log_zinc")]
method = "find_autocorrelation_of (x=, y=)"
for dataset, data, column in MORAN_SETS:
    n = len(data)
    for nearest in (None, 4, 8):
        weights = "1 / distance" if nearest is None else f"nearest={nearest}"
        quantity = f"{column}: Moran's I ({weights})"
        ref(dataset, method, quantity, numpy_moran(data[column], data["x"], data["y"], nearest), "NumPy, the formula")
        ref(dataset, method, f"{column}: value with no pattern ({weights})", -1 / (n - 1), "-1 / (n - 1)")
        if not HAS_R:
            continue
        if nearest is None:
            r = r_run('''w <- 1 / as.matrix(stats::dist(cbind(x, y))); diag(w) <- 0
                         pts <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x", "y"))
                         lw <- spdep::nb2listwdist(spdep::dnearneigh(pts, 0, Inf), pts, type = "idw", alpha = 1, style = "W")
                         set.seed(1); mc <- spdep::moran.mc(z, lw, nsim = 999)
                         list(I = ape::Moran.I(z, w)$observed, p = mc$p.value)''',
                      {"x": data["x"], "y": data["y"], "z": data[column]}, packages=("spdep", "sf", "ape"))
            name = "R ape Moran.I"
        else:
            r = r_run(f'''lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(cbind(x, y), k = {nearest})), style = "W")
                          set.seed(1); mc <- spdep::moran.mc(z, lw, nsim = 999)
                          list(I = spdep::moran(z, lw, length(z), spdep::Szero(lw))$I, p = mc$p.value)''',
                      {"x": data["x"], "y": data["y"], "z": data[column]}, packages=("spdep",))
            name = "R spdep moran (knn2nb)"
        ref(dataset, method, quantity, r["I"], name)
        ref(dataset, method, f"{column}: p-value, 999 shuffles ({weights})", r["p"], "R spdep moran.mc (set.seed(1))")


# ---- variogram
def numpy_variogram(z, x, y, bins=15):
    x, y, z = (np.asarray(v, dtype=float) for v in (x, y, z))
    cutoff = math.hypot(x.max() - x.min(), y.max() - y.min()) / 3
    i, j = np.triu_indices(len(z), k=1)
    d = np.hypot(x[i] - x[j], y[i] - y[j])
    half_squares = (z[i] - z[j]) ** 2 / 2
    width = cutoff / bins
    which = np.ceil(d / width).astype(int) - 1  # bin k holds k*width < d <= (k+1)*width
    out = []
    for k in range(bins):
        keep = which == k
        out.append((keep.sum(), d[keep].mean(), half_squares[keep].mean()))
    return pd.DataFrame(out, columns=["pairs", "distance", "semivariance"])


VARIOGRAMS = [("log_zinc", "spherical", "Sph"), ("log_copper", "exponential", "Exp"), ("log_cadmium", "spherical", "Sph")]
for column, kind, model in VARIOGRAMS:
    dataset, method = "Meuse (sp)", f"find_variogram_of ({kind})"
    mine = numpy_variogram(meuse[column], meuse["x"], meuse["y"])
    for k in range(len(mine)):
        for quantity in ("pairs", "distance", "semivariance"):
            ref(dataset, method, f"{column}: bin {k + 1} {quantity}", mine[quantity][k], "NumPy, the definition")
    if column == "log_zinc":
        for quantity, shown, half_unit in (("nugget", "0.0507", 0.00005), ("partial sill", "0.5906", 0.00005),
                                           ("range", "897", 0.5)):
            printed(dataset, method, f"{column}: {quantity}", shown, "gstat vignette, printed", half_unit, GSTAT)
    if HAS_R:
        r = r_run(f'''v <- gstat::variogram(z ~ 1, locations = ~ x + y, data = data.frame(x = x, y = y, z = z))
                      f <- gstat::fit.variogram(v, gstat::vgm("{model}"))
                      list(np = v$np, dist = v$dist, gamma = v$gamma, nugget = f$psill[1], psill = f$psill[2],
                           range = f$range[2])''', {"x": meuse["x"], "y": meuse["y"], "z": meuse[column]},
                  packages=("gstat", "sf"))
        for k in range(len(r["gamma"])):
            for quantity, key in (("pairs", "np"), ("distance", "dist"), ("semivariance", "gamma")):
                ref(dataset, method, f"{column}: bin {k + 1} {quantity}", r[key][k], "R gstat variogram")
        for quantity, key in (("nugget", "nugget"), ("partial sill", "psill"), ("range", "range")):
            ref(dataset, method, f"{column}: {quantity}", r[key], "R gstat fit.variogram", ITERATIVE)

# ---- GEV and generalised Pareto fits
portpirie = load_csv("portpirie.csv")
fremantle = load_csv("fremantle.csv")
rain = load_csv("coles_rain.csv")
excess = pd.DataFrame({"excess": rain.loc[rain["rain"] > 30, "rain"] - 30})
daily = rain.assign(date=pd.date_range("1914-01-01", periods=len(rain), freq="D").strftime("%Y-%m-%d"))
for name, frame in (("portpirie", portpirie), ("fremantle", fremantle), ("rain_excess", excess), ("rain_daily", daily)):
    checks.data(name, frame)


def tight(func, x0, args=(), disp=0):
    return optimize.fmin(func, x0, args=args, xtol=1e-12, ftol=1e-12, maxiter=200000, maxfun=200000, disp=0)


def r_fevd(values, kind, threshold=None):
    extra = f", threshold = {threshold!r}" if threshold is not None else ""
    code = f'''f <- extRemes::fevd(v, type = "{kind}", method = "MLE"{extra},
                                 optim.args = list(method = "BFGS", control = list(reltol = 1e-14, maxit = 100000)))
               list(par = unname(f$results$par), se = unname(sqrt(diag(solve(f$results$hessian)))))'''
    return r_run(code, {"v": values}, packages=("extRemes",))


PRINTED_GEV = {"Port Pirie (ismev)": {"location": ("3.87", 0.005), "scale": ("0.198", 0.0005),
                                      "shape": ("-0.050", 0.0005), "standard error of location": ("0.028", 0.0005),
                                      "standard error of scale": ("0.020", 0.0005),
                                      "standard error of shape": ("0.098", 0.0005)},
               "Fremantle (ismev)": {"location": ("1.482", 0.0005), "scale": ("0.141", 0.0005),
                                     "shape": ("-0.217", 0.0005)}}
SOURCE = {"Port Pirie (ismev)": "Coles (2001) 3.4.1, printed", "Fremantle (ismev)": "Coles (2001) 6.3.1, printed"}
method = 'fit_distribution_to (distribution="gev")'
for dataset, data in (("Port Pirie (ismev)", portpirie), ("Fremantle (ismev)", fremantle)):
    for quantity, (shown, half_unit) in PRINTED_GEV[dataset].items():
        printed(dataset, method, quantity, shown, SOURCE[dataset], half_unit, COLES)
    c, loc, scale = stats.genextreme.fit(data["sea_level"].to_numpy(dtype=float), optimizer=tight)
    for quantity, value in (("location", loc), ("scale", scale), ("shape", -c)):
        ref(dataset, method, quantity, value, "SciPy genextreme.fit (tight)", ITERATIVE)
    if HAS_R:
        r = r_fevd(data["sea_level"], "GEV")
        for k, quantity in enumerate(("location", "scale", "shape")):
            ref(dataset, method, quantity, r["par"][k], "R extRemes fevd", ITERATIVE)
            ref(dataset, method, f"standard error of {quantity}", r["se"][k], "R extRemes fevd", ITERATIVE)

dataset, method = "Coles rain (ismev)", 'fit_distribution_to (distribution="generalised pareto")'
for quantity, shown, half_unit in (("scale", "7.44", 0.005), ("shape", "0.184", 0.0005)):
    printed(dataset, method, quantity, shown, "Coles (2001) 4.4.1, printed", half_unit, COLES)
c, loc, scale = stats.genpareto.fit(excess["excess"].to_numpy(dtype=float), floc=0, optimizer=tight)
for quantity, value in (("scale", scale), ("shape", c)):
    ref(dataset, method, quantity, value, "SciPy genpareto.fit (floc=0, tight)", ITERATIVE)
if HAS_R:
    r = r_fevd(rain["rain"], "GP", threshold=30.0)
    for k, quantity in enumerate(("scale", "shape")):
        ref(dataset, method, quantity, r["par"][k], "R extRemes fevd (threshold 30)", ITERATIVE)

# ---- return levels
PROFILE_R = '''
nll_gev <- function(z, s, xi, p, x) { mu <- z + s / xi * (1 - (-log(1 - p))^(-xi))
  t <- 1 + xi * (x - mu) / s; if (s <= 0 || any(t <= 0)) return(1e10)
  sum(log(s) + (1 + 1 / xi) * log(t) + t^(-1 / xi)) }
nll_gp <- function(z, xi, m, y, u) { s <- xi * (z - u) / (m^xi - 1)
  t <- 1 + xi * y / s; if (s <= 0 || any(t <= 0)) return(1e10)
  sum(log(s) + (1 + 1 / xi) * log(t)) }
profile_ends <- function(prof, best, zhat, lo, hi) {
  g <- function(z) 2 * (prof(z) - best) - qchisq(0.95, 1)
  c(uniroot(g, c(lo, zhat), tol = 1e-12)$root, uniroot(g, c(zhat, hi), tol = 1e-12)$root) }
'''


def r_return_levels(x, periods, threshold=None):
    if threshold is None:
        fit = 'extRemes::fevd(x, type = "GEV", method = "MLE"'
        level = "p <- 1 / T; zhat <- par[1] - par[2] / par[3] * (1 - (-log(1 - p))^(-par[3]))"
        profile = ("prof <- function(z) optim(par[2:3], function(q) nll_gev(z, q[1], q[2], p, x), "
                   "control = list(reltol = 1e-14, maxit = 100000))$value; lowest <- rl[1] - w")
    else:
        u = float(threshold)
        fit = f'extRemes::fevd(x, threshold = {u!r}, type = "GP", time.units = "days", method = "MLE"'
        level = (f"y <- x[x > {u!r}] - {u!r}; m <- T * 365.25 * mean(x > {u!r}); "
                 f"zhat <- {u!r} + par[1] / par[2] * (m^par[2] - 1)")
        profile = (f"prof <- function(z) optimize(function(xi) nll_gp(z, xi, m, y, {u!r}), c(-0.5, 1.5), "
                   f"tol = 1e-12)$objective; lowest <- max({u!r} + 1e-6, rl[1] - w)")
    code = PROFILE_R + f'''
f0 <- {fit})
f <- {fit}, optim.args = list(method = "BFGS", control = list(reltol = 1e-14, maxit = 100000)))
par <- f$results$par; best <- f$results$value
out <- list()
for (T in c({", ".join(map(str, periods))})) {{
  {level}
  rl <- as.numeric(extRemes::return.level(f, return.period = T, do.ci = TRUE))
  w <- rl[3] - rl[1]; xr <- c(rl[1], rl[3] + 2 * w)
  grid <- as.numeric(extRemes::ci.fevd(f, type = "return.level", return.period = T, method = "proflik",
                                       xrange = xr, nint = 200, verbose = FALSE))
  {profile}
  ends <- profile_ends(prof, best, zhat, lowest, rl[3] + 6 * w)
  default <- as.numeric(extRemes::return.level(f0, return.period = T))
  out[[paste0("T", T)]] <- c(zhat, rl, grid[c(1, 3)], ends, (xr[2] - xr[1]) / 200, default)
}}
out'''
    return r_run(code, {"x": x}, packages=("extRemes",))


def return_level_references(dataset, method, r, periods):
    for T in periods:
        zhat, n_low, level, n_high, g_low, g_high, e_low, e_high, step, default = r[f"T{T}"]
        ref(dataset, method, f"{T}-year level", default, "R extRemes return.level (default stopping rule)", ITERATIVE)
        ref(dataset, method, f"{T}-year level", zhat, "the formula on the R fit", ITERATIVE)
        for quantity, value in (("level", level), ("normal 95 % CI low", n_low), ("normal 95 % CI high", n_high)):
            ref(dataset, method, f"{T}-year {quantity}", value, "R extRemes return.level", ITERATIVE)
        for quantity, grid_value, exact in (("95 % CI low", g_low, e_low), ("95 % CI high", g_high, e_high)):
            ref(dataset, method, f"{T}-year profile {quantity}", grid_value, "R extRemes ci.fevd (same grid)", ITERATIVE)
            ref(dataset, method, f"{T}-year profile {quantity}", exact, "exact profile (one grid step)", step, "absolute")


dataset, method = "Port Pirie (ismev)", "find_return_levels_of"
for T, quantity, shown in ((10, "level", "4.30"), (100, "level", "4.69"), (100, "normal 95 % CI low", "4.38"),
                           (100, "normal 95 % CI high", "5.00")):
    printed(dataset, method, f"{T}-year {quantity}", shown, "Coles (2001) 3.4.1, printed", 0.005, COLES)
if HAS_R:
    return_level_references(dataset, method, r_return_levels(portpirie["sea_level"], [10, 100]), [10, 100])

dataset, method = "Coles rain (ismev)", "find_return_levels_of (values_above=30)"
printed(dataset, method, "100-year level", "106.3", "Coles (2001) 4.4.1, printed", 0.05, COLES)
if HAS_R:
    return_level_references(dataset, method, r_return_levels(rain["rain"], [10, 100], threshold=30), [10, 100])

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
