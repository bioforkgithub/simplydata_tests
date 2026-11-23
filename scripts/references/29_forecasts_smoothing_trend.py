"""Reference values for notebooks/29_forecasts_smoothing_trend.ipynb.

Printed first: the trend package vignette (Mann-Kendall, Sen's slope, seasonal Mann-Kendall and seasonal slope) and
Cobb (1978) for the Nile change point (1898). Where nothing is printed, the help-page example of the R function
simplydata calls or names is run once (HoltWinters, forecast, arima, filter, runmed, signal::sgolayfilt, lowess,
smooth.spline, loess, stl, decompose, lm, acf/pacf/Box.test, spdep moran.mc, urca, tseries, trend, changepoint) and
saved with the link to its page; SciPy for the Savitzky-Golay filter and detrend that simplydata calls.
"""

import numpy as np
import pandas as pd
from scipy import signal

from _tools import References, VIGNETTE, r_available, r_dataset, r_run, r_series, r_version, rman

refs = References("29_forecasts_smoothing_trend")
HAS_R = r_available()
CLOSED, ITERATIVE = 1e-10, 1e-6


def series(name, package="datasets"):
    s = r_series(name, package)
    s["step"] = np.arange(len(s))
    refs.data(name, s)
    return s


co2, air, nile, lh, usacc, presidents = (series(n) for n in ("co2", "AirPassengers", "Nile", "lh", "USAccDeaths", "presidents"))
nottem, lakehuron = series("nottem"), series("LakeHuron")
maxau = r_dataset("maxau", "trend")
maxau["step"] = np.arange(len(maxau))
refs.data("maxau", maxau)
quakes = r_dataset("quakes")
quakes["longitude"] = np.where(quakes["long"] > 180, quakes["long"] - 360, quakes["long"])
quakes = quakes.drop_duplicates(["longitude", "lat"]).reset_index(drop=True)
refs.data("quakes", quakes)

# 1-2. forecasts (R forecast / stats, the help-page series)
ARIMA = "m <- arima(y, order = c({o}){s}, method = 'ML'); p <- predict(m, n.ahead = H); f <- list(mean = p$pred, lower = p$pred - qnorm(0.5 + L / 200) * p$se, upper = p$pred + qnorm(0.5 + L / 200) * p$se)"
CASES = [("co2", co2, "repeats_every=12", "f <- forecast::forecast(HoltWinters(y), h = H, level = L)", "HoltWinters", "stats"),
         ("AirPassengers", air, "repeats_every=12, season='multiplied'", 'f <- forecast::forecast(HoltWinters(y, seasonal = "multiplicative"), h = H, level = L)', "HoltWinters", "stats"),
         ("Nile", nile, "trend=False", "f <- forecast::forecast(HoltWinters(y, beta = FALSE, gamma = FALSE), h = H, level = L)", "HoltWinters", "stats"),
         ("co2", co2, "repeats_every=12, test_method='same season last time'", "f <- forecast::snaive(y, h = H, level = L)", "naive", "forecast"),
         ("Nile", nile, "test_method='last value'", "f <- forecast::naive(y, h = H, level = L)", "naive", "forecast"),
         ("Nile", nile, "test_method='average'", "f <- forecast::meanf(y, h = H, level = L)", "meanf", "forecast"),
         ("Nile", nile, "test_method='straight on'", "f <- forecast::rwf(y, drift = TRUE, h = H, level = L)", "rwf", "forecast"),
         ("lh", lh, "test_method='arima', order=(1, 0, 0)", ARIMA.format(o="1, 0, 0", s=""), "arima", "stats"),
         ("USAccDeaths", usacc, "test_method='arima', order=(0, 1, 1), seasonal_order=(0, 1, 1), repeats_every=12",
          ARIMA.format(o="0, 1, 1", s=", seasonal = list(order = c(0, 1, 1))"), "arima", "stats"),
         ("presidents", presidents, "test_method='arima', order=(1, 0, 0)", ARIMA.format(o="1, 0, 0", s=""), "arima", "stats")]
if HAS_R:
    for name, data, options, code, topic, package in CASES:
        for level, steps in ((95, 12), (80, 6)):
            freq = int(data.attrs.get("frequency", 1)) if hasattr(data, "attrs") else 1
            freq = 12 if name in ("co2", "AirPassengers", "USAccDeaths", "nottem") else (4 if name == "presidents" else 1)
            r = r_run(f'''y <- ts(value, frequency = {freq}); H <- {steps}; L <- {level}
                          {code}
                          list(forecast = as.numeric(f$mean), low = as.numeric(f$lower), high = as.numeric(f$upper))''', data, packages=("forecast",))
            check = f"forecast ({options}, steps_ahead={steps}, interval={level})"
            for key in ("forecast", "low", "high"):
                refs.run_vector(name, check, key, r[key], f"R {topic} + forecast", rman(package, topic), package, ITERATIVE)
                if key == "forecast" and level == 95 and topic == "HoltWinters" and name in ("co2", "Nile"):
                    refs.run_vector(name, check.replace("forecast (", "nums.forecast ("), key, r[key], f"R {topic} + forecast", rman(package, topic), package, ITERATIVE)
    r = r_run("m <- HoltWinters(ts(value, frequency = 12)); list(alpha = m$alpha[[1]], beta = m$beta[[1]], gamma = m$gamma[[1]], sse = m$SSE)", co2)
    refs.run_all("co2", "forecast (repeats_every=12, show='model')", {"alpha": r["alpha"], "beta": r["beta"], "gamma": r["gamma"],
                                                                     "sum of squared errors": r["sse"]}, "R HoltWinters (help page co2)",
                 rman("stats", "HoltWinters"), "stats", ITERATIVE)
    r = r_run('''y <- ts(value, frequency = 12); n <- length(y); tr <- window(y, end = time(y)[n - 12]); te <- window(y, start = time(y)[n - 11])
                 f <- forecast::forecast(HoltWinters(tr), h = 12); list(forecast = as.numeric(f$mean), error = as.numeric(te - f$mean))''', co2, packages=("forecast",))
    refs.run_vector("co2", "forecast (repeats_every=12, check_on_last=12)", "forecast of the held-out year", r["forecast"], "R HoltWinters on the first rows",
                    rman("stats", "HoltWinters"), "stats", ITERATIVE)
    refs.run_vector("co2", "forecast (repeats_every=12, check_on_last=12)", "error (actual minus forecast)", r["error"], "R, actual minus forecast",
                    rman("forecast", "accuracy"), "forecast", ITERATIVE)
    METHODS = {"exponential smoothing": "forecast::forecast(HoltWinters(tr), h = 12)", "same season last time": "forecast::snaive(tr, h = 12)",
               "straight on": "forecast::rwf(tr, drift = TRUE, h = 12)", "last value": "forecast::naive(tr, h = 12)",
               "arima": "forecast::forecast(arima(tr, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1)), method = 'ML'), h = 12)"}
    for name, data, label, methods in (("co2", co2, "compare_forecasts_for ()", ["exponential smoothing", "same season last time", "straight on", "last value"]),
                                       ("USAccDeaths", usacc, "compare_forecasts_for (test_method=['exponential smoothing', 'same season last time', 'arima'], "
                                        "order=(0, 1, 1), seasonal_order=(0, 1, 1))", ["exponential smoothing", "same season last time", "arima"])):
        for m in methods:
            r = r_run(f'''y <- ts(value, frequency = 12); n <- length(y); tr <- window(y, end = time(y)[n - 12]); te <- window(y, start = time(y)[n - 11])
                          a <- forecast::accuracy({METHODS[m]}, te)
                          list(MAE = a[2, "MAE"], RMSE = a[2, "RMSE"], MAPE = a[2, "MAPE"], MASE = a[2, "MASE"])''', data, packages=("forecast",))
            refs.run_all(name, label, {f"{m}: {k}": v for k, v in r.items()}, "R forecast accuracy (each method made directly, test set as ts)",
                         rman("forecast", "accuracy"), "forecast", ITERATIVE)

# 3. smoothing on the Nile flow
x = nile["value"].to_numpy()
if HAS_R:
    r = r_run('''list(ma = as.numeric(stats::filter(value, rep(1 / 5, 5), sides = 2)), w = as.numeric(stats::filter(value, c(1, 2, 1) / 4, sides = 2)),
                      med = as.numeric(runmed(value, 5, endrule = "median")), lw = lowess(seq_along(value), value)$y,
                      sp = predict(smooth.spline(seq_along(value) - 1, value, lambda = 1000, all.knots = TRUE))$y,
                      gcv = predict(smooth.spline(seq_along(value) - 1, value, all.knots = TRUE))$y,
                      sg = as.numeric(signal::sgolayfilt(value, p = 2, n = 11)))''', nile, packages=("signal",))
    ma = [v if v == v else 0.0 for v in r["ma"]]
    refs.run_vector("Nile", 'add_smoothed_column (test_method="moving average", window=5)', "smoothed", r["ma"][2:-2], "R stats filter (sides = 2)",
                    rman("stats", "filter"), "stats", CLOSED)
    refs.run_vector("Nile", 'add_smoothed_column (test_method="moving average", weights=[1, 2, 1])', "smoothed", r["w"][1:-1], "R stats filter (sides = 2)",
                    rman("stats", "filter"), "stats", CLOSED)
    refs.run_vector("Nile", 'nums.smooth (test_method="median", window=5)', "smoothed", r["med"], "R stats runmed (endrule = median)",
                    rman("stats", "runmed"), "stats", CLOSED)
    refs.run_vector("Nile", 'add_smoothed_column (test_method="lowess")', "smoothed", r["lw"], "R stats lowess (f = 2/3)", rman("stats", "lowess"), "stats", ITERATIVE)
    refs.run_vector("Nile", 'add_smoothed_column (test_method="spline", smoothness=1000)', "smoothed", r["sp"], "R smooth.spline (lambda = 1000, all knots)",
                    rman("stats", "smooth.spline"), "stats", ITERATIVE)
    refs.run_vector("Nile", 'nums.smooth (test_method="spline", smoothness="auto")', "smoothed", r["gcv"], "R smooth.spline (GCV, all knots)",
                    rman("stats", "smooth.spline"), "stats", ITERATIVE)
    refs.run_vector("Nile", 'add_smoothed_column (test_method="savitzky-golay", window=11)', "smoothed", r["sg"], "R signal sgolayfilt (degree 2)",
                    rman("signal", "sgolayfilt"), "signal", CLOSED)
    r = r_run('''t <- seq_along(value); p <- predict(loess(value ~ t), se = TRUE); q <- qt(0.975, p$df)
                 list(fit = p$fit, low = p$fit - q * p$se.fit, high = p$fit + q * p$se.fit)''', nile)
    for key, name in (("fit", "smoothed"), ("low", "95 % low"), ("high", "95 % high")):
        refs.run_vector("Nile", 'add_smoothed_column (test_method="loess", interval=95)', name, r[key], "R stats loess + predict(se = TRUE)",
                        rman("stats", "loess"), "stats", ITERATIVE)
refs.run_vector("Nile", 'nums.smooth (test_method="savitzky-golay", window=11)', "smoothed", signal.savgol_filter(x, 11, 2), "SciPy signal.savgol_filter",
                "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.signal.savgol_filter.html", "scipy", CLOSED)
refs.run_vector("Nile", 'add_smoothed_column (test_method="savitzky-golay", window=11, of_degree=2)', "smoothed", signal.savgol_filter(x, 11, 2),
                "SciPy signal.savgol_filter (polyorder 2)", "https://docs.scipy.org/doc/scipy-1.18.1/reference/generated/scipy.signal.savgol_filter.html",
                "scipy", CLOSED)

# 4. trend and season; removing a trend
if HAS_R:
    for name, data, label, call in (("nottem", nottem, "season_window=7", "stl(y, s.window = 7)"),
                                    ("co2", co2, "season_window=7, robust=True", "stl(y, s.window = 7, robust = TRUE)"),
                                    ("nottem", nottem, "season_window='fixed'", 'stl(y, s.window = "periodic")')):
        r = r_run(f"y <- ts(value, frequency = 12); s <- {call}$time.series; list(season = s[, 1], trend = s[, 2])", data)
        for key in ("trend", "season"):
            refs.run_vector(name, f"add_trend_and_season ({label})", key, r[key], f"R stats {call} (help page series)", rman("stats", "stl"), "stats", ITERATIVE)
    for name, data, season, kind in (("co2", co2, "added", "additive"), ("AirPassengers", air, "multiplied", "multiplicative")):
        r = r_run(f'd <- decompose(ts(value, frequency = 12), type = "{kind}"); list(trend = as.numeric(d$trend), season = as.numeric(d$seasonal))', data)
        for key in ("trend", "season"):
            values = r[key][6:-6] if key == "trend" else r[key]
            refs.run_vector(name, f'add_trend_and_season (test_method="moving averages", season="{season}")', key, values,
                            "R stats decompose (help page co2)", rman("stats", "decompose"), "stats", CLOSED)
    r = r_run('''t <- seq_along(value); m <- factor(rep(1:12, length.out = length(value))); f <- lm(value ~ t + m, contrasts = list(m = contr.sum))
                 e <- c(coef(f)[3:13], -sum(coef(f)[3:13])); list(trend = coef(f)[1] + coef(f)[2] * t, season = e[as.integer(m)])''', co2)
    for key in ("trend", "season"):
        refs.run_vector("co2", 'add_trend_and_season (test_method="straight line")', key, r[key], "R stats lm(y ~ t + month), sum-to-zero",
                        rman("stats", "lm"), "stats", CLOSED)
    r = r_run('''a <- value[1:28]; b <- value[29:100]
                 list(r = residuals(lm(value ~ seq_along(value))), b = c(residuals(lm(a ~ seq_along(a))), residuals(lm(b ~ seq_along(b)))))''', nile)
    refs.run_vector("Nile", "nums.remove_trend ()", "detrended", r["r"], "R residuals(lm(x ~ seq_along(x)))", rman("stats", "lm"), "stats", CLOSED)
    refs.run_vector("Nile", "nums.remove_trend (breaks=[28])", "detrended", r["b"], "R residuals(lm()) of each piece", rman("stats", "lm"), "stats", CLOSED)

# 5. autocorrelation; Moran's I on the globe
if HAS_R:
    for name, data in (("lh", lh), ("LakeHuron", lakehuron)):
        for show, call in (("correlation", "acf(value, lag.max = 10, plot = FALSE)$acf[-1]"), ("partial", "pacf(value, lag.max = 10, plot = FALSE)$acf"),
                           ("covariance", 'acf(value, lag.max = 10, type = "covariance", plot = FALSE)$acf[-1]')):
            r = r_run(f"list(v = as.numeric({call}), lb = Box.test(value, lag = 10, type = 'Ljung-Box')$statistic[[1]], n = length(value))", data)
            for word in ("find_autocorrelation_of", "nums.find_autocorrelation"):
                refs.run_vector(name, f'{word} (show="{show}", rows_back=10)', show, r["v"], f"R stats {call.split('(')[0]}", rman("stats", "acf"), "stats", CLOSED)
            if show == "correlation":
                refs.run(name, 'find_autocorrelation_of (show="correlation", rows_back=10)', "Ljung-Box statistic", r["lb"], "R stats Box.test (Ljung-Box)",
                         rman("stats", "Box.test"), "stats", CLOSED)
                acf = np.asarray(r["v"])
                band = 1.959963984540054 * np.sqrt((1 + 2 * np.concatenate([[0], np.cumsum(acf ** 2)[:-1]])) / r["n"])
                refs.run_vector(name, 'find_autocorrelation_of (band="bartlett", rows_back=10)', "band high", band,
                                "R plot.acf (ci.type = 'ma'): qnorm(0.975) sqrt((1 + 2 sum r^2) / n)", rman("stats", "plot.acf"), "stats", CLOSED)
    for seed, permutations in ((1, 999), (9, 499)):
        r = r_run(f'''xy <- cbind(longitude, lat); nb <- spdep::knn2nb(spdep::knearneigh(xy, k = 8, longlat = TRUE)); w <- spdep::nb2listw(nb, style = "W")
                      set.seed({seed}); m <- spdep::moran.mc(mag, w, nsim = {permutations}); list(I = m$statistic[[1]], p = m$p.value)''', quakes, packages=("spdep",))
        check = f"find_autocorrelation_of (on_globe=True, nearest=8, permutations={permutations}, with_seed={seed})"
        refs.run("quakes", check, "Moran's I", r["I"], "R spdep moran.mc (knearneigh longlat = TRUE)", rman("spdep", "moran.mc"), "spdep", CLOSED)
        refs.run("quakes", check, "p-value", r["p"], "R spdep moran.mc (same seed), 4 MCSE", rman("spdep", "moran.mc"), "spdep",
                 4 * np.sqrt(max(r["p"] * (1 - r["p"]), 1 / permutations) / permutations), "absolute")

# 6. steady or wandering (R urca and tseries)
if HAS_R:
    for name, data in (("Nile", nile), ("LakeHuron", lakehuron)):
        for lag in [int(np.trunc((len(data) - 1) ** (1 / 3)))]:   # tseries' default number of lags, trunc((n - 1)^(1/3))
            r = r_run(f"u <- urca::ur.df(value, type = 'drift', lags = {lag}); list(s = u@teststat[1])", data, packages=("urca",))
            refs.run(name, f"check_steady_over_time: ADF statistic with {lag} rows back", "ADF statistic", r["s"], f"R urca ur.df (drift, lag {lag})",
                     rman("urca", "ur.df"), "urca", CLOSED)
        r = r_run("k <- tseries::kpss.test(value, null = 'Level'); list(s = k$statistic[[1]], p = k$p.value)", data, packages=("tseries",))
        for test_method in ("both", "kpss"):
            refs.run_all(name, f'check_steady_over_time (test_method="{test_method}")', {"KPSS statistic": r["s"], "KPSS p-value": r["p"]},
                         "R tseries kpss.test (Level)", rman("tseries", "kpss.test"), "tseries", CLOSED)

# trend: the trend vignette (printed)
TV = "Pohlert, Non-Parametric Trend Tests and Change-Point Detection (the trend package vignette)"
TV_URL = VIGNETTE.format(package="trend", file="trend.pdf")
tv = ("trend", r_version("trend"))
refs.printed_all("maxau Q (trend)", "find_trend_of", {"z": "-1.3989", "p-value": "0.1619", "S": "-144", "Kendall's tau": "-0.1454545"},
                 TV + ", mk.test(maxau[, 'Q'])", TV_URL, "printed", *tv)
refs.printed_all("maxau s (trend)", "find_trend_of", {"Sen's slope": "-0.2876139", "95 % low": "-0.4196477", "95 % high": "-0.1519026"},
                 TV + ", sens.slope(maxau[, 's'])", TV_URL, "printed", *tv)
refs.printed_all("nottem (trend)", "find_trend_of (repeats_every=12)", {"z": "2.0919", "p-value": "0.03645", "seasonal Sen's slope (per year)": "0.05"},
                 TV + ", smk.test(nottem) and sea.sens.slope(nottem)", TV_URL, "printed", *tv)

# change points: Cobb (1978), printed; the PELT searches run once in R changepoint
refs._add("Nile", "find_change_points_in", "year in which the first piece ends", 1898.0, "1898", "published",
          "Cobb (1978), The problem of the Nile: conditional solution to a changepoint problem, Biometrika 65(2), 243-251",
          "https://doi.org/10.1093/biomet/65.2.243", "", "", 0, "absolute")
if HAS_R:
    for change, penalty, fun in (("average", "mbic", "cpt.mean"), ("average", "bic", "cpt.mean"), ("average", "aic", "cpt.mean"),
                                 ("average", "hannan-quinn", "cpt.mean"), ("average and spread", "mbic", "cpt.meanvar")):
        pen = {"mbic": "MBIC", "bic": "BIC", "aic": "AIC", "hannan-quinn": "Hannan-Quinn"}[penalty]
        scale = "x <- value / (mad(diff(value)) / sqrt(2))" if fun == "cpt.mean" else "x <- value"
        r = r_run(f'''{scale}; c <- changepoint::{fun}(x, method = "PELT", penalty = "{pen}"); list(ends = c(changepoint::cpts(c), length(x)))''', nile,
                  packages=("changepoint",))
        ends = np.atleast_1d(r["ends"])
        values = {"number of pieces": len(ends)}
        values.update({f"piece {i + 1} ends at reading": e for i, e in enumerate(ends)})
        refs.run_all("Nile", f'find_change_points_in (change="{change}", penalty="{penalty}")', values, f"R changepoint {fun} (PELT)",
                     rman("changepoint", fun), "changepoint", 0, "absolute")

# 7. the numbers the time charts draw (chart.get_numbers_drawn()): R acf, ccf (lag sign reversed), the readings one row
# back, stl and HoltWinters + forecast on the help-page series
deaths = pd.DataFrame({"step": np.arange(72), "mdeaths": r_series("mdeaths")["value"].to_numpy(),
                       "fdeaths": r_series("fdeaths")["value"].to_numpy()})
refs.data("deaths", deaths)
lh_values = lh["value"].to_numpy(dtype=float)
refs.run_vector("lh", "chart lag_dots (rows_back=1)", "dot across (the reading one row back)", lh_values[:-1],
                "R lh (stats::lag.plot puts each reading against the one before)", rman("stats", "lag.plot"), "stats", CLOSED)
refs.run_vector("lh", "chart lag_dots (rows_back=1)", "dot up (the reading)", lh_values[1:],
                "R lh (stats::lag.plot puts each reading against the one before)", rman("stats", "lag.plot"), "stats", CLOSED)
if HAS_R:
    r = r_run("list(v = as.numeric(acf(value, lag.max = 10, plot = FALSE)$acf[-1]))", lh)
    refs.run_vector("lh", "chart autocorrelation_bars (rows_back=10)", "bar height", r["v"], "R stats acf (plot = FALSE)",
                    rman("stats", "acf"), "stats", CLOSED)
    refs.run_vector("lh", "chart autocorrelation_bars (rows_back=10)", "bar position (rows back)", np.arange(1, 11), "R stats acf, lags 1 to 10",
                    rman("stats", "acf"), "stats", CLOSED)
    r = r_run("c <- ccf(mdeaths, fdeaths, lag.max = 12, plot = FALSE); list(v = rev(as.numeric(c$acf)), s = -rev(as.numeric(c$lag)))", deaths)
    refs.run_vector("mdeaths and fdeaths", "chart cross_correlation_bars (rows_back=12)", "stem height", r["v"],
                    "R stats ccf(mdeaths, fdeaths), lag sign reversed", rman("stats", "acf"), "stats", CLOSED)
    refs.run_vector("mdeaths and fdeaths", "chart cross_correlation_bars (rows_back=12)", "stem position (shift)", r["s"],
                    "R stats ccf(mdeaths, fdeaths), lag sign reversed", rman("stats", "acf"), "stats", CLOSED)
    r = r_run("list(trend = as.numeric(stl(ts(value, frequency = 12), s.window = 7)$time.series[, 2]))", co2)
    refs.run_vector("co2", "chart trend_and_season (repeats_every=12)", "trend panel line", r["trend"], "R stats stl(co2, s.window = 7)",
                    rman("stats", "stl"), "stats", ITERATIVE)
    r = r_run("""f <- forecast::forecast(HoltWinters(ts(value, frequency = 12)), h = 12, level = 95)
                 list(f = as.numeric(f$mean), low = as.numeric(f$lower), high = as.numeric(f$upper))""", co2, packages=("forecast",))
    for key, name in (("f", "forecast line"), ("low", "band low"), ("high", "band high")):
        refs.run_vector("co2", "chart forecast (repeats_every=12)", name, r[key], "R HoltWinters + forecast (h = 12, level = 95)",
                        rman("stats", "HoltWinters"), "stats", ITERATIVE)

refs.save()
