"""An independent converged optimum for drc's four-parameter dose-response fits (LL.4 and W1.4), not drc itself.

drc's relTol is a relative change of the residual sum of squares (RSS); on a flat optimum relTol = 1e-12 fixes the
parameters only to about its square root (1e-6). To tell whether drc's reference or simplydata stopped short, the same
least-squares problem is solved to convergence by two other optimisers, both started from drc's estimates:

* R `optim` (BFGS, reltol 1e-15) on the RSS with its analytic gradient. Its Hessian (optim(hessian = TRUE),
  differentiating the analytic gradient) gives the standard errors as drc defines them for a continuous response:
  vcov = sigma^2 * solve(Hessian / 2), sigma^2 = RSS / (n - 4);
* SciPy `least_squares` (xtol = ftol = gtol = 1e-15) on the residuals with their analytic Jacobian.

The effective doses are drc's formulas written out: LL.4 ED(p) = e * (p / (100 - p))^(1 / b), W1.4 ED(p) =
e * (-log(1 - p / 100))^(1 / b); their standard errors by the delta method and their 95 % intervals as drc's
ED(interval = "delta") gives them for a continuous response, ED -/+ qt(0.975, n - 4) * SE.
`rss` gives the RSS at any set of estimates, evaluated in R, as evidence of which estimates are nearer the optimum.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize, stats

from checks import r_run

MODEL_R = '''
f_of <- function(p, x, kind) { u <- exp(p[1] * (log(x) - log(p[4])))
  if (kind == "LL.4") p[2] + (p[3] - p[2]) / (1 + u) else p[2] + (p[3] - p[2]) * exp(-u) }
'''


def _model(kind, p, x):
    b, c, d, e = p
    with np.errstate(divide="ignore"):
        u = np.exp(b * (np.log(x) - np.log(e)))
    return c + (d - c) / (1 + u) if kind == "LL.4" else c + (d - c) * np.exp(-u)


def _jacobian(kind, p, x):
    """d f / d (b, c, d, e); at x = 0 the curve is d whatever b and e are."""
    b, c, d, e = p
    with np.errstate(divide="ignore", invalid="ignore"):
        lx = np.log(x) - np.log(e)
        u = np.exp(b * lx)
        du_db = np.where(x > 0, u * lx, 0.0)
        du_de = np.where(x > 0, -u * b / e, 0.0)
    if kind == "LL.4":
        g = 1 / (1 + u)
        df_du = -(d - c) * g ** 2
        cols = [df_du * du_db, 1 - g, g, df_du * du_de]
    else:
        h = np.exp(-u)
        df_du = -(d - c) * h
        cols = [df_du * du_db, 1 - h, h, df_du * du_de]
    return np.column_stack(cols)


def _ed(kind, p, levels):
    b, e = p[0], p[3]
    k = np.array([lv / (100 - lv) if kind == "LL.4" else -np.log(1 - lv / 100) for lv in levels])
    ed = e * k ** (1 / b)
    grad = np.column_stack([-ed * np.log(k) / b ** 2, np.zeros_like(ed), np.zeros_like(ed), ed / e])
    return ed, grad


def optim_r(x, y, kind, start, levels=(10, 50, 90)):
    """R optim (BFGS, reltol 1e-15, analytic gradient) from `start`: estimates, standard errors (drc's definition),
    effective doses with their delta standard errors and 95 % intervals, and the RSS."""
    code = MODEL_R + f'''
kind <- "{kind}"
rss <- function(p) sum((y - f_of(p, x, kind))^2)
gr <- function(p) {{ b <- p[1]; c <- p[2]; d <- p[3]; e <- p[4]
  lx <- ifelse(x > 0, log(x) - log(e), 0); u <- ifelse(x > 0, exp(b * lx), 0)
  if (kind == "LL.4") {{ g <- 1 / (1 + u); dfu <- -(d - c) * g^2; dc <- 1 - g; dd <- g }}
  else {{ h <- exp(-u); dfu <- -(d - c) * h; dc <- 1 - h; dd <- h }}
  J <- cbind(dfu * u * lx, dc, dd, ifelse(x > 0, -dfu * u * b / e, 0))
  -2 * colSums((y - f_of(p, x, kind)) * J) }}
o <- optim(start, rss, gr, method = "BFGS", control = list(reltol = 1e-15, maxit = 100000), hessian = TRUE)
n <- length(y); s2 <- o$value / (n - 4); V <- s2 * solve(o$hessian / 2)
list(est = o$par, se = sqrt(diag(V)), V = as.vector(V), rss = o$value, conv = o$convergence, df = n - 4)'''
    r = r_run(code, {"x": np.asarray(x, dtype=float), "y": np.asarray(y, dtype=float),
                     "start": np.asarray(start, dtype=float)})
    if int(r["conv"]) != 0:
        raise SystemExit(f"optim did not converge ({kind}): code {r['conv']}")
    p = np.asarray(r["est"], dtype=float)
    V = np.asarray(r["V"], dtype=float).reshape(4, 4, order="F")
    ed, grad = _ed(kind, p, levels)
    ed_se = np.sqrt(np.einsum("ij,jk,ik->i", grad, V, grad))
    t = stats.t.ppf(0.975, int(r["df"]))
    return {"est": p, "se": np.asarray(r["se"], dtype=float), "ed": ed, "ed_se": ed_se,
            "ed_low": ed - t * ed_se, "ed_high": ed + t * ed_se, "rss": float(r["rss"])}


def least_squares_scipy(x, y, kind, start, levels=(10, 50, 90)):
    """SciPy least_squares (xtol = ftol = gtol = 1e-15, analytic Jacobian) from `start`: estimates and doses."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    fit = optimize.least_squares(lambda p: _model(kind, p, x) - y, np.asarray(start, dtype=float),
                                 jac=lambda p: _jacobian(kind, p, x), method="trf", xtol=1e-15, ftol=1e-15,
                                 gtol=1e-15, max_nfev=100000)
    ed, _ = _ed(kind, fit.x, levels)
    return {"est": fit.x, "ed": ed}


def rss(x, y, kind, estimates: dict) -> dict:
    """The RSS at each named set of estimates, evaluated in R."""
    names = list(estimates)
    code = MODEL_R + f'''kind <- "{kind}"; P <- matrix(p, nrow = 4)
list(rss = apply(P, 2, function(q) sum((y - f_of(q, x, kind))^2)))'''
    r = r_run(code, {"x": np.asarray(x, dtype=float), "y": np.asarray(y, dtype=float),
                     "p": np.concatenate([np.asarray(estimates[n], dtype=float) for n in names])})
    return dict(zip(names, np.atleast_1d(np.asarray(r["rss"], dtype=float))))
