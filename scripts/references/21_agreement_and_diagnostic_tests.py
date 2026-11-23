"""Reference values for notebooks/21_agreement_and_diagnostic_tests.ipynb, saved to references/21_agreement_and_diagnostic_tests.csv.

Computed with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through rpy2) and from the
printed, certified or help-page values cited in the notebook; the data the notebook reads are written to
references/data/21_agreement_and_diagnostic_tests/. The simplydata words also appear here because this script was the notebook's working
source: their values are not stored (Recorder ignores them); the notebook computes them itself. The text of each
section is in the notebook.

    python scripts/references/21_agreement_and_diagnostic_tests.py      (from the repository root; needs R with the packages listed)
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
from checks import Checks, load_csv, quiet, r_available, r_run, read_row, read_rows, versions

versions("pandas", "numpy", "scipy", "statsmodels", "rpy2")
checks = Recorder("21_agreement_and_diagnostic_tests")
HAS_R = r_available()
CLOSED = 1e-10
if HAS_R:
    print(rd.r_version_of("irr", "DescTools", "BlandAltmanLeh", "psych", "pROC", "epiR"))


def table_of(word, *args, **options):
    """Every row of a simplydata result table, each as {column: value} (read with table.get_row)."""
    return read_rows(quiet(word, *args, **options))


def keyed(rows, column="measure"):
    """The rows read by table_of, keyed by one column (the measure)."""
    return {row[column]: row for row in rows}


def columns_of(found):
    """The column names of rows read by table_of (a list of rows, or rows keyed by keyed())."""
    first = found[0] if isinstance(found, list) else next(iter(found.values()))
    return list(first)


def one(word, *args, **options):
    return read_row(quiet(word, *args, **options))


def compare_many(dataset, method, found, reference, reference_name, tolerance=CLOSED, how="relative", prefix=""):
    for quantity, value in reference.items():
        checks.compare(dataset, method, f"{prefix}: {quantity}" if prefix else quantity, found[quantity], value,
                       reference_name, tolerance, how)


diagnoses = rd.r_dataset("diagnoses", "irr")
raters = list(diagnoses.columns)
psych = sdata.table(diagnoses)
psych.agreement_between("rater1", "rater2")

ds = "diagnoses (R irr example; Fleiss 1971)"


def cohen_kappa(a, b, categories, weights=None):
    """Cohen's kappa (weights: None, "linear" or "quadratic"), its standard error and the null standard error."""
    k = len(categories)
    index = {c: i for i, c in enumerate(categories)}
    p = np.zeros((k, k))
    for x, y in zip(a, b):
        p[index[x], index[y]] += 1
    n = p.sum()
    p /= n
    i, j = np.indices((k, k))
    if weights is None:
        w = (i == j).astype(float)
    elif weights == "linear":
        w = 1 - np.abs(i - j) / (k - 1)
    else:
        w = 1 - (i - j) ** 2 / (k - 1) ** 2
    rows, cols = p.sum(axis=1), p.sum(axis=0)
    po, pe = (w * p).sum(), (w * np.outer(rows, cols)).sum()
    kappa = (po - pe) / (1 - pe)
    wbar_i, wbar_j = w @ cols, rows @ w
    var = ((p * (w - (wbar_i[:, None] + wbar_j[None, :]) * (1 - kappa)) ** 2).sum()
           - (kappa - pe * (1 - kappa)) ** 2) / (n * (1 - pe) ** 2)
    var0 = ((np.outer(rows, cols) * (w - (wbar_i[:, None] + wbar_j[None, :])) ** 2).sum() - pe ** 2) / (n * (1 - pe) ** 2)
    return kappa, np.sqrt(var), np.sqrt(var0)


categories = sorted(set(diagnoses.values.ravel()))
kappa, se, se0 = cohen_kappa(diagnoses["rater1"], diagnoses["rater2"], categories)
found = one(psych.agreement_between, "rater1", "rater2")
compare_many(ds, "agreement_between (2 columns of categories: Cohen's kappa)", found,
             {"estimate": kappa, "estimate 95 % CI low": kappa - stats.norm.ppf(0.975) * se,
              "estimate 95 % CI high": kappa + stats.norm.ppf(0.975) * se, "statistic": kappa / se0,
              "p-value": 2 * stats.norm.sf(kappa / se0)}, "Cohen (1960); Fleiss, Cohen & Everitt (1969), NumPy")
if HAS_R:
    r = r_run("""k <- irr::kappa2(cbind(a, b)); list(estimate = k$value, statistic = k$statistic, `p-value` = k$p.value)""",
              {"a": diagnoses["rater1"], "b": diagnoses["rater2"]})
    compare_many(ds, "agreement_between (2 columns of categories: Cohen's kappa)", found, r, "R irr kappa2 (help page example)")

# Fleiss' kappa for the six raters
counts = np.array([[sum(row == c) for c in categories] for row in diagnoses.values])
n_r = counts.sum(axis=1)[0]
p_j = counts.sum(axis=0) / counts.sum()
P_i = ((counts ** 2).sum(axis=1) - n_r) / (n_r * (n_r - 1))
fleiss = (P_i.mean() - (p_j ** 2).sum()) / (1 - (p_j ** 2).sum())
found = one(psych.agreement_between, *raters)
checks.compare(ds, "agreement_between (3 or more columns of categories: Fleiss' kappa)", "estimate", found["estimate"],
               fleiss, "Fleiss (1971), NumPy", CLOSED, "relative")
if HAS_R:
    r = r_run("""k <- irr::kappam.fleiss(data.frame(r1, r2, r3, r4, r5, r6)); list(estimate = k$value, statistic = unname(k$statistic))""",
              {f"r{i + 1}": diagnoses[c] for i, c in enumerate(raters)})
    compare_many(ds, "agreement_between (3 or more columns of categories: Fleiss' kappa)", found, r,
                 "R irr kappam.fleiss (help page example)")


anxiety = rd.r_dataset("anxiety", "irr")
levels = [f"rating {v}" for v in range(1, 7)]
anx_text = "rating " + anxiety.astype(int).astype(str)
ds = "anxiety (R irr example)"
for near, irr_weight in [("linear", "equal"), ("quadratic", "squared")]:
    found = one(sdata.table(anx_text).agreement_between, "rater1", "rater2", near_misses=near, in_order=levels)
    kappa, se, se0 = cohen_kappa(anx_text["rater1"], anx_text["rater2"], levels, near)
    checks.compare(ds, f'agreement_between (near_misses="{near}", in_order=)', "estimate", found["estimate"], kappa,
                   "weighted kappa written out, NumPy", CLOSED, "relative")
    if HAS_R:
        r = r_run(f"""a <- factor(a, levels = 1:6); b <- factor(b, levels = 1:6)
            k <- irr::kappa2(data.frame(a, b), weight = "{irr_weight}"); list(estimate = k$value)""",
                  {"a": anxiety["rater1"], "b": anxiety["rater2"]})
        checks.compare(ds, f'agreement_between (near_misses="{near}", in_order=)', "estimate", found["estimate"],
                       r["estimate"], f'R irr kappa2(weight = "{irr_weight}")', CLOSED, "relative")


X = anxiety.to_numpy(dtype=float)
n, k = X.shape
grand = X.mean()
ssr = k * ((X.mean(axis=1) - grand) ** 2).sum()
ssc = n * ((X.mean(axis=0) - grand) ** 2).sum()
sst = ((X - grand) ** 2).sum()
sse = sst - ssr - ssc
msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
msw = (ssc + sse) / (n * (k - 1))
z = stats.f.ppf


def agreement_interval(icc_value, single):
    """McGraw & Wong (1996), the two-way absolute-agreement interval with approximate degrees of freedom."""
    a = k * icc_value / (n * (1 - icc_value)) if single else icc_value / (n * (1 - icc_value))
    b = 1 + k * icc_value * (n - 1) / (n * (1 - icc_value)) if single else 1 + icc_value * (n - 1) / (n * (1 - icc_value))
    v = (a * msc + b * mse) ** 2 / ((a * msc) ** 2 / (k - 1) + (b * mse) ** 2 / ((n - 1) * (k - 1)))
    fl, fu = z(0.975, n - 1, v), z(0.975, v, n - 1)
    if single:
        low = n * (msr - fl * mse) / (fl * (k * msc + (k * n - k - n) * mse) + n * msr)
        high = n * (fu * msr - mse) / (k * msc + (k * n - k - n) * mse + n * fu * msr)
    else:
        low = n * (msr - fl * mse) / (fl * (msc - mse) + n * msr)
        high = n * (fu * msr - mse) / (msc - mse + n * fu * msr)
    return low, high


forms = {}
F = msr / mse
fl, fu = F / z(0.975, n - 1, (n - 1) * (k - 1)), F * z(0.975, (n - 1) * (k - 1), n - 1)
forms[("a sample", "agreement", "one rater")] = ((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n), F,
                                                 *agreement_interval((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n), True))
forms[("a sample", "agreement", "average of raters")] = ((msr - mse) / (msr + (msc - mse) / n), F,
                                                         *agreement_interval((msr - mse) / (msr + (msc - mse) / n), False))
forms[("a sample", "consistency", "one rater")] = ((msr - mse) / (msr + (k - 1) * mse), F, (fl - 1) / (fl + k - 1),
                                                   (fu - 1) / (fu + k - 1))
forms[("a sample", "consistency", "average of raters")] = ((msr - mse) / msr, F, 1 - 1 / fl, 1 - 1 / fu)
forms[("these only", "consistency", "one rater")] = forms[("a sample", "consistency", "one rater")]
Fw = msr / msw
flw, fuw = Fw / z(0.975, n - 1, n * (k - 1)), Fw * z(0.975, n * (k - 1), n - 1)
forms[("different per row", "agreement", "one rater")] = ((msr - msw) / (msr + (k - 1) * msw), Fw,
                                                          (flw - 1) / (flw + k - 1), (fuw - 1) / (fuw + k - 1))
forms[("different per row", "agreement", "average of raters")] = ((msr - msw) / msr, Fw, 1 - 1 / flw, 1 - 1 / fuw)
anx_t = sdata.table(anxiety)
for (who, checking, of), (estimate, F_value, low, high) in forms.items():
    found = one(anx_t.agreement_between, *anxiety.columns, raters=who, checking_for=checking, of=of)
    label = f'agreement_between (raters="{who}", checking_for="{checking}", of="{of}")'
    compare_many(ds, label, found, {"estimate": estimate, "statistic": F_value, "estimate 95 % CI low": low,
                                    "estimate 95 % CI high": high}, "McGraw & Wong (1996), NumPy")
    if HAS_R:
        model = "oneway" if who == "different per row" else "twoway"
        unit = "single" if of == "one rater" else "average"
        r = r_run(f"""i <- irr::icc(cbind(a, b, c), model = "{model}", type = "{checking}", unit = "{unit}")
            list(estimate = i$value, statistic = i$Fvalue, `p-value` = i$p.value, low = i$lbound, high = i$ubound)""",
                  {"a": anxiety["rater1"], "b": anxiety["rater2"], "c": anxiety["rater3"]})
        compare_many(ds, label, found, {"estimate": r["estimate"], "statistic": r["statistic"], "p-value": r["p-value"],
                                        "estimate 95 % CI low": r["low"], "estimate 95 % CI high": r["high"]},
                     "R irr icc (help page example)")
        psych_row = {("different per row", "agreement", "one rater"): 0, ("a sample", "agreement", "one rater"): 1,
                     ("a sample", "consistency", "one rater"): 2, ("different per row", "agreement", "average of raters"): 3,
                     ("a sample", "agreement", "average of raters"): 4, ("a sample", "consistency", "average of raters"): 5}
        if (who, checking, of) in psych_row:
            i_ = psych_row[(who, checking, of)]
            r = r_run(f"""p <- psych::ICC(cbind(a, b, c), lmer = FALSE)$results[{i_ + 1}, ]
                list(estimate = p$ICC, low = p$`lower bound`, high = p$`upper bound`)""",
                      {"a": anxiety["rater1"], "b": anxiety["rater2"], "c": anxiety["rater3"]})
            compare_many(ds, label, found, {"estimate": r["estimate"], "estimate 95 % CI low": r["low"],
                                            "estimate 95 % CI high": r["high"]}, "R psych ICC")
for interval in (90, 99):
    found = one(anx_t.agreement_between, *anxiety.columns, checking_for="consistency", interval=interval)
    a_ = (100 - interval) / 200
    flx, fux = F / z(1 - a_, n - 1, (n - 1) * (k - 1)), F * z(1 - a_, (n - 1) * (k - 1), n - 1)
    lo_c = [c for c in found if "CI low" in c][0]
    hi_c = [c for c in found if "CI high" in c][0]
    checks.compare(ds, f'agreement_between (checking_for="consistency", interval={interval})', lo_c, found[lo_c],
                   (flx - 1) / (flx + k - 1), "McGraw & Wong (1996), NumPy", CLOSED, "relative")
    checks.compare(ds, f'agreement_between (checking_for="consistency", interval={interval})', hi_c, found[hi_c],
                   (fux - 1) / (fux + k - 1), "McGraw & Wong (1996), NumPy", CLOSED, "relative")


def krippendorff(units, metric):
    values = sorted({v for u in units for v in u})
    index = {v: i for i, v in enumerate(values)}
    o = np.zeros((len(values), len(values)))
    for u in units:
        m = len(u)
        if m < 2:
            continue
        for a_i, a in enumerate(u):
            for b_i, b in enumerate(u):
                if a_i != b_i:
                    o[index[a], index[b]] += 1 / (m - 1)
    vals = np.array(values, dtype=float) if metric == "interval" else None
    i, j = np.indices(o.shape)
    delta = (i != j).astype(float) if metric == "nominal" else (vals[i] - vals[j]) ** 2
    n_c = o.sum(axis=1)
    total = n_c.sum()
    d_o = (o * delta).sum() / total
    d_e = (np.outer(n_c, n_c) * delta).sum() / (total * (total - 1))
    return 1 - d_o / d_e


# The reference implementation reproduces Krippendorff's own worked example (the "C" data of irr's kripp.alpha
# help page, 4 observers and 12 units with gaps; published alpha 0.743 nominal and 0.849 interval, Krippendorff 2011)
made = rd.r_example_objects("irr", "kripp.alpha")
nmm = np.array(r_run("list(v = as.vector(nmm))")["v"]).reshape(-1, 4).T
nmm_units = [[v for v in nmm[:, j] if not np.isnan(v)] for j in range(nmm.shape[1])]
print("reference on Krippendorff's C data: nominal", round(krippendorff(nmm_units, "nominal"), 4), "(published 0.743);",
      "interval", round(krippendorff(nmm_units, "interval"), 4), "(published 0.849)")
assert round(krippendorff(nmm_units, "nominal"), 3) == 0.743 and round(krippendorff(nmm_units, "interval"), 3) == 0.849

for name, frame, metric, columns in [("anxiety (R irr example)", anxiety, "interval", list(anxiety.columns)),
                                     ("diagnoses (R irr example; Fleiss 1971)", diagnoses, "nominal", raters)]:
    found = one(sdata.table(frame).agreement_between, *columns, test_method="krippendorff")
    units = [list(row) for row in frame[columns].itertuples(index=False)]
    checks.compare(name, 'agreement_between (test_method="krippendorff")', "estimate", found["estimate"],
                   krippendorff(units, metric), f"Krippendorff (2011) coincidence matrix, NumPy ({metric})", CLOSED,
                   "relative")
    if HAS_R:
        codes = {v: i + 1 for i, v in enumerate(sorted(set(frame[columns].values.ravel())))}
        cols = {f"r{i}": (frame[c] if metric == "interval" else frame[c].map(codes).astype(float))
                for i, c in enumerate(columns)}
        r = r_run(f"""m <- t(cbind({", ".join(cols)})); list(estimate = irr::kripp.alpha(m, method = "{metric}")$value)""", cols)
        checks.compare(name, 'agreement_between (test_method="krippendorff")', "estimate", found["estimate"],
                       r["estimate"], "R irr kripp.alpha", CLOSED, "relative")

creatinine = load_csv("creatinine.csv")
for name, frame, a, b in [("anxiety (R irr example)", anxiety, "rater1", "rater2"),
                          ("creatinine, serum and plasma (R mcr)", creatinine, "serum", "plasma")]:
    x, y = frame[a].to_numpy(dtype=float), frame[b].to_numpy(dtype=float)
    nn = len(x)
    sx2, sy2 = x.var(), y.var()
    sxy = ((x - x.mean()) * (y - y.mean())).mean()
    rc = 2 * sxy / (sx2 + sy2 + (x.mean() - y.mean()) ** 2)
    r_ = sxy / np.sqrt(sx2 * sy2)
    u = (x.mean() - y.mean()) / (sx2 * sy2) ** 0.25
    var_z = ((1 - r_ ** 2) * rc ** 2 / ((1 - rc ** 2) * r_ ** 2) + 2 * rc ** 3 * (1 - rc) * u ** 2 / (r_ * (1 - rc ** 2) ** 2)
             - rc ** 4 * u ** 4 / (2 * r_ ** 2 * (1 - rc ** 2) ** 2)) / (nn - 2)
    zc = np.arctanh(rc)
    found = one(sdata.table(frame).agreement_between, a, b, test_method="lin")
    compare_many(name, 'agreement_between (test_method="lin")', found,
                 {"estimate": rc, "estimate 95 % CI low": np.tanh(zc - stats.norm.ppf(0.975) * np.sqrt(var_z)),
                  "estimate 95 % CI high": np.tanh(zc + stats.norm.ppf(0.975) * np.sqrt(var_z))},
                 "Lin (1989) written out, NumPy")


for interval in (95, 90):
    found = one(sdata.table(creatinine).limits_of_agreement, "serum", "plasma", interval=interval)
    d = (creatinine["serum"] - creatinine["plasma"]).to_numpy()
    nn, mean_d, s = len(d), d.mean(), d.std(ddof=1)
    t_q = stats.t.ppf(1 - (100 - interval) / 200, nn - 1)
    lower, upper = mean_d - 1.96 * s, mean_d + 1.96 * s
    reference = {"average difference": mean_d, "average difference 95 % CI low": mean_d - t_q * s / np.sqrt(nn),
                 "average difference 95 % CI high": mean_d + t_q * s / np.sqrt(nn), "lower limit": lower,
                 "lower limit 95 % CI low": lower - t_q * s * np.sqrt(3 / nn),
                 "lower limit 95 % CI high": lower + t_q * s * np.sqrt(3 / nn), "upper limit": upper,
                 "upper limit 95 % CI low": upper - t_q * s * np.sqrt(3 / nn),
                 "upper limit 95 % CI high": upper + t_q * s * np.sqrt(3 / nn)}
    if interval != 95:
        reference = {k_.replace("95 %", f"{interval} %"): v for k_, v in reference.items()}
        reference = {k_: v for k_, v in reference.items() if k_ in found}
    compare_many("creatinine, serum and plasma (R mcr)", f"limits_of_agreement (interval={interval})", found, reference,
                 "Bland & Altman (1986, 1999) written out, NumPy")

hs = load_csv("holzinger_swineford_1939.csv")
for name, frame, items in [("Holzinger & Swineford (1939) x4, x5, x6", hs, ["x4", "x5", "x6"]),
                           ("anxiety (R irr example)", anxiety, list(anxiety.columns))]:
    for interval in (95, 90):
        found = table_of(sdata.table(frame).reliability_of, items, interval=interval)
        Xi = frame[items].to_numpy(dtype=float)
        nn, kk = Xi.shape

        def alpha_of(cols):
            return cols.shape[1] / (cols.shape[1] - 1) * (1 - cols.var(axis=0, ddof=1).sum() / cols.sum(axis=1).var(ddof=1))

        alpha = alpha_of(Xi)
        q = (100 - interval) / 200
        low = 1 - (1 - alpha) * stats.f.ppf(1 - q, nn - 1, (nn - 1) * (kk - 1))
        high = 1 - (1 - alpha) * stats.f.ppf(q, nn - 1, (nn - 1) * (kk - 1))
        main = found[0]
        lo_c = [c for c in columns_of(found) if "CI low" in c][0]
        hi_c = [c for c in columns_of(found) if "CI high" in c][0]
        compare_many(name, f"reliability_of (interval={interval})", main,
                     {"estimate": alpha, lo_c: low, hi_c: high}, "Cronbach (1951), Feldt (1965) interval, NumPy")
        if interval == 95:
            for i, item in enumerate(items):
                line = next(row for row in found if row["measure"] == f"alpha without {item}")
                checks.compare(name, "reliability_of", f"alpha without {item}", line["estimate"],
                               alpha_of(np.delete(Xi, i, axis=1)), "Cronbach (1951), NumPy", CLOSED, "relative")


asah = rd.r_dataset("aSAH", "pROC")
patients = sdata.table(asah)
patients.compare_to_truth("s100b", truth="outcome")

def placements(cases, controls):
    """DeLong's structural components: V10 for each case, V01 for each control."""
    psi = (cases[:, None] > controls[None, :]) + 0.5 * (cases[:, None] == controls[None, :])
    return psi.mean(axis=1), psi.mean(axis=0)


def roc_reference(score, truth_is_case):
    cases, controls = score[truth_is_case], score[~truth_is_case]
    v10, v01 = placements(cases, controls)
    auc = v10.mean()
    var = v10.var(ddof=1) / len(cases) + v01.var(ddof=1) / len(controls)
    values = np.unique(score)
    cuts = (values[:-1] + values[1:]) / 2
    youden = [((cases > c).mean() + (controls <= c).mean() - 1) for c in cuts]
    return auc, np.sqrt(var), cuts[int(np.argmax(youden))], v10, v01


case = (asah["outcome"] == "Poor").to_numpy()
ds = "aSAH (R pROC example; Turck et al. 2010)"
for interval in (95, 90):
    found = table_of(patients.compare_to_truth, "s100b", truth="outcome", positive="Poor", interval=interval)
    found = keyed(found)
    lo_c = [c for c in columns_of(found) if "CI low" in c][0]
    hi_c = [c for c in columns_of(found) if "CI high" in c][0]
    auc, se, best, _, _ = roc_reference(asah["s100b"].to_numpy(), case)
    zq = stats.norm.ppf(1 - (100 - interval) / 200)
    label = f"compare_to_truth (a score, interval={interval})"
    row = found["area under the ROC curve (AUC)"]
    compare_many(ds, label, row, {"estimate": auc, lo_c: auc - zq * se, hi_c: auc + zq * se},
                 "Mann-Whitney AUC and DeLong (1988) variance, NumPy", prefix="AUC")
    checks.compare(ds, label, "best cut-off (Youden's index)", found["best cut-off (Youden's index)"]["estimate"],
                   best, "Youden's index over midpoints, NumPy", CLOSED, "relative")
    tp = int((asah["s100b"][case] > best).sum()); fn = int(case.sum() - tp)
    tn = int((asah["s100b"][~case] <= best).sum()); fp = int((~case).sum() - tn)
    if HAS_R:
        for measure, k_, n_ in [("sensitivity", tp, tp + fn), ("specificity", tn, tn + fp),
                                ("positive predictive value", tp, tp + fp), ("negative predictive value", tn, tn + fn)]:
            r = r_run(f"b <- binom.test({k_}, {n_}, conf.level = {interval / 100}); list(e = {k_} / {n_}, low = b$conf.int[1], high = b$conf.int[2])")
            compare_many(ds, label, found[measure], {"estimate": r["e"], lo_c: r["low"], hi_c: r["high"]},
                         "R binom.test", prefix=measure)
        r = r_run(f"""suppressPackageStartupMessages(library(pROC))
            r <- roc(outcome, s100b, levels = c("Good", "Poor"), direction = "<", quiet = TRUE)
            ci <- ci.auc(r, conf.level = {interval / 100}, method = "delong")
            b <- coords(r, "best", best.method = "youden", ret = "threshold")
            list(auc = ci[2], low = ci[1], high = ci[3], best = b$threshold[1])""", asah[["outcome", "s100b"]])
        compare_many(ds, label, row, {"estimate": r["auc"], lo_c: r["low"], hi_c: r["high"]},
                     "R pROC ci.auc (help page example)", prefix="AUC")

found = keyed(table_of(patients.compare_to_truth, "s100b", truth="outcome", positive="Poor", cut_off=0.3,
                 prevalence=10))
print(found)
# a value at the cut-off counts as positive (as pROC's coords: positive when the score is at or above the threshold)
tp = int((asah["s100b"][case] >= 0.3).sum()); fn = int(case.sum() - tp)
tn = int((asah["s100b"][~case] < 0.3).sum()); fp = int((~case).sum() - tn)
sens, spec = tp / (tp + fn), tn / (tn + fp)
label = "compare_to_truth (a score, cut_off=0.3, prevalence=10)"
checks.compare(ds, label, "sensitivity", found["sensitivity"]["estimate"], sens, "counted at the cut-off", CLOSED,
               "relative")
checks.compare(ds, label, "specificity", found["specificity"]["estimate"], spec, "counted at the cut-off", CLOSED,
               "relative")
ppv_rows = [m for m in found if "positive predictive value where" in m]
npv_rows = [m for m in found if "negative predictive value where" in m]
checks.compare(ds, label, ppv_rows[0], found[ppv_rows[0]]["estimate"], sens * 0.1 / (sens * 0.1 + (1 - spec) * 0.9),
               "Bayes' rule, NumPy", CLOSED, "relative")
checks.compare(ds, label, npv_rows[0], found[npv_rows[0]]["estimate"], spec * 0.9 / (spec * 0.9 + (1 - sens) * 0.1),
               "Bayes' rule, NumPy", CLOSED, "relative")


found = keyed(table_of(patients.compare_accuracy, "s100b", "ndka", truth="outcome", positive="Poor"))
print(found)
a1, _, _, v10a, v01a = roc_reference(asah["s100b"].to_numpy(), case)
a2, _, _, v10b, v01b = roc_reference(asah["ndka"].to_numpy(), case)
m_, n_ = case.sum(), (~case).sum()
cov10, cov01 = np.cov(v10a, v10b), np.cov(v01a, v01b)
var_diff = (cov10[0, 0] + cov10[1, 1] - 2 * cov10[0, 1]) / m_ + (cov01[0, 0] + cov01[1, 1] - 2 * cov01[0, 1]) / n_
zz = (a1 - a2) / np.sqrt(var_diff)
row = found["difference in AUC (s100b minus ndka)"]
compare_many(ds, "compare_accuracy", row, {"estimate": a1 - a2, "statistic": zz, "p-value": 2 * stats.norm.sf(abs(zz)),
                                           "estimate 95 % CI low": a1 - a2 - stats.norm.ppf(0.975) * np.sqrt(var_diff),
                                           "estimate 95 % CI high": a1 - a2 + stats.norm.ppf(0.975) * np.sqrt(var_diff)},
             "paired DeLong (1988) written out, NumPy")
for interval in (90,):
    found90 = keyed(table_of(patients.compare_accuracy, "s100b", "ndka", truth="outcome", positive="Poor",
                       interval=interval))
    row = found90["difference in AUC (s100b minus ndka)"]
    lo_c = [c for c in columns_of(found90) if "CI low" in c][0]
    hi_c = [c for c in columns_of(found90) if "CI high" in c][0]
    compare_many(ds, f"compare_accuracy (interval={interval})", row,
                 {lo_c: a1 - a2 - stats.norm.ppf(0.95) * np.sqrt(var_diff),
                  hi_c: a1 - a2 + stats.norm.ppf(0.95) * np.sqrt(var_diff)}, "paired DeLong (1988) written out, NumPy")
if HAS_R:
    r = r_run("""suppressPackageStartupMessages(library(pROC))
        r1 <- roc(outcome, s100b, levels = c("Good", "Poor"), direction = "<", quiet = TRUE)
        r2 <- roc(outcome, ndka, levels = c("Good", "Poor"), direction = "<", quiet = TRUE)
        t <- roc.test(r1, r2, method = "delong", paired = TRUE)
        list(statistic = unname(t$statistic), `p-value` = t$p.value)""", asah[["outcome", "s100b", "ndka"]])
    compare_many(ds, "compare_accuracy", found["difference in AUC (s100b minus ndka)"], r,
                 "R pROC roc.test (help page example)")


made = rd.r_example_objects("epiR", "epi.tests")
tp_, fp_, fn_, tn_ = (int(v) for v in r_run("list(v = dat.v01)")["v"])
text = rd.r_example_code("epiR", "epi.tests")
sens_text = re.search(r"sensitivity is ([\d.]+) \(95% CI ([\d.]+) to ([\d.]+)\)", text).groups()
spec_text = re.search(r"specificity is\s+##\s*([\d.]+) \(95% CI ([\d.]+) to ([\d.]+)\)|specificity is ([\d.]+) \(95% CI ([\d.]+) to ([\d.]+)\)",
                      text.replace("\n", " ")).groups()
spec_text = [v for v in spec_text if v]
lr_text = re.search(r"likelihood ratio of a positive test\s+##\s*is ([\d.]+) \(95% CI ([\d.]+) to ([\d.]+)\)|likelihood ratio of a positive test is ([\d.]+) \(95% CI ([\d.]+) to ([\d.]+)\)",
                    text.replace("\n", " ")).groups()
lr_text = [v for v in lr_text if v]
print("from the help page:", sens_text, spec_text, lr_text)
rows = ([("yes", "yes")] * tp_ + [("yes", "no")] * fp_ + [("no", "yes")] * fn_ + [("no", "no")] * tn_)
scott_rows = pd.DataFrame(rows, columns=["test", "truth"])
scott = sdata.table(scott_rows)
ds = "Scott et al. (2008), R epiR epi.tests example"
try:
    from_rows = keyed(table_of(scott.compare_to_truth, "test", truth="truth", positive="yes"))
except Exception as problem:  # noqa: BLE001
    print("positive='has it' refused:", problem)
    from_rows = None
from_counts = keyed(table_of(sdata.compare_to_truth, found=tp_, missed=fn_, false_alarms=fp_, correct_negatives=tn_,
                       prevalence=5))
lr_pos = (tp_ / (tp_ + fn_)) / (fp_ / (fp_ + tn_))
lr_neg = (fn_ / (tp_ + fn_)) / (tn_ / (fp_ + tn_))
dor = tp_ * tn_ / (fp_ * fn_)
se_pos = np.sqrt(1 / tp_ - 1 / (tp_ + fn_) + 1 / fp_ - 1 / (fp_ + tn_))
se_neg = np.sqrt(1 / fn_ - 1 / (tp_ + fn_) + 1 / tn_ - 1 / (fp_ + tn_))
se_dor = np.sqrt(1 / tp_ + 1 / fp_ + 1 / fn_ + 1 / tn_)
zq = stats.norm.ppf(0.975)
sens, spec = tp_ / (tp_ + fn_), tn_ / (tn_ + fp_)
for label, found in [("sdata.compare_to_truth (four counts, prevalence=5)", from_counts),
                     ("table.compare_to_truth (a yes/no test column)", from_rows)]:
    if found is None:
        checks.skip(ds, label, "every measure", "R / formulas", "the word refused the columns")
        continue
    for measure, (k_, n_) in {"sensitivity": (tp_, tp_ + fn_), "specificity": (tn_, tn_ + fp_),
                              "positive predictive value": (tp_, tp_ + fp_),
                              "negative predictive value": (tn_, tn_ + fn_)}.items():
        if HAS_R:
            r = r_run(f"b <- binom.test({k_}, {n_}); list(estimate = {k_} / {n_}, `estimate 95 % CI low` = b$conf.int[1], `estimate 95 % CI high` = b$conf.int[2])")
            compare_many(ds, label, found[measure], r, "R binom.test", prefix=measure)
    for measure, value, se_log in [("positive likelihood ratio", lr_pos, se_pos), ("negative likelihood ratio", lr_neg, se_neg),
                                   ("diagnostic odds ratio", dor, se_dor)]:
        if measure in found:
            compare_many(ds, label, found[measure], {"estimate": value,
                                                         "estimate 95 % CI low": value * np.exp(-zq * se_log),
                                                         "estimate 95 % CI high": value * np.exp(zq * se_log)},
                         "Simel, Samsa & Matchar (1991) log intervals, NumPy", prefix=measure)
    for measure, shown in [("sensitivity", sens_text), ("specificity", spec_text), ("positive likelihood ratio", lr_text)]:
        for column, s in zip(["estimate", "estimate 95 % CI low", "estimate 95 % CI high"], shown):
            checks.printed(ds, label, f"{measure}: {column}", s, "R epiR epi.tests help page text",
                           title="epiR, epi.tests, example 1 (Scott et al. 2008, Table 1)",
                           url="https://search.r-project.org/CRAN/refmans/epiR/html/epi.tests.html")
ppv5 = [m for m in from_counts if "positive predictive value where" in m][0]
npv5 = [m for m in from_counts if "negative predictive value where" in m][0]
checks.compare(ds, "sdata.compare_to_truth (four counts, prevalence=5)", ppv5, from_counts[ppv5]["estimate"],
               sens * 0.05 / (sens * 0.05 + (1 - spec) * 0.95), "Bayes' rule, NumPy", CLOSED, "relative")
checks.compare(ds, "sdata.compare_to_truth (four counts, prevalence=5)", npv5, from_counts[npv5]["estimate"],
               spec * 0.95 / (spec * 0.95 + (1 - sens) * 0.05), "Bayes' rule, NumPy", CLOSED, "relative")
for interval in (90,):
    found = keyed(table_of(sdata.compare_to_truth, found=tp_, missed=fn_, false_alarms=fp_, correct_negatives=tn_,
                     interval=interval))
    lo_c = [c for c in columns_of(found) if "CI low" in c][0]
    if HAS_R:
        r = r_run(f"b <- binom.test({tp_}, {tp_ + fn_}, conf.level = 0.9); list(low = b$conf.int[1])")
        checks.compare(ds, "sdata.compare_to_truth (interval=90)", "sensitivity 90 % CI low", found["sensitivity"][lo_c],
                       r["low"], "R binom.test", CLOSED, "relative")
found = keyed(table_of(sdata.compare_to_truth, sensitivity=90, specificity=95, prevalence=1))
checks.compare("sensitivity 90 %, specificity 95 %, prevalence 1 %", "sdata.compare_to_truth (sensitivity=, specificity=, prevalence=)",
               "positive predictive value", found["positive predictive value"]["estimate"],
               0.9 * 0.01 / (0.9 * 0.01 + 0.05 * 0.99), "Bayes' rule, NumPy", CLOSED, "relative")
checks.compare("sensitivity 90 %, specificity 95 %, prevalence 1 %", "sdata.compare_to_truth (sensitivity=, specificity=, prevalence=)",
               "negative predictive value", found["negative predictive value"]["estimate"],
               0.95 * 0.99 / (0.95 * 0.99 + 0.1 * 0.01), "Bayes' rule, NumPy", CLOSED, "relative")


checks.data("diagnoses", diagnoses)
checks.data("anxiety", anxiety)
checks.data("anx_text", anx_text)
checks.data("asah", asah)
checks.data("scott_rows", scott_rows)
checks.save()
