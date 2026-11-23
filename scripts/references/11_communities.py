"""Reference values for notebooks/11_communities.ipynb, saved to references/11_communities.csv.

Printed values first: the help pages of vegan `diversity` and `rarefy` (the first three BCI plots) and Oksanen's
vegan tutorial (correspondence analysis of varespec: eigenvalues 0.5249 and 0.3568 of a total inertia 2.083), within
half a unit of the last printed digit. Then the textbook formulas written out (SciPy entropy, pdist, linkage, NumPy):
Shannon and Simpson, rarefaction (Hurlbert 1971), Bray-Curtis, PERMANOVA's pseudo-F and R squared (Anderson 2001),
ANOSIM's R (Clarke 1993), the Mantel correlation, principal coordinates (Gower 1966), the log ratios. Then R
(vegan, zCompositions, compositions, cluster), run once with the same calls; vegan's permutation tests after
set.seed(1), as simplydata runs them, so the permutation p-values must agree exactly.

Tolerances: 1e-12 relative for closed formulas, 1e-10 for decompositions, 1e-6 for the iterative NMDS scores and
the log ratios after the iterative zero replacement; permutation p-values with the same seed 1e-12. Ordination
axes have an arbitrary sign: their scores are compared by absolute value. Whole tables (the 190 Bray-Curtis
distances, the 70 x 35 log ratios) are compared by their largest relative difference.

The data are bci.csv, dune.csv, varespec.csv and mite.csv, written to data/fetched/ from the R package vegan by
scripts/get_data.py (they are not stored in this repository).

    python scripts/references/11_communities.py      (from the repository root; needs R)
"""
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import Recorder  # noqa: E402

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.cluster.hierarchy import linkage  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402

from checks import load_csv, r_available, r_run  # noqa: E402

STEM = "11_communities"
R_PACKAGES = ("vegan", "zCompositions", "compositions", "cluster", "permute")
checks = Recorder(STEM, R_PACKAGES)
HAS_R = r_available()
FORMULA = 1e-12
DECOMPOSED = 1e-10
ITERATIVE = 1e-6
SAME_SEED = 1e-12
VEGAN = "https://search.r-project.org/CRAN/refmans/vegan/html/{}.html"
TUTORIAL = "https://john-quensen.com/wp-content/uploads/2018/10/Oksanen_Jari_vegantutor.pdf"


def rel(dataset, method, quantity, value, source, tolerance):
    checks.compare(dataset, method, quantity, None, value, source, tolerance, "relative")


def absolute(dataset, method, quantity, value, source, tolerance):
    checks.compare(dataset, method, quantity, None, value, source, tolerance)


def printed(dataset, method, quantity, value, source, tolerance, url):
    checks.compare(dataset, method, quantity, None, value, source, tolerance)
    checks.rows[-1].update({"url": url, "kind": "printed", "digits as printed": repr(value)})


bci = load_csv("bci.csv")
dune = load_csv("dune.csv")
DUNE_ENV = ["plot", "A1", "Moisture", "Management", "Use", "Manure"]
dune_species = [c for c in dune.columns if c not in DUNE_ENV]
vare = load_csv("varespec.csv")
VARECHEM = ["N", "P", "K", "Ca", "Mg", "S", "Al", "Fe", "Mn", "Zn", "Mo", "Baresoil", "Humdepth", "pH"]
vare_species = [c for c in vare.columns if c not in VARECHEM + ["site"]]
mite = load_csv("mite.csv")
mite_species = [c for c in mite.columns if c != "Topo"]


def r_with(table, columns, code, extra=None):
    """Run R code with the counts as matrix x (and extra columns as vectors)."""
    data = {"v": table[columns].to_numpy(dtype=float).ravel(order="F")}
    if extra:
        data.update(extra)
    n, k = len(table), len(columns)
    return r_run(f"x <- matrix(v, nrow = {n}, ncol = {k}); colnames(x) <- paste0('s', seq_len({k}))\n" + code, data)


# ---- diversity of BCI
species = list(bci.columns)
counts = bci[species].to_numpy(dtype=float)
share = counts / counts.sum(axis=1, keepdims=True)
formula = {"richness": (counts > 0).sum(axis=1), "shannon": stats.entropy(counts, axis=1),
           "simpson": 1 - (share ** 2).sum(axis=1)}
r = r_with(bci, species, "list(richness = vegan::specnumber(x), shannon = vegan::diversity(x), "
                         "simpson = vegan::diversity(x, 'simpson'))") if HAS_R else None
for i in range(len(bci)):
    for column in ("richness", "shannon", "simpson"):
        rel(f"BCI plot {i + 1}", "add_diversity_columns", column, formula[column][i], "formula (SciPy entropy, NumPy)",
            FORMULA)
        if r is not None:
            rel(f"BCI plot {i + 1}", "add_diversity_columns", column, r[column][i], "R vegan", FORMULA)
for i, (h, d) in enumerate(((4.018, 0.975), (3.848, 0.968), (3.814, 0.965))):
    printed(f"BCI plot {i + 1}", "add_diversity_columns", "shannon", h, "vegan diversity help page, printed", 5e-4,
            VEGAN.format("diversity"))
    printed(f"BCI plot {i + 1}", "add_diversity_columns", "simpson", d, "vegan diversity help page, printed", 5e-4,
            VEGAN.format("diversity"))

hurlbert = []
for plot in counts.astype(int):
    total = int(plot.sum())
    hurlbert.append(float(sum(1 - Fraction(comb(total - k, 20), comb(total, 20)) for k in plot if k > 0)))
r = r_with(bci, species, "list(rarefied = as.numeric(vegan::rarefy(x, sample = 20)))") if HAS_R else None
M = "add_diversity_columns (at_total_count=20)"
for i in range(len(bci)):
    rel(f"BCI plot {i + 1}", M, "expected richness", hurlbert[i], "Hurlbert (1971) formula", FORMULA)
    if r is not None:
        rel(f"BCI plot {i + 1}", M, "expected richness", r["rarefied"][i], "R vegan::rarefy", FORMULA)
for i, value in enumerate((16.361, 15.670, 15.301)):
    printed(f"BCI plot {i + 1}", M, "expected richness", value, "vegan rarefy help page, printed", 5e-4,
            VEGAN.format("rarefy"))

if HAS_R:
    r = r_with(bci, species, "e <- vegan::estimateR(x); h <- vegan::renyi(x, scales = c(1, 2), hill = TRUE)\n"
                             "list(chao = e[2, ], se = e[3, ], q1 = h[, 1], q2 = h[, 2])")
    for i in range(len(bci)):
        dataset = f"BCI plot {i + 1}"
        rel(dataset, "add_diversity_columns (chao1)", "chao1", r["chao"][i], "R vegan::estimateR", FORMULA)
        rel(dataset, "add_diversity_columns (chao1)", "chao1 standard error", r["se"][i], "R vegan::estimateR", FORMULA)
        rel(dataset, "add_diversity_columns (hill numbers)", "hill number 1", r["q1"][i], "R vegan::renyi", FORMULA)
        rel(dataset, "add_diversity_columns (hill numbers)", "hill number 2", r["q2"][i], "R vegan::renyi", FORMULA)
for i in range(len(bci)):
    rel(f"BCI plot {i + 1}", "add_diversity_columns (hill numbers)", "hill number 1", np.exp(stats.entropy(counts[i])),
        "formula exp(Shannon)", FORMULA)
    rel(f"BCI plot {i + 1}", "add_diversity_columns (hill numbers)", "hill number 2", 1 / (share[i] ** 2).sum(),
        "formula 1 / sum p²", FORMULA)

# ---- Bray-Curtis between the dune meadows
x = dune[dune_species].to_numpy(dtype=float)
scipy_d = squareform(pdist(x, "braycurtis"))
upper = np.triu_indices(len(x), 1)
M = "find_distances_between_rows"
checks.vector("dune", M, "distance of each pair", scipy_d[upper], "SciPy pdist (braycurtis)", FORMULA)
for i, j in ((0, 1), (0, 19), (7, 15)):
    rel("dune", M, f"meadow {i + 1} to {j + 1}", scipy_d[i, j], "SciPy pdist (braycurtis)", FORMULA)
if HAS_R:
    r = r_with(dune, dune_species, "list(d = as.numeric(as.matrix(vegan::vegdist(x))))")
    vegan_d = np.asarray(r["d"]).reshape(len(x), len(x), order="F")
    checks.vector("dune", M, "distance of each pair", vegan_d[upper], "R vegan::vegdist", FORMULA)


# ---- PERMANOVA, spread check, ANOSIM
def anderson(d, groups):
    groups = np.asarray(groups)
    n, labels = len(groups), np.unique(groups)
    total = (d[np.triu_indices(n, 1)] ** 2).sum() / n
    within = 0.0
    for g in labels:
        idx = np.flatnonzero(groups == g)
        sub = d[np.ix_(idx, idx)]
        within += (sub[np.triu_indices(len(idx), 1)] ** 2).sum() / len(idx)
    a = len(labels)
    return ((total - within) / (a - 1)) / (within / (n - a)), (total - within) / total


def clarke_r(d, groups):
    groups = np.asarray(groups)
    i, j = np.triu_indices(len(groups), 1)
    ranks = stats.rankdata(d[i, j])
    same = groups[i] == groups[j]
    return (ranks[~same].mean() - ranks[same].mean()) / (len(ranks) / 2)


for name, table, columns, group in (("dune", dune, dune_species, "Management"), ("mite", mite, mite_species, "Topo")):
    M = "compare_profiles_of (PERMANOVA)"
    d = squareform(pdist(table[columns].to_numpy(dtype=float), "braycurtis"))
    f, r2 = anderson(d, table[group])
    rel(name, M, "pseudo-F", f, "Anderson (2001) formula", DECOMPOSED)
    rel(name, M, "R squared", r2, "Anderson (2001) formula", DECOMPOSED)
    if HAS_R:
        r = r_with(table, columns, "set.seed(1); d <- vegan::vegdist(x); g <- factor(g)\n"
                   "a <- vegan::adonis2(d ~ g, permutations = 999)\n"
                   "b <- vegan::betadisper(d, g); t <- vegan::permutest(b, permutations = 999)\n"
                   "list(f = a$F[1], r2 = a$R2[1], p = a[['Pr(>F)']][1], sf = t$tab$F[1], sp = t$tab[['Pr(>F)']][1])",
                   {"g": table[group].astype(str)})
        rel(name, M, "pseudo-F", r["f"], "R vegan::adonis2", DECOMPOSED)
        rel(name, M, "R squared", r["r2"], "R vegan::adonis2", DECOMPOSED)
        absolute(name, M, "p-value (seed 1)", r["p"], "R vegan::adonis2", SAME_SEED)
        rel(name, M, "spread check F", r["sf"], "R vegan::betadisper + permutest", DECOMPOSED)
        absolute(name, M, "spread check p-value (seed 1)", r["sp"], "R vegan::betadisper + permutest", SAME_SEED)
    M = "compare_profiles_of (ANOSIM)"
    rel(name, M, "ANOSIM R", clarke_r(d, table[group]), "Clarke (1993) formula", FORMULA)
    if HAS_R:
        r = r_with(table, columns, "set.seed(1); a <- vegan::anosim(vegan::vegdist(x), factor(g), permutations = 999)\n"
                   "list(stat = a$statistic, p = a$signif)", {"g": table[group].astype(str)})
        rel(name, M, "ANOSIM R", r["stat"], "R vegan::anosim", FORMULA)
        absolute(name, M, "p-value (seed 1)", r["p"], "R vegan::anosim", SAME_SEED)

# ---- principal coordinates, NMDS, correspondence analysis
d = squareform(pdist(dune[dune_species].to_numpy(dtype=float), "braycurtis"))
n = len(d)
centre = np.eye(n) - np.ones((n, n)) / n
values, vectors = np.linalg.eigh(-0.5 * centre @ (d ** 2) @ centre)
order = np.argsort(values)[::-1]
gower = vectors[:, order[:2]] * np.sqrt(values[order[:2]])
for k in (0, 1):
    for i in (0, 9, 19):
        rel("dune", "add_principal_coordinates", f"|PCo{k + 1}| meadow {i + 1}", abs(gower[i, k]),
            "Gower (1966), NumPy eigh", DECOMPOSED)
if HAS_R:
    r = r_with(dune, dune_species, "p <- stats::cmdscale(vegan::vegdist(x), k = 2); list(p1 = p[, 1], p2 = p[, 2])")
    for k in (0, 1):
        for i in (0, 9, 19):
            rel("dune", "add_principal_coordinates", f"|PCo{k + 1}| meadow {i + 1}", abs(r[f"p{k + 1}"][i]),
                "R stats::cmdscale", DECOMPOSED)

    for name, table, columns in (("dune", dune, dune_species), ("varespec", vare, vare_species)):
        r = r_with(table, columns, "set.seed(1); d <- vegan::vegdist(x)\n"
                   "f <- vegan::metaMDS(d, k = 2, try = 20, trymax = 100, autotransform = FALSE, trace = 0)\n"
                   "list(stress = f$stress, p1 = f$points[, 1], p2 = f$points[, 2])")
        absolute(name, "add_nmds", "stress (printed, 3 decimals)", round(r["stress"], 3), "R vegan::metaMDS (seed 1)", 5e-4)
        for k in (1, 2):
            rel(name, "add_nmds", f"|NMDS{k}| row 1", abs(r[f"p{k}"][0]), "R vegan::metaMDS (seed 1)", ITERATIVE)

y = vare[vare_species].to_numpy(dtype=float)
p = y / y.sum()
rows, cols = p.sum(axis=1), p.sum(axis=0)
residuals = (p - np.outer(rows, cols)) / np.sqrt(np.outer(rows, cols))
singular = np.linalg.svd(residuals, compute_uv=False)
inertia = (singular ** 2).sum()
M = "add_correspondence_axes"
for k, value in enumerate((25.2, 17.1)):
    printed("varespec", M, f"CA{k + 1} % (printed)", value, "Oksanen (2015) vegan tutorial, printed", 0.05, TUTORIAL)
printed("varespec", M, "total inertia", 2.083, "Oksanen (2015) vegan tutorial, printed", 5e-4, TUTORIAL)
for k, eigen in enumerate((0.5249, 0.3568)):
    printed("varespec", M, f"CA{k + 1} eigenvalue", eigen, "Oksanen (2015) vegan tutorial, printed", 5e-5, TUTORIAL)
    rel("varespec", M, f"CA{k + 1} %", 100 * singular[k] ** 2 / inertia, "chi-square residuals, NumPy svd", DECOMPOSED)
if HAS_R:
    r = r_with(vare, vare_species, "m <- vegan::cca(x)\n"
               "s <- vegan::scores(m, display = 'sites', scaling = 'sites', choices = 1:2)\n"
               "list(c1 = s[, 1], c2 = s[, 2])")
    for k in (1, 2):
        for i in (0, 11, 23):
            rel("varespec", M, f"|CA{k}| site {i + 1}", abs(r[f"c{k}"][i]), "R vegan::cca (scaling sites)", DECOMPOSED)

# ---- db-RDA, CCA, Mantel
chem = {c: vare[c] for c in VARECHEM}
if HAS_R:
    M = "explain_profiles_of (db-RDA)"
    r = r_with(vare, vare_species, "env <- data.frame(N = N, P = P, K = K)\n"
               "m <- vegan::dbrda(x ~ N + P + K, data = env, distance = 'bray'); q <- vegan::RsquareAdj(m)\n"
               "set.seed(1); a <- anova(m, by = 'margin', permutations = 999)\n"
               "list(r2 = q$r.squared, adj = q$adj.r.squared, fN = a['N', 'F'], fP = a['P', 'F'], fK = a['K', 'F'], "
               "pN = a['N', 4], pP = a['P', 4], pK = a['K', 4])", chem)
    rel("varespec", M, "R squared", r["r2"], "R vegan::dbrda", DECOMPOSED)
    rel("varespec", M, "adjusted R squared", r["adj"], "R vegan::dbrda", DECOMPOSED)
    for term in ("N", "P", "K"):
        rel("varespec", M, f"pseudo-F of {term} (marginal)", r[f"f{term}"], "R vegan::dbrda anova(by = 'margin')",
            DECOMPOSED)
        absolute("varespec", M, f"p-value of {term} (seed 1)", r[f"p{term}"], "R vegan::dbrda anova(by = 'margin')",
                 SAME_SEED)
    M = "explain_profiles_of (CCA)"
    r = r_with(vare, vare_species, "env <- data.frame(Al = Al, P = P, K = K)\n"
               "m <- vegan::cca(x ~ Al + P + K, data = env)\n"
               "list(r2 = vegan::RsquareAdj(m)$r.squared, f = anova(m, permutations = 9)$F[1])", chem)
    rel("varespec", M, "R squared", r["r2"], "R vegan::cca", DECOMPOSED)
    rel("varespec", M, "pseudo-F", r["f"], "R vegan::cca", DECOMPOSED)

M = "find_correlation_between_profiles (Mantel)"
standardised = (vare[VARECHEM] - vare[VARECHEM].mean()) / vare[VARECHEM].std(ddof=1)
pearson = stats.pearsonr(pdist(vare[vare_species].to_numpy(dtype=float), "braycurtis"),
                         pdist(standardised.to_numpy(dtype=float), "euclidean")).statistic
rel("varespec and varechem", M, "Mantel r", pearson, "Pearson r of the distances (SciPy)", FORMULA)
if HAS_R:
    r = r_with(vare, vare_species, "env <- data.frame(" + ", ".join(VARECHEM) + ")\n"
               "set.seed(1); t <- vegan::mantel(vegan::vegdist(x), stats::dist(scale(env)), permutations = 999)\n"
               "list(stat = t$statistic, p = t$signif)", chem)
    rel("varespec and varechem", M, "Mantel r", r["stat"], "R vegan::mantel", FORMULA)
    absolute("varespec and varechem", M, "p-value (seed 1)", r["p"], "R vegan::mantel", SAME_SEED)


# ---- clusters and log ratios
def numbered(labels):
    """Cluster labels renumbered 1, 2, 3 ... in the order they first appear (the numbering itself is arbitrary)."""
    order = {}
    return [order.setdefault(v, len(order) + 1) for v in labels]


x = dune[dune_species].to_numpy(dtype=float)
heights = np.sort(linkage(pdist(x, "braycurtis"), method="average")[:, 2])
T = "draw_cluster_tree (average, Bray-Curtis)"
rel("dune", T, "top join", heights[-1], "SciPy linkage (average)", FORMULA)
rel("dune", T, "lowest join", heights[0], "SciPy linkage (average)", FORMULA)
if HAS_R:
    r = r_with(dune, dune_species, "d <- vegan::vegdist(x); h <- stats::hclust(d, 'average')\n"
               "s <- cluster::silhouette(stats::cutree(h, 3), d)\n"
               "list(top = max(h$height), low = min(h$height), avg = summary(s)$avg.width, w = s[, 3], "
               "g = stats::cutree(h, 3))")
    rel("dune", T, "top join", r["top"], "R stats::hclust", FORMULA)
    M = "add_clusters (hierarchical, 3)"
    rel("dune", M, "average silhouette width", r["avg"], "R cluster::silhouette", FORMULA)
    for i in (0, 9, 19):
        rel("dune", M, f"silhouette width of meadow {i + 1}", r["w"][i], "R cluster::silhouette", FORMULA)
    checks.vector("dune", M, "cluster of each meadow (numbered in order of first appearance)", numbered(r["g"]),
                  "R stats::cutree", 0)

logs = np.log(mite[mite_species].to_numpy(dtype=float) + 1)
clr = logs - logs.mean(axis=1, keepdims=True)
M = "add_log_ratio_columns (adding 1)"
checks.vector("mite", M, "clr of each value (row by row)", clr.ravel(), "formula log(x + 1) - row mean", FORMULA)
if HAS_R:
    r = r_with(mite, mite_species, "list(m = as.numeric(as.matrix(vegan::decostand(x, 'clr', pseudocount = 1))))")
    vegan_clr = np.asarray(r["m"]).reshape(clr.shape, order="F")
    checks.vector("mite", M, "clr of each value (row by row)", vegan_clr.ravel(),
                  "R vegan::decostand(clr, pseudocount = 1)", FORMULA)
    r = r_with(mite, mite_species, "r <- as.matrix(zCompositions::cmultRepl(x, method = 'GBM', output = 'p-counts', "
               "z.delete = FALSE, suppress.print = TRUE))\n"
               "list(m = as.numeric(unclass(compositions::clr(r))))")
    gbm_clr = np.asarray(r["m"]).reshape(clr.shape, order="F")
    for i, j in ((0, 0), (1, 34), (69, 0), (35, 17)):
        rel("mite", "add_log_ratio_columns (estimated small values)", f"clr row {i + 1}, {mite_species[j]}",
            gbm_clr[i, j], "R zCompositions::cmultRepl (GBM) + compositions::clr", ITERATIVE)

from links import add_links  # noqa: E402  (links for the sources common.url_of cannot name)

add_links(checks)
checks.save()
