"""Reference values for notebooks/28_patterns_across_columns.ipynb.

Printed first: Rosseel (2012) for the political democracy model. Elsewhere nothing printed matches these checks, so
the help-page example of the library simplydata calls is run once (R prcomp, kmeans, hclust, cluster, vegan,
zCompositions, compositions, MASS lda, psych, lavaan) and saved with the link to its page. Also saves the data the
notebook reads (the R datasets).
"""

import numpy as np
import pandas as pd

from _tools import References, load_csv, r_available, r_dataset, r_run, r_version, rman

refs = References("28_patterns_across_columns")
HAS_R = r_available()
CLOSED, ITERATIVE = 1e-10, 1e-6

arrests = r_dataset("USArrests")
refs.data("USArrests", arrests)
COLS = ["Murder", "Assault", "UrbanPop", "Rape"]
varespec = r_dataset("varespec", "vegan")
refs.data("varespec", varespec)
SPECIES = [c for c in varespec.columns if c != "row"]
varechem = r_dataset("varechem", "vegan")
both = varespec.merge(varechem, on="row")
refs.data("varespec_and_varechem", both)
dune = r_dataset("dune", "vegan").drop(columns=["row"], errors="ignore")
dune_env = r_dataset("dune.env", "vegan").drop(columns=["row"], errors="ignore")
DUNE_ALL = list(dune.columns)
DUNE = [c for c in DUNE_ALL if (dune[c] > 0).sum() >= 2]
refs.data("dune", dune)
refs.data("dune_and_env", pd.concat([dune, dune_env], axis=1))
# two columns: a one-column file is misread by sdata.table (FINDINGS A24)
refs.data("DUNE_27_species", pd.DataFrame({"number": range(1, len(DUNE) + 1), "species": DUNE}))
hydro = r_dataset("Hydrochem", "compositions")
refs.data("Hydrochem", hydro)
PARTS = ["H", "Na", "K", "Mg", "Ca", "Sr", "Ba", "NH4", "Cl", "NO3", "PO4", "SO4", "HCO3", "TOC"]
politics = r_dataset("PoliticalDemocracy", "lavaan")
refs.data("PoliticalDemocracy", politics)
bci = r_dataset("BCI", "vegan")
refs.data("BCI", bci)
hs = load_csv("holzinger_swineford_1939.csv")
iris = load_csv("iris.csv")

if HAS_R:
    D = "USArrests"
    # 1. principal components: the help-page example prcomp(USArrests, scale. = TRUE)
    r = r_run('''p <- prcomp(data.frame(Murder, Assault, UrbanPop, Rape), scale. = TRUE)
                 list(sdev = p$sdev, share = 100 * p$sdev^2 / sum(p$sdev^2), rot = as.vector(p$rotation), x = as.vector(p$x[, 1:2]))''', arrests)
    url = rman("stats", "prcomp")
    refs.run_vector(D, "find_principal_components_of", "standard deviations", r["sdev"], "R stats prcomp (help page example)", url, "stats", CLOSED)
    refs.run_vector(D, "find_principal_components_of", "spread explained (%)", r["share"], "R stats prcomp (help page example)", url, "stats", CLOSED)
    refs.run_matrix(D, "find_principal_components_of", "weights", np.reshape(r["rot"], (4, 4), order="F"), "R stats prcomp (help page example)",
                    url, "stats", CLOSED, "signs")
    refs.run_matrix(D, "add_principal_components (how_many=2)", "scores", np.reshape(r["x"], (50, 2), order="F"), "R stats prcomp (help page example)",
                    url, "stats", CLOSED, "signs")
    # 2. clusters
    for seed in (1, 2):
        r = r_run(f"z <- scale(data.frame(Murder, Assault, UrbanPop, Rape)); set.seed({seed}); list(cl = kmeans(z, centers = 4, nstart = 25)$cluster)",
                  arrests)
        refs.run_vector(D, f'add_clusters (test_method="k-means", how_many=4, with_seed={seed})', "cluster", r["cl"],
                        "R stats kmeans (nstart = 25, same seed)", rman("stats", "kmeans"), "stats", 0, "groups")
    r = r_run('''z <- scale(data.frame(Murder, Assault, UrbanPop, Rape)); list(cl = cutree(hclust(dist(z), method = "ward.D2"), k = 4))''', arrests)
    refs.run_vector(D, 'add_clusters (test_method="hierarchical", linkage="ward", how_many=4)', "cluster", r["cl"], "R stats hclust (ward.D2) + cutree",
                    rman("stats", "hclust"), "stats", 0, "groups")
    r = r_run('''z <- scale(data.frame(Murder, Assault, UrbanPop, Rape)); set.seed(1); k <- kmeans(z, centers = 4, nstart = 25)
                 s <- cluster::silhouette(k$cluster, dist(z)); list(w = s[, 3])''', arrests, packages=("cluster",))
    refs.run_vector(D, 'add_clusters (how_many=4, show="silhouette")', "silhouette widths", r["w"],
                    "R cluster silhouette (of R kmeans with the same seed)", rman("cluster", "silhouette"), "cluster", CLOSED)
    x = f"x <- data.frame({', '.join(SPECIES)})\n"
    r = r_run(x + 'list(cl = cutree(hclust(vegan::vegdist(x), method = "average"), k = 3))', varespec, packages=("vegan",))
    refs.run_vector("varespec (vegan)", 'add_clusters (test_method="hierarchical", measure="bray-curtis", how_many=3)', "cluster", r["cl"],
                    "R hclust (average) on vegan vegdist", rman("vegan", "vegdist"), "vegan", 0, "groups")
    # 3. distances (vegdist help-page data)
    n = len(varespec)
    for measure, call in (("bray-curtis", "vegan::vegdist(x)"), ("jaccard", 'vegan::vegdist(x, method = "jaccard", binary = TRUE)')):
        r = r_run(x + f"list(d = as.vector(as.matrix({call})))", varespec, packages=("vegan",))
        refs.run_matrix("varespec (vegan)", f'find_distances_between_rows (measure="{measure}")', "distances", np.reshape(r["d"], (n, n), order="F"),
                        "R vegan vegdist (help page data)", rman("vegan", "vegdist"), "vegan", CLOSED)
    r = r_run(f'''x <- data.frame({", ".join(DUNE)})
                  z <- zCompositions::cmultRepl(x, method = "GBM", output = "p-counts", z.warning = 1, suppress.print = TRUE)
                  list(d = as.vector(as.matrix(vegan::vegdist(z, method = "aitchison"))))''', dune, packages=("vegan", "zCompositions"))
    refs.run_matrix("dune (vegan), 27 species", 'find_distances_between_rows (measure="aitchison")', "distances",
                    np.reshape(r["d"], (len(dune), len(dune)), order="F"), "R zCompositions cmultRepl (GBM) + vegan vegdist (aitchison)",
                    rman("zCompositions", "cmultRepl"), "zCompositions", CLOSED)
    # 4. maps
    for seed, tries in ((1, 20), (5, 30)):
        r = r_run(x + f'''set.seed({seed}); m <- vegan::metaMDS(vegan::vegdist(x), k = 2, try = {tries}, trymax = 100, trace = 0)
                          list(p = as.vector(m$points), stress = m$stress)''', varespec, packages=("vegan",))
        refs.run_matrix("varespec (vegan)", f"add_nmds (tries={tries}, with_seed={seed})", "map", np.reshape(r["p"], (n, 2), order="F"),
                        "R vegan metaMDS (same seed)", rman("vegan", "metaMDS"), "vegan", ITERATIVE, "signs")
    r = r_run(x + "list(p = as.vector(cmdscale(vegan::vegdist(x), k = 3)))", varespec, packages=("vegan",))
    refs.run_matrix("varespec (vegan)", 'add_principal_coordinates (measure="bray-curtis", how_many=3)', "coordinates", np.reshape(r["p"], (n, 3), order="F"),
                    "R cmdscale on vegan vegdist", rman("stats", "cmdscale"), "stats", CLOSED, "signs")
    r = r_run(f'''x <- data.frame({", ".join(DUNE_ALL)}); m <- vegan::cca(x)
                  list(s = as.vector(vegan::scores(m, display = "sites", choices = 1:3, scaling = 1)))''', dune, packages=("vegan",))
    refs.run_matrix("dune (vegan)", "add_correspondence_axes (how_many=3)", "row axes", np.reshape(r["s"], (len(dune), 3), order="F"),
                    "R vegan cca without constraints (each axis compared at unit length)", rman("vegan", "cca"), "vegan", CLOSED, "unit")
    # 5. LDA (help-page data iris)
    r = r_run('''d <- data.frame(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, Species = factor(Species))
                 list(x = as.vector(predict(MASS::lda(Species ~ ., data = d))$x))''', iris, packages=("MASS",))
    refs.run_matrix("iris", "add_axes_that_separate_groups (how_many=2)", "LD scores", np.reshape(r["x"], (150, 2), order="F"),
                    "R MASS lda (help page data)", rman("MASS", "lda"), "MASS", CLOSED, "signs")
    # 6. log ratios (clr help-page data Hydrochem; dune counts with zeros)
    L = np.log(hydro[PARTS].to_numpy(float))
    r = r_run(f'''x <- data.frame({", ".join(PARTS)})
                  list(clr = as.vector(unclass(compositions::clr(compositions::acomp(x)))),
                       alr = as.vector(unclass(compositions::alr(compositions::acomp(x[, c(setdiff(names(x), "Cl"), "Cl")])))),
                       ilr = as.vector(unclass(compositions::ilr(compositions::acomp(x)))))''', hydro, packages=("compositions",))
    m = len(hydro)
    refs.run_matrix("Hydrochem (compositions)", 'add_log_ratio_columns (test_method="clr")', "log ratios", np.reshape(r["clr"], (m, 14), order="F"),
                    "R compositions clr (help page data)", rman("compositions", "clr"), "compositions", CLOSED, every=10)
    refs.run_matrix("Hydrochem (compositions)", 'add_log_ratio_columns (test_method="alr", against="Cl")', "log ratios",
                    np.reshape(r["alr"], (m, 13), order="F"), "R compositions alr", rman("compositions", "alr"), "compositions", CLOSED, every=10)
    refs.run_matrix("Hydrochem (compositions)", 'add_log_ratio_columns (test_method="ilr")', "log ratios", np.reshape(r["ilr"], (m, 13), order="F"),
                    "R compositions ilr", rman("compositions", "ilr"), "compositions", CLOSED, every=10)
    r = r_run(f'''x <- data.frame({", ".join(DUNE)})
                  z <- zCompositions::cmultRepl(x, method = "GBM", output = "p-counts", z.warning = 1, suppress.print = TRUE)
                  list(clr = as.vector(unclass(compositions::clr(compositions::acomp(z)))),
                       one = as.vector(unclass(compositions::clr(compositions::acomp(x + 1)))))''', dune, packages=("zCompositions", "compositions"))
    refs.run_matrix("dune (vegan), 27 species", 'add_log_ratio_columns (zeros_replaced_by="estimated small values")', "clr",
                    np.reshape(r["clr"], (len(dune), len(DUNE)), order="F"), "R zCompositions cmultRepl (GBM) + compositions clr",
                    rman("zCompositions", "cmultRepl"), "zCompositions", CLOSED)
    refs.run_matrix("dune (vegan), 27 species", 'add_log_ratio_columns (zeros_replaced_by="adding 1")', "clr",
                    np.reshape(r["one"], (len(dune), len(DUNE)), order="F"), "R compositions clr of the counts plus 1",
                    rman("compositions", "clr"), "compositions", CLOSED)
    # 7. db-RDA, RDA on Hellinger shares, CCA (dbrda help page); p-values within 4 Monte Carlo standard errors
    CASES = [("bray-curtis (default)", "vegan::dbrda(x ~ N + P + K, data = e, distance = 'bray')", ((999, 1), (499, 3))),
             ("bray-curtis, allowing_for Al", "vegan::dbrda(x ~ N + P + K + Condition(Al), data = e, distance = 'bray')", ((999, 1),)),
             ("hellinger", "vegan::rda(vegan::decostand(x, 'hellinger') ~ N + P + K, data = e)", ((999, 1),)),
             ("canonical correspondence", "vegan::cca(x ~ N + P + K, data = e)", ((999, 1),))]
    for label, model, runs in CASES:
        for permutations, seed in runs:
            check = f"explain_profiles_of ({label}, permutations={permutations}, with_seed={seed})"
            r = r_run(x + f'''e <- data.frame(N, P, K, Al); m <- {model}
                          set.seed({seed}); a <- anova(m, permutations = {permutations})
                          set.seed({seed}); b <- anova(m, by = "margin", permutations = {permutations}); r2 <- vegan::RsquareAdj(m)
                          set.seed(11); adj <- replicate(30, vegan::RsquareAdj(m)$adj.r.squared)
                          list(F = a$F[1], p = a[["Pr(>F)"]][1], Fm = b$F[1:3], pm = b[["Pr(>F)"]][1:3], r2 = r2$r.squared,
                               ar2 = r2$adj.r.squared, adjm = mean(adj), adjs = sd(adj))''', both, packages=("vegan",))
            url = rman("vegan", "dbrda")
            refs.run_all("varespec and varechem (vegan)", check, {"all together: pseudo-F": r["F"], "R squared": r["r2"]},
                         "R vegan anova and RsquareAdj", url, "vegan", CLOSED)
            if label == "canonical correspondence":   # vegan's adjusted R squared for CCA is itself a permutation estimate
                refs.run("varespec and varechem (vegan)", check, "adjusted R squared", r["adjm"], "R vegan RsquareAdj (mean of 30 runs), 4 MCSE",
                         url, "vegan", 4 * r["adjs"] * np.sqrt(1 + 1 / 30), "absolute")
            else:
                refs.run("varespec and varechem (vegan)", check, "adjusted R squared", r["ar2"], "R vegan RsquareAdj", url, "vegan", CLOSED)
            mc = lambda p: 4 * np.sqrt(max(p * (1 - p), 1 / permutations) / permutations)  # noqa: E731
            refs.run("varespec and varechem (vegan)", check, "all together: p-value", r["p"], "R vegan anova (same seed), 4 MCSE", url, "vegan",
                     mc(r["p"]), "absolute")
            for i, column in enumerate(["N", "P", "K"]):
                refs.run("varespec and varechem (vegan)", check, f"{column} after the others: pseudo-F", r["Fm"][i], "R vegan anova (by = margin)",
                         url, "vegan", CLOSED)
                refs.run("varespec and varechem (vegan)", check, f"{column} after the others: p-value", r["pm"][i],
                         "R vegan anova (by = margin, same seed), 4 MCSE", url, "vegan", mc(r["pm"][i]), "absolute")
    # 8. PERMANOVA (aitchison), SIMPER (help-page example), Mantel and PROTEST
    dd = pd.concat([dune, dune_env], axis=1)
    r = r_run(f'''x <- data.frame({", ".join(DUNE)}); g <- factor(Management)
                  z <- zCompositions::cmultRepl(x, method = "GBM", output = "p-counts", z.warning = 1, suppress.print = TRUE)
                  d <- vegan::vegdist(z, method = "aitchison"); set.seed(1); a <- vegan::adonis2(d ~ g, permutations = 999)
                  list(F = a$F[1], r2 = a$R2[1], p = a[["Pr(>F)"]][1])''', dd, packages=("vegan", "zCompositions"))
    check = 'compare_profiles_of (measure="aitchison", permutations=999, with_seed=1)'
    refs.run_all("dune and dune.env (vegan)", check, {"pseudo-F": r["F"], "R squared": r["r2"]}, "R vegan adonis2 on the same pipeline",
                 rman("vegan", "adonis2"), "vegan", CLOSED)
    refs.run("dune and dune.env (vegan)", check, "p-value", r["p"], "R vegan adonis2 (same seed), 4 MCSE", rman("vegan", "adonis2"), "vegan",
             4 * np.sqrt(max(r["p"] * (1 - r["p"]), 1 / 999) / 999), "absolute")
    r = r_run(f'''x <- data.frame({", ".join(DUNE_ALL)}); s <- vegan::simper(x, factor(Management))
                  out <- list(); for (n in names(s)) {{ k <- s[[n]]; out[[paste0(n, "|avg")]] <- unname(k$average); out[[paste0(n, "|sp")]] <- names(k$average) }}
                  out''', dd, packages=("vegan",))
    for key in [k for k in r if k.endswith("|avg")]:
        pair = key.split("|")[0]
        refs.run_all("dune and dune.env (vegan)", 'compare_profiles_of (show="contributions")',
                     {f"{pair}: {sp}": v for sp, v in zip(r[pair + "|sp"], r[key])}, "R vegan simper (help page example)",
                     rman("vegan", "simper"), "vegan", CLOSED)
    for test_method, permutations, seed in (("mantel", 999, 1), ("mantel", 499, 4), ("procrustes", 999, 1)):
        test = ("m <- vegan::mantel(d1, d2, permutations = P); list(stat = m$statistic, p = m$signif)" if test_method == "mantel" else
                "m <- vegan::protest(cmdscale(d1, k = 2), cmdscale(d2, k = 2), permutations = P); list(stat = m$t0, p = m$signif)")
        r = r_run(x + f"e <- scale(data.frame(N, P, K)); P <- {permutations}; d1 <- vegan::vegdist(x); d2 <- dist(e); set.seed({seed}); {test}",
                  both, packages=("vegan",))
        check = f'find_correlation_between_profiles (test_method="{test_method}", permutations={permutations}, with_seed={seed})'
        url = rman("vegan", "mantel" if test_method == "mantel" else "procrustes")
        refs.run("varespec and varechem (vegan)", check, "statistic", r["stat"], f"R vegan {'mantel' if test_method == 'mantel' else 'protest'}", url,
                 "vegan", CLOSED)
        refs.run("varespec and varechem (vegan)", check, "p-value", r["p"], "R vegan (same seed), 4 MCSE", url, "vegan",
                 4 * np.sqrt(max(r["p"] * (1 - r["p"]), 1 / permutations) / permutations), "absolute")
    # 9. factor analysis (psych) and the SEM (lavaan)
    TESTS = [f"x{i}" for i in range(1, 10)]
    for test_method, rotation, call in (("minres", "oblimin", 'psych::fa(x, nfactors = 3, fm = "minres", rotate = "oblimin")'),
                                        ("ml", "varimax", 'psych::fa(x, nfactors = 3, fm = "ml", rotate = "varimax")')):
        r = r_run(f"x <- data.frame({', '.join(TESTS)}); f <- {call}; l <- unclass(f$loadings)\n"
                  "list(l = as.vector(l[, 1:3]), u = unname(f$uniquenesses))", hs, packages=("psych", "GPArotation"))
        check = f'find_hidden_factors_behind (test_method="{test_method}", rotation="{rotation}", how_many=3)'
        refs.run_matrix("Holzinger & Swineford (1939)", check, "loadings (factors matched)", np.reshape(r["l"], (9, 3), order="F"),
                        f"R psych fa ({test_method}, {rotation})", rman("psych", "fa"), "psych", ITERATIVE, "factors")
        refs.run_vector("Holzinger & Swineford (1939)", check, "own part (uniqueness)", r["u"], f"R psych fa ({test_method}, {rotation})",
                        rman("psych", "fa"), "psych", ITERATIVE)
    # 10. Chao1 (estimateR help page data)
    trees = ", ".join(f"`{c}`" for c in bci.columns)
    r = r_run(f"x <- data.frame({trees}, check.names = FALSE); e <- vegan::estimateR(x)\nlist(chao = unname(e['S.chao1', ]), se = unname(e['se.chao1', ]))",
              bci, packages=("vegan",))
    refs.run_vector("BCI (vegan)", 'add_diversity_columns (show="chao1")', "Chao1", r["chao"], "R vegan estimateR (help page data)",
                    rman("vegan", "estimateR"), "vegan", CLOSED)
    refs.run_vector("BCI (vegan)", 'add_diversity_columns (show="chao1")', "Chao1 standard error", r["se"], "R vegan estimateR (help page data)",
                    rman("vegan", "estimateR"), "vegan", CLOSED)

# 9. Rosseel (2012), section 5: the political democracy model (printed)
ROS = "Rosseel (2012), lavaan: An R Package for Structural Equation Modeling, J. Stat. Softw. 48(2), section 5"
refs.printed_all("PoliticalDemocracy (lavaan)", "check_hidden_factors (formula=)", {"chi-square": "38.125", "degrees of freedom": "35", "p-value": "0.329"},
                 ROS, "https://doi.org/10.18637/jss.v048.i02", "published", "lavaan", r_version("lavaan"))

refs.save()
