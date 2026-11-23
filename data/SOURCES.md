# Data sources

The datasets the notebooks use are of two kinds.

- **Stored in this folder** (the "Datasets stored here" section and the NIST files): public domain data, US
  Government works, plain published facts (measurements and counts published decades ago, or small tables we typed
  from a book or paper), and data under CC0 or another permissive licence. The files are stored byte for byte as
  obtained (`.gitattributes` stops git from changing their line endings).
- **Fetched when the tests run, not stored here** (the "Datasets fetched when the tests run" section): datasets
  whose copies are under GPL or LGPL, under copyright, or carry no licence. `scripts/get_data.py` writes them to
  `data/fetched/` (ignored by git) from their original source: the installed R package that ships each one (read
  with Rscript, no download) or, for the bioequivalence datasets, the publisher's supplementary files. Each entry
  says exactly how the file is made.

`scripts/get_data.py` checks every file, stored or fetched, against `SHA256SUMS` (`--check` does only that) and
stops with a message naming the file and its source if anything differs, so a changed package or download is
never used without notice. `--refresh` writes the fetched files again and downloads the NIST files again from
NIST. The checksums below hold on every system.

## Datasets stored here

### `sleep.csv`

- **What:** Student's sleep data: extra hours of sleep for 10 patients with each of two drugs (20 rows)
- **Origin and citation:** Cushny, A. R. & Peebles, A. R. (1905). The action of optical isomers. II. Hyoscines. Journal of Physiology 32, 501-510. Analysed by Student (1908). The probable error of a mean. Biometrika 6(1), 1-25.
- **Licence:** Public domain (published 1905 and 1908). Copy from R `datasets::sleep`, GPL-2 | GPL-3 (R's datasets package). Stored here: the measurements were published in 1905 and 1908 and are in the public domain; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/sleep.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `c45604d3160d61f40535c8773b1762e120e47b824d95521950ea54ca17516fd9`

### `anscombe.csv`

- **What:** Anscombe's quartet: four sets of 11 (x, y) pairs
- **Origin and citation:** Anscombe, F. J. (1973). Graphs in statistical analysis. The American Statistician 27(1), 17-21.
- **Licence:** The 88 published numbers. Copy from R `datasets::anscombe`, GPL-2 | GPL-3 (R's datasets package). Stored here: the 88 numbers printed in Anscombe (1973) are plain facts, which copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/anscombe.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `51939a326ba12bf4e312d77ca7102948e9e6b3c7b62cc75208dad185917df2c4`

### `iris.csv`

- **What:** Anderson's / Fisher's iris: 150 flowers of 3 species, 4 measurements
- **Origin and citation:** Anderson, E. (1935). The irises of the Gaspe Peninsula. Bulletin of the American Iris Society 59, 2-5. Fisher, R. A. (1936). The use of multiple measurements in taxonomic problems. Annals of Eugenics 7, 179-188.
- **Licence:** Public domain (1935). Copy from R `datasets::iris`, GPL-2 | GPL-3 (R's datasets package). Stored here: the measurements were published in 1935 and 1936 and are in the public domain; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/iris.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `d440daded18634c1da2f05e6b1a30385f2aca6cd38455b31d263e1657260112a`

### `warpbreaks.csv`

- **What:** Warp breaks per loom for 2 wools at 3 tensions (54 rows)
- **Origin and citation:** Tippett, L. H. C. (1950). Technological Applications of Statistics. Wiley.
- **Licence:** Published facts (1950). Copy from R `datasets::warpbreaks`, GPL-2 | GPL-3 (R's datasets package). Stored here: counts of warp breaks published in 1950, plain facts that copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/warpbreaks.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `9317a8e3f68324caf0dfdf8a98d178255c25c66c019211d2ba67feba043638dc`

### `mtcars.csv`

- **What:** Motor Trend car road tests: 32 cars, 11 columns
- **Origin and citation:** Henderson, H. V. & Velleman, P. F. (1981). Building multiple regression models interactively. Biometrics 37, 391-411 (data from Motor Trend, 1974).
- **Licence:** Published 1974 and 1981. Copy from R `datasets::mtcars`, GPL-2 | GPL-3 (R's datasets package). Stored here: car measurements printed in Motor Trend (1974), plain facts that copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/mtcars.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `48c52ce5303f0abe6e662cde775cecdba87e11982941f94d6f78edb994c36ca5`

### `stackloss.csv`

- **What:** Brownlee's stack loss plant data (21 rows)
- **Origin and citation:** Brownlee, K. A. (1965). Statistical Theory and Methodology in Science and Engineering, 2nd ed. Wiley.
- **Licence:** Published 1960/1965. Copy from R `datasets::stackloss`, GPL-2 | GPL-3 (R's datasets package). Stored here: plant operation readings published in 1960/1965, plain facts that copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/stackloss.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `b511cb73f27e84889a29ac19a19423e9f3a1f13053d925c4f8edd1ffd86f5494`

### `puromycin.csv`

- **What:** Puromycin, treated enzyme only: reaction velocity against substrate concentration (12 rows)
- **Origin and citation:** Treloar, M. A. (1974). Effects of Puromycin on Galactosyltransferase of Golgi Membranes. M.Sc. thesis, University of Toronto. Fit printed in Bates, D. M. & Watts, D. G. (1988). Nonlinear Regression Analysis and Its Applications. Wiley.
- **Licence:** Published facts. Copy from R `datasets::Puromycin` (treated rows), GPL-2 | GPL-3 (R's datasets package). Stored here: enzyme reaction rates measured in 1974 and printed in Bates & Watts (1988), plain facts that copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/Puromycin.html
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67), where it was taken from the R copy named.
- **sha256:** `5c645b859135ec80d93cfdd583fe81c3f72cb106c5da5effa25af5dcaf25e063`

### `bliss.csv`

- **What:** Bliss's flour beetles: killed out of exposed at 8 doses of carbon disulphide (log10 mg/l)
- **Origin and citation:** Bliss, C. I. (1935). The calculation of the dosage-mortality curve. Annals of Applied Biology 22, 134-167. Logistic fit printed in Dobson, A. J. & Barnett, A. G. (2008). An Introduction to Generalized Linear Models, 3rd ed. Chapman & Hall/CRC.
- **Licence:** 24 numbers published in 1935, typed from the paper (facts, not covered by copyright).
- **URL:** none (the 1935 paper in Annals of Applied Biology; the numbers are also printed in Dobson & Barnett 2008)
- **Retrieved:** Copied 2026-11-20 from simplydata's `validation/data/` (github.com/bioforkgithub/simplydata, commit 3a26c67).
- **sha256:** `3323e30cf8bd59a8cf894c3dbf54b369a305620247bb36abb11923adf8263b61`

### `columbus.csv`

- **What:** Columbus, Ohio neighbourhoods: centroids X, Y, CRIME, HOVAL, INC (49 rows)
- **Origin and citation:** Anselin, L. (1988). Spatial Econometrics: Methods and Models. Kluwer. Copy from the R package spData (`columbus`).
- **Licence:** The spData copy is CC0. Stored here with that attribution.
- **URL:** https://cran.r-project.org/package=spData
- **Retrieved:** Exported 2026-11-21 from the installed R package, values unchanged.
- **sha256:** `ef1d2adfaec4733e103e277f2d559da3599b7b612df7ad35ea896f209ca08e12`

### `air_passengers.csv`

- **What:** Monthly international airline passengers, 1949-1960
- **Origin and citation:** Box, G. E. P. & Jenkins, G. M. (1976). Time Series Analysis: Forecasting and Control, series G. Copy from R `datasets::AirPassengers`.
- **Licence:** Published facts; the R copy is GPL-2 | GPL-3. Stored here: monthly airline passenger totals printed by Box & Jenkins (1976), plain facts that copyright does not cover; R's GPL covers its packaging of them, not the numbers.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/AirPassengers.html
- **Retrieved:** Exported 2026-11-21 from the installed R package, values unchanged.
- **sha256:** `ca3e3f11e7d2f373fb322753cd417c986959473fcb9d2ad7c55481c3b582c5f0`

### `rounding_times.csv`

- **What:** Times to round first base by three methods for 22 players
- **Origin and citation:** Hollander, M. & Wolfe, D. A. (1973). Nonparametric Statistical Methods. Wiley, p. 140. Typed from the R help page of `friedman.test` (it gives the printed Friedman statistic 11.143).
- **Licence:** Published facts. Stored here with that attribution.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/stats/html/friedman.test.html
- **Retrieved:** Exported 2026-11-21 from the installed R package, values unchanged.
- **sha256:** `6b43373502be149d1bb59c5661d03fd333bd584df1ff405635b59c60744c0ab2`

### `smokers_fleiss.csv`

- **What:** Smokers among patients in four groups
- **Origin and citation:** Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions, 2nd ed. Wiley, p. 139. From the R help page of `prop.test`.
- **Licence:** Published facts. Stored here with that attribution.
- **URL:** https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prop.test.html
- **Retrieved:** Exported 2026-11-21 from the installed R package, values unchanged.
- **sha256:** `4d5df1c480e12bc91e97d9829b711ab7057b71537d2bfc3c89c085034c8d7d9e`

### `poisons.csv`

- **What:** Survival times of animals under 3 poisons and 4 treatments
- **Origin and citation:** Box, G. E. P. & Cox, D. R. (1964). An analysis of transformations. JRSS B 26, 211-252. Copy from the R package boot (`poisons`).
- **Licence:** Published data (facts); the boot copy is Unlimited. Stored here with that attribution.
- **URL:** https://cran.r-project.org/package=boot
- **Retrieved:** Exported 2026-11-21 from the installed R package, values unchanged.
- **sha256:** `b4fa708856eace7e32c327b30cda14fff5b91dbe4dfdb240bc110c5be539e8b5`

### `flutd.csv`

- **What:** Dry cat food and feline lower urinary tract disease (a 2 x 2 table)
- **Origin and citation:** Willeberg, P. (1977). Animal disease information processing: epidemiologic analyses of the feline urologic syndrome. Acta Veterinaria Scandinavica 18 (Suppl. 64), 1-48. The first example of R epiR `epi.2by2`.
- **Licence:** Published counts (facts).
- **URL:** https://cran.r-project.org/package=epiR
- **Retrieved:** Copied 2026-11-21 byte for byte from simplydata's `validation/data/` (commit 37afb81).
- **sha256:** `d67ab174ba4b61989ce0694f23f2edbeb39a9a19ad3985e9c96a36e189e7db50`

### `down_syndrome.csv`

- **What:** Down syndrome births by birth order and maternal age group
- **Origin and citation:** Selvin, S. (2004). Statistical Analysis of Epidemiologic Data, 3rd ed. Oxford University Press; the example of R epitools `ageadjust.direct` (Fay, M. P. & Feuer, E. J. (1997). Statistics in Medicine 16, 791-801).
- **Licence:** Published counts (facts).
- **URL:** https://cran.r-project.org/package=epitools
- **Retrieved:** Copied 2026-11-21 byte for byte from simplydata's `validation/data/` (commit 37afb81).
- **sha256:** `24daf00427be89aa960d61b92f36e2d5ba17d30c83e42c301d0f661ac6b0f3e4`

## Datasets fetched when the tests run (not stored here)

### `birthwt.csv`

- **What:** Risk factors for low birth weight: 189 births, 10 columns
- **Origin and citation:** Hosmer, D. W. & Lemeshow, S. (1989). Applied Logistic Regression. Wiley. Data from Baystate Medical Center, Springfield, Massachusetts, 1986.
- **Licence:** Copy from the R package MASS (`birthwt`), GPL-2 | GPL-3. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=MASS
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/birthwt.csv` from the installed R package MASS (`get_r("birthwt", "MASS")`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `8daf9b523bb6f4e716f933bda42f2f4c5d2208ffaddab852e0d21102e6d8027e`

### `nmes1988.csv`

- **What:** National Medical Expenditure Survey 1987/88: 4406 people aged 66 or more, 7 of the 19 columns (visits, hospital, health, chronic, gender, school, insurance)
- **Origin and citation:** Deb, P. & Trivedi, P. K. (1997). Demand for medical care by the elderly: a finite mixture approach. Journal of Applied Econometrics 12, 313-336. Analysed in Zeileis, A., Kleiber, C. & Jackman, S. (2008). Regression models for count data in R. Journal of Statistical Software 27(8).
- **Licence:** A US Government survey (public domain). Copy from the R package AER (`NMES1988`, GPL-2 | GPL-3) through Rdatasets; the JSS article is CC BY. Not stored here (see the next lines).
- **URL:** https://vincentarelbundock.github.io/Rdatasets/doc/AER/NMES1988.html
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/nmes1988.csv` from the installed R package AER (`get_r("NMES1988", "AER")[, c("visits", "hospital", "health", "chronic", "gender", "school", "insurance")]`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `035267144c01d23ab38696c89b73804c2c7839719b9c46b93948a4ddd2676c16`

### `lung.csv`

- **What:** NCCTG lung cancer: survival in days of 228 patients with advanced lung cancer, 10 columns
- **Origin and citation:** Loprinzi, C. L. et al. (1994). Prospective evaluation of prognostic variables from patient-completed questionnaires. North Central Cancer Treatment Group. Journal of Clinical Oncology 12(3), 601-607. Shipped with the R package survival (Therneau, T. M.; dataset `lung`).
- **Licence:** LGPL (>= 2) (R package survival). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survival
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/lung.csv` from the installed R package survival (`get_r("lung", "survival")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `06d7fd271a4917256cf8588751bc24bc09320af8286f3bd10816a00c7456a671`

### `dat_bcg.csv`

- **What:** BCG vaccine trials: 13 trials (tpos, tneg, cpos, cneg, ablat, year, alloc)
- **Origin and citation:** Colditz, G. A. et al. (1994). Efficacy of BCG vaccine in the prevention of tuberculosis. JAMA 271, 698-702; analysed in Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software 36(3), doi:10.18637/jss.v036.i03. Copy from the R package metadat (`dat.bcg`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_bcg.csv` from the installed R package metadat (`get_r("dat.bcg", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `d04ccd6b85c8203d99606f95fd9a3e7e1f536ff315d918a5f79357d7bd6336c0`

### `dat_hackshaw1998.csv`

- **What:** Passive smoking and lung cancer in women who never smoked: 37 studies (odds ratios with intervals)
- **Origin and citation:** Hackshaw, A. K., Law, M. R. & Wald, N. J. (1997). The accumulated evidence on lung cancer and environmental tobacco smoke. BMJ 315, 980-988. Copy from metadat (`dat.hackshaw1998`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_hackshaw1998.csv` from the installed R package metadat (`get_r("dat.hackshaw1998", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `eb846c09778dac91034481886064a1e448cdf953cbee0de01175fa834079c66f`

### `dat_normand1999.csv`

- **What:** Length of hospital stay in stroke care: 9 studies (means, SDs, sizes)
- **Origin and citation:** Normand, S.-L. T. (1999). Meta-analysis: formulating, evaluating, combining, and reporting. Statistics in Medicine 18, 321-359. Copy from metadat (`dat.normand1999`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_normand1999.csv` from the installed R package metadat (`get_r("dat.normand1999", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `17c7af5b2650b6c5a46d8ba9485bb5ed70e357e1843237029d9cd0a55b2a3f1f`

### `dat_curtis1998.csv`

- **What:** Plant mass under elevated CO2: 102 observations (means, SDs, sizes)
- **Origin and citation:** Curtis, P. S. & Wang, X. (1998). A meta-analysis of elevated CO2 effects on woody plant mass, form, and physiology. Oecologia 113, 299-313. Copy from metadat (`dat.curtis1998`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_curtis1998.csv` from the installed R package metadat (`get_r("dat.curtis1998", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `c664ec25fc10b74f6f76b2877b5fae2d02eafd420ad4dffa821b6744b8c54baa`

### `dat_hart1999.csv`

- **What:** Warfarin and stroke in atrial fibrillation: 6 trials (events and person-time)
- **Origin and citation:** Hart, R. G., Benavente, O., McBride, R. & Pearce, L. A. (1999). Antithrombotic therapy to prevent stroke in patients with atrial fibrillation. Annals of Internal Medicine 131, 492-501. Copy from metadat (`dat.hart1999`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_hart1999.csv` from the installed R package metadat (`get_r("dat.hart1999", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `6f1df3f6a90fc84b5f96bfc248929abbb69340ea1d95bee3fd8aaf57e792832e`

### `dat_mcdaniel1994.csv`

- **What:** Validity of employment interviews: 160 correlations with sample sizes
- **Origin and citation:** McDaniel, M. A. et al. (1994). The validity of employment interviews. Journal of Applied Psychology 79, 599-616. Copy from metadat (`dat.mcdaniel1994`).
- **Licence:** Published study results (facts); the metadat copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=metadat
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dat_mcdaniel1994.csv` from the installed R package metadat (`get_r("dat.mcdaniel1994", "metadat")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `fa2c67428f72ad6b28977e998ee009e40237c2d55394468dd7d5c4308a1cd2a0`

### `veteran.csv`

- **What:** Veterans' Administration lung cancer trial (137 patients)
- **Origin and citation:** Kalbfleisch, J. D. & Prentice, R. L. (1980). The Statistical Analysis of Failure Time Data. Wiley. Copy from the R package survival (`veteran`).
- **Licence:** The survival copy (Therneau) is LGPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survival
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/veteran.csv` from the installed R package survival (`get_r("veteran", "survival")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `f92e886cb1dce8ef80bcdea3bac19f8f96c681a6601a85a9a34117ce8dfeefd7`

### `mgus2.csv`

- **What:** Monoclonal gammopathy (MGUS): 1384 patients, progression to plasma cell malignancy and death
- **Origin and citation:** Kyle, R. A. et al. (2002). A long-term study of prognosis in monoclonal gammopathy of undetermined significance. New England Journal of Medicine 346, 564-569. Copy from survival (`mgus2`).
- **Licence:** The survival copy (Therneau) is LGPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survival
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/mgus2.csv` from the installed R package survival (`get_r("mgus2", "survival")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a556280f9c87a10d7145d34973461c18f2adf76c5059b9df94c32081ed1cface`

### `heart.csv`

- **What:** Stanford heart transplant data in start-stop form
- **Origin and citation:** Crowley, J. & Hu, M. (1977). Covariance analysis of heart transplant survival data. JASA 72, 27-36. Copy from survival (`heart`).
- **Licence:** The survival copy (Therneau) is LGPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survival
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/heart.csv` from the installed R package survival (`get_r("heart", "survival")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `e37636a3184e98e73580707a5d147213770f1174aaa05ed839d999cbb08b7a5d`

### `udca.csv`

- **What:** Ursodeoxycholic acid trial in primary biliary cirrhosis (dates of entry and of events)
- **Origin and citation:** Lindor, K. D. et al. (1994). Ursodeoxycholic acid in the treatment of primary biliary cirrhosis. Gastroenterology 106, 1284-1290. Copy from survival (`udca`).
- **Licence:** The survival copy (Therneau) is LGPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survival
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/udca.csv` from the installed R package survival (`get_r("udca", "survival")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a1896edfc8c7cd145ccef092a9b78476fe95f40db1849fd55af43837300a8044`

### `bfeed.csv`

- **What:** Duration of breast feeding (927 mothers)
- **Origin and citation:** Klein, J. P. & Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and Truncated Data, 2nd ed. Springer. Copy from the R package KMsurv (`bfeed`).
- **Licence:** The KMsurv copy is GPL (>= 3). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=KMsurv
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/bfeed.csv` from the installed R package KMsurv (`get_r("bfeed", "KMsurv")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `de96962ece40997f12e3b22fb6661751f9b1cb23c8dabf0801ae3aa1ec4b2128`

### `channing.csv`

- **What:** Channing House retirement community: age at entry and at death or leaving (462 residents)
- **Origin and citation:** Klein, J. P. & Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and Truncated Data, 2nd ed. Springer. Copy from the R package KMsurv (`channing`).
- **Licence:** The KMsurv copy is GPL (>= 3). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=KMsurv
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/channing.csv` from the installed R package KMsurv (`get_r("channing", "KMsurv")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `1bf8ad70c8fb4c581f1ef4f8ce7fc2bdca1bb00532e5d278c1d5c6f218b5757f`

### `larynx.csv`

- **What:** Larynx cancer: survival of 90 men by stage
- **Origin and citation:** Klein, J. P. & Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and Truncated Data, 2nd ed. Springer. Copy from the R package KMsurv (`larynx`).
- **Licence:** The KMsurv copy is GPL (>= 3). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=KMsurv
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/larynx.csv` from the installed R package KMsurv (`get_r("larynx", "KMsurv")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `aff49fb728462f68a4749d197fc47d5faf21363eaeda47162c42196ed6236f44`

### `bci.csv`

- **What:** Barro Colorado Island trees: counts of 225 species in 50 one-hectare plots
- **Origin and citation:** Condit, R. et al. (2002). Beta-diversity in tropical forest trees. Science 295, 666-669. Copy from the R package vegan (`BCI`).
- **Licence:** The vegan copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=vegan
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/bci.csv` from the installed R package vegan (`get_r("BCI", "vegan")`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `022a5f8ef1b6dcbdac1130e5298b93d7b945bfa7ebfdd948026b5cbc81eb40bb`

### `dune.csv`

- **What:** Dutch dune meadows: 30 species in 20 sites, with A1, Moisture, Management, Use and Manure (from `dune.env`)
- **Origin and citation:** Jongman, R. H. G., ter Braak, C. J. F. & van Tongeren, O. F. R. (1987). Data Analysis in Community and Landscape Ecology. Pudoc, Wageningen. Copy from the R package vegan (`dune`, `dune.env`).
- **Licence:** The vegan copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=vegan
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/dune.csv` from the installed R package vegan (`local({ d <- get_r("dune", "vegan"); cbind(plot = as.integer(rownames(d)), d, get_r("dune.env", "vegan")) })`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `e3ea44df26e10c60c3ce03f295da4c647287a1784ac94b2c94d30afad11d81be`

### `varespec.csv`

- **What:** Lichen pastures: cover of 44 species at 24 sites, with 14 soil variables (from `varechem`)
- **Origin and citation:** Väre, H., Ohtonen, R. & Oksanen, J. (1995). Effects of reindeer grazing on understorey vegetation in dry Pinus sylvestris forests. Journal of Vegetation Science 6, 523-530. Copy from the R package vegan (`varespec`, `varechem`).
- **Licence:** The vegan copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=vegan
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/varespec.csv` from the installed R package vegan (`local({ d <- get_r("varespec", "vegan"); cbind(site = as.integer(rownames(d)), d, get_r("varechem", "vegan")) })`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `5d1265ae93aad0d558486fd8952d5e17afe17f484ddf9c9c85dc311c527174c6`

### `mite.csv`

- **What:** Oribatid mites: 35 species in 70 soil cores, with Topo (from `mite.env`)
- **Origin and citation:** Borcard, D. & Legendre, P. (1994). Environmental control and spatial structure in ecological communities. Environmental and Ecological Statistics 1, 37-61. Copy from the R package vegan (`mite`, `mite.env`).
- **Licence:** The vegan copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=vegan
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/mite.csv` from the installed R package vegan (`cbind(get_r("mite", "vegan"), Topo = get_r("mite.env", "vegan")$Topo)`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a29a64f042d04c8ecdbb30c8e4bb7c3b46d390225147dbe55ec65a7f211249f9`

### `holzinger_swineford_1939.csv`

- **What:** Mental ability test scores of 301 pupils in two schools (id, sex, age, school, grade, x1-x9)
- **Origin and citation:** Holzinger, K. J. & Swineford, F. (1939). A study in factor analysis: the stability of a bi-factor solution. Supplementary Educational Monographs 48, University of Chicago. Nine-test version and model in Rosseel, Y. (2012). lavaan: an R package for structural equation modeling. Journal of Statistical Software 48(2), 1-36. Copy from the R package lavaan (`HolzingerSwineford1939`).
- **Licence:** Published scores (facts); the lavaan copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://doi.org/10.18637/jss.v048.i02
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/holzinger_swineford_1939.csv` from the installed R package lavaan (`local({ d <- get_r("HolzingerSwineford1939", "lavaan"); d$grade <- as.numeric(d$grade); d$grade[is.na(d$grade)] <- -2147483648; d })`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below). One cell is written as it was in the stored copy: the missing grade of pupil 351 (NA in lavaan) is -2147483648, R's internal code for a missing whole number, which an earlier copy let through. No notebook uses the column grade.
- **sha256:** `fc0cf11342dd793f9f4781c4480efa963bb86d735f8cf7fc1275d53d798a073c`

### `turtles.csv`

- **What:** Homing directions of 10 green sea turtles (degrees)
- **Origin and citation:** Luschi, P. et al. (2001). Behavioral Ecology and Sociobiology 50, 528-534. Copy from the R package circular (`turtles`).
- **Licence:** The circular copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=circular
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/turtles.csv` from the installed R package circular (`data.frame(arrival = get_r("turtles", "circular")$arrival)`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `66d940163cee20390d54b745aa1266118c9ecd1b5d4d14abfbd4b3c28191422b`

### `icu_arrivals.csv`

- **What:** Arrival times of 254 patients at an intensive care unit (hours of the day)
- **Origin and citation:** Cox, D. R. & Lewis, P. A. W. (1966). The Statistical Analysis of Series of Events. Methuen; Fisher, N. I. (1993). Statistical Analysis of Circular Data, data set B.1. Copy from circular (`fisherB1`).
- **Licence:** The circular copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=circular
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/icu_arrivals.csv` from the installed R package circular (`data.frame(hour = as.numeric(get_r("fisherB1", "circular")))`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `39bfeb683d0c15a50eb3b0bba81355189655708816562b16a2889f3daf55b4a9`

### `wind_directions.csv`

- **What:** 310 wind directions (radians)
- **Origin and citation:** Copy from the R package circular (`wind`).
- **Licence:** The circular copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=circular
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/wind_directions.csv` from the installed R package circular (`data.frame(direction_radians = as.numeric(get_r("wind", "circular")))`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `3da9d2081e6b56c4dc2e2cbce44e8a528f213e9c4bb4472cdcd0394b21166844`

### `fisher_b10.csv`

- **What:** Three samples of directions (degrees)
- **Origin and citation:** Fisher, N. I. (1993). Statistical Analysis of Circular Data. Cambridge University Press, data set B.10. Copy from circular (`fisherB10`).
- **Licence:** The circular copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=circular
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/fisher_b10.csv` from the installed R package circular (`local({ b <- get_r("fisherB10", "circular"); do.call(rbind, lapply(names(b), function(s) data.frame(set = s, direction = as.numeric(b[[s]])))) })`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a513757aabf15ee43e3ed671d19b9bb94b3324d467334c5b30c989804155a34a`

### `meuse.csv`

- **What:** Heavy metals in the topsoil of the Meuse flood plain (155 samples)
- **Origin and citation:** Rikken, M. G. J. & Van Rijn, R. P. G. (1993). Soil pollution with heavy metals. Utrecht University. Copy from the R package sp (`meuse`).
- **Licence:** The sp copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=sp
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/meuse.csv` from the installed R package sp (`get_r("meuse", "sp")[, c("x", "y", "cadmium", "copper", "lead", "zinc")]`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `f5d5aea11a05e6f24350cebcb8316caee2866f462e3e305168bdcd08e5128e0b`

### `portpirie.csv`

- **What:** Annual maximum sea levels at Port Pirie, South Australia
- **Origin and citation:** Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values. Springer, section 3.4.1. Copy from the R package ismev (`portpirie`).
- **Licence:** Published data (facts); the ismev copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=ismev
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/portpirie.csv` from the installed R package ismev (`local({ d <- get_r("portpirie", "ismev"); data.frame(year = d$Year, sea_level = d$SeaLevel) })`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `2facadf76c82afcd17af026e474ab11fb40114e4e6d39d9b5bad10f83bef7a23`

### `coles_rain.csv`

- **What:** Daily rainfall at a location in south-west England, 1914-1962
- **Origin and citation:** Coles (2001), section 4.4.1. Copy from ismev (`rain`).
- **Licence:** Published data (facts); the ismev copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=ismev
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/coles_rain.csv` from the installed R package ismev (`data.frame(rain = as.numeric(get_r("rain", "ismev")))`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `0184509c9c36f013d56ef668289e58cd761563b5d19428afb411ec0a8d96fff2`

### `fremantle.csv`

- **What:** Annual maximum sea levels at Fremantle, Western Australia
- **Origin and citation:** Coles (2001), section 6.3.1. Copy from ismev (`fremantle`).
- **Licence:** Published data (facts); the ismev copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=ismev
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/fremantle.csv` from the installed R package ismev (`local({ d <- get_r("fremantle", "ismev"); data.frame(year = d$Year, sea_level = d$SeaLevel, soi = d$SOI) })`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `18c1283d45ea3272b538bb9b5e12746d3d97f5f70456cc2274380c84889f57d1`

### `hedenfalk_p.csv`

- **What:** 3170 p-values from the Hedenfalk breast cancer gene expression study
- **Origin and citation:** Hedenfalk, I. et al. (2001). New England Journal of Medicine 344, 539-548; Storey, J. D. & Tibshirani, R. (2003). PNAS 100, 9440-9445. Copy from the Bioconductor package qvalue (`hedenfalk$p`).
- **Licence:** The qvalue copy is LGPL. Not stored here (see the next lines).
- **URL:** https://bioconductor.org/packages/qvalue
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/hedenfalk_p.csv` from the installed R package qvalue (`data.frame(p = get_r("hedenfalk", "qvalue")$p)`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a11590e1943a9a070799bf1a93eb3b8838e2200043c643e0baf52bfce546af2e`

### `ryegrass.csv`

- **What:** Root length of ryegrass against ferulic acid concentration
- **Origin and citation:** Inderjit, Streibig, J. C. & Olofsdotter, M. (2002). Weed Biology and Management 2, 89-95. Copy from the R package drc (`ryegrass`).
- **Licence:** The drc copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=drc
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/ryegrass.csv` from the installed R package drc (`local({ d <- get_r("ryegrass", "drc"); data.frame(conc = d$conc, root_length = d$rootl) })`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `649ffcd565321b4018e0dd3ee226e1fa6b66548e1d4bccebe40aa42ba523a34d`

### `leaflife.csv`

- **What:** Leaf longevity and leaf mass per area of plant species
- **Origin and citation:** Wright, I. J. et al.; analysed in Warton, D. I. et al. (2006). Biological Reviews 81, 259-291. Copy from the R package smatr (`leaflife`).
- **Licence:** The smatr copy is GPL (>= 3). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=smatr
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/leaflife.csv` from the installed R package smatr (`get_r("leaflife", "smatr")[, c("site", "rain", "soilp", "longev", "lma")]`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `d683106dbfa8b97ad16f4c76a548067751f7870a3393e2c6b9e0485d62cdcf5e`

### `creatinine.csv`

- **What:** Creatinine measured by two methods (108 complete pairs)
- **Origin and citation:** Copy from the R package mcr (`creatinine`, rows with empty cells removed).
- **Licence:** The mcr copy is GPL (>= 3). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=mcr
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/creatinine.csv` from the installed R package mcr (`local({ d <- get_r("creatinine", "mcr"); d <- d[complete.cases(d), ]; data.frame(serum = d$serum.crea, plasma = d$plasma.crea) })`; written with R's `write.csv(row.names = FALSE, na = "")`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `f4cd7ccc1ad534d86723a79a0f6fa02646deeeb237fee737e863b931ea0c3a9d`

### `apistrat.csv`

- **What:** Stratified sample of 200 California schools (cds, stype, dnum, snum, pw, fpc, api00, api99, enroll, meals, ell, sch.wide, awards)
- **Origin and citation:** California Department of Education, Academic Performance Index 2000. Copy from the R package survey (`apistrat`); Lumley, T. (2010). Complex Surveys: A Guide to Analysis Using R. Wiley.
- **Licence:** Public data of the California Department of Education; the survey copy is GPL-2 | GPL-3. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survey
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/apistrat.csv` from the installed R package survey (`get_r("apistrat", "survey", "api")[, c("cds", "stype", "dnum", "snum", "pw", "fpc", "api00", "api99", "enroll", "meals", "ell", "sch.wide", "awards")]`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `3fb5d238813b6c11bf260ba12b3a37bae10c828c8db8af4105928cd8aee00f39`

### `apiclus1.csv`

- **What:** One-stage cluster sample of 183 schools in 15 districts (same columns as apistrat)
- **Origin and citation:** As `apistrat.csv` (R survey `apiclus1`; Lumley 2010).
- **Licence:** As `apistrat.csv`. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=survey
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/apiclus1.csv` from the installed R package survey (`get_r("apiclus1", "survey", "api")[, c("cds", "stype", "dnum", "snum", "pw", "fpc", "api00", "api99", "enroll", "meals", "ell", "sch.wide", "awards")]`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a4e7401c5d1dbf3a094cf4274460912f0f996237218a56107672fd433dc512cf`

### `lalonde.csv`

- **What:** 614 men from the NSW job-training evaluation and comparison groups (treat, age, educ, race, married, nodegree, re74, re75, re78)
- **Origin and citation:** LaLonde, R. J. (1986). Evaluating the econometric evaluations of training programs with experimental data. American Economic Review 76, 604-620; Dehejia, R. H. & Wahba, S. (1999). JASA 94, 1053-1062. Copy from the R package MatchIt (`lalonde`).
- **Licence:** No licence is stated for the data themselves; the MatchIt package is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=MatchIt
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/lalonde.csv` from the installed R package MatchIt (`get_r("lalonde", "MatchIt")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `a297d41d26ce7570b38bdab17d7a2330dc0ea4de4b080ce2c716666a5071f8cf`

### `jobs_ii.csv`

- **What:** 899 participants of the JOBS II job-search intervention (all 17 columns)
- **Origin and citation:** Vinokur, A. D., Price, R. H. & Schul, Y. (1995). Impact of the JOBS intervention on unemployed workers varying in risk for depression. American Journal of Community Psychology 23, 39-74. Copy from the R package mediation (`jobs`); Tingley, D. et al. (2014). Journal of Statistical Software 59(5).
- **Licence:** No licence is stated for the data themselves; the mediation package is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=mediation
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/jobs_ii.csv` from the installed R package mediation (`get_r("jobs", "mediation")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `c59fa202f2134796a366b72cc62c14b4d376357497c5f3159e4710ef9f24623b`

### `bioequivalence_reference.csv`

- **What:** The three 2 x 2 crossover reference datasets A, B and C
- **Origin and citation:** Schütz, H., Labes, D. & Fuglsang, A. (2014). Reference datasets for 2-treatment, 2-sequence, 2-period bioequivalence studies. The AAPS Journal 16, 1292-1297 (datasets A: Sauter et al. 1992; B: Clayton & Leslie 1981; C: B with subjects removed).
- **Licence:** Measured values from published studies (facts); the article itself is (c) AAPS 2014 and is not copied. Not stored here (see the next lines).
- **URL:** https://doi.org/10.1208/s12248-014-9661-0
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/bioequivalence_reference.csv` from the publisher's electronic supplementary material 1, 2 and 3 of the article (datasets A, B and C; https://static-content.springer.com/esm/art%3A10.1208%2Fs12248-014-9661-0/MediaObjects/12248_2014_9661_MOESM1_ESM.txt and MOESM2, MOESM3): subject, period, formulation (R = reference, T = test) and AUC, rows by dataset, subject and period. The article is also on PubMed Central (PMC4389751), but PMC hands its supplementary files only to a web browser, so the same files are taken from the publisher. The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `f7a5a9477ebc4894168933d515f0470ae6e98fc00b6e015043ae54b709f78b98`

### `sleepstudy.csv`

- **What:** Reaction times of 18 sleep-deprived truck drivers over 10 days (180 rows: Reaction, Days, Subject)
- **Origin and citation:** Belenky, G. et al. (2003). Patterns of performance degradation and restoration during sleep restriction and subsequent recovery: a sleep dose-response study. Journal of Sleep Research 12, 1-12. Copy from the R package lme4 (`sleepstudy`); Bates, D., Mächler, M., Bolker, B. & Walker, S. (2015). Journal of Statistical Software 67(1).
- **Licence:** The lme4 copy is GPL (>= 2). Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=lme4
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/sleepstudy.csv` from the installed R package lme4 (`get_r("sleepstudy", "lme4")`; written with R's `write.csv(row.names = FALSE, na = "")`, Windows (CRLF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `46683cc052263b3077cc43af9e5ffe75b46581be098e30f64a15b26d1065783f`

### `race_dolls.csv`

- **What:** Doll choice of 71 white and 89 black children (160 rows: child, doll), the 2 x 2 table 50/21 and 48/41 written out one child per row
- **Origin and citation:** Hraba, J. & Grant, G. (1970). Black is beautiful: a reexamination of racial preference and identification. Journal of Personality and Social Psychology 16(3), 398-402. The table as in the R package BayesFactor (`raceDolls`).
- **Licence:** The published counts (facts from a 1970 paper); BayesFactor's copy is GPL-2. Not stored here (see the next lines).
- **URL:** https://cran.r-project.org/package=BayesFactor
- **Fetched, not stored:** `scripts/get_data.py` writes `data/fetched/race_dolls.csv` from the installed R package BayesFactor (`local({ m <- get_r("raceDolls", "BayesFactor"); do.call(rbind, lapply(colnames(m), function(child) do.call(rbind, lapply(rownames(m), function(doll) data.frame(child = rep(tolower(child), m[doll, child]), doll = rep(tolower(doll), m[doll, child])))))) })`; written with R's `write.csv(row.names = FALSE, na = "", quote = FALSE)`, Unix (LF) line endings). The copy stored in this repository until 2026-11-22 (it stays in the git history) was removed after the fetched file was checked to be the same byte for byte (the sha256 below).
- **sha256:** `3fa017744473a9146f81a4264a5d92a96fe6923c18e6ebb0c8e9a82fcb792bea`

## NIST Statistical Reference Datasets (StRD)

- **Origin and citation:** National Institute of Standards and Technology, Statistical Reference Datasets,
  https://www.itl.nist.gov/div898/strd/ . Cite each by name as NIST StRD, with the original source NIST
  gives in the file header (for example Longley, J. W. (1967), Journal of the American Statistical
  Association 62, 819-841, for Longley; Dorsey, E. N. (1944) for Michelso). The certified values are in
  each file's header.
- **Licence:** public domain in the United States (works of the US Government, 17 U.S.C. 105); stored
  unchanged.
- **Retrieved:** 2026-11-20, downloaded with `scripts/get_data.py` from the URLs below.

### Univariate summary statistics

| file | URL | sha256 |
|---|---|---|
| `nist/PiDigits.dat` | https://www.itl.nist.gov/div898/strd/univ/data/PiDigits.dat | `e9421a2f4c0a31c505099332b4d22e507e808953c3c57bddae9226f25587d7d2` |
| `nist/Lottery.dat` | https://www.itl.nist.gov/div898/strd/univ/data/Lottery.dat | `956cb33392ee28683aa7978326fc356d1461e4a3b0e5a69d317fd481d02d6255` |
| `nist/Lew.dat` | https://www.itl.nist.gov/div898/strd/univ/data/Lew.dat | `081b231f686d01fa00247b8072909206947c8321bd62e060d9aaae5f4866c44e` |
| `nist/Mavro.dat` | https://www.itl.nist.gov/div898/strd/univ/data/Mavro.dat | `6b967538d6a607fdf30df054c6fbf4ca6bfb4fa6dc2ecbe613d4f93237a0b1dd` |
| `nist/Michelso.dat` | https://www.itl.nist.gov/div898/strd/univ/data/Michelso.dat | `68cad435641539ca020adafb4c9a6b7fcc4fefe1d6ab8f2af8fb858d4f1b12ab` |
| `nist/NumAcc1.dat` | https://www.itl.nist.gov/div898/strd/univ/data/NumAcc1.dat | `ce7c7a0fc6d7cedd0d5e1f9911afb9337544be6836730beb9968ca0a4571d375` |
| `nist/NumAcc2.dat` | https://www.itl.nist.gov/div898/strd/univ/data/NumAcc2.dat | `5be29522dfbe299541f0ec98ec1c22a51a5587bf06561ba21a92c209fb7ef0f9` |
| `nist/NumAcc3.dat` | https://www.itl.nist.gov/div898/strd/univ/data/NumAcc3.dat | `66df1b96596735936f37fd889e7c1c6ce998eb2b523bfc46c98ec1737e992aec` |
| `nist/NumAcc4.dat` | https://www.itl.nist.gov/div898/strd/univ/data/NumAcc4.dat | `ca310dc767f5130f980f8280bbe69e26ba414a4d85a89fc11b6c767b828e5fe4` |

### Linear regression

| file | URL | sha256 |
|---|---|---|
| `nist/Norris.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Norris.dat | `5ab6906c68a9c4b2eb1661c0cd7d2a09f5a469a6d4e05b6b691f190f5f7b944f` |
| `nist/Pontius.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Pontius.dat | `de60baa5dc66fd804f9316f895bc5aa2b8f791c09d1ecc1ec39f8b199b8c15c7` |
| `nist/NoInt1.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/NoInt1.dat | `d8234c428dedaa9de23b67df7d04aeee83befae73ceb6b38321b224955a9b667` |
| `nist/NoInt2.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/NoInt2.dat | `a71c2810b1edcef8fec6b864c20c5f9d9908e530282a174772ce610999ca7c69` |
| `nist/Filip.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Filip.dat | `403b34689e401d915cb28bd69311d3ca5c3007d88d3b670b00e3203042315072` |
| `nist/Longley.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Longley.dat | `fc4b0c824f8f0b834eb8e7d1c4c7e0b3cebdd7987bcb085e8f5e1007582480ea` |
| `nist/Wampler1.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler1.dat | `8bcd6f00dbe8c307ebc5c971e12b950eb3cdb75afbc2b38573082b5fcb0d7fbc` |
| `nist/Wampler2.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler2.dat | `01ba5f08b5b9f84431f3754216e006939a80264a1caa69fe08b3cbc543f9ee3e` |
| `nist/Wampler3.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler3.dat | `c5b025f6ebf0365de92d80e478ad826b01a6dfdafb71964ea8bae7eec993ce53` |
| `nist/Wampler4.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler4.dat | `dcff32075a52ecf31af70dad340776237544b1a29b5034a40f8b3e6b607cf9c1` |
| `nist/Wampler5.dat` | https://www.itl.nist.gov/div898/strd/lls/data/LINKS/DATA/Wampler5.dat | `c4e260a6638db61287bfea75a2dd3b324c86c5d40931ddc6d8fab0c0305ffa39` |

### Analysis of variance

| file | URL | sha256 |
|---|---|---|
| `nist/SiRstv.dat` | https://www.itl.nist.gov/div898/strd/anova/SiRstv.dat | `5a4f02d7f256e18c30dcad1b81e0260bb810b5ac113f6107e06dee168d07663e` |
| `nist/SmLs01.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs01.dat | `ad28cac0d4fb79ffea2522b55349083bdea10f12a6f8c46cfbabcccdb4d6b896` |
| `nist/SmLs02.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs02.dat | `708e32eaf36769abe7533e7dcd1b6dc67091b60df938dd37c4f0e65660ef7b0b` |
| `nist/SmLs03.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs03.dat | `e18404d0628b78670580b230a3a3b6cd997d959d00f92a2883a4da942b19377a` |
| `nist/SmLs04.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs04.dat | `1308eb1581f28ce6c9de4ca9520fcc729f916c34f1f9dee10fb6e6581e5b4ce7` |
| `nist/SmLs05.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs05.dat | `7b1333e6e56ced74145096df774cfb1f6305213cb28f9a8222c418ca85d0c296` |
| `nist/SmLs06.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs06.dat | `0c0f09a2df05f2e9db5324a6aaeed85647651e874e9f5c3892b8f29fc4a200f7` |
| `nist/SmLs07.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs07.dat | `fd95c510aa8b3ce7f67ddd32846ccead7d9abae441790336c6fb0fcb0ae3558f` |
| `nist/SmLs08.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs08.dat | `87dd38dfffeb65d39dbe3e8ba73bd16644fc396c1aad03fc0295225f3ba6aef9` |
| `nist/SmLs09.dat` | https://www.itl.nist.gov/div898/strd/anova/SmLs09.dat | `f4688c5a0d5690fc86fa27e0c5a265eef0b24641e74af27aa72d9615084dd13f` |
| `nist/AtmWtAg.dat` | https://www.itl.nist.gov/div898/strd/anova/AtmWtAg.dat | `5a1083c9062d3966f244bffa376fc9b033244abf86a0d0dd016af39cd71fb8a7` |

### Nonlinear regression

| file | URL | sha256 |
|---|---|---|
| `nist/Misra1a.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1a.dat | `9d3194d537d9a6bbf92eb819504070923a1885d0191f40b6b3aaeedd130b98e2` |
| `nist/Chwirut2.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Chwirut2.dat | `c37e5522a34f9a57e18d7d10f8736202eabc0c33b8e7ad0010e618db04b57411` |
| `nist/Chwirut1.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Chwirut1.dat | `d9a055dfe5af71a8754c00f073ef00f8fed2e3fd1c6fd20cea8fd62d7cc3ed84` |
| `nist/Lanczos3.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Lanczos3.dat | `bc998624541c49fedee08346ee2c1d495741c4ea1dd15c052d7c743ab3cf0efe` |
| `nist/Gauss1.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Gauss1.dat | `20d0b85ffc73681ce9bdd69c566e53969a65e0e2192c3de16ad4a1dcfdca6c27` |
| `nist/Gauss2.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Gauss2.dat | `24fe86801c35c4f86140cced68894e31d73b6ed5738cf698474144f94bae3640` |
| `nist/DanWood.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/DanWood.dat | `5eae83de8b34e34c631a5738d640183b0628913e728437ae77f066e832186e24` |
| `nist/Misra1b.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1b.dat | `dd27e310f787dcbc2addc50b9f8e5e3b82155f732926898b1858e4d08eead6a8` |
| `nist/Kirby2.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Kirby2.dat | `cbd288dd9d8586ec97ca408277d1a27ef41b19b543c158d594af189de80ae61c` |
| `nist/Hahn1.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Hahn1.dat | `7e885af48e18a222033a66cb22b491f0f5bc761b652142c0d9b28132daf63956` |
| `nist/Nelson.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Nelson.dat | `82fe77ef1a79a4441e7af5083680fc5046f9eb29b33bc2f2b76c478ecbe976d7` |
| `nist/MGH17.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH17.dat | `75edb5522ee639794ae1f0d0f4a1e3638e963c6854888b853e1f45cea3c2f4d7` |
| `nist/Lanczos1.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Lanczos1.dat | `dc7a419a4d7e7c0891e378dcde953a3e5dd273618b64b2cee8301e51577b428a` |
| `nist/Lanczos2.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Lanczos2.dat | `be1ac732e725ef75298f082bcdedab9e2ee9395ad77d8fad2ef37c95e0421572` |
| `nist/Gauss3.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Gauss3.dat | `24842a9e795db59c2c6cbc4c11047f90faedad014f8146faa6b5555259958d6d` |
| `nist/Misra1c.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1c.dat | `012f75aa0ff112afabdff73880e324ab28e2dc267b678537611c403b959dc941` |
| `nist/Misra1d.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Misra1d.dat | `4dfa39a381a5ffea3e5d9306af1bda9d53bc57226a8e6f0c19af6b9702580481` |
| `nist/Roszman1.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Roszman1.dat | `fb526231ae6b02929e3de54e072e259b3d87424089abd4b1015ffb946b1e0ab0` |
| `nist/ENSO.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/ENSO.dat | `0c644f45b0e586cb1bb0613bad68f63c704efadacd086e768468d33d83f3e5d9` |
| `nist/MGH09.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH09.dat | `6419da83fd9ab34615769413dc2a6d35b11372f512a9401628a1565067c51d10` |
| `nist/Thurber.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Thurber.dat | `806b8d00a5365eff907561db796fec3af8a26a5f353297d4311fd46ed5b0146f` |
| `nist/BoxBOD.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/BoxBOD.dat | `fc44fbed04e3be16df6498e722d4c8aedfd0af210110bc396d9f99acecd77e21` |
| `nist/Rat42.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Rat42.dat | `fc97dd1f94fbcb381495065a053e279f23588a9ceb243276f2a2c36c8225985b` |
| `nist/MGH10.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/MGH10.dat | `a0eb140580fdd6dcf6fba2cbbd81afa598c5226f1e6d702d1313c72a7a6e6fb4` |
| `nist/Eckerle4.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Eckerle4.dat | `8ec2133b29a0a7cecbfccd51f7c61ee432f73c460bf5ce7106d672642e57a61e` |
| `nist/Rat43.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Rat43.dat | `36f08908d0f09b1d075e3701ca567701535c58ac0d795ff68fb3c04e172c2c5d` |
| `nist/Bennett5.dat` | https://www.itl.nist.gov/div898/strd/nls/data/LINKS/DATA/Bennett5.dat | `b34d160605b959c301d13aac6627b693743d32ad98a6c4f830a3f728cdff0f52` |
