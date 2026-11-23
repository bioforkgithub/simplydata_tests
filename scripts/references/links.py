"""Links for the reference rows of notebooks 09 to 17 that common.url_of cannot name.

common.url_of links base R functions, SciPy functions, NIST and "R <package> <function>" sources. The reference
scripts of notebooks 09 to 17 also name R functions with "package::function", R functions of the survival package
by name alone, statsmodels, lifelines, pandas and NumPy, and formulas written out in the script itself. add_links()
fills the empty url of each such row: the version-pinned manual page (or the pinned manual's search page, where the
manual has no page per name), and for a formula written out, the reference script on GitHub. A row whose source is a
published paper or book keeps its citation in the source (the scripts give the DOI where they have one).
"""

from __future__ import annotations

import re
from importlib import metadata

import checks

R_FUNCTION_PACKAGE = {
    "survfit": "survival", "coxph": "survival", "survreg": "survival", "survdiff": "survival",
    "difftime": "base", "arima": "stats", "cov": "stats", "cmdscale": "stats", "factanal": "stats",
    "summary.manova": "stats", "pairwise.prop.test": "stats", "power.t.test": "stats", "power.prop.test": "stats",
    "poisson.test": "stats", "lm": "stats",
}
R_PACKAGE_ONLY = {"vegan", "circular", "survey", "psych", "qvalue", "agricolae", "blockrand", "epitools", "MatchIt",
                  "mediation", "rstanarm", "BayesFactor", "pwr", "lavaan", "extRemes", "gstat", "spdep", "ape"}
BASE_PACKAGES = {"base", "stats", "utils"}
REPOSITORY = "https://github.com/bioforkgithub/simplydata_tests/blob/main/scripts/references/{stem}.py"


def _version(package: str, python: bool = False) -> str:
    if python:
        try:
            return metadata.version(package)
        except metadata.PackageNotFoundError:
            return ""
    if not checks.r_available():
        return ""
    try:
        return checks.r_eval_text(f'as.character(packageVersion("{package}"))')
    except Exception:  # noqa: BLE001 - a package that is not installed gives no version
        return ""


def _r_page(package: str, function: str | None) -> str:
    if package in BASE_PACKAGES:
        return f"https://search.r-project.org/R/refmans/{package}/html/{function}.html"
    version = _version(package)
    tail = f" (version {version})" if version else ""
    if function:
        return f"https://search.r-project.org/CRAN/refmans/{package}/html/{function}.html{tail}"
    return f"https://cran.r-project.org/package={package}{tail}"


def link_for(source: str, stem: str) -> str:
    s = source.strip()
    m = re.match(r"R ([A-Za-z0-9.]+)::([A-Za-z0-9._]+)", s)
    if m:
        return _r_page(m.group(1), m.group(2))
    m = re.match(r"R ([A-Za-z0-9._]+)", s)
    if m:
        name = m.group(1)
        if name in R_FUNCTION_PACKAGE:
            return _r_page(R_FUNCTION_PACKAGE[name], name)
        if name in R_PACKAGE_ONLY:
            return _r_page(name, None)
    m = re.match(r"statsmodels ([A-Za-z_]+)", s)
    if m:
        return f"https://www.statsmodels.org/v{_version('statsmodels', True)}/search.html?q={m.group(1)}"
    m = re.match(r"lifelines ([A-Za-z_]+)", s)
    if m:
        return f"https://lifelines.readthedocs.io/en/v{_version('lifelines', True)}/search.html?q={m.group(1)}"
    m = re.match(r"pandas ([A-Za-z_.]+)", s)
    if m:
        return (f"https://pandas.pydata.org/pandas-docs/version/{_version('pandas', True)}/search.html"
                f"?q={m.group(1).split('(')[0]}")
    if s.startswith("SciPy"):
        return f"https://docs.scipy.org/doc/scipy-{_version('scipy', True)}/reference/stats.html"
    if re.search(r"\b(19|20)\d\d\b", s) and not re.search(r"formula|NumPy|SciPy", s):
        return ""   # a published paper or book: cited in the source itself
    # a formula or definition written out (NumPy, Python, "the formula", "definition of ..."): the script itself
    return REPOSITORY.format(stem=stem)


def add_links(recorder) -> None:
    """Fill the empty url of each reference row of a common.Recorder (before recorder.save())."""
    for row in recorder.rows:
        if not row.get("url"):
            row["url"] = link_for(row["source"], recorder.stem)
