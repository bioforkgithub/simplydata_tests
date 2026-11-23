"""Shared by the reference scripts in scripts/references/: record reference values, and the data the notebook reads.

A reference script computes, with the packages simplydata is compared against (SciPy, statsmodels, NumPy, R through
rpy2) or from certified and published values, every reference value of one notebook, and saves them to
references/<notebook stem>.csv: dataset, method, quantity, value (or values, for a whole curve), source,
tolerance, tolerance type, reason (for a reference that could not be made), versions. The data the notebook needs
(R datasets, docstring examples, derived tables) are written to references/data/<notebook stem>/<name>.csv, so
that the notebook itself reads everything with simplydata (sdata.table) and runs only simplydata words.

Columns of references/<notebook stem>.csv (shared by every notebook with a reference script):
dataset, method, quantity  - the key the notebook uses (checks.ReferenceChecks.check / check_values)
value                      - the reference number (empty when values or reason is filled)
values                     - a whole curve or list, as JSON (compared by its largest relative difference)
digits as printed          - the number exactly as the manual, paper or table prints it (printed references)
source                     - where the number comes from (R function, SciPy docstring, NIST, formula ...)
source title, url          - the title and a version-pinned link of the manual page or paper
kind                       - "printed", "certified" or "run reference"
tolerance, tolerance type  - "absolute", "relative", "lre" or "largest relative difference"
reason                     - filled when the reference could not be made, or the value cannot be read
versions                   - Python, the packages and R (with the R packages used) that made the file

Run from the repository root:  python scripts/references/<notebook stem>.py
"""

from __future__ import annotations

import json
import platform
import sys
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "notebooks"))   # the notebooks' helpers (checks.read_nist, realdata) for reading data

import checks  # noqa: E402

REFERENCES = ROOT / "references"


def versions_line(r_packages=()) -> str:
    parts = [f"Python {platform.python_version()}"]
    for name in ("numpy", "pandas", "scipy", "statsmodels", "rpy2"):
        try:
            parts.append(f"{name} {metadata.version(name)}")
        except metadata.PackageNotFoundError:
            pass
    if checks.r_available():
        parts.append(checks.r_eval_text("R.version.string"))
        for p in r_packages:
            try:
                parts.append(f"{p} {checks.r_eval_text(f'as.character(packageVersion(\"{p}\"))')}")
            except Exception:  # noqa: BLE001
                parts.append(f"{p}: not installed")
    return "; ".join(parts)


SCIPY = "1.18.1"
R_BASE = {"t.test", "wilcox.test", "kruskal.test", "oneway.test", "anova", "chisq.test", "ks.test", "shapiro.test",
          "binom.test", "fisher.test", "prop.test", "prop.trend.test", "poisson.test", "mantelhaen.test", "cor.test",
          "mcnemar.test", "ansari.test", "mood.test", "var.test", "bartlett.test", "fligner.test", "friedman.test",
          "quade.test", "p.adjust", "TukeyHSD", "lowess", "ecdf", "lm", "nls", "smooth.spline", "hclust", "quantile",
          "sd", "mean", "median", "sum", "max", "min", "table", "tapply", "aov", "pairwise.t.test",
          "pairwise.wilcox.test", "loess", "predict", "density", "dnorm", "pnorm", "weighted.mean", "optim"}


def kind_of(source: str) -> str:
    """"certified" (NIST), "printed" (a number printed in a manual, docstring, help page, paper or table) or
    "run reference" (the original library run by the reference script)."""
    low = source.lower()
    if "nist certified" in low:
        return "certified"
    if any(w in low for w in ("printed", "docstring", "help page text", "published")):
        return "printed"
    return "run reference"


def url_of(source: str) -> str:
    """A version-pinned link to the manual page behind a reference, where one can be named from its source."""
    import re

    if "NIST" in source:
        return "https://www.itl.nist.gov/div898/strd/"
    m = re.search(r"SciPy (\w+)", source)
    if m and m.group(1) not in ("quad", "bootstrap", "grid", "written", "formula"):
        name = m.group(1)
        return f"https://docs.scipy.org/doc/scipy-{SCIPY}/reference/generated/scipy.stats.{name}.html"
    if "SciPy bootstrap" in source:
        return f"https://docs.scipy.org/doc/scipy-{SCIPY}/reference/generated/scipy.stats.bootstrap.html"
    m = re.search(r"R ([A-Za-z0-9.]+)(?: ([A-Za-z0-9._]+))?", source)
    if m:
        first, second = m.group(1), m.group(2)
        if first in R_BASE:
            return f"https://search.r-project.org/R/refmans/stats/html/{first}.html"
        if second and checks.r_available():
            try:
                version = checks.r_eval_text(f'as.character(packageVersion("{first}"))')
                return f"https://search.r-project.org/CRAN/refmans/{first}/html/{second}.html (version {version})"
            except Exception:  # noqa: BLE001
                pass
    return ""


class Recorder:
    """The same calls as checks.Checks (compare, compare_all, skip), but it records the reference only: the value
    given as `ours` is not stored (the notebook computes simplydata's value itself, with simplydata words)."""

    def __init__(self, stem: str, r_packages=()):
        self.stem = stem
        self.rows: list[dict] = []
        self.r_packages = tuple(r_packages)
        (REFERENCES / "data" / stem).mkdir(parents=True, exist_ok=True)

    # ---- references
    def compare(self, dataset, method, quantity, ours, reference, reference_name, tolerance, how="absolute"):
        self.rows.append({"dataset": dataset, "method": method, "quantity": quantity, "value": float(reference),
                          "values": "", "digits as printed": "", "source": reference_name, "source title": "",
                          "url": url_of(reference_name), "kind": kind_of(reference_name), "tolerance": tolerance,
                          "tolerance type": how, "reason": ""})
        return True

    def printed(self, dataset, method, quantity, shown: str, reference_name, title="", url=""):
        """A number as a manual prints it (a string such as "0.0968"): compared to half a unit of its last digit."""
        decimals = len(shown.split(".")[1]) if "." in shown else 0
        self.rows.append({"dataset": dataset, "method": method, "quantity": quantity, "value": float(shown),
                          "values": "", "digits as printed": shown, "source": reference_name, "source title": title,
                          "url": url or url_of(reference_name), "kind": "printed",
                          "tolerance": 0.5 * 10 ** (-decimals) * (1 + 1e-9), "tolerance type": "absolute",
                          "reason": ""})

    def compare_all(self, dataset, method, ours, reference, reference_name, tolerance, how="absolute"):
        for quantity, value in reference.items():
            self.compare(dataset, method, quantity, None, value, reference_name, tolerance, how)

    def vector(self, dataset, method, quantity, reference, reference_name, tolerance):
        """A whole curve or list: the notebook compares the largest relative difference with 0."""
        values = [float(v) for v in np.asarray(reference, dtype=float).ravel()]
        self.rows.append({"dataset": dataset, "method": method, "quantity": quantity, "value": np.nan,
                          "values": json.dumps(values), "digits as printed": "", "source": reference_name,
                          "source title": "", "url": url_of(reference_name), "kind": kind_of(reference_name),
                          "tolerance": tolerance, "tolerance type": "largest relative difference", "reason": ""})

    def skip(self, dataset, method, quantity, reference_name, reason):
        self.rows.append({"dataset": dataset, "method": method, "quantity": quantity, "value": np.nan, "values": "",
                          "digits as printed": "", "source": reference_name, "source title": "",
                          "url": url_of(reference_name), "kind": kind_of(reference_name), "tolerance": np.nan,
                          "tolerance type": "", "reason": reason})

    # ---- data for the notebook
    def data(self, name: str, frame) -> Path:
        frame = pd.DataFrame(frame)
        path = REFERENCES / "data" / self.stem / f"{name}.csv"
        frame.to_csv(path, index=False)   # shortest round-trip form of each number (as Python repr)
        return path

    def save(self) -> pd.DataFrame:
        table = pd.DataFrame(self.rows)
        table["versions"] = versions_line(self.r_packages)
        key = ["dataset", "method", "quantity", "source"]
        duplicated = table[table.duplicated(key, keep=False)]
        if len(duplicated):
            print(f"note: {len(duplicated)} rows share a key; the notebook compares with each")
        table.to_csv(REFERENCES / f"{self.stem}.csv", index=False, float_format="%.17g")
        print(f"{self.stem}: {len(table)} reference values saved to references/{self.stem}.csv")
        return table
