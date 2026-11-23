"""Tools for the reference files of notebooks 25 to 32.

The notebooks run simplydata only. What they compare with is written beforehand to references/<notebook>.csv, one row
per value, in this order of preference:

1. printed: a number printed in the documentation of the library simplydata calls (a vignette, a help-page example
   with its output, a docstring), certified (NIST StRD) or published (a paper or textbook table). It is typed in as
   printed (`printed_digits`), with the source and a link; the tolerance is half a unit of the last printed digit.
2. run: only where nothing is printed, the library's own help-page example is run once by the reference script of the
   notebook (scripts/references/<notebook>.py) and its value saved, with the link to the page whose example it
   reproduces and the package version.

Columns: those of scripts/references/common.py (stream A's definition, shared by every notebook with a reference
script): dataset, method, quantity, value, values, digits as printed, source, source title, url, kind (printed |
certified | run reference), tolerance, tolerance type, reason, versions. Notebooks 25 to 32 also use the tolerance
types "signs", "unit", "factors", "groups" and "text" for columns stored value by value ("name [1]", "name [2]" ...),
read by notebooks/reference_checks.py.

The reference script also saves the R datasets the notebook reads (references/data/<notebook>/*.csv; made at run
time, not kept in the repository, so GPL data are never copied into it).
"""

from __future__ import annotations

import sys
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "notebooks"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from checks import load_csv, r_available, r_run, read_nist  # noqa: E402,F401  (shared helpers)

REFERENCES = ROOT / "references"
RMAN = "https://search.r-project.org/CRAN/refmans/{package}/html/{topic}.html"
VIGNETTE = "https://cran.r-project.org/web/packages/{package}/vignettes/{file}"
NIST_NLS = "https://www.itl.nist.gov/div898/strd/nls/data/{name}.shtml"
NIST_LLS = "https://www.itl.nist.gov/div898/strd/lls/data/{name}.shtml"


def half_unit(printed: str) -> float:
    text = printed.strip().lstrip("-")
    return 0.5 * 10 ** (-(len(text.split(".")[1]) if "." in text else 0))


def r_version(package: str) -> str:
    if not r_available():
        return ""
    import checks

    return str(checks._R.r(f'as.character(packageVersion("{package}"))')[0])


def py_version(package: str) -> str:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return ""


def rman(package: str, topic: str) -> str:
    return RMAN.format(package=package, topic=topic)


class References:
    """Collect the reference values and the data of one notebook, then save them."""

    def __init__(self, notebook: str):
        self.notebook = notebook
        self.rows: list[dict] = []
        self.data_dir = REFERENCES / "data" / notebook
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def _add(self, dataset, check, quantity, value, printed, kind, source, url, package, version, tolerance, how, term=""):
        """One row in the columns stream A defined (scripts/references/common.py): dataset, method, quantity, value,
        values, digits as printed, source, source title, url, kind, tolerance, tolerance type, reason, versions."""
        kind = {"published": "printed", "printed": "printed", "certified": "certified"}.get(kind, "run reference")
        title = source if kind != "run reference" else ""
        label = source if kind != "run reference" or not version else f"{source} ({package} {version})"
        self.rows.append({"dataset": dataset, "method": check, "quantity": quantity, "value": value, "values": "",
                          "digits as printed": printed, "source": label, "source title": title, "url": url, "kind": kind,
                          "tolerance": tolerance, "tolerance type": how, "reason": ""})

    def printed(self, dataset, check, quantity, printed: str, source, url, kind="printed", package="", version=""):
        """A value as printed in documentation, a paper or a table: compared to half a unit of its last digit."""
        self._add(dataset, check, quantity, float(printed), printed, kind, source, url, package, version, half_unit(printed),
                  "absolute")

    def printed_all(self, dataset, check, values: dict, source, url, kind="printed", package="", version=""):
        for quantity, printed in values.items():
            self.printed(dataset, check, quantity, printed, source, url, kind, package, version)

    def certified(self, dataset, check, quantity, value, digits, source, url):
        """A certified value (NIST StRD), compared by the number of significant digits that agree."""
        self._add(dataset, check, quantity, float(value), repr(value), "certified", source, url, "NIST StRD", "", digits, "lre")

    def run(self, dataset, check, quantity, value, source, url, package, tolerance, how="relative"):
        """A value from running the library's own help-page example once (only where nothing is printed)."""
        version = r_version(package) if source.startswith("R ") else py_version(package)
        self._add(dataset, check, quantity, float(value), "", "run", source, url, package, version, tolerance, how)

    def run_all(self, dataset, check, values: dict, source, url, package, tolerance, how="relative"):
        for quantity, value in values.items():
            self.run(dataset, check, quantity, value, source, url, package, tolerance, how)

    def run_vector(self, dataset, check, name, values, source, url, package, tolerance, how="relative", every=1):
        """A column of run values (quantity "name [1]", "name [2]" ...). how="signs": compared after giving the
        notebook's column the sign of this one; "unit": after scaling both to length 1 as well; "text": values that must
        be the same text; "groups": labels
        compared up to their numbering. every=100 keeps every 100th value (rows 1, 101, 201 ...) of a long column."""
        version = r_version(package) if source.startswith("R ") else py_version(package)
        values = list(values)
        for i, value in enumerate(values, start=1):
            if (i - 1) % every:
                continue
            self._add(dataset, check, f"{name} [{i}]", value if how in ("groups", "text") else float(value), "", "run", source, url,
                      package, version, tolerance, how)

    def run_matrix(self, dataset, check, name, matrix, source, url, package, tolerance, how="relative", every=1):
        matrix = np.atleast_2d(np.asarray(matrix, float))
        for j in range(matrix.shape[1]):
            self.run_vector(dataset, check, f"{name}, column {j + 1}", matrix[:, j], source, url, package, tolerance, how, every)

    def data(self, name: str, frame: pd.DataFrame) -> None:
        frame.to_csv(self.data_dir / f"{name}.csv", index=False)

    def save(self) -> pd.DataFrame:
        from common import versions_line

        table = pd.DataFrame(self.rows)
        table["versions"] = versions_line()
        REFERENCES.mkdir(exist_ok=True)
        table.to_csv(REFERENCES / f"{self.notebook}.csv", index=False, float_format="%.17g")
        kinds = table["kind"].value_counts().to_dict() if len(table) else {}
        print(f"{self.notebook}: {len(table)} reference values ({kinds}) written to references/{self.notebook}.csv")
        return table


def r_dataset(name: str, package: str = "datasets") -> pd.DataFrame:
    from r_datasets import r_dataset as read

    return read(name, package)


def r_series(name: str, package: str = "datasets") -> pd.DataFrame:
    from r_datasets import r_series as read

    return read(name, package)
