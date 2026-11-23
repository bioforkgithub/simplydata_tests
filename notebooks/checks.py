"""Helpers shared by the notebooks: read the data, run R, compare two numbers and keep a record.

Every comparison is printed as one line and saved to results/<notebook>.csv, which 00_summary.ipynb collects.
Nothing here changes a number; it only reads, compares and records.

The words a reader sees (the column "result") say what each comparison found:
- "agrees": within the stated tolerance;
- "simplydata error – being fixed": a simplydata bug (FINDINGS.md class A), until its fix is released;
- "explained difference – reference": the reference rounds, stops early, uses another definition or has its own
  known bug (class C; the explanation names which);
- "explained difference – machine": the last digits differ between CPUs or systems;
- "explained difference – precision limit": the limit of double precision or of the package simplydata calls,
  the same on every system (class B);
- "difference not yet explained": outside the tolerance and not yet in results/known_fails.csv (CI fails on it);
- "waiting for word" and "skipped": no comparison could be made (the reason is given).
Every difference keeps its finding id, a one-line explanation and a link to FINDINGS.md (results/known_fails.csv).
The machine column "status" keeps PASS, FAIL or SKIP, so nothing is hidden or lost.
"""

from __future__ import annotations

import math
import os
import platform
import re
import sys
import warnings
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
REQUIRE_R = os.environ.get("SIMPLYDATA_TESTS_REQUIRE_R", "") == "1"

KEY = ["notebook", "dataset", "method", "quantity", "reference source"]
AGREES = "agrees"
SIMPLYDATA_ERROR = "simplydata error – being fixed"
REFERENCE = "explained difference – reference"
MACHINE = "explained difference – machine"
PRECISION = "explained difference – precision limit"
NOT_EXPLAINED = "difference not yet explained"
WAITING = "waiting for word"
SKIPPED = "skipped"
RESULT_ORDER = [AGREES, SIMPLYDATA_ERROR, REFERENCE, MACHINE, PRECISION, NOT_EXPLAINED, WAITING, SKIPPED]
_known: dict | None = None


def this_system() -> str:
    return "windows" if sys.platform == "win32" else "ubuntu"


def known_fails() -> dict:
    """The explained FAILs for this system, {key: row of results/known_fails.csv}."""
    global _known
    if _known is None:
        _known = {}
        path = RESULTS / "known_fails.csv"
        if path.exists():
            table = pd.read_csv(path, keep_default_na=False, dtype=str)
            for row in table.to_dict("records"):
                if row["platform"] in ("any", this_system()):
                    _known[tuple(row[c] for c in KEY)] = row
    return _known


def result_words(notebook, dataset, method, quantity, reference_name, status, reason="") -> dict:
    """The reader's words for one comparison: result, finding, explanation and link."""
    if status == "PASS":
        return {"result": AGREES, "finding": "", "explanation": "", "link": ""}
    if status == "SKIP":
        reason = reason if isinstance(reason, str) else ""
        waiting = reason.startswith("waiting for word")
        return {"result": WAITING if waiting else SKIPPED, "finding": "", "explanation": reason, "link": ""}
    row = known_fails().get((notebook, dataset, method, quantity, reference_name))
    if row is None:
        return {"result": NOT_EXPLAINED, "finding": "", "explanation": "", "link": ""}
    return {"result": row.get("result", ""), "finding": row["finding"], "explanation": row.get("explanation", ""),
            "link": row.get("link", "")}


def add_result_words(table: pd.DataFrame) -> pd.DataFrame:
    """Fill the columns result, finding, explanation and link of a results table from its status."""
    if not len(table):
        return table
    words = [result_words(r["notebook"], r["dataset"], r["method"], r["quantity"], r["reference source"],
                          r["status"], r.get("reason", "")) for r in table.to_dict("records")]
    for column in ("result", "finding", "explanation", "link"):
        table[column] = [w[column] for w in words]
    return table


def count_results(table: pd.DataFrame) -> str:
    counts = table["result"].value_counts().to_dict() if len(table) else {}
    return ", ".join(f"{counts[r]} {r}" for r in RESULT_ORDER if counts.get(r))


pd.set_option("display.width", 140)
pd.set_option("display.max_columns", 20)


# ---------------------------------------------------------------- versions
def versions(*packages: str) -> None:
    """Print the Python, package and R versions used, so every result can be reproduced."""
    print(f"Python {platform.python_version()} on {platform.system()} {platform.release()}")
    print(simplydata_version())
    for name in packages:
        try:
            print(f"{name} {metadata.version(name)}")
        except metadata.PackageNotFoundError:
            print(f"{name}: not installed")
    if r_available():
        print(r_eval_text("R.version.string"))


def simplydata_version() -> str:
    """simplydata's version and, when it is not a release, the commit it was taken from.

    The commit comes from pip's record of a `git+https://...@<commit>` install (direct_url.json), or from
    the environment variable SIMPLYDATA_TESTS_SIMPLYDATA_COMMIT when simplydata is read from a checkout."""
    import json

    import simplydata

    text = f"simplydata {simplydata.__version__}"
    commit = os.environ.get("SIMPLYDATA_TESTS_SIMPLYDATA_COMMIT", "")
    if not commit:
        try:
            record = metadata.distribution("simplydata").read_text("direct_url.json")
            commit = json.loads(record).get("vcs_info", {}).get("commit_id", "") if record else ""
        except Exception:  # noqa: BLE001 - no record means a release from PyPI
            commit = ""
    return f"{text} (commit {commit[:7]} of github.com/bioforkgithub/simplydata)" if commit else text


# ---------------------------------------------------------------- R through rpy2
_R = None


def r_available() -> bool:
    """True when R can be reached through rpy2. With SIMPLYDATA_TESTS_REQUIRE_R=1 a missing R is an error."""
    global _R
    if _R is None:
        try:
            import rpy2.robjects as ro

            ro.r("1 + 1")
            _R = ro
        except Exception as problem:  # noqa: BLE001 - any failure means R cannot be used here
            _R = False
            if REQUIRE_R:
                raise RuntimeError(f"R is required (SIMPLYDATA_TESTS_REQUIRE_R=1) but cannot be reached: {problem}")
    return bool(_R)


def r_eval_text(code: str) -> str:
    return str(_R.r(code)[0])


def r_run(code: str, data: pd.DataFrame | dict | None = None, packages: tuple = ()) -> dict:
    """Run R code that ends with a named list of numbers; columns of data come in as R vectors.

    Text columns become R character vectors, numbers R doubles (or integers), True/False R logicals.
    Gives back a dict of floats (or lists of floats for longer results).
    """
    ro = _R
    from rpy2.robjects import vectors

    for package in packages:
        ro.r(f"suppressPackageStartupMessages(library({package}))")
    env = ro.r("new.env()")
    if data is not None:
        columns = data.items() if isinstance(data, dict) else ((c, data[c]) for c in data.columns)
        for name, values in columns:
            values = pd.Series(values)
            if values.dtype == bool:
                env[name] = vectors.BoolVector(values.tolist())
            elif pd.api.types.is_numeric_dtype(values):
                env[name] = vectors.FloatVector(values.astype(float).tolist())
            else:
                env[name] = vectors.StrVector(values.astype(str).tolist())
    found = ro.r["eval"](ro.r["parse"](text=code), env)
    out = {}
    for name, value in zip(found.names, found):
        if isinstance(value, vectors.StrVector):
            out[name] = [str(v) for v in value]
            continue
        numbers = [float(v) for v in value]
        out[name] = numbers[0] if len(numbers) == 1 else numbers
    return out


# ---------------------------------------------------------------- NIST Statistical Reference Datasets
_NUMBER = r"(-?\d+\.?\d*(?:[Ee][-+]?\d+)?)"


@dataclass
class Nist:
    name: str
    data: pd.DataFrame
    certified: dict
    starts: dict          # nonlinear only: parameter -> (start 1, start 2)
    model: str            # the model line(s) of the file


def read_nist(name: str) -> Nist:
    """A NIST StRD file: the data, the certified values and (for nonlinear problems) the starting values."""
    text = (DATA / "nist" / f"{name}.dat").read_text(encoding="ascii").replace("\r", "")
    lines = text.split("\n")
    first, last = map(int, re.search(r"Data\s*[:(]?\s*\(?lines\s+(\d+)\s+to\s+(\d+)", text).groups())
    rows = [line.split() for line in lines[first - 1:last] if line.strip()]
    head = "\n".join(lines[:first - 1])
    certified, starts = {}, {}
    model = ""
    m = re.search(r"Model:\s+(.*?)\n\s*\n", head, flags=re.S)
    if m:
        model = " ".join(m.group(1).split())

    # univariate summary statistics
    for key, pattern in [("mean", r"Sample Mean\s+ybar:\s+"),
                         ("standard deviation", r"Sample Standard Deviation \(denom\. = n-1\)\s+s:\s+"),
                         ("autocorrelation", r"Sample Autocorrelation Coefficient \(lag 1\)\s+r\(1\):\s+")]:
        m = re.search(pattern + _NUMBER, head)
        if m:
            certified[key] = float(m.group(1))

    # linear regression: B0, B1, ... with their standard deviations
    for m in re.finditer(r"\n\s+B(\d+)\s+" + _NUMBER + r"\s+" + _NUMBER, head):
        certified[f"B{m.group(1)}"] = float(m.group(2))
        certified[f"B{m.group(1)} standard error"] = float(m.group(3))
    m = re.search(r"Residual\s+Standard Deviation\s+" + _NUMBER, head)
    if m:
        certified["residual standard deviation"] = float(m.group(1))
    m = re.search(r"\n\s+R-Squared\s+" + _NUMBER, head)
    if m:
        certified["R squared"] = float(m.group(1))

    # analysis of variance
    m = re.search(r"Between \w+\s+(\d+)\s+" + _NUMBER + r"\s+" + _NUMBER + r"\s+" + _NUMBER, head)
    if m:
        certified.update({"between df": float(m.group(1)), "between sum of squares": float(m.group(2)),
                          "between mean square": float(m.group(3)), "F": float(m.group(4))})
        m = re.search(r"Within \w+\s+(\d+)\s+" + _NUMBER + r"\s+" + _NUMBER, head)
        certified.update({"within df": float(m.group(1)), "within sum of squares": float(m.group(2)),
                          "within mean square": float(m.group(3))})
        certified["R squared"] = float(re.search(r"Certified R-Squared\s+" + _NUMBER, head).group(1))
        certified["residual standard deviation"] = float(
            re.search(r"Standard Deviation\s+" + _NUMBER, head).group(1))

    # nonlinear regression: b1 = start1 start2 certified sd
    for m in re.finditer(r"\n\s+(b\d+)\s+=\s+" + _NUMBER + r"\s+" + _NUMBER + r"\s+" + _NUMBER + r"\s+" + _NUMBER,
                         head):
        starts[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        certified[m.group(1)] = float(m.group(4))
        certified[f"{m.group(1)} standard error"] = float(m.group(5))
    m = re.search(r"Residual Sum of Squares:\s+" + _NUMBER, head)
    if m:
        certified["residual sum of squares"] = float(m.group(1))
        certified["residual standard deviation"] = float(
            re.search(r"Residual Standard Deviation:\s+" + _NUMBER, head).group(1))

    # the data: one response y and one or more x columns (ANOVA: the group first, then y)
    if "between df" in certified:
        data = pd.DataFrame({"group": [r[0] for r in rows], "y": [float(r[1]) for r in rows]})
    else:
        width = len(rows[0])
        names = ["y"] + (["x"] if width == 2 else [f"x{i}" for i in range(1, width)])
        data = pd.DataFrame([[float(v) for v in r] for r in rows], columns=names)
    return Nist(name, data, certified, starts, model)


def lre(value: float, certified: float) -> float:
    """The log relative error: how many significant digits of value agree with the certified value.

    -log10(|value - certified| / |certified|), or -log10(|value|) when the certified value is 0.
    Capped at 15 (the digits a double carries); 0 when no digit agrees.
    """
    if not np.isfinite(value):
        return 0.0
    if value == certified:
        return 15.0
    error = abs(value - certified) / abs(certified) if certified != 0 else abs(value)
    return float(min(15.0, max(0.0, -math.log10(error))))


# ---------------------------------------------------------------- comparing and recording
class Checks:
    """Compare numbers, print what each comparison found, and save the record to results/<notebook>.csv."""

    def __init__(self, notebook: str):
        self.notebook = notebook
        self.rows: list[dict] = []

    def compare(self, dataset: str, method: str, quantity: str, ours, reference, reference_name: str,
                tolerance: float, how: str = "absolute") -> bool:
        """how = "absolute" (|ours - reference| <= tolerance), "relative" (divided by |reference|),
        or "lre" (at least `tolerance` significant digits agree, for NIST certified values)."""
        ours, reference = float(ours), float(reference)
        difference = ours - reference
        if how == "lre":
            score = lre(ours, reference)
            passed = score >= tolerance
            shown = f"LRE {score:.1f} (needs >= {tolerance:g})"
        elif how == "relative":
            score = abs(difference) / abs(reference) if reference != 0 else abs(difference)
            passed = bool(score <= tolerance)
            shown = f"relative difference {score:.2e} (allowed {tolerance:.0e})"
        else:
            score = abs(difference)
            passed = bool(score <= tolerance)
            shown = f"difference {score:.2e} (allowed {tolerance:.0e})"
        if not (np.isfinite(ours) and np.isfinite(reference)):
            passed = ours == reference or (np.isnan(ours) and np.isnan(reference))
        status = "PASS" if passed else "FAIL"
        words = result_words(self.notebook, dataset, method, quantity, reference_name, status)
        label = words["result"] + (f" ({words['finding']})" if words["finding"] else "")
        print(f"{label}  {dataset} | {method} | {quantity}: simplydata {ours:.12g}, {reference_name} "
              f"{reference:.12g}; {shown}")
        if words["explanation"]:
            print(f"      why: {words['explanation']} - {words['link']}")
        self.rows.append({"notebook": self.notebook, "dataset": dataset, "method": method, "quantity": quantity,
                          "simplydata": ours, "reference": reference, "reference source": reference_name,
                          "difference": difference, "tolerance": tolerance, "tolerance type": how,
                          "score": score, "status": status, "reason": "", **words})
        return passed

    def compare_all(self, dataset: str, method: str, ours: dict, reference: dict, reference_name: str,
                    tolerance: float, how: str = "absolute") -> None:
        """Compare every quantity the reference has."""
        for quantity, value in reference.items():
            self.compare(dataset, method, quantity, ours[quantity], value, reference_name, tolerance, how)

    def skip(self, dataset: str, method: str, quantity: str, reference_name: str, reason: str) -> None:
        """Record a comparison that could not be made, with the reason (it is listed in the summary)."""
        words = result_words(self.notebook, dataset, method, quantity, reference_name, "SKIP", reason)
        print(f"{words['result']}  {dataset} | {method} | {quantity} against {reference_name}: {reason}")
        self.rows.append({"notebook": self.notebook, "dataset": dataset, "method": method, "quantity": quantity,
                          "simplydata": np.nan, "reference": np.nan, "reference source": reference_name,
                          "difference": np.nan, "tolerance": np.nan, "tolerance type": "", "score": np.nan,
                          "status": "SKIP", "reason": reason, **words})

    def save(self) -> pd.DataFrame:
        RESULTS.mkdir(exist_ok=True)
        table = pd.DataFrame(self.rows)
        table.to_csv(RESULTS / f"{self.notebook}.csv", index=False, float_format="%.17g")
        print(f"{self.notebook}: {len(table)} comparisons - {count_results(table)}")
        return table


def quiet(word, *args, **kwargs):
    """Run a simplydata word without printing its explanation (used inside loops; each notebook
    shows the full printed answer once, as an example)."""
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        return word(*args, **kwargs)


def dataset_path(name: str) -> Path:
    """The path of a dataset: data/<name> when it is stored in this repository, otherwise data/fetched/<name>
    (the datasets this repository does not host, written there by scripts/get_data.py)."""
    for path in (DATA / name, DATA / "fetched" / name):
        if path.exists():
            return path
    raise FileNotFoundError(f"data/fetched/{name} is missing: it is not stored in this repository. "
                            "Run `python scripts/get_data.py` first (it needs R; see README, 'Data').")


def load_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(dataset_path(name))


def _short_warning(message, category, filename, lineno, line=None):
    """Warnings are printed without the file they come from (that path is the reader's own computer)."""
    return f"{category.__name__}: {message}\n"


warnings.formatwarning = _short_warning

if sys.platform == "win32":  # print the minus and other signs on Windows consoles too
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass


# ---------------------------------------------------------------- reading a simplydata result with simplydata words
# The simplydata side of a comparison is read out with simplydata's own words only (get_row, get_column_names,
# count_rows, get_number_at, count); pandas and NumPy are used only on the reference side.
def read_row(result, row=1) -> dict:
    """One row of a simplydata result table as {column: value}: result.get_row(row), counting from 1."""
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        return result.get_row(row)


def read_rows(result) -> list[dict]:
    """Every row of a simplydata result table, each as {column: value} (table.count_rows and table.get_row)."""
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        return [result.get_row(i) for i in range(1, result.count_rows() + 1)]


def read_rows_by(result, column: str) -> dict:
    """The rows of a simplydata result table keyed by the value in one column: {value: {column: value}}."""
    return {row[column]: row for row in read_rows(result)}


def read_numbers(numbers) -> list:
    """The values of sdata.numbers, in order (nums.count and nums.get_number_at)."""
    return [numbers.get_number_at(i) for i in range(1, numbers.count() + 1)]


# ---------------------------------------------------------------- comparing simplydata with a reference file
# The notebooks 18 to 24 run only simplydata words. Their reference values are computed by
# scripts/references/<notebook>.py (SciPy, statsmodels, R, certified and published values) and saved to
# references/<notebook>.csv, which the notebook reads with simplydata (sdata.table). This class does the
# comparison; it is the test harness (it uses plain Python and checks.compare, never simplydata's own numbers
# from anywhere but the notebook).
REFERENCES = ROOT / "references"


def data_file(notebook: str, name: str) -> str:
    """The path of a data file that the notebook's reference script wrote (references/data/<notebook>/<name>.csv)."""
    return str(REFERENCES / "data" / notebook / f"{name}.csv")


class ReferenceChecks(Checks):
    """Checks against references/<notebook>.csv: check() and check_values() compare simplydata's value with every
    reference row of the same dataset, method and quantity; finish() lists the reference rows that no simplydata
    value met (as SKIP, with the reason) and saves the results."""

    def __init__(self, notebook: str):
        super().__init__(notebook)
        import simplydata as sdata

        with __import__("contextlib").redirect_stdout(__import__("io").StringIO()):
            table = sdata.table(str(REFERENCES / f"{notebook}.csv"), separator=",")
            self.references = read_rows(table)
        self.used: set[int] = set()

    def _matching(self, dataset, method, quantity):
        return [(i, r) for i, r in enumerate(self.references)
                if r["dataset"] == dataset and r["method"] == method and r["quantity"] == quantity]

    def check(self, dataset, method, quantity, ours) -> int:
        """Compare one simplydata value with every reference for (dataset, method, quantity)."""
        import json

        if ours is None:   # an empty cell in simplydata's answer
            ours = float("nan")
        elif isinstance(ours, (list, tuple)):
            ours = [float("nan") if v is None else v for v in ours]
        found = self._matching(dataset, method, quantity)
        if not found:
            raise KeyError(f"no reference for {dataset!r} | {method!r} | {quantity!r}")
        for i, ref in found:
            self.used.add(i)
            if ref.get("reason"):
                self.skip(dataset, method, quantity, ref["source"], ref["reason"])
            elif ref.get("values"):
                reference = np.asarray(json.loads(ref["values"]), dtype=float)
                values = np.atleast_1d(np.asarray(ours, dtype=float))
                if values.shape != reference.shape:
                    # a different number of values is a difference, recorded as such (not an error that stops the notebook)
                    self.compare(dataset, method, f"{quantity}: number of values", values.size, reference.size,
                                 ref["source"], 0.0)
                else:
                    score = float(np.max(np.abs(values - reference) / np.maximum(np.abs(reference), 1e-300)))
                    self.compare(dataset, method,
                                 f"{quantity}: largest relative difference over {len(reference)} values",
                                 score, 0.0, ref["source"], float(ref["tolerance"]))
            else:
                reference = float("nan") if ref["value"] is None else ref["value"]
                self.compare(dataset, method, quantity, ours, reference, ref["source"], float(ref["tolerance"]),
                             ref["tolerance type"])
            where = ref.get("url") or ref.get("source title") or ""
            print(f"      Source: {ref['source']} ({ref.get('kind') or ''}){' - ' + where if where else ''}")
            if self.rows:
                self.rows[-1]["kind"] = ref.get("kind") or ""
                self.rows[-1]["url"] = where
        return len(found)

    def check_values(self, dataset, method, values: dict) -> int:
        """Compare each of simplydata's named values with the references for (dataset, method) that name it."""
        n = 0
        for quantity, ours in values.items():
            if self._matching(dataset, method, quantity):
                n += self.check(dataset, method, quantity, ours)
        if n == 0:
            raise KeyError(f"no reference for any of {list(values)} in {dataset!r} | {method!r}")
        return n

    def cannot_read(self, words, reason: str) -> int:
        """Mark the references of these words (method names starting with any of them) as SKIP with the reason, when
        simplydata has no word yet that gives the number back (such as the numbers a chart draws)."""
        n = 0
        for i, ref in enumerate(self.references):
            if i not in self.used and any(ref["method"].startswith(w) for w in words):
                self.used.add(i)
                self.skip(ref["dataset"], ref["method"], ref["quantity"], ref["source"], reason)
                n += 1
        return n

    def finish(self) -> pd.DataFrame:
        for i, ref in enumerate(self.references):
            if i in self.used:
                continue
            reason = ref.get("reason") or "no simplydata value was read for it in the notebook"
            self.skip(ref["dataset"], ref["method"], ref["quantity"], ref["source"], reason)
        table = self.save()
        waiting = table["reason"].astype(str).str.startswith("waiting for word:").sum() if len(table) else 0
        print(f"of the SKIPs, {waiting} are waiting for a simplydata word (see the reason column)")
        where = "GitHub Actions" if os.environ.get("GITHUB_ACTIONS") == "true" else "local"
        cpu = platform.processor() or platform.machine()
        try:
            with open("/proc/cpuinfo", encoding="utf-8") as info:
                cpu = next(line.split(":", 1)[1].strip() for line in info if line.startswith("model name"))
        except (OSError, StopIteration):
            pass
        table["ran on"] = f"{where}, {platform.system()} {platform.release()}, {cpu}"
        table.to_csv(RESULTS / f"{self.notebook}.csv", index=False, float_format="%.17g")
        return table


# ---------------------------------------------------------------- for 00_summary (which reads everything with simplydata)
def result_files() -> list[str]:
    """The results file of every notebook (results/NN_*.csv), in order: 00_summary reads each with sdata.table."""
    return [str(p) for p in sorted(RESULTS.glob("[0-9][0-9]_*.csv"))]

