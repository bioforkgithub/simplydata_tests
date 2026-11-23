"""Compare simplydata results with the reference values of scripts/references/<notebook>.py (notebooks 25 to 32).

The notebook cells hold only simplydata code: they read data with sdata.table(), run a simplydata word, read the
numbers back with simplydata's own words (get_column_names, get_column, get_row, get_cell_at, through `Result` below)
and hand them to `Comparisons.check...`, which looks the reference up in references/<notebook>.csv and records PASS,
FAIL or SKIP in results/<notebook>.csv (the same record as checks.Checks). The arithmetic of a comparison (and the
reading of the reference file) happens here, out of sight of the cells.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import numpy as np
import pandas as pd

from checks import ReferenceChecks, data_file

ROOT = Path(__file__).resolve().parent.parent
REFERENCES = ROOT / "references"


def quiet(word, *args, **kwargs):
    """Run a simplydata word without its printout."""
    with contextlib.redirect_stdout(io.StringIO()):
        return word(*args, **kwargs)


class Result:
    """A simplydata result read with simplydata's own words only (get_column_names, get_column, get_row,
    get_cell_at). result[key, column] is the cell of the one row whose `key` column holds `key`; result.row(1) is the
    first row as a dict; result[name] or result.column(name) a column as a list; `where` keeps only the rows whose
    columns hold the given values (read with get_column, so the table itself is not changed)."""

    def __init__(self, table, key=None, where=None):
        self.table = table
        self.columns = list(quiet(table.get_column_names))
        keep = list(range(1, len(self._raw(self.columns[0])) + 1))
        for name, value in (where or {}).items():
            values = self._raw(name)
            keep = [p for p in keep if values[p - 1] == value]
        self.positions = keep
        self.key = key
        self.index = [self._raw(key)[p - 1] for p in keep] if key else list(range(1, len(keep) + 1))

    def _raw(self, name):
        return list(quiet(self.table.get_column, name))

    def column(self, name):
        values = self._raw(name)
        return [values[p - 1] for p in self.positions]

    def columns_of(self, names):
        """Several columns, row by row (a list of rows)."""
        return [list(values) for values in zip(*(self.column(n) for n in names))]

    def row(self, number):
        """Row `number` (counting from 1, as simplydata does) of the kept rows, as a dict (get_row)."""
        return quiet(self.table.get_row, position=self.positions[number - 1])

    def cell(self, number, name):
        return quiet(self.table.get_cell_at, row=self.positions[number - 1], column=name)

    def __getitem__(self, item):
        if isinstance(item, str):
            return self.column(item)
        key, name = item
        hits = [p for p, k in zip(self.positions, self.index) if k == key]
        if len(hits) != 1:
            raise KeyError(f"{len(hits)} rows with {self.key} = {key!r}")
        return quiet(self.table.get_cell_at, row=hits[0], column=name)

    def __contains__(self, key):
        return key in self.index

    def __len__(self):
        return len(self.positions)


class Comparisons(ReferenceChecks):
    """checks.ReferenceChecks (stream A's comparison with references/<notebook>.csv), with what notebooks 25 to 32 add:
    a missing reference is a SKIP (not an error); check_all; columns stored value by value ("name [1]", "name [2]" ...)
    compared by check_vector / check_matrix, with the tolerance types "signs" (component signs are arbitrary), "unit"
    (also scaled to length 1), "factors" (factor order matched first), "groups" (labels up to their numbering) and
    "text" (the same text); check_equal for a rerun with the same seed; waiting_for_word for a check that waits for a simplydata word."""

    def __init__(self, notebook: str):
        super().__init__(notebook)

    def data_file(self, name: str) -> str:
        return data_file(self.notebook, name)

    def has(self, dataset, method, quantity) -> bool:
        return bool(self._matching(dataset, method, quantity))

    def check(self, dataset, method, quantity, ours) -> int:
        if not self._matching(dataset, method, quantity):
            self.skip(dataset, method, quantity, "references", "no reference value for this quantity")
            return 0
        return super().check(dataset, method, quantity, ours)

    def check_all(self, dataset, method, ours: dict) -> None:
        for quantity, value in ours.items():
            self.check(dataset, method, quantity, value)

    def _column(self, dataset, method, name):
        prefix = f"{name} ["
        return [(i, r) for i, r in enumerate(self.references)
                if r["dataset"] == dataset and r["method"] == method and str(r["quantity"]).startswith(prefix)]

    def _one(self, i, ref, dataset, method, quantity, ours, reference, tolerance, how):
        self.used.add(i)
        self.compare(dataset, method, quantity, ours, reference, ref["source"], tolerance, how)
        where = ref.get("url") or ref.get("source title") or ""
        print(f"      Source: {ref['source']} ({ref.get('kind') or ''}){' - ' + where if where else ''}")
        self.rows[-1]["kind"] = ref.get("kind") or ""
        self.rows[-1]["url"] = where

    def check_vector(self, dataset, method, name, ours) -> None:
        ours = list(ours)
        found = self._column(dataset, method, name)
        if not found:
            self.skip(dataset, method, name, "references", "no reference values for this column")
            return
        by_source: dict = {}
        for i, ref in found:
            by_source.setdefault(ref["source"], []).append((i, ref))
        for source, part in by_source.items():
            how = part[0][1]["tolerance type"]
            positions = [int(str(r["quantity"])[len(name) + 2:-1]) for _, r in part]
            first_i, first = part[0]
            if max(positions) > len(ours) or (how in ("groups", "text") and len(part) != len(ours)):
                self.used.update(i for i, _ in part)
                self._one(first_i, first, dataset, method, f"{name}: number of values", len(ours), max(positions), 0, "absolute")
                continue
            if how in ("groups", "text"):
                self.used.update(i for i, _ in part)
                theirs = [r["value"] for _, r in part]
                same = (int(np.sum(_numbered(ours) == _numbered(theirs))) if how == "groups"
                        else sum(_same_text(a, b) for a, b in zip(ours, theirs)))
                label = "rows in the same group" if how == "groups" else "rows with the same value"
                self._one(first_i, first, dataset, method, f"{name}: {label}", same, len(ours), 0, "absolute")
                continue
            mine = np.asarray([ours[p - 1] for p in positions], float)
            ref = np.asarray([float(r["value"]) for _, r in part], float)
            if how == "unit":
                mine, ref = mine / np.linalg.norm(mine), ref / np.linalg.norm(ref)
            if how in ("signs", "unit", "factors"):
                mine = mine * (np.sign(np.sum(mine * ref)) or 1.0)
            kind = "relative" if how in ("signs", "unit", "factors", "relative") else how
            for (i, r), value, reference in zip(part, mine, ref):
                self._one(i, r, dataset, method, r["quantity"], value, reference, float(r["tolerance"]), kind)

    def check_matrix(self, dataset, method, name, ours_rows) -> None:
        """A matrix given row by row (Result.columns_of), checked column by column ("name, column j"). With tolerance
        type "factors" the columns of ours are first put in the reference's order (largest absolute correlation)."""
        columns = [list(c) for c in zip(*ours_rows)]
        found = [r for _, r in self._column(dataset, method, f"{name}, column 1")]
        if found and found[0]["tolerance type"] == "factors":
            refs = [np.asarray([float(r["value"]) for _, r in self._column(dataset, method, f"{name}, column {j}")], float)
                    for j in range(1, len(columns) + 1)]
            order, free = [], list(range(len(columns)))
            for ref in refs:
                best = max(free, key=lambda k: abs(np.corrcoef(np.asarray(columns[k], float), ref)[0, 1]))
                free.remove(best)
                order.append(best)
            columns = [columns[k] for k in order]
        for j, column in enumerate(columns, start=1):
            self.check_vector(dataset, method, f"{name}, column {j}", column)

    def check_equal(self, dataset, method, quantity, ours, expected, source) -> None:
        """Two simplydata numbers that must be equal (a rerun with the same seed)."""
        self.compare(dataset, method, quantity, float(ours), float(expected), source, 1e-12, "relative")
        print(f"      Source: {source}")

    def waiting_for_word(self, dataset, method, quantity, word: str) -> None:
        """A check written as it will read once simplydata has the word it needs: recorded as SKIP with the reason
        "waiting for word: <word>" (counted apart from the other SKIPs), never as a FAIL."""
        self.skip(dataset, method, quantity, "waiting for a simplydata word", f"waiting for word: {word}")


def _same_text(a, b) -> bool:
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)


def _numbered(labels):
    order: dict = {}
    return np.array([order.setdefault(v, len(order) + 1) for v in labels])
