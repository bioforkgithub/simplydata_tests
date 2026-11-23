"""Fill the reader's columns (result, finding, explanation, link) of results/<notebook>.csv from the machine column
"status" and results/known_fails.csv, without running the notebooks again.

    python scripts/label_results.py [windows|ubuntu]

The words are those of notebooks/checks.py ("agrees", "simplydata error - being fixed", "explained difference -
reference / machine / precision limit", "difference not yet explained", "waiting for word", "skipped"). The system
decides which rows of known_fails.csv apply; by default the one the results came from, read from their "ran on"
column (GitHub Actions runs are Ubuntu), else this computer's. Nothing else in the files changes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "notebooks"))

import checks  # noqa: E402


def system_of(table: pd.DataFrame, default: str) -> str:
    if "ran on" in table.columns and len(table):
        ran = str(table["ran on"].iloc[0])
        if "GitHub Actions" in ran or "Linux" in ran:
            return "ubuntu"
        if "Windows" in ran:
            return "windows"
    return default


def main(args: list[str]) -> int:
    forced = args[0] if args else None
    for path in sorted((ROOT / "results").glob("[0-9][0-9]_*.csv")):
        table = pd.read_csv(path, keep_default_na=False, na_values=[""])
        system = forced or system_of(table, checks.this_system())
        checks.this_system = lambda s=system: s
        checks._known = None
        table = checks.add_result_words(table)
        table.to_csv(path, index=False, float_format="%.17g")
        print(f"{path.name} ({system}): {checks.count_results(table)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
