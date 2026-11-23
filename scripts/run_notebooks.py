"""Run every notebook from top to bottom, save it with its outputs, and fail if any comparison failed.

    python scripts/run_notebooks.py               # all notebooks (01, 02, ... then 00_summary)
    python scripts/run_notebooks.py 03_regression_NIST   # only some (the summary is run again at the end)

First it runs scripts/get_data.py, which writes the datasets this repository does not host to data/fetched/
(from the R packages that ship them, and one publisher download) and checks every data file's checksum.
Each notebook saves its comparisons to results/<notebook>.csv. Every FAIL stays a FAIL in the notebooks.
The exit code is 1 when a notebook stops with an error (with SIMPLYDATA_TESTS_REQUIRE_R=1 a missing R is
such an error) or when a FAIL is NEW: not listed in results/known_fails.csv for this system (the FAILs
explained in FINDINGS.md, written by scripts/make_known_fails.py). A listed FAIL of kind A (worth fixing
in simplydata) that now passes, and that fails on every runner seen, is reported as "fixed - remove from
known_fails" (a warning). A SKIP (a
comparison that cannot be made, with its reason) does not fail the run.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import nbformat
import pandas as pd
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS = ROOT / "notebooks"
RESULTS = ROOT / "results"
sys.path.insert(0, str(NOTEBOOKS))
import checks  # noqa: E402  (the reader's words for each comparison)
KEY = ["notebook", "dataset", "method", "quantity", "reference source"]


def run_reference_script(path: Path) -> bool:
    """Notebooks with a reference script (scripts/references/<notebook>.py) get their reference values and data
    from it: it is run first (it needs R and the packages simplydata is compared against)."""
    import subprocess

    script = ROOT / "scripts" / "references" / f"{path.stem}.py"
    if not script.exists():
        return True
    done = subprocess.run([sys.executable, str(script)], cwd=ROOT)
    if done.returncode:
        print(f"ERROR: {script.name} stopped with an error")
    return done.returncode == 0


def run(path: Path) -> bool:
    if not run_reference_script(path):
        return False
    nb = nbformat.read(path, as_version=4)
    client = NotebookClient(nb, timeout=1200, kernel_name="python3",
                            resources={"metadata": {"path": str(NOTEBOOKS)}})
    started = time.time()
    ok = True
    try:
        client.execute()
    except CellExecutionError as problem:
        print(f"ERROR in {path.name}:\n{str(problem)[:3000]}")
        ok = False
    nbformat.write(nb, path)
    print(f"{path.name}: {'ran' if ok else 'stopped with an error'} in {time.time() - started:.0f} s")
    return ok


def get_data() -> bool:
    """Write the datasets that are not stored here to data/fetched/ (only those missing) and check every checksum."""
    import subprocess

    done = subprocess.run([sys.executable, str(ROOT / "scripts" / "get_data.py")], cwd=ROOT)
    if done.returncode:
        print("ERROR: scripts/get_data.py could not get or check the data (see above)")
    return done.returncode == 0


def main(names: list[str]) -> int:
    if not get_data():
        return 1
    chosen = ([NOTEBOOKS / f"{n.removesuffix('.ipynb')}.ipynb" for n in names] if names
              else sorted(p for p in NOTEBOOKS.glob("[0-9][0-9]_*.ipynb") if not p.name.startswith("00_")))
    ok = all([run(p) for p in chosen])
    summary = NOTEBOOKS / "00_summary.ipynb"
    if summary.exists():
        ok = run(summary) and ok
    files = sorted(RESULTS.glob("[0-9][0-9]_*.csv"))
    if not files:
        print("No results were saved.")
        return 1
    table = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    counts = table["status"].value_counts().to_dict()
    print(f"\n{len(table)} comparisons: " + ", ".join(f"{counts.get(s, 0)} {s}" for s in ("PASS", "FAIL", "SKIP")))
    if "result" in table.columns:
        print("what they found: " + checks.count_results(table))
    # A missing R stops the notebooks themselves when SIMPLYDATA_TESTS_REQUIRE_R=1 (see checks.r_available).
    # With notebooks named, only their results are checked: the stored results of the others may come from another
    # system (run-windows runs one notebook on Windows next to the Ubuntu results stored here).
    ran = {p.stem for p in chosen}
    new = check_against_known(table[table["notebook"].isin(ran)] if names else table)
    if not ok:
        print("\nA notebook stopped with an error.")
    return 0 if ok and not new else 1


def this_platform() -> str:
    return "windows" if sys.platform == "win32" else "ubuntu"


def check_against_known(table: pd.DataFrame) -> int:
    """Compare the FAILs with results/known_fails.csv (the FAILs explained in FINDINGS.md).

    Prints every FAIL with its finding id; gives back the number of NEW FAILs (not listed for this
    platform). A listed A-kind FAIL (worth fixing in simplydata) that now passes is reported as fixed."""
    platform = this_platform()
    known = pd.read_csv(RESULTS / "known_fails.csv")
    here = known[known["platform"].isin(["any", platform])]
    readers = [c for c in ("result", "finding", "explanation", "link") if c in table.columns]
    marked = table.drop(columns=readers).merge(here[KEY + ["finding", "same on every runner"]], on=KEY, how="left")
    failed = marked[marked["status"] == "FAIL"]
    columns = ["notebook", "dataset", "method", "quantity", "simplydata", "reference", "reference source", "finding"]
    explained, unexplained = failed[failed["finding"].notna()], failed[failed["finding"].isna()]
    print(f"\nOn {platform}: {len(failed)} FAIL, {len(explained)} of them explained in FINDINGS.md "
          f"({', '.join(f'{k} {v}' for k, v in explained['finding'].value_counts().sort_index().items())}), "
          f"{len(unexplained)} NEW.")
    if len(unexplained):
        print("\nNEW FAILs (not in results/known_fails.csv; explain them in FINDINGS.md or fix them):")
        print(unexplained[columns].to_string(index=False))
    fixed = marked[(marked["status"] == "PASS") & marked["finding"].str.startswith("A", na=False)
                   & (marked["same on every runner"] == "yes")]
    for _, row in fixed.iterrows():
        message = (f"fixed - remove from known_fails: {row['notebook']} | {row['dataset']} | {row['method']} | "
                   f"{row['quantity']} against {row['reference source']} ({row['finding']}) now passes")
        print(f"::warning::{message}" if platform == "ubuntu" else f"WARNING {message}")
    return len(unexplained)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
