"""A small test of notebooks/checks.py: a whole column of the wrong length is recorded as a difference.

    python tests/test_checks.py      (or pytest tests/)
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "notebooks"))

from checks import ReferenceChecks  # noqa: E402


def make_checks(values):
    """A ReferenceChecks with one reference row holding a whole column (no reference file is read)."""
    checks = ReferenceChecks.__new__(ReferenceChecks)
    checks.notebook, checks.rows, checks.used = "test", [], set()
    checks.references = [{"dataset": "d", "method": "m", "quantity": "q", "values": json.dumps(values), "value": None,
                          "reason": None, "source": "test source", "tolerance": 1e-10, "tolerance type": "relative",
                          "kind": "run reference", "url": ""}]
    return checks


def test_length_mismatch_is_a_difference():
    checks = make_checks([1.0, 2.0, 3.0])
    checks.check("d", "m", "q", [1.0, 2.0])               # one value short: must not raise
    row = checks.rows[-1]
    assert row["quantity"] == "q: number of values"
    assert row["simplydata"] == 2 and row["reference"] == 3 and row["status"] == "FAIL"


def test_same_length_compares_the_values():
    checks = make_checks([1.0, 2.0, 3.0])
    checks.check("d", "m", "q", [1.0, 2.0, 3.0])
    assert checks.rows[-1]["status"] == "PASS"
    assert checks.rows[-1]["quantity"] == "q: largest relative difference over 3 values"


if __name__ == "__main__":
    test_length_mismatch_is_a_difference()
    test_same_length_compares_the_values()
    print("tests/test_checks.py: 2 passed")
