"""Helpers for notebooks 18 to 24: read real datasets from the installed packages at run time, and run the
packages' own manual examples (R help pages, Python docstrings).

Nothing is copied into this repository: each dataset is read from the R package (or Python package) that
ships it, when the notebook runs, and each manual example is read from the installed help page or
docstring. Nothing is sent anywhere; everything runs on this computer.
"""

from __future__ import annotations

import contextlib
import doctest
import inspect
import io
import re

import numpy as np
import pandas as pd

import checks


def r_dataset(name: str, package: str = "datasets", data_set: str | None = None) -> pd.DataFrame:
    """An R dataset as a pandas DataFrame: factors become text, NA becomes an empty cell (NaN or None).

    Row names, when they are not 1, 2, 3 ..., become a first column called "row name". data_set= names the
    data file when it holds several tables (R survey's `api` holds apistrat, apiclus1 ...)."""
    ro = checks._R
    code = f"""local({{
        e <- new.env(); utils::data(list = "{data_set or name}", package = "{package}", envir = e)
        d <- if (exists("{name}", envir = e, inherits = FALSE)) get("{name}", envir = e)
             else getExportedValue("{package}", "{name}")
        d <- as.data.frame(d, stringsAsFactors = FALSE)
        rn <- rownames(d)
        out <- lapply(d, function(v) if (is.factor(v)) as.character(v) else if (inherits(v, "Date")) as.character(v) else as.vector(v))
        list(columns = out, rownames = rn, plain = identical(rn, as.character(seq_len(nrow(d)))))
    }})"""
    found = ro.r(code)
    columns = found.rx2("columns")
    table = {}
    for col_name, values in zip(columns.names, columns):
        r_class = ro.r["class"](values)[0]
        if r_class == "character":
            table[col_name] = [None if (v is ro.NA_Character or str(v) == "NA") else str(v) for v in values]
        elif r_class == "logical":
            table[col_name] = [None if v is ro.NA_Logical else bool(v) for v in values]
        else:
            numbers = np.array(ro.r["as.double"](values), dtype=float)
            numbers[np.array(ro.r["is.na"](values), dtype=bool)] = np.nan
            table[col_name] = numbers
    frame = pd.DataFrame(table)
    if not bool(found.rx2("plain")[0]):
        frame.insert(0, "row name", [str(v) for v in found.rx2("rownames")])
    return frame


def r_example_objects(package: str, topic: str) -> None:
    """Run the assignments (x <- ...) of an R help page's examples, so the objects it builds exist in R.

    Only lines of the form `name <- value` are run (no plots, no printing), and only the first assignment to
    each name (a help page may use x for its first example and again for a later one); each is run in R's
    global environment and an assignment that fails is left out. Gives back the names that were made."""
    ro = checks._R
    code = f"""local({{
        db <- tools::Rd_db("{package}")
        rd <- db[[paste0("{topic}", ".Rd")]]
        f <- tempfile(fileext = ".R"); tools::Rd2ex(rd, f, commentDontrun = TRUE, commentDonttest = FALSE)
        ex <- parse(f, keep.source = FALSE); unlink(f)
        made <- character(0)
        for (e in ex) if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
                          is.name(e[[2]]) && !(as.character(e[[2]]) %in% made)) {{
            ok <- tryCatch({{ eval(e, globalenv()); TRUE }}, error = function(err) FALSE)
            if (ok) made <- c(made, as.character(e[[2]]))
        }}
        made
    }})"""
    return [str(v) for v in ro.r(code)]


def r_example_code(package: str, topic: str) -> str:
    """The text of an R help page's examples section (to print beside the comparison)."""
    ro = checks._R
    code = f"""local({{
        rd <- tools::Rd_db("{package}")[[paste0("{topic}", ".Rd")]]
        f <- tempfile(fileext = ".R"); tools::Rd2ex(rd, f); x <- paste(readLines(f), collapse = "\n"); unlink(f); x
    }})"""
    return str(ro.r(code)[0])


def docstring_example(function, skip_plots: bool = True) -> tuple[dict, list[tuple[str, str]]]:
    """Run the examples of a Python docstring, as its package's own documentation gives them.

    Gives back the namespace the examples made (their data and results) and every example as
    (source, the output the docstring prints). Examples that draw (matplotlib) are left out."""
    import scipy
    from scipy import stats

    text = inspect.getdoc(function) or ""
    examples = doctest.DocTestParser().get_examples(text)
    names = {"np": np, "stats": stats, "scipy": scipy, "pd": pd}
    shown = []
    for example in examples:
        source = example.source
        if skip_plots and re.search(r"\b(plt|ax|fig|matplotlib)\b", source):
            continue
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(source, "<docstring example>", "exec"), names)
        except Exception:  # noqa: BLE001 - an example that needs something not here is left out
            continue
        shown.append((source.strip(), example.want.strip()))
    return names, shown


def printed_number(text: str, name: str) -> float:
    """A number the docstring prints after `name=` (as in `statistic=24.0`), or after `name` alone."""
    m = re.search(re.escape(name) + r"\s*=?\s*(?:np\.float64\()?(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", text)
    if not m:
        raise ValueError(f"{name} not printed in: {text[:200]}")
    return float(m.group(1))


def r_version_of(*packages: str) -> str:
    ro = checks._R
    found = [str(ro.r('as.character(packageVersion("' + p + '"))')[0]) for p in packages]
    return ", ".join(f"{p} {v}" for p, v in zip(packages, found))


def run_printed(word, *args, **options):
    """Run a simplydata word and keep what it prints: (its result, the printed text)."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        result = word(*args, **options)
    return result, buffer.getvalue()


def compare_printed(checks_, dataset, method, quantity, shown: str, reference: float, reference_name: str):
    """Compare a number as simplydata prints it (a string such as "0.108" or "< 0.0001") with a reference:
    within half a unit of the last printed digit, the rule for printed values in these notebooks."""
    shown = shown.strip()
    if shown.startswith("<"):
        bound = float(shown[1:])
        return checks_.compare(dataset, method, f"{quantity} (printed {shown})", 0.0, reference,
                               reference_name + " (below the printed bound)", bound)
    decimals = len(shown.split(".")[1]) if "." in shown and "e" not in shown.lower() else 0
    allowed = 0.5 * 10 ** (-decimals) if "e" not in shown.lower() else 0.5 * 10 ** (
        int(shown.lower().split("e")[1]) - len(shown.lower().split("e")[0].replace("-", "").replace(".", "")) + 1)
    return checks_.compare(dataset, method, f"{quantity} (printed)", float(shown), reference,
                           reference_name + ", half a unit of the printed digit", allowed * (1 + 1e-9))


def printed_p(text: str) -> str:
    """The p-value as a simplydata sentence prints it ("0.108", "0.0026" or "< 0.0001")."""
    m = re.search(r"\bp (=|<) ([0-9.]+(?:e-?\d+)?)", text)
    if not m:
        raise ValueError("no p-value printed in:\n" + text[:500])
    return m.group(2) if m.group(1) == "=" else "< " + m.group(2)


def machine() -> str:
    """Where the notebook ran: GitHub runner or local, with the CPU (some results near a tolerance depend on it)."""
    import os
    import platform

    cpu = platform.processor() or platform.machine()
    try:
        with open("/proc/cpuinfo", encoding="utf-8") as info:
            cpu = next(line.split(":", 1)[1].strip() for line in info if line.startswith("model name"))
    except (OSError, StopIteration):
        pass
    where = "GitHub Actions" if os.environ.get("GITHUB_ACTIONS") == "true" else "local"
    return f"{where}, {platform.system()} {platform.release()}, {cpu}"


def save(checks_) -> pd.DataFrame:
    """checks.save(), with a column "ran on" (GitHub runner or local, and the CPU)."""
    table = checks_.save()
    table["ran on"] = machine()
    table.to_csv(checks.RESULTS / f"{checks_.notebook}.csv", index=False, float_format="%.17g")
    print("ran on:", table["ran on"].iloc[0] if len(table) else machine())
    return table


def docstring_data(function, name: str):
    """The value the docstring's examples first give to `name` (as in `d = [6, 8, 14, ...]`), run alone."""
    import scipy
    from scipy import stats

    text = inspect.getdoc(function) or ""
    for example in doctest.DocTestParser().get_examples(text):
        if re.match(rf"{re.escape(name)}\s*=(?!=)", example.source):
            names = {"np": np, "stats": stats, "scipy": scipy}
            exec(compile(example.source, "<docstring example>", "exec"), names)
            return names[name]
    raise KeyError(f"{name} is not set in the docstring of {function.__name__}")


def docstring_printed(function, source_start: str, occurrence: int = 0) -> str:
    """What the docstring prints for its `occurrence`-th example whose source starts with source_start."""
    text = inspect.getdoc(function) or ""
    found = [e.want.strip() for e in doctest.DocTestParser().get_examples(text) if e.source.startswith(source_start)]
    return found[occurrence]
