"""Read a dataset that ships with an installed R package, at run time, into a pandas table.

The reference scripts of notebooks 25 to 32 use the datasets of the R help pages and vignettes they reproduce. The data are read from
the installed R package each time (nothing is copied into this repository), through rpy2, with no network.
Factors become text columns; their R level order is kept in `levels[column]` so that a notebook can give the
same order to simplydata (in_order=) or to R. Logical columns become True/False.
"""

from __future__ import annotations

import pandas as pd

from checks import r_available


def r_dataset(name: str, package: str = "datasets", levels: dict | None = None) -> pd.DataFrame:
    """The R dataset `package::name` as a pandas table (row names, when they are not 1, 2, ..., in column "row").

    levels: an optional dict that is filled with the R factor levels of each factor column."""
    if not r_available():
        raise RuntimeError("R is needed to read the R datasets of this notebook")
    import checks

    ro = checks._R
    # data() finds most datasets; lazy-loaded ones kept in a shared file (survival's mgus2) come from package::name
    ro.r(f'.sd_tmp <- as.data.frame(tryCatch(get(data("{name}", package = "{package}", envir = environment())), '
         f'error = function(e) {package}::`{name}`, warning = function(w) {package}::`{name}`))')
    columns = list(ro.r("names(.sd_tmp)"))
    out = {}
    for column in columns:
        ro.r(f'.sd_col <- .sd_tmp[["{column}"]]')
        if bool(ro.r("is.factor(.sd_col)")[0]):
            values = [None if v == "NA_character_" else str(v) for v in ro.r("as.character(.sd_col)")]
            na = list(ro.r("is.na(.sd_col)"))
            out[column] = [None if missing else v for v, missing in zip(values, na)]
            if levels is not None:
                levels[column] = [str(v) for v in ro.r("levels(.sd_col)")]
        elif bool(ro.r("is.logical(.sd_col)")[0]):
            na = list(ro.r("is.na(.sd_col)"))
            values = list(ro.r("as.integer(.sd_col)"))
            out[column] = [None if missing else bool(v) for v, missing in zip(values, na)]
        elif bool(ro.r("is.numeric(.sd_col)")[0]):
            na = list(ro.r("is.na(.sd_col)"))
            values = list(ro.r("as.numeric(.sd_col)"))
            out[column] = [float("nan") if missing else float(v) for v, missing in zip(values, na)]
        else:
            na = list(ro.r("is.na(.sd_col)"))
            values = list(ro.r("as.character(.sd_col)"))
            out[column] = [None if missing else str(v) for v, missing in zip(values, na)]
    table = pd.DataFrame(out)
    names = [str(v) for v in ro.r("rownames(.sd_tmp)")]
    if names != [str(i) for i in range(1, len(table) + 1)]:
        table.insert(0, "row", names)
    ro.r("rm(.sd_tmp, .sd_col)")
    return table


def r_series(name: str, package: str = "datasets") -> pd.DataFrame:
    """An R time series (ts) as a table with columns "time" (R time(), as a number) and "value"."""
    if not r_available():
        raise RuntimeError("R is needed to read the R datasets of this notebook")
    import checks

    ro = checks._R
    ro.r(f'.sd_ts <- tryCatch(get(data("{name}", package = "{package}", envir = environment())), '
         f'error = function(e) {package}::`{name}`, warning = function(w) {package}::`{name}`)')
    values = [float(v) for v in ro.r("as.numeric(.sd_ts)")]
    times = [float(v) for v in ro.r("as.numeric(time(.sd_ts))")]
    frequency = float(ro.r("frequency(.sd_ts)")[0])
    ro.r("rm(.sd_ts)")
    table = pd.DataFrame({"time": times, "value": values})
    table.attrs["frequency"] = frequency
    return table
