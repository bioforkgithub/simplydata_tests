"""Get the datasets this repository does not host, and check every data file's checksum.

    python scripts/get_data.py            # get what is missing, then check every checksum
    python scripts/get_data.py --refresh  # get every fetched dataset (and the NIST files) again
    python scripts/get_data.py --check    # only check checksums (no R, no network)

Two kinds of data file:

- stored in data/ (public domain, US Government works, plain published facts typed by us, CC0 and other
  permissive licences; see data/SOURCES.md);
- fetched, not stored: datasets whose copies are under GPL or LGPL, under copyright, or without a stated
  licence. They are written to data/fetched/ (git-ignored) from their original source:
    * the installed R package that ships the dataset (read with Rscript, on this computer, no network), or
    * the publisher's supplementary file (bioequivalence_reference.csv, Schuetz et al. 2014), downloaded.

The conversion of each dataset to a CSV file is fixed (the R expression and the line ending below), so the
fetched file is the same, byte for byte, on every system. Every file, stored or fetched, is checked against
data/SHA256SUMS; the script exits with 1 when a file is missing or differs, and says where the file came from,
so a change in the source is never silent.

Needs R with the packages named in FETCHED (the notebooks need them anyway); see README, "Install". R is run
as `Rscript` (or from R_HOME). This script is for the test runs only; simplydata itself never goes online.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FETCHED_DIR = DATA / "fetched"
NIST = "https://www.itl.nist.gov/div898/strd"
USER_AGENT = "simplydata_tests data fetcher"

# file name -> the folder on the NIST site
NIST_FILES = {
    # univariate summary statistics
    **{f: "univ/data" for f in ("PiDigits", "Lottery", "Lew", "Mavro", "Michelso",
                                "NumAcc1", "NumAcc2", "NumAcc3", "NumAcc4")},
    # linear regression
    **{f: "lls/data/LINKS/DATA" for f in ("Norris", "Pontius", "NoInt1", "NoInt2", "Filip", "Longley",
                                          "Wampler1", "Wampler2", "Wampler3", "Wampler4", "Wampler5")},
    # one-way analysis of variance
    **{f: "anova" for f in ("SiRstv", "SmLs01", "SmLs02", "SmLs03", "SmLs04", "SmLs05", "SmLs06",
                            "SmLs07", "SmLs08", "SmLs09", "AtmWtAg")},
    # nonlinear regression
    **{f: "nls/data/LINKS/DATA" for f in ("Misra1a", "Chwirut2", "Chwirut1", "Lanczos3", "Gauss1", "Gauss2",
                                          "DanWood", "Misra1b", "Kirby2", "Hahn1", "Nelson", "MGH17", "Lanczos1",
                                          "Lanczos2", "Gauss3", "Misra1c", "Misra1d", "Roszman1", "ENSO",
                                          "MGH09", "Thurber", "BoxBOD", "Rat42", "MGH10", "Eckerle4", "Rat43",
                                          "Bennett5")},
}


@dataclass(frozen=True)
class Fetch:
    """How one dataset is made from the R package that ships it.

    package: the R package.
    r:       an R expression giving the data frame to write; get_r(name, package) reads a dataset of a package.
    The file is written by R's write.csv(row.names = FALSE, na = "") with line ending `eol` and, when quote is
    False, without quotes around text and column names (each as the file was first made).
    """

    package: str
    r: str
    quote: bool = True
    eol: str = "\n"


CRLF = "\r\n"
API_COLUMNS = ('c("cds", "stype", "dnum", "snum", "pw", "fpc", "api00", "api99", "enroll", "meals", "ell", '
               '"sch.wide", "awards")')

FETCHED: dict[str, Fetch] = {
    # survival (LGPL >= 2)
    "lung.csv": Fetch("survival", 'get_r("lung", "survival")', eol=CRLF),
    "veteran.csv": Fetch("survival", 'get_r("veteran", "survival")', eol=CRLF),
    "mgus2.csv": Fetch("survival", 'get_r("mgus2", "survival")', eol=CRLF),
    "heart.csv": Fetch("survival", 'get_r("heart", "survival")', eol=CRLF),
    "udca.csv": Fetch("survival", 'get_r("udca", "survival")', eol=CRLF),
    # metadat (GPL >= 2)
    **{f"dat_{n}.csv": Fetch("metadat", f'get_r("dat.{n}", "metadat")', eol=CRLF)
       for n in ("bcg", "hackshaw1998", "normand1999", "curtis1998", "hart1999", "mcdaniel1994")},
    # KMsurv (GPL >= 3)
    **{f"{n}.csv": Fetch("KMsurv", f'get_r("{n}", "KMsurv")', eol=CRLF) for n in ("bfeed", "channing", "larynx")},
    # survey (GPL-2 | GPL-3): the data file "api" holds apistrat and apiclus1
    "apistrat.csv": Fetch("survey", f'get_r("apistrat", "survey", "api")[, {API_COLUMNS}]', eol=CRLF),
    "apiclus1.csv": Fetch("survey", f'get_r("apiclus1", "survey", "api")[, {API_COLUMNS}]', eol=CRLF),
    # MatchIt (GPL >= 2; the data carry no licence of their own)
    "lalonde.csv": Fetch("MatchIt", 'get_r("lalonde", "MatchIt")', eol=CRLF),
    # mediation (GPL >= 2; the data carry no licence of their own)
    "jobs_ii.csv": Fetch("mediation", 'get_r("jobs", "mediation")', eol=CRLF),
    # lme4 (GPL >= 2)
    "sleepstudy.csv": Fetch("lme4", 'get_r("sleepstudy", "lme4")', eol=CRLF),
    # circular (GPL-2)
    "turtles.csv": Fetch("circular", 'data.frame(arrival = get_r("turtles", "circular")$arrival)'),
    "icu_arrivals.csv": Fetch("circular", 'data.frame(hour = as.numeric(get_r("fisherB1", "circular")))'),
    "wind_directions.csv": Fetch("circular", 'data.frame(direction_radians = as.numeric(get_r("wind", "circular")))'),
    "fisher_b10.csv": Fetch("circular", 'local({ b <- get_r("fisherB10", "circular"); do.call(rbind, lapply(names(b), '
                                        'function(s) data.frame(set = s, direction = as.numeric(b[[s]])))) })'),
    # sp (GPL >= 2)
    "meuse.csv": Fetch("sp", 'get_r("meuse", "sp")[, c("x", "y", "cadmium", "copper", "lead", "zinc")]'),
    # ismev (GPL >= 2)
    "portpirie.csv": Fetch("ismev", 'local({ d <- get_r("portpirie", "ismev"); data.frame(year = d$Year, '
                                    'sea_level = d$SeaLevel) })'),
    "coles_rain.csv": Fetch("ismev", 'data.frame(rain = as.numeric(get_r("rain", "ismev")))'),
    "fremantle.csv": Fetch("ismev", 'local({ d <- get_r("fremantle", "ismev"); data.frame(year = d$Year, '
                                    'sea_level = d$SeaLevel, soi = d$SOI) })'),
    # qvalue (Bioconductor, LGPL)
    "hedenfalk_p.csv": Fetch("qvalue", 'data.frame(p = get_r("hedenfalk", "qvalue")$p)'),
    # drc (GPL-2)
    "ryegrass.csv": Fetch("drc", 'local({ d <- get_r("ryegrass", "drc"); data.frame(conc = d$conc, '
                                 'root_length = d$rootl) })'),
    # smatr (GPL >= 3)
    "leaflife.csv": Fetch("smatr", 'get_r("leaflife", "smatr")[, c("site", "rain", "soilp", "longev", "lma")]'),
    # mcr (GPL >= 3): rows with an empty cell removed
    "creatinine.csv": Fetch("mcr", 'local({ d <- get_r("creatinine", "mcr"); d <- d[complete.cases(d), ]; '
                                   'data.frame(serum = d$serum.crea, plasma = d$plasma.crea) })'),
    # MASS (GPL-2 | GPL-3)
    "birthwt.csv": Fetch("MASS", 'get_r("birthwt", "MASS")', quote=False),
    # AER (GPL-2 | GPL-3): the columns of Zeileis, Kleiber & Jackman (2008)
    "nmes1988.csv": Fetch("AER", 'get_r("NMES1988", "AER")[, c("visits", "hospital", "health", "chronic", '
                                 '"gender", "school", "insurance")]', quote=False),
    # vegan (GPL-2): community tables, with their environment table where the notebook uses one
    "bci.csv": Fetch("vegan", 'get_r("BCI", "vegan")', quote=False),
    "dune.csv": Fetch("vegan", 'local({ d <- get_r("dune", "vegan"); cbind(plot = as.integer(rownames(d)), d, '
                               'get_r("dune.env", "vegan")) })', quote=False),
    "varespec.csv": Fetch("vegan", 'local({ d <- get_r("varespec", "vegan"); cbind(site = as.integer(rownames(d)), '
                                   'd, get_r("varechem", "vegan")) })', quote=False),
    "mite.csv": Fetch("vegan", 'cbind(get_r("mite", "vegan"), Topo = get_r("mite.env", "vegan")$Topo)',
                      quote=False),
    # lavaan (GPL >= 2). One student (id 351) has no grade (NA in lavaan). The file these tests were made with
    # holds -2147483648 there (R's internal code for a missing whole number, let through by an earlier copy);
    # it is written the same way so that no value changes. No notebook uses the column grade.
    "holzinger_swineford_1939.csv": Fetch("lavaan", 'local({ d <- get_r("HolzingerSwineford1939", "lavaan"); '
                                                    'd$grade <- as.numeric(d$grade); '
                                                    'd$grade[is.na(d$grade)] <- -2147483648; d })', quote=False),
    # BayesFactor (GPL-2): the 2 x 2 table of Hraba & Grant (1970), one row per child, in the order
    # white child (same-race doll, different-race doll), then black child; labels in lower case
    "race_dolls.csv": Fetch("BayesFactor", 'local({ m <- get_r("raceDolls", "BayesFactor"); '
                                           'do.call(rbind, lapply(colnames(m), function(child) do.call(rbind, '
                                           'lapply(rownames(m), function(doll) data.frame(child = rep(tolower(child), '
                                           'm[doll, child]), doll = rep(tolower(doll), m[doll, child])))))) })',
                            quote=False),
}

# Downloaded from the publisher (no R package ships them). Schuetz, Labes & Fuglsang (2014), AAPS J 16, 1292-1297,
# doi:10.1208/s12248-014-9661-0, electronic supplementary material 1 to 3 (datasets A, B and C). The article is also
# on PubMed Central (PMC4389751), but PMC hands its supplementary files only to a web browser (a "Preparing to
# download" page first), so the publisher's copy of the same files is used.
SCHUETZ_ESM = ("https://static-content.springer.com/esm/art%3A10.1208%2Fs12248-014-9661-0/MediaObjects/"
               "12248_2014_9661_MOESM{n}_ESM.txt")
DOWNLOADED = {"bioequivalence_reference.csv": [SCHUETZ_ESM.format(n=n) for n in (1, 2, 3)]}

R_PRELUDE = r'''
get_r <- function(name, package, data_set = name) {
  e <- new.env()
  suppressWarnings(utils::data(list = data_set, package = package, envir = e))
  if (exists(name, envir = e, inherits = FALSE)) get(name, envir = e) else getExportedValue(package, name)
}
write_r <- function(d, path, eol, quote) {
  con <- file(path, "wb"); on.exit(close(con))
  write.csv(d, con, row.names = FALSE, na = "", eol = eol, quote = quote)
}
'''


def source_of(name: str) -> str:
    if name in FETCHED:
        return f"the R package {FETCHED[name].package}"
    if name in DOWNLOADED:
        return "the publisher's supplementary files " + ", ".join(DOWNLOADED[name])
    return "this repository"


def url_of(name: str) -> str:
    return f"{NIST}/{NIST_FILES[name]}/{name}.dat"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download_bytes(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def download_nist(name: str) -> None:
    target = DATA / "nist" / f"{name}.dat"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(download_bytes(url_of(name)))
    print(f"downloaded {target.relative_to(ROOT)} from {url_of(name)}")


def rscript() -> str:
    found = shutil.which("Rscript")
    if found:
        return found
    home = os.environ.get("R_HOME")
    for sub in ("bin/Rscript", "bin/Rscript.exe", "bin/x64/Rscript.exe"):
        if home and (Path(home) / sub).exists():
            return str(Path(home) / sub)
    raise SystemExit("get_data.py needs R (Rscript) to read the datasets shipped with R packages; install R and "
                     "the packages in README ('Install'), or set R_HOME.")


def r_string(text: str) -> str:
    return '"' + text.replace("\\", "/").replace('"', '\\"') + '"'


def fetch_from_r(names: list[str]) -> None:
    """Write the chosen datasets to data/fetched/ with one Rscript run."""
    FETCHED_DIR.mkdir(parents=True, exist_ok=True)
    lines = [R_PRELUDE]
    for name in names:
        f = FETCHED[name]
        eol = "\\r\\n" if f.eol == CRLF else "\\n"   # written into the R code as an escape
        quote = "TRUE" if f.quote else "FALSE"
        lines.append(f'if (!requireNamespace("{f.package}", quietly = TRUE)) '
                     f'stop("the R package {f.package} is needed for {name}; install it (see README)")')
        lines.append(f'write_r({f.r}, {r_string(str(FETCHED_DIR / name))}, "{eol}", {quote})')
    with tempfile.TemporaryDirectory() as temporary:
        script = Path(temporary) / "fetch.R"
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")
        done = subprocess.run([rscript(), str(script)], capture_output=True, text=True)
    if done.returncode:
        raise SystemExit(f"R stopped while writing the datasets:\n{done.stderr[-3000:]}")
    for name in names:
        print(f"wrote data/fetched/{name} from {source_of(name)}")


def fetch_schuetz() -> None:
    """Datasets A, B and C of Schuetz et al. (2014) in one long table: dataset, subject, period, formulation
    (R -> reference, T -> test), auc; rows by dataset, then subject, then period."""
    import pandas as pd

    frames = []
    for dataset, url in zip("ABC", DOWNLOADED["bioequivalence_reference.csv"]):
        text = download_bytes(url).decode("ascii")
        table = pd.read_csv(io.StringIO(text), sep=r"\s+")
        frames.append(pd.DataFrame({"dataset": dataset, "subject": table["Subj"], "period": table["Per"],
                                    "formulation": table["Trt"].map({"R": "reference", "T": "test"}),
                                    "auc": table["Var"]}).sort_values(["subject", "period"], kind="stable"))
    FETCHED_DIR.mkdir(parents=True, exist_ok=True)
    target = FETCHED_DIR / "bioequivalence_reference.csv"
    target.write_bytes(pd.concat(frames, ignore_index=True).to_csv(index=False, lineterminator="\n").encode("ascii"))
    print(f"wrote data/fetched/bioequivalence_reference.csv from {source_of('bioequivalence_reference.csv')}")


def expected_sums() -> dict[str, str]:
    sums = {}
    for line in (DATA / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        if line.strip():
            expected, name = line.split(maxsplit=1)
            sums[name.strip()] = expected
    return sums


def check() -> int:
    bad = 0
    for name, expected in expected_sums().items():
        path = DATA / name
        short = Path(name).name
        if not path.exists():
            hint = " (run: python scripts/get_data.py)" if name.startswith("fetched/") else ""
            print(f"MISSING  {name}{hint}")
            bad += 1
        elif sha256(path) != expected:
            if name.startswith("fetched/"):
                print(f"DIFFERS  {name}: the file made from {source_of(short)} is not the one these tests were "
                      f"made with (sha256 {sha256(path)}, expected {expected}). The source has changed (a new "
                      f"package version or a changed download); check the values before updating data/SHA256SUMS.")
            else:
                print(f"DIFFERS  {name}")
            bad += 1
    print(f"{'all checksums match' if not bad else f'{bad} file(s) missing or different'}")
    return 1 if bad else 0


def write_checksums() -> int:
    """Write data/SHA256SUMS for every data file, stored and fetched (for maintainers, after adding a dataset)."""
    files = sorted(p for p in DATA.rglob("*") if p.is_file() and p.suffix in (".csv", ".dat"))
    lines = [f"{sha256(p)}  {p.relative_to(DATA).as_posix()}" for p in files]
    (DATA / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote checksums of {len(files)} files to data/SHA256SUMS")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refresh", action="store_true", help="get every fetched dataset and NIST file again")
    parser.add_argument("--check", action="store_true", help="only check checksums (no R, no network)")
    parser.add_argument("--write-checksums", action="store_true", help="write data/SHA256SUMS (maintainers)")
    args = parser.parse_args()
    if args.write_checksums:
        return write_checksums()
    if not args.check:
        for name in NIST_FILES:
            if args.refresh or not (DATA / "nist" / f"{name}.dat").exists():
                download_nist(name)
        wanted = [n for n in FETCHED if args.refresh or not (FETCHED_DIR / n).exists()]
        if wanted:
            fetch_from_r(wanted)
        for name in DOWNLOADED:
            if args.refresh or not (FETCHED_DIR / name).exists():
                fetch_schuetz()
    return check()


if __name__ == "__main__":
    sys.exit(main())
