"""Build the reliability report of simplydata_tests as one PDF, from the files in this repository.

    python scripts/make_pdf_report.py                       # writes reports/simplydata_reliability_report.pdf
    python scripts/make_pdf_report.py --out somewhere.pdf

It reads only what the notebooks left behind: results/NN_*.csv (every comparison and its label),
references/NN_*.csv (where each reference value came from, with the package versions), FINDINGS.md (the evidence
for each explained difference), requirements.txt (the pinned Python packages) and README.md (the run the results
come from). The coverage of words needs simplydata installed (requirements.txt pins the commit the notebooks ran
on); without it that section says so. Nothing is computed again: the PDF reports the stored numbers.

The PDF is made with ReportLab (pure Python, version pinned in .github/workflows/pdf-report.yml) and the DejaVu
fonts that matplotlib ships, so it needs no system library and builds the same way on any machine. ReportLab's
invariant mode leaves the build time out of the file: the same inputs give the same PDF.
"""

from __future__ import annotations

import argparse
import collections
import csv
import glob
import re
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "reports" / "simplydata_reliability_report.pdf"

MOTTO = ("Learn basic Python, and simplydata takes you from your first table to sound statistics, in plain words, "
         "with its results checked against the field's standard tools.")
SECOND_LINE = ("simplydata is a package that is easy for a beginner, honest with a reviewer, and safe for sensitive "
               "data, all at once.")

# the labels of the column "result", in the order the tables show them, with what each means (our own words)
LABELS = [
    ("agrees", "agrees",
     "simplydata and the reference differ by no more than the stated tolerance."),
    ("explained difference – reference", "reference",
     "outside the tolerance because of the reference: it rounds its printed numbers, stops its iterations early, "
     "follows another definition, or has a known fault of its own. FINDINGS.md names which, under a C finding."),
    ("explained difference – machine", "machine",
     "outside the tolerance in the last digits only on some processors or systems; the same comparison agrees "
     "elsewhere (B findings seen on one kind of runner, and B7)."),
    ("explained difference – precision limit", "precision",
     "outside the tolerance on every system because double-precision arithmetic, or the package simplydata calls, "
     "cannot do better on this problem; the reference shares the limit or reaches it differently (B findings, "
     "C29 and C30)."),
    ("simplydata error – being fixed", "being fixed",
     "a fault in simplydata (an A finding) that is waiting for its fix to be released; it becomes \"agrees\" then."),
    ("difference not yet explained", "unexplained",
     "outside the tolerance with no finding to explain it. The notebook workflow fails on any such row, so a "
     "published run has none."),
    ("waiting for word", "waiting",
     "the number exists in the reference, but simplydata had no word that gives it back yet."),
    ("skipped", "skipped",
     "the comparison could not be made; the column \"reason\" says why."),
]
SHORT = {label: short for label, short, _ in LABELS}
GROUPS = [("01 to 08", 1, 8), ("09 to 17", 9, 17), ("18 to 24", 18, 24), ("25 to 32", 25, 32)]
KINDS = {"certified": "a certified value (NIST StRD, NIST/SEMATECH e-Handbook)",
         "printed": "a number printed in a paper, book, manual or help page",
         "run reference": "a run of another program (R and its packages, lifelines, SciPy or statsmodels called "
                          "directly), stored once with its versions",
         "expected": "a value known by construction",
         "": "no kind recorded (skipped rows and rows compared within simplydata)"}


# ---------------------------------------------------------------- reading the repository
def read_csv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def notebooks() -> list[str]:
    return sorted(Path(p).stem for p in glob.glob(str(ROOT / "results" / "[0-9][0-9]_*.csv")))


def results() -> dict[str, list[dict]]:
    return {nb: read_csv(ROOT / "results" / f"{nb}.csv") for nb in notebooks()}


def references() -> list[dict]:
    rows = []
    for p in sorted(glob.glob(str(ROOT / "references" / "[0-9][0-9]_*.csv"))):
        rows += read_csv(Path(p))
    return rows


def word_of(method: str) -> str:
    """The simplydata word a method label starts with ("compare_averages_of (ANOVA)" -> compare_averages_of;
    "chart lag_dots (rows_back=1)" -> lag_dots; "after use_survey_design (strata), ..." -> use_survey_design)."""
    tokens = [x for x in re.split(r"[ (,:]", method.strip()) if x]
    if len(tokens) > 1 and tokens[0] in ("chart", "after"):
        return tokens[1]
    return tokens[0] if tokens else ""


def findings_sections() -> list[tuple[str, str]]:
    """FINDINGS.md cut at its ## and ### headings: (heading, body)."""
    text = (ROOT / "FINDINGS.md").read_text(encoding="utf-8")
    parts = re.split(r"(?m)^(#{2,3} .+)$", text)
    return [(parts[i].lstrip("#").strip(), parts[i + 1].strip()) for i in range(1, len(parts) - 1, 2)]


def evidence_for(fid: str, sections) -> tuple[str, list[str]]:
    """The title of a finding and the pieces of FINDINGS.md about it: its own section ("### B1. ..."), the sections
    whose heading names it among others ("### B6 and B1, more cases"), and the table rows that start with it (the
    table of C1 to C9)."""
    name = re.compile(rf"(?<![A-Za-z0-9]){re.escape(fid)}(?![0-9])")
    own = re.compile(r"^[A-Z]\d+\. ")
    title, pieces = "", []
    for heading, text in sections:
        if heading.startswith(f"{fid}. "):
            title = title or heading[len(fid) + 2:].strip()
            pieces.insert(0, text)
        elif not own.match(heading) and name.search(heading):
            pieces.append(f"**{heading}**\n\n{text}")
        for line in text.splitlines():
            cells = [x.strip() for x in line.strip().strip("|").split("|")]
            if line.startswith("|") and len(cells) > 1 and cells[0] == fid:
                title = title or cells[1]
                pieces.append("; ".join(cells[1:]))
    return title, pieces


def requirements() -> list[str]:
    return [line.strip() for line in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")]


def run_described() -> str:
    """The sentence of README.md that names the run the stored results come from."""
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    found = re.search(r"The results files as stored here \(([^)]*\)[^)]*)\)", text)
    return re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", " ".join(found.group(1).split())) if found else ""


def versions(refs: list[dict]) -> list[tuple[str, str, int]]:
    """(package, version, reference values made with it) from the column "versions" of the reference files."""
    seen: dict[tuple[str, str], int] = collections.Counter()
    for r in refs:
        for part in (r.get("versions") or "").split(";"):
            part = part.strip()
            if not part:
                continue
            if part.startswith("R version"):
                seen[("R", part.removeprefix("R version").strip())] += 1
                continue
            name, _, version = part.rpartition(" ")
            seen[(name or part, version)] += 1
    return sorted(((n, v, c) for (n, v), c in seen.items()), key=lambda x: (x[0] != "Python", x[0] != "R", x[0].lower()))


def simplydata_words():
    """Every simplydata word with its group, from the installed simplydata; None without it."""
    try:
        from simplydata import learn
    except Exception:  # noqa: BLE001 - the section says the words could not be listed
        return None
    words = []
    for group, keys in learn.GROUPS:
        for key in keys:
            if learn._kind(key) != "option":
                words.append((group, key))
    return words


def match_word(head: str, keys: set[str]) -> str | None:
    """The word key a method's first token names: a table word first (table.forecast before the chart's forecast),
    then the token itself (chart words), then the sdata., dist. and nums. words."""
    if "." in head:
        return head if head in keys else None
    for key in (f"table.{head}", head, f"sdata.{head}", f"dist.{head}", f"nums.{head}"):
        if key in keys:
            return key
    return None


# ---------------------------------------------------------------- the PDF
def build(out: Path) -> None:
    import matplotlib
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import mm
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.platypus import (LongTable, PageBreak, Paragraph, Preformatted, SimpleDocTemplate, Spacer,
                                    TableStyle)
    from reportlab.platypus.tableofcontents import TableOfContents

    fonts = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    for name, file in (("Sans", "DejaVuSans.ttf"), ("Sans-Bold", "DejaVuSans-Bold.ttf"),
                       ("Sans-Oblique", "DejaVuSans-Oblique.ttf"), ("Sans-BoldOblique", "DejaVuSans-BoldOblique.ttf"),
                       ("Mono", "DejaVuSansMono.ttf")):
        pdfmetrics.registerFont(TTFont(name, str(fonts / file)))
    pdfmetrics.registerFontFamily("Sans", normal="Sans", bold="Sans-Bold", italic="Sans-Oblique",
                                  boldItalic="Sans-BoldOblique")

    base = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=base["BodyText"], fontName="Sans", fontSize=9, leading=12)
    small = ParagraphStyle("small", parent=body, fontSize=7, leading=8.6)
    cell = ParagraphStyle("cell", parent=body, fontSize=6.6, leading=7.8)
    h1 = ParagraphStyle("h1", parent=base["Heading1"], fontName="Sans-Bold", fontSize=15, leading=19,
                        spaceBefore=6, spaceAfter=8)
    h2 = ParagraphStyle("h2", parent=base["Heading2"], fontName="Sans-Bold", fontSize=11.5, leading=15,
                        spaceBefore=10, spaceAfter=5)
    h3 = ParagraphStyle("h3", parent=base["Heading3"], fontName="Sans-Bold", fontSize=9.5, leading=12.5,
                        spaceBefore=8, spaceAfter=3)
    title_style = ParagraphStyle("title", parent=h1, fontSize=22, leading=27, alignment=TA_CENTER, spaceAfter=14)
    centre = ParagraphStyle("centre", parent=body, alignment=TA_CENTER, fontSize=10.5, leading=14)
    mono = ParagraphStyle("mono", parent=body, fontName="Mono", fontSize=7, leading=8.5)

    def p(text, style=body):
        return Paragraph(text, style)

    def md(text: str) -> str:
        """Inline Markdown of FINDINGS.md to ReportLab's markup: code, bold, italics and links."""
        text = escape(text)
        text = re.sub(r"`([^`]+)`", r'<font name="Mono">\1</font>', text)
        text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
        text = re.sub(r"(?<![*\w])\*([^*]+)\*(?![*\w])", r"<i>\1</i>", text)
        text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)
        return text

    def md_blocks(text: str) -> list:
        """FINDINGS.md text as flowables: paragraphs, bullet lists, tables and code blocks."""
        flow, lines, table, code = [], [], [], None

        def flush_lines():
            item, items = [], []
            for line in lines:
                if line.lstrip().startswith(("- ", "* ")):
                    if item:
                        items.append(item)
                    item = ["- ", line.lstrip()[2:]]
                else:
                    item = item + [line.strip()] if item else ["", line.strip()]
            if item:
                items.append(item)
            for it in items:
                text_ = " ".join(x for x in it[1:] if x)
                if it[0]:
                    flow.append(Paragraph(md(text_), body, bulletText="•"))
                else:
                    flow.append(p(md(text_)))
            lines.clear()

        def flush_table():
            if table:
                rows = [[x.strip() for x in r.strip().strip("|").split("|")] for r in table
                        if not re.match(r"^\|?\s*:?-{2,}", r.strip())]
                width = max(len(r) for r in rows)
                rows = [r + [""] * (width - len(r)) for r in rows]
                flow.append(grid([[Paragraph(md(x), cell) for x in r] for r in rows], header=True))
                table.clear()

        for line in text.splitlines():
            if line.strip().startswith("```"):
                if code is None:
                    flush_lines()
                    flush_table()
                    code = []
                else:
                    flow.append(Preformatted("\n".join(code), mono))
                    code = None
                continue
            if code is not None:
                code.append(line)
            elif line.strip().startswith("|"):
                flush_lines()
                table.append(line)
            elif not line.strip():
                flush_lines()
                flush_table()
            else:
                flush_table()
                if line.lstrip().startswith(("- ", "* ")) and lines and not any(
                        x.lstrip().startswith(("- ", "* ")) for x in lines):
                    flush_lines()  # a list right after a paragraph
                lines.append(line)
        flush_lines()
        flush_table()
        if code:
            flow.append(Preformatted("\n".join(code), mono))
        return flow

    def grid(rows, header=True, widths=None, numeric_from=None):
        t = LongTable(rows, colWidths=widths, repeatRows=1 if header else 0)
        style = [("FONTNAME", (0, 0), (-1, -1), "Sans"), ("FONTSIZE", (0, 0), (-1, -1), 6.6),
                 ("LEADING", (0, 0), (-1, -1), 7.8), ("VALIGN", (0, 0), (-1, -1), "TOP"),
                 ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#9aa3ad")),
                 ("LEFTPADDING", (0, 0), (-1, -1), 2.5), ("RIGHTPADDING", (0, 0), (-1, -1), 2.5),
                 ("TOPPADDING", (0, 0), (-1, -1), 1.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 1.5)]
        if header:
            style += [("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e4e8ec")),
                      ("FONTNAME", (0, 0), (-1, 0), "Sans-Bold")]
        if numeric_from is not None:
            style.append(("ALIGN", (numeric_from, 1), (-1, -1), "RIGHT"))
        t.setStyle(TableStyle(style))
        return t

    def c(text) -> Paragraph:
        return Paragraph(escape(str(text)), cell)

    head_cell = ParagraphStyle("head", parent=cell, fontName="Sans-Bold", fontSize=5.8, leading=7)

    def h(text) -> Paragraph:
        return Paragraph(escape(str(text)), head_cell)

    def num(value: str) -> str:
        """A stored number as it is, shortened to 12 significant digits for the table."""
        try:
            x = float(value)
        except (TypeError, ValueError):
            return value or ""
        return f"{x:.12g}"

    class Doc(SimpleDocTemplate):
        def beforeDocument(self):
            self.seq = 0  # multiBuild makes several passes: the same keys on each, so the contents links hold

        def afterFlowable(self, flowable):
            if isinstance(flowable, Paragraph) and flowable.style.name in ("h1", "h2"):
                level = 0 if flowable.style.name == "h1" else 1
                text = flowable.getPlainText()
                key = f"h{self.seq}"
                self.seq += 1
                self.canv.bookmarkPage(key)
                self.canv.addOutlineEntry(text, key, level=level, closed=level > 0)
                self.notify("TOCEntry", (level, text, self.page, key))

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Sans", 7)
        canvas.drawString(18 * mm, 10 * mm, "simplydata reliability report (simplydata_tests)")
        canvas.drawRightString(A4[0] - 18 * mm, 10 * mm, f"page {doc.page}")
        canvas.restoreState()

    # ---- the numbers
    res = results()
    refs = references()
    allrows = [r for nb in res for r in res[nb]]
    labels_seen = collections.Counter(r["result"] for r in allrows)
    label_order = [label for label, _, _ in LABELS] + sorted(set(labels_seen) - {label for label, _, _ in LABELS})
    sections = findings_sections()
    run = run_described()

    story = []
    story += [Spacer(1, 30 * mm), p("simplydata reliability report", title_style),
              p("Every comparison of the simplydata_tests notebooks, with its reference, tolerance and label",
                centre), Spacer(1, 8 * mm), p(f"<i>{escape(MOTTO)}</i>", centre), Spacer(1, 3 * mm),
              p(f"<i>{escape(SECOND_LINE)}</i>", centre), Spacer(1, 14 * mm)]
    facts = [["comparisons", f"{len(allrows):,}"],
             ["notebooks", f"{len(res)}"],
             ["agree within the tolerance", f"{labels_seen.get('agrees', 0):,}"],
             ["explained differences", f"{sum(n for k, n in labels_seen.items() if k.startswith('explained')):,}"],
             ["differences not yet explained", f"{labels_seen.get('difference not yet explained', 0):,}"],
             ["simplydata errors being fixed", f"{labels_seen.get('simplydata error – being fixed', 0):,}"],
             ["skipped or waiting", f"{labels_seen.get('skipped', 0) + labels_seen.get('waiting for word', 0):,}"]]
    story.append(grid([[c(a), c(b)] for a, b in facts], header=False, widths=[70 * mm, 35 * mm], numeric_from=1))
    story += [Spacer(1, 10 * mm),
              p("Results: " + escape(run or "the results files of this repository") + ".", centre),
              p("Author: Manish Prakash Victor (ORCID 0000-0002-1409-5894). Licence: MIT, as the repository.",
                centre)]
    story.append(p("Built by scripts/make_pdf_report.py (workflow pdf-report) from the files of "
                   "bioforkgithub/simplydata_tests; the same files give the same PDF.", centre))
    story.append(PageBreak())

    toc = TableOfContents()
    toc.levelStyles = [ParagraphStyle("toc1", parent=body, fontName="Sans-Bold", fontSize=9.5, leftIndent=6 * mm,
                                      firstLineIndent=-6 * mm),
                       ParagraphStyle("toc2", parent=body, fontSize=8.5, leftIndent=12 * mm, firstLineIndent=-6 * mm)]
    story += [p("Contents", ParagraphStyle("contents", parent=h1)), toc, PageBreak()]

    # ---- 1. what was compared and how
    story.append(p("1. What was compared, and how", h1))
    story.append(p(
        "simplydata gives plain-language words for statistics in Python and computes nothing itself: each word "
        "calls an established package (pandas, NumPy, SciPy, statsmodels, matplotlib, seaborn, or R through rpy2). "
        "This study asks whether each word makes that call correctly. Every notebook runs simplydata words, and only "
        "simplydata words, on published data and compares each number they give with a reference value that was "
        "made beforehand and stored in <font name=\"Mono\">references/</font>. A comparison has an explicit "
        "tolerance, and its outcome is written with a label (section 2), never dropped."))
    story.append(p(
        "Where a published answer exists it is the reference: values certified by the US National Institute of "
        "Standards and Technology (the Statistical Reference Datasets and the NIST/SEMATECH e-Handbook), or numbers "
        "printed in the paper, book, manual or help page that first analysed the data. Where nothing is printed, the "
        "reference is a run of the field's standard program, usually R and its packages, stored once with the "
        "versions it ran with (section 8). A number that a printed or certified value covers is not compared again "
        "with a program."))
    story.append(p(
        "Certified values are compared by the log relative error, the number of significant digits that agree "
        "(9 for linear problems, 6 for the NIST nonlinear problems, the level used to judge statistical software). "
        "Printed values are compared within their printed rounding. Program runs are compared with a relative "
        "tolerance, mostly 1e-9 to 1e-12 for calculations without iteration, wider where both sides iterate to a "
        "stated stopping rule, and four Monte Carlo standard errors for random answers that do not share a seed."))
    kinds = collections.Counter(r.get("kind", "") for r in allrows)
    story.append(p("Kinds of reference", h3))
    story.append(grid([[c("kind"), c("what it is"), c("comparisons")]] +
                      [[c(k or "(none)"), c(KINDS.get(k, "")), c(f"{n:,}")] for k, n in kinds.most_common()],
                      widths=[28 * mm, 112 * mm, 22 * mm], numeric_from=2))
    sources = collections.Counter(r.get("reference source", "") for r in allrows)
    story.append(p("References most compared with", h3))
    top = sources.most_common(60)
    story.append(grid([[c("reference"), c("comparisons")]] + [[c(s or "(none)"), c(f"{n:,}")] for s, n in top],
                      widths=[140 * mm, 22 * mm], numeric_from=1))
    if len(sources) > len(top):
        rest = sum(n for _, n in sources.most_common()[len(top):])
        story.append(p(f"and {len(sources) - len(top)} other references with {rest:,} comparisons between them "
                       "(each row of every results file names its own).", small))

    # ---- 2. labels
    story.append(p("2. What the labels mean", h1))
    story.append(p("Each comparison carries one label in the column \"result\" of its results file; the column "
                   "\"status\" keeps the machine's PASS, FAIL or SKIP beside it, and the columns \"finding\" and "
                   "\"explanation\" name the entry of FINDINGS.md that explains a difference."))
    story.append(grid([[c("label"), c("meaning")]] + [[c(label), c(meaning)] for label, _, meaning in LABELS],
                      widths=[50 * mm, 112 * mm]))
    story.append(p("FINDINGS.md sorts the differences into three classes: A, a fault in simplydata; B, a limit "
                   "shared with the package underneath or with double-precision numbers; C, a reference that "
                   "differs from the correct answer, not simplydata.", body))

    # ---- 3. totals
    story.append(p("3. Totals", h1))
    shown = [label for label in label_order if labels_seen.get(label) or label in
             ("agrees", "simplydata error – being fixed", "difference not yet explained")]

    def counts_row(name, rows):
        n = collections.Counter(r["result"] for r in rows)
        return [c(name), c(f"{len(rows):,}")] + [c(f"{n.get(label, 0):,}") for label in shown]

    head = [c("notebooks"), c("comparisons")] + [c(SHORT.get(label, label)) for label in shown]
    body_rows = []
    for name, lo, hi in GROUPS:
        rows = [r for nb in res if lo <= int(nb[:2]) <= hi for r in res[nb]]
        if rows:
            body_rows.append(counts_row(name, rows))
    body_rows.append(counts_row("all", allrows))
    w = (162 - 22 - 18) / len(shown)
    story.append(grid([head] + body_rows, widths=[22 * mm, 18 * mm] + [w * mm] * len(shown), numeric_from=1))
    story.append(p("Column names are short forms of the labels of section 2.", small))
    by_kind = [[c(k or "(none)"), c(f"{len(rs):,}")] + [c(f"{collections.Counter(r['result'] for r in rs).get(l, 0):,}")
                                                          for l in shown]
               for k, rs in sorted(((k, [r for r in allrows if r.get("kind", "") == k]) for k in kinds),
                                   key=lambda x: -len(x[1]))]
    story.append(p("By kind of reference", h3))
    story.append(grid([[c("kind"), c("comparisons")] + [c(SHORT.get(l, l)) for l in shown]] + by_kind,
                      widths=[22 * mm, 18 * mm] + [w * mm] * len(shown), numeric_from=1))

    # ---- 4. per notebook
    story.append(PageBreak())
    story.append(p("4. Each notebook", h1))
    story.append(p("One table per notebook: each simplydata word it compares (the first word of the column "
                   "\"method\"), how many comparisons, of which kinds of reference, and their labels. Every single "
                   "comparison, with both numbers, is in results/&lt;notebook&gt;.csv."))
    for nb, rows in res.items():
        story.append(p(f"{escape(nb)}", h2))
        datasets = sorted({r.get("dataset", "") for r in rows if r.get("dataset")})
        n = collections.Counter(r["result"] for r in rows)
        story.append(p(f"{len(rows):,} comparisons on {len(datasets)} datasets; "
                       + "; ".join(f"{SHORT.get(l, l)} {n[l]:,}" for l in label_order if n.get(l)) + ".", small))
        words = collections.OrderedDict()
        for r in rows:
            words.setdefault(word_of(r["method"]), []).append(r)
        table = [[h("word"), h("compared"), h("certified"), h("printed"), h("program")] +
                 [h(SHORT.get(l, l)) for l in shown]]
        for word, rs in words.items():
            k = collections.Counter(r.get("kind", "") for r in rs)
            lab = collections.Counter(r["result"] for r in rs)
            table.append([c(word), c(len(rs)), c(k.get("certified", 0)), c(k.get("printed", 0)),
                          c(k.get("run reference", 0))] + [c(lab.get(l, 0)) for l in shown])
        wl = (166 - 40 - 14 - 13 * 3) / len(shown)
        story.append(grid(table, widths=[40 * mm, 14 * mm, 13 * mm, 13 * mm, 13 * mm] + [wl * mm] * len(shown),
                          numeric_from=1))

    # ---- 5. explained differences by finding
    story.append(PageBreak())
    story.append(p("5. Every explained difference, by finding", h1))
    story.append(p("Each finding below has the evidence FINDINGS.md gives for it, then every comparison it explains "
                   "with both numbers. \"difference\" and \"score\" are as stored: for a certified value the score is "
                   "the number of digits that agree (log relative error), otherwise the relative difference."))
    explained = [r for r in allrows if r["result"].startswith("explained difference")]
    by_finding = collections.defaultdict(list)
    for r in explained:
        by_finding[r.get("finding") or "(none)"].append(r)

    def finding_key(fid):
        m = re.match(r"([A-Z])(\d+)", fid)
        order = {"B": 0, "C": 1, "A": 2}
        return (order.get(m.group(1), 3), int(m.group(2))) if m else (9, 0)

    for fid in sorted(by_finding, key=finding_key):
        rows = by_finding[fid]
        title, pieces = evidence_for(fid, sections)
        labels = collections.Counter(r["result"] for r in rows)
        story.append(p(f"{escape(fid)}. {md(title)}" if title else escape(fid), h2))
        story.append(p(f"{len(rows)} comparisons: " + "; ".join(f"{l} {k}" for l, k in labels.most_common()) + ".",
                       small))
        line = rows[0].get("explanation", "")
        if line and line.strip().rstrip(".") != title.strip().rstrip("."):
            story.append(p(f"<b>In one line:</b> {md(line)}"))
        if pieces:
            story.append(p("Evidence (FINDINGS.md)", h3))
            for piece in pieces:
                story += md_blocks(piece)
        story.append(p("The comparisons", h3))
        table = [[c("notebook"), c("dataset"), c("method"), c("quantity"), c("simplydata"), c("reference"),
                  c("reference source"), c("difference"), c("tolerance")]]
        for r in rows:
            tol = f"{num(r.get('tolerance', ''))} {r.get('tolerance type', '')}".strip()
            table.append([c(r["notebook"][:2]), c(r.get("dataset", "")), c(r["method"]), c(r.get("quantity", "")),
                          c(num(r.get("simplydata", ""))), c(num(r.get("reference", ""))),
                          c(r.get("reference source", "")), c(num(r.get("difference", ""))), c(tol)])
        story.append(grid(table, widths=[mm * x for x in (8, 20, 25, 20, 19, 19, 22, 15, 14)]))

    # ---- 6. coverage of words
    story.append(PageBreak())
    story.append(p("6. Coverage of simplydata's words", h1))
    words = simplydata_words()
    if words is None:
        story.append(p("simplydata is not installed where this report was built, so its words could not be listed. "
                       "Install it from requirements.txt and build again."))
    else:
        keys = {k for _, k in words}
        compared = collections.Counter()
        for r in allrows:
            key = match_word(word_of(r["method"]), keys)
            if key:
                compared[key] += 1
        story.append(p(
            f"simplydata has {len(keys)} words. A word counts as compared here when it is the first word of the "
            f"column \"method\" of at least one comparison: {len(compared)} of the {len(keys)} are. Words used only "
            "inside another word's comparison (for example a distribution word that reads a fit) are not counted. "
            "The words not compared here are checked in simplydata's own reliability study "
            "(validation/report.md in the simplydata repository), or compute nothing (help, lessons, undo)."))
        groups = collections.OrderedDict()
        for group, key in words:
            groups.setdefault(group, []).append(key)
        table = [[c("group of words"), c("words"), c("compared here"), c("comparisons")]]
        for group, ks in groups.items():
            table.append([c(group), c(len(ks)), c(sum(1 for k in ks if compared[k])),
                          c(f"{sum(compared[k] for k in ks):,}")])
        table.append([c("all"), c(len(keys)), c(len(compared)), c(f"{sum(compared.values()):,}")])
        story.append(grid(table, widths=[90 * mm, 20 * mm, 25 * mm, 25 * mm], numeric_from=1))
        story.append(p("Each word and its comparisons", h3))
        table = [[c("group"), c("word"), c("comparisons here")]]
        for group, ks in groups.items():
            for k in ks:
                table.append([c(group), c(k), c(f"{compared[k]:,}" if compared[k] else "–")])
        story.append(grid(table, widths=[70 * mm, 65 * mm, 27 * mm], numeric_from=2))

    # ---- 7. how to rerun
    story.append(PageBreak())
    story.append(p("7. How to run it again", h1))
    story.append(p("Everything runs from a copy of the repository with Python 3.12, R 4.x and the packages below. "
                   "The workflow notebooks.yml does exactly this on GitHub Actions (ubuntu-latest) for every pull "
                   "request and every Monday, and fails if a notebook stops or a comparison falls outside its "
                   "tolerance without an explanation in FINDINGS.md."))
    story.append(Preformatted(
        "git clone https://github.com/bioforkgithub/simplydata_tests\n"
        "cd simplydata_tests\n"
        "python -m pip install -r requirements.txt     # simplydata at the pinned commit, and the references' Python\n"
        "#   R: the packages listed in .github/workflows/notebooks.yml (install.packages(...))\n"
        "python scripts/get_data.py                     # the datasets not stored here, checked against SHA256SUMS\n"
        "python scripts/run_notebooks.py                # every notebook, top to bottom; writes results/*.csv\n"
        "python scripts/references/<notebook>.py        # optional: make a reference file again\n"
        "python -m pip install reportlab==4.4.10\n"
        "python scripts/make_pdf_report.py              # this PDF, from the results\n", mono))
    story.append(p("For the same last digits as the stored results, set OPENBLAS_CORETYPE=Haswell and one thread "
                   "(OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1), as the workflow does: comparisons near a tolerance "
                   "can otherwise move between processors (the label \"machine\")."))

    # ---- 8. versions
    story.append(p("8. Package versions", h1))
    story.append(p("The Python packages are pinned in requirements.txt:"))
    story.append(Preformatted("\n".join(requirements()), mono))
    story.append(p("The programs and packages that made the stored reference values, from the column \"versions\" of "
                   "references/*.csv:"))
    vs = versions(refs)
    story.append(grid([[c("program or package"), c("version"), c("reference values")]] +
                      [[c(n), c(v), c(f"{k:,}")] for n, v, k in vs], widths=[60 * mm, 60 * mm, 30 * mm],
                      numeric_from=2))
    ran_on = collections.Counter(r.get("ran on", "") for r in allrows)
    story.append(p("Where the simplydata side ran (column \"ran on\" of the results files):"))
    story.append(grid([[c("system"), c("comparisons")]] + [[c(k or "(not recorded)"), c(f"{n:,}")]
                                                            for k, n in ran_on.most_common()],
                      widths=[140 * mm, 22 * mm], numeric_from=1))

    out.parent.mkdir(parents=True, exist_ok=True)
    doc = Doc(str(out), pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=16 * mm,
              bottomMargin=16 * mm, title="simplydata reliability report", author="Manish Prakash Victor",
              subject="Every comparison of the simplydata_tests notebooks", creator="scripts/make_pdf_report.py",
              invariant=1)
    doc.seq = 0
    doc.multiBuild(story, onFirstPage=footer, onLaterPages=footer)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(OUT), help="where to write the PDF")
    args = parser.parse_args(argv)
    build(Path(args.out))
    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
