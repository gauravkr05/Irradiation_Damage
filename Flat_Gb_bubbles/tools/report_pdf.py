"""report_pdf.py - render REPORT.md (with its figures) to REPORT.pdf.

Handles the subset of Markdown used in REPORT.md: '#'/'##' headings, paragraphs, pipe tables,
'1.' and '-' lists, images, **bold** and `code`.  Needs reportlab (pip install reportlab) and the
DejaVu Sans fonts (for Δ, θ, √, ≈, ...).
usage: python3 tools/report_pdf.py [REPORT.md] [REPORT.pdf]
"""
import os, re, sys
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (Image, KeepTogether, ListFlowable, ListItem, Paragraph, SimpleDocTemplate,
                                Spacer, Table, TableStyle, CondPageBreak)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "REPORT.md")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "REPORT.pdf")
BASE = os.path.dirname(os.path.abspath(SRC))

FD = "/usr/share/fonts/TTF/"
pdfmetrics.registerFont(TTFont("DV", FD + "DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DV-B", FD + "DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DV-I", FD + "DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFont(TTFont("DV-BI", FD + "DejaVuSans-BoldOblique.ttf"))
pdfmetrics.registerFont(TTFont("DVM", FD + "DejaVuSansMono.ttf"))
pdfmetrics.registerFontFamily("DV", normal="DV", bold="DV-B", italic="DV-I", boldItalic="DV-BI")

INK, MUTED, RULE, HEAD_BG = colors.HexColor("#1b1f24"), colors.HexColor("#5b6470"), colors.HexColor("#c9ced6"), colors.HexColor("#eef1f5")
ST = {
    "title": ParagraphStyle("title", fontName="DV-B", fontSize=19, leading=24, textColor=INK, spaceAfter=4),
    "byline": ParagraphStyle("byline", fontName="DV", fontSize=9, leading=12, textColor=MUTED, spaceAfter=14),
    "h2": ParagraphStyle("h2", fontName="DV-B", fontSize=13, leading=17, textColor=INK, spaceBefore=14, spaceAfter=6),
    "body": ParagraphStyle("body", fontName="DV", fontSize=9.3, leading=13.2, textColor=INK, spaceAfter=6),
    "lead": ParagraphStyle("lead", fontName="DV-B", fontSize=9.3, leading=13.2, textColor=INK, spaceAfter=6),
    "sub": ParagraphStyle("sub", fontName="DV-B", fontSize=9.6, leading=13, textColor=INK, spaceBefore=6, spaceAfter=4),
    "cell": ParagraphStyle("cell", fontName="DV", fontSize=7.6, leading=9.6, textColor=INK),
    "cellh": ParagraphStyle("cellh", fontName="DV-B", fontSize=7.6, leading=9.6, textColor=INK),
    "cap": ParagraphStyle("cap", fontName="DV-I", fontSize=8, leading=10.5, textColor=MUTED, alignment=TA_CENTER, spaceAfter=10),
    "li": ParagraphStyle("li", fontName="DV", fontSize=9.3, leading=13.2, textColor=INK),
}
PAGE_W, PAGE_H = A4
MARGIN = 17 * mm
AVAIL_W = PAGE_W - 2 * MARGIN


def inline(s):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    s = re.sub(r"`([^`]+)`", r'<font name="DVM" size="0.92em">\1</font>', s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    s = s.replace('size="0.92em"', "")
    return s


def table(rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [r for r in cells if not all(re.fullmatch(r":?-{3,}:?", c) for c in r)]
    ncol = max(len(r) for r in cells)
    cells = [r + [""] * (ncol - len(r)) for r in cells]
    pad = 8
    nat, mn = [0.0] * ncol, [0.0] * ncol
    for i, r in enumerate(cells):
        fn = "DV-B" if i == 0 else "DV"
        for j, c in enumerate(r):
            txt = re.sub(r"[*`]", "", c)
            nat[j] = max(nat[j], pdfmetrics.stringWidth(txt, fn, 7.6) + pad)
            mn[j] = max(mn[j], max((pdfmetrics.stringWidth(w, fn, 7.6) for w in txt.split()), default=0) + pad)
    if sum(nat) <= AVAIL_W:
        extra = AVAIL_W - sum(nat)
        w = [n + extra * n / sum(nat) for n in nat]
    else:
        spare = AVAIL_W - sum(mn)
        flex = [n - m for n, m in zip(nat, mn)]
        w = [m + (spare * f / sum(flex) if sum(flex) > 0 else 0) for m, f in zip(mn, flex)] if spare > 0 else [AVAIL_W / ncol] * ncol
    data = [[Paragraph(inline(c), ST["cellh" if i == 0 else "cell"]) for c in r] for i, r in enumerate(cells)]
    t = Table(data, colWidths=w, repeatRows=1, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), HEAD_BG),
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 0), (-1, 0), 0.6, INK),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, INK),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.6), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
        ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
    ]))
    return [t, Spacer(1, 8)]


def image(alt, path):
    p = os.path.join(BASE, path)
    iw, ih = ImageReader(p).getSize()
    w = AVAIL_W; h = w * ih / iw
    hmax = PAGE_H - 2 * MARGIN - 60
    if h > hmax:
        h = hmax; w = h * iw / ih
    return [KeepTogether([Spacer(1, 4), Image(p, width=w, height=h), Spacer(1, 3), Paragraph(inline(alt), ST["cap"])])]


def lst(items, ordered):
    li = [ListItem(Paragraph(inline(x), ST["li"]), leftIndent=14) for x in items]
    return [ListFlowable(li, bulletType="1" if ordered else "bullet", start="1" if ordered else "\u2022", bulletFormat="%s." if ordered else None, leftIndent=14, bulletFontName="DV",
                         bulletFontSize=9 if ordered else 6.5, bulletOffsetY=0 if ordered else -2, spaceBefore=0),
            Spacer(1, 6)]


def build():
    lines = open(SRC, encoding="utf-8").read().splitlines()
    out, i, n = [], 0, len(lines)
    hold = []  # flowables that must stay on the same page as the next block (sub-headings, "Conclusions:")

    class _Story(list):
        def append(self, f):
            self.__iadd__([f])

        def __iadd__(self, fl):
            fl = list(fl)
            if hold:
                fl = [KeepTogether(hold + fl)]; hold.clear()
            elif fl and isinstance(fl[0], Table) and len(fl[0]._cellvalues) <= 16:
                fl = [KeepTogether(fl)]
            out.extend(fl)
            return self

    story = _Story()
    title_seen = False
    while i < n:
        ln = lines[i]
        s = ln.strip()
        if not s:
            i += 1; continue
        if s.startswith("# "):
            story.append(Paragraph(inline(s[2:]), ST["title"])); title_seen = True
            i += 1
            while i < n and not lines[i].strip(): i += 1
            if i < n and not lines[i].startswith("#"):
                story.append(Paragraph(inline(lines[i].strip()), ST["byline"])); i += 1
            continue
        if s.startswith("## "):
            story += [CondPageBreak(170), Paragraph(inline(s[3:]), ST["h2"])]
            i += 1
            # the first sentence after a heading is the section's answer: set it in bold
            while i < n and not lines[i].strip(): i += 1
            if i < n and lines[i].strip() and not lines[i].lstrip().startswith(("|", "!", "-", "1.", "#", "**")):
                story.append(Paragraph(inline(lines[i].strip()), ST["lead"])); i += 1
            continue
        if s.startswith("|"):
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(lines[i]); i += 1
            story += table(rows); continue
        m = re.fullmatch(r"!\[([^\]]*)\]\(([^)]+)\)", s)
        if m:
            story += image(m.group(1), m.group(2)); i += 1; continue
        if re.match(r"^(\d+\.|-)\s", s):
            ordered = s[0].isdigit(); items = []
            while i < n and re.match(r"^(\d+\.|-)\s", lines[i].strip()):
                items.append(re.sub(r"^(\d+\.|-)\s+", "", lines[i].strip())); i += 1
            story += lst(items, ordered); continue
        if re.fullmatch(r"\*\*[^*]+\*\*", s):
            hold.append(Paragraph(inline(s), ST["sub"])); i += 1; continue
        buf = []
        while i < n and lines[i].strip() and not re.match(r"^(#|\||!\[|\d+\.\s|-\s)", lines[i].strip()):
            buf.append(lines[i].strip()); i += 1
        p = Paragraph(inline(" ".join(buf)), ST["body"])
        if " ".join(buf).endswith(":"):
            hold.append(p)
        else:
            story.append(p)
    return out


def footer(c, doc):
    c.saveState()
    c.setFont("DV", 7.5); c.setFillColor(MUTED)
    c.drawString(MARGIN, 10 * mm, "Flat_Gb_bubbles Results Report")
    c.drawRightString(PAGE_W - MARGIN, 10 * mm, f"page {doc.page}")
    c.restoreState()


doc = SimpleDocTemplate(OUT, pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=MARGIN, bottomMargin=MARGIN + 3 * mm,
                        title="Flat_Gb_bubbles Results Report", author="Gaurav Kumar")
doc.build(build(), onFirstPage=footer, onLaterPages=footer)
print("wrote", OUT)
