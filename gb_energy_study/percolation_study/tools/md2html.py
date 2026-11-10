"""Minimal Markdown -> HTML (headings, tables, lists, bold, code) for REPORT.md, figures appended."""
import re, html, sys, os
src = open(sys.argv[1]).read().splitlines(); out = []; i = 0
def inl(t):
    t = html.escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t); t = re.sub(r"\*(.+?)\*", r"<i>\1</i>", t)
    return re.sub(r"`(.+?)`", r"<code>\1</code>", t)
while i < len(src):
    l = src[i]
    if l.startswith("```"):
        j = i + 1; blk = []
        while not src[j].startswith("```"): blk.append(html.escape(src[j])); j += 1
        out.append("<pre>" + "\n".join(blk) + "</pre>"); i = j + 1; continue
    if l.startswith("|"):
        rows = []
        while i < len(src) and src[i].startswith("|"):
            if not re.match(r"^\|[-| ]+\|$", src[i]): rows.append([c.strip() for c in src[i].strip("|").split("|")])
            i += 1
        out.append("<table>" + "".join("<tr>" + "".join(f"<{'th' if k == 0 else 'td'}>{inl(c)}</{'th' if k == 0 else 'td'}>" for c in r) + "</tr>" for k, r in enumerate(rows)) + "</table>"); continue
    m = re.match(r"^(#+) (.*)", l)
    if m: out.append(f"<h{len(m[1])}>{inl(m[2])}</h{len(m[1])}>")
    elif re.match(r"^ *[*\-] ", l):
        out.append("<ul>")
        while i < len(src) and re.match(r"^ *[*\-] ", src[i]):
            ind = len(src[i]) - len(src[i].lstrip())
            out.append(f"<li style='margin-left:{ind*0.8}em'>{inl(src[i].lstrip()[2:])}</li>"); i += 1
        out.append("</ul>"); continue
    elif re.match(r"^\d+\. ", l): out.append(f"<p>{inl(l)}</p>")
    elif l.strip() == "---": out.append("<hr/>")
    elif l.strip(): out.append(f"<p>{inl(l)}</p>")
    i += 1
figs = [("fig1_count_and_coverage.png", "Fig. 1. Bubble count vs coverage (top) and coverage vs time (bottom). Open marker = first merger, filled = single bubble left."),
        ("fig2_events_vs_angle.png", "Fig. 2. Coverage (top) and time (bottom) at first merger, A1/Atot > 0.5 and single bubble, vs measured contact angle."),
        ("fig4_late_stage.png", "Fig. 3. Late stage: smallest gap between islands, coverage and porosity vs time. Dotted = single-bubble time."),
        ("fig3_final_maps.png", "Fig. 4. Column thickness at the single-bubble time (θ = 104°, 120°) or the last map (θ = 128°, 148°; t = 150).")]
out.append("<h2>Figures</h2>")
for f, cap in figs:
    import base64; b64 = base64.b64encode(open('results/' + f, 'rb').read()).decode()
    out.append(f"<p><img src='data:image/png;base64,{b64}' width='660'/></p><p class='cap'>{html.escape(cap)}</p>")
css = "body{font-family:'DejaVu Sans',sans-serif;font-size:9.5pt} table{border-collapse:collapse;margin:6px 0} td,th{border:1px solid #999;padding:2px 5px;font-size:8.5pt} th{background:#eee} h1{font-size:15pt} h2{font-size:12.5pt;margin-top:14px} h3{font-size:11pt} code,pre{font-family:'DejaVu Sans Mono',monospace;font-size:8pt} .cap{font-size:8.5pt;color:#444}"
open(sys.argv[2], "w").write(f"<html><head><meta charset='utf-8'><title>Percolation vs contact angle</title><style>{css}</style></head><body>{''.join(out)}</body></html>")
