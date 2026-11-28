#!/usr/bin/env python
"""Assemble and compile the report.

Generates report/sections/04_classes.typ (class subsections that exist are included, missing ones get a
stub heading so that cross-references resolve) and report/one_loop_gamma5.typ, then compiles to
report/one_loop_gamma5.pdf (or to the path given with --out).

usage: build_report.py [--out PATH] [--strict] [--pages DIR]
  --strict    : missing pieces are an error (default: they are replaced by stub headings).
  --pages DIR : also render every page to DIR/pg-NNN.png (100 dpi) and write DIR/INFO.txt with the page
                count and the first page of every labelled heading.
"""
import json, os, subprocess, sys

REPORT = os.path.dirname(os.path.abspath(__file__))
SEC = os.path.join(REPORT, "sections")

G0 = ["g0_VV", "g0_VVV", "g0_VVVV", "g0_FF", "g0_FFV", "g0_SS", "g0_SSV_SSVV", "g0_SSSS", "g0_FFS", "g0_ghost_sources"]
G1 = ["g1_cV", "g1_cVV", "g1_cVVV_cVVVV", "g1_cFF", "g1_cFFV", "g1_cSS_cSSV", "g1_cFFS", "g1_ccRF"]

HEAD = r'''#import "macros.typ": *

#set document(title: "One-loop renormalisation of the Standard Model in the unbroken phase with a non-anticommuting gamma5", author: "gamma5 working notes")
#set page(paper: "a4", margin: (x: 2.2cm, y: 2.4cm), numbering: "1")
#set text(font: "New Computer Modern", size: 10.5pt, lang: "en")
#show math.equation: set text(font: "New Computer Modern Math")
#set par(justify: true, leading: 0.62em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show heading.where(level: 1): it => { pagebreak(weak: true); v(0.4em); it; v(0.6em) }
#show heading.where(level: 2): it => { v(0.5em); it; v(0.3em) }
#show link: set text(fill: rgb("#1f4e79"))
#show ref: set text(fill: rgb("#1f4e79"))
#show figure.caption: set text(size: 9.5pt)
#set table(inset: 5pt)

#align(center)[
  #text(size: 17pt, weight: "bold")[One-loop renormalisation of the Standard Model \ in the unbroken phase with a non-anticommuting $g5$]
  #v(0.5em)
  #text(size: 11.5pt)[A diagrammatic BMHV calculation with Slavnov--Taylor-restoring counterterms]
  #v(0.8em)
  #text(size: 10pt)[Working report, #datetime(year: 2026, month: 11, day: 27).display("[day] [month repr:long] [year]")]
]
#v(2.5em)
#align(center)[#block(width: 86%)[
  #set par(justify: true)
  #set align(left)
  #text(weight: "bold")[Abstract.] We carry out the complete one-loop renormalisation of the Standard Model with one fermion family in the unbroken phase, in dimensional regularisation with the Breitenlohner--Maison--'t Hooft--Veltman treatment of $g5$, in the Feynman gauge and the $overline("MS")$ scheme. Every chiral fermion receives a sterile partner, so that the only breaking of gauge invariance by the regularisation is the evanescent part of the fermion kinetic term. The calculation is purely diagrammatic. The pole parts of eighteen classes of one-particle-irreducible Green's functions give the divergent counterterms, including their evanescent parts. The one-loop graphs with one insertion of the breaking vertex give the finite counterterms that restore the Slavnov--Taylor identity. All counterterms, their Feynman rules, and the renormalised action in the unbroken and in the broken phase are given explicitly. Each step is stated as an operation on Feynman graphs and as the one-loop case of a construction defined at any loop order, built on the $Rop$-operation. Every class was computed twice independently, the results agree with the literature where it exists, and the restored identity was checked against loop functions with complete finite parts. The last section lists what the present tools lack to automate the recipe.
]]
#v(1em)

#outline(depth: 2, indent: auto)
'''


def exists(rel):
    return os.path.isfile(os.path.join(SEC, rel))


def main():
    out = os.path.join(REPORT, "one_loop_gamma5.pdf")
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    strict = "--strict" in sys.argv
    missing = []

    # ---- section 4 container -------------------------------------------------------------------
    lines = ['#import "../macros.typ": *', "", "= Explicit one-loop calculations <sec:classes>", ""]
    if exists("04a_classes_intro.typ"):
        lines.append('#include "04a_classes_intro.typ"')
    for title, label, keys in (
        ("Ghost number zero: divergent counterterms", "sec:g0", G0),
        ("Ghost number one: insertions of the breaking", "sec:g1", G1),
    ):
        lines += ["", f"== {title} <{label}>", ""]
        intro = f"04b_{label.split(':')[1]}_intro.typ"
        if exists(intro):
            lines.append(f'#include "{intro}"')
        for k in keys:
            if exists(f"classes/{k}.typ"):
                lines.append(f'#include "classes/{k}.typ"')
            else:
                missing.append(f"classes/{k}.typ")
                lines.append(f"=== (pending: {k.replace('_', ' ')}) <sec:{k}>")
    if exists("04z_st_solution.typ"):
        lines += ["", '#include "04z_st_solution.typ"']
    else:
        missing.append("04z_st_solution.typ")
        lines += ["", "== Solution of the Slavnov--Taylor conditions <sec:st-solution>", "", "(pending)"]
    with open(os.path.join(SEC, "04_classes.typ"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # ---- main file -----------------------------------------------------------------------------
    body = [HEAD]
    for rel, stub in (
        ("01_introduction.typ", "= Introduction <sec:intro>"),
        ("02_method.typ", "= Systematics of the method <sec:method>"),
        ("03_setup.typ", "= Setup <sec:setup>"),
        ("03b_operator_bases.typ", "== Counterterm operators <sec:bases>"),
        ("04_classes.typ", None),
        ("05_conclusion.typ", "= Conclusion <sec:conclusion>\n\n(pending)"),
        ("06_automation.typ", "= Towards automation <sec:automation>\n\n(pending)"),
        ("07_supporting_material.typ", "= Supporting material <sec:supporting>\n\n(pending)"),
    ):
        if exists(rel):
            body.append(f'#include "sections/{rel}"')
        else:
            missing.append(rel)
            body.append(stub)
    # the broken-phase subsection is included by the conclusion; stub its label if the conclusion is missing
    if not exists("05_conclusion.typ"):
        if exists("05b_broken_phase.typ"):
            body.append('#include "sections/05b_broken_phase.typ"')
        else:
            body.append("== Subtleties in going to the broken phase <sec:broken>\n\n(pending)")
    body.append('\n#bibliography("refs.bib", style: "american-physics-society", title: "References")')
    body.append('#context [#metadata(query(heading).map(h => (h.location().page(), h.level, if h.has("label") { str(h.label) } else { "" }))) <toc-data>]')
    with open(os.path.join(REPORT, "one_loop_gamma5.typ"), "w") as f:
        f.write("\n".join(body) + "\n")

    if missing:
        print("MISSING (stubbed):", ", ".join(missing))
        if strict:
            sys.exit(2)
    env = {k: v for k, v in os.environ.items() if k not in ("LD_LIBRARY_PATH", "LD_PRELOAD")}
    r = subprocess.run(["nix", "shell", "nixpkgs#typst", "-c", "typst", "compile", "--root", REPORT,
                        os.path.join(REPORT, "one_loop_gamma5.typ"), out], env=env, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    if r.returncode != 0:
        print("TYPST COMPILATION FAILED")
        sys.exit(1)
    p = subprocess.run(["nix", "shell", "nixpkgs#poppler-utils", "-c", "pdfinfo", out], env=env, capture_output=True, text=True)
    pages = [l for l in p.stdout.splitlines() if l.startswith("Pages")]
    print("OK:", out, pages[0] if pages else "")

    if "--pages" in sys.argv:
        d = sys.argv[sys.argv.index("--pages") + 1]
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.startswith("pg-") and f.endswith(".png"):
                os.remove(os.path.join(d, f))
        subprocess.run(["nix", "shell", "nixpkgs#poppler-utils", "-c", "pdftoppm", "-r", "100", "-png", out, os.path.join(d, "pg")], env=env, check=True)
        q = subprocess.run(["nix", "shell", "nixpkgs#typst", "-c", "typst", "query", "--root", REPORT,
                            os.path.join(REPORT, "one_loop_gamma5.typ"), "<toc-data>", "--field", "value", "--one"],
                           env=env, capture_output=True, text=True)
        with open(os.path.join(d, "INFO.txt"), "w") as f:
            f.write(f"PDF: {out}\n{pages[0] if pages else ''}\nPage images: pg-NNN.png (NNN = PDF page number, zero-padded as produced by pdftoppm)\n\n")
            f.write("First page of every labelled heading (level, label):\n")
            try:
                for page, level, label in json.loads(q.stdout):
                    if label:
                        f.write(f"  page {page:4d}  level {level}  <{label}>\n")
            except Exception as e:  # noqa
                f.write(f"  (heading query failed: {e}; {q.stderr[:300]})\n")
        print("pages rendered to", d)


if __name__ == "__main__":
    main()
