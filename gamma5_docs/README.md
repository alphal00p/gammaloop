# gamma5 documentation site

A static, interactive reference to the one-loop renormalised Lagrangian of the Standard Model with a
non-anticommuting γ₅: BMHV scheme with sterile partners, MS-bar, finite Slavnov–Taylor-restoring counterterms,
Feynman gauge, one fermion family. It presents the renormalised action operator by operator, in the unbroken phase
and, after the shift H → H_v + φ, in the broken phase. For every operator it shows the coefficients at tree level,
in the singular counterterm and in the finite restoring counterterm, the Feynman rules (vertex functions V; the rule
is i V), the one-loop graphs of the Green's functions that determine the coefficients, and a theory card taken from
the report. All mathematics is rendered with the Symbolica expression viewer (MathML), all graphs are drawn with
HepKIT.

**What it is not yet.** The page is a picture of what first-class Lagrangian support in HepKIT and a multi-loop
automation are meant to produce directly. Today it is assembled by the build code in `tools/` from the outputs of the
one-loop automation `python/hepkit_renorm` and from the report. It covers one loop only, the Feynman gauge only and
one fermion family; the broken phase is not computed by the automation (the package has no broken-phase
representation) but derived on the site from the unbroken phase by the report's shift of the Higgs doublet.

## Where the content comes from

| content | source |
|---|---|
| operators, coefficients, vertex functions of the unbroken phase | exact decomposition of the order store of the automation run (`site/resources/run_sm_one_family/`), checked during the build against the report's final numbers (`site/resources/final_counterterms.json`); `tools/gdocs/catalog.py` |
| classes of Green's functions, graphs | enumeration of the automation (graph counts checked against the run report), drawn with HepKIT / linnet; `tools/gdocs/graphs.py` |
| class results (pole part, finite breaking) | the order store restricted to the legs of the class; `tools/gdocs/catalog.py` |
| per-graph contribution flags | one evaluation of every graph with the package's one-loop back-end (`tools/gdocs/data/graph_contributions.json`) |
| broken phase | derived from the unbroken-phase operators by the shift H → H_v + φ and the rotation to the report's mass basis, checked against the automation's order store and the report's broken-phase formulas; `tools/gdocs/broken.py` |
| theory cards, group summaries, page texts, type definitions, conventions | written from the report; `tools/gdocs/theory/*.txt` |
| section and page references | the report's Typst sources and `site/resources/report_headings.txt` |
| equation references (`Eq. (N)` with its page) | `tools/gdocs/theory/report_equations.json`: number and page of every labelled equation of the PDF, from `typst query` on the report (see `gdocs.theory.EQUATION_QUERY`), checked against the PDF text |
| report | `site/report/one_loop_gamma5.pdf` |
| downloads | `site/resources/` (see its `README.md`) |

## Layout

```
DATA_FORMAT.md        contract between the build code and the page (JSON formats, interfaces)
serve.py              loopback static server
tools/build_site.py   one command that rebuilds every data file
tools/gdocs/          build code: mathview (notation, Symbolica MathML), graphs, catalog, broken, theory, names (rich
                      term names)
tools/mock_data.py    development aid: writes placeholder data only where no real file exists
tools/tests/          pytest suites and the contract validator check_site.mjs
site/                 the served site: index.html, assets/, data/ (generated JSON), report/, resources/
```

The generated `site/data` is committed, so the site is served without rebuilding.

## Rebuild

Requirements: a Python (3.11 or newer) with Symbolica 3.0.3 and its HepKIT community modules, `numpy` and
`pytest`; Node 18 or newer for the validator. The Symbolica licence is read from the environment variable
`SYMBOLICA_LICENSE` only; pass it on the command line and never write it into a file:

```
cd gamma5_docs
SYMBOLICA_LICENSE=<key> python tools/build_site.py --date 2026-11-28
```

The script puts `python/hepkit_renorm` and `tools` on the Python path itself. It runs, in order: viewer assets,
classes and graphs, unbroken catalogue with class labels and results, group summaries, theory cards, broken
catalogue, resources, manifest, then the contract validator, and prints a table of the files written. A failing step
stops the build with exit status 1. A full build takes about one minute on 32 cores; `--reuse-graphs` keeps the
existing graph files (about 15 s less).

Options: `--date` (else `SOURCE_DATE` from the environment, else today), `--reuse-graphs`, `--graph-cache DIR`
(rendered drawings), `-j N`, `--contributions` (re-evaluates every graph for the contribution flags, about 7 min on 32
cores), `--sources PROJECT_DIR` (refreshes the copies in `site/resources` from the directory that holds `report/` and
`SCRATCHPAD/`; without it the committed copies are kept and only listed), `--galleries` (rewrites the two check pages
`site/data/graphs_gallery.html` and `site/assets/viewer/demo.html`), `--skip STEP` and `--only STEP`.

## Tests

```
cd gamma5_docs/tools
SYMBOLICA_LICENSE=<key> PYTHONPATH=../../python/hepkit_renorm:. python -m pytest tests     # about 90 s
cd ..
node tools/tests/check_site.mjs --no-mock --strict                                         # contract validator
```

## Serve

```
python3 serve.py --port PORT            # serves site/ on 127.0.0.1:PORT (gzip, no directory listing)
```

The server binds to the loopback interface only. To reach it from elsewhere, put a reverse proxy in front, for
example `tailscale serve --bg http://127.0.0.1:PORT` for the machines of a tailnet, or a Cloudflare quick tunnel,
`cloudflared tunnel --url http://127.0.0.1:PORT`, for a public URL (anyone who has the URL can read the site). The
site needs no network access beyond its own files.

## Notation of the page

Field chips, filters and graph labels write ghosts, antighosts, auxiliary fields and the sources of gauge fields and
ghosts with their gauge factor (or, in the broken phase, the field label) as a subscript: ω_G, ω̄_W, b_B, ρ_G, ζ_W,
ω_Z, ρ_γ. In operators and vertex functions a superscript on these symbols is an adjoint index (ω^a), as in the
report. Every term has a plain `name` (document title, tooltips, search) and a `name_html` with the symbols as
MathML (`tools/gdocs/names.py`). The manifest carries the definitions of the five types of terms (`types`) and the
conventions (`conventions`) as rich text.

## Known limitations

- Term order inside sums and the order of factors are Symbolica's canonical order, not always the report's
  arrangement (e.g. `−4 + d`, Yukawa terms in arbitrary order), and may differ between two builds. Some stored
  expressions keep unevaluated group traces of three or four generators or products such as `P_L P_L`; all are
  correct as written.
- Long sums inside narrow table cells and on phones scroll horizontally inside their box. Symbolica omits printed
  fragments above about 4 kB; the build fails if an omission mark reaches a data file, so a few long expressions are
  split or printed in Symbolica's own product order.
- Graphs are drawn at the level of multiplets, as in the report; the sterile partner is not a separate line, and the
  two terms of the Δ₁ vertex are not distinguished in the drawing. A graph's `sign` is the report's closed-loop sign
  (−1 for a closed fermion or ghost loop); `wick_sign` is the package's Wick sign in the leg order of its vertex
  functions, a different quantity.
- Broken-phase terms show the graphs of their unbroken-phase parents; broken-phase graphs are not drawn (the pole parts
  and the finite breaking are the unbroken ones at H = H_v, report Section 5.5.1).
- The poles of the ghost-number-one classes are shown integrated by parts (no derivative on the ghost) or, on the
  classes ωΨ̄Ψ, as a multiple of Δ₁; the build checks numerically that their four-dimensional part vanishes. The
  three-point gauge classes keep traces of three generators unevaluated, so some of their displayed terms cancel
  only through trace identities.
- Contribution flags are per graph, from one evaluation at a random point modulo a prime; graphs of a class can
  cancel in the sum, so the flags are not class results. The evanescent flag is not determined (null) for every
  graph, the pole flag not for the graphs of ghost number one.
- The class `g1:gh_B` (one-point function of the abelian ghost) has no subsection in the report; its summary rests on
  the run report.
- The theory cards, group summaries and overview formulae are hand transcriptions of the report into the site
  notation, checked by eye and, where the build can, against the automation's numbers.
- The report view embeds the PDF; whether `#page=N` is honoured depends on the browser's PDF viewer.
