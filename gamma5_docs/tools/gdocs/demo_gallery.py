"""Gallery of the expression viewer: ``site/assets/viewer/demo.html`` and its data ``demo_data.json``.

Every expression of the gallery is real model content: the operators come from the order store of the run
(``site/resources/run_sm_one_family/sm_one_family_store.json``, read with ``hepkit_renorm.store``), the vertex
functions from ``LocalFunctional.vertex_function`` (through ``hepkit_renorm.rules.vertex_rule``), the coefficients
from the run report. Each item is rendered statically and with the interactive pager side by side.

Run (licensed Symbolica in the environment)::

    PYTHONPATH=python/hepkit_renorm:gamma5_docs/tools python -m gdocs.demo_gallery [SITE_DIR]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from . import mathview as mv

SITE = Path(__file__).resolve().parents[2] / "site"


def _store():
    from hepkit_renorm import builtin_model
    from hepkit_renorm.store import OrderStore

    model = builtin_model("sm_one_family")
    text = (SITE / "resources" / "run_sm_one_family" / "sm_one_family_store.json").read_text(encoding="utf-8")
    return model, OrderStore.from_json(model, text)


def _class_density(functional, cls):
    """The part of ``functional`` with the field content ``cls`` (a Density)."""
    from hepkit_renorm.functional import Density

    monos = [m for m in functional.monomials() if m.fields() == tuple(cls)]
    if not monos:
        raise KeyError(f"no monomial with fields {cls}")
    return Density.from_monomials(functional.fields, monos)


def build(site: Path = SITE) -> dict:
    model, store = _store()
    S0 = store.entry("S0").coefficient(0)
    D1 = store.entry("Delta1").coefficient(0)
    dS = store.entry("deltaS1").coefficient(-1)
    DS = store.entry("DeltaS1").coefficient(0)
    report = json.loads((site / "resources" / "run_sm_one_family" / "run_report.json").read_text(encoding="utf-8"))
    res = report["results"]["1"]
    sections: list[dict] = []

    def item(title: str, expr, note: str = "", interactive: bool = True, page_size: int | None = None,
             budget: int = mv.VIEW_BUDGET) -> dict:
        m_static = mv.M(expr)
        m_inter = mv.M(expr, interactive=True, page_size=page_size, budget=budget) if interactive else None
        return {"title": title, "note": mv.prose(note) if note else "", "static": m_static,
                "interactive": m_inter}

    # ---------------------------------------------------------------- tree level
    tree = []
    tree_classes = [
        ("Gluon kinetic term", ("G", "G")),
        ("W kinetic term", ("W", "W")),
        ("Hypercharge kinetic term", ("B", "B")),
        ("Cubic gluon self-interaction", ("G", "G", "G")),
        ("Quartic gluon self-interaction", ("G", "G", "G", "G")),
        ("Cubic W self-interaction", ("W", "W", "W")),
        ("Gauge fixing (gluon)", ("G", "b_G")),
        ("Auxiliary field (gluon)", ("b_G", "b_G")),
        ("Ghost kinetic term (gluon)", ("gh_G", "agh_G")),
        ("Ghost-gluon vertex", ("G", "gh_G", "agh_G")),
        ("Quark-doublet kinetic term", ("Q", "Qbar")),
        ("Lepton-singlet kinetic term", ("E", "Ebar")),
        ("Quark-gluon vertex", ("G", "Q", "Qbar")),
        ("Lepton-doublet W vertex", ("W", "L", "Lbar")),
        ("Hypercharge vertex of the up singlet", ("B", "U", "Ubar")),
        ("Up-type Yukawa term", ("Q", "Ubar", "H")),
        ("Up-type Yukawa term, conjugate", ("Qbar", "U", "Hdag")),
        ("Down-type Yukawa term", ("Q", "Dbar", "Hdag")),
        ("Electron Yukawa term", ("Lbar", "E", "H")),
        ("Higgs kinetic term", ("H", "Hdag")),
        ("Higgs-W vertex", ("W", "H", "Hdag")),
        ("Higgs-W-W vertex", ("W", "W", "H", "Hdag")),
        ("Higgs-W-B vertex", ("W", "B", "H", "Hdag")),
        ("Quartic Higgs coupling", ("H", "H", "Hdag", "Hdag")),
        ("Source of the gluon (ρ)", ("G", "gh_G", "src_G")),
        ("Source of the ghost (ζ)", ("gh_G", "gh_G", "src_gh_G")),
        ("Source of the quark doublet (R̄)", ("gh_W", "Q", "src_Q")),
        ("Source of the quark doublet (R)", ("gh_W", "Qbar", "src_Qbar")),
        ("Source of the Higgs doublet (σ†)", ("gh_W", "H", "src_H")),
        ("Source of the Higgs doublet (σ)", ("gh_B", "Hdag", "src_Hdag")),
    ]
    for title, cls in tree_classes:
        tree.append(item(title, mv.operator(_class_density(S0, cls)),
                         interactive=len(cls) >= 3))
    sections.append({"id": "tree", "title": "Tree-level operators of S_0",
                     "intro": mv.prose("Operators of the regularised classical action [[ S0 ]] (gauge fixing, "
                                       "ghosts and BRST sources included), one leg class each, with their "
                                       "coefficients as stored in the order store. Indices of the fermion "
                                       "lines are written in matrix notation; the sign of reordering "
                                       "Grassmann-odd fields into the displayed order is included."),
                     "items": tree})

    # ---------------------------------------------------------------- breaking
    brk = []
    for title, cls in [("Breaking, gluon ghost and up singlet", ("gh_G", "U", "Ubar")),
                       ("Breaking, W ghost and lepton doublet", ("gh_W", "L", "Lbar")),
                       ("Breaking, hypercharge ghost and electron singlet", ("gh_B", "E", "Ebar"))]:
        brk.append(item(title, mv.operator(_class_density(D1, cls))))
    brk.append(item("Vertex function of the breaking (gluon ghost, quark doublet)",
                    mv.vertex(D1, ("gh_G", "Qbar", "Q"), model=model, momentum_conservation=True),
                    "Legs [[ leg(omega(A), p(1)) ]], [[ leg(Qbar, p(2)) ]], [[ leg(Q, p(3)) ]], all momenta "
                    "incoming; [[ p(3) ]] is eliminated by momentum conservation."))
    sections.append({"id": "breaking", "title": "The breaking Δ_1",
                     "intro": mv.prose("[[ eqs(Delta1, ord(sop, S0)) ]] couples a ghost, a physical fermion and its "
                                       "sterile partner through an evanescent matrix."),
                     "items": brk})

    # ---------------------------------------------------------------- singular counterterms
    sing = []
    for title, cls in [("Singular counterterm, gluon two-point operators", ("G", "G")),
                       ("Singular counterterm, Higgs two-point operators", ("H", "Hdag")),
                       ("Singular counterterm, quark-gluon operators", ("G", "Q", "Qbar")),
                       ("Singular counterterm, up-type Yukawa operator", ("Q", "Ubar", "H")),
                       ("Singular counterterm, quartic Higgs operators", ("H", "H", "Hdag", "Hdag")),
                       ("Singular counterterm, ghost source operator", ("gh_G", "Q", "src_Q"))]:
        sing.append(item(title, mv.operator(_class_density(dS, cls), loop="pole")))
    const = res["constants"]
    for key, title in [("Z[G]", "δZ of the gluon field"), ("Z[Q]", "δZ of the quark doublet"),
                       ("Z[H]", "δZ of the Higgs doublet"), ("d[g_s]", "δg_s"), ("d[lambda]", "δλ"),
                       ("d[m2]", "δm²"), ("d[y_u]", "δy_u")]:
        sing.append(item(f"Renormalisation constant {title}", mv.coefficient(const[key], loop="pole"),
                         "Value from the run report (units [[ lfe ]])."))
    sections.append({"id": "singular", "title": "Singular counterterms δS^(1)",
                     "intro": mv.prose("Minimal subtraction of the one-loop poles; every renormalisation "
                                       "constant includes [[ lfe ]]. Evanescent operators carry hatted "
                                       "derivatives such as [[ pdhat2(H) ]]."),
                     "items": sing})

    # ---------------------------------------------------------------- finite counterterms
    fin = []
    for title, cls in [("Finite counterterm, gluon two-point operator", ("G", "G")),
                       ("Finite counterterm, quark-gluon operators", ("G", "Q", "Qbar")),
                       ("Finite counterterm, W-Higgs operators", ("W", "H", "Hdag")),
                       ("Finite counterterm, ghost source operators", ("gh_W", "L", "src_L"))]:
        fin.append(item(title, mv.operator(_class_density(DS, cls), loop="finite")))
    rest = res["restoring"]
    for key in list(rest)[:6]:
        fin.append(item(f"Restoring coefficient of {key.split('[')[0]} ({key.split('[')[1].rstrip(']')})",
                        mv.coefficient(rest[key], loop="finite"),
                        "Coefficient of the generated basis element in [[ DeltaS1 ]] (run report); the basis "
                        "names of the automation are not the report's names."))
    sections.append({"id": "finite", "title": "Finite restoring counterterms ΔS^(1)",
                     "intro": mv.prose("[[ eqs(DeltaS1, lf*sumover(alpha, rc(alpha)*Qop(alpha))) ]]; all "
                                       "objects are four-dimensional (barred)."),
                     "items": fin})

    # ---------------------------------------------------------------- vertex functions
    vtx = []
    for title, ent, cls, loop in [
            ("Tree vertex function, three gluons", S0, ("G", "G", "G"), None),
            ("Tree vertex function, four W bosons", S0, ("W", "W", "W", "W"), None),
            ("Tree vertex function, quark-gluon", S0, ("G", "Q", "Qbar"), None),
            ("Tree vertex function, up-type Yukawa", S0, ("Q", "Ubar", "H"), None),
            ("Tree vertex function, Higgs-W", S0, ("W", "H", "Hdag"), None),
            ("Tree vertex function, Higgs-W-W", S0, ("W", "W", "H", "Hdag"), None),
            ("Tree vertex function, fermion two-point", S0, ("Q", "Qbar"), None),
            ("Singular counterterm, gluon two-point function", dS, ("G", "G"), "pole"),
            ("Singular counterterm, Higgs two-point function", dS, ("H", "Hdag"), "pole"),
            ("Singular counterterm, quark-gluon vertex", dS, ("G", "Q", "Qbar"), "pole"),
            ("Singular counterterm, Higgs-B vertex", dS, ("B", "H", "Hdag"), "pole"),
            ("Finite counterterm, gluon two-point function", DS, ("G", "G"), "finite"),
            ("Finite counterterm, three-gluon vertex", DS, ("G", "G", "G"), "finite"),
            ("Finite counterterm, W-Higgs vertex", DS, ("W", "H", "Hdag"), "finite"),
            ("Breaking vertex function, W ghost and lepton doublet", D1, ("gh_W", "Lbar", "L"), None)]:
        info: dict = {}
        val = mv.vertex(ent, cls, model=model, info=info, loop=loop)
        leg_text = ", ".join(f"[[ {mv.M(lg, lead_rationals=False)['src']} ]]" for lg in info["legs"])
        vtx.append(item(title, val, "Legs " + leg_text + "; the Feynman rule is [[ 1𝑖 ]] times the vertex "
                                                           "function."))
    sections.append({"id": "vertex", "title": "Vertex functions",
                     "intro": mv.prose("Vertex functions [[ V ]] in the report's convention (Feynman rule "
                                       "[[ 1𝑖*V ]], momenta incoming, no momentum conservation imposed). "
                                       "Fermion legs are written as matrices between [[ Qbar ]] and [[ Q ]]."),
                     "items": vtx})

    # ---------------------------------------------------------------- long expressions
    long = []
    long.append(item("Pole part of δλ, term by term", mv.coefficient(const["d[lambda]"]),
                     "The sum is paginated: hover or click the marker at its end to see the next terms, "
                     "<em>Inspect</em> opens them in the main view and <em>Back to parent</em> returns.",
                     page_size=5, budget=60_000))
    long.append(item("Singular counterterm, gluon two-point function (pager preview)",
                     mv.vertex(dS, ("G", "G"), model=model, loop="pole"),
                     "The automatic product preview of the pager hides the long sum behind a marker.",
                     budget=60_000))
    long.append(item("Singular counterterm, Higgs two-point function",
                     mv.vertex(dS, ("H", "Hdag"), model=model, loop="pole"), budget=60_000))
    long.append(item("Breaking, three gluons and a ghost",
                     mv.vertex(D1, ("G", "G", "G", "gh_G"), model=model), page_size=4, budget=60_000))
    sections.append({"id": "long", "title": "Long expressions in the pager",
                     "intro": mv.prose("Expressions with many terms, explored with a small page size. The "
                                       "payloads of these items are allowed 60 kB."),
                     "items": long})

    # ---------------------------------------------------------------- prose
    prose_items = [
        "The regulator is [[ ep ]] with [[ eqs(d, 4 - 2*ep) ]]; the Levi-Civita symbol is written "
        "[[ levi(mu, nu, rho, sigma) ]] with the curly glyph, so that it is never confused with the regulator.",
        "Every Lorentz index is split, [[ eqs(gam(mu), gammabar(mu) + gammahat(mu)) ]] and "
        "[[ eqs(g(mu, nu), gbar(mu, nu) + ghat(mu, nu)) ]], with [[ eqs(ghat(mu, lo(mu)), d - 4) ]].",
        "The physical field of a multiplet is [[ eqs(psi(f), ord(Pc(f), Psi(f))) ]] and its sterile partner "
        "[[ eqs(chi(f), ord(Pcbar(f), Psi(f))) ]]; the gauge vertex [[ ord(Psibar(f), gammabar(mu), Tg(a, f), "
        "Pc(f), Psi(f), A(lo(mu), a)) ]] contains only the physical component.",
        "The coupling-weighted generators obey [[ eqs(comm(Tg(a), Tg(b)), 1𝑖*fg(a, b, c)*Tg(c)) ]], "
        "with [[ eqs(fg(I, J, K), g*levi(I, J, K)) ]] for <em>SU</em>(2).",
        "The finite counterterm is [[ eqs(DeltaS1, lf*sumover(alpha, rc(alpha)*Qop(alpha))) ]] and the singular "
        "one carries [[ lfe ]] through the constants [[ dZ(A, G) ]], [[ dZ(Q) ]], [[ dpar(lambda) ]].",
        "Brackets inside a span are allowed: [[ Tr(T(A), acomm(T(C), T(D))) ]]; a literal \\[[ is escaped, and "
        "text such as a &lt; b &amp; c stays text.",
    ]
    sections.append({"id": "prose", "title": "Prose with inline mathematics",
                     "intro": mv.prose("Rich text produced by [[ x ]] spans of <code>mathview.prose</code>."),
                     "items": [{"title": "", "prose": mv.prose(t)} for t in prose_items]})
    count = sum(len(s["items"]) for s in sections)
    return {"title": "Expression viewer gallery", "count": count, "sections": sections}


PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="icon" href="data:,">
<title>Expression viewer gallery</title>
<style>
  :root { color-scheme: light dark; --fg: #1f2933; --bg: #f7f8fa; --card: #ffffff; --line: #dde2e8; --muted: #5b6672;
          --accent: #1f5f8b; }
  @media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --fg: #e6e9ee; --bg: #14171c;
          --card: #1c2027; --line: #2e343d; --muted: #9aa5b1; --accent: #7fb4dc; } }
  :root[data-theme="dark"] { --fg: #e6e9ee; --bg: #14171c; --card: #1c2027; --line: #2e343d; --muted: #9aa5b1;
          --accent: #7fb4dc; }
  :root[data-theme="light"] { color-scheme: light; }
  :root[data-theme="dark"] { color-scheme: dark; }
  body { margin: 0; font: 15px/1.55 system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--fg);
         background: var(--bg); }
  header { padding: 18px 24px 8px; max-width: 1280px; margin: 0 auto; }
  header h1 { font-size: 22px; margin: 0 0 4px; font-weight: 600; }
  header p { margin: 4px 0; color: var(--muted); max-width: 70em; }
  .bar { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; }
  .bar button { font: inherit; color: inherit; background: var(--card); border: 1px solid var(--line);
                border-radius: 6px; padding: 3px 10px; cursor: pointer; }
  main { max-width: 1280px; margin: 0 auto; padding: 0 24px 48px; }
  h2 { font-size: 18px; font-weight: 600; margin: 28px 0 4px; }
  .intro { color: var(--muted); margin: 0 0 10px; max-width: 70em; }
  .item { background: var(--card); border: 1px solid var(--line); border-radius: 8px; padding: 10px 14px 12px;
          margin: 10px 0; --mv-bg: var(--card); --mv-accent: var(--accent); }
  .item h3 { font-size: 14px; font-weight: 600; margin: 0 0 6px; }
  .item .note { font-size: 13px; color: var(--muted); margin: 6px 0 0; }
  .cols { display: grid; grid-template-columns: minmax(0, 1fr) minmax(0, 1fr); gap: 16px; }
  .col { min-width: 0; }
  .col .lab { font-size: 11px; text-transform: uppercase; letter-spacing: .05em; color: var(--muted); }
  .src { font: 11.5px/1.4 ui-monospace, SFMono-Regular, Menlo, monospace; color: var(--muted); margin-top: 4px;
         overflow-wrap: anywhere; }
  .prose { font-size: 15px; }
  @media (max-width: 820px) { .cols { grid-template-columns: minmax(0, 1fr); } }
</style>
</head>
<body>
<header>
  <h1>Expression viewer gallery</h1>
  <p>Expressions of the one-family Standard Model in the BMHV scheme, produced from the automation's order store
  and run report, printed by Symbolica's MathML printer in the site notation (left) and explored with
  Symbolica's pager on a pre-computed payload (right). In the interactive view, hover or click a numbered
  marker to inspect the omitted part; the keys + − 0 zoom when the expression has focus.</p>
  <div class="bar"><span id="count"></span><button type="button" id="theme">Toggle light / dark</button></div>
</header>
<main id="main"></main>
<script type="module">
import { renderMath, ready } from './mathview.js';
const root = document.documentElement;
document.getElementById('theme').onclick = () => {
  const dark = root.dataset.theme ? root.dataset.theme === 'dark' : matchMedia('(prefers-color-scheme: dark)').matches;
  root.dataset.theme = dark ? 'light' : 'dark';
};
const params = new URLSearchParams(location.search);
if (params.get('theme')) root.dataset.theme = params.get('theme');
const data = await (await fetch('./demo_data.json')).json();
await ready();
const main = document.getElementById('main');
document.getElementById('count').textContent = `${data.count} items`;
const el = (tag, cls, html) => { const e = document.createElement(tag); if (cls) e.className = cls;
  if (html !== undefined) e.innerHTML = html; return e; };
for (const sec of data.sections) {
  const h = el('h2'); h.id = sec.id; h.textContent = sec.title; main.append(h);
  main.append(el('p', 'intro', sec.intro));
  for (const it of sec.items) {
    const card = el('section', 'item');
    if (it.title) { const t = el('h3'); t.textContent = it.title; card.append(t); }
    if (it.prose !== undefined) { card.append(el('div', 'prose', it.prose)); main.append(card); continue; }
    const cols = el('div', 'cols');
    const a = el('div', 'col'); a.append(el('div', 'lab', 'Static MathML'));
    const sa = el('div'); a.append(sa); renderMath(sa, it.static);
    const src = el('div', 'src'); src.textContent = it.static.src; a.append(src);
    const b = el('div', 'col');
    const iv = it.interactive && it.interactive.view;
    b.append(el('div', 'lab', iv ? `Interactive pager (${(it.interactive.view.bytes / 1024).toFixed(1)} kB payload` +
      `${it.interactive.view.complete ? '' : ', partial'})` : 'Interactive pager (not pre-computed)'));
    const sb = el('div'); b.append(sb);
    if (iv) renderMath(sb, it.interactive, { interactive: true }); else renderMath(sb, it.static);
    cols.append(a, b); card.append(cols);
    if (it.note) card.append(el('div', 'note', it.note));
    main.append(card);
  }
}
document.body.dataset.ready = '1';
</script>
</body>
</html>
"""


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    site = Path(argv[0]) if argv else SITE
    mv.write_assets(site)
    data = build(site)
    d = site / "assets" / "viewer"
    (d / "demo_data.json").write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    (d / "demo.html").write_text(PAGE, encoding="utf-8")
    print(f"gallery: {data['count']} items, {(d / 'demo_data.json').stat().st_size} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
