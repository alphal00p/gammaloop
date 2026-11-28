"""Green's-function classes of the one-loop calculation and their graphs (``data/classes.json``, ``data/graphs/``).

Formats: sections 5 and 6 of ``DATA_FORMAT.md``. Everything here is derived from the automation package
``hepkit_renorm`` and the run outputs under ``site/resources/run_sm_one_family``; nothing about the graphs is typed in.

Classes
-------
A class is a leg class of the package (:func:`hepkit_renorm.classes.leg_classes`) of ghost number 0 (one-loop graphs
of ``S_0``) or ghost number 1 (one-loop graphs with exactly one vertex of the breaking ``Delta_1 = s S_0``). These are
the classes the run evaluates (``run_report.json``: 85 + 168 classes, 2873 graphs). Their fields are the package's
field ids (``G W B``, ``gh_G`` ghost, ``Q Qbar U Ubar D Dbar L Lbar E Ebar`` Dirac fields including the sterile
partner, ``H Hdag``, BRST sources ``src_<field>``) in the canonical field order of the package.

Class ids (:func:`class_id`)
----------------------------
``g<ghost number>:<field>,<field>,...`` with the fields in the package's canonical order, e.g. ``g0:G,G``,
``g0:gh_G,Q,src_Q``, ``g1:G,gh_G,Q,Qbar``. A ghost-number-1 class always carries one ``Delta_1`` insertion, so the
prefix ``g1`` encodes it; :func:`class_id` checks that the insertion argument agrees. Graph ids are
``<class id>-<n>`` with ``n`` the 1-based position (three digits) in the package's deterministic enumeration order.

Groups
------
Every class is mapped to one of the report's eighteen classes of Green's functions (keys ``g0_VV`` ...
``g1_ccRF``, report tab:classes) by the package's kind signature (:func:`hepkit_renorm.classes.kind_signature`):
classes with a BRST source go to ``g0_ghost_sources`` or ``g1_ccRF``. The ghost-number-1 one-point function of the
abelian ghost (``g1:gh_B``, 5 graphs) has no subsection in the report; it forms the extra group ``g1_c``.

Graphs
------
:meth:`hepkit_renorm.class_results.ClassComputer.enumerate` gives the graphs (canonical keys, symmetry factors), the
graph records give the Wick sign (:meth:`hepkit_renorm.graphs.GraphRecord.grassmann_sign`). The counts agree with
``graphs_per_class`` of the run report class by class (checked by :func:`build_graphs`).

Drawing
-------
Each multiplet-level graph is written as HepKIT physics DOT on a *display model*: a HepKIT model whose particles are
the package fields (antiparticle = conjugate field, ``agh_X`` for the ghost ``gh_X``) and whose vertex rules are the
vertex classes of the families (``S0:<fields>``, ``Delta1:<fields>``; DOT attribute ``int_id``). HepKIT parses it
(:meth:`FeynmanDiagram.from_dot`) and draws it with linnet (:meth:`FeynmanDiagram.render`, static snapshot with
estimated label metrics, no browser; two ImPrEd passes). The generic linnet ``Graph.from_dot`` route and per-edge drawing selectors are
not available in the installed wheel, so the HepKIT renderer serves every graph and :func:`clean_svg` rewrites its
output into the site's visual language, using the geometry linnet computed (edge paths, node positions, label
anchors) and the element metadata of the interactive SVG to identify every edge and vertex:

* geometry from linnet: vertex positions, the base curve of every line and the label anchors; the wavy and curly
  decorations are drawn along linnet's base curves (much lighter than linnet's sampled decorations), and a
  self-loop (tadpole), which linnet draws very small, is redrawn as a circle of radius 9 in linnet's direction;
* gauge bosons ``W``, ``B`` wavy, gluon ``G`` curly; fermions solid with arrow; scalars dashed
  with arrow; ghosts dotted with arrow; BRST sources as double-line external stubs; ordinary vertices dots; the
  ``Delta_1`` vertex a crossed circle. Fermion lines are propagators of the full Dirac field (physical chirality and
  sterile partner); the sterile partner is not a separate line.
* Arrows follow the flow of the field: a line ``<X Xbar>`` points into the vertex where it ends on ``X`` (``Q``,
  ``gh_X``, ``H``); an external ``Q`` leg points into the graph, an external ``Qbar`` leg out of it.
* Every line carries the symbol of its field (MathML in a ``foreignObject``, font ``STIX Two Math``); internal lines
  are labelled by the field ``X`` of the propagator ``<X Xbar>``.
* No ids, no scripts, no fonts, no external references; strokes and fills use ``currentColor``; the background of
  the crossed circle and the inner line of a source stub use ``var(--g5-bg, Canvas)`` (the page sets ``--g5-bg``
  to its background colour); coordinates rounded to 0.1.

Contributions
-------------
``contributes`` of a graph comes from ``gdocs/data/graph_contributions.json`` (:func:`compute_contributions`, one
exact evaluation of every graph with the package's one-loop back-end; ``python -m gdocs.graphs --contributions``,
about ten minutes on 32 cores): ``pole`` for ghost number 0 (``N0 != 0``, ``D = 6``), ``breaking`` for ghost number
1 (``P N1 != 0``, ``D = 4``, physical spinors); the other entries are ``None`` (not determined per graph). Without
the file, ``contributes`` is ``None``. The flags are per graph: the graphs of a class can cancel in the sum (e.g. the
ghost-number-one classes with ``gh_G`` in ``g1_cFFS`` have graphs with ``P N1 != 0``, while the report finds the
breaking of that class for ``gh_W`` only); the class results are those of the report and the run, not these flags.

Command line: ``python -m gdocs.graphs [--site DIR] [-j N] [--cache DIR] [--gallery] [--contributions]``.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

__all__ = ["GROUPS", "build_graphs", "class_id", "clean_svg", "field_mathml", "field_text", "graph_file", "legend"]

RENDER_VERSION = "g5g-1"
_SOURCE_HASH = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]   # cache key: drawing code version
HERE = Path(__file__).resolve().parent
DEFAULT_SITE = HERE.parent.parent / "site"
RUN_DIR = "resources/run_sm_one_family"
SVG_NS = "http://www.w3.org/2000/svg"
XHTML_NS = "http://www.w3.org/1999/xhtml"

# ---------------------------------------------------------------------------------------------- ids


def class_id(fields, ghost_number: int, insertion: str | None = None) -> str:
    """Stable id of a class (module docstring): ``g<ghost number>:<fields in canonical order>``.

    ``fields`` are package field ids in any order; they are sorted into the canonical order of the package's field
    table. ``insertion`` must be ``"Delta1"`` exactly for ghost number 1 (or ``None``, which is then implied)."""
    if ghost_number not in (0, 1):
        raise ValueError("classes of the one-loop calculation have ghost number 0 or 1")
    if insertion not in (None, "Delta1"):
        raise ValueError(f"unknown insertion {insertion!r}")
    if insertion == "Delta1" and ghost_number != 1:
        raise ValueError("a Delta1 insertion means ghost number 1")
    order = _field_order()
    try:
        ordered = sorted(fields, key=order.index)
    except ValueError as exc:
        raise ValueError(f"unknown field in {list(fields)}") from exc
    return f"g{ghost_number}:" + ",".join(ordered)


# canonical field order of the package's one-family Standard Model (FieldTable order); checked against the model
# by build_graphs, so that class_id works without loading the package
FIELD_ORDER = (
    "G", "gh_G", "agh_G", "b_G", "W", "gh_W", "agh_W", "b_W", "B", "gh_B", "agh_B", "b_B",
    "Q", "Qbar", "U", "Ubar", "D", "Dbar", "L", "Lbar", "E", "Ebar", "H", "Hdag",
    "src_G", "src_gh_G", "src_W", "src_gh_W", "src_B", "src_gh_B", "src_Q", "src_Qbar", "src_U", "src_Ubar",
    "src_D", "src_Dbar", "src_L", "src_Lbar", "src_E", "src_Ebar", "src_H", "src_Hdag",
)


def _field_order() -> list[str]:
    return list(FIELD_ORDER)


def graph_file(cid: str) -> str:
    """Path (relative to ``site/``) of the graph file of a class."""
    safe = cid.replace(":", "_").replace(",", "_")
    return f"data/graphs/{safe}.json"


# ---------------------------------------------------------------------------------------------- notation

_GAUGE = ("G", "W", "B")


def _mi(s: str) -> str:
    return f"<mi>{s}</mi>"


def _bar(inner: str) -> str:
    return f'<mover accent="true">{inner}<mo>&#xAF;</mo></mover>'


def _sub(base: str, sub: str) -> str:
    return f'<msub>{base}<mi mathvariant="normal">{sub}</mi></msub>'


def field_mathml(f: str) -> str:
    """Presentation MathML (content of ``<math>``) of the symbol of a package field (report notation:
    ``gh = omega``, ``agh = omegabar``, sources ``rho, zeta, Rbar_f / R_f, sigma^dagger / sigma``)."""
    if f in _GAUGE:
        return _mi(f)
    if f.startswith("gh_"):
        return _sub(_mi("&#x3C9;"), f[3:])
    if f.startswith("agh_"):
        return _sub(_bar(_mi("&#x3C9;")), f[4:])
    if f.startswith("src_gh_"):
        return _sub(_mi("&#x3B6;"), f[7:])
    if f.startswith("src_"):
        x = f[4:]
        if x in _GAUGE:
            return _sub(_mi("&#x3C1;"), x)
        if x == "H":
            return f"<msup>{_mi('&#x3C3;')}<mo>&#x2020;</mo></msup>"
        if x == "Hdag":
            return _mi("&#x3C3;")
        if x.endswith("bar"):
            return _sub(_mi("R"), x[:-3])
        return _sub(_bar(_mi("R")), x)
    if f == "H":
        return _mi("H")
    if f == "Hdag":
        return f"<msup>{_mi('H')}<mo>&#x2020;</mo></msup>"
    if f.endswith("bar"):
        return _bar(_mi(f[:-3]))
    return _mi(f)


def field_text(f: str) -> str:
    """Plain-text form of the symbol of a field (Unicode)."""
    if f in _GAUGE:
        return f
    if f.startswith("gh_"):
        return "ω_" + f[3:]
    if f.startswith("agh_"):
        return "ω̄_" + f[4:]
    if f.startswith("src_gh_"):
        return "ζ_" + f[7:]
    if f.startswith("src_"):
        x = f[4:]
        if x in _GAUGE:
            return "ρ_" + x
        if x == "H":
            return "σ†"
        if x == "Hdag":
            return "σ"
        if x.endswith("bar"):
            return "R_" + x[:-3]
        return "R̄_" + x
    if f == "Hdag":
        return "H†"
    if f.endswith("bar"):
        return f[:-3] + "̄"
    return f


def _math(inner: str, block: bool = False) -> str:
    disp = ' display="block"' if block else ""
    return f'<math xmlns="http://www.w3.org/1998/Math/MathML"{disp}>{inner}</math>'


def _math_fields(fields, sep: str = "<mo>,</mo>") -> str:
    """One ``<math>`` with the symbols of ``fields`` separated by ``sep``."""
    return _math("<mrow>" + sep.join(field_mathml(f) for f in fields) + "</mrow>")


def field_kind(f: str) -> str:
    """Line kind of a field: ``gluon``, ``gauge``, ``fermion``, ``scalar``, ``ghost``, ``source``."""
    if f.startswith("src_"):
        return "source"
    if f == "G":
        return "gluon"
    if f in ("W", "B"):
        return "gauge"
    if f.startswith(("gh_", "agh_")):
        return "ghost"
    if f in ("H", "Hdag"):
        return "scalar"
    return "fermion"


# ---------------------------------------------------------------------------------------------- groups

# (key, title, ghost number, kind signatures, report label, report page, report diagram count)
GROUPS = [
    ("g0_VV", "Gauge two-point functions", 0, ["VV"], "Section 4.1.1", 43, 20),
    ("g0_VVV", "Gauge three-point functions", 0, ["VVV"], "Section 4.1.2", 47, 54),
    ("g0_VVVV", "Gauge four-point functions", 0, ["VVVV"], "Section 4.1.3", 53, 237),
    ("g0_FF", "Fermion two-point functions", 0, ["FF"], "Section 4.1.4", 60, 16),
    ("g0_FFV", "Fermion–gauge vertices", 0, ["FFV"], "Section 4.1.5", 63, 46),
    ("g0_SS", "Higgs two-point function", 0, ["SS"], "Section 4.1.6", 67, 8),
    ("g0_SSV_SSVV", "Higgs–gauge functions", 0, ["SSV", "SSVV"], "Section 4.1.7", 71, 146),
    ("g0_SSSS", "Higgs four-point function", 0, ["SSSS"], "Section 4.1.8", 77, 63),
    ("g0_FFS", "Yukawa vertices", 0, ["FFS"], "Section 4.1.9", 83, 32),
    ("g0_ghost_sources", "Ghost and BRST-source functions", 0, ["sources"], "Section 4.1.10", 87, 20),
    ("g1_cV", "Breaking of ωA", 1, ["cV"], "Section 4.2.1", 92, 10),
    ("g1_cVV", "Breaking of ωAA and the anomaly", 1, ["cVV"], "Section 4.2.2", 96, 40),
    ("g1_cVVV_cVVVV", "Breaking of ωAAA and ωAAAA", 1, ["cVVV", "cVVVV"], "Section 4.2.3", 100, 1320),
    ("g1_cFF", "Breaking of ωΨ̄Ψ", 1, ["cFF"], "Section 4.2.4", 105, 42),
    ("g1_cFFV", "Breaking of ωΨ̄ΨA", 1, ["cFFV"], "Section 4.2.5", 109, 220),
    ("g1_cSS_cSSV", "Breaking of ωH†H, ωH†HA and higher", 1, ["cSS", "cSSV", "cSSVV", "cSSSS"],
     "Section 4.2.6", 114, 525),
    ("g1_cFFS", "Breaking of ωΨ̄ΨH", 1, ["cFFS"], "Section 4.2.7", 119, 122),
    ("g1_ccRF", "Breaking with BRST sources", 1, ["sources"], "Section 4.2.8", 124, 54),
    ("g1_c", "Breaking of the one-point function of the abelian ghost", 1, ["c"], "Section 4.2", 92, None),
]

GROUP_NOTES = {
    "g1_c": "Not one of the report's eighteen classes: the run evaluates the ghost-number-one one-point function of "
            "the abelian ghost, with one insertion of the breaking, and the report has no subsection for it.",
}


_ZERO_BY_GROUP = ("graphs that vanish by group theory; the report's enumeration ran on one representative component "
                  "tuple per assignment and did not generate them, and their values vanish")
COUNT_NOTES = {
    "g0_VVV": "The report counts 54 graphs; the run has 3 more in W W W, " + _ZERO_BY_GROUP + ".",
    "g0_VVVV": "The report counts 237 graphs; the run has 9 more in W W W B, " + _ZERO_BY_GROUP + ".",
    "g0_SSV_SSVV": "The report counts 146 graphs; the run has 1 more in W H H\u2020, " + _ZERO_BY_GROUP + ".",
    "g1_cVVV_cVVVV": "The report counts 1320 graphs. It includes the 120 graphs of the external-field sets "
                     "G G W W \u03c9_W, W W W W \u03c9_W and W W \u03c9_W B B, which the automation does not treat as "
                     "classes: no non-zero, Lorentz-invariant and rigidly invariant local monomial has these fields.",
    "g1_ccRF": "The report and the run both count 54 graphs. The report also lists \u03c9_W \u03c9_W \u03c9_W "
               "\u03b6_W (no graph), which is not a class of the run; the run has seven source classes without "
               "graphs that the report does not list.",
}


def _group_of(sig: str, sector: str, ghost: int) -> str:
    key = "sources" if sector == "sources" else sig
    for g in GROUPS:
        if g[2] == ghost and key in g[3]:
            return g[0]
    raise KeyError(f"no group for signature {sig!r} at ghost number {ghost}")


# ---------------------------------------------------------------------------------------------- package side

_STATE: dict = {}


def _load(site_dir: Path):
    """Model, order store and class computer of the run (``site/resources``)."""
    from hepkit_renorm.class_results import ClassComputer
    from hepkit_renorm.spec import builtin_model
    from hepkit_renorm.store import OrderStore

    model = builtin_model("sm_one_family")
    store_path = site_dir / RUN_DIR / "sm_one_family_store.json"
    store = OrderStore.from_json(model, store_path.read_text())
    computer = ClassComputer(model, store, 1)
    if tuple(model.fields) != FIELD_ORDER:
        raise RuntimeError("the package's field order differs from FIELD_ORDER; update FIELD_ORDER")
    return model, store, computer


def _anti(f: str, table) -> str:
    if f.startswith("gh_"):
        return "agh_" + f[3:]
    if f.startswith("agh_"):
        return "gh_" + f[4:]
    return table[f].conjugate


_PARTICLE_SIDE = {"fermion_bar", "scalar_conj"}


def _is_conjugate(f: str, table) -> bool:
    """True for the field that is written as the antiparticle of the display model."""
    fi = table[f]
    if f.startswith("agh_"):
        return True
    if fi.kind in _PARTICLE_SIDE:
        return True
    return fi.kind == "source" and f.endswith(("bar", "dag"))


def display_model(model, families) -> str:
    """JSON of the HepKIT display model (module docstring)."""
    table = model.fields
    spin = {"gauge": 3, "ghost": -1, "fermion": 2, "scalar": 1, "source": 1}
    parts, pdg = [], 1
    for name, fi in table.items():
        if fi.kind in ("antighost", "auxiliary") or _is_conjugate(name, table):
            continue
        conj = _anti(name, table)
        p = {"pdg_code": pdg, "name": name, "antiname": conj, "spin": spin[fi.kind],
             "color": 8 if name == "G" else 1, "mass": "ZERO", "width": "ZERO", "texname": name,
             "antitexname": conj, "charge": 0.0, "ghost_number": 1 if fi.kind == "ghost" else 0,
             "lepton_number": 0, "y_charge": 0, "propagating": fi.kind != "source"}
        parts.append(p)
        if conj != name:
            q = dict(p)
            q.update(pdg_code=-pdg, name=conj, antiname=name, texname=conj, antitexname=name,
                     ghost_number=-p["ghost_number"])
            parts.append(q)
        pdg += 1
    spins = {p["name"]: p["spin"] for p in parts}
    lorentz: dict[tuple, dict] = {}
    rules = []
    for fam in families:
        for cls in fam.rules:
            key = tuple(spins[f] for f in cls)
            if key not in lorentz:
                lorentz[key] = {"name": f"L{len(lorentz)}", "spins": list(key), "structure": "1"}
            rules.append({"name": f"{fam.name}:{','.join(cls)}", "particles": list(cls), "color_structures": ["1"],
                          "lorentz_structures": [lorentz[key]["name"]], "couplings": [["GC0"]]})
    zero = {"name": "ZERO", "lhablock": None, "lhacode": None, "nature": "internal", "parameter_type": "real",
            "value": [0.0, 0.0], "expression": None}
    doc = {"name": "gamma5_docs_display", "restriction": None, "orders": [], "parameters": [zero],
           "particles": parts, "propagators": [], "lorentz_structures": list(lorentz.values()),
           "couplings": [{"name": "GC0", "expression": "1", "orders": []}], "vertex_rules": rules}
    return json.dumps(doc)


_TOPOLOGY = {1: "tadpole", 2: "bubble", 3: "triangle", 4: "box", 5: "pentagon", 6: "hexagon"}


def _cycle(top) -> list[int]:
    """Line indices in cycle order (one loop: every vertex lies on the cycle)."""
    lines = top.lines
    if len(lines) == 1:
        return [0]
    used, order = set(), []
    v = lines[0][0][0]
    start = v
    while True:
        for e, (a, b, _) in enumerate(lines):
            if e in used:
                continue
            if a[0] == v or b[0] == v:
                used.add(e)
                order.append(e)
                v = b[0] if a[0] == v else a[0]
                break
        else:
            break
        if v == start and len(used) == len(lines):
            break
    return order


def _family_label(name: str) -> str:
    return {"S0": _math("<msub><mi>S</mi><mn>0</mn></msub>"),
            "Delta1": _math("<msub><mi mathvariant=\"normal\">&#x394;</mi><mn>1</mn></msub>")}.get(name, name)


def _graph_data(model, families, top, rec, gid: str) -> tuple[dict, dict]:
    """Metadata (section 6) and drawing specification of one graph."""
    table = model.fields
    fams = [families[f].name for f, _ in top.vertices]
    order = _cycle(top)
    loop = [top.lines[e][2][0] for e in order]
    n_lines = len(top.lines)
    topology = _TOPOLOGY.get(n_lines, f"{n_lines}-gon")
    delta = [v for v, n in enumerate(fams) if families[top.vertices[v][0]].insertion]
    sign = rec.grassmann_sign(table)
    sf = top.symmetry_factor
    ext_slots = set(top.externals)
    # caption
    cap = [f"{topology.capitalize()} with internal lines {_math_fields(loop)}."]
    if delta:
        v = delta[0]
        gh = [j for j, f in enumerate(top.vertices[v][1]) if f.startswith("gh_")]
        internal = all((v, j) not in ext_slots for j in gh)
        cap.append("One " + _family_label("Delta1") + " vertex; its ghost leg is "
                   + ("internal." if internal else "external."))
    # the report's sign column: -1 for a closed fermion loop or a closed ghost loop, +1 otherwise (a loop of
    # fermion and ghost lines is not closed: the fermion line ends on external legs)
    closed = None
    if all(table.odd(x) for x in loop):
        if all(x.startswith("gh_") for x in loop):
            closed = "ghost"
        elif not any(x.startswith("gh_") for x in loop):
            closed = "fermion"
    loop_sign = -1 if closed else 1
    cap.append(f"Symmetry factor {sf}; closed-loop sign {'&#x2212;1 (closed ' + closed + ' loop)' if closed else '+1 (no closed fermion or ghost loop)'}, as in the report's tables. "
               f"The package's Wick sign in the leg order of its vertex functions (Dirac conjugates before fermions, antighosts before ghosts) is "
               f"{'&#x2212;1' if sign < 0 else '+1'}; it is a different quantity.")
    meta = {"id": gid, "svg": None, "loop": loop, "topology": topology,
            "insertion": "Delta1" if delta else None, "symmetry_factor": str(sf), "sign": str(loop_sign),
            "wick_sign": str(sign),
            "caption": "<p>" + " ".join(cap) + "</p>", "contributes": None,
            "vertices": [{"family": fams[v], "fields": list(top.vertices[v][1])} for v in range(len(top.vertices))],
            "superficial_degree": str(rec.superficial_degree(table)),
            "name": top.describe(families)}
    # drawing: DOT + edge specification
    dot = ["digraph g {", "  ext [style=invis];"]
    for v, (fi, cls) in enumerate(top.vertices):
        dot.append(f'  v{v} [int_id="{families[fi].name}:{",".join(cls)}"];')
    # sides of the external legs: the legs in the cyclic order of their vertices around the loop, the first half
    # incoming (left), the rest outgoing (right), so that the two sides are contiguous arcs of the loop and no
    # external leg has to cross it; a leg of a field X with antiparticle Xbar is written as an incoming X or an
    # outgoing Xbar (the arrow then still shows the flow of X)
    vpos = {}
    for e in order:
        for end in top.lines[e][:2]:
            vpos.setdefault(end[0], len(vpos))
    for v in range(len(top.vertices)):
        vpos.setdefault(v, len(vpos))
    legs_cyclic = sorted(range(len(top.legs)), key=lambda k: (vpos[top.externals[k][0]], k))
    groups: list[list[int]] = []
    for k in legs_cyclic:
        if groups and top.externals[groups[-1][0]][0] == top.externals[k][0]:
            groups[-1].append(k)
        else:
            groups.append([k])
    n_legs = len(legs_cyclic)
    if len(groups) == 1:
        cut = max(1, n_legs // 2)
        left, right = legs_cyclic[:cut], legs_cyclic[cut:][::-1]
    else:
        # the left side is an arc of whole vertex groups (legs of one vertex stay on one side), as close to half
        # of the legs as possible, preferably containing the first leg of the class
        best = None
        G = len(groups)
        for i in range(G):
            for j in range(1, G):
                arc = [k for g in (groups[(i + t) % G] for t in range(j)) for k in g]
                score = (abs(len(arc) - n_legs / 2), 0 if 0 in arc else 1, i, j)
                if best is None or score < best[0]:
                    best = (score, arc)
        left = best[1]
        right = [k for k in legs_cyclic if k not in left]
        r0 = legs_cyclic.index(left[-1])
        right = sorted(right, key=lambda k: (legs_cyclic.index(k) - r0) % n_legs)[::-1]
    ext_specs = [None] * len(top.legs)
    for k in left + right:
        v = top.externals[k][0]
        f = top.legs[k]
        a = _anti(f, table)
        undirected = ' dir=none,' if a == f else ''
        if k in left:
            dot.append(f'  ext -> v{v} [{undirected} particle="{f}"];')
            side, particle = "in", f
        else:
            dot.append(f'  v{v} -> ext [{undirected} particle="{a}"];')
            side, particle = "out", a
        ext_specs[k] = {"field": f, "kind": field_kind(f), "vertex": f"v{v}", "particle": particle, "side": side,
                        "flip": _is_conjugate(particle, table)}
    int_specs = []
    for (xa, yb, (x, y)) in top.lines:
        if x == y:
            dot.append(f'  v{xa[0]} -> v{yb[0]} [dir=none, particle="{x}"];')
            int_specs.append({"source": f"v{xa[0]}", "sink": f"v{yb[0]}", "particle": x, "field": x,
                              "kind": field_kind(x), "flip": False})
        else:
            dot.append(f'  v{yb[0]} -> v{xa[0]} [particle="{x}"];')
            int_specs.append({"source": f"v{yb[0]}", "sink": f"v{xa[0]}", "particle": x, "field": x,
                              "kind": field_kind(x), "flip": False})
    dot.append("}")
    spec = {"dot": "\n".join(dot), "externals": ext_specs, "internals": int_specs,
            "delta": [f"v{v}" for v in delta]}
    return meta, spec


# ---------------------------------------------------------------------------------------------- SVG geometry

_NUM = r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_TOK = re.compile(rf"[MLCQZmlcqzHVhv]|{_NUM}")


def _parse_path(d: str) -> list[tuple[str, list[tuple[float, float]]]]:
    """Absolute segments ``(cmd, points)`` of an SVG path (M, L, C, Q, Z, H, V; absolute or relative)."""
    toks = _TOK.findall(d)
    out: list[tuple[str, list[tuple[float, float]]]] = []
    i, cmd = 0, None
    cur = (0.0, 0.0)
    start = cur
    nargs = {"M": 2, "L": 2, "C": 6, "Q": 4, "Z": 0, "H": 1, "V": 1}
    while i < len(toks):
        t = toks[i]
        if t.isalpha():
            cmd = t
            i += 1
            if cmd in "Zz":
                out.append(("L", [start]))
                cur = start
                continue
        n = nargs[cmd.upper()]
        vals = [float(x) for x in toks[i:i + n]]
        i += n
        rel = cmd.islower()
        C = cmd.upper()
        if C == "H":
            pts = [((cur[0] if rel else 0) + vals[0], cur[1])]
            C = "L"
        elif C == "V":
            pts = [(cur[0], (cur[1] if rel else 0) + vals[0])]
            C = "L"
        else:
            pts = [(vals[j] + (cur[0] if rel else 0), vals[j + 1] + (cur[1] if rel else 0))
                   for j in range(0, n, 2)]
        if C == "M":
            start = pts[0]
            out.append(("M", pts))
            cmd = "l" if rel else "L"
        else:
            out.append((C, pts))
        cur = pts[-1]
    return out


def _fmt(x: float) -> str:
    s = f"{x:.1f}"
    if s.endswith(".0"):
        s = s[:-2]
    if s in ("-0",):
        s = "0"
    if s.startswith("0."):
        s = s[1:]
    elif s.startswith("-0."):
        s = "-" + s[2:]
    return s


def _join(nums: list[float]) -> str:
    out = ""
    for x in nums:
        t = _fmt(x)
        if out and not t.startswith("-"):
            out += " "
        out += t
    return out


def _path_d(segs) -> str:
    """Compact path data: absolute coordinates rounded to 0.1, written relative."""
    parts = []
    cur = (0.0, 0.0)
    last = None
    for cmd, pts in segs:
        r = [(round(x, 1), round(y, 1)) for x, y in pts]
        if cmd == "M":
            parts.append("M" + _join([r[0][0], r[0][1]]))
            cur = r[0]
            last = "M"
            continue
        c = cmd.lower()
        nums = []
        for x, y in r:
            nums += [round(x - cur[0], 1), round(y - cur[1], 1)]
        parts.append(("" if c == last else c) + (" " if c == last and not _fmt(nums[0]).startswith("-") else "")
                     + _join(nums))
        cur = r[-1]
        last = c
    return "".join(parts)


def _bez(p0, p1, p2, p3, t):
    u = 1 - t
    return (u ** 3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t ** 3 * p3[0],
            u ** 3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t ** 3 * p3[1])


def _polyline(segs, n: int = 16) -> list[tuple[float, float]]:
    pts: list[tuple[float, float]] = []
    cur = None
    for cmd, p in segs:
        if cmd == "M":
            cur = p[0]
            pts.append(cur)
        elif cmd == "L":
            pts.append(p[0])
            cur = p[0]
        elif cmd == "C":
            for k in range(1, n + 1):
                pts.append(_bez(cur, p[0], p[1], p[2], k / n))
            cur = p[2]
        elif cmd == "Q":
            c1 = (cur[0] + 2 / 3 * (p[0][0] - cur[0]), cur[1] + 2 / 3 * (p[0][1] - cur[1]))
            c2 = (p[1][0] + 2 / 3 * (p[0][0] - p[1][0]), p[1][1] + 2 / 3 * (p[0][1] - p[1][1]))
            for k in range(1, n + 1):
                pts.append(_bez(cur, c1, c2, p[1], k / n))
            cur = p[1]
    return pts


def _midpoint(pts, frac: float = 0.5):
    """Point and unit tangent at the given fraction of the arc length of a polyline."""
    lens = [math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    total = sum(lens)
    if total == 0:
        return pts[0], (1.0, 0.0)
    target, acc = frac * total, 0.0
    for i, L in enumerate(lens):
        if acc + L >= target and L > 0:
            t = (target - acc) / L
            a, b = pts[i], pts[i + 1]
            p = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
            return p, ((b[0] - a[0]) / L, (b[1] - a[1]) / L)
        acc += L
    a, b = pts[-2], pts[-1]
    L = math.dist(a, b) or 1.0
    return b, ((b[0] - a[0]) / L, (b[1] - a[1]) / L)


# ---------------------------------------------------------------------------------------------- SVG rewrite

_ROOT_ATTRS = ('fill="none" stroke="currentColor" stroke-width="1" stroke-linecap="round" stroke-linejoin="round"')
_LABEL_FONT = "'STIX Two Math',math,'Cambria Math',serif"
_BG = "var(--g5-bg,Canvas)"


def _detail(el) -> dict:
    d = el.get("data-linnet-detail")
    return json.loads(d) if d else {}


def _translate(s: str | None) -> tuple[float, float]:
    m = re.match(r"translate\(\s*(" + _NUM + r")[ ,]+(" + _NUM + r")\s*\)", s or "")
    return (float(m.group(1)), float(m.group(2))) if m else (0.0, 0.0)


class _Curve:
    """Arc-length parametrisation of a polyline: point, unit tangent and unit normal at arc length ``s``
    (linear extrapolation beyond the ends)."""

    def __init__(self, pts):
        self.pts = pts
        self.cum = [0.0]
        for i in range(len(pts) - 1):
            self.cum.append(self.cum[-1] + math.dist(pts[i], pts[i + 1]))
        self.length = self.cum[-1]

    def at(self, s: float):
        pts, cum = self.pts, self.cum
        if s <= 0 or s >= self.length:
            i = 0 if s <= 0 else len(pts) - 2
        else:
            lo, hi = 0, len(cum) - 1
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if cum[mid] <= s:
                    lo = mid
                else:
                    hi = mid
            i = lo
        a, b = pts[i], pts[i + 1]
        L = cum[i + 1] - cum[i] or 1e-9
        t = (s - cum[i]) / L
        u = ((b[0] - a[0]) / L, (b[1] - a[1]) / L)
        return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])), u, (-u[1], u[0])


def _wave(curve: _Curve, half: float = 3.1, amp: float = 1.7) -> list:
    """Segments of a wavy line along ``curve`` (one cubic per half wave)."""
    n = max(2, round(curve.length / half))
    h = curve.length / n
    k = 4 * amp / 3
    p0 = curve.at(0)[0]
    segs = [("M", [p0])]
    for i in range(n):
        sg = 1 if i % 2 == 0 else -1
        q1, _, n1 = curve.at(i * h + 0.3 * h)
        q2, _, n2 = curve.at(i * h + 0.7 * h)
        p3 = curve.at((i + 1) * h)[0]
        segs.append(("C", [(q1[0] + sg * k * n1[0], q1[1] + sg * k * n1[1]),
                           (q2[0] + sg * k * n2[0], q2[1] + sg * k * n2[1]), p3]))
    return segs


def _coil(curve: _Curve, period: float = 5.6, amp: float = 2.3, back: float = 1.5) -> list:
    """Segments of a curly (gluon) line along ``curve``: the trochoid ``x = c t + B (1 - cos t)``,
    ``y = A sin t`` (loops for ``2 B > c pi``), one cubic per quarter turn."""
    m = max(1, round(curve.length / period))
    c = curve.length / (2 * math.pi * m)
    B = back * c * math.pi / 2

    def local(t):
        x, y = c * t + B * (1 - math.cos(t)), amp * math.sin(t)
        dx, dy = c + B * math.sin(t), amp * math.cos(t)
        q, u, nrm = curve.at(x)
        return ((q[0] + y * nrm[0], q[1] + y * nrm[1]), (dx * u[0] + dy * nrm[0], dx * u[1] + dy * nrm[1]))

    dt = math.pi / 2
    p, d = local(0.0)
    segs = [("M", [p])]
    for j in range(4 * m):
        t1 = (j + 1) * dt
        q, e = local(t1)
        segs.append(("C", [(p[0] + d[0] * dt / 3, p[1] + d[1] * dt / 3), (q[0] - e[0] * dt / 3, q[1] - e[1] * dt / 3),
                           q]))
        p, d = q, e
    return segs


def clean_svg(raw: str, spec: dict, *, font_size: float = 9.0) -> str:
    """Rewrite a linnet/HepKIT SVG into the site's visual language (module docstring). ``spec`` is the drawing
    specification of :func:`_graph_data`: the field and the vertex/particle/side of every external leg and the
    field of every internal line (by its DOT source, sink and particle), and the names of the ``Delta_1``
    vertices. The geometry is linnet's: vertex positions, the base curve of every line (the half-edge paths) and
    the label anchors; wavy and curly decorations are drawn along the base curves."""
    raw = re.sub(r"<style\b[^>]*>.*?</style>", "", raw, flags=re.S)
    raw = re.sub(r"<script\b[^>]*>.*?</script>", "", raw, flags=re.S)
    root = ET.fromstring(raw)
    vb = [float(x) for x in root.get("viewBox").split()]
    q = lambda t: f"{{{SVG_NS}}}{t}"  # noqa: E731

    nodes_by_index: dict[int, tuple[float, float]] = {}
    node_names: dict[int, str] = {}
    halves: dict[int, list] = defaultdict(list)
    labels: dict[int, tuple[float, float]] = {}
    edges: dict[int, dict] = {}
    pending_label = None

    def walk(parent):
        nonlocal pending_label
        for el in list(parent):
            tag = el.tag
            if tag == q("a"):
                kind = el.get("data-linnet-kind")
                det = _detail(el)
                if kind == "halfedge" and el.get("data-linnet-carrier") == "true":
                    e = det["edge"]
                    p = el.find(q("path"))
                    if p is not None:
                        halves[e].append((det, p.get("d")))
                    edges.setdefault(e, det)
                elif kind == "edge":
                    edges[det["edge"]] = det
                    if pending_label is not None:
                        labels[det["edge"]] = pending_label
                        pending_label = None
                elif kind == "node":
                    node_names[det["node"]] = det.get("name", f"v{det['node']}")
                continue
            if tag == q("circle") and el.get("data-linnet-node") is not None:
                nodes_by_index[int(el.get("data-linnet-node"))] = (float(el.get("cx")), float(el.get("cy")))
            elif tag == q("g"):
                fo = el.find(q("foreignObject"))
                if fo is not None:
                    tx, ty = _translate(el.get("transform"))
                    w, h = float(fo.get("width")), float(fo.get("height"))
                    pending_label = (tx + float(fo.get("x", 0)) + w / 2, ty + float(fo.get("y", 0)) + h / 2)
                else:
                    walk(el)

    walk(root)

    # identify the drawn edges with the specification
    int_pool: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for x in spec["internals"]:
        int_pool[(x["source"], x["sink"], x["particle"])].append(x)
    ext_pool: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for x in spec["externals"]:
        ext_pool[(x["vertex"], x["particle"], x["side"])].append(x)
    name_of = lambda i: node_names.get(i, f"v{i}") if i is not None else None  # noqa: E731
    edge_spec: dict[int, dict] = {}
    for e, det in edges.items():
        if det.get("source") is None or det.get("sink") is None:
            # external leg: identified by its vertex, particle and side (the external index of the drawing is
            # not the order of the DOT statements); equal keys carry equal fields
            side = "in" if det.get("source") is None else "out"
            v = name_of(det["sink"] if side == "in" else det["source"])
            pool = ext_pool.get((v, det["particle"], side))
        else:
            key = (name_of(det["source"]), name_of(det["sink"]), det["particle"])
            pool = int_pool.get(key) or int_pool.get((key[1], key[0], key[2]))
        if not pool:
            raise ValueError(f"edge {det} is not in the specification")
        edge_spec[e] = pool.pop()
    if any(int_pool.values()) or any(ext_pool.values()):
        raise ValueError("lines of the specification were not drawn")
    node_pos = {name_of(i): p for i, p in nodes_by_index.items()}

    def oriented(d: str, start) -> list:
        p = _polyline(_parse_path(d))
        return p[::-1] if math.dist(p[0], start) > math.dist(p[-1], start) else p

    layers: dict[str, list[str]] = defaultdict(list)
    arrows: list[str] = []
    extent: list[tuple[float, float]] = list(node_pos.values())
    for e in sorted(edges):
        sp, det = edge_spec[e], edges[e]
        kind = sp["kind"]
        src, snk = name_of(det.get("source")), name_of(det.get("sink"))
        hv = halves.get(e, [])
        if src and snk:
            # internal line: the source half starts at the source vertex, the sink half ends at the sink vertex
            hs = [h for h in hv if h[0].get("flow") == "source"]
            hk = [h for h in hv if h[0].get("flow") == "sink"]
            a = oriented(hs[0][1], node_pos[src])
            b = oriented(hk[0][1], node_pos[snk])[::-1]
            pts = a + b[1:]
            flow = 1                                   # source -> sink
        else:
            v = src or snk
            pts = oriented(hv[0][1], node_pos[v])     # from the vertex outwards
            flow = 1 if src else -1                    # outgoing: away from the vertex
        if sp.get("flip"):
            flow = -flow
        own_plain = True
        if src and snk and src == snk:
            # self-loop (tadpole): linnet draws it small; redraw it as a circle of radius 9 in the direction of
            # linnet's loop, with the same orientation, and move its label outside
            v = node_pos[src]
            cx = sum(x for x, _ in pts) / len(pts) - v[0]
            cy = sum(y for _, y in pts) / len(pts) - v[1]
            L = math.hypot(cx, cy) or 1.0
            d = (cx / L, cy / L)
            area = sum(pts[i][0] * pts[i + 1][1] - pts[i + 1][0] * pts[i][1] for i in range(len(pts) - 1))
            orient = 1 if area >= 0 else -1
            R = 9.0
            c0 = (v[0] + R * d[0], v[1] + R * d[1])
            phi0 = math.atan2(-d[1], -d[0])
            pts = [(c0[0] + R * math.cos(phi0 + orient * 2 * math.pi * k / 48),
                    c0[1] + R * math.sin(phi0 + orient * 2 * math.pi * k / 48)) for k in range(49)]
            if e in labels:
                labels[e] = (v[0] + d[0] * (2 * R + 0.9 * font_size), v[1] + d[1] * (2 * R + 0.9 * font_size))
            own_plain = False
        curve = _Curve(pts)
        extent.extend(pts)
        if kind == "gauge":
            layers["solid"].append(_path_d(_wave(curve) if own_plain else _wave(curve, 3.4, 1.2)))
        elif kind == "gluon":
            layers["solid"].append(_path_d(_coil(curve) if own_plain else _coil(curve, 6.0, 1.8, 1.4)))
        else:
            # linnet's own base curves (one cubic per half edge for an internal line)
            layers[{"fermion": "solid", "scalar": "dash", "ghost": "dot", "source": "src"}[kind]].append(
                "".join(_path_d(_parse_path(h[1])) for h in hv) if own_plain else
                _path_d(_circle_segs(pts)))
        if kind in ("fermion", "scalar", "ghost"):
            if src and snk:
                mid, u, _ = curve.at(curve.length / 2)
            else:
                mid, u, _ = curve.at(curve.length * 0.55)
            if flow < 0:
                u = (-u[0], -u[1])
            arrows.append(_arrow_d(mid, u))
    out: list[str] = []
    if layers["solid"]:
        out.append(f'<path d="{"".join(layers["solid"])}"/>')
    if layers["dash"]:
        out.append(f'<path stroke-dasharray="4.2 2.6" stroke-linecap="butt" d="{"".join(layers["dash"])}"/>')
    if layers["dot"]:
        out.append(f'<path stroke-dasharray=".1 2.4" stroke-width="1.5" d="{"".join(layers["dot"])}"/>')
    if layers["src"]:
        d = "".join(layers["src"])
        out.append(f'<path stroke-width="3" d="{d}"/><path stroke-width="1.1" style="stroke:{_BG}" d="{d}"/>')
    dots, deltas = [], []
    for i, (x, y) in sorted(nodes_by_index.items()):
        if name_of(i) in spec["delta"]:
            deltas.append(_delta_marker(x, y))
        else:
            dots.append(f'<circle cx="{_fmt(x)}" cy="{_fmt(y)}" r="2.2"/>')
    out.append(f'<g fill="currentColor" stroke="none"><path d="{"".join(arrows)}"/>{"".join(dots)}</g>'
               if arrows else f'<g fill="currentColor" stroke="none">{"".join(dots)}</g>')
    out.extend(deltas)
    fs = font_size
    for e, (cx, cy) in sorted(labels.items()):
        f = edge_spec[e]["field"]
        w, h = 6 * fs, 1.6 * fs
        out.append(f'<foreignObject x="{_fmt(cx - w / 2)}" y="{_fmt(cy - h / 2)}" width="{_fmt(w)}" '
                   f'height="{_fmt(h)}">{_math(field_mathml(f), block=True)}</foreignObject>')
    # bounding box: linnet's, enlarged to the drawn lines (redrawn loops, decorations) and the labels
    for e, (cx, cy) in labels.items():
        w = 0.62 * fs * len(field_text(edge_spec[e]["field"]).replace("\u0304", "")) / 2 + 1
        extent += [(cx - w, cy - 0.75 * fs), (cx + w, cy + 0.75 * fs)]
    xs, ys = [x for x, _ in extent], [y for _, y in extent]
    pad = 3
    xa, ya = min(vb[0], min(xs) - pad), min(vb[1], min(ys) - pad)
    xb, yb = max(vb[0] + vb[2], max(xs) + pad), max(vb[1] + vb[3], max(ys) + pad)
    x0, y0, W, H = xa - 2, ya - 2, xb - xa + 4, yb - ya + 4
    head = (f'<svg xmlns="{SVG_NS}" viewBox="{_fmt(x0)} {_fmt(y0)} {_fmt(W)} {_fmt(H)}" width="{_fmt(W)}" '
            f'height="{_fmt(H)}" {_ROOT_ATTRS} style="font:{_fmt(fs)}px/1.2 {_LABEL_FONT};text-align:center;'
            f'white-space:nowrap" role="img">')
    return head + "".join(out) + "</svg>"


def _circle_segs(pts) -> list:
    """A closed polyline as path segments (every second point, smooth enough at radius 9)."""
    return [("M", [pts[0]])] + [("L", [p]) for p in pts[2::2]]


def _delta_marker(x: float, y: float, r: float = 4.4) -> str:
    c = r * 0.7071
    return (f'<circle cx="{_fmt(x)}" cy="{_fmt(y)}" r="{r}" style="fill:{_BG}"/>'
            f'<path d="M{_join([x - c, y - c])}l{_join([2 * c, 2 * c])}m{_join([0, -2 * c])}'
            f'l{_join([-2 * c, 2 * c])}"/>')


def _arrow_d(p, u, size: float = 4.6) -> str:
    """Path data of a filled arrowhead centred at ``p`` pointing along ``u``."""
    nx, ny = -u[1], u[0]
    tip = (p[0] + 0.55 * size * u[0], p[1] + 0.55 * size * u[1])
    back = (p[0] - 0.45 * size * u[0], p[1] - 0.45 * size * u[1])
    w = 0.42 * size
    a = (back[0] + w * nx, back[1] + w * ny)
    b = (back[0] - w * nx, back[1] - w * ny)
    return _path_d([("M", [tip]), ("L", [a]), ("L", [b])]) + "z"


# ---------------------------------------------------------------------------------------------- rendering

_HK: dict = {}


def _render_raw(model_json: str, dot: str) -> str:
    from symbolica.community import graph as lp
    from symbolica.community import hepkit as hep

    if model_json not in _HK:
        _HK.clear()
        _HK[model_json] = hep.Model.from_json(model_json)
    diagram = hep.FeynmanDiagram.from_dot(_HK[model_json], dot)
    # two ImPrEd passes (the second continues from the first), as in the HepKIT DOT-rendering example
    settings = lp.RenderSettings(layouts=lp.LayoutSettings(seed=1, impred_steps=300, impred_spacing=2.4).then(
        impred_steps=300, impred_spacing=2.4))
    return diagram.render(config=settings, style=hep.DiagramStyle(show_particle=True)).snapshot().to_svg()


def _font_size(svg_width: float) -> float:
    """Label size in user units: at least 10 (linnet lays out with 9-unit estimates; the labels are centred on
    linnet's anchors), larger for wide drawings so that labels stay legible when the drawing is shown 220 px wide."""
    return round(min(12.0, max(10.0, svg_width / 220 * 11.5)), 1)


def render_graph(model_json: str, spec: dict, cache: Path | None = None) -> str:
    """Final SVG of one graph (cached by the hash of the render version, the display model and the spec)."""
    key = hashlib.sha256((RENDER_VERSION + "\0" + _SOURCE_HASH + "\0" + model_json + "\0" + json.dumps(spec, sort_keys=True)).encode()
                         ).hexdigest()
    path = cache / key[:2] / f"{key}.svg" if cache else None
    if path is not None and path.exists():
        return path.read_text()
    raw = _render_raw(model_json, spec["dot"])
    vb = [float(x) for x in re.search(r'viewBox="([^"]+)"', raw).group(1).split()]
    svg = clean_svg(raw, spec, font_size=_font_size(vb[2]))
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(f".{os.getpid()}.tmp")
        tmp.write_text(svg)
        tmp.replace(path)
    return svg


# ---------------------------------------------------------------------------------------------- per class (worker)


def _class_job(args):
    fields, cache = args
    model, computer, model_json = _STATE["model"], _STATE["computer"], _STATE["model_json"]
    from hepkit_renorm.classes import kind_signature, sector

    g = computer.ghost_number(fields)
    cid = class_id(fields, g, "Delta1" if g == 1 else None)
    enum, recs = computer.records(fields)
    graphs = []
    for n, (top, rec) in enumerate(zip(enum.topologies, recs), start=1):
        meta, spec = _graph_data(model, computer.families, top, rec, f"{cid}-{n:03d}")
        meta["svg"] = render_graph(model_json, spec, Path(cache) if cache else None)
        graphs.append(meta)
    sig = kind_signature(model, fields)
    return {"id": cid, "fields": list(fields), "ghost_number": g, "signature": sig,
            "group": _group_of(sig, sector(model, fields), g), "graphs": graphs,
            "vanishing": len(enum.vanishing)}


# ---------------------------------------------------------------------------------------------- contributions

CONTRIBUTIONS = HERE / "data" / "graph_contributions.json"


def compute_contributions(site_dir=DEFAULT_SITE, *, processes: int = 16, out: Path = CONTRIBUTIONS,
                          log=print) -> dict:
    """Per-graph flags from one exact evaluation of every graph (one-loop back-end of the package, over GF(p) at
    one random point of the parameters, Feynman gauge, two random samples of momenta and polarisations):

    * ghost number 0: ``pole`` = ``N0 != 0`` in ``D = 6`` with full spinors (as the run's pole step);
      ``breaking`` = ``False`` (no ``Delta_1`` vertex);
    * ghost number 1: ``breaking`` = ``P N1 != 0`` in ``D = 4`` on physical spinors (as the run's breaking step);
      ``pole`` is not evaluated per graph (``None``): the ``D = 6`` evaluation of the 2218 graphs costs about an
      hour on 32 cores.

    A non-zero value proves the flag; a zero value at the random point means zero up to a probability of order
    ``deg / p`` with ``p ~ 2^31``. ``evanescent`` is not determined per graph (``None``). The result is written to
    ``out`` (a build input kept with the code, so that :func:`build_graphs` runs without loop evaluations)."""
    from hepkit_renorm.class_results import ClassComputer
    from hepkit_renorm.classes import leg_classes
    from hepkit_renorm.exact import default_field
    from hepkit_renorm.fingerprint import Fingerprinter

    site_dir = Path(site_dir)
    model, store, _ = _load(site_dir)
    cc = ClassComputer(model, store, 1, processes=processes)
    fixed = {"xi": 1} if "xi" in model.parameters else {}
    fp6 = Fingerprinter(model, default_field(), 6, seed=11, values=fixed)
    fp4 = Fingerprinter(model, default_field(), 4, seed=11, values=fixed)
    doc: dict = {"method": compute_contributions.__doc__.split("\n\n")[0].replace("\n", " "), "classes": {}}
    t0 = time.time()
    for g in (0, 1):
        for c in leg_classes(model, g):
            fields = c.fields
            if not cc.enumerate(fields).topologies:
                doc["classes"][" ".join(fields)] = {}
                continue
            if g == 0:
                res = cc.compute(fields, fp6, cc.samples(fp6, fields, 2, salt=5))
                flags = {r.name: {"pole": bool(r.laurent.N0.any()), "breaking": False} for r in res.graphs}
            else:
                res = cc.compute(fields, fp4, cc.samples(fp4, fields, 2, salt=5, physical=True))
                flags = {r.name: {"pole": None, "breaking": bool(r.laurent.N1.any())} for r in res.graphs}
            doc["classes"][" ".join(fields)] = flags
            log(f"{' '.join(fields)}: {len(flags)} graphs, {time.time() - t0:.0f} s")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(doc, indent=0, sort_keys=True))
    return doc


def _contributions() -> dict:
    if CONTRIBUTIONS.exists():
        return json.loads(CONTRIBUTIONS.read_text())["classes"]
    return {}


# ---------------------------------------------------------------------------------------------- legend


def legend(model_json: str | None = None) -> list[dict]:
    """The legend: ``[{kind, label, svg}]``, drawn with the same code and style as the graphs."""
    def line(inner: str, w: float = 64, h: float = 18) -> str:
        return (f'<svg xmlns="{SVG_NS}" viewBox="0 0 {w} {h}" width="{w}" height="{h}" {_ROOT_ATTRS} '
                f'role="img">{inner}</svg>')

    base = _Curve([(4.0, 9.0), (60.0, 9.0)])
    straight = _path_d([("M", [(4, 9)]), ("L", [(60, 9)])])
    arrow = f'<path fill="currentColor" stroke="none" d="{_arrow_d((32, 9), (1.0, 0.0))}"/>'
    return [
        {"kind": "gauge", "label": "Electroweak gauge field W or B (wavy line)",
         "svg": line(f'<path d="{_path_d(_wave(base))}"/>')},
        {"kind": "gluon", "label": "Gluon G (curly line)", "svg": line(f'<path d="{_path_d(_coil(base))}"/>')},
        {"kind": "fermion", "label": "Fermion: propagator of the Dirac field of a multiplet, physical chirality and "
                                     "sterile partner on one line (solid); the arrow is the fermion flow",
         "svg": line(f'<path d="{straight}"/>{arrow}')},
        {"kind": "scalar", "label": "Higgs doublet H (dashed); the arrow points from H\u2020 to H",
         "svg": line(f'<path stroke-dasharray="4.2 2.6" stroke-linecap="butt" d="{straight}"/>{arrow}')},
        {"kind": "ghost", "label": "Faddeev\u2013Popov ghost \u03c9 (dotted); the arrow points from \u03c9\u0304 "
                                   "to \u03c9",
         "svg": line(f'<path stroke-dasharray=".1 2.4" stroke-width="1.5" d="{straight}"/>{arrow}')},
        {"kind": "source", "label": "BRST source \u03c1, \u03b6, R, R\u0304, \u03c3, \u03c3\u2020 (double line, "
                                    "external legs only)",
         "svg": line(f'<path stroke-width="3" d="{straight}"/><path stroke-width="1.1" style="stroke:{_BG}" '
                     f'd="{straight}"/>')},
        {"kind": "vertex", "label": "Vertex of the tree-level action S\u2080",
         "svg": line(f'<path d="{straight}"/><circle fill="currentColor" stroke="none" cx="32" cy="9" r="2.2"/>')},
        {"kind": "delta1", "label": "Vertex of the breaking \u0394\u2081 = s S\u2080 (crossed circle); one per "
                                    "graph of ghost number one",
         "svg": line(f'<path d="{straight}"/>{_delta_marker(32, 9)}')},
    ]


#: the legend labels as rich text (symbols as MathML through gdocs.mathview.prose), key "html" of each entry
LEGEND_HTML = {
    "gauge": "Electroweak gauge field [[W]] or [[B]] (wavy line)",
    "gluon": "Gluon [[G]] (curly line)",
    "fermion": "Fermion: propagator of the Dirac field of a multiplet, physical chirality and sterile partner on one "
               "line (solid); the arrow is the fermion flow",
    "scalar": "Higgs doublet [[H]] (dashed); the arrow points from [[Hdag]] to [[H]]",
    "ghost": "Faddeev–Popov ghost [[omega]] (dotted); the arrow points from [[omegabar]] to [[omega]]",
    "source": "BRST source [[rho]], [[zsrc]], [[R]], [[Rbar]], [[sigma]], [[sigmadag]] (double line, external legs "
              "only)",
    "vertex": "Vertex of the tree-level action [[S0]]",
    "delta1": "Vertex of the breaking [[eqs(Delta1, ord(sop, S0))]] (crossed circle); one per graph of ghost number one",
}


def legend_html(entries: list[dict]) -> list[dict]:
    """The legend entries with their rich-text label (needs Symbolica)."""
    from gdocs import mathview as mv

    return [{**e, "html": mv.prose(LEGEND_HTML[e["kind"]])} if e["kind"] in LEGEND_HTML else e for e in entries]


# ---------------------------------------------------------------------------------------------- build


def _report_counts(site_dir: Path) -> dict[str, int]:
    doc = json.loads((site_dir / RUN_DIR / "run_report.json").read_text())
    return doc["results"]["1"]["graphs_per_class"]


def build_graphs(site_dir=DEFAULT_SITE, *, processes: int | None = None, cache=None, log=print) -> dict:
    """Write ``data/classes.json`` and ``data/graphs/*.json`` (module docstring); return a summary.

    ``processes``: forked workers (default: CPU count, at most 64); ``cache``: optional directory of rendered
    SVGs (per graph, keyed by content hash)."""
    import multiprocessing as mp

    from hepkit_renorm.classes import leg_classes

    site_dir = Path(site_dir)
    t0 = time.time()
    model, _store, computer = _load(site_dir)
    model_json = display_model(model, computer.families)
    _STATE.update(model=model, computer=computer, model_json=model_json)
    classes = [c.fields for g in (0, 1) for c in leg_classes(model, g)]
    reported = _report_counts(site_dir)
    if processes is None:
        processes = min(64, os.cpu_count() or 1)
    jobs = [(tuple(f), str(cache) if cache else None) for f in classes]
    # largest classes first for load balance
    jobs.sort(key=lambda j: -reported.get(" ".join(j[0]), 0))
    if processes > 1:
        with mp.get_context("fork").Pool(processes) as pool:
            results = pool.map(_class_job, jobs, chunksize=1)
    else:
        results = [_class_job(j) for j in jobs]
    contrib = _contributions()
    for r in results:
        flags = contrib.get(" ".join(r["fields"]))
        for gr in r["graphs"]:
            f = flags.get(gr["name"]) if flags else None
            gr["contributes"] = ({"pole": f["pole"], "evanescent": None, "breaking": f["breaking"]}
                                 if f is not None else None)
    by_fields = {tuple(r["fields"]): r for r in results}
    results = [by_fields[tuple(f)] for f in classes]
    log(f"graphs: {sum(len(r['graphs']) for r in results)} in {len(results)} classes, {time.time() - t0:.1f} s")

    # checks against the run report
    mismatch = [(" ".join(r["fields"]), len(r["graphs"]), reported.get(" ".join(r["fields"])))
                for r in results if reported.get(" ".join(r["fields"])) != len(r["graphs"])]
    if mismatch or len(results) != len(reported):
        raise RuntimeError(f"graph counts differ from the run report: {mismatch}")

    data = site_dir / "data"
    gdir = data / "graphs"
    gdir.mkdir(parents=True, exist_ok=True)
    keep = set()
    for r in results:
        rel = graph_file(r["id"])
        keep.add(Path(rel).name)
        doc = {"class": r["id"], "graphs": r["graphs"]}
        (site_dir / rel).write_text(json.dumps(doc, ensure_ascii=False, separators=(",", ":")))
    for p in gdir.glob("*.json"):
        if p.name not in keep:
            p.unlink()

    groups = []
    group_classes = defaultdict(list)
    for r in results:
        group_classes[r["group"]].append(r)
    for key, title, ghost, _sigs, label, page, report_n in GROUPS:
        members = group_classes.get(key, [])
        if not members:
            continue
        g = {"key": key, "title": title, "ghost_number": ghost, "classes": [m["id"] for m in members],
             "report": [{"label": label, "page": page}], "summary": None,
             "n_graphs": sum(len(m["graphs"]) for m in members), "report_n_graphs": report_n}
        if key in GROUP_NOTES:
            g["note"] = GROUP_NOTES[key]
        if key in COUNT_NOTES:
            g["count_note"] = COUNT_NOTES[key]
        groups.append(g)
    cls_out = []
    for r in results:
        cls_out.append({"id": r["id"], "group": r["group"], "ghost_number": r["ghost_number"], "fields": r["fields"],
                        "label": None, "insertion": "Delta1" if r["ghost_number"] == 1 else None,
                        "n_graphs": len(r["graphs"]), "graphs": graph_file(r["id"]), "result": None,
                        "summary": None, "signature": r["signature"],
                        "text": "Γ(" + " ".join(field_text(f) for f in r["fields"]) + ")",
                        "symbol": _math_fields(r["fields"], "<mspace width=\"0.2em\"/>")})
    doc = {"groups": groups, "classes": cls_out, "legend": legend_html(legend(model_json)),
           "source": "hepkit_renorm enumeration of the run sm_one_family (order 1); drawings by HepKIT / linnet"}
    (data / "classes.json").write_text(json.dumps(doc, ensure_ascii=False, indent=1))
    size = sum(p.stat().st_size for p in gdir.glob("*.json"))
    summary = {"classes": len(cls_out), "graphs": sum(c["n_graphs"] for c in cls_out),
               "ghost_number_0": sum(1 for c in cls_out if c["ghost_number"] == 0),
               "ghost_number_1": sum(1 for c in cls_out if c["ghost_number"] == 1),
               "groups": {g["key"]: g["n_graphs"] for g in groups}, "bytes": size,
               "vanishing": sum(r["vanishing"] for r in results), "seconds": round(time.time() - t0, 1)}
    log(f"classes.json and {len(cls_out)} graph files written; data/graphs = {size / 1e6:.2f} MB")
    return summary


# ---------------------------------------------------------------------------------------------- gallery


def write_gallery(site_dir=DEFAULT_SITE, per_group: int = 5, out: str = "data/graphs_gallery.html") -> Path:
    """A static page with about ``per_group`` graphs of every group, at thumbnail (220 px) and full size, for
    checking the drawings (light and dark: ``?dark`` in the URL or the system preference)."""
    site_dir = Path(site_dir)
    doc = json.loads((site_dir / "data" / "classes.json").read_text())
    cls = {c["id"]: c for c in doc["classes"]}
    rows = []
    for g in doc["groups"]:
        picked = []
        files = [json.loads((site_dir / cls[c]["graphs"]).read_text()) for c in g["classes"]]
        allg = [(f["class"], x) for f in files for x in f["graphs"]]
        if not allg:
            continue
        step = max(1, len(allg) // per_group)
        picked = allg[::step][:per_group]
        rows.append(f'<h2>{html.escape(g["key"])}: {html.escape(g["title"])} ({g["n_graphs"]} graphs)</h2>')
        rows.append('<div class="thumbs">' + "".join(
            f'<figure class="t">{x["svg"]}<figcaption><code>{html.escape(x["id"])}</code></figcaption></figure>'
            for _, x in picked) + "</div>")
        rows.append('<div class="full">' + "".join(
            f'<figure class="f">{x["svg"]}<figcaption><code>{html.escape(x["id"])}</code> {x["topology"]}, '
            f'S = {x["symmetry_factor"]}, sign {x["sign"]}<br>{x["caption"]}</figcaption></figure>'
            for _, x in picked) + "</div>")
    leg = "".join(f'<div class="lg">{e["svg"]}<span>{html.escape(e["label"])}</span></div>' for e in doc["legend"])
    fonts = sorted(site_dir.glob("assets/**/STIXTwoMath*.woff2"))
    face = (f"@font-face{{font-family:'STIX Two Math';src:url('../{fonts[0].relative_to(site_dir).as_posix()}')}}"
            if fonts else "")
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>Graph gallery (check page)</title>
<style>
{face}
:root{{color-scheme:light;--g5-bg:#fff;--fg:#1d1d1f}}
:root.dark{{color-scheme:dark;--g5-bg:#16181c;--fg:#e8e6e3}}
body{{background:var(--g5-bg);color:var(--fg);font:14px/1.4 sans-serif;margin:16px}}
figure{{margin:0;border:1px solid #8884;padding:6px}}
.thumbs,.full{{display:flex;flex-wrap:wrap;gap:8px;align-items:flex-start;margin-bottom:10px}}
.t svg{{width:220px;height:auto;display:block}} .f svg{{display:block}}
figcaption{{font-size:12px;max-width:340px}} .lg{{display:flex;gap:10px;align-items:center}}
.lg svg{{width:96px;height:auto}}
</style>
<script>if(location.search.includes('dark'))document.documentElement.classList.add('dark')</script>
</head><body><h1>Graph gallery (check page, not part of the site)</h1><section>{leg}</section>{''.join(rows)}
<script>for(const s of document.querySelectorAll('.f svg')){{s.style.width=(1.6*parseFloat(s.getAttribute('width')))+'px';s.style.height='auto'}}</script>
</body></html>"""
    path = site_dir / out
    path.write_text(page)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m gdocs.graphs", description=__doc__.split("\n")[0])
    ap.add_argument("--site", default=str(DEFAULT_SITE))
    ap.add_argument("-j", "--processes", type=int, default=None)
    ap.add_argument("--cache", default=None, help="directory of rendered SVGs (optional)")
    ap.add_argument("--gallery", action="store_true", help="also write data/graphs_gallery.html")
    ap.add_argument("--contributions", action="store_true",
                    help="first evaluate every graph once and write gdocs/data/graph_contributions.json (slow)")
    args = ap.parse_args(argv)
    if args.contributions:
        compute_contributions(args.site, processes=args.processes or 16)
    summary = build_graphs(args.site, processes=args.processes, cache=args.cache)
    print(json.dumps(summary, indent=1))
    if args.gallery:
        print(write_gallery(args.site))
    return 0


if __name__ == "__main__":
    sys.exit(main())
