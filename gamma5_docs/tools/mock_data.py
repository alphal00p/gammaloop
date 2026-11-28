#!/usr/bin/env python3
"""Mock data set for the gamma5 documentation page.

Writes a complete data set that respects DATA_FORMAT.md to ``site/data``: the manifest, both phase files, the
classes of Green's functions with one graph file per class (placeholder drawings), the theory cards and the list of
resources. It exists so that the page can be built and tested before the real build (``build_site.py``) has run.

Every file written here carries the top-level key ``"mock": true`` (``resources.json`` is a list, so there every
entry carries it). The page shows a "sample data" notice for such files. The coefficients, Feynman rules, graphs and
texts of the mock set are placeholders: they have the shape of the real data, not its content.

A file that exists and does not carry the mock marker was written by the real build and is never overwritten.

    python3 tools/mock_data.py [--site PATH]

No dependency beyond the Python standard library; no Symbolica needed (the mathematics is hand-written MathML).
"""

from __future__ import annotations

import argparse
import html
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SITE = HERE.parent / "site"

# --------------------------------------------------------------------------------------------------------------
# Hand-written MathML
# --------------------------------------------------------------------------------------------------------------

NS = "http://www.w3.org/1998/Math/MathML"
LUNATE = "ϵ"  # regulator epsilon
PI = "π"


def esc(s: str) -> str:
    return html.escape(s, quote=False)


def mi(s: str, normal: bool = False) -> str:
    attr = ' mathvariant="normal"' if normal else ""
    return f"<mi{attr}>{esc(s)}</mi>"


def mn(s) -> str:
    return f"<mn>{esc(str(s))}</mn>"


def mo(s: str) -> str:
    return f"<mo>{esc(s)}</mo>"


def mrow(*xs: str) -> str:
    return "<mrow>" + "".join(xs) + "</mrow>"


def msub(b: str, s: str) -> str:
    return f"<msub>{b}{s}</msub>"


def msup(b: str, s: str) -> str:
    return f"<msup>{b}{s}</msup>"


def msubsup(b: str, s: str, p: str) -> str:
    return f"<msubsup>{b}{s}{p}</msubsup>"


def mfrac(a: str, b: str) -> str:
    return f"<mfrac>{a}{b}</mfrac>"


def bar(x: str) -> str:
    return f'<mover accent="true">{x}<mo stretchy="false">‾</mo></mover>'


def hat(x: str) -> str:
    return f'<mover accent="true">{x}<mo stretchy="false">^</mo></mover>'


def paren(*xs: str) -> str:
    return mrow(mo("("), *xs, mo(")"))


def math(*xs: str, block: bool = False) -> str:
    disp = ' display="block"' if block else ""
    return f'<math xmlns="{NS}"{disp}>' + mrow(*xs) + "</math>"


def M(body: str, text: str, block: bool = False) -> dict:
    return {"html": math(body, block=block), "text": text, "src": None, "view": None}


KAPPA_POLE = mfrac(mn(1), mrow(mn(16), msup(mi(PI), mn(2)), mi(LUNATE)))
KAPPA_FIN = mfrac(mn(1), mrow(mn(16), msup(mi(PI), mn(2))))

# couplings: id -> (MathML, text)
CPL = {
    "gs": (msub(mi("g"), mi("s")), "g_s"),
    "g": (mi("g"), "g"),
    "gp": (msup(mi("g"), mo("′")), "g'"),
    "yu": (msub(mi("y"), mi("u")), "y_u"),
    "yd": (msub(mi("y"), mi("d")), "y_d"),
    "ye": (msub(mi("y"), mi("e")), "y_e"),
    "lam": (mi("λ"), "lambda"),
    "m2": (msup(mi("m"), mn(2)), "m^2"),
    "v": (mi("v"), "v"),
}
GROUP_CPL = {"G": "gs", "W": "g", "B": "gp"}


def monomial(cpl: list[tuple[str, int]]) -> tuple[str, str]:
    parts, text = [], []
    for c, p in cpl:
        m, t = CPL[c]
        parts.append(m if p == 1 else msup(m, mn(p)))
        text.append(t if p == 1 else f"{t}^{p}")
    return "".join(parts), " ".join(text)


def poly(terms: list[tuple[int, int, list[tuple[str, int]]]]) -> tuple[str, str]:
    """Sum of rational multiples of coupling monomials, e.g. [(11, 6, [("gs", 2)])]."""
    out, text = [], []
    for k, (num, den, cpl) in enumerate(terms):
        sign = "-" if num < 0 else "+"
        num = abs(num)
        if k > 0 or sign == "-":
            out.append(mo("−" if sign == "-" else "+"))
            text.append(sign)
        mono, mt = monomial(cpl)
        if den != 1:
            out.append(mfrac(mn(num), mn(den)))
            text.append(f"{num}/{den}")
        elif num != 1 or not cpl:
            out.append(mn(num))
            text.append(str(num))
        out.append(mono)
        if mt:
            text.append(mt)
    return mrow(*out), " ".join(text)


def coef(terms, kappa: str | None) -> dict:
    body, text = poly(terms)
    if kappa == "pole":
        return M(mrow(KAPPA_POLE, paren(body) if len(terms) > 1 else body), f"1/(16 pi^2 eps) ({text})")
    if kappa == "finite":
        return M(mrow(KAPPA_FIN, paren(body) if len(terms) > 1 else body), f"1/(16 pi^2) ({text})")
    return M(body, text)


# --------------------------------------------------------------------------------------------------------------
# Fields
# --------------------------------------------------------------------------------------------------------------

ADJ = {"G": "A", "W": "I", "B": "B"}
GROUP_OF = {"G": "SU(3)", "W": "SU(2)", "B": "U(1)"}
FERMIONS = ["Q", "U", "D", "L", "E"]
F_GROUPS = {"Q": ["G", "W", "B"], "U": ["G", "B"], "D": ["G", "B"], "L": ["W", "B"], "E": ["B"]}
F_NAMES = {"Q": "quark doublet", "U": "up-type quark singlet", "D": "down-type quark singlet",
           "L": "lepton doublet", "E": "charged lepton singlet"}


def gauge_sym(r: str, mu: str = "μ", up: bool = False) -> str:
    if r == "B":
        return msup(mi("B"), mi(mu)) if up else msub(mi("B"), mi(mu))
    return msubsup(mi(r), mi(mu), mi(ADJ[r])) if not up else msup(mi(r), mrow(mi(ADJ[r]), mi(mu)))


def fstrength(r: str, up: bool = False) -> str:
    idx = mrow(mi("μ"), mi("ν"))
    if r == "B":
        return msup(mi("B"), idx) if up else msub(mi("B"), idx)
    return msup(mi(r), mrow(mi(ADJ[r]), idx)) if up else msubsup(mi(r), idx, mi(ADJ[r]))


def psi(f: str) -> str:
    return msub(mi("Ψ"), mi(f))


def psibar(f: str) -> str:
    return msub(bar(mi("Ψ")), mi(f))


def chi(f: str) -> str:
    return msub(mi("χ"), mi(f))


H = mi("H")
HD = msup(mi("H"), mo("†"))


def ghost(r: str, anti: bool = False) -> str:
    base = bar(mi("ω")) if anti else mi("ω")
    return msup(base, mi(ADJ[r]))


def aux(r: str) -> str:
    return msup(mi("b"), mi(ADJ[r]))


def unbroken_fields() -> list[dict]:
    out = []
    for r, name in (("G", "gluon"), ("W", "SU(2) gauge field"), ("B", "hypercharge gauge field")):
        out.append({"id": r, "symbol": M(gauge_sym(r), r), "name": name, "kind": "gauge", "group": GROUP_OF[r],
                    "ghost_number": 0, "statistics": "bose", "note": ""})
    for f in FERMIONS:
        out.append({"id": f, "symbol": M(psi(f), f"Psi_{f}"), "name": F_NAMES[f], "kind": "fermion", "group": None,
                    "ghost_number": 0, "statistics": "fermi", "note": "Dirac spinor of the physical field and its sterile partner."})
    for f in FERMIONS:
        out.append({"id": "chi" + f, "symbol": M(chi(f), f"chi_{f}"), "name": f"sterile partner of {f}",
                    "kind": "sterile", "group": None, "ghost_number": 0, "statistics": "fermi", "note": ""})
    out.append({"id": "H", "symbol": M(H, "H"), "name": "Higgs doublet", "kind": "scalar", "group": None,
                "ghost_number": 0, "statistics": "bose", "note": ""})
    out.append({"id": "Hd", "symbol": M(HD, "H^dagger"), "name": "conjugate Higgs doublet", "kind": "scalar",
                "group": None, "ghost_number": 0, "statistics": "bose"})
    for r in "GWB":
        out.append({"id": "c" + r, "symbol": M(ghost(r), f"omega^{ADJ[r]}"), "name": f"{GROUP_OF[r]} ghost",
                    "kind": "ghost", "group": GROUP_OF[r], "ghost_number": 1, "statistics": "fermi", "note": ""})
    for r in "GWB":
        out.append({"id": "cb" + r, "symbol": M(ghost(r, True), f"omegabar^{ADJ[r]}"), "name": f"{GROUP_OF[r]} antighost",
                    "kind": "antighost", "group": GROUP_OF[r], "ghost_number": -1, "statistics": "fermi", "note": ""})
    for r in "GWB":
        out.append({"id": "b" + r, "symbol": M(aux(r), f"b^{ADJ[r]}"), "name": f"{GROUP_OF[r]} auxiliary field",
                    "kind": "auxiliary", "group": GROUP_OF[r], "ghost_number": 0, "statistics": "bose", "note": None})
    for r in "GWB":
        out.append({"id": "rho" + r, "symbol": M(msubsup(mi("ρ"), mi(ADJ[r]), mi("μ")), f"rho^mu_{ADJ[r]}"),
                    "name": f"source of the {GROUP_OF[r]} gauge field", "kind": "source", "group": GROUP_OF[r],
                    "ghost_number": -1, "statistics": "fermi", "note": ""})
        out.append({"id": "zeta" + r, "symbol": M(msub(mi("ζ"), mi(ADJ[r])), f"zeta_{ADJ[r]}"),
                    "name": f"source of the {GROUP_OF[r]} ghost", "kind": "source", "group": GROUP_OF[r],
                    "ghost_number": -2, "statistics": "bose", "note": ""})
    for f in FERMIONS:
        out.append({"id": "R" + f, "symbol": M(msub(mi("R"), mi(f)), f"R_{f}"), "name": f"source of {f}",
                    "kind": "source", "group": None, "ghost_number": -1, "statistics": "bose", "note": ""})
        out.append({"id": "Rb" + f, "symbol": M(msub(bar(mi("R")), mi(f)), f"Rbar_{f}"), "name": f"source of the conjugate of {f}",
                    "kind": "source", "group": None, "ghost_number": -1, "statistics": "bose", "note": ""})
    out.append({"id": "sigma", "symbol": M(mi("σ"), "sigma"), "name": "source of the Higgs doublet",
                "kind": "source", "group": None, "ghost_number": -1, "statistics": "fermi", "note": ""})
    out.append({"id": "src_Hdag", "symbol": M(msup(mi("σ"), mo("†")), "sigma^dagger"), "name": "source of the conjugate Higgs doublet",
                "kind": "source", "group": None, "ghost_number": -1, "statistics": "fermi", "note": ""})
    for f in out:
        f["id"] = ren(f["id"])
    pos = [i for i, f in enumerate(out) if f["id"] == "E"][0] + 1
    for k, f in enumerate(FERMIONS):
        out.insert(pos + k, {"id": f + "bar", "symbol": M(psibar(f), f"Psibar_{f}"), "name": f"conjugate of the {F_NAMES[f]}", "kind": "fermion",
                             "group": None, "ghost_number": 0, "statistics": "fermi", "note": ""})
    return out


import re


def ren(x: str) -> str:
    """Field ids of the mock follow the convention of the real build (gh_G, src_Q, Qbar, Hdag, ...)."""
    if x == "Hd":
        return "Hdag"
    if x == "sigma":
        return "src_H"
    m = re.fullmatch(r"(cb|c|b|rho|zeta|chi|Rb|R)([GWBQUDLE])", x)
    if not m:
        return x
    pre, s = m.groups()
    return {"cb": f"agh_{s}", "c": f"gh_{s}", "b": f"b_{s}", "rho": f"src_{s}", "zeta": f"src_gh_{s}", "chi": f"chi_{s}",
            "Rb": f"src_{s}bar", "R": f"src_{s}"}[pre]


def ren_fields(fields: list[str]) -> list[str]:
    """Rename and mark the conjugate fermion of a bilinear (first occurrence) as FIELDbar."""
    out = [ren(f) for f in fields]
    ferm = [i for i, f in enumerate(out) if f in FERMIONS]
    sterile = any(f.startswith("chi_") for f in out)
    if len(ferm) >= 2 or (len(ferm) == 1 and sterile):
        out[ferm[0]] += "bar"
    return out


class _RenDict(dict):
    def __getitem__(self, k):
        return dict.__getitem__(self, ren(k))

    def get(self, k, d=None):
        return dict.get(self, ren(k), d)

    def __contains__(self, k):
        return dict.__contains__(self, ren(k))


SYM = _RenDict()  # field id -> MathML (filled below)


def fill_sym(fields):
    for f in fields:
        inner = f["symbol"]["html"]
        start = inner.index(">") + 1
        SYM[f["id"]] = inner[start:-len("</math>")]


# --------------------------------------------------------------------------------------------------------------
# Report references (pages from site/resources/report_headings.txt)
# --------------------------------------------------------------------------------------------------------------

REF = {
    "fields": ("Section 2.2", 21), "action": ("Section 2.3", 23), "brst": ("Section 2.4", 24),
    "delta": ("Section 2.5", 26), "rules": ("Section 2.6", 26), "Qbasis": ("Section 2.8.3", 31),
    "evanescent": ("Section 2.8.5", 34), "normalisation": ("Section 2.8.8", 39),
    "g0_VV": ("Section 3.1.1", 43), "g0_VVV": ("Section 3.1.2", 47), "g0_VVVV": ("Section 3.1.3", 53),
    "g0_FF": ("Section 3.1.4", 60), "g0_FFV": ("Section 3.1.5", 63), "g0_SS": ("Section 3.1.6", 67),
    "g0_SSV_SSVV": ("Section 3.1.7", 71), "g0_SSSS": ("Section 3.1.8", 77), "g0_FFS": ("Section 3.1.9", 83),
    "g0_ghost_sources": ("Section 3.1.10", 87), "g1_cV": ("Section 3.2.1", 92), "g1_cVV": ("Section 3.2.2", 96),
    "g1_cVVV_cVVVV": ("Section 3.2.3", 100), "g1_cFF": ("Section 3.2.4", 105), "g1_cFFV": ("Section 3.2.5", 109),
    "g1_cSS_cSSV": ("Section 3.2.6", 114), "g1_cFFS": ("Section 3.2.7", 119), "g1_ccRF": ("Section 3.2.8", 124),
    "st-solve": ("Section 3.3.4", 135), "st-result": ("Section 3.3.5", 141), "Sren1": ("Section 4.1", 146),
    "ct-rules": ("Section 4.2", 151), "Sren-broken": ("Section 4.4", 156), "bl-tree": ("Section 4.4.3", 158),
    "bl-ST": ("Section 4.4.5", 163), "bp-tadpole": ("Section 4.5.4", 171),
}


def ref(*keys):
    return [{"label": REF[k][0], "page": REF[k][1]} for k in keys]


# --------------------------------------------------------------------------------------------------------------
# Classes of Green's functions and their (placeholder) graphs
# --------------------------------------------------------------------------------------------------------------

GROUPS = [
    ("g0_VV", "Gauge two-point functions", 0), ("g0_VVV", "Gauge three-point functions", 0),
    ("g0_VVVV", "Gauge four-point functions", 0), ("g0_FF", "Fermion two-point functions", 0),
    ("g0_FFV", "Fermion-gauge vertices", 0), ("g0_SS", "Higgs two-point function", 0),
    ("g0_SSV_SSVV", "Higgs-gauge functions", 0), ("g0_SSSS", "Higgs four-point function", 0),
    ("g0_FFS", "Yukawa vertices", 0), ("g0_ghost_sources", "Ghost and BRST-source functions", 0),
    ("g1_cV", "Breaking of the ghost-gauge two-point function", 1), ("g1_cVV", "Breaking of ghost-gauge-gauge functions and the anomaly", 1),
    ("g1_cVVV_cVVVV", "Breaking of ghost functions with three and four gauge fields", 1),
    ("g1_cFF", "Breaking of ghost-fermion-fermion functions", 1), ("g1_cFFV", "Breaking of ghost-fermion-fermion-gauge functions", 1),
    ("g1_cSS_cSSV", "Breaking of ghost-Higgs functions", 1), ("g1_cFFS", "Breaking of ghost-Yukawa functions", 1),
    ("g1_ccRF", "Breaking with BRST sources", 1),
]

CLASSES = [  # (group, fields, n_graphs, loop pool)
    ("g0_VV", ["G", "G"], 9, ["G", "cG", "Q", "U", "D"]), ("g0_VV", ["W", "W"], 14, ["W", "cW", "Q", "L", "H"]),
    ("g0_VV", ["B", "B"], 8, ["Q", "U", "D", "L", "E", "H"]),
    ("g0_VVV", ["G", "G", "G"], 46, ["G", "cG", "Q", "U", "D"]), ("g0_VVV", ["W", "W", "W"], 38, ["W", "cW", "Q", "L", "H"]),
    ("g0_VVVV", ["G", "G", "G", "G"], 212, ["G", "cG", "Q", "U", "D"]), ("g0_VVVV", ["W", "W", "W", "W"], 180, ["W", "cW", "Q", "L", "H"]),
    ("g0_VVVV", ["W", "W", "B", "B"], 64, ["W", "B", "Q", "L", "H"]),
    ("g0_FF", ["Q", "Q"], 6, ["Q", "G", "W", "B", "H", "U", "D"]), ("g0_FF", ["U", "U"], 4, ["U", "G", "B", "Q", "H"]),
    ("g0_FF", ["D", "D"], 4, ["D", "G", "B", "Q", "H"]), ("g0_FF", ["L", "L"], 4, ["L", "W", "B", "E", "H"]),
    ("g0_FF", ["E", "E"], 3, ["E", "B", "L", "H"]),
    ("g0_FFV", ["Q", "Q", "G"], 18, ["Q", "G", "W", "B", "H"]), ("g0_FFV", ["Q", "Q", "W"], 16, ["Q", "G", "W", "B", "H"]),
    ("g0_FFV", ["Q", "Q", "B"], 14, ["Q", "G", "W", "B", "H"]), ("g0_FFV", ["L", "L", "W"], 10, ["L", "W", "B", "H"]),
    ("g0_FFV", ["E", "E", "B"], 7, ["E", "B", "L", "H"]),
    ("g0_SS", ["Hd", "H"], 12, ["H", "W", "B", "Q", "U", "D", "L", "E"]),
    ("g0_SSV_SSVV", ["Hd", "H", "W"], 16, ["H", "W", "B", "Q", "L"]), ("g0_SSV_SSVV", ["Hd", "H", "B"], 10, ["H", "W", "B", "Q", "L", "E"]),
    ("g0_SSV_SSVV", ["Hd", "H", "W", "W"], 96, ["H", "W", "B", "cW", "Q", "L"]),
    ("g0_SSSS", ["Hd", "H", "Hd", "H"], 140, ["H", "W", "B", "Q", "U", "D", "L", "E"]),
    ("g0_FFS", ["Q", "U", "H"], 20, ["Q", "U", "G", "B", "W", "H"]), ("g0_FFS", ["Q", "D", "H"], 20, ["Q", "D", "G", "B", "W", "H"]),
    ("g0_FFS", ["L", "E", "H"], 14, ["L", "E", "B", "W", "H"]),
    ("g0_ghost_sources", ["cbG", "cG"], 3, ["cG", "G"]), ("g0_ghost_sources", ["rhoG", "cG"], 3, ["cG", "G"]),
    ("g0_ghost_sources", ["zetaG", "cG", "cG"], 6, ["cG", "G"]),
    ("g1_cV", ["cG", "G"], 10, ["Q", "chiQ", "U", "chiU", "D", "chiD"]), ("g1_cV", ["cW", "W"], 16, ["Q", "chiQ", "L", "chiL"]),
    ("g1_cV", ["cB", "B"], 12, ["Q", "chiQ", "U", "chiU", "E", "chiE"]),
    ("g1_cVV", ["cG", "G", "G"], 54, ["Q", "chiQ", "U", "chiU", "G"]), ("g1_cVV", ["cW", "W", "B"], 48, ["Q", "chiQ", "L", "chiL", "W", "B"]),
    ("g1_cVVV_cVVVV", ["cG", "G", "G", "G"], 388, ["Q", "chiQ", "U", "chiU", "D", "chiD", "G"]),
    ("g1_cVVV_cVVVV", ["cW", "W", "W", "W", "W"], 640, ["Q", "chiQ", "L", "chiL", "W"]),
    ("g1_cFF", ["cG", "Q", "Q"], 8, ["Q", "chiQ", "G", "W", "B"]), ("g1_cFF", ["cB", "E", "E"], 6, ["E", "chiE", "B"]),
    ("g1_cFFV", ["cG", "Q", "Q", "G"], 42, ["Q", "chiQ", "G", "W", "B"]), ("g1_cFFV", ["cB", "L", "L", "W"], 36, ["L", "chiL", "W", "B"]),
    ("g1_cSS_cSSV", ["cW", "Hd", "H"], 18, ["Q", "chiQ", "U", "chiU", "L", "chiL", "H"]),
    ("g1_cSS_cSSV", ["cB", "Hd", "H", "B"], 44, ["Q", "chiQ", "U", "chiU", "L", "chiL", "H", "B"]),
    ("g1_cFFS", ["cB", "Q", "U", "H"], 52, ["Q", "chiQ", "U", "chiU", "H", "B"]),
    ("g1_ccRF", ["cG", "cG", "RbQ", "Q"], 24, ["Q", "chiQ", "G", "cG"]),
]


def class_id(group: str, fields: list[str]) -> str:
    return f"{group[:2]}:{','.join(fields)}"


TOPO = {2: ["bubble", "tadpole"], 3: ["triangle", "bubble"], 4: ["box", "triangle", "bubble"], 5: ["pentagon", "box", "triangle"]}
KIND_OF: dict[str, str] = {}


def svg_graph(ext: list[str], loop: list[str], insertion: bool, rnd: random.Random) -> str:
    """Placeholder drawing: a one-loop polygon with the external legs attached, lines styled by field kind."""
    n = max(len(loop), 1)
    cx, cy, rad = 100, 70, 30 if n > 1 else 22
    import math as _m
    pts = [(cx + rad * _m.cos(2 * _m.pi * k / n - _m.pi / 2), cy + rad * _m.sin(2 * _m.pi * k / n - _m.pi / 2)) for k in range(n)]

    def style(fid):
        kind = KIND_OF.get(fid, "fermion")
        if kind == "gauge":
            return 'stroke-dasharray="1.5 2.5"'
        if kind in ("ghost", "antighost"):
            return 'stroke-dasharray="1 3" stroke-linecap="round"'
        if kind == "scalar":
            return 'stroke-dasharray="6 3"'
        if kind == "sterile":
            return 'stroke-width="2.6"'
        return ""

    lines = []
    if n == 1:
        lines.append(f'<circle cx="{cx}" cy="{cy - rad}" r="{rad}" fill="none" {style(loop[0])}/>')
        pts = [(cx, cy)]
    else:
        for k in range(n):
            (x1, y1), (x2, y2) = pts[k], pts[(k + 1) % n]
            lines.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {style(loop[k])}/>')
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2
            dx, dy = mx - cx, my - cy
            norm = (dx * dx + dy * dy) ** 0.5 or 1
            tx, ty = cx + dx / norm * (rad - 11), cy + dy / norm * (rad - 11)
            lines.append(f'<text x="{tx:.1f}" y="{ty + 3:.1f}" font-size="8" text-anchor="middle" stroke="none" fill="currentColor">{esc(loop[k])}</text>')
    for j, fid in enumerate(ext):
        vx, vy = pts[j % len(pts)]
        ang = 2 * _m.pi * j / len(ext) + (_m.pi if len(ext) == 2 else -_m.pi / 2)
        ex, ey = cx + 62 * _m.cos(ang), cy + 52 * _m.sin(ang)
        lines.append(f'<line x1="{vx:.1f}" y1="{vy:.1f}" x2="{ex:.1f}" y2="{ey:.1f}" {style(fid)}/>')
        lines.append(f'<text x="{ex + (6 if ex >= cx else -6):.1f}" y="{ey + 3:.1f}" font-size="9" text-anchor="{"start" if ex >= cx else "end"}" stroke="none" fill="currentColor">{esc(fid)}</text>')
    for (x, y) in pts:
        lines.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.8" fill="currentColor" stroke="none"/>')
    if insertion:
        x, y = pts[rnd.randrange(len(pts))]
        lines.append(f'<rect x="{x - 4:.1f}" y="{y - 4:.1f}" width="8" height="8" fill="currentColor" stroke="none"/>')
    lines.append('<text x="4" y="136" font-size="7" stroke="none" fill="currentColor" opacity="0.55">placeholder drawing</text>')
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 140" fill="none" stroke="currentColor" stroke-width="1.3">'
            + "".join(lines) + "</svg>")


def make_classes(rnd: random.Random):
    groups, classes, graph_files = [], [], {}
    gmap = {g: {"key": g, "title": t, "ghost_number": gn, "classes": [], "report": ref(g),
                "summary": f"<p>Sample summary of the group <em>{esc(t.lower())}</em>. The real text states which Green's functions "
                           f"the group contains, the number of graphs and what they determine.</p>"} for g, t, gn in GROUPS}
    for group, fields, n, pool in CLASSES:
        fields = ren_fields(fields)
        pool = [ren(x) for x in pool]
        cid = class_id(group, fields)
        g1 = group.startswith("g1")
        label = mrow(mi("Γ"), mo("["), *[SYM[f] for f in fields], *( [mo(";"), msub(mi("Δ"), mn(1))] if g1 else []), mo("]"))
        fname = "data/graphs/" + cid.replace(":", "__").replace(",", "_") + ".json"
        res_pole, _ = poly([(rnd.choice([-1, 1]) * rnd.randint(1, 43), rnd.choice([1, 2, 3, 6, 12, 18]),
                             [(rnd.choice(["gs", "g", "gp"]), 2)])])
        result = {"pole": M(mrow(KAPPA_POLE, res_pole, mo("⋅"), mi("T")), "sample pole")}
        if g1:
            result["breaking"] = M(mrow(KAPPA_FIN, res_pole, mo("⋅"), mi("N")), "sample breaking")
        else:
            result["breaking"] = None
        classes.append({"id": cid, "group": group, "ghost_number": 1 if g1 else 0, "fields": fields,
                        "label": M(label, f"Gamma[{', '.join(fields)}]"), "insertion": "Delta1" if g1 else None,
                        "n_graphs": n, "graphs": fname, "result": result,
                        "summary": f"<p>Sample summary of the class: {n} one-loop graphs with the external fields "
                                   f"{', '.join(fields)}{' and one insertion of the breaking' if g1 else ''}.</p>"})
        gmap[group]["classes"].append(cid)
        graphs = []
        topos = TOPO[min(max(len(fields), 2), 5)]
        for k in range(n):
            topo = topos[k % len(topos)] if rnd.random() < 0.8 else rnd.choice(topos)
            nl = {"tadpole": 1, "bubble": 2, "triangle": 3, "box": 4, "pentagon": 5}[topo]
            loop = [rnd.choice(pool) for _ in range(nl)]
            sf = rnd.choice(["1", "1", "1", "1/2", "1/2", "1/6"])
            ev = rnd.random() < 0.25
            graphs.append({
                "id": f"{cid.replace(':', '_').replace(',', '')}-{k + 1:03d}",
                "svg": svg_graph(fields, loop, g1, rnd),
                "loop": loop, "topology": topo, "insertion": "Delta1" if g1 else None,
                "symmetry_factor": sf, "sign": rnd.choice(["1", "-1"]),
                "caption": f"<p>{topo.capitalize()} with internal lines {', '.join(loop)}; symmetry factor {sf}.</p>",
                "contributes": {"pole": rnd.random() < 0.85, "evanescent": ev, "breaking": g1 and rnd.random() < 0.7},
            })
        graph_files[fname] = {"mock": True, "class": cid, "graphs": graphs}
    groups = [gmap[g] for g, _, _ in GROUPS]
    return groups, classes, graph_files


def find_class(classes, group, fields_subset=None, contains=None):
    out = []
    contains = [ren(x) for x in contains] if contains else None
    for c in classes:
        if c.get("group") != group:
            continue
        if contains and not all(x in c.get("fields", []) for x in contains):
            continue
        out.append(c["id"])
    return out


# --------------------------------------------------------------------------------------------------------------
# Terms of the unbroken phase
# --------------------------------------------------------------------------------------------------------------

def leg(fid: str, p: int) -> dict:
    return {"field": fid, "label": M(mrow(SYM.get(fid, mi(fid)), paren(msub(mi("p"), mn(p)))), f"{fid}(p_{p})")}


def term(tid, name, op, op_text, fields, sector, groups, types, theory, *, ghost_number=0, evanescent=False,
         tree=None, pole=None, finite=None, ev=None, brk=None, structure=None, classes=None, report=None,
         remarks=None, decomposition=None, symbol=None, rules=True, kind=None, parents=None, rnd=None):
    coeffs = {}
    if "tree" in types:
        coeffs["tree"] = {"math": coef(tree or [(1, 1, [])], None)}
    if "pole" in types:
        d = {"math": coef(pole or [(1, 1, [("g", 2)])], "pole")}
        if decomposition:
            d["decomposition"] = M(decomposition, "decomposition")
        d["note"] = None if rnd and rnd.random() < 0.6 else "<p>Sample note on the singular coefficient.</p>"
        coeffs["pole"] = d
    if "evanescent" in types:
        coeffs["evanescent"] = {"math": coef(ev or [(-2, 3, [("g", 2)])], "pole"),
                                "note": "<p>Sample note: the operator vanishes in four dimensions.</p>"}
    if "finite" in types:
        d = {"math": coef(finite or [(1, 6, [("g", 2)])], "finite")}
        if symbol:
            d["symbol"] = M(msub(mi("r"), mi(symbol)), f"r_{symbol}")
        coeffs["finite"] = d
    if "breaking" in types:
        coeffs["breaking"] = {"math": coef(brk or [(1, 1, [])], None)}
    op_m = M(op, op_text)
    total_parts = []
    for k in ("tree", "pole", "evanescent", "finite", "breaking"):
        if k in coeffs:
            inner = coeffs[k]["math"]["html"]
            body = inner[inner.index(">") + 1:-len("</math>")]
            if total_parts:
                total_parts.append(mo("+"))
            total_parts.append(body)
    total = M(mrow(paren(*total_parts) if len(total_parts) > 1 else "".join(total_parts), op), f"(...) {op_text}", block=True)
    t = {"id": tid, "name": name, "operator": op_m, "fields": fields, "legs": len(fields), "ghost_number": ghost_number,
         "sector": sector, "groups": groups, "evanescent": evanescent, "types": types, "coefficients": coeffs, "total": total,
         "classes": classes, "theory": theory, "remarks": remarks, "report": report or []}
    if rules:
        st = structure or mrow(mi("δ"), mo("…"))
        parts = {}
        for k in coeffs:
            inner = coeffs[k]["math"]["html"]
            body = inner[inner.index(">") + 1:-len("</math>")]
            parts[k] = M(mrow(body, mo("⋅"), st), f"{k} part", block=True)
        t["rules"] = [{"title": "Vertex function", "legs": [leg(f, j + 1) for j, f in enumerate(fields)], "parts": parts,
                       "note": "<p>Sample rule: all momenta incoming.</p>" if len(fields) > 2 else None}]
    else:
        t["rules"] = None
    if kind:
        t["kind"] = kind
    if parents is not None:
        t["parents"] = parents
    return t


def unbroken_terms(classes, rnd):
    T = []
    mu, nu = mi("μ"), mi("ν")
    gstruct = mrow(msup(mi("p"), mn(2)), msup(mi("g"), mrow(mu, nu)), mo("−"), msup(mi("p"), mu), msup(mi("p"), nu))
    for r in "GWB":
        grp = [GROUP_OF[r]]
        c = GROUP_CPL[r]
        T.append(term(f"gauge-kin-{r}", f"{GROUP_OF[r]} gauge kinetic term", mrow(fstrength(r), fstrength(r, True)),
                      f"F_{r} F_{r}", [r, r], "gauge", grp, ["tree", "pole", "finite"] if r != "B" else ["tree", "pole"], "gauge-kinetic",
                      tree=[(-1, 4, [])], pole=[(11, 6, [(c, 2)])], finite=[(-1, 3, [(c, 2)])], symbol=f"T{r}",
                      structure=gstruct, classes={"pole": find_class(classes, "g0_VV", contains=[r]), "breaking": find_class(classes, "g1_cV", contains=[r])},
                      report=ref("action", "g0_VV", "Sren1"), decomposition=mrow(mfrac(mn(1), mn(2)), msubsup(mi("δZ"), mi("A"), paren(mi(r)))),
                      remarks="<p>Sample remark: transversality of the singular part is checked in <a data-report-page=\"45\">section 3.1.1</a>.</p>" if r == "G" else None,
                      rnd=rnd))
        T.append(term(f"gauge-ev-{r}", f"Evanescent {GROUP_OF[r]} kinetic operator", mrow(hat(fstrength(r)), hat(fstrength(r, True))),
                      f"Fhat_{r} Fhat_{r}", [r, r], "gauge", grp, ["evanescent"], "ct-evanescent", evanescent=True,
                      ev=[(-1, 3, [(c, 2)])], classes={"pole": find_class(classes, "g0_VV", contains=[r])},
                      report=ref("evanescent", "Sren1"), rnd=rnd))
        T.append(term(f"gauge-mass-{r}", f"Restoring {GROUP_OF[r]} mass-like term", mrow(bar(gauge_sym(r)), bar(gauge_sym(r, up=True))),
                      f"Abar_{r} Abar_{r}", [r, r], "gauge", grp, ["finite"], "st-gauge",
                      finite=[(1, 6, [(c, 2)]), (1, 12, [("yu", 2)])] if r != "G" else [(1, 2, [(c, 2)])], symbol=f"M{r}",
                      classes={"pole": [], "breaking": find_class(classes, "g1_cV", contains=["c" + r])},
                      report=ref("st-solve", "st-result"), rnd=rnd))
        if r != "B":
            T.append(term(f"gauge-3-{r}", f"Cubic {GROUP_OF[r]} self-interaction", mrow(msup(mi("f"), mrow(mi("a"), mi("b"), mi("c"))), mo("("), mi("∂"), gauge_sym(r), mo(")"), gauge_sym(r), gauge_sym(r)),
                          f"f (d A) A A", [r, r, r], "gauge", grp, ["tree", "pole", "finite"], "gauge-self",
                          tree=[(-1, 1, [(c, 1)])], pole=[(-4, 3, [(c, 3)])], finite=[(1, 12, [(c, 3)])], symbol=f"V{r}",
                          classes={"pole": find_class(classes, "g0_VVV", contains=[r]), "breaking": find_class(classes, "g1_cVV", contains=["c" + r])},
                          report=ref("g0_VVV"), rnd=rnd))
            T.append(term(f"gauge-4-{r}", f"Quartic {GROUP_OF[r]} self-interaction", mrow(gauge_sym(r), gauge_sym(r), gauge_sym(r, "ν"), gauge_sym(r, "ν")),
                          "A A A A", [r, r, r, r], "gauge", grp, ["tree", "pole", "finite"], "gauge-self",
                          tree=[(-1, 4, [(c, 2)])], pole=[(2, 3, [(c, 4)])], finite=[(-1, 24, [(c, 4)])], symbol=f"Q{r}",
                          classes={"pole": find_class(classes, "g0_VVVV", contains=[r]), "breaking": find_class(classes, "g1_cVVV_cVVVV", contains=["c" + r])},
                          report=ref("g0_VVVV"), rnd=rnd))
            T.append(term(f"gauge-4-{r}-ev", f"Evanescent quartic {GROUP_OF[r]} operator", mrow(hat(gauge_sym(r)), gauge_sym(r), gauge_sym(r, "ν"), gauge_sym(r, "ν")),
                          "Ahat A A A", [r, r, r, r], "gauge", grp, ["evanescent"], "ct-evanescent", evanescent=True,
                          classes={"pole": find_class(classes, "g0_VVVV", contains=[r])}, rules=False, rnd=rnd))
    T.append(term("gauge-4-WWBB", "Restoring W W B B term", mrow(bar(gauge_sym("W")), bar(gauge_sym("W", up=True)), bar(gauge_sym("B", "ν")), bar(gauge_sym("B", "ν", up=True))),
                  "W W B B", ["W", "W", "B", "B"], "gauge", ["SU(2)", "U(1)"], ["finite"], "st-gauge",
                  finite=[(1, 8, [("g", 2), ("gp", 2)])], symbol="WB", classes={"pole": find_class(classes, "g0_VVVV", contains=["B"]), "breaking": []},
                  report=ref("st-result"), rnd=rnd))
    # gauge fixing and ghosts
    for r in "GWB":
        grp = [GROUP_OF[r]]
        T.append(term(f"gf-b-dA-{r}", f"{GROUP_OF[r]} gauge-fixing term", mrow(aux(r), mi("∂"), msup(mi(""), mi("μ")), gauge_sym(r)),
                      "b d.A", ["b" + r, r], "gauge-fixing", grp, ["tree"], "gauge-fixing", tree=[(1, 1, [])],
                      report=ref("action"), classes=None, rnd=rnd))
        T.append(term(f"gf-bb-{r}", f"{GROUP_OF[r]} auxiliary mass term", mrow(aux(r), aux(r)), "b b", ["b" + r, "b" + r],
                      "gauge-fixing", grp, ["tree"], "gauge-fixing", tree=[(1, 2, [])], report=ref("action"), rnd=rnd))
        T.append(term(f"ghost-kin-{r}", f"{GROUP_OF[r]} ghost kinetic term", mrow(ghost(r, True), mi("∂"), mo("·"), mi("∂"), ghost(r)),
                      "omegabar d.d omega", ["cb" + r, "c" + r], "ghost", grp, ["tree", "pole"] if r != "B" else ["tree"], "ghost",
                      ghost_number=0, tree=[(-1, 1, [])], pole=[(1, 2, [(GROUP_CPL[r], 2)])],
                      classes={"pole": find_class(classes, "g0_ghost_sources", contains=["cb" + r])}, report=ref("action", "g0_ghost_sources"), rnd=rnd))
        if r != "B":
            T.append(term(f"ghost-gauge-{r}", f"{GROUP_OF[r]} ghost-gauge vertex", mrow(ghost(r, True), mi("∂"), gauge_sym(r), ghost(r)),
                          "omegabar d A omega", ["cb" + r, r, "c" + r], "ghost", grp, ["tree", "pole"], "ghost",
                          tree=[(-1, 1, [(GROUP_CPL[r], 1)])], pole=[(1, 2, [(GROUP_CPL[r], 3)])],
                          classes={"pole": find_class(classes, "g0_ghost_sources")}, report=ref("g0_ghost_sources"), rnd=rnd))
    # fermions
    for f in FERMIONS:
        grps = [GROUP_OF[r] for r in F_GROUPS[f]]
        cl_ff = find_class(classes, "g0_FF", contains=[f])
        T.append(term(f"ferm-kin-{f}", f"Kinetic term of {f}", mrow(psibar(f), mi("i"), bar(mi("∂")), mo("/"), psi(f)),
                      f"Psibar_{f} i dslashbar Psi_{f}", [f, f], "fermion", grps, ["tree", "pole", "finite"] if f in "QL" else ["tree", "pole"],
                      "fermion-kinetic", tree=[(1, 1, [])], pole=[(-1, 2, [("g", 2)]), (-1, 8, [("gp", 2)])], finite=[(1, 4, [("gp", 2)])],
                      symbol=f"K{f}", decomposition=mrow(mfrac(mn(1), mn(2)), msub(mi("δZ"), mi(f))),
                      classes={"pole": cl_ff, "breaking": find_class(classes, "g1_cFF", contains=[f])}, report=ref("action", "g0_FF"), rnd=rnd))
        T.append(term(f"ferm-ev-kin-{f}", f"Evanescent kinetic term of {f} and its sterile partner",
                      mrow(psibar(f), mi("i"), hat(mi("∂")), mo("/"), chi(f)), f"Psibar_{f} i dslashhat chi_{f}", [f, "chi" + f],
                      "fermion", [], ["tree"], "fermion-kinetic", evanescent=True, tree=[(1, 1, [])], report=ref("fields", "action"),
                      remarks="<p>Sample remark: this is the only term of the classical action that is not BRST invariant in <i>d</i> dimensions; see <a data-theory=\"breaking\">the breaking</a>.</p>", rnd=rnd))
        T.append(term(f"ferm-sterile-kin-{f}", f"Kinetic term of the sterile partner of {f}", mrow(msub(bar(mi("χ")), mi(f)), mi("i"), bar(mi("∂")), mo("/"), chi(f)),
                      f"chibar i dslashbar chi", ["chi" + f, "chi" + f], "fermion", [], ["tree"], "fermion-kinetic", tree=[(1, 1, [])],
                      report=ref("fields"), rnd=rnd))
        for r in F_GROUPS[f]:
            T.append(term(f"ferm-gauge-{f}-{r}", f"{GROUP_OF[r]} vertex of {f}", mrow(psibar(f), bar(mi("γ")), msup(mi(""), mi("μ")), msub(mi("P"), mi("L" if f in "QL" else "R")), gauge_sym(r), psi(f)),
                          f"Psibar gammabar P A Psi", [f, f, r], "fermion", [GROUP_OF[r]], ["tree", "pole", "finite"], "fermion-gauge",
                          tree=[(1, 1, [(GROUP_CPL[r], 1)])], pole=[(-1, 1, [(GROUP_CPL[r], 3)]), (1, 4, [(GROUP_CPL[r], 1), ("yu", 2)])],
                          finite=[(rnd.choice([-1, 1]) * rnd.randint(1, 17), rnd.choice([6, 12, 18, 36]), [(GROUP_CPL[r], 3)])], symbol=f"F{f}{r}",
                          classes={"pole": find_class(classes, "g0_FFV", contains=[f, r]), "breaking": find_class(classes, "g1_cFFV", contains=[f])},
                          report=ref("g0_FFV", "st-result"), rnd=rnd))
            T.append(term(f"brk-{f}-{r}", f"Breaking vertex of {f} with the {GROUP_OF[r]} ghost",
                          mrow(ghost(r), psibar(f), hat(mi("γ")), msup(mi(""), mi("μ")), mi("∂"), msub(mi("μ"), mi("")), chi(f)),
                          "omega Psibar gammahat d chi", ["c" + r, f, "chi" + f], "breaking", [GROUP_OF[r]], ["breaking"], "breaking",
                          ghost_number=1, evanescent=True, brk=[(1, 1, [(GROUP_CPL[r], 1)])], report=ref("delta"), rnd=rnd))
    for r in "GWB":
        T.append(term(f"ferm-ev-{r}", f"Evanescent fermion operator of the {GROUP_OF[r]} sector", mrow(msubsup(mi("E"), mi("F"), paren(mi(r)))),
                      f"Etilde_F^{r}", ["Q", "Q", r], "fermion", [GROUP_OF[r]], ["evanescent"], "ct-evanescent", evanescent=True,
                      ev=[(-2, 3, [(GROUP_CPL[r], 2)])], classes={"pole": find_class(classes, "g0_FFV", contains=[r])}, report=ref("evanescent"), rnd=rnd))
    # Yukawa
    for f1, f2, y, hs in (("Q", "U", "yu", "Hd"), ("Q", "D", "yd", "H"), ("L", "E", "ye", "H")):
        for conj in (False, True):
            fl = [f2, f1, "H" if hs == "Hd" else "Hd"] if conj else [f1, f2, hs]
            T.append(term(f"yuk-{y}{'-hc' if conj else ''}", f"Yukawa coupling {CPL[y][1]}{' (conjugate)' if conj else ''}",
                          mrow(psibar(fl[0]), SYM[fl[2]], msub(mi("P"), mi("R" if not conj else "L")), psi(fl[1])), "Psibar H P Psi", fl,
                          "yukawa", [], ["tree", "pole"], "yukawa", tree=[(-1, 1, [(y, 1)])],
                          pole=[(-1, 8, [(y, 1), ("g", 2)]), (3, 2, [(y, 3)])], decomposition=msub(mi("δy"), mi(y[1])),
                          classes={"pole": find_class(classes, "g0_FFS", contains=[f1, f2]), "breaking": find_class(classes, "g1_cFFS", contains=[f1, f2])},
                          report=ref("g0_FFS"), rnd=rnd))
    # scalar sector
    T.append(term("higgs-kin", "Higgs kinetic term", mrow(paren(mi("∂"), HD), paren(mi("∂"), H)), "dH^+ dH", ["Hd", "H"], "scalar", [],
                  ["tree", "pole", "finite"], "higgs-kinetic", tree=[(1, 1, [])], pole=[(3, 4, [("g", 2)]), (1, 4, [("gp", 2)]), (-3, 2, [("yu", 2)])],
                  finite=[(1, 2, [("yu", 2)])], symbol="KH", decomposition=msub(mi("δZ"), mi("H")),
                  classes={"pole": find_class(classes, "g0_SS"), "breaking": find_class(classes, "g1_cSS_cSSV", contains=["Hd", "H"])}, report=ref("g0_SS"), rnd=rnd))
    for r in "WB":
        T.append(term(f"higgs-gauge-{r}", f"Higgs-{GROUP_OF[r]} vertex", mrow(HD, gauge_sym(r), mi("∂"), H), "H^+ A dH", ["Hd", "H", r], "scalar", [GROUP_OF[r]],
                      ["tree", "pole", "finite"], "higgs-kinetic", tree=[(-1, 1, [(GROUP_CPL[r], 1)])], pole=[(1, 1, [(GROUP_CPL[r], 3)])], finite=[(1, 8, [(GROUP_CPL[r], 1), ("yu", 2)])],
                      symbol=f"H{r}", classes={"pole": find_class(classes, "g0_SSV_SSVV", contains=[r]), "breaking": find_class(classes, "g1_cSS_cSSV")}, report=ref("g0_SSV_SSVV"), rnd=rnd))
    for r1, r2 in (("W", "W"), ("B", "B"), ("W", "B")):
        T.append(term(f"higgs-gauge2-{r1}{r2}", f"Higgs-gauge contact term {r1}{r2}", mrow(HD, gauge_sym(r1), gauge_sym(r2, up=True), H), "H^+ A A H",
                      ["Hd", "H", r1, r2], "scalar", sorted({GROUP_OF[r1], GROUP_OF[r2]}), ["tree", "pole"] + (["finite"] if r1 != r2 else []), "higgs-kinetic",
                      tree=[(1, 1, [(GROUP_CPL[r1], 1), (GROUP_CPL[r2], 1)])], finite=[(-1, 6, [("g", 1), ("gp", 1), ("yu", 2)])], symbol="HWB",
                      classes={"pole": find_class(classes, "g0_SSV_SSVV", contains=["W", "W"]) if r1 == "W" else [], "breaking": []}, report=ref("g0_SSV_SSVV"), rnd=rnd))
    T.append(term("higgs-mass", "Higgs mass term", mrow(HD, H), "H^+ H", ["Hd", "H"], "scalar", [], ["tree", "pole"], "higgs-potential",
                  tree=[(-1, 1, [("m2", 1)])], pole=[(9, 4, [("m2", 1), ("g", 2)]), (-6, 1, [("m2", 1), ("lam", 1)])], decomposition=msup(mi("δm"), mn(2)),
                  classes={"pole": find_class(classes, "g0_SS")}, report=ref("g0_SS", "Sren1"), rnd=rnd))
    T.append(term("higgs-quartic", "Higgs quartic coupling", msup(paren(HD, H), mn(2)), "(H^+ H)^2", ["Hd", "H", "Hd", "H"], "scalar", [], ["tree", "pole"], "higgs-potential",
                  tree=[(-1, 1, [("lam", 1)])], pole=[(-12, 1, [("lam", 2)]), (-9, 16, [("g", 4)]), (3, 1, [("yu", 4)])], decomposition=mi("δλ"),
                  classes={"pole": find_class(classes, "g0_SSSS")}, report=ref("g0_SSSS"), rnd=rnd))
    for nm, e in (("E_H", mrow(msub(mi("E"), mi("H")))), ("Etilde_H", mrow(msub(mi("Ẽ"), mi("H"))))):
        T.append(term(f"higgs-ev-{nm.replace('_', '')}", f"Evanescent Higgs operator {nm.replace('_', ' ')}", e, nm, ["Hd", "H"], "scalar", [], ["evanescent"], "ct-evanescent",
                      evanescent=True, ev=[(-1, 1, [("yu", 2)])], classes={"pole": find_class(classes, "g0_SS")}, report=ref("evanescent"), rnd=rnd))
    # BRST sources
    for r in "GWB":
        grp = [GROUP_OF[r]]
        T.append(term(f"src-rho-c-{r}", f"Source term of the {GROUP_OF[r]} gauge field, linear part", mrow(SYM["rho" + r], mi("∂"), msub(mi(""), mi("μ")), ghost(r)),
                      "rho d omega", ["rho" + r, "c" + r], "source", grp, ["tree", "pole"] if r != "B" else ["tree"], "brst-sources", ghost_number=0,
                      tree=[(1, 1, [])], pole=[(1, 2, [(GROUP_CPL[r], 2)])], classes={"pole": find_class(classes, "g0_ghost_sources", contains=["rho" + r])}, report=ref("brst"), rnd=rnd))
        if r != "B":
            T.append(term(f"src-rho-Ac-{r}", f"Source term of the {GROUP_OF[r]} gauge field, non-linear part", mrow(SYM["rho" + r], gauge_sym(r), ghost(r)),
                          "rho A omega", ["rho" + r, r, "c" + r], "source", grp, ["tree", "pole"], "brst-sources", tree=[(1, 1, [(GROUP_CPL[r], 1)])],
                          classes={"pole": find_class(classes, "g0_ghost_sources")}, report=ref("brst"), rnd=rnd))
            T.append(term(f"src-zeta-cc-{r}", f"Source term of the {GROUP_OF[r]} ghost", mrow(SYM["zeta" + r], ghost(r), ghost(r)), "zeta omega omega",
                          ["zeta" + r, "c" + r, "c" + r], "source", grp, ["tree", "pole"], "brst-sources", ghost_number=0, tree=[(1, 2, [(GROUP_CPL[r], 1)])],
                          classes={"pole": find_class(classes, "g0_ghost_sources", contains=["zeta" + r])}, report=ref("brst"), rnd=rnd))
    for f in FERMIONS:
        for r in F_GROUPS[f]:
            for barred in (False, True):
                src = ("Rb" if barred else "R") + f
                T.append(term(f"src-{src}-{r}", f"Source term of {'the conjugate of ' if barred else ''}{f}, {GROUP_OF[r]} part",
                              mrow(SYM[src], ghost(r), psi(f) if barred else psibar(f)), "R omega Psi", [src, "c" + r, f], "source", [GROUP_OF[r]],
                              ["tree", "finite"] if r == "G" else ["tree"], "st-source" if r == "G" else "brst-sources",
                              tree=[(1, 1, [(GROUP_CPL[r], 1)])], finite=[(-1, 3, [(GROUP_CPL[r], 3)])], symbol=f"R{f}",
                              classes={"pole": [], "breaking": find_class(classes, "g1_ccRF")} if r == "G" else None,
                              report=ref("brst", "g1_ccRF") if r == "G" else ref("brst"), rnd=rnd))
    T.append(term("src-sigma", "Source term of the Higgs doublet", mrow(SYM["sigma"], ghost("W"), H), "sigma omega H", ["sigma", "cW", "H"], "source", ["SU(2)"],
                  ["tree"], "brst-sources", tree=[(-1, 1, [("g", 1)])], report=ref("brst"), remarks=None, rnd=rnd))
    T.append(term("st-normalisation", "Kernel direction fixed by the normalisation conditions", mrow(mi("Q"), msub(mi(""), mi("N"))), "Q_N", ["W", "W"], "gauge", ["SU(2)"],
                  ["finite"], "normalisation", finite=[(0, 1, [])], report=ref("normalisation"), classes=None, rules=False,
                  remarks="<p>Sample remark: a term whose finite coefficient is fixed by a normalisation condition and not by the breaking.</p>", rnd=rnd))
    T.append(term("anomaly-check", "Would-be anomaly of the ghost-gauge-gauge functions", mrow(ghost("G"), msup(mi("ε"), mrow(mi("μ"), mi("ν"), mi("ρ"), mi("σ"))), mi("∂"), gauge_sym("G"), mi("∂"), gauge_sym("G")),
                  "omega eps dA dA", ["cG", "G", "G"], "gauge", ["SU(3)"], ["breaking"], "anomaly", ghost_number=1,
                  brk=[(0, 1, [])], classes={"pole": [], "breaking": find_class(classes, "g1_cVV")}, report=ref("g1_cVV"), rules=False, rnd=rnd))
    # padding: individual components of the fermion-gauge and Yukawa vertices to reach the expected size
    extra = 0
    for f in FERMIONS:
        for r in F_GROUPS[f]:
            for comp in ("vector", "axial"):
                T.append(term(f"ferm-gauge-{f}-{r}-{comp}", f"{comp.capitalize()} component of the {GROUP_OF[r]} vertex of {f}",
                              mrow(psibar(f), bar(mi("γ")), msup(mi(""), mi("μ")), *( [mi("γ"), msub(mi(""), mn(5))] if comp == "axial" else []), gauge_sym(r), psi(f)),
                              "Psibar gamma (gamma5) A Psi", [f, f, r], "fermion", [GROUP_OF[r]], ["pole", "finite"] if comp == "axial" else ["pole"], "st-fermion",
                              pole=[(1, 4, [(GROUP_CPL[r], 3)])], finite=[(1, 36, [(GROUP_CPL[r], 3)])], symbol=f"A{f}{r}",
                              classes={"pole": find_class(classes, "g0_FFV", contains=[f, r]), "breaking": find_class(classes, "g1_cFFV", contains=[f]) if comp == "axial" else []},
                              report=ref("g0_FFV"), rnd=rnd))
                extra += 1
    for f in FERMIONS:
        T.append(term(f"ferm-ev-ct-{f}", f"Evanescent kinetic counterterm of {f}", mrow(psibar(f), mi("i"), hat(mi("\u2202")), mo("/"), psi(f)),
                      "Psibar i dslashhat Psi", [f, f], "fermion", [], ["evanescent"], "ct-evanescent", evanescent=True,
                      ev=[(-1, 6, [("g", 2)])], classes={"pole": find_class(classes, "g0_FF", contains=[f])}, report=ref("evanescent"), rnd=rnd))
        for r in F_GROUPS[f]:
            T.append(term(f"ferm-gauge-ev-{f}-{r}", f"Evanescent {GROUP_OF[r]} vertex of {f}", mrow(psibar(f), hat(mi("\u03b3")), msup(mi(""), mi("\u03bc")), gauge_sym(r), psi(f)),
                          "Psibar gammahat A Psi", [f, f, r], "fermion", [GROUP_OF[r]], ["evanescent"], "ct-evanescent", evanescent=True,
                          ev=[(-2, 3, [(GROUP_CPL[r], 3)])], classes={"pole": find_class(classes, "g0_FFV", contains=[f, r])}, report=ref("evanescent"),
                          rules=rnd.random() < 0.5, rnd=rnd))
    for t in T:
        t["fields"] = ren_fields(t["fields"])
        t["legs"] = len(t["fields"])
        for r in t.get("rules") or []:
            r["legs"] = [leg(f, j + 1) for j, f in enumerate(t["fields"])]
    return T


# --------------------------------------------------------------------------------------------------------------
# Broken phase
# --------------------------------------------------------------------------------------------------------------

BROKEN_FIELDS = [
    ("A", "A", "photon", "gauge", "U(1)", 0, "bose"), ("Z", "Z", "Z boson", "gauge", None, 0, "bose"),
    ("Wp", "W⁺", "W boson (positive charge)", "gauge", "SU(2)", 0, "bose"), ("Wm", "W⁻", "W boson (negative charge)", "gauge", "SU(2)", 0, "bose"),
    ("G", "G", "gluon", "gauge", "SU(3)", 0, "bose"),
    ("u", "u", "up quark", "fermion", None, 0, "fermi"), ("d", "d", "down quark", "fermion", None, 0, "fermi"),
    ("e", "e", "electron", "fermion", None, 0, "fermi"), ("nu", "ν", "neutrino", "fermion", None, 0, "fermi"),
    ("chiu", "χᵤ", "sterile partner of u", "sterile", None, 0, "fermi"),
    ("h", "h", "Higgs boson", "scalar", None, 0, "bose"), ("G0", "G⁰", "neutral Goldstone boson", "scalar", None, 0, "bose"),
    ("Gp", "G⁺", "charged Goldstone boson", "scalar", None, 0, "bose"), ("Gm", "G⁻", "charged Goldstone boson", "scalar", None, 0, "bose"),
    ("cA", "ωᴀ", "photon ghost", "ghost", None, 1, "fermi"), ("cZ", "ωᶻ", "Z ghost", "ghost", None, 1, "fermi"),
    ("cWp", "ω⁺", "W ghost", "ghost", None, 1, "fermi"), ("cWm", "ω⁻", "W ghost", "ghost", None, 1, "fermi"),
    ("cG", "ω", "gluon ghost", "ghost", "SU(3)", 1, "fermi"),
    ("cbA", "ω̄", "photon antighost", "antighost", None, -1, "fermi"),
    ("bA", "b", "auxiliary field of the photon", "auxiliary", None, 0, "bose"),
]


def broken_data(unbroken, classes, rnd):
    fields = []
    for fid, sym, name, kind, grp, gh, stat in BROKEN_FIELDS:
        fields.append({"id": fid, "symbol": M(mi(sym), sym), "name": name, "kind": kind, "group": grp, "ghost_number": gh, "statistics": stat, "note": ""})
        SYM.setdefault("b:" + fid, mi(sym))
    ub = {t["id"]: t for t in unbroken}
    T = []

    def bt(tid, name, fl, kind, types, theory, parents, sector, groups=(), **kw):
        op = mrow(*[mi(dict((x[0], x[1]) for x in BROKEN_FIELDS)[f]) for f in fl])
        if not groups:
            gs = set()
            for f in fl:
                if f in ("G", "cG"):
                    gs.add("SU(3)")
                if f in ("Wp", "Wm", "Z", "cWp", "cWm", "cZ"):
                    gs.add("SU(2)")
                if f in ("A", "Z", "cA", "bA"):
                    gs.update(("SU(2)", "U(1)"))
            groups = sorted(gs)
        t = term(tid, name, op, " ".join(fl), list(fl), sector, list(groups), types, theory, kind=kind, parents=parents, rnd=rnd,
                 classes=None, **kw)
        # legs refer to broken fields
        if t["rules"]:
            for r in t["rules"]:
                r["legs"] = [{"field": f, "label": M(mrow(mi(dict((x[0], x[1]) for x in BROKEN_FIELDS)[f]), paren(msub(mi("p"), mn(j + 1)))), f)} for j, f in enumerate(fl)]
        T.append(t)

    v2 = [("v", 2)]
    for b, nm, par in (("A", "photon", ["gauge-kin-W", "gauge-kin-B"]), ("Z", "Z", ["gauge-kin-W", "gauge-kin-B", "higgs-kin"]),
                       ("Wp", "W", ["gauge-kin-W", "higgs-kin"]), ("G", "gluon", ["gauge-kin-G"])):
        fl = [b, b] if b != "Wp" else ["Wp", "Wm"]
        bt(f"b-kin-{b}", f"Kinetic term of the {nm}", fl, "kinetic", ["tree", "pole", "finite"] if b != "G" else ["tree", "pole", "finite"], "broken-shift", par, "gauge",
           tree=[(-1, 4, [])], pole=[(11, 6, [("g", 2)])], finite=[(-1, 3, [("g", 2)])], report=ref("bl-tree", "bl-ST"))
        if b in ("Z", "Wp"):
            bt(f"b-mass-{b}", f"Mass term of the {nm}", fl, "mass", ["tree", "pole", "finite"], "broken-masses", par + ["higgs-kin"], "gauge",
               tree=[(1, 8, [("g", 2), ("v", 2)])], pole=[(3, 8, [("g", 4), ("v", 2)])], finite=[(1, 16, [("g", 2), ("yu", 2), ("v", 2)])],
               report=ref("bl-tree", "bl-ST"), remarks="<p>Sample remark: the finite restoring term proportional to <i>v</i>² enters the mass counterterm.</p>")
    bt("b-mass-A", "Photon mass-like restoring term", ["A", "A"], "mass", ["finite"], "broken-st", ["gauge-mass-W", "gauge-mass-B"], "gauge",
       finite=[(0, 1, [])], report=ref("bl-ST"), remarks="<p>Sample remark: the photon remains massless; see <a data-report-page=\"166\">section 4.4.6</a>.</p>")
    bt("b-mix-AZ", "Photon-Z mixing", ["A", "Z"], "mixing", ["pole", "finite"], "broken-mixing", ["gauge-kin-W", "gauge-kin-B", "gauge-mass-W"], "gauge",
       pole=[(1, 6, [("g", 1), ("gp", 1)])], finite=[(1, 12, [("g", 1), ("gp", 1)])], report=ref("bl-ST"))
    bt("b-mix-ZG0", "Z-Goldstone mixing", ["Z", "G0"], "mixing", ["tree", "pole"], "broken-mixing", ["higgs-kin", "higgs-gauge-W", "higgs-gauge-B"], "scalar",
       tree=[(1, 2, [("g", 1), ("v", 1)])], report=ref("bl-tree"))
    bt("b-mix-WG", "W-Goldstone mixing", ["Wp", "Gm"], "mixing", ["tree", "pole"], "broken-mixing", ["higgs-kin", "higgs-gauge-W"], "scalar",
       tree=[(1, 2, [("g", 1), ("v", 1)])], report=ref("bl-tree"))
    bt("b-tad-h", "Higgs tadpole", ["h"], "tadpole", ["tree", "pole"], "broken-tadpole", ["higgs-mass", "higgs-quartic"], "scalar",
       tree=[(-1, 1, [("m2", 1), ("v", 1)]), (-1, 1, [("lam", 1), ("v", 1)])], pole=[(1, 1, [("lam", 1), ("m2", 1), ("v", 1)])], report=ref("bp-tadpole"), rules=True)
    bt("b-kin-h", "Higgs kinetic term", ["h", "h"], "kinetic", ["tree", "pole", "finite"], "broken-shift", ["higgs-kin"], "scalar",
       tree=[(1, 2, [])], pole=[(-3, 2, [("yu", 2)])], finite=[(1, 2, [("yu", 2)])], report=ref("bl-tree"))
    bt("b-mass-h", "Higgs mass term", ["h", "h"], "mass", ["tree", "pole"], "broken-masses", ["higgs-mass", "higgs-quartic"], "scalar",
       tree=[(-1, 1, [("lam", 1), ("v", 2)])], pole=[(3, 1, [("lam", 2), ("v", 2)])], report=ref("bl-tree"))
    for gname, gl in (("G0", ["G0", "G0"]), ("Gpm", ["Gp", "Gm"])):
        bt(f"b-kin-{gname}", f"Kinetic term of the Goldstone boson {gname}", gl, "kinetic", ["tree", "pole"], "broken-shift", ["higgs-kin"], "scalar", tree=[(1, 2, [])], report=ref("bl-tree"))
        bt(f"b-mass-{gname}", f"Goldstone mass term {gname}", gl, "mass", ["tree", "pole"], "broken-masses", ["higgs-mass", "higgs-quartic"], "scalar",
           tree=[(-1, 1, [("m2", 1)]), (-1, 1, [("lam", 1), ("v", 2)])], report=ref("bl-tree"))
    for f, y, par in (("u", "yu", ["yuk-yu", "yuk-yu-hc", "ferm-kin-Q", "ferm-kin-U"]), ("d", "yd", ["yuk-yd", "yuk-yd-hc", "ferm-kin-Q", "ferm-kin-D"]),
                      ("e", "ye", ["yuk-ye", "yuk-ye-hc", "ferm-kin-L", "ferm-kin-E"])):
        bt(f"b-kin-{f}", f"Kinetic term of the {f} quark" if f != "e" else "Kinetic term of the electron", [f, f], "kinetic", ["tree", "pole", "finite"], "broken-shift", par[2:], "fermion",
           tree=[(1, 1, [])], pole=[(-1, 2, [("g", 2)])], finite=[(1, 4, [("gp", 2)])], report=ref("bl-tree"))
        bt(f"b-mass-{f}", f"Mass term of the {f} quark" if f != "e" else "Mass term of the electron", [f, f], "mass", ["tree", "pole"], "broken-masses", par, "fermion",
           tree=[(-1, 2, [(y, 1), ("v", 1)])], pole=[(3, 4, [(y, 3), ("v", 1)])], report=ref("bl-tree"),
           remarks="<p>Sample remark: the mass term couples the physical chiralities only; the sterile partners stay massless.</p>")
        bt(f"b-yuk-{f}", f"Higgs coupling of the {f} quark" if f != "e" else "Higgs coupling of the electron", [f, f, "h"], "interaction", ["tree", "pole"], "broken-interactions", par[:2], "yukawa",
           tree=[(-1, 2, [(y, 1)])], pole=[(3, 4, [(y, 3)])], report=ref("bl-tree"))
        bt(f"b-gold-{f}", f"Neutral Goldstone coupling of the {f} fermion", [f, f, "G0"], "interaction", ["tree", "pole"], "broken-interactions", par[:2], "yukawa",
           tree=[(-1, 2, [(y, 1)])], report=ref("bl-tree"))
        for b, par2 in (("A", ["ferm-gauge-Q-B", "ferm-gauge-Q-W"]), ("Z", ["ferm-gauge-Q-B", "ferm-gauge-Q-W"])):
            bt(f"b-ffv-{f}-{b}", f"{'Photon' if b == 'A' else 'Z'} vertex of the {f} fermion", [f, f, b], "interaction", ["tree", "pole", "finite"], "broken-interactions", par2, "fermion",
               tree=[(-1, 3, [("g", 1)])], pole=[(1, 4, [("g", 3)])], finite=[(1, 12, [("g", 3)])], report=ref("bl-tree", "bl-ST"))
        bt(f"b-brk-{f}", f"Breaking vertex of the {f} fermion with the photon ghost", ["cA", f, "chiu" if f == "u" else f], "interaction", ["breaking"], "breaking",
           [f"brk-{'Q' if f != 'e' else 'L'}-B"], "breaking", ghost_number=1, brk=[(1, 1, [("g", 1)])], report=ref("delta"))
    bt("b-ffv-e-W", "W vertex of the leptons", ["e", "nu", "Wp"], "interaction", ["tree", "pole", "finite"], "broken-interactions", ["ferm-gauge-L-W"], "fermion",
       tree=[(-1, 2, [("g", 1)])], report=ref("bl-tree"))
    bt("b-ffv-ud-W", "W vertex of the quarks", ["u", "d", "Wp"], "interaction", ["tree", "pole", "finite"], "broken-interactions", ["ferm-gauge-Q-W"], "fermion",
       tree=[(-1, 2, [("g", 1)])], report=ref("bl-tree"))
    for k, (fl, nm, par) in enumerate([(["h", "h", "h"], "Cubic Higgs self-coupling", ["higgs-quartic"]), (["h", "h", "h", "h"], "Quartic Higgs self-coupling", ["higgs-quartic"]),
                                        (["h", "G0", "G0"], "Higgs-Goldstone coupling", ["higgs-quartic"]), (["h", "Gp", "Gm"], "Higgs-charged Goldstone coupling", ["higgs-quartic"]),
                                        (["h", "Z", "Z"], "Higgs-Z-Z coupling", ["higgs-gauge2-WW", "higgs-gauge2-BB", "higgs-gauge2-WB"]),
                                        (["h", "Wp", "Wm"], "Higgs-W-W coupling", ["higgs-gauge2-WW"]),
                                        (["h", "h", "Z", "Z"], "Higgs-Higgs-Z-Z coupling", ["higgs-gauge2-WW", "higgs-gauge2-BB"]),
                                        (["h", "h", "Wp", "Wm"], "Higgs-Higgs-W-W coupling", ["higgs-gauge2-WW"]),
                                        (["Wp", "Wm", "A"], "W-W-photon vertex", ["gauge-3-W"]), (["Wp", "Wm", "Z"], "W-W-Z vertex", ["gauge-3-W"]),
                                        (["Wp", "Wm", "A", "A"], "W-W-photon-photon contact term", ["gauge-4-W"]), (["Wp", "Wm", "Z", "Z"], "W-W-Z-Z contact term", ["gauge-4-W"]),
                                        (["Wp", "Wm", "Wp", "Wm"], "Quartic W self-interaction", ["gauge-4-W"]), (["G", "G", "G"], "Cubic gluon self-interaction", ["gauge-3-G"]),
                                        (["G", "G", "G", "G"], "Quartic gluon self-interaction", ["gauge-4-G"]), (["h", "A", "Z"], "Higgs-photon-Z restoring term", ["higgs-gauge2-WB"]),
                                        (["G0", "h", "Z"], "Goldstone-Higgs-Z vertex", ["higgs-gauge-W", "higgs-gauge-B"]), (["Gp", "h", "Wm"], "Charged Goldstone-Higgs-W vertex", ["higgs-gauge-W"]),
                                        (["Gp", "Gm", "A"], "Charged Goldstone-photon vertex", ["higgs-gauge-W", "higgs-gauge-B"]), (["Gp", "Gm", "Z"], "Charged Goldstone-Z vertex", ["higgs-gauge-W", "higgs-gauge-B"])]):
        types = ["tree", "pole"] + (["finite"] if k % 3 == 0 else [])
        if nm.endswith("restoring term"):
            types = ["finite"]
        bt(f"b-int-{k:02d}", nm, fl, "interaction", types, "broken-interactions" if "finite" not in types or nm.endswith("restoring term") is False else "broken-st", par,
           "scalar" if fl[0] in ("h", "G0", "Gp") else "gauge", tree=[(-1, 1, [("lam", 1), ("v", 1)])], finite=[(1, 12, [("g", 2)])], report=ref("bl-tree"))
    for k, (fl, nm, par) in enumerate([(["u", "u", "G"], "Gluon vertex of the u quark", ["ferm-gauge-Q-G", "ferm-gauge-U-G"]),
                                        (["d", "d", "G"], "Gluon vertex of the d quark", ["ferm-gauge-Q-G", "ferm-gauge-D-G"]),
                                        (["nu", "nu", "Z"], "Z vertex of the neutrino", ["ferm-gauge-L-W", "ferm-gauge-L-B"]),
                                        (["u", "d", "Gp"], "Charged Goldstone coupling of the quarks", ["yuk-yu", "yuk-yd"]),
                                        (["e", "nu", "Gm"], "Charged Goldstone coupling of the leptons", ["yuk-ye"]),
                                        (["cbA", "cWp", "Wm"], "Ghost-W vertex", ["ghost-gauge-W"]),
                                        (["cbA", "cG", "G"], "Ghost-gluon vertex", ["ghost-gauge-G"]),
                                        (["h", "h", "G0", "G0"], "Higgs-Higgs-Goldstone contact term", ["higgs-quartic"]),
                                        (["G0", "G0", "Gp", "Gm"], "Goldstone contact term", ["higgs-quartic"]),
                                        (["Gp", "Gm", "A", "A"], "Charged Goldstone-photon contact term", ["higgs-gauge2-WW", "higgs-gauge2-BB"]),
                                        (["Gp", "Gm", "Z", "Z"], "Charged Goldstone-Z contact term", ["higgs-gauge2-WW", "higgs-gauge2-BB"]),
                                        (["h", "Z", "Z"], "Higgs-Z-Z restoring term", ["higgs-gauge2-WB"]),
                                        (["Z", "Z"], "Z mass-like restoring term", ["gauge-mass-W", "gauge-mass-B"]),
                                        (["Wp", "Wm"], "W mass-like restoring term", ["gauge-mass-W"])]):
        types = ["finite"] if "restoring" in nm else ["tree", "pole"]
        bt(f"b-ext-{k:02d}", nm, fl, "mass" if len(fl) == 2 else "interaction", types, "broken-st" if "restoring" in nm else "broken-interactions", par,
           "ghost" if fl[0] == "cbA" else ("gauge" if fl[0] in ("Z", "Wp") else "scalar"), tree=[(-1, 2, [("g", 1)])], finite=[(1, 24, [("g", 2), ("v", 2)])],
           report=ref("bl-ST" if "restoring" in nm else "bl-tree"), ghost_number=0)
    for b, par in (("A", ["ghost-kin-W", "ghost-kin-B"]), ("Z", ["ghost-kin-W", "ghost-kin-B"]), ("Wp", ["ghost-kin-W"]), ("G", ["ghost-kin-G"])):
        gh = {"A": "cA", "Z": "cZ", "Wp": "cWp", "G": "cG"}[b]
        bt(f"b-ghost-kin-{b}", f"Ghost kinetic term ({b})", ["cbA" if b == "A" else "cbA", gh], "kinetic", ["tree", "pole"] if b != "A" else ["tree"], "ghost", par, "ghost", tree=[(-1, 1, [])], report=ref("bl-tree"))
    bt("b-gf-A", "Photon gauge-fixing term", ["bA", "A"], "kinetic", ["tree"], "gauge-fixing", ["gf-b-dA-W", "gf-b-dA-B"], "gauge-fixing", tree=[(1, 1, [])], report=ref("bl-tree"),
       remarks="<p>Sample remark: the gauge fixing is linear; no R<sub>ξ</sub> term is introduced by the shift.</p>")
    # unknown key and missing optional keys, to exercise robustness
    T[-1]["future_key"] = {"ignored": True}
    del T[-2]["remarks"]
    # check parents exist
    for t in T:
        t["parents"] = list(dict.fromkeys(p for p in t["parents"] if p in ub))
    facets_kind = [{"id": k, "label": k.capitalize()} for k in ("kinetic", "mass", "mixing", "tadpole", "interaction")]
    return fields, T, facets_kind


# --------------------------------------------------------------------------------------------------------------
# Theory cards
# --------------------------------------------------------------------------------------------------------------

CARDS = [
    ("gauge-kinetic", "Gauge kinetic terms", "gauge kinetic terms in d dimensions, field renormalisation, transversality", ("action", "g0_VV")),
    ("gauge-self", "Gauge self-interactions", "cubic and quartic gauge self-interactions", ("g0_VVV", "g0_VVVV")),
    ("gauge-fixing", "Gauge fixing", "auxiliary field, linear gauge fixing, Feynman gauge, why it is not renormalised", ("action",)),
    ("ghost", "Ghosts", "ghost kinetic term and ghost-gauge vertex, ghost equation", ("action", "g0_ghost_sources")),
    ("fermion-kinetic", "Fermion kinetic terms and sterile partners", "d-dimensional kinetic term, sterile partners, origin of the breaking", ("fields", "action")),
    ("fermion-gauge", "Fermion-gauge vertices", "four-dimensional (barred) gauge vertex of the physical chirality", ("g0_FFV",)),
    ("yukawa", "Yukawa couplings", "Yukawa couplings", ("g0_FFS",)),
    ("higgs-kinetic", "Higgs kinetic term", "scalar kinetic term and scalar-gauge vertices", ("g0_SS", "g0_SSV_SSVV")),
    ("higgs-potential", "Higgs potential", "mass and quartic terms", ("g0_SS", "g0_SSSS")),
    ("brst-sources", "BRST sources", "source terms of the non-linear BRST variations, their role in the Slavnov-Taylor operator", ("brst",)),
    ("breaking", "The breaking insertion", "the insertion of the tree-level breaking and its Feynman rule", ("delta",)),
    ("ct-fields", "Field renormalisation", "MS-bar field renormalisation constants", ("Sren1",)),
    ("ct-couplings", "Coupling and mass counterterms", "MS-bar coupling and mass counterterms, beta functions", ("Sren1",)),
    ("ct-evanescent", "Evanescent counterterms", "singular counterterms on evanescent operators", ("evanescent", "Sren1")),
    ("st-gauge", "Restoring counterterms: gauge fields", "finite restoring counterterms with gauge fields only (two-, three-, four-point)", ("st-solve", "st-result")),
    ("st-fermion", "Restoring counterterms: fermion-gauge vertices", "finite restoring counterterms of fermion-gauge vertices", ("st-result",)),
    ("st-source", "Restoring counterterms: source terms", "finite counterterms of the source terms (renormalised BRST transformations)", ("st-result",)),
    ("st-scalar", "Restoring counterterms: scalars", "finite restoring counterterms with scalars", ("st-result",)),
    ("normalisation", "Normalisation conditions", "freedom of the finite counterterm (kernel of b) and the normalisation conditions", ("normalisation",)),
    ("anomaly", "The would-be anomaly", "the would-be anomaly and its cancellation", ("g1_cVV",)),
    ("broken-shift", "The shift to the broken phase", "the shift of the Higgs doublet in the linear gauge, no R_xi input", ("Sren-broken",)),
    ("broken-masses", "Masses and mass counterterms", "mass terms and mass counterterms, fermion masses with sterile partners", ("bl-tree",)),
    ("broken-mixing", "Mixing", "gauge-Goldstone and neutral gauge mixing in the linear gauge", ("bl-tree",)),
    ("broken-tadpole", "Tadpoles", "tadpole terms and the definition of the vacuum expectation value", ("bp-tadpole",)),
    ("broken-st", "Restoring counterterms after the shift", "finite restoring counterterms after the shift (e.g. gauge two-point terms proportional to v squared)", ("bl-ST",)),
    ("broken-interactions", "Interactions generated by the shift", "interaction terms generated by the shift", ("bl-tree",)),
]


def theory_cards(class_example: str):
    cards = {}
    for key, title, scope, refs in CARDS:
        body = (f"<p><strong>Sample text.</strong> This card will explain the {esc(scope)}. The real card is written from the report; "
                f"this placeholder only fixes the length and the markup of a card.</p>"
                f"<p>It may contain inline mathematics such as "
                f"{math(mrow(mi('d'), mo('='), mn(4), mo(chr(0x2212)), mn(2), mi(LUNATE)))}, references to the report "
                f"(<a data-report-page=\"{REF[refs[0]][1]}\">{REF[refs[0]][0]}</a>), to another card "
                f"(<a data-theory=\"breaking\">the breaking</a>), to a term (<a data-term=\"unbroken:gauge-kin-G\">the gluon kinetic term</a>) "
                f"and to a class (<a data-class=\"{esc(class_example)}\">a class of Green's functions</a>).</p>")
        cards[key] = {"title": title, "body": body, "report": ref(*refs)}
    return cards


# --------------------------------------------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------------------------------------------

def is_mock_or_absent(path: Path) -> bool:
    if not path.exists():
        return True
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if isinstance(data, dict):
        return data.get("mock") is True
    if isinstance(data, list):
        return bool(data) and all(isinstance(x, dict) and x.get("mock") is True for x in data)
    return False


def write(site: Path, rel: str, data, written: list, skipped: list) -> None:
    path = site / rel
    if not is_mock_or_absent(path):
        skipped.append(rel)
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    written.append(rel)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--site", default=str(DEFAULT_SITE))
    args = ap.parse_args()
    site = Path(args.site).resolve()
    rnd = random.Random(20261128)

    ufields = unbroken_fields()
    for f in ufields:
        KIND_OF[f["id"]] = f["kind"]
    fill_sym(ufields)
    groups, classes, graph_files = make_classes(rnd)
    real_classes = None
    cpath = site / "data/classes.json"
    if cpath.exists() and not is_mock_or_absent(cpath):
        real_classes = json.loads(cpath.read_text(encoding="utf-8"))
        print("classes.json is real: the mock terms refer to its classes; no mock classes or graphs are written")
    uterms = unbroken_terms(real_classes["classes"] if real_classes else classes, rnd)
    bfields, bterms, kinds = broken_data(uterms, classes, rnd)

    sectors = [{"id": s, "label": l} for s, l in (("gauge", "Gauge"), ("gauge-fixing", "Gauge fixing"), ("ghost", "Ghost"), ("fermion", "Fermion"),
                                                   ("yukawa", "Yukawa"), ("scalar", "Scalar"), ("source", "BRST source"), ("breaking", "Breaking"))]
    types = [{"id": "tree", "label": "Tree"}, {"id": "pole", "label": "Singular"}, {"id": "evanescent", "label": "Evanescent"},
             {"id": "finite", "label": "Finite restoring"}, {"id": "breaking", "label": "Breaking"}]
    gfacet = [{"id": g, "label": g} for g in ("SU(3)", "SU(2)", "U(1)")]

    kappa = mfrac(mn(1), mrow(mn(16), msup(mi(PI), mn(2))))
    s_ren = M(mrow(msub(mi("S"), mi("ren")), mo("="), msub(mi("S"), mn(0)), mo("+"), msup(mrow(mi("δ"), mi("S")), paren(mn(1))), mo("+"),
                   msubsup(mrow(mi("Δ"), mi("S")), mi("ST"), paren(mn(1)))), "S_ren = S_0 + delta S^(1) + Delta S^(1)_ST", block=True)
    d_s = M(mrow(msubsup(mrow(mi("Δ"), mi("S")), mi("ST"), paren(mn(1))), mo("="), kappa, mrow(mo("∑"), msub(mi(""), mi("α"))),
                 msub(mi("r"), mi("α")), msub(mi("Q"), mi("α"))), "Delta S = 1/(16 pi^2) sum r_alpha Q_alpha", block=True)
    delta_s = M(mrow(msup(mrow(mi("δ"), mi("S")), paren(mn(1))), mo("="), mo("∑"), mfrac(mrow(mi("δ"), msub(mi("Z"), mi("A"))), mn(2)), msub(mi("O"), mi("A")),
                     mo("+"), mo("…")), "delta S^(1) = ...", block=True)
    unbroken = {
        "mock": True, "phase": "unbroken", "title": "Unbroken phase",
        "intro": "<p><strong>Sample introduction.</strong> The renormalised action of the unbroken phase, operator by operator. Each row is one operator "
                 "with definite field content; its coefficients at tree level, in the singular counterterm and in the finite restoring counterterm "
                 "are listed with the Green's functions that determine them.</p>",
        "overview": [
            {"id": "S-ren", "title": "Renormalised action", "math": s_ren, "text": "<p>Sample text: the renormalised action at one loop, in three parts.</p>"},
            {"id": "dS", "title": "Singular counterterm", "math": delta_s, "text": "<p>Sample text on the singular counterterm.</p>"},
            {"id": "DS", "title": "Finite restoring counterterm", "math": d_s, "text": None},
        ],
        "fields": ufields, "facets": {"sector": sectors, "type": types, "group": gfacet}, "terms": uterms,
    }
    broken = {
        "mock": True, "phase": "broken", "title": "Broken phase",
        "intro": "<p><strong>Sample introduction.</strong> The renormalised action after the shift of the Higgs doublet, in the mass eigenbasis.</p>",
        "overview": [{"id": "S-ren-broken", "title": "Renormalised action after the shift", "math": s_ren, "text": "<p>Sample text.</p>"}],
        "fields": bfields, "facets": {"sector": [x for x in sectors if x["id"] != "source"], "type": [x for x in types if x["id"] != "evanescent"], "group": gfacet, "kind": kinds}, "terms": bterms,
    }
    manifest = {
        "mock": True,
        "title": "One-loop renormalisation of the Standard Model with a non-anticommuting γ₅",
        "subtitle": "BMHV scheme with sterile partners, MS-bar, finite Slavnov-Taylor-restoring counterterms, Feynman gauge, one fermion family.",
        "model": "sm_one_family", "loop_order": 1, "generated": "2026-11-28",
        "report_pdf": "report/one_loop_gamma5.pdf", "report_pages": 186,
        "phases": [{"id": "unbroken", "title": "Unbroken phase", "file": "data/unbroken.json"},
                   {"id": "broken", "title": "Broken phase", "file": "data/broken.json"}],
        "classes": "data/classes.json", "theory": "data/theory.json", "resources": "data/resources.json",
        "status": "<p><strong>Sample status.</strong> This data set is a placeholder generated by <code>tools/mock_data.py</code>; it has the shape of "
                  "the real data and none of its content.</p>",
    }
    resources = []
    res_dir = site / "resources"
    descr = {"final_counterterms.json": "Final counterterm coefficients of the report.", "report_headings.txt": "First page of every labelled heading of the report."}
    if res_dir.is_dir():
        for p in sorted(res_dir.rglob("*")):
            if p.is_file():
                rel = p.relative_to(site).as_posix()
                resources.append({"mock": True, "path": rel, "title": p.name, "description": descr.get(p.name, "Sample description."), "bytes": p.stat().st_size})
    resources.append({"mock": True, "path": "report/one_loop_gamma5.pdf", "title": "Report (PDF)", "description": "The report.",
                      "bytes": (site / "report/one_loop_gamma5.pdf").stat().st_size if (site / "report/one_loop_gamma5.pdf").exists() else None})

    written, skipped = [], []
    write(site, "data/manifest.json", manifest, written, skipped)
    write(site, "data/unbroken.json", unbroken, written, skipped)
    write(site, "data/broken.json", broken, written, skipped)
    write(site, "data/classes.json", {"mock": True, "groups": groups, "classes": classes}, written, skipped)
    cex = (real_classes["classes"] if real_classes else classes)[0]["id"]
    write(site, "data/theory.json", {"mock": True, "cards": theory_cards(cex)}, written, skipped)
    write(site, "data/resources.json", resources, written, skipped)
    removed = 0
    if real_classes is None:
        for rel, data in graph_files.items():
            write(site, rel, data, written, skipped)
    else:
        used = {c.get("graphs") for c in real_classes.get("classes", [])}
        for path in sorted((site / "data/graphs").glob("*.json")):
            rel = path.relative_to(site).as_posix()
            if rel not in used and is_mock_or_absent(path):
                path.unlink()
                removed += 1
        if removed:
            print(f"removed {removed} mock graph files that the real classes.json does not use")
    n_graphs = sum(len(g["graphs"]) for g in graph_files.values())
    print(f"mock data: {len(uterms)} unbroken terms, {len(bterms)} broken terms, {len(classes)} classes in {len(groups)} groups, "
          f"{n_graphs} graphs, {len(CARDS)} theory cards; wrote {len(written)} files, kept {len(skipped)} real files")
    for s in skipped:
        print(f"  kept (not mock): {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
