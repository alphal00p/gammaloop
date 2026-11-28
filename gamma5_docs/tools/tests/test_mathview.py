"""Fast checks of ``gdocs.mathview`` (notation, converters, prose, interactive payload) and of its assets.

Needs a licensed Symbolica (``SYMBOLICA_LICENSE`` in the environment) and ``hepkit_renorm`` on the path. Run from
``gamma5_docs/tools``::

    PYTHONPATH=../../python/hepkit_renorm python -m pytest tests/test_mathview.py
"""

from __future__ import annotations

import json
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parent.parent
SITE = TOOLS.parent / "site"
REPO = TOOLS.parent.parent
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(REPO / "python" / "hepkit_renorm"))

if not os.environ.get("SYMBOLICA_LICENSE"):
    pytest.skip("SYMBOLICA_LICENSE is not set", allow_module_level=True)

from gdocs import mathview as mv  # noqa: E402

MATH = "{http://www.w3.org/1998/Math/MathML}"


def _valid(html: str) -> ET.Element:
    root = ET.fromstring(html)
    assert root.tag == MATH + "math"
    assert "gdocs::" not in html and "hepkit_renorm" not in html
    assert "⋯" not in html, "Symbolica omitted a part of the expression"
    return root


@pytest.fixture(scope="module")
def store():
    from hepkit_renorm import builtin_model
    from hepkit_renorm.store import OrderStore

    model = builtin_model("sm_one_family")
    text = (SITE / "resources" / "run_sm_one_family" / "sm_one_family_store.json").read_text(encoding="utf-8")
    return model, OrderStore.from_json(model, text)


# ------------------------------------------------------------------------------------------------ notation
@pytest.mark.parametrize("entry", mv.registry(), ids=lambda e: e[0])
def test_registry_entry_prints_valid_mathml(entry):
    name, example, doc = entry
    m = mv.M(example)
    _valid(m["html"])
    assert m["text"].strip()
    # the source restores the displayed expression
    again = mv.M(mv.parse(m["src"]), lead_rationals=False)
    assert again["html"] == m["html"]


def test_glyphs_and_accents():
    eps = mv.M("ep")["html"]
    levi = mv.M("levi(mu, nu, rho, sigma)")["html"]
    assert "ϵ" in eps and "ε" not in eps          # regulator: lunate epsilon
    assert "ε" in levi and "ϵ" not in levi        # Levi-Civita: curly epsilon
    bar = _valid(mv.M("gammabar(mu)")["html"])
    assert bar.find(f".//{MATH}mover") is not None
    hat = mv.M("ghat(mu, nu)")["html"]
    assert "<mover" in hat and "ˆ" in hat
    assert "mi>mu<" not in mv.M("G(mu, A)")["html"]           # indices are glyphs, not identifiers
    assert "<msubsup>" in mv.M("G(mu, A)")["html"]
    assert "<msubsup>" in mv.M("g_s^2")["html"]               # power of a subscripted symbol


def test_leading_rationals_and_loop_factor():
    t = mv.M(mv.coefficient("-4/3*g_s^2", loop="pole"))["text"]
    assert t.startswith("-1/(16π²) 1/ϵ 4/3") or t.startswith("-1/(16π²) 1/ϵ (4/3")
    t2 = mv.M(mv.coefficient("-(27*g**2 + 48*g_s**2 + gp**2)/72"))["text"]
    assert "3/8 g²" in t2 and "2/3 g_s²" in t2
    assert "lf" in mv.M(mv.coefficient("3/(16*pi^2)*y_u^2", loop="auto"))["src"]


def test_foreign_symbols_raise():
    from symbolica import E

    with pytest.raises(mv.ConversionError):
        mv.M(E("foreign_ns::lor(1)"))
    with pytest.raises(ValueError):
        mv.parse("g_s +* )")


# ------------------------------------------------------------------------------------------------ converters
def test_operators_of_every_store_functional(store):
    model, st = store
    n = 0
    for entry in st:
        for k, f in entry.laurent.items():
            for mono in f.monomials():
                e = mv.operator(mono, table=f.fields)
                _valid(mv.M(e)["html"])
                n += 1
    assert n > 1000


def test_vertex_functions_of_every_class(store):
    model, st = store
    n = 0
    for name, power in (("S0", 0), ("Delta1", 0), ("deltaS1", -1), ("DeltaS1", 0)):
        f = st.entry(name).coefficient(power)
        for cls in f.leg_classes():
            info = {}
            v = mv.vertex(f, cls, model=model, info=info)
            _valid(mv.M(v)["html"])
            assert len(info["legs"]) == len(cls)
            for lg in info["legs"]:
                _valid(mv.M(lg)["html"])
            n += 1
    assert n > 200


def _class_op(f, cls):
    from hepkit_renorm.functional import Density

    return mv.operator(Density.from_monomials(f.fields, [m for m in f.monomials() if m.fields() == cls]))


def test_tree_operators_match_the_report(store):
    model, st = store
    S0 = st.entry("S0").coefficient(0)
    text = lambda cls: mv.M(_class_op(S0, cls))["text"]  # noqa: E731
    assert text(("Q", "Qbar")) == "𝑖 Q̅ γ^μ ∂_μ Q"                         # Ψ̄ i ∂̸ Ψ
    assert text(("G", "Q", "Qbar")) == "g_s Q̅ γ̅^μ T^A P_L Q G^A_μ"       # Ψ̄ γ̄ 𝒯 P_c Ψ A
    assert text(("Lbar", "E", "H")) == "-y_e L̅ H P_R E"                  # -y_e L̄ H P_R E
    assert text(("H", "Hdag")) == "-m² H^† H+(∂_μ H^†) ∂^μ H"
    assert text(("gh_G", "src_G")) == "ρ^μ_A ∂_μ ω^A"                    # ρ^μ_a ∂_μ ω^a
    assert text(("Q", "Ubar", "H")) == "y_u ε_{ij} U̅ P_L Q_j H_i"        # -y_u Ū H̃† P_L Q


def test_tree_vertex_functions_match_the_report(store):
    model, st = store
    S0 = st.entry("S0").coefficient(0)
    v = lambda cls: mv.M(mv.vertex(S0, cls, model=model))["text"]  # noqa: E731
    assert v(("G", "Q", "Qbar")) == "-g_s γ̅^μ T^A P_L"
    # the order of the two terms follows Symbolica's symbol order, which depends on what ran before in the process
    assert v(("H", "H", "Hdag", "Hdag")) in ("-2 λ (δ_{il} δ_{jk}+δ_{ik} δ_{jl})",
                                             "-2 λ (δ_{ik} δ_{jl}+δ_{il} δ_{jk})")
    three = v(("G", "G", "G"))
    assert three.startswith("𝑖 g_s f^{ACD} (") and "(p_1-p_2)^ρ" in three


# ------------------------------------------------------------------------------------------------ prose
def test_prose_spans_brackets_and_escaping():
    r = mv.prose("A <em>term</em> [[ Tr(T(A), comm(T(C), T(D))) ]] with a[1] and [[ x ]] & y < z.")
    assert r.count("<math") == 2
    assert "<em>term</em>" in r and "a[1]" in r and "&amp; y &lt; z" in r
    assert "[[" not in r
    assert "\\[[" not in mv.prose(r"literal \[[ stays") and "[[" in mv.prose(r"literal \[[ stays")
    r2 = mv.prose('<a data-term="unbroken:x" href="http://example.org" onclick="z()">t</a><script>1</script>')
    assert 'data-term="unbroken:x"' in r2 and "href" not in r2 and "onclick" not in r2
    assert "<script>" not in r2 and "&lt;script&gt;" in r2
    for bad in ("[[ x", "[[ x ] ]]", "[[ ]]", "[[ g_s +* ]]"):
        with pytest.raises(ValueError):
            mv.prose(bad)
    for frag in mv.prose("[[ f(a, b, c) ]] and [[ ord(Qbar, PL, Q) ]]").split("<math")[1:]:
        _valid("<math" + frag.split("</math>")[0] + "</math>")


# ------------------------------------------------------------------------------------------------ interactive payload
def test_view_payload_respects_the_budget(store):
    model, st = store
    dS = st.entry("deltaS1").coefficient(-1)
    exprs = [mv.coefficient("12*lambda^2 - 3*y_u^4 - 3*y_d^4 - y_e^4 + 9/16*g^4 + 3/8*g^2*gp^2 + 3/16*gp^4"),
             mv.vertex(st.entry("S0").coefficient(0), ("W", "W", "W", "W"), model=model),
             mv.vertex(dS, ("G", "Q", "Qbar"), model=model, loop="pole")]
    for e, budget, ps in zip(exprs, (12_000, 30_000, 30_000), (3, None, None)):
        m = mv.M(e, interactive=True, budget=budget, page_size=ps)
        v = m["view"]
        assert v is not None
        assert len(json.dumps(v, separators=(",", ":"), ensure_ascii=False)) <= budget + 64
        for key, ans in v["answers"].items():
            k = json.loads(key)
            assert k[0] in {"0"} | {a[2] for a in v["answers"].values() if a[0] == "p"}
            if ans[0] == "p":
                assert 0 <= ans[1] < len(v["pages"])
            else:
                assert ans[0] == "v" and 0 <= ans[1] < len(v["previews"])
        assert v["pages"][0]["start"] == 0
    # the paginated sum has answers for its marker (hover and pinned preview, inspection)
    m = mv.M(exprs[0], interactive=True, budget=30_000, page_size=3)
    actions = {json.loads(k)[1] for k in m["view"]["answers"]}
    assert {"preview", "open"} <= actions


def test_huge_expression_is_static(store):
    model, st = store
    big = mv.vertex(st.entry("deltaS1").coefficient(-1), ("G", "G", "G", "G"), model=model, loop="pole")
    m = mv.M(big, interactive=True, budget=10_000)
    assert m["view"] is None and m["html"].startswith("<math")


# ------------------------------------------------------------------------------------------------ assets
def test_assets_are_current():
    from symbolica import MathMLViewer

    d = SITE / "assets" / "viewer"
    assert (d / "symbolica_pager.js").read_text(encoding="utf-8") == MathMLViewer.javascript()
    font = (d / "STIXTwoMath-Regular.woff2").read_bytes()
    assert font[:4] == b"wOF2"
    js = (d / "mathview.js").read_text(encoding="utf-8")
    assert "export function renderMath" in js and "export function ready" in js
    assert "symbolica-notebook-math-font" in js     # pre-registered font face: the pager's CDN URL is never used
    assert "http://" not in js and "https://" not in js


def test_pole_factor_is_unambiguous():
    """Interface/physics review: 1/ϵ is a stacked fraction followed by a space, never '1/ϵ g²' in slash form."""
    html = mv.M(mv.coefficient("-g^2", loop="pole"))["html"]
    assert "<mfrac><mn>1</mn><mi>ϵ</mi></mfrac>" in html
    assert '<mo lspace="0" rspace="0">/</mo><mi>ϵ</mi>' not in html
