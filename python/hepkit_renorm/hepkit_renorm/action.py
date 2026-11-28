"""The regularised classical action built from a model specification.

``classical_action(model)`` is ``int d^d x (L_gauge + L_fermion + L_scalar - V + L_Yukawa)`` of report eq:L0
in the BMHV regularisation (report eq:Lgauge, eq:Lferm, eq:LHiggs, eq:LYuk):

* ``L_gauge = -1/4 F^a_{mu nu} F^{a mu nu}`` per factor, ``F^a_{mu nu} = d_mu A^a_nu - d_nu A^a_mu + g f^{abc} A^b_mu
  A^c_nu``, all indices d-dimensional;
* ``L_fermion = sum_f [Psibar_f i gamma^mu d_mu Psi_f + Psibar_f gbar^mu g_r T^a_f P_{c_f} Psi_f A^a_mu]``: the
  d-dimensional kinetic term of the Dirac field ``Psi_f = psi_f + chi_f`` (physical field plus sterile partner)
  and the barred gauge vertex of the physical chirality only;
* ``L_scalar = (D_mu phi)^dagger D^mu phi`` (``1/2 (D phi)^2`` for a real scalar), d-dimensional;
* ``V`` and the Yukawa terms as declared (``L`` contains ``-V``; Yukawa terms with their hermitian conjugates).

The gauge-fixing, ghost and BRST-source terms are not part of the classical action: they are derived from the
BRST structure (``L_gf + L_ext = s(gauge fermion)``) by :mod:`hepkit_renorm.brst`, which builds on
:func:`covariant_derivative` and :func:`field_strength`.

All builders return canonical densities/functionals over ``model.fields``.
"""

from __future__ import annotations

from fractions import Fraction

from symbolica import Expression

from . import tensors as tn
from .functional import Density, LocalFunctional, Occurrence
from .groups import generator_expr, structure_expr
from .spec import Model

__all__ = [
    "classical_action",
    "classical_sectors",
    "covariant_derivative",
    "fermion_terms",
    "field_strength",
    "gauge_kinetic",
    "potential",
    "scalar_kinetic",
    "yukawa_terms",
]


def _num(q: Fraction | int) -> Expression:
    q = Fraction(q)
    return Expression.num(q.numerator) / q.denominator


def covariant_derivative(model: Model, field: str, mu: Expression, labels: tuple[Expression, ...]) -> Density:
    """``(D_mu X)[labels]`` for any field ``X`` of the table: ``d_mu X - i sum_r g_r A^a_{r mu} (T^a_R X)``, with
    ``T_R`` the representation of each group slot of ``X`` (``"F"``, ``"Fbar"``, ``"A"``) and ``q`` for its U(1)
    charges. On a ghost this is ``(D_mu omega)^a = d_mu omega^a + fg^{abc} A^b_mu omega^c``. Free labels:
    ``(mu, *labels)``."""
    fields = model.fields
    info = fields[field]
    out = Density.monomial(fields, 1, Occurrence(field, labels, (mu,)))
    for k, slot in enumerate(info.slots):
        if slot.factor is None:
            continue
        fac = model.factors[slot.factor]
        a = tn.fresh(fac.adjoint_space())
        j = tn.fresh(slot.space)
        inner = list(labels)
        inner[k] = j
        coeff = -tn.I * model.coupling(slot.factor) * generator_expr(fac, slot.rep, a, labels[k], j)
        out = out + Density.monomial(fields, coeff, Occurrence(slot.factor, (mu, a)), Occurrence(field, tuple(inner)))
    for fac_name, q in info.charges.items():
        coeff = -tn.I * model.coupling(fac_name) * _num(q)
        out = out + Density.monomial(fields, coeff, Occurrence(fac_name, (mu,)), Occurrence(field, labels))
    return out


def field_strength(model: Model, factor: str, mu: Expression, nu: Expression, a: Expression | None) -> Density:
    """``F^a_{mu nu}`` of a gauge factor (``a`` is ``None`` for U(1)); free labels ``(mu, nu, a)``."""
    fields = model.fields
    fac = model.factors[factor]
    adj = () if fac.abelian else (a,)
    out = Density.monomial(fields, 1, Occurrence(factor, (nu,) + adj, (mu,))) - \
        Density.monomial(fields, 1, Occurrence(factor, (mu,) + adj, (nu,)))
    if not fac.abelian:
        b, c = tn.fresh(fac.adjoint_space()), tn.fresh(fac.adjoint_space())
        coeff = model.coupling(factor) * structure_expr(fac, a, b, c)
        out = out + Density.monomial(fields, coeff, Occurrence(factor, (mu, b)), Occurrence(factor, (nu, c)))
    return out


def gauge_kinetic(model: Model) -> LocalFunctional:
    """``int -1/4 F^a_{mu nu} F^{a mu nu}`` summed over the factors."""
    total = Density(model.fields)
    for name, fac in model.factors.items():
        mu, nu = tn.fresh(tn.LORENTZ), tn.fresh(tn.LORENTZ)
        a = None if fac.abelian else tn.fresh(fac.adjoint_space())
        F = field_strength(model, name, mu, nu, a)
        total = total + (F * F).scale(Expression.num(-1) / 4)
    return LocalFunctional.integrate(total)


def fermion_terms(model: Model) -> LocalFunctional:
    """d-dimensional kinetic terms and barred physical gauge vertices of all fermion multiplets."""
    fields = model.fields
    total = Density(fields)
    for name, fm in model.fermions.items():
        proj = tn.PL if fm.chirality == "L" else tn.PR
        bar = f"{name}bar"
        info = fields[name]
        mu = tn.fresh(tn.LORENTZ)
        s, t = tn.fresh(tn.SPINOR), tn.fresh(tn.SPINOR)
        rest = tuple(tn.fresh(sl.space) for sl in info.slots[1:])
        total = total + Density.monomial(fields, tn.I * tn.gamma(mu, s, t),
                                         Occurrence(bar, (s,) + rest), Occurrence(name, (t,) + rest, (mu,)))
        for k, slot in enumerate(info.slots):
            if slot.factor is None:
                continue
            fac = model.factors[slot.factor]
            a, j = tn.fresh(fac.adjoint_space()), tn.fresh(slot.space)
            inner = list((t,) + rest)
            inner[k] = j
            u = tn.fresh(tn.SPINOR)
            coeff = (tn.gammabar(mu, s, u) * proj(u, t) * model.coupling(slot.factor)
                     * generator_expr(fac, slot.rep, a, ((s,) + rest)[k], j))
            total = total + Density.monomial(fields, coeff, Occurrence(bar, (s,) + rest),
                                             Occurrence(name, tuple(inner)), Occurrence(slot.factor, (mu, a)))
        for fac_name, q in info.charges.items():
            u = tn.fresh(tn.SPINOR)
            coeff = tn.gammabar(mu, s, u) * proj(u, t) * model.coupling(fac_name) * _num(q)
            total = total + Density.monomial(fields, coeff, Occurrence(bar, (s,) + rest),
                                             Occurrence(name, (t,) + rest), Occurrence(fac_name, (mu,)))
    return LocalFunctional.integrate(total)


def scalar_kinetic(model: Model) -> LocalFunctional:
    """``int (D_mu phi)^dagger D^mu phi`` (complex) or ``int 1/2 D_mu phi D^mu phi`` (real), d-dimensional."""
    total = Density(model.fields)
    for name, sm in model.scalars.items():
        mu = tn.fresh(tn.LORENTZ)
        labels = model.slot_labels(name)
        Dphi = covariant_derivative(model, name, mu, labels)
        if sm.complex:
            total = total + Dphi.conjugate() * Dphi
        else:
            total = total + (Dphi * Dphi).scale(Expression.num(1) / 2)
    return LocalFunctional.integrate(total)


def _bilinear(model: Model, a: str, b: str) -> Density:
    """``a^dagger b`` with all slots contracted by deltas."""
    labels = model.slot_labels(b)
    conj = model.fields[a].conjugate
    return Density.monomial(model.fields, 1, Occurrence(conj, labels), Occurrence(b, labels))


def potential(model: Model) -> LocalFunctional:
    """``-int V``."""
    total = Density(model.fields)
    for term in model.potential:
        prod = Density.monomial(model.fields, 1)
        for a, b in term.bilinears:
            prod = prod * _bilinear(model, a, b)
        total = total + prod.scale(-_num(term.factor) * model.param(term.coupling))
    return LocalFunctional.integrate(total)


def yukawa_terms(model: Model) -> LocalFunctional:
    """``int sum_k [-y_k fbar phi P f + h.c.]`` with the derived group contractions (:meth:`Model.yukawa_structure`)."""
    fields = model.fields
    total = Density(fields)
    for y in model.yukawas:
        proj_name, contractions = model.yukawa_structure(y)
        names = (f"{y.fermion_bar}bar", y.scalar, y.fermion)
        labels = [list(model.slot_labels(n)) for n in names]
        s, t = labels[0][0], labels[2][0]
        coeff = -(tn.PL if proj_name == "PL" else tn.PR)(s, t)
        for _fac, kind, ((k1, j1), (k2, j2)) in contractions:
            l1, l2 = labels[k1][j1], labels[k2][j2]
            coeff = coeff * (tn.delta(l1, l2) if kind == "delta" else tn.epsf(l1, l2))
        fl_bar = _flavour_label(fields, names[0], labels[0])
        fl = _flavour_label(fields, names[2], labels[2])
        coeff = coeff * model.flavour_tensor(y.coupling, fl_bar, fl)
        term = Density.monomial(fields, coeff, *(Occurrence(n, tuple(l)) for n, l in zip(names, labels)))
        total = total + term + term.conjugate()
    return LocalFunctional.integrate(total)


def _flavour_label(fields, name: str, labels: list[Expression]) -> Expression | None:
    for k, slot in enumerate(fields[name].slots):
        if slot.space.kind == "flavour":
            return labels[k]
    return None


def classical_sectors(model: Model) -> dict[str, LocalFunctional]:
    """The sectors of the classical action: ``gauge``, ``fermion``, ``scalar``, ``potential``, ``yukawa``."""
    return {
        "gauge": gauge_kinetic(model),
        "fermion": fermion_terms(model),
        "scalar": scalar_kinetic(model),
        "potential": potential(model),
        "yukawa": yukawa_terms(model),
    }


def classical_action(model: Model) -> LocalFunctional:
    """The regularised classical action (sum of :func:`classical_sectors`)."""
    total = LocalFunctional(model.fields)
    for part in classical_sectors(model).values():
        total = LocalFunctional.integrate(total + part)
    return total
