"""BRST symmetry, gauge fixing, BRST sources and the Slavnov-Taylor operator, derived from a model.

Nothing in this module is typed per model: the BRST transformations follow from the gauge factors and the
representations of the index slots of every field (:class:`~hepkit_renorm.fields.FieldInfo`), the gauge-fixing and
source terms are ``s`` of a gauge fermion, the signs of the Slavnov-Taylor operator are read off the source terms,
and the linearised operators are built from functional derivatives of the action.

Conventions (report 03_setup sec:brst and 03b sec:st-signs; ``SCRATCHPAD/st/st_conventions.md``)
--------------------------------------------------------------------------------------------------
* ``s`` is Grassmann odd, raises the ghost number by one and acts as a *left* graded derivation
  (:class:`~hepkit_renorm.functional.Derivation`); sources are invariant. On the fields (report eq:brst,
  eq:brst-matter), with ``cal(T)^a = g_r T^a`` and ``fg^{abc} = g_r f^{abc}``::

      s A^a_mu = (D_mu omega)^a = d_mu omega^a + fg^{abc} A^b_mu omega^c
      s omega^a = -1/2 fg^{abc} omega^b omega^c          s omegabar^a = b^a        s b^a = 0
      s X = i omega^a (cal(T)^a_R X)                       for every matter field X, the ghost to the left

  ``T_R`` is the generator of the representation of each group slot of ``X`` (:func:`groups.generator_expr`:
  ``Fbar`` gives ``-T^T``, so ``s Psibar = i Psibar P cal(T) omega`` and ``s H^dagger = -i H^dagger cal(T) omega``)
  and the charge for a U(1) factor. With sterile partners only the physical chirality transforms:
  ``s Psi = i omega cal(T) P_c Psi`` and ``s Psibar = i Psibar P_cbar cal(T) omega``.
* Gauge fermion (report eq:gf-exact), Grassmann odd, ghost number -1 and anti-hermitian::

      Psi = sum_F omegabar^a (d^mu A^a_mu + xi/2 b^a) + sum_X kappa_X X^* X

  The source part runs over the primary fields ``X`` (gauge fields, ghosts, fermions, scalars) with
  ``kappa_X = (-1)^{|X^*|}``, so that ``s(kappa_X X^* X) = + X^* sX``; the conjugate partner of a complex field
  enters with minus the hermitian conjugate of the primary term. ``L_gf = s(gauge part)`` and
  ``L_ext = s(source part)``; this reproduces ``- rho A + zeta omega + Rbar Psi - Psibar R - sigma^dag H + H^dag sigma``
  of the report and its eq:Lext.
* Slavnov-Taylor operator (report eq:ST) of a Grassmann-even functional ``F``::

      ST(F) = int [ sum_pairs sign_X (d_side F / d X^*) (dL F / d X) + sum_Y (s Y) (dL F / d Y) ]

  The pairs ``(X, X^*)`` are the fields with a source; ``Y`` runs over the fields without a source whose
  ``s Y`` is not zero (the antighost, giving ``b^a dL F/d omegabar^a``). ``sign_X`` and ``side`` are *derived*:
  the derivative of ``L_ext`` with respect to ``X^*`` that reproduces ``s X`` (:class:`SourcePair`). A left
  derivative is used unless the source is odd and its left derivative gives ``-s X``; then the right derivative
  (``dR = -dL`` on an even functional) is written, as in the report (``sigma``).
* Linearised operator of an even ``F`` (report eq:b-sources, eq:b-action): ``ST(F + kappa H) = ST(F) +
  kappa ST_F H + O(kappa^2)`` with ``ST_F`` the left derivation ``X -> s_X dL F/dX^*``, ``X^* -> s_X dL F/dX``,
  ``Y -> s Y`` (``s_X`` the left-derivative sign of the pair). ``b_d = ST_{S_0}``; ``b`` is its four-dimensional
  limit. ``ST(S_0 + C) = Delta_1 + b_d C + S_q(C)`` (report eq:m-S-exact), with ``S_q`` the source-field part of
  ``ST``.
* ``Delta_1 = s S_0 = ST(S_0)`` (report eq:Delta1).

Zero tests (:func:`vanishes`) work on vertex functions: a local functional is zero up to integration by parts if
every vertex function is zero after the BMHV simplifier, or, failing that, after explicit evaluation of its group
tensors over a prime field (Lie-algebra identities such as ``[T^a, T^b] = i f^{abc} T^c``). A density ``D`` with
free labels is tested through the pairing ``int Z D`` with a field ``Z`` whose slots are contracted with the free
labels (:func:`paired`): it vanishes for all ``Z`` iff ``D`` vanishes.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import cached_property

from symbolica import Expression

from . import tensors as tn
from .action import classical_action, covariant_derivative
from .exact import PrimeField, default_field
from .functional import Density, Derivation, FieldImage, LocalFunctional, Occurrence, standard_legs
from .groups import generator_expr, structure_expr
from .spec import Model

__all__ = [
    "BRST",
    "PRIMARY_KINDS",
    "SourcePair",
    "four_dimensional",
    "paired",
    "vanishes",
]

PRIMARY_KINDS = ("gauge", "ghost", "fermion", "scalar")
"""Field kinds whose source term enters the gauge fermion as ``(-1)^{|X^*|} X^* X``; the conjugate kinds
(``fermion_bar``, ``scalar_conj``) get minus the hermitian conjugate of their partner's term."""


@dataclass(frozen=True)
class SourcePair:
    """A pair (field ``X``, source ``X^*``) of the Slavnov-Taylor operator.

    The term of ``ST(F)`` is ``sign * (d_side F / d X^*) (dL F / d X)``. ``left_sign`` is the equivalent sign with
    a left derivative on ``X^*`` (``s_Phi`` of report eq:m-glue): ``delta' F / delta X^* = left_sign * dL F/dX^*``.
    """

    field: str
    source: str
    sign: int
    side: str
    left_sign: int


class BRST:
    """The BRST structure of a model: ``s``, the gauge fermion, ``S_0``, the Slavnov-Taylor operator, its
    linearisations and the breaking ``Delta_1``.

    ``classical`` is the regularised classical action (default :func:`action.classical_action`); the tree-level
    action is ``S_0 = classical + L_gf + L_ext``. All results are canonical densities / local functionals over
    ``model.fields``; expensive objects are computed once.
    """

    def __init__(self, model: Model, classical: LocalFunctional | None = None) -> None:
        if model.gauge_fixing.kind != "linear_covariant":
            raise NotImplementedError(f"gauge fixing {model.gauge_fixing.kind!r}")
        self.model = model
        self.fields = model.fields
        self._classical = classical

    # ------------------------------------------------------------------ BRST transformations
    @cached_property
    def s(self) -> Derivation:
        """The BRST operator as a left graded derivation of parity 1 (images on every non-source field)."""
        images = {}
        for name, info in self.fields.items():
            if info.is_source:
                continue
            labels = self.model.slot_labels(name)
            image = self._image(name, labels)
            if not image.is_zero():
                images[name] = FieldImage(labels, image)
        return Derivation(1, images)

    def variation(self, name: str, labels: Sequence[Expression]) -> Density:
        """``s X`` of the field ``name`` with free ``labels`` in its slots (zero density for sources and ``b``)."""
        img = self.s.images.get(name)
        if img is None:
            return Density(self.fields)
        return img.density.rename(list(zip(img.labels, labels)))

    def _image(self, name: str, labels: tuple[Expression, ...]) -> Density:
        info, m = self.fields[name], self.model
        if info.kind == "gauge":
            mu, adj = labels[0], labels[1:]
            return covariant_derivative(m, f"gh_{info.factor}", mu, adj)
        if info.kind == "ghost":
            fac = m.factors[info.factor]
            if fac.abelian:
                return Density(self.fields)
            b, c = tn.fresh(fac.adjoint_space()), tn.fresh(fac.adjoint_space())
            coeff = -m.coupling(info.factor) * structure_expr(fac, labels[0], b, c) / 2
            return Density.monomial(self.fields, coeff, Occurrence(name, (b,)), Occurrence(name, (c,)))
        if info.kind == "antighost":
            return Density.monomial(self.fields, 1, Occurrence(f"b_{info.factor}", labels))
        if info.kind in ("fermion", "fermion_bar", "scalar", "scalar_conj"):
            return self._rotation(name, labels)
        return Density(self.fields)  # auxiliary field

    def _rotation(self, name: str, labels: tuple[Expression, ...]) -> Density:
        """``i omega^a (cal(T)^a_R X)`` on all group slots and U(1) charges of a matter field, with the projector on
        the physical chirality of a fermion (on the right of ``Psibar``)."""
        info, m, fields = self.fields[name], self.model, self.fields
        spinor = [k for k, sl in enumerate(info.slots) if sl.space.kind == "spinor"]
        proj = Expression.num(1)
        outer = list(labels)
        if info.kind in ("fermion", "fermion_bar") and info.chirality is not None:
            k = spinor[0]
            u = tn.fresh(tn.SPINOR)
            physical = info.chirality if info.kind == "fermion" else ("R" if info.chirality == "L" else "L")
            head = tn.PL if physical == "L" else tn.PR
            proj = head(labels[k], u) if info.kind == "fermion" else head(u, labels[k])
            outer[k] = u
        out = Density(fields)
        for k, slot in enumerate(info.slots):
            if slot.factor is None:
                continue
            fac = m.factors[slot.factor]
            a, j = tn.fresh(fac.adjoint_space()), tn.fresh(slot.space)
            inner = list(outer)
            inner[k] = j
            coeff = tn.I * m.coupling(slot.factor) * generator_expr(fac, slot.rep, a, labels[k], j) * proj
            out = out + Density.monomial(fields, coeff, Occurrence(f"gh_{slot.factor}", (a,)),
                                         Occurrence(name, tuple(inner)))
        for fac_name, q in info.charges.items():
            coeff = tn.I * m.coupling(fac_name) * Expression.num(q.numerator) / q.denominator * proj
            out = out + Density.monomial(fields, coeff, Occurrence(f"gh_{fac_name}", ()),
                                         Occurrence(name, tuple(outer)))
        return out

    # ------------------------------------------------------------------ gauge fermion and the tree-level action
    @cached_property
    def gauge_fermion_gf(self) -> LocalFunctional:
        """``int sum_F omegabar^a (d^mu A^a_mu + xi/2 b^a)``."""
        m, total = self.model, Density(self.fields)
        for name, fac in m.factors.items():
            mu = tn.fresh(tn.LORENTZ)
            adj = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
            gauge_function = Density.monomial(self.fields, 1, Occurrence(name, (mu,) + adj, (mu,))) + \
                Density.monomial(self.fields, m.gauge_fixing.xi / 2, Occurrence(f"b_{name}", adj))
            total = total + Density.monomial(self.fields, 1, Occurrence(f"agh_{name}", adj)) * gauge_function
        return LocalFunctional.integrate(total)

    @cached_property
    def gauge_fermion_sources(self) -> LocalFunctional:
        """``int sum_X kappa_X X^* X`` (module docstring): primary terms plus the anti-hermitian completion."""
        total = Density(self.fields)
        for name, info in self.fields.items():
            if info.source is None or info.kind not in PRIMARY_KINDS:
                continue
            labels = self.model.slot_labels(name)
            kappa = -1 if self.fields.odd(info.source) else 1
            term = Density.monomial(self.fields, kappa, Occurrence(info.source, labels), Occurrence(name, labels))
            total = total + term
            if info.conjugate != name:
                total = total - term.conjugate()
        return LocalFunctional.integrate(total)

    @cached_property
    def gauge_fermion(self) -> LocalFunctional:
        """The complete gauge fermion ``Psi`` (ghost number -1): ``L_gf + L_ext = s Psi``."""
        return LocalFunctional.integrate(self.gauge_fermion_gf + self.gauge_fermion_sources)

    @cached_property
    def gauge_fixing(self) -> LocalFunctional:
        """``L_gf`` of report eq:Lgf as ``s`` of the gauge part of the gauge fermion:
        ``b^a d.A^a + xi/2 b^a b^a - omegabar^a d^mu (D_mu omega)^a``."""
        return self.s(self.gauge_fermion_gf)

    @cached_property
    def external(self) -> LocalFunctional:
        """``L_ext`` of report eq:Lext: ``s`` of the source part of the gauge fermion."""
        return self.s(self.gauge_fermion_sources)

    @cached_property
    def classical(self) -> LocalFunctional:
        """The regularised classical action."""
        return self._classical if self._classical is not None else classical_action(self.model)

    @cached_property
    def action(self) -> LocalFunctional:
        """The tree-level action ``S_0 = classical + L_gf + L_ext`` (report eq:L0)."""
        return LocalFunctional.integrate(self.classical + self.gauge_fixing + self.external)

    @cached_property
    def breaking(self) -> LocalFunctional:
        """``Delta_1 = s S_0`` (report eq:Delta1); equal to ``ST(S_0)``."""
        return self.s(self.action)

    # ------------------------------------------------------------------ Slavnov-Taylor operator
    @cached_property
    def pairs(self) -> tuple[SourcePair, ...]:
        """The (field, source) pairs with the signs and derivative sides derived from ``L_ext``."""
        out = []
        simp = self.model.simplifier()
        for name, info in self.fields.items():
            if info.source is None:
                continue
            labels = self.model.slot_labels(name)
            derived = self.external.functional_derivative(info.source, labels).simplify(simp)
            target = self.variation(name, labels).simplify(simp)
            if derived == target:  # includes s X = 0 (an abelian ghost): the pair is kept with sign +1
                left = 1
            elif derived == -target:
                left = -1
            else:
                raise ValueError(f"the source term of {name} does not reproduce s {name}")
            if left == -1 and self.fields.odd(info.source):
                out.append(SourcePair(name, info.source, 1, "right", -1))
            else:
                out.append(SourcePair(name, info.source, left, "left", left))
        return tuple(out)

    @cached_property
    def linear_fields(self) -> tuple[str, ...]:
        """Fields without a source with a non-zero ``s`` (their ST term is ``int (sY) dL F/dY``)."""
        return tuple(n for n, img in self.s.images.items() if self.fields[n].source is None
                     and not img.density.is_zero())

    def slavnov_taylor(self, F: LocalFunctional, pairs: Sequence[SourcePair] | None = None) -> LocalFunctional:
        """``ST(F)`` (report eq:ST) of a Grassmann-even functional; ``pairs`` overrides the derived pairs (only for
        negative controls of the sign conventions)."""
        return LocalFunctional.integrate(self.quadratic_part(F, pairs) + self._linear_part(F))

    def quadratic_part(self, F: LocalFunctional, pairs: Sequence[SourcePair] | None = None) -> LocalFunctional:
        """``S_q(F)``: the source-field products of ``ST(F)`` (report sec:m-vi)."""
        if F.parity():
            raise ValueError("the Slavnov-Taylor operator is defined on Grassmann-even functionals")
        total = Density(self.fields)
        for pr in self.pairs if pairs is None else pairs:
            labels = self.model.slot_labels(pr.field)
            first = F.functional_derivative(pr.source, labels, pr.side)
            if first.is_zero():
                continue
            total = total + (first * F.functional_derivative(pr.field, labels)).scale(pr.sign)
        return LocalFunctional.integrate(total)

    def _linear_part(self, F: LocalFunctional) -> Density:
        total = Density(self.fields)
        for name in self.linear_fields:
            labels = self.model.slot_labels(name)
            total = total + self.variation(name, labels) * F.functional_derivative(name, labels)
        return total

    def linearised(self, F: LocalFunctional) -> Derivation:
        """``ST_F``: the left derivation with ``ST(F + kappa H) = ST(F) + kappa ST_F H + O(kappa^2)`` for an even
        ``F`` (module docstring). Its images on the fields are ``s_X dL F/dX^*``, on the sources ``s_X dL F/dX``."""
        if F.parity():
            raise ValueError("linearisation at a Grassmann-odd functional")
        images = {}
        for pr in self.pairs:
            for target, var in ((pr.field, pr.source), (pr.source, pr.field)):
                labels = self.model.slot_labels(target)
                image = F.functional_derivative(var, labels).scale(pr.left_sign)
                if not image.is_zero():
                    images[target] = FieldImage(labels, image)
        for name in self.linear_fields:
            images[name] = self.s.images[name]
        return Derivation(1, images)

    @cached_property
    def b_d(self) -> Derivation:
        """``b_d = ST_{S_0}`` in d dimensions (report eq:b-sources): ``s`` on the fields, the signed equations of
        motion of ``S_0`` on the sources."""
        return self.linearised(self.action)

    @cached_property
    def b(self) -> Derivation:
        """The four-dimensional limit of :attr:`b_d` (images projected with ``P`` and ``d -> 4``)."""
        return four_dimensional(self.b_d, self.model.simplifier())

    # ------------------------------------------------------------------ linear identities
    def gauge_condition(self, F: LocalFunctional, factor: str) -> Density:
        """``dF/db^a - (d^mu A^a_mu + xi b^a)`` (report eq:linear-ids, first identity), free label ``a``."""
        fac = self.model.factors[factor]
        adj = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        mu = tn.fresh(tn.LORENTZ)
        rhs = Density.monomial(self.fields, 1, Occurrence(factor, (mu,) + adj, (mu,))) + \
            Density.monomial(self.fields, self.model.gauge_fixing.xi, Occurrence(f"b_{factor}", adj))
        return F.functional_derivative(f"b_{factor}", adj) - rhs

    def ghost_equation(self, F: LocalFunctional, factor: str) -> Density:
        """``dF/d omegabar^a + d_mu dF/d rho^mu_a`` (report eq:linear-ids, second identity): zero iff ``F`` depends
        on ``omegabar`` and ``rho`` only through ``rhotilde = rho + d omegabar`` (eq:rhotilde). Free label ``a``."""
        fac = self.model.factors[factor]
        adj = () if fac.abelian else (tn.fresh(fac.adjoint_space()),)
        mu = tn.fresh(tn.LORENTZ)
        rho = self.fields[factor].source
        return F.functional_derivative(f"agh_{factor}", adj) + \
            F.functional_derivative(rho, (mu,) + adj).total_derivative(mu)


# ---------------------------------------------------------------------------------------------- helpers
def four_dimensional(derivation: Derivation, simplifier: tn.Simplifier) -> Derivation:
    """The derivation with every image replaced by its four-dimensional limit (``P``, then ``d -> 4``)."""
    images = {}
    for name, img in derivation.images.items():
        projected = img.density.four_dimensional_limit(simplifier)
        if not projected.is_zero():
            images[name] = FieldImage(img.labels, projected)
    return Derivation(derivation.parity, images)


def paired(model: Model, partner: str, labels: Sequence[Expression], density: Density) -> LocalFunctional:
    """``int partner[labels] * density``: the probe functional of a density whose free labels are ``labels``. It
    vanishes (:func:`vanishes`) iff the density does; the partner's slots must lie in the spaces of ``labels``."""
    return LocalFunctional.integrate(Density.monomial(model.fields, 1, Occurrence(partner, tuple(labels))) * density)


def vanishes(model: Model, functional: LocalFunctional, field: PrimeField | None = None,
             simplifier: tn.Simplifier | None = None) -> bool:
    """True if every vertex function of ``functional`` vanishes, after the BMHV simplifier or, if needed, after
    explicit evaluation of the group tensors over ``field`` (``Model.vanishes_modulo_group``). Exact; a ``False``
    can also mean that the structure needs an identity outside both (e.g. reordering of gamma matrices)."""
    simp = simplifier or model.simplifier()
    field = field or default_field()
    for cls in functional.leg_classes():
        legs = standard_legs(model.fields, cls)
        value = functional.vertex_function(legs, simp)
        if value.expand() == Expression.num(0):
            continue
        open_labels = [lab for leg in legs for lab, sl in zip(leg.labels, model.fields[leg.field].slots)
                       if sl.factor is not None]
        if not model.vanishes_modulo_group(value, open_labels, field):
            return False
    return True
