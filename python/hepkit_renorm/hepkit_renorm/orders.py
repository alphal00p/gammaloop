"""The order-``l`` step and what it leaves for order ``l + 1`` (report 02_method sec:m-algorithm steps 5-10,
sec:m-beyond 'What must be kept from lower orders').

:func:`run_order` chains the algebraic steps on the class data of order ``l`` (:mod:`hepkit_renorm.matching`,
:mod:`hepkit_renorm.st_solve`):

5. the ghost-number-0 poles matched onto the singular basis: ``deltaS^(l)`` with its invariant part
   (renormalisation constants) and evanescent part;
6. ``b_d deltaS^(l)`` (a vertex of the ghost-number-1 functional ``Delta^[l]``);
8. the pole identity (ghost-number-1 poles ``= -V_{b_d deltaS}``) and ``X_l`` in the basis ``N`` with ``b X_l = 0``
   and the anomaly analysis;
9. the restoring counterterm ``b DeltaS^(l) = -X_l`` with the normalisation conditions.

:func:`order_entries` turns the result into store entries (step 10), in their d-dimensional form, with the loop
factor ``(1/(16 pi^2))^l`` written out and ``pi`` symbolic:

==================  ==============  ============  =====================================================================
entry               role            bookkeeping   functional
==================  ==============  ============  =====================================================================
``deltaS<l>``       ``ct_singular``  both         ``(1/(16 pi^2))^l C / eps`` (``C`` on the singular basis, d-dimensional)
``DeltaS<l>``       ``ct_finite``    both         ``(1/(16 pi^2))^l sum r_alpha Q_alpha`` (four-dimensional = barred)
``bd_deltaS<l>``    ``breaking``     A            ``b_d`` of ``deltaS<l>``
``bd_DeltaS<l>``    ``breaking``     A            ``b_d`` of ``DeltaS<l>``
``hatDelta<l>``     ``breaking``     B            eq:m-LP: ``Evan[(1/(16 pi^2))^l sum N_1] + Evan b_d DeltaS^(l)``
==================  ==============  ============  =====================================================================

``[S_q(C)]^(L)`` is built at the next order from the counterterm entries of lower orders (it vanishes at
``l = 1``); ``E^(1) = 0`` (no entry). ``hatDelta<l>`` needs the evanescent part of the ``eta/eps`` coefficients of
the ghost-number-1 classes (``remainder``, a local functional in ``(1/(16 pi^2))^l`` units); without it only the
part ``Evan b_d DeltaS^(l)`` is known and the entry is not written.
"""

from __future__ import annotations

from dataclasses import dataclass

from symbolica import Expression

from .fingerprint import is_zero
from .functional import Density, LocalFunctional
from .matching import (
    BreakingResult,
    ClassData,
    OrderContext,
    PoleIdentity,
    SingularResult,
    feynman_gauge,
    match_breaking,
    match_singular,
    pole_identity,
)
from .output import WriterOptions, WrittenModel, write_model
from .st_solve import Normalisation, RestoringResult, solve_restoring
from .store import OrderStore, Role, StoreEntry

__all__ = ["OrderResult", "counterterm_model", "loop_factor", "order_entries", "prune", "run_order",
           "update_store"]


def loop_factor(order: int) -> Expression:
    """``(1/(16 pi^2))^order``."""
    return (Expression.num(16) * Expression.PI ** 2) ** (-order)


@dataclass
class OrderResult:
    """The algebraic results of order ``l`` (module docstring)."""

    order: int
    singular: SingularResult
    restoring: RestoringResult
    breaking: BreakingResult | None = None
    poles: PoleIdentity | None = None
    remainder: LocalFunctional | None = None


def run_order(ctx: OrderContext, poles: ClassData, breaking: ClassData, normalisation: Normalisation | None = None,
              *, analyse_breaking: bool = True, check_poles: bool = True) -> OrderResult:
    """Steps 5-9 of order ``ctx.order`` on the class data ``poles`` (ghost numbers 0 and 1) and ``breaking``."""
    singular = match_singular(ctx, poles)
    identity = pole_identity(ctx, singular, poles) if check_poles else None
    analysis = match_breaking(ctx, breaking) if analyse_breaking else None
    restoring = solve_restoring(ctx, breaking, normalisation)
    return OrderResult(ctx.order, singular, restoring, analysis, identity)


def prune(model, functional: LocalFunctional) -> LocalFunctional:
    """The functional without its vanishing classes (certified zero test :func:`hepkit_renorm.fingerprint.is_zero`
    class by class). The density of a local functional is not a canonical form: a class can be represented by
    terms that cancel (e.g. the auxiliary-field terms of ``b_d`` on the shifted gauge source,
    ``-d b + d b``), and a written model must not contain such classes."""
    keep = Density(model.fields)
    for cls in functional.leg_classes():
        part = LocalFunctional.integrate(Density.from_monomials(
            model.fields, [m for m in functional.monomials() if sorted(m.fields()) == sorted(cls)]))
        if not is_zero(model, part):
            keep = keep + part
    return LocalFunctional.integrate(keep)


def order_entries(result: OrderResult) -> list[StoreEntry]:
    """The store entries of order ``l`` (module docstring); vanishing classes are pruned (:func:`prune`)."""
    ctx = result.singular.context
    model = ctx.model
    l = result.order
    lf = loop_factor(l)
    brst = ctx.brst
    c = prune(model, result.singular.functional())
    ds = prune(model, result.restoring.functional())

    def scaled(f: LocalFunctional) -> LocalFunctional:
        return LocalFunctional.integrate(f.scale(lf))

    entries = [
        StoreEntry(f"deltaS{l}", l, Role.CT_SINGULAR, {-1: scaled(c)},
                   f"singular counterterm deltaS^({l}) (MS-bar poles, invariant and evanescent part)"),
        StoreEntry(f"DeltaS{l}", l, Role.CT_FINITE, {0: scaled(ds)},
                   f"finite Slavnov-Taylor-restoring counterterm DeltaS^({l}) (four-dimensional, barred)"),
        StoreEntry(f"bd_deltaS{l}", l, Role.BREAKING,
                   {-1: scaled(prune(model, LocalFunctional.integrate(brst.b_d(c))))},
                   f"b_d deltaS^({l}) (insertion of bookkeeping A)", "A"),
        StoreEntry(f"bd_DeltaS{l}", l, Role.BREAKING,
                   {0: scaled(prune(model, LocalFunctional.integrate(brst.b_d(ds))))},
                   f"b_d DeltaS^({l}) (insertion of bookkeeping A)", "A"),
    ]
    if result.remainder is not None:
        simp = model.simplifier()
        evan = LocalFunctional.integrate(brst.b_d(ds)).evanescent(simp)
        hat = prune(model, LocalFunctional.integrate(result.remainder + evan))
        entries.append(StoreEntry(f"hatDelta{l}", l, Role.BREAKING, {0: scaled(hat)},
                                  f"evanescent remainder hatDelta^({l}) (eq:m-LP; insertion of bookkeeping B)", "B"))
    return [e for e in entries if not all(f.is_zero() for f in e.laurent.values())]


def update_store(store: OrderStore, result: OrderResult) -> OrderStore:
    """Add the entries of order ``l`` to the store (in place) and return it."""
    for entry in order_entries(result):
        store.add(entry)
    return store


def counterterm_model(model, store: OrderStore, options: WriterOptions | None = None) -> WrittenModel:
    """The renormalised model of the store (tree level, counterterms of every order, breaking insertions) in
    Feynman gauge, through :func:`hepkit_renorm.output.write_model`."""
    if options is None:
        options = WriterOptions(substitutions=feynman_gauge(model))
    return write_model(model, store, options)

