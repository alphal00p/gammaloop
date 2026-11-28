"""Field registry: statistics, quantum numbers, index slots and conjugation of every field and source.

A :class:`FieldTable` is derived by :mod:`hepkit_renorm.spec` from a model specification; nothing in it is typed
by hand. :mod:`hepkit_renorm.functional` needs only the table (parities, slots, conjugates).

Naming scheme of the derived fields (``F`` a gauge factor, ``f`` a fermion multiplet, ``phi`` a scalar):

=================  =========================================================================================
``F``              gauge field ``A^a_mu`` of factor ``F`` (slots: Lorentz, adjoint)
``gh_F``           ghost ``omega^a``                    ``agh_F``  antighost ``omegabar^a``    ``b_F``  auxiliary ``b^a``
``f``, ``fbar``    Dirac field ``Psi_f = psi_f + chi_f`` and ``Psibar_f`` (slots: spinor, group, flavour)
``phi``            scalar; ``phidag`` its conjugate when complex
``src_X``          BRST source ``X^*`` paired with the field ``X`` (report eq:Phistar): ``src_F`` = rho,
                   ``src_gh_F`` = zeta, ``src_f`` = Rbar_f, ``src_fbar`` = R_f, ``src_H`` = sigma^dagger,
                   ``src_Hdag`` = sigma
=================  =========================================================================================

Index slots of dimension one are dropped: a U(1) gauge field has no adjoint slot, a singlet has no group slot
and with one generation there is no flavour slot.

Conjugation convention (``(X_1 ... X_n)^dagger = X_n^dagger ... X_1^dagger``): ``X^dagger = conj_sign *
conjugate``. Gauge fields, ghosts and ``b`` are real (+1), the antighost is imaginary (-1), ``Psi^dagger`` is
represented by ``Psibar`` (the gamma0 factors are absorbed by the Dirac adjoint of the coefficient, see
:func:`hepkit_renorm.tensors.conjugate`), the source of a real field is imaginary (-1) and the sources of a
complex pair are conjugate to each other (+1). These are the conventions of report 03b sec:Qbasis.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from fractions import Fraction

from .tensors import Space

__all__ = ["KINDS", "FieldInfo", "FieldTable", "Slot"]

KINDS = (
    "gauge", "ghost", "antighost", "auxiliary", "fermion", "fermion_bar", "scalar", "scalar_conj", "source",
)


@dataclass(frozen=True)
class Slot:
    """One index slot of a field: its space and, for group slots, the factor and representation with which the
    field transforms under that factor (``"F"``, ``"Fbar"``, ``"A"``)."""

    space: Space
    factor: str | None = None
    rep: str | None = None


@dataclass(frozen=True)
class FieldInfo:
    """Everything the algebra needs to know about one field or source.

    ``odd``: Grassmann parity. ``ghost_number`` and ``dimension`` (mass dimension at d = 4) as in report
    tab:quantum-numbers. ``slots``: index slots in the order Lorentz, spinor, group (factor order), flavour.
    ``charges``: U(1) charges (factor -> rational), already conjugated for conjugate fields and sources.
    ``conjugate``/``conj_sign``: ``X^dagger = conj_sign * conjugate``. ``factor``: gauge factor of gauge-sector
    fields. ``multiplet``: matter multiplet of matter fields and their sources. ``chirality``: physical
    chirality ``"L"``/``"R"`` of fermion multiplets. ``source_of``/``source``: the field of a source / the source
    of a field (``None`` if absent).
    """

    name: str
    kind: str
    odd: bool
    ghost_number: int
    dimension: Fraction
    slots: tuple[Slot, ...]
    charges: Mapping[str, Fraction] = field(default_factory=dict)
    conjugate: str = ""
    conj_sign: int = 1
    factor: str | None = None
    multiplet: str | None = None
    chirality: str | None = None
    source_of: str | None = None
    source: str | None = None

    @property
    def is_source(self) -> bool:
        return self.kind == "source"

    def lorentz_slots(self) -> list[int]:
        """Positions of the Lorentz slots."""
        return [i for i, s in enumerate(self.slots) if s.space.kind == "lorentz"]


class FieldTable(Mapping[str, FieldInfo]):
    """Ordered, immutable mapping ``name -> FieldInfo``. The order is the canonical order of occurrences."""

    def __init__(self, infos: list[FieldInfo]) -> None:
        self._infos = {fi.name: fi for fi in infos}
        if len(self._infos) != len(infos):
            raise ValueError("duplicate field names")
        self._order = {n: k for k, n in enumerate(self._infos)}
        for fi in infos:
            if fi.conjugate not in self._infos:
                raise ValueError(f"conjugate {fi.conjugate!r} of {fi.name!r} is not a field")

    def __getitem__(self, name: str) -> FieldInfo:
        return self._infos[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._infos)

    def __len__(self) -> int:
        return len(self._infos)

    def order(self, name: str) -> int:
        """Position of a field in the canonical order."""
        return self._order[name]

    def odd(self, name: str) -> bool:
        return self._infos[name].odd
