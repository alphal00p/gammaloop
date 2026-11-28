"""The order store: the local functionals that define the model at each loop order, tagged with their role.

At order ``l`` the recipe of report 02_method (sec:m-algorithm, sec:m-beyond) generates diagrams from a model that
contains, besides the tree-level action, the counterterms of all lower orders and the vertices of the breaking.
The store is the single place that holds these functionals; the model writer (:mod:`hepkit_renorm.output`) turns
its entries into vertex rules, and the order loop adds one set of entries per order. At one loop the store holds
only order 0 (``S_0`` and ``Delta_1``); the interfaces are those of every order.

Entries (:class:`StoreEntry`)
-----------------------------
An entry is a local functional with a loop order ``L``, a role (:class:`Role`) and a Laurent series in the
regulator ``epsilon`` (``d = 4 - 2 epsilon``)::

    C = sum_k epsilon^k C_k          (``laurent = {k: C_k}``)

Each coefficient ``C_k`` is an ``epsilon``-free local functional over the model's fields. It may contain the
dimension symbol ``d`` (the writer re-expands ``d = 4 - 2 epsilon``) and the evanescent tag
:data:`EVANESCENT_TAG` (report sec:m-beyond, bookkeeping B), polynomially. Evanescent structures (hatted
metrics, gammas, momenta) are not a role: they are a property of the coefficients, derived where needed.

Roles:

=================  ====================================================================================
``tree``           ``S_0`` (classical action, gauge fixing, ghosts and BRST sources); loop order 0
``breaking``       ghost-number-one insertions used only in Slavnov-Taylor diagrams: ``Delta_1`` at order 0,
                   ``b_d deltaS^(L)``, ``b_d DeltaS^(L)`` and ``S_q`` pieces at order ``L >= 1``
``ct_singular``    ``deltaS^(L)``: the singular (MS-bar) counterterms, invariant and evanescent; ``L >= 1``
``ct_finite``      ``DeltaS^(L)``: the finite Slavnov-Taylor-restoring counterterms; ``L >= 1``
=================  ====================================================================================

The vertices of an entry are its vertex functions (report 03b sec:dictionary): for a ``tree`` entry the classes
with at least three legs or with a BRST source (the other two-point classes are the quadratic form, inverted into
propagators); for every other role all classes, two-point ones included.

Serialisation: :meth:`OrderStore.to_dict` stores every coefficient as the canonical Symbolica string of its
density; :meth:`OrderStore.from_dict` restores it exactly over the fields of the model.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum

from symbolica import Expression, S

from .fields import FieldTable
from .functional import LocalFunctional
from .spec import Model

__all__ = ["EVANESCENT_TAG", "OrderStore", "Role", "StoreEntry", "tree_store"]

EVANESCENT_TAG = S("hepkit_renorm::eta", is_real=True)
"""The evanescent tag ``eta`` of report sec:m-beyond (kept symbolic in counterterm coefficients)."""


class Role(str, Enum):
    """Role of a store entry (module docstring)."""

    TREE = "tree"
    BREAKING = "breaking"
    CT_SINGULAR = "ct_singular"
    CT_FINITE = "ct_finite"


@dataclass(frozen=True)
class StoreEntry:
    """One functional of the store: ``sum_k epsilon^k laurent[k]`` of loop order ``loop_order`` and role ``role``.

    ``name`` identifies the entry (it becomes the ``operator`` of the vertices written from it), ``description``
    is free text (provenance). ``bookkeeping``: ``""`` for an entry of every bookkeeping of the higher orders,
    ``"A"`` or ``"B"`` for one that only that bookkeeping inserts (report sec:m-beyond, tab:m-bookkeeping:
    ``b_d deltaS^(L)`` and ``b_d DeltaS^(L)`` in A, the remainders ``hatDelta^(L)``, ``E^(L)`` in B).
    """

    name: str
    loop_order: int
    role: Role
    laurent: Mapping[int, LocalFunctional]
    description: str = ""
    bookkeeping: str = ""

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("a store entry needs a name")
        if self.loop_order < 0:
            raise ValueError("the loop order is a non-negative integer")
        role = Role(self.role)
        object.__setattr__(self, "role", role)
        if role is Role.TREE and self.loop_order != 0:
            raise ValueError("tree entries have loop order 0")
        if role in (Role.CT_SINGULAR, Role.CT_FINITE) and self.loop_order < 1:
            raise ValueError("counterterm entries have loop order >= 1")
        if role is Role.TREE and set(self.laurent) - {0}:
            raise ValueError("tree entries have no regulator dependence")
        if not self.laurent:
            raise ValueError("an entry needs at least one Laurent coefficient")
        if self.bookkeeping not in ("", "A", "B"):
            raise ValueError("bookkeeping is '', 'A' or 'B'")
        tables = {id(f.fields) for f in self.laurent.values()}
        if len(tables) != 1:
            raise ValueError("all Laurent coefficients of an entry must live over the same field table")
        object.__setattr__(self, "laurent", dict(sorted((int(k), v) for k, v in self.laurent.items())))

    @property
    def fields(self) -> FieldTable:
        return next(iter(self.laurent.values())).fields

    @property
    def powers(self) -> tuple[int, ...]:
        """The powers of ``epsilon`` with a coefficient, ascending."""
        return tuple(self.laurent)

    def coefficient(self, power: int) -> LocalFunctional:
        """``C_power`` (zero if absent)."""
        return self.laurent.get(power, LocalFunctional(self.fields))

    def leg_classes(self) -> list[tuple[str, ...]]:
        """Union of the leg classes of the coefficients, in canonical order."""
        found = {c for f in self.laurent.values() for c in f.leg_classes()}
        order = self.fields.order
        return sorted(found, key=lambda c: [order(x) for x in c])

    def to_dict(self) -> dict:
        out = {"name": self.name, "loop_order": self.loop_order, "role": self.role.value,
               "description": self.description,
               "laurent": {str(k): f.expr.to_canonical_string() for k, f in self.laurent.items()}}
        if self.bookkeeping:
            out["bookkeeping"] = self.bookkeeping
        return out

    @classmethod
    def from_dict(cls, fields: FieldTable, doc: Mapping) -> StoreEntry:
        # the strings are the entries' own canonical expressions: restored as they are (re-normalising a canonical
        # density need not be idempotent, see functional.Density)
        laurent = {int(k): LocalFunctional(fields, Expression.parse(v), _canonical=True)
                   for k, v in doc["laurent"].items()}
        return cls(doc["name"], int(doc["loop_order"]), Role(doc["role"]), laurent, doc.get("description", ""),
                   doc.get("bookkeeping", ""))


@dataclass
class OrderStore:
    """The entries of all orders computed so far, in insertion order (names are unique)."""

    fields: FieldTable
    _entries: dict[str, StoreEntry] = field(default_factory=dict)

    def add(self, entry: StoreEntry) -> None:
        if entry.name in self._entries:
            raise ValueError(f"duplicate store entry {entry.name!r}")
        if entry.fields is not self.fields:
            raise ValueError(f"entry {entry.name!r} lives over another field table")
        self._entries[entry.name] = entry

    def entry(self, name: str) -> StoreEntry:
        return self._entries[name]

    def entries(self, max_loop_order: int | None = None, roles: Iterable[Role] | None = None,
                bookkeeping: str | None = None) -> list[StoreEntry]:
        """Entries with ``loop_order <= max_loop_order`` (all if ``None``), a role in ``roles`` (all if ``None``)
        and, if ``bookkeeping`` is given, an entry of that bookkeeping or of both, in insertion order."""
        wanted = None if roles is None else {Role(r) for r in roles}
        return [e for e in self._entries.values()
                if (max_loop_order is None or e.loop_order <= max_loop_order) and (wanted is None or e.role in wanted)
                and (bookkeeping is None or e.bookkeeping in ("", bookkeeping))]

    def __iter__(self):
        return iter(self._entries.values())

    def __len__(self) -> int:
        return len(self._entries)

    def loop_orders(self) -> list[int]:
        return sorted({e.loop_order for e in self._entries.values()})

    def tree_action(self) -> LocalFunctional:
        """The sum of the ``tree`` entries: the functional whose quadratic part defines the propagators."""
        out = LocalFunctional(self.fields)
        for e in self.entries(roles=[Role.TREE]):
            out = LocalFunctional.integrate(out + e.coefficient(0))
        return out

    def to_dict(self) -> dict:
        return {"entries": [e.to_dict() for e in self._entries.values()]}

    @classmethod
    def from_dict(cls, model: Model, doc: Mapping) -> OrderStore:
        store = cls(model.fields)
        for d in doc["entries"]:
            store.add(StoreEntry.from_dict(model.fields, d))
        return store

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=1)

    @classmethod
    def from_json(cls, model: Model, text: str) -> OrderStore:
        return cls.from_dict(model, json.loads(text))


def tree_store(model: Model, brst=None) -> OrderStore:
    """The order-0 store: ``S0`` (role ``tree``) and ``Delta1`` = ``s S_0`` (role ``breaking``), from
    :class:`hepkit_renorm.brst.BRST` (``brst`` may be passed to reuse its cached objects)."""
    from .brst import BRST

    b = brst if brst is not None else BRST(model)
    store = OrderStore(model.fields)
    store.add(StoreEntry("S0", 0, Role.TREE, {0: b.action}, "tree-level action S_0 = classical + gauge fixing + "
                                                            "ghosts + BRST sources"))
    store.add(StoreEntry("Delta1", 0, Role.BREAKING, {0: b.breaking}, "tree-level breaking Delta_1 = s S_0"))
    return store
