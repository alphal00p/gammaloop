"""hepkit_renorm: algebraic renormalisation of gauge theories with a non-anticommuting gamma5 (BMHV scheme).

Layers (each module documents its conventions):

* :mod:`~hepkit_renorm.exact`      prime fields, square roots, rational reconstruction
* :mod:`~hepkit_renorm.tensors`    index labels, tensor heads, BMHV simplifier, conjugation
* :mod:`~hepkit_renorm.groups`     SU(N) x U(1): explicit generators, invariants, invariant tensors
* :mod:`~hepkit_renorm.fields`     field registry (statistics, quantum numbers, slots, conjugation)
* :mod:`~hepkit_renorm.spec`       model specification and derived field content
* :mod:`~hepkit_renorm.functional` local functionals with Grassmann bookkeeping
* :mod:`~hepkit_renorm.action`     the regularised classical action from a specification
* :mod:`~hepkit_renorm.brst`       BRST operator, gauge fermion, Slavnov-Taylor operator, b_d, b, Delta_1
* :mod:`~hepkit_renorm.rules`      vertex rules, rule sets and propagators
* :mod:`~hepkit_renorm.fingerprint` exact plane-wave evaluation (fingerprints) over GF(p)
* :mod:`~hepkit_renorm.linalg`     exact linear algebra, rational reconstruction, interpolation
"""

from .action import classical_action, classical_sectors, covariant_derivative, field_strength
from .functional import Density, Derivation, FieldImage, Leg, LocalFunctional, Monomial, Occurrence, standard_legs
from .spec import Model, builtin_model, load_model, model_from_dict

__all__ = [
    "Density",
    "Derivation",
    "FieldImage",
    "Leg",
    "LocalFunctional",
    "Model",
    "Monomial",
    "Occurrence",
    "builtin_model",
    "classical_action",
    "classical_sectors",
    "covariant_derivative",
    "field_strength",
    "load_model",
    "model_from_dict",
    "standard_legs",
]
