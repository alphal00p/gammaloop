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
* :mod:`~hepkit_renorm.classes`    classes of external legs (any ghost number)
* :mod:`~hepkit_renorm.bases`      operator bases Q, N, the singular basis; :mod:`~hepkit_renorm.baction` the action of b
* :mod:`~hepkit_renorm.store`      the order store; :mod:`~hepkit_renorm.output` the model writer
* :mod:`~hepkit_renorm.graphs`, :mod:`~hepkit_renorm.diagrams`, :mod:`~hepkit_renorm.hepkit_graphs`  graphs
* :mod:`~hepkit_renorm.loops`      loop back-ends; :mod:`~hepkit_renorm.class_results` per-class results
* :mod:`~hepkit_renorm.parametric` parameter gradings and graded reconstruction
* :mod:`~hepkit_renorm.matching`, :mod:`~hepkit_renorm.st_solve`, :mod:`~hepkit_renorm.orders`  steps 5-10
* :mod:`~hepkit_renorm.loopdata`   class data from the graphs for the algebraic steps
* :mod:`~hepkit_renorm.driver`     the order loop from a specification to the renormalised model
  (command line: ``python -m hepkit_renorm``)
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
