"""Independent review: robustness of the pipeline and of the solver (MS6, MS7).

* Leg order (fast): the class results of the diagram pipeline and the fingerprints of local functionals do not
  depend on the order in which the legs of a class are listed, up to the Grassmann sign of the permutation of the
  odd legs (vertex-function convention, report 03_setup / 02_method eq:m-Gamma). Checked on classes with zero,
  two and three odd legs, with and without a ``Delta_1`` insertion.
* A deliberately wrong sign in one class result (slow, uses the per-graph cache ``HEPKIT_RENORM_CACHE`` if given):
  a flipped ghost-number-0 pole makes the singular counterterm fail its ``b``-closure (or the fit); a flipped
  ghost-number-1 pole fails the pole identity in exactly that class; a flipped ``P N1`` of the class ``(gh_G, G)``
  makes ``b X = 0`` fail. Limit (review finding): a flipped ``P N1`` along a direction that is ``b``-exact within
  one class (``N^(f r)_F = -b Q^(f r)_V`` in ``gh psibar psi``) passes every consistency check and only changes the
  restoring counterterm; it is caught only by an independent computation of the class (gate V6).
* A gauge anomaly (slow): the Standard Model without its leptons (``D_BBB``, ``D_BWW`` non-zero) run end to end
  from the specification stops at the restoring counterterm with ``AnomalyError``; the components of ``X_1`` outside
  ``b(span Q)`` lie on monomials with an ``epsilon`` tensor.

Run the slow part with ``pytest --runslow tests/test_independent_robustness.py`` (about 15 min on a loaded host with
a warm cache, ``HEPKIT_RENORM_PROCESSES`` workers).
"""

from __future__ import annotations

import itertools
import json
import os
import tomllib

import numpy as np
import pytest

from hepkit_renorm import builtin_model, model_from_dict
from hepkit_renorm.bases import breaking_basis
from hepkit_renorm.brst import BRST
from hepkit_renorm.class_results import ClassComputer, ResultCache
from hepkit_renorm.driver import RunOptions, renormalise
from hepkit_renorm.exact import default_field
from hepkit_renorm.fingerprint import Fingerprinter
from hepkit_renorm.functional import LocalFunctional
from hepkit_renorm.loopdata import ResolvedLoopData, localise
from hepkit_renorm.matching import (
    BREAKING,
    POLE,
    FunctionalData,
    MatchingError,
    OrderContext,
    match_breaking,
    match_singular,
    pole_identity,
)
from hepkit_renorm.st_solve import AnomalyError, solve_restoring
from hepkit_renorm.store import tree_store

pytestmark = pytest.mark.filterwarnings("ignore:This process .* is multi-threaded:DeprecationWarning")
PROCESSES = int(os.environ.get("HEPKIT_RENORM_PROCESSES", "4"))


@pytest.fixture(scope="module")
def setup(sm):
    brst = BRST(sm)
    store = tree_store(sm, brst)
    return sm, brst, store


# ---------------------------------------------------------------------------------------------- leg order
def _permuted(fp: Fingerprinter, fields, samples, perm):
    """The same configurations listed in the leg order ``fields[perm[j]]``, and the Grassmann sign."""
    new = [fields[k] for k in perm]
    moms = np.stack([samples.momenta[samples.legs[k].momentum] for k in perm], axis=1)
    pols = [samples.polarisations[k] for k in perm]
    odd = [k for k in perm if fp.model.fields.odd(fields[k])]
    inv = sum(1 for i, j in itertools.combinations(range(len(odd)), 2) if odd[i] > odd[j])
    return tuple(new), fp.explicit_samples(new, moms, pols), (-1) ** inv


@pytest.mark.parametrize("fields", [("G", "G", "G"), ("Qbar", "Q", "G"), ("gh_G", "Qbar", "Q"),
                                    ("gh_W", "H", "Hdag")])
def test_class_results_do_not_depend_on_the_leg_order(setup, fields):
    model, brst, store = setup
    fields = tuple(sorted(fields, key=model.fields.order))
    fp = Fingerprinter(model, default_field(), 4, seed=3, values={"xi": 1})
    comp = ClassComputer(model, store, 1, processes=PROCESSES)
    smp = fp.samples(list(fields), 3, salt=5)
    ref = comp.compute(fields, fp, smp)
    p = fp.field.p
    assert np.any(ref.N0 % p) or np.any(ref.N1 % p)
    tested = 0
    for perm in itertools.permutations(range(len(fields))):
        if list(perm) == sorted(perm):
            continue
        new, smp_p, sign = _permuted(fp, fields, smp, perm)       # identical fields: Bose / Fermi symmetry
        res = comp.compute(new, fp, smp_p)
        for a, b in ((res.N0, ref.N0), (res.N1, ref.N1)):
            assert np.array_equal(np.asarray(a) % p, np.asarray(sign * b) % p), (new, sign)
        # the same relation for the local functionals (the tree action, and Delta_1 at ghost number 1)
        for F in (brst.action, brst.breaking):
            assert np.array_equal(fp.functional(F, smp_p) % p, sign * fp.functional(F, smp) % p)
        tested += 1
        if tested == 2:
            break
    assert tested >= 1


# ---------------------------------------------------------------------------------------------- wrong signs
class _Flipped:
    """Loop data with the sign of one quantity flipped in one class (everything else delegated)."""

    def __init__(self, data: ResolvedLoopData, quantity: str, cls: tuple[str, ...]) -> None:
        self.data, self.quantity, self.cls = data, quantity, tuple(cls)

    def _sign(self, quantity, fields) -> int:
        return -1 if (quantity == self.quantity and tuple(fields) == self.cls) else 1

    def prefetch(self, requests):
        self.data.prefetch(requests)

    def values(self, quantity, fields, fp, samples):
        return self._sign(quantity, fields) * self.data.values(quantity, fields, fp, samples) % fp.field.p

    def parametric(self, quantity, fields, fp, space, samples):
        m = self.data.parametric(quantity, fields, fp, space, samples)
        return m.scale(self._sign(quantity, fields)) if self._sign(quantity, fields) < 0 else m


@pytest.fixture(scope="module")
def loop_setup(setup, tmp_path_factory):
    model, brst, store = setup
    ctx = OrderContext(model, brst=brst)
    cache_dir = os.environ.get("HEPKIT_RENORM_CACHE") or tmp_path_factory.mktemp("cache")
    computer = ClassComputer(model, store, 1)
    data = ResolvedLoopData(computer, ctx.space, processes=PROCESSES, cache=ResultCache(cache_dir))
    return model, brst, ctx, data


def _singular(ctx, brst, data):
    form = localise(ctx, data, POLE, ctx.singular.elements, ctx.classes(0), dimension=6, salt=51)
    C = LocalFunctional.integrate(form.functional(ctx).scale(-1))
    return match_singular(ctx, FunctionalData(brst, counterterm=C))


@pytest.mark.slow
def test_wrong_sign_in_a_ghost_number_0_pole_is_detected(loop_setup):
    _, brst, ctx, data = loop_setup
    assert _singular(ctx, brst, data).closure_rank == 84             # the correct data pass
    with pytest.raises(MatchingError):
        _singular(ctx, brst, _Flipped(data, POLE, ("G", "G")))


@pytest.mark.slow
def test_wrong_sign_in_a_ghost_number_1_pole_fails_the_pole_identity(loop_setup):
    model, brst, ctx, data = loop_setup
    singular = _singular(ctx, brst, data)
    target = ("gh_G", "Qbar", "Q")
    target = tuple(sorted(target, key=model.fields.order))
    others = [c for c in ctx.classes(1) if len(c) == 3][:4]
    classes = [target] + [c for c in others if c != target]
    good = pole_identity(ctx, singular, data, classes=classes)
    assert good.holds
    bad = pole_identity(ctx, singular, _Flipped(data, POLE, target), classes=classes)
    assert bad.failed == [target]


def _breaking(ctx, brst, data):
    form = localise(ctx, data, BREAKING, ctx.N.elements, ctx.classes(1), dimension=4, physical=True,
                    complex_values=True, salt=53)
    X = form.functional(ctx)
    return X, match_breaking(ctx, FunctionalData(brst, breaking=X))


@pytest.mark.slow
def test_wrong_sign_in_one_breaking_class_is_detected(loop_setup):
    """A flipped P N1 of the class (gh_G, G) breaks b X = 0 (N^((G))_3 is not b-closed without its partner
    N^((G))_(f 2) of the class (gh_G, G, G))."""
    model, brst, ctx, data = loop_setup
    target = tuple(sorted(("gh_G", "G"), key=model.fields.order))
    _, res = _breaking(ctx, brst, _Flipped(data, BREAKING, target))
    assert (not res.closed) or res.anomaly.anomalous


@pytest.mark.slow
def test_wrong_sign_along_a_b_exact_direction_is_not_detectable(loop_setup):
    """Limit of the consistency checks (review finding): in the class (gh_G, Qbar, Q) the breaking lies along
    N^(Q G)_F = -b Q^(Q G)_V, which is b-exact on its own. A flipped sign there leaves b X = 0 and the anomaly
    analysis intact; the solver returns a different, equally consistent restoring counterterm. Such an error is
    caught only by an independent computation of the class (gate V6, test_independent_v5_v8_run.py)."""
    model, brst, ctx, data = loop_setup
    target = tuple(sorted(("gh_G", "Qbar", "Q"), key=model.fields.order))
    X_bad, res = _breaking(ctx, brst, _Flipped(data, BREAKING, target))
    assert res.closed and not res.anomaly.anomalous
    X_good, _ = _breaking(ctx, brst, data)
    bad = solve_restoring(ctx, FunctionalData(brst, breaking=X_bad)).nonzero()
    good = solve_restoring(ctx, FunctionalData(brst, breaking=X_good)).nonzero()
    changed = {k for k in set(bad) | set(good) if bad.get(k) != good.get(k)}
    assert changed and all("[G Q Qbar]" in k or "Q Qbar" in k for k in changed), sorted(changed)


# ---------------------------------------------------------------------------------------------- anomaly
NO_LEPTONS = """name = "sm_no_leptons"
generations = 1
[gauge]
G = { group = "SU(3)", coupling = "g_s" }
W = { group = "SU(2)", coupling = "g" }
B = { group = "U(1)", coupling = "gp" }
[fermions]
Q = { chirality = "L", G = "F", W = "F", B = "1/6" }
U = { chirality = "R", G = "F", B = "2/3" }
D = { chirality = "R", G = "F", B = "-1/3" }
[scalars]
H = { complex = true, W = "F", B = "1/2" }
[parameters]
m2 = {}
lambda = {}
y_u = { flavour = true }
y_d = { flavour = true }
[[potential]]
coupling = "m2"
bilinears = [["H", "H"]]
[[potential]]
coupling = "lambda"
bilinears = [["H", "H"], ["H", "H"]]
[[yukawa]]
coupling = "y_u"
fermion_bar = "Q"
scalar = "Hdag"
fermion = "U"
[[yukawa]]
coupling = "y_d"
fermion_bar = "Q"
scalar = "H"
fermion = "D"
[gauge_fixing]
kind = "linear_covariant"
xi = "xi"
"""


@pytest.mark.slow
def test_standard_model_without_leptons_is_anomalous(tmp_path):
    model = model_from_dict(tomllib.loads(NO_LEPTONS))
    # the anomaly coefficients of the content (03b eq:Dabc): the hypercharge ones do not vanish
    D = model.anomaly_tensor(default_field())
    assert np.any(np.asarray(D, dtype=object) % default_field().p)
    out = tmp_path / "run"
    opts = RunOptions(output=out, cache=os.environ.get("HEPKIT_RENORM_CACHE_NOLEPTONS") or None,
                      processes=PROCESSES, pole_identity=False)
    with pytest.raises(AnomalyError) as info:
        renormalise(model, opts)
    comps = info.value.report.components()
    assert comps
    N = {e.name: e for e in breaking_basis_of(model)}
    assert all("eps(" in N[name].functional.expr.to_canonical_string() for name in comps), sorted(comps)
    report = json.loads((out / "run_report.json").read_text())
    assert report["status"] == "failed" and "AnomalyError" in report["error"]
    checks = {c["name"]: c["passed"] for c in report["checks"]}
    assert checks["no anomaly (X in b(span Q))"] is False
    assert checks["four-dimensional limit of deltaS is b-closed"] is True


def breaking_basis_of(model):
    return breaking_basis(model).elements


def test_builtin_model_still_anomaly_free():
    """Control of the anomaly test: the shipped one-family content has a vanishing anomaly tensor."""
    model = builtin_model("sm_one_family")
    D = model.anomaly_tensor(default_field())
    assert not np.any(np.asarray(D, dtype=object) % default_field().p)

