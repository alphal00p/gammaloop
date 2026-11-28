"""The order loop of the recipe from a model specification to the renormalised model (report 02_method
sec:m-algorithm steps 1-10, sec:m-beyond), with a run report.

:func:`renormalise` takes a :class:`~hepkit_renorm.spec.Model` (or a specification file) and a target loop order
and runs, for ``l = 1 .. order``:

====  =====================  ===========================================================================================
step  name                   what is computed (module)
====  =====================  ===========================================================================================
1     ``action``             ``S_0``, ``s``, the sources, ``Delta_1 = s S_0`` (:mod:`~hepkit_renorm.brst`); checks
                             ``ST(S_0) = Delta_1`` and that ``Delta_1`` is evanescent
2     ``tree_model``         the tree-level regularised model file (:func:`~hepkit_renorm.output.tree_model`),
                             loaded in HepKIT
3     ``classes``            the leg classes of ghost number 0, 1, 2 (:mod:`~hepkit_renorm.classes`)
4     ``bases``              ``Q``, ``N`` and the d-dimensional singular basis (:mod:`~hepkit_renorm.bases`)
5     ``graphs``             the graphs of every class of ghost number 0 and 1 from the order store
                             (:mod:`~hepkit_renorm.diagrams`)
6     ``poles``              ghost-number-0 poles from the graphs (``D = 6``, full spinors), resolved in the parameters
                             and matched onto the singular basis (:func:`~hepkit_renorm.loopdata.localise`), then the
                             singular counterterm, renormalisation constants, evanescent part, closure
                             (:func:`~hepkit_renorm.matching.match_singular`); beta functions and anomalous dimensions
7     ``pole_identity``      ghost-number-1 poles from the graphs against ``-V_{b_d deltaS}`` on every class
                             (:func:`~hepkit_renorm.matching.pole_identity`)
8     ``breaking``           ``P N1`` of the ghost-number-1 classes (``D = 4``, physical spinors) matched onto ``N``:
                             ``X_l``; ``b X_l = 0`` and the anomaly analysis (:func:`~hepkit_renorm.matching.
                             match_breaking`)
9     ``restoring``          ``b DeltaS^(l) = -X_l`` with the normalisation conditions (:mod:`~hepkit_renorm.st_solve`)
10    ``store``, ``model``   the order store entries of order ``l`` (:mod:`~hepkit_renorm.orders`) and the renormalised
                             model file (:func:`~hepkit_renorm.orders.counterterm_model`), loaded in HepKIT
====  =====================  ===========================================================================================

Nothing in the loop is specific to a model or to ``l = 1``: the classes, bases, graphs, signs, supports and
normalisation sets are derived; the store holds the entries of every order (the families of the graphs of order
``l`` are the entries of loop order ``< l``). What exists only at one loop is reported by
:func:`missing_for_order` and raises :class:`OrderNotSupported` before any computation.

Outputs (in ``RunOptions.output``)
----------------------------------
* ``<name>_tree.json``: the tree-level regularised model (HepKIT JSON, Feynman gauge);
* ``<name>_renormalised.json``: the renormalised model (tree level, counterterms ``deltaS^(L)``, ``DeltaS^(L)``,
  breaking insertions; :mod:`~hepkit_renorm.output` describes the encoding);
* ``<name>_store.json``: the order store (exact, :meth:`~hepkit_renorm.store.OrderStore.from_json`);
* ``<name>_order_data.json``: per order the matched loop data as exact local functionals (:data:`ORDER_DATA`);
* ``run_report.json`` and ``run_report.txt``: what was computed, every check with its result, the results
  (constants, beta functions, evanescent coefficients, ``X_l`` on ``N``, ``r_alpha`` on ``Q``, existence
  conditions) and the timings (:class:`RunReport`).

The per-graph cache (:class:`~hepkit_renorm.class_results.ResultCache`) is ``RunOptions.cache`` (default
``<output>/graph_cache``): a rerun with the same model reuses every graph value.

Units and conventions of the results: the report's. Renormalisation constants ``delta Z_Phi`` and ``delta p``
without the factor ``(1/(16 pi^2))^l / eps`` (bare ``Phi_0 = Z^(1/2) Phi``, ``p_0 = mu^(w_p eps) (p + delta p)``,
``w_p`` the derived loop weight); ``beta_p = mu dp/dmu = sum_l 2 l delta p^(l) / (16 pi^2)^l`` and ``gamma_Phi =
(1/2) mu d ln Z_Phi / d mu = -sum_l l delta Z_Phi^(l) / (16 pi^2)^l`` (leading poles; the coefficients are given in
units of ``(1/(16 pi^2))^l``); ``X_l`` and ``r_alpha`` in units of ``(1/(16 pi^2))^l``.
"""

from __future__ import annotations

import json
import os
import platform
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from symbolica import Expression

from .brst import BRST
from .class_results import ClassComputer, ResultCache
from .classes import leg_classes
from .fingerprint import is_zero
from .functional import LocalFunctional
from .loopdata import ResolvedLoopData, localise, recorded
from .matching import (
    BREAKING,
    POLE,
    FunctionalData,
    MatchingError,
    OrderContext,
    feynman_gauge,
    match_breaking,
    match_singular,
    pole_identity,
)
from .orders import OrderResult, counterterm_model, update_store
from .output import WriterOptions, tree_model
from .parametric import Polynomial
from .spec import Model, builtin_model, load_model
from .st_solve import AnomalyError, ExistenceError, solve_restoring, standard_normalisation
from .store import OrderStore, tree_store

__all__ = ["ORDER_DATA", "Check", "OrderNotSupported", "RunOptions", "RunReport", "RunResult", "Step", "missing_for_order",
           "renormalise", "resolve_model"]


ORDER_DATA = ("per order l: C with deltaS^(l) = (1/(16 pi^2))^l C / eps (from the ghost-number-0 poles) and X_l in "
              "units of (1/(16 pi^2))^l (from P N1 of the ghost-number-1 classes), as canonical Symbolica strings of "
              "their densities: LocalFunctional(model.fields, Expression.parse(s), _canonical=True)")


class OrderNotSupported(NotImplementedError):
    """The target order needs pieces that are not implemented (:func:`missing_for_order`)."""


def missing_for_order(order: int) -> list[str]:
    """What the order loop lacks at loop order ``order`` (empty for ``order = 1``)."""
    if order <= 1:
        return []
    return [
        ("a loop back-end for L >= 2 loops (OneLoopBackend refuses L != 1; the LoopBackend interface and the "
         "l-loop graph enumeration with counterterm insertions exist)"),
        ("the subtraction of subdivergences (forest of divergent subgraphs, loops.divergent_subgraphs, and the "
         "R-bar operation are not applied)"),
        ("the readout of graphs with L < l loops (counterterm insertions) in units (i/(16 pi^2))^l and of the "
         "general Laurent coefficients c[r, k] (ResolvedLoopData and LoopData read N0 = c[0, 1], N1 = c[1, 1])"),
        ("the quadratic Slavnov-Taylor pieces [S_q(C)]^(L) of the order-l breaking Delta^[l] (not built into the "
         "store)"),
        ("the bookkeeping-B remainder hatDelta^(L) on every class (matching.match_remainder exists, is run on "
         "request only and was validated on two classes)"),
    ]


# ---------------------------------------------------------------------------------------------- options and report
@dataclass
class RunOptions:
    """Options of :func:`renormalise`.

    ``output``: directory of the outputs; ``cache``: per-graph cache directory (default ``<output>/graph_cache``,
    ``False`` for none); ``processes``: forked workers for the graph evaluation; ``normalisation``: ``"full"`` or
    ``"minimal"`` (:func:`~hepkit_renorm.st_solve.standard_normalisation`); ``fermion_kinetic``: the
    representative with fermion kinetic counterterms; ``fallback``: relax a normalisation condition that admits no
    solution; ``pole_identity``: check the ghost-number-1 poles (the most expensive loop pass); ``cohomology``: the
    dimension of the cohomology of ``b`` on ``N`` (slow); ``schema``: of the written model; ``seed``: of every
    random choice (samples, points); ``name``: base name of the files (default: the model name)."""

    output: Path | str = "renormalised"
    order: int = 1
    cache: Path | str | bool | None = None
    processes: int = field(default_factory=lambda: os.cpu_count() or 1)
    normalisation: str = "full"
    fermion_kinetic: bool = False
    fallback: bool = True
    pole_identity: bool = True
    cohomology: bool = False
    schema: str = "v1"
    seed: int = 0
    name: str | None = None
    log: Callable[[str], None] | None = None


@dataclass
class Check:
    """One check: ``passed`` with a short ``detail`` (numbers compared, counts)."""

    name: str
    passed: bool
    detail: str = ""
    order: int = 0


@dataclass
class Step:
    """One step of the run: ``seconds`` and a summary of what it computed."""

    name: str
    seconds: float
    summary: dict[str, Any] = field(default_factory=dict)
    order: int = 0


@dataclass
class RunReport:
    """What was computed, every check and its result, the results and the timings (module docstring)."""

    model: str
    order: int
    options: dict[str, Any]
    steps: list[Step] = field(default_factory=list)
    checks: list[Check] = field(default_factory=list)
    results: dict[str, Any] = field(default_factory=dict)
    outputs: dict[str, str] = field(default_factory=dict)
    status: str = "running"
    error: str = ""
    environment: dict[str, str] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return self.status == "done" and all(c.passed for c in self.checks)

    def to_dict(self) -> dict[str, Any]:
        return {"model": self.model, "order": self.order, "status": self.status, "error": self.error,
                "all_checks_passed": self.passed, "options": self.options, "environment": self.environment,
                "outputs": self.outputs, "steps": [asdict(s) for s in self.steps],
                "checks": [asdict(c) for c in self.checks], "results": self.results}

    def text(self) -> str:
        """A readable summary (the JSON report holds everything)."""
        lines = [f"hepkit_renorm run report: model {self.model}, target order {self.order}",
                 f"status: {self.status}" + (f" ({self.error})" if self.error else ""),
                 f"checks: {sum(c.passed for c in self.checks)} of {len(self.checks)} passed", "", "steps:"]
        for s in self.steps:
            lines.append(f"  [l={s.order}] {s.name:<16} {s.seconds:9.1f} s  " +
                         ", ".join(f"{k}={v}" for k, v in s.summary.items() if not isinstance(v, (dict, list))))
        lines += ["", "checks:"]
        for c in self.checks:
            lines.append(f"  [{'PASS' if c.passed else 'FAIL'}] [l={c.order}] {c.name}: {c.detail}")
        for l_key, res in self.results.items():
            if not isinstance(res, dict):
                continue
            lines += ["", f"results of order {l_key}:"]
            for title, key in (("renormalisation constants (units (1/(16 pi^2))^l / eps)", "constants"),
                               ("beta functions (units (1/(16 pi^2))^l)", "beta"),
                               ("anomalous dimensions (units (1/(16 pi^2))^l)", "gamma"),
                               ("evanescent counterterm coefficients [b_d-closed]", "evanescent"),
                               ("restoring counterterm r_alpha (non-zero)", "restoring"),
                               ("existence conditions", "existence")):
                block = res.get(key)
                if not block:
                    continue
                lines.append(f"  {title}:")
                items = block.items() if isinstance(block, dict) else enumerate(block)
                for k, v in items:
                    lines.append(f"    {k}: {v}")
            if "breaking" in res:
                lines.append(f"  breaking X_l: {res['breaking'].get('nonzero', '?')} non-zero components on N")
        lines += ["", "outputs:"] + [f"  {k}: {v}" for k, v in self.outputs.items()]
        return "\n".join(lines) + "\n"

    def write(self, directory: Path) -> None:
        (directory / "run_report.json").write_text(json.dumps(self.to_dict(), indent=1, default=str))
        (directory / "run_report.txt").write_text(self.text())


@dataclass
class RunResult:
    """The outcome: the model, the order store, the results of every order, the written model and the report."""

    model: Model
    brst: BRST
    store: OrderStore
    orders: list[OrderResult]
    contexts: list[OrderContext]
    written: object
    report: RunReport
    loop_data: list[ResolvedLoopData]


def resolve_model(spec: Model | str | os.PathLike) -> Model:
    """A model from a :class:`Model`, a specification file (``.toml``) or the name of a built-in specification."""
    if isinstance(spec, Model):
        return spec
    path = Path(spec)
    if path.suffix == ".toml" or path.exists():
        return load_model(path)
    return builtin_model(str(spec))


# ---------------------------------------------------------------------------------------------- formatting
def _expr(ctx: OrderContext, poly: Polynomial) -> str:
    return str(ctx.space.expression(poly)) if poly else "0"


def _class_name(c) -> str:
    return " ".join(c)


# ---------------------------------------------------------------------------------------------- the run
class _Run:
    def __init__(self, model: Model, options: RunOptions) -> None:
        self.model = model
        self.opt = options
        self.out = Path(options.output)
        self.out.mkdir(parents=True, exist_ok=True)
        self.name = options.name or model.name
        cache = options.cache
        if cache is None:
            cache = self.out / "graph_cache"
        self.cache = ResultCache(cache) if cache else None
        opts = {k: str(v) for k, v in asdict(options).items() if k != "log"}
        self.report = RunReport(model.name, options.order, opts, environment={
            "python": platform.python_version(), "host": platform.node(), "processes": str(options.processes)})
        self.order = 0
        self.order_data: list[dict[str, Any]] = []

    # ------------------------------------------------------------------ bookkeeping
    def log(self, text: str) -> None:
        if self.opt.log is not None:
            self.opt.log(text)

    def step(self, name: str, fn: Callable[[], dict[str, Any]]):
        self.log(f"[l={self.order}] {name} ...")
        t0 = time.time()
        summary = fn() or {}
        dt = time.time() - t0
        self.report.steps.append(Step(name, round(dt, 2), summary, self.order))
        self.log(f"[l={self.order}] {name} done in {dt:.1f} s")
        self.report.write(self.out)
        return summary

    def check(self, name: str, passed: bool, detail: str = "") -> bool:
        self.report.checks.append(Check(name, bool(passed), detail, self.order))
        self.log(f"    check {name}: {'passed' if passed else 'FAILED'} {detail}")
        return bool(passed)

    def results(self, l: int) -> dict[str, Any]:
        return self.report.results.setdefault(str(l), {})

    # ------------------------------------------------------------------ steps of order 0
    def tree(self) -> None:
        model = self.model

        def action():
            self.brst = BRST(model)
            S0, D1 = self.brst.action, self.brst.breaking
            st = LocalFunctional.integrate(self.brst.slavnov_taylor(S0) - D1)
            self.check("ST(S_0) = Delta_1", is_zero(model, st), "exact zero test (D = 4, 6, 8) of ST(S_0) - s S_0")
            self.check("Delta_1 is non-zero and evanescent",
                       not D1.is_zero() and is_zero(model, LocalFunctional.integrate(
                           D1.four_dimensional_limit(model.simplifier()))),
                       f"{len(D1.leg_classes())} classes; its four-dimensional limit vanishes")
            return {"fields": len(model.fields), "S0_classes": len(S0.leg_classes()),
                    "Delta1_classes": len(D1.leg_classes())}

        self.step("action", action)

        def tree_file():
            self.store = tree_store(model, self.brst)
            written = tree_model(model, WriterOptions(schema=self.opt.schema, substitutions=self._gauge()),
                                 brst=self.brst)
            path = written.write(self.out / f"{self.name}_tree.json")
            loaded = written.load()
            self.report.outputs["tree_model"] = str(path)
            d = written.data
            self.check("tree model loads in HepKIT", len(loaded.vertex_rules) == len(d["vertex_rules"]),
                       f"{len(d['particles'])} particles, {len(d['vertex_rules'])} vertex rules")
            return {"particles": len(d["particles"]), "vertex_rules": len(d["vertex_rules"]),
                    "couplings": len(d["couplings"])}

        self.step("tree_model", tree_file)

    def _gauge(self) -> dict[str, int]:
        return feynman_gauge(self.model)

    # ------------------------------------------------------------------ order l
    def run_order(self, l: int) -> tuple[OrderResult, OrderContext, ResolvedLoopData]:
        model, brst = self.model, self.brst
        self.order = l
        res = self.results(l)
        ctx = OrderContext(model, l, brst=brst, seed=self.opt.seed)

        def classes():
            out = {}
            for g in (0, 1, 2):
                cl = leg_classes(model, g)
                out[f"ghost_number_{g}"] = len(cl)
                out[f"structures_{g}"] = sum(c.n_structures for c in cl)
            self.check("classes derived", out["ghost_number_0"] > 0 and out["ghost_number_1"] > 0,
                       f"{out['ghost_number_0']} / {out['ghost_number_1']} / {out['ghost_number_2']} classes of "
                       "ghost number 0 / 1 / 2")
            return out

        self.step("classes", classes)

        def bases():
            out = {"Q": len(ctx.Q), "N": len(ctx.N), "singular": len(ctx.singular.elements),
                   "evanescent": len(ctx.singular.evanescent)}
            self.check("Q and N complete", len(ctx.Q) == sum(c.n_structures for c in leg_classes(model, 0))
                       and len(ctx.N) == sum(c.n_structures for c in leg_classes(model, 1)),
                       f"|Q| = {len(ctx.Q)}, |N| = {len(ctx.N)} equal the structure counts of their classes")
            return out

        self.step("bases", bases)

        computer = ClassComputer(model, self.store, l)
        data = ResolvedLoopData(computer, ctx.space, processes=self.opt.processes, cache=self.cache,
                                seed=self.opt.seed)
        g0, g1 = ctx.classes(0), ctx.classes(1)

        def graphs():
            per = {}
            vanishing = 0
            for c in g0 + g1:
                enum, recs, _, _ = data.graphs(c)
                per[_class_name(c)] = len(recs)
                vanishing += len(enum.vanishing)
            n0 = sum(per[_class_name(c)] for c in g0)
            n1 = sum(per[_class_name(c)] for c in g1)
            res["graphs_per_class"] = per
            return {"graphs_ghost_number_0": n0, "graphs_ghost_number_1": n1, "vanishing_odd_automorphism": vanishing,
                    "families": [f.name for f in computer.families]}

        self.step("graphs", graphs)

        # 5-6: poles and the singular counterterm
        state: dict[str, Any] = {}

        def poles():
            form = localise(ctx, data, POLE, ctx.singular.elements, g0, dimension=6, salt=51)
            # the poles are N0 = -V_C with deltaS^(l) = (1/(16 pi^2))^l C / eps (matching.FunctionalData)
            C = LocalFunctional.integrate(form.functional(ctx).scale(-1))
            state["C"] = C
            self.check("ghost-number-0 poles in the span of the singular basis", True,
                       f"{form.rows} exact equations on {len(g0)} classes (D = 6, full spinors), coefficients "
                       f"polynomials on their graded supports; reconstructed at {len(form.primes)} prime(s), "
                       "confirmed at a further prime")
            try:
                singular = match_singular(ctx, FunctionalData(brst, counterterm=C))
            except MatchingError as exc:
                self.check("four-dimensional limit of deltaS is b-closed and spanned by the directions", False,
                           str(exc))
                raise
            self.check("four-dimensional limit of deltaS is b-closed", True,
                       f"all {singular.closure_rank} relations of b on the Q coefficients hold")
            self.check("deltaS spanned by field and parameter directions", not singular.unmatched_directions,
                       f"{len(singular.directions)} derived directions; unmatched: {singular.unmatched_directions}")
            state["singular"] = singular
            res["constants"] = {k: _expr(ctx, v) for k, v in singular.constants.items()}
            res["beta"], res["gamma"] = self._beta_gamma(ctx, singular, l)
            res["evanescent"] = {k: f"{_expr(ctx, v)} [{'closed' if cl else 'not closed'}]"
                                 for k, (v, cl) in singular.evanescent.items()}
            res["evanescent_operators"] = {e.name: e.origin for e in ctx.singular.evanescent}
            res["singular_coefficients"] = {k: _expr(ctx, v) for k, v in singular.coefficients.items() if v}
            return {"nonzero_singular_coefficients": len(res["singular_coefficients"]),
                    "evanescent_nonzero": len(singular.evanescent), "closure_rank": singular.closure_rank,
                    "loop_evaluations": data.stats.evaluations, "cache_hits": data.stats.cached}

        self.step("poles", poles)
        singular = state["singular"]

        identity = None
        if self.opt.pole_identity:
            def poles1():
                nonlocal identity
                data.prefetch(recorded(lambda d: pole_identity(ctx, singular, d), ctx.space))
                identity = pole_identity(ctx, singular, data)
                self.check("pole identity (ghost-number-1 poles = -V_{b_d deltaS})", identity.holds,
                           f"{len(identity.classes) - len(identity.failed)} of {len(identity.classes)} classes "
                           f"(D = {identity.dimension}, {identity.points} parameter points); failed: "
                           f"{[_class_name(c) for c in identity.failed]}")
                return {"classes": len(identity.classes), "failed": len(identity.failed),
                        "loop_evaluations": data.stats.evaluations, "cache_hits": data.stats.cached}

            self.step("pole_identity", poles1)
        else:
            self.report.steps.append(Step("pole_identity", 0.0, {"skipped": "RunOptions.pole_identity = False"}, l))

        def breaking():
            form = localise(ctx, data, BREAKING, ctx.N.elements, g1, dimension=4, physical=True,
                            complex_values=True, salt=53)
            X = form.functional(ctx)
            state["X"] = X
            self.check("P N1 of the ghost-number-1 classes in the span of N", True,
                       f"{form.rows} exact equations on {len(g1)} classes (D = 4, physical spinors, both "
                       f"embeddings of i), coefficients polynomials on their graded supports; reconstructed at "
                       f"{len(form.primes)} prime(s), confirmed at a further prime")
            analysis = match_breaking(ctx, FunctionalData(brst, breaking=X), cohomology=self.opt.cohomology)
            state["analysis"] = analysis
            self.check("b X = 0", analysis.closed, f"on {analysis.closure_rows} rows of the ghost-number-2 classes")
            comps = analysis.anomaly.components()
            self.check("no anomaly (X in b(span Q))", not analysis.anomaly.anomalous,
                       f"rank b(span Q) = {analysis.anomaly.image_rank}; components outside: {comps}"
                       + (f"; cohomology dimension {analysis.anomaly.cohomology}"
                          if analysis.anomaly.cohomology is not None else ""))
            nz = form.nonzero()
            res["breaking"] = {"nonzero": len(nz), "coefficients": {k: _expr(ctx, v) for k, v in nz.items()},
                               "anomaly_components": {k: str(v) for k, v in comps.items()},
                               "image_rank": analysis.anomaly.image_rank}
            return {"nonzero_components": len(nz), "closed": analysis.closed,
                    "anomalous": analysis.anomaly.anomalous, "loop_evaluations": data.stats.evaluations,
                    "cache_hits": data.stats.cached}

        self.step("breaking", breaking)

        def restoring():
            norm = standard_normalisation(model, mode=self.opt.normalisation, fallback=self.opt.fallback,
                                         fermion_kinetic=self.opt.fermion_kinetic)
            try:
                sol = solve_restoring(ctx, FunctionalData(brst, breaking=state["X"]), norm)
            except AnomalyError as exc:
                self.check("restoring counterterm exists", False, f"anomaly: {exc}")
                raise
            except ExistenceError as exc:
                self.check("restoring counterterm exists", False, str(exc))
                raise
            state["restoring"] = sol
            relaxed = [c.name for c in sol.conditions if c.relaxed]
            self.check("restoring counterterm unique", sol.unique,
                       f"{sol.unknowns} unknowns, {sol.allowed} allowed by the normalisation, rank {sol.rank_allowed}; "
                       f"rank of b on Q {sol.rank_all}")
            self.check("normalisation conditions applied without relaxation", not relaxed,
                       f"relaxed: {relaxed}" if relaxed else
                       ", ".join(f"{c.name}: {len(c.applied)}/{c.rows}" for c in sol.conditions))
            res["restoring"] = {k: _expr(ctx, v) for k, v in sol.nonzero().items()}
            res["existence"] = [f"{n} (testable in {[_class_name(c) for c in cls]})" for n, cls in sol.existence]
            res["normalisation"] = [asdict(c) for c in sol.conditions]
            return {"unknowns": sol.unknowns, "allowed": sol.allowed, "rank": sol.rank_allowed,
                    "nonzero": len(sol.nonzero()), "existence_conditions": sol.existence_total}

        self.step("restoring", restoring)
        result = OrderResult(l, singular, state["restoring"], state["analysis"], identity)

        def store():
            update_store(self.store, result)
            path = self.out / f"{self.name}_store.json"
            path.write_text(self.store.to_json())
            self.report.outputs["store"] = str(path)
            self.order_data.append({"order": l, "C": state["C"].expr.to_canonical_string(),
                                    "X": state["X"].expr.to_canonical_string()})
            od = self.out / f"{self.name}_order_data.json"
            od.write_text(json.dumps({"description": ORDER_DATA, "orders": self.order_data}, indent=1))
            self.report.outputs["order_data"] = str(od)
            entries = [e for e in self.store if e.loop_order == l]
            return {"entries": [e.name for e in entries],
                    "classes": {e.name: sum(len(f.leg_classes()) for f in e.laurent.values()) for e in entries}}

        self.step("store", store)
        res["loop_statistics"] = asdict(data.stats)
        return result, ctx, data

    def model_file(self):
        def write():
            written = counterterm_model(self.model, self.store,
                                        WriterOptions(schema=self.opt.schema, substitutions=self._gauge()))
            path = written.write(self.out / f"{self.name}_renormalised.json")
            loaded = written.load()
            self.written = written
            self.report.outputs["renormalised_model"] = str(path)
            ops: dict[str, int] = {}
            for v in written.data["vertex_rules"]:
                ops[v.get("operator", "")] = ops.get(v.get("operator", ""), 0) + 1
            self.check("renormalised model loads in HepKIT", len(loaded.vertex_rules) == len(written.data[
                "vertex_rules"]), f"{len(written.data['vertex_rules'])} vertex rules: {ops}")
            return {"vertex_rules": len(written.data["vertex_rules"]), "by_entry": ops,
                    "particles": len(written.data["particles"]), "couplings": len(written.data["couplings"])}

        self.step("model", write)

    @staticmethod
    def _beta_gamma(ctx: OrderContext, singular, l: int) -> tuple[dict[str, str], dict[str, str]]:
        beta, gamma = {}, {}
        for name, poly in singular.constants.items():
            e = ctx.space.expression(poly)
            if name.startswith("d["):
                beta[name[2:-1]] = str((Expression.num(2 * l) * e).expand())
            elif name.startswith("Z["):
                gamma[name[2:-1]] = str((Expression.num(-l) * e).expand())
        return beta, gamma


def renormalise(spec: Model | str | os.PathLike, options: RunOptions | None = None) -> RunResult:
    """Run the order loop up to ``options.order`` from a model specification and write the outputs (module
    docstring). Raises :class:`OrderNotSupported` for an order the loop cannot reach (before any computation);
    a failing step (no anomaly-free solution, data outside the span of a basis) is recorded in the report, which is
    written, and re-raised."""
    options = options or RunOptions()
    missing = missing_for_order(options.order)
    if missing:
        raise OrderNotSupported(f"order {options.order} is not supported; missing: " + "; ".join(missing))
    model = resolve_model(spec)
    run = _Run(model, options)
    orders, contexts, datas = [], [], []
    try:
        run.tree()
        for l in range(1, options.order + 1):
            res, ctx, data = run.run_order(l)
            orders.append(res)
            contexts.append(ctx)
            datas.append(data)
        run.model_file()
        run.report.status = "done"
    except Exception as exc:
        run.report.status = "failed"
        run.report.error = f"{type(exc).__name__}: {exc}"
        run.report.write(run.out)
        raise
    run.report.write(run.out)
    return RunResult(model, run.brst, run.store, orders, contexts, run.written, run.report, datas)

