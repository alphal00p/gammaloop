"""Cases for fermion_ladder.py's consolidation suite; no benchmark scheduler."""

from __future__ import annotations

import hashlib
import re
from math import prod
from pathlib import Path
from time import process_time_ns

from symbolica import E, S, T
from symbolica.community.hep import Symbols
from symbolica.community.spenso import (
    AUTO,
    AliasedTensorExpression,
    CookSettings,
    GammaSimplifySettings,
    Representation,
    TensorExpression,
    TensorName,
    trace,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ORDERS = {
    "original": tuple(range(1, 9)),
    "early": (5, 4, 6, 3, 7, 2, 8, 1),
}


def atom(value):
    if isinstance(value, AliasedTensorExpression):
        value = value.to_expression()
    return value.to_expression() if isinstance(value, TensorExpression) else value


def describe(value):
    expression = atom(value)
    tensor = value.root if isinstance(value, AliasedTensorExpression) else value
    return {
        "terms": sum(1 for _ in expression.terms()),
        "atom_bytes": expression.get_byte_size(),
        "sha256": hashlib.sha256(expression.format_plain().encode()).hexdigest(),
        "rank": tensor.structure.rank if isinstance(tensor, TensorExpression) else None,
    }


def form_batch(declarations, seed, body):
    """Same independent-expression/body-timer boundary as the ladder FORM driver."""
    return f"""Off Statistics;
Format nospaces;
{declarations}
#do j=1,`REPEATS'
  Local F`j' = {seed};
#enddo
.sort
#reset timer
{body}
.sort
#if `DIAGNOSTIC' == 0
  #write "BATCH_CPU_MS=`timer_'"
#else
  #write <`OUTPUT'> "%E",F1
#endif
#$terms = termsin_(F1);
#write "TERMS=%$",$terms
#do j=2,`REPEATS'
  Local Check`j' = F`j'-F1;
#enddo
.sort
#do j=2,`REPEATS'
  #$different = termsin_(Check`j');
  #if `$different' != 0
    #terminate
  #endif
#enddo
.end
"""


class HistoricalLadder:
    """Historical eight-vertex fixture, external k10/k20 polarizations."""

    def __init__(self, order="original"):
        self.mode, self.loops, self.order = "trace4", 4, ORDERS[order]
        self.form_order = order
        lorentz = Representation.mink(4)
        ladder_order = self.order
        (
            ladder_expected_counts,
            ladder_input,
            ladder_indices,
            ladder_metric,
            ladder_momenta,
            ladder_patterns,
            ladder_rhs,
            ladder_vertices,
            ladder_vertex_rule,
        ) = self.geometry(E, S, TensorName, lorentz, prod)
        (
            ladder_native_input,
            ladder_native_patterns,
            ladder_native_rhs,
            ladder_native_rules,
        ) = self.native_fixture(
            S, TensorName, ladder_patterns, ladder_vertex_rule, ladder_vertices, prod
        )
        (
            ladder_outside_absorb,
            ladder_outside_d,
            ladder_outside_input,
            ladder_outside_k,
            ladder_outside_rhs,
        ) = self.scalar_fixture(
            S,
            ladder_indices,
            ladder_input,
            ladder_metric,
            ladder_momenta,
            ladder_vertex_rule,
        )
        (ladder_form_source,) = self.form_fixture(ladder_order)
        self.fixture = {
            "ladder_expected_counts": ladder_expected_counts,
            "ladder_input": ladder_input,
            "ladder_indices": ladder_indices,
            "ladder_metric": ladder_metric,
            "ladder_momenta": ladder_momenta,
            "ladder_patterns": ladder_patterns,
            "ladder_rhs": ladder_rhs,
            "ladder_vertices": ladder_vertices,
            "ladder_vertex_rule": ladder_vertex_rule,
            "ladder_native_input": ladder_native_input,
            "ladder_native_patterns": ladder_native_patterns,
            "ladder_native_rhs": ladder_native_rhs,
            "ladder_native_rules": ladder_native_rules,
            "ladder_outside_absorb": ladder_outside_absorb,
            "ladder_outside_d": ladder_outside_d,
            "ladder_outside_input": ladder_outside_input,
            "ladder_outside_k": ladder_outside_k,
            "ladder_outside_rhs": ladder_outside_rhs,
            "ladder_form_source": ladder_form_source,
        }
        self.source = TensorExpression(ladder_native_input)
        self.rules = ladder_native_rules
        self.strategy = {
            "order": self.order,
            "vertices": 8,
            "internal_edges": 11,
            "scope": "Typed replace, ambient contraction, final typed expansion",
        }

    @staticmethod
    def geometry(E, S, TensorName, lorentz, prod):
        from symbolica import T

        ladder_momenta = {
            i: TensorName.vector(f"gluon_ladder::k{i}").to_expression()(
                lorentz.to_expression()
            )
            for i in (0, 1, 2, 3, 4, 10, 20)
        }
        _k = ladder_momenta
        ladder_indices = {
            i: lorentz(S(f"gluon_ladder::mu{i}")).to_expression() for i in range(1, 12)
        }
        _mu = ladder_indices
        ladder_metric = TensorName.g().to_expression()
        _vx = S("gluon_ladder::vx")
        # All three momenta enter each vertex. The final three arguments are
        # Lorentz slots, or compact momenta for the two external polarizations.
        _routing = [
            (-_k[0], _k[0] - _k[1], _k[1], _k[10], _mu[1], _mu[8]),
            (-_k[1], _k[2], _k[1] - _k[2], _mu[1], _mu[2], _mu[9]),
            (-_k[2], _k[3], _k[2] - _k[3], _mu[2], _mu[3], _mu[10]),
            (-_k[3], _k[4], _k[3] - _k[4], _mu[3], _mu[4], _mu[11]),
            (-_k[4], _k[0], _k[4] - _k[0], _mu[4], _k[20], _mu[5]),
            (-_k[4] + _k[0], -_k[3] + _k[4], _k[3] - _k[0], _mu[5], _mu[11], _mu[6]),
            (-_k[3] + _k[0], -_k[2] + _k[3], _k[2] - _k[0], _mu[6], _mu[10], _mu[7]),
            (-_k[2] + _k[0], -_k[1] + _k[2], _k[1] - _k[0], _mu[7], _mu[9], _mu[8]),
        ]
        assert all(sum(vertex[:3], E("0")) == E("0") for vertex in _routing)
        ladder_vertices = [_vx(i, *vertex) for i, vertex in enumerate(_routing, 1)]
        ladder_input = prod(ladder_vertices)
        _p1, _p2, _p3, _i1, _i2, _i3 = S(
            *(f"gluon_ladder::{name}_" for name in ("p1", "p2", "p3", "i1", "i2", "i3"))
        )
        _g = ladder_metric
        ladder_vertex_rule = (
            -_g(_i1, _i3) * _g(_p1, _i2)
            + _g(_i1, _i2) * _g(_p1, _i3)
            + _g(_i2, _i3) * _g(_p2, _i1)
            - _g(_i1, _i2) * _g(_p2, _i3)
            - _g(_i2, _i3) * _g(_p3, _i1)
            + _g(_i1, _i3) * _g(_p3, _i2)
        )
        ladder_rhs = ladder_vertex_rule.hold(T().expand())
        ladder_patterns = [_vx(i, _p1, _p2, _p3, _i1, _i2, _i3) for i in range(1, 9)]
        ladder_expected_counts = (6, 34, 192, 1084, 6069, 10947, 23937, 9652)
        return (
            ladder_expected_counts,
            ladder_input,
            ladder_indices,
            ladder_metric,
            ladder_momenta,
            ladder_patterns,
            ladder_rhs,
            ladder_vertices,
            ladder_vertex_rule,
        )

    @staticmethod
    def native_fixture(
        S,
        TensorName,
        ladder_patterns,
        ladder_vertex_rule,
        ladder_vertices,
        prod,
    ):
        from symbolica.community.spenso import TensorRule as _TensorRule

        _vertex = TensorName("gluon_ladder_typed::vx").to_expression()
        _routing = S("gluon_ladder_typed::routing", is_scalar=True)

        def _with_ports(value):
            _number, _p1, _p2, _p3, _a, _b, _c = value
            # Momentum routing is opaque metadata. Only the last three arguments
            # are vertex ports; compact tagged vectors consume those ports.
            return _vertex(
                _number, _routing(_p1), _routing(_p2), _routing(_p3), _a, _b, _c
            )

        ladder_native_input = prod(_with_ports(vertex) for vertex in ladder_vertices)
        ladder_native_patterns = [_with_ports(pattern) for pattern in ladder_patterns]
        ladder_native_rhs = ladder_vertex_rule
        ladder_native_rules = [
            _TensorRule(pattern, ladder_native_rhs, rhs_cache_size=1000)
            for pattern in ladder_native_patterns
        ]
        return (
            ladder_native_input,
            ladder_native_patterns,
            ladder_native_rhs,
            ladder_native_rules,
        )

    @staticmethod
    def scalar_fixture(
        S,
        ladder_indices,
        ladder_input,
        ladder_metric,
        ladder_momenta,
        ladder_vertex_rule,
    ):
        from symbolica import T as _T

        # Translate the existing routing and rule; do not define a second graph.
        ladder_outside_d = S(
            "gluon_ladder_outside::d", is_symmetric=True, is_linear=True
        )
        ladder_outside_k = S("gluon_ladder_outside::k")
        _vertex = S("gluon_ladder::vx")
        _a, _b = S("gluon_ladder_outside::a_", "gluon_ladder_outside::b_")
        _index = S("gluon_ladder_outside::index_", tags=["ladder_outside_index"])
        _before, _after = S(
            "gluon_ladder_outside::before___", "gluon_ladder_outside::after___"
        )
        ladder_outside_input = ladder_input
        for _number, _momentum in ladder_momenta.items():
            ladder_outside_input = ladder_outside_input.replace(
                _momentum, ladder_outside_k(_number)
            )
        for _number, _slot in ladder_indices.items():
            _plain = S(
                f"gluon_ladder_outside::mu{_number}", tags=["ladder_outside_index"]
            )
            ladder_outside_input = ladder_outside_input.replace(_slot, _plain)
        _d = ladder_outside_d
        _contract = (
            _T()
            .repeat(
                _T().replace(
                    _d(_index, _a) * _d(_b, _index),
                    _d(_a, _b),
                    min_level=0,
                    max_level=0,
                )
            )
            .replace(_d(_index, _b) ** 2, _d(_b, _b), min_level=0, max_level=0)
            .replace(_d(_index, _index), 4, min_level=0, max_level=0)
        )
        ladder_outside_rhs = ladder_vertex_rule.replace(
            ladder_metric(_a, _b), _d(_a, _b)
        ).hold(_T().expand().chain(_contract))
        ladder_outside_absorb = (
            _d(_index, _a) * _vertex(_before, _index, _after),
            _vertex(_before, _a, _after),
        )
        return (
            ladder_outside_absorb,
            ladder_outside_d,
            ladder_outside_input,
            ladder_outside_k,
            ladder_outside_rhs,
        )

    @staticmethod
    def form_fixture(ladder_order):
        # Preserve the supplied rules and routing; select the same order as Python.
        ladder_form_source = r"""#-
format nospaces;
Dimension 4;
Auto Index mu;
Auto Vector k;
CF vx, f;
S x;

L F = vx(1,-k0, k0-k1, k1, k10, mu1, mu8)*
    vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)*
    vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)*
    vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)*
    vx(5,-k4, k0, k4-k0, mu4, k20, mu5)*
    vx(6,-k4+k0, -k3+k4, k3-k0, mu5, mu11, mu6)*
    vx(7,-k3+k0, -k2+k3, k2-k0, mu6, mu10, mu7)*
    vx(8,-k2+k0, -k1+k2, k1-k0, mu7, mu9, mu8);

#do i=1,8
    id vx(`i', k1?, k2?, k3?, mu1?, mu2?, mu3?) =
        (- d_(mu1, mu3) * d_(k1, mu2)
                    + d_(mu1, mu2) * d_(k1, mu3)
                    + d_(mu2, mu3) * d_(k2, mu1)
                    - d_(mu1, mu2) * d_(k2, mu3)
                    - d_(mu2, mu3) * d_(k3, mu1)
                    + d_(mu1, mu3) * d_(k3, mu2)
                    );

*    Print +s;
    .sort:gluon-`i';
#enddo

*Print +s;
.end
"""
        if ladder_order != tuple(range(1, 9)):
            _prefix, _loop = ladder_form_source.split("#do i=1,8\n")
            _body, _suffix = _loop.split("#enddo", 1)
            ladder_form_source = (
                _prefix
                + "".join(_body.replace("`i'", str(_i)) for _i in ladder_order)
                + _suffix
            )
        return (ladder_form_source,)

    def typed(self, observer=None):
        result = self.source
        for vertex in self.order:
            start = process_time_ns() if observer is not None else 0
            result = result.replace(self.rules[vertex - 1])
            replaced = process_time_ns() if observer is not None else 0
            result = result.contract()
            contracted = process_time_ns() if observer is not None else 0
            if observer is not None:
                observer.append(
                    {
                        "vertex": vertex,
                        "replacement_cpu_ns": replaced - start,
                        "contraction_cpu_ns": contracted - replaced,
                        **describe(result),
                    }
                )
        start = process_time_ns() if observer is not None else 0
        result = result.expand()
        elapsed = process_time_ns() - start if observer is not None else 0
        if observer is not None:
            observer.append(
                {"phase": "materialization", "cpu_ns": elapsed, **describe(result)}
            )
        return result

    def reduce(self):
        return self.typed().to_expression()

    def scalar(self):
        ns = self.fixture
        result = ns["ladder_outside_input"]
        for vertex in self.order:
            result = result.replace(
                ns["ladder_patterns"][vertex - 1],
                ns["ladder_outside_rhs"],
                rhs_cache_size=1000,
            )
            result = result.expand().replace(
                *ns["ladder_outside_absorb"], min_level=0, max_level=0, repeat=True
            )
        return result

    def convert_scalar(self, result):
        ns = self.fixture
        for i, p in ns["ladder_momenta"].items():
            for j, q in ns["ladder_momenta"].items():
                if i <= j:
                    result = result.replace(
                        ns["ladder_outside_d"](
                            ns["ladder_outside_k"](i), ns["ladder_outside_k"](j)
                        ),
                        ns["ladder_metric"](p, q),
                    )
        return result

    def phases(self, expected):
        stages = []
        assert atom(self.typed(stages)) == expected
        return {"stages": stages, "internal_edges": 11}

    def form_source(self):
        source = self.fixture["ladder_form_source"]
        declarations, rest = source.split("L F = ", 1)
        seed, body = rest.split(";", 1)
        return form_batch(
            declarations.replace("#-", ""), seed, body.replace(".end", "")
        )

    def import_form(self, text):
        text = re.sub(r"\s+", "", text)
        ns = self.fixture
        text = re.sub(
            r"k(\d+)\.k(\d+)",
            lambda m: f"r3form::s{m[1]}x{m[2]}",
            text.strip().rstrip(";"),
        )
        result = E(text)
        for i, p in ns["ladder_momenta"].items():
            for j, q in ns["ladder_momenta"].items():
                result = result.replace(
                    S(f"r3form::s{i}x{j}"), ns["ladder_metric"](p, q)
                )
        return result.expand()


class FreeTrace:
    def __init__(self, mode, length):
        self.mode, self.loops, self.length = mode, 3, length
        self.dimension_symbol = S("r3trace::D")
        self.dimension = 4 if mode == "trace4" else self.dimension_symbol
        self.lorentz, self.spin = (
            Representation.mink(self.dimension),
            Representation.bis(4),
        )
        self.indices = [self.lorentz(S(f"r3trace::i{i}")) for i in range(length)]
        gamma = TensorExpression.gamma(self.dimension)
        self.source = trace(self.spin, *(gamma(AUTO, AUTO, i) for i in self.indices))
        self.settings = GammaSimplifySettings()
        self.strategy = {
            "gammas": length,
            "trace_materialization": "explicit alias expansion inside the operation",
            "dimension": str(self.dimension),
        }
        self.form_order = "trace-first"

    def reduce(self):
        return self.source.simplify_gamma(self.settings).expand().to_expression()

    def phases(self, expected):
        return {
            "final_terms": sum(1 for _ in expected.terms()),
            "free_ports": self.length,
        }

    def form_source(self):
        indices = ",".join(f"i{i}" for i in range(self.length))
        dimension = "4" if self.mode == "trace4" else "D"
        return form_batch(
            f"Symbols D;\nDimension {dimension};\nIndices {indices};\nUnitTrace 4;",
            f"g_(1,{indices})",
            f"{self.mode},1;",
        )

    def import_form(self, text):
        text = re.sub(r"\s+", "", text)
        g = TensorName.g().to_expression()
        metrics = {}

        def convert(match):
            pair = tuple(map(int, match.groups()))
            if pair not in metrics:
                metrics[pair] = g(
                    *(self.indices[i].to_expression() for i in pair)
                ).format_plain()
            return metrics[pair]

        text = re.sub(r"\bD\b", "r3trace::D", text.strip().rstrip(";"))
        return E(re.sub(r"d_\(i(\d+),i(\d+)\)", convert, text)).expand()

    def check_form_components(self, native, reference):
        """Finite exact HEP checks; not a symbolic proof of a Schouten identity."""
        import random
        from itertools import permutations

        from symbolica import Expression, Replacement
        from symbolica.community.spenso import Tensor, TensorLibrary, TensorNetwork

        rep = Representation.mink(4)
        vectors = [
            TensorName.vector(f"r3components::v{i}").to_expression()(
                rep.to_expression()
            )
            for i in range(self.length)
        ]
        expressions = {
            "source": self.source.to_expression(),
            "Idenso": native,
            "reference": reference,
        }
        indices = [i.to_expression() for i in self.indices]
        if self.mode == "tracen":
            expressions = {
                key: value.replace(self.dimension_symbol, E("4"))
                for key, value in expressions.items()
            }
            indices = [
                index.replace(self.dimension_symbol, E("4")) for index in indices
            ]
        rules = [
            Replacement(index, vector)
            for index, vector in zip(indices, vectors, strict=True)
        ]
        rules += [
            Replacement(S("r3trace::x"), E("1")),
            Replacement(S("r3trace::y"), E("0")),
        ]
        expressions = {
            key: value.replace_multiple(rules) for key, value in expressions.items()
        }
        rows = []
        for seed in (47, 98, 173):
            rng = random.Random(seed)
            assignments = [
                [rng.randint(-3, 3) for _ in range(4)] for _ in range(self.length)
            ]
            signed_permutations = [
                (
                    permutation,
                    (-1)
                    ** sum(
                        permutation[i] > permutation[j]
                        for i in range(4)
                        for j in range(i + 1, 4)
                    ),
                )
                for permutation in permutations(range(4))
            ]
            determinant = sum(
                sign * prod(assignments[i][permutation[i]] for i in range(4))
                for permutation, sign in signed_permutations
            )
            assert determinant != 0, "HEP vector sample must span four dimensions"
            library = TensorLibrary.hep_lib_atom()
            # Same explicit epsilon components as gamma_simplification.py's HEP check.
            epsilon = Tensor.sparse(
                TensorName(Symbols.levi_civita.get_name())(*[rep] * 4),
                Expression,
            )
            for permutation, sign in signed_permutations:
                epsilon[permutation] = E("-1𝑖" if sign == 1 else "1𝑖")
            library.register(epsilon)
            for vector, components in zip(vectors, assignments, strict=True):
                tensor = Tensor.sparse(TensorExpression(vector), Expression)
                for axis, value in enumerate(components):
                    tensor[axis] = E(str(value))
                library.register(tensor)

            def evaluate_leaf(expression, library=library):
                network = TensorNetwork(expression, library=library)
                network.execute(library=library)
                return network.result_scalar()

            wildcards = S(
                "r3components::a_",
                "r3components::b_",
                "r3components::c_",
                "r3components::d_",
            )
            metric_pattern = TensorName.g().to_expression()(*wildcards[:2])
            epsilon_pattern = Symbols.levi_civita(*wildcards)
            results = {}
            for label, expression in expressions.items():
                if label == "source":
                    value = evaluate_leaf(expression)
                else:
                    # Cache exact HEP evaluations of repeated tensor leaves; the
                    # remaining arithmetic is scalar and needs no network graph.
                    value = expression
                    for pattern in (metric_pattern, epsilon_pattern):
                        value = value.replace(
                            pattern,
                            pattern.hold(T().map(evaluate_leaf)),
                            rhs_cache_size=2048,
                        )
                assert str(value.get_type()).split(".")[-1] == "Num", (
                    "HEP result must be fully numeric"
                )
                results[label] = value.format_plain()
            assert len(set(results.values())) == 1, (seed, results)
            assert results["source"] != "0", "Trace component sample must not vanish"
            rows.append(
                {
                    "seed": seed,
                    "assignments": assignments,
                    "spanning_minor": determinant,
                    "results": results,
                    "exact_at_point": True,
                }
            )
        return {
            "scope": "Three finite exact HEP points: original source network and cached metric/epsilon leaf networks for outputs; no general polynomial identity claim",
            "dimension": 4,
            "rows": rows,
        }


class AxialTrace(FreeTrace):
    def __init__(self):
        super().__init__("trace4", 12)
        gamma = TensorExpression.gamma(4)
        word = trace(
            self.spin,
            TensorExpression.gamma5(4)(AUTO, AUTO),
            *(gamma(AUTO, AUTO, index) for index in self.indices),
        )
        self.spectator = sum(S("r3trace::x", "r3trace::y")) ** 8
        self.source = self.spectator * word
        self.strategy.update(
            gamma5=True,
            scalar_spectator="(x+y)^8 remains factored",
            form_scope="FORM scalar symbol C represents factored (x+y)^8; scalar-factor preparation excluded in both",
        )

    def reduce(self):
        reduced = self.source.simplify_gamma(self.settings)
        # M1 and FORM materialize the trace while retaining the scalar factor.
        # Applying expand to the entire aliased tensor would also distribute
        # (x+y)^8, which is outside this case's declared materialization scope.
        trace = AliasedTensorExpression(reduced.root / self.spectator, reduced.aliases)
        return self.spectator * trace.expand().to_expression()

    def form_source(self):
        indices = ",".join(f"i{i}" for i in range(self.length))
        return form_batch(
            f"Dimension 4; Symbols C; Indices {indices}; UnitTrace 4;",
            f"C*g_(1,5_,{indices})",
            "trace4,1;",
        )

    def import_form(self, text):
        text = re.sub(r"\s+", "", text)
        epsilon = Symbols.levi_civita

        def convert(match):
            return epsilon(
                *(self.indices[int(i)].to_expression() for i in match.groups())
            ).format_plain()

        text = re.sub(r"e_\(i(\d+),i(\d+),i(\d+),i(\d+)\)", convert, text)
        # The repository's epsilon convention already contains FORM's imaginary unit.
        value = super().import_form(text)
        return self.spectator * value.replace(S("C"), E("1"))


class ProductionNumerator:
    """Unmodified captured aa→aa GL16 post-metric numerator; no invented graph."""

    def __init__(self):
        self.path = (
            REPO
            / "crates/idenso/tests/fixtures/aa_aa_2l_gl16_integrated_uv_start_after_simplify_metrics.sym"
        )
        # Q is registered by the combined host, including its display metadata.
        TensorName.gamma()
        TensorName.t()
        Representation.mink(Symbols.dimension)
        Representation.bis(4)
        Representation.cof(3)
        Representation.coad(8)
        self.source = TensorExpression(
            E(self.path.read_text()), cook_indices=CookSettings.indices()
        )
        self.strategy = {
            "fixture": str(self.path.relative_to(REPO)),
            "sha256": hashlib.sha256(self.path.read_bytes()).hexdigest(),
            "scope": "Captured post-metric numerator: gamma, color, then metric simplification",
            "input_adapter": "Existing CookSettings.indices() for nested hedge/edge index payloads; setup excluded",
        }

    def reduce(self):
        return (
            self.source.simplify_gamma()
            .simplify_color()
            .contract()
            .to_expression()
            .to_expression()
        )

    def check_interpreter_components(self, native, reference, output):
        """Exact component checks for this capture, not an arbitrary-momentum proof."""
        import json
        import random
        from itertools import product

        from symbolica import Expression, Replacement
        from symbolica.community.spenso import Tensor, TensorLibrary, TensorNetwork

        dimension = S("gammalooprs::dim")
        rep = Representation.mink(4)
        momenta = [
            S("gammalooprs::Q")(E(str(index)), rep.to_expression())
            for index in range(5, 11)
        ]
        # Fix the component representation while retaining scalar D coefficients.
        rules = [
            Replacement(
                E("spenso::mink(gammalooprs::dim,a_)"), E("spenso::mink(4,a_)")
            ),
            Replacement(E("spenso::mink(gammalooprs::dim)"), E("spenso::mink(4)")),
            Replacement(S("UFO::MT"), E("2")),
            Replacement(S("UFO::GC_2"), E("1")),
            Replacement(S("UFO::GC_11"), E("1")),
            Replacement(E("spenso::idx(2,spenso::cof(3))"), E("1/2")),
        ]
        expressions = {
            label: value.replace_multiple(rules)
            for label, value in {
                "source": atom(self.source),
                "reference": reference,
                "native": native,
            }.items()
        }
        slots = [
            slot.to_expression().replace(dimension, E("4")).format_plain()
            for slot in self.source.structure.slots
        ]
        assert len(slots) == 4 and len(set(slots)) == 4
        report = {
            "passed": False,
            "oracle": "exact HEP components and dimension-only scalar coefficients",
            "scope": (
                "All 256 external 4D components at three fixed six-momentum assignments; "
                "not a proof for arbitrary momenta or a generic-D Clifford algebra"
            ),
            "dimension_scope": (
                "Reference/native scalar D coefficients agree after finite 4D tensor "
                "contraction; degree at most one for this capture"
            ),
            "normalizations": {"MT": 2, "GC_2": 1, "GC_11": 1, "idx(2,cof(3))": "1/2"},
            "external_slots": slots,
            "samples": [],
        }
        output = Path(output)
        output.write_text(json.dumps(report, indent=2) + "\n")
        for seed in (47, 98, 173):
            rng = random.Random(seed)
            assignments = [[rng.randint(-3, 3) for _ in range(4)] for _ in momenta]
            library = TensorLibrary.hep_lib_atom()
            # Generic MixedTensor metrics use f64 constants. Match the existing
            # metric_contraction_validation oracle with explicit Atom entries.
            metric = Tensor.sparse(TensorExpression.g(rep, rep), Expression)
            for axis in range(4):
                metric[axis, axis] = E("1" if axis == 0 else "-1")
            library.register(metric)
            for vector, components in zip(momenta, assignments, strict=True):
                tensor = Tensor.sparse(TensorExpression(vector), Expression)
                for axis, value in enumerate(components):
                    tensor[axis] = E(str(value))
                library.register(tensor)
            row = {"seed": seed, "momenta": assignments, "results": {}}
            report["samples"].append(row)
            for label, expression in expressions.items():
                network = TensorNetwork(expression, library=library)
                network.execute(library=library)
                tensor = network.result_tensor(library=library)
                actual_slots = [
                    slot.to_expression().format_plain()
                    for slot in tensor.structure.slots
                ]
                assert len(actual_slots) == 4 and set(actual_slots) == set(slots)
                axes = [actual_slots.index(slot) for slot in slots]
                tensor = tensor.permute_axes(axes)
                values = [tensor[list(index)] for index in product(range(4), repeat=4)]
                assert all(isinstance(value, Expression) for value in values)
                assert set().union(
                    *(set(value.get_all_symbols(True)) for value in values)
                ) <= {dimension}
                # Tensor execution has finished. Only the remaining scalar D is
                # collected, never a graph numerator or tensor payload.
                coefficients = [value.coefficient_list(dimension) for value in values]
                assert all(
                    monomial in (E("1"), dimension)
                    and re.fullmatch(
                        r"-?\d+(?:/\d+)?(?:𝑖)?", coefficient.format_plain()
                    )
                    for parts in coefficients
                    for monomial, coefficient in parts
                ), (seed, label, "non-rational or higher-degree dimension coefficient")
                components = [
                    value.replace(dimension, E("4")).format_plain() for value in values
                ]
                assert all(
                    re.fullmatch(r"-?\d+(?:/\d+)?(?:𝑖)?", value) for value in components
                )
                row["results"][label] = {
                    "components": components,
                    "dimension_coefficients": [
                        sorted(
                            (monomial.format_plain(), coefficient.format_plain())
                            for monomial, coefficient in parts
                        )
                        for parts in coefficients
                    ],
                    "original_slots": actual_slots,
                    "axes_to_common_order": axes,
                }
                output.write_text(json.dumps(report, indent=2) + "\n")
            results = row["results"]
            assert (
                results["source"]["components"]
                == results["reference"]["components"]
                == results["native"]["components"]
            ), seed
            assert any(value != "0" for value in results["source"]["components"]), seed
            assert (
                results["reference"]["dimension_coefficients"]
                == results["native"]["dimension_coefficients"]
            ), seed
        report["passed"] = True
        output.write_text(json.dumps(report, indent=2) + "\n")
        return {
            "passed": True,
            "oracle": report["oracle"],
            "scope": report["scope"],
            "dimension_scope": report["dimension_scope"],
            "path": str(output),
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        }

    def phases(self, expected):
        rows = []
        value = self.source
        for method in ("simplify_gamma", "simplify_color", "contract"):
            start = process_time_ns()
            value = getattr(value, method)()
            elapsed = process_time_ns() - start
            rows.append({"phase": method, "cpu_ns": elapsed, **describe(value)})
        assert atom(value) == expected
        return {"stages": rows}


def make_case(name):
    if name.startswith("historical-"):
        return HistoricalLadder(name.removeprefix("historical-"))
    if name == "axial12-trace4":
        return AxialTrace()
    if name.startswith("free"):
        length, mode = name[4:].split("-")
        return FreeTrace(mode, int(length))
    if name == "production-aa-aa":
        return ProductionNumerator()
    from fermion_ladder import FermionLadder
    from gluon_ladder import GluonLadder

    particle, loops, mode = name.split("-")
    return (FermionLadder if particle == "fermion" else GluonLadder)(mode, int(loops))


CASES = [
    *(f"fermion-{n}-{m}" for n in (3, 4) for m in ("trace4", "tracen")),
    *(f"gluon-4-{m}" for m in ("trace4", "tracen")),
    "historical-original",
    "historical-early",
    *(f"free{n}-{m}" for n in (4, 12, 14) for m in ("trace4", "tracen")),
    "axial12-trace4",
    "production-aa-aa",
]
