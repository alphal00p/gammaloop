"""Cases for fermion_ladder.py's consolidation suite; no benchmark scheduler."""

from __future__ import annotations

import ast
import hashlib
import re
from math import prod
from pathlib import Path
from time import process_time_ns

from symbolica import E, S, T
from symbolica.community.spenso import (
    AUTO,
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
    return value.to_expression() if isinstance(value, TensorExpression) else value


def describe(value):
    expression = atom(value)
    return {
        "terms": sum(1 for _ in expression.terms()),
        "atom_bytes": expression.get_byte_size(),
        "sha256": hashlib.sha256(expression.format_plain().encode()).hexdigest(),
        "rank": value.structure.rank if isinstance(value, TensorExpression) else None,
    }


def notebook_cell(namespace, export):
    """Load a named pure fixture cell, without importing/running Marimo."""
    source = HERE / "gamma_simplification.py"
    nodes = [
        node
        for node in ast.parse(source.read_text()).body
        if isinstance(node, ast.FunctionDef)
        and isinstance(node.body[-1], ast.Return)
        and isinstance(node.body[-1].value, ast.Tuple)
        and any(
            isinstance(v, ast.Name) and v.id == export for v in node.body[-1].value.elts
        )
    ]
    if len(nodes) != 1:
        raise ValueError(f"Expected one notebook fixture cell exporting {export}")
    node = nodes[0]
    node.decorator_list = []
    exec(  # noqa: S102 -- guarded repository fixture, not user input
        compile(ast.Module(body=[node], type_ignores=[]), str(source), "exec"),
        namespace,
    )
    values = namespace["_"](**{arg.arg: namespace[arg.arg] for arg in node.args.args})
    namespace.update(zip((v.id for v in node.body[-1].value.elts), values, strict=True))


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
    """Exact notebook eight-vertex fixture, external k10/k20 polarizations."""

    def __init__(self, order="original"):
        self.mode, self.loops, self.order = "trace4", 4, ORDERS[order]
        self.form_order = order
        ns = {
            "E": E,
            "S": S,
            "T": T,
            "TensorName": TensorName,
            "TensorExpression": TensorExpression,
            "lorentz": Representation.mink(4),
            "prod": prod,
            "ladder_order": self.order,
        }
        for export in (
            "ladder_vertex_rule",
            "ladder_native_input",
            "ladder_outside_input",
            "ladder_form_source",
        ):
            notebook_cell(ns, export)
        self.fixture = ns
        self.source = TensorExpression(ns["ladder_native_input"])
        self.settings = ns["ladder_native_settings"]
        self.patterns, self.rhs = ns["ladder_native_patterns"], ns["ladder_native_rhs"]
        self.strategy = {
            "order": self.order,
            "vertices": 8,
            "internal_edges": 11,
            "scope": "Notebook typed replace, ambient contraction, final typed expansion",
        }

    def typed(self, observer=None):
        result = self.source
        for vertex in self.order:
            start = process_time_ns() if observer is not None else 0
            result = result.replace_tensor(
                self.patterns[vertex - 1], self.rhs, rhs_cache_size=1000
            )
            replaced = process_time_ns() if observer is not None else 0
            result = result.schoonschip(settings=self.settings)
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
        self.settings = GammaSimplifySettings(expand_traces=True)
        self.strategy = {
            "gammas": length,
            "expand_traces": True,
            "dimension": str(self.dimension),
        }
        self.form_order = "trace-first"

    def reduce(self):
        return self.source.simplify_gamma(self.settings).to_expression()

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
                TensorName("spenso::epsilon", is_antisymmetric=True)(*[rep] * 4),
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
            epsilon_pattern = S("spenso::epsilon")(*wildcards)
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

    def form_source(self):
        indices = ",".join(f"i{i}" for i in range(self.length))
        return form_batch(
            f"Dimension 4; Symbols C; Indices {indices}; UnitTrace 4;",
            f"C*g_(1,5_,{indices})",
            "trace4,1;",
        )

    def import_form(self, text):
        text = re.sub(r"\s+", "", text)
        epsilon = S("spenso::epsilon")

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
        Representation.mink(S("gammalooprs::dim"))
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
            .simplify_metrics()
            .to_expression()
        )

    def phases(self, expected):
        rows = []
        value = self.source
        for method in ("simplify_gamma", "simplify_color", "simplify_metrics"):
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
