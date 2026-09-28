"""M0 diagnostics ported from the plan's Python probes, outside primary clocks."""

from collections import defaultdict
from math import prod
from time import process_time_ns

from symbolica.community.hep import Symbols
from symbolica import E, Replacement, S, T
from symbolica.community.spenso import (
    AUTO,
    Representation,
    TensorExpression,
    TensorName,
    TensorRule,
)
from tensor_benchmark_cases import atom, describe


def observe(call):
    start = process_time_ns()
    try:
        value = call()
        elapsed = process_time_ns() - start
        return value, {"cpu_ns": elapsed, "status": "ok"}
    except Exception as error:  # noqa: BLE001 -- retain diagnostic refusal records
        return None, {
            "cpu_ns": process_time_ns() - start,
            "status": "error",
            "exception": type(error).__name__,
            "message": str(error),
        }


def validation_cost(case):
    """Same fixed inputs for reinference, no-match replacement and rerun probes."""
    absent = TensorName("r3::absent").to_expression()
    old = TensorName("gluon_ladder_typed::vx").to_expression()
    a = S("r3::args___")
    pattern = case.fixture["ladder_native_patterns"][0].replace(old(a), absent(a))
    no_match = TensorRule(
        pattern, case.fixture["ladder_native_rhs"], rhs_cache_size=1000
    )
    stages = {}
    result = case.source
    for step, vertex in enumerate(case.order, 1):
        result = result.replace(case.rules[vertex - 1]).contract()
        stages[step] = result
    stages[8] = result.expand()
    rows = []
    for stage in (5, 7, 8):
        source = stages[stage]
        # Resolving the diagnostic input is outside the measured admission call.
        raw = atom(source)
        for name, call in (
            ("contract", lambda source=source: source.contract()),
            ("expand", lambda source=source: source.expand()),
            (
                "no_match_replace",
                lambda source=source: source.replace(no_match),
            ),
            ("reinfer", lambda raw=raw: TensorExpression(raw)),
        ):
            value, row = observe(call)
            row.update(stage=stage, operation=name, input=describe(source))
            if value is not None:
                row.update(
                    output=describe(value), exact_identity=atom(value) == atom(source)
                )
            rows.append(row)
    return rows


def partial_parse(case, limit=200):
    result, stages = case.source, {}
    for step, vertex in enumerate(case.order, 1):
        result = result.replace(case.rules[vertex - 1]).contract()
        if step in (5, 7, 8):
            stages[step] = result.expand() if step == 8 else result
    rows = []
    for step, source in stages.items():
        terms = [TensorExpression(t) for t in list(atom(source).terms())[:limit]]
        for method in (
            "contract",
            "to_network",
            "canonize",
            "to_dots",
        ):
            values, row = observe(
                lambda method=method, terms=terms: [getattr(t, method)() for t in terms]
            )
            row.update(
                stage=step,
                operation=method,
                terms=len(terms),
                scope="first normalized terms; construction excluded",
            )
            if values is not None:
                row["outputs"] = len(values)
            rows.append(row)
        for count in (50, 100, 200):
            selected = terms[:count]
            if not selected:
                continue
            whole = TensorExpression(sum((atom(t) for t in selected), E("0")))
            singles, individual = observe(
                lambda selected=selected: [t.contract() for t in selected]
            )
            combined, batch = observe(whole.contract)
            row = {
                "stage": step,
                "operation": "contract_identical_termset",
                "terms": len(selected),
                "requested_terms": count,
                "individual": individual,
                "whole_sum": batch,
                "scope": "Same first N normalized terms; term/sum construction and exact checks outside clocks",
            }
            if singles is not None and combined is not None:
                expected = sum((atom(v) for v in singles), E("0"))
                row["exact"] = (atom(combined) - expected).expand() == E("0")
                assert row["exact"]
                row["whole_over_individual"] = batch["cpu_ns"] / individual["cpu_ns"]
            rows.append(row)
    return rows


def keep_rest():
    rep, spin = Representation.mink(4), Representation.bis(4)
    mu, nu, rho = [rep(S(f"r3rest::{n}")).to_expression() for n in ("mu", "nu", "rho")]
    i, j = [spin(S(f"r3rest::{n}")).to_expression() for n in ("i", "j")]
    x, y, z = S("r3rest::x", "r3rest::y", "r3rest::z")
    f = S("r3rest::f", is_scalar=True)
    v, p, q = [
        TensorName.vector(f"r3rest::{n}").to_expression() for n in ("V", "p", "q")
    ]
    gamma = TensorName.gamma().to_expression()
    a = rep(S("r3rest::a_")).to_expression()
    rule = TensorRule(v(a), p(a) + 2 * q(a))
    generator = TensorExpression.t(8, 3)
    color = (
        generator("r3rest::ca", "r3rest::ci", "r3rest::cj")
        * generator("r3rest::ca", "r3rest::cj", "r3rest::ck")
    ).to_expression()
    cases = {
        "plain": v(mu) * p(mu),
        "factored_scalar": (1 + x) * v(mu) * p(mu),
        "scalar_function": f(z) * v(mu) * p(mu),
        "nested_sum": (1 + x) * v(mu) * p(mu) * (y + f(z) * v(nu) * q(nu)),
        "foreign_color": color * v(mu) * p(mu),
        "foreign_gamma": gamma(i, j, rho) * v(mu) * p(mu),
        "connected_gamma": gamma(i, j, mu) * v(mu),
        "vector_power": v(mu) * p(mu) * q(nu) ** 2,
        "float": E("1.5") * v(mu) * p(mu),
        "symbolic_power": x ** S("r3rest::n") * v(mu) * p(mu),
    }
    rows = []
    for name, expression in cases.items():
        source = TensorExpression(expression)
        oracle, oracle_record = observe(
            lambda source=source: source.replace(rule).expand().contract()
        )
        output, row = observe(lambda source=source: source.replace(rule).contract())
        row.update(
            case=name,
            route="replace_then_contract",
            input=expression.format_plain(),
            oracle=oracle_record,
            required_milestone="M3"
            if name in ("vector_power", "foreign_color")
            else "M0",
        )
        if output is not None:
            row.update(
                exact_expanded=(atom(output) - atom(oracle)).expand() == E("0")
                if oracle is not None
                else None,
                output=atom(output).format_plain(),
                rank=output.root.structure.rank,
                preserves_outer_scalar_factor=bool(atom(output).contains(1 + x))
                if name in ("factored_scalar", "nested_sum")
                else None,
                rerun=atom(output.contract()) == atom(output),
            )
        rows.append(row)
    gam = TensorExpression.gamma(4)
    metric = TensorExpression.g(rep)
    ends = {
        "explicit": (spin(S("r3rest::i")), spin(S("r3rest::j"))),
        "AUTO": (AUTO, AUTO),
    }
    for label, indices in ends.items():
        for branches in (1, 2):
            gamma_value = gam(*indices, rep(S("r3rest::nu")))
            if branches == 2:
                gamma_value = gamma_value + gam(*indices, rep(S("r3rest::nu")))
            source = metric(rep(S("r3rest::mu")), rep(S("r3rest::nu"))) * gamma_value
            for method in ("contract", "simplify_gamma"):
                output, row = observe(getattr(source, method))
                row.update(
                    case="metric_into_gamma",
                    spinor_slots=label,
                    branches=branches,
                    method=method,
                    required_milestone="M3"
                    if label == "AUTO" and branches == 2
                    else "M0",
                )
                if output is not None:
                    row["output"] = atom(output).format_plain()
                    expected_gamma = (
                        gamma_value.simplify_gamma()
                        if method == "simplify_gamma"
                        else gamma_value
                    )
                    expected = atom(expected_gamma).replace(
                        rep(S("r3rest::nu")).to_expression(),
                        rep(S("r3rest::mu")).to_expression(),
                    )
                    row["exact_expanded"] = (atom(output) - expected).expand() == E("0")
                    row["rerun"] = atom(getattr(output, method)()) == atom(output)
                rows.append(row)
    return rows


class FactorProbe:
    """The original fixed-fixture DP/DAG prototype, not a tensor engine.

    Its admitted grammar is intentionally limited to the historical 4D ladder.
    One state transition implementation serves Atom and polynomial emission.
    """

    def __init__(self, case):
        self.case, self.ns = case, case.fixture
        self.slots = list(self.ns["ladder_indices"].values())
        self.rep = Representation.mink(4).to_expression()
        self.dot_fix = (
            S("f_")(S("h_")(self.rep)),
            Symbols.metric(S("f_")(self.rep), S("h_")(self.rep)),
        )
        self.internal, self.replaces = 0, 0

    @staticmethod
    def typename(value):
        return str(value.get_type()).split(".")[-1]

    def kind(self, value):
        if self.typename(value) == "Fn":
            return (
                "slot" if value.get_name() == Symbols.lorentz.get_name() else "vector"
            )
        return "other"

    @staticmethod
    def key(state):
        return tuple(sorted(state.items()))

    @staticmethod
    def holder(state, slot):
        return next((i for i, f in state.items() if f.contains(slot)), None)

    def substitute(self, state, index, slot, target, internal=False):
        self.replaces += 1
        value = state[index].replace(slot, target)
        if self.kind(target) == "vector":
            value = value.replace(*self.dot_fix)
        if internal:
            self.internal += 1
            value = atom(TensorExpression(value).contract())
        state[index] = value

    def apply(self, term, state):
        coefficient = E("1")
        for child in list(term) if self.typename(term) == "Mul" else [term]:
            kind = self.typename(child)
            if kind == "Num":
                coefficient *= child
                continue
            if kind == "Pow":
                base = next(iter(child))
                if (
                    self.typename(base) == "Fn"
                    and base.get_name() == Symbols.metric.get_name()
                    and all(self.kind(a) != "slot" for a in base)
                ):
                    coefficient *= child
                    continue
                raise ValueError(f"Prototype unsupported power: {child}")
            if kind != "Fn":
                raise ValueError(f"Prototype unsupported factor: {child}")
            if child.get_name() != Symbols.metric.get_name():
                args = list(child)
                if len(args) != 1 or self.kind(args[0]) != "slot":
                    raise ValueError(f"Prototype unsupported leaf: {child}")
                owner = self.holder(state, args[0])
                if owner is None:
                    coefficient *= child
                else:
                    self.substitute(
                        state, owner, args[0], child.replace(args[0], self.rep)
                    )
                continue
            left, right = list(child)
            lk, rk = self.kind(left), self.kind(right)
            if lk == rk == "slot":
                if left == right:
                    coefficient *= 4
                    continue
                lo, ro = self.holder(state, left), self.holder(state, right)
                if lo is None and ro is None:
                    coefficient *= child
                elif lo is None:
                    self.substitute(state, ro, right, left)
                else:
                    self.substitute(state, lo, left, right, lo == ro)
            elif "slot" in (lk, rk):
                slot, vector = (left, right) if lk == "slot" else (right, left)
                owner = self.holder(state, slot)
                if owner is None:
                    coefficient *= child
                else:
                    self.substitute(state, owner, slot, vector)
            else:
                coefficient *= child
        return coefficient

    def run(self, output="polynomial"):
        ns = self.ns
        rhs = ns["ladder_vertex_rule"].hold(T().expand())
        start = process_time_ns()
        factors = {
            i: vertex.replace(pattern, rhs, rhs_cache_size=1000)
            for i, (vertex, pattern) in enumerate(
                zip(ns["ladder_vertices"], ns["ladder_patterns"], strict=True), 1
            )
        }
        rules_ns = process_time_ns() - start
        start = process_time_ns()
        all_rules = prod(ns["ladder_vertices"]).replace_multiple(
            [Replacement(p, rhs) for p in ns["ladder_patterns"]]
        )
        all_rules_ns = process_time_ns() - start
        start = process_time_ns()
        incidence = [
            [i for i, value in factors.items() if value.contains(s)] for s in self.slots
        ]
        parse_ns = process_time_ns() - start
        root = self.key(factors)
        states, levels, counts = {root: None}, [], []
        edges, coefficient_bytes, expanded_transitions = 0, 0, 0
        start = process_time_ns()
        for vertex in self.case.order:
            current = defaultdict(list)
            for key in states:
                state = dict(key)
                factor = state.pop(vertex)
                indexed = any(factor.contains(s) for s in self.slots)
                terms = factor.terms() if indexed else [factor]
                for term in terms:
                    working = dict(state)
                    coefficient = self.apply(term, working) if indexed else term
                    current[self.key(working)].append((key, coefficient))
                    edges += 1
                    expanded_transitions += int(indexed)
                    coefficient_bytes += coefficient.get_byte_size()
            states = current
            levels.append(current)
            counts.append(
                {
                    "vertex": vertex,
                    "states": len(states),
                    "edges": sum(map(len, current.values())),
                }
            )
        contract_ns = process_time_ns() - start
        start = process_time_ns()
        weights = {root: E("1").to_polynomial() if output == "polynomial" else E("1")}
        for level in levels:
            weights = {
                key: sum(
                    (
                        weights[previous]
                        * (
                            coefficient.to_polynomial()
                            if output == "polynomial"
                            else coefficient
                        )
                        for previous, coefficient in incoming
                    )
                )
                for key, incoming in level.items()
            }
        (value,) = weights.values()
        forward_ns = process_time_ns() - start
        start = process_time_ns()
        value = value.to_expression() if output == "polynomial" else value.expand()
        emission_ns = process_time_ns() - start
        return value, {
            "per_vertex_rules_cpu_ns": rules_ns,
            "all_rules_cpu_ns": all_rules_ns,
            "factored_rules_bytes": all_rules.get_byte_size(),
            "partial_parse_cpu_ns": parse_ns,
            "incidence": incidence,
            "contraction_cpu_ns": contract_ns,
            "levels": counts,
            "edges": edges,
            "expanded_transitions": expanded_transitions,
            "emission_route": output + " forward pass over retained DAG",
            "coefficient_atom_bytes": coefficient_bytes,
            "byte_scope": "sum of edge coefficient Atom sizes; excludes graph/state/container allocation",
            "internal_contractions": self.internal,
            "replacement_calls": self.replaces,
            "forward_cpu_ns": forward_ns,
            "materialization_cpu_ns": emission_ns,
            **describe(value),
        }
