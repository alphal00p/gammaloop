"""Four-loop gluonic propagator numerator with verified physical routing.

The all-incoming vertex rule and its FORM spelling share RULE_TERMS. Integer
routing is checked at every vertex and at both ends of every internal edge.
Only scalar numerator algebra is included: no color, couplings or denominators.
"""

from __future__ import annotations

import re
from math import prod
from pathlib import Path
from time import process_time_ns

from symbolica import E, Replacement, S
from symbolica.community.spenso import (
    Representation,
    TensorExpression,
    TensorName,
    TensorRule,
)

# (sign, metric port pair, momentum leg, remaining component port), zero-based.
RULE_TERMS = (
    (-1, (0, 2), 0, 1),
    (1, (0, 1), 0, 2),
    (1, (1, 2), 1, 0),
    (-1, (0, 1), 1, 2),
    (-1, (1, 2), 2, 0),
    (1, (0, 2), 2, 1),
)
BASIS = ("k1", "k2", "k3", "k4", "q")
RING = (
    (1, 0, 0, 0, 0),
    (1, 1, 0, 0, 0),
    (1, 1, 1, 0, 0),
    (1, 1, 1, 1, 0),
    (1, 1, 1, 1, -1),
    (1, 1, 1, 0, -1),
    (1, 1, 0, 0, -1),
    (1, 0, 0, 0, -1),
)
EXCHANGES = ("mu", "a", "b", "c", "mu", "c", "b", "a")
ORDERS = {
    "original": tuple(range(1, 9)),
    "early": (5, 4, 6, 3, 7, 2, 8, 1),
    "mirror": (1, 2, 8, 3, 7, 4, 6, 5),
}


def graph_data():
    vertices = []
    for i, outgoing in enumerate(RING):
        incoming = RING[(i - 1) % 8]
        exchange = tuple(b - a for a, b in zip(incoming, outgoing, strict=True))
        momenta = (incoming, tuple(-x for x in outgoing), exchange)
        ports = (f"r{(i - 1) % 8 + 1}", f"r{i + 1}", EXCHANGES[i])
        assert all(sum(values) == 0 for values in zip(*momenta, strict=True))
        vertices.append((momenta, ports))
    edges = []
    for i in range(8):
        edges.append((i, 1, (i + 1) % 8, 0, "ring"))
    edges.extend((i, 2, j, 2, "rung") for i, j in ((1, 7), (2, 6), (3, 5)))
    for left, lp, right, rp, _kind in edges:
        assert vertices[left][1][lp] == vertices[right][1][rp]
        assert all(
            a + b == 0
            for a, b in zip(vertices[left][0][lp], vertices[right][0][rp], strict=True)
        )
    assert vertices[0][0][2] == (0, 0, 0, 0, 1)
    assert vertices[4][0][2] == (0, 0, 0, 0, -1)
    assert len(edges) - len(vertices) + 1 == 4
    # The shared external mu is metric projection, not a twelfth internal edge.
    return vertices, edges


class GluonLadder:
    """Physical ring and three rungs, projected with g(mu,nu); no fermion sign."""

    particle = "gluonic"
    form_order = "early"

    def __init__(self, mode="trace4", loops=4, order="mirror"):
        if mode not in ("trace4", "tracen") or loops != 4 or order not in ORDERS:
            raise ValueError(
                "Use trace4/tracen dimension labels, four loops, original/early/mirror"
            )
        self.mode, self.order, self.loops = mode, ORDERS[order], 4
        self.dimension_symbol = S("physical_gluon_ladder::D")
        self.dimension = 4 if mode == "trace4" else self.dimension_symbol
        self.routing = RING
        self.strategy = {
            "vertex_order": self.order,
            "form_vertex_order": ORDERS[self.form_order],
            "local": "compiled factorized tensor rule",
            "ambient": "shared graph contraction with retained aliases",
            "final": "typed expand then scalar notation",
            "trace_instruction": None,
        }
        self.lorentz = Representation.mink(self.dimension)
        self.metric = TensorName.g().to_expression()
        self.vertices, self.edges = graph_data()
        self.momenta = [
            TensorName.vector(f"physical_gluon_ladder::{name}").to_expression()(
                self.lorentz.to_expression()
            )
            for name in BASIS
        ]
        self.vertex = TensorName("physical_gluon_ladder::vx").to_expression()
        metadata = S("physical_gluon_ladder::routing", is_scalar=True)
        slots = {
            name: self.lorentz(S(f"physical_gluon_ladder::{name}")).to_expression()
            for name in sorted(
                {p for _, ports in self.vertices for p in ports} | {"nu"}
            )
        }

        def momentum(row):
            return sum(
                (
                    coefficient * p
                    for coefficient, p in zip(row, self.momenta, strict=True)
                    if coefficient
                ),
                E("0"),
            )

        def node(number, momenta, ports):
            return self.vertex(number, *(metadata(p) for p in momenta), *ports)

        source = []
        explicit = []
        templates = []
        for i, (momenta, ports) in enumerate(self.vertices):
            local = [momentum(row) for row in momenta]
            source.append(node(i + 1, local, [slots[p] for p in ports]))
            explicit_ports = [
                slots["nu"] if i == 4 and leg == 2 else slots[p]
                for leg, p in enumerate(ports)
            ]
            explicit.append(node(i + 1, local, explicit_ports))
            templates.append(node(i + 1, local, [self.lorentz.to_expression()] * 3))
        self.nodes, self.explicit_nodes, self.vertex_templates = (
            source,
            explicit,
            templates,
        )
        self.source = TensorExpression(prod(source))
        self.explicit = TensorExpression(
            self.metric(slots["mu"], slots["nu"]) * prod(explicit)
        )
        p = S(*(f"physical_gluon_ladder::p{i}_" for i in range(3)))
        ports = S(*(f"physical_gluon_ladder::i{i}_" for i in range(3)))
        self.patterns = [node(i, p, ports) for i in range(1, 9)]
        self.rule = sum(
            (
                sign * self.metric(ports[a], ports[b]) * self.metric(p[m], ports[c])
                for sign, (a, b), m, c in RULE_TERMS
            ),
            E("0"),
        )

        self.rules = [
            TensorRule(pattern, self.rule, rhs_cache_size=1000)
            for pattern in self.patterns
        ]
        self.scalars = {
            (i, j): S(f"physical_gluon_ladder::s{i + 1}{j + 1}")
            for i in range(5)
            for j in range(i, 5)
        }
        self.replacements = [
            Replacement(self.metric(self.momenta[i], self.momenta[j]), s)
            for (i, j), s in self.scalars.items()
        ]

    def reduce(self, source=None):
        result = self.source if source is None else source
        for vertex in self.order:
            result = result.replace(self.rules[vertex - 1]).contract()
        # Materialize once, then name scalar products for the exact FORM comparison.
        return result.expand().to_expression().replace_multiple(self.replacements)

    def import_form(self, polynomial):
        positions = {name: i + 1 for i, name in enumerate(BASIS)}

        def dot(match):
            i, j = sorted((positions[match[1]], positions[match[2]]))
            return f"physical_gluon_ladder::s{i}{j}"

        text = re.sub(
            r"\b(k[1-4]|q)\.(k[1-4]|q)\b", dot, polynomial.strip().rstrip(";")
        )
        text = re.sub(r"\bD\b", "physical_gluon_ladder::D", text)
        return E(text).expand(via_poly=True)

    def phases(self, expected):
        result = self.source
        rows = []
        for vertex in self.order:
            start = process_time_ns()
            result = result.replace(self.rules[vertex - 1])
            replaced = process_time_ns()
            result = result.contract()
            contracted = process_time_ns()
            rows.append(
                {
                    "vertex": vertex,
                    "replacement_cpu_ns": replaced - start,
                    "contraction_cpu_ns": contracted - replaced,
                }
            )
        start = process_time_ns()
        scalar = result.expand().to_expression().replace_multiple(self.replacements)
        final = process_time_ns() - start
        assert scalar == expected
        return {
            "stages": rows,
            "final_scalar_cpu_ns": final,
            "final_terms": len(list(scalar.terms())),
        }

    def check_components(self, scalar_results, networks):
        from gluon_ladder_validation import GluonRingComponents

        return GluonRingComponents(self).check(scalar_results, networks=networks)

    def form_source(self):
        """Generate the standalone FORM graph and rule from the Python data."""

        def momentum(row):
            terms = []
            for coefficient, name in zip(row, BASIS, strict=True):
                if coefficient:
                    terms.append(
                        ("+" if coefficient > 0 else "-")
                        + (str(abs(coefficient)) + "*" if abs(coefficient) != 1 else "")
                        + name
                    )
            return "".join(terms).lstrip("+") or "0"

        factors = [
            f"vx({i + 1},"
            + ",".join([*(momentum(row) for row in momenta), *ports])
            + ")"
            for i, (momenta, ports) in enumerate(self.vertices)
        ]
        rule = "".join(
            ("+" if sign > 0 else "-") + f"d_(i{a + 1},i{b + 1})*d_(p{m + 1},i{c + 1})"
            for sign, (a, b), m, c in RULE_TERMS
        ).lstrip("+")
        template = Path(__file__).with_name("gluon_propagator_ladder.frm").read_text()
        graph = "Local F`sample'=" + "*".join(factors) + ";"
        template = re.sub(r"(?m)^Local F`sample.*;$", lambda _: graph, template)
        return re.sub(
            r"(?m)^id vx.*;$",
            lambda _: "id vx(`vertex',p1?,p2?,p3?,i1?,i2?,i3?)=" + rule + ";",
            template,
        )
