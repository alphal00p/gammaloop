"""Reconstruct the LUDY DY/DIS integrands through the pre-IBP boundary.

Run in a current Symbolica community host; see investigation.typ for scope,
conventions, provenance, and the missing integration stages.
"""

import argparse
import importlib
import json
from itertools import combinations
from math import factorial, prod
from pathlib import Path

from symbolica import E, Expression, S
from symbolica.community.idenso import simplify_gamma, simplify_metrics, to_dots
from symbolica.community.spenso import Representation, TensorExpression, TensorName, dot

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--module", choices=("feynkit", "hep"), default="feynkit")
parser.add_argument("--process", choices=("dy", "dis", "both"), default="both")
parser.add_argument(
    "--check", action="store_true", help="compare separate reference oracles"
)
parser.add_argument(
    "--output", type=Path, help="write scalar numerators, cuts and integral indices"
)
args = parser.parse_args()
fk = importlib.import_module(f"symbolica.community.{args.module}")
directory = Path(__file__).resolve().parent
inputs = json.loads((directory / "inputs.json").read_text())
d, kk, kp1, kp2, p1sq, p2sq, p12, msq, x, Qsq = (
    S(f"ludy_port::{name}")
    for name in (
        "d",
        "kk",
        "kp1",
        "kp2",
        "p1_squared",
        "p2_squared",
        "p1_dot_p2",
        "mass_squared",
        "x",
        "Q_squared",
    )
)
lorentz = Representation.mink(dimension=d)
gamma = TensorExpression.gamma(d)
momentum = TensorName.vector("ludy_port::momentum")
k, p1, p2 = (momentum(i, lorentz) for i in range(3))
loop_products = (kk, kp1, kp2)
labels = tuple(S("ludy_port::D")(i) for i in range(1, 7))
report = {
    "source_commit": inputs["source_commit"],
    "stage": "pre-IBP",
    "integrated_coefficients": False,
    "weights_applied": False,
    "processes": {},
    "checks": {},
}
computed = {}

for process in ("dy", "dis") if args.process == "both" else (args.process,):
    # Keep virtualities finite: DIS p.q=(Q²/x-x p²)/2, q²=-Q².
    scalar_products = (kk, kp1, kp2, p1sq, p12, p2sq)
    if process == "dis":
        scalar_products = (kk, kp1, kp2, p1sq, (Qsq / x - x * p1sq) / 2, -Qsq)
    dot_rules = list(
        zip(
            (
                dot(a, b)
                for a, b in ((k, k), (k, p1), (k, p2), (p1, p1), (p1, p2), (p2, p2))
            ),
            scalar_products,
            strict=True,
        )
    )
    process_result = {}
    computed[process] = {}
    for name, spec in inputs[process].items():
        vectors = []
        for coefficients in spec["vectors"]:
            if coefficients is None:
                vectors.append(None)
            else:
                terms = [
                    c * v for c, v in zip(coefficients, (k, p1, p2), strict=True) if c
                ]
                vectors.append(sum(terms[1:], terms[0]))
        source_labels = spec["labels"]
        projections = {"metric": (vectors, source_labels)}
        if process == "dy" and name.startswith("qg "):
            ward_vectors = list(vectors)
            ward_labels = list(source_labels)
            for index, label in enumerate(source_labels):
                if label == "photon":
                    ward_vectors[index] = k + p1 + p2
                    ward_labels[index] = f"ward_{index}"
            projections["ward"] = (ward_vectors, ward_labels)
        elif process == "dis":
            left, right = spec["photon_slots"]
            metric_labels = list(source_labels)
            metric_labels[left] = metric_labels[right] = "photon"
            projections["metric"] = (vectors, metric_labels)
            # Contract the transverse vector directly; Spenso distributes sums.
            transverse = p1 + scalar_products[4] / Qsq * p2
            for label, vector in (("PP", transverse), ("ward", p2)):
                projected = list(vectors)
                projected[left] = projected[right] = vector
                projections[label] = (projected, source_labels)

        numerators = {}
        for projection, (projected, indices) in projections.items():
            chain = E("1")
            for index, (vector, label) in enumerate(
                zip(projected, indices, strict=True)
            ):
                chain *= gamma(
                    index, (index + 1) % len(projected), label
                ).to_expression()
                if vector is not None:
                    chain *= vector(lorentz(label)).to_expression()
            scalar = to_dots(
                simplify_metrics(simplify_gamma(chain).expand()).expand()
            ).expand()
            for source, target in dot_rules:
                scalar = scalar.replace(source, target)
            numerators[projection] = scalar.expand().together().cancel()
        if process == "dis":
            transverse_square = p1sq + scalar_products[4] ** 2 / Qsq
            numerators["longitudinal"] = (
                numerators.pop("PP") / transverse_square
            ).cancel()
            numerators["F1"] = (
                (numerators["metric"] - numerators["longitudinal"]) / (2 - d)
            ).cancel()
            numerators["F2"] = (
                2 * x * (numerators["F1"] + numerators["longitudinal"])
            ).cancel()
            numerators["FL"] = (2 * x * numerators["longitudinal"]).cancel()
            assert (
                numerators["F2"] - 2 * x * numerators["F1"] - numerators["FL"]
            ).expand().cancel() == 0
        computed[process][name] = numerators

        # Ordered momentum/mass inputs determine the scalar propagators.
        if process == "dy":
            propagators = spec["propagators"]
            powers = (1, 1, *spec["powers"])
        else:
            propagators = [
                {"momentum": coefficients, "mass_squared": False}
                for coefficients in (
                    (1, 0, 0),
                    (1, 1, 0),
                    (1, 1, 1),
                    (1, 0, 1),
                    (0, 1, 0),
                    (0, 1, 1),
                )
            ]
            powers = spec["powers"]
        denominators = []
        for propagator in propagators:
            terms = [
                c * v
                for c, v in zip(propagator["momentum"], (k, p1, p2), strict=True)
                if c
            ]
            vector = sum(terms[1:], terms[0])
            denominator = dot(vector, vector) - (
                msq if propagator["mass_squared"] else 0
            )
            for source, target in dot_rules:
                denominator = denominator.replace(source, target)
            denominators.append(denominator.expand())
        loop_lines = tuple(
            i for i, prop in enumerate(propagators) if prop["momentum"][0]
        )
        # An absent fourth DIS propagator need not appear as a redundant coordinate.
        if process == "dis" and powers[3] == 0:
            loop_lines = loop_lines[:3]
        loop_denominators = [denominators[i] for i in loop_lines]
        base_powers = [powers[i] for i in loop_lines]
        external_factor = prod(
            den**power
            for i, (den, power) in enumerate(zip(denominators, powers, strict=True))
            if i not in loop_lines
        )
        if process == "dy":
            # LUDY keeps the two incoming virtuality carriers outside its loop kernel.
            external_factor = (p1sq + 2 * p12 + p2sq) ** spec["divisor_power"]
        kernel = E("1") / prod(
            den**power
            for den, power in zip(loop_denominators, base_powers, strict=True)
        )
        children = [(tuple(range(3)), E("1"), base_powers)]
        if len(loop_lines) == 4:
            apart = kernel.apart(*loop_products)
            assert (apart - kernel).together().cancel() == 0, name
            children = []
            for term in apart.terms():
                matches = []
                for basis in combinations(range(4), 3):
                    coefficient = (
                        (term * prod(loop_denominators[i] for i in basis))
                        .together()
                        .cancel()
                    )
                    if not any(coefficient.contains(v) for v in loop_products):
                        matches.append((basis, coefficient, [1, 1, 1]))
                assert len(matches) == 1, (name, term)
                children.extend(matches)

        reductions = {}
        for projection in ("metric",) if process == "dy" else ("F1", "F2", "FL"):
            reduced_terms = []
            reconstruction = E("0")
            for basis, coefficient, indices in children:
                basis_denominators = [loop_denominators[i] for i in basis]
                local_labels = labels[: len(basis)]
                solutions = Expression.solve(
                    [
                        den - label
                        for den, label in zip(
                            basis_denominators, local_labels, strict=True
                        )
                    ],
                    loop_products,
                )
                assert len(solutions) == 1, name
                reduced = numerators[projection]
                for variable in loop_products:
                    reduced = reduced.replace(variable, solutions[0][variable])
                polynomial = reduced.expand().to_polynomial(vars=list(local_labels))
                for degrees, value in polynomial.coefficient_list(list(local_labels)):
                    shifted = [
                        power - int(degree)
                        for power, degree in zip(indices, degrees, strict=True)
                    ]
                    factor = (
                        coefficient * value.to_expression() / external_factor
                    ).cancel()
                    reconstruction += factor / prod(
                        den**power
                        for den, power in zip(basis_denominators, shifted, strict=True)
                    )
                    reduced_terms.append(
                        {
                            "lines": [loop_lines[i] + 1 for i in basis],
                            "indices": shifted,
                            "coefficient": factor.format_plain(),
                        }
                    )
            assert (
                reconstruction - numerators[projection] * kernel / external_factor
            ).together().cancel() == 0, (name, projection)
            reductions[projection] = reduced_terms

        # Keep the iterated channel memberships; an overlapping line occurs once
        # in the union. Individual CutPropagators do not solve the joint cut.
        sectors = []
        cuts = (
            spec["cuts"].items()
            if process == "dy"
            else ((cut["name"], cut) for cut in spec["cuts"])
        )
        for cut_name, cut in cuts:
            cut_lines = (
                [i + 1 for i, selected in enumerate(cut["mask"]) if selected]
                if process == "dy"
                else sorted(set(cut["initial"] + cut["final"]))
            )
            active = process == "dis" or E(cut["weight"]) != 0
            distributions = []
            if active:
                for line in cut_lines:
                    power = powers[line - 1]
                    assert power > 0, (name, cut_name, line)
                    energy = S("ludy_port::q0")(line)
                    on_shell = S("ludy_port::omega")(line)
                    distributions.append(
                        fk.CutPropagator(energy, on_shell, power=power)
                        .to_expression()
                        .format_plain()
                    )
            supported = {}
            loop_cuts = set(cut_lines) & {i + 1 for i in loop_lines}
            for projection, terms in reductions.items():
                supported[projection] = (
                    [
                        i
                        for i, term in enumerate(terms)
                        if all(
                            dict(zip(term["lines"], term["indices"], strict=True)).get(
                                line, 0
                            )
                            > 0
                            for line in loop_cuts
                        )
                    ]
                    if active
                    else []
                )
            sectors.append(
                {
                    "name": cut_name,
                    "input": cut,
                    "lines": cut_lines,
                    "active": active,
                    "formal_positive_energy_cuts": distributions,
                    "supported_terms": supported,
                    "joint_cut_evaluated": False,
                }
            )
        process_result[name] = {
            "physical_input": spec,
            "numerators": {
                key: value.format_plain() for key, value in numerators.items()
            },
            "denominators": [value.format_plain() for value in denominators],
            "powers": list(powers),
            "reductions": reductions,
            "sectors": sectors,
        }
    report["processes"][process] = process_result

if "dy" in computed:
    # The on-shell qg Ward cancellation holds on the common two-particle cut.
    wards = {}
    ward_integrands = {}
    for piece in ("box", "triangle", "bubble"):
        name = f"qg {piece}"
        value = -computed["dy"][name]["ward"].replace(p1sq, 0).replace(p2sq, 0) / (
            msq * (2 * p12) ** inputs["dy"][name]["divisor_power"]
        )
        wards[piece] = value.together().cancel()
        spec = report["processes"]["dy"][name]
        denominator = (
            E(spec["denominators"][3], default_namespace="ludy_port")
            .replace(p1sq, 0)
            .replace(p2sq, 0)
            ** spec["powers"][3]
        )
        ward_integrands[piece] = (
            (value / denominator)
            .together()
            .cancel()
            .replace(kk, 0)
            .replace(kp2, msq / 2 - p12 - kp1)
        )
    residual = (
        (
            2 * ward_integrands["triangle"]
            - ward_integrands["box"]
            - ward_integrands["bubble"]
        )
        .together()
        .cancel()
    )
    assert residual == 0, residual
    report["checks"]["dy_qg_ward_identity"] = True

# Check the full raised-cut residue, including the uncut opposite-energy factor.
q0, omega, a, b = S(
    "ludy_port::energy", "ludy_port::omega", "ludy_port::a", "ludy_port::b"
)
coefficient = (a + b * q0 + q0**3) / (q0 + 3 * omega)
for power in (1, 2, 3):
    for orientation in (-1, 1):
        residue = coefficient / (q0 + orientation * omega) ** power
        for _ in range(power - 1):
            residue = residue.derivative(q0)
        expected = (
            orientation
            * residue.replace(q0, orientation * omega)
            / factorial(power - 1)
        )
        actual = fk.CutPropagator(
            q0, omega, power=power, orientation=orientation, normalization=1
        ).apply(coefficient, q0)
        assert (actual - expected).together().cancel() == 0, (power, orientation)
report["checks"]["raised_cut_residues"] = 6

if args.check:
    # Reference values are never inputs to the native computation above.
    oracles = json.loads((directory / "oracles.json").read_text())
    assert oracles["source_commit"] == inputs["source_commit"]
    checked = 0
    for process, graphs in computed.items():
        for name, numerators in graphs.items():
            expected = (
                {"metric": oracles[process][name]}
                if process == "dy"
                else oracles[process][name]
            )
            for projection, value in expected.items():
                assert (
                    numerators[projection] - E(value, default_namespace="ludy_port")
                ).together().cancel() == 0, (process, name, projection)
                checked += 1
    if "dy" in computed:
        for piece, value in wards.items():
            assert (
                value - E(oracles["dy_ward"][piece], default_namespace="ludy_port")
            ).together().cancel() == 0, piece
            checked += 1
    report["checks"]["reference_numerators"] = checked

if args.output:
    args.output.write_text(json.dumps(report, indent=2) + "\n")
for process, graphs in report["processes"].items():
    print(
        f"{process.upper()}: {len(graphs)} graphs, {sum(len(g['sectors']) for g in graphs.values())} sectors, exact scalar-integrand reconstruction"
    )
print(json.dumps(report["checks"], sort_keys=True))
print(
    "Pre-IBP only: integrated NLO coefficients and joint iterated cuts remain unimplemented."
)
