"""Generated massive-electron Furry cancellations at generic off-shell momenta.

References:
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Ga
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Ga-GaGa
https://feyncalc.github.io/FeynCalcExamples/QED/OneLoop/Ga-GaGaGaGa

One photon gives a nonzero odd vacuum integrand whose tensor integral vanishes.
Three and five photons cancel in orientation pairs after verified loop-momentum
changes of variables. The independent polarization probes impose no on-shell,
transversality, or special-dimensional relations, so they test the full tensors.
All native graph factors are retained; the tadpole's zero alone does not test
its individual symmetry-factor normalization.
"""

from pathlib import Path
from time import monotonic

from symbolica import E, S
from symbolica.community import hep
from symbolica.community.spenso import TensorExpression

model = hep.Model(Path(__file__).parents[2] / "feynkit-model/tests/fixtures/sm.json")
electron, photon = (model.particle_by_pdg(p) for p in (11, 22))
vertices = [
    v
    for v in model.vertex_rules
    if sorted(v.particles) == sorted([electron.name, electron.antiname, photon.name])
]
assert len(vertices) == 1
D, mass, charge = S("furry::D", "UFO::Me", "UFO::ee")
K, P, mink = S("gammalooprs::K", "gammalooprs::P", "spenso::mink")
index, wave = S("index_", "wave_")
results = {}
for count in (1, 3, 5):
    start = monotonic()
    # A one-point tadpole needs unrestricted bridge counting. Retain selfloops
    # and zero external flow, and disable optional topology/numerator filters.
    result = hep.Process(model, [photon], [photon] * (count - 1)).generate_diagrams(
        loops=1,
        max_vertices=count,
        allow_self_loops=True,
        allow_zero_flow_edges=True,
        maximum_bridges=None,
        vertex_allow=vertices,
        self_energy=None,
        tadpoles=None,
        zero_snails=None,
        numerator_grouping=None,
        progress=None,
    )
    assert len(result.diagrams) == {1: 1, 3: 2, 5: 24}[count]
    # P(100+i) are independent test vectors, distinct from physical P(0..n-2).
    # Contract each external Lorentz slot once; vanishing for arbitrary probes
    # establishes vanishing of the entire amputated tensor.
    probes = [P(100 + i) for i in range(count)]
    kinematics = hep.Kinematics(
        D, momenta=[K(0)] + [P(i) for i in range(count - 1)] + probes
    )
    families, traces, weights = [], [], []
    for diagram_index, diagram in enumerate(result.diagrams):
        numerator = model.expand_couplings(
            diagram.numerator_expression().to_expression()
        )
        for edge in diagram.external_edges:
            matches = list(
                diagram.projector_expression().match(
                    wave(edge.id, mink(4, index)), max_level=0
                )
            )
            assert len(matches) == 1
            port = dict(matches[0])[index]
            numerator *= P(100 + edge.external_index, mink(4, port))
        numerator = numerator.replace(mink(4, index), mink(D, index))
        # Trace compact graph-edge momenta first. Expanding their affine loop
        # routings inside a ten-gamma trace would needlessly multiply its terms.
        tensor = (
            TensorExpression(numerator.expand())
            .simplify_gamma()
            .expand()
            .simplify_metrics()
            .to_dots()
        )
        assert tensor.is_scalar
        weight = (
            diagram.overall_factor_expression(evaluate=True)
            * diagram.numerator_prefactor_expression()
        )
        trace = (
            kinematics.apply(
                diagram.momentum_basis().route_expression(tensor.to_expression())
            )
            * weight
        ).expand()
        assert trace != E("0")
        traces.append(trace)
        weights.append(weight)
        raw = diagram.propagator_family(kinematics=kinematics)
        # Include independent polarization probes among the external vectors so
        # the verified scalar map transports K.e_i as well as K.p_i. No added
        # denominator or kinematic constraint is needed for an affine mapping.
        family = hep.IntegralFamily(
            raw.loop_momenta,
            raw.external_momenta + probes,
            raw.denominators,
            kinematics=kinematics,
        )
        assert len(family.denominators) == count
        families.append(family)
        if count == 5 and (diagram_index + 1) % 4 == 0:
            print(f"Five-photon Dirac traces: {diagram_index + 1}/24", flush=True)
    if count == 1:
        mapping = families[0].mapping_to(families[0], [-K(0)])
        assert mapping is not None
        assert (mapping.apply(traces[0]) + traces[0]).expand() == E("0")
        vacuum = hep.TensorReducer(D).with_integrated_vector(K(0, mink(D)))
        assert vacuum.reduce(traces[0]).expand() == E("0")
        pairs = [(0, 0, mapping)]
    else:
        remaining, pairs = set(range(len(families))), []
        while remaining:
            left = min(remaining)
            remaining.remove(left)
            for right in sorted(remaining):
                mapping = families[right].find_mapping(families[left])
                if mapping is not None and (
                    traces[left] + mapping.apply(traces[right])
                ).expand() == E("0"):
                    pairs.append((left, right, mapping))
                    remaining.remove(right)
                    break
            else:
                raise AssertionError(("unpaired orientation", count, left))
        assert len(pairs) == len(result.diagrams) // 2
        assert sorted(i for left, right, _ in pairs for i in (left, right)) == list(
            range(len(result.diagrams))
        )
    for left, right, mapping in pairs:
        # The map is a reflection with unit absolute Jacobian; its sign does
        # not introduce a second fermion-loop sign. It also preserves all
        # unit propagator powers, so numerator cancellation is integrand cancellation.
        assert dict(mapping.momentum_rules)[K(0)].coefficient(K(0)) == E("-1")
        assert mapping.map_powers([1] * count) == [1] * count
        for source_index, target_index in enumerate(mapping.denominator_map):
            assert (
                mapping.apply(families[right].denominators[source_index])
                - families[left].denominators[target_index]
            ).together() == E("0")
    results[count] = (result, families, traces, pairs)
    print(
        f"{count} photons: {len(result.diagrams)} generated diagrams, "
        f"{len(pairs)} verified reflection pair(s), exact zero in {monotonic() - start:.2f}s",
        flush=True,
    )
print(
    "Full massive off-shell one-, three-, and five-photon Furry cancellations passed.",
    flush=True,
)
