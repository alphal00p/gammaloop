use gammalooprs::{
    graph::{Graph, parse::IntoGraph},
    initialisation::test_initialise,
    processes::Amplitude,
    utils::load_generic_model,
    uv::{ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings},
};
use symbolica::atom::AtomCore;

#[test]
fn massive_fermion_bubble_matches_contracted_numerator_after_uv_integration() {
    test_initialise().unwrap();
    let model = load_generic_model("scalars");
    // Each spin numerator belongs to its own edge or vertex. The independent
    // trace identity uses only the two momenta incident to the same vertex:
    // Tr[(q1_slash+m) gamma_mu (q2_slash+m) gamma^mu]
    // = 4*((2-d)*q1.q2 + d*m^2).
    let spin_edges = [
        "spenso::gamma(spenso::bis(4,0),spenso::bis(4,1),gammalooprs::Q(1,spenso::mink(gammalooprs::dim)))+UFO::mass_scalar_1*spenso::g(spenso::bis(4,0),spenso::bis(4,1))",
        "spenso::gamma(spenso::bis(4,2),spenso::bis(4,3),gammalooprs::Q(2,spenso::mink(gammalooprs::dim)))+UFO::mass_scalar_1*spenso::g(spenso::bis(4,2),spenso::bis(4,3))",
    ];
    let spin_vertices = [
        "spenso::gamma(spenso::bis(4,3),spenso::bis(4,0),spenso::mink(gammalooprs::dim,0))",
        "spenso::gamma(spenso::bis(4,1),spenso::bis(4,2),spenso::mink(gammalooprs::dim,0))",
    ];
    let scalar = "4*((2-gammalooprs::dim)*spenso::dot(gammalooprs::Q(1,spenso::mink(gammalooprs::dim)),gammalooprs::Q(2,spenso::mink(gammalooprs::dim)))+gammalooprs::dim*UFO::mass_scalar_1^2)";
    for prescription in [ApproximationType::PolePart, ApproximationType::MUV] {
        let mut results = Vec::new();
        for (edges, vertices) in [(spin_edges, spin_vertices), (["1", "1"], [scalar, "1"])] {
            let graph: Graph = format!(
                r#"digraph fermion_bubble {{
                    overall_factor="1"; projector="1";
                    edge [particle="scalar_1", num="1", mass="0"];
                    node [num="1", dod="0"];
                    incoming [style=invis]; outgoing [style=invis];
                    incoming -> a [id=0];
                    a -> b [id=1, lmb_id=0, mass="UFO::mass_scalar_1", dod="-1", num="{}"];
                    b -> a [id=2, mass="UFO::mass_scalar_1", dod="-1", num="{}"];
                    b -> outgoing [id=3];
                    a [num="{}"];
                    b [num="{}"];
                }}"#,
                edges[0], edges[1], vertices[0], vertices[1],
            )
            .into_graph(&model)
            .unwrap();
            let mut amplitude = Amplitude::from_graph_list("fermion_bubble", vec![graph]).unwrap();
            let result = amplitude.graphs[0]
                .renormalization_part(&UVgenerationSettings {
                    softct: false,
                    renormalization_prescription: RenormalizationPrescriptionSettings {
                        log_divergent: prescription,
                        massive_power_divergent: prescription,
                        massless_power_divergent: prescription,
                        ..Default::default()
                    },
                    ..Default::default()
                })
                .unwrap();
            assert!(!result.is_zero());
            results.push(result.expression);
        }
        assert_eq!(
            results[0].collect_factors(),
            results[1].collect_factors(),
            "{prescription:?}: the integrated spin numerator must match its exact scalar trace"
        );
    }
}
