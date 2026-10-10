//! Lifecycle and error boundaries of derived diagram metadata.

use super::*;
use crate::tests::{directed_fermion_line, one_loop, scalar_model};
use linnet::half_edge::subgraph::BaseSubgraph;

fn assert_cold(diagram: &FeynmanDiagram) {
    assert!(diagram.id.get().is_none());
    assert!(
        matches!(&diagram.numerator, DiagramNumerator::Derived(value) if value.get().is_none())
    );
}

#[test]
fn validation_and_metadata_preserve_cold_derived_values() {
    let diagram = one_loop();
    assert_cold(&diagram);
    diagram.validate().unwrap();
    assert_cold(&diagram);
    let clone = diagram.clone();
    assert_cold(&clone);
    let mapped = diagram
        .map_data(
            |_, vertex| {
                let mut vertex = vertex.clone();
                vertex.name.push_str("_display");
                vertex
            },
            |_, _, edge| {
                let mut edge = edge.clone();
                if let Some(leg) = &mut edge.external {
                    leg.name.push_str("_display");
                }
                edge
            },
        )
        .unwrap()
        .with_name("")
        .with_overall_factor(Atom::num(7))
        .with_projector(Atom::num(3));
    mapped.validate().unwrap();
    assert_cold(&mapped);
    assert_eq!(mapped.name(), Some(""));

    let id = diagram.id().unwrap();
    assert_eq!(diagram.id.get(), Some(&id));
    assert!(
        matches!(&diagram.numerator, DiagramNumerator::Derived(value) if value.get().is_none())
    );
    assert_eq!(diagram.numerator(), &Atom::one());
    let warm = diagram.clone();
    assert_eq!(warm.id.get(), Some(&id));
    assert!(matches!(&warm.numerator, DiagramNumerator::Derived(value) if value.get().is_some()));
    assert_cold(&clone);
    let mapped = warm
        .map_data(|_, vertex| vertex.clone(), |_, _, edge| edge.clone())
        .unwrap();
    assert_eq!(mapped.id.get(), Some(&id));
    assert!(matches!(&mapped.numerator, DiagramNumerator::Derived(value) if value.get().is_some()));
}

#[test]
fn concurrent_first_reads_agree_without_forcing_each_other() {
    let diagram = one_loop().with_numerator(Atom::num(6)).unwrap();
    assert_cold(&diagram);
    let values = std::thread::scope(|scope| {
        (0..8)
            .map(|_| scope.spawn(|| (diagram.id().unwrap(), diagram.numerator().clone())))
            .collect::<Vec<_>>()
            .into_iter()
            .map(|thread| thread.join().unwrap())
            .collect::<Vec<_>>()
    });
    assert!(values.iter().all(|value| *value == values[0]));
    assert_eq!(values[0].1, Atom::num(6));
    assert_eq!(diagram.id.get(), Some(&values[0].0));
}

#[test]
fn structural_maps_invalidate_only_the_structural_cache() {
    let diagram = one_loop();
    diagram.id().unwrap();
    diagram.numerator();
    let edits: &[fn(&mut DiagramEdge)] = &[
        |edge| edge.directed = !edge.directed,
        |edge| edge.is_dummy = !edge.is_dummy,
        |edge| edge.source_slot.0 += 1,
        |edge| edge.target_slot.0 += 1,
        |edge| edge.external.as_mut().unwrap().index += 2,
        |edge| edge.external.as_mut().unwrap().connection += 2,
        |edge| edge.external.as_mut().unwrap().state = ExternalState::Outgoing,
    ];
    for edit in edits {
        let mapped = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(0) {
                        edit(&mut edge);
                    }
                    edge
                },
            )
            .unwrap();
        assert!(mapped.id.get().is_none());
        assert!(
            matches!(&mapped.numerator, DiagramNumerator::Derived(value) if value.get().is_some())
        );
        assert_eq!(mapped.id().unwrap(), mapped.structural_id().unwrap());
    }
    let mapped = diagram
        .map_data(
            |_, vertex| DiagramVertex {
                interaction: None,
                ..vertex.clone()
            },
            |_, _, edge| edge.clone(),
        )
        .unwrap();
    assert!(mapped.id.get().is_none());
    assert_ne!(mapped.id().unwrap(), diagram.id().unwrap());

    let line = directed_fermion_line();
    let id = line.id().unwrap();
    let anti = line.model.particle_id("f~").unwrap();
    let mapped = line
        .map_data(
            |_, vertex| vertex.clone(),
            |_, _, edge| DiagramEdge {
                particle: anti,
                ..edge.clone()
            },
        )
        .unwrap();
    assert!(mapped.id.get().is_none());
    assert_ne!(mapped.id().unwrap(), id);
}

#[test]
fn fragment_maps_preserve_the_original_aggregate_until_explicit_replacement() {
    let diagram = one_loop();
    assert_cold(&diagram);
    let mapped = diagram
        .map_data(
            |id, vertex| DiagramVertex {
                numerator: if id == VertexId(0) {
                    Atom::num(7)
                } else {
                    vertex.numerator.clone()
                },
                ..vertex.clone()
            },
            |_, _, edge| edge.clone(),
        )
        .unwrap();
    assert!(mapped.id.get().is_none());
    assert!(matches!(&mapped.numerator, DiagramNumerator::Explicit(value) if value.is_one()));
    assert!(matches!(
        mapped.validate(),
        Err(DiagramError::NumeratorFragmentMismatch)
    ));
    let repaired = mapped.with_numerator_fragments([], []).unwrap();
    assert_cold(&repaired);
    assert_eq!(repaired.numerator(), &Atom::num(7));
}

#[test]
fn unknown_model_references_fail_before_an_id_is_requested() {
    let model = scalar_model();
    let missing_particle = ParticleId::from_index(usize::MAX);
    let missing_rule = VertexRuleId::from_index(usize::MAX);
    for directed in [false, true] {
        let mut builder = FeynmanDiagramBuilder::new(model.clone(), None);
        let vertex = builder.add_vertex(DiagramVertex::interaction(
            "v",
            model.vertex_rule_id("V_1").unwrap(),
        ));
        builder
            .add_edge(vertex, vertex, DiagramEdge::new(missing_particle, directed))
            .unwrap();
        assert!(matches!(
            builder.build(),
            Err(DiagramError::Model(ModelError::NotFound { .. }))
        ));
        let diagram = one_loop();
        assert!(matches!(
            diagram.map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    DiagramEdge {
                        particle: missing_particle,
                        directed,
                        ..edge.clone()
                    }
                }
            ),
            Err(DiagramError::Model(ModelError::NotFound { .. }))
        ));
    }
    let mut builder = FeynmanDiagramBuilder::new(model, None);
    builder.add_vertex(DiagramVertex::interaction("v", missing_rule));
    assert!(matches!(
        builder.build(),
        Err(DiagramError::Model(ModelError::NotFound { .. }))
    ));
    assert!(matches!(
        one_loop().map_data(
            |_, vertex| DiagramVertex {
                interaction: Some(missing_rule),
                ..vertex.clone()
            },
            |_, _, edge| edge.clone(),
        ),
        Err(DiagramError::Model(ModelError::NotFound { .. }))
    ));
}

#[test]
fn a_failed_structural_id_is_returned_without_poisoning_the_cache() {
    let mut diagram = one_loop();
    let interaction = diagram.graph[NodeIndex(0)].interaction;
    diagram.graph[NodeIndex(0)].interaction = Some(VertexRuleId::from_index(usize::MAX));
    for _ in 0..2 {
        assert!(matches!(
            diagram.id(),
            Err(DiagramError::Model(ModelError::NotFound { .. }))
        ));
        assert!(diagram.id.get().is_none());
    }
    diagram.graph[NodeIndex(0)].interaction = interaction;
    assert_eq!(diagram.id().unwrap(), one_loop().id().unwrap());
}

#[test]
fn optional_names_and_derived_values_round_trip_at_serialization_boundaries() {
    let unnamed = FeynmanDiagramBuilder::new(scalar_model(), None)
        .build()
        .unwrap();
    assert_eq!(unnamed.name(), None);
    assert_cold(&unnamed);
    for name in [
        None,
        Some(""),
        Some("null"),
        Some("a\"b\\c"),
        Some("trailing\\"),
        Some("line\nbreak"),
    ] {
        let mut diagram = one_loop().with_numerator(Atom::num(6)).unwrap();
        diagram.name = name.map(str::to_owned);
        assert_cold(&diagram);
        let json = diagram.to_json().unwrap();
        let serialized: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(serialized["name"], serde_json::to_value(name).unwrap());
        let dot = diagram.to_dot().unwrap();
        let parsed: DotGraph = DotGraph::from_string(&dot).unwrap();
        assert_eq!(
            parsed.global_data.statements["feynkit_name"],
            serde_json::to_string(&name).unwrap()
        );
        for restored in [
            FeynmanDiagram::from_json(diagram.model_arc(), &json).unwrap(),
            FeynmanDiagram::from_dot(diagram.model_arc(), &dot).unwrap(),
        ] {
            assert_eq!(restored.name(), name);
            assert_eq!(restored.id().unwrap(), diagram.id().unwrap());
            assert_eq!(restored.numerator(), &Atom::num(6));
            assert_eq!(
                restored.loop_momentum_basis(),
                diagram.loop_momentum_basis()
            );
        }
        let mut bad_numerator = diagram.serde_view().unwrap();
        bad_numerator.numerator = "7".into();
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&bad_numerator).unwrap()
            ),
            Err(DiagramError::NumeratorFragmentMismatch)
        ));
        assert!(matches!(
            FeynmanDiagram::from_dot(
                diagram.model_arc(),
                &dot.replacen("numerator=\"6\"", "numerator=\"7\"", 1)
            ),
            Err(DiagramError::NumeratorFragmentMismatch)
        ));
        let mut bad_id = diagram.serde_view().unwrap();
        bad_id.id = DiagramId(0);
        assert!(matches!(
            FeynmanDiagram::from_json(
                diagram.model_arc(),
                &serde_json::to_string(&bad_id).unwrap()
            ),
            Err(DiagramError::DiagramIdMismatch { .. })
        ));
    }
}

#[test]
fn reversal_and_excision_preserve_lazy_or_explicit_aggregate_semantics() {
    let diagram = directed_fermion_line()
        .with_numerator(Atom::num(5))
        .unwrap();
    let reversed = diagram.reverse_edge(EdgeId(0)).unwrap();
    assert_cold(&diagram);
    assert_cold(&reversed);
    let restored = reversed.reverse_edge(EdgeId(0)).unwrap();
    assert_cold(&restored);
    assert_eq!(restored.id().unwrap(), diagram.id().unwrap());
    assert_eq!(restored.numerator(), &Atom::num(5));
    let explicit =
        FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
    let reversed = explicit.reverse_edge(EdgeId(0)).unwrap();
    assert!(
        matches!(&reversed.numerator, DiagramNumerator::Explicit(value) if *value == Atom::num(5))
    );

    let mut bubble = one_loop().with_numerator(Atom::num(3)).unwrap();
    bubble.name = None;
    let full = bubble.graph.full_filter();
    let whole = bubble.excise(&full, &[]).unwrap();
    assert_cold(&bubble);
    assert_cold(&whole);
    bubble.id().unwrap();
    bubble.numerator();
    let whole = bubble.excise(&full, &[]).unwrap();
    assert_eq!(whole.id.get(), bubble.id.get());
    assert!(matches!(&whole.numerator, DiagramNumerator::Derived(value) if value.get().is_some()));
    let crown = SuBitGraph::from_hedge_iter(
        bubble.graph.iter_crown(NodeIndex(0)),
        bubble.graph.n_hedges(),
    );
    let part = bubble.excise(&crown, &[]).unwrap();
    assert_eq!(part.name(), None);
    assert_cold(&part);
    assert_eq!(part.numerator(), &Atom::num(3));
}

#[test]
fn cut_and_threshold_edits_invalidate_identity_without_demanding_aggregates() {
    let diagram = crate::tests::cut_scalar_bubble();
    let id = diagram.id().unwrap();
    let same = diagram.clone().with_cuts(diagram.cuts().to_vec()).unwrap();
    assert_cold(&same);
    assert_eq!(same.id().unwrap(), id);
    let changed = diagram
        .with_topology_threshold_candidates(Vec::new())
        .unwrap();
    assert_cold(&changed);
    assert_ne!(changed.id().unwrap(), id);
}

#[test]
fn native_import_accepts_opposite_registration_order_and_rejects_changed_model() {
    use std::io::{Read, Write};
    use std::process::{Command, Stdio};

    const ROLE: &str = "FEYNKIT_NATIVE_IMPORT_FINGERPRINT_ROLE";
    const EXPRESSION: &str = "FEYNKIT_NATIVE_IMPORT_FINGERPRINT_EXPRESSION";
    const MARKER: &str = "NATIVE_DIAGRAM:";
    if let Ok(role) = std::env::var(ROLE) {
        let mut names = [
            "fingerprint_model::x",
            "fingerprint_model::y",
            "fingerprint_model::z",
        ];
        if role == "import" {
            names.reverse();
        }
        for name in names {
            let _ = Atom::parse(name, "fingerprint_model", ParseSettings::default()).unwrap();
        }
        let cyclic = "fingerprint_model::{cyclesymmetric}::c(fingerprint_model::x,fingerprint_model::y,fingerprint_model::z)";
        let reflected = "fingerprint_model::{cyclesymmetric}::c(fingerprint_model::x,fingerprint_model::z,fingerprint_model::y)";
        let expression = if std::env::var(EXPRESSION).unwrap() == "cyclic" {
            cyclic.to_owned()
        } else {
            format!("fingerprint_model::{{antisymmetric}}::a({cyclic},{reflected})")
        };
        // Always import the same original attributed model input. Public UFO
        // export has a separate attribute-roundtrip contract.
        let mut definition: serde_json::Value = serde_json::from_str(include_str!(
            "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap();
        definition["propagators"][0]["numerator"] = expression.clone().into();
        let model = Arc::new(Model::from_json(&definition.to_string()).unwrap());
        let diagram = FeynmanDiagram::from_dot(
            model.clone(),
            r#"digraph test {
            ext [style=invis];
            ext -> v [particle="scalar_0"];
            v -> ext [particle="scalar_0"];
            v -> ext [particle="scalar_0"];
        }"#,
        )
        .unwrap();
        if role == "export" {
            println!("{MARKER}{}", serde_json::to_string(&diagram).unwrap());
        } else {
            assert_eq!(role, "import");
            let mut serialized = String::new();
            std::io::stdin().read_to_string(&mut serialized).unwrap();
            let imported = FeynmanDiagram::from_json(model, &serialized).unwrap();
            assert_eq!(imported.id().unwrap(), diagram.id().unwrap());
            assert_eq!(
                imported.loop_momentum_basis(),
                diagram.loop_momentum_basis()
            );
            assert_eq!(imported.numerator(), diagram.numerator());
            definition["propagators"][0]["numerator"] = expression
                .replace("fingerprint_model::x", "different_model::x")
                .into();
            let different = Model::from_json(&definition.to_string()).unwrap();
            assert!(matches!(
                FeynmanDiagram::from_json(different, &serialized),
                Err(DiagramError::ModelFingerprintMismatch { .. })
            ));
        }
        return;
    }
    for expression in ["cyclic", "nested"] {
        let command = |role| {
            let mut command = Command::new(std::env::current_exe().unwrap());
            command.args(["--exact", "lazy_tests::native_import_accepts_opposite_registration_order_and_rejects_changed_model", "--nocapture"])
                .env(ROLE, role).env(EXPRESSION, expression);
            command
        };
        let exported = command("export").output().unwrap();
        assert!(
            exported.status.success(),
            "export {expression}: {}\n{}",
            String::from_utf8_lossy(&exported.stdout),
            String::from_utf8_lossy(&exported.stderr)
        );
        let stdout = String::from_utf8(exported.stdout).unwrap();
        let serialized = stdout
            .lines()
            .find_map(|line| line.strip_prefix(MARKER))
            .unwrap();
        let mut importer = command("import")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        importer
            .stdin
            .take()
            .unwrap()
            .write_all(serialized.as_bytes())
            .unwrap();
        let imported = importer.wait_with_output().unwrap();
        assert!(
            imported.status.success(),
            "import {expression}: {}\n{}",
            String::from_utf8_lossy(&imported.stdout),
            String::from_utf8_lossy(&imported.stderr)
        );
    }
}
