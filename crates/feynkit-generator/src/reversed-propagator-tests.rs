// The reversed-propagator fallback of `edge_numerator`: a particle without a
// propagator of its own uses its antiparticle's, read on the reversed line
// (slots swapped, momentum -Q(edge)). A model that lists only one propagator
// of a pair must give the same numerators as one that lists both.

/// `model` without the propagators of `particles`. Particles that name a
/// removed propagator lose that reference, so resolution falls back to the
/// antiparticle.
fn without_propagators(model: &Model, particles: &[&str]) -> Model {
    let mut definition: serde_json::Value =
        serde_json::from_str(&model.to_json().unwrap()).unwrap();
    // Covariant cut multiplets require every member's own propagator.
    definition["covariant_cut_multiplets"] = serde_json::json!({});
    let propagators = definition["propagators"].as_array_mut().unwrap();
    let removed = propagators
        .iter()
        .filter(|propagator| particles.contains(&propagator["particle"].as_str().unwrap()))
        .map(|propagator| propagator["name"].as_str().unwrap().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(removed.len(), particles.len(), "{particles:?}");
    propagators.retain(|propagator| !particles.contains(&propagator["particle"].as_str().unwrap()));
    for particle in definition["particles"].as_array_mut().unwrap() {
        if particle["propagator"]
            .as_str()
            .is_some_and(|name| removed.iter().any(|removed| removed == name))
        {
            particle["propagator"] = serde_json::Value::Null;
        }
    }
    Model::from_json(&definition.to_string()).unwrap()
}

/// `model` with the numerators of the named propagators' particles replaced.
fn with_propagator_numerators(model: &Model, numerators: &[(&str, &str)]) -> Model {
    let mut definition: serde_json::Value =
        serde_json::from_str(&model.to_json().unwrap()).unwrap();
    for (particle, numerator) in numerators {
        let propagator = definition["propagators"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|propagator| propagator["particle"] == *particle)
            .unwrap();
        propagator["numerator"] = serde_json::json!(numerator);
    }
    Model::from_json(&definition.to_string()).unwrap()
}

fn fallback(model: &Model, name: &str) -> bool {
    let particle = model
        .particle_by_id(model.particle_id(name).unwrap())
        .unwrap();
    model.particle_propagator(particle).unwrap().unwrap().1
}

fn edge_numerator_atom(model: &Model, name: &str) -> Atom {
    let generator = Generator::new(model.clone());
    let numerator = generator
        .edge_numerator(model.particle_id(name).unwrap(), 3, [6, 7])
        .unwrap();
    NumeratorInterface::export(&numerator).0
}

/// For each particle pair, the twins without either propagator give the
/// reference model's edge numerator for both labels of the line.
fn assert_twin_edge_numerators(reference: &Model, pairs: &[(&str, &str)]) {
    for (particle, antiparticle) in pairs {
        for removed in [particle, antiparticle] {
            let twin = without_propagators(reference, &[removed]);
            assert!(fallback(&twin, removed), "{removed}");
            for label in [particle, antiparticle] {
                assert_eq!(
                    edge_numerator_atom(&twin, label),
                    edge_numerator_atom(reference, label),
                    "label {label}, removed {removed}"
                );
            }
        }
    }
}

#[test]
fn reversed_fermion_propagator_negates_the_edge_momentum() {
    // Massive spin-1/2: S_e-(p) = i(pslash + m)[sink, source] along the e-
    // flow; on an e- edge the e+ template must read i(-qslash' + m)[slot 1,
    // slot 2] with q' = -Q(edge) and slot 1 = the sink.
    let reference = Model::qed();
    assert!(!fallback(&reference, "e-") && !fallback(&reference, "e+"));
    assert_twin_edge_numerators(&reference, &[("e-", "e+"), ("mu-", "mu+"), ("u", "u~")]);

    // A blanket slot swap without the momentum sign gives the slash term the
    // wrong sign: the twin then differs from the reference.
    let twin = without_propagators(&reference, &["e-"]);
    let fixed = edge_numerator_atom(&twin, "e-");
    let swap_only = {
        let generator = Generator::new(twin.clone());
        let (propagator, reversed) = generator
            .propagator_numerator(
                twin.particle_by_id(twin.particle_id("e-").unwrap())
                    .unwrap(),
            )
            .unwrap();
        assert!(reversed);
        let legs = [Flow::Sink, Flow::Source].map(|flow| NumeratorHalfEdge {
            edge: 3,
            flow,
            hedge: [6, 7][usize::from(flow == Flow::Sink)],
            spin: 2,
            color: 1,
        });
        NumeratorInstantiation {
            owner: NumeratorOwner::Edge(3),
            legs: &legs,
            reversed: false,
        }
        .instantiate(&propagator, NumeratorSector::Spin)
        .unwrap()
    };
    assert_ne!(fixed, swap_only);
    let momentum = symbol!("gammalooprs::Q");
    let negate_q = |atom: &Atom| {
        atom.replace_map(|term, _, out| {
            if matches!(term, AtomView::Fun(f) if f.get_symbol() == momentum) {
                **out = -term.to_owned();
            }
        })
        .expand()
    };
    // The swap-only result is exactly the fixed one at -Q(edge).
    assert_eq!(negate_q(&swap_only), fixed.expand());
}

#[test]
fn reversed_vector_scalar_and_ghost_propagators_are_unchanged() {
    // Momentum-even templates: the fallback equals the explicit propagator.
    let reference = Model::standard_model();
    let pairs = [
        ("W+", "W-"),
        ("G+", "G-"),
        ("ghWp", "ghWp~"),
        ("ghWm", "ghWm~"),
        ("ghG", "ghG~"),
    ];
    assert_twin_edge_numerators(&reference, &pairs);

    // Unitary-gauge W with P(index) and with leg-tagged P(index, leg)
    // momenta: both are even, so the sign of Q(edge) must cancel.
    for numerator in [
        "-1𝑖*Metric(idx(1,1),idx(1,2))+1𝑖*P(idx(1,1))*P(idx(1,2))/MW^2",
        "-1𝑖*Metric(idx(1,1),idx(1,2))+1𝑖*P(idx(1,1),1)*P(idx(1,2),2)/MW^2",
    ] {
        let unitary =
            with_propagator_numerators(&reference, &[("W+", numerator), ("W-", numerator)]);
        assert_twin_edge_numerators(&unitary, &[("W+", "W-")]);
    }
}

fn photon_self_energy_qed(model: &Model, loops: usize) -> GenerationResult {
    Process::new([22_i64], [22_i64])
        .with_particle_filters(None, Some(vec!["a".into(), "e-".into()]), None, vec![])
        .unwrap()
        .generate_diagrams(
            model.clone(),
            &GenerationOptions::default()
                .with_loop_count(loops, loops)
                .unwrap()
                .allow_self_loops(true)
                .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible),
        )
        .unwrap()
}

fn electron_annihilation(model: &Model, loops: usize) -> GenerationResult {
    Process::new([11_i64, -11_i64], [22_i64, 22_i64])
        .with_particle_filters(None, Some(vec!["a".into(), "e-".into()]), None, vec![])
        .unwrap()
        .generate_diagrams(
            model.clone(),
            &GenerationOptions::default()
                .with_loop_count(loops, loops)
                .unwrap(),
        )
        .unwrap()
}

/// Numerical samples plus the canonical exact comparison.
fn sampled_identical() -> NumeratorGrouping {
    NumeratorGrouping::Identical(GraphGroupingOptions {
        number_of_numerical_samples: 3,
        fully_numerical_substitution: true,
        check_canonical_numerator: true,
        ..GraphGroupingOptions::default()
    })
}

/// Whether two diagrams of the same topology fall into one `Identical`
/// group (ratio one) under sampling and canonical comparison.
fn sampled_equal(left: &FeynmanDiagram, right: &FeynmanDiagram) -> bool {
    let grouped = crate::grouping::group_diagrams(
        vec![left.clone(), right.clone()],
        left.model(),
        &sampled_identical(),
        false,
        &GenerationOptions::default().threads(1),
    )
    .unwrap();
    grouped.groups.len() == 1
        && grouped.groups[0]
            .members
            .iter()
            .all(|member| member.ratio == Atom::one())
}

/// Every diagram's fermion lines relabelled: e- from a to b becomes e+ from
/// b to a, the same line. Rebuilding the rules then reads each line through
/// the other species' propagator.
fn relabelled(diagram: &FeynmanDiagram) -> FeynmanDiagram {
    let fermions = diagram
        .edges()
        .filter(|(_, _, edge)| {
            edge.external.is_none()
                && diagram
                    .model()
                    .particle_by_id(edge.particle)
                    .unwrap()
                    .is_fermion()
        })
        .map(|(id, _, _)| id)
        .collect::<Vec<_>>();
    let mut diagram = diagram.clone();
    for edge in fermions {
        diagram = diagram.reverse_edge(edge).unwrap();
    }
    diagram.apply_feynman_rules().unwrap()
}

fn fallback_edge_count(diagram: &FeynmanDiagram) -> usize {
    let model = diagram.model();
    diagram
        .edges()
        .filter(|(_, _, edge)| {
            edge.external.is_none()
                && fallback(model, &model.particle_by_id(edge.particle).unwrap().name)
        })
        .count()
}

/// Exact (term by term) and sampled equality of every twin diagram with the
/// reference, for the stored labels and for the relabelled lines.
fn assert_twin_numerators(
    reference: &GenerationResult,
    twin: &GenerationResult,
    label: &str,
) -> usize {
    assert_eq!(reference.diagrams.len(), twin.diagrams.len(), "{label}");
    assert!(!reference.diagrams.is_empty(), "{label}");
    let mut fallback_edges = 0;
    for (left, right) in reference.diagrams.iter().zip(&twin.diagrams) {
        assert_eq!(left.name(), right.name(), "{label}");
        assert_eq!(
            NumeratorInterface::export(left.numerator()).0,
            NumeratorInterface::export(right.numerator()).0,
            "{label} {} exact",
            left.name()
        );
        assert!(
            sampled_equal(left, right),
            "{label} {} sampled",
            left.name()
        );
        let (left_relabelled, right_relabelled) = (relabelled(left), relabelled(right));
        fallback_edges += fallback_edge_count(right) + fallback_edge_count(&right_relabelled);
        assert_eq!(
            NumeratorInterface::export(left_relabelled.numerator()).0,
            NumeratorInterface::export(right_relabelled.numerator()).0,
            "{label} {} relabelled exact",
            left.name()
        );
        assert!(
            sampled_equal(&left_relabelled, &right_relabelled),
            "{label} {} relabelled sampled",
            left.name()
        );
    }
    fallback_edges
}

#[test]
fn twin_qed_models_give_identical_photon_self_energies() {
    // Generated electron loops of L3 include both orientations of the two
    // Furry pairs; relabelling reads every loop through the other species.
    let reference = Model::qed();
    let twins = [
        ("e+ only", without_propagators(&reference, &["e-"])),
        ("e- only", without_propagators(&reference, &["e+"])),
    ];
    for loops in 1..=3 {
        let expected = photon_self_energy_qed(&reference, loops);
        for (name, twin) in &twins {
            let actual = photon_self_energy_qed(twin, loops);
            let fallback_edges =
                assert_twin_numerators(&expected, &actual, &format!("L{loops} {name}"));
            assert!(fallback_edges > 0, "L{loops} {name}: no fallback edge");
        }
    }
}

#[test]
fn twin_qed_models_give_identical_open_electron_lines() {
    let reference = Model::qed();
    let twins = [
        ("e+ only", without_propagators(&reference, &["e-"])),
        ("e- only", without_propagators(&reference, &["e+"])),
    ];
    for loops in 0..=1 {
        let expected = electron_annihilation(&reference, loops);
        for (name, twin) in &twins {
            let actual = electron_annihilation(twin, loops);
            let fallback_edges =
                assert_twin_numerators(&expected, &actual, &format!("ee->aa L{loops} {name}"));
            assert!(
                fallback_edges > 0,
                "ee->aa L{loops} {name}: no fallback edge"
            );
        }
    }
}

#[test]
fn twin_comparison_detects_a_wrong_slash_sign() {
    // Control: an e- propagator with the slash sign flipped (what the
    // swap-only fallback produced on e- lines) must be told apart, exactly on
    // every electron line and by sampling on the open line.
    let reference = Model::qed();
    let electron = reference
        .propagators()
        .iter()
        .find(|propagator| reference.particle_by_id(propagator.particle).unwrap().name == "e-")
        .unwrap()
        .numerator
        .to_string();
    let wrong = with_propagator_numerators(
        &reference,
        &[(
            "e-",
            "(-1*PSlash(idx(1,2),idx(1,1))+Identity(idx(1,2),idx(1,1))*Me)*1𝑖",
        )],
    );
    assert!(electron.contains("PSlash"));
    let expected = electron_annihilation(&reference, 0);
    let actual = electron_annihilation(&wrong, 0);
    for (left, right) in expected.diagrams.iter().zip(&actual.diagrams) {
        assert_ne!(
            NumeratorInterface::export(left.numerator()).0,
            NumeratorInterface::export(right.numerator()).0
        );
        assert!(!sampled_equal(left, right), "{}", left.name());
    }
}
