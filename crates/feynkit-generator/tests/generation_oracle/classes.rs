//! Full physical equivalence of class refinement and particle-only enumeration.

use feynkit_generator::CutAmplitude;
use feynkit_generator::Orders;
use feynkit_generator::ParticleSelector;
use feynkit_generator::Restrict;
use feynkit_generator::VertexSelector;
use num_bigint::BigUint;
use std::collections::BTreeMap;
use std::env;
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write as _;
use std::sync::Arc;

use feynkit_generator::{
    Classes, Grouping, Loops, Naming, NumeratorOptions, Process, Request, Run, TopologyRules, Veto,
};
use feynkit_model::Model;
use num_rational::BigRational;

use super::cases::{Case, amplitude, plain};
use super::census::{Census, Scope};
use super::threads;

pub fn check(case: &Case) -> Census {
    let model = case.model();
    let request = case.request(&model);
    check_request(
        model,
        request,
        case.numerators(),
        case.grouping(),
        case.name,
    )
}

fn check_request(
    model: Model,
    request: Request,
    numerators: NumeratorOptions,
    grouping: Grouping,
    name: &str,
) -> Census {
    let model = Arc::new(model);
    let mut particle_reference: Option<(Arc<Census>, Arc<Census>, _)> = None;
    for mode in [Classes::None, Classes::Auto] {
        let request = request.clone().with_classes(mode);
        let mut reference: Option<(Arc<Census>, Census, Arc<Census>, _)> = None;
        for (thread_index, count) in threads().into_iter().enumerate() {
            let run = Run::default().threads(count);
            let context = format!("{name}: {mode:?} at {count} threads");
            let prepared = request
                .generate(Arc::clone(&model), &run)
                .unwrap_or_else(|error| panic!("{context}: preparation: {error}"));
            assert!(
                prepared.report().completed,
                "{context}: incomplete preparation"
            );
            let counts = request.count(Arc::clone(&model), &run).unwrap();
            assert!(counts.complete, "{context}: incomplete count");
            // Compare physical totals separately from the intentionally different search work.
            let physical_totals = (
                (
                    counts.diagrams.clone(),
                    counts.sum_inv_aut.clone(),
                    counts.physical_weight.clone(),
                ),
                counts
                    .per_loop
                    .iter()
                    .map(|(&loops, row)| {
                        (
                            loops,
                            (
                                row.diagrams.clone(),
                                row.sum_inv_aut.clone(),
                                row.physical_weight.clone(),
                            ),
                        )
                    })
                    .collect::<BTreeMap<_, _>>(),
            );
            assert_eq!(
                counts.diagrams,
                BigUint::from(prepared.len()),
                "{context}: count"
            );
            let inverse: BigRational = prepared
                .iter()
                .map(|diagram| BigRational::new(1.into(), diagram.aut().into()))
                .sum();
            assert_eq!(counts.sum_inv_aut, inverse, "{context}: exact inverse Aut");

            let build_minimal = || {
                let minimal_set = prepared
                    .build(&numerators.clone().with_naming(Naming::None), &run)
                    .unwrap_or_else(|error| panic!("{context}: Minimal numerators: {error}"))
                    .group(&grouping, &run)
                    .unwrap_or_else(|error| panic!("{context}: Minimal grouping: {error}"));
                Census::of(&minimal_set)
            };
            // Cover both cache orders without reconstructing any identical configuration.
            let minimal = (thread_index == 0).then(&build_minimal);
            let built = prepared
                .build(&numerators.clone().with_naming(Naming::prefix("D")), &run)
                .unwrap_or_else(|error| panic!("{context}: numerators: {error}"));
            let ungrouped = Arc::new(Census::of(&built));
            let named = if matches!(grouping, Grouping::None) {
                drop(built);
                Arc::clone(&ungrouped)
            } else {
                let grouped = built
                    .group(&grouping, &run)
                    .unwrap_or_else(|error| panic!("{context}: grouping: {error}"));
                Arc::new(Census::of(&grouped))
            };
            assert!(named.report.completed, "{context}: incomplete grouping");
            let minimal = minimal.unwrap_or_else(build_minimal);
            named.assert_minimal(&minimal, &context);
            drop(prepared);
            named.assert_group_sources(&ungrouped, &context);
            assert_eq!(
                &counts.diagrams,
                &(&ungrouped.report.zero_numerator_count + ungrouped.diagrams.len()),
                "{context}: retained and zero-color diagrams"
            );
            // The zero-color filter removes physical weight that the count still has.
            if !numerators.filter_zero_color {
                let physical_weight: BigRational = ungrouped
                    .diagrams
                    .iter()
                    .map(|diagram| {
                        BigRational::new(
                            diagram.weight.numerator().to_string().parse().unwrap(),
                            diagram.weight.denominator().to_string().parse().unwrap(),
                        )
                    })
                    .sum();
                assert_eq!(
                    counts.physical_weight, physical_weight,
                    "{context}: exact signed count weight"
                );
            }
            if let Ok(directory) = env::var("FEYNKIT_ORACLE_DUMP") {
                fs::create_dir_all(&directory).unwrap();
                let prefix = format!("{directory}/{name}.classes-{mode:?}-t{count}");
                for (census, suffix) in [(&named, "json"), (&ungrouped, "sources.json")] {
                    if suffix == "sources.json" && Arc::ptr_eq(&named, &ungrouped) {
                        fs::copy(format!("{prefix}.json"), format!("{prefix}.{suffix}")).unwrap();
                        continue;
                    }
                    let file = File::create(format!("{prefix}.{suffix}")).unwrap();
                    let mut output = BufWriter::new(file);
                    census.write_exact_dump(&mut output).unwrap();
                    output.flush().unwrap();
                }
            }
            if let Some((previous, previous_minimal, previous_sources, previous_totals)) =
                &reference
            {
                named.assert_same(previous, Scope::Exact, &context);
                minimal.assert_same(previous_minimal, Scope::Exact, &context);
                ungrouped.assert_same(previous_sources, Scope::Exact, &context);
                assert_eq!(
                    &physical_totals, previous_totals,
                    "{context}: per-loop physical totals"
                );
            } else {
                reference = Some((named, minimal, ungrouped, physical_totals));
            }
        }
        let (named, _, sources, totals) = reference.expect("a thread count");
        if let Some((particles, particle_sources, particle_totals)) = &particle_reference {
            assert_eq!(
                &totals, particle_totals,
                "{name}: Auto vs None per-loop physical totals"
            );
            sources.assert_class_equivalent(
                particle_sources,
                (&sources, particle_sources),
                &format!("{name}: Auto vs None sources"),
            );
            named.assert_class_equivalent(
                particles,
                (&sources, particle_sources),
                &format!("{name}: Auto vs None grouped"),
            );
        } else {
            particle_reference = Some((named, sources, totals));
        }
    }
    let (reference, sources, _) = particle_reference.expect("particle-only reference");
    drop(sources);
    Arc::try_unwrap(reference).unwrap()
}

/// The historical unrestricted SM workload; do not substitute the QED-capped oracle case.
#[test]
fn unrestricted_sm_uudd_retains_every_original_rule() {
    check(&amplitude(
        "sm_uudd_L1",
        Model::standard_model,
        &["u", "u~"],
        &["d", "d~"],
        1,
        plain,
    ));
}

#[test]
fn mass_selective_vetoes_preserve_physical_content() {
    for (name, make_model, particle) in [
        ("qed_ase_L2", Model::qed as fn() -> Model, "a"),
        ("qcd_gse_L2", Model::qcd as fn() -> Model, "g"),
    ] {
        for veto in [
            Veto::default(),
            Veto::massive(),
            Veto::massless(),
            Veto::all(),
        ] {
            for kind in 0..3 {
                let model = make_model();
                let request = Process::new([particle], [particle]).amplitude(Loops::total(2..=2));
                let topology = TopologyRules::new()
                    .with_self_loops(true)
                    .with_zero_flow_edges(true);
                let topology = match kind {
                    0 => topology.with_self_energies(veto),
                    1 => topology.with_tadpoles(veto),
                    _ => topology.with_snails(veto),
                };
                check_request(
                    model,
                    request.with_topology(topology),
                    NumeratorOptions::default(),
                    Grouping::None,
                    &format!("{name}-veto-{kind}-{veto:?}"),
                );
            }
        }
    }
}

#[test]
fn restrictions_preserve_original_ids_and_empty_intersections() {
    let model = Model::qed();
    let ids = [11, 22].map(|pdg| {
        ParticleSelector::by_id(&model, model.particle_id_by_pdg(pdg).unwrap()).unwrap()
    });
    let sparse = (1..model.vertex_rules().len())
        .step_by(2)
        .map(|index| model.vertex_rule_id_at(index).unwrap())
        .collect::<Vec<_>>();
    let vertices = sparse
        .iter()
        .map(|&id| VertexSelector::by_id(&model, id).unwrap())
        .collect::<Vec<_>>();
    for (name, restriction) in [
        ("all", Restrict::new()),
        ("pdg", Restrict::new().only_particles([11_i64, 22_i64])),
        ("id", Restrict::new().only_particles(ids)),
        ("exclude", Restrict::new().exclude_particles([13_i64])),
        (
            "empty-particles",
            Restrict::new().only_particles(Vec::<ParticleSelector>::new()),
        ),
        (
            "empty-vertices",
            Restrict::new().only_vertices(Vec::<VertexSelector>::new()),
        ),
        (
            "external-excluded",
            Restrict::new().exclude_particles([22_i64]),
        ),
        ("sparse-rules", Restrict::new().only_vertices(vertices)),
    ] {
        let request = Process::new(["a"], ["a"])
            .amplitude(Loops::total(1..=1))
            .with_topology(TopologyRules::new().with_self_loops(false))
            .with_restrict(restriction);
        if name == "sparse-rules" {
            for mode in [Classes::None, Classes::Auto] {
                let diagrams = request
                    .clone()
                    .with_classes(mode)
                    .generate(model.clone().into(), &Run::default().threads(1))
                    .unwrap();
                assert!(!diagrams.is_empty());
                assert!(
                    diagrams
                        .iter()
                        .flat_map(|diagram| diagram.rules())
                        .all(|rule| sparse.contains(&rule)),
                    "preserve original rule IDs"
                );
            }
        }
        check_request(
            model.clone(),
            request,
            NumeratorOptions::default(),
            Grouping::None,
            &format!("restriction-{name}"),
        );
    }
}

#[test]
fn cut_amplitude_orders_and_unresolved_membership_preserve_class_content() {
    for name in ["xs_qed_aaee_L2", "xs_qed_eeaa_L2"] {
        let case = Case::named(name);
        for left_right in [false, true] {
            for unresolved in [false, true] {
                let model = case.model();
                let request = case.request(&model);
                let cuts = request
                    .cuts()
                    .clone()
                    .with_symmetrize_left_right(left_right)
                    .with_cut_amplitude(CutAmplitude {
                        orders: Orders::new().min("QED", 2),
                        ..CutAmplitude::default()
                    });
                let orders = if unresolved {
                    Orders::new().perturbative("QED", 1)
                } else {
                    Orders::new()
                };
                check_request(
                    model,
                    request.with_cuts(cuts).with_orders(orders),
                    NumeratorOptions::default(),
                    Grouping::None,
                    &format!("{name}-lr-{left_right}-unresolved-{unresolved}"),
                );
            }
        }
    }
}
