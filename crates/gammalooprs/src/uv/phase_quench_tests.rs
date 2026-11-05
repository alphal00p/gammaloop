//! Algebraic UV checks with local scalar numerators and charged fermion cycles.
//!
//! These diagnostics isolate the weight projection; they are not QCD pressure
//! benchmarks. In particular, the nested fixture has a higher-valence local
//! vertex. Every quark line still belongs to an unambiguous oriented cycle.

use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use linnet::half_edge::subgraph::SubSetOps;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;
use three_dimensional_reps::RepresentationMode;
use three_dimensional_reps::thermal::OCCUPATION_FORM;

use crate::cff::PhaseQuenchedDifference;
use crate::dot;
use crate::graph::Graph;
use crate::graph::parse::from_dot::IntoGraph;
use crate::initialisation;
use crate::numerator::soft_energy_certificate::{EnergyZeroOutcome, SoftEnergyAlgebra};
use crate::numerator::symbolica_ext::NumeratorAtomExt;
use crate::settings::global::{GenerationSettings, MediumMode};
use crate::utils;
use crate::utils::{GS, W_};
use crate::uv::approx::projected_4d::Local4dProjectionContext;
use crate::uv::approx::{CutStructure, OrientationProjection};
use crate::uv::hedge_poset::Wood;

#[test]
fn phase_quenched_difference_commutes_with_local_and_integrated_forests() -> Result<()> {
    initialisation::test_initialise()?;
    let graphs: Vec<Graph> = dot!(
        digraph nested {
            edge [num=1]; node [num=1];
            v1 -> v2 -> v1 [particle="u"];
            v3 -> v4 -> v3 [particle="u"];
            v2 -> v3 [particle="g"];
            v1 -> v4 [particle="g"];
            v1 -> v1 [particle="g"];
        }
        digraph disjoint {
            edge [num=1]; node [num=1];
            v1 -> v2 -> v1 [particle="u"];
            v3 -> v4 -> v3 [particle="u"];
            v2 -> v3 [particle="g"];
            v1 -> v4 [particle="g"];
        }
    )?;
    // Unit graph numerators permit scalar denominator algebra here. Never
    // extend this normalization to distributing a nontrivial graph numerator.
    let occupation_coefficients = |expression: &Atom, context: &str| {
        let scalar = expression
            .replace(function!(GS.thermal_weight_wrapper, W_.a_))
            .with(W_.a_)
            .replace_multiple(&*OCCUPATION_FORM)
            .expand();
        let terms = match scalar.as_view() {
            AtomView::Add(sum) => sum.into_iter().collect::<Vec<_>>(),
            term => vec![term],
        };
        let mut coefficients = BTreeMap::<Atom, Atom>::new();
        // The complete two-loop sign average can exceed the production
        // collector's 128-monomial budget. Each expanded scalar addend has one
        // occupation monomial, so collect it independently and merge exactly.
        for term in terms {
            let mut occupations = BTreeSet::new();
            term.visitor(&mut |view| {
                if let AtomView::Fun(call) = view
                    && call.get_symbol() == GS.thermal_distribution
                {
                    occupations.insert(view.to_owned());
                    return false;
                }
                true
            });
            let local = term
                .coefficient_list_exact(&occupations.into_iter().collect::<Vec<_>>())
                .unwrap_or_else(|| {
                    panic!("{context}: scalar addend must be polynomial in occupations: {term}")
                });
            for (occupation, coefficient) in local {
                *coefficients.entry(occupation).or_insert_with(Atom::zero) += coefficient;
            }
        }
        coefficients.into_iter().collect::<Vec<_>>()
    };
    let algebra = SoftEnergyAlgebra::default();
    for mut graph in graphs {
        let difference = PhaseQuenchedDifference::new(&graph)?;
        for medium in [
            MediumMode::ThermodynamicEquilibrium,
            MediumMode::ZeroTemperatureEquilibrium,
        ] {
            let mut projected_results = Vec::new();
            for generate_integrated in [false, true] {
                let mut settings = GenerationSettings::default();
                settings.medium.mode = medium;
                settings.medium.phase_quenched_difference = true;
                settings.uv.softct = false;
                settings.uv.generate_integrated = generate_integrated;
                let options =
                    graph.production_3d_expression_options(&settings, RepresentationMode::Cff)?;
                let contract =
                    graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
                let production =
                    graph.generate_3d_expression_for_integrand(&contract, &None, &options, None)?;
                let mut results = Vec::new();
                for early_projection in [false, true] {
                    let projection = OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &settings.orientation_pattern,
                        true,
                    )
                    .with_phase_quenched_difference(early_projection.then_some(&difference));
                    // Keep exactly the same proper forests on both routes. The
                    // only difference is the stage at which the weights change.
                    let wood = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
                    assert!(
                        wood.graph
                            .iter_nodes()
                            .any(|(_, _, spinney)| spinney.n_components() > 1),
                        "{} must exercise disconnected UV composition",
                        graph.name
                    );
                    let mut forests = wood.unfold();
                    if graph.name == "nested" {
                        assert!(
                            forests.graph.iter_nodes().any(|(_, _, operation)| {
                                operation.key.iter_levels_top_down().count() > 1
                            }),
                            "nested fixture must contain dependent UV operations"
                        );
                    }
                    forests.compute(
                        &mut graph,
                        utils::vakint()?,
                        projection,
                        &settings.uv,
                        true,
                        &mut Local4dProjectionContext::default(),
                    )?;
                    results.push(forests.orientation_parametric_exprs(&graph, &settings.uv)?);
                }
                assert_eq!(results[0].len(), results[1].len());
                let mut nonzero_projection = false;
                for (late, early) in results[0].iter().zip(&results[1]) {
                    let late = late.integrands.resolved()?;
                    let early = early.integrands.resolved()?;
                    assert_eq!(late.iter().count(), early.iter().count());
                    for ((late_key, raw), (early_key, actual)) in late.iter().zip(early.iter()) {
                        assert_eq!(late_key, early_key);
                        let expected = difference.apply(raw)?;
                        let context = format!(
                            "{}: {medium:?}, integrated={generate_integrated}, residue={late_key:?}",
                            graph.name
                        );
                        for (occupation, coefficient) in
                            occupation_coefficients(&(actual - expected), &context)
                        {
                            let certificate = algebra.scalar_zero(&coefficient)?;
                            assert!(
                                matches!(certificate, EnergyZeroOutcome::Zero(_)),
                                "{}: {medium:?}, integrated={generate_integrated}, residue={late_key:?}, occupation={occupation}: {certificate:?}: {coefficient}",
                                graph.name
                            );
                        }
                        if !nonzero_projection {
                            for (_, coefficient) in occupation_coefficients(actual, &context) {
                                if matches!(
                                    algebra.scalar_zero(&coefficient)?,
                                    EnergyZeroOutcome::Nonzero(_)
                                ) {
                                    nonzero_projection = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                assert!(
                    nonzero_projection,
                    "{}: {medium:?}, integrated={generate_integrated}: projection must be nonzero",
                    graph.name
                );
                projected_results.push(results.pop().expect("early projection result exists"));
            }
            // An integrated-on equality alone could hide a vanishing addback.
            // Require a surviving projected integrated contribution as well.
            let mut nonzero_integrated = false;
            assert_eq!(projected_results[0].len(), projected_results[1].len());
            for (local, integrated) in projected_results[0].iter().zip(&projected_results[1]) {
                let local = local.integrands.resolved()?;
                let integrated = integrated.integrands.resolved()?;
                assert_eq!(local.iter().count(), integrated.iter().count());
                for ((local_key, local), (integrated_key, integrated)) in
                    local.iter().zip(integrated.iter())
                {
                    assert_eq!(local_key, integrated_key);
                    let context = format!(
                        "{}: {medium:?}, integrated addback, residue={local_key:?}",
                        graph.name
                    );
                    for (_, coefficient) in occupation_coefficients(&(integrated - local), &context)
                    {
                        if matches!(
                            algebra.scalar_zero(&coefficient)?,
                            EnergyZeroOutcome::Nonzero(_)
                        ) {
                            nonzero_integrated = true;
                            break;
                        }
                    }
                }
            }
            if graph.name == "nested" {
                // The gluon tadpole leaves both charged cycles in its cograph.
                // The disjoint fixture's bubble counterterms instead remove a
                // whole charged cycle, so their projection may vanish.
                assert!(
                    nonzero_integrated,
                    "{}: {medium:?}: integrated UV must exercise a nonzero projected addback",
                    graph.name
                );
            }
        }
    }
    Ok(())
}
