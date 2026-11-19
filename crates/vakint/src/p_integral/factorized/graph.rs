//! The GRAPH-LEVEL entry (helper VH, numerator assignment, lane
//! rstar-kinds): a FeynKit [`FeynmanDiagram`] (a numerator on every edge and
//! vertex in spenso notation, the model's particles and masses, the diagram's
//! loop-momentum basis) to a value. The caller states nothing beyond the
//! graph; everything else is derived (N/notes/D5-HARD-KIND.md section 3):
//! open indices, external vectors, propagator vs vacuum, the p of a
//! p-integral, masses, routing, colour branches, weights.
//!
//! Dispatch ([`FactorizedNumerator::from_feynkit`] builds the internal
//! representation, one value per colour branch):
//! - PROPAGATOR (an internal edge carries an independent external momentum):
//!   the HAND-OVER PATH, ONE call per graph on the values' own colour branches
//!   (`refpath::direct::lower_entry`, helper VH lane idenso-audit: idenso's
//!   contraction with the caller's projector or the ports open, one
//!   numerator-free denominator match, the direct reader, one reduction).
//!   Each branch enters with a placeholder coefficient that becomes its
//!   `prefactor.total()` (`D → 4 − 2ε`) in the value: colour, n_f, couplings
//!   and open indices stay in the caller's notation. Outside the envelope the
//!   hand-over reproduces (master mode Symbolic or Retain, a massless graph
//!   with one external momentum, a p library) or when it declines, the graph
//!   takes the adapter per branch ([`evaluate_graph_adapter`]);
//! - VACUUM with a massive line: [`super::evaluate_factorized_vacuum`] (idenso,
//!   vakint's vacuum tensor reduction; every line massive: the single-scale
//!   families; ONE massive line among massless ones: the massive closure, cut
//!   and routing derived from the lines; open indices returned);
//! - VACUUM without a massive line: scaleless, exact 0.
//!
//! The internal entries stay the layer under this one. Open indices of a
//! propagator integral are lane d6's (the p path refuses tensor p-integrals
//! today); [`classify`] is the shared detection.

use feynkit_graph::FeynmanDiagram;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};

use super::{
    FactorizedError, FactorizedNumerator, FromFeynkitOptions, ProjectorChoice, Result,
    evaluate_factorized, evaluate_factorized_vacuum,
};
use crate::refpath::{self, ColourBranch, direct};
use crate::utils::vakint_macros::vk_symbol;
use crate::{EvaluationMethod, RustRedMasterMode, Vakint, VakintSettings};

/// What a graph (one colour branch) is, derived from its loop-momentum basis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GraphIntegral {
    /// An internal edge carries an independent external momentum.
    Propagator,
    /// No external momentum through any internal edge, a massive line.
    Vacuum,
    /// No external momentum through any internal edge, no massive line.
    ScalelessVacuum,
}

/// The shared detection (propagator / vacuum) on the internal representation.
pub fn classify(value: &FactorizedNumerator) -> GraphIntegral {
    let loops = value.loops();
    let internal = value.edges.iter().filter(|edge| !edge.is_external());
    let mut massive = false;
    for edge in internal {
        if edge.momentum[loops..].iter().any(|c| *c != 0) {
            return GraphIntegral::Propagator;
        }
        massive |= !edge.mass_squared.is_zero();
    }
    if massive {
        GraphIntegral::Vacuum
    } else {
        GraphIntegral::ScalelessVacuum
    }
}

/// The value of a graph and how it was obtained.
#[derive(Clone, Debug)]
pub struct GraphValue {
    /// Summed over the colour branches, each branch's prefactor applied once
    /// (`Prefactor::total`); open indices in the caller's notation.
    pub value: Atom,
    pub integral: GraphIntegral,
    pub branches: usize,
    pub contract_seconds: f64,
    pub evaluate_seconds: f64,
    /// How the value was obtained.
    pub route: GraphRoute,
}

/// The route [`evaluate_graph`] took.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GraphRoute {
    /// A propagator graph through the hand-over path, one call per graph.
    HandOver,
    /// The adapter per colour branch, and why not the hand-over.
    Adapter(String),
}

/// Evaluate a FeynKit graph (module documentation). `options` keeps the
/// dimension symbol, the colour mode and the projector choice of
/// `from_feynkit`; a graph without external legs is never projected.
pub fn evaluate_graph(
    vakint: &Vakint,
    settings: &VakintSettings,
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &FromFeynkitOptions,
) -> Result<GraphValue> {
    // Lane unknown-masters: strict refusal once, on the sum over the branches.
    if let Some(_boundary) = crate::public_result_boundary() {
        let inner = crate::unknown_master_coefficients_carried(settings);
        let result = evaluate_graph(vakint, &inner, diagram, model, options)?;
        vakint
            .refuse_unknown_master_coefficients(settings, &result.value)
            .map_err(|error| FactorizedError::Evaluation {
                detail: error.to_string(),
            })?;
        return Ok(result);
    }
    let options = entry_options(diagram, options);
    let branches = FactorizedNumerator::from_feynkit(diagram, model, &options)?;
    let why = if branches.is_empty() {
        "every colour branch vanishes".to_string()
    } else if branches
        .iter()
        .any(|branch| classify(branch) != GraphIntegral::Propagator)
    {
        "not a propagator graph".to_string()
    } else {
        match hand_over(vakint, settings, diagram, &branches, options.projector)? {
            Ok(value) => return Ok(value),
            Err(reason) => reason,
        }
    };
    per_branch(vakint, settings, &branches, GraphRoute::Adapter(why))
}

/// The graph through the ADAPTER per colour branch, as [`evaluate_graph`] did
/// before the hand-over path (helper VH): the reference of the tests.
pub fn evaluate_graph_adapter(
    vakint: &Vakint,
    settings: &VakintSettings,
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &FromFeynkitOptions,
) -> Result<GraphValue> {
    // Lane unknown-masters: strict refusal once, on the sum over the branches.
    if let Some(_boundary) = crate::public_result_boundary() {
        let inner = crate::unknown_master_coefficients_carried(settings);
        let result = evaluate_graph_adapter(vakint, &inner, diagram, model, options)?;
        vakint
            .refuse_unknown_master_coefficients(settings, &result.value)
            .map_err(|error| FactorizedError::Evaluation {
                detail: error.to_string(),
            })?;
        return Ok(result);
    }
    let options = entry_options(diagram, options);
    let branches = FactorizedNumerator::from_feynkit(diagram, model, &options)?;
    per_branch(
        vakint,
        settings,
        &branches,
        GraphRoute::Adapter("asked for".into()),
    )
}

fn entry_options(diagram: &FeynmanDiagram, options: &FromFeynkitOptions) -> FromFeynkitOptions {
    let mut options = options.clone();
    if diagram.loop_momentum_basis().external_edges.is_empty() {
        options.projector = ProjectorChoice::None;
    }
    options
}

fn per_branch(
    vakint: &Vakint,
    settings: &VakintSettings,
    branches: &[FactorizedNumerator],
    route: GraphRoute,
) -> Result<GraphValue> {
    let mut value = GraphValue {
        value: Atom::Zero,
        integral: GraphIntegral::ScalelessVacuum,
        branches: branches.len(),
        contract_seconds: 0.0,
        evaluate_seconds: 0.0,
        route,
    };
    for branch in branches {
        let kind = classify(branch);
        value.integral = kind;
        let prefactor = branch.prefactor.total();
        match kind {
            GraphIntegral::Propagator => {
                let got = evaluate_factorized(vakint, settings, branch, &prefactor)?;
                value.value += got.value;
                value.contract_seconds += got.contract_seconds;
                value.evaluate_seconds += got.evaluate_seconds;
            }
            GraphIntegral::Vacuum => {
                let got = evaluate_factorized_vacuum(vakint, settings, branch, &prefactor)?;
                value.value += got.value;
                value.contract_seconds += got.contract_seconds;
                value.evaluate_seconds += got.evaluate_seconds;
            }
            GraphIntegral::ScalelessVacuum => {}
        }
    }
    Ok(value)
}

fn evaluation(detail: impl std::fmt::Display) -> FactorizedError {
    FactorizedError::Evaluation {
        detail: detail.to_string(),
    }
}

/// The placeholder coefficient of colour branch `i` in the hand-over.
fn placeholder(i: usize) -> Atom {
    Atom::var(vk_symbol!(format!("graph_entry_branch_{i}").as_str()))
}

/// A propagator graph's values through the hand-over path (module
/// documentation): `Ok(Err(reason))` when it does not apply (the caller then
/// takes the adapter); errors of the values themselves are the adapter's.
fn hand_over(
    vakint: &Vakint,
    settings: &VakintSettings,
    diagram: &FeynmanDiagram,
    values: &[FactorizedNumerator],
    choice: ProjectorChoice,
) -> Result<std::result::Result<GraphValue, String>> {
    let s = &*refpath::SYMBOLS;
    // The envelope in which the hand-over's reduction and finalisation are
    // those of `Vakint::evaluate` (numerical masters are materialised there
    // after the batch; other backends do not evaluate p-integrals).
    let mode = settings
        .evaluation_order
        .0
        .iter()
        .find_map(|method| match method {
            EvaluationMethod::RustRed(options) => Some(options.master_mode),
            _ => None,
        });
    if !matches!(
        mode,
        Some(RustRedMasterMode::Symbolic | RustRedMasterMode::Retain)
    ) {
        return Ok(Err(format!(
            "RustRed master mode {mode:?}: the hand-over takes Symbolic and Retain"
        )));
    }
    for value in values {
        value.validate()?;
    }
    let first = &values[0];
    if first.externals() != 1 {
        return Ok(Err(format!(
            "{} independent external momenta",
            first.externals()
        )));
    }
    let internal = || first.edges.iter().filter(|edge| !edge.is_external());
    if internal().any(|edge| !edge.mass_squared.is_zero()) {
        return Ok(Err("a massive line".into()));
    }
    let projector = match choice {
        ProjectorChoice::None if first.ports.len() == 2 => direct::EntryProjector::Open,
        ProjectorChoice::None => {
            return Ok(Err(format!("{} open Lorentz ports", first.ports.len())));
        }
        ProjectorChoice::Metric => direct::EntryProjector::Metric,
        ProjectorChoice::Transverse => direct::EntryProjector::Transverse,
    };
    let mut value = GraphValue {
        value: Atom::Zero,
        integral: GraphIntegral::Propagator,
        branches: values.len(),
        contract_seconds: 0.0,
        evaluate_seconds: 0.0,
        route: GraphRoute::HandOver,
    };
    let loops = first.loops();
    if internal().all(|edge| edge.momentum[loops] == 0) {
        // No internal line carries the external momentum: scaleless, exact 0
        // (as the adapter).
        return Ok(Ok(value));
    }
    if vakint.p_library().is_err() {
        return Ok(Err("no p library".into()));
    }
    vakint.validate_settings(settings).map_err(evaluation)?;
    let dimension = Atom::var(first.dimension);
    let terms = values
        .iter()
        .enumerate()
        .map(|(i, value)| ColourBranch {
            coefficient: placeholder(i),
            residuals: vec![
                value
                    .factor_product()
                    .replace(dimension.to_pattern())
                    .with(Atom::var(s.d).to_pattern()),
            ],
        })
        .collect::<Vec<_>>();
    let momentum = (|| -> eyre::Result<Option<Atom>> {
        let kinematics = feynkit_kinematics::Kinematics::in_dimension(&Atom::var(s.d))?;
        let family = diagram.propagator_family(&kinematics)?;
        Ok(match family.external_momenta() {
            [momentum] => Some(momentum.clone()),
            _ => None,
        })
    })();
    let momentum = match momentum {
        Ok(Some(momentum)) => momentum,
        Ok(None) => return Ok(Err("not one independent external momentum".into())),
        Err(error) => return Ok(Err(format!("kinematics: {error}"))),
    };
    let mut seconds = direct::DirectSeconds::default();
    let outcome = settings
        .with_resource_budget(|| {
            Ok(direct::lower_entry(
                vakint,
                settings,
                diagram,
                &terms,
                &momentum,
                projector,
                direct::DirectOptions::default(),
                &mut seconds,
            ))
        })
        .map_err(evaluation)?;
    let projection = match outcome {
        Err(error) => return Ok(Err(format!("hand-over: {error}"))),
        Ok(direct::DirectOutcome::Declined(reason)) => return Ok(Err(reason)),
        Ok(direct::DirectOutcome::Zero(_)) => None,
        Ok(direct::DirectOutcome::Lowered(mut projections)) => Some(projections.remove(0)),
    };
    let mut timed = direct::ReduceSeconds::default();
    if let Some(projection) = projection {
        let outcome = settings
            .with_resource_budget(|| {
                Ok(direct::reduce_projection_timed(
                    vakint,
                    settings,
                    &projection,
                ))
            })
            .map_err(evaluation)?;
        let reduced = match outcome {
            Ok((reduced, seconds)) => {
                timed = seconds;
                reduced
            }
            Err(error) => return Ok(Err(format!("hand-over reduction: {error}"))),
        };
        // The placeholders become the branches' prefactors, `D → 4 − 2ε` as
        // the adapter's input.
        let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let in_epsilon = Atom::num(4) - Atom::num(2) * epsilon;
        let replacements = values
            .iter()
            .enumerate()
            .map(|(i, value)| {
                symbolica::id::Replacement::new(
                    placeholder(i).to_pattern(),
                    value
                        .prefactor
                        .total()
                        .replace(dimension.to_pattern())
                        .with(in_epsilon.to_pattern()),
                )
            })
            .collect::<Vec<_>>();
        value.value = reduced.replace_multiple(&replacements);
        crate::truncation_guard::ensure_known_orders(
            &value.value,
            vk_symbol!(settings.epsilon_symbol.as_str()),
            None,
            "vakint evaluation",
            crate::truncation_guard::plain_name,
        )
        .map_err(evaluation)?;
    }
    value.contract_seconds = seconds.routing + seconds.contract;
    value.evaluate_seconds = seconds.matching
        + seconds.handover
        + seconds.tensor
        + seconds.lower
        + timed.reduce
        + timed.finalise;
    Ok(Ok(value))
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::sync::Arc;

    use super::*;

    /// The graph entry takes the hand-over path on every B -> B fixture of one
    /// to three loops, for the three projector choices (and with merged colour
    /// terms), and its value equals the adapter route's EXACTLY in the
    /// caller's notation (colour invariants, `vakint::nf`, the caller's `D`,
    /// the ports' own index atoms).
    #[test]
    fn the_graph_entry_takes_the_hand_over_and_equals_the_adapter_on_the_bb_fixtures() {
        let mut template = Vakint::new().unwrap();
        if let Some(directory) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
            template
                .load_p_library_directory(std::path::Path::new(&directory), false)
                .unwrap();
        }
        if template.p_library().is_err() {
            eprintln!("skipped: no p library for L1-L3 in this environment");
            return;
        }
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/bb");
        let (mut handed, mut non_zero) = (0, 0);
        for level in ["L1", "L2", "L3"] {
            let directory = root.join(level);
            let model = refpath::load_model(&directory).unwrap();
            let mut files = std::fs::read_dir(&directory)
                .unwrap()
                .filter_map(|e| e.ok().map(|e| e.path()))
                .filter(|p| p.extension().is_some_and(|x| x == "dot"))
                .collect::<Vec<_>>();
            files.sort();
            for file in files {
                let dot = std::fs::read_to_string(&file).unwrap();
                let diagram = FeynmanDiagram::from_dot(Arc::clone(&model), &dot).unwrap();
                let settings = refpath::vakint_settings(diagram.loop_count() as i64 + 1);
                for (projector, merge_colour) in [
                    (ProjectorChoice::Metric, false),
                    (ProjectorChoice::Transverse, false),
                    (ProjectorChoice::None, false),
                    (ProjectorChoice::Metric, true),
                    (ProjectorChoice::None, true),
                ] {
                    let label = format!("{} {projector:?} merge {merge_colour}", file.display());
                    let options = FromFeynkitOptions {
                        projector,
                        merge_colour,
                        ..FromFeynkitOptions::default()
                    };
                    let engine = refpath::fresh_engine(&template).unwrap();
                    let entry = evaluate_graph(&engine, &settings, &diagram, &model, &options)
                        .unwrap_or_else(|e| panic!("{label}: {e}"));
                    let engine = refpath::fresh_engine(&template).unwrap();
                    let adapter =
                        evaluate_graph_adapter(&engine, &settings, &diagram, &model, &options)
                            .unwrap_or_else(|e| panic!("{label}: {e}"));
                    assert_eq!(entry.integral, adapter.integral, "{label}");
                    // Scaleless contact graphs (no internal line carries the
                    // external momentum) are no propagator integral: 0 on both.
                    if adapter.integral == GraphIntegral::Propagator {
                        assert_eq!(entry.route, GraphRoute::HandOver, "{label}");
                        handed += 1;
                    }
                    assert!(
                        (&entry.value - &adapter.value).expand().is_zero(),
                        "{label}: hand-over {} vs adapter {}",
                        entry.value,
                        adapter.value
                    );
                    if !adapter.value.is_zero() {
                        non_zero += 1;
                    }
                }
            }
        }
        eprintln!("graph entry: {handed} hand-overs, {non_zero} non-zero values equal");
        // the 7 fixtures with a non-zero value (as R3's fixture tests) x 5 choices
        assert_eq!((handed, non_zero), (7 * 5, 7 * 5));
    }
}
