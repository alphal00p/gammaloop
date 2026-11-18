//! The GRAPH-LEVEL entry (helper VH, numerator assignment, lane
//! rstar-kinds): a FeynKit [`FeynmanDiagram`] (a numerator on every edge and
//! vertex in spenso notation, the model's particles and masses, the diagram's
//! loop-momentum basis) to a value. The caller states nothing beyond the
//! graph; everything else is derived (N/notes/D5-HARD-KIND.md section 3):
//! open indices, external vectors, propagator vs vacuum, the p of a
//! p-integral, masses, routing, colour branches, weights.
//!
//! Dispatch per colour branch ([`FactorizedNumerator::from_feynkit`] builds
//! the internal representation):
//! - PROPAGATOR (an internal edge carries an independent external momentum):
//!   [`super::evaluate_factorized`] (the adapter: idenso, vakint's p path);
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
use symbolica::atom::Atom;

use super::{
    FactorizedNumerator, FromFeynkitOptions, ProjectorChoice, Result, evaluate_factorized,
    evaluate_factorized_vacuum,
};
use crate::{Vakint, VakintSettings};

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
    let mut options = options.clone();
    if diagram.loop_momentum_basis().external_edges.is_empty() {
        options.projector = ProjectorChoice::None;
    }
    let branches = FactorizedNumerator::from_feynkit(diagram, model, &options)?;
    let mut value = GraphValue {
        value: Atom::Zero,
        integral: GraphIntegral::ScalelessVacuum,
        branches: branches.len(),
        contract_seconds: 0.0,
        evaluate_seconds: 0.0,
    };
    for branch in &branches {
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
