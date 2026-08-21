//! Glue between a gammaloop one-loop graph and the `oneloop` IBP reducer.

use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use oneloop::masters::{MasterBasis, OneLoopMasters};
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use symbolica::atom::{Atom, AtomCore};
use symbolica::symbol;

use crate::feyngen::diagram_generator::evaluate_overall_factor;
use crate::graph::Graph;
use crate::graph::lmb::LMBext;
use crate::utils::{GS, W_};

/// The outcome of attempting to reduce a graph's numerator
pub(crate) enum ReduceOutcome {
    /// Reduced to `Σ coeff × master`, formatted for the dot `reduced_num`.
    Reduced(String),
    /// Not a one-loop graph
    NotOneLoop(usize),
    /// One-loop, but the contracted numerator is identically zero.
    ZeroNumerator,
    /// One-loop with a non-zero numerator, but the reducer produced no terms.
    Unsupported,
}

/// Reduce a one-loop graph's numerator to master integrals.
pub(crate) fn reduce_graph_numerator(graph: &Graph, num: &Atom) -> ReduceOutcome {
    // The reducer only handles a single loop; tree/multi-loop short-circuit.
    let n_loops = graph.loop_momentum_basis.loop_edges.len();
    if n_loops != 1 {
        return ReduceOutcome::NotOneLoop(n_loops);
    }

    // Complete the Dirac traces to a scalar
    let scalar = num
        .collect_gamma_chains()
        .simplify_gamma()
        .expand()
        .simplify_metrics();

    // Collapse the graph grouping / symmetry / sign bookkeeping symbols
    // (`NumeratorDependentGrouping`, `AutG`, `InternalFermionLoopSign`, …) into
    // their numeric values, so grouped-graph duplicate terms sum instead of each
    // monomial carrying a ~95-char tag and appearing once per grouped graph.
    let scalar = evaluate_overall_factor(scalar.as_view());

    // Rewrite the edge momenta `Q(edge, ..)` into the loop basis `K(0,..) + P(j,..)`.
    let reps =
        graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[W_.a___]);
    let num_lmb = scalar
        .replace_multiple(&reps)
        .expand()
        .simplify_metrics()
        .expand();

    // Collect the loop propagators with their masses.
    let mut edges = Vec::new();
    for (_, i, _) in graph.iter_edges() {
        let loop_expr = graph
            .loop_momentum_basis
            .loop_atom(i, GS.loop_mom, &[W_.a___], false);
        if loop_expr.is_zero() {
            continue; // external / tree edge, not a loop propagator
        }
        let ext_expr = graph
            .loop_momentum_basis
            .ext_atom(i, GS.external_mom, &[W_.a___], false);
        let mass = graph[i].particle.mass_atom();
        edges.push(oneloop::bridge::GammaloopEdge {
            lmb_rep: loop_expr + ext_expr,
            mass_sq: &mass * &mass,
        });
    }
    if edges.is_empty() {
        // A one-loop graph should have loop propagators; be defensive.
        return ReduceOutcome::NotOneLoop(n_loops);
    }

    let heads = oneloop::bridge::GammaloopHeads {
        loop_mom: GS.loop_mom,
        external_mom: GS.external_mom,
        index: symbol!("spenso::mink"),
        metric: symbol!("spenso::g"),
    };
    let family = oneloop::bridge::family_from_gammaloop(&num_lmb, &edges, &heads);
    let reduction = oneloop::reduce::reduce(&family);
    if reduction.terms.is_empty() {
        return if num_lmb.is_zero() {
            ReduceOutcome::ZeroNumerator
        } else {
            ReduceOutcome::Unsupported
        };
    }

    let show = |a: &Atom| {
        a.printer(SpensoPrintSettings::typst().typst_symbolica())
            .to_string()
    };

    // Fold the reduction into a single atom `Σ coeff · master`, then
    // `collect_factors` to pull the common coupling / colour / polarization
    // prefactor out front (so it appears once, not once per master) and collect
    // like terms in each coefficient. The master heads (`A0`/`B0`/`C0`/`D0`) are
    // opaque function symbols, so they survive the factoring untouched; strip the
    // `oneloop::` namespace for display.
    let basis = OneLoopMasters;
    let amplitude = reduction
        .terms
        .iter()
        .fold(Atom::Zero, |acc, (coeff, master)| {
            acc + coeff * basis.symbol(master)
        });
    let reduced_num = show(&amplitude.collect_factors())
        .replace("oneloop::", "")
        .replace('"', "\\\"");
    ReduceOutcome::Reduced(reduced_num)
}
