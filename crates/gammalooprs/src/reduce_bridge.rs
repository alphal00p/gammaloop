//! Glue between a gammaloop one-loop graph and the `oneloop` IBP reducer.

use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use oneloop::masters::MasterIntegral;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use symbolica::atom::{Atom, AtomCore};
use symbolica::symbol;

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
    let terms: Vec<String> = reduction
        .terms
        .iter()
        .map(|(coeff, master)| format!("({}) {}", show(coeff), format_master(master, &show)))
        .collect();
    ReduceOutcome::Reduced(terms.join(" + ").replace('"', "\\\""))
}

/// Format a master integral as `A0/B0/C0/D0(args...)`
fn format_master(master: &MasterIntegral, show: &dyn Fn(&Atom) -> String) -> String {
    match master {
        MasterIntegral::Tadpole { m_sq } => format!("A0({})", show(m_sq)),
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => {
            format!("B0({}; {}, {})", show(p_sq), show(m1_sq), show(m2_sq))
        }
        MasterIntegral::Triangle {
            p1_sq,
            p2_sq,
            p12_sq,
            m1_sq,
            m2_sq,
            m3_sq,
        } => format!(
            "C0({}, {}, {}; {}, {}, {})",
            show(p1_sq),
            show(p2_sq),
            show(p12_sq),
            show(m1_sq),
            show(m2_sq),
            show(m3_sq)
        ),
        MasterIntegral::Box {
            p1_sq,
            p2_sq,
            p3_sq,
            p4_sq,
            s,
            t,
            m1_sq,
            m2_sq,
            m3_sq,
            m4_sq,
        } => format!(
            "D0({}, {}, {}, {}, {}, {}; {}, {}, {}, {})",
            show(p1_sq),
            show(p2_sq),
            show(p3_sq),
            show(p4_sq),
            show(s),
            show(t),
            show(m1_sq),
            show(m2_sq),
            show(m3_sq),
            show(m4_sq)
        ),
    }
}
