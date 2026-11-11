//! Bare propagators with their original open tensor boundaries.

use eyre::{Result, WrapErr, ensure, eyre};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use symbolica::atom::Atom;
use vakint::Vakint;

use crate::{
    graph::{FeynmanGraph, Graph},
    model::Model,
    uv::{UVgenerationSettings, UltravioletGraph, uv_graph::UVE},
};

use super::{GraphRegion, expression::RStarIntegrand};

/// Integrate a massless, connected 1PI two-point graph without UV subtractions.
///
/// Local vertex and propagator numerators retain their open Lorentz, color and
/// spinor slots. The result uses the original external edge's `Q` momentum,
/// with the two external momenta identified by momentum conservation. Global
/// graph factors and projectors are excluded, as for [`super::RStarPlan`].
///
/// This invokes the physical p-integral backend directly, without building a
/// forest or computing a counterterm. Include `forcer` in the evaluation
/// methods; it currently supports at most four loops. Propagators must have
/// unit powers. Physical infrared poles, if present, remain in the result.
///
/// `number_of_terms` has Vakint's convention: an L-loop integral is retained
/// through epsilon^(number_of_terms-L-1). Thus L+1 includes the finite term.
/// The normalization, additional per-loop factor and regulator symbols are
/// exactly those of `uv_settings`.
pub fn integrate_propagator(
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv_settings: &UVgenerationSettings,
    number_of_terms: i64,
) -> Result<Atom> {
    ensure!(
        number_of_terms > 0,
        "a propagator expansion needs at least one term"
    );
    let boundary = graph.external_filter::<SuBitGraph>();
    RStarIntegrand::physical_external_routing(graph, &boundary)?;
    let support = graph.no_dummy().subtract(&boundary);
    let loops = graph.n_loops(&support);
    ensure!(
        (1..=4).contains(&loops),
        "direct propagator integration currently supports one through four loops"
    );
    ensure!(
        graph.is_connected(&support) && graph.bridges_of(&support).is_empty(),
        "direct propagator integration requires a connected 1PI graph"
    );
    for (_, _, edge) in graph.iter_edges_of(&support) {
        ensure!(
            edge.data.mass_atom(model).is_zero(),
            "direct propagator integration requires massless propagators"
        );
        ensure!(
            edge.data.extra_data.vakint_edge_power.unwrap_or(1) == 1,
            "direct propagator integration currently requires unit propagator powers"
        );
    }
    let marked = graph
        .iter_edges_of(&support)
        .next()
        .map(|(_, edge, _)| edge)
        .ok_or_else(|| eyre!("a propagator graph contains no internal edges"))?;
    let region = GraphRegion {
        remaining: support.clone(),
        support,
        insertions: Vec::new(),
        boundary,
        loops,
    };
    // The marked edge is only a constructor argument here: its auxiliary
    // mass is zero and no Taylor operation is performed.
    let integrand = RStarIntegrand::new(graph, model, &region, &[], marked, Atom::Zero)?;
    let mut settings = uv_settings.vakint.true_settings();
    settings.number_of_terms_in_epsilon_expansion = number_of_terms;
    integrand
        .integrate_physical(vakint, uv_settings, &settings, 0)
        .wrap_err("while integrating a bare open-tensor propagator")
}

#[cfg(test)]
mod tests {
    use eyre::Result;
    use idenso::representations::Bispinor;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::{
            representation::{Minkowski, RepName},
            slot::{DummyAind, ParseableAind},
        },
    };
    use symbolica::{
        atom::{Atom, Symbol},
        function,
    };
    use vakint::Vakint;

    use crate::{
        finalized_runtime_dot,
        graph::{Autogen, Graph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        model::Model,
        numerator::aind::Aind,
        utils::{GS, load_generic_model},
        uv::UVgenerationSettings,
    };

    use super::integrate_propagator;
    use crate::uv::rstar::{RStarPlan, tests::assert_same_counterterm};

    fn settings() -> UVgenerationSettings {
        let mut settings = crate::uv::rstar::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        settings
    }

    fn make_massless(graph: &mut Graph) {
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
    }

    fn bubble(model: &Model) -> Result<Graph> {
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=7];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=2 sink="{ufo_order:0}"];
                B:1 -> e [id=3 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, model
        )?;
        make_massless(&mut graph);
        graph.overall_factor = Atom::num(11);
        graph.global_prefactor.num = Atom::num(13);
        Ok(graph)
    }

    fn external(graph: &Graph) -> EdgeIndex {
        graph
            .iter_edges()
            .find_map(|(pair, edge, _)| (!pair.is_paired()).then_some(edge))
            .unwrap()
    }

    fn external_square(graph: &Graph) -> Atom {
        let momentum = function!(
            GS.emr_mom,
            usize::from(external(graph)),
            Minkowski {}.new_rep(GS.dim).to_symbolic([])
        );
        function!(SPENSO_TAG.dot, &momentum, &momentum)
    }

    fn scalar_bubble_through_finite(graph: &Graph) -> Atom {
        Atom::one() / Atom::var(GS.dim_epsilon)
            - function!(Symbol::LOG, -external_square(graph) / Atom::var(GS.mu_r_sq))
    }

    #[test]
    fn bare_scalar_bubble_keeps_finite_term_and_excludes_global_factors() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = bubble(&model)?;
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        assert_same_counterterm(&actual, &scalar_bubble_through_finite(&graph))
    }

    #[test]
    fn bare_bubble_retains_open_tensor_and_evanescent_finite_terms() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let nu = minkowski.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let numerator = function!(GS.emr_mom, usize::from(loop_edge), &mu)
            * function!(GS.emr_mom, usize::from(loop_edge), &nu);
        graph.underlying[loop_edge].num = Autogen::explicit(numerator);
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        let scalar = scalar_bubble_through_finite(&graph);
        let q_mu = function!(GS.emr_mom, usize::from(external(&graph)), &mu);
        let q_nu = function!(GS.emr_mom, usize::from(external(&graph)), &nu);
        // Tensor bubble = [d q_mu q_nu - q² g_mu_nu] B / [4(d-1)].
        let expected = q_mu * q_nu * (&scalar / Atom::num(3) + Atom::one() / Atom::num(18))
            - external_square(&graph)
                * ETS.metric(&mu, &nu)
                * (&scalar / Atom::num(12) + Atom::one() / Atom::num(18));
        assert_same_counterterm(&actual, &expected)
    }

    #[test]
    fn bare_bubble_projects_loop_slash_inside_open_clifford_word() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let spinor = Bispinor {}.new_rep(4);
        let left = spinor.to_symbolic([Aind::Hedge(40, 0).to_atom()]);
        let right = spinor.to_symbolic([Aind::Hedge(41, 0).to_atom()]);
        let mu = minkowski.to_symbolic([Aind::new_dummy().to_atom()]);
        let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
        let momentum = function!(
            GS.emr_mom,
            usize::from(loop_edge),
            minkowski.to_symbolic([])
        );
        let gamma = |argument: &Atom| {
            function!(
                idenso::dirac::AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                argument
            )
        };
        graph.underlying[loop_edge].num = Autogen::explicit(spenso::chain!(
            &left,
            &right,
            gamma(&mu),
            gamma(&momentum),
            gamma(&mu)
        ));
        let actual = integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2)?;
        let external = function!(
            GS.emr_mom,
            usize::from(external(&graph)),
            minkowski.to_symbolic([])
        );
        // gamma_mu slash(k) gamma^mu = (2-d) slash(k), and the vector
        // bubble is q_mu B/2. The evanescent factor supplies a finite +1.
        let expected = function!(idenso::dirac::AGS.gamma, &left, &right, &external)
            * (Atom::one() - scalar_bubble_through_finite(&graph));
        assert_same_counterterm(&actual, &expected)
    }

    #[test]
    fn bare_nested_bubble_poles_match_rstar_reconstruction() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1];
                node [num=1];
                edge [particle="scalar_1" num=1];
                e [style=invis];
                e -> A:0 [id=4 sink="{ufo_order:0}"];
                C:1 -> e [id=5 source="{ufo_order:0}"];
                A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
                A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
                A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
                B -> C [id=3 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, &model
        )?;
        make_massless(&mut graph);
        let settings = settings();
        let vakint = Vakint::new()?;
        let actual = integrate_propagator(&graph, &model, &vakint, &settings, 2)?;
        let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        let expected = plan.pole_part(&graph, &model, &vakint, &settings)?;
        assert_same_counterterm(&actual, &expected)
    }

    #[test]
    fn bare_propagator_rejects_ignored_propagator_powers() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph = bubble(&model)?;
        graph.underlying[EdgeIndex::from(0)]
            .extra_data
            .vakint_edge_power = Some(2);
        let error =
            integrate_propagator(&graph, &model, &Vakint::new()?, &settings(), 2).unwrap_err();
        assert!(error.to_string().contains("unit propagator powers"));
        Ok(())
    }
}
