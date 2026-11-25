//! Test the open soft quotient before forming its expensive hard moments.

use eyre::Result;
use linnet::half_edge::involution::HedgePair;
use symbolica::atom::{Atom, AtomCore};

use super::{CountertermCall, Evaluation, Graph, Instant, IrCorrection, top_mass};

/// A primitive two-line quotient around a two-point vertex. Its nominal
/// degree is only a prefilter: q² + m² and q²/ε can both have a nonzero Δ.
fn primitive_pair<'a>(graph: &Graph, correction: &'a IrCorrection) -> Option<&'a CountertermCall> {
    let [forest] = correction.forest.as_slice() else {
        return None;
    };
    let [call] = forest.steps.as_slice() else {
        return None;
    };
    let [slot] = call.argument.insertions.as_slice() else {
        return None;
    };
    if correction.taylor_order != 0 || call.argument.loops != 1 || slot.support != correction.hard {
        return None;
    }
    let mut edges = graph.underlying.iter_edges_of(&call.argument.remaining);
    let [
        (HedgePair::Paired { source: a, sink: b }, _, first),
        (HedgePair::Paired { source: c, sink: d }, _, second),
    ] = [edges.next()?, edges.next()?]
    else {
        return None;
    };
    if edges.next().is_some()
        || first.data.is_dummy
        || second.data.is_dummy
        || [first.data, second.data]
            .iter()
            .any(|e| e.extra_data.vakint_edge_power.is_some_and(|p| p != 1))
    {
        return None;
    }
    let [a, b, c, d] = [a, b, c, d].map(|h| graph.underlying.node_id(h));
    if a == b || c == d {
        return None;
    }
    // Both external ends lie on the connected hard component. The middle
    // vertex has exactly two half-edges in the original graph, so no real
    // external momentum or third carrier enters there, in any loop basis.
    let on_hard = |v| {
        graph
            .underlying
            .iter_nodes_of(&correction.hard)
            .any(|(n, _, _)| n == v)
    };
    let soft_end = |u, v| match (on_hard(u), on_hard(v)) {
        (true, false) => Some(v),
        (false, true) => Some(u),
        _ => None,
    };
    let middle = soft_end(a, b)?;
    if soft_end(c, d)? != middle {
        return None;
    }
    let (_, half_edges, vertex) = graph
        .underlying
        .iter_nodes()
        .find(|(n, _, _)| *n == middle)?;
    (half_edges.count() == 2 && vertex.dod.value > 0).then_some(call)
}

/// Only a regular hard attachment may multiply the completed zero kernel.
/// Scalar Laurent coefficients remain outside K. A singular dimension in a
/// tensor slot or an opaque tensor operation can defeat that separation.
fn regular_hard_attachments(graph: &Graph, correction: &IrCorrection) -> bool {
    use crate::{
        utils::GS,
        uv::rstar::{
            moments::analytic_line_edge,
            regularity::{Factor, has_tensor_structure, normalise},
        },
    };
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{atom::AtomView, function};
    fn regular(view: AtomView<'_>) -> bool {
        // Pure scalar spectators, including arbitrary Laurent coefficients,
        // belong to the hard coefficient and do not enter the quotient K.
        if !has_tensor_structure(&view.to_owned()) {
            return true;
        }
        match view {
            AtomView::Fun(_) => normalise(Factor::HardNumerator, &view.to_owned()).is_ok(),
            AtomView::Add(sum) => sum.iter().all(regular),
            AtomView::Mul(product) => product.iter().all(regular),
            AtomView::Pow(power) => {
                i64::try_from(power.get_exp()).is_ok_and(|n| n >= 0) && regular(power.get_base())
            }
            _ => true,
        }
    }
    let bare = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    let vectors: Vec<_> = graph
        .underlying
        .iter_edges_of(&correction.hard)
        .map(|(_, edge, _)| function!(GS.emr_mom, usize::from(edge), &bare))
        .collect();
    let squares: Vec<_> = vectors
        .iter()
        .map(|q| -function!(SPENSO_TAG.dot, q, q) / Atom::var(GS.mu_r_sq))
        .collect();
    let factors = graph
        .underlying
        .iter_nodes_of(&correction.hard)
        .map(|(_, _, v)| &v.num.value)
        .chain(
            graph
                .underlying
                .iter_edges_of(&correction.hard)
                .map(|(_, _, e)| &e.data.num.value),
        );
    for factor in factors {
        // For this grammar check alone, omit the known scalar invariants of
        // hard lines. Their vector arguments must have the physical dimension
        // and no open slots; the analytic exponents must be -m epsilon.
        let attachments = factor.replace_map(|v, _, out| {
            let scalar = match v {
                AtomView::Fun(f) if f.get_symbol() == SPENSO_TAG.dot && f.get_nargs() == 2 => f
                    .iter()
                    .all(|q| vectors.iter().any(|hard| hard.as_view() == q)),
                AtomView::Pow(p) if squares.iter().any(|sq| sq.as_view() == p.get_base()) => {
                    analytic_line_edge(v).is_some()
                }
                _ => false,
            };
            if scalar {
                **out = Atom::one();
            }
        });
        if !regular(attachments.as_view()) {
            return false;
        }
    }
    true
}

impl<'a> Evaluation<'a> {
    /// Used only by pure counterterm entry points. A physical pole reconstruction
    /// also needs UV children of terms that become scaleless in the rearranged
    /// root and therefore keeps its original child plan.
    pub(super) fn new_counterterm(
        plan: &'a super::RStarPlan,
        graph: &'a Graph,
        model: &'a super::Model,
        vakint: &'a super::Vakint,
        settings: &'a super::UVgenerationSettings,
        projector: Option<&'a super::RootProjector<'a>>,
    ) -> Result<Self> {
        use super::ChildMass;
        use crate::uv::rstar::{expression::dressed_scope, soft_pair_pruning_enabled};
        use std::collections::BTreeSet;
        if !soft_pair_pruning_enabled()
            || !dressed_scope()
            || plan.root_boundary.is_some()
            || !plan.ir.iter().any(|ir| primitive_pair(graph, ir).is_some())
        {
            return Self::new(
                plan,
                graph,
                model,
                vakint,
                settings,
                projector,
                None,
                ChildMass::Reference,
            );
        }
        let mut eval = Self::new_with(
            plan,
            graph,
            model,
            vakint,
            settings,
            projector,
            None,
            ChildMass::Reference,
            false,
        )?;
        let zeros = plan
            .ir
            .iter()
            .map(|ir| {
                if eval.rstar_terms.keeps(true) {
                    eval.zero_soft_pair(ir)
                } else {
                    Ok(false)
                }
            })
            .collect::<Result<Vec<_>>>()?;
        let filtered = if zeros.contains(&true) {
            let needed: BTreeSet<_> = plan
                .ordinary_uv
                .iter()
                .filter(|term| !term.remainder.is_scaleless(graph, plan.marked_edge))
                .chain(
                    plan.ir
                        .iter()
                        .zip(&zeros)
                        .filter(|(_, z)| !**z)
                        .flat_map(|(ir, _)| &ir.preparation),
                )
                .flat_map(|term| term.counterterms.iter())
                .collect();
            (needed.len() < plan.counterterms.len()).then(|| {
                let mut subset = plan.clone();
                subset.counterterms.retain(|call| needed.contains(call));
                subset
            })
        } else {
            None
        };
        if let Some(filtered) = &filtered
            && crate::uv::rstar::plan::uv_children_log()
        {
            eprintln!(
                "IR_SOFT_PAIR_CHILDREN graph={} kept={} original={}",
                graph.name,
                filtered.counterterms.len(),
                plan.counterterms.len()
            );
        }
        let child_plan = filtered.as_ref().unwrap_or(plan);
        if !child_plan.counterterms.is_empty() {
            let start = Instant::now();
            eval.children = child_plan.evaluate_uv_children_for(
                graph,
                model,
                &Atom::var(eval.reference_mass_squared).sqrt(),
                vakint,
                settings,
                ChildMass::Reference,
            )?;
            eval.proper_uv_elapsed += start.elapsed();
        }
        eval.soft_pair_zeros = Some(zeros);
        Ok(eval)
    }
}

impl Evaluation<'_> {
    /// The test kernel is evaluated only once per correction. Canonicalizing
    /// its decorated graph for Δ-cache lookup costs more than the one-loop
    /// reduction. There are no UV or IR subcalls in the admitted quotient.
    fn open_pair_is_zero(&self, call: &CountertermCall) -> Result<bool> {
        use crate::uv::rstar::{RStarPlan, expression, uncontracted};
        let graph = self
            .uv
            .decorated_graph_with_dod(&call.argument, &[Atom::one()], Some(&[0]))?;
        let support = &call.argument.remaining;
        if !RStarPlan::admits_support_root(&graph, support) {
            return Ok(false);
        }
        let plan = RStarPlan::for_support(&graph, self.model, support, self.uv_settings)?;
        if plan.taylor_order <= 0
            || !plan.counterterms.is_empty()
            || !plan.ir.is_empty()
            || plan.ordinary_uv.len() != 1
        {
            return Ok(false);
        }
        let value = uncontracted::scoped(uncontracted::Scope::Quotient, || {
            expression::in_uv_child_scope(plan.taylor_order, || {
                plan.counterterm(&graph, self.model, self.vakint, self.uv_settings)
            })
        });
        Ok(value?.is_zero())
    }

    /// At Taylor grade zero an insertion from the hard component is a
    /// momentum-independent attachment tensor. Its scalar Laurent coefficient
    /// stays *outside* the quotient K. Evaluate the full open Δ of the soft
    /// kernel, with its actual vertex and both attachment indices. If every
    /// channel is regular, contraction with any such attachment still has no
    /// poles. Testing K includes the ε dependence of the entire vertex.
    pub(super) fn zero_soft_pair(&self, correction: &IrCorrection) -> Result<bool> {
        if !crate::uv::rstar::soft_pair_pruning_enabled()
            || self.plan.root_boundary.is_some()
            || !crate::uv::rstar::expression::dressed_scope()
            || top_mass::ir_quotient()? != top_mass::IrQuotient::Rstar
            || super::ir_hard_dryrun()
            || crate::uv::rstar::uncontracted::switch().mode
                == crate::uv::rstar::uncontracted::Mode::Check
            || self.ir_dirac == crate::uv::rstar::IrDiracMode::Check
            || crate::uv::rstar::closed::closed_mode()
                == crate::uv::rstar::closed::ClosedMode::Check
            || crate::uv::rstar::closed::child_mode() == crate::uv::rstar::closed::ChildMode::Check
            || self.plan.massive_routes == crate::uv::rstar::MassiveRoutes::Check
            || top_mass::ir_outer()? == top_mass::IrOuter::Check
        {
            return Ok(false);
        }
        let Some(call) = primitive_pair(self.graph, correction) else {
            return Ok(false);
        };
        if !regular_hard_attachments(self.graph, correction) {
            return Ok(false);
        }
        let start = Instant::now();
        let zero = match self.open_pair_is_zero(call) {
            Ok(zero) => zero,
            Err(error)
                if crate::uv::rstar::ir_dirac::is_diagnostic_failure(&error)
                    || crate::uv::rstar::depth::is_alarm(&error) =>
            {
                return Err(error);
            }
            Err(_) => false,
        };
        if crate::uv::rstar::plan::uv_children_log() {
            eprintln!(
                "IR_SOFT_PAIR graph={} edges={} zero={zero} seconds={:.3}",
                self.graph.name,
                self.edges(&call.argument.remaining),
                start.elapsed().as_secs_f64()
            );
        }
        Ok(zero)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        graph::{Autogen, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        utils::{GS, load_generic_model},
        uv::rstar::{ChildMass, RStarPlan, expression::in_dressed_scope, tests::settings},
    };
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::representation::{Minkowski, RepName},
    };
    use symbolica::{atom::AtomCore, function};
    use vakint::Vakint;

    #[test]
    fn primitive_soft_pair_tests_the_complete_laurent_weight_and_keeps_higher_taylor_orders()
    -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut graph: Graph = r#"digraph G {
            graph [projector=1]; node [num=1]; e [style=invis];
            edge [particle="scalar_1" num=1];
            A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> C [id=2 source="{ufo_order:2}" sink="{ufo_order:0}"];
            C -> B [id=3 lmb_id=1 source="{ufo_order:1}" sink="{ufo_order:2}"];
            e -> A [id=4 sink="{ufo_order:3}"]; B -> e [id=5 source="{ufo_order:3}"];
        }"#
        .into_finalized_runtime_graph(&model)?;
        let edges = graph
            .underlying
            .iter_edges()
            .map(|(_, e, _)| e)
            .collect::<Vec<_>>();
        for e in edges {
            let p = &mut graph.underlying[e].particle;
            *p = p.clone().override_mass(Some(Atom::Zero));
        }
        let c = graph
            .underlying
            .iter_nodes()
            .find(|(_, _, v)| v.name.value == "C")
            .unwrap()
            .0;
        graph.underlying[c].dod = Autogen::explicit(2);
        let q = function!(GS.emr_mom, 2, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let q2 = function!(SPENSO_TAG.dot, &q, &q);
        let eps = Atom::var(GS.dim_epsilon);
        let mass_sq = Atom::var(GS.mu_r_sq);
        let x = (-&q2 / &mass_sq).pow(-&eps);
        let settings = settings();
        let vakint = Vakint::new()?;
        for (rule, should_skip) in [
            (q2.clone(), true),
            (&q2 + &mass_sq, false),
            (&q2 / &eps, false),
            (&q2 * ((x - 1) / &eps).pow(5), true),
        ] {
            let quadratic = rule == q2;
            graph.underlying[c].num = Autogen::explicit(rule);
            let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
            let candidates = plan
                .ir
                .iter()
                .filter(|ir| primitive_pair(&graph, ir).is_some())
                .collect::<Vec<_>>();
            assert_eq!(candidates.len(), 1);
            let mut differentiated = (*candidates[0]).clone();
            differentiated.taylor_order = 1;
            assert!(primitive_pair(&graph, &differentiated).is_none());
            let mut nested = (*candidates[0]).clone();
            let extra = nested.forest[0].steps[0].clone();
            nested.forest[0].steps.push(extra);
            assert!(primitive_pair(&graph, &nested).is_none());
            let evaluator = Evaluation::new_with(
                &plan,
                &graph,
                &model,
                &vakint,
                &settings,
                None,
                None,
                ChildMass::Reference,
                false,
            )?;
            assert!(!evaluator.zero_soft_pair(candidates[0])?);
            assert_eq!(
                in_dressed_scope(|| evaluator.zero_soft_pair(candidates[0]))?,
                should_skip
            );
            if quadratic {
                let reference = plan.counterterm(&graph, &model, &vakint, &settings)?;
                let pruned =
                    in_dressed_scope(|| plan.counterterm(&graph, &model, &vakint, &settings))?;
                crate::uv::rstar::tests::assert_same_counterterm(&pruned, &reference)?;
                let mut isolated = plan.clone();
                isolated.ordinary_uv.clear();
                isolated
                    .ir
                    .retain(|ir| primitive_pair(&graph, ir).is_some());
                let mut eval = in_dressed_scope(|| {
                    Evaluation::new_counterterm(&isolated, &graph, &model, &vakint, &settings, None)
                })?;
                assert!(eval.children.is_empty());
                assert!(eval.counterterm()?.is_zero());
            }
        }
        // A hard tensor has a non-polynomial dimension in its attachment.
        // An open scalar soft kernel is zero, but the ordinary evaluator
        // rejects the malformed hard slot. Pruning must preserve that error.
        use crate::numerator::aind::Aind;
        use spenso::{
            network::library::symbolic::ETS,
            structure::{representation::Euclidean, slot::ParseableAind},
        };
        let hard_vertex = graph
            .underlying
            .iter_nodes()
            .find(|(_, _, v)| v.name.value == "A")
            .unwrap()
            .0;
        let index = |n| {
            Euclidean {}
                .new_rep(GS.dim)
                .to_symbolic([Aind::Hedge(n, 0).to_atom()])
        };
        graph.underlying[hard_vertex].num = Autogen::explicit(
            function!(ETS.metric, index(2), index(3))
                .replace(GS.dim)
                .with(Atom::one() / &eps),
        );
        graph.underlying[c].num = Autogen::explicit(q2);
        let mut plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
        plan.ir.retain(|ir| primitive_pair(&graph, ir).is_some());
        plan.ordinary_uv.clear();
        assert_eq!(plan.ir.len(), 1);
        assert!(!regular_hard_attachments(&graph, &plan.ir[0]));
        let ordinary = plan
            .counterterm(&graph, &model, &vakint, &settings)
            .unwrap_err();
        let dressed =
            in_dressed_scope(|| plan.counterterm(&graph, &model, &vakint, &settings)).unwrap_err();
        assert_eq!(
            dressed.root_cause().to_string(),
            ordinary.root_cause().to_string()
        );
        Ok(())
    }
}
