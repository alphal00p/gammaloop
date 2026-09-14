#[cfg(test)]
use crate::cff::CutCFFIndex;
use crate::{
    cff::expression::OrientationID,
    debug_tags,
    graph::Graph,
    utils::{GS, W_},
    uv::{
        Integrands, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike,
            direct_3d::{Direct3dCts, DirectResidueBranches},
            integrated::IntegratedCts,
            local_3d::Localizer,
            projected_4d::Projected4dCts,
        },
        marker::UvMarker,
    },
};
use color_eyre::Result;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{
    color::{ColorSimplifier, ColorSimplifySettings},
    shorthands::metric::MetricSimplifier,
};
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps};
use spenso::network::parsing::{AtomStructureExt, StrictTensorFilter};
use symbolica::{
    atom::{Atom, AtomCore, AtomType, AtomView, Symbol},
    function,
    id::Replacement,
    symbol,
};
use three_dimensional_reps::CffGenerationContext;
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct FinalIntegrands(Integrands);

impl FinalIntegrands {
    /// Iterate over finalized semantic expressions for diagnostics.
    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (CutCFFIndex, Atom)> {
        self.0.iter().map(|(index, atom)| (*index, atom.clone()))
    }

    pub(crate) fn map(&self, f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(self.0.map(f))
    }

    pub(crate) fn zip_add(self, other: Self) -> Result<Self> {
        Ok(Self(self.0.zip_add([other.0])?))
    }

    /// Recover shared numerator factors after all selected forests are assembled.
    /// Denominator powers and function arguments stay opaque at this boundary.
    pub(crate) fn into_integrands(self) -> Integrands {
        self.0.map(|atom| {
            // Collect only complete factors after Taylor and residue mapping.
            // Opaque powers keep distinct inverse denominators, their owners,
            // and numerator powers intact; functions keep their arguments intact.
            // Wrap original functions too so an input using the temporary head
            // is restored verbatim by the single outer unwrapping pass.
            let opaque = symbol!("gammalooprs::uv::opaque_factor");
            let mut occurrence = 0usize;
            let mut protected = atom.replace_map(|view, context, out| {
                // Sharing a summed vector across contractions can distribute a
                // vanishing contracted factor across separately evaluated terms.
                // Keep these tensor sums local to their original occurrences.
                let tensor_sum = matches!(view, AtomView::Add(_))
                    && context.parent_type == Some(AtomType::Mul)
                    && view.is_tensorial(StrictTensorFilter::ContainsReps);
                if tensor_sum || matches!(view, AtomView::Pow(_) | AtomView::Fun(_)) {
                    let mut branch_local = tensor_sum;
                    view.visitor(&mut |part| {
                        branch_local |= match part {
                            AtomView::Pow(power) => !i64::try_from(power.get_base_exp().1)
                                .is_ok_and(|exponent| exponent >= 0),
                            AtomView::Fun(fun) => [
                                OrientationID::symbol(),
                                GS.theta,
                                GS.orientation_delta,
                                Symbol::IF,
                            ]
                            .contains(&fun.get_symbol()),
                            _ => false,
                        };
                        !branch_local
                    });
                    // Identical inverses must remain inside their branch guards,
                    // including inverses nested in a function or numerator power.
                    // Selectors stay with their contributions so independent
                    // scalar contractions do not become one combined network.
                    **out = if branch_local {
                        occurrence += 1;
                        function!(opaque, view, occurrence)
                    } else {
                        function!(opaque, view)
                    };
                }
            });
            loop {
                let collected = protected.collect_factors();
                if collected == protected {
                    break;
                }
                protected = collected;
            }
            protected.replace_map(|view, _, out| {
                if let AtomView::Fun(fun) = view
                    && fun.get_symbol() == opaque
                {
                    out.set_from_view(
                        &fun.iter()
                            .next()
                            .expect("opaque factor has an expression argument"),
                    );
                }
            })
        })
    }
}

pub(crate) struct FinalIntegrandBuilder<'a> {
    localizer: Localizer<'a>,
    marker: UvMarker,
}

impl<'a> FinalIntegrandBuilder<'a> {
    pub(crate) fn new(localizer: Localizer<'a>, settings: &UVgenerationSettings) -> Self {
        Self {
            localizer,
            marker: UvMarker::new(settings),
        }
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
    )]
    pub(crate) fn build_direct<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &Direct3dCts,
        integrated: &IntegratedCts,
    ) -> Result<FinalIntegrands> {
        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        let full_graph = graph.full_filter();

        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );

        // Resolve the cograph's color algebra before residue mapping repeats its
        // momentum numerator. The final pass still closes attached UV/projector
        // color indices that remain open at this boundary.
        let resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .expect("graph numerator should be available")
            .simplify_color_with(ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            })
            * global_num;
        debug_tags!(#generation, #profile, #uv, #numerator, #dump;
            stage = "final_cograph_numerator_ready",
            graph = %graph.name,
            numerator_bytes = resnum.as_view().get_byte_size(),
            file.atom = %resnum.to_canonical_string(),
            "Cograph numerator before residue mapping"
        );
        let localized_integrated = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?;
        let localized_integrated = DirectResidueBranches::from_transient(&localized_integrated)?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let localized_local = local_terms
            .branches()?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        let final_branches = localized_integrated
            .zip_add(&localized_local)?
            .multiply_key_mapped(self.localizer.orientation, graph, &resnum)?;

        // `DirectResidueBranches` is the sparse, factorization-preserving
        // representation of sum_k sigma(k) I_k while the Taylor forest is
        // built. Materialize that opaque residue-map-key selector only at the
        // evaluator boundary, after every branch-owned numerator map has been
        // applied. An explicit sum simply replaces every sigma by one; physical
        // edge directions are separate sign metadata.
        let selected = final_branches
            .materialize(!self.localizer.orientation.explicit_orientation_sum_only)?;
        Self::simplify_final(graph, &reduced, selected)
    }

    #[debug_instrument(
        graph = %graph.log_display(),
        current = %current.log_display(),
    )]
    pub(crate) fn build_projected<S: ForestNodeLike>(
        &self,
        graph: &mut Graph,
        current: &S,
        local_terms: &Projected4dCts,
        integrated: &IntegratedCts,
    ) -> Result<FinalIntegrands> {
        if current.subgraph().is_empty() {
            return Err(eyre::eyre!(
                "the empty forest must own the complete production CFF, not projected local-4D coefficients"
            ));
        }

        let reduced = graph
            .full_filter()
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        let full_graph = graph.full_filter();

        let global_num = graph.global_atom();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            global_num = %global_num.log_display(),
            "Computed global numerator"
        );
        let resnum = graph
            .numerator(&reduced, current.subgraph())
            .get_single_atom()
            .expect("graph numerator should be available")
            .simplify_color_with(ColorSimplifySettings {
                simplify_non_color: false,
                ..Default::default()
            })
            * global_num;
        let localizer = self.localizer.with_independent_source_sum();
        debug_tags!(#generation, #profile, #uv, #numerator, #dump;
            stage = "final_cograph_numerator_ready",
            graph = %graph.name,
            numerator_bytes = resnum.as_view().get_byte_size(),
            file.atom = %resnum.to_canonical_string(),
            "Cograph numerator before residue mapping"
        );
        // Only the projected local-4D route reaches this assembly boundary.
        // Its child Taylor coefficient deliberately omits the untouched
        // cograph; choose its outer CFF per independent sector here, converting
        // only the selected raw routing proposal into physical surfaces.
        let active_sectors = local_terms.sectors();
        if active_sectors.is_empty() {
            return Err(eyre::eyre!(
                "factorized local term has no active UV sectors"
            ));
        }
        // Restore the production-tree denominators after the outer CFF
        // contracts its loop-energy dependence.  The exact DDx GL0 UV ray
        // verifies that this is the full production tree, including the
        // carrier shared with the factorized self-energy coefficient.
        let fourddenoms = GS.wrap_tree_denoms(
            graph.denominator(&graph.tree_edges.subtract(&graph.initial_state_cut), |_| -1),
        );
        let allowed_zero: Integrands = localizer
            .cutset
            .residue_selector
            .generate_allowed_keys()
            .into_iter()
            .map(|index| (index, Atom::Zero))
            .collect();
        let active_edges = graph
            .iter_edges_of(&reduced)
            .filter_map(|(pair, edge, data)| {
                (pair.is_paired()
                    && !data.data.is_dummy
                    && !graph.tree_edges.includes(&graph[&edge].1))
                .then_some(edge)
            })
            .collect::<Vec<_>>();
        let mut selector_free: Option<Vec<Integrands>> = None;
        for sector in active_sectors {
            if sector.coefficient.is_zero() {
                // A disabled integrated prefix can deliberately retain a
                // typed zero sector for later forest replay. Preserve all
                // allowed cut orders without asking the outer CFF to resolve
                // a map for an identically zero coefficient.
                selector_free.get_or_insert_with(Vec::new);
                continue;
            }
            // Choose the soft Taylor routing once with the untouched cograph
            // numerator present. The child remains factorized, and both the
            // capacity oracle and every outer map consume this same expression.
            // Independent sectors are never summed before deriving their ranks.
            let (numerator, localized) = localizer.projected_cff_from_soft_momentum_proposals(
                graph,
                current.subgraph(),
                &(&sector.coefficient * &resnum),
                active_edges.iter().copied(),
                CffGenerationContext::EmbeddedCffFactor,
            )?;
            let localized = localized
                .map(|atom| atom * &fourddenoms)
                .multiply_mapped(|orientation_id, source_edge_energy_map| {
                    localizer.map_numerator(
                        graph,
                        orientation_id,
                        source_edge_energy_map,
                        &numerator,
                    )
                })?
                .map(|atom| {
                    self.marker.prefix(
                        &full_graph,
                        current.subgraph(),
                        &(atom * &sector.frozen_factor),
                    )
                });
            // Child contours have already been summed. Only actual outer
            // branches carry maps and cut orders; consume each map once before
            // adding its value, checking the complete allowed cut-key shape.
            for (_, _, integrands) in localized.iter_orientations() {
                selector_free
                    .get_or_insert_with(Vec::new)
                    .push(integrands.clone());
            }
        }
        // The integrated addback is shared with the direct route. Its localizer
        // hosts the independent cograph sum in both routes, retaining every
        // source map even when no compatible physical-prefix host survives.
        let localized_integrated = self
            .localizer
            .localize(
                &integrated.physical_finite_counterterm_atom(),
                graph,
                current,
            )?
            .combine()?
            .multiply_mapped(|orientation_id, source_edge_energy_map| {
                self.localizer
                    .map_numerator(graph, orientation_id, source_edge_energy_map, &resnum)
            })?
            .map(|atom| self.marker.prefix(&full_graph, current.subgraph(), atom));
        // Both legitimate projected maps have now consumed every still-unmapped
        // numerator factor. Production hosts are mapping metadata only in this
        // lane: sum them explicitly without ever materializing a selector or
        // traversing another numerator map.
        for (_, _, integrands) in localized_integrated.iter_orientations() {
            selector_free
                .get_or_insert_with(Vec::new)
                .push(integrands.clone());
        }
        let selector_free = selector_free.ok_or_else(|| {
            eyre::eyre!("final 3D UV integrand contains no production energy maps")
        })?;
        Self::simplify_final(graph, &reduced, allowed_zero.zip_add(selector_free)?)
    }

    /// Normalize an already mapped and selector-assembled final integrand. This
    /// tail is deliberately blind to residue maps and cannot map a numerator.
    fn simplify_final(
        graph: &Graph,
        reduced: &SuBitGraph,
        integrands: Integrands,
    ) -> Result<FinalIntegrands> {
        let energy_replacements = graph
            .as_ref()
            .iter_edges_of(reduced)
            .filter(|(pair, _, _)| pair.is_paired())
            .map(|(_, edge_id, _)| {
                let edge_id = usize::from(edge_id) as i64;
                Replacement::new(function!(GS.energy, edge_id), function!(GS.ose, edge_id))
            })
            .collect::<Vec<_>>();
        let simplified = integrands.fallible_map(|atom| {
            let mut atom = atom
                .replace_multiple(&energy_replacements)
                .replace(function!(GS.ose, W_.mass_, W_.prop_))
                .with(W_.prop_);
            // Preserve the sum of CFF denominators after residue mapping, just
            // as the Taylor stage preserves its separate propagator topologies.
            atom = atom.replace(GS.dim).with(4).simplify_metrics();
            debug_tags!(#generation, #profile, #uv, #numerator, #dump;
                stage = "final_integrand_before_color",
                graph = %graph.name,
                numerator_bytes = atom.as_view().get_byte_size(),
                file.atom = %atom.to_canonical_string(),
                "Mapped factorized integrand before final color simplification"
            );
            atom = atom
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                )
                .expand_dots()?;

            // Exact production branches have already mapped every owned
            // numerator fragment. The former coarse export sign replacement
            // has no role here and must not remap these factors a second time.
            Ok(atom
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0))
        })?;
        Ok(FinalIntegrands(simplified))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        cff::esurface::RaisedEsurfaceGroup,
        dot,
        graph::{cuts::CutSet, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        uv::{
            ApproximationType, Spinney,
            approx::{
                OrientationProjection, Rooted, UVCtx,
                local_4d::{Full4dCts, Local4dCts, uv_limit},
                projected_4d::{Projected4dApproximation, Projected4dSector},
            },
            hedge_poset::OwnedForestNode,
        },
    };
    use linnet::half_edge::subgraph::InternalSubGraph;

    #[test]
    fn final_integrand_collects_tensor_factors_without_merging_denominators() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph factorized_final {
            edge [num=1 mass=1];
            node [num=1];
            a -> b [id=0 lmb_id=0];
            a -> b [id=1];
        })?;
        let index = CutCFFIndex::new_all_none();
        let (a, b, d1, d2) = symbol!(
            "final_factor_test::a",
            "final_factor_test::b",
            "final_factor_test::D1",
            "final_factor_test::D2"
        );
        let tensor = spenso::tensor!(final_factor_test, spenso::mink!(4, mu));
        let numerator = &tensor * (Atom::var(a) + b).pow(3);
        let first = Atom::var(a) / Atom::var(d1).pow(2);
        let second = Atom::var(b) / Atom::var(d2).pow(3);
        let distinct = &numerator * &first - &numerator * &second;
        let factored = &numerator * (&first - &second);
        let selector = OrientationID(4).atom();
        // Final expressions expose shared tensor numerators while the original
        // inverse powers remain separate. Opaque arguments and powers retain
        // their exact representation, including an input using our local head.
        // A shared selector belongs to each contribution, not its numerator.
        for (input, expected) in [
            (distinct.clone(), factored.clone()),
            (
                &numerator * &selector * &first - &numerator * &selector * &second,
                &numerator * (&selector * &first - &selector * &second),
            ),
            (
                &distinct + &numerator * (&first + &second).pow(2),
                &numerator * (&first - &second + (&first + &second).pow(2)),
            ),
            (
                function!(symbol!("gammalooprs::uv::opaque_factor"), &distinct),
                function!(symbol!("gammalooprs::uv::opaque_factor"), &distinct),
            ),
            (
                &tensor * (&first + &second).pow(-1),
                &tensor * (&first + &second).pow(-1),
            ),
            (
                &tensor * &first
                    + spenso::tensor!(other_factor_test, spenso::mink!(4, mu)) * &second,
                &tensor * &first
                    + spenso::tensor!(other_factor_test, spenso::mink!(4, mu)) * &second,
            ),
        ] {
            let output = FinalIntegrandBuilder::simplify_final(
                &graph,
                &graph.full_filter(),
                [(index, input)].into_iter().collect(),
            )?
            .into_integrands();
            assert_eq!(output, [(index, expected)].into_iter().collect());
            assert_eq!(
                FinalIntegrandBuilder::simplify_final(
                    &graph,
                    &graph.full_filter(),
                    output.clone()
                )?
                .into_integrands(),
                output
            );
        }
        let first_forest = FinalIntegrandBuilder::simplify_final(
            &graph,
            &graph.full_filter(),
            [(index, &numerator * &first)].into_iter().collect(),
        )?;
        let second_forest = FinalIntegrandBuilder::simplify_final(
            &graph,
            &graph.full_filter(),
            [(index, -&numerator * &second)].into_iter().collect(),
        )?;
        assert_eq!(
            first_forest.zip_add(second_forest)?.into_integrands(),
            [(index, factored)].into_iter().collect()
        );
        Ok(())
    }

    #[test]
    fn final_factor_collection_preserves_inactive_singular_branches() -> Result<()> {
        use crate::{
            cff::expression::OrientationID,
            integrands::process::{
                evaluators::{EvaluatorStack, GenericEvaluatorFloat},
                param_builder::{ParamBuilder, ParamValuePairs},
            },
            processes::EvaluatorSettings,
            settings::global::CompilationOptimizationLevel,
            utils::F,
        };
        use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
        use spenso::algebra::complex::Complex;

        test_initialise()?;
        let x = Atom::var(symbol!("final_factor_test::guarded_x"));
        let wrapper = symbol!("final_factor_test::guarded_wrapper");
        let argument = symbol!("final_factor_test::guarded_argument");
        let numerator = (&x + 1).pow(3);
        let production_ids = [OrientationID(4), OrientationID(9), OrientationID(17)];
        let orientations = vec![EdgeVec::from_iter([Orientation::Default]); 3];
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        builder.pairs.orientations = [GS.sign(EdgeIndex(0))].into_iter().collect();
        builder.pairs.additional_params = [x.clone()].into_iter().collect();
        builder
            .add_function(wrapper, vec![argument], Atom::var(argument))
            .map_err(|error| eyre::eyre!(error))?;
        let parameter_count = builder.pairs.update_ranges();
        builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];

        // This cut has no contribution from map 17. A shared inverse must stay
        // within the guards of maps 4 and 9 even when x vanishes on map 17.
        for (factor, factor_at_two) in [
            (x.pow(-1), 0.5),
            ((x.pow(-1) + 1).pow(2), 2.25),
            (function!(wrapper, x.pow(-1)), 0.5),
        ] {
            let source = &numerator
                * (production_ids[0].atom() * 2 * &factor + production_ids[1].atom() * 3 * &factor);
            let collected = FinalIntegrands(
                [(CutCFFIndex::new_all_none(), source.clone())]
                    .into_iter()
                    .collect(),
            )
            .into_integrands();
            let finalized = collected.iter().next().unwrap().1;
            for atom in [&source, finalized] {
                let (mut stack, _) = EvaluatorStack::new_with_timings(
                    std::slice::from_ref(atom),
                    &builder,
                    &orientations,
                    &production_ids,
                    None,
                    &EvaluatorSettings::default(),
                )?;
                for compiled in [false, true] {
                    if compiled {
                        stack
                            .single_parametric
                            .activate_symjit(CompilationOptimizationLevel::O0)?;
                    }
                    for (map_id, input_x, expected) in [
                        (4, 2.0, 54.0 * factor_at_two),
                        (9, 2.0, 81.0 * factor_at_two),
                        (17, 2.0, 0.0),
                        (17, 0.0, 0.0),
                    ] {
                        let mut values = vec![Complex::new_re(F(1.0)); parameter_count];
                        values[builder.pairs.residue_map_id.value_range.start] =
                            Complex::new_re(F(map_id as f64));
                        values[builder.pairs.additional_params.value_range.start] =
                            Complex::new_re(F(input_x));
                        assert_eq!(
                            <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                                &mut stack.single_parametric
                            )(&values),
                            Complex::new_re(F(expected))
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn final_factor_collection_preserves_closed_tensor_pairs() -> Result<()> {
        use crate::{
            cff::expression::OrientationID,
            integrands::process::{
                evaluators::{EvaluatorStack, GenericEvaluatorFloat},
                param_builder::{ParamBuilder, ParamValuePairs},
            },
            processes::EvaluatorSettings,
            settings::global::CompilationOptimizationLevel,
            utils::F,
        };
        use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
        use spenso::algebra::complex::Complex;
        use symbolica::parse_lit;

        test_initialise()?;
        let index = parse_lit!(spenso::mink(4, 1));
        let temporal = GS.energy_delta(index.as_view());
        let v1 = GS.emr_vec_index(EdgeIndex(1), index.as_view()) + GS.ose(EdgeIndex(1)) * &temporal;
        let v2_minus =
            GS.emr_vec_index(EdgeIndex(2), index.as_view()) - GS.ose(EdgeIndex(2)) * &temporal;
        let v2_plus =
            GS.emr_vec_index(EdgeIndex(2), index.as_view()) + GS.ose(EdgeIndex(2)) * &temporal;
        let d = Atom::var(symbol!("final_factor_test::closed_pair_denominator"));
        let numerator = (&d + 1).pow(3);
        let production_ids = [OrientationID(4), OrientationID(9), OrientationID(17)];
        let orientations = vec![EdgeVec::from_iter([Orientation::Default]); 3];
        let first = production_ids[0].atom() * &v1 * &v2_minus * d.pow(-1);
        let second = production_ids[1].atom() * &v1 * &v2_plus * 3 * d.pow(-2);
        let source = &numerator * &first + &numerator * &second;
        let collected = FinalIntegrands(
            [(CutCFFIndex::new_all_none(), source.clone())]
                .into_iter()
                .collect(),
        )
        .into_integrands();
        let finalized = collected.iter().next().unwrap().1;
        // Recover the shared scalar numerator while each vector pair still
        // closes within its own branch before finite component evaluation.
        assert_eq!(*finalized, &numerator * (&first + &second));

        let mut builder = ParamBuilder::new_empty();
        builder.pairs.residue_map_id = ParamValuePairs::default_from_symbol(GS.residue_map_id);
        builder.pairs.orientations = [GS.sign(EdgeIndex(0))].into_iter().collect();
        builder.pairs.additional_params = std::iter::once(d)
            .chain((1..=2).flat_map(|edge| {
                (1..=3).map(move |component| GS.emr_mom(EdgeIndex(edge), GS.cind(component)))
            }))
            .chain([GS.ose(EdgeIndex(1)), GS.ose(EdgeIndex(2))])
            .collect();
        let parameter_count = builder.pairs.update_ranges();
        builder.values = vec![vec![Complex::new_re(F(0.0)); parameter_count]];
        for atom in [&source, finalized] {
            let (mut stack, _) = EvaluatorStack::new_with_timings(
                std::slice::from_ref(atom),
                &builder,
                &orientations,
                &production_ids,
                None,
                &EvaluatorSettings::default(),
            )?;
            for compiled in [false, true] {
                if compiled {
                    stack
                        .single_parametric
                        .activate_symjit(CompilationOptimizationLevel::O0)?;
                }
                // For q1=(3,4,0), q2=(-3,-4,0), E1=5, the contractions
                // are 25-5*E2 and 25+5*E2. The scalar numerator at D=2 is 27.
                for (map_id, input_d, e2, expected) in [
                    (4, 2.0, 5.0, 0.0),
                    (4, 2.0, 4.0, 67.5),
                    (9, 2.0, 5.0, 1012.5),
                    (17, 2.0, 5.0, 0.0),
                    (17, 0.0, 5.0, 0.0),
                ] {
                    let mut values = vec![Complex::new_re(F(1.0)); parameter_count];
                    values[builder.pairs.residue_map_id.value_range.start] =
                        Complex::new_re(F(map_id as f64));
                    for (slot, value) in values[builder.pairs.additional_params.value_range.clone()]
                        .iter_mut()
                        .zip([input_d, 3.0, 4.0, 0.0, -3.0, -4.0, 0.0, 5.0, e2])
                    {
                        *slot = Complex::new_re(F(value));
                    }
                    assert_eq!(
                        <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                            &mut stack.single_parametric
                        )(&values),
                        Complex::new_re(F(expected))
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn projected_zero_sectors_preserve_cut_orders_without_energy_maps() -> Result<()> {
        test_initialise()?;
        let mut graph: Graph = dot!(digraph projected_zero {
            edge [num=1 mass=1];
            node [num=1];
            a -> b [id=0 lmb_id=0];
            a -> b [id=1];
        })?;
        let current = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )
            .expect("the bubble has a compatible UV spinney"),
            topo_order: 1,
        };
        let mut cutset = CutSet::empty(graph.n_hedges());
        // This is a key-shape diagnostic: a typed zero needs neither actual
        // threshold surfaces nor a production map to preserve two cut orders.
        cutset.residue_selector.left_th_cut = Some(RaisedEsurfaceGroup {
            esurface_ids: Vec::new(),
            max_occurence: 2,
        });
        let production = Default::default();
        let pattern = OrientationPattern::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, true),
        );
        let settings = UVgenerationSettings {
            local_uv_cts_from_expanded_4d_integrands: true,
            ..Default::default()
        };
        let builder = FinalIntegrandBuilder::new(localizer, &settings);
        let local = Projected4dCts::new(vec![Projected4dSector {
            coefficient: Atom::Zero,
            frozen_factor: Atom::var(GS.numerator_sampling_scale),
        }]);
        let finalized =
            builder.build_projected(&mut graph, &current, &local, &IntegratedCts::root())?;
        assert_eq!(
            finalized.into_integrands(),
            cutset
                .residue_selector
                .generate_allowed_keys()
                .into_iter()
                .map(|index| (index, Atom::Zero))
                .collect(),
        );
        // A disabled integrated pole prefix is a legitimate zero recursion
        // input. Its next Taylor operation prunes the zero FourDSector, so
        // exercise the actual producer before the final assembly boundary.
        let zero_prefix = Full4dCts::recursion_input(
            &Local4dCts::root(),
            &IntegratedCts::root(),
            ApproximationType::PolePart,
            false,
            current.lmb(),
        )?;
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let zero_local = uv_limit(
            &zero_prefix,
            &UVCtx::new(&graph, &settings),
            &current,
            &given,
            &current,
            &given,
        )?;
        assert!(zero_local.atom().is_zero());
        let produced = Projected4dApproximation::new(localizer, &mut graph, &settings)
            .project_local_4d(
                &zero_local,
                &mut crate::uv::approx::projected_4d::Local4dProjectionContext::default(),
            )?;
        assert_eq!(
            builder
                .build_projected(&mut graph, &current, &produced, &IntegratedCts::root())?
                .into_integrands(),
            cutset
                .residue_selector
                .generate_allowed_keys()
                .into_iter()
                .map(|index| (index, Atom::Zero))
                .collect(),
            "a pruned local zero must preserve every allowed cut order without energy maps",
        );
        let missing_maps = Projected4dApproximation::new(localizer, &mut graph, &settings)
            .project_local_4d(
                &Local4dCts::root(),
                &mut crate::uv::approx::projected_4d::Local4dProjectionContext::default(),
            )
            .expect_err("a nonzero local source still requires production maps");
        assert!(
            missing_maps
                .to_string()
                .contains("no production energy maps")
        );
        let error = builder
            .build_projected(
                &mut graph,
                &current,
                &Projected4dCts::new(Vec::new()),
                &IntegratedCts::root(),
            )
            .expect_err("absent sectors are different from a deliberate typed zero");
        assert!(error.to_string().contains("no active UV sectors"));
        Ok(())
    }
}
