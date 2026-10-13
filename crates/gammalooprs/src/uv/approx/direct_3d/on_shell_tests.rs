//! Direct OS component and forest checks against the certified 4D route.

use super::*;
use crate::{
    graph::{FeynmanGraph, Graph, cuts::CutSet},
    initialisation::test_initialise,
    numerator::aind::Aind,
    settings::global::OrientationPattern,
    uv::{
        UVgenerationSettings,
        approx::{
            Rooted,
            direct_3d::Direct3dCts,
            local_3d::Localizer,
            local_4d::{Full4dCts, Local4dCts, on_shell_tests::fixture, uv_limit},
            projected_4d::Local4dProjectionContext,
        },
        overlap_numerators::NumeratorJetCache,
    },
};
use eyre::WrapErr;
use idenso::shorthands::UndoShorthands;
use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::abstract_index::AIND_SYMBOLS;
use symbolica::domains::float::Real;

#[test]
fn direct_os_numeric_tensor_oracle_matches_identical_symbolic_gamma_input() -> Result<()> {
    use crate::utils::{ArbPrec, F, FUN_LIB, FloatLike, TENSORLIB};
    use spenso::{
        algebra::complex::Complex,
        iterators::IteratableTensor,
        network::{
            ExecutionResult, MinResultRank, Sequential, SmallestDegree,
            library::{DummyLibrary, panicing::ErroringLibrary, symbolic::ExplicitKey},
            parsing::ShadowedStructure,
            symbolica_interop::Evaluate,
        },
        structure::{
            TensorStructure,
            slot::{DummyAind, IsAbstractSlot},
        },
        tensors::data::DataTensor,
    };
    use symbolica::{
        domains::{float::Complex as SymComplex, rational::Rational},
        evaluate::{ExpressionEvaluator, FunctionMap},
    };
    type Numeric = Complex<F<ArbPrec>>;
    test_initialise()?;
    let left = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    let middle = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    let right = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    let mu = spenso::mink!(4, Atom::from(Aind::new_dummy()));
    let nu = spenso::mink!(4, Atom::from(Aind::new_dummy()));
    let spatial = idenso::gamma!(&right, &left, &mu) * function!(GS.delta_vec, GS.cind(2), &mu);
    let factored = (idenso::gamma0!(&right, &middle)
        + Atom::num(2)
            * idenso::gamma!(&right, &middle, &mu)
            * function!(GS.delta_vec, GS.cind(2), &mu))
        * (idenso::gamma0!(&middle, &left)
            + Atom::i()
                * idenso::gamma!(&middle, &left, &nu)
                * function!(GS.delta_vec, GS.cind(1), &nu));
    let weighted =
        (Atom::num(2) + Atom::num(5).sqrt()).pow(-3) * Atom::num(7).sqrt().pow(-3) * &spatial
            + (Atom::num(3) - Atom::num(11).sqrt()).pow(-2) * &factored;
    for (expression, gaussian) in [(spatial, true), (factored, true), (weighted, false)] {
        let source = expression.parse_into_net()?;
        let symbolic_components = |mut network: crate::numerator::ParsingNet| -> Result<std::collections::BTreeMap<_, _>> {
            network.execute::<Sequential, SmallestDegree, _, _, _>(&*TENSORLIB.read().unwrap(), &*FUN_LIB)?;
            let ExecutionResult::Val(tensor) = network.result_tensor(&*TENSORLIB.read().unwrap())? else {
                return Err(eyre!("the asymmetric gamma control must retain its open matrix"));
            };
            let tensor = tensor.into_owned();
            let slots = tensor.external_structure_iter().map(|slot| slot.to_atom()).collect::<Vec<_>>();
            tensor.iter_flat().map(|(index, value)| {
                let mut key = slots.iter().cloned().zip(tensor.expanded_index(index)?.indices).collect::<Vec<_>>();
                key.sort_by(|left, right| left.0.cmp(&right.0));
                Ok((key, value.to_owned()))
            }).collect()
        };
        let expected = symbolic_components(source.clone())?;
        let mut materialized = source;
        DirectForestFixture::materialize_library(&mut materialized)?;
        assert_eq!(
            symbolic_components(materialized.clone())?,
            expected,
            "library materialization changed the same factored gamma input"
        );
        let mut evaluated = materialized
            .eval_tree(&FunctionMap::default(), &[Atom::var(GS.pi)])?
            .map_coeff(&|coefficient| {
                Numeric::new(F::from(&coefficient.re), F::from(&coefficient.im))
            });
        let zero = F(ArbPrec::default());
        let pi = Numeric::new(zero.clone().pi(), zero);
        let mut numeric = evaluated.evaluate(std::slice::from_ref(&pi));
        let library =
            DummyLibrary::<DataTensor<Numeric, ShadowedStructure<Aind>>, ExplicitKey<Aind>>::new();
        numeric.execute::<Sequential, MinResultRank, _, _, _>(
            &library,
            &ErroringLibrary::<symbolica::atom::Symbol>::new(),
        )?;
        let ExecutionResult::Val(tensor) = numeric.result_tensor(&library)? else {
            return Err(eyre!(
                "the numeric gamma control must retain its open matrix"
            ));
        };
        let tensor = tensor.into_owned();
        let slots = tensor
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        let mut actual = std::collections::BTreeMap::new();
        for (index, value) in tensor.iter_flat() {
            let mut key = slots
                .iter()
                .cloned()
                .zip(tensor.expanded_index(index)?.indices)
                .collect::<Vec<_>>();
            key.sort_by(|left, right| left.0.cmp(&right.0));
            actual.insert(key, value.clone());
        }
        let mut expected_slots = expected
            .keys()
            .next()
            .expect("nonzero gamma control")
            .iter()
            .map(|(slot, _)| slot.clone())
            .collect::<Vec<_>>();
        let mut actual_slots = slots;
        expected_slots.sort();
        actual_slots.sort();
        assert_eq!(
            actual_slots, expected_slots,
            "numeric contraction changed the open slots"
        );
        let zero_atom = Atom::Zero;
        let zero = Numeric::new(F(ArbPrec::default()), F(ArbPrec::default()));
        let keys = expected
            .keys()
            .chain(actual.keys())
            .collect::<std::collections::BTreeSet<_>>();
        for key in keys {
            let rational: ExpressionEvaluator<SymComplex<Rational>> = expected
                .get(key)
                .unwrap_or(&zero_atom)
                .evaluator(&[Atom::var(GS.pi)])
                .build()
                .map_err(|error| eyre!(error))?;
            let mut evaluator = rational.map_coeff(&|coefficient| {
                Numeric::new(F::from(&coefficient.re), F::from(&coefficient.im))
            });
            let actual = actual.get(key).unwrap_or(&zero);
            let expected = evaluator.evaluate_single(std::slice::from_ref(&pi));
            if gaussian {
                assert_eq!(
                    actual, &expected,
                    "numeric contraction changed gamma component {key:?}"
                );
            } else {
                let scale = actual.clone().norm().re.max(expected.clone().norm().re);
                assert!(
                    (actual.clone() - &expected).norm().re
                        <= F(ArbPrec::default().epsilon()).sqrt() * scale,
                    "numeric radical coefficient changed gamma component {key:?}: {actual} versus {expected}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn direct_os_cached_tensor_oracle_preserves_kinematic_parameter_order() -> Result<()> {
    use crate::utils::FloatLike;
    use crate::utils::{ArbPrec, F};
    use crate::uv::approx::local_4d::on_shell_forest_tests::nested_gluon_fixture;
    test_initialise()?;
    let (graph, _, _, _) = nested_gluon_fixture()?;
    let weights = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    let mut expression = Atom::num(31) * Atom::var(GS.pi);
    for (component, weight) in weights[..4].iter().enumerate() {
        // Both physical external-edge names carry the same live momentum.
        let edge = usize::from(component == 3);
        expression += Atom::num(*weight) * function!(GS.emr_mom, edge, GS.cind(component));
    }
    for (coordinate, edge) in graph.loop_momentum_basis.loop_edges.iter().enumerate() {
        for component in 1..=3 {
            expression += Atom::num(weights[3 + 3 * coordinate + component])
                * function!(GS.emr_mom, usize::from(*edge), GS.cind(component));
        }
    }
    let mut evaluator =
        ForestTensorEvaluator::new(DirectForestFixture::OsQuarkInGluon, &graph, 2, &expression)?;
    let denominators = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
    for (numerators, exact_rational_part) in [
        ([1, -2, 3, -5, 7, -11, 13, -17, 19, -23], -15),
        ([-3, 5, -7, 11, -13, 17, -19, 23, -29, 31], 16),
    ] {
        let coordinates = std::array::from_fn::<_, 10, _>(|index| {
            Atom::num(numerators[index]) / Atom::num(denominators[index])
        });
        let external = std::array::from_fn(|index| coordinates[index].clone());
        let loops: [[Atom; 3]; 2] = std::array::from_fn(|loop_index| {
            std::array::from_fn(|component| coordinates[4 + 3 * loop_index + component].clone())
        });
        let values = evaluator.evaluate(&external, &loops)?;
        assert_eq!(values.len(), 1, "the layout control is one scalar");
        let value = values
            .get(&Vec::new())
            .expect("the scalar has no open slot");
        let zero = F(ArbPrec::default());
        let expected = ForestNumber::new(
            F::from(&symbolica::domains::rational::Rational::from(31)) * zero.clone().pi()
                + F::from(&symbolica::domains::rational::Rational::from(
                    exact_rational_part,
                )),
            zero,
        );
        assert!(
            (value.clone() - &expected).norm().re
                <= F(ArbPrec::default().epsilon()).sqrt() * expected.norm().re,
            "cached scalar coordinates must give 31*pi+{exact_rational_part}: {value}"
        );
    }
    let unsampled = function!(
        GS.emr_mom,
        usize::from(*graph.loop_momentum_basis.loop_edges.first().unwrap()),
        GS.cind(0)
    );
    assert!(
        ForestTensorEvaluator::new(DirectForestFixture::OsQuarkInGluon, &graph, 2, &unsampled)
            .is_err(),
        "an unsampled loop energy must not become a zero-valued cache parameter"
    );
    Ok(())
}

#[test]
fn direct_os_spatial_routing_preserves_sampled_temporal_components() -> Result<()> {
    test_initialise()?;
    let (graph, current, _) = fixture(false)?;
    let replacements =
        graph.uv_spatial_wrapped_replacement(current.subgraph(), current.lmb(), &[W_.x___]);
    let temporal = function!(GS.emr_mom, 3, GS.cind(0));
    assert_eq!(temporal.replace_multiple(&replacements), temporal);
    for component in 1..=3 {
        let spatial = function!(GS.emr_mom, 3, GS.cind(component))
            .replace_multiple(&replacements)
            .replace(function!(GS.emr_mom, 0, GS.cind(component)))
            .with(0)
            .replace(function!(GS.emr_mom, 1, GS.cind(component)))
            .with(0);
        // In the rest frame the gluon and fermion have opposite spatial
        // momenta, independently of the chosen orientation of the gluon line.
        assert!(
            (spatial.pow(2) - function!(GS.emr_mom, 2, GS.cind(component)).pow(2))
                .expand_num()
                .is_zero()
        );
    }
    Ok(())
}

#[test]
fn proper_direct_os_child_matches_projected_four_d_in_both_orientation_modes() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check(DirectForestTerms::Child, false)
}

#[test]
fn proper_direct_os_child_in_quark_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInQuark.check(DirectForestTerms::Child, false)
}

#[test]
fn proper_direct_completed_two_loop_os_parent_keeps_live_cograph() -> Result<()> {
    DirectForestFixture::OsQuarkInQuarkWithSpectator.check(DirectForestTerms::ParentOfChild, false)
}

#[test]
fn proper_direct_bare_two_loop_os_parent_keeps_live_cograph() -> Result<()> {
    DirectForestFixture::OsQuarkInQuarkWithSpectator.check(DirectForestTerms::BareParent, false)
}

#[test]
fn proper_direct_two_loop_os_finite_prefix_keeps_live_cograph() -> Result<()> {
    DirectForestFixture::OsQuarkInQuarkWithSpectator.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn proper_direct_soft_child_in_os_quark_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::SoftGluonInQuark.check(DirectForestTerms::Child, false)
}

#[test]
fn proper_direct_disjoint_left_child_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarks.check(DirectForestTerms::Child, false)
}

#[test]
fn proper_direct_disjoint_right_child_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarksReversed.check(DirectForestTerms::Child, false)
}

#[test]
fn direct_os_quark_in_gluon_completed_parent_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check(DirectForestTerms::ParentOfChild, false)
}

#[test]
fn direct_os_quark_in_quark_completed_parent_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInQuark.check(DirectForestTerms::ParentOfChild, false)
}

#[test]
fn direct_soft_gluon_in_os_quark_completed_parent_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::SoftGluonInQuark.check(DirectForestTerms::ParentOfChild, false)
}

#[test]
fn direct_os_quark_finite_prefix_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInQuark.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn direct_disjoint_os_quarks_match_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarks.check(DirectForestTerms::ParentOfChild, false)
}

#[test]
fn direct_disjoint_os_quark_finite_prefix_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarks.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn direct_disjoint_os_quark_reverse_finite_prefix_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarksReversed.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn direct_os_quark_in_gluon_finite_prefix_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn direct_soft_gluon_in_os_quark_finite_prefix_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::SoftGluonInQuark.check(DirectForestTerms::ParentOfChild, true)
}

#[test]
fn direct_os_quark_in_gluon_complete_forest_matches_at_first_point() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check_points(
        DirectForestTerms::CompleteForest,
        false,
        &[[7, 1, 2, 3]],
    )
}

#[test]
fn direct_os_quark_in_gluon_complete_forest_matches_at_second_point() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check_points(
        DirectForestTerms::CompleteForest,
        false,
        &[[11, -2, 1, 4]],
    )
}

#[test]
fn direct_os_quark_in_quark_complete_local_forest_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::OsQuarkInQuark.check(DirectForestTerms::CompleteForest, false)
}

#[test]
fn direct_soft_gluon_in_os_quark_complete_local_forest_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::SoftGluonInQuark.check(DirectForestTerms::CompleteForest, false)
}

#[test]
fn direct_disjoint_os_quarks_complete_finite_product_matches_projected_four_d() -> Result<()> {
    DirectForestFixture::DisjointOsQuarks.check(DirectForestTerms::CompleteForest, true)
}

#[derive(Clone, Copy, Debug)]
enum DirectForestTerms {
    Child,
    BareParent,
    ParentOfChild,
    CompleteForest,
}

#[derive(Clone, Copy, Debug)]
pub(super) enum DirectForestFixture {
    OsQuarkInGluon,
    OsQuarkInQuark,
    OsQuarkInQuarkWithSpectator,
    SoftGluonInQuark,
    DisjointOsQuarks,
    DisjointOsQuarksReversed,
}

impl DirectForestFixture {
    fn materialize_library(network: &mut crate::numerator::ParsingNet) -> Result<()> {
        use crate::utils::TENSORLIB;
        use spenso::iterators::IteratableTensor;
        use spenso::network::{
            graph::{NetworkLeaf, NetworkNode},
            parsing::ShadowedStructure,
            store::TensorScalarStore,
        };
        let library_nodes = network
            .graph
            .graph
            .iter_nodes()
            .filter_map(|(node, _, data)| {
                matches!(data, NetworkNode::Leaf(NetworkLeaf::LibraryKey { .. })).then_some(node)
            })
            .collect::<Vec<_>>();
        for node in library_nodes {
            let tensor = network
                .graph
                .get_lib_data::<ShadowedStructure<Aind>, _, _>(&*TENSORLIB.read().unwrap(), node)?;
            let mut unresolved = std::collections::BTreeSet::new();
            for (_, entry) in tensor.iter_flat() {
                let _ = entry.replace_map(|view, _, _| {
                    if let AtomView::Fun(function) = view
                        && function.get_symbol().get_stripped_name() == "gamma"
                    {
                        unresolved.insert(view.to_owned());
                    }
                });
            }
            if !unresolved.is_empty() {
                return Err(eyre!(
                    "unresolved gamma entries from library node {:?}: {:?}",
                    network.graph.graph[node],
                    unresolved,
                ));
            }
            let index = network.store.add_tensor(tensor);
            network.graph.graph[node] = NetworkNode::Leaf(NetworkLeaf::LocalTensor(index));
        }
        Ok(())
    }

    fn assemble(self, terms: DirectForestTerms, finite_prefix: bool) -> Result<(Graph, Vec<Atom>)> {
        let complete_forest = matches!(terms, DirectForestTerms::CompleteForest);
        let bare_parent = matches!(terms, DirectForestTerms::BareParent);
        let parent_operation = !matches!(terms, DirectForestTerms::Child);
        use crate::uv::{hedge_poset::OwnedForestNode, spinney::Spinney};
        use crate::{
            cff::expression::OrientationID,
            uv::approx::{
                final_integrand::FinalIntegrandBuilder,
                integrated::{Integrated, IntegratedCts},
                local_4d::on_shell_forest_tests::{
                    FIG10_LOGARITHMIC_BOX, disjoint_quark_fixture, nested_gluon_fixture,
                    nested_quark_fixture,
                },
                projected_4d::Projected4dApproximation,
            },
        };
        test_initialise()?;
        let (mut graph, parent, child, root) = match self {
            Self::OsQuarkInGluon => nested_gluon_fixture()?,
            Self::OsQuarkInQuark => nested_quark_fixture(false)?,
            Self::OsQuarkInQuarkWithSpectator => {
                crate::uv::approx::local_4d::on_shell_forest_tests::three_loop_quark_fixture()?
            }
            Self::SoftGluonInQuark => nested_quark_fixture(true)?,
            Self::DisjointOsQuarks => {
                let (graph, left, right, root) = disjoint_quark_fixture()?;
                (graph, right, left, root)
            }
            Self::DisjointOsQuarksReversed => disjoint_quark_fixture()?,
        };
        let disjoint_graph = matches!(
            self,
            Self::DisjointOsQuarks | Self::DisjointOsQuarksReversed
        );
        let disjoint = disjoint_graph && parent_operation && !bare_parent;
        let union = if disjoint {
            Some(OwnedForestNode {
                spinney: Spinney::with_scheme(
                    InternalSubGraph::try_new(
                        parent.subgraph().union(child.subgraph()),
                        graph.as_ref(),
                    )
                    .expect("disjoint self energies form an internal subgraph"),
                    &graph,
                    &graph.loop_momentum_basis,
                    ApproximationType::MUV,
                    parent.dod() + child.dod(),
                )?
                .expect("the disjoint factors are divergent"),
                topo_order: 2,
            })
        } else {
            None
        };
        let settings = UVgenerationSettings {
            generate_integrated: finite_prefix,
            add_marker: false,
            local_uv_cts_from_expanded_4d_integrands: true,
            ..Default::default()
        };
        let full = graph.full_filter();
        let numerator = graph
            .numerator(&full, root.subgraph())
            .get_single_atom()
            .expect("the full physical source numerator is available");
        let options = graph.denominator_only_cff_3d_expression_options();
        let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        // Match production's source preparation: structural tree propagators
        // remain four-dimensional factors owned by final assembly. They are
        // not extra rank-zero energy contours, and their numerator stays in
        // this complete physical capacity analysis.
        let contract_edges =
            graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
        if disjoint_graph {
            assert_eq!(contract_edges, [EdgeIndex(4)]);
        }
        let production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &surfaces,
            &options,
            Some(&numerator),
        )?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let integrated = IntegratedCts::root();
        let child_local = uv_limit(
            &Full4dCts::from_factorized_local(&Local4dCts::root()),
            &UVCtx::new(&graph, &settings),
            &child,
            &root,
            &child,
            &root,
        )?;
        let parent_local = (complete_forest || disjoint || bare_parent)
            .then(|| {
                uv_limit(
                    &Full4dCts::from_factorized_local(&Local4dCts::root()),
                    &UVCtx::new(&graph, &settings),
                    &parent,
                    &root,
                    &parent,
                    &root,
                )
            })
            .transpose()?;
        let box_child = (complete_forest && matches!(self, Self::OsQuarkInGluon))
            .then(|| FIG10_LOGARITHMIC_BOX.node(&graph))
            .transpose()?;
        let box_local = box_child
            .as_ref()
            .map(|component| {
                uv_limit(
                    &Full4dCts::from_factorized_local(&Local4dCts::root()),
                    &UVCtx::new(&graph, &settings),
                    component,
                    &root,
                    component,
                    &root,
                )
            })
            .transpose()?;
        let integrate_child =
            |local: &Local4dCts, component: &OwnedForestNode| -> Result<IntegratedCts> {
                let vakint_settings = vakint::VakintSettings {
                    number_of_terms_in_epsilon_expansion: 4,
                    ..settings.vakint.true_settings()
                };
                let integrated = Integrated::new(crate::utils::vakint()?, &vakint_settings).run(
                    local,
                    &UVCtx::new(&graph, &settings),
                    component,
                    &root,
                    component,
                    &root,
                )?;
                assert!(!integrated.physical_finite_counterterm_atom().is_zero());
                Ok(integrated)
            };
        let child_integrated = if finite_prefix {
            integrate_child(&child_local, &child)?
        } else {
            IntegratedCts::root()
        };
        let parent_integrated = if complete_forest && disjoint && finite_prefix {
            integrate_child(
                parent_local
                    .as_ref()
                    .expect("complete forest has its second local component"),
                &parent,
            )?
        } else {
            IntegratedCts::root()
        };
        let product_integrated = if complete_forest && disjoint && finite_prefix {
            IntegratedCts::factorized_product([&child_integrated, &parent_integrated], 4)?
        } else {
            IntegratedCts::root()
        };
        let current = if let Some(union) = &union {
            union
        } else if parent_operation {
            &parent
        } else {
            &child
        };
        let parent_given = if disjoint || bare_parent {
            &root
        } else {
            &child
        };
        let mut expressions = Vec::new();
        for explicit in [false, true] {
            let orientation =
                OrientationProjection::exact_expression(&production, &options, &pattern, explicit);
            let localizer = Localizer::new(&cutset, orientation);
            let source = Direct3dCts::root(&graph, localizer)?;
            let child_direct = Direct3dApproximation::new(localizer, &mut graph, &settings)
                .run_local(&source, &child, &root, &child, &root)
                .wrap_err_with(|| format!("{self:?}: direct child, explicit={explicit}"))?;
            let direct = if bare_parent {
                Direct3dApproximation::new(localizer, &mut graph, &settings)
                    .run_local(&source, &parent, &root, &parent, &root)?
            } else if disjoint {
                let direct = Direct3dApproximation::new(localizer, &mut graph, &settings)
                    .run_local(&child_direct, &parent, &root, &parent, &root)
                    .wrap_err("disjoint direct local replay")?;
                if finite_prefix {
                    let finite = child_integrated.physical_finite_counterterm_atom();
                    let replay = Direct3dApproximation::new(localizer, &mut graph, &settings)
                        .run_integrated(
                            &[(&finite, &child)],
                            &child,
                            &parent,
                            &root,
                            &parent,
                            &root,
                        )
                        .wrap_err("disjoint direct finite replay")?;
                    assert!(
                        replay.sectors()?.iter().any(|sector| {
                            sector
                                .frozen_integrands
                                .iter()
                                .any(|(_, atom)| !atom.is_one())
                        }),
                        "the finite disjoint prefix must retain its nontrivial frozen kernel"
                    );
                    let mut sectors = direct.sectors()?.to_vec();
                    sectors.extend(replay.sectors()?.iter().cloned());
                    if complete_forest {
                        let finite = parent_integrated.physical_finite_counterterm_atom();
                        let replay = Direct3dApproximation::new(localizer, &mut graph, &settings)
                            .run_integrated(
                                &[(&finite, &parent)],
                                &parent,
                                &child,
                                &root,
                                &child,
                                &root,
                            )
                            .wrap_err("disjoint direct opposite finite replay")?;
                        assert!(
                            replay.sectors()?.iter().any(|sector| {
                                sector
                                    .frozen_integrands
                                    .iter()
                                    .any(|(_, atom)| !atom.is_one())
                            }),
                            "both finite owners retain separate frozen kernels"
                        );
                        sectors.extend(replay.sectors()?.iter().cloned());
                    }
                    Direct3dCts::from_sectors(sectors)?
                } else {
                    direct
                }
            } else if parent_operation {
                Direct3dApproximation::new(localizer, &mut graph, &settings)
                    .run(
                        &child_direct,
                        &child_integrated,
                        &parent,
                        &child,
                        &parent,
                        &child,
                    )
                    .wrap_err_with(|| {
                        format!(
                            "{self:?}: direct parent, explicit={explicit}, finite={finite_prefix}"
                        )
                    })?
            } else {
                child_direct.clone()
            };
            let builder = FinalIntegrandBuilder::new(localizer, &settings);
            let mut completed = builder.build_direct(
                &mut graph,
                current,
                &direct,
                &product_integrated,
            )
            .wrap_err_with(|| format!("{self:?}: direct final assembly, explicit={explicit}, parent={parent_operation}, finite={finite_prefix}"))?;
            if complete_forest {
                let bare = builder.build_direct(&mut graph, &root, &source, &integrated)?;
                let child =
                    builder.build_direct(&mut graph, &child, &child_direct, &child_integrated)?;
                let bare_parent = Direct3dApproximation::new(localizer, &mut graph, &settings)
                    .run_local(&source, &parent, &root, &parent, &root)?;
                let parent_term =
                    builder.build_direct(&mut graph, &parent, &bare_parent, &parent_integrated)?;
                completed = completed
                    .zip_add(bare)?
                    .zip_add(child)?
                    .zip_add(parent_term)?;
                if let Some(box_child) = &box_child {
                    let box_ct = Direct3dApproximation::new(localizer, &mut graph, &settings)
                        .run_local(&source, box_child, &root, box_child, &root)?;
                    let nested_box = Direct3dApproximation::new(localizer, &mut graph, &settings)
                        .run_local(&box_ct, &parent, box_child, &parent, box_child)?;
                    completed = completed
                        .zip_add(builder.build_direct(
                            &mut graph,
                            box_child,
                            &box_ct,
                            &integrated,
                        )?)?
                        .zip_add(builder.build_direct(
                            &mut graph,
                            &parent,
                            &nested_box,
                            &integrated,
                        )?)?;
                }
            }
            crate::debug_tags!(#generation, #profile, #uv, #summary;
                stage = "os_route_oracle_assembly_done", route = usize::from(explicit),
                "Completed direct route assembly before diagnostic materialization"
            );
            let materialization = std::time::Instant::now();
            let value = completed.iter().map(|(_, value)| value).sum::<Atom>();
            crate::debug_tags!(#generation, #profile, #uv, #summary;
                stage = "os_route_oracle_materialized", route = usize::from(explicit),
                elapsed_ms = materialization.elapsed().as_secs_f64() * 1000.0,
                "Materialized factorized diagnostic definitions"
            );
            let selection = std::time::Instant::now();
            let value = if explicit {
                value
            } else {
                (0..production.expression.orientations.len())
                    .map(|selected| OrientationID(selected).select(&value))
                    .sum()
            };
            crate::debug_tags!(#generation, #profile, #uv, #summary;
                stage = "os_route_oracle_selectors_summed", route = usize::from(explicit),
                selector_count = production.expression.orientations.len(),
                elapsed_ms = selection.elapsed().as_secs_f64() * 1000.0,
                "Summed the complete selector truth table"
            );
            expressions.push(value);
        }
        let local = if bare_parent {
            parent_local
                .as_ref()
                .expect("the bare parent owns its local operator")
                .clone()
        } else if disjoint {
            Local4dCts::from_full_product([
                Full4dCts::recursion_input(
                    &child_local,
                    &child_integrated,
                    child.renormalization_scheme(),
                    false,
                    child.lmb(),
                )?,
                Full4dCts::recursion_input(
                    parent_local
                        .as_ref()
                        .expect("complete disjoint forest owns the second component"),
                    &parent_integrated,
                    parent.renormalization_scheme(),
                    false,
                    parent.lmb(),
                )?,
            ])
        } else if parent_operation {
            let completed = Full4dCts::recursion_input(
                &child_local,
                &child_integrated,
                child.renormalization_scheme(),
                false,
                child.lmb(),
            )?;
            uv_limit(
                &completed,
                &UVCtx::new(&graph, &settings),
                &parent,
                parent_given,
                &parent,
                parent_given,
            )?
        } else {
            child_local.clone()
        };
        let orientation =
            OrientationProjection::exact_expression(&production, &options, &pattern, true);
        let localizer = Localizer::new(&cutset, orientation);
        let mut projection_context = Local4dProjectionContext::default();
        let projected = Projected4dApproximation::new(localizer, &mut graph, &settings)
            .project_local_4d(&local, &mut projection_context)
            .wrap_err_with(|| format!("{self:?}: projected coefficient, parent={parent_operation}, finite={finite_prefix}"))?;
        let builder = FinalIntegrandBuilder::new(localizer, &settings);
        let mut projected_final = builder.build_projected(
            &mut graph,
            current,
            &projected,
            &product_integrated,
        )
        .wrap_err_with(|| {
            let mut retained_soft = std::collections::BTreeSet::new();
            for (sector, coefficient) in projected.sectors().iter().enumerate() {
                for entry in coefficient.coefficient.numerators() {
                    let _ = entry.rhs.replace_map(|view, _, _| {
                        if let AtomView::Fun(momentum) = view
                            && momentum.get_symbol() == GS.emr_mom
                            && momentum.get_nargs() >= 1
                            && let Some((owner, role, payload)) =
                                GS.uv_momentum_provenance_data(momentum.get(0))
                            && role == crate::utils::symbols::UvMomentumProvenanceRole::DenominatorDerivedSoft
                        {
                            retained_soft.insert((sector, owner, payload));
                        }
                    });
                }
            }
            format!("{self:?}: projected final assembly, parent={parent_operation}, finite={finite_prefix}; retained soft carriers={retained_soft:?}")
        })?;
        if complete_forest {
            // The empty forest always owns the stored production source.
            // Projected local coefficients begin only at nonempty UV nodes.
            let source = Direct3dCts::root(&graph, localizer)?;
            let bare = builder.build_direct(&mut graph, &root, &source, &integrated)?;
            let child_projected = Projected4dApproximation::new(localizer, &mut graph, &settings)
                .project_local_4d(&child_local, &mut projection_context)?;
            let child_term =
                builder.build_projected(&mut graph, &child, &child_projected, &child_integrated)?;
            let parent_projected = Projected4dApproximation::new(localizer, &mut graph, &settings)
                .project_local_4d(
                    parent_local
                        .as_ref()
                        .expect("complete forest owns the bare parent"),
                    &mut projection_context,
                )?;
            let parent_term = builder.build_projected(
                &mut graph,
                &parent,
                &parent_projected,
                &parent_integrated,
            )?;
            projected_final = projected_final
                .zip_add(bare)?
                .zip_add(child_term)?
                .zip_add(parent_term)?;
            if let (Some(box_child), Some(box_local)) = (&box_child, &box_local) {
                let nested_box = uv_limit(
                    &Full4dCts::from_factorized_local(box_local),
                    &UVCtx::new(&graph, &settings),
                    &parent,
                    box_child,
                    &parent,
                    box_child,
                )?;
                let box_projected = Projected4dApproximation::new(localizer, &mut graph, &settings)
                    .project_local_4d(box_local, &mut projection_context)?;
                let nested_box_projected =
                    Projected4dApproximation::new(localizer, &mut graph, &settings)
                        .project_local_4d(&nested_box, &mut projection_context)?;
                projected_final = projected_final
                    .zip_add(builder.build_projected(
                        &mut graph,
                        box_child,
                        &box_projected,
                        &integrated,
                    )?)?
                    .zip_add(builder.build_projected(
                        &mut graph,
                        &parent,
                        &nested_box_projected,
                        &integrated,
                    )?)?;
            }
        }
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_assembly_done", route = 2,
            "Completed projected route assembly before diagnostic materialization"
        );
        expressions.push(projected_final.iter().map(|(_, value)| value).sum());

        Ok((graph, expressions))
    }

    fn check(self, terms: DirectForestTerms, finite_prefix: bool) -> Result<()> {
        self.check_points(terms, finite_prefix, &[[7, 1, 2, 3], [11, -2, 1, 4]])
    }

    fn check_points(
        self,
        terms: DirectForestTerms,
        finite_prefix: bool,
        points: &[[i32; 4]],
    ) -> Result<()> {
        use crate::utils::{ArbPrec, F, FloatLike};
        use spenso::algebra::complex::Complex;
        let (graph, expressions) = self.assemble(terms, finite_prefix)?;
        let mut evaluators = expressions
            .iter()
            .enumerate()
            .map(|(route, expression)| ForestTensorEvaluator::new(self, &graph, route, expression))
            .collect::<Result<Vec<_>>>()?;
        let loop_points = [[5, 0, 0], [0, 9, 0], [0, 0, 16]];
        let loops = loop_points[..graph.loop_momentum_basis.loop_edges.len()]
            .iter()
            .map(|momentum| momentum.map(Atom::num))
            .collect::<Vec<_>>();
        for &external in points {
            let external_atoms = external.map(Atom::num);
            let values = evaluators
                .iter_mut()
                .map(|evaluator| evaluator.evaluate(&external_atoms, &loops))
                .collect::<Result<Vec<_>>>()?;
            let reference = &values[2];
            let zero = Complex::new(F(ArbPrec::default()), F(ArbPrec::default()));
            assert!(
                reference
                    .values()
                    .any(|value| value.clone().norm().re > F(ArbPrec::default()).zero()),
                "the {self:?} comparison must be nonzero"
            );
            let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
            // First establish selector partition of unity, then compare the
            // independent projected route on the same finite tensor coordinates.
            for (label, actual, expected) in [
                ("selector / explicit", &values[0], &values[1]),
                ("selector / projected", &values[0], reference),
                ("explicit / projected", &values[1], reference),
            ] {
                // Exact zero entries inherit an absolute roundoff bound from
                // the complete tensor, while nonzero entries retain the same
                // relative sensitivity. A tiny cancellation residue must not
                // become its own vanishing accuracy scale.
                let tensor_scale = actual
                    .values()
                    .chain(expected.values())
                    .map(|value| value.clone().norm().re)
                    .fold(F(ArbPrec::default()), |left, right| left.max(right));
                let components = actual
                    .keys()
                    .chain(expected.keys())
                    .collect::<std::collections::BTreeSet<_>>();
                for component in components {
                    let value = actual.get(component).unwrap_or(&zero);
                    let reference = expected.get(component).unwrap_or(&zero);
                    let scale = value.clone().norm().re.max(reference.clone().norm().re);
                    assert!(
                        (value.clone() - reference).norm().re
                            <= &tolerance * (&scale + &tensor_scale),
                        "{self:?} {label}, terms={terms:?}, finite={finite_prefix}, component={component:?}, p={external:?}: {value} versus {reference}"
                    );
                }
            }
        }
        Ok(())
    }

    pub(super) fn full_forest_evaluator(self, route: usize) -> Result<ForestTensorEvaluator> {
        let (graph, expressions) = self.assemble(DirectForestTerms::CompleteForest, false)?;
        if matches!(
            self,
            Self::OsQuarkInGluon | Self::OsQuarkInQuark | Self::SoftGluonInQuark
        ) {
            assert_eq!(
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                vec![EdgeIndex(2), EdgeIndex(3)],
                "the radial fixture rays use e2=p and e3=k"
            );
        } else if matches!(
            self,
            Self::DisjointOsQuarks | Self::DisjointOsQuarksReversed
        ) {
            assert_eq!(
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                [EdgeIndex(2), EdgeIndex(5)],
                "the disjoint radial fixture rays use e2=k_left and e5=k_right"
            );
        }
        let expression = expressions.get(route).ok_or_else(|| {
            eyre!("forest oracle route must be 0 (selector), 1 (explicit), or 2 (projected)")
        })?;
        ForestTensorEvaluator::new(self, &graph, route, expression)
    }
}

type ForestNumber = spenso::algebra::complex::Complex<crate::utils::F<crate::utils::ArbPrec>>;
type ForestComponents = std::collections::BTreeMap<Vec<(Atom, usize)>, ForestNumber>;
type ForestEvaluationNetwork = spenso::network::Network<
    spenso::network::store::NetworkStore<
        spenso::tensors::parametric::EvalTreeTensor<
            ForestNumber,
            spenso::network::parsing::ShadowedStructure<Aind>,
        >,
        symbolica::evaluate::EvalTree<ForestNumber>,
    >,
    spenso::network::library::symbolic::ExplicitKey<Aind>,
    symbolica::atom::Symbol,
    Aind,
>;

/// One factored tensor network, with its source certificates and scalar
/// kinematics prepared once. Every evaluation retains all open components.
pub(super) struct ForestTensorEvaluator {
    fixture: DirectForestFixture,
    route: usize,
    loop_count: usize,
    network: ForestEvaluationNetwork,
}

impl ForestTensorEvaluator {
    fn new(
        fixture: DirectForestFixture,
        graph: &Graph,
        route: usize,
        expression: &Atom,
    ) -> Result<Self> {
        use crate::utils::F;
        use spenso::{
            network::store::TensorScalarStoreMapping, structure::slot::IsAbstractSlot,
            tensors::data::StorageTensor,
        };
        use symbolica::evaluate::FunctionMap;
        type Numeric = ForestNumber;
        let full = graph.full_filter();
        let spatial_replacements =
            graph.uv_spatial_wrapped_replacement(&full, &graph.loop_momentum_basis, &[W_.x___]);
        let tree_replacements =
            graph.uv_wrapped_replacement(&graph.tree_edges, &graph.loop_momentum_basis, &[W_.x___]);
        let loop_count = graph.loop_momentum_basis.loop_edges.len();
        let mut parameters = vec![Atom::var(GS.pi)];
        parameters.extend((0..=3).map(|index| function!(GS.emr_mom, 0, GS.cind(index))));
        for edge in &graph.loop_momentum_basis.loop_edges {
            parameters.extend(
                (1..=3).map(|index| function!(GS.emr_mom, usize::from(*edge), GS.cind(index))),
            );
        }
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_binding_start", route,
            input_bytes = expression.as_view().get_byte_size(),
            file.expression = %expression.to_plain_string(),
            "Binding one factorized diagnostic route"
        );
        let mut value = expression
            .replace(function!(GS.tree_denom_wrapper, W_.x_))
            .with(W_.x_)
            .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
            .with(W_.prop_)
            .replace(GS.m_uv_vacuum)
            .with(20)
            .replace(GS.m_uv_expansion)
            .with(20)
            .replace(GS.numerator_sampling_scale)
            .with(1)
            .replace(GS.dim_epsilon)
            .with(0)
            .replace(GS.mu_r_sq)
            .with(400)
            .replace(function!(GS.localizing_integrand, W_.x_))
            .with(W_.x_)
            .replace(GS.renormalization_localization_scale)
            .with(20);
        // A wavefunction can occupy a spinor slot in Schoonschip
        // notation. Expose that contraction before amputation, then
        // identify its fresh index with the graph's external port.
        value = value.undo_schoonschip::<Aind>()?;
        let mut spinor_ports = std::collections::BTreeMap::new();
        let mut invalid_wavefunction = false;
        let _ = value.replace_map(|view, _, _| {
            let AtomView::Fun(wavefunction) = view else {
                return;
            };
            if ![GS.u, GS.ubar, GS.v, GS.vbar].contains(&wavefunction.get_symbol()) {
                return;
            }
            if wavefunction.get_nargs() != 2
                || !matches!(wavefunction.get(1), AtomView::Fun(index)
                        if index.get_symbol() != AIND_SYMBOLS.cind
                            && index.get_nargs() == 2)
            {
                invalid_wavefunction = true;
                return;
            }
            let Ok(edge) = usize::try_from(wavefunction.get(0)) else {
                invalid_wavefunction = true;
                return;
            };
            let pair = &graph[&EdgeIndex(edge)].1;
            let Some(port) = graph.underlying[pair.any_hedge()].edge_spin_slots().next() else {
                invalid_wavefunction = true;
                return;
            };
            spinor_ports.insert(wavefunction.get(1).to_owned(), port.to_atom());
        });
        if invalid_wavefunction {
            return Err(eyre!(
                "amputation requires one explicit graph-owned external spinor port"
            ));
        }
        value = value.replace_multiple(spinor_ports.into_iter().map(|(source, target)| {
            symbolica::id::Replacement::new(source.to_pattern(), target.to_pattern())
        }));
        value = value.replace_map(|view, _, output| {
            if let AtomView::Fun(wavefunction) = view
                && [GS.u, GS.ubar, GS.v, GS.vbar].contains(&wavefunction.get_symbol())
            {
                **output = Atom::one();
            }
        });
        for (pair, edge, _) in graph.iter_edges_of(&full) {
            if !pair.is_paired() {
                continue;
            }
            let energy_squared = (1..=3).fold(graph[edge].mass_atom().pow(2), |sum, index| {
                let component =
                    graph
                        .loop_momentum_basis
                        .loop_atom(edge, GS.emr_mom, &[GS.cind(index)], true)
                        + graph.loop_momentum_basis.ext_atom(
                            edge,
                            GS.emr_mom,
                            &[GS.cind(index)],
                            true,
                        );
                sum + component.pow(2)
            });
            value = value.replace(GS.ose(edge)).with(energy_squared.sqrt());
            let mass = graph[edge].mass_atom();
            if !mass.is_zero() {
                value = value.replace(mass).with(12);
            }
        }
        value = value.replace_multiple(&spatial_replacements);
        // Actual bridge propagators carry live external four-momenta;
        // they have no sampled temporal pole to overwrite. In particular,
        // the disjoint graph's bridge must retain p rather than either OS point.
        value = value.replace_multiple(&tree_replacements);
        let bind_coordinates = |coefficient: &Atom| {
            if matches!(coefficient.as_view(), AtomView::Num(_)) {
                return coefficient.clone();
            }
            coefficient
                .replace_multiple(&spatial_replacements)
                .replace_multiple(&tree_replacements)
                .replace_map(|view, _, output| {
                    if let AtomView::Var(variable) = view {
                        if variable.get_symbol().get_stripped_name().starts_with("GC_") {
                            **output = Atom::one();
                        }
                        return;
                    }
                    let AtomView::Fun(momentum) = view else {
                        return;
                    };
                    if ![GS.emr_mom, GS.emr_vec].contains(&momentum.get_symbol())
                        || momentum.get_nargs() != 2
                    {
                        return;
                    }
                    let Ok(edge) = usize::try_from(momentum.get(0)) else {
                        return;
                    };
                    let loop_coordinate = graph
                        .loop_momentum_basis
                        .loop_edges
                        .iter()
                        .position(|coordinate| usize::from(*coordinate) == edge);
                    let offset = match loop_coordinate {
                        Some(coordinate) => 4 + 3 * coordinate,
                        None if edge == 0 || edge == 1 => 1,
                        _ => return,
                    };
                    let index = momentum.get(1);
                    if let AtomView::Fun(component) = index
                        && component.get_symbol() == AIND_SYMBOLS.cind
                        && component.get_nargs() == 1
                        && let Ok(component @ 0..=3) = usize::try_from(component.get(0))
                    {
                        if component == 0 && momentum.get_symbol() == GS.emr_vec {
                            **output = Atom::Zero;
                        } else if component != 0 || loop_coordinate.is_none() {
                            // A surviving loop energy has not been sampled;
                            // leave it unresolved so evaluation rejects it.
                            **output = parameters[offset + component].clone();
                        }
                    }
                })
        };
        value = value
            .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
            .with_map(|matched| matched.get(W_.prop_).unwrap().to_atom().sqrt());
        let value = GS.lower_esurfaces(&value).undo_dots::<Aind>()?;
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_tensor_parse_start", route,
            input_bytes = value.as_view().get_byte_size(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            file.expression = %value.to_plain_string(),
            "Bound diagnostic coordinates before finite tensor contraction"
        );
        let mut network = value.parse_into_net()?;
        // Copy the exact library tensors with their canonical
        // index permutations, then evaluate scalar coefficients
        // before contraction. The graph and all open slots stay
        // unchanged; no symbolic radical products are expanded.
        DirectForestFixture::materialize_library(&mut network)?;
        // Bind finite scalar entries after parsing, preserving typed
        // index algebra and the production evaluation order.
        let network = network.map_ref(&bind_coordinates, |tensor| {
            tensor.map_data_ref_self(&bind_coordinates)
        });
        let evaluated = network
            .eval_tree(&FunctionMap::default(), &parameters)
            .wrap_err_with(|| {
                format!("{fixture:?}: numeric tensor coefficients for route {route}")
            })?
            .map_coeff(&|coefficient| {
                Numeric::new(F::from(&coefficient.re), F::from(&coefficient.im))
            });
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_tensor_prepared", route,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Prepared the cached factorized tensor oracle"
        );
        Ok(Self {
            fixture,
            route,
            loop_count,
            network: evaluated,
        })
    }

    pub(super) fn evaluate(
        &mut self,
        external: &[Atom; 4],
        loops: &[[Atom; 3]],
    ) -> Result<ForestComponents> {
        use crate::utils::{ArbPrec, F};
        use spenso::{
            iterators::IteratableTensor,
            network::{
                ExecutionResult, MinResultRank, Sequential,
                library::{DummyLibrary, panicing::ErroringLibrary, symbolic::ExplicitKey},
                parsing::ShadowedStructure,
                symbolica_interop::Evaluate,
            },
            structure::{TensorStructure, slot::IsAbstractSlot},
            tensors::data::DataTensor,
        };
        use symbolica::atom::Symbol;
        type Numeric = ForestNumber;
        if loops.len() != self.loop_count {
            return Err(eyre!(
                "the cached forest oracle requires {} loop momenta, received {}",
                self.loop_count,
                loops.len()
            ));
        }
        let zero = F(ArbPrec::default());
        let pi = Numeric::new(zero.clone().pi(), zero);
        let mut parameters = vec![pi.clone()];
        for atom in external.iter().chain(loops.iter().flatten()) {
            let rational = atom
                .evaluator(&[Atom::var(GS.pi)])
                .build()
                .map_err(|error| eyre!(error))?;
            let mut evaluator = rational.map_coeff(&|coefficient| {
                Numeric::new(F::from(&coefficient.re), F::from(&coefficient.im))
            });
            parameters.push(evaluator.evaluate_single(std::slice::from_ref(&pi)));
        }
        let mut network = self.network.evaluate(&parameters);
        let route = self.route;
        let library =
            DummyLibrary::<DataTensor<Numeric, ShadowedStructure<Aind>>, ExplicitKey<Aind>>::new();
        let contraction = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_tensor_contraction_start", route,
            fixture = ?self.fixture,
            "Parsed the factorized diagnostic tensor network"
        );
        // Use the production batch scheduler. The compatibility
        // SequentialRef path repeatedly extracts the entire remaining
        // network, which dominates these factored forest oracles.
        network.execute::<Sequential, MinResultRank, _, _, _>(
            &library,
            &ErroringLibrary::<Symbol>::new(),
        )?;
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_tensor_contraction_done", route,
            elapsed_ms = contraction.elapsed().as_secs_f64() * 1000.0,
            "Contracted the finite diagnostic tensor network"
        );
        let tensor = match network.result_tensor(&library)? {
            ExecutionResult::One => {
                let zero = F(ArbPrec::default());
                return Ok([(Vec::new(), Numeric::new(zero.one(), zero))]
                    .into_iter()
                    .collect());
            }
            ExecutionResult::Zero => {
                return Ok([(
                    Vec::new(),
                    Numeric::new(F(ArbPrec::default()), F(ArbPrec::default())),
                )]
                .into_iter()
                .collect());
            }
            ExecutionResult::Val(tensor) => tensor.into_owned(),
        };
        let slots = tensor
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        tensor
            .iter_flat()
            .map(|(index, value)| {
                let coordinates = tensor.expanded_index(index)?;
                let mut key = slots
                    .iter()
                    .cloned()
                    .zip(coordinates.indices)
                    .collect::<Vec<_>>();
                key.sort_by(|left, right| left.0.cmp(&right.0));
                Ok((key, value.to_owned()))
            })
            .collect()
    }
}

#[test]
fn primitive_direct_massless_os_equals_soft_refined_operation() -> Result<()> {
    test_initialise()?;
    let (mut graph, current, root) = fixture(true)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        add_marker: false,
        ..Default::default()
    };
    let numerator = graph
        .numerator(current.subgraph(), root.subgraph())
        .get_single_atom()
        .expect("the massless quark numerator is available");
    let options = graph.denominator_only_cff_3d_expression_options();
    let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production =
        graph.generate_3d_expression_for_integrand(&[], &surfaces, &options, Some(&numerator))?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let mut soft_node = crate::uv::hedge_poset::OwnedForestNode {
        spinney: current.spinney.clone(),
        topo_order: current.topo_order,
    };
    soft_node.spinney.renormalization_scheme = ApproximationType::IR;
    for explicit in [false, true] {
        let orientation =
            OrientationProjection::exact_expression(&production, &options, &pattern, explicit);
        let source = Direct3dCts::root(&graph, Localizer::new(&cutset, orientation))?.branches()?;
        let point_source =
            Direct3dApproximation::new(Localizer::new(&cutset, orientation), &mut graph, &settings)
                .on_shell_source(&current, &root, current.lmb(), &source, current.subgraph())?;
        let ctx = UVCtx::new(&graph, &settings);
        let on_shell = apply_taylor(
            &ctx,
            orientation,
            &current,
            &root,
            Some(current.subgraph().clone()),
            current.lmb(),
            &source,
            Some(&point_source),
        )?;
        let soft = apply_taylor(
            &ctx,
            orientation,
            &soft_node,
            &root,
            Some(soft_node.subgraph().clone()),
            soft_node.lmb(),
            &source,
            None,
        )?;
        let difference = on_shell.zip_add(&-soft)?.materialize(false)?.resolved()?;
        for (_, difference) in difference.iter() {
            let value = difference
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .replace(function!(*LOCAL_3D_MASS_SCOPE, W_.x_))
                .with(Atom::one())
                .replace(GS.dim)
                .with(4)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .simplify_color_with(
                    ColorSimplifySettings {
                        simplify_non_color: false,
                        ..Default::default()
                    }
                    .with_cof_dimension_invariants(),
                );
            let certificate = NumeratorJetCache::default()
                .canonicalize(&GS.lower_esurfaces(&value).undo_dots::<Aind>()?)?;
            assert!(
                certificate.is_zero,
                "massless direct OS must equal H1 for explicit={explicit}: {}",
                certificate.key
            );
        }
    }
    Ok(())
}

#[test]
fn primitive_direct_os_matches_projected_four_d_in_both_orientation_modes() -> Result<()> {
    test_initialise()?;
    let (mut graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        add_marker: false,
        ..Default::default()
    };
    let bare_numerator = graph
        .numerator(current.subgraph(), root.subgraph())
        .get_single_atom()
        .expect("the physical quark numerator is available");
    let options = graph.denominator_only_cff_3d_expression_options();
    let canonical_surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &canonical_surfaces,
        &options,
        Some(&bare_numerator),
    )?;
    let pattern = OrientationPattern::default();
    let cutset = CutSet::empty(graph.n_hedges());
    let orientation =
        OrientationProjection::exact_expression(&production, &options, &pattern, false);
    let localizer = Localizer::new(&cutset, orientation);
    let ctx = UVCtx::new(&graph, &settings);
    let local = uv_limit(
        &Full4dCts::from_factorized_local(&Local4dCts::root()),
        &ctx,
        &current,
        &root,
        &current,
        &root,
    )?;
    let mut projection_context = Local4dProjectionContext::default();
    let mut projected = Atom::Zero;
    for sector in local.projection_sectors() {
        let canonical = sector.canonical_projection(&graph)?;
        let coefficient = localizer
            .project_factorized_taylor_sector(&mut graph, &canonical, &mut projection_context)?
            .resolved()?;
        projected += coefficient
            .iter()
            .map(|(_, value)| value.clone())
            .sum::<Atom>();
    }

    let mass = graph[EdgeIndex(2)].mass_atom();
    let mut direct_values = Vec::new();
    let mut direct_host_values = Vec::new();
    for materialize_selectors in [false, true] {
        let orientation = OrientationProjection::exact_expression(
            &production,
            &options,
            &pattern,
            !materialize_selectors,
        );
        let localizer = Localizer::new(&cutset, orientation);
        let source = Direct3dCts::root(&graph, localizer)?.branches()?;
        // Build each representation using its actual orientation setting.
        // Public guards remain while this private component kernel is checked.
        // The forest owns exactly this one subtraction minus.
        let point_source = Direct3dApproximation::new(localizer, &mut graph, &settings)
            .on_shell_source(&current, &root, current.lmb(), &source, current.subgraph())?;
        let direct = -apply_taylor(
            &UVCtx::new(&graph, &settings),
            orientation,
            &current,
            &root,
            Some(current.subgraph().clone()),
            current.lmb(),
            &source,
            Some(&point_source),
        )?;
        // Preserve the complete residue maps until each production host is
        // selected. This checks the two independently constructed modes before
        // summing orientations, including both fixed-point branches per host.
        let selected = direct.materialize(true)?.resolved()?;
        let selected = selected
            .iter()
            .map(|(_, value)| value.clone())
            .sum::<Atom>();
        let hosts = direct
            .iter_keys()
            .map(|(key, _)| key.selector_host)
            .collect::<std::collections::BTreeSet<_>>();
        let host_values = hosts
            .iter()
            .map(|host| (*host, host.select(&selected)))
            .collect::<std::collections::BTreeMap<_, _>>();
        let value = if materialize_selectors {
            host_values.values().cloned().sum()
        } else {
            direct
                .materialize(false)?
                .resolved()?
                .iter()
                .map(|(_, value)| value.clone())
                .sum()
        };
        direct_values.push(value);
        direct_host_values.push(host_values);
    }

    for (spatial, external) in [([9, 12, 0], [11, 2, 3, 4]), ([9, -12, 0], [19, -2, 5, 3])] {
        let specialize = |value: &Atom| {
            // Match the final evaluator boundary on both routes: direct
            // source roots retain tree ownership wrappers, while the typed
            // projection retains symbolic colour Casimirs until finalization.
            // Use the established colour invariant rules, not a numeric CF
            // replacement inferred from the difference under examination.
            let value = GS
                .erase_uv_momentum_provenance(value)
                .replace(function!(GS.tree_denom_wrapper, W_.x_))
                .with(W_.x_)
                .simplify_color_with(
                    ColorSimplifySettings {
                        simplify_non_color: false,
                        ..Default::default()
                    }
                    .with_cof_dimension_invariants(),
                )
                .replace(function!(*LOCAL_3D_MASS_SCOPE, W_.x_))
                .with(Atom::one())
                .replace(GS.dim)
                .with(Atom::num(4))
                .replace(mass.clone())
                .with(Atom::num(8))
                .replace(GS.m_uv_vacuum)
                .with(Atom::num(20))
                .replace(GS.m_uv_expansion)
                .with(Atom::num(20));
            let value = value.replace_multiple(graph.uv_spatial_wrapped_replacement(
                current.subgraph(),
                current.lmb(),
                &[W_.x___],
            ));
            let value = value.replace_map(|view, _, output| {
                let AtomView::Fun(momentum) = view else {
                    return;
                };
                if ![GS.emr_mom, GS.emr_vec].contains(&momentum.get_symbol())
                    || momentum.get_nargs() != 2
                {
                    return;
                }
                let Ok(edge) = usize::try_from(momentum.get(0)) else {
                    return;
                };
                let values = match edge {
                    0 | 1 => external,
                    2 => [0, spatial[0], spatial[1], spatial[2]],
                    _ => return,
                };
                let index = momentum.get(1);
                if let AtomView::Fun(component) = index
                    && component.get_symbol() == AIND_SYMBOLS.cind
                    && component.get_nargs() == 1
                    && let Ok(component @ 0..=3) = usize::try_from(component.get(0))
                {
                    **output = Atom::num(values[component]);
                } else {
                    let first = usize::from(momentum.get_symbol() == GS.emr_vec);
                    **output = (first..=3)
                        .map(|component| {
                            Atom::num(values[component])
                                * function!(GS.delta_vec, GS.cind(component), index)
                        })
                        .sum();
                }
            });
            let value = value
                .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
                .with_map(|matched| matched.get(W_.prop_).unwrap().to_atom().sqrt());
            GS.lower_esurfaces(&value)
        };
        let reference = specialize(&projected);
        assert!(
            !NumeratorJetCache::default()
                .canonicalize(&reference.undo_dots::<Aind>()?)?
                .is_zero,
            "the primitive route comparison must retain a nonzero open-spinor counterterm"
        );
        assert_eq!(
            direct_host_values[0].keys().collect::<Vec<_>>(),
            direct_host_values[1].keys().collect::<Vec<_>>(),
            "the two direct orientation modes must retain the same production hosts"
        );
        for (host, explicit) in &direct_host_values[0] {
            let difference = specialize(explicit) - specialize(&direct_host_values[1][host]);
            let certificate =
                NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
            assert!(
                certificate.is_zero,
                "direct OS orientation host {host:?}, k={spatial:?}, p={external:?}: {}",
                certificate.key
            );
        }
        for (mode, direct) in direct_values.iter().enumerate() {
            let actual = specialize(direct);
            let difference = &actual - &reference;
            let certificate =
                NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
            assert!(
                certificate.is_zero,
                "direct OS selector mode {mode}, k={spatial:?}, p={external:?}: {}",
                certificate.key
            );
        }
    }
    Ok(())
}
