//! Direct OS component and forest checks against the certified 4D route.

use super::*;
use crate::{
    graph::{FeynmanGraph, Graph, cuts::CutSet, parse::NumeratorInterface},
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
use idenso::{
    CookMode, CookSettings,
    shorthands::UndoShorthands,
    tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor},
};
use linnet::half_edge::involution::EdgeIndex;
use spenso::network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing};
use spenso::structure::{IndexedStructure, abstract_index::AIND_SYMBOLS};
use symbolica::domains::float::Real;

#[test]
fn direct_os_numeric_tensor_oracle_matches_identical_symbolic_gamma_input() -> Result<()> {
    use crate::utils::{ArbPrec, F, FUN_LIB, FloatLike, TENSORLIB};
    use spenso::{
        algebra::complex::Complex,
        iterators::IteratableTensor,
        network::{
            ExecutionResult, MinIntermediateCost, MinResultRank, Sequential, SmallestDegree,
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
        let numeric = evaluated.evaluate(std::slice::from_ref(&pi));
        for intermediate_cost in [false, true] {
            for prepare_sums in [false, true] {
                let mut numeric = numeric.clone();
                if prepare_sums {
                    numeric.graph.contract_ready_sum_boundaries();
                }
                let library = DummyLibrary::<
                    DataTensor<Numeric, ShadowedStructure<Aind>>,
                    ExplicitKey<Aind>,
                >::new();
                if intermediate_cost {
                    numeric.execute::<Sequential, MinIntermediateCost, _, _, _>(
                        &library,
                        &ErroringLibrary::<symbolica::atom::Symbol>::new(),
                    )?;
                } else {
                    numeric.execute::<Sequential, MinResultRank, _, _, _>(
                        &library,
                        &ErroringLibrary::<symbolica::atom::Symbol>::new(),
                    )?;
                }
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
    // This proper component has a live surrounding quark line. Nonzero fixed
    // spatial coordinates exercise its attachment boundary, which is absent
    // from the full-domain primitive comparison. This remains pointwise route
    // coverage; the exact source and reconstruction certificates run upstream.
    for velocity in [["1", "0", "0", "0"], ["5/4", "1/4", "1/2", "1/2"]] {
        DirectForestFixture::OsQuarkInQuark.check_points(
            DirectForestTerms::Child,
            false,
            &[[7, 1, 2, 3], [11, -2, 1, 4]],
            velocity,
        )?;
    }
    Ok(())
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
        ["1", "0", "0", "0"],
    )
}

#[test]
fn direct_os_quark_in_gluon_complete_forest_matches_at_second_point() -> Result<()> {
    DirectForestFixture::OsQuarkInGluon.check_points(
        DirectForestTerms::CompleteForest,
        false,
        &[[11, -2, 1, 4]],
        ["1", "0", "0", "0"],
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

    fn assemble(
        self,
        terms: DirectForestTerms,
        finite_prefix: bool,
        velocity: [&str; 4],
    ) -> Result<(Graph, Vec<crate::uv::Integrands>)> {
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
            on_shell_subtraction_velocity: velocity.map(str::to_owned),
            generate_integrated: finite_prefix,
            add_marker: false,
            local_uv_cts_from_expanded_4d_integrands: true,
            ..Default::default()
        };
        let full = graph.full_filter();
        let numerator = graph.numerator(&full, root.subgraph())?.state.expr;
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
            Some(numerator.expression()),
        )?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let integrated = IntegratedCts::root();
        let child_local = uv_limit(
            &Full4dCts::from_factorized_local(&Local4dCts::root()),
            &UVCtx::new(&graph, &graph.model, &settings),
            &child,
            &root,
            &child,
            &root,
        )?;
        let parent_local = (complete_forest || disjoint || bare_parent)
            .then(|| {
                uv_limit(
                    &Full4dCts::from_factorized_local(&Local4dCts::root()),
                    &UVCtx::new(&graph, &graph.model, &settings),
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
                    &UVCtx::new(&graph, &graph.model, &settings),
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
                    &UVCtx::new(&graph, &graph.model, &settings),
                    component,
                    &root,
                    component,
                    &root,
                )?;
                assert!(
                    !integrated
                        .physical_finite_counterterm()?
                        .expression()
                        .is_zero()
                );
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
                    let finite = child_integrated.physical_finite_counterterm()?;
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
                        let finite = parent_integrated.physical_finite_counterterm()?;
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
            let selection = std::time::Instant::now();
            let completed = completed.into_integrands().map(|value| {
                if explicit {
                    value.clone()
                } else {
                    (0..production.expression.orientations.len())
                        .map(|selected| OrientationID(selected).select(value))
                        .sum()
                }
            });
            crate::debug_tags!(#generation, #profile, #uv, #summary;
                stage = "os_route_oracle_selectors_summed", route = usize::from(explicit),
                selector_count = production.expression.orientations.len(),
                numerator_definitions = completed.numerators().len(),
                elapsed_ms = selection.elapsed().as_secs_f64() * 1000.0,
                "Summed selector truth table while retaining shared numerator bodies"
            );
            expressions.push(completed);
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
            ])?
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
                &UVCtx::new(&graph, &graph.model, &settings),
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
                    let _ = entry.rhs.expression().replace_map(|view, _, _| {
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
                    &UVCtx::new(&graph, &graph.model, &settings),
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
        expressions.push(projected_final.into_integrands());

        Ok((graph, expressions))
    }

    fn check(self, terms: DirectForestTerms, finite_prefix: bool) -> Result<()> {
        self.check_points(
            terms,
            finite_prefix,
            &[[7, 1, 2, 3], [11, -2, 1, 4]],
            ["1", "0", "0", "0"],
        )
    }

    fn check_points(
        self,
        terms: DirectForestTerms,
        finite_prefix: bool,
        points: &[[i32; 4]],
        velocity: [&str; 4],
    ) -> Result<()> {
        use crate::utils::{ArbPrec, F, FloatLike};
        use spenso::algebra::complex::Complex;
        let (graph, expressions) = self.assemble(terms, finite_prefix, velocity)?;
        let mut evaluators = expressions
            .iter()
            .enumerate()
            .map(|(route, expression)| match (self, terms) {
                (
                    Self::OsQuarkInGluon,
                    DirectForestTerms::ParentOfChild | DirectForestTerms::CompleteForest,
                ) => {
                    // These few-point checks are cheaper with native numeric
                    // tensor contraction than symbolic family compilation. The
                    // multi-point radial scans retain shared preparation through
                    // full_forest_evaluator. Resolution preserves graph factors.
                    let resolved: Atom = expression.resolved()?.iter().map(|(_, atom)| atom).sum();
                    ForestTensorEvaluator::new(self, &graph, route, &resolved)
                }
                _ => ForestTensorEvaluator::new_shared(self, &graph, route, expression),
            })
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
                        "{self:?} {label}, terms={terms:?}, finite={finite_prefix}, velocity={velocity:?}, component={component:?}, p={external:?}: {value} versus {reference}"
                    );
                }
            }
        }
        Ok(())
    }

    pub(super) fn full_forest_evaluator(self, route: usize) -> Result<ForestTensorEvaluator> {
        if matches!(self, Self::OsQuarkInQuarkWithSpectator) {
            assert_eq!(
                route, 2,
                "complete three-loop radial checks use the validated projected route"
            );
            use crate::{
                settings::global::GenerationSettings,
                uv::{
                    ApproximationType, RenormalizationPrescriptionSettings, UVOrchestrator,
                    UltravioletGraph,
                    approx::{
                        CutStructure, ForestNodeLike, OrientationProjection,
                        local_4d::on_shell_forest_tests::three_loop_quark_fixture,
                    },
                    hedge_poset::Wood,
                },
            };
            use linnet::half_edge::{
                involution::Flow,
                subgraph::{SubSetLike, SubSetOps},
            };
            use std::collections::BTreeSet;

            // The fixture tag selects the same three-loop graph used by the proper
            // two-loop tests. Its manual two-node assembler is deliberately not
            // used here: the production wood must include the outer OS operation.
            let (mut graph, _, _, root) = three_loop_quark_fixture()?;
            assert_eq!(
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                [EdgeIndex(2), EdgeIndex(3), EdgeIndex(4)],
                "the radial chart is e2=p, e3=k, e4=l",
            );
            let model = graph.model.clone();
            let settings = GenerationSettings {
                explicit_orientation_sum_only: true,
                uv: UVgenerationSettings {
                    generate_integrated: false,
                    add_marker: false,
                    local_uv_cts_from_expanded_4d_integrands: true,
                    orchestrator: UVOrchestrator::HedgePoset,
                    renormalization_prescription: RenormalizationPrescriptionSettings {
                        massive_power_divergent: ApproximationType::OS,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                ..Default::default()
            };
            let full = graph.full_filter();
            let spinneys = graph.classified_spinneys(
                &full,
                &model,
                &settings.uv,
                &graph.loop_momentum_basis,
            )?;
            assert!(spinneys.iter().all(|spinney| spinney.n_components() <= 1));
            let components = spinneys
                .iter()
                .filter(|spinney| spinney.n_components() == 1)
                .map(|spinney| {
                    assert_eq!(spinney.dod, 1);
                    assert_eq!(spinney.renormalization_scheme, ApproximationType::OS);
                    graph
                        .paired_edges(spinney.filter())
                        .into_iter()
                        .map(usize::from)
                        .collect::<BTreeSet<_>>()
                })
                .collect::<BTreeSet<_>>();
            assert_eq!(
                components,
                [
                    vec![4, 5],
                    vec![3, 4, 5, 6, 7],
                    vec![2, 3, 4, 5, 6, 7, 8, 9],
                ]
                .into_iter()
                .map(|edges| edges.into_iter().collect::<BTreeSet<_>>())
                .collect::<BTreeSet<_>>()
            );
            // Certify the six soft planes from the stored physical routing and
            // the same fixed-point definition used by on_shell_branches. Only
            // these small linear vector payloads are inspected; no denominator
            // or numerator is modified to create a soft center.
            let velocity = settings.uv.on_shell_velocity()?;
            assert_eq!(velocity, [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
            for (gluon, carrier, external, owned) in [
                (5, 4, 3, vec![4, 5]),
                (7, 3, 2, vec![3, 4, 5, 6, 7]),
                (9, 2, 0, vec![2, 3, 4, 5, 6, 7, 8, 9]),
            ] {
                let gluon = EdgeIndex(gluon);
                assert!(graph[gluon].mass_atom(&model).is_zero());
                let physical =
                    graph
                        .loop_momentum_basis
                        .loop_atom::<Atom>(gluon, GS.emr_mom, &[], true)
                        + graph
                            .loop_momentum_basis
                            .ext_atom::<Atom>(gluon, GS.emr_mom, &[], true);
                let expected = function!(GS.emr_mom, carrier) - function!(GS.emr_mom, external);
                // Expand only this linear routing identity: an overall minus
                // on a two-term sum need not have the same literal AST.
                assert!(
                    (&physical - &expected).expand().is_zero()
                        || (&physical + &expected).expand().is_zero(),
                    "physical massless edge {gluon:?} must carry the stated plane: {physical}",
                );
                let owned = owned.into_iter().collect::<BTreeSet<_>>();
                let component = spinneys
                    .iter()
                    .find(|spinney| {
                        graph
                            .paired_edges(spinney.filter())
                            .into_iter()
                            .map(usize::from)
                            .collect::<BTreeSet<_>>()
                            == owned
                    })
                    .expect("the certified OS component is present");
                let hard = component
                    .lmb
                    .loop_atom::<Atom>(gluon, GS.emr_mom, &[], true);
                let expected_hard = function!(GS.emr_mom, carrier);
                assert!(hard == expected_hard || hard == -&expected_hard);
                let soft = component.lmb.ext_atom::<Atom>(gluon, GS.emr_mom, &[], true);
                let crown = graph.dummy_less_full_crown(component.filter());
                let boundary = crown.included_iter().collect::<Vec<_>>();
                assert_eq!(boundary.len(), 2);
                assert!(
                    boundary.iter().any(|hedge| {
                        let coordinate = function!(GS.emr_mom, usize::from(graph[hedge]));
                        soft == coordinate || soft == -coordinate
                    }),
                    "the gluon soft part must be an actual component crown momentum: {soft}"
                );
                let incoming = boundary
                    .iter()
                    .find(|hedge| graph.flow(**hedge) == Flow::Sink)
                    .expect("the component has one incoming fermion port");
                let owner = graph[incoming];
                let mass = graph[owner].mass_atom(&model);
                // Both signed OS branches set the spatial crown momentum to
                // zero. With the certified hard rows above this gives l=0,
                // k=0 and p=0, separately from the three physical planes.
                for sign in [1, -1] {
                    let point = GS.on_shell_momentum(owner, &(Atom::num(sign) * &mass), &velocity);
                    for axis in 1..4 {
                        let index = GS.cind(axis);
                        assert!(
                            GS.resolve_os_momenta(&GS.indexed_momentum(&point, &[index]))
                                .is_zero()
                        );
                    }
                }
                crate::debug_tags!(#uv, #local, #soft, #summary;
                    stage = "three_loop_soft_plane_inventory",
                    gluon = usize::from(gluon), carrier, point_owner = usize::from(owner),
                    file.physical = %physical, file.hard = %hard, file.soft = %soft,
                    "Native massless physical plane and component-owned rest center"
                );
            }
            assert_eq!(
                Wood::new(CutStructure::empty(&graph), &graph, &model, &settings)?
                    .unfold()
                    .graph
                    .n_nodes(),
                8,
                "the complete local chain has all eight ordered forest families",
            );
            let numerator = graph.numerator(&full, root.subgraph())?.state.expr;
            let options = graph.denominator_only_cff_3d_expression_options();
            let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let contract_edges =
                graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
            let production = graph.generate_3d_expression_for_integrand(
                &contract_edges,
                &surfaces,
                &options,
                Some(numerator.expression()),
            )?;
            let cuts = CutStructure::empty(&graph);
            let mut generated = UVOrchestrator::HedgePoset.parametric_integrands(
                &mut graph,
                &model,
                cuts,
                crate::utils::vakint()?,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &settings.orientation_pattern,
                    true,
                ),
                &settings,
            )?;
            assert_eq!(
                generated.len(),
                1,
                "the uncut amplitude has one complete forest sum"
            );
            let generated = generated.pop().unwrap();
            return ForestTensorEvaluator::new_shared(self, &graph, route, &generated.integrands);
        }
        let (graph, expressions) = self.assemble(
            DirectForestTerms::CompleteForest,
            false,
            ["1", "0", "0", "0"],
        )?;
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
        ForestTensorEvaluator::new_shared(self, &graph, route, expression)
    }
}

type ForestNumber = spenso::algebra::complex::Complex<crate::utils::F<crate::utils::ArbPrec>>;
type ForestComponents = std::collections::BTreeMap<Vec<(Atom, usize)>, ForestNumber>;
// Execution ports remain distinct; graph labels are used only for result keys.
type ForestOutputLabels = std::collections::BTreeMap<Atom, Atom>;
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

enum ForestEvaluationBackend {
    Reference(Box<ForestEvaluationNetwork>),
    Shared {
        evaluator: Box<crate::integrands::process::evaluators::GenericEvaluator>,
        components: Vec<(Vec<(Atom, usize)>, i32)>,
        probe_parameters: usize,
    },
}

/// One factored tensor network, with its source certificates and scalar
/// kinematics prepared once. Every evaluation retains all open components.
pub(super) struct ForestTensorEvaluator {
    fixture: DirectForestFixture,
    route: usize,
    loop_count: usize,
    prepare_sum_boundaries: bool,
    output_labels: ForestOutputLabels,
    backend: ForestEvaluationBackend,
}

impl ForestTensorEvaluator {
    /// Bind the same external/loop chart in compact roots and shared bodies.
    /// Concrete momentum aliases are registered once for components materialized
    /// later by the existing tensor-family owner.
    fn bind(
        graph: &Graph,
        expressions: &crate::uv::Integrands,
    ) -> Result<(
        crate::uv::Integrands,
        crate::integrands::process::ParamBuilder,
        ForestOutputLabels,
    )> {
        use crate::integrands::process::{ParamBuilder, param_builder::FnMapEntry};
        use spenso::{
            network::library::symbolic::ETS,
            structure::{
                abstract_index::AbstractIndex,
                representation::LibrarySlot,
                slot::{IsAbstractSlot, SlotMatcher},
            },
        };
        let full = graph.full_filter();
        let spatial_replacements =
            graph.uv_spatial_wrapped_replacement(&full, &graph.loop_momentum_basis, &[W_.x___]);
        let tree_replacements =
            graph.uv_wrapped_replacement(&graph.tree_edges, &graph.loop_momentum_basis, &[W_.x___]);
        let mut parameters = vec![Atom::var(GS.pi)];
        parameters.extend((0..=3).map(|index| function!(GS.emr_mom, 0, GS.cind(index))));
        for edge in &graph.loop_momentum_basis.loop_edges {
            parameters.extend(
                (1..=3).map(|index| function!(GS.emr_mom, usize::from(*edge), GS.cind(index))),
            );
        }
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
        let prepare = |expression: &Atom| -> Result<Atom> {
            let value = expression
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
            // notation. Expose that contraction before reserving the shared
            // output frame across all roots and retained bodies.
            Ok(value
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode {
                            inner_products: false,
                            expand_inside_chains: true,
                            expand_schoonship: true,
                        },
                        trace: false,
                        chain: false,
                    },
                    ..Default::default()
                })?
                .simple_execute::<()>()?)
        };
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let expressions = expressions
            .resolved_scalars()?
            .fallible_map(&prepare)?
            .map_numerators(|tensor| {
                let source = cooking.uncook(tensor.expression().as_view());
                let prepared = prepare(&source)?;
                Ok(tensor.with_checked_expression(cooking.try_cook(prepared.as_view())?)?)
            })?;
        // Compact CFF roots have an Atom ABI. Admit only their explicit
        // finite factors once here; opaque family calls do not supply an
        // output interface. Keep these carriers until binding completes.
        let prepared_roots = expressions
            .iter()
            .map(|(key, atom)| {
                Ok((
                    *key,
                    NumeratorInterface::admit(atom.clone(), None)?.try_into_indexed()?,
                ))
            })
            .collect::<Result<std::collections::BTreeMap<_, _>>>()?;
        // The native owner reserves written indices and established logical
        // ports, including names carried only by a ranked zero. Allocate the
        // shared output frame only after every root, RHS and projector is known.
        let dummies = SymbolicTensor::<IndexedStructure>::reserved_dummies(
            prepared_roots
                .values()
                .chain(expressions.numerators().iter().map(|entry| &entry.rhs))
                .chain(std::iter::once(&graph.global_prefactor.projector)),
        );
        let mut output_ports = std::collections::BTreeMap::<usize, (Atom, Atom)>::new();
        let mut bind = |source_kind: &'static str,
                        source: &SymbolicTensor<IndexedStructure>|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let mut value = cooking.uncook(source.expression().as_view());
            let mut graph_spinor_ports =
                std::collections::BTreeMap::<usize, LibrarySlot<Aind>>::new();
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
                        if (index.get_symbol() == AIND_SYMBOLS.cind && index.get_nargs() == 1)
                            || (index.get_symbol() != AIND_SYMBOLS.cind && index.get_nargs() == 2))
                {
                    invalid_wavefunction = true;
                    return;
                }
                let Ok(edge) = usize::try_from(wavefunction.get(0)) else {
                    invalid_wavefunction = true;
                    return;
                };
                if edge >= graph.underlying.n_edges() {
                    invalid_wavefunction = true;
                    return;
                }
                let pair = &graph[&EdgeIndex(edge)].1;
                if pair.is_paired() {
                    invalid_wavefunction = true;
                    return;
                }
                let ports = graph.underlying[pair.any_hedge()]
                    .edge_spin_slots()
                    .collect::<Vec<_>>();
                let [port] = ports.as_slice() else {
                    invalid_wavefunction = true;
                    return;
                };
                if port.rep_name() != (idenso::representations::Bispinor {}).into() {
                    invalid_wavefunction = true;
                    return;
                }
                graph_spinor_ports.insert(edge, *port);
                if matches!(wavefunction.get(1), AtomView::Fun(index)
                    if index.get_symbol() == AIND_SYMBOLS.cind)
                {
                    // Parse the fixed coordinate with the native marker owner
                    // and check its range against this graph-owned spinor port.
                    // Tensor inference alone only validates the interface.
                    let component = SlotMatcher::default().concrete_component(wavefunction.get(1));
                    let dimension = usize::try_from(port.dim());
                    if !matches!((component, dimension),
                        (Some(Ok(component)), Ok(dimension)) if component < dimension)
                    {
                        invalid_wavefunction = true;
                    }
                    return;
                }
                // The edge owns this wavefunction, independently of the
                // spelling of its contraction index in another Add scope.
                invalid_wavefunction |= !LibrarySlot::<Aind>::try_from(wavefunction.get(1))
                    .is_ok_and(|source| source.rep() == port.rep());
            });
            if invalid_wavefunction {
                return Err(eyre!(
                    "amputation requires one explicit graph-owned external spinor port ({source_kind})"
                ));
            }
            let opens_spinors = !graph_spinor_ports.is_empty();
            for (edge, port) in graph_spinor_ports {
                output_ports.entry(edge).or_insert_with(|| {
                    (
                        port.rep()
                            .slot::<AbstractIndex, _>(dummies.fresh_index())
                            .to_atom(),
                        port.to_atom(),
                    )
                });
            }
            value = value.replace_map(|view, _, output| {
                if let AtomView::Fun(wavefunction) = view
                    && [GS.u, GS.ubar, GS.v, GS.vbar].contains(&wavefunction.get_symbol())
                {
                    let edge = usize::try_from(wavefunction.get(0))
                        .expect("external wavefunction owner was validated");
                    if matches!(wavefunction.get(1), AtomView::Fun(index)
                        if index.get_symbol() == AIND_SYMBOLS.cind)
                    {
                        **output =
                            function!(GS.delta_vec, wavefunction.get(1), &output_ports[&edge].0);
                    } else {
                        // Open this occurrence locally. A slot that belongs to
                        // a wavefunction in one summand may be contracted inside
                        // a different summand; globally renaming it captures
                        // that private index and changes the matrix coefficient.
                        **output =
                            function!(ETS.metric, wavefunction.get(1), &output_ports[&edge].0);
                    }
                }
            });
            let opened = if opens_spinors {
                // Amputation deliberately changes this retained interface. The
                // native transform owner checks the new explicit ports; only
                // Bispinor metrics are contracted, without re-inferring the RHS.
                source
                    .with_transformed_expression(cooking.try_cook(value.as_view())?)?
                    .0
                    .contract(idenso::tensor::ContractSettings {
                        representations: Some(&[(idenso::representations::Bispinor {}).into()]),
                        rank_one: false,
                        collect_chains: false,
                        collect_traces: false,
                        expand: false,
                        ..Default::default()
                    })?
            } else {
                source.clone()
            };
            // Keep these output indices distinct from the family's live
            // connection ports until all retained bodies are composed.
            value = cooking.uncook(opened.expression().as_view());
            for (pair, edge, _) in graph.iter_edges_of(&full) {
                if !pair.is_paired() {
                    continue;
                }
                let energy_squared =
                    (1..=3).fold(graph[edge].mass_atom(&graph.model).pow(2), |sum, index| {
                        let component = graph.loop_momentum_basis.loop_atom(
                            edge,
                            GS.emr_mom,
                            &[GS.cind(index)],
                            true,
                        ) + graph.loop_momentum_basis.ext_atom(
                            edge,
                            GS.emr_mom,
                            &[GS.cind(index)],
                            true,
                        );
                        sum + component.pow(2)
                    });
                value = value.replace(GS.ose(edge)).with(energy_squared.sqrt());
                let mass = graph[edge].mass_atom(&graph.model);
                if !mass.is_zero() {
                    value = value.replace(mass).with(12);
                }
            }
            value = value.replace_multiple(&spatial_replacements);
            // Actual bridge propagators carry live external four-momenta;
            // they have no sampled temporal pole to overwrite. In particular,
            // the disjoint graph's bridge must retain p rather than either OS point.
            value = value.replace_multiple(&tree_replacements);
            value = value
                .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
                .with_map(|matched| matched.get(W_.prop_).unwrap().to_atom().sqrt());
            let value = GS.lower_esurfaces(&value).undo_dots::<Aind>()?;
            let value = bind_coordinates(&value);
            Ok(opened.with_checked_expression(cooking.try_cook(value.as_view())?)?)
        };
        let bound_roots = prepared_roots
            .iter()
            .map(|(key, source)| {
                let bound = bind("cut-indexed root", source)?;
                Ok((*key, cooking.uncook(bound.expression().as_view())))
            })
            .collect::<Result<Vec<_>>>()?;
        let bound_bodies =
            expressions.map_numerators(|tensor| bind("shared numerator body", tensor))?;
        // resolved_scalars above has consumed the scalar definition store.
        // Rebuild only the root Atom map and retain the mapped typed RHS store.
        let expressions: crate::uv::Integrands = bound_roots.into_iter().collect();
        let expressions = expressions.with_numerators(bound_bodies.numerators().iter().cloned())?;
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = parameters.iter().cloned().collect();
        builder.pairs.update_ranges();
        let mut coordinate_aliases = std::collections::BTreeMap::new();
        for edge in 0..graph.underlying.n_edges() {
            for head in [GS.emr_mom, GS.emr_vec] {
                for component in 0..4 {
                    let lhs = function!(head, edge, GS.cind(component));
                    let rhs = bind_coordinates(&lhs);
                    if lhs != rhs {
                        // Native spatial components can normalize to the same
                        // Q leaf. One canonical component owns one exact chart;
                        // conflicting assignments are never silently replaced.
                        match coordinate_aliases.entry(lhs) {
                            std::collections::btree_map::Entry::Occupied(entry) => {
                                eyre::ensure!(
                                    entry.get() == &rhs,
                                    "canonical momentum component has conflicting bindings"
                                );
                            }
                            std::collections::btree_map::Entry::Vacant(entry) => {
                                entry.insert(rhs);
                            }
                        }
                    }
                }
            }
        }
        for (lhs, rhs) in coordinate_aliases {
            builder
                .add_function_entry(FnMapEntry {
                    lhs,
                    rhs,
                    args: vec![],
                    tags: vec![],
                    inlining: symbolica::evaluate::InliningPolicy::Always,
                    is_alias: true,
                })
                .map_err(|error| eyre!(error))?;
        }
        let output_labels = output_ports.into_values().collect();
        Ok((expressions, builder, output_labels))
    }

    /// Independent resolved-network reference, retained for small controls and
    /// one explicit same-input comparison only.
    fn new(
        fixture: DirectForestFixture,
        graph: &Graph,
        route: usize,
        expression: &Atom,
    ) -> Result<Self> {
        use crate::utils::F;
        let input = [(crate::cff::CutCFFIndex::new_all_none(), expression.clone())]
            .into_iter()
            .collect();
        let (bound, builder, output_labels) = Self::bind(graph, &input).wrap_err_with(|| {
            format!("forest reference binding: fixture={fixture:?}, route={route}, source=resolved input")
        })?;
        let value: Atom = bound.resolved()?.iter().map(|(_, atom)| atom).sum();
        let mut network = value.parse_into_net()?;
        DirectForestFixture::materialize_library(&mut network)?;
        let evaluated = network
            .eval_tree(&builder.fn_map, &builder.pairs.additional_params.params)?
            .map_coeff(&|coefficient| {
                ForestNumber::new(F::from(&coefficient.re), F::from(&coefficient.im))
            });
        Ok(Self {
            fixture,
            route,
            loop_count: graph.loop_momentum_basis.loop_edges.len(),
            prepare_sum_boundaries: true,
            output_labels,
            backend: ForestEvaluationBackend::Reference(Box::new(evaluated)),
        })
    }

    fn new_shared(
        fixture: DirectForestFixture,
        graph: &Graph,
        route: usize,
        input: &crate::uv::Integrands,
    ) -> Result<Self> {
        use crate::integrands::process::evaluators::EvaluatorStack;
        use spenso::{
            iterators::IteratableTensor,
            network::{ExecutionResult, tags::SPENSO_TAG},
            structure::{
                TensorStructure,
                representation::LibraryRep,
                slot::{DualSlotTo, IsAbstractSlot, Slot},
            },
        };
        use symbolica::atom::FunctionBuilder;
        let started = std::time::Instant::now();
        let (bound, mut builder, output_labels) = Self::bind(graph, input).wrap_err_with(|| {
            format!("forest shared binding: fixture={fixture:?}, route={route}, source=retained roots/shared bodies")
        })?;
        // The original graph owns the common external tensor interface for every
        // completed forest term. This small factorized view is used only for
        // admission; evaluator roots remain the original compact family calls.
        let reduced = graph.full_filter().subtract(&graph.initial_state_cut);
        let physical_tensor = graph
            .numerator(&reduced, &graph.empty_subgraph())?
            .state
            .expr
            .multiply(graph.global_numerator()?)?;
        let physical_tag = DirectResidueBranches::numerator_scope().1;
        let physical_lhs = function!(symbol!("gammalooprs::uv::numerator_family"), &physical_tag);
        let physical: crate::uv::Integrands = [(
            crate::cff::CutCFFIndex::new_all_none(),
            physical_lhs.clone(),
        )]
        .into_iter()
        .collect();
        let physical = physical.with_numerators([std::sync::Arc::new(
            crate::integrands::process::param_builder::NumeratorDefinition {
                lhs: physical_lhs,
                rhs: physical_tensor,
                args: vec![],
                tags: vec![physical_tag],
                inlining: symbolica::evaluate::InliningPolicy::Always,
                is_alias: false,
            },
        )])?;
        let (physical, _, physical_labels) = Self::bind(graph, &physical).wrap_err_with(|| {
            format!("forest shared binding: fixture={fixture:?}, route={route}, source=physical output interface")
        })?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut slots = physical.numerators()[0]
            .rhs
            .structure()
            .canonical()
            .external_structure()
            .into_iter()
            .map(|slot| {
                let slot = cooking.uncook(slot.to_atom().as_view());
                // Match the template's graph-owned output to this input's
                // execution frame. Only these probe-port descriptors change;
                // no retained body is renamed onto a live family-link index.
                let label = physical_labels.get(&slot).unwrap_or(&slot);
                let actual = output_labels
                    .iter()
                    .find_map(|(actual, graph_label)| (graph_label == label).then_some(actual))
                    .unwrap_or(label);
                Slot::<LibraryRep, Aind>::try_from(actual.as_view())
            })
            .collect::<std::result::Result<Vec<_>, _>>()?;
        slots.sort_by_key(|slot| slot.to_atom());
        let dimensions = slots
            .iter()
            .map(|slot| usize::try_from(slot.dim()))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let count = dimensions
            .iter()
            .try_fold(1usize, |n, d| n.checked_mul(*d))
            .ok_or_else(|| eyre!("output tensor dimensions overflow"))?;
        eyre::ensure!(
            count <= 4096,
            "the finite oracle output exceeds its component budget"
        );
        let root: Atom = bound.iter().map(|(_, atom)| atom).sum();
        // One symbolic probe closes the common output interface before shared
        // family preparation. Its component values select coordinates only at
        // evaluation time, so the same family bodies serve every output.
        let mut components = Vec::with_capacity(count);
        let mut probe_parameters = Vec::with_capacity(count);
        let probe = if slots.is_empty() {
            components.push((Vec::new(), 1));
            Atom::one()
        } else {
            let scope = DirectResidueBranches::numerator_scope().0;
            let head =
                SPENSO_TAG.tensor_symbol(&format!("gammalooprs::uv::forest_output_probe_{scope}"));
            let probe = FunctionBuilder::new(head)
                .add_args(slots.iter().map(|slot| slot.dual().to_atom()))
                .finish();
            // Materialize only this small probe. Native storage order and
            // dual-color component labels determine the parameter mapping.
            let network = probe.parse_into_net()?;
            let ExecutionResult::Val(tensor) =
                network.result_tensor(&*crate::utils::TENSORLIB.read().unwrap())?
            else {
                return Err(eyre!("the output probe must retain its tensor interface"));
            };
            let tensor = tensor.into_owned();
            let probe_slots = tensor.external_structure_iter().collect::<Vec<_>>();
            for (flat, parameter) in tensor.iter_flat() {
                let coordinates = tensor.expanded_index(flat)?;
                let mut key = probe_slots
                    .iter()
                    .zip(&coordinates.indices)
                    .map(|(slot, coordinate)| {
                        let actual = slot.dual().to_atom();
                        (
                            output_labels.get(&actual).cloned().unwrap_or(actual),
                            *coordinate,
                        )
                    })
                    .collect::<Vec<_>>();
                key.sort_by(|left, right| left.0.cmp(&right.0));
                let sign = probe_slots.iter().zip(&coordinates.indices).fold(
                    1,
                    |sign, (slot, coordinate)| {
                        if slot.rep().is_neg(*coordinate) {
                            -sign
                        } else {
                            sign
                        }
                    },
                );
                probe_parameters.push(parameter.to_owned());
                components.push((key, sign));
            }
            probe
        };
        eyre::ensure!(components.len() == count, "output probe lost components");
        builder
            .pairs
            .additional_params
            .params
            .extend(probe_parameters.iter().cloned());
        builder.pairs.update_ranges();
        let roots = [root.clone() * probe];
        let settings = crate::processes::EvaluatorSettings {
            inline_numerator_functions: false,
            do_algebra: false,
            compile: false,
            store_atom: true,
            ..Default::default()
        };
        let (stack, timings) = EvaluatorStack::new_explicit_sum_with_timings(
            &roots,
            &builder,
            bound.numerators(),
            None,
            &settings,
        )?;
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_shared_evaluator_prepared", route,
            root_bytes = root.as_view().get_byte_size(), numerator_definitions = bound.numerators().len(),
            numerator_bytes = bound.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
            output_components = components.len(), elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            file.timings = ?timings,
            "Prepared shared scalar component functions from factorized numerator families"
        );
        Ok(Self {
            fixture,
            route,
            loop_count: graph.loop_momentum_basis.loop_edges.len(),
            prepare_sum_boundaries: true,
            output_labels,
            backend: ForestEvaluationBackend::Shared {
                evaluator: Box::new(stack.single_parametric),
                components,
                probe_parameters: probe_parameters.len(),
            },
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
        let network = match &mut self.backend {
            ForestEvaluationBackend::Shared {
                evaluator,
                components,
                probe_parameters,
            } => {
                use crate::integrands::process::evaluators::GenericEvaluatorFloat;
                let started = std::time::Instant::now();
                let offset = parameters.len();
                let zero = Numeric::new(F(ArbPrec::default()).zero(), F(ArbPrec::default()).zero());
                parameters.resize(offset + *probe_parameters, zero.clone());
                let mut evaluate =
                    <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(evaluator.as_mut());
                let mut values = ForestComponents::new();
                for (index, (key, sign)) in components.iter().enumerate() {
                    if *probe_parameters != 0 {
                        let one = F(ArbPrec::default()).one();
                        parameters[offset + index] = Numeric::new(
                            if *sign < 0 { -one } else { one },
                            F(ArbPrec::default()).zero(),
                        );
                    }
                    values.insert(key.clone(), evaluate(&parameters));
                    if *probe_parameters != 0 {
                        parameters[offset + index] = zero.clone();
                    }
                }
                crate::debug_tags!(#generation, #profile, #uv, #summary;
                    stage = "os_route_oracle_shared_evaluation_done", route = self.route,
                    elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                    "Evaluated every retained tensor component"
                );
                return Ok(values);
            }
            ForestEvaluationBackend::Reference(network) => network.as_mut(),
        };
        let mut network = network.evaluate(&parameters);
        let route = self.route;
        let preparation = std::time::Instant::now();
        let input_nodes = network.graph.graph.n_nodes();
        let input_hedges = network.graph.graph.n_hedges();
        let moved_leaves = if self.prepare_sum_boundaries {
            network.graph.contract_ready_sum_boundaries()
        } else {
            0
        };
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_sum_boundary_prepared", route,
            prepared = self.prepare_sum_boundaries, moved_leaves, input_nodes, input_hedges,
            nodes = network.graph.graph.n_nodes(), hedges = network.graph.graph.n_hedges(),
            elapsed_ms = preparation.elapsed().as_secs_f64() * 1000.0,
            "Prepared finite sum interfaces without opening scalar factors"
        );
        let library =
            DummyLibrary::<DataTensor<Numeric, ShadowedStructure<Aind>>, ExplicitKey<Aind>>::new();
        let contraction = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #summary;
            stage = "os_route_oracle_tensor_contraction_start", route,
            fixture = ?self.fixture,
            intermediate_cost = false,
            nodes = network.graph.graph.n_nodes(),
            hedges = network.graph.graph.n_hedges(),
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
            intermediate_cost = false,
            nodes = network.graph.graph.n_nodes(),
            hedges = network.graph.graph.n_hedges(),
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
            .map(|slot| {
                let actual = slot.to_atom();
                self.output_labels.get(&actual).cloned().unwrap_or(actual)
            })
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
fn shared_forest_output_probe_recovers_native_components_exactly() -> Result<()> {
    use crate::utils::{FUN_LIB, TENSORLIB};
    use idenso::representations::{Bispinor, ColorFundamental};
    use spenso::{
        iterators::IteratableTensor,
        network::{ExecutionResult, Sequential, SmallestDegree, tags::SPENSO_TAG},
        structure::{
            TensorStructure,
            representation::{LibraryRep, Minkowski, RepName, Representation},
            slot::{DualSlotTo, DummyAind, IsAbstractSlot},
        },
    };
    use std::collections::BTreeMap;
    use symbolica::{atom::FunctionBuilder, id::Replacement};
    test_initialise()?;
    let x = Atom::var(symbolica::symbol!("native_output_probe::x"));
    let y = Atom::var(symbolica::symbol!("native_output_probe::y"));
    let spectator = (&x + &y).pow(11);
    let cases: Vec<(&str, Vec<Representation<LibraryRep>>)> = vec![
        ("bispinor_matrix", vec![Bispinor {}.new_rep(4).to_lib(); 2]),
        (
            "minkowski_matrix",
            vec![Minkowski {}.new_rep(4).to_lib(); 2],
        ),
        (
            "color_matrix",
            vec![
                ColorFundamental {}.new_rep(3).to_lib(),
                ColorFundamental {}.new_rep(3).dual().to_lib(),
            ],
        ),
        (
            "mixed_permuted",
            vec![
                Bispinor {}.new_rep(4).to_lib(),
                ColorFundamental {}.new_rep(3).dual().to_lib(),
                Minkowski {}.new_rep(4).to_lib(),
            ],
        ),
    ];
    for (case, reps) in cases {
        let slots = reps
            .iter()
            .map(|rep| rep.slot::<Aind, _>(Aind::new_dummy()))
            .collect::<Vec<_>>();
        let source_head = SPENSO_TAG.tensor_symbol(&format!("native_output_probe::{case}_source"));
        let probe_head = SPENSO_TAG.tensor_symbol(&format!("native_output_probe::{case}_probe"));
        let source = FunctionBuilder::new(source_head)
            .add_args(slots.iter().map(|slot| slot.to_atom()))
            .finish();
        let probe = FunctionBuilder::new(probe_head)
            .add_args(slots.iter().map(|slot| slot.dual().to_atom()))
            .finish();
        let source_net = source.parse_into_net()?;
        let ExecutionResult::Val(source_tensor) = source_net
            .result_tensor(&*TENSORLIB.read().unwrap())
            .unwrap()
        else {
            panic!("source must retain its interface")
        };
        let source_tensor = source_tensor.into_owned();
        let source_slots = source_tensor.external_structure_iter().collect::<Vec<_>>();
        let mut expected = BTreeMap::<Vec<(Atom, usize)>, Atom>::new();
        let mut source_values = Vec::new();
        for (flat, label) in source_tensor.iter_flat() {
            let coordinates = source_tensor.expanded_index(flat).unwrap();
            let mut key = source_slots
                .iter()
                .zip(coordinates.indices.iter())
                .map(|(slot, index)| (slot.to_atom(), *index))
                .collect::<Vec<_>>();
            key.sort_by(|a, b| a.0.cmp(&b.0));
            // Independent exact dense data with asymmetric real/imaginary entries.
            let weight = key
                .iter()
                .enumerate()
                .map(|(axis, (_, i))| (axis + 2) * (i + 1))
                .sum::<usize>();
            let imaginary = key[0].1 as i64 - key[key.len() - 1].1 as i64;
            let value = Atom::num(weight) * &x + Atom::i() * Atom::num(imaginary) * &y;
            source_values.push(Replacement::new(label.to_pattern(), value.clone()));
            expected.insert(key, value);
        }
        let probe_net = probe.parse_into_net()?;
        let ExecutionResult::Val(probe_tensor) = probe_net
            .result_tensor(&*TENSORLIB.read().unwrap())
            .unwrap()
        else {
            panic!("probe must retain its interface")
        };
        let probe_tensor = probe_tensor.into_owned();
        let probe_slots = probe_tensor.external_structure_iter().collect::<Vec<_>>();
        let mut parameters = Vec::new();
        let mut probe_rows = Vec::new();
        for (flat, label) in probe_tensor.iter_flat() {
            let coordinates = probe_tensor.expanded_index(flat).unwrap();
            let mut key = probe_slots
                .iter()
                .zip(coordinates.indices.iter())
                .map(|(slot, index)| (slot.dual().to_atom(), *index))
                .collect::<Vec<_>>();
            key.sort_by(|a, b| a.0.cmp(&b.0));
            let sign = probe_slots.iter().zip(coordinates.indices.iter()).fold(
                1i32,
                |sign, (slot, index)| {
                    if slot.rep().is_neg(*index) {
                        -sign
                    } else {
                        sign
                    }
                },
            );
            parameters.push(label.to_owned());
            probe_rows.push((key, sign));
        }
        assert_eq!(expected.len(), parameters.len());
        let expression = &spectator * source * probe;
        let mut network = expression.parse_into_net()?;
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&*TENSORLIB.read().unwrap(), &*FUN_LIB)
            .unwrap();
        let scalar = match network.result_scalar().unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        };
        let scalar = scalar.replace_multiple(&source_values);
        let mut all_params = vec![x.clone(), y.clone()];
        all_params.extend(parameters.iter().cloned());
        let _: symbolica::evaluate::ExpressionEvaluator<
            symbolica::domains::float::Complex<symbolica::domains::rational::Rational>,
        > = scalar.evaluator(&all_params).build().unwrap();
        for (row, (key, sign)) in probe_rows.iter().enumerate() {
            let basis = parameters
                .iter()
                .enumerate()
                .map(|(index, param)| {
                    Replacement::new(
                        param.to_pattern(),
                        Atom::num(if index == row { *sign } else { 0 }),
                    )
                })
                .collect::<Vec<_>>();
            let actual = scalar.replace_multiple(&basis);
            assert_eq!(
                actual,
                &spectator * &expected[key],
                "{case} component {key:?}"
            );
            for (a, b) in [(2, 3), (-1, 2)] {
                assert_eq!(
                    actual
                        .replace(x.to_pattern())
                        .with(a)
                        .replace(y.to_pattern())
                        .with(b),
                    (&spectator * &expected[key])
                        .replace(x.to_pattern())
                        .with(a)
                        .replace(y.to_pattern())
                        .with(b)
                );
            }
        }
    }
    Ok(())
}

#[test]
fn shared_forest_oracle_keeps_output_frame_across_matrix_and_state_families() -> Result<()> {
    use crate::{
        integrands::process::param_builder::{FnMapEntry, NumeratorDefinition},
        utils::{ArbPrec, F, FUN_LIB, TENSORLIB},
    };
    use spenso::{
        iterators::IteratableTensor,
        network::{ExecutionResult, Sequential, SmallestDegree},
        structure::{
            TensorStructure,
            slot::{DummyAind, IsAbstractSlot},
        },
        tensors::data::SparseOrDense,
    };
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    let ports = [EdgeIndex(0), EdgeIndex(1)].map(|edge| {
        graph.underlying[graph[&edge].1.any_hedge()]
            .edge_spin_slots()
            .next()
            .expect("one graph-owned spinor port")
            .to_atom()
    });
    let mu = spenso::mink!(4, Atom::from(Aind::new_dummy()));
    let matrix = spenso::g!(&ports[0], &ports[1])
        + (Atom::num(2) + Atom::i()) * idenso::gamma0!(&ports[0], &ports[1])
        + [3, 5, 7]
            .into_iter()
            .enumerate()
            .map(|(index, weight)| {
                Atom::num(weight)
                    * function!(GS.delta_vec, GS.cind(index + 1), &mu)
                    * idenso::gamma!(&ports[0], &ports[1], &mu)
            })
            .sum::<Atom>();
    let states = function!(GS.u, 0, &ports[0]) * function!(GS.ubar, 1, &ports[1]);
    let p0 = GS.emr_mom(EdgeIndex(0), GS.cind(0));
    let p3 = GS.emr_mom(EdgeIndex(0), GS.cind(3));
    let power = (&p0 + &p3).pow(11);
    let spectator = &power / (&p0 - &p3).pow(2);

    // Materialize only the small native matrix, independently of either
    // retained-family layout. Its asymmetric complex entries fix all16 axes.
    let mut network = matrix.parse_into_net()?;
    network
        .execute::<Sequential, SmallestDegree, _, _, _>(&*TENSORLIB.read().unwrap(), &*FUN_LIB)?;
    let ExecutionResult::Val(tensor) = network.result_tensor(&*TENSORLIB.read().unwrap())? else {
        return Err(eyre!("the native matrix must retain its two output ports"));
    };
    // Sparse native iteration omits zero entries. This oracle checks the
    // complete 4×4 matrix, so materialize its sixteen coordinates explicitly.
    let tensor = tensor.into_owned().to_dense();
    let slots = tensor
        .external_structure_iter()
        .map(|slot| slot.to_atom())
        .collect::<Vec<_>>();
    let mut expected_components = std::collections::BTreeMap::new();
    let mut explicit_closed = Atom::Zero;
    for (flat, value) in tensor.iter_flat() {
        let coordinates = tensor.expanded_index(flat)?;
        let mut key = slots
            .iter()
            .cloned()
            .zip(coordinates.indices)
            .collect::<Vec<_>>();
        key.sort_by(|left, right| left.0.cmp(&right.0));
        assert_eq!(
            key.iter().map(|(slot, _)| slot).collect::<Vec<_>>(),
            ports.iter().collect::<Vec<_>>()
        );
        explicit_closed += &spectator
            * value.to_owned()
            * function!(GS.u, 0, GS.cind(key[0].1))
            * function!(GS.ubar, 1, GS.cind(key[1].1));
        expected_components.insert(key, value.to_owned());
    }
    assert_eq!(expected_components.len(), 16);
    let joint = &spectator * &matrix * &states;
    assert!(
        NumeratorJetCache::default()
            .canonicalize(&(&joint - explicit_closed))?
            .is_zero,
        "the single closed joint equals its independent sixteen native matrix coefficients"
    );

    let family = symbol!("gammalooprs::uv::numerator_family");
    let entries = [(&spectator * &matrix), states.clone(), joint.clone()]
        .into_iter()
        .map(|rhs| -> Result<_> {
            let tag = DirectResidueBranches::numerator_scope().1;
            Ok(std::sync::Arc::new(NumeratorDefinition::try_from(
                FnMapEntry {
                    lhs: function!(family, &tag),
                    rhs,
                    args: vec![],
                    tags: vec![tag],
                    inlining: symbolica::evaluate::InliningPolicy::Always,
                    is_alias: false,
                },
            )?))
        })
        .collect::<Result<Vec<_>>>()?;
    assert_eq!(
        entries[0].rhs.expression() * entries[1].rhs.expression(),
        *entries[2].rhs.expression(),
        "one occurrence of each factor is the exact complete source joint"
    );
    let split: crate::uv::Integrands = [(
        crate::cff::CutCFFIndex::new_all_none(),
        &entries[0].lhs * &entries[1].lhs,
    )]
    .into_iter()
    .collect();
    let split = split.with_numerators(entries[..2].iter().cloned())?;
    let complete: crate::uv::Integrands = [(
        crate::cff::CutCFFIndex::new_all_none(),
        entries[2].lhs.clone(),
    )]
    .into_iter()
    .collect();
    let complete = complete.with_numerators([entries[2].clone()])?;
    let (bound, _, output_labels) = ForestTensorEvaluator::bind(&graph, &split)?;
    assert_eq!(output_labels.len(), 2);
    assert!(output_labels.iter().all(|(actual, label)| actual != label));
    let state_body = &bound
        .numerators()
        .iter()
        .find(|entry| entry.tags == entries[1].tags)
        .expect("retained state family")
        .rhs;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let state_ports = state_body
        .structure()
        .canonical()
        .external_structure()
        .into_iter()
        .map(|slot| cooking.uncook(slot.to_atom().as_view()))
        .collect::<Vec<_>>()
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    let expected_ports = ports
        .iter()
        .chain(output_labels.keys())
        .cloned()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        state_ports, expected_ports,
        "the state family retains both live connection ports and two distinct output ports"
    );
    // A zero matrix still owns the original ordered connection ports. The
    // nonzero state family independently supplies the amputation history and
    // therefore the new execution output frame; no zero syntax invents ports.
    let matrix_ports = entries[0].rhs.structure().logical_slots();
    let zero_matrix = entries[0].rhs.with_checked_expression(Atom::Zero)?;
    assert!(!zero_matrix.is_scalar());
    assert_eq!(zero_matrix.structure().logical_slots(), matrix_ports);
    let zero_split = split.clone().with_numerators([
        std::sync::Arc::new(entries[0].with_rhs(zero_matrix)),
        entries[1].clone(),
    ])?;
    let (zero_bound, _, zero_output_labels) = ForestTensorEvaluator::bind(&graph, &zero_split)?;
    let zero_matrix_body = &zero_bound
        .numerators()
        .iter()
        .find(|entry| entry.tags == entries[0].tags)
        .expect("retained zero matrix family")
        .rhs;
    assert!(zero_matrix_body.expression().is_zero());
    assert_eq!(zero_matrix_body.structure().logical_slots(), matrix_ports);
    assert_eq!(zero_output_labels.len(), 2);
    assert!(
        zero_output_labels
            .iter()
            .all(|(actual, label)| actual != label)
    );
    let zero_state_body = &zero_bound
        .numerators()
        .iter()
        .find(|entry| entry.tags == entries[1].tags)
        .expect("retained state family next to zero matrix")
        .rhs;
    let zero_state_ports = zero_state_body
        .structure()
        .canonical()
        .external_structure()
        .into_iter()
        .map(|slot| cooking.uncook(slot.to_atom().as_view()))
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        zero_state_ports,
        ports
            .iter()
            .chain(zero_output_labels.keys())
            .cloned()
            .collect::<std::collections::BTreeSet<_>>(),
        "the actual state family owns the output frame even when its matrix partner is zero"
    );

    // In contrast, zeroing the already closed joint loses its wavefunction
    // history. Preserve that scalar interface; do not assert an amputated
    // matrix or call the physical-output evaluator on this provenance-free zero.
    let closed_zero = entries[2].rhs.with_checked_expression(Atom::Zero)?;
    assert!(closed_zero.is_scalar());
    let closed_zero = complete
        .clone()
        .with_numerators([std::sync::Arc::new(entries[2].with_rhs(closed_zero))])?;
    let (closed_zero, _, closed_zero_labels) = ForestTensorEvaluator::bind(&graph, &closed_zero)?;
    assert!(closed_zero_labels.is_empty());
    assert!(closed_zero.numerators()[0].rhs.expression().is_zero());
    assert!(
        closed_zero.numerators()[0]
            .rhs
            .structure()
            .logical_slots()
            .is_empty()
    );

    let mut zero_split = ForestTensorEvaluator::new_shared(
        DirectForestFixture::OsQuarkInQuark,
        &graph,
        2,
        &zero_split,
    )?;
    let mut split =
        ForestTensorEvaluator::new_shared(DirectForestFixture::OsQuarkInQuark, &graph, 0, &split)?;
    let mut complete = ForestTensorEvaluator::new_shared(
        DirectForestFixture::OsQuarkInQuark,
        &graph,
        1,
        &complete,
    )?;
    for source in [&split, &complete] {
        let ForestEvaluationBackend::Shared { evaluator, .. } = &source.backend else {
            unreachable!()
        };
        assert!(
            evaluator
                .fn_map_entries
                .iter()
                .any(|entry| entry.rhs.contains(power.as_view())),
            "the compiled joint retains the same opaque positive-power factor"
        );
    }
    for external in [[3, 1, 2, 4], [7, -2, 1, 5]] {
        let point = external.map(Atom::num);
        let loops = [[Atom::num(9), Atom::num(12), Atom::Zero]];
        let zeros = zero_split.evaluate(&point, &loops)?;
        assert_eq!(
            zeros.keys().collect::<Vec<_>>(),
            expected_components.keys().collect::<Vec<_>>()
        );
        let zero = ForestNumber::new(F(ArbPrec::default()), F(ArbPrec::default()));
        assert!(zeros.values().all(|value| value == &zero));
        let actual = split.evaluate(&point, &loops)?;
        let same_source = complete.evaluate(&point, &loops)?;
        assert_eq!(actual.len(), 16);
        assert_eq!(
            actual, same_source,
            "separate matrix/state families retain all output coordinates"
        );
        for (key, coefficient) in &expected_components {
            let expected = (&spectator * coefficient)
                .replace(p0.to_pattern())
                .with(external[0])
                .replace(p3.to_pattern())
                .with(external[3]);
            let mut evaluator = expected
                .evaluator::<Atom>(&[])
                .build()
                .map_err(|error| eyre!(error))?
                .map_coeff(&|value| {
                    ForestNumber::new(F::<ArbPrec>::from(&value.re), F::<ArbPrec>::from(&value.im))
                });
            assert_eq!(
                actual[key],
                evaluator.evaluate_single(&[]),
                "exact native matrix component {key:?}"
            );
        }
    }
    Ok(())
}

#[test]
fn shared_forest_oracle_preserves_factored_matrix_functions() -> Result<()> {
    use crate::integrands::process::param_builder::{FnMapEntry, NumeratorDefinition};
    use crate::utils::{ArbPrec, F};
    use spenso::structure::{TensorStructure, slot::IsAbstractSlot};
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    let reduced = graph.full_filter().subtract(&graph.initial_state_cut);
    let physical_tensor = graph
        .numerator(&reduced, &graph.empty_subgraph())?
        .state
        .expr
        .multiply(graph.global_numerator()?)?;
    let physical_tag = DirectResidueBranches::numerator_scope().1;
    let physical_lhs = function!(symbol!("gammalooprs::uv::numerator_family"), &physical_tag);
    let physical: crate::uv::Integrands = [(
        crate::cff::CutCFFIndex::new_all_none(),
        physical_lhs.clone(),
    )]
    .into_iter()
    .collect();
    let physical = physical.with_numerators([std::sync::Arc::new(
        crate::integrands::process::param_builder::NumeratorDefinition {
            lhs: physical_lhs,
            rhs: physical_tensor,
            args: vec![],
            tags: vec![physical_tag],
            inlining: symbolica::evaluate::InliningPolicy::Always,
            is_alias: false,
        },
    )])?;
    let (physical, _, output_labels) = ForestTensorEvaluator::bind(&graph, &physical)?;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let ports = physical.numerators()[0]
        .rhs
        .structure()
        .canonical()
        .external_structure()
        .into_iter()
        .map(|slot| cooking.uncook(slot.to_atom().as_view()))
        .map(|port| output_labels.get(&port).cloned().unwrap_or(port))
        .collect::<Vec<_>>();
    eyre::ensure!(
        ports.len() == 2,
        "the raw fermion fixture has two open spinor ports"
    );
    let argument = Atom::var(symbolica::symbol!("shared_forest_oracle_argument"));
    let tag = DirectResidueBranches::numerator_scope().1;
    let lhs = function!(
        symbolica::symbol!("gammalooprs::uv::numerator_family"),
        &tag,
        &argument
    );
    let rhs = (&argument + Atom::one()).pow(11)
        * (spenso::g!(&ports[0], &ports[1])
            + Atom::i() * &argument * idenso::gamma0!(&ports[0], &ports[1]));
    let entry = std::sync::Arc::new(NumeratorDefinition::try_from(FnMapEntry {
        lhs: lhs.clone(),
        rhs: rhs.clone(),
        args: vec![
            symbolica::atom::Indeterminate::try_from(argument.clone())
                .map_err(|error| eyre!(error))?,
        ],
        tags: vec![tag],
        inlining: symbolica::evaluate::InliningPolicy::Always,
        is_alias: false,
    })?);
    let atom = lhs
        .replace(argument.clone())
        .with(GS.emr_mom(EdgeIndex(0), GS.cind(0)))
        + Atom::num(2)
            * lhs
                .replace(argument.clone())
                .with(GS.emr_mom(EdgeIndex(0), GS.cind(3)));
    let input: crate::uv::Integrands = [(crate::cff::CutCFFIndex::new_all_none(), atom)]
        .into_iter()
        .collect();
    let input = input.with_numerators([entry])?;
    let resolved: Atom = input.resolved()?.iter().map(|(_, atom)| atom).sum();
    let mut reference =
        ForestTensorEvaluator::new(DirectForestFixture::OsQuarkInQuark, &graph, 0, &resolved)?;
    let mut shared =
        ForestTensorEvaluator::new_shared(DirectForestFixture::OsQuarkInQuark, &graph, 0, &input)?;
    assert_eq!(
        NumeratorInterface::export(&input.numerators()[0].rhs).0,
        rhs,
        "preparation must preserve the original power and function body"
    );
    let ForestEvaluationBackend::Shared { evaluator, .. } = &shared.backend else {
        unreachable!()
    };
    let mut retained_power = false;
    for entry in &evaluator.fn_map_entries {
        entry.rhs.visitor(&mut |part| {
            if let AtomView::Pow(power) = part {
                retained_power |= i64::try_from(power.get_exp()).is_ok_and(|power| power == 11)
                    && matches!(power.get_base(), AtomView::Add(_));
            }
            true
        });
    }
    assert!(
        retained_power,
        "compiled shared functions must retain the powered numerator factor"
    );
    for external in [[3, 1, 2, 4], [7, -2, 1, 5]] {
        let external = external.map(Atom::num);
        let loops = [[Atom::num(9), Atom::num(12), Atom::Zero]];
        let expected = reference.evaluate(&external, &loops)?;
        let actual = shared.evaluate(&external, &loops)?;
        let zero = ForestNumber::new(F(ArbPrec::default()).zero(), F(ArbPrec::default()).zero());
        for key in expected
            .keys()
            .chain(actual.keys())
            .collect::<std::collections::BTreeSet<_>>()
        {
            assert_eq!(
                actual.get(key).unwrap_or(&zero),
                expected.get(key).unwrap_or(&zero),
                "exact Gaussian matrix component {key:?}"
            );
        }
    }
    Ok(())
}

#[test]
fn shared_forest_oracle_preserves_kinematic_parameter_order() -> Result<()> {
    use crate::integrands::process::param_builder::{FnMapEntry, NumeratorDefinition};
    use crate::utils::{ArbPrec, F, FloatLike};
    use spenso::structure::{TensorStructure, slot::IsAbstractSlot};
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    let reduced = graph.full_filter().subtract(&graph.initial_state_cut);
    let physical_tensor = graph
        .numerator(&reduced, &graph.empty_subgraph())?
        .state
        .expr
        .multiply(graph.global_numerator()?)?;
    let physical_tag = DirectResidueBranches::numerator_scope().1;
    let physical_lhs = function!(symbol!("gammalooprs::uv::numerator_family"), &physical_tag);
    let physical: crate::uv::Integrands = [(
        crate::cff::CutCFFIndex::new_all_none(),
        physical_lhs.clone(),
    )]
    .into_iter()
    .collect();
    let physical = physical.with_numerators([std::sync::Arc::new(
        crate::integrands::process::param_builder::NumeratorDefinition {
            lhs: physical_lhs,
            rhs: physical_tensor,
            args: vec![],
            tags: vec![physical_tag],
            inlining: symbolica::evaluate::InliningPolicy::Always,
            is_alias: false,
        },
    )])?;
    let (physical, _, output_labels) = ForestTensorEvaluator::bind(&graph, &physical)?;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let ports = physical.numerators()[0]
        .rhs
        .structure()
        .canonical()
        .external_structure()
        .into_iter()
        .map(|slot| cooking.uncook(slot.to_atom().as_view()))
        .map(|port| output_labels.get(&port).cloned().unwrap_or(port))
        .collect::<Vec<_>>();
    eyre::ensure!(
        ports.len() == 2,
        "the raw quark fixture has two spinor ports"
    );
    let argument = Atom::var(symbolica::symbol!("shared_coordinate_order_argument"));
    let tag = DirectResidueBranches::numerator_scope().1;
    let lhs = function!(
        symbolica::symbol!("gammalooprs::uv::numerator_family"),
        &tag,
        &argument
    );
    let weights = [2, 3, 5, 7, 11, 13, 17];
    let mut weighted = Atom::num(31) * Atom::var(GS.pi);
    for (component, weight) in weights[..4].iter().enumerate() {
        weighted += Atom::num(*weight)
            * GS.emr_mom(EdgeIndex(usize::from(component == 3)), GS.cind(component));
    }
    let loop_edge = *graph.loop_momentum_basis.loop_edges.first().unwrap();
    for component in 1..=3 {
        weighted += Atom::num(weights[3 + component]) * GS.emr_mom(loop_edge, GS.cind(component));
    }
    let entry = std::sync::Arc::new(NumeratorDefinition::try_from(FnMapEntry {
        lhs: lhs.clone(),
        rhs: &argument * spenso::g!(&ports[0], &ports[1]),
        args: vec![
            symbolica::atom::Indeterminate::try_from(argument.clone())
                .map_err(|error| eyre!(error))?,
        ],
        tags: vec![tag],
        inlining: symbolica::evaluate::InliningPolicy::Always,
        is_alias: false,
    })?);
    let input: crate::uv::Integrands = [(
        crate::cff::CutCFFIndex::new_all_none(),
        lhs.replace(argument).with(weighted),
    )]
    .into_iter()
    .collect();
    let input = input.with_numerators([entry])?;
    let mut evaluator =
        ForestTensorEvaluator::new_shared(DirectForestFixture::OsQuarkInQuark, &graph, 0, &input)?;
    for (numerators, rational_part) in [
        ([1, -2, 3, -5, 7, -11, 13], 6),
        ([-3, 5, -7, 11, -13, 17, -19], -9),
    ] {
        let coordinates = std::array::from_fn::<_, 7, _>(|index| {
            Atom::num(numerators[index]) / Atom::num(weights[index])
        });
        let external = std::array::from_fn(|index| coordinates[index].clone());
        let loops = [std::array::from_fn(|index| coordinates[4 + index].clone())];
        let values = evaluator.evaluate(&external, &loops)?;
        let zero = F(ArbPrec::default());
        let expected = ForestNumber::new(
            F::from(&symbolica::domains::rational::Rational::from(31)) * zero.clone().pi()
                + F::from(&symbolica::domains::rational::Rational::from(rational_part)),
            zero.clone(),
        );
        for (key, value) in values {
            assert_eq!(key.len(), 2);
            if key[0].1 == key[1].1 {
                assert!(
                    (value.clone() - &expected).norm().re
                        <= F(ArbPrec::default().epsilon()).sqrt() * expected.clone().norm().re,
                    "shared matrix diagonal must give 31*pi+{rational_part}"
                );
            } else {
                assert_eq!(value, ForestNumber::new(zero.clone(), zero.clone()));
            }
        }
    }
    Ok(())
}

#[test]
fn forest_oracle_amputation_requires_external_spinor_ownership() -> Result<()> {
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::slot::{DummyAind, IsAbstractSlot},
    };
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    for (edge, external) in [(EdgeIndex(0), true), (EdgeIndex(2), false)] {
        let pair = &graph[&edge].1;
        let port = graph.underlying[pair.any_hedge()]
            .edge_spin_slots()
            .next()
            .expect("the chosen quark edge has a spinor port");
        let wavefunction = function!(GS.u, usize::from(edge), port.to_atom());
        // Amputation of a free standalone vector is an identity tensor. Use
        // an actual contracted coefficient here to test external ownership.
        let output = idenso::bis!(4, Atom::from(Aind::new_dummy()));
        let matrix = SPENSO_TAG.tensor_symbol("amputation_owner::matrix");
        let input = [(
            crate::cff::CutCFFIndex::new_all_none(),
            function!(matrix, &output, port.to_atom()) * wavefunction,
        )]
        .into_iter()
        .collect();
        let bound = ForestTensorEvaluator::bind(&graph, &input);
        if external {
            let (bound, _, output_labels) = bound?;
            let actual: Atom = bound.iter().map(|(_, atom)| atom).sum();
            let expected = function!(matrix, &output, port.to_atom()).replace_multiple(
                output_labels.iter().map(|(actual, label)| {
                    symbolica::id::Replacement::new(label.to_pattern(), actual.clone())
                }),
            );
            assert!(
                NumeratorJetCache::default()
                    .canonicalize(&(actual - expected))?
                    .is_zero
            );
        } else {
            match bound {
                Err(error) => assert!(error.to_string().contains("external spinor port")),
                Ok(_) => panic!("an internal propagator cannot be amputated as an external state"),
            }
        }
    }
    Ok(())
}

#[test]
fn forest_oracle_amputation_keeps_additive_index_scopes() -> Result<()> {
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::slot::{DummyAind, IsAbstractSlot},
    };
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    let ports = [EdgeIndex(0), EdgeIndex(1)].map(|edge| {
        graph.underlying[graph[&edge].1.any_hedge()]
            .edge_spin_slots()
            .next()
            .expect("one graph-owned spinor port")
            .to_atom()
    });
    let internal = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    let matrix = SPENSO_TAG.tensor_symbol("amputation_scopes::matrix");
    let bridge = SPENSO_TAG.tensor_symbol("amputation_scopes::bridge");
    let other = SPENSO_TAG.tensor_symbol("amputation_scopes::other");
    let mass = Atom::var(symbol!("amputation_scopes::mass"; Scalar));
    // `internal` labels a contracted matrix port in the first Add branch,
    // but the u(0) wavefunction in the second. Globally renaming both it and
    // ports[0] to one output captures the first branch's private contraction.
    let scoped = (function!(matrix, &internal, &ports[0]) * function!(GS.u, 0, &ports[0])
        + &mass * function!(GS.u, 0, &internal))
        * function!(bridge, &ports[1], &internal)
        * function!(GS.ubar, 1, &ports[1]);
    let scoped_expected = function!(matrix, &internal, &ports[0])
        * function!(bridge, &ports[1], &internal)
        + &mass * function!(bridge, &ports[1], &ports[0]);
    // Literal source index names may also belong to different graph-owned
    // external states in different summands. Each occurrence owns its edge.
    let swapped = Atom::num(3)
        * function!(GS.u, 0, &ports[0])
        * function!(matrix, &ports[0], &ports[1])
        * function!(GS.ubar, 1, &ports[1])
        + Atom::num(5)
            * function!(GS.u, 0, &ports[1])
            * function!(other, &ports[0], &ports[1])
            * function!(GS.ubar, 1, &ports[0]);
    let swapped_expected = Atom::num(3) * function!(matrix, &ports[0], &ports[1])
        + Atom::num(5) * function!(other, &ports[1], &ports[0]);
    let contracted_pair = function!(GS.u, 0, &internal) * function!(GS.ubar, 1, &internal);
    let x = Atom::var(symbol!("amputation_scopes::x"; Scalar));
    let y = Atom::var(symbol!("amputation_scopes::y"; Scalar));
    let opaque = symbol!("amputation_scopes::opaque"; Scalar).call((&x + &y).pow(11));
    let spectator = &opaque / (&x - &y).pow(2);
    for (label, source, expected) in [
        ("external versus private slot", scoped, scoped_expected),
        ("swapped external slots", swapped, swapped_expected),
        (
            "contracted abstract states",
            contracted_pair,
            spenso::g!(&ports[0], &ports[1]),
        ),
    ] {
        let input = [(crate::cff::CutCFFIndex::new_all_none(), &spectator * source)]
            .into_iter()
            .collect();
        let (bound, _, output_labels) = ForestTensorEvaluator::bind(&graph, &input)?;
        let actual: Atom = bound.iter().map(|(_, atom)| atom).sum();
        let expected = (&spectator * expected).replace_multiple(output_labels.iter().map(
            |(actual, label)| symbolica::id::Replacement::new(label.to_pattern(), actual.clone()),
        ));
        assert!(
            NumeratorJetCache::default()
                .canonicalize(&(&actual - &expected))?
                .is_zero,
            "{label}: exact open-matrix coefficients"
        );
        assert!(
            !NumeratorJetCache::default().canonicalize(&actual)?.is_zero,
            "{label}: nonzero matrix control"
        );
        assert!(
            actual.contains(opaque.as_view()),
            "opaque powered spectator stays literal"
        );
    }
    Ok(())
}

#[test]
fn forest_oracle_amputation_opens_native_fixed_spinor_components() -> Result<()> {
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{
            representation::LibrarySlot,
            slot::{DummyAind, IsAbstractSlot},
        },
    };
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    let ports = [EdgeIndex(0), EdgeIndex(1)].map(|edge| {
        let pair = &graph[&edge].1;
        graph.underlying[pair.any_hedge()]
            .edge_spin_slots()
            .next()
            .expect("the external quark has one native spinor port")
            .to_atom()
    });
    let diagonal =
        Atom::add_many((0..4).map(|index| {
            function!(GS.u, 0, GS.cind(index)) * function!(GS.ubar, 1, GS.cind(index))
        }));
    let identity = spenso::g!(&ports[0], &ports[1]);
    let off_diagonal = function!(GS.u, 0, GS.cind(1)) * function!(GS.ubar, 1, GS.cind(2));
    let off_diagonal_expected = function!(GS.delta_vec, GS.cind(1), &ports[0])
        * function!(GS.delta_vec, GS.cind(2), &ports[1]);
    let tensor = SPENSO_TAG.tensor_symbol("amputation_exact::matrix");
    // Include a serialized dummy in the input; fresh output allocation must
    // reserve it instead of capturing it as a new graph-output index.
    let right = LibrarySlot::<Aind>::try_from(ports[1].as_view())?
        .rep()
        .slot::<Aind, _>(Aind::new_dummy_at(0))
        .to_atom();
    let abstract_pair = function!(GS.u, 0, &ports[0])
        * function!(tensor, &ports[0], &right)
        * function!(GS.ubar, 1, &right);
    let mixed_pair = function!(GS.u, 0, GS.cind(1))
        * function!(tensor, GS.cind(1), &right)
        * function!(GS.ubar, 1, &right);
    let mixed_expected =
        function!(GS.delta_vec, GS.cind(1), &ports[0]) * function!(tensor, GS.cind(1), &ports[1]);
    let x = Atom::var(symbol!("amputation_exact::x"; Scalar));
    let y = Atom::var(symbol!("amputation_exact::y"; Scalar));
    let opaque = symbol!("amputation_exact::opaque"; Scalar).call((&x + &y).pow(11));
    let spectator = &opaque / (&x - &y).pow(2);
    for (label, source, expected) in [
        (
            "complete component pair",
            diagonal.clone(),
            identity.clone(),
        ),
        ("off-diagonal pair", off_diagonal, off_diagonal_expected),
        ("mixed fixed and abstract pair", mixed_pair, mixed_expected),
        (
            "mixed additive forms",
            Atom::num(3) * diagonal + Atom::num(2) * abstract_pair,
            Atom::num(3) * identity + Atom::num(2) * function!(tensor, &ports[0], &ports[1]),
        ),
    ] {
        for coefficient in [Atom::one(), spectator.clone()] {
            let input = [(
                crate::cff::CutCFFIndex::new_all_none(),
                &coefficient * &source,
            )]
            .into_iter()
            .collect();
            let (bound, _, output_labels) = ForestTensorEvaluator::bind(&graph, &input)?;
            let actual: Atom = bound.iter().map(|(_, atom)| atom).sum();
            let expected = (&coefficient * &expected).replace_multiple(output_labels.iter().map(
                |(actual, label)| {
                    symbolica::id::Replacement::new(label.to_pattern(), actual.clone())
                },
            ));
            let certificate = NumeratorJetCache::default().canonicalize(&(&actual - &expected))?;
            assert!(
                certificate.is_zero,
                "{label}: exact native-component opening"
            );
            assert!(
                !NumeratorJetCache::default().canonicalize(&actual)?.is_zero,
                "{label}: nonzero tensor control"
            );
            if !coefficient.is_one() {
                assert!(
                    actual.contains(opaque.as_view()),
                    "opaque powered payload stays literal"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn forest_oracle_amputation_rejects_invalid_fixed_components_and_owners() -> Result<()> {
    use spenso::structure::slot::DummyAind;
    test_initialise()?;
    let (graph, _, _) = fixture(false)?;
    for (owner, component) in [
        (Atom::Zero, Atom::num(-1)),
        (Atom::Zero, Atom::num(4)),
        (Atom::Zero, Atom::num((1, 2))),
        (Atom::num(2), Atom::Zero), // paired internal quark
        (Atom::num(graph.underlying.n_edges()), Atom::Zero),
        (
            Atom::var(symbol!("amputation_exact::invalid_owner"; Scalar)),
            Atom::Zero,
        ),
    ] {
        let wavefunction = function!(GS.u, &owner, function!(AIND_SYMBOLS.cind, &component));
        let input = [(crate::cff::CutCFFIndex::new_all_none(), wavefunction)]
            .into_iter()
            .collect();
        assert!(
            ForestTensorEvaluator::bind(&graph, &input).is_err(),
            "invalid fixed wavefunction owner={owner}, component={component} must be rejected"
        );
    }
    for port in [
        spenso::mink!(4, Atom::from(Aind::new_dummy())),
        idenso::bis!(3, Atom::from(Aind::new_dummy())),
    ] {
        let input = [(
            crate::cff::CutCFFIndex::new_all_none(),
            function!(GS.u, 0, port),
        )]
        .into_iter()
        .collect();
        assert!(
            ForestTensorEvaluator::bind(&graph, &input).is_err(),
            "an abstract external state must match its graph-owned spinor representation"
        );
    }
    Ok(())
}

#[test]
fn direct_soft_gluon_shared_families_agree_at_identical_point() -> Result<()> {
    use crate::utils::{ArbPrec, F, FloatLike};
    test_initialise()?;
    let fixture = DirectForestFixture::SoftGluonInQuark;
    let (graph, inputs) = fixture.assemble(
        DirectForestTerms::CompleteForest,
        false,
        ["1", "0", "0", "0"],
    )?;
    let mut evaluator = ForestTensorEvaluator::new_shared(fixture, &graph, 1, &inputs[1])?;
    let reference: Atom = inputs[1].resolved()?.iter().map(|(_, atom)| atom).sum();
    let mut reference = ForestTensorEvaluator::new(fixture, &graph, 1, &reference)?;
    let external = [Atom::num(5), Atom::Zero, Atom::Zero, Atom::Zero];
    let loops = [
        [Atom::num((1, 8)), Atom::Zero, Atom::Zero],
        [Atom::Zero, Atom::num(9), Atom::Zero],
    ];
    // Both compilers receive the same factorized route and binding chart;
    // only this independent reference resolves the shared family catalog.
    let baseline = reference.evaluate(&external, &loops)?;
    let actual = evaluator.evaluate(&external, &loops)?;
    let slots = |values: &ForestComponents| {
        values
            .keys()
            .flat_map(|key| key.iter().map(|(slot, _)| slot.clone()))
            .collect::<std::collections::BTreeSet<_>>()
    };
    assert_eq!(
        slots(&baseline),
        slots(&actual),
        "the shared evaluator retains all open slots"
    );
    let tensor_scale = baseline
        .values()
        .chain(actual.values())
        .map(|value| value.clone().norm().re)
        .fold(F(ArbPrec::default()), |left, right| left.max(right));
    assert!(tensor_scale > F(ArbPrec::default()).zero());
    let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
    let zero = ForestNumber::new(F(ArbPrec::default()), F(ArbPrec::default()));
    for key in baseline
        .keys()
        .chain(actual.keys())
        .collect::<std::collections::BTreeSet<_>>()
    {
        let expected = baseline.get(key).unwrap_or(&zero);
        let value = actual.get(key).unwrap_or(&zero);
        let scale = value.clone().norm().re.max(expected.clone().norm().re);
        assert!(
            (value.clone() - expected).norm().re <= &tolerance * (scale + &tensor_scale),
            "shared family evaluator changed tensor component {key:?}"
        );
    }
    Ok(())
}

#[test]
fn boosted_rainbow_cff_outer_soft_limit_distinguishes_fixed_points() -> Result<()> {
    use crate::utils::{ArbPrec, F};
    use symbolica::prelude::SingleFloat;

    test_initialise()?;
    let fixture = DirectForestFixture::OsQuarkInQuark;
    let unit = F(ArbPrec::default()).one();
    let zero = unit.zero();
    for (external, matched_velocity) in [
        ([12, 0, 0, 0], ["1", "0", "0", "0"]),
        ([15, 9, 0, 0], ["5/4", "3/4", "0", "0"]),
        ([15, 3, 6, 6], ["5/4", "1/4", "1/2", "1/2"]),
    ] {
        let external = external.map(Atom::num);
        for (matched, velocity) in [(false, ["1", "0", "0", "0"]), (true, matched_velocity)] {
            if !matched && matched_velocity == ["1", "0", "0", "0"] {
                continue;
            }
            // Include the actual finite child coefficient and its parent operation.
            // An overall analytic addback has a regular localization kernel and
            // cannot contribute to the negative soft powers tested here.
            let (graph, expressions) =
                fixture.assemble(DirectForestTerms::CompleteForest, true, velocity)?;
            assert_eq!(
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                [EdgeIndex(2), EdgeIndex(3)],
                "the soft ray uses the outer quark r=e2 and hard child quark l=e3"
            );
            let mut evaluator =
                ForestTensorEvaluator::new_shared(fixture, &graph, 2, &expressions[2])?;
            for direction in [[1, 0, 0], [0, 1, 0]] {
                let mut weighted = Vec::new();
                for exponent in 3..=7 {
                    let denominator = 1usize << exponent;
                    let radius = Atom::one() / Atom::num(denominator);
                    let loops = [
                        std::array::from_fn(|index| {
                            &external[index + 1] + Atom::num(direction[index]) * &radius
                        }),
                        [Atom::Zero, Atom::num(9), Atom::Zero],
                    ];
                    let values = evaluator.evaluate(&external, &loops)?;
                    let norm = values.values().try_fold(zero.clone(), |norm, value| {
                        let component = value.clone().norm().re;
                        eyre::ensure!(component.is_finite(), "finite complete rainbow tensor");
                        Ok::<_, eyre::Report>(norm.max(component))
                    })?;
                    assert!(norm > zero, "nonzero complete rainbow control");
                    let numeric_radius = &unit / unit.from_usize(denominator);
                    weighted.push(norm * numeric_radius.powi(2));
                }
                crate::debug_tags!(#uv, #local, #soft, #summary;
                    stage = "boosted_rainbow_cff_outer_soft", matched,
                    external = ?external,
                    direction = ?direction, file.weighted = ?weighted,
                    "Complete energy-integrated rainbow radial control"
                );
                let boundary = unit.from_usize(3) / unit.from_usize(2);
                for tail in weighted[2..].windows(2) {
                    if matched {
                        assert!(
                            tail[1] <= &boundary * &tail[0],
                            "matching fixed point removes the logarithmic radial pole: {weighted:?}"
                        );
                    } else {
                        assert!(
                            tail[1] > &boundary * &tail[0]
                                && tail[1] <= unit.from_usize(5) / unit.from_usize(2) * &tail[0],
                            "rest-frame subtraction retains the boosted radial pole: {weighted:?}"
                        );
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn boosted_complete_rainbow_with_finite_child_agrees_across_three_routes() -> Result<()> {
    // One fixed component point, with matching and independently chosen live
    // external momenta. This supplements the projected-route radial matrix.
    DirectForestFixture::OsQuarkInQuark.check_points(
        DirectForestTerms::CompleteForest,
        true,
        &[[15, 3, 6, 6], [15, 9, 0, 0]],
        ["5/4", "1/4", "1/2", "1/2"],
    )
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
        .numerator(current.subgraph(), root.subgraph())?
        .state
        .expr;
    let options = graph.denominator_only_cff_3d_expression_options();
    let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
    let production = graph.generate_3d_expression_for_integrand(
        &[],
        &surfaces,
        &options,
        Some(numerator.expression()),
    )?;
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
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
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
                .with(GS.m_uv_vacuum);
            let cooking = CookSettings::indices()
                .with_mode(CookMode::ReversibleEncoding)
                .with_representation_payloads(true, true);
            let value = numerator
                .with_checked_expression(cooking.try_cook(value.as_view())?)?
                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                })?;
            let value = cooking.uncook(value.expression().as_view());
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
    for velocity in [
        ["1", "0", "0", "0"],
        ["5/4", "3/4", "0", "0"],
        ["5/4", "1/4", "1/2", "1/2"],
    ] {
        let (mut graph, current, root) = fixture(false)?;
        let settings = UVgenerationSettings {
            on_shell_subtraction_velocity: velocity.map(str::to_owned),
            generate_integrated: false,
            add_marker: false,
            ..Default::default()
        };
        let bare_numerator = graph
            .numerator(current.subgraph(), root.subgraph())?
            .state
            .expr;
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonical_surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &canonical_surfaces,
            &options,
            Some(bare_numerator.expression()),
        )?;
        let pattern = OrientationPattern::default();
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation =
            OrientationProjection::exact_expression(&production, &options, &pattern, false);
        let localizer = Localizer::new(&cutset, orientation);
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
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
        for sector in local.projection_sectors()? {
            let canonical = sector.canonical_projection(&graph)?;
            let coefficient = localizer
                .project_factorized_taylor_sector(&mut graph, &canonical, &mut projection_context)?
                .resolved()?;
            projected += coefficient
                .iter()
                .map(|(_, value)| value.clone())
                .sum::<Atom>();
        }

        let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
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
                &UVCtx::new(&graph, &graph.model, &settings),
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
            let specialize = |value: &Atom| -> Result<Atom> {
                // Match the final evaluator boundary on both routes: direct
                // source roots retain tree ownership wrappers, while the typed
                // projection retains symbolic colour Casimirs until finalization.
                // Use the established colour invariant rules, not a numeric CF
                // replacement inferred from the difference under examination.
                let value = GS
                    .erase_uv_momentum_provenance(value)
                    .replace(function!(GS.tree_denom_wrapper, W_.x_))
                    .with(W_.x_);
                let cooking = CookSettings::indices()
                    .with_mode(CookMode::ReversibleEncoding)
                    .with_representation_payloads(true, true);
                let value = bare_numerator
                    .with_checked_expression(cooking.try_cook(value.as_view())?)?
                    .simplify_algebra(&idenso::tensor::AlgebraSettings {
                        color: Some(
                            ColorSimplifySettings::default().with_cof_dimension_invariants(),
                        ),
                        ..Default::default()
                    })?;
                let value = cooking
                    .uncook(value.expression().as_view())
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
                Ok(GS.lower_esurfaces(&value))
            };
            let reference = specialize(&projected)?;
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
                let difference = specialize(explicit)? - specialize(&direct_host_values[1][host])?;
                let certificate =
                    NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
                assert!(
                    certificate.is_zero,
                    "direct OS velocity {velocity:?}, orientation host {host:?}, k={spatial:?}, p={external:?}: {}",
                    certificate.key
                );
            }
            for (mode, direct) in direct_values.iter().enumerate() {
                let actual = specialize(direct)?;
                let difference = &actual - &reference;
                let certificate =
                    NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
                assert!(
                    certificate.is_zero,
                    "direct OS velocity {velocity:?}, selector mode {mode}, k={spatial:?}, p={external:?}: {}",
                    certificate.key
                );
            }
        }
    }
    Ok(())
}
