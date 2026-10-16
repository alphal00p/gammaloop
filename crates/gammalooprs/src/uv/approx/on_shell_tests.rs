//! Independent one-loop OS oracles in explicitly selected homogeneous LMBs.
//! Only denominator polynomials and finite tensor components are expanded.

use super::*;
use crate::{
    graph::parse::NumeratorInterface,
    initialisation::test_initialise,
    model::{InputParamCard, ModelGammaLoopExt},
    utils::{FUN_LIB, TENSORLIB},
    uv::{
        Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode,
        overlap_numerators::NumeratorJetCache,
    },
};
use idenso::{
    color::ColorSimplifySettings,
    dirac::GammaSimplifySettings,
    representations::Bispinor,
    shorthands::UndoShorthands,
    tensor::{AlgebraSettings, SymbolicNetExt, SymbolicNetParse},
};
use linnet::half_edge::involution::{Flow, Hedge};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
        tags::SPENSO_TAG,
    },
    structure::{
        TensorStructure,
        representation::{Minkowski, RepName},
    },
};

pub(crate) fn fixture(massless: bool) -> Result<(Graph, OwnedForestNode, OwnedForestNode)> {
    let source = if massless {
        include_str!("../../../../../tests/resources/graphs/local_ir_quark_self_energy.dot")
    } else {
        include_str!("../../../../../tests/resources/graphs/local_os_top_self_energy.dot")
    };
    let mut model = crate::utils::load_generic_model("sm");
    if massless {
        let mut parameters = InputParamCard::from_str(
            include_str!(
                "../../../../../assets/models/json/sm/restrict_light_quarks_massless.json"
            )
            .to_owned(),
            "json",
        )?;
        model.simplify(&mut parameters)?;
        assert!(
            model.particle("d")?.is_massless(&model),
            "the massless OS fixture requires an exact zero quark mass"
        );
    }
    let mut graphs = Graph::import_dot(source, &model)?;
    assert_eq!(graphs.len(), 1, "the OS fixture contains one graph");
    let graph = graphs.remove(0);
    let external: SuBitGraph = graph.external_filter();
    let component = graph.full_filter().subtract(&external);
    let current = OwnedForestNode {
        spinney: Spinney::with_scheme(
            InternalSubGraph::cleaned_filter_optimist(component, graph.as_ref()),
            &graph,
            &graph.loop_momentum_basis,
            ApproximationType::OS,
            1,
        )?
        .expect("the quark self-energy is UV divergent"),
        topo_order: 1,
    };
    let root = OwnedForestNode {
        spinney: Spinney::empty(&graph),
        topo_order: 0,
    };
    Ok((graph, current, root))
}

pub(super) fn scalar_component(atom: &Atom) -> Result<Atom> {
    // Production intentionally keeps scalar inner products opaque. The
    // independent finite-component oracle must materialize them explicitly,
    // reserving fresh indices on this whole factorized diagnostic copy.
    let explicit_dots = atom.undo_dots::<Aind>()?;
    let mut network = explicit_dots.parse_into_net()?;
    network
        .execute::<Sequential, SmallestDegree, _, _, _>(&*TENSORLIB.read().unwrap(), &*FUN_LIB)?;
    Ok(match network.result_scalar()? {
        ExecutionResult::One => Atom::one(),
        ExecutionResult::Zero => Atom::Zero,
        ExecutionResult::Val(value) => value.into_owned(),
    })
}

fn physical_expression(graph: &Graph, current: &OwnedForestNode, atom: &Atom) -> Atom {
    let expression = GS
        .resolve_os_momenta(
            &GS.erase_uv_momentum_provenance(atom)
                .replace(function!(GS.ct_marker, W_.a_))
                .with(Atom::one())
                .replace_multiple(graph.uv_wrapped_replacement(
                    current.subgraph(),
                    current.lmb(),
                    &[W_.x___],
                )),
        )
        .replace(GS.dim)
        .with(Atom::num(4));
    debug_tags!(#uv, #local, #on_shell, #inspect;
        stage = "os_oracle_physical_interface",
        file.expression = %expression.to_canonical_string(),
        "Factorized OS oracle input before tensor-interface admission"
    );
    // The independent scalar derivatives below act on momentum coefficients.
    // Open only compact vector arguments first: differentiating gamma(..., z*p)
    // as an arbitrary symbolic function would manufacture unevaluated gamma
    // derivatives instead of the linear Clifford coefficient.
    let source = NumeratorInterface::admit(expression, None)
        .expect("the authored OS diagnostic has a valid tensor interface")
        .try_into_indexed()
        .expect("the authored OS diagnostic has explicitly indexed free ports");
    let normalized =
        finalize(source).expect("the physical OS diagnostic has valid tensor contractions");
    NumeratorInterface::export(&normalized)
        .0
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
        })
        .expect("the independent oracle can expose compact vector arguments")
        .simple_execute::<()>()
        .expect("the independent oracle preserves the factorized numerator")
}

/// Collect only identical denominator products. No common denominator or
/// polynomial expansion of a graph numerator is used, even in this diagnostic.
fn rational_coefficients(
    graph: &Graph,
    current: &OwnedForestNode,
    atom: &Atom,
) -> Result<BTreeMap<Vec<Atom>, Atom>> {
    let expression = physical_expression(graph, current, atom);
    let mut coefficients = BTreeMap::<Vec<Atom>, Atom>::new();
    for term in FourDTerm::from_view(expression.as_view())? {
        let mut denominators = term
            .denominators
            .iter()
            .map(|denominator| {
                scalar_component(&denominator.full_expr).map(|polynomial| polynomial.expand())
            })
            .collect::<Result<Vec<_>>>()?;
        denominators.sort();
        *coefficients.entry(denominators).or_insert(Atom::Zero) +=
            NumeratorInterface::export(&term.numerator).0;
    }
    Ok(coefficients)
}

/// Both policies execute the same exact finite-component certificate. Small
/// boosted soft coefficients retain compact matrix products when symbolic
/// gamma-chain collection would expand their fixed-vector sums.
#[derive(Clone, Copy, Debug)]
pub(crate) enum SpinReduction {
    SymbolicFirst,
    FiniteComponents,
}

pub(crate) fn assert_rational_zero(
    label: &str,
    graph: &Graph,
    current: &OwnedForestNode,
    atom: &Atom,
    spin_reduction: SpinReduction,
) -> Result<()> {
    let mut numerators = NumeratorJetCache::default();
    for (class, (denominators, numerator)) in rational_coefficients(graph, current, atom)?
        .into_iter()
        .enumerate()
    {
        // The finite oracle uses repeated-pair reduction for open spinors to
        // retain compact words. Closed coefficients keep canonical trace
        // reduction. The existing typed color pass owns both requested families;
        // finite component execution below supplies the remaining contractions.
        let source = NumeratorInterface::admit(numerator.clone(), None)?;
        let open_spinor = source
            .structure()
            .canonical()
            .external_structure_iter()
            .any(|slot| slot.rep_name() == Bispinor {}.into());
        let gamma_settings = if open_spinor {
            GammaSimplifySettings::repeated_pairs()
        } else {
            GammaSimplifySettings::canonical()
        };
        let reduced = source.simplify_algebra(&AlgebraSettings {
            gamma: matches!(spin_reduction, SpinReduction::SymbolicFirst).then_some(gamma_settings),
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            epsilon: true,
            ..Default::default()
        })?;
        let reduced = NumeratorInterface::export(&reduced).0;
        let certificate = numerators.canonicalize(&reduced.undo_dots::<Aind>()?)?;
        let is_zero = certificate.is_zero;
        if !is_zero {
            debug_tags!(#uv, #local, #four_d, #trace;
                stage = "os_exact_oracle_failure",
                label = label,
                denominator_class = class,
                file.denominators = ?denominators,
                file.numerator = %numerator.to_plain_string(),
                file.contracted_tensor = %certificate.key.to_plain_string(),
                file.rational_expression = %atom.to_plain_string(),
                "Nonzero factorized OS oracle tensor coefficient"
            );
        }
        assert!(
            is_zero,
            "{label}: nonzero tensor coefficient in denominator class {class}"
        );
    }
    Ok(())
}

fn gamma_zero_at_port(graph: &Graph, numerator: &Atom, hedge: usize) -> Atom {
    let port = graph.underlying[Hedge(hedge)]
        .edge_spin_slots()
        .next()
        .expect("the fixture's external quark has one spinor slot")
        .to_atom();
    let dummy = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    idenso::gamma0!(&dummy, &port) * numerator.replace(port).with(dummy)
}

/// Contract a fixed velocity with an ordered spinor port. Concrete Lorentz
/// indices on a gamma remain opaque to the tensor library, so use its delta
/// basis; the Minkowski contraction supplies the spatial minus signs.
pub(crate) fn gamma_velocity_at_port(
    graph: &Graph,
    numerator: &Atom,
    hedge: usize,
    velocity: &[Atom; 4],
) -> Atom {
    let port = graph.underlying[Hedge(hedge)]
        .edge_spin_slots()
        .next()
        .expect("the fixture's external quark has one spinor slot")
        .to_atom();
    let dummy = idenso::bis!(4, Atom::from(Aind::new_dummy()));
    let lorentz = spenso::mink!(4, Atom::from(Aind::new_dummy()));
    let (row, column) = if graph.flow(Hedge(hedge)) == Flow::Source {
        (&port, &dummy)
    } else {
        (&dummy, &port)
    };
    let slash = &velocity[0] * idenso::gamma0!(row.as_view(), column.as_view())
        + Atom::add_many(
            (1..4).map(|mu| &velocity[mu] * function!(GS.delta_vec, GS.cind(mu), &lorentz)),
        ) * idenso::gamma!(row.as_view(), column.as_view(), &lorentz);
    slash * numerator.replace(port).with(dummy)
}

fn dot(left: &Atom, right: &Atom) -> Atom {
    function!(SPENSO_TAG.dot, left, right)
}

/// The fixture routes the *quark* as k and its gluon as k-p. Its numerator
/// N(k,m) has no p dependence. Differentiating the gluon denominator directly
/// gives Eq. 4.24's vacuum term 2[N(k,0) k.p-m gamma0 N(k,0) k0]/D^3.
/// This oracle does not invoke either production Taylor projector.
fn one_loop_reference(
    graph: &Graph,
    current: &OwnedForestNode,
    root: &OwnedForestNode,
    velocity: &[Atom; 4],
) -> Result<(Atom, Atom, Atom)> {
    assert_eq!(
        current.lmb().loop_edges.iter().copied().collect::<Vec<_>>(),
        vec![EdgeIndex(2)]
    );
    let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
    let vector_slot = Minkowski {}.new_rep(4).to_symbolic([]);
    let k = function!(GS.emr_mom, 2, &vector_slot);
    let p = function!(GS.emr_mom, 0, &vector_slot);
    let k_dot_velocity: Atom = (0..4)
        .map(|mu| {
            Atom::num(if mu == 0 { 1 } else { -1 })
                * &velocity[mu]
                * GS.emr_mom(EdgeIndex(2), GS.cind(mu))
        })
        .sum();
    let k_squared = dot(&k, &k);
    let numerator = graph
        .numerator(current.subgraph(), root.subgraph())?
        .to_d_dim(4)?
        .color_simplify()
        .unwrap()
        .get_single_atom()
        .expect("the source quark numerator is a single factor");
    let numerator = physical_expression(graph, current, &numerator);
    let massless_numerator = numerator.replace(mass.to_pattern()).with(Atom::Zero);
    let left_gamma_numerator = gamma_velocity_at_port(graph, &massless_numerator, 1, velocity);
    let denominator = |edge, polynomial: Atom, mass_squared: Atom| {
        GS.den(edge, function!(GS.emr_mom, 2), mass_squared, polynomial)
    };
    let quark = denominator(2, &k_squared - mass.pow(2), mass.pow(2));
    let vacuum = denominator(
        2,
        &k_squared - Atom::var(GS.m_uv_vacuum).pow(2),
        Atom::var(GS.m_uv_vacuum).pow(2),
    );
    let on_shell = [1, -1].into_iter().fold(Atom::Zero, |sum, sigma| {
        let projected = (&numerator
            + Atom::num(sigma) * gamma_velocity_at_port(graph, &numerator, 1, velocity))
            / Atom::num(2);
        let gluon = denominator(
            3,
            &k_squared - Atom::num(2 * sigma) * &mass * &k_dot_velocity + mass.pow(2),
            Atom::Zero,
        );
        sum + projected / (&quark * gluon)
    });
    let remainder = Atom::num(2)
        * (&massless_numerator * dot(&k, &p) - &mass * left_gamma_numerator * &k_dot_velocity)
        / vacuum.pow(3);
    Ok((on_shell, remainder, numerator))
}

#[test]
fn on_shell_zero_operand_needs_no_inferred_spinor_port() -> Result<()> {
    test_initialise()?;
    for massless in [false, true] {
        let (graph, current, _) = fixture(massless)?;
        let zero = SymbolicTensor::infer_indexed(Atom::Zero)?;
        assert_eq!(
            graph.on_shell_branches(
                &graph.model,
                &zero,
                current.subgraph(),
                current.lmb(),
                &UVgenerationSettings::default().on_shell_velocity()?
            )?,
            [zero.clone(), zero]
        );
        let source = graph
            .numerator(current.subgraph(), &graph.empty_subgraph())?
            .state
            .expr;
        let ranked_zero = source.with_checked_expression(Atom::Zero)?;
        for branch in graph.on_shell_branches(
            &graph.model,
            &ranked_zero,
            current.subgraph(),
            current.lmb(),
            &UVgenerationSettings::default().on_shell_velocity()?,
        )? {
            assert_eq!(branch, ranked_zero);
        }
        assert!(
            graph
                .on_shell_branches(
                    &graph.model,
                    &SymbolicTensor::infer_indexed(Atom::one())?,
                    current.subgraph(),
                    current.lmb(),
                    &UVgenerationSettings::default().on_shell_velocity()?,
                )
                .is_err()
        );
    }
    Ok(())
}

#[test]
fn exact_os_coefficients_use_fundamental_color_dimension_identities() -> Result<()> {
    test_initialise()?;
    let spectator = parse_lit!((color_a + color_b) * (color_c - color_d) ^ 11);
    for dimension in [3, 5] {
        let casimir = idenso::color_cas!(2, idenso::cof!(dimension as usize));
        let index = idenso::color_idx!(2, idenso::cof!(dimension as usize));
        let dimension_ratio = Atom::num(dimension * dimension - 1) / Atom::num(dimension);
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            ..Default::default()
        };
        let difference = casimir.pow(2) - &dimension_ratio * &casimir * index;
        assert!(
            SymbolicTensor::infer(&spectator * difference)?
                .simplify_algebra(&settings)?
                .expression()
                .is_zero(),
        );
        assert_eq!(
            SymbolicTensor::infer(&spectator * casimir)?
                .simplify_algebra(&settings)?
                .expression(),
            &(&spectator * dimension_ratio / Atom::num(2)),
        );
    }
    Ok(())
}

#[test]
fn symbolic_spin_reduction_preserves_os_phase_order_and_scalar_factors() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let (_, _, numerator) = one_loop_reference(
        &graph,
        &current,
        &root,
        &UVgenerationSettings::default().on_shell_velocity()?,
    )?;
    let spectator = parse_lit!((spin_a + spin_b) * (spin_c + spin_d) ^ 11);
    let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
    let mut tensors = NumeratorJetCache::default();
    for original in [numerator.clone(), gamma_zero_at_port(&graph, &numerator, 1)] {
        for gamma_settings in [
            GammaSimplifySettings::canonical(),
            GammaSimplifySettings::repeated_pairs(),
        ] {
            // This is the same typed admission used by the finite oracle, with
            // finite component execution as the independent spin-matrix check.
            let reduced = SymbolicTensor::infer(cooking.try_cook(original.as_view())?)?
                .simplify_algebra(&AlgebraSettings {
                    gamma: Some(gamma_settings),
                    epsilon: true,
                    ..Default::default()
                })?;
            let reduced = cooking.uncook(reduced.expression().as_view());
            assert!(
                tensors
                    .canonicalize(&(&original - &reduced).undo_dots::<Aind>()?)?
                    .is_zero,
                "symbolic spin reduction must preserve the exact QCD phase and ordered gamma0 insertion"
            );
            let factored =
                SymbolicTensor::infer(cooking.try_cook((&spectator * &original).as_view())?)?
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(gamma_settings),
                        epsilon: true,
                        ..Default::default()
                    })?;
            assert_eq!(
                cooking.uncook(factored.expression().as_view()),
                &spectator * &reduced
            );
            let momentum =
                |edge| function!(GS.emr_mom, edge, Minkowski {}.new_rep(4).to_symbolic([]));
            let scalar_product = dot(&momentum(2), &momentum(0));
            for power in [2, 3] {
                assert_rational_zero(
                    "finite open-spinor proof preserves powered dots and ordered gamma0",
                    &graph,
                    &current,
                    &(&spectator * scalar_product.pow(power) * (&original - &reduced)),
                    SpinReduction::SymbolicFirst,
                )?;
            }
        }
        // The production analytic boundary retains its canonical policy.
        let source = NumeratorInterface::admit(original.clone(), None)?.try_into_indexed()?;
        let analytic = crate::uv::approx::integrated::simplify(&source)?;
        let analytic = NumeratorInterface::export(&analytic).0;
        assert!(
            tensors
                .canonicalize(&(&original - &analytic).undo_dots::<Aind>()?)?
                .is_zero
        );
        assert_eq!(
            NumeratorInterface::export(&crate::uv::approx::integrated::simplify(
                &NumeratorInterface::admit(&spectator * &original, None)?.try_into_indexed()?,
            )?)
            .0,
            &spectator * analytic
        );
    }
    Ok(())
}

#[test]
fn production_os_matches_equation_4_24_with_ordered_projectors_and_one_forest_sign() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &graph.model, &settings);
    let seed = Full4dCts::root();
    let actual = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
    let (on_shell, remainder, numerator) = one_loop_reference(
        &graph,
        &current,
        &root,
        &UVgenerationSettings::default().on_shell_velocity()?,
    )?;
    assert_rational_zero(
        "production -T_OS equals the independent compact expression",
        &graph,
        &current,
        &(NumeratorInterface::export(actual.atom()).0 + &on_shell + &remainder),
        SpinReduction::SymbolicFirst,
    )?;
    let scheme = actual
        .active_sectors()
        .iter()
        .fold(Atom::Zero, |sum, sector| {
            sum + NumeratorInterface::export(sector.integration_atom()).0
        });
    assert_rational_zero(
        "the integration policy retains only the signed vacuum remainder",
        &graph,
        &current,
        &(scheme + remainder),
        SpinReduction::SymbolicFirst,
    )?;

    // A nonzero exact commutator is a negative control: attaching gamma0 to
    // the other end of this open chain changes the spatial gamma components.
    let commutator =
        gamma_zero_at_port(&graph, &numerator, 1) - gamma_zero_at_port(&graph, &numerator, 0);
    assert!(
        !NumeratorJetCache::default()
            .canonicalize(&commutator)?
            .is_zero,
        "the physical QCD numerator must distinguish left and right projectors"
    );
    Ok(())
}

#[test]
fn fixed_boosted_os_matches_compact_operator_and_local_shell_projection() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    for components in [
        ["1", "0", "0", "0"],
        ["5/4", "3/4", "0", "0"],
        ["5/4", "1/4", "1/2", "1/2"],
    ] {
        let settings = UVgenerationSettings {
            generate_integrated: false,
            on_shell_subtraction_velocity: components.map(str::to_owned),
            ..Default::default()
        };
        let velocity = settings.on_shell_velocity()?;
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
        let actual = uv_limit(&Full4dCts::root(), &ctx, &current, &root, &current, &root)?;
        let (on_shell, remainder, _) = one_loop_reference(&graph, &current, &root, &velocity)?;
        assert_rational_zero(
            "fixed-velocity production operator equals independently differentiated compact form",
            &graph,
            &current,
            &(NumeratorInterface::export(actual.atom()).0 + &on_shell + &remainder),
            SpinReduction::SymbolicFirst,
        )?;
        let scheme = actual
            .active_sectors()
            .iter()
            .fold(Atom::Zero, |sum, sector| {
                sum + NumeratorInterface::export(sector.integration_atom()).0
            });
        assert_rational_zero(
            "fixed-velocity analytic payload is the signed vacuum remainder",
            &graph,
            &current,
            &(scheme + &remainder),
            SpinReduction::SymbolicFirst,
        )?;

        let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
        let complete = physical_expression(
            &graph,
            &current,
            &NumeratorInterface::export(
                &grow_factor(&ctx, &current, &root)?.try_add(actual.atom())?,
            )
            .0,
        );
        for sigma in [1, -1] {
            // Bind the live incoming momentum independently of the production
            // OS-point constructor and preserve all loop components.
            let evaluated = complete.replace_map(|view, _, out| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() == 2
                    && momentum.get(0) == Atom::num(0).as_view()
                {
                    **out = Atom::num(sigma)
                        * &mass
                        * Atom::add_many((0..4).map(|mu| {
                            &velocity[mu] * function!(GS.delta_vec, GS.cind(mu), momentum.get(1))
                        }));
                }
            });
            let projected = (&evaluated
                + Atom::num(sigma) * gamma_velocity_at_port(&graph, &evaluated, 1, &velocity))
                / Atom::num(2);
            assert_rational_zero(
                "matching boosted projector annihilates the complete local one-loop remainder",
                &graph,
                &current,
                &projected,
                SpinReduction::SymbolicFirst,
            )?;
        }
    }
    Ok(())
}

#[test]
fn production_os_massless_specialization_equals_soft_refined_degree_one() -> Result<()> {
    test_initialise()?;
    let (graph, mut current, root) = fixture(true)?;
    assert!(graph[EdgeIndex(2)].mass_atom(&graph.model).is_zero());
    for velocity in [
        ["1", "0", "0", "0"],
        ["5/4", "3/4", "0", "0"],
        ["5/4", "1/4", "1/2", "1/2"],
    ] {
        let settings = UVgenerationSettings {
            generate_integrated: false,
            on_shell_subtraction_velocity: velocity.map(str::to_owned),
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
        let seed = Full4dCts::root();
        current.spinney.renormalization_scheme = ApproximationType::OS;
        let on_shell = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
        current.spinney.renormalization_scheme = ApproximationType::IR;
        let soft = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
        assert_rational_zero(
            "OS at zero physical mass is H1 for every fixed velocity",
            &graph,
            &current,
            &NumeratorInterface::export(&on_shell.atom().try_sub(soft.atom())?).0,
            SpinReduction::SymbolicFirst,
        )?;
    }
    Ok(())
}

#[test]
fn production_os_removes_all_one_loop_ultraviolet_tensor_coefficients() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let ports = [0, 1].map(|hedge| {
        graph.underlying[Hedge(hedge)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom()
    });
    let scale = Atom::var(GS.rescale);
    for velocity in [
        ["1", "0", "0", "0"],
        ["5/4", "3/4", "0", "0"],
        ["5/4", "1/4", "1/2", "1/2"],
    ] {
        let settings = UVgenerationSettings {
            generate_integrated: false,
            on_shell_subtraction_velocity: velocity.map(str::to_owned),
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
        let seed = Full4dCts::root();
        let bare = grow(&seed, &ctx, &current, &root)?;
        let subtraction = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
        let coefficients = rational_coefficients(
            &graph,
            &current,
            &NumeratorInterface::export(&bare.try_add(subtraction.atom())?).0,
        )?;
        let mut nonzero_remainder = false;
        for row in 0..4 {
            for column in 0..4 {
                let mut jet = Atom::Zero;
                for (denominators, numerator) in &coefficients {
                    let component = numerator
                        * function!(GS.delta_vec, GS.cind(row), &ports[0])
                        * function!(GS.delta_vec, GS.cind(column), &ports[1]);
                    let scalar = scalar_component(
                        &(component * NumeratorInterface::export(&graph.global_prefactor.num).0),
                    )? / denominators
                        .iter()
                        .fold(Atom::one(), |product, factor| product * factor);
                    let scaled = scalar
                        .replace(GS.emr_mom(EdgeIndex(2), W_.x___))
                        .with(GS.emr_mom(EdgeIndex(2), W_.x___) / &scale)
                        / scale.pow(4);
                    jet += scaled.series_preserving_factors(
                        GS.rescale,
                        Atom::Zero.as_view(),
                        1,
                        &[],
                    )?;
                }
                let series = jet.series(GS.rescale, Atom::Zero, 1)?;
                for power in [-1, 0] {
                    assert!(
                        series.coefficient(power.into()).unwrap().expand().is_zero(),
                        "unsubtracted inverse-hard coefficient {power} in spinor component ({row},{column}), velocity={velocity:?}"
                    );
                }
                nonzero_remainder |= !series.coefficient(1.into()).unwrap().expand().is_zero();
            }
        }
        assert!(
            nonzero_remainder,
            "the UV check must not pass on an identically zero self-energy"
        );
    }
    Ok(())
}

#[test]
fn production_os_matches_equation_4_24_in_the_gluon_loop_basis() -> Result<()> {
    test_initialise()?;
    let (mut graph, mut current, root) = fixture(false)?;
    // Select a second basis on the untouched source graph, before constructing
    // a counterterm. No completed child or numerator is shifted afterward.
    current.spinney.lmb = graph
        .generate_loop_momentum_bases_of(current.subgraph())
        .into_iter()
        .find(|lmb| lmb.loop_edges.contains(&EdgeIndex(3)))
        .expect("the gluon edge defines a second one-loop basis");
    graph.loop_momentum_basis = current.spinney.lmb.clone();
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &graph.model, &settings);
    let source = graph
        .numerator(current.subgraph(), root.subgraph())?
        .to_d_dim(4)?
        .color_simplify()
        .unwrap()
        .get_single_atom()
        .unwrap();
    let numerator = physical_expression(&graph, &current, &source);
    let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
    let rest = |atom: &Atom, energy: Atom| {
        let replaced = atom.replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() > 1
                && usize::try_from(momentum.get(0))
                    .is_ok_and(|edge| current.lmb().ext_edges.contains(&EdgeIndex(edge)))
            {
                let point = GS.on_shell_momentum(
                    EdgeIndex(0),
                    &energy,
                    &UVgenerationSettings::default().on_shell_velocity().unwrap(),
                );
                let indices = momentum
                    .iter()
                    .skip(1)
                    .map(|index| index.to_owned())
                    .collect::<Vec<_>>();
                **output = GS.indexed_momentum(&point, &indices);
            }
        });
        GS.resolve_os_momenta(&replaced)
    };
    let massless = numerator.replace(mass.to_pattern()).with(Atom::Zero);
    let n_zero = rest(&massless, Atom::Zero);
    // Differentiate only the required linear external-momentum coefficient.
    // Keeping it as N(p,0)-N(0,0) would leave cancelling source factors inside
    // a tensor sum rather than exposing the independent Taylor coefficient.
    let external_scale = symbol!("os_gluon_basis_reference::external_scale");
    let scaled_numerator = massless.replace_map(|view, _, output| {
        if let AtomView::Fun(momentum) = view
            && momentum.get_symbol() == GS.emr_mom
            && momentum.get_nargs() > 1
            && usize::try_from(momentum.get(0))
                .is_ok_and(|edge| current.lmb().ext_edges.contains(&EdgeIndex(edge)))
        {
            **output = view.to_owned() * Atom::var(external_scale);
        }
    });
    let n_p = scaled_numerator
        .derivative(external_scale)
        .replace(external_scale)
        .with(Atom::Zero);
    assert!(
        !n_p.contains_symbol(Symbol::DERIVATIVE),
        "the linear reference coefficient is explicit"
    );
    let n_time = rest(&n_p, Atom::one());
    let vector_slot = Minkowski {}.new_rep(4).to_symbolic([]);
    let k = function!(GS.emr_mom, 3, &vector_slot);
    let p = function!(GS.emr_mom, 0, &vector_slot);
    let k_zero = GS.emr_mom(EdgeIndex(3), GS.cind(0));
    let k_squared = dot(&k, &k);
    let vacuum = GS.den(
        3,
        function!(GS.emr_mom, 3),
        Atom::var(GS.m_uv_vacuum).pow(2),
        &k_squared - Atom::var(GS.m_uv_vacuum).pow(2),
    );
    // The quark carries p-k in the fixture's directed gluon basis. Literal
    // differentiation of 1/[k²((p-k)²-m²)] gives +2 k_mu/D³. Together with
    // dN/dp this is (p-gamma0*pOS)^mu d_mu Sigma|p=m=0, Eq. 4.24.
    let compact = (&n_p - &mass * gamma_zero_at_port(&graph, &n_time, 1)) / vacuum.pow(2)
        + Atom::num(2)
            * (&n_zero * dot(&k, &p) - &mass * gamma_zero_at_port(&graph, &n_zero, 1) * &k_zero)
            / vacuum.pow(3);
    let gluon = GS.den(3, function!(GS.emr_mom, 3), Atom::Zero, k_squared.clone());
    let on_shell = [1, -1].into_iter().fold(Atom::Zero, |sum, sigma| {
        let n_sigma = rest(&numerator, Atom::num(sigma) * &mass);
        let projected =
            (&n_sigma + Atom::num(sigma) * gamma_zero_at_port(&graph, &n_sigma, 1)) / Atom::num(2);
        let quark = GS.den(
            2,
            function!(GS.emr_mom, 2),
            mass.pow(2),
            &k_squared - Atom::num(2 * sigma) * &mass * &k_zero,
        );
        sum + projected / (&gluon * quark)
    });
    let actual = uv_limit(&Full4dCts::root(), &ctx, &current, &root, &current, &root)?;
    assert_rational_zero(
        "composed OS equals the independently differentiated Eq. 4.24 in the gluon basis",
        &graph,
        &current,
        &(NumeratorInterface::export(actual.atom()).0 + on_shell + compact),
        SpinReduction::SymbolicFirst,
    )
}

#[test]
fn production_os_dirac_cff_matches_independent_energy_residues() -> Result<()> {
    use crate::{
        graph::cuts::CutSet,
        settings::global::OrientationPattern,
        uv::approx::{
            OrientationProjection, local_3d::Localizer, projected_4d::Local4dProjectionContext,
        },
    };

    test_initialise()?;
    for (velocity_components, samples) in [
        (
            ["1", "0", "0", "0"],
            [([9, 12, 0], [15, 15]), ([9, -12, 0], [15, 15])],
        ),
        (
            ["5/4", "3/4", "0", "0"],
            [([15, 0, 0], [9, 21]), ([-15, 0, 0], [21, 9])],
        ),
        (
            ["5/4", "1/4", "1/2", "1/2"],
            [([5, 10, 10], [9, 21]), ([-5, -10, -10], [21, 9])],
        ),
    ] {
        let (mut graph, current, root) = fixture(false)?;
        let settings = UVgenerationSettings {
            generate_integrated: false,
            on_shell_subtraction_velocity: velocity_components.map(str::to_owned),
            ..Default::default()
        };
        let velocity = settings.on_shell_velocity()?;
        let local = uv_limit(
            &Full4dCts::root(),
            &UVCtx::new(&graph, &graph.model, &settings),
            &current,
            &root,
            &current,
            &root,
        )?;
        let (_, _, numerator) = one_loop_reference(&graph, &current, &root, &velocity)?;
        let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
        let cutset = CutSet::empty(graph.n_hedges());
        let pattern = OrientationPattern::default();
        let production = Default::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact(&production, &options, &pattern, false),
        );
        let mut context = Local4dProjectionContext::default();
        let mut actual = Atom::Zero;
        for sector in local.projection_sectors()? {
            // These production calls require independent exact numerator and
            // denominator reconstruction certificates before generating CFF.
            let canonical = sector.canonical_projection(&graph)?;
            let projected = localizer
                .project_factorized_taylor_sector(&mut graph, &canonical, &mut context)?
                .resolved()?;
            for (_, expression) in projected.iter() {
                actual += expression;
            }
        }
        assert!(!actual.contains_symbol(GS.os_mom));
        assert!(!actual.contains_symbol(symbol!("gammalooprs::uv::numerator_family")));

        let z_symbol = symbol!("os_dirac_residue_reference::z");
        let z = Atom::var(z_symbol);
        let at_energy = |value: &Atom, energy: &Atom| {
            value.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() == 2
                    && momentum.get(0) == Atom::num(2).as_view()
                {
                    **output = function!(GS.emr_vec, 2, momentum.get(1))
                        + energy * GS.energy_delta(momentum.get(1));
                }
            })
        };
        // Independent positive-energy poles in the quark-k routing are
        // k0=E_m and sigma*m*v0+|k-sigma*m*v|. These exact rational samples
        // have distinct poles and rational square roots in every branch.
        // They supplement the symbolic reconstruction certificates above.
        let fermion_energy = Atom::num(17);
        let vacuum_energy = Atom::num(25);
        let n_zero = at_energy(&numerator.replace(mass.to_pattern()).with(Atom::Zero), &z);
        let k_dot_p = (1..=3).fold(
            &z * GS.emr_mom(EdgeIndex(0), GS.cind(0)),
            |sum, component| {
                sum - GS.emr_mom(EdgeIndex(2), GS.cind(component))
                    * GS.emr_mom(EdgeIndex(0), GS.cind(component))
            },
        );
        let k_dot_velocity = (1..=3).fold(&z * &velocity[0], |sum, component| {
            sum - GS.emr_mom(EdgeIndex(2), GS.cind(component)) * &velocity[component]
        });
        let vacuum_numerator = Atom::num(2)
            * (&n_zero * k_dot_p
                - &mass * gamma_velocity_at_port(&graph, &n_zero, 1, &velocity) * k_dot_velocity);
        let at_vacuum_pole = |value: Atom| value.replace(z_symbol).with(vacuum_energy.clone());
        let twice_energy = Atom::num(2) * &vacuum_energy;
        // Res A(z)/(z²-E²)^3 = 1/2 d²[A(z)/(z+E)^3] at z=E.
        let vacuum_residue =
            at_vacuum_pole(vacuum_numerator.derivative(z_symbol).derivative(z_symbol))
                / (Atom::num(2) * twice_energy.pow(3))
                - Atom::num(3) * at_vacuum_pole(vacuum_numerator.derivative(z_symbol))
                    / twice_energy.pow(4)
                + Atom::num(6) * at_vacuum_pole(vacuum_numerator) / twice_energy.pow(5);
        for (spatial, gluon_energies) in samples {
            assert_eq!(spatial.iter().map(|value| value * value).sum::<i32>(), 225);
            let mut residue = vacuum_residue.clone();
            for (sigma, gluon_energy) in [1, -1].into_iter().zip(gluon_energies.map(Atom::num)) {
                let gluon_spatial_squared = spatial
                    .iter()
                    .enumerate()
                    .map(|(index, value)| {
                        (Atom::num(*value) - Atom::num(8 * sigma) * &velocity[index + 1]).pow(2)
                    })
                    .sum::<Atom>();
                assert_eq!(gluon_energy.pow(2), gluon_spatial_squared);
                let shift = Atom::num(sigma) * &mass * &velocity[0];
                let quark_value = at_energy(&numerator, &fermion_energy);
                let gluon_value = at_energy(&numerator, &(&shift + &gluon_energy));
                let project = |value: Atom| {
                    (&value
                        + Atom::num(sigma) * gamma_velocity_at_port(&graph, &value, 1, &velocity))
                        / Atom::num(2)
                };
                residue += project(quark_value)
                    / (Atom::num(2)
                        * &fermion_energy
                        * ((&fermion_energy - &shift).pow(2) - gluon_energy.pow(2)))
                    + project(gluon_value)
                        / (Atom::num(2)
                            * &gluon_energy
                            * ((&shift + &gluon_energy).pow(2) - fermion_energy.pow(2)));
            }
            // Clockwise contour and the single local forest sign cancel. The
            // remaining physical loop measure is i/(2*pi)^3.
            assert!(
                !residue.contains_symbol(Symbol::DERIVATIVE),
                "all boosted residue derivatives are explicit"
            );
            let expected = Atom::i() * residue / (Atom::num(8) * Atom::var(GS.pi).pow(3));
            let specialize = |value: &Atom| {
                let mut value = GS
                    .erase_uv_momentum_provenance(value)
                    .replace(GS.dim)
                    .with(Atom::num(4))
                    .replace(crate::utils::ose_atom_from_index(EdgeIndex(2)))
                    .with(fermion_energy.clone())
                    .replace(mass.to_pattern())
                    .with(Atom::num(8))
                    .replace(GS.m_uv_vacuum)
                    .with(Atom::num(20))
                    .replace(GS.m_uv_expansion)
                    .with(Atom::num(20));
                // Every abstract vector component is materialized only by the
                // finite tensor contraction below. Bind spatial entries through
                // a represented delta basis, preserving the graph factors.
                value = value.replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && (momentum.get_symbol() == GS.emr_mom
                            || momentum.get_symbol() == GS.emr_vec)
                        && momentum.get_nargs() == 2
                        && momentum.get(0) == Atom::num(2).as_view()
                    {
                        let index = momentum.get(1);
                        if let AtomView::Fun(concrete) = index
                            && concrete.get_symbol()
                                == spenso::structure::abstract_index::AIND_SYMBOLS.cind
                            && concrete.get_nargs() == 1
                            && let Ok(component @ 1..=3) = usize::try_from(concrete.get(0))
                        {
                            **output = Atom::num(spatial[component - 1]);
                        } else if momentum.get_symbol() == GS.emr_vec {
                            **output = (1..=3).fold(Atom::Zero, |sum, component| {
                                sum + Atom::num(spatial[component - 1])
                                    * function!(GS.delta_vec, GS.cind(component), index)
                            });
                        }
                    }
                });
                value
            };
            let actual = specialize(&actual);
            let expected = specialize(&expected);
            debug_tags!(#uv, #local, #on_shell, #inspect;
                stage = "os_oracle_residue_interfaces",
                velocity = ?velocity_components, spatial = ?spatial,
                file.actual = %actual.to_canonical_string(),
                file.expected = %expected.to_canonical_string(),
                "Independent factorized residue expressions before finite tensor contraction"
            );
            let difference = actual - expected;
            let certificate =
                NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
            let is_zero = certificate.is_zero;
            if !is_zero {
                debug_tags!(#uv, #local, #four_d, #cff, #trace;
                    stage = "os_dirac_residue_oracle_failure",
                    velocity = ?velocity_components, spatial = ?spatial,
                    file.difference = %difference,
                    file.contracted_tensor = %certificate.key,
                    "Nonzero factorized OS residue-oracle difference"
                );
            }
            assert!(
                is_zero,
                "physical Dirac CFF differs from independent residues at velocity {velocity_components:?}, spatial momentum {spatial:?}"
            );
        }
    }
    Ok(())
}

#[test]
fn on_shell_taylor_preserves_fixed_point_denominator_contractions() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &graph.model, &settings);
    let source = graph
        .numerator(current.subgraph(), &graph.empty_subgraph())?
        .state
        .expr;
    let vector_slot = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
    let edge = current.lmb().loop_edges[LoopIndex(0)];
    let momentum = function!(GS.emr_mom, usize::from(edge), &vector_slot);
    let time = GS.energy_delta(&vector_slot);
    let compact_time_component = function!(GS.delta_vec, GS.cind(0), &momentum);
    let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
    let prescription = Atom::var(symbol!("os_taylor_prescription"; Scalar));
    let opaque_power = function!(symbol!("os_taylor_spectator"; Scalar), 7).pow(3);
    let opaque_inverse = (Atom::var(symbol!("os_taylor_x"; Scalar)) + Atom::one()).pow(-2);
    let opaque = &opaque_power * &opaque_inverse;
    for sign in [-1, 1] {
        let polynomial = |time_component: Atom| {
            dot(&momentum, &momentum)
                + Atom::num(2 * sign) * &mass * time_component
                + mass.pow(2)
                + &prescription
        };
        let operand = |time_component| {
            // One denominator keeps the prescription inside the retained UV
            // order, so the check below is sensitive to dropping it.
            source.with_checked_expression(
                source.expression() * &opaque
                    / GS.den(
                        usize::from(edge),
                        function!(GS.emr_mom, usize::from(edge)),
                        mass.pow(2),
                        polynomial(time_component),
                    ),
            )
        };
        let compact = operand(compact_time_component.clone())?;
        let explicit = operand(dot(&momentum, &time))?;
        let actual = t_raw(&compact, &ctx, &current, &root, current.lmb())?;
        let expected = t_raw(&explicit, &ctx, &current, &root, current.lmb())?;
        assert_eq!(actual.expression(), expected.expression());
        assert_eq!(
            actual.structure().logical_slots(),
            source.structure().logical_slots()
        );
        assert!(!actual.is_zero());
        assert!(actual.expression().contains(&opaque_power));
        assert!(actual.expression().contains(&opaque_inverse));
        assert!(
            actual
                .expression()
                .contains_symbol(prescription.get_symbol().unwrap())
        );
    }
    Ok(())
}

#[test]
fn on_shell_rest_frame_basis_materializes_scalar_inner_products() -> Result<()> {
    test_initialise()?;
    let vector_slot = Minkowski {}.new_rep(4).to_symbolic([]);
    let time = GS.energy_delta(&vector_slot);
    let space = function!(GS.delta_vec, GS.cind(1), &vector_slot);
    let momentum = function!(GS.emr_mom, 0, &vector_slot);
    assert_eq!(scalar_component(&dot(&time, &time))?, Atom::one());
    assert_eq!(scalar_component(&dot(&time, &space))?, Atom::Zero);
    assert_eq!(scalar_component(&dot(&space, &space))?, Atom::num(-1));
    assert_eq!(
        scalar_component(&dot(&time, &momentum))?,
        GS.emr_mom(EdgeIndex(0), GS.cind(0))
    );
    assert_eq!(
        scalar_component(&dot(&space, &momentum))?,
        -GS.emr_mom(EdgeIndex(0), GS.cind(1))
    );
    Ok(())
}

#[test]
fn actual_native_charge_mates_have_the_same_os_component_conventions() -> Result<()> {
    use crate::uv::approx::integrated::{self, Integrated};
    use idenso::color::ColorSimplifier;
    use linnet::half_edge::subgraph::ModifySubSet;
    use spenso::{network::library::symbolic::ETS, shadowing::TensorCollectFilter};

    test_initialise()?;
    // These are the complete native acceptance graphs, including their exact
    // source slots, local Feynman factors, physical LMB and global metadata.
    // Their opposite arrows represent the charge-conjugate insertions; both
    // components are positive-top matrices in their own flow-oriented chart.
    for (source, name, row, column, incoming_edge) in [
        (
            include_str!(
                "../../../../../tests/resources/graphs/local_os_native_gl1_self_energy.dot"
            ),
            "GL1",
            5,
            8,
            4,
        ),
        (
            include_str!(
                "../../../../../tests/resources/graphs/local_os_native_gl2_self_energy.dot"
            ),
            "GL2",
            7,
            6,
            3,
        ),
    ] {
        let mut graphs = Graph::import_dot(source, &crate::utils::load_generic_model("sm"))?;
        assert_eq!(graphs.len(), 1);
        let graph = graphs.remove(0);
        let mut selected = graph.empty_subgraph::<SuBitGraph>();
        for edge in [1, 2] {
            selected.add(graph[&EdgeIndex(edge)].1);
        }
        let current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(selected, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::OS,
                1,
            )?
            .expect("the actual top self-energy is UV divergent"),
            topo_order: 1,
        };
        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        assert_eq!(
            graph
                .dummy_less_full_crown(current.subgraph())
                .included_iter()
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([Hedge(row), Hedge(column)]),
            "{name}: the native component retains its two physical fermion ports"
        );
        assert_eq!(graph.flow(Hedge(row)), Flow::Source);
        assert_eq!(graph.flow(Hedge(column)), Flow::Sink);
        assert_eq!(graph[&Hedge(column)], EdgeIndex(incoming_edge));
        assert_eq!(
            current.lmb().loop_edges.iter().copied().collect::<Vec<_>>(),
            vec![EdgeIndex(1)],
            "{name}: this independent reference uses the actual fermion loop carrier",
        );
        // The incoming half-edge owns p_OS, while the compatible component
        // chart chooses its independent external carrier. Keep this exact
        // routing certificate separate from ownership for each native graph.
        let incoming_signature = &current.lmb().edge_signatures[EdgeIndex(incoming_edge)];
        assert!(
            incoming_signature
                .internal
                .iter()
                .all(|sign| *sign == crate::momentum::SignOrZero::Zero)
        );
        let independent_external = incoming_signature
            .external
            .iter_enumerated()
            .filter(|(_, sign)| sign.is_sign())
            .collect::<Vec<_>>();
        assert_eq!(independent_external.len(), 1);
        let (external_index, sign) = independent_external[0];
        assert_eq!(*sign, crate::momentum::SignOrZero::Plus);
        let external_coordinate = usize::from(current.lmb().ext_edges[external_index]);
        assert_eq!(
            current.lmb().edge_signatures[EdgeIndex(external_coordinate)],
            current.lmb().edge_signatures[EdgeIndex(incoming_edge)],
            "{name}: the reference external carrier equals the physical incoming momentum"
        );
        assert_eq!(
            current.lmb().edge_signatures[EdgeIndex(3)],
            current.lmb().edge_signatures[EdgeIndex(4)],
            "{name}: the two adjacent top lines carry the same signed momentum"
        );
        let mass = graph[EdgeIndex(1)].mass_atom(&graph.model);
        let numerator = graph
            .numerator(current.subgraph(), root.subgraph())?
            .to_d_dim(GS.dim)?
            .color_simplify()?
            .state
            .expr;
        let row_slot = graph.underlying[Hedge(row)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let column_slot = graph.underlying[Hedge(column)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let identity = function!(ETS.metric, &row_slot, &column_slot);
        // Amputation leaves the VERTEX color interface, dual to the removed
        // external propagator. The edge-to-edge Kronecker helper instead has
        // the opposite variance and cannot supply this component identity.
        let color_slots = [row, column].map(|hedge| {
            graph.underlying[Hedge(hedge)]
                .num_indices
                .color_indices
                .vertex_indices
                .external_structure_iter()
                .next()
                .expect("the actual top boundary has one vertex color slot")
                .to_atom()
        });
        let color = function!(ETS.metric, &color_slots[0], &color_slots[1]);
        let scalar_coefficient =
            Atom::num((4, 3)) * Atom::from(crate::model::UFOSymbol::from("GC_11")).pow(2);
        let coefficient = &scalar_coefficient * color;
        let lorentz = spenso::mink!(GS.dim, Atom::from(Aind::new_dummy()));
        let slash_q =
            function!(GS.emr_mom, 1, &lorentz) * idenso::gamma!(&row_slot, &column_slot, &lorentz);
        let expected_numerator = &coefficient
            * ((Atom::num(2) - Atom::var(GS.dim)) * &slash_q
                + Atom::var(GS.dim) * &mass * &identity);
        let actual_tensor = &numerator;
        let expected_tensor =
            NumeratorInterface::admit(expected_numerator.clone(), None)?.try_into_indexed()?;
        let interfaces = [actual_tensor, &expected_tensor].map(|tensor| {
            tensor
                .structure()
                .canonical()
                .external_structure_iter()
                .collect::<Vec<_>>()
        });
        debug_tags!(#uv, #local, #on_shell, #inspect;
            stage = "native_os_component_bare_interfaces",
            graph = name,
            actual_interface = ?interfaces[0], expected_interface = ?interfaces[1],
            file.actual = %NumeratorInterface::export(&numerator).0.to_canonical_string(),
            file.expected = %expected_numerator.to_canonical_string(),
            "Actual native component and independent bare reference before subtraction"
        );
        assert_eq!(
            interfaces[0], interfaces[1],
            "{name}: the independent reference must preserve the amputated spin/color interface"
        );
        let difference = numerator.try_sub(&expected_tensor)?;
        let difference = difference
            .with_checked_expression(difference.expression().to_cof_dimension_invariants())?;
        let normalized_tensor = integrated::simplify(&difference).wrap_err_with(|| {
            format!("{name}: independent d-dimensional bare numerator comparison")
        })?;
        // Canonical Clifford reduction may retain a factored scalar expression
        // such as -d+(d-2)+2 around a complete color/spin word. Certify only
        // those scalar coefficients, keeping the tensor word and d untouched.
        let normalized_tensor: SymbolicTensor<spenso::structure::partial::PartialStructure> =
            normalized_tensor.into();
        for (_, coefficient) in
            normalized_tensor.coefficient_list(TensorCollectFilter::<0>::TaggedTensors)?
        {
            assert!(coefficient.is_scalar());
            assert!(
                coefficient.expression().collect_compact_factors().is_zero(),
                "{name}: independent d-dimensional bare phase/color/spin reference"
            );
        }

        for components in [["1", "0", "0", "0"], ["5/4", "1/4", "1/2", "1/2"]] {
            let settings = UVgenerationSettings {
                add_marker: false,
                generate_integrated: true,
                on_shell_subtraction_velocity: components.map(str::to_owned),
                ..Default::default()
            };
            let velocity = settings.on_shell_velocity()?;
            let ctx = UVCtx::new(&graph, &graph.model, &settings);
            let local = uv_limit(&Full4dCts::root(), &ctx, &current, &root, &current, &root)?;
            let numerator4 =
                physical_expression(&graph, &current, &NumeratorInterface::export(&numerator).0);
            let massless = numerator4.replace(mass.to_pattern()).with(0);
            let vector_slot = Minkowski {}.new_rep(4).to_symbolic([]);
            let q = function!(GS.emr_mom, 1, &vector_slot);
            let p = function!(GS.emr_mom, external_coordinate, &vector_slot);
            let q_squared = dot(&q, &q);
            let qv: Atom = (0..4)
                .map(|mu| {
                    Atom::num(if mu == 0 { 1 } else { -1 })
                        * &velocity[mu]
                        * GS.emr_mom(EdgeIndex(1), GS.cind(mu))
                })
                .sum();
            let quark = GS.den(
                1,
                function!(GS.emr_mom, 1),
                mass.pow(2),
                &q_squared - mass.pow(2),
            );
            let vacuum = GS.den(
                1,
                function!(GS.emr_mom, 1),
                Atom::var(GS.m_uv_vacuum).pow(2),
                &q_squared - Atom::var(GS.m_uv_vacuum).pow(2),
            );
            let mut on_shell = Atom::Zero;
            for sigma in [1, -1] {
                let point = GS.on_shell_momentum(
                    EdgeIndex(incoming_edge),
                    &(Atom::num(sigma) * &mass),
                    &velocity,
                );
                let gluon = GS.den(
                    2,
                    function!(GS.emr_mom, 1) - &point,
                    Atom::Zero,
                    &q_squared - Atom::num(2 * sigma) * &mass * &qv + mass.pow(2),
                );
                on_shell += (&numerator4
                    + Atom::num(sigma)
                        * gamma_velocity_at_port(&graph, &numerator4, row, &velocity))
                    / (Atom::num(2) * &quark * gluon);
            }
            let remainder = Atom::num(2)
                * (&massless * dot(&q, &p)
                    - &mass * gamma_velocity_at_port(&graph, &massless, row, &velocity) * &qv)
                / vacuum.pow(3);
            assert_rational_zero(
                &format!("{name}: actual native signed local OS"),
                &graph,
                &current,
                &(NumeratorInterface::export(local.atom()).0 + &on_shell + &remainder),
                SpinReduction::SymbolicFirst,
            )?;
            let scheme = local
                .active_sectors()
                .iter()
                .fold(Atom::Zero, |sum, sector| {
                    sum + NumeratorInterface::export(sector.integration_atom()).0
                });
            assert_rational_zero(
                &format!("{name}: actual native signed vacuum input"),
                &graph,
                &current,
                &(scheme + &remainder),
                SpinReduction::SymbolicFirst,
            )?;

            // Independent exact nonzero control: the two ends of the actual
            // source cannot be exchanged merely because the graph is a C mate.
            let commutator = gamma_velocity_at_port(&graph, &numerator4, row, &velocity)
                - gamma_velocity_at_port(&graph, &numerator4, column, &velocity);
            assert!(
                !NumeratorJetCache::default()
                    .canonicalize(&commutator.undo_dots::<Aind>()?)?
                    .is_zero
            );

            let vakint_settings = settings.vakint.true_settings();
            let integrated = Integrated::new(crate::utils::vakint()?, &vakint_settings)
                .run(&local, &ctx, &current, &root, &current, &root)?;
            let field = (&coefficient
                * (function!(GS.emr_mom, external_coordinate, &lorentz)
                    * idenso::gamma!(&row_slot, &column_slot, &lorentz)
                    - &mass * &identity))
                .replace(GS.dim)
                .with(4);
            let normalization = Atom::i() / (Atom::num(16) * Atom::var(Symbol::PI).pow(2));
            let epsilon = Atom::var(GS.dim_epsilon);
            let pole = NumeratorInterface::export(&integrated.physical_pole()?)
                .0
                .replace(GS.dim)
                .with(4);
            assert_rational_zero(
                &format!("{name}: actual native field pole"),
                &graph,
                &current,
                &(pole - &normalization * &field / &epsilon),
                SpinReduction::SymbolicFirst,
            )?;
            // Finite epsilon^0 only; no enclosing pole in this physical NLO component.
            let finite = integrated.physical_finite_counterterm()?;
            let finite = finite.with_checked_expression(
                finite
                    .expression()
                    .replace(GS.dim)
                    .with(4)
                    .replace(GS.dim_epsilon)
                    .with(0)
                    .replace(GS.m_uv_vacuum)
                    .with(1),
            )?;
            let expected_finite = &normalization
                * &field
                * (Atom::one() + function!(Symbol::LOG, Atom::one() / Atom::var(GS.mu_r_sq)));
            assert_rational_zero(
                &format!("{name}: actual native finite scheme addback"),
                &graph,
                &current,
                &(NumeratorInterface::export(&finite).0 - expected_finite),
                SpinReduction::SymbolicFirst,
            )?;

            // Check the actual two adjacent propagator NUMERATORS before
            // energy projection. N=i*(p-slash+m) obeys N*F*N=i*D_phys*N,
            // with F=p-slash-m and D_phys=p^2-m^2. A propagator denominator
            // carrying nonzero i0 cannot be cancelled by this identity as an
            // ordinary rational function; no prescription is dropped here.
            let outgoing_edge = graph[&Hedge(row)];
            let mut adjacent = graph.empty_subgraph::<SuBitGraph>();
            for edge in [outgoing_edge, EdgeIndex(incoming_edge)] {
                adjacent.add(graph[&edge].1);
            }
            let adjacent_numerator = graph[outgoing_edge]
                .num
                .value
                .clone()
                .multiply(graph[EdgeIndex(incoming_edge)].num.value.clone())?;
            let adjacent_numerator = adjacent_numerator.with_checked_expression(
                adjacent_numerator
                    .expression()
                    .replace_multiple(graph.uv_wrapped_replacement(
                        &adjacent,
                        current.lmb(),
                        &[W_.x___],
                    )),
            )?;
            let outer_hedges = [graph.inv(Hedge(row)), graph.inv(Hedge(column))];
            let outer_spin = outer_hedges.map(|hedge| {
                graph.underlying[hedge]
                    .edge_spin_slots()
                    .next()
                    .unwrap()
                    .to_atom()
            });
            let outer_color = outer_hedges.map(|hedge| {
                graph.underlying[hedge]
                    .edge_color_slots()
                    .next()
                    .unwrap()
                    .to_atom()
            });
            let pinch_lorentz = spenso::mink!(4, Atom::from(Aind::new_dummy()));
            let pinched_numerator = Atom::i()
                * (function!(GS.emr_mom, external_coordinate, &pinch_lorentz)
                    * idenso::gamma!(&outer_spin[0], &outer_spin[1], &pinch_lorentz)
                    + &mass * function!(ETS.metric, &outer_spin[0], &outer_spin[1]))
                * function!(ETS.metric, &outer_color[0], &outer_color[1]);
            let physical_polynomial = dot(&p, &p) - mass.pow(2);
            let adjacent_denominators = graph.denominator(&adjacent, &graph.model, |_| -1);
            let terms = FourDTerm::from_view(adjacent_denominators.as_view())?;
            assert_eq!(terms.len(), 1);
            assert_eq!(terms[0].denominators.len(), 2);
            assert_eq!(
                terms[0]
                    .denominators
                    .iter()
                    .map(|denominator| denominator.source_edge)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([outgoing_edge, EdgeIndex(incoming_edge)])
            );
            for denominator in &terms[0].denominators {
                assert_eq!(denominator.mass_squared, mass.pow(2));
                let routed_polynomial =
                    denominator
                        .full_expr
                        .replace_multiple(graph.uv_wrapped_replacement(
                            &adjacent,
                            current.lmb(),
                            &[W_.x___],
                        ));
                assert_rational_zero(
                    &format!("{name}: repeated physical top denominator polynomial"),
                    &graph,
                    &current,
                    &(routed_polynomial - &physical_polynomial),
                    SpinReduction::SymbolicFirst,
                )?;
            }
            let expected_pinched =
                Atom::i() * &scalar_coefficient * &physical_polynomial * &pinched_numerator;
            let field_tensor =
                NumeratorInterface::admit(field.clone(), None)?.try_into_indexed()?;
            assert_rational_zero(
                &format!("{name}: native open spin/color field insertion numerator"),
                &graph,
                &current,
                &(NumeratorInterface::export(&adjacent_numerator.clone().multiply(field_tensor)?)
                    .0
                    - &expected_pinched),
                SpinReduction::SymbolicFirst,
            )?;
            assert!(
                !NumeratorJetCache::default()
                    .canonicalize(&expected_pinched.undo_dots::<Aind>()?)?
                    .is_zero,
                "{name}: the native numerator identity is nontrivial off shell"
            );
            // The integrated coefficient supplies another i. Thus its native
            // pinched coefficient is -C/(16*pi^2)*(1-L), not i*C/(16*pi^2).
            let pinched_finite_coefficient = -&scalar_coefficient
                / (Atom::num(16) * Atom::var(Symbol::PI).pow(2))
                * (Atom::one() + function!(Symbol::LOG, Atom::one() / Atom::var(GS.mu_r_sq)));
            assert_rational_zero(
                &format!("{name}: actual finite insertion has the native pinched phase"),
                &graph,
                &current,
                &(NumeratorInterface::export(&adjacent_numerator.multiply(finite.clone())?).0
                    - pinched_finite_coefficient * &physical_polynomial * &pinched_numerator),
                SpinReduction::SymbolicFirst,
            )?;

            if components != ["1", "0", "0", "0"] {
                // The finite matrix's boost independence was checked above;
                // its physical cograph needs one exact projection comparison.
                continue;
            }
            let source_scope = graph.full_filter().subtract(&graph.initial_state_cut);
            let remaining = source_scope.subtract(current.subgraph());
            // Keep both photon vertices and the other top propagator factored.
            // The fixed external polarization projector multiplies both sides
            // unchanged; checking all open photon components is stronger here.
            let overall_factor = SymbolicTensor::infer_indexed(graph.overall_factor.clone())?;
            assert!(overall_factor.is_scalar());
            let spectator = graph
                .numerator(&remaining.subtract(&adjacent), current.subgraph())?
                .state
                .expr
                .multiply(graph.global_prefactor.num.clone())?
                .multiply(overall_factor)?;
            let pinched_prefix = current
                .subgraph()
                .union(&graph.get_edge_subgraph(EdgeIndex(incoming_edge)));
            let pinched_owners = source_scope.subtract(&pinched_prefix);
            let pinched_finite_coefficient = -&scalar_coefficient
                / (Atom::num(16) * Atom::var(Symbol::PI).pow(2))
                * (Atom::one() + function!(Symbol::LOG, Atom::one() / Atom::var(GS.mu_r_sq)));
            let pinched_source_numerator = pinched_numerator
                .replace(function!(GS.emr_mom, external_coordinate, W_.x___))
                .with(function!(GS.emr_mom, usize::from(outgoing_edge), W_.x___).to_pattern());
            let pinched_source =
                NumeratorInterface::admit(pinched_source_numerator, None)?.try_into_indexed()?;
            let pinched_source = pinched_source.with_checked_expression(
                pinched_source.expression() * &pinched_finite_coefficient,
            )?;
            let inputs = [
                (
                    "raised",
                    // Keep the actual distinct edge owners until exact source
                    // reconstruction, including their raised-pole samples.
                    spectator
                        .clone()
                        .multiply(graph[outgoing_edge].num.value.clone())?
                        .multiply(finite)?
                        .multiply(graph[EdgeIndex(incoming_edge)].num.value.clone())?,
                    remaining,
                    current.subgraph().clone(),
                    3,
                    2,
                ),
                (
                    "pinched",
                    spectator.multiply(pinched_source)?,
                    pinched_owners,
                    pinched_prefix,
                    2,
                    1,
                ),
            ];
            let crown = graph.dummy_stripped_external_flows_of(&source_scope);
            let source_edges = graph.paired_edges(&source_scope);
            // The pinched edge cannot remain a parent loop carrier. Choose
            // the equal-momentum surviving top line in the common parent
            // chart, then certify the signed routing of every active edge.
            let parent_carriers = graph
                .get_edge_subgraph(EdgeIndex(1))
                .union(&graph.get_edge_subgraph(outgoing_edge));
            let parent_lmb = graph.lmb_impl(
                &source_scope,
                &source_scope.subtract(&parent_carriers),
                crown.clone(),
            )?;
            assert_eq!(
                parent_lmb.loop_edges.iter().copied().collect::<Vec<_>>(),
                vec![EdgeIndex(1), outgoing_edge]
            );
            for edge in &source_edges {
                assert_eq!(
                    parent_lmb.edge_signatures[*edge],
                    graph.loop_momentum_basis.edge_signatures[*edge],
                    "{name}: the pinched carrier change preserves the exact physical chart"
                );
            }
            let cutset = crate::graph::cuts::CutSet::empty(graph.n_hedges());
            let options = graph.denominator_only_cff_3d_expression_options();
            let mut projection_graph = graph.clone();
            let mut projected = Vec::new();
            for (label, numerator, owners, prefix, denominator_count, maximum_power) in inputs {
                let prefix = InternalSubGraph::cleaned_filter_optimist(prefix, graph.as_ref());
                let quotient = graph.shrunken_sub_lmb(
                    &source_scope,
                    &prefix,
                    crown.clone(),
                    Some(&parent_lmb),
                )?;
                assert_eq!(quotient.loop_edges.len(), 1);
                let expression = numerator.with_checked_expression(
                    numerator.expression() * graph.denominator(&owners, &graph.model, |_| -1),
                )?;
                let sector = FourDSector::new(
                    expression.clone(),
                    vec![(owners, source_scope.clone(), quotient.clone())],
                    Vec::new(),
                );
                // This first certificate preserves the complete numerator and
                // denominator momenta, masses, powers and source ownership.
                let canonical = sector.canonical_projection(&graph)?;
                assert_eq!(canonical.terms.len(), 1);
                assert_eq!(canonical.terms[0].source_witness.len(), denominator_count);
                assert_eq!(
                    canonical.terms[0].powers.values().copied().max(),
                    Some(maximum_power)
                );
                let terms = FourDTerm::from_tensor(&expression)?;
                assert_eq!(terms.len(), 1);
                // A live photon crown is a physical source, not a Taylor
                // vacuum. The existing source owner certifies each chosen
                // numerator assignment before generating its causal CFF.
                let (cff, _) = projection_graph.cff_from_4d_denominators_in_uv_sub_lmb(
                    &terms[0].denominators,
                    source_edges.iter().copied(),
                    crown.included_iter(),
                    &quotient,
                    crate::graph::ExactUvSubLmbFrame::RetainedPhysicalCrown,
                    &cutset,
                    &options,
                    terms[0].numerator.expression(),
                    None,
                    &[],
                    &BTreeMap::new(),
                )?;
                let mut value = terms[0].numerator.with_checked_expression(Atom::Zero)?;
                let mut source_maps = 0;
                for term in cff.terms.values() {
                    for orientation in &term.orientations {
                        source_maps += 1;
                        let mapped = terms[0].numerator.with_checked_expression(
                            term.map_exact_source_numerator(&orientation.orientation, None)?,
                        )?;
                        let contribution = mapped.with_checked_expression(
                            &orientation.expression * mapped.expression(),
                        )?;
                        value = value.try_add(&contribution)?;
                    }
                }
                assert!(
                    source_maps > 0,
                    "{name}/{label}: the complete source-map sum is nonempty"
                );
                value = value.with_checked_expression(
                    value.expression() * Atom::num(cff.production_prefactor_factor()),
                )?;
                assert!(
                    !value
                        .expression()
                        .contains_symbol(GS.renormalization_localization_scale)
                );
                debug_tags!(#uv, #local, #on_shell, #cff, #inspect;
                    stage = "native_os_pinched_cograph_projection",
                    graph = name, source = label, source_maps,
                    file.value = %NumeratorInterface::export(&value).0.to_canonical_string(),
                    "Certified physical cograph projection with its frozen kernel kept outside"
                );
                projected.push(value);
            }
            // CFF represents the shared causal prescription's boundary value.
            // These exact rational component checks supplement the general
            // reconstruction certificates; they are not a rational identity
            // obtained by cancelling a formal nonzero i0 in a double pole.
            for spatial in [[4, 0, 0], [0, 4, 0]] {
                let values = projected
                    .iter()
                    .map(|expression| {
                        let mut value = NumeratorInterface::export(expression)
                            .0
                            .replace(mass.to_pattern())
                            .with(3)
                            .replace(GS.mu_r_sq)
                            .with(1)
                            .replace(GS.dim)
                            .with(4);
                        for edge in [3, 4, 5] {
                            value = value
                                .replace(crate::utils::ose_atom_from_index(EdgeIndex(edge)))
                                .with(5);
                        }
                        value.replace_map(|view, _, output| {
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
                            if ![0, 3, 4, 5].contains(&edge) {
                                return;
                            }
                            let index = momentum.get(1);
                            if edge == 0 {
                                **output = if momentum.get_symbol() == GS.emr_vec {
                                    Atom::Zero
                                } else {
                                    GS.energy_delta(index)
                                };
                            } else if let AtomView::Fun(component) = index
                                && component.get_symbol()
                                    == spenso::structure::abstract_index::AIND_SYMBOLS.cind
                                && let Ok(component @ 1..=3) = usize::try_from(component.get(0))
                            {
                                **output = Atom::num(spatial[component - 1]);
                            } else if momentum.get_symbol() == GS.emr_vec {
                                **output = Atom::add_many((1..=3).map(|component| {
                                    Atom::num(spatial[component - 1])
                                        * function!(GS.delta_vec, GS.cind(component), index)
                                }));
                            }
                        })
                    })
                    .collect::<Vec<_>>();
                let mut tensors = NumeratorJetCache::default();
                assert!(
                    !tensors
                        .canonicalize(&values[1].undo_dots::<Aind>()?)?
                        .is_zero,
                    "{name}: the pinched physical cograph is nonzero"
                );
                let difference = (&values[0] - &values[1]).to_cof_dimension_invariants();
                let certificate = tensors.canonicalize(&difference.undo_dots::<Aind>()?)?;
                debug_tags!(#uv, #local, #on_shell, #cff, #inspect;
                    stage = "native_os_pinched_cograph_comparison",
                    graph = name, spatial = ?spatial, is_zero = certificate.is_zero,
                    file.difference = %difference.to_canonical_string(),
                    file.tensor = %certificate.key.to_canonical_string(),
                    "Exact finite-component comparison of all raised and pinched source maps"
                );
                assert!(
                    certificate.is_zero,
                    "{name}: raised versus pinched CFF at {spatial:?}"
                );
            }
        }
    }
    Ok(())
}
