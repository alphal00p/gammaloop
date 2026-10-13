//! Independent one-loop OS oracles in explicitly selected homogeneous LMBs.
//! Only denominator polynomials and finite tensor components are expanded.

use super::*;
use crate::{
    graph::parse::IntoGraph,
    initialisation::test_initialise,
    utils::{FUN_LIB, TENSORLIB},
    uv::{
        Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode,
        overlap_numerators::NumeratorJetCache,
    },
};
use idenso::shorthands::UndoShorthands;
use linnet::half_edge::involution::Hedge;
use spenso::{
    network::{ExecutionResult, Sequential, SmallestDegree, tags::SPENSO_TAG},
    structure::representation::{Minkowski, RepName},
};

pub(crate) fn fixture(massless: bool) -> Result<(Graph, OwnedForestNode, OwnedForestNode)> {
    let source = if massless {
        include_str!("../../../../../tests/resources/graphs/local_ir_quark_self_energy.dot")
    } else {
        include_str!("../../../../../tests/resources/graphs/local_os_top_self_energy.dot")
    };
    let graph: Graph = source.into_graph(&crate::utils::load_generic_model("sm"))?;
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

fn scalar_component(atom: &Atom) -> Result<Atom> {
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
    GS.resolve_os_momenta(
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
    .with(Atom::num(4))
    .simplify_metrics()
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
        *coefficients.entry(denominators).or_insert(Atom::Zero) += term.numerator;
    }
    Ok(coefficients)
}

fn assert_rational_zero(
    label: &str,
    graph: &Graph,
    current: &OwnedForestNode,
    atom: &Atom,
) -> Result<()> {
    let mut numerators = NumeratorJetCache::default();
    for (class, (denominators, numerator)) in rational_coefficients(graph, current, atom)?
        .into_iter()
        .enumerate()
    {
        let certificate = numerators.canonicalize(&numerator.undo_dots::<Aind>()?)?;
        let is_zero = certificate.is_zero;
        if !is_zero {
            debug_tags!(#uv, #local, #four_d, #trace;
                stage = "os_exact_oracle_failure",
                label = label,
                denominator_class = class,
                file.denominators = ?denominators,
                file.numerator = %numerator,
                file.contracted_tensor = %certificate.key,
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
) -> Result<(Atom, Atom, Atom)> {
    assert_eq!(
        current.lmb().loop_edges.iter().copied().collect::<Vec<_>>(),
        vec![EdgeIndex(2)]
    );
    let mass = graph[EdgeIndex(2)].mass_atom();
    let vector_slot = Minkowski {}.new_rep(4).to_symbolic([]);
    let k = function!(GS.emr_mom, 2, &vector_slot);
    let p = function!(GS.emr_mom, 0, &vector_slot);
    let k_zero = GS.emr_mom(EdgeIndex(2), GS.cind(0));
    let k_squared = dot(&k, &k);
    let numerator = graph
        .numerator(current.subgraph(), root.subgraph())
        .to_d_dim(4)
        .color_simplify()
        .get_single_atom()
        .expect("the source quark numerator is a single factor");
    let numerator = physical_expression(graph, current, &numerator);
    let massless_numerator = numerator.replace(mass.to_pattern()).with(Atom::Zero);
    let left_gamma_numerator = gamma_zero_at_port(graph, &massless_numerator, 1);
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
        let projected = (&numerator + Atom::num(sigma) * gamma_zero_at_port(graph, &numerator, 1))
            / Atom::num(2);
        let gluon = denominator(
            3,
            &k_squared - Atom::num(2 * sigma) * &mass * &k_zero + mass.pow(2),
            Atom::Zero,
        );
        sum + projected / (&quark * gluon)
    });
    let remainder = Atom::num(2)
        * (&massless_numerator * dot(&k, &p) - &mass * left_gamma_numerator * &k_zero)
        / vacuum.pow(3);
    Ok((on_shell, remainder, numerator))
}

#[test]
fn production_os_matches_equation_4_24_with_ordered_projectors_and_one_forest_sign() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &settings);
    let seed = Full4dCts::root();
    let actual = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
    let (on_shell, remainder, numerator) = one_loop_reference(&graph, &current, &root)?;
    assert_rational_zero(
        "production -T_OS equals the independent compact expression",
        &graph,
        &current,
        &(actual.atom() + &on_shell + &remainder),
    )?;
    let scheme = actual
        .active_sectors()
        .iter()
        .fold(Atom::Zero, |sum, sector| sum + sector.integration_atom());
    assert_rational_zero(
        "the integration policy retains only the signed vacuum remainder",
        &graph,
        &current,
        &(scheme + remainder),
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
fn production_os_massless_specialization_equals_soft_refined_degree_one() -> Result<()> {
    test_initialise()?;
    let (graph, mut current, root) = fixture(true)?;
    assert!(graph[EdgeIndex(2)].mass_atom().is_zero());
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &settings);
    let seed = Full4dCts::root();
    let on_shell = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
    current.spinney.renormalization_scheme = ApproximationType::IR;
    let soft = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
    assert_rational_zero(
        "OS at zero physical mass is H1",
        &graph,
        &current,
        &(on_shell.atom() - soft.atom()),
    )
}

#[test]
fn production_os_removes_all_one_loop_ultraviolet_tensor_coefficients() -> Result<()> {
    test_initialise()?;
    let (graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let ctx = UVCtx::new(&graph, &settings);
    let seed = Full4dCts::root();
    let bare = grow(&seed, &ctx, &current, &root)?;
    let subtraction = uv_limit(&seed, &ctx, &current, &root, &current, &root)?;
    let coefficients = rational_coefficients(&graph, &current, &(bare + subtraction.atom()))?;
    let ports = [0, 1].map(|hedge| {
        graph.underlying[Hedge(hedge)]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom()
    });
    let scale = Atom::var(GS.rescale);
    let mut nonzero_remainder = false;
    for row in 0..4 {
        for column in 0..4 {
            let mut jet = Atom::Zero;
            for (denominators, numerator) in &coefficients {
                let component = numerator
                    .replace(ports[0].to_pattern())
                    .with(GS.cind(row))
                    .replace(ports[1].to_pattern())
                    .with(GS.cind(column));
                let scalar = scalar_component(&(component * &graph.global_prefactor.num))?
                    / denominators
                        .iter()
                        .fold(Atom::one(), |product, factor| product * factor);
                let scaled = scalar
                    .replace(GS.emr_mom(EdgeIndex(2), W_.x___))
                    .with(GS.emr_mom(EdgeIndex(2), W_.x___) / &scale)
                    / scale.pow(4);
                jet +=
                    scaled.series_preserving_factors(GS.rescale, Atom::Zero.as_view(), 1, &[])?;
            }
            let series = jet.series(GS.rescale, Atom::Zero, 1)?;
            for power in [-1, 0] {
                assert!(
                    series.coefficient(power.into()).unwrap().expand().is_zero(),
                    "unsubtracted inverse-hard coefficient {power} in spinor component ({row},{column})"
                );
            }
            nonzero_remainder |= !series.coefficient(1.into()).unwrap().expand().is_zero();
        }
    }
    assert!(
        nonzero_remainder,
        "the UV check must not pass on an identically zero self-energy"
    );
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
    let ctx = UVCtx::new(&graph, &settings);
    let source = graph
        .numerator(current.subgraph(), root.subgraph())
        .to_d_dim(4)
        .color_simplify()
        .get_single_atom()
        .unwrap();
    let numerator = physical_expression(&graph, &current, &source);
    let mass = graph[EdgeIndex(2)].mass_atom();
    let rest = |atom: &Atom, energy: Atom| {
        let replaced = atom.replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() > 1
                && usize::try_from(momentum.get(0))
                    .is_ok_and(|edge| current.lmb().ext_edges.contains(&EdgeIndex(edge)))
            {
                let point = function!(GS.os_mom, 0, &energy);
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
        &(actual.atom() + on_shell + compact),
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
    let (mut graph, current, root) = fixture(false)?;
    let settings = UVgenerationSettings {
        generate_integrated: false,
        ..Default::default()
    };
    let local = uv_limit(
        &Full4dCts::root(),
        &UVCtx::new(&graph, &settings),
        &current,
        &root,
        &current,
        &root,
    )?;
    let (_, _, numerator) = one_loop_reference(&graph, &current, &root)?;
    let mass = graph[EdgeIndex(2)].mass_atom();
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
    for sector in local.projection_sectors() {
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
    // Independent poles in the quark-k routing: k0=E_m, sigma*m+E_0.
    // The fixed values below form two distinct exact rational kinematic
    // samples. Reconstruction above remains a general symbolic certificate;
    // these residue comparisons are additional contour checks, not its proof.
    let fermion_energy = Atom::num(17);
    let gluon_energy = Atom::num(15);
    let vacuum_energy = Atom::num(25);
    let mut residue = Atom::Zero;
    for sigma in [1, -1] {
        let shift = Atom::num(sigma) * &mass;
        let quark_value = at_energy(&numerator, &fermion_energy);
        let gluon_value = at_energy(&numerator, &(&shift + &gluon_energy));
        let project = |value: Atom| {
            (&value + Atom::num(sigma) * gamma_zero_at_port(&graph, &value, 1)) / Atom::num(2)
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
    let n_zero = at_energy(&numerator.replace(mass.to_pattern()).with(Atom::Zero), &z);
    let k_dot_p = (1..=3).fold(
        &z * GS.emr_mom(EdgeIndex(0), GS.cind(0)),
        |sum, component| {
            sum - GS.emr_mom(EdgeIndex(2), GS.cind(component))
                * GS.emr_mom(EdgeIndex(0), GS.cind(component))
        },
    );
    let vacuum_numerator =
        Atom::num(2) * (&n_zero * k_dot_p - &mass * gamma_zero_at_port(&graph, &n_zero, 1) * &z);
    let at_vacuum_pole = |value: Atom| value.replace(z_symbol).with(vacuum_energy.clone());
    let twice_energy = Atom::num(2) * &vacuum_energy;
    // Res A(z)/(z²-E²)^3 = 1/2 d²[A(z)/(z+E)^3] at z=E.
    residue += at_vacuum_pole(vacuum_numerator.derivative(z_symbol).derivative(z_symbol))
        / (Atom::num(2) * twice_energy.pow(3))
        - Atom::num(3) * at_vacuum_pole(vacuum_numerator.derivative(z_symbol))
            / twice_energy.pow(4)
        + Atom::num(6) * at_vacuum_pole(vacuum_numerator) / twice_energy.pow(5);
    // Clockwise contour and the single local forest sign cancel. The
    // remaining physical loop measure is i/(2*pi)^3.
    let expected = Atom::i() * residue / (Atom::num(8) * Atom::var(GS.pi).pow(3));
    for spatial in [[9, 12, 0], [9, -12, 0]] {
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
                    && (momentum.get_symbol() == GS.emr_mom || momentum.get_symbol() == GS.emr_vec)
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
        let difference = specialize(&actual) - specialize(&expected);
        let certificate =
            NumeratorJetCache::default().canonicalize(&difference.undo_dots::<Aind>()?)?;
        let is_zero = certificate.is_zero;
        if !is_zero {
            debug_tags!(#uv, #local, #four_d, #cff, #trace;
                stage = "os_dirac_residue_oracle_failure",
                spatial = ?spatial,
                file.difference = %difference,
                file.contracted_tensor = %certificate.key,
                "Nonzero factorized OS residue-oracle difference"
            );
        }
        assert!(
            is_zero,
            "physical Dirac CFF differs from the independent residues at spatial momentum {spatial:?}"
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
