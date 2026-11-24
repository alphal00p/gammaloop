use std::borrow::Cow;
use std::f64::consts::PI;
use std::slice;
use std::sync::Arc;

use eyre::eyre;
use idenso::color::{ColorSimplifier, ColorSimplifySettings};
use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use itertools::Itertools;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike, SubSetOps};
use symbolica::symbol;

use symbolica::domains::float::{Complex, Float, RealLike};
use symbolica::domains::rational::Rational;
use symbolica::evaluate::ExpressionEvaluator;
use vakint::VakintSettings;

use super::super::test_support::{ThermalCff, empty_forest_node, forest_node, host_node};
use super::super::thermal_soft::soft_branch_scope;
use super::super::{
    Direct3dApproximation, DirectCoordinateFrame, LOCAL_3D_SOFT_SCOPE, ZERO_MODE_ORDER,
    ZERO_MODE_TEMPERATURE, ZeroModeTerm, apply_taylor, class_zero_mode_temperature, coordinate_lmb,
    normalize_local_tags, start, start_family, zero_mode_order, zero_mode_temperature,
};
use super::*;
use crate::cff::expression::OrientationID;
use crate::cff::surface::LinearEnergyExpr;
use crate::dot;
use crate::graph::cuts::CutSet;
use crate::graph::parse::IntoGraph;
use crate::graph::{LMBext, LoopMomentumBasis};
use crate::initialisation::test_initialise;
use crate::integrands::process::param_builder::FnMapEntry;
use crate::settings::global::GenerationSettings;
use crate::utils;
use crate::utils::W_;
use crate::uv::approx::direct_3d::branches::DirectResidueKey;
use crate::uv::approx::direct_3d::{Direct3dCts, DirectSector};
use crate::uv::approx::final_integrand::{FinalIntegrandBuilder, FinalIntegrands};
use crate::uv::approx::integrated::{Integrated, IntegratedCts};
use crate::uv::approx::local_3d::{
    IntegratedLocalizationError, Local3DConceptualBranch, Local3DLoopRescaling,
    Local3DSignedBranch, Localizer,
};
use crate::uv::approx::local_4d::{self, Full4dCts, Local4dCts};
use crate::uv::approx::{CutStructure, ForestNodeLike, OrientationProjection, Rooted, UVCtx};
use crate::uv::hedge_poset::{OwnedForestNode, Wood};
use crate::uv::test_fixtures::{
    Ring, edge_set, local_3d_settings, ring_uv_settings, thermal_medium,
};
use crate::uv::thermal_zero_mode::ThermalZeroModeRegion;
use crate::uv::thermal_zero_mode::tests::static_region_through;
use crate::uv::uv_graph::UVE;
use crate::uv::{ApproximationType, UVgenerationSettings};

/// Distinct energies of the sunset's lines, at `T = 1`.
const SUNSET: [f64; 3] = [0.3, 0.7, 1.1];
/// Ring energies: `E0 = E1 = |q|`, then the two bubbles. No CFF surface
/// combination of these values vanishes.
const RING: [f64; 6] = [0.37, 0.37, 0.81, 0.64, 1.07, 0.52];
/// The quark chemical potential at `mu_B = 3`.
const QUARK_MU: f64 = 1.0;

fn sunset() -> Result<Graph> {
    Ok(dot!(
        digraph zero_mode_sunset {
            edge [particle="scalar_1"];
            node [num=1];
            v1 -> v2 [id=0];
            v1 -> v2 [id=1];
            v1 -> v2 [id=2];
        },
        "scalars"
    )?)
}

/// A sunset whose third line is a raised pair of equal propagators: the
/// generator writes their double pole with the occupation derivative.
fn raised_sunset() -> Result<Graph> {
    Ok(dot!(
        digraph zero_mode_raised_sunset {
            edge [particle="scalar_1"];
            node [num=1];
            v1 -> v2 [id=0];
            v1 -> v2 [id=1];
            v1 -> m [id=2];
            m -> v2 [id=3];
        },
        "scalars"
    )?)
}

/// One loop with two lines of different mass: both carry one frequency.
fn bubble() -> Result<Graph> {
    Ok(dot!(
        digraph zero_mode_bubble {
            node [num=1];
            v1 -> v2 [id=0 particle="scalar_1"];
            v2 -> v1 [id=1 particle="scalar_2"];
        },
        "scalars"
    )?)
}

/// A chain of four equal propagators carrying q closes on a bubble: nu = 4.
fn raised_chain() -> Result<Graph> {
    Ok(dot!(
        digraph zero_mode_chain {
            edge [particle="scalar_1"];
            node [num=1];
            a -> c1 [id=0 lmb_id=0];
            c1 -> c2 [id=1];
            c2 -> c3 [id=2];
            c3 -> b [id=3];
            b -> a [id=4 lmb_id=1];
            a -> b [id=5];
        },
        "scalars"
    )?)
}

impl Ring {
    /// The static second polarization loop `Pi_B(0)` at inverse temperature
    /// `beta`.
    fn static_second_loop(self, beta: f64, mu: f64, vacuum: bool) -> f64 {
        let [_, _, _, _, e4, e5] = RING;
        match (self, vacuum) {
            (Self::Tadpole, false) => bose_tadpole(beta, e4),
            (Self::Tadpole, true) => 1.0 / (2.0 * e4),
            (_, true) => vacuum_bubble(e4, e5),
            (Self::Quark, false) => static_bubble(e4, e5, |e| fermi_tadpole(beta, e, mu)),
            (Self::Mixed | Self::Ghost, false) => static_bubble(e4, e5, |e| bose_tadpole(beta, e)),
        }
    }

    /// The zero-mode region the policy assigns to the host on `edges`.
    fn region(
        self,
        graph: &Graph,
        edges: &[usize],
        vacuum_subtraction: bool,
    ) -> Result<ThermalZeroModeRegion> {
        self.expected_regions(graph, vacuum_subtraction)
            .remove(&edge_set(edges))
            .ok_or_else(|| eyre!("{self:?} has no zero-mode host on {edges:?}"))
    }
}

fn n_bose(beta: f64, energy: f64) -> f64 {
    1.0 / (beta * energy).exp_m1()
}

fn n_fermi(beta: f64, energy: f64) -> f64 {
    1.0 / ((beta * energy).exp() + 1.0)
}

/// `T sum_n 1/(w_n^2 + E^2) = (1/2 + n_B(E))/E` at `T = 1/beta`.
fn bose_tadpole(beta: f64, energy: f64) -> f64 {
    (0.5 + n_bose(beta, energy)) / energy
}

/// The fermionic sum with chemical potential `mu`, at `T = 1/beta`.
fn fermi_tadpole(beta: f64, energy: f64, mu: f64) -> f64 {
    (1.0 - n_fermi(beta, energy - mu) - n_fermi(beta, energy + mu)) / (2.0 * energy)
}

/// `T sum 1/((w^2 + a^2)(w^2 + b^2))` as a divided difference of tadpoles.
fn static_bubble(a: f64, b: f64, tadpole: impl Fn(f64) -> f64) -> f64 {
    (tadpole(a) - tadpole(b)) / (b * b - a * a)
}

/// `int dw/(2 pi) 1/((w^2 + a^2)(w^2 + b^2))`.
fn vacuum_bubble(a: f64, b: f64) -> f64 {
    1.0 / (2.0 * a * b * (a + b))
}

/// Assert that `actual` equals the nonzero reference `expected` to the
/// relative `tolerance`. Any zero would match a zero reference, so a test
/// states an expected zero explicitly instead.
fn assert_close(actual: f64, expected: f64, tolerance: f64, context: &str) {
    assert!(
        expected != 0.0,
        "{context}: the reference vanishes, compare with zero explicitly"
    );
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{context}: {actual} != {expected}"
    );
}

/// The real normalization that, with `i^L`, a root of `graph` multiplies
/// into `T sum prod 1/(w^2+E^2)`: each Minkowski propagator is
/// `-1/(w^2+E^2)` and each loop has the measure `i dw/(2 pi)^4`.
fn normalization(graph: &Graph) -> Result<f64> {
    let loops = i32::try_from(graph.loop_momentum_basis.loop_edges.len())?;
    let propagators = i32::try_from(graph.underlying.n_edges())?;
    Ok(f64::from(-1).powi(propagators) / (2.0 * PI).powi(3 * loops))
}

/// The production residue branches of the thermal CFF of `graph` and their
/// [`normalization`].
fn thermal_sector(
    graph: &mut Graph,
    vacuum_subtraction: bool,
) -> Result<(DirectResidueBranches, f64)> {
    let cff = ThermalCff::new(graph)?;
    let branches = cff.branches(graph, vacuum_subtraction)?;
    Ok((branches, normalization(graph)?))
}

/// One numerical point: on-shell energies by edge, the inverse temperature
/// `beta`, the chemical potential `mu` of the one particle species that has
/// one, and any further symbols.
#[derive(Clone)]
struct Point<'a> {
    energies: &'a [f64],
    beta: f64,
    mu: f64,
    extra: Vec<(Atom, f64)>,
}

impl<'a> Point<'a> {
    /// A point at `beta = 1` and `mu = 0`.
    fn new(energies: &'a [f64]) -> Self {
        Self {
            energies,
            beta: 1.0,
            mu: 0.0,
            extra: Vec::new(),
        }
    }

    fn with_beta(mut self, beta: f64) -> Self {
        self.beta = beta;
        self
    }

    fn with_mu(mut self, mu: f64) -> Self {
        self.mu = mu;
        self
    }

    fn with(mut self, symbol: Atom, value: f64) -> Self {
        self.extra.push((symbol, value));
        self
    }
}

/// Evaluate an atom with retained surfaces and distributions at several
/// points sharing their parameters; the remaining distributions take their
/// explicit equilibrium bodies at the point's `beta`, and the local tags left
/// in the atom take their final values ([`normalize_local_tags`]): the
/// zero-mode temperatures become `1/beta`. The factor `i^L` of the L
/// energy integrations is divided out numerically, and the real part is
/// returned.
fn evaluate_points(graph: &Graph, atom: &Atom, points: &[Point<'_>]) -> Result<Vec<f64>> {
    Ok(evaluate_points_complex(graph, atom, points)?
        .into_iter()
        .map(|value| value.re)
        .collect())
}

/// [`evaluate_points`] with the complex values, for atoms whose numerators
/// carry factors of `i`.
fn evaluate_points_complex(
    graph: &Graph,
    atom: &Atom,
    points: &[Point<'_>],
) -> Result<Vec<Complex<f64>>> {
    let PointEvaluator { evaluator, values } = point_evaluator(graph, atom, points)?;
    let mut evaluator = evaluator
        .map_coeff(&|coefficient| Complex::new(coefficient.re.to_f64(), coefficient.im.to_f64()));
    Ok(values
        .iter()
        .map(|values| {
            let values = values
                .iter()
                .map(|value| Complex::new(*value, 0.0))
                .collect::<Vec<_>>();
            let mut value = evaluator.evaluate_single(&values);
            for _ in &graph.loop_momentum_basis.loop_edges {
                value = Complex::new(value.im, -value.re);
            }
            value
        })
        .collect())
}

/// [`evaluate_points`] in `precision`-bit arithmetic at the same
/// double-precision inputs, rounded to double precision at the end: the
/// reference for the rounding error of the double-precision evaluation of
/// the same expression.
fn evaluate_points_precise(
    graph: &Graph,
    atom: &Atom,
    points: &[Point<'_>],
    precision: u32,
) -> Result<Vec<f64>> {
    let PointEvaluator { evaluator, values } = point_evaluator(graph, atom, points)?;
    let mut evaluator = evaluator.map_coeff_with_prec(
        &|coefficient| {
            Complex::new(
                coefficient.re.to_multi_prec_float(precision),
                coefficient.im.to_multi_prec_float(precision),
            )
        },
        precision,
    );
    Ok(values
        .iter()
        .map(|values| {
            let values = values
                .iter()
                .map(|value| {
                    Complex::new(
                        Float::with_val(precision, *value),
                        Float::with_val(precision, 0),
                    )
                })
                .collect::<Vec<_>>();
            let mut value = evaluator.evaluate_single(&values);
            for _ in &graph.loop_momentum_basis.loop_edges {
                value = Complex::new(value.im.clone(), -value.re.clone());
            }
            value.re.to_f64()
        })
        .collect())
}

/// An evaluator with exact coefficients and the parameter values of the
/// points it evaluates.
struct PointEvaluator {
    evaluator: ExpressionEvaluator<Complex<Rational>>,
    /// The on-shell energies, `beta`, the chemical potential and the further
    /// symbols of each point.
    values: Vec<Vec<f64>>,
}

/// The [`PointEvaluator`] of `atom` at `points`.
fn point_evaluator(graph: &Graph, atom: &Atom, points: &[Point<'_>]) -> Result<PointEvaluator> {
    let atom = normalize_local_tags(&GS.lower_esurfaces(atom))
        .replace(function!(GS.thermal_weight_wrapper, W_.a_))
        .with(W_.a_);
    let edges = (0..graph.underlying.n_edges()).map(EdgeIndex);
    let explicit = graph.make_thermal_distributions_explicit(
        &atom,
        MediumMode::ThermodynamicEquilibrium,
        edges.clone(),
    )?;
    // Antiparticle edges carry the opposite chemical-potential parameter.
    let chemical_potentials = edges
        .filter_map(|edge| {
            let pdg = graph[edge].particle_pdg_code()?;
            Some((graph[edge].chemical_potential_atom()?, pdg))
        })
        .collect::<BTreeMap<_, _>>();
    assert!(
        chemical_potentials
            .values()
            .map(|pdg| pdg.abs())
            .collect::<BTreeSet<_>>()
            .len()
            <= 1,
        "a point has one chemical potential: {chemical_potentials:?}"
    );
    let first = points.first().expect("at least one evaluation point");
    let symbols = (0..first.energies.len())
        .map(|edge| ose_atom_from_index(EdgeIndex(edge)))
        .chain([Atom::var(GS.inverse_temperature)])
        .chain(chemical_potentials.keys().cloned())
        .chain(first.extra.iter().map(|(symbol, _)| symbol.clone()))
        .collect::<Vec<_>>();
    let evaluator = explicit
        .evaluator(&symbols)
        .function_map(graph.param_builder.fn_map.clone())
        .build()
        .map_err(|error| eyre!("{error}"))?;
    let values = points
        .iter()
        .map(|point| {
            point
                .energies
                .iter()
                .copied()
                .chain([point.beta])
                .chain(
                    chemical_potentials
                        .values()
                        .map(|pdg| if *pdg > 0 { point.mu } else { -point.mu }),
                )
                .chain(point.extra.iter().map(|(_, value)| *value))
                .collect()
        })
        .collect();
    Ok(PointEvaluator { evaluator, values })
}

fn evaluate_atom(graph: &Graph, atom: &Atom, point: &Point<'_>) -> Result<f64> {
    Ok(evaluate_points(graph, atom, slice::from_ref(point))?[0])
}

fn evaluate(graph: &Graph, branches: &DirectResidueBranches, point: &Point<'_>) -> Result<f64> {
    evaluate_atom(graph, &branches.factorized_sum(), point)
}

fn project(
    graph: &Graph,
    branches: &DirectResidueBranches,
    classes: Vec<ZeroModeClass>,
    sector: ZeroModeSector,
) -> Result<BTreeMap<Vec<u32>, DirectResidueBranches>> {
    Ok(ZeroModeProjection::new(classes, sector)?.coefficients(graph, branches)?)
}

fn class(graph: &Graph, edges: &[usize], orders: &[u32]) -> Result<ZeroModeClass> {
    Ok(ZeroModeClass::new(
        graph,
        edges.iter().copied().map(EdgeIndex),
        orders.iter().copied(),
    )?)
}

#[test]
fn sunset_zero_mode_is_the_static_bubble() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let [a, b, c] = SUNSET;
    let projected = project(
        &graph,
        &branches,
        vec![class(&graph, &[0], &[2])?],
        ZeroModeSector::Complete,
    )?;
    let projected = evaluate(&graph, &projected[&vec![2]], &Point::new(&SUNSET))? / normalization;
    // [x^-2] F = T (h(b) - h(c))/(c^2 - b^2) with h(E) = (1/2 + n_B(E))/E.
    assert_close(projected * a * a, 1.687959575607, 1e-10, "[x^-2] F");
    assert_close(projected, 18.755106395636, 1e-10, "P_a F");
    assert_close(
        projected,
        static_bubble(b, c, |e| bose_tadpole(1.0, e)) / (a * a),
        1e-12,
        "P_a F against the closed form",
    );
    Ok(())
}

#[test]
fn raised_hard_lines_keep_the_occupation_derivative() -> Result<()> {
    test_initialise()?;
    let mut graph = raised_sunset()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let mut hard_derivative = false;
    branches.factorized_sum().as_view().visitor(&mut |view| {
        hard_derivative |= ThermalDistributionCall::try_from(view)
            .is_ok_and(|call| call.edge != EdgeIndex(0) && call.derivative_order == 1);
        !hard_derivative
    });
    assert!(
        hard_derivative,
        "the raised hard pair must carry an occupation derivative"
    );
    let projected = project(
        &graph,
        &branches,
        vec![class(&graph, &[0], &[2])?],
        ZeroModeSector::Complete,
    )?;
    // P_a F = T sum 1/((w^2+b^2)(w^2+c^2)^2)/a^2 at the exactly equal
    // energies of the raised pair c; without the n(1+n)/T part of h'(c) the
    // value would be 20.257602026309.
    let [a, b, c] = SUNSET;
    let value = evaluate(&graph, &projected[&vec![2]], &Point::new(&[a, b, c, c]))? / normalization;
    assert_close(value * a * a, 1.393935790393, 1e-10, "[x^-2] F");
    assert_close(value, 15.488175448815, 1e-10, "P_a F");
    Ok(())
}

#[test]
fn zero_mode_projections_commute_without_redundant_temperatures() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let [a, b, c] = SUNSET;
    let expected = 1.0 / (a * a * b * b * c * c);
    assert_close(expected, 18.740278480538, 1e-11, "T^2/(a^2 b^2 c^2)");
    for edges in [vec![0, 1], vec![1, 0], vec![0, 1, 2], vec![2, 0, 1]] {
        let classes = edges
            .iter()
            .map(|edge| class(&graph, &[*edge], &[2]))
            .collect::<Result<Vec<_>>>()?;
        let projected = project(&graph, &branches, classes, ZeroModeSector::Complete)?;
        let orders = vec![2; edges.len()];
        assert_eq!(projected.keys().collect::<Vec<_>>(), [&orders]);
        // Once two frequencies are static, the third one is too: its
        // projection restores its propagator without another factor of T.
        let value = evaluate(&graph, &projected[&orders], &Point::new(&SUNSET))? / normalization;
        assert_close(value, expected, 1e-10, &format!("P over {edges:?}"));
    }
    Ok(())
}

#[test]
fn raised_classes_have_only_the_zero_mode_pole_without_vacuum_subtraction() -> Result<()> {
    test_initialise()?;
    let mut graph = raised_chain()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let energies = [0.37, 0.37, 0.37, 0.37, 1.07, 0.52];
    let projected = project(
        &graph,
        &branches,
        vec![class(&graph, &[0, 1, 2, 3], &[1, 2, 3, 4, 5, 6, 7, 8])?],
        ZeroModeSector::Complete,
    )?;
    let q = energies[0];
    let static_value = static_bubble(energies[4], energies[5], |e| bose_tadpole(1.0, e));
    // nu = 4 reaches the third occupation derivative of the class.
    let zero_mode = evaluate(&graph, &projected[&vec![8]], &Point::new(&energies))? / normalization;
    assert_close(
        zero_mode * q.powi(8),
        static_value,
        1e-10,
        "[x^-8] = T G(0)",
    );
    for order in 1..8_u32 {
        let value =
            evaluate(&graph, &projected[&vec![order]], &Point::new(&energies))? / normalization;
        assert!(
            (value * q.powi(order as i32)).abs() <= 1e-10 * static_value,
            "[x^-{order}] of an unsubtracted sector must vanish: {value}"
        );
    }
    Ok(())
}

#[test]
fn quark_ring_zero_mode_factorizes_into_static_polarizations() -> Result<()> {
    assert_ring_zero_mode_factorizes(Ring::Quark)
}

#[test]
fn gluon_ring_zero_mode_factorizes_into_static_polarizations() -> Result<()> {
    assert_ring_zero_mode_factorizes(Ring::Mixed)
}

#[test]
fn ghost_ring_zero_mode_factorizes_into_static_polarizations() -> Result<()> {
    assert_ring_zero_mode_factorizes(Ring::Ghost)
}

#[test]
fn tadpole_ring_zero_mode_factorizes_into_static_polarizations() -> Result<()> {
    assert_ring_zero_mode_factorizes(Ring::Tadpole)
}

/// `[x^-4] = T Pi_A(0) Pi_B(0)` with and without vacuum subtraction, the
/// vacuum 4D-soft `[x^-3] = -Pi_A^vac Pi_B^vac/4` only with it, and no other
/// pole, at unit and at another temperature; only the zero mode carries the
/// class temperature `T_zm(0)`, which the final normalization makes the
/// temperature of the remaining occupations.
fn assert_ring_zero_mode_factorizes(ring: Ring) -> Result<()> {
    test_initialise()?;
    let [q, _, e2, e3, _, _] = RING;
    let zero_mode_temperature_value = Atom::var(symbol!("zero_mode_test_temperature"));
    let media = [(1.0, 0.0), (1.0, QUARK_MU), (0.5, 0.0), (0.5, QUARK_MU)];
    let points = media.map(|(beta, mu)| Point::new(&RING).with_beta(beta).with_mu(mu));
    let expected = media.map(|(beta, mu)| {
        static_bubble(e2, e3, |e| fermi_tadpole(beta, e, mu))
            * ring.static_second_loop(beta, mu, false)
            / beta
    });
    for vacuum_subtraction in [false, true] {
        let mut graph = ring.topology()?;
        let (branches, normalization) = thermal_sector(&mut graph, vacuum_subtraction)?;
        let projected = project(
            &graph,
            &branches,
            vec![class(&graph, &[0, 1], &[1, 2, 3, 4])?],
            ZeroModeSector::Complete,
        )?;
        // [x^-a] q^a in the normalization of T sum prod 1/(w^2+E^2).
        let coefficient = |order: u32| -> Result<Vec<f64>> {
            Ok(
                evaluate_points(&graph, &projected[&vec![order]].factorized_sum(), &points)?
                    .into_iter()
                    .map(|value| value * q.powi(order as i32) / normalization)
                    .collect(),
            )
        };
        let odd_expected = if vacuum_subtraction {
            -vacuum_bubble(e2, e3) * ring.static_second_loop(1.0, 0.0, true) / 4.0
        } else {
            0.0
        };
        let (zero_mode, odd, even, first) = (
            coefficient(4)?,
            coefficient(3)?,
            coefficient(2)?,
            coefficient(1)?,
        );
        for (index, (beta, mu)) in media.iter().enumerate() {
            let context =
                format!("{ring:?}, vacuum subtraction {vacuum_subtraction}, beta {beta}, mu {mu}");
            let scale = expected[index].abs();
            assert_close(zero_mode[index], expected[index], 1e-10, &context);
            let odd = odd[index];
            assert!(
                (odd - odd_expected).abs() <= 1e-10 * scale,
                "{context}: [x^-3] q^3 = {odd} != {odd_expected}"
            );
            let even = even[index];
            assert!(even.abs() <= 1e-10 * scale, "{context}: [x^-2] = {even}");
            let first = first[index];
            // With vacuum subtraction [x^-1] is a further vacuum term.
            assert!(
                vacuum_subtraction || first.abs() <= 1e-10 * scale,
                "{context}: [x^-1] = {first}"
            );
        }

        // The zero mode is linear in its class temperature, which the hard
        // occupations at beta = 1 do not share; the vacuum-soft pole has none.
        let tilted = [1.0, 2.0].map(|temperature| {
            Point::new(&RING)
                .with_mu(QUARK_MU)
                .with(zero_mode_temperature_value.clone(), temperature)
        });
        let vacuum_soft = vacuum_subtraction.then_some((3, 1.0));
        for (order, ratio) in [(4, 2.0)].into_iter().chain(vacuum_soft) {
            let detached = projected[&vec![order]]
                .factorized_sum()
                .replace(class_zero_mode_temperature(EdgeIndex(0)))
                .with(zero_mode_temperature_value.clone());
            let tilted = evaluate_points(&graph, &detached, &tilted)?;
            assert_close(
                tilted[1] / tilted[0],
                ratio,
                1e-12,
                &format!("{ring:?} [x^-{order}] in T_zm"),
            );
        }
    }
    Ok(())
}

#[test]
fn hard_vacuum_sector_keeps_the_vacuum_soft_coefficient() -> Result<()> {
    test_initialise()?;
    for ring in Ring::ALL {
        let mut graph = ring.topology()?;
        let (branches, _) = thermal_sector(&mut graph, true)?;
        let [complete, hard_vacuum] =
            [ZeroModeSector::Complete, ZeroModeSector::HardVacuum].map(|sector| {
                project(
                    &graph,
                    &branches,
                    vec![class(&graph, &[0, 1], &[3])?],
                    sector,
                )
            });
        let (complete, hard_vacuum) = (complete?, hard_vacuum?);
        assert!(
            !hard_vacuum[&vec![3]]
                .factorized_sum()
                .contains_symbol(GS.thermal_distribution),
            "{ring:?}: the hard-vacuum sector retained an occupation"
        );
        for mu in [0.0, QUARK_MU] {
            let point = Point::new(&RING).with_mu(mu);
            assert_close(
                evaluate(&graph, &hard_vacuum[&vec![3]], &point)?,
                evaluate(&graph, &complete[&vec![3]], &point)?,
                1e-10,
                &format!("{ring:?}, mu {mu}"),
            );
        }
    }
    Ok(())
}

/// The summands of a root: its terms, or those of its one sum factor times
/// the remaining factors. They stand in for the CFF variants of the root.
fn summands(root: &Atom) -> Vec<Atom> {
    match root.as_view() {
        AtomView::Add(sum) => sum.iter().map(|term| term.to_owned()).collect(),
        AtomView::Mul(product) => {
            let mut sums = product
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)));
            match (sums.next(), sums.next()) {
                (Some(AtomView::Add(sum)), None) => {
                    let rest = product
                        .iter()
                        .filter(|factor| !matches!(factor, AtomView::Add(_)))
                        .fold(Atom::one(), |rest, factor| rest * factor.to_owned());
                    sum.iter().map(|term| term.to_owned() * &rest).collect()
                }
                _ => vec![root.clone()],
            }
        }
        _ => vec![root.clone()],
    }
}

/// The largest summand of any residue row of `coefficient`, at each point.
fn largest_summands(
    graph: &Graph,
    coefficient: &DirectResidueBranches,
    points: &[Point<'_>],
) -> Result<Vec<f64>> {
    let mut largest = vec![0.0_f64; points.len()];
    for (_, roots) in coefficient.iter_keys() {
        for (_, root) in roots.resolved()?.iter() {
            for term in summands(root) {
                let values = evaluate_points(graph, &term, points)?;
                for (largest, value) in largest.iter_mut().zip(values) {
                    *largest = largest.max(value.abs());
                }
            }
        }
    }
    Ok(largest)
}

#[test]
fn hard_vacuum_summands_stay_bounded_near_a_hard_threshold() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Quark.topology()?;
    let (branches, _) = thermal_sector(&mut graph, true)?;
    // E2 - E3 - E4 + E5 = delta along E5. The odd coefficient is regular
    // there in any sector, but only the hard-vacuum one keeps its summands
    // bounded; the complete sector's summands cancel only in the sum.
    let [q, _, e2, e3, e4, _] = RING;
    let approach = [1e-3, 1e-6].map(|delta| [q, q, e2, e3, e4, e4 - (e2 - e3) + delta]);
    let points = approach
        .each_ref()
        .map(|energies| Point::new(energies.as_slice()));
    let largest = |sector| -> Result<Vec<f64>> {
        let projected = project(
            &graph,
            &branches,
            vec![class(&graph, &[0, 1], &[3])?],
            sector,
        )?;
        largest_summands(&graph, &projected[&vec![3]], &points)
    };
    let complete = largest(ZeroModeSector::Complete)?;
    assert!(
        complete[1] > 100.0 * complete[0],
        "the complete sector's summands must grow towards the threshold: {complete:?}"
    );
    let hard_vacuum = largest(ZeroModeSector::HardVacuum)?;
    assert!(
        hard_vacuum[0].is_finite() && hard_vacuum[1] <= 2.0 * hard_vacuum[0],
        "the hard-vacuum summands grow towards the threshold: {hard_vacuum:?}"
    );
    Ok(())
}

#[test]
fn static_conditions_implied_by_others_add_no_temperature_on_rings() -> Result<()> {
    test_initialise()?;
    let [q, _, e2, e3, e4, e5] = RING;
    for ring in [Ring::Mixed, Ring::Ghost] {
        let mut graph = ring.topology()?;
        let (branches, normalization) = thermal_sector(&mut graph, false)?;
        let expected =
            static_bubble(e2, e3, |e| fermi_tadpole(1.0, e, 0.0)) / (q.powi(4) * e4 * e4 * e5 * e5);
        for classes in [
            vec![class(&graph, &[0, 1], &[4])?, class(&graph, &[4], &[2])?],
            vec![class(&graph, &[4], &[2])?, class(&graph, &[0, 1], &[4])?],
            vec![
                class(&graph, &[0, 1], &[4])?,
                class(&graph, &[4], &[2])?,
                class(&graph, &[5], &[2])?,
            ],
            vec![
                class(&graph, &[5], &[2])?,
                class(&graph, &[4], &[2])?,
                class(&graph, &[0, 1], &[4])?,
            ],
        ] {
            let context = format!("{ring:?} {classes:?}");
            let projected = project(&graph, &branches, classes, ZeroModeSector::Complete)?;
            let (_, coefficient) = projected
                .iter()
                .next()
                .expect("one order tuple was requested");
            let value = evaluate(&graph, coefficient, &Point::new(&RING))? / normalization;
            // T^2 Pi_A(0) / (q^4 E4^2 E5^2): two independent static sums.
            assert_close(value, expected, 1e-9, &context);
        }
    }
    Ok(())
}

#[test]
fn products_more_singular_than_the_zero_mode_keep_exact_coefficients() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Quark.topology()?;
    let (branches, _) = thermal_sector(&mut graph, true)?;
    // [x^-4] (F/x) = [x^-3] F. Dividing by x deepens every product by one
    // order beyond those of F, so the distributions must keep more orders.
    let shifted = branches.map(|root| root / ose_atom_from_index(EdgeIndex(0)));
    let shifted = project(
        &graph,
        &shifted,
        vec![class(&graph, &[0, 1], &[4])?],
        ZeroModeSector::Complete,
    )?;
    let direct = project(
        &graph,
        &branches,
        vec![class(&graph, &[0, 1], &[3])?],
        ZeroModeSector::Complete,
    )?;
    let q = RING[0];
    for mu in [0.0, QUARK_MU] {
        let point = Point::new(&RING).with_mu(mu);
        let direct = evaluate(&graph, &direct[&vec![3]], &point)? * q.powi(3);
        assert!(
            direct.abs() > 0.0,
            "the vacuum-soft coefficient must not vanish"
        );
        assert_close(
            evaluate(&graph, &shifted[&vec![4]], &point)? * q.powi(4),
            direct,
            1e-10,
            &format!("mu {mu}"),
        );
    }
    Ok(())
}

#[test]
fn polynomial_class_energies_shift_the_laurent_coefficients() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let (branches, _) = thermal_sector(&mut graph, false)?;
    let [a, b, _] = SUNSET;
    // A numerator or contact surface x + E1 adds no pole of its own:
    // [x^-2] ((x + E1) F) = E1 [x^-2] F and [x^-1] ((x + E1) F) = [x^-2] F,
    // since the unsubtracted F has no [x^-1] and no deeper pole.
    let surface = ose_atom_from_index(EdgeIndex(0)) + ose_atom_from_index(EdgeIndex(1));
    let contact = branches.map(|root| root * &surface);
    let zero_mode = vec![class(&graph, &[0], &[1, 2])?];
    let plain = project(
        &graph,
        &branches,
        zero_mode.clone(),
        ZeroModeSector::Complete,
    )?;
    let shifted = project(&graph, &contact, zero_mode, ZeroModeSector::Complete)?;
    let point = Point::new(&SUNSET);
    let static_value = evaluate(&graph, &plain[&vec![2]], &point)?;
    assert!(static_value.abs() > 0.0, "the zero mode must not vanish");
    assert!(
        evaluate(&graph, &plain[&vec![1]], &point)?.abs() <= 1e-10 * static_value.abs(),
        "an unsubtracted sector has no [x^-1]"
    );
    // P_a restores OSE(0)^-a = a^-a.
    assert_close(
        evaluate(&graph, &shifted[&vec![2]], &point)?,
        b * static_value,
        1e-10,
        "[x^-2] of (x + E1) F",
    );
    assert_close(
        evaluate(&graph, &shifted[&vec![1]], &point)?,
        a * static_value,
        1e-10,
        "[x^-1] of (x + E1) F",
    );
    Ok(())
}

#[test]
fn retained_surfaces_in_numerators_project_like_their_sums() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let (branches, _) = thermal_sector(&mut graph, false)?;
    // The scalar-pole cancellation of an ordinary U brings retained surfaces
    // over a common denominator, as on the gluon bubble of the mixed ring,
    // which leaves surfaces carrying the class energy in numerator position:
    // (3 S1^2 + 5 S2 + 7)/(S1 S2) with S1 = x - E1 + E2 and S2 = x + E1 - E2.
    // They are linear energy sums, polynomial in x, and project like the
    // lowered sums.
    let ose = |edge| ose_atom_from_index(EdgeIndex(edge));
    let surfaces = [ose(0) - ose(1) + ose(2), ose(0) + ose(1) - ose(2)];
    let factor = |[first, second]: &[Atom; 2]| {
        (Atom::num(3) * first.pow(2) + Atom::num(5) * second + Atom::num(7)) / (first * second)
    };
    let retained = factor(&surfaces.clone().map(|surface| GS.wrap_esurface(&surface)));
    assert!(retained.contains_symbol(GS.esurface));
    let lowered = factor(&surfaces);
    let zero_mode = vec![class(&graph, &[0], &[1, 2])?];
    let retained = project(
        &graph,
        &branches.map(|root| root * &retained),
        zero_mode.clone(),
        ZeroModeSector::Complete,
    )?;
    let lowered = project(
        &graph,
        &branches.map(|root| root * &lowered),
        zero_mode,
        ZeroModeSector::Complete,
    )?;
    let point = Point::new(&SUNSET);
    for order in [1, 2] {
        let expected = evaluate(&graph, &lowered[&vec![order]], &point)?;
        assert!(expected.abs() > 0.0, "[x^-{order}] must not vanish");
        assert_close(
            evaluate(&graph, &retained[&vec![order]], &point)?,
            expected,
            1e-12,
            &format!("[x^-{order}] of the retained surfaces"),
        );
    }
    Ok(())
}

#[test]
fn numerators_follow_the_deformed_pole_map() -> Result<()> {
    test_initialise()?;
    let mut graph = bubble()?;
    let cff = ThermalCff::new(&mut graph)?;
    let orientation = cff.orientation(false);
    let branches = cff.branches(&graph, false)?;
    // N(q0) = M^2 - q0^2, i.e. w^2 + M^2 in Matsubara frequencies. Equal to
    // the inverse a-propagator at M = a, it must not cancel that propagator's
    // zero mode before the deformation: the static value is T M^2/(a^2 b^2).
    let mass = Atom::var(symbol!("zero_mode_numerator_mass"));
    let numerator = mass.pow(2) - GS.emr_mom(EdgeIndex(0), GS.cind(0)).pow(2);
    let with_numerator = branches.multiply_key_mapped(
        orientation,
        &graph,
        &numerator,
        DirectResidueBranches::numerator_scope(),
    )?;
    let zero_mode = class(&graph, &[0], &[2])?;
    let plain = project(
        &graph,
        &branches,
        vec![zero_mode.clone()],
        ZeroModeSector::Complete,
    )?;
    let deformed = project(
        &graph,
        &with_numerator,
        vec![zero_mode],
        ZeroModeSector::Complete,
    )?;
    let energies = [0.3, 0.7];
    for m in [0.3, 0.45] {
        let point = Point::new(&energies).with(mass.clone(), m);
        let ratio = evaluate(&graph, &deformed[&vec![2]], &point)?
            / evaluate(&graph, &plain[&vec![2]], &point)?;
        assert_close(ratio, m * m, 1e-12, &format!("M = {m}"));
    }
    Ok(())
}

#[test]
fn naive_static_limits_do_not_reproduce_the_zero_mode() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let [a, b, c] = SUNSET;
    let correct = static_bubble(b, c, |e| bose_tadpole(1.0, e)) / (a * a);
    let x = symbol!("zero_mode_naive_energy");
    let root = GS
        .lower_esurfaces(&branches.factorized_sum())
        .replace(function!(GS.thermal_weight_wrapper, W_.a_))
        .with(W_.a_);

    // Replacing every Bose occupation by its classical limit T/E also
    // touches the hard lines, whose occupations are not soft.
    let classical = ThermalDistributionCall::replace_selected(
        root.as_view(),
        |_| Ok::<_, Report>(true),
        |call| {
            let energy = symbol!("zero_mode_naive_occupation_energy");
            let mut body = (Atom::one() + &call.thermal_sign) / Atom::num(2)
                + Atom::var(GS.inverse_temperature).pow(-1) / Atom::var(energy);
            if call.derivative_order > 0 {
                body = Atom::var(GS.inverse_temperature).pow(-1) / Atom::var(energy);
            }
            for _ in 0..call.derivative_order {
                body = body.derivative(energy);
            }
            Ok(body.replace(Atom::var(energy)).with(call.energy))
        },
    )?;
    let classical = classical
        .replace(ose_atom_from_index(EdgeIndex(0)))
        .with(Atom::var(x))
        .series(x, Atom::Zero, -2)?
        .coefficient(Rational::from(-2))
        .unwrap_or(Atom::Zero)
        * ose_atom_from_index(EdgeIndex(0)).pow(-2);
    let classical = evaluate_atom(&graph, &classical, &Point::new(&SUNSET))? / normalization;
    assert!(
        (classical / correct - 1.0).abs() > 1e-2,
        "the classical occupations reproduced P_a F: {classical} vs {correct}"
    );

    Ok(())
}

#[test]
fn zero_mode_classes_must_share_one_bosonic_frequency() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Quark.topology()?;
    let refused = |graph: &Graph, edges: &[usize], orders: &[u32]| {
        ZeroModeClass::new(
            graph,
            edges.iter().copied().map(EdgeIndex),
            orders.iter().copied(),
        )
        .expect_err("the class must be rejected")
    };
    let charged = dot!(
        digraph zero_mode_charged {
            node [num=1]; edge [num=1];
            a -> b [id=0 particle="G+"];
            b -> a [id=1 particle="G+"];
        }
    )?;
    let error = refused(&graph, &[2], &[2]);
    assert!(
        matches!(
            error,
            ZeroModeError::NoZeroMode {
                reason: NoZeroMode::Fermion,
                ..
            }
        ),
        "{error}"
    );
    let error = refused(&charged, &[0], &[2]);
    assert!(
        matches!(
            error,
            ZeroModeError::NoZeroMode {
                reason: NoZeroMode::ChemicalPotential,
                ..
            }
        ),
        "{error}"
    );
    // Different frequencies, then one frequency at two masses.
    for (graph, edges) in [(sunset()?, [0, 1]), (bubble()?, [0, 1])] {
        let error = refused(&graph, &edges, &[2]);
        assert!(matches!(error, ZeroModeError::NotRaised { .. }), "{error}");
    }
    let error = refused(&graph, &[0, 1], &[5]);
    assert!(
        matches!(error, ZeroModeError::InvalidOrders { .. }),
        "{error}"
    );
    let error = refused(&graph, &[], &[2]);
    assert!(matches!(error, ZeroModeError::EmptyClass), "{error}");
    let error = refused(&graph, &[9], &[2]);
    assert!(
        matches!(error, ZeroModeError::UnknownEdge { .. }),
        "{error}"
    );

    let zero_mode = class(&graph, &[0, 1], &[4])?;
    let refused = |classes: Vec<ZeroModeClass>, sector| {
        ZeroModeProjection::new(classes, sector).expect_err("the projection must be rejected")
    };
    let error = refused(Vec::new(), ZeroModeSector::Complete);
    assert!(matches!(error, ZeroModeError::NoClass), "{error}");
    let error = refused(
        vec![zero_mode.clone(), class(&graph, &[1], &[2])?],
        ZeroModeSector::Complete,
    );
    assert!(
        matches!(error, ZeroModeError::OverlappingClasses { .. }),
        "{error}"
    );
    // The hard-vacuum sector only has the vacuum-soft order [x^-3].
    let error = refused(vec![zero_mode], ZeroModeSector::HardVacuum);
    assert!(
        matches!(error, ZeroModeError::HardVacuumOrders { .. }),
        "{error}"
    );
    Ok(())
}

#[test]
fn zero_mode_guards_reject_unsupported_sectors() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Quark.topology()?;
    let (branches, _) = thermal_sector(&mut graph, false)?;
    let zero_mode = class(&graph, &[0, 1], &[4])?;
    let refusal = |branches: &DirectResidueBranches| -> Result<ZeroModeError> {
        Ok(
            ZeroModeProjection::new(vec![zero_mode.clone()], ZeroModeSector::Complete)?
                .coefficients(&graph, branches)
                .expect_err("the projection must reject this sector"),
        )
    };
    let ose = |edge| ose_atom_from_index(EdgeIndex(edge));
    let routed = function!(
        GS.on_shell_energy,
        0,
        GS.emr_vec(EdgeIndex(0), GS.cind(1)).pow(2)
    );

    // A class energy in a denominator that is not a linear energy sum.
    let error = refusal(&branches.map(|root| root / (ose(0).pow(2) + ose(2))))?;
    assert!(
        matches!(error, ZeroModeError::FormalArgument { .. }),
        "{error}"
    );
    // A class energy in a formal function argument.
    let error =
        refusal(&branches.map(|root| root * function!(symbol!("zero_mode_test_formal"), ose(0))))?;
    assert!(
        matches!(error, ZeroModeError::FormalArgument { .. }),
        "{error}"
    );
    // A class distribution away from the class energy escapes the
    // substitution, and so does a routed on-shell energy of a class edge.
    let off_shell =
        ThermalDistributionCall::atom_at_energy(0, 0, 1, -1, 1, symbol!("zero_mode_test_energy"));
    let error = refusal(&branches.map(|root| root * &off_shell))?;
    assert!(
        matches!(error, ZeroModeError::DeformedDistribution { .. }),
        "{error}"
    );
    let error = refusal(&branches.map(|root| root * &routed))?;
    assert!(
        matches!(error, ZeroModeError::RoutedEnergy { .. }),
        "{error}"
    );
    // A summand that hides the class energy in a routed one has no pole of
    // its own, but beside summands that still show the class energy it may
    // cancel one of theirs only in the sum. Alone it projects to zero.
    let error = refusal(&branches.map(|root| root + &routed))?;
    assert!(
        matches!(error, ZeroModeError::RoutedEnergy { .. }),
        "{error}"
    );
    let hidden = project(
        &graph,
        &branches.map(|_| routed.clone()),
        vec![zero_mode.clone()],
        ZeroModeSector::Complete,
    )?;
    assert!(hidden[&vec![4]].factorized_sum().is_zero());

    // Shared bodies are polynomial in x and never hide a class energy.
    let x = Atom::var(GS.zero_mode_energy);
    let y = Atom::var(symbol!("zero_mode_test_body"));
    zero_mode.check_numerator((x.pow(2) * &y + &y).as_view())?;
    for (body, rejected) in [
        ((x.pow(2) + &y).pow(-1), "inverse"),
        (x.pow(Atom::num(1) / Atom::num(2)), "root"),
        (&off_shell * &y, "class distribution"),
        (routed * &y, "routed energy"),
        (function!(GS.thermal_weight_wrapper, &x), "medium weight"),
    ] {
        let error = zero_mode
            .check_numerator(body.as_view())
            .expect_err(rejected);
        assert!(
            matches!(
                error,
                ZeroModeError::NotPolynomial { .. }
                    | ZeroModeError::FormalArgument { .. }
                    | ZeroModeError::RoutedEnergy { .. }
            ),
            "{rejected}: {error}"
        );
    }
    Ok(())
}

/// The coordinate frame of a nested operation on `node`.
fn frame(node: &OwnedForestNode) -> DirectCoordinateFrame {
    DirectCoordinateFrame {
        active_subgraph: node.subgraph().clone(),
        lmb: node.lmb().clone(),
    }
}

#[test]
fn coordinate_lmb_full_contraction_keeps_retained_carriers() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Mixed.topology()?;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    // The union host and the factorized join it sits above share one filter.
    let host = forest_node(&graph, &[2, 3, 4, 5], ApproximationType::IR, 4)?;
    let join = forest_node(&graph, &[2, 3, 4, 5], ApproximationType::MUV, 4)?;
    let first = forest_node(&graph, &[2, 3], ApproximationType::MUV, 2)?;
    let second = forest_node(&graph, &[4, 5], ApproximationType::MUV, 2)?;
    assert_eq!(first.lmb().loop_edges.raw, vec![EdgeIndex(2)]);
    assert_eq!(second.lmb().loop_edges.raw, vec![EdgeIndex(4)]);
    let given = InternalSubGraph::try_new(join.subgraph().clone(), graph.as_ref())
        .ok_or_else(|| eyre!("the join is not an internal subgraph"))?;
    let quotient = graph.shrunken_sub_lmb(
        host.subgraph(),
        &given,
        graph.dummy_stripped_external_flows_of(host.subgraph()),
        None,
    )?;
    assert!(
        quotient.loop_edges.is_empty(),
        "contracting the host onto itself leaves no quotient loop: {:?}",
        quotient.loop_edges
    );

    // The whole host is active: both bubbles keep their carriers.
    let frames = [frame(&first), frame(&second)];
    let (canonical, selected) = coordinate_lmb(
        &ctx,
        &host,
        &join,
        Some(host.subgraph()),
        &frames,
        host.subgraph(),
    )?;
    assert_eq!(canonical.loop_edges.raw, vec![EdgeIndex(2), EdgeIndex(4)]);
    assert_eq!(selected.loop_edges.raw, vec![EdgeIndex(2), EdgeIndex(4)]);

    // One bubble was integrated: the active support lies strictly inside
    // the host and keeps only the other bubble's carrier.
    for (active, carrier) in [(&first, EdgeIndex(2)), (&second, EdgeIndex(4))] {
        let (_, selected) = coordinate_lmb(
            &ctx,
            &host,
            &join,
            Some(active.subgraph()),
            &[frame(active)],
            active.subgraph(),
        )?;
        assert_eq!(selected.loop_edges.raw, vec![carrier]);
    }
    Ok(())
}

#[test]
fn rekey_static_zeros_sigma_and_strips_its_energies() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Mixed.topology()?;
    let cff = ThermalCff::new(&mut graph)?;
    let orientation = cff.orientation(true);
    let production = cff.branches(&graph, true)?;
    let edges = |edges: &[usize]| {
        edges
            .iter()
            .copied()
            .map(EdgeIndex)
            .collect::<BTreeSet<_>>()
    };
    // The q class, the static k class and the line e5 they imply.
    let (soft, second_loop) = (edges(&[0, 1]), edges(&[4, 5]));
    let sigma = soft.union(&second_loop).copied().collect::<BTreeSet<_>>();
    let routes_sigma = |map: &[LinearEnergyExpr]| {
        map.iter().enumerate().any(|(slot, energy)| {
            !sigma.contains(&EdgeIndex(slot))
                && energy
                    .internal_terms
                    .iter()
                    .any(|(edge, _)| sigma.contains(edge))
        })
    };
    // The CFF orientation maps put every line on shell. The residue map of
    // an earlier operation may route q0 = OSE(0) through e3 = p - q, and
    // sample a static line at a shifted scale; such source keys compose.
    let orientations = orientation.exact_orientations()?;
    let routed = DirectResidueBranches::from_keyed(production.iter_keys().map(|(key, roots)| {
        let mut map = orientations[key.selector_host].edge_energy_map.clone();
        map[3] = map[3].clone() + LinearEnergyExpr::ose(EdgeIndex(0), -1);
        map[0] = map[0].clone()
            + LinearEnergyExpr {
                uniform_scale_coeff: Rational::from(1),
                constant: Rational::from((1, 2)),
                ..LinearEnergyExpr::zero()
            };
        (
            DirectResidueKey::source(key.selector_host, map),
            roots.clone(),
        )
    }))?;
    let point = Point::new(&RING);
    let sum = evaluate(&graph, &production, &point)?;
    let energy = |edge| GS.emr_mom(EdgeIndex(edge), GS.cind(0));
    let names_sigma = |family: &DirectResidueBranches| {
        family.iter_keys().any(|(_, roots)| {
            roots.numerators().iter().any(|entry| {
                sigma
                    .iter()
                    .any(|edge| entry.rhs.contains(ose_atom_from_index(*edge)))
            })
        })
    };
    for (family, routes_hard_lines) in [(&production, false), (&routed, true)] {
        let rekeyed = family.rekey_static(orientation, &sigma)?;
        for (key, _) in rekeyed.iter_keys() {
            let map = key
                .source_edge_energy_map()
                .expect("re-keyed rows carry their own energy maps");
            for edge in &sigma {
                assert!(
                    map[usize::from(*edge)].is_zero(),
                    "the entry of {edge} must vanish: {:?}",
                    map[usize::from(*edge)]
                );
            }
            assert!(!routes_sigma(map), "a static energy survived: {map:?}");
        }
        assert_eq!(
            rekeyed.rekey_static(orientation, &sigma)?,
            rekeyed,
            "re-keying is idempotent"
        );
        assert_eq!(
            family
                .rekey_static(orientation, &soft)?
                .rekey_static(orientation, &second_loop)?,
            rekeyed,
            "re-keying composes on source maps"
        );
        // Every row has its own selector host, so no keys coincide here.
        assert_eq!(rekeyed.iter_keys().count(), family.iter_keys().count());
        assert_close(
            evaluate(&graph, &rekeyed, &point)?,
            sum,
            1e-12,
            "re-keyed sum",
        );
        // Numerators attached afterwards see the static lines at zero
        // energy, and no static energy in the hard lines.
        for edge in [0, 1, 4, 5] {
            let attached = rekeyed.multiply_key_mapped(
                orientation,
                &graph,
                &energy(edge),
                DirectResidueBranches::numerator_scope(),
            )?;
            assert!(
                attached.factorized_sum().is_zero(),
                "the energy of {edge} must vanish"
            );
        }
        let hard = energy(2) + energy(3);
        for (keyed, expected) in [(family, routes_hard_lines), (&rekeyed, false)] {
            let attached = keyed.multiply_key_mapped(
                orientation,
                &graph,
                &hard,
                DirectResidueBranches::numerator_scope(),
            )?;
            assert_eq!(
                names_sigma(&attached),
                expected,
                "the energies of the hard lines"
            );
        }
    }

    // Rows on one selector host whose maps differ only in a static entry, or
    // in a static energy routed through a hard line, coincide once re-keyed
    // and add.
    let (key, roots) = production
        .iter_keys()
        .next()
        .expect("the family has a residue row");
    let map = orientations[key.selector_host].edge_energy_map.clone();
    let mut sampled = map.clone();
    sampled[4] = sampled[4].clone()
        + LinearEnergyExpr {
            uniform_scale_coeff: Rational::from(2),
            ..LinearEnergyExpr::zero()
        };
    let mut routed_through_hard = map.clone();
    routed_through_hard[2] =
        routed_through_hard[2].clone() + LinearEnergyExpr::ose(EdgeIndex(5), 1);
    let doubled = roots.map(|root| root * Atom::num(2));
    let rows = DirectResidueBranches::from_keyed([
        (
            DirectResidueKey::source(key.selector_host, map),
            roots.clone(),
        ),
        (
            DirectResidueKey::source(key.selector_host, sampled),
            doubled,
        ),
        (
            DirectResidueKey::source(key.selector_host, routed_through_hard),
            roots.clone(),
        ),
    ])?;
    assert_eq!(rows.iter_keys().count(), 3);
    let merged = rows.rekey_static(orientation, &sigma)?;
    let merged_roots = match merged.iter_keys().collect::<Vec<_>>().as_slice() {
        [(_, merged_roots)] => (*merged_roots).clone(),
        keys => {
            return Err(eyre!(
                "the re-keyed rows must coincide: {} keys",
                keys.len()
            ));
        }
    };
    for ((cut, merged_root), (_, root)) in merged_roots.iter().zip(roots.iter()) {
        assert!(
            (merged_root - &(root * Atom::num(4))).expand().is_zero(),
            "the merged root of {cut:?} is the sum of the three rows"
        );
    }
    assert_close(
        evaluate(&graph, &merged, &point)?,
        evaluate(&graph, &rows, &point)?,
        1e-12,
        "merged rows",
    );
    Ok(())
}

#[test]
fn zero_mode_projection_short_circuits_projected_rows() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Mixed.topology()?;
    let cff = ThermalCff::new(&mut graph)?;
    let orientation = cff.orientation(true);
    let branches = cff.branches(&graph, true)?;
    let soft = || class(&graph, &[0, 1], &[4]);
    let static_k = || class(&graph, &[4], &[2]);
    let projected = project(
        &graph,
        &branches,
        vec![soft()?, static_k()?],
        ZeroModeSector::Complete,
    )?;
    let static_a = &projected[&vec![4, 2]];
    let is_zero = |family: &DirectResidueBranches| {
        family
            .iter_keys()
            .all(|(_, roots)| roots.iter().all(|(_, root)| root.is_zero()))
    };
    assert!(!is_zero(static_a), "the static host-A term must not vanish");

    // Re-projecting reads the projected rows off unchanged: no series wraps
    // them again, and they have no other pole.
    for classes in [
        vec![class(&graph, &[0, 1], &[3, 4])?],
        vec![soft()?, static_k()?],
        vec![static_k()?, soft()?],
    ] {
        let orders = classes
            .iter()
            .map(|class| *class.orders.last().expect("an order"))
            .collect::<Vec<_>>();
        let again = project(&graph, static_a, classes, ZeroModeSector::Complete)?;
        assert_eq!(again[&orders], *static_a, "orders {orders:?}");
        if let Some(odd) = again.get(&vec![3]) {
            assert!(is_zero(odd), "a projected row has no [x^-3]");
        }
    }

    // Projected rows never reach the hard-vacuum sector, even when frozen,
    // whereas an unprojected frozen row is an internal error.
    let host_a = graph.edge_subgraph(&[2, 3]);
    let frozen =
        |family: &DirectResidueBranches| family.map(|root| root * soft_branch_scope(&host_a));
    let vacuum_soft = || class(&graph, &[0, 1], &[3]);
    let odd = project(
        &graph,
        &frozen(static_a),
        vec![vacuum_soft()?],
        ZeroModeSector::HardVacuum,
    )?;
    assert!(is_zero(&odd[&vec![3]]));
    let error = ZeroModeProjection::new(vec![vacuum_soft()?], ZeroModeSector::HardVacuum)?
        .coefficients(&graph, &frozen(&branches))
        .expect_err("a frozen row must not be vacuum-projected");
    assert!(
        matches!(error, ZeroModeError::FrozenHardVacuum { .. }),
        "{error}"
    );

    // A frozen numerator body that only projected summands call does not
    // stop the vacuum projection of the unprojected summands of its rows.
    let frozen_body = static_a.multiply_key_mapped(
        orientation,
        &graph,
        &(soft_branch_scope(&host_a) * GS.emr_mom(EdgeIndex(2), GS.cind(0))),
        DirectResidueBranches::numerator_scope(),
    )?;
    assert!(
        frozen_body
            .numerator_definitions()
            .values()
            .any(|entry| entry.rhs.contains_symbol(*LOCAL_3D_SOFT_SCOPE)),
        "the attached body carries the frozen tag"
    );
    let point = Point::new(&RING);
    let hard_vacuum = project(
        &graph,
        &branches,
        vec![vacuum_soft()?],
        ZeroModeSector::HardVacuum,
    )?;
    let beside = project(
        &graph,
        &frozen_body.zip_add(&branches)?,
        vec![vacuum_soft()?],
        ZeroModeSector::HardVacuum,
    )?;
    assert_close(
        evaluate(&graph, &beside[&vec![3]], &point)?,
        evaluate(&graph, &hard_vacuum[&vec![3]], &point)?,
        1e-12,
        "the vacuum-soft term beside a frozen body",
    );

    // Projected summands of the vacuum-soft order are read off as they are:
    // the hard-vacuum sector takes them from an earlier hard-vacuum
    // projection, never from the complete sector.
    let again = project(
        &graph,
        &hard_vacuum[&vec![3]],
        vec![vacuum_soft()?],
        ZeroModeSector::HardVacuum,
    )?;
    assert_close(
        evaluate(&graph, &again[&vec![3]], &point)?,
        evaluate(&graph, &hard_vacuum[&vec![3]], &point)?,
        1e-12,
        "the hard-vacuum term projected again",
    );
    let complete = project(
        &graph,
        &branches,
        vec![vacuum_soft()?],
        ZeroModeSector::Complete,
    )?;
    let error = ZeroModeProjection::new(vec![vacuum_soft()?], ZeroModeSector::HardVacuum)?
        .coefficients(&graph, &complete[&vec![3]])
        .expect_err("the complete sector keeps the hard occupations");
    assert!(
        matches!(error, ZeroModeError::HardVacuumOccupation { .. }),
        "{error}"
    );

    // S_A after U_B vanishes: U_B projects the second loop's occupations to
    // vacuum, and its start routes the class energy OSE(4) away.
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let second = forest_node(&graph, &[4, 5], ApproximationType::MUV, 2)?;
    let started = branches.map_expressions(|atom| {
        let vacuum = Direct3dApproximation::project_unfrozen_occupations(
            &ctx,
            atom,
            second.subgraph(),
            SoftBranches::new(&graph, second.subgraph(), Local3DDeformation::Ordinary),
        )?;
        start(&ctx, &second, &vacuum, &Atom::one(), None, second.lmb())
    })?;
    let started_sum = started.factorized_sum();
    assert!(!started_sum.contains(ose_atom_from_index(EdgeIndex(4))));
    assert!(started_sum.contains_symbol(GS.on_shell_energy));
    // The q zero mode itself survives; only the static k mode is gone.
    let soft_only = project(&graph, &started, vec![soft()?], ZeroModeSector::Complete)?;
    assert!(!is_zero(&soft_only[&vec![4]]));
    let after = project(
        &graph,
        &started,
        vec![soft()?, static_k()?],
        ZeroModeSector::Complete,
    )?;
    assert!(is_zero(&after[&vec![4, 2]]), "S_A after U_B must vanish");
    Ok(())
}

/// A projected summand `c OSE(rep)^-a` keeps its order where a later series
/// carried the class propagators into each term of a sum, as the static U of
/// the gluon bubble does on host A's self-suffix
/// ([`gluon_ring_static_u_on_the_self_suffix_cancels_its_static_mode`]): the
/// complete sector reads `c` off, and the hard-vacuum sector takes nothing
/// from it, also when it is frozen. A sum whose terms carry different powers
/// of the propagator is not of that form and still needs its expansion,
/// which a frozen row must not reach.
#[test]
fn projected_summands_keep_their_order_through_sums() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Mixed.topology()?;
    let cff = ThermalCff::new(&mut graph)?;
    let branches = cff.branches(&graph, true)?;
    let soft = || class(&graph, &[0, 1], &[4]);
    let vacuum_soft = || class(&graph, &[0, 1], &[3]);
    let static_a = project(
        &graph,
        &branches,
        vec![soft()?, class(&graph, &[4], &[2])?],
        ZeroModeSector::Complete,
    )?
    .remove(&vec![4, 2])
    .ok_or_else(|| eyre!("the static host-A term is projected"))?;
    let is_zero = |family: &DirectResidueBranches| {
        family
            .iter_keys()
            .all(|(_, roots)| roots.iter().all(|(_, root)| root.is_zero()))
    };
    assert!(
        !is_zero(&static_a),
        "the static host-A term must not vanish"
    );
    let frozen_scope = soft_branch_scope(&graph.edge_subgraph(&[2, 3]));
    let [a, b] = ["a", "b"].map(|name| Atom::var(symbol!(&format!("zero_mode_test_{name}"))));
    // `(a r + b r) S`: each residue root `r` is one frozen summand whose
    // propagators sit in both terms of a sum.
    let distributed = static_a.map(|root| (root * &a + root * &b) * &frozen_scope);
    assert!(
        !has_top_level_propagators(&distributed, &[(0, 4)]),
        "the class propagators are no top-level factors"
    );
    let read = project(
        &graph,
        &distributed,
        vec![soft()?],
        ZeroModeSector::Complete,
    )?;
    for ((key, expected), (read_key, actual)) in
        static_a.iter_keys().zip(read[&vec![4]].iter_keys())
    {
        assert_eq!(key, read_key);
        for ((_, expected), (_, actual)) in expected.iter().zip(actual.iter()) {
            let difference = actual - expected * (&a + &b) * &frozen_scope;
            assert!(
                difference.expand().is_zero(),
                "the complete sector reads the coefficient off: {difference}"
            );
        }
    }
    let odd = project(
        &graph,
        &distributed,
        vec![vacuum_soft()?],
        ZeroModeSector::HardVacuum,
    )?;
    assert!(is_zero(&odd[&vec![3]]), "a projected row has no [x^-3]");

    let propagator = ose_atom_from_index(EdgeIndex(0));
    let mixed = static_a.map(|root| (root * &a + root * &propagator * &b) * &frozen_scope);
    let error = ZeroModeProjection::new(vec![vacuum_soft()?], ZeroModeSector::HardVacuum)?
        .coefficients(&graph, &mixed)
        .expect_err("a frozen sum of two orders needs its expansion");
    assert!(
        matches!(error, ZeroModeError::FrozenHardVacuum { .. }),
        "{error}"
    );
    Ok(())
}

#[test]
fn mixed_ring_sigma_occupations_cancel_identically_in_s0() -> Result<()> {
    assert_sigma_occupations_cancel(Ring::Mixed)
}

#[test]
fn ghost_ring_sigma_occupations_cancel_identically_in_s0() -> Result<()> {
    assert_sigma_occupations_cancel(Ring::Ghost)
}

/// Host A's static term `[x_q^-4 x_4^-2]` keeps occupations of the line e5
/// that q and k make static. They cancel identically in its sum: neither
/// their vacuum projection nor another temperature of theirs changes it,
/// with or without vacuum subtraction. A second loop vertex numerator
/// `E4 E5 + c^2`, of degree two in that loop's energies like the b0/b1
/// vertices, takes its static value `c^2`.
fn assert_sigma_occupations_cancel(ring: Ring) -> Result<()> {
    test_initialise()?;
    let [q, _, e2, e3, e4, e5] = RING;
    let sigma = [0, 1, 4, 5].map(EdgeIndex);
    let vertex = Atom::var(symbol!("zero_mode_test_vertex_scale"));
    let numerator =
        GS.emr_mom(EdgeIndex(4), GS.cind(0)) * GS.emr_mom(EdgeIndex(5), GS.cind(0)) + vertex.pow(2);
    let c = 0.83;
    let sigma_beta = Atom::var(symbol!("zero_mode_test_sigma_beta"));
    let mut retained = false;
    for vacuum_subtraction in [false, true] {
        let mut graph = ring.topology()?;
        let cff = ThermalCff::new(&mut graph)?;
        let orientation = cff.orientation(vacuum_subtraction);
        let branches = cff
            .branches(&graph, vacuum_subtraction)?
            .multiply_key_mapped(
                orientation,
                &graph,
                &numerator,
                DirectResidueBranches::numerator_scope(),
            )?;
        let normalization = normalization(&graph)?;
        let projected = project(
            &graph,
            &branches,
            vec![class(&graph, &[0, 1], &[4])?, class(&graph, &[4], &[2])?],
            ZeroModeSector::Complete,
        )?;
        let static_a = &projected[&vec![4, 2]];
        let is_sigma = |view: AtomView<'_>| {
            ThermalDistributionCall::try_from(view).is_ok_and(|call| sigma.contains(&call.edge))
        };
        let names_sigma = |family: &DirectResidueBranches| {
            let mut found = false;
            family.factorized_sum().as_view().visitor(&mut |view| {
                found |= is_sigma(view);
                !found
            });
            found
        };
        retained |= names_sigma(static_a);
        let sigma_medium = |medium| {
            static_a.map_expressions(|atom| {
                graph.make_selected_thermal_distributions_explicit(atom, medium, |call| {
                    Ok(sigma.contains(&call.edge))
                })
            })
        };
        let vacuum = sigma_medium(MediumMode::Vacuum)?;
        assert!(!names_sigma(&vacuum));
        // The projection leaves no other explicit temperature to rename.
        assert!(
            !static_a
                .factorized_sum()
                .contains_symbol(GS.inverse_temperature)
        );
        let tilted =
            sigma_medium(MediumMode::ThermodynamicEquilibrium)?.map_expressions(|atom| {
                Ok(atom
                    .replace(GS.inverse_temperature)
                    .with(sigma_beta.clone()))
            })?;
        assert!(!names_sigma(&tilted));
        for mu in [0.0, QUARK_MU] {
            let context = format!("{ring:?}, vacuum subtraction {vacuum_subtraction}, mu {mu}");
            let point = Point::new(&RING).with_mu(mu).with(vertex.clone(), c);
            let complete = evaluate(&graph, static_a, &point)? / normalization;
            // T^2 c^2 Pi_A(0) / (q^4 E4^2 E5^2) at T = 1.
            let expected = c * c * static_bubble(e2, e3, |e| fermi_tadpole(1.0, e, mu))
                / (q.powi(4) * e4 * e4 * e5 * e5);
            assert_close(complete, expected, 1e-9, &context);
            assert_close(
                evaluate(&graph, &vacuum, &point)? / normalization,
                complete,
                1e-10,
                &context,
            );
            for beta in [1.0, 0.4] {
                let point = Point::new(&RING)
                    .with_mu(mu)
                    .with(vertex.clone(), c)
                    .with(sigma_beta.clone(), beta);
                assert_close(
                    evaluate(&graph, &tilted, &point)? / normalization,
                    complete,
                    1e-10,
                    &format!("{context}, static-line beta {beta}"),
                );
            }
        }
    }
    assert!(
        retained,
        "{ring:?}: the static term must keep an occupation of a static line"
    );
    Ok(())
}

#[test]
fn zero_mode_temperatures_belong_to_their_classes() -> Result<()> {
    test_initialise()?;
    let [q, _, e2, e3, e4, e5] = RING;
    let mut graph = Ring::Ghost.topology()?;
    let (branches, normalization) = thermal_sector(&mut graph, false)?;
    let reps = [0, 4, 5];
    let values = reps.map(|rep| Atom::var(symbol!(&format!("zero_mode_test_temperature_{rep}"))));
    let point = |temperatures: [f64; 3]| {
        values
            .iter()
            .zip(temperatures)
            .fold(Point::new(&RING), |point, (value, temperature)| {
                point.with(value.clone(), temperature)
            })
    };
    let points = [
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 1.0],
        [1.0, 3.0, 1.0],
        [2.0, 3.0, 1.0],
        [1.0, 1.0, 5.0],
    ]
    .map(point);
    // T^2 Pi_A(0) / (q^4 E4^2 E5^2): one temperature per independent class.
    let expected =
        static_bubble(e2, e3, |e| fermi_tadpole(1.0, e, 0.0)) / (q.powi(4) * e4 * e4 * e5 * e5);
    for classes in [
        vec![class(&graph, &[0, 1], &[4])?, class(&graph, &[4], &[2])?],
        vec![
            class(&graph, &[0, 1], &[4])?,
            class(&graph, &[4], &[2])?,
            class(&graph, &[5], &[2])?,
        ],
    ] {
        let context = format!("{classes:?}");
        let projected = project(&graph, &branches, classes, ZeroModeSector::Complete)?;
        let (_, coefficient) = projected
            .iter()
            .next()
            .expect("one order tuple was requested");
        let detached =
            reps.iter()
                .zip(&values)
                .fold(coefficient.factorized_sum(), |atom, (rep, value)| {
                    atom.replace(class_zero_mode_temperature(EdgeIndex(*rep)))
                        .with(value.clone())
                });
        assert!(
            !detached.contains_symbol(*ZERO_MODE_TEMPERATURE),
            "{context}: only class representatives own zero-mode temperatures"
        );
        let tilted = evaluate_points(&graph, &detached, &points)?;
        assert_close(tilted[0] / normalization, expected, 1e-9, &context);
        // Linear in T_zm(0) and T_zm(4); the implied static line e5 adds none.
        for (index, ratio) in [(1, 2.0), (2, 3.0), (3, 6.0), (4, 1.0)] {
            assert_close(tilted[index] / tilted[0], ratio, 1e-12, &context);
        }
    }
    Ok(())
}

/// The Feynman-rule numerator that the edges of `subgraph` add to `given`,
/// with its colour algebra closed as the kernel attaches it, for the
/// massless quark of the `sm-thermal` restriction.
fn physical_numerator(graph: &Graph, subgraph: &SuBitGraph, given: &SuBitGraph) -> Atom {
    graph
        .color_reduced_numerator(subgraph, given)
        .expect("the graph numerator is available")
        .replace(graph[EdgeIndex(2)].mass_atom())
        .with(Atom::Zero)
}

/// Close the colour, spinor and Lorentz algebra of a closed physical
/// numerator network and write its scalar products in components
/// `Q(e, cind(mu))`. The trace is expanded before its metrics contract.
fn contract_numerators(atom: &Atom) -> Result<Atom> {
    Ok(atom
        .simplify_color_with(ColorSimplifySettings::default().with_cof_dimension_invariants())
        .simplify_gamma()
        .replace(GS.dim)
        .with(Atom::num(4))
        .expand()
        .simplify_metrics()
        .expand_dots()?
        .expand())
}

/// Every closed numerator body of `branches` contracted, then the summed
/// rows with their bodies and scalar definitions resolved.
fn contracted_sum(branches: &DirectResidueBranches) -> Result<Atom> {
    Ok(branches
        .map_numerators(|entry| contract_numerators(&entry.rhs))?
        .resolved_scalars()?
        .factorized_sum())
}

/// The bodies of a family in which every row calls host A's numerator and
/// the argument-free numerator of the b0/b1 vertices: re-keyed rows give
/// every static energy the same value, zero. Each body of host A, including
/// the copies a projection renames into its own scope, absorbs the
/// contracted product, and the vertex numerator becomes one. `definitions`
/// holds one definition per family.
fn absorb_static_vertices<'a>(
    definitions: impl IntoIterator<Item = &'a FnMapEntry>,
) -> Result<BTreeMap<Atom, Atom>> {
    let (vertices, hosts): (Vec<_>, Vec<_>) = definitions
        .into_iter()
        .partition(|entry| entry.args.is_empty());
    let [vertices] = vertices.as_slice() else {
        return Err(eyre!(
            "expected one argument-free vertex numerator, got {}",
            vertices.len()
        ));
    };
    if hosts.is_empty() {
        return Err(eyre!("expected host A's numerator"));
    }
    let mut bodies = BTreeMap::from([(vertices.rhs.clone(), Atom::one())]);
    for host in hosts {
        bodies.insert(
            host.rhs.clone(),
            contract_numerators(&(&host.rhs * &vertices.rhs))?,
        );
    }
    Ok(bodies)
}

/// The summed rows of a re-keyed static host-A term with its b0/b1
/// numerator attached, contracted through [`absorb_static_vertices`].
fn contracted_static_term(branches: &DirectResidueBranches) -> Result<Atom> {
    let bodies = absorb_static_vertices(branches.numerator_definitions().into_values())?;
    Ok(branches
        .map_numerators(|entry| absorbed_body(&bodies, &entry.rhs))?
        .resolved_scalars()?
        .factorized_sum())
}

/// The body that [`absorb_static_vertices`] gives `body`.
fn absorbed_body(bodies: &BTreeMap<Atom, Atom>, body: &Atom) -> Result<Atom> {
    bodies
        .get(body)
        .cloned()
        .ok_or_else(|| eyre!("an unexpected numerator body {body}"))
}

/// The model parameters `UFO::*` of an atom at unit value: the couplings of
/// the physical numerators, common to both sides of every comparison.
fn unit_couplings(atom: &Atom) -> Vec<(Atom, f64)> {
    let mut couplings = atom
        .get_all_symbols(false)
        .into_iter()
        .filter(|symbol| symbol.get_namespace() == "UFO")
        .map(|symbol| (symbol.get_name().to_owned(), Atom::var(symbol)))
        .collect::<Vec<_>>();
    couplings.sort();
    couplings
        .into_iter()
        .map(|(_, coupling)| (coupling, 1.0))
        .collect()
}

fn assert_close_complex(
    actual: Complex<f64>,
    expected: Complex<f64>,
    tolerance: f64,
    context: &str,
) {
    let difference = (actual.re - expected.re).hypot(actual.im - expected.im);
    let scale = expected.re.hypot(expected.im);
    assert!(
        difference <= tolerance * scale && scale > 0.0,
        "{context}: {actual:?} != {expected:?} (relative {})",
        difference / scale
    );
}

/// Spatial loop momenta `[q, p, k]` of the physical ring points.
const LOOP_MOMENTA: [[[f64; 3]; 3]; 3] = [
    [
        [0.31, -0.22, 0.17],
        [0.52, 0.41, -0.63],
        [-0.44, 0.73, 0.28],
    ],
    [
        [-0.27, 0.35, 0.12],
        [0.18, -0.66, 0.47],
        [0.61, 0.09, -0.38],
    ],
    [
        [0.12, 0.29, -0.41],
        [-0.71, 0.24, 0.33],
        [0.36, -0.52, -0.19],
    ],
];

/// One point of the physical rings: the spatial momenta of e0..e5,
/// `q, q, p, p - q, k, k - q`, and their massless on-shell energies.
struct RingKinematics {
    momenta: [[f64; 3]; 6],
    energies: [f64; 6],
}

impl RingKinematics {
    fn new([q, p, k]: [[f64; 3]; 3]) -> Self {
        let difference = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let momenta = [q, q, p, difference(p, q), k, difference(k, q)];
        let energies = momenta.map(|momentum| momentum.iter().map(|c| c * c).sum::<f64>().sqrt());
        Self { momenta, energies }
    }

    /// The spatial components `Q(e, cind(i))` and their values.
    fn components(&self) -> Vec<(Atom, f64)> {
        self.momenta
            .iter()
            .enumerate()
            .flat_map(|(edge, momentum)| {
                momentum.iter().enumerate().map(move |(axis, value)| {
                    (GS.emr_mom(EdgeIndex(edge), GS.cind(axis + 1)), *value)
                })
            })
            .collect()
    }

    /// The evaluation point at chemical potential `mu`, with the spatial
    /// components and the `parameters` of the evaluated atom.
    fn point(&self, mu: f64, parameters: &[(Atom, f64)]) -> Point<'_> {
        self.point_with_energies(&self.energies, mu, parameters)
    }

    /// [`Self::point`] at other on-shell `energies` of the same momenta.
    fn point_with_energies<'a>(
        &self,
        energies: &'a [f64],
        mu: f64,
        parameters: &[(Atom, f64)],
    ) -> Point<'a> {
        self.components()
            .into_iter()
            .chain(parameters.iter().cloned())
            .fold(Point::new(energies).with_mu(mu), |point, (atom, value)| {
                point.with(atom, value)
            })
    }

    /// `T^2 Pi_A(0) N_0 / (q^4 k^2 (k - q)^2)` at `T = 1`, in the
    /// normalization of a residue root after its factor `i^3`: the
    /// contracted 4D `numerator` of the ring at the static modes
    /// `q0 = k0 = 0`, its fermionic sum over `p0 = (p - q)0 = i w_n + mu`
    /// taken by the residues of `n_F(z - mu)` at `z = +-E2, +-E3`, over the
    /// static propagators and the measure `(2 pi)^9`. Independent of the CFF.
    fn static_host_a(&self, graph: &Graph, numerator: &Atom, mu: f64) -> Result<Complex<f64>> {
        let couplings = unit_couplings(numerator);
        let energies = (0..6).map(|edge| GS.emr_mom(EdgeIndex(edge), GS.cind(0)));
        let components = self.components();
        let parameters = energies
            .chain(components.iter().map(|(atom, _)| atom.clone()))
            .chain(couplings.iter().map(|(atom, _)| atom.clone()))
            .collect::<Vec<_>>();
        let mut evaluator = numerator
            .evaluator(&parameters)
            .function_map(graph.param_builder.fn_map.clone())
            .build()
            .map_err(|error| eyre!("{error}"))?
            .map_coeff(&|coefficient| {
                Complex::new(coefficient.re.to_f64(), coefficient.im.to_f64())
            });
        let mut at = |z: f64| {
            let values = [0.0, 0.0, z, z, 0.0, 0.0]
                .into_iter()
                .chain(components.iter().map(|(_, value)| *value))
                .chain(couplings.iter().map(|(_, value)| *value))
                .map(|value| Complex::new(value, 0.0))
                .collect::<Vec<_>>();
            evaluator.evaluate_single(&values)
        };
        let [e0, e1, e2, e3, e4, e5] = self.energies;
        let (mut real, mut imaginary) = (0.0, 0.0);
        for (pole, other) in [(e2, e3), (e3, e2)] {
            for z in [pole, -pole] {
                let weight = n_fermi(1.0, z - mu) / (2.0 * z * (pole * pole - other * other));
                let value = at(z);
                real += weight * value.re;
                imaginary += weight * value.im;
            }
        }
        let scale = (2.0 * PI).powi(9) * e0 * e0 * e1 * e1 * e4 * e4 * e5 * e5;
        Ok(Complex::new(real / scale, imaginary / scale))
    }

    /// A contracted 4D `numerator` without the quark loop's momenta at the
    /// static modes, every energy zero, with the spatial components and the
    /// further `parameters` of the point.
    fn static_numerator(
        &self,
        graph: &Graph,
        numerator: &Atom,
        parameters: &[(Atom, f64)],
    ) -> Result<Complex<f64>> {
        let values = (0..6)
            .map(|edge| (GS.emr_mom(EdgeIndex(edge), GS.cind(0)), 0.0))
            .chain(self.components())
            .chain(parameters.iter().cloned())
            .collect::<Vec<_>>();
        let symbols = values
            .iter()
            .map(|(symbol, _)| symbol.clone())
            .collect::<Vec<_>>();
        let mut evaluator = numerator
            .evaluator(&symbols)
            .function_map(graph.param_builder.fn_map.clone())
            .build()
            .map_err(|error| eyre!("{error}"))?
            .map_coeff(&|coefficient| {
                Complex::new(coefficient.re.to_f64(), coefficient.im.to_f64())
            });
        let values = values
            .iter()
            .map(|(_, value)| Complex::new(*value, 0.0))
            .collect::<Vec<_>>();
        Ok(evaluator.evaluate_single(&values))
    }
}

/// Host A's numerator, attached before the projection as the kernel does,
/// and the b0/b1 vertex numerator of its cograph.
struct HostANumerators {
    host: Atom,
    vertices: Atom,
}

impl HostANumerators {
    fn new(graph: &Graph) -> Self {
        let host = graph.edge_subgraph(&[2, 3]);
        let empty = graph.empty_subgraph::<SuBitGraph>();
        let cograph = graph.full_filter().subtract(&host);
        Self {
            host: physical_numerator(graph, &host, &empty),
            vertices: physical_numerator(graph, &cograph, &host),
        }
    }

    /// The static term `[x_q^-4 x_4^-2]` of the root with host A's
    /// numerator, for the classes q = {e0, e1} and k = {e4}, still keyed by
    /// the production maps.
    fn projected_static_term(
        &self,
        graph: &Graph,
        orientation: OrientationProjection<'_>,
        root: &DirectResidueBranches,
    ) -> Result<DirectResidueBranches> {
        let with_host = root.multiply_key_mapped(
            orientation,
            graph,
            &self.host,
            DirectResidueBranches::numerator_scope(),
        )?;
        let mut projected = project(
            graph,
            &with_host,
            vec![class(graph, &[0, 1], &[4])?, class(graph, &[4], &[2])?],
            ZeroModeSector::Complete,
        )?;
        projected
            .remove(&vec![4, 2])
            .ok_or_else(|| eyre!("the projection returns the requested order tuple"))
    }
}

/// Sigma of host A on the rings: the classes q = {e0, e1}, k = {e4} and the
/// line e5 = k - q.
fn static_lines() -> BTreeSet<EdgeIndex> {
    [0, 1, 4, 5].map(EdgeIndex).into_iter().collect()
}

#[test]
fn mixed_ring_host_a_rekey_matches_the_zero_mode_on_production_roots() -> Result<()> {
    assert_host_a_rekey_matches_the_zero_mode(Ring::Mixed)
}

#[test]
fn ghost_ring_host_a_rekey_matches_the_zero_mode_on_production_roots() -> Result<()> {
    assert_host_a_rekey_matches_the_zero_mode(Ring::Ghost)
}

/// Host A's static term on the production root, with only A's numerator
/// projected and the b0/b1 numerator attached after re-keying, is
/// `T^2 Pi_A(0) N_0(k, q) / (q^4 k^2 (k - q)^2)` with both numerators at
/// zero static energies: neither the contact rows of the energy-degree
/// bounds nor the key maps leave a static energy in the vertices.
fn assert_host_a_rekey_matches_the_zero_mode(ring: Ring) -> Result<()> {
    test_initialise()?;
    let mut graph = ring.physical_graph()?;
    let numerators = HostANumerators::new(&graph);
    let reference = contract_numerators(&(&numerators.host * &numerators.vertices))?;
    for vacuum_subtraction in [true, false] {
        let cff = ThermalCff::production(&mut graph, vacuum_subtraction)?;
        let orientation = cff.orientation(vacuum_subtraction);
        let root = cff.branches(&graph, vacuum_subtraction)?;
        let projected = numerators.projected_static_term(&graph, orientation, &root)?;
        let attach = |branches: &DirectResidueBranches| {
            branches.multiply_key_mapped(
                orientation,
                &graph,
                &numerators.vertices,
                DirectResidueBranches::numerator_scope(),
            )
        };
        // The production maps would give the vertex numerator the on-shell
        // static energies of each row.
        assert!(
            attach(&projected)?
                .iter_keys()
                .flat_map(|(_, roots)| roots.numerators())
                .all(|entry| !entry.args.is_empty()),
            "{ring:?}: production keys hand the vertices their static energies"
        );
        // Re-keyed, every row gives them zero: absorbing the vertex family
        // into host A's requires it to have no row arguments.
        let actual = contracted_static_term(&attach(
            &projected.rekey_static(orientation, &static_lines())?,
        )?)?;
        let couplings = unit_couplings(&actual);
        let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
        let chemical_potentials = [0.0, QUARK_MU];
        let points = kinematics
            .iter()
            .flat_map(|point| {
                chemical_potentials
                    .iter()
                    .map(|mu| point.point(*mu, &couplings))
            })
            .collect::<Vec<_>>();
        let values = evaluate_points_complex(&graph, &actual, &points)?;
        let mut values = values.into_iter();
        for (index, point) in kinematics.iter().enumerate() {
            for mu in chemical_potentials {
                let context = format!(
                    "{ring:?}, vacuum subtraction {vacuum_subtraction}, point {index}, mu {mu}"
                );
                let actual = values.next().expect("one value per point");
                let expected = point.static_host_a(&graph, &reference, mu)?;
                assert_close_complex(actual, expected, 1e-10, &context);
            }
        }
    }
    Ok(())
}

/// The residue-map selectors `sigma(id)` an atom names.
fn selector_ids(atom: &Atom) -> BTreeSet<Atom> {
    let selector = OrientationID::symbol();
    let mut ids = BTreeSet::new();
    atom.as_view().visitor(&mut |view| {
        if let AtomView::Fun(call) = view
            && call.get_symbol() == selector
        {
            ids.insert(view.to_owned());
            return false;
        }
        true
    });
    ids
}

#[test]
fn mixed_ring_source_keys_survive_build_direct_and_materialize() -> Result<()> {
    assert_source_keys_survive_build_direct_and_materialize(Ring::Mixed)
}

#[test]
fn ghost_ring_source_keys_survive_build_direct_and_materialize() -> Result<()> {
    assert_source_keys_survive_build_direct_and_materialize(Ring::Ghost)
}

/// Host A's re-keyed static term, as the final integrand of its forest
/// node: `build_direct` attaches the cograph numerator through the source
/// maps and `materialize` selects every row by its production host. The
/// result equals the direct evaluation of the same branches and the
/// independent static reference.
fn assert_source_keys_survive_build_direct_and_materialize(ring: Ring) -> Result<()> {
    test_initialise()?;
    let mut graph = ring.physical_graph()?;
    let numerators = HostANumerators::new(&graph);
    let reference =
        contract_numerators(&(&numerators.host * &numerators.vertices * graph.global_atom()))?;
    let cff = ThermalCff::production(&mut graph, true)?;
    let orientation = cff.orientation(true);
    let root = cff.branches(&graph, true)?;
    let projected = numerators.projected_static_term(&graph, orientation, &root)?;
    let sector = projected.rekey_static(orientation, &static_lines())?;

    // Every row keeps its own source map on its production selector host.
    let hosts = sector
        .iter_keys()
        .map(|(key, _)| {
            assert!(
                key.source_edge_energy_map().is_some(),
                "a re-keyed row carries its source map"
            );
            key.selector_host.atom()
        })
        .collect::<BTreeSet<_>>();
    // Materializing selects each row by its host, and the explicit sum
    // replaces every selector by one.
    let selected = sector.materialize(true)?;
    let summed = sector.materialize(false)?;
    for ((cut, selected), (_, summed)) in selected.iter().zip(summed.iter()) {
        assert_eq!(selector_ids(selected), hosts, "selectors of {cut:?}");
        let explicit = selected
            .replace(function!(OrientationID::symbol(), W_.a_))
            .with(Atom::one());
        assert!(
            explicit == *summed || (&explicit - summed).expand().is_zero(),
            "the explicit sum of {cut:?}"
        );
    }

    // The final integrand of host A's node attaches the cograph numerator
    // through the source maps: the static energies vanish on every row, so
    // the vertex family has no row arguments; production keys would give
    // it the class energies.
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(&cutset, orientation);
    let settings = local_3d_settings();
    let host = forest_node(&graph, &[2, 3], ApproximationType::IR, 2)?;
    let mut build = |branches: &DirectResidueBranches| {
        FinalIntegrandBuilder::new(localizer, &settings)
            .build_direct(
                &mut graph,
                &host,
                &Direct3dCts::Root(branches.clone()),
                &IntegratedCts::root(),
            )
            .map(FinalIntegrands::into_integrands)
    };
    let production_keyed = build(&projected)?;
    assert!(
        production_keyed
            .numerators()
            .iter()
            .all(|entry| !entry.args.is_empty()),
        "{ring:?}: production keys hand the vertices their static energies"
    );
    let final_integrand = build(&sector)?;
    let bodies = absorb_static_vertices(final_integrand.numerators().iter().map(Arc::as_ref))?;
    let final_sum = final_integrand
        .map_numerators(|body| absorbed_body(&bodies, body))?
        .resolved()?
        .iter()
        .map(|(_, atom)| atom)
        .sum::<Atom>();
    assert_eq!(selector_ids(&final_sum), hosts);
    let final_sum = final_sum
        .replace(function!(OrientationID::symbol(), W_.a_))
        .with(Atom::one());

    // The direct evaluation of the same branches with the same cograph
    // numerator, and the independent static reference.
    let cograph_numerator = &numerators.vertices * graph.global_atom();
    let direct = contracted_static_term(&sector.multiply_key_mapped(
        orientation,
        &graph,
        &cograph_numerator,
        DirectResidueBranches::numerator_scope(),
    )?)?;
    let couplings = unit_couplings(&direct);
    for (index, kinematics) in LOOP_MOMENTA.map(RingKinematics::new).iter().enumerate() {
        for mu in [0.0, QUARK_MU] {
            let context = format!("{ring:?}, point {index}, mu {mu}");
            let point = kinematics.point(mu, &couplings);
            let [final_value, direct_value] = [&final_sum, &direct].map(|atom| {
                evaluate_points_complex(&graph, atom, slice::from_ref(&point))
                    .map(|values| values[0])
            });
            let (final_value, direct_value) = (final_value?, direct_value?);
            assert_close_complex(final_value, direct_value, 1e-12, &context);
            let expected = kinematics.static_host_a(&graph, &reference, mu)?;
            assert_close_complex(final_value, expected, 1e-10, &context);
        }
    }
    Ok(())
}

/// The value at zero of the polynomial through `(nodes[i], values[i])`, by
/// Neville's scheme.
fn extrapolate_to_zero(nodes: &[f64], values: &[f64]) -> f64 {
    let mut table = values.to_vec();
    for level in 1..nodes.len() {
        for index in 0..nodes.len() - level {
            let (near, far) = (nodes[index], nodes[index + level]);
            table[index] = (far * table[index] - near * table[index + 1]) / (far - near);
        }
    }
    table[0]
}

/// Richardson nodes of the class energies `x_q` and `x_4`.
const RICHARDSON_NODES: [f64; 5] = [0.1, 0.05, 0.025, 0.0125, 0.00625];

#[test]
fn mixed_ring_class_contacts_keep_the_zero_mode() -> Result<()> {
    assert_class_contacts_keep_the_zero_mode(Ring::Mixed)
}

#[test]
fn ghost_ring_class_contacts_keep_the_zero_mode() -> Result<()> {
    assert_class_contacts_keep_the_zero_mode(Ring::Ghost)
}

/// The complete physical numerator, attached to the production root with
/// its maps, depends polynomially on the class energies, also through the
/// contact rows of its energy-degree bounds. `[x_q^-4 x_4^-2]` of that root
/// is the joint limit of `x_q^4 x_4^2 F`, extrapolated in both energies,
/// and the independent static host-A term.
fn assert_class_contacts_keep_the_zero_mode(ring: Ring) -> Result<()> {
    test_initialise()?;
    let mut graph = ring.physical_graph()?;
    let numerators = HostANumerators::new(&graph);
    let reference = contract_numerators(&(&numerators.host * &numerators.vertices))?;
    let complete = physical_numerator(
        &graph,
        &graph.full_filter(),
        &graph.empty_subgraph::<SuBitGraph>(),
    );
    let cff = ThermalCff::production(&mut graph, true)?;
    let orientation = cff.orientation(true);
    let integrand = cff.branches(&graph, true)?.multiply_key_mapped(
        orientation,
        &graph,
        &complete,
        DirectResidueBranches::numerator_scope(),
    )?;
    for &rep in ring.numerator_class_energies() {
        assert!(
            integrand.iter_keys().any(|(_, roots)| roots
                .numerators()
                .iter()
                .any(|entry| entry.rhs.contains(ose_atom_from_index(EdgeIndex(rep))))),
            "{ring:?}: the numerator must depend on the class energy of e{rep}"
        );
    }
    let projected = project(
        &graph,
        &integrand,
        vec![class(&graph, &[0, 1], &[4])?, class(&graph, &[4], &[2])?],
        ZeroModeSector::Complete,
    )?;
    let coefficient = contracted_sum(&projected[&vec![4, 2]])?;
    let complete_sum = contracted_sum(&integrand)?;
    let couplings = unit_couplings(&complete_sum);
    let kinematics = RingKinematics::new(LOOP_MOMENTA[0]);
    let [e0, _, _, _, e4, _] = kinematics.energies;
    let grid = RICHARDSON_NODES
        .iter()
        .flat_map(|x_q| {
            RICHARDSON_NODES.iter().map(|x_4| {
                let mut energies = kinematics.energies;
                energies[0] = *x_q;
                energies[1] = *x_q;
                energies[4] = *x_4;
                energies
            })
        })
        .collect::<Vec<_>>();
    for mu in [0.0, QUARK_MU] {
        let context = format!("{ring:?}, mu {mu}");
        let point = kinematics.point(mu, &couplings);
        let static_value =
            evaluate_points_complex(&graph, &coefficient, slice::from_ref(&point))?[0];
        let expected = kinematics.static_host_a(&graph, &reference, mu)?;
        assert_close_complex(static_value, expected, 1e-10, &context);

        // x_q^4 x_4^2 F on the grid, extrapolated to x_4 = 0 for each x_q
        // and then to x_q = 0.
        let points = grid
            .iter()
            .map(|energies| kinematics.point_with_energies(energies, mu, &couplings))
            .collect::<Vec<_>>();
        let scaled = evaluate_points_complex(&graph, &complete_sum, &points)?
            .into_iter()
            .zip(&grid)
            .map(|(value, energies)| {
                let scale = energies[0].powi(4) * energies[4].powi(2);
                (value.re * scale, value.im * scale)
            })
            .collect::<Vec<_>>();
        let nodes = RICHARDSON_NODES.len();
        let limit = |levels: usize, part: fn(&(f64, f64)) -> f64| {
            let per_x_q = scaled
                .chunks(nodes)
                .map(|row| {
                    let row = row.iter().map(part).collect::<Vec<_>>();
                    extrapolate_to_zero(&RICHARDSON_NODES[..levels], &row[..levels])
                })
                .collect::<Vec<_>>();
            extrapolate_to_zero(&RICHARDSON_NODES[..levels], &per_x_q[..levels])
        };
        let extrapolated = Complex::new(limit(nodes, |v| v.0), limit(nodes, |v| v.1));
        let coarser = Complex::new(limit(nodes - 1, |v| v.0), limit(nodes - 1, |v| v.1));
        let projected_value = Complex::new(
            static_value.re * e0.powi(4) * e4.powi(2),
            static_value.im * e0.powi(4) * e4.powi(2),
        );
        assert_close_complex(
            projected_value,
            extrapolated,
            1e-8,
            &format!("{context}: Richardson, one level fewer {coarser:?}"),
        );
    }
    Ok(())
}

/// The normalization of a local operation's output that numerical
/// evaluation needs besides the final values of its tags, which
/// [`evaluate_points`] sets: one UV mass, explicit on-shell energies and the
/// massless quark.
fn normalize_local_operation(graph: &Graph, atom: &Atom) -> Atom {
    atom.replace(function!(GS.m_uv_vacuum, W_.a_))
        .with(Atom::var(GS.m_uv_vacuum))
        .replace(function!(GS.m_uv_expansion, W_.a_))
        .with(Atom::var(GS.m_uv_vacuum))
        .replace(GS.m_uv_expansion)
        .with(Atom::var(GS.m_uv_vacuum))
        .replace(GS.m_ir_expansion)
        .with(Atom::var(GS.m_ir_vacuum))
        .replace(function!(GS.on_shell_energy, W_.a_, W_.prop_))
        .with(Atom::var(W_.prop_).sqrt())
        .replace(graph[EdgeIndex(2)].mass_atom())
        .with(Atom::Zero)
}

#[test]
fn u_a_outputs_keep_the_zero_mode() -> Result<()> {
    test_initialise()?;
    let mut graph = Ring::Quark.topology()?;
    let cff = ThermalCff::new(&mut graph)?;
    let orientation = cff.orientation(true);
    // Degree two in host A's energies, like its quark trace; the factor
    // E0 + c puts the class energy into the numerator body that U_A keeps.
    let energy = |edge| GS.emr_mom(EdgeIndex(edge), GS.cind(0));
    let vertex = Atom::var(symbol!("zero_mode_test_vertex_scale"));
    let numerator = (energy(2) * energy(3) + vertex.pow(2)) * (energy(0) + &vertex);
    let integrand = cff.branches(&graph, true)?.multiply_key_mapped(
        orientation,
        &graph,
        &numerator,
        DirectResidueBranches::numerator_scope(),
    )?;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let host = forest_node(&graph, &[2, 3], ApproximationType::MUV, 2)?;
    let root = empty_forest_node(&graph);
    let u_a = |branches: &DirectResidueBranches| {
        apply_taylor(&ctx, orientation, &host, &root, None, host.lmb(), branches)
    };
    let orders = [3, 4];
    let zero_mode = |branches: &DirectResidueBranches| {
        project(
            &graph,
            branches,
            vec![class(&graph, &[0, 1], &orders)?],
            ZeroModeSector::Complete,
        )
    };
    let u_output = u_a(&integrand)?;
    assert!(
        u_output.iter_keys().any(|(_, roots)| roots
            .numerators()
            .iter()
            .any(|entry| entry.rhs.contains(ose_atom_from_index(EdgeIndex(0))))),
        "the U_A output keeps the class energy in its numerator bodies"
    );
    let projected_u = zero_mode(&u_output)?;
    let u_projected = zero_mode(&integrand)?
        .into_iter()
        .map(|(order, coefficient)| Ok((order, u_a(&coefficient)?)))
        .collect::<Result<BTreeMap<_, _>>>()?;
    let kinematics = RingKinematics::new(LOOP_MOMENTA[0]);
    let parameters = [(vertex.clone(), 0.83), (Atom::var(GS.m_uv_vacuum), 1.3)];
    let points = [0.0, QUARK_MU].map(|mu| kinematics.point(mu, &parameters));
    for order in orders {
        let [after, before] = [&projected_u, &u_projected].map(|family| {
            let atom = family[&vec![order]].resolved_scalars()?.factorized_sum();
            evaluate_points(&graph, &normalize_local_operation(&graph, &atom), &points)
        });
        let (after, before) = (after?, before?);
        for (index, mu) in [0.0, QUARK_MU].iter().enumerate() {
            let context = format!("[x^-{order}], mu {mu}");
            assert!(before[index].abs() > 0.0, "{context}: U_A P_q F vanishes");
            assert_close(after[index], before[index], 1e-12, &context);
        }
    }
    Ok(())
}

/// The uniform numerator sampling scale `M` on the sunset:
/// the static line e0 and the host line e1 carry numerators of degree two,
/// so the sampled CFF keeps `M` in the energy-map entries of both. The zero
/// mode of e0 with the host numerator attached, re-keyed, with the static
/// vertex numerator attached afterwards, is the static term at every `M`.
#[test]
fn sampled_host_rows_keep_the_static_zero_mode_after_rekey() -> Result<()> {
    test_initialise()?;
    let [a, b, c] = SUNSET;
    let energy = |edge| GS.emr_mom(EdgeIndex(edge), GS.cind(0));
    let host_scale = Atom::var(symbol!("zero_mode_test_host_scale"));
    let vertex_scale = Atom::var(symbol!("zero_mode_test_vertex_scale"));
    let host = energy(1).pow(2) + host_scale.pow(2);
    let vertex = energy(0).pow(2) + vertex_scale.pow(2);
    let (h, d) = (0.83, 1.21);
    // T sum_n (h^2 - w^2)/((w^2 + b^2)(w^2 + c^2)) d^2/a^2 at T = 1: the host
    // numerator at p0 = i w, the vertex numerator at the static q0 = 0.
    let expected = d
        * d
        * ((h * h + b * b) * static_bubble(b, c, |e| bose_tadpole(1.0, e)) - bose_tadpole(1.0, c))
        / (a * a);
    let sigma = BTreeSet::from([EdgeIndex(0)]);
    let sampling_scale = Atom::var(GS.numerator_sampling_scale);
    for vacuum_subtraction in [false, true] {
        let mut graph = sunset()?;
        let normalization = normalization(&graph)?;
        let cff = ThermalCff::sampled(&mut graph, &(&host * &vertex))?;
        let orientation = cff.orientation(vacuum_subtraction);
        let with_host = cff
            .branches(&graph, vacuum_subtraction)?
            .multiply_key_mapped(
                orientation,
                &graph,
                &host,
                DirectResidueBranches::numerator_scope(),
            )?;
        let projected = project(
            &graph,
            &with_host,
            vec![class(&graph, &[0], &[2])?],
            ZeroModeSector::Complete,
        )?
        .remove(&vec![2])
        .ok_or_else(|| eyre!("the projection returns the requested order"))?;
        let rekeyed = projected.rekey_static(orientation, &sigma)?;
        assert!(
            rekeyed.iter_keys().any(|(key, _)| key
                .source_edge_energy_map()
                .is_some_and(|map| !map[1].uniform_scale_coeff.is_zero())),
            "the host line keeps the sampling scale in its entries"
        );
        let static_term = rekeyed.multiply_key_mapped(
            orientation,
            &graph,
            &vertex,
            DirectResidueBranches::numerator_scope(),
        )?;
        for scale in [1.0, 2.5] {
            let point = Point::new(&SUNSET)
                .with(host_scale.clone(), h)
                .with(vertex_scale.clone(), d)
                .with(sampling_scale.clone(), scale);
            assert_close(
                evaluate(&graph, &static_term, &point)? / normalization,
                expected,
                1e-10,
                &format!("vacuum subtraction {vacuum_subtraction}, M = {scale}"),
            );
        }
    }
    Ok(())
}

/// The canonical coordinates of `node`, the chart of a zero-mode host.
fn canonical_lmb(graph: &Graph, node: &OwnedForestNode) -> Result<LoopMomentumBasis> {
    Ok(graph.canonical_component_lmb(node.subgraph())?)
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn difference3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scaled3(a: [f64; 3], s: f64) -> [f64; 3] {
    a.map(|component| component * s)
}

fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// The value and the first derivative at `s = 0` of a function smooth in
/// `s`, from its even and odd parts at small nonzero `s`, where a closed form
/// may be singular, extrapolated in `s^2`. The nodes balance the truncation
/// of five orders in `s^2` against the rounding of divided differences.
fn taylor_at_zero(f: impl Fn(f64) -> f64) -> (f64, f64) {
    let nodes = [0.2, 0.1, 0.05, 0.025, 0.0125];
    let squares = nodes.map(|s| s * s);
    let even = nodes.map(|s| (f(s) + f(-s)) / 2.0);
    let odd = nodes.map(|s| (f(s) - f(-s)) / (2.0 * s));
    (
        extrapolate_to_zero(&squares, &even),
        extrapolate_to_zero(&squares, &odd),
    )
}

/// The b0/b1 vertex numerator `N_B(k, k - q)` of a [`SpatialRing`]'s second
/// loop. The two-loop forms are homogeneous of degree two, one momentum per
/// vertex, so the static second loop has degree `+1`; they
/// differ in the static triangles D1/A and D2/A, whose loop runs through
/// e0, e1 and one line of the second loop.
#[derive(Clone, Copy, Debug)]
enum SecondVertices {
    /// `Q4.Q5 = k.(k - q)`: each vertex carries the momentum of one line, like
    /// the ghost vertices and like `N_A`, so each triangle sees one of them:
    /// static degree -2. The quark ring's second bubble takes this form.
    Ghost,
    /// `(Q4 + Q5)^2 = (2k - q)^2`: each vertex depends on the momenta of all
    /// its lines, like the three-gluon vertices, so each triangle sees both:
    /// the gluon ring's static degree -1.
    ThreeGluon,
    /// One: the four-gluon vertex of the tadpole ring carries no momentum.
    FourGluon,
}

impl SecondVertices {
    /// The vertices of `ring`'s second loop when a test does not choose them.
    fn of(ring: Ring) -> Self {
        match ring {
            Ring::Tadpole => Self::FourGluon,
            Ring::Quark | Ring::Mixed | Ring::Ghost => Self::Ghost,
        }
    }

    fn numerator(self) -> Atom {
        match self {
            Self::Ghost => SpatialRing::vertex(4, 5),
            Self::ThreeGluon => (1..=3)
                .map(|axis| {
                    (GS.emr_mom(EdgeIndex(4), GS.cind(axis))
                        + GS.emr_mom(EdgeIndex(5), GS.cind(axis)))
                    .pow(2)
                })
                .sum(),
            Self::FourGluon => Atom::one(),
        }
    }

    /// `N_B(k, k - q)` at the momenta `k` and `q`.
    fn value(self, k: [f64; 3], q: [f64; 3]) -> f64 {
        match self {
            Self::Ghost => SpatialRing::vertex_value(k, q),
            Self::ThreeGluon => {
                let sum = difference3(scaled3(k, 2.0), q);
                dot3(sum, sum)
            }
            Self::FourGluon => 1.0,
        }
    }
}

/// A ring topology with the spatial vertex numerators `N_A = Q2.Q3` and the
/// second loop's [`SecondVertices`]. Like the physical trace and vertices
/// they are homogeneous of degree two in the momenta of each bubble, so the
/// host degrees `dod_A = 2` and `dod_h = 4` and the static
/// degree `+1` of the second loop are physical; a constant term would be a
/// fixed scale that the soft charts expand like the soft momentum. Free of
/// energies, every static mode keeps them and the denominator-only CFF stays
/// exact.
struct SpatialRing {
    ring: Ring,
    second: SecondVertices,
    graph: Graph,
    cff: ThermalCff,
}

impl SpatialRing {
    fn new(ring: Ring) -> Result<Self> {
        Self::with_second_vertices(ring, SecondVertices::of(ring))
    }

    fn with_second_vertices(ring: Ring, second: SecondVertices) -> Result<Self> {
        Self::in_medium(ring, second, MediumMode::ThermodynamicEquilibrium)
    }

    /// The ring with the denominator-only CFF of `medium_mode`.
    fn in_medium(ring: Ring, second: SecondVertices, medium_mode: MediumMode) -> Result<Self> {
        let mut graph = ring.topology()?;
        let cff = ThermalCff::in_medium(&mut graph, medium_mode)?;
        Ok(Self {
            ring,
            second,
            graph,
            cff,
        })
    }

    /// `Q_a . Q_b` in spatial components.
    fn vertex(a: usize, b: usize) -> Atom {
        (1..=3)
            .map(|axis| {
                GS.emr_mom(EdgeIndex(a), GS.cind(axis)) * GS.emr_mom(EdgeIndex(b), GS.cind(axis))
            })
            .sum()
    }

    /// The thermal root times `N_A N_B`.
    fn root(&self, vacuum_subtraction: bool) -> Result<DirectResidueBranches> {
        self.cff
            .branches(&self.graph, vacuum_subtraction)?
            .multiply_key_mapped(
                self.cff.orientation(vacuum_subtraction),
                &self.graph,
                &(Self::vertex(2, 3) * self.second.numerator()),
                DirectResidueBranches::numerator_scope(),
            )
    }

    /// `N_A(p, p - q)` or `N_B(k, k - q)`, the vertex of the bubble with loop
    /// momentum `loop_momentum`.
    fn vertex_value(loop_momentum: [f64; 3], q: [f64; 3]) -> f64 {
        dot3(loop_momentum, difference3(loop_momentum, q))
    }

    /// The local operation's output `branches` at `points`, in the
    /// normalization of `T sum prod 1/(w^2 + E^2)`.
    fn values(&self, branches: &DirectResidueBranches, points: &[Point<'_>]) -> Result<Vec<f64>> {
        let atom =
            normalize_local_operation(&self.graph, &branches.resolved_scalars()?.factorized_sum());
        let normalization = normalization(&self.graph)?;
        Ok(evaluate_points(&self.graph, &atom, points)?
            .into_iter()
            .map(|value| value / normalization)
            .collect())
    }
}

/// The evaluation point of `kinematics` at `(beta, mu)`.
fn ring_point(kinematics: &RingKinematics, beta: f64, mu: f64) -> Point<'_> {
    kinematics.point(mu, &[]).with_beta(beta)
}

/// `|q s|^4 [x^-4] F = T Pi_A(q s) Pi_B(q s)` of the quark ring's
/// [`SpatialRing`] at the loop momenta `[q, p, k]`, at `(beta, mu)`.
fn quark_ring_static_mode([q, p, k]: [[f64; 3]; 3], s: f64, beta: f64, mu: f64) -> f64 {
    static_quark_bubble(p, scaled3(q, s), beta, mu)
        * static_quark_bubble(k, scaled3(q, s), beta, mu)
        / beta
}

/// A family of a union host's chain on the quark or tadpole ring at
/// `points`, split by the zero-mode temperature `T_zm(0)` of the ring's q
/// class: its static part, linear in `T_zm(0)`, and its vacuum-soft part,
/// free of it.
fn zero_mode_temperature_parts(
    spatial: &SpatialRing,
    family: &DirectResidueBranches,
    points: &[Point<'_>],
) -> Result<Vec<(f64, f64)>> {
    let temperature = Atom::var(symbol!("zero_mode_test_temperature"));
    let detached = family.map(|root| {
        root.replace(class_zero_mode_temperature(EdgeIndex(0)))
            .with(temperature.clone())
    });
    // Each point at the zero-mode temperature T and at zero.
    let points = points
        .iter()
        .flat_map(|point| {
            [1.0 / point.beta, 0.0].map(|value| point.clone().with(temperature.clone(), value))
        })
        .collect::<Vec<_>>();
    Ok(spatial
        .values(&detached, &points)?
        .chunks(2)
        .map(|pair| (pair[0] - pair[1], pair[1]))
        .collect())
}

/// The ring kinematics of `LOOP_MOMENTA[index]` with `q` scaled by `s`.
fn scaled_ring(index: usize, s: f64) -> RingKinematics {
    let [q, p, k] = LOOP_MOMENTA[index];
    RingKinematics::new([scaled3(q, s), p, k])
}

/// The union host of the quark ring expands both bubbles jointly to
/// first order in q. Its static part, linear in the zero-mode temperature
/// `T_zm(0)`, is `|q|^4 S_h F = T (c0 + c1 s)` at `q s`, where `c0` is
/// `T Pi_A(0) Pi_B(0)` and `c1` the slope of the exact zero mode
/// `|q|^4 [x^-4] F = T Pi_A(q) Pi_B(q)`. Vacuum subtraction adds only the
/// vacuum-soft pole `-Pi_A^vac(0) Pi_B^vac(0)/(4 |q|^3)`, free of `T_zm(0)`.
/// One operation under vacuum subtraction serves both parts; the region
/// without it only lacks the hard-vacuum order.
#[test]
fn union_host_s0_is_the_joint_taylor_product() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Quark)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let host = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial.ring.region(graph, &[2, 3, 4, 5], true)?,
    )?;
    let orientation = spatial.cff.orientation(true);
    let union_host = apply_taylor(
        &ctx,
        orientation,
        &host,
        &empty_forest_node(graph),
        None,
        &canonical_lmb(graph, &host)?,
        &spatial.root(true)?,
    )?;
    let sigma = host
        .thermal_zero_mode()
        .ok_or_else(|| eyre!("the union host has a zero-mode region"))?
        .sigma();
    assert_keyed_to_soft_set(&union_host, orientation, sigma, "S_h")?;
    let steps = [1.0, 0.5, 0.25];
    let media = [(1.0, 0.0), (1.0, QUARK_MU), (0.5, 0.0), (0.5, QUARK_MU)];
    let kinematics = (0..LOOP_MOMENTA.len())
        .flat_map(|index| steps.map(|s| scaled_ring(index, s)))
        .collect::<Vec<_>>();
    let points = kinematics
        .iter()
        .flat_map(|kinematics| media.map(|(beta, mu)| ring_point(kinematics, beta, mu)))
        .collect::<Vec<_>>();
    let parts = zero_mode_temperature_parts(&spatial, &union_host, &points)?;
    // |q s|^4 S_h F at (point, step, medium): its static part, linear in
    // T_zm(0), or its vacuum-soft part.
    let at = |index: usize, step: usize, medium: usize, thermal: bool| {
        let (static_part, vacuum_soft) = parts[(index * steps.len() + step) * media.len() + medium];
        let [q, _, _] = LOOP_MOMENTA[index];
        let value = if thermal { static_part } else { vacuum_soft };
        value * (norm3(q) * steps[step]).powi(4)
    };
    for (index, [q, p, k]) in LOOP_MOMENTA.into_iter().enumerate() {
        for (medium, (beta, mu)) in media.into_iter().enumerate() {
            let context = format!("point {index}, beta {beta}, mu {mu}");
            let (c0, c1) = taylor_at_zero(|s| quark_ring_static_mode([q, p, k], s, beta, mu));
            assert!(c1.abs() > 1e-3 * c0.abs(), "{context}: no linear term");
            let v = (0..steps.len())
                .map(|step| at(index, step, medium, true))
                .collect::<Vec<_>>();
            // Linear in s, with the coefficients of the exact zero mode.
            assert!(
                ((v[0] - v[1]) - 2.0 * (v[1] - v[2])).abs() <= 1e-9 * (c0.abs() + c1.abs()),
                "{context}: |q|^4 S_h F is not linear in s: {v:?}"
            );
            assert_close(2.0 * v[1] - v[0], c0, 1e-8, &format!("{context}: s^0"));
            // The slope is a difference of values of the size of c0, which
            // sets the precision of both sides.
            let slope = 2.0 * (v[0] - v[1]);
            assert!(
                (slope - c1).abs() <= 1e-7 * (c0.abs() + c1.abs()),
                "{context}: s^1: {slope} != {c1}"
            );
            // The vacuum-soft pole of the hard-vacuum sector.
            let vacuum = SpatialRing::vertex_value(p, [0.0; 3])
                * vacuum_bubble(norm3(p), norm3(p))
                * SpatialRing::vertex_value(k, [0.0; 3])
                * vacuum_bubble(norm3(k), norm3(k));
            for (step, s) in steps.into_iter().enumerate() {
                assert_close(
                    at(index, step, medium, false),
                    -vacuum * norm3(q) * s / 4.0,
                    1e-8,
                    &format!("{context}, s {s}: vacuum-soft pole"),
                );
            }
        }
    }
    Ok(())
}

/// The edges whose occupations `atom` names, also inside function arguments.
fn occupation_edges(atom: &Atom) -> BTreeSet<EdgeIndex> {
    let mut edges = BTreeSet::new();
    atom.as_view()
        .visitor(&mut |view| match ThermalDistributionCall::try_from(view) {
            Ok(call) => {
                edges.insert(call.edge);
                false
            }
            Err(_) => true,
        });
    edges
}

/// The zero-mode chart of the quark ring's union host writes the
/// quark occupations of both bubbles with one representative edge, so that
/// the coincident energies of the CFF variants cancel symbolically in each
/// Taylor coefficient. The chart's input, the static coefficient, still
/// names the occupation of every quark line.
#[test]
fn zero_mode_chart_unifies_occupation_representatives() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Quark)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(false);
    let host = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial.ring.region(graph, &[2, 3, 4, 5], false)?,
    )?;
    let empty = empty_forest_node(graph);
    let lmb = canonical_lmb(graph, &host)?;
    let started = start_family(
        &ctx,
        orientation,
        &host,
        &empty,
        None,
        &lmb,
        &spatial.root(false)?,
    )?
    .branches;
    let coefficient =
        ZeroModeProjection::new(vec![class(graph, &[0, 1], &[4])?], ZeroModeSector::Complete)?
            .stripped_coefficients(graph, &started)?
            .remove(&vec![4])
            .ok_or_else(|| eyre!("the projection returns the requested order"))?
            .resolved_scalars()?;
    let chart = Local3DDeformation::ZeroMode { endpoint: -3 };
    let (mut before, mut after) = (BTreeSet::new(), BTreeSet::new());
    for (_, roots) in coefficient.iter_keys() {
        for (_, root) in roots.iter() {
            before.append(&mut occupation_edges(root));
            let rescaled = Direct3dApproximation::t_rescale(
                chart,
                &ctx,
                &host,
                &empty,
                root,
                Some(host.subgraph()),
                &lmb,
                true,
            )?;
            after.append(&mut occupation_edges(&rescaled));
        }
    }
    assert_eq!(
        before,
        [2, 3, 4, 5].map(EdgeIndex).into_iter().collect(),
        "the static coefficient names every quark occupation"
    );
    assert_eq!(
        after,
        BTreeSet::from([EdgeIndex(2)]),
        "the zero-mode chart writes every quark occupation of the union host with one \
         representative"
    );
    Ok(())
}

/// The sectors of the integrated prefix `vertex` of `node` on a
/// [`SpatialRing`], from the prefix localizer that every replay enters it
/// with ([`Direct3dApproximation::localize_prefix`]). The vertex gets the
/// factor `i` of the energy integral of each contracted loop, which the
/// evaluation of a residue root divides out with those of the whole ring.
fn localized_prefix(
    spatial: &mut SpatialRing,
    vertex: &Atom,
    node: &OwnedForestNode,
    vacuum_subtraction: bool,
) -> Result<Vec<DirectSector>> {
    let settings = UVgenerationSettings::default();
    let cutset = CutSet::empty(spatial.graph.n_hedges());
    let localizer = Localizer::new(&cutset, spatial.cff.orientation(vacuum_subtraction));
    let loops = i64::try_from(node.lmb().loop_edges.len())?;
    let vertex = vertex * Atom::i().pow(loops);
    Direct3dApproximation::new(localizer, &mut spatial.graph, &settings)
        .localize_prefix(&[(&vertex, node)], node)
}

/// The spatial component `Q(0, cind(axis))` of the ring momentum `q`, and its
/// energy at `axis = 0`.
fn ring_q(axis: usize) -> Atom {
    GS.emr_mom(EdgeIndex(0), GS.cind(axis))
}

/// A symbolic integrated vertex of degree two in the vacuum mass `M` and the
/// ring momentum `q`, like the finite counterterm of a ring bubble: the
/// constant `v0 M^2`, the q-linear `v1 M q_1`, the quadratic `v2 q^2`, the
/// energy term `v3 q_0^2` and the constant `v4 mUVexp^2`, with the rational
/// `v = coefficients`. Like vakint's counterterms, whose numerators keep the
/// expansion mass `mUVexp` of the integrated U beside the terminal mass `mUV`
/// of its propagators, the vertex holds the constant under both symbols,
/// which take one value only in the final integrand.
fn ring_vertex(coefficients: [(i64, i64); 5]) -> Atom {
    let [v0, v1, v2, v3, v4] =
        coefficients.map(|(numerator, denominator)| Atom::num(numerator) / Atom::num(denominator));
    let mass = Atom::var(GS.m_uv_vacuum);
    let q_squared = (1..=3).map(|axis| ring_q(axis).pow(2)).sum::<Atom>();
    v0 * mass.pow(2)
        + v1 * &mass * ring_q(1)
        + v2 * q_squared
        + v3 * ring_q(0).pow(2)
        + v4 * Atom::var(GS.m_uv_expansion).pow(2)
}

/// The vertex coefficients of [`ring_vertex`] in the pins of contracted
/// zero-mode operations.
const RING_VERTEX: [(i64, i64); 5] = [(3, 4), (-5, 4), (3, 5), (7, 5), (-1, 2)];

/// A zero-mode host acts on a sector whose integrated prefix contracted
/// part of it: on the quark ring, the union host on the add-back of its
/// first bubble, with the second bubble's vertices that its U attaches,
/// expands the second bubble in the host's canonical chart contracted by
/// the first, `{e4}`, and refuses the uncontracted chart, for a vanishing
/// sector too; on the add-back of the join, a vertex of the host's degree
/// like the product of the bubbles' counterterms, it expands nothing in the
/// empty chart. The soft set stays outside the support.
#[test]
fn zero_mode_chart_accepts_the_contracted_canonical_chart() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let mut spatial = SpatialRing::new(Ring::Quark)?;
    let first = forest_node(&spatial.graph, &[2, 3], ApproximationType::MUV, 2)?;
    let join = forest_node(&spatial.graph, &[2, 3, 4, 5], ApproximationType::MUV, 4)?;
    let host = host_node(
        &spatial.graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(&spatial.graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let vertex = ring_vertex(RING_VERTEX);
    let [after_first] =
        <[DirectSector; 1]>::try_from(localized_prefix(&mut spatial, &vertex, &first, true)?)
            .map_err(|_| eyre!("an ordinary prefix has one sector"))?;
    let [after_join] =
        <[DirectSector; 1]>::try_from(localized_prefix(&mut spatial, &vertex.pow(2), &join, true)?)
            .map_err(|_| eyre!("a join prefix has one sector"))?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let second = graph.edge_subgraph(&[4, 5]);
    let empty = graph.empty_subgraph::<SuBitGraph>();
    let after_first = after_first.active.multiply_key_mapped(
        orientation,
        graph,
        &SpatialRing::vertex(4, 5),
        DirectResidueBranches::numerator_scope(),
    )?;
    let operation = |given: &OwnedForestNode,
                     support: &SuBitGraph,
                     lmb: &LoopMomentumBasis,
                     sector: &DirectResidueBranches| {
        apply_taylor(
            &ctx,
            orientation,
            &host,
            given,
            Some(support.clone()),
            lmb,
            sector,
        )
    };

    let (_, contracted) = coordinate_lmb(&ctx, &host, &first, None, &[], &second)?;
    assert_eq!(contracted.loop_edges.raw, [EdgeIndex(4)]);
    let static_mode = operation(&first, &second, &contracted, &after_first)?;
    assert!(
        !static_mode.is_zero(),
        "S0_h of the first bubble's add-back"
    );
    assert_keyed_to_soft_set(
        &static_mode,
        orientation,
        &edge_set(&[0, 1]),
        "S0_h on the support of the second bubble",
    )?;
    let uncontracted = canonical_lmb(graph, &host)?;
    let zero = after_first.map(|_| Atom::Zero);
    for sector in [&after_first, &zero] {
        let error = operation(&first, &second, &uncontracted, sector)
            .expect_err("the uncontracted chart is refused");
        assert!(
            format!("{error:?}").contains("must expand in its canonical chart"),
            "{error:?}"
        );
    }
    assert_eq!(
        operation(&first, &second, &contracted, &zero)?,
        zero,
        "a vanishing sector passes the check unchanged"
    );

    let (_, nothing) = coordinate_lmb(&ctx, &host, &join, None, &[], &empty)?;
    assert!(nothing.loop_edges.is_empty());
    let static_vertex = operation(&join, &empty, &nothing, &after_join.active)?;
    assert!(!static_vertex.is_zero(), "S0_h of the join's add-back");
    assert_keyed_to_soft_set(
        &static_vertex,
        orientation,
        &edge_set(&[0, 1]),
        "S0_h on the empty support",
    )?;
    let error = operation(&join, &empty, &contracted, &after_join.active)
        .expect_err("a chart with a loop is refused on the empty support");
    assert!(
        format!("{error:?}").contains("must expand in its canonical chart"),
        "{error:?}"
    );
    Ok(())
}

/// The union host on an integrated factor is the zero mode of the
/// contracted host (L1): on the quark ring, `S0_h` on the add-back of the
/// first bubble with the vertex `V(q)` of [`ring_vertex`] expands the
/// exact static mode `[x^-4](V N_B F_{G/A})` of the contracted ring in `q`
/// through its order `J = 1`: `|q s|^4 S0_h = c0 + c1 s` with the Taylor
/// coefficients of `|q s|^4 [x^-4]` at `s = 0`, while its vacuum-soft order
/// is that limit of `|q s|^3 [x^-3]` in the hard-vacuum sector. Only the
/// static mode at `q0 = 0` and the orders `j <= J` of `V` enter: the
/// quadratic and the energy terms of `V` drop out, the q-linear term stays,
/// and the constant enters with its expansion mass as with its vacuum mass.
#[test]
fn union_host_s0_on_an_integrated_factor_is_the_contracted_static_mode() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let mut spatial = SpatialRing::new(Ring::Quark)?;
    let first = forest_node(&spatial.graph, &[2, 3], ApproximationType::MUV, 2)?;
    let host = host_node(
        &spatial.graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(&spatial.graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let [v0, v1, _, _, v4] = RING_VERTEX;
    let vertices = [
        RING_VERTEX,
        [v0, v1, (0, 1), (0, 1), v4],
        [v0, (0, 1), (0, 1), (0, 1), v4],
    ];
    let mut prefixes = Vec::new();
    for coefficients in vertices {
        let [prefix] = <[DirectSector; 1]>::try_from(localized_prefix(
            &mut spatial,
            &ring_vertex(coefficients),
            &first,
            vacuum_subtraction,
        )?)
        .map_err(|_| eyre!("an ordinary prefix has one sector"))?;
        prefixes.push(prefix.active);
    }
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let second = graph.edge_subgraph(&[4, 5]);
    let (_, chart) = coordinate_lmb(&ctx, &host, &first, None, &[], &second)?;
    // The second bubble's vertices, which its operation would attach.
    let with_second = |prefix: &DirectResidueBranches| {
        prefix.multiply_key_mapped(
            orientation,
            graph,
            &SpatialRing::vertex(4, 5),
            DirectResidueBranches::numerator_scope(),
        )
    };
    let mut images = Vec::new();
    for prefix in &prefixes {
        images.push(apply_taylor(
            &ctx,
            orientation,
            &host,
            &first,
            Some(second.clone()),
            &chart,
            &with_second(prefix)?,
        )?);
    }
    let full = with_second(&prefixes[0])?;
    let mut exact = project(
        graph,
        &full,
        vec![class(graph, &[0, 1], &[4])?],
        ZeroModeSector::Complete,
    )?;
    let exact_static = exact
        .remove(&vec![4])
        .ok_or_else(|| eyre!("the static order is projected"))?;
    let mut exact = project(
        graph,
        &full,
        vec![class(graph, &[0, 1], &[3])?],
        ZeroModeSector::HardVacuum,
    )?;
    let exact_vacuum = exact
        .remove(&vec![3])
        .ok_or_else(|| eyre!("the vacuum-soft order is projected"))?;

    let steps = [1.0, 0.5, 0.25];
    let nodes = [0.2, 0.1, 0.05, 0.025, 0.0125];
    let node_steps = nodes.iter().flat_map(|s| [*s, -*s]).collect::<Vec<_>>();
    let (beta, mu) = (0.5, QUARK_MU);
    for (index, [q, _, _]) in LOOP_MOMENTA.into_iter().enumerate() {
        let context = format!("point {index}");
        let at = |family: &DirectResidueBranches, steps: &[f64]| -> Result<Vec<(f64, f64)>> {
            let kinematics = steps
                .iter()
                .map(|s| scaled_ring(index, *s))
                .collect::<Vec<_>>();
            let points = kinematics
                .iter()
                .map(|kinematics| uv_ring_point(kinematics, beta, mu, IR_MASS))
                .collect::<Vec<_>>();
            Ok(zero_mode_temperature_parts(&spatial, family, &points)?
                .into_iter()
                .zip(steps)
                .map(|((static_part, vacuum_soft), s)| {
                    let q = norm3(q) * s.abs();
                    (static_part * q.powi(4), vacuum_soft * q.powi(3))
                })
                .collect())
        };
        let exact_static = at(&exact_static, &node_steps)?;
        let exact_vacuum = at(&exact_vacuum, &node_steps)?;
        let lookup = |values: &[(f64, f64)], thermal: bool, s: f64| {
            let position = node_steps
                .iter()
                .position(|step| *step == s)
                .expect("a Taylor node");
            let (static_part, vacuum_soft) = values[position];
            if thermal { static_part } else { vacuum_soft }
        };
        let (c0, c1) = taylor_at_zero(|s| lookup(&exact_static, true, s));
        let (vacuum, _) = taylor_at_zero(|s| lookup(&exact_vacuum, false, s));
        assert!(c1.abs() > 1e-3 * c0.abs(), "{context}: no q-linear term");
        let images = images
            .iter()
            .map(|image| at(image, &steps))
            .collect::<Result<Vec<_>>>()?;
        for (step, s) in steps.into_iter().enumerate() {
            let context = format!("{context}, s {s}");
            let (static_part, vacuum_soft) = images[0][step];
            assert!(
                (static_part - (c0 + c1 * s)).abs() <= 1e-6 * (c0.abs() + c1.abs()),
                "{context}: |q|^4 S0_h = {static_part} != {c0} + {c1} s"
            );
            assert_close(
                vacuum_soft,
                vacuum,
                1e-6,
                &format!("{context}: vacuum-soft"),
            );
            // The orders beyond J and the energy term drop out exactly.
            assert_close(
                images[1][step].0,
                static_part,
                1e-10,
                &format!("{context}: V^(2)"),
            );
            assert_close(
                images[1][step].1,
                vacuum_soft,
                1e-10,
                &format!("{context}: V^(2), vacuum-soft"),
            );
            // The q-linear term of the vertex is kept.
            assert!(
                (images[2][step].0 - static_part).abs() > 1e-4 * static_part.abs(),
                "{context}: the q-linear vertex term is dropped"
            );
        }
    }
    Ok(())
}

/// A connected host's own U is part of its integrated add-back, so the
/// prefix localizer subtracts the host's zero mode from it (L2):
/// `ι_A S0_A = S0ᶜ_A ι_A`. On the gluon and ghost rings the add-back of
/// host A comes with a second sector, `-S0ᶜ_A` of it, the zero mode of the
/// contracted host alone: at its endpoint `J - dod(A) = -2` it keeps the
/// vertex at `q = 0`, `V(0) = v0 M^2 + v4 mUVexp^2`, times the static modes of the soft
/// set, exactly the projection of the add-back with that vertex. Its rows
/// name the branch of A by its order tag `lambda_A^-2`, and its step applies
/// `S0` alone.
#[test]
fn own_host_s0_on_its_integrated_prefix_is_the_static_vertex() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    for ring in [Ring::Mixed, Ring::Ghost] {
        let mut spatial = SpatialRing::new(ring)?;
        let host = host_node(
            &spatial.graph,
            &[2, 3],
            2,
            spatial
                .ring
                .region(&spatial.graph, &[2, 3], vacuum_subtraction)?,
        )?;
        let [plain, self_suffix] = <[DirectSector; 2]>::try_from(localized_prefix(
            &mut spatial,
            &ring_vertex(RING_VERTEX),
            &host,
            vacuum_subtraction,
        )?)
        .map_err(|_| eyre!("{ring:?}: the add-back of host A has a self-suffix"))?;
        let constant = localized_prefix(
            &mut spatial,
            &ring_vertex([RING_VERTEX[0], (0, 1), (0, 1), (0, 1), RING_VERTEX[4]]),
            &host,
            vacuum_subtraction,
        )?
        .into_iter()
        .next()
        .ok_or_else(|| eyre!("{ring:?}: the add-back has a plain sector"))?;
        let graph = &spatial.graph;
        assert!(plain.projection_paths.is_empty());
        let [path] = self_suffix.projection_paths.as_slice() else {
            panic!("{ring:?}: the self-suffix has one projection path");
        };
        let [step] = path.steps.as_slice() else {
            panic!("{ring:?}: the self-suffix is one step");
        };
        assert_eq!(
            step.conceptual_branches,
            vec![Local3DSignedBranch {
                branch: Local3DConceptualBranch::S,
                coefficient: 1,
            }],
            "{ring:?}: the self-suffix applies S0 alone"
        );
        assert!(
            step.route_loop_edges.is_empty(),
            "{ring:?}: the contracted host has no chart"
        );
        let order_tag = zero_mode_order(Atom::var(host.subgraph().symbol())).pow(-2);
        assert!(
            self_suffix
                .active
                .factorized_sum()
                .contains(order_tag.as_view()),
            "{ring:?}: the self-suffix rows carry lambda_A^-2"
        );
        assert_eq!(self_suffix.frozen_integrands, plain.frozen_integrands);

        let mut projected = project(
            graph,
            &constant.active,
            vec![class(graph, &[0, 1], &[4])?, class(graph, &[4], &[2])?],
            ZeroModeSector::Complete,
        )?;
        let expected = projected
            .remove(&vec![4, 2])
            .ok_or_else(|| eyre!("the static order is projected"))?;
        let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
        let points = kinematics
            .iter()
            .flat_map(|kinematics| {
                PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, IR_MASS))
            })
            .collect::<Vec<_>>();
        let actual = spatial.values(&-self_suffix.active, &points)?;
        let expected = spatial.values(&expected, &points)?;
        for (index, (actual, expected)) in actual.into_iter().zip(expected).enumerate() {
            assert!(
                expected.abs() > 0.0,
                "{ring:?}, point {index}: no static vertex"
            );
            assert_close(actual, expected, 1e-10, &format!("{ring:?}, point {index}"));
        }
    }
    Ok(())
}

/// A connected host annihilates an integrated prefix that contracted one
/// of its frequency classes: on the gluon ring, `S0_A` on the add-back of
/// the second loop B, which holds the class `{e4}` of the soft set of A, is
/// exactly zero, so the host's operation there is its plain U alone,
/// `H_A ι_B = U_A ι_B`, the replay the union states use
/// (`Forests::replayed_plain_u`).
#[test]
fn zero_mode_host_annihilates_an_integrated_class() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let mut spatial = SpatialRing::new(Ring::Mixed)?;
    let second = forest_node(&spatial.graph, &[4, 5], ApproximationType::MUV, 2)?;
    let host = host_node(
        &spatial.graph,
        &[2, 3],
        2,
        spatial
            .ring
            .region(&spatial.graph, &[2, 3], vacuum_subtraction)?,
    )?;
    let plain_u = forest_node(&spatial.graph, &[2, 3], ApproximationType::MUV, 2)?;
    let [prefix] = <[DirectSector; 1]>::try_from(localized_prefix(
        &mut spatial,
        &ring_vertex(RING_VERTEX),
        &second,
        vacuum_subtraction,
    )?)
    .map_err(|_| eyre!("an ordinary prefix has one sector"))?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let root = empty_forest_node(graph);
    let lmb = canonical_lmb(graph, &host)?;
    let region = host
        .thermal_zero_mode()
        .ok_or_else(|| eyre!("host A has a zero-mode region"))?;
    let started = start_family(&ctx, orientation, &host, &root, None, &lmb, &prefix.active)?;
    let terms = ZeroModeTerm::expand(
        &ctx,
        orientation,
        &host,
        &root,
        region,
        host.subgraph(),
        &lmb,
        &started,
    )?;
    assert!(
        terms.iter().all(|term| term.frozen.is_zero()),
        "S0_A of the add-back of B vanishes"
    );
    let operation = |node: &OwnedForestNode| {
        apply_taylor(&ctx, orientation, node, &root, None, &lmb, &prefix.active)
    };
    let hosted = operation(&host)?;
    let plain = operation(&plain_u)?;
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| {
            PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, IR_MASS))
        })
        .collect::<Vec<_>>();
    let hosted = spatial.values(&hosted, &points)?;
    let plain = spatial.values(&plain, &points)?;
    for (index, (hosted, plain)) in hosted.into_iter().zip(plain).enumerate() {
        assert!(plain.abs() > 0.0, "point {index}: U_A of the add-back of B");
        assert_close(hosted, plain, 1e-12, &format!("point {index}"));
    }
    Ok(())
}

/// A zero-mode chart expands an integrated vertex in the soft momenta at
/// fixed masses: it scales the vertex's algebraic vacuum masses with its
/// loops, but holds its logarithms of the vacuum mass over the MS-bar scale
/// fixed. On the quark ring, `S0_h` on the add-back of the first bubble with
/// a vertex whose constant and q-linear terms, the latter at a kept order,
/// carry `log(M^2/mu^2)` leaves no logarithm depending on the order tag, and
/// equals `S0_h` on the vertex with a fixed symbol `L` in place of the
/// logarithm, at `L = log(M^2/mu^2)`.
#[test]
fn zero_mode_chart_holds_prefix_logarithms() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let mut spatial = SpatialRing::new(Ring::Quark)?;
    let first = forest_node(&spatial.graph, &[2, 3], ApproximationType::MUV, 2)?;
    let host = host_node(
        &spatial.graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(&spatial.graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let mass = Atom::var(GS.m_uv_vacuum);
    let logarithm_symbol = Atom::var(symbol!("zero_mode_test_logarithm"));
    let vertex = |logarithm: &Atom| {
        mass.pow(2) * (Atom::num(2) + Atom::num(3) * logarithm)
            + &mass * ring_q(1) * (Atom::num(-1) + Atom::num(5) * logarithm)
            + ring_vertex([(0, 1), (0, 1), (3, 5), (7, 5), (0, 1)])
    };
    let logarithm = function!(Symbol::LOG, mass.pow(2) / Atom::var(GS.mu_r_sq));
    let mut prefixes = Vec::new();
    for logarithm in [&logarithm, &logarithm_symbol] {
        let [prefix] = <[DirectSector; 1]>::try_from(localized_prefix(
            &mut spatial,
            &vertex(logarithm),
            &first,
            vacuum_subtraction,
        )?)
        .map_err(|_| eyre!("an ordinary prefix has one sector"))?;
        prefixes.push(prefix.active);
    }
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let second = graph.edge_subgraph(&[4, 5]);
    let (_, chart) = coordinate_lmb(&ctx, &host, &first, None, &[], &second)?;
    let images = prefixes
        .iter()
        .map(|prefix| {
            apply_taylor(
                &ctx,
                orientation,
                &host,
                &first,
                Some(second.clone()),
                &chart,
                &prefix.multiply_key_mapped(
                    orientation,
                    graph,
                    &SpatialRing::vertex(4, 5),
                    DirectResidueBranches::numerator_scope(),
                )?,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let with_logarithm = images[0].resolved_scalars()?.factorized_sum();
    let mut logarithms = 0;
    with_logarithm.as_view().visitor(&mut |view| {
        if let AtomView::Fun(call) = view
            && call.get_symbol() == Symbol::LOG
        {
            logarithms += 1;
            assert!(
                !view.contains_symbol(*ZERO_MODE_ORDER) && !view.contains_symbol(GS.rescale),
                "a logarithm depends on the order tag: {view}"
            );
            return false;
        }
        true
    });
    assert!(logarithms > 0, "the logarithm of the vertex survives");
    let mu_r = 2.0_f64;
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| {
            PIN_MEDIA.map(|(beta, mu)| {
                uv_ring_point(kinematics, beta, mu, IR_MASS)
                    .with(Atom::var(GS.mu_r_sq), mu_r * mu_r)
                    .with(
                        logarithm_symbol.clone(),
                        (UV_MASS * UV_MASS / (mu_r * mu_r)).ln(),
                    )
            })
        })
        .collect::<Vec<_>>();
    let actual = spatial.values(&images[0], &points)?;
    let fixed = spatial.values(&images[1], &points)?;
    for (index, (actual, fixed)) in actual.into_iter().zip(fixed).enumerate() {
        assert!(fixed.abs() > 0.0, "point {index}: S0_h of the vertex");
        assert_close(actual, fixed, 1e-10, &format!("point {index}"));
    }
    Ok(())
}

/// Under vacuum subtraction the thermal weight `W - V(W)` vanishes on a
/// cograph without loops, so the add-back of an operation that covers every
/// loop localizes to zero, as its local U does; without vacuum subtraction
/// it is the vacuum weight times the kernel. The localizer needs no rule of
/// its own for it: the zero-temperature lane already gives such an operation
/// no add-back (`Forests::has_zero_temperature_counterpart`), and the
/// localized weight vanishes by itself.
#[test]
fn vacuum_subtracted_loopless_cograph_localizes_to_zero() -> Result<()> {
    test_initialise()?;
    let mut spatial = SpatialRing::new(Ring::Quark)?;
    let whole = whole_ring_node(&spatial.graph)?;
    let vertex = Atom::var(GS.m_uv_vacuum).pow(4);
    for vacuum_subtraction in [true, false] {
        let [localized] = <[DirectSector; 1]>::try_from(localized_prefix(
            &mut spatial,
            &vertex,
            &whole,
            vacuum_subtraction,
        )?)
        .map_err(|_| eyre!("an ordinary prefix has one sector"))?;
        assert_eq!(
            localized.active.is_zero(),
            vacuum_subtraction,
            "vacuum subtraction {vacuum_subtraction}"
        );
    }
    Ok(())
}

/// Every nonzero sector of a zero-mode host expands in the
/// canonical chart its region records, which the operations that cancel
/// against it share.
#[test]
fn host_operation_rejects_a_non_canonical_chart() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Quark)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let host = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial.ring.region(graph, &[2, 3, 4, 5], false)?,
    )?;
    let canonical = canonical_lmb(graph, &host)?;
    assert_eq!(canonical.loop_edges.raw, [2, 4].map(EdgeIndex));
    let other = graph
        .generate_loop_momentum_bases_of(host.subgraph())
        .into_iter()
        .find(|lmb| {
            lmb.loop_edges.iter().copied().collect::<BTreeSet<_>>()
                != canonical.loop_edges.iter().copied().collect()
        })
        .ok_or_else(|| eyre!("the union host has another loop basis"))?;
    let error = apply_taylor(
        &ctx,
        spatial.cff.orientation(false),
        &host,
        &empty_forest_node(graph),
        None,
        &other,
        &spatial.root(false)?,
    )
    .expect_err("a non-canonical chart is refused");
    assert!(
        format!("{error:?}").contains("must expand in its canonical chart"),
        "{error:?}"
    );
    Ok(())
}

/// The connected host A of the mixed ring with a [`SpatialRing`]'s
/// numerators, its canonical chart and the terms of `S0_A F`, without the U
/// terms of its operation.
struct MixedRingHostA {
    node: OwnedForestNode,
    lmb: LoopMomentumBasis,
    terms: Vec<ZeroModeTerm>,
}

impl MixedRingHostA {
    fn new(spatial: &SpatialRing, vacuum_subtraction: bool) -> Result<Self> {
        let graph = &spatial.graph;
        let settings = UVgenerationSettings::default();
        let ctx = UVCtx::new(graph, &settings);
        let orientation = spatial.cff.orientation(vacuum_subtraction);
        let node = host_node(
            graph,
            &[2, 3],
            2,
            spatial.ring.region(graph, &[2, 3], vacuum_subtraction)?,
        )?;
        let root = empty_forest_node(graph);
        let lmb = canonical_lmb(graph, &node)?;
        let started = start_family(
            &ctx,
            orientation,
            &node,
            &root,
            None,
            &lmb,
            &spatial.root(vacuum_subtraction)?,
        )?;
        let region = node
            .thermal_zero_mode()
            .ok_or_else(|| eyre!("host A has a zero-mode region"))?;
        let terms = ZeroModeTerm::expand(
            &ctx,
            orientation,
            &node,
            &root,
            region,
            node.subgraph(),
            &lmb,
            &started,
        )?;
        let sigma = region.sigma();
        for term in &terms {
            assert_keyed_to_soft_set(&term.frozen, orientation, sigma, "S0_A")?;
        }
        Ok(Self { node, lmb, terms })
    }

    /// `S0_A F`.
    fn s0(&self) -> Result<DirectResidueBranches> {
        ZeroModeTerm::restore(&self.terms, |term| Ok(Cow::Borrowed(&term.frozen)))
    }

    /// `U_A S0_A F`, the U of A's own support on its zero-mode branch, as
    /// the operation of host A subtracts it.
    fn own_u_on_s0(&self, ctx: &UVCtx<'_>) -> Result<DirectResidueBranches> {
        let root = empty_forest_node(ctx.graph);
        ZeroModeTerm::restore(&self.terms, |term| {
            Ok(Cow::Owned(Direct3dApproximation::t(
                Local3DDeformation::OrdinaryOnSoftBranch,
                ctx,
                &self.node,
                &root,
                &term.frozen,
                None,
                &self.lmb,
            )?))
        })
    }
}

/// Whether every summand of every root of `family` has the factor
/// `OSE(edge)^-power` for each `(edge, power)` at its top level, where an
/// enclosing projection reads it off.
fn has_top_level_propagators(family: &DirectResidueBranches, propagators: &[(usize, i64)]) -> bool {
    let propagators = propagators
        .iter()
        .map(|(edge, power)| ose_atom_from_index(EdgeIndex(*edge)).pow(Atom::num(-power)))
        .collect::<Vec<_>>();
    family.iter_keys().all(|(_, roots)| {
        roots.iter().all(|(_, root)| {
            let summands = match root.as_view() {
                AtomView::Add(sum) => sum.iter().collect(),
                view => vec![view],
            };
            summands.into_iter().all(|summand| {
                let factors = match summand {
                    AtomView::Mul(product) => product.iter().collect(),
                    view => vec![view],
                };
                propagators
                    .iter()
                    .all(|propagator| factors.contains(&propagator.as_view()))
            })
        })
    })
}

/// Assert that every row of the zero-mode branch `family` carries its
/// static source map and is keyed to the static energies of the soft set
/// `sigma`, so that re-keying it again changes nothing.
fn assert_keyed_to_soft_set(
    family: &DirectResidueBranches,
    orientation: OrientationProjection<'_>,
    sigma: &BTreeSet<EdgeIndex>,
    context: &str,
) -> Result<()> {
    assert!(
        family
            .iter_keys()
            .all(|(key, _)| key.source_edge_energy_map().is_some()),
        "{context}: every row carries its static source map"
    );
    assert_eq!(
        family.rekey_static(orientation, sigma)?,
        *family,
        "{context}: the rows are keyed to the static energies of the soft set"
    );
    Ok(())
}

/// The rows of `family` that a zero-mode operation re-keyed.
fn rekeyed_rows(family: &DirectResidueBranches) -> Result<DirectResidueBranches> {
    DirectResidueBranches::from_keyed(
        family
            .iter_keys()
            .filter(|(key, _)| key.source_edge_energy_map().is_some())
            .map(|(key, roots)| (key.clone(), roots.clone())),
    )
}

/// `Pi(q; l) = N(l, l - q) T sum 1/((w^2 + l^2)(w^2 + (l - q)^2))`, the
/// static quark bubble of a [`SpatialRing`] loop with momentum `l`, at
/// `(beta, mu)`.
fn static_quark_bubble(l: [f64; 3], q: [f64; 3], beta: f64, mu: f64) -> f64 {
    SpatialRing::vertex_value(l, q)
        * static_bubble(norm3(l), norm3(difference3(l, q)), |e| {
            fermi_tadpole(beta, e, mu)
        })
}

/// `Pi_A^(0)(p) = N_A(p, p) T sum 1/(w^2 + p^2)^2`, the static quark bubble
/// of [`SpatialRing`] at `q = 0`, at `(beta, mu)`.
fn static_bubble_at_zero(p: [f64; 3], beta: f64, mu: f64) -> f64 {
    let [q, _, _] = LOOP_MOMENTA[0];
    taylor_at_zero(|s| static_quark_bubble(p, scaled3(q, s), beta, mu)).0
}

/// `S0` of the connected host A keeps the static second loop exact,
/// `S_A F = T^2 Pi_A^(0)(p) N_B(k, k - q)/(q^4 k^2 (k - q)^2)`, with both
/// static class propagators as top-level factors, with and without vacuum
/// subtraction.
#[test]
fn inner_s0_keeps_the_static_bubble_exact() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let media = [(1.0, 0.0), (1.0, QUARK_MU), (0.5, 0.0), (0.5, QUARK_MU)];
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| media.map(|(beta, mu)| ring_point(kinematics, beta, mu)))
        .collect::<Vec<_>>();
    for vacuum_subtraction in [false, true] {
        let s_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?.s0()?;
        assert!(
            has_top_level_propagators(&s_a, &[(0, 4), (4, 2)]),
            "S_A keeps its static class propagators as top-level factors"
        );
        let mut values = spatial.values(&s_a, &points)?.into_iter();
        for (index, [q, p, k]) in LOOP_MOMENTA.into_iter().enumerate() {
            for (beta, mu) in media {
                let context = format!(
                    "vacuum subtraction {vacuum_subtraction}, point {index}, beta {beta}, mu {mu}"
                );
                let static_loop = SpatialRing::vertex_value(k, q)
                    / (norm3(k).powi(2) * norm3(difference3(k, q)).powi(2));
                let expected = static_bubble_at_zero(p, beta, mu) * static_loop
                    / (beta * beta * norm3(q).powi(4));
                let actual = values.next().expect("one value per point");
                assert_close(actual, expected, 1e-9, &context);
            }
        }
    }
    Ok(())
}

/// The union host of the mixed ring on `S_A F` co-scales the
/// static temperature of the second loop inside it, and so keeps that
/// loop's term linear in q: `S S_A F = T^2 Pi_A^(0) (B_0^(0) + B_0^(1))/q^4`
/// with `B_0 = N_B(k, k - q)/(k^2 (k - q)^2)` per unit temperature. Its
/// vacuum-soft order reads the projected rows off as zero.
#[test]
fn nested_s0_keeps_the_q_linear_static_term() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| PIN_MEDIA.map(|(beta, mu)| ring_point(kinematics, beta, mu)))
        .collect::<Vec<_>>();
    for vacuum_subtraction in [false, true] {
        let host_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?;
        let union = host_node(
            graph,
            &[2, 3, 4, 5],
            4,
            spatial
                .ring
                .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
        )?;
        let nested = apply_taylor(
            &ctx,
            spatial.cff.orientation(vacuum_subtraction),
            &union,
            &host_a.node,
            None,
            &canonical_lmb(graph, &union)?,
            &host_a.s0()?,
        )?;
        let mut values = spatial.values(&nested, &points)?.into_iter();
        for (index, [q, p, k]) in LOOP_MOMENTA.into_iter().enumerate() {
            let (b0, b1) = taylor_at_zero(|s| {
                let shifted = difference3(k, scaled3(q, s));
                SpatialRing::vertex_value(k, scaled3(q, s))
                    / (norm3(k).powi(2) * norm3(shifted).powi(2))
            });
            assert!(b1.abs() > 1e-3 * b0.abs(), "point {index}: no linear term");
            for (beta, mu) in PIN_MEDIA {
                let context = format!(
                    "vacuum subtraction {vacuum_subtraction}, point {index}, beta {beta}, mu {mu}"
                );
                let expected = static_bubble_at_zero(p, beta, mu) * (b0 + b1)
                    / (beta * beta * norm3(q).powi(4));
                let actual = values.next().expect("one value per point");
                assert_close(actual, expected, 1e-9, &context);
            }
        }
    }
    Ok(())
}

/// The rows a connected host completes keep their class
/// propagators as top-level factors through its chart, its freeze, its own U
/// and a sibling U. The hard-vacuum sector of the enclosing union host
/// therefore reads them off at the order `2 nu`, which makes its
/// vacuum-soft coefficient zero on them, instead of expanding a frozen branch
/// again (`ZeroModeError::FrozenHardVacuum`).
#[test]
fn frozen_zero_mode_rows_skip_the_hard_vacuum_sector() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(true);
    let root = empty_forest_node(graph);
    let host_a = host_node(
        graph,
        &[2, 3],
        2,
        spatial.ring.region(graph, &[2, 3], true)?,
    )?;
    let h_a = apply_taylor(
        &ctx,
        orientation,
        &host_a,
        &root,
        None,
        &canonical_lmb(graph, &host_a)?,
        &spatial.root(true)?,
    )?;
    let second = forest_node(graph, &[4, 5], ApproximationType::MUV, 2)?;
    let u_b = apply_taylor(
        &ctx,
        orientation,
        &second,
        &root,
        None,
        &canonical_lmb(graph, &second)?,
        &h_a,
    )?;
    let twin = forest_node(graph, &[2, 3, 4, 5], ApproximationType::MUV, 4)?;
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial.ring.region(graph, &[2, 3, 4, 5], true)?,
    )?;
    let union_lmb = canonical_lmb(graph, &union)?;
    let vacuum_soft = ZeroModeProjection::new(
        vec![class(graph, &[0, 1], &[3])?],
        ZeroModeSector::HardVacuum,
    )?;
    for (given, family, context) in [(&host_a, &h_a, "H_A F"), (&twin, &u_b, "U_B H_A F")] {
        let completed = rekeyed_rows(family)?;
        assert!(!completed.is_zero(), "{context}: completed rows");
        assert!(
            has_top_level_propagators(&completed, &[(0, 4)]),
            "{context}: the q-class propagator must stay a top-level factor"
        );
        let started =
            start_family(&ctx, orientation, &union, given, None, &union_lmb, family)?.branches;
        let odd = vacuum_soft.coefficients(graph, &rekeyed_rows(&started)?)?;
        assert!(
            odd[&vec![3]].is_zero(),
            "{context}: projected rows have no vacuum-soft pole"
        );
        vacuum_soft
            .coefficients(graph, &started)
            .map_err(|error| eyre!("{context}: {error}"))?;
    }
    Ok(())
}

/// The final integrand of the residue `branches` of the empty forest node,
/// with every residue selector at one.
fn final_root_integrand(
    graph: &mut Graph,
    orientation: OrientationProjection<'_>,
    branches: &DirectResidueBranches,
) -> Result<Atom> {
    let cutset = CutSet::empty(graph.n_hedges());
    let settings = local_3d_settings();
    let empty = empty_forest_node(graph);
    let integrands = FinalIntegrandBuilder::new(Localizer::new(&cutset, orientation), &settings)
        .build_direct(
            graph,
            &empty,
            &Direct3dCts::Root(branches.clone()),
            &IntegratedCts::root(),
        )?
        .into_integrands()
        .resolved()?;
    Ok(integrands
        .iter()
        .map(|(_, atom)| atom)
        .sum::<Atom>()
        .replace(function!(OrientationID::symbol(), W_.a_))
        .with(Atom::one()))
}

/// The final integrand normalizes the zero-mode tags, `lambda -> 1` and
/// `T_zm -> 1/beta`, checked at a temperature other than one. The U of a
/// zero-mode branch's own support refuses a zero-mode temperature that
/// branch took over, while it keeps the branch's order tag, and an enclosing
/// U keeps and co-scales both.
#[test]
fn zero_mode_tags_normalize_and_reject_own_role() -> Result<()> {
    test_initialise()?;
    let mut graph = sunset()?;
    let cff = ThermalCff::new(&mut graph)?;
    let orientation = cff.orientation(false);
    let root = cff.branches(&graph, false)?;
    let branch = graph.edge_subgraph(&[1, 2]);
    let owner = Atom::var(branch.symbol());
    let order_tag = zero_mode_order(owner.clone());
    let branch_temperature = zero_mode_temperature(owner);
    let tags =
        order_tag.pow(-2) * class_zero_mode_temperature(EdgeIndex(0)) * branch_temperature.pow(2);
    let plain = final_root_integrand(&mut graph, orientation, &root)?;
    let tagged = final_root_integrand(&mut graph, orientation, &root.map(|atom| atom * &tags))?;
    assert!(
        !tagged.contains_symbol(*ZERO_MODE_ORDER)
            && !tagged.contains_symbol(*ZERO_MODE_TEMPERATURE),
        "the final integrand keeps a zero-mode tag"
    );
    // The sunset's final integrand is imaginary at the point: compare complex
    // values, whose scale is nonzero.
    for beta in [0.5, 1.0] {
        let point = Point::new(&SUNSET).with_beta(beta);
        let value = |atom: &Atom| -> Result<Complex<f64>> {
            Ok(evaluate_points_complex(&graph, atom, slice::from_ref(&point))?[0])
        };
        let untagged = value(&plain)?;
        let scale = beta.powi(3);
        assert_close_complex(
            value(&tagged)?,
            Complex::new(untagged.re / scale, untagged.im / scale),
            1e-12,
            &format!("beta {beta}"),
        );
    }

    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(&graph, &settings);
    let empty = empty_forest_node(&graph);
    let child = forest_node(&graph, &[1, 2], ApproximationType::IR, 1)?;
    let parent = forest_node(&graph, &[0, 1, 2], ApproximationType::MUV, 2)?;
    let energy = ose_atom_from_index(EdgeIndex(1)).pow(-3);
    let own_u = |source: &Atom| {
        Local3DLoopRescaling::FullSubgraph.project(
            Local3DDeformation::OrdinaryOnSoftBranch,
            &ctx,
            &child,
            &empty,
            source,
            None,
            child.lmb(),
        )
    };
    let error = own_u(&(&branch_temperature * &energy))
        .expect_err("the own U must refuse a zero-mode temperature its branch owns");
    assert!(
        format!("{error:?}").contains("reached the U of its own branch"),
        "{error:?}"
    );
    let kept = own_u(&(&order_tag * &energy))?;
    assert!(
        kept.contains(order_tag.as_view()),
        "the own U keeps the order tag"
    );
    let enclosing = Local3DLoopRescaling::FullSubgraph.t(
        &ctx,
        &parent,
        &child,
        &(&branch_temperature * &order_tag * &energy),
        Some(parent.subgraph()),
        parent.lmb(),
    )?;
    assert!(
        enclosing.contains(branch_temperature.as_view()) && enclosing.contains(order_tag.as_view()),
        "an enclosing U keeps the frozen tags: {enclosing}"
    );
    // Its chart scales both with its loop momenta, as the soft deformation
    // scaled the frozen vertex: by `1/t` in the chart's inverted Taylor
    // variable `t`.
    for tag in [&branch_temperature, &order_tag] {
        let rescaled = Direct3dApproximation::t_rescale(
            Local3DDeformation::Ordinary,
            &ctx,
            &parent,
            &child,
            tag,
            Some(parent.subgraph()),
            parent.lmb(),
            false,
        )?;
        assert_eq!(
            rescaled,
            tag / Atom::var(GS.rescale),
            "an enclosing U co-scales the frozen tag {tag}"
        );
    }
    Ok(())
}

/// The zero-mode chart probes its residue roots through
/// its endpoint only. On host A of the mixed ring both probes give the same
/// series.
#[test]
fn zero_mode_chart_probe_to_depth_matches_the_probe_through_zero() -> Result<()> {
    test_initialise()?;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(true);
    let host = host_node(
        graph,
        &[2, 3],
        2,
        spatial.ring.region(graph, &[2, 3], true)?,
    )?;
    let empty = empty_forest_node(graph);
    let lmb = canonical_lmb(graph, &host)?;
    let started = start_family(
        &ctx,
        orientation,
        &host,
        &empty,
        None,
        &lmb,
        &spatial.root(true)?,
    )?
    .branches;
    let coefficient = ZeroModeProjection::new(
        vec![class(graph, &[0, 1], &[4])?, class(graph, &[4], &[2])?],
        ZeroModeSector::Complete,
    )?
    .stripped_coefficients(graph, &started)?
    .remove(&vec![4, 2])
    .ok_or_else(|| eyre!("the projection returns the requested order"))?
    .resolved_scalars()?;
    let chart = Local3DDeformation::ZeroMode { endpoint: -2 };
    let rescale = |atom: &Atom, include_measure| {
        Direct3dApproximation::t_rescale(
            chart,
            &ctx,
            &host,
            &empty,
            atom,
            Some(host.subgraph()),
            &lmb,
            include_measure,
        )
    };
    let rescaled = coefficient
        .fallible_map(|_, atom| rescale(atom, true))?
        .map_numerators(|entry| rescale(&entry.rhs, false))?;
    let scope = || DirectResidueBranches::numerator_scope().1;
    let center = Atom::Zero;
    let through_zero =
        rescaled.series_preserving_numerators(GS.rescale, center.as_view(), -2, scope())?;
    let to_depth = rescaled.series_preserving_numerators_to_depth(
        GS.rescale,
        center.as_view(),
        -2,
        scope(),
    )?;
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .map(|kinematics| ring_point(kinematics, 0.5, QUARK_MU))
        .collect::<Vec<_>>();
    let [through_zero, to_depth] = [through_zero, to_depth].map(|series| {
        spatial.values(
            &series.map(|root| root.replace(GS.rescale).with(Atom::one())),
            &points,
        )
    });
    for (index, (through_zero, to_depth)) in through_zero?.into_iter().zip(to_depth?).enumerate() {
        assert!(
            through_zero.abs() > 0.0,
            "point {index}: the chart vanishes"
        );
        assert_close(to_depth, through_zero, 1e-12, &format!("point {index}"));
    }
    Ok(())
}

/// The mass `M` of the U operations in the depth pins of an enclosing U on a
/// frozen zero-mode branch.
const UV_MASS: f64 = 1.3;

/// The static-theory mass `mIR` of a U on a static cograph, where the pins
/// tell it apart from [`UV_MASS`].
const IR_MASS: f64 = 0.45;

/// The static-theory masses of the pins: equal to `M`, where the two masses
/// cannot be told apart, and apart from it.
const STATIC_THEORY_MASSES: [f64; 2] = [UV_MASS, IR_MASS];

/// The evaluation point of `kinematics` at `(beta, mu)`, the UV mass and the
/// static-theory mass `ir_mass`.
fn uv_ring_point(kinematics: &RingKinematics, beta: f64, mu: f64, ir_mass: f64) -> Point<'_> {
    ring_point(kinematics, beta, mu)
        .with(Atom::var(GS.m_uv_vacuum), UV_MASS)
        .with(Atom::var(GS.m_ir_vacuum), ir_mass)
}

/// `1/(l^2 + m^2)`: a static line of momentum `l` after an enclosing U of
/// mass `m`, at its leading order.
fn uv_propagator(l: [f64; 3], mass: f64) -> f64 {
    1.0 / (dot3(l, l) + mass * mass)
}

/// `Pi^vac(0; l) = N(l, l) int dw/(2 pi) 1/(w^2 + l^2)^2 = l^2/(4 |l|^3)`,
/// the vacuum part of a [`SpatialRing`] quark bubble with loop momentum `l`
/// at zero external momentum.
fn vacuum_static_bubble(l: [f64; 3]) -> f64 {
    SpatialRing::vertex_value(l, [0.0; 3]) * vacuum_bubble(norm3(l), norm3(l))
}

/// The 4D U of a [`SpatialRing`] quark bubble with loop momentum `l` at the
/// static external momentum `q`, through first order in `q`: its value at
/// `q = 0` and its term linear in `q`.
///
/// U takes the vacuum limit and expands
/// `f(t) = N(l, l - t q) vb(E_l(t), E_{l-q}(t))` with
/// `vb(a, b) = 1/(2 a b (a + b))` and the rearranged energies
/// `E_l(t)^2 = l^2 + M^2 - t^2 M^2`,
/// `E_{l-q}(t)^2 = (l - t q)^2 + M^2 - t^2 M^2` through `t^2`, its degree
/// two, at `t = 1`. With `E^2 = l^2 + M^2` and `d ln vb/db = -3/(2 E)` at
/// `a = b = E`:
/// - `q^0`: `f(0) + [t^2] f = l^2/(4 E^3) (1 + 3 M^2/(2 E^2))`;
/// - `q^1`: `f'(0) = (l.q)/(4 E^3) (3 l^2/(2 E^2) - 1)`;
///
/// the rest of `[t^2] f` is of order `q^2`.
fn uv_static_bubble(l: [f64; 3], q: [f64; 3]) -> (f64, f64) {
    let energy_squared = dot3(l, l) + UV_MASS * UV_MASS;
    let leading = 1.0 / (4.0 * energy_squared.powf(1.5));
    (
        dot3(l, l) * leading * (1.0 + 1.5 * UV_MASS * UV_MASS / energy_squared),
        dot3(l, q) * leading * (1.5 * dot3(l, l) / energy_squared - 1.0),
    )
}

/// The static second loop `N_B(k, k - q)/((k^2 + m^2)((k - q)^2 + m^2))` of
/// the mixed ring, per unit `T_zm(4)`, through first order in `q` at fixed
/// `m`: `(B^(0), B^(1))`. With `N_B = N_0 + N_1 + N_2` by powers of `q`,
/// `N_1` read off exactly from the quadratic `N_B` at `q` and `-q`, and
/// `(k - q)^2 + m^2 = (k^2 + m^2) - 2 k.q + q^2`:
/// `B^(0) = N_0/(k^2 + m^2)^2`,
/// `B^(1) = N_1/(k^2 + m^2)^2 + 2 (k.q) N_0/(k^2 + m^2)^3`.
/// At `m = 0` these are the Taylor terms `B_0^(0)`, `B_0^(1)` of the static
/// mode; at the static-theory mass `m = m_IR` they are `[U B_0]`, U_B of the
/// static loop at its degree `+1`, whose mass rearrangement would only start
/// at the relative order two.
fn static_second_loop(second: SecondVertices, k: [f64; 3], q: [f64; 3], mass: f64) -> (f64, f64) {
    let constant = second.value(k, [0.0; 3]);
    let linear = (second.value(k, q) - second.value(k, scaled3(q, -1.0))) / 2.0;
    let denominator = dot3(k, k) + mass * mass;
    (
        constant / denominator.powi(2),
        linear / denominator.powi(2) + 2.0 * dot3(k, q) * constant / denominator.powi(3),
    )
}

/// The U node of the whole ring. Every [`SpatialRing`] has the 4D degree
/// four: `12 - 2 * 6 + 2 + 2` with two bubbles, `12 - 2 * 5 + 2` with the
/// tadpole.
fn whole_ring_node(graph: &Graph) -> Result<OwnedForestNode> {
    let edges = (0..graph.underlying.n_edges()).collect::<Vec<_>>();
    forest_node(graph, &edges, ApproximationType::MUV, 4)
}

/// The factors of the mixed ring's closed forms at one point and medium, for
/// a [`SpatialRing`]'s numerators. Derived by hand from the static mode
/// `S_A F = T^2 Pi_A^(0)(p) N_B/(q^4 k^2 (k - q)^2)` of host A. The own U_A
/// of the host introduces the UV mass `M`; U_B and U_G, whose supports hold
/// static lines of `S_A`, introduce the static-theory mass `m_IR` on the lines
/// no inner U made massive.
struct MixedRingFactors {
    /// `T_zm(0) T_zm(4) = T^2`, the temperatures of the static q and k modes.
    temperatures: f64,
    /// The frozen vertex `Pi_A^(0)(p)`, and its U `U_A Pi_A^(0)(p)`.
    bubble: f64,
    uv_bubble: f64,
    /// The static ring `1/q^4`, and `1/(q^2 + m_IR^2)^2` after U_G.
    static_ring: f64,
    uv_ring: f64,
    /// The static second loop `N_B/(k^2 (k - q)^2)` per unit `T_zm(4)`, and
    /// `N_B/((k^2 + m_IR^2)((k - q)^2 + m_IR^2))` after U_G.
    static_second: f64,
    uv_second: f64,
    /// [`static_second_loop`] at `m = 0` and, after U_B, at `m = m_IR`.
    second_taylor: (f64, f64),
    uv_second_taylor: (f64, f64),
}

impl MixedRingFactors {
    fn new(
        second: SecondVertices,
        [q, p, k]: [[f64; 3]; 3],
        beta: f64,
        mu: f64,
        ir_mass: f64,
    ) -> Self {
        let k_q = difference3(k, q);
        let vertices = second.value(k, q);
        Self {
            temperatures: 1.0 / (beta * beta),
            bubble: static_bubble_at_zero(p, beta, mu),
            uv_bubble: uv_static_bubble(p, q).0,
            static_ring: 1.0 / dot3(q, q).powi(2),
            uv_ring: uv_propagator(q, ir_mass).powi(2),
            static_second: vertices / (dot3(k, k) * dot3(k_q, k_q)),
            uv_second: vertices * uv_propagator(k, ir_mass) * uv_propagator(k_q, ir_mass),
            second_taylor: static_second_loop(second, k, q, 0.0),
            uv_second_taylor: static_second_loop(second, k, q, ir_mass),
        }
    }
}

/// One depth pin: an image of the mixed ring's S_A family and its closed
/// form.
struct MixedRingPin {
    name: &'static str,
    image: DirectResidueBranches,
    closed_form: fn(&MixedRingFactors) -> f64,
}

/// The two media of the depth pins and of the static-term and precision checks:
/// `T = 1` at `mu = 0`, `T = 2` at the quark chemical potential of `mu_B = 3`.
const PIN_MEDIA: [(f64, f64); 2] = [(1.0, 0.0), (0.5, QUARK_MU)];

/// An enclosing U acts on a frozen zero-mode
/// branch at the degree of its static cograph. The branch's order tag
/// cancels the degree `dod - j` of the co-scaled frozen vertex, so the U
/// keeps the cograph's own orders only, with the frozen vertex fixed and the
/// static lines massive. Every nonzero depth pin against its closed form,
/// from the fixture's `N_0 = N_B(k, k - q) = k.(k - q)`
/// (`N_0(k, 0) = k^2`) and `Pi_A^(0)` ([`MixedRingFactors`]), at three
/// points, `T = 1` and `T = 2`, `M = 1.3`, and the static-theory mass
/// `m = m_IR` of U_B and U_G equal to `M` and apart from it:
/// - G/A has degree `6 - 8 + 2 = 0`:
///   `U_G S_A = T^2 Pi_A^(0) N_0/((q^2+m^2)^2 (k^2+m^2) ((k-q)^2+m^2))`,
///   and `U_G U_A S_A` is the same with the own `U_A Pi_A^(0)` at `M`;
/// - U_B expands the static loop at its degree `+1`,
///   `U_B S_A = T^2 Pi_A^(0) [U B_0]/q^4` at `m`, after which G/A has degree 0:
///   `U_G U_B S_A = T^2 Pi_A^(0) [U B_0]/(q^2+m^2)^2`;
/// - the ring has degree `j - 1`, so U_G keeps the q-linear terms only:
///   `U_G S S_A = T^2 Pi_A^(0) B_0^(1)/(q^2+m^2)^2`,
///   `U_G S U_B S_A = T^2 Pi_A^(0) [U B_0]^(1)/(q^2+m^2)^2`, and on the
///   quark ring without vacuum subtraction
///   `U_G S_h F = T [Pi_A Pi_B]^(1)/(q^2+m^2)^2`.
///
/// The images of the vacuum-soft order are pinned by
/// `vacuum_soft_term_gets_its_ring_u`, and which images carry `mIR` by
/// `static_theory_mass_only_on_static_cographs`.
#[test]
fn enclosing_u_acts_on_frozen_s0_at_the_static_cograph_degree() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let root = empty_forest_node(graph);
    let host_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?;
    let s_a = host_a.s0()?;
    let u_a_s_a = host_a.own_u_on_s0(&ctx)?;
    let second = forest_node(graph, &[4, 5], ApproximationType::MUV, 2)?;
    let u_b_s_a = apply_taylor(
        &ctx,
        orientation,
        &second,
        &root,
        None,
        &canonical_lmb(graph, &second)?,
        &s_a,
    )?;
    // The factorized join {A, B}, which replays A first, and the union host
    // above it.
    let join = forest_node(graph, &[2, 3, 4, 5], ApproximationType::MUV, 4)?;
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let union_lmb = canonical_lmb(graph, &union)?;
    let s_s_a = apply_taylor(
        &ctx,
        orientation,
        &union,
        &host_a.node,
        None,
        &union_lmb,
        &s_a,
    )?;
    let s_u_b_s_a = apply_taylor(&ctx, orientation, &union, &join, None, &union_lmb, &u_b_s_a)?;
    let whole = whole_ring_node(graph)?;
    let whole_lmb = canonical_lmb(graph, &whole)?;
    let u_g = |given: &OwnedForestNode, sector: &DirectResidueBranches| {
        apply_taylor(&ctx, orientation, &whole, given, None, &whole_lmb, sector)
    };
    let pins = [
        MixedRingPin {
            name: "U_A S_A",
            image: u_a_s_a.clone(),
            closed_form: |f| f.temperatures * f.uv_bubble * f.static_ring * f.static_second,
        },
        MixedRingPin {
            name: "U_B S_A",
            image: u_b_s_a.clone(),
            closed_form: |f| {
                f.temperatures
                    * f.bubble
                    * f.static_ring
                    * (f.uv_second_taylor.0 + f.uv_second_taylor.1)
            },
        },
        MixedRingPin {
            name: "U_G S_A",
            image: u_g(&host_a.node, &s_a)?,
            closed_form: |f| f.temperatures * f.bubble * f.uv_ring * f.uv_second,
        },
        MixedRingPin {
            name: "U_G U_A S_A",
            image: u_g(&host_a.node, &u_a_s_a)?,
            closed_form: |f| f.temperatures * f.uv_bubble * f.uv_ring * f.uv_second,
        },
        MixedRingPin {
            name: "U_G U_B S_A",
            image: u_g(&join, &u_b_s_a)?,
            closed_form: |f| {
                f.temperatures
                    * f.bubble
                    * f.uv_ring
                    * (f.uv_second_taylor.0 + f.uv_second_taylor.1)
            },
        },
        MixedRingPin {
            name: "U_G S S_A",
            image: u_g(&union, &s_s_a)?,
            closed_form: |f| f.temperatures * f.bubble * f.uv_ring * f.second_taylor.1,
        },
        MixedRingPin {
            name: "U_G S U_B S_A",
            image: u_g(&union, &s_u_b_s_a)?,
            closed_form: |f| f.temperatures * f.bubble * f.uv_ring * f.uv_second_taylor.1,
        },
    ];
    // Every mismatch, so that a failure shows the footprint of a regression.
    let mut mismatches = Vec::new();
    // S on the frozen S_A takes A's branch over: no tag of A and no class
    // temperature of the second loop's mode remain, while the union's own
    // tags appear.
    let owner_a = Atom::var(host_a.node.subgraph().symbol());
    let owner_h = Atom::var(union.subgraph().symbol());
    let s_s_a_sum = s_s_a.resolved_scalars()?.factorized_sum();
    for (name, tag, expected) in [
        ("order tag of A", zero_mode_order(owner_a.clone()), false),
        (
            "zero-mode temperature of A",
            zero_mode_temperature(owner_a),
            false,
        ),
        (
            "class temperature T_zm(4)",
            class_zero_mode_temperature(EdgeIndex(4)),
            false,
        ),
        ("union's order tag", zero_mode_order(owner_h.clone()), true),
        (
            "union's zero-mode temperature",
            zero_mode_temperature(owner_h),
            true,
        ),
    ] {
        if s_s_a_sum.contains(tag.as_view()) != expected {
            mismatches.push(format!("S S_A: the {name} is present: {}", !expected));
        }
    }
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = |ir_mass| {
        kinematics
            .iter()
            .flat_map(|kinematics| {
                PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, ir_mass))
            })
            .collect::<Vec<_>>()
    };
    for ir_mass in STATIC_THEORY_MASSES {
        let points = points(ir_mass);
        for MixedRingPin {
            name,
            image,
            closed_form,
        } in &pins
        {
            let mut values = spatial.values(image, &points)?.into_iter();
            for (index, momenta) in LOOP_MOMENTA.into_iter().enumerate() {
                for (beta, mu) in PIN_MEDIA {
                    let context =
                        format!("{name}, point {index}, beta {beta}, mu {mu}, m_IR {ir_mass}");
                    let expected = closed_form(&MixedRingFactors::new(
                        spatial.second,
                        momenta,
                        beta,
                        mu,
                        ir_mass,
                    ));
                    assert!(expected.abs() > 0.0, "{context}: the closed form vanishes");
                    let actual = values.next().expect("one value per point");
                    if (actual - expected).abs() > 1e-8 * expected.abs() {
                        mismatches.push(format!("{context}: {actual} != {expected}"));
                    }
                }
            }
        }
    }

    // The quark ring's union host, without vacuum subtraction.
    let quark = SpatialRing::new(Ring::Quark)?;
    let graph = &quark.graph;
    let ctx = UVCtx::new(graph, &settings);
    let orientation = quark.cff.orientation(false);
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        quark.ring.region(graph, &[2, 3, 4, 5], false)?,
    )?;
    let s_h = apply_taylor(
        &ctx,
        orientation,
        &union,
        &empty_forest_node(graph),
        None,
        &canonical_lmb(graph, &union)?,
        &quark.root(false)?,
    )?;
    let whole = whole_ring_node(graph)?;
    let image = apply_taylor(
        &ctx,
        orientation,
        &whole,
        &union,
        None,
        &canonical_lmb(graph, &whole)?,
        &s_h,
    )?;
    for ir_mass in STATIC_THEORY_MASSES {
        let mut values = quark.values(&image, &points(ir_mass))?.into_iter();
        for (index, [q, p, k]) in LOOP_MOMENTA.into_iter().enumerate() {
            for (beta, mu) in PIN_MEDIA {
                let context =
                    format!("U_G S_h F, point {index}, beta {beta}, mu {mu}, m_IR {ir_mass}");
                let (c0, c1) = taylor_at_zero(|s| quark_ring_static_mode([q, p, k], s, beta, mu));
                let ring = uv_propagator(q, ir_mass).powi(2);
                let actual = values.next().expect("one value per point");
                // The q-linear coefficient is a difference of terms of the size
                // of c0, which sets the precision.
                if (actual - c1 * ring).abs() > 1e-7 * (c0.abs() + c1.abs()) * ring {
                    mismatches.push(format!("{context}: {actual} != {}", c1 * ring));
                }
            }
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} pins of the enclosing U on frozen branches fail:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
    Ok(())
}

/// A U whose static cograph on a frozen zero-mode
/// branch has negative degree annihilates it at the minimal depth, although
/// the co-scaled degree of the frozen vertex would leave a term. With the
/// three-gluon-like vertices `N_0 = (2k - q)^2`, each depending on every loop
/// through it as in the gluon ring, the static triangles D1/A (e0, e1, e4)
/// and D2/A (e0, e1, e5) have degree `3 - 6 + 2 = -1`:
/// `U_D1 S_A = U_D2 S_A = U_D U_A S_A = U_G U_D S_A = 0`, while without its
/// order tag the frozen vertex adds its degree two and U_D keeps a term. The
/// ring's order `j = 0` has degree -1: `U_G S_h = 0` for the quark ring's
/// union host through `J = 0`, and again not without the order tag.
#[test]
fn minimal_depth_annihilates_subcritical_cographs() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let spatial = SpatialRing::with_second_vertices(Ring::Mixed, SecondVertices::ThreeGluon)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let host_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?;
    let s_a = host_a.s0()?;
    let u_a_s_a = host_a.own_u_on_s0(&ctx)?;
    // The frozen vertex without its order tag, as the co-scaled depth of the
    // former finite-temperature S counted it.
    let without_order_tag = |family: &DirectResidueBranches, support: &SuBitGraph| {
        let tag = zero_mode_order(Atom::var(support.symbol()));
        family.map(|root| root.replace(tag.clone()).with(Atom::one()))
    };
    let co_scaled = without_order_tag(&s_a, host_a.node.subgraph());
    let whole = whole_ring_node(graph)?;
    let whole_lmb = canonical_lmb(graph, &whole)?;
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| {
            PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, UV_MASS))
        })
        .collect::<Vec<_>>();
    let assert_nonzero = |values: Vec<f64>, context: &str| {
        assert!(
            values.iter().all(|value| value.is_finite())
                && values.iter().any(|value| value.abs() > 0.0),
            "{context}: {values:?}"
        );
    };
    for (name, edges) in [("D1", [0, 1, 2, 3, 4]), ("D2", [0, 1, 2, 3, 5])] {
        let triangle = forest_node(graph, &edges, ApproximationType::MUV, 2)?;
        let lmb = canonical_lmb(graph, &triangle)?;
        let u_d = |sector: &DirectResidueBranches| {
            apply_taylor(
                &ctx,
                orientation,
                &triangle,
                &host_a.node,
                None,
                &lmb,
                sector,
            )
        };
        let u_d_s_a = u_d(&s_a)?;
        assert!(u_d_s_a.is_zero(), "U_{name} S_A must vanish");
        assert!(u_d(&u_a_s_a)?.is_zero(), "U_{name} U_A S_A must vanish");
        assert!(
            apply_taylor(
                &ctx,
                orientation,
                &whole,
                &triangle,
                None,
                &whole_lmb,
                &u_d_s_a
            )?
            .is_zero(),
            "U_G U_{name} S_A must vanish"
        );
        assert_nonzero(
            spatial.values(&u_d(&co_scaled)?, &points)?,
            &format!("U_{name} keeps a term of S_A without its order tag"),
        );
    }

    // The quark ring's union host through the order j = 0.
    let quark = SpatialRing::new(Ring::Quark)?;
    let graph = &quark.graph;
    let ctx = UVCtx::new(graph, &settings);
    let orientation = quark.cff.orientation(false);
    let region = static_region_through(&quark.ring.region(graph, &[2, 3, 4, 5], false)?, 0);
    let union = host_node(graph, &[2, 3, 4, 5], 4, region)?;
    let s_h = apply_taylor(
        &ctx,
        orientation,
        &union,
        &empty_forest_node(graph),
        None,
        &canonical_lmb(graph, &union)?,
        &quark.root(false)?,
    )?;
    assert!(!s_h.is_zero(), "S_h through j = 0 is nonzero");
    let whole = whole_ring_node(graph)?;
    let whole_lmb = canonical_lmb(graph, &whole)?;
    let u_g = |sector: &DirectResidueBranches| {
        apply_taylor(&ctx, orientation, &whole, &union, None, &whole_lmb, sector)
    };
    assert!(u_g(&s_h)?.is_zero(), "U_G S_h at j = 0 must vanish");
    assert_nonzero(
        quark.values(&u_g(&without_order_tag(&s_h, union.subgraph()))?, &points)?,
        "U_G keeps a term of S_h at j = 0 without its order tag",
    );
    Ok(())
}

/// The zero-mode operation commutes with the hard U's.
/// - On the unfrozen quark ring, `S U_A F` is the closed form of `U_A S F`:
///   U_A of the static mode's Taylor terms
///   `(T/q^4)[Pi_A^(0) (Pi_B^(0) + Pi_B^(1)) + Pi_A^(1) Pi_B^(0)]` takes
///   each `Pi_A^(i)` to its 4D U ([`uv_static_bubble`]). U_A sees q only
///   as `t q`, so its depth two and the order `J = 1` commute. The
///   vacuum-soft order is `-U_A Pi_A^vac(0) Pi_B^vac(0)/(4 |q|^3)`.
/// - On the mixed ring, `S U_B S_A = U_B S_A`: `[U B_0]` has degree one in
///   q, all of which the union host's order `J = 1` keeps, and its
///   vacuum-soft order reads the projected rows off as zero.
#[test]
fn zero_mode_s0_commutes_with_hard_u() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let settings = UVgenerationSettings::default();
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    let points = kinematics
        .iter()
        .flat_map(|kinematics| {
            PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, UV_MASS))
        })
        .collect::<Vec<_>>();

    let quark = SpatialRing::new(Ring::Quark)?;
    let graph = &quark.graph;
    let ctx = UVCtx::new(graph, &settings);
    let orientation = quark.cff.orientation(vacuum_subtraction);
    let first = forest_node(graph, &[2, 3], ApproximationType::MUV, 2)?;
    let u_a_f = apply_taylor(
        &ctx,
        orientation,
        &first,
        &empty_forest_node(graph),
        None,
        &canonical_lmb(graph, &first)?,
        &quark.root(vacuum_subtraction)?,
    )?;
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        quark
            .ring
            .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let s_u_a_f = apply_taylor(
        &ctx,
        orientation,
        &union,
        &first,
        None,
        &canonical_lmb(graph, &union)?,
        &u_a_f,
    )?;
    let mut parts = zero_mode_temperature_parts(&quark, &s_u_a_f, &points)?.into_iter();
    for (index, [q, p, k]) in LOOP_MOMENTA.into_iter().enumerate() {
        let (a0, a1) = uv_static_bubble(p, q);
        for (beta, mu) in PIN_MEDIA {
            let context = format!("S U_A F, point {index}, beta {beta}, mu {mu}");
            let (b0, b1) = taylor_at_zero(|s| static_quark_bubble(k, scaled3(q, s), beta, mu));
            let static_ring = 1.0 / (beta * dot3(q, q).powi(2));
            let terms = [a0 * b0, a0 * b1, a1 * b0];
            let expected = static_ring * terms.iter().sum::<f64>();
            let scale = static_ring * terms.iter().map(|term| term.abs()).sum::<f64>();
            let (static_part, vacuum_soft) = parts.next().expect("one value per point");
            assert!(
                (static_part - expected).abs() <= 1e-7 * scale,
                "{context}: {static_part} != {expected}"
            );
            assert_close(
                vacuum_soft,
                -a0 * vacuum_static_bubble(k) / (4.0 * norm3(q).powi(3)),
                1e-9,
                &format!("{context}: vacuum-soft order"),
            );
        }
    }

    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let s_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?.s0()?;
    let second = forest_node(graph, &[4, 5], ApproximationType::MUV, 2)?;
    let u_b_s_a = apply_taylor(
        &ctx,
        orientation,
        &second,
        &empty_forest_node(graph),
        None,
        &canonical_lmb(graph, &second)?,
        &s_a,
    )?;
    let join = forest_node(graph, &[2, 3, 4, 5], ApproximationType::MUV, 4)?;
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let s_u_b_s_a = apply_taylor(
        &ctx,
        orientation,
        &union,
        &join,
        None,
        &canonical_lmb(graph, &union)?,
        &u_b_s_a,
    )?;
    let after = spatial.values(&s_u_b_s_a, &points)?;
    let before = spatial.values(&u_b_s_a, &points)?;
    for (index, (after, before)) in after.into_iter().zip(before).enumerate() {
        let context = format!("S U_B S_A, point {index}");
        assert!(before.abs() > 0.0, "{context}: U_B S_A vanishes");
        assert_close(after, before, 1e-10, &context);
    }
    Ok(())
}

/// The vacuum-soft order of the quark ring's union host gets the
/// ring's U. The ring `|q|^-3` has degree 0, so U_G keeps its leading term
/// with the frozen vacuum vertex fixed:
/// `U_G S^(3) X = -V_A V_B (q^2 + m^2)^(-3/2)/4`, nonzero, with the
/// static-theory mass `m = m_IR` of U_G equal to `M` and apart from it, for the
/// chains `X = F`, `U_A F` and `U_B F`, with `V = Pi^vac(0)` for a bubble the
/// chain leaves alone and `V = U Pi^vac(0)` ([`uv_static_bubble`]) at `M` for
/// one its hard U replaced. Their static orders keep the q-linear terms only,
/// as in `enclosing_u_acts_on_frozen_s0_at_the_static_cograph_degree`:
/// `T [X_A X_B]^(1)/(q^2 + m^2)^2` with the thermal `Pi(q)` or the 4D U of
/// each bubble through first order in q.
#[test]
fn vacuum_soft_term_gets_its_ring_u() -> Result<()> {
    test_initialise()?;
    let quark = SpatialRing::new(Ring::Quark)?;
    let graph = &quark.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = quark.cff.orientation(true);
    let root = empty_forest_node(graph);
    let integrand = quark.root(true)?;
    let first = forest_node(graph, &[2, 3], ApproximationType::MUV, 2)?;
    let second = forest_node(graph, &[4, 5], ApproximationType::MUV, 2)?;
    let hard_u = |node: &OwnedForestNode| {
        apply_taylor(
            &ctx,
            orientation,
            node,
            &root,
            None,
            &canonical_lmb(graph, node)?,
            &integrand,
        )
    };
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        quark.ring.region(graph, &[2, 3, 4, 5], true)?,
    )?;
    let union_lmb = canonical_lmb(graph, &union)?;
    let whole = whole_ring_node(graph)?;
    let whole_lmb = canonical_lmb(graph, &whole)?;
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    // Each chain with the bubbles its hard U replaced.
    let chains = [
        ("U_G S F", &root, integrand.clone(), [false, false]),
        ("U_G S U_A F", &first, hard_u(&first)?, [true, false]),
        ("U_G S U_B F", &second, hard_u(&second)?, [false, true]),
    ];
    for (name, given, sector, hard) in chains {
        let s = apply_taylor(&ctx, orientation, &union, given, None, &union_lmb, &sector)?;
        let image = apply_taylor(&ctx, orientation, &whole, &union, None, &whole_lmb, &s)?;
        for ir_mass in STATIC_THEORY_MASSES {
            let points = kinematics
                .iter()
                .flat_map(|kinematics| {
                    PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, ir_mass))
                })
                .collect::<Vec<_>>();
            let mut parts = zero_mode_temperature_parts(&quark, &image, &points)?.into_iter();
            for ((index, [q, p, k]), (beta, mu)) in LOOP_MOMENTA
                .into_iter()
                .enumerate()
                .cartesian_product(PIN_MEDIA)
            {
                let context =
                    format!("{name}, point {index}, beta {beta}, mu {mu}, m_IR {ir_mass}");
                // (vacuum vertex, value at q = 0, q-linear term) per bubble.
                let [(vertex_a, a0, a1), (vertex_b, b0, b1)] =
                    [(p, hard[0]), (k, hard[1])].map(|(l, replaced)| {
                        if replaced {
                            let (value, slope) = uv_static_bubble(l, q);
                            (value, value, slope)
                        } else {
                            let (value, slope) =
                                taylor_at_zero(|s| static_quark_bubble(l, scaled3(q, s), beta, mu));
                            (vacuum_static_bubble(l), value, slope)
                        }
                    });
                let ring = uv_propagator(q, ir_mass);
                let (static_part, vacuum_soft) = parts.next().expect("one value per point");
                let expected = -vertex_a * vertex_b * ring.powf(1.5) / 4.0;
                assert!(expected.abs() > 0.0, "{context}: the vertex vanishes");
                assert_close(vacuum_soft, expected, 1e-9, &format!("{context}: a = 3"));
                let expected = (a0 * b1 + a1 * b0) * ring.powi(2) / beta;
                let scale = (a0 * b0).abs() + (a0 * b1).abs() + (a1 * b0).abs();
                assert!(
                    (static_part - expected).abs() <= 1e-7 * scale * ring.powi(2) / beta,
                    "{context}: a = 4: {static_part} != {expected}"
                );
            }
        }
    }
    Ok(())
}

/// Without vacuum subtraction one whole-graph U acts on the
/// unfrozen tadpole ring and on its union host's frozen S in one call, at
/// the full 4D depth on the first and at the minimal depth on the second:
/// `U_G(F + S_h F) = U_G F + U_G S_h F`. `U_G F` is the ordinary 4D U, of
/// degree four: nonzero without vacuum subtraction, independent of `T` and
/// `mu`, and equal to the U of the same ring's vacuum CFF, at the UV mass `M`.
/// And `U_G S_h F = T Pi_A^(1)(q; p) Pi_T(k)/(q^2 + m_IR^2)^2` with the exact
/// tadpole `Pi_T(k) = T sum 1/(w^2 + k^2)` of the four-gluon vertex: only
/// the q-linear term of `S_h F = (T/q^4) T_q^(1)[Pi_A] Pi_T` reaches the
/// ring's degree 0. The one call gives the rows of `S_h F` the static-theory
/// mass `m_IR` and those of `F` the mass `M`, with `m_IR` equal to `M` and
/// apart from it.
#[test]
fn vs_false_sector_mixes_full_and_minimal_depth() -> Result<()> {
    test_initialise()?;
    let tadpole = SpatialRing::new(Ring::Tadpole)?;
    let graph = &tadpole.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = tadpole.cff.orientation(false);
    let root = empty_forest_node(graph);
    let integrand = tadpole.root(false)?;
    let union = host_node(
        graph,
        &[2, 3, 4],
        4,
        tadpole.ring.region(graph, &[2, 3, 4], false)?,
    )?;
    let s_h = apply_taylor(
        &ctx,
        orientation,
        &union,
        &root,
        None,
        &canonical_lmb(graph, &union)?,
        &integrand,
    )?;
    let whole = whole_ring_node(graph)?;
    let whole_lmb = canonical_lmb(graph, &whole)?;
    let u_g = |sector: &DirectResidueBranches| {
        apply_taylor(&ctx, orientation, &whole, &root, None, &whole_lmb, sector)
    };
    let mixed = u_g(&integrand.zip_add(&s_h)?)?;
    let unfrozen = u_g(&integrand)?;
    let frozen = u_g(&s_h)?;
    // The same U of the ring's vacuum CFF.
    let vacuum = SpatialRing::in_medium(Ring::Tadpole, tadpole.second, MediumMode::Vacuum)?;
    let vacuum_whole = whole_ring_node(&vacuum.graph)?;
    let vacuum_u = apply_taylor(
        &UVCtx::new(&vacuum.graph, &settings),
        vacuum.cff.orientation(false),
        &vacuum_whole,
        &empty_forest_node(&vacuum.graph),
        None,
        &canonical_lmb(&vacuum.graph, &vacuum_whole)?,
        &vacuum.root(false)?,
    )?;
    // With the static-theory mass apart from M, U_G F keeps M and U_G S_h F
    // takes m_IR in the same call.
    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    for ir_mass in STATIC_THEORY_MASSES {
        let points = kinematics
            .iter()
            .flat_map(|kinematics| {
                PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, ir_mass))
            })
            .collect::<Vec<_>>();
        let [mixed, unfrozen, frozen] =
            [&mixed, &unfrozen, &frozen].map(|image| tadpole.values(image, &points));
        let mut values = mixed?
            .into_iter()
            .zip(unfrozen?)
            .zip(frozen?)
            .zip(vacuum.values(&vacuum_u, &points)?);
        for ((index, [q, p, k]), (beta, mu)) in LOOP_MOMENTA
            .into_iter()
            .enumerate()
            .cartesian_product(PIN_MEDIA)
        {
            let context = format!("point {index}, beta {beta}, mu {mu}, m_IR {ir_mass}");
            let (((mixed, unfrozen), frozen), vacuum_u) =
                values.next().expect("one value per point");
            assert!(unfrozen.abs() > 0.0, "{context}: U_G F vanishes");
            assert_close(unfrozen, vacuum_u, 1e-10, &format!("{context}: U_G F"));
            assert_close(mixed, unfrozen + frozen, 1e-10, &context);
            let (a0, a1) = taylor_at_zero(|s| static_quark_bubble(p, scaled3(q, s), beta, mu));
            let factor = bose_tadpole(beta, norm3(k)) * uv_propagator(q, ir_mass).powi(2) / beta;
            assert!(
                (frozen - a1 * factor).abs() <= 1e-7 * (a0.abs() + a1.abs()) * factor.abs(),
                "{context}: U_G S_h F: {frozen} != {}",
                a1 * factor
            );
        }
    }
    Ok(())
}

/// Which U images carry the static-theory mass `mIR`, and which the UV mass
/// `mUV`, on the mixed ring without vacuum subtraction, with host A and the
/// union host. A U whose support holds a static line of a completed zero-mode
/// branch is a UV counterterm of the static theory on the branch's rows: U_B
/// on the soft lines `{4,5}` of A gives `U_B S_A` and `U_G U_B S_A` `mIR`
/// alone, and U_G, holding every soft line, gives `mIR` to `U_G S_A`,
/// `U_G U_A S_A` and `U_G S S_A`, with `mUV` only where the own U_A introduced
/// it, as in `U_B U_A S_A`. The own U_A on its branch, and U_A, U_B and U_G on
/// the unfrozen ring (`U_B U_A F`, `U_G U_A F`), keep `mUV` alone. On the
/// host's whole operation `H_A F = S_A + U_A F - U_A S_A`, U_B and U_G give
/// `mIR` to the rows of `S_A` and `U_A S_A` and `mUV` to those of `U_A F`:
/// with `m_IR` equal to `M` and apart from it, each image is the sum of its
/// images of the three parts.
#[test]
fn static_theory_mass_only_on_static_cographs() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = false;
    let spatial = SpatialRing::new(Ring::Mixed)?;
    let graph = &spatial.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    let orientation = spatial.cff.orientation(vacuum_subtraction);
    let root = empty_forest_node(graph);
    let integrand = spatial.root(vacuum_subtraction)?;
    let host_a = MixedRingHostA::new(&spatial, vacuum_subtraction)?;
    let s_a = host_a.s0()?;
    let u_a_s_a = host_a.own_u_on_s0(&ctx)?;
    let first = forest_node(graph, &[2, 3], ApproximationType::MUV, 2)?;
    let second = forest_node(graph, &[4, 5], ApproximationType::MUV, 2)?;
    let union = host_node(
        graph,
        &[2, 3, 4, 5],
        4,
        spatial
            .ring
            .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
    )?;
    let whole = whole_ring_node(graph)?;
    let u = |node: &OwnedForestNode, given: &OwnedForestNode, sector: &DirectResidueBranches| {
        apply_taylor(
            &ctx,
            orientation,
            node,
            given,
            None,
            &canonical_lmb(graph, node)?,
            sector,
        )
    };
    let u_a_f = u(&first, &root, &integrand)?;
    let u_b_s_a = u(&second, &root, &s_a)?;
    let h_a_f = u(&host_a.node, &root, &integrand)?;
    let u_b_u_a_s_a = u(&second, &host_a.node, &u_a_s_a)?;
    let u_b_u_a_f = u(&second, &first, &u_a_f)?;
    let u_b_h_a_f = u(&second, &host_a.node, &h_a_f)?;
    let u_g_s_a = u(&whole, &host_a.node, &s_a)?;
    let u_g_u_a_s_a = u(&whole, &host_a.node, &u_a_s_a)?;
    let u_g_u_a_f = u(&whole, &first, &u_a_f)?;
    let u_g_h_a_f = u(&whole, &host_a.node, &h_a_f)?;
    // (image, carries mIR, carries mUV)
    let images = [
        ("U_A S_A", &u_a_s_a, false, true),
        ("U_B S_A", &u_b_s_a, true, false),
        ("U_B U_A S_A", &u_b_u_a_s_a, true, true),
        ("U_A F", &u_a_f, false, true),
        ("U_B U_A F", &u_b_u_a_f, false, true),
        ("U_B H_A F", &u_b_h_a_f, true, true),
        ("U_G F", &u(&whole, &root, &integrand)?, false, true),
        ("U_G U_A F", &u_g_u_a_f, false, true),
        ("U_G S_A", &u_g_s_a, true, false),
        ("U_G U_A S_A", &u_g_u_a_s_a, true, true),
        ("U_G U_B S_A", &u(&whole, &second, &u_b_s_a)?, true, false),
        (
            "U_G S S_A",
            &u(&whole, &union, &u(&union, &host_a.node, &s_a)?)?,
            true,
            false,
        ),
        ("U_G H_A F", &u_g_h_a_f, true, true),
    ];
    let mut mismatches = Vec::new();
    for (name, image, ir, uv) in images {
        let sum = image.resolved_scalars()?.factorized_sum();
        let carries = |masses: [Symbol; 2]| masses.iter().any(|mass| sum.contains_symbol(*mass));
        let carried = (
            carries([GS.m_ir_vacuum, GS.m_ir_expansion]),
            carries([GS.m_uv_vacuum, GS.m_uv_expansion]),
        );
        if carried != (ir, uv) {
            mismatches.push(format!("{name}: (mIR, mUV) {carried:?} != {:?}", (ir, uv)));
        }
    }
    assert!(
        mismatches.is_empty(),
        "the static-theory mass reaches the wrong U images:\n{}",
        mismatches.join("\n")
    );

    let kinematics = LOOP_MOMENTA.map(RingKinematics::new);
    // (U, its images of H_A F, S_A, U_A S_A and U_A F)
    let operations = [
        ("U_B", [&u_b_h_a_f, &u_b_s_a, &u_b_u_a_s_a, &u_b_u_a_f]),
        ("U_G", [&u_g_h_a_f, &u_g_s_a, &u_g_u_a_s_a, &u_g_u_a_f]),
    ];
    for ir_mass in STATIC_THEORY_MASSES {
        let points = kinematics
            .iter()
            .flat_map(|kinematics| {
                PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, ir_mass))
            })
            .collect::<Vec<_>>();
        for (name, images) in operations {
            let [whole_operation, s, u_on_s, ordinary] =
                images.map(|image| spatial.values(image, &points));
            for (index, (((whole_operation, s), u_on_s), ordinary)) in whole_operation?
                .into_iter()
                .zip(s?)
                .zip(u_on_s?)
                .zip(ordinary?)
                .enumerate()
            {
                assert!(
                    ordinary.abs() > 0.0 && s.abs() > 0.0,
                    "point {index}: {name} U_A F or {name} S_A vanishes"
                );
                assert_close(
                    whole_operation,
                    s + ordinary - u_on_s,
                    1e-10,
                    &format!("{name} H_A F, point {index}, m_IR {ir_mass}"),
                );
            }
        }
    }
    Ok(())
}

/// The ordinary U of a support is a UV counterterm of a static theory on a
/// row iff the support holds a line of the soft set of the row's outermost
/// completed zero-mode branch, the soft set of a branch being the loop lines
/// outside it. On the mixed ring the soft set of host A is `{0,1,4,5}`: the
/// whole graph, the triangle D1 `{0,1,2,3,4}`, the bubble B `{4,5}`, the
/// static pair `{0,1}` and the soft set itself all hold some of it, and only
/// the own branch `{2,3}`, which its own U alone revisits, holds none. A row
/// that also holds the branch of the union host `h = {2,3,4,5}` around A is
/// decided by h alone, whose zero-mode operation expanded A's further soft
/// lines `{4,5}` as hard ones: B holds soft lines of A but not of h. A row
/// with the incomparable branches of A and of B, whose soft set is
/// `{0,1,2,3}`, is static for a support meeting both soft sets, and an error
/// for B, which meets A's soft set but not its own.
#[test]
fn static_cograph_needs_a_soft_line() -> Result<()> {
    test_initialise()?;
    let graph = Ring::Mixed.topology()?;
    let loop_lines = graph.uncut_loop_lines();
    assert_eq!(loop_lines, graph.edge_subgraph(&[0, 1, 2, 3, 4, 5]));
    let host_a = graph.edge_subgraph(&[2, 3]);
    let bubble_b = graph.edge_subgraph(&[4, 5]);
    let union_h = graph.edge_subgraph(&[2, 3, 4, 5]);
    let static_pair = graph.edge_subgraph(&[0, 1]);
    let soft_set_a = graph.edge_subgraph(&[0, 1, 4, 5]);
    let expands = |support: &SuBitGraph, owners: &[&SuBitGraph]| {
        SoftBranches::new(&graph, support, Local3DDeformation::Ordinary).expands_static_cograph(
            &owners
                .iter()
                .map(|owner| (*owner).clone())
                .collect::<BTreeSet<_>>(),
            &loop_lines,
        )
    };
    for support in [
        &loop_lines,
        &graph.edge_subgraph(&[0, 1, 2, 3, 4]),
        &bubble_b,
        &static_pair,
        &soft_set_a,
    ] {
        assert!(expands(support, &[&host_a])?, "{}", support.string_label());
    }
    assert!(!expands(&host_a, &[&host_a])?);
    assert!(!expands(&loop_lines, &[])?);
    assert!(expands(&loop_lines, &[&host_a, &union_h])?);
    assert!(expands(&static_pair, &[&host_a, &union_h])?);
    assert!(!expands(&bubble_b, &[&host_a, &union_h])?);
    assert!(expands(&soft_set_a, &[&host_a, &bubble_b])?);
    let error = expands(&bubble_b, &[&host_a, &bubble_b])
        .expect_err("a U holding soft lines of one branch but none of another is refused")
        .to_string();
    assert!(error.contains("unsupported zero-mode forest"), "{error}");
    Ok(())
}

/// `l` with its component along the unit vector `direction` set to `t`.
fn with_component(l: [f64; 3], direction: [f64; 3], t: f64) -> [f64; 3] {
    let shift = t - dot3(l, direction);
    [0, 1, 2].map(|axis| l[axis] + shift * direction[axis])
}

/// The bits of the reference evaluation in the Landau-pair check, far beyond
/// the cancellation of its per-variant poles at `p.q ~ 1e-8`.
const LANDAU_REFERENCE_BITS: u32 = 300;
/// The double-precision error allowed against that reference; without the
/// symbolic cancellation it was 2.8e-10 already at `p.q ~ 1e-2`.
const MAX_LANDAU_DOUBLE_RELATIVE_ERROR: f64 = 1e-11;

/// A zero-mode chart cancels the per-variant Landau-pair poles of its
/// coefficients symbolically. Expanding an energy difference such as
/// `|p| - |p - q|` divides each Taylor coefficient of a CFF variant by `p.q`,
/// and those of the second bubble by `k.q`; the poles `1/(p.q)^n` and
/// `1/(k.q)^n` cancel only in the sum of the variants. With that sum left to
/// the evaluation, double precision lost two digits per decade of
/// `cos(p, q)` on the quark ring's union host, about `1e-1` at `1e-7`, while
/// the bare integrand stayed flat. Approaching `p.q = 0` and `k.q = 0` at
/// fixed `q`, double precision must now agree with a high-precision
/// evaluation of the same expression: on the union host of the quark ring
/// and on its chain after U_A, where the host's numerator multiplies the
/// jets of an earlier one, with and without vacuum subtraction, and on
/// host A of the mixed ring.
#[test]
fn zero_mode_chart_cancels_landau_pair_poles() -> Result<()> {
    test_initialise()?;
    let [q, p, k] = LOOP_MOMENTA[0];
    let q_hat = scaled3(q, 1.0 / norm3(q));
    let steps = [1e-2, 1e-4, 1e-6, 1e-8];
    // The loci p.q = 0 (j) and k.q = 0 (h).
    let kinematics = steps
        .iter()
        .flat_map(|&t| {
            [
                RingKinematics::new([q, with_component(p, q_hat, t), k]),
                RingKinematics::new([q, p, with_component(k, q_hat, t)]),
            ]
        })
        .collect::<Vec<_>>();
    let points = kinematics
        .iter()
        .flat_map(|kinematics| {
            PIN_MEDIA.map(|(beta, mu)| uv_ring_point(kinematics, beta, mu, UV_MASS))
        })
        .collect::<Vec<_>>();
    let assert_precise =
        |spatial: &SpatialRing, branches: &DirectResidueBranches, name: &str| -> Result<()> {
            let atom = normalize_local_operation(
                &spatial.graph,
                &branches.resolved_scalars()?.factorized_sum(),
            );
            let double = evaluate_points(&spatial.graph, &atom, &points)?;
            let precise =
                evaluate_points_precise(&spatial.graph, &atom, &points, LANDAU_REFERENCE_BITS)?;
            for (index, (double, precise)) in double.into_iter().zip(precise).enumerate() {
                let t = steps[index / (2 * PIN_MEDIA.len())];
                let locus = ["p.q", "k.q"][(index / PIN_MEDIA.len()) % 2];
                assert!(precise != 0.0, "{name}: vanishes at {locus} ~ {t}");
                let error = ((double - precise) / precise).abs();
                assert!(
                    error <= MAX_LANDAU_DOUBLE_RELATIVE_ERROR,
                    "{name}: double precision is off by {error:e} at {locus} ~ {t} \
                 ({double} != {precise})"
                );
            }
            Ok(())
        };
    let quark = SpatialRing::new(Ring::Quark)?;
    let graph = &quark.graph;
    let settings = UVgenerationSettings::default();
    let ctx = UVCtx::new(graph, &settings);
    for vacuum_subtraction in [false, true] {
        let host = host_node(
            graph,
            &[2, 3, 4, 5],
            4,
            quark
                .ring
                .region(graph, &[2, 3, 4, 5], vacuum_subtraction)?,
        )?;
        let orientation = quark.cff.orientation(vacuum_subtraction);
        let root = quark.root(vacuum_subtraction)?;
        let union_host = apply_taylor(
            &ctx,
            orientation,
            &host,
            &empty_forest_node(graph),
            None,
            &canonical_lmb(graph, &host)?,
            &root,
        )?;
        assert_precise(
            &quark,
            &union_host,
            &format!("S_h F, vacuum subtraction {vacuum_subtraction}"),
        )?;
        // After U_A the host multiplies the jets of A's numerator by B's.
        let first = forest_node(graph, &[2, 3], ApproximationType::MUV, 2)?;
        let u_a_f = apply_taylor(
            &ctx,
            orientation,
            &first,
            &empty_forest_node(graph),
            None,
            &canonical_lmb(graph, &first)?,
            &root,
        )?;
        let s_u_a_f = apply_taylor(
            &ctx,
            orientation,
            &host,
            &first,
            None,
            &canonical_lmb(graph, &host)?,
            &u_a_f,
        )?;
        assert_precise(
            &quark,
            &s_u_a_f,
            &format!("S_h U_A F, vacuum subtraction {vacuum_subtraction}"),
        )?;
    }
    let mixed = SpatialRing::new(Ring::Mixed)?;
    assert_precise(&mixed, &MixedRingHostA::new(&mixed, false)?.s0()?, "S0_A F")?;
    Ok(())
}

/// The integrated prefixes of a ring whose add-backs carry the static mode
/// of the q-gluons, which the union host's zero mode on the contracted
/// support cancels: the first bubble A, the second loop B and their
/// factorized join.
#[derive(Clone, Copy, Debug)]
enum RingPrefix {
    First,
    Second,
    Join,
}

impl RingPrefix {
    const ALL: [Self; 3] = [Self::First, Self::Second, Self::Join];
}

/// The soft scalings `s` of the ring momentum `q` at which the add-back
/// pins compare the remainder `|q s|^4 (1 - S0_h) ι_X` with the add-back
/// `|q s|^4 ι_X`. The static mode and its first order in `q` cancel, so
/// their ratio vanishes like `s^2`: it stays below `100 s^2` and falls
/// between the two scalings by more than a remainder linear in `q` would.
const ADD_BACK_SCALINGS: [f64; 2] = [1e-2, 3e-3];

/// Asserts that the add-back `|q s|^4 ι_X` and its zero-mode image,
/// `values(s)`, cancel the static mode, so that their remainder vanishes
/// like `s^2` against the add-back ([`ADD_BACK_SCALINGS`]).
fn assert_static_mode_cancels(
    mut values: impl FnMut(f64) -> Result<(Complex<f64>, Complex<f64>)>,
    context: &str,
) -> Result<()> {
    let mut ratios = Vec::new();
    for s in ADD_BACK_SCALINGS {
        let (add_back, image) = values(s)?;
        let scale = add_back.re.hypot(add_back.im);
        assert!(
            scale > 0.0,
            "{context}, s {s}: the add-back has no static mode"
        );
        let ratio = (add_back.re + image.re).hypot(add_back.im + image.im) / scale;
        assert!(
            ratio <= 100.0 * s * s,
            "{context}, s {s}: the zero-mode image does not cancel the static mode of the \
             add-back, image / add-back = {:?}",
            image / add_back
        );
        ratios.push(ratio);
    }
    let (coarse, fine) = (ratios[0], ratios[1]);
    let linear = ADD_BACK_SCALINGS[0] / ADD_BACK_SCALINGS[1];
    assert!(
        coarse >= 1.5 * linear * fine,
        "{context}: the remainder {coarse:e}, {fine:e} vanishes no faster than q"
    );
    Ok(())
}

/// Every integrated add-back cancels the static mode of the q-gluons with
/// its image under the union host's zero mode on the contracted support,
/// `(1 - S0_h) ι_X`. On the quark ring for the first bubble, the second
/// bubble and their join, and on the tadpole ring for the tadpole, the
/// remainder vanishes like `q^2` against the add-back as `q -> 0`, with the
/// vertex of [`ring_vertex`], whose constant holds the vacuum mass both as
/// `mUV` and as the expansion mass `mUVexp`, like vakint's counterterms. The
/// image is compared with the add-back it subtracts, not with a projection
/// of it: a chart that held `mUVexp` fixed dropped that part of the static
/// mode, which made the image three times the add-back here.
#[test]
fn integrated_add_backs_cancel_their_static_mode() -> Result<()> {
    test_initialise()?;
    let vacuum_subtraction = true;
    let settings = UVgenerationSettings::default();
    let cases = [Ring::Quark]
        .into_iter()
        .flat_map(|ring| RingPrefix::ALL.map(|prefix| (ring, prefix)))
        .chain([(Ring::Tadpole, RingPrefix::Second)]);
    for (ring, prefix) in cases {
        let mut spatial = SpatialRing::new(ring)?;
        let graph = &spatial.graph;
        let union = [[2, 3].as_slice(), ring.second_loop_edges()].concat();
        let host = host_node(
            graph,
            &union,
            4,
            ring.region(graph, &union, vacuum_subtraction)?,
        )?;
        let vertex = ring_vertex(RING_VERTEX);
        // The prefix's node and vertex, the support it leaves to the host and
        // the vertices of that support, which the host's operation attaches.
        let (node, vertex, support, vertices) = match prefix {
            RingPrefix::First => (
                forest_node(graph, &[2, 3], ApproximationType::MUV, 2)?,
                vertex,
                ring.second_loop_edges().to_vec(),
                spatial.second.numerator(),
            ),
            RingPrefix::Second => (
                forest_node(graph, ring.second_loop_edges(), ApproximationType::MUV, 2)?,
                vertex,
                vec![2, 3],
                SpatialRing::vertex(2, 3),
            ),
            RingPrefix::Join => (
                forest_node(graph, &union, ApproximationType::MUV, 4)?,
                vertex.pow(2),
                Vec::new(),
                Atom::one(),
            ),
        };
        let [prefix_sector] = <[DirectSector; 1]>::try_from(localized_prefix(
            &mut spatial,
            &vertex,
            &node,
            vacuum_subtraction,
        )?)
        .map_err(|_| eyre!("{ring:?}: an ordinary prefix has one sector"))?;
        let graph = &spatial.graph;
        let ctx = UVCtx::new(graph, &settings);
        let orientation = spatial.cff.orientation(vacuum_subtraction);
        let support = graph.edge_subgraph(&support);
        let (_, chart) = coordinate_lmb(&ctx, &host, &node, None, &[], &support)?;
        let add_back = prefix_sector.active.multiply_key_mapped(
            orientation,
            graph,
            &vertices,
            DirectResidueBranches::numerator_scope(),
        )?;
        let image = -apply_taylor(
            &ctx,
            orientation,
            &host,
            &node,
            Some(support),
            &chart,
            &add_back,
        )?;
        for (index, [q, _, _]) in LOOP_MOMENTA.into_iter().enumerate() {
            assert_static_mode_cancels(
                |s| {
                    let kinematics = scaled_ring(index, s);
                    let point = uv_ring_point(&kinematics, 1.0, QUARK_MU, IR_MASS);
                    let [add_back, image] = [&add_back, &image]
                        .map(|family| spatial.values(family, slice::from_ref(&point)));
                    let q4 = (norm3(q) * s).powi(4);
                    Ok((
                        Complex::new(q4 * add_back?[0], 0.0),
                        Complex::new(q4 * image?[0], 0.0),
                    ))
                },
                &format!("{ring:?}, {prefix:?}, point {index}"),
            )?;
        }
    }
    Ok(())
}

/// The finite counterterm of the plain U of `edges` with degree `dod`,
/// integrated by vakint as the zero-temperature lane integrates a node on
/// the root (`Forests::compute_4d_for_node`).
fn vakint_integrated(
    graph: &Graph,
    settings: &UVgenerationSettings,
    vakint_settings: &VakintSettings,
    edges: &[usize],
    dod: i32,
) -> Result<IntegratedCts> {
    let node = forest_node(graph, edges, ApproximationType::MUV, dod)?;
    let root = empty_forest_node(graph);
    let ctx = UVCtx::new(graph, settings);
    let full = Full4dCts::recursion_input(
        &Local4dCts::root(),
        &IntegratedCts::root(),
        root.renormalization_scheme(),
        true,
        root.lmb(),
    )?;
    let local = local_4d::uv_limit(&full, &ctx, &node, &root, &node, &root)?;
    Integrated::new(utils::vakint()?, vakint_settings).run(&local, &ctx, &node, &root, &node, &root)
}

/// A final integrand summed over its residue rows as a function of the
/// ring's point: every selector and the localization kernel at one, the
/// closed numerators contracted, the local tags and masses at their final
/// values and the quark massless.
fn final_ring_atom(graph: &Graph, integrands: FinalIntegrands) -> Result<Atom> {
    let atom = integrands
        .into_integrands()
        .resolved()?
        .iter()
        .map(|(_, atom)| atom)
        .sum::<Atom>()
        .replace(function!(OrientationID::symbol(), W_.a_))
        .with(Atom::one())
        .replace(function!(GS.localizing_integrand, W_.a_))
        .with(Atom::one());
    Ok(normalize_local_operation(
        graph,
        &contract_numerators(&atom)?,
    ))
}

/// The add-back of an integrated prefix on a physical ring and its image
/// ([`production_add_back_and_image`]).
struct ProductionAddBack {
    graph: Graph,
    /// The prefix's finite counterterm, as the add-back attaches it.
    finite: Atom,
    add_back: Atom,
    image: Atom,
}

/// The add-back of `prefix` on the physical `ring`, as final assembly builds
/// it for the prefix's node, and its image under the union host's zero mode,
/// the integrated sectors that `run_integrated` replays for the host's node,
/// both as final integrands ([`final_ring_atom`]). The prefix's finite
/// counterterms come from vakint, with the 4D lane's settings; a connected
/// host A integrates its plain U and adds its self-suffix.
fn production_add_back_and_image(
    ring: Ring,
    vacuum_subtraction: bool,
    prefix: RingPrefix,
) -> Result<ProductionAddBack> {
    let mut graph = ring.physical_graph()?;
    let settings = GenerationSettings {
        uv: UVgenerationSettings {
            generate_integrated: true,
            ..ring_uv_settings(ApproximationType::IR)
        },
        medium: thermal_medium(vacuum_subtraction),
        ..Default::default()
    };
    let vakint_settings =
        Wood::new(CutStructure::empty(&graph), &graph, &settings)?.vakint_settings;
    let regions = ring.expected_regions(&graph, vacuum_subtraction);
    let node = |edges: &[usize], dod| match regions.get(&edge_set(edges)) {
        Some(region) => host_node(&graph, edges, dod, region.clone()),
        None => forest_node(&graph, edges, ApproximationType::MUV, dod),
    };
    let first = node(&[2, 3], 2)?;
    let second = node(ring.second_loop_edges(), 2)?;
    let union = [[2, 3].as_slice(), ring.second_loop_edges()].concat();
    let host = node(&union, 4)?;
    let join = forest_node(&graph, &union, ApproximationType::MUV, 4)?;
    let integrated =
        |edges: &[usize]| vakint_integrated(&graph, &settings.uv, &vakint_settings, edges, 2);
    let (integrated_first, integrated_second) =
        (integrated(&[2, 3])?, integrated(ring.second_loop_edges())?);
    let finite = [&integrated_first, &integrated_second]
        .map(IntegratedCts::physical_finite_counterterm_atom);
    let (prefix_node, components, integrated) = match prefix {
        RingPrefix::First => (&first, vec![(&finite[0], &first)], integrated_first),
        RingPrefix::Second => (&second, vec![(&finite[1], &second)], integrated_second),
        RingPrefix::Join => (
            &join,
            vec![(&finite[0], &first), (&finite[1], &second)],
            IntegratedCts::factorized_product([&integrated_first, &integrated_second], 3)?,
        ),
    };
    let cff = ThermalCff::production(&mut graph, vacuum_subtraction)?;
    let orientation = cff.orientation(vacuum_subtraction);
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(&cutset, orientation);
    let image = Direct3dApproximation::new(localizer, &mut graph, &settings.uv).run_integrated(
        &components,
        prefix_node,
        &host,
        prefix_node,
        &host,
        prefix_node,
    )?;
    let builder = FinalIntegrandBuilder::new(localizer, &settings.uv);
    let no_local_terms = image.map(|_| Ok(Atom::Zero))?;
    let add_back = builder.build_direct(&mut graph, prefix_node, &no_local_terms, &integrated)?;
    let image = builder.build_direct(&mut graph, &host, &image, &IntegratedCts::root())?;
    let add_back = final_ring_atom(&graph, add_back)?;
    let image = final_ring_atom(&graph, image)?;
    Ok(ProductionAddBack {
        graph,
        finite: integrated.physical_finite_counterterm_atom(),
        add_back,
        image,
    })
}

/// The production add-backs of `prefixes` on the physical `ring` cancel the
/// static mode of the q-gluons with their zero-mode images
/// ([`production_add_back_and_image`]), with and without
/// `vacuum_subtraction` as given, at the media of the pins.
fn assert_production_add_backs_cancel_their_static_mode(
    ring: Ring,
    prefixes: &[(bool, RingPrefix)],
) -> Result<()> {
    test_initialise()?;
    for &(vacuum_subtraction, prefix) in prefixes {
        let ProductionAddBack {
            graph,
            add_back,
            image,
            ..
        } = production_add_back_and_image(ring, vacuum_subtraction, prefix)?;
        assert_final_static_mode_cancels(
            &graph,
            &add_back,
            &image,
            &format!("{ring:?}, vacuum subtraction {vacuum_subtraction}, {prefix:?}"),
        )?;
    }
    Ok(())
}

/// The final `add_back` of a ring and its zero-mode `image`
/// ([`final_ring_atom`]) cancel the static mode of the q-gluons
/// ([`assert_static_mode_cancels`]) at every point of [`LOOP_MOMENTA`] and
/// every medium of [`PIN_MEDIA`], with the vacuum mass [`UV_MASS`], the
/// static-theory mass [`IR_MASS`] and `mu_r^2 = 4`.
fn assert_final_static_mode_cancels(
    graph: &Graph,
    add_back: &Atom,
    image: &Atom,
    context: &str,
) -> Result<()> {
    // The couplings of each atom: their sum may cancel some of them.
    let mut parameters = [add_back, image]
        .into_iter()
        .flat_map(unit_couplings)
        .unique_by(|(coupling, _)| coupling.clone())
        .collect::<Vec<_>>();
    parameters.extend([
        (Atom::var(GS.m_uv_vacuum), UV_MASS),
        (Atom::var(GS.m_ir_vacuum), IR_MASS),
        (Atom::var(GS.mu_r_sq), 4.0),
    ]);
    for (index, [q, _, _]) in LOOP_MOMENTA.into_iter().enumerate() {
        for (beta, mu) in PIN_MEDIA {
            assert_static_mode_cancels(
                |s| {
                    let kinematics = scaled_ring(index, s);
                    let point = kinematics.point(mu, &parameters).with_beta(beta);
                    let [add_back, image] = [add_back, image]
                        .map(|atom| evaluate_points_complex(graph, atom, slice::from_ref(&point)));
                    let q4 = (norm3(q) * s).powi(4);
                    let scaled = |value: Complex<f64>| Complex::new(q4 * value.re, q4 * value.im);
                    Ok((scaled(add_back?[0]), scaled(image?[0])))
                },
                &format!("{context}, point {index}, beta {beta}, mu {mu}"),
            )?;
        }
    }
    Ok(())
}

/// On the tadpole ring, the add-back of the tadpole, `C_B ∝ M^2` from
/// vakint, held as `-(3/8)` of `mUVexp^2` and `mUV^2` terms whose logarithms
/// cancel only between the two, is cancelled at `q -> 0` by its image
/// `(∫B; S0_h|A)`, with and without vacuum subtraction. Before the zero-mode
/// chart co-scaled the vertex's expansion mass, the image kept only the
/// `mUV^2` terms, `-1.72` times the add-back at `M = 1` and `mu_r = 2`.
#[test]
fn tadpole_ring_production_add_back_cancels_its_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(
        Ring::Tadpole,
        &[(true, RingPrefix::Second), (false, RingPrefix::Second)],
    )
}

/// On the quark ring, the add-backs of the first bubble, the second bubble
/// and their join, from vakint, are cancelled at `q -> 0` by their images
/// `(∫A; S0_h|B)`, `(∫B; S0_h|A)` and `(∫A∫B; S0_h|∅)`.
#[test]
fn quark_ring_production_add_backs_cancel_their_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(
        Ring::Quark,
        &RingPrefix::ALL.map(|prefix| (true, prefix)),
    )
}

/// On the ghost ring, the add-back of host A with its self-suffix, those of
/// the ghost bubble and of the join, from vakint, are cancelled at `q -> 0`
/// by their images under the union host's zero mode, the self-suffix's
/// `(∫A; S0ᶜ_A, S0_h|B)` included.
#[test]
fn ghost_ring_production_add_backs_cancel_their_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(
        Ring::Ghost,
        &RingPrefix::ALL.map(|prefix| (true, prefix)),
    )
}

/// On the gluon ring, with vacuum subtraction, the add-back of host A with its
/// self-suffix, from vakint, is cancelled at `q -> 0` by its images under the
/// union host's zero mode, `(∫A; S0_h|B)` and `(∫A; S0ᶜ_A, S0_h|B)`.
#[test]
fn gluon_ring_production_add_back_of_host_a_cancels_its_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(Ring::Mixed, &[(true, RingPrefix::First)])
}

/// On the gluon ring, with vacuum subtraction, the add-back of the gluon
/// bubble, from vakint, is cancelled at `q -> 0` by its image `(∫B; S0_h|A)`.
#[test]
fn gluon_ring_production_add_back_of_the_gluon_bubble_cancels_its_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(Ring::Mixed, &[(true, RingPrefix::Second)])
}

/// On the gluon ring, with vacuum subtraction, the add-back `C_A C_B` of the
/// factorized join, from vakint, is cancelled at `q -> 0` by its image
/// `(∫A∫B; S0_h|∅)`. It sits on the bare gluon loop, whose three-gluon
/// vertices raise the numerator's energy degree above its propagators', so
/// its localization omits the polynomial energy contact, which vacuum
/// subtraction removes ([`gluon_ring_join_add_back_needs_vacuum_subtraction`]).
#[test]
fn gluon_ring_production_join_add_back_cancels_its_static_mode() -> Result<()> {
    assert_production_add_backs_cancel_their_static_mode(Ring::Mixed, &[(true, RingPrefix::Join)])
}

/// On the gluon ring, with vacuum subtraction, the self-suffix of host A's
/// add-back from vakint replayed through the static U of the gluon bubble at
/// the factorized join, `(∫A; S0ᶜ_A, U_B°)`, is cancelled by its image under
/// the union host's zero mode, as final integrands of the join and of the
/// union host. The sector is the static mode itself, polynomial in q within
/// the host's depth, so the image is minus the whole add-back at every point.
///
/// The static U collects the self-suffix's terms of one Taylor coefficient
/// into a sum, each with the q-class propagators, which the hard-vacuum
/// sector of the union host reads off as projected
/// ([`projected_summands_keep_their_order_through_sums`]). Before, it took
/// that sum for an unprojected frozen row and stopped with
/// `ZeroModeError::FrozenHardVacuum`. The sibling sector `(∫A; U_B)` of the
/// ordinary rows also cancels its static mode, but its add-back takes minutes
/// to assemble, beyond a unit test.
#[test]
fn gluon_ring_static_u_on_the_self_suffix_cancels_its_static_mode() -> Result<()> {
    test_initialise()?;
    let ring = Ring::Mixed;
    let vacuum_subtraction = true;
    let mut graph = ring.physical_graph()?;
    let settings = GenerationSettings {
        uv: UVgenerationSettings {
            generate_integrated: true,
            ..ring_uv_settings(ApproximationType::IR)
        },
        medium: thermal_medium(vacuum_subtraction),
        ..Default::default()
    };
    let vakint_settings =
        Wood::new(CutStructure::empty(&graph), &graph, &settings)?.vakint_settings;
    let mut regions = ring.expected_regions(&graph, vacuum_subtraction);
    let mut host = |edges: &[usize], dod| {
        let region = regions
            .remove(&edge_set(edges))
            .ok_or_else(|| eyre!("{edges:?} is a zero-mode host of the gluon ring"))?;
        host_node(&graph, edges, dod, region)
    };
    let union = [2, 3, 4, 5];
    let (first, union_host) = (host(&[2, 3], 2)?, host(&union, 4)?);
    let second = forest_node(&graph, &[4, 5], ApproximationType::MUV, 2)?;
    let join = forest_node(&graph, &union, ApproximationType::MUV, 4)?;
    let root = empty_forest_node(&graph);
    let finite = vakint_integrated(&graph, &settings.uv, &vakint_settings, &[2, 3], 2)?
        .physical_finite_counterterm_atom();
    let cff = ThermalCff::production(&mut graph, vacuum_subtraction)?;
    let orientation = cff.orientation(vacuum_subtraction);
    let cutset = CutSet::empty(graph.n_hedges());
    let localizer = Localizer::new(&cutset, orientation);
    let replay = Direct3dApproximation::new(localizer, &mut graph, &settings.uv).run_integrated(
        &[(&finite, &first)],
        &first,
        &second,
        &root,
        &second,
        &root,
    )?;
    let [_, self_suffix] = <[DirectSector; 2]>::try_from(replay.sectors()?.to_vec())
        .map_err(|_| eyre!("host A's add-back has one self-suffix"))?;
    assert!(
        !has_top_level_propagators(&self_suffix.active, &[(0, 4)]),
        "the static U carries the q-class propagators of the self-suffix into sums"
    );
    let add_back = Direct3dCts::from_sectors(vec![self_suffix])?;
    let image = Direct3dApproximation::new(localizer, &mut graph, &settings.uv).run_local(
        &add_back,
        &union_host,
        &join,
        &union_host,
        &join,
    )?;
    let builder = FinalIntegrandBuilder::new(localizer, &settings.uv);
    let add_back = builder.build_direct(&mut graph, &join, &add_back, &IntegratedCts::root())?;
    let image = builder.build_direct(&mut graph, &union_host, &image, &IntegratedCts::root())?;
    let (add_back, image) = (
        final_ring_atom(&graph, add_back)?,
        final_ring_atom(&graph, image)?,
    );
    let context = format!("{ring:?}, (∫A; S0ᶜ_A, U_B°)");
    assert_final_static_mode_cancels(&graph, &add_back, &image, &context)?;
    // The sector is the static mode itself, polynomial in q within the
    // union host's depth, so the image is minus the whole add-back.
    let mut parameters = unit_couplings(&add_back);
    parameters.extend([
        (Atom::var(GS.m_uv_vacuum), UV_MASS),
        (Atom::var(GS.m_ir_vacuum), IR_MASS),
        (Atom::var(GS.mu_r_sq), 4.0),
    ]);
    for (index, momenta) in LOOP_MOMENTA.into_iter().enumerate() {
        let kinematics = RingKinematics::new(momenta);
        let point = kinematics.point(QUARK_MU, &parameters);
        let [add_back, image] = [&add_back, &image]
            .map(|atom| evaluate_points_complex(&graph, atom, slice::from_ref(&point)));
        let (add_back, image) = (add_back?[0], image?[0]);
        assert!(
            (add_back.re + image.re).hypot(add_back.im + image.im)
                <= 1e-10 * add_back.re.hypot(add_back.im),
            "{context}, point {index}: the image is not minus the add-back, image / add-back = \
             {:?}",
            image / add_back
        );
    }
    Ok(())
}

/// The vacuum-subtracted Matsubara sum over the q-gluon energy of a physical
/// ring's 4D integrand on the bare gluon loop `{e0, e1}` at `kinematics` and
/// `beta`: `T sum_n F(i w_n) - int dw/(2 pi) F(i w)` of
/// `F(q0) = numerator(q0) / (q0^2 - |q|^2)^2`, with the contracted 4D
/// `numerator` at `q0` on the gluons, `-q0` on e3 = p - q and e5 = k - q, the
/// other energies zero, and the spatial components and further `parameters`
/// of the point. It is `-sum_z Res_z [F(z) h(z)]` over the double poles
/// `z = +-|q|`, with `h(z) = (coth(beta z / 2) -+ 1) / 2`, the Bose
/// occupation minus its vacuum value, each residue a trapezoidal contour
/// integral on a circle of radius `|q| / 2`. A polynomial part of `F` has no
/// residue. Independent of the CFF.
fn gluon_loop_matsubara_sum(
    graph: &Graph,
    numerator: &Atom,
    kinematics: &RingKinematics,
    beta: f64,
    parameters: &[(Atom, f64)],
) -> Result<Complex<f64>> {
    const NODES: usize = 64;
    let components = kinematics.components();
    let symbols = (0..6)
        .map(|edge| GS.emr_mom(EdgeIndex(edge), GS.cind(0)))
        .chain(components.iter().map(|(atom, _)| atom.clone()))
        .chain(parameters.iter().map(|(atom, _)| atom.clone()))
        .collect::<Vec<_>>();
    let mut evaluator = numerator
        .evaluator(&symbols)
        .function_map(graph.param_builder.fn_map.clone())
        .build()
        .map_err(|error| eyre!("{error}"))?
        .map_coeff(&|coefficient| Complex::new(coefficient.re.to_f64(), coefficient.im.to_f64()));
    let real = |value: f64| Complex::new(value, 0.0);
    let zero = real(0.0);
    let exp = |w: Complex<f64>| {
        let modulus = w.re.exp();
        Complex::new(modulus * w.im.cos(), modulus * w.im.sin())
    };
    let energy = kinematics.energies[0];
    let mut sum = zero;
    for pole in [energy, -energy] {
        for node in 0..NODES {
            let angle = 2.0 * PI * node as f64 / NODES as f64;
            let offset = Complex::new(0.5 * energy * angle.cos(), 0.5 * energy * angle.sin());
            let z = real(pole) + offset;
            let values = [z, z, zero, zero - z, zero, zero - z]
                .into_iter()
                .chain(components.iter().map(|(_, value)| real(*value)))
                .chain(parameters.iter().map(|(_, value)| real(*value)))
                .collect::<Vec<_>>();
            // coth(beta z / 2) = 1 + 2 / (exp(beta z) - 1).
            let coth = real(1.0) + real(2.0) / (exp(real(beta) * z) - real(1.0));
            let occupation = (coth - real(pole.signum())) / real(2.0);
            let propagators = (z * z - real(energy * energy)) * (z * z - real(energy * energy));
            sum += evaluator.evaluate_single(&values) * occupation * offset
                / (propagators * real(NODES as f64));
        }
    }
    Ok(zero - sum)
}

/// The production add-back of the factorized join, the product `C_A C_B`
/// of the bubbles' counterterms from vakint on the bare gluon loop
/// `{e0, e1}`, is the vacuum-subtracted Matsubara sum of its 4D integrand
/// over the gluons' energy ([`gluon_loop_matsubara_sum`]), with one
/// normalization on the quark, ghost and gluon rings, at every point and
/// medium of the pins. On the gluon ring the three-gluon vertices raise the
/// numerator's energy degree one above the two propagators', so the
/// localization omits the polynomial energy contact
/// ([`gluon_ring_join_add_back_needs_vacuum_subtraction`]); on the quark and
/// ghost rings it localizes without one. The add-back is a value the image
/// pins cannot see ([`production_add_back_and_image`]): an error that the
/// image repeats cancels there.
///
/// The normalization is `-(2 pi)^-3`: the gluon loop's spatial measure, and
/// its one energy integration's `i` against the `i^3` that the evaluation
/// divides out for the ring's three loops ([`evaluate_points_complex`]);
/// the integrated bubbles' loops carry the localization kernel, here one.
#[test]
fn join_add_backs_are_the_matsubara_sums_of_the_gluon_loop() -> Result<()> {
    test_initialise()?;
    let normalization = -1.0 / (2.0 * PI).powi(3);
    for ring in [Ring::Quark, Ring::Ghost, Ring::Mixed] {
        let ProductionAddBack {
            graph,
            finite,
            add_back,
            ..
        } = production_add_back_and_image(ring, true, RingPrefix::Join)?;
        let union = graph.edge_subgraph(&[[2, 3].as_slice(), ring.second_loop_edges()].concat());
        let cograph = graph.full_filter().subtract(&union);
        let numerator = normalize_local_operation(
            &graph,
            &contract_numerators(
                &(finite * physical_numerator(&graph, &cograph, &union) * graph.global_atom())
                    .replace(GS.dim_epsilon)
                    .with(Atom::Zero),
            )?,
        );
        let mut parameters = [&add_back, &numerator]
            .into_iter()
            .flat_map(unit_couplings)
            .unique_by(|(coupling, _)| coupling.clone())
            .collect::<Vec<_>>();
        parameters.extend([
            (Atom::var(GS.m_uv_vacuum), UV_MASS),
            (Atom::var(GS.m_ir_vacuum), IR_MASS),
            (Atom::var(GS.mu_r_sq), 4.0),
        ]);
        for (index, momenta) in LOOP_MOMENTA.into_iter().enumerate() {
            let kinematics = RingKinematics::new(momenta);
            for (beta, mu) in PIN_MEDIA {
                let point = kinematics.point(mu, &parameters).with_beta(beta);
                let actual =
                    evaluate_points_complex(&graph, &add_back, slice::from_ref(&point))?[0];
                let reference =
                    gluon_loop_matsubara_sum(&graph, &numerator, &kinematics, beta, &parameters)?;
                assert_close_complex(
                    actual,
                    Complex::new(normalization * reference.re, normalization * reference.im),
                    1e-10,
                    &format!("{ring:?}, point {index}, beta {beta}, mu {mu}"),
                );
            }
        }
    }
    Ok(())
}

/// On the physical gluon ring the add-back of the factorized join, the
/// product `C_A C_B` of the bubbles' counterterms from vakint, sits on the
/// bare gluon loop `{e0, e1}`. Its three-gluon vertices raise the numerator's
/// energy degree one above the two propagators', so contracting the last of
/// them leaves a polynomial energy contact of power one. With vacuum
/// subtraction the prefix localizer omits it and keeps a nonzero residue
/// part; without, it stops with the dedicated error, which names the join
/// and its cograph.
#[test]
fn gluon_ring_join_add_back_needs_vacuum_subtraction() -> Result<()> {
    test_initialise()?;
    let ring = Ring::Mixed;
    let mut graph = ring.physical_graph()?;
    let settings = GenerationSettings {
        uv: UVgenerationSettings {
            generate_integrated: true,
            ..ring_uv_settings(ApproximationType::IR)
        },
        medium: thermal_medium(true),
        ..Default::default()
    };
    let vakint_settings =
        Wood::new(CutStructure::empty(&graph), &graph, &settings)?.vakint_settings;
    let first = forest_node(&graph, &[2, 3], ApproximationType::MUV, 2)?;
    let second = forest_node(&graph, ring.second_loop_edges(), ApproximationType::MUV, 2)?;
    let union = [[2, 3].as_slice(), ring.second_loop_edges()].concat();
    let join = forest_node(&graph, &union, ApproximationType::MUV, 4)?;
    let [finite_first, finite_second] =
        [[2, 3].as_slice(), ring.second_loop_edges()].map(|edges| {
            vakint_integrated(&graph, &settings.uv, &vakint_settings, edges, 2)
                .map(|integrated| integrated.physical_finite_counterterm_atom())
        });
    let components = [(&finite_first?, &first), (&finite_second?, &second)];
    let cutset = CutSet::empty(graph.n_hedges());
    let cff = ThermalCff::production(&mut graph, true)?;
    for vacuum_subtraction in [true, false] {
        let localizer = Localizer::new(&cutset, cff.orientation(vacuum_subtraction));
        let localized = Direct3dApproximation::new(localizer, &mut graph, &settings.uv)
            .localize_prefix(&components, &join);
        if vacuum_subtraction {
            let [sector] = <[DirectSector; 1]>::try_from(localized?)
                .map_err(|_| eyre!("the join's prefix has one sector"))?;
            assert!(!sector.active.is_zero(), "the residue part is kept");
            continue;
        }
        let error = localized.expect_err("the join's add-back needs vacuum subtraction");
        match error.downcast_ref::<IntegratedLocalizationError>() {
            Some(IntegratedLocalizationError::EnergyContactRequiresVacuumSubtraction {
                node,
                cograph,
                power,
            }) => {
                assert_eq!(node, &edge_set(&union));
                assert_eq!(cograph, &edge_set(&[0, 1]));
                assert_eq!(*power, 1);
            }
            None => panic!("expected the dedicated localization error, got {error:?}"),
        }
    }
    Ok(())
}

/// The self-suffix of host A on the physical ghost ring, the second sector
/// that the prefix localizer adds to the host's add-back, against references
/// that the zero-mode chart does not build, with vakint's counterterm `C_A`
/// and both sides through the final assembly of the host's node:
/// - it is minus the exact static mode `[x_q^-4 x_4^-2]` of the plain
///   add-back, re-keyed to static energies, up to `O(q^2)` as `q -> 0`,
///   which fixes its presence, sign and multiplicity;
/// - at every point it is `T^2 C_A(0) N_0(k, q) / (q^4 k^2 (k - q)^2)` up to
///   one constant, with `C_A` at zero external momentum, holding the vacuum
///   mass as `mUV` and as `mUVexp`, and the cograph's numerator `N_0` at zero
///   static energies, the same at every point, temperature and chemical
///   potential, with and without vacuum subtraction.
///
/// The image pins of the add-backs ([`production_add_back_and_image`]) cannot
/// see an error of the self-suffix, which the image of the add-back repeats.
#[test]
fn ghost_ring_production_self_suffix_is_the_static_vertex() -> Result<()> {
    test_initialise()?;
    let ring = Ring::Ghost;
    let mut constants = Vec::new();
    for vacuum_subtraction in [true, false] {
        let mut graph = ring.physical_graph()?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: true,
                ..ring_uv_settings(ApproximationType::IR)
            },
            medium: thermal_medium(vacuum_subtraction),
            ..Default::default()
        };
        let vakint_settings =
            Wood::new(CutStructure::empty(&graph), &graph, &settings)?.vakint_settings;
        let region = ring
            .expected_regions(&graph, vacuum_subtraction)
            .remove(&edge_set(&[2, 3]))
            .ok_or_else(|| eyre!("A is a connected host of the ghost ring"))?;
        let host = host_node(&graph, &[2, 3], 2, region)?;
        let finite = vakint_integrated(&graph, &settings.uv, &vakint_settings, &[2, 3], 2)?
            .physical_finite_counterterm_atom();
        // `C_A(0) N_0`: the cograph's numerator holds the ghost momenta k and
        // k - q, so the momentum q of the contracted vertex is its own.
        let static_vertex = normalize_local_operation(
            &graph,
            &contract_numerators(
                &(&finite * HostANumerators::new(&graph).vertices)
                    .replace(GS.dim_epsilon)
                    .with(Atom::Zero),
            )?
            .replace(function!(GS.emr_mom, 0, W_.x___))
            .with(Atom::Zero),
        );
        let cff = ThermalCff::production(&mut graph, vacuum_subtraction)?;
        let orientation = cff.orientation(vacuum_subtraction);
        let cutset = CutSet::empty(graph.n_hedges());
        let localizer = Localizer::new(&cutset, orientation);
        let [plain, self_suffix] = <[DirectSector; 2]>::try_from(
            Direct3dApproximation::new(localizer, &mut graph, &settings.uv)
                .localize_prefix(&[(&finite, &host)], &host)?,
        )
        .map_err(|_| eyre!("the add-back of host A has one self-suffix"))?;
        let static_mode = DirectSector {
            active: project(
                &graph,
                &plain.active,
                vec![class(&graph, &[0, 1], &[4])?, class(&graph, &[4], &[2])?],
                ZeroModeSector::Complete,
            )?
            .remove(&vec![4, 2])
            .ok_or_else(|| eyre!("the static order is projected"))?
            .rekey_static(orientation, &static_lines())?,
            ..self_suffix.clone()
        };
        let builder = FinalIntegrandBuilder::new(localizer, &settings.uv);
        let mut assemble = |sector: DirectSector| -> Result<Atom> {
            let integrands = builder.build_direct(
                &mut graph,
                &host,
                &Direct3dCts::from_sectors(vec![sector])?,
                &IntegratedCts::root(),
            )?;
            final_ring_atom(&graph, integrands)
        };
        let self_suffix = assemble(self_suffix)?;
        let static_mode = assemble(static_mode)?;
        let mut parameters = unit_couplings(&(&self_suffix + &static_mode));
        parameters.extend([
            (Atom::var(GS.m_uv_vacuum), UV_MASS),
            (Atom::var(GS.mu_r_sq), 4.0),
        ]);
        let value = |atom: &Atom, kinematics: &RingKinematics, beta: f64, mu: f64| {
            let point = kinematics.point(mu, &parameters).with_beta(beta);
            Ok::<_, eyre::Report>(
                evaluate_points_complex(&graph, atom, slice::from_ref(&point))?[0],
            )
        };
        for (index, momenta) in LOOP_MOMENTA.into_iter().enumerate() {
            for (beta, mu) in PIN_MEDIA {
                let context = format!(
                    "vacuum subtraction {vacuum_subtraction}, point {index}, beta {beta}, mu {mu}"
                );
                for s in ADD_BACK_SCALINGS {
                    let kinematics = scaled_ring(index, s);
                    let actual = value(&self_suffix, &kinematics, beta, mu)?;
                    let reference = value(&static_mode, &kinematics, beta, mu)?;
                    let scale = reference.re.hypot(reference.im);
                    assert!(scale > 0.0, "{context}, s {s}: no static mode");
                    assert!(
                        (actual.re + reference.re).hypot(actual.im + reference.im)
                            <= 100.0 * s * s * scale,
                        "{context}, s {s}: the self-suffix is not minus the static mode of the \
                         add-back, self-suffix / static mode = {:?}",
                        actual / reference
                    );
                }
                let kinematics = RingKinematics::new(momenta);
                let [e0, e1, _, _, e4, e5] = kinematics.energies;
                let propagators = (beta * e0 * e1 * e4 * e5).powi(2);
                let vertex = kinematics.static_numerator(&graph, &static_vertex, &parameters)?;
                let actual = value(&self_suffix, &kinematics, beta, mu)?;
                constants.push((context, actual * Complex::new(propagators, 0.0) / vertex));
            }
        }
    }
    let (_, first) = constants[0].clone();
    assert!(first.re.hypot(first.im) > 0.0, "no self-suffix");
    for (context, constant) in constants {
        assert_close_complex(
            constant,
            first,
            1e-9,
            &format!("{context}: the self-suffix over T^2 C_A(0) N_0 / (q^4 k^2 (k - q)^2)"),
        );
    }
    Ok(())
}
