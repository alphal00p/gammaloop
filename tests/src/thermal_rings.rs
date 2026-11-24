//! Generated structure, profiles, approach probes and the prescription check
//! shared by the thermal two-bubble ring tests.
//!
//! Each ring joins a d-quark polarization to a second polarization loop by two
//! gluons carrying q, with loop momenta [q, p, k] = [e0, e2, e4].

use std::fmt::{Display, Formatter, Result as FmtResult};
use std::path::PathBuf;

use color_eyre::Result;
use eyre::{Context, eyre};
use gammaloop_api::commands::Profile;
use gammaloop_api::commands::inspect::Inspect;
use gammaloop_api::commands::profile::{InfraRedProfile, UltraVioletProfile};
use gammaloop_api::state::ProcessRef;
use gammalooprs::integrands::process::ir::SingleLimitReport;
use gammalooprs::momentum::ThreeMomentum;
use gammalooprs::processes::ProcessCollection;
use gammalooprs::settings::RuntimeSettings;
use gammalooprs::settings::runtime::{Precision, StabilityLevelSetting};
use gammalooprs::uv::profile::{UV_PROFILE_MAX_DOD, UVLimitSelection, UVProfilePassFail};
use gammalooprs::uv::uv_graph::UVE;
use linnet::half_edge::involution::EdgeIndex;
use spenso::algebra::complex::Complex;

use crate::{CLIState, clean_test, get_test_cli, get_tests_workspace_path};

/// The setting that the hostless error of a freezing finite-temperature ring
/// names: only the IR prescription lets the zero-mode policy host its static
/// soft set.
pub const IR_PRESCRIPTION_SETTING: &str =
    "renormalization_prescription.massless_power_divergent = IR";

/// Distance of a soft or joint soft profile from its expected integer power,
/// with the profiler's 3D soft volume factor included.
pub const SCALING_TOLERANCE: f64 = 0.25;
/// The power-law fit quality required of a soft or joint soft profile.
pub const MIN_R_SQUARED: f64 = 0.98;

/// The seed of every profile, so that separate generations of a ring sample
/// identical rays.
const PROFILE_SEED: u64 = 1337;
/// IR profiles scale the soft momenta by 10^-2 .. 10^-5.
const IR_MIN_SCALE_EXPONENT: f64 = -2.0;
const IR_MAX_SCALE_EXPONENT: f64 = -5.0;
/// Points of every UV profile.
const UV_PROFILE_POINTS: usize = 8;

/// Whether a ring generation subtracts the vacuum part of its distributions.
/// Both settings are supported at finite temperature.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VacuumSubtraction {
    On,
    Off,
}

impl VacuumSubtraction {
    pub fn enabled(self) -> bool {
        matches!(self, Self::On)
    }

    /// The suffix that keeps the states of the two settings apart.
    pub fn state_suffix(self) -> &'static str {
        match self {
            Self::On => "vacuum_subtracted",
            Self::Off => "unsubtracted_vacuum",
        }
    }

    /// The bare ring: no UV or soft subtraction.
    pub fn bare_commands(self) -> &'static [&'static str] {
        match self {
            Self::On => &["run generate_bare"],
            Self::Off => &[
                "set global kv global.generation.medium.vacuum_subtraction=false",
                "run generate_bare",
            ],
        }
    }

    /// The zero-mode subtraction of the run card's IR prescription with
    /// `counterterms`.
    pub fn subtraction_commands(self, counterterms: Counterterms) -> &'static [&'static str] {
        match (counterterms, self) {
            (Counterterms::Local, Self::On) => &["run generate"],
            (Counterterms::Local, Self::Off) => &["run generate_novacsub"],
            (Counterterms::Integrated, Self::On) => &["run generate_integrated"],
            (Counterterms::Integrated, Self::Off) => &["run generate_integrated_novacsub"],
        }
    }
}

impl Display for VacuumSubtraction {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            Self::On => write!(f, "with vacuum subtraction"),
            Self::Off => write!(f, "without vacuum subtraction"),
        }
    }
}

/// The UV counterterms of a ring generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Counterterms {
    /// Local counterterms only: `uv.generate_integrated = false`.
    Local,
    /// Every counterterm with `general.m_uv` also adds back its integrated
    /// finite conversion to MS-bar, inside the zero-mode terms too, localized
    /// at the thermal scale. The zero-mode operations and the static theory's
    /// counterterms, with `general.m_ir`, stay local.
    Integrated,
}

impl Counterterms {
    /// Whether the generation adds back integrated counterterms.
    pub fn integrated(self) -> bool {
        matches!(self, Self::Integrated)
    }
}

impl Display for Counterterms {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        match self {
            Self::Local => write!(f, "with local counterterms"),
            Self::Integrated => write!(f, "with integrated counterterms"),
        }
    }
}

/// One thermal ring: its process and integrand name, its second polarization
/// loop and its superficially divergent UV limits in the generation LMB.
pub struct ThermalRing {
    pub name: &'static str,
    /// PDG code of the second polarization loop's edges.
    pub second_loop_pdg: isize,
    pub second_loop_edges: &'static [usize],
    pub uv_divergent_limits: usize,
}

pub const QUARK_RING: ThermalRing = ThermalRing {
    name: "thermal_quark_ring",
    second_loop_pdg: 1,
    second_loop_edges: &[4, 5],
    uv_divergent_limits: 8,
};

pub const MIXED_RING: ThermalRing = ThermalRing {
    name: "thermal_gluon_quark_ring",
    second_loop_pdg: 21,
    second_loop_edges: &[4, 5],
    uv_divergent_limits: 8,
};

pub const GHOST_RING: ThermalRing = ThermalRing {
    name: "thermal_ghost_quark_ring",
    second_loop_pdg: 9000005,
    second_loop_edges: &[4, 5],
    uv_divergent_limits: 8,
};

/// The gluon tadpole e4 carries k through no other edge.
pub const TADPOLE_RING: ThermalRing = ThermalRing {
    name: "thermal_tadpole_quark_ring",
    second_loop_pdg: 21,
    second_loop_edges: &[4],
    uv_divergent_limits: 7,
};

impl ThermalRing {
    /// [`fresh_soft_ir_cli`] with the ring's run card.
    pub fn fresh_cli(&self, state_name: &str) -> Result<(CLIState, PathBuf)> {
        fresh_soft_ir_cli(self.name, state_name)
    }
}

/// A fresh CLI state at `<tests workspace>/<state_name>`, loaded with the run
/// card `<name>_soft_ir.toml`; nothing is generated yet. Returns the state and
/// its path, which the caller cleans once the state is no longer needed.
pub fn fresh_soft_ir_cli(name: &str, state_name: &str) -> Result<(CLIState, PathBuf)> {
    let state_path = get_tests_workspace_path().join(state_name);
    let cli = get_test_cli(
        Some(format!("{name}_soft_ir.toml").into()),
        &state_path,
        None,
        true,
    )?;
    Ok((cli, state_path))
}

/// Check that `ring` generated as one vacuum graph that keeps the loop
/// momenta [q, p, k] = [e0, e2, e4], so S(e0) scales both gluon carriers,
/// with the same q on both gluons and the requested second loop.
pub fn assert_ring_structure(cli: &CLIState, ring: &ThermalRing) -> Result<()> {
    let name = ring.name;
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_owned())))?;
    let ProcessCollection::Amplitudes(amplitudes) =
        &cli.state.process_list.processes[process_id].collection
    else {
        panic!("{name} must be a vacuum amplitude");
    };
    let amplitude = amplitudes.get(name).expect("the ring was generated");
    assert_eq!(amplitude.graphs.len(), 1);
    let graph = &amplitude.graphs[0].graph;
    let lmb = &graph.loop_momentum_basis;
    assert_eq!(
        lmb.loop_edges
            .iter()
            .copied()
            .map(usize::from)
            .collect::<Vec<_>>(),
        [0, 2, 4],
        "{name} must retain [q, p, k] so S(e0) scales both gluon carriers"
    );
    assert_eq!(
        lmb.edge_signatures[EdgeIndex::from(0)],
        lmb.edge_signatures[EdgeIndex::from(1)],
        "the two ring gluons of {name} must carry the same q"
    );
    for &edge in ring.second_loop_edges {
        assert_eq!(
            graph[EdgeIndex::from(edge)]
                .particle_pdg_code()
                .expect("ring edges retain their model particles"),
            ring.second_loop_pdg,
            "{name} must contain the requested second polarization loop"
        );
    }
    Ok(())
}

/// Profile the IR limits `limits` of the integrand `name` of the process
/// `name` in one call, on the scaling range 1e-2 .. 1e-5, and return their
/// reports in the order of `limits`. The sampled rays depend on the whole
/// select string and on `n_points`, so two generations that profile the same
/// limits with the same `n_points` evaluate identical rays.
pub fn profile_ir_limits(
    cli: &mut CLIState,
    name: &str,
    limits: &[&str],
    n_points: usize,
) -> Result<Vec<SingleLimitReport>> {
    let mut report = Profile::InfraRed(InfraRedProfile {
        process: Some(ProcessRef::Unqualified(name.to_owned())),
        integrand_name: Some(name.to_owned()),
        n_points,
        min_scale_exponent: IR_MIN_SCALE_EXPONENT,
        max_scale_exponent: IR_MAX_SCALE_EXPONENT,
        seed: Some(PROFILE_SEED),
        select: Some(format!("{name} {}", limits.join(" "))),
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?
    .unwrap_ir();
    assert_eq!(report.results_per_graph.len(), 1, "{report}");
    let mut graph = report
        .results_per_graph
        .pop()
        .expect("the ring profiles exactly one graph, as asserted above");
    assert_eq!(graph.graph_name, name);
    assert_eq!(graph.single_limit_reports.len(), limits.len(), "{graph}");
    let mut profiles = Vec::with_capacity(limits.len());
    for limit in limits {
        let index = graph
            .single_limit_reports
            .iter()
            .position(|profile| profile.limit_name == *limit)
            .unwrap_or_else(|| panic!("{name} did not profile {limit}: {graph}"));
        let profile = graph.single_limit_reports.swap_remove(index);
        assert!(profile.orientation_label.is_none());
        profiles.push(profile);
    }
    Ok(profiles)
}

/// [`profile_ir_limits`] of the single limit `limit`.
pub fn profile_ir_limit(
    cli: &mut CLIState,
    name: &str,
    limit: &str,
    n_points: usize,
) -> Result<SingleLimitReport> {
    let mut profiles = profile_ir_limits(cli, name, &[limit], n_points)?;
    Ok(profiles
        .pop()
        .expect("one report per profiled limit, as asserted by profile_ir_limits"))
}

/// One UV profile of a ring: which limits, on which scaling range. The
/// number of points, f128 and the seed are fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UvScan {
    pub selection: UVLimitSelection,
    pub min_scale_exponent: f64,
    pub max_scale_exponent: f64,
}

impl UvScan {
    /// The superficially divergent limits in the generation LMB, scaled by 1e3 .. 1e6.
    pub const ONLY_DIVERGENT: Self = Self {
        selection: UVLimitSelection::OnlyDivergent,
        min_scale_exponent: 3.0,
        max_scale_exponent: 6.0,
    };
}

/// Profile the UV limits `scan` selects of the integrand `name` of the
/// process `name`. The vacuum-subtracted all-hard limit decays exponentially,
/// so a vanishing fit is allowed. Every limit must decay at least as fast as
/// [`UV_PROFILE_MAX_DOD`] for the verdict to pass.
pub fn uv_verdict(cli: &mut CLIState, name: &str, scan: UvScan) -> Result<UVProfilePassFail> {
    let report = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(name.to_owned())),
        integrand_name: Some(name.to_owned()),
        selected_limits: scan.selection,
        n_points: UV_PROFILE_POINTS,
        min_scale_exponent: scan.min_scale_exponent,
        max_scale_exponent: scan.max_scale_exponent,
        seed: Some(PROFILE_SEED),
        use_f128: true,
        allow_vanishing_missing_fits: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?
    .unwrap_uv();
    assert!(!report.stopped_early);
    assert_eq!(report.graphs.len(), 1);
    assert_eq!(report.graphs[0].graph_name, name);
    for table in report.tables_per_graph(UV_PROFILE_MAX_DOD) {
        println!("{name} UV {:?} {scan:?}:\n{table}", scan.selection);
    }
    Ok(report.pass_fail(UV_PROFILE_MAX_DOD))
}

/// Check that the `generate_muv` block of `ring`'s run card fails in a fresh
/// state: under the MUV prescription no zero-mode host subtracts the static
/// soft set of the q-gluons, so the scheme validation rejects the generation
/// with the hostless error, which names [`IR_PRESCRIPTION_SETTING`].
pub fn assert_muv_prescription_is_hostless(ring: &ThermalRing, state_name: &str) -> Result<()> {
    let name = ring.name;
    let (mut cli, state_path) = ring.fresh_cli(state_name)?;
    let error = cli
        .run_command("run generate_muv")
        .expect_err("the MUV prescription must leave the ring without a zero-mode host");
    let message = format!("{error:?}");
    assert!(
        message.contains("no zero-mode host") && message.contains(IR_PRESCRIPTION_SETTING),
        "{name} under MUV must report the missing zero-mode host: {message}"
    );
    clean_test(state_path);
    Ok(())
}

/// Spatial loop momenta [q, p, k] = [e0, e2, e4] of one ring evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RingPoint {
    pub q: ThreeMomentum<f64>,
    pub p: ThreeMomentum<f64>,
    pub k: ThreeMomentum<f64>,
}

impl RingPoint {
    pub fn flattened(&self) -> Vec<f64> {
        [self.q, self.p, self.k]
            .into_iter()
            .flat_map(<[f64; 3]>::from)
            .collect()
    }
}

/// The base point of the approach probes, shared by every ring: loop momenta
/// (Q, P, K), a direction D with |D| = 1, and the unit vectors the probes
/// derive from them.
struct BasePoint {
    q: ThreeMomentum<f64>,
    p: ThreeMomentum<f64>,
    k: ThreeMomentum<f64>,
    direction: ThreeMomentum<f64>,
    q_hat: ThreeMomentum<f64>,
    /// The unit part of D orthogonal to Q.
    normal: ThreeMomentum<f64>,
}

impl BasePoint {
    fn new() -> Self {
        let q = ThreeMomentum::from([0.31, -0.42, 0.25]);
        let direction = ThreeMomentum::from([0.48, -0.36, 0.8]);
        Self {
            q,
            p: ThreeMomentum::from([0.7, 0.1, -0.33]),
            k: ThreeMomentum::from([-0.21, 0.45, 0.38]),
            direction,
            q_hat: unit(q),
            normal: unit(perpendicular(direction, q)),
        }
    }

    fn point(&self) -> RingPoint {
        RingPoint {
            q: self.q,
            p: self.p,
            k: self.k,
        }
    }
}

/// The approach probes (a)-(j) of a ring's zero-mode subtraction: the loci
/// where its counterterms could lose integrability or precision. Each maps a
/// parameter t to a [`RingPoint`] near the base point; t -> 0 approaches the
/// probed locus or limit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Probe {
    /// (a) p̂·q = |k| - |k - q| + t at fixed |p|, where the energy
    /// differences of the two bubbles cancel at linear order.
    CrossBubbleThreshold,
    /// (b) |p| = |p - q|: the quark bubble's two energies coincide.
    EqualQuarkEnergies,
    /// (c) k = p + t D.
    CoincidentLoopMomenta,
    /// (d) p = t D -> 0 at the small |q| = 1e-2 |Q|.
    SoftQuarkAtSoftQ,
    /// (e) q = 1.5 P + t n with n ⊥ P: q ∥ p with |q| > |p|, the collinear
    /// threshold of the quark bubble.
    CollinearQuarkThreshold,
    /// (f) k = t D -> 0 at fixed q, where a static bubble's 1/k³ pair cancels.
    SoftSecondLoop,
    /// (g) |k| = |k - q|: the second bubble's two energies coincide.
    EqualSecondLoopEnergies,
    /// (h) k = K⊥ + t q̂: k ⊥ q, the Landau-pair locus of the second bubble.
    SecondLoopOrthogonalToQ,
    /// (i) k = D / t -> ∞ at fixed q.
    HardSecondLoop,
    /// (j) p = P⊥ + t q̂: p ⊥ q, the Landau-pair locus of the quark bubble.
    QuarkOrthogonalToQ,
}

impl Probe {
    /// The probe's label, (a) to (j).
    pub fn label(self) -> &'static str {
        match self {
            Self::CrossBubbleThreshold => "(a)",
            Self::EqualQuarkEnergies => "(b)",
            Self::CoincidentLoopMomenta => "(c)",
            Self::SoftQuarkAtSoftQ => "(d)",
            Self::CollinearQuarkThreshold => "(e)",
            Self::SoftSecondLoop => "(f)",
            Self::EqualSecondLoopEnergies => "(g)",
            Self::SecondLoopOrthogonalToQ => "(h)",
            Self::HardSecondLoop => "(i)",
            Self::QuarkOrthogonalToQ => "(j)",
        }
    }

    /// The probed parameters, ordered towards the locus or limit.
    pub fn parameters(self) -> &'static [f64] {
        match self {
            Self::SoftQuarkAtSoftQ => &[1e-2, 1e-4, 1e-6],
            Self::SoftSecondLoop | Self::HardSecondLoop => &[1e-1, 1e-3, 1e-5],
            Self::CrossBubbleThreshold
            | Self::EqualQuarkEnergies
            | Self::CoincidentLoopMomenta
            | Self::CollinearQuarkThreshold
            | Self::EqualSecondLoopEnergies
            | Self::SecondLoopOrthogonalToQ
            | Self::QuarkOrthogonalToQ => &[1e-2, 1e-5, 1e-8],
        }
    }

    /// The magnitude of the momentum that `t` scales: |p| for (d), |k| for
    /// (f) and (i), and t itself on a locus.
    pub fn momentum_scale(self, t: f64) -> f64 {
        match self {
            Self::HardSecondLoop => 1.0 / t,
            _ => t,
        }
    }

    pub fn point(self, t: f64) -> RingPoint {
        let base = BasePoint::new();
        let point = base.point();
        let BasePoint {
            q,
            p,
            k,
            direction,
            q_hat,
            normal,
        } = base;
        match self {
            Self::CrossBubbleThreshold => {
                let cosine = (norm(k) - norm(k - q) + t) / norm(q);
                let sine = (1.0 - cosine * cosine).sqrt();
                RingPoint {
                    p: (q_hat * cosine + normal * sine) * norm(p),
                    ..point
                }
            }
            Self::EqualQuarkEnergies => RingPoint {
                p: q * 0.5 + normal * 0.6 + q_hat * t,
                ..point
            },
            Self::CoincidentLoopMomenta => RingPoint {
                k: p + direction * t,
                ..point
            },
            Self::SoftQuarkAtSoftQ => RingPoint {
                q: q * 1e-2,
                p: direction * t,
                ..point
            },
            Self::CollinearQuarkThreshold => RingPoint {
                q: p * 1.5 + unit(perpendicular(direction, p)) * t,
                ..point
            },
            Self::SoftSecondLoop => RingPoint {
                k: direction * t,
                ..point
            },
            Self::EqualSecondLoopEnergies => RingPoint {
                k: q * 0.5 + normal * 0.6 + q_hat * t,
                ..point
            },
            Self::SecondLoopOrthogonalToQ => RingPoint {
                k: perpendicular(k, q) + q_hat * t,
                ..point
            },
            Self::HardSecondLoop => RingPoint {
                k: direction * (1.0 / t),
                ..point
            },
            Self::QuarkOrthogonalToQ => RingPoint {
                p: perpendicular(p, q) + q_hat * t,
                ..point
            },
        }
    }
}

/// Rays from the base point on which the loop momenta are widely separated,
/// where the production stability ladder is checked against arb.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeparatedScalesRay {
    /// q = Q, p = P̂ / t, k = t K̂: a hard quark with a soft second loop.
    HardQuarkSoftSecondLoop,
    /// q = Q, p = t P̂, k = K: a quark of momentum t.
    HardQuark,
    /// q = t D, p = P, k = K̂ / t: a soft q with a hard second loop.
    SoftQHardSecondLoop,
    /// q = Q, p = 100 P̂, k = t K̂: a soft second loop next to a hard quark
    /// at a fixed |q| / Σ|k| of about 6e-3.
    SoftSecondLoopAtHardQuark,
}

impl SeparatedScalesRay {
    pub fn point(self, t: f64) -> RingPoint {
        let base = BasePoint::new();
        let point = base.point();
        match self {
            Self::HardQuarkSoftSecondLoop => RingPoint {
                p: unit(base.p) * (1.0 / t),
                k: unit(base.k) * t,
                ..point
            },
            Self::HardQuark => RingPoint {
                p: unit(base.p) * t,
                ..point
            },
            Self::SoftQHardSecondLoop => RingPoint {
                q: base.direction * t,
                k: unit(base.k) * (1.0 / t),
                ..point
            },
            Self::SoftSecondLoopAtHardQuark => RingPoint {
                p: unit(base.p) * 100.0,
                k: unit(base.k) * t,
                ..point
            },
        }
    }
}

fn norm(vector: ThreeMomentum<f64>) -> f64 {
    (vector * vector).sqrt()
}

fn unit(vector: ThreeMomentum<f64>) -> ThreeMomentum<f64> {
    vector * (1.0 / norm(vector))
}

/// The part of `vector` orthogonal to `axis`.
fn perpendicular(vector: ThreeMomentum<f64>, axis: ThreeMomentum<f64>) -> ThreeMomentum<f64> {
    let axis = unit(axis);
    vector + axis * -(vector * axis)
}

/// How a ring point is evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Evaluation {
    /// Arbitrary precision: the reference.
    Arb,
    /// One stability level in this precision, which accepts whatever it
    /// computes: the raw rounding of that precision.
    Single(Precision),
    /// The integrand's own stability levels, as an integration uses them.
    Production,
}

/// Evaluate the integrand `name` of the process `name` at the flattened
/// spatial loop momenta `point`.
pub fn evaluate_momentum_point(
    cli: &mut CLIState,
    name: &str,
    evaluation: Evaluation,
    point: &[f64],
) -> Result<Complex<f64>> {
    evaluate_as(cli, name, evaluation, |cli, precision| {
        inspect(cli, name, point, precision)
    })
}

/// Evaluate the integrand `name` of the process `name` at each of `points`.
pub fn evaluate_ring_points(
    cli: &mut CLIState,
    name: &str,
    evaluation: Evaluation,
    points: &[RingPoint],
) -> Result<Vec<Complex<f64>>> {
    evaluate_as(cli, name, evaluation, |cli, precision| {
        points
            .iter()
            .map(|point| inspect(cli, name, &point.flattened(), precision))
            .collect()
    })
}

/// What one `inspect` evaluates with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum InspectPrecision {
    /// Arbitrary precision, whatever the stability levels.
    Arb,
    /// The integrand's current stability levels.
    StabilityLevels,
}

fn inspect(
    cli: &mut CLIState,
    name: &str,
    point: &[f64],
    precision: InspectPrecision,
) -> Result<Complex<f64>> {
    let (_, value) = Inspect {
        process: Some(ProcessRef::Unqualified(name.to_owned())),
        integrand_name: Some(name.to_owned()),
        point: point.to_vec(),
        momentum_space: true,
        use_arb_prec: precision == InspectPrecision::Arb,
        ..Default::default()
    }
    .run(&mut cli.state)?;
    Ok(value)
}

/// Run `evaluate` as `evaluation` asks. A [`Evaluation::Single`] evaluation
/// replaces the integrand's stability levels for its duration and restores
/// them afterwards.
fn evaluate_as<T>(
    cli: &mut CLIState,
    name: &str,
    evaluation: Evaluation,
    evaluate: impl FnOnce(&mut CLIState, InspectPrecision) -> Result<T>,
) -> Result<T> {
    match evaluation {
        Evaluation::Arb => evaluate(cli, InspectPrecision::Arb),
        Evaluation::Production => evaluate(cli, InspectPrecision::StabilityLevels),
        Evaluation::Single(precision) => {
            let production_levels = stability_levels(cli, name)?;
            set_stability_levels(
                cli,
                name,
                &[StabilityLevelSetting {
                    precision,
                    three_dimensional_representation: None,
                    // No estimate reaches this tolerance, so the only level
                    // always accepts its own result.
                    required_precision_for_re: 1e99,
                    required_precision_for_im: 1e99,
                    ecm_relative_tolerance_for_re: 0.0,
                    ecm_relative_tolerance_for_im: 0.0,
                    escalate_for_large_weight_threshold: -1.0,
                }],
            )?;
            let result = evaluate(cli, InspectPrecision::StabilityLevels);
            set_stability_levels(cli, name, &production_levels)?;
            result
        }
    }
}

/// The runtime settings of the integrand `name` of the process `name`.
pub fn runtime_settings<'a>(cli: &'a CLIState, name: &str) -> Result<&'a RuntimeSettings> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_owned())))?;
    Ok(cli
        .state
        .process_list
        .get_integrand(process_id, name)?
        .require_generated()?
        .get_settings())
}

/// The stability levels of the integrand `name` of the process `name`.
pub fn stability_levels(cli: &CLIState, name: &str) -> Result<Vec<StabilityLevelSetting>> {
    Ok(runtime_settings(cli, name)?.stability.levels.clone())
}

/// Set the infrared regulator mass `general.m_ir` of the static theory's UV
/// counterterms of the integrand `name` of the process `name`, and check
/// that the integrand applies it.
pub fn set_static_theory_mass(cli: &mut CLIState, name: &str, m_ir: f64) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {name} -i {name} kv general.m_ir={m_ir}"
    ))
    .with_context(|| format!("setting general.m_ir of {name} to {m_ir}"))?;
    let kept = runtime_settings(cli, name)?.general.m_ir;
    if kept != m_ir {
        return Err(eyre!("{name} kept general.m_ir = {kept} instead of {m_ir}"));
    }
    Ok(())
}

/// Check that the run card of `name` sets the static-theory mass
/// `general.m_ir` to its `general.m_uv`, as the gates of every ring were
/// measured with, and return that mass.
pub fn assert_run_card_static_theory_mass(cli: &CLIState, name: &str) -> Result<f64> {
    let general = &runtime_settings(cli, name)?.general;
    if general.m_ir != general.m_uv {
        return Err(eyre!(
            "{name} starts with general.m_ir = {} apart from general.m_uv = {}",
            general.m_ir,
            general.m_uv
        ));
    }
    Ok(general.m_ir)
}

fn set_stability_levels(
    cli: &mut CLIState,
    name: &str,
    levels: &[StabilityLevelSetting],
) -> Result<()> {
    let levels_json = serde_json::to_string(levels)?;
    cli.run_command(&format!(
        "set process -p {name} -i {name} kv stability.levels='{levels_json}'"
    ))
    .with_context(|| format!("setting the stability levels of {name} to {levels_json}"))?;
    let applied = stability_levels(cli, name)?;
    if applied != levels {
        return Err(eyre!(
            "{name} kept the stability levels {applied:?} instead of {levels:?}"
        ));
    }
    Ok(())
}
