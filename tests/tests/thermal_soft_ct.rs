//! Slow lane of the thermal zero-mode soft subtraction (see "Finite-temperature
//! zero-mode subtraction" in `docs/architecture/uv-renormalization.typ`).
//!
//! Each ring of [`thermal_rings`] is generated with the zero-mode subtraction,
//! with vacuum subtraction (`run generate`) and, where measured to be
//! affordable, without it (`run generate_novacsub`). A bare control of the same
//! vacuum-subtraction setting, which keeps the Bose-enhanced divergence of the
//! static q-gluons, evaluates the identical rays:
//! - one `profile bulk` call per limit, with the control's `n_points`, and
//!   equal ray fingerprints;
//! - S(e0) integrable with r² ≥ 0.98, the control divergent, and an
//!   improvement by more than one power;
//! - every other soft or joint soft limit integrable with r² ≥ 0.98 and at its
//!   expected reading; the profiler's pass flag is never used for them;
//! - collinear limits integrable, without a fit-quality requirement;
//! - UV profiles in which every limit decays with a slope of at most -0.9
//!   (`UV_PROFILE_MAX_DOD`);
//! - the approach probes (a)-(j): a bounded remainder or its expected power in
//!   arbitrary precision, and f64 values that agree with it or whose errors
//!   track the control's on the same point.
//!
//! The `*_with_static_theory_mass` tests repeat the IR and UV gates of the
//! quark ring and of the ghost ring with and without vacuum subtraction, with
//! the static theory's UV counterterms at `general.m_ir` apart from the run
//! cards' `general.m_uv`.
//!
//! The `*_integrated_zero_mode_subtraction*` tests repeat every gate and probe
//! of a ring generation with integrated counterterms (`run generate_integrated`
//! and `run generate_integrated_novacsub`): every counterterm with `m_uv` adds
//! back its finite conversion to MS-bar, also inside the zero-mode terms, as
//! integrated prefixes followed by the local operations of their forest words.
//! The localized add-backs are part of the remainder, so it must stay locally
//! integrable with the same readings. The run cards localize them at the
//! thermal scale, where the profiles and probes see them at order one. On the
//! gluon ring the factorized join's add-back leaves a polynomial energy contact
//! on the bare gluon loop, which only vacuum subtraction removes: without it
//! generation must stop with the error that asks for vacuum subtraction. With
//! it the join localizes and the gluon ring runs every gate and probe like the
//! other rings.
//!
//! The `*_default_stability_levels` tests check the production stability
//! levels of the quark ring against arb at widely separated scales and next to
//! soft quarks, on both sides of the escalation thresholds, where a rotation
//! cannot detect the error of the remainder and only the a priori zero-mode
//! escalation tiers skip the lower precisions.
//!
//! The profiler fits the differences of successive samples, so a remainder
//! that tends to a nonzero constant (physical scaling +3) reads +4. The
//! expected readings below use that convention; arb evaluations along these
//! rays confirmed the physical +3 behind every +4.

use std::path::PathBuf;
use std::time::Instant;

use color_eyre::Result;
use gammaloop_integration_tests as integration_tests;
use gammaloop_integration_tests::thermal_rings;
use gammaloop_integration_tests::thermal_rings::{
    Counterterms, Evaluation, GHOST_RING, MIN_R_SQUARED, MIXED_RING, Probe, QUARK_RING, RingPoint,
    SCALING_TOLERANCE, SeparatedScalesRay, TADPOLE_RING, ThermalRing, UvScan, VacuumSubtraction,
};
use gammaloop_integration_tests::{CLIState, SERIAL_TEST_LOCK};
use gammalooprs::integrands::process::ir::SingleLimitReport;
use gammalooprs::settings::runtime::Precision;
use gammalooprs::uv::ApproximationType;
use gammalooprs::uv::profile::UVLimitSelection;
use spenso::algebra::complex::{Complex, FloatDerived};

use Approach::{Bounded, Exponent};
use ChemicalPotentials::{Both, ZeroOnly};
use Expectation::{BoseEnhanced, Collinear, Power};
use PrecisionCheck::{MatchesArb, QuadMatchesArb, TracksControl, TracksControlWithOffset};

/// A logarithmic endpoint can fit a tiny positive exponent on a finite ray.
/// These limits have integer power counting, so require a margin from zero.
const MIN_SOFT_SCALING: f64 = 0.5;
/// The bare control keeps the Bose-enhanced S(e0) divergence.
const MAX_CONTROL_SOFT_SCALING: f64 = -0.5;
/// The subtraction must cross the integrability boundary without an
/// integrated addback.
const MIN_SOFT_IMPROVEMENT: f64 = 1.0;

/// Pure f64 relative error on a probe whose remainder f64 must not degrade
/// towards the locus (measured at most 4.9e-13, quark (j) at muB = 3).
const MAX_F64_RELATIVE_ERROR: f64 = 1e-11;
/// Pure f128 relative error at large momenta (measured at most 6.1e-11).
const MAX_F128_RELATIVE_ERROR: f64 = 1e-8;
/// Bound on the f64 error of the remainder, relative to the larger of |R| and
/// |F|, over the bare control's relative f64 error on the same point
/// (measured at most 39, ghost (e) at muB = 0).
const MAX_CONTROL_ERROR_RATIO: f64 = 100.0;
/// The same bound where the static k-mode of a connected host dominates the
/// remainder: its rounding error grows at the bare's rate, one to three digits
/// above it (measured at most 2.2e2, ghost (f) at |k| = 1e-3 and muB = 0).
const MAX_CONTROL_ERROR_OFFSET: f64 = 1e4;
/// Relative f64 errors of the control below this floor count as exact in
/// those comparisons. Where the control is exact, [`MAX_CONTROL_ERROR_RATIO`]
/// then bounds the remainder's error by [`MAX_F64_RELATIVE_ERROR`], and
/// [`MAX_CONTROL_ERROR_OFFSET`] by a hundred times that.
const CONTROL_ERROR_FLOOR: f64 = MAX_F64_RELATIVE_ERROR / MAX_CONTROL_ERROR_RATIO;
/// Relative change of a bounded remainder between the two innermost probe
/// points (measured at most 1.4e-4); a divergence changes it at O(1).
const MAX_LOCUS_CHANGE: f64 = 1e-2;
/// Distance of a probe's local exponent from its expected power (measured at
/// most 6e-3).
const EXPONENT_TOLERANCE: f64 = 0.05;

/// The static-theory mass `general.m_ir` of the `*_with_static_theory_mass`
/// tests, apart from the run cards' `general.m_uv = 1`. The profiles fit the
/// successive differences of the remainder, so the r² gate fails where the
/// linear approach coefficient on the profiled ray vanishes by accident, and
/// the masses move these zeros: on the ghost ring's S(e3) ray with vacuum
/// subtraction at muB = 0, m_IR = 0.3 on U_G alone reads +3.95 at r² 0.938,
/// and the code before `general.m_ir` fails the same way at m_UV = 0.74
/// (+3.84 at r² 0.93), with a bounded remainder in both. At 0.5, with U_B on
/// the rows of S0_A at m_IR as well, the ray reads +3.996 at r² 0.99999.
const STATIC_THEORY_MASS: f64 = 0.5;

/// The static-theory mass of a slow-lane generation.
#[derive(Clone, Copy, Debug, PartialEq)]
enum StaticTheoryMass {
    /// The run card's `general.m_ir`, equal to its `general.m_uv`.
    RunCard,
    /// `general.m_ir` set to this value.
    Apart(f64),
}

impl StaticTheoryMass {
    /// The purpose that names the generation's state.
    fn purpose(self) -> &'static str {
        match self {
            Self::RunCard => "zero_mode",
            Self::Apart(_) => "zero_mode_static_theory_mass",
        }
    }

    /// The suffix that names this mass in a gate's context: none at the run
    /// card's mass.
    fn context_suffix(self) -> String {
        match self {
            Self::RunCard => String::new(),
            Self::Apart(m_ir) => format!(" (m_ir={m_ir})"),
        }
    }
}

/// The baryon chemical potentials of the slow lane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChemicalPotential {
    Zero,
    Three,
}

impl ChemicalPotential {
    const ALL: [Self; 2] = [Self::Zero, Self::Three];

    fn value(self) -> f64 {
        match self {
            Self::Zero => 0.0,
            Self::Three => 3.0,
        }
    }
}

/// The chemical potentials at which a row is checked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ChemicalPotentials {
    Both,
    ZeroOnly,
}

impl ChemicalPotentials {
    fn cover(self, chemical_potential: ChemicalPotential) -> bool {
        match self {
            Self::Both => true,
            Self::ZeroOnly => chemical_potential == ChemicalPotential::Zero,
        }
    }
}

/// What a profiled IR limit must read.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Expectation {
    /// S(e0), which the static q-gluons Bose-enhance: the remainder reads +1,
    /// the bare control diverges and the subtraction gains more than one power.
    BoseEnhanced,
    /// A soft or joint soft limit with this profiler reading.
    Power(f64),
    /// A collinear limit: integrable, without a power-law requirement.
    Collinear,
}

/// Limits profiled with one `n_points` at the chemical potentials of
/// `chemical_potentials`.
struct LimitGroup {
    n_points: usize,
    chemical_potentials: ChemicalPotentials,
    limits: &'static [(&'static str, Expectation)],
}

/// One profiled limit of a [`LimitGroup`].
#[derive(Clone, Copy, Debug, PartialEq)]
struct Limit {
    name: &'static str,
    expectation: Expectation,
    n_points: usize,
}

/// One UV profile of the subtraction, at every chemical potential.
#[derive(Clone, Copy, Debug, PartialEq)]
enum UvRow {
    /// The superficially divergent limits of the generation LMB, whose count is
    /// the ring's `uv_divergent_limits`.
    OnlyDivergent,
    /// Every loop subset in every LMB, scaled by 10^min .. 10^max.
    All {
        limits: usize,
        min_scale_exponent: f64,
        max_scale_exponent: f64,
    },
}

impl UvRow {
    fn scan(self) -> UvScan {
        match self {
            Self::OnlyDivergent => UvScan::ONLY_DIVERGENT,
            Self::All {
                min_scale_exponent,
                max_scale_exponent,
                ..
            } => UvScan {
                selection: UVLimitSelection::All,
                min_scale_exponent,
                max_scale_exponent,
            },
        }
    }

    fn limits(self, ring: &ThermalRing) -> usize {
        match self {
            Self::OnlyDivergent => ring.uv_divergent_limits,
            Self::All { limits, .. } => limits,
        }
    }
}

/// What the arbitrary-precision remainder does along a probe.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Approach {
    /// It stays bounded across the locus.
    Bounded,
    /// |R| scales with this power of the probe's momentum scale.
    Exponent(f64),
}

/// How the lower-precision remainder must agree with the arb reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PrecisionCheck {
    /// Pure f64 within [`MAX_F64_RELATIVE_ERROR`] at every point: no growth.
    MatchesArb,
    /// Pure f64 loses no more digits than the bare control on the same point
    /// ([`MAX_CONTROL_ERROR_RATIO`]), where the control loses digits itself.
    TracksControl,
    /// Pure f64 loses digits at the bare control's rate on the same points,
    /// at most [`MAX_CONTROL_ERROR_OFFSET`] above it: the static k-mode of a
    /// connected host, whose 1/k³ pair cancels in the remainder.
    TracksControlWithOffset,
    /// Pure f128 within [`MAX_F128_RELATIVE_ERROR`] at every point.
    QuadMatchesArb,
}

impl PrecisionCheck {
    /// The precision of the single stability level that evaluates the
    /// remainder.
    fn precision(self) -> Precision {
        match self {
            Self::QuadMatchesArb => Precision::Quad,
            Self::MatchesArb | Self::TracksControl | Self::TracksControlWithOffset => {
                Precision::Double
            }
        }
    }

    /// Whether the check compares with the bare control's values.
    fn uses_control(self) -> bool {
        matches!(self, Self::TracksControl | Self::TracksControlWithOffset)
    }

    /// The bound on every point's measure.
    fn bound(self) -> f64 {
        match self {
            Self::MatchesArb => MAX_F64_RELATIVE_ERROR,
            Self::TracksControl => MAX_CONTROL_ERROR_RATIO,
            Self::TracksControlWithOffset => MAX_CONTROL_ERROR_OFFSET,
            Self::QuadMatchesArb => MAX_F128_RELATIVE_ERROR,
        }
    }
}

/// One approach probe of a ring generation: what its arb remainder does and how
/// its lower-precision values must agree.
struct ProbeRow {
    probe: Probe,
    approach: Approach,
    precision: PrecisionCheck,
}

/// The power of its loop momentum with which the normalized kernel
/// `rls / (pi^2 (|k|^2 + rls^2)^2)` that localizes an integrated add-back on
/// its cograph decays: integrable over the three spatial components.
const LOCALIZATION_KERNEL_TAIL: f64 = -4.0;

impl ProbeRow {
    /// What the arb remainder does along this probe on a generation with
    /// `counterterms`. The add-backs of the second loop and of the join carry
    /// the localization kernel of the second loop's momentum, so along (i),
    /// `k -> ∞`, the remainder with integrated counterterms decays no faster
    /// than that kernel.
    fn approach(&self, counterterms: Counterterms) -> Approach {
        match (counterterms, self.probe, self.approach) {
            (Counterterms::Integrated, Probe::HardSecondLoop, Exponent(exponent)) => {
                Exponent(exponent.max(LOCALIZATION_KERNEL_TAIL))
            }
            (_, _, approach) => approach,
        }
    }
}

/// One ring generation of the slow lane and everything it is checked against.
struct RingCase {
    ring: &'static ThermalRing,
    vacuum_subtraction: VacuumSubtraction,
    limit_groups: &'static [LimitGroup],
    uv: &'static [UvRow],
    probes: &'static [ProbeRow],
    probe_chemical_potentials: ChemicalPotentials,
}

/// A point at which the default stability levels must return the arb value to
/// within their required precision.
struct ProductionGuard {
    ray: SeparatedScalesRay,
    parameter: f64,
    chemical_potential: ChemicalPotential,
}

/// The production stability levels of one ring generation, checked at
/// widely separated scales.
struct ProductionLadder {
    ring: &'static ThermalRing,
    vacuum_subtraction: VacuumSubtraction,
    guards: &'static [ProductionGuard],
}

const fn probe(probe: Probe, approach: Approach, precision: PrecisionCheck) -> ProbeRow {
    ProbeRow {
        probe,
        approach,
        precision,
    }
}

/// The collinear limits of a ring whose second loop has two edges.
const RING_COLLINEAR: &[(&str, Expectation)] = &[
    ("C[e2,e3]", Collinear),
    ("C[e4,e5]", Collinear),
    ("C[e0,e2]", Collinear),
    ("C[e0,e4]", Collinear),
    ("C[e0,e5]", Collinear),
];

// Quark ring: one union host. With vacuum subtraction the vacuum parts of the
// two bubbles leave S(e2) and S(e4) at +2.
const QUARK_LIMITS: &[LimitGroup] = &[
    LimitGroup {
        n_points: 8,
        chemical_potentials: Both,
        limits: &[
            ("S(e0)", BoseEnhanced),
            ("S(e2)", Power(2.0)),
            ("S(e3)", Power(4.0)),
            ("S(e4)", Power(2.0)),
            ("S(e5)", Power(4.0)),
            ("S(e0)S(e4)", Power(2.0)),
            ("S(e0)S(e5)", Power(2.0)),
            ("S(e4)S(e5)", Power(2.0)),
            ("S(e0)S(e2)", Power(2.0)),
            ("S(e0)S(e3)", Power(2.0)),
        ],
    },
    LimitGroup {
        n_points: 8,
        chemical_potentials: Both,
        limits: RING_COLLINEAR,
    },
];

const QUARK_LIMITS_WITHOUT_VACUUM_SUBTRACTION: &[LimitGroup] = &[LimitGroup {
    n_points: 8,
    chemical_potentials: Both,
    limits: &[
        ("S(e0)", BoseEnhanced),
        ("S(e4)", Power(4.0)),
        ("S(e0)S(e4)", Power(4.0)),
        ("S(e0)S(e5)", Power(4.0)),
        ("S(e4)S(e5)", Power(4.0)),
    ],
}];

const QUARK_PROBES: &[ProbeRow] = &[
    probe(Probe::CrossBubbleThreshold, Bounded, MatchesArb),
    probe(Probe::EqualQuarkEnergies, Bounded, MatchesArb),
    probe(Probe::CoincidentLoopMomenta, Bounded, TracksControl),
    // The vacuum term Π_A^v(0;p)/(4|q|³) of S(e2).
    probe(Probe::SoftQuarkAtSoftQ, Exponent(-1.0), TracksControl),
    probe(Probe::CollinearQuarkThreshold, Bounded, TracksControl),
    // The vacuum term Π_B^v(0;k) of S(e4).
    probe(Probe::SoftSecondLoop, Exponent(-1.0), TracksControl),
    probe(Probe::EqualSecondLoopEnergies, Bounded, MatchesArb),
    probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
    probe(Probe::HardSecondLoop, Exponent(-4.0), QuadMatchesArb),
    probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
];

/// Without vacuum subtraction both quark-ring remainders stay bounded at p -> 0
/// and k -> 0.
const QUARK_PROBES_WITHOUT_VACUUM_SUBTRACTION: &[ProbeRow] = &[
    probe(Probe::SoftQuarkAtSoftQ, Exponent(0.0), TracksControl),
    probe(Probe::SoftSecondLoop, Exponent(0.0), TracksControl),
    probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
    probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
];

const QUARK: RingCase = RingCase {
    ring: &QUARK_RING,
    vacuum_subtraction: VacuumSubtraction::On,
    limit_groups: QUARK_LIMITS,
    uv: &[
        UvRow::OnlyDivergent,
        UvRow::All {
            limits: 84,
            min_scale_exponent: 3.0,
            max_scale_exponent: 6.0,
        },
    ],
    probes: QUARK_PROBES,
    probe_chemical_potentials: Both,
};

const QUARK_WITHOUT_VACUUM_SUBTRACTION: RingCase = RingCase {
    ring: &QUARK_RING,
    vacuum_subtraction: VacuumSubtraction::Off,
    limit_groups: QUARK_LIMITS_WITHOUT_VACUUM_SUBTRACTION,
    uv: &[UvRow::OnlyDivergent],
    probes: QUARK_PROBES_WITHOUT_VACUUM_SUBTRACTION,
    probe_chemical_potentials: Both,
};

// Deterministic points of the quark ring where the rotated Double evaluation
// repeats most of the rounding error of the remainder, so its estimate misses
// it, and only the a priori tiers skip Double. They compare the softest massless
// line E_soft, the static q-gluons (E_zm = |q| = 0.58 T) or a softer quark,
// with Σ|k|: Quad below 5e-3 Σ|k|, the last level below 5e-6 Σ|k|. "Before"
// is the previous rule, E_zm alone against 1e-3 and 1e-6 Σ|k|.
const QUARK_LADDER: ProductionLadder = ProductionLadder {
    ring: &QUARK_RING,
    vacuum_subtraction: VacuumSubtraction::On,
    guards: &[
        // Without any hierarchy tier: Double accepted with estimate 4.3e-6,
        // true error 4.0e8.
        ProductionGuard {
            ray: SeparatedScalesRay::HardQuarkSoftSecondLoop,
            parameter: 1e-5,
            chemical_potential: ChemicalPotential::Zero,
        },
        // Without any hierarchy tier: Double accepted with estimate 4.4e-6,
        // true error 4.5e6.
        ProductionGuard {
            ray: SeparatedScalesRay::HardQuark,
            parameter: 3e5,
            chemical_potential: ChemicalPotential::Three,
        },
        // E_zm = 1.2e-3 Σ|k|. Before: Double accepted with estimate 4.7e-8,
        // true error 4.0e-5.
        ProductionGuard {
            ray: SeparatedScalesRay::HardQuark,
            parameter: 500.0,
            chemical_potential: ChemicalPotential::Zero,
        },
        // E_zm = 1.2e-3 Σ|k| with a quark |k| = 2e-3 T. Before: Double
        // accepted with estimate 8.9e-6, true error 3.0e-3.
        ProductionGuard {
            ray: SeparatedScalesRay::HardQuarkSoftSecondLoop,
            parameter: 2e-3,
            chemical_potential: ChemicalPotential::Zero,
        },
        // E_zm = 1.2e-6 Σ|k| with a quark |k| = 2e-6 T. Before: Quad accepted
        // with estimate 2.0e-9, true error 3.6e-4.
        ProductionGuard {
            ray: SeparatedScalesRay::HardQuarkSoftSecondLoop,
            parameter: 2e-6,
            chemical_potential: ChemicalPotential::Zero,
        },
        // E_zm = 5.8e-3 Σ|k|, above both thresholds, with a quark
        // |k| = 1e-4 T. Before: Double accepted with estimate 2.1e-8, true
        // error 1.3e-4.
        ProductionGuard {
            ray: SeparatedScalesRay::SoftSecondLoopAtHardQuark,
            parameter: 1e-4,
            chemical_potential: ChemicalPotential::Zero,
        },
    ],
};

const QUARK_LADDER_WITHOUT_VACUUM_SUBTRACTION: ProductionLadder = ProductionLadder {
    ring: &QUARK_RING,
    vacuum_subtraction: VacuumSubtraction::Off,
    guards: &[
        // E_zm = 1e-2 T, at the energy tier's threshold, with Σ|k| ≈ 1e2.
        // Without any hierarchy tier: Double accepted with estimates 9.5e-6
        // and 6.5e-6, true errors 5.1e-5 and 2.0e-5.
        ProductionGuard {
            ray: SeparatedScalesRay::SoftQHardSecondLoop,
            parameter: 1e-2,
            chemical_potential: ChemicalPotential::Zero,
        },
        ProductionGuard {
            ray: SeparatedScalesRay::SoftQHardSecondLoop,
            parameter: 1e-2,
            chemical_potential: ChemicalPotential::Three,
        },
        // E_zm = 5.8e-3 Σ|k| with a quark |k| = 1e-4 T. Before: Double
        // accepted with estimate 1.5e-6, true error 0.12.
        ProductionGuard {
            ray: SeparatedScalesRay::SoftSecondLoopAtHardQuark,
            parameter: 1e-4,
            chemical_potential: ChemicalPotential::Zero,
        },
    ],
};

// Ghost ring: union host and connected host A. Its UV `all` profile covers
// q -> ∞ at fixed k - q, which the mixed ring's `only-divergent` profile does
// not; on 1e3 .. 1e6 some of its rays are still pre-asymptotic (-0.91).
const GHOST_UV: &[UvRow] = &[
    UvRow::OnlyDivergent,
    UvRow::All {
        limits: 84,
        min_scale_exponent: 4.0,
        max_scale_exponent: 7.0,
    },
];

const GHOST: RingCase = RingCase {
    ring: &GHOST_RING,
    vacuum_subtraction: VacuumSubtraction::On,
    limit_groups: &[
        LimitGroup {
            n_points: 15,
            chemical_potentials: Both,
            limits: &[
                ("S(e0)", BoseEnhanced),
                ("S(e0)S(e4)", Power(1.0)),
                ("S(e0)S(e5)", Power(1.0)),
                ("S(e4)S(e5)", Power(1.0)),
                ("S(e0)S(e2)", Power(2.0)),
                ("S(e0)S(e3)", Power(2.0)),
            ],
        },
        LimitGroup {
            n_points: 10,
            chemical_potentials: Both,
            limits: &[
                ("S(e2)", Power(2.0)),
                ("S(e3)", Power(4.0)),
                ("S(e4)", Power(1.0)),
                ("S(e5)", Power(2.0)),
            ],
        },
        LimitGroup {
            n_points: 10,
            chemical_potentials: Both,
            limits: RING_COLLINEAR,
        },
    ],
    uv: GHOST_UV,
    probes: &[
        probe(Probe::CrossBubbleThreshold, Bounded, MatchesArb),
        probe(Probe::EqualQuarkEnergies, Bounded, MatchesArb),
        probe(Probe::CoincidentLoopMomenta, Bounded, MatchesArb),
        probe(Probe::SoftQuarkAtSoftQ, Exponent(-1.0), TracksControl),
        probe(Probe::CollinearQuarkThreshold, Bounded, TracksControl),
        // The static k-mode keeps S(e4) at +1; its 1/k³ pair cancels.
        probe(
            Probe::SoftSecondLoop,
            Exponent(-2.0),
            TracksControlWithOffset,
        ),
        probe(Probe::EqualSecondLoopEnergies, Bounded, MatchesArb),
        probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
        probe(Probe::HardSecondLoop, Exponent(-4.0), QuadMatchesArb),
        probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
    ],
    probe_chemical_potentials: Both,
};

// The ghost generation without vacuum subtraction takes about two minutes, so
// besides the quark and tadpole rings the ghost ring covers a connected host
// nested in a union host without vacuum subtraction.
const GHOST_WITHOUT_VACUUM_SUBTRACTION: RingCase = RingCase {
    ring: &GHOST_RING,
    vacuum_subtraction: VacuumSubtraction::Off,
    limit_groups: &[
        LimitGroup {
            n_points: 15,
            chemical_potentials: Both,
            limits: &[("S(e0)", BoseEnhanced), ("S(e0)S(e4)", Power(1.0))],
        },
        LimitGroup {
            n_points: 10,
            chemical_potentials: Both,
            limits: &[("S(e4)", Power(1.0))],
        },
    ],
    uv: GHOST_UV,
    probes: &[
        probe(Probe::SoftQuarkAtSoftQ, Exponent(0.0), TracksControl),
        probe(
            Probe::SoftSecondLoop,
            Exponent(-2.0),
            TracksControlWithOffset,
        ),
        probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
        probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
    ],
    probe_chemical_potentials: Both,
};

// Tadpole ring: union host A ∪ T only. There is no e5, so (a) and (g) do not
// apply: (a) reduces to (b) exactly and to (j) at linear order.
const TADPOLE_COLLINEAR: LimitGroup = LimitGroup {
    n_points: 10,
    chemical_potentials: Both,
    limits: &[
        ("C[e2,e3]", Collinear),
        ("C[e0,e2]", Collinear),
        ("C[e0,e4]", Collinear),
        ("C[e2,e4]", Collinear),
    ],
};

const TADPOLE_UV: &[UvRow] = &[
    UvRow::OnlyDivergent,
    UvRow::All {
        limits: 35,
        min_scale_exponent: 3.0,
        max_scale_exponent: 6.0,
    },
];

const TADPOLE: RingCase = RingCase {
    ring: &TADPOLE_RING,
    vacuum_subtraction: VacuumSubtraction::On,
    limit_groups: &[
        LimitGroup {
            n_points: 15,
            chemical_potentials: Both,
            limits: &[
                ("S(e0)", BoseEnhanced),
                ("S(e4)", Power(1.0)),
                ("S(e0)S(e4)", Power(2.0)),
                ("S(e2)", Power(2.0)),
                ("S(e3)", Power(4.0)),
                ("S(e0)S(e2)", Power(2.0)),
                ("S(e0)S(e3)", Power(2.0)),
            ],
        },
        TADPOLE_COLLINEAR,
    ],
    uv: TADPOLE_UV,
    probes: &[
        probe(Probe::EqualQuarkEnergies, Bounded, MatchesArb),
        probe(Probe::CoincidentLoopMomenta, Bounded, MatchesArb),
        probe(Probe::SoftQuarkAtSoftQ, Exponent(-1.0), TracksControl),
        probe(Probe::CollinearQuarkThreshold, Bounded, TracksControl),
        // The q-independent static tadpole T/k² stays unsubtracted, and its
        // f64 values stay exact.
        probe(Probe::SoftSecondLoop, Exponent(-2.0), MatchesArb),
        probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
        probe(Probe::HardSecondLoop, Exponent(-5.0), QuadMatchesArb),
        probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
    ],
    probe_chemical_potentials: Both,
};

// The tadpole generation without vacuum subtraction takes seconds, so it runs
// the whole matrix.
const TADPOLE_WITHOUT_VACUUM_SUBTRACTION: RingCase = RingCase {
    ring: &TADPOLE_RING,
    vacuum_subtraction: VacuumSubtraction::Off,
    limit_groups: &[
        LimitGroup {
            n_points: 15,
            chemical_potentials: Both,
            limits: &[
                ("S(e0)", BoseEnhanced),
                ("S(e4)", Power(1.0)),
                ("S(e0)S(e4)", Power(2.0)),
                ("S(e2)", Power(4.0)),
                ("S(e3)", Power(4.0)),
                ("S(e0)S(e2)", Power(4.0)),
                ("S(e0)S(e3)", Power(4.0)),
            ],
        },
        TADPOLE_COLLINEAR,
    ],
    uv: TADPOLE_UV,
    probes: &[
        probe(Probe::SoftQuarkAtSoftQ, Exponent(0.0), TracksControl),
        probe(Probe::SoftSecondLoop, Exponent(-2.0), MatchesArb),
        probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
        probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
    ],
    probe_chemical_potentials: Both,
};

// Mixed ring: union host and connected host A, with vacuum subtraction only;
// the quark, ghost and tadpole rings cover the generation without it. Each
// profile evaluates the subtraction in arb (about 10 s per point), so the
// matrix is trimmed: both chemical potentials for S(e0), the joints and UV,
// muB = 0 for the rest, all with 8 points.
const MIXED: RingCase = RingCase {
    ring: &MIXED_RING,
    vacuum_subtraction: VacuumSubtraction::On,
    limit_groups: &[
        LimitGroup {
            n_points: 8,
            chemical_potentials: Both,
            limits: &[
                ("S(e0)", BoseEnhanced),
                ("S(e0)S(e4)", Power(1.0)),
                ("S(e0)S(e5)", Power(1.0)),
                ("S(e4)S(e5)", Power(1.0)),
            ],
        },
        LimitGroup {
            n_points: 8,
            chemical_potentials: ZeroOnly,
            limits: &[
                ("S(e4)", Power(1.0)),
                ("S(e5)", Power(1.0)),
                ("S(e2)", Power(2.0)),
                ("S(e3)", Power(4.0)),
                ("S(e0)S(e2)", Power(2.0)),
                ("S(e0)S(e3)", Power(2.0)),
            ],
        },
        LimitGroup {
            n_points: 8,
            chemical_potentials: ZeroOnly,
            limits: RING_COLLINEAR,
        },
    ],
    uv: &[UvRow::OnlyDivergent],
    probes: &[
        probe(Probe::CrossBubbleThreshold, Bounded, MatchesArb),
        probe(Probe::EqualQuarkEnergies, Bounded, MatchesArb),
        probe(Probe::CoincidentLoopMomenta, Bounded, TracksControl),
        probe(Probe::SoftQuarkAtSoftQ, Exponent(-1.0), TracksControl),
        probe(Probe::CollinearQuarkThreshold, Bounded, TracksControl),
        probe(
            Probe::SoftSecondLoop,
            Exponent(-2.0),
            TracksControlWithOffset,
        ),
        probe(Probe::EqualSecondLoopEnergies, Bounded, MatchesArb),
        probe(Probe::SecondLoopOrthogonalToQ, Bounded, MatchesArb),
        probe(Probe::HardSecondLoop, Exponent(-4.0), QuadMatchesArb),
        probe(Probe::QuarkOrthogonalToQ, Bounded, MatchesArb),
    ],
    probe_chemical_potentials: ZeroOnly,
};

#[test]
fn thermal_quark_bubble_soft_ct_preserves_uv_finiteness() -> Result<()> {
    let _serial = SERIAL_TEST_LOCK.lock();
    let name = "thermal_quark_bubble";
    let (mut cli, state_path) = thermal_rings::fresh_soft_ir_cli(name, name)?;
    assert!(!cli.cli_settings.global.generation.uv.generate_integrated);
    assert_eq!(
        cli.cli_settings
            .global
            .generation
            .uv
            .renormalization_prescription
            .massless_power_divergent,
        ApproximationType::IR,
    );
    cli.run_command("run generate")?;

    for mu_b in [0.0, 3.0] {
        cli.run_command(&format!("set model muB={mu_b}"))?;
        // This two-loop graph isolates the UV regression caused by extending
        // a thermal soft subtraction away from its soft external momenta.
        // All three proper divergent cycle unions and the all-hard limit
        // must remain finite, including the alternative physical LMB.
        let verdict = thermal_rings::uv_verdict(&mut cli, name, UvScan::ONLY_DIVERGENT)?;
        assert_eq!(verdict.total, 4, "{verdict}");
        assert_eq!(
            verdict.failed, 0,
            "the thermal soft CT introduced a UV divergence in {name} at muB={mu_b}: {verdict}"
        );

        let soft = thermal_rings::profile_ir_limit(&mut cli, name, "S(e0)", 15)?;
        assert!(
            soft.passed
                && soft.scaling.is_finite()
                && soft.scaling > MIN_SOFT_SCALING
                && soft.power_law_fit.r_squared() >= MIN_R_SQUARED,
            "the bubble's soft remainder is not integrable at muB={mu_b}: {soft}"
        );
    }

    integration_tests::clean_test(state_path);
    Ok(())
}

#[test]
fn thermal_quark_ring_zero_mode_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&QUARK)
}

#[test]
fn thermal_quark_ring_zero_mode_subtraction_without_vacuum_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&QUARK_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_quark_ring_default_stability_levels() -> Result<()> {
    assert_default_stability_levels(&QUARK_LADDER)
}

#[test]
fn thermal_quark_ring_default_stability_levels_without_vacuum_subtraction() -> Result<()> {
    assert_default_stability_levels(&QUARK_LADDER_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_ghost_quark_ring_zero_mode_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&GHOST)
}

#[test]
fn thermal_ghost_quark_ring_zero_mode_subtraction_without_vacuum_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&GHOST_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_tadpole_quark_ring_zero_mode_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&TADPOLE)
}

#[test]
fn thermal_tadpole_quark_ring_zero_mode_subtraction_without_vacuum_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&TADPOLE_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_gluon_quark_ring_zero_mode_subtraction() -> Result<()> {
    assert_zero_mode_subtraction(&MIXED)
}

#[test]
fn thermal_quark_ring_integrated_zero_mode_subtraction() -> Result<()> {
    assert_integrated_zero_mode_subtraction(&QUARK)
}

#[test]
fn thermal_quark_ring_integrated_zero_mode_subtraction_without_vacuum_subtraction() -> Result<()> {
    assert_integrated_zero_mode_subtraction(&QUARK_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_ghost_quark_ring_integrated_zero_mode_subtraction() -> Result<()> {
    assert_integrated_zero_mode_subtraction(&GHOST)
}

#[test]
fn thermal_ghost_quark_ring_integrated_zero_mode_subtraction_without_vacuum_subtraction()
-> Result<()> {
    assert_integrated_zero_mode_subtraction(&GHOST_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_tadpole_quark_ring_integrated_zero_mode_subtraction() -> Result<()> {
    assert_integrated_zero_mode_subtraction(&TADPOLE)
}

#[test]
fn thermal_tadpole_quark_ring_integrated_zero_mode_subtraction_without_vacuum_subtraction()
-> Result<()> {
    assert_integrated_zero_mode_subtraction(&TADPOLE_WITHOUT_VACUUM_SUBTRACTION)
}

#[test]
fn thermal_gluon_quark_ring_integrated_zero_mode_subtraction() -> Result<()> {
    assert_integrated_zero_mode_subtraction(&MIXED)
}

/// Without vacuum subtraction the gluon ring's generation with integrated counterterms
/// stops at the localization of the factorized join's add-back `C_A C_B`: on the bare gluon
/// loop the three-gluon vertices leave a polynomial energy contact, which no 3D integrand
/// represents. The error names the join and its cograph and asks for vacuum subtraction.
#[test]
fn thermal_gluon_quark_ring_integrated_generation_requires_vacuum_subtraction() -> Result<()> {
    let _serial = SERIAL_TEST_LOCK.lock();
    let (mut cli, state_path) = MIXED_RING.fresh_cli(&format!(
        "{}_{}_integrated_generation",
        MIXED_RING.name,
        VacuumSubtraction::Off.state_suffix()
    ))?;
    let error = cli
        .run_command("run generate_integrated_novacsub")
        .expect_err("the join's add-back needs vacuum subtraction");
    let message = format!("{error:?}");
    assert!(
        message.contains("the integrated counterterm of {e2,e3,e4,e5}")
            && message.contains("on its cograph {e0,e1} without vacuum subtraction")
            && message.contains("`global.generation.medium.vacuum_subtraction = true`"),
        "the gluon ring without vacuum subtraction must report the energy contact of the join: \
         {message}"
    );
    integration_tests::clean_test(state_path);
    Ok(())
}

#[test]
fn thermal_quark_ring_zero_mode_subtraction_with_static_theory_mass() -> Result<()> {
    assert_gates_with_static_theory_mass(&QUARK)
}

#[test]
fn thermal_ghost_quark_ring_zero_mode_subtraction_with_static_theory_mass() -> Result<()> {
    assert_gates_with_static_theory_mass(&GHOST)
}

#[test]
fn thermal_ghost_quark_ring_zero_mode_subtraction_without_vacuum_subtraction_with_static_theory_mass()
-> Result<()> {
    assert_gates_with_static_theory_mass(&GHOST_WITHOUT_VACUUM_SUBTRACTION)
}

/// The bare control's profile of one limit at one chemical potential.
struct ControlProfile {
    chemical_potential: ChemicalPotential,
    limit: &'static str,
    report: SingleLimitReport,
}

/// The bare control's arb and pure f64 values along one probe.
struct ControlProbe {
    chemical_potential: ChemicalPotential,
    probe: Probe,
    arb: Vec<Complex<f64>>,
    double: Vec<Complex<f64>>,
}

/// What the bare control evaluated on the rays and points the subtraction is
/// compared on.
struct Control {
    profiles: Vec<ControlProfile>,
    probes: Vec<ControlProbe>,
}

impl Control {
    fn profile(
        &self,
        chemical_potential: ChemicalPotential,
        limit: &str,
    ) -> Option<&SingleLimitReport> {
        self.profiles
            .iter()
            .find(|profile| {
                profile.chemical_potential == chemical_potential && profile.limit == limit
            })
            .map(|profile| &profile.report)
    }

    fn probe(&self, chemical_potential: ChemicalPotential, probe: Probe) -> Option<&ControlProbe> {
        self.probes.iter().find(|control| {
            control.chemical_potential == chemical_potential && control.probe == probe
        })
    }
}

impl RingCase {
    fn covers(&self, chemical_potential: ChemicalPotential) -> bool {
        self.limit_groups
            .iter()
            .any(|group| group.chemical_potentials.cover(chemical_potential))
    }

    fn limits_at(&self, chemical_potential: ChemicalPotential) -> impl Iterator<Item = Limit> {
        let groups: &'static [LimitGroup] = self.limit_groups;
        groups
            .iter()
            .filter(move |group| group.chemical_potentials.cover(chemical_potential))
            .flat_map(|group| {
                group.limits.iter().map(|&(name, expectation)| Limit {
                    name,
                    expectation,
                    n_points: group.n_points,
                })
            })
    }

    fn probes_at(&self, chemical_potential: ChemicalPotential) -> &'static [ProbeRow] {
        if self.probe_chemical_potentials.cover(chemical_potential) {
            self.probes
        } else {
            &[]
        }
    }
}

/// Generate the zero-mode subtraction of `ring` with `counterterms` in a fresh
/// state named after `purpose`, and check its settings and structure. Returns
/// the state and its path.
fn generate_subtraction(
    ring: &ThermalRing,
    vacuum_subtraction: VacuumSubtraction,
    counterterms: Counterterms,
    purpose: &str,
) -> Result<(CLIState, PathBuf)> {
    let (mut cli, state_path) = ring.fresh_cli(&format!(
        "{}_{}_{purpose}",
        ring.name,
        vacuum_subtraction.state_suffix()
    ))?;
    integration_tests::run_commands(
        &mut cli,
        vacuum_subtraction.subtraction_commands(counterterms),
    )?;
    let generation = &cli.cli_settings.global.generation;
    assert!(generation.uv.subtract_uv);
    assert_eq!(generation.uv.generate_integrated, counterterms.integrated());
    assert_eq!(
        generation
            .uv
            .renormalization_prescription
            .massless_power_divergent,
        ApproximationType::IR,
    );
    assert_eq!(
        generation.medium.vacuum_subtraction,
        vacuum_subtraction.enabled()
    );
    thermal_rings::assert_run_card_static_theory_mass(&cli, ring.name)?;
    if counterterms.integrated() {
        // At the default scale the add-backs would be negligible at momenta
        // of order T, and the UV profiles would cross the kernel's scale.
        let general = &thermal_rings::runtime_settings(&cli, ring.name)?.general;
        assert_eq!(
            general.renormalization_localization_scale * general.inverse_temperature,
            1.0,
            "the run card of {} localizes the integrated add-backs at the thermal scale",
            ring.name
        );
    }
    thermal_rings::assert_ring_structure(&cli, ring)?;
    Ok((cli, state_path))
}

fn assert_zero_mode_subtraction(case: &RingCase) -> Result<()> {
    assert_zero_mode_subtraction_with(case, Counterterms::Local, StaticTheoryMass::RunCard)
}

/// Every gate and probe of `case` on its generation with integrated
/// counterterms: the localized add-backs and their prefix sectors must leave
/// the remainder integrable with the readings of the local generation, but
/// along (i), where it decays like the localization kernel
/// ([`ProbeRow::approach`]).
fn assert_integrated_zero_mode_subtraction(case: &RingCase) -> Result<()> {
    assert_zero_mode_subtraction_with(case, Counterterms::Integrated, StaticTheoryMass::RunCard)
}

/// The IR and UV gates of `case`, without its approach probes, with
/// `general.m_ir` at [`STATIC_THEORY_MASS`]: the static theory's UV
/// counterterms then regulate the static soft lines with their own mass, and
/// every other counterterm keeps `general.m_uv`.
fn assert_gates_with_static_theory_mass(case: &RingCase) -> Result<()> {
    assert_zero_mode_subtraction_with(
        &RingCase {
            probes: &[],
            ..*case
        },
        Counterterms::Local,
        StaticTheoryMass::Apart(STATIC_THEORY_MASS),
    )
}

fn assert_zero_mode_subtraction_with(
    case: &RingCase,
    counterterms: Counterterms,
    static_theory_mass: StaticTheoryMass,
) -> Result<()> {
    let _serial = SERIAL_TEST_LOCK.lock();
    let started = Instant::now();
    let ring = case.ring;
    let name = ring.name;
    let purpose = match counterterms {
        Counterterms::Local => static_theory_mass.purpose().to_owned(),
        Counterterms::Integrated => format!("integrated_{}", static_theory_mass.purpose()),
    };

    let control = run_control(case, &purpose)?;
    let control_time = started.elapsed();

    let generation_started = Instant::now();
    let (mut cli, state_path) =
        generate_subtraction(ring, case.vacuum_subtraction, counterterms, &purpose)?;
    if let StaticTheoryMass::Apart(m_ir) = static_theory_mass {
        thermal_rings::set_static_theory_mass(&mut cli, name, m_ir)?;
    }
    let generation_time = generation_started.elapsed();

    for chemical_potential in ChemicalPotential::ALL {
        if !case.covers(chemical_potential) {
            continue;
        }
        let mu_b = chemical_potential.value();
        let context = format!(
            "{name} {}, {counterterms} at muB={mu_b}{}",
            case.vacuum_subtraction,
            static_theory_mass.context_suffix()
        );
        cli.run_command(&format!("set model muB={mu_b}"))?;

        let profiles_started = Instant::now();
        for limit in case.limits_at(chemical_potential) {
            let profile =
                thermal_rings::profile_ir_limit(&mut cli, name, limit.name, limit.n_points)?;
            let bare = control
                .profile(chemical_potential, limit.name)
                .unwrap_or_else(|| panic!("the control profiled {} for {context}", limit.name));
            assert_limit(&context, limit, &profile, bare);
        }
        for &uv in case.uv {
            let verdict = thermal_rings::uv_verdict(&mut cli, name, uv.scan())?;
            println!("{context}: UV {uv:?}: {verdict}");
            assert_eq!(verdict.total, uv.limits(ring), "{context}: {verdict}");
            assert_eq!(
                verdict.failed, 0,
                "the zero-mode subtraction leaves a UV divergence in {context}: {verdict}"
            );
        }
        let profiles_time = profiles_started.elapsed();

        let probes_started = Instant::now();
        for row in case.probes_at(chemical_potential) {
            let control_probe = control.probe(chemical_potential, row.probe);
            assert_probe(&mut cli, &context, name, row, counterterms, control_probe)?;
        }
        println!(
            "{context}: profiles {:.1} s, probes {:.1} s",
            profiles_time.as_secs_f64(),
            probes_started.elapsed().as_secs_f64(),
        );
    }
    integration_tests::clean_test(state_path);

    println!(
        "{name} {}, {counterterms}: total {:.1} s; control {:.1} s; subtraction generation \
         {:.1} s",
        case.vacuum_subtraction,
        started.elapsed().as_secs_f64(),
        control_time.as_secs_f64(),
        generation_time.as_secs_f64(),
    );
    Ok(())
}

/// Generate the bare control of `case` in a state named after the purpose of
/// the subtraction it controls, and evaluate every limit and every probe point
/// that the subtraction is compared with.
fn run_control(case: &RingCase, purpose: &str) -> Result<Control> {
    let ring = case.ring;
    let name = ring.name;
    let (mut cli, state_path) = ring.fresh_cli(&format!(
        "{name}_{}_{purpose}_control",
        case.vacuum_subtraction.state_suffix()
    ))?;
    integration_tests::run_commands(&mut cli, case.vacuum_subtraction.bare_commands())?;
    let generation = &cli.cli_settings.global.generation;
    assert!(!generation.uv.subtract_uv);
    assert_eq!(
        generation.medium.vacuum_subtraction,
        case.vacuum_subtraction.enabled()
    );
    thermal_rings::assert_ring_structure(&cli, ring)?;

    let mut control = Control {
        profiles: Vec::new(),
        probes: Vec::new(),
    };
    for chemical_potential in ChemicalPotential::ALL {
        if !case.covers(chemical_potential) {
            continue;
        }
        let mu_b = chemical_potential.value();
        cli.run_command(&format!("set model muB={mu_b}"))?;
        for limit in case.limits_at(chemical_potential) {
            control.profiles.push(ControlProfile {
                chemical_potential,
                limit: limit.name,
                report: thermal_rings::profile_ir_limit(
                    &mut cli,
                    name,
                    limit.name,
                    limit.n_points,
                )?,
            });
        }
        for row in case
            .probes_at(chemical_potential)
            .iter()
            .filter(|row| row.precision.uses_control())
        {
            let points = probe_points(row.probe);
            let arb =
                thermal_rings::evaluate_ring_points(&mut cli, name, Evaluation::Arb, &points)?;
            // A vanishing reference would make every relative error infinite
            // and the comparison vacuous.
            assert!(
                arb.iter()
                    .all(|value| value.norm().is_finite() && value.norm() > 0.0),
                "the bare control along {} must stay finite and nonzero for {name} {} at \
                 muB={mu_b}: {arb:?}",
                row.probe.label(),
                case.vacuum_subtraction,
            );
            control.probes.push(ControlProbe {
                chemical_potential,
                probe: row.probe,
                arb,
                double: thermal_rings::evaluate_ring_points(
                    &mut cli,
                    name,
                    Evaluation::Single(Precision::Double),
                    &points,
                )?,
            });
        }
    }
    integration_tests::clean_test(state_path);
    Ok(control)
}

fn assert_limit(
    context: &str,
    limit: Limit,
    profile: &SingleLimitReport,
    bare: &SingleLimitReport,
) {
    let Limit {
        name: limit,
        expectation,
        ..
    } = limit;
    let scaling = profile.scaling;
    let r_squared = profile.power_law_fit.r_squared();
    println!(
        "{context}: {limit}: {scaling:+.4} (r² {r_squared:.6}), expected {expectation:?}; \
         control {:+.4}; ray {:?}",
        bare.scaling, profile.ray_fingerprint,
    );
    assert!(bare.ray_fingerprint.is_some());
    assert_eq!(
        bare.ray_fingerprint, profile.ray_fingerprint,
        "the bare control and the subtraction must evaluate the identical {limit} ray in \
         {context}"
    );
    assert!(
        scaling.is_finite() && scaling > MIN_SOFT_SCALING,
        "{limit} is not integrable in {context}: {profile}"
    );
    let expected = match expectation {
        Collinear => return,
        BoseEnhanced => {
            // The control must retain the divergence, and the local
            // subtraction must cross the integrability boundary.
            assert!(
                bare.scaling < MAX_CONTROL_SOFT_SCALING,
                "the bare control lost the {limit} divergence in {context}: {bare}"
            );
            assert!(
                profile.passed && scaling - bare.scaling > MIN_SOFT_IMPROVEMENT,
                "the zero-mode subtraction did not improve {limit} in {context}: bare={bare}; \
                 subtracted={profile}"
            );
            1.0
        }
        Power(expected) => expected,
    };
    assert!(
        r_squared >= MIN_R_SQUARED && (scaling - expected).abs() <= SCALING_TOLERANCE,
        "{limit} must read {expected:+} with r² >= {MIN_R_SQUARED} in {context}: {profile}"
    );
}

/// The local exponent of the arb remainder `arb` along `probe`, read off its
/// two points closest to the locus or limit, for a remainder expected to scale
/// with the power `expected`. With integrated counterterms a remainder that
/// grows towards a locus at finite momenta is read off the successive
/// differences of its three points instead: the localized add-backs add a part
/// that stays regular there, constant in the approach, which the differences
/// cancel and which would otherwise bias the reading. The probe's parameters
/// are evenly spaced in their logarithm, so the differences scale like the
/// remainder. Towards infinite momentum, along (i), the add-backs vanish with
/// the rest.
fn local_exponent(
    probe: Probe,
    counterterms: Counterterms,
    expected: f64,
    arb: &[Complex<f64>],
) -> f64 {
    let &[first_parameter, outer_parameter, inner_parameter] = probe
        .parameters()
        .last_chunk()
        .expect("every probe has three parameters");
    let &[first, outer, inner] = arb
        .last_chunk()
        .expect("every probe parameter has one arb value");
    let scale =
        (probe.momentum_scale(inner_parameter) / probe.momentum_scale(outer_parameter)).ln();
    if counterterms.integrated() && expected < 0.0 && scale < 0.0 {
        let step = (outer_parameter / first_parameter).ln();
        assert!(
            (step - (inner_parameter / outer_parameter).ln()).abs() <= 1e-9 * step.abs(),
            "the parameters of {} are evenly spaced in their logarithm",
            probe.label()
        );
        ((inner - outer).norm() / (outer - first).norm()).ln() / scale
    } else {
        (inner.norm() / outer.norm()).ln() / scale
    }
}

fn probe_points(probe: Probe) -> Vec<RingPoint> {
    probe
        .parameters()
        .iter()
        .map(|&parameter| probe.point(parameter))
        .collect()
}

fn relative_error(value: Complex<f64>, reference: Complex<f64>) -> f64 {
    (value - reference).norm() / reference.norm()
}

fn scientific(values: &[f64]) -> String {
    let entries = values
        .iter()
        .map(|value| format!("{value:.1e}"))
        .collect::<Vec<_>>();
    format!("[{}]", entries.join(", "))
}

/// Check one probe on the subtraction: its arb shape and its lower-precision
/// agreement, against the control's values on the same points when it tracks
/// the control.
fn assert_probe(
    cli: &mut CLIState,
    context: &str,
    name: &str,
    row: &ProbeRow,
    counterterms: Counterterms,
    control: Option<&ControlProbe>,
) -> Result<()> {
    let probe = row.probe;
    let label = probe.label();
    let points = probe_points(probe);
    let arb = thermal_rings::evaluate_ring_points(cli, name, Evaluation::Arb, &points)?;
    assert!(
        arb.iter()
            .all(|value| value.norm().is_finite() && value.norm() > 0.0),
        "{label} must stay finite and nonzero in {context}: {arb:?}"
    );

    // The shape is read from the two points closest to the locus or limit.
    let &[outer, inner] = arb
        .last_chunk()
        .expect("every probe parameter has one arb value");
    match row.approach(counterterms) {
        Bounded => {
            let change = (inner - outer).norm() / inner.norm().max(outer.norm());
            println!("{context}: {label}: arb {arb:?}, inner change {change:.2e}");
            assert!(
                change <= MAX_LOCUS_CHANGE,
                "the remainder is not bounded across {label} in {context}: {arb:?}"
            );
        }
        Exponent(expected) => {
            let exponent = local_exponent(probe, counterterms, expected, &arb);
            println!("{context}: {label}: arb {arb:?}, local exponent {exponent:+.4}");
            assert!(
                (exponent - expected).abs() <= EXPONENT_TOLERANCE,
                "the remainder along {label} must scale with power {expected:+} in {context}: \
                 {exponent:+.4}"
            );
        }
    }

    // Every check evaluates the remainder with one stability level.
    let precision = row.precision.precision();
    let values =
        thermal_rings::evaluate_ring_points(cli, name, Evaluation::Single(precision), &points)?;
    let (measure, measured) = match row.precision {
        MatchesArb | QuadMatchesArb => {
            let errors = values
                .iter()
                .zip(&arb)
                .map(|(&value, &reference)| relative_error(value, reference))
                .collect::<Vec<_>>();
            ("relative errors", errors)
        }
        TracksControl | TracksControlWithOffset => {
            let control =
                control.unwrap_or_else(|| panic!("the control evaluated {label} for {context}"));
            // Normalizing by the larger of |R| and |F| compares absolute
            // errors where the subtraction cancels F, and relative errors
            // where its counterterms exceed F.
            let ratios = values
                .iter()
                .zip(&arb)
                .zip(control.double.iter().zip(&control.arb))
                .map(|((&value, &reference), (&bare, &bare_reference))| {
                    let error =
                        (value - reference).norm() / reference.norm().max(bare_reference.norm());
                    error / relative_error(bare, bare_reference).max(CONTROL_ERROR_FLOOR)
                })
                .collect::<Vec<_>>();
            ("error ratios to the control", ratios)
        }
    };
    let bound = row.precision.bound();
    let within = measured.iter().all(|&value| value <= bound);
    let measured = scientific(&measured);
    println!("{context}: {label}: {precision} {measure} {measured}");
    assert!(
        within,
        "{precision} along {label} exceeds {bound:e} in {context}: {measure} {measured}"
    );
    Ok(())
}

/// Check the default stability levels of `ladder`'s ring generation against
/// arb at every guard. Every guard is evaluated before any failure is
/// reported.
fn assert_default_stability_levels(ladder: &ProductionLadder) -> Result<()> {
    let _serial = SERIAL_TEST_LOCK.lock();
    let name = ladder.ring.name;
    let (mut cli, state_path) = generate_subtraction(
        ladder.ring,
        ladder.vacuum_subtraction,
        Counterterms::Local,
        "default_stability_levels",
    )?;
    let tolerance = production_tolerance(&cli, name)?;
    let mut failures = Vec::new();
    for guard in ladder.guards {
        let mu_b = guard.chemical_potential.value();
        let context = format!("{name} {} at muB={mu_b}", ladder.vacuum_subtraction);
        cli.run_command(&format!("set model muB={mu_b}"))?;
        failures.extend(production_guard_failure(
            &mut cli, &context, name, guard, tolerance,
        )?);
    }
    assert!(
        failures.is_empty(),
        "the default stability levels accept wrong values: {}",
        failures.join("; ")
    );
    integration_tests::clean_test(state_path);
    Ok(())
}

/// The loosest required precision of the integrand's stability levels: no
/// result they accept may be further from the arb value.
fn production_tolerance(cli: &CLIState, name: &str) -> Result<f64> {
    Ok(thermal_rings::stability_levels(cli, name)?
        .iter()
        .flat_map(|level| {
            [
                level.required_precision_for_re,
                level.required_precision_for_im,
            ]
        })
        .fold(0.0, f64::max))
}

/// Evaluate `guard` with the default stability levels and in arb, and
/// describe the failure when they disagree by more than `tolerance`.
fn production_guard_failure(
    cli: &mut CLIState,
    context: &str,
    name: &str,
    guard: &ProductionGuard,
    tolerance: f64,
) -> Result<Option<String>> {
    let point = guard.ray.point(guard.parameter).flattened();
    let production =
        thermal_rings::evaluate_momentum_point(cli, name, Evaluation::Production, &point)?;
    let reference = thermal_rings::evaluate_momentum_point(cli, name, Evaluation::Arb, &point)?;
    let error = relative_error(production, reference);
    let location = format!("{:?}({:e}) in {context}", guard.ray, guard.parameter);
    println!("{location}: production relative error {error:.1e} (tolerance {tolerance:e})");
    Ok((error > tolerance)
        .then(|| format!("{location}: {production} against arb {reference} ({error:.1e})")))
}
