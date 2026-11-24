//! Bare power counting of the thermal two-bubble rings.
//!
//! Each ring joins a d-quark polarization to a second polarization loop by two
//! gluons carrying q, with loop momenta [q, p, k] = [e0, e2, e4]. Without any
//! UV or soft subtraction the Matsubara zero mode of the q-gluons makes S(e0)
//! non-integrable, while the statistics of the second loop fix its own soft
//! limits. These scalings are the controls for the thermal soft counterterm.
//! Power counting cannot see a ring's overall sign or symmetry factor, so the
//! hand-written ghost and tadpole rings are also compared with the
//! generator's diagrams.
//!
//! With UV subtraction, the MUV prescription of massless power divergences
//! leaves the static soft set of the q-gluons without a zero-mode host, which
//! generation rejects before any forest work. The tadpole ring's zero-mode
//! subtraction without vacuum subtraction, which generates in seconds, checks
//! that the static-theory mass `general.m_ir` is independent of `general.m_uv`.

use color_eyre::Result;
use gammaloop_integration_tests as integration_tests;
use gammaloop_integration_tests::CLIState;
use gammaloop_integration_tests::evaluation_comparison;
use gammaloop_integration_tests::thermal_rings;
use gammaloop_integration_tests::thermal_rings::{
    Counterterms, Evaluation, GHOST_RING, MIN_R_SQUARED, MIXED_RING, QUARK_RING, SCALING_TOLERANCE,
    TADPOLE_RING, ThermalRing, UvScan, VacuumSubtraction,
};

use Expected::{Collinear, Power};

/// Collinear rays approach a smooth bare region; their fits are not exact powers.
const MIN_COLLINEAR_SCALING: f64 = 0.5;
/// Relative change of the subtracted tadpole ring at [`SOFT_Q_POINT`] when
/// `general.m_ir` moves from `m_uv = 1` to 0.3 (measured 4.5e-4 at muB = 0).
const MIN_STATIC_THEORY_MASS_CHANGE: f64 = 1e-5;

#[derive(Clone, Copy)]
enum Expected {
    /// A soft or joint soft limit with this exact power.
    Power(f64),
    /// A collinear limit, integrable without a power-law requirement.
    Collinear,
}

/// The bare limits profiled on one ring and their expected scalings.
struct BareRing {
    ring: &'static ThermalRing,
    limits: &'static [(&'static str, Expected)],
}

// Vacuum subtraction changes no soft or collinear limit of the bare ring.
const QUARK: BareRing = BareRing {
    ring: &QUARK_RING,
    limits: &[
        ("S(e0)", Power(-1.0)),
        ("S(e2)", Power(4.0)),
        ("S(e3)", Power(4.0)),
        ("S(e4)", Power(4.0)),
        ("S(e5)", Power(4.0)),
        ("S(e0)S(e4)", Power(2.0)),
        ("S(e0)S(e5)", Power(2.0)),
        ("S(e4)S(e5)", Power(2.0)),
        ("C[e2,e3]", Collinear),
        ("C[e4,e5]", Collinear),
        ("C[e0,e2]", Collinear),
        ("C[e0,e4]", Collinear),
        ("C[e0,e5]", Collinear),
    ],
};

// The complete 13-limit set of the gluon bubble costs about 80 s per chemical
// potential; this subset keeps both Bose limits and their joint.
const MIXED: BareRing = BareRing {
    ring: &MIXED_RING,
    limits: &[
        ("S(e0)", Power(-1.0)),
        ("S(e4)", Power(1.0)),
        ("S(e0)S(e4)", Power(0.0)),
        ("C[e0,e4]", Collinear),
    ],
};

// Ghosts have periodic Matsubara frequencies: S(e4) = +2 pins Bose statistics,
// whereas antiperiodic ghosts would give the quark-loop value +4.
const GHOST: BareRing = BareRing {
    ring: &GHOST_RING,
    limits: &[
        ("S(e0)", Power(-1.0)),
        ("S(e2)", Power(4.0)),
        ("S(e3)", Power(4.0)),
        ("S(e4)", Power(2.0)),
        ("S(e5)", Power(2.0)),
        ("S(e0)S(e4)", Power(0.0)),
        ("S(e0)S(e5)", Power(0.0)),
        ("S(e4)S(e5)", Power(0.0)),
        ("C[e2,e3]", Collinear),
        ("C[e4,e5]", Collinear),
        ("C[e0,e2]", Collinear),
        ("C[e0,e4]", Collinear),
        ("C[e0,e5]", Collinear),
    ],
};

const TADPOLE: BareRing = BareRing {
    ring: &TADPOLE_RING,
    limits: &[
        ("S(e0)", Power(-1.0)),
        ("S(e2)", Power(4.0)),
        ("S(e3)", Power(4.0)),
        ("S(e4)", Power(1.0)),
        ("S(e0)S(e4)", Power(0.0)),
        ("C[e2,e3]", Collinear),
        ("C[e0,e2]", Collinear),
        ("C[e0,e4]", Collinear),
        ("C[e2,e4]", Collinear),
    ],
};

/// The UV limits that the bare generation leaves convergent: with vacuum
/// subtraction only the all-hard limit decays, without it none does.
fn convergent_bare_uv_limits(vacuum_subtraction: VacuumSubtraction) -> usize {
    match vacuum_subtraction {
        VacuumSubtraction::On => 1,
        VacuumSubtraction::Off => 0,
    }
}

#[test]
fn thermal_quark_ring_bare_power_counting() -> Result<()> {
    assert_bare_power_counting(&QUARK, VacuumSubtraction::On)
}

#[test]
fn thermal_quark_ring_bare_power_counting_without_vacuum_subtraction() -> Result<()> {
    assert_bare_power_counting(&QUARK, VacuumSubtraction::Off)
}

#[test]
fn thermal_gluon_quark_ring_bare_power_counting() -> Result<()> {
    assert_bare_power_counting(&MIXED, VacuumSubtraction::On)
}

#[test]
fn thermal_ghost_quark_ring_bare_power_counting() -> Result<()> {
    assert_bare_power_counting(&GHOST, VacuumSubtraction::On)
}

#[test]
fn thermal_tadpole_quark_ring_bare_power_counting() -> Result<()> {
    assert_bare_power_counting(&TADPOLE, VacuumSubtraction::On)
}

/// The `generate_muv` block of every ring card fails in the scheme validation,
/// before any forest work, with the error that names the IR prescription.
#[test]
fn muv_prescription_reports_the_missing_zero_mode_host() -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    for ring in [&QUARK_RING, &MIXED_RING, &GHOST_RING, &TADPOLE_RING] {
        thermal_rings::assert_muv_prescription_is_hostless(
            ring,
            &format!("{}_muv_prescription", ring.name),
        )?;
    }
    Ok(())
}

/// Loop momenta [q, p, k] = [e0, e2, e4] of the hand-written rings.
const RING_POINT: [f64; 9] = [0.31, -0.12, 0.27, 0.71, 0.23, 0.43, -0.37, 0.61, 0.29];

/// One hand-written ring, the generator's diagram `GL<id>` that equals it, and
/// the diagram's loop momenta at `RING_POINT`.
struct GeneratedRing {
    ring: &'static ThermalRing,
    diagram: usize,
    point: [f64; 9],
}

const GENERATED_RINGS: [GeneratedRing; 2] = [
    // Diagram loop momenta [k - q, k, p - q].
    GeneratedRing {
        ring: &GHOST_RING,
        diagram: 30,
        point: [-0.68, 0.73, 0.02, -0.37, 0.61, 0.29, 0.40, 0.35, 0.16],
    },
    // Diagram loop momenta [p, p - q, k].
    GeneratedRing {
        ring: &TADPOLE_RING,
        diagram: 64,
        point: [0.71, 0.23, 0.43, 0.40, 0.35, 0.16, -0.37, 0.61, 0.29],
    },
];

#[test]
fn thermal_ring_graphs_match_the_generated_diagrams() -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    // Any ring card provides the model and the thermal medium.
    let (mut cli, state_path) =
        GHOST_RING.fresh_cli("thermal_ring_graphs_match_the_generated_diagrams")?;
    let generated = state_path.join("feyngen");
    // All three-loop vacuum diagrams of d, g and ghG. With one worker their
    // numbering is deterministic; GL30 and GL64 are the ghost and tadpole rings.
    integration_tests::run_commands(
        &mut cli,
        &[
            "import model sm-thermal",
            "set model aS=1.0 muB=0.0",
            "set global kv global.generation.uv.subtract_uv=false global.n_cores.feyngen=1",
            concat!(
                "generate amp {} > {} | d g ghG [{3}] --only-diagrams ",
                "--clear-existing-processes --numerator-grouping no_grouping ",
                "--max-n-bridges -1 --number-of-factorized-loop-subtopologies -1 -1 ",
                "--filter-zero-flow-edges false --filter-self-loop false ",
                "--filter-snails false --filter-tadpoles false --filter-selfenergies false ",
                "-p vac3 -i vac3",
            ),
            &format!("save dot {}", generated.display()),
        ],
    )?;
    for generated_ring in &GENERATED_RINGS {
        let name = generated_ring.ring.name;
        let diagram = generated.join(format!(
            "processes/amplitudes/vac3/vac3/GL{}.dot",
            generated_ring.diagram
        ));
        integration_tests::run_commands(
            &mut cli,
            &[
                &format!("import graphs ./tests/resources/graphs/{name}.dot -p {name} -i {name}"),
                &format!("generate existing -p {name} -i {name}"),
                &format!(
                    "import graphs {} -p {name}_generated -i {name}_generated",
                    diagram.display()
                ),
                &format!("generate existing -p {name}_generated -i {name}_generated"),
            ],
        )?;
    }
    for mu_b in [0.0, 3.0] {
        cli.run_command(&format!("set model muB={mu_b}"))?;
        for generated_ring in &GENERATED_RINGS {
            let name = generated_ring.ring.name;
            let ring = thermal_rings::evaluate_momentum_point(
                &mut cli,
                name,
                Evaluation::Production,
                &RING_POINT,
            )?;
            let diagram = thermal_rings::evaluate_momentum_point(
                &mut cli,
                &format!("{name}_generated"),
                Evaluation::Production,
                &generated_ring.point,
            )?;
            let scale = ring.re.abs().max(ring.im.abs());
            assert!(
                scale > 0.0
                    && evaluation_comparison::complex_distance(ring, diagram) <= 1e-12 * scale,
                "{name} differs from GL{} at muB={mu_b}: {ring} vs {diagram}",
                generated_ring.diagram
            );
        }
    }
    integration_tests::clean_test(state_path);
    Ok(())
}

/// Loop momenta [q, p, k] = [e0, e2, e4] with a soft static q, |q| = 0.043,
/// against the hard p and k of [`RING_POINT`].
const SOFT_Q_POINT: [f64; 9] = [0.031, -0.012, 0.027, 0.71, 0.23, 0.43, -0.37, 0.61, 0.29];

/// The static-theory mass `general.m_ir` of the zero-mode subtraction is a
/// setting of its own, independent of `general.m_uv`. The tadpole ring
/// without vacuum subtraction, whose whole-graph U is the UV counterterm of
/// the static theory of its union host, starts at the run card's
/// `m_ir = m_uv` and evaluates differently in arb at `m_ir = 0.3`; moving
/// `m_uv` then keeps `m_ir`, and restoring both masses restores the run
/// card's value exactly.
#[test]
fn thermal_ring_static_theory_mass_is_independent_of_m_uv() -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    let name = TADPOLE_RING.name;
    let (mut cli, state_path) = TADPOLE_RING.fresh_cli(&format!("{name}_static_theory_mass"))?;
    integration_tests::run_commands(
        &mut cli,
        VacuumSubtraction::Off.subtraction_commands(Counterterms::Local),
    )?;
    let m_uv = thermal_rings::assert_run_card_static_theory_mass(&cli, name)?;
    let value = |cli: &mut CLIState| {
        thermal_rings::evaluate_momentum_point(cli, name, Evaluation::Arb, &SOFT_Q_POINT)
    };
    let card = value(&mut cli)?;
    let m_ir = 0.3;
    thermal_rings::set_static_theory_mass(&mut cli, name, m_ir)?;
    let apart = value(&mut cli)?;
    let change =
        evaluation_comparison::complex_distance(apart, card) / card.re.abs().max(card.im.abs());
    println!("{name}: m_ir = m_uv {card}, m_ir = {m_ir} {apart}, relative change {change:e}");
    assert!(
        change > MIN_STATIC_THEORY_MASS_CHANGE,
        "{name} at m_ir = {m_ir} must differ from m_ir = m_uv = {m_uv} at a soft-q point: \
         {apart} vs {card}"
    );
    let set_m_uv = |cli: &mut CLIState, m_uv: f64| {
        cli.run_command(&format!(
            "set process -p {name} -i {name} kv general.m_uv={m_uv}"
        ))
    };
    set_m_uv(&mut cli, 2.0 * m_uv)?;
    let kept = thermal_rings::runtime_settings(&cli, name)?.general.m_ir;
    assert_eq!(kept, m_ir, "{name}: moving m_uv must keep m_ir");
    set_m_uv(&mut cli, m_uv)?;
    thermal_rings::set_static_theory_mass(&mut cli, name, m_uv)?;
    assert_eq!(
        value(&mut cli)?,
        card,
        "{name} with both masses restored must evaluate as the run card"
    );
    integration_tests::clean_test(state_path);
    Ok(())
}

fn assert_bare_power_counting(
    bare: &BareRing,
    vacuum_subtraction: VacuumSubtraction,
) -> Result<()> {
    let _serial = integration_tests::SERIAL_TEST_LOCK.lock();
    let ring = bare.ring;
    let name = ring.name;
    let (mut cli, state_path) = ring.fresh_cli(&format!(
        "{name}_bare_power_counting_{}",
        vacuum_subtraction.state_suffix()
    ))?;
    integration_tests::run_commands(&mut cli, vacuum_subtraction.bare_commands())?;
    assert!(!cli.cli_settings.global.generation.uv.subtract_uv);
    assert_eq!(
        cli.cli_settings.global.generation.medium.vacuum_subtraction,
        vacuum_subtraction.enabled()
    );
    thermal_rings::assert_ring_structure(&cli, ring)?;

    for mu_b in [0.0, 3.0] {
        cli.run_command(&format!("set model muB={mu_b}"))?;
        assert_ir_power_counting(&mut cli, bare, mu_b)?;
        let verdict = thermal_rings::uv_verdict(&mut cli, name, UvScan::ONLY_DIVERGENT)?;
        assert_eq!(
            verdict.total, ring.uv_divergent_limits,
            "{name} at muB={mu_b}: {verdict}"
        );
        // Without UV counterterms every superficially divergent limit
        // diverges, except the all-hard one that vacuum subtraction removes.
        assert_eq!(
            verdict.total - verdict.failed,
            convergent_bare_uv_limits(vacuum_subtraction),
            "the bare {name} {vacuum_subtraction} changed its UV divergences at muB={mu_b}: \
             {verdict}"
        );
    }

    integration_tests::clean_test(state_path);
    Ok(())
}

fn assert_ir_power_counting(cli: &mut CLIState, bare: &BareRing, mu_b: f64) -> Result<()> {
    let name = bare.ring.name;
    let limits = bare
        .limits
        .iter()
        .map(|(limit, _)| *limit)
        .collect::<Vec<_>>();
    // One call profiles every limit: the sampled rays depend on the whole
    // select string, so they differ from the one-call-per-limit run cards and
    // from any subtracted run. That is harmless for these bare-only scalings,
    // which are compared with no other ray.
    let profiles = thermal_rings::profile_ir_limits(cli, name, &limits, 10)?;
    for ((limit, expected), profile) in bare.limits.iter().zip(&profiles) {
        // The profiler's pass flag uses an integrability threshold; joint
        // limits at the logarithmic boundary are asserted by their power.
        match *expected {
            Power(power) => assert!(
                (profile.scaling - power).abs() <= SCALING_TOLERANCE
                    && profile.power_law_fit.r_squared() >= MIN_R_SQUARED,
                "bare {name} {limit} at muB={mu_b} must scale as {power:+}: {profile}"
            ),
            Collinear => assert!(
                profile.scaling > MIN_COLLINEAR_SCALING,
                "bare {name} {limit} at muB={mu_b} must stay integrable: {profile}"
            ),
        }
    }
    Ok(())
}
