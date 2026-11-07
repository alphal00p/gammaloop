//! Shared helpers of the Forcer battery (M2-V5, Tiers B and C): the RustRed p
//! backend against Forcer, through vakint's public API.
//!
//! Include it after its two sibling support files, whose items it uses:
//!
//! ```text
//! #[path = "support/p_integrals.rs"] mod p_integrals;
//! #[path = "support/p_rustred.rs"] mod p_rustred;
//! #[path = "support/p_backends.rs"] mod p_backends;
//! ```
//!
//! Agreement means, for one input and a Laurent depth `terms` (G-scheme, 40
//! digits):
//!
//! 1. `Retain`: RustRed's exact value equals Forcer's in the tag basis
//!    `Π Γ(1+mε)^e × Π PMaster` (`vakint::p_integral::compare`);
//! 2. `Symbolic`: the expanded difference is zero;
//! 3. `Numerical`: both values at `p(i) = (0,0,0,1)` (so `Q² = 1`) and
//!    `mursq = 1`, 40 digits, agree within `1e-30`.
//!
//! The RustRed side is [`p_battery_engine`] with RustRed alone
//! ([`rustred_p_settings`]); the Forcer side is a live Forcer (FORM 5, Tier
//! B) or the frozen fixture outputs (Tier C).
//!
//! STAGE 1 (until M2-V4 ships the built-in 1L-3L library): the RustRed engine
//! holds RustRed's one-loop reference library, [`P_BATTERY_LOOPS`] is `1..=1`
//! and the two- and three-loop tests are ignored. STAGE 2: `Vakint::new()`
//! with its built-in library, `1..=3`.

#![allow(dead_code)]

use std::collections::HashMap;
use std::ops::RangeInclusive;
use std::sync::LazyLock;

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::float::Float;
use symbolica::id::Replacement;
use vakint::p_integral::compare::tag_basis_difference;
use vakint::symbols::S;
use vakint::{
    EvaluationMethod, EvaluationOrder, ForcerOptions, LoopNormalizationFactor, Momentum,
    NumericalEvaluationResult, RustRedEvaluationOptions, RustRedMasterMode, Vakint, VakintError,
    VakintSettings, vakint_symbol,
};

use super::p_integrals::{VAKINT, forcer_form, parse};
pub use super::p_rustred::rustred_p_settings;

/// The loop counts [`p_battery_engine`] reduces with RustRed (STAGE 1: the
/// one-loop reference library).
pub const P_BATTERY_LOOPS: RangeInclusive<usize> = 1..=1;

/// The relative tolerance of the numerical agreement at 40 digits.
pub const NUMERICAL_TOLERANCE: f64 = 1e-30;

static P_BATTERY_ENGINE: LazyLock<Vakint> = LazyLock::new(|| {
    super::p_rustred::rustred_p_engine(
        rustred::propagator::testing::one_loop_reference_library().unwrap(),
    )
});

/// The engine of the RustRed side: RustRed's one-loop reference library
/// injected (STAGE 1); STAGE 2 uses `Vakint::new()` and its built-in library.
pub fn p_battery_engine() -> &'static Vakint {
    &P_BATTERY_ENGINE
}

/// The Laurent depth through `ε^(L+2)` for `loops` loops (the depth of the
/// battery): `2L + 3` terms from `ε^-L`.
pub fn battery_terms(loops: usize) -> i64 {
    2 * loops as i64 + 3
}

/// The Forcer side of [`rustred_p_settings`]: the same settings with FORM at
/// `form` and the order Forcer then RustRed, so that the p-integrals go to
/// Forcer and the vacuum integrals of a mixed input to RustRed, as on the
/// RustRed side.
pub fn forcer_side_settings(
    form: &str,
    mode: RustRedMasterMode,
    normalization: LoopNormalizationFactor,
    terms: i64,
) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![
            EvaluationMethod::Forcer(ForcerOptions { master_mode: mode }),
            EvaluationMethod::RustRed(RustRedEvaluationOptions { master_mode: mode }),
        ]),
        form_exe_path: form.into(),
        ..rustred_p_settings(mode, normalization, terms)
    }
}

/// The modes of an agreement, in the order they are checked.
pub const MODES: [RustRedMasterMode; 3] = [
    RustRedMasterMode::Retain,
    RustRedMasterMode::Symbolic,
    RustRedMasterMode::Numerical,
];

/// The value of `input` on the RustRed side in `mode`, G-scheme, `terms`.
pub fn rustred_value(input: &Atom, mode: RustRedMasterMode, terms: i64) -> Result<Atom, String> {
    let settings = rustred_p_settings(mode, LoopNormalizationFactor::GScheme, terms);
    p_battery_engine()
        .evaluate(&settings, input.as_view())
        .map_err(|error| format!("RustRed ({mode}): {error}"))
}

/// The value of `input` with a live Forcer at `form` in `mode`, G-scheme,
/// `terms`.
pub fn live_forcer_value(
    form: &str,
    input: &Atom,
    mode: RustRedMasterMode,
    terms: i64,
) -> Result<Atom, String> {
    let settings = forcer_side_settings(form, mode, LoopNormalizationFactor::GScheme, terms);
    VAKINT
        .evaluate(&settings, input.as_view())
        .map_err(|error| format!("Forcer ({mode}): {error}"))
}

/// `value` at `p(i) = (0,0,0,1)` for the heads `p(1)..p(4)`, `mursq = 1` and
/// the vacuum mass of the mixed inputs `muvsq = 1`, at the 40 digits of the
/// battery settings.
pub fn at_unit_scale(value: &Atom) -> Result<NumericalEvaluationResult, VakintError> {
    let settings = rustred_p_settings(
        RustRedMasterMode::Numerical,
        LoopNormalizationFactor::GScheme,
        1,
    );
    let bits = settings.get_binary_precision();
    let float = |value: i64| Float::with_val(bits, value);
    let params = HashMap::<String, Float, ahash::RandomState>::from_iter([
        (settings.mu_r_sq_symbol.clone(), float(1)),
        ("muvsq".to_owned(), float(1)),
    ]);
    let externals = (1..=4)
        .map(|head| {
            (
                head,
                Momentum::Real((float(0), float(0), float(0), float(1))),
            )
        })
        .collect::<HashMap<usize, Momentum, ahash::RandomState>>();
    p_battery_engine()
        .numerical_evaluation(
            &settings,
            value.as_view(),
            &params,
            &Default::default(),
            Some(&externals),
        )
        .map(|(value, _)| value)
}

/// Whether `rustred` and `forcer`, two values of one input in `mode`, agree
/// (module documentation); in `Numerical` mode `forcer` may be any exact or
/// numerical form of the value. `Err` describes the first difference.
pub fn values_agree(rustred: &Atom, forcer: &Atom, mode: RustRedMasterMode) -> Result<(), String> {
    match mode {
        RustRedMasterMode::Retain => {
            let epsilon = vakint_symbol!(&VakintSettings::default().epsilon_symbol);
            let difference = tag_basis_difference(rustred, forcer, epsilon)
                .map_err(|error| format!("tag basis: {error}"))?;
            if difference.is_empty() {
                Ok(())
            } else {
                Err(format!("tag-basis difference {difference:?}"))
            }
        }
        RustRedMasterMode::Symbolic => {
            let difference = (rustred - forcer).expand();
            if difference.is_zero() {
                Ok(())
            } else {
                Err(format!("symbolic difference {difference}"))
            }
        }
        RustRedMasterMode::Numerical => {
            let numbers = |value: &Atom| {
                at_unit_scale(value)
                    .map_err(|error| format!("numerical evaluation of {value}: {error}"))
            };
            let (rustred, forcer) = (numbers(rustred)?, numbers(forcer)?);
            let (equal, diagnostic) =
                rustred.does_approx_match(&forcer, None, NUMERICAL_TOLERANCE, 5.0);
            if equal { Ok(()) } else { Err(diagnostic) }
        }
    }
}

/// The agreement of `input` through `terms` Laurent terms, the Forcer side
/// given by `forcer(mode)`: every mode of [`MODES`], and a nonzero value (an
/// agreement of zeros tests nothing). `Err` names the mode and the difference.
pub fn agreement(
    input: &Atom,
    terms: i64,
    forcer: impl Fn(RustRedMasterMode) -> Result<Atom, String>,
) -> Result<(), String> {
    for mode in MODES {
        let expected = forcer(mode)?;
        if expected.is_zero() {
            return Err(format!(
                "Forcer gives 0 ({mode}): the agreement tests nothing"
            ));
        }
        let actual = rustred_value(input, mode, terms)?;
        values_agree(&actual, &expected, mode).map_err(|difference| {
            format!("{mode}: RustRed {actual}\nForcer {expected}\n{difference}")
        })?;
    }
    Ok(())
}

/// `input` (a vakint expression with p-integrals, vacuum terms allowed)
/// evaluates alike with RustRed and with a live Forcer through `terms`
/// Laurent terms (module documentation). Skipped without FORM 5.
pub fn assert_p_backends_agree(input: &str, terms: i64) {
    assert_p_backends_agree_atom(&parse(input), terms);
}

/// [`assert_p_backends_agree`] for a parsed input.
pub fn assert_p_backends_agree_atom(input: &Atom, terms: i64) {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the RustRed/Forcer agreement of {input}");
        return;
    };
    agreement(input, terms, |mode| {
        live_forcer_value(&form, input, mode, terms)
    })
    .unwrap_or_else(|difference| panic!("{input}:\n{difference}"));
}

/// The Laurent depth of [`assert_p_backends_agree_routed`].
pub const ROUTED_TERMS: i64 = 5;

/// Routing invariance under RustRed, and agreement: `input` (a long form
/// `topo(prop(...))` with numerator) agrees with a live Forcer
/// ([`assert_p_backends_agree`], [`ROUTED_TERMS`] terms), and for every seed
/// a random unimodular change of loop momenta `k -> M k + c·p` with random
/// line reversals ([`rerouted`]) evaluates with RustRed to exactly the same
/// `Retain` and `Symbolic` values. The reroutings need no FORM.
pub fn assert_p_backends_agree_routed(input: &str, seeds: std::ops::Range<u64>) {
    assert_p_backends_agree_routed_atom(&parse(input), seeds);
}

/// [`assert_p_backends_agree_routed`] for a parsed input.
pub fn assert_p_backends_agree_routed_atom(original: &Atom, seeds: std::ops::Range<u64>) {
    let input = original;
    assert_p_backends_agree_atom(original, ROUTED_TERMS);
    for mode in [RustRedMasterMode::Retain, RustRedMasterMode::Symbolic] {
        let expected = rustred_value(original, mode, ROUTED_TERMS)
            .unwrap_or_else(|error| panic!("{input}: {error}"));
        assert!(!expected.is_zero(), "{input} ({mode}) vanishes");
        for seed in seeds.clone() {
            let routed = rerouted(original, seed);
            let found = rustred_value(&routed, mode, ROUTED_TERMS)
                .unwrap_or_else(|error| panic!("seed {seed}: {routed}: {error}"));
            values_agree(&found, &expected, mode).unwrap_or_else(|difference| {
                panic!("seed {seed} ({mode}): {routed}\n{difference}")
            });
        }
    }
}

/// A deterministic xorshift generator for reproducible reroutings.
struct Random(u64);

impl Random {
    fn next(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

/// The external heads `p(i)` and the number of loop momenta `k(j)` of `input`.
fn heads_and_loops(input: &Atom) -> (Vec<Atom>, usize) {
    let indices = |pattern: &str| {
        input
            .pattern_match(&parse(pattern).to_pattern(), None, None)
            .map(|matched| matched[&parse("i_").get_symbol().unwrap()].clone())
            .collect::<std::collections::BTreeSet<_>>()
    };
    let heads = indices("p(i_)")
        .into_iter()
        .map(|index| S.p.call((index,)))
        .collect();
    (heads, indices("k(i_)").len())
}

/// `input` (a long form with `L` loop momenta `k(1..=L)` and one external
/// head) rerouted by `seed`: `k(i) -> Σ_j M_ij k(j) + c_i·p` for a random
/// unimodular integer matrix `M` and `c_i` in `{-1, 0, 1}`, numerator and
/// lines alike, and every line reversed with probability 1/2
/// (`prop(n, edge(a,b), q, m, x) -> prop(n, edge(b,a), -q, m, x)`). The
/// integral is unchanged.
pub fn rerouted(input: &Atom, seed: u64) -> Atom {
    let (heads, loops) = heads_and_loops(input);
    assert_eq!(heads.len(), 1, "{input}: one external head");
    let p = &heads[0];
    let mut random = Random(0x9E37_79B9_7F4A_7C15 ^ (seed + 1).wrapping_mul(0x2545_F491_4F6C_DD1D));
    let mut matrix = (0..loops)
        .map(|i| (0..loops).map(|j| i64::from(i == j)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for _ in 0..3 * loops {
        let (i, j) = (random.next(loops), random.next(loops));
        if i != j {
            let sign = if random.next(2) == 0 { 1 } else { -1 };
            let row = matrix[j].clone();
            matrix[i]
                .iter_mut()
                .zip(row)
                .for_each(|(a, b)| *a += sign * b);
        } else {
            matrix[i].iter_mut().for_each(|a| *a = -*a);
        }
    }
    let k = |axis: usize| S.k.call((axis as i64,));
    let rules = (0..loops)
        .map(|i| {
            let shift = random.next(3) as i64 - 1;
            let image = (0..loops).fold(Atom::num(shift) * p, |sum, j| {
                sum + Atom::num(matrix[i][j]) * k(j + 1)
            });
            Replacement::new(k(i + 1).to_pattern(), image.to_pattern())
        })
        .collect::<Vec<_>>();
    let substituted = input.replace_multiple(&rules);
    let flips = (0..64).map(|_| random.next(2) == 1).collect::<Vec<_>>();
    substituted.replace_map(|view, _, out| {
        let AtomView::Fun(prop) = view else {
            return;
        };
        if prop.get_symbol() != S.prop || prop.get_nargs() != 5 {
            return;
        }
        let arguments = prop.iter().collect::<Vec<_>>();
        let id = i64::try_from(arguments[0]).unwrap();
        let AtomView::Fun(edge) = arguments[1] else {
            panic!("{view} has no edge(a,b)");
        };
        let ends = edge.iter().map(|end| end.to_owned()).collect::<Vec<_>>();
        let momentum = arguments[2].to_owned().expand();
        let (edge, momentum) = if flips[id as usize % flips.len()] {
            (S.edge.call((&ends[1], &ends[0])), -momentum)
        } else {
            (S.edge.call((&ends[0], &ends[1])), momentum)
        };
        **out = S.prop.call((
            arguments[0].to_owned(),
            edge,
            momentum.expand(),
            arguments[3].to_owned(),
            arguments[4].to_owned(),
        ));
    })
}

/// `lhs` and `rhs` evaluate with RustRed alone to exactly the same nonzero
/// value (Symbolic, G-scheme, through `ε^(5-L)` like the Forcer identities of
/// `p_integral_forcer_tests.rs`): an analytic identity between p-integrals,
/// FORM-free (M1 review, missing test 1, under RustRed).
pub fn assert_rustred_p_identity(lhs: &str, rhs: &str) {
    let value = |input: &str| {
        rustred_value(&parse(input), RustRedMasterMode::Symbolic, 6)
            .unwrap_or_else(|error| panic!("{input}: {error}"))
    };
    let (left, right) = (value(lhs), value(rhs));
    assert!(
        !left.is_zero(),
        "{lhs} vanishes at this truncation: the identity tests nothing"
    );
    let difference = (&left - &right).expand();
    assert!(difference.is_zero(), "{lhs} - ({rhs}) = {difference}");
}

/// The `Oep` guard under RustRed (M1 review, missing test 5): `input`
/// evaluates through `passing` Laurent terms with RustRed and with a live
/// Forcer to the same value without an `Oep` tail, and one term more is a
/// `MissingMasterOrders` error on both sides (the boundary is the masters'
/// table, shared by both backends). The Forcer side is skipped without FORM
/// 5; the RustRed side always runs.
pub fn assert_p_backends_oep_boundary(input: &str, passing: i64) {
    let atom = parse(input);
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let value = rustred_value(&atom, mode, passing)
            .unwrap_or_else(|error| panic!("{input} ({mode}, {passing} terms): {error}"));
        assert!(
            !value.to_string().contains("Oep"),
            "{input} ({mode}): {value}"
        );
        let settings = rustred_p_settings(mode, LoopNormalizationFactor::GScheme, passing + 1);
        match p_battery_engine().evaluate(&settings, atom.as_view()) {
            Err(VakintError::MissingMasterOrders(_)) => {}
            other => panic!(
                "{input} ({mode}) through {} terms: expected MissingMasterOrders, got {other:?}",
                passing + 1
            ),
        }
        let Some(form) = forcer_form() else {
            eprintln!("FORM 5 is unavailable: skipped the Forcer side of the Oep boundary");
            continue;
        };
        let forcer = live_forcer_value(&form, &atom, mode, passing)
            .unwrap_or_else(|error| panic!("{input}: {error}"));
        values_agree(&value, &forcer, mode)
            .unwrap_or_else(|difference| panic!("{input} ({mode}): {difference}"));
        let settings =
            forcer_side_settings(&form, mode, LoopNormalizationFactor::GScheme, passing + 1);
        match VAKINT.evaluate(&settings, atom.as_view()) {
            Err(VakintError::MissingMasterOrders(_)) => {}
            other => panic!(
                "{input} ({mode}) with Forcer through {} terms: expected a missing master order, got {other:?}",
                passing + 1
            ),
        }
    }
}
