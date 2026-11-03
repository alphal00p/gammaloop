//! The `Oep` truncation guard on real integrals, backend by backend.
//!
//! Every backend substitutes a master integral by its known Laurent
//! coefficients plus an explicit unknown tail `ε^k·Oep(k, master)`. A reduction
//! coefficient with a `1/ε^p` pole moves that tail `p` orders down; here the
//! pole is written into the input as the prefactor `ε^-p`. Each case requests
//!
//! 1. the deepest order its masters support with the pole: this succeeds and,
//!    for exact results, equals `ε^-p` times the pole-free result `p` terms
//!    deeper (the guard changes no correct result);
//! 2. one order more: this must fail with `VakintError::MissingMasterOrders`,
//!    naming the master, its first missing Laurent order (needed exactly one
//!    order too early) and the requested order, instead of silently treating
//!    the unknown coefficient as zero.

mod test_utils;

use symbolica::atom::{Atom, AtomCore};
use test_utils::rustred_support::rustred_only;
use vakint::{
    EvaluationOrder, RustRedMasterMode, VakintError, VakintSettings, vakint_parse, vakint_symbol,
};

/// RustRed never runs FORM; a FORM call would fail on this path.
const NO_FORM: &str = "/vakint-truncation-guard-tests-never-run-form";

const ONE_LOOP_TADPOLE: &str = "topo(prop(1,edge(1,1),k(1),muvsq,1))";
const THREE_LOOP_MERCEDES: &str = "topo(
      prop(1,edge(1,2),k(1),muvsq,1)
     *prop(2,edge(2,3),k(2),muvsq,1)
     *prop(3,edge(3,1),k(3),muvsq,1)
     *prop(4,edge(1,4),k(3)-k(1),muvsq,1)
     *prop(5,edge(2,4),k(1)-k(2),muvsq,1)
     *prop(6,edge(3,4),k(2)-k(3),muvsq,1))";
/// Four one-loop tadpoles at one vertex.
const FOUR_LOOP_CLOVER: &str = "topo(
      prop(1,edge(1,1),k(1),muvsq,1)
     *prop(2,edge(1,1),k(2),muvsq,1)
     *prop(3,edge(1,1),k(3),muvsq,1)
     *prop(4,edge(1,1),k(4),muvsq,1))";

/// One integral, one backend and the depth its masters support.
struct GuardCase {
    name: &'static str,
    settings: VakintSettings,
    loops: i64,
    integral: &'static str,
    /// The order `p` of the prefactor `ε^-p`.
    pole: i64,
    /// The deepest number of Laurent terms (from `ε^-loops`) that succeeds
    /// with the pole; one more term needs one unknown master order.
    deepest_terms: i64,
    /// A fragment of the reported master name.
    master: &'static str,
    /// The first Laurent order of that master that is not available.
    first_missing_order: i64,
    /// Whether the backend's result is exact, so that it can be compared
    /// with the pole-free evaluation.
    exact: bool,
}

impl GuardCase {
    fn evaluate(&self, terms: i64, pole: i64) -> Result<Atom, VakintError> {
        let vakint = test_utils::get_vakint(VakintSettings {
            number_of_terms_in_epsilon_expansion: terms,
            ..self.settings.clone()
        });
        let input = vakint_parse!(format!("ε^(-{pole})*{}", self.integral).as_str())
            .expect("guard-case input parses");
        vakint.evaluate(input.as_view())
    }
}

fn check(case: GuardCase) {
    test_utils::run_multi_lane_acceptance(move || {
        if case.loops >= 4 {
            test_utils::ensure_external_data_loaded();
        }
        let within = case
            .evaluate(case.deepest_terms, case.pole)
            .unwrap_or_else(|error| panic!("{}: within the master depth: {error}", case.name));
        if case.exact {
            let pole_free = case
                .evaluate(case.deepest_terms + case.pole, 0)
                .unwrap_or_else(|error| panic!("{}: pole-free evaluation: {error}", case.name));
            let epsilon = Atom::var(vakint_symbol!("ε"));
            let difference = (within - epsilon.pow(Atom::num(-case.pole)) * pole_free).expand();
            assert!(
                difference.is_zero(),
                "{}: the guard changed a correct result: {difference}",
                case.name
            );
        }

        let missing = test_utils::expect_missing_master_orders(
            case.evaluate(case.deepest_terms + 1, case.pole),
        );
        eprintln!("{}:\n{missing}", case.name);
        // terms + 1 Laurent coefficients from ε^-loops end at ε^(terms - loops).
        assert_eq!(
            missing.requested_order,
            Some(case.deepest_terms - case.loops),
            "{}: {missing}",
            case.name
        );
        let master = missing
            .masters
            .iter()
            .find(|master| master.master.contains(case.master))
            .unwrap_or_else(|| panic!("{}: {} not reported: {missing}", case.name, case.master));
        assert_eq!(
            (master.first_missing_order, master.last_needed_order),
            (case.first_missing_order, Some(case.first_missing_order)),
            "{}: {missing}",
            case.name
        );
    });
}

/// RustRed alone, with masters substituted in `mode`.
fn rustred(mode: RustRedMasterMode) -> VakintSettings {
    rustred_only(mode, NO_FORM)
}

#[test]
fn rustred_one_loop_tadpole_symbolic() {
    // The vendored tadpole is known through ε^4.
    check(GuardCase {
        name: "RustRed 1L tadpole, symbolic masters",
        settings: rustred(RustRedMasterMode::Symbolic),
        loops: 1,
        integral: ONE_LOOP_TADPOLE,
        pole: 1,
        deepest_terms: 5,
        master: "1-loop",
        first_missing_order: 5,
        exact: true,
    });
}

#[test]
fn rustred_one_loop_tadpole_numerical() {
    check(GuardCase {
        name: "RustRed 1L tadpole, numerical masters (default engine)",
        settings: rustred(RustRedMasterMode::Numerical),
        loops: 1,
        integral: ONE_LOOP_TADPOLE,
        pole: 2,
        deepest_terms: 4,
        master: "1-loop",
        first_missing_order: 5,
        exact: false,
    });
}

#[test]
fn rustred_three_loop_mercedes() {
    // The vendored Mercedes master is known through ε^2.
    check(GuardCase {
        name: "RustRed 3L Mercedes, symbolic masters",
        settings: rustred(RustRedMasterMode::Symbolic),
        loops: 3,
        integral: THREE_LOOP_MERCEDES,
        pole: 1,
        deepest_terms: 5,
        master: "3-loop",
        first_missing_order: 3,
        exact: true,
    });
}

#[test]
fn rustred_four_loop_clover() {
    // The clover is the four-loop master PR1 (FMFT's tadpole^4), shipped in
    // the external data through ε^2.
    check(GuardCase {
        name: "RustRed 4L clover (master PR1), numerical masters",
        settings: rustred(RustRedMasterMode::Numerical),
        loops: 4,
        integral: FOUR_LOOP_CLOVER,
        pole: 1,
        deepest_terms: 6,
        master: "4-loop",
        first_missing_order: 3,
        exact: false,
    });
}

#[test]
fn alphaloop_one_loop_tadpole() {
    // alphaLoop's `uvid(1,1) = ε·tadpole` ends in `alarmt1m3·ε^4`.
    check(GuardCase {
        name: "alphaLoop 1L tadpole",
        settings: VakintSettings {
            evaluation_order: EvaluationOrder::alphaloop_only(),
            ..VakintSettings::default()
        },
        loops: 1,
        integral: ONE_LOOP_TADPOLE,
        pole: 1,
        deepest_terms: 3,
        master: "uvid(1,1)",
        first_missing_order: 4,
        exact: false,
    });
}

#[test]
fn matad_three_loop_mercedes() {
    // MATAD's `miD6` ends in `ε^3·Oep(3,D6)`.
    check(GuardCase {
        name: "MATAD 3L Mercedes",
        settings: VakintSettings {
            evaluation_order: EvaluationOrder::matad_only(None),
            ..VakintSettings::default()
        },
        loops: 3,
        integral: THREE_LOOP_MERCEDES,
        pole: 2,
        deepest_terms: 4,
        master: "D6",
        first_missing_order: 3,
        exact: false,
    });
}
