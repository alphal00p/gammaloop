//! The Forcer backend on canonical inputs.
//!
//! The batteries in `tests/fixtures/forcer/p{L}l_battery.json` freeze the
//! exact Forcer output of each case (FORM 5.0.0, Forcer 68a30bf). Without
//! FORM 5 they test the input generation, the output parser and the shared
//! finalizer. With FORM 5 (`VAKINT_TEST_FORM5=/path/to/form`, or the default
//! `form` when it is version 5) every case also runs live and must reproduce
//! the frozen output; `VAKINT_FORCER_FIXTURES=write` regenerates the files.

use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};
use symbolica::transcendental::gamma;

use super::super::finalize::{finalize, pg_value, tests as finalize_tests};
use super::super::{P_FAMILIES, PFamilyData, ReducedPropagatorValue, pg_symbol};
use super::*;
use crate::rustred_evaluation::RustRedMasterMode;
use crate::utils::vakint_macros::vk_parse;
use crate::{EvaluationOrder, LoopNormalizationFactor, VakintError, VakintSettings};

#[derive(Debug, Deserialize, Serialize)]
struct Battery {
    forcer_commit: String,
    form_version: String,
    cases: Vec<Case>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct Case {
    name: String,
    family: String,
    powers: Vec<i64>,
    numerator: String,
    external: String,
    /// The raw `%E` output of each Forcer local, in `ForcerInput::parts` order.
    outputs: Vec<String>,
}

fn family(name: &str) -> &'static PFamilyData {
    P_FAMILIES
        .iter()
        .find(|family| family.name == name)
        .unwrap_or_else(|| panic!("unknown family {name}"))
}

fn parse(text: &str) -> Atom {
    Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap_or_else(|error| panic!("{text}: {error}"))
}

/// FORM 5 for live Forcer runs: `VAKINT_TEST_FORM5`, else the default FORM
/// when it is version 5.
fn forcer_test_form() -> Option<String> {
    std::env::var("VAKINT_TEST_FORM5").ok().or_else(|| {
        let form = crate::default_form_exe_path();
        forcer_form_version(&form).is_ok().then_some(form)
    })
}

fn settings(form: &str) -> VakintSettings {
    VakintSettings {
        form_exe_path: form.into(),
        n_cores: 1,
        ..finalize_tests::settings(LoopNormalizationFactor::pySecDec, 4)
    }
}

fn epsilon() -> Symbol {
    vk_symbol!(&VakintSettings::default().epsilon_symbol)
}

impl Case {
    fn input(&self) -> Option<ForcerInput> {
        forcer_input(
            family(&self.family),
            &self.powers,
            parse(&self.numerator).as_view(),
            &parse(&self.external),
        )
        .unwrap_or_else(|error| panic!("{}: {error}", self.name))
    }

    /// The reduced values of the frozen Forcer output.
    fn frozen(&self) -> Vec<ReducedPropagatorValue> {
        match self.input() {
            None => {
                assert!(self.outputs.is_empty(), "{}: scaleless", self.name);
                Vec::new()
            }
            Some(input) => {
                assert_eq!(input.parts.len(), self.outputs.len(), "{}", self.name);
                combine_parts(&input.parts, &self.outputs, epsilon())
                    .unwrap_or_else(|error| panic!("{}: {error}", self.name))
            }
        }
    }
}

fn fixture_path(file: &str) -> String {
    format!("{}/tests/fixtures/forcer/{file}", env!("CARGO_MANIFEST_DIR"))
}

/// Every case of a battery: frozen values, and the same values from a live
/// Forcer run when FORM 5 is available.
fn forcer_battery(file: &str, frozen: &str) -> Vec<(Case, Vec<ReducedPropagatorValue>)> {
    let mut battery: Battery = serde_json::from_str(frozen).unwrap();
    assert_eq!(battery.forcer_commit, FORCER_UPSTREAM_COMMIT);
    let write = std::env::var("VAKINT_FORCER_FIXTURES").is_ok_and(|mode| mode == "write");
    let live = forcer_test_form();
    if write {
        let form = live.as_deref().expect("regenerating Forcer fixtures needs FORM 5");
        battery.form_version = forcer_form_version(form).unwrap();
        for case in &mut battery.cases {
            case.outputs = match case.input() {
                None => Vec::new(),
                Some(input) => {
                    let names = input.parts.iter().map(|p| p.name.as_str()).collect::<Vec<_>>();
                    run_forcer(&settings(form), &input.program, &names).unwrap()
                }
            };
        }
        std::fs::write(
            fixture_path(file),
            serde_json::to_string_pretty(&battery).unwrap() + "\n",
        )
        .unwrap();
    }
    battery
        .cases
        .into_iter()
        .map(|case| {
            let frozen = case.frozen();
            if let Some(form) = &live {
                let reduced = reduce_p(
                    family(&case.family),
                    &case.powers,
                    parse(&case.numerator).as_view(),
                    &parse(&case.external),
                    &settings(form),
                )
                .unwrap_or_else(|error| panic!("{}: {error}", case.name));
                assert_eq!(reduced, frozen, "{}: live Forcer differs from the fixture", case.name);
            }
            (case, frozen)
        })
        .collect()
}

/// Every value finalizes in each master mode without a missing master order.
fn assert_finalizes(case: &Case, values: &[ReducedPropagatorValue], terms: i64) {
    let loops = family(&case.family).loops;
    let external = parse(&case.external);
    for mode in [RustRedMasterMode::Retain, RustRedMasterMode::Symbolic] {
        let settings = finalize_tests::settings(LoopNormalizationFactor::MSbar, terms);
        let result = finalize(values, loops, &external, &settings, mode, "Forcer")
            .unwrap_or_else(|error| panic!("{} ({mode}): {error}", case.name));
        assert_eq!(result.is_zero(), values.is_empty(), "{} ({mode})", case.name);
    }
}

fn assert_series_zero(expression: Atom, depth: i64, context: &str) {
    let difference = super::super::finalize::gamma_series(&expression, epsilon(), depth).unwrap();
    assert!(difference.is_zero(), "{context}: {difference}");
}

/// The reduced values with every `PG` written in Γ functions.
fn exact(values: &[ReducedPropagatorValue], scale_power: i64) -> Atom {
    let mut total = Atom::Zero;
    for value in values.iter().filter(|value| value.scale_power == scale_power) {
        total += &value.reduced;
    }
    total
        .replace(pg_symbol().call((Atom::var(S.x_), Atom::var(S.y_))).to_pattern())
        .with_map(|matched| {
            let x = get_integer_from_atom(matched.get(S.x_).unwrap().to_atom().as_view()).unwrap();
            let y = get_integer_from_atom(matched.get(S.y_).unwrap().to_atom().as_view()).unwrap();
            pg_value(x, y, epsilon())
        })
}

/// `G(a,b)` at `Q² = 1`, zero for a non-positive integer index.
fn g(a: i64, b: Atom) -> Atom {
    if a <= 0 {
        return Atom::Zero;
    }
    let e = Atom::var(epsilon());
    let half_d = Atom::num(2) - &e;
    let d = Atom::num(4) - Atom::num(2) * &e;
    let a = Atom::num(a);
    gamma().call(&a + &b - &half_d)
        * gamma().call(&half_d - &a)
        * gamma().call(&half_d - &b)
        / (gamma().call(a.clone()) * gamma().call(b.clone()) * gamma().call(d - a - b))
}

fn g_scheme_measure() -> Atom {
    Atom::var(epsilon()) * g(1, Atom::num(1))
}

/// The one-loop closed form of `dot(k(1),p)^n` on `(a, b)`: with
/// `k·q = (k² + q² - (k-q)²)/2` every term is a `G` function.
fn bubble_closed_form(a: i64, b: i64, n: i64) -> Atom {
    let mut total = Atom::Zero;
    let factorial = |m: i64| (1..=m).product::<i64>();
    for i in 0..=n {
        for l in 0..=n - i {
            let j = n - i - l;
            let multinomial = factorial(n) / (factorial(i) * factorial(j) * factorial(l));
            let sign = if l % 2 == 0 { 1 } else { -1 };
            if b - l > 0 {
                total += Atom::num(sign * multinomial) * g(a - i, Atom::num(b - l));
            }
        }
    }
    total / Atom::num(1 << n) / g_scheme_measure()
}

/// `(a, b, n, σ)` of a one-loop case `dot(k(1),p(h))^n` with `external = σ·p(h)`.
fn bubble_case(case: &Case) -> (i64, i64, i64, i64) {
    let n = case
        .numerator
        .split('^')
        .nth(1)
        .map_or(if case.numerator == "1" { 0 } else { 1 }, |n| n.parse().unwrap());
    let sigma = if case.external.starts_with('-') { -1 } else { 1 };
    (case.powers[0], case.powers[1], n, sigma)
}

#[test]
fn forcer_one_loop_battery() {
    let battery = forcer_battery(
        "p1l_battery.json",
        include_str!("../../../tests/fixtures/forcer/p1l_battery.json"),
    );
    for (case, values) in battery {
        let (a, b, n, sigma) = bubble_case(&case);
        if a <= 0 || b <= 0 {
            assert!(values.is_empty(), "{}: scaleless", case.name);
            continue;
        }
        // reduced = (-1)^(a+b+n) σ^n C_G at w = 2 - a - b + n.
        let sign = if (a + b + n) % 2 == 0 { 1 } else { -1 } * if n % 2 == 0 { 1 } else { sigma };
        assert_series_zero(
            exact(&values, 2 - a - b + n) - Atom::num(sign) * bubble_closed_form(a, b, n),
            5,
            &case.name,
        );
        assert_finalizes(&case, &values, 4);
    }
}

#[test]
fn forcer_two_loop_battery() {
    let battery = forcer_battery(
        "p2l_battery.json",
        include_str!("../../../tests/fixtures/forcer/p2l_battery.json"),
    );
    for (case, values) in &battery {
        assert_finalizes(case, values, 4);
    }
    let value = |name: &str| {
        battery
            .iter()
            .find(|(case, _)| case.name == name)
            .unwrap_or_else(|| panic!("missing case {name}"))
            .1
            .clone()
    };
    // The kite: G-scheme value 6ζ3 + O(ε), reduced = -C_G at w = -1.
    let head = parse("p(1)");
    let settings = finalize_tests::settings(LoopNormalizationFactor::GScheme, 3);
    let kite = finalize(
        &value("kite"),
        2,
        &head,
        &settings,
        RustRedMasterMode::Symbolic,
        "Forcer",
    )
    .unwrap();
    let coefficients =
        finalize_tests::coefficients_at_unit_scale(&kite, &settings, &head, -2, 0);
    assert_eq!(
        coefficients,
        [
            Atom::Zero,
            Atom::Zero,
            (Atom::num(-6) * symbolica::transcendental::zeta().call(3)).expand()
        ]
    );
    // The sunset: C_G = G(1,1)G(1,ε)/(εG(1,1))², reduced = -C_G at w = 1.
    let e = Atom::var(epsilon());
    let sunset = g(1, Atom::num(1)) * g(1, e.clone()) / g_scheme_measure().pow(Atom::num(2));
    assert_series_zero(exact(&value("sunset"), 1) + sunset, 5, "sunset");
    // q reversal: the canonical p(1) written as -p(3).
    assert_eq!(value("kite_kq_rank1"), value("kite_kq_rank1_reversed"));
}

#[test]
fn forcer_three_loop_battery() {
    let battery = forcer_battery(
        "p3l_battery.json",
        include_str!("../../../tests/fixtures/forcer/p3l_battery.json"),
    );
    let head = parse("p(1)");
    let settings = finalize_tests::settings(LoopNormalizationFactor::GScheme, 4);
    for (case, values) in &battery {
        assert_finalizes(case, values, 4);
        // The scalar planar and non-planar tops have the G-scheme value 20ζ5.
        if ["la_top", "no_top"].contains(&case.name.as_str()) {
            let result = finalize(values, 3, &head, &settings, RustRedMasterMode::Symbolic, "Forcer")
                .unwrap();
            let coefficients =
                finalize_tests::coefficients_at_unit_scale(&result, &settings, &head, -3, 0);
            assert_eq!(
                coefficients,
                [
                    Atom::Zero,
                    Atom::Zero,
                    Atom::Zero,
                    (Atom::num(20) * symbolica::transcendental::zeta().call(5)).expand()
                ],
                "{}",
                case.name
            );
        }
    }
}

#[test]
fn forcer_four_loop_battery() {
    let battery = forcer_battery(
        "p4l_battery.json",
        include_str!("../../../tests/fixtures/forcer/p4l_battery.json"),
    );
    for (case, values) in &battery {
        assert_finalizes(case, values, 5);
    }
}

#[test]
fn canonical_kite_becomes_the_forcer_graph() {
    let input = forcer_input(
        family("P2L"),
        &[1, 1, 1, 1, 2],
        parse("dot(k(1),p(3))^2*dot(k(2),k(2))").as_view(),
        &parse("-p(3)"),
    )
    .unwrap()
    .unwrap();
    // Incoming momenta count positively; Q enters at vertex 1, leaves at 2.
    for fragment in [
        "(k1.Q)^2",
        "(k2.k2)^1",
        "/(p1.p1)^1/(p2.p2)^1/(p3.p3)^1/(p4.p4)^1/(p5.p5)^2\
         *vx(p3,-p5,Q)*vx(-p1,p4,-Q)*vx(-p2,-p4,p5)*vx(p1,p2,-p3);",
    ] {
        assert!(input.program.contains(fragment), "{fragment}\n{}", input.program);
    }
    assert!(input.program.contains("id k1 = +p1;\nid k2 = +p2;"), "{}", input.program);
    // w = 2L - Σa + r = 4 - 6 + 3, sign (-1)^(6+3).
    assert_eq!(
        input.parts,
        [ForcerPart {
            name: "F1".into(),
            scale_power: 1,
            sign: -1,
            coefficient: Atom::num(1),
        }]
    );
}

#[test]
fn pinched_loop_lines_are_solved_through_the_remaining_lines() {
    for family in &P_FAMILIES {
        let lines = canonical_lines(family).unwrap();
        for pinched_line in 0..lines.len() {
            let remaining = lines
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != pinched_line)
                .map(|(_, line)| line)
                .collect::<Vec<_>>();
            let solution = loop_momenta_in_lines(&remaining, family.loops)
                .unwrap_or_else(|| panic!("{}: line {pinched_line}", family.name));
            // Σ_e X_je (A_e k + b_e Q) + y_j Q = k_j identically.
            for (j, row) in solution.iter().enumerate() {
                for column in 0..=family.loops {
                    let mut total = if column == family.loops {
                        *row.last().unwrap()
                    } else {
                        Fraction::new(0, 1)
                    };
                    for (x, line) in row.iter().zip(&remaining) {
                        let entry = if column == family.loops {
                            line.external
                        } else {
                            line.loops[column]
                        };
                        total = total.sub(Fraction::new(-(entry as i128), 1).mul(*x));
                    }
                    let expected = Fraction::new((column == j) as i128, 1);
                    assert_eq!(total, expected, "{}: k({}) column {column}", family.name, j + 1);
                }
            }
        }
    }
}

#[test]
fn scaleless_sectors_vanish_without_running_forcer() {
    for (name, powers) in [
        ("P1L", vec![0, 1]),
        ("P1L", vec![2, -1]),
        // k(1) and k(2) pinched: line 4 becomes a self-loop.
        ("P2L", vec![0, 0, 1, 1, 1]),
        // Both q lines pinched: the external vertices coincide.
        ("P2L", vec![1, 1, 1, 0, 0]),
    ] {
        let settings = settings("/forcer-scaleless-test-never-runs-form");
        let values = reduce_p(family(name), &powers, Atom::num(1).as_view(), &parse("p(1)"), &settings)
            .unwrap();
        assert!(values.is_empty(), "{name} {powers:?}");
    }
}

#[test]
fn a_vacuum_part_behind_a_cut_vertex_is_scaleless() {
    let line = |from, to| CanonicalLine {
        from,
        to,
        loops: Vec::new(),
        external: 0,
    };
    let bubble = [line(1, 2), line(2, 1)];
    let with_tadpole = [line(1, 2), line(2, 1), line(2, 3), line(3, 2)];
    for (lines, scaleless) in [(&bubble[..], false), (&with_tadpole[..], true)] {
        let graph = ContractedGraph::new(lines, &vec![false; lines.len()]);
        assert_eq!(graph.is_scaleless(&lines.iter().collect::<Vec<_>>()), scaleless);
    }
}

#[test]
fn tensor_numerators_are_rejected() {
    for numerator in ["dot(k(1),p(2))", "k(1,1)*k(2,2)"] {
        let error = forcer_input(
            family("P2L"),
            &[1, 1, 1, 1, 1],
            parse(numerator).as_view(),
            &parse("p(1)"),
        )
        .unwrap_err();
        assert!(
            matches!(&error, VakintError::InvalidNumerator(message) if message.contains("tensor")),
            "{numerator}: {error}"
        );
    }
}

#[test]
fn forcer_output_maps_to_the_p_master_basis() {
    let e = epsilon();
    let parsed = parse_forcer_output(
        "GschemeConstants(0,0)^2*rat(2*ep - 1,2*ep) + GschemeConstants(0,0)*\n      \
         GschemeConstants(1,0)*rat( - 3*ep + 1,2*ep) + Master(no)*rat(1,1)",
        e,
    )
    .unwrap();
    let expected = parse(
        "PG(0,0)^2*(2*ε-1)/(2*ε)+PG(0,0)*PG(1,0)*(1-3*ε)/(2*ε)+PMaster(no)",
    );
    assert_eq!((parsed - expected).expand(), Atom::Zero);
    for unexpected in ["Master(no)*xi", "ALARMno", "GschemeConstants(x,0)"] {
        assert!(
            matches!(
                parse_forcer_output(unexpected, e),
                Err(VakintError::FormOutputParsingError(..))
            ),
            "{unexpected}"
        );
    }
}

#[test]
fn forcer_requires_form_five() {
    assert_eq!(
        forcer_form_version("/vakint-forcer-test/missing-form"),
        Err("could not be run".to_string())
    );
    let directory = tempfile::tempdir().unwrap();
    let old_form = directory.path().join("form");
    std::fs::write(&old_form, "#!/bin/sh\necho 'FORM 4.3.1 (Jan  1 1980)'\n").unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&old_form, std::fs::Permissions::from_mode(0o755)).unwrap();
        assert_eq!(
            forcer_form_version(&old_form.to_string_lossy()),
            Err("is FORM 4.3.1".to_string())
        );
    }
    let missing = "/vakint-forcer-test/missing-form";
    assert!(
        !EvaluationOrder::available(missing, "python3")
            .0
            .iter()
            .any(|method| matches!(method, crate::EvaluationMethod::Forcer(_)))
    );
    let error = Vakint::new()
        .unwrap()
        .validate_settings(&VakintSettings {
            form_exe_path: missing.into(),
            evaluation_order: EvaluationOrder::forcer_only(),
            ..VakintSettings::default()
        })
        .unwrap_err();
    assert!(
        matches!(&error, VakintError::ForcerUnavailable { form_exe_path, .. } if form_exe_path == missing),
        "{error}"
    );
    let message = error.to_string();
    assert!(
        message.contains("needs FORM >= 5.0.0") && message.contains("EvaluationOrder::rustred_only()"),
        "{message}"
    );
}

#[test]
fn vendored_sources_materialize_once_and_completely() {
    let directory = tempfile::tempdir().unwrap();
    let base = directory.path().to_string_lossy().into_owned();
    let first = materialized_sources(Some(&base)).unwrap();
    assert_eq!(first, materialized_sources(Some(&base)).unwrap());
    assert!(first.join(COMPLETE_MARKER).is_file());
    assert_eq!(
        std::fs::read(first.join("forcer").join("Masters.prc")).unwrap(),
        include_bytes!("../../../form_src/forcer/forcer/Masters.prc")
    );
    assert_eq!(vendored::FORCER_FILES.len(), 21);
}
