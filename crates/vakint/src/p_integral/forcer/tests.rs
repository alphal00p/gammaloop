//! The Forcer backend on canonical inputs.
//!
//! The batteries in `tests/fixtures/forcer/p{L}l_battery.json` freeze the
//! exact Forcer output of each case (FORM 5.0.0, Forcer 68a30bf). Without
//! FORM 5 they test the input generation, the output parser and the shared
//! finalizer. With FORM 5 (`VAKINT_TEST_FORM5=/path/to/form`, or the default
//! `form` when it is version 5) every case also runs live and must reproduce
//! the frozen output; `VAKINT_FORCER_FIXTURES=write` regenerates the files.

use symbolica::atom::{Atom, AtomCore};
use symbolica::transcendental::gamma;

use super::super::finalize::{finalize, pg_value, tests as finalize_tests};
use super::super::parts::{
    CanonicalLine, ContractedGraph, Fraction, canonical_lines, loop_momenta_in_lines,
};
use super::super::{P_FAMILIES, PFamilyData, ReducedPropagatorValue, pg_symbol};
use super::fixtures::{ForcerCase as Case, ForcerFixture as Battery};
use super::*;
use crate::rustred_evaluation::RustRedMasterMode;
use crate::utils::vakint_macros::vk_parse;
use crate::{EvaluationOrder, LoopNormalizationFactor, VakintError, VakintSettings};

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

/// The reduced values of the frozen Forcer output of a case.
fn frozen(case: &Case) -> Vec<ReducedPropagatorValue> {
    case.frozen(epsilon())
        .unwrap_or_else(|error| panic!("{}: {error}", case.name))
}

fn fixture_path(file: &str) -> String {
    format!(
        "{}/tests/fixtures/forcer/{file}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// `VAKINT_FORCER_FIXTURES=write`: regenerate the fixture files (FORM 5).
fn writing_fixtures() -> bool {
    std::env::var("VAKINT_FORCER_FIXTURES").is_ok_and(|mode| mode == "write")
}

/// FORM processes run at a time by the fixture writer
/// (`VAKINT_FORCER_FIXTURES_JOBS`, default 1).
fn writer_jobs() -> usize {
    std::env::var("VAKINT_FORCER_FIXTURES_JOBS")
        .ok()
        .and_then(|jobs| jobs.parse().ok())
        .unwrap_or(1)
}

/// Every case of a battery: frozen values, and the same values from a live
/// Forcer run when FORM 5 is available.
fn forcer_battery(file: &str, frozen_json: &str) -> Vec<(Case, Vec<ReducedPropagatorValue>)> {
    let mut battery = Battery::parse(frozen_json).unwrap();
    let live = forcer_test_form();
    if writing_fixtures() {
        let form = live
            .as_deref()
            .expect("regenerating Forcer fixtures needs FORM 5");
        battery.write_outputs(form, writer_jobs()).unwrap();
        std::fs::write(fixture_path(file), battery.to_json()).unwrap();
    }
    battery
        .cases
        .into_iter()
        .map(|case| {
            let frozen = frozen(&case);
            if let Some(form) = &live {
                let reduced = reduce_p(
                    family(&case.family),
                    &case.powers,
                    parse(&case.numerator).as_view(),
                    &parse(&case.external),
                    &settings(form),
                )
                .unwrap_or_else(|error| panic!("{}: {error}", case.name));
                assert_eq!(
                    reduced, frozen,
                    "{}: live Forcer differs from the fixture",
                    case.name
                );
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
        assert_eq!(
            result.is_zero(),
            values.is_empty(),
            "{} ({mode})",
            case.name
        );
    }
}

fn assert_series_zero(expression: Atom, depth: i64, context: &str) {
    let difference = super::super::finalize::gamma_series(&expression, epsilon(), depth).unwrap();
    assert!(difference.is_zero(), "{context}: {difference}");
}

/// The reduced values with every `PG` written in Γ functions.
fn exact(values: &[ReducedPropagatorValue], scale_power: i64) -> Atom {
    let mut total = Atom::Zero;
    for value in values
        .iter()
        .filter(|value| value.scale_power == scale_power)
    {
        total += &value.reduced;
    }
    total
        .replace(
            pg_symbol()
                .call((Atom::var(S.x_), Atom::var(S.y_)))
                .to_pattern(),
        )
        .with_map(|matched| {
            let x = get_integer_from_atom(matched.get(S.x_).unwrap().to_atom().as_view()).unwrap();
            let y = get_integer_from_atom(matched.get(S.y_).unwrap().to_atom().as_view()).unwrap();
            pg_value(x, y, epsilon())
        })
}

/// `G(a,b)` at `Q² = 1`, zero for a non-positive integer index.
fn g(a: i64, b: Atom) -> Atom {
    bubble(&Atom::num(a), &b)
}

/// `G(a,b)` for any indices: the one-loop bubble with line powers `a` and
/// `b` at `Q² = 1`, zero (scaleless) for a non-positive integer index.
fn bubble(a: &Atom, b: &Atom) -> Atom {
    if [a, b]
        .iter()
        .any(|index| get_integer_from_atom(index.as_view()).is_some_and(|n| n <= 0))
    {
        return Atom::Zero;
    }
    let e = Atom::var(epsilon());
    let half_d = Atom::num(2) - &e;
    let d = Atom::num(4) - Atom::num(2) * &e;
    gamma().call(a + b - &half_d) * gamma().call(&half_d - a) * gamma().call(&half_d - b)
        / (gamma().call(a.clone()) * gamma().call(b.clone()) * gamma().call(d - a - b))
}

/// `C_G` of a scalar sector whose graph, pinched lines contracted, is
/// series-parallel between the external vertices: two lines between the same
/// vertices (a bubble) are one line of power `a+b-d/2` times `G(a,b)`, two
/// lines through a vertex of degree two are one line of power `a+b`, and a
/// line of power zero is contracted.
fn series_parallel(family: &PFamilyData, powers: &[i64]) -> Atom {
    let lines = canonical_lines(family).unwrap();
    let pinched = powers.iter().map(|&power| power == 0).collect::<Vec<_>>();
    let graph = ContractedGraph::new(&lines, &pinched);
    let (mut q_in, mut q_out) = (graph.find(1), graph.find(2));
    let mut edges = lines
        .iter()
        .zip(powers)
        .filter(|(_, power)| **power != 0)
        .map(|(line, power)| {
            (
                graph.find(line.from),
                graph.find(line.to),
                Atom::num(*power),
            )
        })
        .collect::<Vec<_>>();
    let half_d = Atom::num(2) - Atom::var(epsilon());
    let mut value = Atom::num(1) / g_scheme_measure().pow(Atom::num(family.loops as i64));
    let context = format!("{} {powers:?}", family.name);
    loop {
        if q_in == q_out {
            return Atom::Zero;
        }
        if let Some(index) = edges
            .iter()
            .position(|(from, to, power)| from == to || power.is_zero())
        {
            // A self-loop or coinciding external vertices are scaleless.
            let (from, to, power) = edges.remove(index);
            if from == to && !power.is_zero() {
                return Atom::Zero;
            }
            let (kept, merged) = (from.min(to), from.max(to));
            for end in edges
                .iter_mut()
                .flat_map(|edge| [&mut edge.0, &mut edge.1])
                .chain([&mut q_in, &mut q_out])
            {
                if *end == merged {
                    *end = kept;
                }
            }
            continue;
        }
        let [(from, to, _)] = edges.as_slice() else {
            let same = |i: usize, j: usize| {
                let ((a, b, _), (c, d, _)) = (&edges[i], &edges[j]);
                (a, b) == (c, d) || (a, b) == (d, c)
            };
            if let Some((i, j)) = (0..edges.len())
                .flat_map(|i| (i + 1..edges.len()).map(move |j| (i, j)))
                .find(|&(i, j)| same(i, j))
            {
                let (_, _, b) = edges.remove(j);
                let a = edges[i].2.clone();
                value *= bubble(&a, &b);
                edges[i].2 = (a + b - &half_d).expand();
                continue;
            }
            let degree = |vertex: i64| {
                edges
                    .iter()
                    .filter(|(from, to, _)| *from == vertex || *to == vertex)
                    .count()
            };
            let vertex = edges
                .iter()
                .flat_map(|(from, to, _)| [*from, *to])
                .find(|&vertex| vertex != q_in && vertex != q_out && degree(vertex) == 2)
                .unwrap_or_else(|| panic!("{context} is not series-parallel"));
            let (first, second) = {
                let mut at = (0..edges.len())
                    .filter(|&i| edges[i].0 == vertex || edges[i].1 == vertex)
                    .collect::<Vec<_>>();
                let second = edges.remove(at.pop().unwrap());
                (at[0], second)
            };
            let far = |(from, to, _): &(i64, i64, Atom)| if *from == vertex { *to } else { *from };
            let end = far(&second);
            let edge = &mut edges[first];
            *edge = (far(edge), end, (&edge.2 + second.2).expand());
            continue;
        };
        assert_eq!(
            [*from.min(to), *from.max(to)],
            [q_in.min(q_out), q_in.max(q_out)],
            "{context}"
        );
        return value;
    }
}

fn g_scheme_measure() -> Atom {
    Atom::var(epsilon()) * g(1, Atom::num(1))
}

/// A scalar battery case on a sector with a vertex of five or more momenta
/// once its pinched lines are contracted (the graphs Forcer cannot match):
/// `reduced = (-1)^(Σa+r)·C_G` at `w = 2L-Σa+r`, with `C_G` the
/// series-parallel `G`-function product, `Σa` the positive powers and `r`
/// the degree of the negative powers (lines without `p(1)`).
fn assert_series_parallel(case: &Case, values: &[ReducedPropagatorValue]) {
    let family = family(&case.family);
    assert_eq!(case.numerator, "1", "{}", case.name);
    assert!(
        max_contracted_valence(family, &case.powers) >= 5,
        "{}",
        case.name
    );
    let positive: i64 = case.powers.iter().filter(|a| **a > 0).sum();
    let degree: i64 = -case.powers.iter().filter(|a| **a < 0).sum::<i64>();
    let sign = if (positive + degree) % 2 == 0 { 1 } else { -1 };
    let w = 2 * family.loops as i64 - positive + degree;
    assert_series_zero(
        exact(values, w) - Atom::num(sign) * series_parallel(family, &case.powers),
        2,
        &case.name,
    );
}

/// The largest number of momenta (Q included) at a vertex of the sector
/// graph once the lines of non-positive power are contracted.
fn max_contracted_valence(family: &PFamilyData, powers: &[i64]) -> usize {
    let lines = canonical_lines(family).unwrap();
    let pinched = powers.iter().map(|&power| power <= 0).collect::<Vec<_>>();
    let graph = ContractedGraph::new(&lines, &pinched);
    let mut momenta = BTreeMap::<i64, usize>::new();
    let ends = lines
        .iter()
        .zip(&pinched)
        .filter(|(_, pinched)| !**pinched)
        .flat_map(|(line, _)| [line.from, line.to])
        .chain([1, 2]);
    for vertex in ends {
        *momenta.entry(graph.find(vertex)).or_default() += 1;
    }
    momenta.into_values().max().unwrap()
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
        .map_or(if case.numerator == "1" { 0 } else { 1 }, |n| {
            n.parse().unwrap()
        });
    let sigma = if case.external.starts_with('-') {
        -1
    } else {
        1
    };
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
    let coefficients = finalize_tests::coefficients_at_unit_scale(&kite, &settings, &head, -2, 0);
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
        if case.name.starts_with("contracted") {
            assert_series_parallel(case, values);
        }
        // The scalar planar and non-planar tops have the G-scheme value 20ζ5.
        if ["la_top", "no_top"].contains(&case.name.as_str()) {
            let result = finalize(
                values,
                3,
                &head,
                &settings,
                RustRedMasterMode::Symbolic,
                "Forcer",
            )
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
    // Master(no) is tabulated through ε^7 and Master(t105) through ε^9; LA
    // multiplies t105 by PG(0,0)(3-6ε)/ε, so two further poles reach ε^9 at
    // nine terms (ε^5).
    let value = |name: &str| {
        &battery
            .iter()
            .find(|(case, _)| case.name == name)
            .unwrap()
            .1
    };
    assert_guard(value("no_top"), 3, 1, 10, ("no", 8));
    assert_guard(value("la_top"), 3, 2, 9, ("t105", 10));
}

/// The `Oep` guard on Forcer results: with an extra `1/ε^pole`, `deepest_terms`
/// Laurent terms succeed and one more needs the unknown order
/// `first_missing_order` of `master`.
fn assert_guard(
    values: &[ReducedPropagatorValue],
    loops: usize,
    pole: i64,
    deepest_terms: i64,
    (master, first_missing_order): (&str, i64),
) {
    let head = parse("p(1)");
    let evaluate = |terms: i64| {
        let settings = finalize_tests::settings(LoopNormalizationFactor::GScheme, terms);
        let e = Atom::var(vk_symbol!(&settings.epsilon_symbol));
        let values = values
            .iter()
            .map(|value| ReducedPropagatorValue {
                scale_power: value.scale_power,
                reduced: &value.reduced * e.pow(Atom::num(-pole)),
            })
            .collect::<Vec<_>>();
        finalize(
            &values,
            loops,
            &head,
            &settings,
            RustRedMasterMode::Symbolic,
            "Forcer",
        )
    };
    let within =
        evaluate(deepest_terms).unwrap_or_else(|error| panic!("within the master depth: {error}"));
    assert!(!within.is_zero());
    match evaluate(deepest_terms + 1) {
        Err(VakintError::MissingMasterOrders(missing)) => {
            assert_eq!(missing.backend, "Forcer");
            assert!(
                missing
                    .masters
                    .iter()
                    .any(|m| m.master.contains(master)
                        && m.first_missing_order == first_missing_order),
                "{missing}"
            );
        }
        other => panic!("expected missing master orders, got {other:?}"),
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
        if case.name.starts_with("contracted") {
            assert_series_parallel(case, values);
        }
    }
    // dot(k(3),p(1)) on lines k(3) (3) and k(3)+p(1) (6) of a sector with a
    // six-momentum vertex: k·Q = ((k+Q)² - k² - Q²)/2 in G-functions.
    let (case, values) = battery
        .iter()
        .find(|(case, _)| case.name == "la4_kq_six_momenta")
        .unwrap();
    let la4 = family("P4L_LA4");
    assert!(max_contracted_valence(la4, &case.powers) >= 6);
    let lowered = |line: usize| {
        let mut powers = case.powers.clone();
        powers[line] -= 1;
        series_parallel(la4, &powers)
    };
    let k_dot_q = (lowered(5) - lowered(2) - series_parallel(la4, &case.powers)) / Atom::num(2);
    let positive: i64 = case.powers.iter().sum();
    assert_series_zero(
        exact(values, 8 - positive + 1)
            - Atom::num(if positive % 2 == 0 { -1 } else { 1 }) * k_dot_q,
        2,
        &case.name,
    );
    // The ladder p-integral is C(2L,L)·ζ(2L-1) in the G-scheme: 70ζ7 at four loops.
    let head = parse("p(1)");
    let settings = finalize_tests::settings(LoopNormalizationFactor::GScheme, 5);
    let (_, ladder) = battery
        .iter()
        .find(|(case, _)| case.name == "la4_top")
        .unwrap();
    let result = finalize(
        ladder,
        4,
        &head,
        &settings,
        RustRedMasterMode::Symbolic,
        "Forcer",
    )
    .unwrap();
    let coefficients = finalize_tests::coefficients_at_unit_scale(&result, &settings, &head, -4, 0);
    assert_eq!(
        coefficients,
        [
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            (Atom::num(-70) * symbolica::transcendental::zeta().call(7)).expand()
        ]
    );
}

#[test]
fn tform_workers_reproduce_the_frozen_output() {
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the TFORM run");
        return;
    };
    let battery = Battery::parse(include_str!(
        "../../../tests/fixtures/forcer/p3l_battery.json"
    ))
    .unwrap();
    let case = battery
        .cases
        .iter()
        .find(|case| case.name == "la_rank4")
        .unwrap();
    let settings = VakintSettings {
        n_cores: 2,
        ..settings(&form)
    };
    let reduced = reduce_p(
        family(&case.family),
        &case.powers,
        parse(&case.numerator).as_view(),
        &parse(&case.external),
        &settings,
    )
    .unwrap();
    assert_eq!(reduced, frozen(case));
}

#[test]
fn tagged_homogeneous_channels_match_independent_reductions() {
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the tagged-channel run");
        return;
    };
    let family = family("P2L");
    let powers = [1, 1, 1, 1, 1];
    let external = parse("-p(3)");
    let numerators = [
        parse("1"),
        parse("(alpha+ε)*dot(k(1),p(3))"),
        parse("2*dot(k(1),k(2))^2"),
        parse("beta*dot(k(2),k(2))"),
    ];
    let numerator = Atom::add_many(&numerators);
    let input = forcer_input(family, &powers, numerator.as_view(), &external)
        .unwrap()
        .unwrap();
    assert!(input.parts.len() >= 3);
    assert!(input.program.contains("Local VakintChannels ="));
    assert_eq!(input.program.matches("#call Forcer()").count(), 1);

    let mut expected = BTreeMap::<i64, Atom>::new();
    for numerator in numerators {
        for value in reduce_p(
            family,
            &powers,
            numerator.as_view(),
            &external,
            &settings(&form),
        )
        .unwrap()
        {
            *expected.entry(value.scale_power).or_default() += value.reduced;
        }
    }
    for cores in [1, 2] {
        let settings = VakintSettings {
            n_cores: cores,
            ..settings(&form)
        };
        let actual = reduce_p(family, &powers, numerator.as_view(), &external, &settings).unwrap();
        assert_eq!(actual.len(), expected.len());
        for value in actual {
            assert!(
                (value.reduced - &expected[&value.scale_power])
                    .expand()
                    .together()
                    .is_zero(),
                "tagged channel changed scale {} with {cores} workers",
                value.scale_power
            );
        }
    }
}

#[test]
fn tagged_four_loop_channels_reproduce_frozen_outputs() {
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the four-loop tagged-channel run");
        return;
    };
    let battery = Battery::parse(include_str!(
        "../../../tests/fixtures/forcer/p4l_battery.json"
    ))
    .unwrap();
    let mut numerator = Atom::Zero;
    let mut expected = BTreeMap::<i64, Atom>::new();
    let mut selected = None;
    for (name, coefficient) in [("la4_top", parse("1")), ("la4_isp", parse("1+ε"))] {
        let case = battery.cases.iter().find(|case| case.name == name).unwrap();
        numerator += &coefficient * parse(&case.numerator);
        for value in frozen(case) {
            *expected.entry(value.scale_power).or_default() += &coefficient * value.reduced;
        }
        selected = Some(case);
    }
    let case = selected.unwrap();
    let actual = reduce_p(
        family(&case.family),
        &case.powers,
        numerator.as_view(),
        &parse(&case.external),
        &settings(&form),
    )
    .unwrap();
    assert_eq!(actual.len(), expected.len());
    for value in actual {
        assert!(
            (value.reduced - &expected[&value.scale_power])
                .expand()
                .together()
                .is_zero(),
            "tagged four-loop channel changed scale {}",
            value.scale_power
        );
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
        assert!(
            input.program.contains(fragment),
            "{fragment}\n{}",
            input.program
        );
    }
    assert!(
        input.program.contains("id k1 = +p1;\nid k2 = +p2;"),
        "{}",
        input.program
    );
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

/// The `vx(...)` factors of each local of a Forcer program.
fn forcer_vertices(input: &ForcerInput) -> Vec<Vec<String>> {
    let vertex = Regex::new(r"vx\([^)]*\)").unwrap();
    input
        .program
        .split_once("#call Forcer()")
        .expect("a Forcer program must call Forcer")
        .0
        .lines()
        .filter(|line| line.starts_with("Local "))
        .map(|local| {
            vertex
                .find_iter(local)
                .map(|found| found.as_str().to_owned())
                .collect()
        })
        .collect()
}

#[test]
fn pinched_sectors_keep_the_cubic_family_graph() {
    // Forcer silently drops a graph that is none of its topologies (at most
    // four momenta per vertex): contracting pinched lines would merge
    // vertices into five or more momenta, so every sector is the family graph
    // with the pinched lines left without propagator.
    for family in &P_FAMILIES {
        let n = family.line_count();
        let input = |powers: &[i64]| {
            forcer_input(family, powers, Atom::num(1).as_view(), &parse("p(1)")).unwrap()
        };
        let [top] = &forcer_vertices(&input(&vec![1; n]).unwrap())[..] else {
            panic!("{}: one local", family.name);
        };
        assert!(top.iter().all(|vertex| vertex.matches(',').count() == 2));
        for first in 0..n {
            for second in first + 1..n {
                let mut powers = vec![1; n];
                (powers[first], powers[second]) = (0, -1);
                if let Some(sector) = input(&powers) {
                    for local in forcer_vertices(&sector) {
                        assert_eq!(&local, top, "{} {powers:?}", family.name);
                    }
                    for pinched in [first, second] {
                        let name = format!("p{}", pinched + 1);
                        assert!(!sector.program.contains(&format!("/({name}.{name})")));
                    }
                }
            }
        }
    }
}

#[test]
fn forcer_warnings_for_unrecognized_graphs_are_errors() {
    let stdout = "Time =       0.29 sec    Generated terms =          1\n\
        WARNING: the following diagram did not match a topology:  + acc(v(p1,p2,p6,p9,\n \
        p10)*v(p1,p6,Q)*EX)\n\nTime =       0.30 sec";
    assert_eq!(
        unrecognized_topology(stdout).as_deref(),
        Some(
            "WARNING: the following diagram did not match a topology: + acc(v(p1,p2,p6,p9, \
             p10)*v(p1,p6,Q)*EX)"
        )
    );
    assert_eq!(unrecognized_topology("Time = 0.3 sec\n  F1 = 0;"), None);
}

#[test]
fn forcer_refuses_its_silent_zero() {
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the unrecognized-graph run");
        return;
    };
    // A four-loop graph with a five-momentum vertex: no Forcer topology, so
    // Forcer alone returns 0 (rules-study verify-crossed, valence_test.frm).
    let program = "#-\n#include- forcer.h\nFormat 255;\nLocal F1 = 1/p1.p1/p2.p2/p3.p3/p4.p4\
        /p6.p6/p8.p8/p5.p5/p9.p9/p10.p10*vx(Q,-p1,p6)*vx(p1,-p2,-p6,p9,p10)*vx(p2,-p3,p8)\
        *vx(p3,p4,-Q)*vx(-p4,-p8,-p5)*vx(p5,-p9,-p10);\n#call Forcer()\n.sort\n\
        #write <F1.out> \"%E\", F1\n.end\n";
    match run_forcer(&settings(&form), program, &["F1"]) {
        Err(VakintError::ForcerUnrecognizedTopology(warning)) => {
            assert!(
                warning.replace(' ', "").contains("v(p1,p2,p6,p9,p10)"),
                "{warning}"
            )
        }
        other => panic!("expected ForcerUnrecognizedTopology, got {other:?}"),
    }
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
                    assert_eq!(
                        total,
                        expected,
                        "{}: k({}) column {column}",
                        family.name,
                        j + 1
                    );
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
        let values = reduce_p(
            family(name),
            &powers,
            Atom::num(1).as_view(),
            &parse("p(1)"),
            &settings,
        )
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
        assert_eq!(
            graph.is_scaleless(&lines.iter().collect::<Vec<_>>()),
            scaleless
        );
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
    let expected = parse("PG(0,0)^2*(2*ε-1)/(2*ε)+PG(0,0)*PG(1,0)*(1-3*ε)/(2*ε)+PMaster(no)");
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
        message.contains("needs FORM >= 5.0.0") && !message.contains("RustRed"),
        "{message}"
    );
}

/// The current user's shared source directory under `root`, which the
/// test owns.
#[cfg(unix)]
fn shared_directory(root: &std::path::Path) -> (u32, PathBuf) {
    use std::os::unix::fs::MetadataExt;
    let owner = std::fs::metadata(root).unwrap().uid();
    (owner, root.join(shared_directory_name(owner)))
}

#[cfg(unix)]
#[test]
fn vendored_sources_materialize_once_and_completely() {
    let directory = tempfile::tempdir().unwrap();
    let base = directory.path().to_string_lossy().into_owned();
    let first = materialized_sources(Some(&base)).unwrap().unwrap();
    assert_eq!(
        Some(&first),
        materialized_sources(Some(&base)).unwrap().as_ref()
    );
    let (owner, shared) = shared_directory(directory.path());
    assert_eq!(first, shared);
    assert!(is_trusted(&first, owner, &vendored_sources().unwrap()));
    assert_eq!(
        std::fs::read(first.join("forcer").join("Masters.prc")).unwrap(),
        include_bytes!("../../../form_src/forcer/forcer/Masters.prc")
    );
    assert_eq!(vendored::FORCER_FILES.len(), 21);
}

/// Whether FORM gets the shared source directory once it was planted (a
/// complete private copy) and then altered by `tamper`, and whether the
/// sources FORM gets are trusted. No private source directory is left beside
/// the run directory.
#[cfg(unix)]
fn sources_used_after(tamper: impl Fn(&std::path::Path)) -> (bool, bool) {
    let root = tempfile::tempdir().unwrap();
    let (owner, shared) = shared_directory(root.path());
    std::fs::create_dir(&shared).unwrap();
    make_private(&shared).unwrap();
    write_sources(&shared, &vendored_sources().unwrap()).unwrap();
    tamper(&shared);
    let base = root.path().to_string_lossy().into_owned();
    let used = materialized_sources(Some(&base)).unwrap();
    let (run, private_sources) = prepare_run_directory(&VakintSettings {
        temporary_directory: Some(base),
        ..VakintSettings::default()
    })
    .unwrap();
    assert_eq!(private_sources, used.is_none());
    let sources = vendored_sources().unwrap();
    let trusted = match &used {
        Some(directory) => is_trusted(directory, owner, &sources),
        // The private copy inside the run directory.
        None => is_trusted(&run, owner, &sources),
    };
    std::fs::remove_dir_all(&run).unwrap();
    assert_eq!(leftovers(root.path()), Vec::<String>::new());
    (used.as_ref() == Some(&shared), trusted)
}

#[cfg(unix)]
#[test]
fn untrusted_source_directories_are_never_used() {
    use std::os::unix::fs::PermissionsExt;
    // A complete private copy of the current user is shared.
    assert_eq!(sources_used_after(|_| ()), (true, true));
    // A forcer.h that runs a shell command when FORM includes it.
    assert_eq!(
        sources_used_after(|shared| {
            std::fs::write(shared.join("forcer.h"), "#system touch PLANTED\n").unwrap()
        }),
        (false, true)
    );
    // A partly deleted copy.
    assert_eq!(
        sources_used_after(|shared| {
            std::fs::remove_file(shared.join("forcer").join("Masters.prc")).unwrap()
        }),
        (false, true)
    );
    // Sources that other users can replace.
    assert_eq!(
        sources_used_after(|shared| {
            std::fs::set_permissions(shared, std::fs::Permissions::from_mode(0o777)).unwrap()
        }),
        (false, true)
    );
    // A link to a complete copy elsewhere.
    assert_eq!(
        sources_used_after(|shared| {
            let elsewhere = shared.with_extension("elsewhere");
            std::fs::rename(shared, &elsewhere).unwrap();
            std::os::unix::fs::symlink(&elsewhere, shared).unwrap();
        }),
        (false, true)
    );
}

/// A shared source directory of the current user that is not trusted (one
/// source file altered), under `root`.
#[cfg(unix)]
fn plant_untrusted_shared_directory(root: &std::path::Path) -> PathBuf {
    let (_, shared) = shared_directory(root);
    std::fs::create_dir(&shared).unwrap();
    make_private(&shared).unwrap();
    write_sources(&shared, &vendored_sources().unwrap()).unwrap();
    std::fs::write(shared.join("forcer.h"), "* altered\n").unwrap();
    shared
}

/// The entries of `root` other than the shared source directory (and its
/// moved-away copies).
#[cfg(unix)]
fn leftovers(root: &std::path::Path) -> Vec<String> {
    let mut names = std::fs::read_dir(root)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| !name.starts_with("vakint-forcer-"))
        .collect::<Vec<_>>();
    names.sort();
    names
}

#[cfg(unix)]
#[test]
fn untrusted_sources_leave_no_private_directory_behind() {
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live run beside untrusted sources");
        return;
    };
    let root = tempfile::tempdir().unwrap();
    let shared = plant_untrusted_shared_directory(root.path());
    let settings = VakintSettings {
        temporary_directory: Some(root.path().to_string_lossy().into_owned()),
        clean_tmp_dir: true,
        ..settings(&form)
    };
    let reduced = reduce_p(
        family("P1L"),
        &[1, 1],
        parse("1").as_view(),
        &parse("p(1)"),
        &settings,
    )
    .unwrap();
    assert_eq!(reduced.len(), 1);
    assert_eq!(leftovers(root.path()), Vec::<String>::new());
    assert_eq!(
        std::fs::read(shared.join("forcer.h")).unwrap(),
        b"* altered\n",
        "the planted directory is never touched"
    );
}

/// A run that fails beside untrusted sources keeps its program for diagnosis
/// but not the private copy of the sources.
#[cfg(unix)]
fn assert_failed_run_removes_private_sources(form: &str, program: &str) {
    let root = tempfile::tempdir().unwrap();
    plant_untrusted_shared_directory(root.path());
    let settings = VakintSettings {
        temporary_directory: Some(root.path().to_string_lossy().into_owned()),
        clean_tmp_dir: true,
        ..settings(form)
    };
    run_forcer(&settings, program, &["F1"]).unwrap_err();
    let run = leftovers(root.path());
    assert_eq!(run.len(), 1, "{run:?}");
    let run = root.path().join(&run[0]);
    assert!(run.join("run_forcer.frm").is_file());
    assert!(std::fs::symlink_metadata(run.join("forcer.h")).is_err());
    assert!(std::fs::symlink_metadata(run.join("forcer")).is_err());
}

#[cfg(unix)]
#[test]
fn a_failed_run_leaves_no_private_sources_behind() {
    // FORM cannot be started.
    assert_failed_run_removes_private_sources("/nonexistent/form", ".end\n");
    let Some(form) = forcer_test_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live failures beside untrusted sources");
        return;
    };
    // FORM fails, and FORM succeeds without the expected output.
    assert_failed_run_removes_private_sources(&form, "Local F1 = ;\n.end\n");
    assert_failed_run_removes_private_sources(&form, ".end\n");
}

#[test]
#[ignore = "profiling"]
fn profile_four_loop_finalize() {
    let battery = Battery::parse(include_str!(
        "../../../tests/fixtures/forcer/p4l_battery.json"
    ))
    .unwrap();
    for case in &battery.cases {
        let values = frozen(case);
        for terms in [4, 5, 6] {
            let settings = finalize_tests::settings(LoopNormalizationFactor::MSbar, terms);
            let start = std::time::Instant::now();
            let result = finalize(
                &values,
                4,
                &parse(&case.external),
                &settings,
                RustRedMasterMode::Symbolic,
                "Forcer",
            );
            eprintln!(
                "{} terms={terms}: {:?}, {}",
                case.name,
                start.elapsed(),
                result.map_or_else(|e| e.to_string(), |r| r.to_string().len().to_string())
            );
        }
    }
}

const BENCHMARK_INPUTS: &str = include_str!("../../../tests/fixtures/benchmark_p_integrals.json");

/// `blake3` of the Forcer programs of every fixture case (name, program).
/// The two- and three-loop battery digests include the tagged-channel
/// rendering; all other cases still render as before the shared numerator
/// split moved to `parts.rs`. The grid and anchor digests are taken over
/// their definitions (`grid_cases`, `anchor_cases`), which
/// `tests/p_integral_forcer_fixtures.rs` ties to the files.
const PROGRAM_DIGESTS: [(&str, &str); 8] = [
    (
        "p1l_battery.json",
        "11175f6ceaf2370815cf4e1095dfd6ffaee092d957ff0de954b0e4eeebba5a6b",
    ),
    (
        "p2l_battery.json",
        "45ae11c400c0c50914c74116cc8791b2061962e3a518dfeec075b4b508db53da",
    ),
    (
        "p3l_battery.json",
        "a0ffe7cdef48253703f950afa73e4875cdcf324fefbdc26e6eccdae192a99ed1",
    ),
    (
        "p4l_battery.json",
        "4db77a0e957fc721fc0b981feab8a2197b43828e4973aecf19520343356e6596",
    ),
    (
        "p_anchors.json",
        "aab2618020f36cce87741607310144342b15e1ff77bdc09bc433153c5d6cc0ba",
    ),
    (
        "p1l_grid.json",
        "5869e41489de001c4f4c897c546e84e2ab8b2b06093b2d2d7dc4ebb2227fd849",
    ),
    (
        "p2l_grid.json",
        "44945f050778bbe6fadce599daeca316f3b73faffb835386a6065b8e1b534f2a",
    ),
    (
        "p3l_grid.json",
        "a084c24715b21fc26c5b92f3cdea441510385b5c5b4ad3f872620e63c01924ae",
    ),
];

/// The cases of a fixture file, without outputs.
fn fixture_definitions(file: &str) -> Battery {
    let grid = |loops| {
        Battery::unwritten(
            super::fixtures::GRID_PROVENANCE,
            super::fixtures::grid_cases(loops, BENCHMARK_INPUTS).unwrap(),
        )
    };
    match file {
        "p_anchors.json" => Battery::unwritten(
            super::fixtures::ANCHOR_PROVENANCE,
            super::fixtures::anchor_cases(),
        ),
        "p1l_grid.json" => grid(1),
        "p2l_grid.json" => grid(2),
        "p3l_grid.json" => grid(3),
        battery => {
            Battery::parse(&std::fs::read_to_string(fixture_path(battery)).unwrap()).unwrap()
        }
    }
}

#[test]
fn forcer_programs_are_byte_identical() {
    let found = PROGRAM_DIGESTS
        .iter()
        .map(|(file, _)| (*file, fixture_definitions(file).program_digest().unwrap()))
        .collect::<Vec<_>>();
    for ((file, expected), (_, digest)) in PROGRAM_DIGESTS.iter().zip(&found) {
        assert_eq!(
            expected, digest,
            "{file}: Forcer programs changed; all digests: {found:#?}"
        );
    }
}

#[test]
fn the_grid_covers_every_nonzero_labelled_sector() {
    let count = |loops: usize| {
        P_FAMILIES
            .iter()
            .filter(|family| family.loops == loops)
            .map(|family| super::fixtures::nonzero_sectors(family.name).unwrap().len())
            .sum::<usize>()
    };
    assert_eq!([count(1), count(2), count(3)], [1, 8, 175]);
    let names = |loops| {
        super::fixtures::grid_cases(loops, BENCHMARK_INPUTS)
            .unwrap()
            .into_iter()
            .map(|case| case.name)
            .collect::<Vec<_>>()
    };
    let (one, two, three) = (names(1), names(2), names(3));
    let kind =
        |names: &[String], prefix: &str| names.iter().filter(|n| n.starts_with(prefix)).count();
    assert_eq!(
        [
            one.len(),
            kind(&one, "bench_"),
            kind(&one, "sector_"),
            kind(&one, "dotted_")
        ],
        [17, 15, 2, 0]
    );
    assert_eq!(
        [
            two.len(),
            kind(&two, "bench_"),
            kind(&two, "sector_"),
            kind(&two, "dotted_")
        ],
        [39, 15, 16, 8]
    );
    assert_eq!(
        [
            three.len(),
            kind(&three, "bench_"),
            kind(&three, "sector_"),
            kind(&three, "dotted_")
        ],
        [417, 45, 350, 22]
    );
    for names in [&one, &two, &three] {
        let mut unique = names.clone();
        unique.sort();
        unique.dedup();
        assert_eq!(unique.len(), names.len());
    }
}

/// `p_anchors.json`, checked against its definition.
fn anchors() -> Battery {
    super::fixtures::write_or_read(
        std::path::Path::new(&fixture_path("p_anchors.json")),
        super::fixtures::ANCHOR_PROVENANCE,
        super::fixtures::anchor_cases(),
        None,
    )
    .unwrap()
}

/// The frozen values of an anchor in the tag basis, per scale power.
fn anchor_basis(anchors: &Battery, name: &str) -> Vec<(i64, super::super::compare::TagBasis)> {
    let case = anchors.cases.iter().find(|case| case.name == name).unwrap();
    frozen(case)
        .into_iter()
        .map(|value| {
            let basis = super::super::compare::tag_basis(&value.reduced, epsilon()).unwrap();
            (value.scale_power, basis)
        })
        .collect()
}

#[test]
fn the_anchors_pin_the_normalization() {
    use super::super::compare::TagBasis;
    let anchors = anchors();
    let e = Atom::var(epsilon());
    let unit = |coefficient: Atom| TagBasis::from([((vec![], vec![]), coefficient.together())]);
    // P1L(1,1) = 1/ε: w = 0, Wick sign (-1)^2.
    assert_eq!(
        anchor_basis(&anchors, "p1l_11"),
        [(0, unit(Atom::num(1) / &e))]
    );
    // k·p on the bubble: w = 1, sign (-1)^3, odd in the sign of p.
    assert_eq!(
        anchor_basis(&anchors, "p1l_kq_plus"),
        [(1, unit(Atom::num(-1) / (Atom::num(2) * &e)))]
    );
    assert_eq!(
        anchor_basis(&anchors, "p1l_kq_minus"),
        [(1, unit(Atom::num(1) / (Atom::num(2) * &e)))]
    );
    // NO(1^8, 0) is the master itself: I_M = +(εG(1,1))^3 Master(no) at w = -2.
    assert_eq!(
        anchor_basis(&anchors, "no"),
        [(
            -2,
            TagBasis::from([((vec![], vec![("no".to_string(), 1)]), Atom::num(1))])
        )]
    );
    // The kite: finalized ε^0 = (-1)^5 · 6ζ3.
    let head = parse("p(1)");
    let settings = finalize_tests::settings(LoopNormalizationFactor::GScheme, 3);
    let kite = anchors
        .cases
        .iter()
        .find(|case| case.name == "kite")
        .unwrap();
    let value = kite
        .finalized(&settings, RustRedMasterMode::Symbolic)
        .unwrap();
    assert_eq!(
        finalize_tests::coefficients_at_unit_scale(&value, &settings, &head, -2, 0),
        [
            Atom::Zero,
            Atom::Zero,
            (Atom::num(-6) * symbolica::transcendental::zeta().call(3)).expand()
        ]
    );
    // The rank-3 kites carry (Q²)^w, w = 2L - Σa + (number of scalar
    // products): (k1·p)²(k2·p) has three, (k1·p)(k1·k2) two.
    for (name, w) in [("kite_rank3_kp", 2), ("kite_rank3_mixed", 1)] {
        let basis = anchor_basis(&anchors, name);
        assert_eq!(basis.len(), 1, "{name}");
        assert_eq!(basis[0].0, w, "{name}");
    }
    // The plain three-loop kite with a doubled central line integrates to
    // t1star05 with its ε-line at n = 0: the master t105 appears, and the LA
    // and NO routings of the same sector agree exactly.
    let la = anchor_basis(&anchors, "t105_parent_la");
    assert_eq!(la, anchor_basis(&anchors, "t105_parent_no"));
    assert_eq!(la.len(), 1);
    assert!(
        la[0]
            .1
            .keys()
            .any(|(_, masters)| masters == &[("t105".to_string(), 1)]),
        "{la:?}"
    );
    // A dot on a bubble line (n = 1) leaves the master alone.
    let dotted = anchor_basis(&anchors, "t105_parent_la_dotted");
    assert!(
        dotted[0]
            .1
            .keys()
            .all(|(_, masters)| masters == &[("t105".to_string(), 1)]),
        "{dotted:?}"
    );
    for case in &anchors.cases {
        let value = case
            .finalized(&settings, RustRedMasterMode::Retain)
            .unwrap_or_else(|error| panic!("{}: {error}", case.name));
        assert!(!value.is_zero(), "{}", case.name);
    }
}

#[test]
fn the_benchmark_inputs_are_sums_of_registry_tops() {
    let inputs = super::fixtures::BenchmarkInputs::parse(BENCHMARK_INPUTS).unwrap();
    let mut cells = Vec::new();
    for cell in &inputs.cells {
        cells.push((cell.loops, cell.dots, cell.rank));
        assert_eq!(
            cell.name,
            format!("L{}_d{}_r{}", cell.loops, cell.dots, cell.rank)
        );
        let tops = P_FAMILIES
            .iter()
            .filter(|family| family.loops == cell.loops)
            .collect::<Vec<_>>();
        assert_eq!(
            cell.tops
                .iter()
                .map(|top| top.family.as_str())
                .collect::<Vec<_>>(),
            tops.iter().map(|family| family.name).collect::<Vec<_>>()
        );
        let lines = 3 * cell.loops - 1;
        let mut powers = vec![1; lines];
        if cell.dots >= 1 {
            powers[0] = 2;
        }
        if cell.dots >= 2 {
            powers[lines - 1] = 2;
        }
        let mut terms = Vec::new();
        for (top, family) in cell.tops.iter().zip(&tops) {
            assert_eq!(top.powers, powers, "{}", cell.name);
            let mut long = family.lines.to_string();
            for (line, power) in powers.iter().enumerate().rev() {
                long = long.replace(&format!("pow({})", line + 1), &power.to_string());
            }
            terms.push(format!("({})*{long}", cell.numerator));
        }
        assert_eq!(cell.input, terms.join("+"), "{}", cell.name);
        // The numerator has rank r: r loop momenta in its scalar products.
        assert_eq!(
            cell.numerator.matches("k(").count(),
            cell.rank,
            "{}",
            cell.name
        );
        assert!(!parse(&cell.input).is_zero(), "{}", cell.name);
    }
    let mut grid = Vec::new();
    for loops in 1..=4 {
        for dots in 0..=2 {
            for rank in [0, 2, 4, 6, 8] {
                grid.push((loops, dots, rank));
            }
        }
    }
    assert_eq!(cells, grid);
}
