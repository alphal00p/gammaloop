//! Equivalent routings of one vacuum integrand evaluate identically.
//!
//! A vacuum integral is fixed by its graph, its line masses and powers and its
//! numerator; the momentum routing written in `topo(prop(...))` only chooses
//! the integration variables. Each test below writes a graph once, numerators
//! in its line momenta `line(e, μ)` and a list of equivalent routings: lines
//! written reversed with a negated momentum, loop momenta permuted or negated,
//! and loop momenta that are combinations of line momenta (a loop momentum need
//! not be the momentum of any propagator). [`assert_routing_invariant`]
//! evaluates every routing with each listed backend, through the one-shot
//! `evaluate` and the staged `evaluate_integral` entry points, and requires
//! exactly the same result. The first routing of every list is vakint's own
//! canonical routing, so the comparison is also an absolute check;
//! `two_loop_sunset_numerator_matches_its_scalar_decomposition` anchors that
//! routing independently of any numerator substitution.
//!
//! The matcher-level contract for every registered topology at one to five
//! loops is the unit test module `topologies::routing_tests::user_routings`
//! (five-loop evaluations are far too slow for this suite).

mod test_utils;

use symbolica::atom::{Atom, AtomCore};
use test_utils::rustred_support::rustred_only;
use vakint::{EvaluationOrder, RustRedMasterMode, VakintSettings, vakint_parse};

/// RustRed never runs FORM; a FORM call would fail on this path.
const NO_FORM: &str = "/vakint-routing-invariance-tests-never-run-form";

/// A line `from -> to` of a vacuum graph and its propagator power. Every line
/// has the mass `muvsq`.
type Line = (usize, usize, i64);

/// The momenta of the lines, in the graph's orientation, for one routing. The
/// lines listed in `reversed` are written as `edge(to, from)` with the negated
/// momentum: the same propagator, read in the opposite direction.
#[derive(Clone, Copy)]
struct Routing {
    momenta: &'static [&'static str],
    reversed: &'static [usize],
}

const fn routing(momenta: &'static [&'static str]) -> Routing {
    Routing {
        momenta,
        reversed: &[],
    }
}

const fn reversed(momenta: &'static [&'static str], reversed: &'static [usize]) -> Routing {
    Routing { momenta, reversed }
}

#[derive(Clone, Copy, Debug)]
enum Backend {
    /// RustRed with analytic masters, no FORM. (Retained masters are not
    /// canonical: equal integrals may keep symmetric-equivalent master labels.)
    RustRed,
    Matad,
    AlphaLoop,
}

/// Every backend of one to three loops.
const LOWER_LOOP_BACKENDS: &[Backend] = &[Backend::RustRed, Backend::Matad, Backend::AlphaLoop];

impl Backend {
    fn settings(self, loops: usize) -> VakintSettings {
        let settings = match self {
            Backend::RustRed => rustred_only(RustRedMasterMode::Symbolic, NO_FORM),
            Backend::Matad => VakintSettings {
                evaluation_order: EvaluationOrder::matad_only(None),
                ..VakintSettings::default()
            },
            Backend::AlphaLoop => VakintSettings {
                evaluation_order: EvaluationOrder::alphaloop_only(),
                ..VakintSettings::default()
            },
        };
        VakintSettings {
            number_of_terms_in_epsilon_expansion: loops as i64 + 1,
            ..settings
        }
    }
}

/// `numerator * topo(...)` in the given routing; `line(e, μ)` in the numerator
/// is the momentum of line `e` in the graph's orientation.
fn integrand(lines: &[Line], routing: Routing, numerator: &str) -> Atom {
    assert_eq!(lines.len(), routing.momenta.len(), "one momentum per line");
    let mut numerator = vakint_parse!(numerator).unwrap();
    let mut propagators = Vec::new();
    for (index, (&(from, to, power), momentum)) in lines.iter().zip(routing.momenta).enumerate() {
        let id = index + 1;
        propagators.push(if routing.reversed.contains(&id) {
            format!("prop({id},edge({to},{from}),-({momentum}),muvsq,{power})")
        } else {
            format!("prop({id},edge({from},{to}),{momentum},muvsq,{power})")
        });
        let indexed = vakint_parse!(momentum)
            .unwrap()
            .replace(vakint_parse!("k(i_)").unwrap().to_pattern())
            .allow_new_wildcards_on_rhs(true)
            .with(vakint_parse!("k(i_,mu_)").unwrap().to_pattern());
        numerator = numerator
            .replace(
                vakint_parse!(format!("line({id},mu_)"))
                    .unwrap()
                    .to_pattern(),
            )
            .with(indexed.to_pattern());
    }
    numerator * vakint_parse!(format!("topo({})", propagators.join("*"))).unwrap()
}

/// Evaluate `input` with the one-shot `evaluate` and, for numerators without
/// external vectors (the staged entry point does not tensor-reduce), with the
/// staged `evaluate_integral`.
fn evaluate_both_ways(backend: Backend, loops: usize, input: &Atom) -> Vec<(&'static str, Atom)> {
    let vakint = test_utils::get_vakint(backend.settings(loops));
    let fail = |entry: &str, error: vakint::VakintError| -> Atom {
        panic!("{backend:?} {entry} failed on {input}: {error}")
    };
    let mut values = vec![(
        "evaluate",
        vakint
            .evaluate(input.as_view())
            .unwrap_or_else(|error| fail("evaluate", error)),
    )];
    if !input.contains_symbol(vakint::vakint_symbol!("p")) {
        values.push((
            "evaluate_integral",
            vakint
                .evaluate_integral(input.as_view())
                .unwrap_or_else(|error| fail("evaluate_integral", error)),
        ));
    }
    values
}

struct Case {
    loops: usize,
    lines: &'static [Line],
    numerators: &'static [&'static str],
    routings: Vec<Routing>,
    backends: &'static [Backend],
}

/// Every routing of `case` gives exactly the result of the first one, for each
/// numerator, backend and entry point.
fn assert_routing_invariant(case: Case) {
    test_utils::run_multi_lane_acceptance(move || {
        if case.loops >= 4 {
            test_utils::ensure_external_data_loaded();
        }
        for &backend in case.backends {
            for numerator in case.numerators {
                let inputs = case
                    .routings
                    .iter()
                    .map(|routing| integrand(case.lines, *routing, numerator))
                    .collect::<Vec<_>>();
                let reference = evaluate_both_ways(backend, case.loops, &inputs[0]);
                for input in &inputs[1..] {
                    let values = evaluate_both_ways(backend, case.loops, input);
                    for ((entry, value), (_, expected)) in values.iter().zip(&reference) {
                        assert!(
                            (value - expected).expand().is_zero(),
                            "{backend:?} {entry}: routings disagree for numerator {numerator}\n  \
                             {}\n    = {expected}\n  {input}\n    = {value}",
                            inputs[0]
                        );
                    }
                }
            }
        }
    });
}

const TADPOLE: &[Line] = &[(1, 1, 2)];

#[test]
fn one_loop_tadpole() {
    assert_routing_invariant(Case {
        loops: 1,
        lines: TADPOLE,
        numerators: &["line(1,1)*line(1,1)", "line(1,1)*p(1,1)*line(1,2)*p(2,2)"],
        routings: vec![
            routing(&["k(1)"]),
            routing(&["-k(1)"]),
            reversed(&["k(1)"], &[1]),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

const SUNSET: &[Line] = &[(1, 2, 1), (1, 2, 1), (2, 1, 1)];
const CANONICAL_SUNSET: &[&str] = &["k(1)", "k(2)", "k(1)+k(2)"];

/// Numerators odd and even in the line momenta, and a tensor with an external vector.
const SUNSET_NUMERATORS: &[&str] = &[
    "line(1,1)*line(2,1)",
    "(line(1,1)-line(2,1))*(line(1,1)-line(2,1))",
    "line(1,1)*p(1,1)*line(3,2)*p(1,2)",
];

#[test]
fn two_loop_sunset() {
    assert_routing_invariant(Case {
        loops: 2,
        lines: SUNSET,
        numerators: SUNSET_NUMERATORS,
        routings: vec![
            routing(CANONICAL_SUNSET),
            // One line carrying a single loop momentum read backwards (-k(2)).
            reversed(CANONICAL_SUNSET, &[2]),
            reversed(CANONICAL_SUNSET, &[1, 3]),
            routing(&["k(2)", "k(1)", "k(1)+k(2)"]),
            routing(&["-k(1)", "-k(2)", "-k(1)-k(2)"]),
            // k(1) is not the momentum of any line.
            routing(&["k(1)+k(2)", "k(2)", "k(1)+2*k(2)"]),
            reversed(&["-k(2)", "k(1)+2*k(2)", "k(1)+k(2)"], &[1, 2]),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

#[test]
fn two_loop_sunset_with_a_dotted_line() {
    assert_routing_invariant(Case {
        loops: 2,
        lines: &[(1, 2, 2), (1, 2, 1), (2, 1, 1)],
        numerators: SUNSET_NUMERATORS,
        routings: vec![
            routing(CANONICAL_SUNSET),
            reversed(CANONICAL_SUNSET, &[2]),
            routing(&["k(1)+k(2)", "k(2)", "k(1)+2*k(2)"]),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

/// Two self-loops at one vertex (the sunset with its third line contracted).
#[test]
fn two_loop_self_loops_sharing_a_loop_momentum() {
    assert_routing_invariant(Case {
        loops: 2,
        lines: &[(1, 1, 1), (1, 1, 2)],
        numerators: &[
            "line(1,1)*line(2,1)",
            "line(1,1)*line(1,1)*line(2,2)*line(2,2)",
        ],
        routings: vec![
            routing(&["k(1)", "k(2)"]),
            routing(&["k(1)+k(2)", "k(2)"]),
            reversed(&["-k(1)", "k(1)+k(2)"], &[2]),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

/// `2 line(1).line(2) = D3 - D1 - D2 - muvsq` with `Di = line(i)^2 - muvsq`: the
/// odd sunset numerator equals scalar integrals with contracted lines, which
/// no numerator substitution touches.
#[test]
fn two_loop_sunset_numerator_matches_its_scalar_decomposition() {
    test_utils::run_multi_lane_acceptance(|| {
        let evaluate = |lines: &[Line], numerator: &str| {
            test_utils::get_vakint(Backend::RustRed.settings(2))
                .evaluate(integrand(lines, routing(CANONICAL_SUNSET), numerator).as_view())
                .unwrap()
        };
        let with_powers = |powers: [i64; 3]| {
            let lines = SUNSET
                .iter()
                .zip(powers)
                .map(|(&(from, to, _), power)| (from, to, power))
                .collect::<Vec<_>>();
            evaluate(&lines, "1")
        };
        let numerator = evaluate(SUNSET, "line(1,1)*line(2,1)");
        let decomposition = (with_powers([1, 1, 0])
            - with_powers([0, 1, 1])
            - with_powers([1, 0, 1])
            - vakint_parse!("muvsq").unwrap() * with_powers([1, 1, 1]))
            / Atom::num(2);
        assert!(
            (&numerator - &decomposition).expand().is_zero(),
            "line(1).line(2) sunset = {numerator}\nscalar decomposition = {decomposition}"
        );
    });
}

const TETRAHEDRON: &[Line] = &[
    (1, 2, 1),
    (2, 3, 1),
    (3, 1, 1),
    (1, 4, 1),
    (2, 4, 1),
    (3, 4, 1),
];
const CANONICAL_TETRAHEDRON: &[&str] = &[
    "k(1)",
    "k(2)",
    "k(3)",
    "k(3)-k(1)",
    "k(1)-k(2)",
    "k(2)-k(3)",
];

#[test]
fn three_loop_tetrahedron() {
    assert_routing_invariant(Case {
        loops: 3,
        lines: TETRAHEDRON,
        numerators: &["line(1,1)*line(4,1)", "line(2,1)*p(1,1)*line(6,2)*p(1,2)"],
        routings: vec![
            routing(CANONICAL_TETRAHEDRON),
            reversed(CANONICAL_TETRAHEDRON, &[1, 4]),
            // k(1)+k(2), k(2), k(2)+k(3) for the old k(1), k(2), k(3).
            routing(&[
                "k(1)+k(2)",
                "k(2)",
                "k(2)+k(3)",
                "-k(1)+k(3)",
                "k(1)",
                "-k(3)",
            ]),
            reversed(
                &[
                    "-k(2)",
                    "k(1)+k(2)+k(3)",
                    "k(3)",
                    "k(2)+k(3)",
                    "-k(1)-2*k(2)-k(3)",
                    "k(1)+k(2)",
                ],
                &[2, 5],
            ),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

/// A factorized product: three self-loops at one vertex (the tetrahedron with
/// lines 1, 3 and 6 contracted).
#[test]
fn three_loop_product_of_tadpoles() {
    assert_routing_invariant(Case {
        loops: 3,
        lines: &[(1, 1, 1), (1, 1, 2), (1, 1, 1)],
        numerators: &[
            "line(1,1)*line(2,1)*line(2,2)*line(3,2)",
            "line(1,1)*line(1,1)*line(2,2)*line(3,2)*line(2,3)*line(3,3)",
        ],
        routings: vec![
            routing(&["k(1)", "k(2)", "k(3)"]),
            routing(&["k(1)+k(2)", "k(2)+k(3)", "-k(3)"]),
            reversed(&["k(3)", "k(1)", "k(1)+k(2)"], &[1]),
        ],
        backends: LOWER_LOOP_BACKENDS,
    });
}

/// A four-loop FG family term with a contracted line (the pole study's F7
/// case): `k(3)` of the family routing is not the momentum of any line.
#[test]
fn four_loop_fg_term_with_a_contracted_line() {
    const FAMILY_ROUTING: &[&str] = &[
        "k(1)",
        "k(2)",
        "k(1)-k(3)",
        "k(4)",
        "k(2)-k(3)",
        "k(1)-k(3)+k(4)",
        "k(1)-k(2)",
    ];
    assert_routing_invariant(Case {
        loops: 4,
        lines: &[
            (1, 2, 1),
            (2, 1, 1),
            (3, 4, 1),
            (1, 4, 2),
            (1, 3, 1),
            (4, 1, 1),
            (2, 3, 1),
        ],
        numerators: &[
            "(line(1,1)-line(3,1))*(line(1,1)-line(3,1))-muvsq",
            "line(1,1)*line(3,1)",
        ],
        routings: vec![
            // A spanning-tree routing: every loop momentum is a line momentum.
            reversed(
                &[
                    "k(1)+k(4)",
                    "k(1)",
                    "k(2)+k(4)",
                    "-k(2)+k(3)-k(4)",
                    "k(2)",
                    "k(3)",
                    "k(4)",
                ],
                &[4],
            ),
            routing(FAMILY_ROUTING),
            reversed(FAMILY_ROUTING, &[2, 3, 7]),
        ],
        backends: &[Backend::RustRed],
    });
}

/// A four-loop FG family term with a self-loop whose loop momentum also flows
/// through other lines in the family routing.
#[test]
fn four_loop_fg_term_with_a_self_loop_sharing_its_momentum() {
    const FAMILY_ROUTING: &[&str] = &[
        "k(1)",
        "k(1)-k(3)",
        "k(4)",
        "k(2)-k(3)",
        "k(1)-k(3)+k(4)",
        "k(1)-k(2)",
    ];
    assert_routing_invariant(Case {
        loops: 4,
        lines: &[
            (1, 1, 1),
            (2, 3, 1),
            (1, 3, 2),
            (1, 2, 1),
            (3, 1, 1),
            (1, 2, 2),
        ],
        numerators: &[
            "(line(1,1)-line(6,1))*(line(1,1)-line(6,1))-muvsq",
            "line(1,1)*line(6,1)",
        ],
        routings: vec![
            reversed(
                &[
                    "k(1)",
                    "k(2)+k(4)",
                    "-k(2)+k(3)-k(4)",
                    "k(2)",
                    "k(3)",
                    "k(4)",
                ],
                &[3],
            ),
            routing(FAMILY_ROUTING),
            reversed(FAMILY_ROUTING, &[1, 4]),
        ],
        backends: &[Backend::RustRed],
    });
}
