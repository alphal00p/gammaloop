//! Colour results must not depend on the order in which labels were created.
//!
//! Symbolica stores the arguments of the antisymmetric structure constant by
//! symbol order, and the shared contractor writes adjoint chains from the
//! stored order. The same network therefore reaches the colour rules in
//! different chain notations, depending on whether its labels or the pool
//! dummies `d_n` were created first, and on the order of the labels.

use itertools::Itertools;
use spenso::network::library::symbolic::ETS;
use symbolica::id::Replacement;

use super::*;
use crate::tensor::AlgebraSettings;

/// SU(3) colour slots of the labels `scope::name`.
struct Labels(String);

impl Labels {
    fn label(&self, name: &str) -> Atom {
        Atom::var(symbol!(&format!("{}::{name}", self.0)))
    }

    fn adjoint(&self, name: &str) -> Atom {
        ColorAdjoint {}.to_symbolic([Atom::num(8), self.label(name)])
    }

    fn structure(&self, [a, b, c]: [&str; 3]) -> Atom {
        color_f!(self.adjoint(a), self.adjoint(b), self.adjoint(c))
    }

    /// chain(s,t,F^a1..F^an) = f(s,i1,a1) f(i1,i2,a2) .. f(i_{n-1},t,an).
    fn adjoint_chain(&self, start: &str, end: &str, generators: &[&str]) -> Atom {
        chain!(self.adjoint(start), self.adjoint(end); generators.iter().map(|generator| {
            color_f!(
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out),
                self.adjoint(generator)
            )
        }))
    }

    fn fundamental_chain(&self, start: &str, end: &str, generators: &[&str]) -> Atom {
        chain!(
            ColorFundamental {}.to_symbolic([Atom::num(3), self.label(start)]),
            spenso::dind!(ColorFundamental {}.to_symbolic([Atom::num(3), self.label(end)]));
            generators.iter().map(|generator| color_t!(self.adjoint(generator)))
        )
    }

    fn fundamental_trace(&self, generators: &[&str]) -> Atom {
        trace!(ColorFundamental {}.to_symbolic([Atom::num(3)]);
            generators.iter().map(|generator| color_t!(self.adjoint(generator))))
    }
}

struct Case {
    name: &'static str,
    labels: &'static [&'static str],
    build: fn(&Labels) -> Atom,
    expected: fn(&Labels) -> Atom,
    /// The contractor writes an open adjoint network as chains chosen by
    /// label order, so only the value of the result is order independent.
    open: bool,
}

fn adjoint_casimir() -> Atom {
    parse_lit!(cas(2, coad(8)), default_namespace = "spenso")
}

fn cases() -> Vec<Case> {
    vec![
        // A structure constant closing an adjoint chain at both ends.
        Case {
            name: "closed_chain_ends",
            labels: &["A", "B", "C", "D", "X"],
            build: |l| l.adjoint_chain("A", "B", &["C", "D"]) * l.structure(["X", "A", "B"]),
            expected: |l| -adjoint_casimir() / 2 * l.structure(["C", "D", "X"]),
            open: false,
        },
        // The same chain, also closed by a second structure constant.
        Case {
            name: "closed_chain_ends_and_generators",
            labels: &["A", "B", "C", "D", "X", "Y"],
            build: |l| {
                l.adjoint_chain("A", "B", &["C", "D"])
                    * l.structure(["D", "C", "Y"])
                    * l.structure(["X", "A", "B"])
            },
            expected: |l| {
                adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // Its explicit structure constants, chain(A,B,F^C,F^D) = f(A,E,C) f(E,B,D).
        Case {
            name: "explicit_structure_constants",
            labels: &["A", "B", "C", "D", "E", "X", "Y"],
            build: |l| {
                l.structure(["A", "E", "C"])
                    * l.structure(["E", "B", "D"])
                    * l.structure(["D", "C", "Y"])
                    * l.structure(["X", "A", "B"])
            },
            expected: |l| {
                adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // A chain from an outer label through two triangles,
        // chain(Y,B,F^G,F^E) = f(Y,i,G) f(i,B,E).
        Case {
            name: "chain_through_triangles",
            labels: &["A", "B", "E", "G", "X", "Y"],
            build: |l| {
                l.adjoint_chain("Y", "B", &["G", "E"])
                    * l.structure(["E", "G", "A"])
                    * l.structure(["X", "A", "B"])
            },
            expected: |l| {
                -adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // A chain whose start is a distant generator label of the same chain:
        // f(A,E,C) f(E,G,D) f(G,X,A) is a triangle.
        Case {
            name: "chain_start_as_generator",
            labels: &["A", "C", "D", "X"],
            build: |l| l.adjoint_chain("A", "X", &["C", "D", "A"]),
            expected: |l| -adjoint_casimir() / 2 * l.structure(["C", "D", "X"]),
            open: false,
        },
        Case {
            name: "chain_start_as_generator_closed",
            labels: &["A", "C", "D", "X", "Y"],
            build: |l| l.adjoint_chain("A", "X", &["C", "D", "A"]) * l.structure(["C", "D", "Y"]),
            expected: |l| {
                -adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // Repeated distant generator labels and the start as a generator.
        Case {
            name: "chain_repeated_generators",
            labels: &["D", "E", "X", "Y"],
            build: |l| l.adjoint_chain("E", "Y", &["D", "X", "E", "D"]),
            expected: |l| {
                adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // Two chains, each with an endpoint among the other's generators.
        Case {
            name: "chains_closing_each_other",
            labels: &["D", "E", "G", "X", "Y"],
            build: |l| {
                l.adjoint_chain("E", "Y", &["D", "G"]) * l.adjoint_chain("X", "G", &["D", "E"])
            },
            expected: |l| {
                -adjoint_casimir().pow(2) / 2 * ETS.metric_literal(l.adjoint("X"), l.adjoint("Y"))
            },
            open: false,
        },
        // f^{aij} f^{bjk} f^{cki} = C_A/2 f^{abc}, closed by f^{abc}.
        Case {
            name: "tetrahedron",
            labels: &["a", "b", "c", "i", "j", "k"],
            build: |l| {
                l.structure(["a", "i", "j"])
                    * l.structure(["b", "j", "k"])
                    * l.structure(["c", "k", "i"])
                    * l.structure(["a", "b", "c"])
            },
            expected: |_| Atom::num(8) * adjoint_casimir().pow(2) / 2,
            open: false,
        },
        // f^{ECD} f^{XCD} = C_A delta^{EX}.
        Case {
            name: "structure_constant_loop_on_tree",
            labels: &["A", "B", "C", "D", "E", "X"],
            build: |l| {
                l.structure(["A", "B", "E"])
                    * l.structure(["E", "C", "D"])
                    * l.structure(["C", "D", "X"])
            },
            expected: |l| adjoint_casimir() * l.structure(["A", "B", "X"]),
            open: false,
        },
        // T^a T^b T^a = (C_F - C_A/2) T^b.
        Case {
            name: "fundamental_line",
            labels: &["i", "j", "a", "b"],
            build: |l| l.fundamental_chain("i", "j", &["a", "b", "a"]),
            expected: |l| {
                let line = l.fundamental_chain("i", "j", &["b"]);
                parse_lit!(cas(2, cof(3)), default_namespace = "spenso") * &line
                    - adjoint_casimir() / 2 * line
            },
            open: false,
        },
        // T^A T^B f^{ABC} = i/2 C_A T^C, with f^{ABC} = chain(A,B,F^C).
        Case {
            name: "fundamental_line_on_adjoint_chain",
            labels: &["i", "j", "A", "B", "C"],
            build: |l| {
                l.fundamental_chain("i", "j", &["A", "B"]) * l.adjoint_chain("A", "B", &["C"])
            },
            expected: |l| Atom::i() / 2 * adjoint_casimir() * l.fundamental_chain("i", "j", &["C"]),
            open: false,
        },
        // Open networks: nothing reduces.
        Case {
            name: "open_chain",
            labels: &["A", "B", "C", "D"],
            build: |l| l.adjoint_chain("A", "B", &["C", "D"]),
            expected: |l| l.adjoint_chain("A", "B", &["C", "D"]),
            open: false,
        },
        Case {
            name: "open_tree",
            labels: &["A", "B", "C", "D", "E"],
            build: |l| l.structure(["A", "B", "E"]) * l.structure(["E", "C", "D"]),
            expected: |l| l.structure(["A", "B", "E"]) * l.structure(["E", "C", "D"]),
            open: true,
        },
        Case {
            name: "open_tree_and_line",
            labels: &["i", "j", "A", "B", "C"],
            build: |l| l.fundamental_chain("i", "j", &["A"]) * l.structure(["A", "B", "C"]),
            expected: |l| l.fundamental_chain("i", "j", &["A"]) * l.structure(["A", "B", "C"]),
            open: false,
        },
        // Tr(T^a T^b T^c) f^{abc} = i/2 T_F C_A N_A.
        Case {
            name: "fundamental_trace",
            labels: &["a", "b", "c"],
            build: |l| l.fundamental_trace(&["a", "b", "c"]) * l.structure(["a", "b", "c"]),
            expected: |_| {
                Atom::i() / 2
                    * Atom::num(8)
                    * parse_lit!(idx(2, cof(3)), default_namespace = "spenso")
                    * adjoint_casimir()
            },
            open: false,
        },
    ]
}

/// The value of a colour expression, with its lines written out and its
/// dummies labelled canonically.
fn canonical(expression: Atom) -> Atom {
    SymbolicTensor::infer(expression)
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .canonize()
        .unwrap()
        .into_expression()
}

/// Every label order of each case, once with the labels created before the
/// pool dummies and once after them, simplifies to the same result.
#[test]
fn colour_results_are_independent_of_label_order() {
    test_initialize();
    let cases = cases();
    let scopes = |position: &str| {
        cases
            .iter()
            .map(|case| {
                case.labels
                    .iter()
                    .permutations(case.labels.len())
                    .enumerate()
                    .map(|(order, labels)| {
                        let scope = format!("label_order_{position}_{}_{order}", case.name);
                        for label in labels {
                            symbol!(&format!("{scope}::{label}"));
                        }
                        Labels(scope)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let before_pool = scopes("before");
    for index in 0..32 {
        let _ = Atom::from(AbstractIndex::Dummy(index));
    }
    let after_pool = scopes("after");
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..Default::default()
    };
    let common = Labels("label_order".into());
    for (index, case) in cases.iter().enumerate() {
        let expected = (case.expected)(&common);
        for labels in before_pool[index].iter().chain(&after_pool[index]) {
            let result = SymbolicTensor::infer((case.build)(labels))
                .unwrap()
                .simplify_algebra(&settings)
                .unwrap()
                .into_expression();
            let relabelled = result.replace_multiple(
                case.labels
                    .iter()
                    .map(|label| {
                        Replacement::new(labels.label(label).to_pattern(), common.label(label))
                    })
                    .collect::<Vec<_>>(),
            );
            if case.open {
                assert_eq!(
                    canonical(relabelled),
                    canonical(expected.clone()),
                    "{} in {}",
                    case.name,
                    labels.0
                );
            } else {
                assert_eq!(relabelled, expected, "{} in {}", case.name, labels.0);
            }
        }
    }
}
