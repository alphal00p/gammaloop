use super::*;
use crate::tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus};
use itertools::Itertools;
use std::collections::HashMap;
use symbolica::coefficient::CoefficientView;

fn color_settings() -> AlgebraSettings {
    AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..Default::default()
    }
}

#[test]
fn fk0032_adjoint_boxes_match_form_and_preserve_spectators() {
    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::adjoint_box_test::hedge");
    let slot = |i| reps.coad_da.to_symbolic([function!(label, i, 1)]);
    let triples = [
        [1, 3, 5],
        [1, 7, 9],
        [3, 19, 23],
        [5, 15, 21],
        [7, 11, 23],
        [9, 13, 17],
        [11, 13, 15],
        [17, 19, 21],
    ];
    let network = Atom::mul_many(triples.map(|[a, b, c]| color_f!(slot(a), slot(b), slot(c))));
    let scalar = parse_lit!((adjoint_box_x + adjoint_box_y) ^ 12);
    // A foreign open port and its scalar metadata survive the local colour
    // decomposition. Its index deliberately uses the old global trace dummy.
    let spectator = function!(
        SPENSO_TAG.tensor_symbol("idenso::adjoint_box_test::spectator"),
        parse_lit!((adjoint_metadata_x + adjoint_metadata_y) ^ 7),
        reps.mink4
            .to_symbolic([Atom::var(symbolica::symbol!("idenso::x"))])
    );
    let source = SymbolicTensor::infer(&network * &scalar * &spectator).unwrap();
    let settings = color_settings();
    let result = source.simplify_algebra(&settings).unwrap();
    // FORM color.h, #call docolor, for these exact ordered vertices:
    // N_A*C_A^4/24 - d44(A,A). The independent Gell-Mann contraction is
    // -108 for SU(3), where d44(A,A)=135 and C_A=3, N_A=8.
    let expected_color = reps.coad_da.dim.to_symbolic() * color_cas!(2, &reps.coad_da).pow(4)
        / Atom::num(24)
        - color_gram!(4, &reps.coad_da, &reps.coad_da);
    let expected = &expected_color * &scalar * &spectator;
    assert!(
        (result.expression.clone() - expected)
            .collect_factors()
            .factor()
            .is_zero()
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    for factor in [&scalar, &spectator] {
        assert!(
            result
                .expression
                .replace(factor.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );
    }
}

#[test]
fn four_adjoint_loop_components_preserve_orientation_and_dummy_scope() {
    test_initialize();
    let label = spenso::index_symbol!("idenso::adjoint_box_components::edge");
    let adjoint = ColorAdjoint {}.new_rep(3);
    let external: [Atom; 4] =
        std::array::from_fn(|i| adjoint.to_symbolic([function!(label, i as i64, 0)]));
    let settings = color_settings();
    // Even and odd permutations, with independently named copies of all
    // internal links. Every one of the 3^4 open components is checked by a
    // finite Levi-Civita contraction, not by another symbolic reduction.
    for (copy, permutation) in [[0, 1, 2], [1, 0, 2], [2, 0, 1], [2, 1, 0]]
        .into_iter()
        .enumerate()
    {
        let internal: [Atom; 4] = std::array::from_fn(|i| {
            adjoint.to_symbolic([function!(label, 10 + i as i64, copy as i64)])
        });
        let first = [&internal[3], &internal[0], &external[0]];
        let mut factors = vec![color_f!(
            first[permutation[0]],
            first[permutation[1]],
            first[permutation[2]]
        )];
        for i in 1..4 {
            factors.push(color_f!(&internal[i - 1], &internal[i], &external[i]));
        }
        let source = SymbolicTensor::infer(Atom::mul_many(factors)).unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_ne!(result.expression, source.expression);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        // Default structural collection may combine f pairs into open chains.
        // Unfold notation for finite components without changing the reduction.
        let explicit = result.undo_chain().unwrap();
        for a in 0..3 {
            for b in 0..3 {
                for c in 0..3 {
                    for d in 0..3 {
                        let values = [a, b, c, d];
                        let mut bindings = external
                            .iter()
                            .cloned()
                            .zip(values)
                            .collect::<HashMap<_, _>>();
                        let actual =
                            Su2Components::evaluate(explicit.expression.as_view(), &mut bindings);
                        let mut expected = 0;
                        for i in 0..3 {
                            for j in 0..3 {
                                for k in 0..3 {
                                    for l in 0..3 {
                                        let first = [l, i, a];
                                        expected += Su2Components::f(
                                            first[permutation[0]],
                                            first[permutation[1]],
                                            first[permutation[2]],
                                        ) * Su2Components::f(i, j, b)
                                            * Su2Components::f(j, k, c)
                                            * Su2Components::f(k, l, d);
                                    }
                                }
                            }
                        }
                        assert!(
                            (actual - expected as f64).abs() < 1e-12,
                            "copy={copy}, component={values:?}: {actual} != {expected}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn adjoint_symmetric_box_trace_contracts_repeated_indices() {
    test_initialize();
    let adjoint = ColorAdjoint {}.new_rep(3);
    let ports = [
        slot!(adjoint, 811).into_atom(),
        slot!(adjoint, 812).into_atom(),
        slot!(adjoint, 813).into_atom(),
    ];
    let source = spenso::trace_sym!(&adjoint; [&ports[0], &ports[0], &ports[1], &ports[2]].map(|port| {
        color_f!(Atom::var(SPENSO_TAG.chain_in), Atom::var(SPENSO_TAG.chain_out), port)
    }));
    let source = SymbolicTensor::infer(source).unwrap();
    let settings = color_settings();
    let result = source.simplify_algebra(&settings).unwrap();
    let expected = Atom::num((5, 6))
        * color_cas!(2, &adjoint).pow(2)
        * function!(
            spenso::network::library::symbolic::ETS.metric,
            &ports[1],
            &ports[2]
        );
    assert_eq!(result.expression, expected);
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    for b in 0..3 {
        for c in 0..3 {
            let mut bindings = [(ports[1].clone(), b), (ports[2].clone(), c)]
                .into_iter()
                .collect();
            let actual = Su2Components::evaluate(result.expression.as_view(), &mut bindings);
            let expected = (0..3)
                .map(|a| Su2Components::symmetric_trace([a, a, b, c]))
                .sum::<f64>();
            assert!((actual - expected).abs() < 1e-12);
        }
    }
}

#[test]
fn adjoint_box_requires_compatible_ports_and_enabled_trace_evaluation() {
    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::adjoint_box_scope::edge");
    let a = ColorAdjoint {}.new_rep(3);
    let ports: [Atom; 8] = std::array::from_fn(|i| a.to_symbolic([function!(label, i as i64, 0)]));
    let network = color_f!(&ports[3], &ports[0], &ports[4])
        * color_f!(&ports[0], &ports[1], &ports[5])
        * color_f!(&ports[1], &ports[2], &ports[6])
        * color_f!(&ports[2], &ports[3], &ports[7]);
    let source = SymbolicTensor::infer(network.clone()).unwrap();
    for color in [
        None,
        Some(ColorSimplifySettings::default().without_trace_evaluation()),
    ] {
        let settings = AlgebraSettings {
            color,
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        assert_eq!(
            source.simplify_algebra(&settings).unwrap().expression,
            source.expression
        );
    }
    // Each f has a valid homogeneous representation. Equal index labels
    // across different dimensions must not be mistaken for connecting edges.
    let incompatible = Atom::mul_many((0..4).map(|i| {
        let representation = ColorAdjoint {}.new_rep(if i % 2 == 0 { 3 } else { 8 });
        let port = |j| representation.to_symbolic([function!(label, j, 0)]);
        color_f!(port((i + 3) % 4), port(i), port(4 + i))
    }));
    let incompatible = SymbolicTensor::infer(incompatible).unwrap();
    assert_eq!(
        incompatible
            .simplify_algebra(&color_settings())
            .unwrap()
            .expression,
        incompatible.expression
    );

    let spin = [slot!(reps.bis4, 701), slot!(reps.bis4, 702)];
    let mu = slot!(reps.mink4, 703);
    let gamma = crate::gamma!(spin[0], spin[1], (mu)) * crate::gamma!(spin[1], spin[0], (mu));
    let foreign = function!(
        SPENSO_TAG.tensor_symbol("idenso::adjoint_box_scope::foreign"),
        a.to_symbolic([Atom::var(symbolica::symbol!("idenso::x"))])
    );
    let tensor = SymbolicTensor::infer(&network * &gamma * &foreign).unwrap();
    let settings = AlgebraSettings {
        contract: AlgebraContraction::None,
        ..color_settings()
    };
    let reduced = tensor.simplify_algebra(&settings).unwrap();
    assert_eq!(reduced.structure, tensor.structure);
    for spectator in [&gamma, &foreign] {
        assert!(
            reduced
                .expression
                .replace(spectator.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );
    }
    assert_eq!(reduced.simplify_algebra(&settings).unwrap(), reduced);
}

/// Numerical SU(2) components with f_abc=epsilon_abc. Only the small open
/// four-vertex identity is evaluated; neither the production graph numerator
/// nor any unrelated scalar sum is distributed by this oracle.
struct Su2Components;

impl Su2Components {
    fn f(a: usize, b: usize, c: usize) -> i64 {
        if a == b || b == c || a == c {
            0
        } else if (a + 1) % 3 == b {
            1
        } else {
            -1
        }
    }

    fn trace([a, b, c, d]: [usize; 4]) -> i64 {
        let mut result = 0;
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        result += Self::f(a, i, j)
                            * Self::f(b, j, k)
                            * Self::f(c, k, l)
                            * Self::f(d, l, i);
                    }
                }
            }
        }
        result
    }

    fn symmetric_trace(args: [usize; 4]) -> f64 {
        args[1..]
            .iter()
            .copied()
            .permutations(3)
            .map(|p| Self::trace([args[0], p[0], p[1], p[2]]) as f64)
            .sum::<f64>()
            / 6.
    }

    fn slots(expression: AtomView<'_>, output: &mut HashMap<Atom, usize>) {
        match expression {
            AtomView::Fun(f) if f.get_symbol() == CS.adjoint_rep && f.get_nargs() == 2 => {
                *output.entry(expression.to_owned()).or_default() += 1;
            }
            AtomView::Fun(f) => {
                for arg in f.iter() {
                    Self::slots(arg, output);
                }
            }
            AtomView::Add(a) => {
                // Each branch owns its internal dummy pairs. Expose only its
                // common external interface to a surrounding product, so an
                // outer minus sign cannot multiply unrelated terms by the
                // dimensions of another branch's dummy sums.
                let mut interface = None;
                for term in a.iter() {
                    let mut branch = HashMap::new();
                    Self::slots(term, &mut branch);
                    assert!(branch.values().all(|count| *count <= 2));
                    branch.retain(|_, count| *count == 1);
                    if let Some(interface) = &interface {
                        assert_eq!(&branch, interface);
                    } else {
                        interface = Some(branch);
                    }
                }
                for (slot, count) in interface.unwrap_or_default() {
                    *output.entry(slot).or_default() += count;
                }
            }
            AtomView::Mul(m) => {
                for arg in m.iter() {
                    Self::slots(arg, output);
                }
            }
            _ => {}
        }
    }

    fn evaluate(expression: AtomView<'_>, bindings: &mut HashMap<Atom, usize>) -> f64 {
        match expression {
            AtomView::Num(number) => {
                let CoefficientView::Natural(n, d, imag, _) = number.get_coeff_view() else {
                    panic!("Unexpected oracle coefficient: {expression}")
                };
                assert_eq!(imag, 0);
                n as f64 / d as f64
            }
            AtomView::Add(add) => add.iter().map(|term| Self::evaluate(term, bindings)).sum(),
            AtomView::Mul(mul) => {
                let mut slots = HashMap::new();
                Self::slots(expression, &mut slots);
                if let Some((dummy, count)) = slots
                    .into_iter()
                    .find(|(slot, _)| !bindings.contains_key(slot))
                {
                    assert_eq!(count, 2);
                    let mut sum = 0.;
                    for value in 0..3 {
                        bindings.insert(dummy.clone(), value);
                        sum += Self::evaluate(expression, bindings);
                    }
                    bindings.remove(&dummy);
                    sum
                } else {
                    mul.iter()
                        .map(|factor| Self::evaluate(factor, bindings))
                        .product()
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::evaluate(base, bindings).powf(Self::evaluate(exponent, bindings))
            }
            AtomView::Fun(f)
                if f.get_symbol() == spenso::network::library::symbolic::ETS.metric =>
            {
                let args = f
                    .iter()
                    .map(|arg| bindings[&arg.to_owned()])
                    .collect::<Vec<_>>();
                if args[0] == args[1] { 1. } else { 0. }
            }
            AtomView::Fun(f) if f.get_symbol() == CS.cas => 2.,
            AtomView::Fun(f) if f.get_symbol() == CS.f => {
                let args = f
                    .iter()
                    .map(|arg| bindings[&arg.to_owned()])
                    .collect::<Vec<_>>();
                Self::f(args[0], args[1], args[2]) as f64
            }
            AtomView::Fun(f) if f.get_symbol() == SPENSO_TAG.trace => {
                let (_, factors) = spenso::shadowing::trace_parts(f).unwrap();
                assert_eq!(factors.len(), 1);
                let AtomView::Fun(sym) = factors[0] else {
                    panic!("Expected symmetric trace")
                };
                assert_eq!(sym.get_symbol(), *spenso::shadowing::SYM);
                let args = sym
                    .iter()
                    .map(|factor| {
                        let AtomView::Fun(t) = factor else {
                            panic!("Expected generator")
                        };
                        assert_eq!(t.get_symbol(), CS.f);
                        bindings[&t.iter().nth(2).unwrap().to_owned()]
                    })
                    .collect::<Vec<_>>();
                assert_eq!(args.len(), 4);
                Self::symmetric_trace(args.try_into().unwrap())
            }
            _ => panic!("Unsupported oracle atom: {expression}"),
        }
    }
}

#[test]
fn six_generator_color_zeros_accept_compound_index_labels() {
    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::color_zero_labels::hedge");
    // Both inputs have the same indexed graph. Only admitted index payloads
    // differ. The existing separated-pair Casimir rule must also cross the
    // cyclic boundary, regardless of the canonical rotation chosen for labels.
    for middle in [[2, 3, 4], [4, 3, 2]] {
        let mut results = Vec::new();
        for compound in [false, true] {
            let slots = [1, 7, 13, 21, 19, 23].map(|i| {
                let index = if compound {
                    function!(label, i, 1)
                } else {
                    Atom::var(symbol!(format!("idenso::color_zero_labels::i{i}")))
                };
                reps.coad_da.to_symbolic([index])
            });
            let [x, a, b, d, e, s] = &slots;
            let word = trace!(&reps.cof_nc;
                [a, s, &slots[middle[0]], &slots[middle[1]], &slots[middle[2]], s]
                    .map(|index| color_t!(index)));
            let expression = color_f!(x, a, d) * color_f!(x, b, e) * word;
            let compact = SymbolicTensor::infer(expression).unwrap();
            let indexed = compact.undo_trace().unwrap().undo_chain().unwrap();
            for source in [compact, indexed] {
                let result = source.simplify_algebra(&color_settings()).unwrap();
                assert!(
                    result.expression.is_zero(),
                    "middle={middle:?}, compound={compound}, source={}, result={}, status={:?}",
                    source.expression,
                    result.expression,
                    result.reduction_status(),
                );
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.simplify_algebra(&color_settings()).unwrap(), result);
                results.push(result);
            }
        }
        assert!(results.windows(2).all(|pair| pair[0] == pair[1]));
    }
}

#[test]
fn color_trace_dummy_cannot_capture_an_external_label() {
    let reps = test_initialize();
    let external = reps
        .coad_da
        .to_symbolic([Atom::var(symbolica::symbol!("idenso::x"))]);
    let expression = trace!(&reps.cof_nc;
        [&external, &reps.coad_da.to_symbolic([Atom::num(98211)]),
         &reps.coad_da.to_symbolic([Atom::num(98212)]), &reps.coad_da.to_symbolic([Atom::num(98213)])]
            .map(|index| color_t!(index)));
    let source = SymbolicTensor::infer(expression).unwrap();
    let reduced = source.simplify_algebra(&color_settings()).unwrap();
    assert!(!reduced.expression.is_zero());
    assert_eq!(reduced.structure, source.structure);
    assert_eq!(
        reduced.simplify_algebra(&color_settings()).unwrap(),
        reduced
    );
}

#[test]
fn color_terminal_metrics_treat_wildcard_spelled_indices_as_literals() {
    let reps = test_initialize();
    let slot = reps
        .coad_da
        .to_symbolic([Atom::var(symbol!("literal_color_index_"))]);
    let source =
        SymbolicTensor::infer(trace!(&reps.cof_nc; [color_t!(&slot), color_t!(&slot)])).unwrap();
    let reduced = source.simplify_algebra(&color_settings()).unwrap();
    assert_eq!(
        reduced.expression,
        reps.coad_da.dim.to_symbolic() * color_idx!(2, &reps.cof_nc)
    );
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
}
