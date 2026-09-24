//! Compare scalar values after contracting actual HEP gamma matrices. Four-
//! dimensional metric/epsilon expressions need not have a unique normal form.

use std::{cell::RefCell, collections::HashMap, sync::Once};

use idenso::{
    dirac::GammaSimplifier,
    epsilon::EPSILON_SYMBOL,
    representations::{Bispinor, initialize},
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey, TensorLibrary},
        parsing::{ParseSettings, StrictTensorFilter},
    },
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
    symbolic_parallelism::{SymbolicParallelism, set_symbolica_rayon_enabled},
    tensors::{
        data::DenseTensor,
        parametric::{MixedTensor, ParamTensor},
    },
};
use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
    function,
};

type Library = TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>;

struct TraceEvaluation {
    library: Library,
    momenta: Vec<Atom>,
    scalar_tensors: RefCell<HashMap<Atom, Atom>>,
}

impl TraceEvaluation {
    fn new(sample: usize) -> Self {
        // Independent directions and signed integer components, including
        // timelike and spacelike vectors. Every sample spans all four axes.
        Self::with_components(|index, axis| {
            if index < 4 {
                if index == axis {
                    1 + sample as i64
                } else {
                    sample as i64
                }
            } else {
                ((index * index * 3 + index * (axis + 1) + axis * 2 + sample * (axis + 1)) % 7)
                    as i64
                    - 3
            }
        })
    }

    fn generic(sample: usize) -> Self {
        // Generic signed integer vectors keep the ordering regressions nonzero;
        // the orthogonal first four vectors above can annihilate these traces.
        Self::with_components(|index, axis| {
            ((index * 19
                + axis * 23
                + index * axis * 11
                + sample * (index + axis * 7 + 3)
                + index * index * (axis + 5))
                % 19) as i64
                - 9
        })
    }

    fn with_components(component: impl Fn(usize, usize) -> i64) -> Self {
        static INITIALIZE: Once = Once::new();
        INITIALIZE.call_once(|| {
            initialize();
            // Tiny component contractions do not benefit from waking the
            // machine-wide Rayon pool for each metric or epsilon tensor.
            set_symbolica_rayon_enabled(SymbolicParallelism::Serial);
        });
        let mut library =
            hep_lib_atom::<AbstractIndex, MixedTensor<f64, ExplicitKey<AbstractIndex>>>();
        let mink = Minkowski {}.new_rep(4);
        let mut momenta = Vec::new();
        for index in 0..14 {
            let name = spenso::network::tags::register_tensor_symbol(
                symbolica::wrap_symbol!(format!("trace_validation_p{index}")),
                Vec::new(),
                true,
            )
            .unwrap();
            let key = ExplicitKey::from_iter([mink.to_lib()], name, None);
            let components = (0..4)
                .map(|axis| Atom::num(component(index, axis)))
                .collect();
            library.insert_explicit(key.map_canonical(|structure| {
                MixedTensor::Param(ParamTensor::param(
                    DenseTensor::from_storage_data(components, structure)
                        .unwrap()
                        .into(),
                ))
            }));
            momenta.push(function!(name, mink.to_symbolic([])));
        }
        // Explicit metrics in the input must also stay exact: the mixed
        // library's generic metric factory otherwise supplies f64 constants.
        let key = ExplicitKey::from_iter([mink.to_lib(); 2], ETS.metric, None);
        let components = (0..16)
            .map(|flat| {
                Atom::num(if flat / 4 != flat % 4 {
                    0
                } else if flat == 0 {
                    1
                } else {
                    -1
                })
            })
            .collect();
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(components, structure)
                    .unwrap()
                    .into(),
            ))
        }));
        // Idenso's epsilon includes -i: epsilon(0,1,2,3) = -i for the
        // (+---) metric and HEP library's gamma5 = i gamma0 gamma1 gamma2 gamma3.
        let key = ExplicitKey::from_iter([mink.to_lib(); 4], *EPSILON_SYMBOL, None);
        let components = (0..256)
            .map(|flat| {
                let indices = [flat / 64, flat / 16 % 4, flat / 4 % 4, flat % 4];
                if (0..4).any(|i| indices[i + 1..].contains(&indices[i])) {
                    Atom::Zero
                } else {
                    let inversions = (0..4)
                        .map(|i| indices[i + 1..].iter().filter(|&&j| indices[i] > j).count())
                        .sum::<usize>();
                    if inversions % 2 == 0 {
                        -Atom::i()
                    } else {
                        Atom::i()
                    }
                }
            })
            .collect();
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(components, structure)
                    .unwrap()
                    .into(),
            ))
        }));
        Self {
            library,
            momenta,
            scalar_tensors: RefCell::default(),
        }
    }

    fn trace(&self, length: usize, axial: bool) -> Atom {
        let mut factors: Vec<_> = self.momenta[..length]
            .iter()
            .map(|p| idenso::gamma!(p))
            .collect();
        if axial {
            factors.insert(0, idenso::gamma5!());
        }
        spenso::trace!(&Bispinor {}.new_rep(4); factors)
    }

    fn ordering_slots() -> [Atom; 4] {
        ["order_a", "order_b", "order_c", "order_d"].map(|name| {
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(symbolica::symbol!(name))
                .to_atom()
        })
    }

    fn ordering_factors(&self) -> Vec<Atom> {
        let [a, b, _, _] = Self::ordering_slots();
        // The a pair has four interior gammas; the b pair has five. Choosing
        // the odd interior first avoids a two-word Chisholm expansion.
        vec![
            idenso::gamma!(&a),
            idenso::gamma!(&self.momenta[0]),
            idenso::gamma!(&b),
            idenso::gamma!(&self.momenta[1]),
            idenso::gamma!(&self.momenta[2]),
            idenso::gamma!(&a),
            idenso::gamma!(&self.momenta[3]),
            idenso::gamma!(&self.momenta[4]),
            idenso::gamma!(&b),
            idenso::gamma!(&self.momenta[5]),
            idenso::gamma!(&self.momenta[6]),
            idenso::gamma!(&self.momenta[7]),
        ]
    }

    fn evaluate(&self, expression: &Atom) -> Atom {
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let mut network =
            HepNet::<AbstractIndex>::try_from_view(expression.as_view(), &self.library, &settings)
                .unwrap();
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&self.library, &*FUN_LIB)
            .unwrap();
        let value = match network.result_scalar().unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        };
        assert!(
            matches!(value.as_view(), AtomView::Num(number)
            if matches!(number.get_coeff_view(), CoefficientView::Natural(_, 1, _, 1))),
            "component contraction must produce an exact Gaussian integer, got {value}"
        );
        value
    }

    #[track_caller]
    fn assert_simplification(&self, expression: &Atom, simplified: &Atom, label: &str) -> Atom {
        let original = self.evaluate(expression);
        // Each closed metric/epsilon subnetwork is independent. Contract it
        // once through HEP and reuse its exact value across the polynomial,
        // avoiding thousands of duplicate network constructions.
        let simplified = simplified.replace_map(|atom, _, output| {
            let AtomView::Fun(function) = atom else {
                return;
            };
            if ![ETS.metric, *EPSILON_SYMBOL].contains(&function.get_symbol()) {
                return;
            }
            let mut tensors = self.scalar_tensors.borrow_mut();
            let value = tensors
                .entry(atom.to_owned())
                .or_insert_with(|| self.evaluate(&atom.to_owned()));
            **output = value.clone();
        });
        let simplified = self.evaluate(&simplified);
        assert_eq!(original, simplified, "HEP component mismatch for {label}");
        original
    }
}

impl TraceEvaluation {
    fn assert_trace_length(length: usize) {
        let evaluations = [0, 1, 2].map(Self::new);
        for axial in [false, true] {
            let expression = evaluations[0].trace(length, axial);
            let simplified = expression.simplify_gamma();
            let mut nonzero_samples = 0;
            for (sample, evaluation) in evaluations.iter().enumerate() {
                let label = format!("sample {sample}, length {length}, gamma5 {axial}");
                let value = evaluation.assert_simplification(&expression, &simplified, &label);
                nonzero_samples += usize::from(!value.is_zero());
                if sample == 0 && axial && length == 4 {
                    assert_eq!(
                        value,
                        Atom::num(4) * Atom::i(),
                        "gamma5 orientation in the (+---) basis"
                    );
                }
            }
            if length.is_multiple_of(2) && (!axial || length >= 4) {
                assert!(
                    nonzero_samples > 0,
                    "all samples were trivial at length {length}, gamma5 {axial}"
                );
            }
        }
    }

    fn assert_repeated_indices([left, right]: [usize; 2]) {
        let evaluations = [0, 1, 2].map(Self::new);
        let repeated = Minkowski {}
            .new_rep(4)
            .slot::<AbstractIndex, _>(symbolica::symbol!("trace_repeated_mu"));
        for gamma5_position in [None, Some(0), Some(3), Some(14)] {
            let mut factors: Vec<_> = evaluations[0]
                .momenta
                .iter()
                .map(|p| idenso::gamma!(p))
                .collect();
            factors[left] = idenso::gamma!(repeated);
            factors[right] = idenso::gamma!(repeated);
            if let Some(position) = gamma5_position {
                factors.insert(position, idenso::gamma5!());
            }
            let expression = spenso::trace!(&Bispinor {}.new_rep(4); factors);
            let simplified = expression.simplify_gamma();
            for (sample, evaluation) in evaluations.iter().enumerate() {
                let label = format!(
                    "sample {sample}, repeated indices {left}/{right}, gamma5 at {gamma5_position:?}"
                );
                let _ = evaluation.assert_simplification(&expression, &simplified, &label);
            }
        }
    }

    fn assert_free_axial_trace(length: usize) {
        let evaluations = [0, 1, 2].map(Self::new);
        let slots: Vec<_> = (0..length)
            .map(|index| {
                Minkowski {}
                    .new_rep(4)
                    .slot::<AbstractIndex, _>(symbolica::symbol!(format!("free_axial_mu{index}")))
                    .to_atom()
            })
            .collect();
        let closing_vectors: Vec<_> = evaluations[0]
            .momenta
            .iter()
            .zip(&slots)
            .map(|(momentum, slot)| {
                let AtomView::Fun(momentum) = momentum.as_view() else {
                    unreachable!("test momenta are vector functions")
                };
                function!(momentum.get_symbol(), slot)
            })
            .collect();
        let mut reference = Vec::new();
        let mut nonzero_samples = 0;
        for position in [0, 1, length / 2, length] {
            let mut factors: Vec<_> = slots.iter().map(|slot| idenso::gamma!(slot)).collect();
            factors.insert(position, idenso::gamma5!());
            let trace = spenso::trace!(&Bispinor {}.new_rep(4); factors);
            // Simplify before attaching spectators: distinct explicit slots
            // select the standalone axial shortcut rather than its full pass.
            let simplified = trace.simplify_gamma();
            assert_ne!(trace, simplified);
            assert_eq!(simplified.simplify_gamma(), simplified);
            // Contract each closing vector at its original gamma before taking
            // the matrix trace. HEP's nested trace boundary would otherwise
            // materialize a tensor with 4^length free Lorentz components before
            // seeing outside vectors. Compact slash arguments use the same
            // registered vector/gamma components and no gamma simplification.
            let mut original_factors: Vec<_> = evaluations[0].momenta[..length]
                .iter()
                .map(|momentum| idenso::gamma!(momentum))
                .collect();
            original_factors.insert(position, idenso::gamma5!());
            let original = spenso::trace!(&Bispinor {}.new_rep(4); original_factors);
            for (sample, evaluation) in evaluations.iter().enumerate() {
                let expected = evaluation.evaluate(&original);
                // Each output slot occurs in one metric or epsilon per term.
                // Close those subnetworks with the same explicit vectors and
                // contract actual HEP components independently. This preserves
                // the polynomial's factorization and avoids a rank-length
                // intermediate tensor; no gamma/epsilon identity is used here.
                let closed = simplified.replace_map(|atom, _, output| {
                    let AtomView::Fun(tensor) = atom else {
                        return;
                    };
                    if ![ETS.metric, *EPSILON_SYMBOL].contains(&tensor.get_symbol()) {
                        return;
                    }
                    let mut tensors = evaluation.scalar_tensors.borrow_mut();
                    let value = tensors.entry(atom.to_owned()).or_insert_with(|| {
                        let vectors = tensor.iter().map(|slot| {
                            &closing_vectors[slots
                                .iter()
                                .position(|expected| expected.as_view() == slot)
                                .expect("kernel tensor arguments retain the explicit slots")]
                        });
                        evaluation.evaluate(&(atom * Atom::mul_many(vectors).as_view()))
                    });
                    **output = value.clone();
                });
                assert_eq!(
                    expected,
                    evaluation.evaluate(&closed),
                    "free axial length {length}, gamma5 position {position}, sample {sample}"
                );
                if position == 0 {
                    nonzero_samples += usize::from(!expected.is_zero());
                    reference.push(expected);
                } else {
                    let sign = if position.is_multiple_of(2) { 1 } else { -1 };
                    assert_eq!(expected, Atom::num(sign) * &reference[sample]);
                }
            }
        }
        assert!(
            nonzero_samples > 0,
            "free axial length {length} is nontrivial"
        );
        if length == 4 {
            assert_eq!(reference[0], Atom::num(4) * Atom::i());
        }
    }
}

// Keep each arity/pattern below the integration runner's timeout, without
// changing the sample count or repeating algebra generation for every sample.
macro_rules! trace_cases {
    ($method:ident; $( $name:ident: $case:expr ),* $(,)?) => {
        $( #[test] fn $name() { TraceEvaluation::$method($case); } )*
    };
}

trace_cases!(assert_trace_length;
    gamma_trace_1: 1, gamma_trace_2: 2, gamma_trace_3: 3, gamma_trace_4: 4,
    gamma_trace_5: 5, gamma_trace_6: 6, gamma_trace_7: 7, gamma_trace_8: 8,
    gamma_trace_9: 9, gamma_trace_10: 10, gamma_trace_11: 11, gamma_trace_12: 12,
    gamma_trace_13: 13, gamma_trace_14: 14,
);
trace_cases!(assert_repeated_indices;
    cyclic_contraction: [0, 13], short_chisholm: [2, 5],
    odd_long_chisholm: [0, 6], even_long_chisholm: [0, 7],
);
trace_cases!(assert_free_axial_trace;
    explicit_free_axial_trace_4: 4,
    explicit_free_axial_trace_8: 8,
    explicit_free_axial_trace_12: 12,
);

#[test]
fn repeated_slashes_match_exact_hep_tensor_contractions() {
    for sample in 0..3 {
        let evaluation = TraceEvaluation::new(sample);
        for order in [
            [0, 0, 1, 2, 3, 4, 5, 6],
            [0, 1, 2, 3, 4, 5, 6, 0],
            [0, 1, 0, 2, 3, 4, 5, 6],
        ] {
            for axial in [false, true] {
                let mut factors: Vec<_> = order
                    .iter()
                    .map(|&index| idenso::gamma!(&evaluation.momenta[index]))
                    .collect();
                if axial {
                    factors.insert(3, idenso::gamma5!());
                }
                let expression = spenso::trace!(&Bispinor {}.new_rep(4); factors);
                let label = format!("sample {sample}, repeated slashes {order:?}, gamma5 {axial}");
                let simplified = expression.simplify_gamma();
                let _ = evaluation.assert_simplification(&expression, &simplified, &label);
            }
        }
    }
}

#[test]
fn repeated_pair_order_and_cyclic_cuts_match_exact_hep_contractions() {
    let evaluations = [0, 1, 2].map(TraceEvaluation::generic);
    let bis = Bispinor {}.new_rep(4);
    let factors = evaluations[0].ordering_factors();
    let original = spenso::trace!(&bis; factors.clone());
    let mut cases = [0, 5, 8]
        .map(|cut| {
            let mut rotated = factors.clone();
            rotated.rotate_left(cut);
            spenso::trace!(&bis; rotated)
        })
        .to_vec();
    // Contract b through its odd interior, then a through its odd interior.
    // Each identity contributes -2, leaving four times this eight-gamma word.
    let odd_first_factors =
        [3, 4, 0, 2, 1, 5, 6, 7].map(|index| idenso::gamma!(&evaluations[0].momenta[index]));
    cases.push(Atom::num(4) * spenso::trace!(&bis; odd_first_factors));
    let simplified: Vec<_> = cases.iter().map(GammaSimplifier::simplify_gamma).collect();
    for (sample, evaluation) in evaluations.iter().enumerate() {
        let expected = evaluation.evaluate(&original);
        assert!(
            !expected.is_zero(),
            "ordering sample {sample} must be nonzero"
        );
        for (case, (expression, simplified)) in cases.iter().zip(&simplified).enumerate() {
            let label = format!("pair order/cyclic cut {case}, sample {sample}");
            assert_eq!(
                evaluation.assert_simplification(expression, simplified, &label),
                expected
            );
        }
    }
}

#[test]
fn external_metric_trace_matches_exact_hep_contractions() {
    let evaluations = [0, 1, 2].map(TraceEvaluation::generic);
    let bis = Bispinor {}.new_rep(4);
    let [a, b, c, d] = TraceEvaluation::ordering_slots();
    let repeated = spenso::trace!(&bis; evaluations[0].ordering_factors());
    let mut factors = evaluations[0].ordering_factors();
    factors[5] = idenso::gamma!(&c);
    factors[8] = idenso::gamma!(&d);
    let original = spenso::g!(&a, &c) * spenso::g!(&b, &d) * spenso::trace!(&bis; factors);
    let precontracted = original.schoonschip_with_settings(
        &SchoonschipSettings::default()
            .with_chain_like_functions()
            .without_rank1_tensors(),
    );
    assert_ne!(
        original, precontracted,
        "external metrics must reach the trace body"
    );
    let cases = [original, precontracted];
    let simplified = cases.each_ref().map(GammaSimplifier::simplify_gamma);
    for (sample, evaluation) in evaluations.iter().enumerate() {
        let expected = evaluation.evaluate(&repeated);
        assert!(
            !expected.is_zero(),
            "metric sample {sample} must be nonzero"
        );
        for (case, (expression, simplified)) in cases.iter().zip(&simplified).enumerate() {
            let label = format!("external metric route {case}, sample {sample}");
            assert_eq!(
                evaluation.assert_simplification(expression, simplified, &label),
                expected
            );
        }
    }
}

#[test]
fn axial_pair_order_and_gamma5_positions_match_exact_hep_contractions() {
    let evaluations = [0, 1, 2].map(TraceEvaluation::generic);
    let bis = Bispinor {}.new_rep(4);
    let positions = [0, 3, 6, 12];
    let cases = positions.map(|position| {
        let mut factors = evaluations[0].ordering_factors();
        factors.insert(position, idenso::gamma5!());
        spenso::trace!(&bis; factors)
    });
    let simplified = cases.each_ref().map(GammaSimplifier::simplify_gamma);
    for (sample, evaluation) in evaluations.iter().enumerate() {
        let expected = evaluation.evaluate(&cases[0]);
        assert!(
            !expected.is_zero(),
            "axial ordering sample {sample} must be nonzero"
        );
        for ((position, expression), simplified) in positions.iter().zip(&cases).zip(&simplified) {
            let label = format!("gamma5 position {position}, sample {sample}");
            let value = evaluation.assert_simplification(expression, simplified, &label);
            // Moving gamma5 past each ordinary gamma changes the sign; this
            // also covers cuts with gamma5 inside a candidate contracted pair.
            let sign = if position % 2 == 0 { 1 } else { -1 };
            assert_eq!(value, Atom::num(sign) * &expected, "{label}");
        }
    }
}
