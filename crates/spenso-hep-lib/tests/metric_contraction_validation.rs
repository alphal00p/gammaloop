//! Check slot substitutions and epsilon rewrites using exact HEP components, without
//! assuming that equivalent index expressions have the same normal form.

use std::sync::Once;

use idenso::{
    epsilon::{EPSILON_SYMBOL, EpsilonSimplifier},
    representations::{Bispinor, ColorFundamental, initialize},
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey, TensorLibrary},
        parsing::{AtomStructureExt, ParseSettings, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    symbolic_parallelism::{SymbolicParallelism, set_symbolica_rayon_enabled},
    tensors::{
        data::DenseTensor,
        parametric::{MixedTensor, ParamTensor},
    },
};
use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
    function, symbol,
};

type Library = TensorLibrary<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>;

struct MetricEvaluation {
    library: Library,
    vectors: [Symbol; 4],
    metric_alias: Symbol,
    antisymmetric: Symbol,
    labelled: Symbol,
}

impl MetricEvaluation {
    fn new(sample: usize) -> Self {
        static INITIALIZE: Once = Once::new();
        INITIALIZE.call_once(|| {
            initialize();
            set_symbolica_rayon_enabled(SymbolicParallelism::Serial);
        });
        let mut library =
            hep_lib_atom::<AbstractIndex, MixedTensor<f64, ExplicitKey<AbstractIndex>>>();
        let mink = Minkowski {}.new_rep(4);
        let bis = Bispinor {}.new_rep(4);
        let cof = ColorFundamental {}.new_rep(3);
        // This tensor has the same exact components as g, without its eager
        // normalization callback, so constructor identities have an oracle.
        let metric_alias = spenso::tensor_symbol!("metric_validation_unreduced_metric"; Symmetric);
        // Generic mixed-tensor metrics use f64 constants; explicit Atom-valued
        // entries keep this oracle exact even when a metric remains a leaf.
        for (representations, dimension, lorentzian) in [
            ([mink.to_lib(); 2], 4, true),
            ([bis.to_lib(); 2], 4, false),
            ([cof.to_lib(), cof.dual().to_lib()], 3, false),
        ] {
            let components: Vec<_> = (0..dimension * dimension)
                .map(|flat| {
                    let (row, column) = (flat / dimension, flat % dimension);
                    Atom::num(if row != column {
                        0
                    } else if lorentzian && row > 0 {
                        -1
                    } else {
                        1
                    })
                })
                .collect();
            for name in [ETS.metric, metric_alias] {
                let key = ExplicitKey::from_iter(representations, name, None);
                library.insert_explicit(key.map_canonical(|structure| {
                    MixedTensor::Param(ParamTensor::param(
                        DenseTensor::from_storage_data(components.clone(), structure)
                            .unwrap()
                            .into(),
                    ))
                }));
            }
        }
        let vectors = std::array::from_fn(|index| {
            let name = SPENSO_TAG.rank_one_tensor_symbol(&format!("metric_validation_p{index}"));
            for (representation, dimension) in [
                (mink.to_lib(), 4),
                (cof.to_lib(), 3),
                (cof.dual().to_lib(), 3),
            ] {
                let key = ExplicitKey::from_iter([representation], name, None);
                let components = (0..dimension)
                    .map(|axis| {
                        Atom::num(((index * 3 + axis * 2 + sample * (axis + 1)) % 7) as i64 - 3)
                    })
                    .collect();
                library.insert_explicit(key.map_canonical(|structure| {
                    MixedTensor::Param(ParamTensor::param(
                        DenseTensor::from_storage_data(components, structure)
                            .unwrap()
                            .into(),
                    ))
                }));
            }
            name
        });
        // With HEP's gamma5 and (+---) metric, Idenso uses epsilon(0,1,2,3) = -i.
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
                    if inversions.is_multiple_of(2) {
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
        let antisymmetric =
            spenso::tensor_symbol!("metric_validation_antisymmetric"; Antisymmetric);
        let key = ExplicitKey::from_iter([mink.to_lib(); 2], antisymmetric, None);
        let components = (0..16)
            .map(|flat| Atom::num((flat / 4) as i64 - (flat % 4) as i64))
            .collect();
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(components, structure)
                    .unwrap()
                    .into(),
            ))
        }));
        let labelled = SPENSO_TAG.rank_one_tensor_symbol("metric_validation_labelled");
        // This scalar tensor parameter deliberately equals an index name.
        let key = ExplicitKey::from_iter(
            [mink.to_lib()],
            labelled,
            Some(vec![Atom::var(symbol!("mu"))]),
        );
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(
                    vec![Atom::num(1), Atom::num(2), Atom::num(3), Atom::num(4)],
                    structure,
                )
                .unwrap()
                .into(),
            ))
        }));
        Self {
            library,
            vectors,
            metric_alias,
            antisymmetric,
            labelled,
        }
    }

    fn vector(&self, index: usize, slot: &Atom) -> Atom {
        function!(self.vectors[index], slot)
    }

    fn compact(&self, index: usize) -> Atom {
        self.vector(index, &Minkowski {}.new_rep(4).to_symbolic([]))
    }

    fn slots() -> [Atom; 4] {
        ["mu", "nu", "rho", "sigma"].map(|name| {
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(symbol!(name))
                .to_atom()
        })
    }

    fn evaluate(&self, expression: &Atom) -> Atom {
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let mut network =
            HepNet::<AbstractIndex>::try_from_view(expression.as_view(), &self.library, &settings)
                .unwrap_or_else(|error| panic!("cannot parse {expression}: {error}"));
        network
            .execute::<Sequential, SmallestDegree, _, _, _>(&self.library, &*FUN_LIB)
            .unwrap_or_else(|error| panic!("cannot contract {expression}: {error}"));
        let value = match network.result_scalar().unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        };
        assert!(
            matches!(value.as_view(), AtomView::Num(number)
            if matches!(number.get_coeff_view(), CoefficientView::Natural(_, 1, _, 1))),
            "expected an exact Gaussian integer for {expression}, got {value}"
        );
        value
    }

    #[track_caller]
    fn assert_same_value(&self, expression: &Atom, simplified: &Atom) -> Atom {
        let before = self.evaluate(expression);
        let after = self.evaluate(simplified);
        assert_eq!(before, after, "HEP mismatch: {expression} -> {simplified}");
        before
    }

    #[track_caller]
    fn assert_rewrite(
        &self,
        expression: Atom,
        settings: &SchoonschipSettings,
        must_change: bool,
    ) -> Atom {
        let simplified = expression.schoonschip_with_settings(settings);
        if must_change {
            assert_ne!(
                expression, simplified,
                "the case must exercise a productive rewrite"
            );
        }
        let before = self.assert_same_value(&expression, &simplified);
        assert_eq!(
            simplified.schoonschip_with_settings(settings),
            simplified,
            "rewrite must reach a fixed point"
        );
        before
    }
}

#[test]
fn metric_chains_and_traces_preserve_exact_components() {
    let settings = SchoonschipSettings::default().without_rank1_tensors();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, rho, _] = MetricEvaluation::slots();
        for (expression, must_change) in [
            (
                spenso::g!(&mu, &nu)
                    * spenso::g!(&nu, &rho)
                    * evaluation.vector(0, &mu)
                    * evaluation.vector(1, &rho),
                true,
            ),
            (
                spenso::g!(&mu, &nu) * spenso::g!(&nu, &rho) * spenso::g!(&rho, &mu),
                true,
            ),
            // Self-traces are already scalar when the metric is constructed.
            (spenso::g!(&mu, &mu), false),
        ] {
            let _ = evaluation.assert_rewrite(expression, &settings, must_change);
        }
    }
}

#[test]
fn long_product_contractions_preserve_independent_sum_scopes() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let slots: Vec<_> = (0..65)
            .map(|index| {
                Minkowski {}
                    .new_rep(4)
                    .slot::<AbstractIndex, _>(symbol!(&format!("local_chain_{index}")))
                    .to_atom()
            })
            .collect();
        for length in [8, 16, 32, 64] {
            let chain: Atom = (0..length)
                .map(|index| spenso::g!(&slots[index], &slots[index + 1]))
                .product();
            // Reusing dummy labels between independent scalar terms must not join
            // their networks. The sum remains factored inside the outer power.
            let first =
                &chain * evaluation.vector(0, &slots[0]) * evaluation.vector(1, &slots[length]);
            let second =
                &chain * evaluation.vector(2, &slots[0]) * evaluation.vector(3, &slots[length]);
            let expression = (first + second * 2 + 101).pow(2);
            let closed: Atom = (0..length)
                .map(|index| spenso::g!(&slots[index], &slots[(index + 1) % length]))
                .product();
            for settings in [
                SchoonschipSettings::default().without_rank1_tensors(),
                SchoonschipSettings::default(),
            ] {
                let value = evaluation.assert_rewrite(expression.clone(), &settings, true);
                assert!(
                    !value.is_zero(),
                    "the component comparison must be nontrivial"
                );
                assert_eq!(
                    evaluation.assert_rewrite(closed.clone(), &settings, true),
                    Atom::num(4),
                    "closed metric chain of length {length} must give the dimension"
                );
            }
        }
    }
}

#[test]
fn metric_substitution_preserves_antisymmetric_orientation_and_scalar_parameters() {
    let settings = SchoonschipSettings::default().without_rank1_tensors();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, rho, sigma] = MetricEvaluation::slots();
        for (source, target) in [(&mu, &rho), (&sigma, &mu)] {
            // Close the free slots only after rewriting, so the metric must
            // act on the antisymmetric tensor rather than a spectator vector.
            let expression =
                spenso::g!(source, target) * function!(evaluation.antisymmetric, source, &nu);
            let simplified = expression.schoonschip_with_settings(&settings);
            assert_ne!(expression, simplified);
            let spectator = evaluation.vector(0, target) * evaluation.vector(1, &nu);
            let before = evaluation.evaluate(&(&expression * &spectator));
            let after = evaluation.evaluate(&(&simplified * &spectator));
            assert_eq!(
                before, after,
                "orientation mismatch: {expression} -> {simplified}"
            );
            assert!(
                !before.is_zero(),
                "orientation test must not be a zero identity"
            );
        }
        let expression = spenso::g!(&mu, &nu) * function!(evaluation.labelled, symbol!("mu"), &mu);
        let simplified = expression.schoonschip_with_settings(&settings);
        assert_ne!(expression, simplified);
        let spectator = evaluation.vector(1, &nu);
        assert_eq!(
            evaluation.evaluate(&(&expression * &spectator)),
            evaluation.evaluate(&(&simplified * &spectator)),
            "scalar parameter changed: {expression} -> {simplified}"
        );
    }
}

#[test]
fn metric_rewrites_preserve_sum_and_power_scopes() {
    let settings = SchoonschipSettings::default().without_rank1_tensors();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, rho, _] = MetricEvaluation::slots();
        let scalar = spenso::g!(&mu, &nu) * evaluation.vector(0, &mu) * evaluation.vector(1, &nu);
        for (expression, must_change) in [
            (
                spenso::g!(&mu, &nu)
                    * (evaluation.vector(0, &mu) + evaluation.vector(1, &mu))
                    * evaluation.vector(2, &nu),
                false,
            ),
            (
                &scalar
                    + spenso::g!(&mu, &rho)
                        * evaluation.vector(2, &mu)
                        * evaluation.vector(3, &rho),
                true,
            ),
            (
                spenso::g!(evaluation.compact(0), evaluation.compact(1)).pow(2),
                false,
            ),
            ((&scalar + 1).pow(2), true),
            (evaluation.vector(0, &mu).pow(2), true),
            (spenso::g!(&mu, &nu).pow(2), true),
        ] {
            let _ = evaluation.assert_rewrite(expression, &settings, must_change);
        }
    }
}

#[test]
fn compact_vectors_and_chain_like_metric_endpoints_preserve_exact_components() {
    let settings = SchoonschipSettings::default().with_chain_like_functions();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, _, _] = MetricEvaluation::slots();
        let bis = Bispinor {}.new_rep(4);
        let [i, j, k] = ["i", "j", "k"].map(|name| bis.slot::<AbstractIndex, _>(symbol!(name)));
        for expression in [
            spenso::g!(&mu, evaluation.compact(0)) * spenso::g!(&mu, evaluation.compact(1)),
            spenso::g!(&mu, &nu)
                * spenso::trace!(
                    &bis,
                    idenso::gamma!(&mu),
                    idenso::gamma!(evaluation.compact(0)),
                    idenso::gamma!(&nu),
                    idenso::gamma!(evaluation.compact(1))
                ),
            spenso::g!(i, j)
                * spenso::chain!(
                    j,
                    k,
                    idenso::gamma!(evaluation.compact(0)),
                    idenso::gamma!(evaluation.compact(1))
                )
                * spenso::g!(k, i),
        ] {
            let _ = evaluation.assert_rewrite(expression, &settings, true);
        }
    }
}

#[test]
fn dual_metric_loops_preserve_exact_components() {
    let evaluation = MetricEvaluation::new(0);
    let cof = ColorFundamental {}.new_rep(3);
    let [i, j, k] =
        ["color_i", "color_j", "color_k"].map(|name| cof.slot::<AbstractIndex, _>(symbol!(name)));
    let settings = SchoonschipSettings::default().without_rank1_tensors();
    for (expression, must_change) in [
        (spenso::g!(i, i.dual()), false),
        (
            spenso::g!(i, j.dual()) * spenso::g!(j, k.dual()) * spenso::g!(k, i.dual()),
            true,
        ),
    ] {
        assert_eq!(
            evaluation.assert_rewrite(expression, &settings, must_change),
            Atom::num(3)
        );
    }
    let slots: Vec<_> = (0..64)
        .map(|index| cof.slot::<AbstractIndex, _>(symbol!(&format!("dual_loop_{index}"))))
        .collect();
    for length in [8, 16, 32, 64] {
        for reverse in [false, true] {
            let expression: Atom = (0..length)
                .map(|index| {
                    let next = (index + 1) % length;
                    if reverse {
                        spenso::g!(slots[index].dual(), slots[next])
                    } else {
                        spenso::g!(slots[index], slots[next].dual())
                    }
                })
                .product();
            assert_eq!(
                evaluation.assert_rewrite(expression, &settings, true),
                Atom::num(3),
                "dual loop of length {length}, reversed={reverse} must give the dimension"
            );
        }
    }
}

#[test]
fn constructed_metric_traces_are_scalar_before_taking_powers() {
    let evaluation = MetricEvaluation::new(0);
    let [mu, _, _, _] = MetricEvaluation::slots();
    let cof = ColorFundamental {}.new_rep(3);
    let color = cof.slot::<AbstractIndex, _>(symbol!("trace_color"));
    let lorentz_trace = function!(evaluation.metric_alias, &mu, &mu);
    let color_trace = function!(
        evaluation.metric_alias,
        color.to_atom(),
        color.dual().to_atom()
    );
    // The alias retains a genuine tensor-network trace. In particular, the
    // square of the scalar g(mu,mu)=D is D², not the old power-before-trace D.
    for (source, oracle, expected) in [
        (
            "g(mink(4,mu),mink(4,mu))^2",
            lorentz_trace.clone().pow(2),
            16,
        ),
        (
            "g(mink(4,mu),mink(4,mu))^3",
            lorentz_trace.clone().pow(3),
            64,
        ),
        (
            "g(cof(3,trace_color),dind(cof(3,trace_color)))^2",
            color_trace.clone().pow(2),
            9,
        ),
        (
            "(-3*g(mink(4,mu),mink(4,mu))+2)^2",
            (&lorentz_trace * -3 + 2).pow(2),
            100,
        ),
        (
            "(g(mink(4,mu),mink(4,mu))+2)^2*(g(cof(3,trace_color),dind(cof(3,trace_color)))+1)",
            (&lorentz_trace + 2).pow(2) * (&color_trace + 1),
            144,
        ),
    ] {
        let constructed = Atom::parse(
            source,
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap();
        assert_eq!(constructed, Atom::num(expected), "source: {source}");
        assert_eq!(constructed.normalize_dots(), constructed);
        assert_eq!(
            evaluation.assert_same_value(&oracle, &constructed),
            Atom::num(expected)
        );
    }
}

#[test]
fn constructed_compact_metric_agrees_with_explicit_component_contraction() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, _, _] = MetricEvaluation::slots();
        for coefficient in [1, -3] {
            let source = format!(
                "g({},{coefficient}*{})",
                mu.to_plain_string(),
                evaluation.compact(0).to_plain_string()
            );
            let constructed = Atom::parse(
                &source,
                "spenso",
                symbolica::parser::ParseSettings::symbolica(),
            )
            .unwrap();
            assert_eq!(constructed, evaluation.vector(0, &mu) * coefficient);
            assert_eq!(constructed.normalize_dots(), constructed);
            let explicit = function!(evaluation.metric_alias, &mu, &nu)
                * evaluation.vector(0, &nu)
                * evaluation.vector(1, &mu)
                * coefficient;
            assert_eq!(
                evaluation.assert_same_value(&explicit, &(constructed * evaluation.vector(1, &mu))),
                Atom::num(coefficient * [8, 1, 1][sample])
            );
        }
    }
}

#[test]
fn compact_metric_and_tagged_vector_substitutions_agree() {
    let settings = SchoonschipSettings::default();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, _, _] = MetricEvaluation::slots();
        let tensor = function!(evaluation.antisymmetric, &mu, &nu);
        let spectator = evaluation.vector(1, &nu);
        let expected = evaluation.evaluate(
            &(function!(evaluation.antisymmetric, evaluation.compact(0), &nu) * &spectator),
        );
        assert!(!expected.is_zero());
        for source in [
            spenso::g!(&mu, evaluation.compact(0)),
            evaluation.vector(0, &mu),
        ] {
            let expression = source * &tensor;
            let simplified = expression.schoonschip_with_settings(&settings);
            assert_ne!(expression, simplified);
            assert!(!simplified.has_repeated_explicit_indices());
            assert_eq!(
                evaluation
                    .assert_same_value(&(&expression * &spectator), &(&simplified * &spectator)),
                expected
            );
        }
        let explicit = evaluation.vector(0, &mu) * &tensor;
        assert_eq!(
            explicit.schoonschip_with_settings(
                &SchoonschipSettings::default().without_rank1_tensors(),
            ),
            explicit,
            "disabling rank-one substitutions must preserve the explicit vector"
        );
    }
}

#[test]
fn tagged_vectors_and_compact_metrics_preserve_epsilon_orientation() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, rho, sigma] = MetricEvaluation::slots();
        let tensor = idenso::epsilon!(&mu, &nu, &rho, &sigma);
        let spectator =
            evaluation.vector(1, &nu) * evaluation.vector(2, &rho) * evaluation.vector(3, &sigma);
        for source in [
            spenso::g!(&mu, evaluation.compact(0)),
            evaluation.vector(0, &mu),
        ] {
            // Close the other slots after substitution, keeping its target
            // and the orientation of the epsilon unambiguous.
            let expression = source * &tensor;
            let simplified = expression.schoonschip();
            assert_ne!(expression, simplified);
            assert!(!simplified.has_repeated_explicit_indices());
            let value = evaluation
                .assert_same_value(&(&expression * &spectator), &(&simplified * &spectator));
            assert!(!value.is_zero());
        }
    }
}

#[test]
fn tagged_vectors_contract_inside_ordered_chains_and_trace_projectors() {
    let settings = SchoonschipSettings::default().with_chain_like_functions();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, _, _, _] = MetricEvaluation::slots();
        let bis = Bispinor {}.new_rep(4);
        let [i, j] =
            ["vector_i", "vector_j"].map(|name| bis.slot::<AbstractIndex, _>(symbol!(name)));
        let gamma_mu = idenso::gamma!(&mu);
        let gamma_q = idenso::gamma!(evaluation.compact(1));
        for word in [
            spenso::trace!(&bis, &gamma_mu, &gamma_q),
            spenso::trace_sym!(&bis, &gamma_mu, &gamma_q),
            spenso::g!(i, j) * spenso::chain!(i, j, &gamma_mu, &gamma_q),
        ] {
            for source in [
                spenso::g!(&mu, evaluation.compact(0)),
                evaluation.vector(0, &mu),
            ] {
                let value = evaluation.assert_rewrite(source * &word, &settings, true);
                // Tr(/p /q) = 4 p.q, also for the normalized symmetric trace.
                assert_eq!(value, Atom::num(4 * [8, 1, 1][sample]));
            }
        }
    }
}

#[test]
#[ignore = "Known compact dual-dot variance loss predates shared slot contraction; see gamma_cleanup.json"]
fn tagged_vector_dots_preserve_dual_slot_orientation() {
    // Both old and new cleanup turn V(cof(i))*W(dind(cof(i))) into
    // g(V(cof),W(cof)); the compact parser requires opposite orientations.
    let settings = SchoonschipSettings::default();
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let cof = ColorFundamental {}.new_rep(3);
        let index = cof.slot::<AbstractIndex, _>(symbol!("vector_color"));
        let (index, dual) = (index.to_atom(), index.dual().to_atom());
        for (left, right) in [(&index, &dual), (&dual, &index)] {
            let expression = evaluation.vector(0, left) * evaluation.vector(1, right);
            assert_eq!(
                evaluation.assert_rewrite(expression, &settings, true),
                Atom::num(-5)
            );
        }
    }
}

#[test]
fn compact_epsilon_powers_simplify_without_repeated_explicit_indices() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let epsilon = idenso::epsilon!(; (0..4).map(|i| evaluation.compact(i)));
        for expression in [
            epsilon.clone().pow(2),
            epsilon.clone().pow(3),
            (epsilon.clone().pow(2) + 1).pow(2),
        ] {
            assert!(!expression.has_repeated_explicit_indices());
            let simplified = expression.simplify_epsilon();
            assert_ne!(expression, simplified);
            assert_eq!(simplified.simplify_epsilon(), simplified);
            let value = evaluation.assert_same_value(&expression, &simplified);
            assert!(!value.is_zero());
        }
    }
}

#[test]
fn disjoint_epsilon_pairs_simplify_before_closing_their_indices() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let left = MetricEvaluation::slots();
        let right = ["epsilon_a", "epsilon_b", "epsilon_c", "epsilon_d"].map(|name| {
            Minkowski {}
                .new_rep(4)
                .slot::<AbstractIndex, _>(symbol!(name))
                .to_atom()
        });
        let expression = idenso::epsilon!(; &left) * idenso::epsilon!(; &right);
        assert!(!expression.has_repeated_explicit_indices());
        let simplified = expression.simplify_epsilon();
        assert_ne!(expression, simplified);
        assert_eq!(simplified.simplify_epsilon(), simplified);
        let spectator = (0..4).fold(Atom::one(), |product, i| {
            product * evaluation.vector(i, &left[i]) * evaluation.vector(i, &right[i])
        });
        let value =
            evaluation.assert_same_value(&(&expression * &spectator), &(&simplified * &spectator));
        assert!(!value.is_zero());
    }
}

#[test]
fn metric_induced_epsilon_zero_preserves_an_unaffected_summand() {
    for sample in 0..3 {
        let evaluation = MetricEvaluation::new(sample);
        let [mu, nu, rho, sigma] = MetricEvaluation::slots();
        let zero = spenso::g!(&mu, &nu)
            * idenso::epsilon!(&mu, &nu, &rho, &sigma)
            * evaluation.vector(2, &rho)
            * evaluation.vector(3, &sigma);
        let unaffected = spenso::g!(evaluation.compact(0), evaluation.compact(1));
        let expression = zero + &unaffected;
        let simplified = expression.simplify_epsilon();
        assert_eq!(simplified, unaffected);
        assert_eq!(
            evaluation.assert_same_value(&expression, &simplified),
            Atom::num([8, 1, 1][sample])
        );
    }
}
