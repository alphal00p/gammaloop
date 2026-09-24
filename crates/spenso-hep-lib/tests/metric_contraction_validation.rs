//! Check metric rewrites by contracting exact HEP tensor components, without
//! assuming that equivalent index expressions have the same normal form.

use std::sync::Once;

use idenso::{
    representations::{Bispinor, ColorFundamental, initialize},
    shorthands::schoonschip::{Schoonschip, SchoonschipSettings},
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey, TensorLibrary},
        parsing::{ParseSettings, StrictTensorFilter},
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
        // Generic mixed-tensor metrics use f64 constants; explicit Atom-valued
        // entries keep this oracle exact even when a metric remains a leaf.
        for (representations, dimension, lorentzian) in [
            ([mink.to_lib(); 2], 4, true),
            ([bis.to_lib(); 2], 4, false),
            ([cof.to_lib(), cof.dual().to_lib()], 3, false),
        ] {
            let key = ExplicitKey::from_iter(representations, ETS.metric, None);
            let components = (0..dimension * dimension)
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
            library.insert_explicit(key.map_canonical(|structure| {
                MixedTensor::Param(ParamTensor::param(
                    DenseTensor::from_storage_data(components, structure)
                        .unwrap()
                        .into(),
                ))
            }));
        }
        let vectors = std::array::from_fn(|index| {
            let name = SPENSO_TAG.rank_one_tensor_symbol(&format!("metric_validation_p{index}"));
            let key = ExplicitKey::from_iter([mink.to_lib()], name, None);
            let components = (0..4)
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
            name
        });
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
        let before = self.evaluate(&expression);
        let after = self.evaluate(&simplified);
        assert_eq!(before, after, "HEP mismatch: {expression} -> {simplified}");
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
        for expression in [
            spenso::g!(&mu, &nu)
                * spenso::g!(&nu, &rho)
                * evaluation.vector(0, &mu)
                * evaluation.vector(1, &rho),
            spenso::g!(&mu, &nu) * spenso::g!(&nu, &rho) * spenso::g!(&rho, &mu),
            spenso::g!(&mu, &mu),
        ] {
            let _ = evaluation.assert_rewrite(expression, &settings, true);
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
    for expression in [
        spenso::g!(i, i.dual()),
        spenso::g!(i, j.dual()) * spenso::g!(j, k.dual()) * spenso::g!(k, i.dual()),
    ] {
        assert_eq!(
            evaluation.assert_rewrite(expression, &settings, true),
            Atom::num(3)
        );
    }
}
