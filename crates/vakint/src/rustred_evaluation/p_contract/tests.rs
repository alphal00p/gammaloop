//! Hermetic tests of the numerator engine on the input type's one-loop
//! fixtures (gluon bubble, quark bubble): pruning is exact (the pruned inputs
//! are the unpruned ones without their zero-sector keys), the two trace
//! engines agree, the contraction order does not change the inputs, and the
//! one-loop quark bubble equals its closed form.

use symbolica::atom::{Atom, AtomCore};

use super::trace::{IdensoTrace, Letter, NaiveTrace, TraceEngine};
use super::*;
use crate::p_integral::factorized::FactorizedNumerator;

const GLUON: &str = include_str!("../../p_integral/factorized/fixtures/l1-gluon-bubble.txt");
const QUARK: &str = include_str!("../../p_integral/factorized/fixtures/l1-quark-bubble.txt");

fn settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        integral_normalization_factor: crate::LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: 3,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

fn value(text: &str) -> FactorizedNumerator {
    crate::Vakint::initialize_vakint_symbols();
    FactorizedNumerator::from_text(text).expect("fixture")
}

fn view(vakint: &Vakint, inputs: &ContractedInputs) -> GateView {
    inputs.gate_view(vakint, &settings()).unwrap().0
}

/// The prefactor used by the tests: the value's own, couplings kept symbolic.
fn prefactor(value: &FactorizedNumerator) -> Atom {
    value.prefactor.total()
}

#[test]
fn the_denominator_match_maps_every_edge_momentum() {
    let vakint = Vakint::new().unwrap();
    for text in [GLUON, QUARK] {
        let value = value(text);
        let (topology, orientation) = r2_topology(&value).unwrap();
        assert!(topology.to_canonical_string().contains("topo"));
        assert_eq!(orientation.abs(), 1);
        let matched = match_denominators(&vakint, &settings(), &value)
            .unwrap()
            .expect("a one-loop bubble is not scaleless");
        assert_eq!(matched.request_family.loops, 1);
        // Momentum conservation survives the linear map: the external legs
        // carry ±head only.
        for edge in value.edges.iter().filter(|edge| edge.is_external()) {
            let vector = &matched.edge_vectors[&edge.id];
            assert!(vector[0].is_zero());
            assert_eq!(vector[1].clone() * vector[1].clone(), Rational::one());
        }
    }
}

#[test]
fn pruning_drops_exactly_the_zero_sector_keys() {
    let vakint = Vakint::new().unwrap();
    for text in [GLUON, QUARK] {
        let value = value(text);
        let on = contract_and_lower(
            &vakint,
            &settings(),
            &value,
            &prefactor(&value),
            &ContractOptions::default(),
        )
        .unwrap();
        let off = contract_and_lower(
            &vakint,
            &settings(),
            &value,
            &prefactor(&value),
            &ContractOptions {
                prune: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(on.stats.oracle.is_empty(), "{}", on.stats.oracle);
        assert!(!view(&vakint, &on).is_empty());
        let difference = p_gate::compare_gate_views(&view(&vakint, &off), &view(&vakint, &on));
        assert!(difference.is_empty(), "{difference:?}");
        assert!(on.entries() <= off.entries());
    }
}

#[test]
fn the_orders_and_the_trace_engines_agree() {
    let vakint = Vakint::new().unwrap();
    for text in [GLUON, QUARK] {
        let value = value(text);
        let reference = contract_and_lower(
            &vakint,
            &settings(),
            &value,
            &prefactor(&value),
            &ContractOptions::default(),
        )
        .unwrap();
        for options in [
            ContractOptions {
                order: ContractionOrder::Natural,
                ..Default::default()
            },
            ContractOptions {
                trace: TraceBackend::Naive,
                ..Default::default()
            },
        ] {
            let other =
                contract_and_lower(&vakint, &settings(), &value, &prefactor(&value), &options)
                    .unwrap();
            let difference =
                p_gate::compare_gate_views(&view(&vakint, &reference), &view(&vakint, &other));
            assert!(difference.is_empty(), "{options:?}: {difference:?}");
        }
    }
}

/// `Tr(p̸ γ^μ q̸ γ_μ) = 4 (2 - D) p·q` from both engines (a contracted pair
/// inside the word).
#[test]
fn both_trace_engines_contract_inside_a_word() {
    crate::Vakint::initialize_vakint_symbols();
    let dimension = FactorizedNumerator::default_dimension();
    let mu = symbolica::function!(
        symbolica::symbol!(names::MINK),
        Atom::var(dimension),
        symbolica::function!(symbolica::symbol!(names::HEDGE), 1, 1)
    );
    let word = [
        Letter::Slash(1),
        Letter::Index(mu.clone()),
        Letter::Slash(2),
        Letter::Index(mu),
    ];
    let naive = NaiveTrace::default().trace(&word, dimension).unwrap();
    let idenso = IdensoTrace.trace(&word, dimension).unwrap();
    // Both, read by the engine's parser, give the same tensor: compare the
    // atoms after renaming every dummy of a contracted pair to one name.
    let canonical = |atom: &Atom| {
        // idenso writes `Q_e · Q_f` as `g(Q(e, mink(D)), Q(f, mink(D)))`, the
        // recursion as `Q(e, mink(D, x)) Q(f, mink(D, x))`.
        let dot = symbolica::parse!(
            "spenso::g(gammalooprs::Q(e_,spenso::mink(d_)),gammalooprs::Q(f_,spenso::mink(d_)))"
        )
        .to_pattern();
        let pattern = symbolica::parse!("gammalooprs::Q(e_,spenso::mink(d_,x_))").to_pattern();
        atom.expand()
            .replace(dot)
            .with(symbolica::parse!("gammalooprs::P(e_)*gammalooprs::P(f_)").to_pattern())
            .replace(pattern)
            .with(symbolica::parse!("gammalooprs::P(e_)").to_pattern())
            .expand()
    };
    let (a, b) = (canonical(&naive), canonical(&idenso));
    assert!((&a - &b).expand().is_zero(), "naive {a}\nidenso {b}");
}
