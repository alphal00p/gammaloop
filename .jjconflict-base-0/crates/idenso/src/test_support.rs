use std::sync::LazyLock;

use spenso::{
    s,
    structure::representation::{Minkowski, RepName, Representation},
    symbol_set,
};

use crate::representations::{Bispinor, ColorAdjoint, ColorFundamental};

#[derive(Clone)]
pub(crate) struct TestReps {
    pub mink4: Representation<Minkowski>,
    pub mink_d: Representation<Minkowski>,
    pub bis4: Representation<Bispinor>,
    pub bis_d: Representation<Bispinor>,
    pub cof_nc: Representation<ColorFundamental>,
    pub coad_da: Representation<ColorAdjoint>,
}

// Generate TestSymbols with all alphabet characters and some multi-character symbols
symbol_set!(TestSymbols, TS;
    a b c d e f h i j k l m n o p q r s t u v w x y z
    A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
    mul nu mu alpha beta rho sigma
);

symbol_set!(SpensoTestSymbols, SPENSO_TS, namespace = "spenso";
    l_1 l_2 l_3 l_4 l_5 l_6 l_7 l_8 l_9 l_10 l_20 l_0
    EMRID G dummy_ss l r dim
    ebar
    edge_1_1 edge_2_1 edge_3_1 edge_4_1 edge_5_1 edge_6_1 edge_7_1 edge_8_1 edge_9_1 edge_10_1 edge_11_1 edge_12_1 edge_13_1 edge_14_1 edge_15_1 edge_16_1 edge_17_1 edge_18_1 edge_19_1 edge_20_1
    hedge0 hedge1 hedge2 hedge3 hedge4 hedge5 hedge6 hedge7 hedge8 hedge9 hedge10 hedge11 hedge12 hedge13 hedge14 hedge15 hedge16 hedge17 hedge18 hedge19 hedge20
    hedge_0 hedge_1 hedge_2 hedge_3 hedge_4 hedge_5 hedge_6 hedge_7 hedge_8 hedge_9 hedge_10 hedge_11 hedge_12 hedge_13 hedge_14 hedge_15 hedge_16 hedge_17 hedge_18 hedge_19 hedge_20
);

pub(crate) static TEST_REPS: LazyLock<TestReps> = LazyLock::new(|| {
    TestReps::initialize_symbols();
    TestReps::build()
});

pub(crate) fn test_initialize() -> &'static TestReps {
    &TEST_REPS
}

impl TestReps {
    pub(crate) fn new() -> Self {
        test_initialize().clone()
    }

    fn build() -> Self {
        Self {
            mink4: Minkowski {}.new_rep(4),
            mink_d: Minkowski {}.new_rep(s!(d)),
            bis4: Bispinor {}.new_rep(4),
            bis_d: Bispinor {}.new_rep(s!(d)),
            cof_nc: ColorFundamental {}.new_rep(s!(Nc)),
            coad_da: ColorAdjoint {}.new_rep(s!(dA)),
        }
    }

    fn initialize_symbols() {
        crate::representations::initialize();

        let tags = &spenso::network::tags::SPENSO_TAG;
        let _ = tags.rank_one_tensor_symbol("P");
        let _ = tags.rank_one_tensor_symbol("Q");
        let _ = tags.rank_one_tensor_symbol("K");

        let _ = spenso::p!();
        let _ = spenso::vector_symbol!(q);
        let _ = spenso::vector_symbol!(P);
        let _ = spenso::vector_symbol!(Q);
        let _ = spenso::vector_symbol!(K);

        crate::color::CS.initialize_tensor_symbols();
        crate::dirac::AGS.initialize_tensor_symbols();
        crate::dirac::PS.initialize_tensor_symbols();
        let _ = *crate::epsilon::EPSILON_SYMBOL;

        let _ = TS.A;
        let _ = SPENSO_TS.G;
    }
}

/// Admit a test expression and resolve the shared contraction result without
/// polynomial expansion. Tests requesting a polynomial call `expand` explicitly.
pub(crate) fn contracted_atom(
    expression: symbolica::atom::AtomView<'_>,
) -> Result<symbolica::atom::Atom, crate::tensor::inference::TensorInferenceError> {
    use crate::tensor::{SymbolicTensor, inference::InterfaceInference};
    let interface = InterfaceInference::replacement_interface(expression)?;
    SymbolicTensor::checked_parts(expression.to_owned(), interface)?
        .contract(Default::default())?
        .resolved()
        .map(SymbolicTensor::into_expression)
}

/// Compare factored snapshot arithmetic exactly, keeping every tensor function
/// opaque. Bounds in each independent leaf and the full interpolation grid prove
/// equality over Q; no numerator is expanded and no tensor identity is assumed.
pub(crate) fn assert_factored_snapshot_eq(actual: &str, expected: &str) {
    use symbolica::{
        atom::{Atom, AtomCore, AtomView},
        coefficient::CoefficientView,
        domains::rational::{Q, Rational},
        parser::{ParseSettings, Token},
    };

    fn opaque_leaves(token: &mut Token, leaves: &mut Vec<Token>) {
        match token {
            Token::ID(_) | Token::Fn(_, _, _) => {
                let position = leaves
                    .iter()
                    .position(|leaf| leaf == token)
                    .unwrap_or_else(|| {
                        leaves.push(token.clone());
                        leaves.len() - 1
                    });
                *token = Token::ID(format!("leaf_{position}").into());
            }
            Token::Op(_, _, _, arguments) => {
                for argument in arguments {
                    opaque_leaves(argument, leaves);
                }
            }
            Token::Number(_, false) => {}
            _ => panic!("factored snapshot requires rational arithmetic: {token}"),
        }
    }

    fn degrees(value: AtomView<'_>, parameters: &[Atom]) -> Vec<usize> {
        let mut bound = vec![0usize; parameters.len()];
        match value {
            AtomView::Num(_) => {
                Rational::try_from(value).expect("snapshot coefficient must be rational");
            }
            AtomView::Var(_) => {
                let position = parameters
                    .iter()
                    .position(|p| p.as_view() == value)
                    .unwrap();
                bound[position] = 1;
            }
            AtomView::Add(sum) => {
                for term in sum {
                    for (degree, next) in bound.iter_mut().zip(degrees(term, parameters)) {
                        *degree = (*degree).max(next);
                    }
                }
            }
            AtomView::Mul(product) => {
                for factor in product {
                    for (degree, next) in bound.iter_mut().zip(degrees(factor, parameters)) {
                        *degree = degree.checked_add(next).expect("snapshot degree overflow");
                    }
                }
            }
            AtomView::Pow(power) => {
                let AtomView::Num(exponent) = power.get_exp() else {
                    panic!("snapshot power must have a nonnegative integer exponent");
                };
                let CoefficientView::Natural(exponent, 1, 0, 1) = exponent.get_coeff_view() else {
                    panic!("snapshot power must have a nonnegative integer exponent");
                };
                let exponent = usize::try_from(exponent).expect("negative snapshot exponent");
                bound = degrees(power.get_base(), parameters)
                    .into_iter()
                    .map(|degree| {
                        degree
                            .checked_mul(exponent)
                            .expect("snapshot degree overflow")
                    })
                    .collect();
            }
            _ => unreachable!("function tokens were replaced before atom construction"),
        }
        bound
    }

    // Replace opaque tokens before constructing any Atom, so builtin names and
    // registered tensor normalizers cannot alter either snapshot reference.
    let settings = ParseSettings::default().convert_mul_to_atom(false);
    let mut tokens = [
        Token::parse(actual, settings.clone()).unwrap(),
        Token::parse(expected, settings).unwrap(),
    ];
    let mut leaves = Vec::new();
    for token in &mut tokens {
        opaque_leaves(token, &mut leaves);
    }
    let parse =
        |text: String| Atom::parse(text, "idenso::factored_snapshot", Default::default()).unwrap();
    let expressions = tokens.map(|token| parse(token.to_string()));
    let parameters = (0..leaves.len())
        .map(|position| parse(format!("leaf_{position}")))
        .collect::<Vec<_>>();
    let bounds = degrees(expressions[0].as_view(), &parameters)
        .into_iter()
        .zip(degrees(expressions[1].as_view(), &parameters))
        .map(|(a, b)| a.max(b))
        .collect::<Vec<_>>();
    let count = bounds
        .iter()
        .try_fold(1usize, |count, degree| {
            count.checked_mul(degree.checked_add(1)?)
        })
        .expect("snapshot interpolation grid overflow");
    assert!(
        count <= 65_536,
        "snapshot interpolation grid is too large: {count}"
    );
    let mut evaluator = Atom::evaluator_multiple(&expressions, &parameters)
        .direct_translation(true)
        .horner_iterations(0)
        .cpe_iterations(Some(0))
        .build()
        .unwrap()
        .map_to_ring(&Q)
        .unwrap();
    for point in 0..count {
        let mut remaining = point;
        let values = bounds
            .iter()
            .map(|degree| {
                let coordinate = remaining % (degree + 1);
                remaining /= degree + 1;
                Rational::from(coordinate as i64)
            })
            .collect::<Vec<_>>();
        let mut output = [Rational::from(0), Rational::from(0)];
        evaluator.evaluate_in_ring(&values, &mut output, &Q);
        assert_eq!(
            output[0], output[1],
            "factored snapshots differ at {values:?}\nactual: {actual}\nexpected: {expected}"
        );
    }
}

#[test]
fn factored_snapshot_equality_keeps_opaque_leaves_and_degree_bounds() {
    assert_factored_snapshot_eq("(d-2)*d*g(a,b)", "(d^2-2*d)*g(a,b)");
    assert_factored_snapshot_eq("(f(a)+2*f(b))/2", "f(b)+f(a)/2");
    assert!(std::panic::catch_unwind(|| assert_factored_snapshot_eq("d^2", "d")).is_err());
    assert!(std::panic::catch_unwind(|| assert_factored_snapshot_eq("f(a)", "f(b)")).is_err());
    assert!(std::panic::catch_unwind(|| assert_factored_snapshot_eq("d^(-1)", "d^(-1)")).is_err());
}
