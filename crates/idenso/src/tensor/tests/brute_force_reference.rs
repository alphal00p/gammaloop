use std::collections::BTreeMap;

use spenso::{
    antisym, slot,
    structure::{abstract_index::AbstractIndex, representation::LibrarySlot, slot::IsAbstractSlot},
    tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
    function, symbol,
};

use crate::{IndexTooling, test_support::test_initialize};

type Slot = LibrarySlot<AbstractIndex>;

/// Test-local copy of the canonicalizer's stable semantic ordering.
///
/// This oracle is graph-first: Graphica owns unsigned canonical numbering. The
/// key orders signed decorations over a fixed small unsigned family and detects
/// when explicit transforms return one rendered representative with both
/// parities; it is deliberately not a second graph-labeling order.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct ReferenceSymbolKey {
    name: String,
    wildcard_level: u8,
    attributes: [bool; 8],
    tags: Vec<String>,
}

impl From<Symbol> for ReferenceSymbolKey {
    fn from(symbol: Symbol) -> Self {
        let mut tags = symbol.get_tags().to_vec();
        tags.sort();
        Self {
            name: symbol.get_name().to_owned(),
            wildcard_level: symbol.get_wildcard_level(),
            attributes: [
                symbol.is_symmetric(),
                symbol.is_antisymmetric(),
                symbol.is_cyclesymmetric(),
                symbol.is_linear(),
                symbol.is_scalar(),
                symbol.is_real(),
                symbol.is_integer(),
                symbol.is_positive(),
            ],
            tags,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum ReferenceAtomKey {
    Number(String),
    Variable(ReferenceSymbolKey),
    Function {
        head: ReferenceSymbolKey,
        arguments: Vec<Self>,
    },
    Power {
        base: Box<Self>,
        exponent: Box<Self>,
    },
    Product(Vec<Self>),
    Sum(Vec<Self>),
}

impl ReferenceAtomKey {
    fn new(value: AtomView<'_>) -> Self {
        match value {
            AtomView::Num(number) => Self::Number(number.as_view().to_owned().to_string()),
            AtomView::Var(variable) => Self::Variable(variable.get_symbol().into()),
            AtomView::Fun(function) => Self::Function {
                head: function.get_symbol().into(),
                arguments: function.iter().map(Self::new).collect(),
            },
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Self::Power {
                    base: Box::new(Self::new(base)),
                    exponent: Box::new(Self::new(exponent)),
                }
            }
            AtomView::Mul(product) => {
                let mut factors = product.iter().map(Self::new).collect::<Vec<_>>();
                factors.sort();
                Self::Product(factors)
            }
            AtomView::Add(sum) => {
                let mut terms = sum.iter().map(Self::new).collect::<Vec<_>>();
                terms.sort();
                Self::Sum(terms)
            }
        }
    }
}

#[derive(Clone)]
struct SignedCandidate {
    unsigned: Atom,
    negative: bool,
    actual: Atom,
}

fn signed(value: Atom, negative: bool) -> Atom {
    if negative { -value } else { value }
}

fn dummies() -> [Slot; 4] {
    let rep = test_initialize().mink4;
    [
        slot!(rep, AbstractIndex::Dummy(0)).cast(),
        slot!(rep, AbstractIndex::Dummy(1)).cast(),
        slot!(rep, AbstractIndex::Dummy(2)).cast(),
        slot!(rep, AbstractIndex::Dummy(3)).cast(),
    ]
}

fn intrinsic_factor(
    tensor: Symbol,
    left: Symbol,
    right: Symbol,
    [first, second]: [Slot; 2],
    rename_flip: bool,
    site_flip: bool,
) -> Atom {
    let [first, second] = if rename_flip {
        [second, first]
    } else {
        [first, second]
    };
    let [tensor_first, tensor_second] = if site_flip {
        [second, first]
    } else {
        [first, second]
    };
    function!(tensor, tensor_first.to_atom(), tensor_second.to_atom())
        * function!(left, second.to_atom())
        * function!(right, first.to_atom())
}

fn reaches_both_orientations(candidates: &[SignedCandidate]) -> bool {
    let mut orientations = BTreeMap::<ReferenceAtomKey, u8>::new();
    for candidate in candidates {
        let bit = if candidate.negative { 2 } else { 1 };
        *orientations
            .entry(ReferenceAtomKey::new(candidate.unsigned.as_view()))
            .or_default() |= bit;
    }
    orientations.values().any(|orientations| *orientations == 3)
}

fn scalar_orbit_minimum(candidates: &[SignedCandidate]) -> Atom {
    if reaches_both_orientations(candidates) {
        return Atom::Zero;
    }

    candidates
        .iter()
        .min_by_key(|candidate| ReferenceAtomKey::new(candidate.actual.as_view()))
        .expect("a bounded oracle orbit is nonempty")
        .actual
        .clone()
}

fn semantic_minimum(candidates: &[Atom]) -> Atom {
    candidates
        .iter()
        .min_by_key(|candidate| ReferenceAtomKey::new(candidate.as_view()))
        .expect("a bounded oracle orbit is nonempty")
        .clone()
}

fn assert_production_orbit(case: &str, candidates: &[Atom], expected: &Atom) {
    for (candidate, input) in candidates.iter().enumerate() {
        assert_eq!(
            input
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            *expected,
            "{case} candidate {candidate}"
        );
    }
}

fn assert_production_closure(case: &str, candidates: &[Atom]) -> Atom {
    let expected = candidates[0]
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    assert_production_orbit(case, candidates, &expected);
    expected
}

#[test]
fn intrinsic_symmetry_factor_exchange_and_dummy_relabeling_match_the_reference_minimum() {
    let slots = dummies();
    let tensor = tensor_symbol!(brute_force_intrinsic_tensor; Antisymmetric);
    let left = tensor_symbol!(brute_force_intrinsic_left);
    let right = tensor_symbol!(brute_force_intrinsic_right);
    let mut candidates = Vec::with_capacity(64);

    for block_swap in [false, true] {
        for first_rename in [false, true] {
            for second_rename in [false, true] {
                for first_site in [false, true] {
                    for second_site in [false, true] {
                        for factor_swap in [false, true] {
                            let blocks = if block_swap {
                                [[slots[2], slots[3]], [slots[0], slots[1]]]
                            } else {
                                [[slots[0], slots[1]], [slots[2], slots[3]]]
                            };
                            let factors = [
                                intrinsic_factor(
                                    tensor,
                                    left,
                                    right,
                                    blocks[0],
                                    first_rename,
                                    first_site,
                                ),
                                intrinsic_factor(
                                    tensor,
                                    left,
                                    right,
                                    blocks[1],
                                    second_rename,
                                    second_site,
                                ),
                            ];
                            let unsigned = if factor_swap {
                                &factors[1] * &factors[0]
                            } else {
                                &factors[0] * &factors[1]
                            };
                            let negative = first_site ^ second_site;
                            candidates.push(SignedCandidate {
                                actual: signed(unsigned.clone(), negative),
                                unsigned,
                                negative,
                            });
                        }
                    }
                }
            }
        }
    }

    assert_eq!(candidates.len(), 64);
    assert!(!reaches_both_orientations(&candidates));
    let inputs = candidates
        .into_iter()
        .map(|candidate| candidate.actual)
        .collect::<Vec<_>>();
    assert!(!assert_production_closure("intrinsic", &inputs).is_zero());
}

#[test]
fn partial_antisymmetry_matches_the_nested_sign_reference_minimum() {
    let [first, second, ..] = dummies();
    let tensor = tensor_symbol!(brute_force_partial_tensor);
    let left = tensor_symbol!(brute_force_partial_left);
    let right = tensor_symbol!(brute_force_partial_right);
    let parameter = Atom::var(symbol!("brute_force_partial_parameter"));
    let mut candidates = Vec::with_capacity(4);

    for rename_flip in [false, true] {
        for site_flip in [false, true] {
            let [first, second] = if rename_flip {
                [second, first]
            } else {
                [first, second]
            };
            let group = if site_flip {
                -antisym!(second, first)
            } else {
                antisym!(first, second)
            };
            candidates.push(
                FunctionBuilder::new(tensor)
                    .add_arg(parameter.clone())
                    .add_arg(group)
                    .finish()
                    * function!(left, second.to_atom())
                    * function!(right, first.to_atom()),
            );
        }
    }

    assert_eq!(candidates.len(), 4);
    let expected = semantic_minimum(&candidates);
    assert_production_orbit("partial", &candidates, &expected);
}

#[test]
fn declared_unary_function_semantics_match_the_reference_minimum() {
    let [first, second, ..] = dummies();
    let tensor = tensor_symbol!(brute_force_function_tensor; Antisymmetric);
    let left = tensor_symbol!(brute_force_function_left);
    let right = tensor_symbol!(brute_force_function_right);
    let nonlinear = symbol!("brute_force_function_nonlinear");
    let linear = symbol!("brute_force_function_linear"; Linear);

    let mut nonlinear_candidates = Vec::with_capacity(4);
    let mut linear_candidates = Vec::with_capacity(4);
    for rename_flip in [false, true] {
        for site_flip in [false, true] {
            let body =
                intrinsic_factor(tensor, left, right, [first, second], rename_flip, site_flip);
            nonlinear_candidates.push(function!(nonlinear, signed(body.clone(), site_flip)));
            let unsigned = function!(linear, body);
            linear_candidates.push(SignedCandidate {
                actual: signed(unsigned.clone(), site_flip),
                unsigned,
                negative: site_flip,
            });
        }
    }

    assert_eq!(nonlinear_candidates.len(), 4);
    let expected = semantic_minimum(&nonlinear_candidates);
    assert_production_orbit("nonlinear function", &nonlinear_candidates, &expected);

    assert_eq!(linear_candidates.len(), 4);
    assert!(!reaches_both_orientations(&linear_candidates));
    let inputs = linear_candidates
        .into_iter()
        .map(|candidate| candidate.actual)
        .collect::<Vec<_>>();
    assert!(!assert_production_closure("linear function", &inputs).is_zero());
}

#[test]
fn signed_sum_branch_action_detects_both_orientations_of_one_representative() {
    let [first, second, ..] = dummies();
    let tensor = tensor_symbol!(brute_force_sum_tensor; Antisymmetric);
    let left = tensor_symbol!(brute_force_sum_left);
    let right = tensor_symbol!(brute_force_sum_right);
    let mut candidates = Vec::with_capacity(8);

    for rename_flip in [false, true] {
        for site_flip in [false, true] {
            for branch_swap in [false, true] {
                let [first, second] = if rename_flip {
                    [second, first]
                } else {
                    [first, second]
                };
                let [tensor_first, tensor_second] = if site_flip {
                    [second, first]
                } else {
                    [first, second]
                };
                let branches = [
                    function!(left, first.to_atom()) * function!(right, second.to_atom()),
                    function!(left, second.to_atom()) * function!(right, first.to_atom()),
                ];
                let sum = if branch_swap {
                    &branches[1] + &branches[0]
                } else {
                    &branches[0] + &branches[1]
                };
                let unsigned =
                    function!(tensor, tensor_first.to_atom(), tensor_second.to_atom()) * sum;
                candidates.push(SignedCandidate {
                    actual: signed(unsigned.clone(), site_flip),
                    unsigned,
                    negative: site_flip,
                });
            }
        }
    }

    assert_eq!(candidates.len(), 8);
    let expected = scalar_orbit_minimum(&candidates);
    assert!(expected.is_zero());
    let inputs = candidates
        .into_iter()
        .map(|candidate| candidate.actual)
        .collect::<Vec<_>>();
    assert_production_orbit("signed Sum", &inputs, &expected);
}
