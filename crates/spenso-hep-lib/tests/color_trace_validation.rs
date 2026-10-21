//! Exact finite-component checks, independent of the colour identity kernels.
//!
//! The input oracle multiplies the HEP library's matrices in the supplied order.
//! The output oracle evaluates only components, including normalized symmetric
//! matrix averages. Neither route expands a symbolic tensor or a projector.

use std::collections::{BTreeMap, HashMap};

use idenso::{
    color::{CS, ColorSimplifySettings},
    representations::{ColorAdjoint, ColorFundamental, initialize},
    tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus, SymbolicTensor},
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing::{self, SYM},
    structure::{abstract_index::AbstractIndex, representation::RepName, slot::ParseableAind},
    tensors::data::GetTensorData,
};
use spenso_hep_lib::{su3_generator_data_atom, su3_structure_f_data_atom};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, parse_lit, symbol,
};

type Matrix = Vec<Vec<Atom>>;
type Bindings = BTreeMap<Atom, usize>;

// Fully distributed term count, without expanding: sums add, products multiply.
fn monomials(view: AtomView<'_>) -> u128 {
    match view {
        AtomView::Add(sum) => sum.iter().map(monomials).sum(),
        AtomView::Mul(product) => product.iter().map(monomials).product(),
        _ => 1,
    }
}

// A certified colour fixed point must survive a fresh admission unchanged;
// the policy cache answers the ordinary rerun without running a kernel.
fn assert_fresh_fixed_point(case: &str, reduced: &Atom) {
    let fresh = SymbolicTensor::infer(reduced.clone()).unwrap();
    let again = fresh.simplify_algebra(&settings()).unwrap();
    assert_eq!(
        again.expression(),
        reduced,
        "{case}: re-admitted colour output is not a fixed point"
    );
}

/// Reduce arithmetic of an already selected scalar matrix component in
/// Q(i)[sqrt(3)]. This never receives an unevaluated tensor expression.
fn scalar(value: Atom) -> Atom {
    value
        .expand()
        .replace(parse_lit!(sqrt(3)).pow(2).to_pattern())
        .with(Atom::num(3))
}

fn identity(size: usize) -> Matrix {
    (0..size)
        .map(|i| (0..size).map(|j| Atom::num(i64::from(i == j))).collect())
        .collect()
}

fn multiply(left: &Matrix, right: &Matrix) -> Matrix {
    let size = left.len();
    (0..size)
        .map(|i| {
            (0..size)
                .map(|j| {
                    scalar(Atom::add_many(
                        (0..size)
                            .filter(|&k| !left[i][k].is_zero() && !right[k][j].is_zero())
                            .map(|k| &left[i][k] * &right[k][j]),
                    ))
                })
                .collect()
        })
        .collect()
}

struct Components {
    fundamental: Vec<Matrix>,
    adjoint: Vec<Matrix>,
    traces: HashMap<(bool, bool, Vec<usize>), Atom>,
}

impl Components {
    fn new() -> Self {
        initialize();
        let t = su3_generator_data_atom(CS.t_strct::<AbstractIndex>(3, 8));
        let f = su3_structure_f_data_atom(CS.f_strct::<AbstractIndex>(8));
        let fundamental = (0..8)
            .map(|a| {
                (0..3)
                    .map(|i| {
                        (0..3)
                            .map(|j| {
                                t.canonical()
                                    .get_owned(t.layout().logical_to_canonical(&[a, i, j]))
                                    .unwrap_or(Atom::Zero)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        // These are the actual raw f(in,out,a) matrices, not Hermitian
        // adjoint generators. Their odd/even trace phases must not be guessed.
        let adjoint = (0..8)
            .map(|a| {
                (0..8)
                    .map(|i| {
                        (0..8)
                            .map(|j| {
                                f.canonical()
                                    .get_owned(f.layout().logical_to_canonical(&[i, j, a]))
                                    .unwrap_or(Atom::Zero)
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        Self {
            fundamental,
            adjoint,
            traces: HashMap::new(),
        }
    }

    fn trace(&mut self, adjoint: bool, symmetric: bool, mut word: Vec<usize>) -> Atom {
        if symmetric {
            word.sort_unstable();
        }
        let key = (adjoint, symmetric, word);
        if let Some(value) = self.traces.get(&key) {
            return value.clone();
        }
        let matrices = if adjoint {
            &self.adjoint
        } else {
            &self.fundamental
        };
        let size = matrices[0].len();
        let matrix = if symmetric {
            // S(mask) = sum_j S(mask\j) M_j / |mask|. Occurrence positions
            // remain distinct even when their component values coincide.
            let mut averages = vec![identity(size)];
            for mask in 1usize..(1 << key.2.len()) {
                let mut average = vec![vec![Atom::Zero; size]; size];
                for (position, &a) in key.2.iter().enumerate() {
                    if mask & (1 << position) == 0 {
                        continue;
                    }
                    let term = multiply(&averages[mask ^ (1 << position)], &matrices[a]);
                    for (row, source) in average.iter_mut().zip(term) {
                        for (entry, value) in row.iter_mut().zip(source) {
                            *entry += value;
                        }
                    }
                }
                for row in &mut average {
                    for entry in row {
                        *entry = scalar(&*entry / Atom::num(mask.count_ones()));
                    }
                }
                averages.push(average);
            }
            averages.pop().unwrap()
        } else {
            key.2.iter().fold(identity(size), |product, &a| {
                multiply(&product, &matrices[a])
            })
        };
        let value = scalar(Atom::add_many((0..size).map(|i| matrix[i][i].clone())));
        self.traces.insert(key, value.clone());
        value
    }

    /// gram(n,R,S): the contraction of the two normalized symmetric
    /// invariants. Raw adjoint words carry i^n relative to that normalization.
    fn gram(&mut self, degree: usize, adjoint: [bool; 2]) -> Atom {
        let mut total = Atom::Zero;
        for tuple in 0..8usize.pow(degree as u32) {
            let word = (0..degree)
                .map(|position| tuple / 8usize.pow(position as u32) % 8)
                .collect::<Vec<_>>();
            let mut product = Atom::num(1);
            for side in adjoint {
                let value = self.trace(side, true, word.clone());
                product *= if side {
                    value / Atom::i().pow(Atom::num(degree as i64))
                } else {
                    value
                };
            }
            total = scalar(total + product);
        }
        total
    }

    fn compile(&mut self, value: AtomView<'_>) -> Component {
        let kind = match value {
            AtomView::Num(_) => Kind::Scalar(value.to_owned()),
            AtomView::Var(variable) => {
                let name = variable.get_symbol();
                Kind::Scalar(if name == symbol!("color_trace_validation::Nc") {
                    Atom::num(3)
                } else if name == symbol!("color_trace_validation::Na") {
                    Atom::num(8)
                } else if let Some(value) = spectator_scalar(name.get_name()) {
                    value
                } else {
                    panic!("unassigned scalar {value}")
                })
            }
            AtomView::Add(sum) => Kind::Sum(sum.iter().map(|x| self.compile(x)).collect()),
            AtomView::Mul(product) => {
                Kind::Product(product.iter().map(|x| self.compile(x)).collect())
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                Kind::Power(
                    Box::new(self.compile(base)),
                    i64::try_from(exponent).unwrap(),
                )
            }
            AtomView::Fun(function) => {
                let head = function.get_symbol();
                let args = function.iter().collect::<Vec<_>>();
                if head == CS.f {
                    Kind::Structure(
                        args.into_iter()
                            .map(index)
                            .collect::<Vec<_>>()
                            .try_into()
                            .unwrap(),
                    )
                } else if head == ETS.metric {
                    Kind::Metric(
                        args.into_iter()
                            .map(index)
                            .collect::<Vec<_>>()
                            .try_into()
                            .unwrap(),
                    )
                } else if head == SPENSO_TAG.trace {
                    let (rep, mut factors) = shadowing::trace_parts(function).unwrap();
                    let symmetric = if let [AtomView::Fun(projector)] = factors.as_slice()
                        && projector.get_symbol() == *SYM
                    {
                        factors = projector.iter().collect();
                        true
                    } else {
                        false
                    };
                    let adjoint = rep.get_symbol().unwrap() == CS.adjoint_rep;
                    let labels = factors
                        .into_iter()
                        .map(|factor| {
                            let AtomView::Fun(generator) = factor else {
                                panic!("invalid word {factor}")
                            };
                            assert_eq!(generator.get_symbol(), if adjoint { CS.f } else { CS.t });
                            let args = generator.iter().collect::<Vec<_>>();
                            let position = args.iter().position(|arg| {
                                matches!(arg, AtomView::Fun(slot) if slot.get_symbol() == CS.adjoint_rep)
                            }).expect("a generator has an indexed adjoint argument");
                            if adjoint {
                                // A cyclic placement of the external argument
                                // preserves f(in,out,a), hence its matrix sign.
                                assert_eq!(args[(position + 1) % 3], Atom::var(SPENSO_TAG.chain_in).as_view());
                                assert_eq!(args[(position + 2) % 3], Atom::var(SPENSO_TAG.chain_out).as_view());
                            }
                            index(args[position])
                        })
                        .collect();
                    Kind::Trace {
                        adjoint,
                        symmetric,
                        labels,
                    }
                } else if head == CS.gram {
                    let degree = usize::try_from(args[0]).unwrap();
                    let adjoint =
                        [args[1], args[2]].map(|rep| rep.get_symbol().unwrap() == CS.adjoint_rep);
                    Kind::Scalar(self.gram(degree, adjoint))
                } else if head == CS.idx || head == CS.cas {
                    assert_eq!(usize::try_from(args[0]).unwrap(), 2);
                    let adjoint = args[1].get_symbol().unwrap() == CS.adjoint_rep;
                    let sign = if adjoint { -1 } else { 1 };
                    // Derive quadratic invariants from the exact matrices too.
                    let trace = if head == CS.idx {
                        self.trace(adjoint, false, vec![0, 0])
                    } else {
                        Atom::add_many((0..8).map(|a| self.trace(adjoint, false, vec![a, a])))
                            / Atom::num(if adjoint { 8 } else { 3 })
                    };
                    Kind::Scalar(scalar(Atom::num(sign) * trace))
                } else if head.get_name().starts_with("color_trace_validation::")
                    && args.iter().all(|arg| matches!(arg, AtomView::Fun(slot) if slot.get_symbol() == CS.adjoint_rep))
                {
                    Kind::Foreign(head.get_name().to_string(), args.into_iter().map(index).collect())
                } else {
                    panic!("unsupported component node {value}")
                }
            }
        };
        Component::new(kind)
    }
}

/// Values of the named scalar coefficients and spectators of the
/// keep-the-rest cases.
fn spectator_scalar(name: &str) -> Option<Atom> {
    Some(match name.strip_prefix("color_trace_validation::keep_")? {
        "v0" => Atom::num(2),
        "v1" => Atom::num(-3),
        "v2" => Atom::num(7) / Atom::num(2),
        "x" => Atom::num(5),
        _ => return None,
    })
}

fn index(slot: AtomView<'_>) -> Atom {
    let AtomView::Fun(slot) = slot else {
        panic!("expected adjoint slot")
    };
    assert_eq!(slot.get_symbol(), CS.adjoint_rep);
    assert_eq!(slot.get_nargs(), 2);
    slot.iter().nth(1).unwrap().to_owned()
}

enum Kind {
    Scalar(Atom),
    Sum(Vec<Component>),
    Product(Vec<Component>),
    Power(Box<Component>, i64),
    Structure([Atom; 3]),
    Metric([Atom; 2]),
    Trace {
        adjoint: bool,
        symmetric: bool,
        labels: Vec<Atom>,
    },
    /// A foreign tensor with adjoint slots: fixed small integer components.
    Foreign(String, Vec<Atom>),
}

struct Component {
    kind: Kind,
    open: Vec<Atom>,
    summed: Vec<Atom>,
}

impl Component {
    fn new(kind: Kind) -> Self {
        let labels = match &kind {
            Kind::Scalar(_) => vec![],
            Kind::Sum(terms) => {
                let open = terms[0].open.clone();
                assert!(terms.iter().all(|term| term.open == open));
                return Self {
                    kind,
                    open,
                    summed: vec![],
                };
            }
            Kind::Product(factors) => factors.iter().flat_map(|x| x.open.clone()).collect(),
            Kind::Power(base, exponent) => {
                if base.open.is_empty() {
                    vec![]
                } else {
                    assert_eq!(*exponent, 2);
                    base.open.iter().chain(&base.open).cloned().collect()
                }
            }
            Kind::Structure(labels) => labels.to_vec(),
            Kind::Metric(labels) => labels.to_vec(),
            Kind::Trace { labels, .. } => labels.clone(),
            Kind::Foreign(_, labels) => labels.clone(),
        };
        let mut counts = BTreeMap::<_, usize>::new();
        for label in labels {
            *counts.entry(label).or_default() += 1;
        }
        assert!(counts.values().all(|&count| count <= 2));
        let (mut open, mut summed) = (vec![], vec![]);
        for (label, count) in counts {
            if count == 2 {
                summed.push(label);
            } else {
                open.push(label);
            }
        }
        Self { kind, open, summed }
    }

    fn evaluate(&self, data: &mut Components, bindings: &mut Bindings) -> Atom {
        self.sum_indices(data, bindings, 0)
    }

    fn bound_zero(&self, data: &Components, bindings: &Bindings) -> bool {
        match &self.kind {
            Kind::Scalar(value) => value.is_zero(),
            Kind::Structure([a, b, c]) => {
                if let (Some(&a), Some(&b), Some(&c)) =
                    (bindings.get(a), bindings.get(b), bindings.get(c))
                {
                    data.adjoint[c][a][b].is_zero()
                } else {
                    false
                }
            }
            Kind::Metric([a, b]) => {
                matches!((bindings.get(a), bindings.get(b)), (Some(a), Some(b)) if a != b)
            }
            Kind::Product(factors) => factors
                .iter()
                .any(|factor| factor.bound_zero(data, bindings)),
            Kind::Power(base, exponent) => *exponent > 0 && base.bound_zero(data, bindings),
            _ => false,
        }
    }

    fn inexpensive(&self) -> bool {
        match &self.kind {
            Kind::Scalar(_) | Kind::Structure(_) | Kind::Metric(_) | Kind::Foreign(..) => true,
            Kind::Power(base, _) => base.inexpensive(),
            _ => false,
        }
    }

    fn sum_indices(&self, data: &mut Components, bindings: &mut Bindings, position: usize) -> Atom {
        // Sparse HEP components can settle a product before the remaining
        // coordinates are bound. This is exact component lookup, not a colour
        // identity or an approximation to the remaining contraction.
        if self.bound_zero(data, bindings) {
            return Atom::Zero;
        }
        if let Some(label) = self.summed.get(position) {
            assert!(!bindings.contains_key(label));
            let result = Atom::add_many((0..8).map(|value| {
                bindings.insert(label.clone(), value);
                self.sum_indices(data, bindings, position + 1)
            }));
            bindings.remove(label);
            return scalar(result);
        }
        match &self.kind {
            Kind::Scalar(value) => value.clone(),
            Kind::Sum(terms) => scalar(Atom::add_many(
                terms.iter().map(|x| x.evaluate(data, bindings)),
            )),
            Kind::Product(factors) => {
                let mut product = Atom::num(1);
                for factor in factors
                    .iter()
                    .filter(|factor| factor.inexpensive())
                    .chain(factors.iter().filter(|factor| !factor.inexpensive()))
                {
                    let value = factor.evaluate(data, bindings);
                    if value.is_zero() {
                        return Atom::Zero;
                    }
                    product = scalar(product * value);
                }
                product
            }
            Kind::Power(base, exponent) => scalar(base.evaluate(data, bindings).pow(*exponent)),
            Kind::Structure([a, b, c]) => {
                data.adjoint[bindings[c]][bindings[a]][bindings[b]].clone()
            }
            Kind::Metric([a, b]) => Atom::num(i64::from(bindings[a] == bindings[b])),
            Kind::Trace {
                adjoint,
                symmetric,
                labels,
            } => data.trace(
                *adjoint,
                *symmetric,
                labels.iter().map(|label| bindings[label]).collect(),
            ),
            Kind::Foreign(name, labels) => {
                let seed = name
                    .bytes()
                    .map(usize::from)
                    .chain(labels.iter().map(|label| bindings[label]))
                    .fold(17usize, |hash, value| {
                        hash.wrapping_mul(31).wrapping_add(value)
                    });
                Atom::num((seed % 7) as i64 - 3)
            }
        }
    }
}

fn settings() -> AlgebraSettings {
    AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    }
}

fn word(adjoint: bool, labels: &[Atom]) -> Atom {
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let representation = if adjoint {
        coad.to_symbolic([])
    } else {
        ColorFundamental {}
            .new_rep(symbol!("color_trace_validation::Nc"))
            .to_symbolic([])
    };
    spenso::trace!(representation; labels.iter().map(|label| {
        let slot = coad.to_symbolic([label.clone()]);
        if adjoint {
            idenso::color_f!(Atom::var(SPENSO_TAG.chain_in), Atom::var(SPENSO_TAG.chain_out), slot)
        } else {
            idenso::color_t!(slot)
        }
    }))
}

fn assert_trace_case(adjoint: bool, length: usize, reverse: bool) {
    let mut components = Components::new();
    let label = spenso::index_symbol!("color_trace_validation::edge");
    let labels = (0..length)
        .map(|i| function!(label, i, 1))
        .collect::<Vec<_>>();
    let original = word(adjoint, &labels);
    let mut rotated = labels.clone();
    rotated.rotate_left(2);
    assert_eq!(word(adjoint, &rotated), original);
    let case = format!("trace/{adjoint}/{length}/{reverse}");
    let order = if reverse {
        labels.iter().rev().cloned().collect()
    } else {
        labels.clone()
    };
    let source = SymbolicTensor::infer(word(adjoint, &order)).unwrap();
    let reduced = source.simplify_algebra(&settings()).unwrap();
    assert_fresh_fixed_point(&case, reduced.expression());
    assert_ne!(source.expression(), reduced.expression());
    assert_eq!(source.structure(), reduced.structure());
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.simplify_algebra(&settings()).unwrap(), reduced);
    let explicit = reduced.undo_chain().unwrap();
    let actual = components.compile(explicit.expression().as_view());
    let mut nonzero = false;
    for sample in 0..3 {
        let values = match sample {
            0 => vec![7; length],
            1 => (0..length)
                .map(|i| if i + 2 < length { 0 } else { i + 3 - length })
                .collect(),
            _ => (0..length).map(|i| (i * 3 + length) % 8).collect(),
        };
        let mut bindings: Bindings = labels.iter().cloned().zip(values).collect();
        let ordered_values = order.iter().map(|label| bindings[label]).collect();
        // This is a direct matrix product, not the symbolic trace parser.
        let expected = components.trace(adjoint, false, ordered_values);
        nonzero |= !expected.is_zero();
        assert_eq!(
            actual.evaluate(&mut components, &mut bindings),
            expected,
            "adjoint={adjoint}, length={length}, reverse={reverse}, sample={sample}"
        );
    }
    assert!(
        nonzero,
        "all component checks were trivial at length {length}"
    );
}

fn assert_spectator_case(adjoint: bool) {
    initialize();
    let label = spenso::index_symbol!("color_trace_validation::spectator_edge");
    let labels = (0..8).map(|i| function!(label, i, 2)).collect::<Vec<_>>();
    let spectator = parse_lit!((color_trace_spectator_x + color_trace_spectator_y) ^ 17);
    let source = SymbolicTensor::infer(word(adjoint, &labels) * &spectator).unwrap();
    for color in [
        None,
        Some(ColorSimplifySettings::default().without_trace_evaluation()),
        Some(
            ColorSimplifySettings::default()
                .without_trace_evaluation()
                .without_cross_chain_fierz_expansion(),
        ),
    ] {
        let disabled = AlgebraSettings {
            color,
            ..settings()
        };
        assert_eq!(source.simplify_algebra(&disabled).unwrap(), source);
    }
    let reduced = source.simplify_algebra(&settings()).unwrap();
    assert_eq!(source.structure(), reduced.structure());
    assert!(
        reduced
            .expression()
            .replace(spectator.to_pattern())
            .with(Atom::Zero)
            .is_zero()
    );
    let mut ordinary_trace = false;
    reduced.expression().visitor(&mut |node| {
            if let AtomView::Fun(function) = node
                && let Some((_, factors)) = shadowing::trace_parts(function)
            {
                ordinary_trace |= !matches!(factors.as_slice(), [AtomView::Fun(projector)] if projector.get_symbol() == *SYM);
            }
            if let AtomView::Fun(function) = node {
                let expected_dimension = if function.get_symbol() == CS.adjoint_rep {
                    Some(symbol!("color_trace_validation::Na"))
                } else if function.get_symbol() == CS.fundamental_rep {
                    Some(symbol!("color_trace_validation::Nc"))
                } else { None };
                if let Some(dimension) = expected_dimension {
                    assert_eq!(function.iter().next().unwrap(), Atom::var(dimension).as_view());
                }
            }
            true
        });
    assert!(
        !ordinary_trace,
        "an eligible ordered colour trace was left behind"
    );
}

fn assert_two_plus_two_prefix(adjoint: bool) {
    let mut components = Components::new();
    let label = spenso::index_symbol!("color_trace_validation::prefix_edge");
    let labels = [0, 1, 2, 3].map(|axis| function!(label, Atom::num(axis), Atom::num(2)));
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let representation = if adjoint {
        coad.to_symbolic([])
    } else {
        ColorFundamental {}
            .new_rep(symbol!("color_trace_validation::Nc"))
            .to_symbolic([])
    };
    let factors = labels.each_ref().map(|label| {
        let slot = coad.to_symbolic([label.clone()]);
        if adjoint {
            idenso::color_f!(
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out),
                slot
            )
        } else {
            idenso::color_t!(slot)
        }
    });
    let expression = spenso::trace!(
        representation,
        shadowing::sym(factors[..2].iter().cloned()),
        &factors[2],
        &factors[3]
    );
    let source = SymbolicTensor::infer(expression).unwrap();
    let result = source.simplify_algebra(&settings()).unwrap();
    assert_eq!(source.structure(), result.structure());
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings()).unwrap(), result);
    assert!(!result.expression().contains_symbol(ETS.metric));
    assert!(
        result
            .expression()
            .contains_symbol(if adjoint { CS.cas } else { CS.idx })
    );
    assert!(result.expression().contains_symbol(SPENSO_TAG.trace));
    let actual = components.compile(result.expression().as_view());
    let points = [[7, 7, 7, 7], [0, 0, 1, 1], [0, 1, 2, 7], [3, 4, 5, 6]]
        .into_iter()
        .chain((0..32usize).map(|sample| {
            [
                sample % 8,
                (sample / 4) % 8,
                (sample * 3 + 1) % 8,
                (sample * 5 + 2) % 8,
            ]
        }));
    for values in points {
        let first = components.trace(adjoint, false, values.to_vec());
        let second = components.trace(
            adjoint,
            false,
            vec![values[1], values[0], values[2], values[3]],
        );
        let expected = scalar((first + second) / Atom::num(2));
        let mut bindings = labels.iter().cloned().zip(values).collect();
        assert_eq!(
            actual.evaluate(&mut components, &mut bindings),
            expected,
            "two-plus-two adjoint={adjoint}, components={values:?}"
        );
    }
}

fn assert_repeated_case(adjoint: bool) {
    let mut components = Components::new();
    let label = spenso::index_symbol!("color_trace_validation::repeated_edge");
    let labels = (0..7).map(|i| function!(label, i, 3)).collect::<Vec<_>>();
    let sequence = [0, 1, 2, 3, 0, 4, 5, 6];
    let word_labels = sequence.map(|i| labels[i].clone());
    let case = format!("repeated/{adjoint}/8");
    let source = SymbolicTensor::infer(word(adjoint, &word_labels)).unwrap();
    let reduced = source.simplify_algebra(&settings()).unwrap();
    assert_fresh_fixed_point(&case, reduced.expression());
    assert_eq!(source.structure(), reduced.structure());
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.simplify_algebra(&settings()).unwrap(), reduced);
    let explicit = reduced.undo_chain().unwrap();
    let actual = components.compile(explicit.expression().as_view());
    let mut bindings = labels[1..]
        .iter()
        .cloned()
        .map(|label| (label, 7))
        .collect();
    // Sum the shared colour coordinate explicitly in the original matrix
    // product, independent of both the tensor parser and the output oracle.
    let expected =
        scalar(Atom::add_many((0..8).map(|a| {
            components.trace(adjoint, false, vec![a, 7, 7, 7, a, 7, 7, 7])
        })));
    assert!(!expected.is_zero());
    assert_eq!(actual.evaluate(&mut components, &mut bindings), expected);
}

fn assert_cycle_case(length: usize, flipped: bool) {
    let mut components = Components::new();
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let label = spenso::index_symbol!("color_trace_validation::cycle_external");
    let spectator = parse_lit!((cycle_spectator_x + cycle_spectator_y) ^ 19);
    let external = (0..length)
        .map(|i| function!(label, i, 4))
        .collect::<Vec<_>>();
    let case = format!("cycle/{length}/{flipped}");
    let scope = symbol!(format!(
        "color_trace_validation::cycle_scope_{length}_{flipped}"
    ));
    let internal = (0..length)
        .map(|i| coad.to_symbolic([Atom::from(AbstractIndex::Normal(i).scoped(scope))]))
        .collect::<Vec<_>>();
    let cycle = Atom::mul_many((0..length).map(|i| {
        let left = &internal[i];
        let right = &internal[(i + 1) % length];
        let external = coad.to_symbolic([external[i].clone()]);
        if flipped && i == 0 {
            idenso::color_f!(right, left, external)
        } else {
            idenso::color_f!(left, right, external)
        }
    }));
    let source = SymbolicTensor::infer(cycle * &spectator).unwrap();
    let reduced = source.simplify_algebra(&settings()).unwrap();
    assert_fresh_fixed_point(&case, reduced.expression());
    assert_ne!(source.expression(), reduced.expression());
    assert_eq!(source.structure(), reduced.structure());
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.simplify_algebra(&settings()).unwrap(), reduced);
    assert!(
        reduced
            .expression()
            .replace(spectator.to_pattern())
            .with(Atom::Zero)
            .is_zero()
    );
    if length == 5 {
        // color.h's five-loop table: nine terms, against fifteen from the
        // symmetric-prefix recursion.
        let colour = reduced
            .expression()
            .replace(spectator.to_pattern())
            .with(Atom::num(1));
        assert!(monomials(colour.as_view()) <= 9, "{colour}");
    }
    let explicit = reduced.undo_chain().unwrap();
    let colour = explicit
        .expression()
        .replace(spectator.to_pattern())
        .with(Atom::num(1));
    let actual = components.compile(colour.as_view());
    let values = if length.is_multiple_of(2) {
        vec![7; length]
    } else {
        (0..length)
            .map(|i| if i + 2 < length { 0 } else { i + 3 - length })
            .collect()
    };
    let expected =
        components.trace(true, false, values.clone()) * Atom::num(if flipped { -1 } else { 1 });
    assert!(!expected.is_zero(), "trivial cycle at length {length}");
    let mut bindings = external.iter().cloned().zip(values).collect();
    assert_eq!(
        actual.evaluate(&mut components, &mut bindings),
        expected,
        "indexed adjoint cycle length={length}, flipped={flipped}"
    );
}

macro_rules! component_cases {
    ($helper:ident; $($name:ident: ($($argument:expr),*)),* $(,)?) => {
        $(#[test]
        fn $name() { $helper($($argument),*); })*
    };
}

component_cases!(assert_trace_case;
    fundamental_trace_5: (false, 5, false),
    fundamental_trace_5_reversed: (false, 5, true),
    fundamental_trace_6: (false, 6, false),
    fundamental_trace_6_reversed: (false, 6, true),
    fundamental_trace_7: (false, 7, false),
    fundamental_trace_7_reversed: (false, 7, true),
    fundamental_trace_8: (false, 8, false),
    fundamental_trace_8_reversed: (false, 8, true),
    adjoint_trace_5: (true, 5, false),
    adjoint_trace_5_reversed: (true, 5, true),
    adjoint_trace_6: (true, 6, false),
    adjoint_trace_6_reversed: (true, 6, true),
    adjoint_trace_7: (true, 7, false),
    adjoint_trace_7_reversed: (true, 7, true),
    adjoint_trace_8: (true, 8, false),
    adjoint_trace_8_reversed: (true, 8, true),
);

component_cases!(assert_spectator_case;
    fundamental_trace_preserves_capabilities_and_factorization: (false),
    adjoint_trace_preserves_capabilities_and_factorization: (true),
);

component_cases!(assert_two_plus_two_prefix;
    fundamental_two_plus_two_prefix_matches_exact_matrices: (false),
    adjoint_two_plus_two_prefix_matches_exact_matrices: (true),
);

component_cases!(assert_repeated_case;
    fundamental_trace_8_with_separated_repeated_index: (false),
    adjoint_trace_8_with_separated_repeated_index: (true),
);

component_cases!(assert_cycle_case;
    indexed_adjoint_cycle_5: (5, false),
    indexed_adjoint_cycle_5_flipped: (5, true),
    indexed_adjoint_cycle_6: (6, false),
    indexed_adjoint_cycle_6_flipped: (6, true),
    indexed_adjoint_cycle_7: (7, false),
    indexed_adjoint_cycle_7_flipped: (7, true),
    indexed_adjoint_cycle_8: (8, false),
    indexed_adjoint_cycle_8_flipped: (8, true),
);

fn adjoint_label(name: &str) -> Atom {
    ColorAdjoint {}
        .new_rep(symbol!("color_trace_validation::Na"))
        .to_symbolic([Atom::var(symbol!(&format!(
            "color_trace_validation::{name}"
        )))])
}

fn reduced_scalar(case: &str, source: &Atom, settings: &AlgebraSettings) -> Atom {
    let source = SymbolicTensor::infer(source.clone()).unwrap();
    let reduced = source.simplify_algebra(settings).unwrap();
    assert_eq!(
        reduced.reduction_status(),
        ReductionStatus::Complete,
        "{case}"
    );
    assert_eq!(source.structure(), reduced.structure(), "{case}");
    let mut colour = false;
    reduced.expression().visitor(&mut |node| {
        colour |= matches!(node, AtomView::Fun(function)
            if [CS.f, CS.t, SPENSO_TAG.trace, SPENSO_TAG.chain, ETS.metric].contains(&function.get_symbol()));
        !colour
    });
    assert!(!colour, "{case}: not a scalar: {}", reduced.expression());
    reduced.expression().clone()
}

/// K_{3,3} with a structure constant at each vertex has an automorphism of
/// odd sign. It vanishes; the kernel reaches that zero algebraically.
#[test]
fn signed_graph_zero_reduces_algebraically() {
    let mut components = Components::new();
    let edge = |i: usize, j: usize| adjoint_label(&format!("bipartite_{i}{j}"));
    let graph = Atom::mul_many(
        (0..3)
            .map(|i| idenso::color_f!(edge(i, 0), edge(i, 1), edge(i, 2)))
            .chain((0..3).map(|j| idenso::color_f!(edge(0, j), edge(1, j), edge(2, j)))),
    );
    let expected = components
        .compile(graph.as_view())
        .evaluate(&mut components, &mut Bindings::new());
    assert!(expected.is_zero());
    assert_eq!(
        reduced_scalar("signed_zero", &graph, &settings()),
        Atom::Zero
    );
}

/// Colour labels spelled like Symbolica wildcards: the reduction runs on
/// private aliases, so normalization zeroes f(x_,x_,y) and sorts f. The
/// closed four-loop graph FK0522 is FORM's N_A C_A^4/4 with either spelling.
#[test]
fn wildcard_spelled_labels_reach_structural_zeros() {
    let mut components = Components::new();
    let [a, b, c] = ["a_", "b_", "c_"].map(adjoint_label);
    let source = idenso::color_f!(&a, &b, &c) * spenso::g!(&a, &b);
    let oracle = components.compile(source.as_view());
    for value in 0..8 {
        let mut bindings = Bindings::from([(index_label(&c), value)]);
        assert!(oracle.evaluate(&mut components, &mut bindings).is_zero());
    }
    let reduced = SymbolicTensor::infer(source)
        .unwrap()
        .simplify_algebra(&settings())
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert!(reduced.expression().is_zero(), "{}", reduced.expression());

    let edges = [
        [1, 2, 3],
        [1, 4, 5],
        [2, 6, 7],
        [3, 8, 9],
        [10, 11, 4],
        [10, 8, 12],
        [12, 9, 5],
        [11, 6, 7],
    ];
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let expected = Atom::var(symbol!("color_trace_validation::Na"))
        * idenso::color_cas!(2, coad).pow(Atom::num(4))
        / Atom::num(4);
    for suffix in ["", "_"] {
        let graph = Atom::mul_many(edges.map(|labels| {
            let [x, y, z] = labels.map(|label| adjoint_label(&format!("fk0522_a{label}{suffix}")));
            idenso::color_f!(x, y, z)
        }));
        let reduced = reduced_scalar("fk0522", &graph, &settings());
        assert!(
            (&reduced - &expected).expand().is_zero(),
            "suffix {suffix:?}: {reduced}"
        );
    }
}

/// Two copies of one five-cycle with different internal labels, behind a
/// factored spectator: both cuts give the same canonical states, which
/// merge, and the spectator stays factored.
#[test]
fn equal_generated_states_merge_behind_a_spectator() {
    let mut components = Components::new();
    let spectator = parse_lit!((merge_spectator_x + merge_spectator_y) ^ 7);
    let external = ["a_", "b_", "c_", "d_", "e_"].map(adjoint_label);
    let cycle = |prefix: &str| {
        let internal = (0..5)
            .map(|i| adjoint_label(&format!("{prefix}{i}")))
            .collect::<Vec<_>>();
        Atom::mul_many(
            (0..5).map(|i| idenso::color_f!(&internal[i], &internal[(i + 1) % 5], &external[i])),
        )
    };
    let source = (cycle("merge_x") + cycle("merge_y")) * &spectator;
    let reduced = SymbolicTensor::infer(source.clone())
        .unwrap()
        .simplify_algebra(&settings())
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    let colour = reduced
        .expression()
        .replace(spectator.to_pattern())
        .with(Atom::num(1));
    assert!(
        reduced
            .expression()
            .replace(spectator.to_pattern())
            .with(Atom::Zero)
            .is_zero(),
        "the spectator is not kept factored: {}",
        reduced.expression()
    );
    assert!(monomials(colour.as_view()) <= 9, "{colour}");
    let expected = components.compile(
        source
            .replace(spectator.to_pattern())
            .with(Atom::num(1))
            .as_view(),
    );
    let actual = components.compile(colour.as_view());
    let mut nonzero = false;
    for values in [
        [0, 0, 0, 1, 2],
        [7, 7, 7, 7, 7],
        [0, 1, 2, 3, 4],
        [3, 4, 5, 6, 7],
    ] {
        let mut bindings = external.iter().map(index_label).zip(values).collect();
        let value = expected.evaluate(&mut components, &mut bindings);
        nonzero |= !value.is_zero();
        assert_eq!(
            actual.evaluate(&mut components, &mut bindings),
            value,
            "{values:?}"
        );
    }
    assert!(nonzero);
}

fn index_label(slot: &Atom) -> Atom {
    index(slot.as_view())
}

/// An f on two legs of one line at any distance is absorbed (color.h):
/// X^a W X^c R f^{ace} = κ (C_A/2) X^e W R + κ Σ_k f^{w_k c y_k} f^{ace} X^a W[w_k→y_k] R.
/// Traces in both representations and an open chain, plain and wildcard
/// labels, behind a factored spectator. An open chain M is probed by the
/// independent closures Tr(M) and Tr(M T^z).
#[test]
fn structure_constant_on_distant_legs_is_absorbed() {
    let mut components = Components::new();
    let spectator = parse_lit!((absorb_spectator_x + absorb_spectator_y) ^ 3);
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let fundamental = ColorFundamental {}
        .new_rep(symbol!("color_trace_validation::Nc"))
        .to_symbolic([]);
    let line_in = ColorFundamental {}
        .new_rep(symbol!("color_trace_validation::Nc"))
        .to_symbolic([Atom::var(symbol!("color_trace_validation::absorb_i"))]);
    let line_out = spenso::dind!(
        ColorFundamental {}
            .new_rep(symbol!("color_trace_validation::Nc"))
            .to_symbolic([Atom::var(symbol!("color_trace_validation::absorb_j"))])
    );
    for suffix in ["", "_"] {
        let labels = (0..6)
            .map(|i| {
                Atom::var(symbol!(&format!(
                    "color_trace_validation::absorb_a{i}{suffix}"
                )))
            })
            .collect::<Vec<_>>();
        let x = Atom::var(symbol!(&format!(
            "color_trace_validation::absorb_x{suffix}"
        )));
        let z = Atom::var(symbol!("color_trace_validation::absorb_z"));
        let slot = |label: &Atom| coad.to_symbolic([label.clone()]);
        for (first, second) in [(0, 2), (0, 3), (3, 0), (4, 1)] {
            let f = idenso::color_f!(slot(&labels[first]), slot(&labels[second]), slot(&x));
            let mut sources = [false, true]
                .map(|adjoint| (word(adjoint, &labels), Vec::<Atom>::new()))
                .to_vec();
            if first < second {
                // Close the open chain M with T^z and with the identity.
                let chain = spenso::chain!(&line_in, &line_out; labels.iter().map(|label| idenso::color_t!(slot(label))));
                let probes = [
                    spenso::trace!(&fundamental; labels.iter().map(|label| idenso::color_t!(slot(label))).chain([idenso::color_t!(slot(&z))])),
                    spenso::trace!(&fundamental; labels.iter().map(|label| idenso::color_t!(slot(label)))),
                ];
                sources.push((chain, probes.to_vec()));
            }
            for (line, probes) in sources {
                let source = &line * &f * &spectator;
                let reduced = SymbolicTensor::infer(source.clone())
                    .unwrap()
                    .simplify_algebra(&settings())
                    .unwrap();
                assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
                assert!(
                    reduced
                        .expression()
                        .replace(spectator.to_pattern())
                        .with(Atom::Zero)
                        .is_zero(),
                    "the spectator is not kept factored: {}",
                    reduced.expression()
                );
                let colour = reduced
                    .expression()
                    .replace(spectator.to_pattern())
                    .with(Atom::num(1));
                // Close an open chain on both sides alike; a trace needs none.
                let closings: Vec<(Atom, Atom)> = if probes.is_empty() {
                    vec![(&line * &f, colour.clone())]
                } else {
                    let chain_head = symbol!("spenso::chain");
                    let closed = |expression: &Atom, extra: Option<&Atom>| {
                        expression.replace_map(|node, _, out| {
                            if let AtomView::Fun(function) = node
                                && function.get_symbol() == chain_head
                            {
                                let args = function.iter().skip(2).map(|arg| arg.to_owned());
                                **out = spenso::trace!(&fundamental; args.chain(extra.cloned()));
                            }
                        })
                    };
                    let probe = idenso::color_t!(slot(&z));
                    vec![
                        (&probes[0] * &f, closed(&colour, Some(&probe))),
                        (&probes[1] * &f, closed(&colour, None)),
                    ]
                };
                for (input, output) in closings {
                    let expected = components.compile(input.as_view());
                    let actual = components.compile(output.as_view());
                    let mut open = expected.open.clone();
                    open.sort();
                    let mut nonzero = false;
                    for sample in 0..40usize {
                        let mut bindings: Bindings = open
                            .iter()
                            .cloned()
                            .zip(sample_point(open.len(), sample))
                            .collect();
                        let value = expected.evaluate(&mut components, &mut bindings);
                        nonzero |= !value.is_zero();
                        assert_eq!(
                            actual.evaluate(&mut components, &mut bindings),
                            value,
                            "f({first},{second}) {suffix:?}: {output}"
                        );
                    }
                    assert!(
                        nonzero,
                        "f({first},{second}) {suffix:?}: trivial samples for {input}"
                    );
                }
            }
        }
    }
}

/// Tr(T^a T^b T^c T^d) f^{abx} f^{cdx}: each structure constant absorbs a
/// generator pair; the term-local repeat finishes it in one pass. Only
/// one-shot traces certify the result and skip the confirming no-op round.
#[test]
fn trace_with_two_structure_pairs_reduces_in_one_pass() {
    let mut components = Components::new();
    let names = [["a", "b", "c", "d", "x"], ["a_", "b_", "c_", "d_", "x_"]];
    for (names, one_shot) in names
        .into_iter()
        .flat_map(|names| [(names, true), (names, false)])
    {
        let [a, b, c, d, x] =
            names.map(|name| Atom::var(symbol!(&format!("color_trace_validation::pair_{name}"))));
        let slot = |label: &Atom| {
            ColorAdjoint {}
                .new_rep(symbol!("color_trace_validation::Na"))
                .to_symbolic([label.clone()])
        };
        let source = word(false, &[a.clone(), b.clone(), c.clone(), d.clone()])
            * idenso::color_f!(slot(&a), slot(&b), slot(&x))
            * idenso::color_f!(slot(&c), slot(&d), slot(&x));
        let reduced = reduced_scalar(
            "structure_pairs",
            &source,
            &AlgebraSettings {
                color: Some(if one_shot {
                    ColorSimplifySettings::default()
                } else {
                    ColorSimplifySettings::default().without_one_shot_traces()
                }),
                max_passes: one_shot.then_some(1),
                ..settings()
            },
        );
        let expected = components
            .compile(source.as_view())
            .evaluate(&mut components, &mut Bindings::new());
        assert!(!expected.is_zero());
        assert_eq!(
            components
                .compile(reduced.as_view())
                .evaluate(&mut components, &mut Bindings::new()),
            expected,
            "{names:?} one_shot={one_shot}: {reduced}"
        );
    }
}

/// d_R^{abx} f^{aik} f^{bjk} = (C_A/2) d_R^{ijx}, in both orientations.
#[test]
fn symmetric_trace_with_two_linked_structure_legs() {
    let mut components = Components::new();
    for (names, flipped) in [
        (["a", "b", "x", "i", "j", "k"], false),
        (["a", "b", "x", "i", "j", "k"], true),
        (["a_", "b_", "x_", "i_", "j_", "k_"], false),
    ] {
        let [a, b, x, i, j, k] =
            names.map(|name| Atom::var(symbol!(&format!("color_trace_validation::dff_{name}"))));
        let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
        let slot = |label: &Atom| coad.to_symbolic([label.clone()]);
        let symmetric = spenso::trace!(
            ColorFundamental {}
                .new_rep(symbol!("color_trace_validation::Nc"))
                .to_symbolic([]),
            shadowing::sym([&a, &b, &x].map(|label| idenso::color_t!(slot(label))))
        );
        let left = if flipped {
            idenso::color_f!(slot(&i), slot(&a), slot(&k))
        } else {
            idenso::color_f!(slot(&a), slot(&i), slot(&k))
        };
        let source = symmetric * left * idenso::color_f!(slot(&b), slot(&j), slot(&k));
        let reduced = SymbolicTensor::infer(source.clone())
            .unwrap()
            .simplify_algebra(&settings())
            .unwrap();
        assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
        assert!(
            !reduced.expression().contains_symbol(CS.f),
            "{}",
            reduced.expression()
        );
        let expected = components.compile(source.as_view());
        let actual = components.compile(reduced.expression().as_view());
        let mut nonzero = false;
        for values in [[0, 0, 7], [0, 1, 7], [3, 4, 7], [2, 5, 6], [7, 7, 7]] {
            let mut bindings = [&i, &j, &x].into_iter().cloned().zip(values).collect();
            let value = expected.evaluate(&mut components, &mut bindings);
            nonzero |= !value.is_zero();
            assert_eq!(
                actual.evaluate(&mut components, &mut bindings),
                value,
                "{names:?} flipped={flipped} {values:?}"
            );
        }
        assert!(nonzero);
    }
}

/// X^2 = sum_i x_i X for a colour sum X with open ports. The second case
/// has internal dummies in a summand, which the copies must not capture.
#[test]
fn colour_power_of_a_sum_reduces_to_a_scalar() {
    let mut components = Components::new();
    let [a, b, c, d, e, g] = ["a", "b", "c", "d", "e", "g"]
        .map(|name| Atom::var(symbol!(&format!("color_trace_validation::power_{name}"))));
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let slot = |label: &Atom| coad.to_symbolic([label.clone()]);
    let trace = word(false, &[a.clone(), b.clone(), c.clone()]);
    let triangle = idenso::color_f!(slot(&a), slot(&d), slot(&e))
        * idenso::color_f!(slot(&b), slot(&e), slot(&g))
        * idenso::color_f!(slot(&c), slot(&g), slot(&d));
    for (case, source) in [
        ("power_trace", &trace * &trace),
        ("power_internal", (&trace + &triangle).pow(Atom::num(2))),
    ] {
        assert!(matches!(source.as_view(), AtomView::Pow(_)));
        let reduced = reduced_scalar(case, &source, &settings());
        let expected = components
            .compile(source.as_view())
            .evaluate(&mut components, &mut Bindings::new());
        assert!(!expected.is_zero());
        assert_eq!(
            components
                .compile(reduced.as_view())
                .evaluate(&mut components, &mut Bindings::new()),
            expected,
            "{case}: {reduced}"
        );
    }
}

/// Powers of colour sums whose summands contract their own dummies, also in
/// metrics and nested sums (gate reproducers). The planner hands the kernel
/// the copies of such a power as a product of equal sums; distributing one
/// copy must not capture the other copy's dummies.
#[test]
fn colour_power_with_summand_local_dummies() {
    let mut components = Components::new();
    let [a, b, c, e, x, y] = ["a", "b", "c", "e", "x", "y"]
        .map(|name| Atom::var(symbol!(&format!("color_trace_validation::local_{name}"))));
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let slot = |label: &Atom| coad.to_symbolic([label.clone()]);
    let tr = |labels: &[&Atom]| {
        word(
            false,
            &labels
                .iter()
                .map(|label| (*label).clone())
                .collect::<Vec<_>>(),
        )
    };
    let g = |left: &Atom, right: &Atom| spenso::g!(slot(left), slot(right));
    let f = |p: &Atom, q: &Atom, r: &Atom| idenso::color_f!(slot(p), slot(q), slot(r));
    let cases = [
        (
            "trg_distributed",
            tr(&[&a, &c]) * tr(&[&b, &c]) + tr(&[&a, &c]) * g(&b, &c) + g(&a, &b),
        ),
        (
            "trg_noouter",
            tr(&[&a, &c]) * (tr(&[&b, &c]) + g(&b, &c)) + Atom::num(2) * tr(&[&a, &b]),
        ),
        (
            "nested_trg",
            tr(&[&a, &c]) * (tr(&[&b, &c]) + g(&b, &c)) + g(&a, &b),
        ),
        (
            "trg_nometric_inner",
            tr(&[&a, &c]) * (tr(&[&b, &c]) + Atom::num(2) * tr(&[&c, &b])) + g(&a, &b),
        ),
        (
            "nested",
            f(&a, &c, &e) * (f(&b, &c, &e) + tr(&[&b, &c, &e])) + g(&a, &b),
        ),
        (
            "nested_ffpartner",
            f(&a, &c, &e) * (f(&b, &c, &e) + tr(&[&b, &c, &e])) + f(&a, &x, &y) * f(&b, &x, &y),
        ),
    ];
    for (case, base) in cases {
        let source = base.pow(Atom::num(2));
        assert!(matches!(source.as_view(), AtomView::Pow(_)), "{case}");
        let reduced = reduced_scalar(case, &source, &settings());
        let expected = components
            .compile(source.as_view())
            .evaluate(&mut components, &mut Bindings::new());
        assert!(!expected.is_zero(), "{case}");
        assert_eq!(
            components
                .compile(reduced.as_view())
                .evaluate(&mut components, &mut Bindings::new()),
            expected,
            "{case}: {reduced}"
        );
    }
}

/// A summand-local dummy spelled like a label outside its sum is scoped to
/// the summand. Distributing the sum or giving the generated states canonical
/// labels must not capture it: the result equals that of the renamed input
/// and is itself admissible. Foreign vectors and pure colour.
#[test]
fn summand_local_dummies_shadowing_outer_labels() {
    let mut components = Components::new();
    let label = |name: &str| Atom::var(symbol!(&format!("color_trace_validation::shadow_{name}")));
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let slot = |name: &str| coad.to_symbolic([label(name)]);
    let p = |name: &str| {
        function!(
            spenso::tensor_symbol!("color_trace_validation::shadow_P"),
            slot(name)
        )
    };
    let q = |name: &str| {
        function!(
            spenso::tensor_symbol!("color_trace_validation::shadow_Q"),
            slot(name)
        )
    };
    let f = |x: &str, y: &str, z: &str| idenso::color_f!(slot(x), slot(y), slot(z));
    let cycle = Atom::mul_many([
        f("e0", "x1", "x0"),
        f("e1", "x1", "x2"),
        f("x0", "x4", "e4"),
        f("x2", "e2", "x3"),
        f("x3", "e3", "x4"),
    ]);
    let five = word(true, &["a", "b", "c", "e", "g"].map(label));
    let cases = [
        (
            "foreign_vectors",
            p("a") + p("b") * q("a") * q("b"),
            p("a") + p("w") * q("a") * q("w"),
            five,
        ),
        (
            "colour_bubble",
            spenso::g!(slot("e0"), slot("z")) + f("e0", "e1", "k") * f("z", "e1", "k"),
            spenso::g!(slot("e0"), slot("z")) + f("e0", "w", "k") * f("z", "w", "k"),
            cycle,
        ),
    ];
    for (case, shadowed, renamed, line) in cases {
        let reduced = SymbolicTensor::infer(&shadowed * &line)
            .unwrap()
            .simplify_algebra(&settings())
            .unwrap();
        assert_eq!(
            reduced.reduction_status(),
            ReductionStatus::Complete,
            "{case}"
        );
        SymbolicTensor::infer(reduced.expression().clone())
            .unwrap_or_else(|error| panic!("{case}: inadmissible output: {error}"));
        // The oracle scopes no dummy to a summand: a sum the output keeps
        // factored is compared in its renamed spelling.
        let output = reduced
            .expression()
            .replace(shadowed.to_pattern())
            .with(renamed.to_pattern());
        let expected = components.compile((&renamed * &line).as_view());
        let actual = components.compile(output.as_view());
        let mut open = expected.open.clone();
        open.sort();
        let mut nonzero = false;
        for sample in 0..24usize {
            let mut bindings: Bindings = open
                .iter()
                .cloned()
                .zip(sample_point(open.len(), sample))
                .collect();
            let value = expected.evaluate(&mut components, &mut bindings);
            nonzero |= !value.is_zero();
            assert_eq!(
                actual.evaluate(&mut components, &mut bindings),
                value,
                "{case} sample={sample}: {output}"
            );
        }
        assert!(nonzero, "{case}: all samples trivial");
    }
}

/// A power of a tensor contracts its copies: a closed scope. Its labels may
/// also be spelled by a port of the summand holding it, which a factor the
/// sum is distributed into contracts. The reduction keeps the scopes apart.
/// The oracle evaluates the input with the power's labels renamed; plain and
/// wildcard-spelled labels.
#[test]
fn power_scopes_stay_apart_from_ports() {
    let mut components = Components::new();
    for suffix in ["", "_"] {
        let label = |name: &str| {
            Atom::var(symbol!(&format!(
                "color_trace_validation::scope_{name}{suffix}"
            )))
        };
        let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
        let slot = |name: &str| coad.to_symbolic([label(name)]);
        let tr = |names: &[&str]| {
            word(
                false,
                &names.iter().map(|name| label(name)).collect::<Vec<_>>(),
            )
        };
        let f = |x: &str, y: &str, z: &str| idenso::color_f!(slot(x), slot(y), slot(z));
        let g = |x: &str, y: &str| spenso::g!(slot(x), slot(y));
        let nc = Atom::var(symbol!("color_trace_validation::Nc"));
        let na = Atom::var(symbol!("color_trace_validation::Na"));
        let two = Atom::num(2);
        let powers = |c: &str, d: &str| tr(&[c, "x"]).pow(&two) + &nc * tr(&[d, "y"]).pow(&two);
        let cases = [
            (
                "two_powers_f",
                powers("c", "c") * tr(&["c", "e", "g"]) * f("c", "g", "e"),
                powers("u", "w") * tr(&["c", "e", "g"]) * f("c", "g", "e"),
            ),
            (
                "two_powers_tr",
                powers("c", "c") * tr(&["c", "e", "g"]) * tr(&["c", "g", "e"]),
                powers("u", "w") * tr(&["c", "e", "g"]) * tr(&["c", "g", "e"]),
            ),
            (
                "power_beside_port",
                (tr(&["c", "x"]).pow(&two) + &na) * tr(&["c", "e", "g"]) * f("c", "g", "e"),
                (tr(&["u", "x"]).pow(&two) + &na) * tr(&["c", "e", "g"]) * f("c", "g", "e"),
            ),
            (
                "power_of_sum_beside_port",
                ((tr(&["c", "x"]) + g("c", "x")).pow(&two) + &na)
                    * tr(&["c", "e", "g"])
                    * tr(&["c", "g", "e"]),
                ((tr(&["u", "x"]) + g("u", "x")).pow(&two) + &na)
                    * tr(&["c", "e", "g"])
                    * tr(&["c", "g", "e"]),
            ),
        ];
        for (case, shadowed, renamed) in cases {
            let case = format!("{case}{suffix}");
            let reduced = reduced_scalar(&case, &shadowed, &settings());
            let expected = components
                .compile(renamed.as_view())
                .evaluate(&mut components, &mut Bindings::new());
            assert!(!expected.is_zero(), "{case}");
            assert_eq!(
                components
                    .compile(reduced.as_view())
                    .evaluate(&mut components, &mut Bindings::new()),
                expected,
                "{case}: {reduced}"
            );
        }
    }
}

/// A foreign tensor with a free adjoint index beside a decomposed trace
/// stays a top-level factor of the result in both trace modes: no colour
/// identity reaches it, so it is never part of a colour row. Exact at SU(3).
#[test]
fn foreign_free_adjoint_factor_stays_factored() {
    let mut components = Components::new();
    let label = |name: &str| Atom::var(symbol!(&format!("color_trace_validation::free_{name}")));
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let vector = function!(
        spenso::tensor_symbol!("color_trace_validation::shadow_P"),
        coad.to_symbolic([label("z")])
    );
    for (case, line) in [
        (
            "adjoint_five",
            word(true, &["a", "b", "c", "e", "g"].map(label)),
        ),
        (
            "fundamental_block",
            word(false, &["a", "b", "c", "a", "b", "e", "g"].map(label)),
        ),
    ] {
        let source = &vector * &line;
        let reduced = SymbolicTensor::infer(source.clone())
            .unwrap()
            .simplify_algebra(&settings())
            .unwrap();
        assert_eq!(
            reduced.reduction_status(),
            ReductionStatus::Complete,
            "{case}"
        );
        assert_kept_as_factors(case, reduced.expression(), &[&vector]);
        let expected = components.compile(source.as_view());
        let actual = components.compile(reduced.expression().as_view());
        let mut open = expected.open.clone();
        open.sort();
        let mut nonzero = false;
        for sample in 0..24usize {
            let mut bindings: Bindings = open
                .iter()
                .cloned()
                .zip(sample_point(open.len(), sample))
                .collect();
            let value = expected.evaluate(&mut components, &mut bindings);
            nonzero |= !value.is_zero();
            assert_eq!(
                actual.evaluate(&mut components, &mut bindings),
                value,
                "{case} sample={sample}"
            );
        }
        assert!(nonzero, "{case}: all samples trivial");
    }
}

/// Each kept tensor occurs once in `reduced`, as a factor of the whole
/// product, of a term of the sum, or of a term of a sum factor: it was
/// neither distributed nor rewritten.
fn assert_kept_as_factors(case: &str, reduced: &Atom, kept: &[&Atom]) {
    fn reaches(expression: AtomView<'_>, kept: AtomView<'_>) -> bool {
        match expression {
            AtomView::Add(sum) => sum.iter().any(|term| reaches(term, kept)),
            AtomView::Mul(product) => product.iter().any(|factor| {
                factor == kept || matches!(factor, AtomView::Add(_)) && reaches(factor, kept)
            }),
            _ => expression == kept,
        }
    }
    for tensor in kept {
        let AtomView::Fun(function) = tensor.as_view() else {
            panic!("{case}: kept {tensor} is not a tensor")
        };
        let mut occurrences = 0;
        reduced.visitor(&mut |node| {
            occurrences += usize::from(
                matches!(node, AtomView::Fun(other) if other.get_symbol() == function.get_symbol()),
            );
            true
        });
        assert_eq!(
            occurrences, 1,
            "{case}: {tensor} distributed or rewritten: {reduced}"
        );
        assert!(
            reaches(reduced.as_view(), tensor.as_view()),
            "{case}: {tensor} is not a factor: {reduced}"
        );
    }
}

/// Foreign tensors with their own free adjoint index beside closed colour
/// cores whose sums the reduction distributes stay factors as found, in
/// settings (i), (ii) and (iii) and in both trace modes: a top-level factor,
/// or a factor of its own term of a sum. The cores are f(a,b,c) times a sum
/// of f's, FK0945's four-gluon colour sum, an eight-f vacuum graph, Tr_F Tr_A
/// of length four, and a power of a trace beside a port. The input values
/// multiply the cores' exact SU(3) components by the foreign tensors'.
#[test]
fn foreign_factors_stay_factored_beside_closed_colour() {
    let mut components = Components::new();
    let label = |name: &str| Atom::var(symbol!(&format!("color_trace_validation::keep_{name}")));
    let foreign_v = spenso::tensor_symbol!("color_trace_validation::keep_V");
    let foreign_u = spenso::tensor_symbol!("color_trace_validation::keep_U");
    let two = Atom::num(2);
    for setting in ["i", "ii", "iii"] {
        let cof = setting == "iii";
        let (adjoint_dimension, fundamental_dimension): (
            spenso::structure::dimension::Dimension,
            spenso::structure::dimension::Dimension,
        ) = if cof {
            (8usize.into(), 3usize.into())
        } else {
            (
                symbol!("color_trace_validation::Na").into(),
                symbol!("color_trace_validation::Nc").into(),
            )
        };
        let coad = ColorAdjoint {}.new_rep(adjoint_dimension);
        let cof_rep = ColorFundamental {}.new_rep(fundamental_dimension);
        let na = if cof {
            Atom::num(8)
        } else {
            Atom::var(symbol!("color_trace_validation::Na"))
        };
        let slot = |name: &str| coad.to_symbolic([label(name)]);
        let f = |a: &str, b: &str, c: &str| idenso::color_f!(slot(a), slot(b), slot(c));
        let line = |adjoint: bool, names: &[&str]| {
            let representation = if adjoint {
                coad.to_symbolic([])
            } else {
                cof_rep.to_symbolic([])
            };
            spenso::trace!(representation; names.iter().map(|name| {
                if adjoint {
                    idenso::color_f!(Atom::var(SPENSO_TAG.chain_in), Atom::var(SPENSO_TAG.chain_out), slot(name))
                } else {
                    idenso::color_t!(slot(name))
                }
            }))
        };
        let v = |k: usize| label(&format!("v{k}"));
        let product = |factors: Vec<Atom>| Atom::mul_many(factors.iter());
        // f(a,b,c) * (v0 f(a,b,c) + v1 f(a,c,b)): the minimal reproducer.
        let m1 = f("a", "b", "c") * (v(0) * f("a", "b", "c") + v(1) * f("a", "c", "b"));
        // FK0945: labels in an order whose prefixes close structure constants
        // early, which keeps the input oracle sparse.
        let g4 = product(vec![
            f("g01", "g02", "g03"),
            f("g01", "g04", "g05"),
            f("g02", "g06", "g07"),
            f("g03", "g08", "g09"),
            f("g10", "g11", "g06"),
            f("g10", "g11", "g08"),
        ]) * (v(0) * f("g12", "g07", "g05") * f("g12", "g09", "g04")
            + v(1) * f("g12", "g07", "g04") * f("g12", "g09", "g05")
            + v(2) * f("g12", "g07", "g09") * f("g12", "g04", "g05"));
        let closed8 = product(
            [
                [1, 3, 5],
                [1, 7, 9],
                [3, 19, 23],
                [5, 15, 21],
                [7, 11, 23],
                [9, 13, 17],
                [11, 13, 15],
                [17, 19, 21],
            ]
            .map(|[a, b, c]| {
                f(
                    &format!("h{a:02}"),
                    &format!("h{b:02}"),
                    &format!("h{c:02}"),
                )
            })
            .to_vec(),
        );
        let mixed44 = line(false, &["a", "b", "c", "d"]) * line(true, &["a", "b", "c", "d"]);
        // The power contracts its copies: its labels are a closed scope, which
        // the oracle evaluates renamed.
        let powport = (line(false, &["c", "x"]).pow(&two) + &na)
            * line(false, &["c", "e", "g"])
            * f("c", "g", "e");
        let powport_renamed = (line(false, &["u", "x"]).pow(&two) + &na)
            * line(false, &["c", "e", "g"])
            * f("c", "g", "e");
        let kept_v = function!(foreign_v, slot("zz"));
        let kept_u = function!(foreign_u, slot("zz"));
        let spectator = Atom::num(1) + label("x");
        let value = |components: &mut Components, core: &Atom| {
            components
                .compile(core.as_view())
                .evaluate(components, &mut Bindings::new())
        };
        let cores = [
            ("m1", m1.clone(), m1.clone()),
            ("fk0945", g4.clone(), g4),
            ("closed8", closed8.clone(), closed8),
            ("mixed44", mixed44.clone(), mixed44),
            ("powport", powport, powport_renamed),
        ];
        let ff = f("a", "b", "c") * f("a", "b", "c");
        let mut cases = Vec::new();
        for (name, core, renamed) in cores {
            let expected = value(&mut components, &renamed);
            assert!(!expected.is_zero(), "{name}");
            cases.push((
                format!("{name}xV"),
                &kept_v * &core,
                vec![(kept_v.clone(), expected)],
            ));
        }
        // Each term of a sum keeps its own foreign factor, also behind a
        // scalar spectator of the sum.
        let terms = &kept_v * &m1 + &kept_u * &ff;
        let term_values = vec![
            (kept_v.clone(), value(&mut components, &m1)),
            (kept_u.clone(), value(&mut components, &ff)),
        ];
        cases.push(("m1xV+ffxU".into(), terms.clone(), term_values.clone()));
        cases.push((
            "(1+x)(m1xV+ffxU)".into(),
            &spectator * &terms,
            term_values
                .into_iter()
                .map(|(tensor, core)| (tensor, scalar(core * Atom::num(6))))
                .collect(),
        ));
        let settings = match setting {
            "i" => settings(),
            "ii" => AlgebraSettings {
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            },
            _ => AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            },
        };
        for (case, source, kept) in cases {
            let case = format!("{case}/{setting}");
            let admitted = SymbolicTensor::infer(source.clone()).unwrap();
            let reduced = admitted.simplify_algebra(&settings).unwrap();
            assert_eq!(
                reduced.reduction_status(),
                ReductionStatus::Complete,
                "{case}: {}",
                reduced.expression()
            );
            assert_eq!(admitted.structure(), reduced.structure(), "{case}");
            assert_kept_as_factors(
                &case,
                reduced.expression(),
                &kept.iter().map(|(tensor, _)| tensor).collect::<Vec<_>>(),
            );
            let actual = components.compile(reduced.expression().as_view());
            assert_eq!(actual.open, vec![label("zz")], "{case}");
            let mut nonzero = false;
            for zz in 0..8 {
                let mut bindings = Bindings::from([(label("zz"), zz)]);
                let expected = scalar(Atom::add_many(kept.iter().map(|(tensor, core)| {
                    core * components
                        .compile(tensor.as_view())
                        .evaluate(&mut components, &mut bindings.clone())
                })));
                nonzero |= !expected.is_zero();
                assert_eq!(
                    actual.evaluate(&mut components, &mut bindings),
                    expected,
                    "{case} zz={zz}: {}",
                    reduced.expression()
                );
            }
            assert!(nonzero, "{case}: all samples trivial");
        }
        // A foreign factor spelled with a counter label d_0, alone or as a
        // contracted pair, beside an open adjoint trace whose states take
        // canonical dummy labels: those labels avoid d_0. Open cores beside
        // a foreign tensor stay Deferred under full contraction, as before.
        let counter = coad.to_symbolic([AbstractIndex::Dummy(0).to_atom()]);
        let open_core = line(true, &["a0", "a1", "a2", "a3", "a4"]);
        for (case, kept) in [
            ("adj5xV(d_0)", vec![function!(foreign_v, &counter)]),
            (
                "adj5xV(d_0)U(d_0)",
                vec![
                    function!(foreign_v, &counter),
                    function!(foreign_u, &counter),
                ],
            ),
        ] {
            let case = format!("{case}/{setting}");
            let source = product(kept.clone()) * &open_core;
            let admitted = SymbolicTensor::infer(source.clone()).unwrap();
            let reduced = admitted.simplify_algebra(&settings).unwrap();
            let status = reduced.reduction_status();
            assert!(
                status == ReductionStatus::Complete
                    || (setting != "i" && status == ReductionStatus::Deferred),
                "{case}: {status:?}"
            );
            assert_eq!(admitted.structure(), reduced.structure(), "{case}");
            assert_kept_as_factors(
                &case,
                reduced.expression(),
                &kept.iter().collect::<Vec<_>>(),
            );
            let explicit = reduced.undo_chain().unwrap();
            let actual = components.compile(explicit.expression().as_view());
            let expected = components.compile(source.as_view());
            let mut open = expected.open.clone();
            open.sort();
            let mut actual_open = actual.open.clone();
            actual_open.sort();
            assert_eq!(actual_open, open, "{case}");
            let mut nonzero = false;
            for sample in 0..16usize {
                let mut bindings: Bindings = open
                    .iter()
                    .cloned()
                    .zip(sample_point(open.len(), sample))
                    .collect();
                let value = expected.evaluate(&mut components, &mut bindings);
                nonzero |= !value.is_zero();
                assert_eq!(
                    actual.evaluate(&mut components, &mut bindings),
                    value,
                    "{case} sample={sample}: {}",
                    reduced.expression()
                );
            }
            assert!(nonzero, "{case}: all samples trivial");
        }
    }
}

/// Wildcard-spelled labels are reduced on private aliases. An alias never
/// names a label the input already uses, whichever spelling that label has.
#[test]
fn wildcard_aliases_never_collide_with_input_labels() {
    let mut components = Components::new();
    let coad = ColorAdjoint {}.new_rep(symbol!("color_trace_validation::Na"));
    let wild = symbol!("color_trace_validation::alias_a_");
    for alias in [
        format!("idenso::wildcard_label_{}", wild.get_id()),
        format!("idenso::wildcard_label_{}_0", wild.get_id()),
    ] {
        let slot = |label: Atom| coad.to_symbolic([label]);
        let [w, c] = ["alias_w", "alias_c"].map(|name| {
            slot(Atom::var(symbol!(&format!(
                "color_trace_validation::{name}"
            ))))
        });
        let a = slot(Atom::var(wild));
        let l = slot(Atom::var(symbol!(&alias)));
        let source = idenso::color_f!(&a, &w, &c)
            * idenso::color_f!(&l, &w, &c)
            * function!(
                spenso::tensor_symbol!("color_trace_validation::alias_P"),
                &a
            )
            * function!(
                spenso::tensor_symbol!("color_trace_validation::alias_Q"),
                &l
            );
        let reduced = SymbolicTensor::infer(source.clone())
            .unwrap()
            .simplify_algebra(&settings())
            .unwrap();
        assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
        let expected = components
            .compile(source.as_view())
            .evaluate(&mut components, &mut Bindings::new());
        assert!(!expected.is_zero());
        assert_eq!(
            components
                .compile(reduced.expression().as_view())
                .evaluate(&mut components, &mut Bindings::new()),
            expected,
            "{alias}: {}",
            reduced.expression()
        );
    }
}

/// d_A^{abcd} d_F^{abcd} = N (N^2-1)(N^2+6)/48: the mixed quartic Gram
/// invariant substitutes like d44(A,A) and d44(F,F) (SU(3): 15/2).
#[test]
fn mixed_quartic_gram_substitutes_at_fixed_dimension() {
    let mut components = Components::new();
    let coad = ColorAdjoint {}.new_rep(8);
    let cof = ColorFundamental {}.new_rep(3);
    let labels = ["a", "b", "c", "d"]
        .map(|name| Atom::var(symbol!(&format!("color_trace_validation::gram_{name}"))));
    let slots = labels.clone().map(|label| coad.to_symbolic([label]));
    let source = spenso::trace!(
        coad.to_symbolic([]),
        shadowing::sym(slots.iter().map(|slot| idenso::color_f!(
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            slot
        )))
    ) * spenso::trace!(
        cof.to_symbolic([]),
        shadowing::sym(slots.iter().map(|slot| idenso::color_t!(slot)))
    );
    let reduced = SymbolicTensor::infer(source.clone())
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.expression(), &(Atom::num(15) / Atom::num(2)));
    let expected = components
        .compile(source.as_view())
        .evaluate(&mut components, &mut Bindings::new());
    assert_eq!(reduced.expression(), &expected);
}

/// color.h's two-block identity, Tr(X^a X^b W X^a X^b R) =
/// Tr(X^a X^b W X^b X^a R) + κ² (C_A/2) Tr(X^a W X^a R), for T and raw F,
/// plain and wildcard-spelled labels.
#[test]
fn repeated_generator_block_is_reversed() {
    let mut components = Components::new();
    for adjoint in [false, true] {
        for names in [["a", "b", "c", "d", "e"], ["a_", "b_", "c_", "d_", "e_"]] {
            let [a, b, c, d, e] = names
                .map(|name| Atom::var(symbol!(&format!("color_trace_validation::block_{name}"))));
            // Raw adjoint words need four open indices for nonzero samples.
            let g = Atom::var(symbol!(&format!(
                "color_trace_validation::block_g{}",
                &names[0][1..]
            )));
            let open = |word: &[&Atom]| {
                let mut word = word
                    .iter()
                    .map(|label| (*label).clone())
                    .collect::<Vec<_>>();
                if adjoint {
                    word.extend([e.clone(), g.clone()]);
                }
                word
            };
            // Raw adjoint F^a F^b X F^a F^b W summed over a, b vanishes for
            // the first word; the second is nonzero.
            let mut nonzero = false;
            for source in [
                word(adjoint, &open(&[&a, &b, &c, &a, &b, &d])),
                word(adjoint, &open(&[&a, &b, &c, &d, &a, &b])),
            ] {
                let reduced = SymbolicTensor::infer(source.clone())
                    .unwrap()
                    .simplify_algebra(&settings())
                    .unwrap();
                assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
                let expected = components.compile(source.as_view());
                let actual = components.compile(reduced.expression().as_view());
                let mut open = expected.open.clone();
                open.sort();
                for sample in 0..24usize {
                    let mut bindings: Bindings = open
                        .iter()
                        .cloned()
                        .zip(sample_point(open.len(), sample))
                        .collect();
                    let value = expected.evaluate(&mut components, &mut bindings);
                    nonzero |= !value.is_zero();
                    assert_eq!(
                        actual.evaluate(&mut components, &mut bindings),
                        value,
                        "adjoint={adjoint} {names:?} sample={sample}: {}",
                        reduced.expression()
                    );
                }
            }
            assert!(nonzero, "adjoint={adjoint} {names:?}");
        }
    }
}

/// Structured and hashed points over the 8 colours. SU(3) components are
/// sparse; repeating each colour in pairs (shuffled) keeps more of them
/// nonzero.
fn sample_point(count: usize, sample: usize) -> Vec<usize> {
    let hash = |seed: u64| {
        let mut z = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 31)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 29)
    };
    let mut values = (0..count)
        .map(|i| match sample {
            0 => 7,
            1 => [1, 2][i % 2],
            2 => [0, 1, 2][i % 3],
            3 => [3, 4, 5, 6][i % 4],
            4..=7 => (hash((sample * 64 + i / 2) as u64) % 8) as usize,
            _ => (hash((sample * 64 + i) as u64) % 8) as usize,
        })
        .collect::<Vec<_>>();
    for i in (1..values.len()).rev() {
        let j = (hash((sample * 4096 + i) as u64) % (i as u64 + 1)) as usize;
        values.swap(i, j);
    }
    values
}

/// The open cases of `color_trace_scaling matrix`, reduced in settings (i)
/// and (iii) and compared with the input oracle at sample components.
/// Slow; run explicitly with `--ignored`.
#[test]
#[ignore = "slow component oracle over every matrix case; run with --ignored"]
fn matrix_cases_match_component_oracle() {
    let mut components = Components::new();
    let label = |name: &str| Atom::var(symbol!(&format!("color_trace_validation::m_{name}")));
    let names = |prefix: &str, length: usize| {
        (0..length)
            .map(|i| format!("{prefix}{i}"))
            .collect::<Vec<_>>()
    };
    let repeated = ["r0", "r1", "r2", "r3", "r0", "r4", "r5", "r6"].map(String::from);
    // (name, adjoint line, line word, attached f, second line sharing ports)
    type MatrixCase = (
        &'static str,
        bool,
        Vec<String>,
        Option<[&'static str; 3]>,
        Option<[&'static str; 4]>,
    );
    let cases: [MatrixCase; 10] = [
        ("fund7", false, names("a", 7), None, None),
        ("adj7", true, names("a", 7), None, None),
        ("fund8", false, names("a", 8), None, None),
        ("adj8", true, names("a", 8), None, None),
        ("fundrep8", false, repeated.to_vec(), None, None),
        ("adjrep8", true, repeated.to_vec(), None, None),
        (
            "emb_tf_fund7",
            false,
            names("a", 7),
            Some(["a0", "a3", "x"]),
            None,
        ),
        (
            "emb_tf_adj7",
            true,
            names("a", 7),
            Some(["a0", "a3", "x"]),
            None,
        ),
        (
            "emb_tt_fund",
            false,
            names("a", 6),
            None,
            Some(["a0", "a3", "b0", "b1"]),
        ),
        (
            "emb_tt_adj",
            true,
            names("a", 6),
            None,
            Some(["a0", "a3", "b0", "b1"]),
        ),
    ];
    let only = std::env::var("COLOR_MATRIX_CASES").ok();
    for (case, adjoint, word_labels, structure, second) in cases {
        if only
            .as_ref()
            .is_some_and(|only| !only.split(',').any(|name| name == case))
        {
            continue;
        }
        for cof in [false, true] {
            let (adjoint_dimension, fundamental_dimension): (
                spenso::structure::dimension::Dimension,
                spenso::structure::dimension::Dimension,
            ) = if cof {
                (8usize.into(), 3usize.into())
            } else {
                (
                    symbol!("color_trace_validation::Na").into(),
                    symbol!("color_trace_validation::Nc").into(),
                )
            };
            let coad = ColorAdjoint {}.new_rep(adjoint_dimension);
            let slot = |name: &str| coad.to_symbolic([label(name)]);
            let line = |labels: &[String]| {
                let representation = if adjoint {
                    coad.to_symbolic([])
                } else {
                    ColorFundamental {}
                        .new_rep(fundamental_dimension)
                        .to_symbolic([])
                };
                spenso::trace!(representation; labels.iter().map(|name| {
                    if adjoint {
                        idenso::color_f!(Atom::var(SPENSO_TAG.chain_in), Atom::var(SPENSO_TAG.chain_out), slot(name))
                    } else {
                        idenso::color_t!(slot(name))
                    }
                }))
            };
            let mut source = line(&word_labels);
            if let Some([a, b, c]) = structure {
                source *= idenso::color_f!(slot(a), slot(b), slot(c));
            }
            if let Some(second) = second {
                source *= line(&second.map(String::from));
            }
            let settings = if cof {
                AlgebraSettings {
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                }
            } else {
                settings()
            };
            let case_name = format!("matrix/{case}/{}", if cof { "iii" } else { "i" });
            let admitted = SymbolicTensor::infer(source.clone()).unwrap();
            let reduced = admitted.simplify_algebra(&settings).unwrap();
            assert_eq!(
                reduced.reduction_status(),
                ReductionStatus::Complete,
                "{case_name}"
            );
            let explicit = reduced.undo_chain().unwrap();
            let actual = components.compile(explicit.expression().as_view());
            let expected = components.compile(source.as_view());
            let mut open = actual.open.clone();
            open.sort();
            let mut source_open = expected.open.clone();
            source_open.sort();
            assert_eq!(open, source_open, "{case_name}");
            let mut nonzero = 0;
            for sample in 0..12usize {
                let mut bindings: Bindings = open
                    .iter()
                    .cloned()
                    .zip(sample_point(open.len(), sample))
                    .collect();
                let value = expected.evaluate(&mut components, &mut bindings);
                nonzero += usize::from(!value.is_zero());
                assert_eq!(
                    actual.evaluate(&mut components, &mut bindings),
                    value,
                    "{case_name} sample={sample}"
                );
            }
            assert!(nonzero > 0, "{case_name}: all samples trivial");
        }
    }
}
