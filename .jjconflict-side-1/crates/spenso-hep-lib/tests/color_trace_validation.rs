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
    structure::{abstract_index::AbstractIndex, representation::RepName},
    tensors::data::GetTensorData,
};
use spenso_hep_lib::{su3_generator_data_atom, su3_structure_f_data_atom};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, parse_lit, symbol,
};

type Matrix = Vec<Vec<Atom>>;
type Bindings = BTreeMap<Atom, usize>;

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

    fn compile(&mut self, value: AtomView<'_>) -> Component {
        let kind = match value {
            AtomView::Num(_) => Kind::Scalar(value.to_owned()),
            AtomView::Var(variable) => {
                let name = variable.get_symbol();
                Kind::Scalar(Atom::num(
                    if name == symbol!("color_trace_validation::Nc") {
                        3
                    } else if name == symbol!("color_trace_validation::Na") {
                        8
                    } else {
                        panic!("unassigned scalar {value}")
                    },
                ))
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
                } else {
                    panic!("unsupported component node {value}")
                }
            }
        };
        Component::new(kind)
    }
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
            Kind::Scalar(_) | Kind::Structure(_) | Kind::Metric(_) => true,
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
    let order = if reverse {
        labels.iter().rev().cloned().collect()
    } else {
        labels.clone()
    };
    let source = SymbolicTensor::infer(word(adjoint, &order)).unwrap();
    let reduced = source.simplify_algebra(&settings()).unwrap();
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
    let source = SymbolicTensor::infer(word(adjoint, &word_labels)).unwrap();
    let reduced = source.simplify_algebra(&settings()).unwrap();
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
    let expected = {
        scalar(Atom::add_many((0..8).map(|a| {
            components.trace(adjoint, false, vec![a, 7, 7, 7, a, 7, 7, 7])
        })))
    };
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
        { components.trace(true, false, values.clone()) * Atom::num(if flipped { -1 } else { 1 }) };
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
