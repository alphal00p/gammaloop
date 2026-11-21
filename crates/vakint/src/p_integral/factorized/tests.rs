//! Hermetic contract tests: the two fixtures are records written by the text
//! writer from the beta thread's L1 DOT files (`diagram-000000.dot`, gluon
//! loop; `diagram-000002.dot`, quark loop; colour mode `Invariants`, metric
//! projector). Every negative test mutates one of them and expects its own
//! error.
use super::*;
use symbolica::function;

const GLUON: &str = include_str!("fixtures/l1-gluon-bubble.txt");
const QUARK: &str = include_str!("fixtures/l1-quark-bubble.txt");

fn read(text: &str) -> FactorizedNumerator {
    crate::Vakint::initialize_vakint_symbols();
    FactorizedNumerator::from_text(text).unwrap_or_else(|e| panic!("{e}\n{text}"))
}

fn slot(rep: &str, dim: Atom, index: Atom) -> Atom {
    use spenso::structure::representation::RepName;
    match rep {
        names::MINK => spenso::structure::representation::Minkowski {}.to_symbolic([dim, index]),
        names::BIS => idenso::representations::Bispinor {}.to_symbolic([dim, index]),
        names::COAD => idenso::representations::ColorAdjoint {}.to_symbolic([dim, index]),
        _ => unreachable!(),
    }
}

fn hedge(h: i64) -> Atom {
    function!(feynkit_graph::symbols::hedge_index(), h, 1)
}

fn mink(value: &FactorizedNumerator, h: i64) -> Atom {
    slot(names::MINK, Atom::var(value.dimension), hedge(h))
}

#[test]
fn the_fixtures_validate_and_round_trip_byte_identically() {
    for text in [GLUON, QUARK] {
        let value = read(text);
        value.validate().unwrap();
        assert_eq!(value.to_text().unwrap(), text);
        let again = FactorizedNumerator::from_text(&value.to_text().unwrap()).unwrap();
        assert_eq!(again, value);
    }
    let both = format!("{GLUON}{QUARK}");
    assert_eq!(parse_factorized_records(&both).unwrap().len(), 2);
    let quark = read(QUARK);
    assert_eq!(quark.closed_loops().unwrap().len(), 1);
    assert_eq!(quark.prefactor.loop_sign, -1);
    assert_eq!(quark.prefactor.nf_power, 1);
    assert!(read(GLUON).closed_loops().unwrap().is_empty());
}

#[test]
fn an_index_contracted_twice_is_rejected() {
    // Three occurrences in one factor.
    let mut value = read(GLUON);
    let edge = value
        .edges
        .iter()
        .find(|e| !e.is_external())
        .unwrap()
        .clone();
    let h = edge.source_hedge.unwrap() as i64;
    let dimension = value.dimension;
    let q = move |e: usize, h: i64| {
        function!(
            feynkit_momentum(),
            e as i64,
            slot(names::MINK, Atom::var(dimension), hedge(h))
        )
    };
    let extra = q(edge.id, h) * q(edge.id, h);
    value
        .edges
        .iter_mut()
        .find(|e| e.id == edge.id)
        .unwrap()
        .factor = &edge.factor * extra;
    assert!(
        matches!(
            value.validate(),
            Err(FactorizedError::IndexOverContracted { .. })
        ),
        "{:?}",
        value.validate()
    );
    // Contracted inside the edge factor AND with the vertex.
    let mut value = read(GLUON);
    let extra = q(edge.id, h);
    value
        .edges
        .iter_mut()
        .find(|e| e.id == edge.id)
        .unwrap()
        .factor = &edge.factor * extra;
    assert!(
        matches!(
            value.validate(),
            Err(FactorizedError::DummyIndexShared { .. })
        ),
        "{:?}",
        value.validate()
    );
}

fn feynkit_momentum() -> Symbol {
    feynkit_graph::momentum_symbol()
}

#[test]
fn an_uncontracted_index_is_rejected() {
    // The gluon propagator g(μ_s, μ_t) replaced by Q(e, μ_s): μ_t dangles.
    let mut value = read(GLUON);
    let edge = value
        .edges
        .iter()
        .find(|e| !e.is_external())
        .unwrap()
        .clone();
    let h = edge.source_hedge.unwrap() as i64;
    let replacement = function!(feynkit_momentum(), edge.id as i64, mink(&value, h));
    value
        .edges
        .iter_mut()
        .find(|e| e.id == edge.id)
        .unwrap()
        .factor = replacement;
    assert!(
        matches!(
            value.validate(),
            Err(FactorizedError::IndexUncontracted { .. })
        ),
        "{:?}",
        value.validate()
    );
}

#[test]
fn broken_momentum_conservation_is_rejected() {
    let mut value = read(QUARK);
    let edge = value.edges.iter_mut().find(|e| !e.is_external()).unwrap();
    edge.momentum[0] += 1;
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::MomentumNotConserved { .. } | FactorizedError::MomentumBasis { .. })
    ));
    // Conservation broken without touching a basis edge: the external leg.
    let mut value = read(QUARK);
    let external = value
        .edges
        .iter()
        .filter(|e| e.is_external())
        .map(|e| e.id)
        .max()
        .unwrap();
    let edge = value.edges.iter_mut().find(|e| e.id == external).unwrap();
    edge.momentum[1] = -edge.momentum[1];
    assert!(
        matches!(
            value.validate(),
            Err(FactorizedError::MomentumNotConserved { .. })
        ),
        "{:?}",
        value.validate()
    );
}

#[test]
fn a_colour_tensor_left_in_a_factor_is_rejected() {
    let mut value = read(GLUON);
    let adjoint = |h: i64| slot(names::COAD, Atom::num(8), hedge(h));
    let colour = function!(
        idenso::color::CS.f,
        adjoint(100),
        adjoint(101),
        adjoint(102)
    ) * function!(
        idenso::color::CS.f,
        adjoint(100),
        adjoint(101),
        adjoint(102)
    );
    value.vertices[0].factor = &value.vertices[0].factor * colour;
    assert!(
        matches!(
            value.validate(),
            Err(FactorizedError::ColourInFactor { .. })
        ),
        "{:?}",
        value.validate()
    );
}

#[test]
fn a_gamma5_is_rejected() {
    let mut value = read(QUARK);
    let spinor = |h: i64| slot(names::BIS, Atom::num(4), hedge(h));
    let edge = value
        .edges
        .iter()
        .find(|e| e.line == LineKind::Fermion)
        .unwrap()
        .clone();
    // γ5 between the propagator and its source half-edge: (S γ5)_{t s}.
    let s = edge.source_hedge.unwrap() as i64;
    let factor = edge
        .factor
        .replace(spinor(s).to_pattern())
        .with(spinor(900).to_pattern());
    let factor = factor * function!(idenso::dirac::AGS.gamma5, spinor(900), spinor(s));
    value
        .edges
        .iter_mut()
        .find(|e| e.id == edge.id)
        .unwrap()
        .factor = factor;
    assert!(
        matches!(value.validate(), Err(FactorizedError::Gamma5 { .. })),
        "{:?}",
        value.validate()
    );
}

#[test]
fn four_dimensional_slots_and_signs_are_rejected() {
    // mink(4, ..) left in a factor.
    let four = GLUON.replace(
        &mink(&read(GLUON), 2).to_canonical_string(),
        &slot(names::MINK, Atom::num(4), hedge(2)).to_canonical_string(),
    );
    assert_ne!(four, GLUON);
    assert!(matches!(
        FactorizedNumerator::from_text(&four),
        Err(FactorizedError::WrongDimension { .. })
    ));
    // The closed-loop sign is the prefactor's, exactly once.
    let mut value = read(QUARK);
    value.prefactor.loop_sign = 1;
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::LoopSign { .. })
    ));
    let mut value = read(QUARK);
    value.prefactor.nf_power = 2;
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::NfPower { .. })
    ));
}

#[test]
fn the_reserved_extensions_are_rejected() {
    // Insertions are accepted (insertion_tests.rs); one on a vertex that
    // keeps its own factor is refused (the factor is stated once).
    let mut value = read(GLUON);
    let hedges = value
        .local_hedges(Location::Vertex(0))
        .into_iter()
        .collect();
    value.insertions.push(Insertion {
        vertex: 0,
        hedges,
        cycles: (0, 0),
        coefficient: Atom::num(1),
        factor: Atom::num(1),
    });
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::Insertion { vertex: 0, .. })
    ));
    let mut value = read(GLUON);
    value.kind = FactorizedKind::Closure(ClosureSpec {
        edge: 2,
        mass_squared: Atom::num(1),
        auxiliary_mass_squared: Atom::num(1),
        power: 1,
        independent_externals: false,
        taylor: None,
        softs: Vec::new(),
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    });
    // the closure kind is accepted now: this one is refused by its own rules
    // (the massless edge 2 is not a closing line)
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::Closure { .. })
    ));
    let closure = GLUON.replace("kind\tdirect", "kind\tclosure");
    assert!(matches!(
        FactorizedNumerator::parse_text(&closure),
        Err(FactorizedError::Reserved { .. })
    ));
    let insertions = GLUON.replace("insertions\t0", "insertions\t1");
    assert!(matches!(
        FactorizedNumerator::parse_text(&insertions),
        Err(FactorizedError::Reserved { .. })
    ));
    // a v3 record may carry insertions: the declared one is missing here
    let mut marked = read(GLUON);
    marked.record = Some(DirectRecord {
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    });
    let v3 = marked.to_text().unwrap();
    assert!(v3.starts_with(&format!("schema\t{FACTORIZED_SCHEMA_V3}\n")));
    let insertions_v3 = v3.replace("insertions\t0", "insertions\t1");
    assert_ne!(insertions_v3, v3);
    assert!(matches!(
        FactorizedNumerator::parse_text(&insertions_v3),
        Err(FactorizedError::TextSectionCount { .. })
    ));
}

#[test]
fn text_assertions_reject_damaged_records() {
    // An unknown version.
    let v2 = GLUON.replace(FACTORIZED_SCHEMA, "pint6.factorized-numerator.v2");
    assert!(matches!(
        FactorizedNumerator::parse_text(&v2),
        Err(FactorizedError::TextSchema { .. })
    ));
    let v1_boundary = GLUON.replace(FACTORIZED_SCHEMA, "pint6.remainder-boundary.v1");
    assert!(matches!(
        FactorizedNumerator::parse_text(&v1_boundary),
        Err(FactorizedError::TextSchema { .. })
    ));
    // A wrong line-count assertion.
    let lines = GLUON
        .lines()
        .find(|l| l.starts_with("lines\t"))
        .unwrap()
        .to_string();
    let n: usize = lines["lines\t".len()..].parse().unwrap();
    let wrong = GLUON.replace(&lines, &format!("lines\t{}", n + 1));
    assert!(matches!(
        FactorizedNumerator::parse_text(&wrong),
        Err(FactorizedError::TextLineCount { .. })
    ));
    // A section count that does not match.
    let edges = GLUON
        .lines()
        .find(|l| l.starts_with("edges\t"))
        .unwrap()
        .to_string();
    let n: usize = edges["edges\t".len()..].parse().unwrap();
    let wrong = GLUON.replace(&edges, &format!("edges\t{}", n + 1));
    assert!(matches!(
        FactorizedNumerator::parse_text(&wrong),
        Err(FactorizedError::TextSectionCount { .. })
    ));
    // A dropped vertex line: caught by the section count.
    let vertex = GLUON
        .lines()
        .find(|l| l.starts_with("vertex\t"))
        .unwrap()
        .to_string();
    let dropped = GLUON.replace(&format!("{vertex}\n"), "");
    assert!(matches!(
        FactorizedNumerator::parse_text(&dropped),
        Err(FactorizedError::TextSectionCount { .. })
    ));
    // An unknown key.
    let unknown = GLUON.replace("insertions\t0\n", "insertions\t0\nbogus\t1\n");
    assert!(matches!(
        FactorizedNumerator::parse_text(&unknown),
        Err(FactorizedError::TextSyntax { .. })
    ));
    // Text after `end`, and a truncated record.
    assert!(FactorizedNumerator::parse_text(&format!("{GLUON}x\n")).is_err());
    let truncated = GLUON.lines().take(10).collect::<Vec<_>>().join("\n");
    assert!(matches!(
        FactorizedNumerator::parse_text(&truncated),
        Err(FactorizedError::TextSectionCount { .. })
    ));
    let truncated = GLUON.lines().take(30).collect::<Vec<_>>().join("\n");
    assert!(matches!(
        FactorizedNumerator::parse_text(&truncated),
        Err(FactorizedError::TextSyntax { .. })
    ));
}

#[test]
fn from_parts_checks_what_it_is_given() {
    let value = read(QUARK);
    let parts = FactorizedParts {
        label: value.label.clone(),
        dimension: value.dimension,
        loop_edges: value.loop_edges.clone(),
        external_edges: value.external_edges.clone(),
        edges: value.edges.iter().rev().cloned().collect(),
        vertices: value.vertices.iter().rev().cloned().collect(),
        ports: value.ports.clone(),
        projector: value.projector.clone(),
        prefactor: value.prefactor.clone(),
    };
    let rebuilt = FactorizedNumerator::from_parts(parts.clone()).unwrap();
    assert_eq!(rebuilt.edges, value.edges);
    assert_eq!(rebuilt.vertices, value.vertices);
    let mut wrong = parts;
    wrong.ports.clear();
    assert!(matches!(
        FactorizedNumerator::from_parts(wrong),
        Err(FactorizedError::Ports { .. })
    ));
}

/// The refpath convention of a value's prefactor (lane refpath, REFPATH.md
/// section 5): couplings at `gs = 1` (`gs_minus = -1`, `gs_i = i`, `gs2_i = i`,
/// `gs2_minus_i = -i`), colour invariants in refpath's symbols, `refpath::nf`.
#[cfg(feature = "feynkit-ingress")]
fn refpath_prefactor(value: &FactorizedNumerator) -> Atom {
    use idenso::color::CS;
    let i = Atom::i();
    let mut total = value.prefactor.total();
    for (name, image) in [
        ("UFO::gs_minus", Atom::num(-1)),
        ("UFO::gs_i", i.clone()),
        ("UFO::gs2_i", i.clone()),
        ("UFO::gs2_minus_i", -i.clone()),
    ] {
        total = total
            .replace(Atom::var(symbolica::symbol!(name)))
            .with(image);
    }
    let total = total
        .replace(Atom::var(symbolica::symbol!(names::NF)))
        .with(Atom::var(crate::refpath::SYMBOLS.nf));
    crate::refpath::map_invariants(
        &total,
        &Atom::var(symbolica::symbol!("spenso::dA")),
        &Atom::var(CS.nc),
    )
    .expand()
    .together()
}

/// The adapter (`evaluate_factorized`: idenso contraction, vakint's existing
/// p evaluation) on the two fixtures with lane refpath's metric normalisation
/// `1/(i q² (D-1))` gives exactly lane refpath's archived L1 values of the same
/// graphs (B/num-refpath/reference/L1-metric, graphs 000000 and 000002, two
/// Laurent terms, G-scheme).
#[cfg(feature = "feynkit-ingress")]
#[test]
fn the_adapter_reproduces_lane_refpaths_archived_l1_values() {
    let vakint = crate::Vakint::new().unwrap();
    let settings = crate::refpath::vakint_settings(2);
    let gluon = "(-10𝑖/3)*refpath::{}::CA*vakint::{}::ε^(-1)+(10𝑖/3)*log(-1*vakint::{symmetric,linear}::dot(vakint::{spenso::tensor,spenso::rank1}::p(1),vakint::{spenso::tensor,spenso::rank1}::p(1))*vakint::{}::mursq^(-1))*refpath::{}::CA+(𝑖/9)*refpath::{}::CA";
    let quark = "(-4𝑖/3)*log(-1*vakint::{symmetric,linear}::dot(vakint::{spenso::tensor,spenso::rank1}::p(1),vakint::{spenso::tensor,spenso::rank1}::p(1))*vakint::{}::mursq^(-1))*refpath::{}::TF*refpath::{}::nf+(-4𝑖/9)*refpath::{}::TF*refpath::{}::nf+(4𝑖/3)*refpath::{}::TF*refpath::{}::nf*vakint::{}::ε^(-1)";
    for (text, expected) in [(GLUON, gluon), (QUARK, quark)] {
        let mut value = read(text);
        value.projector.as_mut().unwrap().p2_power = -1;
        let d = Atom::var(value.dimension);
        let prefactor = refpath_prefactor(&value) / (Atom::i() * (&d - Atom::num(1)));
        let got = evaluate_factorized(&vakint, &settings, &value, &prefactor).unwrap();
        let expected = crate::utils::vakint_macros::vk_parse!(expected).unwrap();
        assert!(
            (&got.value - &expected).expand().is_zero(),
            "{}: {} vs {}",
            value.label,
            got.value,
            expected
        );
    }
    // The closure kind stays reserved for the adapter too.
    let mut value = read(GLUON);
    value.kind = FactorizedKind::Closure(ClosureSpec {
        edge: 2,
        mass_squared: Atom::num(1),
        auxiliary_mass_squared: Atom::num(1),
        power: 1,
        independent_externals: false,
        taylor: None,
        softs: Vec::new(),
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    });
    assert!(matches!(
        evaluate_factorized(&vakint, &settings, &value, &Atom::num(1)),
        Err(FactorizedError::Reserved { .. })
    ));
}

/// The flavour symbol is the caller's: in the beta thread's namespace the
/// quark bubble's weight has coefficient -1 on `beta::nf` (loop sign -1,
/// n_f^1) and nothing is left on `vakint::nf`; with the default total the
/// opposite (review O3/interface.md section 2).
#[test]
fn the_flavour_symbol_is_the_callers() {
    let quark = read(QUARK);
    let beta_nf = Atom::var(symbolica::symbol!("beta::nf"));
    let vakint_nf = Atom::var(symbolica::symbol!(names::NF));
    let weight = |total: Atom| {
        (&total / &quark.prefactor.scalars / &quark.prefactor.colour / &quark.prefactor.overall)
            .expand()
    };
    let in_beta = weight(quark.prefactor.total_with_flavour(&beta_nf));
    assert_eq!(in_beta.coefficient(beta_nf.as_view()), Atom::num(-1));
    assert_eq!(in_beta.coefficient(vakint_nf.as_view()), Atom::num(0));
    let in_vakint = weight(quark.prefactor.total());
    assert_eq!(in_vakint.coefficient(vakint_nf.as_view()), Atom::num(-1));
    assert_eq!(in_vakint.coefficient(beta_nf.as_view()), Atom::num(0));
}

/// `Σ_t coefficient_t Π_v w_{t,v}(k_v)` at one channel.
fn merged_value(merged: &[MergedColourTerm], key: &[usize]) -> Atom {
    let mut sum = Atom::Zero;
    for term in merged {
        let mut product = term.coefficient.clone();
        for (vector, row) in term.weights.iter().zip(key) {
            product *= Atom::num(vector[*row].clone());
        }
        sum += product;
    }
    sum.expand()
}

fn all_keys(dims: &[usize]) -> Vec<Vec<usize>> {
    let mut keys = vec![Vec::new()];
    for &rows in dims {
        keys = keys
            .into_iter()
            .flat_map(|key| {
                (0..rows).map(move |row| {
                    let mut next = key.clone();
                    next.push(row);
                    next
                })
            })
            .collect();
    }
    keys
}

/// Lane share (A1): the merged colour terms reproduce every channel exactly
/// (including the absent, i.e. zero, ones), numeric and symbolic, and give
/// the expected number of terms on tensors of known rank.
#[test]
fn merged_colour_terms_reconstruct_every_channel() {
    let ca = Atom::var(symbolica::symbol!("spenso::share_test_CA"));
    let cf = Atom::var(symbolica::symbol!("spenso::share_test_CF"));
    let n = |v: i64| Atom::num(v);
    let q = |a: i64, b: i64| Atom::num(symbolica::domains::rational::Rational::from((a, b)));
    // (dims, channels (choices, coefficient), expected number of merged terms)
    type Case = (Vec<usize>, Vec<(Vec<usize>, Atom)>, usize);
    let cases: Vec<Case> = vec![
        // A rank-one 3x3 tensor (a_i b_j): one term.
        (
            vec![3, 3],
            all_keys(&[3, 3])
                .into_iter()
                .map(|k| {
                    let (a, b) = ([1, -2, 3][k[0]], [2, 1, -1][k[1]]);
                    (k, n(a * b))
                })
                .collect(),
            1,
        ),
        // One multi-row local (three rows): any vector is one term.
        (
            vec![3],
            vec![(vec![0], n(1)), (vec![1], n(1)), (vec![2], n(2))],
            1,
        ),
        // Symbolic: c(k) = CA^2 x_k + CA CF y_k, x and y independent: two terms.
        (
            vec![3],
            vec![
                (vec![0], &ca * &ca * n(2) + &ca * &cf),
                (vec![1], &ca * &ca * q(1, 2)),
                (vec![2], &ca * &cf * n(-3)),
            ],
            2,
        ),
        // Symbolic, proportional channels: one term.
        (
            vec![3, 3],
            all_keys(&[3, 3])
                .into_iter()
                .filter(|k| k != &vec![1, 1])
                .map(|k| (k.clone(), (&ca - &cf) * n(1 + k[0] as i64)))
                .collect(),
            0, // checked by reconstruction only
        ),
    ];
    for (dims, channels, expected) in cases {
        let merged = merge_colour_channels(&dims, &channels).expect("mergeable");
        if expected > 0 {
            assert_eq!(merged.len(), expected, "{dims:?}");
        }
        for term in &merged {
            for vector in &term.weights {
                let first = vector.iter().find(|w| !w.is_zero()).unwrap();
                assert!(first.is_one(), "weights normalised to a leading 1");
            }
        }
        let table = channels
            .iter()
            .cloned()
            .collect::<std::collections::BTreeMap<_, _>>();
        for key in all_keys(&dims) {
            let expected = table.get(&key).cloned().unwrap_or(Atom::Zero).expand();
            let found = merged_value(&merged, &key);
            assert!(
                (&expected - &found).expand().is_zero(),
                "{dims:?} {key:?}: {expected} vs {found}"
            );
        }
    }
    // Imaginary channel coefficients (a quark loop with three gluons) merge
    // with real weights: the imaginary unit is part of the colour monomial.
    let i = Atom::i();
    let imaginary = vec![
        (vec![0], &i * q(-3, 4)),
        (vec![1], &i * q(3, 8)),
        (vec![2], &i * q(3, 8)),
    ];
    let merged = merge_colour_channels(&[3], &imaginary).expect("mergeable");
    assert_eq!(merged.len(), 1);
    for (key, expected) in &imaginary {
        assert!((expected - merged_value(&merged, key)).expand().is_zero());
    }
    // A non-rational number refuses the merge.
    assert!(merge_colour_channels(&[2], &[(vec![0], Atom::num(1.5)), (vec![1], n(1))]).is_none());
}
