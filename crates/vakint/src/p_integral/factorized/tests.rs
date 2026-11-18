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
    let mut value = read(GLUON);
    value.insertions.push(Insertion {
        vertex: 0,
        coefficient: Atom::num(1),
        factor: Atom::num(1),
    });
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::Reserved { .. })
    ));
    assert!(matches!(
        value.to_text(),
        Err(FactorizedError::Reserved { .. })
    ));
    let mut value = read(GLUON);
    value.kind = FactorizedKind::Closure(ClosureSpec {
        edge: 2,
        mass_squared: Atom::num(1),
        power: 1,
        taylor: None,
    });
    assert!(matches!(
        value.validate(),
        Err(FactorizedError::Reserved { .. })
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
