//! Offline bridge compiled separately against the old and current Symbolica.
//! Only trusted generated RustRed assets are supported; this is not a file parser
//! for untrusted input. RustRed's structural sections never pass through algebra.
use serde::{Deserialize, Serialize};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
    sync::Arc,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::{Coefficient, CoefficientView},
    domains::{
        integer::{Integer, IntegerRing, Z},
        rational_polynomial::RationalPolynomial,
    },
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
    state::State,
};

type Poly = MultivariatePolynomial<IntegerRing, u16>;
type Terms = Vec<(String, Vec<u16>)>;

// Preserve sparse terms and the complete ordered variable map, including unused
// variables and zero coefficients. Printing only the rational function would
// lose that map. General catalog atoms use Symbolica's exact canonical syntax.
#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
enum Record {
    Rational {
        variables: Vec<String>,
        numerator: Terms,
        denominator: Terms,
    },
    Expression {
        canonical: String,
    },
}

impl Record {
    fn from_atom(atom: &Atom) -> Self {
        if let AtomView::Num(n) = atom.as_view() {
            if let CoefficientView::RationalPolynomial(p) = n.get_coeff_view() {
                let p = p.deserialize();
                let variables = p
                    .numerator
                    .variables()
                    .iter()
                    .map(|v| match v {
                        PolyVariable::Symbol(s) => s.to_atom().to_canonical_string(),
                        _ => panic!("unsupported non-symbol polynomial variable"),
                    })
                    .collect();
                assert_eq!(p.numerator.variables(), p.denominator.variables());
                let terms = |poly: &Poly| {
                    poly.into_iter()
                        .map(|m| (m.coefficient.to_string(), m.exponents.to_vec()))
                        .collect()
                };
                return Self::Rational {
                    variables,
                    numerator: terms(&p.numerator),
                    denominator: terms(&p.denominator),
                };
            }
        }
        Self::Expression {
            canonical: atom.to_canonical_string(),
        }
    }

    fn to_atom(&self) -> Atom {
        match self {
            Self::Rational {
                variables,
                numerator,
                denominator,
            } => {
                let variables = Arc::new(
                    variables
                        .iter()
                        .map(|name| {
                            let atom = Atom::parse(name, "migration", Default::default()).unwrap();
                            let AtomView::Var(var) = atom.as_view() else {
                                panic!("not a symbol")
                            };
                            PolyVariable::Symbol(var.get_symbol())
                        })
                        .collect::<Vec<_>>(),
                );
                let poly = |terms: &Terms| {
                    let mut poly = Poly::new(&Z, Some(terms.len()), variables.clone());
                    for (coefficient, exponents) in terms {
                        poly.append_monomial_back(
                            coefficient.parse::<Integer>().unwrap(),
                            exponents,
                        );
                    }
                    poly
                };
                let p = RationalPolynomial {
                    numerator: poly(numerator),
                    denominator: poly(denominator),
                };
                let mut atom = Atom::new();
                atom.to_num(Coefficient::RationalPolynomial(p));
                assert_eq!(
                    &Self::from_atom(&atom),
                    self,
                    "polynomial transport must preserve every term and variable"
                );
                atom
            }
            Self::Expression { canonical } => {
                Atom::parse(canonical, "migration", Default::default()).unwrap()
            }
        }
    }
}

fn read_file(path: &Path) -> Vec<u8> {
    let file = File::open(path).unwrap();
    let mut reader: Box<dyn Read> = if path.extension().is_some_and(|x| x == "gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut bytes = vec![];
    reader.read_to_end(&mut bytes).unwrap();
    bytes
}

fn length(source: &mut &[u8]) -> usize {
    let n = u64::from_le_bytes(source[..8].try_into().unwrap()) as usize;
    *source = &source[8..];
    n
}

fn sections(bytes: &[u8]) -> Vec<(u16, &[u8])> {
    assert_eq!(&bytes[..8], b"RRPBIN\r\n");
    assert_eq!(u32::from_le_bytes(bytes[8..12].try_into().unwrap()), 1);
    let n = u32::from_le_bytes(bytes[16..20].try_into().unwrap());
    let mut source = &bytes[20..];
    let mut sections = vec![];
    for _ in 0..n {
        let tag = u16::from_le_bytes(source[..2].try_into().unwrap());
        assert_eq!(&source[2..4], &[0, 0]);
        source = &source[4..];
        let size = length(&mut source);
        sections.push((tag, &source[..size]));
        source = &source[size..];
    }
    assert!(source.is_empty());
    sections
}

fn main() {
    let _noise = Atom::parse(
        "native_asset_migration_noise::f(a,b)",
        "native_asset_migration_noise",
        Default::default(),
    )
    .unwrap();
    let args: Vec<_> = std::env::args().collect();
    let bytes = read_file(Path::new(&args[2]));
    let sections = sections(&bytes);
    let &(atom_tag, table) = sections
        .iter()
        .find(|(tag, _)| *tag == 2 || *tag == 7)
        .unwrap();
    if args[1] == "decode" {
        let mut state = sections.iter().find(|(tag, _)| *tag == 1).unwrap().1;
        let map = State::import(&mut state, None).unwrap();
        assert!(state.is_empty());
        let mut table = table;
        let count = length(&mut table);
        let mut output = BufWriter::new(flate2::write::GzEncoder::new(
            File::create(&args[3]).unwrap(),
            flate2::Compression::fast(),
        ));
        for _ in 0..count {
            let size = length(&mut table);
            let mut frame = &table[..size];
            let atom = Atom::import_with_map(&mut frame, &map).unwrap();
            assert!(frame.is_empty());
            let record = Record::from_atom(&atom);
            assert_eq!(
                record.to_atom(),
                atom,
                "native decoder text/term round trip failed"
            );
            serde_json::to_writer(&mut output, &record).unwrap();
            output.write_all(b"\n").unwrap();
            table = &table[size..];
        }
        assert!(table.is_empty());
        output.flush().unwrap();
        output.into_inner().unwrap().finish().unwrap();
        eprintln!("decoded {count} records from {}", args[2]);
    } else if args[1] == "encode" {
        let mut old_table = table;
        let count = length(&mut old_table);
        let input = BufReader::new(flate2::read::GzDecoder::new(File::open(&args[3]).unwrap()));
        let mut table = Vec::new();
        table.extend_from_slice(&(count as u64).to_le_bytes());
        let mut symbols = Atom::new().get_all_symbols(true);
        let mut n = 0;
        for line in input.lines() {
            let record: Record = serde_json::from_str(&line.unwrap()).unwrap();
            let atom = record.to_atom();
            symbols.extend(atom.get_all_symbols(true));
            let mut frame = vec![];
            atom.as_view().write(&mut frame).unwrap();
            table.extend_from_slice(&(frame.len() as u64).to_le_bytes());
            table.extend_from_slice(&frame);
            n += 1;
        }
        assert_eq!(n, count);
        let mut state = vec![];
        State::export_partial(&mut state, symbols).unwrap();
        let mut output = bytes[..20].to_vec();
        for (tag, original) in sections {
            let payload = if tag == 1 {
                &state[..]
            } else if tag == atom_tag {
                &table[..]
            } else {
                original
            };
            output.extend_from_slice(&tag.to_le_bytes());
            output.extend_from_slice(&0u16.to_le_bytes());
            output.extend_from_slice(&(payload.len() as u64).to_le_bytes());
            output.extend_from_slice(payload);
        }
        std::fs::write(&args[4], output).unwrap();
        eprintln!("encoded {count} records to {}", args[4]);
    } else {
        panic!("expected decode or encode");
    }
}
