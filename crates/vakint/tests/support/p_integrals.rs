//! Shared helpers of the p-integral (massless propagator-type) tests.
//!
//! Include it with `#[path = "support/p_integrals.rs"] mod p_integrals;`.
//! Integrals are written as strings in vakint's `topo(prop(...))` notation,
//! with an optional numerator factor in front.

use std::sync::LazyLock;

use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::id::Replacement;
use vakint::diagram_integral::DiagramIntegralOptions;
use vakint::symbols::S;
use vakint::{Vakint, VakintExpression, VakintSettings, p_integral, vakint_parse};

static VAKINT: LazyLock<Vakint> = LazyLock::new(|| Vakint::new().unwrap());

/// The scalar `phi` of the DOT tests: massless, vertices of valence 3 to 10.
pub static PHI: LazyLock<Model> = LazyLock::new(|| {
    Model::from_json(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/phi_p_graphs.json"
    )))
    .unwrap()
});

pub fn parse(input: &str) -> Atom {
    vakint_parse!(input).unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// Vakint's canonical long form (family, matched powers, `σ·p(i)`) of `input`.
pub fn canonical(input: &Atom) -> Atom {
    VAKINT
        .to_canonical(&VakintSettings::default(), input.as_view(), false)
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// The canonical short form `topo(NAME(σ·p(i), a_1, ..., a_N))` of `input`.
pub fn short_form(input: &str) -> Atom {
    VAKINT
        .to_canonical(&VakintSettings::default(), parse(input).as_view(), true)
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// `input` is scaleless: it canonicalizes and evaluates to an exact zero.
pub fn assert_p_scaleless(input: &str) {
    let atom = parse(input);
    assert!(canonical(&atom).is_zero(), "{input} must vanish");
    assert!(
        VAKINT
            .evaluate(&VakintSettings::default(), atom.as_view())
            .unwrap_or_else(|error| panic!("{input}: {error}"))
            .is_zero(),
        "{input} must evaluate to 0"
    );
}

/// `input` is rejected with an error mentioning `fragment`.
pub fn assert_p_error(input: &str, fragment: &str) {
    match VAKINT.to_canonical(&VakintSettings::default(), parse(input).as_view(), false) {
        Ok(result) => panic!("{input} was accepted as {result}"),
        Err(error) => assert!(
            error.to_string().contains(fragment),
            "{input}: '{error}' does not mention '{fragment}'"
        ),
    }
}

/// The short form of `input` is a sector of `family` with `pinched` lines.
pub fn assert_p_family(input: &str, family: &str, pinched: usize) {
    let short = short_form(input);
    let factors = match short.as_view() {
        AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
        view => vec![view],
    };
    let head = factors
        .into_iter()
        .find_map(|factor| match factor {
            AtomView::Fun(topo) if topo.get_symbol() == S.topo => topo.iter().next(),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{input}: {short} has no topo(...)"));
    let AtomView::Fun(head) = head else {
        panic!("{input}: {short} is not a short form");
    };
    assert_eq!(
        head.get_symbol().get_stripped_name(),
        family,
        "{input}: {short}"
    );
    let zeros = head.iter().skip(1).filter(|power| power.is_zero()).count();
    assert_eq!(zeros, pinched, "{input}: {short}");
}

/// The top of `family` with line powers `1, 2, ..., N` (no graph symmetry
/// survives distinct powers) in the family routing, times `numerator`.
pub fn distinct_power_top(family: &str, lines: usize, numerator: &str) -> Atom {
    let powers = (1..=lines)
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",");
    parse(numerator) * canonical(&parse(&format!("topo({family}(p(1),{powers}))")))
}

/// Two inputs have the same canonical form.
pub fn assert_p_same_canonical(left: &str, right: &str) {
    assert_eq!(
        canonical(&parse(left)),
        canonical(&parse(right)),
        "{left} vs {right}"
    );
}

/// A deterministic xorshift generator for reproducible random routings.
struct Random(u64);

impl Random {
    fn next(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }

    fn sign(&mut self) -> i64 {
        if self.next(2) == 0 { 1 } else { -1 }
    }
}

/// Number of distinct loop momenta `k(j)` of `input`, which uses `k(1..=L)`.
fn loops(input: &Atom) -> usize {
    input
        .pattern_match(&parse("k(i_)").to_pattern(), None, None)
        .map(|matched| matched[&wildcard("i_")].clone())
        .collect::<std::collections::BTreeSet<_>>()
        .len()
}

fn wildcard(name: &str) -> symbolica::atom::Symbol {
    parse(name).get_symbol().unwrap()
}

/// Rewrite every `prop(...)` of `input` with `rewrite(id, from, to, momentum)`.
fn map_props(input: &Atom, rewrite: impl Fn(i64, i64, i64, Atom) -> (i64, i64, Atom)) -> Atom {
    input.replace_map(|view, _, out| {
        let AtomView::Fun(prop) = view else {
            return;
        };
        if prop.get_symbol() != S.prop || prop.get_nargs() != 5 {
            return;
        }
        let arguments = prop.iter().collect::<Vec<_>>();
        let integer = |view: AtomView| i64::try_from(view).unwrap();
        let AtomView::Fun(edge) = arguments[1] else {
            panic!("{view} has no edge(a,b)");
        };
        let mut ends = edge.iter().map(integer);
        let id = integer(arguments[0]);
        let (from, to, momentum) = rewrite(
            id,
            ends.next().unwrap(),
            ends.next().unwrap(),
            arguments[2].to_owned(),
        );
        **out = S.prop.call((
            id,
            S.edge.call((from, to)),
            momentum.expand(),
            arguments[3].to_owned(),
            arguments[4].to_owned(),
        ));
    })
}

/// The whole integrand (numerator and props) with `rules` applied at once.
fn substitute(input: &Atom, rules: &[(Atom, Atom)]) -> Atom {
    let rules = rules
        .iter()
        .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
        .collect::<Vec<_>>();
    map_props(&input.replace_multiple(&rules), |_, a, b, q| (a, b, q))
}

/// Routing invariance (section 7): for every seed, a random unimodular change
/// of loop momenta `k -> M k + c·p(1)`, random edge flips and a random vertex
/// relabeling, then the reversal `p(1) -> -p(1)` and the renaming
/// `p(1) -> p(3)` of the whole integrand, all give exactly the canonical form
/// of `input` (with the reversal and the renaming undone). Choose powers that
/// break the graph symmetries: a vertex relabeling is only exact up to the
/// automorphisms that preserve the powers.
pub fn assert_p_routing_invariant(input: &str, seeds: std::ops::Range<u64>) {
    routing_invariance(&parse(input), seeds, true);
}

/// [`assert_p_routing_invariant`] for a parsed integral.
pub fn assert_p_atom_routing_invariant(original: &Atom, seeds: std::ops::Range<u64>) {
    routing_invariance(original, seeds, true);
}

/// [`assert_p_routing_invariant`] without vertex relabelings, for graphs with
/// two vertices: exchanging them maps every line to itself, so the reversal
/// `q -> -q` is a symmetry and a relabeling may pick either orientation.
pub fn assert_p_reroutings_invariant(input: &str, seeds: std::ops::Range<u64>) {
    routing_invariance(&parse(input), seeds, false);
}

fn routing_invariance(original: &Atom, seeds: std::ops::Range<u64>, relabel: bool) {
    let original = original.clone();
    let expected = canonical(&original);
    let loops = loops(&original);
    let (k, p, p3) = (
        |axis: usize| vakint_parse!(format!("k({axis})")).unwrap(),
        parse("p(1)"),
        parse("p(3)"),
    );
    for seed in seeds {
        let mut random =
            Random(0x9E37_79B9_7F4A_7C15 ^ (seed + 1).wrapping_mul(0x2545_F491_4F6C_DD1D));
        // A unimodular integer matrix from random elementary operations.
        let mut matrix = (0..loops)
            .map(|i| (0..loops).map(|j| i64::from(i == j)).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        for _ in 0..3 * loops {
            let (i, j) = (random.next(loops), random.next(loops));
            if i != j {
                let sign = random.sign();
                let row = matrix[j].clone();
                matrix[i]
                    .iter_mut()
                    .zip(row)
                    .for_each(|(a, b)| *a += sign * b);
            } else {
                matrix[i].iter_mut().for_each(|a| *a = -*a);
            }
        }
        let images = (0..loops)
            .map(|i| {
                let shift = random.next(3) as i64 - 1;
                let image = (0..loops).fold(Atom::num(shift) * &p, |sum, j| {
                    sum + Atom::num(matrix[i][j]) * k(j + 1)
                });
                (k(i + 1), image)
            })
            .collect::<Vec<_>>();
        let rerouted = substitute(&original, &images);
        let flips = (0..64).map(|_| random.next(2) == 1).collect::<Vec<_>>();
        let shift = 10 + random.next(5) as i64;
        let mut order = (0..32).collect::<Vec<i64>>();
        for i in (1..order.len()).rev() {
            order.swap(i, random.next(i + 1));
        }
        let relabeled = map_props(&rerouted, |id, a, b, q| {
            let (a, b) = if relabel {
                (order[a as usize] + shift, order[b as usize] + shift)
            } else {
                (a, b)
            };
            if flips[id as usize % flips.len()] {
                (b, a, -q)
            } else {
                (a, b, q)
            }
        });
        let found = canonical(&relabeled);
        assert_eq!(found, expected, "seed {seed}: {relabeled}");
        let reversed = substitute(&relabeled, &[(p.clone(), -p.clone())]);
        let found = substitute(&canonical(&reversed), &[(p.clone(), -p.clone())]);
        assert_eq!(found, expected, "seed {seed}, p(1) -> -p(1): {reversed}");
        let renamed = substitute(&relabeled, &[(p.clone(), p3.clone())]);
        let found = substitute(&canonical(&renamed), &[(p3.clone(), p.clone())]);
        assert_eq!(found, expected, "seed {seed}, p(1) -> p(3): {renamed}");
    }
}

/// Read a compact DOT p-graph back into vakint (`from_diagram`).
pub fn from_dot(dot: &str) -> Atom {
    let diagram =
        FeynmanDiagram::from_dot(PHI.clone(), dot).unwrap_or_else(|error| panic!("{dot}: {error}"));
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    Atom::from(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions::default(),
        )
        .unwrap_or_else(|error| panic!("{dot}: {error}")),
    )
}

/// DOT round trip (section 3): `atom -> to_dot -> from_diagram` keeps the
/// canonical form, and `dot -> atom -> dot -> atom` is a fixed point after
/// the first export.
pub fn assert_p_dot_round_trip(input: &str) {
    let atom = parse(input);
    let dot = p_integral::to_dot(atom.as_view(), "phi").unwrap();
    let imported = from_dot(&dot);
    assert_eq!(canonical(&imported), canonical(&atom), "{input}\n{dot}");
    let again = p_integral::to_dot(imported.as_view(), "phi").unwrap();
    assert_eq!(from_dot(&again), imported, "{dot}\n{again}");
}
