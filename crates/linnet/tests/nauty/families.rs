//! Families from the Traces benchmarks that are hard for the search rather than for refinement:
//! incidence graphs of the Desarguesian projective planes PG(2,q), whose many equal-trace leaves
//! need a target cell that splits as many cells as possible, and Latin-square graphs, which are
//! strongly regular with (for random squares) trivial groups, so every leaf has to be compared.
//! Group orders are known in closed form or recorded from nauty.

use rand::{rngs::SmallRng, SeedableRng};

use super::{graph::Graph, runs::full_size};

fn rng() -> SmallRng {
    SmallRng::seed_from_u64(2026)
}

/// The field GF(p^k) on `0..p^k`, an element being the base-`p` digits of its polynomial.
struct Field {
    q: usize,
    add: Vec<usize>,
    mul: Vec<usize>,
}

impl Field {
    /// `modulus` holds the `k` low coefficients of a monic irreducible polynomial of degree `k`
    /// over GF(p); `[0]` gives the prime field.
    fn new(p: usize, modulus: &[usize]) -> Field {
        let k = modulus.len();
        let q = p.pow(k as u32);
        let digits = |a: usize| (0..k).map(move |i| a / p.pow(i as u32) % p);
        let pack = |d: &[usize]| d.iter().rev().fold(0, |a, &x| a * p + x);
        let mut add = vec![0; q * q];
        let mut mul = vec![0; q * q];
        for a in 0..q {
            for b in 0..q {
                let sum: Vec<usize> = digits(a).zip(digits(b)).map(|(x, y)| (x + y) % p).collect();
                add[a * q + b] = pack(&sum);
                let mut product = vec![0; 2 * k - 1];
                for (i, x) in digits(a).enumerate() {
                    for (j, y) in digits(b).enumerate() {
                        product[i + j] = (product[i + j] + x * y) % p;
                    }
                }
                // x^k = -(modulus): reduce from the top degree down
                for d in (k..2 * k - 1).rev() {
                    let t = std::mem::take(&mut product[d]);
                    for (i, &c) in modulus.iter().enumerate() {
                        product[d - k + i] = (product[d - k + i] + (p - c) * t) % p;
                    }
                }
                mul[a * q + b] = pack(&product[..k]);
            }
        }
        Field { q, add, mul }
    }

    fn dot(&self, a: [usize; 3], b: [usize; 3]) -> usize {
        a.iter().zip(b).fold(0, |s, (&x, y)| {
            self.add[s * self.q + self.mul[x * self.q + y]]
        })
    }
}

/// The point-line incidence graph of PG(2,q): points first, then lines, both as normalised
/// homogeneous coordinates; a line is incident to the points orthogonal to it.
fn projective_plane(field: &Field) -> Graph {
    let q = field.q;
    let points: Vec<[usize; 3]> = (0..q)
        .flat_map(|x| (0..q).map(move |y| [x, y, 1]))
        .chain((0..q).map(|x| [x, 1, 0]))
        .chain([[1, 0, 0]])
        .collect();
    let m = points.len();
    let mut edges = vec![];
    for (l, &line) in points.iter().enumerate() {
        for (v, &point) in points.iter().enumerate() {
            if field.dot(line, point) == 0 {
                edges.push((v, m + l));
            }
        }
    }
    assert_eq!(edges.len(), m * (q + 1), "every line has q + 1 points");
    Graph::undirected(2 * m, edges)
}

#[test]
fn projective_plane_incidence_graphs() {
    // (p, modulus) of GF(p^e); larger planes only at full size
    let small: &[(usize, &[usize])] = &[(2, &[0]), (3, &[0]), (2, &[1, 1]), (5, &[0]), (7, &[0])];
    let large: &[(usize, &[usize])] = &[(2, &[1, 1, 0]), (3, &[1, 0]), (11, &[0]), (13, &[0])];
    let mut rng = rng();
    for &(p, modulus) in small.iter().chain(if full_size() { large } else { &[] }) {
        let field = Field::new(p, modulus);
        let (q, e) = (field.q as u128, modulus.len() as u128);
        let g = projective_plane(&field);
        let symmetry = g.symmetry();
        // PGammaL(3, q) and the dualities: 2 e q^3 (q^3 - 1) (q^2 - 1)
        let order = 2 * e * q.pow(3) * (q.pow(3) - 1) * (q * q - 1);
        assert_eq!(symmetry.order, order, "PG(2,{q})");
        assert_eq!(symmetry.node_orbits.len(), 1, "PG(2,{q}) is self-dual");
        g.check_invariance(&symmetry.canonical_form, 2, &mut rng);
    }
}

/// The Latin-square graph: cells adjacent when they share a row, a column or a symbol.
fn latin_square_graph(rows: &[&str]) -> Graph {
    let n = rows.len();
    let symbol = |c: usize| rows[c / n].as_bytes()[c % n];
    let edges = (0..n * n)
        .flat_map(|a| (a + 1..n * n).map(move |b| (a, b)))
        .filter(|&(a, b)| a / n == b / n || a % n == b % n || symbol(a) == symbol(b))
        .collect();
    Graph::undirected(n * n, edges)
}

#[test]
fn latin_square_graphs() {
    // random squares (Jacobson-Matthews), whose graphs have trivial groups per nauty
    let order_8 = [
        [
            "01536274", "72304651", "60172345", "14720563", "47265130", "35417026", "23651407",
            "56043712",
        ],
        [
            "47012356", "62150734", "71534062", "04673125", "56207413", "35461270", "23746501",
            "10325647",
        ],
    ];
    let order_12 = [
        "90381ba72564",
        "27b658034a91",
        "6419a572380b",
        "a853924071b6",
        "0a723489b615",
        "529046ba1378",
        "394701586b2a",
        "136b7095a482",
        "8ba127649053",
        "7584631b02a9",
        "b62a89315740",
        "4105ba268937",
    ];
    let cyclic: Vec<String> = (0..8)
        .map(|r| (0..8).map(|c| char::from(b'0' + (r + c) % 8)).collect())
        .collect();
    let cyclic: Vec<&str> = cyclic.iter().map(String::as_str).collect();
    let mut cases: Vec<(&[&str], u128)> = vec![
        (&order_8[0], 1),
        (&order_8[1], 1),
        // the Cayley table of Z_n: 6 n^2 phi(n)
        (&cyclic, 6 * 64 * 4),
    ];
    if full_size() {
        cases.push((&order_12, 1));
    }
    let mut rng = rng();
    let mut forms = vec![];
    for (rows, order) in cases {
        let g = latin_square_graph(rows);
        let symmetry = g.symmetry();
        assert_eq!(symmetry.order, order, "{rows:?}");
        g.check_invariance(&symmetry.canonical_form, 2, &mut rng);
        forms.push(symmetry.canonical_form);
    }
    // different main classes: nauty finds the two graphs non-isomorphic
    assert_ne!(forms[0], forms[1]);
}
