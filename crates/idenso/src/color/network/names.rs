//! Names of fully contracted networks of symmetric tensors (color.h
//! `makeds`), and the relations color.h applies between them.

use std::{collections::HashMap, sync::LazyLock};

use super::graph::{Net, Obj, Rep};
use super::poly::{C, CA, CR, FIRST_NAMED, I2R, Mono, NA, NR, Poly};

/// A topology: one index list per symmetric tensor, in argument order.
struct Template {
    name: &'static str,
    nodes: &'static [&'static [u16]],
}

macro_rules! t {
    ($name:literal, $($node:expr),+ $(,)?) => {
        Template { name: $name, nodes: &[$(&$node),+] }
    };
}

/// color.h's topologies (two typos in its d5333 and d5533a patterns, which
/// can never match there, are corrected).
const TEMPLATES: &[Template] = &[
    t!("d33", [1, 2, 3], [1, 2, 3]),
    t!("d44", [1, 2, 3, 4], [1, 2, 3, 4]),
    t!("d55", [1, 2, 3, 4, 5], [1, 2, 3, 4, 5]),
    t!("d66", [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6]),
    t!("d77", [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 6, 7]),
    t!("d88", [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5, 6, 7, 8]),
    t!("d433", [1, 2, 3, 4], [1, 2, 5], [3, 4, 5]),
    t!("d444", [1, 2, 3, 4], [1, 2, 5, 6], [3, 4, 5, 6]),
    t!("d633", [1, 2, 3, 4, 5, 6], [1, 2, 3], [4, 5, 6]),
    t!("d3333", [1, 2, 3], [1, 4, 5], [2, 4, 6], [3, 5, 6]),
    t!("d543", [1, 2, 3, 4, 5], [1, 2, 3, 6], [4, 5, 6]),
    t!("d743", [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4], [5, 6, 7]),
    t!("d644", [1, 2, 3, 4, 5, 6], [1, 2, 3, 7], [4, 5, 6, 7]),
    t!("d554", [1, 2, 3, 4, 5], [1, 2, 3, 6, 7], [4, 5, 6, 7]),
    t!("d653", [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 7], [5, 6, 7]),
    t!("d5333", [1, 2, 3, 4, 5], [1, 2, 6], [3, 4, 7], [5, 6, 7]),
    t!("d4433a", [1, 2, 3, 4], [1, 5, 6, 7], [2, 3, 5], [4, 6, 7]),
    t!("d4433b", [1, 2, 3, 4], [1, 2, 5, 6], [3, 5, 7], [4, 6, 7]),
    t!("d4433c", [1, 2, 3, 4], [1, 2, 5, 6], [3, 4, 7], [5, 6, 7]),
    t!("d844", [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4], [5, 6, 7, 8]),
    t!("d853", [1, 2, 3, 4, 5, 6, 7, 8], [1, 2, 3, 4, 5], [6, 7, 8]),
    t!("d763", [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 5, 8], [6, 7, 8]),
    t!("d754", [1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4, 8], [5, 6, 7, 8]),
    t!("d7333", [1, 2, 3, 4, 5, 6, 7], [1, 2, 3], [4, 5, 8], [6, 7, 8]),
    t!("d664", [1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 7, 8], [5, 6, 7, 8]),
    t!("d655", [1, 2, 3, 4, 5, 6], [1, 2, 3, 7, 8], [4, 5, 6, 7, 8]),
    t!("d6433a", [1, 2, 3, 4, 5, 6], [1, 2, 7, 8], [3, 4, 5], [6, 7, 8]),
    t!("d6433b", [1, 2, 3, 4, 5, 6], [1, 2, 7, 8], [3, 4, 7], [5, 6, 8]),
    t!("d6433c", [1, 2, 3, 4, 5, 6], [1, 2, 3, 7], [4, 5, 8], [6, 7, 8]),
    t!("d5533a", [1, 2, 3, 4, 5], [1, 2, 6, 7, 8], [3, 4, 6], [5, 7, 8]),
    t!("d5533b", [1, 2, 3, 4, 5], [1, 2, 6, 7, 8], [3, 4, 5], [6, 7, 8]),
    t!("d5533c", [1, 2, 3, 4, 5], [1, 2, 3, 6, 7], [4, 6, 8], [5, 7, 8]),
    t!("d5533d", [1, 2, 3, 4, 5], [1, 2, 3, 6, 7], [4, 5, 8], [6, 7, 8]),
    t!("d5443a", [1, 2, 3, 4, 5], [1, 2, 6, 7], [3, 6, 7, 8], [4, 5, 8]),
    t!("d5443b", [1, 2, 3, 4, 5], [1, 2, 3, 6], [4, 6, 7, 8], [5, 7, 8]),
    t!("d5443c", [1, 2, 3, 4, 5], [1, 2, 6, 7], [3, 4, 6, 8], [5, 7, 8]),
    t!("d5443d", [1, 2, 3, 4, 5], [1, 2, 3, 6], [4, 5, 7, 8], [6, 7, 8]),
    t!("d4444a", [1, 2, 3, 4], [1, 2, 5, 6], [5, 6, 7, 8], [3, 4, 7, 8]),
    t!("d4444b", [1, 2, 3, 4], [1, 2, 5, 6], [3, 5, 7, 8], [4, 6, 7, 8]),
    t!("d43333a", [1, 2, 3, 4], [1, 5, 6], [2, 5, 7], [3, 6, 8], [4, 7, 8]),
    t!("d43333b", [1, 2, 3, 4], [1, 2, 5], [3, 6, 7], [4, 6, 8], [5, 7, 8]),
];

/// A named invariant: a topology and the representation of each tensor.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub(crate) struct Invariant {
    pub(crate) name: &'static str,
    pub(crate) reps: Vec<Rep>,
}

/// Every topology with every representation assignment, by certificate.
struct Table {
    by_certificate: HashMap<Vec<u32>, u16>,
    entries: Vec<Invariant>,
}

static TABLE: LazyLock<Table> = LazyLock::new(build);

fn template_net(template: &Template, reps: &[Rep]) -> Net {
    Net::new(
        template
            .nodes
            .iter()
            .zip(reps)
            .map(|(node, rep)| Obj::D(*rep, node.iter().map(|&i| i - 1).collect()))
            .collect(),
    )
}

fn build() -> Table {
    let mut by_certificate = HashMap::new();
    let mut entries = Vec::new();
    for template in TEMPLATES {
        let n = template.nodes.len();
        // R before A, as color.h orders its labels.
        for mask in 0..(1u32 << n) {
            let reps: Vec<Rep> = (0..n).map(|k| if mask >> (n - 1 - k) & 1 == 1 { Rep::A } else { Rep::R }).collect();
            if template.nodes.iter().zip(&reps).any(|(node, rep)| *rep == Rep::A && node.len() % 2 == 1) {
                continue;
            }
            let (certificate, sign, _) = template_net(template, &reps).canonicalize();
            debug_assert_eq!(sign, 1);
            by_certificate.entry(certificate).or_insert_with(|| {
                entries.push(Invariant { name: template.name, reps });
                (entries.len() - 1) as u16
            });
        }
    }
    // color.h's dddff object, kept as an invariant where it has no relation:
    // symmetric in its three tensors (`symmetrize cOldddff`), R before A.
    for reps in [vec![Rep::R, Rep::R, Rep::R], vec![Rep::R, Rep::R, Rep::A]] {
        entries.push(Invariant { name: "dddff", reps });
    }
    Table { by_certificate, entries }
}

/// The invariant of a canonical pure symmetric-tensor network, if named.
pub(crate) fn lookup(certificate: &[u32]) -> Option<u16> {
    TABLE.by_certificate.get(certificate).copied()
}

/// All invariant names.
pub(crate) fn template_names() -> impl Iterator<Item = &'static str> {
    TEMPLATES.iter().map(|t| t.name).chain(["dddff"])
}

/// The canonical entry of a topology with the given representations.
pub(crate) fn invariant(name: &str, reps: &[Rep]) -> Option<u16> {
    if name == "dddff" {
        let mut sorted = reps.to_vec();
        sorted.sort();
        return TABLE.entries.iter().position(|e| e.name == "dddff" && e.reps == sorted).map(|p| p as u16);
    }
    let template = TEMPLATES.iter().find(|t| t.name == name)?;
    if template.nodes.len() != reps.len() {
        return None;
    }
    let (certificate, _, _) = template_net(template, reps).canonicalize();
    lookup(&certificate)
}

pub(crate) fn variable_of(id: u16) -> Option<&'static Invariant> {
    TABLE.entries.get((id - FIRST_NAMED) as usize)
}

fn variable(entry: u16) -> Poly {
    Poly::var(FIRST_NAMED + entry)
}

fn v(name: &str, reps: &[Rep]) -> Poly {
    value(invariant(name, reps).expect("dependency names a known topology"))
}

/// The value of a canonical entry after color.h's dependencies.
pub(crate) fn value(entry: u16) -> Poly {
    value_of(entry, &TABLE.entries[entry as usize])
}

fn ca(e: i8) -> Poly {
    Poly::mono(Mono::var(CA, e), C::one())
}

fn value_of(entry: u16, inv: &Invariant) -> Poly {
    use Rep::{A, R};
    let r = inv.reps.clone();
    let lin = |terms: &[(i64, i64, Poly)]| -> Poly {
        let mut out = Poly::zero();
        for (n, d, p) in terms {
            out.add_scaled(p, &C::rat(*n, *d));
        }
        out
    };
    match inv.name {
        "d433" if r[0] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d33", &[r[1], r[2]])),
        "d4433a" | "d4433c" if r[0] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d433", &[r[1], r[2], r[3]])),
        "d4433a" | "d4433c" if r[1] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d433", &[r[0], r[2], r[3]])),
        // A d_A^{(4)} with two legs on a d_R^{(3)} is cA^2/6 d_R^{(3)}, as for
        // d433(A, ..). color.h writes these two dependencies without the
        // factor, which its makeds has already taken out at that point.
        "d5443b" | "d5443d" if r[2] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d543", &[r[0], r[1], r[3]])),
        // Not in color.h, whose reductions never leave it: two d_A^{(4)}
        // sharing two legs, joined to two d_R^{(3)}. Exact for SU(N), where all
        // d_R^{abc} are proportional, hence for every group with a d^{(3)}.
        "d4433b" if r[0] == A && r[1] == A => lin(&[
            (1, 2, Poly::mono(Mono::var(NA, -1), C::one()).mul(&v("d33", &[r[2], r[3]])).mul(&v("d44", &[A, A]))),
            (-1, 144, ca(4).mul(&v("d33", &[r[2], r[3]]))),
        ]),
        // Not applied by color.h, which can leave d4433b(R,A,R,R) and
        // d653(A,R,R) where the reducer reaches cA d543, and dddff(R,R,A):
        // checked numerically, slot-wise and for representations spanning
        // every space of symmetric invariants, in SU(3) to SU(7) (d4433b),
        // SU(3) to SU(5) (d653; odd-rank tensors vanish outside SU(N)), and
        // SU(3..5), SO(7), SO(8), Sp(4), Sp(6) and G2 (dddff).
        "d4433b" if r[1] == A => lin(&[
            (1, 2, Poly::mono(Mono::var(NA, -1), C::one()).mul(&v("d33", &[r[2], r[3]])).mul(&v("d44", &[A, r[0]]))),
            (-1, 24, ca(2).mul(&v("d433", &[r[0], r[2], r[3]]))),
        ]),
        "d653" if r[0] == A => ca(1).scale(&C::rat(1, 10)).mul(&v("d543", &[r[1], A, r[2]])),
        "dddff" if r.iter().filter(|x| **x == R).count() == 2 => dddff_pair(R, R),
        "d6433a" if r[1] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d633", &[r[0], r[2], r[3]])),
        // color.h's makeds has -1/4 for the d433 terms of both d6433a
        // relations; the d_R^{(6)} sub-reduction is the one of d633, -1/40.
        "d6433a" if r[0] == R => lin(&[
            (-1, 8, ca(1).mul(&v("d4433c", &r))),
            (-1, 4, ca(1).mul(&v("d4433a", &r))),
            (-1, 40, Poly::var(I2R).mul(&ca(2)).mul(&v("d433", &[r[1], r[2], r[3]]))),
            (1, 1, Poly::mono(Mono::var(NR, -1), C::one()).mul(&v("d33", &[r[0], r[2]])).mul(&v("d433", &[r[1], r[0], r[3]]))),
        ]),
        "d6433a" => lin(&[
            (-1, 8, ca(1).mul(&v("d4433c", &r))),
            (-1, 4, ca(1).mul(&v("d4433a", &r))),
            (-1, 40, ca(3).mul(&v("d433", &[r[1], r[2], r[3]]))),
        ]),
        "d43333b" if r[0] == A => ca(2).scale(&C::rat(1, 6)).mul(&v("d3333", &r[1..])),
        "d5533b" => {
            let x = Poly::var(CR) + &ca(1).scale(&C::rat(-1, 4));
            v("d33", &[r[0], r[2]])
                .mul(&v("d33", &[r[1], r[3]]))
                .mul(&Poly::mono(Mono::var(NA, -1), C::one()))
                .mul(&x.mul(&x))
        }
        "d633" if r[0] == R => lin(&[
            (-3, 8, ca(1).mul(&v("d433", &r))),
            (-1, 40, Poly::var(I2R).mul(&ca(2)).mul(&v("d33", &[r[1], r[2]]))),
            (1, 1, Poly::mono(Mono::var(NR, -1), C::one()).mul(&v("d33", &[r[0], r[1]])).mul(&v("d33", &[R, r[2]]))),
        ]),
        "d633" => lin(&[
            (-3, 8, ca(1).mul(&v("d433", &r))),
            (-1, 40, ca(3).mul(&v("d33", &[r[1], r[2]]))),
        ]),
        "d743" if r[0] == R => lin(&[
            (-1, 2, ca(1).mul(&v("d543", &r))),
            (1, 1, Poly::mono(Mono::var(NR, -1), C::one()).mul(&v("d33", &[R, r[2]])).mul(&v("d44", &[r[0], r[1]]))),
            (-43, 120, Poly::mono(Mono::var(NA, -1), C::one()).mul(&v("d33", &[r[0], r[2]])).mul(&v("d44", &[A, r[1]]))),
            (-23, 1440, ca(2).mul(&v("d433", &[r[1], r[0], r[2]]))),
            (49, 60, v("d4433b", &[A, r[1], r[0], r[2]])),
        ]),
        "d66" if r == [A, A] => lin(&[
            (5, 8, v("d444", &[A, A, A])),
            (-7, 240, ca(2).mul(&v("d44", &[A, A]))),
            (-1, 864, ca(6).mul(&Poly::var(NA))),
        ]),
        _ => variable(entry),
    }
}

/// color.h's `dddff` object with representations `reps` (in any order): its
/// relation if at least one tensor is adjoint, else the invariant dddff(R,R,R).
pub(crate) fn dddff_value(reps: &[Rep]) -> Poly {
    match reps.iter().filter(|r| **r == Rep::R).count() {
        0 => dddff_pair(Rep::A, Rep::A),
        1 => dddff_pair(Rep::R, Rep::A),
        2 => dddff_pair(Rep::R, Rep::R),
        _ => value(invariant("dddff", reps).expect("dddff entries")),
    }
}

/// color.h's `dddff`, d_X^{a1a2a5a6} d_Y^{a3a4a7a8} d_Z^{a5a6a7a8}
/// f^{a1a3a9} f^{a2a4a9} (f(p1,p2,i) f(p1,p2,i) p1.p^2 p2.p^2 in color.h),
/// which is symmetric in X, Y, Z, with Z adjoint. color.h has this relation
/// for X = A; it holds for X, Y in R as well (see value_of).
pub(crate) fn dddff_pair(x: Rep, y: Rep) -> Poly {
    use Rep::A;
    let mut out = ca(3).mul(&v("d44", &[x, y])).scale(&C::rat(-2, 27));
    out.add_scaled(&ca(1).mul(&v("d444", &[x, y, A])), &C::rat(19, 15));
    out.add_scaled(&v("d644", &[A, x, y]), &C::rat(-8, 9));
    out
}
