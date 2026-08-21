//! Emit oneloop reductions for the cross-engine benchmark (gitignored, local-only).
//!
//! For a battery of scalar + dotted integral families over a fixed spacelike massive
//! geometry, run reduce() and print each reduction as
//!     coeff (a rational in `d`)  *  master(numeric args)
//! in a machine-readable form for `crosscheck.py`, which evaluates the masters with
//! OneLOop (avh_olo) and compares  sum_i c_i M_i  against a direct scipy integration.
//!
//! Geometry: Euclidean offsets in centi-units (r/100); the reducer's invariants are the
//! spacelike s_ij = -(r_i - r_j)^2, which feed OneLOop's three/four_point directly.
//!   cargo run --release --example emit_reductions -p oneloop > /tmp/oneloop_reductions.txt

use oneloop::masters::MasterIntegral;
use oneloop::symbols::S;
use oneloop::{Integral, IntegralFamily, Kinematics, Propagator, reduce};
use symbolica::atom::{Atom, AtomCore};
use symbolica::function;

// Two spacelike massive geometries (centi-unit integer offsets, centi-unit masses^2).
// Chosen so every C0/D0 is finite and OneLOop-evaluable (all separations spacelike).
const OFF_A: [[i64; 4]; 7] = [
    [0, 0, 0, 0],
    [30, 0, 0, 0],
    [10, 35, 0, 0],
    [0, 20, 40, 0],
    [5, 0, 15, 45],
    [22, 12, 0, 18],
    [8, 25, 30, 10],
];
const MSQ_A: [i64; 7] = [100, 121, 81, 144, 80, 105, 95];

const OFF_B: [[i64; 4]; 7] = [
    [0, 0, 0, 0],
    [40, 0, 0, 0],
    [15, 25, 0, 0],
    [5, 30, 20, 0],
    [0, 10, 35, 25],
    [18, 8, 12, 30],
    [12, 20, 5, 15],
];
const MSQ_B: [i64; 7] = [90, 110, 130, 75, 160, 100, 115];
// Massless internal line (line 0), rest massive; used with spacelike offsets -> IR-finite.
const MSQ_M: [i64; 7] = [0, 121, 81, 144, 80, 105, 95];

// NEAR-DEGENERATE geometry: lines 1,2 nearly collinear with line 0 (r0=0) -> the Gram of
// {r1,r2} is tiny (900 centi^4 vs |r1|^2|r2|^2 ~ 3.2e6), so 1/Gram numerator/dot coefficients
// blow up ~1e5x.  Stresses float catastrophic-cancellation in sum_i c_i M_i near a threshold.
// NB: line 3 carries a tiny 3rd-dim component (z=1) so the box's 3 external momenta span 3D
// with a *small* Gram (near-singular, valid) rather than being exactly rank-deficient. With
// [90,2,0,0] instead, q2==q3 exactly -> singular Gram -> gram_solve's Matrix::solve().unwrap()
// PANICS (Err(Inconsistent), reduce.rs:614) on any rank-2 numerator: a robustness gap to fix.
const OFF_D: [[i64; 4]; 7] = [
    [0, 0, 0, 0],
    [30, 0, 0, 0],
    [60, 1, 0, 0], // ~2*r1, tiny perpendicular -> near-collinear
    [90, 2, 1, 0],
    [45, 1, 3, 0],
    [22, 12, 0, 18],
    [8, 25, 30, 10],
];

fn mass_atom(msq_centi: i64) -> Atom {
    Atom::num(msq_centi) / Atom::num(100)
}

// Lexicographic pairwise invariants s_ij = -(r_i - r_j)^2 for lines 0..n, in centi^2 units
// (offsets are r/100, so (r_i-r_j)^2 = dist2_centi / 10000).
fn invariants(off: &[[i64; 4]; 7], n: usize) -> Vec<Atom> {
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let d: i64 = (0..4).map(|k| (off[i][k] - off[j][k]).pow(2)).sum();
            out.push(Atom::num(-d) / Atom::num(10000));
        }
    }
    out
}

fn family(
    off: &[[i64; 4]; 7],
    msq: &[i64; 7],
    n: usize,
    exps: Vec<i32>,
    numerator: Atom,
) -> IntegralFamily {
    IntegralFamily {
        propagators: (0..n)
            .map(|i| Propagator {
                momentum: Atom::Zero,
                mass_sq: mass_atom(msq[i]),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics {
            invariants: invariants(off, n),
        },
        targets: vec![Integral {
            propagator_exponents: exps,
            isp_exponents: vec![],
        }],
        numerator,
    }
}

// Family with EXPLICIT (lexicographic) invariants -- lets us inject timelike (positive) s_ij
// that no Euclidean offset geometry can produce.  invs are real values as (num, den) fractions.
fn family_inv(
    msq: &[i64; 7],
    n: usize,
    exps: Vec<i32>,
    invs: Vec<Atom>,
    numerator: Atom,
) -> IntegralFamily {
    IntegralFamily {
        propagators: (0..n)
            .map(|i| Propagator {
                momentum: Atom::Zero,
                mass_sq: mass_atom(msq[i]),
            })
            .collect(),
        isps: vec![],
        kinematics: Kinematics { invariants: invs },
        targets: vec![Integral {
            propagator_exponents: exps,
            isp_exponents: vec![],
        }],
        numerator,
    }
}

fn dot_ll() -> Atom {
    function!(S.dot, Atom::var(S.k), Atom::var(S.k))
}

fn dot_lq(j: usize) -> Atom {
    let q = symbolica::symbol!(format!("oneloop::q{}", j + 1));
    function!(S.dot, Atom::var(S.k), Atom::var(q))
}

fn master_line(m: &MasterIntegral) -> String {
    let s = |x: &Atom| x.expand().to_string();
    match m {
        MasterIntegral::Tadpole { m_sq } => format!("A0 {}", s(m_sq)),
        MasterIntegral::Bubble { p_sq, m1_sq, m2_sq } => {
            format!("B0 {} {} {}", s(p_sq), s(m1_sq), s(m2_sq))
        }
        MasterIntegral::Triangle {
            p1_sq,
            p2_sq,
            p12_sq,
            m1_sq,
            m2_sq,
            m3_sq,
        } => format!(
            "C0 {} {} {} {} {} {}",
            s(p1_sq),
            s(p2_sq),
            s(p12_sq),
            s(m1_sq),
            s(m2_sq),
            s(m3_sq)
        ),
        MasterIntegral::Box {
            p1_sq,
            p2_sq,
            p3_sq,
            p4_sq,
            s: si,
            t,
            m1_sq,
            m2_sq,
            m3_sq,
            m4_sq,
        } => format!(
            "D0 {} {} {} {} {} {} {} {} {} {}",
            s(p1_sq),
            s(p2_sq),
            s(p3_sq),
            s(p4_sq),
            s(si),
            s(t),
            s(m1_sq),
            s(m2_sq),
            s(m3_sq),
            s(m4_sq)
        ),
    }
}

fn emit(name: &str, off: &[[i64; 4]; 7], msq: &[i64; 7], n: usize, exps: Vec<i32>) {
    emit_num(name, off, msq, n, exps, Atom::num(1), "1");
}

fn emit_num(
    name: &str,
    off: &[[i64; 4]; 7],
    msq: &[i64; 7],
    n: usize,
    exps: Vec<i32>,
    numerator: Atom,
    num_label: &str,
) {
    let cfg = if std::ptr::eq(off, &OFF_A) { 0 } else { 1 };
    let fam = family(off, msq, n, exps.clone(), numerator);
    let r = reduce(&fam);
    let exps_s: Vec<String> = exps.iter().map(|e| e.to_string()).collect();
    println!(
        "PROCESS name={name} cfg={cfg} n={n} exps={} num={num_label}",
        exps_s.join(",")
    );
    // geometry (real = centi/100) so the driver can integrate the same family directly
    let off_s: Vec<String> = (0..n)
        .map(|i| {
            (0..4)
                .map(|k| format!("{}", off[i][k] as f64 / 100.0))
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    println!("OFFSETS {}", off_s.join(";"));
    let m_s: Vec<String> = (0..n)
        .map(|i| format!("{}", msq[i] as f64 / 100.0))
        .collect();
    println!("MASSES {}", m_s.join(","));
    for (c, m) in &r.terms {
        // coeff is a rational function of d only (kinematics are numeric); print for sympy.
        println!("TERM coeff=( {} ) {}", c.expand(), master_line(m));
    }
    println!("ENDPROCESS");
}

// Emit a family with explicit lex invariants sinv = [(num,den),...] (real units).  A positive
// s_ij is TIMELIKE; kept below threshold ((sqrt m_i + sqrt m_j)^2) so the integral stays real
// and scipy can integrate it with the Minkowski F = sum x m^2 - sum_{i<j} x_i x_j s_ij.
fn emit_tl(
    name: &str,
    msq: &[i64; 7],
    n: usize,
    exps: Vec<i32>,
    sinv: &[(i64, i64)],
    numerator: Atom,
    num_label: &str,
) {
    let invs: Vec<Atom> = sinv
        .iter()
        .map(|(a, b)| Atom::num(*a) / Atom::num(*b))
        .collect();
    let fam = family_inv(msq, n, exps.clone(), invs, numerator);
    let r = reduce(&fam);
    let exps_s: Vec<String> = exps.iter().map(|e| e.to_string()).collect();
    println!(
        "PROCESS name={name} cfg=9 n={n} exps={} num={num_label}",
        exps_s.join(",")
    );
    println!("OFFSETS {}", vec!["0,0,0,0"; n].join(";")); // dummy; direct value uses SINV
    let m_s: Vec<String> = (0..n)
        .map(|i| format!("{}", msq[i] as f64 / 100.0))
        .collect();
    println!("MASSES {}", m_s.join(","));
    let sv: Vec<String> = sinv
        .iter()
        .map(|(a, b)| format!("{}", *a as f64 / *b as f64))
        .collect();
    println!("SINV {}", sv.join(",")); // signed lex invariants (real)
    for (c, m) in &r.terms {
        println!("TERM coeff=( {} ) {}", c.expand(), master_line(m));
    }
    println!("ENDPROCESS");
}

fn main() {
    gammalooprs::initialisation::initialise().expect("activate Symbolica license");

    for (off, msq) in [(&OFF_A, &MSQ_A), (&OFF_B, &MSQ_B)] {
        // dotted bubble (N=2): the minimal FJT index-lowering leaf
        emit("dot_N2", off, msq, 2, vec![2, 1]);
        emit("dot2_N2", off, msq, 2, vec![3, 1]);
        // scalar N-point: triangle .. heptagon (finite; vNV / degenerate-Cayley path for N>4)
        for n in 3..=7 {
            emit(&format!("scalar_N{n}"), off, msq, n, vec![1; n]);
        }
        // dotted (one raised power) N=3..6: the FJT dimension-shift / index-lowering path
        for n in 3..=6 {
            let mut exps = vec![1; n];
            exps[0] = 2;
            emit(&format!("dot_N{n}"), off, msq, n, exps);
        }
        // doubly dotted triangle + box
        emit("dot2_N3", off, msq, 3, vec![2, 2, 1]);
        emit("dot2_N4", off, msq, 4, vec![2, 1, 2, 1]);
        // numerator k^2 = dot(k,k): convention-free (r0=0, so k^2 = D0 + m0^2). Tests the
        // separate triangle_numerator / box_numerator / ngon_numerator path.
        for n in 3..=6 {
            emit_num(
                &format!("numll_N{n}"),
                off,
                msq,
                n,
                vec![1; n],
                dot_ll(),
                "ll",
            );
        }
        // (k^2)^2: rank-4 convention-free numerator; finite (and OneLOop-oracle-able) on
        // pentagon/hexagon (A >= 5).  Exercises the irreducible tensor / ISP machinery.
        let ll2 = &dot_ll() * &dot_ll();
        emit_num("numll2_N4", off, msq, 4, vec![1; 4], ll2.clone(), "ll2d"); // rank-4 box: UV-divergent
        emit_num("numll2_N5", off, msq, 5, vec![1; 5], ll2.clone(), "ll2");
        emit_num("numll2_N6", off, msq, 6, vec![1; 6], ll2.clone(), "ll2");
        // rank-6 (k^2)^3: pentagon (UV-div, A=5) + hexagon (finite, A=6)
        let ll3 = &ll2 * &dot_ll();
        emit_num("numll3_N5", off, msq, 5, vec![1; 5], ll3.clone(), "ll3d");
        emit_num("numll3_N6", off, msq, 6, vec![1; 6], ll3, "ll3");
        // dot(k, q_j) numerators: q_j = r_j - r_{j-1} (leg differences). Rank-1 and mixed.
        for n in 3..=5 {
            emit_num(
                &format!("numq1_N{n}"),
                off,
                msq,
                n,
                vec![1; n],
                dot_lq(0),
                "q1",
            );
        }
        emit_num("numq2_N4", off, msq, 4, vec![1; 4], dot_lq(1), "q2");
        emit_num(
            "numq1q2_N4",
            off,
            msq,
            4,
            vec![1; 4],
            &dot_lq(0) * &dot_lq(1),
            "q1q2",
        );
        emit_num(
            "numllq1_N4",
            off,
            msq,
            4,
            vec![1; 4],
            &dot_ll() * &dot_lq(0),
            "llq1",
        );
        // N>4 NUMERATOR PATH (ngon_numerator): dot(k,q_j) for j>=2 + mixed products, on a
        // pentagon/hexagon.  Highest-risk untested area (the box rule_q2 bug's N>4 analog).
        emit_num("numq2_N5", off, msq, 5, vec![1; 5], dot_lq(1), "q2");
        emit_num("numq3_N5", off, msq, 5, vec![1; 5], dot_lq(2), "q3");
        emit_num("numq4_N5", off, msq, 5, vec![1; 5], dot_lq(3), "q4");
        emit_num(
            "numq1q2_N5",
            off,
            msq,
            5,
            vec![1; 5],
            &dot_lq(0) * &dot_lq(1),
            "q1q2",
        );
        emit_num(
            "numq2q3_N5",
            off,
            msq,
            5,
            vec![1; 5],
            &dot_lq(1) * &dot_lq(2),
            "q2q3",
        );
        emit_num(
            "numllq2_N5",
            off,
            msq,
            5,
            vec![1; 5],
            &dot_ll() * &dot_lq(1),
            "llq2",
        );
        emit_num("numq3_N6", off, msq, 6, vec![1; 6], dot_lq(2), "q3");
        // MIXED dotted + numerator (box_numerator with raised powers): finite.
        emit_num(
            "numq1_dotN4",
            off,
            msq,
            4,
            vec![2, 1, 1, 1],
            dot_lq(0),
            "q1",
        );
        emit_num(
            "numq1q2_dotN4",
            off,
            msq,
            4,
            vec![2, 1, 1, 1],
            &dot_lq(0) * &dot_lq(1),
            "q1q2",
        );
        emit_num(
            "numll_dotN4",
            off,
            msq,
            4,
            vec![2, 1, 1, 1],
            dot_ll(),
            "lld",
        );
        // rank-4 on a TRIANGLE (k^2)^2: UV-divergent (quadratic -> single 1/eps in dim reg).
        emit_num(
            "numll2_N3",
            off,
            msq,
            3,
            vec![1; 3],
            &dot_ll() * &dot_ll(),
            "ll2d",
        );
        // triple dot (higher exponent) scalar + a dotted numerator combo.
        emit("dot3_N3", off, msq, 3, vec![3, 1, 1]);
        emit_num("numq1_dot2N3", off, msq, 3, vec![2, 2, 1], dot_lq(0), "q1");
        // rank-2 EXTERNAL numerator (k.q3)^2 on a BOX: genuine 3x3 Gram-inverse reduction,
        // UV-finite (rank-2 on a box).  Baseline for the near-degenerate stress below.
        emit_num(
            "numq3q3_N4",
            off,
            msq,
            4,
            vec![1; 4],
            &dot_lq(2) * &dot_lq(2),
            "q3q3",
        );
    }

    // ---- KINEMATIC regime: massless internal line (line 0 = 0), spacelike externals ----
    // With off-shell externals a single massless line stays IR-finite (the full integral AND
    // the individual massless B0/C0/D0 masters), so scipy + pole-cancellation still apply.
    // Exercises the m=0 reduction path: tadpole A0(0)=0, massless masters.
    for n in 3..=6 {
        emit(&format!("mless_N{n}"), &OFF_A, &MSQ_M, n, vec![1; n]);
    }
    emit("mless_dot_N4", &OFF_A, &MSQ_M, 4, vec![1, 2, 1, 1]); // dot a MASSIVE line
    emit_num(
        "mless_numll_N4",
        &OFF_A,
        &MSQ_M,
        4,
        vec![1; 4],
        dot_ll(),
        "ll",
    );
    emit_num(
        "mless_numq1_N4",
        &OFF_A,
        &MSQ_M,
        4,
        vec![1; 4],
        dot_lq(0),
        "q1",
    );

    // ---- NEAR-DEGENERATE kinematics (tiny Gram): float catastrophic-cancellation stress ----
    // Coefficients of the dot/numerator reductions carry 1/Gram ~ 1e5, so sum_i c_i M_i is a
    // large near-cancellation.  scalar (coeff 1) is exact; the real test is dot/numerator.
    emit("dg_scalar_N3", &OFF_D, &MSQ_A, 3, vec![1; 3]);
    emit("dg_dot_N3", &OFF_D, &MSQ_A, 3, vec![2, 1, 1]);
    emit_num("dg_numll_N3", &OFF_D, &MSQ_A, 3, vec![1; 3], dot_ll(), "ll");
    emit_num(
        "dg_numq1_N3",
        &OFF_D,
        &MSQ_A,
        3,
        vec![1; 3],
        dot_lq(0),
        "q1",
    );
    emit_num(
        "dg_numll2_N3",
        &OFF_D,
        &MSQ_A,
        3,
        vec![1; 3],
        &dot_ll() * &dot_ll(),
        "ll2d",
    );
    // near-degenerate BOX too (lines 1..3 nearly collinear): dot + k^2 numerators.
    emit("dg_scalar_N4", &OFF_D, &MSQ_A, 4, vec![1; 4]);
    emit_num("dg_numll_N4", &OFF_D, &MSQ_A, 4, vec![1; 4], dot_ll(), "ll");
    emit_num(
        "dg_numq1_N4",
        &OFF_D,
        &MSQ_A,
        4,
        vec![1; 4],
        dot_lq(0),
        "q1",
    );
    // genuine 1/Gram rank-2 external numerator (k.q1)^2 on the near-degenerate triangle+box:
    // THE catastrophic-cancellation stress (coeffs ~1/Gram vs the numq1q1_N3 baseline).
    // (k.q3)^2 and (k.q2)(k.q3) on the box WEIGHT the near-collinear q2~=q3 direction -> hit
    // the small 3x3 Gram(q1,q2,q3), a genuine 1/Gram reduction, and stay UV-finite (box).
    // THE catastrophic-cancellation stress vs the numq3q3_N4 baseline.
    emit_num(
        "dg_numq1q1_N4",
        &OFF_D,
        &MSQ_A,
        4,
        vec![1; 4],
        &dot_lq(0) * &dot_lq(0),
        "q1q1",
    );
    emit_num(
        "dg_numq3q3_N4",
        &OFF_D,
        &MSQ_A,
        4,
        vec![1; 4],
        &dot_lq(2) * &dot_lq(2),
        "q3q3",
    );
    emit_num(
        "dg_numq2q3_N4",
        &OFF_D,
        &MSQ_A,
        4,
        vec![1; 4],
        &dot_lq(1) * &dot_lq(2),
        "q2q3",
    );
    // NB: an EXACTLY-singular box (OFF_DS, q2==q3) now surfaces a clear "singular Gram matrix"
    // error from reduce() rather than reducing, so it is intentionally not emitted here.

    // ---- TIMELIKE (positive s_ij) but BELOW threshold: integral stays REAL, scipy-checkable ----
    // All masses^2 = 2.0 (threshold (sqrt2+sqrt2)^2 = 8); one invariant is timelike (+), rest
    // spacelike (-).  Tests master args + coefficient evaluation at timelike kinematics.
    // DISTINCT masses (so the crosscheck can identify each pinched sub-topology's lines).
    let m2: [i64; 7] = [200, 210, 190, 220, 180, 205, 195];
    let tri = [(1, 2), (-3, 10), (-1, 5)]; // lex s01=+0.5 (timelike), s02=-0.3, s12=-0.2
    emit_tl("tl_scalar_N3", &m2, 3, vec![1; 3], &tri, Atom::num(1), "1");
    emit_tl("tl_dot_N3", &m2, 3, vec![2, 1, 1], &tri, Atom::num(1), "1");
    let box_inv = [(3, 5), (-3, 10), (-1, 5), (-2, 5), (-1, 4), (-3, 10)]; // s01=+0.6 timelike
    emit_tl(
        "tl_scalar_N4",
        &m2,
        4,
        vec![1; 4],
        &box_inv,
        Atom::num(1),
        "1",
    );
    emit_tl(
        "tl_dot_N4",
        &m2,
        4,
        vec![2, 1, 1, 1],
        &box_inv,
        Atom::num(1),
        "1",
    );
    // NEAR-threshold: masses^2~0.5 (threshold ~2.1), s01=+1.5 -> C0 real but close to the cut.
    let mlt: [i64; 7] = [50, 55, 45, 52, 48, 53, 47];
    let near = [(3, 2), (-3, 10), (-1, 5)]; // s01=+1.5 (below threshold ~2.1)
    emit_tl("tl_near_N3", &mlt, 3, vec![1; 3], &near, Atom::num(1), "1");
    emit_tl(
        "tl_neardot_N3",
        &mlt,
        3,
        vec![2, 1, 1],
        &near,
        Atom::num(1),
        "1",
    );

    // ---- ABOVE threshold (complex, absorptive): checked via the IBP mass-derivative ----
    // I[a0=2] = d/dm0^2 I[scalar] (line 0 doubled), the RHS finite-differenced from OneLOop's
    // COMPLEX master.  "tlx_" names route to that check (scipy can't cross F=0 above threshold).
    let abv3 = [(10, 1), (-3, 10), (-1, 5)]; // s01=+10 > threshold ~8.2 -> C0 complex
    emit_tl(
        "tlx_dot_N3",
        &m2,
        3,
        vec![2, 1, 1],
        &abv3,
        Atom::num(1),
        "1",
    );
    let abv4 = [(-3, 10), (12, 1), (-1, 5), (-2, 5), (-1, 4), (-3, 10)]; // s02(diag)=+12 -> D0 complex
    emit_tl(
        "tlx_dot_N4",
        &m2,
        4,
        vec![2, 1, 1, 1],
        &abv4,
        Atom::num(1),
        "1",
    );

    // ---- RANK >= 3 MIXED tensor numerators (the open generalization) ----
    // High-rank numerators with DISTINCT external directions (not just isotropic (k^2)^p).
    // Chosen UV-finite (rank r < 2(N-2)) so the moment/scipy oracle validates the finite
    // part directly: box r=3; pentagon r=3,4; hexagon r=3,5.  Two spacelike geometries.
    // Placed last so any panic in this frontier area cannot truncate the families above.
    for (off, msq) in [(&OFF_A, &MSQ_A), (&OFF_B, &MSQ_B)] {
        let q1q2q3 = &(&dot_lq(0) * &dot_lq(1)) * &dot_lq(2);
        // box (N=4), rank-3
        emit_num("mix_q1q2q3_N4", off, msq, 4, vec![1; 4], q1q2q3.clone(), "q1q2q3");
        emit_num(
            "mix_q1q1q2_N4",
            off,
            msq,
            4,
            vec![1; 4],
            &(&dot_lq(0) * &dot_lq(0)) * &dot_lq(1),
            "q1q1q2",
        );
        // pentagon (N=5), rank-3 and rank-4
        emit_num("mix_q1q2q3_N5", off, msq, 5, vec![1; 5], q1q2q3.clone(), "q1q2q3");
        emit_num(
            "mix_q1q2q3q4_N5",
            off,
            msq,
            5,
            vec![1; 5],
            &(&(&dot_lq(0) * &dot_lq(1)) * &dot_lq(2)) * &dot_lq(3),
            "q1q2q3q4",
        );
        emit_num(
            "mix_llq1q2_N5",
            off,
            msq,
            5,
            vec![1; 5],
            &(&dot_ll() * &dot_lq(0)) * &dot_lq(1),
            "llq1q2",
        );
        // hexagon (N=6), rank-3 and rank-5
        emit_num("mix_q1q2q3_N6", off, msq, 6, vec![1; 6], q1q2q3.clone(), "q1q2q3");
        emit_num(
            "mix_llq1q2q3_N6",
            off,
            msq,
            6,
            vec![1; 6],
            &(&(&dot_ll() * &dot_lq(0)) * &dot_lq(1)) * &dot_lq(2),
            "llq1q2q3",
        );
        // dotted (raised power) x rank-3 mixed: box [2,1,1,1], UV-finite (r < 6).
        emit_num(
            "mix_q1q2q3_dotN4",
            off,
            msq,
            4,
            vec![2, 1, 1, 1],
            q1q2q3.clone(),
            "q1q2q3",
        );
        emit_num(
            "mix_llq1q2_dotN4",
            off,
            msq,
            4,
            vec![2, 1, 1, 1],
            &(&dot_ll() * &dot_lq(0)) * &dot_lq(1),
            "llq1q2",
        );
        // rank-6 mixed on hexagon: k^2 q1 q2 q3 q4.
        emit_num(
            "mix_llq1q2q3q4_N6",
            off,
            msq,
            6,
            vec![1; 6],
            &(&(&(&dot_ll() * &dot_lq(0)) * &dot_lq(1)) * &dot_lq(2)) * &dot_lq(3),
            "llq1q2q3q4",
        );
        // heptagon (N=7) MIXED tensor: N=7's 6 external momenta are linearly dependent in 4D,
        // so the reducible-direction Gram is rank-deficient -> handled by gram_solve's subset
        // pseudo-inverse (validated here vs the moment oracle). See docs/04-frontier.md.
        emit_num("mix_q1q2q3_N7", off, msq, 7, vec![1; 7], q1q2q3.clone(), "q1q2q3");
    }

    // (on-shell / SINV scaffolding — ms_/ggh_/ggho_/reg_ — removed: the offset-based
    //  scipy/tensor oracle can't validate SINV kinematics. The regularization fix is
    //  unit-tested in src (on_shell_massless_triangle_rank2_regularizes); results in
    //  benchmarks/madloop_reference.md.)
}
