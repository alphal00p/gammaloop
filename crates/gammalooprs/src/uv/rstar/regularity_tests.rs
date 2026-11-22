//! IRH-DEPTH round 5 ([`super::super::regularity`]): the review's witnesses
//! of 2026-11-21 (Astra, IRH-DEPTH/native/inputs, byte for byte) and one
//! pole witness per entry of the closed factor list, on B_self_energy_L3_006
//! as the round-2/3 tests ([`super`]).
//!
//! The independent value is −117/64 C_A³/ε (the review's subtraction on the
//! same tensor denominator, and the grid's three methods at k3j1); the
//! round-3 child route gave +117/32 on the hidden tensor and on the
//! numerator prefactor. Each witness runs the round-3 planner
//! ([`regularity::as_round3`]) next to the fixed one and pins its wrong value,
//! so it fails without the certificate.

use eyre::{Result, ensure};
use symbolica::atom::{Atom, AtomCore};

use super::super::{
    closed::with_ir_hard,
    depth::{self, DepthCounts, DepthTest},
    regularity::{self, Factor},
};
use super::{Case, EPS, FIXED, L3_006, ensure_same, fnv, pinned, q, run, scaled};

const DIM: &str = "gammalooprs::{}::dim";
const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";

/// The review's closed metric trace with fresh Lorentz slots (make_inputs.py):
/// `−2/(g_ab g^ab − 4)`, equal to 1/ε in d = 4−2ε.
fn tensor() -> String {
    let slot = |n: usize| format!("{MINK}(4,{HEDGE}({n},1))");
    format!("(-2)/(({METRIC}({},{}))^2-4)", slot(81), slot(82))
}

const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";

/// A closed Dirac trace with fresh slots: `−8/(tr(γ_a γ^a) − 16)`, with
/// `tr(γ_a γ^a) = d tr 1 = 4d`, equal to 1/ε in d = 4−2ε (MATRIX-R45,
/// testdata/irh_matrix/gen_matrix.py).
fn dirac_trace() -> String {
    format!("(-8)/({}-16)", dirac_product())
}

/// `tr(γ_a γ^a)` with fresh slots (bispinor 83, 84; Lorentz 85).
fn dirac_product() -> String {
    let bis = |n: usize| format!("{BIS}(4,{HEDGE}({n},1))");
    let mink = format!("{MINK}(4,{HEDGE}(85,1))");
    format!(
        "{GAMMA}({},{},{mink})*{GAMMA}({},{},{mink})",
        bis(83),
        bis(84),
        bis(84),
        bis(83)
    )
}

/// ε⁻¹ in the spellings of the mechanical matrix: explicit, a `1/(d−4)`
/// form, a closed metric trace (the review's), a closed Dirac trace
/// (MATRIX-R45).
#[derive(Clone, Copy, Debug)]
enum Spelling {
    Epsilon,
    DimensionMinusFour,
    Trace,
    DiracTrace,
}

impl Spelling {
    fn text(self) -> String {
        match self {
            Self::Epsilon => format!("({EPS})^(-1)"),
            Self::DimensionMinusFour => format!("(-2)/({DIM}-4)"),
            Self::Trace => tensor(),
            Self::DiracTrace => dirac_trace(),
        }
    }
}

const SPELLINGS: [Spelling; 4] = [
    Spelling::Epsilon,
    Spelling::DimensionMinusFour,
    Spelling::Trace,
    Spelling::DiracTrace,
];

/// The span of the first `field="…"` of `line` not preceded by a lowercase
/// letter or `_` (the fixture scripts' `(?<![a-z_])field="([^"]+)"`).
fn field_span(line: &str, field: &str) -> Option<(usize, usize)> {
    let key = format!("{field}=\"");
    let mut from = 0;
    while let Some(found) = line[from..].find(&key) {
        let start = from + found;
        let before = line[..start].chars().last();
        if before.is_some_and(|c| c.is_ascii_lowercase() || c == '_') {
            from = start + key.len();
            continue;
        }
        let body = start + key.len();
        let end = body + line[body..].find('"')?;
        return Some((body, end));
    }
    None
}

/// `dot` with `field` of the line starting with `  {term} [` (or
/// `  {term} -> …` for an edge) rewritten by `transform` (make_inputs.py's
/// `modify`).
fn modify(dot: &str, term: &str, field: &str, transform: impl Fn(&str) -> String) -> String {
    let mut count = 0;
    let out = dot
        .split_inclusive('\n')
        .map(|line| {
            if line.starts_with(&format!("  {term} ")) && count == 0 {
                let (body, end) = field_span(line, field).expect("the field");
                count += 1;
                format!("{}{}{}", &line[..body], transform(&line[body..end]), &line[end..])
            } else {
                line.to_owned()
            }
        })
        .collect();
    assert_eq!(count, 1, "{term} {field}");
    out
}

/// Replace the tail `tail` of a field by `with` (make_inputs.py).
fn retail(tail: &str, with: &str) -> impl Fn(&str) -> String {
    let (tail, with) = (tail.to_owned(), with.to_owned());
    move |x: &str| {
        assert!(x.ends_with(&tail), "{x} ends with {tail}");
        format!("{}{with}", &x[..x.len() - tail.len()])
    }
}

/// `({x})*({factor})`.
fn times(factor: &str) -> impl Fn(&str) -> String {
    let factor = factor.to_owned();
    move |x: &str| format!("({x})*({factor})")
}

fn hidden() -> String {
    scaled(L3_006, &[("v0", 3), ("v1", -1)])
}
fn hard_v0() -> String {
    scaled(L3_006, &[("v0", 3)])
}
fn quotient_v1() -> String {
    scaled(L3_006, &[("v1", -1)])
}

/// The review's tensor-hidden: hidden with v1's ε⁻¹ spelled as the trace.
fn tensor_hidden() -> String {
    let t = format!("*({})", tensor());
    let dot = modify(&hidden(), "v1", "numerator", retail(&format!("*({EPS})^(-1)"), &t));
    modify(
        &dot,
        "graph",
        "numerator",
        retail(&format!("*({EPS})^(2)"), &format!("*({EPS})^(3){t}")),
    )
}

/// The review's visible neighbour: quotient-v1 with the trace spelling.
fn tensor_visible() -> String {
    let t = format!("*({})", tensor());
    let dot = modify(&quotient_v1(), "v1", "numerator", retail(&format!("*({EPS})^(-1)"), &t));
    modify(&dot, "graph", "numerator", retail(&format!("*({EPS})^(-1)"), &t))
}

/// The review's dim-quotient: hidden with v1's ε⁻¹ spelled −2/(d−4).
fn dim_quotient() -> String {
    let d = format!("*(-2)/({DIM}-4)");
    let dot = modify(&hidden(), "v1", "numerator", retail(&format!("*({EPS})^(-1)"), &d));
    modify(
        &dot,
        "graph",
        "numerator",
        retail(&format!("*({EPS})^(2)"), &format!("*({EPS})^(3){d}")),
    )
}

/// The review's hard-pole-<field>: hard-v0 with the graph's `field` over ε.
fn hard_pole(field: &str) -> String {
    modify(&hard_v0(), "graph", field, |x| format!("({x})/({EPS})"))
}

#[test]
fn the_witnesses_are_the_reviewed_inputs() {
    for (name, fixture, hash, length) in [
        ("tensor-hidden", tensor_hidden(), 0x968d_72f7_fd27_8a4c_u64, 47874),
        ("tensor-visible", tensor_visible(), 0x2bf5_cd73_73d1_a349, 47820),
        ("dim-quotient", dim_quotient(), 0xf653_c052_3a4a_6190, 47394),
        (
            "hard-pole-numerator_prefactor",
            hard_pole("numerator_prefactor"),
            0x7cdf_2feb_2653_d8c1,
            47356,
        ),
        (
            "hard-pole-overall_factor",
            hard_pole("overall_factor"),
            0x01c7_042b_872f_dcfd,
            47356,
        ),
    ] {
        assert_eq!((fnv(&fixture), fixture.len()), (hash, length), "{name}");
    }
}

/// The child route of `case` with the graph's numerator prefactor applied in
/// the root projector, as beta's callback does
/// (`project_physical_root_mode`: `tensor · σ(global_prefactor.num)`), on
/// `route` ("child" or "closed").
fn projected(case: &Case, route: &'static str, test: DepthTest) -> Result<(Atom, DepthCounts)> {
    let prefactor = case.graph.global_prefactor.num.clone();
    let projector = move |tensor: &Atom| -> Result<Atom> { Ok(tensor * &prefactor) };
    let (value, counts) = depth::with_depth_test(test, || {
        depth::collect(|| {
            with_ir_hard(route, || {
                case.plan.projected_counterterm(
                    &case.graph,
                    case.model(),
                    &case.vakint,
                    &case.settings,
                    &projector,
                )
            })
        })
    });
    eprintln!("IRH_R5_PROJECTED route={route} counts={counts:?}");
    Ok((value?, counts))
}

/// A child run certified without a refusal, the bound raised to the hidden
/// pole: no retry, no exceeded alarm.
fn ensure_raised(label: &str, counts: &DepthCounts) -> Result<()> {
    ensure!(
        counts.corrections >= 1
            && counts.bounded >= 1
            && counts.exceeded == 0
            && counts.refused == 0
            && counts.retried == 0,
        "{label}: {counts:?}"
    );
    Ok(())
}

/// Finding 1, the hidden tensor (review native/inputs/tensor-hidden): the
/// quotient vertex v1 carries −2/(g_ab g^ab − 4) = 1/ε with fresh indices
/// and v0·ε³ hides every row. The trace is contracted, B = 2, the child
/// route equals today's −117/64 C_A³/ε; round 3 read B = 1 and gave
/// +117/32 with no refusal.
///
/// Every witness runs through a root projector (here the graph prefactor,
/// 1), as beta does: the projected path regularises `d` after the
/// projection, before the final pole map (`ir_outer_integrated`), so a
/// symbolic `d` from a numerator (a contracted trace, `1/(d−4)`) is read in
/// `4−2ε`; the unprojected `counterterm` keeps it symbolic through its `K`
/// (both routes, consistently: 0 here), outside beta's use.
#[test]
fn the_hidden_tensor_pole_is_read() -> Result<()> {
    let case = Case::new(&tensor_hidden(), None)?;
    let (today, _) = projected(&case, "closed", FIXED)?;
    pinned("tensor-hidden today", &today, &[(-1, q(-117, 64))], true)?;
    let (child, counts) = projected(&case, "child", FIXED)?;
    ensure_same("tensor-hidden", &child, &today)?;
    ensure_raised("tensor-hidden", &counts)?;
    let (round3, counts3) = regularity::as_round3(|| projected(&case, "child", FIXED))?;
    ensure!(
        counts3.refused == 0 && counts3.exceeded == 0,
        "round 3 accepted it silently: {counts3:?}"
    );
    pinned("tensor-hidden round3", &round3, &[(-1, q(117, 32))], true)
}

/// The visible neighbour (review native/inputs/tensor-visible): v1's trace
/// alone. Round 3 measured a pole of order 2 over its bound 1 (exceeded,
/// declined to today's path); now the bound is 2: no alarm, 117/32 ε⁻⁴.
#[test]
fn the_visible_tensor_pole_is_bounded() -> Result<()> {
    let case = Case::new(&tensor_visible(), None)?;
    let (today, _) = projected(&case, "closed", FIXED)?;
    pinned("tensor-visible", &today, &[(-4, q(117, 32))], false)?;
    let (child, counts) = projected(&case, "child", FIXED)?;
    ensure_same("tensor-visible", &child, &today)?;
    ensure!(
        counts.bounded >= 1 && counts.exceeded == 0 && counts.refused == 0,
        "{counts:?}"
    );
    let (round3, counts3) = regularity::as_round3(|| projected(&case, "child", FIXED))?;
    ensure_same("tensor-visible round3", &round3, &today)?;
    ensure!(
        counts3.exceeded >= 1 && counts3.reason == Some("bound_exceeded"),
        "{counts3:?}"
    );
    Ok(())
}

/// The scalar-dimension neighbour (review native/inputs/dim-quotient):
/// v1's 1/ε as −2/(d−4); read on both planners, −117/64.
#[test]
fn the_dimension_spelling_is_read() -> Result<()> {
    let case = Case::new(&dim_quotient(), None)?;
    let (today, _) = projected(&case, "closed", FIXED)?;
    let (child, counts) = projected(&case, "child", FIXED)?;
    ensure_same("dim-quotient", &child, &today)?;
    ensure_raised("dim-quotient", &counts)?;
    pinned("dim-quotient", &child, &[(-1, q(-117, 64))], true)
}

/// Finding 2, the public prefactor (review
/// native/inputs/hard-pole-numerator_prefactor): hard-v0 with
/// `numerator_prefactor` = 1/ε, applied in the root projector as beta does.
/// Its excess is in B (2): −117/64 on both routes; round 3 did not read
/// it (B = 1) and gave +117/32.
#[test]
fn the_numerator_prefactor_pole_is_read() -> Result<()> {
    let case = Case::new(&hard_pole("numerator_prefactor"), None)?;
    let (today, _) = projected(&case, "closed", FIXED)?;
    pinned("prefactor today", &today, &[(-1, q(-117, 64))], true)?;
    let (child, counts) = projected(&case, "child", FIXED)?;
    ensure_same("prefactor", &child, &today)?;
    ensure_raised("prefactor", &counts)?;
    let (round3, counts3) = regularity::as_round3(|| projected(&case, "child", FIXED))?;
    ensure!(counts3.refused == 0, "{counts3:?}");
    pinned("prefactor round3", &round3, &[(-1, q(117, 32))], true)
}

/// One matrix cell (projected as beta): the child route equals today's path
/// exactly, with no exceeded alarm (a refusal is allowed and counted: then
/// the value is today's by construction).
fn cell(label: &str, dot: &str) -> Result<Option<(Atom, DepthCounts)>> {
    let case = Case::new(dot, None)?;
    let today = projected(&case, "closed", FIXED);
    let child = projected(&case, "child", FIXED);
    let (today, child, counts) = match (today, child) {
        (Ok((today, _)), Ok((child, counts))) => (today, child, counts),
        (Err(today), Err(child)) => {
            // Loud on both routes: no value is returned.
            eprintln!("IRH_R5_CELL {label} loud_error today={today:#} child={child:#}");
            return Ok(None);
        }
        (Err(today), Ok(_)) => {
            eyre::bail!("{label}: today's path fails ({today:#}) but the child route returns a value")
        }
        (Ok(_), Err(child)) => return Err(child),
    };
    if !super::same(&child, &today)? {
        eprintln!(
            "IRH_R5_CELL_MISMATCH {label} child_pole={} today_pole={} counts={counts:?}",
            super::laurent(&super::record(&child)?, -1)?,
            super::laurent(&super::record(&today)?, -1)?,
        );
    }
    ensure_same(label, &child, &today)?;
    ensure!(counts.exceeded == 0, "{label}: {counts:?}");
    eprintln!("IRH_R5_CELL {label} counts={counts:?}");
    // Before round 5 (not asserted, the matrix's "before" column).
    let (round3, counts3) = regularity::as_round3(|| projected(&case, "child", FIXED))?;
    eprintln!(
        "IRH_R5_CELL_ROUND3 {label} same={} refused={} exceeded={}",
        super::same(&round3, &today)?,
        counts3.refused,
        counts3.exceeded
    );
    Ok(Some((today, counts)))
}

/// The matrix for one attribute: the three spellings give the same pole part
/// (in record units) with the child route equal to today's path; `pin` is
/// the explicit spelling's ε⁻¹ value when known independently.
fn matrix(
    label: &str,
    fixture: impl Fn(Spelling) -> String,
    pin: Option<Atom>,
) -> Result<Vec<Option<DepthCounts>>> {
    matrix_over(label, &SPELLINGS, fixture, pin)
}

/// [`matrix`] over `spellings`.
fn matrix_over(
    label: &str,
    spellings: &[Spelling],
    fixture: impl Fn(Spelling) -> String,
    pin: Option<Atom>,
) -> Result<Vec<Option<DepthCounts>>> {
    let mut first = None::<Atom>;
    let mut all = Vec::new();
    for &spelling in spellings {
        let cell_label = format!("{label} {spelling:?}");
        let Some((today, counts)) = cell(&cell_label, &fixture(spelling))? else {
            all.push(None);
            continue;
        };
        all.push(Some(counts));
        let poles = super::record(&today)?
            .series(crate::utils::GS.dim_epsilon, 0, -1)
            .map_err(|error| eyre::eyre!("{error}"))?
            .to_atom()
            .expand();
        if let Some(pin) = &pin {
            pinned(&cell_label, &today, &[(-1, pin.clone())], false)?;
        }
        match &first {
            None => first = Some(poles),
            Some(first) => ensure!(
                (first - &poles).expand().together().is_zero(),
                "{cell_label}: poles {poles} differ from the explicit spelling's {first}"
            ),
        }
    }
    Ok(all)
}

/// Entry quotient_vertex_numerator: v1 · ε⁻¹ (three spellings) with v0·ε³.
#[test]
fn the_quotient_vertex_pole_is_read_in_every_spelling() -> Result<()> {
    matrix(
        "v1",
        |spelling| {
            let dot = modify(&hard_v0(), "v1", "numerator", times(&spelling.text()));
            modify(&dot, "graph", "numerator", times(&spelling.text()))
        },
        Some(q(-117, 64)),
    )
    .map(|_| ())
}

/// Entry quotient_edge_numerator: the quotient edge 4 (v1 → v3) · ε⁻¹ with
/// v0·ε³.
#[test]
fn the_quotient_edge_pole_is_read_in_every_spelling() -> Result<()> {
    matrix(
        "edge4",
        |spelling| {
            let dot = modify(&hard_v0(), "v1 ->", "numerator", times(&spelling.text()));
            modify(&dot, "graph", "numerator", times(&spelling.text()))
        },
        None,
    )
    .map(|_| ())
}

/// Entry hard_numerator: v0 · ε⁴ · ε⁻¹ (every spelling) with v1·ε⁻¹, the
/// k3j1 point: −117/64. The hard part is integrated and truncated before
/// d = 4−2ε is substituted in it: −2/(d−4) and the trace hide a pole of h_j
/// from the truncation (round 5 before this rule: +117/32 against today's
/// −117/64). They are refused (symbol_in_denominator, unresolved_tensor)
/// into today's value; the explicit spelling is read.
#[test]
fn the_hard_vertex_spellings_agree() -> Result<()> {
    let counts = matrix(
        "v0",
        |spelling| {
            let base = scaled(L3_006, &[("v0", 4), ("v1", -1)]);
            let dot = modify(&base, "v0", "numerator", times(&spelling.text()));
            modify(&dot, "graph", "numerator", times(&spelling.text()))
        },
        Some(q(-117, 64)),
    )?;
    let reasons = counts
        .iter()
        .map(|counts| counts.as_ref().and_then(|counts| counts.reason))
        .collect::<Vec<_>>();
    ensure!(
        reasons
            == [
                None,
                Some("symbol_in_denominator"),
                Some("unresolved_tensor"),
                Some("unresolved_tensor")
            ],
        "{reasons:?}"
    );
    Ok(())
}

/// Entries numerator_prefactor and root_projector: hard-v0 with the graph
/// prefactor ε⁻¹ (explicit, 1/(d−4), metric trace), applied in the root
/// projector. The Dirac-trace spelling is not in this pole comparison: this
/// test's projector (`tensor · prefactor`) leaves a bispinor trace
/// unevaluated, so both routes give 0 there (MATRIX-R45
/// `matrix_r45_numerator_prefactor`, child = today); beta refuses every
/// non-regular prefactor before R*.
#[test]
fn a_numerator_prefactor_pole_is_read_in_every_spelling() -> Result<()> {
    matrix_over(
        "prefactor",
        &SPELLINGS[..3],
        |spelling| modify(&hard_v0(), "graph", "numerator_prefactor", times(&spelling.text())),
        Some(q(-117, 64)),
    )
    .map(|_| ())
}

/// A refusal end to end: v1 · 1/(1 + y) (an unknown symbol in a
/// denominator) with v0·ε³: refused before H, counted, today's value.
#[test]
fn a_symbol_in_a_denominator_refuses_the_correction() -> Result<()> {
    let y = "irh_r5_test::y";
    let factor = format!("1/(1+{y})");
    let dot = modify(&hard_v0(), "v1", "numerator", times(&factor));
    let dot = modify(&dot, "graph", "numerator", times(&factor));
    let (_, counts) = cell("v1 symbol", &dot)?.ok_or_else(|| eyre::eyre!("a value"))?;
    ensure!(
        counts.refused >= 1
            && counts.symbol_in_denominator >= 1
            && counts.reason == Some("symbol_in_denominator"),
        "{counts:?}"
    );
    Ok(())
}

/// Entry projector: the graph projector is read (B) and refused as a graph
/// weight when poleful; R* and beta's R* path do not multiply it.
#[test]
fn a_projector_pole_is_read() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
    for (expression, excess) in [
        (Atom::one(), 0),
        (Atom::one() / &epsilon, 1),
        (Atom::num(-2) / (Atom::var(crate::utils::GS.dim) - Atom::num(4)), 1),
    ] {
        assert_eq!(
            regularity::certified_excess(Factor::GraphProjector, &expression),
            Ok(Some(excess))
        );
        assert_eq!(
            regularity::ensure_regular_weight("projector", &expression).is_ok(),
            excess == 0
        );
    }
    Ok(())
}

/// Entries overall_factor (and the weights at beta's boundary): every
/// spelling of a pole is refused, the regular neighbours and the
/// production factors pass.
#[test]
fn the_graph_weights_are_certified() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let parse = |text: &str| -> Result<Atom> {
        Ok(symbolica::try_parse!(text).map_err(|error| eyre::eyre!("{error}"))?)
    };
    for spelling in SPELLINGS {
        let pole = parse(&spelling.text())?;
        for name in ["overall_factor", "numerator_prefactor", "projector"] {
            let refused = regularity::ensure_regular_weight(name, &pole);
            ensure!(refused.is_err(), "{name} {spelling:?} must be refused");
            eprintln!("IRH_R5_WEIGHT {name} {spelling:?} {}", refused.unwrap_err());
        }
    }
    let neighbour = parse(&tensor().replace("-4)", "-3)"))?;
    for regular in [Atom::one(), Atom::num(1) / Atom::num(2), neighbour] {
        regularity::ensure_regular_weight("overall_factor", &regular)
            .map_err(|error| eyre::eyre!("{error}"))?;
    }
    // The production overall factor, as beta evaluates it.
    let production = parse(
        "(feynkit_generator_factor::{}::AutG(1))^(-1)*feynkit_generator_factor::{}::ExternalFermionOrderingSign(1)",
    )?;
    let evaluated = crate::graph::global::evaluate_overall_factor(production.as_view());
    regularity::ensure_regular_weight("overall_factor", &evaluated)
        .map_err(|error| eyre::eyre!("{error}"))?;
    ensure!(
        regularity::ensure_regular_weight("vertex_weight", &Atom::one()).is_err(),
        "a factor outside the list is refused"
    );
    Ok(())
}

/// Entry propagator_power: a power other than 1 never reaches R* (the
/// runtime graph rejects it), so no ε can enter through a power.
#[test]
fn a_propagator_power_cannot_carry_a_pole() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let dot = L3_006.replacen(
        "v1 -> v3 [id=4, particle_id=1, pdg=21, particle=\"Q\", directed=false, source_slot=1, target_slot=0, power=1,",
        "v1 -> v3 [id=4, particle_id=1, pdg=21, particle=\"Q\", directed=false, source_slot=1, target_slot=0, power=2,",
        1,
    );
    ensure!(dot != L3_006, "the fixture changes edge 4's power");
    use crate::feyngen::feynkit::FeynmanDiagramGammaLoopExt;
    let refused = super::super::recursive_tests::load(&dot)
        .and_then(|diagram| super::super::recursive_tests::colour_factored(&diagram))
        .and_then(|diagram| {
            diagram
                .to_gamma_loop_graph(None, true)
                .map(|_| ())
                .map_err(|error| eyre::eyre!("{error}"))
        });
    let error = format!("{:#}", refused.expect_err("a power of 2 is rejected"));
    ensure!(error.contains("propagator power"), "{error}");
    Ok(())
}

/// Entries compensation, root_normalisation, moment_reduction: the fixed
/// factors have no excess (the compensation polynomial, 1/(d−1), 1/(i q²)
/// with an ε-free invariant, the rank-two vacuum reduction 1/d); a closed
/// trace contracts to d = 4−2ε.
#[test]
fn the_fixed_normalisations_are_regular() -> Result<()> {
    crate::initialisation::test_initialise()?;
    use crate::utils::GS;
    let d = Atom::var(GS.dim);
    for a in [-2_i64, 0, 1, 3] {
        for (loops, order, degree) in [(1, 0, 2), (2, 1, 3), (3, 0, 4)] {
            let factor = super::super::top_mass::compensation_factor(
                super::super::top_mass::TopMassTerms::Keep,
                a,
                loops,
                order,
                degree,
                GS.dim_epsilon,
            );
            assert_eq!(
                regularity::certified_excess(Factor::Compensation, &factor),
                Ok(Some(0)),
                "{factor}"
            );
        }
    }
    let q2 = symbolica::parse!("irh_r5_test::q2");
    for regular in [
        Atom::one() / (&d - Atom::num(1)),
        Atom::one() / &d,
        Atom::one() / (&d * (&d + Atom::num(2))),
    ] {
        assert_eq!(
            regularity::certified_excess(Factor::MomentReduction, &regular),
            Ok(Some(0))
        );
    }
    // 1/(i q²): q² is a symbol in a denominator, refused as a graph factor;
    // beta's normalisation is fixed and builds it from the ε-free external
    // invariant (the entry's certificate is by construction).
    assert!(regularity::certified_excess(Factor::RootNormalisation, &(Atom::one() / &q2)).is_err());
    let trace = super::super::regularity::contract_closed(
        &symbolica::try_parse!(&format!(
            "({METRIC}({MINK}(4,{HEDGE}(81,1)),{MINK}(4,{HEDGE}(82,1))))^2"
        ))
        .map_err(|error| eyre::eyre!("{error}"))?,
    );
    let expected = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
    ensure!(
        trace.as_ref().is_some_and(|trace| (trace - &expected).expand().is_zero()),
        "g_ab g^ab = {trace:?}"
    );
    Ok(())
}

/// Integration with round 4: a regularity refusal is an ordinary fallback
/// under every guard setting with STRICT on: never the fatal depth alarm,
/// never lost by the guarded parent (counted, today's value); a certified
/// raised bound (the hidden tensor) raises no exceeded alarm under STRICT.
#[test]
fn a_refusal_is_a_fallback_under_the_guard_and_strict() -> Result<()> {
    let y = "irh_r5_test::y";
    let factor = format!("1/(1+{y})");
    let symbol = modify(&hard_v0(), "v1", "numerator", times(&factor));
    let symbol = modify(&symbol, "graph", "numerator", times(&factor));
    for (label, dot) in [("symbol", symbol), ("tensor-hidden", tensor_hidden())] {
        let case = Case::new(&dot, None)?;
        let (today, _) = projected(&case, "closed", FIXED)?;
        for guard in [false, true] {
            let cell = format!("{label} guard={guard} strict");
            let test = DepthTest {
                guard: Some(guard),
                strict: Some(true),
                ..FIXED
            };
            let (child, counts) = projected(&case, "child", test)?;
            eprintln!("IRH_R6_STRICT {cell} counts={counts:?}");
            ensure_same(&cell, &child, &today)?;
            if label == "symbol" {
                ensure!(
                    counts.exceeded == 0
                        && counts.refused >= 1
                        && counts.symbol_in_denominator >= 1
                        && counts.reason == Some("symbol_in_denominator"),
                    "{cell}: {counts:?}"
                );
            } else {
                ensure_raised(&cell, &counts)?;
            }
        }
    }
    Ok(())
}

/// U7b: the uncontracted K3 path (`GAMMALOOP_UNCONTRACTED_NUMERATOR`, kind
/// `IrHard`) integrates the remainders of H(γ) only after the round-5
/// certificate, so it honours the same refusals and the same raised bounds:
/// Astra's three −117/64 witnesses and the hard-vertex matrix row, with K3
/// on and in check, the uncontracted strict switch on (a K3 error fails the
/// run instead of falling back).
mod k3 {
    use super::super::super::uncontracted::{
        Kind, Mode, Switch, with_capture, with_strict, with_switch,
    };
    use super::*;

    /// The child route with K3 `mode`: the value, the depth counts and the
    /// number of K3 attempts reported.
    fn child_k3(case: &Case, mode: Mode) -> Result<(Atom, DepthCounts, usize)> {
        let (result, captured) = with_switch(Switch::new(mode, &[Kind::IrHard]), || {
            with_strict(true, || with_capture(|| projected(case, "child", FIXED)))
        });
        let (value, counts) = result?;
        let k3 = captured
            .iter()
            .filter(|captured| captured.kind == Kind::IrHard)
            .count();
        Ok((value, counts, k3))
    }

    #[test]
    fn the_uncontracted_k3_path_reads_the_certified_bound() -> Result<()> {
        for (label, dot) in [
            ("tensor-hidden", tensor_hidden()),
            ("prefactor", hard_pole("numerator_prefactor")),
            ("dim-quotient", dim_quotient()),
        ] {
            let case = Case::new(&dot, None)?;
            let (today, _) = projected(&case, "closed", FIXED)?;
            pinned(label, &today, &[(-1, q(-117, 64))], true)?;
            for mode in [Mode::On, Mode::Check] {
                let cell = format!("{label} K3 {mode:?}");
                let (child, counts, k3) = child_k3(&case, mode)?;
                eprintln!("IRH_R6_K3 {cell} k3_attempts={k3} counts={counts:?}");
                ensure_same(&cell, &child, &today)?;
                ensure_raised(&cell, &counts)?;
                ensure!(k3 >= 1, "{cell}: the K3 path was not exercised");
            }
        }
        Ok(())
    }

    #[test]
    fn the_uncontracted_k3_path_honours_the_refusals() -> Result<()> {
        let expected = [
            None,
            Some("symbol_in_denominator"),
            Some("unresolved_tensor"),
            Some("unresolved_tensor"),
        ];
        for (spelling, reason) in SPELLINGS.into_iter().zip(expected) {
            let base = scaled(L3_006, &[("v0", 4), ("v1", -1)]);
            let dot = modify(&base, "v0", "numerator", times(&spelling.text()));
            let dot = modify(&dot, "graph", "numerator", times(&spelling.text()));
            let case = Case::new(&dot, None)?;
            let (today, _) = projected(&case, "closed", FIXED)?;
            for mode in [Mode::On, Mode::Check] {
                let cell = format!("v0 {spelling:?} K3 {mode:?}");
                let (child, counts, k3) = child_k3(&case, mode)?;
                eprintln!("IRH_R6_K3 {cell} k3_attempts={k3} counts={counts:?}");
                ensure_same(&cell, &child, &today)?;
                ensure!(counts.exceeded == 0 && counts.reason == reason, "{cell}: {counts:?}");
                if reason.is_some() {
                    ensure!(counts.refused >= 1, "{cell}: {counts:?}");
                }
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// MATRIX-R45: the mechanical variant matrix (coordinator requirement, 0
// WRONG). Every attribute that enters the product × ε⁻¹ in every spelling
// ([`SPELLINGS`]) × the child route under the production planner, STRICT and
// the guard, against today's (generic) path. A cell is RIGHT (the child route
// itself succeeded, no refusal, and equals today's value exactly), REFUSED
// (loud and counted, today's value; the class is pinned), LOUD (both routes
// fail with an error, no value) or WRONG (a value different from today's: the
// test fails). The native half (beta binaries, the DOT attributes) is
// testdata/irh_matrix/gen_matrix.py with its expected table
// testdata/irh_matrix/expected.tsv.
// ---------------------------------------------------------------------------

/// A MATRIX-R45 disposition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Disposition {
    Right,
    Refused(&'static str),
    Loud,
}

use Disposition::{Loud, Refused, Right};

const SID: Disposition = Refused("symbol_in_denominator");
const UT: Disposition = Refused("unresolved_tensor");

/// The three child-route settings of every cell.
const MODES: [(&str, DepthTest); 3] = [
    ("child", FIXED),
    (
        "strict",
        DepthTest {
            strict: Some(true),
            ..FIXED
        },
    ),
    (
        "guard",
        DepthTest {
            guard: Some(true),
            ..FIXED
        },
    ),
];

/// Classify one child run against today's; a WRONG value, a child-only
/// failure or a child route that is not taken fail the test.
fn classify(
    label: &str,
    today: &Result<Atom>,
    child: Result<(Atom, DepthCounts)>,
    nested: bool,
) -> Result<Disposition> {
    let (today, child, counts) = match (today, child) {
        (Err(today), Err(child)) => {
            eprintln!("IRH_R45 {label} LOUD today={today:#} child={child:#}");
            return Ok(Loud);
        }
        (Err(today), Ok((child, _))) => {
            eyre::bail!("{label}: today's path fails ({today:#}) but the child route returns {child}")
        }
        (Ok(_), Err(child)) => eyre::bail!("{label}: the child route alone fails: {child:#}"),
        (Ok(today), Ok((child, counts))) => (today, child, counts),
    };
    if !super::same(&child, today)? {
        eyre::bail!(
            "{label}: WRONG child_pole={} today_pole={} counts={counts:?}",
            super::laurent(&super::record(&child)?, -1)?,
            super::laurent(&super::record(today)?, -1)?,
        );
    }
    ensure!(counts.exceeded == 0, "{label}: exceeded {counts:?}");
    let disposition = if counts.refused >= 1 {
        Refused(counts.reason.unwrap_or("unknown"))
    } else {
        ensure!(
            counts.corrections >= 1 && (!nested || counts.nested >= 1),
            "{label}: the child route was not taken (nested={nested}): {counts:?}"
        );
        Right
    };
    eprintln!(
        "IRH_R45 {label} {disposition:?} today_pole={} counts={counts:?}",
        super::laurent(&super::record(today)?, -1)?
    );
    Ok(disposition)
}

/// One attribute of the root (projected as beta): every spelling under every
/// mode, against `expected` (one disposition per spelling, the same under
/// every mode).
fn r45_root(label: &str, fixture: impl Fn(Spelling) -> String, expected: [Disposition; 4]) -> Result<()> {
    let mut wrong = Vec::new();
    for (spelling, expected) in SPELLINGS.into_iter().zip(expected) {
        let case = Case::new(&fixture(spelling), None)?;
        let today = projected(&case, "closed", FIXED).map(|(value, _)| value);
        for (mode, test) in MODES {
            let cell = format!("{label} {spelling:?} {mode}");
            match classify(&cell, &today, projected(&case, "child", test), false) {
                Ok(found) if found == expected => {}
                Ok(found) => wrong.push(format!("{cell}: {found:?}, expected {expected:?}")),
                Err(error) => wrong.push(format!("{error:#}")),
            }
        }
    }
    ensure!(wrong.is_empty(), "{label}: {wrong:#?}");
    Ok(())
}

/// Entry quotient_vertex_numerator on the base hard-v0 (v0·ε³).
#[test]
fn matrix_r45_quotient_vertex() -> Result<()> {
    r45_root(
        "v1",
        |s| {
            let dot = modify(&hard_v0(), "v1", "numerator", times(&s.text()));
            modify(&dot, "graph", "numerator", times(&s.text()))
        },
        [Right, Right, Right, Right],
    )
}

/// Entry quotient_edge_numerator: the quotient edge 4 (v1 → v3).
#[test]
fn matrix_r45_quotient_edge() -> Result<()> {
    r45_root(
        "e4",
        |s| {
            let dot = modify(&hard_v0(), "v1 ->", "numerator", times(&s.text()));
            modify(&dot, "graph", "numerator", times(&s.text()))
        },
        [Right, Right, Right, Right],
    )
}

/// Entry hard_numerator: the hard vertex v0 on the base (v0·ε³·S).
#[test]
fn matrix_r45_hard_vertex() -> Result<()> {
    r45_root(
        "v0",
        |s| {
            let dot = modify(&hard_v0(), "v0", "numerator", times(&s.text()));
            modify(&dot, "graph", "numerator", times(&s.text()))
        },
        [Right, SID, UT, UT],
    )
}

/// Entry hard_numerator: the hard edge 7 (v2 → v4) on the base.
#[test]
fn matrix_r45_hard_edge() -> Result<()> {
    r45_root(
        "e7",
        |s| {
            let dot = modify(&hard_v0(), "v2 ->", "numerator", times(&s.text()));
            modify(&dot, "graph", "numerator", times(&s.text()))
        },
        [Right, SID, UT, UT],
    )
}

/// Entry hard_numerator with every row hidden: v0·ε⁴·S, v1·ε⁻¹ (k3j1).
#[test]
fn matrix_r45_hard_vertex_hidden_rows() -> Result<()> {
    r45_root(
        "v0k3j1",
        |s| {
            let base = scaled(L3_006, &[("v0", 4), ("v1", -1)]);
            let dot = modify(&base, "v0", "numerator", times(&s.text()));
            modify(&dot, "graph", "numerator", times(&s.text()))
        },
        [Right, SID, UT, UT],
    )
}

/// Entries numerator_prefactor, root_projector (beta's callback multiplies
/// σ(prefactor)) and root_normalisation: the graph prefactor.
#[test]
fn matrix_r45_numerator_prefactor() -> Result<()> {
    r45_root(
        "pref",
        |s| modify(&hard_v0(), "graph", "numerator_prefactor", times(&s.text())),
        // The trace: vakint's conversion fails on both routes ('¿'); the
        // Dirac trace: this test's projector leaves the bispinor trace
        // unevaluated, 0 on both routes (beta refuses both before R*).
        [Right, Right, Loud, Right],
    )
}

/// Entry projector: the graph projector (R* does not multiply it; read in
/// B).
#[test]
fn matrix_r45_projector() -> Result<()> {
    r45_root(
        "proj",
        |s| modify(&hard_v0(), "graph", "projector", times(&s.text())),
        [Right, Right, Right, Right],
    )
}

/// The nested cells: c2-m9 (k2j1 with the marked line 9; its only
/// child-route condition is in the child calls nested in the root, the R*
/// of its decorated quotients, IRH-DEPTH.md 2.6) with v1's ε⁻¹ in every
/// spelling, on the unprojected path as the nested tests (round 2-4). RIGHT
/// requires nested child calls (`nested ≥ 1`); the refusals are counted in
/// the root's collect scope, the nested evaluations included.
#[test]
fn matrix_r45_nested_decorated_quotient() -> Result<()> {
    let mut wrong = Vec::new();
    for spelling in SPELLINGS {
        let base = scaled(L3_006, &[("v0", 2)]);
        let dot = modify(&base, "v1", "numerator", times(&spelling.text()));
        let dot = modify(&dot, "graph", "numerator", times(&spelling.text()));
        let case = Case::new(&dot, Some(9))?;
        let today = case.today();
        if matches!(spelling, Spelling::Epsilon) {
            pinned("c2-m9 r45", today.as_ref().map_err(|e| eyre::eyre!("{e:#}"))?, &[(-2, q(-117, 64))], false)?;
        }
        // v1 lies in the hard part of some correction of the nested
        // evaluations: a 1/(d−4) or a trace there is refused (hard_numerator).
        let expected = match spelling {
            Spelling::Epsilon => Right,
            Spelling::DimensionMinusFour => SID,
            Spelling::Trace | Spelling::DiracTrace => UT,
        };
        for (mode, test) in MODES {
            let cell = format!("nested-v1 {spelling:?} {mode}");
            let child = case.child(run(test)).map(|run| (run.value, run.root));
            match classify(&cell, &today, child, true) {
                Ok(found) if found == expected => {}
                Ok(found) => wrong.push(format!("{cell}: {found:?}, expected {expected:?}")),
                Err(error) => wrong.push(format!("{error:#}")),
            }
        }
    }
    ensure!(wrong.is_empty(), "{wrong:#?}");
    Ok(())
}

/// Every entry of the closed list × every spelling, at the certificate: the
/// normalised factor shows the pole (excess ≥ 1) or is refused; never read
/// as regular. The graph weights are refused at beta's boundary in every
/// spelling.
#[test]
fn matrix_r45_every_listed_factor_reads_every_spelling() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let x = symbolica::parse!("irh_r45_test::x");
    let mut wrong = Vec::new();
    for entry in regularity::FACTORS {
        for spelling in SPELLINGS {
            let pole = symbolica::try_parse!(&spelling.text()).map_err(|error| eyre::eyre!("{error}"))?;
            let outcome = regularity::certified_excess(entry.factor, &(&x * &pole));
            eprintln!("IRH_R45_UNIT {} {spelling:?} {outcome:?}", entry.name);
            match outcome {
                Ok(Some(n)) if n >= 1 => {}
                Err(_) => {}
                other => wrong.push(format!("{} {spelling:?}: {other:?}", entry.name)),
            }
        }
    }
    for spelling in SPELLINGS {
        let pole = symbolica::try_parse!(&spelling.text()).map_err(|error| eyre::eyre!("{error}"))?;
        for name in ["overall_factor", "numerator_prefactor", "projector"] {
            if regularity::ensure_regular_weight(name, &pole).is_ok() {
                wrong.push(format!("weight {name} {spelling:?} accepted"));
            }
        }
    }
    // The Dirac trace contracts to 4d.
    let dirac = symbolica::try_parse!(&dirac_product()).map_err(|error| eyre::eyre!("{error}"))?;
    let contracted = regularity::contract_closed(&dirac);
    let expected = Atom::num(16) - Atom::num(8) * Atom::var(crate::utils::GS.dim_epsilon);
    ensure!(
        contracted.as_ref().is_some_and(|trace| (trace - &expected).expand().is_zero()),
        "tr(γ_a γ^a) = {contracted:?}"
    );
    ensure!(wrong.is_empty(), "{wrong:#?}");
    Ok(())
}

/// The root projector callback (entry root_projector) carrying the pole
/// itself. The callback is opaque to R*: it is certified by CONTRACT
/// (`RStarPlan::projected_counterterm`: regular in ε), not at run time, and
/// beta's callback multiplies only σ(numerator_prefactor), refused by beta
/// when poleful, and fixed regular normalisations. A caller that breaks the
/// contract gets a wrong child value in the explicit and 1/(d−4) spellings
/// (+117/32 against today's −117/64, also under STRICT and the guard: no
/// row measures the callback). This test PINS that known gap of the public
/// API (MATRIX-R45 open item), so a change of it is noticed; it is not
/// reachable from beta.
#[test]
fn matrix_r45_a_poleful_root_projector_callback_is_outside_the_contract() -> Result<()> {
    let case = Case::new(&hard_v0(), None)?;
    let mut unexpected = Vec::new();
    for spelling in SPELLINGS {
        let pole = symbolica::try_parse!(&spelling.text()).map_err(|error| eyre::eyre!("{error}"))?;
        let run = |route: &'static str, test: DepthTest| -> Result<(Atom, DepthCounts)> {
            let pole = pole.clone();
            let projector = move |tensor: &Atom| -> Result<Atom> { Ok(tensor * &pole) };
            let (value, counts) = depth::with_depth_test(test, || {
                depth::collect(|| {
                    with_ir_hard(route, || {
                        case.plan.projected_counterterm(
                            &case.graph,
                            case.model(),
                            &case.vakint,
                            &case.settings,
                            &projector,
                        )
                    })
                })
            });
            Ok((value?, counts))
        };
        let today = run("closed", FIXED).map(|(value, _)| value);
        for (mode, test) in MODES {
            let cell = format!("callback {spelling:?} {mode}");
            let found = classify(&cell, &today, run("child", test), false);
            eprintln!("IRH_R45_CALLBACK {cell} {:?}", found.as_ref().map_err(|e| format!("{e:#}")));
            let pinned = match (spelling, &found) {
                (Spelling::Epsilon | Spelling::DimensionMinusFour, Err(error)) => {
                    let error = format!("{error:#}");
                    error.contains("WRONG child_pole=117/32 today_pole=-117/64")
                }
                (Spelling::Trace, Ok(Loud)) => true,
                (Spelling::DiracTrace, Ok(Right)) => true,
                _ => false,
            };
            if !pinned {
                unexpected.push(format!("{cell}: {:?}", found.map_err(|e| format!("{e:#}"))));
            }
        }
    }
    ensure!(unexpected.is_empty(), "the pinned callback gap changed: {unexpected:#?}");
    Ok(())
}

/// The DOT overall factor of gen_matrix.py's `overall-<spelling>` variants:
/// the production factor times the spelling (and its regular neighbour).
fn overall_variants() -> Vec<(String, String)> {
    let production = "(feynkit_generator_factor::{}::AutG(1))^(-1)*feynkit_generator_factor::{}::ExternalFermionOrderingSign(1)";
    let mut variants = vec![
        ("eps", Spelling::Epsilon.text()),
        ("dm4", Spelling::DimensionMinusFour.text()),
        ("gtr", Spelling::Trace.text()),
        ("dtr", Spelling::DiracTrace.text()),
        ("gtr3", tensor().replace("-4)", "-3)")),
    ]
    .into_iter()
    .map(|(name, factor)| (format!("overall-{name}"), format!("({production})*({factor})")))
    .collect::<Vec<_>>();
    variants.push(("overall-none".into(), production.into()));
    variants
}

/// gen_matrix.py writes the manifest's overall_factor binding from
/// testdata/irh_matrix/overall_canonical.tsv; it must be the canonical
/// string beta compares with (`diagram.overall_factor()`, the parsed DOT
/// attribute).
#[test]
fn matrix_r45_the_overall_canonical_table_is_current() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let table = include_str!("testdata/irh_matrix/overall_canonical.tsv");
    let mut wrong = Vec::new();
    for (name, text) in overall_variants() {
        let canonical = symbolica::try_parse!(&text)
            .map_err(|error| eyre::eyre!("{error}"))?
            .to_canonical_string();
        let line = format!("{name}\t{canonical}");
        if !table.lines().any(|row| row == line) {
            wrong.push(line);
        }
    }
    ensure!(
        wrong.is_empty(),
        "IRH_R45_CANONICAL_TABLE missing or stale rows:\n{}",
        wrong.join("\n")
    );
    Ok(())
}
