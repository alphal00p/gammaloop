//! IRH-DEPTH round 8b: the reviewer's NO GO on rounds 4-7
//! (review-20261122/IRH-DEPTH/r47) and the topological bound
//! (results/consolidation/irh-depth/TOPOLOGICAL-BOUND.md), on
//! B_self_energy_L3_006 as the earlier rounds ([`super`]).
//!
//! * C1 (finding 1 and the quotient refutation): a scalar `d` is
//!   `4 - 2ε` before ANY expansion in `ε`. A factor equal to −1 at
//!   `d = 4 - 2ε`, spelled `1/(1+(d−4)/ε)` or `ε/(ε+d−4)`, gives the value of
//!   the literal −1 on every route (child with STRICT off and on, the closed
//!   path, the projected path), on the hard vertex, a quotient vertex, a
//!   quotient edge, a hard edge and inside a propagator insertion.
//! * Finding 4: the graph's stored projector is not applied by R*, so 0, ε
//!   and ε² there move neither the bound nor the value; a regular numerator
//!   prefactor ε is no STRICT alarm.
//! * R1: an opaque root projector callback is never probed; the child route
//!   declines (counted, `opaque_projector`).
//! * Finding 5: an alarm a callback raises is never swallowed.

use eyre::{Result, ensure, eyre};
use symbolica::atom::{Atom, AtomCore};

use super::super::{
    ClosedRootNormalization,
    closed::with_ir_hard,
    depth::{self, DepthCounts, DepthTest},
    expression::RStarIntegrand,
};
use super::{Case, EPS, FIXED, L3_006, k, q, same};
use crate::{graph::Graph, utils::GS, uv::approx::integrated::Integrated};
use linnet::half_edge::subgraph::SuBitGraph;

const DIM: &str = "gammalooprs::{}::dim";

/// The three spellings of −1 at `d = 4 − 2ε` (the review's analytic inputs):
/// the literal, `1/(1+(d−4)/ε)` (the cross form, ε/(d−4) at fixed d) and
/// `ε/(ε+d−4)` (the quotient form).
fn spellings() -> [(&'static str, String); 3] {
    [
        ("literal", "-1".to_owned()),
        ("cross", format!("1/(1+({DIM}-4)/{EPS})")),
        ("form", format!("{EPS}/({EPS}+{DIM}-4)")),
    ]
}

/// The span of the first `field="…"` of `line` not preceded by a lowercase
/// letter or `_`.
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

/// `dot` with `field` of the first line starting with `  {term} ` multiplied
/// by `factor`.
fn times(dot: &str, term: &str, field: &str, factor: &str) -> String {
    let mut count = 0;
    let out = dot
        .split_inclusive('\n')
        .map(|line| {
            if line.starts_with(&format!("  {term} ")) && count == 0 {
                let (body, end) = field_span(line, field).expect("the field");
                count += 1;
                format!(
                    "{}({})*({factor}){}",
                    &line[..body],
                    &line[body..end],
                    &line[end..]
                )
            } else {
                line.to_owned()
            }
        })
        .collect();
    assert_eq!(count, 1, "{term} {field}");
    out
}

/// `base` with the numerator of `term` (a vertex `v0`, or an edge `v1 ->`)
/// and the aggregate graph numerator multiplied by `factor`, as the
/// review's inputs.
fn spelled(base: &str, term: &str, factor: &str) -> String {
    let dot = times(base, term, "numerator", factor);
    times(&dot, "graph", "numerator", factor)
}

/// One run of the unweighted public API (`counterterm_with_children`, no
/// root projector) on `route` (`child`, `closed`, `off`), STRICT `strict`.
fn unweighted(case: &Case, route: &'static str, strict: bool) -> Result<(Atom, DepthCounts)> {
    let test = DepthTest {
        strict: Some(strict),
        ..FIXED
    };
    let (value, counts) = depth::with_depth_test(test, || {
        depth::collect(|| {
            with_ir_hard(route, || {
                case.plan.counterterm_with_children(
                    &case.graph,
                    case.model(),
                    &case.vakint,
                    &case.settings,
                    case.children.clone(),
                )
            })
        })
    });
    eprintln!("IRH_R8B_UNWEIGHTED route={route} strict={strict} counts={counts:?}");
    Ok((value?, counts))
}

/// One run through an opaque root projector callback (the weighted public
/// API, `projected_counterterm`): `projector` applied to every complete outer
/// Taylor polynomial.
fn opaque(
    case: &Case,
    route: &'static str,
    strict: bool,
    projector: &dyn Fn(&Atom) -> Result<Atom>,
) -> Result<(Atom, DepthCounts)> {
    let test = DepthTest {
        strict: Some(strict),
        ..FIXED
    };
    let (value, counts) = depth::with_depth_test(test, || {
        depth::collect(|| {
            with_ir_hard(route, || {
                case.plan.projected_counterterm(
                    &case.graph,
                    case.model(),
                    &case.vakint,
                    &case.settings,
                    projector,
                )
            })
        })
    });
    eprintln!("IRH_R8B_OPAQUE route={route} strict={strict} counts={counts:?}");
    Ok((value?, counts))
}

/// The closed metric root with the graph's numerator prefactor as the
/// declared scalar normalisation (`c = σ(prefactor)`), as beta's typed
/// projector: the callback is `c · σ(X)` on closed scalars, and the same data
/// is passed as [`ClosedRootNormalization`]
/// (`projected_counterterm_closing_root_scalar`).
fn typed_root(graph: &Graph) -> Result<ClosedRootNormalization> {
    let (external, external_routing) =
        RStarIntegrand::physical_external_routing(graph, &graph.external_filter::<SuBitGraph>())?;
    Ok(ClosedRootNormalization {
        external,
        factor: graph
            .global_prefactor
            .num
            .replace_multiple(&external_routing),
        external_routing,
    })
}

/// One run of the closed metric root through the typed normalisation.
fn typed(case: &Case, route: &'static str, strict: bool) -> Result<(Atom, DepthCounts)> {
    let closed_root = typed_root(&case.graph)?;
    let factor = closed_root.factor.clone();
    let routing = closed_root.external_routing.clone();
    let callback =
        move |tensor: &Atom| -> Result<Atom> { Ok(&factor * tensor.replace_multiple(&routing)) };
    let test = DepthTest {
        strict: Some(strict),
        ..FIXED
    };
    // The closed root linearises momenta factor by factor: the streamed and
    // the whole maps may differ in spelling only (compared algebraically).
    let ((value, counts), _) = crate::uv::approx::streamed::with_algebraic_check(|| {
        depth::with_depth_test(test, || {
            depth::collect(|| {
                with_ir_hard(route, || {
                    case.plan.projected_counterterm_closing_root_scalar(
                        &case.graph,
                        case.model(),
                        &case.vakint,
                        &case.settings,
                        &callback,
                        &closed_root,
                    )
                })
            })
        })
    });
    eprintln!("IRH_R8B_TYPED route={route} strict={strict} counts={counts:?}");
    Ok((value?, counts))
}

/// The pole parts of two closed (scalar) values agree at `d = 4 − 2ε`.
fn same_scalar(a: &Atom, b: &Atom) -> Result<bool> {
    let difference = Integrated::dimensionally_regularized(&(a - b));
    Ok(difference
        .series(GS.dim_epsilon, 0, -1)
        .map_err(|error| eyre!("{error}"))?
        .to_atom()
        .expand()
        .together()
        .is_zero())
}

/// The ε⁻¹..ε⁻⁴ part of a value in the record units, for the messages.
fn label_poles(value: &Atom) -> String {
    super::record(value)
        .and_then(|record| {
            Ok(record
                .series(GS.dim_epsilon, 0, -1)
                .map_err(|error| eyre!("{error}"))?
                .to_atom()
                .expand())
        })
        .map_or_else(
            |error| format!("unreadable ({error})"),
            |poles| poles.to_string(),
        )
}

/// The spelling matrix of one position (`term` of `base`): for every
/// spelling of −1 and every route (child STRICT off and on, closed, off)
/// the unweighted counterterm equals the literal's on the projected path
/// (`off`); the literal's equals −1 times `base`'s on that path.
/// `taken`: the child route takes the correction (every declared order is
/// 0: the spellings normalise to −1); otherwise it declines it
/// (`declared_order`, the k3j1 carrier's v1·ε⁻¹).
fn unweighted_matrix(label: &str, base: &str, term: &str, taken: bool) -> Result<()> {
    let base_case = Case::new(base, None)?;
    let (base_value, _) = unweighted(&base_case, "off", false)?;
    drop(base_case);
    let mut reference = None::<Atom>;
    let mut failures = Vec::new();
    for (name, factor) in spellings() {
        let case = Case::new(&spelled(base, term, &factor), None)?;
        for (route, strict) in [
            ("off", false),
            ("closed", false),
            ("child", false),
            ("child", true),
        ] {
            let cell = format!("{label} {name} {route} strict={strict}");
            let value = match unweighted(&case, route, strict) {
                Ok((value, counts)) => {
                    if route == "child" {
                        // The quotient form ε/(ε+d−4) on a HARD factor keeps
                        // round 5's conservative refusal (a d-denominator of
                        // the hard part that vanishes at d = 4), counted; its
                        // value is checked like every other cell's.
                        let hard_form = name == "form" && (term == "v0" || term == "v2 ->");
                        let disposition = if taken {
                            (counts.taken >= 1 && counts.declined == 0 && counts.exceeded == 0)
                                || (hard_form
                                    && counts.taken == 0
                                    && counts.reason == Some("symbol_in_denominator"))
                        } else {
                            counts.taken == 0
                                && ((counts.declined >= 1
                                    && counts.reason == Some("declared_order"))
                                    || (hard_form
                                        && counts.reason == Some("symbol_in_denominator")))
                        };
                        if !disposition {
                            failures.push(format!("{cell}: disposition {counts:?}"));
                        }
                    }
                    value
                }
                Err(error) => {
                    failures.push(format!("{cell}: error {error:#}"));
                    continue;
                }
            };
            eprintln!("IRH_R8B_CELL {cell} poles={}", label_poles(&value));
            match &reference {
                None => {
                    ensure!(
                        name == "literal" && route == "off",
                        "the reference is the literal on the projected path"
                    );
                    ensure!(
                        same(&value, &(-&base_value))?,
                        "{cell}: the literal −1 is not −1 times the base: {} against {}",
                        label_poles(&value),
                        label_poles(&base_value)
                    );
                    reference = Some(value);
                }
                Some(reference) => {
                    if !same(&value, reference)? {
                        failures.push(format!(
                            "{cell}: {} against the literal's {}",
                            label_poles(&value),
                            label_poles(reference)
                        ));
                    }
                }
            }
        }
    }
    ensure!(failures.is_empty(), "{label}: {failures:#?}");
    Ok(())
}

/// (a) Finding 1, the hard numerator (vertex v0): the child route took the
/// correction and gave −117/32 against +117/64 on the review's k3j1 input
/// (`1/(1+(d−4)/ε)` at fixed d is ε/(d−4) − ε²/(d−4)² + …, every term of
/// order 0 at d = 4 − 2ε: no finite truncation of it is the value).
#[test]
fn round8b_a_hard_vertex_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("hard v0", L3_006, "v0", true)
}

/// (a) the review's own carrier: k3j1 (v0·ε³, v1·ε⁻¹).
#[test]
fn round8b_a_hard_vertex_spellings_on_k3j1() -> Result<()> {
    unweighted_matrix("hard v0 k3j1", &k(3, 1), "v0", false)
}

/// (b) The refutation (topo/refute): the same spellings on the quotient
/// vertex v1 gave 117/256 on the child AND on the generic route.
#[test]
fn round8b_b_quotient_vertex_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("quotient v1", L3_006, "v1", true)
}

/// (b) the refutation's carrier: k3j1.
#[test]
fn round8b_b_quotient_vertex_spellings_on_k3j1() -> Result<()> {
    unweighted_matrix("quotient v1 k3j1", &k(3, 1), "v1", false)
}

/// (c) A quotient edge (v1 → v3, edge 4) and a hard edge (v2 → v3).
#[test]
fn round8b_c_quotient_edge_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("quotient edge v1->v3", L3_006, "v1 ->", true)
}

#[test]
fn round8b_c_hard_edge_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("hard edge v2->v3", L3_006, "v2 ->", true)
}

/// (c)/(a) through the typed closed root (beta's production form): the
/// spellings on the hard vertex and on the quotient vertex equal the literal
/// on every route.
#[test]
fn round8b_c_typed_root_spellings_equal_the_literal() -> Result<()> {
    let mut failures = Vec::new();
    for term in ["v0", "v1"] {
        let mut reference = None::<Atom>;
        for (name, factor) in spellings() {
            let case = Case::new(&spelled(L3_006, term, &factor), None)?;
            for (route, strict) in [
                ("off", false),
                ("closed", false),
                ("child", false),
                ("child", true),
            ] {
                let cell = format!("typed {term} {name} {route} strict={strict}");
                let value = match typed(&case, route, strict) {
                    Ok((value, counts)) => {
                        eprintln!("IRH_R8B_CELL {cell} counts={counts:?}");
                        let hard_form = name == "form" && term == "v0";
                        if route == "child"
                            && !(counts.taken >= 1 && counts.declined == 0 && counts.exceeded == 0)
                            && !(hard_form && counts.reason == Some("symbol_in_denominator"))
                        {
                            failures
                                .push(format!("{cell}: the typed root must be taken: {counts:?}"));
                        }
                        value
                    }
                    Err(error) => {
                        failures.push(format!("{cell}: error {error:#}"));
                        continue;
                    }
                };
                match &reference {
                    None => reference = Some(value),
                    Some(reference) => {
                        if !same_scalar(&value, reference)? {
                            failures.push(format!("{cell}: differs from the literal on off"));
                        }
                    }
                }
            }
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// (d) Finding 4: the graph's stored projector is not applied by R* (nor by
/// beta's R* path), so 0, ε and ε² there must not move B or the value, on
/// the unweighted and on the weighted (projected) public API. On 969447edcf
/// a zero projector made the quotient certificate 0 (B = 1) and the child
/// route took k3j1 at depth 1: +117/32 against −117/64.
#[test]
fn round8b_d_the_stored_projector_moves_nothing() -> Result<()> {
    let case = Case::new(&k(3, 1), None)?;
    let (unit, unit_counts) = unweighted(&case, "child", false)?;
    let today = case.today()?;
    ensure!(same(&unit, &today)?, "k3j1 with projector 1");
    let prefactor = case.graph.global_prefactor.num.clone();
    let weighted = move |tensor: &Atom| -> Result<Atom> { Ok(tensor * &prefactor) };
    let (weighted_unit, weighted_counts) = opaque(&case, "child", false, &weighted)?;
    let (weighted_today, _) = opaque(&case, "closed", false, &weighted)?;
    ensure!(
        same(&weighted_unit, &weighted_today)?,
        "k3j1 weighted with projector 1"
    );
    let epsilon = Atom::var(GS.dim_epsilon);
    let mut failures = Vec::new();
    for (name, projector) in [
        ("0", Atom::Zero),
        ("eps", epsilon.clone()),
        ("eps^2", epsilon.pow(2)),
    ] {
        let mut moved = Case::new(&k(3, 1), None)?;
        moved.graph.global_prefactor.projector = projector;
        let (value, counts) = unweighted(&moved, "child", false)?;
        if !same(&value, &today)? || counts.bounded != unit_counts.bounded {
            failures.push(format!(
                "unweighted projector={name}: {} (bounded {}) against {} (bounded {})",
                label_poles(&value),
                counts.bounded,
                label_poles(&today),
                unit_counts.bounded
            ));
        }
        let (value, counts) = opaque(&moved, "child", false, &weighted)?;
        if !same(&value, &weighted_today)? || counts.bounded != weighted_counts.bounded {
            failures.push(format!(
                "weighted projector={name}: {} (bounded {}) against {} (bounded {})",
                label_poles(&value),
                counts.bounded,
                label_poles(&weighted_today),
                weighted_counts.bounded
            ));
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// (d) The values finding: k2j1 (v0·ε², v1·ε⁻¹) with the regular numerator
/// prefactor ε, applied by the root projector as beta does. The prefactor's
/// order must not offset the quotient pole (the old certificate multiplied
/// them: bound 1, a measured need 2, a false STRICT alarm on both release
/// binaries). Child with STRICT on and off equals today's value.
#[test]
fn round8b_d_a_regular_prefactor_is_no_strict_alarm() -> Result<()> {
    let dot = times(&k(2, 1), "graph", "numerator_prefactor", EPS);
    let case = Case::new(&dot, None)?;
    let prefactor = case.graph.global_prefactor.num.clone();
    let weighted = move |tensor: &Atom| -> Result<Atom> { Ok(tensor * &prefactor) };
    let (today, _) = opaque(&case, "closed", false, &weighted)?;
    super::pinned("k2j1 prefactor eps", &today, &[(-1, q(-117, 64))], true)?;
    for strict in [false, true] {
        let (child, counts) = opaque(&case, "child", strict, &weighted)
            .map_err(|error| eyre!("strict={strict}: a false alarm: {error:#}"))?;
        ensure!(
            same(&child, &today)?,
            "strict={strict}: {} against {}",
            label_poles(&child),
            label_poles(&today)
        );
        ensure!(counts.exceeded == 0, "strict={strict}: {counts:?}");
        let (typed_child, typed_counts) = typed(&case, "child", strict)
            .map_err(|error| eyre!("typed strict={strict}: a false alarm: {error:#}"))?;
        let (typed_today, _) = typed(&case, "closed", false)?;
        ensure!(
            same_scalar(&typed_child, &typed_today)?,
            "typed strict={strict}: the child differs from today's"
        );
        ensure!(
            typed_counts.exceeded == 0,
            "typed strict={strict}: {typed_counts:?}"
        );
    }
    Ok(())
}

/// The symbol the old unit-insertion probe renamed ε to: an input that holds
/// it reached the callback from the probe.
fn probe_symbol() -> symbolica::atom::Symbol {
    symbolica::symbol!("irh_callback_probe_epsilon")
}

/// (e) R1: an opaque root projector callback is not a certificate; the child
/// route declines to the generic path, loudly and counted
/// (`opaque_projector`), and never probes it. On L3_006 (every row shown)
/// and on k3j1 (every group vanishes at the old first depth, where the old
/// planner probed the callback on the unit insertion).
#[test]
fn round8b_e_an_opaque_callback_declines_without_a_probe() -> Result<()> {
    for (label, dot) in [("k0j0", L3_006.to_owned()), ("k3j1", k(3, 1))] {
        let case = Case::new(&dot, None)?;
        let probed = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let seen = probed.clone();
        let identity = move |tensor: &Atom| -> Result<Atom> {
            if tensor.contains_symbol(probe_symbol()) {
                seen.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            }
            Ok(tensor.clone())
        };
        let (today, _) = opaque(&case, "closed", false, &identity)?;
        let (child, counts) = opaque(&case, "child", false, &identity)?;
        ensure!(same(&child, &today)?, "{label}: the value");
        ensure!(
            counts.taken == 0 && counts.refused >= 1 && counts.reason == Some("opaque_projector"),
            "{label}: an opaque callback must decline (counted, opaque_projector): {counts:?}"
        );
        let probes = probed.load(std::sync::atomic::Ordering::Relaxed);
        ensure!(
            probes == 0,
            "{label}: the callback was probed {probes} time(s)"
        );
    }
    Ok(())
}

/// (f) Finding 5: every alarm a root projector callback raises is fatal on
/// every route, never swallowed into a decline or a value. (Repair RC-S3:
/// the half with a callback forwarding an alarm only on the old probe's
/// input could not fail once the probe was removed, and is dropped.)
#[test]
fn round8b_f_a_callback_alarm_is_never_swallowed() -> Result<()> {
    let case = Case::new(&k(3, 1), None)?;
    // A callback that always raises the alarm: fatal on every route.
    let always = |_: &Atom| -> Result<Atom> {
        Err(depth::alarm(
            "IRH-DEPTH: an alarm raised by the root projector".to_owned(),
        ))
    };
    for route in ["closed", "child", "off"] {
        match opaque(&case, route, false, &always) {
            Ok(_) => eyre::bail!("{route}: the alarm was swallowed"),
            Err(error) => ensure!(depth::is_alarm(&error), "{route}: {error:#}"),
        }
    }
    Ok(())
}

/// C1 by construction: the input of every ε expansion
/// (`Integrated::regulated_for_expansion`) has no scalar d (the cross form
/// becomes −1 exactly), a Lorentz `mink(d, i)` keeps its dimension, and an
/// expression without a scalar d is not copied.
#[test]
fn round8b_c1_the_expansion_input_is_dimension_free() -> Result<()> {
    crate::initialisation::test_initialise()?;
    use spenso::structure::representation::{Minkowski, RepName};
    let d = Atom::var(GS.dim);
    let epsilon = Atom::var(GS.dim_epsilon);
    let cross = Atom::one() / (Atom::one() + (&d - Atom::num(4)) / &epsilon);
    let regulated = Integrated::regulated_for_expansion(&cross, "test")?;
    ensure!(
        (&*regulated + Atom::one()).expand().together().is_zero(),
        "{regulated}"
    );
    ensure!(crate::uv::approx::integrated::scalar_dimension(&regulated).is_none());
    let mink = Minkowski {}.new_rep(GS.dim);
    let slot = mink.to_symbolic([Atom::num(1)]);
    let vector = symbolica::function!(GS.emr_mom, 1, &slot);
    let tensor = &vector / (&d - Atom::num(1));
    let regulated = Integrated::regulated_for_expansion(&tensor, "test")?;
    ensure!(
        regulated.contains_symbol(GS.dim),
        "the Lorentz dimension is kept: {regulated}"
    );
    ensure!(
        crate::uv::approx::integrated::scalar_dimension(&regulated).is_none(),
        "{regulated}"
    );
    let free = &vector * &epsilon;
    ensure!(matches!(
        Integrated::regulated_for_expansion(&free, "test")?,
        std::borrow::Cow::Borrowed(_)
    ));
    Ok(())
}

// ---------------------------------------------------------------------------
// Round 8b repair (review of 81a8ceab3b: C1-5, RC-B1, RC-S1, RC-S3).
// ---------------------------------------------------------------------------

const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
const BIS: &str = "spenso::{spenso::representation,spenso::self_dual}::bis";
const GAMMA: &str = "spenso::{linear,spenso::tensor}::gamma";

/// The review's closed Dirac trace −8/(tr(γ_a γ^a) − 16), 1/ε at d = 4 − 2ε.
fn dirac_trace_pole() -> String {
    let bis = |n: usize| format!("{BIS}(4,{HEDGE}({n},1))");
    let mink = format!("{MINK}(4,{HEDGE}(85,1))");
    format!(
        "(-8)/({GAMMA}({},{},{mink})*{GAMMA}({},{},{mink})-16)",
        bis(83),
        bis(84),
        bis(84),
        bis(83)
    )
}

/// The closed metric trace T = g_ab g^ab (T − 4 = −2ε at d = 4 − 2ε).
fn metric_trace() -> String {
    let slot = |n: usize| format!("{MINK}(4,{HEDGE}({n},1))");
    format!("({METRIC}({},{}))^2", slot(81), slot(82))
}

/// One run of the typed closed root with the declared data `data` and a
/// SEPARATE callback `X ↦ callback_factor · σ(X)`; also returns how many
/// times the callback was called on a closed scalar (`closed_dot_scalar`).
fn typed_pair(
    case: &Case,
    route: &'static str,
    strict: bool,
    data: &ClosedRootNormalization,
    callback_factor: &Atom,
) -> Result<(Atom, DepthCounts, usize)> {
    let factor = callback_factor.clone();
    let routing = data.external_routing.clone();
    let closed_calls = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
    let seen = closed_calls.clone();
    let callback = move |tensor: &Atom| -> Result<Atom> {
        if super::super::closed_dot_scalar(tensor.as_view()) {
            seen.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
        Ok(&factor * tensor.replace_multiple(&routing))
    };
    let test = DepthTest {
        strict: Some(strict),
        ..FIXED
    };
    let ((value, counts), _) = crate::uv::approx::streamed::with_algebraic_check(|| {
        depth::with_depth_test(test, || {
            depth::collect(|| {
                with_ir_hard(route, || {
                    case.plan.projected_counterterm_closing_root_scalar(
                        &case.graph,
                        case.model(),
                        &case.vakint,
                        &case.settings,
                        &callback,
                        data,
                    )
                })
            })
        })
    });
    let calls = closed_calls.load(std::sync::atomic::Ordering::Relaxed);
    eprintln!(
        "IRH_R8B_REPAIR_TYPED route={route} strict={strict} closed_callback_calls={calls} counts={counts:?}"
    );
    Ok((value?, counts, calls))
}

/// `typed_root` with the declared factor multiplied by `extra`.
fn typed_root_times(graph: &Graph, extra: &Atom) -> Result<ClosedRootNormalization> {
    let mut root = typed_root(graph)?;
    root.factor = &root.factor * extra;
    Ok(root)
}

/// Repair C1-5: `k` of a factor is the ε order of THE FACTOR (its Laurent
/// valuation at d = 4 − 2ε), not the sum of its multiplicative pieces' pole
/// orders: a regular factor spelled as a product (2ε/(4−d) = 1,
/// (d²−16)/(ε(d+4)) = −2) on a quotient vertex or edge has k = 0, so the
/// child route TAKES it (no decline) and gives the literal's value on every
/// route. On 81a8ceab3b each was `declared_order` (k = 1).
#[test]
fn repair_c1_5_a_regular_product_on_the_quotient_side_has_k_zero() -> Result<()> {
    let base_case = Case::new(L3_006, None)?;
    let (base_value, _) = unweighted(&base_case, "off", false)?;
    drop(base_case);
    let mut failures = Vec::new();
    for (name, factor, literal) in [
        ("half", format!("2*{EPS}/(4-{DIM})"), Atom::one()),
        (
            "dsq",
            format!("({DIM}^2-16)/({EPS}*({DIM}+4))"),
            Atom::num(-2),
        ),
    ] {
        for term in ["v1", "v1 ->"] {
            let case = Case::new(&spelled(L3_006, term, &factor), None)?;
            for (route, strict) in [("off", false), ("child", false), ("child", true)] {
                let cell = format!("{name} {term} {route} strict={strict}");
                match unweighted(&case, route, strict) {
                    Ok((value, counts)) => {
                        eprintln!(
                            "IRH_R8B_REPAIR_CELL {cell} poles={} counts={counts:?}",
                            label_poles(&value)
                        );
                        if route == "child"
                            && !(counts.taken >= 1 && counts.declined == 0 && counts.refused == 0)
                        {
                            failures.push(format!("{cell}: not taken: {counts:?}"));
                        }
                        if !same(&value, &(&literal * &base_value))? {
                            failures.push(format!(
                                "{cell}: {} against the literal's {}",
                                label_poles(&value),
                                label_poles(&(&literal * &base_value))
                            ));
                        }
                    }
                    Err(error) => failures.push(format!("{cell}: error {error:#}")),
                }
            }
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// Repair RC-B1: a typed root normalisation that holds a tensor structure
/// (here a closed Dirac trace equal to 1/ε, whose value depends on the
/// dimension the trace is taken in) is not a typed regular form: the child
/// route never certifies or multiplies it; it declines (`opaque_projector`)
/// and equals the generic route. On 81a8ceab3b `excess_poles` read k = 0 and
/// the child took it (taken=1).
#[test]
fn repair_rc_b1_a_traced_typed_factor_is_never_taken() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let spelled = symbolica::try_parse!(&dirac_trace_pole()).map_err(|error| eyre!("{error}"))?;
    let data = typed_root_times(&case.graph, &spelled)?;
    let (off, _, _) = typed_pair(&case, "off", false, &data, &data.factor)?;
    for strict in [false, true] {
        let (child, counts, _) = typed_pair(&case, "child", strict, &data, &data.factor)?;
        ensure!(
            counts.taken == 0 && counts.declined >= 1 && counts.reason == Some("opaque_projector"),
            "strict={strict}: a traced typed factor must decline: {counts:?}"
        );
        ensure!(
            same_scalar(&child, &off)?,
            "strict={strict}: the child differs from off"
        );
    }
    Ok(())
}

/// Repair RC-B1 (the independent gate's cell, helper VH 2): the graph
/// numerator prefactor spelled with a closed metric trace,
/// 1/(1+(T−4)/ε) = −1, and the typed data built from the RAW prefactor (as
/// beta builds it), while the callback multiplies the contracted prefactor.
/// The checker and the multiplier must be one object: the child route
/// declines the raw traced data and every route gives the literal −1's
/// value. On 81a8ceab3b the child applied the raw data.
#[test]
fn repair_rc_b1_a_traced_prefactor_in_the_typed_data_declines() -> Result<()> {
    let literal_dot = times(L3_006, "graph", "numerator_prefactor", "-1");
    let literal_case = Case::new(&literal_dot, None)?;
    let literal_data = typed_root(&literal_case.graph)?;
    let (literal, _, _) = typed_pair(
        &literal_case,
        "off",
        false,
        &literal_data,
        &literal_data.factor,
    )?;
    drop(literal_case);
    let spelling = format!("1/(1+(({})-4)/{EPS})", metric_trace());
    let dot = times(L3_006, "graph", "numerator_prefactor", &spelling);
    let case = Case::new(&dot, None)?;
    let data = typed_root(&case.graph)?;
    let contracted = super::super::regularity::normalise(
        super::super::regularity::Factor::GraphNumeratorPrefactor,
        &case.graph.global_prefactor.num,
    )
    .map_err(|refused| eyre!("{refused}"))?;
    let callback_factor = contracted.replace_multiple(&data.external_routing);
    let mut failures = Vec::new();
    for (route, strict) in [("off", false), ("child", false), ("child", true)] {
        match typed_pair(&case, route, strict, &data, &callback_factor) {
            Ok((value, counts, _)) => {
                if !same_scalar(&value, &literal)? {
                    failures.push(format!(
                        "{route} strict={strict}: differs from the literal −1"
                    ));
                }
                if route == "child"
                    && (counts.taken > 0 || counts.reason != Some("opaque_projector"))
                {
                    failures.push(format!(
                        "{route} strict={strict}: raw traced data taken: {counts:?}"
                    ));
                }
            }
            Err(error) => failures.push(format!("{route} strict={strict}: error {error:#}")),
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// Repair RC-S1: the typed data and the callback must not be two sources of
/// truth. On a closed scalar the evaluation's ONE projection is the typed
/// form `c · σ(X)` (the callback is never called there), on every route, so
/// the child route's certified and multiplied scalar is the generic route's:
/// with data c and an inconsistent callback 2c, child equals off. On
/// 81a8ceab3b: child_equals_off=false, taken=1, declined=0.
#[test]
fn repair_rc_s1_one_projection_object() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let data = typed_root(&case.graph)?;
    let twice = &data.factor * Atom::num(2);
    let (off, _, off_calls) = typed_pair(&case, "off", false, &data, &twice)?;
    let (consistent, _, _) = typed_pair(&case, "off", false, &data, &data.factor)?;
    for strict in [false, true] {
        let (child, counts, calls) = typed_pair(&case, "child", strict, &data, &twice)?;
        ensure!(
            counts.taken >= 1 && counts.declined == 0,
            "strict={strict}: {counts:?}"
        );
        ensure!(
            same_scalar(&child, &off)?,
            "strict={strict}: the child differs from off"
        );
        ensure!(
            calls == 0,
            "strict={strict}: the callback was called on {calls} closed scalar(s)"
        );
    }
    ensure!(
        off_calls == 0,
        "off: the callback was called on {off_calls} closed scalar(s)"
    );
    ensure!(
        same_scalar(&off, &consistent)?,
        "off: the closed scalars took the callback"
    );
    Ok(())
}

// ---- IRH-ATTACK (helper VH 2, 2026-11-22, public patch irh-attack-gltest.patch of
// 20261122T1304Z-gltest-typed-root, base 19ddf454fc): the four adversarial cells as regression tests in their
// PRESENT form on the G-G stack (round 8 integrated + the round-8b repair). T = g_ab g^ab with fresh 4-dimensional
// slots (as every DOT writes its Lorentz slots). Every cell prints an `IRH_ATTACK_CELL` line.

/// T = g_ab g^ab, a closed metric trace with fresh slots.
fn attack_trace() -> String {
    let g = "spenso::{symmetric,linear,real,spenso::tensor}::g";
    let m = |h: u32| {
        format!(
            "spenso::{{spenso::representation,spenso::self_dual}}::mink(4,gammalooprs::{{spenso::index}}::hedge({h},1))"
        )
    };
    format!("({g}({},{}))^2", m(81), m(82))
}

/// One run of `projected_counterterm_scalar` (the typed scalar root projection given as data).
fn attack_scalar(case: &Case, route: &'static str, factor: &Atom) -> Result<(Atom, DepthCounts)> {
    let test = DepthTest {
        strict: Some(false),
        ..FIXED
    };
    let (value, counts) = depth::with_depth_test(test, || {
        depth::collect(|| {
            with_ir_hard(route, || {
                case.plan.projected_counterterm_scalar(
                    &case.graph,
                    case.model(),
                    &case.vakint,
                    &case.settings,
                    factor,
                )
            })
        })
    });
    eprintln!("IRH_ATTACK_SCALAR route={route} counts={counts:?}");
    Ok((value?, counts))
}

/// IRH-ATTACK F-B (typed root from the RAW traced prefactor): the child route never takes a traced typed factor; it
/// DECLINES (counted, `opaque_projector`, repair RC-B1) to the generic route, or the cell is refused on both routes
/// (the 0/0 spelling −2ε/(T−4): vakint's "invalid input integral" on both). The child equals the generic route.
/// Present form, a listed limit (BATCH2-ON-STACK.md "Listed limits"; helper VH 2's PRODUCIBLE.md: no producer writes
/// a trace into a weight): this API multiplies the factor it is GIVEN, so a raw factor whose trace has 4-dimensional
/// slots is read at T = 4 on both routes (the +1 literal's value for the −1 spelling 1/(1+(T−4)/ε)). Beta normalises
/// every weight at entry (round 8a, the trace contracted in d) and passes that one object; with the normalised data
/// the cell equals the literal (`repair_rc_b1_a_traced_prefactor_in_the_typed_data_declines`; helper VH 2's beta
/// cells y-pt-*/a-tr-pf-*, 20261122T1343Z-attack-r8b-repair).
#[test]
fn irh_attack_typed_root_with_a_traced_prefactor_declines() -> Result<()> {
    let t = attack_trace();
    let cells = [
        ("literal+1", "1".to_owned()),
        ("trace-cross", format!("1/(1+({t}-4)/{EPS})")),
        ("trace-one", format!("-2*{EPS}/({t}-4)")),
    ];
    let mut literal_plus_one = None::<Atom>;
    let mut failures = Vec::new();
    for (name, factor) in cells.iter() {
        let dot = times(L3_006, "graph", "numerator_prefactor", factor);
        let case = Case::new(&dot, None)?;
        let off = typed(&case, "off", false);
        let child = typed(&case, "child", false);
        eprintln!(
            "IRH_ATTACK_CELL typed prefactor {name} off={:?} child={:?}",
            off.as_ref()
                .map(|(v, _)| v.to_string())
                .map_err(|e| format!("{e:#}")),
            child
                .as_ref()
                .map(|(v, c)| (v.to_string(), format!("{c:?}")))
                .map_err(|e| format!("{e:#}")),
        );
        match (off, child) {
            (Ok((off, _)), Ok((child, counts))) => {
                if *name == "literal+1" {
                    if counts.taken == 0 {
                        failures.push(format!("{name}: the literal is not taken: {counts:?}"));
                    }
                    literal_plus_one = Some(off.clone());
                } else {
                    if !(counts.taken == 0
                        && counts.declined >= 1
                        && counts.reason == Some("opaque_projector"))
                    {
                        failures.push(format!(
                            "{name}: a traced typed factor must decline: {counts:?}"
                        ));
                    }
                    if *name == "trace-cross"
                        && let Some(reference) = &literal_plus_one
                        && !same_scalar(&off, reference)?
                    {
                        failures.push(format!(
                            "{name}: off is not the T = 4 reading (the +1 literal's value)"
                        ));
                    }
                }
                if !same_scalar(&child, &off)? {
                    failures.push(format!("{name}: the child differs from off"));
                }
            }
            (Err(_), Err(_)) if *name == "trace-one" => {}
            (off, child) => failures.push(format!(
                "{name}: off ok={} child ok={}: the routes disagree on refusing",
                off.is_ok(),
                child.is_ok()
            )),
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// IRH-ATTACK F-B at the public typed scalar entry (`projected_counterterm_scalar`): a factor −1 spelled through a
/// closed metric trace. Present form: the child DECLINES it (`opaque_projector`, repair RC-B1: a traced factor is not
/// a closed scalar) and equals the generic route; both multiply the factor as given (T = 4 with 4-dimensional slots:
/// the +1 value for both spellings). Listed limit as above: the caller passes a normalised factor (beta does).
#[test]
fn irh_attack_typed_scalar_api_with_a_traced_factor_declines() -> Result<()> {
    let t = attack_trace();
    let case = Case::new(L3_006, None)?;
    let parse = |text: &str| -> Result<Atom> {
        symbolica::try_parse!(text).map_err(|error| eyre!("{text}: {error}"))
    };
    let (plus_one, _) = attack_scalar(&case, "off", &Atom::num(1))?;
    let mut failures = Vec::new();
    for (name, text) in [
        ("trace-cross", format!("1/(1+({t}-4)/{EPS})")),
        ("trace-form", format!("{EPS}/({EPS}+{t}-4)")),
    ] {
        let factor = parse(&text)?;
        let (off, _) = attack_scalar(&case, "off", &factor)?;
        let (child, counts) = attack_scalar(&case, "child", &factor)?;
        eprintln!("IRH_ATTACK_CELL scalar api {name} off={off} child={child} counts={counts:?}");
        if !(counts.taken == 0 && counts.declined >= 1 && counts.reason == Some("opaque_projector"))
        {
            failures.push(format!(
                "{name}: a traced typed factor must decline: {counts:?}"
            ));
        }
        if !same(&child, &off)? {
            failures.push(format!("{name}: the child differs from off"));
        }
        if !same(&off, &plus_one)? {
            failures.push(format!(
                "{name}: off is not the T = 4 reading (the +1 value)"
            ));
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// IRH-ATTACK F-A at the source, EXPECTED CURRENT BEHAVIOUR (not the spec): (d−4+2ε)/(2d−8+4ε), 1/2 as a rational
/// function of d and ε and 0/0 at d = 4 − 2ε, is NOT refused by the C1 regulation or the certificate's normalisation:
/// both return 0. A LISTED constructed-input limit (BATCH2-ON-STACK.md "Listed limits"; helper VH 2's PRODUCIBLE.md:
/// no producer writes d or ε into a numerator or weight); not fixed in this batch. This test fails when the
/// behaviour changes (then the listing and this test are updated together).
#[test]
fn irh_attack_a_vanishing_denominator_is_zero_today_listed_limit() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let d = Atom::var(GS.dim);
    let epsilon = Atom::var(GS.dim_epsilon);
    let zz = (&d - Atom::num(4) + Atom::num(2) * &epsilon)
        / (Atom::num(2) * &d - Atom::num(8) + Atom::num(4) * &epsilon);
    let regulated = Integrated::regulated_for_expansion(&zz, "irh_attack").map_err(|error| {
        eyre!("regulated_for_expansion refused (the listed limit changed): {error:#}")
    })?;
    eprintln!(
        "IRH_ATTACK_CELL regulated_for_expansion accepted as {}",
        &*regulated
    );
    ensure!(
        regulated.expand().is_zero(),
        "regulated_for_expansion: {} (expected 0 today)",
        &*regulated
    );
    let normalised = super::super::regularity::normalise(
        super::super::regularity::Factor::QuotientVertexNumerator,
        &zz,
    )
    .map_err(|refused| eyre!("normalise refused (the listed limit changed): {refused}"))?;
    eprintln!("IRH_ATTACK_CELL normalise accepted as {normalised}");
    ensure!(
        normalised.expand().is_zero(),
        "normalise: {normalised} (expected 0 today)"
    );
    Ok(())
}

/// IRH-ATTACK F-A end to end on the unweighted API, EXPECTED CURRENT BEHAVIOUR: quotient vertex v1 of L3_006 times
/// the 0/0 spelling of 1/2 gives 0 on the generic, closed and child routes (the literal 1/2 gives a nonzero
/// counterterm). LISTED constructed-input limit (BATCH2-ON-STACK.md "Listed limits"; helper VH 2's PRODUCIBLE.md);
/// not fixed in this batch.
#[test]
fn irh_attack_a_vanishing_denominator_on_a_quotient_vertex_is_zero_today_listed_limit() -> Result<()>
{
    let zz = format!("({DIM}-4+2*{EPS})/(2*{DIM}-8+4*{EPS})");
    let literal = Case::new(&spelled(L3_006, "v1", "1/2"), None)?;
    let (reference, _) = unweighted(&literal, "off", false)?;
    ensure!(
        !reference.is_zero(),
        "the literal 1/2's counterterm is zero"
    );
    let case = Case::new(&spelled(L3_006, "v1", &zz), None)?;
    for route in ["off", "closed", "child"] {
        let (value, counts) = unweighted(&case, route, false)?;
        eprintln!("IRH_ATTACK_CELL zz v1 {route} value={value} counts={counts:?}");
        ensure!(
            value.is_zero(),
            "{route}: {value} (expected 0 today: the listed limit changed)"
        );
    }
    Ok(())
}
