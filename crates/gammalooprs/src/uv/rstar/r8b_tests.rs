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
use crate::{
    graph::Graph,
    utils::GS,
    uv::{UltravioletGraph as _, approx::integrated::Integrated},
};
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
    let (external, external_routing) = RStarIntegrand::physical_external_routing(
        graph,
        &graph.external_filter::<SuBitGraph>(),
    )?;
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
        .map_or_else(|error| format!("unreadable ({error})"), |poles| poles.to_string())
}

/// The spelling matrix of one position (`term` of `base`): for every
/// spelling of −1 and every route (child STRICT off and on, closed, off)
/// the unweighted counterterm equals the literal's on the projected path
/// (`off`); the literal's equals −1 times `base`'s on that path.
fn unweighted_matrix(label: &str, base: &str, term: &str) -> Result<()> {
    let base_case = Case::new(base, None)?;
    let (base_value, _) = unweighted(&base_case, "off", false)?;
    drop(base_case);
    let mut reference = None::<Atom>;
    let mut failures = Vec::new();
    for (name, factor) in spellings() {
        let case = Case::new(&spelled(base, term, &factor), None)?;
        for (route, strict) in [("off", false), ("closed", false), ("child", false), ("child", true)] {
            let cell = format!("{label} {name} {route} strict={strict}");
            let value = match unweighted(&case, route, strict) {
                Ok((value, _)) => value,
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
    unweighted_matrix("hard v0", L3_006, "v0")
}

/// (a) the review's own carrier: k3j1 (v0·ε³, v1·ε⁻¹).
#[test]
fn round8b_a_hard_vertex_spellings_on_k3j1() -> Result<()> {
    unweighted_matrix("hard v0 k3j1", &k(3, 1), "v0")
}

/// (b) The refutation (topo/refute): the same spellings on the quotient
/// vertex v1 gave 117/256 on the child AND on the generic route.
#[test]
fn round8b_b_quotient_vertex_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("quotient v1", L3_006, "v1")
}

/// (b) the refutation's carrier: k3j1.
#[test]
fn round8b_b_quotient_vertex_spellings_on_k3j1() -> Result<()> {
    unweighted_matrix("quotient v1 k3j1", &k(3, 1), "v1")
}

/// (c) A quotient edge (v1 → v3, edge 4) and a hard edge (v2 → v3).
#[test]
fn round8b_c_quotient_edge_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("quotient edge v1->v3", L3_006, "v1 ->")
}

#[test]
fn round8b_c_hard_edge_spellings_equal_the_literal() -> Result<()> {
    unweighted_matrix("hard edge v2->v3", L3_006, "v2 ->")
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
            for (route, strict) in [("off", false), ("closed", false), ("child", false), ("child", true)] {
                let cell = format!("typed {term} {name} {route} strict={strict}");
                let value = match typed(&case, route, strict) {
                    Ok((value, counts)) => {
                        eprintln!("IRH_R8B_CELL {cell} counts={counts:?}");
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
    let weighted =
        move |tensor: &Atom| -> Result<Atom> { Ok(tensor * &prefactor) };
    let (weighted_unit, weighted_counts) = opaque(&case, "child", false, &weighted)?;
    let (weighted_today, _) = opaque(&case, "closed", false, &weighted)?;
    ensure!(same(&weighted_unit, &weighted_today)?, "k3j1 weighted with projector 1");
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
        ensure!(typed_counts.exceeded == 0, "typed strict={strict}: {typed_counts:?}");
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
        ensure!(probes == 0, "{label}: the callback was probed {probes} time(s)");
    }
    Ok(())
}

/// (f) Finding 5: a callback that forwards a typed alarm (here only on the
/// probe's input, as a callback delegating to a nested R* evaluation could)
/// was skipped by the probe and the correction admitted. Every alarm a
/// callback raises must be fatal; without a probe none is raised here, and
/// the child route declines the opaque callback.
#[test]
fn round8b_f_a_callback_alarm_is_never_swallowed() -> Result<()> {
    let case = Case::new(&k(3, 1), None)?;
    let forwarding = |tensor: &Atom| -> Result<Atom> {
        if tensor.contains_symbol(probe_symbol()) {
            return Err(depth::alarm(
                "IRH-DEPTH: a nested alarm forwarded by the root projector".to_owned(),
            ));
        }
        Ok(tensor.clone())
    };
    let (today, _) = opaque(&case, "closed", false, &forwarding)?;
    match opaque(&case, "child", false, &forwarding) {
        Ok((child, counts)) => {
            ensure!(same(&child, &today)?, "the value");
            ensure!(
                counts.callback_unprobed == 0 && counts.taken == 0,
                "the alarm was swallowed and the correction admitted: {counts:?}"
            );
        }
        Err(error) => ensure!(
            depth::is_alarm(&error),
            "an error that lost the alarm class: {error:#}"
        ),
    }
    // A callback that always raises the alarm: fatal on every route.
    let always = |_: &Atom| -> Result<Atom> {
        Err(depth::alarm("IRH-DEPTH: an alarm raised by the root projector".to_owned()))
    };
    for route in ["closed", "child", "off"] {
        match opaque(&case, route, false, &always) {
            Ok(_) => eyre::bail!("{route}: the alarm was swallowed"),
            Err(error) => ensure!(depth::is_alarm(&error), "{route}: {error:#}"),
        }
    }
    Ok(())
}
