//! IRH-DEPTH on B_self_energy_L3_006 (results/prepared/L3, as in the
//! review of U7b, depth/REVIEW.md): its correction 1 has a two-loop hard part
//! (vertex v0) and a one-loop quotient (vertex v1, an external vertex), and
//! takes the child route (omega 1). The variations multiply the numerator of
//! v0 or v1 by a power of epsilon (and the aggregate graph numerator by
//! their product), exactly as the review's fixtures do. Every test compares
//! the complete counterterm of the graph on the child route with today's
//! path (`GAMMALOOP_IR_HARD=closed`) on the same proper UV children.

use eyre::{Result, ensure};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{feyngen::feynkit::FeynmanDiagramGammaLoopExt, initialisation::test_initialise};

use super::{
    RStarPlan,
    closed::with_ir_hard,
    depth::{self, HookSite, RowHook},
    recursive_tests::{colour_factored, load, settings},
    tests::canonical_tensor,
};

const L3_006: &str = include_str!("testdata/l3-d6.dot");

fn epsilon() -> Atom {
    Atom::var(crate::utils::GS.dim_epsilon)
}

/// The review's fixture: the first `numerator="…"` of the `v<k>` line and of
/// the `graph` line (not one of `…_numerator=`) times `ε^power`.
fn scaled(dot: &str, powers: &[(&str, i64)]) -> String {
    let scale = |line: &str, power: i64| -> String {
        if power == 0 {
            return line.to_owned();
        }
        let key = "numerator=\"";
        let mut from = 0;
        while let Some(found) = line[from..].find(key) {
            let start = from + found;
            let before = line[..start].chars().last();
            if before.is_some_and(|c| c.is_ascii_lowercase() || c == '_') {
                from = start + key.len();
                continue;
            }
            let body = start + key.len();
            let end = body + line[body..].find('"').expect("a closed numerator string");
            return format!(
                "{}({})*(gammalooprs::{{}}::ε)^({power}){}",
                &line[..body],
                &line[body..end],
                &line[end..]
            );
        }
        panic!("no numerator on {line}");
    };
    let total = powers.iter().map(|(_, power)| power).sum();
    dot.split_inclusive('\n')
        .map(|line| {
            if line.starts_with("  graph [") {
                return scale(line, total);
            }
            for (vertex, power) in powers {
                if line.starts_with(&format!("  {vertex} [")) {
                    return scale(line, *power);
                }
            }
            line.to_owned()
        })
        .collect()
}

/// FNV-1a (64 bits) of a fixture.
fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

/// The variations are byte for byte the review's inputs
/// (`depth/<run>/input/L3/diagram-000006.dot`).
#[test]
fn the_variations_are_the_reviewed_inputs() {
    for (powers, hash, length) in [
        (&[("v0", 3), ("v1", -1)][..], 0xa9ef_30bf_146a_1514, 47361),
        (&[("v1", -1)][..], 0xebad_e5c4_aa6f_0719, 47334),
        (&[("v0", 3)][..], 0x7805_e917_fd4a_9a6f, 47332),
        (&[("v0", -1)][..], 0xd343_de22_357a_640d, 47334),
        (&[("v0", 2), ("v1", -2)][..], 0xfc1d_eae0_8537_647c, 47333),
        (&[("v0", 3), ("v1", -3)][..], 0x9b0a_64e8_c3d4_3d18, 47333),
    ] {
        let fixture = scaled(L3_006, powers);
        assert_eq!((fnv(&fixture), fixture.len()), (hash, length), "{powers:?}");
    }
}

/// Whether two counterterms are the same tensor.
fn same(actual: &Atom, expected: &Atom) -> Result<bool> {
    let difference = canonical_tensor(&(canonical_tensor(actual)? - canonical_tensor(expected)?))?;
    Ok(difference.together().is_zero())
}

struct Routes {
    today: Atom,
    child: Atom,
    /// IRH-DEPTH counts of the child-route evaluation.
    counts: Counts,
}

#[derive(Clone, Copy, Debug)]
struct Counts {
    corrections: usize,
    retried: usize,
    /// Retries to the bound beyond the measured need.
    raised: usize,
    refused: usize,
}

/// The counterterm of `dot` on today's path and on the child route (with
/// `hook` rewriting its rows and its first depth at least `start`), on the
/// same proper UV children (computed on today's path).
fn routes(dot: &str, hook: Option<RowHook>, start: usize) -> Result<Routes> {
    test_initialise()?;
    let vakint = Vakint::new()?;
    let diagram = colour_factored(&load(dot)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let settings = settings();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
    let children = with_ir_hard("closed", || {
        plan.evaluate_uv_children(&graph, model, &mass, &vakint, &settings)
    })?;
    let today = with_ir_hard("closed", || {
        plan.counterterm_with_children(&graph, model, &vakint, &settings, children.clone())
    })?;
    let before = depth::totals();
    let run = || {
        depth::with_start_extra(start, || {
            with_ir_hard("child", || {
                plan.counterterm_with_children(&graph, model, &vakint, &settings, children.clone())
            })
        })
    };
    let child = match hook {
        Some(hook) => depth::with_row_hook(hook, run),
        None => run(),
    }?;
    let after = depth::totals();
    let [corrections, retried, raised, refused] = [0, 1, 2, 3].map(|k| after[k] - before[k]);
    let counts = Counts {
        corrections,
        retried,
        raised,
        refused,
    };
    eprintln!("IRH_DEPTH_TEST {counts:?}");
    Ok(Routes {
        today,
        child,
        counts,
    })
}

/// The review's counterexample: v0·ε³ makes every group of the hard part
/// vanish at the first depth (its valuation is 1), v1·ε⁻¹ gives the quotient
/// a double pole (bound 2). The pole of their product was lost on the child
/// route (+117/32 against -117/64 C_A³/ε in the review's units).
#[test]
fn a_hard_part_vanishing_at_the_first_depth_keeps_its_quotient_pole() -> Result<()> {
    let routes = routes(&scaled(L3_006, &[("v0", 3), ("v1", -1)]), None, 0)?;
    ensure!(!routes.today.is_zero(), "the variation must have poles");
    ensure!(
        same(&routes.child, &routes.today)?,
        "child route {} differs from today's {}",
        routes.child,
        routes.today
    );
    // The vanished groups show no insertion: the bound stands for their
    // valuation and raises the depth to 2, where the hard part is shown and
    // the child route itself gets the pole; the guard never fires.
    let Counts {
        corrections,
        retried,
        raised,
        refused,
    } = routes.counts;
    ensure!(
        corrections >= 1 && retried >= 1 && raised >= 1 && refused == 0,
        "{:?}",
        routes.counts
    );
    Ok(())
}

/// Positive control (review: quotient-v1): v1·ε⁻¹ alone. The shown rows
/// measure the double quotient pole and the planner retries at depth 2.
#[test]
fn a_measured_quotient_pole_still_retries() -> Result<()> {
    let routes = routes(&scaled(L3_006, &[("v1", -1)]), None, 0)?;
    ensure!(
        same(&routes.child, &routes.today)?,
        "child route {} differs from today's {}",
        routes.child,
        routes.today
    );
    let Counts {
        corrections,
        retried,
        raised,
        refused,
    } = routes.counts;
    ensure!(
        corrections >= 1 && retried >= 1 && raised == 0 && refused == 0,
        "{:?}",
        routes.counts
    );
    Ok(())
}

/// Hard-only control (review: single-softer): v0·ε³ alone. The hard part
/// vanishes at the first depth, but the quotient bound is its single loop:
/// nothing can hide, no retry, no refusal.
#[test]
fn a_vanishing_hard_part_under_a_simple_pole_needs_nothing() -> Result<()> {
    let routes = routes(&scaled(L3_006, &[("v0", 3)]), None, 0)?;
    ensure!(
        same(&routes.child, &routes.today)?,
        "child route {} differs from today's {}",
        routes.child,
        routes.today
    );
    let Counts {
        corrections,
        retried,
        raised,
        refused,
    } = routes.counts;
    ensure!(
        corrections >= 1 && retried == 0 && raised == 0 && refused == 0,
        "{:?}",
        routes.counts
    );
    Ok(())
}

/// The trailing ε order of `value` (`None` when zero).
fn order(value: &Atom) -> Option<i64> {
    let series = value
        .series(
            crate::utils::GS.dim_epsilon,
            0,
            symbolica::poly::series::SeriesDepth::relative(1),
        )
        .ok()?;
    if series.is_zero() {
        return None;
    }
    series.get_trailing_exponent().numerator().to_i64()
}

/// A row of a nonzero group whose coefficient vanishes at the first depth
/// and whose quotient factor alone has the deeper pole (the review's
/// row-dependent valuation): on the unmodified graph, the hook turns the row
/// `(h, E)` of the top group (remainder loops = hard loops) of the top-level
/// correction with insertion `target` into
/// `(ε³ h, E/ε²)`, its coefficient kept only within the attempt's depth
/// (`h` is exact through `ε^(extra-1)`, so `ε³ h` is known there exactly
/// and is zero below its own order; `shown = false` keeps it zero at every
/// depth). With `h` from `ε⁻²`, `ε³ h` starts at `ε`: zero at depth 1. Its
/// quotient factor `Z̃(E/ε²)` starts at `ε⁻³` (bound 3); every other row's
/// at `ε⁻¹` (bound 1).
fn vanished_row(target: Atom, shown: bool) -> RowHook {
    Box::new(
        move |site: &HookSite, rows: &mut Vec<(Atom, Atom, usize)>| {
            if site.hard_loops != 2 || site.quotient_loops != 1 {
                return;
            }
            for (h, insertion, loops) in rows.iter_mut() {
                if *loops != site.hard_loops || *insertion != target {
                    continue;
                }
                let within = order(h).is_some_and(|order| order + 3 <= site.extra as i64 - 1);
                *h = if shown && within {
                    &*h * epsilon().pow(3)
                } else {
                    Atom::Zero
                };
                *insertion = &*insertion / epsilon().pow(2);
            }
        },
    )
}

/// The rows of the top-level correction at the first depth.
fn first_rows() -> Result<Vec<(Atom, Atom, usize)>> {
    let rows = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
    let record = rows.clone();
    let hook: RowHook = Box::new(
        move |site: &HookSite, shown: &mut Vec<(Atom, Atom, usize)>| {
            if site.hard_loops == 2 && site.quotient_loops == 1 && site.extra == 1 {
                *record.borrow_mut() = shown.clone();
            }
        },
    );
    routes(L3_006, Some(hook), 0)?;
    Ok(std::rc::Rc::try_unwrap(rows)
        .map_err(|_| eyre::eyre!("the hook outlived its run"))?
        .into_inner())
}

#[test]
fn a_vanished_row_of_a_nonzero_group_keeps_its_quotient_pole() -> Result<()> {
    let rows = first_rows()?;
    let top = rows.iter().map(|(_, _, loops)| *loops).max().unwrap_or(0);
    let group = rows
        .iter()
        .filter(|(h, _, loops)| *loops == top && !h.is_zero())
        .collect::<Vec<_>>();
    ensure!(group.len() >= 2, "the top group needs two rows: {rows:?}");
    for (h, target, _) in group {
        if order(h) != Some(-2) {
            continue;
        }
        // The value of the hooked correction: from depth 3 on, every row's
        // coefficient is shown through the order its quotient needs.
        let reference = routes(L3_006, Some(vanished_row(target.clone(), true)), 3)?;
        // The same without the row's contribution: the row must reach the
        // pole part.
        let without = routes(L3_006, Some(vanished_row(target.clone(), false)), 3)?;
        if same(&reference.child, &without.child)? {
            continue;
        }
        let routes = routes(L3_006, Some(vanished_row(target.clone(), true)), 0)?;
        ensure!(
            same(&routes.child, &reference.child)?,
            "child route {} differs from the reference {}",
            routes.child,
            reference.child
        );
        // The vanished row's quotient factor is measured (no raise to a
        // bound): the planner retries at depth 3 and the child route gets the
        // reference; the guard never fires.
        let Counts {
            corrections,
            retried,
            raised,
            refused,
        } = routes.counts;
        ensure!(
            corrections >= 1 && retried >= 1 && raised == 0 && refused == 0,
            "{:?}",
            routes.counts
        );
        return Ok(());
    }
    Err(eyre::eyre!(
        "no row of the top group reaches the pole part: {rows:?}"
    ))
}
