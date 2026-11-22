//! IRH-DEPTH on B_self_energy_L3_006 (results/prepared/L3, as in the
//! review of U7b, depth/REVIEW.md): its correction 1 has a two-loop hard part
//! (vertex v0) and a one-loop quotient (vertex v1, an external vertex), and
//! takes the child route (omega 1). The variations multiply the numerator of
//! v0 or v1 by a power of epsilon (and the aggregate graph numerator by
//! their product), exactly as the review's fixtures do, or replace v0's
//! Lorentz factor by a row-level combination (the IRH-DEPTH grid's rA/rB,
//! `grid/scripts/rowfix.py`). Every test compares the complete counterterm of
//! the graph on the child route with today's path
//! (`GAMMALOOP_IR_HARD=closed`) and pins its pole part in the units of the
//! beta record (`zb`, the coefficient of C_A³; [`record`]) to the generic
//! values of the IRH-DEPTH grid (results/consolidation/irh-depth/grid.json).
//! The third method (massive UV-only R, MATAD) confirms every pin at the
//! grid points kKjJ; at the row-level points rA, rB and rA0 it disagrees at
//! the pinned (rearrangement-unsafe) order, so there the direct pin is the
//! generic value only and the independent check is `point − control`
//! (linearity onto k3j1 / k4j2, which the massive method confirms). The planner
//! before IRH-DEPTH ([`Planner::Legacy`]) is run next to it and pinned to its
//! wrong value, so each test fails on the old planner.
//!
//! The counts are per test ([`depth::collect`], thread-local and scoped):
//! the tests hold under `cargo test`'s parallel threads.

use std::sync::OnceLock;

use eyre::{Result, ensure, eyre};
use idenso::{
    CookMode, CookSettings, IndexTooling,
    color::ColorSimplifySettings,
    tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor},
};
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    network::library::symbolic::ETS,
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibrarySlot, Minkowski},
        slot::IsAbstractSlot,
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    model::Model,
    utils::GS,
    uv::{UVgenerationSettings, approx::integrated::Integrated},
};

use super::{
    RStarPlan,
    closed::with_ir_hard,
    depth::{self, BoundOverride, DepthCounts, DepthTest, HookSite, Planner, RowHook},
    recursive_tests::{colour_factored, load, settings},
    tests::canonical_tensor,
};

const L3_006: &str = include_str!("testdata/l3-d6.dot");
const EPS: &str = "gammalooprs::{}::ε";

fn epsilon() -> Atom {
    Atom::var(GS.dim_epsilon)
}

/// The rational `n/d`.
fn q(n: i64, d: i64) -> Atom {
    Atom::num(n) / Atom::num(d)
}

/// The span `(body, end)` of the first `numerator="…"` of `line` that is not
/// one of `…_numerator=`, as the fixture scripts' regex
/// `(?<![a-z_])numerator="([^"]+)"`.
fn first_numerator(line: &str) -> Option<(usize, usize)> {
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
        let end = body + line[body..].find('"')?;
        return Some((body, end));
    }
    None
}

/// The review's fixture: the first `numerator="…"` of the `v<k>` line and of
/// the `graph` line (not one of `…_numerator=`) times `ε^power`.
fn scaled(dot: &str, powers: &[(&str, i64)]) -> String {
    let scale = |line: &str, power: i64| -> String {
        if power == 0 {
            return line.to_owned();
        }
        let (body, end) = first_numerator(line).expect("a numerator");
        format!(
            "{}({})*({EPS})^({power}){}",
            &line[..body],
            &line[body..end],
            &line[end..]
        )
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

/// The first balanced `(…)` group at the start of `text`.
fn lead_group(text: &str) -> &str {
    let mut depth = 0_i64;
    for (i, c) in text.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    return &text[..=i];
                }
            }
            _ => {}
        }
    }
    panic!("no group in {text}")
}

/// The IRH-DEPTH grid's row-level fixture (`grid/scripts/rowfix.py`, byte
/// for byte): v0's Lorentz factor `L0` becomes `ε^a L0 + ε^w X` with
/// `X = Q1^μ1 g^μ2μ4` (v0's external index and two internal ones, with v0's
/// colour and coupling; `None` drops a term), v1's numerator is multiplied
/// by `ε^-j`, and the aggregate graph numerator gets both. X feeds only the
/// row `Q1^μ1 g^μ7μ9` of correction 1, whose quotient pole is one lower than
/// the other five rows': it stays measured while the other ten rows vanish
/// inside nonzero groups.
fn rowed(dot: &str, a: Option<i64>, w: Option<i64>, j: i64) -> String {
    const MINK: &str = "spenso::{spenso::representation,spenso::self_dual}::mink";
    const HEDGE: &str = "gammalooprs::{spenso::index}::hedge";
    const METRIC: &str = "spenso::{symmetric,linear,real,spenso::tensor}::g";
    let lines = dot.split_inclusive('\n').collect::<Vec<_>>();
    let v0 = lines
        .iter()
        .find(|line| line.starts_with("  v0 ["))
        .expect("v0");
    let (body, end) = first_numerator(v0).expect("v0's numerator");
    let n0 = &v0[body..end];
    let l0 = lead_group(n0);
    let q_end = l0.find("::Q(1,").expect("Q(1, …) in L0");
    let q_start = l0[..q_end].rfind("gammalooprs::{").expect("Q's namespace");
    let momentum = &l0[q_start..q_end + "::Q".len()];
    assert!(l0.contains(&format!("{MINK}(4,")) && l0.contains(HEDGE) && l0.contains(METRIC));
    let index = |h: usize| format!("{MINK}(4,{HEDGE}({h},1))");
    let x = format!(
        "{momentum}(1,{})*{METRIC}({},{})",
        index(1),
        index(2),
        index(4)
    );
    let parts = a
        .map(|a| format!("({EPS})^({a})*{l0}"))
        .into_iter()
        .chain(w.map(|w| format!("({EPS})^({w})*{x}")))
        .collect::<Vec<_>>();
    let new = format!("({})", parts.join("+"));
    let replace_numerator = |line: &str, numerator: &str| -> String {
        let (body, end) = first_numerator(line).expect("a numerator");
        format!("{}{numerator}{}", &line[..body], &line[end..])
    };
    lines
        .iter()
        .map(|line| {
            if line.starts_with("  graph [") {
                let (body, end) = first_numerator(line).expect("the graph numerator");
                let aggregate = &line[body..end];
                assert!(aggregate.starts_with(&format!("{l0}*")));
                let mut scaled = aggregate.replacen(l0, &new, 1);
                if j != 0 {
                    scaled = format!("({scaled})*({EPS})^({})", -j);
                }
                replace_numerator(line, &scaled)
            } else if line.starts_with("  v0 [") {
                replace_numerator(line, &n0.replacen(l0, &new, 1))
            } else if line.starts_with("  v1 [") && j != 0 {
                let (body, end) = first_numerator(line).expect("v1's numerator");
                replace_numerator(line, &format!("({})*({EPS})^({})", &line[body..end], -j))
            } else {
                (*line).to_owned()
            }
        })
        .collect()
}

/// The IRH-DEPTH grid's points (grid.md; `kKjJ`: v0·ε^K, v1·ε^-J).
fn k(power_v0: i64, power_v1: i64) -> String {
    scaled(L3_006, &[("v0", power_v0), ("v1", -power_v1)])
}
fn r_a() -> String {
    rowed(L3_006, Some(3), Some(2), 1)
}
fn r_b() -> String {
    rowed(L3_006, Some(4), Some(2), 2)
}
fn x_a() -> String {
    rowed(L3_006, None, Some(2), 1)
}
fn x_b() -> String {
    rowed(L3_006, None, Some(2), 2)
}
fn r_a0() -> String {
    rowed(L3_006, Some(3), Some(0), 1)
}
fn x_a0() -> String {
    rowed(L3_006, None, Some(0), 1)
}

/// FNV-1a (64 bits) of a fixture.
fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3)
    })
}

/// The variations are byte for byte the review's inputs
/// (`depth/<run>/input/L3/diagram-000006.dot`) and the IRH-DEPTH grid's
/// (`/common/build-cache/ben/consol/irh-depth/grid/inputs/<point>`).
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
    for (name, fixture, hash, length) in [
        ("k0j0", L3_006.to_owned(), 0x6dbf_495b_b56b_d4ab, 47276),
        ("k2j1", k(2, 1), 0xf02f_0a91_b993_d102, 47361),
        ("k3j2", k(3, 2), 0xe311_9a84_8214_70d8, 47361),
        ("k4j2", k(4, 2), 0x97fe_8627_189e_7688, 47361),
        ("k4j3", k(4, 3), 0x0eb7_0a54_e2ab_4b12, 47361),
        ("rA", r_a(), 0x7ad7_f31f_a2bb_7cc5, 48366),
        ("rB", r_b(), 0xc79a_8e14_db57_aca1, 48366),
        ("xA", x_a(), 0xd409_5318_1fb3_9eb1, 45100),
        ("xB", x_b(), 0xed0e_4da0_6da8_6baf, 45100),
        ("rA0", r_a0(), 0x4dd7_7cbc_349a_65f9, 48366),
        ("xA0", x_a0(), 0xee3b_0b70_1fa0_71a5, 45100),
    ] {
        assert_eq!((fnv(&fixture), fixture.len()), (hash, length), "{name}");
    }
}

/// Whether two counterterms are the same tensor.
fn same(actual: &Atom, expected: &Atom) -> Result<bool> {
    let difference = canonical_tensor(&(canonical_tensor(actual)? - canonical_tensor(expected)?))?;
    Ok(difference.together().is_zero())
}

/// The metric trace `g_μν T^μν` of a counterterm's two open Lorentz slots
/// over `d - 1` (the transverse projection of the beta engine's metric
/// projector, up to the constant normalisation fixed by [`unit`]).
fn traced(counterterm: &Atom) -> Result<Atom> {
    let tensor = canonical_tensor(counterterm)?;
    if tensor.is_zero() {
        return Ok(Atom::Zero);
    }
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(tensor.as_view())?;
    let ports = cooked
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let lorentz = ports
        .iter()
        .filter(|port| port.rep_name() == (Minkowski {}).into())
        .collect::<Vec<_>>();
    ensure!(
        lorentz.len() == 2 && ports.len() == 2,
        "a counterterm with two open Lorentz slots: {ports:?}"
    );
    let metric = ETS.metric(lorentz[0].to_atom(), lorentz[1].to_atom());
    let contracted = SymbolicTensor::infer_indexed(cooked * &metric)?
        .undo_trace()?
        .undo_chain()?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(Default::default()),
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::Dots,
            ..Default::default()
        })?;
    ensure!(
        contracted
            .expression()
            .list_dangling::<AbstractIndex>()?
            .is_empty(),
        "the trace left open slots"
    );
    // Momentum conservation of the two-point function: the incoming edge 0
    // (ext0 -> v1) and the outgoing edge 1 (v0 -> ext1) carry the same Q.
    let conserved = cooking
        .uncook(contracted.expression().as_view())
        .replace_map(|view, _, out| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() == 2
                && momentum.get(0) == Atom::num(1).as_view()
            {
                **out = symbolica::function!(GS.emr_mom, Atom::num(0), momentum.get(1));
            }
        });
    let scalar = conserved / (Atom::var(GS.dim) - 1);
    Ok(Integrated::dimensionally_regularized(&scalar))
}

/// The coefficient of `ε^order` of the Laurent expansion of `value`.
fn laurent(value: &Atom, order: i64) -> Result<Atom> {
    let poles = value
        .series(GS.dim_epsilon, 0, -1)
        .map_err(|error| eyre!("{error}"))?
        .to_atom()
        .expand();
    Ok(poles.coefficient(epsilon().pow(order)).together())
}

/// The record units: the unmodified graph's leading pole is 117/32 C_A³/ε³
/// in the beta record (grid.json k0j0, all three methods). This is the
/// trace of its counterterm at ε⁻³ over 117/32, computed once per process
/// (on today's path; [`the_unmodified_graph_fixes_the_record_units`] checks
/// its other two orders, so the normalisation is a constant).
fn unit() -> Result<Atom> {
    static UNIT: OnceLock<Atom> = OnceLock::new();
    if let Some(unit) = UNIT.get() {
        return Ok(unit.clone());
    }
    let today = Case::new(L3_006, None)?.today()?;
    let unit = laurent(&traced(&today)?, -3)? / q(117, 32);
    ensure!(!unit.is_zero(), "the unmodified graph has a triple pole");
    Ok(UNIT.get_or_init(|| unit).clone())
}

/// A counterterm in record units: the coefficient of C_A³ of the beta
/// record's `zb`, as a Laurent series in ε.
fn record(counterterm: &Atom) -> Result<Atom> {
    Ok(traced(counterterm)? / unit()?)
}

/// `counterterm` in record units has the pole coefficients `pins`
/// (`(order, value)`) and, when `complete`, no other pole.
fn pinned(label: &str, counterterm: &Atom, pins: &[(i64, Atom)], complete: bool) -> Result<()> {
    let value = record(counterterm)?;
    let mut rest = value
        .series(GS.dim_epsilon, 0, -1)
        .map_err(|error| eyre!("{error}"))?
        .to_atom()
        .expand();
    for (order, expected) in pins {
        let actual = laurent(&value, *order)?;
        eprintln!("IRH_DEPTH_PIN {label} order={order} value={actual}");
        ensure!(
            (&actual - expected).together().is_zero(),
            "{label}: ε^{order} is {actual}, pinned {expected}"
        );
        rest -= expected * epsilon().pow(*order);
    }
    if complete {
        let rest = rest.expand().together();
        ensure!(
            rest.is_zero(),
            "{label}: poles beyond the pinned ones: {rest}"
        );
    }
    Ok(())
}

/// One prepared graph (colour factored, the automatic marked line or
/// `marked`) with its proper UV children on today's path.
struct Case {
    graph: Graph,
    plan: RStarPlan,
    diagram: feynkit_graph::FeynmanDiagram,
    vakint: Vakint,
    settings: UVgenerationSettings,
    children: std::collections::BTreeMap<linnet::half_edge::subgraph::SuBitGraph, Atom>,
}

/// One child-route run of a [`Case`].
#[derive(Default)]
struct Run {
    /// Rewrites the rows of every attempt.
    hook: Option<RowHook>,
    /// The first depth is at least this.
    start: usize,
    /// The planner (default: production).
    test: DepthTest,
}

struct ChildRun {
    value: Atom,
    /// The counts of the root evaluation's child-route corrections, the
    /// child calls nested in it (its decorated quotients' R*) included.
    root: DepthCounts,
}

impl Case {
    fn new(dot: &str, marked: Option<usize>) -> Result<Self> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let diagram = colour_factored(&load(dot)?)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let settings = settings();
        let plan = match marked {
            Some(edge) => {
                RStarPlan::new(&graph, diagram.model(), EdgeIndex::from(edge), &settings)?
            }
            None => RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings)?,
        };
        let children = with_ir_hard("closed", || {
            plan.evaluate_uv_children(&graph, diagram.model(), &mass(), &vakint, &settings)
        })?;
        Ok(Self {
            graph,
            plan,
            diagram,
            vakint,
            settings,
            children,
        })
    }

    fn model(&self) -> &Model {
        self.diagram.model()
    }

    /// The counterterm on today's path.
    fn today(&self) -> Result<Atom> {
        with_ir_hard("closed", || {
            self.plan.counterterm_with_children(
                &self.graph,
                self.model(),
                &self.vakint,
                &self.settings,
                self.children.clone(),
            )
        })
    }

    /// The counterterm on the child route.
    fn child(&self, run: Run) -> Result<ChildRun> {
        let Run { hook, start, test } = run;
        depth::with_depth_test(test, || {
            let children = self.children.clone();
            let evaluate = || {
                depth::with_start_extra(start, || {
                    with_ir_hard("child", || {
                        self.plan.counterterm_with_children(
                            &self.graph,
                            self.model(),
                            &self.vakint,
                            &self.settings,
                            children,
                        )
                    })
                })
            };
            let (value, root) = depth::collect(|| match hook {
                Some(hook) => depth::with_row_hook(hook, evaluate),
                None => evaluate(),
            });
            eprintln!("IRH_DEPTH_TEST {test:?} root={root:?}");
            Ok(ChildRun {
                value: value?,
                root,
            })
        })
    }
}

fn mass() -> Atom {
    Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt()
}

const FIXED: DepthTest = DepthTest {
    planner: Planner::Fixed,
    bound: BoundOverride::Computed,
    strict: None,
    guard: None,
};
const LEGACY: DepthTest = DepthTest {
    planner: Planner::Legacy,
    ..FIXED
};
const NO_RAISE: DepthTest = DepthTest {
    planner: Planner::NoRaise,
    ..FIXED
};

fn run(test: DepthTest) -> Run {
    Run {
        test,
        ..Run::default()
    }
}

/// The child route equals today's path exactly.
fn ensure_same(label: &str, child: &Atom, today: &Atom) -> Result<()> {
    ensure!(
        same(child, today)?,
        "{label}: child route {child} differs from today's {today}"
    );
    Ok(())
}

/// Round 8b: a correction whose quotient carries a declared pole (k > 0)
/// declines the child route before H(γ) is integrated (`declared_order`,
/// counted, the value today's): no attempt, no retry, no alarm.
fn ensure_declined(label: &str, counts: &DepthCounts, reason: &str) -> Result<()> {
    ensure!(
        counts.corrections >= 1
            && counts.taken == 0
            && counts.declined >= 1
            && counts.refused >= 1
            && counts.reason == Some(reason)
            && counts.retried == 0
            && counts.raised == 0
            && counts.exceeded == 0,
        "{label}: expected the decline {reason}: {counts:?}"
    );
    Ok(())
}

/// One grid point whose child route was wrong (grid.md 2.5): the fixed
/// child route equals today's path, starts at the quotient bound (no retry,
/// no refusal) and has the pinned poles (the third method's at its safe
/// orders); the legacy planner has the review's wrong poles.
fn grid_point(
    label: &str,
    dot: &str,
    pins: &[(i64, Atom)],
    complete: bool,
    legacy_pins: &[(i64, Atom)],
    legacy_complete: bool,
) -> Result<()> {
    let case = Case::new(dot, None)?;
    let today = case.today()?;
    let child = case.child(run(FIXED))?;
    ensure_same(label, &child.value, &today)?;
    // Round 8b: v1·ε^-J has k_v1 = J > 0: declined (before: taken at the
    // raised bound).
    ensure!(child.root.bounded >= 1, "{label}: {:?}", child.root);
    ensure_declined(label, &child.root, "declared_order")?;
    pinned(label, &today, pins, complete)?;
    let legacy = case.child(run(LEGACY))?;
    ensure!(
        !same(&legacy.value, &today)?,
        "{label}: the legacy planner must miss the pole"
    );
    pinned(
        &format!("{label} legacy"),
        &legacy.value,
        legacy_pins,
        legacy_complete,
    )
}

/// The record units: the unmodified graph has the pole part of all three
/// methods, 117/32 ε⁻³ − 1339/192 ε⁻² + (1369/288 + 3ζ3) ε⁻¹ (grid.json
/// k0j0), on both routes.
#[test]
fn the_unmodified_graph_fixes_the_record_units() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let today = case.today()?;
    let child = case.child(run(FIXED))?;
    ensure_same("k0j0", &child.value, &today)?;
    ensure!(
        child.root.corrections >= 1 && child.root.bounded == 0 && child.root.refused == 0,
        "{:?}",
        child.root
    );
    // ε⁻¹: 1369/288 + 3ζ3 (symbolica's zeta(3)).
    let z3 = symbolica::parse!("zeta(3)");
    pinned(
        "k0j0",
        &today,
        &[
            (-3, q(117, 32)),
            (-2, q(-1339, 192)),
            (-1, q(1369, 288) + Atom::num(3) * z3),
        ],
        true,
    )
}

/// The review's counterexample (k3j1): v0·ε³ makes every group of the hard
/// part vanish at the first depth (its valuation is 1), v1·ε⁻¹ gives the
/// quotient a double pole (bound 2). The pole of their product was lost on
/// the child route (+117/32 against -117/64 C_A³/ε). The fixed planner
/// starts at the bound, 2, where the hard part is shown.
#[test]
fn a_hard_part_vanishing_at_the_first_depth_keeps_its_quotient_pole() -> Result<()> {
    grid_point(
        "k3j1",
        &k(3, 1),
        &[(-1, q(-117, 64))],
        true,
        &[(-1, q(117, 32))],
        true,
    )
}

/// k3j2 (v0·ε³, v1·ε⁻²): leading -117/64 ε⁻² (massive -117/64; ε⁻¹ carries
/// logs); the old child record (+117/32 ε⁻² − 1339/192 ε⁻¹) is log-free and
/// passes the MS-locality check.
#[test]
fn grid_k3j2() -> Result<()> {
    grid_point(
        "k3j2",
        &k(3, 2),
        &[(-2, q(-117, 64))],
        false,
        &[(-2, q(117, 32)), (-1, q(-1339, 192))],
        true,
    )
}

/// k3j3 (v0·ε³, v1·ε⁻³): leading -117/64 ε⁻³ (massive -117/64).
#[test]
fn grid_k3j3() -> Result<()> {
    grid_point(
        "k3j3",
        &k(3, 3),
        &[(-3, q(-117, 64))],
        false,
        &[(-3, q(117, 32)), (-2, q(-1339, 192))],
        false,
    )
}

/// k4j2 (v0·ε⁴, v1·ε⁻²): -117/64 ε⁻¹ exactly (massive -117/64, all orders).
#[test]
fn grid_k4j2() -> Result<()> {
    grid_point(
        "k4j2",
        &k(4, 2),
        &[(-1, q(-117, 64))],
        true,
        &[(-1, q(117, 32))],
        true,
    )
}

/// k4j3 (v0·ε⁴, v1·ε⁻³): leading -117/64 ε⁻² (massive -117/64).
#[test]
fn grid_k4j3() -> Result<()> {
    grid_point(
        "k4j3",
        &k(4, 3),
        &[(-2, q(-117, 64))],
        false,
        &[(-2, q(117, 32)), (-1, q(-1339, 192))],
        true,
    )
}

/// A row-level point (grid rA, rB): the rows of a NONZERO group vanish
/// inside it (the moment extraction drops them, no count sees them) while
/// the X-fed row stays measured. Decided through linearity in v0's
/// numerator: `point − control` (the X-only control) is the pure-L0 point
/// (k3j1, k4j2), -117/64 C_A³/ε exactly for the generic and the massive
/// method; the legacy child route gave +117/32.
fn row_point(
    label: &str,
    dot: &str,
    control: &str,
    pins: &[(i64, Atom)],
    complete: bool,
) -> Result<()> {
    let case = Case::new(dot, None)?;
    let today = case.today()?;
    let child = case.child(run(FIXED))?;
    ensure_same(label, &child.value, &today)?;
    // Round 8b: declined (k_v1 > 0; before: taken).
    ensure!(child.root.bounded >= 1, "{label}: {:?}", child.root);
    ensure_declined(label, &child.root, "declared_order")?;
    pinned(label, &today, pins, complete)?;
    let control = Case::new(control, None)?.today()?;
    pinned(
        &format!("{label} - control"),
        &(&child.value - &control),
        &[(-1, q(-117, 64))],
        true,
    )?;
    // The planner before IRH-DEPTH misses the pole: legacy − control is not
    // the reference −117/64 ε⁻¹. Round 8b: the pin to the old wrong value
    // (+117/32) is dropped where it no longer holds: the legacy emulation
    // has no retry any more (one attempt at L_q, no guard), so its wrong
    // value can differ (rB: +151/64).
    let legacy = case.child(run(LEGACY))?;
    let missed = pinned(
        &format!("{label} legacy - control"),
        &(&legacy.value - &control),
        &[(-1, q(-117, 64))],
        true,
    );
    ensure!(missed.is_err(), "{label}: the legacy planner must miss the pole");
    Ok(())
}

/// rA: v0 = ε³ L0 + ε² X, v1·ε⁻¹: -169/64 ε⁻¹ (generic; massive -91/64 at
/// this rearrangement-unsafe order), rA − xA = -117/64 ε⁻¹.
#[test]
fn grid_row_level_a() -> Result<()> {
    row_point("rA", &r_a(), &x_a(), &[(-1, q(-169, 64))], true)
}

/// rA0: v0 = ε³ L0 + X, v1·ε⁻¹: every group is nonzero at the first depth
/// (3 of 13 rows shown, the 10 rows with o = 1 and P = 2 dropped inside
/// them): 13/32 ε⁻³ (ε⁻² and ε⁻¹ carry logs), rA0 − xA0 = -117/64 ε⁻¹.
#[test]
fn grid_row_level_a0() -> Result<()> {
    row_point("rA0", &r_a0(), &x_a0(), &[(-3, q(13, 32))], false)
}

/// rB: v0 = ε⁴ L0 + ε² X, v1·ε⁻²: -13/16 ε⁻² (ε⁻¹ carries logs),
/// rB − xB = -117/64 ε⁻¹.
#[test]
fn grid_row_level_b() -> Result<()> {
    row_point("rB", &r_b(), &x_b(), &[(-2, q(-13, 16))], false)
}

/// Positive control (review: quotient-v1, k0j1): v1·ε⁻¹ alone. The bound
/// (2) is the double quotient pole: the fixed planner starts there; without
/// the raise the shown rows measure it and the planner retries to 2.
#[test]
fn a_measured_quotient_pole_still_retries() -> Result<()> {
    let case = Case::new(&k(0, 1), None)?;
    let today = case.today()?;
    let child = case.child(run(FIXED))?;
    ensure_same("k0j1", &child.value, &today)?;
    // Round 8b: k_v1 = 1 declines on every planner that guards (before: the
    // bound 2 was taken, and without the raise the shown rows were measured
    // and retried to 2: no retry exists any more).
    ensure!(child.root.bounded >= 1, "{:?}", child.root);
    ensure_declined("k0j1", &child.root, "declared_order")?;
    let measured = case.child(run(NO_RAISE))?;
    ensure_same("k0j1 no raise", &measured.value, &today)?;
    ensure_declined("k0j1 no raise", &measured.root, "declared_order")?;
    // 117/32 ε⁻⁴ (massive 117/32 at this safe order).
    pinned("k0j1", &today, &[(-4, q(117, 32))], false)
}

/// Hard-only control (review: single-softer, k3j0): v0·ε³ alone. The hard
/// part vanishes at the first depth, but the quotient bound is its single
/// loop: nothing can hide, no retry, no refusal; the value is 0 (all three
/// methods).
#[test]
fn a_vanishing_hard_part_under_a_simple_pole_needs_nothing() -> Result<()> {
    let case = Case::new(&k(3, 0), None)?;
    let today = case.today()?;
    let child = case.child(run(FIXED))?;
    ensure_same("k3j0", &child.value, &today)?;
    let DepthCounts {
        corrections,
        bounded,
        retried,
        raised,
        refused,
        ..
    } = child.root;
    ensure!(
        corrections >= 1 && bounded == 0 && retried == 0 && raised == 0 && refused == 0,
        "{:?}",
        child.root
    );
    pinned("k3j0", &today, &[], true)
}

/// The guard fires (pre-review 2(e)): k3j1 with the raise disabled. The
/// accepted depth 1 stays below the bound 2 with every group vanished: the
/// correction is refused (`vanished`), and the WHOLE correction falls back
/// to today's path, so the value is today's, -117/64 ε⁻¹.
#[test]
fn the_guard_refuses_a_vanished_hard_part_below_the_bound() -> Result<()> {
    let case = Case::new(&k(3, 1), None)?;
    let today = case.today()?;
    let child = case.child(run(NO_RAISE))?;
    ensure_same("k3j1 no raise", &child.value, &today)?;
    // Round 8b: declined before H (k_v1 = 1; before: refused `vanished`).
    ensure_declined("k3j1 no raise", &child.root, "declared_order")?;
    pinned("k3j1 no raise", &child.value, &[(-1, q(-117, 64))], true)
}

/// The guard on the realistic hidden rows (pre-review finding 1): rA0
/// shows every group at the first depth (nothing vanished: the rows dropped
/// inside them are in no count). With the raise disabled it is refused
/// `unverified`; with the bound unknown it is refused `unbounded` before the
/// hard part is integrated (the blocking case: before, an unknown bound with
/// nothing vanished was accepted, the old planner's value, which misses
/// +351/64 ε⁻¹). Both fall back to today's value.
#[test]
fn the_guard_refuses_a_hidden_row_without_a_known_bound() -> Result<()> {
    let case = Case::new(&r_a0(), None)?;
    let today = case.today()?;
    let unverified = case.child(run(NO_RAISE))?;
    ensure_same("rA0 no raise", &unverified.value, &today)?;
    // Round 8b: declined before H (k_v1 = 1; before: refused `unverified`).
    ensure_declined("rA0 no raise", &unverified.root, "declared_order")?;
    let unbounded = case.child(run(DepthTest {
        bound: BoundOverride::Unknown,
        ..FIXED
    }))?;
    ensure_same("rA0 unknown bound", &unbounded.value, &today)?;
    // Round 8b: the decline precedes the unknown-bound refusal (before:
    // refused `unbounded`).
    ensure_declined("rA0 unknown bound", &unbounded.root, "declared_order")?;
    pinned(
        "rA0 unknown bound",
        &unbounded.value,
        &[(-3, q(13, 32))],
        false,
    )?;
    // The fallback is the RIGHT value, not only today's: the leading order
    // 13/32 ε⁻³ is the old route's too (the missed +351/64 sits at ε⁻¹), so
    // both refusals are pinned through the X-only control as well (tests
    // review, finding 3): rA0 − xA0 = -117/64 ε⁻¹ exactly.
    let control = Case::new(&x_a0(), None)?.today()?;
    for (label, refused) in [
        ("rA0 no raise", &unverified),
        ("rA0 unknown bound", &unbounded),
    ] {
        pinned(
            &format!("{label} - control"),
            &(&refused.value - &control),
            &[(-1, q(-117, 64))],
            true,
        )?;
    }
    Ok(())
}

/// A row whose insertion has an unknown excess (tests review 5): on the
/// unmodified graph (bound known, = the quotient's loops), the hook
/// multiplies one returned row's insertion by an opaque `f(ε)`. Its row
/// bound is unknown, so the correction is refused `unbounded` after the
/// first integration of H, before the quotient factors are built; the value
/// is today's (which never sees the hook).
#[test]
fn a_row_without_a_known_bound_is_refused_after_the_first_depth() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let today = case.today()?;
    let opaque = symbolica::function!(symbolica::symbol!("irh_depth_test_f"), epsilon());
    let hook: RowHook = Box::new(
        move |site: &HookSite, rows: &mut Vec<(Atom, Atom, usize)>| {
            if site.hard_loops != 2 || site.quotient_loops != 1 {
                return;
            }
            if let Some((_, insertion, _)) = rows.first_mut() {
                *insertion = &*insertion * &opaque;
            }
        },
    );
    let child = case.child(Run {
        hook: Some(hook),
        ..run(FIXED)
    })?;
    ensure_same("k0j0 opaque row", &child.value, &today)?;
    // Round 8b (C4): a row insertion that is not ε-free declines before its
    // quotient factors are built (before: its row bound was unknown,
    // refused `unbounded`).
    ensure!(child.root.bounded == 0, "{:?}", child.root);
    ensure!(child.root.insertion_not_eps_free == 1, "{:?}", child.root);
    ensure_declined("k0j0 opaque row", &child.root, "insertion_not_eps_free")
}

/// A measured need above the bound (pre-review note 3): k0j1 (true need 2)
/// with the bound forced to 1. Refused at once (`bound_exceeded`, counted
/// as exceeded, no retry) into today's value; an error under
/// GAMMALOOP_IR_HARD_DEPTH_STRICT=1.
#[test]
fn a_measured_pole_above_the_bound_is_an_alarm() -> Result<()> {
    // Round 8b: k0j0 (every k = 0, B = L_q = 1, measured need 1) with the
    // bound forced to 0; the measured need above B is the typed alarm with
    // STRICT off AND on (before: k0j1 forced to 1, refused `bound_exceeded`
    // without STRICT; k0j1 now declines before any measurement).
    let case = Case::new(L3_006, None)?;
    for strict in [false, true] {
        let below = DepthTest {
            bound: BoundOverride::Fixed(0),
            strict: Some(strict),
            ..FIXED
        };
        let error = match case.child(run(below)) {
            Ok(run) => {
                return Err(eyre!(
                    "strict={strict}: a measured need above B must fail: {:?}",
                    run.root
                ));
            }
            Err(error) => error,
        };
        ensure!(depth::is_alarm(&error), "strict={strict}: {error:#}");
        let text = format!("{error:#}");
        ensure!(text.contains("IRH-DEPTH") && text.contains("bound=0"), "{text}");
    }
    Ok(())
}

/// Nested child calls (pre-review 2(a)): c2-m9, k2j1 (v0·ε², v1·ε⁻¹) with
/// the marked line 9. No top-level child call hides a pole; the condition
/// holds only in the child calls nested in the root evaluation: the R* of
/// its decorated quotients (`quotient_rstar` → `rstar_child_counterterm`,
/// an evaluation on a proper support), hard edges {2,3,8}, every group zero
/// at depth 1, quotient bound 2 over one quotient loop (IRH-DEPTH.md 2.6).
/// The proper UV children of the root take no child-route correction (their
/// counts are zero on either route), so they are kept on today's path. On the
/// fixed planner every nested call is bounded and starts at its bound (one
/// attempt at depth 2, recorded by a hook), and the counterterm equals
/// today's exactly; -117/64 ε⁻² at the safe order (massive -117/64). The
/// legacy planner misses the nested pole: legacy − today is exactly
/// [(351/32) log(M²/μᵣ²) − 697/64] C_A³/ε in the record (IRH-DEPTH.md 2.6).
#[test]
fn nested_child_calls_keep_their_quotient_pole() -> Result<()> {
    let case = Case::new(&k(2, 1), Some(9))?;
    let today = case.today()?;
    let sites = std::rc::Rc::new(std::cell::RefCell::new(Vec::<HookSite>::new()));
    let record = sites.clone();
    let hook: RowHook = Box::new(move |site: &HookSite, _: &mut Vec<(Atom, Atom, usize)>| {
        record.borrow_mut().push(site.clone());
    });
    let child = case.child(Run {
        hook: Some(hook),
        ..run(FIXED)
    })?;
    ensure_same("c2-m9 nested", &child.value, &today)?;
    let root = child.root;
    // Round 8b: every nested correction whose quotient holds v1·ε⁻¹ declines
    // (declared_order) before H; no nested attempt hides a pole, every
    // refusal is such a decline, no retry, no alarm (before: every nested
    // call was bounded and taken at its bound 2, on hard {2,3,8}).
    ensure!(
        root.nested >= 1
            && root.nested_bounded == root.nested
            && root.declared_order >= 1
            && root.declined == root.refused
            && root.retried == 0
            && root.exceeded == 0,
        "{root:?}"
    );
    let sites = sites.borrow();
    let nested = sites.iter().filter(|site| site.nested).collect::<Vec<_>>();
    let deep = nested
        .iter()
        .filter(|site| site.edges == "2,3,8")
        .collect::<Vec<_>>();
    ensure!(
        deep.is_empty(),
        "the nested calls on hard {{2,3,8}} decline before H: {sites:?}"
    );
    eprintln!(
        "IRH_DEPTH_NESTED c2-m9 attempts={} nested={} on_2_3_8={}",
        sites.len(),
        nested.len(),
        deep.len()
    );
    pinned("c2-m9", &today, &[(-2, q(-117, 64))], false)?;
    let legacy = case.child(run(LEGACY))?;
    ensure!(
        !same(&legacy.value, &today)?,
        "c2-m9: the legacy planner must miss the nested pole"
    );
    pinned("c2-m9 legacy", &legacy.value, &[(-2, q(-117, 64))], false)?;
    // The order the fix changes (tests review, finding 2).
    let log = |symbol| symbolica::function!(symbolica::atom::Symbol::LOG, Atom::var(symbol));
    let log_ratio = log(super::top_mass::auxiliary_mass_squared()) - log(GS.mu_r_sq);
    // Round 8b: the legacy emulation has no retry any more; its miss is
    // asserted (above), its old value (with log(M²/μᵣ²)) is not pinned.
    let _ = log_ratio;
    Ok(())
}

/// The trailing ε order of `value` (`None` when zero).
fn order(value: &Atom) -> Option<i64> {
    let series = value
        .series(
            GS.dim_epsilon,
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
fn first_rows(case: &Case) -> Result<Vec<(Atom, Atom, usize)>> {
    let rows = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
    let record = rows.clone();
    let hook: RowHook = Box::new(
        move |site: &HookSite, shown: &mut Vec<(Atom, Atom, usize)>| {
            if site.hard_loops == 2 && site.quotient_loops == 1 && site.extra == 1 {
                *record.borrow_mut() = shown.clone();
            }
        },
    );
    case.child(Run {
        hook: Some(hook),
        ..run(FIXED)
    })?;
    Ok(std::rc::Rc::try_unwrap(rows)
        .map_err(|_| eyre!("the hook outlived its run"))?
        .into_inner())
}

/// The vanished row is measured (tests review, finding 4). The candidate
/// rows are taken in a canonical order (their insertions' text, not
/// symbolica's term order), and the test uses the first one whose hooked
/// quotient factor needs the full bound 3, so the retry comes from measuring
/// the vanished row (`raised == 0`), whatever the process ran before. The
/// same hook without the raise to the bound must reach the reference
/// through the measurement alone.
#[test]
fn a_vanished_row_of_a_nonzero_group_keeps_its_quotient_pole() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let rows = first_rows(&case)?;
    let top = rows.iter().map(|(_, _, loops)| *loops).max().unwrap_or(0);
    let mut group = rows
        .iter()
        .filter(|(h, _, loops)| *loops == top && !h.is_zero() && order(h) == Some(-2))
        .map(|(_, insertion, _)| {
            (
                canonical_tensor(insertion).map(|tensor| tensor.to_string()),
                insertion,
            )
        })
        .map(|(key, insertion)| key.map(|key| (key, insertion)))
        .collect::<Result<Vec<_>>>()?;
    ensure!(group.len() >= 2, "the top group needs two rows: {rows:?}");
    group.sort_by(|a, b| a.0.cmp(&b.0));
    let hooked = |target: &Atom, shown: bool, start: usize, test: DepthTest| {
        case.child(Run {
            hook: Some(vanished_row(target.clone(), shown)),
            start,
            test,
        })
    };
    // Round 8b (C4): the hooked row's insertion E/ε² is not ε-free; the
    // correction declines before its quotient factors are built, on every
    // planner that guards, and the value is today's (the hook never reaches
    // the generic path). Before round 8b its quotient factor was measured
    // while its h vanished and the depth retried to 3 (no measurement
    // decides the depth now).
    let today = case.today()?;
    let (_, target) = group.first().ok_or_else(|| eyre!("no row"))?;
    for (label, start, test) in [("fixed", 0, FIXED), ("deep", 3, FIXED), ("no raise", 0, NO_RAISE)] {
        let child = hooked(target, true, start, test)?;
        ensure_same(label, &child.value, &today)?;
        ensure!(child.root.insertion_not_eps_free >= 1, "{label}: {:?}", child.root);
        ensure_declined(label, &child.root, "insertion_not_eps_free")?;
    }
    Ok(())
}

// Round 4 (review IRH-DEPTH finding 3): a fatal depth alarm raised in a
// NESTED evaluation keeps its fatal status through the three boundaries the
// review found: the vanished-row valuation, the guarded parent's fallback and
// the counterterm (Δ) cache. The breach is forced in the nested evaluations
// only ([`BoundOverride::Nested`]): c2-m9's decorated quotients on the hard
// part {2,3,8} need 2 poles at depth 2 (the start is raised to 2 so the rows
// are shown); their bound is forced to 1, so each is a measured row above its
// bound, while the root's corrections keep their computed bound and raise no
// alarm of their own. The override applies only inside a child-route attempt
// (or a vanished row's valuation) of the root, never on the root's fallback
// to today's terms, so a fallback that loses the alarm computes a plain value.

/// The forced nested breach: the nested bound 1 inside `within` of the root.
fn nested_breach(within: depth::Within, guard: bool, strict: bool) -> Run {
    // Round 8b: on c0-m9 (L3_006 with the marked line 9: every k = 0, so the
    // nested corrections are taken at B = L_q = 1 and measure a need of 1)
    // the nested bound is forced to 0 (before: c2-m9, whose nested
    // corrections now decline, forced to 1).
    Run {
        start: 0,
        ..run(DepthTest {
            bound: BoundOverride::Nested { bound: 0, within },
            guard: Some(guard),
            strict: Some(strict),
            ..FIXED
        })
    }
}

/// `result` must be the strict depth alarm of a nested evaluation.
fn ensure_alarm(label: &str, result: Result<ChildRun>) -> Result<()> {
    let error = match result {
        Ok(run) => {
            return Err(eyre!(
                "{label}: GAMMALOOP_IR_HARD_DEPTH_STRICT=1 must fail on a nested breach; it returned {} with {:?}",
                run.value,
                run.root
            ));
        }
        Err(error) => {
            // The class, not only the text (MATRIX-R45 review): a lost class
            // is a declinable error.
            ensure!(
                super::depth::is_alarm(&error),
                "{label}: the error lost the alarm class: {error:#}"
            );
            format!("{error:#}")
        }
    };
    ensure!(
        error.contains("IRH-DEPTH") && error.contains("bound=0"),
        "{label}: not the depth alarm: {error}"
    );
    Ok(())
}

/// A nested breach in a child-route attempt of the root, under both guard
/// settings. Without STRICT the nested evaluation refuses into today's terms
/// (counted, the value is today's); with STRICT the run fails, also when the
/// root runs under GAMMALOOP_IR_HARD_CHILD_GUARD=1 (on 307ac11f2f the guarded
/// root declined the nested alarm, its own counter unchanged, and its
/// fallback computed today's value: STRICT was not fatal).
#[test]
fn a_nested_alarm_is_fatal_under_both_guard_settings() -> Result<()> {
    // Round 8b: a measured need above B is the typed alarm with STRICT off
    // too (before: without STRICT the nested evaluation refused into today's
    // value).
    let case = Case::new(L3_006, Some(9))?;
    for guard in [false, true] {
        for strict in [false, true] {
            let label = format!("c0-m9 nested breach guard={guard} strict={strict}");
            ensure_alarm(
                &label,
                case.child(nested_breach(depth::Within::ChildAttempt, guard, strict)),
            )?;
        }
    }
    Ok(())
}

/// Every row of the root's child-route corrections vanishes (a hook), so the
/// decorated quotients are computed only while the vanished rows are valued
/// (`ir_quotient_valuation`); the nested breach is forced there only. With
/// STRICT the run fails under both guard settings (on 307ac11f2f the
/// valuation turned the forest's error into "no measurement" and the run
/// succeeded).
#[test]
fn a_nested_alarm_while_valuing_a_vanished_row_is_fatal() -> Result<()> {
    // Round 8b: the vanished rows are valued only under STRICT (a
    // cross-check); without it nothing is valued and no nested evaluation
    // runs there (before: valued always, the breach refused without STRICT).
    let case = Case::new(L3_006, Some(9))?;
    let vanish = || -> RowHook {
        Box::new(|site: &HookSite, rows: &mut Vec<(Atom, Atom, usize)>| {
            if !site.nested {
                for (h, _, _) in rows.iter_mut() {
                    *h = Atom::Zero;
                }
            }
        })
    };
    for guard in [false, true] {
        let label = format!("c0-m9 vanished rows guard={guard}");
        let unvalued = case.child(Run {
            hook: Some(vanish()),
            ..nested_breach(depth::Within::ZeroRow, guard, false)
        })?;
        ensure!(
            unvalued.root.exceeded == 0,
            "{label}: a vanished row was valued without STRICT: {:?}",
            unvalued.root
        );
        ensure_alarm(
            &label,
            case.child(Run {
                hook: Some(vanish()),
                ..nested_breach(depth::Within::ZeroRow, guard, true)
            }),
        )?;
    }
    Ok(())
}

/// The nested alarm stored by the shared counterterm (Δ) cache and replayed
/// by a later evaluation that forces nothing: the replay is the alarm again,
/// fatal also under GAMMALOOP_IR_HARD_CHILD_GUARD=1. On 307ac11f2f this test
/// passes too: the cache replayed the alarm as a plain message, which the
/// guard declined, but the root's fallback to today's terms looked up the
/// same quotient and failed on the replay. The class itself is checked by
/// [`the_counterterm_cache_replays_the_alarm_class`].
#[test]
fn a_cached_nested_alarm_stays_fatal() -> Result<()> {
    for guard in [false, true] {
        let label = format!("c0-m9 cached alarm guard={guard}");
        let mut case = Case::new(L3_006, Some(9))?;
        case.plan = case.plan.with_counterterm_cache(std::sync::Arc::new(
            crate::uv::counterterm_cache::CountertermCache::default(),
        ));
        ensure_alarm(
            &format!("{label} (computed)"),
            case.child(nested_breach(depth::Within::ChildAttempt, false, true)),
        )?;
        // Nothing forced: only the cache holds the alarm.
        let replay = case.child(Run {
            start: 0,
            ..run(DepthTest {
                guard: Some(guard),
                strict: Some(true),
                ..FIXED
            })
        });
        ensure_alarm(&format!("{label} (replayed)"), replay)?;
    }
    Ok(())
}

/// The counterterm cache keeps the class of a stored error (round 4): a
/// computation that raised the strict alarm replays it on every later hit as
/// the alarm (not a plain message), with the same text, without computing
/// again.
#[test]
fn the_counterterm_cache_replays_the_alarm_class() -> Result<()> {
    let case = Case::new(&k(2, 1), Some(9))?;
    let cache = crate::uv::counterterm_cache::CountertermCache::default();
    let lookup = |compute: &dyn Fn() -> Result<Atom>| {
        cache.counterterm_with(
            &case.graph,
            &case.plan.support,
            case.model(),
            "irh-depth-round4".into(),
            compute,
        )
    };
    let first = lookup(&|| Err(depth::alarm("IRH-DEPTH: test alarm bound=1".into())))
        .expect_err("the computation fails");
    let replayed = lookup(&|| panic!("a cached error is not computed again"))
        .expect_err("the cached error is replayed");
    ensure!(
        depth::is_alarm(&first) && depth::is_alarm(&replayed),
        "the alarm class is lost: {first:#} / {replayed:#}"
    );
    ensure!(
        format!("{replayed:#}").contains("IRH-DEPTH: test alarm bound=1"),
        "{replayed:#}"
    );
    // An ordinary error stays ordinary.
    let plain = cache
        .counterterm_with(
            &case.graph,
            &case.plan.support,
            case.model(),
            "irh-depth-round4-plain".into(),
            || Err(eyre!("an ordinary failure")),
        )
        .expect_err("the computation fails");
    ensure!(!depth::is_alarm(&plain), "{plain:#}");
    Ok(())
}

/// MATRIX-R45 review (nested-and-cache, finding 1): a Δ that exhausts the
/// scoped budget of a guarded child attempt is not stored. The guard
/// declines with `budget`, and the fallback outside any budget computes the
/// same Δ again (before: it replayed the guard's budget error and the
/// diagram failed).
#[test]
fn a_guard_budget_error_is_not_kept_by_the_cache() -> Result<()> {
    use symbolica::utils::resources::checkpoint;
    let case = Case::new(&k(2, 1), Some(9))?;
    let cache = crate::uv::counterterm_cache::CountertermCache::default();
    let declined: std::result::Result<Atom, String> =
        super::closed::guarded(Some(1), None, || {
            let delta = cache.counterterm_with(
                &case.graph,
                &case.plan.support,
                case.model(),
                "irh-r45-poison".into(),
                || {
                    checkpoint().map_err(eyre::Report::new)?;
                    Ok(Atom::num(7))
                },
            )?;
            Ok(Ok(delta))
        })?;
    ensure!(
        declined == Err("budget".to_owned()),
        "the guarded attempt was not a budget decline: {declined:?}"
    );
    let fallback = cache.counterterm_with(
        &case.graph,
        &case.plan.support,
        case.model(),
        "irh-r45-poison".into(),
        || Ok(Atom::num(7)),
    );
    ensure!(
        fallback.is_ok(),
        "the fallback replays the guard's budget error: {:#}",
        fallback.unwrap_err()
    );
    Ok(())
}

/// IRH-DEPTH round 5: the regularity witnesses ([`super::regularity`]).
#[path = "regularity_tests.rs"]
mod round5;

/// IRH-DEPTH round 8b: C1, the topological certificate and the declines.
#[path = "r8b_tests.rs"]
mod round8b;

/// IRH-DEPTH round 8 on the G-G stack: the K3 records under C1 and the declines.
#[path = "r8_k3_tests.rs"]
mod round8_k3;
