//! Glue between a gammaloop one-loop graph and the `oneloop` IBP reducer.

use idenso::dirac::GammaSimplifier;
use idenso::shorthands::metric::MetricSimplifier;
use oneloop::masters::{MasterBasis, OneLoopMasters};
use spenso::network::library::symbolic::ETS;
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::structure::representation::{LibraryRep, Minkowski};
use symbolica::atom::{Atom, AtomCore};
use tracing::warn;

use crate::feyngen::diagram_generator::evaluate_overall_factor;
use crate::graph::Graph;
use crate::graph::lmb::LMBext;
use crate::utils::{GS, W_};
use crate::uv::uv_graph::UVE;

/// The outcome of attempting to reduce a graph's numerator
pub(crate) enum ReduceOutcome {
    /// Reduced to `Σ coeff × master`, formatted for the dot `reduced_num`.
    Reduced(String),
    /// Not a one-loop graph
    NotOneLoop(usize),
    /// One-loop, but the contracted numerator is identically zero.
    ZeroNumerator,
    /// One-loop with a non-zero numerator, but the reducer could not produce it.
    ///
    /// Covers both "the reducer returned no terms" and "the graph could not be faithfully
    /// translated into the reducer's conventions" (routing, external labelling, an
    /// uncontracted Lorentz index, ...). The second case logs the reason at WARN; it
    /// deliberately does *not* get its own status string, because the app-side consumer
    /// regex-matches the `not_one_loop` / `zero_numerator` / `unsupported` vocabulary.
    /// Surfacing the reason as a `reduce_error` dot attribute would be a strict improvement
    /// but needs a matching arm in `graph::parse::serialization`.
    Unsupported,
}

/// The gammaloop tensor heads the `oneloop` bridge pattern-matches on.
///
/// These come from the canonical accessors rather than string literals, so a spenso rename
/// is a compile error instead of a silent degradation of every numerator to a constant.
/// (`GS` has no `mink`/`g` fields; `LibraryRep::symbol()` and `ETS.metric` are the
/// authoritative sources -- see `spenso::structure::representation` and
/// `spenso::network::library::symbolic`.)
pub(crate) fn gammaloop_heads() -> oneloop::bridge::GammaloopHeads {
    oneloop::bridge::GammaloopHeads {
        loop_mom: GS.loop_mom,
        external_mom: GS.external_mom,
        index: LibraryRep::from(Minkowski {}).symbol(),
        metric: ETS.metric,
    }
}

/// The loop propagators of `graph`, with their masses, in `iter_edges()` order.
///
/// Split out from [`reduce_graph_numerator`] so the edge filter and the mass extraction can
/// be asserted on directly.
pub(crate) fn graph_edges(graph: &Graph) -> Vec<oneloop::bridge::GammaloopEdge> {
    let mut edges = Vec::new();
    for (_, i, _) in graph.iter_edges() {
        let loop_expr = graph
            .loop_momentum_basis
            .loop_atom(i, GS.loop_mom, &[W_.a___], false);
        if loop_expr.is_zero() {
            continue; // external / tree edge, not a loop propagator
        }
        let ext_expr = graph
            .loop_momentum_basis
            .ext_atom(i, GS.external_mom, &[W_.a___], false);
        // `UVE::mass_atom`, not the inherent `PossibleParticle::mass_atom`: only the former
        // maps the UFO `ZERO` symbol to `Atom::Zero`. Without it a massless propagator
        // arrives as `mass_sq = ZERO^2`, which is not `Atom::Zero`, so `tadpole_coefficient`'s
        // scaleless short-circuit never fires and every master argument carries a `ZERO^2`.
        let mass = UVE::mass_atom(&graph[i]);
        edges.push(oneloop::bridge::GammaloopEdge {
            lmb_rep: loop_expr + ext_expr,
            mass_sq: &mass * &mass,
        });
    }
    edges
}

/// Reduce a one-loop graph's numerator to master integrals.
pub(crate) fn reduce_graph_numerator(graph: &Graph, num: &Atom) -> ReduceOutcome {
    // The reducer only handles a single loop; tree/multi-loop short-circuit.
    let n_loops = graph.loop_momentum_basis.loop_edges.len();
    if n_loops != 1 {
        return ReduceOutcome::NotOneLoop(n_loops);
    }

    // Complete the Dirac traces to a scalar. `simplify_gamma` already collects
    // gamma chains itself, so no separate `collect_gamma_chains()` is needed
    // here (and the caller in `serialization` has gamma-simplified already).
    let scalar = num.simplify_gamma().expand().simplify_metrics();

    // Collapse the graph grouping / symmetry / sign bookkeeping symbols
    // (`NumeratorDependentGrouping`, `AutG`, `InternalFermionLoopSign`, …) into
    // their numeric values, so grouped-graph duplicate terms sum instead of each
    // monomial carrying a ~95-char tag and appearing once per grouped graph.
    let scalar = evaluate_overall_factor(scalar.as_view());

    // Rewrite the edge momenta `Q(edge, ..)` into the loop basis `K(0,..) + P(j,..)`.
    let reps =
        graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[W_.a___]);
    let num_lmb = scalar
        .replace_multiple(&reps)
        .expand()
        .simplify_metrics()
        .expand();

    // Collect the loop propagators with their masses.
    let edges = graph_edges(graph);
    if edges.is_empty() {
        // A one-loop graph should have loop propagators; be defensive. Reporting this as
        // `NotOneLoop(1)` would emit the contradictory pair
        // `reduce_status=not_one_loop, reduce_loops=1`.
        warn!("reduce: graph reports one loop but no loop propagators were found");
        return ReduceOutcome::Unsupported;
    }

    let family = match oneloop::bridge::family_from_gammaloop(&num_lmb, &edges, &gammaloop_heads())
    {
        Ok(family) => family,
        Err(e) => {
            warn!("reduce: cannot translate graph into the one-loop reducer's basis: {e}");
            return ReduceOutcome::Unsupported;
        }
    };
    let reduction = oneloop::reduce::reduce(&family);
    if reduction.terms.is_empty() {
        return if num_lmb.is_zero() {
            ReduceOutcome::ZeroNumerator
        } else {
            ReduceOutcome::Unsupported
        };
    }

    let show = |a: &Atom| a.printer(SpensoPrintSettings::typst_options()).to_string();

    // Fold the reduction into a single atom `Σ coeff · master`, then
    // `collect_factors` to pull the common coupling / colour / polarization
    // prefactor out front (so it appears once, not once per master) and collect
    // like terms in each coefficient. The master heads (`A0`/`B0`/`C0`/`D0`) are
    // opaque function symbols, so they survive the factoring untouched; strip the
    // `oneloop::` namespace for display.
    let basis = OneLoopMasters;
    let amplitude = reduction
        .terms
        .iter()
        .fold(Atom::Zero, |acc, (coeff, master)| {
            acc + coeff * basis.symbol(master)
        });
    let reduced_num = show(&amplitude.collect_factors())
        .replace("oneloop::", "")
        .replace('"', "\\\"");
    ReduceOutcome::Reduced(reduced_num)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dot;
    use crate::graph::parse::IntoGraph;
    use crate::initialisation::test_initialise;
    use crate::model::Model;
    use crate::utils::load_generic_model;
    use std::sync::LazyLock;
    use symbolica::atom::Symbol;
    use symbolica::{function, symbol};

    /// One shared model instance for the whole test binary. `load_generic_model` re-registers
    /// the UFO parameter symbols (with fresh print closures) on every call, which Symbolica
    /// rejects as a redefinition once another test thread has already done it.
    static SCALARS: LazyLock<Model> = LazyLock::new(|| load_generic_model("scalars"));

    /// A one-loop triangle of the massless (`mass = "ZERO"`) pdg-1000 scalar, the same shape
    /// `graph::parse::tests::test_load` prints. Its LMB routes the loop propagators as
    /// `K(0)`, `-P(0)+K(0)`, `-P(0)-P(1)+K(0)`.
    fn one_loop_triangle() -> Graph {
        dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v4
            ext -> v5
            v6 -> ext
            v5 -> v4;
            v6 -> v5;
            v4 -> v6;
        }, &*SCALARS)
        .unwrap()
    }

    #[test]
    fn heads_resolve_to_the_spenso_symbols_the_bridge_expects() {
        test_initialise().unwrap();
        // `reduce_bridge` used to build these from the raw string literals
        // `symbol!("spenso::mink")` / `symbol!("spenso::g")`. Pinning the canonical accessors
        // against those names means a spenso rename breaks *here* instead of silently
        // degrading every numerator to a constant.
        let heads = gammaloop_heads();
        assert_eq!(heads.index, symbol!("spenso::mink"));
        assert_eq!(heads.metric, symbol!("spenso::g"));
        assert_eq!(heads.loop_mom, GS.loop_mom);
        assert_eq!(heads.external_mom, GS.external_mom);
    }

    #[test]
    fn graph_edges_keeps_only_loop_propagators() {
        test_initialise().unwrap();
        let g = one_loop_triangle();
        let edges = graph_edges(&g);
        // 6 edges in the graph, 3 of them external/tree.
        assert_eq!(g.iter_edges().count(), 6);
        assert_eq!(edges.len(), 3);
        let reps: Vec<String> = edges.iter().map(|e| e.lmb_rep.to_string()).collect();
        assert!(
            reps.iter().all(|r| r.contains("K(0")),
            "every kept edge must carry the loop momentum: {reps:?}"
        );
    }

    #[test]
    fn massless_propagators_reach_the_reducer_as_exact_zero() {
        test_initialise().unwrap();
        // pdg 1000 has `"mass": "ZERO"`. The inherent `PossibleParticle::mass_atom` returns
        // `Atom::var(ZERO)`; only `UVE::mass_atom` maps it to `Atom::Zero`. With the wrong one
        // every propagator arrives as `mass_sq = ZERO^2`, which is not `Atom::Zero`, so
        // `tadpole_coefficient`'s scaleless short-circuit never fires and every master
        // argument carries a spurious `ZERO^2`.
        for e in graph_edges(&one_loop_triangle()) {
            assert_eq!(
                e.mass_sq,
                Atom::Zero,
                "massless propagator reached the reducer as `{}`",
                e.mass_sq
            );
        }
    }

    #[test]
    fn one_loop_triangle_reduces_to_a_triangle_master() {
        test_initialise().unwrap();
        let g = one_loop_triangle();
        match reduce_graph_numerator(&g, &Atom::num(1)) {
            ReduceOutcome::Reduced(s) => {
                assert!(s.contains("C0"), "expected a C0 master, got `{s}`");
                assert!(
                    !s.contains("ZERO"),
                    "the UFO `ZERO` symbol leaked into the masters: `{s}`"
                );
            }
            other => panic!("expected Reduced, got {}", outcome_name(&other)),
        }
    }

    #[test]
    fn a_zero_numerator_is_reported_as_such() {
        test_initialise().unwrap();
        let g = one_loop_triangle();
        assert!(matches!(
            reduce_graph_numerator(&g, &Atom::Zero),
            ReduceOutcome::ZeroNumerator
        ));
    }

    #[test]
    fn a_tree_graph_is_not_one_loop() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v1
            ext -> v1
            v1 -> ext
        }, &*SCALARS)
        .unwrap();
        match reduce_graph_numerator(&g, &Atom::num(1)) {
            ReduceOutcome::NotOneLoop(n) => assert_eq!(n, 0),
            other => panic!("expected NotOneLoop(0), got {}", outcome_name(&other)),
        }
    }

    #[test]
    fn a_two_loop_graph_is_not_one_loop() {
        test_initialise().unwrap();
        // Kite: four 3-valent vertices, five internal lines => 5 - 4 + 1 = 2 loops.
        let g: Graph = dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v1
            v3 -> ext
            v1 -> v2;
            v2 -> v3;
            v3 -> v4;
            v4 -> v1;
            v2 -> v4;
        }, &*SCALARS)
        .unwrap();
        match reduce_graph_numerator(&g, &Atom::num(1)) {
            ReduceOutcome::NotOneLoop(n) => assert_eq!(n, 2),
            other => panic!("expected NotOneLoop(2), got {}", outcome_name(&other)),
        }
    }

    #[test]
    fn dot_serialize_with_reduce_emits_a_reduced_numerator() {
        test_initialise().unwrap();
        // The actual user entry point: `save dot ... --reduce`. Nothing exercised this path
        // end to end before, so `reduce_status` / `reduced_num` had no test at all.
        use crate::processes::DotExportSettings;
        let g = one_loop_triangle();
        let dot = g.dot_serialize(&DotExportSettings {
            output_full_numerator: true,
            reduce: true,
            ..Default::default()
        });
        assert!(
            dot.contains("reduced_num"),
            "no `reduced_num` in the emitted dot:\n{dot}"
        );
        assert!(
            !dot.contains("reduce_status"),
            "a one-loop scalar triangle must reduce, but got a status:\n{dot}"
        );
        assert!(dot.contains("C0"), "expected a C0 master:\n{dot}");
        assert!(
            !dot.contains("ZERO"),
            "the UFO `ZERO` symbol leaked into the emitted dot:\n{dot}"
        );
    }

    #[test]
    fn dot_serialize_with_reduce_reports_a_two_loop_graph() {
        test_initialise().unwrap();
        use crate::processes::DotExportSettings;
        let g: Graph = dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v1
            v3 -> ext
            v1 -> v2;
            v2 -> v3;
            v3 -> v4;
            v4 -> v1;
            v2 -> v4;
        }, &*SCALARS)
        .unwrap();
        let dot = g.dot_serialize(&DotExportSettings {
            output_full_numerator: true,
            reduce: true,
            ..Default::default()
        });
        assert!(dot.contains("reduce_status"), "no status:\n{dot}");
        assert!(dot.contains("not_one_loop"), "wrong status:\n{dot}");
        assert!(dot.contains("reduce_loops"), "no loop count:\n{dot}");
        assert!(!dot.contains("reduced_num"), "unexpected reduction:\n{dot}");
    }

    /// `K(0, mink(4, idx))` / `P(j, mink(4, idx))` in the exact form gammaloop emits.
    fn mink_idx(i: i64) -> Atom {
        function!(symbol!("spenso::mink"), Atom::num(4), Atom::num(i))
    }

    #[test]
    fn graph_edges_route_the_triangle_as_the_reducers_chain_with_flipped_signs() {
        test_initialise().unwrap();
        // The premise the whole q-basis fix rests on, pinned against the *real* LMB rather
        // than against a hand-built `GammaloopEdge`. gammaloop routes the loop propagators as
        // `k`, `k - P(0)`, `k - P(0) - P(1)`, i.e. `r_a - r_{a-1} = -q_a`, whereas the reducer
        // hard-codes `r_i = +q1 + ... + q_{i-1}` (`reduce.rs`'s `triangle_topo` r_coeffs
        // `[[0,0,0],[1,0,0],[1,1,0]]` and its RSP rule `k.q1 = (D2-D1-m1+m2-p1)/2`). If the
        // LMB ever re-routes, this test says so instead of the bridge silently rejecting the
        // graph as `unsupported`.
        let a = Atom::var(W_.a___);
        let k = function!(GS.loop_mom, Atom::num(0), a.clone());
        let p = |j: i64| function!(GS.external_mom, Atom::num(j), a.clone());
        let want = [k.clone(), &k - &p(0), &k - &p(0) - &p(1)];
        let got: Vec<Atom> = graph_edges(&one_loop_triangle())
            .iter()
            .map(|e| e.lmb_rep.clone())
            .collect();
        assert_eq!(got, want.to_vec(), "gammaloop's one-loop routing changed");
    }

    #[test]
    fn an_odd_rank_numerator_on_the_real_graph_picks_up_the_routing_sign() {
        test_initialise().unwrap();
        // The joint the oneloop-side tests cannot close on their own: *real* graph edges +
        // an odd-rank numerator. `k.P(0)` is `-dot(k, q1)` in the reducer's chain basis,
        // because gammaloop's chain direction is `-P(0)`. Before the relabelling fix this
        // came out as `+dot(k, q1)`, i.e. every odd-rank monomial was negated -- invisible to
        // the scalar and rank-2 validation because both are sign-blind.
        let edges = graph_edges(&one_loop_triangle());
        let num = &function!(GS.loop_mom, Atom::num(0), mink_idx(77))
            * &function!(GS.external_mom, Atom::num(0), mink_idx(77));
        let fam = oneloop::bridge::family_from_gammaloop(&num, &edges, &gammaloop_heads())
            .expect("the real one-loop triangle must translate");
        let dot =
            |a: Symbol, b: Symbol| function!(oneloop::symbols::S.dot, Atom::var(a), Atom::var(b));
        let s = &oneloop::symbols::S;
        assert_eq!(fam.numerator, -dot(s.k, s.q1));
        // and the C(3,2) invariants in the lexicographic order `reduce_core` reads them in:
        // `(r1-r2)^2 = q1^2`, `(r1-r3)^2 = (q1+q2)^2`, `(r2-r3)^2 = q2^2`.
        assert_eq!(
            fam.kinematics
                .invariants
                .iter()
                .map(|i| i.expand())
                .collect::<Vec<_>>(),
            vec![
                dot(s.q1, s.q1),
                (dot(s.q1, s.q1) + Atom::num(2) * dot(s.q1, s.q2) + dot(s.q2, s.q2)).expand(),
                dot(s.q2, s.q2),
            ]
        );
    }

    /// A one-loop box and a one-loop pentagon of the massless pdg-1000 scalar.
    fn one_loop_box() -> Graph {
        dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v1
            ext -> v2
            v3 -> ext
            v4 -> ext
            v1 -> v2;
            v2 -> v3;
            v3 -> v4;
            v4 -> v1;
        }, &*SCALARS)
        .unwrap()
    }

    fn one_loop_pentagon() -> Graph {
        dot!(digraph {
            graph [ overall_factor = 1; ]
            edge [ pdg=1000 ]
            ext [style=invis]
            ext -> v1
            ext -> v2
            ext -> v3
            v4 -> ext
            v5 -> ext
            v1 -> v2;
            v2 -> v3;
            v3 -> v4;
            v4 -> v5;
            v5 -> v1;
        }, &*SCALARS)
        .unwrap()
    }

    #[test]
    fn a_real_box_and_pentagon_still_reduce_despite_their_lmb_routing() {
        test_initialise().unwrap();
        // REGRESSION. gammaloop's LMB does not hand the loop propagators over in loop order,
        // and its offsets are not the reducer's chain even up to sign: a real box arrives as
        //
        //     r = [ -q2 + q3,  q3,  0,  -q1 - q2 + q3 ]
        //
        // because exactly one propagator is the LMB basis edge and one polygon leg is the
        // *dependent* external, expanded by momentum conservation into a sum of the others.
        // Only the triangle happens to come out already in chain form. Demanding the chain
        // outright therefore sent every real box and pentagon to `unsupported`, even though
        // their scalar reduction is routing-independent and was correct before the guard.
        // `oneloop::bridge::chain_order` reorders instead.
        for (name, g, master) in [
            ("box", one_loop_box(), "D0"),
            ("pentagon", one_loop_pentagon(), "D0"),
        ] {
            match reduce_graph_numerator(&g, &Atom::num(1)) {
                ReduceOutcome::Reduced(s) => assert!(
                    s.contains(master),
                    "{name}: expected a {master} master, got `{s}`"
                ),
                other => panic!("{name}: expected Reduced, got {}", outcome_name(&other)),
            }
        }
    }

    #[test]
    fn a_reordered_box_keeps_the_invariants_the_reducer_reads_lexicographically() {
        test_initialise().unwrap();
        // The reorder is only faithful if the pairwise invariants are recomputed from the
        // *reordered* offsets: `reduce_core` reads `invariants` as the C(4,2) lexicographic
        // pairs of the chain `r_i = q1 + ... + q_{i-1}`, so a permutation of the propagators
        // has to permute them too. Chain order for this box is `[0, q3, -q2+q3, -q1-q2+q3]`,
        // whose lexicographic pairwise squares are the six entries below.
        let fam = oneloop::bridge::family_from_gammaloop(
            &Atom::num(1),
            &graph_edges(&one_loop_box()),
            &gammaloop_heads(),
        )
        .expect("a real one-loop box must translate");
        let s = &oneloop::symbols::S;
        let d = |a: Symbol, b: Symbol| function!(s.dot, Atom::var(a), Atom::var(b));
        let two = Atom::num(2);
        // r1=0, r2=q3, r3=-q2+q3, r4=-q1-q2+q3
        let want = [
            d(s.q3, s.q3),                                                   // (r1-r2)^2 = q3^2
            (d(s.q2, s.q2) - &two * d(s.q2, s.q3) + d(s.q3, s.q3)).expand(), // (r1-r3)^2
            (d(s.q1, s.q1) + &two * d(s.q1, s.q2) - &two * d(s.q1, s.q3) + d(s.q2, s.q2)
                - &two * d(s.q2, s.q3)
                + d(s.q3, s.q3))
            .expand(), // (r1-r4)^2
            d(s.q2, s.q2),                                                   // (r2-r3)^2 = q2^2
            (d(s.q1, s.q1) + &two * d(s.q1, s.q2) + d(s.q2, s.q2)).expand(), // (r2-r4)^2
            d(s.q1, s.q1),                                                   // (r3-r4)^2 = q1^2
        ];
        assert_eq!(
            fam.kinematics
                .invariants
                .iter()
                .map(|i| i.expand())
                .collect::<Vec<_>>(),
            want.to_vec()
        );
        assert_eq!(fam.propagators.len(), 4);
    }

    fn outcome_name(o: &ReduceOutcome) -> String {
        match o {
            ReduceOutcome::Reduced(s) => format!("Reduced({s})"),
            ReduceOutcome::NotOneLoop(n) => format!("NotOneLoop({n})"),
            ReduceOutcome::ZeroNumerator => "ZeroNumerator".into(),
            ReduceOutcome::Unsupported => "Unsupported".into(),
        }
    }
}
