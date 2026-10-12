use std::{
    ops::ControlFlow,
    sync::{
        atomic::{AtomicBool, AtomicUsize, Ordering},
        Arc,
    },
};

use linnet::half_edge::generate::{
    batch::{Event, Hook, HookMask, Output},
    Bridges, GenerationError, Legs, Plan, Query, Run, Structure, Stub, Theory,
};
use num_bigint::BigUint;
use num_rational::BigRational;
use num_traits::Zero;

type TestPlan = Plan<char, u8>;

fn phi(rules: &[(char, usize)]) -> Theory<char> {
    rules
        .iter()
        .fold(Theory::builder(), |builder, &(field, degree)| {
            builder.vertex(vec![Stub::undirected(field); degree], [])
        })
        .build()
        .unwrap()
}

#[test]
fn deferred_hooks_preserve_diagrams_counts_and_dispatch_order_across_threads() {
    let theory = phi(&[('p', 3), ('p', 3)]);
    let query = Query::<char, u8>::new(&Legs::new())
        .loops(3..=3)
        .structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        });
    let plain = Plan::new(&theory, &query).unwrap();
    let scalar = Run::default();
    let expected: Vec<_> = plain
        .diagrams(&scalar)
        .unwrap()
        .map(|d| d.canonical_form().clone())
        .collect();
    let expected_counts = plain.count(&scalar).unwrap();
    let plan = Plan::new(
        &theory,
        &query.deferred_hooks(HookMask {
            skeletons: true,
            dressings: true,
            select: true,
        }),
    )
    .unwrap();
    assert!(matches!(
        plan.diagrams(&scalar),
        Err(GenerationError::DeferredHooks)
    ));
    let mut reference = None;
    for threads in [1, 2, 4] {
        let observed = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let run = Run::default().parallel(true);
                let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
                let mut forms = Vec::new();
                let mut dispatch = Vec::new();
                let mut widest = 0;
                loop {
                    match driver.poll().unwrap() {
                        Event::Batch(requests) => {
                            widest = widest.max(requests.len());
                            let mut wave = Vec::new();
                            for request in &requests {
                                let shape = match request.hook() {
                                    Hook::Skeleton => {
                                        let v = request.skeleton().unwrap();
                                        (v.n_vertices(), v.edges().count(), v.is_complete())
                                    }
                                    Hook::Dressing => {
                                        let v = request.dressing().unwrap();
                                        (
                                            v.n_vertices(),
                                            v.edges().filter(|e| e.stub.is_some()).count(),
                                            v.is_complete(),
                                        )
                                    }
                                    Hook::Select => {
                                        let v = request.diagram().unwrap();
                                        (v.n_vertices(), v.edges().count(), v.is_complete())
                                    }
                                };
                                wave.push((request.hook(), shape));
                            }
                            dispatch.push(wave);
                            driver
                                .resolve(
                                    &requests.iter().map(|r| r.reply(true)).collect::<Vec<_>>(),
                                )
                                .unwrap();
                        }
                        Event::Diagram(diagram) => forms.push(diagram.canonical_form().clone()),
                        Event::Yield => {}
                        Event::Done => break,
                        _ => panic!("unexpected output mode"),
                    }
                }
                assert_eq!(forms, expected);
                assert_eq!(driver.counts().unwrap(), &expected_counts);
                assert!(
                    widest > 1,
                    "one worker must also supply multiple independent hook requests"
                );
                dispatch
            });
        if let Some(reference) = &reference {
            assert_eq!(&observed, reference);
        } else {
            reference = Some(observed);
        }
    }
    let mut driver = plan.batched(&scalar, Output::Count).unwrap();
    loop {
        match driver.poll().unwrap() {
            Event::Batch(requests) => driver
                .resolve(&requests.iter().map(|r| r.reply(true)).collect::<Vec<_>>())
                .unwrap(),
            Event::Yield => {}
            Event::Done => break,
            _ => panic!("counting must not construct output diagrams"),
        }
    }
    assert_eq!(driver.counts().unwrap(), &expected_counts);
}

#[test]
fn deferred_progress_stops_count_at_first_output_and_rejects_stale_tickets() {
    let theory = phi(&[('p', 3), ('p', 3)]);
    let plan = vacuum(&theory, 3);
    let run = Run::default().parallel(true).deferred_progress(true);
    let mut driver = plan.batched(&run, Output::Count).unwrap();
    loop {
        match driver.poll().unwrap() {
            Event::Progress(p) if !p.diagrams.is_zero() => {
                assert_eq!(p.diagrams, BigUint::from(1u8));
                driver.stop_after_progress();
                break;
            }
            Event::Progress(_) => {}
            Event::Yield => {}
            _ => panic!("first accepted output must notify progress"),
        }
    }
    assert_eq!(driver.counts().unwrap().diagrams, BigUint::from(1u8));
    assert!(!driver.counts().unwrap().complete);
    assert!(matches!(driver.poll().unwrap(), Event::Yield));
    assert!(matches!(driver.poll().unwrap(), Event::Done));

    let query = Query::<char, u8>::new(&Legs::new())
        .loops(3..=3)
        .structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        })
        .deferred_hooks(HookMask {
            skeletons: true,
            ..HookMask::default()
        });
    let plan = Plan::new(&theory, &query).unwrap();
    let run = Run::default();
    let mut old = plan.batched(&run, Output::Count).unwrap();
    let old_replies = loop {
        match old.poll().unwrap() {
            Event::Batch(requests) => {
                break requests.iter().map(|r| r.reply(true)).collect::<Vec<_>>()
            }
            Event::Yield => {}
            _ => panic!("expected a callback batch"),
        }
    };
    old.abort();
    let mut fresh = plan.batched(&run, Output::Count).unwrap();
    loop {
        match fresh.poll().unwrap() {
            Event::Batch(_) => break,
            Event::Yield => {}
            _ => panic!("expected callback"),
        }
    }
    assert!(matches!(
        fresh.resolve(&old_replies),
        Err(GenerationError::InvalidReply)
    ));
    assert!(fresh.aborted());
    let mut unstarted = plan.batched(&run, Output::Count).unwrap();
    assert!(matches!(
        unstarted.resolve(&[]),
        Err(GenerationError::InvalidReply)
    ));
}

#[test]
fn batch_progress_observes_rejected_outputs_and_delivers_its_last_accepted_item() {
    let theory = phi(&[('p', 3), ('p', 3)]);
    let query = Query::<char, u8>::new(&Legs::new())
        .loops(3..=3)
        .structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        });
    let plan = Plan::new(&theory, &query.clone().select(|_| false)).unwrap();
    let run = Run::default().progress(|p| {
        assert!(p.diagrams.is_zero());
        assert!(!p.topologies.is_zero());
        ControlFlow::Break(())
    });
    let mut driver = plan.batched(&run, Output::Count).unwrap();
    while !matches!(driver.poll().unwrap(), Event::Done) {}
    assert_eq!(driver.counts().unwrap().skeletons, BigUint::from(1u8));
    assert!(driver.counts().unwrap().diagrams.is_zero());
    assert!(!driver.counts().unwrap().complete);

    let plan = Plan::new(&theory, &query).unwrap();
    let run = Run::default().progress(|p| {
        if p.diagrams == BigUint::from(3u8) {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    });
    let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
    let mut diagrams = 0;
    loop {
        match driver.poll().unwrap() {
            Event::Diagram(_) => diagrams += 1,
            Event::Done => break,
            Event::Yield => {}
            _ => panic!("unexpected event"),
        }
    }
    assert_eq!(diagrams, 3);
    assert_eq!(driver.counts().unwrap().diagrams, BigUint::from(3u8));
    assert!(!driver.counts().unwrap().complete);
}

fn vacuum(theory: &Theory<char>, loops: usize) -> TestPlan {
    Plan::new(
        theory,
        &Query::new(&Legs::new())
            .loops(loops..=loops)
            .structure(Structure {
                self_loops: true,
                zero_flow_edges: true,
                ..Structure::default()
            }),
    )
    .unwrap()
}

#[test]
fn count_and_skeleton_dressings_keep_distinct_rules_and_native_certificates() {
    // Two identical cubic declarations give 2^V times the scalar inverse-aut sum.
    // This exact result is also asserted by the pre-redesign scalar oracle.
    let plan = vacuum(&phi(&[('p', 3), ('p', 3)]), 3);
    let mut reference = None;
    for threads in [1, 4] {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                let run = Run::default().parallel(true);
                let diagrams: Vec<_> = plan.diagrams(&run).unwrap().collect();
                let counts = plan.count(&run).unwrap();
                assert!(counts.complete);
                assert_eq!(counts.diagrams, BigUint::from(diagrams.len()));
                assert_eq!(counts.sum_inv_aut, BigRational::from_integer(5.into()));
                assert_eq!(counts.per_loop[&3].sum_inv_aut, counts.sum_inv_aut);
                let skeletons: Vec<_> = plan.skeletons(&run).unwrap().collect();
                assert_eq!(counts.skeletons, BigUint::from(skeletons.len()));
                for skeleton in &skeletons {
                    let native =
                        skeleton
                            .graph()
                            .canonize_by(Clone::clone, Clone::clone, Clone::clone);
                    assert_eq!(skeleton.automorphisms(), native.automorphism_count());
                    assert_eq!(skeleton.id(), native.canonical_form());
                    let first: Vec<_> = skeleton
                        .diagrams()
                        .map(|d| d.canonical_form().clone())
                        .collect();
                    let second: Vec<_> = skeleton
                        .diagrams()
                        .map(|d| d.canonical_form().clone())
                        .collect();
                    assert_eq!(first, second);
                }
                let flattened: Vec<_> = skeletons
                    .iter()
                    .flat_map(|s| s.diagrams())
                    .map(|d| d.canonical_form().clone())
                    .collect();
                let forms: Vec<_> = diagrams
                    .iter()
                    .map(|d| d.canonical_form().clone())
                    .collect();
                assert_eq!(forms, flattened);
                if let Some(reference) = &reference {
                    assert_eq!(&forms, reference);
                } else {
                    reference = Some(forms);
                }
            });
    }
}

#[test]
fn skeletons_without_dressings_and_arbitrary_selection_remain_visible() {
    let theory = phi(&[('p', 3), ('x', 4)]);
    let query = Query::<char, u8>::new(&Legs::new())
        .loops(3..=3)
        .structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        });
    let plan = Plan::new(&theory, &query).unwrap();
    let selected = Plan::new(&theory, &query.clone().select(|_| false)).unwrap();
    let run = Run::default();
    let skeletons: Vec<_> = plan.skeletons(&run).unwrap().collect();
    assert!(skeletons.iter().any(|s| s.diagrams().next().is_none()));
    let counts = plan.count(&run).unwrap();
    let rejected = selected.count(&run).unwrap();
    assert_eq!(counts.skeletons, BigUint::from(skeletons.len()));
    assert_eq!(rejected.skeletons, counts.skeletons);
    assert!(rejected.diagrams.is_zero());
    assert!(rejected.sum_inv_aut.is_zero());
    assert!(rejected.complete);
    assert_eq!(selected.skeletons(&run).unwrap().count(), skeletons.len());
}

#[test]
fn retained_skeleton_dresses_after_its_execution_and_rejects_another_plan() {
    let theory = phi(&[('p', 3), ('p', 3)]);
    let plan = vacuum(&theory, 3);
    let (snapshot, reference) = {
        let run = Run::default();
        let skeleton = plan.skeletons(&run).unwrap().next().unwrap();
        let reference: Vec<_> = skeleton
            .diagrams()
            .map(|d| d.canonical_form().clone())
            .collect();
        (skeleton.snapshot(), reference)
    };
    let run = Run::default();
    let retained: Vec<_> = plan
        .dressings(&snapshot, &run)
        .unwrap()
        .map(|d| d.canonical_form().clone())
        .collect();
    assert_eq!(retained, reference);
    let other = vacuum(&theory, 3);
    assert!(matches!(
        other.dressings(&snapshot, &run),
        Err(GenerationError::ForeignSkeleton)
    ));
    let id = snapshot.id().clone();
    drop(plan);
    drop(theory);
    assert_eq!(snapshot.id(), &id);
    assert!(snapshot.automorphisms().is_some());
}

#[test]
fn abort_stops_decorating_and_buffered_outputs_after_the_first_item() {
    let plan = vacuum(&phi(&[('p', 3), ('p', 3)]), 3);
    for parallel in [false, true] {
        let abort = AtomicBool::new(false);
        let run = Run::default().parallel(parallel).abort(&abort);
        let mut diagrams = plan.diagrams(&run).unwrap();
        assert!(diagrams.next().is_some());
        abort.store(true, Ordering::Relaxed);
        assert!(diagrams.next().is_none());
        assert!(diagrams.aborted());
        let counts = plan.count(&run).unwrap();
        assert!(!counts.complete);
        assert!(counts.diagrams.is_zero());
    }
    let run = Run::default().progress(|progress| {
        if progress.diagrams >= BigUint::from(1u8) {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    });
    let counts = plan.count(&run).unwrap();
    assert_eq!(counts.diagrams, BigUint::from(1u8));
    assert!(!counts.complete);
}

#[test]
fn abort_inside_skeleton_hook_prevents_later_calls_and_outputs() {
    // One worker keeps cancellation deterministic while exercising both native routes;
    // wasm uses the sequential route even when parallel execution is requested.
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            let theory = phi(&[('p', 3), ('p', 3)]);
            for parallel in [false, true] {
                for keep in [false, true] {
                    let abort = Arc::new(AtomicBool::new(false));
                    let invoked = Arc::new(AtomicUsize::new(0));
                    let flag = Arc::clone(&abort);
                    let calls = Arc::clone(&invoked);
                    let query = Query::<char, u8>::new(&Legs::new())
                        .loops(3..=3)
                        .structure(Structure {
                            self_loops: true,
                            zero_flow_edges: true,
                            ..Structure::default()
                        })
                        .prune_skeletons(move |_| {
                            calls.fetch_add(1, Ordering::Relaxed);
                            flag.store(true, Ordering::Relaxed);
                            keep
                        });
                    let plan = Plan::new(&theory, &query).unwrap();
                    let run = Run::default().parallel(parallel).abort(&abort);
                    let counts = plan.count(&run).unwrap();
                    assert!(!counts.complete);
                    assert!(counts.diagrams.is_zero());
                    assert_eq!(invoked.load(Ordering::Relaxed), 1);

                    abort.store(false, Ordering::Relaxed);
                    invoked.store(0, Ordering::Relaxed);
                    let mut diagrams = plan.diagrams(&run).unwrap();
                    assert!(diagrams.next().is_none());
                    assert!(diagrams.aborted());
                    assert_eq!(invoked.load(Ordering::Relaxed), 1);
                }
            }
        });
}

// Each entry point must observe cancellation inside the predicate before another call or output.
fn assert_predicate_abort(plan: &TestPlan, abort: &AtomicBool, calls: &AtomicUsize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| {
            for parallel in [false, true] {
                let run = Run::default()
                    .parallel(parallel)
                    .check_heredity(true)
                    .abort(abort);
                abort.store(false, Ordering::Relaxed);
                calls.store(0, Ordering::Relaxed);
                let counts = plan.count(&run).unwrap();
                assert!(!counts.complete);
                assert!(counts.diagrams.is_zero());
                assert_eq!(calls.load(Ordering::Relaxed), 1);

                abort.store(false, Ordering::Relaxed);
                calls.store(0, Ordering::Relaxed);
                let mut diagrams = plan.diagrams(&run).unwrap();
                assert!(diagrams.next().is_none());
                assert!(diagrams.aborted());
                assert_eq!(calls.load(Ordering::Relaxed), 1);
            }
        });
    for owned in [false, true] {
        abort.store(false, Ordering::Relaxed);
        calls.store(0, Ordering::Relaxed);
        let run = Run::default().check_heredity(true).abort(abort);
        let skeleton = plan.skeletons(&run).unwrap().next().unwrap();
        let snapshot = skeleton.snapshot();
        let mut dressings = if owned {
            plan.dressings(&snapshot, &run).unwrap()
        } else {
            skeleton.diagrams()
        };
        assert!(dressings.next().is_none());
        assert!(dressings.aborted());
        assert_eq!(calls.load(Ordering::Relaxed), 1);
    }
}

#[test]
fn abort_inside_partial_dressing_stops_each_public_output_route() {
    let theory = phi(&[('p', 3), ('q', 3)]);
    let bare = Query::<char, u8>::new(&Legs::new()).loops(2..=2);
    let control = Plan::new(&theory, &bare)
        .unwrap()
        .count(&Run::default().parallel(false))
        .unwrap();
    assert!(control.complete && !control.diagrams.is_zero());
    for keep in [false, true] {
        let abort = Arc::new(AtomicBool::new(false));
        let calls = Arc::new(AtomicUsize::new(0));
        let flag = Arc::clone(&abort);
        let invoked = Arc::clone(&calls);
        let query = bare.clone().prune_dressings(move |view| {
            assert!(!view.is_complete());
            invoked.fetch_add(1, Ordering::Relaxed);
            flag.store(true, Ordering::Relaxed);
            keep
        });
        assert_predicate_abort(&Plan::new(&theory, &query).unwrap(), &abort, &calls);
    }
}

#[test]
fn abort_inside_complete_dressing_or_select_prevents_output() {
    // The three-leg tree has no internal edges; the one-loop bubble takes Single's fast path.
    let theory = phi(&[('p', 3)]);
    for n_legs in [2, 3] {
        let legs = (0..n_legs).fold(Legs::new(), |legs, label| {
            legs.leg(label, Stub::undirected('p'))
        });
        let loops = usize::from(n_legs == 2);
        let bare = Query::new(&legs).loops(loops..=loops);
        assert!(!Plan::new(&theory, &bare)
            .unwrap()
            .count(&Run::default())
            .unwrap()
            .diagrams
            .is_zero());
        for in_select in [false, true] {
            for keep in [false, true] {
                let abort = Arc::new(AtomicBool::new(false));
                let calls = Arc::new(AtomicUsize::new(0));
                let selections = Arc::new(AtomicUsize::new(0));
                let flag = Arc::clone(&abort);
                let invoked = Arc::clone(&calls);
                let query = bare.clone().prune_dressings(move |view| {
                    assert!(view.is_complete());
                    if in_select {
                        return true;
                    }
                    invoked.fetch_add(1, Ordering::Relaxed);
                    flag.store(true, Ordering::Relaxed);
                    keep
                });
                let flag = Arc::clone(&abort);
                let invoked = Arc::clone(&calls);
                let selected = Arc::clone(&selections);
                let query = query.select(move |_| {
                    selected.fetch_add(1, Ordering::Relaxed);
                    if in_select {
                        invoked.fetch_add(1, Ordering::Relaxed);
                        flag.store(true, Ordering::Relaxed);
                        keep
                    } else {
                        true
                    }
                });
                assert_predicate_abort(&Plan::new(&theory, &query).unwrap(), &abort, &calls);
                if !in_select {
                    assert_eq!(selections.load(Ordering::Relaxed), 0);
                }
            }
        }
    }
}

#[test]
fn retained_dressing_polls_abort_during_rejected_rule_search_without_hooks() {
    let p = Stub::undirected('p');
    let theory = (0..16)
        .fold(Theory::builder().weights(["order"]), |builder, rule| {
            builder.vertex(vec![p.clone(); 3], [2 * (rule % 2)])
        })
        .vertex(vec![Stub::undirected('q'); 4], [1])
        .build()
        .unwrap();
    // Planning admits one q vertex, but a connected graph with p legs uses only p rules.
    // Their state bounds admit one while every exact rule tuple has even total weight.
    let legs = Legs::new().leg(1u8, p.clone()).leg(2, p);
    let query = Query::new(&legs).loops(2..=2).weight("order", 1..=1);
    let plan = Plan::new(&theory, &query).unwrap();
    let run = Run::default();
    let snapshot = plan
        .skeletons(&run)
        .unwrap()
        .find(|skeleton| skeleton.view().n_vertices() == 4)
        .unwrap()
        .snapshot();
    assert!((0..4).all(|vertex| snapshot.view().degree(vertex) == 3));
    let polls = AtomicUsize::new(0);
    let run = Run::default().abort_check(|| polls.fetch_add(1, Ordering::Relaxed) >= 3);
    let mut dressings = plan.dressings(&snapshot, &run).unwrap();
    assert!(dressings.next().is_none());
    assert!(dressings.aborted());
}

#[test]
fn abort_from_rejected_last_dressing_is_reported_without_output() {
    let abort = Arc::new(AtomicBool::new(false));
    let flag = Arc::clone(&abort);
    let p = Stub::undirected('p');
    let legs = Legs::new().leg(1u8, p.clone()).leg(2, p.clone()).leg(3, p);
    let query = Query::new(&legs).loops(0..=0).select(move |_| {
        flag.store(true, Ordering::Relaxed);
        false
    });
    let plan = Plan::new(&phi(&[('p', 3)]), &query).unwrap();
    let run = Run::default().abort(&abort);
    let skeleton = plan.skeletons(&run).unwrap().next().unwrap();
    let mut dressings = skeleton.diagrams();
    assert!(dressings.next().is_none());
    assert!(dressings.aborted());
    abort.store(false, Ordering::Relaxed);
    let counts = plan.count(&run).unwrap();
    assert!(!counts.complete);
    assert!(counts.diagrams.is_zero());
}

#[test]
fn final_selection_cannot_hide_dressing_heredity_violations() {
    let theory = phi(&[('p', 3), ('p', 3)]);
    let query = Query::<char, u8>::new(&Legs::new())
        .loops(3..=3)
        .structure(Structure {
            self_loops: true,
            zero_flow_edges: true,
            ..Structure::default()
        })
        .prune_dressings(|view| view.is_complete())
        .select(|_| false);
    let plan = Plan::new(&theory, &query).unwrap();
    for parallel in [false, true] {
        let run = Run::default().check_heredity(true).parallel(parallel);
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| plan.count(&run))).is_err()
        );
    }
}

#[test]
fn one_theory_serves_separate_orders_and_owned_diagram_views() {
    let p = Stub::undirected('p');
    let theory = Theory::builder()
        .weights(["tag"])
        .vertex(vec![p.clone(); 3], [1])
        .vertex(vec![p.clone(); 3], [-1])
        .build()
        .unwrap();
    let legs = Legs::new().leg(1u8, p.clone()).leg(2, p);
    let query = Query::new(&legs).loops(1..=1);
    let positive = Plan::new(&theory, &query.clone().weight("tag", 2..=2)).unwrap();
    let negative = Plan::new(&theory, &query.weight("tag", -2..=-2)).unwrap();
    let run = Run::default();
    assert_eq!(positive.count(&run).unwrap().diagrams, BigUint::from(1u8));
    assert_eq!(negative.count(&run).unwrap().diagrams, BigUint::from(1u8));
    let owned = positive.diagrams(&run).unwrap().next().unwrap();
    assert_eq!(owned.rules(), &[0, 0]);
    let expected_id = owned.skeleton().skeleton_id().unwrap().clone();
    let expected_domains: Vec<_> = (0..2).map(|v| owned.skeleton().domains(v)).collect();
    assert_eq!(expected_domains, [vec![0, 1], vec![0, 1]]);
    drop((positive, negative, theory));
    assert_eq!(owned.skeleton().skeleton_id().unwrap(), &expected_id);
    assert_eq!(owned.skeleton().order_bounds().unwrap()[0].range, [-2, 2]);
    assert_eq!(owned.skeleton().domains(0), vec![0, 1]);
    assert_eq!(owned.graph().n_live(), 2);
}

#[test]
fn qed_explanation_and_shared_plan_keep_known_self_energy_count() {
    let (a, e) = (Stub::undirected('a'), Stub::incoming('e'));
    let theory = Theory::builder()
        .weights(["QED"])
        .vertex([e.clone(), e.flipped(), a.clone()], [1])
        .build()
        .unwrap();
    let legs = Legs::new().leg(1u8, a.clone()).leg(2, a);
    let seen = Arc::new(AtomicBool::new(false));
    let capture = Arc::clone(&seen);
    let query = Query::new(&legs)
        .loops(2..=2)
        .weight("QED", 4..=4)
        .structure(Structure {
            bridges: Bridges::OnePI,
            ..Structure::default()
        })
        .select(move |view| {
            assert_eq!(view.order_bounds().unwrap()[0].range, [4, 4]);
            assert!(view.skeleton_id().is_some());
            assert!((0..view.n_vertices()).all(|v| view.domains(v) == [0]));
            capture.store(true, Ordering::Relaxed);
            true
        });
    let plan = Plan::new(&theory, &query).unwrap();
    assert!(plan
        .explain()
        .criteria
        .iter()
        .any(|entry| entry.criterion == "weight QED"));
    let run = Run::default();
    assert_eq!(plan.count(&run).unwrap().diagrams, BigUint::from(3u8));
    assert_eq!(plan.diagrams(&run).unwrap().count(), 3);
    assert!(seen.load(Ordering::Relaxed));
}

#[test]
fn planning_reports_infeasible_bounds_and_preserves_low_degree_vertex_bound() {
    let theory = phi(&[('p', 2)]);
    let legs = Legs::<char, u8>::new();
    assert_eq!(
        Plan::new(&theory, &Query::new(&legs).loops(1..=1)).err(),
        Some(GenerationError::Unbounded)
    );
    let query = Query::new(&legs).loops(1..=1).structure(Structure {
        max_vertices: Some(3),
        self_loops: true,
        ..Structure::default()
    });
    assert!(Plan::new(&theory, &query).is_ok());
    assert!(matches!(
        Plan::new(&theory, &query.clone().loops(..0)),
        Err(GenerationError::Infeasible(_))
    ));
    let absent = Legs::new().leg(1u8, Stub::undirected('x'));
    assert!(matches!(
        Plan::new(&theory, &Query::new(&absent).loops(1..=1)),
        Err(GenerationError::Infeasible(_))
    ));
}
