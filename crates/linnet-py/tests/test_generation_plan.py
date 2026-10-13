"""Reusable plans preserve exact counts, rule identities and Python object lifetimes."""

import gc
import re
import sys
import threading
import unittest
import weakref
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

from linnet_under_test import field
from linnet_under_test import graph as g


class _Owner:
    pass


class _Field:
    def __init__(self, owner):
        self.owner = owner

    def __repr__(self):
        return "GenerationField('phi')"


def _two_flavour_qed(keys=None):
    # Class members straddle the fixed photon in key order; each key is itself a tuple.
    if keys is None:
        keys = (("a", ("lepton", 1)), ("b", ("photon",)), ("c", ("lepton", 2)))
    first, photon_key, second = keys
    photon = g.EdgeSignature(photon_key)
    charged = [g.EdgeSignature(key, direction=False) for key in (first, second)]
    theory = g.Theory([[line, line.flip(), photon] for line in charged])
    return theory, g.Legs([(1, photon), (2, photon)]), (first, second)


class GenerationPlanTests(unittest.TestCase):
    def test_named_signed_weights_keep_duplicate_rules_distinct(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3, [phi] * 3], weights={"tag": [1, -1]})
        legs = g.Legs([(1, phi), (2, phi)])
        self.assertEqual(theory.n_rules, 2)
        self.assertEqual(theory.weight_names, ["tag"])

        forms = []
        for weight, rule in [(2, 0), (-2, 1)]:
            with self.subTest(weight=weight):
                plan = theory.query(legs, loops=1, weights={"tag": (weight, weight)})
                self.assertEqual(plan.count().diagrams, 1)
                (diagram,) = plan.diagrams()
                self.assertEqual(diagram.rules(), [rule, rule])
                forms.append(diagram.canonical_form())
        self.assertNotEqual(*forms)

    def test_repeated_counts_and_streams_match_exact_vacuum_reference(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3, [phi] * 3])
        plan = theory.query(
            g.Legs(),
            loops=3,
            structure=g.Structure(self_loops=True, zero_flow_edges=True),
        )
        reference = None
        for parallel in (False, True):
            with self.subTest(parallel=parallel):
                run = g.Run(parallel=parallel)
                first = plan.count(run)
                stream = plan.diagrams(run)
                self.assertIs(iter(stream), stream)
                diagrams = list(stream)
                self.assertEqual(list(stream), [])
                self.assertFalse(stream.aborted)
                second = plan.count(run)

                # Independent native scalar oracle: two cubic rules at L=3 sum to 5.
                for counts in (first, second):
                    self.assertIs(type(counts), g.Counts)
                    self.assertIs(type(counts.diagrams), int)
                    self.assertIs(type(counts.skeletons), int)
                    self.assertIs(type(counts.sum_inv_aut), Fraction)
                    self.assertTrue(counts.complete)
                    self.assertGreater(counts.skeletons, 0)
                    self.assertEqual(counts.diagrams, len(diagrams))
                    self.assertEqual(counts.sum_inv_aut, Fraction(5))
                    self.assertEqual(set(counts.per_loop), {3})
                    row = counts.per_loop[3]
                    self.assertIs(type(row), g.LoopCounts)
                    self.assertEqual(row.diagrams, counts.diagrams)
                    self.assertEqual(row.skeletons, counts.skeletons)
                    self.assertEqual(row.sum_inv_aut, counts.sum_inv_aut)
                self.assertEqual(first.skeletons, second.skeletons)
                self.assertTrue(all(diagram.loops == 3 for diagram in diagrams))
                self.assertEqual(
                    sum(
                        (Fraction(1, diagram.symmetry_factor()) for diagram in diagrams),
                        Fraction(0),
                    ),
                    Fraction(5),
                )
                forms = [diagram.canonical_form() for diagram in diagrams]
                self.assertEqual(len(set(forms)), len(forms))
                skeletons = plan.skeletons(run)
                self.assertIs(iter(skeletons), skeletons)
                dressed = [
                    diagram.canonical_form()
                    for skeleton in skeletons
                    for diagram in skeleton.diagrams()
                ]
                self.assertCountEqual(dressed, forms)
                self.assertEqual(list(skeletons), [])
                self.assertFalse(skeletons.aborted)
                if reference is None:
                    reference = forms
                else:
                    self.assertEqual(forms, reference)

    def test_query_snapshots_legs_before_later_mutation(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3])
        legs = g.Legs([(1, phi), (2, phi), (3, phi)])
        query = g.Query(legs, loops=0)
        self.assertIs(legs.leg(4, phi), legs)
        self.assertEqual(len(legs), 4)

        original = g.Plan(theory, query)
        (diagram,) = original.diagrams()
        self.assertEqual(original.count().diagrams, 1)
        self.assertEqual(sorted(label for _, label, _ in diagram.legs()), [1, 2, 3])
        self.assertEqual(theory.query(legs, loops=0).count().diagrams, 3)

    def test_vertex_counts_select_the_declared_rule_index(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3, [phi] * 3])
        plan = theory.query(
            g.Legs([(1, phi), (2, phi)]),
            loops=1,
            counts=[(g.Countable.vertices([1]), (2, 2))],
        )
        self.assertEqual(plan.count().diagrams, 1)
        (diagram,) = plan.diagrams()
        self.assertEqual(diagram.rules(), [1, 1])

    def test_qed_plan_explains_the_known_three_diagram_count(self):
        photon = g.EdgeSignature(field("a"))
        electron = g.EdgeSignature(field("e"), direction=False)
        theory = g.Theory(
            [[electron, electron.flip(), photon]], weights={"QED": [1]}
        )
        plan = theory.query(
            g.Legs([(1, photon), (2, photon)]),
            loops=2,
            weights={"QED": (4, 4)},
            structure=g.Structure.one_pi(),
        )
        report = plan.explain()
        self.assertIs(type(report), g.PlanReport)
        self.assertIsNone(report.infeasible)
        self.assertTrue(all(type(item) is g.PlanCriterion for item in report.criteria))
        self.assertIn("weight QED", [item.criterion for item in report.criteria])
        self.assertIn("weight QED", str(report))
        self.assertEqual(plan.count().diagrams, 3)
        self.assertEqual(len(list(plan.diagrams())), 3)

    def test_finite_vertex_bounds_and_infeasible_reports(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 2])
        legs = g.Legs()
        with self.assertRaisesRegex(ValueError, "bound"):
            theory.query(legs)
        with self.assertRaisesRegex(ValueError, "vertex bound"):
            theory.query(legs, loops=1)

        structure = g.Structure(max_vertices=3, self_loops=True)
        # Generic Linnet may omit loops when a vertex bound makes the search finite.
        plan = g.Plan(theory, g.Query(legs, structure=structure))
        counts = plan.count()
        self.assertTrue(counts.complete)
        self.assertGreater(counts.diagrams, 0)

        with self.assertRaises(ValueError) as raised:
            theory.query(legs, loops=range(0), structure=structure)
        report = raised.exception.report
        self.assertIs(type(report), g.PlanReport)
        self.assertIsInstance(report.infeasible, str)
        self.assertTrue(report.infeasible)
        self.assertIn(report.infeasible, str(report))
        self.assertIsInstance(report.empty_partitions, list)
        loops = next(item for item in report.criteria if item.criterion == "loops")
        self.assertEqual(loops.stage, "planning")
        self.assertEqual(loops.role, "restriction")
        self.assertTrue(loops.hereditary)

    def test_compact_diagram_survives_its_plan_and_theory(self):
        def retain_diagram():
            phi = g.EdgeSignature(field("phi"))
            theory = g.Theory([[phi] * 3])
            plan = theory.query(g.Legs([(1, phi), (2, phi)]), loops=1)
            return next(plan.diagrams())

        diagram = retain_diagram()
        gc.collect()
        self.assertIs(type(diagram), g.Diagram)
        self.assertEqual(diagram.loops, 1)
        self.assertEqual(diagram.rules(), [0, 0])
        self.assertEqual(len(diagram.edges()), 2)
        self.assertEqual(sorted(label for _, label, _ in diagram.legs()), [1, 2])
        identity = diagram.canonical_form()
        hash(identity)
        graph = diagram.graph()
        self.assertEqual(graph.n_nodes, 2)
        self.assertEqual(graph.n_edges, 4)
        self.assertEqual(graph.canonize()[2], diagram.symmetry_factor())
        graph.add_node(g.Node("later"))
        self.assertEqual(graph.n_nodes, 3)
        self.assertEqual(diagram.graph().n_nodes, 2)
        self.assertEqual(diagram.canonical_form(), identity)

    def test_retained_skeleton_repeats_dressing_after_its_source_is_deleted(self):
        def retain_skeleton():
            phi = g.EdgeSignature(field("phi"))
            theory = g.Theory([[phi] * 3, [phi] * 3])
            plan = theory.query(g.Legs([(1, phi), (2, phi)]), loops=1)
            stream = plan.skeletons()
            return next(stream)

        skeleton = retain_skeleton()
        gc.collect()
        self.assertIs(type(skeleton), g.Skeleton)
        self.assertEqual(skeleton.loops, 1)
        self.assertGreater(skeleton.automorphisms(), 0)
        identity = skeleton.id()
        hash(identity)
        self.assertEqual(skeleton.graph().n_nodes, 2)
        first = [diagram.canonical_form() for diagram in skeleton.diagrams()]
        self.assertTrue(first)
        second = [diagram.canonical_form() for diagram in skeleton.diagrams()]
        self.assertEqual(first, second)
        self.assertEqual(skeleton.id(), identity)

    def test_custom_field_cycles_through_theory_and_diagram_are_collected(self):
        def cycle(output):
            owner = _Owner()
            phi = g.EdgeSignature(_Field(owner))
            theory = g.Theory([[phi] * 3])
            if output == "theory":
                owner.output = theory
            else:
                plan = theory.query(
                    g.Legs([(1, phi), (2, phi), (3, phi)]), loops=0
                )
                owner.output = next(plan.diagrams())
            # owner -> output -> field tuple -> custom field -> owner.
            return weakref.ref(owner)

        for output in ("theory", "diagram"):
            with self.subTest(output=output):
                probe = cycle(output)
                gc.collect()
                self.assertIsNone(probe())

    @unittest.skipIf(sys.platform in {"emscripten", "wasi"}, "requires native Python threads")
    def test_parallel_callbacks_use_the_invoking_thread_and_keep_exact_counts(self):
        self._check_callbacks_use_invoking_thread(worker=True)

    def test_inline_callbacks_use_the_invoking_thread_and_keep_exact_counts(self):
        self._check_callbacks_use_invoking_thread(worker=False)

    def _check_callbacks_use_invoking_thread(self, worker):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3, [phi] * 3])
        structure = g.Structure(self_loops=True, zero_flow_edges=True)
        expected = theory.query(g.Legs(), loops=3, structure=structure).count()
        calls = []

        def keep(value):
            calls.append((type(value), threading.get_ident()))
            return True

        plan = theory.query(
            g.Legs(),
            loops=3,
            structure=structure,
            prune_skeletons=keep,
            prune_dressings=keep,
            select=keep,
        )
        run = g.Run(parallel=True, progress=keep)

        def count():
            return threading.get_ident(), plan.count(run)

        if worker:
            # Construct the plan here, but advance it from another Python thread.
            with ThreadPoolExecutor(max_workers=1) as executor:
                invoking_thread, actual = executor.submit(count).result()
            self.assertNotEqual(invoking_thread, threading.get_ident())
        else:
            invoking_thread, actual = count()
            self.assertEqual(invoking_thread, threading.get_ident())
        self.assertEqual({thread for _, thread in calls}, {invoking_thread})
        self.assertEqual(
            {kind for kind, _ in calls},
            {g.SkeletonView, g.DressingView, g.DiagramView, g.Progress},
        )
        self.assertTrue(actual.complete)
        self.assertEqual(actual.skeletons, expected.skeletons)
        self.assertEqual(actual.diagrams, expected.diagrams)
        self.assertEqual(actual.sum_inv_aut, expected.sum_inv_aut)
        self.assertEqual(actual.sum_inv_aut, Fraction(5))
        self.assertEqual(actual.per_loop[3].sum_inv_aut, expected.per_loop[3].sum_inv_aut)
        self.assertGreater(run.callback_round_trips, 0)
        self.assertLessEqual(
            run.callback_round_trips, sum(kind is not g.Progress for kind, _ in calls)
        )

    def test_final_selection_sees_typed_complete_diagrams(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3, [phi] * 3])
        legs = g.Legs([(1, phi), (2, phi)])
        baseline = list(theory.query(legs, loops=1).diagrams())
        expected = [
            diagram.canonical_form()
            for diagram in baseline
            if diagram.rules().count(1) == 1
        ]
        self.assertGreater(len(expected), 0)
        self.assertLess(len(expected), len(baseline))

        def exactly_one_second_rule(view):
            self.assertIs(type(view), g.DiagramView)
            self.assertEqual(view.stage, g.GenerationStage.Complete)
            self.assertTrue(view.is_complete)
            return sum(view.rules(v) == [1] for v in range(view.n_vertices)) == 1

        plan = theory.query(legs, loops=1, select=exactly_one_second_rule)
        actual = [diagram.canonical_form() for diagram in plan.diagrams(parallel=True)]
        self.assertCountEqual(actual, expected)
        self.assertEqual(plan.count().diagrams, len(expected))

    def test_final_rejection_does_not_hide_nonhereditary_dressing(self):
        phi = g.EdgeSignature(field("phi"))
        plan = g.Theory([[phi] * 3, [phi] * 3]).query(
            g.Legs(),
            loops=3,
            structure=g.Structure(self_loops=True, zero_flow_edges=True),
            prune_dressings=lambda view: view.is_complete,
            select=lambda view: False,
        )
        for parallel in (False, True):
            with self.subTest(parallel=parallel):
                with self.assertRaisesRegex(RuntimeError, "hereditary"):
                    plan.count(g.Run(parallel=parallel, check_heredity=True))

    def test_saved_callback_views_expire_on_every_exit_path(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3])
        legs = g.Legs([(i, phi) for i in range(4)])
        for hook, view_type in (
            ("prune_skeletons", g.SkeletonView),
            ("prune_dressings", g.DressingView),
            ("select", g.DiagramView),
        ):
            for outcome in ("keep", "reject", "type_error", "exception"):
                with self.subTest(hook=hook, outcome=outcome):
                    saved = []

                    def predicate(view):
                        self.assertIs(type(view), view_type)
                        saved.append(view)
                        if outcome == "exception":
                            raise RuntimeError("saved-view sentinel")
                        if outcome == "type_error":
                            return object()
                        return outcome == "keep"

                    stream = theory.query(legs, loops=0, **{hook: predicate}).diagrams(
                        parallel=True
                    )
                    if outcome in ("keep", "reject"):
                        self.assertEqual(len(list(stream)), 3 if outcome == "keep" else 0)
                    else:
                        error = TypeError if outcome == "type_error" else RuntimeError
                        message = (
                            "bool" if outcome == "type_error" else "saved-view sentinel"
                        )
                        with self.assertRaisesRegex(error, message):
                            next(stream)
                        self.assertEqual(len(saved), 1)
                        self.assertTrue(stream.aborted)
                        self.assertEqual(list(stream), [])
                    self.assertTrue(saved)
                    for view in saved:
                        with self.assertRaises(ReferenceError):
                            _ = view.stage
                        with self.assertRaises(ReferenceError):
                            view.graph()
                        with self.assertRaises(ReferenceError):
                            view.field(0)
                        with self.assertRaises(ReferenceError):
                            view.flow(0)

    def test_expired_callback_views_release_their_python_keys(self):
        saved = []

        def keep(view):
            saved.append(view)
            return True

        def execute():
            owner = _Owner()
            phi = g.EdgeSignature(_Field(owner))
            plan = g.Theory([[phi] * 3]).query(
                g.Legs([(i, phi) for i in range(3)]), loops=0, select=keep
            )
            self.assertEqual(plan.count().diagrams, 1)
            return weakref.ref(owner)

        probe = execute()
        gc.collect()
        self.assertEqual(len(saved), 1)
        self.assertIsNone(probe())
        with self.assertRaises(ReferenceError):
            saved[0].graph()

    @unittest.skipIf(sys.platform in {"emscripten", "wasi"}, "requires native Python threads")
    def test_owned_and_expired_views_can_cross_python_threads(self):
        phi = g.EdgeSignature(field("phi"))
        saved = []

        def keep(view):
            saved.append(view)
            return True

        plan = g.Theory([[phi] * 3]).query(
            g.Legs([(i, phi) for i in range(4)]), loops=0, select=keep
        )
        diagrams = list(plan.diagrams())
        owned = diagrams[0].skeleton()

        def inspect_views():
            for view in saved:
                with self.assertRaises(ReferenceError):
                    _ = view.stage
            return owned.stage, owned.n_vertices, len(list(owned.graph().nodes()))

        with ThreadPoolExecutor(max_workers=1) as executor:
            observed = executor.submit(inspect_views).result()
        self.assertEqual(observed, (g.GenerationStage.Topology, 2, 2))

    def test_type_conversion_failure_prevents_the_next_predicate(self):
        phi = g.EdgeSignature(field("phi"))
        calls = []

        def invalid(view):
            calls.append("skeleton")
            return 1

        def later(view):
            calls.append("later")
            return True

        plan = g.Theory([[phi] * 3]).query(
            g.Legs([(i, phi) for i in range(4)]),
            loops=0,
            prune_skeletons=invalid,
            prune_dressings=later,
            select=later,
        )
        stream = plan.diagrams(parallel=True)
        with self.assertRaisesRegex(TypeError, "bool"):
            next(stream)
        self.assertEqual(calls, ["skeleton"])
        self.assertEqual(list(stream), [])
        self.assertEqual(calls, ["skeleton"])

    def test_cancellation_inside_a_predicate_stops_later_calls(self):
        phi = g.EdgeSignature(field("phi"))
        run = g.Run(parallel=True)
        saved = []

        def cancel(view):
            saved.append(view)
            run.cancel()
            return True

        plan = g.Theory([[phi] * 3]).query(
            g.Legs([(i, phi) for i in range(4)]),
            loops=0,
            prune_skeletons=cancel,
            prune_dressings=cancel,
            select=cancel,
        )
        stream = plan.diagrams(run)
        self.assertEqual(list(stream), [])
        self.assertTrue(run.cancelled)
        self.assertTrue(stream.aborted)
        self.assertEqual(len(saved), 1)
        with self.assertRaises(ReferenceError):
            saved[0].degree(0)

    def test_failed_and_cancelled_runs_do_not_poison_another_execution(self):
        phi = g.EdgeSignature(field("phi"))
        theory = g.Theory([[phi] * 3])
        legs = g.Legs([(i, phi) for i in range(4)])
        expected = [d.canonical_form() for d in theory.query(legs, loops=0).diagrams()]
        cancel_run = g.Run(parallel=True)
        mode = "error"

        def keep(view):
            if mode == "error":
                raise RuntimeError("independent-run sentinel")
            if mode == "cancel":
                cancel_run.cancel()
            return True

        plan = theory.query(legs, loops=0, prune_skeletons=keep)
        healthy_run = g.Run(parallel=True)
        healthy = plan.diagrams(healthy_run)
        failed = plan.diagrams(g.Run(parallel=True))
        cancelled = plan.diagrams(cancel_run)
        with self.assertRaisesRegex(RuntimeError, "independent-run sentinel"):
            next(failed)
        self.assertTrue(failed.aborted)
        self.assertEqual(list(failed), [])
        mode = "cancel"
        self.assertEqual(list(cancelled), [])
        self.assertTrue(cancelled.aborted)
        mode = "keep"
        self.assertFalse(healthy_run.cancelled)
        self.assertCountEqual([d.canonical_form() for d in healthy], expected)
        self.assertFalse(healthy.aborted)
        self.assertEqual(plan.count(g.Run()).diagrams, 3)

    def test_progress_veto_delivers_accepted_item_but_cancel_discards_it(self):
        phi = g.EdgeSignature(field("phi"))
        plan = g.Theory([[phi] * 3]).query(
            g.Legs([(i, phi) for i in range(4)]), loops=0
        )
        for cancel in (False, True):
            with self.subTest(cancel=cancel):
                reports = []

                def progress(report):
                    self.assertIs(type(report), g.Progress)
                    reports.append(report.diagrams)
                    if report.diagrams == 0:
                        return True
                    if cancel:
                        run.cancel()
                        return True
                    return False

                run = g.Run(parallel=True, progress=progress)
                stream = plan.diagrams(run)
                self.assertEqual(len(list(stream)), 0 if cancel else 1)
                self.assertTrue(stream.aborted)
                self.assertEqual(reports[-1], 1)
                self.assertEqual(run.cancelled, cancel)

    def test_flavour_class_modes_preserve_canonical_maps_and_original_rules(self):
        theory, legs, _ = _two_flavour_qed()
        reference = None
        for classes in (g.Classes.none(), g.Classes.auto(), None):
            for direct in (False, True):
                with self.subTest(classes=classes, direct=direct):
                    options = {} if classes is None else {"classes": classes}
                    plan = (
                        theory.query(legs, loops=1, **options)
                        if direct
                        else g.Plan(theory, g.Query(legs, loops=1, **options))
                    )
                    diagrams = list(plan.diagrams())
                    actual = {
                        diagram.canonical_form(): (
                            diagram.symmetry_factor(),
                            tuple(sorted(diagram.rules())),
                        )
                        for diagram in diagrams
                    }
                    self.assertEqual(len(actual), len(diagrams))
                    self.assertEqual(len(actual), 2)
                    self.assertEqual(set(actual.values()), {(1, (0, 0)), (1, (1, 1))})
                    counts = plan.count()
                    self.assertTrue(counts.complete)
                    self.assertEqual(counts.diagrams, 2)
                    self.assertEqual(counts.sum_inv_aut, Fraction(2))
                    dressed = {
                        diagram.canonical_form(): diagram.symmetry_factor()
                        for skeleton in plan.skeletons()
                        for diagram in skeleton.diagrams()
                    }
                    self.assertEqual(
                        dressed, {key: value[0] for key, value in actual.items()}
                    )
                    if reference is None:
                        reference = actual
                    else:
                        self.assertEqual(actual, reference)

    def test_class_fields_keep_original_tuple_keys_and_directed_flow(self):
        def execute():
            theory, legs, keys = _two_flavour_qed()
            saved_views, owned = [], []
            seen = set()

            def inspect(view):
                saved_views.append(view)
                for index, (_, _, stub) in enumerate(view.edges()):
                    value = view.field(index)
                    flow = view.flow(index)
                    if value is None:
                        self.assertIsNone(stub)
                        self.assertIsNone(flow)
                        continue
                    self.assertIs(type(value), g.Dressed)
                    seen.add(value.kind)
                    owned.append(value)
                    if value.kind == "class":
                        self.assertIsNone(stub)
                        self.assertIsNone(value.particle)
                        self.assertIs(type(value.members), tuple)
                        self.assertEqual(len(value.members), 2)
                        self.assertIs(value.members[0], keys[0])
                        self.assertIs(value.members[1], keys[1])
                        self.assertIn(flow, (g.Flow.Source, g.Flow.Sink))
                    else:
                        self.assertEqual(value.kind, "particle")
                        self.assertIsNone(value.members)
                        self.assertTrue(any(value.particle is key for key in keys))
                        self.assertIs(stub.data, value.particle)
                        self.assertEqual(
                            flow, g.Flow.Source if stub.direction else g.Flow.Sink
                        )
                    if type(view) is g.DiagramView:
                        self.assertEqual(value.kind, "particle")
                return True

            # Omitting classes must take the Auto route and expose a real class state.
            plan = theory.query(legs, loops=1, prune_dressings=inspect, select=inspect)
            self.assertEqual(plan.count(g.Run(parallel=True)).diagrams, 2)
            self.assertEqual(seen, {"class", "particle"})
            return saved_views, owned, keys

        views, owned, keys = execute()
        gc.collect()
        for view in views:
            with self.assertRaises(ReferenceError):
                view.field(0)
            with self.assertRaises(ReferenceError):
                view.flow(0)
        for value in owned:
            if value.kind == "class":
                self.assertIs(value.members[0], keys[0])
                self.assertIs(value.members[1], keys[1])
                self.assertIsNone(value.particle)
            else:
                # A concrete tuple-valued key is not a tuple of class members.
                self.assertIs(type(value.particle), tuple)
                self.assertTrue(any(value.particle is key for key in keys))
                self.assertIsNone(value.members)

    def test_owned_dressed_values_participate_in_key_cycles(self):
        def cycle(kind):
            owner = _Owner()
            payload = _Field(owner)
            theory, legs, _ = _two_flavour_qed(
                (("a", payload), ("b", "photon"), ("c", payload))
            )

            def capture(view):
                for index in range(len(view.edges())):
                    value = view.field(index)
                    if value is not None and value.kind == kind:
                        owner.output = value
                return True

            self.assertEqual(
                theory.query(legs, loops=1, prune_dressings=capture).count().diagrams, 2
            )
            self.assertIs(type(owner.output), g.Dressed)
            # owner -> Dressed -> original tuple key -> custom payload -> owner.
            return weakref.ref(owner)

        for kind in ("class", "particle"):
            with self.subTest(kind=kind):
                probe = cycle(kind)
                gc.collect()
                self.assertIsNone(probe())

    def test_class_callback_failure_and_cancel_expire_views_without_poisoning_plan(
        self,
    ):
        theory, legs, keys = _two_flavour_qed()
        for outcome in ("type_error", "exception", "cancel"):
            with self.subTest(outcome=outcome):
                run = g.Run(parallel=True)
                saved = []
                stopped = False
                mode = outcome
                sentinel = RuntimeError("class callback sentinel")

                def inspect(view):
                    nonlocal stopped
                    self.assertFalse(
                        stopped, "a callback ran after failure or cancellation"
                    )
                    for index in range(len(view.edges())):
                        value = view.field(index)
                        if value is None or value.kind != "class" or mode == "keep":
                            continue
                        saved.append((view, value))
                        stopped = True
                        if mode == "exception":
                            raise sentinel
                        if mode == "type_error":
                            return 1
                        run.cancel()
                        return True
                    return True

                plan = theory.query(legs, loops=1, prune_dressings=inspect)
                stream = plan.diagrams(run)
                if outcome == "cancel":
                    self.assertEqual(list(stream), [])
                    self.assertTrue(run.cancelled)
                else:
                    error = TypeError if outcome == "type_error" else RuntimeError
                    with self.assertRaises(error) as caught:
                        list(stream)
                    if outcome == "exception":
                        self.assertIs(caught.exception, sentinel)
                self.assertTrue(stream.aborted)
                self.assertEqual(len(saved), 1)
                view, value = saved[0]
                self.assertEqual(value.kind, "class")
                self.assertIs(value.members[0], keys[0])
                self.assertIs(value.members[1], keys[1])
                with self.assertRaises(ReferenceError):
                    view.field(0)
                with self.assertRaises(ReferenceError):
                    view.flow(0)
                mode, stopped = "keep", False
                self.assertEqual(plan.count(g.Run()).diagrams, 2)

    def test_rejecting_a_class_but_accepting_its_particles_is_not_hereditary(self):
        theory, legs, _ = _two_flavour_qed()

        def only_particles(view):
            return all(
                view.field(index) is None or view.field(index).kind != "class"
                for index in range(len(view.edges()))
            )

        plan = theory.query(legs, loops=1, prune_dressings=only_particles)
        for parallel in (False, True):
            with self.subTest(parallel=parallel):
                with self.assertRaisesRegex(RuntimeError, "hereditary"):
                    plan.count(g.Run(parallel=parallel, check_heredity=True))

    def test_published_generation_example_runs(self):
        source_path = (
            Path(__file__).resolve().parents[3] / "docs/products/linnet/content/api.typ"
        )
        examples = [
            code
            for code in re.findall(
                r"```python\n(.*?)\n```",
                source_path.read_text(encoding="utf-8"),
                re.DOTALL,
            )
            if "theory = graph.Theory(vertices=" in code
        ]
        self.assertEqual(len(examples), 1)
        code = examples[0]
        module_import = "import linnet_py as graph\n"
        self.assertEqual(code.count(module_import), 1)
        code = code.replace(module_import, "\n", 1)
        exec(compile(code, str(source_path), "exec"), {"graph": g})


if __name__ == "__main__":
    unittest.main()
