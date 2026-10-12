"""Shared native/installed contract for Request hook registration."""

import gc
import json
import sys
import threading
import weakref
from pathlib import Path


def check_request_decorators(fk, model_json):
    unraisable = []
    previous = sys.unraisablehook

    def record_unraisable(event):
        unraisable.append((event.exc_type.__name__, str(event.exc_value)))

    try:
        sys.unraisablehook = record_unraisable
        _check_request_decorators(fk, model_json)
    finally:
        # The worker has exited; retained exception tracebacks must collect here.
        gc.collect()
        gc.collect()
        sys.unraisablehook = previous
    assert not unraisable, unraisable


def _check_request_decorators(fk, model_json):
    model = fk.Model.from_json(model_json)
    process = model.process([1000], [1000, 1000])
    criteria = {
        "loops": fk.Loops.total(1),
        "restrict": fk.Restrict().only_vertices(["V_3_SCALAR_000"]),
        "topology": fk.TopologyRules(max_vertices=3),
    }

    def request():
        return process.amplitude(**criteria)

    def effective(spec):
        text = repr(spec)
        assert text.startswith("Request(") and text.endswith(")")
        return json.loads(text[8:-1])

    caller = threading.get_ident()
    baseline = request().count().diagrams
    assert baseline > 0
    spec = request()
    initial = effective(spec)
    assert not any(
        initial[name] for name in ("prune_skeletons", "prune_dressings", "select")
    )
    assert not initial["check_heredity"]
    threads = {name: set() for name in ("SkeletonView", "DressingView", "DiagramView")}
    retained = []

    def record(view):
        kind = type(view).__name__
        threads[kind].add(threading.get_ident())
        retained.append(view)
        assert view.n_vertices > 0
        assert all(
            stub is None or stub.data.pdg_code == 1000 for _, _, stub in view.edges()
        )
        return True

    @spec.prune_skeletons
    def topology(view):
        assert type(view).__name__ == "SkeletonView"
        return record(view)

    @spec.prune_dressings(thread_safe=True)
    def dressing(view):
        assert type(view).__name__ == "DressingView"
        return record(view)

    @spec.select()
    def complete(view):
        assert type(view).__name__ == "DiagramView" and view.is_complete
        return record(view)

    assert spec.prune_skeletons(topology) is topology
    assert spec.prune_dressings(thread_safe=True)(dressing) is dressing
    assert spec.select(complete) is complete
    state = effective(spec)
    assert all(state[name] for name in ("prune_skeletons", "prune_dressings", "select"))
    assert state["thread_safe_hooks"] == {
        "prune_skeletons": False,
        "prune_dressings": True,
        "select": False,
    }
    progress_threads = set()
    run = fk.Run(
        threads=2, progress=lambda p: progress_threads.add(threading.get_ident())
    )
    assert "progress=callable" in repr(run) and "thread_safe=false" in repr(run)
    assert spec.count(run=run).diagrams == baseline
    assert progress_threads == {caller}
    assert threads["SkeletonView"] == threads["DiagramView"] == {caller}
    assert len(threads["DressingView"]) == 1
    assert (
        threads["DressingView"] == {caller}
        if sys.platform in ("emscripten", "wasi")
        else caller not in threads["DressingView"]
    )
    for view in retained:
        try:
            _ = view.n_vertices
        except ReferenceError:
            pass
        else:
            raise AssertionError("decorated callback view escaped its scope")
    retained.clear()
    del view
    # Replacement updates both callback and permission; invalid replacement is atomic.
    for hook in ("prune_skeletons", "prune_dressings", "select"):
        current = request()
        previous = lambda view: False
        getattr(current, hook)(previous, thread_safe=True)
        replacement_threads = set()

        def replacement(view, seen=replacement_threads):
            seen.add(threading.get_ident())
            return True

        assert getattr(current, hook)(replacement) is replacement
        before = effective(current)
        try:
            getattr(current, hook)(object())
        except TypeError as error:
            assert "callable" in str(error)
        else:
            raise AssertionError("non-callable decorator accepted")
        assert effective(current) == before
        assert current.count().diagrams == baseline
        assert replacement_threads == {caller}
        assert not effective(current)["thread_safe_hooks"][hook]

    # The same exception/cancellation boundary applies to every registration and thread route.
    for hook in ("prune_skeletons", "prune_dressings", "select"):
        for worker in (False, True):
            current = request()
            expected = RuntimeError("decorated callback failed")
            calls = []

            def fail(view, calls=calls, expected=expected):
                calls.append(1)
                raise expected

            getattr(current, hook)(thread_safe=worker)(fail)
            try:
                current.generate(run=fk.Run(threads=2))
            except RuntimeError as error:
                assert error is expected
            else:
                raise AssertionError("decorator exception was lost")
            assert calls == [1]
            current = request()
            calls.clear()
            token = fk.CancellationToken()

            def cancel(view, calls=calls, token=token):
                calls.append(1)
                token.cancel()
                return True

            getattr(current, hook)(thread_safe=worker)(cancel)
            cancelled = current.count(run=fk.Run(threads=2, cancellation=token))
            assert calls == [1] and not cancelled.complete

    # Re-registering from an active callback reports a Python borrow error, never a panic.
    current = request()

    @current.prune_skeletons
    def reenter(view):
        current.prune_dressings(lambda inner: True)
        return True

    try:
        current.count()
    except RuntimeError as error:
        assert "borrow" in str(error).lower()
    else:
        raise AssertionError("active request registration did not fail")
    current.prune_skeletons(lambda view: True)
    assert current.count().diagrams == baseline

    # functools.partial owns parameterized registration and is collected with existing Request GC.
    def cycle():
        current = request()
        decorator = current.prune_dressings(thread_safe=True)

        def keep(view, held=(current, decorator)):
            return True

        decorator(keep)
        return weakref.ref(keep), weakref.ref(decorator)

    references = cycle()
    gc.collect()
    assert all(reference() is None for reference in references)

    # Request and Run preferences monotonically enable the same native heredity check.
    for request_check, run_check in (
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ):
        # Two possible internal fields create actual partial assignments. A singleton
        # rule has only complete decisions and cannot witness non-hereditary pruning.
        current = process.amplitude(
            **(
                criteria
                | {
                    "restrict": fk.Restrict().only_vertices(
                        ["V_3_SCALAR_000", "V_3_SCALAR_011"]
                    )
                }
            )
        )
        decisions = []
        assert current.check_heredity is False
        current.check_heredity = request_check
        assert current.check_heredity is request_check
        assert effective(current)["check_heredity"] is request_check

        @current.prune_dressings
        def nonhereditary(view, decisions=decisions):
            decisions.append(view.is_complete)
            return view.is_complete

        @current.select
        def reject_complete(view):
            return False

        if not (request_check or run_check):
            counts = current.count(run=fk.Run(check_heredity=False))
            assert counts.complete and counts.diagrams == 0
            assert decisions and not any(decisions)
            continue
        try:
            current.count(run=fk.Run(check_heredity=run_check))
        except fk.GenerationError as error:
            assert "hereditary" in str(error)
            assert False in decisions and True in decisions
        else:
            raise AssertionError("request/Run heredity preference was ignored")


if __name__ == "__main__":
    from symbolica.community import hepkit as fk

    check_request_decorators(
        fk, (Path(__file__).parent / "fixtures/scalars_2p_3p.json").read_text()
    )
    print(
        "Request decorators: dispatch, replacement, GC, exceptions, cancellation and heredity passed"
    )
