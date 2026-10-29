import threading
import time

fk = globals()["fk"]  # Supplied by the native test or standalone probe.

request = (
    fk.Model.qcd().process(["u", "u~"], ["g", "g"]).amplitude(loops=fk.Loops.total(2))
)
caller = threading.get_ident()
seen = []
stopped = False


def progress(update):
    global stopped
    assert threading.get_ident() == caller
    seen.append((update.stage, update.completed, stopped))
    if not stopped:
        assert update.stage == fk.Stage.Skeletons
        # Let the native worker advance while the caller handles this update.
        # The Rust test separately controls the cross-batch boundary without time.
        time.sleep(0.25)
        stopped = True
        return False
    return None


result = request.generate(run=fk.Run(threads=1, progress=progress))
print(
    [(str(stage), count, after_false) for stage, count, after_false in seen], flush=True
)
assert not result.report.completed
assert seen[0][0] == fk.Stage.Skeletons, seen
assert seen[-1] == (fk.Stage.Done, result.report.retained_count, True), seen
assert seen[-1][1] == len(result), seen
assert sum(stage == fk.Stage.Done for stage, _, _ in seen) == 1, seen
assert all(stage == fk.Stage.Done for stage, _, after_false in seen if after_false), (
    seen
)
assert len(seen) == 2, seen
delivered = len(seen)
time.sleep(0.02)
assert len(seen) == delivered
