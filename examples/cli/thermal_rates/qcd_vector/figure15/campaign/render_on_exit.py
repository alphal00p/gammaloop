"""Wait on kernel process descriptors, then render completed scan artifacts."""
import argparse
import os
from pathlib import Path
import select
import subprocess
import sys

from common import process_identity, read_json, same_process, utc_now, write_json

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("run_root", type=Path)
args = parser.parse_args()
root = args.run_root.resolve()
controller = read_json(root / "campaign.json")["controller_identity"]
watchdog = read_json(root / "watchdog.json")["identity"]
poller = select.poll()
descriptors = set()
for identity in (controller, watchdog):
    if not same_process(identity):
        raise SystemExit("Cannot arm renderer: a registered process already exited")
    descriptor = os.pidfd_open(identity["pid"])
    descriptors.add(descriptor)
    poller.register(descriptor, select.POLLIN)
status = dict(status="armed", armed_utc=utc_now(), identity=process_identity(os.getpid()),
              controller=controller, watchdog=watchdog,
              notification="Local report and plots only; no automatic chat notification")
write_json(root / "completion-renderer.json", status)
while descriptors:
    for descriptor, _ in poller.poll():
        poller.unregister(descriptor)
        descriptors.remove(descriptor)
        os.close(descriptor)
status.update(status="rendering", compute_exited_utc=utc_now())
write_json(root / "completion-renderer.json", status)
with (root / "completion-renderer.log").open("w") as log:
    result = subprocess.run([sys.executable, str(Path(__file__).with_name("plot.py")), str(root)],
                            stdout=log, stderr=subprocess.STDOUT)
status.update(status="complete" if result.returncode == 0 else "failed",
              exit_code=result.returncode, finished_utc=utc_now())
write_json(root / "completion-renderer.json", status)
