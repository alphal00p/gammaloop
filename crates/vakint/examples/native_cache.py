"""Small one-loop prepare/load/full-snapshot example for a coherent HepKit host.

No numerical master computation or implicit network access. This example's
runtime validation is recorded separately; it does not certify four-loop prep.
"""

import argparse
import json
from pathlib import Path

from symbolica import E
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prepared", type=Path)
    parser.add_argument("snapshot", type=Path)
    args = parser.parse_args()
    if args.prepared.exists() or args.snapshot.exists():
        parser.error("use two new directories; existing caches are never overwritten")
    obj = Vakint(
        evaluation_order=[
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        form_exe_path="/native-cache-example-must-not-use-form",
        integral_normalization_factor="FMFTandMATAD",
        use_dot_product_notation=True,
    )

    def progress(report):
        print(
            f"{report['status']}: {report['completed_roots']}/{report['total_roots']} durable roots"
        )

    report = obj.prepare_cache(1, str(args.prepared), progress=progress)
    print(
        json.dumps(
            {
                key: report[key]
                for key in ("status", "completed_roots", "charged_disk_bytes")
            }
        )
    )
    if report["status"] != "complete":
        raise RuntimeError(f"partial preparation: {report['failure']}")
    obj.load_cache(1, str(args.prepared))
    prepared = E("topo(I1L(4,3))", default_namespace="vakint")
    first = obj.evaluate_integral(prepared)
    assert (first - obj.evaluate_integral(prepared)).expand().together() == E("0")
    held_out = E("topo(I1L(4,5))", default_namespace="vakint")
    expected = obj.evaluate_integral(held_out)
    saved = obj.save_cache(1, str(args.snapshot), mode="full-snapshot")
    print(
        json.dumps(
            {
                key: saved[key]
                for key in ("complete", "durable_records", "charged_disk_bytes")
            }
        )
    )
    if not saved["complete"]:
        raise RuntimeError(
            f"partial snapshot: {saved['failure']} / {saved['checkpoint_failure']}"
        )
    obj.clear_cache(1)
    obj.load_cache(1, str(args.snapshot))
    assert (expected - obj.evaluate_integral(held_out)).expand().together() == E("0")
    print("exact_parity=passed cache=", obj.cache_statistics(1)["native"])


if __name__ == "__main__":
    main()
