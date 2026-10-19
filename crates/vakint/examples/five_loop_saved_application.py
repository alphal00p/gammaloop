"""Apply an externally hosted five-loop package to a native scalar vacuum DOT.

Requires one Symbolica assembly containing the matching FeynKit and Vakint crates.
Input powers must belong to the selected package's recorded application union.
"""

import argparse
from pathlib import Path
from time import perf_counter

from symbolica import E
from symbolica.community.feynkit import FeynmanDiagram, Kinematics, Model
from symbolica.community.hepkit.vakint import (
    Vakint,
    VakintEvaluationMethod,
    integral_from_diagram,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("application", type=Path)
    parser.add_argument("dot", type=Path)
    parser.add_argument(
        "--powers", help="comma-separated powers in native family order"
    )
    parser.add_argument("--cache-gib", type=int, default=1)
    parser.add_argument(
        "--cache-policy", choices=("evict-lru", "fail-on-limit"), default="evict-lru"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="print the complete guarded result"
    )
    args = parser.parse_args()

    model = Model(Path(__file__).resolve().parents[1] / "tests/fixtures/phi_3_4.json")
    diagram = FeynmanDiagram.from_topology_dot(model, args.dot.read_text())
    family = diagram.propagator_family(kinematics=Kinematics())
    powers = None if args.powers is None else [int(p) for p in args.powers.split(",")]
    integral = integral_from_diagram(diagram, family, E("1"), powers=powers)
    engine = Vakint(
        evaluation_order=[
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        form_exe_path="/five-loop-example-must-not-use-form",
        integral_normalization_factor="1",
        rustred_cache_bytes=args.cache_gib * 1024**3,
        rustred_cache_policy=args.cache_policy,
    )
    start = perf_counter()
    engine.load_rustred_application(str(args.application))
    print(f"load_seconds={perf_counter() - start:.6f}")
    previous = None
    for phase in ("first", "warm"):
        start = perf_counter()
        result = engine.evaluate_rustred_symbolic(integral.to_expression())
        print(f"{phase}_application_seconds={perf_counter() - start:.6f}")
        if previous is not None and result != previous:
            raise AssertionError("warm application changed the exact result")
        previous = result
        print("cache:", engine.rustred_application_cache_statistics())
    value, conditions, families = result
    print(f"nonzero_conditions={len(conditions)} families={len(families)}")
    if args.verbose:
        print("value:", value)
        print("nonzero_conditions:", conditions)
        print("families:", families)
    else:
        preview = str(value)
        print("value_preview:", preview[:240] + ("..." if len(preview) > 240 else ""))


if __name__ == "__main__":
    main()
