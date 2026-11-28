"""Command line: ``python -m hepkit_renorm SPEC [options]`` runs :func:`hepkit_renorm.driver.renormalise`.

``SPEC`` is a specification file (``.toml``) or the name of a built-in specification (``sm_one_family``). The
outputs (tree-level and renormalised model files, the order store, ``run_report.json`` / ``run_report.txt``) are
written to ``--output``. The exit status is 0 when the run finished and every check passed, 1 when a check failed,
2 when the run stopped with an error (the report says where).
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from .driver import OrderNotSupported, RunOptions, renormalise


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m hepkit_renorm", description=__doc__.split("\n\n")[0])
    ap.add_argument("spec", help="specification file (.toml) or built-in model name (e.g. sm_one_family)")
    ap.add_argument("-o", "--output", default="renormalised", help="output directory (default: ./renormalised)")
    ap.add_argument("-l", "--order", type=int, default=1, help="target loop order (default 1; only 1 is supported)")
    ap.add_argument("-j", "--processes", type=int, default=os.cpu_count() or 1,
                    help="worker processes for the graph evaluation (default: all cores)")
    ap.add_argument("--cache", default=None, help="per-graph cache directory (default: <output>/graph_cache)")
    ap.add_argument("--no-cache", action="store_true", help="do not cache graph values")
    ap.add_argument("--normalisation", choices=("full", "minimal"), default="full",
                    help="normalisation conditions of the restoring counterterm (default: full)")
    ap.add_argument("--fermion-kinetic", action="store_true",
                    help="representative with fermion kinetic restoring counterterms")
    ap.add_argument("--no-fallback", action="store_true", help="do not relax a normalisation condition")
    ap.add_argument("--no-pole-identity", action="store_true",
                    help="skip the check of the ghost-number-1 poles (the most expensive loop pass)")
    ap.add_argument("--cohomology", action="store_true", help="also compute the cohomology dimension of b (slow)")
    ap.add_argument("--schema", choices=("v1", "v2"), default="v1", help="schema of the written models")
    ap.add_argument("--seed", type=int, default=0, help="seed of every random choice (default 0)")
    ap.add_argument("--name", default=None, help="base name of the output files (default: the model name)")
    ap.add_argument("-q", "--quiet", action="store_true", help="no progress lines")
    args = ap.parse_args(argv)
    start = time.time()

    def log(text: str) -> None:
        print(f"[{time.time() - start:8.1f} s] {text}", flush=True)

    opts = RunOptions(output=args.output, order=args.order, cache=False if args.no_cache else args.cache,
                      processes=args.processes, normalisation=args.normalisation,
                      fermion_kinetic=args.fermion_kinetic, fallback=not args.no_fallback,
                      pole_identity=not args.no_pole_identity, cohomology=args.cohomology, schema=args.schema,
                      seed=args.seed, name=args.name, log=None if args.quiet else log)
    try:
        result = renormalise(args.spec, opts)
    except OrderNotSupported as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:   # noqa: BLE001 -- any failure: the report is written, say where it is
        print(f"error: {type(exc).__name__}: {exc}\nsee {os.path.join(args.output, 'run_report.txt')}",
              file=sys.stderr)
        return 2
    print(result.report.text())
    return 0 if result.report.passed else 1


if __name__ == "__main__":
    sys.exit(main())
