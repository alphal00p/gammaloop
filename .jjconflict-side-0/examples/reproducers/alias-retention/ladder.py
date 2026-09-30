"""Compare full spacetime reductions with and without generated alias definitions."""

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("case", choices=("gluon", "quark"))
parser.add_argument("--hep-examples", type=Path, required=True)
parser.add_argument("--save", type=Path, required=True)
parser.add_argument("--compare", type=Path)
parser.add_argument("--no-aliases", action="store_true")
parser.add_argument(
    "--symbol-indices",
    action="store_true",
    help="Diagnose graph-label eligibility after trace collection",
)
parser.add_argument("--max-steps-per-domain", type=int, default=16)
args = parser.parse_args()
sys.path.insert(0, str(args.hep_examples))

from generated_ladder_helpers import ladder_filter, project_and_split
from symbolica import E, Expression, S, core
from symbolica.community import hep
from symbolica.community import tensor as sp

phases = []


def measured(name, operation):
    wall, cpu = time.perf_counter(), time.process_time()
    value = operation()
    phases.append(
        {
            "name": name,
            "wall": time.perf_counter() - wall,
            "cpu": time.process_time() - cpu,
        }
    )
    print(json.dumps(phases[-1]), flush=True)
    return value


quark = args.case == "quark"
model = hep.Model.standard_model()
dimension = S("ladder_measure::D")
diagram = measured(
    "generation",
    lambda: (
        model.process(["g"], ["g"], vertex_allow=["V_76" if quark else "V_36"])
        .generate_diagrams(
            loops=4,
            max_vertices=8,
            allow_self_loops=False,
            self_energy=None,
            tadpoles=None,
            zero_snails=None,
            numerator_grouping=None,
            progress=None,
            filter=ladder_filter(quark),
        )
        .diagrams[0]
    ),
)
raw = measured(
    "numerator",
    lambda: model.expand_couplings(
        diagram.numerator_expression(in_lmb=not quark).to_expression()
    ).replace(model.particle("b").mass, E("0")),
)
numerator = measured(
    "admission", lambda: sp.TensorExpression(raw).with_lorentz_dimension(dimension)
)
_, _, spacetime = measured(
    "projection_split", lambda: project_and_split(numerator, dimension)
)
if args.symbol_indices:
    assert quark, "this diagnostic selects the quark trace's four dummy labels"
    collected = measured("collect_trace", spacetime.contract)
    raw = collected.root.to_expression()
    labels = [
        x
        for x in raw.get_all_indeterminates()
        if x.get_name().split("::")[-1] == "hedge"
    ]
    assert len(labels) == 4, labels
    for i, label in enumerate(labels):
        raw = raw.replace(label, S(f"alias_probe::mu{i}"))
    spacetime = sp.TensorExpression(raw)
settings = {
    "gamma": quark,
    "color": False,
    "epsilon": False,
    "contract": "dots",
    "max_steps_per_domain": args.max_steps_per_domain,
}
result = measured("simplify_algebra", lambda: spacetime.simplify_algebra(**settings))
aliases = len(result.aliases)
print(
    f"status={result.reduction_status} aliases={aliases} bytes={result.get_byte_size()}",
    flush=True,
)
assert result.reduction_status == sp.ReductionStatus.Complete, str(
    result.reduction_status
)
if args.no_aliases:
    assert aliases == 0
measured("completed_rerun", lambda: result.simplify_algebra(**settings))
scalar = measured("expand", lambda: result.expand().to_expression())
if args.compare:
    assert scalar == Expression.load(str(args.compare)), "unequal exact polynomials"
scalar.save(str(args.save))
record = {
    "case": args.case,
    "symbol_indices": args.symbol_indices,
    "settings": settings,
    "phases": phases,
    "complete": True,
    "aliases": aliases,
    "bytes": result.get_byte_size(),
    "terms": len(list(scalar.terms())),
    "exact_comparison": str(args.compare) if args.compare else None,
    "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    "extension": str(core.__file__),
    "extension_sha256": hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
}
args.save.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record), flush=True)
