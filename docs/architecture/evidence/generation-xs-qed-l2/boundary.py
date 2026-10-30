"""New-domain inspection and source-qualified finite checker loading."""

import importlib.util
import sys

import common as M

ROUND3 = M.HERE.parent / "generation-oracle-family-diagnosis/round3"
ROOTS = {r: M.OUT / "runtime" / r / "captures" for r in M.REVISIONS}
FILES = (
    "model.json",
    "request.json",
    "options.json",
    "grouping.json",
    "payload.json",
    "typed-owners.json",
    "capture.json",
    "prepared-report.txt",
    "built-report.txt",
)


def checker():
    # Load unchanged approved mechanisms with their own import owner. No source
    # file is patched, nor is its four-case conclusion imported into this domain.
    previous = sys.modules["common"]
    sys.modules["common"] = M.C
    spec = importlib.util.spec_from_file_location(
        "xs_qed_finite_design", ROUND3 / "validate.py"
    )
    v = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(v)
    finally:
        sys.modules["common"] = previous
    v.B.CASES = M.CASES
    v.B.ROOTS = ROOTS
    v.SCHEMA = M.load(M.HERE / "schema.json")
    return v


def hashes(roots=ROOTS):
    return {
        str(p): M.sha(p)
        for root in roots.values()
        for p in root.rglob("*")
        if p.is_file()
    }


def schema():
    import collections

    types, keys, sources = (
        collections.defaultdict(set),
        collections.defaultdict(set),
        {},
    )

    def walk(value, path):
        types[path].add(type(value).__name__)
        if isinstance(value, dict):
            numeric = path.endswith("/edge_signatures")
            keys[path].add("numeric-map" if numeric else tuple(sorted(value)))
            for k, val in value.items():
                walk(val, path + "/" + ("#" if numeric else k))
        elif isinstance(value, list):
            for i, val in enumerate(value):
                walk(
                    val,
                    path + "/" + (str(i) if path.endswith("/payload/edges/*") else "*"),
                )

    for root in ROOTS.values():
        for mode in ("None", "Auto"):
            folder = root / M.CASES[0] / f"classes-{mode}-t1"
            for name in FILES[:6]:
                p = folder / name
                sources[str(p)] = M.sha(p)
                walk(M.load(p), name)
    result = {
        "scope": M.CASES[0],
        "types": {k: sorted(v) for k, v in sorted(types.items())},
        "object_keys": {
            k: sorted([list(v) if isinstance(v, tuple) else v for v in vals], key=str)
            for k, vals in sorted(keys.items())
        },
        "derivation_inputs": sources,
        "qualification": "independently observed fixed producer paths/kinds, not copied four-family schema; exact edge tuples and symmetric domains checked separately",
    }
    M.write(M.HERE / "schema.json", result)
    print("new supported shape", len(types), len(keys), "input files", len(sources))


if __name__ == "__main__":
    schema()
