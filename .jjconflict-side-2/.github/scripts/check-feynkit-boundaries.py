#!/usr/bin/env python3
"""Reject dependencies from standalone FeynKit crates into GammaLoop."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN_PACKAGES = {"gammaloop-api", "gammalooprs"}
FORBIDDEN_SOURCE_REFERENCES = (
    "gammaloop-api",
    "gammaloop_api",
    "gammalooprs",
)

# GammaLoop may re-export and extend the canonical model types, but it must not
# redefine them or hide them behind a local alias/newtype.  Include every
# identifier owned by ``feynkit-model`` so a second ID namespace cannot emerge
# on either side of the GammaLoop API boundary.
FEYNKIT_MODEL_OWNER_TYPE = (
    r"(?:Model|Particle|OrderId|ParameterId|ParticleId|PropagatorId|"
    r"LorentzStructureId|CouplingId|VertexRuleId|ModelFunctionId|"
    r"ModelFormFactorId)"
)
FEYNKIT_MODEL_OWNER_PATH = (
    rf"(?:[A-Za-z_][A-Za-z0-9_]*::)*{FEYNKIT_MODEL_OWNER_TYPE}"
)
FORBIDDEN_MODEL_OWNER_WRAPPERS = (
    rf"\b(?:pub(?:\([^)]*\))?\s+)?(?:struct|enum|type)\s+{FEYNKIT_MODEL_OWNER_TYPE}\b",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w+\s*\(\s*(?:pub(?:\([^)]*\))?\s+)?{FEYNKIT_MODEL_OWNER_PATH}\s*\)\s*;",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w+(?:\s*<[^>{{;]*>)?\s*\{{\s*(?:pub(?:\([^)]*\))?\s+)?\w+\s*:\s*{FEYNKIT_MODEL_OWNER_PATH}\s*,?\s*\}}",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?enum\s+\w+(?:\s*<[^>{{;]*>)?\s*\{{\s*\w+\s*\(\s*{FEYNKIT_MODEL_OWNER_PATH}\s*\)\s*,?\s*\}}",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?type\s+\w+\s*=\s*{FEYNKIT_MODEL_OWNER_PATH}\s*;",
)

FEYNKIT_FINALIZED_OWNER_TYPE = (
    r"(?:FeynmanDiagram|GenerationResult|GenerationOptions|CffGraph|"
    r"CffExpression|CffResult|SurfaceCache|EnergySurface|HSurface)"
)
FEYNKIT_FINALIZED_OWNER_PATH = (
    rf"(?:[A-Za-z_][A-Za-z0-9_]*::)*{FEYNKIT_FINALIZED_OWNER_TYPE}"
)
FORBIDDEN_FINALIZED_OWNER_WRAPPERS = (
    rf"\b(?:pub(?:\([^)]*\))?\s+)?(?:struct|enum|type)\s+{FEYNKIT_FINALIZED_OWNER_TYPE}\b",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w+\s*\(\s*(?:pub(?:\([^)]*\))?\s+)?{FEYNKIT_FINALIZED_OWNER_PATH}\s*\)\s*;",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w+(?:\s*<[^>{{;]*>)?\s*\{{\s*(?:pub(?:\([^)]*\))?\s+)?\w+\s*:\s*{FEYNKIT_FINALIZED_OWNER_PATH}\s*,?\s*\}}",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?enum\s+\w+(?:\s*<[^>{{;]*>)?\s*\{{\s*\w+\s*\(\s*{FEYNKIT_FINALIZED_OWNER_PATH}\s*\)\s*,?\s*\}}",
    rf"\b(?:pub(?:\([^)]*\))?\s+)?type\s+\w+\s*=\s*{FEYNKIT_FINALIZED_OWNER_PATH}\s*;",
)

# GammaLoop consumes these semantic objects from FeynKit.  Defining local
# stand-ins makes the dependency graph look correct while silently restoring a
# second model, diagram-generation pipeline, or CFF IR.  Keep this list narrow:
# numerical evaluator data and the final runtime graph intentionally remain in
# ``gammalooprs``.
FORBIDDEN_GAMMALOOP_DEFINITIONS = {
    "crates/gammalooprs/src": (
        *FORBIDDEN_FINALIZED_OWNER_WRAPPERS,
        *FORBIDDEN_MODEL_OWNER_WRAPPERS,
    ),
    "crates/gammaloop-api/src": (
        *FORBIDDEN_FINALIZED_OWNER_WRAPPERS,
        *FORBIDDEN_MODEL_OWNER_WRAPPERS,
    ),
    "crates/gammalooprs/src/model": (
        r"\b(?:pub\s+)?struct\s+SerializableModel\b",
        r"\b(?:pub\s+)?type\s+Arc(?:Particle|VertexRule|Propagator)\b",
        r"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w*(?:Adapter|Wrapper)\b",
    ),
    "crates/gammalooprs/src/feyngen": (
        r"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w*(?:Adapter|Wrapper)\b",
    ),
    "crates/gammalooprs/src/cff": (
        r"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+(?:CFFExpression|CffExpression|CffGraph|CffResult|SurfaceCache|Esurface|EnergySurface|Hsurface|HSurface)\b",
        r"\b(?:pub(?:\([^)]*\))?\s+)?type\s+(?:CFFExpression|CffExpression|CffGraph|CffResult|OrientationID|EsurfaceID|EnergySurfaceId|HsurfaceID|HSurfaceId|HybridSurfaceID)\b",
        r"\b(?:pub(?:\([^)]*\))?\s+)?struct\s+\w*(?:Adapter|Wrapper|CffInput)\b",
    ),
}

FORBIDDEN_GAMMALOOP_GENERATION_PATTERNS = (
    r"\b(?:fn\s+)?generate_numerators\b",
    r"\b(?:fn\s+)?setup_loop_momentum_basis\b",
    r"\bfn\s+fix_cp_vertex_rules\b",
)

# Surface construction and shift canonicalization are part of FeynKit's CFF
# engine. GammaLoop may evaluate canonical surfaces through extension traits,
# but it must not grow a second structural implementation.
FORBIDDEN_GAMMALOOP_CFF_PATTERNS = (
    r"\bfn\s+new_from_subgraph\b",
    r"\bfn\s+canonicalize_shift\b",
    r"\bfn\s+(?:new_)?from_cut_(?:left|right|side)\b",
)

# Keep the public extension-trait entry point thin. The implementation it calls
# is checked separately below because it lives alongside the finalized DOT
# runtime-artifact importer.
FORBIDDEN_FEYNKIT_ENTRYPOINT_PATTERNS = (
    r"\bstruct\s+FeynkitGeneratorAdapter\b",
    r"\bstruct\s+FeynkitAdapterResult\b",
    r"\bfn\s+prepare_diagram\b",
    r"\bGraph::from_parsed\s*\(",
)

FEYNKIT_BRIDGE_UNITS = (
    "from_feynkit_diagram",
    "feynkit_external_edges",
    "feynkit_external_connections",
    "feynkit_tag_index",
    "feynkit_runtime_map",
    "feynkit_runtime_half_edge",
    "feynkit_runtime_cuts",
    "translate_feynkit_atom",
    "add_feynkit_sign",
    "feynkit_runtime_lmb",
    # These generic helpers are shared with the finalized-runtime importer, but
    # are part of the transitive path used by ``Graph::from_feynkit``. Keep
    # them under the same no-recanonicalization rule so physics work cannot be
    # hidden one call below the explicitly FeynKit-named bridge functions.
    "extract_initial_data",
    "process_cut_edges",
    "build_underlying_graph",
    "from_feynkit",
)
FEYNKIT_BRIDGE_HELPER_PATTERN = (
    r"(?m)^[ \t]*(?:pub(?:\([^)]*\))?[ \t]+)?fn[ \t]+(feynkit_[a-z0-9_]+)\b"
)

# The canonical bridge may translate IDs, sew explicitly paired cuts, and copy
# FeynKit's chosen routing into runtime storage. It must not call the adjacent
# finalized-DOT importer or any routine that chooses rules, numerators, or an
# alternative loop-momentum basis. Scanning named bridge functions rather than
# all of graph/parse intentionally leaves that strict runtime-artifact importer
# available for fully specified GammaLoop artifacts.
FORBIDDEN_MECHANICAL_BRIDGE_PATTERNS = (
    r"\b(?:Self|Graph|ParseGraph)::from_parsed(?:_with_validation)?\s*\(",
    r"\b(?:generate|setup)_numerators?\s*\(",
    r"\bsetup_loop_momentum_basis\s*\(",
    r"\bmaterialize_explicit_loop_momentum_basis\s*\(",
    r"\bgenerate_loop_momentum_bases(?:_of)?\s*\(",
    r"\blmb_impl\s*\(",
    r"\bfix_cp_vertex_rules\s*\(",
    r"\bextract_vertex_particles\s*\(",
    r"\bprocess_vertex_hedges\s*\(",
)


def cargo_metadata() -> dict[str, object]:
    command = [
        "cargo",
        "metadata",
        "--all-features",
        "--format-version=1",
        "--locked",
    ]
    result = subprocess.run(
        command,
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        sys.stderr.write(result.stderr)
        raise SystemExit(result.returncode)
    return json.loads(result.stdout)


def is_feynkit(package: dict[str, object]) -> bool:
    name = str(package["name"])
    return name == "feynkit" or name.startswith("feynkit-")


def dependency_violations(metadata: dict[str, object]) -> list[str]:
    packages = {package["id"]: package for package in metadata["packages"]}
    nodes = {node["id"]: node for node in metadata["resolve"]["nodes"]}
    roots = [package for package in packages.values() if is_feynkit(package)]
    if not roots:
        return ["cargo metadata did not contain any FeynKit packages"]

    violations = []
    for root in roots:
        pending = [(root["id"], [root["name"]])]
        visited = set()
        while pending:
            package_id, path = pending.pop()
            if package_id in visited:
                continue
            visited.add(package_id)
            for dependency_id in nodes[package_id]["dependencies"]:
                dependency = packages[dependency_id]
                dependency_path = [*path, dependency["name"]]
                if dependency["name"] in FORBIDDEN_PACKAGES:
                    violations.append(" -> ".join(dependency_path))
                else:
                    pending.append((dependency_id, dependency_path))
    return sorted(set(violations))


def source_violations(metadata: dict[str, object]) -> list[str]:
    violations = []
    for package in metadata["packages"]:
        if not is_feynkit(package):
            continue
        package_root = Path(package["manifest_path"]).parent
        for source in package_root.rglob("*.rs"):
            text = source.read_text(encoding="utf-8")
            for reference in FORBIDDEN_SOURCE_REFERENCES:
                if reference in text:
                    relative = source.relative_to(REPOSITORY_ROOT)
                    violations.append(f"{relative}: references {reference}")
    return sorted(violations)


def rust_function(source: str, name: str) -> str | None:
    """Return one Rust function, excluding adjacent functions in the same file."""
    declaration = re.search(
        rf"(?m)^[ \t]*(?:pub(?:\([^)]*\))?[ \t]+)?fn[ \t]+{re.escape(name)}\b",
        source,
    )
    if declaration is None:
        return None
    opening = source.find("{", declaration.end())
    if opening < 0:
        return None

    depth = 0
    index = opening
    while index < len(source):
        if source.startswith("//", index):
            newline = source.find("\n", index + 2)
            index = len(source) if newline < 0 else newline + 1
            continue
        if source.startswith("/*", index):
            comment_depth = 1
            index += 2
            while index < len(source) and comment_depth:
                if source.startswith("/*", index):
                    comment_depth += 1
                    index += 2
                elif source.startswith("*/", index):
                    comment_depth -= 1
                    index += 2
                else:
                    index += 1
            continue
        if source[index] == '"':
            index += 1
            while index < len(source):
                if source[index] == "\\":
                    index += 2
                elif source[index] == '"':
                    index += 1
                    break
                else:
                    index += 1
            continue
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[declaration.start() : index + 1]
        index += 1
    return None


def mechanical_bridge_violations(source: str, relative: Path) -> list[str]:
    violations = []
    discovered_units = re.findall(FEYNKIT_BRIDGE_HELPER_PATTERN, source)
    units = dict.fromkeys((*FEYNKIT_BRIDGE_UNITS, *discovered_units))
    for unit in units:
        function = rust_function(source, unit)
        if function is None:
            violations.append(f"{relative}: canonical bridge unit `{unit}` is missing")
            continue
        for pattern in FORBIDDEN_MECHANICAL_BRIDGE_PATTERNS:
            if re.search(pattern, function):
                violations.append(
                    f"{relative}:{unit}: repeats finalized FeynKit work ({pattern})"
                )
    return violations


def duplicate_ownership_violations() -> list[str]:
    violations = []
    for relative_root, patterns in FORBIDDEN_GAMMALOOP_DEFINITIONS.items():
        root = REPOSITORY_ROOT / relative_root
        for source in root.rglob("*.rs"):
            text = source.read_text(encoding="utf-8")
            for pattern in patterns:
                if re.search(pattern, text):
                    relative = source.relative_to(REPOSITORY_ROOT)
                    violations.append(
                        f"{relative}: duplicates a FeynKit-owned semantic type ({pattern})"
                    )

    # Generation is a semantic ownership boundary, not merely a constraint on
    # one adapter file.  Search the whole GammaLoop domain crate so moving a
    # duplicate implementation to a DOT parser or helper module cannot evade
    # this guard.
    gammaloop_root = REPOSITORY_ROOT / "crates/gammalooprs/src"
    for source in gammaloop_root.rglob("*.rs"):
        text = source.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_GAMMALOOP_GENERATION_PATTERNS:
            if re.search(pattern, text):
                relative = source.relative_to(REPOSITORY_ROOT)
                violations.append(
                    f"{relative}: repeats FeynKit-owned generation work ({pattern})"
                )

    cff_root = REPOSITORY_ROOT / "crates/gammalooprs/src/cff"
    for source in cff_root.rglob("*.rs"):
        text = source.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_GAMMALOOP_CFF_PATTERNS:
            if re.search(pattern, text):
                relative = source.relative_to(REPOSITORY_ROOT)
                violations.append(
                    f"{relative}: repeats FeynKit-owned structural CFF work ({pattern})"
                )

    bridge = REPOSITORY_ROOT / "crates/gammalooprs/src/feyngen/feynkit.rs"
    if bridge.exists():
        text = bridge.read_text(encoding="utf-8")
        for pattern in FORBIDDEN_FEYNKIT_ENTRYPOINT_PATTERNS:
            if re.search(pattern, text):
                relative = bridge.relative_to(REPOSITORY_ROOT)
                violations.append(
                    f"{relative}: repeats finalized FeynKit work ({pattern})"
                )

    mechanical_bridge = REPOSITORY_ROOT / "crates/gammalooprs/src/graph/parse/mod.rs"
    if not mechanical_bridge.exists():
        violations.append(
            "crates/gammalooprs/src/graph/parse/mod.rs: canonical FeynKit bridge is missing"
        )
    else:
        relative = mechanical_bridge.relative_to(REPOSITORY_ROOT)
        violations.extend(
            mechanical_bridge_violations(
                mechanical_bridge.read_text(encoding="utf-8"), relative
            )
        )
    return sorted(violations)


def self_test() -> None:
    source = """
fn from_feynkit() {
    copy_finalized_data();
}

fn from_parsed() {
    materialize_explicit_loop_momentum_basis();
}

fn feynkit_new_bridge_helper() {
    copy_finalized_data();
}

impl Graph {
    pub(crate) fn feynkit_associated_helper(&self) {
        copy_finalized_data();
    }
}
"""
    extracted = rust_function(source, "from_feynkit")
    if extracted is None or "copy_finalized_data" not in extracted:
        raise AssertionError("failed to extract the mechanical bridge function")
    if "materialize_explicit_loop_momentum_basis" in extracted:
        raise AssertionError("bridge extraction leaked into the finalized importer")
    if "feynkit_new_bridge_helper" not in re.findall(
        FEYNKIT_BRIDGE_HELPER_PATTERN, source
    ):
        raise AssertionError(
            "new FeynKit bridge helpers are not discovered automatically"
        )
    if "feynkit_associated_helper" not in re.findall(
        FEYNKIT_BRIDGE_HELPER_PATTERN, source
    ):
        raise AssertionError(
            "associated or visibility-qualified FeynKit bridge helpers are not discovered"
        )

    for root in ("crates/gammalooprs/src", "crates/gammaloop-api/src"):
        patterns = FORBIDDEN_GAMMALOOP_DEFINITIONS.get(root, ())
        for owner_wrapper in (
            "struct LocalModel(feynkit_model::Model);",
            "struct LocalParticle(Particle);",
            "struct NamedModelWrapper {\n    pub(crate) inner: feynkit_model::Model,\n}",
            "enum ParticleWrapper {\n    Canonical(Particle),\n}",
            "type LocalParticleId = gammalooprs::model::ParticleId;",
            "enum VertexRuleId { Local }",
            "struct DiagramAdapter(feynkit_graph::FeynmanDiagram);",
            "struct NamedCffWrapper {\n    inner: feynkit_cff::CffResult,\n}",
            "enum DiagramWrapper {\n    Finalized(FeynmanDiagram),\n}",
            "struct CffAdapter(CffResult);",
            "type LocalGenerationResult = feynkit_generator::GenerationResult;",
            "enum GenerationOptions { Local }",
        ):
            if not any(re.search(pattern, owner_wrapper) for pattern in patterns):
                raise AssertionError(
                    f"failed to reject `{owner_wrapper}` in the {root} boundary"
                )

    duplicate_surface_builder = "fn new_from_subgraph() {}"
    if not any(
        re.search(pattern, duplicate_surface_builder)
        for pattern in FORBIDDEN_GAMMALOOP_CFF_PATTERNS
    ):
        raise AssertionError("failed to reject structural CFF construction in GammaLoop")

    duplicate_cut_surface_builder = "fn new_from_cut_left() {}"
    if not any(
        re.search(pattern, duplicate_cut_surface_builder)
        for pattern in FORBIDDEN_GAMMALOOP_CFF_PATTERNS
    ):
        raise AssertionError("failed to reject cut-surface construction in GammaLoop")

    forbidden_source = source.replace(
        "copy_finalized_data();", "Graph::from_parsed(finalized_dot);"
    )
    function = rust_function(forbidden_source, "from_feynkit")
    if function is None or not any(
        re.search(pattern, function) for pattern in FORBIDDEN_MECHANICAL_BRIDGE_PATTERNS
    ):
        raise AssertionError("failed to reject physics work in the mechanical bridge")

    delegated_helpers = (
        "extract_initial_data",
        "process_cut_edges",
        "build_underlying_graph",
    )
    if not all(helper in FEYNKIT_BRIDGE_UNITS for helper in delegated_helpers):
        raise AssertionError("delegated bridge helpers are not guarded transitively")
    delegated_source = "\n".join(
        f"fn {unit}() {{ copy_finalized_data(); }}"
        for unit in dict.fromkeys(FEYNKIT_BRIDGE_UNITS)
    ).replace(
        "fn process_cut_edges() { copy_finalized_data(); }",
        "fn process_cut_edges() { materialize_explicit_loop_momentum_basis(); }",
    )
    delegated_errors = mechanical_bridge_violations(
        delegated_source, Path("synthetic/graph/parse/mod.rs")
    )
    if not any(
        ":process_cut_edges:" in error and "repeats finalized FeynKit work" in error
        for error in delegated_errors
    ):
        raise AssertionError("failed to reject physics work in a delegated bridge helper")


def main() -> None:
    metadata = cargo_metadata()
    dependency_errors = dependency_violations(metadata)
    source_errors = source_violations(metadata)
    ownership_errors = duplicate_ownership_violations()
    if dependency_errors or source_errors or ownership_errors:
        print(
            "FeynKit must remain independent of GammaLoop application crates.",
            file=sys.stderr,
        )
        for violation in [*dependency_errors, *source_errors, *ownership_errors]:
            print(f"  - {violation}", file=sys.stderr)
        raise SystemExit(1)
    print("FeynKit dependency boundary is valid.")


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-test"]:
        self_test()
        print("FeynKit boundary guard self-test passed.")
    elif sys.argv[1:]:
        raise SystemExit("usage: check-feynkit-boundaries.py [--self-test]")
    else:
        main()
