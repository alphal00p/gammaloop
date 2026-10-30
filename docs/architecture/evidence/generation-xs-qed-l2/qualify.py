"""Source, build and cut/owner-domain qualifications derived from actual captures."""

import collections
import hashlib

import boundary as B
import common as M


def main():
    source_rows = []
    for route, revision in M.REVISIONS.items():
        for path in (
            "crates/feynkit-model/src/model/builtins.rs",
            "crates/feynkit-graph/src/expressions.rs",
            "crates/feynkit-generator/src/generation.rs",
            "crates/feynkit-graph/src/symbols.rs",
        ):
            text, command = M.immutable(revision, path)
            source_rows.append(
                {
                    "route": route,
                    "revision": revision,
                    "path": path,
                    "sha256": hashlib.sha256(text.encode()).hexdigest(),
                    "command": command,
                }
            )
    # Exact production delta: current3432 consists of originalb9 production
    # with the four approved expectation edits plus docs; none affects this Row.
    p, command = M.run(
        [
            "jj",
            "diff",
            "--ignore-working-copy",
            "--from",
            "b9f788922341e7bde2e801a054423ae70ed7e9ae",
            "--to",
            M.APPROVED,
            "--summary",
            "--no-pager",
        ],
        M.OUT / "qualification/current-base-delta",
    )
    paths = [line[2:] for line in p.stdout.splitlines() if line]
    M.require(
        p.returncode == 0
        and all(path.startswith("docs/") or path == M.MAIN for path in paths),
        "production-delta",
        "current3432vsb9",
    )
    projection = lambda q: sorted(
        (
            r["unit"].rsplit("-", 1)[0],
            r["kind"],
            r["features"],
            r["profile"],
            r["rustc"],
            str(r["rustflags"]),
        )
        for r in q["closure"]
    )
    builds = {}
    for route in M.REVISIONS:
        a, b = (
            M.load(M.OUT / "build" / route / name / "qualification.json")
            for name in ("generation_oracle", "oracle_xs_qed_l2_capture")
        )
        M.require(
            projection(a) == projection(b),
            "dependency-features",
            route + ": actual within-route full closure",
        )
        M.require(a["root"]["rustc"] == b["root"]["rustc"], "compiler", route)
        builds[route] = {
            "dependency_units": a["closure_units"],
            "test_profile": a["root"]["profile"],
            "observer_profile": b["root"]["profile"],
            "full_dependency_projection_equal": True,
            "compiler": a["root"]["rustc"],
            "test_elf": a["artifact_sha256"],
            "observer_elf": b["artifact_sha256"],
            "kernel": "c78cffc20f318df1fce50539a587cbe8ffc717c0"
            if route == "historical"
            else "1a1829b1f3d0a88b53b44085b3b9f34ef77945b3",
        }
    domains = []
    for route, root in B.ROOTS.items():
        for mode in ("None", "Auto"):
            folder = root / M.CASES[0] / f"classes-{mode}-t1"
            model = M.load(folder / "model.json")
            p = M.load(folder / "payload.json")
            typed = M.load(folder / "typed-owners.json")
            interaction = collections.Counter()
            cut_sizes = collections.Counter()
            interface = 0
            for graph in p["diagrams"]:
                payload = graph["payload"]
                for vi, vertex in enumerate(payload["vertices"]):
                    interaction[
                        model["vertex_rules"][vertex["interaction"]]["name"]
                    ] += 1
                for edge in payload["edges"]:
                    M.require(len(edge) == 2, "edge-tuple", "positional")
                    for side in ("source", "target"):
                        index = edge[0][side]
                        M.require(
                            index is None
                            or type(index) is int
                            and 0 <= index < len(payload["vertices"]),
                            "incidence",
                            route,
                        )
                for cut in payload["cuts"]:
                    cut_sizes[len(cut["cut"])] += 1
                    for port in cut["cut"]:
                        M.require(
                            type(port["edge"]) is int
                            and 0 <= port["edge"] < len(payload["edges"])
                            and port["endpoint"] in ("source", "target"),
                            "cut-owner",
                            route,
                        )
                    for side in ("left", "right"):
                        for port in cut[side]["half_edges"]:
                            M.require(
                                type(port["edge"]) is int
                                and 0 <= port["edge"] < len(payload["edges"])
                                and port["endpoint"] in ("source", "target"),
                                "cut-side-owner",
                                route,
                            )
                interface += sum(
                    "num_interface" in v for v in payload["vertices"]
                ) + sum("num_interface" in e[1] for e in payload["edges"])
            domains.append(
                {
                    "route": route,
                    "mode": mode,
                    "graphs": len(p["diagrams"]),
                    "owners": sum(len(g["owners"]) for g in typed),
                    "axes": sum(
                        len(o["ordered_ports"]) for g in typed for o in g["owners"]
                    ),
                    "cut_size_histogram": dict(cut_sizes),
                    "rule_occurrences": dict(interaction),
                    "local_interfaces_present": interface,
                    "provenance_limit": "physical cut/rule values are observed and cross-route compared, not independently derived physical correctness",
                }
            )
    M.write(
        M.OUT / "qualification/source-domain.json",
        {
            "sources": source_rows,
            "current_base_delta": command,
            "builds": builds,
            "domains": domains,
            "estimate_correction": "pre-capture20 owners estimate corrected by observed forward6vertices9edges+3globals=18 owners/graph,540 total; no budget inferred",
            "within_route_binding": "separate observer joined to pristine assertion ELF through byte-identical complete census payload/source companions; typed port/expression typing and locality remain capture facts",
        },
    )
    print(builds)


if __name__ == "__main__":
    main()
