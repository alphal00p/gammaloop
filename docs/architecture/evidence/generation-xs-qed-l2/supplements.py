"""Case-specific cut/schema controls and explicit infrastructure separation."""

import sys

import boundary as B
import common as M


def main():
    rows = []
    specs = {
        "cut-endpoint": "cross-physical",
        "cut-side": "cross-physical",
        "bool-incidence": "schema",
        "unknown-owner-path": "schema",
        "malformed-json": "json",
        "request-both": "selected-mode",
        "interface-value": "interface-live-link",
        "missing-input": "infrastructure",
        "missing-license": "infrastructure",
    }
    for mode in ("None", "Auto"):
        for name, stage in specs.items():
            root = M.OUT / "supplements" / mode / name
            roots = {}
            for route, source in B.ROOTS.items():
                roots[route] = root / route
                d = roots[route] / M.CASES[0] / f"classes-{mode}-t1"
                d.mkdir(parents=True, exist_ok=False)
                for file in (source / M.CASES[0] / f"classes-{mode}-t1").iterdir():
                    (d / file.name).symlink_to(file)
            current = roots["current"] / M.CASES[0] / f"classes-{mode}-t1"
            if name in (
                "cut-endpoint",
                "cut-side",
                "bool-incidence",
                "unknown-owner-path",
                "malformed-json",
                "interface-value",
            ):
                p = M.load(current / "payload.json")
                g = p["diagrams"][0]["payload"]
                if name == "cut-endpoint":
                    g["cuts"][0]["cut"][0]["endpoint"] = (
                        "source"
                        if g["cuts"][0]["cut"][0]["endpoint"] == "target"
                        else "target"
                    )
                elif name == "cut-side":
                    g["cuts"][0]["left"]["loop_count"] += 1
                elif name == "bool-incidence":
                    g["edges"][0][0]["source"] = False
                elif name == "unknown-owner-path":
                    for value in roots.values():
                        path = value / M.CASES[0] / f"classes-{mode}-t1/payload.json"
                        other = M.load(path)
                        other["diagrams"][0]["payload"]["hidden_interface"] = []
                        path.unlink()
                        M.write(path, other)
                elif name == "interface-value":
                    next(v for v in g["vertices"] if "num_interface" in v)[
                        "num_interface"
                    ].reverse()
                if name != "unknown-owner-path":
                    (current / "payload.json").unlink()
                    M.write(current / "payload.json", p)
                    if name == "malformed-json":
                        with (current / "payload.json").open("a") as stream:
                            stream.write("{}")
            if name == "request-both":
                for value in roots.values():
                    path = value / M.CASES[0] / f"classes-{mode}-t1/request.json"
                    request = M.load(path)
                    request["classes"] = "Auto" if mode == "None" else "None"
                    path.unlink()
                    M.write(path, request)
            if name == "missing-input":
                roots["current"] = root / "absent"
            argv = [
                sys.executable,
                "-B",
                M.HERE / "validate.py",
                "--mode",
                mode,
                "--out",
                root / "result",
            ]
            for route, value in roots.items():
                argv += ["--" + route, value, "--anchors-" + route, value]
            env = {}
            if name == "missing-license":
                (root / "empty-config").mkdir()
                env["XDG_CONFIG_HOME"] = str(root / "empty-config")
            p, command = M.run(argv, root / "execution", env_extra=env)
            verdict = M.load(root / "result/verdict.json")
            expected = "harness-error" if stage == "infrastructure" else "reject"
            matched = (
                verdict["stage"] == stage
                and verdict["disposition"] == expected
                and p.returncode == (3 if stage == "infrastructure" else 2)
            )
            rows.append(
                {
                    "name": name,
                    "mode": mode,
                    "expected_stage": stage,
                    "observed": verdict,
                    "matched": matched,
                    "command": command,
                }
            )
            M.write(M.OUT / "supplements/results.json", rows)
            print(mode, name, matched, verdict["stage"], flush=True)
    M.require(
        all(r["matched"] for r in rows), "supplement", "unexpected predicate; no credit"
    )


if __name__ == "__main__":
    main()
