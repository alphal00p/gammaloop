"""Source-specific one-literal migration preparation, never a physics repair."""

import argparse
import collections
import datetime
import json
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import common as M

sys.path.append(str(M.FROZEN))
import bind_runtime as BIND
import build as FBUILD

CASE = M.CASES[0]
OLD = M.WORKSPACE / "target/xs-qed-l2"
OLD_ELF = OLD / "bin/current/generation_oracle"
NEW_ELF = M.OUT / "bin/generation_oracle"
OBSERVER = M.OUT / "bin/oracle_xs_qed_l2_migration_capture"
EXEC = M.HERE.parent / "generation-oracle-family-migration/runtime_exec.py"
OLD_HEX, NEW_HEX = "0xf36c3111407873cc", "0xf889b37c1fd93143"
REVIEW = {
    "A": (
        ".local/review-xs-qed-l2-a-opus-20261030T063905Z",
        "f5173b1d86a44c18d11179b13a0542eb0282097cb303079d709ab75d2184c554",
    ),
    "B": (
        ".local/review-xs-qed-l2-b-opus-20261030T063831Z",
        "5c97912c4b6fa9660c55b3bf03b785b37b36f2d3a042d677acd388f12ec8346b",
    ),
}


def manifest(paths):
    return {
        str(p): {"sha256": M.sha(p), "bytes": p.stat().st_size}
        for p in sorted(set(paths))
    }


def immutable(rev, path):
    p, _ = M.run(
        ["jj", "file", "show", "--ignore-working-copy", "-r", rev, path],
        M.OUT / "sources" / rev / path,
    )
    M.require(p.returncode == 0, "source", path)
    return p.stdout


def preserve(phase):
    if phase == "before":
        paths = set(M.load(OLD / "preservation/after.json")["files"])
        paths.update(M.load(OLD / "inventory/artifacts.json")["files"])
        for root in [M.FROZEN, *[M.WORKSPACE / r[0] for r in REVIEW.values()]]:
            paths.update(
                str(p)
                for p in root.rglob("*")
                if p.is_file()
                and not any(x in p.parts for x in ("cache", "tmp", "state"))
            )
        paths.update(
            str(p)
            for p in (
                OLD / "freeze.json",
                OLD / "inventory/sources.json",
                OLD / "inventory/artifacts.json",
                OLD / "preservation/after.json",
            )
        )
        result = manifest(Path(p) for p in paths)
    else:
        result = manifest(
            Path(p) for p in M.load(M.OUT / "preservation/before.json")["files"]
        )
        M.require(
            result == M.load(M.OUT / "preservation/before.json")["files"],
            "preservation",
            "all prior relied-on bytes unchanged",
        )
    M.write(
        M.OUT / f"preservation/{phase}.json",
        {"utc": datetime.datetime.now(datetime.UTC).isoformat(), "files": result},
    )
    print(phase, len(result), flush=True)


def derive():
    d = BIND.D
    d.weight = BIND.weight
    source = immutable(M.APPROVED, M.MAIN)
    historical = immutable("2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94", M.MAIN)
    expected = d.row(source, CASE)
    M.require(
        expected == d.row(historical, CASE), "old-row", "complete immutable Rows equal"
    )
    rows = {}
    for route in ("historical", "current"):
        root = OLD / f"runtime/{route}"
        p = M.load(root / f"captures/{CASE}/classes-None-t1/payload.json")
        rows[route] = d.summary(p, (root / f"dumps/{CASE}.census").read_text())
    M.require(rows["historical"] == expected, "derive", "historical complete Row")
    M.require(
        {k: v for k, v in rows["current"].items() if k != "numerators"}
        == {k: v for k, v in expected.items() if k != "numerators"},
        "derive",
        "all other fields",
    )
    M.require(
        f"0x{expected['numerators']:016x}" == OLD_HEX
        and f"0x{rows['current']['numerators']:016x}" == NEW_HEX,
        "derive",
        "exact independent integers",
    )
    M.write(
        M.OUT / "derivation.json",
        {
            "expected": expected,
            "rows": rows,
            "source_sha256": M.sha(
                M.OUT / "sources" / M.APPROVED / M.MAIN / "stdout.txt"
            ),
            "dependencies": manifest(
                [
                    M.FROZEN / "bind_runtime.py",
                    M.HERE.parent / "generation-oracle-family-migration/derive.py",
                ]
            ),
        },
    )
    print(rows, flush=True)


def audit():
    old = immutable(M.APPROVED, M.MAIN)
    new = (M.WORKSPACE / M.MAIN).read_text()
    M.require(
        old.count(OLD_HEX) == 1 and new == old.replace(OLD_HEX, NEW_HEX, 1),
        "textual-diff",
        "exact single substitution; every other byte identical",
    )
    M.write(
        M.OUT / "code-audit.json",
        {
            "base": M.APPROVED,
            "path": M.MAIN,
            "old_sha256": M.sha(M.OUT / "sources" / M.APPROVED / M.MAIN / "stdout.txt"),
            "new_sha256": M.sha(M.WORKSPACE / M.MAIN),
            "old": OLD_HEX,
            "new": NEW_HEX,
            "exact_substitution": True,
        },
    )


def prepare():
    for label, (root, expected) in REVIEW.items():
        file = M.WORKSPACE / root / "review.json"
        M.require(M.sha(file) == expected, "approval", label)
        M.write(
            M.OUT / "approvals" / (label + ".json"),
            {
                "path": str(file),
                "sha256": expected,
                "verdict": M.load(file)["verdict"],
                "notes": M.load(file)["nonblocking"],
            },
        )
    M.require(
        M.sha(OLD_ELF)
        == "7da2b7619edb99cc215e298b9dcb4616ce7eab5fea3195174a9692ad705dfd3c",
        "old-elf",
        "approved exact hash",
    )
    root = M.OUT / "source"
    root.mkdir(parents=True, exist_ok=False)
    archive = M.OUT / "source.tar"
    argv = [
        "git",
        "--git-dir=/shared/lcnbr/dev/gammaloop/main/.git",
        "archive",
        M.APPROVED,
    ]
    # Binary archive is not decoded through the text command recorder.
    start = datetime.datetime.now(datetime.UTC).isoformat()
    with archive.open("wb") as out:
        p = subprocess.run(
            argv, stdout=out, stderr=subprocess.PIPE, env=M.M.clean_env(), check=False
        )
    M.require(p.returncode == 0, "archive", "immutable source")
    M.write(
        M.OUT / "archive-command.json",
        {
            "argv": argv,
            "started_utc": start,
            "ended_utc": datetime.datetime.now(datetime.UTC).isoformat(),
            "exit": p.returncode,
            "stderr": p.stderr.decode(),
            "archive_sha256": M.sha(archive),
            "both_license_variables": "unset",
        },
    )
    with tarfile.open(archive) as t:
        t.extractall(root, filter="data")
    target = root / M.MAIN
    M.require(
        target.read_text().replace(OLD_HEX, NEW_HEX, 1)
        == (M.WORKSPACE / M.MAIN).read_text(),
        "overlay",
        "one-literal source",
    )
    shutil.copy2(M.WORKSPACE / M.MAIN, target)
    # Observer is a separate process. Only pool selection changes, from a
    # checked explicit environment value shared with the assertion harness.
    raw = (M.FROZEN / "observer.rs.txt").read_text()
    adapted = (
        raw.replace("pool1 only", "pools 1/8/32 observational only")
        .replace(
            'let folder = root.join(case.name).join(format!("classes-{classes:?}-t1"));',
            'let pool: usize = env::var("FEYNKIT_ORACLE_THREADS").unwrap().parse().unwrap();\n    assert!([1, 8, 32].contains(&pool));\n    let folder = root.join(case.name).join(format!("classes-{classes:?}-t{pool}"));',
        )
        .replace("Run::default().threads(1)", "Run::default().threads(pool)")
    )
    M.require(
        adapted != raw and "-t1" not in adapted and "threads(1)" not in adapted,
        "observer-overlay",
        "exact pool parameterization",
    )
    overlay = M.OUT / "overlay/oracle_xs_qed_l2_migration_capture.rs.txt"
    overlay.parent.mkdir(parents=True)
    overlay.write_text(adapted)
    M.write(
        M.OUT / "overlay/source.json",
        {
            "base": M.APPROVED,
            "source": str(M.FROZEN / "observer.rs.txt"),
            "source_sha256": M.sha(M.FROZEN / "observer.rs.txt"),
            "overlay_sha256": M.sha(overlay),
            "limit": "separate source-qualified observer, no assertion instrumentation",
        },
    )


def build():
    root = M.OUT / "source"
    for name, example in (
        ("generation_oracle", False),
        ("oracle_xs_qed_l2_migration_capture", True),
    ):
        if example:
            shutil.copy2(
                M.OUT / "overlay" / (name + ".rs.txt"),
                root / "crates/feynkit-generator/examples" / (name + ".rs"),
            )
        inputs = [
            p
            for p in root.rglob("*")
            if p.is_file() and "target" not in p.relative_to(root).parts
        ]
        dest = M.OUT / "build" / name
        M.write(dest / "sources-before.json", manifest(inputs))
        cargo = [
            "cargo",
            "build" if example else "test",
            "--locked",
            "--profile",
            "dev-optim",
            "-p",
            "feynkit-generator",
            "--example" if example else "--test",
            name,
            *([] if example else ["--no-run"]),
            "--message-format",
            "json-render-diagnostics",
        ]
        p, command = M.run(
            ["nix", "develop", M.URI, "-c", *cargo],
            dest / "execution",
            cwd=root,
            inputs=inputs,
        )
        M.require(p.returncode == 0, "build", name)
        messages = [json.loads(s) for s in p.stdout.splitlines() if s.startswith("{")]
        items = [
            r
            for r in messages
            if r.get("reason") == "compiler-artifact"
            and r["target"]["name"] == name
            and r.get("executable")
        ]
        M.require(
            len(items) == 1 and not items[0]["fresh"],
            "build",
            "actual new emitted executable",
        )
        item = items[0]
        elf = Path(item["executable"])
        saved = M.OUT / "bin" / name
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(elf, saved)
        kind = "example" if example else "test-integration-test"
        data = FBUILD.fingerprints(saved, name, kind)
        data.update(
            elf=str(saved),
            elf_sha256=M.sha(saved),
            command=command,
            emitted=item,
            source_base=M.APPROVED,
        )
        kernel_msgs = [
            r
            for r in messages
            if r.get("reason") == "compiler-artifact"
            and r["target"]["name"] == "symbolica"
        ]
        M.require(
            kernel_msgs and all("1a1829b" in r["manifest_path"] for r in kernel_msgs),
            "kernel",
            "actual source checkout from compiler outputs",
        )
        data["kernel_build_messages"] = kernel_msgs
        M.write(dest / "qualification.json", data)
        M.write(dest / "sources-after.json", manifest(inputs))
        M.require(
            M.load(dest / "sources-before.json") == M.load(dest / "sources-after.json"),
            "source",
            "build input stability",
        )
        print(name, M.sha(saved), data["closure_units"], flush=True)
    old = M.load(OLD / "build/current/generation_oracle/qualification.json")
    new = M.load(M.OUT / "build/generation_oracle/qualification.json")
    obs = M.load(M.OUT / "build/oracle_xs_qed_l2_migration_capture/qualification.json")
    projection = lambda q: sorted(
        (
            x["unit"].rsplit("-", 1)[0],
            x["kind"],
            x["features"],
            x["profile"],
            x["rustc"],
            str(x["rustflags"]),
        )
        for x in q["closure"]
    )
    M.require(
        projection(old) == projection(new) == projection(obs),
        "closure",
        "actual complete dependency feature/profile/compiler units",
    )
    M.require(
        new["closure"] == obs["closure"],
        "closure",
        "identical linked observer/test dependency units",
    )
    M.write(
        M.OUT / "build/closure-comparison.json",
        {
            "old_units": old["closure_units"],
            "new_units": new["closure_units"],
            "observer_units": obs["closure_units"],
            "full_projection_equal": True,
            "new_observer_exact_units_equal": True,
            "old_exact_units_equal": old["closure"] == new["closure"],
            "inputs": manifest(
                [
                    OLD / "build/current/generation_oracle/qualification.json",
                    M.OUT / "build/generation_oracle/qualification.json",
                    M.OUT
                    / "build/oracle_xs_qed_l2_migration_capture/qualification.json",
                ]
            ),
        },
    )
    for name, args in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "--version"])):
        p, _ = M.run(
            ["nix", "develop", M.URI, "-c", *args],
            M.OUT / "build" / name,
            cwd=Path("/tmp"),
        )
        M.require(p.returncode == 0, "toolchain", name)


def runpath(elf, label):
    p, _ = M.run(
        [
            "/nix/store/31gnwsgag5qa14hfblqi5hyn71pyx7py-binutils-2.46/bin/readelf",
            "-d",
            elf,
        ],
        M.OUT / "runpath" / label,
        inputs=[elf],
    )
    M.require(p.returncode == 0, "runpath", str(elf))
    paths = [
        x
        for line in p.stdout.splitlines()
        if "RUNPATH" in line or "RPATH" in line
        for x in re.search(r"\[(.*?)\]", line)[1].split(":")
    ]
    private = [x for x in paths if not x.startswith("/nix/store/")]
    M.require(
        private and all(not Path(x).exists() for x in private), "runpath", str(private)
    )
    for workspace in (
        M.WORKSPACE,
        M.WORKSPACE.parent / "generation-oracle-family-historical",
    ):
        M.require(
            not (workspace / "outputs/out/lib").exists(), "runpath", str(workspace)
        )
    return {
        "elf": str(elf),
        "sha256": M.sha(elf),
        "runpath": paths,
        "non_store_candidates_absent": private,
    }


def qualify():
    """Qualify retained successful Cargo outputs after the first helper's
    example-name fingerprint lookup failed; never call that helper a pass.
    """
    for name, kind in (
        ("generation_oracle", "test-integration-test"),
        ("oracle_xs_qed_l2_migration_capture", "example"),
    ):
        dest = M.OUT / "build" / name
        command = M.load(dest / "execution/command.json")
        M.require(
            command["exit"] == 0
            and command["inputs_before"] == command["inputs_after"],
            "retained-build",
            name,
        )
        before = M.load(dest / "sources-before.json")
        after = manifest(Path(p) for p in before)
        M.require(before == after, "retained-build", "actual full sources invariant")
        M.write(dest / "sources-after-qualification.json", after)
        messages = [
            json.loads(s)
            for s in (dest / "execution/stdout.txt").read_text().splitlines()
            if s.startswith("{")
        ]
        item = next(
            r
            for r in messages
            if r.get("reason") == "compiler-artifact"
            and r["target"]["name"] == name
            and r.get("executable")
        )
        saved = M.OUT / "bin" / name
        M.require(
            not item["fresh"] and M.sha(saved) == M.sha(item["executable"]),
            "retained-build",
            "actual compiled ELF",
        )
        data = FBUILD.fingerprints(saved, name, kind)
        kernels = [
            r
            for r in messages
            if r.get("reason") == "compiler-artifact"
            and r["target"]["name"] == "symbolica"
        ]
        M.require(
            kernels and all("1a1829b" in r["manifest_path"] for r in kernels),
            "kernel",
            "compiler output source",
        )
        data.update(
            elf=str(saved),
            elf_sha256=M.sha(saved),
            command=command,
            emitted=item,
            kernel_build_messages=kernels,
            source_base=M.APPROVED,
        )
        raw = {}
        for path, digest in data["fingerprint_files"].items():
            M.require(M.sha(path) == digest, "fingerprints", path)
            file = Path(path)
            to = M.OUT / "fingerprints" / name / file.parent.name / file.name
            to.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(file, to)
            raw[str(to)] = M.sha(to)
        data["retained_raw_fingerprints"] = raw
        M.write(dest / "qualification-v2.json", data)
    old = M.load(OLD / "build/current/generation_oracle/qualification.json")
    M.require(
        M.sha(OLD_ELF) == old["artifact_sha256"], "old-provenance", "exact approved ELF"
    )
    for path, digest in old["fingerprint_files"].items():
        M.require(
            M.sha(path) == digest,
            "old-provenance",
            "actual inherited fingerprint unchanged",
        )
    new = M.load(M.OUT / "build/generation_oracle/qualification-v2.json")
    obs = M.load(
        M.OUT / "build/oracle_xs_qed_l2_migration_capture/qualification-v2.json"
    )
    projection = lambda q: sorted(
        (
            x["unit"].rsplit("-", 1)[0],
            x["kind"],
            x["features"],
            x["profile"],
            x["rustc"],
            str(x["rustflags"]),
        )
        for x in q["closure"]
    )
    M.require(
        projection(old) == projection(new) == projection(obs),
        "closure",
        "complete actual dependency-unit projections",
    )
    M.require(
        new["closure"] == obs["closure"],
        "closure",
        "identical linked compiled dependency units",
    )
    for field in (
        "rustc",
        "profile",
        "features",
        "rustflags",
        "compile_kind",
        "target",
        "config",
    ):
        M.require(old["root"][field] == new["root"][field], "assertion-root", field)
    M.write(
        M.OUT / "build/closure-comparison.json",
        {
            "old_units": old["closure_units"],
            "new_units": new["closure_units"],
            "observer_units": obs["closure_units"],
            "full_projection_equal": True,
            "new_observer_exact_units_equal": True,
            "old_exact_units_equal": old["closure"] == new["closure"],
            "inputs": manifest(
                [
                    OLD / "build/current/generation_oracle/qualification.json",
                    M.OUT / "build/generation_oracle/qualification-v2.json",
                    M.OUT
                    / "build/oracle_xs_qed_l2_migration_capture/qualification-v2.json",
                ]
            ),
            "limits": "root test/build profiles differ; full dependency-unit equality is separately checked, not inferred from root features",
        },
    )
    for name, argv in (("rustc", ["rustc", "-vV"]), ("cargo", ["cargo", "--version"])):
        p, _ = M.run(
            ["nix", "develop", M.URI, "-c", *argv],
            M.OUT / "build" / name,
            cwd=Path("/tmp"),
        )
        M.require(p.returncode == 0, "toolchain", name)
    print(
        "qualified retained successful Cargo builds; original helper failure not credited",
        flush=True,
    )


def sandbox(elf, args, root, flags):
    root.mkdir(parents=True, exist_ok=False)
    for name in ("cache", "tmp", "state", "cwd", "dumps", "captures"):
        (root / name).mkdir()
    guard = runpath(elf, root.name + "-" + root.parent.name)
    argv = [
        shutil.which("bwrap"),
        "--ro-bind",
        "/",
        "/",
        "--dev",
        "/dev",
        "--proc",
        "/proc",
        "--bind",
        root,
        root,
        "--bind",
        root / "tmp",
        "/tmp",
        "--bind",
        root / "cache",
        "/home/lcnbr/.cache",
        "--unshare-net",
        "--die-with-parent",
        "--chdir",
        root / "cwd",
        "--",
        sys.executable,
        "-B",
        EXEC,
        elf,
        *args,
    ]
    p, command = M.run(
        argv,
        root / "execution",
        env_extra={
            "TMPDIR": str(root / "tmp"),
            "XDG_CACHE_HOME": str(root / "cache"),
            "GL_TEST_LOG_DIR": str(root / "state"),
            **flags,
        },
        inputs=[
            elf,
            EXEC,
            M.HERE / "migrate.py",
            M.FROZEN / "bind_runtime.py",
            M.HERE.parent / "generation-oracle-family-migration/derive.py",
        ],
    )
    M.write(root / "runpath.json", guard)
    return p, command


def selected(p):
    rows = re.findall(
        r"test result: (ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out",
        p.stdout,
    )
    M.require(len(rows) == 1, "selection", "one actual top-level summary")
    return rows[0]


def runtime():
    records = []
    for route, elf in (("old", OLD_ELF), ("new", NEW_ELF)):
        for pool in (1, 8, 32):
            root = M.OUT / f"runtime/{route}/pool{pool}"
            p, cmd = sandbox(
                elf,
                [CASE, "--exact", "--test-threads", "1", "--nocapture"],
                root,
                {
                    "FEYNKIT_ORACLE_THREADS": str(pool),
                    "FEYNKIT_ORACLE_DUMP": str(root / "dumps"),
                },
            )
            row = selected(p)
            M.require(
                list(map(int, row[1:5]))
                == ([0, 1, 0, 0] if route == "old" else [1, 0, 0, 0])
                and p.returncode == (101 if route == "old" else 0),
                "runtime",
                f"{route}/{pool}: {row}",
            )
            if route == "old":
                M.require(
                    f"{CASE}: summary differs from its oracle row" in p.stderr,
                    "old-boundary",
                    "Row numerator assertion",
                )
            records.append(
                {
                    "route": route,
                    "pool": pool,
                    "summary": row,
                    "exit": p.returncode,
                    "command": cmd,
                }
            )
            M.write(M.OUT / "runtime/results.json", records)
            obsroot = M.OUT / f"observers/{route}/pool{pool}"
            q, _ = sandbox(
                OBSERVER,
                [obsroot / "captures", CASE],
                obsroot,
                {"FEYNKIT_ORACLE_THREADS": str(pool)},
            )
            M.require(q.returncode == 0, "observer", f"{route}/{pool}")
            print(route, pool, row, flush=True)


def binding():
    consumed = [M.OUT / "derivation.json", OLD / f"runtime/current/dumps/{CASE}.census"]
    for mode in ("None", "Auto"):
        consumed.extend(
            p
            for p in (
                OLD / f"runtime/current/captures/{CASE}/classes-{mode}-t1"
            ).iterdir()
            if p.is_file()
        )
    for route in ("old", "new"):
        for pool in (1, 8, 32):
            consumed.extend(
                p
                for p in (M.OUT / f"runtime/{route}/pool{pool}/dumps").iterdir()
                if p.is_file()
            )
            for mode in ("None", "Auto"):
                folder = (
                    M.OUT
                    / f"observers/{route}/pool{pool}/captures/{CASE}/classes-{mode}-t{pool}"
                )
                consumed.extend(
                    folder / n
                    for n in (
                        "payload.json",
                        "typed-owners.json",
                        "model.json",
                        "request.json",
                        "options.json",
                        "grouping.json",
                        "capture.json",
                        "prepared-report.txt",
                        "built-report.txt",
                    )
                )
    inputs_before = manifest(consumed)
    d = BIND.D
    d.weight = BIND.weight
    rows, modes = [], []
    expected = M.load(M.OUT / "derivation.json")["rows"]["current"]
    for route in ("old", "new"):
        for pool in (1, 8, 32):
            compact_file = M.OUT / f"runtime/{route}/pool{pool}/dumps/{CASE}.census"
            M.require(
                compact_file.read_bytes()
                == (OLD / f"runtime/current/dumps/{CASE}.census").read_bytes(),
                "complete-census",
                "all actual compact entries/report equal approved current census",
            )
            payloads = {}
            for mode in ("None", "Auto"):
                base = (
                    M.OUT
                    / f"observers/{route}/pool{pool}/captures/{CASE}/classes-{mode}-t{pool}"
                )
                dump = (
                    M.OUT
                    / f"runtime/{route}/pool{pool}/dumps/{CASE}.classes-{mode}-t{pool}.json"
                )
                sources = dump.with_name(dump.name[:-5] + ".sources.json")
                M.require(
                    dump.read_bytes()
                    == sources.read_bytes()
                    == (base / "payload.json").read_bytes(),
                    "observer-test",
                    "complete bytes",
                )
                p = M.load(dump)
                typed = M.load(base / "typed-owners.json")
                owners = sum(len(g["owners"]) for g in typed)
                axes = sum(len(o["ordered_ports"]) for g in typed for o in g["owners"])
                M.require(
                    owners == 540
                    and axes == 1296
                    and [g["id"] for g in typed] == [g["id"] for g in p["diagrams"]],
                    "actual-owner-domain",
                    "ordered complete capture",
                )
                payloads[mode] = p
                approved = OLD / f"runtime/current/captures/{CASE}/classes-{mode}-t1"
                # No mode normalization; each full independent mode binds to
                # its own certified domain, including all physical/port fields.
                for file in (
                    "payload.json",
                    "typed-owners.json",
                    "model.json",
                    "request.json",
                    "options.json",
                    "grouping.json",
                    "prepared-report.txt",
                    "built-report.txt",
                ):
                    M.require(
                        (base / file).read_bytes() == (approved / file).read_bytes(),
                        "approved-domain",
                        f"{route}/{pool}/{mode}/{file}",
                    )
                meta = M.load(base / "capture.json")
                M.require(
                    meta["case"] == CASE
                    and meta["classes"] == mode
                    and meta["grouping"] == "None"
                    and meta["prepared_len"] == meta["built_len"] == 30
                    and meta["run"]
                    == f"Run {{ threads: Some({pool}), parallel: true, check_heredity: false, cancelled: false }}",
                    "observer-Run",
                    str(base),
                )
                M.require(
                    all(
                        g["overall_factor"] == g["payload"]["overall_factor"]
                        for g in p["diagrams"]
                    ),
                    "within-factor",
                    str(base),
                )
                compact = (
                    M.OUT / f"runtime/{route}/pool{pool}/dumps/{CASE}.census"
                ).read_text()
                summary = d.summary(p, compact)
                target = (
                    expected
                    if mode == "None"
                    else {**expected, "named": 14666857465192919060}
                )
                M.require(summary == target, "full-Row", str(base))
                weights = collections.Counter(
                    str(d.weight(g["overall_factor"])) for g in p["diagrams"]
                )
                M.require(
                    weights == {"1": 17, "-1": 13},
                    "weights",
                    "exact signed rational distribution",
                )
                rows.append(
                    {
                        "route": route,
                        "pool": pool,
                        "mode": mode,
                        "owners": owners,
                        "axes": axes,
                        "Row": summary,
                        "weights": dict(weights),
                        "capture": meta,
                        "inputs": manifest(
                            [
                                dump,
                                sources,
                                *[
                                    base / n
                                    for n in (
                                        "payload.json",
                                        "typed-owners.json",
                                        "model.json",
                                        "request.json",
                                        "options.json",
                                        "grouping.json",
                                        "capture.json",
                                        "built-report.txt",
                                        "prepared-report.txt",
                                    )
                                ],
                            ]
                        ),
                    }
                )
            key = lambda g: json.dumps(g["key"], sort_keys=True, separators=(",", ":"))
            a, b = payloads["None"], payloads["Auto"]
            bk = [key(g) for g in b["diagrams"]]
            ak = [key(g) for g in a["diagrams"]]
            M.require(
                collections.Counter(ak) == collections.Counter(bk)
                and len(set(ak)) == 30,
                "mode-domain",
                "full key multiplicity",
            )
            mapping = []
            for i, g in enumerate(a["diagrams"]):
                j = bk.index(key(g))
                other = b["diagrams"][j]
                M.require(
                    M.C.strict_equal(
                        {k: v for k, v in g.items() if k != "name"},
                        {k: v for k, v in other.items() if k != "name"},
                    ),
                    "mode-payload",
                    str(i),
                )
                for payload, graph, index in ((a, g, i), (b, other, j)):
                    M.require(
                        payload["groups"][index]
                        == {
                            "master": index,
                            "members": [
                                {
                                    "id": graph["id"],
                                    "name": graph["name"],
                                    "ratio": "1",
                                    "overall_factor": graph["overall_factor"],
                                    "source_index": index,
                                }
                            ],
                        },
                        "singleton-chain",
                        str(index),
                    )
                mapping.append(
                    {
                        "None": i,
                        "Auto": j,
                        "id": g["id"],
                        "key": g["key"],
                        "names": [g["name"], other["name"]],
                        "groups": [a["groups"][i], b["groups"][j]],
                    }
                )
            moved = {
                str(x["None"]): x["Auto"] for x in mapping if x["None"] != x["Auto"]
            }
            M.require(
                moved
                == {
                    "12": 13,
                    "13": 14,
                    "14": 15,
                    "15": 16,
                    "16": 17,
                    "17": 18,
                    "18": 19,
                    "19": 21,
                    "21": 12,
                },
                "mode-cycle",
                str(moved),
            )
            modes.append(
                {
                    "route": route,
                    "pool": pool,
                    "mapping": mapping,
                    "moved": moved,
                    "byte_equal": False,
                }
            )
    # Full Summary in each old assertion, not just panic text or digest hashes.
    for pool in (1, 8, 32):
        stderr = (M.OUT / f"runtime/old/pool{pool}/execution/stderr.txt").read_text()
        for prefix, target in (
            ("  left", expected),
            (" right", M.load(M.OUT / "derivation.json")["expected"]),
        ):
            text = re.search(
                "^" + prefix + r": (Summary \{.*\})$", stderr, re.MULTILINE
            )[1]
            text = (
                re.sub(
                    r"Fraction \{ numerator: (-?\d+), denominator: (\d+) \}",
                    r"(\1, \2)",
                    text,
                )
                .replace("symmetry_factors: [", "symmetry_factors: &[")
                .replace("Summary {", "Row {", 1)
            )
            M.require(
                d.row("observed: " + text, "observed") == target,
                "assertion-summary",
                str(pool) + prefix,
            )
    M.write(
        M.OUT / "binding.json",
        {
            "domains": rows,
            "mode_correspondence": modes,
            "inputs_before": inputs_before,
            "inputs_after": manifest(consumed),
            "limit": "pool1 certified; 8/32 current runtime observations only; separate observer/test join, not in-memory assertion metadata",
        },
    )
    M.require(
        inputs_before == manifest(consumed),
        "binding-inputs",
        "all consumed captures/dumps unchanged",
    )


def negatives():
    rows = []
    for label, elf, name, flags in (
        ("zero-selection", NEW_ELF, "xs_qed_l2_absent_selection", {}),
        ("record-mode", OLD_ELF, CASE, {"FEYNKIT_ORACLE_RECORD": "1"}),
    ):
        p, command = sandbox(
            elf,
            [name, "--exact", "--test-threads", "1", "--nocapture"],
            M.OUT / "negatives" / label,
            {"FEYNKIT_ORACLE_THREADS": "1", **flags},
        )
        row = selected(p)
        M.require(
            p.returncode == 0
            and int(row[1]) == (0 if label == "zero-selection" else 1),
            "negative-setup",
            label,
        )
        if label == "record-mode":
            M.require("Row {" in p.stdout, "record-mode", "actual bypass observed")
        rows.append(
            {
                "name": label,
                "summary": row,
                "command": command,
                "admitted": False,
                "credit": "none: zero selection or explicit record-mode bypass",
            }
        )
    M.write(M.OUT / "negatives/results.json", rows)


def replay():
    # Immutable pins/required capture closure are checked before portable calls;
    # no registry/inventory guard is disabled or represented as passing.
    pins = M.load(OLD / "inventory/sources.json")
    mismatches = [
        p for p, rec in pins["files"].items() if M.sha(M.WORKSPACE / p) != rec["sha256"]
    ]
    M.write(
        M.OUT / "replay/outer-source-inventory-check.json",
        {
            "mismatches": mismatches,
            "disposition": "reject; zero credit" if mismatches else "match",
            "qualification": "source inventory check, not full frozen outer validator; portable finite entrypoint below pins mathematical inputs separately",
        },
    )
    M.require(
        mismatches == ["docs/developers.toml"],
        "outer-inventory",
        "only authorized new registry insertion",
    )
    for path, rec in pins["files"].items():
        if path == "docs/developers.toml":
            continue  # Not a portable semantic input; full inventory rejected above.
        M.require(
            M.sha(M.WORKSPACE / path) == rec["sha256"], "certificate-source-pin", path
        )
    M.write(
        M.OUT / "replay/source-inventory-binding.json",
        {
            "path": str(OLD / "inventory/sources.json"),
            "sha256": M.sha(OLD / "inventory/sources.json"),
            "inventory": pins,
        },
    )
    for mode in ("None", "Auto"):
        dest = M.OUT / "replay" / mode
        root = M.OUT / "replay-sandbox" / mode
        root.mkdir(parents=True)
        (root / "cache").mkdir()
        (root / "tmp").mkdir()
        guards = [
            runpath(elf, "native-" + mode + "-" + route)
            for route, elf in M.C.B.ELFS.items()
        ]
        inputs = [
            p
            for r in ("historical", "current")
            for p in (OLD / f"runtime/{r}/captures/{CASE}/classes-{mode}-t1").iterdir()
            if p.is_file()
        ]
        inputs += (
            list(M.C.B.ELFS.values())
            + [p for p in M.FROZEN.iterdir() if p.is_file()]
            + list(M.C.source_hashes())
        )
        argv = [
            "bwrap",
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--bind",
            M.OUT,
            M.OUT,
            "--bind",
            root / "tmp",
            "/tmp",
            "--bind",
            root / "cache",
            "/home/lcnbr/.cache",
            "--unshare-net",
            "--die-with-parent",
            "--",
            sys.executable,
            "-B",
            M.FROZEN / "validate.py",
            "--mode",
            mode,
            "--out",
            dest,
        ]
        p, _ = M.run(argv, root / "execution", inputs=inputs)
        M.require(
            p.returncode == 0
            and M.load(dest / "verdict.json")["disposition"] == "accept",
            "replay",
            mode,
        )
        for name in (
            "cross-pairs.json",
            "aggregate-pairs.json",
            "faithfulness.json",
            "anchor-pairs.json",
            "cross-ports.json",
            "definitions.json",
        ):
            M.require(
                (dest / name).read_bytes()
                == (OLD / "proof" / mode / name).read_bytes(),
                "replay-pairs",
                name,
            )
        M.write(dest / "runpath-guards.json", guards)
        print("replay", mode, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "stage",
        choices=(
            "before",
            "derive",
            "audit",
            "prepare",
            "build",
            "qualify",
            "runtime",
            "binding",
            "negatives",
            "replay",
            "after",
        ),
    )
    ap.add_argument("--execute", action="store_true")
    args = ap.parse_args()
    stage = args.stage
    if not args.execute:
        stamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%S%f")
        p, _ = M.run(
            [sys.executable, "-B", Path(__file__), stage, "--execute"],
            M.OUT / "commands" / (stamp + "-" + stage),
        )
        print(p.stdout, end="", flush=True)
        print(p.stderr, end="", file=sys.stderr, flush=True)
        sys.exit(p.returncode)
    if stage in ("before", "after"):
        preserve(stage)
    else:
        globals()[stage]()


if __name__ == "__main__":
    main()
