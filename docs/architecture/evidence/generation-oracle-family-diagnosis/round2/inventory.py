"""Supplemental evidence closure; hashes never decide semantic equality."""

import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess

SPEC = importlib.util.spec_from_file_location(
    "certificate", Path(__file__).with_name("validate.py")
)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
OUT = Path(
    os.environ.get(
        "FAMILY_ROUND2_OUT", str(V.WORKSPACE / "target/oracle-numerator-family-round2")
    )
)


def sha(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main():
    roots = [
        V.BASE,
        V.HIST_WORKSPACE / "target/oracle-numerator-family-diagnosis",
        V.WORKSPACE / "target/oracle-numerator-families-preflight",
        V.WORKSPACE / "target/oracle-family-independent-review-r1",
    ]
    rows = []
    for root in roots:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                rows.append(
                    {
                        "path": str(path),
                        "sha256": sha(path),
                        "bytes": path.stat().st_size,
                        "symlink": str(path.readlink()) if path.is_symlink() else None,
                    }
                )
    destination = OUT / "relied-on-frozen-evidence.json"
    if destination.exists():
        V.require(
            V.load(destination) == rows, "frozen-inventory", "relied-on files changed"
        )
    else:
        V.write(destination, rows)
    guards = []
    for root, revision in (
        (V.WORKSPACE, "e814b6b4c3463d9aa1c018939365b624619af24b"),
        (V.HIST_WORKSPACE, "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94"),
    ):
        paths = [
            "Cargo.lock",
            "Cargo.toml",
            "nix-ci.nix",
            "crates/feynkit-generator/tests/generation_oracle/main.rs",
            "crates/feynkit-generator/tests/generation_oracle/classes.rs",
            "crates/feynkit-generator/tests/generation_oracle/census.rs",
            "crates/feynkit-generator/tests/generation_oracle/cases.rs",
        ]
        for path in paths:
            command = [
                "jj",
                "file",
                "show",
                "--ignore-working-copy",
                "-r",
                revision,
                path,
            ]
            proc = subprocess.run(command, cwd=root, capture_output=True, check=True)
            V.require(
                (root / path).read_bytes() == proc.stdout,
                "source-invariance",
                str(root / path),
            )
            guards.append(
                {
                    "workspace": str(root),
                    "revision": revision,
                    "path": path,
                    "sha256": sha(root / path),
                }
            )
    V.write(OUT / "source-invariance.json", guards)
    print("frozen evidence files", len(rows), "source guards", len(guards))


if __name__ == "__main__":
    main()
