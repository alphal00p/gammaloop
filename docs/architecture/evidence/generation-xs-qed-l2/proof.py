"""No-network proof subprocesses with actual input and checker ELF bindings."""

import argparse
import shutil
import sys

import boundary as B
import common as M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kind", choices=("baseline", "controls"))
    args = ap.parse_args()
    root = M.OUT / ("proof-sandbox-" + args.kind)
    root.mkdir(parents=True, exist_ok=False)
    (root / "cache").mkdir()
    (root / "tmp").mkdir()
    commands = []
    if args.kind == "baseline":
        commands = [
            (
                mode,
                [
                    sys.executable,
                    "-B",
                    M.HERE / "validate.py",
                    "--mode",
                    mode,
                    "--out",
                    M.OUT / "proof" / mode,
                ],
            )
            for mode in ("None", "Auto")
        ]
    else:
        commands = [("controls", [sys.executable, "-B", M.HERE / "controls.py"])]
    for label, command in commands:
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
            M.OUT,
            M.OUT,
            "--bind",
            root / "cache",
            "/home/lcnbr/.cache",
            "--bind",
            root / "tmp",
            "/tmp",
            "--unshare-net",
            "--die-with-parent",
            "--",
            *command,
        ]
        inputs = [
            *B.hashes(),
            *M.C.B.ELFS.values(),
            M.HERE / "schema.json",
            M.HERE / "validate.py",
            M.HERE / "boundary.py",
            *M.C.source_hashes(),
        ]
        p, _ = M.run(argv, root / label, inputs=inputs)
        M.require(
            p.returncode == 0,
            "proof",
            label + ": inspect exact verdict and stop if new obstruction",
        )
        print(args.kind, label, p.returncode, flush=True)


if __name__ == "__main__":
    main()
