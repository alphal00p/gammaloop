"""Native-parser layer controls, separate from full-validator syntax controls."""

import shutil
import sys

import boundary as B
import common as M


def main():
    rows = []
    source = B.ROOTS["current"] / M.CASES[0] / "classes-None-t1/typed-owners.json"
    expression = M.load(source)[0]["owners"][0]["producer_decoded_expression"]
    for label, malformed in (
        ("trailing-semicolon", expression + ";"),
        ("trailing-term-junk", expression + "+)"),
    ):
        for route, elf in M.C.B.ELFS.items():
            d = M.OUT / "native-parser-controls" / label / route
            d.mkdir(parents=True, exist_ok=False)
            (d / "cache").mkdir()
            (d / "tmp").mkdir()
            pairs = d / "pairs.json"
            M.write(pairs, [{"historical": malformed, "current": expression}])
            runner = (
                M.HERE.parent / "generation-oracle-family-migration/runtime_exec.py"
            )
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
                d,
                d,
                "--bind",
                d / "cache",
                "/home/lcnbr/.cache",
                "--bind",
                d / "tmp",
                "/tmp",
                "--unshare-net",
                "--die-with-parent",
                "--",
                sys.executable,
                "-B",
                runner,
                elf,
                pairs,
                d / "result.json",
            ]
            p, command = M.run(
                argv, d / "execution", inputs=[source, elf, pairs, runner]
            )
            M.require(
                p.returncode == 101
                and "oracle_family_exact.rs:35:" in p.stderr
                and "Err" in p.stderr,
                "native-parser-control",
                label + route + ": unexpected failure (not inequality credit)",
            )
            rows.append(
                {
                    "label": label,
                    "kernel": route,
                    "expected": "native parse rejection",
                    "stage": "direct-native-parser",
                    "exit": p.returncode,
                    "command": command,
                    "pair_sha256": M.sha(pairs),
                    "limit": "direct native-layer control only; full validator rejects these strings earlier at syntax, separately tested",
                }
            )
    M.write(M.OUT / "native-parser-controls/results.json", rows)


if __name__ == "__main__":
    main()
