"""Run the retained native MathML producer and its Chromium client together."""

import os
import subprocess
import sys
import time
from pathlib import Path


def main():
    output = (
        Path(
            os.environ.get(
                "LINNET_MATHML_BROWSER_ROOT",
                Path.cwd() / ".local/graph-mathml-labels/browser",
            )
        )
        / "ci"
    )
    output.mkdir(parents=True)
    environment = dict(os.environ, LINNET_MATHML_BROWSER_DIR=str(output))
    with (output / "producer.log").open("w", encoding="utf-8") as log:
        producer = subprocess.Popen(
            [
                "cargo",
                "nextest",
                "run",
                "--archive-file",
                sys.argv[1],
                "--extract-to",
                ".",
                "--extract-overwrite",
                "--workspace-remap",
                ".",
                "--profile",
                "ci_gammaloop",
                "--run-ignored",
                "only",
                "-E",
                "test(=mathml_view_browser_producer)",
            ],
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            ready = output / "ready.json"
            deadline = time.monotonic() + 30
            while not ready.exists():
                if producer.poll() is not None or time.monotonic() >= deadline:
                    raise RuntimeError("native browser producer failed to become ready")
                time.sleep(0.05)
            subprocess.run(
                [
                    sys.executable,
                    "crates/linnet-render-py/tests/mathml_view_browser.py",
                    "--ready",
                    str(ready),
                    "--output",
                    str(output / "client"),
                    "--engines",
                    "chromium",
                ],
                check=True,
                timeout=230,
                env=environment,
            )
            if producer.wait(timeout=15) != 0:
                raise RuntimeError("native browser producer rejected the round trip")
            subprocess.run(
                [
                    sys.executable,
                    "crates/linnet-render-py/tests/iframe_display_browser.py",
                    "--html",
                    str(output / "chromium-350-a.html"),
                    "--html",
                    str(output / "chromium-350-b.html"),
                    "--output",
                    str(output / "iframe"),
                ],
                check=True,
                timeout=120,
                env=environment,
            )
        finally:
            if producer.poll() is None:
                producer.terminate()
                try:
                    producer.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    producer.kill()
                    producer.wait()
            print(log.name)
            print(Path(log.name).read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
