"""Ship the selected Symbolica kernel's authoritative Python stub unchanged.

Run with --check in documentation validation, or without it after updating the
locked kernel revision. Cargo uses the project's configured source for locked inputs.
"""

import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    host = Path(__file__).resolve().parent
    metadata = json.loads(
        subprocess.check_output(
            [
                "cargo",
                "metadata",
                "--locked",
                "--format-version=1",
                "--manifest-path",
                str(host / "Cargo.toml"),
            ],
            text=True,
        )
    )
    kernels = [
        package for package in metadata["packages"] if package["name"] == "symbolica"
    ]
    if len(kernels) != 1:
        raise SystemExit("the Community host must select exactly one Symbolica kernel")
    kernel = kernels[0]
    source = Path(kernel["manifest_path"]).parent / "symbolica.pyi"
    files = {
        host / "python/symbolica/core.pyi": source.read_bytes(),
        host / "python/symbolica/py.typed": b"",
    }
    for path, expected in files.items():
        if args.check:
            if not path.is_file() or path.read_bytes() != expected:
                raise SystemExit(
                    f"{path} differs from the locked Symbolica stub; run {__file__}"
                )
        else:
            path.write_bytes(expected)
    print(f"Community core typing matches {kernel['source']}")


if __name__ == "__main__":
    main()
