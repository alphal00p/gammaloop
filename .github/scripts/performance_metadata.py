"""Record public CI and Nix identities without serializing the environment."""

import argparse
import hashlib
import json
import os
import subprocess
import tomllib
from datetime import UTC, datetime
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("identity", "builds", "finish"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--identity", type=Path)
    args = parser.parse_args()
    if args.command == "identity":
        value = {
            "repository": os.environ["GITHUB_REPOSITORY"],
            "repository_url": "https://github.com/" + os.environ["GITHUB_REPOSITORY"],
            "framework_revision": os.environ["FRAMEWORK_REVISION"],
            "candidate": os.environ["CANDIDATE"],
            "baseline": os.environ["BASELINE"],
            "event": os.environ["SOURCE_EVENT"],
            "run_id": int(os.environ["GITHUB_RUN_ID"]),
            "run_attempt": int(os.environ["GITHUB_RUN_ATTEMPT"]),
            "pr": int(os.environ["PR_NUMBER"]) if os.environ.get("PR_NUMBER") else None,
        }
    elif args.command == "builds":
        value = {}
        for side in ("baseline", "candidate"):
            package = os.environ[side.upper() + "_PACKAGE"]
            info = json.loads(
                subprocess.check_output(["nix", "path-info", "--json", package])
            )
            manifest = tomllib.loads((Path(side) / "Cargo.toml").read_text())
            value[side] = {
                "package": package,
                "store": info,
                "context": {
                    "flake_lock_sha256": digest(Path(side) / "flake.lock"),
                    "cargo_lock_sha256": digest(Path(side) / "Cargo.lock"),
                    "profiles": manifest.get("profile", {}),
                    "system": "x86_64-linux",
                },
            }
    else:
        if args.output.exists():
            return
        if not args.identity or not args.identity.exists():
            raise ValueError("Source identity was not recorded")
        value = {
            "schema_version": 1,
            "ci": json.loads(args.identity.read_text()),
            "cases": [],
            "status": 2,
            "created_at": datetime.now(UTC).isoformat(),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
