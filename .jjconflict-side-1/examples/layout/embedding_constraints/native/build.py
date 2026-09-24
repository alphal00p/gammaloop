"""Build the pinned OGDF SPQR primitive; downloaded sources stay under target/."""

import argparse
import hashlib
import io
import os
import shutil
import ssl
import subprocess
import tarfile
import urllib.request
from pathlib import Path

OGDF_COMMIT = "1a40505ac00677a58ab5ce532e1dcc3bd63116e5"
OGDF_SHA256 = "f777924563d31813103e4b9df55bb081885a8b65b7c6868d13d162f892c72d84"
JSON_VERSION = "v3.12.0"
JSON_SHA256 = "aaf127c04cb31c406e5b04a63f1ae89369fccde6d8fa7cdda1ed4f32dfc5de63"


class NativeBuild:
    def __init__(self, output):
        self.output = Path(output).resolve()
        self.output.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def fetch(url, expected):
        # Some Nix Python builds do not discover the system's CA bundle.
        system_ca = Path("/etc/ssl/certs/ca-certificates.crt")
        context = ssl.create_default_context(
            cafile=os.environ.get("SSL_CERT_FILE")
            or (str(system_ca) if system_ca.exists() else None)
        )
        with urllib.request.urlopen(url, context=context, timeout=60) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != expected:
            raise RuntimeError(f"SHA256 mismatch for {url}")
        return data

    def run(self, jobs):
        source = self.output / f"ogdf-{OGDF_COMMIT}"
        if not source.exists():
            archive = self.fetch(
                f"https://codeload.github.com/ogdf/ogdf/tar.gz/{OGDF_COMMIT}",
                OGDF_SHA256,
            )
            with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as handle:
                handle.extractall(self.output, filter="data")
        json_include = self.output / "json-include"
        header = json_include / "nlohmann/json.hpp"
        if (
            not header.exists()
            or hashlib.sha256(header.read_bytes()).hexdigest() != JSON_SHA256
        ):
            header.parent.mkdir(parents=True, exist_ok=True)
            header.write_bytes(
                self.fetch(
                    f"https://raw.githubusercontent.com/nlohmann/json/{JSON_VERSION}/"
                    "single_include/nlohmann/json.hpp",
                    JSON_SHA256,
                )
            )
        cmake = os.environ.get("CMAKE") or shutil.which("cmake")
        if not cmake:
            raise RuntimeError(
                "CMake is required; use a development shell or set CMAKE"
            )
        subprocess.run(
            [
                cmake,
                "-S",
                str(Path(__file__).parent),
                "-B",
                str(self.output / "build"),
                "-DCMAKE_BUILD_TYPE=Release",
                f"-DOGDF_SOURCE_DIR={source}",
                f"-DJSON_INCLUDE_DIR={json_include}",
            ],
            check=True,
        )
        subprocess.run(
            [
                cmake,
                "--build",
                str(self.output / "build"),
                "--target",
                "ec-spqr",
                "-j",
                str(jobs),
            ],
            check=True,
        )
        binary = self.output / "ec-spqr"
        shutil.copy2(self.output / "build/ec-spqr", binary)
        print(binary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[4] / "target/ec-planarity-native",
    )
    parser.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    NativeBuild(args.output).run(args.jobs)
