"""Build the native EC layout engine with pinned OGDF; downloaded sources stay under target/."""

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

    def run(self, jobs, wasm=False):
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
        configure = [cmake]
        if wasm:
            emcmake = os.environ.get("EMCMAKE") or shutil.which("emcmake")
            if not emcmake:
                raise RuntimeError("Emscripten is required for --wasm; set EMCMAKE")
            configure.insert(0, emcmake)
            # Toolchain adaptation only; OGDF graph algorithms are unchanged.
            compiler = source / "cmake/compiler-specifics.cmake"
            original = 'if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NOT ${CMAKE_SYSTEM_PROCESSOR} MATCHES "^arm")'
            patched = original.replace(" AND NOT ${", " AND NOT EMSCRIPTEN AND NOT ${")
            text = compiler.read_text()
            if original not in text and patched not in text:
                raise RuntimeError("Pinned OGDF compiler configuration has changed")
            compiler.write_text(text.replace(original, patched))
        build = self.output / ("build-wasm" if wasm else "build")
        target = "ec-layout" if wasm else "ec-spqr"
        subprocess.run(
            configure
            + [
                "-S",
                str(Path(__file__).parent),
                "-B",
                str(build),
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
                str(build),
                "--target",
                target,
                "-j",
                str(jobs),
            ],
            check=True,
        )
        filename = "ec-layout.wasm" if wasm else "ec-spqr"
        binary = self.output / filename
        shutil.copy2(build / filename, binary)
        print(binary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        help="Build cache directory (default: target/ec-planarity-native in the workspace)",
    )
    parser.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument(
        "--wasm", action="store_true", help="Build the standalone Typst plugin"
    )
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    output = args.output
    if output is None:
        output = Path(__file__).resolve().parents[4] / "target/ec-planarity-native"
    NativeBuild(output).run(args.jobs, wasm=args.wasm)
