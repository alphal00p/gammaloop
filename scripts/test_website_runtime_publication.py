"""Shared-file publication must verify bytes and retain historical identities."""

import gzip
import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path


class RuntimePublicationTest(unittest.TestCase):
    def test_integrity_collisions_and_historical_retention(self):
        script = Path(__file__).with_name("update-docs-pages.sh")
        functions = script.read_text().split('[ "$#" -ge 3 ]', 1)[0]
        command = (
            functions
            + "\nbuild_root=$1\npages_root=$2\ncheck_shared_runtime_assets\npublish_shared_runtime_assets\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build, pages = root / "build", root / "pages"
            payload = b"verified runtime payload"
            digest = hashlib.sha256(payload).hexdigest()
            source = build / "assets/notebook-runtime-files" / digest
            source.mkdir(parents=True)
            (source / "native.wasm").write_bytes(payload)
            (source / "native.wasm.gz").write_bytes(gzip.compress(payload))
            old = pages / "assets/notebook-runtime-files" / ("a" * 64) / "old.zip"
            old.parent.mkdir(parents=True)
            old.write_bytes(b"historical archive")

            def run():
                return subprocess.run(
                    ["bash", "-c", command, str(script), str(build), str(pages)],
                    capture_output=True,
                    check=False,
                )

            self.assertEqual(run().returncode, 0)
            installed = pages / source.relative_to(build) / "native.wasm"
            self.assertEqual(installed.read_bytes(), payload)
            self.assertEqual(old.read_bytes(), b"historical archive")
            (source / "native.wasm.gz").write_bytes(gzip.compress(b"wrong bytes"))
            self.assertNotEqual(run().returncode, 0)
            self.assertEqual(installed.read_bytes(), payload)
            (source / "native.wasm.gz").write_bytes(gzip.compress(payload))
            installed.write_bytes(b"collision")
            self.assertNotEqual(run().returncode, 0)
            self.assertEqual(installed.read_bytes(), b"collision")


if __name__ == "__main__":
    unittest.main()
