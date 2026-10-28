"""Typed cache identities on one coherently built standalone or Community host.

These schema and forwarding checks need no external campaign or cache payload.
Set TMPDIR to an owned workspace directory. The generated-file integration
control lives in running_cache_smoke.py and authenticates actual native output.
"""

import os
import tempfile
import unittest
from pathlib import Path

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    from vakint_hepkit_standalone import CacheFileIdentity, Vakint
else:
    from symbolica.community.hepkit.vakint import CacheFileIdentity, Vakint


def document():
    return {
        "compressed_bytes": (1 << 63) + 7,
        "compressed_blake3": list(range(32)),
        "decoded_bytes": (1 << 64) - 1,
        "decoded_blake3": list(reversed(range(32))),
    }


class CacheFileIdentitySmoke(unittest.TestCase):
    def test_native_document_round_trip_and_keyword_constructor(self):
        expected = document()
        for identity in (
            CacheFileIdentity.from_document(expected),
            CacheFileIdentity(**expected),
        ):
            with self.subTest(identity=identity):
                self.assertEqual(identity.to_document(), expected)
                self.assertEqual(
                    identity.compressed_bytes, expected["compressed_bytes"]
                )
                self.assertEqual(identity.decoded_bytes, expected["decoded_bytes"])
                self.assertIsInstance(identity.compressed_blake3, bytes)
                self.assertIsInstance(identity.decoded_blake3, bytes)
                self.assertEqual(
                    identity.compressed_blake3, bytes(expected["compressed_blake3"])
                )
                self.assertEqual(
                    identity.decoded_blake3, bytes(expected["decoded_blake3"])
                )

    def test_identity_and_exported_document_are_independent(self):
        original = document()
        identity = CacheFileIdentity.from_document(original)
        exported = identity.to_document()
        original["compressed_blake3"][0] = 255
        exported["decoded_blake3"][0] = 255
        exported["compressed_bytes"] = 0
        self.assertEqual(identity.to_document(), document())
        for field in document():
            with self.subTest(field=field), self.assertRaises(AttributeError):
                setattr(identity, field, None)

    def test_native_schema_rejects_malformed_documents(self):
        cases = []
        for field in document():
            missing = document()
            del missing[field]
            cases.append(missing)
        extra = document()
        extra["unexpected"] = 1
        cases.append(extra)
        for field in ("compressed_blake3", "decoded_blake3"):
            for value in (
                [],
                [0] * 31,
                [0] * 33,
                [-1] + [0] * 31,
                [256] + [0] * 31,
                [True] + [0] * 31,
                [1.5] + [0] * 31,
                "00" * 32,
                None,
            ):
                invalid = document()
                invalid[field] = value
                cases.append(invalid)
        for field in ("compressed_bytes", "decoded_bytes"):
            for value in (-1, 1 << 64, True, 1.5, "1", None):
                invalid = document()
                invalid[field] = value
                cases.append(invalid)
        for invalid in cases:
            with self.subTest(document=invalid), self.assertRaises(ValueError):
                CacheFileIdentity.from_document(invalid)
        for invalid in (None, [], "{}"):
            with self.subTest(document=invalid), self.assertRaises(TypeError):
                CacheFileIdentity.from_document(invalid)

    def test_keyword_constructor_rejects_invalid_digests_and_counts(self):
        for field in ("compressed_blake3", "decoded_blake3"):
            for length in (31, 33):
                invalid = document()
                invalid[field] = [0] * length
                with (
                    self.subTest(field=field, length=length),
                    self.assertRaises(ValueError),
                ):
                    CacheFileIdentity(**invalid)
        for field in ("compressed_bytes", "decoded_bytes"):
            for value in (-1, 1 << 64):
                invalid = document()
                invalid[field] = value
                with (
                    self.subTest(field=field, value=value),
                    self.assertRaises(OverflowError),
                ):
                    CacheFileIdentity(**invalid)

    def test_load_cache_forwards_typed_identity_before_lower_loop_path_access(self):
        engine = Vakint()
        identity = CacheFileIdentity.from_document(document())
        with tempfile.TemporaryDirectory(prefix="vakint-file-identity-") as scratch:
            absent = Path(scratch) / "must-not-be-read-or-created"
            with self.assertRaisesRegex(ValueError, "single-file identity/options"):
                engine.load_cache(1, str(absent), file_identity=identity)
            self.assertFalse(absent.exists())
            with self.assertRaises(TypeError):
                engine.load_cache(1, str(absent), file_identity=document())
            self.assertFalse(absent.exists())


if __name__ == "__main__":
    unittest.main()
