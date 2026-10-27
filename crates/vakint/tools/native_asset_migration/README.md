# Native asset codec migration

This offline bridge transports Vakint's 15 trusted RustRed assets from
Symbolica state format 6 to 11. It does not generate rules or recompute masters.
The runtime contains only the migrated assets, not this bridge.

`bridge.rs` is compiled separately against the two Symbolica revisions. The old
decoder exports exact integer coefficients, exponent vectors and complete
ordered variable maps for rational polynomials; general catalog atoms use
canonical text. Both decoders verify their own text/term round trip. A fresh
current decoder must produce identical portable records after encoding. Every
process initializes unrelated symbols before importing the asset.

`migrate.py` preserves the envelope header, dictionary positions and every
non-algebra section byte for byte. It checks the source hashes against the
October 5 receipt and records before/after hashes, sizes, frame counts and exact
portable-record hashes. Gzip output uses `gzip -n -9`.

Use Python 3, Cargo/Rust compatible with the selected Symbolica revisions,
GMP/MPFR build prerequisites, `gzip`, and an available Symbolica license. Obtain
clean source checkouts of these revisions from the Symbolica repository:

- Decoder: `a19c760dd567c239f30d87e4e924ca2f8b8457ab`.
- Encoder: `f514a7fe1b19e851207f9421ef9a4a4157d5f0e5`.

Run from the repository root, providing the original format-6 data directory:

```sh
python3 crates/vakint/tools/native_asset_migration/migrate.py \
  /path/to/symbolica-state6 /path/to/symbolica-state11 \
  /path/to/original/data/rustred /tmp/migrated-rustred /tmp/codec-work
```

The output directory must not exist. The work directory retains both Cargo
manifests/lockfiles, binaries and portable records for inspection. The tool
checks the source Git revisions; use unmodified checkouts. Generated outputs
must also pass Vakint's `rustred_asset_compatibility_tests`, which cold-load the
programs, replay the normalizers and check exact catalog key coverage. These
checks do not establish arbitrary-index candidate closure or numerical accuracy.
