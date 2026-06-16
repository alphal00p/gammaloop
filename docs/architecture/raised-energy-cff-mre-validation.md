# Standalone Symbolica MRE evidence

Executed revision: `a3af927afc633c6e6c2474d3f03fa02e0242eaa4`. The standalone package passes locked all-target check/build/Clippy and two controls whose assertions require **39** and **4**. The tiny writer exits 0; the reader exits the expected **101**, reproducing `f(y,x)` where `f(x,y)` is required. This is **two positive controls and one reproduced defect**, contributing **zero workspace tests**.

Observed C7: `6d7a9220129c0fac37bdc0e04ebefccb1688c7be`. Source relationship: **Source-equivalent standalone package; no additional execution**. Package tree `b917afeca167a3b114c6f8bc0d22be24dc888292` and the exact manifest/lock/Cargo-config/toolchain-file comparison are recorded in the JSON sidecar. The executed revision and original receipt/binary hashes are retained.

| Stage | Classification | Exit | Seconds |
| --- | --- | --- | --- |
| check | build_gate | 0 | 42.099 |
| build | build_gate | 0 | 315.154 |
| clippy | build_gate | 0 | 0.978 |
| control | positive_control | 0 | 0.498 |
| abstract-control | positive_control | 0 | 0.501 |
| import-write | diagnostic_input_creation | 0 | 0.495 |
| import-read | reproduced_upstream_import_diagnostic | 101 | 0.493 |

Separately executed standalone formatting: **PASS** at `bbb8ff1db499a78b3b25c41f774c9a034acb5b4d` in 0.479s, with complete tracked-source guards. Package/build-input identity with observed C7: **True**. This is a formatting operation, not a control or additional runtime execution.


No GL262 physical replay or large evaluator-construction panic reproduction was attempted. The import argument-order defect does not establish the cause of the GL262 panic. Controls build small synthetic inputs; they do not decode or validate the three serialized Symbolica blobs.

The original unpinned formatting receipt has no revision, time or source guard and is not a fresh pinned MRE gate. The separate pinned formatting record, when present, carries its own revision and source guards.

Exact result/output fingerprints and complete tracked-source guards are verified against the original immutable revision. Binary and tiny-input hashes are retained from the execution driver; this summarizer does not read/hash the binaries or read/decode the serialized input. The original toolchain/environment and execution date remain authoritative.
