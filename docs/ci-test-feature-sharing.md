# Shared features for test builds

The integration and Python API groups build the same integration-test crate.
`python-api-tests` has no dependencies; it only enables the `test_python_api`
executable. Both groups now use that feature through `testFeatures` in
`nix/ci.nix`. Their archive filters select disjoint binaries, so Python tests
still run once and both groups report separately.

This removes the second integration dependency build, test-binary build, and
their two input merges. Public output names remain available; the Python
context aliases resolve to the shared producers. Regenerate `nix-ci.nix` with
`just ci-update` after changing the workspace or CI definitions.

The feature is applied only to package test contexts. Production artifacts,
the Python module, other test groups, Clippy, and doctests keep their existing
derivation identities. Avoid a workspace-wide `--all-features`: Python linking
and ABI options are distinct configurations, and `no_pyo3` conflicts with Python
support. The remaining explicit feature differences in Linnet (`rkyv`) and
Spenso (`python`) need their own coverage and completion-time measurements.

## Local validation, 9 September 2026

| Check | Main layout | FeynKit layout |
|---|---:|---:|
| Selected integration tests, before and after | 104 | 106 |
| Selected Python API tests, before and after | 5 | 5 |
| Both runtime groups | Passed | Passed |
| Prepare archives and graph/config checks from cached producers | 12.68s | 10.53s |
| Cargo compilation during that preparation or runtime | 0 | 0 |
| Exact preparation repeat | 0.22s, no builders | 0.22s, no builders |
| Explicit NixCI targets, before → after | 55 → 53 | 75 → 73 |

Inventories match by binary, test name, and ignored status; no extra runtime
filter was supplied to make them match. Artifact identities were compared,
including Clippy, doctest, and formatting derivations. Production
crate identities and FeynKit's eighth group remain unchanged. No tests were
edited or retried. Runtime checks used the authorized local license privately.

Identity comparison uses main CI revision `e40357da7182061b3045ca6ee9062c76dc972b52`
and FeynKit revision `19a6460c658f411e7362060603c87d1d1e990af4`, each with this
patch. Main runtime validation applies the same patch to the frozen measured
revision `7e6dc57327fbd6ab3e6181151d66963e8e93675c` to reuse its existing cached
producers; it is not a full build of the newer main revision. FeynKit runtime
validation uses `19a6460…` plus the patch. Disposable stores used eight assigned
CPUs, eight build cores and four concurrent builders, without cache uploads.

The four eliminated steps account for 475s of combined job time and seven Cargo
compilation messages in the previous FeynKit preparation trace. Their durations
overlap: this is removed work, **not a measured 475s completion-time speedup**.
The archive timings above are warm preparation, not clean compilation results.

The shared context currently depends on the Python module even for integration
tests. If that module is rebuilding, integration can start later. Before rollout,
compare individual-group and all-test completion on NixCI, particularly after
Python-only edits. No NixCI suite has been submitted for this feature change;
the earlier CI speedup measurements do not measure it. This also does not
provide native Cargo's persistent incremental compilation within a changed crate.
