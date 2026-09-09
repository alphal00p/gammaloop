# Python inputs for CI tests

The shared integration/Python Rust binaries do not import the Python extension
while compiling. Supplying the extension in their build-time `PYTHONPATH` added
a dependency on its entire production build. The ordinary integration group
then waited for that path even though its filter excludes the Python API tests.

Compilation and archive creation now use the same Python interpreter and Python
support packages without the extension. The Python API runtime group explicitly
requests the extension. Both groups still share one test-binary producer, with
unchanged features, profile and disjoint filters. No project test was modified.

## Local validation, 9 September 2026

On FeynKit, preparing the two archives with cached upstream dependencies took
254.64s, with seven workspace Cargo compilation messages and no external ones.
The Python extension was physically absent from the disposable store throughout
that build. An unchanged repeat took 0.18s without builders.

After copying the unchanged extension for runtime use, both groups passed in
151.32s without compilation: 106 integration tests passed (67 skipped), and all
five Python API tests passed. The archive inventories match the previous
candidate: 173 integration entries and five Python entries, with identical names,
test kinds and ignored flags. This is correctness validation and preparation
measurement, not a demonstrated remote speedup.

Derivation inspection on main and FeynKit confirms that the extension has been
removed from both archive build graphs and the ordinary integration runner, and
retained by the Python API runner. Its underlying derivation is unchanged;
revision-specific scheduling wrappers differ as intended. Clippy, doctests,
formatting and all other test-group derivations are unchanged. Generated CI graph
checks pass on both layouts. New runtime builds were executed on FeynKit only.

The first disposable-store attempt was stopped after 219s: its seed retained
final archives but lacked intermediate dependency outputs, so it began rebuilding
them. The successful attempt uses the complete earlier validation store. Both
attempts and the interruption reason are retained in the measurement directory.

Measurements use eight fixed CPUs (8–15), eight build cores, up to four builders,
a shared host and no uploads. Evidence is under
`/tmp/gammaloop-ci-validation/python-runtime-only/`, including inventories,
`identity-summary.json`, build phases and runtime summaries. The original working
checkout and shared application branches are unchanged.
