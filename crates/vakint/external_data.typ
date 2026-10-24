= External reduction data and master capabilities

Vakint ships its 1L–3L RustRed programs, master coefficients and graphs without
prepared caches. The separate #link("https://github.com/alphal00p/vakint-external-data")[external-data repository]
contains independently packaged loop data and optional complete native caches.
FORM reduction rules and their hard-coded master formulas remain in Vakint.
See #link("external_data_validation.typ")[the migration validation report] for
release identities, numerical evidence, completed checks and remaining gates.

== Explicit loading

The official-download API and options below are available on native Unix targets.
They are not exposed on Windows or WebAssembly targets.

```rust
use vakint::{OfficialExternalDataOptions, Vakint};
let mut engine = Vakint::new()?;
let report = engine.load_official_external_data(4, &OfficialExternalDataOptions {
    local_artifacts_path: Some("./data".into()),
    ..Default::default()
})?;
```

The matching Python method is:

```python
report = engine.load_official_external_data(
    4, version="latest", local_artifacts_path="./data", refresh=False
)
```

The loop argument is always a maximum: four loads 1L, 2L, 3L and 4L. All packs
come from one immutable release. The default path is the current working
directory, with managed installations below `.vakint-external-data/<version>`.
The report identifies the version, per-loop paths, available caches, and whether
each pack was local, reused, or downloaded. Prepared caches are loaded by default:
`OfficialExternalDataOptions.load_prepared_caches` is `true`, and the Python
keyword `load_prepared_caches` defaults to `True`. Set the Rust option to `false`
or pass `load_prepared_caches=False` in Python to load the same authenticated data
without native cache decoding or import. Checksums still cover every listed
file, including optional cache bytes; this is an activation policy, not a change
to the dataset. Use this option with `local_artifacts_path` for a local checkout.

`cache_available` reports the manifest declaration. Per-loop `cache_loaded`
reports actual automatic attachment, and the top-level `load_prepared_caches`
field records the requested policy. Availability may therefore be true while
loading is false. Explicit cache load methods remain available after data loading.

If the supplied directory contains the repository manifest, Vakint validates and
loads it directly. It never downloads into or modifies this checkout. `latest`
means the checkout's current content; an explicit version must match its declared
identity. A missing file or incorrect checksum is an error.
Loop-pack manifests without a valid root release manifest also identify damaged
local data and cause an error without downloading into that checkout. An
unrelated valid `manifest.json` is preserved and does not block a managed install;
empty loop-named directories alone do not identify an external-data checkout.

Managed installations reuse validated files without downloading or extracting
again. `refresh=True` checks the `data-latest` release alias; pinned versions remain
immutable. A failed refresh reports an error and leaves the active engine intact.
Checksums, relative paths and schemas are checked before activation. Native
program/cache bindings are checked when their native payloads are activated. Loading a new snapshot replaces the prior external selection;
old clones keep their existing immutable data and cache sessions.

`load_external_data(path, maximum)` loads an existing checkout directly. The
`load_rustred_artifacts(loop_count, path)` method accepts a local pack/IBP path
or a sealed 1L–3L rules file; bare programs retain the active master values and
discard the preceding program's cache attachment. The
explicit native cache load/import/download methods and the saved-application
path interface remain available. No ordinary evaluation initiates a download.

== Four-loop acceptance tests

The four-loop numerical acceptance gates and native asset tests call the official
loader with maximum loop count four, the compatible immutable release declared
by the vendored manifest, and `refresh=false`. This avoids a stale cached `latest`
marker while preserving offline reuse.
Tests that inspect the five-loop application or graphs explicitly request five.
On first use it downloads
missing packs into the test process's current directory; subsequent processes
reuse the validated installation without network access. Cargo normally runs a
test binary from its package directory. Set `VAKINT_TEST_EXTERNAL_DATA` to an
existing checkout to override that location.
Managed `.vakint-external-data` directories are ignored by Git.

The scientific gates retain their explicit execution flags and FORM requirements.
For example, this existing gate uses the RustRed backend and an archived
independent reference without invoking FORM:

```sh
cargo test -p vakint --test rustred_four_loop_tests \
  h_rank_four_matches_reference_without_form -- --ignored --nocapture
```

Ordinary tests using only the vendored snapshot continue to run offline.

== Master representation and numerical capability

Rust uses `RustRedEvaluationOptions { master_mode: RustRedMasterMode::Symbolic }`.
Python uses `VakintEvaluationMethod.new_rustred_method(master_mode="symbolic")`.
The other values are `retain` and `numerical`; numerical substitution is default.

Retain mode returns family-qualified master integrals. Symbolic mode expands
available exact Laurent coefficients and leaves numerical-only coefficients as
snapshot-qualified symbols. Evaluate those through the same engine's numerical
evaluation method so it can resolve their data. Switching to a different data
snapshot does not silently reinterpret old coefficient symbols.

Numerical mode evaluates exact formulas at the requested precision. Numerical-only
coefficients enforce their validated significant-digit limit. Required epsilon
powers are computed from actual reduced coefficients, including reduction poles,
products and normalization, after exact cancellation. Missing entries are not zero.
Errors name the integral, coefficient power, requested capability, available data,
and package version. This source-data check is not a final-result error estimate.

HPLs at `exp(±i*pi/3)` use Symbolica's arbitrary-precision numerical algorithms.
The same evaluator supplies FORM-backed HPL substitution. Finite constants without
an evaluatable exact definition retain explicit numerical precision; resizing a
stored decimal does not provide more accurate digits.

== Partial campaigns

Recorded coverage predicates describe the actual campaign scope. Summary maxima
and completed cache workloads do not enlarge that scope. Outside certified
coverage, Vakint warns and applies only rules whose mathematical conditions hold.
Unsupported integrals remain exact, family-qualified residuals. Corruption, unsafe
guards, invalid transports and operational-limit failures remain errors.

After cancellation, warnings list remaining undeclared masters and parent-topology
DOT paths. Numerical evaluation reports these masters again, together with declared
terminals lacking coefficients. The initial 5L release has 905 declared campaign
terminals and no computed 5L coefficients; it does not certify a complete basis.

== Maintaining the vendored snapshot

Run `python3 tools/vendor_external_data.py /path/to/vakint-external-data` from the
Vakint crate. The tool validates input hashes and copies all available 1L–3L
coefficients, graphs and programs, excluding caches. Native cache compatibility is
determined by codec and authenticated program identities, not producer revision
strings alone. Rebuild only caches invalidated by actual program/format changes.

Mathematica may be used to develop independent numerical reference values. It is
not a build or runtime dependency of Vakint or the official-data loader.
