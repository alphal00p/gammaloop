= External reduction data and master capabilities

Vakint ships its 1L–3L RustRed programs, master coefficients and graphs without
prepared caches. The separate #link("https://github.com/alphal00p/vakint-external-data")[external-data repository]
contains independently packaged loop data and optional complete native caches.
FORM reduction rules and their hard-coded master formulas remain in Vakint.
See #link("external_data_validation.typ")[the migration validation report] for
release identities, numerical evidence, completed checks and remaining gates.

== Explicit loading

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
each pack was local, reused, or downloaded.

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
Checksums, relative paths, schemas and native program/cache bindings are checked
before activation. Loading a new snapshot replaces the prior external selection;
old clones keep their existing immutable data and cache sessions.

`load_external_data(path, maximum)` loads an existing checkout directly. The
`load_rustred_artifacts(loop_count, path)` method accepts a local pack/IBP path
or a sealed 1L–3L rules file; bare programs retain the active master values and
discard the preceding program's cache attachment. The
explicit native cache load/import/download methods and the saved-application
path interface remain available. No ordinary evaluation initiates a download.

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
