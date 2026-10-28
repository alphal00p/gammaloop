= Staged signed sampling reconciliation

This preparation ports the cumulative approved sampling production and regressions
at `c82916079b5cc643151e0213499ca9875739a4e3` onto reviewed staged sources 22–26 at
`f2e32c44c986058ede1d7e585b9e2e4fefe1110c`. The original sampling bookmark and
independent production and P3 reviews remain preserved. This branch omits sources
19–21. Ordered main replay must apply 19–21, then 22–26, then 27–29 and rerun
combined acceptance. A NEW independent Astra review is required before that replay.

== Sampling and physical ownership

Source 27 co-orients parallel bundles against their first stored edge. Selected
three-momentum input blocks are reflected before the existing affine routing:
`A R k + b`. Fixed external shifts are unchanged. The parent channel uses an
identity frame when reflection is required. Forward mapping, inverse density,
absolute Jacobian and OSE weights share the same frame. Evaluation still receives
signed master momenta; physical edge incidence and numerator ownership remain
with the current graph. Denominator evenness never supplies numerator signs.

Source 28 omits only exact numeric-domain zero affine coefficients. Small nonzero
coefficients remain at active precision; no approximate threshold or conversion
through binary64 is introduced. Existing finite-input, NaN and underflow policies
remain in their current owners. Source 29 retains the specifically accepted CoolQM
three-sigma gate, seed 1337, 750000 samples, original target and error semantics.
Other tolerances and random streams are unchanged.

== Staged interfaces and graph regression

The patch is applied as cumulative hunks, without replacing old files wholesale.
Graph tests use the current native compact DOT parser with scalar particles;
obsolete half-edge `ufo_order` metadata is removed. The existing indexed numerator
interface exports the authored edge factors. Current NodeForest incidence, signed
routing, graph finalization, hidden-history admission and source-error owners are
retained. No compatibility aliases or older storage are restored.

The graph-derived mixed two-loop regression retains three graph-generated bases,
both map-density and OSE weighting, two cube points per selected channel and each
parent basis: 36 generating cases and 108 inverse checks. Its independent routing
boundary uses the physical edge signatures, including fixed external shifts. The
odd factor belongs only to the reversed physical edge, and factors stay separate.
Channel identities supply score lookup; existing catalogue and deduplication
controls supply channel coverage. The weighted numerator sum tests partition
unity, not independent full-evaluator multiplicity.

== Evidence and limits

#link("evidence/generation-sampling-staged/receipt.json")[The receipt] records exact
source mappings, staged code identity, files, commands, producer identities and
fresh results. Old 74ce-based physics evidence remains historical and supplies no
new staged acceptance credit. Logs and review inputs remain in the task workspace.

Validation uses the Nix development environment, silently preserved or loaded
license, two Cargo jobs and the configured shared cache without copying, manually
seeding, deleting or overriding it. NixCI remains disabled. The protected kernel,
root manifests and locks, Tydenso, uv and notebook host are unchanged. Inherited
source31 documentation scope failures are retained; this lane does not refresh
unrelated digests or duplicate the staged22 audit correction.

BRST, full WP3, historical multi-seed studies and full unrelated heavy physics are
outside this preparation. Stage19 grouping is absent and final combined acceptance
is mandatory. No main integration, publication, upload, push or PR is authorized
by this result.
