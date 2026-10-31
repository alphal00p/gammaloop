= Exact RustRed private host prerequisite proposal

The missing `IBPFamily` belongs to RustRed's native FeynKit bridge. The f47
minimal notebook host registers FeynKit alone. The committed full Community
host calls both FeynKit initialization and RustRed registration on its HEP
module. The first difference is native host composition before the attribute
lookup, not demonstrated failure of a shared family algorithm. Keep the usual
notebook host minimal; use the private validation superset supplied here.

The root-selected candidate is `alphal00p/rustred` commit
`7c1ed03722b8c05daf60c89ba4ecc79457ed2ada`, tree
`c78d15aec5648706f066e939650aecdced550bbd`. Despite older checkout names, exact
local Cargo Git objects contained it. All 2,955 archived file blobs were checked
against the commit and the immutable store snapshot. The source is
`/nix/store/5hr22rycs6dpj77afgq0maqd8qbipk6a-rustred-source-7c1ed037`.
No network retrieval or donor-tree modification occurred. The unmaterialized
vendored Symbolica gitlink is not selected. Community's own FeynKit and kernel
pins were not adopted: actual inputs remain f47 and 8eb.

`proposal/private-host.patch` supplies exactly three private host changes:
manifest, standalone lock and registration. The existing registration macro
returns its child module; `rustred_feynkit::register_hep_module` adds the genuine
IBP classes to that same module. The original public wrapper remains unchanged.
RustRed source is unchanged and loaded from its exact content-addressed source.
The whole host is one extension with one Symbolica/PyO3 graph.

The bridge has no direct Cargo dependency on another FeynKit revision. It reads
the existing independent/complete flags, momenta, denominators, kinematics,
dimension and scalar-product API, extracting native PythonExpression values.
Those properties exist at f47. The inspected 8eb atom and PyO3 boundary APIs
also exist. This source inspection is not a whole-crate type-check certificate.
No concrete missing method, type or links conflict was found at the inspected
boundaries. Compilation and scientific output remain unqualified.

Disable RustRed defaults and campaign-api. Only rustred-order, rustred and the
bridge are selected; reconstruction and the Python/application adapters remain
disabled. Their optional lock entries are retained. The upstream build selector
`RUSTRED_RUNTIME_ARITIES=1` supplies the one-denominator QCD family. This does
not qualify multi-denominator reduction or complete Community coverage.

The candidate lock retains all 464 existing package records except the host's
new dependency edge and adds five RustRed path records. Registry versions and
8eb kernel sources do not change. Static references are unambiguous and closed.
No Cargo resolver was run. The future producer first runs locked offline
metadata and verifies the unique kernel, PyO3 versions, source paths and empty
RustRed feature lists before compilation. No root lock is edited. A production
adoption would require the ordinary coordinated lock/subset and CI workflow.

`proposal/private-producer.nix` reuses the qualified bounded wheel pattern,
FeynKit artifact owner and existing immutable vendor set. Enforced pure offline
Nix evaluation, with import-from-derivation disabled, succeeded from a sealed
content-addressed proposal. It binds both exact source snapshots. The historical
graph has 2,996 derivations; pruning already available outputs leaves 28 reused
frontier outputs and one unbuilt private wheel output. This is derivation and
local-store evidence, not compilation success. Cargo may rebuild additional
units when final feature fingerprints require it.

`next-build.sh` contains the exact proposed build command. It remains unexecuted
and requires root authorization for the next native build phase. The producer retains metadata, host locks,
arity selection, fingerprints and depfiles for the later provenance gate.
`next-commands.json` specifies the pure-to-runtime boundary: inspect actual new
wheel/NAR/DSO hashes and ELF dependencies, prepare a fresh private environment,
install the entire wheel offline, seal a compiler-free namespace, and only then
enter separately authorized licensed product execution. Old native hashes and
old passing tests cannot qualify this new host.

Run the original f47 installed QCD script at its original full source/fixture
layout. Its SHA256 remains
`2622f4f8f53624e177e6ae059eb12953f1b16f6cf08ebe85b70ff6155423dfc6`.
All 23 assertion sites, targets `[[2], [3]]`, depth 2, residual assertion `[[1]]`,
scientific reference inputs and original completion marker remain unchanged.
The stub-only correction is a separate packet. No alias, analytical replacement,
skip, native pass, new-host credit or review approval is claimed.

All new workspace files stay inside this ignored packet. No tracked file,
reference, index, JJ state, dirty donor manifest, shared installation, license,
credential or user configuration was changed or read for credential material.
No owned process remains. Exact snapshots, source errors and sealed evidence
are retained for root review and continuation.
