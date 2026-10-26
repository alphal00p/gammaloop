= Native cache-union census expectations checkpoint

This new supplement preserves the exact private v4-census helper source and its
reproduction controls. Earlier source supplements and all cache payloads remain
unchanged. The new mode delegates an arbitrary ordered number of committed
sources to RustRed's existing `copy_base_union_to_new`; it adds no codec,
coefficient import, reduction system, dependency pin or cache-selection policy.

Frozen predecessor v4 source and its #raw("ad6f13ee…") executable remain
unchanged. This revision requires generation, stored-record count and index
length, while prior run/chunk-reference/referenced-byte expectations may be
explicit null. Every supplied value must match Native. Full output arithmetic
uses actual authenticated Native censuses exclusively; unknowns are never
reported as zero. Missing fields, mismatches and overflow are tested.

The actual private compiled source identity is
#raw("87a1a0a48c80deaabbe82a31ab9ad59a372f50a7c07d997571f48472fc6225bc").
Its separate local optimized executable is
#raw("ba4a7b9c6b26cbed29ff9a00ef17c1d44424ccf06829467c6897494e142f967d").
Neither identity is the later source-publication commit. This utility links
Native #raw("b209adb1e774242bcbe0417fed223317566561bd") (PR14), official Symbolica
#raw("9cd57678d8817ec33c05cacba4925162f1c86b0b") and its Graphica/Numerica
subtrees, with native state codec 6. It does not link Vakint Main.

The actual isolated release gates passed: 38 synthetic tests, zero failures,
one deliberately ignored real-admission test; Clippy, release build and
#raw("build-info") all exited zero. Own code was compiled at optimization 3,
fat LTO and one codegen unit; matching already-qualified dependency artifacts
were independently warm-copied. The bounded compile unit used CPU149, 16 GiB
maximum RAM, no swap and a one-hour maximum. Its actual wall time was 402.799
seconds. These gates opened no real Native cache or application and performed
no copy, import, reduction or selected export. No fresh-host rebuild, real
N-source union, final HEAVY/LIGHT cache or release acceptance is claimed here.

After extraction, verify #raw("manifest.json") and #raw("SHA256SUMS"). The
archive contains all ten frozen helper source/documentation files, including
#raw("Cargo.lock"), #raw("build.rs"), the shared admission module and new
#raw("src/union.rs"). It also preserves the actual build launcher/recorder,
existing Python admission/guard source, dependency pins, package inventory and
reused deterministic packaging source. Regular files, sorted paths, normalized
tar metadata and gzip timestamp zero were checked; two archive builds match.
No executable, target tree, cache/index/chunk, proof, environment, development
license, credentials or scientific result is included.

#raw("helper/PRODUCER.typ") documents the request, ordering, validation and
Native API boundary. It retains earlier design notes at their original
inspection times; later acceptance is limited to the four gates stated above.
The original absolute paths in #raw("Cargo.toml"), #raw("build.rs") and controls
are historical build facts. To relocate, check out the exact public Native and
Symbolica revisions, adjust all source roots deliberately, preserve lock and
features, establish a new source freeze, and run the documented gate commands
with one permitted CPU, #raw("RUSTRED_RUNTIME_ARITIES=15") and the stated
release profile. A changed path/source creates a new cohort and executable;
never reuse the old executable or attestation identity. Registry dependencies
are recorded with locked checksums and must be supplied normally.

The archived launcher is an original control example, not a portable installer:
its private #raw("env.sh") and former proof-root inputs are intentionally absent.
Recreate environment/toolchain and a truthful local source freeze before using
it. The compile gates deliberately unset all five runtime-admission variables.
Real data modes require the unchanged shared before-Native admission, an actual
registered per-job guard, a fresh memory-floor check and explicit closed-source
runtime controls. Publishing or extracting this checkpoint supplies none of
that authorization. No live q0 source or cache was touched.
