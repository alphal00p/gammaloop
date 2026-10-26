= Private source-alias helper checkpoint

This auxiliary source-only archive preserves the reviewed private route
qualifier and its unchanged resource-admission module, build launcher, and
three small source/linkage/readiness records. It is a diagnostic checkpoint,
not a production API, registered Cargo example, benchmark result, or cache.
The archive is
#raw("vakint-5l-source-alias-helper-20261026.tar.gz");
#raw("source_alias_checkpoint.json") records its checksum and validation.

== Source identity and linked libraries

The helper's own frozen source identity is recorded by
#raw("helper-first3/source-identity.json"). Its qualifier source has SHA256
#raw("8c5ed3e5447a5c4109acb4c9ada2f4e6263afbdbd22bcc013203c9bdd59418d8");
the unchanged admission module has SHA256
#raw("bb3223d27584b42069d3813ebdcaefa6598a0935ca4531714a50248b334870cb").
The helper has no assigned Main revision. This publication commit records
source preservation; it does not retroactively identify the helper as code
compiled from that commit.

The linked libraries separately belong to the authenticated consumer cohort:
Main #raw("03ad1ecfd75b776723562ab2a17e1393b17ccd84"),
RustRed #raw("4190f2c6f3158e3fc3cbcb1937e12e5971da9175"), and official
Symbolica #raw("9cd57678d8817ec33c05cacba4925162f1c86b0b").
The cohort receipt's SHA256 is
#raw("ecde038c079ab65012d02c15bf598f16246277eaf614233a9c9b39b9037ce5ca").
The archived linkage audit inventories exact compiler, library and native-link
inputs; none of their binaries or build directories is included. Paths and
hashes in these frozen records describe the original private layout.

== State and validation limits

The source/build snapshot sealed at 2026-10-26 17:44:40 UTC records
#raw("SEALED_NOT_EXECUTED"), with zero Native executions at that point and no
recorded input drift. The sealed executable's recorded SHA256 was
#raw("43861e6e373128c5dea9445ebd4faf732969f568cdfb7722b10b564428f3a12d").
Runtime admission, source-owner qualification, and protected-alias exclusion
remain #raw("PENDING") in this source/build snapshot. No later runtime controls,
proofs or acceptance are included by this backup. The earlier readiness note
describes the inspection before this compilation; it must not be read as a
later runtime result.
No executable, RLIB, scientific proof, raw trace, runtime output, cache payload,
private environment, or license configuration is included.

The checkpoint checks exact original bytes and per-member SHA256 values, two
byte-identical archive passes, normalized regular-file members, safe paths,
secret scans, and Rust formatting. These checks do not constitute a fresh
compile, Native execution, reduction, complete CI run, runtime qualification,
cache-independence proof, or final-release acceptance. Existing canonical
benchmark results and dependency pins are unchanged.

== Recovery and relocation

Extract only into a new private workspace and verify #raw("SHA256SUMS") and
#raw("manifest.json"). Extraction starts no process. The archived
#raw("build.py") is the exact original launcher, not a portable installer;
its absolute paths and evidence bindings must not be silently reused elsewhere.

Any relocation, source edit, rebuild, library change, or new registration needs
fresh source hashes, an authenticated library/cohort binding, and fresh explicit
registration and resource admission before execution. Recheck the compiler and
link inputs, seal the new executable, and preserve its new identity. Do not
reuse historical PIDs, guard tickets, executable hashes or compiler attestations.
Keep the helper's source identity distinct from its linked library cohort;
publishing or extracting this archive authorizes no Native operation.
