= Post-consolidation V2 suite source controls

#strong[DRAFT / INACTIVE: source preservation only.]
This capsule preserves the later V2 input identities, inactive qualification
controls, source-exclusion inventory, and historical collector normalization.
It does not establish a final runtime freeze, authorize execution, or supply
new qualification, mapper, cache-selection, or performance results.

The archive contains 28 explicitly whitelisted source/control files: 21 JSON,
four Python, one shell, one Typst, and one source patch. The source directory
had 29 files; its generated #raw("readiness.pdf") is excluded while canonical
#raw("readiness.typ") is retained. No external payload referenced by a control
is copied into the capsule. The earlier proposal capsule from commit
#raw("7decbbab8d7b8ca59417c770deceb5030d2404eb") remains unchanged.

The original collector cohort is #raw("03ad1ecf"), RustRed #raw("4190f2c6"),
and Symbolica #raw("9cd57678"). Its collector and mapper ELF hashes begin
#raw("319ddef4") and #raw("6e5e1b8d"). These historical identities must not be
relabeled as the new consolidated source. Workspace paths and missing receipt
references are preserved verbatim; this is not a portable executable package.

Before any future activation, seal the actual consolidated build and cohort,
materialize new controls, obtain a fresh coordinator GO with the existing
PIDFD and host-memory protections, and produce genuine complete uncached
qualification and mapper receipts. Keep known-seen historical controls outside
eligible new tuning and held-out arrays. No such jobs ran to create this archive.

#raw("suite-v2-post-consolidation-source-20261026.tar.gz") is 79,340 bytes;
its SHA-256 is
#raw("0a34c0f84a08661a9da95e3c648fa0e29e86235e6ee31aa155646b838601bb5e").
#raw("source-supplement-manifest.json") records every member's length and
SHA-256, exclusions, cohort boundaries, and validation. Archive paths are safe
regular files; order, modes, ownership, timestamps, and gzip headers are
normalized. Two encodings matched, and all extracted bytes matched the source.
The checksum list covers the archive and each extracted member.
