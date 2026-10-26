= Five-loop producer source supplement

This source-only supplement preserves the actual private helper that produced
the qualified optimized executable with SHA256
#raw("48917fd9be0abc900f8a511057b5eb0e3e6f3f275c0c7b74cf318ec47cd7264a").
It contains all nine frozen helper source/document files, its resolved Cargo
lock, the actual optimized build launcher/compiler recorder, and the private
Python admission/guard/controller sources and historical configuration examples.
Qualification JSON is evidence, not a compiled binary. No cache, numerical data,
license configuration, credential, private environment, target directory or
compiler-output log is included. The historical #raw("a906414e…") provenance
archive remains separate and unchanged.

== Dependencies and build

Check out #link("https://github.com/alphal00p/rustred")[RustRed] at
#raw("b209adb1e774242bcbe0417fed223317566561bd"), with its Symbolica submodule at
#raw("9cd57678d8817ec33c05cacba4925162f1c86b0b") from the official
#link("https://github.com/symbolica-dev/symbolica")[Symbolica repository].
Graphica and Numerica are trees inside that Symbolica commit; #raw("rustred-order")
comes from the RustRed checkout. #raw("dependencies/public-sources.json") records
these facts, and #raw("dependencies/locked-dependencies.json") inventories all
143 resolved packages and their registry checksums. Registry downloads, Rust,
GMP/MPFR and any selected dynamic libraries remain external dependencies.

The preserved #raw("Cargo.toml") uses absolute paths into the original RustRed
checkout. #raw("build.rs") checks that checkout's clean Git revision and derives
the actual helper/Native/kernel identity. Either restore the original layout
recorded in the manifest, or explicitly relocate those Cargo paths, the
#raw("NATIVE") constant in #raw("build.rs"), and the output-parent
#raw("WORKSPACE") constant in #raw("src/main.rs"). Relocation changes source identity:
record a *new cohort and executable hash*. Never relabel a rebuilt executable
with the original hash or reuse its compiler attestations.

The recorded build used Rust 1.98.1, release optimization 3, fat LTO and one
codegen unit. The included #raw("run.sh") and #raw("rustc_record.py") preserve its
exact flags and source checks. Their original workspace paths and private
#raw("env.sh") reference must be replaced with the reproducer's explicit local
environment. Supply a valid Symbolica license externally; no value is included.
Use a fresh build-output directory and the archived Cargo lock. Run
#raw("build-info") afterward and authenticate the new binary and source cohort.
No fresh-host rebuild or bit-identical rebuild has been tested by this packaging.

== Runtime and provenance

The helper requires its reviewed before-Native admission tuple before reading
five-loop inputs. The archived controller and Python leaf/policy/guard sources
provide the actual implementation: one worker, an exact observed PID/start/unit
identity, a post-spawn ticket and a live PIDFD guard at the 50 GB host-memory
floor. They use Python's standard library and Linux/systemd user services.
The original 1200-second/16 GiB examples are historical pilot envelopes, not
internal reduction limits or ready-to-use permission for another run.

After relocation, author fresh controls with the actual UID, paths, CPU/unit
names, source hashes, executable and loader environment. Preserve all guard
checks; observe process identities instead of copying historical PIDs. The
original controllers contain fixed path/hash bindings, so unchanged examples
are not portable launch scripts. No service or Native operation is started by
extracting or inspecting this archive.

Supply separately the authenticated application, binding, closed Native source
and exact source inventory, plus an admitted usage profile and selected keys.
The frozen schema sources are included. #raw("export-selected") uses the existing
one-gzip/shared-table Native codec and validates exact records and ordered guards
after reimport. Public source admission lives in external-data commit
#raw("cf87c035a2d3a7b3a83126f67c00c0d034944665"); finite-source metadata composition
lives in GammaLoop commit #raw("ab57b4d0ba069295d187793126809024744b098e").
These public workflows are not linked into the private Native-only helper.

Check #raw("SHA256SUMS") and #raw("manifest.json") after extraction. Archive paths,
order, timestamps, ownership and file modes are normalized. This supplement
claims source preservation only: no new build, Native import/export, reduction,
final cache selection, release publication or final-release inclusion.
