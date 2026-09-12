# Agent Guidance

Agents and human contributors follow the same repository guidance. Before making
changes, read [CONTRIBUTING.typ](CONTRIBUTING.typ) and treat it as authoritative
for this repository.

If instructions conflict, prefer the more specific local guidance and preserve
the user's current work unless explicitly asked to change it.

Use descriptive branch names without an agent prefix, for example
`ci-final-review-readiness`. Do not add `codex/` unless explicitly requested.

During implementation, set top-level `enable = false` in `nix-ci.nix` unless
explicitly instructed otherwise. Read-only investigations do not change it.
Before final review, follow the enable, validate, upload, push, and `final-review`
label sequence in [CONTRIBUTING.typ](CONTRIBUTING.typ#ci-readiness). On `itphlies`,
reuse/download matching NixCI cache outputs and finish `just ci-checks-and-upload`
before pushing CI-enabled work, so NixCI can reuse the results.

Ask for approval before changing a test expectation only when the change is
non-trivial or deep. Fix routine fixture mistakes and straightforward obsolete
expectations without asking, preserving the test's validation purpose.

Prefer completion notifications over actively polling checks, uploads, or remote
CI. Preserve the run/commit identity and logs, and resume on completion, failure,
or required input. Never claim a completion wake-up is configured unless it is;
see [CONTRIBUTING.typ](CONTRIBUTING.typ#ci-completion).

Search existing abstractions before adding helpers. Preserve useful comments
and test coverage; follow the shared guidance when their intent is unclear.
Before finishing, review the diff for concise, idiomatic code and duplication.
