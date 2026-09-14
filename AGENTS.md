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

Floating-point underflow in exponentially suppressed tails may round to zero,
including when this makes stability comparisons trivially pass. Apply routine
underflow-to-zero corrections and update their test expectations without asking
again. Preserve meaningful contributions by combining numerical factors before
rounding; overflow and invalid numerical operations remain separate issues.

Before adding helper functions, structs, or methods, check the codebase for
similar use cases, and whether the functionality is already provided by the
existing code or only needs a small adjustment/API change. When adding a new
helper, confirm with the codebase maintainers that the functionality is not
already provided by an existing helper.

Search existing abstractions before adding helpers. Preserve useful comments
and test coverage; follow the shared guidance when their intent is unclear.
Before finishing, review the diff for concise, idiomatic code and duplication.
