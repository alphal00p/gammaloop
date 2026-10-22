# Agent Guidance

Agents and human contributors follow the same repository guidance. Before making
changes, read [CONTRIBUTING.typ](CONTRIBUTING.typ) and treat it as authoritative
for this repository.

Before editing Rust, read the [Rust policy](CONTRIBUTING.typ#rust-style) and
applicable local guidance. Compliance is mandatory for new and substantially
modified code. Check the final diff against the policy, including design rules
that formatting and Clippy cannot enforce, and report validation and any
departures from the documented allowances.

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
non-trivial or deep and has not already been authorized. Fix routine fixture
mistakes and straightforward obsolete expectations without asking, preserving
the test's validation purpose.

Prefer completion notifications over actively polling checks, uploads, or remote
CI. Preserve the run/commit identity and logs, and resume on completion, failure,
or required input. Never claim a completion wake-up is configured unless it is;
see [CONTRIBUTING.typ](CONTRIBUTING.typ#ci-completion).

Follow the [numeric precision and underflow policy](CONTRIBUTING.typ#numeric-precision).
Routine underflow-to-zero corrections and their test expectation updates do not
require renewed approval. Preserve meaningful contributions; overflow and
invalid numerical operations remain separate issues.

Before adding helpers, structs, methods, traits, or adapters, search existing
abstractions and similar use cases. Reuse or adjust the existing owner when
appropriate. Ask maintainers only if ownership or duplication remains unclear
after that search.

Preserve useful comments and test coverage; follow the shared guidance when
their intent is unclear. Before finishing, review the diff for concise,
idiomatic code and duplication.
