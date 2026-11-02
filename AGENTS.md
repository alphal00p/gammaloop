# Agent Guidance

Agents and human contributors follow the same repository guidance. Before making
changes, read [CONTRIBUTING.typ](CONTRIBUTING.typ) and treat it as authoritative
for this repository.

Before editing Rust, read the [Rust policy](CONTRIBUTING.typ#rust-style) and
applicable local guidance. Compliance is mandatory for new and substantially
modified code. Check the final diff against the policy, including design rules
that formatting and Clippy cannot enforce, and report validation and any
departures from the documented allowances.

Before editing Python, read the [Python guidance](CONTRIBUTING.typ#python-style)
and use the Ruff formatter. For native tests, builds and PyO3 packaging, use the
[development commands](CONTRIBUTING.typ#development-commands): `just test CRATE`,
`just build CRATE`, and `just python-wheel` / `just python-develop`.
The default `backend=auto` selects Nix when available and Cargo/Maturin otherwise;
a Nix failure remains a failure. Rust command construction stays in the shared
shell dispatcher. Explicit `backend=dyndrv` uses `cargo-dyndrv-worker` from the
fork pinned in `flake.lock` at `c1935d6`; auto never selects it. See the
[dyndrv workflow](CONTRIBUTING.typ#dyndrv-development) for qualification limits,
session evidence and cleanup.
Native commands and wheels share `.config/build-profile`; use
`just profile=dev COMMAND` or `just profile=release COMMAND` for explicit modes.
Use `backend=cargo` to explicitly select iteration with host tools. Preserve
`SYMBOLICA_LICENSE` and the documented environment and wheel-lifetime rules.
Nix prepares reusable archives and executable runners; licensed tests, doctests
and documentation run as the invoking user on every invocation. Building or
substituting a runner is preparation only, never evidence that its tests passed.
Keep licenses out of Nix evaluation and derivations; use the runtime license
loading described in the development commands.

If instructions conflict, prefer the more specific local guidance and preserve
the user's current work unless explicitly asked to change it.

Before task edits, builds, tests, or experiments, create or resume your own
feature workspace at `../<feature-name>`, relative to `main`. Keep `main`'s
checked-out files read-only for implementation. Reuse only a workspace owned by
the same task; see [CONTRIBUTING.typ](CONTRIBUTING.typ#agent-workspaces) for setup
and resumption rules.

Clean up your own temporary files, abandoned edits, and local processes as work
continues. Preserve other contributors' work, shared caches, useful tests, and
artifacts still needed for debugging, CI, or review. Follow
[CONTRIBUTING.typ](CONTRIBUTING.typ#agent-cleanup) before removing a workspace.

Before creating any new jj change, always ask the user which author name and
email to use, and wait for their answer. Set that identity explicitly when
creating the change.

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
