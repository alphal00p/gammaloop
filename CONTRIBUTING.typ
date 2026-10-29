= Contributing

These guidelines are shared by human contributors and coding agents. They cover
the project workflow, the engineering style we want, and the subsystem-specific
constraints that are easy to miss during refactors.

== Core Engineering Principles

These rules are intentionally broad and should shape most code changes.

- Prefer one semantic boundary per concept. Before adding a trait, helper, mode,
  or adapter, look for the existing abstraction that already owns the behavior
  and extend that instead. Before adding helpers, structs, or methods, search
  similar use cases; ask maintainers only if ownership or duplication remains
  unclear after that search.
- Collapse duplicate abstractions aggressively. If two names describe the same
  responsibility, merge them and update call sites rather than keeping
  compatibility shims.
- Prefer methods, trait methods, or impl-associated functions over bare free
  functions. A helper should live on the type or trait whose invariant it uses.
- Do not preserve internal API compatibility unless explicitly requested.
  Breaking internal Rust/Python APIs is acceptable when it removes duplication
  or clarifies ownership.
- Avoid compatibility fallbacks and "old path plus new path" designs. Pick the
  clearer model and migrate usages.
- Keep changes minimal but complete: remove obsolete code paths, imports, docs,
  and call sites in the same change. Preserve useful comments and rationale
  when moving code; remove demonstrably obsolete comments with the behavior
  they describe, and ask if their relevance is unclear.
- Prefer concise, idiomatic code and remove duplication. Review the diff for
  readability before finishing; line counts are not a quality target.
- Prefer concrete return types over sentinel-style APIs. Avoid `Option` where
  `None` means "nothing happened" if the caller always needs a usable result.
- Make transformations explicit and idempotent. If recursive rewriting is
  involved, include a fixed-point check or another clear loop guard.
- Avoid relying on fragile errors for control flow when a direct predicate or
  structural check would express the intent better.
- Do not add broad helper layers just to avoid touching call sites. Update the
  call sites when the API should change.


== Investigation Discipline

- Before fixing a suspected bug, find the existing code path that should own the
  behavior.
- When asked to identify failure modes, do not patch yet; report the concrete
  failing cases and where they enter the code.
- When asked for a plan, do not make code changes until the plan is accepted.
- If a requested change looks like only a rename, stop and check whether there
  is a deeper design change implied.
- Ignore unrelated untracked local artifacts by default (editor swap files,
  scratch docs, local example edits, profiling outputs, etc.) unless the task
  clearly requires them.

=== Discrepancy Triage
<discrepancy-triage>
Before investigating an implementation mismatch at a low level, first
perform a bird\'s-eye comparison of the two complete pipelines. This
step is mandatory: it prevents a difference in a shallow representation
boundary from being debugged as a problem in shared algebra or physics.

+ Write down the inputs, requested modes, and ordered processing stages
  for both routes. Mark every stage as shared or different, and identify
  the last shared boundary and the first boundary that can produce
  different state.
+ Record the smallest result matrix that distinguishes the routes.
  Toggle one independent feature at a time, such as empty versus
  nonempty UV forest, local versus integrated counterterms, threshold
  subtraction, cut, residue order, LMB channel, orientation, evaluator
  mode, and numeric precision.
+ Apply logical exclusions before proposing causes. A mechanism shared
  by two routes cannot explain a difference between them unless the
  inputs reaching that mechanism have already diverged. An independently
  validated downstream engine is not a candidate until its two actual
  inputs are shown to differ or an identical-input A/B test fails.
+ Inspect the representation at the first differing boundary before
  tracing deeper code. Prefer durable artifacts such as
  `save standalone --json`, UV forest exports, generation reports, and
  structured debug logs. Compare branch counts, keys, maps, selectors,
  and factorized coefficients before comparing only final floating-point
  totals.
+ Reduce the mismatch hierarchically: total, graph, forest, cut, residue
  order, LMB, orientation, then individual term. Stop as soon as the
  first unequal pair is found and make that pair the reproducer.
+ For selector-local versus explicit-sum representations, compare
  `sum(selector * body)` with `sum(body)` directly. Evaluate the
  complete selector truth table and prove that each explicit branch is
  selected exactly once. Do not investigate contour, residue,
  reconstruction, or numerator machinery until this shallow
  partition-of-unity comparison passes.

Diagnostic normalization must preserve the factorization of graph
numerators, including test-only copies. Denominator-only algebra and
finite tensor-component contraction are distinct from distributing a
graph numerator. Each progress report for a discrepancy should state:
what is shared, what first differs, what has been excluded, the smallest
current reproducer, and the single next comparison that will reduce it.

==== Projected local-4D UV reconstruction stop rule
<projected-local-4d-uv-reconstruction-stop-rule>
For a projected local-4D UV-to-CFF mismatch, reconstruction must be
certified #emph[before] investigating CFF recursion, contour signs, or
residue aggregation. This is a mandatory correctness boundary, not an
optional diagnostic:

+ Retain raw post-Taylor sectors in the compatible hard sub-LMB,
  including original edge owners and provenance roles for subsequent
  outer Taylor operations. Build canonical algebra only in a separate
  projection view.
+ Resolve denominator classes by exact signed routing, mass, full
  polynomial (including its prescription), and component domain. Keep
  physical source occurrences as a separate incidence witness. Construct
  the UV skeleton from that witness; derivatives may add serial copies
  of retained lines. Contract equivalent channels only with a certified
  serial-path or pure-cycle incidence. Nonadjacent equivalent channels
  retain their original incidence; do not infer a graph from a signature
  matrix.
+ Completed hard provenance roles zero and one may transfer to a
  certified denominator class and use occurrences of that class,
  including original numerator factors. Preserve the literal
  odd-momentum sign independently of denominator evenness. A hard
  carrier without an unambiguous surviving pole class stays a fixed
  affine carrier; physical-source and soft provenance keep their
  existing meanings. New soft factors retain explicit provenance until
  outer-CFF assembly; only those factors may use an exactly certified
  off-shell routing through active cograph edges, including every fixed
  external shift.
+ Apply the exact immutable production assignment plan, including the
  signed hard/raw/parsed conversion `H = h R`, `P = r R`, and
  `H^0 = h r P^0`.
+ Substitute the source post-Taylor numerator and the reconstructed
  UV-EMR numerator into one neutral set of formal loop four-momenta (and
  the same fixed external data) and require their exact symbolic
  difference to vanish.
+ Independently require equality of denominator momentum, mass,
  multiplicity, and component domain, using `D(Q) = D(-Q)` only on the
  denominator side.

Only after both exact certificates pass may a discrepancy be attributed
to generalized-CFF input normalization, CFF generation, component
composition, or residue aggregation. The common LMB is only a coordinate
chart for this proof; it never supplies EMR ownership or CFF rank
capacity. Correctness of the EMR rewrite precedes dispatch optimization.
Production compares at most three certified proposals by their actual
native generated source-map row count, before surface conversion or
cut/host selection. Rank orders the proposals and breaks count ties; it
is not a substitute for that count or a proof of global optimality. The
selected payload and its exact assignment must stay together. Graph
numerators remain factorized, including diagnostic copies. Establish the
identity with factor-preserving rewrites and exact cancellation;
evaluations at explicit coordinates provide additional checks, not a
general symbolic proof. The worked GL04 `1zs/T2` certificate is
maintained in
#link("docs/architecture/exact-powered-denominator-cff-lifting.typ#worked-live-reproducer-gl04-temporal-square-1zst2")[`docs/architecture/exact-powered-denominator-cff-lifting.typ`];.

== Debug Logging Pattern

- Prefer `debug_tags!` plus the log filter environment variables over ad hoc
  debug files, one-off `println!`, or custom dump env vars. This keeps
  investigative output routable to display and JSONL log sinks with the same
  filtering syntax.
- Add semantic tags that describe the pipeline and subsystem, for example:
  ```rust
  debug_tags!(#uv, #integrated, #vakint, #trace;
      stage = "to_vakint_integrand_after_den_to_prop",
      reduced = %reduced_label,
      display.preview = %expr.log_print(Some(120)),
      file.expr = %expr.to_plain_string(),
      "Vakint trace"
  );
  ```
- Use a `stage` field for fine-grained ordering within a code path. Use
  concise display fields for terminal output and `file.*` fields for large or
  exact payloads. `file.*` fields are hidden from the display sink and written
  as structured JSONL logfile fields with the routing prefix stripped
  (`file.expr` is stored as `expr`).
- To capture a focused test trace to file while keeping stderr quiet, run for
  example:
  ```bash
  GL_DISPLAY_FILTER=off \
  GL_LOGFILE_FILTER='[{#uv,#integrated,#vakint,#trace}]=debug' \
  GL_TEST_LOG_DIR=/tmp/gammaloop-test-logs \
  cargo test -p gammalooprs --test test_renormalization --profile dev-optim \
    finite_part_ghost_3loop_renormalization_sum_origin_mwe -- --ignored --nocapture
  ```
- Use `GL_ALL_LOG_FILTER` only when display and logfile should receive the same
  filter. Otherwise prefer `GL_DISPLAY_FILTER` and `GL_LOGFILE_FILTER` so large
  `file.*` payloads can be kept out of terminal output.

== Code Quality Baseline

=== Rust
<rust-style>

==== Scope And Authority

This document defines GammaLoop's mandatory, self-contained Rust policy for
agents and human contributors. The
#link("https://google.github.io/styleguide/rustguide.html")[Google Rust Style Guide]
and #link("https://rust-lang.github.io/api-guidelines/")[Rust API Guidelines]
are source material for this policy; their rules and future updates are adopted
only through explicit changes to this document. Read this policy and applicable
local guidance before editing Rust.

Apply the policy to new and substantially modified code. Existing violations
are not precedents for new code. Keep unrelated cleanup outside the task; a
small edit does not require rewriting its surrounding module.

State exceptions alongside the affected rule, with their scope and rationale.
Numerical precision, mathematical notation, generated code, and interoperability
may need specific allowances; none is a blanket exemption. Explain any departure
outside the documented allowances in the change's review or handoff, identifying
the affected rule and why the departure is necessary.

Enforce mechanical requirements through formatting checks and selected Clippy
lints. Review the final diff against the design rules as well: passing automated
checks does not establish compliance with the whole policy. Report the relevant
validation performed and any checks that could not be completed.

==== Crates And Modules

- Prefer small, cohesive crates with narrow dependencies. When adding a coherent
  capability or substantially restructuring a subsystem, actively consider
  extracting a focused crate instead of growing an already broad crate. Use
  ownership, reuse, and dependency boundaries to guide the split. Reducing
  recompilation, enabling parallel compilation, and improving reuse of unchanged
  Nix build outputs are explicit design goals; crate boundaries and build-source
  filtering must support those goals. Existing large crates are not a target
  size for new code.
- Use conventional Cargo layout with nested modules. Prefer `foo/mod.rs` for
  modules with child files, keeping those files together under `foo/`. Leaf
  modules may use `foo.rs`. Avoid unrelated layout-only renames when extending
  an existing module tree.
- Keep parent modules focused on declarations and re-exports. Put substantial
  algorithms in child modules. Small shared types and coordination code may
  live in the parent when they connect its children.
- Keep implementation modules private. Use `pub(crate)` for items that need
  access elsewhere in the crate, and make modules public only when callers need
  that organization. Prefer explicit re-exports of the public items callers
  need over exposing implementation structure.

==== Imports And Paths

- Put imports at module scope, including trait imports. Do not place `use`
  statements inside functions or other local blocks. Imports for a test module
  belong at that module's scope.
- Use explicit imports. `use symbolica::prelude::*;` is allowed. Test modules
  may also use `use super::*;`. Other dependency preludes do not inherit the
  Symbolica exception. Import enum variants explicitly if they are needed as
  unqualified names; other glob imports are not allowed.
- Each `use` statement imports items from one source module. Group names from
  that module, as in `use std::fmt::{Debug, Display};`, and use separate
  statements for different modules instead of nested import trees.
- Outside the macro-hygiene exception below, always import referenced items and
  use their imported names instead of fully or partially qualified module or
  crate paths. Names already in scope through a language prelude or an allowed
  glob need no additional import. Resolve name collisions with import aliases,
  such as `use std::io::Error as IoError;`,
  rather than using `io::Error` or `std::io::Error` at the point of use.
- In ordinary code, the only exception to module-path qualification is calling
  a free function from another module: import its module and call
  `module::function(...)`. For example, use `mem::take(...)` after
  `use std::mem;`, rather than writing
  `std::mem::take(...)` or importing `take` directly. Alias the imported module
  if its name collides with another import.
- Use ordinary method-call syntax. Associated items on an imported type or
  `Self`, such as `Type::new()` or `Self::Item`, are distinct from module or
  crate qualification and remain valid.

==== Formatting

- Format Rust with the repository's pinned stable `rustfmt`, available through
  the Nix development environment. Use `just fmt` and the standard formatting
  with a 100-column target. Follow the formatter's output.
- Allow `#[rustfmt::skip]` only for meaningful mathematical alignment, such as
  matrices or tables, or for generated code. Keep the exemption narrowly scoped
  and explain its reason alongside the attribute or in the generator.
- Use lowercase hexadecimal digits in numeric literals. Write invisible
  characters in literals using explicit escapes, such as `\u{200d}` for a
  zero-width joiner, so their presence is visible in source.
- Keep formatting compatible with the pinned stable toolchain. Rules that its
  formatter cannot enforce, including import grouping where unsupported, remain
  review requirements. Do not require nightly-only formatting options.

==== Ownership And Borrowing

- Make ownership intentional. Prefer uniquely owned values and borrowing. Move
  values when transferring ownership; borrow when a function only needs
  temporary access.
- Prefer borrowing over cloning. First consider whether borrowing or changing
  the data flow avoids a copy. Clone when an independent value is required, or
  when borrowing would make the design substantially harder to maintain. The
  simplicity exception requires measurements showing that the cloning cost is
  acceptable for the relevant workload. Use representative benchmarks or
  profiles and record the workload and result in the review or handoff; an
  assumption that the copy is cheap is not sufficient. Check the actual cloning
  cost, especially for symbolic expressions, tensors, and evaluation loops.
- Use `Rc` and `Arc` for deliberate shared ownership when multiple owners need
  the same data and borrowing would be awkward. Cloning a reference-counted
  handle is expected and does not copy its underlying data, but still has a
  runtime cost. Apply the same measured-performance condition when choosing
  shared ownership for simplicity.

==== Function Inputs And Outputs

- Accept flexible borrowed inputs when only their contents matter: prefer
  `&str`, `&[T]`, `&T`, and `Option<&T>` over references to their owning
  containers.
- When iteration is the contract, reusable library APIs should accept
  `impl IntoIterator`. Use slices when indexing, contiguous storage, or a
  simpler internal interface matters.
- Return iterators when results naturally stream. Return a collection when
  the implementation already builds one, callers need indexing, or producing
  an owned collection is the operation's purpose.
- Use `impl Into<T>` selectively where a public API benefits from accepting
  several input types. Prefer concrete inputs internally; weigh call-site
  convenience against compilation cost, code size, and inference complexity.
- Use named enums for ambiguous boolean flags. A clear boolean operation such
  as `set_enabled(bool)` remains appropriate. Prefer a configuration object or
  a small builder on the owning type for several independent options. Named
  constructors are sufficient when there are only a few modes.

==== Errors And Diagnostics

- Use `Option` for meaningful absence. Use `Result` when callers need to
  understand why an operation failed. Avoid `None` as an ambiguous signal that
  nothing happened.
- Return errors for expected failures, including invalid user input,
  unsupported operations, I/O failures, and anticipated numerical failures.
  Propagate them with useful context.
- Use `thiserror` for typed errors in reusable Rust components, including their
  public APIs. Expose the distinctions callers need to handle failures.
- Use `eyre` results and contextual reports in application orchestration, such
  as CLI commands and Python-facing workflows. This boundary depends on the
  component's responsibility, not whether its crate name contains "api".
  Reusable library APIs still use typed errors.
- Preserve underlying errors as sources when adding context. Convert errors
  explicitly to CLI diagnostics or Python exceptions at the presentation
  boundary, preserving distinctions callers need to act on. Do not flatten
  structured errors into message strings while propagating them internally.
- At an operation boundary where an internal failure can be contained and the
  operation cleanly abandoned, prefer returning an explicit internal-error
  report. Identify it as a bug and preserve the context needed to investigate.
  Retain panics for broken internal invariants; never continue using potentially
  invalid state or turn an internal failure into a successful fallback result.
- Prefer pattern matching when it proves a value is present. When a panic is
  justified by an invariant, use `expect` with an explanation of that invariant
  and give assertions useful failure context. Tests may use straightforward
  assertions and `unwrap`.
- Use `catch_unwind` only at explicit outer boundaries for cleanup and reporting.
  The affected operation must fail; catching a panic does not establish that
  its state is safe to reuse. Do not use panic catching as ordinary error
  handling or resume work with potentially invalid state.

==== Type Design

- Use newtypes for distinct meanings or invariants, such as index domains,
  units, or validated values that could otherwise be confused. Reuse existing
  domain types and typed indices first.
- Use type aliases for readability when semantic identity stays the same.
  An alias does not establish a distinct type or enforce a new invariant.
- Prefer named fields for multi-field data. Single-field newtypes and small,
  obvious local tuples remain appropriate. Avoid APIs that make callers
  remember the meanings of tuple positions.
- Keep fields private when construction or mutation requires validation.
  Plain data containers may expose fields. Avoid getters and setters that add
  no useful contract.
- Keep mode and strategy enums private unless callers need to inspect or match
  on them.

==== Traits And Dispatch

- Prefer inherent methods for behavior belonging to one type. Introduce traits
  for meaningful generic contracts, interchangeable implementations, or
  extension methods, rather than merely organizing code.
- Implement standard traits such as `From`, `TryFrom`, `FromStr`, and `Iterator`
  when their contracts fit. Preserve GammaLoop's precision-aware semantics
  when choosing conversions.
- Eagerly implement applicable common traits for the types a crate defines,
  following the adopted
  #link("https://rust-lang.github.io/api-guidelines/interoperability.html#c-common-traits")[C-COMMON-TRAITS]
  guidance. Consider `Copy`, `Clone`, `Eq`, `PartialEq`, `Ord`, `PartialOrd`,
  `Hash`, `Debug`, `Display`, and `Default` before a current caller needs them:
  downstream crates often cannot add these implementations themselves.
  Applicability must respect the numerical and standard-trait contracts below.
  Prefer derives when they produce the correct behavior, and avoid unnecessary
  generic bounds.
- Prefer generics for reusable libraries and performance-sensitive kernels.
  Use `dyn Trait` for runtime-selected implementations or when reducing
  compilation cost and code size justifies it. Measure performance-driven
  tradeoffs.
- Use extension traits deliberately and name them `SomethingExt`. For types
  owned within the workspace, prefer adding functionality to the owning type.
- Keep shared behavior behind the primary trait method when a type implements a
  trait. Put mode selection in the implementing type's configuration rather
  than exposing parallel helper methods that bypass the trait.

==== Standard Trait Contracts

- Keep approximate comparisons in `ApproxEq` or explicitly named tolerance
  checks. `PartialEq` must respect symmetry and transitivity; implement `Eq`
  only when reflexivity also holds.
- Implement `Ord` only with a defined total ordering consistent with equality
  and `PartialOrd`. Numerical wrappers must define their NaN and signed-zero
  handling explicitly. Do not silently treat incomparable values as equal.
- Keep hashing consistent with equality: equal keys must hash identically.
  Review custom equality and hashing together, particularly for numerical
  values and cache keys.
- Use `Copy` for small, inexpensive value types, and `Clone` when duplication
  has meaningful semantics. Implementing either trait does not remove the
  performance requirements on copying. Do not add `Copy` bounds to
  precision-generic APIs merely because fixed-precision backends support them;
  arbitrary-precision values may own storage and have different cloning costs.
- Implement `Default` only for a valid, unambiguous default. Configuration
  objects should support it when sensible. If construction requires a precision
  or model context, use a context-aware constructor instead of inventing a
  default. Keep settings defaults aligned with serialization and default
  suppression as required by the settings policy below.

==== Concurrency

- Reuse existing parallel infrastructure. Prefer the established Rayon pools
  for numerical work and respect the configured core budget. Avoid creating
  threads or independent pools inside frequently called operations.
- Prefer worker-owned state and share immutable data where appropriate. Keep
  shared mutation and lock durations small. Avoid holding locks across
  callbacks, blocking I/O, or expensive computation.
- Treat poisoned locks as failures. Do not discard the poisoning error and
  continue computing with potentially invalid state.
- Give background workers an explicit lifecycle: define ownership,
  shutdown/cancellation, and error propagation. Join workers or otherwise
  ensure their completion is accounted for.
- Make reproducibility requirements explicit. Parallel changes must preserve
  the promised random-stream assignment and numerical reproducibility.
  Validate changes to reduction order against those requirements, alongside
  performance.

==== Macros And Code Generation

- Prefer ordinary functions, types, and traits when they express the
  abstraction clearly.
- Use `macro_rules!` for useful structural repetition, such as repeated
  numerical implementations.
- Introduce procedural macros only when their capabilities justify the
  complexity. Established derives such as Serde and `thiserror` remain normal
  tools.
- Keep generated code reproducible. Change its generator or source inputs;
  inspect representative expansions and test meaningful transformation
  behavior.
- Macro expansion templates may use qualified paths required for hygiene,
  such as `$crate::LogMessage`, so an expansion resolves the intended definition
  independently of the caller's imports. This exception is limited to the
  generated expansion; ordinary implementation code and macro call sites still
  follow the import rules above.

==== Naming

- Use standard Rust casing: `snake_case` for functions, methods, variables, and
  modules; `UpperCamelCase` for types, traits, and enum variants; and
  `SCREAMING_SNAKE_CASE` for constants and statics. Name operations for their
  semantic result. Avoid `*_for_<caller>` names unless the result has a
  caller-specific invariant.
- Established mathematical notation and domain abbreviations such as CFF, UV,
  and LMB are appropriate where their meaning is clear. Explain them at public
  boundaries.
- Use `as_` for cheap views, `to_` for conversions that preserve the source, and
  `into_` for consuming conversions.
- Use `new` and `try_new` for ordinary constructors, and meaningful names for
  alternatives. Add an explicit constructor alongside `From` only when it
  materially improves discoverability.
- Receiver-based construction such as `from_usize(&self, ...)` and constant
  methods such as `PI()` are allowed where they preserve the active numerical
  context. Document the reason for associated lint exemptions.
- Keep generic and lifetime names concise but informative. Single letters are
  appropriate for conventional mathematics and simple relationships; use
  descriptive names when several roles would otherwise be hard to distinguish.

==== Rust Documentation

- Document public API behavior and non-obvious internal invariants. Numerical
  API contracts should explain relevant units, conventions, precision
  assumptions, and data shapes. Avoid comments that merely repeat names.
- Include `Errors`, `Panics`, and `Safety` sections where applicable. Unsafe
  APIs must state what callers must guarantee.
- Cover every public API item with an example exercising its functionality,
  following the adopted
  #link("https://rust-lang.github.io/api-guidelines/documentation.html#c-example")[C-EXAMPLE]
  guidance. This includes public modules, traits, structs, enums, functions,
  methods, macros, and type definitions. A link to a relevant example on another
  item can provide that coverage without duplication. Apply the requirement
  within reason: examples should show useful behavior or why an API is needed,
  rather than mechanically illustrating basic Rust syntax. Prefer runnable
  examples and doctests where practical, and demonstrate normal error
  propagation with `?`.
- Keep Rust API documentation beside the code in `///` and `//!` comments using
  Rustdoc Markdown. Standalone manuals, architecture notes, and contributor
  guidance remain Typst. Use ordinary comments for implementation rationale.

==== Unsafe Code And Toolchain

- Expect unsafe code to be rare. GammaLoop's numerical algorithms and
  application logic should normally use safe Rust. Introduce `unsafe` only for
  a concrete need in low-level operations or FFI, or when representative
  benchmarks justify an optimization that safe alternatives cannot adequately
  provide. Explain why the safe alternatives are insufficient.
- Keep unsafe operations small and behind an abstraction that enforces their
  requirements. Separate them from numerical algorithms and application logic
  where practical.
- Each unsafe block or implementation needs a `SAFETY:` comment explaining why
  its obligations hold. Unsafe APIs must document caller obligations.
- Validate relevant invariants with targeted tests and, where applicable,
  Miri or sanitizers. These checks supplement the safety argument; passing
  tests does not establish soundness.
- Keep production code on the pinned stable toolchain, avoiding nightly-only
  language and library features. Specialized validation tools may use a
  separate toolchain without changing production requirements.

==== API Conventions And Initialization

- Use predictable accessors: `foo()`, `foo_mut()`, `set_foo()`, `take_foo()`, and
  `replace_foo()` according to their behavior. Provide expected companions such
  as `is_empty()` with `len()`.
- Overload operators only for meaningful mathematical or logical operations.
  Document conventions where the meaning could be ambiguous.
- Use `#[must_use]` when ignoring an output is probably a mistake, including
  pure transformations and lazy computations. Avoid redundant annotations
  when the return type already provides the warning.
- Make `const fn` a deliberate API choice for meaningful compile-time
  evaluation. Account for callers when changing that capability.
- Use `#[non_exhaustive]` selectively at external compatibility boundaries.
  Do not apply it mechanically to internal types; GammaLoop's internal APIs
  remain freely refactorable.
- Treat dependencies' `#[doc(hidden)]` items as private.
- Keep initialization explicit or controlled on first use; avoid hidden
  execution before `main`. Symbolica's managed global state, thread-local
  workspaces, and library-symbol registration are an explicit allowance.
  Use its module-scope `initialize!` mechanism when symbols must be registered
  as the global state initializes, before expressions use them. Preserve
  initializer dependencies and the existing initializer-aware lazy wrappers
  that prevent re-entry and lock-order problems. Retain the `StateMap` context
  needed to translate imported identifiers into the current state. This
  allowance covers Symbolica integration, not unrelated application startup
  work in registration callbacks.

==== Implementation Readability

- Make inferred types understandable. Add annotations where `.parse()`,
  `.into()`, or complex expressions obscure the result type. Prefer annotating
  a binding when that reads better than a turbofish.
- Choose loops or iterator chains for clarity. Keep concise transformations as
  chains; use loops for complicated control flow or side effects. Name long
  pipelines before consuming them.
- Make matching intentional. Handle variants explicitly when each needs
  consideration. Use `if let` or a deliberate wildcard when other variants are
  irrelevant; do not silently swallow newly added cases.
- Use `format!` and `write!` for interpolation and multi-part text assembly.
  Reuse buffers where allocation cost matters.
- Prefer explicit function composition over function nesting: when a value
  crosses a semantic boundary, bind it or use `From`/`TryFrom` so the
  conversion is visible at the call site.

==== Clippy And Validation

- Run checks appropriate to the affected code. Use `cargo check` for quick
  feedback when useful; it is not a prerequisite for every build or test run.
- Follow Clippy diligently. Run the repository's Clippy checks for relevant
  Rust changes with the pinned toolchain and configured lint set. Address
  warnings in new and substantially modified code before finishing. Report
  unrelated existing warnings and any checks that could not be completed.
- Apply suggestions with attention to semantics, numerical precision, and
  measured performance. Review automated fixes; a lint suggestion does not
  override the contracts in this guide.
- A lint exemption needs a concrete reason, such as an intentional numerical
  convention or a false positive. Name the specific lint, use the smallest
  applicable scope, and document the reason. Do not add broad suppressions or
  weaken lint settings merely to make checks pass.

=== Python

- Use the ruff formatter.

=== Tests

- Never weaken tests merely to make failures disappear. Fix routine fixture
  mistakes and straightforward obsolete expectations without asking,
  preserving the test's validation purpose. Ask before non-trivial or deep
  expectation changes unless they are already authorized, including changes
  required by an explicitly requested behavior change.
- During development, run checks relevant to the change. Before final review,
  prepare and upload the selected CI artifacts; NixCI runs the full suite.
- Install `cargo-nextest` 0.9.137 or newer. The repository configuration
  enforces this minimum; update an existing installation with
  `cargo nextest self update`.
- Rust integration tests live in `tests/` with shared fixtures in
  `tests/resources/`.
- Use `cargo nextest` via `just test TEST_NAME` for targeted integration runs.
- Pass emulated libtest arguments after `--`, for example
  `just test_gammaloop -- --skip TEST_NAME`.
- Put targeted unit tests next to the implementation they exercise.
- Use broader integration tests only for cross-module behavior.
- Add tests for edge cases that motivated the change, especially when collapsing
  duplicated logic.
- Numerator oracles must retain physical locality: edge factors depend only on
  their edge momentum, and vertex factors only on incident momenta. Non-local
  algebraic diagnostics do not certify UV subtraction or graph reconstruction.

=== Docs

- Documentation should explain conventions and algorithmic intent, not restate
  function names. Avoid tautological docs like "apply the convention"; say what
  the convention does.
- Keep docs synchronized with API refactors; remove references to deleted
  traits, functions, and modes.

==== Documentation Source Format

- Typst (`.typ`) is the canonical source format for standalone authored
  repository documentation, including contributor guidance, manuals,
  architecture notes, investigations, plans, and changelogs. Never add a new
  Markdown/HTML prose source or a parallel copy.
- Files named exactly `AGENTS.md` and `README.md` are the only Markdown
  compatibility exceptions: agent tooling discovers the former, while GitHub,
  Cargo, and package indexes discover the latter. Keep those files concise and
  link them to the canonical Typst documentation.
- Rustdoc comments embedded in Rust source use Markdown as required by
  Rustdoc; the standalone-document format rule does not apply to them.
- Do not add a separately maintained `README.typ`; the compatibility README is
  its own deliberately small source. This applies in every directory. Changing
  extension spelling or case (`.MD`, `.mdx`, `.HTML`, `.TYP`, and similar) does
  not create another exception; canonical Typst sources use lowercase `.typ`.
- Existing Markdown and HTML prose is legacy migration debt, not precedent for
  new files. Any change to its authored content must migrate the canonical
  source to Typst and update every registry and link in the same change instead
  of retaining both formats.
- Build outputs may use a downstream-required format, but generated output is
  not a second source and must not be checked in as hand-maintained prose.
- Architecture docs live in `docs/architecture/` and are classified by
  `docs/developers.toml`. New entries must point directly to `.typ` sources.
  Each source begins with exactly one `= Title`, matching its registry title;
  the builder supplies the page-level heading and verifies that contract.
- `docs/legacy-prose.toml` records that the pre-policy prose migration is
  complete. Its inventory must remain empty; the documentation check rejects
  additions and parallel editable copies.

==== Typst Graph Debugging

When debugging Linnest layout output from Typst, expose the graph records through
Typst metadata and inspect them with `typst query`. The metadata element must be
wrapped in markup, because Typst labels can only be attached in markup mode.

```typst
[#metadata(graph.nodes(g)) <linnest-tree-nodes>]
[#metadata(graph.edges(g)) <linnest-tree-edges>]
```

Then query the values from the same root used for compilation:

```bash
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-nodes>' --field value --pretty
typst query --root crates crates/linnest/typst/examples/tree.typ '<linnest-tree-edges>' --field value --pretty
```

This is useful for checking generated `pos`, `label-pos`, `layout-width`,
`layout-height`, `label-width`, and `label-height` values without inferring them
from the rendered PDF.

== Repository Map

- `crates/gammalooprs/src/` holds the core Rust implementation.
- `crates/gammaloop-api/` is the workspace member for the Rust API and Python
  bindings (`crates/gammaloop-api/src` and
  `crates/gammaloop-api/python/gammaloop`).
- `tests/` contains the integration-test crate and shared fixtures.
- `examples/` includes command cards and notebooks for end-to-end runs.
- `bin/` contains the process-tree memory watchdog used by guarded test runs.
  The repository-root `gammaloop` script launches the CLI.
- `assets/` stores schemas/data and model files (`assets/models/`).
- Rust benches live under `crates/gammalooprs/benches/`.

== Common Commands

- Build Rust CLI in dev-optim:
  ```bash
  just build-cli
  ```
- Build Rust CLI with stable Python ABI:
  ```bash
  just build-cli-abi
  ```
- Build Python API bindings:
  ```bash
  just build-api
  ```
- Build Python API bindings with stable ABI:
  ```bash
  just build-api-abi
  ```
- Format and lint:
  ```bash
  just fmt
  just clippy
  ```
- If you run into a macOS linking issue complaining about a missing `__emul...`
  symbol, try building with `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T`; see
  #link("crates/gammaloop-api/build.rs")[`crates/gammaloop-api/build.rs`]
  for the impact of this setting.

The root `justfile` keeps build and lint commands. Test, NixCI and drawing recipes
are imported from `just/tests.just`, `just/ci.just` and `just/drawing.just`; run all
commands from the repository as before. Use `just --list` to see them.

== CI Readiness
<ci-readiness>

The `final-review` label identifies PRs ready for final review. Merging requires
an open, non-draft PR with this label, explicit top-level `enable = true` in
`nix-ci.nix`, and successful required CI checks. Draft and unlabeled PRs are
valid work in progress; their readiness gate intentionally fails to block merging.

During implementation, agents set `enable = false` unless explicitly instructed
otherwise. Human contributors may enable CI earlier. Read-only investigations
do not change the toggle. Keep `enable = true` on `main`.

Before requesting final review:

+ Set top-level `enable = true` in `nix-ci.nix`.
+ On `itphlies`, reuse/download matching cache outputs through Nix substitution.
  Run static checks and prepare/upload the selected runtime producers against
  the final code when credentials are available, following the cache workflow
  below. Runtime tests run in NixCI after the push.
+ Push the prepared code, mark the PR non-draft, then apply `final-review`.

When returning to development, remove `final-review` and commit
`enable = false`. The label gates merging; the committed toggle controls NixCI
work independently, and changing either does not cancel already-running NixCI
jobs. Full GitHub Actions build workflows run automatically on pushes to `main`, not
on PRs, merge groups, or feature-branch pushes. Existing manual runs, schedules,
and release-tag publishing remain available. The lightweight NixCI readiness
check still runs on PR and merge-group events because it is required to merge.
After configuration validation, it runs a paired performance comparison for
relevant changes on same-repository final-review PRs and merge-group commits.
Documentation-only PR changes skip that comparison. The separate main-push
workflow records fresh gammabench measurements for every merge. See #link("docs/architecture/ci.typ")[CI maintenance and
measurements] for the workload, coverage limits and trigger policy.

Edit only the top-level toggle manually: `just ci-update` preserves it while
regenerating the scheduling configuration. Missing or nonboolean toggles and
other generated-content drift fail validation. `dependency-discovery.enable`
is independent and must not be used to disable CI. Local checks and uploads
remain available with NixCI disabled.

== NixCI Cache
<nixci-cache>

See #link("docs/architecture/ci.typ")[CI maintenance and measurements] for the build design, measured
results, branch migration and `just ci-report` usage.

Before final review, contributors and agents must run the selected static checks
locally and upload compiled runtime producers when cache credentials are available.
NixCI runs the selected runtime suite after the push.
Ordinary development pushes require checks relevant to the change, not the full
suite and upload cycle.
Work from the repository root; enter `nix develop` if you need the pinned Just
and build tools. Stage newly added files first so the Git-backed flake sees them.

`just check` runs Cargo checking only. For the selected CI suite without uploading:

```sh
just ci-checks
```

This runs the selected Rust and Python tests, Clippy, doctests, formatting,
workspace/CI graph checks, and the full five-product documentation check. The
documentation check validates generated inputs, renders latest and snapshot
sites, checks their HTML, and exercises publication behavior. Licensed checks,
including documentation, require `SYMBOLICA_LICENSE` in your local environment.
An unchanged successful Nix check can be reused without executing its tests again.
`nix flake check --impure` additionally builds the CLI, standalone Rustdoc,
persistent Typst renderer and WASM checks exported by the flake.

On `itphlies`, run the following before pushing CI-enabled work, including
updates to a PR already in final review. The host's Nix daemon is configured to
download from the NixCI cache using its signing key and machine credentials.
Leave substitution enabled: the Nix builds in this command automatically reuse
local outputs or download matching cached outputs, build/check what is missing,
then upload the static-check outputs, runtime launchers and producer closures.
There is no need to download the entire cache or force rebuilds. Completing this locally
before the push minimizes work left for NixCI.

Use the same command on other hosts for the final-review pre-push workflow when
cache credentials are available:

```sh
just ci-checks-and-upload
```

This runs static checks and compiles runtime producers without running Rust or
Python tests, doctests or the documentation publication suite. Use `just ci-checks`
explicitly when you also want to run that suite locally.
Wait for the command to succeed and print `CI upload completed` before pushing.
Run it against the final code you intend to push; rerun it if that code changes.
If you lack upload credentials, run checks relevant to the change and mention
in your handoff or PR that the artifacts were not uploaded. NixCI builds any
missing artifacts before running its tests.

After changing dependencies or features, run `cargo hakari generate` and
`cargo hakari verify` in `nix develop`. Commit the updated workspace-hack manifest
and lockfile. After changing workspace manifests or test groups, run
`just ci-update` before the checks and commit the regenerated
`nix/ci-workspace-graph.json` and `nix-ci.nix` with the change.

Put your NixCI token in `~/.netrc`, or point `NIXCI_NETRC` at an existing file:

```sh
export NIXCI_NETRC=/path/to/existing/netrc
```

On Clan-managed hosts, use the configured `NIXCI_NETRC`; a fresh login may be
needed after upload access is enabled. Keep tokens and license values outside
Git. See #link("https://nix-ci.com/documentation/nix-ci-cache")[the NixCI cache documentation].
Uploads run as your user and need no additional Nix privileges. Entering the dev
shell enables no upload hook. Bare Cargo builds do not populate this cache.
An untrusted daemon may warn that it ignores `netrc-file`; the upload client
still uses that file to authenticate with the destination cache.

The upload command realizes the build jobs selected by `nix/ci.nix` and publishes
their runtime closures. These jobs include static checks, lightweight runtime
launchers and the producers that retain test binaries and compiler artifacts.
The command does not select runtime test-result checks. Existing matching build
outputs are published too. It does not select packaging, WASM, dev shells or
unrelated repositories.
Required shared dependencies can still be part of the uploaded closures.

Static checks and artifact preparation share one timed build phase; upload time
is reported separately.
An upload failure returns a nonzero status; fix the reported problem and rerun
`just ci-checks-and-upload`. Prepared outputs remain cached, so retrying does not
require recompiling matching inputs. Use plain `just ci-checks` for local
benchmarks without publishing. Remote reuse also requires matching outputs: a macOS build,
for example, does not replace a Linux build.

NixCI starts each selected runtime launcher even when that launcher is cached.
The launcher requests the corresponding Nix test-result check; a matching cached
result can skip test execution. Uploading compiled producers prepares the inputs
for uncached runtime checks, rather than certifying their results locally.

The flake retains NixCI's substituter and signing-key settings for downloads.
On a shared daemon, an administrator must configure the cache's trust and access
credentials for substitution; enabling user uploads does not configure downloads.

== CI Completion
<ci-completion>

Prefer completion notifications over an agent repeatedly polling local checks,
cache transfers, or remote CI. When the execution environment supports it, run
long commands with captured logs and an exit-status notification, and resume
only on completion, failure, or required input. Local checks and uploads must
still succeed before the dependent push; background execution is not evidence
of success.

After a push, record the commit and CI run or PR URL, then yield instead of
actively tracking unchanged status unless explicitly asked to monitor. Use an
available event-driven completion wake-up to resume the agent. Verify that the
notification is actually registered; if no such mechanism is available, state
that limitation in the handoff instead of promising an automatic follow-up.
A scheduled monitor still polls; do not substitute it for a completion event
without making that distinction explicit. On resuming, verify that the result
belongs to the intended commit before reporting success or investigating failure.

== Version Control Workflow

=== jj

- If the repo is colocated with jj metadata, prefer `jj` over `git`.
- Start new tasks with `jj new` and immediately
  `jj describe -m "<plan>"`.
- Use `jj log --no-pager` and `jj diff --no-pager` to avoid paging.
- Check the current change before editing; if the working copy has changes but
  no description, run `jj diff` then `jj describe -m "<description>"`.
- Update the current change description with `jj describe -m "<description>"`
  as the plan evolves.
- Use `jj diff` to see changed files.
- Reference: #link("https://jj-vcs.github.io/jj/prerelease/cli-reference/")[Jujutsu CLI reference].

=== Commits And PRs

- Use descriptive branch names without agent prefixes, for example
  `ci-final-review-readiness`. Only use `codex/` when explicitly requested.
- Commit messages are short and descriptive, typically lowercase without scopes
  (for example, `remove edge quotes`).
- For user-requested changes, describe the requested intent and resulting
  behavior; do not automatically copy the full user prompt.
- If there is already a commit, prefer making a new one.
- PRs should include a clear summary, test command(s) run, and any example
  command card used to validate behavior.
- If changes touch CLI state, mention `gammaloop_state/` impacts or migration
  steps.

== Internal API Policy

- API stability is currently a non-goal for internal development: prioritize
  maintainability and structure over preserving external-facing module paths.
- Breaking Rust/Python API changes are acceptable when they simplify
  architecture; document major moves in the change description.
- Backward compatibility with old on-disk GammaLoop states is also currently a
  non-goal: it is acceptable to rename/remove legacy state files, layouts, and
  loaders when that simplifies the implementation.
- Do not keep compatibility fallbacks for old state formats/names unless the
  task explicitly asks for migration support.
- Do not hesitate to modify existing functions when that produces a cleaner and
  more concise design than preserving compatibility or leaving old code
  untouched.

== Project-Specific Notes

The following sections are narrower invariants. They are not general style
preferences; they protect behavior that has been easy to regress.

=== Numeric Precision
<numeric-precision>

- In generic or arbitrary-precision code, do not introduce constants or
  intermediate values through `from_f64`, `std::f64::consts::*`, or other lossy
  `f64` routes unless that exact location is an explicit `f64` boundary by
  design, such as persisted settings that are already `f64` or
  histogram/output accumulation that is intentionally `f64`.
- `f64` values originating from user-supplied settings are an allowed boundary:
  converting those setting values into `F<T>` with `from_f64` is acceptable. Do
  not extend that exception to internally generated constants or intermediate
  values.
- When working with `F<T>` or other precision-generic numeric code, build
  constants from an in-scope representative value of the correct type using
  helpers such as `.zero()`, `.one()`, `.epsilon()`, `.PI()`, `.TAU()`,
  `.from_usize()`, `.from_isize()`, etc., so the active precision is preserved
  exactly.
- If a needed constant or operation cannot be expressed without a dubious
  precision-losing conversion, stop and ask for clarification instead of
  guessing.
- Do not explicitly bring `symbolica::domains::float::FloatLike` into scope in
  GammaLoop code. Prefer the GammaLoop wrappers in `crate::utils`, and if
  `F<T>` is missing an operation needed for precision-safe generic code, add or
  use a helper there instead. If it genuinely looks unavoidable to use
  Symbolica's trait directly, stop and ask first.
- Floating-point underflow in exponentially suppressed tails may round to
  zero, including when this makes stability comparisons trivially pass. Apply
  routine underflow-to-zero corrections and update their test expectations
  without asking again. Preserve meaningful contributions by combining
  numerical factors before rounding; overflow and invalid numerical operations
  remain separate issues.
- For stability of an averaged component, its mean absolute fully weighted
  probe value may establish underflow below binary64's smallest normal value.
  Keep this bound separate for each component and observable; never borrow an
  imaginary component's scale to accept a meaningful real discrepancy.

=== Settings Serialization

- The settings serialization model relies on the `SHOWDEFAULTS` escape hatch in
  `gammalooprs::utils::serde_utils` when writing persisted state defaults and
  when building completion catalogs from serialized settings.
- When adding or changing settings fields, verify that `Default::default()` and
  the corresponding `skip_serializing_if` helper stay aligned.
- Custom helpers must route through `show_defaults_helper(...)` or
  `IsDefault::is_default`; otherwise the field will disappear from saved
  `global_settings.toml` / `default_runtime_settings.toml` and from completion.

=== State And CLI

- Treat saved-state detection as manifest-based only. Empty folders or folders
  containing only transient runtime artifacts such as `logs/` are scratch state,
  not legacy state, and must be treated as blank state.
- Do not introduce writes into the state folder outside explicit `save state` /
  `quit -o`, except for logfile tracing when that logger is actually enabled.
- Honor `--read-only-state` consistently. When it is enabled, do not write into
  the state folder and prefer cwd-based fallbacks for transient artifacts that
  would otherwise default into the state.
- The run-card field is `command_blocks` with singular `command`, not
  `commands_blocks`. No compatibility alias is kept.
- Integration workspace resume must restore the exact state at the end of the
  last completed iteration. The authoritative resume data lives in the
  workspace manifest/internal state plus per-slot settings and internal
  observable checkpoints, not in user-facing `integration_result.json`.
- Treat workspace snapshots as completed-iteration checkpoints only: do not
  persist partial-iteration state, and keep user-facing latest snapshots plus
  optional `_iter_0001` history separate from the internal resume state.
- CLI runs create a `gammaloop_state/` directory by default; keep it out of
  commits unless intentionally sharing a reproducible state.
- Use `./gammaloop -s <state_dir>` for isolated experiments.

=== Differential LU / Event Processing

- The GammaLoop current-architecture record registered in
  `docs/developers.toml` is the implemented summary for the differential LU
  stack. Keep it synchronized when changing selectors, observables, event
  grouping, or sample-evaluation output.
- Event grouping semantics are by graph-group, not just by graph: if multiple
  LU graphs share the same `group_id`, all accepted cuts from all of those
  graphs belong in the same retained `EventGroup`.
- Histogram error propagation is based on Monte Carlo samples, not on individual
  observable entries. Treat the full retained event-group list from one
  evaluation as one statistical sample for each histogram.
- Histogram accumulation is sparse. Do not explicitly zero-fill untouched bins.
  The sparse state is:
  - histogram-level `sample_count`;
  - per-bin raw stats `entry_count`, `sum_weights`, `sum_weights_squared`,
    `mitigated_fill_count`;
  - underflow and overflow as full bins with the same raw-stat structure.
- Histogram snapshots are intentionally raw-stat snapshots, not presentation-only
  views. They must remain mergeable and reconstructible back into live histogram
  state. Rust-side helpers already exist for `merge(...)`, `merge_in_place(...)`,
  `rebin(...)`, and reconstruction into accumulator state.
- Sample-evaluation results retain the full event-group list when observables
  are configured, without applying the outer sampler weight or accumulating
  statistics. `general.generate_events` additionally requests diagnostic events
  when no observables are configured. Selectors may force temporary events.
  Accumulate each sample explicitly with its complete outer coefficient; a zero
  coefficient still advances the sample count.
- `general.store_additional_weights_in_event` controls whether
  `additional_weights` are stored on events. The map is a `BTreeMap` keyed by
  lightweight identifiers such as `Original`,
  `ThresholdCounterterm { subset_index }`, and `FullMultiplicativeFactor`.

=== API / Python Interop

- `gammaloop-api` defaults to `cli`, which owns the executable, interactive
  editor, and shell completion generators. Python packaging disables default
  features; reusable state loading and command sessions remain available.
- `ufo-import` enables Python embedding for UFO model imports; JSON model
  imports remain available without it.
- `python` enables the Python API and implies `ufo-import`. `python-stubs`
  adds stub metadata and implies `python`.
- `extension-module` enables importable-extension linking and implies `python`.
  Keep this mode separate from CLI and Rust test builds that embed Python.
- `abi3` selects the stable Python ABI without activating Python support.
  Use `ufo-import,abi3` for a compatible embedding CLI or
  `extension-module,abi3` for a distribution extension. Stable-ABI builds may
  have more overhead than version-specific builds.
- The Rust API has precise endpoints `evaluate_sample_precise` and
  `evaluate_samples_precise`; the Python API intentionally stays `f64`-only.
- Python API end-to-end tests are subprocess-based and assume the contributor
  already built/installed the Python extension in the active environment
  (`maturin develop` or `just build-api`). Do not try to embed Python into the
  Rust test binary for these tests.
- The top-level `.venv` in the repo is a reasonable default development
  environment for running the Python API examples/tests after building the
  extension there.
- If the extension module environment is an issue when compiling the Python API,
  try `build-api-abi`.

=== Environment And Build Quirks

- The `gammaloop` CLI activates its crate-bound application (OEM) license with
  `symbolica::set_application_key!` before calling the reusable API.
  `NO_SYMBOLICA_OEM_LICENSE` is read there via `option_env!`, so
  opting out requires setting it when compiling the binary, not at runtime.
  Rust/Python libraries and their tests do not activate this OEM license; use a
  user license through `SYMBOLICA_LICENSE` for those entry points.
