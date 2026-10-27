= Naming and routing implementation preparation

Read-only source assessment for the already approved decisions 2, 3 and 6.
Baseline: xspzosrl, commit 217c0695ddd936756019ef3b4b57e0685e0e09b9.
The combined workspace remains the empty change swztlkyxkrry, commit
64a27bf078436103a3f548c756d69c53bca7a573. No repository source was edited
and no tests were run for this assessment. Recheck these boundaries after
integration: active A, P2, A6 and S lanes can change their implementation.

== Established owners

- `Request::generate` in `crates/feynkit-generator/src/generation/staged.rs`
  materializes preparations, merges mirror orbits, applies the resolved squared
  order relation, and retains the execution-free search context.
- `PreparedDiagram.compact.canonical_form()` already supplies the ordering used
  by `BuiltDiagramSet::group`. Reuse that cache and the existing ordered map/pool
  infrastructure; do not add a second representation or an enumeration stream.
- `DiagramSet::build` currently assigns names by prepared index, constructs
  numerators through `Generator::diagram_pair`, then applies named or ID-keyed
  LMB overrides. Unknown keys currently disappear without a diagnostic.
- `BuiltDiagramSet::select` currently filters names/IDs after construction.
  It separately tracks forced-cut matches and interrupted progress.
- `SignedRoutedDiagram` already wraps finalized geometry, exact sign and
  multiplicity, exposes structural records, IDs and momentum bases, and protects
  the unset fragments behind a private geometry owner. `DiagramSet::routed`
  currently exists only under `cfg(test)`.
- `Request.classes` is persisted through serde, equality and bincode. `Run`
  already owns ephemeral execution policy; move class policy there in the named
  change, update callers and remove the serde fallback rather than preserving
  two policy owners.
- Current process-definition magic is `GLPROC` version 3; the state manifest is
  version 22. Bump the relevant persisted format once when the request encoding
  changes, using the then-current integrated version. No legacy-format fallback.

== Approved next changes

+ After enumeration/mirror and squared-order filtering, establish canonical
  output order using the existing compact canonical forms in parallel. Named
  output must use that order before grouping, selection or routing. Avoid eager
  DiagramId calculation for compact unnamed consumers. Account for sorting cost
  in the count/materialized/named benchmark boundaries.
+ Resolve name and optional DiagramId-prefix selectors and LMB keys through
  one owner. A moved pin must identify the new name of its diagram; unknown or
  ambiguous references must identify the offending selector. Do not turn
  cancellation into a misleading missing-reference error: distinguish incomplete
  observed output from a completed request with no match. Preserve union
  semantics for selected names/IDs and the meaning of an explicitly empty
  selection. Keep forced-cut matching and squared-order pairing intact.
+ Expose `DiagramSet::route` and `RoutedSet::numerators` with `build` as their
  composition. Use the existing signed-routing geometry and preserve the mapping
  back to prepared inputs for comparison-frame grouping. Route and apply LMBs
  before local numerator construction, so selection can avoid building discarded
  numerators. Check prefactor/projector/zero-colour policy ownership explicitly;
  avoid leaving numerator-only options on an allegedly numerator-free contract.
+ Keep Classes Auto as default only after mode-independent names, the inherited
  dressing-hook rule, and the measured gate pass. Custom particle predicates
  must not reject class-valued states before refinement unless they opt in and
  satisfy the heredity contract.

== Discriminating validation

- Compare name-to-structural-ID maps across Auto/None and 1/8 threads, before and
  after grouping, using the saved P6 representatives and currently failing named
  scalar/gg-hhh cards. Do not accept snapshot rewrites as a naming fix.
- Exercise valid pins, moved pins, unknown names/IDs/LMB keys, ambiguous prefixes,
  an explicitly empty selection, and cancellation during each public boundary.
- Compare `build` with `route` then `numerators` for weight, slots, cuts, routing,
  names/IDs and exact factorized fragments. Include amplitude and cross-section
  requests, mirror orbits, squared-order pairing, explicit LMB choices, projector
  overrides, zero-colour filtering and grouped comparison frames.
- Confirm rejected selections never construct their local numerators and that
  routed output has no local/projector construction. Verify one caller-thread
  execution path remains serial and ordered.
- Update Python wrappers, module registration, generated stubs, persisted
  request/state documentation and gammabench adapters in the corresponding
  change; do not add compatibility aliases.
- Run the approved Classes gate on reference and tadpole/self-energy/vacuum
  workloads: Auto/None Ir at most 1.005 in generation/full and allocation bytes
  no more than 5 percent above None. The heavy tier runs in both modes. Reuse A6's
  completed physical-frame proof instead of introducing another canonical pass.

== Sequencing

This preparation changes no API and does not race active fix/review work. Begin
implementation after integration fixes the exact baseline. Keep names,
route/numerators, and the Classes gate as separate described changes, each with
its own implementation, meaningful checks, adversarial review and fix pass.
Refresh developer-note digests against the resulting source, and record
instruction/allocation evidence with gammabench. M1 remains local evidence;
M2 is recorded at the eventually published commit.
