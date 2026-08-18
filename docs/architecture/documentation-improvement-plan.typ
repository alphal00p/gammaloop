= Documentation quality review and improvement plan

#quote(block: true)[
  *Status:* In progress

  *Review date:* 2026-08-18

  *Reviewed baseline:* `c9f4e32acd2c`

  *Scope:* The public portal, all five product manuals, generated API
  references, developer notes, authoring workflow, validation, and Pages
  publication.
]

== Executive assessment

The documentation delivery system is substantially stronger than the content
programme built on top of it. It has a typed product registry, deterministic
generation, immutable release snapshots, generated Rust and Python references,
compiled examples, local-link validation, search, and a shared responsive
shell. Those are valuable foundations and should be retained.

The principal weakness is usefulness at the user boundary. The generated
Python and CLI inventories are valuable correctness artifacts, but generation
has been treated too much as the end of the documentation job: long, flat
pages expose hundreds of entries with weak hierarchy, incomplete descriptions,
and presentation that is difficult to scan. At least one first-user tutorial
is not runnable in the environment it prescribes, and the uniform five-page
manuals stop before several products' headline workflows. Substantial current
notes and examples also remain outside the rendered documentation system.

The right next step is not a new documentation framework. It is to make the
existing system clear, task-complete, integrated, owned, and continuously
verified. Generated provenance is a strength, not an excuse for a lower bar:
Python, CLI, settings, and other generated reference pages should be as
readable, navigable, and thorough as authored pages.

This review uses the following priorities:

- *P0 — integrity:* a defect can publish stale material, produce a broken
  artifact, or give a new user instructions that fail.
- *P1 — user value and maintainability:* the information is incomplete,
  difficult to find, or likely to drift.
- *P2 — quality and reach:* the information exists, but presentation,
  accessibility, scientific context, or discovery should improve.

== What is already strong

The improvement work should preserve these properties:

- `docs/products/registry.toml` gives five products one explicit content and
  component model. The builder rejects missing sources, invalid routes, and an
  incomplete minimum manual shape.
- Rust catalogs and Python stubs are generated from implementation-owned
  metadata and compared with checked-in snapshots. GammaLoop has deep runtime
  signature/default parity, and Linnet has a real import-surface check.
- Authored Rust examples compile, shell examples are parsed, Python examples
  compile, data blocks are parsed, and selected side-effect-free catalog Rust
  examples execute.
- Full-site builds validate local links, fragments, and search-index targets.
  Nix builds the site twice and compares the results.
- Tagged snapshots are immutable, contain source and component provenance, and
  coexist with the moving `latest` channel.
- Developer notes under `docs/architecture/` must be classified in
  `docs/developers.toml`, and their rendered pages expose status, kind, pinned
  source links, search, and navigation.
- The generated shell already includes skip links, labelled navigation,
  breadcrumbs, theme controls, responsive layout, and readable no-script
  fallbacks for publication filtering.

== Review method and limits

The review inspected every registered product page, both content registries,
the documentation builder and example harness, the generated Python, CLI, and
settings output, the authoring commands, the Nix derivations, the Pages
workflow, and the persistent-history update script. The Vakint sibling-project
tutorial was also attempted from a clean temporary project, where its
documented dependency setup failed before compilation.

Every published developer note was also checked for lifecycle accuracy,
resolvable evidence, and agreement with its mapped implementation. That audit
found false present-tense claims in three pages labelled implemented, obsolete
fixtures and controls in performance records, an implemented drawing proposal
still presented as a future design, and a chronological build investigation
whose historical and current contracts were not clearly separated. This branch
corrects or explicitly archives the verified cases; the proposed freshness gate
is intended to prevent their recurrence.

This revision records the agreed publication model: `docs` is a temporary
branch used to prepare and review a pull request; after merge, `main` is the
sole authoritative source of mutable `latest` documentation. That is an
implementation constraint, not an open governance decision. Tagged snapshots
remain independent and immutable.

This was a source and build-contract review, not user research. It did not run
every licensed scientific backend, conduct a screen-reader session, or measure
search behavior with production analytics. Recommendations about common tasks
are inferred from each product's stated headline capabilities and should be
confirmed with maintainers and users before expanding the manuals. The absence
of an automated accessibility or browser gate is a verified process gap; it is
not, by itself, a claim that every rendered page currently violates a standard.

== Baseline

#table(
  columns: (1.5fr, 3fr, 3fr),
  table.header([*Surface*], [*Current baseline*], [*Consequence*]),
  [Product manuals],
  [5 products, 25 registered authored pages: one overview, one tutorial, one subject chapter, one interface chapter, and one release page per product],
  [Consistent, but too shallow for several headline capabilities],
  [Python reference],
  [5 generated module references with 95 exported classes/functions and downloadable stubs],
  [Inventory parity is strong, but long module pages, raw docstring presentation, implementation-detail noise, and missing member prose make the result hard to use],
  [GammaLoop CLI/settings reference],
  [285 commands, 412 arguments, and 166 settings on one generated page],
  [The data is exhaustive, but 132 commands, 117 arguments, and 163 settings lack explanatory text, and the flat presentation hides hierarchy and task context],
  [Developer area],
  [12 classified notes at the reviewed baseline; entries contain only free-form kind/status metadata],
  [Classification is strong, but “Implemented” has no owner, review date, verified code scope, expiry, or enforceable meaning],
  [Unintegrated corpus],
  [At the reviewed `c9f4e32acd2c` baseline, 4,203 lines under legacy `docs/*.md` and `docs/idenso/*` sat outside the product and developer registries],
  [Useful current guidance, examples, dated evidence, plans, and raw output were neither canonical Typst sources nor part of one navigable documentation system],
  [Manual examples],
  [Rust compiles; shell uses `bash -n`; Python uses `compile()`; only selected catalog Rust examples execute],
  [Syntax errors are caught, but most promised outcomes are not],
  [Release guidance],
  [Every product family records at least one missing or fragmented changelog surface],
  [Users cannot reliably derive upgrade risk from one canonical place],
  [Link checks],
  [Full-site local HTML links/fragments/search entries are checked, excluding Rustdoc and external URLs],
  [Important broken targets can still pass CI],
  [Publication policy],
  [`docs` prepares a reviewed change; `main` is authoritative after merge],
  [Preview and validation should happen before merge; only the merged main revision should promote `latest`],
)

The baseline is a snapshot, not a permanent score. Counts should be generated
by the quality job once the proposed inventory exists.

== Implementation progress

The first implementation slice is
#link("https://github.com/alphal00p/gammaloop/commit/6feaba98e48a2cedb4688d3d0cc9a214716a8677")[`6feaba98e48a`]
(2026-08-18). It delivers foundations and selected content improvements; it
does not close the full acceptance criteria in this plan.

The follow-up implementation is tracked in
#link("https://github.com/alphal00p/gammaloop/pull/96")[pull request 96]. The
table below describes that branch, not the older production site. It is kept
deliberately explicit about acceptance debt: a generated inventory or a new
page counts as progress, but not as completion when its browser, execution, or
coverage gate is still absent.

#table(
  columns: (1.5fr, 3fr, 3fr),
  table.header([*Area*], [*Implemented on the documentation PR*], [*Still open*]),
  [Python and CLI reference],
  [Exact Clap actions, arity, aliases, conflicts, globals, and usage; compact family and symbol indexes; nested disclosures; stable deep links; responsive setting cards; structured Python docstrings; and a generation gate that now rejects empty descriptions for all 79 public commands, 412 public arguments, and 166 settings],
  [The non-empty gate does not yet measure substantive prose or task backlinks. Python coverage remains uneven outside the GammaLoop vertical slice, and responsive/light/dark/accessibility browser fixtures are not yet enforced],
  [Manuals and examples],
  [Every product overview now offers a three-way task chooser. The GammaLoop process-generation, Linnet algorithms, Spenso networks, Idenso algebra, and Vakint evaluation chapters now include realistic source/configuration examples, expected invariants, failure interpretation, and targeted next references. A structurally validated registry records one beginner and one real-value canonical journey per product with its source, route, prerequisites, expected invariant, verification tier, owner group, and test command],
  [Only three authored examples currently execute; clean external beginner journeys, remaining repository/scientific examples, named individual ownership, and the scheduled licensed/heavy tier remain open],
  [Integrated corpus],
  [Legacy current material is classified and routed, authored prose is Typst-only apart from exact `README.md` and `AGENTS.md` compatibility files, companion package changelogs now have manual routes, and the portal routes calculations, graphs, tensors, identities, and integrals to the appropriate product with accurate ecosystem roles],
  [Version-level headings included from companion changelogs are not yet individually indexed, and release coverage remains incomplete],
  [Developer architecture],
  [The reviewed current notes were corrected against implementation and now carry lifecycle, review, trigger, and optional verified-scope metadata],
  [Named ownership, verified scopes for eight of the ten current records, trigger-aware owner review, and review automation remain incomplete],
  [Cross-run build reuse],
  [The Pages derivation is split into a Cargo producer and documentation consumer. The cold PR run rebuilt the producer and completed the five-site step in 1,005 seconds. A following docs-only run restored the producer, and its retry found the producer already present and completed the five-site step in 6 seconds],
  [Collect at least five comparable canonical warm runs before enforcing the proposed percentile objective],
)

The cold-run evidence is
#link("https://github.com/alphal00p/gammaloop/actions/runs/32131249528")[Documentation Pages run 32131249528].
The quality summary classified the producer as `rebuilt`; this avoids treating
a successful cold build as proof of cross-run reuse.

The cross-run evidence is
#link("https://github.com/alphal00p/gammaloop/actions/runs/32134216164")[Documentation Pages run 32134216164].
Its first attempt restored the reusable producer; the retry classified it as
`already-present` and measured the five-site consumer at 6 seconds. The two
attempts are reported separately because queueing, downloading, and a first
consumer build must not be disguised as compilation time.

The expanded source-owned reference and manual implementation at revision
`9b4525426311` provides a second honest cold point:
#link("https://github.com/alphal00p/gammaloop/actions/runs/32144386817")[run 32144386817]
classified the Cargo producer as `rebuilt` and completed the five-site step in
1,337 seconds. The same run passed the deterministic semantic/no-JavaScript
audit on 19 representative pages and retained the complete preview artifact.
This longer result is recorded rather than averaged into the warm evidence;
the follow-up content-only run is the cross-run reuse comparison.

The public site is not a preview of this branch. At the review date its latest
deployment was
#link("https://github.com/alphal00p/gammaloop/actions/runs/32073830002")[run 32073830002]
from revision `f46892d216f1`; the pull-request workflow correctly retained the
newer complete preview as an artifact and skipped publication. The visible
`latest` site changes only after the reviewed branch is merged to `main`.

== Acceptance status at this revision

#table(
  columns: (auto, auto, 4fr),
  table.header([*ID*], [*State*], [*Evidence or remaining boundary*]),
  [`DOC-001`], [Partial], [Compiled usage semantics and representative parser checks exist; the exhaustive invocation matrix remains open],
  [`DOC-002`], [Partial], [The misleading beginner promises were corrected, but five clean external runs and the scheduled tier are not complete],
  [`DOC-003`], [Complete], [Clean single-product previews rewrite ecosystem links and pass local-link validation],
  [`DOC-004`], [Partial], [Shared compact reference components are implemented; browser visual and accessibility fixtures remain open],
  [`DOC-005`], [Partial], [GammaLoop's principal Python workflow is documented from implementation-owned docstrings; all-module member coverage remains open],
  [`DOC-006`], [Partial], [Catalog generation rejects empty public CLI, argument, or settings prose; current non-empty coverage is 79/79, 412/412, and 166/166, but substantive-quality and task-backlink gates remain open],
  [`DOC-007`], [Partial], [Five task choosers and five deeper subject chapters are present; checked beginner and real-value paths do not yet cover every headline capability],
  [`DOC-008`], [Partial], [Ten canonical product journeys now have validated owner-group, route, tier, invariant, and command metadata; remaining maintained repository/scientific examples and scheduled verification are not yet registered],
  [`DOC-009`], [Partial], [Canonical prose is Typst-only and the declared legacy corpus is routed; a generated zero-orphan content graph remains open],
  [`DOC-010`], [Partial], [Lifecycle, review triggers, and selected digests are checked and changed selected scopes fail validation; named owners, scopes for eight current records, and trigger-aware owner review require maintainer assignment],
  [`DOC-011`], [Partial], [Companion release routes exist, but component-to-release coverage is not complete],
  [`DOC-012`], [Partial], [The portal offers five task-first product routes and each overview has task-to-guide/reference links; a registry-owned shared task model and complete reference backlinks remain open],
  [`DOC-013`], [Complete], [Pull requests retain complete previews without publishing, while `main` and immutable tags own publication],
  [`DOC-014`], [Partial], [Local links, fragments, containment, source paths, and search targets are checked; scheduled external checking remains open],
  [`DOC-015`], [Partial], [A required deterministic audit checks semantics and no-JavaScript navigation on 19 representative generated pages; standards HTML validation, browser/axe, spelling, keyboard, viewport, and visual gates remain open],
  [`DOC-016`], [Partial], [Search is federated across five products and developer notes with deterministic client scoring, and product/developer pages have pinned source and issue routes; top-five ranking execution and feedback/search on every portal and Rustdoc page remain open],
  [`DOC-017`], [Partial], [The source-backed Idenso convention specification is deep; equivalent auditable anchors remain uneven across products],
)

== Developer architecture audit at this baseline

Recent file dates did not guarantee correctness. The audit checked present-tense
claims against code rather than treating the publication commit as proof:

#table(
  columns: (1.6fr, 3fr, 2.8fr),
  table.header([*Area*], [*Verified drift*], [*Disposition in this branch*]),
  [GammaLoop current architecture],
  [Repository scope was overstated; unmanifested state loading, external evaluator activation, additional weights, integration-workspace persistence, and the built-in quantity list had drifted],
  [Corrected the current architecture against the state, evaluator, event, integration, and quantity implementations],
  [Spenso parsing diagram],
  [The `ShorthandParsing` shape, strict tensor filtering, opaque-function boundary, compact-vector chain endpoints, and trace closure no longer matched the parser/materializer],
  [Updated the diagram and its examples to the current parser contract],
  [Schoonschip network diagram],
  [`schoonschip_net` was said to run a final pattern cleanup that the current wrapper no longer performs],
  [Corrected the entry-point and comparison text],
  [Drawing migration],
  [A completed and subsequently evolved migration remained under proposals as a future callback sketch],
  [Reclassified it as an implemented record and pointed to the canonical generator, style module, and templates; retained the sketch as history],
  [Spenso performance records],
  [Two root fixtures, tests, environment controls, and unresolved-failure claims no longer exist or are contradicted by current regression coverage],
  [Marked the records historical/superseded and non-runnable instead of presenting their measurements as current],
  [Schoonschip investigation and archived baseline],
  [One removed setting and two removed trace controls were advertised; the archived baseline has a non-resolving revision and removed tests],
  [Corrected live controls and labelled the baseline non-reproducible],
  [Nix/Crane record],
  [A 1,500-line chronology mixes historical “now” statements with the maintained cache contract],
  [Added an explicit reading boundary now; splitting the current contract from experiment history remains improvement work],
)

The UV-renormalization architecture still matched its implementation in the
scoped audit. The CFF surface-cache proposal also remains genuinely open:
`Graph::generate_cff` still owns and mutates the graph cache. These checks are
point-in-time evidence, not a substitute for the continuous contract below.

== Critical findings

=== P0. The generated CLI reference can emit false invocation syntax

GammaLoop's compiled Clap tree is a sound source of truth, but the exported
schema does not preserve enough of that truth. It records argument names,
defaults, and value names without the action and arity needed to distinguish a
flag from a value-taking option. The renderer consequently documents the valid
boolean flag `--clean-state` as `--clean-state <CLEAN_STATE>`. A generated page
that looks exhaustive while showing an invocation the parser rejects is more
dangerous than an obviously incomplete page.

The 285-node command export also mixes user operations with generated help and
choice nodes without explaining the distinction. Generation prevents manual
drift only when the generated representation and renderer preserve the
semantics of the compiled interface.

*Required change*

- Export exact Clap usage, action, arity, value requirements, conflicts, and
  inherited/global option semantics.
- Render flags, value-taking options, positionals, aliases, and generated help
  nodes according to those semantics.
- Parity-test documented invocations against the compiled parser rather than
  snapshot-testing JSON shape alone.
- Version the generated schema when adding the missing semantics.

*Complete when*

- no documented option falsely requires or omits a value;
- generated command examples and a representative matrix of flag, option,
  positional, alias, and nested-subcommand invocations parse successfully;
- rendered usage agrees with compiled `--help`; and
- generated helper/choice nodes are omitted from primary navigation or clearly
  labelled and grouped.

=== P0. First-user promises are not consistently executable

The manuals are careful and polished, but three visible promises do not match
the verified behavior:

- Vakint tells a user to create a sibling project and add
  `../gammaloop/crates/vakint`. In a clean sibling project, dependency
  resolution fails because the required Symbolica override is owned by the
  workspace root and is not inherited by a path consumer.
- Linnet's registry promises a first subgraph-aware operation and its tutorial
  says the builder returns node and half-edge identifiers. The example retains
  only node identifiers and calls `cycle_basis()` on the full graph.
- GammaLoop calls its tutorial “First calculation”, but the first-success path
  stops after process generation and explicitly defers integration.

These are more serious than missing prose: they teach users to distrust the
rest of the manual.

*Required change*

- Run every beginner tutorial verbatim in a clean temporary environment outside
  the workspace.
- Fix the dependency contract or documented installation path for Vakint.
- Either demonstrate an actual Linnet subgraph operation or change the promise
  to the operation shown.
- Add a bounded GammaLoop sample-evaluation result, or rename the tutorial to
  “Create your first state”.
- Give every procedural example prerequisites, expected output or invariant,
  expected runtime/resource class, failure modes, and a next step.

*Complete when*

- all five beginner journeys pass in CI from empty temporary directories;
- every fence is declared `run`, `compile`, or `syntax` and CI rejects an
  unclassified fence; and
- heavy, licensed, or external-tool examples are visibly classified and run in
  a provisioned scheduled tier rather than silently skipped.

=== P0. The advertised partial-build contract can emit broken navigation

`just docs-site PRODUCT` selects and builds one product, then renders a portal
containing links to all five products. Generated-link validation runs only when
the selected product is `all`. A clean single-product output can therefore
finish successfully with four sets of broken portal links.

*Required change*

Define one honest partial-build contract:

- render a product-scoped preview without the multi-product portal;
- emit safe links to the canonical deployed products; or
- build the minimal redirect/stub set needed for a complete local portal.

*Complete when*

- every documented `docs-site` form is tested from an empty output directory;
- local link and fragment validation always runs for the artifact it emits; and
- a single-product build contains zero broken navigation or search targets.

=== P1. Generated Python and CLI data is not yet a good reference experience

The reference pipeline has strong implementation provenance, but the rendered
result is mechanically broad and editorially shallow. Autogeneration is a
sourcing mechanism; it does not reduce the need for hierarchy, explanation,
examples, or visual design.

The Python reference currently exposes 95 top-level exports across five
modules. Each module is rendered as one long sequence of generic cards, with
no local symbol index or filter and no per-symbol route. NumPy-style docstrings
are escaped and split into paragraphs instead of being parsed into semantic
Parameters, Returns, Raises, Notes, and Examples sections. Method parameters
become recursively nested member blocks, internal enum labels such as
`PythonClass` and `Getter` leak into the UI, and internal build features receive
user-facing badges. All 95 automatically attached examples amount to importing
the symbol and calling `help(...)`, rather than demonstrating a task.

Completeness varies sharply. GammaLoop renders 40 exports, but only four are in
the curated supported surface; 36 implementation-detail entries receive almost
equal visual and search prominence. All 360 rendered GammaLoop members lack
member documentation. Linnet has 189 undocumented members out of 313, while
Spenso has 69 out of 244. Inventory parity is therefore not yet documentation
parity.

The CLI/settings page has a related shape problem. It puts 285 commands, 412
arguments, and 166 settings on one route. Beyond the truthfulness issue above,
132 command nodes and 117 arguments have no explanatory text, while only three
settings have descriptions. Full command paths and a flat settings table expose
the data, but do not help a user understand command families, configuration
groups, precedence, or which settings matter to a particular task.

*Required change — Python*

- Build an information hierarchy of module, class/function, member, and
  parameter with stable URLs or stable indexed routes, breadcrumbs, a local
  symbol index, kind filters, and search.
- Parse supported docstring structure into semantic sections and human labels;
  render normalized readable signatures while retaining the exact signature
  and downloadable stub.
- Put the curated supported surface first. Collapse, separate, and de-rank
  implementation details unless users need them.
- Define a completeness contract for supported exports and public members:
  substantive summary, signature, parameters/defaults, return behavior,
  raises/failure modes where applicable, provenance, and a guide or meaningful
  tested example.
- Link source provenance to the binding implementation where possible, not
  only to the generated `.pyi` file.

*Required change — CLI and settings*

- Group commands by real command hierarchy and settings by conceptual path;
  provide compact indexes, filtering, and focused detail routes or sections.
- Show exact invocation, aliases, inherited options, defaults, allowed values,
  conflicts, environment/config precedence, and copyable examples.
- Require substantive descriptions for public commands, arguments, and
  settings; generated placeholders and empty cells do not count.
- Link each primary command family and settings group to the manual tasks that
  use it, and link those guides back to exact reference anchors.

*Required change — presentation*

- Use the shared manual shell and design tokens for every reference route, with
  readable long signatures, responsive parameter/settings layouts, useful
  local navigation, and restrained visual density.
- Add light/dark visual regression coverage and browser/accessibility checks at
  375 px, 768 px, and a desktop width.

*Complete when*

- a user can reach a known supported Python symbol, CLI command, or setting in
  at most two navigation interactions from its reference hub and in the first
  five search results;
- 100% of supported Python exports and public members meet the completeness
  contract or are explicitly marked not applicable;
- 100% of public commands, arguments, and settings have truthful syntax and
  substantive descriptions;
- no raw docstring markup, schema debug label, or internal-only feature name is
  exposed as primary user documentation;
- each module and primary command family links to meaningful, tested examples
  for its principal tasks; and
- representative pages have no viewport-level horizontal scroll at 375 px and
  no serious or critical accessibility findings.

=== P1. The manuals cover a template, not each product's real task surface

The five manuals each have exactly five short authored pages. Their common
shape is useful, but one tutorial and one subject chapter cannot carry all of
the distinct workflows:

#table(
  columns: (1.2fr, 2.2fr, 3.5fr),
  table.header([*Product*], [*Strong current entry point*], [*Highest-value missing worked path*]),
  [GammaLoop], [Stateful CLI generation and resume], [A bounded calculation, sample evaluation, observables, output interpretation, precision/stability, and validation],
  [Linnet], [Half-edge construction and cycle basis], [Subgraph traversal, cuts, mutation mappings, and one drawing workflow],
  [Spenso], [Pairwise tensor contraction], [Network construction, planning, execution, and inspection],
  [Idenso], [One Python metric contraction], [Dummy-index handling, reversible cooking, Dirac/color passes, and a Rust path],
  [Vakint], [Dependency-free topology matching], [One analytic or numerical evaluation with normalization and expected Laurent coefficients],
)

API hubs are broad, but tutorials rarely deep-link the generated item anchors
for the symbols they use. Some Python snippets are explanatory fragments with
undefined inputs while looking like runnable examples. Conceptual physics
claims also need stronger anchors: Local Unitarity, normalization/sign
conventions, and precision-sensitive backend behavior should point to defining
references and regression-backed examples.

*Required change*

- Replace the five-page minimum as a content goal with a per-product task map
  covering the three principal user goals, prerequisites, beginner journeys,
  real-value workflows, troubleshooting, and result interpretation.
- Build vertical journeys: a manual task should lead to a canonical executable
  example, explain the important choices, link every used Python/Rust/CLI item
  to its exact reference entry, and show the expected result or invariant.
- Write shared cross-product concepts once and link them contextually from each
  product rather than copying abbreviated explanations into five manuals.
- Treat operational details—precision, resource expectations, persistence,
  diagnostics, failure recovery, and upgrade effects—as part of the workflow,
  not optional appendices.

*Complete when*

- each product has a maintained task matrix derived from its three most common
  user goals;
- each headline capability has at least one tested beginner path and one tested
  real-value path;
- every runnable snippet is self-contained, and every fragment is labelled as
  such;
- tutorial symbols deep-link to the generated anchored API entry; and
- convention-sensitive scientific claims include a citation, an explicit
  convention or equation, and a checked example.

=== P1. Examples are shallow, scattered, and weakly connected to the manuals

The example harness catches useful classes of error, but it primarily validates
syntax and compilation. Only selected catalog Rust examples execute, and many
manual snippets are fragments rather than complete programs. Repository tests,
examples, command cards, and scientific reference cases can contain the real
workflow and expected behavior while the public manual shows a smaller,
independently maintained approximation.

There is no common inventory recording what task an example teaches, whether it
is canonical, its prerequisites, expected output, runtime and infrastructure
tier, owner, test command, or public documentation route. Consequently an
example can remain green without being discoverable, or remain visible while
no longer demonstrating a meaningful result.

*Required change*

- Create one example registry covering manual fences, repository examples,
  command cards, reference cases, and any maintained notebooks. Record product,
  task, audience, prerequisites, canonical inputs, expected result/invariant,
  runtime tier, owner, test command, and documentation route.
- Give each maintained example one canonical source. Include or extract it into
  manuals and reference pages so copied snippets cannot drift independently.
- Provide both a fast beginner example and a representative real-value example
  for each principal product task, including output interpretation and common
  failure modes.
- Execute hermetic examples on pull requests. Run licensed, networked, or heavy
  examples in an explicit scheduled tier and publish their last verified
  revision and result.

*Complete when*

- every maintained example has exactly one canonical documentation route and
  verification command;
- the three principal tasks for each product each have a beginner and a
  real-value example;
- every displayed result is checked by value, tolerance, snapshot, or explicit
  invariant rather than compilation alone;
- manual snippets and repository examples cannot silently diverge; and
- stale scheduled examples are visibly marked instead of looking current.

=== P1. Release and freshness information is visible but not trustworthy

The manuals honestly disclose their debt: Spenso 0.6.0 notes stop at 0.5.6,
Idenso 0.3.0 notes stop at 0.2.5, Vakint has no standalone changelog,
GammaLoop has no consolidated changelog, and Linnet's Python distribution has
no separate history. Several prose pages also repeat manifest versions or
external-tool minima that the builder already knows from generated metadata.

The publication cache has a resilience fallback, but the rendered portal does
not expose a stale-data state and the builder does not enforce cache age.

*Required change*

- Render package versions and generated dependency minima from their canonical
  metadata rather than handwritten literals.
- Establish one release-history location for every user-installable component.
- Compare each manifest version with the latest documented release entry.
- Add `reviewed_at` to current manuals and expose stale publication data rather
  than silently presenting it as current.

*Complete when*

- no current component version exceeds its documented release coverage;
- generated version claims have no handwritten duplicate;
- current user pages have an owner and a review no older than the agreed
  service level; and
- publication data older than 14 days is visibly marked and causes a scheduled
  health failure or issue while offline builds remain reproducible.

=== P1. The documentation corpus is fragmented instead of integrated

Classification is enforced only for files directly under
`docs/architecture/`. At least 4,203 lines elsewhere under `docs/` include
current logging and drawing guidance, Idenso references, living plans, dated
CI/benchmark evidence, rebase notes, and a raw terminal dump. Some are linked
only as GitHub source; most have no public route, audience, owner, status, or
supersession relationship. Classification alone would not solve this: useful
material must be migrated into the manuals, shared guides, reference, and
search rather than merely labelled in place.

The configured ecosystem model is also underused. `pillar`, `affiliation`,
`featured`, and product `related` metadata is validated or deserialized but is
not turned into task-oriented navigation. All five portal cards are called a
“Research project” even though one is the application and the others are
specialized libraries or layers. Product manuals link registered developer
notes to raw GitHub blobs instead of their rendered, searchable pages.

*Required change*

- Make Typst the canonical source format for all maintained authored
  documentation. Do not create Markdown or HTML mirrors: rendered HTML/PDF is
  build output, not a second editable source. Files named exactly `AGENTS.md`
  and `README.md` remain the only compatibility exceptions because external
  tooling discovers those names.
- Treat every other existing Markdown or hand-authored HTML document as
  explicit migration debt. New registry entries must use `.typ`; a substantive
  revision to a legacy source must migrate it to Typst, update its routes and
  links, and remove the old source in the same change.
- Establish one content graph with multiple views. Each concept or example has
  one canonical source; product manuals, generated reference, developer notes,
  examples, navigation, and search provide contextual routes to it.
- Classify every authored source as product documentation, current developer
  guidance, proposal, investigation, archive, source-only reference, fixture,
  or draft.
- Extend registry metadata with a constrained lifecycle status, named owner,
  `reviewed_at`, applicable products/topics, and optional
  `supersedes`/`superseded_by`.
- Publish current logging and diagnostics guidance and link it from the
  GammaLoop CLI/API and developer material. Publish drawing guidance with clear
  GammaLoop/Linnet/Linnest ownership. Migrate `docs/idenso/*` into Idenso manual
  reference or appendices with Spenso syntax cross-links.
- Turn living plans into lifecycle-labelled proposals or current-status pages.
  Move CI summaries, benchmarks, raw terminal output, and rebase notes to dated
  investigation/archive routes with a current summary or explicit supersession
  link; preserve the evidence rather than deleting it.
- Use product relationships and ecosystem roles to let users choose by task.
- Link registered developer notes to their rendered routes, keeping the pinned
  source link as secondary provenance.
- Keep README files as contributor entry points that link to rendered canonical
  guidance rather than acting as parallel manuals.

*Complete when*

- 100% of meaningful authored documentation has one canonical Typst source and
  a rendered route; only fixtures and raw evidence may be explicitly
  source-only;
- no new Markdown/HTML source or parallel format can enter the content graph,
  and the legacy-source count decreases monotonically to zero;
- CI fails an orphan-document fixture;
- all published pages have constrained status, owner, review date, and scope;
- superseded material cannot look current in navigation or search; and
- every current page is in navigation or contextually linked from a maintained
  guide, and a generated link-graph report has zero orphan current pages;
- no current user guidance exists only as a GitHub blob; and
- a new user can choose “calculations”, “graphs”, “tensors”, “identities”, or
  “integrals” from the portal and reach the correct first task.

=== P1. “Current” developer architecture has no freshness contract

`docs/developers.toml` schema 1 gives a note an ID, title, summary, source,
kind, and free-form status. The builder checks that those strings are non-empty,
that the source exists, and that every top-level architecture file is
classified. It does not know who owns a claim, when it was reviewed, which code
proves it, whether that code changed, or whether a proposal was implemented or
superseded. A note can therefore say “Implemented” forever while remaining
fully green.

The rendered source link makes this worse by pinning the note and its
hand-authored code links to the current Pages build SHA. That proves which prose
was published, not that the architecture was verified at that revision. Fixed
line anchors in the legacy HTML diagrams can move to unrelated code while their paths
continue to pass link validation. The concrete drift corrected by this audit
shows that recent modification dates and working links are insufficient.

*Required change*

- Introduce a developer-registry schema with a constrained lifecycle separate
  from freshness: `current`, `proposal`, `investigation`, `archived`, or
  `superseded`. Drafts remain unpublished. An implemented proposal is not a
  lifecycle of its own: migrate it to `current` when the page is maintained as
  the active contract, or to `archived` when it is only a design record, and
  require `implemented_by` evidence in either case.
- Migrate the existing free-form values explicitly:
  #table(
    columns: (2fr, 4fr),
    table.header([*Existing status*], [*Schema-2 lifecycle*]),
    [`Implemented`], [`current`],
    [`Implemented record`], [`current` if it retains an active contract; otherwise `archived`, with `implemented_by`],
    [`Proposal`], [`proposal`],
    [`Investigation record`], [`investigation`],
    [`Superseded investigation`], [`superseded`, with `superseded_by`],
    [`Historical experiment` and `Archived · non-reproducible`], [`archived`],
    [`Living engineering record`], [split into a `current` contract and `archived` evidence before migration],
  )

- Resolve owners through one owner table. Require `reviewed_at`, a review/PR
  record, applicable products/topics, stable code paths and symbols, and broad
  review-trigger globs for current and proposal notes.
- Resolve stable Rust symbols to their current declarations and hash the
  canonical verified scopes. A changed or missing digest is a required review,
  not an automatic claim that prose is still correct. Do not use stored line
  numbers or a self-referential HEAD SHA as semantic evidence.
- Give investigations and archives an immutable evidence revision, capture
  date, command, fixtures, toolchain/environment, and expected artifact or
  result. They remain valid historical records when HEAD changes, but must not
  look current.
- Require `superseded_by` for superseded material, validate that lifecycle
  graph for missing targets/cycles, and remove superseded pages from primary
  navigation and search while preserving redirects and provenance.
- Map changed code through review triggers and the reverse code-reference graph.
  Require the note owner to refresh its scope digest or explicitly attest that
  the change has no documentation impact. Protect the note and registry through
  CODEOWNERS.
- Render owner, lifecycle, review date, evidence/verified scope, and freshness
  state. Keep the build revision distinct from the review/evidence revision.
- Start with proposed defaults—warn current architecture at 60 days, block
  publication at 90 days unless reviewed, and require proposal disposition at
  180 days—then ratify or adjust them with the named owners before enforcement.
  Do not age-fail frozen investigation/archive evidence; fail it only for
  invalid provenance or misleading lifecycle presentation.
- Split living engineering logs into a short canonical current contract and
  dated experiment/incident records.

*Complete when*

- 100% of current developer notes have a resolving owner, review record,
  verified scopes/digests, review triggers, and a contextual inbound link beyond
  the developer index;
- no current code path or symbol reference is missing or ambiguous;
- no changed verified scope can publish until its owner reviews or explicitly
  attests the impact;
- no overdue current note can publish, while historical records remain visibly
  pinned to their evidence revision;
- every investigation's declared fixture exists and command parses at capture,
  or the page is explicitly marked non-reproducible; and
- every implemented/superseded proposal names the implementing or replacement
  evidence.

=== P1. Production is doing work that belongs in a review gate

The accepted delivery flow is now clear: `docs` is temporary PR staging and
`main` is authoritative after merge. The automation should encode that flow.
Updates to `docs` should build the complete preview and quality report through
its pull request without promoting `latest`; the merged `main` revision should
rebuild and deploy the reviewed sources. This is a workflow-enforcement task,
not a choice between two long-lived documentation sources.

The Pages workflow currently has no pull-request trigger, and its validation
set differs from the full unit-test group for the six documentation crates. A
documentation edit can therefore reach the deployment workflow before the
reviewer has a complete rendered artifact and quality report.

Link validation has additional blind spots:

- Rustdoc trees and external URLs are skipped;
- repository `source-link` targets are constructed without checking that their
  paths exist;
- local path normalization does not explicitly require the result to stay under
  the generated output, allowing a false pass against a host file; and
- HTML validity, spelling, browser behavior, and accessibility are not checked
  by standards-based tools.

*Required change*

- Add one required, path-routed “Documentation quality” pull-request check.
- Build and retain the complete rendered preview artifact on pull requests,
  including the temporary `docs` staging pull request, but do not let those
  revisions mutate production `latest`.
- Promote `latest` only from merged `main`; keep scheduled mutable builds on
  `main` and preserve tag-triggered immutable snapshots independently.
- Run content checks for content-only changes and the full six-crate docs test
  group for renderer, schema, exporter, script, or workflow changes.
- Keep local checks deterministic; run flaky external HTTP checks separately on
  a schedule with retries and an allowlist.

*Complete when*

- no production deployment is the first full validation of its commit;
- an update to the `docs` staging pull request produces a preview but cannot
  change `latest`, while the merged `main` commit can publish it;
- local links, fragments, Rustdoc, repository source paths, search entries,
  redirects, and asset references have zero failures;
- link resolution cannot escape the generated root;
- every preview is retained long enough for review; and
- the job summary reports example modes, page counts, link counts, cache reuse,
  and step timings.

=== P2. Accessibility, progressive enhancement, search, and feedback are unmeasured

The shell has good semantic foundations, but mobile manual navigation is moved
off-screen and depends on JavaScript to open. Search is scoped separately to a
single product or the developer area, indexes headings and API entries rather
than full prose, uses unranked substring matching, and caps results. The portal
has no federated search. There is also no direct “report a documentation issue”
path, CODEOWNERS file, HTML validator, browser smoke test, accessibility test,
canonical metadata, or sitemap.

*Complete when*

- all routes remain navigable at 375 px with JavaScript disabled;
- keyboard focus, search, menu, theme, and dialog behavior pass browser smoke
  tests;
- HTML validation has zero errors and axe reports zero serious or critical
  findings on the portal, a tutorial, a manual chapter, a reference page, and a
  developer note;
- representative task queries return a current relevant page in the first five
  federated results;
- every page offers a revision-pinned edit/source link and a pre-labelled issue
  path; and
- canonical URLs and a sitemap distinguish mutable latest pages from immutable
  snapshots.

== Target documentation model

The target is one content graph presented through five user-facing views:

+ *Portal:* choose a task or product; find people, publications, citations,
   and developer material.
+ *Tutorial:* one maintained path from an empty environment to a verified
   result.
+ *How-to and explanation:* task recipes and scientific/software concepts,
   linked rather than duplicated across products.
+ *Reference:* generated, version-specific Rust, Python, CLI, settings, and
   topology facts.
+ *Developer record:* implemented architecture with verified code scopes and
   review triggers, plus proposals, investigations, and archives with explicit
   lifecycle and evidence provenance.

Compatibility READMEs should orient contributors and link to the rendered
material, not act as second user manuals. Generated reference is a first-class product
surface paired with task documentation; generation is a sourcing mechanism,
not an exemption from editorial or visual quality. Dated evidence should be
preserved, but never presented as the current contract without an explicit
current summary.

The views must be bidirectionally connected. A manual task links to the exact
symbols, commands, settings, concepts, and canonical examples it uses. A
reference item links back to the guides and examples that give it purpose.
Shared concepts have one maintained explanation, while contextual links make
them visible from every relevant product.

== Reference-quality design requirements

=== Python

Keep the runtime inventory and checked stubs as the source of signatures and
availability. Add editorial information at the implementation boundary through
binding docstrings, structured annotations, and validated grouping metadata;
do not create a second handwritten signature catalog.

The renderer should provide:

- a module index followed by stable class/function and member routes or indexed
  sections;
- fully qualified names, human-readable kinds, normalized and exact
  signatures, source provenance, availability, and version information;
- semantic Parameters, Attributes, Returns, Raises, Notes, and Examples
  sections rather than escaped raw docstring layout;
- a curated supported surface in primary navigation, with implementation
  details separated, collapsed, and de-ranked in search;
- local symbol filtering, breadcrumbs, a useful table of contents, copyable
  signatures/examples, and responsive tables; and
- guide/example backlinks for the principal tasks of every module.

Completeness reporting must distinguish inventory parity from documentation
parity. An export can be present in the stub and still fail the documentation
gate because its members, parameters, return semantics, or examples are empty.

=== CLI and settings

Keep compiled Clap and settings schemas authoritative, but extend the export so
the renderer does not have to infer semantics. The reference model should
preserve command hierarchy, exact usage, action, arity, aliases,
global/inherited options, conflicts, defaults, allowed values, configuration
precedence, and source provenance.

The top-level page should explain command families and configuration groups,
not print the full database. Focused detail routes or sections should show
copyable invocations, expected effects, persistence behavior, related settings,
failure modes, and links to task guides. A filter should search command paths,
aliases, arguments, and setting paths without adding hundreds of entries to the
page table of contents.

=== Shared visual contract

All generated reference uses the same shell, tokens, navigation behaviors,
feedback routes, and accessibility targets as the manuals. Reference-specific
components should be designed for dense technical data without becoming a wall
of equally weighted cards: compact indexes, progressive disclosure, readable
line wrapping, responsive field layouts, clear status badges, and restrained
use of borders and pills.

Visual regression fixtures should cover a module index, a large Python class,
a function with rich docstrings, a command family, and a settings group in
light and dark modes at 375 px, 768 px, and desktop width.

== Prioritized work register

#table(
  columns: (auto, auto, 3fr, 3fr),
  table.header([*ID*], [*Priority*], [*Deliverable*], [*Verification*]),
  [`DOC-001`], [P0], [Make generated CLI syntax semantically faithful to compiled Clap], [Zero false syntaxes; rendered invocation matrix parses and matches `--help`],
  [`DOC-002`], [P0], [Repair or rename the Vakint, Linnet, and GammaLoop beginner promises], [Five clean, verbatim tutorial runs pass],
  [`DOC-003`], [P0], [Define the single-product preview contract], [Every clean partial build has zero broken links],
  [`DOC-004`], [P1], [Build the shared reference information architecture and visual components], [Python/CLI fixtures pass navigation, responsive, dark/light, and accessibility checks],
  [`DOC-005`], [P1], [Make Python reference complete, structured, and task-linked], [Runtime/export/member/signature/default parity plus substantive supported-symbol documentation reaches 100%],
  [`DOC-006`], [P1], [Make CLI and settings reference complete, grouped, and task-linked], [100% of public commands, arguments, and settings are truthful and substantively described],
  [`DOC-007`], [P1], [Create a top-task matrix and deepen each product manual], [One beginner and one real-value checked path per headline capability],
  [`DOC-008`], [P1], [Register canonical examples and execute them by verification tier], [Every maintained example has one route, owner, expected result, and test command],
  [`DOC-009`], [P1], [Integrate the full corpus with lifecycle metadata, canonical Typst sources, and a generated link graph], [100% of meaningful sources have one Typst source and a canonical route; zero orphan current pages],
  [`DOC-010`], [P1], [Enforce developer-architecture freshness and verified code scopes], [Every current note is owned, recently reviewed, scope-clean, symbol-valid, and contextually linked; changed scopes block publication pending owner review],
  [`DOC-011`], [P1], [Make release/version/tool claims generated and complete], [Manifest-to-release coverage is current for every component],
  [`DOC-012`], [P1], [Use ecosystem metadata for task navigation and bidirectional manual/reference links], [Task chooser and contextual cross-link fixtures pass],
  [`DOC-013`], [P1], [Add the PR quality gate and enforce docs-preview/main-promotion flow], [The temporary `docs` pull request retains a complete preview but cannot mutate `latest`; merged `main` deploys],
  [`DOC-014`], [P1], [Complete link, source-path, containment, and redirect validation], [Zero internal failures, including Rustdoc; external report is scheduled],
  [`DOC-015`], [P2], [Add browser, HTML, accessibility, spelling, and no-JS checks], [Zero HTML errors and serious/critical axe findings on representative pages],
  [`DOC-016`], [P2], [Federate search and add edit/report feedback routes], [Task and symbol fixtures rank the correct page in the first five results],
  [`DOC-017`], [P2], [Add scientific citations and regression-backed conventions], [Each precision/sign/normalization-sensitive claim has an auditable anchor],
)

== Delivery sequence

=== Phase 0 — stop integrity regressions

Complete DOC-001 through DOC-003 first. They protect production, preview
artifacts, and first-user trust. Do not expand a generated reference whose
rendered invocation syntax can be false, or a content surface whose documented
local build can silently emit broken navigation. Enable the source-format guard
in this phase as an immediate invariant: CI must reject newly added authored
`.md`/`.html` files outside the two filename compatibility exceptions and reject
parallel editable copies. DOC-009 then migrates the declared legacy inventory.

=== Phase 1 — productize reference and establish one content graph

Complete DOC-004 through DOC-006, establish the inventory/schema portion of
DOC-009, establish the developer lifecycle/scope schema and initial audit under
DOC-010, and complete DOC-013. Prototype the shared reference design with one
large GammaLoop Python class, one command family, and one settings group before
mechanically applying it to all modules. Publish coverage and developer-scope
reports as warnings, fill the source-owned documentation gaps, then promote the
substantive coverage and freshness thresholds to required checks.

=== Phase 2 — make the five product journeys valuable

Complete DOC-007, DOC-008, DOC-012, and DOC-017 product by product. Prefer a
vertical slice—tutorial, canonical example, explanation, exact Python/CLI/Rust
deep links, troubleshooting, and result interpretation—for one product over
adding disconnected pages to all five. Keep licensed and expensive scientific
validation in a provisioned scheduled tier while preserving a hermetic
pull-request smoke path.

=== Phase 3 — migrate and connect the remaining corpus

Finish DOC-009 through DOC-011. Publish current guidance in manuals or shared
guides and move evidence into explicit investigation/archive routes. Apply the
developer freshness contract to every record, and split large living records
into a short current contract plus dated evidence. Add supersession links and
redirects before moving a source; do not silently delete historical
investigation material. Require the generated link graph to show how every
current page is reached.

=== Phase 4 — improve reach and continuous quality

Complete DOC-014 through DOC-016. Run deterministic checks on pull requests and
network/licensed/full-browser health checks on a schedule. Publish the scorecard
so degradation is visible without making every transient external failure block
an author.

== First vertical implementation slice

The first implementation pull requests should prove the model end to end:

+ Use the canonical `docs/architecture/architecture-current.typ` source,
   introduce the developer lifecycle/review schema, register stable code
   scopes, and turn the corrected state-loading,
   evaluator, event-weight, integration-workspace, and quantity claims into the
   first scope-digest fixtures.
+ Extend the CLI schema with action, arity, and exact usage; correct
   `--clean-state`; and add compiled-parser parity fixtures.
+ Build shared reference index/detail components and apply them to the
   GammaLoop Python module, one command family, and one settings group.
+ Replace raw Python docstring rendering with structured sections and add a
   coverage report that separates runtime parity from documentation parity.
+ Deepen one GammaLoop task from clean setup through a checked result, using a
   canonical example and bidirectional links to every referenced command,
   setting, and Python symbol.
+ Migrate the legacy `docs/logging.md` content into a canonical Typst diagnostics
   journey and use the
   content graph to prove it is discoverable from the manual, CLI reference,
   developer area, and search without duplicating its prose.

This slice exercises generated truth, visual design, substantive content,
examples, and integration before the same structure is replicated across five
products.

== Scorecard

The quality job should publish these measures on every canonical build:

#table(
  columns: (4fr, 1.4fr),
  align: (left, right),
  table.header([*Measure*], [*Target*]),
  [Non-main revisions able to mutate `latest`], [0],
  [False generated CLI invocation syntaxes], [0],
  [Public CLI commands, arguments, and settings with substantive descriptions], [100%],
  [Supported Python exports and public members meeting the documentation contract], [100%],
  [Generated `help(...)` placeholders counted as meaningful examples], [0],
  [Reference fixtures passing responsive light/dark visual checks], [100%],
  [Broken local links, fragments, source paths, redirects, or search targets], [0],
  [Clean beginner journeys executed], [5 / 5 products],
  [Principal product tasks with beginner and real-value checked examples], [100%],
  [Maintained examples with one canonical route and verification tier], [100%],
  [Meaningful authored sources with one canonical Typst source and rendered route], [100%],
  [Orphan current pages], [0],
  [Published pages classified, owned, scoped, and reviewed], [100%],
  [Current developer notes meeting the freshness/verified-scope contract], [100%],
  [Missing or ambiguous developer code references], [0],
  [Changed verified scopes awaiting owner disposition at publication], [0],
  [Overdue current developer reviews], [0],
  [Current developer notes reachable only from the developer index], [0],
  [Current components with release-coverage gaps], [0],
  [Runtime-capable Python surfaces with full parity], [100%],
  [Serious or critical axe findings on representative pages], [0],
  [Publication cache age], [at most 14 days, or visibly stale],
  [Representative task searches returning a correct result in the first five], [100%],
  [Warm docs-only builds that rebuild the reusable Cargo producer], [0],
)

For build performance, collect at least five canonical runs before enforcing a
percentile. The initial service objective is a warm p95 of at most 10 minutes
and a cold p95 of at most 25 minutes, with publication/deployment time reported
separately. A content-only edit should never run unrelated workspace tests.

== Decisions required before implementation

The publication source is no longer an open choice: `docs` is PR staging and
`main` alone is authoritative for mutable `latest`. Maintainers still need to
record two ownership and execution choices:

+ Who owns the portal/toolchain, generated Python and CLI reference quality,
   each of the five product manuals, and the GammaLoop runtime/UV,
   Spenso/Idenso parsing, build-system, and documentation-system architecture
   records? Owners must also approve the current/proposal review service levels.
+ Which scientific examples are hermetic pull-request checks, and which require
   licensed or scheduled infrastructure?

Once those decisions are recorded, the work register can be split into small
pull requests with one acceptance criterion per change. The reference prototype
and corpus inventory can begin while named owners and expensive-example tiers
are assigned; publication and product behavior changes still require their
ordinary review.
