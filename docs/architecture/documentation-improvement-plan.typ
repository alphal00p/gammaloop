= Documentation quality review and improvement plan

#quote(block: true)[
  *Status:* In progress

  *Review date:* 2026-08-19

  *Reviewed baseline:* `c9f4e32acd2c`

  *Scope:* The public portal, all five product manuals, generated API
  references, developer notes, authoring workflow, validation, and Pages
  publication.
]

== Executive assessment

The documentation delivery system is substantially stronger than the content
programme built on top of it. It has a typed product registry, deterministic
generation, immutable release snapshots, revision-specific native Rustdoc,
generated Python references, compiled examples, local-link validation, search,
and a shared responsive shell. Those are valuable foundations and should be
retained.

The principal weakness is usefulness at the user boundary. The generated
Python and CLI inventories are valuable correctness artifacts, but generation
has been treated too much as the end of the documentation job: long, flat
pages expose hundreds of entries with weak hierarchy, incomplete descriptions,
and presentation that is difficult to scan. At least one first-user tutorial
is not runnable in the environment it prescribes, and the expanded manuals
still stop before several products' headline results. Substantial maintained
examples also remain outside the rendered documentation system.

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
- Internal Rust support catalogs and Python stubs are generated from
  implementation-owned metadata. The Rust catalogs guard a deliberately named
  support boundary without becoming a second API browser; GammaLoop has deep
  Python runtime signature/default parity, and Linnet has a real import-surface
  check.
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

This revision records the agreed publication model: while pull request 96 is
open, pushes to its in-repository `docs` branch publish mutable `latest` as a
development review site. The workflow checks the pull-request state before
granting that temporary exception. Ordinary pull requests remain artifact-only;
after pull request 96 closes or merges, `main` is the sole authoritative source
of mutable `latest` documentation. That is an implementation constraint, not an
open governance decision. Tagged snapshots remain independent and immutable.

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
  [`docs` temporarily publishes the development site while pull request 96 is open; `main` is authoritative after merge],
  [The temporary publisher must retire automatically with the pull request; other pull requests remain artifact-only],
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
  [Python, CLI, and Rust reference boundary],
  [Exact Clap actions, arity, aliases, conflicts, globals, and usage; a flat command index with 79 manpage-style command routes; 45 hierarchical settings-namespace routes; Python module indexes with one flat page per supported class or function; stable item/member/legacy links; structured Python docstrings; normalized signatures; documented enum variants; coalesced overloads; composite settings types/defaults; linked internal Python types; and no inline raw-stub dump. Rust instead uses concise crate/version/feature orientation pages that lead directly to native, revision-specific Rustdoc; source-backed Rust catalogs remain internal drift checks],
  [The non-empty gate does not yet measure substantive prose or task backlinks. Linnet Python member and parameter prose remains uneven; GammaLoop's 12 implementation-detail exports remain deliberately de-ranked, including eleven currently unreachable integration-result records; and responsive/light/dark/accessibility browser fixtures are not yet enforced],
  [Manuals and examples],
  [Every product overview now offers a task chooser. The hierarchy is consistently Start here, Tutorial, Guides, Reference, and Version history. GammaLoop state/integration workflows, Linnet algorithms and the Clinnet rendering/cache guide, Spenso Rust and Python networks, Idenso algebra, and Vakint evaluation include realistic source/configuration examples, expected invariants, failure interpretation, and exact reference links. The Vakint tutorial now follows the checked-out API rather than an incompatible published-source example],
  [Only three authored examples currently execute; Clinnet's documented commands do not yet come from a generated parser-parity inventory, and clean external beginner journeys, remaining repository/scientific examples, named individual ownership, and the scheduled licensed/heavy tier remain open],
  [Integrated corpus],
  [Legacy current material is classified and routed, authored prose is Typst-only apart from exact `README.md` and `AGENTS.md` compatibility files, long README manuals were reduced to concise compatibility indexes, every product has a common Version history route, and search is derived from rendered HTML so imported Linnest/Kurvst manuals and included companion versions are indexed],
  [Release coverage remains incomplete, repository examples are not yet fully registered, and a generated graph proving that every current concept has a canonical route and contextual inbound link remains open],
  [Developer architecture],
  [The reviewed current notes were corrected against implementation and carry constrained lifecycle, review, trigger, and verified-scope metadata. The builder now rejects every current note without a scope, all 14 current records are tied to reviewed code digests, and current high-level architecture is registered for all five products. Long imported pages are also checked for a single semantic page heading],
  [All 31 records still resolve to the placeholder owner; review aging, trigger-aware owner approval, contextual inbound-link coverage, and browser-level architecture-page checks remain incomplete],
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

Before the temporary publisher was enabled, the public site was not a preview
of this branch. At that review point its latest deployment was
#link("https://github.com/alphal00p/gammaloop/actions/runs/32073830002")[run 32073830002]
from revision `f46892d216f1`. Pushes to `docs` now build, validate, and publish
the complete development site while pull request 96 remains open. Once it
closes or merges, that branch loses publication authority automatically and
only `main` can update `latest`.

== Cross-site audit at this revision

The examples that initiated this review were symptoms rather than the audit
boundary. The follow-up review compared navigation, rendered output, search,
reference generation, task journeys, release coverage, source ownership, and
developer freshness across the whole site. It found these broader conditions:

The current registries contain 39 authored product routes, 12 canonical
journeys, and 31 developer records, of which 14 are current. The counts below
are derived from `docs/products/registry.toml`, `docs/examples.toml`, and
`docs/developers.toml`; they are not estimates from the rendered navigation.

#table(
  columns: (1.5fr, 3fr, 3fr),
  table.header([*System*], [*What this revision establishes*], [*What is still lacking*]),
  [Information architecture],
  [39 authored product routes now use one visible hierarchy: 5 Start pages, 5 Tutorials, 11 Guides, 6 authored Reference pages, and 12 Version history pages. Product and task navigation no longer mixes guides, interfaces, and releases under a generic Manual label, while all 23 moved `manual/*` routes remain as non-navigated redirects. The builder also rejects a PDF chapter set or order that differs from the HTML page registry],
  [Depth is still unbalanced: GammaLoop has 8 authored routes, Linnet 9, Spenso 10, Idenso 7, and Vakint 5. A shared installation/compatibility matrix and a maintained three-task map per product are absent],
  [Reference generation],
  [Python signatures drop `typing.`/`builtins.` noise; module indexes link one page per supported class/function; parameters render as tables; members and overloads use flat stable sections; adjacent and cross-symbol links replace nested disclosure cards. GammaLoop commands use one flat manpage each, while settings follow all inferred schema namespaces with parent/child navigation. Rust no longer duplicates declarations and members in a custom browser: each product has a concise orientation page, and authored tasks link to exact pages in the native Rustdoc built for that revision],
  [GammaLoop's 191 supported reachable Python members have prose, but parameters, meaningful task examples, and backlinks remain uneven across Python modules. Vakint's 27-row topology inventory exposes only name, loop count, and slot count—not the match expression, mass/propagator pattern, backend applicability, or epsilon-depth limits needed to compare an unrecognized integral. Clinnet is not in the generated CLI inventory. Linnest and Kurvst now expose their generated Tidy module APIs in HTML, but compact symbol indexes and browser-level presentation checks remain absent. Rustdoc federation and visual/browser regression coverage are not enforced],
  [Search and inclusion],
  [Product search is built from rendered HTML body prose and headings, so content imported or included by Typst—most visibly Linnest, Kurvst, and companion version histories—is searchable. Representative task-query fixtures cover installation, integration, drawing, and version lookups; the portal federates all products and developer notes; and five registry-backed citation entries plus per-product `Cite` links expose the software citation hub],
  [Synonym coverage, Rustdoc federation, browser-level search behavior, a generated orphan-content report, and product-to-method/publication mappings beyond the software citation remain open],
  [Product journeys],
  [GammaLoop now maps its supported Rust/Python workflow and integration persistence, Spenso adds a complete Python tensor/network guide, Linnet covers cuts and drawing boundaries, Idenso preserves its source-backed identity specification, and Vakint's setup/API contradiction is corrected],
  [GammaLoop's events and observables guide still needs a bounded executed integration through `observables_final.json` and HWU output with a checked result; Linnet needs a complete mutation/drawing result; Spenso's HEP library needs a worked tensor/projector workflow and convention table; Idenso needs a worked Rust pipeline; Vakint needs evaluated Laurent coefficients; scientific conventions and citations remain uneven],
  [Examples and outcomes],
  [12 canonical journeys record prerequisites, expected invariants, routes, owner groups, verification tiers, commands, and exact source markers],
  [Only 3 execute; 6 compile and 3 are syntax-only. The wider maintained example corpus, licensed/heavy scheduled tier, result-value checks, and named owners are not integrated],
  [Version history],
  [Every product has the same Version history entry point; the canonical core Linnet, Spenso, and Idenso changelogs now have rendered nested routes alongside the companion histories, and each PDF keeps that complete history together before generated API appendices],
  [GammaLoop and Vakint lack consolidated changelogs; current Idenso/Spenso versions exceed their latest narrative entries; Linnet's Python distribution lacks an independent history; and readers have no generated selector from `latest` to available immutable snapshots],
  [Developer architecture],
  [All 14 current notes require resolving anchors and reviewed digests; their registered scopes were rechecked, high-level system architecture is now included for all five products, and the rendered lifecycle/evidence boundary distinguishes current guidance from 17 proposal, investigation, superseded, or archived records],
  [Every one of the 31 records still has placeholder ownership; nine historical/investigation records lack immutable evidence revisions. Review expiry, CODEOWNERS integration, contextual product links, and automatic changed-trigger review remain open],
  [Semantic quality],
  [The deterministic HTML gate discovers and checks every non-redirect page owned by the shared portal/manual shell. The redesigned full build contains 305 such pages, including generated symbol, command, and namespace routes; invalid multiple-page-heading structures and nine headerless data tables in current architecture notes were corrected],
  [Standards HTML validation, axe/browser execution, keyboard/viewport/visual regression, spelling, and scheduled external-link health remain open],
)

This table is the current prioritization boundary. Uniform labels and generated
inventories are not counted as task completeness, and acknowledged missing
release notes are not counted as complete history.

== Acceptance status at this revision

#table(
  columns: (auto, auto, 4fr),
  table.header([*ID*], [*State*], [*Evidence or remaining boundary*]),
  [`DOC-001`], [Partial], [Compiled usage semantics and representative parser checks exist; the exhaustive invocation matrix remains open],
  [`DOC-002`], [Partial], [The misleading beginner promises were corrected, but five clean external runs and the scheduled tier are not complete],
  [`DOC-003`], [Complete], [Clean single-product previews rewrite ecosystem links and pass local-link validation],
  [`DOC-004`], [Partial], [Python and CLI now use interface-native indexes and flat detail pages rather than nested disclosure boxes; the 305-page deterministic semantic/no-JS audit passes, while browser visual and axe fixtures remain open],
  [`DOC-005`], [Partial], [GammaLoop's curated Python surface expands from 4 to 28 of 40 exports, all 191 supported reachable members and 96 of 103 supported parameter rows have source-owned documentation, and Spynso3's 14 top-level exports are promoted; implementation-detail prose, generated setter parameters, and task coverage remain open],
  [`DOC-006`], [Partial], [Catalog generation rejects empty public CLI, argument, or settings prose; current non-empty coverage is 79/79, 412/412, and 181/181, including nested array and tagged-union settings, but substantive-quality and task-backlink gates remain open],
  [`DOC-007`], [Partial], [Five task choosers, a normalized Guides hierarchy, deeper subject chapters, and a new Spenso Python workflow are present; checked beginner and real-value paths do not yet cover every headline capability],
  [`DOC-008`], [Partial], [12 canonical product journeys now have validated owner-group, route, tier, invariant, command, and exact source-marker metadata; their modes are 3 run, 6 compile, and 3 syntax. Remaining maintained repository/scientific examples and scheduled verification are not yet registered],
  [`DOC-009`], [Partial], [Canonical prose is Typst-only and the declared legacy corpus is routed; a generated zero-orphan content graph remains open],
  [`DOC-010`], [Partial], [Lifecycle, review triggers, anchors, and digests are checked; every current record must declare a scope and all 14 do. Named owners, review aging, and trigger-aware owner approval still require maintainer assignment],
  [`DOC-011`], [Partial], [Every product has a consistently named Version history route; canonical core Linnet, Spenso, and Idenso changelogs and the companion histories are nested below it. GammaLoop/Vakint histories, component-to-release coverage, and a generated selector for available immutable snapshots remain open],
  [`DOC-012`], [Partial], [The portal offers five task-first product routes and each overview has task-to-guide/reference links; a registry-owned shared task model and complete reference backlinks remain open],
  [`DOC-013`], [Partial], [Ordinary pull requests retain complete previews, and the temporary `docs` publisher is gated on pull request 96's exact open head revision; branch protection does not yet enforce an always-created path-routed documentation check],
  [`DOC-014`], [Partial], [Local links, fragments, containment, source paths, and search targets are checked; scheduled external checking remains open],
  [`DOC-015`], [Partial], [A required deterministic audit discovers all generated shell pages and checks their headings, table headers, control labels, landmarks, internal ARIA references, and no-JavaScript navigation; the current full build covers 305 pages. Standards HTML validation, browser/axe, spelling, keyboard, viewport, and visual gates remain open],
  [`DOC-016`], [Partial], [Search is federated across five products, developer notes, and registry-backed software citations; it indexes rendered body prose and imported/included headings and has representative task-query ranking fixtures. Product pages link their exact citation fragment, while product/developer pages have pinned source and issue routes. Synonyms, Rustdoc federation, and browser-level search checks remain open],
  [`DOC-017`], [Partial], [The source-backed Idenso convention specification is deep; equivalent auditable anchors remain uneven across products],
  [`DOC-018`], [Partial], [Native Rustdoc is the canonical version-specific API browser. Concise product orientation pages identify crates, versions, features, support expectations, and authored workflow guides; task pages link directly to exact Rustdoc items. Source-backed support catalogs remain internal coverage and drift checks, and the former custom `reference/rust/supported/*` browser routes are compatibility redirects rather than navigation or search targets. Redirect-fragment and Rustdoc-federation coverage remain to be enforced],
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

=== P0. Generated CLI syntax needs exhaustive parity evidence

At baseline the exported schema omitted action and arity, so the renderer
falsely documented the boolean `--clean-state` flag as requiring a value and
mixed user operations with generated help/choice nodes. The current schema and
renderer preserve exact usage, actions, value requirements, aliases, conflicts,
globals, and arity; generated internals are separated from primary navigation,
and representative parser fixtures cover the corrected semantics.

The remaining integrity boundary is exhaustive parity. A generated page that
looks complete while one less-common nested/alias/positional combination is
wrong is still more dangerous than an obviously partial page.

*Implemented contract and remaining gate*

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

The baseline exposed three misleading promises: a Vakint sibling-path setup
whose workspace override was not inherited, a Linnet description that
overstated the identifiers/subgraph operation shown, and a GammaLoop “First
calculation” that stopped after state generation. This branch corrects the
Linnet promise, renames GammaLoop's entry journey to “Create your first state”,
and makes Vakint's `latest` tutorial follow the checked-out source API and Nix
workspace instead of combining incompatible source and release contracts.

The residual P0 is execution evidence. Only three canonical journeys run;
compile/syntax checks do not prove the five empty-environment tutorials or
their promised scientific results.

*Required change*

- Run every beginner tutorial verbatim in a clean temporary environment outside
  the workspace.
- Keep mutable `latest` source instructions distinct from immutable released
  package instructions, and parity-test both supported installation boundaries.
- Add a bounded GammaLoop integration/result journey and equivalent
  result-bearing paths where the other products still stop at compilation.
- Give every procedural example prerequisites, expected output or invariant,
  expected runtime/resource class, failure modes, and a next step.

*Complete when*

- all five beginner journeys pass in CI from empty temporary directories;
- every fence is declared `run`, `compile`, or `syntax` and CI rejects an
  unclassified fence; and
- heavy, licensed, or external-tool examples are visibly classified and run in
  a provisioned scheduled tier rather than silently skipped.

=== Resolved P0. Product-scoped previews have an honest navigation contract

At baseline `just docs-site PRODUCT` emitted a five-product portal while only
building one product and skipped generated-link validation, leaving four sets
of broken links. Product-scoped builds now rewrite ecosystem routes for the
preview they actually contain and always validate local links, fragments, and
search targets. Clean preview fixtures cover each product. Keep this as a
required regression contract; it no longer belongs in the open work queue.

=== P1. Generated references are not yet a uniformly good experience

The reference pipeline has strong implementation provenance, but the rendered
result is mechanically broad and editorially shallow. Autogeneration is a
sourcing mechanism; it does not reduce the need for hierarchy, explanation,
examples, or visual design.

At the reviewed baseline, the Python reference exposed 95 top-level exports
across five modules as long sequences of generic cards. It lacked a local
symbol index and rendered NumPy-style sections, recursive parameter objects,
internal kind labels, feature flags, and generated `help(...)` placeholders as
though they were user documentation. GammaLoop rendered 40 exports but curated
only four, all of its rendered members lacked prose, and Linnet and Spenso also
had large undocumented member counts. Inventory parity was not documentation
parity.

This branch fixes the Python rendering defects: supported items come first;
signatures and kind labels are normalized; docstrings become Parameters,
Returns, Raises, Notes, and Examples; parameters use tables instead of nested
disclosures; each supported class/function owns a flat page with stable member
anchors and adjacent-symbol navigation; internal types and exact source items
are linked; and the downloadable stub no longer occupies the human page.
GammaLoop now curates 28 of 40 exports and documents all 191 supported
reachable members and 96 of 103 supported parameter rows; Spynso3 promotes all
14 top-level exports.
Implementation-detail prose, generated setter parameters, task examples, and
backlinks remain uneven, so the quality finding remains open even though the
visual pathology is corrected.

The rendered artifact makes the remaining editorial debt measurable:

#table(
  columns: (1.1fr, 2.8fr, 2.8fr),
  table.header([*Surface*], [*Missing substantive prose*], [*Missing task evidence*]),
  [GammaLoop Python], [All 191 supported reachable members are documented; 7 of 103 supported parameter rows are generated setter values without distinct prose. Across the deliberately separated implementation-detail exports, 65 of 71 members and the only parameter lack descriptions], [Only 13 genuine Examples sections; automatic `help(...)` placeholders are no longer rendered as examples],
  [Linnet Python], [56 of 180 members and all 133 parameter rows lack descriptions], [No genuine Examples section; automatic introspection placeholders are no longer added],
  [Vakint Python], [All 20 public members are documented; 41 of 42 parameter rows have prose], [Task examples do not yet cover the complete backend and result-interpretation lifecycle],
)

These are usefulness measures, not a demand to promote every implementation
detail. Python's supported inventory should be deliberate, and every promoted
item should satisfy the stronger prose/example contract. Rust support catalogs
retain source-backed entry-point and drift checks internally, but their former
member and parameter counts are no longer public-renderer acceptance targets;
native Rustdoc owns exhaustive API presentation.

At baseline the CLI/settings page put 285 command-tree nodes, 412 arguments,
and 166 settings on one flat route, with false value syntax and extensive empty
prose. The current renderer preserves compiled action/arity semantics, separates
public commands from help-tree internals, gives all 79 public commands a flat
manpage-style route, organizes 181 current settings into 45 schema-derived
namespace routes, adds filtering and stable legacy relocation links, and
rejects empty public descriptions. It still needs a substantive-quality/backlink
contract and browser evidence that the densest pages remain usable at narrow
widths.

The imported Linnest and Kurvst manuals exposed a different but equally visible
reference gap at baseline: their prose rendered in HTML, while each generated
Tidy module reference was guarded by `target() == "paged"`. This branch now
renders those inventories in both HTML and PDF and documents the actual bundled
source-package import paths. The remaining work is to prove stable per-symbol
linking, provide compact indexes for the long inventories, and cover their web
presentation in a browser rather than merely proving compilation.

*Required change — Python*

- Build an information hierarchy of module, class/function, member, and
  parameter with stable URLs or stable indexed routes, breadcrumbs, a local
  symbol index, kind filters, and search.
- Parse supported docstring structure into semantic sections and human labels;
  render normalized readable signatures while retaining the exact signature
  and downloadable stub.
- Put the curated supported surface first. Collapse, separate, and de-rank
  implementation details unless users need them; when they are not part of the
  compatibility surface, keep them out of public symbol routes and search.
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
- Link each primary command page and settings namespace to the manual tasks that
  use it, and link those guides back to exact reference anchors.

*Required change — presentation*

- Use the shared manual shell and design tokens for every reference route, with
  readable long signatures, responsive parameter/settings layouts, useful
  local navigation, and restrained visual density.
- Add compact Linnest and Kurvst symbol indexes, stable per-symbol anchors, and
  focused search entries for their now-rendered HTML inventories.
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
- each module and primary command family or namespace links to meaningful,
  tested examples for its principal tasks;
- Linnest and Kurvst symbols have stable direct links and focused search entries;
  and
- representative pages have no viewport-level horizontal scroll at 375 px and
  no serious or critical accessibility findings.

=== P1. The manuals cover a template, not each product's real task surface

At the reviewed baseline each product had exactly five short authored pages.
The current branch replaces that quota with a common Start, Tutorial, Guides,
Reference, and Version history hierarchy and grows the products to between five
and ten authored routes. The shape is easier to navigate, but route count is
not task completeness and the remaining depth is still uneven:

#table(
  columns: (1.2fr, 2.2fr, 3.5fr),
  table.header([*Product*], [*Strong current entry point*], [*Highest-value missing worked path*]),
  [GammaLoop], [State generation, inspection, sample APIs, integration slots, workspace resume, and source-level event/observable inspection], [A bounded executed integration carrying the existing observables/selectors guide through `observables_final.json` and HWU interpretation, precision/stability, and validation],
  [Linnet], [Half-edge construction, subgraph-aware traversal, cuts, and drawing boundaries], [A checked mutation-mapping workflow and one complete layout/rendered drawing result],
  [Spenso], [Rust contraction/network execution and Python tensor/network ownership], [Planning diagnostics, compiled evaluator lifecycle, a checked representative network result across both interfaces, and an owned HEP conventions/workflow chapter covering gamma/projector and SU(3) data],
  [Idenso], [Python metric/dummy-index passes and a deep source-backed FORM rule specification], [A complete Rust identity pipeline with reversible cooking and checked Dirac/color output],
  [Vakint], [Current-source topology matching, normalization, and explicit backend policy], [One analytic or numerical evaluation with expected Laurent coefficients and backend reproducibility, plus topology entries that expose match/mass patterns and backend/epsilon limits rather than only counts],
)

GammaLoop and the new Spenso Python guide now deep-link the exact generated
items they use, but equivalent backlinks remain sparse in other tutorials and
on reference items themselves. Authored item links currently number 27 Python,
11 Rust, and 12 CLI for GammaLoop; Linnet has 0 Python and 4 Rust, Spenso 9 and
5, Idenso 0 and 3, and Vakint 0 and 4. Some snippets are still compile-only fragments
where users need an observed result. Conceptual physics claims also need
stronger anchors: Local Unitarity, normalization/sign conventions, and
precision-sensitive backend behavior should point to defining references and
regression-backed examples.

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
no separate history. The core Linnet, Spenso, and Idenso changelog sources are
now rendered directly under their products' Version history sections, so those
records are discoverable without maintaining a duplicate. Their presence does
not fill the missing Spenso/Idenso releases or create histories for GammaLoop
and Vakint. Several prose pages also repeat manifest versions or external-tool
minima that the builder already knows from generated metadata.

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

At the reviewed baseline, at least 4,203 lines of current guidance, references,
plans, dated evidence, and raw output lived outside the product and developer
registries. The current branch completes that declared migration: maintained
prose is canonical Typst, architecture records are lifecycle-classified,
logging/drawing/Idenso material has rendered routes, compatibility READMEs point
to canonical pages, and the source-format gate rejects a new parallel Markdown
or HTML manual.

The portal now turns the ecosystem into task routes for calculations, graphs,
tensors, identities, and integrals and distinguishes the application from its
supporting layers. The remaining integration problem is finer grained: there is
no generated content graph proving that repository examples, shared concepts,
generated items, developer records, and every current page have one canonical
route plus contextual inbound links. Companion tools are also only partly in
that graph: Kurvst, Linnest, and Clinnet have manuals, but they are not all
first-class versioned components in the product registry, and Clinnet lacks the
compiled-parser reference contract used by GammaLoop. Migration is complete;
connection, component provenance, and ownership are not.

*Established invariant*

- Maintained prose is Typst; only files named exactly `AGENTS.md` and
  `README.md` are compatibility entry points, and the source-format check
  rejects parallel Markdown/HTML manuals.
- Product, developer, proposal, investigation, archive, and supersession
  lifecycles are registry-owned. Current diagnostics, drawing, and Idenso
  reference material has rendered routes; READMEs link to those routes.
- Portal task routes and federated rendered-content search provide the first
  integrated views without copying the underlying prose.

*Still required*

- Generate the content/reference/example graph and fail orphan current sources,
  missing canonical routes, duplicate editable concepts, and missing
  contextual inbound links.
- Extend ownership/review metadata from developer records to current product
  pages and canonical examples.
- Split long living logs into current contracts plus dated evidence and finish
  contextual links from product guides to relevant developer records.
- Register independently versioned companion tools as components, expose their
  compatibility and provenance, and generate Clinnet's command reference from
  its compiled parser.

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

=== P1. “Current” developer architecture lacks owned freshness enforcement

Schema 3 now separates lifecycle from display status, requires review metadata
and triggers, resolves stable anchors, hashes verified code scopes, rejects a
current record without a scope, validates supersession/evidence, and renders
those facts separately from the build revision. All 14 current records pass
that contract at this revision, including high-level architecture for all five
products. This closes the earlier failure mode where a
recent-looking “Implemented” page could remain green without any code evidence.

The remaining weakness is governance and ecosystem composition. All 31
records resolve to the `unassigned` placeholder, scope digests are refreshed by
the documentation author rather than an accountable owner, trigger globs do not
yet map a changed pull request to a required review, and age warnings do not
escalate into owner review or publication failures. A digest proves that the
reviewed bytes did not change; it does not prove that the chosen scope is
sufficient or that the prose
is architecturally complete. The five product notes intentionally describe
their own cores; no current contributor map composes their dependency,
language/binding, code-generation, release, and ownership boundaries across
GammaLoop, Linnet/Clinnet/Linnest, Spenso/Spynso3, Idenso, Vakint, and Kurvst.

*Still required*

- Assign the portal/toolchain and each product/runtime architecture area to a
  resolving owner, protect those mappings through CODEOWNERS, and require owner
  disposition when a review trigger or verified scope changes.
- Ratify review service levels, then warn or block overdue current records while
  leaving frozen investigation/archive evidence governed by provenance rather
  than age.
- Generate the reverse reference graph and require contextual inbound links,
  not only a developer-index entry.
- Give investigations complete immutable capture provenance and split long
  living engineering logs into a short current contract plus dated records.
- Keep the five product system notes connected to product workflows
  and expand their verified scopes when architecture ownership moves, rather
  than treating registration itself as permanent freshness.
- Add one source-backed ecosystem map that composes the five product notes and
  companion packages without duplicating their internal descriptions.

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

The accepted delivery flow is now clear: `docs` is the temporary live
development publisher while pull request 96 is open, and `main` is
authoritative after merge. The automation should encode that lifecycle rather
than rely on a remembered cleanup. Updates to `docs` should build, validate,
and publish the complete site once, without duplicating the pull-request build.
Ordinary pull requests remain artifact-only. When pull request 96 closes or
merges, the workflow must reject further `docs` publication and the merged
`main` revision should rebuild and deploy the reviewed sources. This is not a
choice between two long-lived documentation sources.

At the audited baseline, the Pages workflow had no pull-request trigger and its
validation set differed from the documentation-crate checks. The current
workflow now produces a complete rendered quality result before publication,
and the temporary `docs` path uses that same result as its publication input.

The current validation still has bounded blind spots:

- local generated pages, fragments, source paths, containment, redirects,
  search targets, and included Rustdoc trees are checked; external URLs remain
  outside the deterministic link gate;
- the full-corpus semantic/no-JavaScript audit is deterministic but is not a
  standards HTML validator or browser accessibility engine; and
- spelling, keyboard behavior, viewport rendering, and visual regression are
  not yet checked.

*Required change*

- Add one required, path-routed “Documentation quality” pull-request check.
- Build and retain the complete rendered preview artifact on ordinary pull
  requests without letting those revisions mutate `latest`.
- Let pushes to the in-repository `docs` branch publish the development
  `latest` site only while pull request 96 is open, and skip the duplicate
  pull-request event for the same head revision.
- Retire that exception automatically when pull request 96 closes or merges;
  then promote `latest` only from `main`. Keep scheduled mutable builds on
  `main` and preserve tag-triggered immutable snapshots independently.
- Run content checks for content-only changes and the full six-crate docs test
  group for renderer, schema, exporter, script, or workflow changes.
- Keep local checks deterministic; run flaky external HTTP checks separately on
  a schedule with retries and an allowlist.

*Complete when*

- no production deployment is the first full validation of its commit;
- a push to `docs` publishes the validated development site exactly once while
  pull request 96 is open, cannot publish after it closes or merges, and the
  merged `main` commit can publish the same content;
- local links, fragments, Rustdoc, repository source paths, search entries,
  redirects, and asset references have zero failures;
- link resolution cannot escape the generated root;
- every preview is retained long enough for review; and
- the job summary reports example modes, page counts, link counts, cache reuse,
  and step timings.

=== P2. Accessibility, progressive enhancement, and search are only partially measured

The shell has good semantic foundations, revision-pinned source/issue feedback,
federated search, and deterministic no-JavaScript navigation checks. Search now
indexes both body prose and headings from the rendered page rather than only
the wrapper Typst, so imported manuals and included release sections are
visible. Representative task-query fixtures exercise installation,
integration, drawing, and version discovery. The semantic gate discovers every
non-redirect page owned by the shared shell, rather than sampling a fixed route
list, and rejects malformed headings, headerless data tables, broken internal
ARIA references, and no-JavaScript navigation failures across the full corpus.

Those are structural checks, not browser evidence. Search still lacks synonym
coverage and Rustdoc federation, and its maintained scorer checks do not
exercise the dialog, keyboard interaction, or actual layout in a browser.
There is no standards HTML validator, browser/axe run, keyboard or viewport
smoke suite, spelling or visual regression gate, CODEOWNERS file, canonical
metadata, or sitemap. Typst's experimental HTML export can still ignore layout
behavior inside third-party drawing components, so representative visual
fixtures remain necessary even when compilation is clean.

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
+ *Reference:* concise Rust orientation plus native version-specific Rustdoc,
   and generated Python, CLI, settings, and topology facts.
+ *Developer record:* implemented architecture with verified code scopes and
   review triggers, plus proposals, investigations, and archives with explicit
   lifecycle and evidence provenance.

Compatibility READMEs should orient contributors and link to the rendered
material, not act as second user manuals. Generated Python, CLI, settings, and
topology reference is a first-class product surface paired with task
documentation; generation is a sourcing mechanism, not an exemption from
editorial or visual quality. Rustdoc remains native, while authored orientation
and task pages explain support boundaries and workflows. Dated evidence should
be preserved, but never presented as the current contract without an explicit
current summary.

The views must be connected. A manual task links to the exact symbols,
commands, settings, concepts, and canonical examples it uses. Generated
Python/CLI items link back to the guides and examples that give them purpose;
Rust orientation pages connect authored workflows to exact Rustdoc items
without wrapping Rustdoc in a second member browser. Shared concepts have one
maintained explanation, while contextual links make them visible from every
relevant product.

== Reference-quality design requirements

=== Rust

Use native Rustdoc as the canonical exhaustive API browser for the exact
revision and feature matrix being documented. Do not reproduce its signatures,
members, implementations, trait relationships, or search in a custom public
renderer. Each product's `reference/rust/` route should instead be a concise
orientation page that identifies crate roles, versions, enabled features,
support/stability expectations, and the authored workflow guide, then links
directly to each native crate root.

Authored tutorials and guides should link exact Rustdoc item or method URLs when
they name an API. Public items not selected by a maintained workflow are not an
implicit compatibility promise merely because Rustdoc exposes them. Existing
`reference/rust/supported/*` URLs remain non-indexed compatibility redirects;
they must not appear in navigation, search, the PDF, or canonical route lists.

Keep source-owned comments and curated catalog annotations as internal checks
where they detect a deliberate support boundary, missing source provenance, or
drift in key entry points. They should not create a second public member browser
or force maintainers to duplicate parameter and example prose that belongs in
Rustdoc and task documentation. Complete when every configured Rust crate has a
validated orientation link and native crate index, authored deep links resolve,
legacy redirects have a safe fallback, and federated search points to Rustdoc
rather than retired catalog fragments.

=== Python

Keep the runtime inventory and checked stubs as the source of signatures and
availability. Add editorial information at the implementation boundary through
binding docstrings, structured annotations, and validated grouping metadata;
do not create a second handwritten signature catalog.

The renderer should provide:

- a filterable module index followed by one stable page per supported public
  class or function, with stable member fragments on the owning symbol page;
- fully qualified names, human-readable kinds, normalized and exact
  signatures, source provenance, availability, and version information;
- semantic Parameters, Attributes, Returns, Raises, Notes, and Examples
  sections rather than escaped raw docstring layout;
- a curated supported surface in primary navigation, with implementation-only
  exports retained as drift evidence but omitted from public symbol pages and
  search;
- local symbol filtering, breadcrumbs, a useful table of contents, copyable
  signatures/examples, adjacent-symbol navigation, and responsive tables; and
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

The top-level page should be a filterable command index, not print the full
database. Every public command should have a flat manpage-style route with its
synopsis, positional arguments, local options, inherited/global options, and
parent/direct-child links. Settings should use a separate namespace index plus
one page per actual schema namespace, with parent/child navigation and only the
settings directly owned there. Legacy fragments should lead to these canonical
routes without adding hundreds of entries to a page table of contents.

=== Shared visual contract

All generated reference uses the same shell, tokens, navigation behaviors,
feedback routes, and accessibility targets as the manuals. Reference-specific
components should be designed for dense technical data without becoming a wall
of equally weighted cards: compact indexes, flat semantic sections, readable
line wrapping, responsive field layouts, clear status badges, and restrained
use of borders and pills. Disclosures are reserved for genuinely secondary
content, not the primary class, member, command, or settings hierarchy.

Visual regression fixtures should cover a module index, a large Python class,
a function with rich docstrings, a command manpage, and a settings namespace in
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
  [`DOC-011`], [P1], [Make release/version/tool claims generated and complete and expose available immutable snapshots], [Manifest-to-release coverage is current for every component and readers can switch versions],
  [`DOC-012`], [P1], [Use ecosystem metadata for task navigation and bidirectional manual/reference links], [Task chooser and contextual cross-link fixtures pass],
  [`DOC-013`], [P1], [Add the PR quality gate and enforce the temporary docs-publisher/main-promotion flow], [While pull request 96 is open, a `docs` push validates and publishes once; after it closes or merges, only `main` can update `latest`],
  [`DOC-014`], [P1], [Complete link, source-path, containment, and redirect validation], [Zero internal failures, including Rustdoc; external report is scheduled],
  [`DOC-015`], [P2], [Add browser, HTML, accessibility, spelling, and no-JS checks], [Zero HTML errors and serious/critical axe findings on representative pages],
  [`DOC-016`], [P2], [Federate search and add edit/report feedback routes], [Task and symbol fixtures rank the correct page in the first five results],
  [`DOC-017`], [P2], [Add scientific citations and regression-backed conventions], [Each precision/sign/normalization-sensitive claim has an auditable anchor],
  [`DOC-018`], [P1], [Make native Rustdoc canonical and keep the custom Rust layer orientational], [Every configured crate has a version/feature-aware orientation link and validated Rustdoc root; authored deep links target Rustdoc; internal support audits pass; retired browser URLs redirect outside navigation, search, and PDF output],
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

Complete DOC-004 through DOC-006 and DOC-018, establish the inventory/schema
portion of DOC-009, establish the developer lifecycle/scope schema and initial audit under
DOC-010, and complete DOC-013. The shared reference design is now applied to
all supported Python symbols, GammaLoop commands, and settings namespaces;
finish this phase with representative browser fixtures rather than another
mechanical renderer expansion. Publish coverage and developer-scope reports as
warnings, fill the source-owned documentation gaps, then promote the substantive
coverage and freshness thresholds to required checks.

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

== Next vertical implementation slice

The completed slice established the registry, reference, source-format,
freshness, and build foundations. The next slice should prove task completeness
where the main product is still weakest:

+ Add a bounded GammaLoop integration card whose runtime is suitable for a
  provisioned check, and assert a stable result or tolerance instead of only
  compiling its commands.
+ Extend the existing events and observables guide through an executed
  integration that produces `observables_final.json` and HWU output; explain
  raw statistics and merge ownership, and link every CLI setting plus
  Rust/Python record it uses.
+ Register that journey as the canonical source, classify its fast and
  licensed/heavy tiers, and publish the last verified revision when the heavy
  tier cannot run on every pull request.
+ Add representative browser/axe and 375/768/desktop visual fixtures for the
  GammaLoop Python and CLI reference, plus search queries for the new workflow.
+ Use the same acceptance shape—not copied prose—to complete one result-bearing
  journey in Linnet, Spenso, Idenso, and Vakint.

This slice tests whether the improved structure actually carries a user from
setup to an interpreted scientific result. More generated entries alone do not
satisfy it.

== Scorecard

The quality job should publish these measures on every canonical build:

#table(
  columns: (4fr, 1.4fr),
  align: (left, right),
  table.header([*Measure*], [*Target*]),
  [Non-main revisions able to mutate `latest` outside the open pull request 96 exception], [0],
  [False generated CLI invocation syntaxes], [0],
  [Registered public CLI commands, arguments, and settings with substantive descriptions], [100%],
  [User-facing CLIs inventoried with generated parser parity, including Clinnet], [100%],
  [Configured Rust crates with a validated native Rustdoc root and orientation entry], [100%],
  [Authored Rust API links targeting canonical Rustdoc rather than the retired browser], [100%],
  [Retired Rust browser URLs present in navigation, search, or PDF output], [0],
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

== Remaining maintainer decisions

The publication source is no longer an open choice: `docs` is the temporary
development publisher only while pull request 96 remains open, and `main`
alone is authoritative for mutable `latest` afterward. A long-lived
per-branch selector is deferred unless multiple simultaneous documentation
branches make its routing, retention, and cleanup costs worthwhile.
Maintainers still need to record two ownership and execution choices:

+ Who owns the portal/toolchain, generated Python and CLI reference quality,
   each of the five product manuals, and the GammaLoop runtime/UV,
   Spenso/Idenso parsing, build-system, and documentation-system architecture
   records? Owners must also approve the current/proposal review service levels.
+ Which scientific examples are hermetic pull-request checks, and which require
   licensed or scheduled infrastructure?

Implementation can continue in small acceptance-focused changes while those
assignments are made, but ownership-dependent freshness enforcement and the
scheduled scientific tier cannot be declared complete until accountable people
and infrastructure are recorded.
