= Manual content debt register
<manual-content-debt-register>

#quote(block: true)[
  *Status:* Open register

  *Review date:* 2026-08-21

  *Reviewed baseline:* `4e77d8804aa3394992572575dab5d3f9dedf65a9`

  *Scope:* Authored product manuals, their canonical examples, scientific
  terminology, method citations, result interpretation, release narrative, and
  product handoffs. Visual presentation and item-level API documentation are out
  of scope.
]

== Purpose

This is the canonical backlog for known prose-manual content deficiencies. The
#link("documentation-improvement-plan.typ")[documentation improvement plan]
sets the cross-product strategy and quality gates. The companion
#link("api-documentation-debt-register.typ")[API documentation debt register]
tracks signatures, item contracts, support boundaries, and local API examples.
This register retains findings about audience, scientific explanation, complete
workflows, terminology, and the interpretation of results.

Page count and compilation alone do not close a manual-content item. Resolution
requires canonical Typst prose for an explicitly named audience and task,
scientific review and citations where relevant, honest prerequisites and
runtime/license limits, and a verified workflow that reaches an interpretable
result or invariant at its stated verification tier.

== Audience contract

The documentation should become progressively more computational as readers move
deeper into the stack:

1. GammaLoop starts from a collider-physics question, states the supported physics
   and conventions, produces a result with an uncertainty, and teaches validation.
2. Vakint starts from a mathematical integral and Idenso from a convention-explicit
   tensor identity.
3. Spenso then maps indexed notation into representations, slots, storage, and
   networks.
4. Linnet and Linnest may foreground graph models, invariants, algorithms, and
   programmable interfaces after showing the scientific operation they support.
5. Native Rustdoc and developer records own detailed types, storage, wire formats,
   logging internals, complexity, and implementation architecture.

Every overview should state its intended reader, primary task and interface,
expected programming fluency, and “use this directly when / use GammaLoop when”
boundary. A GammaLoop reader should not need to infer the crate architecture to
complete or assess a calculation.

== Audit summary

#table(
  columns: (1fr, 1.8fr, 2.8fr),
  table.header([*Product*], [*Content to preserve*], [*Critical remaining gap*]),
  [GammaLoop], [Careful state, persistence, command, and operational boundaries; new scope and conventions chapters establish the current scientific boundary.], [The primary path still stops at state creation or point inspection. No maintained cross section or distribution is carried through uncertainty, convergence, normalization, and independent validation.],
  [Vakint], [Clear matching, reduction, evaluation, backend, and citation stages.], [Software notation is not mapped to a mathematical integral; normalization and checked Laurent coefficients are absent.],
  [Idenso], [Deep source-backed rule material and strong rewrite-stage cautions.], [Shipped scientific behavior, target specifications, implementation evidence, and historical validation are mixed in one public reference.],
  [Spenso], [Strong structure/data/execution model and useful failure ownership.], [The bridge from indexed physics to representations, slots, duality, networks, and HEP conventions is missing.],
  [Linnet], [The strongest progressive conceptual introduction, including the physical value of half-edge boundaries.], [The maintained strategy still needs the checked mutation and complete graph-result journey recorded under `DOC-007`.],
  [Linnest], [Broad graph, layout, generic drawing, and application-template composition.], [Archive/plugin payload mechanics and advanced callbacks precede the smallest useful drawing result.],
)

== Maintenance contract

Identifiers are stable. Use `Open`, `In progress`, `Blocked`, or `Resolved`. Add or
refine an entry in the same change that discovers a prose deficiency. Mark it
resolved only when its completion condition is backed by the canonical prose,
scientific review where required, and the relevant example-verification record.
Keep resolved entries through the next documentation review with the resolving
change ID, then remove them from the active register.

Do not duplicate item-level deficiencies from `APIDOC-*`. A manual item may link an
API item when a complete user journey depends on it. Broad delivery work already
owned by `DOC-*` remains in the strategy register and is referenced below with its
concrete manual manifestation.

== GammaLoop manual

=== MANUAL-101 · Physics scope and Local Unitarity contract

*Priority:* P1 · *Status:* In progress

The overview names Local Unitarity and immediately switches to software subsystems;
the process guide gives only one sentence about forward-scattering graphs and cuts.
The manual does not let a physicist determine whether a proposed calculation is
partonic or hadronic, supported or experimental, or dependent on an unavailable
model, contribution, backend, tool, or license.

*Source-audit findings (2026-08-21):* the demonstrable runtime path uses fixed partonic external
momenta and a PDF factor of one; no PDF convolution or factorization scale was found. Automatic
polarization sums cover scalar, spinor, and vector states, while complex masses are unsupported.
Built-in JSON models, optional UFO import, observables, clustering, and multiple evaluator/UV
backends have narrower feature and tool boundaries than the overview implies. One evaluation can
produce multiple correlated event groups which together form a single statistical sample.
Compiled loop-count ceilings are implementation bounds, not evidence of complete physics coverage
at those orders.

*Progress (`nlymvowy`):* a new physicist-facing scope chapter now explains the real/virtual
reorganization at method level, cites the defining Local Unitarity papers, separates their
all-order formal claim from release capability, explains correlated event groups, and publishes
a conservative matrix for initial states/PDFs, generated objects, filters, models, spins, masses,
subtraction, integration components, observables, tools, and licenses. A complete account of the
supported subtraction envelope and a maintained end-to-end cross-section validation remain open.

*Complete when:* a physicist-facing chapter explains the forward-graph/cut relation,
real/virtual combination, common integration variables, local cancellation, the role
of correlated events and UV/IR/threshold subtraction, and applicable observable
conditions. A maintained matrix classifies supported, experimental, and unsupported
initial states, orders/contributions, models, masses/spins, PDFs, scales and schemes,
observables, external tools, and license requirements, with defining method citations.

=== MANUAL-102 · Scientific conventions and GammaLoop glossary

*Priority:* P1 · *Status:* In progress

Settings expose convention-sensitive choices, but no authored chapter assembles the
scientific contract. Terms such as state, event group, graph group, raised cut,
orientation, loop-momentum basis, E-surface, integration coordinate, and target are
introduced without one canonical physical definition.

*Source-audit findings (2026-08-21):* the implemented metric is `diag(+,-,-,-)` and stored external
momenta are positive-energy four-vectors with separate incoming/final signatures. Scalar helicity
zero is the physical card convention but is not enforced by the evaluator. The one- and
two-particle flux formulas and the sample × parameterization-Jacobian × integrator-weight relation
are source-verifiable, while reporting-unit conversion is applied only in the two-incoming flux
branch. Powered and unpowered constraints populate different filter containers; the powered
exponent is discarded, and existing differential regression cards therefore do not establish a
reviewed perturbative-order definition. A universal phase-space/amplitude normalization,
automatic color average, complete symmetry-factor contract, perturbative-order mapping, and the
energy-unit premise behind picobarn conversion remain unverified and must not be inferred.

*Progress (`nlymvowy`):* a canonical conventions chapter now records four-vector order, metric,
physical momentum signatures, dependent-momentum conservation, helicity/spin averaging, explicit
color ownership, one/two-particle flux, reporting units, the $I J w_"MC"$ sample relation,
correlated event reconstruction, normalization ownership, scale distinctions, and overloaded
workspace/event terms. The process guide now distinguishes amplitude-side, cross-section-side,
and unresolved-cut selectors and explicitly warns that the repository's current differential
cards are regression fixtures rather than audited order/normalization exemplars. Complete
phase-space/amplitude normalization and method-author reviews of scheme and threshold conventions
remain explicitly open. Parser/runtime correction remains tracked by
#link("api-documentation-debt-register.typ#apidoc-007-gammaloop-accepted-but-unsupported-calculation-options")[APIDOC-007].

*Complete when:* the manual states metric and momentum flow, integration and
phase-space measures, flux, symmetry, color and helicity factors, units, amplitude
normalization, scale/scheme conventions, event-weight composition, and the relation
between integrand value, Jacobian, integrator weight, and result. It maps `{n}`,
`{{n}}`, `QCD=n`, powered coupling constraints, amplitude-loop count, forward-loop
count, and LO/NLO terminology, and links every overloaded term at first use.

=== MANUAL-103 · Result interpretation and scientific validation

*Priority:* P1 · *Status:* Open

The process guide explains where result files live, but not how their estimates,
errors, components, slots, or observable bins should be interpreted or trusted. The
scientifically useful inspection, limiting-profile, precision, and stability commands
are not assembled into a validation workflow.

*Complete when:* one executed audit explains the authoritative checkpoint versus
summary JSON and observable JSON/HwU outputs; total and differential estimates,
units, real/imaginary components, Monte Carlo error, bin normalization and
correlations; target versus stopping accuracy; convergence and chi-squared;
zero/invalid samples and precision escalation; and seed, parameterization,
subtraction-parameter, scale, and independent-reference checks with explicit pass/fail
criteria. The result-bearing calculation becomes the default continuation from the
physicist tutorial rather than an advanced audit. Link `APIDOC-201` and `APIDOC-202`
for Python result/event contracts.

=== MANUAL-104 · Public troubleshooting and companion-tool boundary

*Priority:* P1 · *Status:* Open

GammaLoop's public diagnostics chapter becomes a contributor tracing-DSL manual with
callsite metadata, filter implementation, routing prefixes, and tag-naming guidance.
The GammaLoop manual also embeds the complete Kurvst geometry API despite drawing
being a companion workflow.

*Ownership finding (2026-08-21):* the GammaLoop page calls Kurvst its own layer, while Linnest
links back to that route as canonical. Resolution first requires assigning Kurvst one canonical
companion-documentation owner and route; GammaLoop and Linnest should then retain task-specific
handoffs rather than competing or embedded copies.

*Complete when:* public troubleshooting is organized by physicist-facing symptoms,
diagnoses, safe actions, and useful bug-report evidence. Logging internals and naming
rules live in developer material. GammaLoop retains a focused diagram-customization
handoff and links to independently owned Linnest/Kurvst documentation instead of
embedding their complete manuals. Beginner and physics workflows require no tracing
DSL, Rust type, archive format, orientation table, or curve-API knowledge.

== Domain and graph manuals

=== MANUAL-201 · Vakint mathematical contract

*Priority:* P1 · *Status:* Open

The first `topo(prop(...))` expression is never translated into a mathematical
integral, and the evaluation guide names normalization components without providing
the convention. The topology inventory cannot be used to compare a mathematical
integral with entries whose meaning is reduced to a name and counts.

*Result-contract findings (2026-08-21):* the worked evaluation selects
`LoopNormalizationFactor::MSbar` and then calls the output an “MS-bar Laurent series,” without
distinguishing a loop-measure normalization from pole subtraction or renormalization. The optional
backend uncertainty is not classified as absolute or relative, deterministic or statistical, or
defined coefficient by coefficient. The worked AlphaLoop path is absent from the adjacent methods
and citations list.

*Complete when:* a representative expression is mapped argument by argument to its
integral, including dimension/epsilon, metric, loop measure, (+i0), mass and scale,
momentum routing, propagator power, and external-momentum conventions. Topology
entries expose denominators or diagrams, routing, masses, active propagators, backend
support, and epsilon limits. One checked Laurent series records coefficients,
tolerance, backend/tool versions, and normalization. It distinguishes normalization from
renormalization or subtraction, defines every reported uncertainty, and cites every selected
method. API-local work remains under `APIDOC-103` and `APIDOC-206`.

=== MANUAL-202 · Spenso physicist-to-type bridge

*Priority:* P0 · *Status:* In progress

The first tutorial uses a dimension-two `Lorentz` representation for a generic scalar
contraction, while the HEP library receives only a brief mention. The network guide
also describes a reusable stored plan more strongly than the documented public
strategy interface appears to support.

*Progress (`umvyopkm`):* the toy contraction now identifies its two-dimensional representation
as Euclidean linear algebra rather than Lorentz physics, and the network guide accurately
describes execution-time strategy selection instead of a cached public plan. The complete
Einstein-notation-to-type bridge and convention-explicit HEP workflow remain open.

*Example finding (2026-08-21):* the authored “network” example computes the affine combination
`2 a + 3 b`; it has no shared dummy index and demonstrates no tensor contraction. The Python
network example checks only result length, so incorrect components or index structure could pass.

*Complete when:* the planning/reuse claim matches the public implementation and the
toy representation is physically honest. A familiar Einstein-notation expression is
mapped through free/dummy indices, variance/duality, representations, slots,
structure/data, storage, contraction, and result structure. A checked network and
gamma/projector or SU(3) workflow states the relevant metric, basis, index-order, and
normalization conventions. The connected multi-tensor contraction checks its free indices and
components and agrees across at least the Rust and Python paths. Link `APIDOC-104`,
`APIDOC-105`, and `APIDOC-204` for item contracts.

=== MANUAL-203 · Idenso user contract versus implementation evidence

*Priority:* P0 · *Status:* In progress

The FORM color/Dirac reference combines mathematical rules, shipped behavior, target
patterns, source mappings, implementation excerpts, historical Typst and search
commands, validation logs, and unresolved uncertainties. A public reader cannot
reliably distinguish current behavior from a proposed implementation contract.

*Progress (`umvyopkm`):* the public page is now a concise shipped-rule contract with explicit
dimension, color-normalization, and partial-normal-form boundaries. Four-dimensional gating and
the implemented gamma-five trace relation are explicit, but the gamma-five definition and epsilon
orientation/index convention remain open. FORM sources, target patterns, historical commands, and
unresolved provenance remain available in a separate classified developer source map. The
beginner metric example is explicitly four-dimensional Minkowski, but it checks only free-index
counts rather than the claimed expression, coefficient, or surviving index.

*Complete when:* every public rule is classified as shipped and scientifically
reviewed, with conventions and checked identities. Target designs, source maps,
historical commands, validation mechanics, and unresolved work move to appropriately
classified developer records without losing provenance. The beginner path uses a
convention-explicit Minkowski/Dirac/color workflow and checks structural equality to an expected
transformed expression, coefficient, and free indices. Nontrivial public identities cite their
mathematical sources next to the claim. Link `APIDOC-004`, `APIDOC-106`, and `APIDOC-205` for
API-local claims.

=== MANUAL-204 · Linnest draw-first progression

*Priority:* P1 · *Status:* Open

The Linnest manual explains archived graph bytes and opaque payloads before the first
useful drawing. Its “minimal” example immediately adds subgraph filtering, DOT export,
metadata styling, and multiple callbacks, then duplicates much of the construction.

*Output and convention findings (2026-08-21):* the HTML path substitutes prose for the rendered
drawing, so a site reader cannot compare copied source with its expected visual result. The
physics-style prose distinguishes graph orientation, fermion flow, and momentum arrows
mechanically, but does not establish incoming/outgoing conventions or teach readers how to avoid
a physically misleading arrow.

*Complete when:* one canonical build or parse → layout → draw example produces a
verified small physics diagram before customization. Physics styling, subgraphs,
layout constraints, and callbacks are introduced progressively; duplicated source is
removed; algorithm-selection advice covers every supported layout; and archive,
plugin, payload, and wire-format details live in an advanced or developer boundary. The HTML
path shows the checked diagram, and an annotated example distinguishes topology direction,
fermion flow, and momentum flow.

=== MANUAL-205 · Linnet physical cuts, flow, and mutation mappings

*Priority:* P1 · *Status:* Open

Linnet explains half-edge storage well but defers source/sink flow and orientation semantics to
API variants. Its cut example checks only the number of returned partitions, not their members,
boundary half-edges, directions, or physical interpretation. Storage-changing mutation is
described without a checked mapping journey.

*Complete when:* an annotated Feynman- or tensor-network example maps physical incoming/outgoing
legs to half-edges, `Flow`, orientation, subgraphs, and returned cut partitions; distinguishes a
graph-separating cut from an on-shell/Cutkosky cut; and verifies the actual cut members rather than
only their count. One extraction, excision, sewing, or contraction journey preserves and checks
the returned old-to-new mappings instead of reusing stale indices.

== Cross-product findings

=== MANUAL-301 · Audience ladder and product handoffs

*Priority:* P1 · *Status:* Open

The desired progression from physicist-facing application to increasingly technical
domain and algorithm crates is implicit rather than stated. Portal tasks can therefore
land on software-mechanical tutorials when the reader expects a scientific outcome.

*Complete when:* all five product overviews plus Linnest state the intended reader,
scientific task, required programming fluency, primary interface, and direct-use versus
GammaLoop boundary. One ecosystem map and the portal task routes follow the audience
contract above.

=== MANUAL-302 · Canonical terminology ownership

*Priority:* P1 · *Status:* Open

Terms including state/workspace, event/event group, graph/hedge/subgraph, tensor
structure/network, rewrite/cooking, topology/canonicalization/evaluation, and
precision/error have overlapping physics and software meanings across products.

*Complete when:* each overloaded term has one canonical owner and definition,
product-specific prose links that definition at first use, and a terminology audit
finds no conflicting current definitions.

== Concrete findings already owned by the strategy register

#table(
  columns: (auto, 1.8fr, 3fr),
  table.header([*Work item*], [*Existing owner*], [*Concrete manifestation retained by this audit*]),
  [Beginner and real-value paths], [`DOC-002` and `DOC-007`], [Reclassify the selected `gg → hhh` contribution as a state/workspace or benchmark example; make the portal's calculation path reach a bounded checked value or distribution; give the events workflow meaningful cuts and bins; complete the Linnet mutation/drawing, Spenso HEP, Idenso Rust, and Vakint Laurent-series outcomes.],
  [Canonical scientific examples], [`DOC-008`], [Classify examples as beginner, validated result, regression, benchmark, experiment, or historical; record setup, conventions, prerequisites/license/backend, canonical source, expected result/tolerance, runtime tier, verification command, and last verified revision. Current Vakint matching, Idenso identity, and Spenso network prose exceeds its compile-only evidence, while Linnest HTML omits the rendered result.],
  [Release and reproducibility narrative], [`DOC-011`], [Fill GammaLoop/Vakint histories, current Spenso/Idenso gaps, and Linnet-Python history; distinguish API, physics/numerical/default, and compatibility changes; pin revision-sensitive dependencies and sources.],
  [Methods and citations], [`DOC-017`], [Place defining sources next to Local Unitarity and subtraction, Vakint backends, Idenso identities, Spenso HEP conventions, and nontrivial graph algorithms; a broad publication feed is not a method map.],
)
