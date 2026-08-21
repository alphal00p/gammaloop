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
2. Vakint, Idenso, and Spenso start from mathematical or HEP notation, then explain
   the computational representation and domain workflow.
3. Linnet and Linnest may foreground graph models, invariants, algorithms, and
   programmable interfaces after showing the scientific operation they support.
4. Native Rustdoc and developer records own detailed types, storage, wire formats,
   logging internals, complexity, and implementation architecture.

Every overview should state its intended reader, primary task and interface,
expected programming fluency, and “use this directly when / use GammaLoop when”
boundary. A GammaLoop reader should not need to infer the crate architecture to
complete or assess a calculation.

== Audit summary

#table(
  columns: (1fr, 1.8fr, 2.8fr),
  table.header([*Product*], [*Content to preserve*], [*Critical remaining gap*]),
  [GammaLoop], [Careful state, persistence, command, and operational boundaries.], [The method, capability envelope, conventions, first result, uncertainty, and validation are missing from the primary physicist path.],
  [Vakint], [Clear matching, reduction, evaluation, backend, and citation stages.], [Software notation is not mapped to a mathematical integral; normalization and checked Laurent coefficients are absent.],
  [Idenso], [Deep source-backed rule material and strong rewrite-stage cautions.], [Shipped scientific behavior, target specifications, implementation evidence, and historical validation are mixed in one public reference.],
  [Spenso], [Strong structure/data/execution model and useful failure ownership.], [The bridge from indexed physics to representations, slots, duality, networks, and HEP conventions is missing.],
  [Linnet], [The strongest progressive conceptual introduction, including the physical value of half-edge boundaries.], [The maintained strategy still needs the checked mutation and complete graph-result journey recorded under `DOC-007`.],
  [Linnest], [Broad graph, layout, drawing, and physics-style capability.], [Archive/plugin payload mechanics and advanced callbacks precede the smallest useful drawing result.],
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

*Priority:* P1 · *Status:* Open

The overview names Local Unitarity and immediately switches to software subsystems;
the process guide gives only one sentence about forward-scattering graphs and cuts.
The manual does not let a physicist determine whether a proposed calculation is
partonic or hadronic, supported or experimental, or dependent on an unavailable
model, contribution, backend, tool, or license.

*Complete when:* a physicist-facing chapter explains the forward-graph/cut relation,
real/virtual combination, common integration variables, local cancellation, the role
of correlated events and UV/IR/threshold subtraction, and applicable observable
conditions. A maintained matrix classifies supported, experimental, and unsupported
initial states, orders/contributions, models, masses/spins, PDFs, scales and schemes,
observables, external tools, and license requirements, with defining method citations.

=== MANUAL-102 · Scientific conventions and GammaLoop glossary

*Priority:* P1 · *Status:* Open

Settings expose convention-sensitive choices, but no authored chapter assembles the
scientific contract. Terms such as state, event group, graph group, raised cut,
orientation, loop-momentum basis, E-surface, integration coordinate, and target are
introduced without one canonical physical definition.

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
criteria. Link `APIDOC-201` and `APIDOC-202` for Python result/event contracts.

=== MANUAL-104 · Public troubleshooting and companion-tool boundary

*Priority:* P1 · *Status:* Open

GammaLoop's public diagnostics chapter becomes a contributor tracing-DSL manual with
callsite metadata, filter implementation, routing prefixes, and tag-naming guidance.
The GammaLoop manual also embeds the complete Kurvst geometry API despite drawing
being a companion workflow.

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

*Complete when:* a representative expression is mapped argument by argument to its
integral, including dimension/epsilon, metric, loop measure, (+i0), mass and scale,
momentum routing, propagator power, and external-momentum conventions. Topology
entries expose denominators or diagrams, routing, masses, active propagators, backend
support, and epsilon limits. One checked Laurent series records coefficients,
tolerance, backend/tool versions, and normalization. API-local work remains under
`APIDOC-103` and `APIDOC-206`.

=== MANUAL-202 · Spenso physicist-to-type bridge

*Priority:* P0 · *Status:* Open

The first tutorial uses a dimension-two `Lorentz` representation for a generic scalar
contraction, while the HEP library receives only a brief mention. The network guide
also describes a reusable stored plan more strongly than the documented public
strategy interface appears to support.

*Complete when:* the planning/reuse claim matches the public implementation and the
toy representation is physically honest. A familiar Einstein-notation expression is
mapped through free/dummy indices, variance/duality, representations, slots,
structure/data, storage, contraction, and result structure. A checked network and
gamma/projector or SU(3) workflow states the relevant metric, basis, index-order, and
normalization conventions. Link `APIDOC-104`, `APIDOC-105`, and `APIDOC-204` for item
contracts.

=== MANUAL-203 · Idenso user contract versus implementation evidence

*Priority:* P0 · *Status:* Open

The FORM color/Dirac reference combines mathematical rules, shipped behavior, target
patterns, source mappings, implementation excerpts, historical Typst and search
commands, validation logs, and unresolved uncertainties. A public reader cannot
reliably distinguish current behavior from a proposed implementation contract.

*Complete when:* every public rule is classified as shipped and scientifically
reviewed, with conventions and checked identities. Target designs, source maps,
historical commands, validation mechanics, and unresolved work move to appropriately
classified developer records without losing provenance. The beginner path uses a
convention-explicit Minkowski/Dirac/color workflow with expected identity and free
indices. Link `APIDOC-004`, `APIDOC-106`, and `APIDOC-205` for API-local claims.

=== MANUAL-204 · Linnest draw-first progression

*Priority:* P1 · *Status:* Open

The Linnest manual explains archived graph bytes and opaque payloads before the first
useful drawing. Its “minimal” example immediately adds subgraph filtering, DOT export,
metadata styling, and multiple callbacks, then duplicates much of the construction.

*Complete when:* one canonical build or parse → layout → draw example produces a
verified small physics diagram before customization. Physics styling, subgraphs,
layout constraints, and callbacks are introduced progressively; duplicated source is
removed; algorithm-selection advice covers every supported layout; and archive,
plugin, payload, and wire-format details live in an advanced or developer boundary.

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
  [Canonical scientific examples], [`DOC-008`], [Classify examples as beginner, validated result, regression, benchmark, experiment, or historical; record setup, conventions, prerequisites/license/backend, canonical source, expected result/tolerance, runtime tier, verification command, and last verified revision.],
  [Release and reproducibility narrative], [`DOC-011`], [Fill GammaLoop/Vakint histories, current Spenso/Idenso gaps, and Linnet-Python history; distinguish API, physics/numerical/default, and compatibility changes; pin revision-sensitive dependencies and sources.],
  [Methods and citations], [`DOC-017`], [Place defining sources next to Local Unitarity and subtraction, Vakint backends, Idenso identities, Spenso HEP conventions, and nontrivial graph algorithms; a broad publication feed is not a method map.],
)

