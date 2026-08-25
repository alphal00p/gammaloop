#import "../../shared.typ": boundary, callout, source-link

#let physics-scope = [
= Physics scope and Local Unitarity

GammaLoop is alpha-stage research software for perturbative quantum-field-theory
calculations. Its present, demonstrable calculation path is partonic: it constructs
amplitude integrands and uses Local Unitarity forward-scattering representations for
differential cross sections, then evaluates or integrates selected contributions. It is not
currently a general hadron-level event generator or a promise that every process allowed by
a model is numerically supported.

#callout("Intended reader", [
  This chapter is for physicists deciding whether GammaLoop is appropriate for a
  calculation. The command-line workflow requires no Rust knowledge. Read the
  #link("guides/conventions/")[kinematics and normalization contract] before
  interpreting a value, and treat compiled loop-count limits as software bounds rather
  than demonstrated perturbative reach.
])

== What Local Unitarity changes

At fixed perturbative order, a conventional calculation separates real-emission and
virtual contributions even though their infrared singularities cancel only in the sum.
The Local Unitarity representation uses Loop--Tree Duality for each forward-scattering
graph to express those degrees of freedom on common integration variables. Contributions
that approach the same soft or collinear configuration can then cancel locally, before a
Monte Carlo integral relies on cancellations between separately integrated numbers.

The defining #link("https://arxiv.org/abs/2010.01068")[Local Unitarity paper] proves this
representation for processes without initial-state collinear singularities and introduces
local ultraviolet counterterms. The
#link("https://arxiv.org/abs/2203.11038")[raised-propagator and local-renormalization paper]
extends the formalism with distributional cutting rules and a local R-operation treatment.
Those are statements about the method. A method valid at arbitrary perturbative order is
not evidence that this release implements every order, process, initial state, or
renormalization choice.

For one Monte Carlo point, GammaLoop can produce several cut or counterterm contributions in
one or more correlated event groups. The observable layer commits the complete group list as
one statistical sample. Treating the members as independent events would discard the
local-cancellation relation and underestimate correlations.

== Capability envelope in this release

The labels below are deliberately conservative:

- *Supported* means an explicit implementation path was identified; it does not certify
  every process or numerical regime.
- *Conditional or experimental* means feature-, backend-, model-, or workflow-dependent.
- *Unsupported* means the accepted value reaches a placeholder, explicit error, or panic.
- *Unverified* means the source does not justify a public scientific guarantee.

#table(
  columns: (1.05fr, 1fr, 2.8fr),
  table.header([*Area*], [*Classification*], [*Current boundary*]),
  [Calculation objects], [Supported, scoped], [Amplitude integrands and forward-scattering cross-section integrands are generated. Numerical support must still be established for the selected topology and contribution.],
  [Initial states and PDFs], [Partonic only], [The current external kinematics are fixed momenta and helicities, with a PDF factor of one. There is no PDF convolution or factorization scale. The built-in flux path supports exactly one or two incoming particles.],
  [Orders and selectors], [Conditional filters], [`{n}` constrains an amplitude graph's loops or the sum on two cut sides, while `{{n}}` constrains forward-graph loops. Outside the bracketed block, unpowered constraints populate the amplitude filter and powered constraints populate the cross-section filter, whose exponent is currently discarded. Bracketed named orders bound model-resolved additional cut content on the cross-section path. These mechanisms are not interchangeable and do not provide a universal LO/NLO/NNLO mapping.],
  [Models], [Conditional], [Built-in Standard Model, scalar, and scalar-gravity JSON models are available; other JSON models can be imported. UFO import requires the optional Python-backed feature and its loader.],
  [External polarizations], [Scoped], [Automatic sums are implemented for scalar, spinor, and vector states, with Feynman or light-like axial gauge for vectors. Higher-spin polarization replacement is unsupported.],
  [Masses and widths], [Real masses only], [Real masses enter evaluation. Complex masses are unsupported. Width names are model metadata, but no finite-width or complex-mass scheme was verified.],
  [Subtraction and renormalization], [Conditional or experimental], [Local UV settings and several vacuum-integral backends exist. Threshold subtraction is experimental, and the physical meaning of every exposed prescription still needs a method-author review.],
  [Integration components], [Scoped], [Real and imaginary components can be integrated separately. The accepted `both` setting is not implemented.],
  [Observables], [Supported; requires user validation], [Particle and jet quantities, selectors, histograms, and kT-family clustering are available. GammaLoop does not establish the infrared and collinear safety of an arbitrary user-defined observable.],
  [Licenses and tools], [Conditional], [The Rust crates use MIT or Apache-2.0 terms, while Symbolica has separate license conditions. FORM, vacuum-integral backends, UFO import, and drawing tools have independent availability requirements.],
)

== Decide whether a calculation is in scope

Before generating graphs, write down these physics choices outside the run card:

+ Is the requested quantity partonic with fixed external kinematics, or does it require
  PDFs, a factorization scale, beam structure, or more than two incoming particles?
+ Is it an amplitude or a differential cross section, and which exact graph, loop, cut,
  and coupling-order filters define the requested contribution?
+ Does the active model use only supported external spins and real masses? If a width is
  physically essential, do not infer a finite-width scheme from imported metadata.
+ Which spin and color sums or averages, projectors, flux, symmetry factors, and units
  belong in the normalization? Record them using the
  #link("guides/conventions/")[conventions checklist].
+ Is the observable infrared and collinear safe for the selected real/virtual
  combination? This is a physics precondition, not an automatic validation step.
+ What maintained limit, independent implementation, or published number will be used
  to validate the result and its Monte Carlo uncertainty?

#boundary("Method scope versus software evidence", [
  The first-state tutorial proves that a model, generated integrand, settings, and
  persisted state fit together; it does not validate a physical cross section. Use the
  #source-link("examples/cli/epem_a_ddx/LO/epem_a_ddx.toml", label: "partonic cross-section regression card")
  and the scalar-topology targets as concrete implementation sources to inspect. The former
  still uses a powered selector whose exponent is not enforced, so it is not an audited
  perturbative-order or normalization exemplar. Publish a number only after recording its
  contribution definition, normalization, convergence, and an independent comparison.
])
]
