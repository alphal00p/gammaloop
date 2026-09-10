#import "../../shared.typ": boundary, callout, source-link

#let conventions = [
= Kinematics, normalization, and weights

This page assembles the source-verifiable scientific conventions needed to interpret a
GammaLoop sample. It intentionally separates factors supplied by the runtime from factors
that remain the physicist's responsibility. Where the implementation does not establish a
universal convention, the text says so instead of filling the gap from customary notation.

== Four-vectors and momentum flow

External four-vectors use component order $(E, p_x, p_y, p_z)$ and metric signature
$g = "diag"(+1, -1, -1, -1)$, so

$ p · q = E_p E_q - (p_x q_x + p_y q_y + p_z q_z), quad p^2 = E^2 - p_x^2 - p_y^2 - p_z^2. $

Stored external momenta are physical, normally positive-energy vectors. A separate process
signature marks a particle as initial or final; GammaLoop does not require the user-facing
list to store every final momentum with a minus sign. When one momentum is dependent, it is
reconstructed to satisfy

$ sum_i^"initial" p_i - sum_j^"final" p_j = 0. $

The runtime field `e_cm` is a construction and normalization scale. If it is omitted for
fixed momenta, the current fallback is the mean absolute value of their nonzero components,
not a physical reconstruction of $sqrt(s)$. Supply `e_cm` explicitly whenever its physical
meaning matters.

== Helicity and polarization

Fixed external helicities use `-1`, `0`, or `+1`. Zero is the physical and documented
convention for a scalar card, but the current evaluator ignores rather than enforces a scalar's
helicity. The automatic external-momentum generator currently assigns `+1` to every external
entry; that is a fixed-polarization default, not an unpolarized calculation.

Use `summed` for a polarization sum and `summed_averaged` for the same replacement with the
implemented spin-state average. The setting acts on every external entry explicitly marked
that way; incoming particles are not averaged merely because they are incoming. Implemented
averaging factors are 1 for a scalar, $1 / 2$ for a spinor, $1 / 2$ for a massless vector,
and $1 / 3$ for a massive vector. Automatic higher-spin polarization sums are unsupported.

#callout("Color is a separate normalization decision", [
  No general automatic incoming-color average was identified. A maintained calculation may
  supply a color projector or global prefactor explicitly. Record that choice rather than
  assuming that `summed_averaged` includes color.
])

== Flux and reporting units

With the flux factor enabled, the implemented one-incoming-particle factor is

$ F_"1 in" = 1 / (2 E), $

and the implemented two-incoming-particle factor is

$ F_"2 in" = 1 / (4 sqrt((p_1 · p_2)^2 - m_1^2 m_2^2)). $

Zero or three-or-more incoming particles are not implemented by this default path. Reporting
units resolve to `pb`, `fb`, `ab`, `mb`, or `none`; `auto` selects picobarns for more than one
incoming particle and no conversion otherwise. The fixed factor
$3.89379372171859 times 10^8$ relative to picobarns is applied only by the two-incoming-particle
flux branch. The one-incoming-particle branch returns $F_"1 in"$ without a reporting-unit
conversion. Disabling the flux factor also bypasses the conversion branch.

The code describes input scales as model energy units, while the fixed picobarn conversion
is the usual natural-unit form associated with GeV. Until the energy-unit premise is made an
explicit model contract, do not claim that arbitrary model units convert correctly to pb.

== From an integrand value to a Monte Carlo sample

For the public sample-evaluation contract, the exact weighting relation is

$ w_"sample" = I_"returned" J_"parameterization" w_"MC". $

`integrand_result` is the returned integrand before the parameterization Jacobian, and
`integrator_weight` excludes that Jacobian. During integration, GammaLoop first multiplies by
the parameterization Jacobian and then applies the Monte Carlo sample weight. For a cross
section, the returned integrand already contains the enabled flux and reporting-unit factor.

An event group may retain an original contribution, threshold counterterms, and one full
multiplicative factor. Its implemented reconstruction is

$ w_"event" = (w_"Original" + sum_i w_"CT,i") w_"FullMultiplicativeFactor". $

All groups originating from one Monte Carlo point are combined before one statistical sample
is committed to a histogram bin. They are correlated contributions, not independent events.
The #link("guides/events-and-observables/")[events and observables guide] shows the corresponding
Python and Rust records.

== Normalization ownership checklist

#table(
  columns: (1.35fr, 1.15fr, 2.4fr),
  table.header([*Factor*], [*Usual owner here*], [*What to record*]),
  [Polarization sum and spin average], [Run card], [The explicit helicity mode for every external particle and the vector gauge choice.],
  [Color sum or average], [Physicist / projector], [The projector or prefactor and its normalization; do not infer a universal automatic average.],
  [Graph multiplicity and signs], [Generated graph factor], [Automorphism, fermion-loop, external-ordering, and grouping contributions visible in the generated state.],
  [Flux and output units], [Runtime when enabled], [Incoming momenta and masses, the selected unit, and whether the flux branch was disabled.],
  [Parameterization Jacobian], [Integration driver], [The parameterization and whether a reported value is before or after its Jacobian.],
  [Monte Carlo weight], [Integrator], [The sample weight, seed, component, and correlation boundary.],
  [Subtraction contributions], [Correlated event group], [Original and counterterm weights, signs, and the common multiplicative factor.],
  [Overall physical normalization], [Calculation definition], [Any remaining coupling, projector, color, symmetry, phase-space, or conventional factor used for comparison.],
)

Generated graph factors retain several combinatorial and fermionic contributions, but the
source audit did not establish one universal statement that these plus a chosen projector form
the complete physical symmetry and averaging factor for every process. Inspect the generated
factor and write the full normalization used by the comparison calculation.

== Scales and scheme names

Keep these settings distinct:

- `mu_r` is the renormalization scale and enters the parameter set as its square;
- `m_uv` is a UV reference mass;
- `renormalization_localization_scale` controls localization of local UV counterterms;
- MUV, PolePart, OS, IR, Unsubtracted, and VaccuumLimit are accepted local-counterterm
  prescription names whose full physical definitions still require a method-author contract;
  OS, IR, and VaccuumLimit reach explicit unimplemented paths in at least one current
  counterterm evaluator; and
- Vakint's default `MSbar` choice describes its vacuum-integral normalization and must not be
  used as a synonym for GammaLoop's MUV prescription.

No factorization scale belongs to the current fixed-partonic external-kinematics path.

== Terms that are easy to overload

- A *GammaLoop state* is the persisted calculation workspace, not a quantum state and not an
  integration checkpoint.
- An *integration workspace* owns resumable Monte Carlo iterations for selected
  `(process, integrand)` slots; its summary JSON is not the authoritative checkpoint.
- An *integration coordinate* is a parameterization coordinate unless `momentum_space` is
  explicitly selected, in which case each loop momentum is a flattened `(p_x, p_y, p_z)`
  triplet.
- A *target* is a comparison value used to report a delta. It is not the relative or absolute
  accuracy that stops integration.
- An *event group* contains correlated contributions from one Monte Carlo point. A *graph
  group* is a generation-time grouping of related graphs and multiplicity/numerator factors;
  the two are not interchangeable.

#boundary("Conventions not yet established as universal contracts", [
  Before publishing a result, obtain an explicit calculation-level statement for all
  $2 pi$, $i$, sign, contour, and phase-space-measure factors; overall amplitude and
  cross-section normalization; color sums and averages; remaining symmetry factors;
  powered coupling constraints; the LO/NLO meaning of graph filters; finite-width treatment;
  and the selected renormalization prescription. The current source does not justify silently
  supplying any of these from convention alone.
])

The
#source-link("examples/cli/epem_a_ddx/LO/epem_a_ddx.toml", label: "e+ e- partonic regression card")
is useful for tracing implemented factors, but it is not yet a normalization exemplar: its
powered selector exponent is not enforced. For any audit, record the model, contribution
filters, helicities, projectors, runtime settings, observable definition, and independently
reviewed perturbative-order definition.
]
