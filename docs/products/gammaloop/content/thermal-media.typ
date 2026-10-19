#import "../../shared.typ": boundary, developer-link, source-link

#let thermal-media = [
= Thermal and dense media

Medium modes evaluate amplitudes in an equilibrium medium at finite temperature or density.
The experimental cross-section path evaluates fully inclusive thermal rates.
The medium is at rest in the frame of the supplied momenta.

== Scope and limitations

Connected vacuum graphs, which give the pressure, are the only amplitudes validated against
literature values. Higher-point amplitudes are supported but remain experimental.

Thermal `xs` reuses the ordinary inclusive cross-section machinery for particles with spin,
masses, chemical potentials, virtual propagators, and raised cuts. There is no separate
scalar-only or one-probe whitelist. Existing vacuum kinematic capabilities, model-domain
checks, and the Fermi-surface limitations below still apply. Removing a scope guard does
not establish correctness or finiteness for every process. A spacelike spectral kernel is
not an on-shell particle’s lifetime.

Local 3D UV subtraction and integrated finite counterterms are shared with vacuum rates.
Their physical renormalization-scale dependence is retained. Subtracting the vacuum rate
alone does not renormalize vacuum subgraphs multiplying medium-dependent outer loops.
Thermal threshold subtraction must remain disabled at generation and runtime; choose
infrared-finite inclusive observables whose complete cut sum cancels the relevant internal
singularities, or kinematics without such singularities. This does not make an isolated
singular cut integrable.

The radial solver retains the inexpensive positive-energy interior chart. General signed
or exterior charts use directed value and derivative enclosures over the compactified
positive ray, without a fixed energy-count cap. They run in binary64 first and again at
2048 bits whenever binary64 cannot decide, so binary64 uncertainty never closes a channel.
A certified closed channel contributes zero;
unresolved tangencies, cusps, and exhausted work or precision budgets are reported explicitly.
Original affine routing is retained throughout. Positive-energy surface sampling is available
when its host is strictly interior; the complete signed rate is evaluated at the mapped point.
Mixed-sign surface targeting is not implemented. Compact maps require a full-support companion.

Optional vacuum subtraction uses the same distribution-weight transformation as medium amplitudes.
Thermal crossings and their occupation factors pass through the shared residue, normalization,
selector, event, and precision machinery. The benchmark coverage described below is narrower
than the shared architecture’s intended inclusive scope.

#boundary("Infrared safety and thresholds", [
  GammaLoop has no thermal infrared subtraction or resummation, so fixed-order amplitudes
  and rates must be infrared finite. At finite temperature, soft bosonic modes, whose
  occupation grows like `T/E`, eventually require resummation that GammaLoop does not perform.
  For amplitudes, threshold subtraction is also unavailable, so choose kinematics without
  thresholds, including the thermal (Landau-damping) cuts that a medium adds to the vacuum ones.
])

The following are also unsupported:

- local UV counterterms from expanded 4D integrands
  (#link("reference/cli/settings/cli/global/generation/uv/#setting-cli-global-generation-uv-local-uv-cts-from-expanded-4d-integrands-7def650172e108c5")[`local_uv_cts_from_expanded_4d_integrands`]);
  medium modes use 3D local subtraction;
- external momenta or distinct mass parameters with values equal up to sign that make otherwise
  distinct propagator poles coincide, leaving unresolved `0/0` thermal factors. Generation uses the
  derivative treatment below only when propagators share one mass parameter and the same
  momentum routing, including external shifts;
- at zero temperature, a Fermi-surface term in which another fermion's occupation step
  depends on the localized loop momentum, so that the two Fermi surfaces intersect; generation
  reports it;
- standalone evaluator export of processes with Fermi-surface localization; keep the saved
  process state instead;
- evaluating medium expressions with the standalone evaluator of the diagnostic `3Drep`
  command.

== Select the medium

Set
#link("reference/cli/settings/cli/global/generation/medium/#setting-cli-global-generation-medium-mode-aa98d46a871caf94")[`global.generation.medium.mode`]
before `generate`:

- `thermodynamic_equilibrium` evaluates the distributions at the inverse temperature `1/T` given
  by
  #link("reference/cli/settings/runtime/general/#setting-runtime-general-inverse-temperature-fce360cba2a48bb2")[`runtime.general.inverse_temperature`]
  (default `1`, in inverse model energy units);
- `zero_temperature_equilibrium` takes their `T -> 0` limit: the inverse temperature is unused,
  and the distributions become step functions at the chemical potentials.

The mode is fixed at generation; changing it requires regenerating the integrand. The inverse
temperature and chemical potentials are runtime values that can change per integrand without
regeneration.
#link("reference/cli/settings/cli/global/generation/threshold-subtraction/#setting-cli-global-generation-threshold-subtraction-enable-thresholds-39fd8091d536c278")[`enable_thresholds`]
defaults to `true` and is rejected in medium modes, so disable it explicitly:

// docs-example: syntax
```toml
[cli_settings.global.generation.medium]
mode = "thermodynamic_equilibrium"
vacuum_subtraction = true

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = false

[default_runtime_settings.general]
inverse_temperature = 10.0

[default_runtime_settings.model]
muB = [0.5, 0.0]
```

For thermal `xs`, disable threshold subtraction at both generation and runtime:

// docs-example: syntax
```toml
[cli_settings.global.generation.medium]
mode = "thermodynamic_equilibrium"
vacuum_subtraction = false

[cli_settings.global.generation.uv]
subtract_uv = true
generate_integrated = true
local_uv_cts_from_expanded_4d_integrands = false

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = false

[default_runtime_settings.subtraction]
disable_threshold_subtraction = true
```

The experimental scalar sunset #source-link("examples/cli/thermal_rates/scalar_sunset/run.toml", label: "run card")
and #source-link("examples/cli/thermal_rates/scalar_sunset/benchmark.typ", label: "benchmark notes")
provide a concrete example and an independent phase-space reference. The reference alone does
not validate native generation or integration.
The #source-link("examples/cli/thermal_rates/scalar_bubble/run.toml", label: "massive scalar bubble card")
and #source-link("examples/cli/thermal_rates/scalar_bubble/benchmark.typ", label: "bubble benchmark")
cover decay, absorption, a closed channel, spacelike support, and raised-propagator mass
derivatives. The benchmark gives the normalization conversions for
#link("https://arxiv.org/abs/1503.03534")[Ho and Scherrer (2015)] and
#link("https://arxiv.org/abs/hep-ph/0302098")[Nishikawa, Morimatsu, and Hidaka (2003)].

The restored #source-link("examples/cli/thermal_rates/massive_double_triangle/run.toml", label: "massive double-triangle card")
checks inclusive virtual and real cuts below and above internal thresholds. The
#source-link("examples/cli/thermal_rates/massive_self_energy/run.toml", label: "raised scalar insertion card")
includes ordinary local and integrated UV subtraction and an auxiliary-scale control.
Their benchmark notes explicitly separate historical native integrations from current
acceptance runs. The #source-link("examples/cli/thermal_rates/qcd_vector/run.toml", label: "QCD vector-current card")
exercises the shared spinning and gauge-theory path. All 116 selected native regressions
pass, including restored scalar cut/event closure, UV-addback controls, saved-state reload,
and QCD fixtures. Current integrated scalar and QCD comparisons agree with their
independent references within the reported Monte Carlo errors. Precision varies by case:
the UV-renormalized scalar insertion has 2.50% full-rate error, while the above-threshold
double-triangle inclusive error is 28.8% and its smaller virtual mean remains unresolved
relative to its own error. Direct pole controls separately check its local cancellation.
The benchmark records preserve these limits and distinguish rescued precision attempts
from final discarded samples.

== Chemical potentials

Chemical potentials are model parameters. Each particle may name one (`chemical_potential` in
JSON and UFO models); a missing assignment or the literal `ZERO` means zero chemical potential.
Restriction-fixed chemical potentials are canonicalized to `ZERO`, including derived expressions.
A runtime value of zero does not erase an adjustable assignment. UFO import keeps these
assignments and therefore requires `ufo-model-loader>=1.0.0`. In the bundled `sm` model,
antiparticles use the negated `minus_<name>` parameter, and every particle's value is derived
from the external parameters `muB`, `muQ`, `muLe`, `muLmu`, and `muLtau` (LHA block
`CHEMICALPOTENTIAL`) using its baryon number, electric charge, and lepton flavor. Assignments use
explicit particle metadata, not the `mu` prefix: `muH` belongs to the Higgs potential. The other
bundled models declare no chemical potentials.

For medium evaluation, a particle and its antiparticle must have equal and opposite
chemical-potential values; a self-conjugate particle must have zero chemical potential.
Warm-up checks these relations at the current model values, including after parameter updates
or state reloads. Independent external parameters are allowed when their values satisfy the
relations.

With the default `--simplify-model=true`, a restriction card turns every parameter it sets to
zero into a constant. `sm-default` therefore fixes all five chemical potentials at zero. Import
`sm-thermal` instead, which sets `muB = 3` and keeps it adjustable. To vary other chemical
potentials, create a custom restriction card based on `thermal` and give the desired parameters
nonzero initial values. Keep `MU = MD = MS = ymu = ymd = yms = 0.0` in that card when the
light quarks should remain massless with vanishing Yukawa couplings.
`sm-thermal_s_mass` matches `sm-thermal` but also keeps the
strange-quark mass `MS` and Yukawa input `yms` adjustable, retaining the corresponding Higgs
and Goldstone interactions. After import, set the desired runtime values, including zero,
for example with `set model muB=0.5`.

Finite-temperature evaluation requires a finite, positive inverse temperature. Chemical
potentials must be finite and real. Bosons require a finite, nonnegative real mass and
`|mu| < m`, with the massless `m = mu = 0` case also supported. Massive saturation
(`|mu| = m > 0`) and boson condensation are unsupported. Fermions may have `|mu| >= m`.
Finite-temperature fermions require a finite real mass; at zero temperature the same check
applies when a fermion has a chemical potential and therefore an occupation step. These
checks use the current model values and any DOT mass overrides when the integrand is warmed up,
including after model updates or state reloads.

Amplitude warm-up also checks mixed-sign denominators in the selected generated orientations
and rejects those that vanish for every loop momentum at the current masses and external
momenta. For coincident poles with the same routing, regenerate
with a shared mass parameter so generation can use the derivative treatment. The check uses
exact equality; nearly coincident poles remain a numerical issue.

== Vacuum subtraction and the pressure

With
#link("reference/cli/settings/cli/global/generation/medium/#setting-cli-global-generation-medium-vacuum-subtraction-622cb5e219356a31")[`vacuum_subtraction = true`],
only the medium-dependent contribution is generated: each distribution weight is replaced by its
difference from the vacuum limit, and UV counterterms covering all loops are omitted,
including when bridges remain in the graph.
This option requires a medium mode.

For supported thermal `xs` processes, the result is the excess spectral kernel
$Gamma_rho(P) - Gamma_("vac")(P)$. The scalar sunset card calls this process `thermal_excess`,
alongside the full `equilibrium` kernel and the `vacuum` control. To obtain the production
source for an initially empty bosonic probe through KMS, use the full kernel:
$S(P) = n_B(E) Gamma_rho(P)$. Add the vacuum contribution back before applying this relation
if starting from `thermal_excess`.

GammaLoop keeps equilibrium vacuum amplitudes in its Minkowski convention. For the complete
connected vacuum-graph sum $A$, including model factors and counterterms, convert to the pressure
contribution by hand: $delta p = -i A$. This holds at every loop order and also after vacuum
subtraction. The
#developer-link("phase-conventions", "phase-conventions.typ", "phase-convention derivation")
follows from Wick rotating the single overall spacetime-volume factor; stripped integrals and
arbitrary numerator replacements require their own conversion. Keep
`integrator.integrated_phase = "imag"` and comparison targets in the Minkowski convention; the
real pressure contribution is then the reported imaginary part of $A$, with its sign and
uncertainty preserved. The maintained
#source-link("examples/cli/eos/cool_qm/NLO/cool_qm_eos_NLO.toml", label: "NLO cold quark matter run card")
integrates such a contribution.

== Fermi surfaces at zero temperature

Several propagators sharing one mass parameter and the same momentum up to an overall sign,
including external shifts, produce energy derivatives of a distribution, as in self-energy
insertions. In `zero_temperature_equilibrium`, those
derivatives of a fermion with a nonzero chemical potential are delta functions on its Fermi
surface. Bosons and fermions at zero chemical potential have no Fermi surface, so their
derivatives vanish. The pointwise occupation is undefined (NaN) exactly on the Fermi surface,
`E = sigma*mu`; deterministic inspect points should lie off that discontinuity. Localized
delta terms are evaluated separately and do not assign a pointwise value to `theta(0)`.
GammaLoop localizes each delta function by rescaling the loop momentum and
inserting the normalized positive-scale profile configured in
#link("reference/cli/settings/runtime/h-function/")[`runtime.h_function`]. The profile needs a
positive `sigma` and a supported `power`; `exponential_ct` is rejected. Integrated results do not
depend on that profile, but pointwise values do, so fix it explicitly when comparing sampled
values. Positive-order Fermi-surface terms at the degenerate onset
`sigma*mu = |m| > 0` are unsupported unless their product already has empty support.
This is checked when the selected term is evaluated, at that evaluation's precision;
ordinary order-zero occupation steps are not rejected by this check. The
#source-link("tests/resources/run_cards/fermi_surface_1l_integration.toml", label: "one-loop Fermi-surface run card")
illustrates regular-shell localization.

Rates use the same Fermi-surface sector extraction and coefficient evaluators as amplitudes.
Independent virtual spectator shells retain the existing derivative localization. A cut
depending on the shell momentum is currently supported only for one regular Fermi delta
and a two-energy cut with a simple LU coefficient. That chart fixes the shell radius and
localizes the cut in its polar angle; its prepared event passes through the usual selectors
and event-weight accounting. Higher derivatives of a delta on a dependent cut, multiple
dependent Fermi shells, tangent intersections, and a vanishing angular axis are rejected.
For a double LU pole, generation also retains the boundary delta produced when its radial
derivative acts on an order-zero Fermi occupation. This symbolic boundary is passed to the
same simple-cut localizer; differentiating the step numerically would omit it. Cold occupation
terms at LU pole order above two and pre-existing Fermi distributions inside a raised LU
coefficient are rejected until their higher boundary intersections are implemented.
This does not provide thermal threshold or infrared subtraction. The finite-density
boundary terms and their pole-derivative origin are discussed by
#link("https://arxiv.org/abs/2208.14479")[Gorda, Österman, and Säppi (2022)].
]
