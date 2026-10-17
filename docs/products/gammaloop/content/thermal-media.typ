#import "../../shared.typ": boundary, developer-link, source-link

#let thermal-media = [
= Thermal and dense media

Medium modes evaluate amplitudes in an equilibrium medium at finite temperature or density.
An initial finite-temperature cross-section path covers the restricted probe cuts below.
The medium is at rest in the frame of the supplied momenta.

== Scope and limitations

Connected vacuum graphs, which give the pressure, are the only amplitudes validated against
literature values. Higher-point amplitudes are supported but remain experimental.

The currently enabled thermal `xs` path is restricted to one future-directed timelike probe and simple cuts
through massless scalar particles, with no virtual propagators and no chemical-potential
assignments on the cut particles. Choose a generation loop-momentum basis in which at most
one cut-particle momentum contains both loop and external momentum; the signed radial solver
supports at most one distinct cut energy that scales nonhomogeneously with the radial variable.
A runtime chemical potential of zero does not remove an adjustable assignment. This path remains
experimental. Spacelike probes, massive or spinning
cut particles, raised cuts, cold modes, vacuum subtraction, UV subtraction, threshold subtraction,
and surface-targeted sampling are unavailable for thermal `xs`.
These are temporary validation limits. Thermal crossings and their occupation factors use the
shared cross-section residue, normalization, selector, event, and precision machinery.
The internal representation of signed cuts does not by itself establish support for the
excluded processes.

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

For the restricted thermal `xs` path, disable threshold subtraction at both generation and
runtime, and use ordinary loop-momentum sampling:

// docs-example: syntax
```toml
[cli_settings.global.generation.medium]
mode = "thermodynamic_equilibrium"
vacuum_subtraction = false

[cli_settings.global.generation.uv]
subtract_uv = false

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = false

[default_runtime_settings.subtraction]
disable_threshold_subtraction = true
```

The experimental scalar sunset #source-link("examples/cli/thermal_rates/scalar_sunset/run.toml", label: "run card")
and #source-link("examples/cli/thermal_rates/scalar_sunset/benchmark.typ", label: "benchmark notes")
provide a concrete example and an independent phase-space reference. The reference alone does
not validate native generation or integration.

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
only the medium-dependent amplitude is generated: each distribution weight is replaced by its
difference from the vacuum limit, and UV counterterms covering all loops are omitted,
including when bridges remain in the graph.
This option requires a medium mode.

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
]
