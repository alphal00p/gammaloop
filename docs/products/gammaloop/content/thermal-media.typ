#import "../../shared.typ": boundary, developer-link, source-link

#let thermal-media = [
= Thermal and dense media

Medium modes evaluate amplitudes in an equilibrium medium at finite temperature or density.
The medium is at rest in the frame of the supplied momenta.

== Scope and limitations

Medium modes generate fixed-order amplitudes only; cross sections are rejected. Connected vacuum
graphs, which give the pressure, are the only amplitudes validated against literature values.
Higher-point amplitudes are supported but remain experimental.

#boundary("Infrared safety and thresholds", [
  GammaLoop has no thermal infrared subtraction or resummation, so only infrared-safe
  fixed-order amplitudes can be evaluated. At finite temperature, soft bosonic modes, whose
  occupation grows like `T/E`, eventually require resummation that GammaLoop does not perform.
  Threshold subtraction is also unavailable, so choose kinematics without thresholds, including
  the thermal (Landau-damping) cuts that a medium adds to the vacuum ones.
])

The following are also unsupported:

- local UV counterterms from expanded 4D integrands
  (#link("reference/cli/settings/cli/global/generation/uv/#setting-cli-global-generation-uv-local-uv-cts-from-expanded-4d-integrands-7def650172e108c5")[`local_uv_cts_from_expanded_4d_integrands`]);
  medium modes use 3D local subtraction;
- external momenta that make otherwise distinct propagator poles coincide, such as zero
  four-momentum transfer between equal-mass propagators, which can leave unresolved `0/0`
  thermal factors (repeated poles found during generation use the derivative treatment below);
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

== Chemical potentials

Chemical potentials are model parameters. Each particle may name one (`chemical_potential` in
JSON and UFO models); a particle without one has zero chemical potential. UFO import keeps these
assignments and therefore requires `ufo-model-loader>=1.0.0`. In the bundled `sm` model,
antiparticles use the negated `minus_<name>` parameter, and every particle's value is derived
from the external parameters `muB`, `muQ`, `muLe`, `muLmu`, and `muLtau` (LHA block
`CHEMICALPOTENTIAL`) using its baryon number, electric charge, and lepton flavor. Assignments use
explicit particle metadata, not the `mu` prefix: `muH` belongs to the Higgs potential. The other
bundled models declare no chemical potentials.

With the default `--simplify-model=true`, a restriction card turns every parameter it sets to
zero into a constant. `sm-default` therefore fixes all five chemical potentials at zero. Import
`sm-thermal` instead, which sets `muB = 3` and keeps it adjustable, or `sm-full` when other
chemical potentials must vary. Then set the values, for example with `set model muB=0.5`.

Finite-temperature evaluation requires a finite, positive inverse temperature. Chemical
potentials must be finite and real. Bosons require a finite, nonnegative real mass and
`|mu| < m`, with the massless `m = mu = 0` case also supported. Massive saturation
(`|mu| = m > 0`) and boson condensation are unsupported. Fermions may have `|mu| >= m`; at zero
temperature, a fermion with a chemical potential needs a finite real mass. These
checks use the current model values and any DOT mass overrides when the integrand is warmed up,
including after model updates or state reloads.

== Vacuum subtraction and the pressure

With
#link("reference/cli/settings/cli/global/generation/medium/#setting-cli-global-generation-medium-vacuum-subtraction-622cb5e219356a31")[`vacuum_subtraction = true`],
only the medium-dependent part is generated: each distribution weight is replaced by its
difference from the vacuum limit, and the overall UV counterterm of the full graph is omitted.
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

Several equal-mass propagators carrying the same loop momentum, as in self-energy insertions,
produce energy derivatives of a distribution. In `zero_temperature_equilibrium`, those
derivatives of a fermion with a nonzero chemical potential are delta functions on its Fermi
surface. Bosons and fermions at zero chemical potential have no Fermi surface, so their
derivatives vanish. GammaLoop localizes each delta function by rescaling the loop momentum and
inserting the normalized positive-scale profile configured in
#link("reference/cli/settings/runtime/h-function/")[`runtime.h_function`]. The profile needs a
positive `sigma` and a supported `power`; `exponential_ct` is rejected. Integrated results do not
depend on that profile, but pointwise values do, so fix it explicitly when comparing sampled
values. A Fermi surface that shrinks to a point, `|mu| = m > 0`, is unsupported. The
#source-link("tests/resources/run_cards/fermi_surface_1l_integration.toml", label: "one-loop Fermi-surface run card")
integrates such a contribution.
]
