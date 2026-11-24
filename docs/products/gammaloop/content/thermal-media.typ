#import "../../shared.typ": boundary, developer-link, source-link

#let thermal-media = [
= Thermal and dense media

Medium modes evaluate amplitudes in an equilibrium medium at finite temperature or density.
The medium is at rest in the frame of the supplied momenta.

== Scope and limitations

Medium modes generate fixed-order amplitudes only; cross sections are rejected. Connected vacuum
graphs, which give the pressure, are the only amplitudes validated against literature values.
Higher-point amplitudes are supported but remain experimental.

Medium modes require `global.generation.three_dimensional_representations = ["cff"]`.
LTD supports vacuum integrands only; requests containing LTD, including mixed CFF/LTD
generation, are rejected.

#boundary("Infrared safety and thresholds", [
  Thermal zero-mode subtraction is experimental: local finiteness is validated only on
  two-bubble ring graphs, and near soft static modes the evaluation needs higher precision
  (see Local soft subtraction).
  GammaLoop has no thermal resummation, so only infrared-safe
  fixed-order amplitudes can be evaluated. At finite temperature, soft bosonic modes, whose
  occupation grows like `T/E`, eventually require resummation that GammaLoop does not perform.
  Threshold subtraction is also unavailable, so choose kinematics without thresholds, including
  the thermal (Landau-damping) cuts that a medium adds to the vacuum ones.
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

== Local soft subtraction

At finite temperature the Matsubara zero mode of massless bosonic lines without a chemical
potential makes soft regions such as a static gluon between two polarization loops
divergent. GammaLoop subtracts them locally on the default hedge-poset direct-3D route:
each divergent static soft set is assigned to a zero-mode
host, the subgraph that stays hard; its counterterm keeps the zero-frequency mode of the soft
lines and expands the hard part in their momenta. Hosts are `IR` subgraphs, so this needs
`renormalization_prescription.massless_power_divergent = IR`; every other `IR` subgraph
receives a plain UV counterterm. On this route a divergent static soft set without a supported
host is an error that names this setting, so with the default `MUV` prescription such graphs
are rejected at finite temperature, with and without integrated counterterms. With the default
`generate_integrated = true`, every UV counterterm with the mass `runtime.general.m_uv` gets its
integrated counterpart, which converts it to 4D MS-bar, also where a zero-mode subtraction acts
on it, so that results do not depend on `m_uv`. The zero-mode subtractions themselves and the
static UV counterterms described next stay local: their integrated add-back is not generated yet.
The integrated counterparts are spread over the loop momenta of their subgraphs with the width
`runtime.general.renormalization_localization_scale`, which should be of the order of the
temperature. An integrated counterterm is localized on the rest of the graph. Where the
counterterm's numerator raises that rest's loop-energy degree above its propagators', as for a
gluon polarization loop and a second polarization loop joined by two gluons, a part polynomial
in the loop energy is left over. It has no thermal part, and vacuum subtraction removes it
exactly, so such graphs need `global.generation.medium.vacuum_subtraction = true` with
integrated counterterms; without it generation stops with an error that names the counterterm
and asks for vacuum subtraction. With vacuum subtraction such graphs generate like the others.
Where the left-over part multiplies propagators of that rest, vacuum subtraction does not remove
it, and generation stops with an error either way.
`generate_integrated = false` generates the local subtraction alone. The subtraction works with
and without vacuum subtraction. A UV counterterm applied to a
zero-mode subtraction whose subgraph holds static soft lines of that subtraction is, on that
term, a UV counterterm of the static theory: it regulates those lines with the mass
`runtime.general.m_ir`. That setting is a parameter of the static theory's scheme, independent
of `runtime.general.m_uv`: it has its own default, 1000, and does not follow `m_uv`, so set it
explicitly. The choice is made per term: the same
counterterm keeps `m_uv` on the terms without such a zero-mode subtraction, and so do the
host's own counterterm and all other counterterms. Since the add-back of the zero-mode
subtractions and static counterterms is deferred, results can depend on
`runtime.general.m_ir`, and with `generate_integrated = false` also on `runtime.general.m_uv`.
Near a soft static mode the
remainder's rounding error cannot be
estimated by rotations, so in graphs with a zero-mode
subtraction evaluation starts at higher precision when the softest static line is soft
against the temperature
(`runtime.stability.zero_mode_energy_escalation_factor`), or when the softest massless line of
any species, such as a soft quark, is soft against the loop momenta
(`runtime.stability.zero_mode_hierarchy_quad_escalation_factor` for quadruple precision,
`runtime.stability.zero_mode_hierarchy_escalation_factor` for the last stability level). The
construction, its orders and
its known limitations are described in the
#developer-link("uv-renormalization", "uv-renormalization.typ", "UV renormalization architecture").

== Vacuum subtraction and the pressure

With
#link("reference/cli/settings/cli/global/generation/medium/#setting-cli-global-generation-medium-vacuum-subtraction-622cb5e219356a31")[`vacuum_subtraction = true`],
only the medium-dependent part is generated: each distribution weight is replaced by its
difference from the vacuum limit, and UV counterterms covering all loops are omitted,
including when bridges remain in the graph. When a zero-mode host subtracts a static soft mode
(see above), the all-loop counterterm is instead kept as a plain UV operation: it subtracts
the UV divergence of each host's static reduced cograph, whose zero-mode coefficient
enclosing operations treat as a fixed vertex.
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
