= Thermal rates development handoff

== Current milestone

Branch: `thermal_rates`, based on `ff533fcd9` from main. The implementation was rebuilt
against this base; the old `thermal_rate_investigation` branch is background material,
not a patch source. The first supported example is the scalar sunset with simple thermal
cuts, using the existing signed CFF representations and thermal weights.

#link("scalar_sunset/run.toml")[The run card],
#link("scalar_sunset/model.json")[the portable scalar model],
#link("scalar_sunset/reference.py")[the independent phase-space reference], and
#link("scalar_sunset/benchmark.typ")[the benchmark notes] are the entry points.
The native observable is the spectral relaxation kernel
$Gamma_rho(P) = -op("Im") Pi_R(P)/E$; the existing one-particle flux supplies
$1/(2E)$. Optional vacuum subtraction uses the shared thermal weights and returns
$Gamma_rho - Gamma_("vac")$. The card names this process `thermal_excess`; `equilibrium`
retains the full kernel and `vacuum` supplies the control. The KMS production-source conversion
requires the full $Gamma_rho$, as explained in the benchmark.

== Shared design and code owners

- `crates/three-dimensional-reps/src/expression.rs`: `RaisedSurfaceGroup` carries
  E/H aliases and effective pole order. `select_surface_residues` owns common alias
  normalization, numerator cancellation, denominator-sign handling, and residue extraction.
  Thermal weights remain attached to the expression.
- `crates/gammalooprs/src/cff/mod.rs`: signed equations are grouped after raised-edge
  normalization. Physical cut support is selected before the common CFF residue operation.
- `crates/gammalooprs/src/processes/cross_section.rs`: cut groups include matching physical
  support from all aliases, without duplicates; E and H groups are collected before building
  derivative shapes and evaluators. Zero-order E groups retain the physical-cut partition.
- `crates/gammalooprs/src/cff/radial.rs`: `EnergySurfaceRay` retains ordered signed energy
  occurrences, affine momenta, masses, and the external shift. It owns native/enclosed/dual
  evaluation, precision transfer, and shared root refinement. Algebraic collapsing is only
  a root-isolation view; residue jets retain the represented occurrences.
- `crates/gammalooprs/src/cff/esurface.rs` and `hsurface.rs`: positive-E and signed-H
  interfaces construct the common ray and retain their respective root-isolation policies.
  Preserve `Esurface` as positive-energy geometry: its threshold and sampling algorithms
  depend on that invariant.
- `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`: each localized root
  enters the common LU normalization, residue, event, selector, and precision pipeline.
  Absorbed bath particles become incoming particles before selectors run; simple cuts use
  the absolute radial Jacobian. The original probe stays incoming.

Keep the thermal addition centered on signed cuts and occupation factors. Extend existing
owners when adding capabilities; avoid introducing a second thermal residue or event pipeline.
Shared signed algebra and derivative jets do not by themselves validate raised thermal cuts.

== Supported scope and deferred work

The enabled thermal `xs` path has one future-directed timelike probe, finite temperature,
massless scalar cut particles, no chemical-potential assignments, and simple cuts without
virtual propagators. The generation loop-momentum basis must route external momentum through
at most one cut edge that also carries loop momentum.

Signed radial isolation currently covers
$A r + B + s sqrt(a r^2 + 2 b r + c)$ after combining identical affine energies.
It enumerates simple positive roots within this chart and distinguishes absent roots from
unsupported equations and unresolved tangencies. General mixed-sign root isolation needs
further algorithms; keep that work separate from positive-E geometric assumptions.

Thermal threshold subtraction remains out of scope. Disable it both at generation and runtime.
Cold modes, massive or spinning cut particles, raised thermal cuts,
UV subtraction, and surface-targeted sampling remain outside the enabled thermal-xs scope.
These are current validation limits, not a reason to specialize the common residue machinery.

Focus on infrared-finite cases. Thermal infrared subtraction and screening/resummation are
not implemented, and no new infrared-finiteness checks were added. Spacelike and near-light-cone
QCD vector benchmarks and their possible additional subtraction requirements remain later work.

Saved-state format is 16. Regenerate older states; the persisted cut layouts changed.

== Reproduction and verification

From the repository root:

```sh
just build-cli-no-pyo3
./gammaloop -n examples/cli/thermal_rates/scalar_sunset/run.toml run generate rest inspect
python3 examples/cli/thermal_rates/scalar_sunset/reference.py --orders 64 128 256 512 --output /tmp/thermal-rates-sunset-reference.json
```

Use `boosted` or `fast` instead of `rest` for probe momenta 1 or 3. The card also supplies
an `integrate` block for `equilibrium`, `vacuum`, and `thermal_excess`; use separate scratch
workspaces for distinct kinematics.
#link("scalar_sunset/validation.json")[The validation record] retains exact commands,
settings and executable hashes, numerical values, and log paths. Its raw `/tmp` artifacts
are local to the validation host and are not included in this commit.

=== Vacuum-subtraction continuation

The clean `3e73e6b07` checkout passed 195 selected tests before this extension
(run `caaaf13b-1dd3-43bc-88a2-472b21e0b4a2`). After enabling the existing subtraction
option, all 12 targeted tests passed (run `bced47f4-4695-41e5-a347-a5c772800dc6`).
The scalar sunset regression checks total and per-event subtraction at momenta 0, 1,
and 3, saved-state reload, conservation, closure, selectors, and an independently computed
positive thermal excess of approximately `1.7071821495e-30` at inverse temperature 200.
This contribution is too small to recover by subtracting full binary64 evaluations.

Relevant Clippy checks passed with warnings treated as errors, changed Rust sources
passed formatting, and the four edited Typst documents rendered. The fresh CLI rest
smoke gives `thermal_excess = 5.248854262631212e-05` with four events, exactly equal
to equilibrium minus vacuum, including each physical channel. Event closure is exact;
the maximum four-momentum residual is approximately `2.22e-16`.
These checks ran on `itphlies`; commands, hashes, run IDs, and logs are recorded under
`vacuum_subtraction_continuation` in the validation record. Saved-state format remains 16.
No new integration runs were performed.

=== Earlier milestone verification

- The shared refactor passed 195 focused tests across five binaries with `just test`,
  `ci-optim`, and the `ci_gammaloop` profile. Run ID:
  `e8fd0805-f408-427a-bc7d-6c883d997ae4`. The suite covers residue signs and aliases,
  E/H rays, existing raised vacuum behavior, sampling, saved-state reload, boosted sunset
  kinematics, event conservation and closure, and selectors with event retention on and off.
- After setting the final manifest number, all five manifest compatibility tests passed.
  Run ID: `2a70b1e2-a8b6-40ba-8622-d0e415a251d7`;
  log: `/tmp/thermal-rates-state16-tests.log`.
- Relevant Clippy checks passed without warnings, changed Rust files passed formatting checks,
  and the architecture and benchmark Typst documents compiled. Final manifest-edit Clippy log:
  `/tmp/thermal-rates-state16-clippy.log`.
- The rebuilt CLI smoke reproduced the prior rest-frame point values exactly:
  equilibrium `5.2497393862529534e-05` with four events and vacuum
  `8.851236217417727e-09` with one event. Event weights sum to the integrand;
  maximum four-momentum residual is approximately `2.22e-16`.
- Historical pointwise validation checked all twelve signed events at momenta 0, 1, and 3
  against independent phase space within `3.4e-13` relative. The earlier integrations at
  momenta 0 and 1 agree with the reference. A discrepant momentum-3 pilot remains recorded;
  its independent million-sample repeat agrees, with approximately 2.85% uncertainty.
  The historical 0.5% integration target was not reached; reaching it is not a prerequisite
  for extending validation to other infrared-finite cases. Those integrations were not rerun
  after the structural refactor.

These are focused checks, not the full repository suite. Top-level `enable = false`
in `nix-ci.nix` is intentional; no NixCI artifacts were uploaded and no final-review label
was requested. No automatic CI completion follow-up is registered.

== Suggested continuation

+ Continue from the shared residue and radial architecture, using the scalar sunset card
  as a regression. For each additional infrared-finite case, check which existing machinery
  applies and identify concrete missing functionality before adding code. Preserve
  phase-space checks, vacuum controls, event conservation, and selector behavior.
+ Broaden signed root isolation and validation case by case. Retain complete root accounting
  and explicit uncertainty for unsupported charts or tangencies.
+ Validate each additional physical capability before relaxing its scope guard, including
  the Jacobian and derivative conventions for raised thermal cuts.
+ Keep thermal threshold subtraction and thermal infrared treatment separate from this
  milestone. Before requesting final review, follow the enable/check/upload/push/label
  workflow in `CONTRIBUTING.typ`.
