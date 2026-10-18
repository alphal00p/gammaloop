= Thermal rates development handoff

== Current milestone

Branch: `thermal_rates`, based on `ff533fcd9` from main. The implementation was rebuilt
against this base; the old `thermal_rate_investigation` branch is background material,
not a patch source. The first supported example is the scalar sunset with simple thermal
cuts, using the existing signed CFF representations and thermal weights.
The 2026-10-17–18 continuation broadens root geometry and raised-cut normalization, reuses
positive-E sampling maps, and connects cold Fermi sectors to the common rate evaluation.
The initial continuation’s #link("scalar_bubble/run.toml")[massive scalar bubble card] and
#link("scalar_bubble/benchmark.typ")[literature benchmark] provide IR-finite decay,
absorption, spacelike, and raised-propagator controls.
The subsequent inclusive-parity continuation removes the process-specific thermal
whitelist and reuses the full vacuum graph, numerator, virtual-line, and local/integrated
UV machinery with signed cuts and thermal weights. The restored
#link("massive_double_triangle/benchmark.typ")[double triangle] and
#link("massive_self_energy/benchmark.typ")[UV-renormalized insertion] distinguish
historical evidence from current validation and its numerical precision limits.

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
  Signed derivative jets enter the Laurent residue unchanged; one root-orientation factor
  converts it to the absolute measure at every pole order. Absorbed bath particles become
  incoming particles before selectors run. The original probe stays incoming.
- `crates/gammalooprs/src/integrands/process/fermi_surface/`: amplitudes and rates share
  orientation-aware sector extraction, shell parameters, derivative profiles, and coefficient
  evaluation. Independent shells and the regular single-shell/two-energy angular chart
  prepare contributions for the common physical branch and event loop.
- `crates/linnet/src/half_edge.rs`: cut enumeration separates visited source growth from
  accepted partitions. Connectivity bounds apply to accepted cuts; an intermediate disconnected
  complement may become valid after absorbing a source-attached spectator cycle.

Keep the thermal addition centered on signed cuts and occupation factors. Extend existing
owners when adding capabilities; avoid introducing a second thermal residue or event pipeline.
Shared signed algebra and derivative jets do not by themselves validate every raised process.

== Supported scope and deferred work

Thermal `xs` uses the ordinary inclusive cross-section scope: the shared vacuum
machinery handles graph topology, particle spin, masses, chemical-potential assignments,
virtual propagators, and raised cuts. There is no scalar-only, one-probe, or virtual-spectator
whitelist. The usual model-domain checks and supported vacuum kinematics still apply;
removing a whitelist does not certify every process or cancel singularities in an incomplete
cut sum. Infrared-finite inclusive observables are the present target.

Signed radial isolation now uses directed value and derivative enclosures on the entire
compactified positive ray. Certified empty cells and simple crossing brackets account for
all roots without a fixed cap on distinct energies or an elimination-polynomial degree cap.
Positive-E rays starting inside their convex surface retain the existing inexpensive native
solve. Exterior rays can have two roots; certified closed channels contribute zero.
Unresolved tangencies, vanishing-energy cusps, and exhausted subdivision or precision budgets
remain explicit errors, never evidence that a channel is absent. Ordered signed energy
occurrences and their original affine routing are retained for every residue jet.

Existing positive-E sampling maps are available for thermal rates. Registration, host matching,
and cached-host adoption retain positive-E ownership even when H groups share physical cut IDs.
H cuts use their own radial equations on the resulting mapped sample. Strict-interior map
requirements and full-support mixture checks are unchanged; signed-H targeting is not provided.

Cold rate coefficients use the shared Fermi-sector evaluator. If a cut is independent of
the localized shell coordinates, existing shell derivative products apply with constant LU
parameters. A dependent two-energy cut can instead use one regular Fermi delta and a simple
LU coefficient: the shell fixes the radius, the cut fixes the polar angle, and the azimuth is
retained. Its angular Jacobian supplies the cut measure without a second LU profile. The
prepared contribution and physical event enter the common branch loop. Higher delta derivatives
on a dependent cut, multiple dependent shells, zero angular axes, tangent intersections, and
occupation-step intersections remain unsupported. These restrictions are distinct from the
existing higher-derivative support for independent shells.

Double LU poles also differentiate order-zero cold occupations. Generation now extracts
that distributional derivative symbolically after selecting each retained orientation and
resolving scalar definitions, preserving factorized graph numerators. Contour offsets are
canonicalized before differentiation: they describe the same occupation plus a constant,
so their derivatives must combine before checking for dependent steps. The occupation's first
energy derivative carries its original affine radial energy derivative and the first-derivative
coefficient from the shared Laurent-residue builder. Converted to a simple-cut coefficient,
it enters the same Fermi-sector localizer and physical event loop. The smooth raised residue
is retained. Active cold occupations at LU order above two, or existing positive-order Fermi
distributions in a raised LU coefficient, remain explicit higher-intersection errors.

Thermal threshold subtraction remains out of scope. Disable it both at generation and runtime.
Thermal rates reuse ordinary local 3D UV subtraction and integrated counterterms, including
their finite renormalization-scale dependence. Vacuum subtraction of the whole medium rate
does not replace subgraph renormalization. Expanded-4D local medium counterterms remain
unsupported. Keep algorithmic limitations separate from validation coverage.

Focus on infrared-finite cases. Thermal infrared subtraction and screening/resummation are
not implemented, and no new infrared-finiteness checks were added. Spacelike and near-light-cone
QCD vector benchmarks and their possible additional subtraction requirements remain later work.

Saved-state format is 17. Regenerate older states; cross-section cut groups now persist
Fermi-surface sectors and their cut-localization classification. Composite edge masses
also retain the shared native-precision evaluator, matching shell and coefficient arithmetic.

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

=== Inclusive vacuum-parity continuation: implementation and validation

The comparison with `origin/thermal_rate_investigation` at
`a0b2d4191cad17d8dc2f20891bb0f245f0df18b7` identified previously demonstrated
finite-temperature scalar virtual corrections and local-plus-integrated UV renormalization
that the rebuilt branch’s temporary guards excluded. The current continuation removes those
policy restrictions through the existing shared owners; it does not import a parallel thermal
evaluator or the old dependency-specific hyperbolic workaround.

Restored entry points are #link("massive_double_triangle/run.toml")[the double-triangle card]
for below/above-threshold real–virtual sums, #link("massive_self_energy/run.toml")[the exact
raised insertion card] with both UV pieces and an auxiliary-scale control, and
#link("qcd_vector/run.toml")[the inclusive QCD vector-current card]. Thermal threshold
subtraction remains disabled. The above-threshold scalar comparison tests cancellation
inside the inclusive cut sum rather than adding a new thermal subtraction prescription.

#link("inclusive_validation.json")[The current acceptance manifest] records source and binary
fingerprints, the test selection, checks, benchmark records, and numerical limitations.

The restored independent NumPy quadratures reproduce `3.057651940618648e-4` below
threshold, approximately `5.2820877369156e-6` above threshold, and
`1.4011350428792406e-4` for the MSbar insertion. Their scripts, source formulas, and
historical comparison records are retained in the corresponding example directories.
The restored Typst benchmark/derivation documents render; the input JSON/TOML and Python
syntax checks and Ruff pass. All 116 selected native regressions passed in run
`07460a80-6bc8-4df6-b20a-d00fc0f77ac0` (64.499 seconds), including restored scalar
point/channel/event closure, UV finite-addback and auxiliary-scale controls, state reload,
and the QCD fixtures. Current scalar integrations total 1.1 million process evaluations,
with zero final NaNs or unresolved stability failures. The below-threshold inclusive
estimate is `3.01438e-4 +/- 6.50100e-6`; the independent-seed UV-renormalized estimate is
`1.35990e-4 +/- 3.40253e-6`. Both agree with their references within reported Monte Carlo
uncertainties, as do the auxiliary-scale and channel controls. Above threshold, the
independent-seed estimate `4.77231e-6 +/- 1.37286e-6` agrees with `5.28209e-6`; its
28.8% error and unresolved virtual relative precision remain documented. Four direct
pole-approach controls independently reproduce the historical joint cancellation within
`1.54e-16` relative. The QCD card records its separate current integrated comparison.

Historical million-sample runs, interrupted pilots, and accepted current estimates remain
separate in the validation records. Rejected low-precision attempts are rescued by higher
precision and must not be counted as discarded samples. For the scalar insertion,
different default sampling LMBs explain why correlated cube samples do not enforce exact
closure of the final channel means; selecting the same basis restores the checked
pointwise closure to `3.39e-21`.

=== 2026-10-17–18 implementation and validation at `4a6519c64`

The continuation from `9bc647ef2fb7ded5b6590f2e112a2b14a27d9c22` implements the root,
raised-cut, positive-E sampling, and cold-sector connections in commit `4a6519c64`. Its root
isolation used a six-energy algebraic elimination chart, superseded by the interval solver
in the inclusive continuation. All 111 selected
tests passed in run `ff7869bf-c0be-49ae-8e97-6881e41a0363` on `itphlies`, using the
`ci-optim` Cargo profile and `ci_gammaloop` Nextest profile. Runtime was 171.028 seconds;
log: `/tmp/thermal-rates-verified-tests-3.log`. This is a targeted regression run, not the full suite.

#link("scalar_bubble/validation.json")[The continuation validation record] preserves the
source hashes, exact test selection, run identities, native comparison values, and log paths.
The retained CLI artifacts are under `/tmp/thermal-rates-verified-artifacts-3`; core test logs
are under `/tmp/thermal-rates-verified-logs-3`. These paths are local to the validation host.

The primary physics controls are documented in
#link("scalar_bubble/benchmark.typ")[the scalar bubble and cold Yukawa benchmark]:
#link("https://arxiv.org/abs/1503.03534")[Ho and Scherrer (2015)] for boosted timelike decay,
#link("https://arxiv.org/abs/hep-ph/0302098")[Nishikawa, Morimatsu, and Hidaka (2003)] for
massive bubble cut domains, #link("https://arxiv.org/abs/2007.01224")[Lundberg and Pasechnik
(2021)] for finite-density fermionic phase space, and
#link("https://arxiv.org/abs/2208.14479")[Gorda, Österman, and Säppi (2022)] for cold
pole-derivative boundary terms. The independent Python scalar reference uses no native
roots, residues, or sampling profiles.

The native scalar checks cover rest and boosted decay, closed channels, absorption, and
spacelike support. Maximum relative rate error is `2.65e-5`, below the `1e-4` requirement;
the closed channel is exactly zero. Both raised scalar mass-squared derivatives, including
negative-slope absorption, agree within `8.2e-15` relative, against a `1e-8` requirement.
The physical Yukawa checks retain native Dirac and color numerators. At both chemical-potential
signs, ordinary rates agree within `1.0e-7`, the individual mass insertions within `6.6e-6`,
and their sums within `2.0e-6` relative. All use an unchanged `1e-4` acceptance bound.
The Gauss angular and radial refinements were checked independently before selecting them.
These are deterministic quadrature controls, not new Monte Carlo convergence measurements.

The run also covers the existing sunset, positive-E thermal sampling with signed channels,
amplitude Fermi sectors, virtual-spectator localization, event conservation and weight closure,
saved-state reload, root isolation, cut enumeration, and physical-basis validation.
Composite-mass checks compare exact `sqrt(2)` and `1/3` expressions in Quad, Arb, and fixed256
lanes after graph reload. A rotated near-axis angular regression retains the same strict
cut and shell checks in Double, Quad, and Arb.

The physical mass-insertion comparison exposed a missing boundary from LU differentiation
of cold occupations. Before the correction, its discrepancy matched the independently
predicted missing term to approximately `4e-11` absolute. The corrected native comparisons
now pass. Exact contour-offset cancellation and the stable transverse construction were
necessary to retain valid crossed coefficients and translated-shell sampling points.

Relevant Clippy checks pass with warnings treated as errors; changed Rust files pass formatting,
the Python reference passes Ruff, and the authored Typst notes render. The independent
128-point scalar occupation quadrature agrees with the closed formulas within `1.1e-17`
absolute. Full-suite checks, CI artifact upload, and final-review readiness are not claimed.

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
+ Extend the literature comparisons to further infrared-finite processes, retaining event
  conservation, selector behavior, saved-state reload, and independent mass-derivative controls.
  Assess each remaining scope guard against existing shared functionality before adding code.
+ Broaden the remaining signed-root and joint-surface charts case by case. Retain complete
  root accounting and explicit uncertainty for unsupported charts or tangencies; validate each
  capability before relaxing its guard.
+ Keep thermal threshold subtraction and thermal infrared treatment separate from this
  milestone. Before requesting final review, follow the enable/check/upload/push/label
  workflow in `CONTRIBUTING.typ`.
