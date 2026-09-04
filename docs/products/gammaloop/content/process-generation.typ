#import "../../shared.typ": callout, boundary, developer-link, source-link

#let process-generation = [
= Process generation and state workflows

GammaLoop describes a calculation with a TOML run card and stores its reusable state in a
directory. Start the command-line interface through the repository wrapper (`./gammaloop`) or
the Cargo-built executable, then generate processes, integrands, and integration jobs within
that state.

== Choose a generation mode

The interactive command tree separates cross sections, amplitudes, and operations on processes
already stored in the active state:

// docs-example: syntax
```text
generate xs INITIAL > FINAL ...
generate amp INITIAL > FINAL ...
generate existing ...
generate help xs
```

`xs` constructs forward-scattering graphs and filters compatible Cutkosky cuts. `amp` constructs
ordinary amplitudes. Their final-state lists therefore have different meanings even when the
printed particle syntax is identical. Use `generate help <mode>` inside the stateful CLI to see
the flags supported by your installed version.

== Process specification grammar

Process specifications accept `to` or `>` between space-separated initial and final states.
Braces express alternatives in a cross-section final state, while `{}` denotes an empty state:

// docs-example: syntax
```text
e+ e- > mu+ mu-
e+ e- > { Z Z, a a, H H }
{} > {}
```

The empty-initial-state form is generation grammar, not evidence that the default
cross-section flux evaluator supports a zero-incoming-particle calculation. Its implemented
flux branches require exactly one or two incoming particles.

A slash vetoes particles and a vertical bar selects an allow-list. The words `veto` and `only`
are equivalent forms. Outside the bracketed perturbative block, coupling constraints use the
active model's coupling-order name:

// docs-example: syntax
```text
e+ e- > d d~ g / u c QED==2 QCD>=2 QCD<=4
e+ e- > mu+ mu- | g ghG QED==2 QCD>=2
```

#callout("Powered and unpowered constraints are not interchangeable", [
  An unpowered spelling such as `QCD==2` builds an amplitude-side coupling-order filter. A
  powered spelling such as `QCD^2==2` builds a cross-section-side coupling-order filter, but
  the parsed exponent is currently discarded. Substituting one spelling for the other changes
  which filter container is used and still does not define a squared amplitude or a complete
  cross-section order. Production work therefore needs a scientifically reviewed selector set;
  do not infer one from the spelling alone.
])

The bracketed perturbative block is intentionally distinct from ordinary flags. For `amp`,
`{n}` fixes the amplitude graph's loop count; for `xs`, it fixes the sum of loop counts on the
two cut sides. `{{n}}` fixes the forward graph's loop count. On a cross-section process,
bracketed `QCD=n` or the shorthand `QCD` resolves the model's unresolved-particle set and permits
at most the summed requested count of additional cut particles from those sets; it is not an
exact graph coupling power. No corresponding amplitude-side selection effect has yet been
established for that bracketed named-order spelling:

// docs-example: syntax
```text
e+ e- > Z [ {1} {{2}} QCD=2 QED=1 ]
e+ e- > Z [ QCD ]
```

These selectors are generation filters, not universal definitions of LO, NLO, or NNLO. Record
the amplitude-side loops, forward loops, unresolved-cut allowance, and coupling filters
explicitly for the chosen generation mode; use the #link("physics/")[physics-scope guide] to
distinguish that calculation definition from the formal perturbative reach of Local Unitarity.

#callout("Resolve names through the active model", [
  Particle names, antiparticles, coupling-order names, and allowed interactions come from the
  model loaded in the current state. A syntactically valid specification can still be rejected
  when that model does not define a requested name.
])

== Filters, grouping, and multiplicity

Generation flags control topology filters, loop ranges, cut blobs and spectators,
symmetrization, and numerator-aware grouping. These choices are physics-significant. In
particular, grouping can combine graph multiplicities and numerator ratios; external-fermion
symmetrization requires the caller to retain the associated ordering signs. Run
`generate help <mode>` before choosing filters to see the available aliases, defaults, values,
and argument counts.

The generated graph's overall factor keeps separate contributions for automorphisms, internal
fermion loops, external-fermion ordering, numerator-independent symmetry grouping, and
numerator-dependent grouping. It does not establish a universal incoming-color average.
Persist the state and inspect that factor together with any explicit projector when auditing
diagram normalization instead of replacing it with an assumed factorial; the
#link("guides/conventions/")[normalization checklist] assigns each factor to its owner.

== A maintainable generation card

The following reduced card uses the repository's scalar model and the same one-loop bubble
command exercised by the maintained scalar-topology example. Save it as `manual-bubble.toml`
and run `./gammaloop --clean-state manual-bubble.toml run generate -c "quit -o"` from the
repository root.

// docs-example: syntax gammaloop-process-generation
```toml
commands = []

[cli_settings.state]
folder = "./gammaloop_state/manual-bubble"
name = "manual_bubble"

[[command_blocks]]
name = "generate"
commands = [
  "import model scalars-default.json",
  "generate amp scalar_1 > scalar_1 [{1}] --allowed-vertex-interactions V_3_SCALAR_122 -p bubble -i scalar_bubble",
  "generate",
  "save state -o",
]
```

The acceptance invariant is a persisted state in which `display processes` identifies process
`bubble` and integrand `scalar_bubble`; replaying the card with a clean destination must produce
the same named objects. The exact command and setting semantics are linked from the
#link("reference/cli/commands/gammaloop/generate/amp/#command-gammaloop-generate-amp-23b125bb12cb4529")[amplitude-generation reference],
#link("reference/cli/commands/gammaloop/save/state/#command-gammaloop-save-state-6caa72db99ca51d0")[state-persistence reference], and
#link("reference/cli/settings/cli/state/#setting-cli-state-folder-ea5760beb855f86c")[state-folder setting]. The full maintained source is
the #source-link("examples/cli/scalar_topologies/bubble.toml", label: "scalar bubble run card").

#callout("Interpret the first failing stage", [
  An unknown particle or `V_3_SCALAR_122` error means the scalar model was not imported or no
  longer exposes that interaction. A process with no retained diagram points to the process
  specification or generation filters. A generated process that disappears after restart means
  `save state -o` did not complete or the resumed `-s` path differs from `cli_settings.state.folder`.
])

== From generation to evaluation

A reproducible run keeps its stages explicit:

- import or select the model before resolving a process specification;
- generate the process and inspect the retained diagrams;
- generate the requested integrand representation;
- set kinematics, sampling, stability, and subtraction settings;
- evaluate samples or run an integration command block;
- exit with persistence enabled when the state should be resumed.

== Select integration slots and targets

An integration slot is one `(process, integrand)` pair. Repeat `-p/--process` and
`-i/--integrand-name` in matching order to integrate several slots together. By default those
slots train and sample a shared grid, so their sample shapes must be compatible; add
`--uncorrelated` when each slot needs an independent grid.

A single `--target re im` (or `--target re,im`) applies the same known comparison value to every
selected slot. Multi-slot runs can instead use one keyed value per slot:

// docs-example: syntax
```text
integrate \
  -p bubble_no_integrated_UV -i scalar_bubble_below_thres \
  -p bubble -i scalar_bubble_below_thres \
  --target bubble_no_integrated_UV@scalar_bubble_below_thres=1.471664021721597e-2,0.0 \
           bubble@scalar_bubble_below_thres=2.7029875684552542e-3,0.0 \
  --workspace-path ./integration/bubble-below-threshold \
  --n-cores 1
```

Targets annotate result deltas; they do not replace `integrator.target_relative_accuracy`,
`integrator.target_absolute_accuracy`, or the sample-budget settings that control convergence.
The maintained
#source-link("examples/cli/scalar_topologies/bubble.toml", label: "scalar bubble run card")
contains the complete settings and both below- and above-threshold target sets. The generated
#link("reference/cli/commands/gammaloop/integrate/#command-gammaloop-integrate-eac2ec37911d2fab")[integration reference] is authoritative for
selectors, renderers, batching controls, and allowed values.

== Resume or replace an integration workspace

An integration workspace is a completed-iteration checkpoint, not a second GammaLoop state.
Without `--workspace-path`, writable sessions use `STATE_FOLDER/integration_workspace`;
read-only sessions use `./integration_workspace_STATE_NAME` (or `./integration_workspace`) so
they do not write into the active state.

- Without `--restart`, GammaLoop resumes when the workspace already contains a compatible
  checkpoint. Slots, shared-versus-uncorrelated sampling, effective model parameters, and
  generated-integrand fingerprints must match. Saved targets and non-model runtime settings win
  over changed values and produce a warning.
- With `--restart`, GammaLoop removes that specific workspace and starts the integration from
  scratch. Use it when changing slots, sampling correlation, generated integrands, or settings
  that should take effect. It does not mean “resume”.
- `--show-summary-only` reads the last completed result without evaluating more samples and
  cannot be combined with `--restart`.

The workspace root contains `manifest.json`, the latest user-facing `integration_result.json`,
per-slot settings and results under `integrands/`, and the internal resume checkpoint under
`state/integration_state.bin`. Observable resume snapshots live at
`state/observables/<process@integrand>/latest.json`. Configured user-facing output is written per
slot as `integrands/<process@integrand>/observables_final.json` or
`integrands/<process@integrand>/observables_final.hwu`. `--write-results-for-each-iteration`
additionally retains numbered result and observable history at the corresponding root or slot
boundary.

#callout("Archive the right boundary", [
  `integration_result.json` is a presentation snapshot, not the authoritative resume state.
  Keep the complete integration workspace when a run must be resumed, and keep the GammaLoop
  state/run card that produced the referenced integrands. A copied result JSON alone is suitable
  for analysis, not continuation.
])

#boundary("Keep the complete state together", [
  `run.toml` records boot settings, reusable command blocks, and top-level commands. Generated
  process and integrand artifacts live beside it. Resume with `-s STATE_FOLDER`; replay the run
  card to reconstruct the state from scratch. Copy or archive the entire state directory rather
  than moving a generated subdirectory on its own, because settings and fingerprints are
  checked together.
])

Contributor-facing ownership and data flow are documented in the
#developer-link("gammaloop-architecture", "architecture-current.typ", "GammaLoop architecture")
and the
#developer-link("uv-renormalization", "uv-renormalization.typ", "UV renormalization architecture").
]
