#import "../../shared.typ": callout, boundary, source-link

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

A slash vetoes particles and a vertical bar selects an allow-list. The words `veto` and `only`
are equivalent forms. Amplitude coupling constraints use the coupling name directly; powered
constraints apply to the complete cross section:

// docs-example: syntax
```text
e+ e- > d d~ g / u c QED==2 QCD>=2 QCD<=4
e+ e- > mu+ mu- | g ghG QED^2==2 QCD^2>=2
```

The bracketed perturbative block is intentionally distinct from ordinary flags. `{n}` fixes the
amplitude loop count (or the sum across the two cut sides), `{{n}}` fixes the forward-graph loop
count, and `QCD=n` or the shorthand `QCD` selects a relative perturbative order:

// docs-example: syntax
```text
e+ e- > Z [ {1} {{2}} QCD=2 QED=1 ]
e+ e- > Z [ QCD ]
```

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
numerator-dependent grouping. Persist the state and inspect that factor when auditing diagram
normalization instead of replacing it with an assumed factorial.

== A maintainable generation card

The following reduced card uses the repository's scalar model and the same one-loop bubble
command exercised by the maintained scalar-topology example. Save it as `manual-bubble.toml`
and run `./gammaloop --clean-state manual-bubble.toml run generate -c "quit -o"` from the
repository root.

// docs-example: syntax
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
#link("reference/cli/?q=generate%20amp")[amplitude-generation reference],
#link("reference/cli/?q=save%20state")[state-persistence reference], and
#link("reference/cli/?q=cli.state.folder")[state-folder setting]. The full maintained source is
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

#boundary("Keep the complete state together", [
  `run.toml` records boot settings, reusable command blocks, and top-level commands. Generated
  process and integrand artifacts live beside it. Resume with `-s STATE_FOLDER`; replay the run
  card to reconstruct the state from scratch. Copy or archive the entire state directory rather
  than moving a generated subdirectory on its own, because settings and fingerprints are
  checked together.
])
]
