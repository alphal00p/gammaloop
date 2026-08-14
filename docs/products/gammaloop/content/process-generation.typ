#import "../../shared.typ": callout, boundary

#let process-generation = [
= Process generation and state workflows

GammaLoop describes a calculation with a TOML run card and stores its reusable state in a
directory. Start the command-line interface through the repository wrapper (`./gammaloop`) or
the Cargo-built executable, then generate processes, integrands, and integration jobs within
that state.

== Choose a generation mode

The interactive command tree separates cross sections, amplitudes, and operations on processes
already stored in the active state:

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

```text
e+ e- > mu+ mu-
e+ e- > { Z Z, a a, H H }
{} > {}
```

A slash vetoes particles and a vertical bar selects an allow-list. The words `veto` and `only`
are equivalent forms. Amplitude coupling constraints use the coupling name directly; powered
constraints apply to the complete cross section:

```text
e+ e- > d d~ g / u c QED==2 QCD>=2 QCD<=4
e+ e- > mu+ mu- | g ghG QED^2==2 QCD^2>=2
```

The bracketed perturbative block is intentionally distinct from ordinary flags. `{n}` fixes the
amplitude loop count (or the sum across the two cut sides), `{{n}}` fixes the forward-graph loop
count, and `QCD=n` or the shorthand `QCD` selects a relative perturbative order:

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
