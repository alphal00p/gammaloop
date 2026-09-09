# Process and graph generation

GammaLoop has two generation modes. Cross-section generation (`xs`) builds
forward-scattering graphs and compatible cuts. Amplitude generation (`amp`,
also accepted as `amplitude`) builds ordinary outgoing amplitudes. `existing`
reuses a process already present in the state.

The command grammar is:

```text
generate xs  INITIAL to FINAL [process-options] [generation-options]
generate amp INITIAL >  FINAL [process-options] [generation-options]
generate existing --process PROCESS [--integrand-name NAME]
```

`>` and `to` are equivalent separators. Particle names come from the active
model. In XS mode, `FINAL` selects compatible Cutkosky cuts; unresolved
massless particles may be added according to the perturbative order. A brace
list selects several allowed final-state channels:

```text
generate xs e+ e- > {z z, a a}
```

Empty braces are useful for vacuum graphs. In XS mode, an empty final-state set
also disables final-state cut filtering.

## Process filters

```text
generate xs e+ e- > d d~ g / u d
generate xs e+ e- > d d~ g | u d c s
generate xs e+ e- > d d~ g QED==2 QCD>=2 QCD<=4
generate xs e+ e- > d d~ g QED^2==4 QCD^2<=2
generate xs e+ e- > d d~ g [{1} {{2}} QCD=2 QED=1]
```

- `/` (or `veto`) removes the listed particle species before graph generation.
- `|` (or `only`) keeps only the listed particle species.
- Bare coupling names and comparisons such as `QED==2` constrain each generated
  amplitude. In XS mode, both sides of a cut must satisfy them.
- Squared coupling names such as `QED^2==4` constrain the complete
  forward-scattering graph and are available only in XS mode.
- The `[...]` block selects perturbative terms. `{n}` fixes the amplitude loop
  count, `{{n}}` fixes the forward-graph loop count, and `QCD=n` selects the
  corresponding perturbative order. A bare `QCD` means `QCD=1`.

## Generation options

Use `./gammaloop help generate` for the complete command-generated list. The
options below change graph construction or grouping.

### State and topology

- `-a` adds results to existing processes.
- `--clear-existing-processes` (`-c`, `--clear`) removes existing processes
  before generating.
- `--process-name` (`-p`) and `--integrand-name` (`-i`) name generated
  resources. `--only-diagrams` stops after graph generation. `--keep-sources`
  retains generated C++ source files after compilation.
- `--filter-self-loop` removes one-edge self-loops.
- `--filter-selfenergies`, `--filter-snails`, and `--filter-tadpoles` remove
  the corresponding special topologies. The `--veto-*` and
  `--veto-only-scaleless-*` options select whether massive, massless, or only
  scaleless instances are vetoed.
- `--max-n-bridges` limits graph bridges. `-1` disables this limit.
- `--number-of-factorized-loop-subtopologies MIN MAX` filters by the number of
  factorized loop subtopologies; negative bounds disable it.
- `--number-of-anticommutating-loops MIN MAX` filters closed fermion or ghost
  loops; negative bounds disable the filter.
- `--veto-vertex-interactions` and `--allowed-vertex-interactions` restrict
  model interaction identities.
- XS-only cut-shape filters are `--n-cut-blobs`, `--n-cut-spectators`, and
  `--filter-cross-section-tadpoles`. `--filter-zero-flow-edges` removes edges
  with zero flow, and `--max-multiplicity-for-fast-cut-filter` limits the fast
  cut-filter path.

Defaults for special-topology filters depend on the generation mode and whether
the process is a vacuum process. Set them explicitly when reproducibility is
important.

### Symmetrization and grouping

External symmetrization is disabled unless requested:

- `--symmetrize-initial-states` and `--symmetrize-final-states` group graphs
  differing only by permutations of those external labels.
- `--symmetrize-left-right-states` identifies left/right-related graphs and
  requires the selected theory and process to support that CP identification.
- `--allow-symmetrization-of-external-fermions-in-amplitudes` (`--symferm`)
  enables amplitude external-fermion permutations. Use it only when the
  associated fermion-flow signs and helicity routing are understood.

Numerator-aware grouping is selected with `--numerator-grouping`:

```text
--numerator-grouping no_grouping
--numerator-grouping only_detect_zeroes
--numerator-grouping group_identical_graphs_up_to_sign
--numerator-grouping group_identical_graphs_up_to_scalar_rescaling
```

The comparison strategy can be tuned with
`--consider-internal-masses-only-in-numerator-isomorphisms`,
`--compare-canonized-numerator`,
`--number-of-samples-for-numerator-comparisons`,
`--numerical-samples-seed`, and
`--fully-numerical-substitution-when-comparing-numerators`, and
`--symmetric-left-right-polarizations`. Numerical grouping
uses exact integer arithmetic on sampled components; sampled momenta are not
physical kinematics, and fully numerical substitution fixes `N_c=3` while
replacing functions with random primes.

### Selecting graphs and loop bases

```text
--graph-prefix partonic_channel
--select-graphs GL_12 GL_13
--veto-graphs GL_11 GL_15
--loop-momentum-bases GL_12=7,10 GL_77=4,2
```

The loop-basis map names a graph and lists its edge indices. GammaLoop does not
use this option to prove that a supplied basis is topologically valid; invalid
choices fail when the basis is used.

Global numerator/projector overrides are available through
`--global-prefactor-num` and `--global-prefactor-projector`. Quote Symbolica
expressions when the shell would otherwise split or expand their punctuation.

Generation writes processes into the active state. After changing generation
options, regenerate dependent integrands and compiled evaluators so that the
saved graph set and its grouping metadata agree with the new process.
