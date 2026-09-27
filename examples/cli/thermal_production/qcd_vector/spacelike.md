# Spacelike vector-current diagnostic

At `T=1`, `K=(1,0,0,4)` has `K^2=-15`. It probes the same massless QCD
vector-current correlator as the timelike benchmark, with the same couplings,
UV settings, and `medium.vacuum_subtraction=false`.

The vacuum spectral function vanishes here: a spacelike momentum cannot equal
a sum of positive-energy on-shell final-state momenta. In a bath, mixed-sign
cuts describe scattering, for example `K + p = p'`, with
`k0 = E(p') - E(p)`. They have spacelike support and give Landau damping.
This is a current response to spacelike momentum transfer, not the production
of an on-shell particle with spacelike four-momentum. See
[arXiv:1910.07552v2, Section 4 and Appendix A](https://arxiv.org/pdf/1910.07552v2).

## Numerical comparison

Values below are `rho_V/T^2`, including `Nc=3`, with the electromagnetic
charge removed and `g=1`. NLO means the correction alone. The native-to-current
conversion is now `rho_V = 9*Gamma_rho`, since `k0=1`; using the timelike
factor 36 would be incorrect.

| Contribution | GammaLoop | Reference |
| --- | ---: | ---: |
| Thermal LO, 100,000 samples | 0.21817041 ± 0.00224614 | 0.21937756 |
| Thermal NLO, 100,000 samples | -0.02465424 ± 0.00061311 | -0.02206600 |
| Thermal NLO, independent seed, **150,000-sample checkpoint** | -0.02430015 ± 0.00056276 | -0.02206600 |
| Vacuum LO and NLO, 1,000 samples each | exactly zero | zero |

LO agrees within 0.54 Monte Carlo errors. The completed NLO run differs by
4.22 quoted errors (11.7% of the reference correction's magnitude). The
independent seed trends in the same direction, but **did not finish** its
requested 300,000 samples. Its table entry is the last completed iteration,
not a completed integration result.

The author's reference driver gives NLO `-0.022065903635688392` at tolerance
`1e-4` and `-0.022066002519195806` at `1e-5`. Their difference, `9.89e-8`, is
well below the native errors. The completed 100,000-sample run and the saved
150,000-sample checkpoint report no nonfinite or unstable samples.

## Concrete failure and limits of the diagnosis

The independent run exited with status 1 during its fourth 50,000-sample
iteration:

```text
Cannot isolate simple radial roots in
[0.9999998078938658, 0.9999998078938659]; increase precision
```

The error originates in `all_roots_in_interval` at
`crates/gammalooprs/src/utils/newton_solver.rs:773`, called by
`positive_radius_roots` from the LU evaluator at
`crates/gammalooprs/src/integrands/process/cross_section/mod.rs:1577`.
The bracket is in the compact coordinate `u=t/(1+t)`, near `t=5.2e6`.
Root-isolation errors propagate through the stability-level call at
`crates/gammalooprs/src/integrands/process/mod.rs:3783`; they do not return an
unstable value that the ordinary precision escalation can retry.
The failing sample's full coordinates were not retained in the display log.

The ordered comparison is: same model, couplings, projector, graph generation,
thermal CFF construction, residue selection, UV subtraction, and evaluator
construction; then different external momentum at numerical evaluation, which
changes which radial roots are present. The spacelike standalone JSON is
**exactly identical** to the saved timelike thermal standalone JSON, including
both graphs, all expressions, orientation tables, function maps, and residue
orders. Thus the change did not generate a different set of graphs or drop
cut groups during preprocessing. LO and vacuum checks also exclude a simple
normalization or incorrectly applied external-momentum setting.

The NLO mismatch is not yet reduced to a graph, cut, or residue. In particular,
its relation to the separate root-isolation failure is unproven. An understated
Monte Carlo uncertainty, an evaluation/representation issue, or a reference
matching issue has not been excluded. This does not demonstrate a physical
spacelike divergence or establish that threshold subtraction is required.
The paper has a finite fixed-order answer at this hard point. Near the light
cone, however, its Eq. (5.10) exhibits a logarithmic singularity and the
Introduction explains the need for thermal resummation. That is a separate
physical limitation.

## Reproduction

The optional `spacelike` command block installs the external momentum before
integrand preprocessing. The positive `e_cm=sqrt(15)` remains a sampling scale;
it is not a rest-frame energy for this spacelike momentum.

```sh
./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-spacelike-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model spacelike lo nlo export inspect integrate \
  -D samples=100000 -D batch=25000 \
  -D output=/tmp/gammaloop-qcd-vector-spacelike

./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-spacelike-seed-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model spacelike nlo integrate \
  -D samples=300000 -D batch=50000 -D seed=271828 \
  -D output=/tmp/gammaloop-qcd-vector-spacelike-seed

./gammaloop --dev-optim --read-only-state \
  -s /tmp/gammaloop-qcd-vector-spacelike-vacuum-state -n -L off \
  examples/cli/thermal_production/qcd_vector/run.toml \
  run model spacelike lo nlo integrate -D mode=vacuum \
  -D lo_samples=1000 -D samples=1000 -D batch=1000 \
  -D output=/tmp/gammaloop-qcd-vector-spacelike-vacuum

/tmp/gammaloop-gauge-reference/reference 1 4 0.0001
/tmp/gammaloop-gauge-reference/reference 1 4 0.00001
```

Build the author-code driver as described in [README.md](README.md).
[spacelike_results.json](spacelike_results.json) preserves the complete result
records, explicit incomplete-run status, reference values, and solver error.
No core implementation or tests were changed for this diagnostic.
