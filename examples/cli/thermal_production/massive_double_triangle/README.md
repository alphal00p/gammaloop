# Massive thermal double triangle

This benchmark integrates an isolated two-loop scalar self-energy graph with
both virtual two-particle cuts and real three-particle cuts, including their
thermal crossings. It uses the double-triangle integral in Blaizot and Reinosa,
[hep-ph/0406109, Fig. 5 and Eqs. (58)–(59)](https://arxiv.org/pdf/hep-ph/0406109).
The comparison numbers are our independent evaluation of the published
arbitrary-mass representation, rather than numbers tabulated in that paper.

The three real scalars have masses `MX=1`, `MH=0.6`, `ML=0.25`. The interaction
is `-lam/2 * (X*H^2 + L*H^2 + X*L^2)`, with `lam=1`; each cubic vertex is
`-i*lam`. All chemical potentials vanish. The imported graph has factor `1/2`
and edges `a-b:H`, `a-c:H`, `b-c:H`, `b-d:L`, `c-d:L`, with X entering at a
and leaving at d. This is one labelled graph contribution, not the complete
order-lam^4 correction of the scalar model.

At external momentum `(1,0,0,0)`, the LL virtual cut is open while HH and the
all-outgoing HHL cuts are closed. The real thermal channels are `X+L -> H+H`
and `X+H -> H+L`, with distinct interference terms from the two three-line cuts.
Their summed contribution is negative and cancels part of the virtual result.
These are interference contributions, so the separate real part need not be
positive.

All surviving cuts are free of physical threshold poles for this choice of
masses and momentum. The virtual triangle is below its thresholds. The real
cut's uncut heavy propagators have denominator magnitude at least `0.2`, and
the potentially dangerous light denominators at least `0.0625`.
[reference_derivation.md](reference_derivation.md) gives the bounds and the
independent cut formulas. All masses are nonzero. Both threshold and UV
subtraction are disabled; no experimental subtraction is needed here.

## Reproduce

From the repository root, after building the CLI:

```sh
./gammaloop --read-only-state -s /tmp/gammaloop-massive-double-triangle-state -n -L off \
  examples/cli/thermal_production/massive_double_triangle/run.toml \
  run generate export inspect integrate lower_temperature

python3 examples/cli/thermal_production/massive_double_triangle/reference.py
python3 examples/cli/thermal_production/massive_double_triangle/reference.py \
  --temperature 0.25 \
  --output examples/cli/thermal_production/massive_double_triangle/reference_t025.json
```

The native card uses the `SummedFunctionMap` evaluator, two integration workers,
one million samples per process, and seed `271828`. It integrates vacuum,
zero-temperature equilibrium, and the inclusive/virtual/real equilibrium
pieces at T=0.5, then the three thermal pieces at T=0.25. All use the same
external momentum. Process filters select cut supports; thermal crossing is
handled by the equilibrium LU implementation. No bath-averaging driver is used.
The independent reference requires NumPy and uses orders 32–512 of tensor
Gauss-Legendre quadrature after analytic cut-delta integration.

The `export` block saves generated graph/CFF diagnostics under `/tmp`; the
`inspect` block saves three point evaluations there. Integration workspaces are
`/tmp/gammaloop-massive-double-triangle-final` and
`/tmp/gammaloop-massive-double-triangle-t025`. The run card uses `--restart`.

## Normalization and results

GammaLoop integrates the graph's spectral relaxation contribution
`Gamma_rho = -Im(Pi_R)/MX`. The production source for an initially empty X
population is `S = n_B(MX) * Gamma_rho`; its rate per `d^3p` is `S/(2*pi)^3`.
This conversion is applied after the native integration. It is not an extra
bath phase-space integral.

At **T=0.5**, one million samples per process:

| Contribution | GammaLoop | Independent reference | Difference / MC error |
|---|---:|---:|---:|
| vacuum | 0.00010978005 ± 9.268e-08 | 0.00010980937553 | -0.316 |
| thermal | 0.0003057930053 ± 6.797e-07 | 0.000305765194062 | +0.041 |
| virtual | 0.0003265925453 ± 2.857e-07 | 0.000326690187024 | -0.342 |
| real | -2.079954004e-05 ± 6.323e-07 | -2.09249929626e-05 | +0.198 |
| zero_temperature | 0.00010978005 ± 9.268e-08 | 0.00010980937553 | -0.316 |

At **T=0.25**, one million samples per process:

| Contribution | GammaLoop | Independent reference | Difference / MC error |
|---|---:|---:|---:|
| thermal | 0.0001463104601 ± 1.474e-07 | 0.000146476852203 | -1.129 |
| virtual | 0.0001533661725 ± 1.336e-07 | 0.000153438271538 | -0.540 |
| real | -7.055712461e-06 ± 6.515e-08 | -6.96141933505e-06 | -1.447 |

Errors are native Monte Carlo errors; reference quadrature convergence is
recorded separately. All eight comparisons are within 1.45 MC standard errors.
The T=0.5 real contribution cancels 6.37% of the virtual contribution (reference:
6.40515%). Vacuum and zero-temperature equilibrium coincide exactly. The
inclusive integral agrees with virtual plus real to better than `8e-19` at
both temperatures; these slots share correlated sampling.

The production source at T=0.5 is `(4.78620 ± 0.01064)e-5`, compared with
`4.78576474095e-5`. Across all eight million native evaluations there were zero
NaNs and zero unresolved stability failures. The T=0.5 and T=0.25 runs used
quadruple precision for 0.2073% and 0.4148% of evaluations, respectively, and
never required arbitrary precision.

`native_result.json` and `native_result_t025.json` retain the complete native
results. `comparison.json` contains the comparisons, production-source and
`d^3p` conversions, closure checks, and unrounded stability statistics.

## Closed-cut solver correction

The first probe exposed an existing runtime bug: unsigned cuts used the final
Newton iterate even when it did not solve the cut equation. Closed massive
channels could therefore produce events that violated energy conservation.
The runtime now uses the existing positive-root routine for both unsigned and
mixed-sign cuts. A cut with no positive root contributes no event. No new
solver helper was added.

`closed_cut_reproducer.json` retains the original pre-fix diagnostic with
`MH=1` (different from the final benchmark's `MH=0.6`): it includes spurious HH
and HHL decay events with energy residuals exceeding one. The new unit test
checks a closed massive two-particle cut at rest and after a boost.

`point_check.json` records the corrected final benchmark. It retains one
virtual and six real events, with maximum four-momentum residual
`4.20e-14`. Inclusive minus virtual minus real is `4.24e-21`; all three
evaluations pass the rotation stability check in double precision.

## Validation

The following checks passed, without changing existing test expectations:

- 30 Esurface/Newton-solver unit tests, including the new closed-cut regression.
- Existing scalar virtual-triangle interference regression, covering both
  virtual mirrors below and above threshold and causal threshold counterterms.
- Equilibrium LU residue-selection regression.
- Cargo format/check, clippy for the core and API crates, and Python ruff checks.

The broader vacuum and root tests can be reproduced with:

```sh
cargo nextest run -p gammalooprs --lib --cargo-profile dev-optim \
  -P test_gammaloop --test-threads 2 \
  -E 'test(cff::esurface::tests::) or test(utils::newton_solver::tests::)'
cargo nextest run -p gammaloop-integration-tests --test test_runs \
  --cargo-profile dev-optim -P test_gammaloop --test-threads 1 \
  -E 'test(=scalar_virtual_phase_conventions::scalar_virtual_triangle_mirrors_match_signed_born_interference)'
```

The previously reported unrelated thermal failures were outside this test
selection. No existing failing tests were modified.


## Above the thresholds

The [above-threshold continuation](above_threshold/README.md) keeps the same
internal masses and temperature and raises MX to 1.6. It opens the HH and HHL
cuts, verifies local cancellation of the virtual threshold poles, and compares
two-million-sample native integrations with independent principal-value
evaluations of the literature formulas.
