# Powered scalar-dot freshening is required by whole-numerator FORM reduction

Current behavior and complete-value reproducer · 11 September 2026

The final freshening pass in `Vakint::convert_from_dot_notation` has a concrete correctness requirement. Disabling only that pass changes the public tensor-reduction result for a factorized scalar square when `project_onto_tensor_integrals=false`. Keep the pass. Its measured conversion overhead is real, but removing it is unsafe for the supported whole-numerator route.

## Independent oracle and observed result

Use the numerator `(dot(k(1),p(1)) + dot(k(1),p(2)))^2` with one-loop vacuum topology `I1L(muvsq,1)`. Write the two external vectors as `p` and `q`. Vacuum isotropy gives `<k_mu*k_nu> = g_mu_nu*k^2/D`, hence the complete angular average is

```
k^2*(p^2 + 2*p·q + q^2)/D,  D = 4 - 2*eps.
```

Disabling freshening gives

```
k^2*(p^2 + q^2) + 2*k^2*p·q/D.
```

Both diagonal terms have lost `1/D`. The complete residual is `k^2*(p^2+q^2)*(D-1)/D`. The example is a supported vacuum-projection algebra input; it is not a full physical cross-section acceptance.

| Public tensor-reduction mode | Current conversion | Freshening disabled |
| --- | --- | --- |
| Projected tensor integrals | Exact oracle match | Exact oracle match |
| Whole numerator | Exact oracle match | Nonzero residual shown above |

All four public API executions complete successfully. The wrong result is silent. Both direct runs of the unchanged production FORM template corroborate this matrix at the backend boundary.

## Why the scalar round trips miss it

The Rust Lorentz parser contracts a scalar power's base before applying its exponent. Reusing the already contracted base therefore gives correct scalar round trips.

Whole-numerator reduction instead sends the indexed representation to FORM. Without freshening, the expression contains a power such as `(k(mu)*p(mu) + k(nu)*q(nu))^2`. FORM's production expansion exposes repeated copies with the same index. Its vector-pair reduction can contract the loop vectors with each other before the required angular average. Giving each copy independent dummy indices prevents that incorrect contraction.

Projected mode decomposes scalar monomials before this handoff. Direct powers of individual dots already receive separate dummy contractions through the converter's existing `dot_pow` loop. It consequently passes this example even with the final sum-power pass disabled. The complete pipeline comparison locates the failure at the indexed FORM handoff; topology, Symbolica dependency and angular-average oracle are shared.

## Behavioral regression

`tensor_reduction_of_powered_scalar_sum_matches_angular_average` in [tensor_reduction_tests.rs](../../crates/vakint/tests/tensor_reduction_tests.rs) exercises both public modes against the same independent formula. It compares the complete rational residual after `together().cancel()`. It does not constrain dummy names, copy counts, intermediate storage or expression formatting, and does not expand the input numerator.

The original diagnostic's literal collected-factor comparison was false even for equivalent denominator signs. Its raw outputs remain unchanged in the evidence. A separate checker maps the exact scalar invariants and common topology, then cancels the complete rational residual; it proves three matches and the one wrong whole-numerator result. No native observation was rerun to repair that representation comparison.

## Reproduction and scope

The [evidence bundle](raised-energy-cff-powered-dot-evidence.tar.gz) retains the public probe, exact three-line intervention, source and binary hashes, commands, original outputs, FORM inputs/resources, and the independent residual checker. The executable is based on `d55a0f3e9d25e5c64bdb9ef9fb57749ecf680390`; the intervention adds an opt-in early return immediately before the final freshening traversal in an isolated checkout. Production logic is unchanged. The workspace adds the outward regression and a comment explaining the FORM requirement.

The public probe runs all four combinations concurrently under an 8 GB process-tree cap, with a 60-second per-process limit and no retries. Compilation is guarded at 15 GB. Native benchmark timings from the motivation audit remain unchanged; these new executions establish correctness, not a workflow speed claim. No PySecDec integration is launched or certified.

Validation on the unchanged production implementation: all 44 selected Vakint unit and tensor-reduction tests pass with four workers, zero retries and disabled snapshot updates. The new regression checks both modes. Formatting, the locked test-target check and Clippy on all Vakint targets pass. The retained test log distinguishes build time from the test execution time.

```
cargo fmt --all -- --check
cargo check -p vakint --test tensor_reduction_tests --profile dev-optim --locked
cargo nextest run -p vakint --lib --test tensor_reduction_tests \
  --cargo-profile dev-optim --locked --test-threads 4 --retries 0
cargo clippy -p vakint --all-targets --profile dev-optim --locked
```
