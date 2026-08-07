# GL297 correlated soft/threshold cure

This analysis direct-imports the directive-dressed `GL297.dot` together with a
local directive-free baseline, generates only the verified orientation for
each, and retains all nine process-valid cuts.  Its four explicit records
replace only the genuinely two-loop legacy threshold CTs:
`(5,7)` and `(5,11,13)` under cut `(2,4,9)` are solved in `[5]`, while `(3,9)`
and `(3,11,12)` under cut `(2,6,7)` are solved in `[3]`.  Every other
association uses the implicit legacy default.

Run from the repository root:

```bash
MallocNanoZone=1 CARGO_INCREMENTAL=0 cargo run -p gammaloop-api --bin gammaloop -- \
  --clean-state -o -l info -L off -p none \
  examples/cli/epem_a_ttxh/NNLO/analyses/gl297_correlated_soft_threshold/analysis.toml \
  run generate approach_e12 approach_e13

python3 examples/cli/epem_a_ttxh/NNLO/analyses/gl297_correlated_soft_threshold/plot.py
```

`inspect_correlated_positive` reproduces the seven positive-λ momentum-space
samples used by the Rust cure regression.  `approach_e12` and `approach_e13`
center the affine paths exactly at the simultaneous soft/threshold point,
scale the displayed parameter as `λ=10^-2 t`, and use `--skip-midpoint` so that
the undefined `λ=0` point is not evaluated.  Each branch therefore samples
`|λ|` logarithmically from `10^-6` through `10^-2` for the threshold-off,
legacy, and forced prescriptions.

Here “threshold off” means the complete selected GL297 graph and orientation,
summed over all nine required LU cuts with the numerator, flux, local UV,
integrated UV, and soft counterterms unchanged, but with only threshold
subtraction disabled.  It is neither a bare Feynman graph nor the sum over all
166 process graphs.  Older generated files and plots called this curve
“bare”; the less ambiguous label is used here.

With `MT=173`, `K* = sqrt(500^2-MT^2) (3,-4,12)/13`.  The two generation-LMB
paths `(K0,K1,K2,K3)=(e4,e5,e6,e11)` are:

```text
e13: K1=K*+λ a, K3=K1+λ d, q13=-K1+K3=λ d
e12: K3=K*+λ a, K0=K2+K3-K1+λ d, q12=K0+K1-K2-K3=λ d
a=(47,29,-31), d=(211,-157,193).
```

Both relevant threshold equations and the soft momentum scale linearly.  The
validated fitted powers are threshold-off approximately `λ^-1`, legacy
approximately `λ^-3`, and forced approximately `λ^-1` in each direction.
`plot.py` produces `e12_cure.pdf` and `e13_cure.pdf`.  Each multipage report
opens with a mirrored double-sided logarithmic comparison of the threshold-off,
legacy, and forced prescriptions, then shows signed event closure, exact
one-sided L+I variants, and ordered LL/LI/IL/II components.  This makes the
legacy excess two inverse powers, the forced restoration, and the resolved
threshold decomposition directly visible without duplicating the original
integrand.
