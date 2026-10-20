= Independent analytic LO targets

All values use $T=1$, $N_c=3$, one massless unit-charge flavor,
with $e^2 Q_d^2$ factored out. The observables obey
$rho_V=2rho_T+rho_L$. They are leading order and have no strong-coupling
multiplier. The 35 benchmark points have 70 targets in `targets.json`;
`rho_lo` and `value` are aliases. No NLO master or NLO reference evaluation
is used.

The source is Jackson, arXiv:1910.07552v2, and the pinned author code
`v0.9-spf`, Zenodo DOI 10.5281/zenodo.3478144, commit
`dca510d2d2c8bee2288345e972ee3a54b37cdbd5`.
The archive and individual source-file hashes are recorded in
`provenance.json`. The archive SHA256 was checked against its local bytes.
The vector and density LO expressions are `examples/dileptons.cpp:76,164`,
and the thermal moments are `src/therm.cpp:61-132`.
The density conversion is
$rho_L=-(E^2-k^2)/k^2 times "author_Rho00.lo"$.

Instead of the author's truncated trilogarithm routine, the implementation
integrates equivalent smooth one-dimensional Fermi phase-space moments.
Let $n_F(p)=1/(exp(p)+1)$ and $K^2=E^2-k^2$.
For $E>k$, use $r=(E+k x)/2$ and
$W(x)=1-n_F(r)-n_F(E-r)$:
$
rho_T = N_c K^2/(32 pi) integral_(-1)^1 (1+x^2)W(x) dif x,
quad
rho_L = N_c K^2/(16 pi) integral_(-1)^1 (1-x^2)W(x) dif x.
$
For $0<E<k$, let $a=(k-E)/2$, $D(p)=n_F(p)-n_F(p+E)$,
and $t=2p+E$:
$
rho_T = -N_c K^2/(8 pi k^3) integral_a^infinity (t^2+k^2)D(p) dif p,
quad
rho_L = N_c K^2/(4 pi k^3) integral_a^infinity (t^2-k^2)D(p) dif p.
$
In particular, spacelike transverse targets are positive and longitudinal
targets are negative. At rest the exact isotropic result is
$rho_T=rho_L=N_c E^2 tanh(E/4)/(12 pi)$.
At $E=0$ and on the exact light cone the LO functions vanish.

Gauss–Legendre rules with 64 and 96 nodes agree to binary64 precision.
The spacelike integral is split into smooth panels up to $p=a+64$;
`validation.json` includes an explicit bound on the omitted exponentially
suppressed tail. Per-target error estimates include refinement differences,
16 ulps, and the tail bound; they are numerical estimates rather than
rigorous total roundoff bounds.

Validation checks all 64 previously saved author LO targets, the closed
analytic vector trace at all 35 points, spacelike signs, and all eight rest
points against an independent small-$k^2$ extrapolation of the phase-space
formulas. Differences from old author values are below $4.4 times 10^(-6)$
relative and reflect the pinned special-function approximation. Full
per-point comparisons are retained in `validation.json`.

Reproduce using a Python 3 interpreter with only the standard library:
```sh
python3 /tmp/qcd-fig15-lo-20261020/reference/analytic_lo.py
```
The generator reads the pinned author files and archive, plus the earlier
64-row LO target table at the explicit scratch paths recorded in the script
and provenance. Generated targets are independent of the native benchmark
binary and its sampling statistics. The figure renderer uses the original
extracted published LO paths; these analytic targets are for precision and
residual comparisons.
