= Figure 15 author reference

The independent target is the unmodified author implementation accompanying
G. Jackson, _Two-loop thermal spectral functions with general kinematics_,
#link("https://arxiv.org/abs/1910.07552v2")[arXiv:1910.07552v2].
Use #link("https://doi.org/10.5281/zenodo.3478144")[Zenodo3478144],
version `v0.9-spf`, source commit
`dca510d2d2c8bee2288345e972ee3a54b37cdbd5`. The archive SHA256 is
`207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb`.
`build.sh` downloads that archive, verifies it, and compiles `driver.cpp`
against the original `examples/dileptons.cpp` and masters. It changes
neither master integrands nor integration domains.

The driver selects the existing `Rho_V` or `Rho_00` assembly. Output is
`E,k,tolerance,lo,nlo,seconds`; `nlo` is the full vacuum-plus-thermal
coefficient of $g^2$, with $N_c=3$, $C_F=4/3$, one unit-charge flavor,
$T=1$, and zero chemical potential. The class fields are
$op("Rho_V")=-op("Im")Pi^mu_mu$ and
$op("Rho_00")=-op("Im")Pi^(00)$. Thus
$rho_L=-(K^2/k^2)op("Rho_00")$ and
$rho_T=(op("Rho_V")-rho_L)/2$. Figure15 shows LO plus nine times NLO,
with $g=3$ and the electromagnetic coupling factored out.

== Build and one-point evaluation

From the repository root, enter an environment with a C++11 compiler,
GSL (including `gsl-config`), curl, unzip, Python3, and Bash. For example:

```sh
nix shell nixpkgs#gcc nixpkgs#gsl.dev nixpkgs#curl nixpkgs#unzip nixpkgs#python3 --command bash
bash examples/cli/thermal_rates/qcd_vector/figure15/reference/build.sh /tmp/figure15-author
python3 examples/cli/thermal_rates/qcd_vector/figure15/reference/reference_cli.py --binary /tmp/figure15-author/reference --energy 1.05 --momentum 1 --tolerance 1e-5 --basis both --timeout 1800 --cache-directory /tmp/figure15-author/cache --output /tmp/figure15-author/point.json
```

The build directory retains the archive, exact build command, compiler and
GSL versions, and source/binary hashes. The upstream renamed example main
emits a missing-return warning; this entry point is never invoked.
`QCD_FIG15_REFERENCE_BINARY` and `QCD_FIG15_REFERENCE_CACHE` can supply
CLI defaults. No scratch path is embedded in the source or required at runtime.

The CLI supports `both`, `vector`, and `density` bases. At $k=0$, upstream
coordinates are singular, so the CLI evaluates only the vector assembly
at $k/E=0.05,0.025,0.0125$, extrapolates quadratically in $k^2$, and uses
$rho_T=rho_L=rho_V/3$. It reports the adjacent-fit difference and an
additional quartic extrapolation. Direct tiny-$k$ density division is
unsuitable because special-function cancellation is amplified by $1/k^2$.
The rest LO result uses its exact analytic limit.

The cache stores raw class fields and identifies each record by the exact
floating-point inputs, tolerance, component, and verified archive hash.
Cached components remain available after interruption. The total timeout
covers all source calls for a point. `--cpu` sets inherited CPU affinity;
OpenMP and BLAS thread counts are restricted to one. SIGTERM/SIGHUP kills
and reaps the active source child before emitting an incomplete result.

== Grid scheduler contract

`reference_grid.py` covers 35 points: $k/T=0,0.5,1,1.5$ and
$E/T=0.25,0.5,...,2$, excluding $E=k$, plus $E=k plus.minus 0.05$
for each nonzero momentum. It uses at most six single-CPU workers on
CPUs236--241 by default. `--workers`, `--cpus` (comma-separated IDs or
half-open ranges), and `--deadline-epoch` control these bounds.
All subprocesses inherit the controller's process group; none starts a new
session or process group. The full controller must terminate that whole
group at its hard deadline.

```sh
python3 examples/cli/thermal_rates/qcd_vector/figure15/reference/reference_grid.py --binary /tmp/figure15-author/reference --cache-directory /tmp/figure15-author/cache --output-directory /tmp/figure15-author/grid --workers 6 --cpus 236:242 --deadline-epoch DEADLINE_EPOCH
```

`--dry-run` writes initial metadata and known pilot targets without
starting any source calculations. `targets.json` is an atomically replaced
list with `E`, `k`, `component` (`T` or `L`), `rho_lo`, `rho_nlo`,
`rho_g3`, and reference provenance. `status.json`, per-point result and
stdout/stderr files, and `jobs.jsonl` retain progress and failures. A failed
reference preserves an existing pilot target and never blocks native work.
The checked-in `pilot-targets.json` seeds available targets before any
new reference calculation finishes; it is a generated numerical record.

== Reference accuracy

The author disables GSL's fatal handler and does not expose quadrature
return codes. Tolerance variation and momentum extrapolation are empirical
convergence checks rather than rigorous error bounds. The pilot compared
absolute tolerances $10^(-4)$ and $10^(-5)$; the spacelike point
$(E,k)=(0.5,1)$ reached a 900-second cap at the finer tolerance, so its
retained target explicitly records $10^(-4)$. Its vector-only result was
independently repeated. Other Figure15 pilot changes are below 0.0051%
relative in either NLO component. The hard $(4,1)$ control is outside
Figure15's energy range and received a separate $10^(-6)$ refinement
because density projection amplifies the coarser tolerance change.

== Published vector curves

The overlay's smooth lines are extracted directly from Figure 15 on page 33
of #link("https://arxiv.org/pdf/1910.07552v2")[the published arXiv PDF].
They are separate from the independently computed numerical targets above.
The calibrated PDF SHA256 is
`b92fd775f7b5ebb15ea382925acc6e89b6ff1f1f9eec1898432a319c259a7d65`.

`extract_published_curves.py PAPER.pdf published-curves.json --mutool MUTOOL`
uses MuPDF's SVG export, interprets the original move, line, and horizontal
path commands, and maps the resulting vertices with the printed tick marks.
It extracts fourteen curves into fifteen paths: eight LO-plus-NLO curves
for four momenta in two panels, and six gray LO curves for nonzero momenta.
The blue transverse curve retains its two original clipped segments.
There are 4,187 vertices; no new interpolation or smoothing is applied.
The JSON records source hashes, axis anchors, affine fits, path hashes,
and the original path segmentation.

These are graphical references, with the precision and clipping of the
published plot. They neither replace the numerical targets nor change
Monte Carlo errors, residuals, precision acceptance, or the native marker
normalization (analytic LO plus nine times native NLO at $g=3$).
