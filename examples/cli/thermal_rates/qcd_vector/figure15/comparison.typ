= QCD vector polarizations for Figure 15

This card projects the native thermal QCD calculation onto the transverse
and longitudinal spectral functions in Figure 15 of
#link("https://arxiv.org/pdf/1910.07552")[G. Jackson, _Two-loop thermal
spectral functions with general kinematics_]. The figure uses
$k/T=0,0.5,1,1.5$, $g=3$, horizontal coordinate $k_0/T$, and vertical
coordinate $rho_(T,L)/(e T)^2$. The electromagnetic coupling is factored
out. The reference includes $N_c=3$ and one unit-charge massless flavor;
there is no additional color or $2pi$ conversion.

== Native projections and normalization

The model and subtraction settings follow the parent
#link("../benchmark.typ")[hard timelike vector benchmark]: one down-quark
flavor, $e=g=1$, $Q_d=-1/3$, thermal equilibrium at $T=1$, local
three-dimensional UV subtraction, integrated MSbar counterterms, and
physical renormalization scale $mu_r=4$. Auxiliary UV and localization
scales are also four. Threshold and soft counterterms remain disabled.
The reference scale must match when comparing the isolated NLO correction.

Write $K=(E,0,0,k)$ in the bath rest frame, $n=(1,0,0,0)$, and
$s=K-E n$, so $s^2=-k^2$. The existing global-prefactor projector
interface supplies

$ P_V^(mu nu) = -g^(mu nu), $
$ P_T^(mu nu) = 1/2(-g^(mu nu)+n^mu n^nu+s^mu s^nu/s^2), $
$ P_L^(mu nu) = -n^mu n^nu-s^mu s^nu/s^2. $

The card expresses $n$ through the existing temporal delta tensor and
$s$ through `Q3(0,...)`. Their sum satisfies $2P_T+P_L=P_V$ before
generation, including for spacelike $K$. Each projected numerator remains
factorized through the shared generation and subtraction pipeline.
Longitudinal results are integrated directly; no subtraction of separately
estimated Monte Carlo integrals is required.

For a native result $I_X$ from any LO or NLO slot,

$ rho_X = 9 E I_X, quad X in {T,L,V}. $

This removes $e^2 Q_d^2$ and the one-particle flux normalization.
The factor 36 belongs only to $E=4$. The extra diagnostic `00` slots
have projector $n^mu n^nu$ and give the positive native charge-density
convention. After integration, current conservation gives
$rho_L=(K^2/k^2)rho_(00)$; the author's `Rho_00` class instead returns
the opposite sign. At NLO the density relation need not hold pointwise,
whereas the direct projector trace relation does.

At $k=0$, rotational symmetry gives $rho_T=rho_L=rho_V/3$.
Generate and integrate only the trace slots there: the direction-dependent
projectors contain $1/s^2$ and are undefined at rest. For $k>0$ the
projectors contain no $1/K^2$ denominator. This does not supply light-cone
resummation to the underlying unresummed two-loop calculation.

Native NLO slots are the isolated $g^2$ coefficient because generation
uses $g=1$. Figure 15 compares the full fixed-order expression
$rho_X^"LO"+9rho_X^"NLO"$, divided by $T^2$ for the plotted units.
Both that sum and the isolated coefficient should be retained when
assessing agreement; LO+NLO is not an NLO-correction acceptance target.

The external current half-edge indices are `0,5` at LO and `0,11` at NLO,
for both NLO masters `GL0` and `GL2`. These are the actual sewn
cross-section indices of this process, not the `0,1` amplitude indices.
The `export` block makes those assignments inspectable in DOT output.

== Reproduction

Run from the repository root using the pinned development environment.
The explicit state directory keeps generated state outside the checkout:

```sh
nix develop --command ./gammaloop -n \
  -s /tmp/gammaloop-qcd-figure15-state \
  examples/cli/thermal_rates/qcd_vector/figure15/direct-projectors.toml \
  run model LO_T LO_L LO_00 LO_V NLO_T NLO_L NLO_00 NLO_V save
```

The blocks create separate processes `vector_T`, `vector_L`, `vector_00`,
and `vector_V` with integrand `LO`, and corresponding `_nlo` processes
with integrand `NLO`. The `00` blocks are optional diagnostics. For rest
points generate only `LO_V` and `NLO_V` instead.

The `point` block changes one existing slot, for example:

```sh
./gammaloop -n -s /tmp/gammaloop-qcd-figure15-state \
  examples/cli/thermal_rates/qcd_vector/figure15/direct-projectors.toml \
  run point save -D process=vector_L_nlo -D energy=1.0 -D momentum=0.5 \
  -D e_cm=0.8660254037844386
```

The explicit `save` persists this change. Apply the same kinematics to
each selected slot; use `order=LO` for LO.
Here `e_cm` is a positive sampling scale, chosen as $sqrt(abs(K^2))$
away from the light cone. Changing runtime kinematics reuses the generated
tensor expressions.

The `validate_anchor`, `validate_timelike`, and `validate_spacelike`
blocks reproduce the 20,000-sample LO projection controls at $(E,k)=$
$(4,1)$, $(1,0.5)$, and $(0.5,1)$. They require all four LO processes.
The `inspect_nlo` block reproduces the fixed NLO point, with all four NLO
processes still at their default $(4,1)$ kinematics. Output paths use
`-D output=/tmp/my-figure15-check`; the LO integration worker count uses
`-D cores=8` by default.

The `pilot` block integrates the direct NLO $T,L,V$ slots together at
their current kinematics, initially with 10,000 samples per slot. Options
`samples`, `batch`, `seed`, `cores`, and `output` select a finite pilot.
Its workspace uses `--restart`; choose a distinct output prefix for each
independent point or seed. The card requests 0.5% relative error and keeps
the default precision ladder and generic Euler rotation checks. An
external campaign controller must enforce the shared worker cap and
wall-time deadline; this card alone is not a deadline controller.
The agreed campaign ceiling is 60 workers in total. Pilot results and
cost estimates are reported before starting the full scan, whose hard
wall-time limit is twelve hours.

== Preserved validation

#link("projector-validation.json")[The numerical certificate] records the
accepted native projections on source commit
`116767bd3ed30491a2b15559d228171d20082abe` using the immutable `itphlies`
CLI with SHA256
`9f4f7bc9b5f5d1caffa8235cb88b93429b90677aabe9e05420dd155e03e8de24`.
It preserves the original scratch-card identity and build command;
the durable card parameterizes output paths and retains the same projector
expressions. No Rust changes or relaxed precision settings were needed.

At each of the three LO points, 20,000 samples per slot finished with
zero unresolved nonfinite or unstable samples. Direct longitudinal
results agree with independent massless phase-space integration within
0.64 quoted errors. The LO direct/density relation agrees pointwise to
$5.5 times 10^(-20)$ and after integration to $5.6 times 10^(-17)$ in
native units. Direct trace closure agrees to $1.7 times 10^(-19)$ at
the LO points and $2.9 times 10^(-19)$ at the preserved NLO point.
These certify projection and normalization; integrated NLO convergence
requires the separately recorded pilots and subsequent scan.

Earlier fixed Cartesian projector trials were excluded because a fixed
spatial axis does not follow the default rotated external momentum.
The accepted covariant projectors preserve the original generic rotation
checks and retain no axis-specific stability overrides.
