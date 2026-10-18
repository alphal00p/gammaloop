= UV-renormalized massive scalar self-energy insertion

This example restores the exact raised-line, virtual, and ordinary UV-renormalized
thermal benchmark recorded at `a0b2d4191cad17d8dc2f20891bb0f245f0df18b7`.
It uses the #link("../massive_double_triangle/model.json")[same scalar model] with
`M=1.6`, `H=0.6`, `L=0.25`, `T=0.5`, unit cubic couplings, zero chemical potentials,
and `P=(1.6,0,0,0)`. The labelled graph multiplier is one.

```text
D_H(Q)^2 D_H(P-Q) D_H(K) D_L(Q-K), D_m(Q)=1/(Q^2-m^2+i0).
```

The two serial `H` edges share momentum `Q`. The repeated shell is outside the
inserted bubble's spectral support:
`(H-L)^2=0.1225 < H^2=0.36 < (H+L)^2=0.7225`.
The probe also avoids the `HH` and `HHL` thresholds, 1.2 and 1.45.
The virtual `HH` contribution and real `HHL` continuum, including thermal crossings,
are evaluated through the shared cut and Laurent machinery.

Only the inserted bubble is UV divergent; the full graph has degree -2.
The ordinary local subtraction and integrated finite addback implement the MUV
prescription with Vakint MSbar normalization, at physical scale `mu_r=1.6`:

```text
B_MSbar(s)=-1/(16*pi^2) integral_0^1 dx
 log([x*H^2+(1-x)*L^2-x*(1-x)*s-i0]/mu_r^2).
```

Medium vacuum subtraction is explicitly false. Subtracting the vacuum *rate* would
not remove this subdivergence: the vacuum bubble divergence multiplies a thermal outer
loop. No thermal counterterm is introduced. Threshold subtraction remains disabled;
thermal infrared treatment and resummation remain out of scope.

== Reproduction

```sh
./gammaloop -n examples/cli/thermal_rates/massive_self_energy/run.toml run generate inspect export save
./gammaloop -n -s /tmp/gammaloop-thermal-self-energy-state run integrate -D samples=250000 -D workspace=/tmp/thermal-self-energy-integral
python3 examples/cli/thermal_rates/massive_self_energy/thermal_reference.py --output /tmp/thermal-self-energy-reference.json
python3 examples/cli/thermal_rates/massive_self_energy/vacuum_reference.py
```

The `save` block persists the generated processes and command history for later
integrations without regeneration. Top-level `-n` disables the additional save on exit;
it does not disable the explicit `save` block. The card builds full thermal, vacuum, cold
zero-density, raised `HH`, continuum `HHL`,
and auxiliary-scale control processes. `thermal_aux` changes `m_uv=1` to `1.4` and
localization scale `1` to `0.7`, keeping every physical parameter fixed. Its integral
must agree with `thermal`; its pointwise integrand need not. `thermal_bare` is generated
only as a pointwise control and is deliberately omitted from the integration block.
The reference scripts require NumPy. The vacuum script writes its adjacent reference JSON;
the thermal script accepts `--output`, `--temperature`, and `--mu-r`.

`export` retains the standalone coefficients, including the integrated localization
kernel and nonzero `log(mu_r^2/m_uv^2)` coefficient. This checks that the integrated
addback is present. `samples`, `batch`, `cores`, `workspace`, and `output` are card
substitutions. The default run cap is 250,000 per process and seed `271828`; its 0.5%
target is a stopping criterion rather than a development prerequisite.

== Independent reference and historical results

The reference uses #link("https://arxiv.org/pdf/hep-ph/9409388")[Bauberger et al.,
Eqs. (3), (13)–(15)] for mass differentiation and bubble dispersion, and
#link("https://arxiv.org/pdf/hep-ph/0302098")[Nishikawa, Morimatsu and Hidaka,
section III] for the thermal bubble. The insertion reduction is derived in
#link("reference_derivation.typ")[the reference notes]. These numerical values are
independent evaluations of the literature ingredients, not published tabulated numbers.

#table(
  columns: 3,
  table.header([Contribution], [Vacuum kernel], [Full kernel at `T=0.5`]),
  [Raised `HH`], [`6.118053358288014e-5`], [`-8.825872261661233e-5`],
  [`HHL` continuum], [`7.930999160990028e-6`], [`2.283722269045362e-4`],
  [Inclusive], [`6.911153274387017e-5`], [`1.401135042879238e-4`],
)

The old native thermal estimate was `1.4041890505e-4 +/- 7.137e-7` from two million
samples, a 0.43-error difference. The six runs totalled 7.5 million evaluations;
all comparisons were within 1.13 errors, with zero NaNs or unresolved stability failures.
Changing both auxiliary scales shifted the central value by only `3.99e-8`.
These are historical results, explicitly identified in
#link("historical_validation.json")[the retained comparison record]. Current validation
must be recorded independently; historical evidence is not a new checkout result.

The old pinned symbolic evaluator initially could not vectorize native hyperbolic
functions in raised residues and used an exponential workaround. That workaround is
not imported: current shared hyperbolic derivatives are already exercised by raised
thermal tests. The exact bare thermal point was `2.0918100443146576e-5`; agreement with
its split-mass limit provided a separate control. Literal/composite mass precision and
shared routing must agree before a pointwise comparison is interpreted physically.

== Current checkout validation

All 116 selected native regressions pass, including serialized-state reload, exact
raised-line and UV settings, physical event closure, and the finite integrated-addback
scale response. The current CLI (`a7752162aa0b917e…`) also completes all six 25,000-point
integration controls and a separate 100,000-point thermal comparison with seed `314159`:

#table(
  columns: 4,
  table.header([Contribution], [Samples], [Native estimate and error], [Independent target]),
  [Full thermal], [100,000], [`1.35989609e-4 +/- 3.40253e-6`], [`1.40113504e-4`],
  [Auxiliary scales], [100,000], [`1.35848099e-4 +/- 3.53726e-6`], [`1.40113504e-4`],
  [Raised `HH`], [100,000], [`-9.11693078e-5 +/- 3.35659e-6`], [`-8.82587226e-5`],
  [`HHL` continuum], [100,000], [`2.28239495e-4 +/- 2.04348e-6`], [`2.28372227e-4`],
  [Vacuum / cold zero density], [25,000 each], [`6.60382933e-5 +/- 3.62448e-6`], [`6.91115327e-5`],
)

The four 100,000-point comparisons are within 1.22 reported standard errors; changing
the two auxiliary UV scales shifts the central value by `1.42e-7`. The full-rate Monte
Carlo error is 2.50%, so this run does not claim the card's 0.5% stopping target.
Both completed runs have zero final NaNs and zero unresolved stability failures.
Warnings about unsuccessful low-precision attempts precede successful precision retries
and are recorded separately from these final failure counters.

#link("validation.json")[The validation record] retains exact commands, executable and
result hashes, all channel statistics, the earlier pilots, and the independent references.
The two seeds' estimates are reported separately. The 100,000-point run trains on the
full thermal slot and samples all four contributions with correlated coordinates;
pointwise channel closure is also checked directly by the native regression.

Correlated cube samples need not describe identical physical momenta in each slot.
The parent LMB is `[5,1]` throughout, but the default sampling channel is `lmb(2,5)` for
the full and `HH` slots and `lmb(4,2)` for `HHL`. Selecting `lmb(2,5)` in all three
restores closure at the checked cube point to `3.39e-21`. The reported integrated means
use ordinary sample averages; their small closure residual is not an inverse-variance
averaging effect or evidence of a physical coefficient discrepancy.
