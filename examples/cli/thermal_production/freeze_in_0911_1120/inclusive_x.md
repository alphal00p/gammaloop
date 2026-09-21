# Forward-X thermal production benchmark

The forward-X construction is the appropriate way to obtain a differential
production rate without an external bath-averaging driver. The current native
cross-section implementation generates the necessary thermal CFF, but its LU
residue selection drops the mixed-sign scattering surfaces. The integral it
returns is the thermal decay spectral rate, not the complete inclusive rate.
No Rust implementation or existing test was changed in this investigation.

## Setup and production-rate convention

Use the custom `model.json`: four distinct real scalars with interaction
`lambda X B1 B2 B3`, massless bath fields, and `MX = T = lambda = 1`.
The forward graph at order lambda squared is the two-loop sunset with three
internal bath lines. Its only bath-species content is B1 B2 B3, so generating

```
generate xs X > B1 B2 B3 FI^2==2 [{{2}}] -p thermal -i LO
```

produces the complete forward graph at this order. The explicit two-loop
specification is necessary: the CLI defaults to one forward loop. There is no
internal X line and no need to assign X an equilibrium population.

Define `n_X = integral d^3p/(2*pi)^3 f_X(p)` and
`S(p) = (df_X(p)/dt)|_{f_X=0}`. For an equilibrium bath, detailed balance gives

\[
S(p)=n_B(E_X)\,\Gamma_\rho(p),\qquad
\Gamma_\rho=\Gamma_{\rm abs}-\Gamma_{\rm prod}
=-\frac{\operatorname{Im}\Pi_R}{E_X},\qquad
\frac{d\gamma}{d^3p}=\frac{S(p)}{(2\pi)^3}.
\]

Here the Bose factor at E_X follows from KMS for the bath correlator; it does
not assume an actual thermal X population. See the general relation in
[Redondo and Raffelt, equations (2.2)-(2.3)](https://arxiv.org/pdf/1305.2920).
GammaLoop's one-incoming-particle flux already supplies `1/(2 E_X)`;
no additional flux factor should be applied to its result.

## Where the native calculation loses scattering

The common inputs are the model, graph, constant squared matrix element,
external X momentum, and thermal CFF. The last complete representation is the
saved `inclusive_cff.json`, which has eight orientations: two E-surfaces and
six H-surfaces. For positive on-shell E_X, the physically relevant surfaces are

- `E1 + E2 + E3 - E_X = 0`: decay/inverse decay;
- `Ei + Ej - Ek - E_X = 0`, for the three choices of k: scattering.

The native cross-section path constructs only positive-energy E-surface cuts
in `crates/gammalooprs/src/processes/cross_section.rs:1700`, using
`Esurface::new_from_cut_left` in `crates/gammalooprs/src/cff/esurface.rs:873`.
It then calls `select_esurface_residue` in
`crates/gammalooprs/src/cff/mod.rs:362`.
This is the first boundary that excludes the required scattering surfaces.

The durable before/after evidence is:

| Orientation ID | Canonical denominator | Kind | Nonzero native residue |
|---:|---|---|---|
| 0 | E1+E2+E3+E_X | E | no |
| 1 | -E1+E2+E3+E_X | H | no |
| 2 | E1-E2+E3+E_X | H | no |
| 3 | E1+E2-E3-E_X | H | no |
| 4 | E1+E2-E3+E_X | H | no |
| 5 | E1-E2+E3-E_X | H | no |
| 6 | -E1+E2+E3-E_X | H | no |
| 7 | E1+E2+E3-E_X | E | yes |

`inclusive_standalone.json` contains a single cut group with an explicit
`if(-7+residue_map_id,0,...)` selector. Thus summing all eight orientations still
keeps only orientation 7. Its statistical numerator is exactly

\[
(1+n_1)(1+n_2)(1+n_3)-n_1n_2n_3.
\]

Multiplication by n_B(E_X), on this cut, gives `n1*n2*n3`: production by
`B1+B2+B3 -> X`. The paper's three `Bi+Bj -> Bk+X` channels instead need
orientations 3, 5, and 6. This discrepancy is established before integration;
UV subtraction, threshold counterterms, the numerator, and Monte Carlo
sampling do not explain these missing contributions.

## Native momentum scan

The saved `inclusive_x.toml` runs vacuum and equilibrium modes at
`|p_X| = 0, 1, 3`, each with one million samples and two workers. The table
reports the source S(p), obtained by multiplying the native equilibrium
integral by n_B(E_X). All entries use mass units set by MX=1.

| p_X | Native S(p) | Bose inverse-decay reference | Bose scattering contribution omitted |
|---:|---:|---:|---:|
| 0 | (1.48103 +/- 0.01420)e-3 | 1.4839167951e-3 | 2.7451682566e-3 |
| 1 | (4.27167 +/- 0.04710)e-4 | 4.3583698394e-4 | 1.0654035432e-3 |
| 3 | (1.15407 +/- 0.05459)e-5 | 1.2671830764e-5 | 5.7220506593e-5 |

The inverse-decay comparisons differ by -0.20, -1.84, and -2.07 Monte Carlo
standard errors. Precision deteriorates at the boosted point. No samples were
reported as NaN or unstable. The paired vacuum integrals agree with the known
three-body decay rate `lambda^2 MX^2/(512*pi^3 E_X)` within 1.1 standard errors.
At rest the native unconverted thermal integral is
`0.002544820884 +/- 0.000024397616`, compared with `0.002549787264`.
These checks establish the contribution being calculated; they do not validate
complete thermal production.

The inverse-decay reference integrates the flat massless Dalitz phase space,
weighted by the product of three bath Bose functions. The energy triangle is
split according to its softest particle so that the endpoint singularity is
cancelled by the phase-space width. Its orientation is averaged after boosting
the three energies into the bath frame. Quadrature orders 32 and 64 agree to
better than 1e-9 relative at every quoted momentum.

The independent Bose scattering reference sums three identical channels.
Writing k for the absorbed bath particle's energy in the inverse process,
`q0 = E_X+k`, `Q = |p_X+k|`, and `E_+- = (q0+-Q)/2`, it uses

\[
S_{2\to2}^{\rm BE}(p)=\frac{3\lambda^2 e^{-E_X/T}}{128\pi^3 E_X}
\int_0^\infty dk\,k n_B(k)\int_{-1}^{1}dc\,
[1+n_B(q^0)]\left[1+\frac{2T}{Q}
\log\frac{1-e^{-E_+/T}}{1-e^{-E_-/T}}\right].
\]

This follows by integrating the outgoing two-body phase space analytically in
`X+Bk -> Bi+Bj` and using detailed balance. Orders 64, 128, and 256 agree at the
displayed precision. Adding this column to the inverse-decay reference gives
the complete Bose production source at this order.

## Differential reference for section 6.3

[Hall et al., section 6.3](https://arxiv.org/pdf/0911.1120) counts only the three
scattering channels, with Maxwell-Boltzmann bath occupations and no final-state
Bose enhancement. Under those approximations the differential source is

\[
S_{2\to2}^{\rm MB}(p)=
\frac{3\lambda^2 T^2}{64\pi^3 E_X}\,e^{-E_X/T}.
\]

One way to derive this is to use the inverse process. The massless outgoing
pair has integrated phase space `1/(8*pi)` and the absorbed bath integral is
`integral d^3k/[(2*pi)^3 2k] exp(-k/T) = T^2/(4*pi^2)`.
The flux `1/(2 E_X)`, detailed-balance factor `exp(-E_X/T)`, and three species
channels give the expression above. Integrating only the produced momentum
then recovers the paper's equation (6.24):

\[
\gamma_{2\to2}^{\rm MB}=
\frac{3\lambda^2 m_X T^3}{128\pi^5}K_1(m_X/T)
=4.6098997578035\times10^{-5}\quad (m_X=T=\lambda=1).
\]

The corresponding S(p) values at p=0,1,3 are respectively
`5.5615670926e-4`, `2.5989116652e-4`, and `2.0236337997e-5`.
Divide by `(2*pi)^3` for `d gamma/d^3p`, or multiply S(p) by
`p^2/(2*pi^2)` for the radial distribution `d gamma/dp`.

Two adjustments are needed for an exact comparison to that paper: retain the
mixed-sign scattering cuts, and match its Maxwell-Boltzmann approximation.
A fully inclusive calculation for this massless-bath model also contains
inverse three-body decay, which section 6.3's scattering formula does not
include. In the MB approximation that extra source is
`lambda^2 MX^2 exp(-E_X/T)/(512*pi^3 E_X)`.

## Reproduce and artifacts

From the repository root, using the existing built binary:

```sh
./gammaloop --read-only-state -s /tmp/gammaloop-freeze-in-inclusive-state \
  -n -L off \
  examples/cli/thermal_production/freeze_in_0911_1120/inclusive_x.toml \
  run generate export integrate p1 p3
```

`inclusive_results.json` preserves all three raw native results, and
`inclusive_comparison.json` contains normalizations, references, convergence,
selector evidence, and provenance. `inclusive_cff.json` and
`inclusive_standalone.json` preserve the representation before and after cut
selection; `inclusive_graph.dot` preserves the generated graph. The run uses
the interpreted summed function-map evaluator, with UV and threshold
subtraction disabled. No separate bath-averaging script is used by this card.
