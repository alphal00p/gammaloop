# Independent spectral reference for the self-energy insertion

[thermal_reference.py](thermal_reference.py) computes the renormalized labelled
insertion graph without reading any GammaLoop output. The numerical values
are new evaluations for this benchmark, not numbers quoted from a paper.

The bubble ingredients are supported by
[Nishikawa, Morimatsu and Hidaka, hep-ph/0302098](https://arxiv.org/pdf/hep-ph/0302098):
section III treats arbitrary unequal masses; Eqs. (III.30)-(III.34) give the
vacuum bubble, Eqs. (III.35)-(III.42) its thermal discontinuities, and
Eq. (III.44) the one-dimensional thermal real part. Their calculation uses the
11 component; Eq. (II.14) supplies its relation to the physical discontinuity.
The formulas below instead use the retarded continuation of the positive
Euclidean bubble directly. The insertion reduction below is our derivation.

[Bauberger et al., hep-ph/9409388](https://arxiv.org/pdf/hep-ph/9409388)
provide the vacuum mass-derivative relation in Eq. (3) and bubble dispersion
construction in Eqs. (13)-(15). The independently evaluated vacuum result is
also saved by [vacuum_reference.py](vacuum_reference.py).

## Conventions

Write the positive Euclidean propagator as
`D_a(Q) = 1/(Q_E^2+a)` and the MSbar bubble with fixed internal masses H,L as
`B(Q)`. The integral is

```text
I(P) = sum-integral_Q D_a(Q)^2 D_b(P-Q) B(Q),
a=b=H^2,  P=(M,0),  M=1.6, H=0.6, L=0.25, T=0.5.
```

There is one labelled insertion and graph overall factor 1. With the model's
cubic couplings equal to 1, `Pi_R=-I_R`, so the spectral relaxation-rate
contribution is `Gamma_rho=Im(I_R)/M`. The production source for an initially
empty X mode is `n_B(M)*Gamma_rho`. These are contributions of this graph,
not the complete model's rate at this perturbative order.

The vacuum bubble, including its finite MSbar term, is

```text
B_0(s) = -1/(16*pi^2) integral_0^1 dx
         log([x*H^2+(1-x)*L^2-x*(1-x)*s-i0]/mu_r^2).
```

We use `mu_r=1.6`. The vacuum contribution is retained throughout. Only the
bubble is UV divergent; this subtraction supplies its ordinary vacuum mass
counterterm. No thermal counterterm is introduced.

## Raised two-particle cut

For a trial repeated-line squared mass a, with the bubble and spectator masses
held fixed, define

```text
E_a = (M^2+a-b)/(2*M),  E_b=M-E_a,
q_a = sqrt(E_a^2-a),   phi(a,b)=q_a/(8*pi*M).
```

The bubble has no spectral support on this shell because
`(H-L)^2 < a < (H+L)^2`. The repeated propagator identity
`D_a^2=-dD_a/da` then gives the pole contribution

```text
Im I_HH = -d/da { phi(a,b) * [1+n_B(E_a)+n_B(E_b)]
                 * Re B_R(E_a,q_a) } at a=b=H^2.
```

The derivative acts on the kinematics and the thermal functions. At equal
a and b, the two outer occupation-factor derivatives happen to cancel each
other; the energy/momentum dependence of the internal thermal bubble still
contributes. Dropping that dependence gives the wrong result.

For `s=omega^2-q^2`, let `E_i=sqrt(k^2+m_i^2)` and
`A_i=m_j^2-m_i^2-s`. The finite thermal bubble is

```text
B_T(omega,q) = sum_(i,j) 1/(16*pi^2*q) integral_0^infinity dk
  k*n_B(E_i)/E_i * sum_(sigma=+-1)
  log |[A_i+2*sigma*omega*E_i+2*q*k]
       /[A_i+2*sigma*omega*E_i-2*q*k]|,
(i,j) = (H,L),(L,H).
```

No principal-value singularity is encountered near the stable H shell.
The script therefore evaluates the derivative by complex steps in a, using
the analytic logarithm branches continued from that shell. It compares steps
`1e-4`, `1e-8`, and `1e-20`. The last two agree to floating-point precision.

## Three-particle continuum and thermal crossings

Define `E=sqrt(q^2+b)`, `omega_-=M-E`, `omega_+=M+E`, and
`s_+-=omega_+-^2-q^2`. The continuum part of the insertion's spectral function
is simply the bubble spectral function divided by `(s-a)^2`. Convolution
with the spectator gives

```text
Im I_HHL = integral_0^infinity dq q^2/(4*pi^2*E) * {
  [1+n_B(omega_-)+n_B(E)] * Im B_R(omega_-,q)/(s_--a)^2
  + [n_B(E)-n_B(omega_+)] * Im B_R(omega_+,q)/(s_+-a)^2 }.
```

Here `n_B(-omega)=-1-n_B(omega)` and
`Im B_R(-omega,q)=-Im B_R(omega,q)`. The second term and negative-energy
part of the first term are essential thermal crossings. Both denominators
stay away from zero on the bubble's spectral support for these masses.

For completeness, the script's bubble discontinuity can be reconstructed
without a Matsubara sum. Take positive omega, `s=omega^2-q^2`,
`A=s+H^2-L^2`, and `lambda=A^2-4*s*H^2`. Angular integration of the cut delta
function yields the following one-dimensional energy integrals:

- If `s>(H+L)^2`, set `E_+-=(omega*A +- q*sqrt(lambda))/(2*s)`.
  Then `16*pi*q*Im B_R = integral_(E_-)^(E_+) de
  [1+n_B(e)+n_B(omega-e)]`.
- If `0<s<(H-L)^2`, use the same endpoints. Then
  `16*pi*q*Im B_R = integral_(E_-)^(E_+) de
  [n_B(e-omega)-n_B(e)]`.
- If `s<0`, set `L_+=(omega*A-q*sqrt(lambda))/(2*s)` and
  `L_-=(-omega*A-q*sqrt(lambda))/(2*s)`. Then
  `16*pi*q*Im B_R = integral_(L_+)^infinity de
  [n_B(e-omega)-n_B(e)] + integral_(L_-)^infinity de
  [n_B(e)-n_B(e+omega)]`.
- The bubble discontinuity vanishes between the two timelike thresholds.

These expressions specialize to H>L, as in this benchmark. All energy
integrals are evaluated analytically using the antiderivative
`T*log(1-exp(-e/T))`. Only the spectator integral remains numerical.

## Numerical convergence and independent checks

The continuum is integrated in the spectator energy E, split at
`H`, `(M^2+H^2-(H+L)^2)/(2*M)`,
`(M^2+H^2-(H-L)^2)/(2*M)`, `(M^2+H^2)/(2*M)`, `M`, and infinity.
A sine-squared coordinate change removes square-root endpoints on finite
intervals. The last interval and the virtual thermal momentum integral use
rational maps to infinity.

| Gauss order | Vacuum Gamma_rho | Thermal Gamma_rho |
|---:|---:|---:|
| 32 | 6.911153274387024e-5 | 1.401135057160701e-4 |
| 64 | 6.911153274387015e-5 | 1.401135042879249e-4 |
| 128 | 6.911153274387078e-5 | 1.401135042879250e-4 |
| 256 | 6.911153274387031e-5 | 1.401135042879244e-4 |
| 512 | 6.911153274387017e-5 | 1.401135042879238e-4 |

The converged contributions are

| Contribution | Vacuum | T=0.5 |
|---|---:|---:|
| HH raised cut | 6.118053358288014e-5 | -8.825872261661233e-5 |
| HHL continuum, all thermal crossings | 7.930999160990028e-6 | 2.283722269045362e-4 |
| Inclusive | 6.911153274387017e-5 | 1.401135042879238e-4 |

The thermal production source is `5.954034823893899e-6`. The negative virtual
contribution is permitted for an individual perturbative correction and
cancels about 39% of the positive continuum contribution.

The vacuum limit agrees with the separate dispersion calculation. An
independent analytic derivative of the thermal real-bubble logarithms also
reproduces the HH result. Quadrature and derivative-step convergence test
numerical accuracy; they do not replace the native implementation comparison.
The saved script and [JSON results](thermal_reference.json) do not incorporate
any native value or fitted normalization.

```sh
python3 examples/cli/thermal_production/massive_self_energy/thermal_reference.py
```
