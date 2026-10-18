= Independent spectral reference for the scalar insertion

The scripts compute the labelled renormalized insertion without using native outputs.
#link("https://arxiv.org/pdf/hep-ph/9409388")[Bauberger et al., Eq. (3) and
Eqs. (13)–(15)] supplies the vacuum mass-derivative/dispersion construction.
#link("https://arxiv.org/pdf/hep-ph/0302098")[Nishikawa, Morimatsu and Hidaka,
section III] supplies the unequal-mass thermal bubble: Eqs. (III.30)–(III.34) are the
vacuum bubble, (III.35)–(III.42) its thermal discontinuities, and (III.44) the thermal
real part. Their 11-component convention is related to the physical discontinuity by
Eq. (II.14); the formulas below use the retarded positive Euclidean bubble directly.

Let `D_a(Q)=1/(Q_E^2+a)` and `I(P)=sum-integral_Q D_a(Q)^2 D_b(P-Q) B(Q)`.
Here `a=b=H^2`, `P=(M,0)`, `M=1.6`, `H=0.6`, `L=0.25`, `T=0.5`, and
`mu_r=1.6`. With unit cubic couplings and graph factor one, `Pi_R=-I_R`, hence
`Gamma_rho=Im(I_R)/M`. This is a graph contribution, not the complete theory at this order.

== Raised two-body pole

The renormalized vacuum bubble is

```text
B_0(s)=-1/(16*pi^2) integral_0^1 dx
 log([x*H^2+(1-x)*L^2-x*(1-x)*s-i0]/mu_r^2).
E_a=(M^2+a-b)/(2*M), E_b=M-E_a,
q_a=sqrt(E_a^2-a), phi(a,b)=q_a/(8*pi*M).
Im I_HH=-d/da { phi(a,b)*[1+n_B(E_a)+n_B(E_b)]
                *Re B_R(E_a,q_a) } at a=b=H^2.
```

Only the repeated propagator's squared mass is differentiated; bubble and spectator
masses remain fixed. This follows `D_a^2=-dD_a/da` in Euclidean conventions.
The repeated shell lies in the bubble's spectral gap. The script therefore uses
analytic continuation and complex mass steps `1e-4`, `1e-8`, and `1e-20`;
the last two agree to floating-point precision. At equal outer masses, derivatives
of their two occupations cancel, but the internal thermal bubble still depends on
the external shell energy and momentum and must be differentiated.

For `s=omega^2-q^2`, `E_i=sqrt(k^2+m_i^2)`, and `A_i=m_j^2-m_i^2-s`,

```text
B_T(omega,q)=sum_(i,j) 1/(16*pi^2*q) integral_0^infinity dk
 k*n_B(E_i)/E_i * sum_(sigma=+-1)
 log |[A_i+2*sigma*omega*E_i+2*q*k]/[A_i+2*sigma*omega*E_i-2*q*k]|,
(i,j)=(H,L),(L,H).
```

There is no principal-value pole near the stable repeated shell.

== Continuum and thermal crossings

Set `E=sqrt(q^2+b)`, `omega_-=M-E`, `omega_+=M+E`, and `s_+-=omega_+-^2-q^2`.
The bubble continuum divided by `(s-a)^2` convolves with the spectator:

```text
Im I_HHL=integral_0^infinity dq q^2/(4*pi^2*E) * {
 [1+n_B(omega_-)+n_B(E)]*Im B_R(omega_-,q)/(s_--a)^2
 +[n_B(E)-n_B(omega_+)]*Im B_R(omega_+,q)/(s_+-a)^2 }.
```

Use `n_B(-omega)=-1-n_B(omega)` and `Im B_R(-omega,q)=-Im B_R(omega,q)`.
Both crossed terms are required. The denominators remain separated from zero on
bubble spectral support for these stable masses.

For positive `omega`, define `s=omega^2-q^2`, `A=s+H^2-L^2`,
`lambda=A^2-4*s*H^2`, and `E_+-=(omega*A+-q*sqrt(lambda))/(2*s)`.
Angular integration gives the following `16*pi*q*Im B_R`:

- Above `(H+L)^2`: integrate `1+n_B(e)+n_B(omega-e)` from `E_-` to `E_+`.
- Between zero and `(H-L)^2`: integrate `n_B(e-omega)-n_B(e)` over the same endpoints.
- For spacelike `s`: sum integrals from
  `L_+=(omega*A-q*sqrt(lambda))/(2*s)` to infinity of `n_B(e-omega)-n_B(e)`,
  and from `L_-=(-omega*A-q*sqrt(lambda))/(2*s)` to infinity of `n_B(e)-n_B(e+omega)`.
- Between the two timelike thresholds: the discontinuity vanishes.

These expressions specialize to `H>L`. The antiderivative
`T*log(1-exp(-e/T))` performs every energy integral analytically.
Only the spectator energy integral remains numerical, split at
`H`, `(M^2+H^2-(H+L)^2)/(2*M)`, `(M^2+H^2-(H-L)^2)/(2*M)`,
`(M^2+H^2)/(2*M)`, `M`, and infinity. Sine-squared maps resolve finite endpoints;
rational maps cover the tails.

Gauss orders 64 through 512 yield the thermal kernel
`1.40113504287924e-4` to approximately 14 significant figures. A separate vacuum
dispersion calculation gives `6.91115327438702e-5`, and an analytic thermal derivative
checks the raised `HH` term. Quadrature convergence establishes the oracle's numerical
accuracy; it does not replace comparison with the native implementation.
