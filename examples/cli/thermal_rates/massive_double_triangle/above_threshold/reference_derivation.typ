= Independent double-triangle reference above threshold

The routing and normalization follow #link("../reference_derivation.typ")[the parent
reference]: `H=0.6`, `L=0.25`, `T=0.5`, `kappa=1/2`, now `M=1.6`.
#link("https://arxiv.org/pdf/hep-ph/0406109")[Blaizot and Reinosa, Eq. (58)] supplies
the thermal representation; #link("https://arxiv.org/pdf/hep-ph/9501201")[Bauberger and
Böhm, Eqs. (4)–(12)] supplies the independent vacuum four-cut decomposition.
No generated GammaLoop expressions or fitted normalizations enter the script.

The newly open `HH` cut and existing `LL` cut multiply the real parts of their virtual
triangles. Their remaining poles carry Cauchy principal values, not an imaginary
double-delta contribution:

```text
Im I_B = coth(M/(4*T))*sqrt(1-4*L^2/M^2)*Re C_HHH(M^2,L^2,L^2)/(16*pi).
Im I_A = coth(M/(4*T))*sqrt(1-4*H^2/M^2)*Re C_LLH(M^2,H^2,H^2)/(16*pi).
```

== Vacuum triangles

For equal side masses `a`, cut-particle mass `b`, and central mass `H`,

```text
C_vac = 1/(16*pi^2) integral dxdy /
        [a^2*(x+y)+H^2*z-M^2*x*y-b^2*z*(x+y)-i0], z=1-x-y.
A(u)=a^2*u+H^2*(1-u)-b^2*u*(1-u), Delta=4*A(u)-M^2*u^2.
```

Integrate one Feynman parameter analytically. The real kernel is
`4*atan(M*u/sqrt(Delta))/(M*sqrt(Delta))` for positive `Delta`, and
`2*log(abs((M*u-sqrt(-Delta))/(M*u+sqrt(-Delta))))/(M*sqrt(-Delta))` otherwise.
The remaining integral is split at the discriminant zero with a squared endpoint
coordinate; its vanishing factor is evaluated directly in that coordinate.
The resulting real triangles are `C_HHH=0.009074680552697728` and
`C_LLH=-0.0030958640063356822`. Their imaginary parts are diagnostics only.

== Thermal triangles and principal values

Put each thermal internal line on shell and perform the angular integral. With
`q=sqrt(M^2/4-b^2)`, `p=sqrt(E^2-a^2)`, and `a0=b^2+a^2-H^2`,

```text
L_s(E)=log(abs((a0-s*M*E+2*q*p)/(a0-s*M*E-2*q*p))).
C_T,side = 1/(8*pi^2*q) PV integral_a^infinity dE n_B(E)
           sum_(s=+-1) L_s(E)/(M^2-2*M*s*E).
```

The `s=+` pole lies at `E0=M/2`. Set `f(E)=n_B(E)*L_+(E)` and `d=E0-a`:

```text
PV integral_a^infinity f(E)/(E-E0) dE
 = integral_0^d [f(E0+t)-f(E0-t)]/t dt
 + integral_(E0+d)^infinity f(E)/(E-E0) dE.
```

This pairing removes the pole before quadrature. If `E0<a`, use the ordinary integral,
which also reproduces the parent's `M=1` result. For the central thermal line,

```text
r=sqrt(E^2-H^2), B_E=2*q*r, c0=b^2+H^2-a^2,
J(E)=log(abs(((c0+B_E-M*E)*(c0-B_E+M*E)) /
             ((c0+B_E+M*E)*(c0-B_E-M*E))))/(2*M*E*B_E),
C_T,central=1/(4*pi^2) integral_H^infinity dE r*n_B(E)*J(E).
```

The stable internal vertices make this integral regular. The thermal corrections are
`C_T,HHH=-9.424040868346835e-5` and `C_T,LLH=-0.009538272001851462`.

== Real cuts and convergence

The all-outgoing `U` cut has signs `(+,-,+)`. The light energy spans
`L<=F<=(M^2+L^2-4*H^2)/(2*M)=0.36953125`; the heavy energy lies between
`(M-F+-w)/2`, where `w=sqrt(F^2-L^2)*sqrt(1-4*H^2/(M^2+L^2-2*M*F))`.
Its statistical weight is `1+n2+n3+n5+n2*n3+n2*n5+n3*n5`.
A sine-squared coordinate removes both Dalitz endpoints. The crossed-cut domains
remain those of the parent derivation. Real uncut denominators stay away from zero:
heavy magnitude at least `0.2375`, light magnitude at least `0.0625`.

At `T=0.5`, virtual `LL=7.986779004617656e-5` and
`HH=-7.823887809088843e-5` leave `1.628911955288129e-6`.
Real cuts contribute `3.653175781627507e-6`, so the inclusive kernel is
`5.282087736915636e-6`. The vacuum kernel is `4.328445495072546e-5`.
Gauss orders 128 through 512 agree within approximately `2e-19` absolute.
Running the same script with `--mass 1` reproduces the independent below-threshold
kernel `3.0576519406186564e-4`.
