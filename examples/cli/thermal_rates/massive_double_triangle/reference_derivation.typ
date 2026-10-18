= Independent double-triangle spectral reference

The scripts evaluate the retarded discontinuity of the arbitrary-mass double-triangle
sum-integral in #link("https://arxiv.org/pdf/hep-ph/0406109")[Blaizot and Reinosa,
Eq. (58)]. The routing is `p1=p`, `p2=q`, `p3=p-K`, `p4=q-K`, `p5=p-q`, with
`m1=m3=m5=H=0.6`, `m2=m4=L=0.25`, and `K=(M,0)`. For the positive Euclidean
master `I`, the graph multiplier `kappa=1/2` gives `Pi_R=-kappa*I_R` and
`Gamma_rho=kappa*Im(I_R)/M`. All cubic couplings are one and chemical potentials vanish.
The formulas are independent of generated CFF expressions and LU residues.

== Cut decomposition

For signed on-shell energies `ei=si*Ei`, let `ni=n_B(ei)` and `z=M+i0`:

```text
A=e3-e1+z, B=e4-e2+z, C=e5+e2-e1, D=C+B-A.
R_B=(n2-n4)[n1/(A*C)-n3/(A*D)+n5/(C*D)].
R_U=[n3*n2+n3*(1+n5)-n2*n5]/(A*B), U=C-A.
R_V=[-n1*n4-n1*(1+n5)+n4*n5]/(A*B), V=C+B.
```

Free spectral densities supply `product(si)/(32*product(Ei))` and spatial measure
`d^3p*d^3q/(2*pi)^6`; every retarded pole uses
`Im[1/(M-a+i0)]=-pi*delta(M-a)`. At `M=1`, the virtual `A` cut is closed and
`B` fixes `Eq=M/2`. Reflection gives equal integrated `U` and `V` contributions.
After summing the two uncut signs and integrating the cut delta over the relative angle,

```text
Im I_U = 1/(64*pi^3) sum_(s2,s3,s5) integral dE dF
         [s2*s3*s5*N]/[(M^2+2*M*s3*E)*(M^2-2*M*s2*F)],
N=n3*n2+n3*(1+n5)-n2*n5,
E5=s5*(M-s2*F+s3*E)>0.
```

Below threshold the surviving signs are `(+,+,+)`, `(+,-,-)`, and `(-,-,+)`.
For `s2=+`, integrate `F` from `(M^2+L^2)/(2*M)` to infinity; writing
`a=M-F`, `q=sqrt(F^2-L^2)`, and `Q2=M^2+L^2-2*M*F`, the lower `E` limit is
`[-s3*a+q*sqrt(1-4*H^2/Q2)]/2`. For `s2=-`, start at
`max(L,(4*H^2-M^2-L^2)/(2*M))`, use `a=M+F` and `Q2=M^2+L^2+2*M*F`,
and integrate `E` between `[a+-q*sqrt(1-4*H^2/Q2)]/2`.
Rational maps compactify the unbounded energies; a squared coordinate resolves the pair endpoint.

== Vacuum and thermal triangle checks

The independent Feynman-parameter vacuum triangle is

```text
C0 = 1/(16*pi^2) integral_(x,y>=0,x+y<=1) dx dy /
     [H^2-M^2*x*y-L^2*z*(x+y)], z=1-x-y.
Im I_vac = sqrt(1-4*L^2/M^2)*C0/(16*pi).
```

At `M=1`, `C0=0.012747019233442522`. Direct integration of the vacuum `B` residue
agrees to approximately `1e-17`. Replacing each occupation's medium excess by
`n(sE)-n(sE)|T=0=s*n(E)` gives `C_T=0.004777955558031722` at `T=0.5`.
The virtual result is obtained by multiplying `C0+C_T` by
`coth(M/(4*T))*sqrt(1-4*L^2/M^2)/(16*pi)`, then `kappa/M`.

The real interference is negative below threshold. Every internal `HHL` vertex is
kinematically stable: `H<H+L`, `L<2*H`. Heavy uncut denominators have magnitude
at least `2*M*H-M^2=0.2`; potentially dangerous light denominators have magnitude
at least `L^2=0.0625`. The below-threshold comparison therefore uses no physical
threshold subtraction or thermal infrared treatment.

The script retains Gauss orders 32, 64, 128, 256, and 512. At `T=0.5`, the inclusive
result progresses from `3.057651932923972e-4` to `3.057651940618655e-4`.
#link("above_threshold/reference_derivation.typ")[The continuation] treats the two
virtual principal values and the newly open three-particle decay explicitly.
