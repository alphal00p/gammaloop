# Independent massive double-triangle reference

The benchmark is the massive specialization of the scalar double-triangle sum-integral in Blaizot and Reinosa, [hep-ph/0406109, Eq. (58)](https://arxiv.org/pdf/hep-ph/0406109). The published formula is an exact integral representation, not a tabulated number for these particular masses. The values below are independent numerical evaluations of its discontinuity. No generated GammaLoop CFF or LU residue is used.

Use external K=(M,0), M=1; line masses m1=m3=m5=H=0.6 and m2=m4=L=0.25; all Bose statistics with zero chemical potentials. The momentum routing is p1=p, p2=q, p3=p−K, p4=q−K, p5=p−q. All four cubic couplings equal one, and the labelled graph multiplier is κ=1/2. For the positive Euclidean master I of Eq. (58), Π_R=−κ I_R and Γρ=−ImΠ_R/M=κ Im I_R/M. Every value is for this graph contribution; it is not a complete scalar theory's order-g^4 correction.

| T | Γρ virtual | Γρ real | Γρ inclusive |
|---|---:|---:|---:|
| 0 | 0.00010980937552959695 | 0 | 0.00010980937552959695 |
| 0.25 | 0.0001534382715382693 | −0.000006961419335046441 | 0.00014647685220322285 |
| 0.5 | 0.0003266901870244345 | −0.000020924992962569012 | 0.0003057651940618655 |

At T=0.5 the real contribution cancels 6.40515% of the virtual one. These are interference corrections; individual real channels need not be positive.

The separate production source for an initially empty X population is S=n_B(M)Γρ. It equals 0.00004785764740945725 at T=0.5 and 0.0000027328714165067995 at T=0.25. The differential rate per d³p is S/(2π)³. The tabulated Γρ values above contain no external Bose factor.

Reproduce both quadratures with Python and NumPy:

```sh
python3 examples/cli/thermal_production/massive_double_triangle/reference.py
python3 examples/cli/thermal_production/massive_double_triangle/reference.py --temperature 0.25 --output examples/cli/thermal_production/massive_double_triangle/reference_t025.json
```

Each output retains orders 32,64,128,256,512. The default output is reference.json beside the script. Both runs include an independent vacuum check.

## Discontinuity before integration

Let e_i=s_i E_i, n_i=n_B(e_i), z=M+i0, A=e3−e1+z, B=e4−e2+z, C=e5+e2−e1, D=e5+e4−e3=C+B−A. Insert the free spectral densities in Eq. (58), giving the weight product(s_i)/(32 product(E_i)) summed over all 32 sign assignments and spatial measure d³p d³q/(2π)^6. Take each simple pole using Im[1/(M−a+i0)]=−πδ(M−a).

The A cut is kinematically closed because M<2H. The virtual B cut fixes E_q=M/2 and s2=+,s4=−. Its residue is

R_B=(n2−n4)[n1/(A C)−n3/(A D)+n5/(C D)].

The three-line U=C−A cut, M=e5+e2−e3, has residue

R_U=[n3 n2+n3(1+n5)−n2 n5]/(A B).

The other three-line V=C+B cut has residue

R_V=[−n1 n4−n1(1+n5)+n4 n5]/(A B).

Reflecting the graph gives equal integrated U and V contributions for these equal pairs of masses. The two uncut signs in U can be summed analytically. Integrating the cut delta over the relative angle and changing to positive energies E=E_p and F=E_q gives

Im I_U = (64π³)^−1 sum_(s2,s3,s5) ∫dE dF [s2 s3 s5 N]/[(M²+2M s3 E)(M²−2M s2 F)],

where N=n3 n2+n3(1+n5)−n2 n5, E5=s5(M−s2 F+s3 E)>0, and the angular condition is imposed exactly. Only (+++), (+−−), and (−−+) survive. At T=0.5 their separate Im I_U values are −0.000038541409621711654, +0.00008943154965106708, and −0.00007181513299192444. Multiply their sum by two for U+V, and by κ/M to get Γρ real.

For s2=+ define F0=(M²+L²)/(2M), Q²=M²+L²−2MF<0, a=M−F, q=sqrt(F²−L²). Integrate F from F0 to infinity and E from [−s3 a+q sqrt(1−4H²/Q²)]/2 to infinity. For s2=− use F from max[L,(4H²−M²−L²)/(2M)] to infinity, Q²=M²+L²+2MF, a=M+F, and E between [a±q sqrt(1−4H²/Q²)]/2. The reference script maps these domains to finite intervals and uses tensor Gauss-Legendre quadrature; a squared coordinate removes the timelike pair-production endpoint.

## Virtual check and convergence

The vacuum triangle is independently integrated over Feynman parameters:

C0=(16π²)^−1 ∫_(x,y≥0,x+y≤1) dxdy/[H²−M²xy−L²z(x+y)], z=1−x−y.

It gives C0=0.012747019233442522 and Im I_vac=sqrt(1−4L²/M²) C0/(16π)=0.0002196187510591939. Directly integrating the zero-temperature B-cut residue gives the same C0 to about 10^−17. The thermal triangle correction is obtained by replacing n(sE)−n(sE)|_(T=0)=s n(E) in R_B; it is 0.004777955558031722 at T=0.5. Multiply C0+C_T by coth(M/(4T)) sqrt(1−4L²/M²)/(16π) for the full virtual Im I.

For T=0.5, quadrature orders 32,64,128,256,512 yield Γρ inclusive:

- 0.0003057651932923972
- 0.00030576519405596476
- 0.0003057651940618624
- 0.00030576519406186494
- 0.0003057651940618655

## Why no threshold subtraction is needed

All HHL vertices are stable: H<H+L and L<2H. The virtual HHH triangle is below all of its two-particle thresholds because M<2H and L<2H. In a real cut the heavy uncut denominator obeys |M²±2ME_H|≥2MH−M²=0.2. For the potentially dangerous light uncut denominator, crossed heavy scattering requires F≥(M²+L²)/(2M), so |M²−2MF|≥L²=0.0625; the other light sign gives M²+2MF>0. Thus the surviving cuts have no physical threshold poles. All bath particles are massive, so no soft or collinear resummation is invoked.
