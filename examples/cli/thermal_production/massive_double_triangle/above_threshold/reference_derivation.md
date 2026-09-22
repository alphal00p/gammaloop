# Independent reference above the thresholds

This evaluates the arbitrary-mass scalar double-triangle representation in [Blaizot and Reinosa, hep-ph/0406109, Eq. (58)](https://arxiv.org/pdf/hep-ph/0406109). The numbers are our numerical evaluation at the specified masses, not a published numerical table. No GammaLoop expression enters the reference. The vacuum decomposition into two two-particle cuts and two three-particle cuts also follows [Bauberger and Böhm, hep-ph/9501201, Eqs. (4)–(12)](https://arxiv.org/pdf/hep-ph/9501201).

The routing and normalization agree with the parent example: p1=p, p2=q, p3=p−K, p4=q−K, p5=p−q; masses m1=m3=m5=H=0.6 and m2=m4=L=0.25; external K=(M,0), now M=1.6; T=0.5; zero chemical potentials; four unit cubic couplings; labelled graph factor κ=1/2. For the positive Euclidean master I, Π_R=−κ I_R and Γρ=κ Im I_R/M. The production source is S=n_B(M)Γρ, with an additional (2π)^−3 for the rate per d³p.

| Contribution | Vacuum Γρ | T=0.5 Γρ |
|---|---:|---:|
| LL virtual cut B | +0.00005359169825801313 | +0.00007986779004617656 |
| HH virtual cut A | −0.00001273066440986632 | −0.00007823887809088843 |
| Sum of virtual cuts | +0.00004086103384814681 | +0.000001628911955288129 |
| Sum of real cuts | +0.000002423421102578654 | +0.000003653175781627507 |
| Inclusive | **0.00004328445495072546** | **0.000005282087736915636** |

Here both aggregate virtual and aggregate real contributions are positive. Large cancellations occur between the two virtual cuts and between signed real channels. The below-threshold parent example has opposite signs for its aggregate virtual and real contributions. Individual terms of this interference graph need not be positive. The thermal production source is S=0.00000022445898051219618.

## Newly open cuts and prescriptions

Use e_i=s_i E_i, n_i=n_B(e_i), z=M+i0 and

- A=e3−e1+z,
- B=e4−e2+z,
- C=e5+e2−e1,
- D=e5+e4−e3=C+B−A.

At M=1.6 the HH two-particle cut A opens because M>2H=1.2. The previously open LL cut B remains. Both U=C−A and V=C+B now also have all-outgoing three-particle cuts because M>2H+L=1.45, in addition to their three crossed thermal cuts. Reflection gives equal integrated U and V contributions.

After taking the retarded discontinuity, the two virtual cuts multiply the **real part** of the respective triangle. Their remaining threshold poles therefore carry Cauchy principal values. There is no imaginary double-delta contribution from their overlapping poles. The virtual terms are

Im I_B = coth[M/(4T)] sqrt(1−4L²/M²) Re C_HHH(M²,L²,L²)/(16π),

Im I_A = coth[M/(4T)] sqrt(1−4H²/M²) Re C_LLH(M²,H²,H²)/(16π).

The triangle includes its vacuum and thermal parts. Both Bose prefactors coincide since each on-shell pair has energy M/2. The reference evaluates each principal value explicitly; it does not use a finite iε extrapolation or discard the threshold pole.

## Vacuum triangle

Write a for the two equal side masses and b for the external cut-particle mass. The central mass is H. Thus (a,b)=(H,L) for B and (L,H) for A. A Feynman-parameter representation is

C_vac=(16π²)^−1 ∫dxdy/[a²(x+y)+H²z−M²xy−b²z(x+y)−i0], z=1−x−y.

Set u=x+y and integrate the other parameter analytically. Define A_u=a²u+H²(1−u)−b²u(1−u) and Δ=4A_u−M²u². The real one-dimensional kernel is

- 4 atan(Mu/sqrt(Δ))/(M sqrt(Δ)) for Δ>0;
- 2 log|(Mu−sqrt(−Δ))/(Mu+sqrt(−Δ))|/(M sqrt(−Δ)) for Δ<0.

The script splits at the zero of Δ and squares the endpoint coordinate. It evaluates the vanishing quadratic factor directly in that coordinate, avoiding loss of precision from subtracting nearby floating-point numbers. The resulting vacuum real parts are

Re C_HHH=0.009074680552697728,

Re C_LLH=−0.0030958640063356822.

The triangle imaginary parts are recorded as diagnostics; the two-particle cut rate uses their real parts. The newly open vacuum three-particle contribution is Im I_U=Im I_V=0.000003877473764125846. Together these give the vacuum rate in the table.

## Thermal triangle and explicit principal value

Put each thermal internal line on shell and integrate its relative angle analytically. With q=sqrt(M²/4−b²), p=sqrt(E²−a²), a0=b²+a²−H², define

L_s(E)=log|[a0−sME+2qp]/[a0−sME−2qp]|.

The two equivalent side-line contributions together are

C_T,side=(8π²q)^−1 PV ∫_a^∞ dE n_B(E) sum_(s=±1) L_s(E)/(M²−2MsE).

Only s=+ has a pole, at E0=M/2. For f(E)=n_B(E)L_+(E), d=E0−a, its principal value is evaluated exactly as

PV ∫_a^∞ f(E)/(E−E0) dE

= ∫_0^d [f(E0+t)−f(E0−t)]/t dt + ∫_(E0+d)^∞ f(E)/(E−E0) dE.

The paired numerator removes the pole before numerical integration. A squared coordinate also removes the lower-energy square-root endpoint. When E0<a, as in the M=1 regression, the ordinary unpaired integral is used.

For the central thermal line use r=sqrt(E²−H²), B_E=2qr and c0=b²+H²−a². Its angular integral is

J(E)=log|[(c0+B_E−ME)(c0−B_E+ME)]/[(c0+B_E+ME)(c0−B_E−ME)]|/(2ME B_E),

and C_T,central=(4π²)^−1 ∫_H^∞ dE r n_B(E) J(E). It is regular because every HHL vertex is stable. The resulting thermal triangle corrections are C_T,HHH=−0.00009424040868346835 and C_T,LLH=−0.009538272001851462.

## Real cuts

The crossed-cut energy integrals and signs are derived in the parent reference notes. The new all-outgoing U cut has (s2,s3,s5)=(+,−,+). If F is the light energy and E one heavy energy, its integration domain is

L≤F≤(M²+L²−4H²)/(2M)=0.36953125,

E∈[(M−F−w)/2,(M−F+w)/2],

w=sqrt(F²−L²) sqrt[1−4H²/(M²+L²−2MF)].

Its positive statistical weight is 1+n2+n3+n5+n2n3+n2n5+n3n5, divided by (M²−2ME)(M²−2MF) and multiplied by (64π³)^−1. The vacuum weight is one. A sin² coordinate for F removes both Dalitz endpoints. At T=0.5 its Im I_U is 0.00001542040978726106.

All physical real cuts remain free of uncut-propagator threshold poles. An on-shell uncut line would require a forbidden H↔H+L or L↔H+H process at an internal vertex. More explicitly, the heavy uncut denominator has magnitude at least min[2HL−L²,2HL+L²]=0.2375; the light denominator has magnitude at least min[L²,4H²−L²]=0.0625. The new threshold issue is confined to the virtual triangles.

## Reproduction and convergence

Run with Python and NumPy:

```sh
python3 examples/cli/thermal_production/massive_double_triangle/above_threshold/reference.py
python3 examples/cli/thermal_production/massive_double_triangle/above_threshold/reference.py --mass 1 --output /tmp/double_triangle_below_check.json
```

The script also accepts --temperature. The default output is reference.json beside the script and retains every quadrature order, triangle term, cut contribution, and normalization convention.

| Order | Inclusive thermal Γρ |
|---:|---:|
| 32 | 0.000005282088185841872 |
| 64 | 0.000005282087737202180 |
| 128 | 0.000005282087736915817 |
| 256 | 0.000005282087736915668 |
| 512 | 0.000005282087736915636 |

With --mass 1, the analytically angle-integrated thermal triangle reproduces the parent's direct two-dimensional virtual-cut integral: C_T,HHH=0.004777955558031728 and inclusive Γρ=0.00030576519406186564. This checks the new triangle representation against the unchanged below-threshold calculation, in addition to the independent vacuum Feynman-parameter checks.
