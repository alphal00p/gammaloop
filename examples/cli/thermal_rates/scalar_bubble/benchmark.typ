= Massive scalar bubble benchmarks

The portable model defines three distinct real scalars with
$cal(L)_("int")=-g X B_1 B_2$. Its parameters `MX`, `M1`, `M2`, and `g`
are adjustable without changing the generated graph. The command
`generate xs X > B1 B2 FI^2==2 [{{1}}]` produces the one-loop forward
bubble. Both internal masses are positive, chemical potentials vanish,
and threshold subtraction remains disabled at generation and evaluation.

The observable is $Gamma_rho(P)=-op("Im") Pi_R(P)/E$, with
$P=(E,0,0,p)$ and $E>0$. For spacelike $P$ this is a normalized spectral
kernel, not the lifetime of an on-shell particle. The `equilibrium`,
`vacuum`, and `thermal_excess` processes respectively return $Gamma_rho$,
$Gamma_("vac")$, and their difference. The native one-particle flux
supplies $1/(2E)$; the distinct daughters introduce no factorial.

== Literature and normalization

#link("https://arxiv.org/abs/1503.03534")[Ho and Scherrer,
_Cosmological Particle Decays at Finite Temperature_, Phys. Rev. D 92,
025019 (2015)] give boosted equal-mass scalar and fermionic results.
Their Eqs. (2.2)--(2.4) distinguish spectral relaxation from forward
loss: for the neutral bosonic probe,
$Gamma_rho=(1-exp(-E/T)) Gamma_D$.
Their scalar Eqs. (3.2)--(3.3), or directly (A17)--(A18), give a
timelike control. Their interaction is $g Phi phi^2$; for the present
distinct-species interaction with vertex magnitude $g$, replace their
overall $g^2/(8pi)$ by $g^2/(16pi)$. The thermal enhancement ratio after
the spectral conversion is unchanged. In particular, at rest with equal
daughter masses the spectral factor is $1+2n_B(E/2)$, whereas the
forward-loss factor is $(1+n_B(E/2))^2$.

#link("https://arxiv.org/abs/hep-ph/0302098")[Nishikawa, Morimatsu, and
Hidaka, _Thermal sunset diagram for scalar field theories_, Phys. Rev. D
68, 076002 (2003)] provide the arbitrary-mass bubble used as a sunset
subdiagram. In the arXiv version, Eqs. (III.35)--(III.42) give the thermal
pieces and their energy endpoints. They distinguish decay above
$(m_1+m_2)^2$, a closed interval between the two thresholds, absorption
below $(m_1-m_2)^2$, and spacelike support. These expressions are for the
real-time 11 component; Eq. (II.14) supplies the factor
$tanh(E/(2T))$ when recovering the positive-frequency retarded imaginary
part. The reference below instead integrates the corresponding physical
phase-space occupations directly, making this conversion explicit in
the observable convention.

== Independent phase-space reference

`reference.py` uses no native roots, residues, radial profiles, or generated
expressions. Define $s=E^2-p^2$,
$Lambda=[s-(m_1+m_2)^2][s-(m_1-m_2)^2]$, and

$ B(a,b)=integral_a^b dif e n_B(e)
  =T ln((1-exp(-b/T))/(1-exp(-a/T))). $

The script evaluates this logarithm with `expm1` and `log1p` to retain
small thermal contributions. Its optional Gauss--Legendre calculation
integrates the occupation functions themselves; the spacelike half-line
uses a rational coordinate. Agreement of quadrature orders checks the
reference evaluation, not the native implementation.

For $s>(m_1+m_2)^2$, set
$e_i^plus.minus=[E(s+m_i^2-m_j^2) plus.minus p sqrt(Lambda)]/(2s)$.
The integrated result is

$ Gamma_rho = g^2/(16pi E p) [p sqrt(Lambda)/s
  + B(e_1^-,e_1^+) + B(e_2^-,e_2^+)]. $

At rest this becomes
$g^2 sqrt(Lambda)/(16pi E^3)[1+n_B(e_1)+n_B(e_2)]$.
The vacuum contribution is the first term; `thermal_excess` keeps the
occupation terms directly instead of subtracting two rounded totals.

For $0<s<(m_h-m_l)^2$, the only channel is
$X+B_l -> B_h$. The incoming bath energy has endpoints
$a^plus.minus=[E(m_h^2-m_l^2-s) plus.minus p sqrt(Lambda)]/(2s)$ and

$ Gamma_rho=g^2/(16pi E p)[B(a^-,a^+)-B(a^-+E,a^++E)]. $

At rest, $a=(m_h^2-m_l^2-E^2)/(2E)$ and the result is
$g^2 sqrt(Lambda)/(16pi E^3)[n_B(a)-n_B(a+E)]$.
The rest-frame radial equation
$sqrt(k^2+m_h^2)-sqrt(k^2+m_l^2)-E=0$ has a negative derivative
at its physical root. It therefore tests the orientation of the radial
normalization independently of the positive derivative in decay.

For equal-mass spacelike kinematics, both absorption channels contribute.
With $a=[p sqrt(1+4m^2/(-s))-E]/2$,

$ Gamma_rho=g^2/(8pi E p) B(a,a+E). $

The script also handles unequal masses in this domain by integrating each
incoming species separately. Its energy threshold follows directly from
$abs(cos theta)<=1$ in $P+K_i=K_j$; it is not a radial-root calculation.
All three processes vanish in the closed interval
$(m_1-m_2)^2<s<(m_1+m_2)^2$.

== Raised-propagator controls

At fixed external momentum, temperature, and coupling,
$dif/(dif m_i^2) (q_i^2-m_i^2+i 0)^(-1)
  =(q_i^2-m_i^2+i 0)^(-2)$.
Consequently the daughter-mass-squared derivatives of the scalar reference
are independent integrated oracles for a single raised propagator.
The script returns analytic rest-frame derivatives, including the energy
dependence of the occupations. Neither a derivative of the external mass
nor an extra sign is included. For a fermion benchmark the numerator's
mass derivative must be accounted for separately.

The related cold contour issue is treated directly by
#link("https://arxiv.org/abs/2208.14479")[Gorda, Österman, and Säppi,
_Augmenting the residue theorem with boundary terms in finite-density
calculations_, Phys. Rev. D 106, 105026 (2022)]. Their Eqs. (12)--(14)
show the boundary delta term for a double pole and its higher-order
generalization; Eqs. (50)--(51) extend energy-squared differentiation to
bubbles. In the repository's vacuum-subtracted clockwise contour
convention the massive double-pole check is

$ integral (dif^3 k)/(2pi)^3 [-Theta(mu-E_k)/(4E_k^3)
  -delta(E_k-mu)/(4E_k^2)]
  =-1/(8pi^2) ln((mu+sqrt(mu^2-m^2))/m), quad mu>m. $

This last massive expression follows by differentiating the occupied-ball
integral at fixed chemical potential. It is a useful independent factor
for a cold virtual spectator in a rate; it does not establish support
for overlapping Fermi and LU localization constraints.

== Cold Yukawa comparison

The native regression `tests/tests/test_runs/cold_thermal_rates.rs` uses
the Standard Model $H b overline(b)$ vertex to test spin, Pauli blocking,
and a Fermi surface that intersects the physical cut. This is a separate
fixture from the scalar command card below. It retains the model's Dirac
and color tensors. A test-only two-point vertex reuses their identities
with coupling $-i$, so that
$S(q) (-i bold(1)) S(q)=dif S(q)/(dif m)$ for
$S(q)=i(slash(q)+m)/(q^2-m^2+i 0)$.
The insertion differentiates the full propagator, including its numerator.
Both insertion positions are required for the derivative in the common
fermion mass. The physical parameter `MB` varies independently of the
Yukawa parameter `ymb`.

Ho and Scherrer's scalar-to-fermion rate, Eqs. (3.15)--(3.16) and
(C2)--(C3), fixes the zero-density normalization. The finite-chemical-
potential extension is provided by
#link("https://arxiv.org/abs/2007.01224")[Lundberg and Pasechnik,
_Thermal Field Theory in real-time formalism: concepts and applications
for particle decays_, Eur. Phys. J. A 57, 71 (2021)], in particular
Eqs. (VIII.3) and (VIII.9) of the arXiv version.
The expressions here are derived by integrating the corresponding
two-body phase space and then taking its cold limit. Gorda et al.
establish the distributional boundary prescription described above;
their paper is not the source of this Yukawa decay formula.

For $s=E^2-p^2>4m^2$, define $v=sqrt(1-4m^2/s)$ and
$e^plus.minus=(E plus.minus p v)/2$. With the native vertex
$-i y_b/sqrt(2)$ and $N_c$ colors, the spin/color sum is
$N_c y_b^2(s-4m^2)$. Therefore

$ Gamma_("vac") = (N_c y_b^2 s v^3)/(16pi E), quad
  Gamma_rho = Gamma_("vac") integral_(-1)^1 (dif c)/2
  [1-f_F(e_b-mu)-f_F(e_(overline(b))+mu)]. $

Here $c$ is the rest-frame decay angle and
$e_b=(E+p v c)/2$, $e_(overline(b))=E-e_b$.
At zero temperature the occupation functions become steps. The result
is even in $mu$:

$ Gamma_rho = Gamma_("vac")
  min(1,max(0,(e^+ - abs(mu))/(p v))). $

Inside the partially blocked interval $e^-<abs(mu)<e^+$, set
$a=E/2-abs(mu)$. At fixed $y_b$, external momentum, and chemical
potential, the common-mass derivative is

$ (dif Gamma_rho)/(dif m)
  = -(N_c y_b^2 m)/(8pi E) [3v + (4a)/p]. $

The benchmark chooses
$(M,m,E,p,abs(mu))=(5,3/2,25/4,15/4,5/2)$.
Then $v=4/5$, $(e^-,e^+)=(13/8,37/8)$, and the unblocked fraction is
$17/24$. The Fermi momentum is exactly two. On the simultaneous Fermi
and cut support, the occupied outgoing species has bath-frame energy
$5/2$ and polar cosine $5/12$ relative to the probe momentum.
For $N_c=3$ the independent targets are

#table(
  columns: (1fr, 1fr),
  [Quantity divided by $y_b^2$], [Target],
  [$Gamma_("vac")$], [$48/(125pi)$],
  [$Gamma_rho$], [$34/(125pi)$],
  [$Gamma_rho-Gamma_("vac")$], [$-14/(125pi)$],
  [$dif Gamma_rho / dif m$], [$-69/(250pi)$],
  [Insertion on the outgoing $b$ line, $mu>0$], [$-9/(100pi)$],
  [Insertion on the outgoing $overline(b)$ line, $mu>0$], [$-93/(500pi)$],
)

The last two entries exchange under $mu -> -mu$; they are not equal
halves. They follow by giving the two propagators independent auxiliary
masses, differentiating
$N_c y_b^2[s-(m_1+m_2)^2](e_1^+-mu)/(16pi E p)$ in the partially
blocked region, and then setting $m_1=m_2=m$.
Independent angular quadrature using boosted four-momenta and the Dirac
trace agrees with these integrated rates to $5 times 10^(-16)$ absolute.
Mass central differences through step $10^(-5)$ agree with the analytic
derivatives to $2 times 10^(-12)$ absolute.
The native regression compares the rates and both insertions with these
targets separately for positive and negative chemical potentials.

The distributional derivative has to include every differentiation
stage. A zero-order occupation in a double LU pole produces a Fermi
delta when the radial residue differentiates its coefficient, even when
the original coefficient contains no explicit distribution derivative.
Differentiating an already evaluated step function pointwise cannot
recover this contribution. Graph numerator factors remain intact when
separating such distributional support from the regular coefficient.

== Command card and reference targets

From the repository root:

```sh
./gammaloop -n examples/cli/thermal_rates/scalar_bubble/run.toml run generate rest_decay inspect
python3 examples/cli/thermal_rates/scalar_bubble/reference.py --orders 32 64 128 --output /tmp/thermal-rates-bubble-reference.json
```

Replace `rest_decay` by any case below. Each case sets both model masses
and runtime kinematics/temperature for all three processes. The `inspect`
block uses one loop momentum $(0.6,0,0.8)$; pointwise values include the LU
profile and must not be compared directly with the integrated targets.
The `integrate` block uses a shared scratch workspace and `--restart`;
preserve its output before running another case. For an arbitrary reference
point, supply `--energy`, `--momentum`, `--masses m1 m2`, and `--temperature`.
Omit `--orders` to use only the Python standard library.

The following values have $g=1$ and are analytic reference targets. They
are not a native validation record.

#table(
  columns: (auto, auto, auto, auto, 1fr),
  [Case], [$(E,p)$], [$(m_1,m_2)$], [$T$], [$Gamma_rho$],
  [`rest_decay`], [$(2,0)$], [$(0.4,0.7)$], [$0.6$], [$0.012104501639179106$],
  [`boosted_decay`], [$(2.5,1.5)$], [$(0.4,0.7)$], [$0.6$], [$0.009094016657615578$],
  [`closed`], [$(1,0)$], [$(0.6,0.8)$], [$0.6$], [$0$],
  [`absorption`], [$(1,0)$], [$(0.3,1.8)$], [$0.6$], [$0.006880508783133370$],
  [`spacelike`], [$(1,2)$], [$(0.5,0.5)$], [$1$], [$0.010364773879410214$],
)

The two rest-decay mass-squared derivatives are
$-0.005206172268741891$ and $-0.003670617610479006$.
For absorption they are $0.00036810356346982456$ and
$-0.0035964947813623143$, in the order `M1`, `M2`.
Mass-squared central differences with decreasing steps provide a separate
check of these derivative formulas. Thermal threshold subtraction,
infrared subtraction, and resummation are outside this example.

== Native verification

The final focused run passed all 111 selected tests in 171.028 seconds:
nextest `ff7869bf-c0be-49ae-8e97-6881e41a0363`.
The scalar checks cover all five domains, independent pointwise phase
space, event conservation and weight closure, saved-state reload, and
raised mass-squared derivatives for decay and negative-slope absorption.
The physical Yukawa check covers both signs of the chemical potential,
full/vacuum/subtracted rates, the two distinct mass insertions and their
sum, and events on the simultaneous Fermi and physical-cut support.

Nineteen tagged comparisons retain coarse and fine deterministic
quadratures, Richardson extrapolations, analytic targets, and measured
errors. These quadrature differences are numerical checks, not rigorous
integration-error bounds. The largest relative discrepancies were:

#table(
  columns: (1fr, auto, auto),
  [Comparison], [Measured maximum], [Required tolerance],
  [Scalar spectral rates], [$2.65 times 10^(-5)$], [$10^(-4)$],
  [Raised scalar mass-squared derivatives], [$8.17 times 10^(-15)$], [$10^(-8)$],
  [Cold Yukawa full/vacuum/subtracted rates], [$9.97 times 10^(-8)$], [$10^(-4)$],
  [Individual fermion mass insertions], [$6.53 times 10^(-6)$], [$10^(-4)$],
  [Sum of fermion mass insertions], [$1.98 times 10^(-6)$], [$10^(-4)$],
)

The closed scalar channel is exactly zero. The scalar rate tolerance also
contains a $10^(-12)$ absolute floor. At negative chemical potential the
two Yukawa insertions are $-0.05920565268708278$ and
$-0.02864770278166487$, with relative discrepancies
$2.35 times 10^(-7)$ and $6.53 times 10^(-6)$.
Their sum is $-0.08785335546874765$; the positive-chemical-potential sum
agrees within $1.06 times 10^(-9)$ relative. The small difference from
the displayed $y_b=1$ targets uses the resolved native Yukawa coupling
consistently in every comparison.

The negative-chemical-potential shell is translated in the original loop
basis. Gauss--Legendre angular rules of orders 8 and 16 on intervals
split at the exact Pauli boundaries, with radial midpoint refinement,
resolve this contribution while keeping the test below its existing
timeout. Independent integration of the translated shell profile and
the complete mass-derivative density certified this choice before the
native comparison. No target or acceptance tolerance was relaxed.

#link("validation.json")[The validation record] retains the run identity,
commands, source and log provenance, and all nineteen numerical comparisons.
This is focused validation of the stated configurations, not the full
repository test suite or a general claim about higher Fermi derivatives,
arbitrary simultaneous localization constraints, or singularity treatment.
