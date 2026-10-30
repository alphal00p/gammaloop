= Master definitions, data and numerical capabilities

RustRed master data lives in the
#link("https://github.com/alphal00p/vakint-external-data")[external-data repository].
Vakint vendors its 1L–3L programs, coefficients and graphs under `data/`, without
caches. Higher-loop data and prepared caches require an explicit load. See
#link("external_data.typ")[external-data loading] for releases, local checkouts,
snapshot replacement and the `retain`, `symbolic` and `numerical` modes.
FORM reductions and their hard-coded master formulas remain in Vakint.

This document describes the data contract and existing finite coefficients.
It does not prescribe new five-loop numerical calculations or campaign generation.
Use the checked-in dependency revisions and lockfile when reproducing results.

== Definitions, relations and values

A terminal is identified by its full family fingerprint and integer powers.
A friendly name or displayed hash alone is not its mathematical definition.
The family specifies ordered denominators, loop momenta, dimension, scalar-product
coordinates and power shifts. Zero powers remove a denominator; negative powers
put the full denominator polynomial, including its mass term, in the numerator.

Inspect three layers independently:

+ Native programs, their guards and declared terminal keys.
+ Exact symmetry, factorization, terminal-normalization and conventional-basis maps.
+ Laurent coefficients, their exactness and validated numerical precision.

Multiple keys may share one value or project onto a smaller conventional basis.
Neither 74 four-loop output keys nor 905 five-loop campaign terminals is an
independence claim. Per-loop manifests bind each key to its graph, coefficient
metadata and hashed payloads. `symbolic.json` may contain historical placeholder
expressions; the manifest's per-coefficient `symbolic` flag determines whether
an evaluatable exact formula is available.

== One through three loops

- The one-loop K1 terminal `[1]` is the tadpole
  `-Gam(1,1)/(ep*(ep-1))`. K counts denominator coordinates, not loops.
- The two-loop K3 terminals are `[1,1,1] -> -miT111` and
  `[0,1,1] -> (Gam(1,1)/(ep*(ep-1)))^2`.
- The three-loop K6 family has 38 terminal keys. Its exact projections use
  tadpole products, sunset-times-tadpole, `miBN`, `miD6`, and the five-line
  expression involving `miD5`. Retain that complete rational mapping and its
  outer minus sign; the five-line integral is not simply `miD5`.

The K6 denominator ordering is
`k1^2-1, k2^2-1, k3^2-1, (k3-k1)^2-1, (k1-k2)^2-1, (k2-k3)^2-1`.
Its terminals occupy six sector orbits and five distinct mapped expressions.
Definitions and historical mapping implementations are in
#link("src/rustred_evaluation/artifact.rs")[artifact.rs],
#link("src/rustred_evaluation/terminal/k6.rs")[terminal/k6.rs] and
#link("form_src/matad/matad-ng.hh")[matad-ng.hh]. Runtime RustRed evaluation
loads coefficients from the active data registry, independently of FORM tables.
Native program files are `<loop>L/ibps/rules.rrbin` in the data repository.

The retained MATAD authority includes `MASTERS_EXPANSION` and
`MASTERS_SUBSTITUTION` in #link("src/matad.rs")[matad.rs], and
`DIRECT_SUBSTITUTIONS` in #link("src/matad_numerics.rs")[matad_numerics.rs].
The broader legacy basis also includes `miD4`, `miD3`, `miDN`, `miDM`, `miE3`
and `miBN1`. Explicit table tails begin at epsilon^5 for `miT111`, `miBN`
and `miBN1`; epsilon^3 for the `miD*` entries; and epsilon^4 for `miE3`.
These source tails must be propagated through rational prefactors before
assigning a finite composite coefficient range.

The current packaged finite maxima are:

#table(
  columns: (1fr, 3fr, 1fr),
  [Loops], [Master class], [Last power],
  [1], [Tadpole], [`epsilon^4`],
  [2], [Tadpole product], [`epsilon^3`],
  [2], [Sunset], [`epsilon^4`],
  [3], [Three tadpoles], [`epsilon^2`],
  [3], [Sunset times tadpole; five-line integral], [`epsilon^3`],
  [3], [`miBN`], [`epsilon^4`],
  [3], [`miD6`], [`epsilon^2`],
)

The migration preserved all original coefficients and added 25 previously omitted
2L/3L coefficients from existing formulas. The new two-loop numerical coefficient
has 20,000 validated digits; the new three-loop numerical coefficients have 200.
All have exact formulas. Gamma-only closed forms remain recorded, but loading
does not synthesize an unlimited Laurent series beyond the declared finite rows.

== Four-loop programs and finite coefficient coverage

The external repository's
#link("https://github.com/alphal00p/vakint-external-data/tree/main/4L/ibps")[4L/ibps directory]
contains the H, X, BMW and FG parent definitions and saved native data:

- `.toml` and `.csv` describe denominator and momentum ordering.
- `.candidates.rrbin.gz` contain candidate programs and raw terminal declarations.
- `.rrnorm.bin` carry exact weighted terminal normalization.
- `.rrcat.bin.gz` retain exact projections onto FMFT PR expressions.

#table(
  columns: (1fr, 1fr, 1fr),
  [Parent], [Raw terminals], [Retained keys],
  [H], [386], [22],
  [X], [445], [19],
  [BMW], [179], [17],
  [FG], [145], [16],
  [Family-local total], [1,155], [74],
)

Use #link("src/rustred_evaluation/four_loop.rs")[the four-loop loader] and
the native catalog codec to read these files. Native Atom/state serialization
is not a text format. A finite catalog's completeness describes its declared
terminal set; it does not certify arbitrary-index IBP closure.

The conventional basis includes `PR0` through `PR15` and auxiliary/dotted names
`PR4d`, `PR9d`, `PR9x`, `PR11d`. Their definitions and exact relations are in
#link("form_src/fmft/fmft.frm")[fmft.frm]; in particular, `PR9x` is eliminated
by an exact relation. The count “16 masters” does not make every auxiliary
symbol independent.

The complete migration audit inspected all 38 finite FMFT numerical table entries
and all 74 master mappings. Graph factorization certified 22 keys in six classes:
`PR1=tadpole^4`, `PR2=tadpole*miBN`, `PR3=-tadpole^2*miT111`,
`PR6=tadpole*five-line`, `PR8=miT111^2`, `PR13=tadpole*miD6`.
Existing lower-loop formulas establish `PR2` through epsilon^3;
`PR1`, `PR3`, `PR6`, `PR8` through epsilon^2; and `PR13` through epsilon^1.
The `PR1` extension uses the already stored finite order, without introducing
an unbounded Gamma-only expansion policy.

This supplies 24 additional exact and 200-digit numerical coefficients across
20 keys. All 132 existing factorized coefficients were checked numerically at
200 digits, including their already available symbolic formulas. Original
payloads and reproducible validation are retained in the external repository's
`provenance/` and `validation/` directories. Existing consumer `Oep` cutoffs in
FORM tables do not override the independently verified package capabilities:
for example, the stored `PR3ep2` dependency resolves from `miT111` epsilon^4,
and the archived `PR13` epsilon^1 formula is exact despite an older FORM sentinel.

The remaining connected coefficients retain their individually documented
capabilities. Legacy FMFT numerical entries often provide only 26–50 digits;
neither working precision nor a nearby exact coefficient increases those limits.
Per-coefficient manifest metadata, rather than one uniform loop-level range,
is the runtime authority.

== Shared constants and numerical precision

#link("src/hpl.rs")[The shared constant bridge] converts historical `HR`/`HI`
notation to Symbolica's exact sixth-root HPL constants. Words over `{-1,0,1}`
through weight six at `exp(±i*pi/3)` are evaluated at requested precision in
both FORM and RustRed paths. Every one of the 464 removed legacy HPL projections
was checked at 200 and 20,000 digits; the original decimal strings are archived
with the external data. There is no runtime finite HPL substitution table.
Zeta and square-root constants use native exact functions; Clausen and `S2`
use the same HPL evaluator. Finite opaque constants retain honest precision caps.

=== Defining integrals of constants

Every `VakintConstant` record carries `definition()` (Python: `definition`), an
exact Symbolica atom of the multiple integral that defines the constant. It is
documentation, not an evaluation route: the heads below are inert `vakint`
functions without evaluation rules.

#table(
  columns: 2,
  [*Head*], [*Meaning*],
  [`Integral(f, Interval(x1,a1,b1), …, Interval(xn,an,bn))`],
  [Iterated definite integral with the first interval outermost, as in
   Mathematica's `Integrate`. Inner bounds may depend on outer variables;
   bounds may be `∞`.],
  [`LoopIntegral(D, f, k(1), …, k(L))`],
  [`∏_i ∫ d^D k_i / π^(D/2) f` over Euclidean momenta, with `f` rational in
   `EuclideanDot` and analytic continuation in `D` where needed.],
  [`SeriesCoefficient(f, x, x0, n)`],
  [The coefficient of `(x-x0)^n`, not the `n`-th derivative.],
  [`EuclideanDot(p, q)`],
  [Symmetric, linear Euclidean scalar product, distinct from `dot`.],
)

Bound variables are `t`, `x(i)`, `k(i)`, `ep` (`d=4-2*ep`) and `delta`
(`D=2-2*delta`); the only free symbols are `𝜋`, `𝑒` and `γ`. Admitted native
functions are `real_part`, `imaginary_part`, `exp`, `log`, `sin`, `gamma`,
`bessel_k` and `bessel_i`. A definition never references another named
constant: no `vakint` constant, HPL, GPL or master symbol appears inside it,
and `zeta`/`polylog` are admitted only as the complete definition of that
standard function. Every defining integral and every normalization factor is
written out explicitly outside the heads:

- A public FMFT/MATAD master with Minkowski rows `D_M = c_0 + Σ c_j q_j` is
  `(-1)^(Σ a) exp(L γ ep) LoopIntegral(4-2 ep, ∏ D_E^(-a))` with
  `D_E = Σ c_j q_j^E - c_0`. The one-loop tadpole is `-exp(γ ep) Γ(-1+ep)`.
- MATAD internal masters (`miT111`, `miD6`, …) are Euclidean without a sign:
  `exp(L γ ep) LoopIntegral(4-2 ep, …)`. `D6`, `T1ep` and `T1ep2` are their
  coefficients; the external sunset is `-miT111`.
- MATAD bars invert the expansions of `matad.rs`: `D5bar1`, `D5bar2`, `D6bar1`
  and `D6bar2` are `SeriesCoefficient((1-2 ep)(1-ep) exp(3 γ ep) LoopIntegral(4-2 ep, …), ep, 0, n)`,
  `BNbar4` uses `(1-ep)(1-2 ep)(1-3 ep)(2-3 ep)/Γ(1+ep)^3` and `T111bar4` uses
  `√3 (1-2 ep)(ep-1)/Γ(1+ep)^2`, without an `exp(L γ ep)` factor; their values
  are the tabulated HPL combinations. `BNbar4` is a convenience consistent with
  `matad.rs`, not a basis generator: it equals a combination of
  `Li_6(1/2)`, `s6weak` and classical constants.
  The weight-seven bars `D5bar3`, `T111bar5` and `BNbar5` come from external data
  with the same conventions.
- `D=2` periods are `LoopIntegral(2, ∏ 1/(q_E²+1))` and their `D=2-2δ` jets are
  `SeriesCoefficient(exp(L γ delta) LoopIntegral(2-2 delta, …), delta, 0, n)`.
- `K5Moment(k) = Integral(t^k bessel_k(0,t)^5, Interval(t,0,∞))` and
  `K5LogMoment(k) = 3 Integral(t^k bessel_k(0,t)^5 (log t + γ - log 2), …)`.
- `π`, `γ`, `ζ(n)`, `Li_n`, `log(q) = Integral(1/t, Interval(t,1,q))`,
  `Cl_2(π/3)`, `Cl_4(π/3) = √3/12 ∫_0^1 (-log t)^3/(1-t+t^2) dt` (`cl4`) and
  `S2 = 4/(9√3) Cl_2(π/3)` use one-dimensional integrals; `𝑒` and square roots
  are their own elementary definitions.
- `HR`/`HI(a1,…,an, z)` are `real_part`/`imaginary_part` of the
  radial-path HPL at `z=(1±𝑖√3)/2`; each letter contributes `z f_a(z s_i)` over
  `1 > s_1 > … > 0`, trailing zeros become the innermost `log(z s_k)^m/m!`.
  `GPLReal`/`GPLImag` use the kernel `1/(t-a)` on `[0,1]`.
- Recorded period aliases such as `PR4fin` are the stated coefficient of their
  bound master. Historical FMFT and MATAD `NumericalSource` leaves are the complete
  Laurent coefficient they tabulate. Runtime inputs, pySecDec estimates,
  AlphaLoop literals and incomplete table entries (`PR3ep2`, `PR13ep1`) have no
  defining integral: `get_constant` still returns their value, `definition()` is
  `None` and `definition_unavailable_reason()` states why.

Every record also carries `normalisation()`, a plain-language statement of the
conventions of its definition class, shown as "Normalisation" on the constant
card:

#table(
  columns: (1fr, 1.3fr, 2.4fr),
  [*Class*], [*Examples*], [*Normalisation*],
  [Elementary and polylogarithmic numbers],
  [`𝜋`, `γ`, `zeta(n)`, `log(q)`, `polylog(n,q)`, `s6weak`],
  [Plain real number given by the convergent integral shown, first interval
   outermost; no loop measure and no further factor.],
  [Clausen values],
  [`cl2`, `cl4`, `S2`],
  [`Cl_n(π/3) = Im Li_n(exp(i π/3))`, with `Cl_2(π/3) = -∫_0^(π/3) log(2 sin(t/2)) dt`
   and `S2 = 4/(9√3) Cl_2(π/3)`; plain real numbers.],
  [HPL projections],
  [`HR(w, z)`, `HI(w, z)`],
  [Real or imaginary part of `H(w; z)` at `z = (1±i√3)/2` on the straight path
   `t = z s`, `1 > s_1 > … > 0`; letters `1/t`, `1/(1-t)`, `1/(1+t)`; trailing zeros
   become `log(z s)^m/m!` on the principal branch.],
  [GPL and polylogarithm projections],
  [`GPLReal`, `GPLImag`, `PolylogReal`, `PolylogImag`],
  [Real or imaginary part of `G(a_1,…,a_n; 1) = ∫_(1>t_1>…>t_n>0) ∏ dt_i/(t_i-a_i)`,
   resp. `Li_n(z) = z/(n-1)! ∫_0^1 (-log t)^(n-1)/(1-z t) dt`.],
  [Bessel moments],
  [`K5Moment(k)`, `K5LogMoment(k)`, `K5MomentTaylor(k, j)`],
  [`∫_0^∞ t^k K_0(t)^5 dt`; `3 ∫_0^∞ t^k K_0(t)^5 (log t + γ - log 2) dt`; the
   coefficient of `ε^j` of the Taylor expansion (not the derivative) of
   `(e^γ/2)^(3ε) ∫_0^∞ t^(k+3ε) K_ε(t)^5 dt`.],
  [MATAD masters],
  [`D6`, `T1ep`, `T1ep2`, MATAD `NumericalSource` leaves],
  [Euclidean measure `d^d k/π^(d/2)` per loop, `d = 4-2ε`, unit masses, no sign,
   factor `e^(γε)` per loop; the coefficient of `ε^n` of the Laurent expansion,
   not the derivative.],
  [MATAD bars],
  [`D5bar1`, `D5bar2`, `D5bar3`, `D6bar1`, `D6bar2`, `T111bar4`, `T111bar5`,
   `BNbar4`, `BNbar5`],
  [The coefficient of `ε^n` (not the derivative) of the `matad.rs` prefactor
   times `LoopIntegral(4-2ε, …)` (measure `d^d k/π^(d/2)` per loop, no sign):
   `(1-2ε)(1-ε) e^(3γε)` for D5 and D6; `(1-ε)(1-2ε)(1-3ε)(2-3ε) Γ(1+ε)^(-3)`
   for BN and `√3 (1-2ε)(ε-1) Γ(1+ε)^(-2)` for T111, with no `e^(γε)` factor
   because the explicit `e^(-Lγε)` of `matad.rs` cancels the per-loop `e^(γε)`
   of the MATAD master.],
  [Public four-loop masters],
  [`PR4fin`, `PR4dfin`, `PR11dep0`, `PR4ep1`, `PR4dep1`, `PR14ep2`, `PR14ep3`,
   FMFT `NumericalSource` leaves],
  [Minkowski measure `(i π^(d/2) e^(-γε))^(-1)` per loop written as
   `(-1)^(Σa) e^(4γε) LoopIntegral(4-2ε, ∏ D_E^(-a))`, unit masses; the
   coefficient of `ε^n` of the Laurent expansion, not the derivative.],
  [`D = 2` periods],
  [`PR7D2`, `PR10D2`, `PR15D2`],
  [`LoopIntegral(2, ∏ 1/(q_E²+1))`: measure `d^2 k/π` per loop, unit masses,
   positive (hat) sign convention without `(-1)^(Σa)`, no exponential factor.],
  [`D = 2-2δ` jets],
  [`PR7D2Eps1`, …, `PR11dD2Eps4`, `PR0D2Eps2`],
  [The coefficient of `δ^n` of the Taylor expansion (not the derivative) of
   `e^(Lγδ) LoopIntegral(2-2δ, ∏ 1/(q_E²+1))`: measure `d^D k/π^(D/2)` per
   loop, unit masses, positive (hat) sign convention, `e^(4γδ)` for four loops.],
  [No definition],
  [runtime inputs, pySecDec estimates, AlphaLoop literals, `PR3ep2`, `PR13ep1`],
  [`definition` is `None`; the value is the recorded source value and
   `definition_unavailable_reason` states why.],
)

For a definition built from `LoopIntegral` the statement is not looked up by
name: it is derived from the definition itself. It names the coefficient taken,
the measure and dimension, the number of loops, the propagators with their
masses and `Σa`, and every factor in front of the integral (sign, `e^(Lγε)`,
`Γ` factors, polynomials, added terms), and ends with "No other factor". Its
label states the convention those factors realise (public Vakint/FMFT master,
MATAD master, MATAD bar with the matching `matad.rs` prefactor, `D = 2` period
or `D = 2-2δ` jet), or "Loop integral with explicit factors" when none applies;
a statement therefore never contradicts its definition.

The momentum routing transcribes the master's family rows in family order.
External records supply the atom in `symbolic_definition`; records of older
packs are reconstructed from their definition metadata
(#link("src/inspection/definitions.rs")[definitions.rs]).

Numerical mode prefers an available exact formula and evaluates it at the user's
requested precision. A numerical-only coefficient enforces its validated digit
limit. Symbolic mode retains snapshot-qualified symbols for numerical-only
coefficients so later evaluation uses the same data snapshot. These checks
describe source capability, not a certified error estimate for the final result.

== Normalization and required orders

All public `RustRedMaster` symbols, including five-loop residuals and lower-loop
product factors, denote unit-mass integrals in the FMFT/MATAD convention with
`d=4-2*epsilon`. They do not denote bare unnormalized five-loop measures.
The adapter restores the mass dimension and configured loop measure once for
the complete expression; factorization does not multiply that conversion twice.
Power shifts and numerator mass powers remain part of the exact lowering.

The convention uses
`Gam(a,b)=exp(b*epsilon*EulerGamma)*Gamma(a+b*epsilon)`.
`LoopNormalizationFactor::FMFTandMATAD` supplies
`(i*pi^(d/2)*exp(-EulerGamma*epsilon))^(-L)`;
the default pySecDec normalization supplies `(i*pi^(d/2))^(-L)`.
Custom and MSbar settings differ. See the finalizers and `LoopNormalizationFactor`
in #link("src/lib.rs")[lib.rs] for the exact conversion.
Euclidean denominator signs and other programs' tadpole-normalized conventions
must not be treated as interchangeable with these definitions.

For `c(epsilon)*product_i M_i(epsilon)`, let `v_c` be the coefficient valuation
and `ell_j` the leading powers of the other factors. A safe required depth for
factor i through output power T is `K_i >= T-v_c-sum_(j!=i) ell_j`.
The runtime coalesces exact coefficients and cancellations before checking
which powers actually survive. Its output target is
`number_of_terms_in_epsilon_expansion-L-1`; there is no historical fixed
four-loop finite-part ceiling. Missing coefficients produce an error, not zero.

== Five-loop campaign inventory and historical provenance

The initial external 5L pack contains native IBPs, graphs and 905 declared
family-qualified campaign terminals, with no computed five-loop coefficients.
Coverage metadata records the certified scope. Inputs outside it produce a
warning; guarded reduction can retain additional unknown terminals. Final
warnings list unknown keys and parent DOT paths. Numerical evaluation gracefully
reports unresolved keys and declared terminals without coefficients.

The source package was the historical THE_ONE generation 54 publication, scope
`f7ee869008fccc6ed6cf2f4e3342c78d2b8607f59ac4d523dcae504595c1b0f8`:
67 owners, 9,982 rules, 804 recorded requests, 939 source terminal keys and
905 keys after exact finite substitutions. Its original `application.json`
SHA256 was `01425e1cecc71893d2e3a592a162600100326f285914926df61e43dd0c870e4b`.
These are provenance facts, not a statement about the moving D<=10 campaign
or a proof of a minimal basis. The source inventory is distinct from 813 test
witnesses and 804 query regions. The package is now distributed through the
external repository; obsolete machine-specific campaign paths are not load APIs.

RustRed's `load_saved_application::<15>` and `PreparedSavedApplication` expose
`remaining_terminals`, `terminal_sources`, `output_family` and family definitions.
The readable TIDE geometry has 15 shifted coordinates `q_i^2-1`, including the
intentional `D14=(k1+k2-k4)^2-1`. Auxiliary coordinates are not a contiguous
trailing block, and 15 coordinates do not mean 15 physical graph propagators.
The runtime applies exact terminal substitutions once. Native product services
can return lower-family factors; these are resolved by their own family registry,
not mistaken for missing rows in the flat five-loop inventory.

The 13 October 2026 handoff recorded 119 passing native unit tests and four
independently reproduced integration failures at baseline `6203c6cb`. Those are
historical observations, not the acceptance status of this migration. Current
validation must identify the tested revision, data release and scope. See
#link("five_loop_application.typ")[five-loop application details] and
#link("tests/experimental_rustred_4l/README.md")[four-loop acceptance history].
