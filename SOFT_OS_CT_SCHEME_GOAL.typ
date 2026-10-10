#let ink = rgb("#192334")
#let muted = rgb("#667085")
#let accent = rgb("#176B87")
#let accent-soft = rgb("#EAF5F8")
#let rule-color = rgb("#CAD5E0")
#let code-fill = rgb("#F4F7FA")

#let horizontalrule = line(
  length: 100%,
  stroke: 0.6pt + rule-color,
)

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
    ])
    .join()
}

#set document(
  title: "Soft and On-Shell Counterterm Scheme Goal",
  author: "GammaLoop Development Notes",
  keywords: ("GammaLoop", "soft counterterms", "renormalization", "CFF"),
)

#set page(
  paper: "a4",
  margin: (x: 20mm, top: 19mm, bottom: 20mm),
  numbering: "1",
  footer: context [
    #set text(font: "DejaVu Sans", size: 7.5pt, fill: muted)
    #grid(
      columns: (1fr, auto),
      column-gutter: 1em,
      [GammaLoop · Soft counterterm implementation record],
      [#counter(page).display("1")],
    )
  ],
)

#set text(
  font: "Libertinus Serif",
  size: 10pt,
  fill: ink,
  lang: "en",
)
#set par(justify: true, leading: 0.68em, spacing: 0.72em)
#set heading(numbering: none, outlined: true)
#set list(indent: 1.25em, body-indent: 0.55em, spacing: 0.28em)
#set enum(indent: 1.25em, body-indent: 0.55em, spacing: 0.28em)
#set table(
  inset: 6pt,
  stroke: 0.45pt + rule-color,
)

#show heading: set text(font: "DejaVu Sans", fill: ink)
#show heading.where(level: 1): set text(size: 20pt, weight: "bold")
#show heading.where(level: 2): set text(size: 15pt, weight: "bold", fill: accent)
#show heading.where(level: 3): set text(size: 11.5pt, weight: "bold")
#show link: set text(fill: accent)
#show raw.where(block: false): set text(
  font: "DejaVu Sans Mono",
  size: 0.86em,
  fill: rgb("#324A5F"),
)
#show raw.where(block: true): it => block(
  width: 100%,
  breakable: true,
  fill: code-fill,
  stroke: 0.45pt + rule-color,
  radius: 3pt,
  inset: 8pt,
  above: 0.8em,
  below: 0.8em,
)[
  #set text(font: "DejaVu Sans Mono", size: 7.6pt)
  #it
]
#show quote.where(block: true): it => block(
  width: 100%,
  breakable: true,
  fill: accent-soft,
  stroke: (left: 2pt + accent),
  inset: (x: 11pt, y: 7pt),
  above: 0.8em,
  below: 0.8em,
)[#it]
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.where(kind: image): set figure.caption(position: bottom)

#align(center + horizon)[
  #block(width: 86%)[
    #text(font: "DejaVu Sans", size: 9pt, weight: "bold", fill: accent)[
      GAMMALOOP · RENORMALIZATION
    ]
    #v(16mm)
    #text(font: "DejaVu Sans", size: 28pt, weight: "bold", fill: ink)[
      Soft and On-Shell\
      Counterterm Scheme Goal
    ]
    #v(7mm)
    #line(length: 42mm, stroke: 1.4pt + accent)
    #v(7mm)
    #text(font: "DejaVu Sans", size: 13pt, fill: muted)[
      Phase 1 · Local soft counterterms\
      with on-shell counterterms deferred
    ]
    #v(16mm)
    #text(size: 10pt, fill: muted)[
      Complete implementation plan, derivations, investigation record,\
      authoritative corrections, and final acceptance reconciliation
    ]
    #v(24mm)
    #text(font: "DejaVu Sans", size: 8.5pt, fill: muted)[
      Source snapshot · 7 September 2026
    ]
  ]
]

#pagebreak()

#outline(
  title: [Contents],
  depth: 3,
  indent: auto,
)

#pagebreak()

#text(size: 8pt, fill: muted)[#raw("<proposed_plan>")]

= Phase 1 --- Local Soft Counterterms, with On-Shell Counterterms Deferred
<phase-1-local-soft-counterterms-with-on-shell-counterterms-deferred>
== Summary
<summary>
Implement the local soft-refined counterterm operator of
#link("https://arxiv.org/pdf/2203.11038v1")[arXiv:2203.11038v1];,
especially §§4.2--4.4 and Appendix B, in both GammaLoop representations:

$ H_d = U_d + S_(d - 1) - U_d S_(d - 1) \, #h(2em) S_(- 1) = 0 . $

The implementation will:

- Cleanly defer all on-shell counterterm work behind explicit
  `unimplemented!` panics. No rest-frame substitution, $gamma^0$, or OS
  projectors will remain as WIP behavior.
- Complete local soft counterterms in the 4D recursive forest and the 3D
  CFF/residue representation.
- Apply every parent operator to the complete nested child atom,
  preserving the exact forest cross terms when parent and child use
  different expansion points.
- Keep integrated soft and OS counterterms disabled until the local
  implementation and nesting tests are complete.
- Add the paper's Appendix-B topology, a fully derived toy nested
  example, inspect-level structural checks, actual
  `profile ultra-violet` CLI tests, and a physical $q_g arrow.r 0$ test
  using a massive-top bubble in a two-loop
  $gamma^star.op arrow.r d macron(d)$ vertex.

== 1. Mandatory implementation startup
<mandatory-implementation-startup>
Before implementation, perform read-only preflight checks: reread
`CONTRIBUTING.typ`, inspect the current worktree, and verify that no
conflicting active goal exists.

The first repository mutation must write this entire
#raw("<proposed_plan>...</proposed_plan>") block
verbatim, including its tags, to:

`/common/dev/gammaloop/lcnbr/SOFT_OS_CT_SCHEME_GOAL.typ`

The literal next tool action, with no intervening inspection or
implementation action, must call `create_goal` with the objective:

#quote(block: true)[
Implement `/common/dev/gammaloop/lcnbr/SOFT_OS_CT_SCHEME_GOAL.typ`: defer
local OS counterterms and complete and verify the local soft-counterterm
scheme.
]

No implementation work starts before that goal exists.

== 2. Scope, dispatch, and failure policy
<scope-dispatch-and-failure-policy>
=== Connected components and schemes
<connected-components-and-schemes>
- Leave `CTIdentifier` and `ApproximationType` unchanged. They already
  select the approximation policy at the correct level: one connected
  1PI component.
- Classify every connected component first. A disconnected spinney is
  retained only when all of its connected components are selected; it
  has no independent approximation scheme.
- Apply the selected operator to each connected component and multiply
  completed results for disjoint components.
- Use HedgePoset as the reference path. Compare with the legacy path
  only on topologies it already supports.
- Keep the existing inert `softct` setting inert; it receives no
  dispatch or precedence semantics.
- Allow hybrid woods such as a soft-refined child and ordinary-MUV
  parent.
- Reject an IR/soft component if `generate_integrated=true`.
- Reject an IR/soft and `PolePart` mixture in the same wood until an
  explicit integrated policy exists.

=== OS cleanup
<os-cleanup>
In both
#link("/common/dev/gammaloop/lcnbr/crates/gammalooprs/src/uv/approx/local_4d.rs")[local\_4d.rs]
and
#link("/common/dev/gammaloop/lcnbr/crates/gammalooprs/src/uv/approx/local_3d.rs")[local\_3d.rs];:

- Route `ApproximationType::OS` directly to:

```rust
unimplemented!(
    "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
)
```

- Remove incomplete OS transformations, rest-frame graph rewriting,
  $gamma^0$ construction, and $P_plus.minus$ projector logic.
- Retain the enum and configuration syntax for compatibility.
- Move useful WIP comments into the goal document or the adjacent
  deferred-OS explanation rather than discarding them.
- Replace malformed OS graph fixtures with one minimal valid
  component-dispatch test which reaches the intended panic.
- Do not retain a massless-OS shortcut to the soft operator: OS is
  wholly deferred in this phase.

=== Color closure
<color-closure>
For bare or amplitude-level fixtures with open fundamental color
indices, multiply by the color-singlet projector

$ frac(delta_i zws^j, N_c) . $

Existing Lorentz projector machinery may close Lorentz indices
automatically; color closure must be explicit. Tests must fail if an
unsaturated color index remains.

== 3. Operator definitions and implementation contract
<operator-definitions-and-implementation-contract>
For a connected component $delta$ of superficial UV degree $d$, define
the ordinary UV deformation

$ Phi_delta^U \( u \) X = X #h(-0.1667em) ({ u p } \, { u m } \, { k } ; D_(e \, u)) \, $

with

$ D_(e \, u) = k_e^2 - M^2 + 2 u thin k_e #h(-0.1667em) dot.op p_e + u^2 p_e^2 - u^2 \( m_e^2 - M^2 \) . $

Define the soft deformation

$ Phi_delta^S \( s \) X = X #h(-0.1667em) ({ s p } \, { m } \, { k }) . $

Thus:

- $U$ scales component-external four-momenta and physical masses and
  introduces the UV mass $M$ through the existing UV propagator
  expansion.
- $S$ scales only component-external four-momenta.
- $S$ leaves loop momenta and every physical mass unchanged.
- $S$ never introduces $M$.
- A terminal denominator already containing $M$ is not reinterpreted as
  a physical-mass denominator by an outer operation.
- A free physical mass in a numerator remains graded by an outer $U$,
  while a literal $M$ remains fixed.

With

$ J_z^(lt.eq n) f = sum_(j = 0)^n frac(1, j !) partial_z^j f\|_(z = 0) \, #h(2em) J_z^(lt.eq - 1) = 0 \, $

the operators are

$ U_d = "ev"_(u = 1) J_u^(lt.eq d) Phi^U \, #h(2em) S_(d - 1) = "ev"_(s = 1) J_s^(lt.eq d - 1) Phi^S \, $

and

$ H_d = U_d + S_(d - 1) - U_d S_(d - 1) . $

`U_d S_{d-1}` means: construct and combine the complete soft Taylor
polynomial first, then apply a fresh ordinary UV operation to that
completed expression.

For $d = 0$,

$ H_0 = U_0 . $

=== 4D construction
<d-construction>
Refactor only as much as needed to reuse the existing compatible-LMB
selection, raw Taylor projection, and final tensor/gamma normalization:

+ Construct raw $U_d X$.
+ For $d > 0$, construct raw $S_(d - 1) X$.
+ Apply ordinary $U_d$ to the completed soft result, producing
  $U_d S_(d - 1) X$.
+ Form $U_d X + S_(d - 1) X - U_d S_(d - 1) X$.
+ Only then apply the forest subtraction sign, CT marker, and canonical
  simplification.

There must be one completed local atom per component branch, not
separately marked UV, soft, and overlap counterterms.

=== 3D CFF construction
<d-cff-construction>
Apply the same algebra to the full component-local CFF expression:

- Determine surface ownership from the connected component's canonical
  route and topology shore.
- Rewrite only E/H-surface occurrences owned wholly by the active
  component.
- Leave co-graph and mixed-boundary surfaces unchanged.
- Scale both external spatial dependence and external-energy shifts.
- Keep loop energies, loop spatial momenta, physical masses, and
  existing $M$-denominators fixed in the soft branch.
- Perform the replacement before Symbolica atom collapse whenever
  structural equality could otherwise merge distinct topology
  occurrences.
- Preserve structural equality of identical CFF denominators; do not
  manufacture different mathematical cache IDs merely to encode
  ownership.
- Carry an occurrence-provenance side table keyed by residue/forest
  history when equal denominators on different shores need different
  ownership decisions. This plan explicitly authorizes that minimal
  internal metadata if the existing cache provenance cannot be extended.
- Error with component, residue, surface, and route context if ownership
  is ambiguous. Never fall back to global replacement.

For $d > 0$, retain the fixed Laurent endpoint $- 1$. In the dressed CFF
expression, Taylor degree $j$ appears at Laurent power $j - d$, so
$j = 0 \, dots.h \, d - 1$ is exactly the strictly negative part
required for $S_(d - 1)$. Bypass the soft branch for $d = 0$.

Apply $U_d$ to the completed CFF soft branch, combine $U + S - U S$, and
return one combined atom per residue. Preserve the export identity

$ \( upright("forest index") \, upright("node index") \, upright("node key") \, upright("term index") \, upright("residue index") \) $

through normalization.

== 4. Why differently centred nested jets define a valid R-operation
<why-differently-centred-nested-jets-define-a-valid-r-operation>
The implementation must be guided by these algebraic invariants rather
than by assuming that all Taylor operators share a centre.

First,

$ 1 - H_d = \( 1 - U_d \) \( 1 - S_(d - 1) \) $

holds exactly, without assuming that $U$ and $S$ commute. Therefore:

- $S$ removes the unwanted low-order soft behavior.
- $U$ subsequently removes the UV behavior of the soft-completed branch.
- The difference from ordinary UV subtraction is

$ H_d - U_d = \( 1 - U_d \) S_(d - 1) \, $

which is UV-improved by construction.

Do not assume that $H$ is idempotent or that $U S = S U$.

For a nested child $gamma subset Gamma$, let $K_gamma$ and $K_Gamma$
denote whichever local operator, $U$ or $H$, is selected for each
connected component. Define

$ P_Gamma = \( 1 - K_gamma \) I \, #h(2em) Z_Gamma = - K_Gamma P_Gamma . $

The nested forest result is

$ R = \( 1 - K_Gamma \) \( 1 - K_gamma \) I = I - K_gamma I - K_Gamma I + K_Gamma K_gamma I . $

The pure outer terms combine as

$ - K_Gamma I + K_Gamma K_gamma I = - K_Gamma \( 1 - K_gamma \) I \, $

so the outer approximation acts on the child-subtracted integrand, not
on the original integrand independently. Conversely,

$ - K_gamma I + K_Gamma K_gamma I = - \( 1 - K_Gamma \) K_gamma I . $

The implementation and tests must establish

$ A_gamma #h(-0.1667em) [K_Gamma \( 1 - K_gamma \) I] = 0 $

for the child hard limit, and

$ A_(Gamma \/ gamma) #h(-0.1667em) [\( 1 - K_Gamma \) K_gamma I] = 0 $

for the reduced outer hard limit.

A sufficient routing condition is

$ A_gamma K_Gamma = K_Gamma A_gamma . $

Enforce this by using a single canonical component route based on
homogeneous, unimodular loop-basis transformations. Do not use a
path-dependent “first compatible LMB”, and reject an affine shift by an
external momentum when it would alter which variable an outer operator
regards as fixed.

The practical nested rule is:

- Finish and combine the child operator first.
- Treat the child boundary momentum as an internal loop variable when it
  belongs to the parent loop system.
- A parent $S$ scales only the parent's external momentum.
- A parent $U$ scales the parent's external momentum and physical
  masses.
- Both parent operators act on the actual completed child expression,
  including its polynomial dependence on the child boundary momentum.
- Existing terminal $M$-denominators remain terminal; physical masses
  still participate in an outer $U$.

This is mathematically sufficient for nested subdivergence cancellation
under the paper's renormalizable, nonexceptional-momentum hypotheses.
Genuine overlapping, non-nested components remain separate compatible
forests in HedgePoset; incompatible components are never multiplied in
one forest.

The 3D implementation must additionally verify, rather than merely
assume, residue compatibility:

$ sum_rho "Res"_rho \( B_(4 d) I \) = B_(3 d) #h(-0.1667em) (sum_rho "Res"_rho I) \, #h(2em) B in { U \, S \, U S \, H } \, $

including nested compositions. Per-residue equality is required only
after matching the canonical forest/residue history.

== 5. Explicit nested toy derivation
<explicit-nested-toy-derivation>
This derivation must be copied into the goal file and converted into
term-by-term tests. Work in a $1 + 1$-dimensional kinematic slice,

$ ell^2 = ell_0^2 - ell_1^2 \, $

while retaining the original four-dimensional superficial degrees for
the Taylor orders.

Define

$ V_ell = ell^2 - m^2 \, #h(2em) U_ell = ell^2 - M^2 \, #h(2em) Delta = M^2 - m^2 \, $

and

$ r = k - p \, #h(2em) K = k #h(-0.1667em) dot.op p \, #h(2em) P = p #h(-0.1667em) dot.op q . $

Take

$ A \( k \, p \) = frac(k_0 + m, V_k V_(k - p)) \, #h(2em) C \( p \, q \) = frac(p_0, V_p V_(p + q)) \, #h(2em) X = A thin C \, $

with

$ d_gamma = 1 \, #h(2em) d_Gamma = 2 . $

=== Child operator
<child-operator>
The ordinary child projection is

$ U_(gamma \, 1) A = frac(k_0 + m, U_k^2) + frac(2 k_0 K, U_k^3) . $

The child soft projection is

$ S_(gamma \, 0) A = frac(k_0 + m, V_k^2) . $

Applying the ordinary UV projection to that completed soft term gives

$ U_(gamma \, 1) S_(gamma \, 0) A = frac(k_0 + m, U_k^2) . $

Therefore

$ H_gamma A = frac(k_0 + m, V_k^2) + frac(2 k_0 K, U_k^3) . $

This explicitly shows that the soft branch retains $m$, while $M$ occurs
only through an ordinary UV projection.

=== Generic parent operator
<generic-parent-operator>
For

$ Z = J_(upright(p h y s)) / V_(p + q) \, #h(2em) J_U \( u \) = j_0 + u thin j_1 + u^2 j_2 \, $

the parent UV projection is

$ U_(Gamma \, 2) Z =  & j_0 / U_p + (j_1 / U_p - frac(2 P j_0, U_p^2))\
 & + [j_2 / U_p - frac(2 P j_1, U_p^2) + j_0 (frac(4 P^2, U_p^3) - frac(q^2 + Delta, U_p^2))] . $

The parent soft projection is

$ S_(Gamma \, 1) Z = J_(upright(p h y s)) (1 / V_p - frac(2 P, V_p^2)) . $

The overlap is

$ U_(Gamma \, 2) S_(Gamma \, 1) Z =  & j_0 / U_p + (j_1 / U_p - frac(2 P j_0, U_p^2))\
 & + [j_2 / U_p - frac(2 P j_1, U_p^2) - frac(Delta j_0, U_p^2)] . $

Consequently,

$ H_Gamma Z = J_(upright(p h y s)) (1 / V_p - frac(2 P, V_p^2)) + j_0 (frac(4 P^2, U_p^3) - q^2 / U_p^2) . $

Define

$ L_m = 1 / V_p^2 - frac(2 P, V_p^3) \, #h(2em) Q_U = frac(4 P^2, U_p^4) - q^2 / U_p^3 . $

=== Parent acting on the bare child
<parent-acting-on-the-bare-child>
For the bare expression,

$ A_0 = frac(k_0, U_k U_r) \, $

and

$ j_0^X = frac(p_0 k_0, U_k U_r U_p) \, #h(2em) j_1^X = frac(p_0 m, U_k U_r U_p) \, $

$ j_2^X = - Delta j_0^X (1 / U_k + 1 / U_r + 1 / U_p) . $

Thus

$ H_Gamma X = p_0 A thin L_m + p_0 A_0 thin Q_U . $

=== Parent acting on the completed child
<parent-acting-on-the-completed-child>
Let

$ hat(A) = H_gamma A . $

Under the parent UV grading,

$ g_0 = k_0 / U_k^2 + frac(2 k_0 K, U_k^3) \, #h(2em) g_1 = m / U_k^2 \, #h(2em) g_2 = - frac(2 k_0 Delta, U_k^3) . $

Then

$ H_Gamma \( hat(A) C \) = p_0 hat(A) thin L_m + p_0 g_0 thin Q_U . $

=== Four distinct forest terms
<four-distinct-forest-terms>
GammaLoop must retain these four signed forest terms separately:

$ F_nothing = + frac(p_0 A, V_p V_(p + q)) \, $

$ F_gamma = - frac(p_0 hat(A), V_p V_(p + q)) \, $

$ F_Gamma = - p_0 A thin L_m - p_0 A_0 thin Q_U \, $

$ F_(gamma ; Gamma) = + p_0 hat(A) thin L_m + p_0 g_0 thin Q_U . $

Only after each term has been matched to its forest provenance may tests
regroup them:

$ F_nothing + F_gamma = frac(p_0 \( A - hat(A) \), V_p V_(p + q)) \, $

$ F_Gamma + F_(gamma ; Gamma) = - p_0 \( A - hat(A) \) L_m + p_0 \( g_0 - A_0 \) Q_U \, $

and therefore

$ R = p_0 \( A - hat(A) \) [frac(1, V_p V_(p + q)) - L_m] + p_0 \( g_0 - A_0 \) Q_U . $

The required asymptotic checks are

$ A - hat(A) = O \( k^(- 5) \) \, #h(2em) A_0 - g_0 = O \( k^(- 5) \) \, $

the pure parent hard behavior is $O \( p^(- 5) \)$, the simultaneous
two-loop hard behavior is $O \( Lambda^(- 9) \)$, and

$ frac(1, V_p V_(p + q)) - L_m = O \( q^2 \) \, #h(2em) Q_U = O \( q^2 \) . $

Run this derivation at generic $m eq.not 0$ to prove mass preservation
and at $m = 0$ to provide a simple soft-profile oracle.

Paper formulas may combine several forest contributions. GammaLoop tests
must first compare each forest family separately; in 3D, residues
belonging to the same forest family may be summed, but raw terms from
different families may not be pre-combined merely to match the paper's
presentation.

== 6. Tests and fixtures
<tests-and-fixtures>
=== Algebraic operator tests
<algebraic-operator-tests>
Add active normalized-atom tests for the soft portions of Eqs. 4.15,
4.16, and 4.19:

- $d = 0$: $H_0 = U_0$.
- $d = 1$: $S_0$ retains the constant soft term.
- $d = 2$: $S_1$ retains constant and linear soft terms.
- Physical masses survive in $S$.
- $M$ appears only in $U$ and $U S$.
- No auxiliary Taylor parameter survives evaluation.
- Verify $1 - H = \( 1 - U \) \( 1 - S \)$ and $H - U = \( 1 - U \) S$.
- Include a noncommuting example demonstrating that $U S eq.not S U$.
- Assert the marker and subtraction sign occur once, after branch
  combination.

=== Nested and disconnected tests
<nested-and-disconnected-tests>
Cover the operator matrix

$ \( U \/ U \) \, quad \( H \/ U \) \, quad \( U \/ H \) \, quad \( H \/ H \) $

for child/parent choices:

- Retain all four signed forest terms.
- Check the child remainder in the child-hard limit.
- Check the paired outer terms in the child-hard limit.
- Check the complementary pair in the reduced-parent limit.
- Check simultaneous hard rays.
- Claim soft improvement only for a complete soft wood, not for an
  arbitrarily incomplete subset.
- Verify all branches use the same canonical route.
- Add an affine-routing negative test which produces a contextual error.
- Test disconnected soft components with HedgePoset and verify
  factorization into completed component results.

=== 3D micro-oracle and ownership tests
<d-micro-oracle-and-ownership-tests>
For a one-time/one-space CFF surface, use

$ E_m \( ell \) = sqrt(ell^2 + m^2) \, $

$ eta \( s \) = E_m \( ell \) + E_m \( ell - s r_1 \) + s r_0 . $

Its soft coefficients are

$ s_0 = 2 E_m \, #h(2em) s_1 = r_0 - frac(ell r_1, E_m) \, #h(2em) s_2 = frac(m^2 r_1^2, 2 E_m^3) . $

For the ordinary UV deformation, with $E_U = sqrt(ell^2 + M^2)$,

$ u_0 = 2 E_U \, #h(2em) u_1 = r_0 - frac(ell r_1, E_U) \, $

$ u_2 = frac(M^2 r_1^2, 2 E_U^3) - Delta / E_U . $

For

$ R \( z \) = frac(N_0 + z N_1 + z^2 N_2, S_0 + z S_1 + z^2 S_2) \, $

use the independent oracle

$ R_0 = N_0 / S_0 \, $

$ R_1 = N_1 / S_0 - frac(N_0 S_1, S_0^2) \, $

$ R_2 = N_2 / S_0 - frac(N_1 S_1, S_0^2) + N_0 (S_1^2 / S_0^3 - S_2 / S_0^2) . $

Test:

- Equal CFF atoms originating on different shores.
- Nested ownership where the child surface changes but its co-graph does
  not.
- Ambiguous ownership producing an error rather than global rewriting.
- $U$, $S$, $U S$, and $H$ against residues of the corresponding 4D
  expressions.
- Nested 4D-to-3D residue compatibility.
- Exactly one final combined $H$ atom per residue.

=== Paper and topology fixtures
<paper-and-topology-fixtures>
Keep separately named fixtures for:

- The exact massive nested self-energy topology presented in Appendix
  B.1.
- The graph labelled B.1 in the paper's figure material, so it cannot be
  confused with the Appendix equation.
- A double triangle with a massless-quark loop.
- Its soft-IR variant.
- Bare massless fermion and gluon self-energies.
- The DGSE/figure-8-style soft self-energy example.
- Consecutive, nested, and disconnected soft components.

The exact massive Appendix-B.1 graph must parse, validate, identify its
connected components, expose the expected forest provenance through
inspect-level output, and reach the intentional OS panic if configured
with the paper's OS assignment. A separately configured IR
specialization must exercise the active soft implementation.

After canonical tensor/gamma normalization, compare the applicable local
soft expressions with Appendix B.5--B.7, the massless B.10--B.12
specialization, and the nested B.16 result. Do not claim coverage of
integrated expressions B.8--B.9, B.13--B.14, B.17 onward, or OS-only
formulas.

Inspect-level tests must display or serialize:

- Selected connected components and their approximation types.
- Canonical loop routes.
- Parent/child forest keys.
- Separate $U$, $S$, and $U S$ provenance before their final
  combination.
- Physical versus UV-mass denominators.
- Residue and CFF-surface ownership.

No profile test may pass vacuously: selected and resolved limit counts
must both be nonzero.

=== Actual UV-profile CLI acceptance
<actual-uv-profile-cli-acceptance>
Every representative graph-level fixture must include a run-card command
that is parsed and executed through the real CLI path:

```text
profile ultra-violet
```

Constructing `Profile::UltraViolet` directly is insufficient as the only
integration test.

Use summed-orientation mode unless a test explicitly targets orientation
bookkeeping. For example:

```bash
profile ultra-violet \
  -p gamma_star_ddbar_top_bubble \
  -i soft_ir \
  --min-scaling 8 \
  --max-scaling 12 \
  --n-points 25 \
  --seed 1337
```

The command writes a profile output directory; acceptance must parse its
`uv_profile.json`, not infer success from the process exit code.

For each locally subtracted fixture require:

- Total selected limits $> 0$.
- Resolved fits $> 0$.
- Finite fitted exponents.
- `pass_fail(-0.9).failed == 0`.

Run a matched unsubtracted control with the same topology, routing,
seed, and scaling window. It must have selected and resolved limits
$> 0$ and fail at least one expected UV limit.

Lightweight primitive and self-energy CLI profiles belong in the normal
suite. Larger Appendix and vertex profiles may be marked serial/slow,
but they remain mandatory in the project's acceptance CI and may not be
silently ignored.

=== Physical $q_g arrow.r 0$ massive-top-bubble vertex test
<physical-q_gto0-massive-top-bubble-vertex-test>
Add a two-loop vertex fixture named `gamma_star_ddbar_top_bubble`:

- External process: $gamma^star.op arrow.r d macron(d)$, with massless
  external $d \, macron(d)$.
- Use $q_g$ for the internal soft gluon momentum to avoid confusing it
  with the external quark label.
- Dress the exchanged gluon with a closed massive top-quark self-energy.
- The two massless gluon propagators on either side of the top bubble
  carry the same $q_g$.
- Give the corrected gluon a stable selector edge, fixed in the fixture
  as `e6`.
- Close external light-quark color with $delta_i zws^j \/ N_c$.
- Use nonexceptional fixed kinematics below the top-pair threshold, with
  $Q = sqrt(p_(gamma^star.op)^2) = 300 med upright(G e V)$,
  $m_t = 173 med upright(G e V)$, and zero top width.

Generate three otherwise identical local-only integrands:

+ `bare`: no local counterterms.
+ `uv_only`: ordinary $U$ on every selected connected component.
+ `soft_ir`: $H$ on the connected massive-top gluon self-energy and
  ordinary $U$ on the overall vertex and other components.

The soft component has $d = 2$, so

$ H_2 = U_2 + S_1 - U_2 S_1 . $

At fixed top-loop momentum $ell$,

$ \( 1 - S_1 \) Pi_(mu nu) \( q_g \, ell ; m_t \) = O \( q_g^2 \) . $

Consequently, the repeated-gluon insertion changes from

$ I_U \( q_g \) tilde.op q_g^(- 6) $

to

$ I_H \( q_g \) tilde.op q_g^(- 4) \, #h(2em) I_H / I_U tilde.op q_g^2 . $

This removes the spurious two-power enhancement but does not remove the
genuine logarithmic soft behavior of the full massless vertex.

Run `profile ultra-violet` on all three variants. The bare variant must
fail at least one hard limit; both `uv_only` and `soft_ir` must pass,
proving that the soft refinement preserves ordinary UV subtraction.

Then run matched soft profiles:

```bash
profile bulk \
  -p gamma_star_ddbar_top_bubble \
  -i uv_only \
  --select 'top_bubble_vertex S(e6)' \
  --min-scaling -2 \
  --max-scaling -5 \
  --n-points 25 \
  --seed 1337
```

```bash
profile bulk \
  -p gamma_star_ddbar_top_bubble \
  -i soft_ir \
  --select 'top_bubble_vertex S(e6)' \
  --min-scaling -2 \
  --max-scaling -5 \
  --n-points 25 \
  --seed 1337
```

Run the same soft profile for `bare` as a diagnostic control.

Acceptance requires:

- Exactly one requested soft report per variant.
- A resolved, finite fit with $R^2 gt.eq 0.98$.
- Identical route, random direction, seed, points, and fit window
  between variants.
- The UV-only fitted scaling is less than $- 1.0$, with theoretical
  expectation approximately $- 2$.
- The soft-refined fitted scaling is greater than $- 0.5$, with
  theoretical expectation approximately $0$.
- Most importantly,

$ upright("scaling")_H - upright("scaling")_U gt.eq 1.5 \, $

with theoretical target $+ 2$.

GammaLoop's bulk profile includes the three-dimensional soft integration
measure, which is why $q_g^(- 6) arrow.r q_g^(- 4)$ appears
approximately as $- 2 arrow.r 0$. Do not require `soft_ir.all_passed`:
the remaining marginal behavior is physical.

This fixture is also the principal mixed-nesting test: the outer
ordinary $U$ must re-expand the complete child $H$ result without
destroying its $O \( q_g^2 \)$ soft remainder.

== 7. Completion criteria and deferred work
<completion-criteria-and-deferred-work>
Phase 1 is complete only when:

- OS dispatch is cleanly and intentionally unimplemented in both
  representations.
- No partial $gamma^0$, $P_plus.minus$, or rest-frame graph rewrite
  remains active.
- IR/soft plus integrated generation fails before expression generation.
- The 4D and 3D implementations satisfy $H = U + S - U S$ for
  $d = 0 \, 1 \, 2$.
- Physical masses remain in soft branches and $M$ arises only through
  ordinary UV projection.
- Nested $U \/ H$ forests pass separate child, parent, and simultaneous
  asymptotic checks.
- The toy example's four forest terms match independently before
  regrouping.
- The Appendix-B.1 topology is a valid, non-vacuous GammaLoop fixture.
- CFF replacement is component-local and occurrence-safe, with ambiguous
  ownership rejected.
- Lightweight profile tests run in the normal suite and
  disconnected-union tests use HedgePoset.
- Representative fixtures execute the actual `profile ultra-violet` CLI
  command with nonzero selected/resolved limits and matched failing bare
  controls.
- The massive-top-bubble vertex passes ordinary UV profiling in both
  `uv_only` and `soft_ir` modes and demonstrates the expected two-power
  $q_g arrow.r 0$ improvement.
- Formatting, targeted tests, and the applicable full test suite pass
  without weakening unrelated tests.

Phase 2 will integrate the single completed local soft atom under a
scheme-specific policy and test $\[ S \] = 0$ and the integrated
Appendix-B formulas. No soft branch may be dropped or replaced by an
integrated identity before nesting is complete.

Phase 3 will revisit OS counterterms after GammaLoop can derive local
counterterms directly from the expanded 4D representation, where the
nonzero rest-frame expansion point and $gamma^0$ projectors can be
introduced without rewriting the 3D CFF graph in place.

#text(size: 8pt, fill: muted)[#raw("</proposed_plan>")]

== Implementation finding: contracted-child grading in hybrid soft/MUV woods
<implementation-finding-contracted-child-grading-in-hybrid-softmuv-woods>
The first numerical DGSE profile falsified an implementation inference
made while realizing the ordered forest recurrence. Applying the same
simultaneous full-parent Laurent grading to both the homogeneous all-MUV
child channel and the finite soft refinement violates the
contracted-graph semantics of the forest formula.

For a degree-$d_gamma$ soft-refined child, write

$ D_gamma = H_(d_gamma) - U_(d_gamma) = \( 1 - U_(d_gamma) \) S_(d_gamma - 1) . $

Write the completed child as $H_gamma = U_gamma + D_gamma$. The soft jet
decomposes $D_gamma$ into child-boundary polynomials of degree
$j lt.eq d_gamma - 1$, with child-UV-finite coefficients after the
factor $\( 1 - U_(d_gamma) \)$. When such an insertion is viewed in a
containing component of superficial degree $d_Gamma$, the reduced-parent
degree of its order-$j$ term is

$ d_Gamma - d_gamma + j lt.eq d_Gamma - 1 . $

Equation 3.3 of the paper represents the inner counterterm as an
effective vertex of the contracted graph $Gamma \/ gamma$. The
nested-cancellation conditions in Eqs. 3.7--3.8 then scale only the
quotient loops while holding the completed child's internal loops fixed.
The discussion following Eq. 4.5 makes the same requirement for the
cancellation of a nested forest pair. In three dimensions the
corresponding Laurent measure is therefore
$lambda^(3 L_(Gamma \/ gamma))$, not $lambda^(3 L_Gamma)$, on the
$D_gamma$ channel.

Linearity gives the required ordinary-parent action

$ U_Gamma \( H_gamma \* Gamma \/ gamma \) = U_Gamma^(upright(f u l l)) \( U_gamma \* Gamma \/ gamma \) + U_Gamma^(upright(q u o t)) \( D_gamma \* Gamma \/ gamma \) . $

The superscripts distinguish only loop-coordinate and loop-measure
grading. The quotient operation still acts on the complete $D_gamma$
expression: it scales a child boundary carrier when that carrier is a
quotient-loop momentum, grades physical masses and free $M$-numerator
factors under the outer $U$, and preserves terminal-$M$ denominator
provenance. Only the child's internal loop coordinates and their measure
are held fixed. For an $H$ parent the same split is used in its $U$ and
$U S$ branches, while its $S$ branch acts on the exact completed sum.

This is not a replacement of $H_gamma$ by $U_gamma$, nor a dropped
forest term. The implementation carries both channels explicitly and
recombines them into the single signed parent atom. A logarithmic parent
has no quotient projection of $D_gamma$, because every term has degree
at most $- 1$. For a positive-degree parent,
$U_Gamma^(upright(q u o t))$ subtracts any non-negative reduced-parent
coefficient. A degree-zero soft child remains trivial because
$H_0 = U_0$. Focused tests must verify both the reduced-parent and
simultaneous hard limits before the positive-degree hybrid is accepted.

The exact DGSE failure exposed the distinction: with the child loop held
fixed, the old two-loop three-dimensional grading inserted one extra
loop measure and shifted the predicted degree from $- 1$ to $+ 2$. The
all-MUV reference channel already had the correct cancellation and
remains unchanged.

== Implementation finding: free versus terminal UV-mass factors
<implementation-finding-free-versus-terminal-uv-mass-factors>
The sentence in the verbatim plan saying that a literal $M$ remains
fixed is too broad. Section 4.1 of arXiv:2203.11038v1 states after Eq.
(4.4) that an outer Taylor operation must also expand the
$m_(upright(U V))$ factors which an inner Taylor operation introduced in
a numerator. Consequently, a free inner-counterterm numerator factor
$m_(upright(U V))^r$ has outer-$U$ Taylor degree $r$, just like any
other polynomial mass factor.

This is distinct from a terminal UV-vacuum denominator. Such a
denominator retains its terminal provenance and vacuum-denominator shape
under an outer operation: it must not be reinterpreted as a
physical-mass propagator and subjected to a second UV rearrangement. The
implementation preserves this role distinction in both representations:

- physical masses and free numerator mass factors are graded once by an
  ordinary outer $U$;
- a soft $S$ operation leaves all masses unchanged;
- terminal $M$-denominators remain terminal and are graded only through
  the homogeneous scaling of the completed local expression.

The corresponding direct oracle must therefore distinguish a free
inner-$U$ numerator factor from the mass embedded in a terminal
denominator; using one symbol at Taylor degree zero to represent both
roles would contradict the paper.

== Implementation finding: mixed CFF occurrences require a partial image
<implementation-finding-mixed-cff-occurrences-require-a-partial-image>
The verbatim plan's instruction to leave a mixed-boundary surface
unchanged must be understood as forbidding a cache-global replacement or
a deformation of its pure co-graph contribution. It cannot mean leaving
the concrete mixed occurrence algebraically unchanged. Taking the nested
energy kernel

$ I_s = frac(1, \( k_0^2 - A^2 \) \( \( k_0 - s ell_0 \)^2 - B \( s \)^2 \) \( ell_0^2 - C^2 \)) $

at the matched poles $k_0 = A$ and $ell_0 = C$ produces

$ "Res" I_s = frac(1, 4 A C thin \[ A - s C - B \( s \) \] thin \[ A - s C + B \( s \) \]) . $

The standalone co-graph residue factor $1 \/ \( 2 C \)$ is fixed. Inside
each mixed CFF denominator, however, the same $C$ represents the
temporal momentum entering the active child and must therefore acquire
the soft parameter. If

$ B \( s \) = sqrt(\( k_1 - s ell_1 \)^2 + m^2) \, $

then at $s = 0$

$ frac(d, d s) \[ A - s C - B \( s \) \] = - C + frac(k_1 ell_1, B \( 0 \)) \, #h(2em) frac(d, d s) \[ A - s C + B \( s \) \] = - C - frac(k_1 ell_1, B \( 0 \)) . $

Freezing the mixed occurrence would miss the $- C$ term, while globally
scaling its cached denominator would also deform pure co-graph
occurrences and energy prefactors. The occurrence-local `Partial` image
is therefore the required residue map: prepare the component-internal
energy dependence and scale only the non-component-energy plus
external-shift combination which represents the active component's
boundary energy. Pure co-graph occurrences remain unchanged, and equal
algebraic denominators on different shores retain separate occurrence
histories.

== Appendix-B v1 traceability notes
<appendix-b-v1-traceability-notes>
The implementation audit found two internal inconsistencies in the
literal v1 Appendix which must not be mistaken for GammaLoop routing or
operator choices.

First, B.3 prints

$ gamma_2 = { e_1 \, e_2 \, e_4 \, e_5 } \, $

and the surrounding B.5 text describes the quotient as a massive
tadpole. That printed set is not the 1PI fermion cycle in the displayed
topology. The topology itself and the surviving $\( k - p \)^2$ cograph
in B.20 instead select ${ e_1 \, e_2 \, e_3 \, e_4 }$, leaving the
massless paper edge $e_5$. In the checked-in GammaLoop graph the
internal-edge map is

$ e_2^(upright(G L)) mapsto e_2^(upright(p a p e r)) \, quad e_3^(upright(G L)) mapsto e_4^(upright(p a p e r)) \, quad e_4^(upright(G L)) mapsto e_5^(upright(p a p e r)) \, quad e_5^(upright(G L)) mapsto e_3^(upright(p a p e r)) \, quad e_6^(upright(G L)) mapsto e_1^(upright(p a p e r)) . $

GammaLoop therefore correctly uses the 1PI child
${ e_2^(upright(G L)) \, e_3^(upright(G L)) \, e_5^(upright(G L)) \, e_6^(upright(G L)) }$
and leaves $e_4^(upright(G L)) = e_5^(upright(p a p e r))$ in the
cograph. Fixture comments and test names record this explicitly instead
of claiming an unqualified literal B.5 edge assignment.

Second, the first quadratic numerator in the TeX of B.16 has a sign
which, read literally, conflicts with both B.7 and direct second
differentiation of the B.1 integrand. The implementation uses the
direct-expansion sign $- q^2 N_0 \/ U_p^4$, and its test records the v1
discrepancy before comparing the physical and UV products as separate
forest families.

== Corrected implementation finding: the outer forest operator remains the full parent operator
<corrected-implementation-finding-the-outer-forest-operator-remains-the-full-parent-operator>
The earlier "contracted-child grading" note above records an
intermediate diagnosis which is superseded by the explicit DGSE power
count and CFF trace below. Contracted-graph scaling is a useful
#emph[diagnostic limit] in the nested cancellation proof, but it is not
a second production definition of $K_Gamma$. The forest recurrence
remains

$ K_Gamma \( 1 - K_gamma \) I \, $

with the one full parent operator acting linearly on the completed child
remainder. In particular, there is no scheme-dependent split in which
the all-$U$ part is graded with all parent loops while
$H_gamma - U_gamma$ is graded with quotient loops only.

For the DGSE fixture, the massless two-edge child has $d_gamma = 1$,
while the containing graph is logarithmic, $d_Gamma = 0$. At fixed child
loop momentum $k$,

$ D_gamma = H_(gamma \, 1) - U_(gamma \, 1) = S_(gamma \, 0) - U_(gamma \, 1) S_(gamma \, 0) $

is independent of the quotient loop momentum $q$. The reduced cograph
has four propagators and a numerator of rank three in $q$, so

$ D_gamma C_(Gamma \/ gamma) = O \( q^(- 5) \) \, $

or degree $- 1$ after the four-dimensional loop measure. Its quotient
hard projection is therefore zero, exactly as the nested forest
recurrence requires for a logarithmic parent. Consequently the DGSE
$H \/ U$ and $U \/ U$ outer-node atoms agree, even though the child
$H - U$ atom is nonzero. Tests must establish both facts non-vacuously.

== Corrected implementation finding: the soft CFF jet must use factorised residue blocks
<corrected-implementation-finding-the-soft-cff-jet-must-use-factorised-residue-blocks>
The erroneous DGSE three-dimensional scaling was not caused by a missing
outer counterterm. For one contributing orientation the unexpanded full
CFF contains the four mixed causal factors

$ \[ 0 \, 5 \, 6 \] \, quad \[ 4 \, 5 \, 6 \] \, quad \[ 1 \, 5 \, 6 \] \, quad \[ 5 \, 6 \, 8 \] . $

Applying the occurrence-local partial image at the exceptional
child-soft point sends all four denominators to the same child energy
sum $E_5 + E_6$. The result contains $\( E_5 + E_6 \)^(- 4)$, whereas
the residue of the Taylor-expanded four-dimensional graph contains only
one child causal factor and three causal denominators of the contracted
cograph. Losing those three quotient denominators shifts the observed
hard degree by precisely three powers, from the required $- 1$ to $+ 2$.

Thus the earlier mixed-occurrence note is locally correct as a
derivative of one generic, nonexceptional denominator, but it is
insufficient as a definition of the #emph[completed] CFF soft jet:
taking the exceptional soft Taylor jet and resolving the energy residues
is non-uniform because poles coalesce. The production source for a soft
branch must instead be the fixed-orientation factorisation

$ "CFF" #h(-0.1667em) (delta \/ { upright("immediate soft children") }) product_(gamma prec delta) "CFF" #h(-0.1667em) (gamma \/ { upright("its immediate soft children") }) "CFF" \( Gamma \/ delta \) \, $

with the corresponding numerator kernels and with the full-orientation
selector inserted exactly once. A soft map acts on the active component
block and its numerator while the pure cograph block remains literal; an
outer operation subsequently acts on the actual completed product.
Nested and disjoint soft components therefore form a laminar block
forest. Equal denominators retain mathematical equality, while
block/occurrence provenance records where each factor came from.

The same `CutCFFIndex` and canonical full orientation must be retained
while the block inventory is allowed to differ between $U$, $S$, and
$U S$. Branch materialisation must consequently validate a union of
branch-tagged occurrences rather than require byte-for-byte identical
occurrence tables. If a selected residue cannot be transported to the
factorised block forest unambiguously, generation must fail with graph,
component, route, residue, orientation, and block context; it must never
revert to a global surface replacement.

== Corrected implementation finding: Phase-1 hybrid closure is degree dependent
<corrected-implementation-finding-phase-1-hybrid-closure-is-degree-dependent>
The unconditional hybrid-wood allowance in the verbatim plan is too
broad. Let a positive-degree child use (H\_gamma), and isolate its
finite scheme refinement (D\_gamma=H\_gamma-U\_gamma). A boundary Taylor
monomial of order (jd\_gamma-1) has reduced-cograph degree

$ d_(G a m m a \/ gamma) \( D_g a m m a \)_j = d_G a m m a - d_g a m m a + j . $

The worst retained monomial therefore has degree (d\_Gamma-1). A
logarithmic ordinary parent ((d\_Gamma=0)) is safe: its quotient
projection vanishes, as the explicit DGSE power count above shows. For
an ordinary MUV parent with (d\_Gamma\>0), however, the same hybrid
leaves a nonnegative reduced-cograph limit. Two independent 4D fixtures
expose the limiting degree zero for (d\_Gamma=1), even though the
full-parent recurrence is implemented exactly.

Replacing the full parent operation by a quotient-graded operation on
(D\_gamma) is not a valid repair: it contradicts the Appendix-B.16
completed-child expression and makes the production operator depend on
an arbitrary algebraic split. The missing object is a finite
scheme-change counterterm after the child has been integrated to a local
vertex, which is part of the deferred integrated policy. Phase 1 must
consequently allow (H/U) only for (d\_Gamma=0) and reject a
positive-degree ordinary-MUV ancestor of a retained soft refinement with
a contextual error. Same-scheme (H/H) nesting remains active, as does
the logarithmic DGSE (H/U) case.

== Corrected implementation finding: ordinary production and soft replay are separate 3D channels
<corrected-implementation-finding-ordinary-production-and-soft-replay-are-separate-3d-channels>
Ordinary U continues to act on the established full CFF production atom.
Replacing that atom by a product of generated component/cograph blocks
at a generic hard point changes the residue representation and spoils
ordinary UV cancellation. In parallel, each deferred term retains a
laminar factorised replay state for a possible later S operation.
Ordinary U updates that replay with the same Taylor parameter and route
but does not publish it as the production atom. Only S consumes the
replay and reconstructs the production atom from its blocks.

Each replay block therefore carries its own transformed-occurrence
history. An outer S must freeze the actual image left by an inner U, S,
or US, rather than cancel its tag with the original physical
denominator; otherwise a spurious ratio such as (D/U(D)) survives and
the ordered jet is lost. Replay provenance becomes public production
provenance only when S consumes the factorisation. Splitting a
previously deformed container is rejected if its ordered occurrence
images cannot be transported unambiguously.

== Corrected implementation finding: only soft components create CFF block boundaries
<corrected-implementation-finding-only-soft-components-create-cff-block-boundaries>
The fixed-orientation laminar factorisation is an exceptional-point
replay program, not an alternative value for a generic ordinary
counterterm. In particular, for an ordinary child $gamma subset Gamma$
there is no orientation-wise identity

$ B_Gamma^(\( o \)) = B_gamma^(\( o \)) B_(Gamma \/ gamma)^(\( o \)) . $

Consequently an ordinary $U_gamma$ operation must deform the smallest
existing causal container in place and record its unprojected
deformation history; it must not introduce a new component/cograph
split. Only a soft operation creates the block boundary required when
its exceptional limit coalesces poles. The replay root is the maximal
internal causal region: external half-edges are not part of that region
and must not force an overall soft component through a spurious split.

If a later proper soft component splits a container carrying earlier
operations, occurrence images from the discarded causal tree cannot
simply be copied to the regenerated trees. There is generally no
occurrence bijection, and an operation may have changed inverse-energy
or tree factors even when it emitted no surface event. Instead, retain
the exact ordered, unprojected deformation program. For each regenerated
block, replay the inherited stages in chronological order before
appending the new soft stage, while leaving the already-deformed
numerator/measure kernel untouched. A stage represented by an already
contracted child block is not replayed on its cograph.

This ordering follows from the ring-homomorphism property of the
deformation,

$ Phi \( K product_i C_i \) = Phi \( K \) product_i Phi \( C_i \) \, $

whereas the Taylor projector itself is not multiplicative:

$ J \( f g \) eq.not J \( f \) J \( g \) . $

The formal Taylor parameters therefore remain shared across the
regenerated factors and are projected only after their product has been
recombined. The same replay representation is required for inherited
ordinary and soft stages; special-casing only $U \/ H$ would silently
mis-order $H \/ H$ woods. Any non-laminar scope or route which cannot be
replayed exactly remains a contextual generation error rather than
falling back to global replacement.

== Implementation finding: selected CFF residues require local count provenance
<implementation-finding-selected-cff-residues-require-local-count-provenance>
A selected `CutCFFIndex` records a residue order, but it does not
identify the raised surface group which was selected. The factorised
replay therefore retains the originating `CutSet` and assigns a local
right-threshold, left-threshold, and LU count to every causal block.
`Some(0)` is meaningful internally: the selector exists globally, while
that particular factor must contain no occurrence. Branch selection is
performed in the same right-to-left-to-LU order as `Graph::cff`, and a
selected denominator is mapped to the residue unit only after the
matching branch count has been enforced.

Pruning a causal tree must not renumber away occurrence provenance. A
`NormalizedCFFBlock` consequently stores the complete occurrence object
in each tree node, rather than recovering it from the post-pruning node
index. This preserves the original occurrence identity, shore, and
support even when equal denominators occur on different branches.

The currently proven transport is deliberately conservative. An already
selected complete block can be transformed directly. When a block is
split around a component, every regenerated occurrence is first aligned
with the orientation-filtered occurrence table of its source container.
The key is the block-local generated shore, compared with the source
shore after restriction to the visible region and modulo exchange with
its complementary component shore, together with signed energy support.
A source on the same restricted shore is preferred; use the
complementary shore only when no same-shore source matches. Within that
class, when several global shores have the same component restriction,
compare their sets of source-shore vertices outside the component by
inclusion. Discard only strict supersets; retain incomparable minimal
extensions, which therefore remain an ambiguity error unless they give
the same exact copied image. Cardinality is not a valid topological
order for branching cographs. For the freshly generated component block,
the source energy support is projected onto the exact visible
uncontracted domain `paired_edges(region - generator_contract)` before
comparison. This is only a matching key: after the unique match, its
effective surface ID, complete denominator, and unprojected support are
copied from the source while the generated surface ID, local shore, and
stable node ID remain block-local provenance. For example, the
top-bubble block generates support $E_7 + E_8$, while the
fixed-orientation full-CFF source shore carries $E_5 + E_7 + E_8$;
projecting the source key identifies the component shore, but retaining
its full support lets the soft map form $E_7 + E_8 + s E_5$. Cograph
matching remains exact and unprojected, since projecting a cograph
source could copy a contracted child's energies back into the quotient.
After topological ranking, source contexts which give exactly the same
copied image -- equal surface ID, full support, and algebraically equal
denominator -- are one operational remap candidate: source context is
not copied, and the generated tree keeps its own occurrence node and
multiplicity. Two distinct full images at the same best rank remain an
ambiguity error. In particular, regenerated and source node indices are
not used as a tie-break because their causal trees can have different
local numbering. The block atom is then rebuilt before any residue
selection. This shore lift is essential: contraction can erase a literal
external-energy shift such as $plus.minus Q_5^0$, and no later soft
boundary map can reconstruct a symbol which is already absent. Copying
the unique source occurrence leaves a pure cograph denominator literal
while still allowing the active component's soft map to scale its own
copy of the boundary shift. The freshly aligned component/cograph blocks
are the residue representation of the corresponding #emph[soft forest
branch];; they are not an undeformed identity factorisation of the
generic root CFF at (s=1). The top-bubble orientation zero is an
explicit counterexample: its full graph orientation is acyclic and has a
nonzero root-CFF atom, but contracting the top bubble makes the
inherited cograph orientation cyclic. The component block is nonzero and
the contracted cograph block is therefore zero. Requiring their product
to reproduce the generic root term would compare two different forest
graphs and reject the mathematically correct residue representation. The
top regression records this counterexample explicitly, while the exact
operation-before-residue tests compare the (S), (US), and (H) branches
with their corresponding 4D residues. The mandatory summed-orientation
top-bubble profile supplies the full-graph physical acceptance.

After this lift, a selected count may be assigned wholly to the cograph
only if the prospective zero-count component block contains no
occurrence of the selected raised surface. Exact support/denominator
matching is preferred; an energy-support fallback is allowed only when
no global surface outside the selected raised group has the same energy
multiset. A positive local count with no match is an error, as is a zero
local count with a match, because the latter requires residue-count
convolution. All physical cut edges are passed to the regenerated block.
If the selected cut intersects the component, a source shore has zero or
multiple lifts, or a selected-equivalent occurrence appears in a block
provisionally assigned count zero, Phase 1 rejects the operation with
graph, component/block, residue, orientation, cut, and route context
rather than generating an unselected or shift-erased CFF.

== Superseded implementation finding: fixed-orientation provenance versus numerical profiles
<superseded-implementation-finding-fixed-orientation-provenance-versus-numerical-profiles>
An earlier diagnostic proposed selecting one stored orientation inside
the computed-inspect exporter and remapping a dense occurrence index
back to that stored index. That selector was not retained: graph terms
already contain the orientations selected during generation, and no
export consumer selected a further orientation. Computed export
therefore consumes the graph term's existing orientation set without
another filtering or indexing convention. Actual orientation-local
validation selects orientations through the normal generation/profile
path, while summed-orientation UV and bulk profiles remain the
additional physical cancellation checks.

The large top-bubble fixture uses a computed HedgePoset 4D forest export
for its graph-specific parent/child, route, and U/S/US/H provenance
assertions. That export integrates only the recursive 4D forest, emits
one all-none residue atom per node, and deliberately bypasses every
three-dimensional $\( 2 pi \)^(3 L)$, LU, and CFF post-factor. Its
actual generated integrands and all UV and $q_g$-soft CLI profiles
remain three-dimensional and summed over orientations. Smaller exact 3D
residue oracles retain the detailed CFF route, mass-wrapper,
occurrence-ownership, and nested-replay checks. This separation avoids
recomputing and serializing a giant exact 3D provenance atom merely to
repeat structure already checked by those oracles.

Generation RAM telemetry polls at one-second cadence. `sysinfo` 0.33
scans the host `/proc` inventory even for a single requested PID, so the
previous 100-ms cadence consumed a substantial fraction of one CPU
during exact Symbolica work. The monitor parks with a timeout and is
explicitly unparked before joining, preserving immediate generation
shutdown despite the longer sampling interval.

Profile JSON is itself part of the acceptance boundary. The serialized
form therefore retains `allow_vanishing_missing_fits`, aggregate
inspect-fit status, and per-orientation inspect-fit status. `pass_fail`
reads those public serialized entries after a round trip, so parsing
`uv_profile.json` has the same missing-fit and vanishing-limit semantics
as the in-memory result.

== Latest implementation status: production deforms the summed unsplit CFF
<latest-implementation-status-production-deforms-the-summed-unsplit-cff>
This finding supersedes the later claims above that a soft branch should
use a fixed-orientation component/cograph factorisation, that `S` should
publish such a replay as its production atom, or that blockwise
fixed-orientation identities are acceptance criteria. Those
constructions may remain useful for topology and occurrence-provenance
diagnostics, but they are not the mathematical CFF integrand on which
the production Taylor operator acts.

The physical CFF is the sum over compatible orientations. Contracting a
soft component does not preserve acyclicity orientation by orientation:
an acyclic full-graph orientation can induce a cyclic contracted
cograph, even though its term in the unsplit full CFF is nonzero.
Contributions therefore redistribute and cancel between orientations,
and in general there is no identity

$ "CFF"^(\( o \)) \( Gamma \) = "CFF"^(\( o \)) \( delta \) "CFF"^(\( o \)) \( Gamma \/ delta \) $

for each fixed orientation $o$. Replacing an unsplit term by that
product changes the residue representation before the Taylor jet is
taken. Production must instead apply the component-local soft
deformation to the complete summed, unsplit CFF atom, retain occurrence
provenance only to distinguish equal owned and unowned factors, and
project the Taylor series after that complete atom has been deformed.
Fixed-orientation block tests must consequently be labelled as
diagnostic construction tests; they must not require root/block
equality, branch cancellation, or operation-before-residue identities
orientation by orientation. Physical residue comparisons and profile
acceptance are summed over the matched orientations.

The massive-top-bubble child profile is the numerical acceptance for
this ordering. On the same route, seed, direction, and fit window,
ordinary $U_2$ scales approximately as $- 2$, while the summed unsplit
$H_2 = U_2 + S_1 - U_2 S_1$ result scales approximately as $0$. The
observed two-power improvement is the required cancellation; differing
individual orientation scalings are not a replacement for this summed
physical result.

The canonical stored spatial-vector form is
`Q3(edge, minkowski_representation)`, constructed as
`GS.emr_vec(edge, Minkowski {}.new_rep(4).to_symbolic([]))`. An
explicitly indexed component may use `Q3(edge, cind(mu))`, but a typed
bare `Q3(edge)` is not a second convention. Replacement and routing code
must match `Q3(edge, W_.x___)` and preserve its second argument so
represented tensors and explicit components follow one deformation
without losing their structure.

== Authoritative correction: exceptional residues require conditional causal factors
<authoritative-correction-exceptional-residues-require-conditional-causal-factors>
This section supersedes the immediately preceding “Latest implementation
status” section. It does not alter the verbatim implementation plan at
the start of this file.

The summed, unsplit CFF remains the correct production representation
for an ordinary $U$ operation at generic kinematics. It is not, however,
the correct object on which to take the component-local exceptional soft
jet after the residue sum has flattened distinct child and co-graph
poles. In the massive-top-bubble test, three mixed child/co-graph causal
denominators coalesce with the child denominator at the child soft
point. Taylor expanding that collapsed atom produces a fourth power of
the child denominator and loses the three quotient denominators. This is
precisely the observed spurious $+ 2$ scaling in five reduced-parent
hard rays. Occurrence tags can distinguish equal factors, but cannot
recover causal factors which have already coalesced algebraically.

The correct exceptional residue is a #emph[conditional] causal
factorisation; it is not the rejected construction in which one original
full-graph orientation is imposed on every factor. For each quotient
orientation $o_Q$, generate the contracted quotient first. Copy its
directed boundary edges into the split edges of the child, independently
sum every acyclic child-internal orientation compatible with that
boundary, apply the child operation to that completed child sum, and
only then multiply the untouched quotient factor:

$ "Res" \( S_delta I \) = sum_(o_Q) Theta_Q \( o_Q \) thin C_Q \( o_Q \) S_delta #h(-0.1667em) [sum_(o_delta divides partial o_delta = partial o_Q) Theta_delta \( o_delta \) thin C_delta \( o_delta \) thin N_delta] . $

The same ordering defines $U S$: finish the soft Taylor polynomial of
the child factor first, apply the ordinary child $U$ operation to that
completed factor, and multiply the quotient only afterwards. The
ordinary $U$ branch continues to use the established unsplit CFF path. A
containing operator must replay the complete factorwise $S$ or $U S$
branch from retained block and forest history; collapsing it and
globally guessing ownership is invalid.

Selectors are owned exactly once. Original paired edges internal to the
child contribute child selectors; quotient edges, including the child's
boundary carriers, contribute quotient selectors. Split boundary
directions condition child generation but do not contribute a second
selector. A conditionally valid child sector need not lift to any single
globally acyclic root orientation. Such a sector must be retained with
generated factor-local occurrence provenance; source-CFF orientation
matching is used only to align a unique occurrence image when one
exists, never to pair or prune causal factors.

The top-bubble oracle establishes the concrete structure: the one-loop
quotient has $2^4 - 2 = 14$ acyclic orientations, each has two
conditioned child orientations, and every product contains exactly three
quotient causal denominators and one child causal denominator. The
inverse-energy factors partition edges $e 3 \, dots.h \, e 6$ and
$e 7 \, e 8$ exactly once; there are no internal bridge/tree factors.
Reversing the two boundary causal arrows alone does not send $q_g$ to
$- q_g$. The physical reversal combines the matched boundary-energy
sectors with reversal of the canonical spatial carrier
$Q 3 \( e 6 \, upright(m i n k) \( 4 \) \)$. Under that correlated map
the selector-free $S_0$ coefficient is even and $S_1$ is odd.

Selected residues are convolved across factors. For every total
`CutCFFIndex` order, enumerate non-negative child and quotient
allocations whose componentwise sum is the requested right-threshold,
left-threshold, and LU order. Select occurrences inside each factor
before multiplication, retain selected images and factor-path occurrence
namespaces, and reject ambiguous source alignment or overlapping raised
groups contextually.

Finally, the single spatial-momentum convention is a two-argument
function: `Q3(edge,mink(4))` for a represented spatial vector,
`Q3(edge,mink(4,index))` for its automatically indexed tensor-network
image, or `Q3(edge,cind(i))` for an explicit component. Bare `Q3(edge)`,
arbitrary second arguments, and non-four-dimensional Minkowski
representations are invalid. Replacement patterns use the exact one-slot
form `Q3(edge,W_.x_)`; `W_.x___` is not a supported compatibility
convention.

== Authoritative implementation-state addendum: completed conditional stages and one momentum representation
<authoritative-implementation-state-addendum-completed-conditional-stages-and-one-momentum-representation>
This section supplements the “Authoritative correction: exceptional
residues require conditional causal factors” section and has precedence
over every earlier post-plan implementation finding where they disagree.
Historical findings remain in this file for traceability.

The ordinary $U$ branch at generic kinematics continues to use the
established summed, unsplit output CFF. Every branch containing an
exceptional soft operation $S$, including $U S$ and a containing
operation replaying a completed child $S$ or $U S$ branch, uses the
quotient-first conditional causal-factor tree. For each quotient
orientation, first sum all acyclic child orientations compatible with
its directed boundary, apply the connected child operation to that
complete child sum, and only then multiply the untouched quotient
factor. No single full-graph orientation is imposed on both factors.

The output `CutCFF` determines the emitted residue-key domain. The
separate unpruned, default-orientation alignment `CutCFF` is used only
for occurrence provenance and raised-surface residue classification. A
generated conditional factor and its generated denominator are always
the mathematical production object; an aligned source occurrence never
replaces that denominator, supplies a causal branch, or prunes a
conditionally valid child sector. Occurrence identity includes the
output residue, local residue allocation, factor path, and generated
local occurrence namespace. Equal mathematical denominators remain
equal. Requested right-threshold, left-threshold, and LU orders are
convolved additively over child and quotient factors.

Every connected deformation stage is a completed sequential operator. If
$A_0$ is the complete active-factor expression, then

$ A_i = "ev"_(z_i = 1) J_(z_i)^(lt.eq n_i) Phi_i \( z_i \) A_(i - 1) \, $

where $n_i = d_i$ for $U$ and $n_i = d_i - 1$ for $S$. The parameter
$z_i$ is shared by the numerator, inverse-energy factors, causal
denominators, tree factors, and every compatible-orientation term
belonging to that active factor. The Taylor projector therefore acts on
their complete sum and product, not on individual denominators or
orientations.

After that projection, $A_i$ contains no live $z_i$ before stage $i + 1$
begins. In particular, $U S X$ means that the complete $S$ Taylor
polynomial is constructed first and a fresh $U$ acts on that polynomial;
it is not a joint bivariate Taylor projection. Stored jet entries are
chronological replay recipes only. A retained pre-series occurrence
image may be kept as diagnostic provenance, but it is not the physical
atom on which a later operator acts. Conditional reconstruction replays
one recipe, completes its Taylor jet at the active factor boundary, and
only then advances to the next recipe or multiplies the factor into its
quotient.

There is one canonical homogeneous full-route chart for the CFF kernel
and its numerator. It is selected from the active connected component
and its independent soft boundary carriers; an incoming child chart may
be validated and converted, but it does not select a different parent
chart. If $cal(C)_r \[ N \]$ is the numerator prepared in that
ordinary-$U$ full-route chart and $q_a$ are the independent soft
carriers, then

$ Phi^S \( s \) N = cal(C)_r \[ N \]\|_(q_a mapsto s q_a) . $

Dependent internal routes already expressed through $cal(C)_r$ are not
solved or substituted a second time. Carrier substitutions must not
revisit the right-hand sides of those routed momenta, since that would
introduce a spurious second power of $s$. Consequently,

$ Phi^S \( 1 \) N = cal(C)_r \[ N \] $

exactly. The massive-top-bubble numerator verifies this before
tensor/gamma normalization in the canonical full route `[e6,e8]`, with
`e6` as its independent carrier. This equality is a chart invariant
only: it does not by itself prove $S U = U S$. A complete comparison
must first route every dependent edge representative into the same
homogeneous chart. In the top-bubble factor, `e5` and `e6` are opposite
representatives of the same boundary momentum; comparing them as
independent symbols produces a spurious nonzero $S U - U S$ remainder.
After routing both to `[e6,e8]`, the production soft remainder, the
physical numerator, and the independent $k^2$, $k #h(-0.1667em) dot.op q$,
and $q^2$ oracles vanish exactly through the requested soft order. This
canonical equality must not be replaced by a second numerator-routing
convention.

`Q3` likewise has one structural convention: it is always a function of
exactly two arguments. Its second argument is one representation object:

- `Q3(edge,mink(4))` for the represented unindexed spatial vector;
- `Q3(edge,mink(4,index))` for its automatically indexed tensor-network
  image;
- `Q3(edge,cind(i))` for an explicit Cartesian component.

Thus these are different values of the same second slot, not different
`Q3` arities or conventions. Bare `Q3(edge)`, arbitrary second
arguments, and non-four-dimensional Minkowski representations are
invalid. Replacement code uses exactly `Q3(edge,W_.x_)`; `W_.x___` is
not a compatibility path.

Finally, `ct_marker` is unit-valued forest-history metadata, not part of
the physical deformation algebra. Each deferred sum or active sector
retains one accumulated marker history beside its marker-free physical
terms. Neither $Phi^U$, $Phi^S$, a Taylor projector, conditional
reconstruction, nor occurrence replay may transform that marker history.

Materialization first forms the signed marker-free $U + S - U S$ result.
Only after branch combination and the forest subtraction sign are
complete is the current approximation marker attached once and the
accumulated history updated. Public selected integrands and their
ordinary-reference channel carry the same history, while conditional
sources, jet recipes, and occurrence images remain marker-free.
Disconnected products multiply their accumulated histories. Stripping
markers must recover exactly the marker-disabled physical atom.

== Authoritative implementation-state addendum: tree-bridge boundary flow
<authoritative-implementation-state-addendum-tree-bridge-boundary-flow>
Conditional factorization distinguishes causal selector ownership from
the direction needed to interpret a split child boundary. If a factor
owns an ordinary tree bridge, that bridge carries no causal denominator
and its two selector sectors have already combined through

$ theta \( q^0 \) + theta \( - q^0 \) = 1 . $

It is therefore contracted during that factor's CFF generation and
contributes its ordinary propagator exactly once. The contraction
records the bridge as undirected in the generated block, but a child
meeting the same edge still requires a direction for its external-energy
shift. Before recursing, restore the graph's deterministic
`Orientation::Default` flow only on tree edges owned by the current
factor. Both disjoint children then see the same global bridge flow,
while opposite half-edge incidence supplies the appropriate local sign.
This routing convention introduces neither a selector nor a second
orientation sum.

Do not default arbitrary undirected split edges and do not relax the
normalized CFF generator's missing-boundary-direction error. A non-tree
split edge must still inherit its direction from the containing causal
quotient; otherwise the conditional factorization is incomplete and must
fail contextually.

== Authoritative implementation-state addendum: conditional evaluator domain
<authoritative-implementation-state-addendum-conditional-evaluator-domain>
This section has precedence over any earlier implementation note
suggesting that exceptional (S) and (US) branches can remain on the
flattened summed CFF through numerical evaluation. The ordinary (U)
branch remains on that established representation. A branch containing
(S) must instead use the quotient-first conditional residue
representation described above.

The distinction is forced by the DGSE reduced-parent hard limits. For
the degree-one child with internal edges (e5,e6), the four mixed
supports

$ \[ e 0 \, e 5 \, e 6 \] \, quad \[ e 4 \, e 5 \, e 6 \] \, quad \[ e 1 \, e 5 \, e 6 \] \, quad \[ e 5 \, e 6 \, e 8 \] $

have the correct partial soft images

$ D_i \( s \) = E_5 + E_6 + s B_i . $

Applying (S\_0) only after those factors have been flattened produces
$\( E_5 + E_6 \)^(- 4)$, removes three quotient denominators, and
changes five reduced-parent hard rays from approximately (-1) to exactly
(+2). Freezing the mixed factors is not a valid repair: it omits the
boundary derivative required by the matched four-dimensional residue and
destroys the top-bubble soft cancellation. The conditional residue
instead contains one child denominator and three untouched quotient
denominators.

A conditionally acyclic product need not lift to a globally acyclic
root-CFF orientation. Such a product is nevertheless a physical term of
the local counterterm and must be included in the numerical orientation
sum. The conditional factor tree therefore exposes a separate
evaluator-support sideband. For each complete product, merge only the
directions on every factor's `owned_selector_edges`. These edge sets
must be disjoint and must partition all paired non-tree, non-initial-cut
selector edges. In the full signature, initial-state-cut edges use their
fixed default direction and all other non-selector, tree, bridge, and
unpaired edges remain undirected. Exact signatures are deduplicated.
Occurrence namespaces are never interpreted as evaluator orientation
IDs.

The root CFF orientation list remains the sole basis for CFF generation,
threshold geometry, tropical sampling, public orientation labels, and
explicit root-orientation selection. Only the parametric evaluator for
the affected local atom receives the stable union

$ cal(O)_(upright(e v a l)) = cal(O)_(upright(r o o t)) union cal(O)_(upright(c o n d i t i o n a l)) . $

Root orientations retain their existing indices as the prefix of that
union. Threshold evaluators which do not contain an exceptional local
branch continue to use the root list. A per-residue local or threshold
atom which does contain such a branch receives its own conditional
support. Summed evaluator modes must include the full union. Sampling or
requesting one root orientation must fail contextually while
non-liftable supplemental sectors are active, because there is no
canonical one-to-one attribution of those sectors to a root orientation.

The massive-top-bubble factor supplies the structural oracle: its
fourteen quotient alternatives and two conditioned child alternatives
yield exactly twenty-eight complete selector signatures, with (e3,,e8)
directed and every other edge undirected. Numerical acceptance requires
both sides of the representation issue: the top-bubble (q\_g) fit
improves by two powers, and all five DGSE reduced-parent hard rays
return below the UV pass threshold.

== Authoritative correction: complete conditional products use the root domain
<authoritative-correction-complete-conditional-products-use-the-root-domain>
This section supersedes the preceding claim that a completed conditional
factor product can require a non-root evaluator orientation, and
therefore supersedes the residue-wide evaluator-domain union proposed
there.

The public CFF contracts only graph-theoretic tree bridges, not every
edge outside the chosen loop-momentum basis. The massive-top-bubble
topology is bridgeless: its ordinary root CFF consequently has thirty
fully directed, globally acyclic orientations on (e3,,e8). The
conditional factor tree has fourteen acyclic quotient choices and two
acyclic child choices. Its twenty-eight completed products are a strict
subset of those thirty root orientations, rather than supplemental
orientations.

More generally, under the validated laminar factor-tree construction,
every complete conditional product lifts to a globally acyclic root
orientation. If its merged directions contained a directed cycle, choose
the lowest factor containing that cycle. The cycle either lies wholly in
one child, contradicting the inductive hypothesis, or survives
contraction as a directed cycle in that factor's quotient, contradicting
the quotient CFF's acyclicity. Boundary directions are inherited
consistently, and graph-theoretic tree bridges remain undirected in both
representations. A quotient orientation may fail to define a
#emph[unique] root lift before child choices are supplied; this does not
mean that a completed product lacks a root lift.

Conditional orientation signatures remain useful provenance and a
structural validation oracle, but they must not enlarge the numerical
evaluator domain. For every residue, verify that each completed
signature is contained in the ordinary root domain and fail contextually
otherwise. Numerical evaluation, threshold geometry, tropical sampling,
public labels, and explicit orientation selection all continue to use
the ordinary root orientations alone. Remove the staged runtime and
standalone-schema plumbing that treated conditional signatures as
supplemental evaluator sectors unless a future topology provides a
concrete counterexample to the containment proof.

For the top-bubble acceptance failure, evaluator support is therefore
not the current suspect: the measured child-only fits are

$ a_U = - 2.0011044489 \, #h(2em) a_H = - 2.0011035828 . $

The implementation must next compare the separately retained (U), (S),
and (US) terms on the identical numerical soft ray and locate where the
required leading cancellation in (S-US) is lost. A symbolic assertion
that (H-U) is nonzero is necessary but insufficient; the branchwise
leading powers and coefficients must demonstrate the expected two-power
improvement before the physical profile can pass.

The conditional factor tree is not required to reproduce the undeformed,
selector-valued root CFF term by term, or even to define a separately
proved raw CFF identity. Write the orientation-summed deformed root and
conditional representations as (R(s)) and (C(s)). For the degree-two top
bubble, the operator-level compatibility condition is precisely

$ J_s^(lt.eq 1) #scale(x: 120%, y: 120%)[\[] R \( s \) - C \( s \) #scale(x: 120%, y: 120%)[\]] = 0 . $

The two root-only sectors may therefore be absent from (C(s)) provided
their constant and linear soft content is redistributed among the
changed conditional residues. A difference beginning at order (s^2) is
immaterial to (S\_1). Tests must compare the completed,
orientation-summed (S) and (US) jets with residues of their 4D-expanded
counterparts; a raw pre-jet equality is neither a substitute nor an
acceptance requirement.

== Authoritative diagnostic: quotient-first rebuilding deletes the required soft jet
<authoritative-diagnostic-quotient-first-rebuilding-deletes-the-required-soft-jet>
The branchwise massive-top-bubble profile locates the first numerical
failure more precisely than the structural factor-tree tests. On one
identical summed-orientation soft ray, the undeformed root and ordinary
branch have the expected reported scaling near (-2), while the completed
occurrence-aware soft jet #emph[before] conditional rebuilding has the
same leading coefficient with the opposite forest sign:

$ a_(upright(r o o t)) = - 2.00111449 \, #h(2em) a_(S \, thin upright(p r e - r e b u i l d)) = - 2.00110944 . $

At the first and last sampled points, their imaginary parts combine as

$ 1.7181377 thin 10^(- 7) - 1.7184258 thin 10^(- 7) = - 2.88 thin 10^(- 11) \, $

$ 1.7462234254 thin 10^8 - 1.7462234257 thin 10^8 = - 2.93 thin 10^(- 2) . $

The remainder is two powers softer. Thus the occurrence-aware
deformation, its physical-mass handling, its numerator route, and the
(S\_1) Taylor projection have already produced the low-order polynomial
required by the physical acceptance test.

The final `rebuild_conditional` call changes that same completed branch
to

$ a_(S \, thin upright(c o n d i t i o n a l)) = + 0.00135784 \, $

and changes the completed (US) branch similarly to approximately zero
soft scaling. The reconstructed branches are about two powers too soft
to cancel the root and ordinary contributions. The failure therefore
begins at the 30-orientation root-to-28-product reconstruction, not in
CFF occurrence ownership, the canonical momentum route, the Taylor
endpoint, or evaluator orientation support.

Concretely, the current conditional tree computes

$ sum_(alpha = 1)^14 Q_alpha \( q \) thin J_s^(lt.eq 1) sum_(beta = plus.minus) B_(beta divides alpha) \( s q \) \, $

but implements no redistribution which would make it equal to

$ J_s^(lt.eq 1) sum_(omega = 1)^30 R_omega \( q \, s q \) . $

The two expressions agree only on the undeformed physical diagonal if
the appropriate causal identities are used; that does not imply equality
after the child occurrence of the boundary carrier is scaled
independently of its co-graph occurrence. In particular, root
orientations 0 and 29 have no conditional product, and the current 28
changed products do not inherit their constant or linear jet content.

The immediate implementation audit must therefore treat the completed
full-root occurrence-aware jet as the reference production value. Before
removing conditional rebuilding, re-run the DGSE reduced-parent and
simultaneous-hard tests with both (S) and (US) retained on that same
representation. If those limits fail, repair how a #emph[containing]
operator acts on the completed child polynomial; do not replace the
child polynomial by the already-falsified quotient-first atom. The
conditional tree may remain a topology/provenance diagnostic only if it
is not used as a physical branch.

The same-component overlap must remain on the same representation as the
soft branch. Preserving only (S) while rebuilding (US) still leaves the
top-bubble fit at approximately (-2). Preserving both completed
full-root branches makes the existing non-vacuous child-only acceptance
test pass: the ordinary control retains raw behavior near (q\_g^{-5}),
while the complete soft-refined result behaves near (q\_g^{-3}),
i.e.~the reported bulk scaling improves from approximately (-2) to
approximately (0). This is also a direct numerical check of the ordered
definition (US=U(SX)); (S) and its overlap may not be evaluated in
different residue representations.

== Containing-parent discriminator: a full-root child is necessary but not sufficient
<containing-parent-discriminator-a-full-root-child-is-necessary-but-not-sufficient>
Keeping both the child (S) branch and its same-component (US) overlap on
the completed full-root CFF fixes the massive-top-bubble soft limit, but
applying the same representation indiscriminately through every
containing operation fails the DGSE hard limits. In the matched DGSE
profile, the full-root result has

$ a_(upright(s o f t)) = + 4.06730619 \, #h(2em) R^2 = 0.99965276 \, $

while five of 27 resolved UV limits fail. Each failed ray fixes the
child edge (e5) and scales a containing or complementary loop direction;
its fitted slope is (+2). The otherwise identical all-(U) control fails
zero of 27 limits. This rules out the universal two-line bypass of
conditional rebuilding: the full-root representation is the correct
completed child soft polynomial, but it is not automatically a valid
coordinate representation for every strict containing hard projection.

The next admissible implementation is therefore operation-boundary
aware:

- retain the physical child (S) and same-component (US) values on the
  occurrence-aware, orientation-resolved root CFF;
- apply a strict containing parent (U) directly to each completed child
  term, using the containing component's homogeneous route and
  active-loop scope;
- retain the chronological deformation program so that a later operator
  sees the actual child polynomial and does not infer its centre from
  its final algebraic form;
- never reconstruct the child (S) polynomial from quotient-first
  factors, since the top-bubble result has already disproved that
  equality.

For each completed child history (B), any proposed conditional parent
chart must satisfy, coefficient by coefficient,

$ J_u^(lt.eq d_Gamma) Phi_Gamma^U {N_(Gamma \/ gamma) [C_(upright(c o n d)) \( B \) - R_(upright(r o o t)) \( B \)]} = 0 \, $

without summing distinct orientations and with the same canonical route,
mass provenance, selector, and tensor contraction. For the available
top-bubble and DGSE strict parents, (d\_), so the constant coefficient
alone is a sharp first discriminator. A passing hard profile is not
proof of this identity; the branchwise recurrence and the matched
profile must both pass.

== Corrected representation verdict: retain the orientation-resolved CFF
<corrected-representation-verdict-retain-the-orientation-resolved-cff>
The branchwise and hard-limit diagnostics rule out the current
conditional reconstruction, but they do #strong[not] imply that exact
raised residues are a new requirement of the soft scheme:

- quotient-first conditional rebuilding deletes the constant and linear
  soft jet needed by the massive-top-bubble limit;
- retaining the completed root child jet restores that soft limit, while
  the current containing-parent route/scope then violates strict hard
  projections in DGSE;
- combining a full-root child with a conditional strict parent fails
  both the exact outer-pair recurrence and fourteen of twenty-seven DGSE
  UV rays;
- discarding the strict degree-zero parent image of the refinement,
  although its exact four-dimensional value is zero, still fails the
  same fourteen three-dimensional rays.

The exact four-dimensional DGSE topology independently establishes

$ U_(Gamma \, 0) [\( H_gamma - U_gamma \) C_(Gamma \/ gamma)] = 0 . $

Ordinary UV counterterms already generate raised four-dimensional
propagator powers. The established 3D implementation does not construct
their residues explicitly. Instead it generates the ordinary CFF per
orientation, substitutes the corresponding on-shell energies in the
numerator, applies the complete UV deformation to that moving 3D
expression, and only then extracts its Laurent jet. Derivatives of pole
locations, numerator energies, causal surfaces, and the 3D measure are
therefore supplied by the chain rule within each orientation. Raised
propagator powers alone cannot explain a failure first seen for soft
counterterms.

The genuinely new requirement is the ordered composition of deformations
with different centres and grading rules. For every orientation (rho),
the production path must implement

$ H_(delta \, rho) C_rho = U_(delta \, rho) C_rho + S_(delta \, rho) C_rho - U_(delta \, rho) #h(-0.1667em) (S_(delta \, rho) C_rho) \, $

and a containing operator must act on each complete signed child term on
that same orientation-resolved representation. The orientation selector
itself is fixed; only the component-owned momenta, energy shifts,
masses, and surfaces follow the selected deformation. No counterterm may
rely on cancellation between different orientations to become locally
subtracted.

The failed full-root experiment therefore diagnoses an incorrect
containing route or active-variable scope in the retained replay, not
missing raised-residue data. The failed conditional experiment diagnoses
a change of mathematical CFF representation under an exceptional
component-local deformation. Repair those operation scopes on the
established full CFF path; do not convert the result to quotient-first
denominators after the child jet has been formed.

An independent positive-energy residue oracle with a linear energy
numerator now tests this mechanism explicitly. For

$ D = k_0^2 - E^2 \, #h(2em) A = 2 \( k_0 p_0 - ell p_1 \) \, #h(2em) N_0 = a + b k_0 \, #h(2em) N_1 = c + d k_0 \, $

the first two coefficients of

$ frac(N_0 + s N_1, D^2 + s A D) $

at the positive pole give

$ R_0 = - frac(a, 4 E^3) \, #h(2em) R_1 = - frac(c, 4 E^3) + frac(p_0 b, 8 E^3) + frac(3 ell p_1 a, 8 E^5) . $

Multiplication by a strict parent's positive-pole factor yields

$ P_0 = - frac(a, 8 F E^3) \, #h(2em) P_1 = - frac(c, 8 F E^3) + frac(b, 16 E^3) + frac(3 ell p_1 a, 16 F E^5) . $

The terms proportional to derivatives of the numerator and cograph
remain a useful oracle for a future implementation which derives local
3D counterterms from expanded 4D quantities. They are not a Phase-1
production gate and do not justify replacing the per-orientation CFF
construction.

The current conditional factor tree remains useful as topology and
ownership provenance, but it is not a valid production replacement for a
nested soft atom. Phase 1 continues on the established
orientation-resolved CFF path.

== Canonical spatial-momentum representation
<canonical-spatial-momentum-representation>
There is one persistent abstract spatial-vector convention:

$ Q 3 \( e \, "mink" \( 4 \) \) quad upright("or") quad Q 3 \( e \, "mink" \( 4 \, mu \) \) . $

The edge argument is a non-negative integer graph-edge ID. A bare
`Q3(e)` is invalid. Explicit components are stored as `Q(e,cind(i))`,
not as a second `Q3` convention. Internally projecting `Q3` onto
`cind(0)` normalizes to zero, while projecting onto a spatial `cind(i)`
normalizes to `Q(e,cind(i))`. Malformed edge arguments must remain
visible long enough for graph validation and may not be normalized into
`Q` first.

== Rejected investigation: summed exact 4D-to-3D bridge
<rejected-investigation-summed-exact-4d-to-3d-bridge>
#strong[Do not implement this bridge.] The following records the
investigated approach so that its failure mode remains explicit. It was
rejected because a summed residue hosted on one orientation does not
subtract any individual orientation and therefore cannot be used for
local sampling, profiling, contour deformation, or residue selection.

The experiment was rejected: an orientation-summed residue cannot serve
as an orientation-local counterterm. Hosting it on one selector only
prevents duplicate counting; it does not subtract the other orientation
channels. Consequently it is unsuitable for local sampling, profiles,
contour deformation, selected cuts, or residue-specific evaluation. No
part of this bridge is connected to production.

== Corrected Phase-1 implementation boundary
<corrected-phase-1-implementation-boundary>
Reuse the existing per-orientation CFF representation and extend its
retained deformation history so that:

- #block[
  #set enum(numbering: "(A)", start: 19)
  + deforms only occurrences owned by its connected component and
    preserves physical masses;
  ]
- (US) applies a fresh ordinary deformation to the complete (S)
  polynomial on the same orientation;
- a parent applies its selected deformation to every complete signed
  child term, with the child boundary carrier reclassified according to
  the parent loop system;
- selectors and unowned co-graph occurrences stay fixed;
- normalization combines orientations only after each one has received
  its complete forest subtraction.

Because the counterterm is constructed directly on the CFF, cancellation
is an orientation-local requirement, not merely a property of the summed
orientation result. UV profiling must therefore be run for each
orientation individually, using generation-time orientation selection
where useful to shorten the diagnostic loop. A failure report must
retain the orientation signature/selector together with the routed hard
ray so that the faulty surface ownership or expansion route can be
isolated. The summed-orientation profile is an additional end-to-end
check and may not replace these per-orientation checks.

The positive-energy raised-pole test remains an independent algebraic
oracle for future 4D-derived counterterms. It is not wired into Phase-1
production.

This branch does not yet contain the separate work supporting original
numerators with more-than-linear dependence on an individual EMR energy.
Phase-1 therefore selects fixtures whose original numerator is already
at most linear in each EMR energy (in particular, quark-loop examples
are valid, whereas the gluon-loop contribution to a gluon self-energy is
excluded). This is a temporary fixture precondition only: do not add a
degree analyser, an early rejection, or any other production machinery
that would immediately be removed when that branch is merged.

== Authoritative mixed-surface and orientation verdict
<authoritative-mixed-surface-and-orientation-verdict>
This section supersedes every earlier research-log passage proposing
partial mixed-surface deformation, conditional-CFF regeneration,
supplemental orientations, summed-only validation, or a quotient-only
containing operator. Those passages are retained only as the record of
rejected investigations.

The connected soft operator owns a CFF denominator occurrence only when
both its stored topology shore and its energy/shift support are local to
the active connected component (including its crown carriers). Such an
occurrence receives the complete component-local spatial and temporal
deformation. A pure co-graph or mixed/spanning parent-CFF shore receives
no deformation and is restored literally from its occurrence tag. There
is no production `Partial` ownership class.

The rejected partial rule mapped a mixed denominator schematically as

$ sum E_(upright(c h i l d)) + sum E_(upright(c o g r a p h)) + q_0 mapsto "prepare" #h(-0.1667em) (sum E_(upright(c h i l d))) + s (sum E_(upright(c o g r a p h)) + q_0) . $

At $s = 0$ this erased the co-graph causal suppression. In the selected
DGSE orientation every causal denominator was mixed, so the pure child
atom lost all quotient causal denominators and five quotient-hard rays
incorrectly scaled as $+ 2$. This was not evidence for a quotient-only
outer operator. Production keeps the full containing `FullSubgraph`
ordinary operator required by the forest recurrence.

CFF generation records each concrete occurrence's generating shore and
support before algebraic cache deduplication. The occurrence lift
freezes every tagged inverse denominator before applying component-local
substitutions. Therefore an unowned occurrence remains exact even when
its mathematical denominator equals an owned occurrence; only the owned
occurrence receives the transformed image. Ambiguous or missing
provenance must still produce a contextual error rather than a global
replacement.

The direct 3D validation for this correction is orientation-local. For
DGSE, generation-time filtering of orientation `(-,-,0,0,+,+,+,0,+)`
leaves exactly `--00+++0+`; all 27 routed UV limits then pass
individually, including the five former $+ 2$ failures. Unfiltered
profiling covers 30 orientations times 27 limits. Some individual fits
have pre-asymptotic zero crossings in the original $10^8$-to-$10^12$
window; over a 33-point $10^8$-to-$10^16$ window all 810 resolve below
`-0.9`, while the isolated problematic channels converge to $- 1$ over
$10^12$-to-$10^16$. Acceptance must retain orientation labels and hard
rays, and must choose a window which actually reaches the orientation's
asymptotic regime rather than classifying a zero crossing as a
subtraction failure. Summed profiling remains an additional check only.

A pre-residue 4D deformation can give a different image for a selected
mixed pole. That fact is a negative algebraic oracle showing that
operation-before- residue does not define direct CFF ownership; it must
not be connected to the Phase-1 production operator. Exact 3D-versus-4D
local evaluation remains a later merge gate once local counterterms
derived from expanded 4D quantities are available.

== Tier-1 cleanup of the rejected conditional replay
<tier-1-cleanup-of-the-rejected-conditional-replay>
The quotient-first conditional factor tree and its deferred replay
program were an experimental implementation of the representation change
rejected above. They are not a dormant alternative production path:
rebuilding a completed soft jet on those factors deletes the constant
and linear content needed for the top-bubble cancellation and changes
the evaluator domain. The module, replay source, replay scopes,
factor-local deformation callback, completed- factor direct-map branch,
and tests whose only oracle was that representation are therefore
removed.

The retained jet history records only the connected component, operator
kind, Taylor parameter, and diagnostic label needed by the active
full-CFF forest recurrence. A containing operator continues to act
directly on every completed signed child term in the same
orientation-resolved CFF. The serialized `supplemental_orientations`
field and the broader CFF-generation APIs remain in place for schema
compatibility and possible future topology work; Phase 1 does not
populate them through conditional replay.

== Current direct-CFF refinement: contraction-transition ownership
<current-direct-cff-refinement-contraction-transition-ownership>
The later orientation-local audit supersedes the preceding claim that a
mixed parent-CFF occurrence must always remain literal. That rule
preserves the hard behavior but omits the child soft jet. Conversely,
deforming every mixed factor makes all four DGSE denominators coalesce
at the exceptional point and removes three quotient suppressions. The
missing distinction is not the cached surface, shore containment, or
edge support alone; it is the role of a factor in one concrete CFF
contraction branch.

For a root-to-leaf CFF contraction history with current vertex partition
$Pi$, define

$ b_delta \( Pi \) = \# { B in Pi divides B inter V_delta eq.not nothing } . $

A transition which merges blocks $A$ and $B$ belongs to the component
contraction tree exactly when

$ A inter V_delta eq.not nothing \, #h(2em) B inter V_delta eq.not nothing . $

It then decreases $b_delta$ by one. Restricting the full contraction
dendrogram to $V_delta$, deleting outside leaves, and suppressing unary
merges leaves precisely these transitions. Hence every path contains
$b_delta \( Pi_(upright(s e e d)) \) - 1$ component factors, or
$\| V_delta \| - 1$ for the singleton seed partition of the root CFF.
This is topology neutral and also states the correct generalization when
a future factor starts from precontracted blocks.

For DGSE and $delta = { v_2 \, v_3 ; e_5 \, e_6 }$, the selected
orientation `--00+++0+` has the chain

```text
{v2}+{v1}
{v1,v2}+{v5}
{v3}+{v4}
{v3,v4}+{v1,v2,v5}
```

Only the terminal transition joins two blocks which both meet the child.
Its mixed denominator is the component-reducing causal factor and
receives the complete deformation

$ D_i \( s \) = I_i \( k + s q \) + s B_i \( q \) \, $

where $I_i$ is its component-internal energy current and $B_i$ its
complementary boundary flux. This occurrence is classified as
`Complete`. The other three mixed denominators are spanning/quotient
occurrences. When one contains a component-internal energy, that
internal dependence still receives $I_i \( k + s q \)$, because the
parent CFF is not factorized into separate component and co-graph CFFs,
while its complementary boundary flux $B_i \( q \)$ remains literal.
Such an occurrence is classified as `InternalOnly`. An occurrence
containing no component-internal energy remains entirely literal and is
classified as `None`. Across all 30 acyclic DGSE orientations and 70
root-to-leaf monomials there is exactly one `Complete` component factor
per monomial. Shore restriction alone is ambiguous, and 24 shared CFF
nodes have different roles on different outgoing branches.

The minimal sufficient provenance is therefore the selected block and
its actual adjacent contraction partner on each parent-to-child
transition, plus the remaining partner at a terminal leaf. A node-only
occurrence tag is insufficient because its inverse denominator is
factored outside the child sum. The direct-CFF representation must
distribute that factor over outgoing branches and give each transition a
distinct occurrence identity before Symbolica collapse. Equal
mathematical denominators remain equal after the unit-valued transition
tags are stripped.

The active implementation records only this topology-neutral transition
metadata during the existing CFF generation, branch-expands the
occurrence-tagged form, applies the three deformation modes above, and
keeps the ordinary parent operator as `FullSubgraph`. The focused
orientation's 27 hard rays, all 810 unfiltered orientation/ray pairs,
and the matched soft profile pass. Failures must still be diagnosed on
the same orientation before any summed result is considered.

== Nested-jet representation audit
<nested-jet-representation-audit>
The contraction-transition rule solves #emph[which] causal factor a soft
operation owns, but a unit-valued tag on a flattened inverse is not by
itself sufficient to compose differently centred completed jets. After
an ordinary projection,

$ J_u \[ f \( u \) g \( u \) \] $

contains product-rule cross terms. The tag originally attached to
$f^(- 1)$ has factorized away from its coefficient, so replacing that
tag cannot freeze or transform the exact completed occurrence
independently of $g$. In particular, using the soft `Complete`/`None`
classification as a record of which factors a previous ordinary $U$
changed is invalid: $U$ acts on every factor carrying the active loop
variables, whereas $S$ owns only the component-reducing contraction
transitions. Silently restoring the cached pre-jet denominator would
corrupt nested $U \/ H$ terms.

The mathematically clean direct-CFF representation is therefore to
retain an opaque occurrence factor

$ cal(D)_omega \[ D_omega \] quad arrow.l.r quad D_omega^(- 1) $

for each orientation/branch occurrence $omega$ until every nested
deformation has been composed. Ordinary $U$ transforms the stored
argument $D_omega$ wherever its variables require it. Soft $S$
transforms that argument only for component-reducing transitions and
holds every other occurrence wrapper fixed while applying the component
map to the numerator, tree factors, and other exposed structures.

Each retained forest term then stores the fully composed multivariate
deformation together with its ordered jet list, but does not project a
jet during composition. At materialization the occurrence wrappers are
lowered back to inverse causal denominators and the jets are evaluated
in their original order, from the innermost child outwards. For two
stages this uses

$ J_s thin Phi^S \( s \) thin J_u thin Phi^U \( u \) X = J_s J_u thin Phi^S \( s \) Phi^U \( u \) X \, $

where the second equality only moves coefficient extraction past a later
substitution which is independent of $u$; it does not exchange the two
deformation maps or assume $U S = S U$. Thus all product-rule cross
terms are generated by the final multivariate expression, while the
occurrence identity is still available when each deformation is
composed. This remains a direct, per-orientation operation on the parent
CFF and introduces no 4D reconstruction.

Until this representation is implemented, a completed prior jet followed
by a soft operation which must leave one of its transformed causal
occurrences literal must fail closed rather than use the cached
denominator as a surrogate.

== Authoritative implementation-state addendum: immutable edge scope, materialization boundary, and canonical routes
<authoritative-implementation-state-addendum-immutable-edge-scope-materialization-boundary-and-canonical-routes>
For the direct 3D CFF operator, ownership of every exposed noncausal
inverse-energy or tree-propagator factor is determined only by its
immutable physical graph-edge ID. In
`LOCAL_3D_OSE_EDGE(e_phys,e_basis,mass_kind,D)`, `Den(e_phys,...)`, and
raw `OSE(e_phys)`, `e_phys` is the ownership key; the mutable
expansion-basis edge is never an ownership key. Before a connected
deformation, freeze occurrence-wrapped causal inverses first and then
freeze every remaining edge-local noncausal factor whose physical edge
is outside the active component. Route preparation and the selected
deformation act only on the exposed owned factors. Restore each frozen
factor literally before Taylor projection, so it is constant in that
stage but still participates in the complete product. The same factor is
exposed when a containing component owns its physical edge. A prepared
`OSE(basis,D)` without an immutable physical-edge ID, a malformed or
out-of-range ID, or a private freeze wrapper surviving to
materialization is a contextual error.

Materialization is a one-way boundary for occurrence-local composition.
It lowers the exact-current occurrence wrappers and projects the stored
jets; the resulting coefficients contain numerator/denominator
product-rule cross terms which cannot in general be refactorized into
the original owned occurrences. Unit-valued occurrence tags or cached
pre-jet denominators do not restore that lost factorization.
Consequently, applying a fresh eager `T` mapper to a materialized atom
globally reinterprets the completed CFF and is not an oracle for the
orientation-local action of a containing operator. Production composes
every later deformation on the retained deferred term before
materialization and projects the jets in chronological child-to-parent
order. Structural tests inspect that retained history; hard-limit
acceptance is supplied by the real per-orientation
`profile ultra-violet` path, with the summed profile as an additional
check.

A forest node's nominal `current.lmb()` is enumeration-path metadata,
not the semantic chart of its counterterm. The same connected component
may otherwise appear with equivalent representatives such as `[e2,e3]`
and `[e2,e4]`. Both 3D and 4D operations therefore first use the
graph-canonical full route or the compatible component subroute induced
from it, and every accepted loop coordinate must be homogeneous with
respect to graph-external momenta. The 3D full chart is derived
deterministically from that component chart and the operation's
genuinely independent boundary carriers; an incoming child chart is
validated and converted but does not choose a different parent chart. In
4D, the completed child records the coordinates on which either its
selected atom or its all-ordinary reference actually depends. Physical
membership of an otherwise absent coordinate representative in the
contracted child is not an additional chart constraint. If the canonical
candidate cannot retain a recorded coordinate, enumerate only compatible
homogeneous bases and choose the edge-lexicographic minimum. This
deterministic fallback is permitted only for such genuine retained
coordinates. Never select the nominal forest LMB or a
traversal-dependent first compatible basis, and reject any required
affine external-momentum shift contextually.

== Authoritative implementation correction: reuse the orientation-local MUV map
<authoritative-implementation-correction-reuse-the-orientation-local-muv-map>
The occurrence/shore and deferred-multijet designs above are superseded
by a direct audit of the existing, working 3D MUV pipeline. The ordinary
operator does not identify or replace individual E/H-surface
occurrences. It starts from the complete orientation-resolved parent
CFF, rewrites edge momenta into the active connected component's
compatible local LMB, and applies its deformation to that expression as
a whole. The soft operator must use the same scope and route; only its
deformation differs.

For every edge momentum resolved in the operation LMB as

$ Q_e = L_e + P_e \, $

ordinary large-loop scaling uses the existing representation
`L_e/t + P_e`, whereas soft scaling uses `L_e + s P_e`. The latter is
applied to both `Q3(e,mink(...))` and explicit `Q(e,cind(...))`
occurrences, including the temporal external-energy shifts already
present in the parent CFF. Loop variables, on-shell energies, and
physical masses remain fixed in the soft deformation, and no UV mass is
introduced. The Taylor parameter is projected before the completed child
atom is returned, so a parent applies its fresh graph/LMB deformation to
that actual completed expression just as it already does for nested
ordinary counterterms.

Accordingly, production must not attach topology-shore identities,
freeze individual causal denominators, retain deferred jets, generate
component CFFs, or reconstruct a 3D counterterm from a 4D one. Equal
mathematical CFF factors remain equal. Component locality is supplied by
the graph subgraph and the compatible local LMB replacements already
used by MUV. Validation is per-orientation structural inspection and
per-orientation UV/soft profiling; summed-orientation checks are
additional rather than substitutes.

== Superseded implementation-state addendum: one parent-CFF convention
<superseded-implementation-state-addendum-one-parent-cff-convention>
The temporary normalized component/cograph CFF-block builders are not
part of the direct 3D counterterm design and have been removed. Their
proposed per-block orientation inventory, residue-count transport, and
generated-occurrence namespace encoded the rejected factorized replay
approach described in the superseded findings above. No
occurrence-orientation indexing convention remains in the exporter: it
consumes the already-selected graph-term orientations directly, and no
generated-block namespace or second orientation filter exists.

`generate_uv_cff` remains an unnormalized fixed-orientation graph
diagnostic. It may demonstrate, for example, that contraction turns a
root-acyclic orientation into a cyclic cograph, but local soft
counterterms never use its result as a component atom or multiply it
into the parent expression. Production starts from `Graph::cff` or
`Graph::cff_for_orientations`, attaches occurrence identity to that
parent CFF, and composes every U/S deformation directly on its retained
orientation-resolved terms. Contracted undirected edges receive a
temporary direction only inside the diagnostic builder before the edge
is removed; this is not a second orientation or momentum convention.

== Final direct-CFF rule: stage-relative pole energies
<final-direct-cff-rule-stage-relative-pole-energies>
The existing MUV inspection trace fixes the precise meaning of the
preceding direct-CFF rule. For a connected component $delta$ and a
completed orientation atom, first route only the reduced component
$delta \\ gamma$ through its compatible local LMB. If

$ Q_e = L_e + P_e \, $

then the soft deformation is applied to the exposed component dependence
as

$ L_e + P_e mapsto L_e + s P_e . $

The scaled generators are `Q3(e,mink(...))` and explicit
`Q(e,cind(...))` for physical edges in `dummy_less_full_crown(delta)`.
Although a compatible sub-LMB also retains the original graph-external
generators, those are spectators unless their physical edges belong to
this component crown. The crown is that of the complete current
component, not the reduced graph (whose holes are internal child
boundaries) and not an active disconnected-sector union.

A remaining one-argument `OSE(e)` is held literal during this stage. It
is a pole energy of the already retained parent residue, not an
independent component-external Taylor generator. Scaling it globally
would also deform the universal parent-CFF factors and co-graph causal
denominators. For the massive-top bubble this changes the intended
Laurent form from $s^(- 2) F \( s \)$ to $s^(- 4) F' \( s \)$, so the
fixed endpoint $- 1$ no longer selects the component's $j = 0 \, 1$
Taylor coefficients. In contrast, an active internal energy has already
been materialized by `start` as `OSE(lmb_id,D[Q3])`; its argument is
transformed through the routed `Q3` dependence.

This distinction is stage-relative rather than permanent ownership
metadata. An `OSE(e)` held fixed by a child becomes an internal energy
when a containing component owns edge $e$. The containing operation's
fresh `start` then materializes it as `OSE(lmb_id,D[Q3])`, after which
the containing deformation acts on its argument. This is how eager child
completion composes with a parent operation without surface-occurrence
IDs or deferred jets.

The direct 3D prescription is therefore:

+ retain the complete parent CFF separately for every orientation;
+ materialize the active internal energies with the existing `start`
  step;
+ route the reduced component with the same compatible local LMB
  machinery as MUV;
+ scale only the current component-crown `Q3/Q` generators, keeping raw
  parent-residue `OSE`, masses, selectors, and graph-external spectators
  literal;
+ multiply by $s^(- d)$, keep Laurent powers through $- 1$, and set
  $s = 1$;
+ apply the existing ordinary 3D MUV operator to this completed soft
  atom for the $U S$ branch, then combine $H = U + S - U S$.

This is an orientation-local CFF jet at fixed retained pole energies. It
is not obtained by taking a residue of the implemented 4D soft
counterterm, and Phase 1 does not claim per-term equality to such an
operation. Its current acceptance oracle is cancellation in every
orientation: in the top-bubble fixture the eight problematic ordinary-U
orientations scale near $- 2$, while all corresponding H orientations
scale near zero and improve by two powers; a containing degree-zero U
changes the result only at numerical-noise level. The nested H/H fixture
has 315 resolved orientation-level UV fits with no failure, while its
matched bare control fails 183.

Finally, the local-only forest algebra does not admit a positive-degree
ordinary-U parent above a positive-degree H child in general. The finite
$\( 1 - U \) S$ child refinement needs the containing soft completion
appearing in the paper's nested formulas; a degree-zero parent is
harmless because $H_0 = U_0$. Until the deferred
scheme-change/integrated policy exists, production therefore accepts
U/U, U/H, H/H, and H/U only for a degree-zero U parent, and rejects
positive-degree H-child/U-parent woods contextually.

== Authoritative correction: hard-dual soft deformation in the parent CFF
<authoritative-correction-hard-dual-soft-deformation-in-the-parent-cff>
The preceding stage-relative crown-scaling rule is superseded by the
hard-dual construction below. Its restriction on a positive-degree
ordinary parent is independent of that representation choice and remains
necessary, as shown by the mixed-wood counterexample at the end of this
addendum. Ordinary UV Taylor expansion already admits two exactly
equivalent homogeneous charts: scale external momenta and masses to
zero, or scale the component-internal loops to infinity. The soft Taylor
expansion has the same freedom, with a different mass grading.

Let a connected component integrand without its loop measure have
homogeneous engineering degree $h$, with internal loop momenta $k$,
component-external momenta $p$, and all dimension-one scales $m$ and any
existing terminal UV scale $M$ held fixed by the paper's soft limit.
Then

$ X \( k \, s p \, m \, M \) = s^h X \( k \/ s \, p \, m \/ s \, M \/ s \) . $

Consequently the two hard representations are

$  & upright("direct small-variable chart") & upright("dual hard chart")\
U & \( p \, m \) mapsto \( s p \, s m \) \, #h(0em) M upright(" fixed") & \( k \, M \) mapsto \( k \/ s \, M \/ s \) \,\
S & p mapsto s p \, #h(0em) \( m \, M \) upright(" fixed") & \( k \, m \, M \) mapsto \( k \/ s \, m \/ s \, M \/ s \) . $

The hard chart does not change the definition of $S$, expand a physical
mass, or introduce $M_(upright(U V))$. Co-scaling $m$ is precisely the
homogeneous representation of leaving that physical mass fixed in
$X \( k \, s p \, m \)$.

For an edge routed in the component LMB as

$ q_e = L_e + P_e \, $

the soft hard deformation is

$ q_e mapsto P_e + frac(q_e - P_e, s) = P_e + L_e / s \, #h(2em) m_e mapsto m_e / s . $

Thus it is a dilation of the internal part around the fixed
external-flow point $P_e$, not an affine change of loop basis. In the 3D
representation,

$ E_(m \/ s) \( L \/ s + P \) = s^(- 1) E_m \( L + s P \) . $

An E- or H-surface containing active energies therefore generates the
correct soft scaling of its spatial dependence and temporal boundary
shift automatically. No component crown, raw boundary `OSE`, tree
denominator, CFF surface, shore, or contraction occurrence is globally
rewritten. Mixed surfaces are handled in exactly the same
orientation-local manner as by the existing MUV hard-loop map: active
energies are deformed and spectator energies remain fixed.

With $L$ active residue loops and 3D homogeneous degree $h_(3 d)$,

$ s^(- 3 L) X \( k \/ s \, p \, m \/ s \) = s^(- d) X \( k \, s p \, m \) \, #h(2em) d = 3 L + h_(3 d) . $

If the direct soft expansion is $sum_(j gt.eq 0) X_j s^j$, its hard form
is $sum_(j gt.eq 0) X_j s^(j - d)$. The fixed Laurent endpoint `-1`
therefore selects exactly $j = 0 \, dots.h \, d - 1$, while $d = 0$ has
no soft branch. In the code's pre-inversion variable $x = 1 \/ s$, scale
active loop generators and component-local fixed mass scales by $x$,
apply the existing $x^(3 L)$ measure grading, replace $x$ by $1 \/ s$,
and retain strictly negative Laurent powers. There is no additional
manual factor $s^(- d)$.

The active OSE normalization before inversion is

$ "OSE" \( e \, D \) mapsto x^2 "OSE" #h(-0.1667em) (e \, D / x^2) $

inside the outer square root. This extracts the expected hard energy
factor and keeps the OSE argument regular for Symbolica. It differs from
ordinary MUV only in the finite argument: MUV retains its

$ x^2 "OSE" #h(-0.1667em) (e \, frac(M^2 x^2 + D - M^2, x^2)) $

rearrangement, whereas soft introduces no $M$.

Loop and OSE edge identities already provide all required CFF locality.
The only information lost by an eager Taylor projection is the support
of a free mass factor. For example, a fermion numerator may leave a
standalone $m$, and an inner ordinary projection may leave a coefficient
proportional to $m^2 - M^2$. Two disconnected components can use the
same model mass symbol, so a global mass substitution would act on the
wrong sibling. The minimal persistent metadata is therefore a
unit-valued mass-support tag anchored to a representative physical graph
edge of the connected component. It is attached only to masses in the
newly introduced component numerator and materialized active OSEs.
Ordinary U-generated $M$ receives the same support tag. A containing
operation rebases contained tags to its representative; a disjoint
operation leaves them unchanged. Soft scales only its active tag. The
tag is replaced by one only on the cloned final-evaluator expression,
after all forest operations have composed.

This mass-support tag is not CFF occurrence provenance. Equal causal
denominators remain mathematically equal, and there are no shore IDs,
surface ownership tables, deferred multijets, component CFF products, or
4D-to-3D counterterm reconstruction.

For disconnected components with

$ I = T \( p \) A \( k_A \, p \, m \) B \( k_B \, p \, m \) \, $

the two hard maps use independent loop supports and independently tagged
mass occurrences. Hence

$ S_B S_A I = T \( p \) \( a_0 + a_1 p \) \( b_0 + b_1 p \) \, $

including the $a_1 b_1 p^2$ term. The shared bridge $T \( p \)$ is
fixed. This is precisely the term lost by one global total-degree crown
Taylor series. For nested components, the child projection is completed
first; the parent then scales every loop and mass support contained in
the parent. The maps need not commute, and their eager chronological
order is exactly the forest product
$K_(upright(p a r e n t)) K_(upright(c h i l d))$. This establishes the
implementation order, but it does not by itself make every mixed
operator assignment locally finite. For $D_gamma = H_gamma - U_gamma$,
the change from a U/U forest to an H-child/U-parent forest is

$ R_(H \/ U) - R_(U \/ U) = - \( 1 - U_Gamma \) D_gamma C_(Gamma \/ gamma) . $

A child boundary monomial of degree $j lt.eq d_gamma - 1$ has
reduced-parent degree $d_Gamma - d_gamma + j$, which can reach
$d_Gamma - 1$. Thus a positive-degree ordinary parent can leave a
nonnegative reduced-parent hard limit even though $D_gamma$ is finite in
both the child-hard and simultaneous-hard limits. The existing exact
algebraic and graph-level counterexamples verify this obstruction.
Production therefore continues to allow a degree-zero U parent above H,
but rejects a positive-degree U parent until either the ancestor is also
assigned H or an explicit finite scheme-change/effective-vertex
counterterm is implemented. This is separate from, and in addition to,
the ban on integrated soft/PolePart mixtures.

The active implementation follows this rule directly on every parent-CFF
orientation. Its first validation set establishes:

- homogeneous $d = 0 \, 1 \, 2$ hard/soft micro-oracles, including a
  materialized OSE and a fixed tree/co-graph denominator;
- a nested child whose spatial and temporal boundary dependence is
  softened while graph-external and raw parent-pole factors remain
  fixed;
- the exact two-component degree-two scalar spectacles fixture,
  including its summed and per-orientation UV profile and matched
  bare/MUV controls;
- consecutive soft components;
- a nested massive H/H top self-energy with a non-vacuous
  per-orientation UV profile; and
- all 30 orientations of the massive-top bubble child, including the
  eight ordinary-U orientations with the expected approximately
  two-power soft improvement under H.

These checks derive S from the 3D CFF itself. The 4D identity above is
an algebraic proof and test oracle only; it is never used to construct
or map the 3D local counterterm.

== Verified mixed-wood and terminal-mass correction
<verified-mixed-wood-and-terminal-mass-correction>
The degree-qualified mixed-wood restriction above has now been tested by
temporarily enabling all four production rows (U/U), (H/U), (U/H), and
(H/H) in the same two-loop scalar fixture, retaining the four forest
families separately and measuring the child, reduced-parent, and
simultaneous hard rays. Exactly one obstruction occurs:

$ H \/ U : #h(2em) - \( 1 - U_Gamma \) H_gamma I = O \( lambda^0 \) $

in the reduced-parent ray, and the complete (H/U) forest has the same
logarithmic remainder. The other three rows are locally improved in
every required ray. This is direct production algebra, not a 4D-to-3D
comparison, and confirms that the original unqualified request to allow
a soft child below an ordinary positive-degree parent cannot be part of
the local-only Phase-1 acceptance. A logarithmic parent remains allowed
because (H\_0=U\_0); this is the mixed nesting exercised by the
top-bubble vertex.

The same audit corrected an earlier malformed terminal-mass oracle.
GammaLoop uses `mUVexp` as denominator provenance and `mUV` as the
actual vacuum mass:

$ "den" \( e \, q \, m_(upright(U V \, e x p))^2 \, q^2 - m_(upright(U V))^2 \) . $

An outer inverse-hard (U) keeps a free numerator `mUVexp` fixed, which
gives it the required direct-Taylor mass grading, while it sends the
actual terminal `mUV` to `mUV/s`. The terminal propagator then
reproduces the same (q#super[2-m\_{}];2) coefficient and is not
infrared-rearranged twice. The active regression now obtains this
structure from a genuine inner `t_raw` projection before applying the
outer operation; it no longer constructs an unphysical denominator with
`mUVexp` as its actual mass.

== Exact factorized materialization of the local soft operator
<exact-factorized-materialization-of-the-local-soft-operator>
The separate $U \( X \)$, $S \( X \)$, and $U \( S \( X \) \)$
expressions remain the definition and the small-oracle provenance of the
local operator. Production 3D materialization uses the exactly
equivalent linear form

$ H \( X \) = S \( X \) + U #h(-0.1667em) (X - S \( X \)) . $

This does not exchange the order of $U$ and $S$, and assumes neither
commutativity nor idempotence. The existing ordinary 3D projector is
linear for a fixed component, active sector, and compatible LMB: its
routing and mass scope substitutions are fixed ring homomorphisms, and
multiplication by the measure factor, Laurent truncation, and final
parameter evaluation are linear. Consequently,

$ U #h(-0.1667em) (X - S \( X \)) = U \( X \) - U \( S \( X \) \) $

with exactly the same component route and mass support as the two
separately materialized branches. A parent still receives the completed
child atom and applies its own factorized operation afterwards.
HedgePoset forest nodes and their signs remain separate; only the
internal representation of one completed $H$ atom is changed.

This factorization supersedes earlier implementation wording which
required production to retain two complete, nearly cancelling
$U \( X \)$ and $U \( S \( X \) \)$ atoms simultaneously. Focused
algebraic and inspect tests must still compare the conceptual
$U + S - U S$ branches independently on expressions small enough to
materialize them safely.

The immediate engineering reason is also mathematically aligned with the
subtraction: on the one-orientation Figure-B.1 fixture, separately
retaining the two UV branches produced a contracted scalar Add of
5,384,386,060 bytes. Symbolica stores Add lengths as 64-bit values but
stores a Mul/Fun child length as 32 bits, so wrapping that scalar
truncated the size by exactly $2^32$ and a later scalar-alias traversal
panicked. Applying $U$ once to the soft remainder allows equal Laurent
coefficients to cancel during the series projection rather than after
tensor contraction. Compact byte-size logging is retained at the local
projection, forest aggregation, final-numerator, and evaluator
boundaries. If a future graph still reaches the representation ceiling,
the next clean boundary is to keep signed forest-node atoms as separate
evaluator inputs and sum their numerical values per orientation; it is
not to reconstruct the 3D counterterm from a 4D expression or to hide
large branches in tagged function-map bodies.

== Signed forest-node evaluator decomposition
<signed-forest-node-evaluator-decomposition>
The factorized local operator reduced, but did not remove, the
representation ceiling in the complete one-orientation Figure-B.1
fixture. Its second contracted scalar was 5,241,104,034 bytes; the
enclosing Symbolica child length again retained only 946,136,738 bytes,
differing by exactly (2^{32}). The failure occurred after the direct-CFF
forest had been completed, in tensor network scalar reconstruction. It
is therefore neither evidence for a different soft deformation nor
permission to reconstruct the counterterm from the four-dimensional
representation.

HedgePoset now retains every signed final node contribution as a named
parametric evaluator term before its residue-wise analytical sum. The
latter remains available only for consumers which explicitly require one
symbolic expression, including legacy comparison and analytical
inspection. Amplitude evaluator construction consumes the named terms
directly. The existing multi-output evaluator machinery contracts and
optimizes each term independently; after the requested orientation or
orientation sum has been evaluated, the amplitude wrapper adds every
returned value. No output may be dropped, and an empty output list is an
error.

This changes only the association of the exact forest sum. It does not
split one connected approximation into independently marked (U), (S),
and (US) counterterms: a completed (H) node remains one signed evaluator
channel. It also retains the node key needed to diagnose any channel
which remains too large. If one node alone exceeds the representation
ceiling, that node may be partitioned further only along its existing
top-level additive terms while retaining the same node identity;
products are never split and no tagged function-map surrogate is
introduced.

Because exact symbolic cancellation between different forest nodes is
then performed numerically, acceptance requires the existing stability
escalation and matched per-orientation UV profiles. The Figure-B.1
regression must both complete evaluator construction and pass its
per-orientation profile; merely avoiding the 32-bit panic is not
sufficient.

== Superseding implementation-state addendum: forest-wide runtime aggregation
<superseding-implementation-state-addendum-forest-wide-runtime-aggregation>
This section supersedes the proposed runtime behavior in #strong[Signed
forest-node evaluator decomposition] above, while preserving that
section as the record of the investigated and rejected design. Signed
forest nodes must remain separately available through the computed
inspect/export path, with their node keys, parent keys, residue indices,
local projection paths, and 4D branch provenance. They are not
independently evaluator-ready runtime channels.

The evaluator-facing representation is the historical residue-wise sum
of the complete signed forest. HedgePoset first applies per-node
normalization and color collection, adds all compatible signed nodes,
and only then performs the final orientation-local four-vector split.
This ordering is required because temporal components in localized
integrated-counterterm polynomials can cancel between forest nodes
before an internal four-momentum is replaced by its orientation-resolved
energy and spatial parts. Moving that split across the forest sum
exposed spurious internal energy parameters and failed two sunrise
subgraph profiles by two powers.

Consequently, `ParametricIntegrands` again contains one `Integrands`
residue map, `AmplitudeDerivedData` again contains one completed
`all_mighty_integrand`, and evaluator construction receives exactly that
one atom. Standalone archives therefore retain their version-5
one-output layout; there is no additive-term count, logical-result
post-reduction, or pairwise forest reduction. The rejected multi-output
design first summed each term's contiguous orientation group and then
used an adjacent-pair tree across terms; that association was
deterministic, but it could not reproduce symbolic cross-node
cancellation and therefore is not a valid replacement for the
forest-wide sum.

This runtime aggregation does not weaken the requirement to test forest
terms one by one. Small algebraic oracles and computed UV-forest exports
continue to retain the four signed nested families independently before
any test regroups them. Production evaluation alone consumes their
completed symbolic sum.

== Canonical 3D route transitivity audit
<canonical-3d-route-transitivity-audit>
No additional retained-carrier or CFF-occurrence metadata is required
for nested local-LMB replacement. Let $B$ be the ordered cotree-edge set
of the graph-canonical full LMB and let $T = G \\ B$ be its spanning
tree. For any (possibly disconnected) subgraph $A$,
`try_compatible_sub_lmb` chooses the first $r \( A \)$-element subset
$R_A subset.eq B inter A$ whose deletion leaves every component of $A$
connected. Equivalently, start from $T inter A$, scan the edges of
$B inter A$ in reverse canonical order, and add an edge exactly when it
joins two current components. The rejected edges are $R_A$.

For nested subgraphs $C subset.eq P$, connectivity in the working forest
for $P$ is, at every point of that reverse scan, a coarsening of
connectivity in the working forest for $C$. Therefore, whenever an edge
is rejected from the child forest because its endpoints are already
connected, its endpoints are already connected in the parent forest as
well. Hence

$ R_C subset.eq R_P . $

The argument applies componentwise and includes self-loops. Since
`Spinney::with_scheme` induces every production component LMB from the
same graph-canonical full LMB, every Q3/Q carrier retained by a
completed child is necessarily a loop carrier of a containing parent.
The parent hard map thus scales it through the existing MSbar-style
local-LMB replacement. LMB generator ordering may differ, but the
physical carrier-edge set used by the replacement is nested, which is
the required invariant. This rules out the suspected
frozen-child-carrier failure without introducing shore IDs, occurrence
IDs, or a 4D reconstruction path.

== Final acceptance reconciliation
<final-acceptance-reconciliation>
The occurrence-ownership and ambiguity-rejection bullets in the original
plan, including the corresponding completion bullet, are superseded by
#strong[Authoritative implementation correction: reuse the
orientation-local MUV map] and #strong[Authoritative correction:
hard-dual soft deformation in the parent CFF];. They are retained above
only as the investigation record. Phase 1 does not add shore IDs,
occurrence IDs, frozen causal factors, component CFFs, or a fallback
global surface rewrite. Its final 3D locality criterion is that the
existing orientation-local MUV route acts on the complete parent CFF,
with the soft hard-dual changing only the deformation and Laurent
endpoint. The per-orientation profiles are the production acceptance
test for that rule.

Likewise, the one-second `sysinfo` telemetry paragraph belongs to the
explicitly superseded fixed-orientation investigation. It is not a
Phase-1 acceptance criterion and no process-monitor dependency or
scheduler change is introduced for the soft operator.

The two UV-mass roles are now kept distinct during 3D forest
composition. For a fresh active on-shell energy, the ordinary hard chart
uses

$ frac(m_(upright(U V))^2 r^2 + D - m_(upright(U V \, e x p))^2, r^2) \, $

so the terminal energy contains the actual vacuum mass while the Taylor
coefficient retains the expansion-mass provenance. If an active OSE
already contains `mUV`, it is an inner-U terminal energy: `mUV` and the
active loop are scaled together and its argument is divided by $r^2$,
without applying the MUV rearrangement again. This recognition is
OSE-local and therefore cannot rescale a disconnected sibling. It relies
on the current reserved-symbol invariant that `mUV` inside an active OSE
argument denotes a terminal vacuum mass; a future nonterminal use there
would require explicit provenance.

Finally, a `--per-orientation` UV profile must not construct its summed
curve by adding the already exported `f64` orientation values. Local
forest terms can cancel by substantially more than machine precision
even when every orientation was evaluated with arbitrary precision. The
summed curve is therefore evaluated independently with
`orientation=None`, which performs the orientation sum before the final
downcast. If any individual orientation required the arbitrary-precision
retry, the independent summed evaluation is also forced through that
precision. Per-orientation diagnostics remain unchanged and retain their
own stability metadata.

== Single completed 4D counterterm state
<single-completed-4d-counterterm-state>
This addendum supersedes the historical statements that a selected 4D atom
travels with an all-ordinary reference and that both determine its retained
expansion coordinates. Those statements remain above as the investigation
record. `Local4dCts` and `Full4dCts` now retain one completed atom; each parent
applies its selected operator to that atom, and expansion coordinates are
derived only from its actual momentum dependence. Only the selected result is
finalized and marked.

The `has_soft_ancestry` flag records whether a nonzero S branch has been
executed in the recursive contribution. It propagates through later operators
and disconnected products. A degree-zero H or a zero S branch does not
introduce soft ancestry, and later algebraic cancellation does not erase
existing ancestry. This history supports routing and contextual scheme checks
without constructing or comparing a second counterterm expression.

The orchestrator owns early validation of unsupported scheme assignments.
Direct `uv_limit` calls retain the corresponding inexpensive invariant errors
using `has_soft_ancestry`. Ordinary-scheme baselines needed by algebraic and
graph comparisons are constructed explicitly in tests. The selected forest
operator and its U, S, US, and combined branch provenance are unchanged.
