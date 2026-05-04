#set document(title: "FORM gamma and color simplification rules")
#set page(paper: "a4", margin: (x: 19mm, y: 18mm))
#set text(font: "Libertinus Serif", size: 10pt)
#set raw(tab-size: 2)
#set math.equation(numbering: "(1)")
#show heading.where(level: 1): it => {
  pagebreak(weak: true)
  it
}

#let source(url) = link(url)[#url]
#let form-manual = "https://www.nikhef.nl/~form/maindir/documentation/reference/html/manual.html#dx1-85012"
#let color-doc = "https://www.nikhef.nl/~form/maindir/packages/color/color.html"
#let color-src = "https://www.nikhef.nl/~form/maindir/packages/color/color.h"
#let color-examples = "https://www.nikhef.nl/~form/maindir/packages/color/color.tar.gz"
#let opera = "https://github.com/form-dev/form/blob/master/sources/opera.c"
#let gamma-rs = "https://github.com/alphal00p/gammaloop/blob/main/crates/idenso/src/gamma.rs"
#let color-rs = "https://github.com/alphal00p/gammaloop/blob/main/crates/idenso/src/color.rs"
#let sym-pattern = "https://symbolica.io/docs/pattern_matching.html"
#let sym-rust-id = "https://docs.rs/symbolica/latest/symbolica/id/index.html"
#let sym-rust-replace = "https://docs.rs/symbolica/latest/symbolica/id/struct.ReplaceBuilder.html"
#let spenso-doc = "https://docs.rs/spenso/latest/spenso/"
#let spenso-symbolic = "https://docs.rs/spenso/latest/src/spenso/network/library/symbolic.rs.html"
#let typst-syntax = "https://typst.app/docs/reference/syntax/"
#let typst-math = "https://typst.app/docs/reference/math/"

= FORM Gamma- and Color-Algebra Simplification Rules

This specification records source-backed replacement rules for FORM gamma-algebra simplification and the FORM `color.h` color-algebra package. The mathematical notation is independent of FORM token names except where a source symbol is documented literally. Open gamma strings always carry explicit bispinor indices; closed gamma traces are written as explicit bispinor contractions, not as a separate trace-index placeholder.

The Symbolica + spenso section uses the approved `chain` and `trace` tensor-pattern notation:
```rs
chain(bis(4,a),bis(4,c),
  gamma(in,out,p(2,mink(4))),
  gamma(in,out,mink(4,mu)))
```
and
```rs
chain(cof(Nc,i),dind(cof(Nc,j)),
  t(coad(Nc^2-1,a),in,out),
  t(coad(Nc^2-1,b),in,out))
```
This is a specification-level pattern notation. GammaLoop currently has internal helpers named `spenso::gamma_chain` and `spenso::gamma_trace`; those names are mentioned only when documenting the implementation source.

== Source Inventory

- FORM manual gamma-algebra section: #source(form-manual).
- FORM current source repository, gamma trace implementation: #source(opera).
- FORM package page for `color.h`: #source(color-doc).
- FORM package source `color.h`: #source(color-src). This source is not currently present in the `form-dev/form` GitHub tree; line references to `color.h` therefore use package-source line numbers recorded in this document and not GitHub anchors.
- FORM color example programs: #source(color-examples). The package page describes these as simple FORM programs that use `color.h`.
- GammaLoop/idenso current tensor simplifiers: #source(gamma-rs), #source(color-rs).
- Symbolica pattern matching documentation: #source(sym-pattern).
- Symbolica Rust `id` API documentation: #source(sym-rust-id), #source(sym-rust-replace).
- spenso Rust crate documentation: #source(spenso-doc), symbolic tensor-library source: #source(spenso-symbolic).
- Typst syntax and math documentation: #source(typst-syntax), #source(typst-math).

== Notational Conventions

Lorentz indices are written as $μ, ν, ρ, σ$ and are raised and lowered with $η^(μ ν)$ and $η_(μ ν)$. The dimension is denoted by $D$ in dimension-generic identities and fixed to $4$ in the four-dimensional FORM `trace4` rules. Repeated Lorentz labels are contracted.

Bispinor indices are $α, β, χ, δ$. An open gamma chain is written as $(Γ^(μ_1 ... μ_n))_α^β$. A closed trace is an explicit contraction $(Γ^(μ_1 ... μ_n))_α^α$.

Fundamental color indices are $i, j, k, l$. Adjoint color indices are $a, b, c, d, e$. Fundamental generators are $(T^a)_i^j$. The package constants are mapped as $C_A ↔$ `cA`, $C_F ↔$ `cR`, $T_R ↔$ `I2R`, and $N_c ↔$ `NR`; `NA` is the adjoint dimension. For ordinary $"SU"(N_c)$, $C_A = N_c$, $T_R = 1 / 2$, and $C_F = (N_c^2 - 1) / (2 N_c)$, but `color.h` is written in terms of invariants.

The word `product` denotes an ordered noncommutative string inside a chain. It is not used as FORM syntax.

= Gamma-Algebra Replacement Rules

== Overview

#table(
  columns: (2.4cm, 4.4cm, 6.8cm),
  [Rule], [Dimension], [Purpose],
  [Gamma collection], [all], [Collect adjacent bispinor-contracted gamma factors into an ordered chain.],
  [Chain joining],
  [all],
  [Join two open chains when the outgoing bispinor index of one equals the incoming index of the next.],

  [Metric contraction], [all or $D$], [Replace $γ^μ γ_μ$ and related repeated Lorentz slots by the active dimension.],
  [Trace closure], [all], [Turn an open chain with equal endpoints into an explicit bispinor contraction.],
  [Trace recursion],
  [$D$ or $4$],
  [Evaluate even traces recursively and annihilate odd traces without gamma-five insertions.],

  [Chisholm reduction], [$4$ only], [Reduce $γ^μ Γ γ_μ$ using four-dimensional identities.],
  [Gamma-five and epsilon], [$4$ only], [Use Levi-Civita and gamma-five branch logic in FORM `trace4`.],
)

== G1. Gamma-chain product and joining <rule-gamma-chain>

Mathematical identity:
$ γ(v)_α^β γ(w)_β^χ = (Γ^(v w))_α^χ $ <eq-gamma-chain-join>

For longer strings, the same rule constructs an ordered product inside the chain:
$
  (Γ^(v_1 ... v_m))_α^β (Γ^(w_1 ... w_n))_β^χ
  = (Γ^(v_1 ... v_m w_1 ... w_n))_α^χ
$ <eq-gamma-chain-product>

FORM source context:
```c
/* Trace4 receives a string of gamma matrices in 'instring'. */
number = *params - 6;
...
for ( i = 0; i < number; i++ ) *p++ = *m++;
```
Source: #source(opera + "#L670-L686"), #source(opera + "#L715-L745"). Manual: #source(form-manual).

GammaLoop/idenso implementation currently collects adjacent gamma matrices by replacing two contracted `spenso::gamma` calls with an internal chain:
```rs
AGS.gamma_pattern(RS.a__, RS.b__, RS.c__) * AGS.gamma_pattern(RS.d__, RS.c__, RS.e__),
GS.chain_pattern(RS.b__, RS.e__, [RS.a__, RS.c__, RS.d__]),
```
Source: #source(gamma-rs + "#L604-L613").

Symbolica + spenso pattern:
```rs
gamma(bis(4,a),bis(4,b),p(2,mink(4)))
* gamma(bis(4,b),bis(4,c),mink(4,mu))
  -> chain(bis(4,a),bis(4,c),
       gamma(in,out,p(2,mink(4))),
       gamma(in,out,mink(4,mu)))
```
Assumptions: `in` and `out` are local slots of each factor in the chain. Dummy bispinor labels consumed during collection are not retained as free indices.

== G2. Chain orientation and reversal <rule-gamma-orientation>

Mathematical identity:
$ (Γ^(v_1 ... v_n))_α^β = ((Γ^(v_n ... v_1))_β^α)^"rev" $ <eq-gamma-orientation>

Here `rev` means that the implementation may represent the same contracted tensor network with opposite local slot orientation. This is not a gamma-algebra identity by itself; it is a tensor-pattern canonicalization convention.

Pattern examples:
```rs
A(1,j,la,i,mu) * B(j,k) * C(k,l)
  -> chain(rep,i,l,A(1,out,la,in,mu),B(in,out),C(in,out))

B(j,k) * C(l,k)
  -> chain(rep,j,l,B(in,out),C(out,in))
```

Source context: Symbolica matches syntactically with wildcard variables and respects function-argument ordering for non-symmetric functions; see #source(sym-pattern). spenso provides symbolic tensor structures and contraction support; see #source(spenso-doc).

== G3. Metric contraction in gamma chains <rule-gamma-metric-contract>

Dimension-generic terminal identity:
$ γ^μ_α^β γ_(μ,β)^χ = D δ_α^χ $ <eq-gamma-metric-contract>

In the approved chain notation:
```rs
chain(bis(D,a),bis(D,c),
  gamma(in,out,mink(D,mu)),
  gamma(in,out,mink(D,mu)))
  -> D * g(bis(D,a),bis(D,c))
```

FORM source context: `Trace4Gen` tests repeated adjacent objects and puts the repeated pair into the accumulated delta list before recursing.
```c
if ( *p == p[1] ) {
  *(t->accup)++ = *p;
  *(t->accup)++ = *p;
  ...
  Trace4Gen(..., number-2)
}
```
Source: #source(opera + "#L958-L972"). Manual: #source(form-manual).

GammaLoop/idenso records the dimension contraction in its repeated-Lorentz chain patterns:
```rs
function!(GS.gamma_chain, ..., mink(d_,a_), ..., mink(d_,a_), ...)
  -> function!(GS.gamma_chain, ...) * d_
```
Source: #source(gamma-rs + "#L784-L829").

Assumptions: $D$ is the dimension attached to the Lorentz representation. Four-dimensional Chisholm rules below require $D = 4$.

== G4. Trace closure with explicit bispinor contraction <rule-gamma-trace-closure>

Mathematical identity:
$ (Γ^(v_1 ... v_n))_α^α = γ(v_1)_α^β ... γ(v_n)_χ^α $ <eq-gamma-trace-closure>

Pattern:
```rs
chain(bis(4,a),bis(4,a),
  gamma(in,out,p(2,mink(4))),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,p(3,mink(4))))
  -> trace(bis(4),
       gamma(in,out,p(2,mink(4))),
       gamma(in,out,mink(4,mu)),
       gamma(in,out,p(3,mink(4))))
```

FORM source context: FORM `trace4` and `tracen` operate on a spin line, represented internally as a string of gamma matrices. `Trace4` copies that string into `t->inlist` and then dispatches to `Trace4Gen`.
```c
t->inlist = AT.WorkPointer;
...
for ( i = 0; i < number; i++ ) *p++ = *m++;
...
ret = Trace4Gen(...);
```
Source: #source(opera + "#L730-L745"), #source(opera + "#L820-L823"). Manual: #source(form-manual).

GammaLoop/idenso closes chains when the endpoints agree:
```rs
replace(function!(GS.gamma_chain, RS.a__, RS.x_, RS.x_).to_pattern())
  .repeat()
  .with(function!(GS.gamma_trace, RS.a__).to_pattern())
```
Source: #source(gamma-rs + "#L877-L879").

== G5. Odd and even trace recursion <rule-gamma-trace-recursion>

Odd trace without gamma-five:
$ (Γ^(v_1 ... v_(2 k + 1)))_α^α = 0 $ <eq-gamma-odd-trace>

Two-gamma trace:
$ γ^μ_α^β γ^ν_β^α = 4 η^(μ ν) $ <eq-gamma-two-trace>

Recursive even trace:
$
  (Γ^(v_1 ... v_(2 k)))_α^α
  = sum_(r=2)^(2 k) (-1)^r η^(v_1 v_r)
  (Γ^(v_2 ... hat(v_r) ... v_(2 k)))_β^β
$ <eq-gamma-even-trace>

FORM source context:
```c
if ( ( number < 0 ) || ( number & 1 ) ) return(0);
...
if ( number == 2 ) { *kron++ = *t->inlist; *kron++ = t->inlist[1]; }
```
Source: #source(opera + "#L424-L432"). The recursive trace generator comments describe zero- and two-gamma terminal cases before general recursion: #source(opera + "#L830-L856"). Manual: #source(form-manual).

Symbolica + spenso pattern:
```rs
trace(bis(4),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)))
  -> 4 * g(mink(4,mu),mink(4,nu))
```
Dummy-index freshness: recursive trace replacement removes the paired Lorentz argument and must not reuse bound labels introduced by metric contractions.

== G6. Four-dimensional Chisholm reductions <rule-gamma-chisholm>

Odd interior:
$
  γ^μ_α^β (Γ^(ν_1 ... ν_(2 k + 1)))_β^χ γ_(μ,χ)^δ
  = -2 (Γ^(ν_(2 k + 1) ... ν_1))_α^δ
$ <eq-gamma-chisholm-odd>

Even interior:
$
  γ^μ_α^β (Γ^(ν_1 ... ν_(2 k)))_β^χ γ_(μ,χ)^δ
  = 2 (Γ^(ν_(2 k) ν_1 ... ν_(2 k - 1)))_α^δ
  + 2 (Γ^(ν_(2 k - 1) ... ν_1 ν_(2 k)))_α^δ
$ <eq-gamma-chisholm-even>

Special two-interior case:
$
  γ^μ_α^β γ^ν_β^χ γ^ρ_χ^λ γ_(μ,λ)^δ
  = 4 η^(ν ρ) δ_α^δ
$ <eq-gamma-chisholm-two>

FORM source comment:
```c
g(mu)*g(a1)*...*g(an)*g(mu)=
n=odd:  -2*g(an)*...*g(a1)
n=even: 2*g(an)*g(a1)*...*g(a(n-1))
    +2*g(a(n-1))*...*g(a1)*g(an)
There is a special case for n=2 : 4*d(a1,a2)*gi
```
Source: #source(opera + "#L670-L680"). Implementation odd contraction: #source(opera + "#L978-L1012"). Implementation even contraction: #source(opera + "#L1021-L1096"), #source(opera + "#L1118-L1157"). Manual: #source(form-manual).

Pattern:
```rs
chain(bis(4,a),bis(4,d),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu1)),
  gamma(in,out,mink(4,nu2)),
  gamma(in,out,mink(4,nu3)),
  gamma(in,out,mink(4,mu)))
  -> -2 * chain(bis(4,a),bis(4,d),
       gamma(in,out,mink(4,nu3)),
       gamma(in,out,mink(4,nu2)),
       gamma(in,out,mink(4,nu1)))
```
Assumptions: The contracted Lorentz label must belong to a four-dimensional Lorentz representation. The source code checks dimension equality with `4` before applying the Chisholm branches.

== G7. Gamma-five and epsilon rules <rule-gamma-five>

FORM `trace4` uses the four-dimensional identity:
$
  γ^μ γ^ν γ^ρ
  = ε^(μ ν ρ σ) γ_5 γ_σ
  + η^(μ ν) γ^ρ - η^(μ ρ) γ^ν + η^(ν ρ) γ^μ
$ <eq-gamma-trick>

Source comment:
```c
g_(j,mu)*g_(j,nu)*g_(j,ro)=e_(mu,nu,ro,si)*g5_(j)*g_(j,si)
+d_(mu,nu)*g_(j,ro)-d_(mu,ro)*g_(j,nu)+d_(nu,ro)*g_(j,mu)
which is for 4 dimensions only!
```
Source: #source(opera + "#L320-L329"). The gamma-five branch logic in trace generation is at #source(opera + "#L480-L489") and #source(opera + "#L808-L813"). Manual: #source(form-manual).

Pattern:
```rs
chain(bis(4,a),bis(4,b),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)),
  gamma(in,out,mink(4,rho)))
  -> epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))
     * chain(bis(4,a),bis(4,b),
         gamma5(in,out),
         gamma(in,out,mink(4,sigma)))
     + g(mink(4,mu),mink(4,nu)) * chain(... gamma(rho) ...)
     - g(mink(4,mu),mink(4,rho)) * chain(... gamma(nu) ...)
     + g(mink(4,nu),mink(4,rho)) * chain(... gamma(mu) ...)
```
This is explicitly four-dimensional and requires a fresh dummy Lorentz label `sigma`.

= Color-Algebra Replacement Rules from `color.h`

== Overview

#table(
  columns: (2.8cm, 4cm, 6.8cm),
  [Rule], [Convention], [Purpose],
  [Line joining], [`T(i,j,?a)`], [Join open noncommutative color lines.],
  [Trace closure], [`cOlTr`, `cOlTt`], [Convert closed lines to trace blocks and choose trace blocks to simplify.],
  [Casimir reduction], [`cR`, `cA`], [Remove repeated adjacent or separated generators in traces.],
  [Trace normalization], [`I2R`], [Evaluate two-generator trace.],
  [Structure constants], [`f` antisymmetric], [Contract small f loops and rewrite larger loops through invariants.],
  [Symmetric invariants], [`cOldR`, `cOldA`], [Represent and contract generalized d tensors.],
  [`simpli` strategy], [procedure], [Eliminate f tensors in environments of invariants using generalized Jacobi moves.],
)

== C1. Open color-line joining <rule-color-line-join>

Mathematical identity:
$
  (T^(a_1) ... T^(a_m))_i^j (T^(b_1) ... T^(b_n))_j^k
  = (T^(a_1) ... T^(a_m) T^(b_1) ... T^(b_n))_i^k
$ <eq-color-line-join>

FORM source:
```form
repeat id T(cOli1?,cOli2?,?a)*T(cOli2?,cOli3?,?b) = T(cOli1,cOli3,?a,?b);
id  T(cOli1?,cOli1?,?a) = cOlTr(?a);
```
Source package lines: `color.h` lines 91--93 at #source(color-src). Package documentation: #source(color-doc). No official GitHub line anchor was found for `color.h`; see the validation appendix.

Pattern:
```rs
chain(cof(Nc,i),dind(cof(Nc,j)),
  t(coad(Nc^2-1,a1),in,out),
  t(coad(Nc^2-1,a2),in,out))
*
chain(cof(Nc,j),dind(cof(Nc,k)),
  t(coad(Nc^2-1,b1),in,out))
  -> chain(cof(Nc,i),dind(cof(Nc,k)),
       t(coad(Nc^2-1,a1),in,out),
       t(coad(Nc^2-1,a2),in,out),
       t(coad(Nc^2-1,b1),in,out))
```

== C2. Trace closure and trace normalization <rule-color-trace>

Closed line:
$ (T^(a_1) ... T^(a_n))_i^i = op("tr")_R (T^(a_1) ... T^(a_n)) $ <eq-color-trace-closure>

One generator:
$ op("tr")_R (T^a) = 0 $ <eq-color-trace-one>

Two generators:
$ op("tr")_R (T^a T^b) = T_R δ^(a b) $ <eq-color-trace-two>

FORM source:
```form
id  cOlTt(cOli1?) = 0;
id  cOlTt(cOli1?,cOli2?) = I2R*d_(cOli1,cOli2);
id  cOlTt(cOli1?,cOli2?,cOli3?) = cOldR(...) + i_/2*f(...)*I2R;
```
Source package lines: `color.h` lines 459--461 at #source(color-src). Package documentation: #source(color-doc).

Pattern:
```rs
chain(cof(Nc,i),dind(cof(Nc,i)),
  t(coad(Nc^2-1,a),in,out),
  t(coad(Nc^2-1,b),in,out))
  -> trace(cof(Nc),
       t(coad(Nc^2-1,a),in,out),
       t(coad(Nc^2-1,b),in,out))
  -> TR * g(coad(Nc^2-1,a),coad(Nc^2-1,b))
```

== C3. Fundamental Casimir and adjacent contractions <rule-color-casimir>

Mathematical identities:
$ (T^a)_i^j (T^a)_j^k = C_F δ_i^k $ <eq-color-casimir-fund>

$
  op("tr")_R (T^a T^b T^a op("product")(X))
  = (C_F - C_A / 2) op("tr")_R (T^b op("product")(X))
$ <eq-color-casimir-inside-trace>

FORM source:
```form
id  cOlTr(cOli1?,cOli1?,?a) = cR*cOlTr(?a);
id  cOlTr(cOli1?,cOli2?,cOli1?,?a) = [cR-cA/2]*cOlTr(cOli2,?a);
```
Source package lines: `color.h` lines 108--111 at #source(color-src), repeated for `cOlTt` at lines 173--175. Package documentation: #source(color-doc).

Pattern:
```rs
chain(cof(Nc,i),dind(cof(Nc,k)),
  t(coad(Nc^2-1,a),in,out),
  t(coad(Nc^2-1,a),in,out))
  -> C_F * g(cof(Nc,i),dind(cof(Nc,k)))
```
Assumptions: `color.h` keeps `C_F` as `cR`, not automatically as $(N_c^2 - 1) / (2 N_c)$.

== C4. Fundamental Fierz generator contraction <rule-color-fierz>

Mathematical identity:
$
  (T^a)_i^j (T^a)_k^l
  = T_R (δ_i^l δ_k^j - N_c^(-1) δ_i^j δ_k^l)
$ <eq-color-fierz>

This identity is not written in this explicit form in `color.h`, whose main algorithm uses trace joining and invariant reductions. It is implemented directly in GammaLoop/idenso:
```rs
t(e,a,b) * t(e,c,d)
  -> TR * (id(a,d) * id(c,b) - id(a,b) * id(c,d) / Nc)
```
Source: #source(color-rs + "#L407-L414"). Symbolica documentation for replacement mechanics: #source(sym-pattern).

Pattern:
```rs
t(coad(Nc^2-1,a),cof(Nc,i),dind(cof(Nc,j)))
* t(coad(Nc^2-1,a),cof(Nc,k),dind(cof(Nc,l)))
  -> TR * (
       g(cof(Nc,i),dind(cof(Nc,l)))
       * g(cof(Nc,k),dind(cof(Nc,j)))
       - g(cof(Nc,i),dind(cof(Nc,j)))
       * g(cof(Nc,k),dind(cof(Nc,l))) / Nc)
```

== C5. Structure-constant loop contractions <rule-color-ff>

Mathematical identities:
$ f^(a c d) f^(b c d) = C_A δ^(a b) $ <eq-color-ff-two>

$ f^(a b e) f^(b c f) f^(c a g) = C_A / 2 f^(e f g) $ <eq-color-fff-three>

FORM source:
```form
ReplaceLoop,f,a=3,l<4,outfun=cOlff;
id  cOlff(cOli1?,cOli2?) = -cA*d_(cOli1,cOli2);
id  cOlff(cOli1?,cOli2?,cOli3?) = cA/2*f(cOli1,cOli2,cOli3);
```
Source package lines: `color.h` lines 148--150, 208--210, 505--507, and 538--541 at #source(color-src). Package documentation: #source(color-doc).

Pattern:
```rs
f(coad(Nc^2-1,a),coad(Nc^2-1,c),coad(Nc^2-1,d))
* f(coad(Nc^2-1,b),coad(Nc^2-1,c),coad(Nc^2-1,d))
  -> ca * g(coad(Nc^2-1,a),coad(Nc^2-1,b))
```
Sign note: FORM's `ReplaceLoop` loop orientation produces `-cA*d_(...)` for the two-argument cyclic helper. The standard identity above assumes the displayed orientation $f^(a c d) f^(b c d)$.

== C6. Mixed trace--structure contraction <rule-color-trace-f>

Mathematical identity:
$
  op("tr")_R (T^a T^b op("product")(X)) f^(a b c)
  = i C_A / 2 op("tr")_R (T^c op("product")(X))
$ <eq-color-trace-f>

FORM source:
```form
id  cOlTr(cOli1?,cOli2?,?a)*f(cOli1?,cOli2?,cOli3?) =
  i_*cA/2*cOlTr(cOli3,?a);
```
Source package lines: `color.h` lines 108--111 and 173--175 at #source(color-src). Package documentation: #source(color-doc).

Pattern:
```rs
trace(cof(Nc),
  t(coad(Nc^2-1,a),in,out),
  t(coad(Nc^2-1,b),in,out),
  rest__)
* f(coad(Nc^2-1,a),coad(Nc^2-1,b),coad(Nc^2-1,c))
  -> i_ * ca / 2 * trace(cof(Nc),
       t(coad(Nc^2-1,c),in,out),
       rest__)
```

== C7. Symmetric invariant contractions <rule-color-dd>

`color.h` uses `cOldR` and `cOldA` for symmetric invariant tensors and reduces contractions to named invariants such as `d33`, `d44`, and `d55`.

Representative identities:
$ d_R^(a b c) d_R^(a b d) = d_(33)(R,R) δ^(c d) / N_A $ <eq-color-d33-partial>

$ d_X^(a b c) d_Y^(a b c) = d_(33)(X,Y) $ <eq-color-d33-full>

FORM source:
```form
id cOldr1?cOldar[cOlx3](cOli1?,cOli2?,cOli3?)
 *cOldr2?cOldar[cOlx4](cOli1?,cOli2?,cOli4?)
 = d33(...)*cOlx1*d_(cOli3,cOli4)/NA;
...
id cOldr1?cOldar[cOlx3](cOli1?,...,cOli3?)
 *cOldr2?cOldar[cOlx4](cOli1?,...,cOli3?) = d33(...);
```
Source package lines: `color.h` lines 1287--1294 and 1316--1335 at #source(color-src). Package documentation: #source(color-doc).

Pattern:
```rs
d(sym_rep_x,coad(Nc^2-1,a),coad(Nc^2-1,b),coad(Nc^2-1,c))
* d(sym_rep_y,coad(Nc^2-1,a),coad(Nc^2-1,b),coad(Nc^2-1,d))
  -> d33(sym_rep_x,sym_rep_y)
     * g(coad(Nc^2-1,c),coad(Nc^2-1,d)) / (Nc^2 - 1)
```
Assumptions: The package treats these as group-invariant placeholders rather than expanding them into $N_c$ polynomials.

== C8. `simpli` and generalized Jacobi strategy <rule-color-simpli>

`simpli` is an implementation strategy, not one scalar identity. It eliminates pairs of structure constants in environments of symmetric invariants and applies generalized Jacobi transformations.

FORM source:
```form
* Procedure tries to eliminate f tensors in an environment of invariants.
* Apply the generalized Jacobi identity ...
if ( count(f,1) != multipleof(2) ) discard;
#call contract
...
id f(cOlpR`isimpli',cOli1?,cOli3?)*f(cOlpR`isimpli',cOli2?,cOli3?) =
  cOlpR`isimpli'(cOli1)*cOlpR`isimpli'(cOli2)*cA/2;
```
Source package lines: `color.h` lines 688--713 at #source(color-src). Package documentation: #source(color-doc).

Pattern family:
```rs
f(coad(Nc^2-1,x),coad(Nc^2-1,a),coad(Nc^2-1,c))
* f(coad(Nc^2-1,x),coad(Nc^2-1,b),coad(Nc^2-1,c))
* invariant_environment__
  -> ca / 2 * projected_invariant(a,b,invariant_environment__)
```
This rule class requires canonicalization of symmetric invariants and freshness of all auxiliary adjoint labels introduced during generalized Jacobi expansion.

= Symbolica + spenso Rewrite Patterns

This section is written for a Rust implementation with Symbolica as a dependency. Symbolica's Rust `id` module exposes pattern matching through `AtomCore::replace`, `AtomCore::pattern_match`, `Replacement`, `MatchSettings`, and `ReplaceBuilder`; complex right-hand sides can be built with `ReplaceBuilder::with_map` or with `replace_map` callbacks. See #source(sym-rust-id) and #source(sym-rust-replace). spenso represents tensor slots through symbolic structures and representation slots; the crate documents tensor structures, symbolic tensors, and contraction support at #source(spenso-doc), with symbolic tensor-library keys in #source(spenso-symbolic).

Symbolica wildcard convention: `x_` is one atom, `x__` is one-or-more atoms, and `x___` is zero-or-more atoms. Function argument order is preserved for non-symmetric functions. Repeated wildcard names require exact structural equality; this is the right mechanism for repeated dummy labels such as `mu_` or `a_`.

The public target syntax for this specification is `chain(...)` and `trace(...)`. If implemented inside GammaLoop/idenso, the current internal Rust symbols are `spenso::gamma_chain` and `spenso::gamma_trace`; those should be adapter names only, not the public pattern vocabulary. The approved surface form remains:
```rs
chain(bis(4,a),bis(4,c),
  gamma(in,out,p(2,mink(4))),
  gamma(in,out,mink(4,mu)))
```

== Rust Pattern Construction Guidelines

Use concrete spenso builders for tensor slots where available and Symbolica wildcards only for the abstract labels. This keeps representation constraints attached to the pattern.

```rs
use symbolica::{
  atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
  function, symbol,
  id::{Match, MatchSettings, Pattern, Replacement},
  utils::Settable,
};
```

Implementation rules:
- Use `Replacement::new(lhs.to_pattern(), rhs.to_pattern())` for local structural rules whose right side is another static pattern.
- Use `replace_multiple_into` in a loop for FORM-style short-circuit passes where rule order matters and each pass should be followed by expansion and metric simplification.
- Use `replace_map` or `ReplaceBuilder::with_map` for Chisholm parity, sequence reversal, trace recursion, `ReplaceLoop`-style loop detection, cyclic color traces, and fresh dummy allocation.
- Add `MatchSettings { rhs_cache_size: 0, ..Default::default() }` for replacement maps whose right side depends on match-local fresh symbols.
- Keep `chain` and `trace` non-symmetric. Their argument order is algebraically meaningful.

== Gamma Collection Patterns

Mathematical identity: @eq-gamma-chain-product. Source: #source(gamma-rs + "#L604-L613").

Raw gamma pair:
```rs
gamma(bis(d_,a_),bis(d_,b_),x_)
* gamma(bis(d_,b_),bis(d_,c_),y_)
  -> chain(bis(d_,a_),bis(d_,c_),
       gamma(in,out,x_),
       gamma(in,out,y_))
```

Append and prepend raw gamma factors:
```rs
chain(bis(d_,a_),bis(d_,b_),xs___)
* gamma(bis(d_,b_),bis(d_,c_),y_)
  -> chain(bis(d_,a_),bis(d_,c_),xs___,gamma(in,out,y_))

gamma(bis(d_,a_),bis(d_,b_),x_)
* chain(bis(d_,b_),bis(d_,c_),ys___)
  -> chain(bis(d_,a_),bis(d_,c_),gamma(in,out,x_),ys___)
```

Join two chains:
```rs
chain(bis(d_,a_),bis(d_,b_),xs___)
* chain(bis(d_,b_),bis(d_,c_),ys___)
  -> chain(bis(d_,a_),bis(d_,c_),xs___,ys___)
```

Orientation case:
```rs
chain(bis(d_,a_),bis(d_,b_),xs___)
* chain(bis(d_,c_),bis(d_,b_),ys___)
  -> chain(bis(d_,a_),bis(d_,c_),xs___,reverse_flip(ys___))
```

`reverse_flip` is not a plain Symbolica replacement. Build it in Rust by reading the matched `ys___` argument list, reversing the entries, and swapping `in` and `out` in each tensor factor.

== Gamma Terminal and Trace Patterns

Mathematical identity: @eq-gamma-trace-closure.

```rs
chain(bis(d_,a_),bis(d_,a_),xs___)
  -> trace(bis(d_),xs___)

trace(bis(d_),gamma(in,out,x_))
  -> 0

trace(bis(4),
  gamma(in,out,mink(4,mu_)),
  gamma(in,out,mink(4,nu_)))
  -> 4 * g(mink(4,mu_),mink(4,nu_))
```

Adjacent repeated Lorentz object:
```rs
chain(bis(d_,a_),bis(d_,b_),
  xs___,
  gamma(in,out,mink(d_,mu_)),
  gamma(in,out,mink(d_,mu_)),
  ys___)
  -> d_ * chain(bis(d_,a_),bis(d_,b_),xs___,ys___)
```

GammaLoop/idenso currently has the equivalent repeated-object contraction in Rust at #source(gamma-rs + "#L784-L829") and trace closure at #source(gamma-rs + "#L877-L879"). FORM's corresponding short-circuit in `Trace4Gen` removes adjacent equal objects before deeper recursion; source #source(opera + "#L958-L972").

== Gamma Chisholm and Trace Recursion Builders

Mathematical identities: @eq-gamma-chisholm-odd, @eq-gamma-chisholm-even, and @eq-gamma-even-trace. Source: #source(opera + "#L670-L680"), #source(opera + "#L978-L1012"), #source(opera + "#L1021-L1096"), #source(opera + "#L1118-L1157").

Match contracted endpoints in four dimensions:
```rs
chain(bis(4,a_),bis(4,b_),
  gamma(in,out,mink(4,mu_)),
  middle___,
  gamma(in,out,mink(4,mu_)))
```

Rust-side builder:
```rs
fn build_chisholm(middle: &[Atom], a: Atom, b: Atom) -> Option<Atom> {
  if !all_gamma_factors(middle) {
    return None;
  }
  match middle.len() {
    n if n % 2 == 1 => {
      Some(-Atom::num(2) * chain(a, b, reverse(middle)))
    }
    2 => {
      let (x, y) = lorentz_pair(&middle[0], &middle[1])?;
      Some(Atom::num(4) * g(x, y) * bis_id(a, b))
    }
    n if n % 2 == 0 => {
      Some(Atom::num(2) * chain(a.clone(), b.clone(), move_last_to_front(middle))
        + Atom::num(2) * chain(a, b, reverse_except_last_then_last(middle)))
    }
    _ => None,
  }
}
```

Trace recursion should also be a Rust-side builder, not a static RHS:
```rs
trace(bis(4),first_,rest___)
```

Builder rule: if the trace length is odd, return zero. If the length is two, return `4 * g(first, second)`. Otherwise sum over all pairings of `first_` with one later gamma, alternating signs and recursively rebuilding `trace(bis(4), ...)`. This mirrors the source branch that uses zero/two-gamma terminals and then recursive generation; see #source(opera + "#L424-L432") and #source(opera + "#L830-L856").

== Color Collection Patterns

Mathematical identity: @eq-color-line-join.

Raw generator pair:
```rs
t(coad(na_,a_),cof(nc_,i_),dind(cof(nc_,j_)))
* t(coad(na_,b_),cof(nc_,j_),dind(cof(nc_,k_)))
  -> chain(cof(nc_,i_),dind(cof(nc_,k_)),
       t(coad(na_,a_),in,out),
       t(coad(na_,b_),in,out))
```

Join chains:
```rs
chain(cof(nc_,i_),dind(cof(nc_,j_)),xs___)
* chain(cof(nc_,j_),dind(cof(nc_,k_)),ys___)
  -> chain(cof(nc_,i_),dind(cof(nc_,k_)),xs___,ys___)
```

Opposite endpoint orientation:
```rs
chain(cof(nc_,i_),dind(cof(nc_,j_)),xs___)
* chain(cof(nc_,k_),dind(cof(nc_,j_)),ys___)
  -> chain(cof(nc_,i_),dind(cof(nc_,k_)),xs___,reverse_flip(ys___))
```

The first two rules correspond to `color.h` line joining and trace closure at package lines 91--93. GammaLoop/idenso's current direct color simplifier starts from raw `spenso::t` and `spenso::f` patterns; source #source(color-rs + "#L360-L414").

== Color Trace and Casimir Patterns

Mathematical identity: @eq-color-trace-closure and @eq-color-trace-two.

```rs
chain(cof(nc_,i_),dind(cof(nc_,i_)),xs___)
  -> trace(cof(nc_),xs___)

trace(cof(nc_),t(coad(na_,a_),in,out))
  -> 0

trace(cof(nc_),
  t(coad(na_,a_),in,out),
  t(coad(na_,b_),in,out))
  -> TR * g(coad(na_,a_),coad(na_,b_))

trace(cof(nc_),
  t(coad(na_,a_),in,out),
  t(coad(na_,b_),in,out),
  t(coad(na_,c_),in,out))
  -> dR(coad(na_,a_),coad(na_,b_),coad(na_,c_))
     + i_ / 2 * TR * f(coad(na_,a_),coad(na_,b_),coad(na_,c_))
```

Casimir and separated-generator shortcuts:
```rs
chain(cof(nc_,i_),dind(cof(nc_,k_)),
  xs___,
  t(coad(na_,a_),in,out),
  t(coad(na_,a_),in,out),
  ys___)
  -> CF * chain(cof(nc_,i_),dind(cof(nc_,k_)),xs___,ys___)

trace(cof(nc_),
  t(coad(na_,a_),in,out),
  t(coad(na_,b_),in,out),
  t(coad(na_,a_),in,out),
  xs___)
  -> (CF - CA / 2) * trace(cof(nc_),t(coad(na_,b_),in,out),xs___)
```

Source: `color.h` selected-trace terminal rules at lines 459--464 and easy contraction rules at lines 108--111 and 173--175.

== Color Structure-Constant and Invariant Patterns

Structure-constant contractions:
```rs
f(coad(na_,a_),coad(na_,c_),coad(na_,d_))
* f(coad(na_,b_),coad(na_,c_),coad(na_,d_))
  -> CA * g(coad(na_,a_),coad(na_,b_))

f(coad(na_,a_),coad(na_,b_),coad(na_,e_))
* f(coad(na_,b_),coad(na_,c_),coad(na_,f_))
* f(coad(na_,c_),coad(na_,a_),coad(na_,g_))
  -> CA / 2 * f(coad(na_,e_),coad(na_,f_),coad(na_,g_))
```

Trace--structure shortcut:
```rs
trace(cof(nc_),
  t(coad(na_,a_),in,out),
  t(coad(na_,b_),in,out),
  xs___)
* f(coad(na_,a_),coad(na_,b_),coad(na_,c_))
  -> i_ * CA / 2 * trace(cof(nc_),t(coad(na_,c_),in,out),xs___)
```

Symmetric-invariant contraction:
```rs
d(sym_x_,coad(na_,a_),coad(na_,b_),coad(na_,c_))
* d(sym_y_,coad(na_,a_),coad(na_,b_),coad(na_,d_))
  -> d33(sym_x_,sym_y_) * g(coad(na_,c_),coad(na_,d_)) / na_
```

`f` must be registered or canonicalized as antisymmetric. Symbolica's documentation warns that antisymmetric functions canonicalize written patterns; implementers should either use canonical argument order in the pattern or absorb the sign in a pre-canonicalization pass. Source for small `f` loop reductions: `color.h` lines 148--150, 208--210, 505--507, and 537--541. Source for `d33` and higher invariant contractions: `color.h` lines 1287--1348.

== FORM Short-Circuit Rewrite Order

These are the high-priority short-circuit patterns FORM uses before falling back to expensive generic recursion.

Gamma `trace4` short-circuit order:
- Reject impossible input: negative number, odd number, or missing spin-line parameter returns zero. Source #source(opera + "#L715-L717").
- Resolve zero- and two-gamma traces before recursion; gamma-five can force zero. Source #source(opera + "#L424-L432") and #source(opera + "#L867-L875").
- Remove adjacent equal objects before Chisholm scans. Source #source(opera + "#L958-L972").
- Try four-dimensional odd Chisholm contraction first. Source #source(opera + "#L978-L1012").
- Try four-dimensional even Chisholm contraction, with a special two-interior case before the longer case. Source #source(opera + "#L1021-L1096") and #source(opera + "#L1118-L1157").
- Search same objects by shortest cyclic distance and emit alternating recursive terms. Source #source(opera + "#L1168-L1232").
- If all objects are different, fall back to `Trace4no`; the three-gamma epsilon trick is part of the four-dimensional branch. Source #source(opera + "#L320-L329") and #source(opera + "#L1240-L1244").

Color `color.h` short-circuit order:
- Join open color lines and close equal endpoints before any trace selection. Source `color.h` lines 91--93.
- Skip the main trace route unless at least one closed trace exists. Source `color.h` lines 95--103.
- Apply easy trace contractions first: adjacent equal generator, separated `a b a`, and trace times `f(a,b,c)`. Source `color.h` lines 108--111 and 173--175.
- Rewrite duplicated trace blocks of lengths four, three, and two before general expansion. Source `color.h` lines 112--147 and 178--207.
- Collapse small `f` loops via `ReplaceLoop`: two-leg and three-leg loops are terminal. Source `color.h` lines 148--150, 208--210, 505--507, 537--541, and 672--676.
- Select the trace with the fewest different indices by wrapping traces in `cOlTT`. Source `color.h` lines 156--163.
- Use terminal selected-trace rules for one, two, three, and four generators before converting longer traces back to `T`. Source `color.h` lines 459--466.
- In `simpli`, discard odd `f` counts, contract matching invariant-vector `f f` pairs, then run generalized Jacobi/invariant passes. Source `color.h` lines 688--724.
- Convert completed invariant-vector products into `d33`, `d44`, `d55`, and higher invariants, and discard impossible overlap topologies. Source `color.h` lines 1258--1294 and 1316--1348.

Rust implementation note: implement the short-circuit order as ordered passes, not as one unordered `replace_multiple` containing every rule. FORM's performance relies on reducing the expression before the broader searches are attempted.

= FORM-Side Test Cases in Symbolica/spenso Syntax

This section translates the compact FORM-side rule tests and the package example cases into the `chain`/`trace` notation used above. Expected values are written in the same symbolic convention: `TR` for `I2R`, `CA` for `cA`, `CF` for `cR`, and `NA` for adjoint dimension. The local validation used FORM 4.3.1 for the explicitly listed compact values.

== Gamma Trace Tests

Source basis: FORM manual gamma algebra and `trace4` implementation; see #source(form-manual), #source(opera + "#L424-L432"), and #source(opera + "#L867-L875").

```rs
// FORM: g_(1,mu,nu); trace4,1;
trace(bis(4),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)))
  -> 4 * g(mink(4,mu),mink(4,nu))

// FORM: g_(1,mu,nu,rho); trace4,1;
trace(bis(4),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)),
  gamma(in,out,mink(4,rho)))
  -> 0

// FORM: g_(1,mu,nu,rho,sigma); trace4,1;
trace(bis(4),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)),
  gamma(in,out,mink(4,rho)),
  gamma(in,out,mink(4,sigma)))
  -> 4 * g(mink(4,mu),mink(4,nu)) * g(mink(4,rho),mink(4,sigma))
   - 4 * g(mink(4,mu),mink(4,rho)) * g(mink(4,nu),mink(4,sigma))
   + 4 * g(mink(4,mu),mink(4,sigma)) * g(mink(4,nu),mink(4,rho))
```

Chisholm and adjacent-contraction chain tests:
```rs
// FORM manual: g_(1,mu,mu) = gi_(1)*d_(mu,mu).
chain(bis(d,a),bis(d,b),
  gamma(in,out,mink(d,mu)),
  gamma(in,out,mink(d,mu)))
  -> d * g(bis(d,a),bis(d,b))

// FORM manual, four-dimensional odd interior.
chain(bis(4,a),bis(4,b),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu1)),
  gamma(in,out,mink(4,nu2)),
  gamma(in,out,mink(4,nu3)),
  gamma(in,out,mink(4,mu)))
  -> -2 * chain(bis(4,a),bis(4,b),
       gamma(in,out,mink(4,nu3)),
       gamma(in,out,mink(4,nu2)),
       gamma(in,out,mink(4,nu1)))

// FORM manual, four-dimensional two-interior special case.
chain(bis(4,a),bis(4,b),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)),
  gamma(in,out,mink(4,rho)),
  gamma(in,out,mink(4,mu)))
  -> 4 * g(mink(4,nu),mink(4,rho)) * g(bis(4,a),bis(4,b))
```

Gamma-five epsilon test:
```rs
// FORM Trick identity in opera.c.
chain(bis(4,a),bis(4,b),
  gamma(in,out,mink(4,mu)),
  gamma(in,out,mink(4,nu)),
  gamma(in,out,mink(4,rho)))
  -> epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))
     * chain(bis(4,a),bis(4,b),
         gamma5(in,out),
         gamma(in,out,mink(4,sigma)))
     + g(mink(4,mu),mink(4,nu))
       * chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,rho)))
     - g(mink(4,mu),mink(4,rho))
       * chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu)))
     + g(mink(4,nu),mink(4,rho))
       * chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))
```

FORM manual gamma-five stress examples:
```rs
// FORM: symmetric trace of gamma5 and 12 regular matrices via distrib_ and tracen.
g5_trace_symmetric_12(m1,...,m12)
  -> 51975 terms

// FORM: regular trace g_(1,5_,m1,...,m12); trace4,1.
trace(bis(4),gamma5(in,out),gamma(in,out,mink(4,m1)),...,gamma(in,out,mink(4,m12)))
  -> 1029 output terms after FORM trace4
```
These are term-count regression tests from the manual, not compact algebraic values.

== Compact Color Rule Tests

Source basis: `color.h` line joining, selected trace terminals, and small `f` loops; see #source(color-src), especially lines 91--93, 459--464, and 505--507.

```rs
// FORM: T(i1,i1,a)
chain(cof(Nc,i),dind(cof(Nc,i)),
  t(coad(NA,a),in,out))
  -> 0

// FORM: T(i1,i1,a,b)
chain(cof(Nc,i),dind(cof(Nc,i)),
  t(coad(NA,a),in,out),
  t(coad(NA,b),in,out))
  -> TR * g(coad(NA,a),coad(NA,b))

// FORM: T(i1,i1,a,b,c)
chain(cof(Nc,i),dind(cof(Nc,i)),
  t(coad(NA,a),in,out),
  t(coad(NA,b),in,out),
  t(coad(NA,c),in,out))
  -> dR(coad(NA,a),coad(NA,b),coad(NA,c))
     + i_ / 2 * TR * f(coad(NA,a),coad(NA,b),coad(NA,c))
```

```rs
// FORM: f(a,c,d)*f(b,c,d)
f(coad(NA,a),coad(NA,c),coad(NA,d))
* f(coad(NA,b),coad(NA,c),coad(NA,d))
  -> CA * g(coad(NA,a),coad(NA,b))

// FORM: f(a,b,e)*f(b,c,f)*f(c,a,g)
f(coad(NA,a),coad(NA,b),coad(NA,e))
* f(coad(NA,b),coad(NA,c),coad(NA,f))
* f(coad(NA,c),coad(NA,a),coad(NA,g))
  -> CA / 2 * f(coad(NA,e),coad(NA,f),coad(NA,g))
```

Four-generator trace terminal:
```rs
trace(cof(Nc),
  t(coad(NA,a),in,out),
  t(coad(NA,b),in,out),
  t(coad(NA,c),in,out),
  t(coad(NA,d),in,out))
  -> dR(coad(NA,a),coad(NA,b),coad(NA,c),coad(NA,d))
     + i_ / 2 * (
         dR(coad(NA,a),coad(NA,b),coad(NA,x))
         * f(coad(NA,c),coad(NA,d),coad(NA,x))
         + dR(coad(NA,c),coad(NA,d),coad(NA,x))
         * f(coad(NA,a),coad(NA,b),coad(NA,x)))
     - TR / 6 * f(coad(NA,a),coad(NA,c),coad(NA,x))
       * f(coad(NA,b),coad(NA,d),coad(NA,x))
     + TR / 3 * f(coad(NA,a),coad(NA,d),coad(NA,x))
       * f(coad(NA,b),coad(NA,c),coad(NA,x))
```

== FORM `color.tar.gz` Example Families

The color package page links `color.tar.gz` as "some FORM programs that use color.h"; see #source(color-doc) and #source(color-examples). The example programs define source-side graph families. The following schematics preserve their topology in Symbolica/spenso notation.

`tloop.frm` cases:
```rs
// TYPE = qloop, SIZE = n.
chain(cof(Nc,i1),dind(cof(Nc,i1)),
  t(coad(NA,j1),in,out),...,t(coad(NA,jn),in,out),
  t(coad(NA,j1),in,out),...,t(coad(NA,jn),in,out))

// For SIZE = 3, FORM color+simpli expected:
Q6 -> NA * TR * CF^2 - 3/2 * NA * TR * CA * CF + 1/2 * NA * TR * CA^2

// TYPE = gloop, SIZE = n.
f(coad(NA,k1),coad(NA,k2),coad(NA,j1)) * ...
* f(coad(NA,k(2 n)),coad(NA,k1),coad(NA,jn))

// For SIZE = 3:
G6 -> 0

// TYPE = qqloop, SIZE = n: two fundamental loops linked by adjoint labels.
trace(cof(Nc),t(j1),...,t(jn)) * trace(cof(Nc),t(j1),...,t(jn))

// For SIZE = 3:
QQ3 -> -1/4 * NA * TR^2 * CA + d33(R1,R2)

// TYPE = qgloop, SIZE = n: one fundamental loop and one adjoint loop.
trace(cof(Nc),t(j1),...,t(jn)) * adjoint_loop(j1,...,jn)

// For SIZE = 3:
QG3 -> i_ / 4 * NA * TR * CA^2

// TYPE = ggloop, SIZE = n: two linked adjoint loops.
adjoint_loop(j1,...,jn) * adjoint_loop(j1,...,jn)

// For SIZE = 3:
GG3 -> 1/4 * NA * CA^3
```

Fixed graph examples:
```rs
// TYPE = g14: Coxeter graph with 14 trivalent adjoint vertices.
g14 -> 1/648 * NA * CA^7
    - 8/15 * d444(A1,A2,A3) * CA
    + 16/9 * d644(A1,A2,A3)

// TYPE = fiveq: five linked quark loops.
fiveq -> 1/192 * NA * TR^5 * CA^3
      + 1/4 * d33(R1,R2) * TR^3 * CA^2
      + 5/48 * d33(R1,R2) * d33(R3,R4) * NA^-1 * TR * CA
      + 5/48 * d33(R1,R3) * d33(R2,R4) * NA^-1 * TR * CA
      + 1/8 * d33(R1,R4) * d33(R2,R3) * NA^-1 * TR * CA
      + 1/16 * d44(R1,A1) * TR^4
      + 3/8 * d433(R3,R1,R2) * TR^2 * CA
      + 1/2 * d3333(R1,R2,R3,R4) * TR * CA
      + d43333a(R5,R2,R1,R4,R3)

// TYPE = g24 is a high-cost stress example in tloop.frm.
// No compact expected value is included in the package source; it is not a short unit test.
```

`su.frm` examples using `#call SUn`:
```rs
// CASE = F3F3.
F3F3 -> 2*a^3*NF^-1 - 5/2*a^3*NF + 1/2*a^3*NF^3

// CASE = F4F4.
F4F4 -> -3*a^4*nf^2*NF^-2 + 4*a^4*nf^2
        - 7/6*a^4*nf^2*NF^2 + 1/6*a^4*nf^2*NF^4

// CASE = F4A4.
F4A4 -> -2*a^4*nf*NF + 5/3*a^4*nf*NF^3 + 1/3*a^4*nf*NF^5

// CASE = A4A4.
A4A4 -> -24*a^4*NF^2 + 70/3*a^4*NF^4 + 2/3*a^4*NF^6

// CASE = FnFn with n = 3.
F3F3 -> 2*a^3*nf^2*NF^-1 - 2*a^3*nf^2*NF

// CASE = FnAn with n = 3.
F3A3 -> -i_*a^3*nf*NF^2 + i_*a^3*nf*NF^4

// CASE = AnAn with n = 3.
A3A3 -> -2*a^3*NF^3 + 2*a^3*NF^5

// CASE = qloop with n = 3.
Q6 -> -a^3*nf*NF^-2 + a^3*nf*NF^2

// CASE = gloop with n = 3.
G6 -> 0
```

The SU checker procedure uses $T_R = a$, $C_F = a (N_F^2 - 1) / N_F$, $C_A = 2 a N_F$, and $N_A = N_F^2 - 1$; source `SUn.prc` lines 26--58.

= Implementation Caveats

- FORM's four-dimensional gamma simplifications depend on index dimension checks in `Trace4Gen`; applying them in dimension-generic algebra is wrong.
- `trace4` accepts gamma-five branches with special sign changes and gamma-six/gamma-seven toggles; `tracen` is dimension-generic and does not tolerate gamma-five variants.
- `color.h` is group- and representation-invariant. It intentionally keeps symbols such as `cR`, `cA`, `I2R`, `NA`, and higher `d` invariants rather than reducing all cases to $N_c$.
- The `simpli` procedure is a rewrite strategy over invariant environments. Its output depends on canonicalization choices and on the `RANK` bound.
- Symbolica matching is syntactic. Sequence reversal, parity tests, fresh dummy generation, and antisymmetric signs require Rust-side construction or explicit pattern restrictions.
- The approved `chain` and `trace` syntax in this document is a target specification for tensor-pattern rules, not a claim that Symbolica ships built-in functions with those names.

= Appendix: Source Map

#table(
  columns: (2.4cm, 4.2cm, 3.8cm, 4.4cm),
  [Rule], [Identity], [Equation label], [Primary source],
  [G1],
  [$γ(v)_α^β γ(w)_β^χ = (Γ^(v w))_α^χ$],
  [@eq-gamma-chain-join],
  [#source(opera + "#L715-L745"); #source(gamma-rs + "#L604-L613")],

  [G2], [chain orientation], [@eq-gamma-orientation], [#source(sym-pattern); #source(spenso-doc)],
  [G3],
  [$γ^μ_α^β γ_(μ,β)^χ = D δ_α^χ$],
  [@eq-gamma-metric-contract],
  [#source(opera + "#L958-L972"); #source(gamma-rs + "#L784-L829")],

  [G4],
  [$(Γ^(v_1 ... v_n))_α^α$],
  [@eq-gamma-trace-closure],
  [#source(opera + "#L730-L745"); #source(gamma-rs + "#L877-L879")],

  [G5],
  [odd trace zero, two trace],
  [@eq-gamma-odd-trace, @eq-gamma-two-trace],
  [#source(opera + "#L424-L432"); #source(opera + "#L830-L856")],

  [G6],
  [Chisholm odd/even],
  [@eq-gamma-chisholm-odd, @eq-gamma-chisholm-even],
  [#source(opera + "#L670-L680"); #source(opera + "#L978-L1012"); #source(opera + "#L1118-L1157")],

  [G7],
  [gamma-five/equation epsilon trick],
  [@eq-gamma-trick],
  [#source(opera + "#L320-L329"); #source(opera + "#L480-L489")],

  [C1],
  [open color-line joining],
  [@eq-color-line-join],
  [`color.h` lines 91--93, #source(color-src); #source(color-doc)],

  [C2],
  [trace closure and $T_R$],
  [@eq-color-trace-closure, @eq-color-trace-two],
  [`color.h` lines 459--461, #source(color-src)],

  [C3],
  [Casimir inside trace],
  [@eq-color-casimir-fund, @eq-color-casimir-inside-trace],
  [`color.h` lines 108--111, #source(color-src)],

  [C4], [fundamental Fierz], [@eq-color-fierz], [#source(color-rs + "#L407-L414")],
  [C5],
  [$f f$ and $f f f$],
  [@eq-color-ff-two, @eq-color-fff-three],
  [`color.h` lines 148--150 and 505--507, #source(color-src)],

  [C6], [trace--$f$ contraction], [@eq-color-trace-f], [`color.h` lines 108--111 and 173--175, #source(color-src)],
  [C7],
  [`d33` contractions],
  [@eq-color-d33-partial, @eq-color-d33-full],
  [`color.h` lines 1287--1294 and 1316--1335, #source(color-src)],

  [C8], [`simpli` strategy], [`<rule-color-simpli>`], [`color.h` lines 688--713, #source(color-src)],
  [S1], [Rust Symbolica pattern API], [Symbolica + spenso section], [#source(sym-rust-id); #source(sym-rust-replace)],
  [S2],
  [spenso symbolic tensor structures],
  [Symbolica + spenso section],
  [#source(spenso-doc); #source(spenso-symbolic)],

  [S3],
  [FORM gamma short-circuit order],
  [FORM short-circuit section],
  [#source(opera + "#L424-L432"); #source(opera + "#L958-L972"); #source(opera + "#L1168-L1232")],

  [S4],
  [FORM color short-circuit order],
  [FORM short-circuit section],
  [`color.h` lines 91--163, 459--507, 688--724, and 1258--1348, #source(color-src)],

  [T1], [FORM gamma test values], [FORM-side test section], [#source(form-manual); #source(opera + "#L424-L432")],
  [T2],
  [FORM compact color test values],
  [FORM-side test section],
  [`color.h` lines 459--464 and 505--507, #source(color-src)],

  [T3],
  [FORM color example families],
  [FORM-side test section],
  [`color.tar.gz`: `tloop.frm`, `su.frm`, `SUn.prc`, #source(color-examples)],
)

Source snippets appear next to each detailed rule above. Symbolica/spenso pattern references are the code blocks in the corresponding detailed rule and in the Symbolica + spenso rewrite-pattern section.

= Appendix: Validation

== Typst compile result

Typst version:
```sh
typst 0.14.2 (b33de9de)
```

Compilation command attempted:
```sh
typst compile form_gamma_color_rules_from_scratch.typ form_gamma_color_rules_from_scratch.pdf
```

Result: succeeded. The PDF file `form_gamma_color_rules_from_scratch.pdf` was produced.

== Static checks

Checks run against this Typst source:
```sh
rg -n '\b\x70rod\b' form_gamma_color_rules_from_scratch.typ
rg -n 'tr_[A-Za-z]+\(|Tr_[A-Za-z]+\(|op\("tr"\)_[A-Za-z]+\(' form_gamma_color_rules_from_scratch.typ
rg -n '#L[0-9]+' form_gamma_color_rules_from_scratch.typ
rg -n 'gamma_chain|gamma_trace' form_gamma_color_rules_from_scratch.typ
rg -n 'with_map|replace_map|FORM Short-Circuit' form_gamma_color_rules_from_scratch.typ
rg -n 'FORM-Side Test Cases|Q6 ->|F3F3 ->|g14 ->|fiveq ->' form_gamma_color_rules_from_scratch.typ
```

Results:
- No word token `\x70rod` was found; intended occurrences use `product`.
- No subscripted trace operator was immediately attached to a parenthesized argument; trace calls use explicit spacing.
- GitHub line-anchor fragments were found for the FORM `opera.c`, GammaLoop `gamma.rs`, and GammaLoop `color.rs` source links. `color.h` rules cite exact package-source line numbers because no official GitHub-hosted `color.h` was found in `form-dev/form`.
- The internal names `spenso::gamma_chain` and `spenso::gamma_trace` occur only in source-code snippets or explanatory text documenting GammaLoop internals; target patterns use `chain` and `trace`.
- Rust-only operations that cannot be static Symbolica RHS patterns, including sequence reversal, parity checks, cyclic loop detection, and fresh-index creation, are assigned to `with_map`, `replace_map`, or explicit Rust builders.
- The FORM-side test-case section is present and contains expected values for compact gamma tests, compact color tests, and short `color.tar.gz` example instantiations.

Additional FORM checks used to populate the test-case section:
```sh
form -q /tmp/gamma-test.frm
form -q /tmp/color-trace.frm
form -q /tmp/color-short.frm
cd /tmp/form-color && form -q tloop-qloop.frm
cd /tmp/form-color && form -q tloop-gloop.frm
cd /tmp/form-color && form -q tloop-qqloop.frm
cd /tmp/form-color && form -q tloop-qgloop.frm
cd /tmp/form-color && form -q tloop-ggloop.frm
cd /tmp/form-color && form -q tloop-g14.frm
cd /tmp/form-color && form -q tloop-fiveq.frm
cd /tmp/form-color && form -q su-F3F3.frm
cd /tmp/form-color && form -q su-F4F4.frm
cd /tmp/form-color && form -q su-F4A4.frm
cd /tmp/form-color && form -q su-A4A4.frm
cd /tmp/form-color && form -q su-FnFn.frm
cd /tmp/form-color && form -q su-FnAn.frm
cd /tmp/form-color && form -q su-AnAn.frm
cd /tmp/form-color && form -q su-qloop.frm
cd /tmp/form-color && form -q su-gloop.frm
```

== Unresolved uncertainties

- The current official FORM GitHub repository contains `sources/opera.c` but not the legacy `color.h` package file. The package page links the actual source. Therefore color rules cite exact package-source line numbers and the package source URL, but they do not have official GitHub line anchors.
- The approved `chain` and `trace` notation is a specification-level Symbolica/spenso tensor-pattern notation. Exact Rust constructors for sequence reversal and parity-filtered Chisholm replacement must be implemented with Rust-side match processing.

= References

- FORM manual gamma algebra: #source(form-manual).
- FORM `opera.c`: #source(opera).
- FORM `color.h` package page: #source(color-doc).
- FORM `color.h` source: #source(color-src).
- FORM `color.h` example programs: #source(color-examples).
- GammaLoop/idenso gamma implementation: #source(gamma-rs).
- GammaLoop/idenso color implementation: #source(color-rs).
- Symbolica pattern matching: #source(sym-pattern).
- Symbolica Rust `id` API: #source(sym-rust-id).
- Symbolica Rust `ReplaceBuilder`: #source(sym-rust-replace).
- spenso crate documentation: #source(spenso-doc).
- spenso symbolic tensor-library source: #source(spenso-symbolic).
- Typst syntax: #source(typst-syntax).
- Typst math: #source(typst-math).
