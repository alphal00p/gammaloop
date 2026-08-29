# FeynKit tensor reduction

`feynkit-tensor` reduces Lorentz-covariant vacuum tensor numerators to scalar
Spenso invariants. It is designed for the Taylor-expanded vacuum graphs that
enter the R-operation, but the projector is independent of a particular
integral family.

The public boundary follows the rest of the FeynKit ecosystem:

- vectors carry a final `spenso::mink(D,index)` argument;
- scalar contractions are emitted as `spenso::dot`;
- residual free-index metric tensors are emitted as `spenso::g` with complete
  Minkowski slots;
- ordinary Spenso metrics are simplified before the projector is constructed;
- `FeynmanDiagramTensorExt` adds the operation to the native graph type without
  introducing a wrapper graph.

## Why high rank remains tractable

At rank `2k`, a naive Lorentz projector contains `(2k-1)!!` metric
pairings. The Gram matrix of two pairings only depends on the alternating
cycles formed when they are superimposed. Those cycle lengths define a coset
type, or equivalently an integer partition of `k`. Orthogonal Weingarten
coefficients are constant on each type, reducing the coefficient problem to
`p(k)` unknowns:

| rank | labeled pairings | coefficient classes |
| ---: | ---: | ---: |
| 8 | 105 | 5 |
| 12 | 10,395 | 11 |
| 16 | 2,027,025 | 22 |
| 20 | 654,729,075 | 42 |

The coefficient engine implements the Collins--Matsumoto orthogonality
recurrence and caches exact symbolic tables. When all loop momenta are equal,
it uses the closed isotropic coefficient directly,

\[
  \frac{1}{D(D+2)\cdots(D+2k-2)}.
\]

The universal `OrthogonalWeingarten` API retains those `p(k)` coefficient
tables. The fully contracted reducer goes one step further when vectors repeat
on either side. Let `S_A` and `S_B` be the raw zero-or-one incidence matrices
from labeled matchings to internal and projector contraction orbits, let `G` be
the pairing Gram matrix, and let `W = G^-1`. The requested coefficient matrix
is

\[
  C = S_A^T W S_B.
\]

Symmetry gives `G S_A = S_A H_A` on the much smaller internal orbit space (and
similarly for the projector orbit space). FeynKit therefore evaluates the
exact identity

\[
  C = H_A^{-T} N = N H_B^{-1}, \qquad N = S_A^T S_B.
\]

It constructs `N` with a memoized joint-color matching recurrence and `H` from
memoized alternating cycles, then solves whichever side has fewer orbits. No
Cartesian product of labeled pairings is generated. For example, multiplicity
patterns `[2,18]` on both sides of a rank-20 projector require only a `2 x 2`
orbit solve and produce four compact terms.

Each contraction orbit is encoded by exponents `alpha_ij` of
`prod_(i <= j) dot(p_i,p_j)^alpha_ij`, together with the number of labeled
pairings it represents. This preserves the raw-incidence normalization used in
the formulas above.

Completely unsymmetrized free-index output can itself contain `(2k-1)!!`
distinct metric terms; no representation can materialize that cheaply. That
path retains explicit pairing, pairing-product, and output-term budgets so such
requests become a clear error instead of an accidental memory blow-up.

## Rust example

```rust,no_run
use feynkit_tensor::TensorReducer;
use symbolica::parse;

let dimension = parse!("D");
let numerator = parse!("K(1,spenso::mink(D,mu))*K(1,spenso::mink(D,nu))*P(spenso::mink(D,mu))*P(spenso::mink(D,nu))");
let reducer = TensorReducer::new(dimension).with_integrated_head_name("K");
let scalar = reducer.reduce(numerator.as_view())?.into_expression();
println!("{scalar}");
# Ok::<(), feynkit_tensor::TensorReductionError>(())
```

The Python binding mirrors this API and returns one concrete Symbolica
`Expression`:

```python
from symbolica import E
import symbolica.community.feynkit as fk

reducer = fk.TensorReducer.feynkit(E("4"))
scalar_numerator = vacuum_diagram.reduce_tensor_numerator(reducer)
scalar_graphs = vacuum_diagram.reduce_tensor_graphs(reducer)
```

`reduce_tensor_graphs` requires a fully contracted numerator and creates one
deterministically named contribution per compact scalar term while preserving
the topology ID. Use `reduce_tensor_numerator` when residual free Lorentz
indices—and therefore metric-tensor output—must be retained.

The reducer dimension must match the dimension stored in the input
`spenso::mink` slots. Native FeynKit rules currently carry `4`; an
R-operation expression prepared with symbolic `D` should use `E("D")`
instead.

For mixed high-rank tensors, keep the dimension symbolic (for example `D` or
`4-2*eps`) until after reduction. At exceptional fixed integer dimensions the
labeled metric Gram matrix can be singular because dimension-specific tensor
identities make its spanning set redundant. The all-equal isotropic fast path
does not require that inverse and remains well defined at positive physical
dimensions.

For a graph with external `FeynKit::Momentum` tensors, select the internal
compact momenta explicitly with `with_integrated_vector`; the `feynkit`
constructor intentionally selects the whole momentum head and is intended for
pure vacuum numerators.

## Method and validation sources

- B. Collins and S. Matsumoto,
  [On some properties of orthogonal Weingarten functions](https://arxiv.org/abs/0903.5143).
- B. Collins and S. Matsumoto,
  [Weingarten calculus via orthogonality relations](https://arxiv.org/abs/1701.04493).
- F. Herzog *et al.*,
  [The five-loop beta function of Yang--Mills theory with fermions](https://arxiv.org/abs/1801.06084),
  for symmetry-reduced tensor projectors in the R* operation.
- M. Goode *et al.*,
  [Tensor reduction of vacuum integrals](https://arxiv.org/abs/2408.05137),
  for the modern orbit-partition formulation and high-rank projectors.
- J. Davies *et al.*,
  [OPITeR](https://arxiv.org/abs/2411.02233), for integrand-symmetry-aware
  tensor reduction.

Exact coefficients at ranks 2, 4, 6, 8, and 10 are tested against VAKINT's
FORM implementation. At rank 20, the generated table is checked to contain all
42 integer-partition classes. Goode *et al.* publish coefficients in a
transformed symmetric-`dd` orbit basis, so that ancillary table is not a
direct coefficient-by-coefficient oracle for the raw orthogonal-Weingarten
table used internally here.
