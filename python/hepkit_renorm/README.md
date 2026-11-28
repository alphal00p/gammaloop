# hepkit_renorm

Algebraic renormalisation of gauge theories with a non-anticommuting gamma5 (BMHV scheme, sterile partners,
finite Slavnov-Taylor-restoring counterterms), built on the Symbolica / HepKIT Python API. Input: a declarative
model specification (gauge group, matter, classical four-dimensional Lagrangian, gauge fixing, regularisation).
Target output: a JSON model with the renormalised theory. Physics reference: the report *One-loop renormalisation
of the Standard Model with a non-anticommuting gamma5* (sections 02_method, 03_setup, 03b_operator_bases).

Status (milestones MS0 and MS1): specification, groups, local functionals, the regularised classical action,
BRST / Slavnov-Taylor (`brst`), Feynman rules (`rules`), exact fingerprints (`fingerprint`) and exact linear
algebra (`linalg`). Gates V1 (tree-level Feynman rules) and V2 (Slavnov-Taylor identities of `S_0`) pass, also in
the independent review tests `tests/test_independent_*.py`. Pure Python; no Rust change is needed.

## Layout

```
python/hepkit_renorm/
  pyproject.toml              package metadata, pytest and ruff settings
  hepkit_renorm/
    exact.py                  prime fields GF(p), square roots, rational reconstruction, exact linear algebra
    tensors.py                index labels, tensor heads, BMHV simplifier, conjugation, canonical relabelling
    groups.py                 SU(N) x U(1): explicit generators, structure constants, invariants, invariant tensors
    fields.py                 field registry (statistics, quantum numbers, slots, conjugation)
    spec.py                   model specification -> validated Model with derived field content
    functional.py             monomials, densities, local functionals, derivations, vertex functions
    action.py                 regularised classical action from a Model
    brst.py                   BRST operator, gauge fermion, S_0, Slavnov-Taylor operator, b_d, b, Delta_1
    rules.py                  vertex rules, rule sets, propagators
    fingerprint.py            exact plane-wave evaluation of functionals in D = 4, 6, 8 over GF(p)
    linalg.py                 exact rank / kernel / solve, rational reconstruction (CRT), interpolation
    models/sm_one_family.toml the one-family Standard Model of the report
  tests/                      pytest suite; tests/data holds frozen golden data (see tests/data/README.md);
                              test_independent_*.py + independent_support.py: the independent review of MS0/MS1
                              (own evaluator: Weyl representation, own generators, 62-bit prime)
```

## Running the tests

```
cd python/hepkit_renorm
SYMBOLICA_LICENSE=... /common/build-cache/vhirschi/gamma5/venv/bin/python -m pytest   # 192 tests, about 2 min
```

The runner is pytest (installed into the venv with `uv pip install --python <venv>/bin/python pytest`; the
repository's own dev group also uses pytest). `SYMBOLICA_LICENSE` is read from the environment only (pass it
inline, never write it to a file). The tests also run in Symbolica's unlicensed mode, which is limited to one
Symbolica process per user: with other Symbolica processes running, an unlicensed run aborts at import. A test
marked `requires_licence` is skipped with a message when the variable is absent. Tests never read files outside this directory. Lint: `uvx ruff check hepkit_renorm tests`.

## Conventions (report 03_setup, 03b; binding)

1. `d = 4 - 2 eps`, metric `diag(+,-,-,-,...)`; every index explicit; a label occurring twice is contracted.
2. `gbar`/`ghat` are the 4- and (d-4)-dimensional metrics; `gammabar`/`gammahat` the split gammas;
   `{gammabar, gamma5} = 0`, `[gammahat, gamma5] = 0`; `ghat(mu,mu) = d - 4` (tree level: no eta tag).
3. `eps^{0123} = +1` (Idenso's symbol is `-i` times it: convert at the boundary).
4. `D_mu = d_mu - i g_r T^a A^a_mu`; `cal(T)^a = g_r T^a`, `fg^{abc} = g_r f^{abc}`; `Tr(T^a T^b) = delta/2`.
5. Grassmann order = order of occurrences; derivations act from the left; left/right derivatives as in eq:ST.
6. Vertex function `V` = coefficient of `theta_1..theta_n` in `Y[sum theta_k e_k e^{-i p_k x}]`, momenta incoming,
   `d_mu -> -i p_mu`; Feynman rule `i V`; exchanging two odd legs flips the sign.
7. Hermitian conjugation reverses the occurrences; `gh` real, `agh^dag = -agh`, the source of a real field is
   imaginary, `Psi^dag -> Psibar` with the Dirac adjoint on coefficients.
8. Fields: `F` (gauge field of factor F), `gh_F`, `agh_F`, `b_F`, fermion `f`/`fbar` (Dirac field incl. the
   sterile partner), scalar `phi`/`phidag`, source of X `src_X` (= X* of eq:Phistar).
9. Slots of dimension one are dropped (U(1) adjoint, singlets, one generation).
10. Exact arithmetic only: Symbolica rationals and `i` symbolically; group tensors over GF(p).

## Public API

Import paths are `hepkit_renorm.<module>`; the names re-exported by `hepkit_renorm` are marked (top).

### spec

- `builtin_model(name) -> Model` (top): load `models/<name>.toml` (e.g. `"sm_one_family"`).
- `load_model(path) -> Model` (top), `model_from_dict(doc) -> Model` (top): TOML format in the `spec` module
  docstring and in `models/sm_one_family.toml`. Unknown keys, unsupported options, non-invariant Yukawa or
  potential terms raise.
- `Model` (top), a dataclass with
  - `name`, `generations`, `factors: dict[str, GaugeFactor]` (ordered), `fermions: dict[str, FermionMultiplet]`,
    `scalars: dict[str, ScalarMultiplet]`, `parameters: dict[str, Parameter]`, `potential: tuple[PotentialTerm]`,
    `yukawas: tuple[YukawaTerm]`, `gauge_fixing: GaugeFixing` (`.kind == "linear_covariant"`, `.xi` an
    Expression), `regularisation: Regularisation`, `fields: FieldTable` (derived);
  - `param(name) -> Expression`: symbol of a parameter or coupling (real unless declared complex);
    `coupling(factor) -> Expression`; `flavour_tensor(name, i, j)`: `y(i, j)` or the plain symbol without
    flavour slots;
  - `simplifier(expand_projectors=True) -> tensors.Simplifier` with the model's dimensions;
    `group_evaluator(field) -> groups.GroupEvaluator`; `slot_labels(field_name) -> tuple[Expression]` (fresh);
  - `yukawa_structure(y) -> (projector, [(factor, "delta" | "epsf", [(field_position, slot), (..)])])`;
  - `multiplet_generators(multiplet, field, couplings=None) -> ndarray[a, i, j]` (combined adjoint index over
    all factors, the report's 12 values for the SM); `group_invariants(field) -> {"S": {(f, r): expr},
    "C": {f: expr}, "C_ad": {r: expr}}` (tab:invariants); `anomaly_tensor(field) -> ndarray[a, b, c]`
    (eq:Dabc, unit couplings, times generations);
  - `rigid_defects(density, field) -> list[Monomial]`: monomials that are not rigidly invariant (empty list:
    invariant; the test is monomial by monomial, which is sufficient);
  - `vanishes_modulo_group(expr, open_labels, field) -> bool`: zero test of a simplified expression (e.g. a
    vertex function) after explicit evaluation of its group tensors, terms collected by their non-group part;
    decides Lie-algebra identities such as `[T^a, T^b] = i f^{abc} T^c` that the simplifier does not apply
    (needed for `s S_0`: its physical four-point parts cancel only by the Lie algebra).
- Data classes: `FermionMultiplet(name, chirality, reps, flavoured=True)`, `ScalarMultiplet(name, complex, reps,
  flavoured=False)`, `Parameter(name, real=True, flavour=False)`, `PotentialTerm(coupling, bilinears,
  factor=1)`, `YukawaTerm(coupling, fermion_bar, scalar, fermion)`, `GaugeFixing(kind, xi)`, `Regularisation`.

### fields

- `FieldTable`: ordered immutable mapping `name -> FieldInfo`; `order(name)` (canonical occurrence order),
  `odd(name)`.
- `FieldInfo(name, kind, odd, ghost_number, dimension, slots, charges, conjugate, conj_sign, factor, multiplet,
  chirality, source_of, source)`; `kind` in `gauge, ghost, antighost, auxiliary, fermion, fermion_bar, scalar,
  scalar_conj, source`; `is_source`; `lorentz_slots()`.
- `Slot(space, factor=None, rep=None)`: `rep` in `"F"`, `"Fbar"`, `"A"` for group slots (already conjugated for
  conjugate fields and sources); U(1) charges are in `FieldInfo.charges`.
- Slot order of a field: Lorentz, spinor, group slots in factor order, flavour. Statistics, ghost number and
  dimension of the sources follow from their field (opposite statistics, `-1 - gh`, `3 - dim`).

### tensors

- Label constructors: `lor(n)`, `spin(n)`, `adj(factor_symbol(F), n)`, `fund(factor_symbol(F), n)`, `flav(n)`;
  `Space(kind, factor, dim)` with `.label(n)`, `.key`; `LORENTZ`, `SPINOR`, `FLAVOUR(n)`, `adjoint_space`,
  `fundamental_space`; `fresh(space) -> label` (process-wide, deterministic, numbers >= 10^6).
- Heads (Symbolica symbols, also in `HEADS`): `g, gbar, ghat, eps, p, sp, spbar, sphat, gamma, gammabar, gammahat,
  gamma5, PL, PR, delta, T, f, dsym, epsf, d`; `I` is the imaginary unit. Meaning in the module docstring.
- `Simplifier(group_dims, n_flavours, expand_projectors=True)`; `__call__(expr) -> expr`: expanded, rewritten to
  fixpoint, dummies relabelled canonically. Canonical mode: projectors as `(1 -/+ gamma5)/2`, full gammas split,
  `gamma5` moved right, free full metrics and `sp` split into barred + hatted parts. Gamma matrices are not
  reordered among themselves: general equality needs exact numerical fingerprints.
- `conjugate(expr)`: complex conjugate with the Dirac adjoint of spinor heads; `absorb(expr)`: absorb `g` and
  `delta`; `canonical_relabel(term, keep=())`, `canonical_form(expr, keep=())`; `labels_in(term)`,
  `free_labels(term)`, `label_space(label)`, `rename(expr, pairs)`, `contains_head(expr, head)`,
  `is_evanescent(term)` (contains `ghat`, `gammahat` or `sphat`).

### groups

- `GaugeFactor(name, group, coupling)`: `group` is `"SU(N)"` or `"U(1)"`; `abelian`, `n`, `dim_adjoint`,
  `rep_dim(rep)`, `adjoint_space()`, `rep_space(rep)`. Representations: `"1"`, `"F"`, `"Fbar"`, `"A"`, or a
  `Fraction` charge for U(1); `conjugate_rep(factor, rep)`.
- Explicit tensors over a `PrimeField`: `su_generators(n, field)` (generalised Gell-Mann / 2, Gell-Mann order),
  `su_structure_constants`, `su_d_tensor`, `rep_generators(factor, rep, field)`, `product_generators(factors,
  reps, field, couplings=None)`.
- Invariants: `dynkin_index(factor, rep, field)`, `casimir_value(factor, rep, field)` (Fractions, exact);
  `invariant_tensors(factor, reps, field) -> list[ndarray]` (null space of `sum_k T^a_(k)`, dense, for products of
  up to a few thousand components).
- Symbolic builders: `generator_expr(factor, rep, a, i, j)` (`T(a,i,j)`, `-T(a,j,i)`, `-i f(a,i,j)` or the
  charge), `structure_expr(factor, a, b, c)`.
- `GroupEvaluator(factors, field, n_flavours)`: `evaluate(term, open_labels, values=None) -> ndarray` for a
  product term of group heads, numbers and symbols (values for the symbols); `number(expr)`;
  `complex_rational(num) -> (re, im)`.

### exact

- `PrimeField(p)`: `field(x)` image of an int/Fraction, `inv`, `sqrt` (smaller root, fixes the embedding of
  `Q(i, sqrt n)`), `i`, `reconstruct(x) -> Fraction`, `array`, `zeros`, `reduce`, `einsum`, `matmul`, `is_zero`,
  `row_reduce`, `rank`, `nullspace`. Elements are Python ints, arrays are `numpy` object arrays.
- `default_field()` (31-bit prime with `sqrt` of -1, 2, 3, 5, 6, 7, 10), `second_field()` (independent prime,
  for confirming reconstructions), `find_prime(start, squares)`, `is_prime`.

### functional

- `Occurrence(field, slots=(), derivatives=())` (top): `d_{derivatives} field[slots]`.
- `Monomial(coefficient, occurrences)` (top): ordered Grassmann product; `fields()`, `separated()`.
- `Density(fields, expr=0)` (top), immutable, canonical, possibly with free labels:
  - constructors `Density.monomial(fields, coefficient, *occurrences)`, `Density.from_monomials(fields, monos)`;
  - views `monomials()`, iteration, `len`, `is_zero()`, `==` (exact canonical equality), `free_labels()`,
    `parity()`, `field_content() -> Counter`;
  - algebra `+`, `-`, unary `-`, `scale(c)`, `*` (graded product, left factor first, equal free labels
    contract), `total_derivative(mu)`, `rename(pairs)`, `map_coefficients(fn)`, `apply(derivation)`,
    `conjugate()`, `simplify(simplifier)`;
  - BMHV: `split_lorentz(simplifier)`, `physical(simplifier)` (P), `evanescent(simplifier)` (E = 1 - P),
    `four_dimensional_limit(simplifier)` (P, then d -> 4).
- `LocalFunctional(fields, expr=0)` (top): a `Density` without free labels, meaning `int d^d x`;
  `LocalFunctional.integrate(density)`, `density()`;
  `functional_derivative(field, labels, side="left" | "right") -> Density` (free `labels` in the slots of the
  field); `integrate_by_parts(field)`; `vertex_function(legs, simplifier=None, momentum_conservation=True) ->
  Expression` (the single code path for V: Feynman rules, counterterm rules and basis vertex functions);
  `leg_classes()`, `vertex_functions(simplifier) -> {class: V}`, `vanishes(simplifier)` (zero up to integration by
  parts). The product of two local functionals is not defined (multiply densities).
- `Leg(field, labels, momentum)` (top), `standard_legs(fields, names) -> list[Leg]` (top): deterministic labels
  `space.label(10000 + 10 k + j)`, momenta `1..n`.
- `Derivation(parity, images)` (top), `FieldImage(labels, density)` (top): a graded left derivation defined by its
  images on fields (missing fields are annihilated, images are checked for parity); `derivation(density)`;
  commutes with partial derivatives.
- `field_symbol(name)`: the symbol naming a field in the canonical expression
  `sum coeff * occ(pos, field, slots..) * dv(pos, mu)` held in `Density.expr`.

### action

- `classical_action(model) -> LocalFunctional` (top): `int (L_gauge + L_fermion + L_scalar - V + L_Yuk)` in the
  BMHV regularisation (eq:Lgauge, eq:Lferm, eq:LHiggs, eq:LYuk); `classical_sectors(model) -> {"gauge",
  "fermion", "scalar", "potential", "yukawa"}` (top).
- Building blocks: `covariant_derivative(model, field, mu, labels) -> Density` (top) for any field, including
  `(D_mu omega)^a` on ghosts; `field_strength(model, factor, mu, nu, a) -> Density` (top); `gauge_kinetic`,
  `fermion_terms`, `scalar_kinetic`, `potential`, `yukawa_terms`.

## Interfaces for the next modules

- `brst` builds `L_gf + L_ext = s(gauge fermion)` (eq:gf-exact), the BRST derivation `s` (`Derivation(1, ...)`
  with images built from `covariant_derivative`, `generator_expr`, `structure_expr` and the slot reps of
  `FieldInfo`), the Slavnov-Taylor operator (products of `functional_derivative`s over the pairs
  `(src_X, X)`), `b_d`, and `Delta_1 = s S_0`. `S_0 = classical_action(model) + L_gf + L_ghost + L_ext`. The
  gauge-fixing parameter is `model.gauge_fixing.xi`; the ghost sector fields and sources are already in
  `model.fields`.
- `rules` wraps `LocalFunctional.vertex_function` (Feynman rule `i V`), inverts the quadratic part for the
  propagators and writes HepKIT JSON; the readable form uses `model.simplifier(expand_projectors=False)`.
- `fingerprint` evaluates coefficients exactly: the group part with `GroupEvaluator`, Lorentz and Dirac heads
  with explicit components in d = 4 and d > 4 (to be added there), over `exact.PrimeField`.

## Design decisions

- Coefficients are Symbolica expressions over the package's own head vocabulary, not spenso/Idenso tensors:
  Idenso's tensor inference rejects metrics between 4- and (d-4)-dimensional slots and Python accepts only
  single-symbol dimensions (scout notes on the algebra), while the functional layer needs explicit
  bar/hat projectors and exact Grassmann bookkeeping. The writer translates to UFO/spenso structures.
- Group tensors are explicit over a prime field: an orthonormal basis of su(3) has no rational form (the trace
  form on the Cartan part has discriminant 3 modulo squares), and SU(N) needs sqrt(2k(k+1)) for every k < N. In
  GF(p) with the required roots every identity is exact; rational invariants are reconstructed.
- Grassmann order lives in the ordered occurrence list, outside Symbolica's commutative product; a density is
  stored as one canonical Symbolica expression, so equal monomials merge and equality is exact.
