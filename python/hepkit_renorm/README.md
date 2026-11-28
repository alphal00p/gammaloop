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
    graphs.py                 graph records (input of the loop back-ends): routing, Wick sign, power counting
    loops.py                  loop back-end interface and the one-loop back-end: N0, N1 (Laurent data) on samples
    store.py                  order store: functionals per loop order and role (tree, breaking, ct_singular,
                              ct_finite), Laurent coefficients in eps
    output.py                 model writer: order store -> HepKIT JSON (today's format + v2 metadata) and
                              the reverse map to package fields, components and vertex functions
    diagrams.py               graphs of a leg class from the store's rule sets (any loop order, counterterm and
                              breaking insertions): enumeration, canonical keys, symmetry factors, graph records
    hepkit_graphs.py          the generator route: HepKIT on the written model, mapped back to multiplet-level graphs
    class_results.py          per-class local results (summed Laurent data of all graphs on fingerprint samples),
                              physical samples, on-disk per-graph cache
    parametric.py             parameter gradings (derived), parametric fingerprints, graded polynomial
                              reconstruction of exact results in the couplings
    matching.py               class data -> singular counterterm (singular basis, renormalisation constants,
                              evanescent part), pole identity, X_l in N, b X = 0, anomaly analysis, bookkeeping-B
                              remainder
    st_solve.py               b DeltaS^(l) = -X_l with model-independent normalisation conditions (N1)-(N4),
                              existence conditions, automatic fallback (N4'), uniqueness as a rank statement
    orders.py                 the order-l step (steps 5-9) and the store entries it leaves for order l + 1;
                              the renormalised model file
    loopdata.py               loop data for the algebraic steps: classes from the graphs, resolved in the
                              parameters, batched and cached; their local form on operator bases
    driver.py                 the order loop from a specification to the renormalised model, with the run report
    __main__.py               command line: python -m hepkit_renorm SPEC
    models/sm_one_family.toml the one-family Standard Model of the report
  tests/                      pytest suite; tests/data holds frozen golden data (see tests/data/README.md);
                              test_independent_*.py + independent_support.py: the independent reviews of MS0/MS1
                              and of MS2-MS7 (own evaluator: Weyl representation, own generators, 62-bit prime;
                              test_independent_v5_v8_run.py: the end-to-end run against the report, slow;
                              test_independent_robustness.py: leg order, wrong signs, an anomalous content);
                              test_review_fixes.py: the fixes of the MS2-MS7 review
```

## Running the tests

```
cd python/hepkit_renorm
SYMBOLICA_LICENSE=... python -m pytest   # 377 passed, 79 slow skipped; 8-12 min
```

The runner is pytest (installed into the venv with `uv pip install --python <venv>/bin/python pytest`; the
repository's own dev group also uses pytest). `SYMBOLICA_LICENSE` is read from the environment only (pass it
inline, never write it to a file). The tests also run in Symbolica's unlicensed mode, which is limited to one
Symbolica process per user: with other Symbolica processes running, an unlicensed run aborts at import. A test
marked `requires_licence` is skipped with a message when the variable is absent. Tests never read files outside this directory. Lint: `uvx ruff check hepkit_renorm tests`.
Long end-to-end tests are marked `slow` and skipped by default; run them with `--runslow` (e.g.
`python -m pytest --runslow tests/test_sm_gates.py`). The whole suite with `--runslow` takes about 66 min on a loaded
host with `HEPKIT_RENORM_PROCESSES=96` and a warm per-graph cache `HEPKIT_RENORM_CACHE=<dir outside the repository>`
(`test_class_results.py::test_all_classes`, which never uses a cache, takes 29 min of it).

## User guide: from a model specification to the renormalised model

### What it does

`python -m hepkit_renorm SPEC` (or `hepkit_renorm.driver.renormalise(SPEC, RunOptions(...))`) renormalises a gauge
theory at one loop in the BMHV scheme with finite Slavnov-Taylor-restoring counterterms, starting from the
classical four-dimensional Lagrangian in the unbroken phase, and writes the renormalised model as a HepKIT JSON file.
Everything between the specification and the file is derived: sterile partners, gauge fixing, ghosts, BRST
sources and transformations, `Delta_1 = s S_0`, the classes of Green's functions, the operator bases, the action of
`b`, the diagrams with their signs and symmetry factors, the loop integrals, the normalisation conditions. The
steps are those of the report's recipe (02_method sec:m-algorithm) at `l = 1`; the order loop and the order store are
those of every order (sec:m-beyond).

### Input

A specification file in TOML (format: the `spec` module docstring; example: `hepkit_renorm/models/sm_one_family.toml`):
gauge factors (`SU(N)` or `U(1)` with a coupling), chiral fermion multiplets with their representations and
hypercharges, scalar multiplets, parameters, the scalar potential, the Yukawa terms, the gauge fixing (linear
covariant; the run uses Feynman gauge) and the regularisation choices. A built-in name (`sm_one_family`) also works.

### Running the one-family Standard Model

```
cd python/hepkit_renorm
SYMBOLICA_LICENSE=... python -m hepkit_renorm sm_one_family \
    -o /path/outside/the/repository/run_sm -j 128
```

Options: `-l/--order` (only 1 is supported; a higher order is refused before any computation with the list of what
is missing), `-j/--processes` (forked workers for the graphs), `--cache DIR` (per-graph cache, default
`<output>/graph_cache`; a rerun with a warm cache evaluates no graph), `--no-cache`, `--normalisation full|minimal`,
`--fermion-kinetic` (the representative with fermion kinetic restoring counterterms), `--no-fallback`,
`--no-pole-identity` (skips the most expensive loop pass, the check of the ghost-number-1 poles in `D = 6`),
`--cohomology`, `--schema v1|v2`, `--seed`, `--name`, `-q`. Exit status 0: finished and every check passed; 1: a check
failed; 2: stopped with an error (anomaly, data outside a basis, unsupported order).

Cost (one-family SM, 128 worker processes on a shared, fully loaded host, no cache): about 25 min, of which the
pole identity (ghost-number-1 poles in `D = 6`) takes 17 min, `P N1` in `D = 4` (both embeddings of `i` and the
confirmation at a second prime) 5 min and the ghost-number-0 poles 1 min; with a warm cache the run takes a few
minutes (the algebra). `run_report.json` records every step's time. The per-graph cache of the SM takes about
50 MB.

### Outputs (in the output directory)

| file | content |
|---|---|
| `<name>_tree.json` | tree-level regularised model (HepKIT JSON, Feynman gauge): `S_0` with ghosts and sources, `Delta_1` (order `DELTA`) |
| `<name>_renormalised.json` | the renormalised model: `S_0`, `Delta_1`, `deltaS1` (MS-bar poles, invariant and evanescent; orders `CT = 1`, `CTPOLE`), `DeltaS1` (finite restoring counterterm; `CT = 1`, `CTFIN`), `bd_deltaS1`, `bd_DeltaS1` (breaking insertions for order 2). Encoding: `output` module docstring; the reverse map to package fields is embedded under `hepkit_renorm.reverse_map` |
| `<name>_store.json` | the order store (exact; `OrderStore.from_json(model, text)`), the input of the next order |
| `<name>_order_data.json` | per order the matched loop data as exact local functionals: `C` (`deltaS^(l) = (1/(16 pi^2))^l C/eps`) and `X_l` |
| `run_report.json`, `run_report.txt` | the run report (below) |

### Reading the run report

`run_report.txt` is a summary; `run_report.json` holds everything:

- `status` (`done` or `failed` with `error`), `all_checks_passed`;
- `steps`: one entry per step (`action`, `tree_model`, `classes`, `bases`, `graphs`, `poles`, `pole_identity`,
  `breaking`, `restoring`, `store`, `model`) with its order `l`, wall time in seconds and what it computed (counts,
  loop evaluations, cache hits);
- `checks`: every check with `passed` and the numbers behind it: `ST(S_0) = Delta_1`; `Delta_1` evanescent; the models
  load in HepKIT; `Q` and `N` complete; the poles lie in the span of the singular basis (exact fit on every class,
  coefficients polynomials on their graded supports, reconstruction confirmed at a further prime); the
  four-dimensional limit of `deltaS` is `b`-closed and spanned by the field and parameter directions; the pole identity
  on every ghost-number-1 class; `P N1` lies in the span of `N`; `b X = 0`; no anomaly; the restoring counterterm is
  unique and the normalisation conditions were applied without relaxation;
- `results["<l>"]`: `constants` (`Z[<field>]` = `delta Z`, `d[<parameter>]` = `delta p`, units
  `(1/(16 pi^2))^l / eps`; bare `Phi_0 = Z^(1/2) Phi`, `p_0 = mu^(w_p eps)(p + delta p)`), `beta` (`beta_p = 2 l delta
  p^(l)`, units `(1/(16 pi^2))^l`), `gamma` (`gamma_Phi = -l delta Z_Phi^(l)`), `evanescent` (coefficients of the
  evanescent operators, with whether `b_d` annihilates them; `evanescent_operators` gives their origin),
  `singular_coefficients` (the singular counterterm on the generated singular basis), `breaking` (`X_l` on the
  generated basis `N`, anomaly components), `restoring` (non-zero `r_alpha` on the generated basis `Q`), `existence`
  (the existence conditions with the classes in which each can be tested), `normalisation` (the conditions with their
  applied and relaxed equations), `graphs_per_class`, `loop_statistics`.

Names of basis elements are generated (`Q<k>[class]`, `N<k>[class]`, `E<k>[class]`) and deterministic; they are not
the report's names. Comparisons with the report go by functionals (fingerprints), see `tests/driver_golden.py`.

### What is supported

- One loop (`l = 1`), any gauge group that is a product of `SU(N)` and `U(1)` factors, chiral fermions and scalars
  as in the specification format, one generation; Feynman gauge; MS-bar; BMHV with sterile partners; the report's
  normalisation conditions (N1)-(N4) stated for any model, with the automatic fallback (N4').
- Exact results: every decision over `GF(p)` at random points with rational reconstruction; symbolic coefficients in
  all free parameters.
- Anomalies are detected and reported (the run stops at the restoring counterterm with `AnomalyError`).

### What is not supported (yet)

- Loop order 2 and higher (`driver.missing_for_order(2)`): an `L >= 2` loop back-end, the subtraction of
  subdivergences, the readout of graphs with counterterm insertions, the `S_q` pieces of `Delta^[l]`, the bookkeeping-B
  remainder on every class.
- Several generations, flavour matrices, complex parameters (CKM; stage D); gauges other than Feynman gauge in the
  graphs; four-fermion vertices.
- The v2 model schema is written as extra keys (`--schema v2` writes it natively, but today's HepKIT loader drops it).

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

### store

- `Role` (`tree`, `breaking`, `ct_singular`, `ct_finite`); `StoreEntry(name, loop_order, role, laurent={k: LocalFunctional},
  description)`: the functional `sum_k eps^k C_k` (coefficients eps-free, may contain `d` and `EVANESCENT_TAG`);
  `powers`, `coefficient(k)`, `leg_classes()`, `to_dict()` / `from_dict(fields, doc)`.
- `OrderStore(fields)`: `add(entry)`, `entry(name)`, `entries(max_loop_order=None, roles=None)`, `loop_orders()`,
  `tree_action()` (sum of the tree entries, the quadratic form of the propagators), `to_json()` /
  `from_json(model, text)` (exact round trip). `tree_store(model, brst=None)`: entries `S0` and `Delta1`.

### output

- `write_model(model, store, options=WriterOptions()) -> WrittenModel`; `tree_model(model, options=None, brst=None)`:
  the tree-level regularised model (`S_0` with gauge fixing, ghosts and sources; `Delta_1` with order `DELTA`) in
  Feynman gauge. Encoding conventions (components for non-native groups, commuting ghost encoding, markers
  `bmhv::bar`/`bmhv::hat`, coupling orders, v2 metadata keys and `v2_overrides`) in the module docstring.
- `WriterOptions(schema="v1" | "v2", name, native_colour="auto", substitutions={"xi": 1}, max_loop_order, orders=
  OrderNames(), regulator="eps", tag="eta", symbols, dimensions=(6, 8, 4), samples, seed)`; `OrderNames` (defaults
  `G<factor>`, `YUK`, `POT`, `DELTA`, `SRC`, `CT`, `CTPOLE`, `CTFIN`).
- `WrittenModel`: `data` (JSON document), `reverse_map`, `json_string()`, `write(path)`, `load()` (HepKIT model),
  `package_rule(vertex, power=0)` (the package `VertexRule` behind a written vertex, legs in its order).
- `ReverseMap`: `particles` (`ParticleRecord`: field, `{slot: component}`, antiname, pdg), `vertices`
  (`VertexRecord`: store entry, loop order, role, package leg class, components per leg, orders, Laurent powers,
  evanescent), `propagators`, `parameters` (UFO name -> package kind and reference), `orders`; `particle(name)`,
  `vertex(name)`, `symbol(model, ufo_name)`, `to_dict()` / `from_dict()` (also embedded in the JSON under
  `hepkit_renorm.reverse_map`). `component_samples(fingerprinter, fields, components, n)`: fingerprint samples with
  the component slots fixed (evaluates a package vertex function at the components of a written vertex).

### graphs

- `GraphRecord(legs, vertices, propagators, symmetry_factor=1, name="", tags=())`: one graph of `Gamma_{legs}`
  (value with Feynman rules `i V`, propagators, symmetry factor and Wick sign); `legs` are the class legs
  (`standard_legs`, Grassmann order, incoming momenta). Derived: `loops`, `routing()` (propagator momenta from momentum
  conservation, loop momenta `p(LOOP_MOMENTUM + l)` on the chords), `leg_momentum(v, j)`, `grassmann_sign(table)`
  (global Wick rule, module docstring), `superficial_degree(table)`, `is_one_particle_irreducible()`, `attachment(end)`.
- `GraphVertex(legs, value, ends, outgoing=frozenset(), dimension=4, order=0, origin="")`, `GraphVertex.from_rule(rule,
  ends, dimension=..)`: any vertex function in the package vocabulary (S_0, Delta_1 with `dimension=5`, counterterms);
  `ends[j]` is `ExternalEnd(leg)` or `PropagatorEnd(propagator, side)` (side 0: the `X` end, momentum flows in).
- `GraphPropagator(fields, labels, numerator, denominators)`: `numerator / prod (k^2 - m^2)^n`, `k = p(0, mu)`;
  `GraphPropagator.from_rule(rules.Propagator, {xi: 1})` splits a rule (Feynman gauge by substitution);
  `split_denominators(value)`.

### loops

- `OneLoopBackend()`: `supports(graph)`, `prepare(graph, model) -> PreparedGraph` (symbolic, dimension-independent),
  `evaluate(prepared, fingerprinter, samples) -> Laurent`. `local_uv(graph, model, fp, samples)`; `evaluate_graphs(
  prepared_list, fp, {leg fields: Samples}, processes=N)` (forked workers, identical results).
- `Laurent(loops, field, n, coefficients)`: `T G = (i/(16 pi^2))^L sum c[r, k] eta^r / eps^k`, `c[r, k]` as values on
  the samples (contracted with the polarisations like `Fingerprinter.evaluate`, any `D`); `N0 = c[0,1]`, `N1 =
  c[1,1]`; pole operators `K()`, `Kf()`, `Ke()` (sec:m-K); `+`, `scale`. Compare `sum_G N0` with the fingerprint of
  `-16 pi^2 eps V_{delta S}` and `P N1` (D = 4) of the ghost-number-one graphs with `X_1`.
- `LoopBackend` (Protocol) for later l-loop back-ends; `divergent_subgraphs(graph, table)` (forest elements);
  `embed_samples(samples, fp, fp_larger)`.
- Method (module docstring): massless Taylor expansion to the superficial degree (logarithmic terms only), tensor
  reduction by pairings of loop-momentum slots, dense Lorentz-Dirac tensors per factor part with separate group
  tensors, exact network contraction; the tag `eta` from the polynomial dependence on the internal hatted dimension,
  obtained by evaluating with internal indices in `D_ext + m` dimensions (odd ones included) and solving the
  Vandermonde system. Validated against the toolkit (`tests/data/loops_toolkit_golden.json`) and the class results
  g0_VV, g1_cV, g1_cFF, g1_cVVV_cVVVV (`tests/test_loops.py`, about 20 s).

### diagrams

- `VertexFamily(name, rules, dimension=4, loop_order=0, power=0, insertion=False)`: the rules of one Laurent
  coefficient of one store entry; `store_families(model, store, order)` (entries of loop order `< order`; tree:
  `tree_rules`, other roles: all classes, two-point ones included; dimension `4 + g dim s`; role `breaking` is an
  insertion); `graph_propagators(model, store, substitutions=None)` (`{(X, Y): GraphPropagator}`, Feynman gauge).
- `enumerate_graphs(model, legs, families, propagators, *, order=1, insertions=0, min_loops=1, loops=None,
  one_particle_irreducible=True) -> Enumeration`: connected (1PI) graphs with `L + sum_v L_v = order`, exactly
  `insertions` insertion vertices; legs in the given (Grassmann) order. `Enumeration.topologies` (kept),
  `.vanishing` (odd automorphism), `.records(model, propagators)`.
- `Topology(legs, vertices, externals, lines, key, symmetry_factor, loops, order, power)`, `.describe(families)`;
  `record(topology, families, model, propagators) -> GraphRecord`; `topology_key(vertices, externals, lines)`
  (canonical key, equal iff isomorphic); `automorphisms(topology)`; `line_keys(propagators)` (one orientation per
  line), `line_fields(a, b, keys)`.

### hepkit_graphs

- `generated_keys(written, model, fields, families, propagators, *, loops=1, order=None, insertions=0,
  tuples=None) -> {key: [Generated(particles, lines, symmetry_factor)]}`: HepKIT's one- or l-loop 1PI diagrams on all
  component tuples (`component_tuples(written, fields)`) mapped to multiplet-level keys through the reverse map
  (workarounds of the generation scout: vertex bound from Euler, exact DELTA/SRC/CT orders, no internal sources).

### class_results

- `ClassComputer(model, store, order=1, *, backend=None, cache=None, processes=1, families=None,
  propagators=None)`: `enumerate(fields)`, `records(fields)`, `samples(fp, fields, n, salt=0, physical=False)`,
  `compute(fields, fp, samples) -> ClassResult`; insertions = ghost number of the class (0 or 1).
- `ClassResult(fields, ghost_number, order, samples, fingerprinter, graphs, vanishing, seconds, cached)`:
  `laurent` (sum; `by_loops()` at `l > 1`), `N0`, `N1`, `compare(values, r=0, k=1)`; `GraphResult(name,
  symmetry_factor, sign, power, laurent)`.
- `physical_samples(fp, samples)` (`Pphys` on the fermion legs); `ResultCache(directory)` (per-graph JSON files keyed
  by SHA-256 of the graph record, back-end, prime, dimension, `i` embedding, parameter values and samples).
- Gates (`tests/test_class_results.py`): V5 (pole, D = 6, full spinors: `V_G`, `V_{b_d G}`) and V6 (`P N1`, D = 4:
  `V_{X_1}`) for every class of ghost number 0 and 1 against the reference engines (`--runslow`; a selection covering
  every sector by default).

### parametric

- `ParameterSpace(model, action, fixed=None)`: free parameters (`names`, `symbols`), fixed ones (rational values,
  e.g. `{"xi": 1}` or every parameter for a computation at one point), the derived gradings `weights` (every
  monomial of `S_0` with `n` legs has weight `n - 2`) and `dimensions`; `expand(F) -> {exponents: F_e}`,
  `grade(F) -> (grade, excess)`, `support(weight, dimension)`, `values(point)`, `random_points`, `expression(poly)`.
- `ParametricMatrix(field, shape, terms)`: `at(point)`, `hstack`, `vstack`, `rows`, `columns`, `scale`,
  `from_values`; `parametric_rows(fp, space, functionals, blocks)`: fingerprints as polynomials in the parameters.
- `reconstruct(make_evaluator, supports, space, complex_values=False)`: exact polynomials from evaluations at random
  points over several primes (CRT, confirmed by one more prime). `Polynomial` = `{exponents: Fraction | (re, im)}`.

### matching

- Class data (units of `class_results`): `ClassData` protocol (`values(quantity, fields, fp, samples)`), quantities
  `POLE` (`N0`) and `BREAKING` (`P N1`); `FunctionalData(brst, counterterm=C, breaking=X, counterterm_variation=None)`
  (golden data or symbolic results; also parametric), `LoopData(ClassComputer)` (values at one point: fix every
  parameter), `CombinedData(default, {class: provider})`.
- `OrderContext(model, order=1, brst=None, fixed=None (Feynman gauge), Q=None, N=None, singular=None, seed=0)`: bases
  (built on first use), fingerprinters per prime, sample blocks, supports.
- `match_singular(ctx, data, dimensions=(6,)) -> SingularResult`: `coefficients` on the singular basis,
  `four_dimensional`, `constants` (`Z[<multiplet or factor>]` = `delta Z`, `d[<parameter>]` = `delta p`; report
  convention, without `1/(16 pi^2 eps)`), `relative(p)`, `evanescent` (`{name: (poly, b_d-closed)}`),
  `closure_rank`, `functional()`. `field_directions(model, action)`, `parameter_directions(model, action)`.
- `pole_identity(ctx, singular, data) -> PoleIdentity(classes, failed, ...)`, `.holds`.
- `match_breaking(ctx, data, cohomology=False, reconstruct_coordinates=False) -> BreakingResult` (`coordinates`,
  `closed`, `anomaly: AnomalyReport(complement, coefficients, image_rank, cohomology)`, `.anomalous`).
- `match_remainder(ctx, data, restoring_functional, classes) -> RemainderResult` (eq:m-LP, bookkeeping B);
  `evanescent_basis(ctx, classes)`.

### st_solve

- `standard_normalisation(model, mode="full" | "minimal", fallback=True, fermion_kinetic=False) -> Normalisation`:
  (N1) no pure-matter classes, (N2) longitudinal two-point gauge terms, (N3) no `rhotilde.d omega` and no
  structure-constant part of `(rhotilde, A, omega)`, (N4) no abelian ghost; `fermion_kinetic` keeps the fermion
  kinetic terms (N1k: no gauge vertex of the first factor of each multiplet).
- Conditions: `ClassCondition(name, description, predicate)`, `ElementCondition(name, description, names)`,
  `ProjectedCondition`, `transverse_two_point(model, factor)`, `antisymmetric_vertex(model, factor)`.
- `solve_restoring(ctx, data, normalisation=None, basis=None) -> RestoringResult`: `coefficients`, `functional()`,
  `nonzero()`, `unique`, `unknowns`, `allowed`, `rank_allowed`, `rank_all`, `kernel`, `existence` (`[(direction,
  classes in which it can be tested)]`), `existence_total`, `conditions` (`ConditionReport`: applied, relaxed,
  dropped). Raises `AnomalyError` (with the `AnomalyReport`) or `ExistenceError` (fallback off).

### orders

- `run_order(ctx, poles, breaking, normalisation=None) -> OrderResult(order, singular, restoring, breaking, poles,
  remainder)`; `order_entries(result)` / `update_store(store, result)`: `deltaS<l>` (`ct_singular`, `eps^-1`),
  `DeltaS<l>` (`ct_finite`), `bd_deltaS<l>`, `bd_DeltaS<l>` (`breaking`, bookkeeping A), `hatDelta<l>` (`breaking`,
  bookkeeping B, with a remainder), loop factor `(1/(16 pi^2))^l` written out; `prune(model, F)` (drop vanishing
  classes); `counterterm_model(model, store)` (the renormalised model file, Feynman gauge).
- Store entries carry `bookkeeping` (`""`, `"A"`, `"B"`); `OrderStore.entries(..., bookkeeping=...)` filters.
- Gates: `tests/test_sm_gates.py` (V7, V8 on the report's golden data, `--runslow`, about 10 min),
  `tests/test_sm_pipeline.py` (pipeline data through the same interface at a rational point, `--runslow`); fast
  toy-model tests `test_parametric.py`, `test_matching.py`, `test_st_solve.py`, `test_orders.py`.

### loopdata

- `ResolvedLoopData(computer, space, processes=1, cache=None, seed=0)`: the `ClassData` protocol from the graphs of a
  `ClassComputer` (order 1): `values(quantity, fields, fp, samples)` at the point of `fp`, `parametric(...)` as a
  polynomial in the free parameters (each graph evaluated at as many points as its derived support has monomials:
  `graph_support(graph, space, omega)`), `prefetch([Request])` (all graphs, points and classes in one pool of
  forked workers, per-graph cache), `graphs(fields)`, `stats` (`LoopStatistics`).
- `RequestRecorder` / `recorded(consumer, space)`: dry run of a consumer of class data (its requests, answered with
  zeros), so that they can be evaluated in one batch.
- `localise(ctx, data, quantity, elements, classes, *, dimension, physical=False, complex_values=False, salt=51)
  -> LocalForm`: the data as a local functional on the elements, coefficients exact polynomials on their graded
  supports, reconstructed at one prime (both embeddings of `i` if complex) and confirmed at a further prime
  (`MatchingError` with the failing classes outside the span). `LocalForm.functional(ctx)`, `.nonzero()`.
- `loops.evaluate_jobs([(prepared, fingerprinter, samples)], processes=N)`: the batch evaluator behind
  `evaluate_graphs`; `parametric.interpolate_outputs` (one prime of `reconstruct`); `linalg.reconstruct_all`.

### driver

- `renormalise(spec, RunOptions(...)) -> RunResult(model, brst, store, orders, contexts, written, report,
  loop_data)`; `RunOptions(output, order=1, cache=None, processes, normalisation="full", fermion_kinetic=False,
  fallback=True, pole_identity=True, cohomology=False, schema="v1", seed=0, name=None, log=None)`;
  `RunReport` (`steps`, `checks`, `results`, `outputs`, `status`, `passed`, `to_dict()`, `text()`); `Check`, `Step`;
  `missing_for_order(l)`, `OrderNotSupported`; `resolve_model(spec)` (model, `.toml` file or built-in name);
  `ORDER_DATA`. Command line: `python -m hepkit_renorm --help`.
- Gates: `tests/test_driver.py` (chiral QED end to end, fast; `--runslow`: its pole identity and an anomalous model),
  `tests/test_driver_sm.py` (the one-family SM from its specification compared with the report coefficient by
  coefficient through `tests/driver_golden.py`, `--runslow`).

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
