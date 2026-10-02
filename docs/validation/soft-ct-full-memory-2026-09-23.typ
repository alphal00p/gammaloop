= Why full U/H exceed 12 GiB after the main rebase

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

A subsequent #link("soft-ct-main-comparison-2026-09-23.typ")[matched main control]
also exceeds 12 GiB on full ordinary U. The allocation diagnosis below therefore
does not establish a new memory regression caused by the soft-counterterm stack.

A further #link("soft-ct-root-sharing-2026-09-23.typ")[structural-size audit]
measures a 75.40-MiB root before protection and demonstrates substantial exact
subtree repetition. Retained numerator-body sizes quoted below do not measure
that assembled root.

Both full direct-3D runs exhaust the memory bound inside the global
`collect_factors()` call in `FinalIntegrands::into_integrands`, after the signed
forest sum has been assembled and before evaluator preprocessing starts. Six
live debugger samples, taken at approximately 3, 6 and 10 GiB for each
prescription, identify the same insertion in Symbolica 3.0.0's
`collect_factors_impl`. This is temporary factor-map growth, including deep
copies of owned expression keys. It is not evidence that a 12-GiB expanded
physical numerator has been constructed.

No production source, test, assertion, algebraic setting or CI setting was
changed for this investigation. The tested source is the prior validation
change `2c2694036e93`; the investigation is a separate Jujutsu change. The
source/method manifest
and binary fingerprint
identify the exact input implementation.

== Comparison and exclusions

Shared input is the three-loop `paper_figure_b1_double_triangle_soft_ir` graph:
102 native CFF source maps, the same physical numerator/projector, graph LMB and
kinematics, and disabled integrated counterterms and threshold subtraction.
Full U/H use explicit orientation sums. The orientation-58 control selects one
localized source pattern; both the source inventory and eventual selector
handling differ. Numerical equivalence of these unequal inventories is not
assumed.

The shared pipeline is CFF/residue maps, classified forest, component numerator
families, Taylor projection, branch normalization, source and forest sums,
global factor collection, finite tensor preprocessing, and evaluator
construction. U/H first differ at Taylor: U(X), versus S(X)+U(X-S(X)) for
positive-degree H. Degree-zero H uses U. Bare removes the forest.

The previous uninstrumented
#link("soft-ct-main-2026-09-23.typ")[measurement report] supplies the comparison:

#table(
  columns: (2fr, 1fr, 1fr, 1.5fr),
  table.header([Case], [Wall (s)], [Peak RSS (GiB)], [Outcome]),
  [Full bare], [1.806], [0.722], [Completed],
  [Orientation 58 U], [7.420], [0.890], [Completed],
  [Orientation 58 H], [33.486], [1.482], [Completed],
  [Full U], [132.806], [12.266], [12 GiB stop],
  [Full H], [448.838], [12.139], [12 GiB stop],
)

Full source inventory alone is insufficient to explain failure: bare completes.
A soft-only mechanism is insufficient: ordinary U also fails. The live stacks
now further exclude branch normalization and numerical evaluator construction
as the active allocation boundary for these stops. The failing owner is the
completed forest's scalar root entering global factor collection. The exact
number of top-level summands and distinct factor keys at this point was not
extracted from the live heap; no exact byte accounting is claimed.

== Live evidence and allocation mechanism

The six-sample stack summary
records the measured RSS and relevant source frames. Full all-thread stacks,
resource traces, commands and compact generation logs accompany it. The active
path is:

```text
AmplitudeGraph::build_integrands
  -> Forests::orientation_parametric_exprs
  -> ParametricIntegrands::from_final
  -> FinalIntegrands::into_integrands
  -> protected.collect_factors()                final_integrand.rs:112
  -> AtomView::collect_factors_impl             Symbolica collect.rs:1343
  -> HashMap insertion / owned Atom clone / Vec<u8> copy
```

The dependency source is
`/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/symbolica-3.0.0/src/collect.rs`.
Its hash is preserved in the manifest. In the Add case, lines 1308-1345:

+ Build a factor map for every summand.
+ Build the union of all factor keys and their candidate common exponents.
+ For every union key, visit every summand map.
+ Subtract the common exponent from present factors; insert a cloned key into
  every map where it is absent, even when that exponent is zero.

The decisive operation is `f.insert(ff.clone(), -*p)`. For the many factors
that are not common to every summand, `p` is zero. An absent factor is already
implicitly raised to zero, but the implementation stores that entry explicitly.
Only the outer common-factor map is pruned immediately afterwards; the per-term
maps retain the inserted zero entries until reconstruction.

For T summands and F distinct factor keys, this can create O(T F) temporary map
entries. When F grows with T, storage is quadratic in the number of terms.
`AtomOrView` distinguishes borrowed views from owned atoms; cloning an owned
key copies its expression buffer. The 6/10-GiB U samples and all H samples show
these byte-buffer clones at the insertion itself. The first U sample catches
hash-map growth and rehashing at the same source line.

This call is inside a fixed-point loop in GammaLoop. The evidence does not show
an infinite loop or identify the iteration number. Repetition is not needed to
reproduce the growth: the independent example below uses one call.

The factor pass and its opaque-wrapper rules are also present in
`main@origin` at `4c84f3cdc68a`. These samples identify the current failure, not
a claim that Symbolica 3 introduced the algorithm or that the soft prescription
introduced the factor pass.

== Independent reproducer

The #link("soft-ct-full-memory-2026-09-23/collect-reproducer.rs")[standalone Rust probe]
links the already-built Symbolica 3 library and forms this factored sum:

```text
E_n = sum_i c(i) * (a(i) + b(i))
```

One `collect_factors()` call returns an expression structurally identical to
its input in every run. There is no CFF, Taylor operation, tensor contraction,
explicit expansion or fixed-point loop. Fresh-process
measurements are:

#table(
  columns: (1fr, 1.4fr, 1.4fr, 1fr),
  table.header([Terms], [Input/output bytes], [Peak RSS (MiB)], [Time (s)]),
  [128], [5,253], [6.000], [0.0064],
  [256], [10,502], [15.012], [0.0265],
  [512], [21,766], [57.008], [0.1416],
  [1024], [44,294], [219.020], [0.7051],
)

The near-fourfold RSS increase on doubling the larger inputs supports the
factor-map mechanism seen in the actual runs. These synthetic sizes are not
measurements of the QCD expression. They isolate the collector's allocation
behavior without changing or expanding a graph numerator.

== Expression walkthrough and expansion audit

+ #strong[Construct and retain physical factors.] The graph numerator builder
  multiplies the selected edge and vertex numerator factors. This is ordinary
  symbolic product construction, not a whole-numerator `.expand()`. The CFF
  supplies scalar denominator expressions and exact residue energy maps.
  #link("../../crates/gammalooprs/src/numerator/uninit.rs")[Numerator construction].

+ #strong[Attach a shared numerator family.] For each current UV component,
  color algebra is simplified with `simplify_non_color=false`. Exact affine
  energy-map coefficients become formal arguments of one retained numerator
  body; each residue row supplies its numeric argument values. Scalar roots
  multiply a `numerator_family(...)` call rather than an inlined copy of the
  body. The flat definition store shares entries through `Arc`.
  #link("../../crates/gammalooprs/src/uv/approx/direct_3d/branches.rs")[Branch/family implementation].

+ #strong[Prepare the hard/soft scaling.] Physical masses retain component
  ownership tags. Active energies are replaced by their explicit on-shell
  forms; momentum splitting and the selected affine loop chart introduce
  sums inside those factors. U applies its MUV energy deformation; S co-scales
  the component's physical masses. The scalar root receives the loop-measure
  power exactly once; retained numerator bodies receive the momentum/mass
  substitutions without another measure. The inverse hard-scaling variable
  is then used for the Laurent series.
  #link("../../crates/gammalooprs/src/uv/approx/direct_3d/kernel.rs")[Taylor preparation].

+ #strong[Taylor-expand the dependent algebra.] Known multilinear argument
  slots containing the Taylor variable are normalized with temporary linear
  heads, then their original heads are restored. This distributes linear
  argument sums where required; it does not expand inert projector
  permutations. Independent factors of a root product are kept outside the
  native series. The dependent product undergoes genuine series coefficient
  convolution: for example, the first-order coefficient of
  `(A0+t A1)(B0+t B1)` is `A0 B1+A1 B0`. The common precision budget accounts for
  numerator and scalar-root valuations. Coefficients are stored as new shared
  families, and roots are Taylor-expanded using their opaque coefficient calls.
  This is not a promise of zero expansion; it is a specific boundary around
  the expansion. U keeps through order zero and S through order minus one in
  this inverse-scaling convention, including the root measure. Positive-degree
  H is assembled as S(X)+U(X-S(X)); the forest supplies subtraction signs.
  #link("../../crates/gammalooprs/src/numerator/symbolica_ext.rs")[Factor-preserving scalar series].

+ #strong[Normalize each finalized branch and body.] Attach the remaining
  cograph numerator as another family. Remove the unit-valued mass-support
  tags, put the remaining energies on shell, replace dimension symbols, run
  `simplify_metrics`, perform color-only simplification, and call
  `expand_dots`. These names do not mean every operation is a harmless rewrite:
  color collection uses polynomial collection in selected color factors while
  making unselected coefficients opaque; Schoonschip metric simplification
  parses and contracts a symbolic tensor network; `expand_dots` normalizes dot
  notation and executes the corresponding two-argument metric networks.
  Schoonschip's default has `expand_contracted_sums=false`; its residual-boundary
  fallback nevertheless contains a local `.expand()`, guarded by 512 estimated
  terms and a 100,000-byte input-product limit. Such local contraction/expansion
  is distinct from distributing the entire graph numerator. All these operations
  have completed before the sampled global collector.
  #link("../../crates/gammalooprs/src/uv/approx/final_integrand.rs")[Final normalization],
  #link("../../crates/spenso/src/shadowing/collect.rs")[selective collection],
  #link("../../crates/idenso/src/shorthands/schoonschip/contraction.rs")[bounded residual expansion],
  #link("../../crates/idenso/src/shorthands/metric.rs")[dot expansion].

+ #strong[Sum sources and forests.] Source selectors are materialized for the
  localized route; full explicit mode uses selector one. `zip_add` assembles
  sums with `Atom::add_many` while merging the retained definitions. The forest
  assembly collects color factors, sums signed forest contributions, splits
  the remaining momenta, and removes the metadata wrapper `den(..., D)` while
  retaining D itself. Addition/normalization does not mean a blanket distributive
  expansion of every product. Numerator bodies are still separate definitions.
  #link("../../crates/gammalooprs/src/uv/hedge_poset.rs")[Forest assembly],
  #link("../../crates/gammalooprs/src/uv/mod.rs")[root/definition storage].

+ #strong[Attempt global factor recovery: observed failure.] Before returning
  the complete expression, `into_integrands` protects powers, functions and
  tensor sums with opaque wrappers. Inverses, selectors and numerator-family
  occurrences receive occurrence-local wrappers to keep them inside their
  original branches and contraction contexts. This intentionally prevents
  unsafe sharing across guards and preserves factor scopes. The completed
  scalar root is then sent to `collect_factors()` until unchanged. The retained
  numerator-definition store is not inlined here. Nevertheless, the collector
  densifies its temporary factor maps over the large sum and clones owned
  keys. Opaque arguments protect algebraic scope, not the complexity of this
  global bookkeeping. This pass is unconditional at this boundary; the later
  evaluator's `do_algebra=false` setting does not disable it.

+ #strong[Later, unreachable in the full failing runs.] If factoring returns,
  wrappers are removed and evaluator construction can combine compatible
  numerator interfaces, prune unused definitions, perform finite tensor
  component contractions, and build numerical evaluators. Component contractions
  really create finite sums, but no full-U/H sample in this investigation reaches
  that stage. The successful orientation-58 runs do. Their smaller source
  inventory does not certify the memory behavior of the full sum.
  #link("../../crates/gammalooprs/src/integrands/process/evaluators.rs")[Evaluator preparation].

The compact Taylor counters report 20 completed coefficient-family events for
U and 33 for H. The maximum sum of numerator coefficient-body bytes within one
such event is 18,599 bytes for U and 54,910 bytes for H. These exclude scalar
roots and other retained state and are not total-expression or heap sizes.
They support keeping numerator-body size separate from the collector's
multi-GiB temporary allocations.

== Repair target and limits

The first repair target is the missing-factor/zero-exponent insertion in the
collector and the scope of the unconditional whole-forest factoring pass.
Avoiding zero entries would preserve their implicit identity; a bounded or
sparsely implemented factor-recovery pass should also preserve the opaque
branch and numerator scopes. Removing protections or flattening physical
numerators to work around this would change the wrong representation boundary.
No proposed repair has been applied or validated here.

The ordinary-U case already reproduces the failure without the soft Taylor
terms. The same source line explains the H memory-growth samples. This does
not prove that full generation will finish within 12 GiB after a repair:
remaining evaluator costs have not yet been reached and must be measured.

Debugger attachment paused each target for roughly 9-11 seconds per sample.
The sampled runs remain capped at 12 GiB target-process-tree RSS, with 0.2-second
polling; debugger memory is separate. Their wall times include instrumentation
and are not replacement benchmark timings. The small standalone probe is a
separate uninstrumented measurement. No full expression dump or diagnostic
numerator expansion was used, no license value is retained in the artifacts,
and all investigation processes have stopped.
