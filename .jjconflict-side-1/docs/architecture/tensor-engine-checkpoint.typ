= Tensor engine correctness checkpoint

This checkpoint consolidates symbolic tensor operations in Idenso's existing
`SymbolicTensor`. Composition, index substitution, interface bookkeeping,
contraction, aliases and tensor rules share that owner. Spynso performs Python
conversion and dispatch. The duplicate structured-atom type and symbolic-network
Schoonschip engine are removed; production numerator callers use the shared
operations and explicit materialization boundaries.

The final fixes preserve dummy scopes in powers, logical ports, typed zeros,
callback-dependent interfaces and incomplete contraction status. In particular,
a callback may change rank after metric substitution, and an exhausted budget
must not acquire a completion certificate. A basis-selector argument in the
registered delta tensor is metadata, not another open port.

== Validation scope

Qualification combines saved broad runs with targeted runs after the final
corrections. It is not a claim that one monolithic CI run passed on this commit.

- The saved Idenso/Spenso broad run passed 1,140 of 1,142 tests, with 22 skipped.
  Its two failures were outdated capped-frontier and uniformly transposed-trace
  expectations. Their corrected cases pass targeted checks.
- The final contraction checks passed 79 focused regressions, including 72
  callback tests, and 49 HEP integration checks, with one skipped. The GL03
  comparison now uses exact components instead of comparing independently
  rounded floating-point evaluations. Original momenta and graph inputs remain.
- Native Spynso passed 127 unit tests and ten integration tests. Symmetric
  metric cooking may normalize argument order; its declared layout is restored
  explicitly. Nonsymmetric cooking preserves exact arguments and logical ports
  in both tested orders.
- The ghost-crossing fixture retains its UFO coupling, momentum and colour
  oracles. Its external Grassmann-order signs follow the declared convention;
  final arithmetic cancellation uses `collect_factors` and `collect_num`, preserving opaque
  tensor leaves and the literal-zero assertions.
- The production checks include 110 selected GammaLoop/Vakint tests and a rebuilt
  GL000 three-loop aa-to-aa smoke run with 20 finite samples. This smoke run is
  not a precision-integral or whole-amplitude parity certificate.
- The refreshed Wasm passes 21 actual Typst runtime checks. Relevant strict
  Clippy, formatting and generated Python-stub checks pass.
- Sixteen saved ladder, trace and production cases pass the FORM/HEP comparison
  replay. The production comparison uses all 256 external four-dimensional
  components at three fixed momentum assignments; it does not establish a
  general kinematic or arbitrary-dimensional Clifford identity.

The installed Python run retains one known baseline failure: Graphica's directed
fermion tadpole normalization gives an extra factor of one half. It reproduces
with both the earlier and current Symbolica revisions. The original assertion
is retained, with a
#link("../../examples/reproducers/symbolica-graph/report.typ")[standalone reproducer].

== Saved measurements and notebook

Symbolica is pinned to `6a96c9d77b217ea394770524c119e0613d4d9ef5`. Its
non-polynomial expansion improves markedly on overlapping factors, but is not
uniformly faster than polynomial expansion. The
#link("../../examples/reproducers/symbolica-expansion/performance.typ")[measurement report]
retains the 19-case comparison. No new end-to-end timing claim is inferred from
those microbenchmarks or from the final correctness replay. Additional
performance work is deferred at the user's requested stopping point.

The ladder and gamma notebook, FORM programs and saved comparisons are in
#link("https://github.com/alphal00p/symbolica-community/blob/3da847a9c8f9453ead70c2328860cf5e68a74298/examples/hep/gamma_simplification.py")[the HEP examples checkpoint].
