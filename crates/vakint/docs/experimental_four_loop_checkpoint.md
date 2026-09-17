# Experimental four-loop RustRed checkpoint

This checkpoint records finite-target development evidence. It is **not** a
four-loop production artifact, a proof that the candidate rules close the full
family, or a production `EvaluationMethod::RustRed` four-loop backend. The
experimental adapter remains feature-gated and absent from the normal method
registry.

## Numerical status

The unchanged comparative harness has passed parent, dotted-parent and pinch
inputs for four supplied families:

- H, with 26,956 cold candidate-rule applications for the dotted target;
- FG, with 3,362 cold candidate-rule applications for the dotted target;
- X, with 82,637 cold candidate-rule applications for the dotted target;
- BMW with an explicit first-propagator power of three, with 16,113 cold
  candidate-rule applications for the dotted target.

The BMW power-two probe is itself a declared finite residual, so it is not used
as recurrence evidence. The power-three target changes only the requested
integer powers, not the family, solver, or adapter. In each accepted experiment
the terminal catalog was built offline from the actual residual keys reached by
RustRed, using FMFT as an oracle. After clearing RustRed's point cache, the
candidate scalar tail used an invalid FORM path and performed the positive rule
counts above. FORM was still used by the separately labelled offline catalog
preparation and legacy FMFT comparison lane.

These are twelve finite-target comparisons. Separately, all fifteen existing
legacy four-loop FMFT reference entrypoints pass after the shared exact-before-
approximate finalizer correction. Those fifteen legacy passes are not RustRed
candidate-mode passes. No four-loop artifact or numerical terminal catalog is
currently shipped.

## Frozen timing run

On 2026-09-17 the already-built test driver
`dd4e435825c6baa182d311068ce99d2cbe93b7ebad1f529577ade51fb009e528`
was run with affinity restricted to CPUs 50–51, two candidate-generation
workers and three timing repeats. Its linked Vakint and RustRed release-library
hashes were respectively
`0af8ab8e75e55742e822038a5d50b7745f82ed525e01c6b4fad159dc9308cc94`
and
`d8356baaa569e50496149650e619c22d836f8ba175ac71aa95ca3c3ac328c971`.
The immediately preceding C++ canary on CPU 50 had ended before these runs
started. Other work elsewhere on the shared host was not controlled.

Generation, offline terminal evaluation and catalog construction are outside
every interval below. The scalar boundary starts from an already tensor-reduced
expression. The full boundary reruns the FeynKit prepass inside each interval.
FMFT includes FORM process startup and I/O. `RustRed cold` clears pointwise
memoization immediately before evaluation; `RustRed warm` is the following
cache-reusing evaluation. “Cold” is therefore logical-cold candidate
application in the same live process: it does not reset Symbolica or allocator
state, Vakint objects, OS page caches, or the already prepared offline catalog.
Every timed result was compared numerically with FMFT outside the timer and
passed. Times are medians of three wall-clock samples.

| Family and target | Boundary | FMFT | RustRed cold | Cold/FMFT | RustRed warm | Cold rule applications |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| H dotted | scalar | 3.168 s | 38.271 s | 12.08× | 0.0336 s | 26,956 |
| H dotted | full FeynKit | 1.934 s | 37.749 s | 19.51× | 0.0349 s | 26,956 |
| FG dotted | scalar | 0.968 s | 1.873 s | 1.93× | 0.0406 s | 3,362 |
| FG dotted | full FeynKit | 0.969 s | 1.868 s | 1.93× | 0.294 s | 3,362 |
| BMW power-three dotted | scalar | 2.628 s | 16.481 s | 6.27× | 0.0586 s | 16,113 |
| BMW power-three dotted | full FeynKit | 0.701 s | 4.477 s | 6.39× | 0.0612 s | 16,113 |

The H, FG and BMW parent and pinch targets are declared residuals and therefore
show zero rule applications. Their fast native timings are terminal/catalog
lookups, not recurrence performance, and are intentionally excluded from the
comparison table.

The complete timing processes also reran numerical acceptance successfully.
Their candidate-generation observations were 62.277 s for H, 145.288 s for FG
and 301.241 s for BMW; these setup times are not paired solver benchmarks and
were strongly affected by shared-host load. Process-resident-memory readings
during the timing phases were approximately 935,000 KiB for H, 307,000 KiB for
FG and 592,000–596,000 KiB for BMW. These are snapshots, not phase peaks and
not FORM-child memory.

## Interpretation

The frozen scalar applier does not yet meet the requested “within 75% of FMFT”
objective on a genuine recurrence. FG is the closest at about 1.93× FMFT,
while H and BMW expose substantially larger cold-application costs. Immediate
warm reuse is fast, showing that memoization works, but it does not excuse the
cold gap. Profiling and optimization remain required after the checkpoint.

The measurements also show why terminal-only inputs cannot stand in for scalar
reduction benchmarks. Future performance work should retain a mandatory
positive rule-application count, paired lane order, fresh-cache intervals and
the same post-timing numerical comparison. Production four-loop acceptance
still additionally requires exact original-source replay for every rule used,
a finite-target or whole-domain coverage claim with its precise scope, and a
shipped terminal catalog independent of runtime FORM.
