# Current goal for the raised-energy CFF handoff

This file is temporary handoff material, intentionally committed under
`HANDOFF_REMOVE_FOR_PR/` for the next agent. Remove this directory before the
merge-ready PR is finalized.

## Immediate first task: resolve the GL04 `T2` discrepancy

Start by reading
`HANDOFF_REMOVE_FOR_PR/HANDOFF_T2_DISCREPANCY.md` in full. Reproduce its
isolated GL04 `{1zs}` / Taylor-`T2` mismatch before changing production code. Finish the
complete one-call generalized-CFF oracle described there, use it to decide
whether the direct local-3D or projected local-4D construction is wrong, and
then locate the first algebraically different boundary. Apply only a
first-principles, topology- and process-independent fix. Do not fit a sign,
special-case GL04, add a process-specific factor, dispatch the original
numerator, or algebraically cancel numerator and denominator factors upstream
of generalized CFF.

Treat local-3D as the stronger current oracle, not as infallible: despite its
simpler construction it could still contain an implementation defect. A
successful UV profile proves cancellation of the profiled divergent limits,
not equality of UV-finite terms. If both local-3D and projected-local-4D pass
their UV profiles while their complete residue sums still differ, classify the
remaining mismatch as potentially UV-finite and resolve it by the complete
local comparison and identical-input oracle rather than declaring either route
correct from the profiles alone.

The non-negotiable architecture is:

1. **Local UV counterterms from 3D:** generate the complete generalized CFF for
   the original graph, write it as a sum over complete residue-map keys (with
   an opaque selector per key when orientation localization is requested), and
   apply the UV Taylor operator directly to each complete keyed expression.
   This lane performs no UV graph reconstruction. Its explicit-orientation-sum
   mode is the same construction with the selectors omitted.
2. **Local UV counterterms from 4D:** apply the Taylor operator to the
   factorized 4D integrand, retain original-edge provenance, construct the UV
   skeleton directly from the original graph and derivative-created serial
   occurrences (never from Kirchhoff/incidence reconstruction), rewrite the
   Taylor numerator in the EMR variables of the resulting factorized graphs,
   call generalized CFF on those factors, and sum all residue maps explicitly.
   This lane is not required to subtract orientation by orientation.
3. Keep the numerator factorized at every production boundary. The original
   numerator remains on its original local graph elements. The EMR dispatcher
   may act only on new temporal factors produced by denominator derivatives;
   every new hard factor should be placed on a corresponding new serial edge
   of that same degenerate line. Retained soft/cograph factors stay with the
   factorized outer graph.
4. Compare local-3D and projected-local-4D only after summing every generalized
   residue-map piece belonging to one complete physical Cutkosky residue.
   Individual contact/remainder pieces and artificial moments need not match.
5. Do not pull LTD support into this branch. Keep the representation interfaces
   clear enough that future LTD residue keys, whose energy rules may be affine
   momentum combinations, can be added without redesigning the UV lanes.

## Establish the all-green UV/CFF milestone

After fixing the `T2` discrepancy:

- Run the focused generalized-CFF, exact-source, EMR/LMB-equivalence, direct
  local-3D, projected local-4D, raised-threshold, and integrated-UV tests.
- Run the complete scalar three-loop LU acceptance matrix in all three modes:
  orientation-local local UV from 3D, explicit-sum local UV from 3D, and
  explicit-sum local UV from 4D. Integrated UV counterterms and threshold
  counterterms must be enabled in the production acceptance variants.
- Run UV profiles on the selected scalar cases. Local UV from 3D must be UV
  finite separately for every complete residue-map key; projected local UV
  from 4D is tested after the complete residue sum.
- Ensure every starting integral is energy convergent. The energy DOD check is
  the ordinary UV DOD calculation with measure contribution one per loop
  rather than four per loop; unsupported energy-divergent inputs are errors.
- Run `just test_gammaloop` and follow the same selection as CI. Ignore the
  obsolete `tests::failing` namespace.

When this line is fully green, commit and push it immediately. The commit
subject must clearly contain `ALL GREEN`. Use the required identity for every
git operation:

```bash
git -c user.name=ValentinHirschi \
    -c user.email=valentin.hirschi@gmail.com ...
```

## Critical simplification/refactoring review

Before returning to physical DDx/TTx validation, critically review the entire
local-UV implementation from both 3D and 4D. Refactor it to be as small,
readable, and explicit as possible while preserving the green tests. In
particular:

- keep the direct and projected lanes in clearly separated modules and types;
- remove obsolete reverse-engineered-local-3D and diagnostic-only production
  paths;
- use names that describe semantic roles rather than contour folklore;
- document all sign, normalization, energy-map, and factor-placement
  conventions at their point of use;
- preserve provenance without introducing a graph-reconstruction algorithm;
- avoid helpers already provided elsewhere in the codebase;
- after every material refactor, rerun the focused route-equivalence and UV
  profile tests before making the next change.

## Physical acceptance and scale validation

Once the all-green/refactored UV core is stable, complete the fast physical
acceptance coverage:

### `e+ e- -> a -> d d~` at NLO

- Use the existing `epem_a_ddx` NLO example with integrated UV counterterms
  enabled and a target of `(alpha_s/pi) * LO` for the NLO correction. Keep the
  phase/sign caveat isolated only where still physically unresolved; do not
  weaken the ratio test otherwise.
- Make the example card run in approximately three minutes on five cores.
- Verify the result in all three CFF/UV modes above.
- Verify each graph is independent of `M_UV`, the complete NLO sum is
  independent of `mu_r`, and the localizing scale cancels.
- Compare the GL0 and GL2 `mu_r` logarithms analytically and numerically,
  including the GL0 multiplicity of two and the distinction between `mu_r`
  and `mu_r^2`.
- Verify the EMR energy bounds: GL2 at most degree one and the GL0
  self-energy-carrying structure degree two where required.
- The paper's direct `gamma -> d d~` values are not numerically identical to
  the `e+e-` process; use only the appropriate constant normalization or the
  inclusive NLO/LO ratio when comparing them.

### Fully-MSbar top-pair production at NLO

- Add the corresponding fast `epem_a_ttx` acceptance in all three modes, with
  local 3D and local 4D UV counterterms.
- Also reproduce the published direct `gamma -> t t~` setup by replacing the
  off-shell photon spin projector with `-g^(mu nu)` and disabling picobarn
  units. Do not compare unconverted `e+e-` numbers directly with the published
  direct-photon numbers.
- On-shell counterterms are still out of scope; use the fully-MSbar scheme.

Use explicit orientation sums for numerical integration, explicit summed LMB
channels with OSE multichannel weights, a few iterations, weak continuous-grid
learning, and no discrete-grid learning. Acceptance statistics must be as low
as practical while remaining conclusive and fast enough for
`just test_gammaloop`; retain higher-statistics runs only as manual validation.

## UV-profile CLI and acceptance coverage

Finish and validate the UV-profile command for amplitudes and LU cross
sections:

- optional graph selection and optional physical Cutkosky-cut selection by
  edge-ID list;
- default `--selected-limits only-divergent`, covering every cycle combination
  with expected UV DOD `>= 0`;
- use the generation LMB when suitable, otherwise the first suitable LMB in a
  deterministic sorted list;
- retain exhaustive all-cycle/all-LMB profiling as an opt-in mode;
- support per-residue-map-key profiling for orientation-local local UV from
  3D, and complete-residue-sum profiling for local UV from 4D;
- `--fail-fast` and a clear colored summary table containing every failure and
  enough details to reproduce it;
- concise reusable test harness and fast acceptance tests for the selected
  scalar LU cases plus DDx and TTx in both UV constructions.

## Merge-ready finish

Record generation and evaluation performance for the three CFF/UV modes,
remove temporary diagnostics, run the complete CI-equivalent gate, synchronize
with the latest `main` without losing user work, resolve all in-scope failures,
commit and push `raised_energy_cff_wip`, and open a clean PR to `main` with the
tests, numerical results, known limitations, and non-LTD scope documented.
