= Reimplementing the LUDY DY and DIS calculations

== Outcome and scope

The executable `compute.py` is a partial reimplementation through scalar
integrands and pre-IBP integral indices. It uses the current FeynKit community
stack: Spenso for dimension-dependent tensor expressions, Idenso for Dirac
algebra, Symbolica for exact scalar algebra, and FeynKit for raised cut
distributions. Most of the present implementation is in the first three
packages; it is not an end-to-end FeynKit computation of the integrated NLO
coefficients.

The physical inputs cover all eight DY graph records and all six DIS rooted
graphs, including the two inequivalent quark self-energy roots. They retain
17 DY cut-sector records (one has zero weight) and 15 DIS iterated cut-sector
records. Explicit physical input is necessary: the Python graph surface does
not yet ingest these LUDY propagator/cut records directly.

The program computes:

- Eight transverse DY traces and three longitudinal qg traces, with finite
  incoming virtualities and symbolic Lorentz dimension.
- The metric, transverse-vector and Ward projections of all six DIS traces,
  followed by their F1, F2 and FL projections.
- Propagator denominators from momentum vectors and mass squares.
- Exact multivariate partial fractions for both crossed boxes.
- Affine denominator coordinates and scalar integral-index shifts for eight
  DY numerators and eighteen DIS F1/F2/FL numerators. Every resulting scalar
  integrand is reconstructed exactly before its cut metadata is attached.
- The common two-particle-cut DY qg Ward cancellation. This requires both
  incoming legs on shell, `k²=0`, and `(k+p1+p2)²=m²`; it is not a vanishing
  uncut rational integrand.
- Formal individual FeynKit cut distributions, preserving propagator powers,
  and cut-support filtering of the scalar terms. Six independent residue
  checks exercise powers one through three, both energy orientations, and a
  nonconstant rational coefficient.

The JSON output explicitly states `stage="pre-IBP"`,
`integrated_coefficients=false`, and `joint_cut_evaluated=false`.
It contains no fabricated integrated coefficient or successful pole-cancellation
claim. Coupling/color normalization, channel phases, and scheme weights must
still be applied by the downstream physical assembly.

== Reproduction and provenance

The reference is the private
#link("https://github.com/ZenoCapatti/LUDY-resources/tree/e5554019e398bbbc114ae3dafedc88bf09cb6eeb")[LUDY-resources revision e5554019e398bbbc114ae3dafedc88bf09cb6eeb].
Authenticated GitHub access succeeded. Reference source files are:

- `calculations/dy/python/notebook.py`: physical line lists, eight raw chains,
  three Ward chains, and the independent expanded-polynomial DY oracle appendix.
- `calculations/dy/python/physical.py`: graph/cut semantics and the downstream
  native assembly.
- `calculations/dis/python/physical.py`: five topology inputs, six rooted
  traces, finite-virtuality kinematics, DIS projectors, and fifteen cut records.
- `notebooks/ludy.py`: the existing reference trace, apart, numerator reduction,
  and Kira interfaces.
- DY `masters.py`, DIS `masters.py`, `assembly.py`, and `native_scheme.py`:
  integration, master conventions, raised external cuts, and physical assembly.

`inputs.json` contains physical inputs only. A vector is the integer triple
of coefficients in `(k,p1,p2)` for DY or `(k,p,q)` for DIS; null denotes an
unslashed gamma matrix. Repeated labels contract Lorentz indices. The dimension
is symbolic while the spinor trace has dimension four.

`oracles.json` is a separate validation artifact. DY transverse values come
from the reference's expanded-polynomial appendix, independently checked
against its raw-chain implementation during export. DY Ward and DIS values
come from the reference implementations, so they validate the migration of
inputs and projectors, not an independent Dirac-algebra engine. The port
reads the oracle only after completing its own computation, with `--check`.

Use a current combined Symbolica/FeynKit/Spenso/Idenso host:

```sh
python examples/ludy/compute.py --check --output /tmp/ludy.json
python examples/ludy/compute.py --process dy --check
python examples/ludy/compute.py --process dis --check
```

For the community host that calls the FeynKit module `hep`, add `--module hep`.
There is no automatic import fallback. The port uses
`TensorExpression.gamma(d)` and the current `Expression.solve` API. The pinned
reference uses the older `TensorName.gamma` API, so regenerating its fixtures
requires its own locked environment, separate from the environment running
the port. Neither execution nor validation starts Kira or downloads data.
To regenerate the physical inputs and independent reference values, run this
from the repository root with the pinned LUDY environment's interpreter:

```sh
/path/to/LUDY-resources/.venv/bin/python examples/ludy/export_reference.py \
  /path/to/LUDY-resources /tmp/ludy-reference-export
```

The exporter checks the reference commit and rejects changed tracked calculation
sources. It executes only the pinned input/oracle cells, not the full Marimo app.
Inspect the generated JSON before replacing source-controlled fixtures.

== Validation status

Both process branches passed `--check` in the installed `hep` host: fourteen
physical graph records, thirty-two sector records, twenty-six exact scalar
integrand reconstructions, the DY cut Ward identity, and six raised-cut residue
checks. All forty-seven reference comparisons pass: eight DY transverse
numerators, three DY Ward projections, and six projections for each of the
six DIS graphs. Ruff formatting/linting passed and this Typst document compiled.
No Rust code was changed by this example.

The fixture exporter replaces dimension-dependent scalar products before
renaming the dimension. Reversing that order changes the expressions before
the scalar-product patterns can match. The regenerated eleven DY expected
expressions now use the same scalar coordinates as the port; physical inputs
and DIS reference expressions are unchanged. The integrated NLO calculation
remains outside the implemented pre-IBP boundary.

== Conventions that must survive a complete port

The DIS kinematics retain `p²`, `q²=-Q²`, and
`p.q=(Q²/x-x p²)/2`. The projection vector is
`P=p+(p.q/Q²)q`, with `P²=p²+(p.q)²/Q²`. The code contracts this vector directly;
Spenso distributes its components, replacing the reference's four separate
pp/pq/qp/qq traces. It verifies `F2=2x F1+FL` for each rooted graph.

DY loop kernels exclude the incoming virtuality carriers, exactly as in
the reference numerator-reduction boundary. The qg triangle/bubble also retain
the explicit inverse first/second powers of `(p1+p2)²`. These kernels are not
complete normalized cross sections. DIS external D5/D6 denominator powers
remain in the scalar coefficients; their derivative cuts still need physical
assembly.

The cut records preserve initial and final memberships separately. Their union
contains each shared line once. Multiplying independent cut deltas for both
memberships would double count a shared line. The displayed positive-energy
distributions are formal local factors with independent energy placeholders;
they do not establish physical energy support, routings, or joint Jacobians.
DY reference master weights and observable phases are separate conventions
from FeynKit's default `-2 pi i` cut normalization.

`supported_terms` identifies terms with every cut loop line present at positive
power. This is a pre-IBP bookkeeping result, not an evaluated discontinuity.
The DY crossed-box `00` sector requires its separate analytic maximal-cut
contribution. No three-denominator apart child supports all four cut lines;
an empty supported-term list must not be interpreted as a zero physical sector.

The reference documents differences among historical Mathematica output,
earlier manuscript expressions, and its corrected native coefficients.
Those older closed formulas are not substituted for a calculation here.

== Missing functionality and existing owners

=== Physical input to native graph objects

Existing: `Model.generate_diagrams`, `FeynmanDiagram` import, momentum bases,
subgraph selection, numerator/denominator extraction, and ordinary generated
final-state cuts. Rust already has `FeynmanDiagram::builder`; serialized Python
imports expect finalized model identities and graph annotations.

Missing for this port: a supported Python construction path from ordered
physical lines, momentum routing, local tensors, and iterated cut membership.
This should expose or extend the existing graph builder, rather than create
another graph schema. A minimal acceptance case is the DIS `q_lsz` graph with
a squared external carrier; the second is the DY crossed box with its four
distinct cut sectors. The current script keeps reference inputs as data and
does not pretend its traces were generated by the Standard Model generator.

Native generation probes in the available installed host encountered
Symbolica's restricted thread allowance. The generator enters a Rayon pool
even when `threads=1`. This is a runtime validation limitation, distinct from
the missing LUDY input bridge; ordinary generation support already exists.

=== Iterated discontinuities and raised external cuts

Existing: `DiagramCut`, `CutPropagator`, cut orientations, raised powers, and
complete energy-residue derivatives. CFF also has surface/pole/residue tools.
It would be incorrect to request a new generic raised-cut primitive.

Missing: composition into two ordered discontinuity channels, shared-line
ownership, physical support, relative phases, external-virtuality derivative
actions, and their Jacobians/measure factors. The reference's D5/D6 derivative
acts in a virtuality coordinate; applying an energy cut to only the numerator
does not implement it. Acceptance cases: DIS `q_box/triple123`,
`q_bubble/external_external`, and `q_lsz/external_external`.

=== Cut-aware scalar-family reduction and IBP

Existing: Symbolica provides `apart`, exact solving, and polynomial
coefficients. The small explicit pipeline here uses those operations rather
than reimplementing elimination. FeynKit's tensor reducer performs vacuum
isotropic reduction; applying it to these scattering denominators would erase
physical dependence on the external momenta.

Missing: an integral-family object tied to native graph lines and propagator
powers, cut-preserving numerator-index reduction, cut specialization, and an
IBP backend boundary. The LUDY `kirapy` package already owns deterministic Kira
projects, exact expression transport, cache manifests, and master selection.
Reuse that interface or agree an adapter before adding a second job runner.
The reference DY computation needs eighteen projects; DIS needs sixteen
ordinary projects and two additional raised-external projects. Acceptance
should compare target indices and selected-master coefficients, not just
final finite formulas.

=== Analytic masters, branches, and normalization

Existing: LUDY has independently derived bubble/triangle cut masters, massive
variants, and the crossed-box maximal-cut master. Vakint supplies vacuum
integral machinery; this does not cover the scattering cut-master basis.

Missing in FeynKit: analytic evaluation for these cut families, a branch and
phase contract, the positive phase-space normalization, and the explicit
external `dp²/(2 pi)` factors. An initial integration boundary should accept
the reference's exact master records without requiring its Marimo UI. The
analytic four-line DY `00` contribution must remain an explicit separate
route when an apart decomposition has no supporting child.

=== Virtuality expansion and endpoint distributions

Existing: Symbolica series/differentiation and FeynKit's local UV Taylor
expansion. UV scaling is a different limit from the incoming-virtuality
expansion used here.

Missing: leading-power expansion with physical virtuality boundaries,
epsilon-dependent virtuality integration, delta/plus/regular endpoint algebra,
convolution and endpoint normalization, and distribution-valued Laurent
assembly. Generalized energy deltas in the cut API do not implement
`[log(1-z)^n/(1-z)]_+`. Acceptance cases should include endpoint moments,
double/simple pole cancellation, and the independent DIS/DY one-leg scheme
comparison, retaining the incoming-gluon polarization conversion.

=== Dimension-dependent incoming spin averages

Existing: the current source includes particle spin sums, including covariant
and axial vector projectors, and Spenso can promote Lorentz dimensions.

Missing for the reference normalization: the spin-sum averaging API is explicitly
four-dimensional; DIS requires the gluon factor `1/(d-2)` and its associated
finite scheme conversion. Promoting slots alone does not promote a previously
applied factor `1/2`. The script leaves channel prefactors for the downstream
assembly rather than silently equating these conventions.

== Suggested implementation order

First expose the existing graph builder and define iterated cut metadata,
using the fourteen graphs and thirty-two sector records here as acceptance
inputs. Then connect graph-owned scalar families to kirapy and the existing
LUDY analytic-master boundary. Finally add virtuality and endpoint assembly,
requiring cancellation of the physical Laurent poles and the DIS/DY scheme
comparison. D-dimensional incoming spin normalization belongs in that same
validated assembly path. These are proposed extensions, not new APIs added by
this exploratory port.
