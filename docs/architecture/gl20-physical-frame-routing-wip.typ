= GL20 physical frame and strict regression

The physical-frame and routing fixes pass the original Main GL20 test with its
exact JSON assertion. The snapshot migration changes exactly 78 numeric leaves;
the other 669 leaves, header, test source, model, settings, semantic topology
selector, event order, cuts, channels, and momenta remain unchanged. The complete
original test passed once with `INSTA_UPDATE=no`: one passed, zero failed, in
7.58 seconds. This is focused validation; consolidated full CI remains separate.

== Physical ownership boundaries

For cross sections with left/right symmetry and incoming particles, generation
restores the paired external momentum assignment after each skeleton pass.
Changing the incoming occurrence conjugates the all-incoming species while
preserving the physical particle, attachments, and vertex slots. Amplitude,
vacuum, and left/right-disabled paths retain their scope.

Physical representative selection retains the historical stored graph order:
external tags, vertex rule names, incident edges, endpoints, and signed PDGs.
It does not re-canonicalize the representative or change comparison numerators,
ratios, or grouping keys. The retained witness selected current member 4 instead
of its physical counterpart 7; historical physical ordering chooses the latter.
Those member numbers are evidence, not production constants.

For genuinely sewn external carriers with both endpoints, automatic loop routing
starts at the owner of the first included native half-edge, obtained through the
existing graph query. Unpaired, outgoing, vacuum, and disconnected fallbacks
retain their contracts. The code contains no fixture node, basis, or helper-sign
override.

Four Graph and seven Generator contracts passed. They check paired roles and
mirrors, physical ordering, particle and slot integrity, complete cut metadata,
vertex conservation, basis transformations, and scope gates. Coordinate-specific
cut shifts remain tested in their explicit alternate basis. An earlier synthetic
fixture omitted required cut metadata; its failure is retained separately, and
the corrected fixture keeps all previous assertions.

== Numerical diagnosis and snapshot migration

The corrected saved state matches the historical physical frame exactly,
including all 35 loop bases and their integer routing matrices, default basis,
and four active channels. All 16 events match the historical and Main golden
event identities, cut metadata, and incoming and outgoing momenta. The earlier
meaningful counterterm imaginary-sign discrepancy is resolved by these physical
ownership fixes.

Matched retained evaluators were compared on the same recorded input bits,
including historical override input 16 equal to one and the actual GC_94 binding.
The active model inputs agree. The old Arb program contains binary64 pi widened
to its working type; the current program uses native-precision pi. Their explicit
normalization is pi to the power minus eight. The f64 pi bits agree, so this
coefficient distinction does not explain native f64 rounding differences.

One further native-gated comparison observed actual 1000-bit inner storage in
both owners at the inputs, roots, and selected energy/coefficient checkpoints.
All eight on-shell energies agree exactly across owners. The retained rounded
input has cut residual `-5.25147579242245e-14` even when lifted exactly to that
precision. A separate formal diagnostic replaced only the dependent E11 energy
by the cut relation. After accounting for the independently proved pi factor,
the old/current root difference fell from `1.6123e-23` to `7.4175e-301`.
This establishes the mechanism for that retained case; the intervention is
neither a native-equivalent trace, a new physical point, nor a general identity
proof.

A separate forced-Arb full-point reference retains the existing decimal lift of
the public f64 coordinates. Its 1000-bit metadata is configured/wrapper-reported;
the single identity rotation returns no accuracy estimate. It is not an interval
or convergence certificate. The reference and retained f64 results agree closely
in complete complex norms; tiny imaginary cancellation residues are not used as
independent relative-error scales. No numerical tolerance is introduced.

The Main fixture now records the actual corrected f64 output. Its old SHA256 is
`1fcc97bcb754b8aefd936530108d4f7d8281b8ea2e4c6f17e199776363ac7a0b`;
the new SHA256 is
`32acac110c9e3859753585e17787bbd559bd3f5547dea7c26c4092a59e8151a2`.
The historical golden remains a distinct reference. The Main fixture is not
transplanted into the separate maintained FeynKit owner fixture.

The strict-run completion SHA256 is
`d1fe6e0fa2dd62f21b913cd5c5ee31d9fd2af9e2e2e27b0ec7a5e105ca7e5435`;
independent final review is
`69c02e89358ced25ec7cdde33a3e35d65f8ee91f85b95fff69125867f8a2e0fb`.
One test harness was compiled against nine authenticated unchanged libraries.
Compilation took 62.8 seconds and is separate from the test's 7.58-second runtime.
A preceding controller import failure occurred before compilation or evaluation
and is preserved; it is not a scientific retry or a test pass.

== Source checkpoints and integration

The source supplements in `docs/validation/gl20-2026-10-26/` preserve completed
diagnostic methods and controllers without state payloads, binaries, environments,
licence values, or bound activation records. They are source references rather
than standalone launch bundles.

The focused owner fixes remain public on
`gl20-physical-frame-routing-restoration`. Their owner integration belongs in the
maintained FeynKit publication work, including its existing PR 220 review path;
the Main-specific snapshot migration remains on the Main-based GL20 branch.
No protected FeynKit branch, Native dependency pin, five-loop worker, cache, model,
runtime counterterm sign, or numerical acceptance policy is changed here.
