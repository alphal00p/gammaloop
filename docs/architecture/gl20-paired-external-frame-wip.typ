= WIP: paired external momentum frame

This branch checkpoints a reviewed generation candidate. Three focused tests
pass; the complete original GL20 regression is still pending at this source
checkpoint. NixCI remains disabled. This is not a final-review or delivered-fix
claim.

== Source boundary

The historical generator calls `canonize_external_momenta_assignment` after each
of its two skeleton-canonicalization passes. For a sewn incoming/outgoing pair,
that operation assigns the first occurrence in rebuilt edge order to the
incoming momentum slot, then assigns the other occurrence to the outgoing slot.
Changing the role conjugates the edge's all-incoming particle species. The
physical species attached to the paired external momentum stays unchanged.

The current generator omitted that operation after rebuilding the graph. This
candidate restores it on `ResolvedProcess` for cross sections with left/right
symmetry. It runs after each pass, preserves stored node and edge order, and
validates all proposed updates before committing them. Amplitude, vacuum, and
left/right-disabled paths retain their existing behavior.

== Evidence and exclusions

The GL20 regression selects its topology by colored connectivity; its current
native label is GL04. The retained witness was matched against the entire
ordered CI graph, including endpoints, particles, external metadata, slots, and
half-edge order. A literal replay of the historical external-assignment step on
that witness moves the incoming attachment from vertex 6 to the historical
vertex 7 and conjugates the four external half-edge particle species.

Earlier matched original-point diagnostics located the significant imaginary
sign difference at the side-specific integrated counterterm helper: six
geometry inputs matched bit for bit while the selected side differed. Runtime
helper signs are unchanged by this candidate.

A separate raw-singleton frame-preservation experiment passed its focused
contracts but left the complete failing GL20 snapshot byte-identical. That
withdrawn experiment is not included in this branch.

== Validation at this checkpoint

The private matching-dependency Rust 1.98.1 harness passed all three focused
`paired_external_momenta_` tests. They cover role assignment and conditional
particle conjugation, mirrored and idempotent counterparts, amplitude/vacuum
scope guards, atomic rejection of invalid metadata, and the existing UFO slot
matcher with half-edge IDs, flow, spin, and color.

The complete original `gl20_multichannel_local_inspect_event_snapshot` test is
being rebuilt and will run once with `INSTA_UPDATE=no`, its original semantic
selector, fixture, model settings, point, golden snapshot, and assertions.
Full regression success, Clippy, and broader CI are not claimed here. No
production dependency pin is changed.
