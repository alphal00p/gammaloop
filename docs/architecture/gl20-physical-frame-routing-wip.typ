= WIP: physical graph frame and paired routing root

This branch checkpoints the combined, source-reviewed GL20 compatibility
candidate. The corrected cohort passed four Graph contracts and seven Generator
contracts, but the complete original GL20 regression failed its unchanged golden
snapshot. All 16 events now match the golden event indices, cut metadata, and
incoming and outgoing momenta exactly. Numeric weight differences remain, so this
candidate is not a sufficient GL20 fix. NixCI remains disabled; no production
dependency pin, final-review status, Clippy, or full CI success is claimed.

== Generation and routing boundaries

For cross sections with left/right symmetry and incoming particles, the
candidate restores the historical paired external momentum assignment after
each skeleton pass. The incoming occurrence is assigned in rebuilt edge order;
a role change conjugates the all-incoming species while preserving its physical
particle. Amplitude, vacuum, and left/right-disabled paths retain their scope.

The same path retains the raw physical representative's historical graph order
before finalization. Its key uses the stored node and edge order, external tags,
vertex rule names, incident edge IDs, endpoints, and signed PDGs. It does not
re-canonicalize the representative or change comparison numerators, ratios, or
grouping keys.

For genuinely sewn external carriers with both endpoints, automatic loop routing
starts at the owner of the first included native half-edge. The existing graph
query determines that owner. Existing roots and the disconnected fallback remain
available; the code contains no fixture node or loop-basis override.

== Retained witness

The original test selects its physical topology by colored connectivity; its
native label in the retained current population is GL04. In the previous paired
assignment-only candidate, the actual selected member 4 matched historical
member 21, while its discarded counterpart 7 matched historical master 20.
The historical raw graph ordering chooses the latter counterpart.

Both counterparts previously used default loop edges 3, 4, and 7. Exact integer
DFS replay showed that the native first-half-edge owner restores the historical
automatic loop edges 4, 5, and 9 on both role frames. Those integers are witness
outputs rather than production constants.

The previous paired assignment-only branch passed three focused tests but
failed the complete original golden snapshot; it remains a separate public WIP
checkpoint. A withdrawn singleton shortcut is excluded from this candidate.
Runtime counterterm helper signs, model inputs, original fixture, golden
snapshot, tolerances, and scientific assertions are unchanged.

== Validation results

The sealed private plan first runs four Graph contracts and seven Generator
contracts. Graph contracts include a genuinely paired sewn box: the automatic
basis is checked for vertex and paired-cut conservation, while the previous
coordinate-specific cut-shift assertions remain on an explicit alternate basis.
Both bases retain graph and numerator identity and exact basis-change checks.
Generator contracts cover paired role conjugation, mirrored counterparts,
physical member ordering, scope guards, and half-edge and vertex-slot integrity.

Only after both focused suites pass does the plan compile the remaining matching
dependencies and run the complete original
`gl20_multichannel_local_inspect_event_snapshot` once, with `INSTA_UPDATE=no`,
the unchanged semantic selector, model, settings, point, golden, and assertions.
The first focused Graph harness compiled with Rust 1.98.1 and ran once. The
paired-box conservation and paired-mirror tests passed. Two other new tests
failed when full diagram validation reported `MissingCrossSectionCuts`; those
synthetic paired builders did not install the cut metadata required by that
validator. Their preceding routing and basis assertions passed up to that
validation call. This does not establish both complete flow contracts.

The controller stopped at this gate. It did not compile the remaining owners,
run Generator contracts, generate the original fixture, or evaluate its point.
No retry ran. Source, fixture, and golden hashes passed the post-failure guard.
A separate corrected-fixture cohort preserves that failed run and changes only
the new Graph test child. It installs complete cut partitions through the
existing `with_cut_partitions` constructor. The paired fixture uses retained
historical cut 0, including both attachment crowns; its mirrored counterpart
reverses the sides. The disconnected fixture uses an explicitly synthetic cut.
All 18 prior assertions remain, with eight completeness and crossing assertions
added. Production Graph and Generator sources are byte-identical to the first
cohort. Independent source review passed. The corrected cohort compiled all 12 planned
artifacts, passed all four Graph and seven Generator contracts, and then ran the
complete original GL20 test once with `INSTA_UPDATE=no`. Setup took 1104.7 seconds;
the original test ran for 7.55 seconds and returned 101 after a snapshot failure.
Source, fixture, and golden hashes passed the post-failure guard. No retry ran.

All 16 reconstructed events now equal both current and historical golden data
for event order, complete cut metadata, and incoming and outgoing momenta.
The complete retained diff changes 78 numeric rows. The current aggregate is
`(-1.4716101189747047e-29, -4.5025928118207965e-30)` in `(re, im)` order;
the current golden is
`(-1.4716101189741638e-29, -4.502592811820771e-30)`.
The large imaginary-sign discrepancy is resolved, but these strict numerical
mismatches remain. The historical golden is a separate numeric reference and is
not substituted for the current one. Stored numerator, evaluator IR, and
operation ordering remain the next comparison boundary; no tolerance or golden
change is proposed. No five-loop physical
execution is involved.
