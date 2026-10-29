= Five certified ungrouped source44 snapshot migrations

The migration is approved at `7a9f24095862a1edd98fd929800fa6665c65f16b`
and reconciled into the ordered MAIN44–45 candidate. Its qualified later failures
remain limits; a NEW independent combined MAIN44–45 review is required before46.

The independent combined Astra review of `121ac233d5ab014f4fadefccce702041bbbdf487` approved the fixed five-row scalar-weight certificate only. This prepared migration changes exactly five inline expectations in `slow` tests. It does not change a request, model, parameter card, filter, grouping policy, worker selection, state lifecycle, production physics, lock, or NixCI setting. The earlier approved `failing::test_generate_aa_ttx_cross_section_slow_filter` row stays at `874 | -186`; the grouped sunset and other grouped rows remain untouched.

The #link("evidence/generation-integration/wp3/ghost-weight-migration/receipt.json")[migration receipt] gives the exact source/commit/feature/kernel/binary/log identities, old and new assertion lines, all checks and retained artifacts. The #link("evidence/generation-integration/wp3/ghost-weight-migration/check-five-row-diff.mjs")[line inventory] verifies that only source lines 2334, 2366, 2402, 2404 and 2407 changed, and each changed line differs only in its literal. All other 2,584 lines of `test_feyngen.rs` are byte identical to the approved certificate source.

#table(
  columns: (1.8fr, 0.9fr, 0.9fr),
  table.header([Module-qualified case], [Prior literal], [Migrated literal]),
  [`slow::test_generate_epem_a_qqh` QCD2], [`594 | -261`], [`594 | -201`],
  [`slow::test_generate_aa_ttx_cross_section` QCD2], [`874 | -266`], [`874 | -186`],
  [`slow::test_vacuum_amplitude_generation` nonfactorizable L2], [`5 | -11/12`], [`5 | -23/12`],
  [same L3], [`30 | -17/12`], [`30 | -13/12`],
  [same L4], [`200 | -199/48`], [`200 | 457/48`],
)

== Actual assertion behavior

An immutable exact-source old-literal test executable failed the named electron and photon slow tests at the old QCD2 expectations, with actual results `594 | -201` and `874 | -186`. The freshly built final Cargo test executable reaches and passes those QCD2 cells after their original earlier cells. Each complete test then fails at a later *unchanged grouped* assertion: electron row 7 has different group IDs at the same total, and photon row 4 has a different grouped representation at the same total. These test functions are not claimed to pass as wholes.

The original vacuum test stops at an earlier unchanged grouped two-loop assertion, before the three target nonfactorizable cells. A separately source-bound scoped harness copies the original scalar generation helpers, the exact `import model sm.json` and `global.n_cores.feyngen=16` setup, and each of the three original command and inline assertion lines into independent tests. All three old-literal tests fail with their certified actual values; all three migrated-literal tests pass with one real test each. The harness validates those cells only. It neither executes nor waives the grouped prefix, and it does not establish a full pristine vacuum-test pass.

The final executable also ran the unchanged ignored diagnostic. Its 4,239 DOT files, 15 graph TSVs, 15 process definitions, and 63 live input files match the approved producer capture; only the isolated `state.folder` value was normalized in 21 settings files. All five sums and physical graph witnesses therefore remain bound to the independently approved certificate. The original round1/round2/round3 joins, source43 baseline, binaries, logs, reviews, and all approved refs remain intact.

== Review boundary

This migration is prepared for a *new independent Astra migration/test review*. It does not certify the grouped source44 rows, source46 combined color-cycle behavior, the complete test module, main integration, or the generation gate. No graph numerator was expanded and no snapshot was broadly rerecorded.
