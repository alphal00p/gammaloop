= Grouping-frame numerator digest migration

== MAIN39–43 integration

This approved preparation record is imported into the ordered MAIN39–43 candidate.
The original source, proof and test bindings below remain historical; current combined
checks and approvals are recorded in #link("evidence/generation-integration/wp3/main39-43/receipt.json")[the integration receipt].
A NEW independent Astra combined review is required before44.

== Historical preparation record

This record proposes exactly three native-schema numerator digest updates in `crates/feynkit-generator/tests/generation_oracle/frames.rs`. The source rows, diagrams, topologies, assignments, inverse symmetry, weights, content, names, ids, grouping options, prefactor and named loop-momentum bases are unchanged. The candidate awaits a new independent Astra migration and test review. No integration or publication has occurred.

== Entitlement and preservation

Independent Astra Scope A (`target/astra-r2-scope-a-20261028T200300Z/review.json`) and Scope B (`target/astra-frame-scope-b-round2-20261028T200020Z/review.json`) both approve the bounded certificate at `7dd0b0c160e9c1b62533b54217f1d62275dde157`. That exact commit is preserved at the Jujutsu bookmark `generation-frame-certificate-approved`. The shared Git store at `/shared/lcnbr/dev/gammaloop/main/.git`, as resolved by `jj git root`, resolves the three older tags to `c5a19031d8a527ca61903728b6465187f3e87928`, `67f36295d41cdf872d978692c8aff93d3f3b855e` and `3718714c791d40fb746dd9abc6eeb61fe86548c4`. The reported visibility discrepancy was a review presentation issue; no tag was missing.

The approved certificate and all six original captures remain historical. Its proof is an exact factor-preserving, signed affine structural replay over qualified historical/current producer captures. Only the D65 typed admission helper ran inside both pinned kernels. This migration makes no private-cache or general grouping claim.

== Exact literal map

#table(
  columns: (2fr, 1fr, 1fr),
  table.header([Test], [Old digest], [Reviewed native digest]),
  [`comparison_frames_with_fermion_flow_reversals`], [`0xbbb6f06b55b1dd03`], [`0xf2c181dfbe49930f`],
  [`representative_frames_with_fermion_flow_reversals`], [`0xbd1d483fc596a496`], [`0x196ad47f7c361abc`],
  [`representative_frames_with_left_right_symmetrization`], [`0x87988a044c993bcf`], [`0x9a59879459492f0e`],
)

The reviewed native decimal values are 17492405225653900047, 1831509842489383612 and 11122069825991028494. Exact Python integer checks join these values to the old literals, both approved reviews and frozen diagnosis.

== Hygiene and test evidence

A separate hygiene change `8c4011890ba7cef3edc8be5d17ce3cc893c25aa6` formats the twelve proof scripts with pinned Ruff 0.16.6 and the diagnostic Rust files with pinned stable rustfmt. The Python ASTs were identical immediately after formatting. Two explicit lint rewrites bind an immediately used lambda's loop variables and replace generator-to-set expressions; an intentional diagnostic `exec` has a scoped lint annotation. Eleven substantive proof/corruption scripts replayed into `target/frame-certificate/migration-r1/formatted-replay/`; every mathematical summary equals the frozen approved summary. The approved certificate JSON and raw captures were not rewritten.

The hygiene binary with old literals ran the 1/8/32 thread sequence and failed only the three numerator assertions. Its observed actual values equal the reviewed native digests, with every other summary field equal. The new qualified `ci-optim`, default plus native oracle binary passes all three named tests separately at 1, 8 and 32 threads, with record, dump and force-pass modes unset. The test output confirms three executions and zero failures at each setting. The current D65 typed admission and effective-model fingerprint diagnostic tests also pass on that binary. Scoped Clippy with warnings denied, pinned Ruff lint and format, and pinned rustfmt checks pass.

The structured `docs/architecture/evidence/generation-frame-migration/receipt.json` records exact source and binary SHA-256 identities, logs, ref resolutions, row values, and preserved artifact inventory. It does not hash itself. The original 2102 pair captures, 18 model-hook captures and 64 non-script approved artifact hashes were verified after migration; the twelve original script bytes remain available at the approved commit, while their formatted successors have separate hashes. No broad docs or unrelated physics build is claimed.
