= Three one-diagram QED oracle numerator diagnosis correction

== Correction scope

This round-two record repairs the certificate and provenance defects found by two
new independent Opus reviews of round-one commit `e1d44235`. The round-one
17-file candidate, its raw JSON bytes, receipt, manifest, validator, report,
ELFs, logs and review directories remain unchanged. The correction still covers
only `xs_qed_eedd_L1`, `xs_qed_eemm_L1`, and `xs_qed_eeee_1pi_L1` at the
single-worker pool. It changes no production algorithm, expectation, request, generation option,
lock, kernel, manifest or NixCI setting. It does not extend to the other 78 rows
or pools 8 and 32.

The physical conclusion remains bounded and unchanged: no mathematical
numerator defect was found in these three factorized graphs. Their historical
and current capture digests differ at a representation/export boundary. This is
not a general grouping theorem, a proof for another model family, or authority
to migrate snapshots.

== Asserted route and fresh provenance

The row asserted by `main.rs::check` is the particle-only reference returned by
`classes::check`: `Classes::None` at the first requested worker count.
`Classes::Auto` is generated afterwards and compared with that reference. The
round-one report incorrectly called the asserted row Auto. Fresh pristine runs
now retain both `<case>.classes-None-t1.json` and
`<case>.classes-Auto-t1.json`; they are byte-identical for every route and case.
The Auto live-metadata capture is therefore bridged exactly to the asserted None
payload, rather than assumed to be the asserting execution.

Opus provenance review A ran, in each locked source workspace,
`FEYNKIT_ORACLE_THREADS=1 FEYNKIT_ORACLE_DUMP=<review dump>` with record mode
unset and `nix develop -c cargo test --locked --profile dev-optim -p
feynkit-generator --test generation_oracle -- --test-threads 1 --nocapture` plus
the three exact names. Historical source `2b1ef4f0`, Symbolica `c78cffc2`, and
ELF SHA-256 `32d17aaa...0695` pass three of three. Current source `b9f78892`,
Symbolica `1a1829b1`, and ELF `c2ef6aea...be84` fail all three only at the
unchanged numerator assertions. Runs without `FEYNKIT_ORACLE_DUMP` reproduce the
same results, establishing that dump capture does not record or force-pass.
Both shells use Rust 1.98.1 and the native default features.

The round-one raw historical/current payloads are byte-identical to the newly
retained fresh None and Auto dumps. The corresponding `.census` files, exact
environment-bearing logs, producer paths and hashes are now in
#link("evidence/generation-oracle-numerator-diagnosis/round2/provenance")[round-two provenance].
The original round-one dump directories were not retained and are not claimed
to exist. The pristine ELFs saved in round one came from a later post-restore
rebuild, not the earlier assertion process; the fresh Opus executions ran those
exact bytes without recompilation and now qualify them.

Earlier nextest attempts were prebuild attempts without `--locked`; the current
one exited 100 and neither establishes an assertion result nor a worker/dump
environment. The adapted diagnostic went through two compile failures, an r2
run and the final decoded run. Only the final patch survived round one; no
original adapted ELF was retained. Review A independently reapplied the archived
patches to clean source exports and reproduced all six metadata files byte for
byte. Its retained adapted ELFs are historical `31e9ea39...b01ef` and current
`1a76aa51...66c50`. The correction preserves the patched sources, exact logs,
features and producer hashes and does not reconstruct the lost intermediate
ELFs.

The diagnostic did not serialize `Request` or `Run`. Equality of those inputs
is consequently qualified as source-level identity: the bound case definitions,
`cross_section`, `plain`, `loop_one_pi`, particle constants and request-building
functions are unchanged for these three definitions; the captured full QED
model, numerator options and grouping settings agree, and the live reports give
17/34, 17/34 and 3/4. The command and t1 filenames bind `Run` to one worker.

== Complete semantic boundary

Capture identity and semantic equivalence are separate entrypoints. The capture
checker verifies the immutable manifests, fresh None/Auto bridge, six exact raw
payload hashes and the six historical/current FNV numerator digests. It gives no
semantic credit merely because a fixed digest matches.

The semantic checker parses every expression without expansion. It admits only
the complete symbol-head inventories frozen from the qualified captures and a
closed Q/K/P/edge/vertex mapping bound to the authoritative historical and
current `symbols.rs` hashes. Only Q and edge occur in these payloads. Mapping is
on whole parsed heads; prefixes, foreign namespaces, trailing syntax and
unsupported scalar grammar reject.

For each side it requires one exact diagram ID and a bijection between payload
vertices, edges, aggregate numerator, prefactor and projector and the complete
live owner table. Current live expressions must equal their payload fragments.
Historical live expressions must decode to their payload fragments through an
injective map of the actual 14 cooked symbols per row. The metadata ports were
decoded after generation from live `logical_slots` using each producer's actual
`CookSettings`; they are observations, not native admission metadata.

Every interface-like field is rejected outside its exact known owner. A present
interface must be an array equal in order to that owner's live ports before it
can be omitted from the cross-route payload comparison. Missing interfaces are
allowed only at owners whose exact producer-observed inferred order is frozen in
the semantic contract. Null, malformed, missing, extra, wrong-owner,
dimension/variance and mirrored corruptions reject. Metric/composite flags and
historical logical order are likewise bound as observed native producer outputs,
not promoted to independent theorems.

The checker derives every free typed port from the parsed factorized expression,
retaining owner, representation, dimension and variance, and validates graph
incidence. There are exactly 92 axes. The six nonidentity historical-to-current
transports are eedd `vertex:0` and `vertex:1` `[4,3,0,1,2]`, eedd `edge:2` and
`edge:3` `[3,2,0,1]`, and eeee-1pi `vertex:1` and `vertex:2` `[1,0,2]`.
Consistent explicit-interface transposes remain semantically valid even though
they correctly change capture bytes and FNV.

Without distributing a numerator, the aggregate factor multiset and exact
Gaussian-integer coefficient equal the product of every local vertex and edge
factor. The prefactor, numerator and projector form a completely contracted
scalar. Both routes' sole group has master 0, one member, source index 0 and
ratio 1; diagram/member IDs, names, source factors and payload/diagram overall
factors bind exactly, with scalar syntax parsed or rejected. No trace,
whole-graph expansion, parameter-at-one evaluation or numerical sample is used.

== Controls, evidence and limits

The original review scripts and results remain unchanged under `review-inputs`:
provenance review SHA-256 `747af1ff...4c9b8`, independent checker
`885eb2bd...b7e1`, mutation driver `aa27ced2...f9a7`, results
`e8618262...e3159`, and independent FNV source `099bae69...3395`. The corrected
full semantic entrypoint exercises all 43 records actually present in the
retained review matrix, including valid identity, key order, owner order and
consistent transpose controls, plus six held-outs for a mirrored foreign head,
missing-interface order, mirrored model, unknown historical interface, mirrored
diagram ID and a second valid transpose. All 49 outcomes match their intent
without SHA or fixed-digest shortcuts.

The #link("evidence/generation-oracle-numerator-diagnosis/round2/receipt.json")[correction receipt],
#link("evidence/generation-oracle-numerator-diagnosis/round2/validate.mjs")[combined validator],
separate capture and semantic checkers, semantic contract, control results,
qualified provenance and inventories form the correction bundle. Direct Typst
PDF compilation and semantic HTML bundle compilation are source- and hash-bound
there. In tracked scope the full documentation builder remains blocked by the
inherited eight mark-geometry HTML inventory entries. In this working workspace
it additionally reports the retained untracked Opus review HTML; that evidence
was not deleted to manufacture a pass.

No expectation migration, integration, publication, upload, push, fetch or CI
activation occurred. The next comparison outside this scope remains a new
full-key/local-owner typed-axis transport plus exact group master, member source
factor and ratio binding for one separately selected model family or grouped
case. This correction stops for new independent Opus A/B reviews, with Astra as
fallback if Opus is unavailable.
