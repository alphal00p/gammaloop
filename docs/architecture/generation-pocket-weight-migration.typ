= Photon to top-pair QCD=2 single-row snapshot migration

== MAIN39–43 integration

This approved preparation record is imported into the ordered MAIN39–43 candidate.
The original source, proof and test bindings below remain historical; current combined
checks and approvals are recorded in #link("evidence/generation-integration/wp3/main39-43/receipt.json")[the integration receipt].
A NEW independent Astra combined review is required before44.

== Historical preparation record

This prepared migration changes one inline expectation in `failing::test_generate_aa_ttx_cross_section_slow_filter`: `874 | -266 = -266` becomes `874 | -186 = -186`. The test name, module, imported `sm.json`, process command and options, clean-state lifecycle, QCD=0/1 expectations, and every other assertion are unchanged. The separate `slow::test_generate_aa_ttx_cross_section` row is untouched. No production generator or physics code changes.

The exact old-fail control uses the retained `ci-optim`/default-feature binary from approved certificate `5f7d96ba675f76dcb7e814e085c5f66ece13cb2b`. It passes the first two inline snapshots and fails the named QCD=2 row: expected `-266`, produced `-186`. The newly compiled binary from the one-row edit runs the same full module-qualified test with snapshot update modes unset. One test runs and passes all three cells: QCD=0 `4 | -4`, QCD=1 `40 | -40`, and QCD=2 `874 | -186`. The #link("evidence/generation-integration/wp3/pocket-weight-migration/migration.json")[migration receipt] binds both source and binary hashes, exact commands, logs, model and lock inputs, and the one-line diff.

The #link("evidence/generation-integration/wp3/pocket-weight-diagnosis/diagnosis.json")[frozen diagnosis] received independent Astra scope A approval at `5f7d96ba…`, recorded in `target/pocket-weight-diagnosis/review-astra-scope-a-round2-JOlXh9/review.json`. Independent scope B approved every retained graph, cut, automorphism, Grassmann sign, and weight at the earlier `a8c5a11e…` candidate, recorded in `target/pocket-weight-diagnosis/review-astra-scope-b-YmoQie/review.json`. Scope B's raw bindings were rechecked unchanged before this migration. Its 40-ghost-sign counterfactual reaches `-266`, but does not establish why the historical literal was written. The first observed incompatibility remains the named snapshot versus its own pristine producer.

The named failing row remains `-266` through `087bda82ce5a14dab092776ea470e66278dda12d`, which changes the distinct `slow::test_generate_aa_ttx_cross_section` QCD=2 row. `dfe6387606e44a96f110c1ce03f52d92824dbdb1` later changes the named failing row to `-186`. These are source-history observations, and no v4 execution is claimed.

Before the expectation edit, a separate import-only test-helper change replaced qualified `File`, `Value`, and `json!` use sites with module-scope imports. Its isolated diagnostic rerun reproduced all 874 filtered and 1,106 unfiltered DOT files, graph TSVs, resolved process definitions, and live model/parameter captures byte for byte. This cleanup has no numerical effect and leaves the approved certificate artifacts frozen.

This is only the named aa to top-pair QCD=2 row. It does not resolve other source44/45 snapshots, older frame certificates, drawing38, or the larger ordered generation integration. A fresh independent Astra migration/test review must approve this concrete edit before any gate recertification. No main integration, upload, push, publication, or CI enablement is authorized by this record.
