= FeynKit community refresh and shared qualification

This investigation records the non-destructive reconciliation of approved
integration #raw("fdff0d18953348b84a0f1b7e6ab5ef703a5090de") from upstream
#raw("3f0df2a870806076820786994d4883decec8e9e3") onto exact upstream
#raw("9102f885722e856e30227ed46e685e5c465905f7"). The kernel changes from
Symbolica main #raw("1a1829b1") to community
#raw("8eb7c3729d6abefc916decc1bdf0c06536f8af80"). Root, Tydenso and notebook
locks retain one shared Symbolica/Numerica/Graphica 3.0.1 identity; the uv,
Hakari and workspace-hack inputs stay aligned. This is a source reconciliation
and bounded qualification record, not integration approval.

The account originally added at #raw("c3449988839d9ad3389dd9fd4c4b8bb9ec447055")
is preserved in that commit. The additive review aggregation below updates its
pending-review language only within each reviewer's actual scope. This updated
note and receipt remain an unreviewed metadata candidate. The registry uses
#raw("evidence_revision") at #raw("cacd5182704082df929a066548e3501be22c9533"),
without #raw("review_ref") or #raw("reviewed_at"): XS-1's paired-field defect is
repaired without inventing approval of the new note.

== Source history and corrections

The original 106 reviewed changes and their evidence remain immutable. A
separate 106-change replay onto the new upstream is conflict-free. Round-one
code #raw("5ea149ffd88ae3afc4072a38c0871c1329be8424") and evidence
#raw("0afe4d3d940e8b181d97243cf78508e854bdef6f") remain named checkpoints.
The first independent source review requested changes. A clean final tree and
clean replay history did not establish semantic compatibility: textually clean
upstream hunks still contradicted approved model and API owners.

The explicit five-file follow-up
#raw("f8422fcce0b3e9580fcb05ef6aa908b93b70521e") restores the canonical UUV1
sum, three ghost crossings and negative-PDG branch while preserving upstream
API conversions and the scalar-prefactor reuse. It migrates the upstream
squared-amplitude display test to the existing staged Request API with the
same typed V_98 restriction, tree level and two-vertex bound. Its 21 assertion
ASTs and display suffix are unchanged. Both sector-decomposition examples
return to the staged API; the linked ExtensionPoint body and tests remain
unchanged, and the existing stubs already contain these examples.

The browser test derives natural and fitted scales from measured drawing
dimensions, using the current renderer contract. Percentage rounding matches
JavaScript's positive half ties; label ratios and zoom follow those scales.
The approved pinned drawings, stable node locator, tolerances and assertion
coverage remain. Added regime guards still require actual browser execution;
the worker's seven arithmetic regimes are source checks only.

The next bounded follow-up
#raw("5ba04864ca506e259d9eef4a54f65acdd8323fc2") addresses two representation
boundaries identified by the second source review. It converts both authored
axes to float32 before far viewBox containment checks, retaining the 1e-6
margin. For the tall outer SVG only, root expressly authorized equality within
one CSS layout unit (1/64 pixel), because the new renderer removed toolbar
slack. CTM delta=1e-4, label places=2, origins, drawings and all other controls
remain. These are source-derived precision allowances, not recorded browser
outputs or a universal engine proof. The
#link("evidence/generation-integration/feynkit-refresh/browser-precision/receipt.json")[precision receipt]
binds the one-file delta and root's requested fresh focused source review.
The prior CODEf842/EVIDENCE594c and the independent reviews stay unchanged.

The additive
#link("evidence/generation-integration/feynkit-refresh/corrections/receipt.json")[correction receipt]
binds the review, worker patches and exact five-file delta. At c344 a new
independent source review was still required; source R4 below supplies that
bounded verdict. No runtime, installed-host, browser or generated-stub pass
follows from these source checks. The original handoff prose
#raw("reconciliation.txt") remains historical round-one evidence. Its incomplete
semantic judgments are superseded by this correction record, not silently
rewritten. The old hash handoff combines CODE5ea and EVIDENCE0afe inputs; the
additive receipt explicitly names both identities. Original replay52's crossing
prefactor conversion was already present upstream, independently of the ghost
expectation corrections.

== Outer SVG normalization and source R4

Code #raw("f47f27213a1c55291239c9bc0f334a7ea1cd1cc1") implements the selected
Option A correction in the browser test alone. Camera CTM and node-label
dimensions are divided by the validated positive uniform outer SVG screen
scale. Independent viewport/viewBox, meet-matrix, one-layout-unit and
source-derived gamma5 bounds constrain that normalization; resize rereads the
outer scale. Camera delta 1e-4, label places 2, binary32 authored containment
1e-6 and the R2 tall allowance 1/64 remain unchanged. Additive evidence
#raw("cacd5182704082df929a066548e3501be22c9533") retains the exact worker
rationale and explicitly identifies the root-authorized comment-only departure.

Independent Opus source R4 approves source reconciliation through f47f and
the truthfulness of cacd's evidence. It expressly excludes c344's shared
metadata; its approval cannot be extended to this new note. The reviewer
executed the exact final test method with static current-JS/SVG fixtures in
Chromium 153, including near-maximal layout residuals, and seven negative
controls. This is isolated fixture evidence, not an installed product-renderer
or WebKit pass. Worker headful DPR 2 snapped at 1/128 CSS pixel; the reviewer
also exercised the DPR 1, 1/64 regime. Crop and label control attribution was
checked in the actual assertion order.

R4-N3 is an #emph[executed] pre-existing residual: Chromium's preserved tall
assertion fails for 48.9831 by 620.3672 pt by 1.67e-5 pixel beyond its 1/64
allowance, while the new helper's representation-aware height check passes.
The reviewer classifies it as nonblocking because R2 was explicitly preserved.
The model's 21 failures in 300,000 random tall geometries are simulation,
not additional browser executions. Installed-host qualification remains open;
WebKit parsing is untested and Gecko's 1/60 unit lies outside this allowance.
R4's formatting and verbatim-rationale Typst nits remain disclosed in the
copied review; the historical rationale is not imported into this note.

== One shared producer boundary

The one authorized pure build selects the package-specific
#raw("gammaloop-nextest-feynkit-feynkit-generator.tar.zst") producer, not the
all-package link farm. Its immutable CODE5ea Git URI has an explicit full
revision, shallow and allRefs parameters. The fetched tree matches all 4580
tracked files and symlinks. The archive and package-test sources contain 670
and 669 byte-identical selected source files. The recorded recursive chain
includes dependency-filtered sources, synthetic dev-dependency manifests and
actual phase/hook commands. The selected native compiler is rustc 1.98.1,
commit #raw("48a229ceaefd4985c50990b14116b6d856af0985"), LLVM 22.1.8.
The archive profile is #raw("ci-optim"), distinct from old development observers.

The exact corrected CODEf842 and CODE5ba probes resolve to the identical package-specific
archive derivation, with unchanged source, dependency artifacts, synthetic
inputs, phases, profile and Rust flags. All 47 immutable filtered source inputs
remain the same. The broader FeynKit group changes because other package tests
and examples changed; that group is not built. The metadata-only round-one
EVIDENCE0afe probe also reproduces its producer. These are explicit producer
invariance checks, not claims that every source revision was independently
runtime-built or that the corrected Python/browser paths were executed.

Nix may fetch locked public dependencies or use normal substitution during the
authorized pure build. This is not an offline claim. Every Nix child receives
neither license variable and uses a validated immutable Git URI or exact store
path from a neutral working directory. No private path flake, licensed Nix
result, cache transplant, unrelated dependency upgrade or license/budget policy
is introduced. The source build request remains CODE5ea even when identical
producer inputs permit bounded transfer to CODEf842 and CODE5ba.

== Helpers and runtime limits

The canonical native #raw("oracle_family_exact.rs") preserves its pair/result
contract and sign, coefficient and power controls. Its negative controls reject
zero and idempotent units; existing additive-witness protocols remain necessary
for those faithfulness cases. Atom admission alone is not a graph-membership,
definition, schema or physical-equivalence certificate.

The retained public capture entry point supports exactly its four named
ungrouped cases at pool1 for Classes None and Auto. The selector-capable copied
harness observes its existing named raw build stage; the observational library
module is a separate source owner, not an interchangeable standalone program.
The private diagnostic source manifest retains original inputs and the exact
observation-only overlay. Helpers must bind actual compiler arguments and
qualified current8eb dependency units. No old1a objects, moving grouped adapter,
production manifest change or fabricated universal observer is permitted.

Python validators and their local import/dynamic-loader sources are additive
lane proof dependencies, outside Nix and Rust compilation inputs. Merely
hashing or compiling adjacent sources does not execute their proofs. Existing
replay scripts still name old-base executables until explicitly rebound.

Runtime characterization is limited to the eight previously migrated rows and
four new raw rows, pool1 only, with exact nonempty selections and record/update
flags absent. Both class modes require actual public payload, census and source
companions. Observational public dumps must join the pristine producer before
semantic use. All actual passes, failures and full-Row changes must be retained;
no eight-pass/four-fail outcome is presumed and no literal is re-recorded.
A minimal-mount unlicensed preflight verifies a distinct network namespace.
Private runtime license delivery occurs only through the sandbox's anonymous
pipe followed by direct execution, with no network or Nix operation in that
licensed runtime boundary.

The sole pure archive build completed successfully. Its archive SHA256 is
#raw("9688999459ceb9e3747f35e5ed599437fbdf64ef03a4564ca537b3254c0dac1e").
The layered compiled test closure resolves 452 actual Cargo dependency units
through three separately hashed artifact layers. Actual archive Cargo metadata
resolves all three kernel crates to community8eb7. Compiler/profile/features,
synthetic inputs, phase commands, NAR metadata and actual helper Rustc arguments
are retained; no working-copy cargo tree substitutes for executed provenance.
All four helper artifacts compiled once using those qualified8eb units. The
native admission executable's SHA256 is
#raw("5141d8fddacf40a81f5452e212cd788daa15c3cb3e6d2290fbfc6b65f76c48e1").
The bounded public observer, private observational library and real-harness
observer have separate artifact and command identities; no ELF identity is
invented. The archive ELF and layered test ELF have distinct bytes: only their
loader sections differ, reflecting the archive's added own-output RUNPATH.
All non-loader section contents and decoded dynamic symbols agree. Both ELF
identities remain bound, and no binary is patched or relinked for transfer.

All twelve exact pool1 characterizations completed: the eight previously
migrated rows passed their existing assertions, and the four new raw rows
failed with only the asserted numerator digest differing in each full Row.
No literal or expected output changed. Each pristine run selected exactly one
test with exit0 or101 as actually observed. The observed harness reproduced
that result and byte-joined every public payload, census and source companion.
Both actual None/Auto cells retain their own model/request/options/Run,
typed-owner/port ledger and three existing raw observation stages. The table
assertion uses the harness-returned particle/None reference; Auto counts and
payloads remain separately captured, never copied from None.

The
#link("evidence/generation-integration/feynkit-refresh/qualification/capture-index.json")[capture index]
and
#link("evidence/generation-integration/feynkit-refresh/qualification/full-row-characterization.json")[full-Row characterization]
point to actual completed records without altering them. These are current8eb
runtime characterization and source joins only. Independent lane schema,
definition, source-factor, exact-pair, mode-correspondence and adversarial
replays were still pending at c344; the separately bounded reviews below now
record their actual dispositions. Migration approval remains separate.
Old-capture admission on the new checker is not a fresh producer certificate.
Those c344 characterizations added no extra rows, pools8/32, grouped adapter,
installed-host or browser execution.

The c344 metadata-only documentation/evidence head was not separately
runtime-built. The actual build source remains CODE5ea, with explicit identical
selected-producer transfer to CODE5ba; the new metadata files are outside the
recorded generator input set. Source R4 also read cacd's recorded selected-
producer invariance through f47f; it did not re-evaluate or rebuild it.
Full installed-host, generated-stub parity,
Chromium/WebKit including regime guards, browser/widget and full-site checks
remain pending. The 452-unit build does not establish whole-workspace kernel
compatibility. Source review, arithmetic controls, compilation, licensed
characterization and physical proof are separate evidence categories.

== Preserved evidence and review boundary

All 1204 historical evidence files and 243 oracle/test_feyngen paths remain
byte-identical. Generator, model, amplitude, Python generation bridge and
physical fixture owners are unchanged from approved fdff; all eight authorized
literal migrations remain untouched. The 118 scope keys remain ordered and
valid. The two earlier Nix scope bindings describe upstream's explicit iframe
and compile-time test delivery plus widget browser execution: added source and
runtime coverage, not a no-behavior-change claim or inherited browser pass.
The correction itself changes no covered region or digest.

The task's initial inventory scratch commit #raw("883184320f00") is retained at
#raw("generation-feynkit-refresh-inventory-scratch"). Earlier attempts, source
copies, failed preparation/selection attempts, build logs and qualification
outputs remain available for review. No shared cache or other contributor's
work is removed. The physical main index guard remains unchanged, the approved
integration reference stays at fdff, and top-level NixCI enable remains false.
This record and its registry note are a separate unreviewed evidence candidate.
A new combined review remains necessary before integration approval.

== Separately bounded review dispositions

The
#link("evidence/generation-integration/feynkit-refresh/review-aggregation-20261030/receipt.json")[review receipt]
and its
#link("evidence/generation-integration/feynkit-refresh/review-aggregation-20261030/receipt.typ")[scope record]
bind six verbatim independent review JSON records, their original hashes,
source locations and available manifests. Referenced scratch, captures, helper
sources and native results remain retained in their original locations; they
are not all vendored here. Earlier reviews, attempts and evidence files remain
unchanged. No review here approves the conjunction of all lanes.

- #strong[Shared Scope A:] Opus approves the twelve-case pool1 provenance,
  runtime and evidence chain and the shared provenance for both new raw lanes.
  It excludes source reconciliation, mathematics, Scope B, migration and
  combined approval. Its nonblockers retain the historical shared-target and
  crates.io Graphica qualification, lane-internal QED input-version history,
  diagnostic comment/unused-helper nits, and the pre-existing Git-colocated
  integration ref at 863f506c; the authoritative jj bookmark remains fdff.
- #strong[Prior eight:] Opus approves finite application of the three already
  approved certificates to fresh 8eb pool1, Grouping None, Classes None and
  Auto independently (16 cells). It does not re-review the original physics or
  approve new literals, grouped modes or pools 8/32.
- #strong[Four new raw rows:] Opus approves corrected Scope B preparations
  only. QED's four newly explicit interfaces per cell and QCD's 32 permuted
  explicit interfaces per payload/mode remain capture-bound and uncertified;
  migration must explicitly accept the representation changes. SM has zero
  explicit interfaces on both routes. Raw-four wrapper findings R1–R8 remain
  attributed to that review; separate errata work at
  #raw("target/raw-four-wrapper-errata-20261030-r5") in the numerator-families
  workspace is pending, with no general-runner or native-adapter approval here.
- #strong[Grouped None provenance:] the review remains
  #raw("REQUEST_CHANGES"). Lossy JSON DiagramIds, observer side effects and
  overstated producer-validator coverage require additive disposition. The
  correction packet
  #raw(".local/feynkit-refresh/grouped-none-provenance-disposition-20261030-r2")
  remains pending and unapproved. Shared Scope A approval covers its shared
  dependencies only; it does not discharge this rejection.
- #strong[Grouped None mathematics:] Opus approves the finite capture-premised
  decoded-numerator certificate for #raw("xs_qed_eemm_identical_L2"), NamingD,
  Classes None, Grouping Identical with the captured five-sample/seed-3 options,
  pool1, proof source 5ba, dependency archive 5ea and kernel 8eb. It proves the
  signed sewing, D22/D23 relation and observed 30-to-29 collapse within those
  premises. Aggregate certification still requires grouped provenance approval
  and root aggregation. Raw-cook/decoded identity, ordered interfaces and
  historical FinalFrame provenance remain qualified. S1–S7, including the
  separate one-shot transport/spinor checks, chart-invertibility assertion gap
  and unsnapshotted helper sources, remain disclosed. Auto and migration remain
  unproved and unapproved.

== Prior-eight reporting dispositions

F1–F7 reporting is closed additively here, using the independent prior-eight
Opus review's checks and attribution; this metadata change reruns no proof or
control and changes no sealed packet.

- #strong[F1:] the historical c78 kernel leg carries by exact pair/input
  invariance; the old 1a current-kernel leg is superseded and gets no fresh
  credit; the 8eb leg is fresh in the reviewed packet and reviewer replay.
- #strong[F2:] L2's constant #raw("None") maps asserted live #raw("none");
  #raw("built_len") is the length of C's raw built-diagram list, equal to the
  old value. This is an explicit presentation map.
- #strong[F3:] the reviewer verified fresh QED None/Auto and old-current dump
  byte equality, matching source companions and completed D0 census. FNV
  equality follows from byte identity; the frozen outer validator still rejects
  the changed registry and receives no new pass.
- #strong[F4:] with the reserved additive witness, sign/coefficient/power
  booleans may be satisfied by the witness alone. Mathematical content is
  parsed/canonical equality with additive cancellation, not five independent
  mathematical checks.
- #strong[F5:] the original adapter used old typed-owner order. The reviewer
  independently reproduced it with a C-only fixed tail rotation and accepted
  raw C order with identical obligation multisets.
- #strong[F6:] the original sandbox exposed the whole runtime root writable,
  including command/stdout/stderr records. Independent reproduction supports
  the results; the original containment caveat remains.
- #strong[F7:] the inherited #raw("anchor_self_comparisons") label now counts
  C-versus-B cross-capture anchors; it remains unchanged in the sealed outputs.

The reviewer carries 428 pure rejections and 115 QED controls only by complete
source/input invariance. The 194 historical-only control records are not
upgraded to current credit. The reviewer executed named held-outs; this change
claims no fresh full-validator mutant suite. All decoded/capture, definition,
ordered-interface and old QED Request/Run source-constructor limitations remain.
No full documentation, runtime, installed-host, browser, migration, combined,
WP3 or integration gate is closed by this metadata aggregation.
