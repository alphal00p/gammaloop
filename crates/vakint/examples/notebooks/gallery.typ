= Interactive Vakint field notes

These five Marimo applications share the native HEPKit, tensor, graph and
Symbolica types from one coherent host. They demonstrate graph dressing and
efficient gamma/color contractions, precision-aware constant/master inspection,
an explore-only view of the zero-momentum four-rung (five-loop) vacuum fixture,
exact numerator/pinch closure, and end-to-end Ordered reductions of the three-
and four-loop ladders in cancellable background jobs.
The source follows the native workflow of the
#link("https://www.hepkit.org/gallery/notebooks/four_loop_numerator.py")[HEPKit numerator tutorial].

== Run the gallery

Use the Python environment containing the coordinated Community host, Marimo
0.25.2, Uvicorn and the host's native display dependencies. Do not install a
second Symbolica wheel into that environment. From this directory run:

#raw("python serve_gallery.py", block: true)

The server binds only to #raw("127.0.0.1:8417"). Configure
#raw("VAKINT_GALLERY_PORT") to select another local port. A persistent user
service and tailnet-only HTTPS proxy can expose that port; the notebook process
does not change the tailnet or existing services itself. The five routes are
#raw("/graphs/"), #raw("/inspect/"), #raw("/five-loop/"), #raw("/closure/") and
#raw("/reduce/").
The server uses application mode, with source editing disabled.

The dedicated tailnet deployment is
#link("https://itphlies.tailb3264.ts.net:8453/")[itphlies.tailb3264.ts.net:8453].
Its user service is #raw("vakint-five-loop-gallery.service"); restart or stop it
with #raw("systemctl --user restart vakint-five-loop-gallery.service") or
#raw("systemctl --user stop vakint-five-loop-gallery.service"). Keep the final
deployment's sealed operations receipt for its interpreter, logs and rollback
configuration. Restarting the app does not reconfigure Tailscale.

Opening a page does not load five-loop owners or start a physical reduction.
The inspection app constructs the offline vendored 1L–3L engine. Its explicit
load button (1L–4L) validates #raw("VAKINT_GALLERY_ARTIFACTS") when
set; otherwise the official latest packs install once into
#raw("VAKINT_GALLERY_DATA_CACHE") (default: the user cache directory,
#raw("~/.cache/vakint-gallery")), never into the server's working directory. A loaded engine
is reused for that session's subsequent lookups. The 1L–5L catalogue is offered
only when the operator sets #raw("VAKINT_GALLERY_ALLOW_5L=1"): its load runs in
the page kernel for minutes, holds several GB per session and cannot be
cancelled. Every master remains a native
snapshot-bound inspection object, and #raw("master.graph") is an ordinary
HEPKit graph. The 5L catalogue contains declared, uncomputed terminals with empty
coefficient maps, never fabricated zeros.

== Recorded evidence

Each optional producer receipt is listed in #raw("data/manifest.json") with its
relative path and SHA-256. #raw("VAKINT_GALLERY_EVIDENCE") can select another
directory containing the same manifest. The app checks the bytes before display
and rejects traversal or corruption. Receipt names are #raw("campaign"),
#raw("zero-momentum"), #raw("closure"), #raw("cache-build"),
#raw("held-out") and #raw("warning-controls"). The last records a real tuning
case which triggered a lowered-entry coverage advisory and still reduced to
declared catalogue masters. Its provisional owner mapping does not establish
global physical coverage, a speedup, or an example of unknown residuals.
An absent receipt makes no timing or successful-reduction
claim: the five-loop and closure pages omit cards whose receipt is absent.
The producer remains responsible for scientific validation and source/data provenance; a local hash does not certify those facts.

The #raw("cache-build") and #raw("held-out") receipts, and the older
#raw("zero-momentum-running-cache"), #raw("-light-cache") and
#raw("-heavy-cache") ones, describe the removed depth-first completed-result
caches. They are historical: Vakint no longer builds or reads those caches, and
no page shows them unless an evidence directory provides them.

Browser previews have bounded depth, item counts and string lengths. A receipt
above 8 MiB is offered as a download instead of being parsed into the page.
Complete guard lists and unknown-integral identities remain in their original
scientific artifacts. Each receipt entry can list additional
#raw("artifacts: [{path, sha256}]") below the same evidence root. The
#raw("/evidence/") routes verify these listed hashes and stream downloads without
embedding large arrays in the browser. Unlisted paths and checksum mismatches
are rejected. A returned expression with unknown terminals is labeled a partial
retained result, never a completed reduction to catalogue masters.

The primary physical fixture is the complete zero-momentum contraction with both
current vertices, unaveraged #raw("g_mu_nu delta_a_b"), massless quark slash
numerators, and the auxiliary mass only in vacuum denominators. It has
#raw("R=5, A=14, D=9"). The distinct quadratic fixture and its transverse
normalization remain deferred; no completed or cached quadratic result is claimed.
Timings distinguish numerator preparation, loading,
routing, IBP reduction, cache attachment and total runtime. Final comparisons
use medians/ranges of three fresh processes and exact result/guard checks.

== Background jobs: bounded and cancellable

Every reduction started from a page is a background job, never an
in-kernel computation. #raw("gallery_jobs.py") starts the supervisor
#raw("benchmark_job.py") in a new session, so a job is exactly one process
group. The supervisor runs the job's command with
#raw("VAKINT_N_CORES") and #raw("RAYON_NUM_THREADS") set to the job's
#raw("n_cores"), pinned to that many of the server's CPUs, at niceness 10, and
with a wall-clock deadline. It hands the command the work root's lock, so one
job runs at a time and the lock is released only when every process of the job
has exited. The command dies with its supervisor (#raw("PR_SET_PDEATHSIG")):
no worker outlives its job.

Cancel sends SIGTERM to the job's process group and SIGKILL after five seconds,
then confirms that no member survives; recorded #raw("/proc") start times guard
every signal against PID reuse, and the server's own process group is never
signalled. The supervisor records #raw("cancelled"), #raw("timed_out"),
#raw("failed") or #raw("complete"); a supervisor that vanished without a final
status is reported as #raw("lost"). Jobs live under #raw("VAKINT_GALLERY_WORK")
(default: a per-user directory under the temporary directory); a fresh page
re-attaches to the running job and can cancel it.

Server limits: #raw("VAKINT_GALLERY_MAX_CORES") caps every engine's and job's
#raw("n_cores") (default 4, never more than the server's CPUs), and
#raw("VAKINT_GALLERY_MAX_SECONDS") caps a job's deadline (default 1800 s).
#raw("VAKINT_GALLERY_N_CORES") and #raw("VAKINT_GALLERY_ARITHMETIC")
(#raw("exact"), #raw("modular_sampled") or the default
#raw("modular_polynomial")) select the in-page engines' defaults.

== Three- and four-loop reductions (#raw("/reduce/"))

The page reduces the massless-quark QCD ladder self-energy with two or three
gluon rungs to retained masters with the Ordered engine (MSbar, dot-product
notation, L+1 epsilon orders, as in #raw("tests/ordered_reduction_rungs.rs")).
The prepared scalar inputs ship in #raw("data/ladder/") with
#raw("cases.json"): 3L zero momentum and quadratic, 4L zero momentum. Each
input's sha256 equals the rung reference's #raw("source_input_sha256"); the
worker #raw("reduce_worker.py") checks it, binds the dimension, checks the
recorded evaluated-input sha256, reduces, and compares the retained result's
sha256 with the recorded fixture (identical for every arithmetic and core
count). The 4L cases first load the 4L external-data pack. The 4L quadratic
input (1.7 MB, about 25 s of reduction) is not shipped; a server provides it
with #raw("VAKINT_GALLERY_LADDER_INPUTS"). An optional persistent Ordered step
cache under the work root makes reruns warm.

The page keeps showing its last run, finished or not, when a control changes;
the URL carries #raw("?job=<id>"), so a reload shows the same run, and the
Recent jobs list links every finished ladder run of the server back to its
result.

Run a case headlessly with the page itself, which then waits for the job:

#raw("python 05_ladder_reductions.py --case 3L-zero-momentum --n-cores 2", block: true)

(optional #raw("--arithmetic"), #raw("--max-seconds") and
#raw("--cancel-after SECONDS"), which cancels exactly like the button). It
prints one JSON line: the job directory, its state, the retained sha256 with
#raw("matches_recorded"), the catalogue closure, the seconds and #raw("ok"). The
exit status is 0 only if the run completed with the recorded fixture's sha256
and a closed result, or was cancelled as #raw("--cancel-after") asked with no
surviving process; 1 for any other outcome (failed, timed out, lost, a
different hash); 2 if nothing started (unknown case or option, missing input,
another run holding the server); 130 after Ctrl-C, which cancels the job. A
headless run watches and cancels only the job it started itself.

== Five loops are explore-only (#raw("/five-loop/"))

The five-loop page shows the fixture's definition, the Ordered reduction plan
(owners and supports, read from #raw("5L/ibps/selection.json") of a configured
or installed 5L pack; no IBP payload is decoded) and the saved exact
zero-momentum result (#raw("examples/data/five_loop_q0"), checked against its
recorded sha256). It never starts a five-loop reduction by default. An operator
can enable one zero-momentum run by setting #raw("VAKINT_GALLERY_ALLOW_5L=1")
and configuring the compiled #raw("VAKINT_GALLERY_BENCHMARK") executable
(#raw("five_loop_qcd_ladder")), the #raw("VAKINT_GALLERY_ARTIFACTS") checkout,
#raw("VAKINT_GALLERY_WORK") and #raw("VAKINT_GALLERY_PREPARED_INPUT_ROOT")
(a #raw("zero-momentum/") bundle with #raw("prepared-input.json"),
#raw("vakint-input.symbolica") and #raw("projector-conditions.symbolica"); the
native CLI validates it and never falls back to fresh preparation). That run is
a background job like the others: bounded cores, deadline and Cancel. The
same opt-in offers the inspection page's 1L–5L catalogue load. The
quadratic projection has no Run action. The graph page contracts the bubble and
the three- and four-loop ladders (well under a second each); the five-loop
numerator (about 16 s that the page kernel could not interrupt) is shown as
geometry only.

The ten small teaching formulae use their objects' existing native HTML protocol
through the shared gallery presentation method. This retains native MathML and
any bounded-preview omission notes without a view-owned widget model. The
underlying tensor/Expression objects, native gamma/color/dot contraction, and
interactive FeynKit graphs are unchanged. Verify every default teaching formula
is complete, and both directions of every formula-tab switch work, in the actual
browser. Never strip native omission markers to make a preview appear complete.

== Validation

Run #raw("marimo check 0*.py") and the repository Ruff formatter/checker in the
coherent environment. Run the launcher, cancellation and configuration tests with
#raw("python -m unittest discover -s . -p test_shared_gallery.py"); they use fake
commands or short #raw("sleep") processes and perform no graph preparation or
reduction. With the native host, run the 3L and 4L cases headlessly as above
and once with #raw("--cancel-after"); each must exit 0 with #raw("\"ok\": true"). Inspect all five apps in a browser, including native SVG,
formula rendering, lookup forms, source-limit warnings, exact closure controls
and background status. A notebook source passing static checks does not by itself
establish native runtime interoperability. Preserve the actual browser and
producer receipts with the release validation evidence.
