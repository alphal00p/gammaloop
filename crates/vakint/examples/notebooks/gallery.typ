= Interactive Vakint field notes

These four Marimo applications share the native HEPKit, tensor, graph and
Symbolica types from one coherent host. They demonstrate graph dressing and
efficient gamma/color contractions, precision-aware constant/master inspection,
the zero-momentum four-rung vacuum fixture, and exact numerator/pinch closure.
The quadratic Taylor fixture is retained as deferred history. Its reduction,
retesting and cache comparison require a later user decision and have no gallery
Run action.
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
does not change the tailnet or existing services itself. The four routes are
#raw("/graphs/"), #raw("/inspect/"), #raw("/five-loop/") and #raw("/closure/").
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
load button validates #raw("VAKINT_GALLERY_ARTIFACTS") when set, otherwise calls
the official latest loader in the caller's working directory. A loaded engine
is reused for that session's subsequent lookups. Every master remains a native
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
An absent receipt is visibly pending and makes no timing
or successful-reduction claim. The producer remains responsible for scientific
validation and source/data provenance; a local hash does not certify those facts.

The optional #raw("zero-momentum-running-cache") receipt is an immutable
observation of an optimized producer and its cumulative native cache. Display
enqueued and acknowledged-durable record counts separately. Neither count
establishes completion of the full physical reduction. A fresh process reusing
that same input is a seen-input check, distinct from the held-out suite.
Keep final reduction times, speedups and estimates absent until their completed
authenticated receipts are available. Refreshing this evidence does not require
changing the deployed native host or starting a numerical calculation.
The separate #raw("zero-momentum-light-cache") and
#raw("zero-momentum-heavy-cache") receipts contain completed comparisons for
their named tier. Keep them pending until exact results, cache identities and
controlled measurements have been validated.

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

== Explicit expensive runs

To enable the saved-input reduction button, configure the qualified compiled
#raw("VAKINT_GALLERY_BENCHMARK") executable, #raw("VAKINT_GALLERY_ARTIFACTS")
checkout, #raw("VAKINT_GALLERY_WORK") development directory, and
#raw("VAKINT_GALLERY_PREPARED_INPUT_ROOT"). This last directory contains
#raw("zero-momentum/") bundle with
#raw("prepared-input.json"), #raw("vakint-input.symbolica") and its complete
#raw("projector-conditions.symbolica"). The selected bundle must exist; the native
CLI validates its exact payload, projection, scientific identity and producer
provenance. Missing or invalid configured input reports an error and never falls
back to fresh preparation. Historical quadratic bundles can be retained,
but the gallery launcher rejects that projection while it is deferred.
Bind the executable's bytes and capabilities to the reviewed build and
consumer-validation receipts before enabling the button.
These paths are deployment configuration, not browser-controlled commands.

A separate worker runs the same Rust example with #raw("--projection"),
#raw("--prepared-input"), #raw("--artifacts"), #raw("--work-dir") and
#raw("--progress-snapshot"). The reduction uses RustRed's ordered scheduler,
which reads and writes no persistent cache. Do not reuse a failed job's
directory for another attempt.

Each fresh process uses the same complete scalar input without
repeating native graph generation or trace contraction. The report separates
current verification/loading from original
producer preparation; historical preparation is excluded from current end-to-end
time. Validate fresh scientific generation independently with
#raw("--fresh-prepare --prepare-only") in a new directory before installing a
bundle. Reusing scalar input does not cache its IBP result.

The browser stays responsive and reads atomic status records. A process lock
permits one gallery benchmark at a time. Physical reductions have no default
time cutoff.
The worker refreshes elapsed time and worker/machine RAM observations every five
seconds while waiting for the collector. Its native #raw("WORK/monitor.json")
snapshot is displayed through a bounded JSON preview with its actual file age.
Before the first Native event, progress remains explicitly pending; repeated
heartbeats do not imply scientific progress or a new Native observation.
Requests, logs and resulting receipts remain below the configured workspace.
The depth-first cache evidence shown by the closure page is historical: Vakint
no longer prepares, captures or reads depth-first caches. The graph tutorial's
explicit tensor-contraction button is a separate, interactive algebra
demonstration, not a controlled timing sample of the physical zero-momentum fixture.

The ten small teaching formulae use their objects' existing native HTML protocol
through the shared gallery presentation method. This retains native MathML and
any bounded-preview omission notes without a view-owned widget model. The
underlying tensor/Expression objects, native gamma/color/dot contraction, and
interactive FeynKit graphs are unchanged. Verify every default teaching formula
is complete, and both directions of every formula-tab switch work, in the actual
browser. Never strip native omission markers to make a preview appear complete.

== Validation

Run #raw("marimo check 0*.py") and the repository Ruff formatter/checker in the
coherent environment. Run the cheap launcher/configuration tests with
#raw("python -m unittest discover -s . -p test_shared_gallery.py"); they use fake
commands and perform no graph preparation or reduction. Inspect all four apps in a browser, including native SVG,
formula rendering, lookup forms, source-limit warnings, exact closure controls
and background status. A notebook source passing static checks does not by itself
establish native runtime interoperability. Preserve the actual browser and
producer receipts with the release validation evidence.
