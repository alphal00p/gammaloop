= Installed host browser correction

The unchanged f47 readable-scale method passes in pinned Chromium 153.0.8010.36
and WebKit 26.5, revision 2311, after correcting the private runtime capsule.
The final run executed one method, four actual native pinned drawings, two
successful engine subtests, two launches and two pages at DPR 1. There were zero
skips, failures or errors. Initial camera/crop/label checks, plus, resize and reset
all ran with the original formulas, fixtures, locators, regime guards and
precision bounds. Root independent review remains pending.

== Cause and correction

Chromium initialized both the camera and DOM label correctly. The actual native
HTML replay had textContent and output value `100%`, but innerText was empty and
the rendered label had a zero-sized box. This persisted through fonts.ready and
an animation frame, without page errors. The empty-root capsule had the DejaVu
font bytes mounted but no fontconfig directory list. An identical-input replay,
changing only FONTCONFIG_FILE to the private configuration enumerating that
existing verified font, rendered `100%` with positive text dimensions; camera
and outer geometry stayed identical. No readiness wait or oracle change is
needed. This matches the explicit fontconfig input in the repository browser
recipe. The earlier R4 fixture method and JS/CSS bytes match, but its static
producer and inherited filesystem under unshare differ from the installed-host
empty-root capsule. No old fixture proof was rerun.

WebKit independently created its browser context, then its first page target
crashed with “Could not create WPE EGL display: EGL_SUCCESS. Aborting...”. The
wrapper stayed in page creation until the diagnostic watchdog interrupted it.
The pinned libglvnd EGL loader searches absent system vendor directories. A
source-identified, existing Mesa 26.2.1 vendor and exact software-driver/native
library files, plus its small drirc data directory, were supplied read-only.
Explicit EGL vendor and driver discovery paths then allowed the same pinned
WebKit to create a page and finish loading about:blank at DPR 1. The successful
blank diagnostic's automation and Mesa probe warnings are retained. No dbus,
display server, host driver tree, GPU, host root or credentials were mounted.
The extra runtime is fully inventoried in `egl-runtime.json`: 71 native/data
files or aliases and 72 additional file/data mounts, including the drirc
directory. The original 479 store mounts and all browser/Python pins remain.

== Evidence and source identity

The immutable source is f47f27213a1c55291239c9bc0f334a7ea1cd1cc1 and the
kernel is 8eb7c3729d6abefc916decc1bdf0c06536f8af80. All 4602 source entries,
479 store NARs, installed venv files, selected producer outputs, native library
resolutions, and the supplied host/WebKit report identities were verified before
use. The exact unchanged test SHA256 is
`ee700e0629fad62358e1f0a778353757a9a305bdbc131f2ae2675b7b74d82c6c`.
Native HTML and SVG captures, exact embedded JS/CSS bindings, console/page-error
records, timing, DOM states and geometry are retained. Both engines read initial
labels `100%`, `100%`, `100%`, `43%`, then `105%` and reset `100%`.

`report.json` binds the originals, hashes, real counts, per-engine outcomes and
all seven capsule attempts. The original host failures remain untouched.
`setup-failures.json` retains the pure command and ELF-reader setup failures.
The diagnostic four-drawing capture and final four-drawing native method total
eight native drawings in this phase. Static replay is diagnostic evidence and
is not substituted for the successful full native method.

`capsule.py` and `helpers/selected_browser_runner.py` are the executed private
controller and observation helper. `capsule-correction.patch` shows the minimal
review-only change to the historical writable-state owner: font configuration
and explicit read-only EGL inputs. Its independent font mountpoint avoids
creating a file inside the historical read-only helper mount. Port the proposal
into the final owner's fresh capsule; never execute it in sealed H. No tracked
product/test patch or fresh host build is required for this correction.

== Containment, preservation and limits

A separate unlicensed containment audit passed: fresh network namespace,
loopback only, no routes, read-only source/helpers/venv/libraries, no Nix socket,
credentials, VCS or compiler executables. Only fresh private work was writable;
evidence remained outside the child mounts. The three private license variables
were absent from pure children and stripped from the Playwright driver/browser
environments. The license entered only the restricted immediate native-product
exec through an anonymous pipe; no value, hash, argv entry, file or receipt was
created. No network, Nix, Cargo, compiler or producer-build command ran.

The baseline WebKit watchdog interruption is not a test pass or skip. Diagnostic
watchdogs do not alter test timeout/retry/memory policy; 20 GiB was not adopted.
No product physics/model, expectations or tolerances changed. The original
passing Rust tests, UFO crossing, squared display and signature checks were
preserved without reruns. API/stubs and the later token owner remain separate.
Main, sealed packets, tracked source/index/refs and JJ were not changed.

This is bounded installed-host camera validation using the existing fixed
external metrics snapshot. It does not qualify browser font measurement, a full
widget/site, future rebuilt host, combined WP3 work or the recorded synthetic
R4-N3 residual. The additional Mesa system runtime is explicit review input;
the eventual combined runtime must supply and verify it or independently
qualify an equivalent. No independent approval is claimed.

All owned capsule processes were reaped. Disposable browser state and caches
were removed. The new ignored correction directory retains only review/debugging
sources and evidence; the reused installed venv, store inputs, original packets
and workspace remain available to their owners.
