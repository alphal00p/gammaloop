= R3 outer SVG scale correction

Frozen patch base: `5ba04864ca506e259d9eef4a54f65acdd8323fc2`.
The sole proposed path is `crates/linnet-py/tests/test_svg_browser.py`.
The new review matches SHA256
`41d01ed9f2fda6b85eefc7d5bea0877e18da7a2a29c196e36cec2dd001349953`.
All prior packets are preserved untouched. No tracked patch or ref was changed.

== Contract and normalization

The authored SVG header in output.rs has matching numerical viewBox sizes
and pt dimensions. The established natural and fitted camera formulas remain
unchanged. The two representation boundaries repaired in R2 remain unchanged:
containment uses binary32 authored extents with the original 1e-6 margin;
tall outer height keeps its separate delta=1/64 CSS pixel equality check.

The missing stage is between the outer viewBox and its CSS viewport.
interactive.js sets both the outer viewBox and CSS height using graphHeight;
the engine snaps CSS layout while its SVGRect viewBox read is binary32.
The default meet transform is uniform, with scale
`min(viewport_width/viewBox_width, viewport_height/viewBox_height)`.
A camera screen CTM contains that outer scale. The fixed-size node-0 target's
screen width and height contain it too. Divide each by the measured outer
scale to express the original camera and target-size assertions in the outer
SVG's viewBox units. Translations do not affect dimensions. No ratio direction,
natural/fitted formula, positive percent rounding, fixture, guard, locator,
delta=1e-4 or places=2 changes. After plus/resize the helper rereads the current
outer transform; the initial value is not reused.

Normalization is conditional on independent assertions, not an unchecked
cancellation. The helper requires a positive, uniform, unrotated outer
matrix (a=d, b=c=0), a positive viewBox height, the exact expected column
viewBox width, a screen width within one layout unit of that column, and
viewBox/CSS heights consistent with the independently computed fitted height.
It checks matrix scale against the actual viewport/viewBox meet relationship.
Finally `abs(a-1)*viewBox_height` must stay within one layout unit plus the
representation term. This last quantity is in CSS pixels, not pt. The tests
therefore cannot absorb an arbitrary CSS zoom, wrong outer transform or
capped preview merely by dividing it away.

The bounded helper is local to the existing test. A bounded frozen-source
search found no existing outer-scale assertion helper in linnet-py tests or
the documented browser fixture harness. No product abstraction was added.

== Binary32 term

For positive normal dimensions the relative rounding unit is u=2^-24.
The source natural formula can be written as
`max(A, min(2.25*A, 450/W32, 300/H32))`, where A is the authored CSS width
read divided by W32. Thus graphHeight is a min/max of positive expressions
with at most three uncertain binary32 factors: H32, authored CSS width and
the reciprocal of W32. The outer viewBox adds one binary32 read; the screen
matrix adds one possible binary32 rounding. The conservative bound is
`gamma5 = 5*u/(1-5*u)`. The dimensioned representation term is
`expected_height*gamma5`. Double arithmetic is much finer and covered by the
unused slack in this bound. This term is source-derived, not calibrated from
measurements. Matrix-versus-meet comparison uses 2^-23, one binary32 ulp near
unity. The allowance for layout remains one 1/64 CSS pixel unit; neither
camera delta nor target places is relaxed.

The tall R2 assertion is deliberately preserved verbatim, including its
existing one-unit allowance. The review's rare combined-representation
residuals there and possible decimal parser rounding remain limitations to
diagnose on the installed host; this correction does not broaden that repair.

== Actual isolated browser evidence

The preferred T3 tools were used: preview_status, preview_open with
reuseExistingTab=false/open=false and the incognito profile, then the owned
blank tab only. No user tab, login or external site was used. Exact frozen
interactive.js and interactive.css bytes were embedded in the retained local
fixture document. A CSP denies external resources; the engine reported zero
resource requests. Source and HTML hashes are in source-bindings.json.

The real browser reports Chromium Chrome/154.0.0.0 on Linux, viewport
1280x800 CSS pixels and devicePixelRatio=2. The high-entropy build-version API
was unavailable on about:blank; the recorded version is the actual reported
UA version, not an inferred patch build. WebKit is not exposed by these
collaborative tab tools and was not exercised. The own tab was closed through
the normal lifecycle tool. No browser was built, installed or downloaded.

Nine initial fixtures were exercised: the four original integer fixtures,
the review's four fractional regime examples, and a nearby fractional
boundary. Their benign static node-0 square is 36.0156 authored units. This
is a fixture input and not a measured golden. Three plus/300px-resize/reset
sequences used real DOM keyboard-event handlers from the exact JavaScript.
The fractional fixtures satisfy the product test's small-reset guard and
reset to 100%. The integer 120x80 fixture does not satisfy that guard
(120*3 > 300), so its independently predicted reset is 83%, which the engine
reported. Its height changes from 240 to 200px after narrowing; this exercises
the actual viewBox/viewport relationship rather than assuming height stays
constant. These fixture semantics do not alter the product guard or checks.

All nine normalized camera scales, percentages, target dimensions and
containment bounds passed. The two tall fixtures passed the preserved R2
height limit. Relational width and height label checks passed at places=2.
The fractional fitted fixture's unnormalized CTM differed from its formula
by more than 1e-4 while its normalized value passed at delta=1e-4, confirming
that this is a real outer-scale boundary in this engine.

Six real DOM controls were rejected: a 1.01 wrong camera scale, a 1.01 outer
transform, CSS zoom 1.01, a 1.02 target-width factor, a 300px tall cap, and a
10x10 camera crop. The outer-transform and cap controls retain the correct
relative camera scale; their rejection comes from the independent outer
bounds. The target-factor control is rejected at the unchanged places=2.
No observed result was copied into an expectation.

The Python validator executes the exact proposed local outer helper on the
recorded engine geometry, and applies the unchanged source-derived camera,
ratio and containment assertions. Other tests and module items have identical
ASTs. The R2 assertions, fixed tolerances and locator are checked for
preservation, and exact unified-patch replay reproduces the supplied after
bytes. Separate Fraction probes demonstrate the dimensional one-unit plus
representation bound; those probes are explicitly mathematical simulations,
not additional browser evidence.

== Supersession and limits

The previous packet's statement that a 1e-5 crop on either axis is rejected
was overbroad. It described selected probes only. An unfitted column has
intentional horizontal slack; even fitted camera widths may have binary32
slack. A tiny width reduction does not universally crop authored content.
That prose is superseded additively here, leaving the original evidence
untouched. The new control removes enough extent to be a meaningful crop.

This earns isolated exact-JavaScript/static-SVG Chromium fixture evidence
only. It is not product rendering, an installed-current-host test, or full
Chromium/WebKit test_svg_browser credit. The product fixture geometry and
all behavior checks still require the exact final installed host with both
engines and Playwright present. An import-based unittest skip is not a pass.
No kernel/product build, licensed execution, Nix, compiler, producer, secret,
license/config read or network was used. Ruff was not on PATH, so formatting
was inspected but no formatter run is claimed.

All writes are within this assigned ignored packet. Frozen source/review
copies, patch, hashes, HTML/input/probe, actual results, separate validation,
preparation drafts and bounded tool-error disclosures are retained for the
integration owner and fresh independent review. Initial transport limits
(truncated internal input, oversized evaluation expression and trailing
whitespace) were corrected before successful execution. No failed browser
assertion was hidden and no skip was counted as success. No background
process remains, and the owned browser tab is closed.
