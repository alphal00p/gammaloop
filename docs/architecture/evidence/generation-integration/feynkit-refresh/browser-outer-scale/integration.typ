= Outer SVG screen-scale correction: R4 integration record

This is additive evidence for the existing community-refresh investigation,
not source approval or a new runtime qualification. The code follow-up is
`f47f27213a1c55291239c9bc0f334a7ea1cd1cc1`, on preserved shared metadata
`c3449988839d9ad3389dd9fd4c4b8bb9ec447055`. Earlier code `5ba04864` and
evidence `2c7dbbb0` remain unchanged.

The sole source delta applies the root-selected Option A for R3-1. Camera
screen matrices and stable node-label dimensions are divided by the outer SVG
screen scale. Positive uniform unrotated geometry, viewport and viewBox sizes,
default meet and bounded departure from unity are checked independently before
normalization. The source-derived gamma5 representation term and one 1/64 CSS
pixel layout allowance do not relax camera delta=1e-4 or label places=2. R2
containment and tall-height repairs remain. The worker's complete derivation is
retained verbatim in `worker-rationale.typ`.

The root separately authorized one comment correction: “one layout unit is a
CSS pixel” becomes “one layout unit is 1/64 of a CSS pixel”, wrapped on two
lines. This changes no AST, code or allowance. The final browser SHA256 is
`ee700e0629fad62358e1f0a778353757a9a305bdbc131f2ae2675b7b74d82c6c`;
it intentionally differs from the worker's after hash. The worker packet is
unchanged and fully hash-bound, including its actual fixture evidence.

The worker exercised nine initial static-SVG/current-JavaScript Chromium
fixtures, three plus/resize/reset sequences and six rejected real DOM controls.
This is real isolated fixture evidence. It supplies no installed-current-host,
WebKit or product-test pass. Fraction simulations are separate mathematical
evidence. The older either-axis tiny-crop claim is qualified additively by the
worker's explanation of horizontal slack, without altering the historical note.

Integration checks parse the final Python source, compare its AST with the
worker, verify the single source-path diff and whitespace, preserve all other
4,600 previously tracked files, and verify all 118 registry scope bindings.
No scope or registry bytes change. Explicit immutable pure evaluation, with
both license variables absent, selects the identical generator archive drv
`/nix/store/5pxm7dg3gnz25a30iqzmyvh7kwcrwbsq-gammaloop-nextest-binaries-feynkit-feynkit-generator-0.1.0.drv`.
Its package, test and all 47 dependency-filtered source inputs remain identical;
the browser path is absent from each. The new full source is not called newly
runtime-built. No archive, helper, host, browser or site build/capture was run.

The receipt binds the R3 review, worker packet, final patch and source,
producer-invariance record and pending obligations. Existing c344 shared
qualification, runtime12 and sealed prior-eight records remain immutable.
NixCI remains false. Fresh root-assigned R4 source review, installed-host
Chromium/WebKit coverage and independent combined review remain pending.
