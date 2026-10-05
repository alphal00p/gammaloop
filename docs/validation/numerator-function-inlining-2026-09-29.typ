= Numerator function inlining validation

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

#set text(size: 10pt)

The current compact-expression implementation recovers main-like evaluator
speed when generated numerator functions are inlined during numerical lowering.
This experiment covers the bubble-dressed double triangle in explicit 3D with
ordinary MUV subtraction and the complete sum of 102 orientations.

The experiment used `inline_numerator_functions = false` by default.
The current default is `true`; setting it to `false` retains shared evaluator
bodies. Inlining registers generated component and coefficient functions
with Symbolica's `Always` inlining policy. Tensor contraction, compact energies,
scalar-factor collection and root aliases retain their existing behavior.
Function-map replacements remain disabled, so this introduces no earlier
symbolic body expansion. Dual evaluators continue to inline unconditionally.

The experimental executable was built from commit
`47f3432f2d6034b0b72f6b524e916ad2add26845`, based on the frozen retained-function
implementation `6f192251c9dacd557d82eaac45e4425948934bf7`. Both use Symbolica
`70375b9ef06411d8bac3e14d1f57382384547eef`. The generation settings differ only
in the new flag: direct translation remains enabled, compilation disabled,
Horner iterations one, and the UV prescription and physical inputs unchanged.
Fresh states were generated with the matching executable.

#table(
  columns: (2fr, 1fr, 1fr),
  table.header([Measurement], [Retained functions], [Inlined functions]),
  [Complete generation], [275.403 s], [244.594 s],
  [Peak generation RSS], [1.307 GiB], [1.351 GiB],
  [Symbolica evaluator build], [1.026 s], [5.827 s],
  [Double kernel median], [4.779 ms], [1.303 ms],
  [1000-bit kernel], [480.635 ms], [82.945 ms],
)

Generation timings come from separate runs on a shared host. The numerical
timings use successive paired runs after generation: one discarded Double
evaluation, three measured Double evaluations and one forced 1000-bit
evaluation. Each evaluation uses one identity rotation and the full orientation
sum at spatial coordinates `(17, 23, 31, 37, 41, 43, 47, 53, 59)` in ordered
basis `(6, 5, 8)`. Kernel timings exclude state loading and command setup.
The observed improvements are 3.67 times for Double and 5.79 times for 1000-bit
arithmetic. Historical main timings at this point were 1.277 ms and 87.114 ms.

All six recorded source-size and count metrics agree before numerical
construction: the outer atom and alias bodies occupy 17,269,327 bytes, with a
13-byte root and 214 scalar aliases; 107 recorded function entries contain
14,149,571 body bytes. Generated functions registered for retention decrease
from 59 to zero. Root additions increase from 14,256 to 124,516 and
multiplications from 19,300 to 132,236 as their bodies enter the root program.
The retained root counts omit work inside its callees.
The aggregate function-call counter decreases from 1,280 to 12. That counter
also includes noninteger powers and built-in mathematical functions; it does
not count only retained numerator calls. The exact heads of those 12
instructions were not decoded in this experiment.

The inlined and retained results pass the signed complex comparison. All three
Double outputs agree to relative error `1.70e-14`. Both 1000-bit calculations
export the same imaginary part, `-2.2199969069050517e-22`; their exported real
parts differ only at the subnormal boundary. Native JSON narrows these outputs
to Double, so this checks those exports rather than every high-precision bit.
The focused scalar regression additionally compares the two paths directly in
1000-bit arithmetic within absolute tolerance `1e-290`.

Inlining leaves the current/main value ratio `2.2553906100326313` unchanged.
This rules out retaining numerator calls as the cause of that discrepancy at
the tested point. The main/current normalization difference remains unresolved.

Validation includes the focused retained-body regression, Clippy with only the
existing `uv/profile.rs` argument-count warning, and independent checks of
saved settings, model values, all 102 ordered orientations, signed momentum
routing, native precision and output metadata, and state/executable hashes.
No numerical or generation jobs remain running. The adjacent
The measurements above summarize the original run; its raw receipts and generated states are not distributed.
