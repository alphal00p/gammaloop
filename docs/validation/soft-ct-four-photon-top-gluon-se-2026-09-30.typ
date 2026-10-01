= Four-photon top loop with a gluon or ghost self-energy

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

#set text(size: 10pt)

Diagnostic date: September 30, 2026. The previous overlap-cancellation audit was
put on hold before creating jj change `rqokynns` in a separate workspace.
The tested source is `032969abd45c53e0a5bb10d883d7b1d3a4bb149b`;
Symbolica is pinned to `578dfcb55fb0662d871456ea8ecc57c1390caa77`.
The final change adds diagnostic inputs and results only. It does not change
production Rust code, dependencies, evaluator settings, or the original graph files.

== What was tested

The physical three-loop graphs are GL262 (gluon bubble) and GL256 (its ghost
counterpart), with four external photons and a massive top loop. The original
projector, momentum-dependent numerator, vertices, edge directions and overall
factors are retained. Their graph factors are respectively −2 and −4; an extra
ghost sign is not added during evaluation.

The Full-H prescription uses IR subtraction for the quadratic gluon/ghost bubble
and the linearly divergent top self-energy, with MUV for logarithmic regions.
Integrated counterterms and threshold subtraction are disabled. All physical
orientation contributions are summed. These are individual graph and local
counterterm tests, not an integrated four-photon amplitude calculation.

The evaluator settings are the existing inlined configuration: numerator-function
inlining enabled, retained scalar aliases, single parametric evaluation, one
Horner iteration, and native compilation disabled. The trace field `compile=true`
is a logging tag, not the value of the compilation setting. Symbolica must still
build its arithmetic evaluator. Generation uses one worker and eight Rayon
threads. The run cap is 900 seconds; the successful 4D retries raise the initial
12 GiB RSS allowance to 100 GiB. Concurrent cases ran on a shared host, so the
wall times are diagnostic measurements rather than isolated benchmarks.

Physical inputs are MT = 173, zero top width, two 86.5 GeV incoming photons,
helicities ++−−, mUV = 9.188, and muR = 91.188. The cards preserve all other
model parameters and external momenta explicitly.

== Routing boundary

The original basis `[4, 9, 12]` fails projected 4D before evaluator generation:
the nested top-self-energy component `12nI` needs child expansion carriers whose
canonical routes contain external momentum. The implementation rejects their
promotion as an affine shift. This happens after 1.118 seconds for GL262 and
0.567 seconds for GL256; it is not a memory cap or a failed numerical scaling test.

Using the basis `[5, 9, 6]` makes both nested boundaries explicit loop coordinates:
`q = Q5 = -Q7`, `b = Q9`, `h = Q6 = Q8`. This valid input rerouting has absolute
Jacobian one and passes the guard. It does not establish pointwise equality of
counterterms regenerated in different bases. Any route comparison here must use
the same basis. The adjacent `routing-notes.typ` gives the transformation,
source references and ghost momentum mapping.

== Generation and numerical results

#table(
  columns: (1fr, auto, auto, auto),
  table.header([Graph], [4D generation], [Peak RSS], [Worst Double/Arb difference]),
  [GL262 — gluon], [606.07 s], [27.53 GiB], [2.20 × 10^-14],
  [GL256 — ghost], [308.96 s], [15.49 GiB], [1.06 × 10^-14],
)

Both generation settings and saved-state/binary hashes were verified. The six
point evaluations plus process startup took 4.872 s for GL262 and 2.269 s for
GL256; these are whole-process times, not isolated evaluator kernel timings.

The native central soft fits give powers 2.00038996 (GL262, R-squared
0.99999968275) and 2.00000751 (GL256, R-squared 0.999999999883). Those scans took
36.009 s and 13.334 s respectively. Their rays have distinct fingerprints.

A separate common canonical ray uses `q = lambda*(30,30,114.5)`,
`b_gluon = (37,41,43)`, `b_ghost = -b_gluon`, and `h = (47,53,145.5)`.
The signed sum, with the original graph factors already included, gives a
measure-weighted slope of 1.99991827 and R-squared 0.999999999415 over the last
nine points (lambda from 10^-6 to 10^-8). The full 17-point fit is 1.99699231,
with R-squared 0.999997675. The minimum ratio of the sum magnitude to the sum
of individual magnitudes is 0.99999208; this ray does not rely on a delicate
additional gluon/ghost cancellation. Raw complex values are retained.

The original-basis explicit 3D attempt reached its 900.28 s cap at 1.899 GiB.
The aligned localized 3D attempt reached 900.34 s at 1.778 GiB. Both remained
in symbolic Taylor construction, before evaluator building; no comparison of
all three routes is established. The failed 12 GiB projected-4D attempts are
also retained separately from the successful larger-memory retries.

Both graphs pass one native representative for every divergent cycle union:
four rays and 68 samples per graph in the aligned basis. GL262 takes 71.576 s;
GL256 takes 25.806 s. The measure-weighted UV powers are:

#table(
  columns: (1fr, auto, auto),
  table.header([Region], [GL262], [GL256]),
  [Bubble, edges 9 and 10], [−1], [−1],
  [Top self-energy, edges 4, 5, 7, 9, 10], [−1], [−1],
  [Disconnected top hexagon and bubble], [−2], [−2],
  [Overall graph], [−1], [−1],
)

All eight fits use all 17 points, with R-squared 1 at exported precision;
the fitted powers differ from the displayed integers by less than 1.4 × 10^-9.
This is the native `only-divergent` selection. It covers each divergent cycle
union once, not every possible basis/subset ray.

The exhaustive `all` attempts hit the 900-second cap for both graphs, at
900.620 s for GL262 and 900.513 s for GL256. The native profiler only writes
its numerical report after completing the sweep, so those attempts yield no
partial numerical certificate. The 280-ray coverage is not validated and the
timeouts are not classified as failed scaling laws.

No individual-orientation scaling, complete integrated amplitude, integrated
counterterm validation, or algebraic cancellation proof is claimed here.


The three generic points are evaluated in Double and forced Arb precision, with
the complete orientation sum. Their relative comparison is a complex-norm
comparison; a small real component consistent with roundoff is not assessed by
its own tiny denominator. Native inspect exports its final values as f64. These
checks are numerical evidence and do not establish algebraic identities.

The native central soft check is `S(e5)` with 17 points from lambda = 10^-4 to
10^-8, seed 1337, and Arb evaluation. Both central gluon lines soften together,
so there is one independent soft momentum. The native profiler fits adjacent
differences of the norm and adds three powers for the soft momentum measure.
Its native pass condition is a positive resulting exponent; the measured exponent
and fit quality are reported separately. The native gluon and ghost profiles
use graph-specific rays and are not silently added together.

The all-limit UV scan requests all 40 loop-momentum bases and all seven nonempty
loop subsets, with 17 scales from 10^8 to 10^12: 280 rays and 4,760 evaluations
per graph. A full-pass claim requires all of that coverage, native slope at most
−0.9, and R-squared at least 0.99. A cap does not count as a physical failure.

== Reproduction and artifacts

Run the commands below from the repository root with a valid Symbolica
environment and a CLI built from the recorded source. Replace the binary path
with the matching build. The recorded runs additionally use the RSS supervisor
and immutable binary/state hash checks retained in the campaign directory.

```sh
gl_bin=/path/to/gammaloop
data=docs/validation/soft-ct-four-photon-top-gluon-se-2026-09-30
state=/tmp/four-photon-top-gl262-check

"$gl_bin" "$data/gl262-aligned.toml" -s "$state" run -c \
  "set global file $data/h-projected-4d.json; run generate; save state -o; quit -n"

"$gl_bin" -s "$state" --read-only-state --no-save-state run -c \
  "inspect -p four_photon_top_self_energy -i GL262 --momentum-space --graph-id 0 --point 30 30 114.5 37 41 43 47 53 145.5 --use_arb_prec; quit -n"

"$gl_bin" -s "$state" --read-only-state --no-save-state run -c \
  "profile bulk -p four_photon_top_self_energy -i GL262 --select 'GL262 S(e5)' --min-scaling=-4 --max-scaling=-8 --n-points 17 --seed 1337 --output /tmp/gl262-soft.json; quit -n"
```

The UV check uses the same saved state:

```sh
"$gl_bin" -s "$state" --read-only-state --no-save-state run -c \
  "profile ultra-violet -p four_photon_top_self_energy -i GL262 --graph GL262 --selected-limits only-divergent --min-scaling 8 --max-scaling 12 --n-points 17 --seed 1337 --use_f128 --output /tmp/gl262-uv; quit -n"
```

Changing the selection to `all` reproduces the more expensive exhaustive scan.
Use the GL256 card and integrand name for the ghost contribution. The explicit
and localized 3D overlays are also retained. Use a fresh state directory for
every generation and keep the exact input/binary identity with each result.

The adjacent directory retains portable run cards, aligned DOT graphs, and
routing and numerical-contract notes. These are sufficient to repeat the
documented generation and profiling commands. The historical receipts,
raw numerical outputs, and evidence archive are not distributed.
