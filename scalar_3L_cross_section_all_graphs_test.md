# Scalar 3L Cross-Section All-Graphs Inspect Coverage

## Current Scope

- Test file: `tests/tests/test_runs/scalar_3L_cross_section_inspects.rs`.
- The current process generation produces graph labels `GL00` through `GL48`.
  This is 49 labels in the current branch, so the harness covers all present
  labels even though the original planning note referred to 45 graphs.
- The exhaustive checks live under the `slow` test module and are excluded from
  the base `just test_gammaloop` profile. They are visible through
  `just test_gammaloop slow`.
- Six timed rich-inspect anchors have been moved out of the `slow` module and
  are now part of the base `just test_gammaloop` profile:
  `GL02` baseline, `GL00` with `Q(7)^2`, `GL16` with `Q(1)^2`, and quartic
  `Q(1)^4` probes for `GL02`, `GL16`, and `GL24`.
- `CONSIDER_COMPARISON_WITH_LTD = true`. The retained slow sweep compares:
  - `CFF + local UV from 3D expansions`,
  - `CFF + local UV from expanded 4D integrands`,
  - `LTD + local UV from expanded 4D integrands`.
- `ALLOW_DISABLING_UV_SUBTRACTION = true` and
  `ALLOW_DISABLING_THRESHOLD_SUBTRACTION = true` keep known broader limitations
  expressible per graph. Turning either flag to `false` promotes the
  corresponding fallback to a hard test failure.

## Common Configuration

- Process generation command:
  `generate xs scalar_1 > {scalar_0 scalar_0, scalar_0 scalar_0 scalar_1} | scalar_0 scalar_1 [{{3}}] --allowed-vertex-interactions V_3_SCALAR_001 V_3_SCALAR_000`.
- Default higher-energy numerator probe:
  `(Q(1,mink(4,1))*Q(2,mink(4,1)))*(Q(5,mink(4,2))*Q(6,mink(4,2)))`.
  This is the lowered default retained for the all-graph slow sweep; heavier
  four-factor probes should stay in targeted tests where their runtime is
  acceptable.
- Each graph test evaluates both:
  - no extra numerator,
  - the configured higher-energy numerator probe.
- Each retained case performs a rich local-inspect comparison of the two CFF UV
  construction modes, including total weight, event grouping, event weights, and
  additional event weights.
- The CFF local-3D result is the benchmark source for the total-weight snapshot
  once the rich CFF local-3D/local-4D comparison has passed.
- In explicit-orientation-sum CFF, the nominal local-3D setting is allowed to
  route a UV forest term through the expanded-4D bridge when the UV-rescaled
  subgraph contains duplicate leading loop-denominator signatures. This is a
  structural local-UV boundary, not a graph-name fallback: after the UV
  expansion those denominators become confluent `mUV` denominators, and the
  representation-independent 4D bridge is the canonical way to build that local
  counterterm.

## Current Per-Graph Configuration

| status | graphs | notes |
| --- | --- | --- |
| three-way parity passes with UV and thresholds enabled | `GL02`, `GL05`, `GL06`, `GL09`, `GL13`, `GL16`, `GL20`, `GL22`, `GL25`, `GL27`, `GL28`, `GL29`, `GL31`, `GL32`, `GL34`, `GL39`, `GL40`, `GL41`, `GL42`, `GL43`, `GL47`, `GL48` | `GL22`, `GL25`, `GL28`, `GL32`, `GL39`, `GL43`, and `GL47` use the lighter single-dot probe `Q(1,mink(4,1))*Q(2,mink(4,1))`; `GL48` uses final states `scalar_0 scalar_0` and numerator `Q(7,mink(4,1))*Q(8,mink(4,1))`. |
| three-way parity passes with UV disabled | `GL00`, `GL07`, `GL08`, `GL10`, `GL11`, `GL12`, `GL14`, `GL15`, `GL17`, `GL18`, `GL19`, `GL21`, `GL23`, `GL26`, `GL33`, `GL36`, `GL37`, `GL38`, `GL44`, `GL45`, `GL46` | `GL23` uses the lighter single-dot probe `Q(1,mink(4,1))*Q(2,mink(4,1))`; the others use the default two-dot numerator probe and keep thresholds enabled. Disabling UV is guarded by `ALLOW_DISABLING_UV_SUBTRACTION`. |
| three-way parity passes with UV and thresholds disabled | `GL01`, `GL03`, `GL04`, `GL30` | `GL03` uses `Q(7,mink(4,1))*Q(8,mink(4,1))`; `GL30` uses the lighter `Q(1,mink(4,1))*Q(2,mink(4,1))`. `GL04` is kept in this fallback because the LTD threshold-enabled all-graph case is above the intended slow-test runtime budget. `GL30` hit a threshold singular sample and a UV local-3D/local-4D mismatch with UV enabled. |
| three-way parity passes with UV enabled and thresholds disabled | `GL24`, `GL35` | `GL24` uses the lighter `Q(1,mink(4,1))*Q(2,mink(4,1))` numerator; its UV forest triggers the duplicate-leading-denominator expanded-4D bridge. `GL35` produced a matching NaN threshold counterterm at the fixed inspect point with thresholds enabled, so the exhaustive benchmark disables thresholds for this graph. |
| retained scalar rich-inspect parity currently fails | none | The current retained scalar cross-section suite has no known rich local-inspect parity failures among the three representations. |

## Current Validation State

- The harness is implemented for all current labels with LTD included for every
  retained graph/numerator configuration.
- The selected full GammaLoop suite
  `just test_gammaloop` passed 1090/1090 tests on May 6, 2026 after the
  zero-loop energy-bound probe was short-circuited for empty bounds, the
  imported pure-tree external-tree inspect coverage was restored, and the six
  fast scalar cross-section rich-inspect anchors were promoted into the default
  suite.
- Latest full pure-CFF scalar sweep:
  `env RUST_MIN_STACK=33554432 INSTA_UPDATE=always cargo nextest run --cargo-profile dev-optim -P test_gammaloop -p gammaloop-integration-tests --run-ignored all --ignore-default-filter -E 'test(/scalar_3l_cross_section_/)' --no-fail-fast`
  passed on May 6, 2026. This includes the all-graph baseline, the retained
  quadratic numerator probes, and the retained quartic subset.
- Latest LTD-enabled scalar validation:
  `env RUST_MIN_STACK=33554432 cargo nextest run --cargo-profile dev-optim -P test_gammaloop -p gammaloop-integration-tests --run-ignored all --ignore-default-filter -j4 -E 'test(/scalar_3l_cross_section_/)' --no-fail-fast`
  passed 149/149 retained scalar tests in 1586.599s on May 6, 2026, after the
  fast-anchor promotion. This
  covers the all-graph baseline, the retained quadratic probes, and the
  retained quartic subset with the rich three-way inspect comparison enabled.
- The May 6 LTD-enabled rerun refreshed three benchmark snapshots
  (`GL26`, `GL33`, and `GL44`, all in the `quadratic_pair_numerator` case with
  UV disabled and thresholds enabled) after fixing the repeated-sector
  lower-graph construction and finite-pole contact completion. These snapshot
  updates were accepted only after the rich three-way local-inspect comparison
  had already succeeded.
- The expanded-4D UV bridge now shape-expands a zero local residue vector to the
  finite integrated-UV residue shape. This removed one false residue-count
  boundary and exposed the genuine `UV forest union nodes are not supported`
  boundary on `GL10`, which remains a guarded UV fallback.
- Existing snapshots cover every currently retained graph/numerator pair in the
  scalar sweep. Benchmark snapshots were only accepted after the rich
  local-inspect comparison between CFF local-3D, CFF local-4D, and LTD local-4D
  succeeded.

## Quadratic Energy Numerator Probe

- A second slow section,
  `slow::quadratic_energy_numerators`, reuses the same graph list and the same
  per-graph UV/threshold configuration as the all-graph pure-CFF sweep.
- Most graphs are split into two independent tests:
  - `q1_squared`: numerator `Q(1,mink(4,1))*Q(1,mink(4,1))`;
  - `q7_squared`: numerator `Q(7,mink(4,1))*Q(7,mink(4,1))`.
- The split is deliberate: it keeps coverage for one squared-edge probe from
  being masked when the other probe hits a generation boundary.
- The q7-squared probe is omitted for `GL02`, `GL05`, `GL27`, `GL29`, `GL31`,
  `GL34`, `GL35`, and `GL42`: these graph/edge combinations enter expensive
  generalized repeated-sector paths and remained CPU-bound past the intended
  slow-test budget. The q1-squared probe is still retained for these graphs.
- The retained `q7_squared` cases pass with the three-way comparison enabled;
  the removed q7 cases are runtime exclusions, not known parity mismatches.
- Latest `q1_squared` sweep:
  `env RUST_MIN_STACK=33554432 INSTA_UPDATE=always cargo nextest run --cargo-profile dev-optim -P test_gammaloop -p gammaloop-integration-tests --run-ignored all --ignore-default-filter -E 'test(/quadratic_energy_numerators::scalar_3l_cross_section_gl[0-9][0-9]_quadratic_energy_inspects_match::q1_squared/)' --no-fail-fast`
  passed 49/49 tests on May 6, 2026.
- The former `q1_squared` lower-sector failures on `GL16`, `GL20`, `GL28`,
  `GL32`, and `GL43` were caused by reconstructing auxiliary lower-sector
  graphs from momentum signatures. The CFF lower-sector builder now projects
  actual subgraphs from the parsed topology, so graph reconstruction from
  signatures remains a CLI/testing convenience and is not used in generation
  logic.
- Snapshot coverage currently contains 49 `q1_squared` benchmarks and retained
  `q7_squared` benchmarks for every graph except the runtime exclusions listed
  above. No `.snap.new` files are left pending.

## Quartic Energy Numerator Probe

- A targeted slow section, `slow::quartic_energy_numerators`, checks the
  numerator
  `Q(1,mink(4,1))*Q(1,mink(4,1))*Q(1,mink(4,2))*Q(1,mink(4,2))`.
- The retained subset is `GL00`, `GL02`, `GL06`, `GL10`, `GL16`, `GL20`,
  `GL24`, `GL30`, `GL41`, and `GL42`. This includes UV-on, UV-off,
  threshold-on, threshold-off, and the lower-sector cases that exposed the
  signature-reconstruction boundary.
- Latest targeted quartic sweep:
  `env RUST_MIN_STACK=33554432 INSTA_UPDATE=always cargo nextest run --cargo-profile dev-optim -P test_gammaloop -p gammaloop-integration-tests --run-ignored all --ignore-default-filter -E 'test(/quartic_energy_numerators::scalar_3l_cross_section_gl(00|02|06|10|16|20|24|30|41|42)_quartic_energy_inspects_match::q1_quartic/)' --no-fail-fast`
  passed 10/10 tests with LTD enabled on May 6, 2026.
- A trial all-graph quartic sweep was not retained: the heaviest tail had
  several individual graph tests still running after more than 11 minutes.
  The all-graph CFF anchor remains the baseline plus the two quadratic numerator
  sweeps above.

## Default-Suite Fast Anchors

- The default suite now includes a curated sub-20s-rich-inspect subset while the
  exhaustive scalar sweep remains under `slow`.
- Timed anchors:
  - `GL02` baseline, 7.7s wall, UV and thresholds enabled, no numerator plus
    the default higher-energy probe.
  - `GL00` with `Q(7,mink(4,1))^2`, 2.8s wall, UV disabled and thresholds
    enabled, covering a repeated-sector scalar graph.
  - `GL16` with `Q(1,mink(4,1))^2`, 4.3s wall, UV and thresholds enabled,
    covering a lower-sector energy-bound case.
  - Quartic `Q(1)^4` probes for `GL02`, `GL16`, and `GL24`, 8--10s wall each;
    `GL24` keeps thresholds disabled and exercises the duplicate-leading-
    denominator expanded-4D UV bridge.
- Slower but important anchors such as `GL06`, `GL09`, and `GL48` remain in the
  slow module. They continue to be covered by the 149-test scalar sweep.
