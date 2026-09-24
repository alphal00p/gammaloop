= Build5 finite Monte Carlo execution
<build5-finite-monte-carlo-execution>
Prepared only; no new numerical acceptance is claimed. Every command uses frozen b53 `gate-physics-build5` SHA `2fd035e87812264837a6cb0feb752a6abb63f827a90d1b461285442fd0dd311d`. Its `physics-pilot` stage requires the exact successful build5 preflight, source/state identity, and per-card correctness settings. The preflight is still running at preparation time. No old partial acceptance is substituted.

The companion `screen-build5-commands.json` contains exact argv and resource/log paths. Original cards, physical settings, N and seeds are unchanged. GNU time is available at `/run/current-system/sw/bin/time`. Peak RSS covers the whole driver process, including load/warmup; use `pilots.json[0].elapsed_seconds` for the existing Integrate interval and `2048 / elapsed_seconds` for throughput. The whole-process wall time is separate. Root observes resident memory during warmup and stays below the authorized 300 GB; no unmeasured RAM-fit assumption is made.

```bash
source /tmp/gammaloop-rebase-dev-env.sh
gl_gate_packet=/tmp/gl638-hosted-joint-gate
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-screen-throughput20-build5.resources.txt" \
  "$gl_gate_packet/drivers/gate-physics-build5" \
  "$gl_gate_packet/manifest-screen-throughput20-build5.json" \
  "$gl_gate_packet/results-screen-throughput20-build5" all physics-pilot 20 1 1 1337 \
  > "$gl_gate_packet/results-screen-throughput20-build5.log" 2>&1
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-screen-throughput30-build5.resources.txt" \
  "$gl_gate_packet/drivers/gate-physics-build5" \
  "$gl_gate_packet/manifest-screen-throughput30-build5.json" \
  "$gl_gate_packet/results-screen-throughput30-build5" all physics-pilot 30 1 1 1337 \
  > "$gl_gate_packet/results-screen-throughput30-build5.log" 2>&1
```

Choose workers before examining advanced screen outcomes: use 30 only if the valid, equal-N throughput ratio (30 workers / 20 workers) is at least 1.35, which is 90% of ideal 1.5 scaling, and peak resident memory is at most 300 GB (300,000,000,000 bytes); otherwise use 20. This criterion is fixed before results. The throughput interval includes Integrate warmup; whole-process wall is reported separately. Retain both control results, cards, source hashes and RSS records. Record the choice and chosen baseline output directory. The same seed pairs RNG schedules, not a bit-equality claim for parallel floating-point accumulation.

```bash
# Set to exactly 20 or 30 from the two completed throughput controls.
: "${gl_screen_workers:?Set to 20 or 30 from the recorded throughput-only decision}"
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-screen-advanced-firstseed-build5.resources.txt" \
  "$gl_gate_packet/drivers/gate-physics-build5" \
  "$gl_gate_packet/manifest-screen-advanced-firstseed-build5.json" \
  "$gl_gate_packet/results-screen-advanced-firstseed-build5" all physics-pilot "$gl_screen_workers" 1 1 1337 \
  > "$gl_gate_packet/results-screen-advanced-firstseed-build5.log" 2>&1
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-screen-all-remaining-seeds-build5.resources.txt" \
  "$gl_gate_packet/drivers/gate-physics-build5" \
  "$gl_gate_packet/manifest-screen-all-remaining-seeds-build5.json" \
  "$gl_gate_packet/results-screen-all-remaining-seeds-build5" all physics-pilot "$gl_screen_workers" 1 1 1337 \
  > "$gl_gate_packet/results-screen-all-remaining-seeds-build5.log" 2>&1
```

The selected screen is1 chosen baseline + 6 advanced at 11113 + 14 rows at 22229/33331 = 21 runs × 2048 = 43,008 draws. Four cold loads suffice. The other 2,048-draw baseline is an excluded throughput control, never a fourth seed or duplicated screen row. `manifest-screen-full-build5.json` is the complete analysis declaration, not an additional executable pilot request. The existing analyzer now reads these three real directories directly, preserves per-batch acceptance and enforces unique(mode,seed), equal N and chosen workers. Its formulas remain unchanged; no synthetic pilot collection or state is created. The losing throughput control is retained explicitly and excluded from pooling.

For screen maximum replay, call the frozen max29b `drivers/max-replay-build5` separately on each chosen real pilot directory using that batch\'s matching manifest and a new output directory. This keeps complete original Samples/checkpoints and their source/card lineage. Do not replay the excluded worker control as screen evidence. Per-batch top 2 controls are labelled as such; the later primary confirmation is one complete Cartesian 3 methods × 5 seeds collection, so its own global retained-extremum ranking stays intact.

After all three selected screen batches finish successfully, analyze the actual directories directly:

```bash
# The earlier throughput-only choice fixes these two directory names.
gl_chosen_baseline="$gl_gate_packet/results-screen-throughput${gl_screen_workers}-build5"
if [ "$gl_screen_workers" = 20 ]; then
  gl_other_baseline="$gl_gate_packet/results-screen-throughput30-build5"
else
  gl_other_baseline="$gl_gate_packet/results-screen-throughput20-build5"
fi
python "$gl_gate_packet/analyze_physics_pilots.py" \
  "$gl_gate_packet/manifest-screen-full-build5.json" "$gl_chosen_baseline" \
  --pilot-directory "$gl_gate_packet/results-screen-advanced-firstseed-build5" \
  --pilot-directory "$gl_gate_packet/results-screen-all-remaining-seeds-build5" \
  --workers "$gl_screen_workers" --exclude-pilot-directory "$gl_other_baseline"
```

The script writes `physics_pilot_analysis.json` in the chosen baseline directory and records every absolute input origin. Inspect its `passed` field and retained failures before freezing a/b/c. Do not substitute the other worker result based on its numerical estimate or variance.
