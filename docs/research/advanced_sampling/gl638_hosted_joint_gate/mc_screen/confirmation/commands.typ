= Confirmation and retained maxima after the fixed screen
<confirmation-and-retained-maxima-after-the-fixed-screen>
Prepared recipe only. No a/b/c selection or confirmation has been performed by this preparation. The existing screen analyzer owns selection. Workers are already frozen at 20 by the preregistered throughput/RSS decision; the excluded 30-worker baseline never enters the screen or its maximum replay.

`confirmation-build5-template.json` records the fixed overrides and source hashes; it is not a runnable Gamma manifest. `maximum-replay-build5-commands.json` contains exact arguments for the four later maximum commands. All commands reuse optimized build5, the frozen b53 physical client and max29b client. The original seven cards, saved state, full 936 orientations/six cuts, UV/threshold settings, physical tolerances and accepted seven-catalogue preflight remain unchanged. Only confirmation N/seeds/mode selection and execution metadata differ.

== Freeze the screen before maximum replay
<freeze-the-screen-before-maximum-replay>
After all three included batches finish, run the existing screen analysis exactly once:

```bash
source /tmp/gammaloop-rebase-dev-env.sh
gl_gate_packet=/tmp/gl638-hosted-joint-gate
python "$gl_gate_packet/analyze_physics_pilots.py" \
  "$gl_gate_packet/manifest-screen-full-build5.json" \
  "$gl_gate_packet/results-screen-throughput20-build5" \
  --pilot-directory "$gl_gate_packet/results-screen-advanced-firstseed-build5" \
  --pilot-directory "$gl_gate_packet/results-screen-all-remaining-seeds-build5" \
  --workers 20 --exclude-pilot-directory "$gl_gate_packet/results-screen-throughput30-build5"
```

Stop on a failed analysis and retain its failures. The included matrix is exactly baseline 20 at 11113, the six advanced modes at 11113, and all seven at 22229/33331:21 runs, 43,008 draws. There is no duplicated first baseline or synthetic combined workspace.

The following metadata-only command takes the analyzer\'s selected order unchanged and writes new files exclusively. It neither selects a new candidate nor evaluates physics. Run it only after the complete screen analysis passes:

```bash
python - <<'PY'
from pathlib import Path
import hashlib
import json
p = Path('/tmp/gl638-hosted-joint-gate')
t = json.loads((p / 'confirmation-build5-template.json').read_text())
base = Path(t['base_manifest'])
assert hashlib.sha256(base.read_bytes()).hexdigest() == t['base_manifest_sha256']
for name, expected in t['required_input_identities'].items():
    assert hashlib.sha256((p / name).read_bytes()).hexdigest() == expected, name
assert json.loads((p / 'screen-worker-choice-build5.json').read_text())['chosen_workers'] ==20
report_path = Path(t['required_screen_report'])
report_bytes = report_path.read_bytes()
screen = json.loads(report_bytes)
assert screen['passed'] and not screen['failures'] and screen['workers'] ==20
assert screen['manifest_sha256'] == t['base_manifest_sha256']
selection = screen['selection']
modes = selection['confirmation_modes']
assert len(modes) == len(set(modes)) ==3 and modes[0] == 'optimized_lmb'
assert modes[-1] == selection['primary_estimate_mode']
manifest = json.loads(base.read_text())
assert selection['candidate_driver_sha256'] == manifest['candidate_driver_sha256']
for mode in modes:
    card = Path(manifest['cards'][mode])
    assert hashlib.sha256(card.read_bytes()).hexdigest() == selection['card_sha256'][mode]
    assert manifest['card_sha256'][str(card)] == selection['card_sha256'][mode]
frozen_path = Path(t['frozen_screen_copy'])
out = Path(t['output_manifest'])
assert not frozen_path.exists() and not out.exists()
manifest.update(t['overrides'])
manifest.update({
    'status': 'prepared_confirmation_not_run',
    'proposal_order': modes,
    'batch_purpose': 'one complete frozen3method x5seed independent confirmation; no new winner selection',
    'execution_prerequisite': 'Complete valid frozen screen plus same accepted build5 source/state/card preflight; physics-pilot enforces preflight reuse',
    'gates_pending': ['Independent fixed confirmation and retained maximum replay/accuracy audit'],
    'confirmation_screen_report': str(frozen_path),
    'confirmation_screen_report_sha256': hashlib.sha256(report_bytes).hexdigest(),
    'confirmation_frozen_selection': selection,
})
with frozen_path.open('xb') as f:
    f.write(report_bytes)
with out.open('x') as f:
    f.write(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'manifest': str(out), 'sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
                  'frozen_screen': str(frozen_path), 'modes': modes, 'workers':20,
                  'samples_per_method':163840, 'total_samples':491520}))
PY
```

Review and retain the printed hashes before running confirmation. A failure to select c does not authorize choosing it from local-ray plateaus or from a partial screen.

== One complete confirmation collection
<one-complete-confirmation-collection>
```bash
RAYON_NUM_THREADS=20 GL_DISPLAY_FILTER=off \
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-confirmation-build5.resources.txt" \
  "$gl_gate_packet/drivers/gate-physics-build5" \
  "$gl_gate_packet/manifest-confirmation-build5.json" \
  "$gl_gate_packet/results-confirmation-build5" all physics-pilot 20 1 1 1337 \
  > "$gl_gate_packet/results-confirmation-build5.log" 2>&1
```

The driver applies N=32768 to each of seeds 10007,20011,30013,40009,50021, one fresh untrained iteration and batch\_size=256. The collection contains 15 runs and 491,520 draws, 163,840 per method. Physical failures remain failures; no diagnostic-pilot mode or selective extra samples are used. The primary result is the screen-selected c even if confirmation does not favor it.

After the terminal successful confirmation, write its statistics before hashing the pilot for maximum replay:

```bash
python "$gl_gate_packet/analyze_physics_pilots.py" \
  "$gl_gate_packet/manifest-confirmation-build5.json" \
  "$gl_gate_packet/results-confirmation-build5" \
  --workers 20 --screen "$gl_gate_packet/screen-selection-frozen-build5.json"
```

Require `physics_pilot_analysis.json` to pass. The existing analyzer authenticates the frozen screen/card/driver/workers and uses its absolute-phase scales unchanged; it does not choose a confirmation winner.

== Exact maximum replay and independent accuracy audit
<exact-maximum-replay-and-independent-accuracy-audit>
`max-replay-build5` accepts exactly `MANIFEST PILOT_DIRECTORY NEW_OUTPUT_DIRECTORY`. It reads the original integration\_state.bin complete Samples, preserves all phase/sign aliases and uses first-iteration max\_eval0. Its single 128MB-stack evaluator worker may run with `RAYON_NUM_THREADS=1`. Native returned events and detailed CT terms already include the outer Sample factor.

The three screen inputs and the later confirmation have these exact mappings:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,auto,auto,),
    table.header([Manifest], [Immutable pilot directory], [New maximum output],),
    table.hline(),
    [manifest-screen-throughput20-build5.json], [results-screen-throughput20-build5], [results-max-screen-throughput20-build5],
    [manifest-screen-advanced-firstseed-build5.json], [results-screen-advanced-firstseed-build5], [results-max-screen-advanced-firstseed-build5],
    [manifest-screen-all-remaining-seeds-build5.json], [results-screen-all-remaining-seeds-build5], [results-max-screen-all-remaining-seeds-build5],
    [manifest-confirmation-build5.json], [results-confirmation-build5], [results-max-confirmation-build5],
  )]
  , kind: table
  )

All paths are under `/tmp/gl638-hosted-joint-gate`; exact unambiguous argv arrays are in `maximum-replay-build5-commands.json`. For example, the first batch is:

```bash
RAYON_NUM_THREADS=1 GL_DISPLAY_FILTER=off \
/run/current-system/sw/bin/time -v -o "$gl_gate_packet/results-max-screen-throughput20-build5.resources.txt" \
  "$gl_gate_packet/drivers/max-replay-build5" \
  "$gl_gate_packet/manifest-screen-throughput20-build5.json" \
  "$gl_gate_packet/results-screen-throughput20-build5" \
  "$gl_gate_packet/results-max-screen-throughput20-build5" \
  > "$gl_gate_packet/results-max-screen-throughput20-build5.log" 2>&1
python "$gl_gate_packet/audit_maximum_results.py" \
  "$gl_gate_packet/results-max-screen-throughput20-build5"
```

Require maximum `summary.json` status completed, no error, unchanged state and pilot hashes; then require `maximum_accuracy_audit.json` passed. The audit reuses the existing physical complex-norm/cut metric, checks complete-factor/event/decomposition identities and ordinary-versus-forced-Arb controls. The driver alone establishes reporting replay, validity and six-cut coverage; it does not replace the independent numerical comparison.

Screen top-two controls are per method within each real batch, among retained signed-extremum Samples; they are not a global complex-norm ranking across the three directories. Confirmation\'s complete 3×5 collection provides a single retained-extremum ranking across its five seeds. Neither is the unrecorded maximum complex norm over all samples or a global boundedness certificate. Detailed returned CT decomposition and the existing metadata registry support component attribution; complete native star geometry can require separately retained existing traces.

== Safe scheduling and immutable directory hashes
<safe-scheduling-and-immutable-directory-hashes>
Maximum replay hashes its entire pilot input directory before and after evaluation. Therefore statistics, selection copies and any other writes inside that pilot directory must finish first. Do not rerun its analyzer, add logs or write audit products there while maximum replay is active. Logs, resource records, maximum outputs and their accuracy audits use the separate sibling paths above. Archiving may read these immutable inputs and write elsewhere.

After the screen analysis and selection are frozen, one screen-max process can overlap the 20-worker confirmation because they use separate read-only loaded-state owners and disjoint writable outputs. Run only one maximum replay at a time and observe combined peak resident memory≤300GB. Nominal evaluation parallelism is 21, within 30. Keep the shared generated state immutable and retain each driver\'s before/after checks. This overlap is for physics throughput only; its timings do not support an isolated scaling or performance comparison. Confirmation maximum replay starts only after confirmation statistics have been written and frozen.
