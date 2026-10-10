from pathlib import Path
import gzip, hashlib, json, shutil

root = Path('/tmp/soft-ct-rebase-2026-09-23')
repo = Path('/common/dev/gammaloop/lcnbr')
dest = root / 'evidence'
dest.mkdir(parents=True, exist_ok=True)
labels = ['single-u', 'single-h', 'orientation-58-u', 'orientation-58-h', 'all-bare', 'all-u', 'all-h', 'all-h-projected']
summary = []
for label in labels:
    folder = root / label
    if not (folder / 'result.json').exists():
        continue
    receipt = json.loads((folder / 'result.json').read_text())
    trace = folder / 'state/logs/gammalog-trace.jsonl'
    rows = [json.loads(line) for line in trace.read_text().splitlines()] if trace.exists() else []
    summary.append({
        **{k: receipt[k] for k in ['label', 'elapsed_seconds', 'exit_code', 'termination_reason', 'peak_sampled_tree_rss_kib', 'peak_observed_process_hwm_kib', 'state_saved']},
        'native_source_maps': [row['native_source_maps'] for row in rows if row.get('stage') == 'raw_cff_generation'],
        'evaluator_stages': [row for row in rows if row.get('stage') in ['evaluator_stack_numerator_components_reachable', 'evaluator_stack_parse_atoms_done', 'evaluator_stack_new_single_parametric_done', 'evaluator_stack_new_done']],
        'last_stage': next((row for row in reversed(rows) if row.get('stage')), None),
    })
    for name in ['result.json', 'command.json', 'resources.jsonl']:
        shutil.copyfile(folder / name, dest / (label + '-' + name))
    shutil.copyfile(folder / 'stdout.log', dest / (label + '-stdout.log.txt'))
    if trace.exists():
        with trace.open('rb') as source, gzip.open(dest / (label + '-trace.jsonl.gz'), 'wb') as target:
            shutil.copyfileobj(source, target)
profiles = []
for label in labels:
    folder = root / (label + '-profile')
    if not (folder / 'result.json').exists():
        continue
    receipt = json.loads((folder / 'result.json').read_text())
    profiles.append(receipt)
    for name in ['result.json', 'command.json', 'resources.jsonl']:
        shutil.copyfile(folder / name, dest / (label + '-profile-' + name))
    shutil.copyfile(folder / 'stdout.log', dest / (label + '-profile-stdout.log.txt'))
    profile = folder / 'profile/uv_profile.json'
    if profile.exists():
        with profile.open('rb') as source, gzip.open(dest / (label + '-uv-profile.json.gz'), 'wb') as target:
            shutil.copyfileobj(source, target)
(dest / 'profile-results.json').write_text(json.dumps(profiles, indent=2) + '\n')
for name in ['single.toml', 'orientation-58.toml', 'all.toml', 'all-projected.toml', 'muv.json', 'bare.json', 'rerun-generation.py', 'profile-generated.py', 'collect-evidence.py', 'binary-runtime-license-sha256.txt']:
    shutil.copyfile(root / name, dest / name)
shutil.copyfile(root / 'build-runtime-license.log', dest / 'build-runtime-license.log.txt')
manifest = json.loads((root / 'source-environment.json').read_text())
manifest['runtime_license_build'] = {'NO_SYMBOLICA_OEM_LICENSE': '1', 'SYMBOLICA_LICENSE': 'set in ignored local .envrc; value omitted'}
for file in manifest['source_sha256']:
    assert hashlib.sha256((repo / file).read_bytes()).hexdigest() == manifest['source_sha256'][file], file
(dest / 'source-environment.json').write_text(json.dumps(manifest, indent=2) + '\n')
(dest / 'generation-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
for result in summary:
    print(result['label'], round(result['elapsed_seconds'], 3), round(result['peak_observed_process_hwm_kib'] / 1024**2, 3), result['exit_code'], result['termination_reason'])
