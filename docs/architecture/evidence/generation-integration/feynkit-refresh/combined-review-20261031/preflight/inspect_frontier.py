"""Pure, local demand-frontier inspection; never realize or execute an artifact."""
from pathlib import Path
exec(Path(__file__).with_name('evidence_io.py').read_text())

def norm(value):
    return value if value.startswith('/') else '/nix/store/' + value

core = json.loads((P/'code2cf-core-only-eval.json').read_text())['stdout']
groups_record = json.loads((P/'code2cf-scoped-groups-corrected-eval.json').read_text())
groups = json.loads(groups_record['stdout']) if groups_record.get('exit_code') == 0 else {}
group_drvs = {'core': core, **groups}
raw_groups = {}
for group, drv in group_drvs.items():
    r = run('frontier-group-' + group, ['nix', '--offline', '--option', 'allow-import-from-derivation', 'false', 'derivation', 'show', drv])
    assert r['exit_code'] == 0
    obj = json.loads(r['stdout']); obj = obj.get('derivations', obj)
    raw_groups[group] = (drv, next(iter(obj.values())))

# Select normal package-specific owners; never build entire FeynKit/Spenso/docs groups.
wanted = {
    'core': ['gammaloop-api', 'gammalooprs', 'gammaloop-tracing-filter', 'gammaloop-tracing-filter-macros', 'three-dimensional-reps'],
    'integration': ['gammaloop-integration-tests'],
    'feynkit': ['feynkit-graph', 'feynkit-model', 'feynkit-tensor'],
    'spenso': ['idenso', 'spenso', 'spenso-hep-lib'],
    'docs': ['alphal00p-docs-builder'],
}
roots = {}
for group, packages in wanted.items():
    if group not in raw_groups:
        continue
    drv, d = raw_groups[group]
    deps = d.get('inputs', {}).get('drvs', d.get('inputDrvs', {}))
    for package in packages:
        matches = [norm(k) for k in deps if k.endswith('-gammaloop-nextest-binaries-' + group + '-' + package + '-0.1.0.drv')]
        assert len(matches) == 1, (group, package, matches)
        roots[group + '/' + package] = matches[0]

graph = {}
required = {v: {'out'} for v in roots.values()}
todo = set(required)
missing, reused, unknown = {}, {}, {}
step = 0
allowed = ['name','pname','src','cargoArtifacts','cargoExtraArgs','buildPhase','buildPhaseCargoCommand','installPhase','installPhaseCommand','doCheck','doInstallCargoArtifacts','CARGO_BUILD_INCREMENTAL','RUSTFLAGS','CARGO_PROFILE','cargoProfile','nativeBuildInputs','buildInputs','out','stdenv','NO_SYMBOLICA_OEM_LICENSE','PYO3_PYTHON']
while todo:
    batch = sorted(todo); todo = set(); step += 1
    assert step <= 40, 'bounded demand traversal exceeded 40 rounds'
    r = run('frontier-demand-' + str(step), ['nix','--offline','--option','allow-import-from-derivation','false','derivation','show',*batch])
    assert r['exit_code'] == 0
    obj = json.loads(r['stdout']); obj = obj.get('derivations', obj)
    for key, d in obj.items():
        key = norm(key)
        e = d.get('env', {})
        assert not any(k in e for k in ['SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'])
        outputs = {k: norm(v['path']) for k,v in d['outputs'].items() if v.get('path')}
        deps = {norm(k): set(v.get('outputs', []) if isinstance(v, dict) else v) for k,v in d.get('inputs',{}).get('drvs',d.get('inputDrvs',{})).items()}
        graph[key] = {'outputs':outputs,'inputDrvs':{k:sorted(v) for k,v in deps.items()},'inputSrcs':[norm(s) for s in d.get('inputs',{}).get('srcs',d.get('inputSrcs',[]))], 'builder':d.get('builder'),'args':d.get('args'),'env':{k:e[k] for k in allowed if k in e},'drv_sha256':hashlib.sha256(Path(key).read_bytes()).hexdigest()}
        absent = False
        for output in required[key]:
            slot = key + '!' + output
            path = outputs.get(output)
            if path is None:
                unknown[slot] = {'drv':key,'output':output}; absent = True
            elif Path(path).exists():
                reused[slot] = {'drv':key,'output':output,'path':path}
            else:
                missing[slot] = {'drv':key,'output':output,'path':path}; absent = True
        if absent:
            for dep, outs in deps.items():
                old = required.setdefault(dep,set())
                if outs - old:
                    old.update(outs); todo.add(dep)

paths = sorted({v['path'] for v in reused.values()})
r = run('frontier-reused-path-info', ['nix','--offline','--option','allow-import-from-derivation','false','path-info','--json',*paths])
assert r['exit_code'] == 0
root_records = {}
for name, drv in roots.items():
    d = graph[drv]; root_records[name] = {'drv':drv,'output':d['outputs']['out'],'available':Path(d['outputs']['out']).exists(),'env':d['env'],'drv_sha256':d['drv_sha256']}
(P/'producer-roots.json').write_text(json.dumps(root_records,indent=2)+'\n')
(P/'demand-frontier.json').write_text(json.dumps({'classification':'PURE_LOCAL_STORE_METADATA_ONLY','roots':roots,'examined_derivations':len(graph),'missing_required_outputs':list(missing.values()),'reused_frontier_outputs':list(reused.values()),'unknown_required_outputs':list(unknown.values()),'reused_outputs_registered_valid':True,'network_cache_availability':'not queried; no network authorized','build_time_disk_memory':'not measured; graph size is not a measurement','graph':graph},indent=2)+'\n')
print(json.dumps({'roots':len(roots),'missing_roots':sum(not x['available'] for x in root_records.values()),'examined':len(graph),'missing_outputs':len(missing),'reused_outputs':len(reused),'unknown':len(unknown)}))
