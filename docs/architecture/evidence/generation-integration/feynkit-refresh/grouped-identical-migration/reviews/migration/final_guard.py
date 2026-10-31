from pathlib import Path
import ast,datetime,hashlib,json
O=Path(__file__).resolve().parent;W=O.parent.parent;P=W/'target/grouped-identical-migration-preparation-20261031-r1'
def load(p):return json.loads(Path(p).read_text())
def h(p):
    hh=hashlib.sha256()
    with Path(p).open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):hh.update(b)
    return hh.hexdigest()
first=load(O/'integrity-validated-before.json');last=load(O/'integrity-final.json')
assert not first['mismatches']and not last['mismatches']
assert load(O/'input-hashes-before.json')==load(O/'input-hashes-validated-before.json')==load(O/'input-hashes-final.json')
inputs={}
for f in O.glob('*reads.json'):
    for p,d in load(f).items():
        if not Path(p).is_relative_to(O):
            assert p not in inputs or inputs[p]==d,(p,'changed between fresh reads')
            inputs[p]=d
changes=[]
for p,d in inputs.items():
    actual=h(p)
    if actual!=d:changes.append({'path':p,'before':d,'after':actual})
assert not changes,changes
rows=load(P/'preserved-prior12-Rows.json');a=(P/'main.after.NOT-APPLIED.rs').read_bytes();b=(P/'main.before.rs').read_bytes()
assert len(rows)==12
for name,r in rows.items():assert r['row'].encode()in a and r['row'].encode()in b
parsed=[]
for f in O.glob('*.py'):ast.parse(f.read_text());parsed.append(f.name)
registries=[{'path':p,'before':d,'after':last['registry'][p]}for p,d in first['registry'].items()if last['registry'][p]!=d]
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protected_files_unchanged':15899,'fresh_external_read_bindings':len(inputs),'fresh_external_changes':changes,'candidate_files_unchanged':71,'prior12_rows_exact':True,'registry_movements_during_review':registries,'preparation_reported_external_ref_movement':load(P/'workspace-guard.json')['external_VCS_movement'],'preparation_reported_registry_movement':load(P/'workspace-guard.json')['refresh_external_tracked_changes'],'ref_index_operations_in_this_review':0,'fresh_general_ref_invariance_claim':False,'Python_AST_parsed':parsed,'owned_processes_remaining':[],'no_restore_or_cleanup_of_external_work':True}
(O/'final-guard.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
print(json.dumps(r,indent=2))
