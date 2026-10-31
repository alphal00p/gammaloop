"""Final read-only guards and noncircular packet seal; writes only this directory."""
import ast
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

O=Path(__file__).resolve().parent
W=O.parents[1]
assert O==W/'target/grouped-identical-migration-preparation-20261031-r1'
assert not any('LICENSE' in k for k in os.environ)

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def load(n):return json.loads((O/n).read_text())
def save(n,x):(O/n).write_text(json.dumps(x,indent=2)+'\n')

def git(*args):
    result=subprocess.run(['/home/lcnbr/.nix-profile/bin/git','--no-optional-locks','--git-dir='+str(W.parent/'main/.git'),*args],capture_output=True,env={'LANG':'C.UTF-8','GIT_OPTIONAL_LOCKS':'0','GIT_NO_LAZY_FETCH':'1'})
    with (O/'vcs-reads.jsonl').open('a') as f:f.write(json.dumps({'args':args,'exit':result.returncode,'stderr':result.stderr.decode(),'stdout_sha256':hashlib.sha256(result.stdout).hexdigest()})+'\n')
    assert result.returncode==0
    return result.stdout

# Newly consumed accepted characterization review: separate observation-only guard.
D=W.parent/'generation-feynkit-refresh/.local/review-full-generator-characterization-astra-20261031-022430Z-320628f7'
s=json.loads((D/'SEAL.json').read_text())
assert sha(D/'review.json')=='fb688585bc7efd0b494fef7648df5abad398f88bfe8c439e8b42c0dd520eaa37'
assert sha(D/'SEAL.json')=='5e4792d8c2a2e678b21f92b3e3e4f28348a9ca8e30d7b8b6edf4d751dc769776'
assert sha(D/'SHA256SUMS')==s['manifest_sha256']
checked=[]
for line in (D/'SHA256SUMS').read_text().splitlines():
    h,path=line.split('  ',1)
    p=D/path
    assert p.resolve().is_relative_to(D)
    assert sha(p)==h
    checked.append({'path':str(p),'sha256':h})
assert len(checked)==44
save('accepted-characterization-review-guard.json',{'status':'PASS','files':checked,'approval_scope':'Root accepted asserted-summary/pool1 characterization only, with deleted payloads and no pool8/32/mathematical/migration credit.'})

# Snapshot-free working-file and shared metadata comparison, recording outside movement.
before=load('workspace-before.json')
g=load('workspace-guard.json')
for key,root in [('feature_tracked_changes',W),('refresh_external_tracked_changes',W.parent/'generation-feynkit-refresh')]:
    base=before['feature' if root==W else 'refresh'];changes=[]
    for path,expected in base.items():
        p=root/path
        actual={'symlink':os.readlink(p)} if p.is_symlink() else sha(p) if p.is_file() else None
        if expected!=actual:changes.append({'path':path,'before':expected,'after':actual})
    g[key]=changes
assert not g['feature_tracked_changes']
now={'main_head':git('rev-parse','HEAD').decode().strip(),'refs_sha256':hashlib.sha256(git('show-ref')).hexdigest(),'index_sha256':sha(W.parent/'main/.git/index')}
g['shared_metadata_after']=now
g['external_VCS_movement']={k:before[k]!=v for k,v in now.items()}
g['final_guard_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
save('workspace-guard.json',g)

# Preserve the exact patch and all byte guards, without git apply or object writes.
b=(O/'main.before.rs').read_bytes();a=(O/'main.after.NOT-APPLIED.rs').read_bytes()
p=load('patch-binding.json')
blob=lambda data:hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
assert blob(b)==p['before_blob'] and blob(a)==p['after_blob']
assert hashlib.sha256(b).hexdigest()==p['before_sha256']
assert hashlib.sha256(a).hexdigest()==p['after_sha256']
assert sha(O/'proposed-NOT-APPLIED.patch')==p['patch_sha256']
assert a.replace(('numerators: '+p['new']).encode(),('numerators: '+p['old']).encode(),1)==b
assert git('rev-parse',p['immutable_base']+':'+p['path']).decode().strip()==p['before_blob']
assert len(re.findall(rb'^\s+(?:fast|heavy) \w+: Row \{',a,re.M))==149
for rec in load('preserved-prior12-Rows.json').values():assert rec['row'].encode() in a
assert load('dependency-gap-matrix.json')['Auto_aggregation_gate']=='ACCEPTED_PINNED_CONDITIONALLY_BY_ROOT'
assert load('dependency-gap-matrix.json')['status']=='NEW_MIGRATION_REVIEW_REQUIRED'

# New Python is parsed, Typst links checked. No application/compiler/native import.
parsed=[]
for f in O.glob('*.py'):ast.parse(f.read_text(),filename=str(f));parsed.append(f.name)
t=(O/'handoff.typ').read_text()
assert len(re.findall(r'^= ',t,re.M))==1
links=re.findall(r'#link\("([^\"]+)"\)',t)
assert all((O/x).exists() or x in ['manifest.json','seal.json'] for x in links)
save('static-checks.json',{'Python_AST_parsed':parsed,'Typst_single_title':True,'Typst_relative_links_checked':links,'unrendered':True,'no_linter_compiler_or_native_checks':True})
h=load('handoff.json');h['external_VCS_movement']=g['external_VCS_movement'];h['refresh_external_changes']=g['refresh_external_tracked_changes'];h['characterization_review_files_rehashed']=44;h['guard_utc']=g['final_guard_utc'];save('handoff.json',h)

# No further packet writes after the seal and its own sidecar.
excluded={'manifest.json','SHA256SUMS','seal.json','seal.sha256'}
files={str(f.relative_to(O)):{'sha256':sha(f),'bytes':f.stat().st_size} for f in sorted(O.rglob('*')) if f.is_file() and f.name not in excluded}
save('manifest.json',{'schema':'grouped-identical-preparation-manifest/1','files':files,'self_exclusions':sorted(excluded),'status':'NEW_MIGRATION_REVIEW_REQUIRED','not_applied':True})
manifest_sha=sha(O/'manifest.json')
(O/'SHA256SUMS').write_text(''.join(r['sha256']+'  '+f+'\n' for f,r in files.items())+manifest_sha+'  manifest.json\n')
save('seal.json',{'schema':'grouped-identical-preparation-seal/1','status':'NEW_MIGRATION_REVIEW_REQUIRED','manifest_sha256':manifest_sha,'SHA256SUMS_sha256':sha(O/'SHA256SUMS'),'handoff_sha256':sha(O/'handoff.json'),'patch_sha256':p['patch_sha256'],'files':len(files),'before_blob':p['before_blob'],'after_blob':p['after_blob'],'approval':False,'applied':False,'Auto_accepted_conditionally_by_root':True,'None_accepted_separately_conditionally_by_root':True,'self_exclusions':'manifest.json and SHA256SUMS bind ordinary files; seal.json binds both. seal.sha256 binds seal.json; no circular self-hash.','owned_live_processes':[]})
seal_sha=sha(O/'seal.json')
(O/'seal.sha256').write_text(seal_sha+'  seal.json\n')
assert all(sha(O/f)==r['sha256'] for f,r in files.items())
print(json.dumps({'status':'SEALED_NOT_APPLIED','packet':str(O),'files':len(files),'seal_sha256':seal_sha,'new':p['new'],'before_blob':p['before_blob'],'after_blob':p['after_blob'],'remaining':'NEW_MIGRATION_REVIEW_REQUIRED','external_refresh_changes':g['refresh_external_tracked_changes'],'external_VCS_movement':g['external_VCS_movement'],'owned_live_processes':[]}))
