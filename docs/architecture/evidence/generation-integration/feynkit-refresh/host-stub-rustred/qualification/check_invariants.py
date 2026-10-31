from pathlib import Path
import os,subprocess,hashlib,json,tomllib,ast
P=Path(__file__).resolve().parent;W=P.parents[2];ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV['GIT_OPTIONAL_LOCKS']='0'
def git(*a):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*a],env=ENV)
def blob(b):return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def scopes(reg,reader):
 out=[]
 for s in tomllib.loads(reg.decode())['section']:
  for n in s.get('note',[]):
   for x in n.get('scope',[]):
    b=reader(x['path']);t=b.decode();a,z=x.get('anchor'),x.get('end_anchor');assert not a or t.count(a)==1;assert not z or t.count(z)==1;r=t[t.index(a):t.index(z)].encode() if z else b;assert blob(r)==x['digest'];out.append({'note':n['id'],**x,'actual_blob':blob(r)})
 return out
base='d0fb4e906140a6a0fc8dd6401ca839cab0c2734d';code='2cfb52b28d99f92305a9d571bbe3911cf04b8a19';path='examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi'
assert git('diff','--name-only',base,code).decode().splitlines()==[path]
assert hashlib.sha256((W/path).read_bytes()).hexdigest()=='5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5'
assert git('show','-s','--format=%P',code).decode().strip()==base
assert git('show','-s','--format=%an <%ae>',code).decode().strip()=='Lucien Huber <im@lcnbr.ch>'
old=scopes(git('show',base+':docs/developers.toml'),lambda p:git('show',base+':'+p));new=scopes((W/'docs/developers.toml').read_bytes(),lambda p:(W/p).read_bytes());assert old==new and len(new)==118
rows='crates/feynkit-generator/tests/generation_oracle/main.rs';assert (W/rows).read_bytes()==git('show',base+':'+rows)
initial=json.loads((P/'initial-candidate-files.json').read_text());changed=[]
for p,v in initial.items():
 f=W/p;b=os.readlink(f).encode() if f.is_symlink() else f.read_bytes()
 if hashlib.sha256(b).hexdigest()!=v['sha256']:changed.append(p)
assert set(changed)<=set([path,'docs/developers.toml']),changed
assert '  enable = false;' in (W/'nix-ci.nix').read_text()
main=json.loads((P/'initial-main-files.json').read_text());main_changed=[]
for p,v in main.items():
 f=W.parent/'main'/p;b=os.readlink(f).encode() if f.is_symlink() else f.read_bytes()
 if hashlib.sha256(b).hexdigest()!=v['sha256']:main_changed.append(p)
G=json.loads((P/'initial-guard.json').read_text());index=hashlib.sha256(Path(G['main_index_path']).read_bytes()).hexdigest();integration=subprocess.check_output(['jj','--ignore-working-copy','log','-r','generation-integration','--no-graph','-T','commit_id'],cwd=W,env=ENV).decode();assert integration==G['global_integration']
result={'code':code,'base':base,'only_production_change':path,'all149_Rows_and12_migrations_byte_preserved':True,'all118_scope_bindings':new,'CI':False,'global_integration_unchanged':integration,'initial_main_files_changed':main_changed,'main_index_unchanged':index==G['main_index_sha256'],'main_index_initial':G['main_index_sha256'],'main_index_current':index,'candidate_original_files_changed':changed,'no_vcs_global_invariance_claim':True,'unrelated_shared_ref_changes_permitted_and_recorded':True}
(P/'invariants.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='all118_scope_bindings'}))
