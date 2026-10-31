"""Pure exact applicability guard; no license reads and no tracked writes."""
import ast,datetime,hashlib,json,os,re,subprocess,tomllib
from pathlib import Path
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
F=W.parent/'generation-oracle-numerator-families';P=F/'target/grouped-identical-migration-preparation-20261031-r1'
R=F/'.local/review-grouped-identical-migration-astra-20261031-024225Z-1469f5e4'
BASE='457e361d7d39e0a3daae801651f8c1dfabbf47aa';D0='d0fb4e906140a6a0fc8dd6401ca839cab0c2734d';FFA='ffa76ec420642c03349ab756d2d7720ab93c0057'
PATH='crates/feynkit-generator/tests/generation_oracle/main.rs';KEYS=('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')
assert not any(k in os.environ for k in KEYS)
ENV=dict(os.environ,GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def write(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
def git(*args):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*args],env=ENV)
def tree(rev):
 out={}
 for row in git('ls-tree','-rz','--full-tree',rev).split(b'\0'):
  if row:
   a,p=row.split(b'\t',1);mode,kind,oid=a.decode().split();assert kind=='blob';out[os.fsdecode(p)]={'mode':mode,'blob':oid}
 return out
def scopes(reg):
 out=[]
 for sec in tomllib.loads(reg.decode())['section']:
  for n in sec.get('note',[]):
   for s in n.get('scope',[]):
    b=(W/s['path']).read_bytes();t=b.decode();a,z=s.get('anchor'),s.get('end_anchor');assert not a or t.count(a)==1;assert not z or t.count(z)==1
    r=t[t.index(a):t.index(z)].encode() if z else b;assert blob(r)==s['digest'];out.append({'note':n['id'],**s})
 return out
assert subprocess.check_output(['jj','--ignore-working-copy','log','--no-graph','-r','@','-T','commit_id'],cwd=W,env=ENV,text=True).strip()==BASE
base=tree(BASE)
for p,r in base.items():
 f=W/p;b=os.fsencode(os.readlink(f)) if r['mode']=='120000' else f.read_bytes();assert blob(b)==r['blob'],p
assert re.search(r'^  enable = false;$',(W/'nix-ci.nix').read_text(),re.M)
assert sha(R/'review.json')=='7d96505dcb700eeae42b5a4e87676551ba910cb660df586cafe34d6b90e944f1'
assert sha(R/'seal.json')=='e4f574d2346fd4c32c27fbcfbe5b027acd1c8772185371d83bc701429218a26f'
for p,r in load(R/'seal.json')['files'].items():assert sha(R/p)==r['sha256'] and (R/p).stat().st_size==r['bytes'],p
assert sha(P/'handoff.json')=='6279084309a00acfadede1c160cee1a2715f5b45a2e16e82dfba3470194b099c'
assert sha(P/'seal.json')=='49bfcf3c5b9bf9dd67552b5811498f2c3a80d12c1813c9aba973497ad0635b42'
for line in (P/'SHA256SUMS').read_text().splitlines():
 h,p=line.split('  ',1);assert not Path(p).is_absolute() and '..' not in Path(p).parts;assert sha(P/p)==h,p
assert sha(P/'proposed-NOT-APPLIED.patch')=='2450f56adfe23ea3b157b4d687c4ad6edd39866ff841f6ba8e85fa6f6d864fb7'
old=(W/PATH).read_bytes();assert old==git('show',D0+':'+PATH)==git('show',FFA+':'+PATH)==(P/'main.before.rs').read_bytes()
assert blob(old)=='9ef9f1ba005f20149716b22346d45c5501b62e87' and hashlib.sha256(old).hexdigest()=='6dacf61baca9f13f3e4a6d719f89eaced429fb3cb14c7904dcb41fe4c86b66df'
new=old.replace(b'numerators: 0x00fac78da0e0d159',b'numerators: 0xa60df38a6b8079e3');assert old.count(b'numerators: 0x00fac78da0e0d159')==1
assert new==(P/'main.after.NOT-APPLIED.rs').read_bytes()
assert blob(new)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f' and hashlib.sha256(new).hexdigest()=='66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2'
rows=lambda b:{m[1].decode():m[0] for m in re.finditer(rb'(\w+): Row \{[^\n]*?\}',b)}
a,b=rows(old),rows(new);assert len(a)==len(b)==149
assert [n for n in a if a[n]!=b[n]]==['xs_qed_eemm_identical_L2']
for n in a:assert len(re.findall(rb'\w+:',a[n]))==11 and re.sub(rb'numerators: 0x[0-9a-f]+',b'numerators: TOKEN',a[n])==re.sub(rb'numerators: 0x[0-9a-f]+',b'numerators: TOKEN',b[n])
prior=load(W/'docs/architecture/evidence/generation-integration/feynkit-refresh/raw-three-migration/receipt.json')['runtime']['results'][0]['outcomes'];prior=[r['case'] for r in prior];assert len(prior)==12 and all(a[n]==b[n] for n in prior)
# Complete path/blob equality, including every manifest, lock, filter and source owner.
comparisons=[]
for rev in [D0,FFA]:
 t=tree(rev);assert set(t)==set(base) or all(p.startswith('docs/') for p in set(t)^set(base))
 delta=[p for p in sorted(set(t)&set(base)) if t[p]!=base[p]]
 production=[p for p in delta if not p.startswith('docs/')]
 assert production==['examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi'],production
 comparisons.append({'revision':rev,'path_count':len(t),'compared_paths':len(set(t)&set(base)),'all_deltas':delta,'all_non_docs_deltas':production,'complete_generator_producer_source_invariant':True})
assert sha(W/'examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi')=='5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5'
# Rehash reviewed original source reconstruction inputs and its detailed source owners.
supp=load(R/'supplementary-v3.json');bound=[]
for r in supp['source_bindings']:
 p=r['path'];assert sha(W/p)==r['refresh_current_sha256'];bound.append(r)
# Reconstruct original observational overlay in memory using the review's exact function only.
Q=W/'.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
ns={'Q':Q,'read':lambda p:Path(p).read_bytes(),'need':lambda x,s:None if x else (_ for _ in ()).throw(AssertionError(s)),'hashlib':hashlib,'re':re}
asttree=ast.parse((R/'supplementary-v3.py').read_bytes());fn=next(n for n in asttree.body if isinstance(n,ast.FunctionDef) and n.name=='uninstrumented_sources');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(R/'supplementary-v3.py'),'exec'),ns);originals=ns['uninstrumented_sources']()
domain=['crates/feynkit-generator/','crates/feynkit-graph/','crates/feynkit-model/','crates/idenso/','crates/spenso/','assets/'];orig=[]
for f in sorted((Q/'source').rglob('*')):
 if not f.is_file() or f.is_symlink():continue
 p=f.relative_to(Q/'source').as_posix()
 if not(any(p.startswith(d) for d in domain) or p in ['Cargo.toml','Cargo.lock']):continue
 o=originals[p] if p in originals else f.read_bytes()
 if o is None:continue
 c=(W/p).read_bytes()
 if p==PATH:assert re.sub(rb'numerators: 0x[0-9a-f]+',b'numerators: TOKEN',o)==re.sub(rb'numerators: 0x[0-9a-f]+',b'numerators: TOKEN',c)
 else:assert o==c,p
 orig.append({'path':p,'original_sha256':hashlib.sha256(o).hexdigest(),'accepted_host_sha256':hashlib.sha256(c).hexdigest(),'only_prior_literal_delta':p==PATH})
sc=scopes((W/'docs/developers.toml').read_bytes());assert len(sc)==118
write('source-invariance.json',{'base':BASE,'exact_before_blob':blob(old),'exact_after_blob':blob(new),'all149_Rows_10fields_preserved_except_single_numerator':True,'prior12':prior,'prior12_rows_sha256':{n:hashlib.sha256(a[n]).hexdigest() for n in prior},'complete_tree_comparisons':comparisons,'detailed_review_source_owners':bound,'original_source_domain':orig,'scope_bindings':sc,'review_files_rehashed':len(load(R/'seal.json')['files']),'preparation_SHA256SUMS_rehashed':len((P/'SHA256SUMS').read_text().splitlines()),'stub_sha256':sha(W/'examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi')})
write('start.json',{'time':datetime.datetime.now(datetime.UTC).isoformat(),'base':BASE,'base_tree':base,'refs':git('show-ref').decode(),'root_acceptance_sha256':sha(W/'.local/feynkit-refresh/root-grouped-identical-literal-acceptance-20261031.json'),'tracking_limit':'1MiB','author':'Lucien Huber <im@lcnbr.ch>'})
print('PASS exact before/after,149Rows,prior12,118scopes,complete d0fb/ffa producer invariance and',len(orig),'original source bindings; no tracked edits')
