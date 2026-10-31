"""Exact authorized applicability and preserved-input checks. No private config access."""
import datetime,hashlib,json,os,re,subprocess,tomllib
from pathlib import Path
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh'); O=Path(__file__).resolve().parent
F=W.parent/'generation-oracle-numerator-families'; P=F/'target/raw-next-qcd-migration-preparation-20261031-r1'
BASE='83276dc047815d1dcc7f003aeb512c204b26758e'; START='bd2642ac1a5cfd9d9bae35f5ada68360b62a2535'
FILE='crates/feynkit-generator/tests/generation_oracle/main.rs'
KEYS=('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')
ENV={k:v for k,v in os.environ.items() if k not in KEYS};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def write(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
def git(*a):return subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),*a],env=ENV)
def jj(*a):return subprocess.check_output(['jj','--ignore-working-copy',*a],cwd=W,env=ENV,text=True).strip()
def tree(rev):
 out={}
 for row in git('ls-tree','-rz','--full-tree',rev).split(b'\0'):
  if row:
   a,p=row.split(b'\t',1);mode,kind,oid=a.decode().split();assert kind=='blob';out[os.fsdecode(p)]={'mode':mode,'blob':oid}
 return out
pins={
 W/'.local/feynkit-refresh/root-next-qcd-migration-acceptance-and-application-authorization-20261031.json':'91b7937ed1bc94e55f56195656f7b26865b6bc508cb82e8d7047e01e3abb1ed1',
 W/'.local/feynkit-refresh/root-next-qcd-scientific-interface-acceptance-20261031.json':'225b2c244ba02f16ca3ee73a59d23df3fff12204deef44f3aa405397247a0ea7',
 W/'.local/feynkit-refresh/root-next-qcd-migration-preparation-receipt-20261031.json':'4c66d289489823726da2bf9355c43a937830a970e0e1e53e18f34e4a530f8804',
 P/'SEAL.json':'b4e75df006b68563219947ecbdef5cb8fd289f9b4c7ba840d79c0ce599377324',
 P/'handoff.json':'1b69607db0184093de3b16d32541d6f7a5a3f00968b85694685a6d17fb1506b3',
 P/'proposed-not-applied.patch':'1df84cbea910df5ff092203ff6050e785bccf0027bd627b91a601d1f7b98cf84',
 F/'.local/review-next-qcd-migration-astra-20261031-efvezb1p/review.json':'132a73c7a6711a7d4f74bbe47b0427e5e8fc2ddcb58079ec4f686b78bfc39f74',
 F/'.local/review-next-qcd-migration-astra-20261031-efvezb1p/SEAL.json':'c66f4a1039b9c9fad744633c99d99d05f3b57c30db82a3e101f7fc04b9ee8e25',
 F/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x/review.json':'6d892bd9263a7c093bba5fd372fb4663a5a4e81dfaad4bfc709bbba052b96d21',
 F/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x/SEAL.json':'888b82c478697e917cb243587653cc16d09a484f299a671c5ffd8981e58c1d99',
 F/'target/raw-next-qcd-proof-native-execution-20261031-r1/report.json':'3baf87d1d031f790b753bb944d2a4c180f0620297b7808ad348dfa7ea65c76f0',
 F/'target/raw-next-qcd-proof-native-execution-20261031-r1/COMPLETE-SEAL.json':'22b11020019c4b9f8e058a2951de177f5a6f111fbb4f6976932f98ef0c399626',
 F/'target/raw-next-qcd-captures-20261031-r2/handoff.json':'9da40bc2d34381c7c6b67fa81a6d91032ff699222349f6910f9787a18cc7e660',
 F/'target/raw-next-qcd-captures-20261031-r2/SEAL.json':'166dd40ee056d9bf67da1a84eef200dc5bc9ff0fea55e15a88e9feb8f4e98d4b',
 F/'.local/review-next-qcd-captures-astra-20261031-IYGt5Jrv/review.json':'7ccf8f133bff073e5f428d1fabc0ebd6552aa1e04a8b7edd149f526b37330aa7',
 W.parent/'generation-integration-plan/docs/architecture/generation-integration-plan.typ':'4489dc7383468cf2daf7f8c177d3ea87fe9cad45c99af1afc0f160b202ff06f6',
}
checked={}
def check(p,h):
 p=Path(p);actual=sha(p);assert actual==h,(str(p),actual,h);checked[str(p)]={'sha256':actual,'bytes':p.stat().st_size}
for p,h in pins.items():check(p,h)
manifests=[(P/'SHA256SUMS',122),(F/'.local/review-next-qcd-migration-astra-20261031-efvezb1p/SHA256SUMS',74),(F/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x/SHA256SUMS',57),(F/'target/raw-next-qcd-proof-native-execution-20261031-r1/COMPLETE-SHA256SUMS',232),(F/'target/raw-next-qcd-proof-native-execution-20261031-r1/EXECUTION-SHA256SUMS',135),(F/'target/raw-next-qcd-captures-20261031-r2/SHA256SUMS',1976)]
manifest_records=[]
for m,count in manifests:
 lines=m.read_text().splitlines();assert len(lines)==count
 for line in lines:
  h,p=line.split('  ',1);p=Path(p);check(p if p.is_absolute() else m.parent/p,h)
 manifest_records.append({'path':str(m),'sha256':sha(m),'entries':count})
external=load(P/'complete-input-manifest.json')['external_files'];assert len(external)==3869
for p,r in external.items():check(p,r['sha256']);assert Path(p).stat().st_size==r['bytes']
assert jj('log','--no-graph','-r','@','-T','commit_id')==START
assert jj('log','--no-graph','-r','@-','-T','commit_id')==BASE
assert jj('log','--no-graph','-r','generation-integration','-T','commit_id')==BASE
base=tree(BASE);assert git('rev-parse',BASE+'^{tree}').decode().strip()=='c9de5c3fd962bd5d7c9375ab1f86f75699e8a69b'
assert base==tree(START)
for p,r in base.items():
 f=W/p;b=os.fsencode(os.readlink(f)) if r['mode']=='120000' else f.read_bytes();assert blob(b)==r['blob'],p
old=(W/FILE).read_bytes();new=(P/'after.NOT-APPLIED.rs').read_bytes();assert old==git('show',BASE+':'+FILE)==(P/'before.rs').read_bytes()
assert blob(old)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f' and blob(new)=='cc4505a6c4f6fc1687c075ed79e154a8e19f05a6'
rows=lambda b:{m[1].decode():m[0] for m in re.finditer(rb'(\w+): Row \{[^\n]*?\}',b)}
a,b=rows(old),rows(new);assert len(a)==len(b)==149
changes=load(P/'handoff.json')['candidates'];assert set(n for n in a if a[n]!=b[n])==set(changes)
expected=old
for n,v in changes.items():
 assert a[n].count(('numerators: '+v['old']).encode())==1
 replacement=a[n].replace(('numerators: '+v['old']).encode(),('numerators: '+v['new']).encode());assert replacement==b[n]
 expected=expected.replace(a[n],replacement)
assert expected==new
prior=['xs_qed_eemm_L1','xs_qed_eedd_L1','xs_qed_eeee_1pi_L1','qed_eemm_L1','ctorder_qcd_uugg_L1','sm_eemm_L0','ym_g6_L0','xs_qed_eemm_L2','sm_eeww_L0','qed_ee4f_L0','xs_qed_eeee_L1','xs_qcd_gguu_1pi_L1','xs_qed_eemm_identical_L2']
assert all(a[n]==b[n] for n in prior);assert a['xs_qcd_gggg_sym_L1']==b['xs_qcd_gggg_sym_L1']
sc=[]
for section in tomllib.loads((W/'docs/developers.toml').read_text())['section']:
 for note in section.get('note',[]):
  for scope in note.get('scope',[]):
   t=(W/scope['path']).read_text();start,end=scope.get('anchor'),scope.get('end_anchor');assert not start or t.count(start)==1;assert not end or t.count(end)==1
   data=t[t.index(start):t.index(end)].encode() if end else (W/scope['path']).read_bytes();assert blob(data)==scope['digest'];sc.append({'note':note['id'],**scope})
assert len(sc)==118
assert re.search(r'^  enable = false;$',(W/'nix-ci.nix').read_text(),re.M)
write('input-verification.json',{'pins':{str(p):h for p,h in pins.items()},'manifests':manifest_records,'external_files':len(external),'unique_verified_files':checked,'no_private_config_read':True})
write('source-invariance.json',{'base':BASE,'start':START,'base_tree':base,'exact_before_blob':blob(old),'exact_after_blob':blob(new),'changes':changes,'all149_Rows':True,'unchanged_full_Rows':146,'prior13':prior,'prior13_rows_sha256':{n:hashlib.sha256(a[n]).hexdigest() for n in prior},'same_hash_sibling_unchanged':True,'scope_bindings':sc,'registry_before':{'blob':base['docs/developers.toml']['blob'],'sha256':sha(W/'docs/developers.toml')},'stub_before':base['examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi'],'all_other_source_bytes_equal':True})
main=W.parent/'main';index=main/'.git/index'
write('start.json',{'utc':datetime.datetime.now(datetime.UTC).isoformat(),'JJ_refs':jj('bookmark','list','--all'),'Git_refs':git('show-ref').decode(),'main_head':git('rev-parse','HEAD').decode().strip(),'main_index_sha256':sha(index),'main_status':subprocess.check_output(['git','-C',str(main),'status','--porcelain=v1','--untracked-files=no'],env=ENV,text=True),'main_tracked_files':tree('HEAD'),'author':'Lucien Huber <im@lcnbr.ch>','tracking_cap':'1MiB; unchanged','resource_policy':'existing two-heavy-slot /16worker /240s; generator one-job two-core schedule; core union untouched'})
print('PASS pins, manifests,',len(checked),'unique files,149 Rows/13 prior/118 scopes; exact authorized3-token bytes; no tracked edit')
