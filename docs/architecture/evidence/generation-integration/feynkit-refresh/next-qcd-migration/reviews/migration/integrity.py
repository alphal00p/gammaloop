import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import common as C
mode=sys.argv[1]
pins={C.P/'SEAL.json':'b4e75df006b68563219947ecbdef5cb8fd289f9b4c7ba840d79c0ce599377324',C.P/'SHA256SUMS':'7c8c48341828a16d00b4f692f6ce3941d35b5679e73eac2fd4176f139942b503',C.P/'handoff.json':'1b69607db0184093de3b16d32541d6f7a5a3f00968b85694685a6d17fb1506b3',C.P/'proposed-not-applied.patch':'1df84cbea910df5ff092203ff6050e785bccf0027bd627b91a601d1f7b98cf84',C.P/'complete-input-manifest.json':'600ced69f0172b82a4043aa3d2a255e270c9782a79726f4d745ab93d3cc40e6e',C.P/'before.rs':'66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2',C.P/'after.NOT-APPLIED.rs':'17e0abe46c76c6b6922e5b33782925020e1412e1a6ed8cd8dc7e23de438d4aa2',C.S/'review.json':'6d892bd9263a7c093bba5fd372fb4663a5a4e81dfaad4bfc709bbba052b96d21',C.S/'SEAL.json':'888b82c478697e917cb243587653cc16d09a484f299a671c5ffd8981e58c1d99',C.F/'root-next-qcd-scientific-interface-acceptance-20261031.json':'225b2c244ba02f16ca3ee73a59d23df3fff12204deef44f3aa405397247a0ea7',C.F/'root-next-qcd-migration-preparation-receipt-20261031.json':'4c66d289489823726da2bf9355c43a937830a970e0e1e53e18f34e4a530f8804'}
expected={str(p):h for p,h in pins.items()};counts={}
def add(p,h):
 p=str(p);assert p not in expected or expected[p]==h,(p,expected.get(p),h);expected[p]=h
for folder,name in [(C.P,'SHA256SUMS'),(C.S,'SHA256SUMS'),(C.CAP,'SHA256SUMS'),(C.R1,'SHA256SUMS'),(C.R2,'SHA256SUMS'),(C.R3,'SHA256SUMS'),(C.N,'COMPLETE-SHA256SUMS'),(C.N,'EXECUTION-SHA256SUMS')]:
 m=folder/name;rows=C.manifest(m);counts[str(m)]=len(rows);add(m,C.sha(C.read(m)))
 for p,h in rows:add(p,h)
x=C.load(C.P/'complete-input-manifest.json')
for p,r in x['external_files'].items():add(p,r['sha256'])
for p,r in x['source_config_references_not_read'].items():add(p,r['recorded_sha256'])
y=C.load(C.S/'integrity-final.json')
for r in y['checks']:
 # Exact schema checked below in initial run.
 p=r['path'];h=r.get('expected',r.get('expected_sha256',r.get('sha256')))
 if isinstance(h,dict):h=h['sha256']
 assert isinstance(h,str) and len(h)==64,r
 add(p,h)
for p in [C.W/'AGENTS.md',C.W/'CONTRIBUTING.typ',C.W/'.gitignore',C.W.parent/'main/AGENTS.md',C.W.parent/'main/CONTRIBUTING.typ',C.W.parent/'generation-integration-plan/docs/architecture/generation-integration-plan.typ']:
 add(p,C.sha(C.read(p)))
if mode=='final':
 for f in C.O.glob('*-consumed.json'):
  for p,r in C.load(f).items():
   if not Path(p).is_relative_to(C.O):add(p,r['sha256'])
results={};failures=[]
for p,h in sorted(expected.items()):
 try:
  r=C.fingerprint(p);results[p]=r
  if r['sha256']!=h:failures.append({'path':p,'expected':h,'actual':r['sha256']})
 except Exception as e:failures.append({'path':p,'error':repr(e)})
report={'phase':mode,'unique_files':len(results),'manifest_entries':counts,'pinned_files':{str(p):h for p,h in pins.items()},'checks':results,'failures':failures,'source_configs_freshly_rehashed':list(x['source_config_references_not_read']),'refs':C.refs()}
C.save('integrity-'+mode+'.json',report)
print({k:v for k,v in report.items() if k not in ('checks','refs','pinned_files','manifest_entries')})
assert not failures,failures
if mode=='final':
 old=C.load(C.O/'integrity-before.json');assert all(report['checks'][p]==r for p,r in old['checks'].items())
 C.save('refs-final-disclosure.json',{'before':old['refs'],'after':report['refs'],'equal':old['refs']==report['refs'],'action':'Read-only disclosure; no restoration or mutation.'})
