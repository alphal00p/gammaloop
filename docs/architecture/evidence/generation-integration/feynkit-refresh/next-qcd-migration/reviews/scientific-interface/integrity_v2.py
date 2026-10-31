import hashlib,json,pathlib,sys,datetime
W=pathlib.Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families'); O=pathlib.Path(__file__).resolve().parent
E=W/'target/raw-next-qcd-proof-native-execution-20261031-r1'; R=W/'target/raw-next-qcd-proof-preparation-20261031-r3'
checks=[];cache={}
def check(p,h,kind,n=None):
 p=pathlib.Path(p); assert not any(x in str(p) for x in ['/home/lcnbr/.config/','/license','/.netrc']),str(p)
 if str(p) not in cache:
  d=hashlib.sha256();size=0
  with p.open('rb') as f:
   while b:=f.read(1024*1024):d.update(b);size+=len(b)
  cache[str(p)]=(d.hexdigest(),size)
 v,s=cache[str(p)];checks.append(dict(path=str(p),kind=kind,expected=h,actual=v,bytes=s,ok=v==h and (n is None or n==s)))
def manifest(p,kind):
 for line in p.read_text().splitlines():
  h,name=line.split('  ',1);check(p.parent/name,h,kind)
pins={E/'report.json':'3baf87d1d031f790b753bb944d2a4c180f0620297b7808ad348dfa7ea65c76f0', E/'COMPLETE-SEAL.json':'22b11020019c4b9f8e058a2951de177f5a6f111fbb4f6976932f98ef0c399626', E/'COMPLETE-SHA256SUMS':'844465ea17e9142778839984035d1851b73866aa371a0aee75a278f9f7d5fc7a', E/'EXECUTION-SHA256SUMS':'a08f0940cc6e282d98540ac0dc10779fa747a0c5c37e438b53fd87b5d6638fc5', E/'SEAL.json':'522817b3b339b7ce0d8b82cdb013c03721da8b99700ee98f4be63ff3d309d0ef', W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/root-next-qcd-proof-native-execution-assignment-20261031-r3.json':'650f597bb727c64e53898d676f2884ee2125a5e73be8a8c0f0edf28376ba1382',R/'SEAL.json':'a7e3a1e8e4e8e6d979d8cf92f8a9de5878e88d6ac1f105bf19a86611e4381cc1', W/'.local/review-next-qcd-proof-preparation-astra-20261031-6xpiqaAK/review.json':'62dade9062105a077ae6089d7e261742d57d63d9d07737b020c9b3098a71db95', W/'.local/review-next-qcd-proof-r2-supervisor-astra-20261031-mtFpBwVu/review.json':'64e9bf551fc949878e0c6647238fad4daf9eb69d82e0622939c89b138ca8d1d3', W/'.local/review-next-qcd-supervisor-r3-astra-20261031-22d61lkb/review.json':'6ad6242a78d52d0b371fb39d3942e66ed32bbbafec8878a91e6a72d5e3a9208b'}
for p,h in pins.items():check(p,h,'user pin')
for f in ['COMPLETE-SHA256SUMS','EXECUTION-SHA256SUMS','SHA256SUMS']:manifest(E/f,f)
for d in [R,W/'target/raw-next-qcd-proof-preparation-20261031-r2',W/'target/raw-next-qcd-proof-preparation-20261031-r1',W/'target/raw-next-qcd-captures-20261031-r2',W/'.local/review-next-qcd-proof-preparation-astra-20261031-6xpiqaAK',W/'.local/review-next-qcd-proof-r2-supervisor-astra-20261031-mtFpBwVu',W/'.local/review-next-qcd-supervisor-r3-astra-20261031-22d61lkb']:
 seal=json.loads((d/'SEAL.json').read_text());check(d/'SHA256SUMS',seal['manifest_sha256'],'historical seal manifest');manifest(d/'SHA256SUMS',d.name)
before=json.loads((E/'inputs-before.json').read_text());after=json.loads((E/'inputs-after.json').read_text());assert before==after
for p,m in before.items():check(p,m['sha256'],'before/after input',m['bytes'])
for p,m in json.loads((E/'execution-artifacts.json').read_text()).items():assert p==m['path'];check(p,m['sha256'],'execution artifact',m['bytes'])
report={'timestamp':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'failures':[x for x in checks if not x['ok']],'unique_files':len(cache),'input_count':len(before)}
(O/('integrity-'+(sys.argv[1] if len(sys.argv)>1 else 'before')+'.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}));assert not report['failures']
