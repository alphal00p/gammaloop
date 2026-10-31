from pathlib import Path
import hashlib,json,os,subprocess,datetime
O=Path(__file__).resolve().parent
W=O.parents[1]
P=W/'target/raw-next-qcd-migration-preparation-20261031-r1'
S=W/'.local/review-next-qcd-scientific-interface-astra-20261031-5a3KbP5x'
T=W/'target'; CAP=T/'raw-next-qcd-captures-20261031-r2'
R1=T/'raw-next-qcd-proof-preparation-20261031-r1'; R2=T/'raw-next-qcd-proof-preparation-20261031-r2'; R3=T/'raw-next-qcd-proof-preparation-20261031-r3'
N=T/'raw-next-qcd-proof-native-execution-20261031-r1'
F=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh'
REV='83276dc047815d1dcc7f003aeb512c204b26758e'; ACT='ffa76ec420642c03349ab756d2d7720ab93c0057'; CODE='6f41807b064bcf2f579883134daef65c8ac279d3'
MAIN='crates/feynkit-generator/tests/generation_oracle/main.rs'
CASES=('xs_qcd_gguu_sym_L1','xs_qcd_gguu_L1','xs_qcd_gggg_L1')
MODES=('None','Auto')
used={}
def sha(b): return hashlib.sha256(b).hexdigest()
def blob(b): return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def safe(p):
 p=Path(p).absolute()
 assert str(p).startswith(('/shared/lcnbr/dev/gammaloop/','/nix/store/')),str(p)
 assert p.name not in ('.envrc','.env','license'),str(p)
 return p
def read(p):
 p=safe(p);b=p.read_bytes();used[str(p)]={'sha256':sha(b),'bytes':len(b)};return b
def unique(pairs):
 d={}
 for k,v in pairs:
  assert k not in d,k
  d[k]=v
 return d
def load(p): return json.loads(read(p),object_pairs_hook=unique,parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
def put(name,b):
 p=(O/name).resolve();assert p.is_relative_to(O)
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
def save(name,v): put(name,(json.dumps(v,indent=2,ensure_ascii=False)+'\n').encode())
def flush(name): save(name,used)
def git(*args):
 assert args[0] in ('show','ls-tree','cat-file','rev-parse','for-each-ref','diff')
 return subprocess.check_output(['git','--git-dir=/shared/lcnbr/dev/gammaloop/main/.git',*args],cwd=W,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0','GIT_CONFIG_NOSYSTEM':'1','GIT_CONFIG_GLOBAL':'/dev/null'})
def refs(): return git('for-each-ref','--format=%(refname) %(objectname)').decode()
def fingerprint(p):
 p=safe(p);h=hashlib.sha256();n=0
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b);n+=len(b)
 return {'sha256':h.hexdigest(),'bytes':n}
def manifest(p):
 return [(Path(n) if Path(n).is_absolute() else p.parent/n,h) for h,n in (line.split('  ',1) for line in read(p).decode().splitlines())]
