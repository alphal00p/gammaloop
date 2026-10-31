"""Bounded pure producer qualification. Every subprocess has all license keys removed."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
O = Path(__file__).resolve().parent
OLD = W / '.local/feynkit-refresh/qualification'
REV = '432b6544f3d2d884f6efffb8be104ae5820183b0'
URI = 'git+file:///shared/lcnbr/dev/gammaloop/main?rev=' + REV + '&shallow=1&allRefs=1'
KEYS = ('SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED', 'SYMBOLICA_LICENSE_FILE')
ALLOWED = ('PATH','HOME','USER','LOGNAME','LANG','LC_ALL','TMPDIR','NIX_REMOTE','NIX_BUILD_CORES','SSL_CERT_FILE','NIX_SSL_CERT_FILE')
ENV = {k:os.environ[k] for k in ALLOWED if k in os.environ}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')
assert not any(k in ENV for k in KEYS)

def write(name, data):
    (O/name).write_text(json.dumps(data, indent=2) + '\n')

def load(name):
    return json.loads((O/name).read_text())

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def call(name, args, parse=True):
    assert not any(k in ENV for k in KEYS)
    r={'argv':args, 'cwd':str(W), 'license_keys_removed':list(KEYS), 'start':datetime.datetime.now(datetime.UTC).isoformat()}
    write(name+'.command.json',r)
    with (O/(name+'.stdout')).open('wb') as out, (O/(name+'.stderr')).open('wb') as err:
        p=subprocess.run(args, cwd=str(W), env=ENV, stdout=out, stderr=err)
    r.update(exit=p.returncode, end=datetime.datetime.now(datetime.UTC).isoformat(),stdout_sha256=sha(O/(name+'.stdout')),stderr_sha256=sha(O/(name+'.stderr')))
    write(name+'.command.json',r)
    assert p.returncode==0, (name,p.returncode)
    return json.loads((O/(name+'.stdout')).read_text()) if parse else (O/(name+'.stdout')).read_bytes()

def nix(name,args):
    r=call(name,['nix',*(['--offline'] if args[0]!='build' else []),'--option','allow-import-from-derivation','false',*args])
    if args[:2]==['derivation','show']:
        r=r.get('derivations',r)
        r={k if k.startswith('/') else '/nix/store/'+k:v for k,v in r.items()}
        assert all(k not in d.get('env',{}) for d in r.values() for k in KEYS)
    return r

def evaluate():
    assert re.fullmatch(r'git\+file:///shared/lcnbr/dev/gammaloop/main\?rev=[0-9a-f]{40}&shallow=1&allRefs=1',URI)
    m=nix('evaluate',['eval','--json','--expr','let f=builtins.getFlake '+json.dumps(URI)+'; in {source=f.outPath; drvPath=f.checks.x86_64-linux.gammaloop-nextest-binaries-feynkit.drvPath;}'])
    write('flake.json',dict(uri=URI,rev=REV,**m))
    d=next(iter(nix('group-derivation',['derivation','show',m['drvPath']]).values()))
    ids=list(d.get('inputDrvs',d.get('inputs',{}).get('drvs',{})))
    ids=[p if p.startswith('/') else '/nix/store/'+p for p in ids]
    ds=nix('direct-derivations',['derivation','show',*ids])
    chosen={p:d for p,d in ds.items() if d.get('env',{}).get('pname')=='gammaloop-nextest-binaries-feynkit-feynkit-generator'}
    assert len(chosen)==1
    write('archive-derivation.json',chosen)
    (O/'selected-drv.txt').write_text(next(iter(chosen))+'\n')
    old=json.loads((OLD/'archive-derivation.json').read_text()); od=next(iter(old.values())); nd=next(iter(chosen.values()))
    invariant={k:od['env'].get(k)==nd['env'].get(k) for k in ['RUSTFLAGS','cargoExtraArgs','buildPhase','installPhase','nativeBuildInputs','buildInputs','stdenv']}
    write('producer-delta.json',{'old_drv':next(iter(old)),'new_drv':next(iter(chosen)),'old_source':od['env']['src'],'new_source':nd['env']['src'],'old_artifacts':od['env']['cargoArtifacts'],'new_artifacts':nd['env']['cargoArtifacts'],'invariant_env_fields':invariant})
    print('selected new archive',next(iter(chosen)),flush=True)

def build():
    drv=(O/'selected-drv.txt').read_text().strip()
    assert re.fullmatch(r'/nix/store/[a-z0-9]{32}-gammaloop-nextest-binaries-feynkit-feynkit-generator-0\.1\.0\.drv',drv)
    assert not (O/'build.command.json').exists(), 'one new producer build only'
    nix('build',['build','--max-jobs','1','--cores','2','--no-link','--json',drv+'^out'])
    print('selected archive realization PASS',flush=True)

if __name__=='__main__':
    {'evaluate':evaluate,'build':build}[sys.argv[1]]()
