"""Read-only evidence access; all authoring is bounded to this private packet."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

O = Path(__file__).resolve().parent
W = O.parents[1]
T = W / 'target'
CAP = T / 'raw-next-qcd-captures-20261031-r2'
PRE = T / 'raw-next-qcd-preflight-20261031-r1'
R1 = T / 'raw-next-qcd-proof-preparation-20261031-r1'
R2 = T / 'raw-next-qcd-proof-preparation-20261031-r2'
R3 = T / 'raw-next-qcd-proof-preparation-20261031-r3'
N = T / 'raw-next-qcd-proof-native-execution-20261031-r1'
REFRESH = W.parent / 'generation-feynkit-refresh'
TARGET = '83276dc047815d1dcc7f003aeb512c204b26758e'
CODE = '6f41807b064bcf2f579883134daef65c8ac279d3'
ACTUAL = 'ffa76ec420642c03349ab756d2d7720ab93c0057'
HISTORICAL = '2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94'
MAIN = 'crates/feynkit-generator/tests/generation_oracle/main.rs'
CASES = ('xs_qcd_gguu_sym_L1', 'xs_qcd_gguu_L1', 'xs_qcd_gggg_L1')
MODES = ('None', 'Auto')
used = {}

class Reject(Exception):
    def __init__(self, predicate, detail=''):
        self.predicate = predicate
        super().__init__(predicate + ': ' + str(detail))

def need(ok, predicate, detail=''):
    if not ok:
        raise Reject(predicate, detail)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()

def restricted(path):
    s = str(path)
    return '/.config/' in s or '/.cargo/config' in s or s.endswith('/.envrc') or s.endswith('/.env') or s.endswith('/license')

def read(path):
    p = Path(path).absolute()
    need(not restricted(p), 'no-config-read', p)
    b = p.read_bytes()
    used[str(p)] = {'sha256': sha(b), 'bytes': len(b)}
    return b

def unique(pairs):
    d = {}
    for k,v in pairs:
        need(k not in d, 'unique-json-key', k)
        d[k] = v
    return d

def load(path):
    return json.loads(read(path), object_pairs_hook=unique, parse_constant=lambda x: need(False, 'finite-json', x))

def put(name, data):
    p = (O / name).resolve()
    need(p.is_relative_to(O), 'private-write-boundary', p)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.exists():
        need(p.read_bytes()==data, 'existing-private-output-byte-equal', p)
        return
    with p.open('xb') as f:
        f.write(data)

def save(name, data):
    put(name, (json.dumps(data, indent=2, ensure_ascii=False) + '\n').encode())

def ref(path):
    p = Path(path).absolute()
    read(p)
    return {'path': str(p), **used[str(p)]}

def git(*args):
    need(args[0] in ('cat-file', 'ls-tree', 'diff', 'show', 'rev-parse', 'for-each-ref', 'check-ignore'), 'read-only-git', args)
    env = dict(os.environ, GIT_OPTIONAL_LOCKS='0')
    return subprocess.check_output(['git', '--git-dir=/shared/lcnbr/dev/gammaloop/main/.git', *args], env=env, cwd=O)

def jj_state():
    return subprocess.check_output(['jj', '--ignore-working-copy', 'log', '--no-pager', '--no-graph', '-r', 'generation-integration | @', '-T', 'commit_id ++ " " ++ change_id ++ "\\n"'], cwd=W).decode()

def fingerprint(path):
    p = Path(path)
    need(not restricted(p), 'no-config-read', p)
    h = hashlib.sha256()
    count = 0
    with p.open('rb') as f:
        while data := f.read(1024*1024):
            h.update(data)
            count += len(data)
    return {'sha256': h.hexdigest(), 'bytes': count}

def manifest(path):
    p = Path(path)
    rows = []
    for line in read(p).decode().splitlines():
        h, name = line.split('  ', 1)
        q = Path(name)
        rows.append((q if q.is_absolute() else p.parent/q, h))
    return rows
