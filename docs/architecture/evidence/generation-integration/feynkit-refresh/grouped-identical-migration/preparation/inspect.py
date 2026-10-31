from pathlib import Path
import hashlib
import json
import os
import sys

W = Path('/shared/lcnbr/dev/gammaloop/generation-oracle-numerator-families')
O = W / 'target/grouped-identical-migration-preparation-20261031-r1'
N = W / 'docs/architecture/evidence/generation-xs-qed-eemm-identical-aggregation-r1'
A = W / 'docs/architecture/evidence/generation-xs-qed-eemm-identical-auto-aggregation-r1'
def audit(event, args):
    if event == 'open':
        p, mode, flags = args
        if isinstance(p, (str, bytes)):
            path = Path(os.fsdecode(p)).resolve()
            if flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND):
                assert path.is_relative_to(O), ('write outside output', str(path))
            assert not str(path).startswith('/home/lcnbr/.config'), 'private config read'
    if event.startswith('socket.') or event in ('subprocess.Popen', 'os.system'):
        raise RuntimeError(('forbidden process/network', event))
sys.addaudithook(audit)
assert not any('LICENSE' in k for k in os.environ)
def read(p):
    p = Path(p)
    b = p.read_bytes()
    with (O / 'reads.jsonl').open('a') as f:
        f.write(json.dumps({'path':str(p), 'sha256':hashlib.sha256(b).hexdigest(), 'bytes':len(b)})+'\n')
    return b
def j(p):
    return json.loads(read(p))
def save(name, x):
    (O / name).write_text(json.dumps(x, indent=2)+'\n')
if __name__ == '__main__':
    for s in sys.argv[1:]:
        p = W / s
        if p.is_dir():
            print(s, sorted(x.name for x in p.iterdir()))
        else:
            x = j(p)
            if isinstance(x, dict):
                print(s, json.dumps({k:(v if not isinstance(v,(dict,list)) else {'kind':type(v).__name__,'len':len(v),'keys':list(v)[:12] if isinstance(v,dict) else None}) for k,v in x.items()},indent=2))
            else:
                print(s,'list',len(x),json.dumps(x[:1])[:5000])
