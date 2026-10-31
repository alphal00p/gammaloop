import hashlib
import json
from pathlib import Path
import runpy

O=Path(__file__).resolve().parent
U=runpy.run_path(str(O/'inspect.py'))
W,A,j,save=(U[k] for k in ['W','A','j','save'])
B=runpy.run_path(str(O/'bind_evidence.py'))
sha=B['sha']
expected=j(A/'protected-inputs.json')['files']
expected.update({p:{'sha256':h} for p,h in j(O/'protected-inputs-before.json').items()})
D=W/'.local/review-grouped-auto-aggregation-astra-20261031-independent-r1-1793413266607828149'
seal=j(D/'seal.json')
for p,r in seal['files'].items():expected[str(D/p)]=r
cli=Path(j(O/'consumed-Auto-independent-review-copy.json')['fresh_full_protected_CLI']['output'])
for p,r in seal['fresh_protected_CLI_output_files'].items():expected[str(cli/p)]=r
root=W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/root-grouped-auto-acceptance-20261031.json'
for p in [D/'review.json',D/'seal.json',root]:expected[str(p)]={'sha256':sha(p)}
observed={}
failures=[]
for path,r in expected.items():
    p=Path(path)
    try:
        h=sha(p)
        observed[path]=h
        if h!=r['sha256'] or ('bytes' in r and p.stat().st_size!=r['bytes']):failures.append({'path':path,'expected':r,'actual_sha256':h})
    except OSError as e:failures.append({'path':path,'error':str(e)})
save('complete-input-guard.json',{'expected_files':len(expected),'verified_files':len(observed),'failures':failures,'status':'PASS' if not failures else 'FAIL_EXTERNAL_INPUT_MOVEMENT','scope':'Complete accepted Auto protected file closure plus original producer/None candidate seals and newly supplied Auto review/CLI. Inert hashing only; no helper execution.'})
save('complete-input-hashes.json',observed)
print(json.dumps({'expected':len(expected),'failures':failures[:5],'status':'PASS' if not failures else 'FAIL'}))
if failures:raise SystemExit(1)
