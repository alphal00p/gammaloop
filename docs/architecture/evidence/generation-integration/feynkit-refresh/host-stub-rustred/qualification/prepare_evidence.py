from pathlib import Path
import hashlib,json,os,shutil,tomllib,subprocess
P=Path(__file__).resolve().parent;W=P.parents[2];D=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/host-stub-rustred'
D.mkdir(parents=True,exist_ok=False)
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def copy(src,dst):dst=D/dst;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
R=W/'.local/review-host-browser-stub-astra-20261031-013106Z-5655cc61';T=W/'.local/review-raw-three-integration-astra-20261031-012811Z-323f4486'
for name,src in [('host-browser-stub',R),('bounded-twelve-row',T)]:
 for f in ('review.json','SEAL.json','SHA256SUMS'):copy(src/f,Path('reviews')/name/f)
for n in ('host-api-stub-correction-20261031-r1','host-browser-correction-20261031-r1','rustred-host-preflight-20261031-r1'):
 src=P.parent/n
 for f in ('report.json','report.typ','SHA256SUMS','SEAL.json','seal.json'):
  if (src/f).exists():copy(src/f,Path('premises')/n/f)
for src,dst in [(P.parent/'host-api-stub-correction-20261031-r1/candidate.patch','premises/approved-stub.patch'),(P.parent/'rustred-host-preflight-20261031-r1/proposal/private-host.patch','premises/original-private-host.patch'),(P.parent/'rustred-host-preflight-20261031-r1/proposal/private-producer.nix','premises/original-private-producer.nix'),(P.parent/'host-browser-correction-20261031-r1/egl-runtime.json','premises/egl-runtime.json'),(P.parent/'host-browser-correction-20261031-r1/extra-runtime-data-hashes.json','premises/extra-runtime-data-hashes.json'),(P.parent/'host-api-stub-correction-20261031-r1/inspection/formatter-cache-correction.json','premises/formatter-cache-correction.json'),(R/'inspection/combined-static-symbol-difference.json','premises/combined-static-symbol-difference.json')]:copy(src,dst)
(P/'evidence-directory.json').write_text(json.dumps({'path':str(D),'purpose':'Additive exact prior review and seal bindings plus small new qualification evidence; no independent review of this candidate yet.'},indent=2)+'\n')
print(D)
