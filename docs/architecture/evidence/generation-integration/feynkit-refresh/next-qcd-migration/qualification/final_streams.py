"""Bind closed finalizer streams. No active stream is listed as a completed input."""
import ast,gzip,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration'
def load(p):return json.loads(Path(p).read_text())
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings=load(E/'qualification-copy-bindings.json');ns={'O':O,'W':W,'E':E,'gzip':gzip,'hashlib':hashlib,'json':json,'Path':Path,'bindings':bindings,'sha':lambda b:hashlib.sha256(b).hexdigest()};t=ast.parse((O/'archive_inputs.py').read_bytes());fn=next(x for x in t.body if isinstance(x,ast.FunctionDef) and x.name=='cp');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(O/'archive_inputs.py'),'exec'),ns)
for name in ['close_packet.stdout','close_packet.stderr','final_streams.py']:ns['cp'](O/name,Path('qualification')/name)
save(E/'qualification-copy-bindings.json',bindings);r=load(E/'receipt.json');r['retention']['committed_qualification_records']=len(bindings);r['assembly_stream_closure']['finalizer_streams']='Closed complete finalizer stdout/stderr are committed in qualification/close_packet.{stdout,stderr}; last sealing helper is silent';save(E/'receipt.json',r)
files=sorted([p for p in E.rglob('*') if p.is_file() and p not in [E/'SHA256SUMS',E/'SEAL.json']]+[W/'docs/developers.toml',W/'docs/architecture/generation-next-qcd-migration.typ',W/'crates/feynkit-generator/tests/generation_oracle/main.rs']);assert all(p.stat().st_size<=1048576 for p in files)
(E/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(W))+'\n' for p in files));s=load(E/'SEAL.json');s.update(manifest_sha256=sha(E/'SHA256SUMS'),receipt_sha256=sha(E/'receipt.json'),files=len(files),bytes=sum(p.stat().st_size for p in files));save(E/'SEAL.json',s)
