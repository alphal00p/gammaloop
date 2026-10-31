"""Close assembler streams before the primary immutable freeze; preserve provisional seals."""
import ast,gzip,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration'
def load(p):return json.loads(Path(p).read_text())
def save(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name in ['SHA256SUMS','SEAL.json']:(O/('assembly-provisional-'+name)).write_bytes((E/name).read_bytes())
bindings=load(E/'qualification-copy-bindings.json');ns={'O':O,'W':W,'E':E,'gzip':gzip,'hashlib':hashlib,'json':json,'Path':Path,'bindings':bindings,'sha':lambda b:hashlib.sha256(b).hexdigest()}
t=ast.parse((O/'archive_inputs.py').read_bytes());fn=next(x for x in t.body if isinstance(x,ast.FunctionDef) and x.name=='cp');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(O/'archive_inputs.py'),'exec'),ns)
for name in ['assemble_packet.stdout','assemble_packet.stderr','close_packet.py','assembly-provisional-SHA256SUMS','assembly-provisional-SEAL.json']:
 p=O/name;bindings[:]=[r for r in bindings if r['original']!=str(p)];ns['cp'](p,Path('qualification')/name)
save(E/'qualification-copy-bindings.json',bindings)
r=load(E/'receipt.json');r['retention']['committed_qualification_records']=len(bindings);r['assembly_stream_closure']={'status':'Closed complete assembler streams rebound before immutable freeze','provisional_seal':'qualification/assembly-provisional-SEAL.json','provisional_manifest':'qualification/assembly-provisional-SHA256SUMS','provisional_scope':'In-progress assembly snapshot only; not used as complete closure or an immutable candidate','finalizer_streams':'Closed finalizer diagnostics are preserved in the subordinate durability packet before its freeze'};save(E/'receipt.json',r)
for p in O.glob('*.py'):ast.parse(p.read_bytes())
for r in bindings:
 pieces=[]
 for part in r['copies']:
  p=W/part['path'];assert sha(p)==part['sha256'] and p.stat().st_size==part['bytes'];pieces.append(p.read_bytes())
 b=b''.join(pieces);b=b if r['encoding']=='identity' else gzip.decompress(b);assert hashlib.sha256(b).hexdigest()==r['sha256'] and len(b)==r['bytes']
files=sorted([p for p in E.rglob('*') if p.is_file() and p not in [E/'SHA256SUMS',E/'SEAL.json']]+[W/'docs/developers.toml',W/'docs/architecture/generation-next-qcd-migration.typ',W/'crates/feynkit-generator/tests/generation_oracle/main.rs']);assert all(p.stat().st_size<=1048576 for p in files)
(E/'SHA256SUMS').write_text(''.join(sha(p)+'  '+str(p.relative_to(W))+'\n' for p in files));s=load(E/'SEAL.json');s.update(manifest_sha256=sha(E/'SHA256SUMS'),receipt_sha256=sha(E/'receipt.json'),files=len(files),bytes=sum(p.stat().st_size for p in files));save(E/'SEAL.json',s)
print('PASS closed primary evidence',len(files),'entries,',len(bindings),'exact lossless qualification bindings; seal',sha(E/'SEAL.json'))
