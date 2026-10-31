"""Retain the exact accepted plan bytes rather than relying on a mutable shared path."""
import ast,gzip,hashlib,json
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/next-qcd-migration';P=W.parent/'generation-integration-plan/docs/architecture/generation-integration-plan.typ'
assert hashlib.sha256(P.read_bytes()).hexdigest()=='4489dc7383468cf2daf7f8c177d3ea87fe9cad45c99af1afc0f160b202ff06f6'
bindings=json.loads((O/'input-copy-bindings.json').read_text());ns={'O':O,'W':W,'E':E,'gzip':gzip,'hashlib':hashlib,'json':json,'Path':Path,'bindings':bindings,'sha':lambda b:hashlib.sha256(b).hexdigest()};t=ast.parse((O/'archive_inputs.py').read_bytes());fn=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='cp');exec(compile(ast.Module(body=[fn],type_ignores=[]),str(O/'archive_inputs.py'),'exec'),ns);ns['cp'](P,'roots/accepted-integration-plan.typ');(O/'input-copy-bindings.json').write_text(json.dumps(bindings,indent=2)+'\n')
print('Preserved exact accepted plan as input record',len(bindings))
