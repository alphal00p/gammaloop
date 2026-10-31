import hashlib,json,sys,os
from pathlib import Path
import symbolica.core as core
from symbolica.community import hepkit,feynkit_native
e=json.loads(Path('/proposal/expected.json').read_text());p=Path(core.__file__);assert str(p).startswith('/work/venv/');assert hashlib.sha256(p.read_bytes()).hexdigest()==e['native_sha256'];assert hepkit.FeynmanDiagram is feynkit_native.FeynmanDiagram;stub=Path(hepkit.__file__).with_suffix('.pyi');assert hashlib.sha256(stub.read_bytes()).hexdigest()==e['stub_sha256'];assert hasattr(hepkit,'IBPFamily')==e['private'];assert __debug__
Path('/work/identity.json').write_text(json.dumps({'core_path':str(p),'core_sha256':e['native_sha256'],'stub_path':str(stub),'stub_sha256':e['stub_sha256'],'IBPFamily':hasattr(hepkit,'IBPFamily'),'native_class_identity':True,'assertions_enabled':__debug__},indent=2)+'\n');print('Exact installed wheel identity passed')
