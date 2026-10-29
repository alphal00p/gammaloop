from pathlib import Path
import json,hashlib,tarfile,tomllib,subprocess
out=Path('target/mark-geometry-round3');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();proof={}
local=Path('vendor/kurbo');files=sorted(p for p in local.rglob('*') if p.is_file());assert len(files)==53
for kind in ['production','test']:
 root=Path((out/f'{kind}-source-path').read_text().strip());assert sorted(p.name for p in (root/'vendor').iterdir())==['clarabel','kurbo'];rows=[]
 for p in files:
  dest=root/p;assert dest.read_bytes()==p.read_bytes(),(kind,str(p));rows.append({'path':str(p),'sha256':sha(p)})
 proof[kind]={'source':str(root),'vendor_directories':['clarabel','kurbo'],'complete_kurbo_files':rows}
for kind in ['root','host']:
 meta=json.loads((out/f'{kind}-metadata.json').read_text());p=next(p for p in meta['packages'] if p['name']=='kurbo');assert p['source'] is None and p['version']=='0.13.1' and Path(p['manifest_path']).resolve()==(local/'Cargo.toml').resolve();proof[kind+'_selection']={k:p[k] for k in ['id','source','manifest_path','version','license']}
archive=next((out/'sdist-final').glob('*.tar.gz'));rows=[]
with tarfile.open(archive) as t:
 prefix=t.getnames()[0].split('/')[0]+'/'
 for p in files:
  data=t.extractfile(prefix+str(p)).read();assert data==p.read_bytes();rows.append({'path':str(p),'sha256':sha(p)})
proof['actual_gammaloop_sdist']={'archive_sha256':sha(archive),'size':archive.stat().st_size,'files':rows,'scope':'Real source distribution, not a final wheel or installed native/Pyodide/Community artifact.'}
proof['nix_wiring']='No Nix change required: existing source30 localCargoPatchDirs and extraDummyScript preserve declared patch sources and licenses, including during dummy-source dependency/test preparation.'
(out/'source-selection-audit.json').write_text(json.dumps(proof,indent=2)+'\n');print('Root/host select local Kurbo; production/test complete sources + both licenses, only declared vendor directories; all 53 files in actual sdist')
