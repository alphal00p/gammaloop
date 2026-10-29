set -euo pipefail
source target/mark-geometry-round3/environment.sh
out=target/mark-geometry-round3
python3 - <<'PY'
from pathlib import Path
import json,hashlib
out=Path('target/mark-geometry-round3');old=(out/'clean-result').resolve();new=(out/'final-result').resolve();oldfiles=sorted(str(p.relative_to(old)) for p in old.rglob('*') if p.is_file());newfiles=sorted(str(p.relative_to(new)) for p in new.rglob('*') if p.is_file());assert oldfiles==newfiles;rows=[]
for rel in oldfiles:
 assert (old/rel).read_bytes()==(new/rel).read_bytes();rows.append({'path':rel,'sha256':hashlib.sha256((new/rel).read_bytes()).hexdigest()})
(out/'final-payload-bindings.json').write_text(json.dumps({'runtime_asset_tested_output':str(old),'final_source_qualified_output':str(new),'files':rows,'all29_byte_identical':len(rows)==29},indent=2)+'\n')
p=Path('result');assert p.is_symlink();p.unlink();p.symlink_to(new)
p=Path('.local/mark-geometry-round3-producers/wasm-bundle')
if p.is_symlink():p.unlink()
paths=json.loads((out/'final-source-paths.json').read_text())
for kind in ['production','test']:(out/f'{kind}-source-path').write_text(paths[kind]+'\n')
PY
python3 "$out/source-audit.py" > "$out/source-audit.log"
python3 "$out/audit-dependency.py" > "$out/dependency-audit.log"
python3 "$out/retain.py"
python3 "$out/provenance.py" > "$out/provenance.log"
python3 "$out/package-evidence.py" > "$out/package-evidence.log"
cat "$out/provenance.log" "$out/package-evidence.log"
