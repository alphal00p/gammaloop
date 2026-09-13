#!/usr/bin/env bash
set -euo pipefail
: "${X2_TARGET:?Set the exact gated target directory containing deps}"
: "${X2_API_FINGERPRINT:?Set the current built API fingerprint JSON}"
cd "$(dirname "$0")"
python3 - <<'PYCODE'
import hashlib,json,os,subprocess
from pathlib import Path
root=Path(os.environ['X2_TARGET']); fingerprint=Path(os.environ['X2_API_FINGERPRINT'])
api=json.loads(fingerprint.read_text()); wanted={name:dep_hash for _,name,_,dep_hash in api['deps']}
libs={'gammaloop_api':root/'deps/libgammaloop_api.rlib'}
for name in ['gammalooprs','figment','toml','color_eyre','serde_json','symbolica']:
    matches=[]
    for p in (root/'.fingerprint').glob('*/lib-'+name):
        if int.from_bytes(bytes.fromhex(p.read_text().strip()),'little')==wanted[name]:
            suffix=p.parent.name.rsplit('-',1)[1]
            lib=root/f'deps/lib{name}-{suffix}.rlib'
            if lib.is_file(): matches.append(lib)
    assert len(matches)==1,(name,matches)
    libs[name]=matches[0]
Path('reference_libraries.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p)+'\n' for p in libs.values()))
args=['rustc','--edition=2021','-C','opt-level=1','-C','debuginfo=0','-L','dependency='+str(root/'deps')]
for name,path in libs.items(): args+=['--extern',f'{name}={path}']
args+=['gl638_x2_direct_h_driver.rs','-o','gl638_x2_direct_h_driver']
Path('rustc-command.json').write_text(json.dumps(args,indent=2)+'\n')
subprocess.run(args,check=True)
Path('reference_driver.sha256').write_text(''.join(hashlib.sha256(Path(name).read_bytes()).hexdigest()+'  '+name+'\n' for name in ['gl638_x2_direct_h_driver','gl638_x2_direct_h_driver.rs']))
PYCODE
