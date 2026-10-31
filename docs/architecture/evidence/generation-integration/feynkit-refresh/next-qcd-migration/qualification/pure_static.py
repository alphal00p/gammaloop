"""Pinned unlicensed static checks and existing Typst tool resolution; no build."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent));import pure as P
O=P.O;W=P.W
compiler=json.loads((W/'.local/feynkit-refresh/grouped-identical-migration-20261031/archive-qualification.json').read_text())['compiler']
P.call('rustfmt',[compiler['rustfmt']['path'],'--check','--edition','2024','crates/feynkit-generator/tests/generation_oracle/main.rs'],False)
P.call('code-diff-check',['git','--git-dir='+str(W.parent/'main/.git'),'diff','--check','83276dc047815d1dcc7f003aeb512c204b26758e',P.REV,'--','crates/feynkit-generator/tests/generation_oracle/main.rs'],False)
expr='let f=builtins.getFlake '+json.dumps(P.URI)+'; in map (p: { name=p.name or ""; path=toString p; drv=p.drvPath; }) (builtins.filter (p: builtins.match "typst.*" (p.name or "") != null) f.devShells.x86_64-linux.default.buildInputs)'
tools=P.nix('typst-resolution',['eval','--json','--expr',expr]);P.write('typst-tools.json',tools)
print('PASS pinned rustfmt and exact immutable code diff; Typst candidates',tools)
