from pathlib import Path
exec((Path(__file__).resolve().parent/'pure_stage.py').read_text().split('# Bind the exact immutable new code')[0])
import tomllib
b=json.loads((P/'source-bindings.json').read_text());d=json.loads((P/'lock-reconciliation.json').read_text());assert d['added']==[] and d['removed']==[['rustred-app','0.1.0',None],['rustred-python','0.1.0',None]];assert len(d['changed'])==1 and d['changed'][0]['key']==['rustred-feynkit','0.1.0',None]
m=json.loads((P/'lock-resolver.stdout').read_text());nodes={n['id']:n for n in m['resolve']['nodes']}
for name in ('rustred','rustred-order','rustred-feynkit'):
 p,=[p for p in m['packages'] if p['name']==name];assert nodes[p['id']]['features']==[]
assert not any(p['name'] in ('rustred-app','rustred-python') for p in m['packages'])
scratch=P/'private-lock-work';a=manifest(b['private_source']);c=manifest(scratch);delta0=delta(a,c);assert delta0=={'changed':['examples/notebooks/symbolica-host/Cargo.lock'],'removed':[],'added':[]},delta0
shutil.copytree(P/'proposal',P/'proposal-attempt1')
priv=run('r2-seal-private-source',[NIX,'store','add-path',str(scratch),'--name','source']).decode().strip();h=nar('r2-private-source-nar',priv)
s=(P/'proposal/private-producer.nix').read_text().replace(b['private_source'],priv).replace(b['private_nar'],h);(P/'proposal/private-producer.nix').write_text(s)
sealed=run('r2-seal-proposal',[NIX,'store','add-path',str(P/'proposal'),'--name','host-stub-rustred-proposal']).decode().strip();nh=nar('r2-proposal-nar',sealed)
expr='let s=builtins.path {path='+sealed+';name="host-stub-rustred-proposal";sha256="'+nh+'";}; o=import (s+"/ordinary-producer.nix"); p=import (s+"/private-producer.nix"); in {ordinary={drv=o.hostWheel.drvPath;out=o.hostWheel.outPath;};private={drv=p.hostWheel.drvPath;out=p.hostWheel.outPath;};}'
(P/'eval-sealed-r2.nix').write_text(expr+'\n');outputs=ev('r2-corrected-pure-evaluation',expr);assert outputs['ordinary']==b['producer_outputs']['ordinary'];save('producer-outputs-r2.json',outputs)
graph=json.loads(run('r2-private-graph',[NIX,'derivation','show','-r',outputs['private']['drv']]));save('pure/r2-private-graph.json',graph);save('pure/r2-private-root-derivation.json',graph['derivations'][Path(outputs['private']['drv']).name]);run('r2-private-dry-run',[NIX,'build','--offline','--dry-run','--no-link','--json',outputs['private']['drv']+'^out'],cwd='/tmp')
save('source-bindings-r2.json',{**b,'private_source':priv,'private_nar':h,'proposal_source':sealed,'proposal_nar':nh,'producer_outputs':outputs,'lock_correction':{'reason':'Cargo canonical lock prunes two disabled optional path adapters and their unused bridge dependency; no version/source changed','delta':d,'locked_gate_preserved':True}})
shutil.rmtree(scratch);shutil.rmtree(P/'private-cargo-home');shutil.rmtree(P/'empty-home')
print(json.dumps(outputs),flush=True)
