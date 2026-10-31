from pathlib import Path
exec((Path(__file__).resolve().parent/'pure_stage.py').read_text().split('# Bind the exact immutable new code')[0])
outputs=json.loads((P/'producer-outputs.json').read_text())
for name,item in outputs.items():
 if name=='ordinary':graph=json.loads((P/'pure/ordinary-graph.json').read_text())
 else:graph=json.loads(run(name+'-graph',[NIX,'derivation','show','-r',item['drv']]));save('pure/'+name+'-graph.json',graph)
 save('pure/'+name+'-root-derivation.json',graph['derivations'][Path(item['drv']).name])
 run(name+'-dry-run',[NIX,'build','--offline','--dry-run','--no-link','--json',item['drv']+'^out'],cwd='/tmp')
get=lambda n:(P/'pure'/n/'stdout').read_text().strip()
info=json.loads((P/'pure/source-evaluation.json').read_text())
full=ev('full-source-evaluation','let f=builtins.getFlake "git+file:///shared/lcnbr/dev/gammaloop/main?rev='+CODE+'&shallow=1&allRefs=1"; in toString f.outPath')
save('source-bindings.json',{'new_code':CODE,'old_host_code':OLD,'full_source':full,'full_source_nar':nar('full-source-nar',full),'ordinary_source':info['source'],'ordinary_nar':get('ordinary-source-nar'),'private_source':get('seal-private-source'),'private_nar':get('private-source-nar'),'rustred_source':'/nix/store/5hr22rycs6dpj77afgq0maqd8qbipk6a-rustred-source-7c1ed037','rustred_nar':get('rustred-source-nar'),'proposal_source':get('seal-proposal'),'proposal_nar':get('proposal-nar'),'source_filter':str(P/'proposal/ordinary-base.nix'),'producer_outputs':outputs,'metadata_gate':'offline locked before each compilation','no_runtime_execution':True,'secret_variables_omitted':True,'retained_setup_failure':'Nix v4 graph has derivations wrapper and basename keys; adapter corrected without changing sources or producer inputs'})
print(json.dumps(outputs),flush=True)
