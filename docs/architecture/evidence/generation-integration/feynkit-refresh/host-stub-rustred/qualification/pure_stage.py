from pathlib import Path
import os, subprocess, json, hashlib, shutil, time, datetime
W=Path.cwd(); P=W/'.local/feynkit-refresh/host-stub-rustred-integration-20261031-r1'; PRE=W/'.local/feynkit-refresh/rustred-host-preflight-20261031-r1'
ENV={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')};ENV['GIT_OPTIONAL_LOCKS']='0'
NIX='/run/current-system/sw/bin/nix'; CODE='2cfb52b28d99f92305a9d571bbe3911cf04b8a19'; OLD='f47f27213a1c55291239c9bc0f334a7ea1cd1cc1'
(P/'pure').mkdir(exist_ok=True);(P/'proposal').mkdir(exist_ok=True)
def save(n,o): p=P/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2)+'\n')
def run(label,args,cwd=W):
 t=time.monotonic();r=subprocess.run(args,cwd=cwd,env=ENV,capture_output=True);d=P/'pure'/label;d.mkdir(exist_ok=False);(d/'stdout').write_bytes(r.stdout);(d/'stderr').write_bytes(r.stderr);save(str(d.relative_to(P)/'command.json'),{'argv':args,'cwd':str(cwd),'exit_code':r.returncode,'elapsed_s':time.monotonic()-t,'private_variables_absent':True});assert r.returncode==0,(label,r.stderr.decode(errors='replace')[-4000:]);return r.stdout

def ev(label,expr): return json.loads(run(label,[NIX,'eval','--offline','--option','pure-eval','true','--option','allow-import-from-derivation','false','--json','--expr',expr]))
def nar(label,path): return run(label,[NIX,'hash','path',str(path)]).decode().strip()
def manifest(root):
 root=Path(root); return {str(p.relative_to(root)):{'kind':'symlink' if p.is_symlink() else 'file','sha256':hashlib.sha256(os.readlink(p).encode() if p.is_symlink() else p.read_bytes()).hexdigest(),'executable':bool(p.stat().st_mode&0o111)} for p in sorted(root.rglob('*')) if p.is_file() or p.is_symlink()}
def delta(a,b): return {'changed':[k for k in a.keys()&b.keys() if a[k]!=b[k]],'removed':sorted(a.keys()-b.keys()),'added':sorted(b.keys()-a.keys())}
# Bind the exact immutable new code before selecting producers.
changed=run('code-diff',['git','-C',str(W.parent/'main'),'diff','--name-only','d0fb4e906140a6a0fc8dd6401ca839cab0c2734d',CODE]).decode().splitlines();assert changed==['examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi']
orig=(W/'.local/feynkit-refresh/installed-host-preflight-20261030-r1/bounded-producers.nix').read_text();new=orig.replace(OLD,CODE);(P/'proposal/ordinary-base.nix').write_text(new)
expr='let p = import '+str(P/'proposal/ordinary-base.nix')+'; in {source=toString p.selectedSource; normal=toString p.hostWheel.src; base=p.hostWheel.cargoArtifacts; vendor=p.hostWheel.cargoVendorDir;}'
# Seal expression sources so pure evaluator never reads mutable workspace paths.
sealedbase=run('seal-base',[NIX,'store','add-path',str(P/'proposal/ordinary-base.nix'),'--name','host-stub-ordinary-base.nix']).decode().strip();h=nar('base-nar',sealedbase)
expr='let p = import (builtins.path {path='+sealedbase+'; sha256="'+h+'"; name="host-stub-ordinary-base.nix";}); in {source=toString p.selectedSource; base=toString p.hostWheel.cargoArtifacts; vendor=toString p.hostWheel.cargoVendorDir;}'
info=ev('source-evaluation',expr);save('pure/source-evaluation.json',info);source=Path(info['source']);srch=nar('ordinary-source-nar',source)
assert hashlib.sha256((source/'examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi').read_bytes()).hexdigest()=='5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5'
oldsrc=Path('/nix/store/0q772x7ilixakyvvd43h9z4119409f9h-source');a=manifest(oldsrc);b=manifest(source);d=delta(a,b);assert sorted(d['changed'])==changed and not d['removed'] and not d['added'],d;save('ordinary-source-manifest.json',b);save('ordinary-vs-f47-source-delta.json',{'old_source':str(oldsrc),'old_nar':nar('old-source-nar',oldsrc),'new_source':str(source),'new_nar':srch,**d,'native_and_generator_sources_byte_invariant':True})
rr=Path('/nix/store/5hr22rycs6dpj77afgq0maqd8qbipk6a-rustred-source-7c1ed037');rrh=nar('rustred-source-nar',rr);assert rrh=='sha256-6Sgm1lGDDGCe2ch784rGp0krsqlOiZtKAEuafvNnKlA='
scratch=P/'proposal/private-source-input';shutil.copytree(source,scratch,symlinks=True)
for p in scratch.rglob('*'):
 if not p.is_symlink():p.chmod(0o755 if p.is_dir() else 0o755 if p.stat().st_mode&0o111 else 0o644)
scratch.chmod(0o755)
patch=(PRE/'proposal/private-host.patch').read_bytes();r=subprocess.run(['patch','--batch','-p1','--directory',str(scratch)],input=patch,env=ENV,cwd=W,capture_output=True);(P/'pure/private-patch.stdout').write_bytes(r.stdout);(P/'pure/private-patch.stderr').write_bytes(r.stderr);assert r.returncode==0
c=manifest(scratch);d=delta(b,c);assert sorted(d['changed'])==['examples/notebooks/symbolica-host/Cargo.lock','examples/notebooks/symbolica-host/Cargo.toml','examples/notebooks/symbolica-host/src/lib.rs'] and not d['removed'] and not d['added'];save('private-source-manifest.json',c)
priv=run('seal-private-source',[NIX,'store','add-path',str(scratch),'--name','source']).decode().strip();privh=nar('private-source-nar',priv);shutil.rmtree(scratch);save('private-source-delta.json',{'base_code':CODE,'ordinary_source':str(source),'private_source':priv,'private_nar':privh,'reviewed_patch_sha256':hashlib.sha256(patch).hexdigest(),**d})
producer=(PRE/'proposal/private-producer.nix').read_text().replace(OLD,CODE).replace('/nix/store/3rgbf8hy7dfaarqah5s1dm1blhf1a37w-source',priv).replace('sha256-Xab4epwGb9q004L3u9QukLcIwKfpy5d3UYWTpRUOEUU=',privh)
# Retain final feature-context fingerprints and build compiler identity, no runtime.
producer=producer.replace('runHook preBuild\n      cargo metadata','runHook preBuild\n      rustc -Vv > producer-rustc.txt\n      cargo metadata').replace('cp private-host-metadata.json "$out/provenance/"','cp private-host-metadata.json producer-rustc.txt "$out/provenance/"')
(P/'proposal/private-producer.nix').write_text(producer)
ordinary=producer.replace(priv,str(source)).replace(privh,srch).replace('generation-refresh-rustred-validation-host-wheel','generation-refresh-corrected-stub-host-wheel').replace('    RUSTRED_SOURCE = rustredSource;\n','').replace('    RUSTRED_RUNTIME_ARITIES = "1";\n','')
start=ordinary.index('      for name in ("rustred-feynkit", "rustred", "rustred-order"):');end=ordinary.index('\n      METADATA',start);ordinary=ordinary[:start]+'      assert not any(p["name"].startswith("rustred") for p in packages)'+ordinary[end:]
ordinary=ordinary.replace('      printf \'%s\\n\' "$RUSTRED_RUNTIME_ARITIES" > "$out/provenance/runtime-arities.txt"\n','')
(P/'proposal/ordinary-producer.nix').write_text(ordinary)
# Only two wheel attributes are selected. No exporter/archive/site derivations realized.
sealed=run('seal-proposal',[NIX,'store','add-path',str(P/'proposal'),'--name','host-stub-rustred-proposal']).decode().strip();sealedh=nar('proposal-nar',sealed)
expr='let s=builtins.path {path='+sealed+';name="host-stub-rustred-proposal";sha256="'+sealedh+'";}; o=import (s+"/ordinary-producer.nix"); p=import (s+"/private-producer.nix"); in {ordinary={drv=o.hostWheel.drvPath;out=o.hostWheel.outPath;};private={drv=p.hostWheel.drvPath;out=p.hostWheel.outPath;};}'
(P/'eval-sealed.nix').write_text(expr+'\n');outputs=ev('corrected-pure-evaluation',expr);save('producer-outputs.json',outputs)
for name,item in outputs.items():
 graph=json.loads(run(name+'-graph',[NIX,'derivation','show','-r',item['drv']]));save('pure/'+name+'-graph.json',graph)
 save('pure/'+name+'-root-derivation.json',graph[item['drv']])
 run(name+'-dry-run',[NIX,'build','--offline','--dry-run','--no-link','--json',item['drv']+'^out'],cwd='/tmp')
save('source-bindings.json',{'new_code':CODE,'old_host_code':OLD,'ordinary_source':str(source),'ordinary_nar':srch,'private_source':priv,'private_nar':privh,'rustred_source':str(rr),'rustred_nar':rrh,'proposal_source':sealed,'proposal_nar':sealedh,'source_filter':str(P/'proposal/ordinary-base.nix'),'producer_outputs':outputs,'metadata_gate':'offline locked before each compilation','no_runtime_execution':True,'secret_variables_omitted':True})
print(json.dumps(outputs),flush=True)
