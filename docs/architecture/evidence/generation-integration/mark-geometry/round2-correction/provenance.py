from pathlib import Path
import hashlib,json,subprocess,tomllib,re
root=Path.cwd();out=root/'target/mark-geometry-round3';base='e426ff5fcbe1487c9bf2cf1ed5aaaa03cabf8722';code='f35643bae1a67fb6ec0c30151e3624f661cf5f96';read=lambda rev,p:subprocess.check_output(['jj','--ignore-working-copy','file','show','-r',rev,p]);sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();digest=lambda b:hashlib.sha256(b).hexdigest()
drvname,drv=next(iter(json.loads((out/'clean-derivation.json').read_text())['derivations'].items()));src=Path(drv['env']['src']);bundle=(out/'final-result').resolve();sources=[]
for p in sorted(src.rglob('*')):
 if p.is_file():
  rel=str(p.relative_to(src));assert (root/rel).read_bytes()==p.read_bytes()==read(code,rel),rel;sources.append({'path':rel,'sha256':sha(p)})
payload=[]
for p in sorted(bundle.rglob('*')):
 if p.is_file():
  rel=p.relative_to(bundle).as_posix()
  owner=(root/'crates'/p.stem/'typst'/p.name) if rel in ['linnest.wasm','kurvst.wasm'] else root/(rel.removeprefix('templates/') if rel.startswith('templates/crates/') else 'crates/clinnet/'+rel)
  assert p.read_bytes()==owner.read_bytes()==read(code,str(owner.relative_to(root))),rel;payload.append({'path':rel,'owner':str(owner.relative_to(root)),'sha256':sha(p)})
assert len(payload)==29 and len(sources)==1455,(len(payload),len(sources));assert not any('tiptoe' in x['path'].lower() for x in payload)
assert tomllib.loads((bundle/'templates/crates/linnest/typst/typst.toml').read_text())['package']['license']=='MIT AND Apache-2.0 AND BSD-3-Clause'
notices=['vendor/clarabel/LICENSE.md','crates/linnest/layout/native/LICENSE.networkx','crates/linnest/typst/LICENSE.clarabel','crates/linnest/typst/LICENSE.networkx','crates/kurvst/typst/LICENSE','crates/linnest/typst/LICENSE','examples/notebooks/symbolica-host/licenses/clarabel-LICENSE','examples/notebooks/symbolica-host/licenses/networkx-LICENSE','crates/linnest/build_assets.rs','just/drawing.just']
for p in notices:assert (root/p).read_bytes()==read(base,p),p
allfiles=subprocess.check_output(['jj','--ignore-working-copy','file','list','-r',base],text=True).splitlines();protected=[p for p in allfiles if (Path(p).name in ['Cargo.lock','uv.lock','flake.lock','nix-ci.nix'] or p.endswith('.pyi')) and p not in ['Cargo.lock','examples/notebooks/symbolica-host/Cargo.lock']]
for p in protected:assert (root/p).read_bytes()==read(base,p),p
assert re.search(r'\benable\s*=\s*false',Path('nix-ci.nix').read_text())
for p in ['nix/ci-workspace-graph.json','crates/gammaloop-workspace-hack/Cargo.toml','nix-ci.nix']:assert (root/p).read_bytes()==read(base,p)
for ref,expected in [('generation-drawing-ports','98cbf5b8f984e702f4699b1b176fbebfee059d3a'),('generation-mark-geometry-round1','cc7bc03ab97edade07532a586f5830f933451b07'),('generation-mark-geometry-round2',base)]:assert subprocess.check_output(['jj','log','-r',ref,'--no-graph','-T','commit_id'],text=True).strip()==expected
binaries=[]
for filename in ['native-test-binaries.json','clinnet-binaries.json']:
 for key,v in json.loads((out/filename).read_text())['rust-binaries'].items():binaries.append({**v,'sha256':sha(Path(v['binary-path']))})
assets=[]
paths=sorted([p.relative_to(root).as_posix() for p in (root/'docs/assets').rglob('*.svg')]+[p.relative_to(root).as_posix() for p in (root/'assets').glob('gammalooplogo*.svg')]+['crates/kurvst/typst/kurvst.wasm','crates/linnest/typst/linnest.wasm','crates/kurvst/typst/docs/manual.pdf']+[f'crates/linnest/typst/examples/{p}.pdf' for p in ['soft-triangle','thrown-out','xbox-sum']])
for p in paths:
 old=read(base,p);assets.append({'path':p,'before_sha256':digest(old),'sha256':sha(root/p),'changed':old!=(root/p).read_bytes()})
assert len(assets)==40;assert [p['path'] for p in assets if p['changed']]==['crates/kurvst/typst/kurvst.wasm','crates/linnest/typst/linnest.wasm']
(out/'asset-inventory.json').write_text(json.dumps(assets,indent=2)+'\n')
galleries=[]
for gallery in ['marks','boundaries']:
 for theme in ['light','dark']:
  for fmt in ['svg','png','pdf']:
   p=out/f'{gallery}-{theme}.{fmt}';old=root/'docs/architecture/evidence/generation-integration/mark-geometry/round1-correction'/p.name;reproduced=out/f'{gallery}-{theme}-reproduced.{fmt}';assert p.read_bytes()==reproduced.read_bytes();galleries.append({'path':p.name,'before_sha256':sha(old),'sha256':sha(p),'changed':p.read_bytes()!=old.read_bytes(),'reproduced':True})
(out/'gallery-inventory.json').write_text(json.dumps(galleries,indent=2)+'\n')
changed=subprocess.check_output(['jj','--ignore-working-copy','diff','--from',base,'--to',code,'--name-only'],text=True).splitlines();inventory=[{'path':p,'sha256':sha(root/p),'new':p not in allfiles} for p in changed];(out/'changed-code-files.json').write_text(json.dumps(inventory,indent=2)+'\n')
report={'base':base,'corrected_code':code,'code_change':'vvwwnloulqpnmrmxpkyspwunkwsxsspu','source38':'15dd7e1b188dad112d7887bdc6d769125b2697a1','clean_build_input_commit':(out/'build-source-commit').read_text().strip(),'clean_build_flake':(out/'clean-flake').read_text().strip(),'clean_source_files_match_final_code':True,'nix_derivation':'/nix/store/'+drvname,'nix_source':str(src),'nix_output':str(bundle),'nix_package':json.loads((out/'clean-package.json').read_text()),'cargo_artifacts':drv['env'].get('cargoArtifacts'),'source_count':len(sources),'source_files':sources,'payload_count':len(payload),'payload':payload,'native_test_binaries':binaries,'retained_notices_and_selector':{p:sha(root/p) for p in notices},'protected_files':{p:sha(root/p) for p in protected},'generated_hakari_and_ci_files_unchanged':True,'independent_final_producer_rebuild_exit':int((out/'clean-reproduction.exit').read_text()),'reproducibility_limit':'Only a bounded exact-source final-producer rebuild, with normally reused source-qualified dependency artifacts. No cold dependency, historical WASM determinism, browser/full-site, final native/Pyodide/Community wheels or physics certification.'}
(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n');print(f'{len(sources)} source files and {len(payload)} delivered files match exact final code; 40 assets with only 2 WASMs changed; protected refs/files/notices unchanged')
