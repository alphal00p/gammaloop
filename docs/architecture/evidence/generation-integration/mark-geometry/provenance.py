from pathlib import Path
import json,subprocess,hashlib,tomllib,re
root=Path.cwd();out=root/'target/mark-geometry';base='100da20603a083d031c31bd4e8c6970fc44dda2f';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();read=lambda rev,path:subprocess.check_output(['jj','--ignore-working-copy','file','show','-r',rev,path])
assert read(base,'crates/kurvst/src/marks.rs')==read('15dd7e1b-','crates/kurvst/src/marks.rs')
derivations=json.loads((out/'derivation.json').read_text())['derivations'];name,drv=next(iter(derivations.items()));source=Path(drv['env']['src']);bundle=(root/'result').resolve();sources=[]
for p in sorted(source.rglob('*')):
 if p.is_file():
  rel=p.relative_to(source).as_posix();assert (root/rel).is_file(),rel;assert p.read_bytes()==(root/rel).read_bytes(),rel;sources.append({'path':rel,'sha256':sha(p)})
payload=[]
for p in sorted(bundle.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(bundle).as_posix()
 if rel in ['linnest.wasm','kurvst.wasm']:owner=root/'crates'/p.stem/'typst'/p.name
 elif rel.startswith('templates/crates/'):owner=root/rel.removeprefix('templates/')
 else:owner=root/'crates/clinnet'/rel
 assert p.read_bytes()==owner.read_bytes(),rel;payload.append({'path':rel,'source':owner.relative_to(root).as_posix(),'sha256':sha(p)})
assert not any('tiptoe' in x['path'].lower() for x in payload)
assert tomllib.loads((bundle/'templates/crates/linnest/typst/typst.toml').read_text())['package']['license']=='MIT AND Apache-2.0 AND BSD-3-Clause'
notice_paths=['vendor/clarabel/LICENSE.md','crates/linnest/layout/native/LICENSE.networkx','crates/linnest/typst/LICENSE.clarabel','crates/linnest/typst/LICENSE.networkx','crates/kurvst/typst/LICENSE','crates/linnest/typst/LICENSE','examples/notebooks/symbolica-host/licenses/clarabel-LICENSE','examples/notebooks/symbolica-host/licenses/networkx-LICENSE','crates/linnest/build_assets.rs']
for path in notice_paths:assert (root/path).read_bytes()==read(base,path),path
allfiles=subprocess.check_output(['jj','file','list','-r',base],text=True).splitlines();protected=[p for p in allfiles if Path(p).name in ['Cargo.lock','uv.lock','flake.lock','nix-ci.nix'] or p.endswith('.pyi')]
for path in protected:assert (root/path).read_bytes()==read(base,path),path
oldbookmark=subprocess.check_output(['jj','log','-r','generation-drawing-ports','--no-graph','-T','commit_id'],text=True).strip();assert oldbookmark=='98cbf5b8f984e702f4699b1b176fbebfee059d3a'
metadata=json.loads((out/'native-test-binaries.json').read_text());binaries=[]
for key,v in metadata['rust-binaries'].items():binaries.append({**v,'sha256':sha(Path(v['binary-path']))})
def varint(data,at):
 value=0;shift=0
 while True:
  b=data[at];at+=1;value|=(b&127)<<shift
  if b<128:return value,at
  shift+=7
sections={}
for package in ['linnest','kurvst']:
 p=bundle/f'{package}.wasm';data=p.read_bytes();assert data[:8]==b'\0asm\x01\0\0\0';at=8;items=[]
 while at<len(data):
  ident=data[at];at+=1;size,start=varint(data,at);body=data[start:start+size];name=''
  if ident==0:n,first=varint(body,0);name=body[first:first+n].decode(errors='replace')
  items.append({'id':ident,'name':name,'size':size,'sha256':hashlib.sha256(body).hexdigest()});at=start+size
 sections[package]={'sha256':sha(p),'sections':items}
package_info=json.loads((out/'package.json').read_text())
report={'base':base,'source38':'15dd7e1b188dad112d7887bdc6d769125b2697a1','source38_parent_marks_equals_base':True,'nix_derivation':'/nix/store/'+next(iter(derivations)),'nix_source':str(source),'nix_output':str(bundle),'nix_package':package_info,'cargo_dependency_artifacts':drv['env'].get('cargoArtifacts'),'wasm_native_build_inputs':drv['env'].get('nativeBuildInputs'),'cargo_build_phase':drv['env'].get('buildPhase'),'source_files':sources,'source_count':len(sources),'payload':payload,'payload_count':len(payload),'wasm':sections,'native_tests':binaries,'notice_hashes':{p:sha(root/p) for p in notice_paths},'protected_files':{p:sha(root/p) for p in protected},'old_drawing_bookmark':oldbookmark,'producer_rebuild_exit':int((out/'wasm-reproduction.exit').read_text()),'reproduction_scope':'Independent Nix final producer rebuild matched both complete WASM/bundle outputs. Source-qualified dependency artifacts were reused normally. This is not a cold-dependency reproducibility proof, nor a retroactive repair/certification of historical Linnest data/name differences.'}
(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'{len(sources)} exact Nix sources; {len(payload)} exact delivered files; {len(protected)} locks/stubs/configs preserved; notice selector and all retained notices unchanged')
# All output changes are observed, not imported historical blobs.
asset_paths=[p.relative_to(root).as_posix() for p in (root/'docs/assets').rglob('*.svg')]+[p.relative_to(root).as_posix() for p in (root/'assets').glob('gammalooplogo*.svg')]+['crates/kurvst/typst/kurvst.wasm','crates/linnest/typst/linnest.wasm','crates/kurvst/typst/docs/manual.pdf']+[f'crates/linnest/typst/examples/{p}.pdf' for p in ['soft-triangle','thrown-out','xbox-sum']]
assets=[]
for path in sorted(asset_paths):
 try:before=read(base,path)
 except subprocess.CalledProcessError:continue
 assets.append({'path':path,'before_sha256':hashlib.sha256(before).hexdigest(),'sha256':sha(root/path),'changed':before!=(root/path).read_bytes()})
(out/'asset-inventory.json').write_text(json.dumps(assets,indent=2)+'\n')
tiptoe=Path('/nix/store/v7rw00hg1wwrr7gcl128yh9nw4v3dgrg-typst-package-tiptoe-0.4.0/lib/typst-packages/tiptoe/0.4.0')
reference={p.relative_to(tiptoe).as_posix():sha(p) for p in sorted(tiptoe.rglob('*')) if p.is_file()}
retired=[]
for path in ['crates/kurvst/typst/LICENSE.tiptoe','examples/notebooks/symbolica-host/licenses/tiptoe-LICENSE','crates/linnet-py/licenses/tiptoe-LICENSE']:
 prior=subprocess.run(['jj','--ignore-working-copy','file','show','-r',base,path],capture_output=True)
 retired.append({'path':path,'present_in_base':prior.returncode==0,'prior_sha256':hashlib.sha256(prior.stdout).hexdigest() if prior.returncode==0 else None,'disposition':'removed with replaced geometry' if prior.returncode==0 else 'already absent; standalone packaging not restored'})
references=[]
for path in subprocess.check_output(['jj','file','list'],text=True).splitlines():
 p=root/path
 if p.suffix not in ['.rs','.typ','.py','.nix','.toml','.md','.just'] or '/evidence/' in path:continue
 for no,line in enumerate(p.read_text().splitlines(),1):
  if 'tiptoe' in line.lower():references.append({'path':path,'line':no,'text':line.strip(),'kind':'negative payload/upgrade control' if path.startswith(('crates/','tests/','nix/','just/')) else 'historical provenance, source-bound record or current audit explanation'})
audit={'reference_package':str(tiptoe),'reference_package_files':reference,'reference_license_sha256':sha(tiptoe/'LICENSE.txt'),'actual_provenance':'Ported original source38, whose parent marks.rs exactly equals this approved base. Inspected original/current mark code and pinned Tiptoe marks/line/curve/utility/license sources. No Tiptoe source was downloaded into or vendored with the package.','replacement_constructions':{'triangle_stealth':'Generic offset-line intersections/corner bisectors and edge-direction tests replace Tiptoe trigonometric x1/x2/x3 special cases.','straight':'Tangent from the tip to round endpoint caps replaces old fixed-angle inset.','round_tikz':'Inner-box dart and two Kurbo quarter arcs replace old L-thickness special case and 0.42/0.32/0.1 control points.','barb_hooks_circle':'Kurbo Arc/Ellipse geometry from local centre/radius definitions and stroke growth; no Tiptoe path/arc renderer retained. Elementary shape meaning and public arc fields remain.','bar_bracket_square_diamond':'Local symmetric lines/polygon vertex definitions, generic offsets/miter and exact sqrt(2) diamond default. Elementary endpoints retain shape semantics, not a retained material Tiptoe implementation.','rays_composite':'Evenly spaced vector spokes and typed Rust composite accumulation through existing MarkSpec owners; explicit public phase and signed gaps retained.','transport_fitting':'Existing Rust CBOR schemas, packed paths, Kurbo carrier/chord/refinement/stroke footprint owners preserved. They are not copied Tiptoe rendering/resource code.'},'retained_compatibility':'Shape names, option vocabulary and dimensional default data retain the earlier Tiptoe-inspired public API. This is not a clean-room claim. Engineering source audit found no material Tiptoe package implementation/resource remaining in runtime geometry; old source/credits remain in immutable history and historical notes for independent review.','retired':retired,'retained_notice_sha256':report['notice_hashes'],'remaining_text_references':references,'source_owner_reconciliation':['Current Clinnet/API embed selectors are generic; only obsolete positive assertions become absence controls.','Common build_assets.rs root-license selector unchanged byte-for-byte.','Both actual sync recipe helpers remove only known obsolete generated license and preserve missing-license upgrade, idempotence and unrelated user notice controls.','Existing host license-files glob and root sdist owner retained; absent standalone wheel not revived.','Nix source filters drop only obsolete Tiptoe paths; check asserts absence plus exact retained notice bytes and approved expression.','Current API prose/README already reflects Spec/Community owners and contains no obsolete attribution, so not overwritten by historical patch.','Historical frozen architecture/receipt text remains historical; current prefaces point at this separate implementation.']}
(out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
