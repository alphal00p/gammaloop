from pathlib import Path,PurePosixPath
import os,subprocess,json,hashlib,zipfile,csv,io,base64,struct,re,tomllib,ast
P=Path(__file__).resolve().parent;W=P.parents[2];H=P.parent/'installed-host-qualification-20261030-r1';R=W/'.local/review-host-browser-stub-astra-20261031-013106Z-5655cc61';B=P.parent/'host-browser-correction-20261031-r1'
E={k:v for k,v in os.environ.items() if k not in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')}
def save(n,v):p=P/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def sha(p):return hashlib.file_digest(Path(p).open('rb'),'sha256').hexdigest()
def run(n,a,allow=False):
 p=P/'qualification'/n;p.mkdir(parents=True,exist_ok=False);r=subprocess.run(a,cwd='/tmp',env=E,capture_output=True);(p/'stdout').write_bytes(r.stdout);(p/'stderr').write_bytes(r.stderr);save(str(p.relative_to(P)/'command.json'),{'argv':a,'cwd':'/tmp','exit_code':r.returncode,'private_variables_absent':True});assert allow or r.returncode==0,(n,r.stderr.decode()[-3000:]);return r.stdout
src=(H/'inspect_artifacts.py').read_text();exec(src[src.index('def elf('):src.index('def scan_tar(')])
b=json.loads((P/'source-bindings-r2.json').read_text());outputs=b['producer_outputs'];nix='/run/current-system/sw/bin/nix';records={}
for attempt in ('','-r2'):
 for k,v in json.loads((P/('producer-outputs'+attempt+'.json')).read_text()).items():run('log'+attempt+'-'+k,[nix,'log',v['drv']],True)
for name,item in outputs.items():
 if not Path(item['out']).exists():continue
 wheel,=list((Path(item['out'])/'dist').glob('*.whl'));members=[];native=[];metadata={};extracted=P/'artifacts'/name
 with zipfile.ZipFile(wheel) as z:
  assert len(z.namelist())==len(set(z.namelist()))
  for f in z.namelist():
   q=PurePosixPath(f);assert not q.is_absolute() and '..' not in q.parts;data=z.read(f);members.append({'path':f,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
   if f.endswith(('/METADATA','/WHEEL','/RECORD')):metadata[f]=data.decode()
   if f.endswith('.so'):
    p=extracted/q;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data);native.append(elf(p))
  record,=[v for k,v in metadata.items() if k.endswith('/RECORD')]
  for f,h,size in csv.reader(io.StringIO(record)):
   data=z.read(f)
   if h:assert h=='sha256='+base64.urlsafe_b64encode(hashlib.sha256(data).digest()).decode().rstrip('=')
   if size:assert len(data)==int(size)
  stub='symbolica/community/hepkit/__init__.pyi';assert hashlib.sha256(z.read(stub)).hexdigest()=='5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5'
 assert len(native)==1 and native[0]['path'].endswith('/symbolica/core.abi3.so') and not native[0]['nonstore_library_paths']
 provenance=Path(item['out'])/'provenance';m=json.loads((provenance/'private-host-metadata.json').read_text());nodes={n['id']:n for n in m['resolve']['nodes']};packages=m['packages'];features={}
 for n in ('symbolica','numerica','graphica','pyo3','pyo3-ffi','feynkit-py','spynso3','linnet-py','rustred','rustred-order','rustred-feynkit'):
  ps=[p for p in packages if p['name']==n]
  if not ps:continue
  assert len(ps)==1;n0=ps[0];features[n]={'id':n0['id'],'version':n0['version'],'source':n0['source'],'manifest_path':n0['manifest_path'],'features':nodes[n0['id']]['features']}
  if n in ('symbolica','numerica','graphica'):assert n0['version']=='3.0.1' and n0['source'].endswith('#8eb7c3729d6abefc916decc1bdf0c06536f8af80')
  if n.startswith('pyo3'):assert n0['version']=='0.28.3'
  if n.startswith('rustred'):assert nodes[n0['id']]['features']==[] and n0['manifest_path'].startswith(b['rustred_source']+'/')
 assert not any(p['name'] in ('rustred-app','rustred-python') for p in packages)
 if name=='private':assert {'rustred','rustred-order','rustred-feynkit'}<=features.keys() and (provenance/'runtime-arities.txt').read_text()=='1\n'
 else:assert not any(p['name'].startswith('rustred') for p in packages)
 prov=[{'path':str(p.relative_to(provenance)),'sha256':sha(p),'bytes':p.stat().st_size} for p in provenance.rglob('*') if p.is_file()]
 fingerprint_data=[]
 for f in (provenance/'cargo-fingerprints').rglob('*.json'):
  data=json.loads(f.read_text());fingerprint_data.append({'path':str(f.relative_to(provenance)),'data':data})
 assert any('alphal00p-notebook-host' in x['path'] for x in fingerprint_data)
 if name=='private':
  for stem in ('rustred-','rustred-feynkit-','rustred-order-'):assert any(stem in x['path'] for x in fingerprint_data),stem
 nar=run('nar-'+name,[nix,'hash','path',item['out']]).decode().strip();pathinfo=json.loads(run('path-info-'+name,[nix,'path-info','--json',item['out']]));pi=pathinfo[item['out']] if isinstance(pathinfo,dict) else pathinfo[0];assert nar==pi['narHash']
 records[name]={'drv':item['drv'],'output':item['out'],'output_nar':nar,'wheel':{'path':str(wheel),'sha256':sha(wheel),'bytes':wheel.stat().st_size},'native':native[0],'members':members,'metadata':metadata,'features':features,'provenance_files':prov,'fingerprints':fingerprint_data,'compiler':(provenance/'producer-rustc.txt').read_text(),'stub_sha256':'5e2dfa31a10d9c545d38f4257b006daa7b1b53bcf0d96aafb060ee71a68f77e5','whole_wheel_RECORD_verified':True}
 assert 'rustc 1.98.1' in records[name]['compiler']
 depfiles={f.name:f.read_text() for f in (provenance/'depfiles').glob('*.d')}
 for stem in ('symbolica','numerica','graphica'):
  ds=[v for k,v in depfiles.items() if k.startswith(stem+'-')];assert ds,stem
  assert any('/nix/store/m0kqy74hq5qsbmhx1zcgs7y21was4k8m-vendor-cargo-deps/' in v or '8eb7c3729d6abefc916decc1bdf0c06536f8af80' in v for v in ds),stem
 if name=='private':
  for stem in ('rustred','rustred_order','rustred_feynkit'):
   ds=[v for k,v in depfiles.items() if k.startswith(stem+'-')];assert ds and any(b['rustred_source']+'/' in v for v in ds),stem
  for stem in ('rustred','rustred-order','rustred-feynkit'):
   fs=[f for f in fingerprint_data if f['path'].split('/')[1].startswith(stem+'-') and '/lib-' in f['path']]
   assert fs and all(f['data']['features']=='[]' for f in fs),stem
 save('artifacts/'+name+'/compiled-source-bindings.json',{'actual_depfiles':depfiles,'source_roots_verified':True,'features_from_retained_compiler_fingerprints':fingerprint_data})
 save('artifacts/'+name+'/qualification.json',records[name])
# Same exporter binary, with exact code-owner source proof, not a new build.
exporter=Path('/nix/store/mn3vzg4c35q658jgg078hahj3hcic2f4-generation-refresh-feynkit-stub-exporter-0.1.0/bin/alphal00p-docs-python-exporter');assert sha(exporter)=='ea6b689912dd12eb1d40fc9f65f6e52ded87a33a755b94ceb42363a6f5480e9e'
full=Path(b['full_source']);oldfull=Path('/nix/store/0vsa28jjcrwyrna4iflavxm0gc8msqxm-source');proof=[]
for root in ('crates/alphal00p-docs-python-exporter','crates/feynkit-py','crates/spynso3','crates/linnet-py'):
 for p in (full/root).rglob('*'):
  if p.is_file() and not p.is_symlink():
   rel=p.relative_to(full)
   if str(rel)=='crates/feynkit-py/python/symbolica/community/hepkit/__init__.pyi':continue
   assert p.read_bytes()==(oldfull/rel).read_bytes(),rel;proof.append({'path':str(rel),'sha256':sha(p)})
save('exporter-reuse-source-proof.json',{'exporter':str(exporter),'sha256':sha(exporter),'source_owners_verified':proof,'new_build':False,'canonical_check_pending':True})
# Original QCD source/fixture and tests remain byte-exact.
checks={'crates/feynkit-py/tests/installed_feyncalc_qcd_ghost_vertex.py':'2622f4f8f53624e177e6ae059eb12953f1b16f6cf08ebe85b70ff6155423dfc6','crates/feynkit-model/tests/fixtures/sm.json':'a0193ae09a6f920c8fcfbe8231c6a70e820b82b573257e07c216bdcd8e57acd9','crates/feynkit-py/tests/installed_ufo_ghost_crossing.py':'22231e8bd4f950b386fd72a917c6468afa802ed24c2218b5c48e6d29ee44f314','crates/feynkit-py/tests/installed_squared_amplitude_display.py':'ee598de2e554ef35f5ca2166e6789df5d8d3dbb7743bcaa77ce5dd994a7587d8','crates/linnet-py/tests/test_svg_browser.py':'ee700e0629fad62358e1f0a778353757a9a305bdbc131f2ae2675b7b74d82c6c'}
for p,h in checks.items():assert sha(full/p)==h,p
q=ast.parse((full/'crates/feynkit-py/tests/installed_feyncalc_qcd_ghost_vertex.py').read_text());assert len([n for n in ast.walk(q) if isinstance(n,ast.Assert)])==23
save('original-test-bindings.json',{'source':str(full),'sha256':checks,'QCD_assertion_sites':23,'test_expectations_modified':False})
# Native resolution through immutable paths; no ambient loader search.
mounts=set(json.loads((R/'inspection/verified-inputs.json').read_text())['mounts']);mounts.discard('/nix/store/yg65yzb4ckljb4bg5w673j8ap56yrxfb-generation-refresh-installed-host-wheel-0.1.0');metadata={};edges=[];queue=[Path(x['native']['path']) for x in records.values()]+[exporter]
def carrier(p):return str(Path('/nix/store')/Path(p).parts[3])
while queue:
 p=queue.pop().resolve()
 if str(p) in metadata:continue
 m=elf(p);metadata[str(p)]=m;bad=m['nonstore_library_paths'];assert not bad or (bad==[''] and str(p)=='/nix/store/wkiqjb44n5k0fwhmc523dvflspd47q3b-gcc-15.3.0-libgcc/lib/libgcc_s.so.1'),(str(p),bad)
 if str(p).startswith('/nix/store/'):mounts.add(carrier(p))
 search=[p.parent]
 for t in m['dynamic']:
  if t['tag'] in ('RUNPATH','RPATH'):
   for v in t['value'].split(':'):
    if not v:continue
    expanded=v.replace('$ORIGIN',str(p.parent));assert expanded.startswith('/nix/store/'),expanded;search.append(Path(expanded));mounts.add(carrier(expanded))
 if m['interpreter']:queue.append(Path(m['interpreter']));mounts.add(carrier(m['interpreter']))
 for t in m['dynamic']:
  if t['tag']!='NEEDED':continue
  hits=[s/t['value'] for s in search if (s/t['value']).is_file()]
  if not hits:hits=[Path(s)/'lib'/t['value'] for s in sorted(mounts) if (Path(s)/'lib'/t['value']).is_file()]
  assert hits,(p,t);resolved=hits[0].resolve();assert str(resolved).startswith('/nix/store/');mounts.add(carrier(resolved));edges.append({'consumer':str(p),'needed':t['value'],'resolved':str(resolved),'sha256':sha(resolved)});queue.append(resolved)
for name,item in records.items():mounts.add(item['output'])
for root in mounts:
 for exe in ('git','nix','cargo','rustc','cc','gcc','g++','c++','clang','clang++','tcc'):assert not (Path(root)/'bin'/exe).exists(),(root,exe)
save('native-library-resolution.json',{'elf_metadata':metadata,'resolved_NEEDED':edges,'only_store_resolution':True,'native_loader_not_yet_executed':True})
# All four secret variables absent from every Nix child. Last Nix call precedes runtime preparation/activation.
info=json.loads(run('runtime-mount-path-info',[nix,'path-info','--json',*sorted(mounts)]));save('runtime-mount-path-info.json',info)
compiler=Path('/nix/store/c3bpv60f5wgp4pmgsifaykxk214h6463-rust-stable-with-components-2026-09-03/bin/rustc');save('compiler-binding.json',{'path':str(compiler),'sha256':sha(compiler),'version':next(iter(records.values()))['compiler'],'root_derivations':'pure/*root-derivation.json','fingerprints':'artifacts/*/qualification.json'})
save('artifact-qualification.json',{'records':records,'source_bindings':b,'exporter':{'path':str(exporter),'sha256':sha(exporter)},'ordinary_old_native_byte_equal':records.get('ordinary',{}).get('native',{}).get('sha256')=='53c5efb6344b233dfba09612835c33940f4daf64b215462f4292335067c420a2','native_compositions_distinct':len(records)==2 and records['ordinary']['native']['sha256']!=records['private']['native']['sha256'],'pure_artifact_qualification_complete':True,'store_mounts':sorted(mounts),'extra_browser_mounts':json.loads((B/'egl-runtime.json').read_text())['extra_mounts']})
print(json.dumps({'qualified_wheels':list(records),'native_hashes':{n:r['native']['sha256'] for n,r in records.items()},'mounts':len(mounts),'ELFs':len(metadata),'NEEDED_edges':len(edges)}),flush=True)
