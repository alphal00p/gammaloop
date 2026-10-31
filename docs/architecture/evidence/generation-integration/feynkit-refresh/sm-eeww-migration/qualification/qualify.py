"""Actual new source/archive/units; reuse old metadata only for identical inputs."""
import importlib.util
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import re
import subprocess
import tarfile

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('pure',HERE/'pure.py');P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
O,W,OLD=P.O,P.W,P.OLD
sha,write,load,nix=P.sha,P.write,P.load,P.nix
assert load('build.command.json')['exit']==0
meta=load('flake.json'); full=Path(meta['source']); selected=(O/'selected-drv.txt').read_text().strip()
queue=[selected];chain={};sources={}
while queue:
    batch=sorted(set(queue)-set(chain));queue=[]
    if not batch:break
    ds=nix('chain-'+str(len(chain)),['derivation','show',*batch])
    for drv,d in ds.items():
        chain[drv]=d;e=d.get('env',{});src=e.get('src')
        if src and Path(src).is_dir() and src not in sources:
            files=[]
            for f in sorted(Path(src).rglob('*')):
                if not f.is_file() and not f.is_symlink():continue
                rel=str(f.relative_to(src));data=os.fsencode(os.readlink(f)) if f.is_symlink() else f.read_bytes();base=full/rel
                eq=base.exists() and data==(os.fsencode(os.readlink(base)) if base.is_symlink() else base.read_bytes())
                files.append({'path':rel,'sha256':hashlib.sha256(data).hexdigest(),'full_tree_byte_equal':eq})
            sources[src]={'used_by':drv,'files':files}
        for k in d.get('inputDrvs',d.get('inputs',{}).get('drvs',{})):
            p=k if k.startswith('/') else '/nix/store/'+k
            if re.search(r'-(gammaloop-|vendor-cargo|rust-stable-with-components)',p):queue.append(p)
oldchain=json.loads((OLD/'producer-chain.json').read_text())
shared=sorted(set(chain)&set(oldchain['derivations']))
write('producer-chain.json',{'derivations':chain,'filtered_sources':sources,'prior_chain_sha256':sha(OLD/'producer-chain.json'),'identical_derivation_ids_reused':shared,'new_derivation_ids':sorted(set(chain)-set(oldchain['derivations']))})
# Exact entire fetched tree, including symlinks and path set.
raw=P.call('git-source-tree',['git','--git-dir='+str(W.parent/'main/.git'),'ls-tree','-rz','--full-tree',P.REV],False)
expected={};bad=[]
for row in raw.split(b'\0'):
    if not row:continue
    header,name=row.split(b'\t',1);mode,kind,oid=header.decode().split();name=os.fsdecode(name);f=full/name
    assert kind=='blob'
    data=os.fsencode(os.readlink(f)) if mode=='120000' else f.read_bytes();actual=hashlib.sha1(b'blob %d\0'%len(data)+data).hexdigest();expected[name]=oid
    if actual!=oid:bad.append(name)
actual=set()
for d,ds,fs in os.walk(full):
    for n in fs+[n for n in ds if (Path(d)/n).is_symlink()]:actual.add(str((Path(d)/n).relative_to(full)))
assert not bad and actual==set(expected)
write('fetched-tree-check.json',{'revision':P.REV,'source':str(full),'tracked_paths':len(expected),'blobs':expected,'mismatches':bad,'extra':sorted(actual-set(expected)),'missing':sorted(set(expected)-actual)})
src=next(iter(load('archive-derivation.json').values()))['env']['src']; oldsrc=next(iter(json.loads((OLD/'archive-derivation.json').read_text()).values()))['env']['src']
a={r['path']:r['sha256'] for r in sources[src]['files']};b={r['path']:r['sha256'] for r in oldchain['filtered_sources'][oldsrc]['files']}
assert a.keys()==b.keys();delta=[p for p in a if a[p]!=b[p]];assert delta==['crates/feynkit-generator/tests/generation_oracle/main.rs'];assert all(r['full_tree_byte_equal'] for r in sources[src]['files'])
write('selected-source-delta.json',{'old':oldsrc,'new':src,'paths':len(a),'changed_paths':delta,'new_files':a})
build=json.loads((O/'build.stdout').read_text());assert len(build)==1;store=Path(build[0]['outputs']['out']);archive=store/'gammaloop-nextest-feynkit-feynkit-generator.tar.zst'
dest=O/'archive-extracted';dest.mkdir(exist_ok=False)
listing=P.call('archive-list',['tar','--zstd','-tf',str(archive)],False).decode().splitlines()
assert all(not PurePosixPath(p).is_absolute() and '..' not in PurePosixPath(p).parts for p in listing)
P.call('archive-extract',['tar','--zstd','--no-same-owner','--no-same-permissions','-xf',str(archive),'-C',str(dest)],False)
oldq=json.loads((OLD/'archive-qualification.json').read_text());art=Path(next(iter(load('archive-derivation.json').values()))['env']['cargoArtifacts']);layer=art/'target.tar.zst'
layers=[];p=layer
while p.exists():
    layers.append({'path':str(p.resolve()),'sha256':sha(p)})
    prev=Path(str(p)+'.prev')
    if not prev.exists():break
    p=prev.resolve()
old_layers=oldq['compiled_layers'];assert len(layers)==len(old_layers)
assert [l['sha256'] for l in layers[1:]]==[l['sha256'] for l in reversed(old_layers[:-1])],'dependency layers changed; do not reuse unqualified units'
# Read new terminal layer fingerprints only; never copy a donor target or relink.
fpdir=O/'new-fingerprints';fpdir.mkdir(exist_ok=False)
flist=P.call('artifact-layer-list',['tar','--zstd','-tf',str(layer)],False).decode().splitlines()
fps=[s for s in flist if '/.fingerprint/' in s and not s.endswith('/')]
assert fps and all(not PurePosixPath(p).is_absolute() and '..' not in PurePosixPath(p).parts for p in fps)
P.call('artifact-fingerprints',['tar','--zstd','--no-same-owner','--no-same-permissions','-xf',str(layer),'-C',str(fpdir),*fps],False)
units={k:list(v) for k,v in oldq['all_fingerprint_index'].items()};new=[]
for f in fpdir.rglob('*.json'):
    marker=f.with_suffix('')
    if not marker.is_file():continue
    try:k=marker.read_text().strip()
    except UnicodeError:continue
    if not re.fullmatch('[0-9a-f]{16}',k):continue
    r={'unit':f.parent.name,'kind':f.stem,'path':str(f),'sha256':sha(f),'stored_fingerprint':k,'data':json.loads(f.read_text()),'source_layer':str(layer)}
    units.setdefault(k,[]).append(r);new.append(r)
roots=[r for r in new if r['kind']=='test-integration-test-generation_oracle'];assert len(roots)==1;root=roots[0];seen=set();closure=[];todo=[root]
while todo:
    r=todo.pop();k=r['stored_fingerprint']
    if k in seen:continue
    seen.add(k);closure.append(r)
    for dep in r['data']['deps']:
        key=dep[3].to_bytes(8,'little').hex();assert key in units,(dep,key)
        choices=units[key];assert all(x['data']==choices[0]['data'] for x in choices);todo.append(choices[-1])
metadata=json.loads((dest/'target/nextest/cargo-metadata.json').read_text());kernels=[{k:p[k] for k in ['name','version','source','manifest_path']} for p in metadata['packages'] if p['name'] in ['symbolica','numerica','graphica']]
assert len(kernels)==3 and all(p['source']=='git+https://github.com/symbolica-dev/symbolica?branch=community#8eb7c3729d6abefc916decc1bdf0c06536f8af80' for p in kernels)
compiler=json.loads((OLD/'compiler.json').read_text());P.call('compiler',[compiler['rustc']['path'],'-vV'],False)
assert (O/'compiler.stdout').read_text().strip()==compiler['rustc']['version']
elfs=[]
for f in dest.rglob('*'):
    if f.is_file():
        with f.open('rb') as stream:magic=stream.read(4)
        if magic==b'\x7fELF':elfs.append({'path':str(f),'sha256':sha(f),'size':f.stat().st_size})
oracle=[p for p in elfs if Path(p['path']).name.startswith('generation_oracle-')];assert len(oracle)==1;assert oracle[0]['sha256']!=next(p['sha256'] for p in oldq['elfs'] if Path(p['path']).name.startswith('generation_oracle-'))
nar=nix('archive-path-info',['path-info','--json',str(store)])
artnar=nix('artifact-path-info',['path-info','--json',str(art)])
nextest=Path('/nix/store/kwm3g9vh36d9zd171v0nz64fyan1hvfz-cargo-nextest-0.9.143/bin/cargo-nextest');python=Path('/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin/python3');bwrap=Path('/nix/store/dk9qhjgg469lv6mriys7v4c59igarmvx-bubblewrap-0.11.1/bin/bwrap')
nix('runtime-tool-path-info',['path-info','--json',str(nextest.parents[1]),str(python.parents[1]),str(bwrap.parents[1])])
write('archive-qualification.json',{'code':P.REV,'uri':P.URI,'source':str(full),'source_paths':len(expected),'archive_drv':selected,'archive':str(archive),'archive_sha256':sha(archive),'archive_store':str(store),'archive_nar':nar,'artifact_nar':artnar,'selected_source':src,'selected_source_delta':delta,'compiled_layers':layers,'unchanged_dependency_layers':True,'actual_root_fingerprint':root,'closure_units':len(closure),'closure':closure,'new_terminal_layer_fingerprints':new,'kernels':kernels,'compiler':compiler,'profile':'ci-optim','elfs':elfs,'oracle_elf':oracle[0],'runtime_tools':{'nextest':str(nextest),'python':str(python),'bwrap':str(bwrap)},'old_452_unit_metadata_reused':'Only unchanged dependency layers with exact archive bytes and recursive dependency fingerprints. New terminal root and archive ELFs are actual new outputs.','runtime_credit':'NONE YET'})
print('qualified',len(expected),'source paths;',len(closure),'actual closure units; new ELF',oracle[0]['sha256'],flush=True)
