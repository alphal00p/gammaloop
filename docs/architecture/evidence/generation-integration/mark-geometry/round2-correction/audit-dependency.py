from pathlib import Path
import hashlib,tarfile,json,subprocess,tomllib,difflib
out=Path('target/mark-geometry-round3');src=Path('/home/lcnbr/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/kurbo-0.13.1');crate=Path('/home/lcnbr/.cargo/registry/cache/index.crates.io-1949cf8c6b5b557f/kurbo-0.13.1.crate');sha=lambda b:hashlib.sha256(b).hexdigest();assert sha(crate.read_bytes())=='4b60dfc32f652b926df6192e55525b16d186c69d47876c3ead4da5cc9f8450e2';rows=[]
with tarfile.open(crate) as t:
 for m in t.getmembers():
  if not m.isfile():continue
  rel=Path(m.name).relative_to('kurbo-0.13.1');data=t.extractfile(m).read();assert (src/rel).read_bytes()==data;local=Path('vendor/kurbo')/rel;assert local.is_file();rows.append({'path':str(rel),'upstream_sha256':sha(data),'local_sha256':sha(local.read_bytes()),'changed':local.read_bytes()!=data})
assert [r['path'] for r in rows if r['changed']]==['src/stroke.rs']
patch=''.join(difflib.unified_diff((src/'src/stroke.rs').read_text().splitlines(True),Path('vendor/kurbo/src/stroke.rs').read_text().splitlines(True),fromfile='upstream/kurbo-0.13.1/src/stroke.rs',tofile='local/kurbo-0.13.1/src/stroke.rs'));(out/'kurbo.patch').write_text(patch)
base='e426ff5fcbe1487c9bf2cf1ed5aaaa03cabf8722';locks={}
for path in ['Cargo.lock','examples/notebooks/symbolica-host/Cargo.lock','tydenso/Cargo.lock','uv.lock','flake.lock','nix/ufo-runtime/uv.lock','crates/linnest/layout/uv.lock']:
 old=subprocess.check_output(['jj','--ignore-working-copy','file','show','-r',base,path]);now=Path(path).read_bytes();item={'before_sha256':sha(old),'sha256':sha(now)}
 if path.endswith('Cargo.lock'):
  a=tomllib.loads(old.decode())['package'];b=tomllib.loads(now.decode())['package'];key=lambda p:(p['name'],p['version']);aa={key(p):p for p in a};bb={key(p):p for p in b};assert aa.keys()==bb.keys();diff=[{'before':aa[k],'after':bb[k]} for k in aa if aa[k]!=bb[k]];item['package_changes']=diff
  if path.startswith('tydenso'):assert not diff
  else:
   assert len(diff)==1 and diff[0]['before']['name']=='kurbo';before=diff[0]['before'];after=diff[0]['after'];assert {k:v for k,v in before.items() if k not in ['source','checksum']}==after
 else:assert old==now
 locks[path]=item
(out/'dependency-audit.json').write_text(json.dumps({'crate_archive':str(crate),'archive_sha256':sha(crate.read_bytes()),'upstream_git':'838b69ef35f3381cefc007d86aa8364f5a6dc673','version':'0.13.1','files':rows,'local_patch_sha256':sha(patch.encode()),'locks':locks,'licenses':'Both exact upstream license files retained; existing shipped package license/notice owners unchanged. No broader legal clearance claimed.'},indent=2)+'\n')
print(f'{len(rows)} published files verified; stroke.rs sole local source delta; exactly two Kurbo selection entries changed')
