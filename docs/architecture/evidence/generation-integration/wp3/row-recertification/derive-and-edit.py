from pathlib import Path
import json,struct,hashlib,re,shutil
O=Path('target/wp3-units19-21-recertification');F=Path('target/wp3-units19-21/frozen-row-diagnosis');M=(1<<64)-1
# Exact native FNV64 encoding copied from the reviewed Census::digest_json/Digest owner.
def fnv(b,h=0xcbf29ce484222325):
 for byte in b:h=((h^byte)*0x100000001b3)&M
 return h
def enc(v):return json.dumps(v,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()
def digest(ds):
 result=[]
 for d in ds:
  p=d['payload'];h=0xcbf29ce484222325
  for k in ['numerator','numerator_prefactor','projector']:h=fnv(enc(p[k]),h)
  for k in ['num_interface','numerator_prefactor_interface','projector_interface']:
   if k in p:h=fnv(enc(p[k]),h)
  for frag in p['vertices']+[e[1] for e in p['edges']]:
   h=fnv(enc(frag['numerator']),h)
   if 'num_interface' in frag:h=fnv(enc(frag['num_interface']),h)
  result.append(fnv(struct.pack('<QQ',fnv(enc(d['key'])),h)))
 return fnv(b''.join(struct.pack('<Q',n) for n in sorted(result)))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
proof=json.loads((F/'producer-row-reproduction.json').read_text());rows=[];p=Path('crates/feynkit-generator/tests/generation_oracle/main.rs');before=p.read_text();after=before
for r in proof:
 name=r['case'];src=F/'current-decoded/captures'/name/'grouped.json';old=F/'baseline-decoded/captures'/name/'grouped.json';value=digest(json.loads(src.read_text())['diagrams']);prior=digest(json.loads(old.read_text())['diagrams']);assert(value,prior)==(r['current_native_digest'],r['baseline_native_digest'])
 pattern=r'(^    fast '+re.escape(name)+r': Row \{[^\n]*?numerators: )(0x[0-9a-f]+)';m=re.search(pattern,after,re.M);assert m and int(m[2],16)==prior
 after,n=re.subn(pattern,lambda m:m[1]+f'0x{value:016x}',after,flags=re.M);assert n==1
 rows.append({'case':name,'old_integer':prior,'old_hex':f'0x{prior:016x}','new_integer':value,'new_hex':f'0x{value:016x}','current_native_capture':str(src),'current_sha256':sha(src),'baseline_native_capture':str(old),'baseline_sha256':sha(old),'hashing':'Every native interface retained; no head mapping, payload normalization or parameter substitution in the digest derivation.'})
assert before!=after
# Only the two authorized literal tokens may differ.
expected=before
for r in rows:expected=expected.replace(r['old_hex'],r['new_hex'],1)
assert expected==after
p.write_text(after);(O/'derivation.json').write_text(json.dumps({'input_code':'82c1eb5b40c4911e59dfc25c9ba1786e282e111c','frozen_candidate':'df55185059f31bea548a5de27c95828656fbdc52','before_file_sha256':hashlib.sha256(before.encode()).hexdigest(),'after_file_sha256':sha(p),'rows':rows,'precision':'Python exact arbitrary-precision integers; u64 wrapping is explicit. No JS-number or floating conversion.'},indent=2)+'\n')
print([(r['case'],r['old_hex'],r['new_hex']) for r in rows])
