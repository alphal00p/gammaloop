from pathlib import Path
import json,hashlib,subprocess
root=Path.cwd();e=root/'docs/architecture/evidence/generation-integration/mark-geometry/round2-correction';old='User explicitly allowed narrow exact Kurbo vendor patch and necessary root/host selection entries.';new="The root coordinator authorized the narrow exact Kurbo 0.13.1 vendor patch and necessary root/notebook-host patch-selection entries within the accepted correctness task, preserving all unrelated dependency versions, pins and locks."
p=e/'receipt.json';d=json.loads(p.read_text());assert d['dependency_exception']['authorization']==old;d['dependency_exception']['authorization']=new;p.write_text(json.dumps(d,indent=2)+'\n')
p=e/'package-evidence.py';s=p.read_text();assert s.count(old)==1;p.write_text(s.replace(old,new))
p=root/'docs/architecture/generation-mark-geometry.typ';s=p.read_text();oldnote='both the root and notebook-host manifests. The user authorized only the necessary\npatch-selection lock changes; dependency versions, kernel and other locks stay\nfixed.';newnote='both the root and notebook-host manifests. Within the accepted correctness task,\nthe root coordinator authorized the narrow exact Kurbo 0.13.1 vendor patch and\nnecessary root/notebook-host patch-selection entries; all unrelated dependency\nversions, pins and locks stay fixed.';assert oldnote in s;p.write_text(s.replace(oldnote,newnote))
p=e/'evidence-inventory.json';rows=json.loads(p.read_text());changed=[]
for r in rows:
 f=e/r['path'];sha=hashlib.sha256(f.read_bytes()).hexdigest()
 if sha!=r['sha256']:
  assert r['path'] in ['receipt.json','package-evidence.py'],r['path'];changed.append(r['path']);r.update(sha256=sha,size=f.stat().st_size)
assert sorted(changed)==['package-evidence.py','receipt.json'];p.write_text(json.dumps(rows,indent=2)+'\n')
