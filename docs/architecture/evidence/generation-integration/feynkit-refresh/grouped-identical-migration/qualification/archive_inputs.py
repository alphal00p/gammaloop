"""Archive exact approved packet bytes and declare compressed/external bindings."""
from pathlib import Path
import gzip,hashlib,json,shutil
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
E=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
F=W.parent/'generation-oracle-numerator-families';P=F/'target/grouped-identical-migration-preparation-20261031-r1'
R=F/'.local/review-grouped-identical-migration-astra-20261031-024225Z-1469f5e4'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
bindings=[]
def cp(source,dest,reference):
 source=Path(source);dest=E/dest;b=source.read_bytes();dest.parent.mkdir(parents=True,exist_ok=True);compressed=len(b)>1048576
 if compressed:dest=Path(str(dest)+'.gz');dest.write_bytes(gzip.compress(b,compresslevel=9,mtime=0));assert gzip.decompress(dest.read_bytes())==b
 else:dest.write_bytes(b)
 assert dest.stat().st_size<=1048576
 bindings.append({'original':str(source),'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'durable_copy':str(dest.relative_to(W)),'copy_sha256':sha(dest),'copy_bytes':dest.stat().st_size,'encoding':'gzip; exact original bytes' if compressed else 'identity','durable_reference':reference})
# Exact reviewed migration outputs, including source proof, failed attempts and seal.
seal=load(R/'seal.json')
for rel,r in seal['files'].items():
 source=R/rel;assert sha(source)==r['sha256'] and source.stat().st_size==r['bytes'];cp(source,'reviews/migration/'+rel,'generation-grouped-identical-migration-candidate-20261031 (new candidate, no approval)')
for rel in ['seal.json','seal.sha256']:cp(R/rel,'reviews/migration/'+rel,'new candidate; original review packet path retained')
for line in (P/'SHA256SUMS').read_text().splitlines():
 digest,rel=line.split('  ',1);assert sha(P/rel)==digest;cp(P/rel,'preparation/'+rel,'new candidate; original preparation exact path retained')
for rel in ['SHA256SUMS','seal.json','seal.sha256']:
 cp(P/rel,'preparation/'+rel,'new candidate; original preparation exact path retained')
for name in ['root-grouped-identical-literal-acceptance-20261031.json','root-host-stub-rustred-acceptance-20261031.json','root-grouped-none-acceptance-20261031.json','root-grouped-auto-acceptance-20261031.json','root-full-generator-pool1-acceptance-20261031.json','root-three-integration-host-acceptance-and-next-20261031.json']:
 cp(W/'.local/feynkit-refresh'/name,'roots/'+name,'exact retained root decision; no new self approval')
accepted=[]
for mode,root,ref in [('None',F/'.local/review-grouped-none-aggregation-astra-20261030-independent-b1','ad7ccbe7cee2373ed54fe16d17971491eb581e80'),('Auto',F/'.local/review-grouped-auto-aggregation-astra-20261031-independent-r1-1793413266607828149','113608b3917d903d0fcc730aa3760f114fbc3f04')]:
 for rel in ['review.json','seal.json','seal.sha256']:cp(root/rel,'reviews/accepted-'+mode+'/'+rel,ref+' and its retained independent review')
 rec=load(root/'review.json');accepted.append({'mode':mode,'candidate':ref,'review':str(root/'review.json'),'review_sha256':sha(root/'review.json'),'seal':str(root/'seal.json'),'seal_sha256':sha(root/'seal.json'),'source_conditions':rec.get('review_composition',rec.get('premise_conjunction')),'full_review_copied':True,'fresh_replay_credit':False})
(E/'input-bindings.json').write_text(json.dumps({'bindings':bindings,'accepted_modes':accepted,'source_input_roots_retained':True,'compressed_copies_are_exact_uncompressed_bytes':True,'tracking_limit_bytes':1048576},indent=2)+'\n')
print('Archived',len(bindings),'exact approved input bindings, including full new source/interface proof and separate accepted mode reviews')
