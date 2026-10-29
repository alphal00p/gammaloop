from pathlib import Path
import json,hashlib,subprocess
r=Path.cwd();out=r/'target/wp3-source38';base='100da20603a083d031c31bd4e8c6970fc44dda2f';code='0c65bf7c519764460a15a922340d958a9dfa7429';prep='861c19e3e7afc067851b13f3e7548baa923286ee'
def git(*a):return subprocess.check_output(['git','--git-dir=../main/.git',*a])
def sha(b):return hashlib.sha256(b).hexdigest()
paths=git('diff','--name-only',base,code).decode().splitlines();rows=[]
for p in paths:
 try:b=git('show',code+':'+p);s=sha(b)
 except subprocess.CalledProcessError:s=None
 category='evidence' if '/evidence/' in p else 'generated-assets' if p.endswith(('.wasm','.svg','.png','.pdf')) else 'documentation' if p.startswith('docs/') or p.endswith('.typ') else 'runtime-and-build'
 rows.append({'path':p,'category':category,'sha256_at_frozen_qualified_code':s})
prep_paths=git('diff','--name-only',base,prep).decode().splitlines();reconciled=[]
for p in prep_paths:
 try:a=git('show',prep+':'+p)
 except subprocess.CalledProcessError:a=None
 try:b=git('show',code+':'+p)
 except subprocess.CalledProcessError:b=None
 if a!=b:reconciled.append(p)
assert reconciled==['docs/architecture/generation-mark-geometry.typ'],reconciled
refs=json.loads((r/'docs/architecture/evidence/generation-integration/wp3/source36/integrity.json').read_text())['protected_refs'];refs.update({'generation-integration-source36':base,'generation-integration-source38-code':'925f739799572d23f46a762c895ebb85e5bd88d6','generation-integration-source38-qualified':code})
refs.update({name:subprocess.check_output(['jj','--ignore-working-copy','log','-r',name,'--no-graph','-T','commit_id'],text=True).strip() for name in ['generation-mark-geometry','generation-mark-geometry-round1','generation-mark-geometry-round2','generation-mark-geometry-round3','generation-mark-geometry-round3-code','generation-mark-geometry-round3-build-input']})
for name,want in refs.items():
 actual=subprocess.check_output(['jj','--ignore-working-copy','log','-r',name,'--no-graph','-T','commit_id'],text=True).strip();assert actual==want,(name,actual,want)
for rev in ['925f739799572d23f46a762c895ebb85e5bd88d6',code]:assert subprocess.check_output(['jj','--ignore-working-copy','log','-r',rev,'--no-graph','-T','conflict'],text=True).strip()=='false'
report={'base':base,'frozen_qualified_code':code,'runtime_code':'925f739799572d23f46a762c895ebb85e5bd88d6','approved_preparation':prep,'preparation_changed_paths':len(prep_paths),'runtime_and_selected_source_bytes_match_preparation':True,'intentional_reconciliations':reconciled,'new_plan_and_main_review_bookkeeping':'Canonical plan and main38 review/archive files additionally update MAIN state; docs/developers.toml remains exact approved preparation bytes.','inventory':rows,'counts':{k:sum(x['category']==k for x in rows) for k in sorted({x['category'] for x in rows})},'protected_refs':refs,'NixCI':False,'kernel':'1a1829b1f3d0a88b53b44085b3b9f34ef77945b3','ignored_log_import':'Initial runtime commit925 omitted153 historical ignored .log paths. Explicit include-ignored tracking restored them byte-for-byte in0c65 evidence/metadata commit. No source or runtime difference; all original logs historical, not fresh passes.','self_hash_boundary':'Inventory is frozen0c65 Git blobs, excludes later validation receipts; no mutable receipt hashes used as code identity.'}
(out/'integrity.json').write_text(json.dumps(report,indent=2)+'\n');print(report['counts'],len(refs))
