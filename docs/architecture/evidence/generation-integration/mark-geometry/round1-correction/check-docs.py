from pathlib import Path
import runpy,json,tomllib,subprocess,hashlib
root=Path.cwd();out=root/'target/mark-geometry-corrected';owner=root/'docs/architecture/evidence/generation-integration/wp3/oracle-registry-audit.py';check=runpy.run_path(str(owner))['check']
blob=lambda text:subprocess.check_output(['git','hash-object','--stdin'],input=text.encode()).decode().strip()
registry=tomllib.loads((root/'docs/developers.toml').read_text());rows=[];titles=[]
for section in registry['section']:
 for note in section.get('note',[]):
  source=(root/note['source']).read_text();assert [l[2:] for l in source.splitlines() if l.startswith('= ')]==[note['title']],note['id'];titles.append(note['id'])
  if note['lifecycle']=='current':assert note.get('scope'),note['id']
  for scope in note.get('scope',[]):
   source=(root/scope['path']).read_text();r=check(scope,source);assert r['status']=='valid',(note['id'],scope,r)
   payload=source[source.index(scope['anchor']):source.index(scope['end_anchor'])] if 'end_anchor' in scope else source
   assert blob(payload)==scope['digest'];rows.append({'note':note['id'],**scope,**r})
text='outside\nBEGIN\ninside\nEND\nafter\n';full={'anchor':'BEGIN','digest':blob(text)};part={'anchor':'BEGIN','end_anchor':'END','digest':blob('BEGIN\ninside\n')}
cases=[(full,text,'valid'),({**full,'anchor':'missing'},text,'missing-or-ambiguous-anchor'),(full,text+'BEGIN','missing-or-ambiguous-anchor'),({**full,'anchor':''},text,'missing-or-ambiguous-anchor'),(part,text,'valid'),({**part,'end_anchor':''},text,'missing-or-ambiguous-end'),({**part,'end_anchor':'absent'},text,'missing-or-ambiguous-end'),(part,text+'END','missing-or-ambiguous-end'),({**part,'anchor':'END','end_anchor':'BEGIN'},text,'invalid-range'),({'anchor':'BEGIN','end_anchor':'GIN','digest':'x'},text,'invalid-range'),(full,text+'changed','stale-digest'),(full,text+'#[cfg(test)]','stale-digest'),(part,'outside changed'+text,'valid'),(part,text.replace('inside','changed'),'stale-digest'),({'anchor':'α','digest':blob('αβγ')},'αβγ','valid'),(full,None,'missing-source')]
for s,t,expected in cases:assert check(s,t)['status']==expected,(s,t,expected)
base=subprocess.check_output(['jj','file','show','-r','100da206','crates/kurvst/typst/docs/manual.typ'],text=True);current=(root/'crates/kurvst/typst/docs/manual.typ').read_text();assert current.count('table.header(')==base.count('table.header(')
(out/'registry.json').write_text(json.dumps({'source_owner':str(owner.relative_to(root)),'owner_sha256':hashlib.sha256(owner.read_bytes()).hexdigest(),'valid_scopes':len(rows),'valid_titles':len(titles),'negative_controls':len(cases),'manual_semantic_table_headers_preserved':True,'scopes':rows},indent=2)+'\n')
print(f'{len(rows)} scopes, {len(titles)} titles, {len(cases)} canonical controls pass; manual table headers preserved')
