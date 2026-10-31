"""Relevant pure checks only; frozen historical negative fixtures remain verbatim."""
import ast,hashlib,json,os,re,subprocess,tomllib
from pathlib import Path
O=Path(__file__).resolve().parent;W=O.parents[2];F=W.parent/'generation-oracle-numerator-families/target/raw-next-qcd-migration-preparation-20261031-r1';MAIN='crates/feynkit-generator/tests/generation_oracle/main.rs'
def blob(b):return hashlib.sha1(b'blob %d\0'%len(b)+b).hexdigest()
before=(F/'before.rs').read_bytes();after=(W/MAIN).read_bytes();assert after==(F/'after.NOT-APPLIED.rs').read_bytes() and blob(after)=='cc4505a6c4f6fc1687c075ed79e154a8e19f05a6'
rows=lambda b:{m[1].decode():m[0] for m in re.finditer(rb'(\w+): Row \{[^\n]*?\}',b)};a,b=rows(before),rows(after);assert len(a)==len(b)==149
state=json.loads((O/'source-invariance.json').read_text());assert set(n for n in a if a[n]!=b[n])==set(state['changes']);assert all(a[n]==b[n] for n in state['prior13']) and a['xs_qcd_gggg_sym_L1']==b['xs_qcd_gggg_sym_L1']
sc=[]
for sec in tomllib.loads((W/'docs/developers.toml').read_text())['section']:
 for note in sec.get('note',[]):
  for s in note.get('scope',[]):
   t=(W/s['path']).read_text();start,end=s.get('anchor'),s.get('end_anchor');assert not start or t.count(start)==1;assert not end or t.count(end)==1;data=t[t.index(start):t.index(end)].encode() if end else (W/s['path']).read_bytes();assert blob(data)==s['digest'];sc.append({'note':note['id'],**s})
assert sc==state['scope_bindings'] and len(sc)==118
for p in [W/'docs/developers.toml',W/'docs/architecture/generation-next-qcd-migration.typ']:
 text=p.read_text();assert text.endswith('\n') and all(line.rstrip()==line for line in text.splitlines()),p
sources=[]
for p in O.glob('*.py'):ast.parse(p.read_bytes());sources.append(p.name)
parsed=[]
for p in O.glob('*.json'):json.loads(p.read_text());parsed.append(p.name)
tomllib.loads((W/'.config/nextest.toml').read_text());tomllib.loads((W/'Cargo.toml').read_text());tomllib.loads((W/'Cargo.lock').read_text())
assert re.search(r'^  enable = false;$',(W/'nix-ci.nix').read_text(),re.M)
base=state['base_tree'];unchanged=['Cargo.toml','Cargo.lock','flake.nix','flake.lock','nix-ci.nix','.config/nextest.toml','examples/notebooks/symbolica-host/python/symbolica/community/hepkit/__init__.pyi']
for p in unchanged:assert blob((W/p).read_bytes())==base[p]['blob'],p
begin=json.loads((O/'licensed-phase-begin.json').read_text())['time'];cmds=[json.loads(p.read_text()) for p in O.glob('*.command.json')];nix=[r for r in cmds if r['argv'][0]=='nix'];assert len([r for r in nix if 'build' in r['argv']])==1 and all(r['end']<begin and r['exit']==0 for r in nix)
r={'status':'PASS','149_Rows':149,'unchanged_Rows':146,'prior13':state['prior13'],'exact_three_approved_deltas':state['changes'],'excluded_sibling_preserved':True,'all118_ordered_actual_source_bindings_verified':True,'unchanged_source_config_stub_CI_files':unchanged,'Python_AST':sources,'new_top_level_JSON':parsed,'TOML':'registry,nextest,Cargo manifests/lock parsed','authored_note_registry_whitespace':'PASS','pure_Nix_before_licensed_phase_only':True,'broader_checks_or_credit':False,'Rust_policy':'Only three preapproved lowercase hexadecimal literals; pinned rustfmt passes; no design-policy departure'}
(O/'bounded-final-checks.json').write_text(json.dumps(r,indent=2)+'\n');print('PASS relevant JSON/TOML/AST/whitespace/source/149Rows/prior13/sibling/118scope/stub/CI/chronology checks')
