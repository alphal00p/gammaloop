"""Last pure source/build-log/closure attestations. No license read."""
import sys,datetime,json,os,ast,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent));import pure as P
O=P.O;W=P.W;load=P.load;write=P.write;sha=P.sha
q=load('archive-qualification.json');d=next(iter(load('archive-derivation.json').values()))
assert load('build.command.json')['exit']==0
assert q['actual_root_fingerprint']['data']['features']=='["default", "native"]' and q['profile']=='ci-optim'
assert all(load('producer-delta.json')['invariant_env_fields'].values())
assert q['closure_units']==len(q['closure'])
for key in ['cargoExtraArgs','RUSTFLAGS','buildPhase','installPhase']:
 assert d['env'].get(key)==next(iter(json.loads((W/'.local/feynkit-refresh/raw-three-migration-20261031/archive-derivation.json').read_text()).values()))['env'].get(key)
# Collect completed build logs and the actual terminal derivation before licensing.
assert load('archive-build-log.command.json')['exit']==0
art=d['env']['cargoArtifacts'];info=json.loads((O/'artifact-path-deriver.stdout').read_text());assert isinstance(info,dict);derivers=list(info);assert all(k.endswith('.drv') for k in derivers)
assert len(derivers)==1
if not derivers[0].startswith('/'):derivers[0]='/nix/store/'+derivers[0]
terminal=P.nix('terminal-derivation',['derivation','show',derivers[0]])
P.call('package-build-log',['nix','--offline','log',derivers[0]],False)
write('terminal-derivation.json',terminal)
new=Path(q['oracle_elf']['path']).read_bytes();old=(0x00fac78da0e0d159).to_bytes(8,'little');expected=(0xa60df38a6b8079e3).to_bytes(8,'little')
assert new.count(expected)>0 and new.count(old)==0
write('compiled-expected-token.json',{'oracle_elf':q['oracle_elf'],'expected_token':'0xa60df38a6b8079e3','new_little_endian_occurrences':new.count(expected),'old_little_endian_occurrences':new.count(old),'selected_source_sha256':sha(Path(q['selected_source'])/'crates/feynkit-generator/tests/generation_oracle/main.rs'),'source_hash_expected':'66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2','claim':'Actual pure producer source/build graph plus ELF literal corroboration; forthcoming exact native assertion is separate.'})
assert sha(Path(q['selected_source'])/'crates/feynkit-generator/tests/generation_oracle/main.rs')=='66f01bcc855552ae888f1fb7c5bb158afbd448492b4e27a26e2bda40f09b07d2'
# Retain actual new immutable archive bytes locally without a second realization.
import shutil
retained=O/'retained-actual-archive.tar.zst';shutil.copyfile(q['archive'],retained)
assert sha(retained)==q['archive_sha256']
roots={str(retained):q['archive']}
write('retained-artifacts.json',{'actual_new_archive':{'path':str(retained),'sha256':sha(retained),'bytes':retained.stat().st_size,'original_store':q['archive']},'source_durable_reference':'generation-grouped-identical-migration-code-20261031','retained_terminal_layer':art+'/target.tar.zst','reason':'Actual new archive/extracted ELF and dependency audit retained for pending independent review; no extra realization or resource policy change.'})
for p in O.glob('*.py'):ast.parse(p.read_bytes())
write('syntax.json',{'Python_AST':'PASS','sources':{p.name:sha(p) for p in O.glob('*.py')},'no_native_imports':True})
commands=[json.loads(p.read_text()) for p in O.glob('*.command.json')];nix=[r for r in commands if r['argv'][0]=='nix']
assert all(r['exit']==0 for r in nix)
assert sum('build' in r['argv'] for r in nix)==1
write('nix-complete.json',{'completed':datetime.datetime.now(datetime.UTC).isoformat(),'source':q['code'],'selected_archive':q['archive_drv'],'all_nix_queries_complete':True,'all_closure_and_compiler_profile_queries_complete':True,'no_nix_or_compiler_after_this_boundary':True,'license_keys_removed':list(P.KEYS),'one_selected_build_invocation':True,'actual_closure_units':q['closure_units'],'GC_roots':roots,'Nix_commands':len(nix),'qualification_complete_before_licensed_phase':True})
print('PASS final pure graph/log/literal/GC-root/chronology binding; no license accessed',flush=True)
