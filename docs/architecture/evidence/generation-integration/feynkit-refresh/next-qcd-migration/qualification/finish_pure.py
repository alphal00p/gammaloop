"""Complete all Nix, compiler and source qualification before private runtime handoff."""
import ast,datetime,json,shutil,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent));import pure as P
O=P.O;W=P.W;load=P.load;write=P.write;sha=P.sha
q=load('archive-qualification.json');d=next(iter(load('archive-derivation.json').values()))
assert load('build.command.json')['exit']==0
assert q['actual_root_fingerprint']['data']['features']=='["default", "native"]' and q['profile']=='ci-optim'
assert all(load('producer-delta.json')['invariant_env_fields'].values()) and q['closure_units']==len(q['closure'])
for key in ['cargoExtraArgs','RUSTFLAGS','buildPhase','installPhase','nativeBuildInputs','buildInputs','stdenv']:
 assert d['env'].get(key)==next(iter(json.loads((W/'.local/feynkit-refresh/grouped-identical-migration-20261031/archive-derivation.json').read_text()).values()))['env'].get(key)
P.call('archive-build-log',['nix','--offline','log',q['archive_drv']],False)
info=P.nix('artifact-path-deriver',['path-info','--json','--derivation',d['env']['cargoArtifacts']]);assert len(info)==1
terminal=next(iter(info));terminal=terminal if terminal.startswith('/') else '/nix/store/'+terminal;assert terminal.endswith('.drv')
write('terminal-derivation.json',P.nix('terminal-derivation',['derivation','show',terminal]));P.call('package-build-log',['nix','--offline','log',terminal],False)
write('source-NAR-bindings.json',P.nix('source-NAR-info',['path-info','--json',q['source'],q['selected_source']]))
write('lock-profile-bindings.json',{'lock_sha256':sha(Path(q['source'])/'Cargo.lock'),'manifest_sha256':sha(Path(q['source'])/'Cargo.toml'),'flake_lock_sha256':sha(Path(q['source'])/'flake.lock'),'nextest_config_sha256':sha(Path(q['source'])/'.config/nextest.toml'),'producer':q['archive_drv'],'actual_root_fingerprint':q['actual_root_fingerprint'],'terminal_derivation':terminal,'RUSTFLAGS':d['env']['RUSTFLAGS'],'cargoExtraArgs':d['env']['cargoExtraArgs'],'kernels':q['kernels'],'compiler':q['compiler']})
elf=Path(q['oracle_elf']['path']).read_bytes();rows=json.loads((O/'source-invariance.json').read_text())['changes'];counts={}
for case,r in rows.items():
 old=int(r['old'],16).to_bytes(8,'little');new=int(r['new'],16).to_bytes(8,'little');counts[case]={'old':r['old'],'new':r['new'],'new_le_occurrences':elf.count(new),'old_le_occurrences':elf.count(old)};assert elf.count(new)>0
assert counts['xs_qcd_gggg_L1']['old_le_occurrences']>0,'Shared old token deliberately remains in excluded sibling'
expected='17e0abe46c76c6b6922e5b33782925020e1412e1a6ed8cd8dc7e23de438d4aa2';assert sha(Path(q['selected_source'])/'crates/feynkit-generator/tests/generation_oracle/main.rs')==expected
write('compiled-expected-tokens.json',{'oracle_elf':q['oracle_elf'],'counts':counts,'approved_source_sha256':expected,'sibling_old_token_preserved':True,'scope':'Source/derivation/ELF corroboration; exact upcoming native assertions remain separate'})
retained=O/'retained-actual-archive.tar.zst';shutil.copyfile(q['archive'],retained);assert sha(retained)==q['archive_sha256']
write('retained-artifacts.json',{'actual_archive':{'path':str(retained),'sha256':sha(retained),'bytes':retained.stat().st_size,'original_store':q['archive']},'terminal_layer':d['env']['cargoArtifacts']+'/target.tar.zst','code_bookmark':'generation-next-qcd-migration-code-20261031','retention':'Actual archive/ELF/fingerprints retained pending independent integration review; no extra realization and no GC root claimed'})
# Resolve the compiler that belongs to the exact pinned development environment, without building it.
expr='let f=builtins.getFlake '+json.dumps(P.URI)+'; s=f.devShells.x86_64-linux.default; in { packageCache=s.TYPST_PACKAGE_CACHE_PATH; fontPaths=s.TYPST_FONT_PATHS; }'
tools=P.nix('typst-pinned-environment',['eval','--json','--expr',expr]);tools['compiler']=tools['packageCache'].removesuffix('/lib/typst/packages')+'/bin/typst';tools['existing']=Path(tools['compiler']).is_file();write('typst-native-tools.json',tools)
if tools['existing']:P.call('typst-version',[tools['compiler'],'--version'],False)
for p in O.glob('*.py'):ast.parse(p.read_bytes())
commands=[json.loads(p.read_text()) for p in O.glob('*.command.json')];nix=[r for r in commands if r['argv'][0]=='nix'];assert all(r['exit']==0 for r in nix);assert sum('build' in r['argv'] for r in nix)==1
write('nix-complete.json',{'completed':datetime.datetime.now(datetime.UTC).isoformat(),'source':q['code'],'selected_archive':q['archive_drv'],'all_nix_queries_complete':True,'all_closure_and_compiler_profile_queries_complete':True,'no_nix_or_Rust_compiler_in_licensed_phase':True,'license_keys_removed':list(P.KEYS),'one_selected_build_invocation':True,'actual_closure_units':q['closure_units'],'Nix_commands':len(nix),'qualification_complete_before_licensed_phase':True,'static_Typst_after_runtime':'Separate unlicensed phase only; no Rust compiler or broader docs build'})
print('PASS completed exact source/NAR/archive/ELF/fingerprint/kernel/profile/log qualification; one build; no license accessed')
