"""Validate and seal the finite implementation candidate; no Nix or licensing."""
import ast
import datetime
import importlib.util
import re
import shutil
import sys
import tomllib
import xml.etree.ElementTree as ET
sys.path.insert(0, str(__file__).rsplit('/', 1)[0])
from runtime_common import *

assert not any(k in os.environ for k in ['SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED'])
BASE = '5ea4c93b79dd33094072d7717e912aa4b05d6689'
CODE = 'ffa76ec420642c03349ab756d2d7720ab93c0057'
E = W / 'docs/architecture/evidence/generation-integration/feynkit-refresh/raw-three-migration'
NOTE = W / 'docs/architecture/generation-raw-three-migration.typ'
ENV = dict(os.environ, GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')
spec = importlib.util.spec_from_file_location('guard', O / 'guard_apply.py')
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)
git, blob = G.git, G.blob

r = load(O / 'runtime-results.json')
assert r['code'] == Q['code'] == CODE
assert r['actual_executions'] == r['passes'] == 36
assert r['failures'] == r['retries'] == r['skipped_credit'] == 0
assert len(r['results']) == 3
assert load(O / 'selected-inventory.json')['count'] == 12
assert load(O / 'capture-joins.json')['same_physical_payloads']
assert load(O / 'capture-joins.json')['affected'] == 80
logs = []
boundary = datetime.datetime.fromisoformat(load(O / 'nix-complete.json')['completed'])
for pool in [1, 8, 32]:
    d = O / 'runs' / ('pool' + str(pool))
    c = load(d / 'command.json')
    assert c['exit'] == 0 and c['process_reaped'] and c['inputs_before'] == c['inputs_after']
    assert datetime.datetime.fromisoformat(c['start']) > boundary
    assert c['command_records_outside_writable_capsule']
    assert not d.is_relative_to(O / 'work')
    assert c['environment_before_exec']['FEYNKIT_ORACLE_THREADS'] == str(pool)
    assert c['environment_before_exec']['NEXTEST_TEST_THREADS'] == '16'
    assert not any(k in c['environment_before_exec'] for k in ['FEYNKIT_ORACLE_RECORD', 'SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED'])
    assert ('FEYNKIT_ORACLE_DUMP' in c['environment_before_exec']) == (pool == 1)
    assert '--retries' in c['argv'] and c['argv'][c['argv'].index('--retries') + 1] == '0'
    assert '--no-tests' in c['argv'] and c['argv'][c['argv'].index('--no-tests') + 1] == 'fail'
    assert c['argv'][c['argv'].index('--profile') + 1] == 'ci_gammaloop'
    assert sha(d / 'stdout') == c['stdout_sha256'] and sha(d / 'stderr') == c['stderr_sha256']
    tests = ET.parse(d / 'junit.xml').findall('.//testcase')
    assert len(tests) == 12 and {t.attrib['name'] for t in tests} == set(CASES)
    assert all(not any(t.find(k) is not None for k in ['failure', 'error', 'skipped', 'rerun']) for t in tests)
    textlog = (d / 'stderr').read_text()
    m = re.search(r'^\s*Summary \[[^\]]+\] 12 tests run: 12 passed, 171 skipped$', textlog, re.M)
    assert m and 'Starting 12 tests across 1 binary (171 tests and 3 binaries skipped)' in textlog
    passes = re.findall(r'^\s*PASS \[[^\]]+\] \(\s*\d+/12\) feynkit-generator::generation_oracle (\w+)$', textlog[:m.start()], re.M)
    assert len(passes) == 12 and set(passes) == set(CASES)
    assert not re.search(r'^\s*(RETRY|FAIL|TIMEOUT|FLAKY)\b', textlog, re.M)
    logs.append({'pool': pool, 'JUnit_selected': len(tests), 'unique_before_summary_PASS_records': len(passes), 'summary_selected_run': 12, 'summary_selected_pass': 12, 'selected_skips': 0, 'unselected_filtered_tests': 171, 'unselected_filtered_binaries': 3, 'final_output_repeats_not_new_executions': True})
write(O / 'runtime-log-audit.json', {'records': logs, 'actual_selected_executions': 36, 'retries': 0, 'selected_skips': 0, 'filtered_skip_credit': 0})

nix = []
for f in sorted(O.glob('*.command.json')):
    c = load(f)
    assert c['exit'] == 0
    assert sha(O / f.name.replace('.command.json', '.stdout')) == c['stdout_sha256']
    assert sha(O / f.name.replace('.command.json', '.stderr')) == c['stderr_sha256']
    if c['argv'][0] == 'nix':
        assert datetime.datetime.fromisoformat(c['end']) <= boundary
        assert set(c['license_keys_removed']) == {'SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED'}
        nix.append(c)
assert len([c for c in nix if c['argv'][1] == 'build']) == 1
write(O / 'pure-chronology.json', {'Nix_commands': len(nix), 'selected_build_calls': 1, 'all_exit_zero_before_boundary': True, 'boundary': boundary.isoformat(), 'Nix_calls_after_licensed_phase': 0, 'license_keys_removed': ['SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED']})

assert sha(Q['archive']) == Q['archive_sha256'] and sha(Q['oracle_elf']['path']) == Q['oracle_elf']['sha256']
prep = load(O / 'runtime-preparation.json')
for target, source in prep['file_mounts'].items():
    assert sha(source) == prep['file_hashes'][target]
for root in prep['stdlib_mounts']:
    assert {str(f.relative_to(root)): sha(f) for f in sorted(Path(root).rglob('*')) if f.is_file()} == prep['stdlib_files']
for p, oid in load(O / 'fetched-tree-check.json')['blobs'].items():
    f = Path(Q['source']) / p
    data = os.fsencode(os.readlink(f)) if f.is_symlink() else f.read_bytes()
    assert blob(data) == oid
for c in load(O / 'copied-evidence.json')['copies']:
    assert sha(c['source']) == sha(E / c['copy']) == c['sha256']
for j in load(O / 'capture-joins.json')['files']:
    assert sha(j['fresh']) == sha(j['approved_original']) == sha(j['approved_observer']) == sha(E / j['copy']) == j['fresh_sha256']

start = load(O / 'start.json')
path = 'crates/feynkit-generator/tests/generation_oracle/main.rs'
assert blob((W / path).read_bytes()) == '9ef9f1ba005f20149716b22346d45c5501b62e87'
assert git('show', CODE + ':' + path) == (W / path).read_bytes()
assert git('diff-tree', '--no-commit-id', '--name-only', '-r', BASE, CODE).decode().splitlines() == [path]
assert git('show', '-s', '--format=%P', CODE).decode().strip() == BASE
assert git('show', '-s', '--format=%an <%ae>', CODE).decode().strip() == 'Lucien Huber <im@lcnbr.ch>'
for p, record in start['tracked_paths'].items():
    if p in [path, 'docs/developers.toml']:
        continue
    f = W / p
    data = os.fsencode(os.readlink(f)) if f.is_symlink() else f.read_bytes()
    assert blob(data) == record['blob'], ('unapproved base drift', p)
oldreg = git('show', BASE + ':docs/developers.toml')
newreg = (W / 'docs/developers.toml').read_bytes()
oldsc = G.scopes(oldreg, lambda p: git('show', BASE + ':' + p))
newsc = G.scopes(newreg, lambda p: (W / p).read_bytes())
assert len(oldsc) == 118 and oldsc == newsc == load(O / 'source-check.json')['scope_bindings']
oldnotes = [n for s in tomllib.loads(oldreg.decode())['section'] for n in s.get('note', [])]
newnotes = [n for s in tomllib.loads(newreg.decode())['section'] for n in s.get('note', [])]
assert [n for n in newnotes if n['id'] != 'generation-raw-three-migration'] == oldnotes
note = next(n for n in newnotes if n['id'] == 'generation-raw-three-migration')
assert 'reviewed_at' not in note and 'review_ref' not in note and 'scope' not in note
assert NOTE.read_text().splitlines()[0] == '= ' + note['title']
assert sum(l.startswith('= ') for l in NOTE.read_text().splitlines()) == 1
assert re.search(r'^  enable = false;$', (W / 'nix-ci.nix').read_text(), re.M)
main = W.parent / 'main/.git'
assert {p: sha(main / p) if (main / p).exists() else None for p in start['main_git']} == start['main_git']
refs = lambda s: {r: h for h, r in (l.split() for l in s.splitlines())}
initialrefs, currentrefs = refs(start['refs']), refs(git('show-ref').decode())
assert all(currentrefs.get(k) == v for k, v in initialrefs.items())
integration = subprocess.check_output(['jj', '--ignore-working-copy', 'log', '--no-graph', '--no-pager', '-r', 'generation-integration', '-T', 'commit_id'], cwd=W, env=ENV, text=True).strip()
assert integration == 'fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
guard = {'base': BASE, 'code': CODE, 'all_base_files_except_exact_main_tokens_and_additive_registry_unchanged': True, 'all_host_API_browser_files_byte_pinned': True, 'all118_scopes_anchors_and_digests_equal': True, 'all_prior_registry_notes_equal': True, 'main_HEAD_index_packed_refs_equal': True, 'all_initial_git_refs_equal': True, 'generation_integration': integration, 'CI': False, 'all_archive_source_runtime_files_rehashed': True, 'new_independent_review_required': True}
write(O / 'final-guard.json', guard)

syntax = {}
for f in sorted(O.glob('*.py')):
    ast.parse(f.read_text(), filename=str(f))
    syntax[f.name] = sha(f)
assert load(O / 'format.json')['source_sha256'] == sha(W / path)
write(O / 'syntax.json', {'Python_AST': syntax, 'registry_toml': 'PASS', 'one_matching_Typst_title': True, 'pinned_rustfmt_record': load(O / 'format.json'), 'Clippy_full_docs_credit': False})
typst = Path('/nix/store/wzpxqpc3v2xxah69rqifq5kac67hx36g-typst-0.15.1/bin/typst')
args = [str(typst), 'compile', '--root', str(W), str(NOTE), str(O / 'canonical-note.pdf')]
p = subprocess.run(args, cwd=W, env=ENV, capture_output=True, text=True)
write(O / 'typst-check.json', {'argv': args, 'exit': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr, 'tool_sha256': sha(typst), 'source_sha256': sha(NOTE), 'output_sha256': sha(O / 'canonical-note.pdf') if p.returncode == 0 else None, 'scope': 'Only this standalone note. No full docs/site suite credit.'})
assert p.returncode == 0, p.stderr

# Evidence files are never placed in a child-writable root. The actual archive,
# extraction, fingerprint files and other twelve-row dumps stay private for review.
qualification = E / 'qualification'
qualification.mkdir(exist_ok=False)
for f in sorted(O.iterdir()):
    if not f.is_file() or f.name in ['finish.stdout', 'finish.stderr', 'canonical-note.pdf', 'local-artifacts.json']:
        continue
    shutil.copyfile(f, qualification / f.name)
for d in sorted((O / 'runs').iterdir()):
    dest = E / 'runtime' / d.name
    dest.mkdir(parents=True, exist_ok=False)
    for f in d.iterdir():
        if f.is_file():
            shutil.copyfile(f, dest / f.name)
local = {str(f.relative_to(O)): sha(f) for f in sorted(O.rglob('*')) if f.is_file() and str(f.relative_to(O)) not in ['finish.stdout', 'finish.stderr', 'local-artifacts.json']}
write(O / 'local-artifacts.json', {'root': str(O), 'files': local, 'exclusions': ['local-artifacts.json', 'finish.stdout', 'finish.stderr'], 'reason': 'Self manifest and live seal-driver streams. All completed build/runtime/capture evidence copied and sealed separately.', 'retained_for': 'NEW independent integration review: actual archive extraction, fingerprints, full twelve-row dumps, failed setup and outputs. No owned background process remains.'})
shutil.copyfile(O / 'local-artifacts.json', qualification / 'local-artifacts.json')
receipt = {'schema': 'raw-three-token-integration-candidate/1', 'status': 'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED', 'created_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'author': 'Lucien Huber <im@lcnbr.ch>', 'workspace': str(W), 'accepted_SM_base': BASE, 'new_code_commit': CODE, 'code_bookmark': 'generation-raw-three-migration-code-20261031', 'expected_candidate_bookmark': 'generation-raw-three-migration-candidate-20261031', 'literal_approval': {'review_sha256': sha(E / 'reviews/migration/review.json'), 'verdict': load(E / 'reviews/migration/review.json')['verdict'], 'root_explicitly_adopted': True, 'new_integration_self_approval': False, 'before_blob': 'e6af4ec644688e1333064c50e3b157ae631473ee', 'after_blob': '9ef9f1ba005f20149716b22346d45c5501b62e87', 'patch_sha256': sha(E / 'preparation/proposed-not-applied.patch'), 'exact_three_tokens_only': True, 'all149_literal_Rows_10_fields_prior9_every_other_byte_preserved': True}, 'actual_producer': {'uri': Q['uri'], 'source_paths': Q['source_paths'], 'selected_source_files': 670, 'selected_source_delta': Q['selected_source_delta'], 'archive_drv': Q['archive_drv'], 'archive': Q['archive'], 'archive_sha256': Q['archive_sha256'], 'oracle_elf': Q['oracle_elf'], 'profile': Q['profile'], 'compiler': Q['compiler'], 'kernel_packages': Q['kernels'], 'actual_root_fingerprint': Q['actual_root_fingerprint'], 'compiled_dependency_units': Q['closure_units'], 'compiled_layers': Q['compiled_layers'], 'producer_derivations': len(load(O / 'producer-chain.json')['derivations']), 'one_pure_package_archive_realization': True, 'actual_new_ELF_not_SM605': True, 'reuse_policy': 'Only exact unchanged dependency-layer bytes/fingerprints. No donor cache mount, relink, feature/pin/dependency/budget change.', 'pure_chronology': load(O / 'pure-chronology.json')}, 'runtime': r, 'runtime_log_audit': load(O / 'runtime-log-audit.json'), 'capture_joins': {'path': 'capture-joins.json', 'sha256': sha(O / 'capture-joins.json'), 'actual_new_files': 15, 'mode_cells': 6, 'owners': 400, 'affected': 80, 'unaffected': 320, 'all_physical_payloads_exact_byte_equal': True, 'normalizations': [], 'Auto_native_Row_assertions': 0, 'Auto_report_counter_limit': 'Existing public hook emits the returned None report only. Auto topology/assignment counters retain the explicit approved original-capture premise; actual Auto JSON/source payloads, names, identities, interfaces and all other serialized physics are fresh exact joins. ee4f7 vs42/name permutation preserved.'}, 'capsule': {'new_command_stdout_stderr_outside_child_writable_mount': True, 'SMINT_N1_historical_caveat_preserved': True, 'empty_root_bwrap_unshare_all': True, 'source_archive_input_runtime_files_readonly': True, 'host_roots_home_credentials_Nix_socket_compiler_donorcache_absent': True, 'unlicensed_loader_netns_input_preflight': 'PASS', 'license': 'Exported value preserved else current XDG config; private memory -> anonymous pipe -> immediate native exec only. No value in argv/logs/files/hash/Nix/pure child environments.', 'Nix_after_licensed_phase': 0}, 'finite_semantics': {'review': 'd9303d91d922c994f87b07e10cce895b44dba96e294d64be231a3ecdf921b14c', 'scope': 'Three separate finite cases, None/Auto independently, GroupingNone, math pool1. Pooled8/32 are runtime only.', 'operand_and_output_maps': 'H[i]=C[p[i]], b[p[i]]=a[i]. Recompute named plans or transport every positional operand AND retain output-axis map. No historical positional output identity or arbitrary stale ProductPlan interchangeability.', 'M2': '12 flag controls bind claims to actual imported owners, not serialized flags.', 'M3': 'Separate formatted control-generator repair is not assumed approved; neither generator used here.', 'native_reuse': '40 successful algebra batches/3088 obligations only; SM32 excluded. No old observer/algebra/helper/producer reruns.', 'source_roles': {'historical': '2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94', 'proof': '5ba04864ca506e259d9eef4a54f65acdd8323fc2', 'reconciled': 'f47f27213a1c55291239c9bc0f334a7ea1cd1cc1', 'prior_built': '5ea149ffd88ae3afc4072a38c0871c1329be8424', 'actual_new_built': CODE, 'kernel': '8eb7c3729d6abefc916decc1bdf0c06536f8af80'}, 'premises': 'Exact original raw/cooked/decoded capture association, registration state, complete symbol definitions, local/aggregate/projector meaning and usual compiler/build honesty remain finite premises. No universal decoder/locality/history/component-physics theorem.'}, 'validation': {'guard': guard, 'pinned_rustfmt': 'PASS', 'source_scope_118': 'PASS unchanged', 'standalone_Typst': load(O / 'typst-check.json'), 'Python_and_registry_syntax': 'PASS', 'Rust_policy_departures': None, 'new_Clippy_or_full_docs_credit': False, 'existing_build_warnings': 'Clarabel unused import and elided lifetimes remain in retained build log.'}, 'failures_retained': load(O / 'setup-failures.json'), 'remaining_gates': ['NEW independent integration review of this actual candidate', 'Auto mathematical review and root aggregation; only Auto ScopeA approved', 'host API/stubs/browser failed qualification and separate corrective work', 'other371/full oracle/heavy, performance/memory, wholeWP3/M-rebuild and combined integration', 'full docs/site, distribution and publication'], 'process': {'CI': False, 'generation_integration_advanced': False, 'main_index_frozen_refs_or_history_changed': False, 'delegations_push_fetch_PR_upload_publication': False, '20_GiB_policy_adopted': False, 'owned_background_processes_remaining': []}, 'retention': {'tracked': str(E), 'ignored_holder': str(O), 'purpose': 'Actual archive extraction, new fingerprints, complete queries/results, full selected dumps and failed setup retained for independent review. All prior packets/shared caches/other workers artifacts untouched.'}}
write(E / 'capture-joins.json', load(O / 'capture-joins.json'))
write(E / 'receipt.json', receipt)
files = {str(f.relative_to(E)): sha(f) for f in sorted(E.rglob('*')) if f.is_file() and str(f.relative_to(E)) not in ['SHA256SUMS', 'SEAL.json']}
(E / 'SHA256SUMS').write_text(''.join(h + '  ' + p + '\n' for p, h in files.items()))
write(E / 'SEAL.json', {'schema': 'raw-three-migration-evidence-seal/1', 'status': 'NEW_INDEPENDENT_INTEGRATION_REVIEW_REQUIRED', 'code_commit': CODE, 'accepted_base': BASE, 'manifest_sha256': sha(E / 'SHA256SUMS'), 'receipt_sha256': sha(E / 'receipt.json'), 'capture_joins_sha256': sha(E / 'capture-joins.json'), 'canonical_note_sha256': sha(NOTE), 'files': len(files), 'bytes': sum((E / p).stat().st_size for p in files), 'actual_selected_executions': 36, 'selected_retries_or_skips': 0, 'new_self_or_combined_approval': False, 'exclusions': ['SHA256SUMS (bound by this seal)', 'SEAL.json (external hash returned in committed-seal record)']})
print('SEALED', len(files), 'files; 36 actual passes; 15 exact capture joins; all118 scopes unchanged.', flush=True)
print('receipt', sha(E / 'receipt.json'), 'seal', sha(E / 'SEAL.json'), flush=True)
