"""Pure loader, namespace and compiled-selection preflight before license access."""
import importlib.util
import sys
sys.path.insert(0, str(__file__).rsplit('/', 1)[0])
from runtime_common import *

assert not any(k in os.environ for k in ['SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED', 'SYMBOLICA_LICENSE_FILE'])
assert load(O / 'nix-complete.json')['all_nix_queries_complete']
SM = O.parent / 'sm-eeww-migration-20261031'
oldq = load(SM / 'archive-qualification.json')
assert Q['oracle_elf']['sha256'] != oldq['oracle_elf']['sha256']
assert Q['archive_sha256'] != oldq['archive_sha256']
assert Q['closure_units'] == len(Q['closure']) and Q['closure_units'] > 0
spec = importlib.util.spec_from_file_location('elfparser', O / 'elf_helpers.py')
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)
INPUT.mkdir(exist_ok=False)
config = (SOURCE / '.config/nextest.toml').read_bytes()
assert b'[store]' not in config
(INPUT / 'nextest.toml').write_bytes(config + b'\n[store]\ndir = "/work/nextest"\n')
(INPUT / 'exec.py').write_text('''import os,sys
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'))
for p in ('/source','/archive','/input'):assert os.statvfs(p).f_flag & os.ST_RDONLY
parts=[]
while b:=os.read(0,8192):parts.append(b)
fd=os.open('/dev/null',os.O_RDONLY);os.dup2(fd,0);os.close(fd)
value=b''.join(parts).decode();assert value.strip()
e=dict(os.environ);e['SYMBOLICA_LICENSE']=value
os.execve(sys.argv[1],sys.argv[1:],e)
''')
(INPUT / 'preflight.py').write_text('''import os,json
assert not any(k in os.environ for k in ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE'))
for p in ('/source','/archive','/input'):assert os.statvfs(p).f_flag & os.ST_RDONLY
for p in ('/shared','/root','/home/lcnbr','/run/nix/daemon-socket/socket'):assert not os.path.exists(p)
print(json.dumps({'netns':os.readlink('/proc/self/ns/net'),'readonly':['/source','/archive','/input'],'host_roots_credentials_nix_socket_absent':True,'license_absent':True}))
''')
old = load(SM / 'runtime-closure.json')
files = {}
records = {}
def add(target, source):
    source = Path(source).resolve()
    assert source.is_file()
    assert str(target).startswith('/nix/store/')
    assert target not in files or files[target] == str(source)
    files[str(target)] = str(source)

for path, record in old['ELF_RUNPATH_NEEDED'].items():
    p = Path(path)
    if not str(p).startswith('/nix/store/'):
        continue
    assert sha(p) == record['sha256']
    info = C.elf(p)
    assert all(info[k] == v for k, v in record.items() if k != 'sha256')
    records[str(p)] = {'sha256': sha(p), **info}
    add(str(p), p)
for r in old['resolutions']:
    p = Path(r['resolved'])
    assert sha(p) == r['sha256']
    if not str(p).startswith('/nix/store/'):
        continue
    add(str(p), p)
    requester = old['ELF_RUNPATH_NEEDED'][r['requester']]
    dirs = [Path(d) for d in requester.get('runpath', []) if d]
    if requester.get('interpreter'):
        dirs.append(Path(requester['interpreter']).parent)
    for d in dirs:
        alias = d / r['needed']
        if str(alias).startswith('/nix/store/') and alias.is_file() and alias.resolve() == p.resolve():
            add(str(alias), p)
add(TOOLS['python'], TOOLS['python'])
add(TOOLS['nextest'], TOOLS['nextest'])
new = C.elf(Q['oracle_elf']['path'])
prior = C.elf(oldq['oracle_elf']['path'])
assert new['needed'] == prior['needed'] and new['interpreter'] == prior['interpreter']
assert new['runpath'][1:] == prior['runpath'][1:]
assert new['runpath'][0] == Q['archive_store'] + '/lib' and not Path(new['runpath'][0]).exists()
assert not any('rust-stable' in p or '/bin/gcc' in p or '/bin/nix' in p or oldq['archive_store'] in p for p in files)
std = str(Path(TOOLS['python']).parents[1] / 'lib/python3.13')
assert Path(std).is_dir()
write(O / 'runtime-preparation.json', {'file_mounts': files, 'file_hashes': {p: sha(s) for p, s in files.items()}, 'stdlib_mounts': [std], 'stdlib_files': {str(f.relative_to(std)): sha(f) for f in sorted(Path(std).rglob('*')) if f.is_file()}, 'ELF_RUNPATH_NEEDED': records, 'new_oracle': {'sha256': sha(Q['oracle_elf']['path']), **new}, 'prior_SM_oracle_comparison': {'needed_equal': True, 'runpath_tail_equal': True, 'new_archive_only_first_path': True}, 'source_config_sha256': sha(SOURCE / '.config/nextest.toml'), 'input_hashes': {p.name: sha(p) for p in INPUT.iterdir()}, 'license_access': False, 'exact_file_library_mounts_no_compiler_packages': True})
host = os.readlink('/proc/self/ns/net')
root, e, rc = execute('preflight', [TOOLS['python'], '-I', '-S', '-B', '/input/preflight.py'])
assert rc == 0
r = load(e / 'stdout')
assert r['netns'] != host
write(O / 'network-preflight.json', {'host_netns': host, 'child': r, 'distinct_network_namespace': True})
root, e, rc = execute('loader-list', [new['interpreter'], '--list', '/archive/target/ci-optim/deps/' + Path(Q['oracle_elf']['path']).name])
assert rc == 0
root, e, rc = execute('inventory', [TOOLS['nextest'], 'nextest', 'list', *reuse(), '--message-format', 'json', '-E', FILTER])
assert rc == 0
inv = load(e / 'stdout')
found = [{'binary_id': bid, 'name': n, **t} for bid, s in inv['rust-suites'].items() for n, t in s['testcases'].items() if t['filter-match']['status'] == 'matches']
assert len(found) == 16 and {r['name'] for r in found} == set(CASES)
assert all(r['binary_id'] == 'feynkit-generator::generation_oracle' for r in found)
assert inv['rust-suites']['feynkit-generator::generation_oracle']['binary-path'] == '/archive/target/ci-optim/deps/' + Path(Q['oracle_elf']['path']).name
write(O / 'discovery-scope.json', {'authorization':'Bounded inventory discovery only','test_executions':0,'pass_credit':0,'selected_names':CASES,'command_record':str(e/'command.json')})
write(O / 'selected-inventory.json', {'status': 'READY_FOR_AUTHORIZED_RUNTIME', 'actual_compiled_selected': found, 'count': len(found), 'filter': FILTER, 'inventory_sha256': sha(e / 'stdout')})
print('PASS unlicensed exact-file loader/netns/read-only preflight and 16 exact compiled identities.', flush=True)
