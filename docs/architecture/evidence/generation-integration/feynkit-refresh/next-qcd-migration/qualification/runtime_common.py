"""Task-local capsule orchestration; contains no Nix or license lookup."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import resource

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
O = Path(__file__).resolve().parent
Q = json.loads((O / 'archive-qualification.json').read_text())
ARCH = O / 'archive-extracted'
SOURCE = Path(Q['source'])
INPUT = O / 'capsule-input'
TOOLS = Q['runtime_tools']
PRIOR13 = json.loads((O / 'source-invariance.json').read_text())['prior13']
NEW3 = ['xs_qcd_gguu_sym_L1', 'xs_qcd_gguu_L1', 'xs_qcd_gggg_L1']
CASES = PRIOR13 + NEW3
FILTER = 'package(=feynkit-generator) and binary(=generation_oracle) and (' + ' or '.join('test(=' + n + ')' for n in CASES) + ')'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def load(p):
    return json.loads(Path(p).read_text())

def write(p, data):
    Path(p).write_text(json.dumps(data, indent=2) + '\n')

def reuse():
    return ['--cargo-metadata', '/archive/target/nextest/cargo-metadata.json', '--binaries-metadata', '/archive/target/nextest/binaries-metadata.json', '--workspace-remap', '/source', '--target-dir-remap', '/archive/target', '--config-file', '/input/nextest.toml', '--user-config-file', 'none']

def execute(name, command, pool=1, value=None, dump=False):
    assert not any(k in os.environ for k in ['SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED', 'SYMBOLICA_LICENSE_FILE'])
    assert load(O / 'nix-complete.json')['all_nix_queries_complete']
    root = O / 'work' / name
    evidence = O / 'runs' / name
    root.mkdir(parents=True, exist_ok=False)
    evidence.mkdir(parents=True, exist_ok=False)
    for n in ['tmp', 'home', 'cache', 'state', 'dumps']:
        (root / n).mkdir()
    workers = min(16, len(os.sched_getaffinity(0)))
    cores = int(os.environ.get('NIX_BUILD_CORES', '0'))
    if cores > 0:
        workers = min(workers, cores)
    e = {'PATH': str(Path(TOOLS['nextest']).parent), 'HOME': '/home/runtime', 'USER': 'runtime', 'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8', 'TMPDIR': '/work/tmp', 'XDG_CACHE_HOME': '/work/cache', 'GL_TEST_LOG_DIR': '/work/state', 'NEXTEST_TEST_THREADS': str(workers), 'NEXTEST_NO_INPUT_HANDLER': '1', 'CARGO_NET_OFFLINE': 'true', 'FEYNKIT_ORACLE_THREADS': str(pool)}
    if dump:
        e['FEYNKIT_ORACLE_DUMP'] = '/work/dumps'
    argv = [TOOLS['bwrap'], '--unshare-all', '--die-with-parent', '--new-session', '--proc', '/proc', '--dev', '/dev', '--bind', str(root), '/work', '--bind', str(root / 'tmp'), '/tmp', '--dir', '/home', '--bind', str(root / 'home'), '/home/runtime', '--ro-bind', str(INPUT), '/input', '--ro-bind', str(ARCH), '/archive', '--ro-bind', Q['archive'], '/archive.tar.zst', '--ro-bind', str(SOURCE), '/source']
    for target, source in sorted(load(O / 'runtime-preparation.json')['file_mounts'].items()):
        argv += ['--ro-bind', source, target]
    for p in load(O / 'runtime-preparation.json')['stdlib_mounts']:
        argv += ['--ro-bind', p, p]
    argv += ['--chdir', '/source', '--', *command]
    before = {str(p): sha(p) for p in [Path(Q['archive']), Path(Q['oracle_elf']['path']), *INPUT.iterdir()]}
    record = {'argv': argv, 'environment_before_exec': e, 'start': datetime.datetime.now(datetime.UTC).isoformat(), 'licensed': value is not None, 'license_delivery': 'anonymous stdin pipe -> isolated Python -> immediate native exec; no value recorded' if value is not None else 'UNLICENSED', 'command_records_outside_writable_capsule': True, 'writable_mount': str(root), 'evidence_root': str(evidence), 'inputs_before': before, 'empty_root_no_network': True}
    write(evidence / 'command.json', record)
    start = time.monotonic()
    ru_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (evidence / 'stdout').open('wb') as out, (evidence / 'stderr').open('wb') as err:
        p = subprocess.Popen(argv, cwd=W, env=e, stdin=subprocess.PIPE if value is not None else subprocess.DEVNULL, stdout=out, stderr=err)
        try:
            if value is not None:
                p.stdin.write(value.encode())
                p.stdin.close()
            rc = p.wait(timeout=900)
        except BaseException:
            p.kill()
            p.wait()
            raise
    ru_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    record['resource_observation'] = {'child_user_seconds_delta':ru_after.ru_utime-ru_before.ru_utime,'child_system_seconds_delta':ru_after.ru_stime-ru_before.ru_stime,'controller_cumulative_child_maxrss_KiB':ru_after.ru_maxrss,'maxrss_scope':'Cumulative reaped children of this controller; not a per-test isolated peak; no memory ceiling adopted'}
    after = {p: sha(p) for p in before}
    record.update(exit=rc, end=datetime.datetime.now(datetime.UTC).isoformat(), seconds=time.monotonic() - start, inputs_after=after, stdout_sha256=sha(evidence / 'stdout'), stderr_sha256=sha(evidence / 'stderr'), process_reaped=True)
    write(evidence / 'command.json', record)
    assert before == after
    return root, evidence, rc
