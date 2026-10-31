"""Only authorized fresh runtime. No Nix operations, compiler, or algebra helper."""
import os
import sys
EXPORTED = os.environ.pop('SYMBOLICA_LICENSE', None)
os.environ.pop('SYMBOLICA_LICENSE_PATH', None)
os.environ.pop('SYMBOLICA_LICENSE_SIGNED', None)
sys.path.insert(0, str(__file__).rsplit('/', 1)[0])
from runtime_common import *
import xml.etree.ElementTree as ET

assert load(O / 'selected-inventory.json')['status'] == 'READY_FOR_AUTHORIZED_RUNTIME'
value = EXPORTED
if not value:
    config = Path(os.environ.get('XDG_CONFIG_HOME', str(Path.home() / '.config')))
    value = (config / 'symbolica/license').read_text().strip()
assert value
results = []
for pool in [1, 8, 32]:
    cmd = [TOOLS['python'], '-I', '-S', '-B', '/input/exec.py', TOOLS['nextest'], 'nextest', 'run', *reuse(), '--profile', 'ci_gammaloop', '--run-ignored', 'all', '--no-tests', 'fail', '--no-fail-fast', '--retries', '0', '--success-output', 'final', '--status-level', 'all', '--final-status-level', 'all', '-E', FILTER]
    root, e, rc = execute('pool' + str(pool), cmd, pool=pool, value=value, dump=(pool == 1))
    f = root / 'nextest/ci_gammaloop/junit.xml'
    assert f.is_file()
    tests = ET.parse(f).findall('.//testcase')
    assert len(tests) == 12 and {t.attrib['name'] for t in tests} == set(CASES)
    statuses = [{'case': t.attrib['name'], 'status': 'skipped' if t.find('skipped') is not None else 'failed' if any(t.find(k) is not None for k in ['failure', 'error']) else 'passed', 'time': t.attrib.get('time')} for t in tests]
    (e / 'junit.xml').write_bytes(f.read_bytes())
    result = {'pool': pool, 'exit': rc, 'outcomes': statuses, 'actual_executions': sum(t['status'] != 'skipped' for t in statuses), 'junit_sha256': sha(f), 'command': str(e / 'command.json'), 'public_dump_hook_enabled': pool == 1}
    results.append(result)
    write(O / 'runtime-results.json', {'code': Q['code'], 'archive_sha256': Q['archive_sha256'], 'oracle_elf': Q['oracle_elf'], 'results': results, 'expected_selected_executions': 36, 'actual_executions': sum(r['actual_executions'] for r in results), 'passes': sum(t['status'] == 'passed' for r in results for t in r['outcomes']), 'failures': sum(t['status'] == 'failed' for r in results for t in r['outcomes']), 'retries': 0, 'skipped_credit': 0, 'finite_math_scope': 'Existing three separately approved pool1 None/Auto proofs only; pools8/32 runtime only.'})
    print('pool', pool, ':', [(t['case'], t['status']) for t in statuses], flush=True)
    assert rc == 0 and all(t['status'] == 'passed' for t in statuses), 'Actual divergence retained; expectations must not be weakened.'
value = EXPORTED = None
print('36 actual full-Row passes; zero retries/skipped credit. Existing hook captured pool1 modes.', flush=True)
