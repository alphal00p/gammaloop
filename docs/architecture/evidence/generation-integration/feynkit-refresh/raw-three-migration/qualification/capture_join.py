"""Join actual production dumps to approved physical inputs, without native proof runs."""
import copy
from collections import Counter
from fractions import Fraction
import re
import shutil
import sys
sys.path.insert(0, str(__file__).rsplit('/', 1)[0])
from runtime_common import *

assert not any(k in os.environ for k in ['SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED'])
F = W.parent / 'generation-oracle-numerator-families'
PROOF = F / 'target/raw-three-interface-proof-20261030-r1'
E = W / 'docs/architecture/evidence/generation-integration/feynkit-refresh/raw-three-migration'
review = load(E / 'reviews/migration/review.json')
proof = load(PROOF / 'inputs.json')
assert sha(PROOF / 'validation-receipt.json') == '7c443ba689cd2e04cddf4d6a2d1f60c4e41512cdad9cf64237e770b86d74349f'
refs = load(E / 'preparation/public-census-source-joins.json')
DUMPS = O / 'work/pool1/dumps'

def enc(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()

def fnv(data):
    v = 0xcbf29ce484222325
    for b in data:
        v = ((v ^ b) * 0x100000001b3) & 0xffffffffffffffff
    return v

def le(n):
    return n.to_bytes(8, 'little')

def text(s):
    return s.encode() + b'\xff'

def name(s):
    return b'\xfe' if s is None else text(s)

def bag(xs):
    return fnv(b''.join(le(x) for x in sorted(xs)))

def weight(s):
    # Exact restricted bookkeeping grammar present in these three finite cells.
    # Reject every other factor; never substitute a physical parameter.
    v = Fraction(1)
    for part in s.split('*'):
        m = re.fullmatch(r'\(feynkit_generator_factor::\{\}::AutG\((\d+)\)\)\^\(-1\)', part)
        if m:
            v /= int(m[1])
            continue
        m = re.fullmatch(r'feynkit_generator_factor::\{\}::(AntiFermionSpinSumSign|ExternalFermionOrderingSign|InternalFermionLoopSign)\((-?1)\)', part)
        assert m, ('unsupported physical factor', part)
        v *= int(m[2])
    return v

def digests(p):
    p = copy.deepcopy(p)
    parts = [enc(p.pop(k)) for k in ['numerator', 'numerator_prefactor', 'projector']]
    for k in ['num_interface', 'numerator_prefactor_interface', 'projector_interface']:
        if k in p:
            parts.append(enc(p.pop(k)))
    for o in [*p['vertices'], *[e[1] for e in p['edges']]]:
        parts.append(enc(o['numerator']))
        o['numerator'] = None
        if 'num_interface' in o:
            parts.append(enc(o.pop('num_interface')))
    for k in ['model', 'symmetry_factor', 'overall_factor']:
        p.pop(k)
    return fnv(enc(p)), fnv(b''.join(parts))

def derive(payload, report):
    field = lambda k: int(re.search(r'\b' + k + r': (\d+)', report)[1])
    assert 'completed: true' in report
    assert field('retained_count') == len(payload['diagrams'])
    syms = Counter()
    inv = Fraction()
    ws = Fraction()
    content, nums, ids = [], [], []
    named = b''
    for g in payload['diagrams']:
        k = fnv(enc(g['key']))
        frame, numerator = digests(g['payload'])
        s = int(g['aut'])
        assert s > 0
        syms[s] += 1
        inv += Fraction(1, s)
        ws += weight(g['overall_factor'])
        content.append(fnv(le(k) + le(s) + text(g['overall_factor']) + le(frame)))
        nums.append(fnv(le(k) + le(numerator)))
        ids.append(fnv(int(g['id'], 16).to_bytes(16, 'little')))
        named += name(g['name']) + le(k)
    for group in payload['groups']:
        named += le(group['master'])
        for member in group['members']:
            named += name(member['name']) + text(member['ratio']) + text(member['overall_factor'])
    return {'diagrams': len(payload['diagrams']), 'topologies': field('topology_count'), 'assignments': field('interaction_assignment_count'), 'inverse_symmetry': [inv.numerator, inv.denominator], 'weights': [ws.numerator, ws.denominator], 'symmetry_factors': [list(x) for x in sorted(syms.items())], 'content': bag(content), 'numerators': bag(nums), 'named': fnv(named), 'ids': bag(ids)}

def row(case):
    s = re.search(r'\b' + case + r': Row \{([^\n]+)\}', (W / 'crates/feynkit-generator/tests/generation_oracle/main.rs').read_text())[1]
    d = {k: int(re.search(r'\b' + k + r': (0x[0-9a-f]+|\d+)', s)[1], 0) for k in ['diagrams', 'topologies', 'assignments', 'content', 'numerators', 'named', 'ids']}
    for k in ['inverse_symmetry', 'weights']:
        d[k] = list(map(int, re.search(k + r': \((-?\d+), (\d+)\)', s).groups()))
    d['symmetry_factors'] = [list(map(int, t)) for t in re.findall(r'\((\d+), (\d+)\)', re.search(r'symmetry_factors: &\[([^\]]+)\]', s)[1])]
    return d

records = []
joins = []
for approved in review['per_row']:
    case = approved['case']
    binding = next(r for r in refs if r['case'] == case)
    for j in binding['public_joins']:
        fresh = DUMPS / j['file']
        assert sha(j['pristine']) == sha(j['observer']) == j['sha256']
        assert fresh.read_bytes() == Path(j['pristine']).read_bytes(), ('first unequal production boundary', case, j['file'])
        dest = E / 'captures' / case / j['file']
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(fresh, dest)
        joins.append({'case': case, 'file': j['file'], 'fresh': str(fresh), 'copy': str(dest.relative_to(E)), 'fresh_sha256': sha(fresh), 'approved_original': j['pristine'], 'approved_observer': j['observer'], 'approved_sha256': j['sha256'], 'byte_equal': True, 'normalizations': []})
    census = (DUMPS / (case + '.census')).read_text()
    # The existing public hook emits the returned None report only. Its Auto
    # topology/assignment work counters are retained proof premises, not fresh
    # Auto full-Row assertions. Every physical Auto payload byte is newly joined.
    for mode in ['None', 'Auto']:
        p = DUMPS / f'{case}.classes-{mode}-t1.json'
        d = load(p)
        assert p.read_bytes() == p.with_name(p.stem + '.sources.json').read_bytes()
        cell = next(c for c in proof['cells'] if c['case'] == case and c['mode'] == mode)
        assert cell['actual_request']['classes'] == mode and cell['grouping'] == 'none' and 'threads: Some(1)' in cell['Run']
        assert len(d['diagrams']) == len(cell['graphs'])
        owners = 0
        affected = 0
        for fresh, g in zip(d['diagrams'], cell['graphs']):
            assert fresh['id'] == g['id'] and fresh['name'] == g['name']
            assert fresh['payload'] == g['current_payload'], ('proof payload drift', case, mode, g['id'])
            assert fresh['payload']['model'] == g['current_payload']['model']
            for owner in g['owners']:
                n = owner['owner']
                if n.startswith('vertex:'):
                    actual = fresh['payload']['vertices'][int(n.split(':')[1])]
                elif n.startswith('edge:'):
                    actual = fresh['payload']['edges'][int(n.split(':')[1])][1]
                else:
                    key = {'numerator': 'numerator', 'numerator_prefactor': 'numerator_prefactor', 'projector': 'projector'}[n]
                    actual = {'numerator': fresh['payload'][key]}
                    interface = {'numerator': 'num_interface', 'numerator_prefactor': 'numerator_prefactor_interface', 'projector': 'projector_interface'}[n]
                    if interface in fresh['payload']:
                        actual['num_interface'] = fresh['payload'][interface]
                assert actual['numerator'] == owner['expression']
                assert ('num_interface' in actual) == owner['explicit']
                if owner['explicit']:
                    assert actual['num_interface'] == owner['current']
                    affected += 1
                owners += 1
        expected = approved['after_Row'] if mode == 'None' else approved['Auto_current_Row']
        report = census.splitlines()[0] if mode == 'None' else f"GenerationReport {{ topology_count: {expected['topologies']}, interaction_assignment_count: {expected['assignments']}, retained_count: {expected['diagrams']}, completed: true }}"
        actual = derive(d, report)
        assert actual == expected
        if mode == 'None':
            assert actual == row(case)
        old = approved['before_Row'] if mode == 'None' else approved['Auto_historical_Row']
        assert [k for k in actual if actual[k] != old[k]] == ['numerators']
        records.append({'case': case, 'mode': mode, 'grouping': 'None', 'pool': 1, 'physical_graphs': len(d['diagrams']), 'owners': owners, 'explicit_affected_owners': affected, 'fresh_payload_sha256': sha(p), 'proof_input_sha256': sha(PROOF / 'inputs.json'), 'all_payloads_full_keys_models_names_IDs_groups_sources_ordered_interfaces_byte_bound': True, 'Row': actual, 'other_nine_fields_value_equal_to_approved_historical_mode': True, 'actual_native_full_Row_assertion': mode == 'None', 'report_counter_origin': 'fresh None census report' if mode == 'None' else 'approved Auto report premise from original capture; public hook emits no Auto report', 'normalizations': []})
assert len(records) == 6 and sum(r['owners'] for r in records) == 400 and sum(r['explicit_affected_owners'] for r in records) == 80
assert len(joins) == 15
shutil.copyfile(PROOF / 'inputs.json', E / 'preparation/retained-proof/inputs.json')
write(O / 'capture-joins.json', {'schema': 'raw-three-actual-production-capture-join/1', 'code': Q['code'], 'archive_sha256': Q['archive_sha256'], 'oracle_elf': Q['oracle_elf'], 'production_hook': {'main_rs': 'check/FEYNKIT_ORACLE_DUMP -> census text', 'classes_rs': 'check_request/FEYNKIT_ORACLE_DUMP -> per-mode JSON and source companion', 'census_rs': 'Census::write_exact_dump -> exact ordered payload schema', 'changed': False, 'new_observer_or_algebra_runs': 0}, 'files': joins, 'cells': records, 'same_physical_payloads': True, 'owners': 400, 'affected': 80, 'unaffected': 320, 'original_math_only_pool1': True, 'byte_normalizations': [], 'ee4f_mode_limits': 'None literal Row topologies7/name order; independently fresh Auto payload with original report42 and actual mode name permutation. No duplicate Auto Row assertion.', 'new_native_proof_runs': 0, 'old_algebra_reuse': '40 successful batches/3088 obligations only; SM32 excluded', 'arbitrary_cached_positional_plans': 'UNSAFE; recompute named plans or transport each operand AND preserve output-axis map.'})
print('PASS 15 exact production capture joins, six separate mode cells, 400 owners/80 affected; no normalization or proof rerun.', flush=True)
