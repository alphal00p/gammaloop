"""Pure byte qualification and unlicensed trusted runtime probes; no private reader."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import struct
import sys

W = Path(__file__).resolve().parents[2]
P = Path(__file__).resolve().parent
O = W / 'target/raw-next-qcd-proof-preparation-20261031-r3'
R1 = W / 'target/raw-next-qcd-proof-preparation-20261031-r1'
ASSIGN = W.parent / 'generation-feynkit-refresh/.local/feynkit-refresh/root-next-qcd-proof-native-execution-assignment-20261031-r3.json'
PYTHON = '/nix/store/4168xqs4l5jbfc7xyi9zrhdhszx0pnlb-python3-3.13.16/bin/python3'
PRIVATE = {'SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED', 'SYMBOLICA_LICENSE_FILE'}
sys.path.insert(0, str(O))
import supervisor as S
import supervisor_guard as K
from lifecycle import bounded

def load(p):
    return json.loads(Path(p).read_bytes(), object_pairs_hook=S.A.unique)

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def save(name, value):
    with (P / name).open('x') as f:
        json.dump(value, f, indent=2)
        f.write('\n')

def elf(path):
    """Read ELF64 little-endian headers only; never execute the object."""
    p = Path(path)
    with p.open('rb') as f:
        h = f.read(64)
        assert h[:6] == b'\x7fELF\x02\x01', str(p)
        phoff = struct.unpack_from('<Q', h, 32)[0]
        entsize, count = struct.unpack_from('<HH', h, 54)
        f.seek(phoff)
        headers = [struct.unpack('<IIQQQQQQ', f.read(entsize)) for _ in range(count)]
        interp = None
        dynamic = []
        for typ, flags, off, va, pa, filesz, memsz, align in headers:
            if typ == 3:
                f.seek(off); interp = f.read(filesz).rstrip(b'\0').decode()
            if typ == 2:
                f.seek(off)
                dynamic = [struct.unpack('<QQ', f.read(16)) for _ in range(filesz // 16)]
        tags = {}
        for tag, value in dynamic:
            if tag == 0: break
            tags.setdefault(tag, []).append(value)
        strings = b''
        if 5 in tags:
            address = tags[5][0]
            segment = next(h for h in headers if h[0] == 1 and h[3] <= address < h[3] + h[5])
            f.seek(segment[2] + address - segment[3])
            strings = f.read(tags[10][0])
        def string(n): return strings[n:strings.index(b'\0', n)].decode()
        return {'path': str(p), 'resolved': str(p.resolve()), 'sha256': sha(p), 'interpreter': interp,
                'needed': [string(v) for v in tags.get(1, [])],
                'rpath': [v for s in tags.get(15, []) for v in string(s).split(':')],
                'runpath': [v for s in tags.get(29, []) for v in string(s).split(':')]}

def main():
    assert not PRIVATE.intersection(os.environ)
    assert not (O / 'execution').exists()
    assert sha(ASSIGN) == '650f597bb727c64e53898d676f2884ee2125a5e73be8a8c0f0edf28376ba1382'
    assignment = load(ASSIGN)
    review = Path(assignment['supervisor_review']['path'])
    assert sha(review) == '6ad6242a78d52d0b371fb39d3942e66ed32bbbafec8878a91e6a72d5e3a9208b'
    assert sha(review.parent / 'SEAL.json') == assignment['supervisor_review_seal_sha256']
    manifests = [O, O.parent / 'raw-next-qcd-proof-preparation-20261031-r2', R1, review.parent]
    originals = load(O / 'originals-before.json')
    inventory = {}
    # Preserve the entire declared old evidence domain, including prior reviews,
    # captures and reference inputs; do not import any scientific implementation.
    for path, value in originals.items():
        expected = value['sha256'] if isinstance(value, dict) else value
        assert sha(path) == expected, path
        inventory[path] = {'sha256': expected, 'bytes': Path(path).stat().st_size}
    manifest_counts = {}
    for directory in manifests:
        seal = load(directory / 'SEAL.json')
        assert sha(directory / 'SHA256SUMS') == seal['manifest_sha256']
        lines = (directory / 'SHA256SUMS').read_text().splitlines()
        assert len(lines) == seal['bound_files']
        for line in lines:
            expected, rel = line.split('  ', 1)
            path = (directory / rel).resolve()
            assert path.is_relative_to(directory)
            assert sha(path) == expected, str(path)
            inventory[str(path)] = {'sha256': expected, 'bytes': path.stat().st_size}
        for name in ['SEAL.json', 'SHA256SUMS']:
            path = directory / name
            inventory[str(path)] = {'sha256': sha(path), 'bytes': path.stat().st_size}
        manifest_counts[str(directory)] = len(lines)
    cards = load(O / 'native-commands.json')['commands']
    assert assignment['command_ids'] == [r['id'] for r in cards]
    assert len(cards) == 21
    for field, name in [('supervisor_source_sha256','supervisor.py'), ('supervisor_guard_sha256','supervisor_guard.py'), ('native_commands_sha256','native-commands.json'), ('result_contracts_sha256','result-contracts.json'), ('ordered_input_bindings_sha256','ordered-input-bindings.json')]:
        assert sha(O / name) == assignment[field]
    bundles = []
    for index, (card, binding, prospective) in enumerate(zip(cards, assignment['program_input_order_result_bindings'], load(O / 'prospective-supervisor-commands.json'))):
        b = K._check_bundle(card['id'])
        a = K.verify_assignment(ASSIGN, b['supervisor_seal_sha256'])
        expected = {'id': card['id'], 'program': card['program'],
                    'input': {'path': card['input'], 'sha256': card['input_sha256']},
                    'card_sha256': S.A.sha(json.dumps(card, sort_keys=True, separators=(',', ':')).encode()),
                    'ordered_input_binding': b['definition_order'], 'result_contract': b['contract']}
        assert binding == expected
        argv = prospective['FUTURE_NOT_AUTHORIZED']
        assert argv == [PYTHON, '-I', '-B', str(O / 'supervisor.py'), '--id', card['id'], '--execute-licensed', '--assignment', str(ASSIGN)]
        for path, digest in b['protected_inputs'].items():
            inventory[path] = {'sha256': digest, 'bytes': Path(path).stat().st_size}
        bundles.append({'index': index, 'id': card['id'], 'bundle_and_assignment': 'PASS', 'ordered_input': load(b['definition_order']['path']), 'contract': b['contract']})
        print('PURE', index + 1, card['id'], 'PASS', flush=True)
    for path in [ASSIGN, review, Path(assignment['root_supervisor_acceptance'])]:
        inventory[str(path)] = {'sha256': sha(path), 'bytes': path.stat().st_size}
    save('inputs-before.json', inventory)
    save('pure-bundles.json', bundles)
    save('planned-invocations.json', [{'index': i, 'id': r['id'], 'argv': [PYTHON,'-I','-S','-B',str(O / 'supervisor.py'),'--id',r['id'],'--execute-licensed','--assignment',str(ASSIGN)], 'cwd': str(W)} for i,r in enumerate(cards)])
    save('pure-qualification.json', {'status':'PASS', 'utc':datetime.datetime.now(datetime.UTC).isoformat(), 'manifests':manifest_counts, 'input_files':len(inventory), 'all_four_private_names_absent':True, 'private_loader_imported':False, 'execution_root_absent':not (O/'execution').exists(), 'assignment':a, 'S_flag':'User-requested -S adds stdlib-only isolation to reviewed prospective argv; pinned Python, script, ID and assignment are identical.'})
    # Static transitive ELF resolution in the exact native mount namespace.
    selected = [cards[0], cards[1], cards[4]]
    configs = []
    for card in selected:
        argv = card['argv']; divider = argv.index('--')
        mounts = [(argv[i+1], argv[i+2]) for i,v in enumerate(argv[:divider]) if v == '--ro-bind']
        def exposed(path):
            resolved = Path(path).resolve()
            return any(resolved == Path(src).resolve() or (Path(src).is_dir() and resolved.is_relative_to(Path(src).resolve())) for src,dest in mounts if src == dest)
        seen = {}
        missing = []
        def walk(path, inherited=()):
            path = str(Path(path))
            if path in seen: return
            info = elf(path); seen[path] = info
            info['resolved_needed'] = {}
            origin = str(Path(path).parent)
            expand = lambda s:s.replace('${ORIGIN}',origin).replace('$ORIGIN',origin)
            dirs = [expand(d) for d in (info['rpath'] or info['runpath'])] + list(inherited)
            if info['interpreter']:
                interp = info['interpreter']
                if not exposed(interp): missing.append({'object':path,'interpreter':interp,'reason':'not mounted'})
                elif not Path(interp).is_file(): missing.append({'object':path,'interpreter':interp,'reason':'absent'})
                else: walk(interp)
            child_inherited = tuple([expand(d) for d in info['rpath']] + list(inherited))
            for name in info['needed']:
                candidates = [str(Path(d) / name) for d in dirs]
                found = next((p for p in candidates if Path(p).is_file() and exposed(p)), None)
                if not found: missing.append({'object':path,'needed':name,'search':candidates,'reason':'unqualified or missing from exact readonly mounts'})
                else:
                    info['resolved_needed'][name] = found
                    walk(found, child_inherited)
        walk(card['program']['path'])
        walk(PYTHON)
        config = {'id':card['id'], 'program':card['program'], 'readonly_mounts':mounts, 'objects':seen, 'missing':missing}
        configs.append(config)
    save('elf-qualification.json', {'configurations':configs, 'status':'FAIL' if any(c['missing'] for c in configs) else 'PASS'})
    assert not any(c['missing'] for c in configs), 'exact mounted ELF dependency boundary; no native/private phase allowed'
    # Run the actual default pure CLI without the four private names, saving
    # its unmodified complete output outside all selected native leaf paths.
    card = cards[0]
    with (P/'default-pure.stdout').open('xb') as out, (P/'default-pure.stderr').open('xb') as err:
        result = bounded([PYTHON,'-I','-S','-B',str(O/'supervisor.py'),'--id',card['id']], dict(os.environ), W, out, err, seconds=240)
    save('default-pure-lifecycle.json',result)
    assert result['exit'] == 0 and result['process_reaped'] and not result['timed_out']
    assert not (O/'execution').exists()
    print('PURE_AND_STATIC_ELF_PASS', flush=True)

if __name__ == '__main__':
    main()
