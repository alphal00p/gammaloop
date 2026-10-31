"""Independent streaming byte verification; stdlib only, no candidate execution."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PACKETS = [
    ('raw-next-qcd-preflight-20261031-r1', 936, '3e31a23aeff6de506bf6fb98ff3d0961a6d8bc49cb0576c57ad4ec9ca17cdb6d', 'a8cf51bb72bb650604532452cc560952d8067f31acb75a78366054e8106f5304'),
    ('raw-next-qcd-captures-20261031-r1', 1764, 'd73651512a36aad066ca52734bd7b560c0372f92245b1f6e038864cd7e454666', 'a5383b9bd408f6bff2a0afdce5e2e50b537bf84c84fe7e14093187c0bd27aebf'),
    ('raw-next-qcd-captures-20261031-r2', 1976, '9da40bc2d34381c7c6b67fa81a6d91032ff699222349f6910f9787a18cc7e660', '166dd40ee056d9bf67da1a84eef200dc5bc9ff0fea55e15a88e9feb8f4e98d4b'),
]
cache = {}
def sha(path):
    path = Path(path)
    key = str(path.absolute())
    if key not in cache:
        h = hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
                h.update(block)
        cache[key] = h.hexdigest()
    return cache[key]

report = {'packets': [], 'errors': [], 'native_executions': 0, 'candidate_imports': 0}
for name, count, handoff, seal_hash in PACKETS:
    base = ROOT / 'target' / name
    seal = json.loads((base / 'SEAL.json').read_text())
    assert sha(base / 'SEAL.json') == seal_hash
    assert sha(base / 'handoff.json') == handoff == seal['handoff_json_sha256']
    assert sha(base / 'SHA256SUMS') == seal['manifest_sha256']
    lines = (base / 'SHA256SUMS').read_text().splitlines()
    assert len(lines) == count == seal['files']
    paths = set()
    for line in lines:
        expected, text = line.split('  ', 1)
        path = Path(text)
        if not path.is_absolute():
            path = base / path
        # Only sealed repository evidence and explicit store objects, never user credentials.
        assert str(path).startswith(('/shared/lcnbr/dev/gammaloop/', '/nix/store/'))
        assert str(path) not in paths
        paths.add(str(path))
        try:
            actual = sha(path)
            if actual != expected:
                report['errors'].append({'path': str(path), 'expected': expected, 'actual': actual})
        except OSError as error:
            report['errors'].append({'path': str(path), 'error': str(error)})
    report['packets'].append({'name': name, 'bindings': count, 'manifest_sha256': seal['manifest_sha256'], 'handoff_sha256': handoff, 'seal_sha256': seal_hash})
report['unique_files_hashed'] = len(cache)
report['bindings_verified'] = sum(p['bindings'] for p in report['packets'])
if sys.argv[1] == 'after':
    original = json.loads((OUT / 'frozen-before-files.json').read_text())
    assert all(cache[p] == expected for p, expected in original.items())
    supplemental = json.loads((OUT / 'supplemental-bindings.json').read_text())
    for path, expected in supplemental.items():
        assert sha(path) == expected, path
    report['original_before_after_equal'] = True
    report['supplemental_bindings_unchanged'] = len(supplemental)
    report['total_unique_files_hashed_with_supplemental'] = len(cache)
report['pass'] = not report['errors']
(OUT / ('frozen-' + sys.argv[1] + '.json')).write_text(json.dumps(report, indent=2) + '\n')
(OUT / ('frozen-' + sys.argv[1] + '-files.json')).write_text(json.dumps(cache, indent=2) + '\n')
print(json.dumps(report))
