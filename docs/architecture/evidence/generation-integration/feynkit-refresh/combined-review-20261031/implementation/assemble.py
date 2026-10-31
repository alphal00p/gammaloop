from pathlib import Path
import gzip
import hashlib
import json
import os
import subprocess
import tomllib

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
P = W / '.local/feynkit-refresh/combined-candidate-20261031-r1'
F = W / '.local/feynkit-refresh/combined-integration-preflight-20261031-r1'
E = 'docs/architecture/evidence/generation-integration/feynkit-refresh/combined-review-20261031'
BASE = '7f03ab8fd5834d7c78164429e2a44dab8d01ef74'
CODE = '6f41807b064bcf2f579883134daef65c8ac279d3'
ENV = {k: os.environ[k] for k in ['PATH', 'HOME', 'USER', 'LOGNAME', 'LANG', 'LC_ALL'] if k in os.environ}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def git(*args, input=None):
    return subprocess.check_output(['git', '--git-dir=' + str(W.parent / 'main/.git'), *args], env=ENV, input=input)

def emit(name, value):
    path = W / E / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')

def tree(rev):
    out = {}
    for line in git('ls-tree', '-rz', '--full-tree', rev).split(b'\0'):
        if line:
            header, path = line.split(b'\t', 1)
            mode, kind, oid = header.decode().split()
            assert kind == 'blob'
            out[os.fsdecode(path)] = {'mode': mode, 'blob': oid}
    return out

bindings = []
def copy(src, name, role):
    src = Path(src).resolve()
    raw = src.read_bytes()
    target = W / E / name
    record = {'source': str(src), 'logical_copy': E + '/' + name, 'sha256': sha(raw), 'bytes': len(raw), 'role': role}
    if len(raw) <= 1048576:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        record.update(encoding='identity', parts=[{'path': str(target.relative_to(W)), 'sha256': sha(raw), 'bytes': len(raw)}])
    else:
        encoded = gzip.compress(raw, mtime=0)
        parts = []
        for i, start in enumerate(range(0, len(encoded), 524288)):
            data = encoded[start:start + 524288]
            part = W / E / 'chunks' / (sha(raw) + '.' + str(i).zfill(4) + '.part')
            part.parent.mkdir(parents=True, exist_ok=True)
            part.write_bytes(data)
            parts.append({'path': str(part.relative_to(W)), 'sha256': sha(data), 'bytes': len(data)})
        record.update(encoding='gzip-chunks', parts=parts)
    bindings.append(record)
    return record

def manifest_entries(path):
    result = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        h, name = line.split(maxsplit=1)
        name = name.lstrip('*')
        assert len(h) == 64 and all(c in '0123456789abcdef' for c in h), (path, line)
        result.append((h, name))
    return result

def main():
    assert not (W / E).exists()
    before = tree(BASE)
    for f in sorted(F.iterdir()):
        if f.is_file():
            copy(f, 'preflight/' + f.name, 'verbatim prior preflight, NOT-APPLIED proposals remain historical')
    decisions = json.loads((F / 'decision-index.json').read_text())
    extra = [
        ('grouped-applied-r3', W / '.local/review-grouped-identical-B1-r3-astra-20261031-DSnqH402/review.json', '66db3d6dd9886e072fabdd78f3dab82d669b17408e5e05c10709f50e5eba9d7c'),
        ('grouped-applied-original', W / '.local/review-grouped-identical-integration-astra-20261031-033116Z-c0749b38/review.json', '0f0851e068653202a7981e71b4d5ccf54316d1bf544caa0725f99cc1794066f4'),
    ]
    names = ['host-integration', 'full-pool1', 'raw-three', 'host-prior', 'grouped-none', 'grouped-auto', 'finite-runner-r3', 'source-r4', 'shared-scope-a', 'metadata-48d0', 'source46', 'xs-original-a', 'xs-original-b', 'grouped-preparation']
    reviews = [(names[i], (W / d['review']).resolve(), d['expected_sha256']) for i, d in enumerate(decisions)] + extra
    # The prerequisite prior-eight review is already durably attributed by the original aggregation.
    old_index = json.loads((W / 'docs/architecture/evidence/generation-integration/feynkit-refresh/review-aggregation-20261030/reviews/manifest.json').read_text())
    prior = next(x for x in old_index['reviews'] if x['id'] == 'prior-eight')
    reviews.append(('prior-eight', Path(prior['source_directory']) / 'review.json', prior['review_sha256']))
    external, seals, review_index = [], [], []
    for name, review, expected in reviews:
        assert sha(review.read_bytes()) == expected, review
        parsed = json.loads(review.read_text())
        root_files = {review}
        for n in ['SHA256SUMS', 'MANIFEST.sha256', 'SEAL.json', 'seal.json', 'seal.sha256', 'review.json.sha256']:
            q = review.parent / n
            if q.is_file():
                root_files.add(q)
        for q in sorted(root_files):
            copy(q, 'reviews/' + name + '/' + q.name, 'verbatim independent review/seal; original path semantics retained')
        closure = {}
        for m in [q for q in root_files if q.name in ['SHA256SUMS', 'MANIFEST.sha256']]:
            entries = manifest_entries(m)
            seals.append({'review': name, 'source_manifest': str(m), 'sha256': sha(m.read_bytes()), 'entries': len(entries)})
            for h, n in entries:
                q = Path(n) if Path(n).is_absolute() else m.parent / n
                q = q.resolve()
                assert q.is_file(), (m, n)
                raw = q.read_bytes()
                assert sha(raw) == h, (m, n)
                closure[str(q)] = {'source': str(q), 'sha256': h, 'bytes': len(raw), 'manifest': str(m), 'manifest_name': n}
        for s in [q for q in root_files if q.name in ['seal.json', 'SEAL.json']]:
            spec = json.loads(s.read_text())
            fs = spec.get('files')
            if isinstance(fs, dict):
                for n, item in fs.items():
                    if not isinstance(item, dict) or 'sha256' not in item:
                        continue
                    q = (Path(n) if Path(n).is_absolute() else s.parent / n).resolve()
                    raw = q.read_bytes()
                    assert sha(raw) == item['sha256'], (s, n)
                    assert len(raw) == item.get('bytes', len(raw))
                    closure[str(q)] = {'source': str(q), 'sha256': sha(raw), 'bytes': len(raw), 'manifest': str(s), 'manifest_name': n}
        # All sealed inputs are explicitly retained external unless copied above. This is not a self-contained native/proof archive.
        external.extend(dict(r, review=name, retention='exact original retained review input; not claimed committed') for r in closure.values() if Path(r['source']) not in root_files)
        review_index.append({'id': name, 'original_review': str(review), 'review': E + '/reviews/' + name + '/review.json', 'sha256': expected, 'verdict': parsed.get('verdict', parsed.get('decision')), 'original_scope': parsed.get('scope', parsed.get('approved_scope')), 'sealed_external_entries': len(closure)})
    roots = sorted({d['acceptance'] for d in decisions if 'acceptance' in d} | {'.local/feynkit-refresh/root-grouped-applied-integration-acceptance-20261031.json'})
    for name in roots:
        copy(W / name, 'root-dispositions/' + Path(name).name, 'verbatim coordinator disposition, no new approval')
    emit('review-index.json', {'reviews': review_index, 'manifest_verifications': seals, 'current_combined_review': 'NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED'})
    # Historical global snapshots. Resolve each original manifest at its immutable evidence revision, not at today's mutable metadata.
    changed = ['docs/developers.toml', 'docs/architecture/generation-feynkit-community-refresh.typ', 'docs/architecture/evidence/generation-integration/execution.json']
    history = []
    for path in changed:
        raw = git('show', BASE + ':' + path)
        snapshot = 'historical/7f03ab8f/' + path.replace('/', '__') + '.verbatim'
        rec = copy(W / path, snapshot, 'exact pre-edit immutable evidence snapshot')
        assert rec['sha256'] == sha(raw)
        history.append({'original_path': path, 'immutable_evidence_revision': BASE, 'git_blob': before[path]['blob'], 'sha256': sha(raw), 'bytes': len(raw), 'durable_snapshot': rec['logical_copy'], 'meaning': 'pre-edit evidence metadata, not compiled/runtime source attribution'})
    emit('historical-snapshots.json', {'snapshots': history, 'rule': 'Historical SHA checks use the original immutable evidence revision or exact retained snapshot; current canonical metadata is a new revision. Never rewrite old seals or reviews.'})
    # All 425 and 29 entries are hash-verified now and again from the frozen candidate, with the exact registry snapshot substituted explicitly.
    historical_checks = []
    for folder in ['grouped-identical-migration', 'grouped-identical-migration-b1-r3']:
        m = W / 'docs/architecture/evidence/generation-integration/feynkit-refresh' / folder / 'SHA256SUMS'
        rows = []
        for expected, path in manifest_entries(m):
            raw = git('show', BASE + ':' + path)
            assert sha(raw) == expected, path
            rows.append({'path': path, 'sha256': expected, 'bytes': len(raw), 'git_blob': before[path]['blob'], 'immutable_evidence_revision': BASE})
        historical_checks.append({'manifest': str(m.relative_to(W)), 'manifest_sha256': sha(m.read_bytes()), 'evidence_revision': BASE, 'entries': rows, 'source_producer_revision': CODE, 'source_and_evidence_roles_distinct': True})
    emit('historical-manifest-bindings.json', {'manifests': historical_checks})
    emit('retained-review-inputs.json', {'self_contained': False, 'entries': external, 'all_rehashed': True, 'scope': 'All entries of the copied review SHA256 manifests and file-map seals. Other external scientific/runtime inputs stay under their original receipts and are not newly claimed committed.'})
    emit('copy-bindings.json', {'tracking_limit': 1048576, 'encoding': 'identity or lossless gzip split into at most 512KiB parts; logical names are not tracking claims', 'bindings': bindings})
    print(json.dumps({'copied_logical_files': len(bindings), 'review_count': len(review_index), 'external_sealed_entries': len(external), 'historical_manifest_entries': sum(len(m['entries']) for m in historical_checks)}))

if __name__ == '__main__':
    main()
