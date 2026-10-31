"""Pure audit of exact Git blobs, original seals, lossless copies and guards."""
from pathlib import Path
from collections import Counter
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import tomllib

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
P = W / '.local/feynkit-refresh/combined-candidate-20261031-r1'
E = 'docs/architecture/evidence/generation-integration/feynkit-refresh/combined-review-20261031'
BASE = '7f03ab8fd5834d7c78164429e2a44dab8d01ef74'
ENV = {k: os.environ[k] for k in ['PATH', 'HOME', 'USER', 'LOGNAME', 'LANG', 'LC_ALL'] if k in os.environ}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')

def sha(b):
    return hashlib.sha256(b).hexdigest()

def git(*args, input=None):
    return subprocess.check_output(['git', '--git-dir=' + str(W.parent / 'main/.git'), *args], env=ENV, input=input)

def tree(rev):
    result = {}
    for row in git('ls-tree', '-rz', '--full-tree', rev).split(b'\0'):
        if row:
            header, path = row.split(b'\t', 1)
            mode, kind, oid = header.decode().split()
            assert kind == 'blob'
            result[os.fsdecode(path)] = {'mode': mode, 'blob': oid}
    return result

def batch(oids):
    order = list(dict.fromkeys(oids))
    raw = git('cat-file', '--batch', input=''.join(x + '\n' for x in order).encode())
    result, pos = {}, 0
    for wanted in order:
        end = raw.index(b'\n', pos)
        oid, kind, size = raw[pos:end].decode().split()
        assert oid == wanted and kind == 'blob'
        size = int(size); pos = end + 1
        data = raw[pos:pos + size]; pos += size + 1
        assert hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest() == oid
        result[oid] = data
    assert pos == len(raw)
    return result

def main():
    rev = sys.argv[1]
    t, old = tree(rev), tree(BASE)
    added = sorted(set(t) - set(old))
    removed = sorted(set(old) - set(t))
    changed = sorted(p for p in set(t) & set(old) if t[p] != old[p])
    expected = ['docs/architecture/evidence/generation-integration/execution.json', 'docs/architecture/generation-feynkit-community-refresh.typ', 'docs/developers.toml']
    assert not removed and changed == expected, (removed, changed)
    assert all(p.startswith(E + '/') or p == 'docs/architecture/generation-combined-refresh-review.typ' for p in added)
    cache = batch([entry['blob'] for entry in t.values()])
    def read(p):
        assert p in t, ('untracked claim', p)
        oid = t[p]['blob']
        if oid not in cache:
            cache.update(batch([oid]))
        return cache[oid]
    for p in added:
        assert len(read(p)) <= 1048576, ('tracking limit', p, len(read(p)))
    for p in t:
        assert t[p]['mode'] in ['100644', '100755', '120000']
        if t[p]['mode'] != '120000':
            assert (W / p).read_bytes() == read(p), ('working copy differs from frozen Git', p)
    manifest_checks = []
    for mp in [E + '/SHA256SUMS', E + '/freeze/SHA256SUMS']:
        if mp not in t:
            continue
        entries = []
        for line in read(mp).decode().splitlines():
            h, p = line.split('  ', 1)
            data = read(p)
            assert sha(data) == h, (mp, p)
            entries.append({'path': p, 'sha256': h, 'bytes': len(data), **t[p]})
        spec = json.loads(read(mp.rsplit('/', 1)[0] + '/SEAL.json'))
        assert sha(read(mp)) == spec['manifest_sha256']
        assert len(entries) == spec['files']
        manifest_checks.append({'manifest': mp, 'sha256': sha(read(mp)), 'entries': entries})
    binding_index = json.loads(read(E + '/copy-bindings.json'))
    copied_records = []
    for r in binding_index['bindings']:
        chunks = []
        for part in r['parts']:
            data = read(part['path'])
            assert len(data) == part['bytes'] and sha(data) == part['sha256']
            chunks.append(data)
        raw = b''.join(chunks)
        if r['encoding'] == 'gzip-chunks':
            raw = gzip.decompress(raw)
        assert len(raw) == r['bytes'] and sha(raw) == r['sha256']
        copied_records.append({'logical_copy': r['logical_copy'], 'sha256': sha(raw), 'bytes': len(raw), 'tracked_parts': len(chunks)})
    histories = json.loads(read(E + '/historical-snapshots.json'))['snapshots']
    snapshots = {}
    for r in histories:
        raw = read(r['durable_snapshot'])
        historic = git('show', r['immutable_evidence_revision'] + ':' + r['original_path'])
        assert raw == historic and sha(raw) == r['sha256']
        assert old[r['original_path']]['blob'] == r['git_blob']
        snapshots[r['original_path']] = raw
    old_manifests = []
    for m in json.loads(read(E + '/historical-manifest-bindings.json'))['manifests']:
        assert sha(read(m['manifest'])) == m['manifest_sha256']
        results = []
        for r in m['entries']:
            p = r['path']
            # Revision-aware original global metadata; never compare today's new registry to the old seal.
            raw = snapshots[p] if p in snapshots else read(p)
            assert sha(raw) == r['sha256'] and len(raw) == r['bytes'], (m['manifest'], p)
            assert old[p]['blob'] == r['git_blob']
            results.append({'path': p, 'sha256': sha(raw), 'bytes': len(raw), 'verified_at': BASE if p in snapshots else rev, 'via_snapshot': p in snapshots})
        old_manifests.append({'manifest': m['manifest'], 'entries': results})
    # Reconstruct the existing sealed grouped copied closure from actual tracked blobs too.
    ge = 'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration'
    logical, reconstructed = {}, []
    for filename, key in [('input-bindings.json', 'bindings'), ('large-evidence-bindings.json', 'all_own_evidence_bindings')]:
        for r in json.loads(read(ge + '/' + filename))[key]:
            path = r['durable_copy']
            data = read(path)
            assert sha(data) == r['copy_sha256'] and len(data) == r['copy_bytes']
            if r['encoding'] == 'gzip-chunks':
                spec = json.loads(data)
                assert spec['parts'] == r['parts']
                pieces = []
                for part in r['parts']:
                    b = read(part['path'])
                    assert sha(b) == part['sha256'] and len(b) == part['bytes']
                    pieces.append(b)
                data = gzip.decompress(b''.join(pieces))
                target = path[:-len('.gz.index.json')]
            elif r['encoding'].startswith('gzip'):
                data = gzip.decompress(data)
                target = path[:-3]
            else:
                assert r['encoding'] == 'identity'
                target = path
            assert sha(data) == r['sha256'] and len(data) == r['bytes']
            assert target not in logical or logical[target] == data
            logical[target] = data
            reconstructed.append({'path': target, 'sha256': sha(data), 'bytes': len(data)})
    def copied(path):
        return logical[path] if path in logical else read(path)
    migration = []
    for n, r in json.loads(copied(ge + '/reviews/migration/seal.json'))['files'].items():
        data = copied(ge + '/reviews/migration/' + n)
        assert sha(data) == r['sha256'] and len(data) == r['bytes']
        migration.append({'path': n, 'sha256': sha(data), 'bytes': len(data)})
    preparation = []
    for line in copied(ge + '/preparation/SHA256SUMS').decode().splitlines():
        h, n = line.split('  ', 1)
        data = copied(ge + '/preparation/' + n)
        assert sha(data) == h
        preparation.append({'path': n, 'sha256': sha(data), 'bytes': len(data)})
    external = json.loads(read(E + '/retained-review-inputs.json'))['entries']
    for r in external:
        data = Path(r['source']).read_bytes()
        assert sha(data) == r['sha256'] and len(data) == r['bytes'], ('external sealed input changed', r['source'])
    reg = tomllib.loads(read('docs/developers.toml').decode())
    old_reg = tomllib.loads(git('show', BASE + ':docs/developers.toml').decode())
    def scopes(registry):
        return [(n['id'], s) for section in registry['section'] for n in section.get('note', []) for s in n.get('scope', [])]
    assert scopes(reg) == scopes(old_reg) and len(scopes(reg)) == 118
    scope_results = []
    for ident, s in scopes(reg):
        source = read(s['path']).decode()
        assert source.count(s['anchor']) == 1
        if 'end_anchor' in s:
            assert source.count(s['end_anchor']) == 1
            source = source[source.index(s['anchor']):source.index(s['end_anchor'])]
        data = source.encode()
        digest = hashlib.sha1(b'blob %d\0' % len(data) + data).hexdigest()
        assert digest == s['digest']
        scope_results.append({'note': ident, 'path': s['path'], 'digest': digest})
    for note_id in ['generation-feynkit-community-refresh', 'generation-combined-refresh-review']:
        note = next(n for section in reg['section'] for n in section.get('note', []) if n['id'] == note_id)
        source = read(note['source']).decode()
        titles = re.findall(r'^= (.+)$', source, re.M)
        assert titles == [note['title']]
        assert not any(k in note for k in ['reviewed_at', 'review_ref'])
        for target in re.findall(r'#link\("([^"#]+)(?:#[^"]*)?"\)', source):
            if '://' not in target:
                assert str(Path(note['source']).parent / target) in t, ('broken local link', note_id, target)
    for p in added + changed:
        if p.endswith('.json'):
            json.loads(read(p))
    guards = json.loads((P / 'guards-before.json').read_text())
    main = W.parent / 'main'
    assert sha((main / '.git/index').read_bytes()) == guards['index_sha256']
    for p, h in guards['main_files'].items():
        assert sha((main / p).read_bytes()) == h, ('main changed', p)
    bookmarks = subprocess.check_output(['jj', '--ignore-working-copy', 'bookmark', 'list', '-T', 'name ++ " " ++ normal_target.commit_id() ++ "\n"'], cwd=W).decode()
    # The original template intentionally retains local AND remote (@git) entries.
    # Preserve the complete multiset rather than losing one in a name-keyed dict.
    before_refs, now_refs = Counter(guards['bookmarks'].splitlines()), Counter(bookmarks.splitlines())
    assert not before_refs - now_refs, ('old bookmark entries changed', before_refs - now_refs)
    local_global = subprocess.check_output(['jj', '--ignore-working-copy', 'log', '-r', 'generation-integration', '--no-graph', '-T', 'commit_id'], cwd=W).decode().strip()
    assert local_global == 'fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
    assert git('rev-parse', rev + ':crates/feynkit-generator/tests/generation_oracle/main.rs').decode().strip() == '7283198bf0152c82c1cd5f0173703cddfa7f5c0f'
    assert read('nix-ci.nix') == git('show', BASE + ':nix-ci.nix')
    author = git('show', '-s', '--format=%an <%ae>', rev).decode().strip()
    assert author == 'Lucien Huber <im@lcnbr.ch>', author
    # Enumerate the whole immutable tree and SHA256 every blob (including ignored logs and global paths).
    vector = [{'path': p, **t[p], 'bytes': len(read(p)), 'sha256': sha(read(p))} for p in sorted(t)]
    result = {'schema': 'combined-refresh-immutable-audit/1', 'candidate': rev, 'tree': git('rev-parse', rev + '^{tree}').decode().strip(), 'base': BASE, 'author': author, 'full_tree': vector, 'full_tree_vector_sha256': sha(json.dumps(vector, sort_keys=True, separators=(',', ':')).encode()), 'added_paths': added, 'changed_paths': changed, 'removed_paths': removed, 'new_manifest_checks': manifest_checks, 'copied_logical_closure': copied_records, 'historical_manifest_checks': old_manifests, 'historical_snapshots': histories, 'original_grouped_lossless_bindings': reconstructed, 'original_grouped_migration_review_closure': migration, 'original_grouped_preparation_closure': preparation, 'external_sealed_inputs_rehashed': len(external), 'scope_digests': scope_results, 'old_bookmarks_main_index_and_production_preserved': True, 'build_native_project_runtime_license_config_network_or_delegation': 0, 'owned_background_processes': [], 'review': 'NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED'}
    path = P / ('immutable-' + rev + '-audit.json')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'candidate': rev, 'tree': result['tree'], 'full_tree_entries': len(vector), 'full_tree_vector_sha256': result['full_tree_vector_sha256'], 'new_manifests': [(x['manifest'], len(x['entries'])) for x in manifest_checks], 'historical_manifests': [(x['manifest'], len(x['entries'])) for x in old_manifests], 'copied_logical_files': len(copied_records), 'external_sealed_inputs_rehashed': len(external), 'audit': str(path), 'audit_sha256': sha(path.read_bytes()), 'audit_bytes': path.stat().st_size}))

if __name__ == '__main__':
    main()
