from pathlib import Path
import hashlib
import json
import os
import subprocess
import tomllib

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
P = W / '.local/feynkit-refresh/combined-candidate-20261031-r1'
E = 'docs/architecture/evidence/generation-integration/feynkit-refresh/combined-review-20261031'
F = W / '.local/feynkit-refresh/combined-integration-preflight-20261031-r1'
BASE = '7f03ab8fd5834d7c78164429e2a44dab8d01ef74'
CODE = '6f41807b064bcf2f579883134daef65c8ac279d3'
ENV = {k: os.environ[k] for k in ['PATH', 'HOME', 'USER', 'LOGNAME', 'LANG', 'LC_ALL'] if k in os.environ}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')

def sha(b):
    return hashlib.sha256(b).hexdigest()

def git(*args, input=None):
    return subprocess.check_output(['git', '--git-dir=' + str(W.parent / 'main/.git'), *args], env=ENV, input=input)

def emit(name, data):
    (W / E / name).write_text(json.dumps(data, indent=2) + '\n')

def blobs(rev, paths):
    raw = git('cat-file', '--batch', input=''.join(rev + ':' + p + '\n' for p in paths).encode())
    pos, result = 0, {}
    for p in paths:
        end = raw.index(b'\n', pos)
        header = raw[pos:end].decode().split()
        assert len(header) == 3 and header[1] == 'blob', (rev, p, header)
        size = int(header[2]); pos = end + 1
        data = raw[pos:pos + size]; pos += size + 1
        assert hashlib.sha1(b'blob %d\0' % size + data).hexdigest() == header[0]
        result[p] = data
    assert pos == len(raw)
    return result

def main():
    filters = json.loads((F / 'source-filter-bindings.json').read_text())
    paths = sorted({r['path'] for f in filters for r in f['manifest']})
    old, code, base = (blobs(rev, paths) for rev in ['2cfb52b28d99f92305a9d571bbe3911cf04b8a19', CODE, BASE])
    rows = []
    for f in filters:
        bad, now = [], []
        for r in f['manifest']:
            p = r['path']
            if sha(old[p]) != r['sha256'] or old[p] != code[p] or code[p] != base[p] or base[p] != (W / p).read_bytes():
                bad.append(p)
            now.append({'path': p, 'sha256': sha(base[p])})
        assert not bad, (f['producer'], bad)
        rows.append({'producer': f['producer'], 'filtered_source': f['src'], 'files': len(now), 'all_recorded_compile_file_hashes_equal_2cf_6f_7f_and_candidate': True, 'manifest_sha256': sha(json.dumps(now, sort_keys=True, separators=(',', ':')).encode()), 'source_owner': f['source_owner'], 'phase_inputs': f['phases_and_synthetic_inputs']})
    delta = git('diff-tree', '--no-commit-id', '--name-status', '-r', '2cfb52b28d99f92305a9d571bbe3911cf04b8a19', CODE).decode().splitlines()
    actual_non_docs = [x for x in delta if not x.split('\t')[-1].startswith('docs/')]
    assert actual_non_docs == ['M\tcrates/feynkit-generator/tests/generation_oracle/main.rs'], actual_non_docs
    mainrs = 'crates/feynkit-generator/tests/generation_oracle/main.rs'
    a, b = (git('show', r + ':' + mainrs) for r in ['2cfb52b28d99f92305a9d571bbe3911cf04b8a19', CODE])
    assert a.replace(b'0x00fac78da0e0d159', b'0xa60df38a6b8079e3') == b
    assert (W / mainrs).read_bytes() == b
    assert git('rev-parse', CODE + ':' + mainrs).decode().strip() == '7283198bf0152c82c1cd5f0173703cddfa7f5c0f'
    emit('compiled-input-invariance.json', {'schema': 'pure-compiled-input-invariance/1', 'source_roles': {'prior_core_source': '2cfb52b28d99f92305a9d571bbe3911cf04b8a19', 'actual_grouped_source': CODE, 'accepted_evidence_base': BASE}, 'filters': rows, 'unique_files': len(paths), 'non_documentation_delta_2cf_to6f': actual_non_docs, 'only_generator_literal_outside_13_filters': True, 'candidate_production_delta_from7f': [], 'new_Nix_evaluation_build_runtime_or_native_runs': 0, 'scope': 'Exact recorded compile source membership and bytes only. Does not re-evaluate derivations, compiler synthetic manifests, runtime discovery/tool closure or docs runtime catalogs. Producer sibling ACTIVE; zero positive runtime credit.'})
    reg = tomllib.loads((W / 'docs/developers.toml').read_text())
    scopes = []
    for section in reg['section']:
        for note in section.get('note', []):
            for s in note.get('scope', []):
                source = (W / s['path']).read_text()
                assert source.count(s['anchor']) == 1
                if 'end_anchor' in s:
                    assert source.count(s['end_anchor']) == 1
                    start, end = source.index(s['anchor']), source.index(s['end_anchor'])
                    assert end >= start + len(s['anchor'])
                    reviewed = source[start:end].encode()
                else:
                    reviewed = source.encode()
                actual = hashlib.sha1(b'blob %d\0' % len(reviewed) + reviewed).hexdigest()
                assert actual == s['digest'], (note['id'], s['path'], actual)
                scopes.append({'note': note['id'], **s, 'actual_digest': actual})
    assert len(scopes) == 118
    emit('source-scope-validation.json', {'canonical_owner': 'crates/alphal00p-docs-builder/src/lib.rs: DeveloperScope extraction; unique start/end anchors, Git blob digest of exact range or full file', 'ordered_scopes': scopes, 'count': 118, 'project_code_executed': False})
    # Retain original failed-name data as data, never infer a new package run.
    full = W / '.local/feynkit-refresh/full-generator-characterization-20261031-r1'
    for name in ['characterization.json', 'failure-matrix.json']:
        data = (full / name).read_bytes()
        assert len(data) <= 1048576
        (W / E / ('prior-full-' + name)).write_bytes(data)
    hist = json.loads((W / E / 'historical-manifest-bindings.json').read_text())
    for m in hist['manifests']:
        m['original_manifest_evidence_revision'] = '7c0e13a08786416a9f24b3bec6df8a50e4592346' if m['manifest'].endswith('/grouped-identical-migration/SHA256SUMS') else BASE
        m['repaired_complete_evidence_revision'] = BASE
        m['old7c_complete_claim'] = False
    emit('historical-manifest-bindings.json', hist)
    print(json.dumps({'filters_verified': len(rows), 'unique_compile_files': len(paths), 'scope_digests': len(scopes), 'new_runtime_credit': 0}))

if __name__ == '__main__':
    main()
