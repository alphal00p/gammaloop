from pathlib import Path
import gzip
import hashlib
import json
import subprocess

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
P = W / '.local/feynkit-refresh/combined-candidate-20261031-r1'
E = W / 'docs/architecture/evidence/generation-integration/feynkit-refresh/combined-review-20261031'
META = '6a12f775c8671194dfb05f4e382c3064575c9186'

def sha(b):
    return hashlib.sha256(b).hexdigest()

def main():
    freeze = E / 'freeze'
    freeze.mkdir(exist_ok=False)
    # A historical rejected review is retained as evidence, never called accepted.
    execution = W / 'docs/architecture/evidence/generation-integration/execution.json'
    for path in [execution, E / 'receipt.json']:
        data = json.loads(path.read_text())
        state = data['current_status_reconciliation_20261031'] if path == execution else data
        state['review_records'] = state.pop('accepted_reviews')
        state['historical_rejected_review_ids'] = ['grouped-applied-original']
        state['review_record_semantics'] = 'Verbatim review index includes approvals within their scopes AND original0f0851e REQUEST_CHANGES. Only latest66db3d6d/root disposition accepts the corrected entire bounded applied grouped conjunction at7f. No combined review approval.'
        path.write_text(json.dumps(data, indent=2) + '\n')
    source = P / ('immutable-' + META + '-audit.json')
    raw = source.read_bytes()
    encoded = gzip.compress(raw, mtime=0)
    parts = []
    for i, start in enumerate(range(0, len(encoded), 524288)):
        chunk = encoded[start:start + 524288]
        target = freeze / ('audited-metadata-tree.' + str(i).zfill(4) + '.json.gz.part')
        target.write_bytes(chunk)
        parts.append({'path': str(target.relative_to(W)), 'sha256': sha(chunk), 'bytes': len(chunk)})
    record = {
        'schema': 'combined-refresh-freeze-record/1',
        'audited_metadata_revision': META,
        'audited_metadata_tree': json.loads(raw)['tree'],
        'audited_metadata_vector_sha256': json.loads(raw)['full_tree_vector_sha256'],
        'audit_encoding': 'lossless gzip chunks, concatenated in order',
        'audit_sha256': sha(raw), 'audit_bytes': len(raw), 'audit_parts': parts,
        'new_candidate_delta_from_audited_metadata': 'This freeze packet, receipt/execution schema clarification distinguishing the historical REQUEST_CHANGES record from scoped accepted reviews, and corresponding combined seal refresh. No other canonical or production change.',
        'status': 'NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED',
        'source': '6f41807b064bcf2f579883134daef65c8ac279d3',
        'accepted_evidence_base': '7f03ab8fd5834d7c78164429e2a44dab8d01ef74',
        'final_candidate_audit': 'Exact final immutable candidate/tree audit is performed after snapshot and retained as a detached sealed handoff packet; this parent audit does not claim to include future tree bytes.',
        'main_index_old_refs_CI_production_preserved': True,
        'global_integration': 'fdff0d18953348b84a0f1b7e6ab5ef703a5090de',
        'new_native_build_or_runtime_credit': 0,
        'independent_review': 'root to assign; no selfapproval',
        'remaining_gates': '../current-status-matrix.json',
    }
    (freeze / 'record.json').write_text(json.dumps(record, indent=2) + '\n')
    (freeze / 'freeze.py').write_bytes(Path(__file__).read_bytes())
    (freeze / 'metadata-audit.stdout').write_bytes((P / 'audit-metadata-corrected.stdout').read_bytes())
    (freeze / 'metadata-audit.stderr').write_bytes((P / 'audit-metadata-corrected.stderr').read_bytes())
    (freeze / 'guards-before.json').write_bytes((P / 'guards-before.json').read_bytes())
    # Preserve the original assembly and corrected metadata revisions under separate bookmarks.
    paths = [str(x.relative_to(W)) for x in E.rglob('*') if x.is_file() and not x.is_relative_to(freeze) and x not in [E / 'SHA256SUMS', E / 'SEAL.json']]
    paths += ['docs/architecture/generation-combined-refresh-review.typ', 'docs/architecture/generation-feynkit-community-refresh.typ', 'docs/architecture/evidence/generation-integration/execution.json', 'docs/developers.toml']
    text = ''.join(sha((W / name).read_bytes()) + '  ' + name + '\n' for name in sorted(paths))
    (E / 'SHA256SUMS').write_text(text)
    seal = json.loads((E / 'SEAL.json').read_text())
    seal.update(manifest_sha256=sha(text.encode()), files=len(paths))
    (E / 'SEAL.json').write_text(json.dumps(seal, indent=2) + '\n')
    files = sorted(x for x in freeze.rglob('*') if x.is_file())
    assert all(x.stat().st_size <= 1048576 for x in files)
    text = ''.join(sha(x.read_bytes()) + '  ' + str(x.relative_to(W)) + '\n' for x in files)
    (freeze / 'SHA256SUMS').write_text(text)
    seal = {'schema': 'combined-refresh-freeze-seal/1', 'manifest': 'SHA256SUMS', 'manifest_sha256': sha(text.encode()), 'files': len(files), 'status': 'NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED', 'audited_parent_metadata': META, 'approval': False}
    (freeze / 'SEAL.json').write_text(json.dumps(seal, indent=2) + '\n')
    files += [freeze / 'SHA256SUMS', freeze / 'SEAL.json']
    subprocess.run(['jj', '--config', 'snapshot.auto-track="none()"', 'file', 'track', '--include-ignored', *[str(x.relative_to(W)) for x in files]], cwd=W, check=True)
    print(json.dumps({'new_freeze_files': len(files), 'combined_seal_sha256': sha((E / 'SEAL.json').read_bytes()), 'freeze_seal_sha256': sha((freeze / 'SEAL.json').read_bytes())}))

if __name__ == '__main__':
    main()
