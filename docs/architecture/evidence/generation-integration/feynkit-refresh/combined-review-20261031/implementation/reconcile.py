from pathlib import Path
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
GLOBAL = 'fdff0d18953348b84a0f1b7e6ab5ef703a5090de'
ENV = {k: os.environ[k] for k in ['PATH', 'HOME', 'USER', 'LOGNAME', 'LANG', 'LC_ALL'] if k in os.environ}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')

def sha(b):
    return hashlib.sha256(b).hexdigest()

def git(*args):
    return subprocess.check_output(['git', '--git-dir=' + str(W.parent / 'main/.git'), *args], env=ENV)

def emit(path, value):
    (W / path).write_text(json.dumps(value, indent=2) + '\n')

def main():
    execution = 'docs/architecture/evidence/generation-integration/execution.json'
    doc = 'docs/architecture/generation-feynkit-community-refresh.typ'
    registry = 'docs/developers.toml'
    original = json.loads(git('show', BASE + ':' + execution))
    current = json.loads((W / execution).read_text())
    assert current == original
    patch = json.loads((F / 'execution-reconciliation.NOT-APPLIED.json-patch').read_text())
    for op in patch:
        if op['path'] == '/current_status_reconciliation_20261031':
            continue
        parts = op['path'].strip('/').split('/')
        owner = current
        for part in parts[:-1]:
            owner = owner[part]
        if op['op'] == 'test':
            assert owner[parts[-1]] == op['value']
        else:
            assert op['op'] == 'replace'
            owner[parts[-1]] = op['value']
    current['packages']['WP3']['current_slice_review'] = 'Source/evidence frozen on accepted grouped7f/source6f for NEW independent combined refresh review. Constituent acceptances retained; no combined approval or WP3/M-rebuild closure.'
    current['active_xs_qed_l2']['status'] = 'Original fdff/eae bounded integration accepted through independent scopes A/B;298/73/49 remains old1a history. PRIOR8eb ffa full-pool1 characterization is302/69/49. Later accepted grouped migration leaves68 other failure identities pending case proof/disposition; no fresh303/68/49 run.'
    accepted = json.loads((W / E / 'review-index.json').read_text())
    state = {
        'schema': 'generation-combined-refresh-candidate/1',
        'status': 'NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED',
        'accepted_evidence_base': BASE, 'actual_grouped_source': CODE,
        'earlier_host_source': '2cfb52b28d99f92305a9d571bbe3911cf04b8a19',
        'earlier_host_evidence': '457e361d7d39e0a3daae801651f8c1dfabbf47aa',
        'original_approved_global_integration': GLOBAL,
        'source_upstream': '9102f885722e856e30227ed46e685e5c465905f7',
        'kernel': '8eb7c3729d6abefc916decc1bdf0c06536f8af80',
        'author': 'Lucien Huber <im@lcnbr.ch>',
        'CI_enabled': False, 'publication_authorized': False, 'global_advance_performed': False,
        'accepted_reviews': accepted['reviews'],
        'combined_review': 'Unassigned: root must obtain a NEW independent review of this exact conjunction; constituent approvals alone do not approve it.',
        'bounded_grouped': {
            'root_disposition': E + '/root-dispositions/root-grouped-applied-integration-acceptance-20261031.json',
            'scope': 'ENTIRE corrected bounded applied grouped conjunction at7f, not durability alone',
            'archive_sha256': 'ba02f91fc47a35b03cc9d51dded39836a7635bdf423c6025ab6e4775c504fc35',
            'oracle_elf_sha256': '7f00172ace12f04ef53d51a70e6ac4a12248e05e3d61edbde42e5e02079f5c81',
            'migrated_rows': 13, 'actual_runtime_passes': 39, 'pools': [1, 8, 32],
            'public_files': 5, 'owner_associations': 2124,
            'full_Row': 'Only literal Classes::None; Auto physical captures and captured search-counter premises separately retained. Pools8/32 are runtime evidence only.',
            'proof': 'Independent current raw30 plus separate None/Auto P1-P5. Both operand AND output-axis maps; finite current named logical consumers with validated canonical p/p-inverse. No historical grouped producer upgrade, dense/positional/general decoder/import/locality/full-integrand theorem.',
            'historical_B1': 'Old7c review0f0851e remains REQUEST_CHANGES:412/425 closure and13 missing logs. Repaired7f425/425 plus29 subordinate entries accepted by66db3d6d. Prior positive39freshchecks reused by correction; zero additional native runs.'
        },
        'prior_full_generator_pool1': {
            'source': 'ffa76ec420642c03349ab756d2d7720ab93c0057',
            'inventory': 420, 'selected': 371, 'passed': 302, 'failed': 69, 'ignored': 49,
            'fresh': False, 'meaning': 'Actual prior characterization; numerator-only asserted Summary pairs are not physical proof.',
            'derived_pending_other_identities': 68,
            'derivation': 'One of prior69 failure identities now has an accepted bounded grouped migration. This is a worklist subtraction, never a newly executed303/68/49 result.',
            'complete_pool8_or32': False
        },
        'prior_twelve_rows': 'Accepted raw-three integration535b7926 on SM/prior-eight prerequisites:36 actual full-Row None checks at1/8/32 and15 pool1 physical joins; mode-specific physical Auto and observer/whole-definition/complete operand-output maps retained.',
        'host_scope': 'Accepted applied stub2cf, corrected ordinary wheel and separate private RustRed7c1/kernel8eb host with8fresh independent passes. Ordinary lacksIBPFamily. Private arity1 only; tracing_max_level_info differs ordinary. Loader, wheel-vs-loose-companion and external closure qualifications retained. No production RustRed/fullCommunity.',
        'finite_runner': 'Exact finite r3 F1/F2 guards accepted under b370752/root disposition. Old verbatim pending wording does not reopen them. No generic runner/adapter/native/math approval.',
        'source46': 'Cumulative79dfeaab approval is real and source/profile/kernel-bound. No repeat implementation. New8eb consumer qualification is separate.',
        'production_core_API_evaluator': {
            'status': 'ACTIVE_SIBLING_NO_COMPLETION_OR_RUNTIME_RESULTS_SUPPLIED',
            'source': '2cfb52b28d99f92305a9d571bbe3911cf04b8a19',
            'producer_roots': 6,
            'source_filter_transition': E + '/compiled-input-invariance.json',
            'positive_runtime_credit_here': 0,
            'limits': 'No producer completion, current native discovery counts or runtime passes inferred. Source/filter byte comparison does not replace exact derivation or runtime qualification.'
        },
        'next_three_QCD': {
            'capture_review_sha256': '7ccf8f133bff073e5f428d1fabc0ebd6552aa1e04a8b7edd149f526b37330aa7',
            'capture_scope': '12invocations;6oldpass/6current101;24worlds only',
            'preparation_review_sha256': '62dade9062105a077ae6089d7e261742d57d63d9d07737b020c9b3098a71db95',
            'preparation_status': 'REQUEST_CHANGES: PREP1 fold endpoints; PREP2 live case/mode. r2 repair/supervisor ACTIVE, no completion supplied.',
            'math_interface_native21_migration_approval': False
        },
        'remaining_gates': E + '/current-status-matrix.json',
        'WP3': 'incomplete', 'M_rebuild': 'not claimed',
        'all_numerator_migrations_blanket_prerequisite': False,
        'proposed_local_bookmark_advance': {'bookmark': 'generation-integration', 'from': GLOBAL, 'to': 'this exact frozen combined candidate, only after NEW independent combined review and root scoped disposition', 'performed': False, 'publication': 'separate explicit user-go gate; no push/fetch/PR/upload/CI enable'},
        'historical_metadata': E + '/historical-snapshots.json',
        'scope_bindings_preserved': 118,
        'new_runtime_build_or_native_executions': 0,
        'no_self_contained_native_proof_archive_claim': True
    }
    current['current_status_reconciliation_20261031'] = state
    emit(execution, current)
    emit(E + '/receipt.json', state)
    # Preserve all ordered scope records, including every grouped-writer addition.
    reg = (W / registry).read_text()
    old_scopes = [(n['id'], s) for section in tomllib.loads(reg)['section'] for n in section.get('note', []) for s in n.get('scope', [])]
    replacements = {
        'generation-grouped-identical-migration': {
            'kind': 'Bounded migration qualification record',
            'status': 'Accepted entire bounded applied grouped conjunction at7f; broader combined refresh review required',
            'review_ref': BASE, 'reviewed_at': '2026-10-31', 'evidence_revision': BASE
        },
        'generation-host-stub-rustred': {'status': 'Accepted bounded applied stub, corrected ordinary distribution and private RustRed integration; no production RustRed/fullCommunity'},
        'generation-raw-three-migration': {'status': 'Accepted bounded twelve-row/raw-three integration; broader combined refresh review required'},
        'generation-sm-eeww-migration': {'status': 'Accepted bounded SM migration prerequisite; no whole-WP3 approval'},
        'generation-feynkit-community-refresh': {
            'summary': 'Separately accepted source, provenance, finite migration, pool1 and host evidence; new combined refresh review remains required.',
            'status': 'Bounded constituents accepted; NEW independent combined review and broader WP3 gates remain open',
            'captured_at': '2026-10-31', 'evidence_revision': BASE
        }
    }
    blocks = reg.split('[[section.note]]')
    for i, block in enumerate(blocks):
        for ident, fields in replacements.items():
            if ('\nid = "' + ident + '"\n') not in block:
                continue
            for key, value in fields.items():
                line = key + ' = ' + json.dumps(value)
                lines = block.splitlines()
                matches = [j for j, old in enumerate(lines) if old.startswith(key + ' = ')]
                if matches:
                    assert len(matches) == 1
                    lines[matches[0]] = line
                else:
                    lines.insert(2, line)
                block = '\n'.join(lines) + '\n'
            blocks[i] = block
    reg = '[[section.note]]'.join(blocks)
    marker = '[[section.note]]\nid = "generation-feynkit-community-refresh"'
    new = '''[[section.note]]
id = "generation-combined-refresh-review"
title = "Combined FeynKit refresh review candidate"
summary = "Exact7f source/evidence conjunction with accepted bounded constituents, revision-aware historical closure and explicit remaining gates."
source = "docs/architecture/generation-combined-refresh-review.typ"
kind = "Combined integration review candidate"
status = "NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED"
lifecycle = "investigation"
owner = "unassigned"
captured_at = "2026-10-31"
evidence_revision = "7f03ab8fd5834d7c78164429e2a44dab8d01ef74"
products = ["gammaloop", "feynkit", "linnet", "spenso"]
topics = ["generation", "integration", "community", "validation"]
review_triggers = ["docs/architecture/evidence/generation-integration/execution.json", "docs/architecture/generation-feynkit-community-refresh.typ", "docs/architecture/generation-combined-refresh-review.typ", "docs/architecture/evidence/generation-integration/feynkit-refresh/combined-review-20261031/**"]

'''
    assert marker in reg
    reg = reg.replace(marker, new + marker, 1)
    after_scopes = [(n['id'], s) for section in tomllib.loads(reg)['section'] for n in section.get('note', []) for s in n.get('scope', [])]
    assert len(old_scopes) == 118 and old_scopes == after_scopes
    (W / registry).write_text(reg)
    text = (W / doc).read_text()
    # Label the old account as history; preserve all original review/proof qualifications.
    text = text.replace('This updated\nnote and receipt remain an unreviewed metadata candidate. The registry uses\n#raw("evidence_revision") at #raw("cacd5182704082df929a066548e3501be22c9533"),', 'That attribution was independently accepted at metadata48d0; the new\ncombined candidate below still requires independent review. The historical\nregistry used #raw("evidence_revision") at #raw("cacd5182704082df929a066548e3501be22c9533"),')
    text = text.replace('This record and its registry note are a separate unreviewed evidence candidate.\nA new combined review remains necessary before integration approval.', 'The following review accounts describe their original bounded state; current\ndispositions are reconciled below without rewriting the sealed originals.\nA new combined review remains necessary before integration approval.')
    text = text.replace('== Separately bounded review dispositions', '== Historical separately bounded review dispositions', 1)
    addition = '''

== Current combined review candidate, 2026-10-31

#strong[NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED.] The
#link("generation-combined-refresh-review.typ")[combined candidate record]
and #link("evidence/generation-integration/feynkit-refresh/combined-review-20261031/receipt.json")[machine-readable receipt]
reconcile the exact accepted evidence base
`7f03ab8fd5834d7c78164429e2a44dab8d01ef74`, actual grouped source
`6f41807b064bcf2f579883134daef65c8ac279d3`, and current metadata.
Metadata48d0 acceptance concerns its original attribution only. SourceR4,
shared Scope A, the SM/prior-eight and raw-three prerequisites, separate
None/Auto aggregates, finite runner r3, host distributions and the entire
corrected bounded applied grouped conjunction retain their own scopes.

The latest grouped review66db3d6d accepts the whole bounded conjunction at7f:
13 migrated Rows,39 actual passes at1/8/32,5 public files and2124 owner
associations under the retained finite premises. Its correction reuses the
original independent39 checks and adds zero native runs. Old7c remains
412/425 with13 absent logs and review0f0851e remains REQUEST_CHANGES.
Exact7f closes425/425 plus29 subordinate entries. Both operand AND output maps,
whole definitions and observer premises remain mandatory in each mode; no
historical grouped-producer or universal layout/decoder/locality theorem follows.

The full pool1 run remains PRIOR ffa evidence:420 inventory,371 selected,
302 passes,69 failures,49 ignored. One failed identity now has an accepted
grouped migration;68 other identities remain pending case proof/disposition.
This derived worklist is not a fresh303/68/49 result. Only None asserts the
literal full Row; physical Auto and captured counter premises remain separate.
Ordinary corrected host still lacks IBPFamily. Private RustRed7c1/8eb arity1
and tracing_max_level_info differ from ordinary; loader, wheel/loose-companion
and external closure limits remain. No production RustRed/full Community.

Source46 cumulative79dfeaab approval stands within its original source/profile;
8eb consumer qualification is separate. Exact finite runner r3 F1/F2 guards
are accepted, despite old verbatim pending phrases. The six-producer2cf sibling
and source-filter transition qualification remain active with no supplied
completion/runtime results and zero positive runtime credit here. Next-three
QCD captures are accepted only as captures; preparation62dade90 requested
PREP1 fold-endpoint and PREP2 live-case/mode repairs, now active. No math,
interface, native21 or migration approval is inferred.

Constituent approvals do not approve this conjunction. Outstanding numerator
migrations are not a blanket prerequisite for each bounded refresh review.
WP3/M-rebuild, full fast/oracle1/8/32,3150/4664/306, fresh a_ddx/BRST/NLO/CoolQM,
R/F, fullDocs/Clippy, correctness/performance/memory and later plan stages stay
open as indexed. No20GiB or1.2x cap was adopted. A proposed local
generation-integration advance fromfdff requires NEW combined review and root
disposition; it is distinct from publication, which requires explicit user go.
Neither advance nor CI enable/publication is performed here.

Historical shared metadata is retained at exact revision-aware snapshots;
old seals continue to verify their original evidence revisions, never today's
new registry/note bytes. All118 ordered source bindings remain unchanged.
'''
    (W / doc).write_text(text + addition)
    matrix = json.loads((F / 'status-matrix.json').read_text())
    overrides = {
        'NEXT bounded combined refresh integration': ('NEW_INDEPENDENT_COMBINED_REVIEW_REQUIRED on accepted7f/source6f plus this metadata', 'No conjunction approval. Proposed local advance is conditional on NEW review/root disposition, distinct from explicit user-go publication.'),
        'current12 rows / full generator': ('Prior12 accepted; latest13-row grouped conjunction accepted separately', 'Prior302/69/49 characterization stays unchanged. One accepted grouped migration leaves68 other identities as a derived worklist, never a fresh303/68/49 run. No completepool8/32 or physics from Summary.'),
        'grouped None / Auto and finite runner': ('Separate None/Auto aggregates, finite r3 F1/F2 guards and entire corrected applied grouped conjunction at7f ACCEPTED within exact domains', 'Latest66db3d6d reuses original0f0851e positive39checks but old7c remains REQUEST_CHANGES and412/425.7f425/425+29. No historical grouped producer upgrade; both operand/output maps and P1-P5 remain mandatory.'),
        'source46 integration': ('Accepted cumulative79dfeaab, original source/profile-bound', 'Implementation accepted; no repeat needed.8eb consumer qualification separate.'),
    }
    for row in matrix:
        row['preflight_status'] = row['status']
        row['preflight_limits'] = row['limits']
        if row['gate'] in overrides:
            row['status'], row['limits'] = overrides[row['gate']]
    emit(E + '/current-status-matrix.json', {'preflight_verbatim': E + '/preflight/status-matrix.json', 'rows': matrix, 'unresolved_in_progress': {'production': state['production_core_API_evaluator'], 'QCD': state['next_three_QCD']}})
    emit(E + '/metadata-reconciliation.json', {'base': BASE, 'proposal_base': '457e361d7d39e0a3daae801651f8c1dfabbf47aa', 'method': 'Test original execution replacements against7f; update latest grouped/current scope explicitly; preserve additive grouped-writer registry and all original ordered scopes.', 'ordered_source_bindings': after_scopes, 'count': len(after_scopes), 'review_dates': 'Only actual grouped review completed2026-10-31 used for grouped registry. New combined/refresh metadata has no invented review date or approval.', 'original_execution_histories_preserved_except_explicit_status_reconciliation': True})
    print(json.dumps({'status': state['status'], 'ordered_scopes': len(after_scopes), 'tracked_metadata_paths': [execution, doc, registry]}))

if __name__ == '__main__':
    main()
