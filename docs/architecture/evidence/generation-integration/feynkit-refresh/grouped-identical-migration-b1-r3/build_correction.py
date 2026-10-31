"""Build a subordinate B1 record; never modifies original packet members."""
from pathlib import Path
import hashlib,json,os,subprocess,shutil,sys
W=Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh');O=Path(__file__).resolve().parent
C=W/'docs/architecture/evidence/generation-integration/feynkit-refresh/grouped-identical-migration-b1-r3'
R=W/'.local/review-grouped-identical-integration-astra-20261031-033116Z-c0749b38';OLD='7c0e13a08786416a9f24b3bec6df8a50e4592346';FIX='ed669e2beeb39f36b3e5753ab9c746770ba7653d'
ENV={k:os.environ[k] for k in ['PATH','HOME','USER','LOGNAME','LANG','LC_ALL'] if k in os.environ};ENV.update(GIT_OPTIONAL_LOCKS='0',GIT_CONFIG_NOSYSTEM='1',GIT_CONFIG_GLOBAL='/dev/null')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,d):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(d,indent=2)+'\n')
def cp(source,target):
 dest=C/target;assert Path(source).stat().st_size<=1048576;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,dest);assert sha(source)==sha(dest)
initial=load(O/'initial.json');review=load(R/'review.json');aud=load(O/('immutable-'+FIX+'-audit.json'));assert aud['matched']==425 and not aud['missing_paths'];assert sha(W/'.gitignore')==initial['gitignore_sha256'];assert sha(W/'docs/developers.toml')==initial['registry_sha256']
assert subprocess.check_output(['git','--git-dir='+str(W.parent/'main/.git'),'show','-s','--format=%an <%ae>',FIX],env=ENV,text=True).strip()=='Lucien Huber <im@lcnbr.ch>'
for r in initial['missing_paths']:assert sha(W/r['path'])==r['sha256'] and (W/r['path']).stat().st_size==r['local_bytes']
C.mkdir(parents=True,exist_ok=True)
cp(R/'review.json','review/review.json');cp(R/'seal.json','review/seal.json');cp(R/'git-manifest-missing.json','review/git-manifest-missing.json')
external=[]
for name,r in review['evidence'].items():
 path=Path(r['path']);assert sha(path)==r['sha256'] and path.stat().st_size==r['bytes'];external.append({'name':name,**r,'retained_external':True,'new_execution':False})
archive=Path(review['producer']['archive']);elf=Path(review['producer']['executed_ELF']['path']);assert sha(archive)=='ba02f91fc47a35b03cc9d51dded39836a7635bdf423c6025ab6e4775c504fc35';assert sha(elf)=='7f00172ace12f04ef53d51a70e6ac4a12248e05e3d61edbde42e5e02079f5c81'
runtime=load(R/'runtime-results.json');assert review['fresh_runtime']['passes']==39 and review['fresh_proof_join']['owner_associations']==2124
for f in review['fresh_proof_join']['public_files']:
 for role in ['fresh_path','approved_path','implementation_path']:assert sha(f[role])==f['sha256'] and Path(f[role]).stat().st_size==f['bytes']
cp(O/('immutable-'+FIX+'-audit.json'),'manifest-audit.json');cp(O/'historical-7c-audit.json','historical-7c-audit.json');cp(O/'audit.py','audit.py');cp(O/'audit-before-snapshot.py','audit-before-snapshot.py');cp(O/'build_correction.py','build_correction.py')
for p in sorted(O.iterdir()):
 if p.is_file() and (p.name.startswith(('track-exact13.','jj-track-help.','repaired-audit.','historical-audit.','pre-snapshot-audit.')) or p.name=='initial.json'):cp(p,'tracking/'+p.name)
errors=[{'attempt':'optional packet-closure.json read','kind':'read-only setup','detail':'Filename absent; reviewer-listed git-manifest-missing.json consumed instead. No tracking or data modification.'},{'attempt':'pre-snapshot immutable audit','source':'audit-before-snapshot.py','stderr':'tracking/pre-snapshot-audit.stderr','detail':'jj file track --include-ignored returned0 and queued exact files. Initial --ignore-working-copy read saw the unsnapshotted change;425-completeness assertion stopped. Explicit snapshot then produced ed669e2b; same audit source passes425 actual Git blobs. No audit condition weakened.'}]
save(C/'tracking/errors-and-source-versions.json',{'tracking_exit':load(O/'track-exact13.command.json')['exit'],'tracking_stderr_sha256':sha(O/'track-exact13.stderr'),'tracking_warning_bytes':(O/'track-exact13.stderr').stat().st_size,'errors':errors,'audit_before_snapshot_sha256':sha(O/'audit-before-snapshot.py'),'audit_corrected_workflow_sha256':sha(O/'audit.py'),'source_bytes_unchanged_workflow_corrected':sha(O/'audit-before-snapshot.py')==sha(O/'audit.py'),'temporary_ignore_exceptions':False,'gitignore_before_after_sha256':initial['gitignore_sha256'],'tracking_limit':'1MiB','include_ignored_override_bounded_to_exact13_all_under_limit':True})
save(C/'review/external-artifacts.json',{'review':initial['review'],'lowercase_seal':initial['seal'],'artifacts':external,'reviewer_runtime_credit':39,'reviewer_five_capture_credit':5,'reviewer_owner_associations':2124,'fresh_repair_runtime_credit':0,'external_runtime_archive_libraries_proof_inputs_required':True,'selfcontained_claim':False})
save(C/'missing-logs.json',{'source_review_sha256':initial['review']['sha256'],'finding':'B1','exact_paths':initial['missing_paths'],'repair':'Same path and exact bytes, explicitly tracked using supported JJ --include-ignored; no rename/reformat/truncation.'})
receipt={'schema':'grouped-identical-B1-durability-correction/1','status':'NEW_INDEPENDENT_FOCUSED_REREVIEW_REQUIRED','finding':'B1 only','author':'Lucien Huber <im@lcnbr.ch>','old_candidate':OLD,'old_closure_claim':'FALSE: local425 hashes and working-file equality did not establish committed closure. Immutable7c contains412/425 and omits13 ignored logs. Old claims remain historical and are corrected here; old7c is not relabelled complete.','review_sha256':initial['review']['sha256'],'review_lowercase_seal_sha256':initial['seal']['sha256'],'log_repair_commit':FIX,'log_repair_bookmark':'generation-grouped-identical-b1-log-repair-20261031','candidate_bookmark':'generation-grouped-identical-migration-candidate-r3-20261031','evidence_bookmark':'generation-grouped-identical-migration-evidence-r3-20261031','manifest_entries':425,'complete_immutable_log_repair_audit':{'path':'manifest-audit.json','commit':FIX,'sha256':sha(C/'manifest-audit.json'),'matched':425,'missing':[],'mismatched':[],'unexpected':[],'scoped_manifest_tree_sha256':aud['scoped_manifest_tree_sha256'],'method':aud['method']},'final_candidate_audit':'Final candidate adds only this subordinate packet on the complete log-repair commit. After freezing, repeat the same425 actual Git blob audit against its exact final commit and require the identical scoped-manifest vector. Exact final commit/audit/seals are retained in task scratch handoff.json and handoff.typ, avoiding a circular self-commit hash.','unchanged_old_manifest_sha256':'f9a8dec1335e5d23f90561c919b6e59136da023f344ee5f120217afd7a1fffdf','unchanged_old_receipt_sha256':'1f5048596ee6c2cb2f2fcbc50f0c6f64ef258e203715251538ac7e82fb9c8331','unchanged_old_seal_sha256':'5a1d9886a2dc0156267b39818ef0c2d3c3fec72d87739ac295acaa848d5e559f','source_archive_runtime_proof_invariance':{'code':'6f41807b064bcf2f579883134daef65c8ac279d3','archive':str(archive),'archive_sha256':sha(archive),'ELF':str(elf),'ELF_sha256':sha(elf),'all_preexisting_Git_path_modes_and_blobs_preserved':True,'old425_members_unchanged':True,'canonical_note_and_registry_unchanged':True,'all118_bindings_preserved':True,'accepted_math_review':'7d96505dcb700eeae42b5a4e87676551ba910cb660df586cafe34d6b90e944f1','both_modes_P1_P5_and_operand_output_maps_limits_retained':True,'old_implementation_39_and_new_review_39_remain_separate_prior_runs':True,'reviewer_passes':39,'public_captures':5,'owner_associations':2124,'repair_build_Nix_compiler_native_or_test_runs':0,'review_and_external_runtime_bytes_rehashed':True},'reconstructed_committed_lossless_bindings':aud['lossless_bindings_count'],'copied_migration_review_seal_files':len(aud['complete_copied_migration_review_seal']),'copied_preparation_manifest_entries':len(aud['complete_copied_preparation_manifest']),'ignored_scratch':str(O),'preserved_refs':initial['preserved_refs'],'review_verdict_remains':'REQUEST_CHANGES; root has not approved applied integration','self_approval':False,'Typst_rendered':False,'delegation':0,'license_or_private_config_reads':0,'new_runtime_executions':0,'CI_or_budget_or_global_integration_or_publication_change':False,'remaining_gate':'NEW independent focused B1 corrected-candidate review; all unrelated existing gates remain.'}
save(C/'receipt.json',receipt)
(C/'record.typ').write_text('''= Grouped Identical evidence durability correction
<grouped-identical-evidence-durability-correction>

#strong[NEW_INDEPENDENT_FOCUSED_REREVIEW_REQUIRED.] This subordinate record
repairs only B1 from the complete independent applied-integration review
`0f0851e068653202a7981e71b4d5ccf54316d1bf544caa0725f99cc1794066f4`.
The actual lowercase `seal.json` is bound by
`121a41e327b0c353b5fe24fae1582e239925b86024218719072dcbfd08e7d51a`.
The verdict remains REQUEST_CHANGES; root has granted no applied integration
approval. The exact review and seal are retained verbatim under `review/`.

The old claim of complete Git closure at `7c0e13a0…` was false. Its425 entries
matched local files, but only412 existed in that immutable Git tree;13 ignored
`.log` paths were absent. The original receipt, note, manifest and seal remain
unchanged historical bytes. This correction does not call old7c complete.
`historical-7c-audit.json` reproduces every existing blob and all13 omissions.

`ed669e2beeb39f36b3e5753ab9c746770ba7653d` adds exactly those13 existing logs,
at their original manifest paths with exact hashes/sizes. Documented JJ
`file track --include-ignored` selects only these paths after checking every
size against the unchanged1MiB cap. No ignore exception, main index or global
Git configuration edit is used. The first audit preceded the tracking snapshot
and correctly stopped; its source/stderr are retained. After an explicit
snapshot, the same audit verifies425 actual immutable Git blobs. Tracking
commands, stdout/stderr and corrected workflow/source versions are retained.

`manifest-audit.json` enumerates all425 workspace-relative entries, including
the canonical note and `docs/developers.toml`. It checks exact tree membership,
Git blob identity and SHA-256 from `git cat-file --batch`, independently
rejecting missing, mismatched or unexpected scoped paths. It also reconstructs
all417 lossless bindings and the complete copied migration-review/preparation
closure from committed blobs. Current-file equality and `git ls-files` are not
its membership evidence.

The final r3 candidate adds only this subordinate packet after the immutable
log repair. A post-freeze audit addresses the exact final commit and requires
its complete425-entry path/mode/blob/hash vector to equal the committed audit.
Its exact commit, audit hash and packet seal are recorded in the retained
ignored scratch `handoff.json` and `handoff.typ`; the packet does not assert a
circular self-commit hash. The r3 bookmarks identify the final candidate.

Code6f41807b, archiveba02f91f and executed ELF7f00172a are unchanged. The
independent review's39 passed Rows, five captures and2124 owner associations
remain that review's prior runtime evidence. No native/test/build/Nix/compiler
or observer/math work is rerun for this repair. Separate None/Auto P1–P5,
finite named-consumer source argument, both operand/output maps, all observer/
definition/tag-subset limits and every prior scope remain intact. Canonical
registry,118 source bindings and all already-sealed members are preserved.
Large review/runtime/proof artifacts retain their exact external paths/hashes/
sizes; no self-contained archive claim is made.

This Typst record is unrendered. Only source/JSON/hash/tree/link/diff integrity
checks apply. All original refs and partials stay; CI, resource policy, global
integration and publication remain unchanged. A new independent focused
re-review is required before root can dispose B1 or grant integration approval.
''')
# Prepare tracking records before requesting tracking; later snapshot captures final results.
for n in ['track-packet.stdout','track-packet.stderr']:(C/'tracking'/n).write_bytes(b'')
command=['jj','--config','snapshot.auto-track=none()','file','track',str(C.relative_to(W))]
save(C/'tracking/track-packet.command.json',{'argv':command,'cwd':str(W),'status':'pending'})
assert all(p.stat().st_size<=1048576 for p in C.rglob('*') if p.is_file())
r=subprocess.run(command,cwd=W,env=ENV,capture_output=True);(C/'tracking/track-packet.stdout').write_bytes(r.stdout);(C/'tracking/track-packet.stderr').write_bytes(r.stderr);save(C/'tracking/track-packet.command.json',{'argv':command,'cwd':str(W),'exit':r.returncode,'stdout_sha256':hashlib.sha256(r.stdout).hexdigest(),'stderr_sha256':hashlib.sha256(r.stderr).hexdigest(),'all_new_files_checked_before_tracking_bytes':1048576});assert r.returncode==0,(r.stdout,r.stderr)
print('PASS subordinate B1 packet prepared; original425 bytes untouched; reviewer artifacts rehashed; tracking exit0',flush=True)
