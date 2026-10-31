"""Prove exact finite source applicability to832, without relabelling producer executions."""
import collections
import difflib
import context as C
import reconstruct as R

def tree(revision,*paths):
    return {line.split('\t')[1]:line.split('\t')[0].split() for line in C.git('ls-tree','-r',revision,*paths).decode().splitlines()}

def main():
    a=tree(C.ACTUAL,'crates','assets','Cargo.toml','Cargo.lock')
    b=tree(C.TARGET,'crates','assets','Cargo.toml','Cargo.lock')
    C.need(a.keys()==b.keys(),'source-file-inventory')
    changed=[p for p in a if a[p]!=b[p]]
    C.need(changed==[C.MAIN],'actual-pristine-target-only-main',changed)
    before=C.read(C.O/'source-blobs/actual-pristine'/C.MAIN).decode();after=C.read(C.O/'source-blobs/target'/C.MAIN).decode()
    ra,rb=R.rows(before),R.rows(after);C.need(ra.keys()==rb.keys(),'source-row-domain')
    changed_rows=[];masked=after
    for name in ra:
        if ra[name]['line']==rb[name]['line']: continue
        C.need([k for k in R.FIELDS if ra[name]['fields'][k]!=rb[name]['fields'][k]]==['numerators'],'prior-delta-only-numerator',name)
        C.need(name not in C.CASES,'target-does-not-preapply-next-three',name)
        C.need(ra[name]['line'].replace(f"numerators: 0x{ra[name]['fields']['numerators']:016x}",'numerators: TOKEN')==rb[name]['line'].replace(f"numerators: 0x{rb[name]['fields']['numerators']:016x}",'numerators: TOKEN'),'prior-delta-token-only')
        changed_rows.append({'case':name,'actual_pristine_Row':ra[name]['fields'],'target_Row':rb[name]['fields']})
        masked=masked.replace(rb[name]['line'],ra[name]['line'],1)
    C.need(masked==before and len(changed_rows)==4,'all-main-bytes-four-prior-tokens')
    C.put('actual-pristine-to-target-main.patch',''.join(difflib.unified_diff(before.splitlines(keepends=True),after.splitlines(keepends=True),fromfile=C.ACTUAL+'/'+C.MAIN,tofile=C.TARGET+'/'+C.MAIN)).encode())
    selected=['crates/feynkit-generator/src/lib.rs','crates/feynkit-generator/src/grouping.rs','crates/feynkit-graph/src/lib.rs','crates/feynkit-graph/src/numerator_interface.rs','crates/feynkit-model/src/lib.rs','crates/idenso/src/tensor/composition.rs','crates/idenso/src/tensor/composition/product.rs','crates/spenso/src/structure/partial.rs','Cargo.toml','Cargo.lock']
    binds=[]
    for path in selected:
        data=C.git('show',C.TARGET+':'+path);C.put('source-blobs/target/'+path,data)
        C.need(C.blob(data)==b[path][2]==a[path][2],'selected-source-byte-identity',path)
        binds.append({'path':path,'actual_pristine':C.ACTUAL,'target':C.TARGET,'git_blob':C.blob(data),'sha256':C.sha(data),'bytes':len(data),'byte_identical':True})
    # Full path inventories bind all files/models/manifests, including inert source config
    # object IDs, without reading configuration files or executing scientific code.
    C.save('complete-source-object-bindings.json',{'actual_pristine':C.ACTUAL,'target':C.TARGET,'object_format':'Git mode/type/blob; exact immutable byte identity','paths':{p:{'actual':a[p],'target':b[p]} for p in a},'only_changed_path':C.MAIN})
    full_delta=C.git('diff','--name-status',C.ACTUAL,C.TARGET).decode()
    C.put('actual-pristine-to-target-all-paths.txt',full_delta.encode())
    code_delta=C.git('diff','--name-status',C.CODE,C.TARGET).decode()
    C.put('relevant-code-to-target-all-paths.txt',code_delta.encode())
    C.need(C.git('show',C.CODE+':'+C.MAIN).decode()==after,'6f-target-expectation-byte-identity')
    C.save('source-applicability.json',{'target_revision':C.TARGET,'relevant_code':C.CODE,'actual_capture_producer':C.ACTUAL,'actual_capture_producer_relabelled':False,'historical_capture_producer':C.HISTORICAL,'historical_kernel':'c78cffc20f318df1fce50539a587cbe8ffc717c0','current_kernel':'8eb7c3729d6abefc916decc1bdf0c06536f8af80','observer_built':'5ea149ffd88ae3afc4072a38c0871c1329be8424 + qualified overlay','proof_source':'5ba04864ca506e259d9eef4a54f65acdd8323fc2','source_files':len(a),'byte_identical_source_object_files':len(a)-1,'delta_fields':['four unrelated prior numerator literals in149-Row table'],'prior_rows_differing':changed_rows,'all_nonrow_main_bytes_identical':True,'selected_frozen_sources':binds,'model_assets':'Entire assets subtree identical by exact Git objects; complete captured model/request/options/Run and physical ownership bound separately in physical-input-bindings.json.','features':'Cargo manifests/lock exact832==actualffa; actual historical dev-optim/default+native, current ci-optim/default+native and direct-rustc observer qualifications retained from producer-source-bindings.json. No new build/profile equivalence.','external_producer_qualification':C.ref(C.T/'raw-next-qcd-captures-20261031-r1/producer-source-bindings.json'),'other_changes':'Complete all-path deltas retained. Docs/evidence and notebook stub are outside the identical crates/assets/manifests/lock domain; they are not generator or case changes.','runtime_at_target':'NOT_RUN; source applicability only. No new compile/runtime pass.'})
    C.save('applicability-consumed-inputs.json',C.used)
    print('PASS source applicability:',len(a),'files;',len(changed_rows),'unrelated tokens; no new runtime claim')

if __name__=='__main__': main()
