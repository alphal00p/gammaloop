"""Pure retained-data assembly and seals; no new native/runtime/compiler work."""
import collections
import datetime
import functools
import json
import os
from pathlib import Path
import sys
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
from qualification import W,O,R1,ASSIGN,PRIVATE,load,save,sha

@functools.lru_cache(maxsize=None)
def file_reference(path):
    return {'path':str(path),'sha256':sha(path),'bytes':Path(path).stat().st_size}

def reference(path,pointer=None):
    result=dict(file_reference(str(path)))
    if pointer is not None:result['json_pointer']=pointer
    return result

assert not PRIVATE.intersection(os.environ)
before=load(P/'inputs-before.json')
after={p:{'sha256':sha(p),'bytes':Path(p).stat().st_size} for p in before}
save('inputs-after.json',after)
unchanged=before==after
actual=load(P/'actual-invocations.json') if (P/'actual-invocations.json').is_file() else []
plan=load(P/'planned-invocations.json')
outcomes=[]
for row in actual:
    path=O/'execution/outcomes'/(row['id']+'.json')
    receipt=load(path) if path.is_file() else {}
    outcomes.append({'id':row['id'],'status':receipt.get('status','NO_PROTECTED_OUTCOME'),
                     'native_launched':receipt.get('native_launched',False),
                     'authenticated_execution':receipt.get('authenticated_execution',False),
                     'validated_predicate':receipt.get('validated_predicate',False),
                     'semantic_negative_credit':receipt.get('semantic_negative_credit',False),
                     'orchestration_reaudit':row['validated'],
                     'error':receipt.get('error'),
                     'outcome':reference(path) if path.is_file() else None,
                     'protected_receipt':reference(O/'execution/records'/row['id']/'result.json') if (O/'execution/records'/row['id']/'result.json').is_file() else None})
save('outcome-index.json',outcomes)
counts={'planned':21,'attempted_controller_commands':len(actual),
        'native_launches':sum(r['native_launched'] for r in outcomes),
        'authenticated_executions':sum(r['authenticated_execution'] for r in outcomes),
        'baseline_passes':sum(r['status']=='BASELINE_PASS' and r['validated_predicate'] and r['orchestration_reaudit'] for r in outcomes),
        'expected_semantic_rejections':sum(r['semantic_negative_credit'] and r['orchestration_reaudit'] for r in outcomes),
        'infrastructure_or_unvalidated_outcomes':sum(not r['validated_predicate'] or not r['orchestration_reaudit'] for r in outcomes),
        'unattempted':21-len(actual),'retries':0,
        'analysis_tool_errors_before_native':1}
execution_files={str(p):reference(p) for p in sorted((O/'execution').rglob('*')) if p.is_file()} if (O/'execution').exists() else {}
save('execution-artifacts.json',execution_files)

# Every result index is an actual protected native array member. Existing
# definition/owner associations are retained as premises, never new obligations.
inverse=[]
for batch in load(R1/'prepared/native-math-batches.json'):
    source=R1/'prepared/math'/batch['batch']/'obligations.json'
    obligations=load(source)
    for kernel in batch['kernels']:
        ident=batch['batch']+'--'+kernel
        result=O/'execution/records'/ident/'results.json'
        if not result.is_file():continue
        rows=load(result)
        assert len(rows)==len(obligations)==batch['unique_pairs']
        for index,(obligation,native) in enumerate(zip(obligations,rows)):
            for ai,association in enumerate(obligation['obligations']):
                for member in association.get('owners',[None]):
                    inverse.append({'batch':batch['batch'],'kernel':kernel,'row':index,
                                    'result':reference(result,'/'+str(index)),
                                    'obligation':reference(source,'/'+str(index)+'/obligations/'+str(ai)),
                                    'member':member,'obligation_kind':association.get('obligation'),
                                    'native_predicates':native,'definition_order':reference(R1/'prepared/math'/batch['batch']/'definition-order.json'),
                                    'credit':'Association to one native unique-pair row, not an additional fresh scientific obligation.'})
save('math-member-to-native-row.json',inverse)

typed=O/'execution/records/all-six-cells--typed8eb/results.json'
typed_associations=[]
if typed.is_file():
    native=load(typed)
    inputs=load(R1/'prepared/interface-inputs.json')
    for ci,(cell,source_cell) in enumerate(zip(native['cells'],inputs['cells'])):
        for gi,(graph,source_graph) in enumerate(zip(cell.get('graphs',[]),source_cell['graphs'])):
            base=f'/cells/{ci}/graphs/{gi}'
            typed_associations.append({'case':cell['case'],'mode':cell['mode'],'id':graph['id'],
                'full_native_graph':reference(typed,base),'exact_supplied_graph':reference(R1/'prepared/interface-inputs.json',base),
                'owners':[{'owner':owner['owner'],'actual_import_export_flags_operand_maps_and_lowering':reference(typed,base+'/owners/'+str(oi))} for oi,owner in enumerate(graph['owners'])],
                'actual_pair_plans':reference(typed,base+'/plans'),
                'actual_fold_plans_survivors_canonical_layout_and_output_axes':reference(typed,base+'/actual_fold'),
                'public_full_or_empty_numerator_of':reference(typed,base+'/public_numerator_of'),
                'prepared_operand_and_output_axis_maps':reference(R1/'prepared/all-layout-cells.json'),
                'prepared_finite_locality':reference(R1/'prepared/finite-locality-and-mass-inputs.json')})
save('typed-native-associations.json',typed_associations)
premises={name:reference(R1/'prepared'/name) for name in ['all-local-relations.json','all-layout-cells.json','all-aggregate-relations.json','consumer-obligations.json','finite-locality-and-mass-inputs.json','pristine-observer-census-associations.json','complete-key-multiplicity-domains.json','gggg-even-cycle-witnesses.json','native-math-batches.json','interface-definition-order.json']}
save('prepared-data-associations.json',premises)
limits=[
 'NEW independent scientific/finite-interface review is required; this execution packet supplies no final approval or migration permission.',
 'Exactly three QCD cases; None/Auto separately; GroupingNone; captured pool1. OnlyNone historical Row assertion; Auto has own inputs. Factor2 remains4of5 sym diagrams.',
 '600relations=578literal+22even-f-cycle/144rewrites;96unexpandedproducts;48cells;1890ownerpairinputs;912finite-locality records are finite premises, not a general decoder/locality/UV/fullintegrand/component-physics theorem.',
 '166/118 unique-pair whole-definition/ordered-side partitions per c78/8eb kernel. Repeated definitions and literal identities are not fresh scientific obligations. No old onePI cross-credit.',
 'Typed8eb imports historical/current data; no new c78 consumer execution. Both operand AND output-axis maps remain. Rank8num+rank4projector leaves FOURcoad8 outputs; no groupedrank-zero transfer.',
 '210 changed-owner lowerings require actual output;390 unchanged cases remain conditional identical-input deterministic-source arguments.26built-inbis no-op controls get zero rejection credit.',
 'Historical2b/c78Symbolica/Numerica+registryGraphica3.0.1; currentpristineffa/archive8993/kernel8eb; observerbuilt5ea+overlay differs from proofsource5ba. Dev/CI/profile/observer/callbackuserdata gaps retain identity. Raw/cooked/decoded association remains a captured premise.',
 'Authentication is exact cooperative native supervisor/receipt consistency, not cryptographic/hostile attestation. Outer PID is bwrap launcher; immediate bridge/native private delivery has source and reached-native evidence, no secret or inner PID dump.',
 'No source/expectation/metadata/ref/index changes, generator/capture/full371/host/old-batch tests, Nix/build/compiler/network-host calls, publication or delegation. Globalintegration832 and all old evidence remain untouched.'
]
completion=load(P/'phase-completion.json') if (P/'phase-completion.json').exists() else {'status':'NOT_STARTED'}
errors=[]
for row in outcomes:
    if row['error']:errors.append({'id':row['id'],'error':row['error']})
errors.append({'phase':'pre_native_static_analysis_tool','type':'false missing dependency scan','detail':'Initial depth-first scan omitted already loaded library reuse; original source/stdout/stderr/FAIL report preserved. Corrected separate-exec breadth-first scan and actual unlicensed probes provide qualification; no card/supervisor or mounts were changed.'})
save('source-and-setup-errors.json',errors)
report={'schema':'next-qcd-exact21-native-execution-candidate/1',
        'status':'READY_FOR_NEW_INDEPENDENT_SCIENTIFIC_FINITE_INTERFACE_REVIEW' if counts['baseline_passes']==5 and counts['expected_semantic_rejections']==16 and unchanged else 'STOPPED_BOUNDARY_FOR_ROOT_DISPOSITION',
        'utc':datetime.datetime.now(datetime.UTC).isoformat(),'workspace':str(W),'packet':str(P),
        'root_assignment':reference(ASSIGN),'root_acceptance':reference(load(ASSIGN)['root_supervisor_acceptance']),
        'R3_seal':reference(O/'SEAL.json'),'independent_R3_review':reference(load(ASSIGN)['supervisor_review']['path']),
        'counts':counts,'completion':completion,'all_preserved_input_bytes_unchanged':unchanged,
        'input_files':len(before),'execution_artifact_files':len(execution_files),
        'genuine_preflight':reference(P/'runtime-qualification.json') if (P/'runtime-qualification.json').exists() else None,
        'actual_outcomes':'outcome-index.json','actual_invocations':'actual-invocations.json',
        'full_execution_artifacts':'execution-artifacts.json','full_math_inverse_references':'math-member-to-native-row.json',
        'full_typed_associations':'typed-native-associations.json','prepared_associations':'prepared-data-associations.json',
        'full_before_after':['inputs-before.json','inputs-after.json'],'errors':'source-and-setup-errors.json',
        'unresolved_obligations':limits,'final_scientific_approval':False,
        'retention':'Retain NEW ignored native execution packet and actual R3 execution trees at their actual locations for independent review/root continuation. No owned background process retained.'}
save('report.json',report)
text='= Exact next-QCD native execution candidate\n\n'
text+='Status: `'+report['status']+'`. Final independent scientific/finite-interface approval remains separate.\n\n'
text+=f"Actual controller attempts: {counts['attempted_controller_commands']}/21; authenticated native executions: {counts['authenticated_executions']}; validated baselines: {counts['baseline_passes']}/5; exact semantic rejections: {counts['expected_semantic_rejections']}/16; infrastructure/unvalidated outcomes: {counts['infrastructure_or_unvalidated_outcomes']}; unattempted: {counts['unattempted']}; retries:0.\n\n"
text+=f"All {len(before)} protected old/source/input/review/assignment files unchanged: {str(unchanged).lower()}. Full output arrays/streams and partial work stay at actual R3 execution locations, externally bound by execution-artifacts.json.\n\n"
text+='The initial static analysis tool false alarm is retained and separately accounted in source-and-setup-errors.json. Actual loader and stdlib probes ran without private names or scientific helper execution.\n\n'
for limit in limits:text+='- '+limit+'\n'
text+='\nRetain this ignored packet and actual R3 execution trees for review/continuation; no source, tracked, JJ, ref or index operation was performed.\n'
with (P/'report.typ').open('x') as stream:stream.write(text)
# Seal actual locations separately, then bind that manifest in the packet seal.
with (P/'EXECUTION-SHA256SUMS').open('x') as stream:
    for path,meta in execution_files.items():stream.write(meta['sha256']+'  '+path+'\n')
files=[p for p in sorted(P.rglob('*')) if p.is_file() and p.name not in ['SHA256SUMS','SEAL.json'] and p.name not in ['assembly.stdout','assembly.stderr']]
with (P/'SHA256SUMS').open('x') as stream:
    for path in files:stream.write(sha(path)+'  '+str(path.relative_to(P))+'\n')
save('SEAL.json',{'schema':'next-qcd-native-execution-evidence-seal/1','manifest_sha256':sha(P/'SHA256SUMS'),'bound_files':len(files),'execution_manifest_sha256':sha(P/'EXECUTION-SHA256SUMS'),'execution_files':len(execution_files),'report_sha256':sha(P/'report.json'),'final_scientific_approval':False,'excluded_active_assembly_streams':['assembly.stdout','assembly.stderr']})
print(json.dumps({'status':report['status'],'counts':counts,'inputs_unchanged':unchanged,'seal_sha256':sha(P/'SEAL.json')}),flush=True)
