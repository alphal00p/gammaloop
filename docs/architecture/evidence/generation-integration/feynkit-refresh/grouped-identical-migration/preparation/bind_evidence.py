"""Prepare finite proof bindings and holdouts. No historical helper execution."""
import collections
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import runpy
import traceback

U = runpy.run_path(str(Path(__file__).with_name('inspect.py')))
W, O, N, A = (U[k] for k in ('W', 'O', 'N', 'A'))
read, j, save = (U[k] for k in ('read', 'j', 'save'))
P = W.parent/'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
T = W/'target/xs-qed-eemm-identical-fresh8eb-r3'
RAW = W/'target/xs-qed-l2/runtime/current/captures/xs_qed_eemm_L2/classes-None-t1'
CASE = 'xs_qed_eemm_identical_L2'

def need(ok, p):
    if not ok:
        raise AssertionError(p)

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def canon(x):
    return json.dumps(x, sort_keys=True, ensure_ascii=False, separators=(',', ':'))

def digest(x):
    return hashlib.sha256(canon(x).encode()).hexdigest()

def inventory():
    pins = [
        (N/'receipt.json','4a91e2bd0036c7e32933069dfa147d3141ab4438a706bb7d96dc81a56b8994ac'),
        (N/'seal.json','c2e4448d0930b5422e8a40bcda6d0a9ab81cc41a1ef9e1be23c3fbdf177f7710'),
        (W/'.local/review-grouped-none-aggregation-astra-20261030-independent-b1/review.json','45368ac2dc5783a6ba31daa9b7e8bd7f03a80746bf3b07b6d3e34cc909adcb81'),
        (A/'receipt.json','d7ad071e60426eafe5662b9f5b57fb8484508e99e9967252a342e1abc3f3fce3'),
        (A/'seal.json','a2446a215715c2e4757dd08b6f4cccdcc953b2543932671e3c91cd1e8283fb3c'),
        (A/'reviews/Scope-A.json','80e87bc1dfc062183d1d5bb62b2ca889b1ebc6f82461fc9e1b033259c79e3d5a'),
        (A/'reviews/Scope-B.json','ac52b96f58e87e9180c5c7083d9dccdf973de38fa66aee549bcb7341d01dd76a'),
    ]
    result = {}
    for p,h in pins:
        need(sha(p)==h, 'supplied-root-pin '+str(p))
        result[str(p)] = h
    for root, seal, base in [(N,N/'seal.json',W),(A,A/'seal.json',A),(P,P/'seal.json',P)]:
        spec = j(seal)
        for path, rec in spec['files'].items():
            p = base/path
            need(sha(p)==rec['sha256'] and p.stat().st_size==rec['bytes'], 'sealed-file '+str(p))
            result[str(p)] = rec['sha256']
    # Direct load-bearing original records: pin before any finite projection.
    for p in [T/'maps.json',T/'raw-bridge.json',T/'fresh-triangle.json',T/'spinor-adoption.json',T/'final-seal.json',N/'constituents.json',N/'native-reuse.json',N/'producer-binding.json',N/'qualifications.json',A/'premise-matrix.json',A/'guard-receipt.json',A/'qualifications.json']:
        result[str(p)] = sha(p)
    save('protected-inputs-before.json',result)
    return result

def one(rows, stage):
    vals = [x['value'] for x in rows if x['stage']==stage]
    need(len(vals)==1, 'single-stage '+stage)
    return vals[0]

def owner_payload(p, owner):
    if owner.startswith('vertex:'):
        return p['vertices'][int(owner.split(':')[1])], 'numerator', 'num_interface'
    if owner.startswith('edge:'):
        return p['edges'][int(owner.split(':')[1])][1], 'numerator', 'num_interface'
    return p, owner, {'numerator':'num_interface','projector':'projector_interface','numerator_prefactor':'numerator_prefactor_interface'}[owner]

def owners(g):
    out = {x['owner']:x['tensor'] for x in g['owners']}
    need(len(out)==len(g['owners'])==18, 'unique-18-named-owners')
    for owner,t in out.items():
        p, field, interface = owner_payload(g['json'],owner)
        need(p[field]==t['decoded']==t['raw'], 'owner-raw-decoded-source '+owner)
        if interface in p:
            need(p[interface]==t['ports'], 'declared-interface-ordered '+owner)
        need(len(t['ports'])==len(set(t['ports'])), 'repeat-free-named-ports '+owner)
        for port in t['ports']:
            explicit = r'gammalooprs::\{spenso::index\}::(?:hedge|edge)\([0-9]+,1\)'
            spacetime = r'spenso::\{spenso::representation,spenso::self_dual\}::(?:bis|mink)\(4,' + explicit + r'\)'
            color = r'spenso::\{spenso::representation,spenso::dualizable\}::cof\(3,' + explicit + r'\)'
            need(any(re.fullmatch(pat,port) is not None for pat in [spacetime,color,r'spenso::\{spenso::lower\}::dind\('+color+r'\)']), 'finite-explicit-port-domain')
    return out

def project(g):
    out = copy.deepcopy(g['json'])
    out.pop('name',None)
    out.pop('id',None)  # Rounded json.id is never identity evidence in this new join.
    return out

def bind():
    rows = [json.loads(x) for x in read(P/'observer-runtime/captures/ledger.jsonl').splitlines()]
    live = one(rows,'live-inputs')
    raw = one(rows,'raw-diagrams')
    final = one(rows,'final-diagrams')
    rawc, finalc = one(rows,'raw-census'),one(rows,'final-census')
    pristine = j(P/f'pristine-runtime/dumps/{CASE}.classes-None-t1.json')
    pristine_raw = j(P/f'pristine-runtime/dumps/{CASE}.classes-None-t1.sources.json')
    need(finalc==pristine and rawc==pristine_raw, 'observer-pristine-complete-census')
    need(read(P/f'pristine-runtime/dumps/{CASE}.classes-None-t1.sources.json')==read(RAW/'payload.json'),'raw30-approved-complete-payload-bytes')
    model, request, options = (j(RAW/(k+'.json')) for k in ['model','request','options'])
    need(all(live[k]==v for k,v in [('model',model),('request',request),('options',options)]),'exact-None-build-inputs')
    need(live['request']['classes']=='None','live-None-request-mode')
    meta=j(RAW/'capture.json')
    need(live['Run']==meta['run'],'raw30-Run-exact')
    need(j(RAW/'grouping.json')=='none','raw30-ungrouped-boundary')
    need(live['grouping']=={'identical':{'check_canonical_numerator':False,'differentiate_particle_masses_only':True,'fully_numerical_substitution':False,'number_of_numerical_samples':5,'numerical_sample_seed':3,'symmetric_polarizations':False}},'exact-Identical-options')
    approved = j(RAW/'typed-owners.json')
    source_records=[]
    for idx,(g,c,old) in enumerate(zip(raw,rawc['diagrams'],approved,strict=True)):
        need(g['id']==c['id']==old['id'],'exact-raw-native-ID')
        need(project(g)==c['payload'],'raw-physical-census-binding')
        actual=owners(g)
        oldmap={o['owner']:o for o in old['owners']}
        need(set(actual)==set(oldmap),'raw30-owner-domain')
        for owner,t in actual.items():
            expected={'owner':owner,'expression':t['raw'],'producer_decoded_expression':t['decoded'],'ordered_ports':t['ports'],'is_metric':t['is_metric'],'is_composite':t['is_composite']}
            need(oldmap[owner]==expected,'raw30-complete-owner-binding')
        source_records.append({'raw_index':idx,'native_ID':g['id'],'full_key_sha256':digest(c['key']),'payload_sha256':digest(c['payload']),'source_factor':c['overall_factor'],'owners':{k:digest(v) for k,v in actual.items()}})
    need(len(raw)==len(rawc['diagrams'])==len(approved)==30,'all-raw30')
    pairs = [x['value'] for x in rows if x['stage']=='staged-pair-final']
    maps=j(T/'maps.json')['current']
    need(len(pairs)==30 and len(maps)==60,'sixty-source-phase-maps')
    map_bindings=[]
    for p in pairs:
        idx=p['raw_index']
        need(p['source_index']==idx,'raw-source-index')
        for side in ['comparison','representative']:
            target=p['pair'][side]
            owner_target=owners(target)
            m,=[x for x in maps if x['raw_index']==idx and x['side']==side]
            # Map typed labels, then record the actual operand/output-axis permutation.
            pm={}
            for name,t in owners(raw[idx]).items():
                if name.startswith('vertex:'):
                    mapped='vertex:'+str(m['vertices'][int(name.split(':')[1])])
                elif name.startswith('edge:'):
                    mapped='edge:'+str(m['edges'][name.split(':')[1]][0])
                else:
                    mapped=name
                ports=[re.sub(r'::hedge\(([0-9]+),1\)',lambda q:'::hedge('+str(m['hedges'][q[1]])+',1)',v) for v in t['ports']]
                dest=owner_target[mapped]['ports']
                need(collections.Counter(ports)==collections.Counter(dest),'mapped-owner-port-domain '+str(idx)+' '+side+' '+name)
                permutation=[dest.index(v) for v in ports]
                pm[name]={'target_owner':mapped,'source_ordered_ports':t['ports'],'mapped_ordered_ports':ports,'actual_target_ordered_ports':dest,'source_axis_to_target_axis':permutation,'role':'output-axis map' if name in ['numerator','projector','numerator_prefactor'] else 'operand-axis map'}
            map_bindings.append({'raw_index':idx,'source_native_ID':raw[idx]['id'],'side':side,'target_native_ID':target['id'],'full_map':m,'owners':pm,'expression_covariance':'Reused accepted None finite owner-covariance theorem; not established from ports alone.'})
    need(len(final)==len(pristine['diagrams'])==29,'29-final-payloads')
    final_owners=[]
    for g,c in zip(final,pristine['diagrams'],strict=True):
        need(g['id']==c['id'] and project(g)==c['payload'],'final-native-ID-full-payload')
        om=owners(g)
        final_owners.append({'id':g['id'],'key_sha256':digest(c['key']),'physical_payload_sha256':digest(c['payload']),'owners':om})
    partitions=[]
    for plan,g,c in zip([x['value'] for x in rows if x['stage']=='collapse-plan'],pristine['groups'],pristine['diagrams'],strict=True):
        member_records=[]
        for member in g['members']:
            idx=member['source_index']
            need(member['id']==raw[idx]['id'] and member['name']==raw[idx]['name'],'member-source-identity')
            need(member['overall_factor']==rawc['diagrams'][idx]['overall_factor'],'member-source-factor')
            need(type(member['ratio']) is str and member['ratio']=='1','ratio-one-membership')
            member_records.append({'raw_index':idx,'native_ID':raw[idx]['id'],'source_factor':member['overall_factor'],'ratio':member['ratio'],'map_records':[i for i,m in enumerate(map_bindings) if m['raw_index']==idx]})
        need(plan['members']==[m['name'] for m in g['members']],'collapse-members-match')
        need(plan['factor']==c['overall_factor'],'collapse-final-factor')
        partitions.append({'final_native_ID':c['id'],'full_key_sha256':digest(c['key']),'master_retained_position':g['master'],'members':member_records,'final_factor':c['overall_factor'],'cuts_sha256':digest(c['payload']['cuts']),'thresholds_sha256':digest(c['payload']['topology_threshold_candidates']),'collapse_plan_sha256':digest(plan),'payload_sha256':digest(c['payload'])})
    need(sorted(m['raw_index'] for g in partitions for m in g['members'])==list(range(30)),'no-omission-or-duplicate')
    need([sorted(m['raw_index'] for m in g['members']) for g in partitions if len(g['members'])>1]==[[22,23]],'exact-sole-merge-membership')
    save('raw30-current-owner-bindings.json',{'input_equality':{'model':True,'request':True,'options':True,'Run':True,'raw_grouping':'None','grouped_grouping':'Identical','last_shared_boundary':'raw built source payloads before grouping'},'sources':source_records,'raw30_native_results':'Reused exact approved c95 native inputs/results only; no kernel or program relabel.'})
    save('source-phase-owner-axis-maps.json',map_bindings)
    save('final-physical-owner-bindings.json',final_owners)
    save('membership-weight-payload-chain.json',{'groups':partitions,'no_omission_duplicate':True,'sole_merge':[22,23],'ratio':'1','current_correctness_argument':'Accepted independent raw30 correctness + exact current raw30/input/owner rebind + accepted conditional None signed sewing/chart/owner-covariance/closed-spinor relation + complete partition, source-factor, no-zero/cancellation, cut/threshold and final physical joins. No historical grouped frame is a correctness premise.','limits':'The new migration requires independent review of this composition and its material interface conditions; ports alone and scalar total4 do not prove it.'})
    save('consumer-interface-boundary.json',{'finite_matching_serialized_declarations':sum('num_interface' in obj for d in pristine['diagrams'] for obj in [d['payload']]+d['payload']['vertices']+[e[1] for e in d['payload']['edges']]),'owner_ports_checked':{'raw':540,'source_phase':1080,'final':522},'all_typed_ports':'Explicit repeat-free self-dual bis/mink dimension4 and cof/dual-cof dimension3; complete tags and dual orientation retained, including the quark-loop sources; exact ordered capture ports retained. No compact/unresolved port observed.','source_read':['crates/feynkit-graph/src/numerator_interface.rs:admit/export/canonical_ports','crates/feynkit-graph/src/lib.rs:numerator/fragment_numerator','crates/feynkit-generator/src/grouping.rs:complete.multiply(prefactor).multiply(projector)','crates/idenso/src/tensor/composition.rs:interface_is_named/product_plan/named/free_interface/multiply','crates/idenso/src/tensor/composition/product.rs:multiply_many/merge_named','crates/spenso/src/structure/partial.rs:from_logical_slots/logical_slots'],'inference':'Declared serialization is bound to actual current logical owner order. Current consumers select named plans afresh, with canonical survivors; no arbitrary historical positional plan is reused. This is finite source/data evidence offered to the migration reviewer, not imported three-QED/QCD interface approval.','material_boundary':'P4/accepted raw30 certificates certify named axis multisets and captured interfaces, excluding a positional/free-index theorem. Direct native grouped product-plan/layout/whole-import results are not present for this owner set. Source inspection and ordered declarations narrow this boundary but are not a new independent interface approval.','smallest_remaining_witness_if_source_argument_not_accepted':'For the exact current raw30, 60 comparison/representative graphs and 29 final graphs in this packet, bind logical/canonical layouts and named product operand/output-axis maps (including prefactor/projector and 54 declared final owner interfaces) to the already qualified actual consumer; show finite commutation and reject an order-only mutation. It need not prove arbitrary tensors or historical positional identity. Fresh native work or an additional theorem requires a separate assignment.','status':'NEW_MIGRATION_REVIEW_MUST_DISPOSE_INTERFACE_BOUNDARY','borrowed_three_QED_QCD_approval':False})
    return rows,pristine,live,raw,approved,partitions

def main():
    pins=inventory()
    rows,pristine,live,raw,approved,partitions=bind()
    need(all(sha(Path(p))==h for p,h in pins.items()),'protected-inputs-after')
    save('protected-inputs-after.json',pins)
    save('evidence-binding-result.json',{'status':'FINITE_BINDINGS_PREPARED_FOR_REVIEW','sealed_files_verified':len(pins),'raw_sources':30,'groups':29,'raw_owner_bindings':540,'phase_owner_axis_maps':1080,'final_owner_bindings':522,'historical_grouped_equivalence':False,'Auto_aggregation':'AUTO_AGGREGATION_REVIEW_PENDING','migration':'NEW_MIGRATION_REVIEW_REQUIRED','native_executions':0})
    print(json.dumps({'status':'FINITE_BINDINGS_PREPARED_FOR_REVIEW','sealed_files':len(pins),'raw':30,'groups':29,'phase_owner_maps':1080}))

if __name__=='__main__':
    try:
        main()
    except BaseException:
        with (O/'evidence-failed-attempts.log').open('a') as f:f.write(traceback.format_exc())
        raise
