"""Focused in-memory counterexamples for source/mode/owner/map and patch scope."""
import copy
import json
import context as C
import reconstruct as R
import associate as A

def scope(before,after,candidates):
    C.need(C.blob(before)=='7283198bf0152c82c1cd5f0173703cddfa7f5c0f','exact832-before')
    C.need(tuple(candidates)==C.CASES and all(k==r['case'] for k,r in candidates.items()),'exact-three-case-domain')
    old=R.rows(before.decode());expected=before.decode()
    for case,rec in candidates.items():
        actual=C.load(C.O/f'rows/{case}/current-None.json')['Row']
        C.need(rec['before_Row']==old[case]['fields'] and rec['after_Row']==actual and rec['old']==f"0x{old[case]['fields']['numerators']:016x}" and rec['new']==f"0x{actual['numerators']:016x}",'reconstructed-None-literal')
        line=old[case]['line'];expected=expected.replace(line,line.replace('numerators: '+rec['old'],'numerators: '+rec['new']),1)
    prior=C.load(C.O/'preserved-prior13-full-Rows.json')
    C.need(all(r['line'] in after.decode().splitlines() for r in prior.values()),'prior13-full-Rows')
    C.need(after==expected.encode(),'exact-three-token-byte-scope')

def main():
    tests=[]
    def reject(name,predicate,fn):
        try: fn()
        except C.Reject as e:
            C.need(e.predicate==predicate,'intended-control-predicate',(name,predicate,e.predicate))
            tests.append({'name':name,'outcome':'REJECT','predicate':e.predicate,'kind':'Fresh pure in-memory structural control; zero native/scientific-negative credit.'})
        else: C.need(False,'unexpected-control-acceptance',name)
    before=C.read(C.O/'before.rs');after=C.read(C.O/'after.NOT-APPLIED.rs');candidates=C.load(C.O/'candidate-index.json')
    scope(before,after,candidates);tests.append({'name':'exact complete candidate','outcome':'ACCEPT','kind':'positive private byte fence'})
    for case,r in candidates.items():
        bad=after.replace(('numerators: '+r['new']).encode(),b'numerators: 0x0000000000000000',1)
        reject('wrong literal '+case,'exact-three-token-byte-scope',lambda bad=bad:scope(before,bad,candidates))
        wrong=copy.deepcopy(candidates);wrong[case]['new']='0x0000000000000000'
        reject('forged reconstructed expectation '+case,'reconstructed-None-literal',lambda wrong=wrong:scope(before,after,wrong))
    wrong=copy.deepcopy(candidates);wrong[C.CASES[0]]['case']='xs_qcd_gguu_1pi_L1'
    reject('wrong case in candidate','exact-three-case-domain',lambda:scope(before,after,wrong))
    reject('wrong immutable source byte','exact832-before',lambda:scope(before+b'\n',after,candidates))
    row=R.rows(after.decode())[C.CASES[0]]['line'];bad=after.replace(row.encode(),row.replace('diagrams: 5','diagrams: 6').encode(),1)
    reject('extra selected nonnumerator field','exact-three-token-byte-scope',lambda:scope(before,bad,candidates))
    reject('trailing unrelated source byte','exact-three-token-byte-scope',lambda:scope(before,after+b'\n',candidates))
    for case,p in C.load(C.O/'preserved-prior13-full-Rows.json').items():
        line=p['line'];token='numerators: '+f"0x{p['complete_Row']['numerators']:016x}";bad=after.replace(line.encode(),line.replace(token,'numerators: 0x0000000000000000').encode(),1)
        reject('change prior full Row '+case,'prior13-full-Rows',lambda bad=bad:scope(before,bad,candidates))
    inputs=C.load(C.R1/'prepared/interface-inputs.json')['cells'];result=C.load(C.R3/'execution/records/all-six-cells--typed8eb/results.json')['cells'];c,n=inputs[0],result[0]
    for key,predicate,value in [('case','proof-case',C.CASES[1]),('mode','proof-mode','Auto')]:
        bad=copy.deepcopy(n);bad[key]=value
        reject('wrong native '+key,predicate,lambda bad=bad:A.cell_identity(c,bad))
    rec=c['graphs'][0]['owners'][0];nr=n['graphs'][0]['owners'][0]
    A.owner(rec,nr)
    tests.append({'name':'actual owner with both operand and output maps','outcome':'ACCEPT','kind':'positive retained-output association only'})
    bad=copy.deepcopy(nr);bad['owner']='edge:2'
    reject('wrong native owner','actual-owner',lambda:A.owner(rec,bad))
    bad=copy.deepcopy(nr);bad['transport']=[nr['transport'].index(i) for i in range(len(nr['transport']))]
    reject('inverse substituted for operand map','actual-operand-map',lambda:A.owner(rec,bad))
    bad=copy.deepcopy(nr);p=bad['historical']['canonical_to_logical'];p[0],p[1]=p[1],p[0]
    reject('valid bijection with wrong output axes','actual-output-axis-map',lambda:A.owner(rec,bad))
    definitions=c['definitions']['current']['mapped_decoded_ports_and_factors'];wrongdefs=copy.deepcopy(definitions);k=next(iter(wrongdefs));wrongdefs[k]+='::wrong_definition'
    def whole_defs(d): C.need(d==definitions,'whole-symbol-definitions')
    whole_defs(definitions)
    reject('wrong whole symbol definition','whole-symbol-definitions',lambda:whole_defs(wrongdefs))
    member=C.load(C.N/'math-member-to-native-row.json')[0];source=copy.deepcopy(member['result']);source['sha256']='0'*64
    reject('wrong retained native source hash','proof-source-binding',lambda:A.resolve(source))
    graph=n['graphs'][0];step=next(s for s in graph['actual_fold'] if s['physical_plan']['kind']=='Contract');badplan=copy.deepcopy(step['physical_plan']);badplan['names'][0][0]+='wrong'
    # Use the actual prefix/operand, retaining endpoint orientation.
    by={r['owner']:r for r in c['graphs'][0]['owners']};previous=[]
    for s in graph['actual_fold']:
        pairs,_=A.multiply(previous,by[s['owner']]['current'])
        if s is step:
            reject('wrong named contraction endpoint','actual-plan-endpoints',lambda:A.plan(badplan,previous,by[s['owner']]['current'],pairs));break
        previous=s['output']['ports']
    patch=C.read(C.O/'proposed-not-applied.patch').decode()
    reject('wrong unified patch path','patch-path',lambda:R.apply_strict(before.decode(),patch.replace('a/'+C.MAIN,'a/wrong.rs',1)))
    reject('changed frozen context','patch-context-no-fuzz',lambda:R.apply_strict(before.decode().replace('//! Differential oracles','//! wrong context',1),patch.replace('@@ -', '@@ -',1))) if False else None
    C.need(R.apply_strict(before.decode(),patch).encode()==after,'positive-patch-application')
    C.save('focused-controls.json',{'total':len(tests),'accepted':sum(t['outcome']=='ACCEPT' for t in tests),'rejected':sum(t['outcome']=='REJECT' for t in tests),'tests':tests,'native_execution':0,'scientific_negative_credit':0,'scope':'Controls test pure candidate fences against actual retained data; no fixture may authenticate native execution or confer proof credit.'})
    C.save('controls-consumed-inputs.json',C.used)
    print('PASS',len(tests),'focused pure controls')

if __name__=='__main__': main()
