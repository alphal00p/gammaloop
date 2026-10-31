"""Pure saved-record probes. No module/candidate/library loading or native execution."""
from pathlib import Path
import ast
import copy
import hashlib
import json
import re

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
P = ROOT / 'target/raw-next-qcd-preflight-20261031-r1'
R1 = ROOT / 'target/raw-next-qcd-captures-20261031-r1'
R2 = ROOT / 'target/raw-next-qcd-captures-20261031-r2'
cached = json.loads((OUT/'frozen-before-files.json').read_text())
def digest(path):
    key = str(Path(path).absolute())
    if key not in cached:
        with Path(path).open('rb') as stream:
            cached[key] = hashlib.file_digest(stream, 'sha256').hexdigest()
    return cached[key]
def need(ok, detail):
    if not ok:
        raise ValueError(detail)
def load(path):
    return json.loads(Path(path).read_text())
def functions(path, names, namespace):
    nodes = [n for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == set(names)
    assert not any(isinstance(n, (ast.Import, ast.ImportFrom)) for node in nodes for n in ast.walk(node))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), namespace)

# Independent complete token grammar for Rust Debug summaries and literal Rows.
# This does not share the repaired delimiter splitting or the old regex parser.
class RecordParser:
    def __init__(self, text):
        self.tokens = re.findall(r'0x[0-9a-fA-F]+|[A-Za-z_][A-Za-z_0-9]*|-?[0-9]+|[{}\[\](),:&]', text)
        assert ''.join(self.tokens) == re.sub(r'\s+', '', text), 'unsupported record syntax'
        self.pos = 0
    def take(self, expected=None):
        t = self.tokens[self.pos]
        self.pos += 1
        assert expected is None or expected == t, (expected,t)
        return t
    def value(self):
        t = self.take()
        if t == '&':
            return self.value()
        if re.fullmatch(r'-?[0-9]+|0x[0-9a-fA-F]+', t):
            return int(t, 16 if t.startswith('0x') else 10)
        if t in ['[','(']:
            end = ']' if t == '[' else ')'
            result = []
            while self.tokens[self.pos] != end:
                result.append(self.value())
                if self.tokens[self.pos] != ',':
                    break
                self.take(',')
            self.take(end)
            return result
        assert t in ['Summary','Row','Fraction']
        self.take('{')
        result = {}
        while self.tokens[self.pos] != '}':
            name = self.take()
            assert name not in result, 'duplicate record field'
            self.take(':')
            result[name] = self.value()
            if self.tokens[self.pos] != ',':
                break
            self.take(',')
        self.take('}')
        if t == 'Fraction':
            assert set(result)=={'numerator','denominator'} and result['denominator'] > 0
            return [result['numerator'],result['denominator']]
        return result
    def parse(self):
        result = self.value()
        assert self.pos == len(self.tokens)
        return result

owner = P/'capsule/helpers/gauge/4-row.py'
semantic = {'re': re}
functions(owner, ['parsed_summary','expected_row'], semantic)
fields = ('diagrams','topologies','assignments','inverse_symmetry','weights','symmetry_factors','content','numerators','named','ids')
namespace = {'Path':Path,'ast':ast,'json':json,'re':re,'need':need,'P':P,'FIELDS':fields,'_scope':semantic,
             'MAIN':P/'capsule/source/original_reference/crates/feynkit-generator/tests/generation_oracle/main.rs'}
functions(R2/'row_boundary.py', ['fields','summary_fields','expected_fields','assertion_boundary'], namespace)
rows = load(P/'capture-command-plans.json')['commands']
versions = load(P/'expected-rows.json')
combined = load(R2/'combined-twelve-commands.json')
records = []
last_end = ''
for index,row in enumerate(rows):
    root = Path(row['records_root'])
    command,result = load(root/'command.json'),load(root/'result.json')
    out,err=(root/'stdout').read_text(),(root/'stderr').read_text()
    assert command['plan']==row==combined[index]['plan']
    assert command['deadline_seconds']==240 and command['retries']==0
    assert command['controller_sha256']==digest((R1 if index<3 else R2)/'controller.py')
    assert last_end < result['start'] and command['start']==result['start']
    last_end=result['end']
    assert result['id']==row['id'] and result['process_reaped'] and not result['timed_out']
    assert result['native_launched'] and result['licensed_phase_started'] and not result['cleanup']['remaining_pids']
    assert digest(root/'stdout')==result['stdout_sha256'] and digest(root/'stderr')==result['stderr_sha256']
    assert command['inputs_before']=={p:digest(p) for p in command['inputs_before']}
    if index==2:
        assert 'inputs_after' not in result and not result['validated_assertion_boundary']
        assert result['error']=={'phase':'assertion_boundary','type':'ValueError','predicate':'ten Summary fields'}
    else:
        assert result['inputs_after']==command['inputs_before'] and result['validated_assertion_boundary']
        assert 'error' not in result
    boundary = namespace['assertion_boundary'](row,result['exit'],out,err,Path(row['writable_root']))
    literal_text=versions[row['case']]['original_reference']['text']
    literal=RecordParser(literal_text[literal_text.index('Row {'):].rstrip(';')).parse()
    assert set(literal)==set(fields) and literal==namespace['expected_fields'](row['case'])
    actual=None
    if row['producer'].startswith('current'):
        texts=[next(line.split(': ',1)[1] for line in err.splitlines() if line.startswith(prefix)) for prefix in ['  left: ',' right: ']]
        actual,expected=map(lambda t:RecordParser(t).parse(),texts)
        assert expected==literal
        assert [f for f in fields if actual[f]!=expected[f]]==['numerators']
        assert namespace['summary_fields'](texts[0])==actual and namespace['summary_fields'](texts[1])==expected
    records.append({'id':row['id'],'exit':result['exit'],'native_pass':boundary['native_test_pass'],'actual':actual,'original_Row':literal,'original_controller_error':result.get('error')})

control=[]
def reject(name, call, contains):
    try:
        call()
    except ValueError as error:
        assert contains in str(error), (name,str(error),contains)
        control.append({'name':name,'rejected_by':str(error),'setup_failure_credit':False})
    else:
        raise AssertionError('negative accepted: '+name)

row=rows[2];root=Path(row['records_root']);work=Path(row['writable_root'])
out,err=(root/'stdout').read_text(),(root/'stderr').read_text()
left=next(line.removeprefix('  left: ') for line in err.splitlines() if line.startswith('  left: '))
old={'need':need,'re':re}
functions(R1/'controller.py',['summary_fields'],old)
reject('original-parser-real-Fraction-regression',lambda:old['summary_fields'](left),'ten Summary fields')
assert len(re.findall(r'(\w+): (.*?)(?=, \w+: | \}$)',left))==12
parser_cases=[
 ('duplicate-top-field',left.replace('diagrams: 5','diagrams: 5, diagrams: 5',1),'duplicate'),
 ('missing-top-field',left.replace('diagrams: 5, ','',1),'exact ten'),
 ('extra-top-field',left[:-2]+', unknown: 0 }','exact ten'),
 ('duplicate-Fraction-field',left.replace('numerator: 5','numerator: 5, numerator: 5',1),'duplicate'),
 ('missing-Fraction-field',left.replace('numerator: 5, ','',1),'Fraction fields'),
 ('extra-Fraction-field',left.replace('numerator: 5','unknown: 0, numerator: 5',1),'Fraction fields'),
 ('zero-denominator',left.replace('denominator: 1','denominator: 0',1),'positive'),
 ('negative-denominator',left.replace('denominator: 1','denominator: -1',1),'integer syntax'),
 ('scalar-type-boolean',left.replace('diagrams: 5','diagrams: True',1),'unsigned'),
 ('scalar-trailing-token',left.replace('diagrams: 5','diagrams: 5x',1),'unsigned'),
 ('mismatched-nesting',left.replace('[(1, 5)]','[(1, 5)}'),'balanced'),
 ('unclosed-nesting',left.replace('denominator: 1 }','denominator: 1',1),'unclosed'),
 ('truncated-record',left[:-2],'framing'),
 ('trailing-record',left+' x','framing'),
 ('nested-wrong-scalar',left.replace('diagrams: 5','diagrams: (5)',1),'unsigned'),
 ('bool-in-symmetry-tuple',left.replace('[(1, 5)]','[(True, 5)]'),'exact tuple'),
]
for name,text,predicate in parser_cases:
    reject(name,lambda text=text:namespace['summary_fields'](text),predicate)
stream_cases=[
 ('wrong-stdout-ID',out.replace('test xs_qcd_gguu_sym_L1 ...','test wrong ...'),err,101,'named test'),
 ('wrong-panic-ID',out,err.replace("thread 'xs_qcd_gguu_sym_L1'","thread 'wrong'"),101,'named assertion'),
 ('wrong-failure-list',out.replace('    xs_qcd_gguu_sym_L1','    wrong'),err,101,'named failed'),
 ('wrong-nonnumerator',out,err.replace('left: Summary { diagrams: 5','left: Summary { diagrams: 6'),101,'sole numerator'),
 ('wrong-original-Row',out,err.replace('right: Summary { diagrams: 5','right: Summary { diagrams: 6'),101,'immutable expected'),
 ('zero-tests',out.replace('running 1 test','running 0 tests'),err,101,'running'),
 ('extra-failure',out.replace('1 failed','2 failed'),err,101,'counts'),
 ('wrong-filtered',out.replace('182 filtered','181 filtered'),err,101,'counts'),
 ('duplicate-result',out+out,err,101,'running'),
 ('wrong-exit',out,err,0,'native exit'),
 ('missing-right-Summary',out,'\n'.join(x for x in err.splitlines() if not x.startswith(' right: ')),101,'Summary pair'),
 ('duplicate-left-Summary',out,err+'\n  left: '+left+'\n',101,'Summary pair'),
 ('second-panic',out,err+"\nthread 'xs_qcd_gguu_sym_L1' panicked at fake:1:\n",101,'named assertion'),
]
for name,stdout,stderr,code,predicate in stream_cases:
    reject(name,lambda stdout=stdout,stderr=stderr,code=code:namespace['assertion_boundary'](row,code,stdout,stderr,work),predicate)

v={'Path':Path,'json':json,'need':need,'digest':digest,'O':R2,'R1':R1,'P':P,
   'IDENT':row['id'],'H':digest(R1/'handoff.json'),'SEAL':digest(R1/'SEAL.json'),
   'assertion_boundary':namespace['assertion_boundary']}
functions(R2/'revalidation.py',['r1_manifest','bound_original','observation','verify_receipt','effective_result'],v)
receipt=load(R2/'invocation3-revalidation.json')
assert v['verify_receipt'](receipt)==receipt
assert v['effective_result'](row)==receipt
assert receipt['inputs_after_pure_revalidation']==load(root/'command.json')['inputs_before']
mutations=[
 ('receipt-ID','id','xs_qcd_gguu_L1--current-pristine','exact invocation3'),
 ('receipt-path','source_result_path',str(root/'other.json'),'path/hash'),
 ('receipt-hash','source_result_sha256','0'*64,'path/hash'),
 ('receipt-r1-hash','r1_handoff_sha256','0'*64,'exact r1 seal'),
 ('receipt-false-pass','native_test_pass',True,'named failure'),
 ('receipt-rerun','native_rerun',True,'named failure'),
 ('receipt-amended','original_result_amended',True,'named failure'),
 ('receipt-parser-hash','parser_sha256','0'*64,'parser source'),
 ('receipt-parser-path','parser_source',str(R1/'controller.py'),'parser source'),
 ('receipt-exit','native_exit',0,'named failure'),
]
for name,key,value,predicate in mutations:
    changed=copy.deepcopy(receipt);changed[key]=value
    reject(name,lambda changed=changed:v['verify_receipt'](changed),predicate)
for name,key,predicate in [('receipt-artifact','bindings','artifact hashes'),('receipt-before','inputs_before','before/after'),('receipt-after','inputs_after_pure_revalidation','before/after')]:
    changed=copy.deepcopy(receipt);changed[key][next(iter(changed[key]))]='0'*64
    reject(name,lambda changed=changed:v['verify_receipt'](changed),predicate)
changed=copy.deepcopy(receipt);changed['boundary']['row_fields']['actual']['diagrams']=6
reject('receipt-nonnumerator-boundary',lambda:v['verify_receipt'](changed),'named failure')
# Every other prerequisite really uses its own result, including the first two.
for other in rows:
    effective=v['effective_result'](other)
    assert effective['id']==other['id'] and effective['validated_assertion_boundary'] and effective['process_reaped']
for name in ['write','subreaper','descendants','reap','bounded']:
    def fn(path):
        return next(ast.dump(n,include_attributes=False) for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    assert fn(R1/'controller.py')==fn(R2/'controller.py')

report={'status':'PASS','record_replays':records,'independent_current_Summary_pairs':6,
        'independent_original_Rows':3,'negative_controls':control,'negative_control_count':len(control),
        'historical_native_passes_retained':6,'current_native_failures_retained':6,
        'old_regex_entries':12,'repaired_top_level_fields':10,'receipt_original_after_record_absent':True,
        'candidate_reported_controls':{'positive':4,'negative':31},
        'reviewer_native_invocations':0,'candidate_module_imports':0,
        'method':'Only named pure function ASTs extracted; no controller capture, top-level imports or helper execution. Independent record grammar and original record joins separately check results.'}
(OUT/'parser-probe.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['record_replays','negative_controls']}))
