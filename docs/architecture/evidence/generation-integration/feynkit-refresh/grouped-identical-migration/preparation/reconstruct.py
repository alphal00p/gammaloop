"""Source-equivalent finite Census reconstruction. No producer or native imports."""
import copy
import difflib
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
P = W.parent / 'generation-feynkit-refresh/.local/feynkit-refresh/grouped-none-qualification-20261030T191919'
CASE = 'xs_qed_eemm_identical_L2'

def need(ok, predicate):
    if not ok:
        raise AssertionError(predicate)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def encode(x):
    return json.dumps(x, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False).encode()

def fnv(b, h=0xcbf29ce484222325):
    for v in b:
        h = ((h ^ v) * 0x100000001b3) & ((1 << 64) - 1)
    return h

def number(x):
    return x.to_bytes(8, 'little')

def text(s):
    return s.encode() + b'\xff'

def multiset(xs):
    return fnv(b''.join(number(x) for x in sorted(xs)))

class Scalar:
    """Only observed bookkeeping heads; never normalize a graph numerator."""
    def __init__(self, s):
        self.s, self.i = s, 0

    def expr(self):
        x = self.product()
        while self.i < len(self.s) and self.s[self.i] == '+':
            self.i += 1
            x += self.product()
        return x

    def product(self):
        x = self.power()
        while self.i < len(self.s) and self.s[self.i] == '*':
            self.i += 1
            x *= self.power()
        return x

    def power(self):
        x = self.atom()
        if self.i < len(self.s) and self.s[self.i] == '^':
            self.i += 1
            n = self.atom()
            need(n.denominator == 1, 'scalar-integer-power')
            x **= n.numerator
        return x

    def atom(self):
        if self.s[self.i] == '(':
            self.i += 1
            x = self.expr()
            need(self.s[self.i] == ')', 'scalar-closing-parenthesis')
            self.i += 1
            return x
        m = re.match(r'-?[0-9]+', self.s[self.i:])
        if m:
            self.i += len(m[0])
            return Fraction(int(m[0]))
        m = re.match(r'feynkit_generator_factor::\{\}::([A-Za-z]+)\(', self.s[self.i:])
        need(m is not None, 'scalar-known-qualified-head')
        self.i += len(m[0])
        args = [self.expr()]
        while self.s[self.i] == ',':
            self.i += 1
            args.append(self.expr())
        need(self.s[self.i] == ')', 'scalar-head-close')
        self.i += 1
        if m[1] == 'NumeratorDependentGrouping':
            need(len(args) == 3, 'scalar-group-arity')
            return args[1] * args[2]
        need(m[1] in ['AutG', 'AntiFermionSpinSumSign', 'ExternalFermionOrderingSign', 'InternalFermionLoopSign'] and len(args) == 1, 'scalar-source-sign-head')
        return args[0]

def weight(s):
    p = Scalar(s)
    x = p.expr()
    need(p.i == len(s), 'scalar-complete-parse')
    return x

def payload_digests(payload):
    p = copy.deepcopy(payload)
    chunks = [encode(p.pop(k)) for k in ['numerator', 'numerator_prefactor', 'projector']]
    for k in ['num_interface', 'numerator_prefactor_interface', 'projector_interface']:
        if k in p:
            chunks.append(encode(p.pop(k)))
    for owner in p['vertices'] + [e[1] for e in p['edges']]:
        chunks.append(encode(owner['numerator']))
        owner['numerator'] = None
        if 'num_interface' in owner:
            chunks.append(encode(owner.pop('num_interface')))
    for k in ['model', 'symmetry_factor', 'overall_factor']:
        p.pop(k)
    return fnv(encode(p)), fnv(b''.join(chunks))

def census(data, report):
    content, nums, ids, per = [], [], [], []
    named = b''
    inv, weights = Fraction(), Fraction()
    histogram = {}
    for d in data['diagrams']:
        key = fnv(encode(d['key']))
        frame, numerator = payload_digests(d['payload'])
        aut = int(d['aut'])
        need(aut == d['payload']['symmetry_factor'], 'aut-payload-binding')
        need(d['overall_factor'] == d['payload']['overall_factor'], 'factor-payload-binding')
        inv += Fraction(1, aut)
        weights += weight(d['overall_factor'])
        histogram[aut] = histogram.get(aut, 0) + 1
        content.append(fnv(number(key) + number(aut) + text(d['overall_factor']) + number(frame)))
        nums.append(fnv(number(key) + number(numerator)))
        named += (text(d['name']) if d['name'] is not None else b'\xfe') + number(key)
        ids.append(fnv(int(d['id'], 16).to_bytes(16, 'little')))
        per.append({'id': d['id'], 'key_digest': f'{key:016x}', 'payload_sha256': sha(encode(d['payload'])), 'numerator_digest': f'{numerator:016x}', 'key_bound_numerator_digest': f'{nums[-1]:016x}'})
    for g in data['groups']:
        named += number(g['master'])
        for m in g['members']:
            named += (text(m['name']) if m['name'] is not None else b'\xfe') + text(m['ratio']) + text(m['overall_factor'])
    pair = lambda f: [f.numerator, f.denominator]
    row = {'diagrams': len(per), 'topologies': int(re.search(r'topology_count: ([0-9]+)', report)[1]), 'assignments': int(re.search(r'interaction_assignment_count: ([0-9]+)', report)[1]), 'inverse_symmetry': pair(inv), 'weights': pair(weights), 'symmetry_factors': [list(x) for x in sorted(histogram.items())], 'content': f'0x{multiset(content):016x}', 'numerators': f'0x{multiset(nums):016x}', 'named': f'0x{fnv(named):016x}', 'ids': f'0x{multiset(ids):016x}'}
    return row, per

def main():
    pristine = P / f'pristine-runtime/dumps/{CASE}.classes-None-t1.json'
    data = j(pristine)
    report = read(P / f'pristine-runtime/dumps/{CASE}.census').decode().splitlines()[0]
    row, per = census(data, report)
    expected = {'diagrams':29,'topologies':157,'assignments':714,'inverse_symmetry':[29,1],'weights':[4,1],'symmetry_factors':[[1,29]],'content':'0xc6e7d902ad8d316f','numerators':'0x00fac78da0e0d159','named':'0x3d72c09c50112e41','ids':'0xcd1abf37940d48cf'}
    need([k for k in row if row[k] != expected[k]] == ['numerators'], 'exact-nine-other-fields')
    save('reconstructed-None-row.json', {'derivation':'Exact pristine physical payload, canonical full key, source Digest/Entry/Summary algorithm. Recorded/native hash only corroborates after derivation.', 'source':str(pristine), 'actual':row,'expected':expected,'per_diagram':per,'nine_other_fields_match':True,'report':report})
    before = read(O / 'main.before.rs')
    lines = before.decode().splitlines(keepends=True)
    rows = [x for x in lines if re.match(r'\s+(fast|heavy) \w+: Row \{', x)]
    need(len(rows) == 149, 'current-149-Rows')
    oldline, = [x for x in rows if re.match(r'\s+(fast|heavy) '+CASE+r':', x)]
    old, new = expected['numerators'], row['numerators']
    need(oldline.count('numerators: '+old) == 1, 'one-old-row-token')
    newline = oldline.replace('numerators: '+old, 'numerators: '+new)
    after = before.replace(oldline.encode(), newline.encode(), 1)
    need(after.replace(newline.encode(),oldline.encode(),1) == before, 'all-other-bytes-preserved')
    need(len(after) == len(before), 'equal-length-single-token')
    (O/'main.after.NOT-APPLIED.rs').write_bytes(after)
    path = 'crates/feynkit-generator/tests/generation_oracle/main.rs'
    blob = lambda b: hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
    diff = ''.join(difflib.unified_diff(before.decode().splitlines(keepends=True), after.decode().splitlines(keepends=True), fromfile='a/'+path,tofile='b/'+path))
    patch = 'diff --git a/'+path+' b/'+path+'\nindex '+blob(before)+'..'+blob(after)+' 100644\n'+diff
    (O/'proposed-NOT-APPLIED.patch').write_text(patch)
    # Validate patch reconstruction in memory, without git apply or index access.
    body = patch.splitlines(keepends=True)
    removed = [x[1:] for x in body if x.startswith('-') and not x.startswith('---')]
    added = [x[1:] for x in body if x.startswith('+') and not x.startswith('+++')]
    need(removed == [oldline] and added == [newline], 'one-removed-one-added-line')
    need(before.replace(removed[0].encode(), added[0].encode(), 1) == after, 'patch-reconstruction')
    save('patch-binding.json', {'status':'NOT_APPLIED','approval':False,'immutable_base':'d0fb4e906140a6a0fc8dd6401ca839cab0c2734d','producer_applicability_code':'ffa76ec420642c03349ab756d2d7720ab93c0057','path':path,'old':old,'new':new,'before_blob':blob(before),'after_blob':blob(after),'before_sha256':sha(before),'after_sha256':sha(after),'patch_sha256':sha(patch.encode()),'Rows':149,'changed_tokens':1,'other_148_Rows_unchanged':True,'all_other_bytes_unchanged':True,'all_Rows_before_sha256':{re.search(r'(\w+): Row',x)[1]:sha(x.encode()) for x in rows},'preserves_all_prior_12_migrations':True,'gates':['AUTO_AGGREGATION_REVIEW_PENDING','NEW_MIGRATION_REVIEW_REQUIRED'],'runtime_tested_after_patch':False})
    print(json.dumps({'status':'RECONSTRUCTED_NOT_APPLIED','old':old,'new':new,'new_decimal':int(new,16),'nine_other_fields':True,'Rows':149,'before_blob':blob(before),'after_blob':blob(after)}))

if __name__ == '__main__':
    try:
        main()
    except BaseException:
        with (O/'reconstruction-failed-attempts.log').open('a') as f:
            f.write(traceback.format_exc())
        raise
