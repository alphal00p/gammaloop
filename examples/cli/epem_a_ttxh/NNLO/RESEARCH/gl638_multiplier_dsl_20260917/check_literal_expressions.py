"""Evaluate the six literal DOT expressions against their short energy formulas.

Standard-library diagnostic; no GammaLoop execution or source mutation.
"""

import ast
import hashlib
import json
import operator
import random
import re
import tomllib
from decimal import Decimal as D, getcontext
from pathlib import Path

getcontext().prec = 90
ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent / 'Cargo.toml').exists())
DOT = ROOT / 'examples/cli/epem_a_ttxh/NNLO/graphs/GL638.dot'
text = DOT.read_text()
payload = re.search(r'threshold_counterterms = "(.*?)";', text, re.S)[1]
data = tomllib.loads(payload.replace('\\"', '"'))
expressions = []
graph_edges = re.findall(r'(\w+)(?::\d+)?\s*->\s*(\w+)(?::\d+)?\s*\[id=(\d+)', text)
assert len(graph_edges) == 17
for threshold in data['cuts'][0]['thresholds']:
    for variant in threshold['counterterms']:
        if 'multiplier' not in variant:
            continue
        source = variant['multiplier']['expression']
        translated = re.sub(r'Q3\(star,(\d+),cind\((\d+)\)\)', r'q_\1_\2', source)
        translated = translated.replace('P(0,cind(0))', '(Q/2)').replace('P(1,cind(0))', '(Q/2)')
        translated = translated.replace('UFO::MT', 'mt').replace('^', '**').replace('if(', 'choose(')
        expressions.append((threshold['edges'], variant['name'], ast.parse(translated, mode='eval').body))
assert len(expressions) == 6


def evaluate(node, values):
    if isinstance(node, ast.Constant):
        return D(node.value)
    if isinstance(node, ast.Name):
        return values[node.id]
    if isinstance(node, ast.BinOp):
        operation = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
                     ast.Div: operator.truediv, ast.Pow: operator.pow}[type(node.op)]
        return operation(evaluate(node.left, values), evaluate(node.right, values))
    if isinstance(node, ast.Call) and node.func.id == 'sqrt':
        return evaluate(node.args[0], values).sqrt()
    if isinstance(node, ast.Call) and node.func.id == 'choose':
        return evaluate(node.args[1 if evaluate(node.args[0], values) else 2], values)
    raise ValueError(ast.dump(node))


rng = random.Random(638)
rows = []
maximum = D(0)
for index in range(43):
    p, u, s, t = [[D(rng.randint(-350, 350)) for _ in range(3)] for _ in range(4)]
    if index >= 40:
        p = [t[k] + D(10) ** (-10 * (index - 39)) * D(k + 1) for k in range(3)]
    q = {3: p, 6: u, 7: s, 10: t, 12: p, 8: s, 5: t,
         9: [D(0)] * 3, 11: [D(0)] * 3}
    q.update({4: [u[k] + p[k] - t[k] for k in range(3)],
              2: [u[k] - t[k] for k in range(3)],
              13: [p[k] - t[k] for k in range(3)],
              14: [s[k] - p[k] for k in range(3)]})
    energy = {e: (sum(x*x for x in vector) + D(125 if e == 2 else 0 if e in (9,11,13,14) else 173)**2).sqrt()
              for e, vector in q.items()}
    Q = energy[2] + energy[6] + energy[10]
    q.update({0: [D(0), D(0), Q/2], 1: [D(0), D(0), -Q/2],
              15: [D(0), D(0), Q/2], 16: [D(0), D(0), -Q/2]})
    balance = [[D(0)] * 3 for _ in range(10)]
    for tail, head, edge in graph_edges:
        for vertex, sign in [(tail, -1), (head, 1)]:
            if vertex.isdigit():
                for k in range(3):
                    balance[int(vertex)][k] += sign*q[int(edge)][k]
    assert all(x == 0 for vector in balance for x in vector), balance
    values = {f'q_{e}_{k+1}': x for e, vector in q.items() for k, x in enumerate(vector)}
    values.update(Q=Q, mt=D(173))
    H = energy[2] + energy[4] + energy[12] - Q
    P = energy[3] + energy[12] - Q
    Z = energy[3] + energy[10] + energy[13] - Q
    G0 = energy[2] + energy[6] + energy[12] + energy[13] - Q
    G4 = energy[2] + energy[4] + energy[10] + energy[13] - Q
    WH = H*H/(H*H+(P*Z/Q)**2)
    WF = (G0*G4)**2/((G0*G4)**2+(P*Z)**2)
    errors = []
    for edges, name, expression in expressions:
        expected = WF if edges == [8,10,13,14] else WH
        if name == 'shared_1l':
            expected = 1-expected
        actual = evaluate(expression, values)
        errors.append(abs(actual-expected))
    maximum = max(maximum, *errors)
    rows.append({'case': index, 'Q_GeV': str(Q), 'max_absolute_weight_error': str(max(errors)),
                 'soft_q13': index >= 40})
assert maximum < D('1e-70'), maximum
# Unphysical all-zero spatial input with Q=2*mt gives H=P=0: only the
# guarded WH pair has an explicit assigned value there. Do not evaluate WF.
zero_values = dict(values, Q=D(346))
zero_values.update({key: D(0) for key in zero_values if key.startswith('q_')})
common_zero = [str(evaluate(expr, zero_values)) for _, _, expr in expressions[:4]]
assert common_zero == ['1', '0', '1', '0']
output = {'dot': str(DOT.relative_to(ROOT)), 'dot_sha256': hashlib.sha256(DOT.read_bytes()).hexdigest(),
          'literal_expressions': len(expressions), 'decimal_digits': getcontext().prec,
          'cases': rows, 'max_absolute_weight_error': str(maximum),
          'all_vertex_spatial_momentum_balances_exact': True,
          'WH_common_zero_shared_native_A_then_U': common_zero,
          'scope': '43 finite host-shell points including 3 soft-q13 points; exact graph incidence checked; not a runtime/variance test'}
Path(__file__).with_name('results.json').write_text(json.dumps(output, indent=2)+'\n')
print(json.dumps({key: value for key, value in output.items() if key != 'cases'}, indent=2))
