#!/usr/bin/env python3
"""IRH-DEPTH mechanical variant matrix (MATRIX-R45): native inputs.

For every attribute of B_self_energy_L3_006 that enters the product and can be
written in a DOT input, multiply it by eps^-1 written four ways:

  eps  explicit                     (gammalooprs::{}::ε)^(-1)
  dm4  a 1/(D-4) form               (-2)/(gammalooprs::{}::dim-4)
  gtr  a closed metric trace = d    (-2)/((g(mink(4,h81),mink(4,h82)))^2-4)
  dtr  a closed Dirac trace = 4 d   (-8)/(gamma(bis(4,h83),bis(4,h84),mink(4,h85))
                                         *gamma(bis(4,h84),bis(4,h83),mink(4,h85))-16)

on top of the base "hard-v0" (L3_006 with the hard vertex v0 x eps^3; fixtures.py
of the IRH native witnesses, byte for byte), and write one beta input root per
variant: OUT/<variant>/L3/{diagram-000006.dot,generation.json,model.json}
(manifest reduced to that diagram, index 0, as fixtures.py).

Attributes (closed list entry of regularity::FACTORS in brackets):
  v1        quotient vertex numerator          [quotient_vertex_numerator]
  e4        quotient edge 4 (v1 -> v3)         [quotient_edge_numerator]
  v0        hard vertex v0                     [hard_numerator]
  e7        hard edge 7 (v2 -> v4)             [hard_numerator]
  v0k3j1    v0 x eps^4 x S with v1 x eps^-1    [hard_numerator, rows hidden]
  pref      graph numerator_prefactor          [numerator_prefactor, root_projector, root_normalisation]
  proj      graph projector                    [projector]
  overall   graph overall_factor               [overall_factor]  (the manifest's
            overall_factor binding is beta's canonical string of the DOT factor, read
            from overall_canonical.tsv, which the in-tree test
            matrix_r45_the_overall_canonical_table_is_current keeps current)
  power     edge 4 power=<S>                   [propagator_power]
For vertex and edge attributes the graph's aggregate numerator is multiplied by the
same factor (it stays the product of its vertex and edge factors).
The entries without a DOT spelling (insertion_coefficient, root_closure, row_insertion,
compensation, moment_reduction) and the nested cells are in the in-tree test
depth_tests::round5::matrix_r45.

usage: gen_matrix.py PREPARED_L3_DIR OUT_DIR   (prints one TSV line per variant)
"""
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

EPS = 'gammalooprs::{}::ε'
DIM = 'gammalooprs::{}::dim'
MINK = 'spenso::{spenso::representation,spenso::self_dual}::mink'
BIS = 'spenso::{spenso::representation,spenso::self_dual}::bis'
HEDGE = 'gammalooprs::{spenso::index}::hedge'
METRIC = 'spenso::{symmetric,linear,real,spenso::tensor}::g'
GAMMA = 'spenso::{linear,spenso::tensor}::gamma'


def slot(rep, n):
    return f'{rep}(4,{HEDGE}({n},1))'


SPELLINGS = {
    'eps': f'({EPS})^(-1)',
    'dm4': f'(-2)/({DIM}-4)',
    # the review's spelling (make_inputs.py): g_ab g^ab as a squared metric
    'gtr': f'(-2)/(({METRIC}({slot(MINK, 81)},{slot(MINK, 82)}))^2-4)',
    # tr(gamma_a gamma^a) = d tr(1) = 4 d
    'dtr': (f'(-8)/({GAMMA}({slot(BIS, 83)},{slot(BIS, 84)},{slot(MINK, 85)})'
            f'*{GAMMA}({slot(BIS, 84)},{slot(BIS, 83)},{slot(MINK, 85)})-16)'),
}

# attribute -> (closed-list entries, base powers, [(line prefix, field)], also the aggregate)
ATTRIBUTES = {
    'v1': ('quotient_vertex_numerator', {'v0': 3}, [('v1 [', 'numerator')], True),
    'e4': ('quotient_edge_numerator', {'v0': 3}, [('v1 -> v3 [', 'numerator')], True),
    'v0': ('hard_numerator', {'v0': 3}, [('v0 [', 'numerator')], True),
    'e7': ('hard_numerator', {'v0': 3}, [('v2 -> v4 [', 'numerator')], True),
    'v0k3j1': ('hard_numerator', {'v0': 4, 'v1': -1}, [('v0 [', 'numerator')], True),
    'pref': ('numerator_prefactor,root_projector,root_normalisation', {'v0': 3},
             [('graph [', 'numerator_prefactor')], False),
    'proj': ('projector', {'v0': 3}, [('graph [', 'projector')], False),
    'overall': ('overall_factor', {'v0': 3}, [('graph [', 'overall_factor')], False),
    'power': ('propagator_power', {'v0': 3}, [('v1 -> v3 [', 'power')], False),
}


def scale(line, power):
    # fixtures.py
    if not power:
        return line
    return re.sub(r'(?<![a-z_])numerator="([^"]+)"',
                  lambda m: f'numerator="({m[1]})*({EPS})^({power})"', line, count=1)


def base_lines(dot, powers):
    lines = dot.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if line.startswith('  graph ['):
            lines[i] = scale(line, sum(powers.values()))
        for term, power in powers.items():
            if line.startswith('  ' + term + ' ['):
                lines[i] = scale(line, power)
    return lines


def times(lines, prefix, field, factor):
    idx = [i for i, v in enumerate(lines) if v.startswith('  ' + prefix)]
    assert len(idx) == 1, (prefix, idx)
    i = idx[0]
    if field == 'power':
        new, n = re.subn(r'(?<![a-z_])power=[^,\]]+', lambda m: f'power="{factor}"', lines[i], count=1)
    else:
        new, n = re.subn(r'(?<![a-z_])' + field + r'="([^"]+)"',
                         lambda m: f'{field}="({m[1]})*({factor})"', lines[i], count=1)
    assert n == 1, (prefix, field)
    lines[i] = new


def canonical_overall():
    table = {}
    for line in (Path(__file__).resolve().parent / 'overall_canonical.tsv').read_text().splitlines():
        if line and not line.startswith('#'):
            name, canonical = line.split('\t', 1)
            table[name] = canonical
    return table


def write(out, name, prepared, lines, overall=None):
    dst = out / name / 'L3'
    if dst.parent.exists():
        shutil.rmtree(dst.parent)
    dst.mkdir(parents=True)
    shutil.copyfile(prepared / 'model.json', dst / 'model.json')
    manifest = json.loads((prepared / 'generation.json').read_bytes())
    spec = dict(manifest['diagrams'][6])
    assert spec['name'] == 'B_self_energy_L3_006'
    spec['index'] = 0
    if overall is not None:
        spec['overall_factor'] = overall
    manifest['diagrams'] = [spec]
    manifest['report']['retained_count'] = 1
    (dst / 'generation.json').write_text(json.dumps(manifest))
    text = ''.join(lines)
    (dst / spec['dot']).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def main():
    prepared, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    dot = (prepared / 'diagram-000006.dot').read_text()
    rows = []
    for name, powers in [('base-hard-v0', {'v0': 3}), ('base-hidden', {'v0': 3, 'v1': -1})]:
        rows.append((name, 'control', '-', write(out, name, prepared, base_lines(dot, powers))))
    for attribute, (entries, powers, places, aggregate) in ATTRIBUTES.items():
        for spelling, factor in SPELLINGS.items():
            lines = base_lines(dot, powers)
            for prefix, field in places:
                times(lines, prefix, field, factor)
            if aggregate:
                times(lines, 'graph [', 'numerator', factor)
            name = f'{attribute}-{spelling}'
            overall = canonical_overall()[name] if attribute == 'overall' else None
            rows.append((name, entries, spelling, write(out, name, prepared, lines, overall)))
    # the overall-factor controls: the base (its manifest binding re-derived) and a regular neighbour
    lines = base_lines(dot, {'v0': 3})
    rows.append(('overall-none', 'overall_factor', '-', write(out, 'overall-none', prepared, lines, canonical_overall()['overall-none'])))
    lines = base_lines(dot, {'v0': 3})
    times(lines, 'graph [', 'overall_factor', SPELLINGS['gtr'].replace('-4)', '-3)'))
    rows.append(('overall-gtr3', 'overall_factor', 'regular', write(out, 'overall-gtr3', prepared, lines, canonical_overall()['overall-gtr3'])))
    # the projector control: a regular projector (1+eps)
    lines = base_lines(dot, {'v0': 3})
    times(lines, 'graph [', 'projector', f'1+{EPS}')
    rows.append(('proj-reg', 'projector', 'regular', write(out, 'proj-reg', prepared, lines)))
    for row in rows:
        print('\t'.join(row))


if __name__ == '__main__':
    main()
