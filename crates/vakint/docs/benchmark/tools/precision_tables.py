#!/usr/bin/env python3
"""Precision tables of vakint's RustRed results, generated from a vakint-external-data tree and nothing else.

  precision_tables.py --data-root ROOT --json OUT.json [--md OUT.md] [--status provisional|final] [--label TEXT]
                      [--constants-rs FILE] [--pole-depth FILE]
  precision_tables.py --data-root ROOT --verify TABLES.json [--constants-rs FILE] [--pole-depth FILE]

Table 1: every master of the tree (1L-5L) x every Laurent order the data hold; table 2: for each loop count L, every order
eps^k of the result of an arbitrary L-loop single-scale vacuum input (RustRed backend), worst case over the inputs.
Cell classes (worst first):  N(X)  numerical only, X validated digits (rounded down) of the master coefficient(s);
                             S*(X) exact expression containing finite-digit periods; X = the fewest digits among them;
                             S     exact expression in periods vakint evaluates at any requested precision.
Rule (vakint inspection.rs PreparedExpression::available_precision): a coefficient's available precision is the minimum
significant_digits of its finite sources; none -> arbitrary precision. Finite sources = the records of the tree's
constant_bindings (4L/constants.json) and numerical-only master coefficients (manifest significant_digits). Everything
inspection/constants.rs builtin_definition() resolves is evaluated at the requested precision; that list is read from
--constants-rs (default: this crate's src/inspection/constants.rs) and an unknown symbol makes the generator fail.
Table 2 (rustred_evaluation/masters.rs finalize): the order eps^k needs every master m through eps^(k + p_m), p_m = spurious-
pole order of m's reduction coefficient; vakint tries the basis-alias route first (aliased masters rewritten onto their
representatives, rebuilt here from the manifests' exact_basis_expression as masters.rs load_basis_aliases does) and falls back
to the direct route; if neither has the orders it returns MissingMasterOrders (no value). p_m: 1L-4L from --pole-depth
(measured, pinned to the reduction tables by digest; default docs/benchmark/data/pole-depth.json), 5L 0 (pole-free basis,
read from 5L/provenance/refinement-layer.json). Worst case = any combination of masters can occur in an input.
--verify regenerates and compares the content hash with TABLES.json (exit 1 if they differ): stale tables are refused.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent                      # crates/vakint/docs/benchmark
CRATE = BENCH.parent.parent              # crates/vakint
SCHEMA = 'vakint.precision-tables.v1'
MARK = {'S': '🟢', 'S*': '🟠', 'N': '🔴'}
INF = 10 ** 9


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def git(root, *a):
    r = subprocess.run(['git', '-C', str(root), *a], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else None


# ------------------------------------------------------------------------------------------------ vakint's built-ins --

def read_builtins(path: Path):
    """The constants inspection/constants.rs builtin_definition() resolves exactly (evaluated at any precision)."""
    src = path.read_text()
    try:
        body = src[src.index('pub(crate) fn builtin_definition'):src.index('pub(crate) fn builtin_symbols')]
        var_branch = body[body.index('AtomView::Var(variable)'):body.index('AtomView::Fun(function)')]
        fun_branch = body[body.index('AtomView::Fun(function)'):]
        named = re.findall(r'"([A-Za-z0-9_]+)"', re.search(r'const NAMED_CONSTANTS: &\[&str\] = &\[(.*?)\];', src, re.S).group(1))
        bars = re.findall(r'\("([A-Za-z0-9_]+)",\s*"[A-Za-z0-9_]+",\s*-?\d+\)',
                          re.search(r'MATAD_BARS: \[[^\]]*\] = \[(.*?)\];', src, re.S).group(1))
        projections = re.findall(r'"([A-Za-z]+)"', re.search(r'if \[\s*((?:"[A-Za-z]+",\s*)+)\]\s*\.iter\(\)\s*\.any', fun_branch).group(1))
    except (ValueError, AttributeError) as e:
        raise SystemExit(f'{path}: builtin_definition layout not recognised ({e}); update precision_tables.py')
    explicit = set(re.findall(r'symbol == vk_symbol!\("([A-Za-z0-9_]+)"\)', var_branch))
    if 'NAMED_CONSTANTS' not in var_branch or 'MATAD_BARS' not in var_branch:
        raise SystemExit(f'{path}: the Var branch no longer reads NAMED_CONSTANTS / MATAD_BARS; update precision_tables.py')
    need = {'zeta': 'symbol == zeta()', 'sqrt': 'Symbol::SQRT', 'log': 'Symbol::LOG', 'polylog': 'symbol == polylog()'}
    missing = [k for k, v in need.items() if v not in fun_branch]
    if missing:
        raise SystemExit(f'{path}: builtin_definition no longer handles {missing}; update precision_tables.py')
    return {'variables': sorted(set(named) | set(bars) | explicit), 'projections': sorted(projections),
            'functions': ['pi', 'e', 'zeta(n), n >= 2', 'sqrt(n), n >= 0', 'log(q), q > 0 rational, q != 1',
                          'polylog(n, q), q rational < 1 (or q = 1, n >= 2)'],
            'source': str(path), 'sha256': sha256_file(path)}


# ------------------------------------------------------------------------------------------- expression scanning --

IDENT = re.compile(r'(?<![A-Za-z0-9_:{}])((?:[A-Za-z_][A-Za-z0-9_]*::(?:\{\}::)?)*[A-Za-z_][A-Za-z0-9_]*)')
FLOAT = re.compile(r'(?<![A-Za-z0-9_])\d+\.\d|\d[eE][-+]?\d|`')
RATIONAL = re.compile(r'^\s*-?\d+(?:/\d+)?\s*$')
KNOWN_UNICODE = {'𝜋', '𝑖'}                 # Symbolica's pi and imaginary unit


def balanced(expr, i):
    """expr[i] == '(': index after the matching ')'"""
    depth = 0
    for j in range(i, len(expr)):
        if expr[j] == '(':
            depth += 1
        elif expr[j] == ')':
            depth -= 1
            if depth == 0:
                return j + 1
    raise ValueError('unbalanced parentheses')


def finite_sources(expr, B, consts, where):
    """{symbol: digits} of the finite sources in an exact row; fails on a symbol vakint would not resolve"""
    found = {}
    odd = {ch for ch in expr if ord(ch) > 127} - KNOWN_UNICODE
    if odd:
        raise SystemExit(f'{where}: unknown symbols {sorted(odd)}')
    i = 0
    while True:
        m = IDENT.search(expr, i)
        if not m:
            break
        name, j = m.group(1), m.end()
        args = ''
        if j < len(expr) and expr[j] == '(':
            k = balanced(expr, j)
            args, j = expr[j:k], k
        short = name.split('::')[-1]
        if name.startswith('vakint::'):
            if args and short in B['projections']:
                if IDENT.search(args[1:-1]):
                    raise SystemExit(f'{where}: {name}{args[:80]} has symbolic arguments')
            elif not args and short in B['variables']:
                pass
            else:
                key = f'vakint::{short}{args.replace(" ", "")}'
                if key not in consts:
                    raise SystemExit(f'{where}: {key!r} is neither a vakint built-in nor a finite constant of the tree')
                found[key] = consts[key]['digits']
        elif name == 'symbolica::zeta':
            if not (args and RATIONAL.match(args[1:-1]) and '/' not in args and int(args[1:-1]) >= 2):  # zeta(n), n >= 2
                raise SystemExit(f'{where}: zeta{args} is not a built-in')
        elif name == 'symbolica::polylog':
            n, _, q = args[1:-1].partition(',')
            if not (RATIONAL.match(n) and RATIONAL.match(q) and '/' not in n and int(n) >= 1
                    and (eval_rational(q) < 1 or (eval_rational(q) == 1 and int(n) >= 2))):
                raise SystemExit(f'{where}: polylog{args} is not a built-in')
        elif name == 'log':
            if not (args and RATIONAL.match(args[1:-1]) and eval_rational(args[1:-1]) > 0 and eval_rational(args[1:-1]) != 1):
                raise SystemExit(f'{where}: log{args} is not a built-in')
        else:
            raise SystemExit(f'{where}: unknown symbol {name!r}')
        i = j
    return found


def eval_rational(s):
    a, _, b = s.strip().partition('/')
    return int(a) / (int(b) if b else 1)


# ------------------------------------------------------------------------------------------------------- the tree --

def load_tree(root: Path, B):
    consts, labels = {}, {}
    for L in range(1, 6):
        man = json.loads((root / f'{L}L' / 'manifest.json').read_text())
        cb = man.get('constant_bindings')
        if cb:
            p = root / f'{L}L' / cb['path']
            if sha256_file(p) != cb['sha256']:
                raise SystemExit(f'{p}: sha256 differs from {L}L/manifest.json constant_bindings')
            for r in json.loads(p.read_text()):
                v = (r.get('provenance') or {}).get('validation') or {}
                consts[r['symbol'].replace(' ', '')] = {'digits': r['significant_digits'], 'loop': L,
                                                        'error_quality': r.get('error_quality')}
                lab = v.get('digits_label') or ''
                if 'two profiles of one code' in lab:
                    if v.get('independent') is False and 'NOT independent' in lab:
                        labels[r['symbol']] = {'strict': None, 'independent': v.get('independent_method_digits'),
                                               'method': v.get('independent_method')}
                    elif v.get('independent') == 'to prec_validated_digits_strict only':
                        labels[r['symbol']] = {'strict': v.get('prec_validated_digits_strict'),
                                               'independent': v.get('independent_method_digits'), 'method': v.get('independent_method')}
    masters = []
    for L in range(1, 6):
        man = json.loads((root / f'{L}L' / 'manifest.json').read_text())
        for m in man['masters']:
            sym = {}
            f = (m.get('files') or {}).get('symbolic.json')
            if f:
                p = root / f'{L}L' / f['path']
                if sha256_file(p) != f['sha256']:
                    raise SystemExit(f'{p}: sha256 differs from the manifest')
                sym = {o: e for o, e in json.loads(p.read_text())}
            orders = {}
            for c in m['coefficients']:
                o, where = c['epsilon_power'], f"{L}L/{m['name']} eps^{c['epsilon_power']}"
                if c.get('symbolic'):
                    if o not in sym:
                        raise SystemExit(f'{where}: symbolic coefficient without a row')
                    src = finite_sources(sym[o], B, consts, where)
                    if FLOAT.search(sym[o]):        # bound as a finite literal at the coefficient's digits (bind_finite_literals)
                        src[f'literal in {where}'] = c['significant_digits']
                    cell = {'class': 'S*', 'X': min(src.values()), 'constants': src} if src else {'class': 'S'}
                else:
                    if not c.get('numerical') or not c.get('significant_digits'):
                        raise SystemExit(f'{where}: coefficient with neither a row nor validated digits')
                    cell = {'class': 'N', 'X': int(c['significant_digits'])}
                cell['payload_digits'] = c.get('significant_digits') if c.get('numerical') else None
                orders[o] = cell
            masters.append({'loop': L, 'name': m['name'], 'master_id': m['master_id'], 'expr': m.get('exact_basis_expression'),
                            'orders': orders, 'family': m['family_sha256'], 'powers': m['powers']})
    return consts, labels, masters


def rank(c):
    """larger = worse"""
    if c is None:
        return (3, 0)
    return (0, 0) if c['class'] == 'S' else (1 if c['class'] == 'S*' else 2, -c['X'])


def label(c):
    return '—' if c is None else c['class'] if c['class'] == 'S' else f"{c['class']}({c['X']})"


# ------------------------------------------------------------------------------------- basis aliases (masters.rs) --

VAR = re.compile(r'^(?:[A-Za-z_][A-Za-z0-9_]*::(?:\{\}::)?)*[A-Za-z_][A-Za-z0-9_]*$')


def aliases(masters):
    """masters.rs load_basis_aliases: identical complete expressions -> first master; expressions rational in ep and in
    primitive symbols (bare-variable expressions) -> that relation on the primitives' first masters. Returns
    {master_id: [representative master_ids]} for aliased masters."""
    out, identical, prim = {}, {}, {}
    for m in masters:
        e = m['expr']
        if not e or not m['orders']:
            continue
        rep = identical.setdefault((m['loop'], e), m['master_id'])
        if rep != m['master_id']:
            out[m['master_id']] = [rep]
        bare = e.strip()
        if VAR.match(bare) and bare.split('::')[-1] != 'ep':
            prim.setdefault((m['loop'], bare), m['master_id'])
    for m in masters:
        e = m['expr']
        if not e or re.search(r'[A-Za-z_][A-Za-z0-9_:{}]*\s*\(', e) or FLOAT.search(e):
            continue                                  # a function call or a float: not an exact rational relation
        names = set(IDENT.findall(e))
        if not names or not all(n.split('::')[-1] == 'ep' or (m['loop'], n) in prim for n in names):
            continue
        reps = sorted({prim[(m['loop'], n)] for n in names if n.split('::')[-1] != 'ep'})
        if reps != [m['master_id']]:
            out[m['master_id']] = reps
    return out


# ---------------------------------------------------------------------------------------------- spurious poles --

def framed_4l(root, n):
    prog = gzip.decompress((root / '4L' / 'ibps' / f'{n}.candidates.rrbin.gz').read_bytes())
    norm = (root / '4L' / 'ibps' / f'{n}.rrnorm.bin').read_bytes()
    h = hashlib.sha256(b'vakint-native-candidate-and-normalization-v1\0')
    for b in (prog, norm):
        h.update(struct.pack('<Q', len(b)))
        h.update(b)
    return h.hexdigest()


def load_depths(root: Path, path: Path, masters, alias_map):
    d = json.loads(path.read_text())
    if d.get('schema') != 'vakint.pole-depth.v1':
        raise SystemExit(f'{path}: unknown schema')
    for L, rec in d['measured_on'].items():
        if 'rules_rrbin_sha256' in rec:
            now = sha256_file(root / L / 'ibps' / 'rules.rrbin')
            if now != rec['rules_rrbin_sha256']:
                raise SystemExit(f'{L}/ibps/rules.rrbin differs from the reduction table whose poles were measured ({path})')
        for n, want in (rec.get('framed_program_sha256') or {}).items():
            if framed_4l(root, n) != want:
                raise SystemExit(f'4L/ibps/{n} differs from the reduction program whose poles were measured ({path})')
    depth = {}
    for m in masters:
        if m['loop'] == 5:
            continue
        r = d['masters'].get(m['master_id'])
        if r is None:
            raise SystemExit(f"{m['loop']}L/{m['name']}: no measured pole order in {path}")
        depth[m['master_id']] = r
        want = sorted(r.get('aliased_to') or [])
        have = sorted(alias_map.get(m['master_id'], []))
        if want != have:          # the measured alias route must be the one vakint takes on this tree
            raise SystemExit(f"{m['loop']}L/{m['name']}: alias representatives {have} differ from the measured ones {want}")
    rl = json.loads((root / '5L' / 'provenance' / 'refinement-layer.json').read_text())
    bc = rl.get('basis_change') or {}
    gp = (bc.get('checks') or {}).get('g_pole') or {}
    if not (gp.get('passed') is True and gp.get('worst_valuation') == 0 and gp.get('n_polar') == 0
            and 'regular at d = 4' in (bc.get('design') or '')):
        raise SystemExit('5L/provenance/refinement-layer.json does not record a pole-free basis (basis_change.checks.g_pole)')
    proof = re.search(r'exact proof: (\d+) / (\d+) generators regular at d = 4', bc.get('design') or '')
    five = {'tag': bc.get('tag'), 'g_pole_coefficients': gp.get('coefficients'),
            'proof': f'{proof.group(1)} / {proof.group(2)} generators' if proof else None}
    return depth, d, five


# ---------------------------------------------------------------------------------------------------- tables --

def needs(m, through):
    """the cells of master m for orders first..through, or None if one is not held"""
    o = m['orders']
    if not o:
        return None
    lo = min(o)
    if through < lo:
        return []
    cells = []
    for j in range(lo, through + 1):
        if j not in o:
            return None
        cells.append((m, j, o[j]))
    return cells


def worst_of(entries):
    w, who = None, []
    for m, j, c in entries:
        if w is None or rank(c) > rank(w):
            w, who = c, [(m, j, c)]
        elif rank(c) == rank(w) and c['class'] != 'S':
            who.append((m, j, c))
    return w, who


def table2_cell(L, k, ms, by_id, alias_map, depth):
    """worst case over inputs of the L-loop result at eps^k (see the module docstring)"""
    info = []
    for m in ms:
        mid = m['master_id']
        if L == 5:
            kd = ka_self = 0
        else:
            r = depth[mid]
            kd = r['direct']
            ka_self = r.get('alias_representative', kd)
        aliased = mid in alias_map
        if aliased:
            A = []
            for rid in alias_map[mid]:
                rr = depth[rid]
                n = needs(by_id[rid], k + rr.get('alias_representative', rr['direct']))
                if n is None:
                    A = None
                    break
                A += [(e[0], e[1], e[2], 'alias') for e in n]
        else:
            n = needs(m, k + ka_self)          # unchanged by the normalization: the record's own (alias-measured) order
            A = None if n is None else [(e[0], e[1], e[2], 'route') for e in n]
        n = needs(m, k + kd)
        D = None if n is None else [(e[0], e[1], e[2], 'direct') for e in n]
        info.append({'m': m, 'aliased': aliased, 'A': A, 'D': D})
    # alone: the master by itself gets no value; breakers: aliased masters whose alias route lacks orders but whose direct
    # route has them (an input containing one falls back to the direct route for ALL its masters); combo: masters whose direct
    # route lacks orders, so an input combining them with a breaker gets no value
    alone = [i for i in info if i['A'] is None and (not i['aliased'] or i['D'] is None)]
    breakers = [i for i in info if i['aliased'] and i['A'] is None and i['D'] is not None]
    combo = [i for i in info if i['D'] is None and i not in alone] if breakers else []
    entries = [e for i in info if i['A'] is not None for e in i['A']]
    if breakers:
        entries += [e for i in info if i['D'] is not None for e in i['D']]
    w, who = worst_of([(m, j, c) for m, j, c, _ in entries]) if entries else (None, [])
    routes = {}
    for m, j, c, r in entries:
        routes.setdefault((m['master_id'], j), r)
    return {'missing': [i['m'] for i in alone], 'combo': [i['m'] for i in combo], 'breakers': [i['m'] for i in breakers],
            'returns': bool(entries), 'worst': w,
            'who': [(m, j, c, routes.get((m['master_id'], j))) for m, j, c in who]}


def build(root: Path, B, depth_path: Path, status: str, label_text: str | None):
    consts, labels, masters = load_tree(root, B)
    by_id = {m['master_id']: m for m in masters}
    checks = tier_check(root, by_id)
    alias_map = aliases(masters)
    depth, depth_doc, five = load_depths(root, depth_path, masters, alias_map)
    all_orders = sorted({o for m in masters for o in m['orders']})
    orders = list(range(all_orders[0], all_orders[-1] + 1))

    # coefficient groups (registry.rs coefficient_group): same loop + identical exact_basis_expression; Numerical mode takes
    # the best of the group's exact formulas (cap = their finite digits) and numerical payloads (their digits)
    groups = {}
    for m in masters:
        if m['expr']:
            groups.setdefault((m['loop'], m['expr']), []).append(m)

    def numerical_digits(m, o):
        best = 0
        for g in groups.get((m['loop'], m['expr']), [m]) if m['expr'] else [m]:
            c = g['orders'].get(o)
            if c is None:
                continue
            if c['class'] == 'S':
                return INF
            if c['class'] == 'S*':
                best = max(best, c['X'])
            if c.get('payload_digits'):
                best = max(best, int(c['payload_digits']))
        return best

    def confirmed(s):
        lab = labels[s]
        return lab['strict'] if lab['strict'] is not None else lab['independent']

    def two_profile(cell):
        """constants of an S* cell whose digits up to the cell's X rest (in part) on two profiles of one code: their
        independently confirmed digits are below X"""
        return sorted(s for s in (cell.get('constants') or {}) if s in labels and (confirmed(s) or 0) < cell['X'])

    # ------------------------------------------------------------------------------------------------ table 1 --
    t1 = []
    tp_users, nm_users = {}, []
    for m in masters:
        cells = {}
        for o in orders:
            c = m['orders'].get(o)
            if c is None:
                cells[str(o)] = None
                continue
            e = {'class': c['class'], 'X': c.get('X'), 'marks': []}
            if c['class'] == 'S*':
                e['constants'] = sorted(c['constants'], key=lambda s: (c['constants'][s], s))[:6]
                tp = two_profile(c)
                if tp:
                    e['marks'].append('p')
                    for s in tp:
                        tp_users.setdefault(s, []).append(f"{m['loop']}L/{m['name']} ε^{o}")
                nd = numerical_digits(m, o)
                if nd > c['X']:
                    e['marks'].append('n')
                    e['numerical_mode_digits'] = nd
                    nm_users.append((f"{m['loop']}L/{m['name']} ε^{o}", c['X'], nd))
            cells[str(o)] = e
        t1.append({'loop': m['loop'], 'name': m['name'], 'master_id': m['master_id'], 'cells': cells,
                   **({'name_640': m['name_640']} if m.get('name_640') else {})})

    # ------------------------------------------------------------------------------------------------ table 2 --
    # columns: up to the highest order any L-loop master holds (that master itself, as an input, returns it)
    t2, notes2 = [], []
    k_max = {}
    for L in range(1, 6):
        ms = [m for m in masters if m['loop'] == L]
        k_max[L] = max(o for m in ms for o in m['orders'])
        t2.append((L, {k: table2_cell(L, k, ms, by_id, alias_map, depth) for k in range(-L, k_max[L] + 1)}))
    orders2 = list(range(-5, max(k_max.values()) + 1))

    def note(text):
        notes2.append(text)
        return chr(ord('a') + len(notes2) - 1)

    def masters_text(items, limit=6):
        names = []
        for x in items:
            s = f"{x['loop']}L/{x['name']}"
            if s not in names:
                names.append(s)
        return ', '.join(names[:limit]) + (f' and {len(names) - limit} more' if len(names) > limit else '')

    t2_rows = []
    for L, row in t2:
        cells = {}
        for k in orders2:
            if k < -L:
                cells[str(k)] = {'class': None, 'na': 'below'}
                continue
            if k not in row:
                cells[str(k)] = {'class': None, 'na': 'beyond'}
                continue
            ms = [x for x in masters if x['loop'] == L]
            c = row[k]
            w, who = c['worst'], c['who']
            if not c['returns'] and not any(k in x['orders'] for x in masters if x['loop'] == L):
                cells[str(k)] = {'class': None, 'na': 'beyond'}
                continue
            partial = bool(c['missing'] or c['combo'] or w is None)
            e = {'class': None if partial else w['class'], 'X': None if partial else w.get('X'), 'marks': []}
            binding = sorted({f"{m['loop']}L/{m['name']}" + (f" ({m['name_640']})" if m.get('name_640') else '') + f" ε^{j}"
                              + (' (alias route)' if r == 'alias' else '') for m, j, cc, r in who})
            if partial:
                e['na'] = 'partial'
                e['returned'] = label(w)
                e['missing_masters'] = [f"{x['loop']}L/{x['name']}" for x in c['missing']]
                e['missing_in_combination'] = [f"{x['loop']}L/{x['name']}" for x in c['combo']]
                txt = f"{L}L ε^{k}: not available for every input. "
                if c['missing']:
                    txt += f"Inputs reaching {masters_text(c['missing'])} get MissingMasterOrders (no value; vakint fails closed). "
                if c['combo']:
                    txt += (f"So do inputs that combine {masters_text(c['combo'])} with {masters_text(c['breakers'])} (whose alias "
                            f"route lacks orders, so vakint takes the direct route, which lacks the former's orders). ")
                if w is not None:
                    txt += f"The other inputs return {label(w)}" + (f" (set by {', '.join(binding[:4])})" if w['class'] != 'S' else '') + '.'
                else:
                    own = [(x, x['orders'][k]) for x in ms if k in x['orders']]
                    ow, _ = worst_of([(x, j, x['orders'][j]) for x, _ in own for j in x['orders'] if j <= k])
                    txt += (f"Only inputs whose reduction coefficients have lower spurious-pole orders return it, e.g. the "
                            f"{len(own)} masters that hold ε^{k} ({masters_text([x for x, _ in own], 4)}) themselves: {label(ow)}.")
                    e['returned'] = label(ow)
                e['marks'].append(note(txt))
            elif w['class'] != 'S':
                pole = [x for x in who if x[1] > k]
                why = f"set by {', '.join(binding[:4])}{' and more' if len(binding) > 4 else ''}"
                if w['class'] == 'S*':
                    cs = sorted({s for m, j, cc, r in who for s, dg in (cc.get('constants') or {}).items() if dg == w['X']})
                    why += f"; constant {', '.join(cs)}"
                    tp = sorted({s for m, j, cc, r in who for s in two_profile(cc)})
                    if tp:
                        why += ('; ' + '; '.join(
                            f"{s} rests on two profiles of one code ({consts[s]['digits']} d; independent confirmation "
                            f"{confirmed(s)} d)" for s in tp))
                if pole:
                    why += (f"; spurious poles up to 1/ε^{max(j - k for m, j, cc, r in pole)} make ε^{k} need ε^{max(j for m, j, cc, r in pole)} "
                            "of these masters (the physical pole coefficients are polylogarithmic by the R-operation, so these "
                            "constants should cancel; the data do not guarantee it)")
                e['binding'] = binding[:12]
                e['marks'].append(note(f"{L}L ε^{k} {label(w)}: {why}."))
            cells[str(k)] = e
        t2_rows.append({'loop': L, 'cells': cells, 'last_order': k_max[L]})

    # ------------------------------------------------------------------------------------------------ notes --
    tp_list = []
    for s, users in sorted(tp_users.items()):
        lab = labels[s]
        conf = (f"independent implementations agree to {lab['strict']} d" if lab['strict'] is not None else
                f"independent confirmation {lab['independent']} d" if lab['independent'] is not None else 'no independent confirmation')
        tp_list.append(f"{s} ({consts[s]['digits']} d; {conf}): {len(users)} order{'s' if len(users) != 1 else ''}")
    counts = {}
    for r in t1:
        for o, c in r['cells'].items():
            if c:
                key = f"{r['loop']}L"
                d = counts.setdefault(key, {'S': 0, 'S*': 0, 'N': 0})
                d[c['class']] += 1
    sstar_min = {}
    for r in t1:
        for o, c in r['cells'].items():
            if c and c['class'] in ('S*', 'N'):
                k = (f"{r['loop']}L", c['class'])
                sstar_min[k] = min(sstar_min.get(k, INF), c['X'])

    # identities behind 'S' (data provenance of the tree)
    ident = identity_notes(root)
    pslq = pslq_notes(root)

    rm = root / 'manifest.json'
    sources = {
        'data_root_manifest_sha256': sha256_file(rm),
        'loop_manifest_sha256': {f'{L}L': sha256_file(root / f'{L}L' / 'manifest.json') for L in range(1, 6)},
        'constants_sha256': {},
        'five_loop_values': five_loop_build(root),
        'version': json.loads(rm.read_text()).get('version'),
        'vakint_constants_rs_sha256': B['sha256'],
        'pole_depth_sha256': sha256_file(depth_path),
    }
    for L in range(1, 6):
        cb = json.loads((root / f'{L}L' / 'manifest.json').read_text()).get('constant_bindings')
        if cb:
            sources['constants_sha256'][f'{L}L'] = cb['sha256']
    content = {
        'schema': SCHEMA, 'sources': sources, 'orders_table1': orders, 'orders_table2': orders2,
        'table1': t1, 'table2': t2_rows, 'notes_table2': notes2,
        'two_profile_constants': tp_list, 'numerical_mode_more_digits': [list(x) for x in nm_users],
        'counts_table1': counts, 'worst_table1': {f'{a} {b}': v for (a, b), v in sorted(sstar_min.items())},
        'rule': {'builtins': {k: B[k] for k in ('variables', 'projections', 'functions')},
                 'finite_constants': {s: v['digits'] for s, v in sorted(consts.items())}},
        'pole_depth': {'method': depth_doc['method'], 'five_loop': five,
                       'orders': {f'{L}L': {'direct': sorted({x['direct'] for x in depth_doc['masters'].values() if x['loop'] == L}),
                                            'alias_representative': sorted({x['alias_representative'] for x in depth_doc['masters'].values()
                                                                            if x['loop'] == L and 'alias_representative' in x})}
                                  for L in range(1, 5)},
                       'not_reduced_samples': {L: s for L, v in depth_doc['sources'].items() for f, s in v.items() if 'not reduced' in s}},
        'identities': ident, 'pslq': pslq, 'checks': checks,
    }
    content['method'] = method_notes(content)
    dg = None
    if (root / '.git').exists():
        dg = {'commit': git(root, 'rev-parse', 'HEAD'), 'tree': git(root, 'rev-parse', 'HEAD^{tree}'),
              'dirty_files': len([l for l in (git(root, 'status', '--porcelain', '--untracked-files=normal') or '').splitlines() if l])}
    if status == 'final' and not (dg and dg['dirty_files'] == 0):
        raise SystemExit('--status final needs a clean git checkout of the published data tree')
    return {**content, 'content_sha256': content_hash(content), 'status': status,
            'generated_utc': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            'data_root': str(root), 'data_git': dg, 'label': label_text,
            'generator': {'path': 'crates/vakint/docs/benchmark/tools/precision_tables.py', 'sha256': sha256_file(Path(__file__))}}


def tier_check(root, by_id):
    """5L: the release values' tiers (5L/producer-update.json value_provenance) must map onto the classes read from the rows:
    S -> S, S* -> S* (same finite constants), numerical tiers (A, B, F, ...) -> N with X <= the tier's digits"""
    d = json.loads((root / '5L' / 'producer-update.json').read_text())
    stale, tiers = [], {}
    for mid, v in d['masters'].items():
        m = by_id.get(mid)
        if m is None:
            raise SystemExit(f'5L/producer-update.json: master {mid[:40]} is not in the 5L manifest')
        m['name_640'] = v.get('name_640')
        for o, r in (v.get('value_provenance') or {}).items():
            c, t = m['orders'].get(int(o)), r.get('tier')
            tiers[t] = tiers.get(t, 0) + 1
            where = f"5L/{m['name']} ε^{o}"
            if c is None:
                raise SystemExit(f'{where}: tier {t} in producer-update.json but no coefficient in the manifest')
            want = 'S' if t == 'S' else 'S*' if t == 'S*' else 'N'
            if c['class'] != want:
                raise SystemExit(f"{where}: release tier {t} but the row is {label(c)}")
            if want == 'N' and c['X'] > float(r['digits']):
                raise SystemExit(f"{where}: manifest digits {c['X']} exceed the release tier's {r['digits']}")
            if want == 'S*' and r.get('sstar_constants'):
                rec = dict((a.replace(' ', ''), b) for a, b in (json.loads(r['sstar_constants'].replace("'", '"'))
                                                                if isinstance(r['sstar_constants'], str) else r['sstar_constants']))
                if set(rec) != set(c['constants']):
                    raise SystemExit(f"{where}: S* constants {sorted(rec)} of the release record differ from the row's {sorted(c['constants'])}")
                if any(rec[k] != c['constants'][k] for k in rec):
                    stale.append(where)
    return {'five_loop_tiers': dict(sorted(tiers.items())), 'five_loop_sstar_digit_records_differing_from_4L_constants': stale}


NOT_CONTENT = {'content_sha256', 'status', 'generated_utc', 'data_root', 'data_git', 'label', 'generator'}


def content_hash(t):
    """sha256 of everything the tables say (not where, when or by which file they were written)"""
    c = {k: v for k, v in t.items() if k not in NOT_CONTENT}
    return hashlib.sha256(json.dumps(c, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def five_loop_build(root):
    p = root / '5L' / 'producer-update.json'
    d = json.loads(p.read_text())
    prod = (d.get('provenance') or {}).get('producer', '')
    m = re.findall(r'(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)', prod)
    return {'producer_update_sha256': sha256_file(p), 'build_utc': m[0] if m else None,
            'basis': (d.get('provenance') or {}).get('basis')}


def identity_notes(root):
    """numerical (not proven) identities behind basis-v4 rows: provenance/takeover-basis-v4-rewrite-v1/receipt.json (status,
    checked digits), updated by later checks the tree records (precision-lift identity-checks.json, independence audit)"""
    p = root / 'provenance' / 'takeover-basis-v4-rewrite-v1' / 'receipt.json'
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    rel = {r['id']: r for r in d['equivalence']['relations_used']}
    numerical = {i for i, r in rel.items()
                 if r['status'].startswith(('COMPUTED', 'ACCEPTED-')) or re.search(r'PSLQ|lindep|integer relation', r['status'])}
    use = {}
    for mid, x in d['rewrites'].items():
        for o, row in x['rows'].items():
            for rid in row.get('rules') or []:
                if rid in numerical:
                    use.setdefault(rid, {}).setdefault(f"{x['loop_count']}L", set()).add((x['name'], int(o)))
    later = {}
    ic = root / 'provenance' / 'precision-lift-4l-runtime-v4' / 'identity-checks.json'
    if ic.exists():
        for rid, c in (json.loads(ic.read_text()).get('checks') or {}).items():
            if c.get('verdict') == 'PASS':
                later[rid] = (c['agreement_digits'], 'at the shipped constants (precision lift)')
    for f in sorted((root / 'provenance' / 'precision-lift-4l-runtime-v4' / 'objects').glob('*.md.gz')):
        t = gzip.decompress(f.read_bytes()).decode('utf-8', 'replace')
        mm = re.search(r'X3 and X4 \(the Im G swaps\) re-found blind by two codes at (\d+) d', t)
        if mm:
            for rid in ('X3_TI3', 'X4_TI4'):
                later[rid] = (float(mm.group(1)), 'independence audit archived with the precision lift')
    out = {}
    for rid, u in sorted(use.items()):
        st = rel[rid]['status']
        dig = [float(x) for x in re.findall(r'(\d{3,5}(?:\.\d+)?) digits', st)]
        v = {'status': st.split(' (')[0], 'rows': {L: len(v) for L, v in sorted(u.items())},
             'masters': {L: len({n for n, _ in v}) for L, v in sorted(u.items())},
             'lowest_order': {L: min(o for _, o in v) for L, v in sorted(u.items())},
             'checked_digits': max(dig, default=None)}
        if rid in later and (v['checked_digits'] is None or later[rid][0] > v['checked_digits']):
            v['rechecked_digits'], v['rechecked_by'] = later[rid]
        out[rid] = v
    return out


def pslq_notes(root):
    p = root / '5L' / 'producer-update.json'
    d = json.loads(p.read_text())
    n, exact = 0, 0
    for m in (d['masters'].values() if isinstance(d['masters'], dict) else d['masters']):
        for o, v in (m.get('value_provenance') or {}).items():
            if v.get('tier') in ('S', 'S*'):
                exact += 1
                if (v.get('bar') or {}).get('kind') in ('pslq', 'symp-v4'):
                    n += 1
    return {'exact_orders': exact, 'pslq_orders': n}


# ---------------------------------------------------------------------------------------------------- markdown --

def cell_md(c):
    if c is None:
        return '—'
    if c.get('class') is None:
        if c.get('na') == 'partial':
            return '— [' + ','.join(c['marks']) + ']'
        return ' ' if c.get('na') == 'below' else '—'
    t = f"{MARK[c['class']]} {c['class'] if c['class'] == 'S' else c['class'] + '(' + str(c['X']) + ')'}"
    return t + (' [' + ','.join(c['marks']) + ']' if c.get('marks') else '')


def eps(o):
    return f'ε^{o}'


def markdown(T):
    s = T['sources']
    f5 = s['five_loop_values']
    head = (f"**{T['status'].upper()}** · data tree {T.get('label') or T['data_root']}"
            + (f" (git {T['data_git']['commit'][:12]}, tree {T['data_git']['tree'][:12]}"
               f"{', ' + str(T['data_git']['dirty_files']) + ' dirty files' if T['data_git']['dirty_files'] else ''})" if T.get('data_git') else '')
            + f" · {s['version']} root manifest `{s['data_root_manifest_sha256'][:16]}` · 4L constants `{s['constants_sha256'].get('4L', '')[:16]}`"
            + f" · 5L values build {f5['build_utc']} (producer-update `{f5['producer_update_sha256'][:16]}`)"
            + f" · generated {T['generated_utc']} · tables `{T['content_sha256'][:16]}`")
    out = ['# Precision of vakint’s master coefficients and results', '', head, '',
           '🟢 **S** exact, in periods vakint evaluates at any requested precision · '
           '🟠 **S\\*(X)** exact, but contains periods vakint has to at most X digits · '
           '🔴 **N(X)** numerical only, at most X digits · — not available. X = worst-case digits, rounded down.', '']
    o1 = T['orders_table1']
    out += ['## Table 1: master coefficients (rows: masters, columns: Laurent orders the data hold)', '',
            '| master | ' + ' | '.join(eps(o) for o in o1) + ' |', '|---|' + '---|' * len(o1)]
    last = None
    for r in T['table1']:
        if r['loop'] != last:
            c = T['counts_table1'][f"{r['loop']}L"]
            n = sum(1 for x in T['table1'] if x['loop'] == r['loop'])
            out.append(f"| **{r['loop']}L** ({n} master{'s' if n != 1 else ''}: "
                       f"{c['S']} S, {c['S*']} S\\*, {c['N']} N) |" + ' |' * len(o1))
            last = r['loop']
        out.append(f"| {r['loop']}L/{r['name']}{' (' + r['name_640'] + ')' if r.get('name_640') else ''} | "
                   + ' | '.join(cell_md(r['cells'][str(o)]) for o in o1) + ' |')
    out += ['', '[p] Part of the X digits rests on a constant computed with two profiles of one code beyond its independent '
            'confirmation (constant, its digits, its independent confirmation, orders affected): '
            + '; '.join(T['two_profile_constants']) + '.',
            '[n] vakint’s Numerical mode returns more digits for these orders: it uses the best payload of the master’s coefficient '
            'group (masters with the same exact basis expression), here a validated numerical value: '
            + '; '.join(f'{w} S*({x}) -> {d} d' for w, x, d in T['numerical_mode_more_digits']) + '.', '']
    out = [x.replace(' eps^', ' ε^') for x in out]
    o2 = T['orders_table2']
    out += ['## Table 2: L-loop result of any single-scale vacuum input (rows: L, columns: Laurent orders of the result)', '',
            '| L | ' + ' | '.join(eps(o) for o in o2) + ' |', '|---|' + '---|' * len(o2)]
    for r in T['table2']:
        out.append(f"| {r['loop']}L | " + ' | '.join(cell_md(r['cells'][str(o)]) for o in o2) + ' |')
    out += [''] + [f"[{chr(ord('a') + i)}] {n}" for i, n in enumerate(T['notes_table2'])] + ['']
    out += ['## How the tables are made', ''] + [f'- {x}' for x in T['method']] + ['']
    return '\n'.join(out)


def method_notes(T):
    B = T['rule']['builtins']
    pd = T['pole_depth']
    five = pd['five_loop']
    ids = T.get('identities') or {}
    fc = T['rule']['finite_constants']
    lo = min(fc.values()) if fc else None
    notes = [
        'Rule (vakint PreparedExpression::available_precision): a coefficient is limited by the fewest significant digits of its '
        'finite sources; an expression without finite sources is evaluated at any requested precision. Built-ins (evaluated at '
        f"any precision, read from inspection/constants.rs): {', '.join(B['variables'])}; the projections "
        f"{', '.join(B['projections'])} (HPL/GPL words); {'; '.join(B['functions'])}. Finite sources: the {len(fc)} records of "
        f"4L/constants.json ({lo}-{max(fc.values())} d) and numerical-only master coefficients.",
        'Finite-digit periods (significant digits as vakint reads them): '
        + ', '.join(f"{k.replace('vakint::', '')} {v}" for k, v in sorted(fc.items(), key=lambda kv: (kv[1], kv[0]))) + '.',
        'Table 1 classifies each stored order: an exact row without finite sources is S, with finite sources S*(min digits), a '
        'row without an exact form N(validated digits). Symbolic-mode class; numerical-only values come back with X-1 to X-3 '
        'tracked digits.',
        'Table 2: the result at ε^k needs every master m through ε^(k+p_m), p_m = spurious-pole order of m’s reduction '
        'coefficient. vakint uses its basis-alias route first (masters with the same exact value rewritten onto one representative) '
        'and falls back to the direct route; the class is the worst over all master orders a route needs, over all inputs (any '
        'combination of masters). A needed order that is not held gives MissingMasterOrders: no value, never a wrong one.',
        f"Spurious-pole orders 1L-4L (direct route {', '.join(f'{L} {max(v['direct'])}' for L, v in pd['orders'].items())}; alias route "
        f"{', '.join(f'{L} {max(v['alias_representative'])}' for L, v in pd['orders'].items() if v['alias_representative'])}) are "
        "measured on every record of the shipped caches and on exact samples outside them (up to 6 dots and 3 numerators at 4L; "
        "no growth with dots or numerators), not proven; they are pinned to the reduction tables by digest. 5L: pole-free basis "
        f"({five.get('tag')}; exact proof on {five.get('proof')}, {five.get('g_pole_coefficients')} sampled coefficients without a pole), "
        "p = 0 for every 5L master.",
    ]
    if ids:
        parts = []
        for rid, v in ids.items():
            parts.append(f"{rid} ({', '.join(f'{L} {n} row{chr(115) if n != 1 else chr(32)[:0]} from ε^{v['lowest_order'][L]}' for L, n in v['rows'].items())}"
                         + (f"; checked to {v['checked_digits']:.0f} d" if v.get('checked_digits') else '')
                         + (f", re-checked to {v['rechecked_digits']:.0f} d ({v['rechecked_by']})" if v.get('rechecked_digits') else '') + ')')
        notes.append('Some 1L-4L S rows (and the 5L product rows built on them) rest on numerical identities of the constant basis v4, '
                     'checked to thousands of digits but not proven: ' + '; '.join(parts) + '. Above those digits such an S is exact '
                     'only modulo the identity; vakint gives no warning for it.')
    notes.append(f"5L S includes {T['pslq']['pslq_orders']} of {T['pslq']['exact_orders']} exact orders identified by PSLQ that pass "
                 'the 40-excess-digit bar (numerical identification, not proofs).')
    notes.append('Digits lost to cancellations in the reduction coefficients are tracked and reported by vakint (PrecisionWarning), not '
                 'included in X.')
    ch = T.get('checks') or {}
    if ch.get('five_loop_tiers'):
        st = ch.get('five_loop_sstar_digit_records_differing_from_4L_constants') or []
        notes.append('5L release tiers (5L/producer-update.json) checked against the rows: '
                     + ', '.join(f'{k} {v}' for k, v in ch['five_loop_tiers'].items())
                     + ' (S -> S, S* -> S* with the same constants, numerical tiers -> N)'
                     + (f'; {len(st)} S* records there still carry older digits of their constants than 4L/constants.json, which is '
                        'what vakint uses and what the tables use' if st else '') + '.')
    notes.append('RustRed backend only (vakint’s first backend); MATAD/FMFT/alphaLoop fallbacks are not part of these tables.')
    return notes


# ---------------------------------------------------------------------------------------------------- main --

def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--json', type=Path)
    ap.add_argument('--md', type=Path)
    ap.add_argument('--status', choices=['provisional', 'final'], default='provisional')
    ap.add_argument('--label')
    ap.add_argument('--constants-rs', type=Path, default=CRATE / 'src' / 'inspection' / 'constants.rs')
    ap.add_argument('--pole-depth', type=Path, default=BENCH / 'data' / 'pole-depth.json')
    ap.add_argument('--verify', type=Path, help='compare a committed tables JSON with a regeneration from --data-root')
    a = ap.parse_args(argv)
    root = a.data_root.resolve()
    B = read_builtins(a.constants_rs)
    if a.verify:
        old = json.loads(a.verify.read_text())
        new = build(root, B, a.pole_depth, 'provisional', None)
        if content_hash(old) != old.get('content_sha256'):
            print(f"precision tables MODIFIED: {a.verify} content differs from its own content_sha256")
            return 1
        if old.get('content_sha256') != new['content_sha256']:
            print(f"precision tables STALE: {a.verify} content {str(old.get('content_sha256'))[:16]} != regenerated "
                  f"{new['content_sha256'][:16]} from {root}")
            return 1
        print(f"precision tables OK: {a.verify} = regeneration from {root} (content {new['content_sha256'][:16]}, {old.get('status')})")
        return 0
    T = build(root, B, a.pole_depth, a.status, a.label)
    if a.json:
        a.json.write_text(json.dumps(T, indent=1, ensure_ascii=False) + '\n')
    if a.md:
        a.md.write_text(markdown(T) + '\n')
    print(f"precision tables {T['status']}: content {T['content_sha256'][:16]}; "
          + '; '.join(f"{r['loop']}L " + ' '.join(cell_md(r['cells'][str(o)]).replace(' ', '') for o in T['orders_table2'] if o >= -r['loop'])
                      for r in T['table2']))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
