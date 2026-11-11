#!/usr/bin/env python3
"""Overview of what vakint computes today, symbolically or numerically, and with what precision (1L-5L single-scale vacuum).

  overview.py --data-root ROOT --tables TABLES.json --json OUT.json [--md OUT.md] [--overlay OVERLAY.json] [--names NAMES.json]
              [--vakint-src DIR] [--pole-depth FILE]

Input modes (the PDF says which one it shows):
  RELEASED  (no --overlay): only the data tree ROOT (TABLES.json = tools/precision_tables.py on the same tree; checked by digest).
  PREVIEW   (--overlay): the released tree plus derived results that no released tree holds yet (schema vakint.precision-overlay.v1:
            master, order, status, exact form, digits, digits of numerical agreement, provenance, review state; optionally a basis
            change at one loop order: kept / added / dropped masters). Overlay results are 'pending': shown with their own marker and
            counted apart from the released ones, never as released.
Sections (rendered by docs/overview.typ, also written as markdown):
  (a) summary per loop and Laurent order: coefficients S / S* / N / structural zero (exact 0 row), digits of S* and N;
  (b) per-master table (table 1 of the precision tables + zero class + pending cells);
  (c) the result of any single-scale vacuum input at L loops, per order (table 2 of the precision tables; PREVIEW adds the row of the
      overlay's basis and the dropped former masters, which fail closed at an order while numerical coefficients enter it);
  (d) named constants: standard definition, digits available, evaluated at any precision or not, how many coefficients use them;
      NAMES.json (optional) supplies traditional names / formulas: {code: {display?, traditional, terms?, status, source}};
  (e) how the symbolic forms were obtained and checked.
Classes: S = exact, periods at any precision; S*(X) = exact, contains periods known to X digits; N(X) = numerical only, X digits;
0 = structural zero (exact zero row; counted apart from S).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import precision_tables as pt  # noqa: E402  (same directory: the tree loader and route model of the precision tables)

SCHEMA = 'vakint.benchmark-overview.v1'
OVERLAY_SCHEMA = 'vakint.precision-overlay.v1'
RANK = {'S': 0, 'Z': 0, 'S*': 1, 'N': 2}


def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rank(c):
    """larger = worse (None = not available)"""
    if c is None:
        return (3, 0)
    return (RANK[c['class']], -(c.get('X') or 0) if c['class'] in ('S*', 'N') else 0)


def label(c):
    if c is None:
        return '—'
    return {'S': 'S', 'Z': '0'}.get(c['class']) or f"{c['class']}({c['X']})"


# ------------------------------------------------------------------------------------------------------ zero rows --

def zero_rows(root: Path):
    """{(loop, name, order)} of the coefficients whose exact row is the number 0 (structural zeros)"""
    out = set()
    for L in range(1, 6):
        man = json.loads((root / f'{L}L' / 'manifest.json').read_text())
        for m in man['masters']:
            f = (m.get('files') or {}).get('symbolic.json')
            if not f:
                continue
            for o, e in json.loads((root / f'{L}L' / f['path']).read_text()):
                if re.fullmatch(r'\s*[-+]?\s*0\s*', e):
                    out.add((L, m['name'], int(o)))
    return out


# ----------------------------------------------------------------------------------------- overlay exact forms --

RRC = re.compile(r'RRC\(([A-Za-z0-9_]+(?:\([0-9, ]+\))?)\)')
KNOWN = {'zeta', 'pi', 'log', 'polylog', 'sqrt', 'HI', 'HR', 'om1', 'cl2', 'cl4', 'Cl2', 'Cl4', 's6', 'D5bar1', 'D5bar2', 'D6bar1',
         'D6bar2', 'T111bar4', 'BNbar4', 'GPLReal', 'GPLImag', 'PolylogImag', 'PolylogReal', 'ImLi3', 'ImLi4', 'i', 'z2', 'z3', 'z4',
         'z5', 'z6', 'sqrt3'}


def classify_form(form, finite, where):
    """(class, X, constants) of an overlay exact form: RRC(X) = a finite constant of the released tree (digits from the tree)"""
    found = {}
    for name in RRC.findall(form):
        key = 'vakint::' + name.replace(' ', '')
        if key not in finite:
            raise SystemExit(f'{where}: RRC({name}) is not a finite constant of the data tree')
        found[key] = finite[key]
    rest = RRC.sub('', form).replace('𝑖', ' i ').replace('𝜋', ' pi ')
    unknown = sorted({w for w in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', rest)} - KNOWN)
    if unknown:
        raise SystemExit(f'{where}: unknown symbols {unknown} in the exact form (add them to KNOWN only if vakint evaluates them '
                         'at any precision)')
    if form.strip() in ('0', '-0'):
        return 'Z', None, {}
    return ('S*', min(found.values()), found) if found else ('S', None, {})


def load_overlay(path: Path, finite, tree_by_name):
    ov = json.loads(path.read_text())
    if ov.get('schema') != OVERLAY_SCHEMA:
        raise SystemExit(f'{path}: schema {ov.get("schema")!r} is not {OVERLAY_SCHEMA}')
    cells = {}
    for e in ov['entries']:
        where = f"overlay {e['loop']}L/{e['master']} ε^{e['order']}"
        if e.get('exact_form'):
            cls, X, consts = classify_form(e['exact_form'], finite, where)
            if e.get('status') and e['status'] != 'N' and {e['status'], cls} not in ({'S'}, {'S*'}, {'Z'}, {'S', 'Z'}):
                raise SystemExit(f"{where}: status {e['status']} but the form classifies as {cls}")
        elif e.get('status') == 'N':
            cls, X, consts = 'N', int(e['digits']), {}
        else:
            raise SystemExit(f'{where}: neither an exact form nor numerical digits')
        if e['kind'] == 'override':
            rel = tree_by_name.get((e['loop'], e['master']))
            if rel is None or rel['master_id'] != e['master_id']:
                raise SystemExit(f'{where}: not a master of the data tree with this master_id (overlay for another tree?)')
        cells[(e['loop'], e['master'], int(e['order']))] = {
            'class': cls, 'X': X, 'constants': sorted(consts), 'state': e['review_state'], 'agreement_digits': e.get('agreement_digits'),
            'claimed_digits': (e.get('released') or {}).get('digits'), 'exact_form': e.get('exact_form'), 'method': e.get('method'),
            'fit': e.get('fit'), 'provenance': e.get('provenance'), 'flags': e.get('flags') or [], 'also': e.get('also') or [],
            'kind': e['kind']}
    return ov, cells


# ---------------------------------------------------------------------------------------------------- constants --

BUILTIN_DEF = {
    'pi': 'π', 'EulerGamma': 'Euler–Mascheroni constant γ_E', 'sqrt3': '√3',
    'z2': 'ζ(2) = π²/6', 'z3': 'ζ(3)', 'z4': 'ζ(4) = π⁴/90', 'z5': 'ζ(5)', 'z6': 'ζ(6) = π⁶/945',
    'S2': 'S₂ = 4/(9√3) Cl₂(π/3)', 'cl2': 'Cl₂(π/3) = Im Li₂(e^{iπ/3})', 'cl4': 'Cl₄(π/3) = Im Li₄(e^{iπ/3})',
    's6': 's₆ = Σ_{m≥k≥1} (−1)^{m+k}/(m⁵ k) (alternating Euler sum, weight 6)',
    'D6': 'MATAD D6 (three-loop tetrahedron, six massive lines): D6bar(0) + 6 ζ(3)',
    'T1ep': 'ε¹ coefficient of the two-loop equal-mass sunset (MATAD T1)', 'T1ep2': 'ε² coefficient of the two-loop equal-mass sunset (MATAD T1)',
    'D5bar1': 'MATAD D5bar(1): ε¹ coefficient of the three-loop tetrahedron with five massive lines',
    'D5bar2': 'MATAD D5bar(2): ε² coefficient of the three-loop tetrahedron with five massive lines',
    'D6bar1': 'MATAD D6bar(1): ε¹ coefficient of the three-loop tetrahedron with six massive lines',
    'D6bar2': 'MATAD D6bar(2): ε² coefficient of the three-loop tetrahedron with six massive lines',
    'T111bar4': 'MATAD T111bar(4): ε⁴ coefficient of the two-loop sunset with three massive lines',
    'BNbar4': 'MATAD BNbar(4): ε⁴ coefficient of the three-loop equal-mass banana',
    'HR': 'Re H(w; e^{iπ/3}): harmonic polylogarithms, letters {−1, 0, 1}, weight ≤ 6',
    'HI': 'Im H(w; e^{iπ/3}): harmonic polylogarithms, letters {−1, 0, 1}, weight ≤ 6',
    'GPLReal': 'Re G(a₁,…,aₙ; 1): multiple polylogarithms, letters 0, ±1, sixth roots of unity, weight ≤ 6',
    'GPLImag': 'Im G(a₁,…,aₙ; 1): multiple polylogarithms, letters 0, ±1, sixth roots of unity, weight ≤ 6',
    'PolylogReal': 'Re Li_n(z), z = e^{iπ/3} (n ≤ 6) or i/√3', 'PolylogImag': 'Im Li_n(z), z = e^{iπ/3} (n ≤ 6) or i/√3',
    'zeta': 'ζ(n), n ≥ 2 (Riemann zeta values)', 'polylog': 'Li_n(q), q rational (classical polylogarithm)',
    'log': 'ln q, q rational', 'MZV': 'multiple zeta values ζ(n₁,…,n_k)',
}


def finite_definition(name):
    """standard definition of a finite-digit record, by its name pattern (REPORT: periods/l6 LIFTED-DEFINITIONS.md)"""
    m = re.fullmatch(r'K5Moment\((\d)\)', name)
    if m:
        return f'c₅,{m.group(1)} = ∫₀^∞ t^{m.group(1)} K₀(t)⁵ dt (Bessel moment, Broadhurst notation)'
    m = re.fullmatch(r'K5LogMoment\((\d)\)', name)
    if m:
        return f'3 ∫₀^∞ t^{m.group(1)} K₀(t)⁵ (ln t + γ_E − ln 2) dt (log-weighted Bessel moment)'
    m = re.fullmatch(r'K5MomentTaylor\((\d),(\d)\)', name)
    if m:
        return (f'[ε^{m.group(2)}] e^{{3(γ_E−ln2)ε}} ∫₀^∞ t^{{{m.group(1)}+3ε}} K_ε(t)⁵ dt (ε-jet of the D = 2−2ε four-loop banana)')
    m = re.fullmatch(r'PR(\d+)(p|d)?D2(?:Eps(\d))?', name)
    if m:
        top = f"FMFT four-loop topology PR{m.group(1)}{m.group(2) or ''}" + (' (one dotted line)' if m.group(2) else '')
        j = m.group(3)
        return (f'{top} in D = 2, unit masses, Euclidean' if j is None else
                f'[δ^{j}] e^{{4γ_E δ}} × {top} in D = 2−2δ, unit masses, Euclidean')
    m = re.fullmatch(r'PR14ep(\d)', name)
    if m:
        return f'[ε^{m.group(1)}] of the FMFT four-loop master PR14 in D = 4−2ε (unit masses)'
    m = re.fullmatch(r'(T111bar|BNbar|D5bar)(\d)', name)
    if m:
        what = {'T111bar': 'two-loop sunset with three massive lines', 'BNbar': 'three-loop equal-mass banana',
                'D5bar': 'three-loop tetrahedron with five massive lines'}[m.group(1)]
        return f'MATAD {m.group(1)}({m.group(2)}): ε^{m.group(2)} coefficient of the {what}'
    if name == 's6weak':
        return 'Σ_{m≥k≥1} (−1)^{m+k}/(m⁵ k) (alternating Euler sum; PARI convention)'
    return None


TOKEN = re.compile(r'vakint::([A-Za-z0-9_]+)(\([^()]*(?:\([^()]*\))?[^()]*\))?|symbolica::(zeta|polylog)\(([^()]*)\)|(?<![A-Za-z_:])log\(([^()]*)\)|𝜋')


def constant_usage(root: Path, B, finite):
    """{key: number of coefficients of the tree that contain it}"""
    use = {}
    for L in range(1, 6):
        man = json.loads((root / f'{L}L' / 'manifest.json').read_text())
        for m in man['masters']:
            f = (m.get('files') or {}).get('symbolic.json')
            if not f:
                continue
            for o, e in json.loads((root / f'{L}L' / f['path']).read_text()):
                keys = set()
                for t in TOKEN.finditer(e):
                    if t.group(0) == '𝜋':
                        keys.add('pi')
                    elif t.group(1):
                        name, args = t.group(1), (t.group(2) or '').replace(' ', '')
                        if f'vakint::{name}{args}' in finite:
                            keys.add(f'vakint::{name}{args}')
                        elif name in B['projections']:
                            keys.add(name)
                        else:
                            keys.add(name)
                    elif t.group(3) == 'zeta':
                        keys.add(f'zeta({t.group(4)})')
                    elif t.group(3) == 'polylog':
                        keys.add(f'polylog({t.group(4).replace(" ", "")})')
                    else:
                        keys.add(f'log({t.group(5)})')
                for k in keys:
                    use.setdefault(k, {}).setdefault(f'{L}L', 0)
                    use[k][f'{L}L'] += 1
    return use


def constants_table(root, B, finite, names, overlay_cells):
    use = constant_usage(root, B, finite)
    rows = []
    for key, by_loop in use.items():
        if key == 'om1' or key == 'ep':
            continue
        if key.startswith('vakint::'):
            name = key[len('vakint::'):]
            d = finite_definition(name)
            kind, digits, arb = 'finite', finite[key], False
        else:
            name = key
            base = re.sub(r'\(.*', '', key)
            d = BUILTIN_DEF.get(key) or BUILTIN_DEF.get(base)
            if base == 'zeta':
                d = f'ζ({key[5:-1]})'
            elif base == 'polylog':
                n, q = key[8:-1].split(',')
                d = f'Li_{n}({q})'
            elif base == 'log':
                d = f'ln {key[4:-1]}'
            kind, digits, arb = 'built-in', None, True
        if d is None:
            raise SystemExit(f'constant {key}: no definition known to overview.py (add one)')
        nm = (names or {}).get(name) or {}
        rows.append({'code': name, 'display': nm.get('display'), 'kind': kind, 'definition': d, 'digits': digits,
                     'arbitrary_precision': arb, 'used_in': dict(sorted(by_loop.items())), 'used_total': sum(by_loop.values()),
                     'traditional': nm.get('traditional'), 'traditional_terms': nm.get('terms'),
                     'traditional_status': nm.get('status'), 'traditional_check': nm.get('check')})
    # constants that only overlay (pending) forms use
    pend = {}
    for (L, n, o), c in overlay_cells.items():
        for k in c['constants']:
            pend[k] = pend.get(k, 0) + 1
    for r in rows:
        r['pending_uses'] = pend.get('vakint::' + r['code'], 0) if r['kind'] == 'finite' else 0

    def order(r):
        group = 0 if r['kind'] == 'built-in' else 1
        return (group, -(r['digits'] or 10 ** 9) if group else 0, r['code'].lower())
    rows.sort(key=order)
    return rows


# ------------------------------------------------------------------------------------------------------- build --

def build(a):
    root = a.data_root.resolve()
    T = json.loads(a.tables.read_text())
    if T['sources']['data_root_manifest_sha256'] != sha256_file(root / 'manifest.json'):
        raise SystemExit(f'{a.tables}: made from another data tree than {root} (root manifest digest differs)')
    if pt.content_hash(T) != T.get('content_sha256'):
        raise SystemExit(f'{a.tables}: content differs from its own content_sha256')
    finite = T['rule']['finite_constants']
    B = pt.read_builtins(a.vakint_src / 'inspection' / 'constants.rs')
    B['route_sources'] = pt.read_route_sources(a.vakint_src)
    consts, labels, masters = pt.load_tree(root, B)
    by_id = {m['master_id']: m for m in masters}
    pt.tier_check(root, by_id)                       # sets name_640 on the 5L masters
    tree_by_name = {(m['loop'], m['name']): m for m in masters}
    zeros = zero_rows(root)
    names = json.loads(a.names.read_text()) if a.names else None
    if names is not None and 'names' in names:          # vakint.constant-names.v1 document
        names = names['names']

    ov, ocells = (load_overlay(a.overlay, finite, tree_by_name) if a.overlay else (None, {}))
    mode = 'preview' if ov else 'released'
    basis = (ov or {}).get('basis')

    # ---------------------------------------------------------------- per-master table (b), released cells
    def rel_cell(m, o):
        c = m['orders'].get(o)
        if c is None:
            return None
        cls = 'Z' if (m['loop'], m['name'], o) in zeros else c['class']
        if cls == 'Z' and c['class'] != 'S':
            raise SystemExit(f"{m['loop']}L/{m['name']} ε^{o}: zero row classified {c['class']}")
        return {'class': cls, 'X': c.get('X') if cls in ('S*', 'N') else None}

    t1_rel = {(r['loop'], r['name']): r for r in T['table1']}
    orders1 = T['orders_table1']

    def row_of(m, pending_only=False):
        cells = {}
        base = t1_rel.get((m['loop'], m['name']))
        for o in orders1:
            rc = None if pending_only else rel_cell(m, o)
            if rc is not None and base is not None:
                bc = base['cells'][str(o)]
                if bc is None or bc['class'] != (rc['class'] if rc['class'] != 'Z' else 'S') or bc.get('X') != rc.get('X'):
                    raise SystemExit(f"{m['loop']}L/{m['name']} ε^{o}: tables.json says {bc}, the tree {rc}")
                rc['marks'] = bc.get('marks') or []
            pc = ocells.get((m['loop'], m['name'], o))
            if pc is not None:
                e = dict(rc) if rc else {'class': None, 'X': None, 'marks': []}
                e['pending'] = {'class': pc['class'], 'X': pc['X'], 'state': pc['state'], 'agreement_digits': pc['agreement_digits'],
                                'claimed_digits': pc['claimed_digits'], 'kind': pc['kind']}
                if rc is not None and rank(pc) > rank(rc):
                    e['pending']['worse_than_released'] = True
                cells[str(o)] = e
            else:
                cells[str(o)] = rc
        return cells

    rows = []
    loops_basis = basis['loop'] if basis else None
    dropped = {d['name'] for d in basis['dropped']} if basis else set()
    for m in masters:
        if m['loop'] == loops_basis and m['name'] in dropped:
            continue
        if basis and m['loop'] == loops_basis and basis['kept'].get(m['name']) != m['master_id']:
            raise SystemExit(f"overlay basis: {m['loop']}L/{m['name']} is neither dropped nor kept with this master_id")
        rows.append({'loop': m['loop'], 'name': m['name'], 'master_id': m['master_id'], 'name_640': m.get('name_640'),
                     'cells': row_of(m)})
    if basis:
        for x in basis['added']:
            if (loops_basis, x['name']) in tree_by_name:
                raise SystemExit(f"overlay basis: added master {x['name']} already in the tree")
            fake = {'loop': loops_basis, 'name': x['name'], 'orders': {}}
            rows.append({'loop': loops_basis, 'name': x['name'], 'master_id': x['master_id'], 'design_name': x['design_name'],
                         'replaces': x.get('replaces'), 'pending_master': True, 'cells': row_of(fake, pending_only=True)})
        rows.sort(key=lambda r: (r['loop'], int(re.sub(r'\D', '', r['name']) or 0)))

    def eff(c):
        """the class a cell has with the overlay applied (pending if it comes from the overlay)"""
        if c is None:
            return None, False
        p = c.get('pending')
        if p and (c['class'] is None or rank(p) <= rank(c)):
            return {'class': p['class'], 'X': p['X']}, True
        return ({'class': c['class'], 'X': c['X']} if c['class'] else None), False

    # ---------------------------------------------------------------- (a) summary per loop and order
    def summarize(rws, use_overlay):
        out = {}
        for o in orders1:
            acc = {'S': 0, 'S*': 0, 'N': 0, 'Z': 0, 'pending': {'S': 0, 'S*': 0, 'N': 0, 'Z': 0}, 'sx': [], 'nx': [], 'held': 0}
            for r in rws:
                c = r['cells'][str(o)]
                if use_overlay:
                    e, pend = eff(c)
                else:
                    e, pend = (({'class': c['class'], 'X': c['X']} if c and c['class'] else None), False)
                if e is None:
                    continue
                acc['held'] += 1
                (acc['pending'] if pend else acc)[e['class']] += 1
                if e['class'] == 'S*':
                    acc['sx'].append(e['X'])
                if e['class'] == 'N':
                    acc['nx'].append(e['X'])
            if acc['held'] == 0:
                continue
            out[str(o)] = {'held': acc['held'], 'S': acc['S'], 'S*': acc['S*'], 'N': acc['N'], 'Z': acc['Z'], 'pending': acc['pending'],
                           'sstar_digits': [min(acc['sx']), max(acc['sx'])] if acc['sx'] else None,
                           'n_digits_min': min(acc['nx']) if acc['nx'] else None,
                           'n_digits_median': statistics.median(acc['nx']) if acc['nx'] else None}
        return out

    summary = []
    for L in range(1, 6):
        tree_rows = [{'cells': row_of(m)} for m in masters if m['loop'] == L]
        s = summarize([{'cells': {k: (None if v is None or v.get('class') is None else v) for k, v in r['cells'].items()}}
                       for r in tree_rows], False)
        summary.append({'loop': L, 'kind': 'released', 'label': f'{L}L', 'masters': len(tree_rows), 'orders': s,
                        'basis': (T['pole_depth']['five_loop']['tag'] if L == 5 else None)})
        if ov and L == loops_basis:
            br = [r for r in rows if r['loop'] == L]
            summary.append({'loop': L, 'kind': 'preview', 'label': f"{L}L {basis['name']}", 'masters': len(br),
                            'orders': summarize(br, True), 'basis': basis['name']})
        elif ov and any(L == k[0] for k in ocells):
            summary.append({'loop': L, 'kind': 'preview', 'label': f'{L}L + pending', 'masters': len(tree_rows),
                            'orders': summarize([r for r in rows if r['loop'] == L], True), 'basis': None})

    # ---------------------------------------------------------------- (c) result of any L-loop input
    orders2 = T['orders_table2']
    res_rows = []
    for r in T['table2']:
        cells = {}
        for o in orders2:
            c = r['cells'][str(o)]
            if c.get('na') == 'below':
                cells[str(o)] = {'state': 'below'}
            elif c.get('na') == 'beyond':
                cells[str(o)] = {'state': 'none'}
            elif c.get('na') == 'partial':
                cells[str(o)] = {'state': 'partial', 'returned': c.get('returned'), 'marks': c.get('marks') or []}
            else:
                cells[str(o)] = {'state': 'ok', 'class': c['class'], 'X': c.get('X'), 'marks': c.get('marks') or [],
                                 'others': c.get('others')}
        res_rows.append({'label': f"{r['loop']}L", 'loop': r['loop'], 'kind': 'released', 'cells': cells, 'marks': r.get('marks') or []})

    former = []
    if basis:
        # the overlay's basis (p = 0: pole-free basis) + the lower-loop product factors of the precision tables' route model
        alias_map = pt.aliases(masters)
        depth, depth_doc, five, factors, measured = pt.load_depths(root, a.pole_depth, masters, alias_map)
        br = [r for r in rows if r['loop'] == loops_basis]
        cells = {}
        k_max = max(int(o) for r in br for o, c in r['cells'].items() if eff(c)[0] is not None)
        for o in orders2:
            if o < -loops_basis:
                cells[str(o)] = {'state': 'below'}
                continue
            if o > k_max:
                cells[str(o)] = {'state': 'none'}
                continue
            ent, missing = [], []
            for r in br:
                hold = [int(j) for j, c in r['cells'].items() if eff(c)[0] is not None]
                lo = min(hold)
                for j in range(lo, o + 1):
                    e, pend = eff(r['cells'].get(str(j)))
                    if e is None:
                        missing.append(r['name'])
                        break
                    ent.append((e, pend))
            fmiss = []
            for f in factors:
                for rid in alias_map.get(f['m']['master_id'], [f['m']['master_id']]):
                    n = pt.needs(by_id[rid], o + f['offset'])
                    if n is None:
                        fmiss.append(f"{f['m']['loop']}L/{f['m']['name']}")
                        break
                    ent += [({'class': 'Z' if (mm['loop'], mm['name'], j) in zeros else cc['class'], 'X': cc.get('X')}, False)
                            for mm, j, cc in n]
            w = max((e for e, p in ent), key=rank) if ent else None
            wp = any(p for e, p in ent if rank(e) == rank(w))
            if missing or fmiss:
                cells[str(o)] = {'state': 'partial', 'returned': label(w), 'missing': sorted(set(missing)) + fmiss, 'marks': []}
            else:
                cells[str(o)] = {'state': 'ok', 'class': 'S' if w['class'] == 'Z' else w['class'], 'X': w.get('X'), 'marks': [],
                                 'pending': True, 'set_by_pending': wp}
        res_rows.append({'label': f"{loops_basis}L {basis['name']}", 'loop': loops_basis, 'kind': 'preview', 'cells': cells,
                         'marks': [], 'note': f"{basis['name']}: {basis['status']}; spurious-pole order p = 0 for its masters (pole-free "
                                              'basis), lower-loop product factors as in the released row'})
        # former masters: exact where every numerical coefficient entering the order is exact, else fail closed
        eff_class = {}
        for r in br:
            for j, c in r['cells'].items():
                e, p = eff(c)
                rel = None if c is None or c.get('class') is None else {'class': c['class'], 'X': c['X']}
                eff_class[(r['name'], int(j))] = (rel, e)
        for d in basis['dropped']:
            cells = {}
            for o in orders2:
                if o < -loops_basis:
                    cells[str(o)] = {'state': 'below'}
                    continue
                x = d['orders'].get(str(o))
                if x is None:
                    cells[str(o)] = {'state': 'none'}
                    continue
                if x['contributors'] is None:
                    cells[str(o)] = {'state': 'numerical', 'v02x': f"N({int(x['v02x_digits'])})"}
                    continue
                open_rel = [f'{n}|{j}' for n, j in x['contributors']
                            if (eff_class.get((n, j)) or (None, None))[0] is None or eff_class[(n, j)][0]['class'] in ('N',)]
                open_pend = [f'{n}|{j}' for n, j in x['contributors']
                             if (eff_class.get((n, j)) or (None, None))[1] is None or eff_class[(n, j)][1]['class'] in ('N',)]
                sstar = [eff_class[(n, j)][1]['X'] for n, j in x['contributors']
                         if (eff_class.get((n, j)) or (None, None))[1] is not None and eff_class[(n, j)][1]['class'] == 'S*']
                cells[str(o)] = {'state': 'fails-closed' if open_pend else 'exact', 'open_released': len(open_rel),
                                 'open_pending': len(open_pend), 'n_contributors': len(x['contributors']),
                                 'open_list': open_pend[:12], 'v02x': x['v02x'], 'sstar_X': min(sstar) if sstar else None}
            former.append({'label': f"{d['name']} (former)", 'name': d['name'], 'replaced_by': d.get('replaced_by'), 'cells': cells})

    # ---------------------------------------------------------------- (d) constants, (e) method
    ctab = constants_table(root, B, finite, names, ocells)
    method = method_paragraphs(T, ov, ocells)

    pend_list = []
    for (L, n, o), c in sorted(ocells.items(), key=lambda kv: (kv[0][0], int(re.sub(r'\D', '', kv[0][1])), kv[0][2])):
        pend_list.append({'loop': L, 'master': n, 'order': o, 'label': label(c), 'state': c['state'], 'kind': c['kind'],
                          'agreement_digits': c['agreement_digits'], 'claimed_digits': c['claimed_digits'], 'method': c['method'],
                          'fit': c['fit'], 'exact_form': c['exact_form'], 'provenance': c['provenance'], 'flags': c['flags'],
                          'also': c['also'],
                          'below_claim': (c['agreement_digits'] is not None and c['claimed_digits'] is not None
                                          and c['agreement_digits'] < c['claimed_digits'])})

    src = T['sources']
    out = {
        'schema': SCHEMA, 'mode': mode,
        'tree': {'label': T.get('label'), 'data_root': str(root), 'version': src['version'], 'git': T.get('data_git'),
                 'root_manifest_sha256': src['data_root_manifest_sha256'], 'tables_content_sha256': T['content_sha256'],
                 'tables_status': T['status'], 'five_loop_basis': T['pole_depth']['five_loop']['tag'],
                 'five_loop_build': src['five_loop_values']},
        'overlay': None if not ov else {
            'path': str(a.overlay.resolve()), 'sha256': sha256_file(a.overlay), 'title': ov.get('title'), 'created_utc': ov.get('created_utc'),
            'basis': None if not basis else {k: basis[k] for k in ('loop', 'name', 'status', 'rule')} | {
                'added': basis['added'], 'dropped': [d['name'] for d in basis['dropped']]},
            'states': {s: sum(1 for c in ocells.values() if c['state'] == s) for s in sorted({c['state'] for c in ocells.values()})},
            'sources': ov.get('sources')},
        'orders1': orders1, 'orders2': orders2,
        'summary': summary, 'table1': rows,
        'table1_counts': table1_counts(rows, eff),
        'result': res_rows, 'former': former,
        'constants': ctab, 'names': None if not a.names else {'path': str(a.names.resolve()), 'sha256': sha256_file(a.names)},
        'method': method, 'pending': pend_list,
        'generated_utc': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'generator': {'path': 'crates/vakint/docs/benchmark/tools/overview.py', 'sha256': sha256_file(Path(__file__))},
    }
    return out


def table1_counts(rows, eff):
    out = {}
    for r in rows:
        d = out.setdefault(f"{r['loop']}L", {'masters': 0, 'S': 0, 'S*': 0, 'N': 0, 'Z': 0, 'pending': 0, 'pending_masters': 0})
        d['masters'] += 1
        d['pending_masters'] += bool(r.get('pending_master'))
        for c in r['cells'].values():
            e, p = eff(c)
            if e is None:
                continue
            if p:
                d['pending'] += 1
            else:
                d[e['class']] += 1
    return out


def method_paragraphs(T, ov, ocells):
    q = T['pslq']
    ids = T.get('identities') or {}
    fc = T['rule']['finite_constants']
    paras = [
        ('Closed forms', '1L-3L and the poles of 4L: closed forms from the literature and the MATAD/FMFT expansions (harmonic '
         'polylogarithms at the sixth root of unity, the MATAD constants), carried exactly by vakint. 4L finite parts: exact in '
         f'{len(fc)} lifted four-loop constants (Bessel moments, two-dimensional ε-jets) computed to {min(fc.values())}-{max(fc.values())} '
         'digits (S*). Some S rows rest on numerical identities of the constant basis checked to thousands of digits '
         f"({len(ids)} identities; PRECISION-TABLES.md)."),
        ('Five loops', f"Structural zeros from the UV forest bound and exact rows from pole cancellation; {q['pslq_orders']} of the "
         f"{q['exact_orders']} exact 5L orders are PSLQ identifications ({q['pslq_sstar_orders']} S*) and {q['measured_zero_orders']} "
         'measured zeros, each passing the 40-excess-digit bar (numerical identifications, not proofs). The other 5L orders are '
         'numerical (high-precision grids and dispersive gluing; digits = validated significant digits).'),
        ('R operation (preview)', 'Pole coefficient K G = (A) − Δ: part (A) = −K Σ_S Δ(S) G/S from the subdivergences (ordinary massive R '
         'operation, exact, lower-loop coefficients of vakint); Δ = overall counterterm, a pure number: absent for ω < 0 (complete '
         'result), otherwise fitted as an integer relation of [numeric − (A)] in the basis 1, ζ3, π⁴, ζ5, π⁶, ζ3², ζ7, ζ3π⁴ '
         '(accepted at a margin ≥ 20 digits; to be confirmed by R*).'),
        ('Numerical check', 'Every symbolic form is compared with the numerical value it replaces, to that value’s full claimed digits '
         '(user ruling 2026-11-10), at adoption and again when the release tooling integrates it; a mismatch blocks the substitution. '
         'Pending entries list their digits of agreement; where they fall short of the claimed digits the entry is flagged.'),
    ]
    if ov:
        n = {s: sum(1 for c in ocells.values() if c['state'] == s) for s in {c['state'] for c in ocells.values()}}
        paras.append(('Pending (this PDF)', '; '.join(f'{v} {k}' for k, v in sorted(n.items())) + '. None of them is in the released '
                      'data tree; they are counted apart and never as released.'))
    return [{'title': t, 'text': x} for t, x in paras]


# ------------------------------------------------------------------------------------------------------ markdown --

def cell_md(c, mode):
    if c is None:
        return '—'
    p = c.get('pending')
    base = '—' if c.get('class') is None else label(c)
    if p and mode == 'preview' and not p.get('worse_than_released'):
        return f"{base} → **{label(p)}†**" if c.get('class') else f"**{label(p)}†**"
    return base


def markdown(O):
    m = O['mode']
    out = [f"# What vakint computes today (1L-5L single-scale vacuum): {m.upper()}", '',
           f"Data tree {O['tree']['label'] or O['tree']['data_root']} ({O['tree']['version']}, 5L basis {O['tree']['five_loop_basis']}) · "
           f"tables `{O['tree']['tables_content_sha256'][:16]}` · generated {O['generated_utc']}", '']
    if O['overlay']:
        ov = O['overlay']
        out += [f"**PREVIEW**: + overlay `{ov['path']}` (`{ov['sha256'][:16]}`, {ov['created_utc']}): {ov['title']}. Pending entries "
                f"({', '.join(f'{v} {k}' for k, v in ov['states'].items())}) are marked † and counted apart; none is released.", '']
    out += ['S = exact (periods at any precision) · S*(X) = exact, contains periods known to X digits · N(X) = numerical only, X digits '
            '· 0 = structural zero · — = not held / not available.', '', '## (a) Summary per loop and Laurent order', '',
            '| row | order | held | S | S* | S* digits | N | N digits min / median | 0 | pending (S/S*/N/0) |',
            '|---|---|---|---|---|---|---|---|---|---|']
    for s in O['summary']:
        for o, x in s['orders'].items():
            p = x['pending']
            out.append(f"| {s['label']}{' (' + str(s['masters']) + ' masters)' if o == next(iter(s['orders'])) else ''} | ε^{o} | {x['held']} | "
                       f"{x['S']} | {x['S*']} | {'-'.join(map(str, x['sstar_digits'])) if x['sstar_digits'] else ''} | {x['N']} | "
                       f"{'' if x['n_digits_min'] is None else str(x['n_digits_min']) + ' / ' + str(x['n_digits_median'])} | {x['Z']} | "
                       f"{'' if s['kind'] == 'released' else '/'.join(str(p[k]) for k in ('S', 'S*', 'N', 'Z'))} |")
    o1 = O['orders1']
    out += ['', '## (b) Master coefficients', '', '| master | ' + ' | '.join(f'ε^{o}' for o in o1) + ' |', '|---|' + '---|' * len(o1)]
    last = None
    for r in O['table1']:
        if r['loop'] != last:
            c = O['table1_counts'][f"{r['loop']}L"]
            out.append(f"| **{r['loop']}L** ({c['masters']} masters: {c['S']} S, {c['S*']} S*, {c['N']} N, {c['Z']} zero"
                       + (f"; {c['pending']} pending†" if c['pending'] else '') + ') |' + ' |' * len(o1))
            last = r['loop']
        nm = r['name'] + (f" ({r['name_640']})" if r.get('name_640') else '') + (f" {r['design_name']} (new, replaces {r['replaces']})†"
                                                                                  if r.get('design_name') else '')
        out.append(f"| {r['loop']}L/{nm} | " + ' | '.join(cell_md(r['cells'][str(o)], m) for o in o1) + ' |')
    o2 = O['orders2']
    out += ['', '## (c) Result of any single-scale vacuum input at L loops', '', '| L | ' + ' | '.join(f'ε^{o}' for o in o2) + ' |',
            '|---|' + '---|' * len(o2)]

    def rc(c):
        s = c['state']
        if s == 'below':
            return ' '
        if s == 'none':
            return 'not available'
        if s == 'partial':
            return f"not for every input (others {c['returned']})"
        if s == 'numerical':
            return f"numerical (v0.2.x {c['v02x']})"
        if s == 'exact':
            return 'exact' + (f" (S*({c['sstar_X']}))" if c.get('sstar_X') else '') + '†'
        if s == 'fails-closed':
            return f"fails closed ({c['open_pending']} of {c['n_contributors']} coefficients numerical; released data {c['open_released']})"
        t = {'S': 'exact', 'S*': f"exact, X = {c['X']} d", 'N': f"numerical, ≥ {c['X']} d"}[c['class']]
        return t + ('†' if c.get('pending') else '')
    for r in O['result'] + O['former']:
        out.append(f"| {r['label']} | " + ' | '.join(rc(r['cells'][str(o)]) for o in o2) + ' |')
    out += ['', '## (d) Named constants', '', '| constant | definition | digits | any precision | coefficients (released) | traditional form |',
            '|---|---|---|---|---|---|']
    for c in O['constants']:
        out.append(f"| {c['display'] or c['code']}{' (' + c['code'] + ')' if c['display'] else ''} | {c['definition']} | "
                   f"{'any' if c['arbitrary_precision'] else c['digits']} | {'yes' if c['arbitrary_precision'] else 'no'} | "
                   f"{c['used_total']} ({', '.join(f'{k} {v}' for k, v in c['used_in'].items())})"
                   + (f" + {c['pending_uses']} pending" if c.get('pending_uses') else '') + ' | '
                   + (f"{c['traditional']}" + (f" ({c['traditional_status']})" if c.get('traditional_status') else '') if c.get('traditional') else '')
                   + ' |')
    out += ['', '## (e) How the symbolic forms were obtained and checked', ''] + [f"- **{p['title']}.** {p['text']}" for p in O['method']]
    if O['pending']:
        out += ['', '## Pending entries (preview overlay)', '', '| coefficient | class | review state | method | agreement / claimed digits | flags |',
                '|---|---|---|---|---|---|']
        for p in O['pending']:
            out.append(f"| {p['loop']}L/{p['master']} ε^{p['order']} | {p['label']} | {p['state']} | {p['method'] or p['kind']}"
                       f"{' · ' + p['fit'] if p.get('fit') else ''} | "
                       f"{p['agreement_digits'] if p['agreement_digits'] is not None else '—'} / {p['claimed_digits'] if p['claimed_digits'] is not None else '—'} | "
                       + '; '.join((['agreement below the claimed digits: check to full digits before adoption'] if p['below_claim'] else [])
                                   + p['flags']) + ' |')
    return '\n'.join(out) + '\n'


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--tables', type=Path, required=True)
    ap.add_argument('--json', type=Path, required=True)
    ap.add_argument('--md', type=Path)
    ap.add_argument('--overlay', type=Path)
    ap.add_argument('--names', type=Path)
    ap.add_argument('--vakint-src', type=Path, default=pt.CRATE / 'src')
    ap.add_argument('--pole-depth', type=Path, default=pt.BENCH / 'data' / 'pole-depth.json')
    a = ap.parse_args(argv)
    O = build(a)
    a.json.write_text(json.dumps(O, indent=1, ensure_ascii=False) + '\n')
    if a.md:
        a.md.write_text(markdown(O))
    print(f"overview {O['mode']}: " + '; '.join(
        f"{s['label']} " + ' '.join(f"{o}:{x['S']}/{x['S*']}/{x['N']}/{x['Z']}" + (f"+{sum(x['pending'].values())}p" if s['kind'] != 'released' else '')
                                    for o, x in s['orders'].items()) for s in O['summary']))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
