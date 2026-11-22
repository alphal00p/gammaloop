#!/usr/bin/env python3
"""IRH-DEPTH mechanical variant matrix (MATRIX-R45): classify the native cells.

Input: RUNS/<profile>/<variant>/<mode>/{log,meta.txt,out/L3/integral-000000.json}
for mode in generic (GAMMALOOP_IR_HARD=off), child, strict
(GAMMALOOP_IR_HARD_DEPTH_STRICT=1), guard (GAMMALOOP_IR_HARD_CHILD_GUARD=1).

Each child-route cell is compared with the generic cell of the same profile and
variant:
  RIGHT    the child record equals generic (zb, status, longitudinal) and the
           IR-hard correction itself took the child route (an IR_HARD_CHILD line
           with outcome=child, no IR_HARD_DEPTH_REFUSED line);
  REFUSED  equal to generic, the correction refused loudly and counted
           (IR_HARD_DEPTH_REFUSED reason=<class>, IR_HARD_DEPTH_SUMMARY refused>=1);
  BOUNDARY both routes fail with beta's public-boundary refusal of a graph weight
           (status failed, "IRH-DEPTH regularity: <factor> ..."): loud, no value;
  LOUD     both routes fail with the same kind of error (no value; e.g. the input
           cannot be read);
  DECLINED equal to generic but the child route was neither taken nor refused by
           the certificate (another decline, e.g. the guard); not RIGHT;
  WRONG    a successful child record whose value differs from generic, or a
           child success where generic fails.
usage: classify.py RUNS EXPECTED.tsv [OUT.json]   (exit 1 on any WRONG or any cell
off its expected disposition)
"""
import json
import re
import sys
from pathlib import Path


def record(cell):
    path = cell / 'out/L3/integral-000000.json'
    if not path.exists():
        return None
    rec = json.loads(path.read_bytes())
    return {k: rec.get(k) for k in ('zb', 'status', 'longitudinal', 'zero_reason', 'locality_error')}


def meta_rc(cell):
    text = (cell / 'meta.txt').read_text() if (cell / 'meta.txt').exists() else ''
    m = re.search(r'rc=(\d+)', text)
    return int(m[1]) if m else None


def log_facts(cell):
    log = (cell / 'log').read_text(errors='replace') if (cell / 'log').exists() else ''
    outcomes = re.findall(r'^IR_HARD_CHILD graph=\S+ corr=\d+ hard_loops=\d+ omega=\d+ outcome=(\S+)', log, re.M)
    refused = re.findall(r'IR_HARD_DEPTH_REFUSED .*?reason=(\S+)', log)
    regularity = re.findall(r'^IR_HARD_DEPTH_REGULARITY .*', log, re.M)
    summary = [dict(kv.split('=', 1) for kv in line.split()[1:] if '=' in kv)
               for line in re.findall(r'^IR_HARD_DEPTH_SUMMARY .*', log, re.M)]
    boundary = re.findall(r'IRH-DEPTH regularity: (\S+) (?:has a pole of order (\d+)|the ε order)', log)
    errors = [l for l in log.splitlines() if l.startswith(('Error', 'error', 'thread ')) or 'panicked' in l][:3]
    # U7b: K3 records built on the uncontracted path (on or check), and the
    # checks that compared equal
    k3 = len(re.findall(r'^UNCONTRACTED kind=ir-hard \S+ id=\d+ outcome=uncontracted', log, re.M))
    k3_equal = len(re.findall(r'^UNCONTRACTED_CHECK kind=ir-hard .*equal=true', log, re.M))
    return dict(outcomes=outcomes, refused=refused, regularity=regularity[:2], summary=summary,
                boundary=boundary, errors=errors, k3=k3, k3_equal=k3_equal)


def pole(zb):
    if not zb:
        return None
    m = re.findall(r'([+-]?[0-9/]*)\*?beta::\{\}::cA\^3\*gammalooprs::\{\}::ε\^\(-1\)', zb.replace(' ', ''))
    return m[0] if m else zb[:80]


def classify(gen, child):
    g, c = gen['record'], child['record']
    gok = g is not None and g['status'] == 'computed'
    cok = c is not None and c['status'] == 'computed'
    f = child['facts']
    if cok and not gok:
        return 'WRONG', 'child succeeds where generic fails'
    if cok and gok:
        if c != g:
            return 'WRONG', f"child {c['zb']} != generic {g['zb']}"
        if f['refused']:
            return 'REFUSED', ','.join(sorted(set(f['refused'])))
        if 'child' in f['outcomes']:
            return 'RIGHT', 'outcome=child'
        return 'DECLINED', ','.join(sorted(set(f['outcomes']))) or 'no IR_HARD_CHILD line'
    # neither computed
    gb, cb = gen['facts']['boundary'], f['boundary']
    if gb and cb:
        return 'BOUNDARY', ','.join(sorted({b[0] for b in cb}))
    if (g is None) == (c is None) and (g is None or g['status'] == c['status']):
        return 'LOUD', (f['errors'][:1] or [str(c and c['status'])])[0][:160]
    return 'WRONG?', f'generic {g} child {c}'


def main():
    runs, expected_path = Path(sys.argv[1]), Path(sys.argv[2])
    expected = {}
    for line in expected_path.read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        variant, profiles, disposition, detail = (line.split('\t') + [''] * 4)[:4]
        expected[(variant, profiles)] = (disposition, detail)
    table, bad = [], []
    for profile in sorted(p.name for p in runs.iterdir() if p.is_dir()):
        for variant in sorted(v.name for v in (runs / profile).iterdir() if v.is_dir()):
            base = runs / profile / variant
            cells = {m: dict(record=record(base / m), facts=log_facts(base / m), rc=meta_rc(base / m))
                     for m in ('generic', 'child', 'strict', 'guard') if (base / m).exists()}
            if 'generic' not in cells:
                continue
            for mode in ('child', 'strict', 'guard'):
                if mode not in cells:
                    continue
                disposition, detail = classify(cells['generic'], cells[mode])
                key = next((k for k in [(variant, profile), (variant, '*')] if k in expected), None)
                want = expected.get(key)
                ok = want is not None and want[0] == disposition and (not want[1] or want[1] == detail or disposition in ('LOUD', 'RIGHT'))
                row = dict(profile=profile, variant=variant, mode=mode, disposition=disposition, detail=detail,
                           child=pole((cells[mode]['record'] or {}).get('zb')),
                           generic=pole((cells['generic']['record'] or {}).get('zb')),
                           child_status=(cells[mode]['record'] or {}).get('status'),
                           generic_status=(cells['generic']['record'] or {}).get('status'),
                           rc=cells[mode]['rc'], generic_rc=cells['generic']['rc'],
                           k3=cells[mode]['facts']['k3'], k3_equal=cells[mode]['facts']['k3_equal'],
                           summary=cells[mode]['facts']['summary'],
                           expected=want, as_expected=ok)
                table.append(row)
                if disposition.startswith('WRONG') or not ok:
                    bad.append(row)
    counts = {}
    for row in table:
        counts[row['disposition']] = counts.get(row['disposition'], 0) + 1
    out = dict(counts=counts, cells=len(table), off_expected=len(bad), table=table)
    if len(sys.argv) > 3:
        Path(sys.argv[3]).write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    for row in table:
        print('\t'.join(str(row[k]) for k in ('profile', 'variant', 'mode', 'disposition', 'detail', 'child', 'generic', 'k3', 'k3_equal', 'as_expected')))
    print('COUNTS', json.dumps(counts), 'cells', len(table), 'off_expected', len(bad))
    for row in bad:
        print('BAD', json.dumps(row, ensure_ascii=False)[:1500])
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
