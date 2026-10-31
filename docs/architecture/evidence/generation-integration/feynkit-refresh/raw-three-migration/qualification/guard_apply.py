"""Verify approved immutable inputs and restrict the new code to three tokens."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib

W = Path('/shared/lcnbr/dev/gammaloop/generation-feynkit-refresh')
O = Path(__file__).resolve().parent
F = W.parent / 'generation-oracle-numerator-families'
P = F / 'target/raw-three-migration-preparation-20261031-r2'
BASE = '5ea4c93b79dd33094072d7717e912aa4b05d6689'
PATH = 'crates/feynkit-generator/tests/generation_oracle/main.rs'
ENV = {k: v for k, v in os.environ.items() if k not in ('SYMBOLICA_LICENSE', 'SYMBOLICA_LICENSE_PATH', 'SYMBOLICA_LICENSE_SIGNED')}
ENV.update(GIT_OPTIONAL_LOCKS='0', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL='/dev/null')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def blob(b):
    return hashlib.sha1(b'blob %d\0' % len(b) + b).hexdigest()

def write(name, value):
    (O / name).write_text(json.dumps(value, indent=2) + '\n')

def git(*args):
    return subprocess.check_output(['git', '--git-dir=' + str(W.parent / 'main/.git'), *args], env=ENV)

def scopes(registry, reader):
    out = []
    for sec in tomllib.loads(registry.decode())['section']:
        for n in sec.get('note', []):
            for s in n.get('scope', []):
                b = reader(s['path'])
                t = b.decode()
                a, z = s.get('anchor'), s.get('end_anchor')
                assert not a or t.count(a) == 1
                assert not z or t.count(z) == 1
                r = t[t.index(a):t.index(z)].encode() if z else b
                assert blob(r) == s['digest'], (n['id'], s)
                out.append({'note': n['id'], **s, 'actual_blob': blob(r)})
    return out

if __name__ == '__main__':
    assert subprocess.check_output(['jj', '--ignore-working-copy', 'log', '--no-graph', '--no-pager', '-r', '@', '-T', 'commit_id'], cwd=W, env=ENV, text=True).strip() == BASE
    tree = git('ls-tree', '-rz', '--full-tree', BASE)
    paths = {}
    for row in tree.split(b'\0'):
        if not row:
            continue
        h, p = row.split(b'\t', 1)
        mode, kind, oid = h.decode().split()
        p = os.fsdecode(p)
        f = W / p
        data = os.fsencode(os.readlink(f)) if mode == '120000' else f.read_bytes()
        assert kind == 'blob' and blob(data) == oid, p
        paths[p] = {'blob': oid, 'mode': mode}
    assert re.search(r'^  enable = false;$', (W / 'nix-ci.nix').read_text(), re.M)
    specs = [
        (F / '.local/review-raw-three-migration-astra-20261031-004323Z-a3d128ed', 'c0692180c797034f923ae4bcf0a91712816a67bcf726d088787914c9a4d211e6', 'a7eb2b7d437e7c761681d3163bc7dad6bc3e67a71ee72be68be972274d7b2d59'),
        (F / '.local/review-raw-three-interface-astra-20261031-002103Z-8a5d8001', 'd9303d91d922c994f87b07e10cce895b44dba96e294d64be231a3ecdf921b14c', '0ce63d1e764a9c5a36d8639174738eecc3a0c921aa800ed7e089e355991318f8'),
        (W / '.local/review-sm-eeww-integration-astra-20261031-independent-r1', '5564bed2f642ba72a9208858d86292f8b9326fb1c55f24ef23a07d027a93624b', '740a537002c817054b2cf97315509298d144dd0c2c46d53cfbf332d939d48abc'),
    ]
    packets = []
    for root, review, seal in specs:
        assert sha(root / 'review.json') == review
        assert sha(root / 'SEAL.json') == seal
        packets.append(root)
    assert sha(P / 'handoff.json') == '067f4aa62546e10fdf35d8a99b40812ae4ac5c0e927daddcf2190ae353eb23ad'
    assert sha(P / 'SEAL.json') == '6a6ee263fa4d69d48023b896a607675af0ec88b3993ccfa3fd5ce8b17fddd6c2'
    assert sha(P / 'proposed-not-applied.patch') == '872261e0211f57960a1eb0544dbdb303934cd95e4c44b47c4e9cccd358d09fe6'
    packets += [P, F / 'target/raw-three-interface-proof-20261030-r1']
    sealed = []
    for root in packets:
        seal = json.loads((root / 'SEAL.json').read_text())
        assert sha(root / 'SHA256SUMS') == seal['manifest_sha256']
        count = 0
        for line in (root / 'SHA256SUMS').read_text().splitlines():
            digest, rel = line.split('  ', 1)
            assert not Path(rel).is_absolute() and '..' not in Path(rel).parts
            assert sha(root / rel) == digest, (root, rel)
            count += 1
        sealed.append({'root': str(root), 'seal_sha256': sha(root / 'SEAL.json'), 'files_reverified': count})
    old = (W / PATH).read_bytes()
    assert old == git('show', BASE + ':' + PATH) == git('show', '605178a65620b104805b012f1edc299710e8a176:' + PATH) == (P / 'before-main.rs.txt').read_bytes()
    assert blob(old) == 'e6af4ec644688e1333064c50e3b157ae631473ee'
    swaps = {'qed_ee4f_L0': ('7f2a34ae141a6f7c', '1d107b2170d73cd2'), 'xs_qed_eeee_L1': ('8eb2c3186607a1a2', '7734a2de85dd54c9'), 'xs_qcd_gguu_1pi_L1': ('75fbcf2f108240fb', '631d8c73b7db146b')}
    pat = rb'(\w+): Row \{[^\n]*?\}'
    rows = {m[1].decode(): m[0] for m in re.finditer(pat, old)}
    assert len(rows) == 149
    new = old
    for case, (a, b) in swaps.items():
        r = rows[case]
        a, b = ('numerators: 0x' + a).encode(), ('numerators: 0x' + b).encode()
        assert r.count(a) == 1
        new = new.replace(r, r.replace(a, b), 1)
    assert new == (P / 'proposed-main.rs.txt').read_bytes()
    assert blob(new) == '9ef9f1ba005f20149716b22346d45c5501b62e87'
    newrows = {m[1].decode(): m[0] for m in re.finditer(pat, new)}
    assert rows.keys() == newrows.keys()
    changed = [n for n in rows if rows[n] != newrows[n]]
    assert set(changed) == set(swaps)
    for n in rows:
        assert len(re.findall(rb'\w+:', rows[n])) == 11
        assert re.sub(rb'numerators: 0x[0-9a-f]+', b'numerators: TOKEN', rows[n]) == re.sub(rb'numerators: 0x[0-9a-f]+', b'numerators: TOKEN', newrows[n])
    prior = ['xs_qed_eemm_L1', 'xs_qed_eedd_L1', 'xs_qed_eeee_1pi_L1', 'qed_eemm_L1', 'ctorder_qcd_uugg_L1', 'sm_eemm_L0', 'ym_g6_L0', 'xs_qed_eemm_L2', 'sm_eeww_L0']
    assert all(rows[n] == newrows[n] for n in prior)
    sc = scopes((W / 'docs/developers.toml').read_bytes(), lambda p: (W / p).read_bytes())
    assert len(sc) == 118
    main = W.parent / 'main/.git'
    write('start.json', {'base': BASE, 'time': datetime.datetime.now(datetime.UTC).isoformat(), 'author': 'Lucien Huber <im@lcnbr.ch>', 'tracked_paths': paths, 'main_git': {p: sha(main / p) if (main / p).exists() else None for p in ['HEAD', 'index', 'packed-refs']}, 'refs': git('show-ref').decode(), 'sealed_packets': sealed, 'CI': False})
    write('source-check.json', {'base': BASE, 'before_blob': blob(old), 'after_blob': blob(new), 'changed_rows': changed, 'tokens': 3, 'literal_Rows': 149, 'fields_each': 10, 'other_nine_fields_equal': True, 'every_other_byte_equal': True, 'prior_nine': {n: hashlib.sha256(rows[n]).hexdigest() for n in prior}, 'scope_bindings': sc, 'scope_anchors_or_digests_changed': 0, 'host_API_browser_and_all_other_base_files_pinned': True})
    print('Guard PASS:', len(paths), 'base files;', len(sc), 'scopes;', len(rows), 'complete Rows; exact three-token result verified. No tracked edit yet.')
