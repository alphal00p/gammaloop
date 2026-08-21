#!/usr/bin/env python3
"""Deployed-app coverage sweep: does the full gammaloop -> oneloop reduce pipeline
handle real Standard-Model processes at scale without breaking?

This is a *coverage / robustness* benchmark, complementary to the cross-engine
*correctness* checks (`crosscheck.py`, `verify_*.py`). It drives the live
FeynmanEngine API (https://gammaloop.hirschi.lu/feynmangraph) exactly as a user
would: for each process it POSTs to `/api/generate-amp` (one loop), then reduces
an evenly-spaced sample of the emitted diagrams via `/api/reduce`, and tallies:

  - reduced OK (a coeff x master result was returned),
  - graceful `zero_numerator` (the diagram's contracted numerator vanishes),
  - `unsupported` "walls" (a reducer degenerate-Cayley panic surfaced gracefully),
  - HTTP errors / timeouts,
  - which master classes (A0/B0/C0/D0) each process exercises.

It needs no external engines (Python stdlib only) and talks only to the public
deployment, so it is safe to run from any machine (network I/O, negligible local
memory). Requests are issued strictly sequentially to avoid contending on the
server's per-request gammaloop subprocess.

Usage:
    python3 app_process_sweep.py [diagrams_per_process=20] [base_url]

Example (reduce up to 40 diagrams per process):
    python3 app_process_sweep.py 40
"""
import gzip
import json
import re
import sys
import time
import urllib.error
import urllib.request

BASE = "https://gammaloop.hirschi.lu/feynmangraph/api"

# Standard-Model one-loop processes spanning every topology the reducer supports:
# self-energies (bubble), V -> f f-bar (triangle), Higgs, 4-fermion + QCD (box),
# and loop-induced / high-rank (light-by-light is a rank-4 box). Each tuple is
# (label, initial_state, final_state, [coupling_orders to try in order]).
Q2 = [{"QED": 2}]
QCD = [{"QED": 2}, {"QED": 2, "QCD": 2}]
PROCESSES = [
    # --- self-energies (2-point / bubble) ---
    ("a>a", ["a"], ["a"], Q2), ("z>z", ["z"], ["z"], Q2), ("w+>w+", ["w+"], ["w+"], Q2),
    ("h>h", ["h"], ["h"], Q2), ("g>g", ["g"], ["g"], Q2), ("e->e-", ["e-"], ["e-"], Q2),
    ("u>u", ["u"], ["u"], Q2), ("d>d", ["d"], ["d"], Q2), ("b>b", ["b"], ["b"], Q2),
    ("t>t", ["t"], ["t"], Q2), ("mu>mu", ["mu-"], ["mu-"], Q2), ("ta>ta", ["ta-"], ["ta-"], Q2),
    # --- V -> f f-bar (3-point / triangle) ---
    ("z>e+e-", ["z"], ["e+", "e-"], Q2), ("z>dd~", ["z"], ["d", "d~"], Q2),
    ("z>uu~", ["z"], ["u", "u~"], Q2), ("z>bb~", ["z"], ["b", "b~"], Q2),
    ("w+>ud~", ["w+"], ["u", "d~"], Q2), ("a>e+e-", ["a"], ["e+", "e-"], Q2),
    ("a>dd~", ["a"], ["d", "d~"], [{"QED": 1, "QCD": 2}, {"QED": 3}]),
    # --- Higgs ---
    ("h>aa", ["h"], ["a", "a"], Q2), ("h>gg", ["h"], ["g", "g"], Q2),
    ("h>bb~", ["h"], ["b", "b~"], Q2), ("h>ta+ta-", ["h"], ["ta+", "ta-"], Q2),
    ("h>mu+mu-", ["h"], ["mu+", "mu-"], Q2),
    # --- 4-fermion (4-point / box) ---
    ("e+e->mu+mu-", ["e+", "e-"], ["mu+", "mu-"], Q2),
    ("e+e->ta+ta-", ["e+", "e-"], ["ta+", "ta-"], Q2),
    ("e+e->bb~", ["e+", "e-"], ["b", "b~"], Q2), ("uu~>dd~", ["u", "u~"], ["d", "d~"], Q2),
    ("dd~>dd~", ["d", "d~"], ["d", "d~"], Q2), ("uu~>tt~", ["u", "u~"], ["t", "t~"], Q2),
    # --- QCD multi-parton (box) ---
    ("gg>gg", ["g", "g"], ["g", "g"], Q2), ("gg>ddx", ["g", "g"], ["d", "d~"], Q2),
    ("uux>gg", ["u", "u~"], ["g", "g"], Q2), ("gu>gu", ["g", "u"], ["g", "u"], Q2),
    ("ddx>gg", ["d", "d~"], ["g", "g"], Q2), ("uux>ddx", ["u", "u~"], ["d", "d~"], Q2),
    # --- loop-induced / high-rank tensor ---
    ("gg>h", ["g", "g"], ["h"], Q2), ("gg>hh", ["g", "g"], ["h", "h"], Q2),
    ("aa>aa(LbL)", ["a", "a"], ["a", "a"], [{"QED": 4}]),
    ("aa>e+e-", ["a", "a"], ["e+", "e-"], Q2),
]

_MASTER = re.compile(r"\b(A0|B0|C0|D0)\(")


def _post(path, body, timeout=140):
    data = json.dumps(body).encode()
    req = urllib.request.Request(
        BASE + path, data=data,
        headers={"Content-Type": "application/json", "Accept-Encoding": "gzip"},
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read()
            if resp.headers.get("Content-Encoding") == "gzip":
                raw = gzip.decompress(raw)
            return resp.status, json.loads(raw.decode())
    except urllib.error.HTTPError as exc:
        return exc.code, {"_err": exc.read().decode(errors="replace")[:150]}
    except Exception as exc:  # noqa: BLE001 (network sweep: report, don't crash)
        return -1, {"_err": f"{type(exc).__name__}: {exc}"[:150]}


def _masters(raw):
    return sorted(set(_MASTER.findall(raw)))


def _generate(initial, final, orders_list):
    for orders in orders_list:
        st, body = _post("/generate-amp", dict(
            initial_state=initial, final_state=final, coupling_orders=orders,
            loop_count=1, model_id="sm", theory_id="sm", max_diagrams=200))
        if st == 200 and body.get("count", 0) > 0:
            return orders, body["diagrams"]
    return None, []


def _even_sample(n, k):
    if n <= k:
        return list(range(n))
    return sorted({round(i * (n - 1) / (k - 1)) for i in range(k)})


def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    global BASE
    if len(sys.argv) > 2:
        BASE = sys.argv[2].rstrip("/")

    t0 = time.time()
    tot_diag = tot_ok = tot_err = tot_wall = tot_zero = 0
    d0_procs, anomalies = [], []
    print(f"# deployed-app coverage sweep  ({time.strftime('%Y-%m-%d %H:%M:%S')})")
    print(f"# base={BASE}  diagrams_per_process<={budget}\n")

    for label, initial, final, orders_list in PROCESSES:
        orders, diagrams = _generate(initial, final, orders_list)
        if not diagrams:
            print(f"{label:16} GEN-FAIL")
            anomalies.append((label, "generate"))
            continue
        idxs = _even_sample(len(diagrams), budget)
        mset, ok, err, wall, zero, maxc = set(), 0, 0, 0, 0, 0
        for i in idxs:
            status, resp = _post("/reduce", diagrams[i])
            tot_diag += 1
            if status != 200:
                err += 1
                tot_err += 1
                continue
            reason, raw = resp.get("reason"), resp.get("raw", "") or ""
            if reason == "unsupported":
                wall += 1
                tot_wall += 1
            elif reason == "zero_numerator":
                zero += 1
                tot_zero += 1
            elif reason is None:
                ok += 1
                tot_ok += 1
                mset |= set(_masters(raw))
                maxc = max(maxc, len(raw))
        if "D0" in mset:
            d0_procs.append(label)
        if err or wall:
            anomalies.append((label, f"err={err} wall={wall}"))
        flag = (f" ERR={err}" if err else "") + (f" WALL={wall}" if wall else "")
        print(f"{label:16} n={len(diagrams):3} sampled={len(idxs):2} | "
              f"ok={ok:2} masters={''.join(sorted(mset)) or '-':8} zero={zero} "
              f"maxchars={maxc}{flag}")

    print(f"\n=== SUMMARY ({time.time() - t0:.0f}s) ===")
    print(f"processes={len(PROCESSES)}  diagrams_reduced={tot_diag}  ok={tot_ok}  "
          f"zero_numerator={tot_zero}  http_err={tot_err}  walls(unsupported)={tot_wall}")
    print(f"D0/box coverage: {len(d0_procs)} processes -> {d0_procs}")
    print(f"anomalies: {anomalies if anomalies else 'NONE'}")


if __name__ == "__main__":
    main()
