"""Check the iterated-residue Matsubara sums of matsubara.py against direct sums.

At random spatial momenta (zero and non-zero external momentum) and non-zero bosonic
external frequencies omega_K = 2 pi T m, R(i omega_K) must equal the Euclidean sum.
The direct sum truncated at |n_l| <= nmax has a tail ~ nmax^-p (p ~ 5); three
truncations nmax/4, nmax/2, nmax estimate p and remove it (Richardson). The tolerances
are about ten times the largest deviation measured with this seed: 4.4e-15 (DT,
roundoff) and 3.4e-10 (TBT, remaining truncation error of the triple sum; 1.6e-11 at
nmax=80, which is slow on a loaded machine).

Usage:  python check_matsubara.py      (exits nonzero on failure)
"""

import sys

import numpy as np

# isort: split
from matsubara import Compiled, brute_force, build_terms, get_graph, line_energies

# graph: (nmax, relative tolerance)
CASES = {"dt": (300, 5e-14), "tbt": (48, 3e-9)}


def main():
    rng = np.random.default_rng(1)
    failed = False
    for name, (nmax, tol) in CASES.items():
        g = get_graph(name)
        C = Compiled(g, build_terms(g))
        print(f"{name}: {C.n_terms} terms, {len(C.thresholds)} threshold functions")
        rels = []
        for trial in range(3):
            loops = rng.normal(size=(1, g.n_loops, 3)) * 0.7
            kvec = rng.normal(size=3) * 0.5 if trial else np.zeros(3)
            E = line_energies(g, loops, kvec)
            for T in (0.5, 0.3):
                for m in (1, 2):
                    val = C.evaluate(E, [1j * 2 * np.pi * T * m], T)[0, 0]
                    b0, b1, b2 = (
                        brute_force(g, E[0], T, m, n)
                        for n in (nmax // 4, nmax // 2, nmax)
                    )
                    with np.errstate(divide="ignore", invalid="ignore"):
                        p = np.log2(abs((b0 - b1) / (b1 - b2)))
                    # without a resolved tail the truncated sum has already converged
                    if np.isfinite(p) and p > 1:
                        direct = b2 + (b2 - b1) / (2**p - 1)
                    else:
                        direct = b2
                    rel = abs(val - direct) / abs(direct)
                    rels.append(rel)
                    print(
                        f"  trial {trial} T={T} m={m}: residues {val.real:.13e} "
                        f"direct {direct:.13e} (tail p={p:.2f}) rel {rel:.1e}"
                    )
        worst = np.max(rels)  # nan propagates and fails
        ok = worst < tol
        failed |= not ok
        print(f"{name}: worst rel {worst:.1e} {'<' if ok else '>='} {tol:.0e}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
