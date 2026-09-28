= Directed self-loop symmetry reproducer

The installed tensor-consolidation checks found a discrepancy before tensor
algebra: a one-vertex closed Dirac-fermion tadpole has overall factor $-1/2$,
where the independent Wick count in
`crates/feynkit-py/tests/installed_tadpole_normalization.py` requires $-1$.
The generated graph carries `AutG(2)` and one closed-fermion-loop minus sign.
No numerator simplification or expansion is involved in this comparison.

#link("directed_self_loop.py")[The minimal Python program] uses only Symbolica's
public graph API. Both a directed self-loop and an undirected self-loop return
an automorphism group size of two, with or without a uniquely marked external
vertex. An oriented fermion line distinguishes its incoming and outgoing ends;
its one allowed Wick pairing has no extra factor of two. The undirected scalar
control has two interchangeable ends and does require that factor.

This was measured with Symbolica
`6a96c9d77b217ea394770524c119e0613d4d9ef5`, whose installed lockfile uses
registry `graphica 3.0.0` (checksum
`afa7211a7cfdba21893f634868e19459e5d1882e5957f98d7c084875a6e70ab6`).
The same four outputs were reproduced with the preserved `939c4de` build,
using the identical program. This issue therefore predates the update to
`6a96c9d`; it is not a regression from the measured expansion changes.
The tensor implementation and its assertions are unchanged by this reproducer;
no generator compensation for the upstream edge symmetry is applied.

#link("observation.json")[The observation record] binds the exact tested program,
installed core, commands, dependency versions and output. The portable program
has the same four graph constructions, omitting only JSON/provenance bookkeeping.
