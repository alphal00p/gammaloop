# nauty test-suite fixtures

These files come from, or were produced with, **nauty and Traces 2.9.3** by Brendan McKay and
Adolfo Piperno (<https://pallini.di.uniroma1.it>, <https://users.cecs.anu.edu.au/~bdm/nauty/>).
nauty is licensed under the Apache License 2.0; its notice is in [`COPYRIGHT`](COPYRIGHT) and the
licence text in [`LICENSE-2.0.txt`](LICENSE-2.0.txt). They are used only by Linnet's
`tests/nauty` integration test, which checks Linnet's canonizer and graph generator against nauty's
`runalltests` suite, and are excluded from the published `linnet` crate (Linnet itself is MIT).

## Copied unchanged from the nauty 2.9.3 distribution

- `runalltests.in`: the test driver; the tests parse its 98 runs.
- `nautest1.dre`, `nautest2.dre`, `nautest3.dre`: the dreadnaut scripts.
- `nautest1a.ans` ... `nautest3c.ans`: the dreadnaut transcripts (dense, sparse, Traces).
- `nautest{a,b1,b2,b3,d,e,h,i,j,k,m,n,q,s,y,ff,ii}.ans`: answers of the utility pipelines.

## Recorded from a nauty 2.9.3 build

Inputs of pipelines whose answers are too large to keep (or are labelled), written by the first
stage of the pipeline:

| file | command | for |
|---|---|---|
| `directg-bases.g6` | `geng -c -qp 8 9` | `nautestf.ans` |
| `vcolg-bases.g6` | `genbg -q 4 3 10` | `nautestg.ans` |
| `watercluster2-bases.g6` | `geng -Cq 8 12` | `nautestff.ans` |
| `distgraphg-bases.s6` | `genspecialg -q -X12 -X44 -X23 -X68 -X100` | `nautestii.ans` |
| `multig-bases.txt` | `geng -qd3D4 8`, each line followed by the line count of `multig -qTm2` on it | `nautestgg.ans` |

The expected group sizes, orbit counts and per-base class counts written into the tests were
read from the same build with `countg -V -1 --a,o,O` and from the answers above.

Optimised builds (release, or CI's `ci-optim`) run every family at full size; unoptimised
`target/debug` builds run reduced sizes, listed in `tests/nauty/main.rs`.

The tests that need nauty itself (see `tests/nauty/oracle.rs`) are ignored by default and fail
without a nauty build; run them with
`NAUTY_DIR=<nauty build> cargo nextest run --release -p linnet --test nauty --run-ignored only oracle::`.
