# Compression during cache publication

On 9 September 2026, a local single-artifact Nix binary-cache export measured:

| Outer NAR compression | Wall time | CPU time | Stored bytes |
|---|---:|---:|---:|
| XZ, default single thread | 217.35s | 214.69s | 472,995,704 |
| XZ, parallel | 34.56s | 212.43s | 474,632,864 |
| Zstd, parallel | 1.74s | 2.22s | 480,617,203 |

The input was the same 483,715,824-byte NAR containing one Crane
`target.tar.zst`. The Zstd export was 1.61% larger than single-threaded XZ and
1.26% larger than parallel XZ. This used eight fixed CPUs (8–15), a shared host,
local file caches and no network. Reference objects were present before timing;
`--no-recursive` limited each timed export to the artifact. The variants used
separate destination caches and preserved the source bytes.

This establishes a potential CPU cost inside publication intervals. It does not
establish which compression settings NixCI workers currently use, or how much of
a remote interval is compression, reference checking, storage or network time.
Ask Syd whether worker publication uses XZ, whether parallel compression is enabled,
and whether Zstd can be used for these already-compressed Crane archives.
The [Nix binary-cache settings](https://nix.dev/manual/nix/2.35/store/types/local-binary-cache-store.html)
document both methods and parallel compression.

Two earlier harness attempts are retained separately: a recursive copy included
598 closure paths and was stopped; an unseeded single-path copy failed because
its references were absent from the destination. Neither is used in this timing
table. The successful JSON and CSV results are under
`/tmp/gammaloop-ci-validation/cache-compression/seeded-single-path/`.

An unauthenticated read of NixCI cache metadata returned HTTP 401, so it did not
reveal current cache visibility or compression settings. This is an access
limitation of the probe, not evidence of a worker authentication problem.
