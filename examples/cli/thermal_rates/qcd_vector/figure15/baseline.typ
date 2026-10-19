= Figure 15 historical baseline

Every native integration, timing measurement, and CPU profile in these dated result
directories was produced with source commit
`116767bd3ed30491a2b15559d228171d20082abe` on `itphlies`, using native
binary SHA256
`9f4f7bc9b5f5d1caffa8235cb88b93429b90677aabe9e05420dd155e03e8de24`.
The measurements predate `75a0ced4e`, which changes signed radial
certification to try binary64 bounds first. A later rebase of these result
commits does not make the measurements results for that optimization.
#link("baseline.json")[The baseline record] makes this attribution explicit.

The #link("comparison.typ")[run card and normalization],
#link("pilot.typ")[pilot],
#link("full-scan-20261019/report.typ")[full-scan report],
#link("vacuum-timing-20261018/report.typ")[vacuum comparison], and
#link("profiling-20261018/report.typ")[CPU profile] describe the investigation.
The full scan ended within its 60-worker, twelve-hour limit. Seven rest-frame
points reached the declared reference-normalized 0.5% statistical precision;
moving-frame convergence remains incomplete. Near-light-cone results have
large errors, and the exact $E=3k$ author-code target remains suspect.
Published vector curves are graphical references, not numerical targets.

The committed files retain the summaries, final raw numerical results,
reference values, reports, and selected figures. Final native and author
point objects are consolidated without changing their numerical values.
The compact full-scan manifest retains seeds, runtime settings, input hashes,
stopping rules, and final statuses. The original path-bearing manifests,
executed source snapshots, all checkpoints and logs, raw/recovered perf
recordings, generated states, source PDF, and baseline executable are in the
verified archive identified by #link("archive.json")[the archive record].

The archive is stored outside Git on `itphlies` under
`/shared/kaapo/gammaloop/artifacts/figure15-baseline-116767bd3-20261019T094021Z/`.
Its `SHA256SUMS` authenticates `artifacts.tar.zst` and the per-file content
manifest. Zstandard integrity testing and a full tar comparison with the
original files passed before curation. An unpacked copy of the original
checkout artifacts is retained there as well.

From that archive directory, verify and restore into a separate destination:

```sh
sha256sum -c SHA256SUMS
mkdir restored
tar --zstd -xf artifacts.tar.zst -C restored
```

Archive members preserve the original absolute paths with the leading slash
removed. For example, scratch runs appear below `restored/tmp/` and the
original checkout artifacts below
`restored/shared/kaapo/gammaloop/gammaloop-thermal-rates/examples/cli/thermal_rates/qcd_vector/figure15/`.
Original absolute paths inside historical records are provenance, not
requirements for the portable reproduction tools.
