# Drawing benchmark

See [the canonical instructions and provenance](instructions.typ) and
[the fixture manifest](manifest.typ). Run from the repository root with
`python3 tests/drawing-bench/run.py --help`.

The archives and original acceptance measurements come from branch
`983da1d9310bf885fcebe5ba4cb8b47166f82fa0`. They are historical evidence,
not refreshed performance results. `native.py` now uses typed `RenderSettings`,
`DiagramStyle` and `DiagramView.to_svg()`; this timing boundary is recorded in
new metadata and has not been benchmarked as part of the refresh.
