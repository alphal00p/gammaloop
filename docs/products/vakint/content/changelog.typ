#import "../../shared.typ": release-note, source-link

#let changelog = [
= Version history

#release-note([
  Vakint does not yet have a standalone changelog. The package metadata, tagged repository
  history, and API reference are the available version record.
])

== Components and version sources

- Rust package: #source-link("crates/vakint/Cargo.toml", label: "vakint package metadata")
- Release record: #link("https://github.com/alphal00p/vakint/tags")[Vakint repository tags]
- Usage and external tools: #link("guides/evaluation/")[matching and evaluation guide]

The Rust crate and the Python interface included with Symbolica can be released on different
schedules. Record both versions when sharing a calculation or reporting a problem.

== Upgrade checklist

Pin the Vakint version together with FORM, pySecDec, MATAD, and FMFT versions used for a
calculation. Then compare:

- supported topologies and canonical momentum routing;
- normalization conventions and epsilon-expansion defaults;
- tensor reduction and analytic backend behavior;
- pySecDec integration and numerical precision;
- Rust and Python method signatures and feature requirements.

== Reproducibility record

Record the Rust crate, Symbolica/Python module, FORM, pySecDec, MATAD, and FMFT versions together
with the selected topology and normalization settings.
]
