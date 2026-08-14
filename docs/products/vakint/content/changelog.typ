#import "../../shared.typ": callout, source-link

#let changelog = [
= Versions and compatibility

#callout("Before upgrading", [
  Pin the Vakint version together with FORM, pySecDec, MATAD, and FMFT versions used for a
  calculation. Changes to topology matching, normalization, or a backend can change a result even
  when the high-level method call remains the same.
])

The Rust crate and the Python interface included with Symbolica can be released on different
schedules. Record both versions when sharing a calculation or reporting a problem.

There is not yet a standalone Vakint changelog. When evaluating an upgrade, check the selected
release's API reference and
#source-link("crates/vakint/README.md", label: "Vakint README"), then compare the following areas:

- supported topologies and canonical momentum routing;
- normalization conventions and epsilon-expansion defaults;
- tensor reduction and analytic backend behavior;
- pySecDec integration and numerical precision;
- Rust and Python method signatures and feature requirements.
]
