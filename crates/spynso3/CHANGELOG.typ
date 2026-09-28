= Changelog

All notable changes to this project will be documented in this file.

The format is based on #link("https://keepachangelog.com/en/1.0.0/")[Keep a Changelog],
and this project adheres to #link("https://semver.org/spec/v2.0.0.html")[Semantic Versioning].

== Unreleased

- MathML output loads the full STIX Two Math font from a pinned CDN URL when
  it is not installed locally. Expressions share the font download and contain
  no embedded font bytes.
- `symbolica.community.spenso.load_math_font()` returns a style element with
  the bundled font for offline notebooks and standalone HTML exports. Display
  it once before expressions in the same HTML document:

  ```python
  from IPython.display import HTML, display
  from symbolica.community import spenso

  display(HTML(spenso.load_math_font()))
  ```

  Keep this setup output when saving a notebook. Isolated output frames require
  their own setup; normal online output needs no setup cell.

== #link("https://github.com/alphal00p/spenso/compare/spynso3-v0.1.1...spynso3-v0.1.2")[0.1.2] - 2026-01-14

=== Other

- update versions

== #link("https://github.com/alphal00p/spenso/compare/spynso3-v0.1.0...spynso3-v0.1.1")[0.1.1] - 2025-12-15

=== Fixed

- fix sum normalization in network parsing
- fix warnings
- fix `add_assign` with sparse tensors
- fix projp
- fix param docs and iter return type
- fix `pyo3_stubs`

=== Other

- Update implementation to pyo3-stub-gen 0.17
- Update symbolica to 1.1.0
- update linnet
- add readmes
- update to 0.17 `pyo3_stub_gen`
- update linnet and make classattr to classmethods
- move initialize to api
- Release spenso-macros 0.3
- update to symbolica 1.0
- update sympolica to 0.20
- Properly precontract scalars
- Add back pyo3-stub-gen-derive
- update to current dev symbolica and add pattern for pslash + m conjugate
- first try conj
- update python api
- add function generic key
- again
- return any because python is dumb
- update `__iter__` method
- use numpy doc format
- add overloaded pyi stub gen
