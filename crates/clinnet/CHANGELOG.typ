= Changelog

All notable changes to this project will be documented in this file.

The format is based on #link("https://keepachangelog.com/en/1.0.0/")[Keep a Changelog],
and this project adheres to #link("https://semver.org/spec/v2.0.0.html")[Semantic Versioning].

== Unreleased

=== Changed

- The bundled Linnest package declares `MIT AND Apache-2.0 AND BSD-3-Clause`. Its new
  `LICENSE.networkx` holds the BSD-3-Clause notice for the NetworkX code that the
  `linnest.wasm` layout seed translates, previously only part of `LICENSE.ec-layout`.

=== Removed

- The bundled Linnest package no longer carries the unused `ec-layout.wasm` plugin, which
  linked the GPL-licensed OGDF, or its `LICENSE.ec-layout`. Layout seeds come from
  `linnest.wasm` as before.

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.7...clinnet-v0.1.8")[0.1.8] - 2025-12-04

=== Other

- redo cli to have subcommands

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.6...clinnet-v0.1.7")[0.1.7] - 2025-12-04

=== Other

- properly forward inputs for a single figure

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.5...clinnet-v0.1.6")[0.1.6] - 2025-12-03

=== Other

- Add parallel processing and directional constraints
- add inputs

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.4...clinnet-v0.1.5")[0.1.5] - 2025-11-28

=== Other

- change output name and make it short, and add directional groups

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.3...clinnet-v0.1.4")[0.1.4] - 2025-11-26

=== Other

- Make figure and grid templates optional with defaults

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.2...clinnet-v0.1.3")[0.1.3] - 2025-11-26

=== Other

- don't overwrite existing templates
- strip quotes strings for typst

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.1...clinnet-v0.1.2")[0.1.2] - 2025-11-24

=== Other

- update symbolica

== #link("https://github.com/alphal00p/linnet/compare/clinnet-v0.1.0...clinnet-v0.1.1")[0.1.1] - 2025-11-12

=== Other

- remove unnecessary outputs
- better templates
