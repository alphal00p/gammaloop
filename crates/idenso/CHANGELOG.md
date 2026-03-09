# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.5](https://github.com/alphal00p/spenso/compare/idenso-v0.2.4...idenso-v0.2.5) - 2026-01-08

### Fixed

- fix metric initialization and update to symbolica 1.2
- fix initialize to also load ETS

### Other

- update linnet

## [0.2.4](https://github.com/alphal00p/spenso/compare/idenso-v0.2.3...idenso-v0.2.4) - 2025-12-15

### Fixed

- fix canonization in presence of duals
- fix gamma normalisation
- fix warnings
- fix add_assign with sparse tensors
- fix projp

### Other

- Import initialize function
- Update symbolica to 1.1.0
- update linnet
- align to atom in argument of cof for Nc
- update linnet
- add readmes
- formatting and add ref for parametric
- update to 0.17 pyo3_stub_gen
- update linnet and make classattr to classmethods
- proper imports and feature flags for python
- move initialize to api
- Release spenso-macros 0.3
- update to symbolica 1.0
- update sympolica to 0.20
- Properly precontract scalars
- Add back pyo3-stub-gen-derive
- update to current dev symbolica and add pattern for pslash + m conjugate
- update to linnet 0.13.1
- working canonize (up to functions) and dirac adjoint
- spenso canonize buggy
- robust_conj but with lots of gamma_0s
- first try conj
- working pow execution with wrapping
- add function generic key

## [0.2.2](https://github.com/alphal00p/spenso/compare/idenso-v0.2.1...idenso-v0.2.2) - 2025-08-18

### Fixed

- fix cook indices
- fix clippy errors

### Other

- update linnet to crates
- all tests pass
- Fix Multiple Heads Bug
- update linnet
- update linnet
- update linnet
- update parse problem and add spenso:: to all symbols
- make all symbols take spenso namespace
- operate directly on coefficient list
- remove expand

## [0.2.1](https://github.com/alphal00p/spenso/compare/idenso-v0.2.0...idenso-v0.2.1) - 2025-07-15

### Fixed

- fix gamma algebra again

### Other

- update versions

## [0.2.0](https://github.com/alphal00p/spenso/compare/idenso-v0.1.0...idenso-v0.2.0) - 2025-07-15

### Fixed

- fix wrap_indices
- fix gamma algebra≈

### Other

- allow arbitrary arguments for to dots
