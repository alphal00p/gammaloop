# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.6](https://github.com/alphal00p/spenso/compare/spenso-v0.5.5...spenso-v0.5.6) - 2026-01-08

### Fixed

- fix python_stub_gen feature
- fix metric initialization and update to symbolica 1.2

### Other

- update linnet

## [0.5.5](https://github.com/alphal00p/spenso/compare/spenso-v0.5.4...spenso-v0.5.5) - 2025-12-15

### Fixed

- fix canonization in presence of duals
- fix sum normalization in network parsing
- fix warnings
- fix overflowing slot order
- fix add_assign with sparse tensors

### Other

- Update symbolica to 1.1.0
- update linnet
- add depth to parser, to stop at first mul
- remove more prints
- remove prints
- update linnet
- add readmes
- formatting and add ref for parametric
- update to 0.17 pyo3_stub_gen
- update linnet and make classattr to classmethods
- Release spenso-macros 0.3
- update to symbolica 1.0
- update sympolica to 0.20
- Properly precontract scalars
- Add back pyo3-stub-gen-derive
- Updated spenso with minor canonical string update for symbolica f2eb0c1ddcc61c342f2ff9616c2d129141d433dc
- Refactor coefficient conversions to support generic float types
- update to current dev symbolica and add pattern for pslash + m conjugate
- update to linnet 0.13.1
- working canonize (up to functions) and dirac adjoint
- spenso canonize buggy
- robust_conj but with lots of gamma_0s
- first try conj
- update python api
- working pow and fun execution
- working pow execution with wrapping
- support for executing pow and functions on tensors
- implement function library trait
- add function generic key

## [0.5.3](https://github.com/alphal00p/spenso/compare/spenso-v0.5.2...spenso-v0.5.3) - 2025-08-18

### Other

- update linnet to crates
- all tests pass
- Fix Multiple Heads Bug
- implement ExportNumber for spenso::Complex
- to Atom for ExpandedIndex
- remove printouts
- update linnet
- update linnet
- update linnet
- add infinite loop error test
- update parse problem and add spenso:: to all symbols
- remove expand

## [0.5.2](https://github.com/alphal00p/spenso/compare/spenso-v0.5.1...spenso-v0.5.2) - 2025-07-15

### Other

- remove printout

## [0.5.1](https://github.com/alphal00p/spenso/compare/spenso-v0.5.0...spenso-v0.5.1) - 2025-07-15

### Fixed

- fix gamma algebra≈

### Other

- allow arbitrary arguments for to dots
