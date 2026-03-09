<div align="center">

![logo](https://raw.githubusercontent.com/alphal00p/spenso/master/spensologo.svg)


[![Documentation](https://img.shields.io/docsrs/spenso/latest)](https://docs.rs/spenso/latest/spenso/)
[![crates.io](https://img.shields.io/crates/v/spenso.svg)](https://crates.io/crates/spenso)
[![Build Status](https://github.com/alphal00p/spenso/actions/workflows/ci.yml/badge.svg)](https://github.com/alphal00p/spenso/actions/workflows/ci.yml)
[![codecov](https://codecov.io/github/alphal00p/spenso/graph/badge.svg?token=ST0XA54QSF)](https://codecov.io/github/alphal00p/spenso)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/alphal00p/spenso)
</div>


The spenso crate is the foundational tensor manipulation library that provides tensor network representation, automatic index matching, data storage, and execution capabilities. It handles the core mathematical operations on tensors with abstract indices, enabling automatic contraction and simplification of tensor expressions.

The spenso library provides:

1. Network System: Computation graph representation with Network, NetworkGraph, and execution strategies
1. Data Storage: Flexible DataTensor enum over DenseTensor and SparseTensor with complex number support
1. `symbolica` Integration: Parsing from symbolica expressions with automatic structure recognition
1. Execution Engine: Multiple strategies for network evaluation with automatic index matching
1. Type-Safe Abstractions: Rich trait system separating structure from data

One of the nicest features of the network, is that it can be parsed from a `symbolica` expression directly, allowing for seamless integration between symbolic manipulation and concrete tensor evaluation.
