= Standalone WASM timing adaptation

This directory vendors Clarabel 0.11.1 from its crates.io distribution under the
upstream Apache-2.0 license. Numerical algorithms and native timing are unchanged.
Standalone WASM builds have no JavaScript clock: timers report zero elapsed time,
and settings validation rejects every finite time limit. ImPrEd uses an infinite
time limit and its existing iteration limit. The WASM `web-time` dependency is
removed; native builds retain `std::time::Instant`.
