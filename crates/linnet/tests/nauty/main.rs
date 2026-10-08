//! Linnet's canonizer and graph generator against nauty 2.9.3's test suite.
//!
//! Optimised builds, including CI's `ci-optim` profile, run every family at full size. Plain
//! `cargo test` builds (opt-level 0, binaries under `target/debug`) run reduced sizes, see
//! [`runs::full_size`]: tournaments on 6 vertices only, the radius-6 trees of run 68, generated
//! trees up to 9 vertices, the first base of run 70, the radius-3 bases of run 95 and the cubic
//! bases of run 96. The tests in [`oracle`] are ignored and need a nauty build in `NAUTY_DIR`.

mod dreadnaut;
mod formats;
mod graph;
mod oracle;
mod runs;
mod utilities;
