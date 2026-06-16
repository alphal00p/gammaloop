- If running into a linking issue complaining about missing `__emul...` symbol, then it may be a macos specific issue, that can be fixed by building with the following environment variable change `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` (see impact of this in file `build.rs`.
- Physical PySecDec integration tests can be slow and remain manual only. They
  are ignored and excluded from automatic local and CI selection. See the
  [README](README.md) for the explicit four-test manual command. Pure Rust
  adapter tests remain automatic because they do not launch PySecDec.
- A normal Vakint run retaining temporary outputs for diagnosis is:

RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo test --package vakint --no-default-features
