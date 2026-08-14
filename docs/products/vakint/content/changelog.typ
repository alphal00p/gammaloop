#import "../../shared.typ": release-note, source-link

#let changelog = [
= Releases and change history

#release-note([
  Vakint does not currently have a canonical `CHANGELOG.md` in this repository. The selected
  release is identified by `crates/vakint/Cargo.toml`, repository tags, and the exact source
  revision. This manual must not synthesize historical entries from commits or GammaLoop's
  release history.
])

A future curated changelog should distinguish topology-library changes, canonicalization and
tensor-reduction changes, analytic backend changes, pySecDec integration, normalization and
precision changes, and Rust/Python API migrations. These categories can affect reproducibility
even when a high-level method name is unchanged.

Documentation may publish immutable snapshots for Vakint tags and a moving `latest` channel.
Every snapshot records the crate version, feature set, FORM and pySecDec prerequisites, and
repository revision. Until curated release notes exist, verify behavior against the selected
source and tests beginning at
#source-link("crates/vakint/README.md", label: "the maintained README").
]
