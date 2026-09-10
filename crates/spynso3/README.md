# Spenso Python adapter

`spynso3` registers Spenso's native Python interface as
`symbolica.community.spenso`. It is an adapter crate for the Symbolica community
module assembly, not a standalone `spenso` Python distribution.

The canonical user guide is the rendered
[Spenso Python workflow](https://alphal00p.github.io/gammaloop/products/spenso/latest/guides/python/).
It covers tensor construction, reusable libraries, network execution, symbolic
evaluation, parallelism, failure boundaries, and links to the generated API.
Exact signatures and defaults are in the
[generated Python reference](https://alphal00p.github.io/gammaloop/products/spenso/latest/reference/python/spynso3/);
the checked [`spynso3.pyi`](../../docs/api/python/spynso3.pyi) remains a tooling artifact.

## Availability

An installed Symbolica package exposes Spenso only when its assembly includes
this community module. Verify the actual environment with:

```bash
python -c "import symbolica.community.spenso as spenso; print(spenso.__name__)"
```

For a custom build, add this crate to the
[`symbolica-community`](https://github.com/benruijl/symbolica-community)
assembly and invoke its `SymbolicaCommunityModule` registration. Building this
Rust crate alone does not add the module to an existing Python installation.

Contributor policy and documentation commands are in
[`CONTRIBUTING.typ`](../../CONTRIBUTING.typ).
