# FeynKit Python community module

`feynkit-py` provides the Rust side of `symbolica.community.feynkit`. It is an
`rlib` implementing `SymbolicaCommunityModule`; it deliberately does not export
a PyO3 extension module or build a wheel by itself. Linking it as a separate
extension would create a second Symbolica kernel, whose expression types and
global state would not be interchangeable with the installed `symbolica`
package.

## Symbolica community host integration

The combined
[`symbolica-community`](https://github.com/symbolica-dev/symbolica-community)
host owns the single `symbolica.core` extension. Integrate FeynKit there in the
same way as `spynso3`:

1. Add the Rust library to the host dependencies. Track `main` while integrating
   the change, then replace `branch` with a pinned `rev` for a release build:

   ```toml
   [dependencies]
   feynkit-py = { git = "https://github.com/alphal00p/gammaloop", branch = "main" }
   ```

2. Include FeynKit's inventory when the host generates stubs:

   ```toml
   [features]
   python_stubgen = [
     # existing community modules...
     "feynkit-py/python_stubgen",
   ]
   ```

3. Register the module in the host's `#[pymodule] fn core(...)`, after
   `create_symbolica_module(m)?`:

   ```rust
   register_module!(m, feynkit_py::FeynkitModule);
   ```

   The host's existing `register_module!` creates `feynkit_native`, inserts it
   as `symbolica.community.feynkit_native`, and exposes an
   `initialize_module()` callback. Do not add `#[pymodule]` to `feynkit-py`.

4. Copy the checked-in `python/symbolica/community/feynkit` directory into the
   host's `python/symbolica/community` tree. Its `__init__.py` imports the native
   child and runs the initialization callback; `__init__.pyi` is the generated
   public stub.

Run the stub generator after changing the Python API:

```console
cargo run -p feynkit-py --features python_stubgen --bin feynkit-stubgen
```

## Model inspection and recomputation

`Model` exposes immutable `Parameter`, `Coupling`, `VertexRule`,
`LorentzStructure`, `Propagator`, `ModelFunction`, and `FormFactor` values
through both collections and name lookups. Process selectors similarly retain
whether they were specified by name or by PDG code through `ParticleSelector`.

Expression evaluation stays outside the model crate. Pass a Python callable to
`Model.recompute_with()` or as the optional `evaluator` argument to
`Model.with_parameter_card()`. The callable receives an `EvaluationRequest`
and returns `EvaluatedValues`; incomplete or unexpected results raise
`ModelError`, and neither the input model nor a partially updated model is
published on failure. Omitting `evaluator` from `with_parameter_card()` keeps
the invalidation-only behavior.

## Native particle-physics workflows and notebook figures

The common notebook workflows are native PyO3 methods, not a Python facade.
Each operation lives on the object that owns its configuration:

```python
import symbolica.community.feynkit as fk

result = model.generate_diagrams(["e-", "e+"], ["mu-", "mu+"], loops=1)
diagram = result.diagrams[0]
cff = diagram.build_cff()
jets = fk.JetDefinition.anti_kt(0.4).cluster(final_state_momenta)
loaded = fk.UfoLoader(restriction_name="massless").load("models/sm")
```

`FeynmanDiagram.to_linnest()` emits the complete Typst/Linnest source used for
figures. `to_svg()`, `to_html()`, `_repr_svg_()`, and `_repr_html_()` compile
that source with the Python package `typst>=0.15,<0.16`; Linnest and Kurvst are
embedded in the native module. The community host should declare `typst` as a
Python dependency so diagrams render automatically in Jupyter and Marimo.

## Installed import smoke test

This repository can validate registration in-process, but cannot perform an
installed import by itself because it intentionally does not ship the combined
wheel. After building and installing the host's `symbolica` wheel, run:

```console
python crates/feynkit-py/tests/installed_import_smoke.py
```

The host's release CI should run the same script after installing the wheel,
alongside its existing Spenso, Idenso, and Vakint import tests.
