"""Executable tour of every public Spenso type and its display protocols."""

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="full")


@app.function
def example_specification():
    # Each snippet is executed as shown, and can also be extracted for static checking.
    EXAMPLES = [
        (
            "Representation",
            "The index space: its dimension and which dual can contract with it. Built-ins include euc, mink, bis, cof, coad and cos.",
            "rep = sp.Representation.euc(2)",
            "rep",
        ),
        (
            "RepresentationName",
            "Dimension-independent representation identity, duality and canonical metric rule.",
            "rep_name = rep.name",
            "rep_name",
        ),
        (
            "Slot",
            "One labelled port in an index space. dual() changes its variance; representation exposes its space.",
            'slot = rep("mu")',
            "slot",
        ),
        (
            "TensorName",
            "A symbolic tensor head. Calling it with representations leaves ports open; slots label them immediately. print customizes its displayed name.",
            'A = sp.TensorName("spenso_api_tour::A", print={"typst": "A", "latex": "A"})',
            "A",
        ),
        (
            "TensorExpression",
            "An immutable symbolic tensor and ordered interface. Calling it assigns indices. Tensor-aware operations preserve that interface; to_expression() is the explicit boundary for general Symbolica algebra.",
            'expression = A(rep, rep)\nindexed = expression("mu", "nu")',
            "indexed",
        ),
        (
            "TensorRule",
            "A reusable tensor-safe replacement: construction validates wildcard closure, and each application checks the actual interface. Ordered rule lists share the same matching owner.",
            "rule = sp.TensorRule(indexed, 2 * indexed)\nreplaced = indexed.replace(rule)",
            "rule",
        ),
        (
            "TensorStructure",
            "Immutable metadata: optional name, scalar arguments, ordered slots, rank and shape. Unresolved ports retain their Representation.",
            "structure = indexed.structure",
            "structure",
        ),
        (
            "_AutoIndex",
            "AUTO leaves a port unresolved when other ports are indexed. The underscore alias denotes the same singleton.",
            'placeholder = sp.AUTO\npartially_indexed = expression(placeholder, "nu")',
            "placeholder",
        ),
        (
            "Tensor",
            "Actual dense or sparse component data. Integer indexing is flat row-major indexing; a coordinate sequence follows axes. expression() gives its symbolic descriptor; with_name gives stored data a library identity.",
            'tensor = sp.Tensor.dense(expression, [x, E("1"), E("0"), x + 1]).with_name("spenso_api_tour::M")',
            "tensor",
        ),
        (
            "TensorLibrary",
            "Resolves named tensors to stored data. Empty, constructed and HEP libraries have different contents; hep_lib_atom keeps exact symbolic entries.",
            'library = sp.TensorLibrary()\nlibrary.register(tensor)\nstored = library["spenso_api_tour::M"]\nreference = stored.expression()',
            "library",
        ),
        (
            "TensorNetwork",
            "The default display draws the current executable graph through Linnest. render and to_linnest accept graph.RenderSettings; expression() retains the source formula. execute evaluates the graph and result_tensor retrieves its data.",
            'network = tensor("i", "j") * tensor("j", "k")\nexecuted_network = tensor("i", "j") * tensor("j", "k")\nexecuted_network.execute()\nnetwork_result = executed_network.result_tensor()',
            "network",
        ),
        (
            "ExecutionStatus",
            "Immutable progress snapshot: pending graph operations, contractions and ready operations. step() returns an independent intermediate network.",
            "status = network.status\nnext_network = network.step()",
            "status",
        ),
        (
            "ReductionStatus",
            "Complete means no eligible work remains under the selected settings. Deferred and Capped retain an exact expression for subsequent work.",
            "contracted = indexed.contract(representations=[rep], collect_traces=False)\nreduction_status = contracted.reduction_status",
            "reduction_status",
        ),
        (
            "TensorEvaluator",
            "Optimizes symbolic components for repeated numerical evaluation. Each input row follows the params order; each output is a Tensor with the original interface.",
            "evaluator = tensor.evaluator([x], iterations=1, n_cores=1)\nevaluated = evaluator.evaluate([[2.0]])[0]",
            "evaluator",
        ),
        (
            "CompiledTensorEvaluator",
            "The same component evaluator compiled to a C++ shared library. The number_type selects real or complex evaluation; evaluate returns tensors.",
            'compiled = evaluator.compile("spenso_api_demo", str(build_directory / "tensor.cpp"), str(build_directory / "tensor.so"), number_type="complex", inline_asm="none", optimization_level=0)\ncompiled_result = compiled.evaluate([[2.0 + 0j]])[0]',
            "compiled",
        ),
        (
            "BroadcastFunction",
            "An elementwise function head. Symbolic input stays symbolic; concrete input creates a network that uses a registered numerical callback.",
            'broadcast = sp.BroadcastFunction("spenso_api_tour::square")',
            "broadcast",
        ),
        (
            "TensorFunctionLibrary",
            "Owns elementwise numerical callbacks used by network execution. It is separate from the registry of named tensor data.",
            "function_library = sp.TensorFunctionLibrary()\nfunction_library.register(broadcast, lambda value: value * value)",
            "function_library",
        ),
        (
            "PortPattern",
            "A structural port matcher, including wildcard representations and dimensions. It is a Symbolica pattern expression, not a concrete Slot.",
            'port_pattern = sp.PortPattern.exact(rep, S("spenso_api_tour::i_"))',
            "port_pattern",
        ),
        (
            "TensorPattern",
            "A tensor-aware replacement pattern. args matches ordinary scalar arguments; ports matches its tensor slots. any/vector and the HEP factories construct common patterns.",
            'pattern = sp.TensorPattern(A, ports=[port_pattern, sp.PortPattern.exact(rep, S("spenso_api_tour::j_"))])',
            "pattern",
        ),
        (
            "FactorProjector",
            "A normalized group of matrix factors for chain or trace; symbolic groups stay compact and component groups retain their data.",
            'factor_projector = sp.FactorProjector.symmetric(sp.TensorExpression.dirac_gamma(4)(sp.AUTO, sp.AUTO, "mu"), sp.TensorExpression.dirac_gamma(4)(sp.AUTO, sp.AUTO, "nu"))',
            "factor_projector",
        ),
        (
            "DisplaySettings",
            "Controls presentation: ports/schoonschip/call layouts, alphabet/graph/raw index labels, and superscript/array component coordinates. It does not change the tensor.",
            'display_settings = sp.DisplaySettings(component_style="array", show_dimensions=True)',
            "display_settings",
        ),
        (
            "ExecutionMode",
            "Network execution policy: All performs every possible contraction, Single selects one smallest-degree rewrite per step, Scalar only contracts scalar operations.",
            "execution_mode = sp.ExecutionMode.All",
            "execution_mode",
        ),
        (
            "SymbolicParallelism",
            "Selects how Symbolica work is scheduled through Rayon. set_symbolica_rayon_enabled sets this process-wide policy.",
            "parallelism = sp.SymbolicParallelism.Auto",
            "parallelism",
        ),
    ]

    for _name, _description in [
        ("CanonicalizationError", "A tensor network cannot be canonically relabelled."),
        ("CookingError", "A function or index payload cannot be cooked or restored."),
        (
            "DiracAdjointError",
            "A Dirac adjoint cannot be formed for the supplied tensor structure.",
        ),
        ("NetworkToolingError", "A symbolic network cannot be parsed or processed."),
    ]:
        EXAMPLES.append(
            (
                _name,
                _description
                + " Catch this specific exception, or its standard Python base class. The displayed instance illustrates its message.",
                f'error_{_name} = sp.{_name}("example diagnostic")',
                f"error_{_name}",
            )
        )

    PRELUDE = """from pathlib import Path
    from tempfile import TemporaryDirectory
    from symbolica import E, S
    from symbolica.community import tensor as sp

    _build = TemporaryDirectory(prefix="spenso-api-tour-")
    build_directory = Path(_build.name)
    x = S("spenso_api_tour::x")
    """
    from inspect import cleandoc

    return EXAMPLES, cleandoc(PRELUDE) + "\n"


@app.function
def run_examples():
    """Execute the displayed examples in one shared namespace."""
    EXAMPLES, PRELUDE = example_specification()
    namespace = {}
    exec(PRELUDE, namespace)  # noqa: S102 -- execute the owned examples shown below.
    examples = []
    for name, explanation, code, variable in EXAMPLES:
        exec(code, namespace)  # noqa: S102 -- code is an owned literal from the tour.
        examples.append((name, explanation, code, namespace[variable]))
    return examples, namespace


@app.function
def stub_catalog():
    """Include the complete declared member list beside each example."""
    import ast
    from pathlib import Path

    stub = Path(__file__).resolve().parents[2] / "docs/api/python/spynso3.pyi"
    source = stub.read_text()
    catalog = {}
    for node in ast.parse(source).body:
        if isinstance(node, ast.ClassDef):
            members = {}
            for member in node.body:
                if isinstance(member, ast.FunctionDef):
                    doc = ast.get_docstring(member) or ""
                    members.setdefault(member.name, doc.split("\n")[0])
                elif isinstance(member, ast.AnnAssign) and isinstance(
                    member.target, ast.Name
                ):
                    members[member.target.id] = ast.unparse(member.annotation)
            catalog[node.name] = members
    return catalog


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # Spenso: the complete Python API tour

    Start with **Representation → Slot → TensorName → TensorExpression**.
    Add data with **Tensor**, execute with **TensorNetwork**, and evaluate repeatedly
    with **TensorEvaluator**. The remaining entries cover registries, patterns,
    display and simplification settings, policies, and exceptions.

    Every example below has been executed. Expand an entry for its code, live
    display, text representation and complete declared member list. All examples
    share `from symbolica.community import tensor as sp` and the scalar `x = S("spenso_api_tour::x")`.
    """)
    return


@app.cell
def _(mo):
    examples, namespace = run_examples()
    catalog = stub_catalog()
    panels = {}
    for name, explanation, code, value in examples:
        member_rows = [
            {"Member": key, "Purpose": description}
            for key, description in catalog[name].items()
        ]
        panels[name] = mo.vstack(
            [
                mo.md(explanation),
                mo.md(f"```python\n{code}\n```"),
                mo.as_html(value),
                mo.md(f"**Python `repr`**\n\n```text\n{value!r}\n```"),
                mo.accordion(
                    {"All declared members": mo.ui.table(member_rows, selection=None)}
                ),
            ]
        )
    mo.accordion(panels)
    return (namespace,)


@app.cell
def _(mo, namespace):
    mo.vstack(
        [
            mo.md("## Results: symbolic composition, evaluated data, compiled data"),
            mo.as_html(namespace["network_result"]),
            mo.as_html(namespace["evaluated"]),
            mo.as_html(namespace["compiled_result"]),
            mo.md(
                "`network.expression()` keeps the source tensor expression; execution results live in `result_tensor()` or `result_scalar()`. The two numerical results above agree."
            ),
        ]
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Operations and simplification methods

    - `dot(p, q)`, `chain(start, end, *factors)` and `trace(rep, *factors)` construct tensor operations.
      Symbolic operands return `TensorExpression`; a concrete `Tensor` or `TensorNetwork` operand gives a `TensorNetwork`.
    - `outer`, `contract_ports`, `compose`, arithmetic and `trace` are methods on tensor expressions. `contract` performs structural contraction and returns a `TensorExpression`.
    - `replace` accepts Symbolica patterns and replacements, a reusable `TensorRule`, or an ordered rule list, preserving the tensor interface.
    - Materialization is explicit: `expand` emits a polynomial result and `evaluator` builds numerical evaluation.
    - Symbolic tensor operations return a `TensorExpression`. Its `to_expression()` crosses to ordinary Symbolica algebra. Reconstruct with `TensorExpression(raw, structure=source.structure)` to validate and retain a logical layout.
    - Dirac algebra: `simplify_algebra`, with `gamma=True, gamma_output="chains"` for collected words, gamma conjugation and `dirac_adjoint`.
    - Color algebra: `simplify_algebra`, with dimension-invariant substitution selected by `color_substitute_cof_dimension_invariants=True`.
    - Index algebra: `contract`, `simplify_algebra`, `with_lorentz_dimension`, cooking and canonicalization.
    - Compact notation: `to_dots`, `undo_dots`, `undo_chain`, and `undo_trace`. Shared `contract` handles indexed contractions and chain or trace collection.
    - `format_tensor`, `to_typst`, `to_html`, `to_svg` and `formatted` also accept ordinary Symbolica expressions.
      `TensorExpression(raw)` validates and wraps an ordinary expression.

    `TensorExpression` subclasses `Expression` and retains the shared typed value.
    Inspection, matching, conversions and evaluators are inherited; tensor-aware
    algebra and rewrite methods retain the tensor interface.
    `TensorPattern` remains a Symbolica pattern. Tensor gamma factories and index
    operations have tensor-specific meanings; a scalar gamma special function uses
    `Expression.gamma()`.
    """)
    return


if __name__ == "__main__":
    app.run()
