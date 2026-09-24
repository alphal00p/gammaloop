"""Executable tour of every public Spenso type and its display protocols."""

# ruff: noqa: PLR1711 -- marimo cells retain explicit return statements.

import marimo

__generated_with = "0.21.1"
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
            "Symbolic algebra plus an ordered external interface. Calling it assigns indices; repeated compatible indices contract. The simplification methods preserve and validate that interface.",
            'expression = A(rep, rep)\nindexed = expression("mu", "nu")',
            "indexed",
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
            "Actual dense or sparse component data. Integer indexing is flat row-major indexing; a coordinate sequence follows structure.slots. expression() gives its symbolic descriptor; with_name gives stored data a library identity.",
            'tensor = sp.Tensor.dense(expression, [x, E("1"), E("0"), x + 1]).with_name("spenso_api_tour::M")',
            "tensor",
        ),
        (
            "TensorLibrary",
            "Resolves named tensors to stored data. Empty, constructed and HEP libraries have different contents; hep_lib_atom keeps exact symbolic entries.",
            'library = sp.TensorLibrary()\nlibrary.register(tensor)\nreference = library["spenso_api_tour::M"]',
            "library",
        ),
        (
            "TensorNetwork",
            "A graph of operations with a semantic source expression. execute evaluates the graph; result_tensor/result_scalar retrieve the answer. to_dot exports the graph itself.",
            'network = tensor("i", "j") * tensor("j", "k")\nnetwork.execute()\nnetwork_result = network.result_tensor()',
            "network",
        ),
        (
            "TensorEvaluator",
            "Optimizes symbolic components for repeated numerical evaluation. Each input row follows the params order; each output is a Tensor with the original interface.",
            "evaluator = tensor.evaluator({}, {}, [x], iterations=1, n_cores=1)\nevaluated = evaluator.evaluate([[2.0]])[0]",
            "evaluator",
        ),
        (
            "CompiledTensorEvaluator",
            "The same component evaluator compiled to a C++ shared library. Complex evaluation is available for both real and complex inputs.",
            'compiled = evaluator.compile("spenso_api_demo", str(build_directory / "tensor.cpp"), str(build_directory / "tensor.so"), inline_asm="none", optimization_level=0)\ncompiled_result = compiled.evaluate_complex([[2.0 + 0j]])[0]',
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
            "DisplaySettings",
            "Controls presentation: ports/schoonschip/call layouts, alphabet/graph/raw index labels, and superscript/array component coordinates. It does not change the tensor.",
            'display_settings = sp.DisplaySettings(component_style="array", show_dimensions=True)',
            "display_settings",
        ),
        (
            "GammaSimplifySettings",
            "Controls gamma-chain ordering, closed traces and the optional four-dimensional three-gamma epsilon identity.",
            "gamma_settings = sp.GammaSimplifySettings.canonical()",
            "gamma_settings",
        ),
        (
            "ColorSimplifySettings",
            "Controls closed color traces, Fierz expansion across chains and substitution of fundamental-dimension invariants.",
            "color_settings = sp.ColorSimplifySettings(expand_cross_chain_fierz=False)",
            "color_settings",
        ),
        (
            "ColorCasimirSettings",
            "Controls conversion to a Casimir basis, including the fundamental dimension and Dynkin-index normalization.",
            "casimir_settings = sp.ColorCasimirSettings(substitute_fundamental_index=True)",
            "casimir_settings",
        ),
        (
            "CookTagFilter",
            "Selects function heads by tags: any, all, or the output tags of the associated CookSettings.",
            'tag_filter = sp.CookTagFilter.any(["spenso_api_tour::edge"])',
            "tag_filter",
        ),
        (
            "CookSourceFilter",
            "Selects whole function occurrences or only payloads nested inside representation indices.",
            "source_filter = sp.CookSourceFilter.representation_index_payload(tag_filter)",
            "source_filter",
        ),
        (
            "CookSettings",
            "Configures flattening or reversible encoding of function payloads. indices() is the preset for readable graph-based tensor indices; reversible() supports uncook.",
            "cook_settings = sp.CookSettings.indices()",
            "cook_settings",
        ),
        (
            "SchoonschipSettings",
            "Controls compact vector, dot and chain notation, recursion, sum expansion and contraction order. Expansion is opt-in.",
            "schoonschip_settings = sp.SchoonschipSettings.depth_first(depth_limit=2)",
            "schoonschip_settings",
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
        (
            "GammaChainOrdering",
            "RepeatedPairs moves matching gamma matrices together; Canonical orders the whole open chain.",
            "gamma_order = sp.GammaChainOrdering.Canonical",
            "gamma_order",
        ),
        (
            "CookMode",
            "FlattenedSymbol gives readable cooked names; ReversibleEncoding retains the information required by uncook.",
            "cook_mode = sp.CookMode.ReversibleEncoding",
            "cook_mode",
        ),
        (
            "SchoonschipMode",
            "Selects one simplification pass or recursive processing.",
            "schoonschip_mode = sp.SchoonschipMode.Recursive",
            "schoonschip_mode",
        ),
        (
            "SchoonschipTraversal",
            "Selects depth-first or breadth-first recursive processing.",
            "traversal = sp.SchoonschipTraversal.DepthFirst",
            "traversal",
        ),
        (
            "SchoonschipContractionOrder",
            "Chooses which contraction is performed first. This controls the execution strategy, not the mathematical expression.",
            "contraction_order = sp.SchoonschipContractionOrder.SmallestDegree",
            "contraction_order",
        ),
    ]

    for _name, _description in [
        ("CanonicalizationError", "A tensor network cannot be canonically relabelled."),
        ("CookingError", "A function or index payload cannot be cooked or restored."),
        (
            "DiracAdjointError",
            "A Dirac adjoint cannot be formed for the supplied tensor structure.",
        ),
        ("DotExpansionError", "Dot notation cannot be expanded consistently."),
        (
            "GammaConjugationError",
            "Conjugated gamma matrices cannot be rewritten consistently.",
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
    from symbolica.community import spenso as sp

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
    share `from symbolica.community import spenso as sp` and the scalar `x = S("spenso_api_tour::x")`.
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
    return examples, namespace


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
    - `outer`, `contract`, `compose`, arithmetic and `trace` are also methods on the relevant tensor types.
    - Scalar algebra: `expand`, `factor`, `collect`, `collect_num`, `collect_factors`, `together`, `cancel`, `apart` and Horner collection.
    - Dirac algebra: `simplify_gamma`, `collect_gamma_chains`, gamma conjugation and `dirac_adjoint`.
    - Color algebra: `simplify_color`, `collect_color`, `to_color_casimir` and dimension-invariant substitution.
    - Index algebra: `simplify_metrics`, `simplify_epsilon`, `with_lorentz_dimension`, cooking and canonicalization.
    - Compact notation: `schoonschip`, `schoonschip_net`, `to_dots`, chain collection and their expansion/undo methods.
    - `format_tensor`, `to_typst`, `to_html`, `to_svg` and `formatted` also accept ordinary Symbolica expressions.
      `as_tensor` restores a structured tensor after an ordinary expression operation.

    `TensorExpression` and `TensorPattern` inherit Symbolica algebra, but their gamma
    factories and tensor indexing have tensor-specific meanings. A scalar gamma
    special function remains `Expression.gamma()`.
    """)
    return


if __name__ == "__main__":
    app.run()
