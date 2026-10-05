"""Export the executable current-kernel notebook for the pyAmpliCol paper.

Run with the Symbolica build used by the notebook and rsvg-convert on PATH:
    python examples/notebooks/export_spenso_current_kernel.py --output target/spenso-paper

The generated section uses the paper's software macros, pythoncode listing
style and cross-reference labels. It belongs in sections/ with the generated
network PDF in assets/. The notebook remains the prose and code source.
"""

import argparse
import ast
import hashlib
import json
import re
import subprocess
import textwrap
from fractions import Fraction
from functools import reduce
from operator import mul
from pathlib import Path

import numpy as np
import symbolica
import typst
from symbolica import AtomType, E, Replacement, Symbol
from symbolica.community.tensor import to_typst


def markdown_to_latex(markdown):
    """Convert the notebook's prose, preserving mathematics and paper references."""
    markdown = textwrap.dedent(markdown).strip()
    markdown = markdown.replace(
        "## Deriving current kernels with Spenso",
        r"\subsection{Deriving current kernels with \spenso{}}"
        "\n" + r"\label{sec:spenso}",
    )
    heading = ""
    if markdown.startswith(r"\subsection"):
        heading, markdown = markdown.split("\n\n", 1)
        heading += "\n\n"
    software = {
        "pyAmpliCol": "pyac",
        "Symbolica": "symbolica",
        "SymJIT": "symjit",
        "Spenso": "spenso",
        "Idenso": "idenso",
        "ALOHA": "aloha",
        "Comix": "comix",
    }
    references = {
        "equation (3.6)": r"eq.~\eqref{eq:zgg-currents}",
        "the amplitude expression following (3.6)": r"eq.~\eqref{eq:zgg-closure}",
        "the section on general contact\ndecomposition": r"section~\ref{sec:general-contact-decomposition}",
        "the section on general\ncontact decomposition": r"section~\ref{sec:general-contact-decomposition}",
        "The section on general contact\ndecomposition": r"Section~\ref{sec:general-contact-decomposition}",
        "the UFO-coverage table": r"table~\ref{tab:ufo-coverage}",
    }
    parts = re.split(r"(\$\$[\s\S]*?\$\$|\$[^$]*\$|`[^`]*`|\[@[^\]]+\])", markdown)
    result = []
    for part in parts:
        if part.startswith("$$"):
            result.append("\\[\n" + part[2:-2].strip() + "\n\\]")
        elif part.startswith("$"):
            result.append(part)
        elif part.startswith("`"):
            code = part[1:-1].replace("_", r"\_")
            result.append(r"\texttt{" + code + "}")
        elif part.startswith("[@"):
            result.append(r"\cite{" + part[2:-1] + "}")
        else:
            part = part.replace("–", "--").replace("’", "'")
            for name, macro in software.items():
                part = re.sub(
                    r"\b" + name + r"\b",
                    lambda _, macro=macro: "\\" + macro + "{}",
                    part,
                )
            for phrase, reference in references.items():
                part = part.replace(phrase, reference)
            result.append(part)
    return heading + "".join(result) + "\n\n"


def exact_number(value):
    """Recover the small rational inputs used in this pedagogical example."""
    value = complex(value)
    real = Fraction(value.real).limit_denominator(1000000)
    imag = Fraction(value.imag).limit_denominator(1000000)
    if abs(complex(real, imag) - value) > 1e-14:
        raise ValueError(
            "Input cannot be represented by the displayed rational notation"
        )
    return E(str(real)) + Symbol.I * E(str(imag))


def wrapped_component_svg(component):
    """Line-wrap the displayed component without introducing algebraic aliases."""
    factors = list(component)
    total = next(x for x in factors if x.get_type() == AtomType.Add)
    inverse = next(
        x for x in factors if x.get_type() == AtomType.Pow and list(x)[1] == -1
    )
    denominator = next(iter(inverse))
    coefficient = reduce(
        mul, (x for x in factors if x != total and x != inverse), E("1")
    )
    terms = list(total)
    if (coefficient * sum(terms) / denominator - component).expand().cancel() != 0:
        raise ValueError("Line-wrapped output differs from the original component")
    lines = [to_typst(term) for term in terms]
    rows = lines[0] + "".join(
        " \\ &" + ("" if row.startswith("-") else "+") + row for row in lines[1:]
    )
    sign = "-" if coefficient == -1 else to_typst(coefficient)
    source = "frac(" + sign + "lr((&" + rows + "))," + to_typst(denominator) + ")"
    document = (
        "#set page(width:auto,height:auto,margin:4pt)\n#set text(size:11pt)\n$ "
        + source
        + " $"
    )
    return typst.compile(document.encode(), format="svg").decode()


def numerical_check(namespace):
    inputs = [*namespace["quark"], *namespace["gluon"], *namespace["momentum"]]
    substitutions = [
        Replacement(p, exact_number(x))
        for p, x in zip(namespace["parameters"], inputs, strict=True)
    ]
    exact = [
        component.replace_multiple(substitutions).expand().cancel()
        for component in namespace["Q24"][:]
    ]
    numeric = [complex(component.evaluate({})) for component in exact]
    np.testing.assert_allclose(namespace["value"][:], numeric, atol=1e-14, rtol=1e-14)
    expected = [(7 + 9 * Symbol.I) / 10, (5 - 3 * Symbol.I) / 10, E("0"), E("0")]
    if any((a - b).expand() != 0 for a, b in zip(exact, expected, strict=True)):
        raise ValueError("Update the numerical-output equation for the new inputs")
    return exact


def output_block(values, slug, assets, rsvg_convert):
    """Capture the displayed value before later cells mutate any network."""
    graphics = []
    for index, value in enumerate(values):
        if slug == "network":
            svg = (
                value.render().to_svg().replace("<svg ", '<svg data-theme="light" ', 1)
            )
        elif slug == "components":
            svg = wrapped_component_svg(value)
        else:
            svg = value.to_svg()
        basename = f"{slug}-{index}"
        (assets / f"{basename}.svg").write_text(svg)
        subprocess.run(
            [
                rsvg_convert,
                "--format",
                "pdf",
                "--output",
                str(assets / f"{basename}.pdf"),
                str(assets / f"{basename}.svg"),
            ],
            check=True,
        )
        width = float(re.search(r'width="([0-9.]+)pt"', svg)[1])
        # Fit tensor expressions to the paper column without enlarging them.
        width = min(width, 395)
        if slug == "network":
            width *= 0.7  # Keep the diagram compact beside the other cell outputs.
        graphics.append(
            rf"\includegraphics[width={width:.3f}pt]{{spenso-current/{basename}.pdf}}"
        )
    note = " (line-wrapped)" if slug == "components" else ""
    return (
        "\\tcblower\n"
        + rf"{{\small\sffamily\color{{codecomment}} Output{note}}}\par\smallskip"
        + "\n\\centering\n"
        + "\\qquad\n".join(graphics)
        + "\n"
    )


def export(notebook, output, rsvg_convert):
    source = notebook.read_text()
    namespace = {}
    blocks = []
    listing_names = {
        "Q2": ("inputs", "Symbolic input currents and momentum."),
        "join": ("vertex", "The quark--gluon vertex contraction."),
        "propagator": ("propagator", "The massless quark propagator."),
        "current": ("current", "Applying the propagator to the vertex contribution."),
        "network": ("network", "Parsing the current into a tensor network."),
        "Q24": (
            "components",
            "Executing the network and displaying the first current component.",
        ),
        "value": ("evaluation", "Numerical evaluation of the current."),
    }
    output.mkdir(parents=True, exist_ok=True)
    assets = output / "assets/spenso-current"
    assets.mkdir(parents=True, exist_ok=True)
    (output / "sections").mkdir(exist_ok=True)
    for cell in ast.parse(source).body:
        if not isinstance(cell, ast.FunctionDef):
            continue
        if any(arg.arg == "mo" for arg in cell.args.args):
            blocks.append(("prose", cell.body[0].value.args[0].value))
            continue
        body = [
            statement
            for statement in cell.body
            if not isinstance(statement, ast.Return)
        ]
        if not body or isinstance(body[0], ast.Import):
            continue  # Hidden marimo setup.
        if not isinstance(body[-1], ast.Expr):
            raise TypeError("Each listing must end with the expression to display")
        module = ast.fix_missing_locations(ast.Module(body=body[:-1], type_ignores=[]))
        exec(compile(module, str(notebook), "exec"), namespace)  # noqa: S102 -- Execute the requested local notebook.
        assignments = {
            node.id
            for stmt in body
            for node in ast.walk(stmt)
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store)
        }
        names = assignments & listing_names.keys()
        if len(names) != 1:
            raise ValueError(f"Update the paper export for this changed cell: {names}")
        name = names.pop()
        slug, caption = listing_names[name]

        lines = source.splitlines()[body[0].lineno - 1 : body[-1].end_lineno]
        code = textwrap.dedent("\n".join(lines))
        # Capture the notebook's actual final expression.
        displayed = eval(
            compile(ast.Expression(body[-1].value), str(notebook), "eval"), namespace
        )
        values = list(displayed) if isinstance(displayed, tuple) else [displayed]
        listing = (
            rf"\begin{{notebookcell}}{{{caption}}}{{lst:spenso-{slug}}}"
            + "\n\\begin{lstlisting}[style=pythoncode]"
            + "\n"
            + code
            + "\n\\end{lstlisting}\n"
            + output_block(values, slug, assets, rsvg_convert)
            + "\\end{notebookcell}\n\n"
        )
        blocks.append((slug, listing))
    exact = numerical_check(namespace)
    section = "% Generated by export_spenso_current_kernel.py; edit the notebook, not this file.\n"
    for kind, content in blocks:
        section += markdown_to_latex(content) if kind == "prose" else content
    (output / "sections/spenso-current-kernel.tex").write_text(section.rstrip() + "\n")
    receipt = {
        "notebook": str(notebook),
        "notebook_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "symbolica_module": symbolica.__file__,
        "exact_output": [x.format(color_builtin_symbols=False) for x in exact],
        "checks": [
            "native component rendering preserved under line wrapping",
            "exact rational substitution agrees with the tensor evaluator",
        ],
    }
    (output / "validation.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Exported and verified {output / 'sections/spenso-current-kernel.tex'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--notebook",
        type=Path,
        default=Path(__file__).with_name("spenso_current_kernel.py"),
    )
    parser.add_argument("--output", type=Path, default=Path("target/spenso-paper"))
    parser.add_argument("--rsvg-convert", default="rsvg-convert")
    args = parser.parse_args()
    export(args.notebook, args.output, args.rsvg_convert)
