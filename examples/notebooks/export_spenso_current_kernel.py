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
from pathlib import Path

import numpy as np
import symbolica
from symbolica import E, PrintMode, Replacement, S, Symbol


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


def latex(expression):
    result = expression.format(mode=PrintMode.Latex, max_line_length=10000)
    for prefix, notation in [
        ("u", r"u_{%s}"),
        ("e", r"e^{%s}"),
        ("p", r"p^{%s}"),
        ("B", r"B_{%s}"),
    ]:
        result = re.sub(
            r"\b" + prefix + r"([0-3])\b",
            lambda m, notation=notation: notation % m[1],
            result,
        )
    return result


def aligned(rows):
    return "\\begin{align*}\n" + " \\\\\n".join(rows) + "\n\\end{align*}\n\n"


def renamed(expression, parameters, symbols):
    return expression.replace_multiple(
        [
            Replacement(parameter, symbol)
            for parameter, symbol in zip(parameters, symbols, strict=True)
        ]
    )


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


def component_outputs(namespace):
    """Derive compact output and prove it equals the contracted tensor."""
    u = list(S(*(f"spenso_paper::u{k}" for k in range(4))))
    e = list(S(*(f"spenso_paper::e{k}" for k in range(4))))
    p = list(S(*(f"spenso_paper::p{k}" for k in range(4))))
    b = list(S(*(f"spenso_paper::B{k}" for k in range(4))))
    parameters = namespace["parameters"]
    symbols = u + e + p
    vertex = [
        renamed(component / Symbol.I, parameters, symbols).expand()
        for component in namespace["join"].to_tensor()[:]
    ]
    for symbol in u:
        vertex = [component.collect(symbol) for component in vertex]
    i = Symbol.I
    # Weyl slash(p), independently checked against the contracted current below.
    slash = [
        [0, 0, p[0] - p[3], -p[1] + i * p[2]],
        [0, 0, -p[1] - i * p[2], p[0] + p[3]],
        [p[0] + p[3], p[1] - i * p[2], 0, 0],
        [p[1] + i * p[2], p[0] - p[3], 0, 0],
    ]
    denominator = p[0] ** 2 - sum(x**2 for x in p[1:])
    numerators = [
        -sum(b[row] * slash[row][col] for row in range(4)) for col in range(4)
    ]
    actual = [
        renamed(component, parameters, symbols) for component in namespace["Q24"][:]
    ]
    for component, numerator in zip(actual, numerators, strict=True):
        expanded = renamed(numerator, b, vertex) / denominator
        if (component - expanded).expand().cancel() != 0:
            raise ValueError("Displayed component does not equal the notebook output")
    inputs = [*namespace["quark"], *namespace["gluon"], *namespace["momentum"]]
    exact = [
        renamed(component, symbols, [exact_number(x) for x in inputs]).expand().cancel()
        for component in actual
    ]
    numeric = [complex(component.evaluate({})) for component in exact]
    np.testing.assert_allclose(namespace["value"][:], numeric, atol=1e-14, rtol=1e-14)
    vertex_latex = (
        "For the component output, write $u_a=Q(2)_a$, $e^\\mu=G(4)^\\mu$ "
        "and $p^\\mu=q_{24}^\\mu$. Defining $B_b=\\sum_{a,\\mu}u_a"
        "(\\gamma^\\mu)_{ab}e_\\mu$, the vertex contribution is $iB_b$, with\n"
        + aligned(
            [f"B_{{{k}}} &= {latex(component)}" for k, component in enumerate(vertex)]
        )
    )
    result_latex = (
        "In terms of the vertex components above, the output can be written as\n"
        + aligned(
            [
                rf"Q(2,4)_{{{k}}} &= -\frac{{{latex(-component)}}}{{s_{{24}}}}"
                for k, component in enumerate(numerators)
            ]
        )
        + r"where $s_{24}=(p^0)^2-(p^1)^2-(p^2)^2-(p^3)^2$."
        + "\n\n"
    )
    return vertex_latex, result_latex, exact


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
            "Executing the network and extracting the current components.",
        ),
        "value": ("evaluation", "Numerical evaluation of the current."),
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "assets").mkdir(exist_ok=True)
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
        module = ast.fix_missing_locations(ast.Module(body=body, type_ignores=[]))
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
        statements = body[:-1] if isinstance(body[-1], ast.Expr) else body
        lines = source.splitlines()[
            statements[0].lineno - 1 : statements[-1].end_lineno
        ]
        code = textwrap.dedent("\n".join(lines))
        blocks.append(
            (
                slug,
                "\\noindent\\begin{minipage}{\\linewidth}\n"
                + rf"\begin{{lstlisting}}[style=pythoncode,caption={{{caption}}},label={{lst:spenso-{slug}}}]"
                + "\n"
                + code
                + "\n\\end{lstlisting}\n\\end{minipage}\n\n",
            )
        )
        if name == "network":
            svg = namespace["network"].render()
            svg = svg.replace("<svg ", '<svg data-theme="light" ', 1)
            (output / "assets/spenso-current-network.svg").write_text(svg)
    vertex, result, exact = component_outputs(namespace)
    extra = {
        "inputs": r"\[Q(2)_a,\qquad G(4)^\mu,\qquad q_{24}^\mu.\]" + "\n\n",
        "vertex": vertex,
        "propagator": r"\[P_d=\frac{i}{s_{24}}\begin{pmatrix}"
        r"0&0&p^0-p^3&-p^1+ip^2\\0&0&-p^1-ip^2&p^0+p^3\\"
        r"p^0+p^3&p^1-ip^2&0&0\\p^1+ip^2&p^0-p^3&0&0"
        r"\end{pmatrix}.\]" + "\n\n",
        "current": r"\[Q(2,4)_c=-\frac{1}{s_{24}}\sum_b B_b(\not p)_{bc}.\]" + "\n\n",
        "network": r"\begin{figure}[H]"
        + "\n"
        + r"\centering"
        + "\n"
        + r"\includegraphics[width=.95\linewidth]{spenso-current-network.pdf}"
        + "\n"
        + r"\caption{The parsed current before execution. Tensor contractions and scalar operations, including the inverse propagator denominator, remain explicit.}"
        + "\n"
        + r"\label{fig:spenso-current-network}"
        + "\n"
        + r"\end{figure}"
        + "\n\n",
        "components": result,
    }
    # Verify the exact numerical equation already present in the source prose.
    expected = [(7 + 9 * Symbol.I) / 10, (5 - 3 * Symbol.I) / 10, E("0"), E("0")]
    if any((a - b).expand() != 0 for a, b in zip(exact, expected, strict=True)):
        raise ValueError(
            "Update the notebook's numerical-output equation for the new inputs"
        )
    section = "% Generated by export_spenso_current_kernel.py; edit the notebook, not this file.\n"
    for kind, content in blocks:
        section += (
            markdown_to_latex(content)
            if kind == "prose"
            else content + extra.get(kind, "")
        )
    (output / "sections/spenso-current-kernel.tex").write_text(section)
    subprocess.run(
        [
            rsvg_convert,
            "--format",
            "pdf",
            "--output",
            str(output / "assets/spenso-current-network.pdf"),
            str(output / "assets/spenso-current-network.svg"),
        ],
        check=True,
    )
    receipt = {
        "notebook": str(notebook),
        "notebook_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "symbolica_module": symbolica.__file__,
        "exact_output": [x.format(color_builtin_symbols=False) for x in exact],
        "checks": [
            "all four symbolic components equal the displayed abbreviated expressions",
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
