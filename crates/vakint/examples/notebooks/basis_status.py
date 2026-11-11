"""State of the single-scale vacuum basis (1L-5L) for page 06.

Two parts, both data-driven:

* The status tables come from the benchmark generator itself
  (``crates/vakint/docs/benchmark/tools/{precision_tables,overview}.py``), the
  code that writes the precision tables of ``benchmark.pdf``. The page never
  retypes a number: it reads the generator's JSON (``overview.json`` and
  ``tables.json``) or runs the generator on the configured data tree.
* The live examples build their inputs from shipped master graphs
  (``data/basis-status``, SHA-256 checked) and evaluate them with Vakint on
  whatever data the engine loaded. Each order is classified from Vakint's own
  output: exact, exact in named constants (digits capped by those constants),
  numerical (digits Vakint vouches for), or failed closed.

Sources of the status tables, in this order:

1. ``VAKINT_GALLERY_STATUS``: an ``overview.json`` written by the generator
   (``tables.json`` next to it is read too), for example the PREVIEW edition
   of ``regen.sh`` that adds results not yet in a released tree.
2. ``VAKINT_GALLERY_ARTIFACTS``: the operator's data checkout. The generator
   runs on it in-process (about three seconds; cached by the tree's manifest
   digest). A tree the generator refuses is reported with its reason and the
   page falls back to 3.
3. The committed generator output ``docs/benchmark/data/{overview,tables}.json``
   (the released data the committed ``benchmark.pdf`` describes).

``VAKINT_GALLERY_BENCHMARK_DOCS`` names ``crates/vakint/docs/benchmark`` when
the notebooks are served from a copy outside the checkout.
"""

from __future__ import annotations

import contextlib
import hashlib
import html
import importlib
import io
import json
import os
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

DIRECTORY = Path(__file__).resolve().parent
INPUTS = DIRECTORY / "data" / "basis-status"
ORDERS = tuple(range(-5, 6))
CLASS_TEXT = {
    "S": "exact",
    "S*": "exact in named constants",
    "N": "numerical",
    "Z": "exact zero",
}
STATE_TEXT = {
    "below": "",
    "none": "not available",
    "partial": "not for every input",
}


def sup(power):
    return str(power).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def eps(power):
    return f"ε{sup(power)}"


# -- Where the generator and its output live ------------------------------------
def benchmark_docs():
    """``crates/vakint/docs/benchmark`` (generator tools + committed output), or None."""
    configured = os.environ.get("VAKINT_GALLERY_BENCHMARK_DOCS")
    candidates = [Path(configured)] if configured else []
    candidates.append(DIRECTORY.parents[1] / "docs" / "benchmark")
    for candidate in candidates:
        if (candidate / "tools" / "overview.py").is_file():
            return candidate.resolve()
    return None


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _generator():
    """The generator modules (precision_tables, overview), imported from the docs tools."""
    docs = benchmark_docs()
    if docs is None:
        raise FileNotFoundError(
            "the benchmark generator (crates/vakint/docs/benchmark/tools) is not "
            "next to these notebooks; set VAKINT_GALLERY_BENCHMARK_DOCS"
        )
    tools = str(docs / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return (
        importlib.import_module("precision_tables"),
        importlib.import_module("overview"),
        docs,
    )


def _cache_directory():
    configured = os.environ.get("VAKINT_GALLERY_WORK")
    base = Path(configured) if configured else Path.home() / ".cache" / "vakint-gallery"
    directory = base / "basis-status"
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def generate(root):
    """Run the benchmark generator on the data tree ``root`` (RELEASED mode):
    (overview, tables). Cached by the tree's root manifest and the generator's
    sources. Raises ValueError with the generator's own reason when it refuses
    the tree."""
    pt, ov, docs = _generator()
    root = Path(root).resolve()
    key = hashlib.sha256(
        "|".join(
            [
                _sha256(root / "manifest.json"),
                _sha256(docs / "tools" / "precision_tables.py"),
                _sha256(docs / "tools" / "overview.py"),
                *(
                    _sha256(path)
                    for path in sorted((docs / "data").glob("*.json"))
                    if path.name
                    in {
                        "pole-depth.json",
                        "r-operation-poles.json",
                        "constant-names.json",
                    }
                ),
            ]
        ).encode()
    ).hexdigest()[:24]
    cache = _cache_directory()
    tables_path = cache / f"tables-{key}.json"
    overview_path = cache / f"overview-{key}.json"
    if tables_path.is_file() and overview_path.is_file():
        overview = json.loads(overview_path.read_text())
        overview["cached"] = True
        return overview, json.loads(tables_path.read_text())
    src = docs.parents[1] / "src"
    names = docs / "data" / "constant-names.json"
    # The generator reports refusals as SystemExit(reason) and prints progress.
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            builtins = pt.read_builtins(src / "inspection" / "constants.rs")
            builtins["route_sources"] = pt.read_route_sources(src)
            tables = pt.build(
                root,
                builtins,
                docs / "data" / "pole-depth.json",
                "provisional",
                None,
                docs / "data" / "r-operation-poles.json",
            )
            tables_path.write_text(json.dumps(tables, ensure_ascii=False))
            overview = ov.build(
                SimpleNamespace(
                    data_root=root,
                    tables=tables_path,
                    overlay=None,
                    names=names if names.is_file() else None,
                    vakint_src=src,
                    pole_depth=docs / "data" / "pole-depth.json",
                )
            )
    except SystemExit as refusal:
        tables_path.unlink(missing_ok=True)
        raise ValueError(
            str(refusal.code if refusal.code is not None else refusal)
        ) from None
    overview_path.write_text(json.dumps(overview, ensure_ascii=False))
    return overview, tables


def load_status(checkout=None):
    """The status tables and where they come from: a dict with ``overview``,
    ``tables`` (may be None), ``origin`` (file | generated | committed),
    ``source`` (a short description) and ``notes`` (messages for the reader,
    e.g. why a configured tree could not be tabulated)."""
    notes = []
    configured = os.environ.get("VAKINT_GALLERY_STATUS")
    if configured:
        path = Path(configured)
        try:
            overview = json.loads(path.read_text())
            tables_path = path.with_name("tables.json")
            tables = (
                json.loads(tables_path.read_text()) if tables_path.is_file() else None
            )
            return {
                "overview": overview,
                "tables": tables,
                "origin": "file",
                "source": f"generator output {path.parent.name}/{path.name}",
                "notes": notes,
            }
        except (OSError, ValueError) as error:
            notes.append(
                f"VAKINT_GALLERY_STATUS could not be read ({error}); using the next source."
            )
    if checkout is not None:
        started = time.perf_counter()
        try:
            overview, tables = generate(checkout)
            return {
                "overview": overview,
                "tables": tables,
                "origin": "generated",
                "source": (
                    "the benchmark generator run on the configured data tree "
                    + (
                        "(cached)"
                        if overview.get("cached")
                        else f"({time.perf_counter() - started:.1f} s)"
                    )
                ),
                "notes": notes,
            }
        except (OSError, ValueError, KeyError) as error:
            notes.append(
                "The benchmark generator does not tabulate the configured data tree "
                f"{Path(checkout).name}: {error}. Its tables below are those of the "
                "committed benchmark (released data); the live examples further down "
                "still evaluate on the configured tree."
            )
    docs = benchmark_docs()
    if docs is None:
        notes.append(
            "No benchmark generator output is available (set VAKINT_GALLERY_BENCHMARK_DOCS)."
        )
        return {
            "overview": None,
            "tables": None,
            "origin": None,
            "source": None,
            "notes": notes,
        }
    overview = json.loads((docs / "data" / "overview.json").read_text())
    tables_path = docs / "data" / "tables.json"
    tables = json.loads(tables_path.read_text()) if tables_path.is_file() else None
    return {
        "overview": overview,
        "tables": tables,
        "origin": "committed",
        "source": "the committed benchmark tables (docs/benchmark/data/overview.json)",
        "notes": notes,
    }


# -- Status tables as plain rows --------------------------------------------------
def class_label(cls, digits=None):
    """S, S*(X), N(X), 0 as short text."""
    if cls in ("S*", "N") and digits is not None:
        return f"{cls}({digits})"
    return {"Z": "0"}.get(cls, cls)


def result_cell(cell):
    """(short label, kind) of one cell of the overview's result table."""
    state = cell.get("state")
    if state == "ok":
        label = class_label(cell.get("class"), cell.get("X"))
        if cell.get("marks"):
            label += "".join(f"[{m}]" for m in cell["marks"])
        return label, cell.get("class")
    if state == "partial":
        return "not every input" + "".join(
            f"[{m}]" for m in cell.get("marks") or []
        ), "partial"
    if state == "none":
        return "—", "none"
    return "", "below"


def result_rows(overview):
    """The per-loop, per-order table (the benchmark's table 2 / overview (c))."""
    orders = overview["orders2"]
    rows = []
    for row in overview["result"]:
        entry = {
            "loops": row["label"]
            + "".join(f"[{m}]" for m in row.get("marks") or [])
            + (" (preview)" if row["kind"] != "released" else "")
        }
        for order in orders:
            entry[eps(order)] = result_cell(row["cells"][str(order)])[0]
        rows.append(entry)
    for former in overview.get("former") or []:
        entry = {"loops": former["label"]}
        for order in orders:
            cell = former["cells"][str(order)]
            state = cell.get("state")
            if state == "exact":
                value = "0" if cell.get("value") == "zero" else "S"
                if cell.get("sstar_X"):
                    value = f"S*({cell['sstar_X']})"
                entry[eps(order)] = value
            elif state == "numerical":
                entry[eps(order)] = cell.get("v02x") or "N"
            elif state == "fails-closed":
                entry[eps(order)] = "fails closed"
            elif state == "none":
                entry[eps(order)] = "—"
            else:
                entry[eps(order)] = ""
        rows.append(entry)
    columns = ["loops"] + [
        eps(o) for o in orders if any(r.get(eps(o)) not in ("", None) for r in rows)
    ]
    return rows, columns


def result_notes(tables, rows):
    """Footnotes [a]... used in ``rows``, from the generator's tables.json."""
    if not tables:
        return []
    used = sorted(
        {
            m
            for r in rows
            for v in r.values()
            for m in re.findall(r"\[([a-z])\]", str(v))
        }
    )
    notes = tables.get("notes_table2") or []
    if isinstance(notes, list):  # note i is mark chr(ord("a") + i), as in the PDF
        notes = {chr(ord("a") + i): text for i, text in enumerate(notes)}
    return [(m, notes.get(m)) for m in used if notes.get(m)]


def summary_rows(overview):
    """Coefficient counts per loop and order (overview (a)): S / S* / N / 0."""
    rows = []
    for item in overview["summary"]:
        for order, x in item["orders"].items():
            pending = sum((x.get("pending") or {}).values())
            row = {
                "loops": item["label"],
                "order": eps(int(order)),
                "masters holding it": x["held"],
                "exact": x["S"],
                "exact in named constants": x["S*"],
                "numerical": x["N"],
                "exact zero": x["Z"],
                "S* digits": "–".join(str(v) for v in x["sstar_digits"])
                if x.get("sstar_digits")
                else "",
                "N digits (min / median)": (
                    f"{x['n_digits_min']} / {x['n_digits_median']:g}"
                    if x.get("n_digits_min") is not None
                    else ""
                ),
            }
            if item["kind"] != "released":
                row["of which pending"] = pending
            rows.append(row)
    return rows


def master_rows(overview, loop=5, status=None, order=None):
    """Per-master rows at ``loop`` (overview (b)); ``status`` filters on the class
    at ``order`` (any order when ``order`` is None): exact | named | numerical |
    zero | missing."""
    want = {
        "exact": {"S"},
        "named": {"S*"},
        "numerical": {"N"},
        "zero": {"Z"},
        "missing": {None},
        None: None,
    }[status]
    orders = [o for o in overview["orders1"] if -loop <= o <= 1]
    rows = []
    for r in overview["table1"]:
        if r["loop"] != loop:
            continue
        cells = {}
        classes = {}
        for o in orders:
            c = r["cells"].get(str(o))
            if c is None:
                cells[eps(o)], classes[o] = "—", None
                continue
            released = c.get("class")
            pending = c.get("pending")
            if pending and (released is None or _rank(pending) <= _rank(c)):
                cells[eps(o)] = (
                    class_label(pending["class"], pending.get("X")) + " (pending)"
                )
                classes[o] = pending["class"]
            elif released:
                cells[eps(o)] = class_label(released, c.get("X"))
                classes[o] = released
            else:
                cells[eps(o)], classes[o] = "—", None
        if want is not None:
            look = [order] if order is not None else orders
            if not any(classes.get(o) in want for o in look):
                continue
        name = r["name"]
        if r.get("design_name"):
            name += f" ({r['design_name']})"
        rows.append({"master": name, **cells})
    return rows, ["master", *[eps(o) for o in orders]]


def _rank(c):
    return (
        {"S": 0, "Z": 0, "S*": 1, "N": 2}.get(c.get("class"), 3),
        -(c.get("X") or 0),
    )


def constants_rows(overview):
    """Named constants (overview (d)): definition, digits, traditional form."""
    rows = []
    for c in overview.get("constants") or []:
        rows.append(
            {
                "constant": c["code"],
                "definition": c.get("definition_full") or c.get("definition") or "",
                "digits": "any precision"
                if c.get("arbitrary_precision")
                else (c.get("digits") or ""),
                "traditional form": c.get("traditional") or "",
                "status of the traditional form": c.get("traditional_status") or "",
                "used in": ", ".join(
                    f"{k} {v}" for k, v in (c.get("used_in") or {}).items()
                ),
            }
        )
    return rows


def method_paragraphs(overview):
    return [(m["title"], m["text"]) for m in overview.get("method") or []]


def tree_label(overview):
    """The data tree the tables describe, e.g. "vakint-external-data v0.2.2 (0faa1c795664)"."""
    tree = (overview or {}).get("tree") or {}
    if tree.get("label"):
        return tree["label"]
    commit = ((tree.get("git") or {}).get("commit") or "")[:12]
    return f"vakint-external-data {tree.get('version')}" + (
        f" ({commit})" if commit else ""
    )


def five_loop_basis(overview):
    """The 5L basis the tables describe ("BASIS-v2", "BASIS-v3" from a preview overlay, ...)."""
    basis = ((overview.get("overlay") or {}).get("basis") or {}).get("name")
    return basis or (overview.get("tree") or {}).get("five_loop_basis")


# -- Live examples ----------------------------------------------------------------
def input_manifest():
    return json.loads((INPUTS / "manifest.json").read_text())


EXAMPLES = {
    "4L-spurious": {
        "title": "4L: a dotted input whose reduction has spurious poles",
        "graph": "4L-P15",
        "dots": [0],
        "loops": 4,
        "point": (
            "The reduction coefficients of this input carry 1/ε poles on masters "
            "whose finite parts hold four-loop Bessel-moment constants "
            "(K5Moment, PR7D2, ...). By the R-operation the 1/ε pole of any "
            "four-loop vacuum integral needs only three-loop finite parts, so "
            "those constants must cancel there: Vakint returns ε⁻¹ exact (no "
            "named constant), and ε⁰ keeps them."
        ),
    },
    "5L-master": {
        "title": "5L master P1 (a numerator master)",
        "graph": "5L-P1",
        "dots": [],
        "loops": 5,
        "point": (
            "On BASIS-v3 with the adopted pole forms: ε⁻⁵..ε⁻² exact, ε⁻¹ exact "
            "in named four-loop constants, ε⁰ numerical."
        ),
    },
    "5L-dotted": {
        "title": "5L input: master P1 with one dot",
        "graph": "5L-P1",
        "dots": [0],
        "loops": 5,
        "heavy": "23-32 minutes on 4 cores (measured on the v0.3.0 RC)",
        "point": (
            "A non-master input is reduced first. On BASIS-v3 its reduction "
            "coefficients have no pole at d = 4, so ε⁻⁵..ε⁻² stay exact, ε⁻¹ "
            "exact in named constants and ε⁰ numerical."
        ),
    },
    "5L-P108": {
        "title": "5L former master P108 as an input",
        "graph": "5L-P108",
        "dots": [],
        "loops": 5,
        "needs": "BASIS-v3",
        "heavy": (
            "long: its power difference (11) lies outside the recorded 5L entry "
            "predicates (9 and 10), so RustRed reduces it without a recorded "
            "route; unfinished after 2 h 15 min on the v0.3.0 RC"
        ),
        "point": (
            "P108 was a master of BASIS-v2. On BASIS-v3 it reduces to the new "
            "basis, and its double pole comes back exact: ε⁻² = −ζ3/2."
        ),
    },
    "5L-T1-8175": {
        "title": "5L T1-8175 (master P135 of BASIS-v3)",
        "graph": "5L-P135",
        "dots": [],
        "loops": 5,
        "needs": "BASIS-v3",
        "point": (
            "Its double pole is ε⁻² = 4 ζ5 (part A of the R-operation 5 ζ5, "
            "overall counterterm −ζ5, derived by R*)."
        ),
    },
}


def example_graph(key):
    """(dot source, model JSON, manifest record) of a shipped example input,
    SHA-256 checked."""
    spec = EXAMPLES[key]
    manifest = input_manifest()
    record = manifest["graphs"][spec["graph"]]
    dot_path = INPUTS / record["path"]
    model_path = INPUTS / manifest["models"][record["model"]]["path"]
    for path, expected in (
        (dot_path, record["sha256"]),
        (model_path, manifest["models"][record["model"]]["sha256"]),
    ):
        if _sha256(path) != expected:
            raise ValueError(f"{path.name}: checksum mismatch")
    return dot_path.read_text(), model_path.read_text(), record


def example_input(key):
    """The Vakint integral of an example: the master graph (with its numerator
    and powers; the example's dots added), built with HEPKit and
    ``integral_from_diagram`` as Vakint's own acceptance tests do."""
    from symbolica import E
    from symbolica.community.hepkit import FeynmanDiagram, Kinematics, Model
    from symbolica.community.hepkit.vakint import integral_from_diagram

    dot, model_json, _record = example_graph(key)
    diagram = FeynmanDiagram.from_dot(Model.from_json(model_json), dot)
    family = diagram.propagator_family(kinematics=Kinematics(E("gammalooprs::dim")))
    powers = {}
    for line in dot.splitlines():
        edge = re.search(r"\[id=(\d+),", line)
        if "->" in line and edge:
            powers[int(edge.group(1))] = int(re.search(r"power=(-?\d+)", line).group(1))
    for edge in EXAMPLES[key]["dots"]:
        powers[edge] += 1
    integral = integral_from_diagram(
        diagram,
        family,
        E(diagram.numerator),
        powers=[powers[k] for k in sorted(powers)],
    )
    return integral.to_expression(), powers


_engines = {}


def example_engine(
    gallery, loops, *, precision, mode="symbolic", terms=None, progress=None
):
    """A FMFTandMATAD engine (the data's own normalization, so master
    coefficients read as in the data) with the ``loops``-loop data loaded.
    At most two engines stay alive (a 5L load holds over 10 GB)."""
    key = (loops, int(precision), mode, terms or loops + 1)
    if key in _engines:
        return _engines[key]
    while len(_engines) >= 2:
        _engines.pop(next(iter(_engines)))
    engine = gallery.engine(
        precision=int(precision),
        mode=mode,
        normalization="FMFTandMATAD",
        use_dot_product_notation=True,
        number_of_terms_in_epsilon_expansion=terms or loops + 1,
    )
    report = gallery.load_data(engine, loops, progress=progress)
    _engines[key] = (engine, report)
    return engine, report


_PLACEHOLDER = re.compile(r"RustRedCoefficient\(")
_NAMED = re.compile(
    r"RustRedConstant\([^,]+,\s*(?:[A-Za-z_]+::)?(?:\{\}::)?([A-Za-z0-9_]+(?:\([0-9, ]+\))?)\)"
)
_FLOAT = re.compile(r"(?<![A-Za-z0-9_^])\d+\.\d+(?:e[-+]?\d+)?")


def _unit(value):
    from symbolica import E

    for symbol in ("UFO::mUV", "vakint::mursq", "vakint::muvsq"):
        value = value.replace(E(symbol), E("1"))
    return value.expand()


def classify(text):
    """Class of one exact-output coefficient: (class, named constants)."""
    named = sorted(set(_NAMED.findall(text)))
    if _PLACEHOLDER.search(text) or _FLOAT.search(text):
        return "N", named
    if named:
        return "S*", named
    if text.strip() in ("0", "-0"):
        return "Z", named
    return "S", named


def readable(text):
    """Vakint's canonical text with namespaces and data tags dropped."""
    text = re.sub(r"RustRedConstant\([^,]+,\s*", "RRC(", text)
    text = re.sub(r"(?:vakint|symbolica|UFO)::(?:\{\}::)?", "", text)
    return text.replace("\x1b[35m", "").replace("\x1b[0m", "")


def evaluate_example(gallery, key, *, precision=60, progress=None):
    """Evaluate one example with symbolic masters and report per ε order:
    class, digits, value text, the named constants, and the exact expression
    (a Symbolica Expression, for the pager). Orders Vakint refuses are reported
    as failed closed with Vakint's reason. Plain data except ``expression``."""
    from symbolica import E
    from symbolica.community.hepkit.vakint import MissingMasterOrdersError

    spec = EXAMPLES[key]
    loops = spec["loops"]
    started = time.perf_counter()
    integral, powers = example_input(key)
    engine, report = example_engine(
        gallery, loops, precision=precision, progress=progress
    )
    loaded = time.perf_counter()
    out = {
        "key": key,
        "title": spec["title"],
        "loops": loops,
        "powers": powers,
        "data": (report or {}).get("version"),
        "precision": int(precision),
        "orders": [],
        "missing": None,
    }
    try:
        value = engine.evaluate(
            integral, **gallery.progress_options(engine.evaluate, progress)
        )
    except MissingMasterOrdersError as error:
        out["missing"] = gallery.missing_orders(error)
        out["seconds"] = {
            "load": loaded - started,
            "evaluate": time.perf_counter() - loaded,
        }
        return out
    evaluated = time.perf_counter()
    value = _unit(value)
    epsilon = E("vakint::ε")
    exact = {}
    for power, coefficient in value.coefficient_list(epsilon):
        exact[epsilon_order(power)] = coefficient
    numerical, _ = engine.numerical_evaluation(value, {})
    rows = {}
    for order, coefficient in sorted(numerical.coefficients.items()):
        rows[int(order)] = gallery.laurent_row(order, coefficient)
    # Orders absent from Vakint's output are exact zeros (from the most singular
    # one the loop count allows).
    for order in sorted(set(range(-loops, max([*exact, *rows], default=0) + 1))):
        text = exact[order].to_canonical_string() if order in exact else "0"
        cls, named = classify(text)
        cap = constant_digits(engine, named) if cls == "S*" else {}
        row = rows.get(order)
        digits = row["real"]["digits"] if row else None
        out["orders"].append(
            {
                "order": order,
                "class": cls,
                "named": named,
                "cap": min((v for v in cap.values() if v), default=None),
                "digits": digits,
                "exact_zero": bool(row and row["real"]["exact_zero"]),
                "value": row["real"]["text"] if row else "0",
                "exact": (
                    readable(str(exact[order]))
                    if cls in ("S", "Z", "S*") and order in exact and len(text) < 4000
                    else None
                ),
                "expression": exact.get(order),
            }
        )
    out["seconds"] = {
        "load": loaded - started,
        "evaluate": evaluated - loaded,
        "numerical": time.perf_counter() - evaluated,
    }
    return out


def constant_digits(engine, names):
    """Stored digits of named finite constants (the S* cap), by name."""
    wanted = set(names)
    out = {}
    for symbol in engine.constant_symbols():
        name = str(symbol.expression).rsplit("::", 1)[-1]
        if name not in wanted:
            continue
        try:
            precision = engine.get_constant(symbol, prec=2).available_precision
            out[name] = getattr(precision, "max_decimal_precision", None)
        except (ValueError, TypeError):
            out[name] = None
    return out


def retain_poles(gallery, key, *, progress=None):
    """Spurious poles of an example's reduction, from Vakint's Retain output
    (masters kept as symbols): {master: lowest ε power of its reduction
    coefficient} for the masters whose coefficient has a pole, and the number
    of masters."""
    from symbolica import E

    spec = EXAMPLES[key]
    integral, _ = example_input(key)
    engine, _ = example_engine(gallery, spec["loops"], precision=32, mode="retain")
    value = _unit(
        engine.evaluate(integral, **gallery.progress_options(engine.evaluate, progress))
    )
    epsilon = E("vakint::ε")
    masters = sorted(set(_MASTER.findall(value.to_canonical_string())))
    poles = {}
    for master in masters:
        coefficient = value.coefficient(E(master))
        lowest = coefficient.series(epsilon, 0, 0).get_trailing_exponent()
        if isinstance(lowest, tuple):
            lowest = lowest[0] / lowest[1]
        if lowest < 0:
            poles[readable(master)] = int(lowest)
    return poles, len(masters)


_MASTER = re.compile(r"vakint::(?:\{\}::)?RustRedMaster\([^()]*\)")


def epsilon_order(key):
    """The ε power of a ``coefficient_list`` key (1, ε, ε^2, 1/ε, 1/ε^2, ε^(-2), ...)."""
    text = re.sub(r"(?:vakint::)?(?:\{\}::)?", "", key.to_canonical_string()).replace(
        " ", ""
    )
    if text == "1":
        return 0
    match = re.fullmatch(r"(1/)?ε(?:\^\(?(-?\d+)\)?)?", text)
    if not match:
        raise ValueError(f"not a power of ε: {text}")
    power = int(match.group(2) or 1)
    return -power if match.group(1) else power


def class_badge(cls):
    return {
        "S": "exact",
        "Z": "exact 0",
        "S*": "exact, named constants",
        "N": "numerical",
    }.get(cls, cls)


def example_table_rows(result):
    """Rows of one evaluated example: order, class, digits, value, exact form."""
    rows = []
    for o in result["orders"]:
        if o["class"] in ("S", "Z"):
            digits = "exact"
        elif o["class"] == "S*":
            digits = f"{o['digits']} here; the exact form is good to {o['cap']} digits"
            digits += f" ({', '.join(o['named'])})"
        else:
            digits = o["digits"]
        form = html.unescape(o["exact"] or "")
        rows.append(
            {
                "order": eps(o["order"]),
                "what Vakint returns": class_badge(o["class"]),
                "value": "0" if o["class"] == "Z" else _short(o["value"]),
                "digits vouched for": digits,
                "exact form": form if len(form) <= SHORT else form[:SHORT] + " …",
            }
        )
    return rows


SHORT = 160


def _short(text, digits=24):
    text = str(text)
    m = re.match(rf"^(-?\d+\.\d{{{digits}}})\d+(e[-+]?\d+)?$", text)
    return f"{m.group(1)}…{m.group(2) or ''}" if m else text


# -- What the loaded data holds -----------------------------------------------------
def data_basis(root):
    """("BASIS-v3" | "BASIS-v2" | None, version) of a data tree's 5L pack: BASIS-v3
    holds T1-8175 as master P135 and no longer P108 (a former master)."""
    if root is None:
        return None, None
    root = Path(root)
    try:
        version = json.loads((root / "manifest.json").read_text()).get("version")
        five = json.loads((root / "5L" / "manifest.json").read_text())
    except (OSError, ValueError):
        return None, locals().get("version")
    names = {m.get("name") for m in five.get("masters") or []}
    if "P135" in names and "P108" not in names:
        return "BASIS-v3", version
    if "P108" in names:
        return "BASIS-v2", version
    return None, version


def master_names(root, loops):
    """{RustRedMaster text: catalogue name} for the ``loops``-loop masters of a
    data tree (master_id = L<loops>_<family>_<signed powers>)."""
    out = {}
    if root is None:
        return out
    try:
        manifest = json.loads((Path(root) / f"{loops}L" / "manifest.json").read_text())
    except (OSError, ValueError):
        return out
    for master in manifest.get("masters") or []:
        powers = ",".join(str(p) for p in master.get("powers") or [])
        out[f"RustRedMaster(rustred_family_{master.get('family_sha256')},{powers})"] = (
            master.get("name")
        )
    return out


def _run_setting():
    return os.environ.get("VAKINT_GALLERY_RUN_EXAMPLES", "").strip()


def autorun():
    """Evaluate the live examples without a button press (headless runs and
    static exports): VAKINT_GALLERY_RUN_EXAMPLES = 1 (the quick ones), all (the
    heavy five-loop ones too) or a comma-separated list of example keys."""
    return _run_setting() not in ("", "0")


def default_examples():
    """The five-loop examples selected by default: the quick ones, or as
    VAKINT_GALLERY_RUN_EXAMPLES says (all | key,key,...)."""
    setting = _run_setting()
    five = [key for key, spec in EXAMPLES.items() if spec["loops"] == 5]
    if setting == "all":
        return five
    chosen = [key for key in setting.split(",") if key in five]
    return chosen or [key for key in five if not EXAMPLES[key].get("heavy")]
