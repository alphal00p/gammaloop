use std::ffi::CString;

use pyo3::{
    prelude::*,
    types::{PyCFunction, PyDict, PyList, PyModule},
};
use spenso::{
    portable_payload::{
        MATH_DISPLAY_ATTACHMENT_SCHEMA, MathDisplayDeclaration, MathDisplayDeclarations,
        PortableRepresentationClass, REPRESENTATION_ATTACHMENT_SCHEMA, RepresentationDeclaration,
        RepresentationDeclarations, attachments_for_atom, canonical_math_display_symbol_name,
        register_math_display_symbol,
    },
    structure::representation::{IndexDisplay, IndexPalette, IndexRow, LibraryRep},
};
use spynso3::SpensoModule;
use symbolica::{
    api::python::{SymbolicaCommunityModule, create_symbolica_module},
    atom::{Atom, FunctionBuilder, NamespacedSymbol, SymbolBuilder},
    parse,
};
use tymbolica_atom_payload::{encode_atom_from_set, parse_payload};

const SPENSO_WRAPPER: &str = r#"
from ..spenso_native import *

initialize_module()
"#;

fn install_package<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyModule>> {
    let package = PyModule::new(py, name)?;
    package.setattr("__package__", name)?;
    package.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(name, &package)?;
    Ok(package)
}

fn register_spenso<'py>(core: &Bound<'py, PyModule>) -> PyResult<Bound<'py, PyModule>> {
    let py = core.py();
    let native = PyModule::new(py, "spenso_native")?;
    let initialize_module =
        PyCFunction::new_closure(py, Some(c"initialize_module"), None, |args, _kwargs| {
            SpensoModule::initialize(args.py())
        })?;
    native.add("initialize_module", initialize_module)?;
    SpensoModule::register_module(&native)?;
    core.add_submodule(&native)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("symbolica.community.spenso_native", &native)?;
    Ok(native)
}

fn import_spenso_wrapper<'py>(
    py: Python<'py>,
    community: &Bound<'py, PyModule>,
) -> PyResult<Bound<'py, PyModule>> {
    let name = "symbolica.community.spenso";
    let wrapper = PyModule::new(py, name)?;
    wrapper.setattr("__package__", name)?;
    wrapper.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(name, &wrapper)?;
    let source = CString::new(SPENSO_WRAPPER).expect("Spenso wrapper contains a null byte");
    py.run(&source, Some(&wrapper.dict()), Some(&wrapper.dict()))?;
    community.add("spenso", &wrapper)?;
    PyModule::import(py, name)
}

#[test]
fn combined_core_exposes_rich_display_with_one_expression_type() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;

        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_spenso(&core)?;
        let spenso = import_spenso_wrapper(py, &community)?;

        let locals = PyDict::new(py);
        locals.set_item("core", &core)?;
        locals.set_item("spenso", &spenso)?;
        let assertions = CString::new(
            r##"
for name in (
    "DisplaySettings",
    "format_tensor",
    "formatted",
    "to_html",
    "to_svg",
    "to_typst",
):
    assert hasattr(spenso, name), name

mink = spenso.Representation.mink(4)
tensor = spenso.TensorName("T")(mink)
expression = tensor.to_expression()

assert type(expression) is core.Expression
assert type(spenso.as_tensor(expression).to_expression()) is core.Expression
assert isinstance(tensor.to_typst(), str)
assert isinstance(spenso.to_typst(expression), str)
assert type(spenso.formatted(expression)) is core.FormattedOutput
assert isinstance(tensor.to_typst(True), str)
assert isinstance(spenso.to_typst(expression, True), str)
assert type(tensor.formatted(True)) is core.FormattedOutput
assert type(spenso.formatted(expression, True)) is core.FormattedOutput
for name in ("formatted", "to_html", "to_svg", "to_typst"):
    assert hasattr(tensor, name), name

# Explicit rendering reports the optional dependency, while notebook and
# FormattedOutput paths silently retain their existing LaTeX/text fallback.
import builtins
import sys

previous_typst = sys.modules.pop("typst", None)
original_import = builtins.__import__

def import_without_typst(name, globals=None, locals=None, fromlist=(), level=0):
    if name == "typst":
        raise ImportError("typst deliberately unavailable in this test")
    return original_import(name, globals, locals, fromlist, level)

builtins.__import__ = import_without_typst
try:
    try:
        spenso.to_html(expression)
    except ImportError as error:
        assert "gammaloop[typst-display]" in str(error)
    else:
        raise AssertionError("explicit HTML rendering must require typst-display")

    assert tensor._repr_html_() is None
    assert type(spenso.formatted(expression, True)) is core.FormattedOutput
finally:
    builtins.__import__ = original_import
    if previous_typst is not None:
        sys.modules["typst"] = previous_typst

# A fake compiler lets this test verify the trusted notation source and virtual
# file routing without installing typst-py or pretending to test its compiler.
import types

typst = types.ModuleType("typst")
typst.calls = []

def compile_typst(files, *, format, pretty):
    typst.calls.append((files, format, pretty))
    if format == "html":
        return b"<html><head><style>.math{}</style></head><body><math><mi>T</mi></math></body></html>"
    if format == "svg":
        return b"<svg xmlns='http://www.w3.org/2000/svg'></svg>"
    raise AssertionError(format)

typst.compile = compile_typst
sys.modules["typst"] = typst
notation_source = "#let trusted-notation-sentinel = 42"

html = spenso.to_html(expression, notation_source=notation_source)
assert "<math" in html
files, format, pretty = typst.calls[-1]
assert files["notation.typ"] == notation_source.encode()
assert files["render.typ"]
assert files["tree.cbor"]
assert format == "html"
assert pretty is True

svg = tensor.to_svg(notation_source=notation_source)
assert svg.startswith("<svg")
files, format, pretty = typst.calls[-1]
assert files["notation.typ"] == notation_source.encode()
assert format == "svg"
assert pretty is True

del sys.modules["typst"]
if previous_typst is not None:
    sys.modules["typst"] = previous_typst
"##,
        )
        .expect("Python assertions contain a null byte");
        py.run(&assertions, Some(&locals), Some(&locals))
    })
    .unwrap();
}

#[test]
fn native_payload_round_trip_preserves_nested_calls_and_spenso_attachments() {
    let namespace = "spynso_native_payload_round_trip";
    let representation_name = format!("{namespace}::M");
    let index_palette = IndexPalette::cyclic(
        1,
        [
            IndexDisplay::symbol("mu").unwrap(),
            IndexDisplay::symbol("nu").unwrap(),
        ],
    )
    .unwrap();
    let representation = LibraryRep::new_self_dual_with_index_palette_and_row(
        &representation_name,
        index_palette.clone(),
        IndexRow::Bottom,
    )
    .unwrap();
    let representation_declaration = RepresentationDeclaration::new(
        PortableRepresentationClass::SelfDual,
        index_palette,
        IndexRow::Bottom,
    );

    let manual_display = IndexDisplay::symbol("rho")
        .unwrap()
        .with_bottom(IndexDisplay::Number(2));
    let math_display_declaration = MathDisplayDeclaration::new(manual_display.clone()).unwrap();
    let display_name = canonical_math_display_symbol_name(&manual_display, namespace).unwrap();
    let display_symbol = register_math_display_symbol(&manual_display, namespace).unwrap();

    let nested = parse!("f(2*g(r-x))");
    let envelope_name = format!("{namespace}::payload");
    let envelope = SymbolBuilder::new(NamespacedSymbol::parse(&envelope_name))
        .build()
        .unwrap();
    let atom = FunctionBuilder::new(envelope)
        .add_arg(nested)
        .add_arg(Atom::var(representation.symbol()))
        .add_arg(Atom::var(display_symbol))
        .finish();

    let attachments = attachments_for_atom(&atom).unwrap();
    assert_eq!(attachments.len(), 2);
    assert!(
        attachments
            .iter()
            .any(|attachment| attachment.schema() == REPRESENTATION_ATTACHMENT_SCHEMA)
    );
    assert!(
        attachments
            .iter()
            .any(|attachment| attachment.schema() == MATH_DISPLAY_ATTACHMENT_SCHEMA)
    );

    let payload = encode_atom_from_set(&atom, &attachments).unwrap();
    let parsed = parse_payload(&payload).unwrap();
    let imported_attachments = parsed.attachment_set();
    assert_eq!(imported_attachments, attachments);
    let representations =
        RepresentationDeclarations::from_attachment_set(&imported_attachments).unwrap();
    let math_displays =
        MathDisplayDeclarations::from_attachment_set(&imported_attachments).unwrap();

    assert_eq!(
        representations.get(&representation_name),
        Some(&representation_declaration)
    );
    assert_eq!(
        math_displays.get(&display_name),
        Some(&math_display_declaration)
    );
    representations.preflight_registration().unwrap();
    math_displays.preflight_registration().unwrap();
    representations.register_before_atom_import().unwrap();
    math_displays.register_before_atom_import().unwrap();

    let imported = parsed.import_atom().unwrap();
    assert_eq!(imported, atom);
    assert_eq!(attachments_for_atom(&imported).unwrap(), attachments);
}
