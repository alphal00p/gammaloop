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
use symbolica_typst_plugin::payload::{encode_atom_from_set, parse_payload};

const SPENSO_WRAPPER: &str = r#"
from ..tensor_native import *

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
    assert_eq!(SpensoModule::get_name(), "tensor");
    let native_name = format!("{}_native", SpensoModule::get_name());
    let native = PyModule::new(py, &native_name)?;
    let initialize_module =
        PyCFunction::new_closure(py, Some(c"initialize_module"), None, |args, _kwargs| {
            SpensoModule::initialize(args.py())
        })?;
    native.add("initialize_module", initialize_module)?;
    SpensoModule::register_module(&native)?;
    core.add_submodule(&native)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("symbolica.community.{native_name}"), &native)?;
    Ok(native)
}

fn import_spenso_wrapper<'py>(
    py: Python<'py>,
    community: &Bound<'py, PyModule>,
) -> PyResult<Bound<'py, PyModule>> {
    let name = "symbolica.community.tensor";
    let wrapper = PyModule::new(py, name)?;
    wrapper.setattr("__package__", name)?;
    wrapper.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(name, &wrapper)?;
    let source = CString::new(SPENSO_WRAPPER).expect("Spenso wrapper contains a null byte");
    py.run(&source, Some(&wrapper.dict()), Some(&wrapper.dict()))?;
    community.add("tensor", &wrapper)?;
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
        assert!(!community.hasattr("spenso")?);
        for name in ["Tensor", "TensorExpression", "Representation"] {
            assert_eq!(
                spenso
                    .getattr(name)?
                    .getattr("__module__")?
                    .extract::<String>()?,
                "symbolica.community.tensor"
            );
        }

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

# Rendering is self-contained even when the unrelated PyPI Linnet and
# typst-py cannot be imported. Trusted notation still crosses the real compiler.
import builtins
import sys

original_import = builtins.__import__
def import_without_renderers(name, globals=None, locals=None, fromlist=(), level=0):
    if name in ("typst", "linnet"):
        raise ImportError(f"{name} deliberately unavailable in this test")
    return original_import(name, globals, locals, fromlist, level)

builtins.__import__ = import_without_renderers
try:
    html = spenso.to_html(expression)
    assert "<math" in html
    assert "data:font" not in html
    assert len(html.encode()) < 10_000
    assert "<math" in tensor._repr_html_()
    assert tensor.to_svg().startswith("<svg")
    assert type(spenso.formatted(expression, True)) is core.FormattedOutput
    for render in (spenso.to_html, spenso.to_svg):
        try:
            render(expression, notation_source='#panic("trusted-notation-sentinel")')
        except RuntimeError as error:
            assert "trusted-notation-sentinel" in str(error), str(error)
        else:
            raise AssertionError("custom notation was not compiled")
finally:
    builtins.__import__ = original_import

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
