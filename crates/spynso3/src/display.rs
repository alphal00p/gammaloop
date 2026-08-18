use pyo3::prelude::*;
use spenso::{
    algebra::complex::RealOrComplexRef,
    iterators::IteratableTensor,
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        HasStructure, TensorStructure,
        abstract_index::AbstractIndex,
        permuted::PermuteTensor,
        representation::LibraryRep,
        slot::{IsAbstractSlot, Slot},
    },
    tensors::{data::SparseOrDense, parametric::AtomViewOrConcrete},
};
use symbolica::{
    api::python::{PythonExpression, PythonFormattedOutput},
    atom::Atom,
    domains::SelfRing,
    printer::{PrintOptions, PrintState},
};

use crate::{Spensor, composition::StructuredAtom};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::gen_stub_pyfunction;

#[derive(Clone, Copy)]
enum TensorDisplayMode {
    Plain,
    Latex,
    Typst,
}

fn display_options(mode: TensorDisplayMode, show_dimensions: bool) -> PrintOptions {
    let mut presentation = match mode {
        TensorDisplayMode::Plain => SpensoPrintSettings::compact(),
        TensorDisplayMode::Latex | TensorDisplayMode::Typst => SpensoPrintSettings::typst(),
    };
    presentation.with_dim = show_dimensions;

    let mut options = match mode {
        TensorDisplayMode::Plain => presentation.nice_symbolica(),
        TensorDisplayMode::Latex => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::latex()
        },
        TensorDisplayMode::Typst => PrintOptions {
            custom_print_mode: presentation.into(),
            ..PrintOptions::typst()
        },
    };
    options.color_builtin_symbols = false;
    options.color_top_level_sum = false;
    options.terms_on_new_line = false;
    options.max_line_length = None;
    options
}

fn format_atom_with_mode(atom: &Atom, mode: TensorDisplayMode, show_dimensions: bool) -> String {
    atom.format_string(&display_options(mode, show_dimensions), PrintState::new())
}

fn format_structured_with_mode(
    value: &StructuredAtom,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    format_atom_with_mode(&value.presentation_atom(), mode, show_dimensions)
}

pub(crate) fn format_atom(atom: &Atom, show_dimensions: bool) -> String {
    format_atom_with_mode(atom, TensorDisplayMode::Plain, show_dimensions)
}

pub(crate) fn atom_to_typst(atom: &Atom, show_dimensions: bool) -> String {
    format_atom_with_mode(atom, TensorDisplayMode::Typst, show_dimensions)
}

pub(crate) fn atom_to_latex(atom: &Atom, show_dimensions: bool) -> String {
    let body = format_atom_with_mode(atom, TensorDisplayMode::Latex, show_dimensions);
    format!("$${body}$$")
}

pub(crate) fn format_atom_output(atom: &Atom, show_dimensions: bool) -> PythonFormattedOutput {
    PythonFormattedOutput {
        text: format_atom(atom, show_dimensions),
        // Symbolica does not yet expose its Typst-to-SVG renderer to community
        // modules. Keep rich construction here so that renderer is one call to
        // wire once the public dependency API is available.
        html: None,
        latex: Some(atom_to_latex(atom, show_dimensions)),
    }
}

pub(crate) fn format_structured(value: &StructuredAtom, show_dimensions: bool) -> String {
    format_structured_with_mode(value, TensorDisplayMode::Plain, show_dimensions)
}

pub(crate) fn structured_to_typst(value: &StructuredAtom, show_dimensions: bool) -> String {
    format_structured_with_mode(value, TensorDisplayMode::Typst, show_dimensions)
}

pub(crate) fn structured_to_latex(value: &StructuredAtom, show_dimensions: bool) -> String {
    let body = format_structured_with_mode(value, TensorDisplayMode::Latex, show_dimensions);
    format!("$${body}$$")
}

pub(crate) fn format_structured_output(
    value: &StructuredAtom,
    show_dimensions: bool,
) -> PythonFormattedOutput {
    PythonFormattedOutput {
        text: format_structured(value, show_dimensions),
        // Symbolica does not yet expose its Typst-to-SVG renderer to community
        // modules. Keep rich construction here so that renderer is one call to
        // wire once the public dependency API is available.
        html: None,
        latex: Some(structured_to_latex(value, show_dimensions)),
    }
}

fn format_tensor_value(
    value: AtomViewOrConcrete<'_, RealOrComplexRef<'_, f64>>,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    match value {
        AtomViewOrConcrete::Atom(atom) => {
            format_atom_with_mode(&atom.to_owned(), mode, show_dimensions)
        }
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Real(value)) => value.to_string(),
        AtomViewOrConcrete::Concrete(RealOrComplexRef::Complex(value)) => value.to_string(),
    }
}

fn format_tensor_interface(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    logical_tensor_slots(tensor)
        .into_iter()
        .map(|slot| format_atom_with_mode(&slot.to_atom(), mode, show_dimensions))
        .collect::<Vec<_>>()
        .join(" ")
}

fn logical_tensor_slots(tensor: &Spensor) -> Vec<Slot<LibraryRep, AbstractIndex>> {
    let canonical = tensor
        .tensor
        .structure
        .structure()
        .external_structure_iter()
        .collect::<Vec<_>>();
    let representation_sorted = tensor.tensor.index_permutation.apply_slice_inv(&canonical);
    tensor
        .tensor
        .rep_permutation
        .apply_slice_inv(&representation_sorted)
}

fn concrete_tensor_body(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    let dense = tensor
        .tensor
        .structure
        .clone()
        .to_dense()
        .permute_inds(&tensor.tensor.index_permutation.inverse())
        .permute_inds(&tensor.tensor.rep_permutation.inverse());
    let values = dense
        .iter_flat()
        .map(|(_, value)| format_tensor_value(value, mode, show_dimensions))
        .collect::<Vec<_>>();
    let shape = logical_tensor_slots(tensor)
        .into_iter()
        .map(|slot| usize::try_from(slot.dim()))
        .collect::<Result<Vec<_>, _>>()
        .ok();

    match (mode, shape.as_deref()) {
        (_, Some([])) => values.first().cloned().unwrap_or_else(|| "0".to_string()),
        (TensorDisplayMode::Plain, Some([_])) => format!("[{}]", values.join(", ")),
        (TensorDisplayMode::Plain, Some([_, columns])) if *columns > 0 => values
            .chunks(*columns)
            .map(|row| format!("[{}]", row.join(", ")))
            .collect::<Vec<_>>()
            .join("\n"),
        (TensorDisplayMode::Plain, _) => format!("[{}]", values.join(", ")),
        (TensorDisplayMode::Typst, Some([_])) => format!("vec({})", values.join(",")),
        (TensorDisplayMode::Typst, Some([_, columns])) if *columns > 0 => {
            let rows = values
                .chunks(*columns)
                .map(|row| row.join(","))
                .collect::<Vec<_>>();
            format!("mat({})", rows.join(";"))
        }
        (TensorDisplayMode::Typst, _) => {
            format!(r#"op("Tensor")({})"#, values.join(","))
        }
        (TensorDisplayMode::Latex, Some([_])) => {
            format!(
                r#"\begin{{pmatrix}}{}\end{{pmatrix}}"#,
                values.join(r#" \\ "#)
            )
        }
        (TensorDisplayMode::Latex, Some([_, columns])) if *columns > 0 => {
            let rows = values
                .chunks(*columns)
                .map(|row| row.join(" & "))
                .collect::<Vec<_>>();
            format!(
                r#"\begin{{pmatrix}}{}\end{{pmatrix}}"#,
                rows.join(r#" \\ "#)
            )
        }
        (TensorDisplayMode::Latex, _) => {
            format!(
                r#"\operatorname{{Tensor}}\left({}\right)"#,
                values.join(",")
            )
        }
    }
}

fn format_concrete_with_mode(
    tensor: &Spensor,
    mode: TensorDisplayMode,
    show_dimensions: bool,
) -> String {
    let interface = format_tensor_interface(tensor, mode, show_dimensions);
    let body = concrete_tensor_body(tensor, mode, show_dimensions);
    if interface.is_empty() {
        return body;
    }

    match mode {
        TensorDisplayMode::Plain => format!("Tensor_({interface})\n{body}"),
        TensorDisplayMode::Typst => format!("attach({body},b:({interface}))"),
        TensorDisplayMode::Latex => format!("$${body}_{{{interface}}}$$"),
    }
}

pub(crate) fn format_concrete_tensor(tensor: &Spensor, show_dimensions: bool) -> String {
    format_concrete_with_mode(tensor, TensorDisplayMode::Plain, show_dimensions)
}

pub(crate) fn concrete_tensor_to_typst(tensor: &Spensor, show_dimensions: bool) -> String {
    format_concrete_with_mode(tensor, TensorDisplayMode::Typst, show_dimensions)
}

pub(crate) fn format_concrete_tensor_output(
    tensor: &Spensor,
    show_dimensions: bool,
) -> PythonFormattedOutput {
    PythonFormattedOutput {
        text: format_concrete_tensor(tensor, show_dimensions),
        html: None,
        latex: Some(format_concrete_with_mode(
            tensor,
            TensorDisplayMode::Latex,
            show_dimensions,
        )),
    }
}

/// Format a tensor expression using compact Spenso notation.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = false))]
fn format_tensor(expression: &PythonExpression, show_dimensions: bool) -> String {
    format_atom(&expression.expr, show_dimensions)
}

/// Format a tensor expression as Typst math source.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = false))]
fn to_typst(expression: &PythonExpression, show_dimensions: bool) -> String {
    atom_to_typst(&expression.expr, show_dimensions)
}

/// Build Symbolica's rich display wrapper for a tensor expression.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.spenso")
)]
#[pyfunction]
#[pyo3(signature = (expression, show_dimensions = false))]
fn formatted(expression: &PythonExpression, show_dimensions: bool) -> PythonFormattedOutput {
    format_atom_output(&expression.expr, show_dimensions)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(format_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(to_typst, m)?)?;
    m.add_function(wrap_pyfunction!(formatted, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::{
        network::tags::SPENSO_TAG,
        structure::{
            abstract_index::AbstractIndex,
            dimension::Dimension,
            representation::{LibraryRep, Representation},
        },
        vector_symbol,
    };
    use symbolica::{atom::AtomCore, function, symbol};

    fn vector() -> Atom {
        let representation = Representation {
            rep: LibraryRep::SelfDual(1),
            dim: Dimension::Concrete(4),
        };
        function!(
            vector_symbol!("display_test_vector"),
            Atom::num(1),
            representation
                .slot::<AbstractIndex, _>(AbstractIndex::from(symbol!("mu")))
                .to_atom()
        )
    }

    #[test]
    fn dimensions_are_opt_in_and_formatting_does_not_change_the_atom() {
        let atom = vector();
        let canonical = atom.to_canonical_string();

        let compact = format_atom(&atom, false);
        let dimensioned = format_atom(&atom, true);

        assert_ne!(compact, dimensioned);
        assert!(!compact.contains('4'));
        assert!(dimensioned.contains('4'));
        assert_eq!(atom.to_canonical_string(), canonical);
    }

    #[test]
    fn typst_automatically_uses_spenso_printers() {
        let atom = vector();
        let automatic = atom.printer(PrintOptions::typst()).to_string();

        assert!(!automatic.contains("spenso::"));
        assert!(symbol!("display_test_vector").has_tag(&SPENSO_TAG.rank1));
    }
}
