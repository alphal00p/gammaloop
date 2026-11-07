//! Tensor notation leaves in Symbolica's responsive MathML renderer.
use crate::display;
use crate::display::{DisplaySettings, IndexAliases, TensorDisplayMode};
use idenso::tensor::{SymbolicTensor, composition};
use linnet_render_py::PyErrExt;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyAnyMethods, PyDict, PyDictMethods};
use pyo3::{Bound, Py, PyAny, PyErr, PyResult, Python};
use spenso::network::tags;
use spenso::network::tags::SPENSO_TAG;
use spenso::network::tags::mathml;
use spenso::network::tags::mathml::ExportError;
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use std::cell::RefCell;
use std::collections::VecDeque;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::api::python::MathMLViewer;
use symbolica::api::python::PythonExpression;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::domains::SelfRing;
use symbolica::printer::{PrintMode, PrintOptions, PrintState, PrintUserData};
use symbolica_notation::{MathMl, Notation};

const MAX_LEAF_NODES: usize = 256;
const MAX_LEAF_BYTES: usize = 4096;
const MAX_ALIASES: usize = 4096;

fn fits(
    atom: AtomView<'_>,
    budget: &mut usize,
    depth: usize,
    max_bytes: usize,
    max_terms: usize,
) -> bool {
    if *budget == 0
        || depth > 32
        || matches!(atom, AtomView::Add(sum) if sum.get_nargs() > max_terms)
        || atom.get_byte_size() > max_bytes
        || atom
            .get_symbol()
            .is_some_and(|symbol| symbol.get_name().len() > max_bytes)
    {
        return false;
    }
    *budget -= 1;
    atom.children()
        .all(|child| fits(child, budget, depth + 1, max_bytes, max_terms))
}

pub(crate) fn is_large(value: &SymbolicTensor<PartialStructure>) -> bool {
    // Text keeps the base limits, independently of notebook layout budgets.
    !fits(value.expression().as_view(), &mut 2000, 0, 64 * 1024, 100)
}

/// Notebook previews retain their smaller rendering threshold.
pub(crate) fn needs_paging(value: &SymbolicTensor<PartialStructure>) -> bool {
    !fits(
        value.expression().as_view(),
        &mut 1000,
        0,
        MAX_LEAF_BYTES,
        usize::MAX,
    )
}

type LeafCacheKey = (Atom, [u16; 10]);

impl DisplaySettings {
    /// Tensor presentation options shared by static exports and scoped paging.
    pub(crate) fn mathml_options(&self) -> PrintOptions {
        let mut options = display::display_options_with_settings(TensorDisplayMode::Typst, self);
        options.mode = PrintMode::MathML;
        for (key, value) in [
            ("spenso.tensor_layout", &self.tensor_layout),
            ("spenso.index_gap", &self.index_gap),
            ("spenso.factor_gap", &self.factor_gap),
        ] {
            options
                .custom_print_mode
                .insert(key.into(), PrintUserData::String(value.clone()));
        }
        options
    }
}

// Complete export and scoped paging share the tensor notation owner.
pub(super) fn html(
    py: Python<'_>,
    atom: &Atom,
    settings: &DisplaySettings,
    aliases: IndexAliases,
) -> PyResult<String> {
    let options = settings.mathml_options();
    // Validate settings even for a scalar root, before optional leaf callbacks.
    mathml::validate_options(&options).map_err(|error| export_error(py, error))?;
    let aliases = Arc::new(aliases);
    let mut body = None;
    // Avoid allocating a viewer for inputs beyond its hard node/depth limits.
    if fits(atom.as_view(), &mut 1000, 0, 64 * 1024, 500) {
        let leaf_aliases = Arc::clone(&aliases);
        let offered = Arc::new(AtomicUsize::new(0));
        let leaf_offers = Arc::clone(&offered);
        let leaf_options = settings.mathml_options();
        // The Python viewer requires an independent immutable expression owner.
        // This copy is bounded by the same small-input gate as its rendering.
        let expression = Py::new(py, PythonExpression::from(atom.clone()))?;
        let viewer = MathMLViewer::expression(
            expression,
            options
                .clone()
                .multiplication_operator('·')
                .num_exp_as_superscript(true),
            Some(Box::new(move |py, view, _, state| {
                if !TensorLeaves::semantic(view) {
                    return Ok(None);
                }
                let fragment = complete_body(view, &leaf_options, state, &leaf_aliases)
                    .map_err(|error| export_error(py, error))?;
                let serial = leaf_offers.fetch_add(1, Ordering::Relaxed);
                // A declined semantic offer can become a complete native call.
                // Track acceptance without changing the final leaf or its layout.
                Ok(Some(format!(
                    "<!--spenso-complete-leaf-{serial}-->{fragment}"
                )))
            })),
        )
        .page_size(500)
        .into_python(py)?;
        let result = viewer.call_method1(py, "_render_page", (0,));
        // Close the temporary source on both success and failure.
        let closed = viewer.call_method0(py, "close");
        let page = result?;
        closed?;
        body = complete_page(page.bind(py), offered.load(Ordering::Relaxed));
    }
    let body = match body {
        Some(body) => body,
        None => complete_body(atom.as_view(), &options, &PrintState::new(), &aliases)
            .map_err(|error| export_error(py, error))?,
    };
    Ok(display::notebook_html_fragment(&format!(
        "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"><mrow>{body}</mrow></math>"
    )))
}

// Symbolica's notebook page protocol has no public Rust page API at the pinned
// commit. Arithmetic tests guard against silently falling back after a change.
fn complete_page(page: &Bound<'_, PyAny>, offered: usize) -> Option<String> {
    let page = page.cast::<PyDict>().ok()?;
    let start: usize = page.get_item("start").ok()??.extract().ok()?;
    let end: usize = page.get_item("end").ok()??.extract().ok()?;
    let total: usize = page.get_item("total").ok()??.extract().ok()?;
    let markers: Vec<Py<PyAny>> = page.get_item("markers").ok()??.extract().ok()?;
    if start != 0 || end != total || !markers.is_empty() {
        return None;
    }
    let mut html: String = page.get_item("html").ok()??.extract().ok()?;
    for serial in 0..offered {
        let marker = format!("<!--spenso-complete-leaf-{serial}-->");
        if !html.contains(&marker) {
            return None;
        }
        html = html.replace(&marker, "");
    }
    if html.contains("<mtext data-symbolica-marker=") {
        return None;
    }
    let (root, rest) = html.split_once('>')?;
    if root != "<math" && !root.starts_with("<math ") {
        return None;
    }
    let fragment = rest.strip_suffix("</math>")?;
    MathMl::fragment(fragment.to_owned(), 256 * 1024)
}

fn export_error(py: Python<'_>, error: ExportError) -> PyErr {
    let exception = PyValueError::new_err(format!(
        "complete MathML export failed: {error}; use paged() for bounded inspection"
    ));
    if let ExportError::InvalidSettings {
        source: Some(source),
        ..
    } = error
    {
        exception.set_cause(
            py,
            Some(PyValueError::new_err(source.to_string()).with_rust_cause(&source)),
        );
    }
    exception
}

/// Complete fallback beyond Symbolica's printer bounds. This walk has no
/// sign/fraction layout; arithmetic layout belongs to Symbolica when it fits.
fn complete_body(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
    aliases: &IndexAliases,
) -> Result<String, ExportError> {
    let label_error = RefCell::new(None);
    let body = mathml::expression(atom, options, state, &|representation, index| {
        let key = (representation, index.to_owned());
        if let Some((_, display)) = aliases.entries.get(&key) {
            let rendered = display.to_mathml();
            if rendered.is_none() {
                label_error
                    .borrow_mut()
                    .get_or_insert(ExportError::IncompleteOutput);
            }
            return rendered;
        }
        aliases.source_indices.get(&key).and_then(|source| {
            let rendered =
                mathml::expression(source.as_view(), options, &PrintState::new(), &|_, _| None);
            match rendered {
                Ok(rendered) => Some(rendered),
                Err(error) => {
                    label_error.borrow_mut().get_or_insert(error);
                    None
                }
            }
        })
    });
    // Preserve callback failures instead of silently displaying encoded labels.
    label_error.into_inner().map_or(body, Err)
}

struct TensorLeaves {
    settings: DisplaySettings,
    aliases: IndexAliases,
    // Responsive alternatives share direct fragments, with bounded retention.
    // A cached fragment and the returned viewer response need independent owned
    // strings. Tensor metadata itself remains shared through the viewer's Arc.
    cache: VecDeque<(LeafCacheKey, Option<String>)>,
}

impl TensorLeaves {
    fn semantic(view: AtomView<'_>) -> bool {
        let Some(symbol) = view.get_symbol() else {
            return false;
        };
        symbol.has_tag(&SPENSO_TAG.tensor)
            || symbol.has_tag("tensor")
            || symbol.has_tag(&SPENSO_TAG.index)
            || symbol.has_tag(&SPENSO_TAG.representation)
            || tags::tensor_has_custom_print(symbol)
            || Notation::lookup(symbol)
                .expect("notation registry must not be poisoned")
                .is_some()
            || matches!(
                symbol.get_name(),
                "spenso::bracket"
                    | "spenso::cas"
                    | "spenso::idx"
                    | "spenso::gram"
                    | "spenso::gamma"
                    | "spenso::gamma0"
                    | "spenso::gamma5"
                    | "spenso::projp"
                    | "spenso::projm"
                    | "spenso::dot"
                    | "spenso::chain"
                    | "spenso::trace"
                    | "spenso::Nc"
            )
    }

    fn render(
        &mut self,
        value: &SymbolicTensor<PartialStructure>,
        view: AtomView<'_>,
        state: &PrintState,
    ) -> Option<String> {
        if !fits(view, &mut { MAX_LEAF_NODES }, 0, MAX_LEAF_BYTES, usize::MAX) {
            return None;
        }
        // Arithmetic and ordinary calls belong to Symbolica. Only semantic
        // tensor leaves enter the scoped semantic printer.
        if !Self::semantic(view) {
            return None;
        }
        let root = value.expression().as_view().get_data();
        // Custom printers receive the complete state. Reusing a fragment must
        // preserve every observable flag and level, including exponent context.
        let context = [
            u16::from(state.in_sum),
            u16::from(state.in_product),
            u16::from(state.suppress_one),
            u16::from(state.in_exp),
            u16::from(state.in_exp_base),
            u16::from(state.top_level_add_child),
            u16::from(state.superscript),
            state.level,
            state.bracket_level,
            state.indentation_level,
        ];
        if let Some(index) = self
            .cache
            .iter()
            .position(|((atom, flags), _)| atom.as_view() == view && *flags == context)
        {
            let entry = self
                .cache
                .remove(index)
                .expect("the cache position was found without intervening mutation");
            let result = entry.1.clone();
            self.cache.push_back(entry);
            return result;
        }
        // The inherited Python expression may own a byte-for-byte Atom clone.
        // Equality is bounded by fits above and does not scan a large root.
        let atom = if view.get_data() == root {
            value.presentation_atom()
        } else {
            value.presentation_fragment(view)
        };
        // Reserve the public free ports without walking the original expression.
        let mut bundle = FunctionBuilder::new(AIND_SYMBOLS.aind).add_arg(&atom);
        for slot in value.structure().logical_slots() {
            bundle = bundle.add_arg(composition::port_atom(slot));
        }
        // Resolve base labels before decorating scoped indices. Reallocating each
        // leaf independently can call two different bases a′, then discard the
        // second alias because scoped displays have no numeric palette position.
        let mut aliases = IndexAliases::for_atom_with_aliases(
            &bundle.finish(),
            &self.settings.index_style,
            &self.aliases,
        );
        for (representation, labels) in &aliases.reserved_labels {
            self.aliases
                .reserved_labels
                .entry(*representation)
                .or_default()
                .extend(labels.iter().cloned());
        }
        aliases.entries.retain(|key, alias| {
            if let Some(existing) = self.aliases.entries.get(key) {
                // A later explicit label must not acquire an earlier generated
                // alias. Preserve the source spelling in that unusual case.
                return !aliases
                    .reserved_labels
                    .get(&key.0)
                    .is_some_and(|labels| labels.contains(&IndexAliases::label_key(&existing.1)));
            }
            if self.aliases.entries.len() >= MAX_ALIASES {
                return false;
            }
            self.aliases
                .reserved_labels
                .entry(key.0)
                .or_default()
                .insert(IndexAliases::label_key(&alias.1));
            self.aliases.entries.insert(key.clone(), alias.clone());
            true
        });
        let options = self.settings.mathml_options();
        // Aliases remain borrowed viewer state, rather than serialized options
        // or symbols inserted into the expression being paged.
        let result = mathml::render(atom.as_view(), &options, state, &|representation, index| {
            let key = (representation, index.to_owned());
            aliases
                .entries
                .get(&key)
                .and_then(|(_, display)| display.to_mathml())
                .or_else(|| {
                    let source = aliases.source_indices.get(&key)?;
                    if !fits(
                        source.as_view(),
                        &mut { MAX_LEAF_NODES },
                        0,
                        MAX_LEAF_BYTES,
                        usize::MAX,
                    ) {
                        return None;
                    }
                    Some(source.format_string(&options, PrintState::new()))
                })
        });
        // Oversized calls defer to structural rendering so their arguments can
        // still be paged and inspected, rather than becoming an atomic loop.
        let result = result.filter(|html| html.len() <= MAX_LEAF_BYTES);
        self.cache
            .push_back(((view.to_owned(), context), result.clone()));
        if self.cache.len() > 128 {
            self.cache.pop_front();
        }
        result
    }
}

pub(crate) fn viewer(
    py: Python<'_>,
    expression: Py<PythonExpression>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
    page_size: Option<usize>,
) -> PyResult<Py<PyAny>> {
    if notation.is_some() {
        return Err(PyValueError::new_err(
            "The direct MathML pager does not accept Typst notation_source; use to_typst_project() or to_svg() for custom Typst exports",
        ));
    }
    // Reserve explicit labels across the complete descriptor once, before any
    // responsive layout can allocate a generated alias. Retain only that
    // inventory; generated aliases and decoded sources remain leaf-local/bounded.
    let aliases = IndexAliases {
        reserved_labels: IndexAliases::reserved_for_descriptor(&value, &settings.index_style),
        ..IndexAliases::default()
    };
    let mut leaves = TensorLeaves {
        settings,
        aliases,
        cache: VecDeque::new(),
    };
    let viewer = MathMLViewer::expression(
        expression,
        PrintOptions::new()
            .multiplication_operator('·')
            .num_exp_as_superscript(true),
        Some(Box::new(move |_py, atom, _options, state| {
            Ok(leaves.render(&value, atom, state))
        })),
    );
    match page_size {
        Some(size) => viewer.page_size(size),
        None => viewer.automatic_product_preview(),
    }
    .into_python(py)
}

pub(crate) fn static_html(
    py: Python<'_>,
    expression: Py<PythonExpression>,
    value: Arc<SymbolicTensor<PartialStructure>>,
    settings: DisplaySettings,
    notation: Option<String>,
) -> PyResult<String> {
    viewer(py, expression, value, settings, notation, None)?
        .call_method0(py, "_repr_html_")?
        .extract(py)
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
    use symbolica::atom::{NamespacedSymbol, SymbolBuilder};
    use symbolica::{parse, symbol};

    #[test]
    fn semantic_leaf_work_is_bounded_before_notation_callbacks() {
        let small = parse!("f(x,y)");
        assert!(fits(
            small.as_view(),
            &mut { MAX_LEAF_NODES },
            0,
            MAX_LEAF_BYTES,
            usize::MAX
        ));
        assert!(!fits(
            small.as_view(),
            &mut 1,
            0,
            MAX_LEAF_BYTES,
            usize::MAX
        ));
        assert!(!fits(
            small.as_view(),
            &mut { MAX_LEAF_NODES },
            33,
            MAX_LEAF_BYTES,
            usize::MAX
        ));
    }

    #[test]
    fn text_limits_are_independent_of_notebook_leaf_limits() {
        let head = symbol!("paging_text_tests::f");
        let wide = FunctionBuilder::new(head)
            .add_args((0..1024).map(Atom::num))
            .finish();
        assert!(fits(wide.as_view(), &mut 2000, 0, 64 * 1024, 100));
        assert!(!fits(
            wide.as_view(),
            &mut 1000,
            0,
            MAX_LEAF_BYTES,
            usize::MAX
        ));
        let sum = Atom::add_many(
            (0..101).map(|i| FunctionBuilder::new(head).add_arg(Atom::num(i)).finish()),
        );
        assert!(!fits(sum.as_view(), &mut 2000, 0, 64 * 1024, 100));
        assert!(fits(sum.as_view(), &mut 2000, 0, 64 * 1024, usize::MAX));
        let too_wide = FunctionBuilder::new(head)
            .add_args((0..2000).map(Atom::num))
            .finish();
        assert!(!fits(too_wide.as_view(), &mut 2000, 0, 64 * 1024, 100));
    }

    #[test]
    fn complete_arithmetic_uses_native_signs_and_fractions() {
        Python::initialize();
        Python::attach(|py| {
            let mut cases = [
                ("a-b", vec!["a", "−", "b"]),
                ("x-3", vec!["−", "3", "+", "x"]),
                ("-x", vec!["−", "x"]),
                ("a-2*b", vec!["a", "−", "2", "b"]),
                ("x/y", vec!["x", "y"]),
            ]
            .into_iter()
            .map(|(input, tokens)| (parse!(input), tokens))
            .collect::<Vec<_>>();
            let representation = LibraryRep::from(Minkowski {});
            let port = representation
                .to_symbolic([Atom::num(4), Atom::var(symbol!("arithmetic_tests::i"))]);
            let left = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("arithmetic_tests::A"))
                .add_arg(&port)
                .finish();
            let right = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("arithmetic_tests::B"))
                .add_arg(port)
                .finish();
            cases.push((left - right, vec!["A", "i", "−", "B", "i"]));
            for (atom, tokens) in cases {
                let rendered = html(
                    py,
                    &atom,
                    &DisplaySettings::default(),
                    IndexAliases::default(),
                )
                .unwrap();
                let math = rendered
                    .split_once("<math ")
                    .unwrap()
                    .1
                    .split_once("</math>")
                    .unwrap()
                    .0;
                let visible = math
                    .split('<')
                    .filter_map(|part| {
                        let (tag, text) = part.split_once('>')?;
                        (matches!(tag.split_whitespace().next(), Some("mi" | "mn" | "mo"))
                            && text != "&#x2062;")
                            .then_some(text)
                    })
                    .collect::<Vec<_>>();
                assert_eq!(visible, tokens, "{atom}: {rendered}");
                if tokens == ["x", "y"] {
                    assert!(rendered.contains("<mfrac>"));
                }
            }
        });
    }

    #[test]
    fn complete_fallback_preserves_every_term_when_viewer_cannot_fit() {
        Python::initialize();
        Python::attach(|py| {
            let atom = Atom::add_many((0..600).map(|i| {
                FunctionBuilder::new(symbol!("fallback_tests::f"))
                    .add_arg(Atom::num(i))
                    .finish()
            }));
            let rendered = html(
                py,
                &atom,
                &DisplaySettings::default(),
                IndexAliases::default(),
            )
            .unwrap();
            assert!(!fits(atom.as_view(), &mut 1000, 0, 64 * 1024, 500));
            for i in 0..600 {
                assert!(rendered.contains(&format!("<mn>{i}</mn>")));
            }
            assert!(!rendered.contains("⋯"));
            assert!(!rendered.contains("data-symbolica-marker"));
        });
    }

    #[test]
    fn changed_notebook_page_protocol_declines_without_panicking() {
        Python::initialize();
        Python::attach(|py| {
            assert!(complete_page(py.None().bind(py), 0).is_none());
            let page = PyDict::new(py);
            assert!(complete_page(page.as_any(), 0).is_none());
            page.set_item("start", 0).unwrap();
            page.set_item("end", 1).unwrap();
            page.set_item("total", 1).unwrap();
            page.set_item("markers", Vec::<usize>::new()).unwrap();
            page.set_item("html", "<math><mi>x</mi></math>").unwrap();
            assert!(complete_page(page.as_any(), 0).is_some());
            for key in ["start", "end", "total", "markers", "html"] {
                let saved = page.get_item(key).unwrap().unwrap();
                page.del_item(key).unwrap();
                assert!(complete_page(page.as_any(), 0).is_none(), "{key}");
                page.set_item(key, py.None()).unwrap();
                assert!(complete_page(page.as_any(), 0).is_none(), "{key}");
                page.set_item(key, saved).unwrap();
            }
            for root in ["<div><mi>x</mi></math>", "<math><mi>x</mi>", "garbage"] {
                page.set_item("html", root).unwrap();
                assert!(complete_page(page.as_any(), 0).is_none(), "{root}");
            }
        });
    }

    #[test]
    fn complete_export_preserves_decoded_label_errors() {
        let representation = LibraryRep::from(Minkowski {});
        let index = Atom::var(symbol!("paging_mathml_tests::i"));
        let tensor = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("paging_mathml_tests::T"))
            .add_arg(representation.to_symbolic([Atom::num(4), index.clone()]))
            .finish();
        let source = Atom::var(
            SymbolBuilder::new(NamespacedSymbol::parse(&format!(
                "paging_mathml_tests::{}",
                "x".repeat(5000)
            )))
            .build()
            .unwrap(),
        );
        let mut aliases = IndexAliases::default();
        aliases
            .source_indices
            .insert((representation.symbol(), index), source);
        Python::initialize();
        Python::attach(|py| {
            let error = html(py, &tensor, &DisplaySettings::default(), aliases).unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            assert!(error.to_string().contains("source exceeds"), "{error}");
        });
    }
}
