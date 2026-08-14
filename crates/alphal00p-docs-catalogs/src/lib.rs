//! Explicit, ordered supported-API catalogs for the five documentation sites.
//!
//! These adapters intentionally do not use a process-global inventory. Each
//! component exporter constructs its own scope, which keeps registration
//! deterministic and avoids pulling incompatible PyO3 feature sets together.

mod annotated_items;
#[cfg(feature = "gammaloop-reference")]
pub mod gammaloop_reference;
pub mod generated;
mod python_stub;
#[cfg(feature = "vakint-reference")]
pub mod vakint_reference;

use std::{
    fs,
    path::{Path, PathBuf},
};

use alphal00p_docs_schema::{
    ApiLanguage, DocCatalog, DocComponent, DocExample, DocFormat, DocItem, DocItemKind, DocProduct,
    DocScope, DocText, SourceLocation,
};
use eyre::{ContextCompat, Result, WrapErr, bail, ensure};
use python_stub::python_declarations;

// Keep the proc-macro integration in this internal adapter crate so the
// independently publishable product crates do not acquire documentation-only
// dependencies. Each descriptor is constructed explicitly below, mirroring
// Typst's ordered runtime scope instead of relying on a global inventory.
mod annotated_scopes {
    /// GammaLoop's supported Rust entry points. The complete compiled surface
    /// remains available in the #link("reference/rust/")[Rustdoc sidecar].
    #[alphal00p_docs::scope(
        id = "gammalooprs",
        title = "gammalooprs supported API",
        format = "typst"
    )]
    mod gammalooprs {}

    /// GammaLoop's stable facade and command API. Use this curated boundary for
    /// applications while consulting the #link("reference/rust/")[full Rustdoc]
    /// for exhaustive details.
    #[alphal00p_docs::scope(
        id = "gammaloop-api",
        title = "gammaloop-api supported API",
        format = "typst"
    )]
    mod gammaloop_api {}

    /// Linnet's supported half-edge graph API. Algorithms are grouped in the
    /// order used by the manual and #link("reference/rust/")[full Rustdoc].
    #[alphal00p_docs::scope(id = "linnet", title = "linnet supported API", format = "typst")]
    mod linnet {}

    /// Spenso's supported tensor and network API. The scope is #emph[intentionally]
    /// ordered so foundational tensor types precede planning and execution.
    #[alphal00p_docs::scope(id = "spenso", title = "spenso supported API", format = "typst")]
    mod spenso {}

    /// Spenso's exported procedural macros. Expansion details remain in the
    /// #link("reference/rust/")[exhaustive Rustdoc sidecar].
    #[alphal00p_docs::scope(
        id = "spenso-macros",
        title = "spenso-macros supported API",
        format = "typst"
    )]
    mod spenso_macros {}

    /// Spenso's high-energy-physics companion library. These entries #emph[complement]
    /// the tensor concepts owned by the main Spenso documentation.
    #[alphal00p_docs::scope(
        id = "spenso-hep-lib",
        title = "spenso-hep-lib supported API",
        format = "typst"
    )]
    mod spenso_hep_lib {}

    /// Idenso's supported representation and algebra API. The manual owns the
    /// representation conventions; this scope provides #emph[stable] source links.
    #[alphal00p_docs::scope(id = "idenso", title = "idenso supported API", format = "typst")]
    mod idenso {}

    /// Vakint's supported topology and evaluation API. Generated topology and
    /// backend tables are #emph[validated] against runtime sources.
    #[alphal00p_docs::scope(id = "vakint", title = "vakint supported API", format = "typst")]
    mod vakint {}

    pub(super) fn for_component(component: &str) -> Option<alphal00p_docs_schema::DocScope> {
        match component {
            "gammalooprs" => Some(__alphal00p_docs_scope_gammalooprs()),
            "gammaloop-api" => Some(__alphal00p_docs_scope_gammaloop_api()),
            "linnet" => Some(__alphal00p_docs_scope_linnet()),
            "spenso" => Some(__alphal00p_docs_scope_spenso()),
            "spenso-macros" => Some(__alphal00p_docs_scope_spenso_macros()),
            "spenso-hep-lib" => Some(__alphal00p_docs_scope_spenso_hep_lib()),
            "idenso" => Some(__alphal00p_docs_scope_idenso()),
            "vakint" => Some(__alphal00p_docs_scope_vakint()),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CatalogRequest {
    pub product_id: String,
    pub product_title: String,
    pub component_id: String,
    pub package: String,
    pub component_title: String,
    pub version: String,
    pub language: ApiLanguage,
    pub module: Option<String>,
    pub features: Vec<String>,
}

pub fn export_catalog(
    request: &CatalogRequest,
    workspace_root: &Path,
    python_stub: Option<&Path>,
) -> Result<DocCatalog> {
    let root = match request.language {
        ApiLanguage::Rust => rust_scope(&request.component_id, workspace_root)?,
        ApiLanguage::Python => {
            let stub = python_stub.context("a checked-in Python stub is required")?;
            python_scope(request, workspace_root, stub)?
        }
        language => bail!("unsupported catalog exporter language {language}"),
    };
    let product = DocProduct::new(&request.product_id, &request.product_title);
    let mut component = DocComponent::new(
        &request.component_id,
        &request.package,
        &request.component_title,
        &request.version,
        request.language,
    );
    component.features.clone_from(&request.features);
    let catalog = DocCatalog::new(product, component, root);
    catalog.validate_supported()?;
    Ok(catalog)
}

#[derive(Clone, Copy)]
struct RustExampleSpec {
    id: &'static str,
    language: &'static str,
    code: &'static str,
    signature_override: Option<&'static str>,
}

macro_rules! example {
    ($id:literal, $language:literal, $code:literal) => {
        RustExampleSpec {
            id: $id,
            language: $language,
            code: $code,
            signature_override: None,
        }
    };
    ($id:literal, $language:literal, $code:literal, signature = $signature:literal) => {
        RustExampleSpec {
            id: $id,
            language: $language,
            code: $code,
            signature_override: Some($signature),
        }
    };
}

fn rust_scope(component: &str, workspace_root: &Path) -> Result<DocScope> {
    let examples = rust_examples(component)
        .wrap_err_with(|| format!("no supported Rust examples registered for {component}"))?;
    let items = annotated_items::for_component(component)
        .wrap_err_with(|| format!("no source-backed items registered for {component}"))?;
    ensure!(
        examples.len() == items.len(),
        "{component} has {} example records but {} source-backed descriptors",
        examples.len(),
        items.len()
    );
    let mut root = annotated_scopes::for_component(component)
        .wrap_err_with(|| format!("no annotated documentation scope registered for {component}"))?;
    let mut supported = DocScope::new("supported", "Supported entry points");
    for (example, mut item) in examples.iter().zip(items) {
        ensure!(
            item.id == example.id,
            "{component} example {} is paired with source-backed item {}",
            example.id,
            item.id
        );
        let source = item
            .source
            .as_ref()
            .context("source-backed descriptor has no source location")?;
        ensure!(
            !source.identifier.is_empty()
                && source.line > 0
                && workspace_root.join(&source.file).is_file(),
            "{component} item {} has an invalid source-backed location",
            item.id
        );
        ensure!(
            item.signature
                .as_ref()
                .is_some_and(|value| !value.is_empty()),
            "{component} source-backed item {} has no signature",
            item.id
        );
        if !item
            .examples
            .iter()
            .any(|candidate| candidate.code == example.code)
        {
            item.examples.push(DocExample::new(
                "Supported usage",
                example.language,
                example.code,
            ));
        }
        if let Some(signature) = example.signature_override {
            item.signature = Some(signature.to_owned());
        }
        item.required_features = rust_item_required_features(component, example.id)
            .iter()
            .map(|feature| (*feature).to_owned())
            .collect();
        supported.define_item(item)?;
    }
    root.define_scope(supported)?;
    Ok(root)
}

fn rust_item_required_features(component: &str, item: &str) -> &'static [&'static str] {
    match (component, item) {
        // The entire parametric module is behind Spenso's Symbolica-backed
        // `shadowing` feature, not merely extra implementations on the type.
        ("spenso", "ParamTensor") => &["shadowing"],
        _ => &[],
    }
}

fn rust_examples(component: &str) -> Result<&'static [RustExampleSpec]> {
    match component {
        "gammalooprs" => Ok(GAMMALOOPRS),
        "gammaloop-api" => Ok(GAMMALOOP_API),
        "linnet" => Ok(LINNET),
        "spenso" => Ok(SPENSO),
        "spenso-macros" => Ok(SPENSO_MACROS),
        "spenso-hep-lib" => Ok(SPENSO_HEP_LIB),
        "idenso" => Ok(IDENSO),
        "vakint" => Ok(VAKINT),
        _ => bail!("unknown Rust component {component}"),
    }
}

const GAMMALOOPRS: &[RustExampleSpec] = &[
    example!(
        "GammaLoopContext",
        "rust",
        "fn accepts_context<C: gammalooprs::GammaLoopContext>(_context: &C) {}"
    ),
    example!(
        "HasModel",
        "rust",
        "fn model_owner<C: gammalooprs::HasModel>(_context: &C) {}"
    ),
    example!(
        "set_interrupt_handler",
        "rust",
        "gammalooprs::set_interrupt_handler();"
    ),
];

const GAMMALOOP_API: &[RustExampleSpec] = &[
    example!(
        "StateLoadOption",
        "rust",
        "let options = gammaloop_api::StateLoadOption::default();"
    ),
    example!(
        "StateLoadOption::load",
        "rust",
        "let loaded = gammaloop_api::StateLoadOption::default().load()?;"
    ),
    example!(
        "LoadedState::cli_session",
        "rust",
        "let mut loaded = gammaloop_api::StateLoadOption::default().load()?; let session = loaded.cli_session();"
    ),
    example!(
        "CLISettings",
        "rust",
        "let settings = gammaloop_api::CLISettings::default();"
    ),
    example!(
        "gammaloop",
        "console",
        "gammaloop --help",
        signature = "gammaloop [OPTIONS] [COMMAND]"
    ),
];

const LINNET: &[RustExampleSpec] = &[
    example!("HedgeGraph", "rust", "use linnet::half_edge::HedgeGraph;"),
    example!(
        "HedgeGraphBuilder",
        "rust",
        "use linnet::half_edge::builder::HedgeGraphBuilder;"
    ),
    example!(
        "SuBitGraph",
        "rust",
        "use linnet::half_edge::subgraph::SuBitGraph;"
    ),
    example!(
        "SimpleTraversalTree",
        "rust",
        "use linnet::half_edge::tree::SimpleTraversalTree;"
    ),
    example!("DotGraph", "rust", "use linnet::parser::DotGraph;"),
];

const SPENSO: &[RustExampleSpec] = &[
    example!(
        "TensorStructure",
        "rust",
        "fn accepts_structure<S: spenso::structure::TensorStructure>(_structure: &S) {}"
    ),
    example!(
        "DenseTensor",
        "rust",
        "use spenso::tensors::data::DenseTensor;"
    ),
    example!(
        "SparseTensor",
        "rust",
        "use spenso::tensors::data::SparseTensor;"
    ),
    example!(
        "Contract",
        "rust",
        "fn can_contract<T: spenso::contraction::Contract>(_tensor: &T) {}"
    ),
    example!("Network", "rust", "use spenso::network::Network;"),
    example!(
        "ParamTensor",
        "rust",
        "use spenso::tensors::parametric::ParamTensor;"
    ),
];

const SPENSO_MACROS: &[RustExampleSpec] = &[example!(
    "SimpleRepresentation",
    "rust",
    "#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default, spenso_macros::SimpleRepresentation)] #[representation(name = \"flavor\", self_dual)] struct Flavor;"
)];

const SPENSO_HEP_LIB: &[RustExampleSpec] = &[
    example!("hep_lib", "rust", "use spenso_hep_lib::hep_lib;"),
    example!(
        "gamma_data_dirac",
        "rust",
        "use spenso_hep_lib::gamma_data_dirac;"
    ),
    example!(
        "su3_generator_data",
        "rust",
        "use spenso_hep_lib::su3_generator_data;"
    ),
];

const IDENSO: &[RustExampleSpec] = &[
    example!(
        "IndexTooling",
        "rust",
        "fn accepts_index_tools<T: idenso::IndexTooling>(_expression: &T) {}"
    ),
    example!(
        "CookSettings",
        "rust",
        "let settings = idenso::CookSettings::reversible();"
    ),
    example!(
        "Cookable",
        "rust",
        "fn accepts_cookable<T: idenso::Cookable>(_expression: &T) {}"
    ),
    example!(
        "SelectiveExpand",
        "rust",
        "fn accepts_expansion<T: idenso::selective_expand::SelectiveExpand>(_expression: &T) {}"
    ),
    example!("bis", "rust", "let representation = idenso::bis!(4);"),
    example!("cof", "rust", "let representation = idenso::cof!(Nc);"),
    example!("coad", "rust", "let representation = idenso::coad!(Na);"),
    example!("gamma", "rust", "let tensor = idenso::gamma!(1, 2, 3);"),
    example!("gamma0", "rust", "let tensor = idenso::gamma0!(1, 2);"),
    example!("gamma5", "rust", "let tensor = idenso::gamma5!(1, 2);"),
    example!(
        "epsilon",
        "rust",
        "let tensor = idenso::epsilon!(1, 2, 3, 4);"
    ),
    example!("color_t", "rust", "let tensor = idenso::color_t!(1, 2, 3);"),
    example!("color_f", "rust", "let tensor = idenso::color_f!(1, 2, 3);"),
    example!(
        "representations::initialize",
        "rust",
        "idenso::representations::initialize();"
    ),
];

const VAKINT: &[RustExampleSpec] = &[
    example!(
        "vakint_parse",
        "rust",
        "let expression = vakint::vakint_parse!(\"x\").expect(\"valid Vakint expression\"); let explicit = vakint::vakint_parse!(\"x\", vakint::NAMESPACE).expect(\"valid explicit namespace\");"
    ),
    example!(
        "vakint_symbol",
        "rust",
        "let symbol = vakint::vakint_symbol!(\"mass\");"
    ),
    example!("Vakint", "rust", "let engine = vakint::Vakint::new()?;"),
    example!(
        "VakintSettings",
        "rust",
        "let settings = vakint::VakintSettings::default();"
    ),
    example!("VakintExpression", "rust", "use vakint::VakintExpression;"),
    example!(
        "EvaluationOrder",
        "rust",
        "let order = vakint::EvaluationOrder::empty();"
    ),
    example!(
        "VakintExpression::canonicalize",
        "rust",
        "use vakint::VakintExpression;"
    ),
    example!(
        "VakintExpression::tensor_reduce",
        "rust",
        "use vakint::VakintExpression;"
    ),
    example!(
        "VakintExpression::evaluate_integral",
        "rust",
        "use vakint::VakintExpression;"
    ),
    example!(
        "NumericalEvaluationResult",
        "rust",
        "use vakint::NumericalEvaluationResult;"
    ),
];

fn python_scope(request: &CatalogRequest, workspace_root: &Path, stub: &Path) -> Result<DocScope> {
    let stub_path = if stub.is_absolute() {
        stub.to_path_buf()
    } else {
        workspace_root.join(stub)
    };
    let source = fs::read_to_string(&stub_path)
        .wrap_err_with(|| format!("failed to read {}", stub_path.display()))?;
    let declarations = python_declarations(&source)?;
    ensure!(
        !declarations.is_empty(),
        "{} exports no Python declarations",
        request.component_id
    );
    for required in python_required_exports(&request.component_id)? {
        ensure!(
            declarations
                .iter()
                .any(|declaration| declaration.name == *required),
            "{} no longer exports supported entry {}",
            request.component_id,
            required
        );
    }
    let source_file = stub_path
        .strip_prefix(workspace_root)
        .unwrap_or(&stub_path)
        .to_string_lossy()
        .replace('\\', "/");

    let mut root = DocScope::new(
        &request.component_id,
        format!("{} Python API", request.product_title),
    );
    root.docs = Some(DocText::new(
        DocFormat::PythonDocstring,
        format!("Exports registered by {}.", request.package),
    ));
    root.summary = Some("Ordered Python module exports".to_owned());
    let mut exports = DocScope::new("exports", "Module exports");
    for declaration in declarations {
        let supported = python_export_is_supported(&request.component_id, &declaration.name)?;
        if supported {
            ensure!(
                !declaration.docs.trim().is_empty(),
                "supported export {} is missing a Python docstring",
                declaration.name
            );
        }
        let docs = if declaration.docs.trim().is_empty() {
            format!(
                "Reachable generated Python {} {}. This export is outside the curated support catalog.",
                match declaration.kind {
                    DocItemKind::PythonClass => "class",
                    _ => "function",
                },
                declaration.name
            )
        } else {
            declaration.docs
        };
        let summary = first_sentence(&docs);
        let mut item = DocItem::new(
            &declaration.name,
            &declaration.name,
            &declaration.name,
            declaration.kind,
        );
        item.supported = supported;
        // A Python binding export exists only when its isolated exporter feature
        // pair is enabled. Record that requirement per item as well as on the
        // component so an all-feature reference cannot imply default presence.
        item.required_features.clone_from(&request.features);
        item.docs = Some(DocText::new(DocFormat::PythonDocstring, docs));
        item.summary = Some(summary);
        item.signature = Some(declaration.signature);
        item.members = declaration.members;
        item.examples.push(DocExample::new(
            "Inspect the registered export",
            "python",
            format!(
                "from {} import {}\nhelp({})",
                request.module.as_deref().unwrap_or(&request.package),
                declaration.name,
                declaration.name
            ),
        ));
        item.source = Some(SourceLocation::new(
            format!("{}::{}", request.component_id, declaration.name),
            &source_file,
            declaration.line,
            1,
        ));
        exports.define_item(item)?;
    }
    root.define_scope(exports)?;
    Ok(root)
}

fn python_required_exports(component: &str) -> Result<&'static [&'static str]> {
    match component {
        "gammaloop-python" => Ok(&[
            "GammaLoopAPI",
            "atom_to_canonical_string",
            "evaluate_graph_overall_factor",
            "to_dots",
        ]),
        "spynso3" => Ok(&[
            "Representation",
            "Slot",
            "Tensor",
            "TensorIndices",
            "TensorStructure",
            "TensorNetwork",
            "TensorLibrary",
            "ExecutionMode",
            "SymbolicParallelism",
        ]),
        "linnet-py" | "idenso-community" | "vakint-community" => Ok(&[]),
        _ => bail!("unknown Python component {component}"),
    }
}

fn python_export_is_supported(component: &str, name: &str) -> Result<bool> {
    match component {
        "linnet-py" | "idenso-community" | "vakint-community" => Ok(true),
        "gammaloop-python" | "spynso3" => Ok(python_required_exports(component)?.contains(&name)),
        _ => bail!("unknown Python component {component}"),
    }
}

fn first_sentence(docs: &str) -> String {
    let text = docs.split_whitespace().collect::<Vec<_>>().join(" ");
    text.char_indices()
        .find(|(index, character)| {
            *character == '.'
                && text[*index + 1..]
                    .chars()
                    .next()
                    .is_none_or(char::is_whitespace)
        })
        .map(|(index, _)| text[..=index].to_owned())
        .unwrap_or(text)
}

pub fn write_catalog(path: &Path, catalog: &DocCatalog) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, serde_json::to_vec_pretty(catalog)?)
        .wrap_err_with(|| format!("failed to write {}", path.display()))
}

pub fn resolve_stub(workspace_root: &Path, component: &str) -> Option<PathBuf> {
    let checked = workspace_root
        .join("docs/api/python")
        .join(format!("{component}.pyi"));
    checked.is_file().then_some(checked)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_rust_component_has_a_strict_catalog() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        for component in [
            "gammalooprs",
            "gammaloop-api",
            "linnet",
            "spenso",
            "spenso-macros",
            "spenso-hep-lib",
            "idenso",
            "vakint",
        ] {
            let request = CatalogRequest {
                product_id: "test".to_owned(),
                product_title: "Test".to_owned(),
                component_id: component.to_owned(),
                package: component.to_owned(),
                component_title: component.to_owned(),
                version: "0.1.0".to_owned(),
                language: ApiLanguage::Rust,
                module: None,
                features: vec![],
            };
            let catalog = export_catalog(&request, root, None).unwrap();
            let docs = catalog.root.docs.expect("annotated component scope docs");
            assert_eq!(docs.format, DocFormat::TypstMarkup);
            assert!(
                docs.body.contains("#link(") || docs.body.contains("#emph["),
                "{component} scope should preserve raw Typst markup"
            );
            assert!(catalog.root.source.is_some());
        }
    }

    #[test]
    fn source_backed_descriptors_use_definition_owning_modules() {
        for (component, id, expected) in [
            (
                "linnet",
                "SuBitGraph",
                "linnet::half_edge::subgraph::subset::SuBitGraph",
            ),
            (
                "spenso",
                "DenseTensor",
                "spenso::tensors::data::dense::DenseTensor",
            ),
            (
                "spenso",
                "SparseTensor",
                "spenso::tensors::data::sparse::SparseTensor",
            ),
            ("idenso", "CookSettings", "idenso::cook::CookSettings"),
            ("idenso", "Cookable", "idenso::cook::Cookable"),
            (
                "vakint",
                "VakintExpression::canonicalize",
                "vakint::VakintExpression::canonicalize",
            ),
        ] {
            let item = annotated_items::for_component(component)
                .unwrap()
                .into_iter()
                .find(|item| item.id == id)
                .unwrap();
            let source = item.source.unwrap();
            assert_eq!(source.identifier, expected);
            assert!(source.line > 0);
            assert!(
                item.signature
                    .is_some_and(|signature| !signature.is_empty())
            );
        }
    }

    #[test]
    fn parses_python_exports_in_declared_order() {
        let stub = r#"
__all__ = [
    "run",
    "Engine",
]

class Engine:
    r"""
    Runs work. More detail.
    """

def run(value: int) -> int:
    r"""
    Return one value.
    """
"#;
        let declarations = python_declarations(stub).unwrap();
        assert_eq!(declarations[0].name, "run");
        assert_eq!(declarations[1].name, "Engine");
        assert_eq!(first_sentence(&declarations[1].docs), "Runs work.");
    }

    #[test]
    fn parses_one_line_python_docstrings() {
        let declarations =
            python_declarations("class Engine:\n    \"\"\"Runs work without setup.\"\"\"\n")
                .unwrap();
        assert_eq!(declarations[0].docs, "Runs work without setup.");
    }
}
