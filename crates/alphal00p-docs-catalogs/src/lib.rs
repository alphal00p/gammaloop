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
    match request.language {
        ApiLanguage::Python => validate_python_catalog(&catalog)?,
        _ => catalog.validate_supported()?,
    }
    Ok(catalog)
}

fn validate_python_catalog(catalog: &DocCatalog) -> Result<()> {
    fn validate_scope(scope: &DocScope, supported: &mut usize) -> Result<()> {
        for item in scope.items.values().filter(|item| item.supported) {
            *supported += 1;
            ensure!(
                item.docs
                    .as_ref()
                    .is_some_and(|docs| !docs.body.trim().is_empty()),
                "supported Python export {} is missing docs",
                item.id
            );
            ensure!(
                item.summary
                    .as_ref()
                    .is_some_and(|summary| !summary.trim().is_empty()),
                "supported Python export {} is missing a summary",
                item.id
            );
            ensure!(
                item.signature
                    .as_ref()
                    .is_some_and(|signature| !signature.trim().is_empty()),
                "supported Python export {} is missing a signature",
                item.id
            );
            ensure!(
                item.source.as_ref().is_some_and(|source| {
                    !source.identifier.trim().is_empty() && !source.file.trim().is_empty()
                }),
                "supported Python export {} is missing a source location",
                item.id
            );
        }
        for child in scope.scopes.values() {
            validate_scope(child, supported)?;
        }
        Ok(())
    }

    catalog.validate()?;
    let mut supported = 0;
    validate_scope(&catalog.root, &mut supported)?;
    ensure!(supported > 0, "Python catalog has no supported exports");
    Ok(())
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
        ("spenso", "NetworkParse" | "ParamTensor" | "ParseSettings" | "SymbolicParallelism") => {
            &["shadowing"]
        }
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
        "GammaLoopContextContainer",
        "rust",
        "use gammalooprs::GammaLoopContextContainer;"
    ),
    example!(
        "HasIntegrand",
        "rust",
        "fn accepts_integrand<I: gammalooprs::integrands::HasIntegrand>(_integrand: &I) {}"
    ),
    example!(
        "Integrand",
        "rust",
        "use gammalooprs::integrands::Integrand;"
    ),
    example!(
        "GlobalSettings",
        "rust",
        "let global = gammalooprs::settings::GlobalSettings::default();"
    ),
    example!(
        "RuntimeSettings",
        "rust",
        "let runtime = gammalooprs::settings::RuntimeSettings::default();"
    ),
    example!(
        "SingleSampleEvaluationResult",
        "rust",
        "use gammalooprs::integrands::evaluation::SingleSampleEvaluationResult;"
    ),
    example!(
        "BatchSampleEvaluationResult",
        "rust",
        "use gammalooprs::integrands::evaluation::BatchSampleEvaluationResult;"
    ),
    example!(
        "PreciseSingleSampleEvaluationResult",
        "rust",
        "use gammalooprs::integrands::evaluation::PreciseSingleSampleEvaluationResult;"
    ),
    example!(
        "PreciseBatchSampleEvaluationResult",
        "rust",
        "use gammalooprs::integrands::evaluation::PreciseBatchSampleEvaluationResult;"
    ),
    example!(
        "HistogramSnapshot",
        "rust",
        "use gammalooprs::observables::HistogramSnapshot;"
    ),
    example!(
        "HistogramSnapshot::merge",
        "rust",
        "let left = gammalooprs::observables::HistogramAccumulatorState::default().snapshot(); let right = gammalooprs::observables::HistogramAccumulatorState::default().snapshot(); let merged = left.merge(&right)?;"
    ),
    example!(
        "HistogramSnapshot::rebin",
        "rust",
        "let accumulator = gammalooprs::observables::HistogramAccumulatorState::continuous(\"energy\".into(), \"weighted events\".into(), gammalooprs::observables::ObservablePhase::Real, gammalooprs::observables::ObservableValueTransform::Identity, 0.0, 4.0, false, false, 4); let coarser = accumulator.snapshot().rebin(2)?; assert_eq!(coarser.bins.len(), 2);"
    ),
    example!(
        "HistogramAccumulatorState",
        "rust",
        "let accumulator = gammalooprs::observables::HistogramAccumulatorState::default();"
    ),
    example!(
        "HistogramAccumulatorState::snapshot",
        "rust",
        "let accumulator = gammalooprs::observables::HistogramAccumulatorState::default(); let snapshot = accumulator.snapshot();"
    ),
    example!(
        "set_interrupt_handler",
        "rust",
        "gammalooprs::set_interrupt_handler();"
    ),
    example!(
        "request_interrupt",
        "rust",
        "gammalooprs::request_interrupt();"
    ),
    example!(
        "is_interrupt_requested",
        "rust",
        "let interrupted = gammalooprs::is_interrupt_requested();"
    ),
    example!(
        "clear_interrupt_request",
        "rust",
        "gammalooprs::clear_interrupt_request();"
    ),
];

const GAMMALOOP_API: &[RustExampleSpec] = &[
    example!(
        "StateLoadOption",
        "rust",
        "let options = gammaloop_api::StateLoadOption::default();"
    ),
    example!("LoadedState", "rust", "use gammaloop_api::LoadedState;"),
    example!(
        "StateLoadOption::load",
        "rust",
        "let loaded = gammaloop_api::StateLoadOption::default().load()?;"
    ),
    example!(
        "LoadedState::cli_session",
        "rust",
        "fn open_session(loaded: &mut gammaloop_api::LoadedState) -> gammaloop_api::session::CliSession<'_> { loaded.cli_session() }"
    ),
    example!(
        "CliSession",
        "rust",
        "use gammaloop_api::session::CliSession;"
    ),
    example!(
        "CliSession::execute_command",
        "rust",
        "fn execute(session: &mut gammaloop_api::session::CliSession<'_>, command: gammaloop_api::state::CommandHistory) -> eyre::Result<()> { session.execute_command(command)?; Ok(()) }"
    ),
    example!(
        "CommandHistory",
        "rust",
        "use gammaloop_api::state::CommandHistory;"
    ),
    example!(
        "CommandHistory::from_raw_string",
        "rust",
        "let command = gammaloop_api::state::CommandHistory::from_raw_string(\"display processes\")?;"
    ),
    example!(
        "RunHistory",
        "rust",
        "use gammaloop_api::state::RunHistory;"
    ),
    example!(
        "RunHistory::load",
        "rust",
        "let history = gammaloop_api::state::RunHistory::load(\"run.toml\")?;"
    ),
    example!("State", "rust", "use gammaloop_api::state::State;"),
    example!(
        "State::get_integrand_info",
        "rust",
        "fn inspect(state: &gammaloop_api::state::State) -> eyre::Result<gammaloop_api::integrand_info::IntegrandInfo> { state.get_integrand_info(None, None) }"
    ),
    example!(
        "IntegrandInfo",
        "rust",
        "use gammaloop_api::integrand_info::IntegrandInfo;"
    ),
    example!(
        "ImportModel",
        "rust",
        "use gammaloop_api::commands::import::model::ImportModel;"
    ),
    example!(
        "ImportModel::run",
        "rust",
        "fn import(request: &gammaloop_api::commands::import::model::ImportModel, state: &mut gammaloop_api::state::State) -> eyre::Result<()> { request.run(state) }"
    ),
    example!("Generate", "rust", "use gammaloop_api::commands::Generate;"),
    example!(
        "Generate::run",
        "rust",
        "fn generate(request: &gammaloop_api::commands::Generate, state: &mut gammaloop_api::state::State, global: &gammalooprs::settings::GlobalSettings, runtime: &gammalooprs::settings::RuntimeSettings) -> eyre::Result<()> { request.run(state, \"compiled\", false, global, runtime) }"
    ),
    example!(
        "Integrate",
        "rust",
        "use gammaloop_api::commands::Integrate;"
    ),
    example!(
        "Integrate::run",
        "rust",
        "fn integrate(request: &gammaloop_api::commands::Integrate, state: &mut gammaloop_api::state::State, settings: &gammaloop_api::CLISettings) -> eyre::Result<gammaloop_api::commands::IntegrationOutput> { request.run(state, settings) }"
    ),
    example!(
        "IntegrationOutput",
        "rust",
        "use gammaloop_api::commands::IntegrationOutput;"
    ),
    example!(
        "EvaluateSamples",
        "rust",
        "use gammaloop_api::commands::evaluate_samples::EvaluateSamples;"
    ),
    example!(
        "evaluate_sample",
        "rust",
        "fn evaluate(state: &mut gammaloop_api::state::State, request: &gammaloop_api::commands::evaluate_samples::EvaluateSamples<'_>) -> eyre::Result<()> { let _result = gammaloop_api::commands::evaluate_samples::evaluate_sample(state, request)?; Ok(()) }"
    ),
    example!(
        "evaluate_samples",
        "rust",
        "fn evaluate(state: &mut gammaloop_api::state::State, request: &gammaloop_api::commands::evaluate_samples::EvaluateSamples<'_>) -> eyre::Result<()> { let _batch = gammaloop_api::commands::evaluate_samples::evaluate_samples(state, request)?; Ok(()) }"
    ),
    example!(
        "EvaluateSamplesPrecise",
        "rust",
        "use gammaloop_api::commands::evaluate_samples::EvaluateSamplesPrecise;"
    ),
    example!(
        "evaluate_sample_precise",
        "rust",
        "fn evaluate(state: &mut gammaloop_api::state::State, request: &gammaloop_api::commands::evaluate_samples::EvaluateSamplesPrecise<'_>) -> eyre::Result<()> { let _result = gammaloop_api::commands::evaluate_samples::evaluate_sample_precise(state, request)?; Ok(()) }"
    ),
    example!(
        "evaluate_samples_precise",
        "rust",
        "fn evaluate(state: &mut gammaloop_api::state::State, request: &gammaloop_api::commands::evaluate_samples::EvaluateSamplesPrecise<'_>) -> eyre::Result<()> { let _batch = gammaloop_api::commands::evaluate_samples::evaluate_samples_precise(state, request)?; Ok(()) }"
    ),
    example!(
        "StateFolderKind",
        "rust",
        "use gammaloop_api::state::StateFolderKind;"
    ),
    example!(
        "classify_state_folder",
        "rust",
        "let kind = gammaloop_api::state::classify_state_folder(std::path::Path::new(\"./state\"))?;"
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
    example!(
        "HedgePair",
        "rust",
        "use linnet::half_edge::involution::HedgePair;"
    ),
    example!("Flow", "rust", "use linnet::half_edge::involution::Flow;"),
    example!(
        "Orientation",
        "rust",
        "use linnet::half_edge::involution::Orientation;"
    ),
    example!(
        "SubSetLike",
        "rust",
        "use linnet::half_edge::subgraph::SubSetLike;"
    ),
    example!(
        "SubSetOps",
        "rust",
        "use linnet::half_edge::subgraph::SubSetOps;"
    ),
    example!(
        "Inclusion",
        "rust",
        "use linnet::half_edge::subgraph::Inclusion;"
    ),
    example!(
        "BaseSubgraph",
        "rust",
        "use linnet::half_edge::subgraph::BaseSubgraph;"
    ),
    example!(
        "SubGraphLike",
        "rust",
        "use linnet::half_edge::subgraph::SubGraphLike;"
    ),
    example!(
        "NodeStorage",
        "rust",
        "use linnet::half_edge::nodestore::NodeStorage;"
    ),
    example!(
        "NodeStorageOps",
        "rust",
        "use linnet::half_edge::nodestore::NodeStorageOps;"
    ),
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
        "ContractionStrategy",
        "rust",
        "use spenso::network::ContractionStrategy;"
    ),
    example!(
        "ExecutionStrategy",
        "rust",
        "use spenso::network::ExecutionStrategy;"
    ),
    example!(
        "ExecutionMode",
        "rust",
        "use spenso::network::ExecutionMode;"
    ),
    example!(
        "ExecutionResult",
        "rust",
        "use spenso::network::ExecutionResult;"
    ),
    example!(
        "TensorNetworkError",
        "rust",
        "use spenso::network::TensorNetworkError;"
    ),
    example!(
        "ParseSettings",
        "rust",
        "let settings = spenso::network::parsing::ParseSettings::default();"
    ),
    example!(
        "NetworkParse",
        "rust",
        "use spenso::network::parsing::NetworkParse;"
    ),
    example!(
        "SymbolicParallelism",
        "rust",
        "use spenso::symbolic_parallelism::SymbolicParallelism;"
    ),
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
            "AdditionalWeight",
            "BatchEvaluationResult",
            "CutInfo",
            "DotExportSettings",
            "EvaluationResult",
            "Event",
            "EventGroup",
            "FourMomentum",
            "GammaLoopAPI",
            "HistogramAccumulator",
            "HistogramBin",
            "HistogramSnapshot",
            "HistogramStats",
            "IntegrandActiveThresholdCut",
            "IntegrandCut",
            "IntegrandCutThreshold",
            "IntegrandGraph",
            "IntegrandGraphGroup",
            "IntegrandInfo",
            "IntegrandLoopMomentumBasis",
            "IntegrandOrientation",
            "IntegrandThresholdEsurface",
            "SampleEvaluationResult",
            "SettingsValue",
            "StabilityResult",
            "atom_to_canonical_string",
            "evaluate_graph_overall_factor",
            "to_dots",
        ]),
        "spynso3" => Ok(&[
            "CompiledTensorEvaluator",
            "LibraryTensor",
            "Representation",
            "Slot",
            "Tensor",
            "TensorEvaluator",
            "TensorIndices",
            "TensorStructure",
            "TensorNetwork",
            "TensorLibrary",
            "TensorName",
            "ExecutionMode",
            "SymbolicParallelism",
            "set_symbolica_rayon_enabled",
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
                "linnet",
                "HedgePair",
                "linnet::half_edge::involution::HedgePair",
            ),
            (
                "linnet",
                "SubSetLike",
                "linnet::half_edge::subgraph::SubSetLike",
            ),
            (
                "linnet",
                "NodeStorageOps",
                "linnet::half_edge::nodestore::NodeStorageOps",
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
            (
                "spenso",
                "ContractionStrategy",
                "spenso::network::contract::ContractionStrategy",
            ),
            (
                "spenso",
                "NetworkParse",
                "spenso::network::parsing::NetworkParse",
            ),
            (
                "spenso",
                "SymbolicParallelism",
                "spenso::symbolic_parallelism::SymbolicParallelism",
            ),
            ("idenso", "CookSettings", "idenso::cook::CookSettings"),
            ("idenso", "Cookable", "idenso::cook::Cookable"),
            (
                "vakint",
                "VakintExpression::canonicalize",
                "vakint::VakintExpression::canonicalize",
            ),
            (
                "gammaloop-api",
                "ImportModel",
                "gammaloop_api::commands::import::model::ImportModel",
            ),
            (
                "gammaloop-api",
                "Generate::run",
                "gammaloop_api::commands::generate::Generate::run",
            ),
            (
                "gammaloop-api",
                "IntegrationOutput",
                "gammaloop_api::commands::integrate::IntegrationOutput",
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
    fn gammaloop_supported_surfaces_have_regression_guarded_coverage() {
        assert_eq!(rust_examples("gammalooprs").unwrap().len(), 20);
        let rust = rust_examples("gammaloop-api").unwrap();
        assert_eq!(rust.len(), 30);
        for entry in [
            "ImportModel",
            "ImportModel::run",
            "Generate",
            "Generate::run",
            "Integrate",
            "Integrate::run",
            "IntegrationOutput",
        ] {
            assert!(
                rust.iter().any(|example| example.id == entry),
                "missing supported Rust entry {entry}"
            );
        }

        let python = python_required_exports("gammaloop-python").unwrap();
        assert_eq!(python.len(), 28);
        for entry in [
            "GammaLoopAPI",
            "IntegrandInfo",
            "HistogramAccumulator",
            "EventGroup",
            "SettingsValue",
        ] {
            assert!(
                python.contains(&entry),
                "missing supported Python entry {entry}"
            );
        }
    }

    #[test]
    fn gammalooprs_curated_direct_members_are_documented() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let request = CatalogRequest {
            product_id: "gammaloop".to_owned(),
            product_title: "GammaLoop".to_owned(),
            component_id: "gammalooprs".to_owned(),
            package: "gammalooprs".to_owned(),
            component_title: "GammaLoop Rust core".to_owned(),
            version: "0.1.0".to_owned(),
            language: ApiLanguage::Rust,
            module: None,
            features: vec![],
        };
        let catalog = export_catalog(&request, root, None).unwrap();
        let supported = &catalog.root.scopes["supported"];
        assert_eq!(supported.items.len(), 20);

        let members = supported
            .items
            .values()
            .flat_map(|item| item.members.iter().map(move |member| (item, member)))
            .collect::<Vec<_>>();
        assert_eq!(members.len(), 74);
        let undocumented = members
            .iter()
            .filter(|(_, member)| {
                member
                    .docs
                    .as_ref()
                    .is_none_or(|docs| docs.body.trim().is_empty())
            })
            .map(|(item, member)| format!("{}.{}", item.id, member.name))
            .collect::<Vec<_>>();
        assert!(
            undocumented.is_empty(),
            "undocumented curated gammalooprs members: {undocumented:?}"
        );
    }

    #[test]
    fn gammaloop_lifecycle_source_docs_are_complete() {
        let items = annotated_items::for_component("gammaloop-api").unwrap();
        let import_model = items.iter().find(|item| item.id == "ImportModel").unwrap();
        assert_eq!(
            import_model.docs.as_ref().map(|docs| docs.body.as_str()),
            Some(
                "Imports a built-in JSON model, a JSON model file, or a UFO model directory into the active state."
            )
        );

        for (item_id, member_name, expected) in [
            (
                "Generate",
                "mode",
                "Selects cross-section, amplitude, or existing-process generation; when\nomitted, generates all integrands in the active state.",
            ),
            (
                "Integrate",
                "workspace_path",
                "Workspace directory used for resumable integration state and persisted results",
            ),
            (
                "IntegrationOutput",
                "result",
                "Final estimates and convergence data for every selected process/integrand slot.",
            ),
            (
                "IntegrationOutput",
                "observables",
                "Observable snapshots keyed by each slot's `process@integrand` identifier.",
            ),
            (
                "IntegrationOutput",
                "workspace_path",
                "Workspace that stores resumable integration state, results, and observable snapshots.",
            ),
        ] {
            let item = items.iter().find(|item| item.id == item_id).unwrap();
            let member = item
                .members
                .iter()
                .find(|member| member.name == member_name)
                .unwrap();
            assert_eq!(
                member.docs.as_ref().map(|docs| docs.body.as_str()),
                Some(expected),
                "unexpected source docs for {item_id}.{member_name}"
            );
        }
    }

    #[test]
    fn spenso_supported_rust_surface_covers_execution_and_parsing() {
        let rust = rust_examples("spenso").unwrap();
        assert_eq!(rust.len(), 14);
        for entry in [
            "ContractionStrategy",
            "ExecutionStrategy",
            "ExecutionMode",
            "ExecutionResult",
            "TensorNetworkError",
            "ParseSettings",
            "NetworkParse",
            "SymbolicParallelism",
        ] {
            assert!(
                rust.iter().any(|example| example.id == entry),
                "missing supported Rust entry {entry}"
            );
        }
        for entry in ["ParseSettings", "NetworkParse", "SymbolicParallelism"] {
            assert_eq!(rust_item_required_features("spenso", entry), &["shadowing"]);
        }
    }

    #[test]
    fn linnet_supported_rust_surface_covers_foundational_graph_contracts() {
        let rust = rust_examples("linnet").unwrap();
        assert_eq!(rust.len(), 15);
        for entry in [
            "HedgePair",
            "Flow",
            "Orientation",
            "SubSetLike",
            "SubSetOps",
            "Inclusion",
            "BaseSubgraph",
            "SubGraphLike",
            "NodeStorage",
            "NodeStorageOps",
        ] {
            assert!(
                rust.iter().any(|example| example.id == entry),
                "missing supported Rust entry {entry}"
            );
        }
    }

    #[test]
    fn every_supported_gammaloop_python_export_has_docs() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let source = fs::read_to_string(root.join("docs/api/python/gammaloop-python.pyi"))
            .expect("checked-in GammaLoop Python stub");
        let declarations = python_declarations(&source).unwrap();
        for entry in python_required_exports("gammaloop-python").unwrap() {
            let declaration = declarations
                .iter()
                .find(|declaration| declaration.name == *entry)
                .unwrap_or_else(|| panic!("missing supported Python declaration {entry}"));
            assert!(
                !declaration.docs.trim().is_empty(),
                "supported Python declaration {entry} has no docs"
            );
        }
    }

    #[test]
    fn spynso_supported_surface_covers_the_documented_workflow_types() {
        let required = python_required_exports("spynso3").unwrap();
        assert_eq!(required.len(), 14);
        for entry in [
            "CompiledTensorEvaluator",
            "LibraryTensor",
            "TensorEvaluator",
            "TensorName",
            "set_symbolica_rayon_enabled",
        ] {
            assert!(
                required.contains(&entry),
                "missing supported Python entry {entry}"
            );
        }

        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let source = fs::read_to_string(root.join("docs/api/python/spynso3.pyi"))
            .expect("checked-in spynso3 Python stub");
        let declarations = python_declarations(&source).unwrap();
        let promoted = declarations
            .iter()
            .filter(|declaration| {
                matches!(
                    declaration.name.as_str(),
                    "CompiledTensorEvaluator"
                        | "LibraryTensor"
                        | "TensorEvaluator"
                        | "TensorName"
                        | "set_symbolica_rayon_enabled"
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(promoted.len(), 5);
        assert!(
            promoted
                .iter()
                .all(|declaration| !declaration.docs.trim().is_empty())
        );
        let members = promoted
            .iter()
            .flat_map(|declaration| &declaration.members)
            .collect::<Vec<_>>();
        assert_eq!(members.len(), 31);
        let overload_groups = members
            .iter()
            .filter(|member| {
                member
                    .members
                    .iter()
                    .any(|member| member.kind == alphal00p_docs_schema::DocMemberKind::Overload)
            })
            .count();
        let overloads = members
            .iter()
            .flat_map(|member| &member.members)
            .filter(|member| member.kind == alphal00p_docs_schema::DocMemberKind::Overload)
            .collect::<Vec<_>>();
        assert_eq!(members.len() - overload_groups + overloads.len(), 35);
        assert_eq!(
            members
                .iter()
                .filter(|member| {
                    !member
                        .members
                        .iter()
                        .any(|member| member.kind == alphal00p_docs_schema::DocMemberKind::Overload)
                        && member.docs.is_some()
                })
                .count()
                + overloads
                    .iter()
                    .filter(|member| member.docs.is_some())
                    .count(),
            28
        );
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

    #[test]
    fn python_catalogs_keep_source_examples_without_introspection_scaffolds() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let request = CatalogRequest {
            product_id: "test".to_owned(),
            product_title: "Test".to_owned(),
            component_id: "idenso-community".to_owned(),
            package: "idenso".to_owned(),
            component_title: "Idenso Python".to_owned(),
            version: "0.1.0".to_owned(),
            language: ApiLanguage::Python,
            module: Some("symbolica.community.idenso".to_owned()),
            features: vec![],
        };
        let catalog = export_catalog(
            &request,
            root,
            Some(&root.join("docs/api/python/idenso-community.pyi")),
        )
        .unwrap();
        let items = catalog.root.scopes["exports"]
            .items
            .values()
            .collect::<Vec<_>>();

        assert!(items.iter().all(|item| item.examples.is_empty()));
        assert!(items.iter().any(|item| {
            item.docs
                .as_ref()
                .is_some_and(|docs| docs.body.contains("Examples"))
        }));
        assert!(items.iter().all(|item| {
            item.examples.iter().all(|example| {
                !example.title.contains("Introspection scaffold") && !example.code.contains("help(")
            })
        }));
    }
}
