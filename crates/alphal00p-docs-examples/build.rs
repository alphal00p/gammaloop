use std::{
    collections::BTreeSet,
    env,
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

use alphal00p_docs_catalogs::{CatalogRequest, export_catalog};
use alphal00p_docs_schema::{ApiLanguage, DocCatalog, DocItem, DocScope};
use eyre::{Context, ContextCompat, Result, ensure};
use serde::Deserialize;

#[path = "src/example_markers.rs"]
mod example_markers;

use example_markers::{FencedBlock, fenced_blocks, named_blocks, require_named_block};

#[derive(Debug, Deserialize)]
struct Registry {
    product: Vec<Product>,
}

#[derive(Debug, Deserialize)]
struct Product {
    id: String,
    title: String,
    #[serde(default)]
    pages: Vec<Page>,
    #[serde(default)]
    rust_components: Vec<Component>,
    #[serde(default)]
    python_components: Vec<Component>,
}

#[derive(Debug, Deserialize)]
struct Page {
    source: PathBuf,
    route: String,
}

#[derive(Debug, Deserialize)]
struct ExampleRegistry {
    schema: u32,
    example: Vec<RegisteredExample>,
}

#[derive(Debug, Deserialize)]
struct RegisteredExample {
    id: String,
    marker: String,
    product: String,
    task: String,
    audience: String,
    source: PathBuf,
    route: String,
    prerequisites: Vec<String>,
    expected: String,
    verification: String,
    tier: String,
    owner: String,
    test_command: String,
}

#[derive(Clone, Debug, Deserialize)]
struct Component {
    id: String,
    package: String,
    version_source: PathBuf,
    module: Option<String>,
    #[serde(default)]
    features: Vec<String>,
}

#[derive(Default)]
struct ManualExamples {
    rust: Vec<(String, String, String)>,
    python: Vec<(String, String, String)>,
    shell: Vec<(String, String)>,
    data_blocks: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RustExampleMode {
    Run,
    RunIsolated,
    CompileOnly,
    ShellSyntax,
    RustdocIgnored,
}

fn main() -> Result<()> {
    let crate_root =
        PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").context("missing CARGO_MANIFEST_DIR")?);
    let workspace_root = crate_root
        .parent()
        .and_then(Path::parent)
        .context("documentation examples crate is not two levels below the workspace root")?;
    let registry_path = workspace_root.join("docs/products/registry.toml");
    println!("cargo:rerun-if-changed={}", registry_path.display());
    let registry: Registry = toml::from_str(
        &fs::read_to_string(&registry_path)
            .wrap_err_with(|| format!("failed to read {}", registry_path.display()))?,
    )
    .wrap_err_with(|| format!("failed to parse {}", registry_path.display()))?;
    let examples_path = workspace_root.join("docs/examples.toml");
    println!("cargo:rerun-if-changed={}", examples_path.display());
    let example_registry: ExampleRegistry = toml::from_str(
        &fs::read_to_string(&examples_path)
            .wrap_err_with(|| format!("failed to read {}", examples_path.display()))?,
    )
    .wrap_err_with(|| format!("failed to parse {}", examples_path.display()))?;
    validate_example_registry(workspace_root, &registry, &example_registry)?;

    let output = PathBuf::from(env::var_os("OUT_DIR").context("missing OUT_DIR")?);
    let mut rust_source = String::from(
        "static ISOLATED_RUST_EXAMPLE: std::sync::Mutex<()> = std::sync::Mutex::new(());\n",
    );
    let mut shell_source = String::from("#!/usr/bin/env bash\nset -euo pipefail\n");
    let mut python_source = String::from("EXAMPLES = [\n");
    let mut rust_count = 0;
    let mut run_count = 0;
    let mut shell_count = 0;
    let mut python_count = 0;
    let mut rust_products = BTreeSet::new();
    let mut python_products = BTreeSet::new();

    for product in &registry.product {
        for component in &product.rust_components {
            watch_component_inputs(workspace_root, component, None);
            let catalog =
                component_catalog(workspace_root, product, component, ApiLanguage::Rust, None)?;
            for (item_path, item) in supported_items(&catalog) {
                for (example_index, example) in item.examples.iter().enumerate() {
                    let mode = rust_example_mode(&component.id, &item.id, &example.language)?;
                    let case = format!("{}::{item_path}::{example_index}", component.id);
                    match mode {
                        RustExampleMode::Run
                        | RustExampleMode::RunIsolated
                        | RustExampleMode::CompileOnly => {
                            rust_products.insert(product.id.as_str());
                            let function = rust_identifier(&case);
                            writeln!(
                                rust_source,
                                "#[allow(dead_code)]\nfn {function}() -> eyre::Result<()> {{\n{}\n    Ok(())\n}}",
                                indent(&example.code, 4)
                            )?;
                            rust_count += 1;
                            match mode {
                                RustExampleMode::Run => writeln!(
                                    rust_source,
                                    "#[test]\nfn run_{function}() -> eyre::Result<()> {{ {function}() }}"
                                )?,
                                RustExampleMode::RunIsolated => writeln!(
                                    rust_source,
                                    "#[test]\nfn run_{function}() -> eyre::Result<()> {{\n    const MARKER: &str = \"{function}\";\n    if std::env::var(\"ALPHAL00P_DOCS_ISOLATED_EXAMPLE\").as_deref() == Ok(MARKER) {{\n        return {function}();\n    }}\n    let _guard = ISOLATED_RUST_EXAMPLE.lock().expect(\"isolated example mutex poisoned\");\n    let status = std::process::Command::new(std::env::current_exe()?)\n        .arg(\"--exact\")\n        .arg(\"catalog_examples::run_{function}\")\n        .env(\"ALPHAL00P_DOCS_ISOLATED_EXAMPLE\", MARKER)\n        .status()?;\n    eyre::ensure!(status.success(), \"isolated documentation example {case} failed\");\n    Ok(())\n}}"
                                )?,
                                RustExampleMode::CompileOnly => {}
                                RustExampleMode::ShellSyntax | RustExampleMode::RustdocIgnored => {
                                    unreachable!()
                                }
                            }
                            if matches!(mode, RustExampleMode::Run | RustExampleMode::RunIsolated) {
                                run_count += 1;
                            }
                        }
                        RustExampleMode::ShellSyntax => {
                            writeln!(shell_source, "# {case}\n{}", example.code)?;
                            shell_count += 1;
                        }
                        RustExampleMode::RustdocIgnored => {}
                    }
                }
            }
        }

        for component in &product.python_components {
            let stub = workspace_root
                .join("docs/api/python")
                .join(format!("{}.pyi", component.id));
            watch_component_inputs(workspace_root, component, Some(&stub));
            let catalog = component_catalog(
                workspace_root,
                product,
                component,
                ApiLanguage::Python,
                Some(&stub),
            )?;
            for (item_path, item) in supported_items(&catalog) {
                for (example_index, example) in item.examples.iter().enumerate() {
                    ensure!(
                        example.language == "python",
                        "Python catalog example {}::{item_path}::{example_index} uses unsupported language {}",
                        component.id,
                        example.language
                    );
                    python_products.insert(product.id.as_str());
                    let case = format!("{}::{item_path}::{example_index}", component.id);
                    writeln!(
                        python_source,
                        "    ({}, {}),",
                        serde_json::to_string(&case)?,
                        serde_json::to_string(&example.code)?
                    )?;
                    python_count += 1;
                }
            }
        }
    }

    let manual = manual_examples(workspace_root, &registry)?;
    for (case, code, mode) in &manual.rust {
        let function = rust_identifier(case);
        if mode == "syntax" {
            rust_count += 1;
            continue;
        }
        if syn::parse_file(code).is_ok() {
            let run_test = if mode == "run" {
                run_count += 1;
                format!(
                    "\n    #[test]\n    fn run_documented_main() -> eyre::Result<()> {{\n        const MARKER: &str = \"{function}\";\n        if std::env::var(\"ALPHAL00P_DOCS_ISOLATED_EXAMPLE\").as_deref() == Ok(MARKER) {{\n            main();\n            return Ok(());\n        }}\n        let _guard = super::ISOLATED_RUST_EXAMPLE.lock().expect(\"isolated example mutex poisoned\");\n        let status = std::process::Command::new(std::env::current_exe()?)\n            .arg(\"--exact\")\n            .arg(\"catalog_examples::{function}::run_documented_main\")\n            .env(\"ALPHAL00P_DOCS_ISOLATED_EXAMPLE\", MARKER)\n            .status()?;\n        eyre::ensure!(status.success(), \"isolated documentation example {function} failed\");\n        Ok(())\n    }}"
                )
            } else {
                String::new()
            };
            writeln!(
                rust_source,
                "#[allow(dead_code, unused_imports)]\nmod {function} {{\n{}{}\n}}",
                indent(code, 4),
                run_test
            )?;
        } else {
            ensure!(
                mode != "run",
                "runnable Rust documentation example {case} must be a complete source file"
            );
            writeln!(
                rust_source,
                "#[allow(dead_code)]\nfn {function}() -> eyre::Result<()> {{\n{}\n    Ok(())\n}}",
                indent(code, 4)
            )?;
        }
        rust_count += 1;
    }
    for (product, case, code) in &manual.python {
        python_products.insert(product);
        writeln!(
            python_source,
            "    ({}, {}),",
            serde_json::to_string(case)?,
            serde_json::to_string(code)?
        )?;
        python_count += 1;
    }
    for (case, code) in &manual.shell {
        writeln!(shell_source, "# {case}\n{code}")?;
        shell_count += 1;
    }

    let expected_products = registry
        .product
        .iter()
        .map(|product| product.id.as_str())
        .collect::<BTreeSet<_>>();
    ensure!(
        rust_products == expected_products,
        "Rust examples do not cover every registered product: expected {expected_products:?}, found {rust_products:?}"
    );
    ensure!(
        python_products == expected_products,
        "Python examples do not cover every registered product: expected {expected_products:?}, found {python_products:?}"
    );
    ensure!(rust_count > 0, "no Rust catalog examples were generated");
    ensure!(
        run_count > 0,
        "no safe Rust catalog examples were selected to run"
    );
    ensure!(shell_count > 0, "no shell catalog examples were generated");
    ensure!(
        python_count > 0,
        "no Python catalog examples were generated"
    );
    ensure!(
        !manual.rust.is_empty() && !manual.python.is_empty() && !manual.shell.is_empty(),
        "manual examples must include Rust, Python, and shell blocks"
    );
    ensure!(
        manual.data_blocks > 0,
        "manual examples must include validated TOML or text blocks"
    );

    python_source.push_str(
        r#"]

for name, source in EXAMPLES:
    compile(source, name, "exec")

print(f"compiled {len(EXAMPLES)} Python documentation examples")
"#,
    );

    fs::write(output.join("rust_catalog_examples.rs"), rust_source)?;
    fs::write(output.join("shell_catalog_examples.sh"), shell_source)?;
    fs::write(output.join("python_catalog_examples.py"), python_source)?;
    Ok(())
}

fn validate_example_registry(
    workspace_root: &Path,
    products: &Registry,
    registry: &ExampleRegistry,
) -> Result<()> {
    ensure!(registry.schema == 2, "unsupported example registry schema");
    let page_sources = products
        .product
        .iter()
        .flat_map(|product| &product.pages)
        .map(|page| page.source.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|path| {
            fs::read_to_string(workspace_root.join(&path))
                .map(|source| (path.clone(), source))
                .wrap_err_with(|| format!("failed to read {}", path.display()))
        })
        .collect::<Result<Vec<_>>>()?;
    let named = named_blocks(&page_sources)?;
    let mut ids = BTreeSet::new();
    let mut markers = BTreeSet::new();
    let mut coverage = BTreeSet::new();
    for example in &registry.example {
        ensure!(
            ids.insert(&example.id),
            "duplicate example id {}",
            example.id
        );
        ensure!(
            example.marker == example.id,
            "example {} marker must equal its example id",
            example.id
        );
        ensure!(
            markers.insert(&example.marker),
            "duplicate example marker {}",
            example.marker
        );
        let product = products
            .product
            .iter()
            .find(|product| product.id == example.product)
            .wrap_err_with(|| format!("example {} has unknown product", example.id))?;
        let page = product
            .pages
            .iter()
            .find(|page| page.source == example.source)
            .wrap_err_with(|| format!("example {} has an unregistered source", example.id))?;
        ensure!(
            page.route == example.route,
            "example {} route does not match its registered page",
            example.id
        );
        ensure!(
            matches!(example.audience.as_str(), "beginner" | "real-value"),
            "example {} has unsupported audience {}",
            example.id,
            example.audience
        );
        ensure!(
            matches!(example.verification.as_str(), "run" | "compile" | "syntax"),
            "example {} has unsupported verification mode {}",
            example.id,
            example.verification
        );
        ensure!(
            matches!(example.tier.as_str(), "pull-request" | "scheduled"),
            "example {} has unsupported verification tier {}",
            example.id,
            example.tier
        );
        ensure!(
            !example.task.trim().is_empty()
                && !example.expected.trim().is_empty()
                && !example.owner.trim().is_empty()
                && !example.test_command.trim().is_empty()
                && !example.prerequisites.is_empty()
                && example
                    .prerequisites
                    .iter()
                    .all(|prerequisite| !prerequisite.trim().is_empty()),
            "example {} has incomplete maintenance metadata",
            example.id
        );
        require_named_block(
            &named,
            &example.marker,
            &example.source,
            &example.verification,
        )?;
        coverage.insert((example.product.as_str(), example.audience.as_str()));
    }
    let unregistered = named
        .keys()
        .filter(|marker| !markers.contains(*marker))
        .collect::<Vec<_>>();
    ensure!(
        unregistered.is_empty(),
        "named documentation example markers are missing registry records: {unregistered:?}"
    );
    let expected = products
        .product
        .iter()
        .flat_map(|product| {
            ["beginner", "real-value"]
                .into_iter()
                .map(move |audience| (product.id.as_str(), audience))
        })
        .collect::<BTreeSet<_>>();
    ensure!(
        coverage == expected,
        "canonical examples must cover beginner and real-value journeys for every product"
    );
    Ok(())
}

fn manual_examples(workspace_root: &Path, registry: &Registry) -> Result<ManualExamples> {
    let mut examples = ManualExamples::default();
    for product in &registry.product {
        let content_root = workspace_root
            .join("docs/products")
            .join(&product.id)
            .join("content");
        println!("cargo:rerun-if-changed={}", content_root.display());
        let mut files = fs::read_dir(&content_root)
            .wrap_err_with(|| format!("failed to read {}", content_root.display()))?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?;
        files.retain(|path| path.extension().is_some_and(|extension| extension == "typ"));
        files.sort();

        for path in files {
            let source = fs::read_to_string(&path)
                .wrap_err_with(|| format!("failed to read {}", path.display()))?;
            let relative = path
                .strip_prefix(workspace_root)
                .expect("manual content is below the workspace");
            let is_tutorial = relative
                .file_name()
                .is_some_and(|name| name == "tutorial.typ");
            for FencedBlock {
                language,
                marker,
                line,
                code,
            } in fenced_blocks(&source, relative)?
            {
                let case = format!("manual::{}::{}:{line}", product.id, relative.display());
                ensure!(
                    !is_tutorial || marker.is_some(),
                    "tutorial example {case} needs an immediately preceding // docs-example: MODE marker"
                );
                let mode = marker
                    .as_ref()
                    .map(|marker| marker.mode.as_str())
                    .unwrap_or(match language.as_str() {
                        "rust" | "python" => "compile",
                        _ => "syntax",
                    });
                match (language.as_str(), mode) {
                    ("rust", mode @ ("run" | "compile" | "syntax")) => {
                        if matches!(mode, "run" | "syntax") {
                            syn::parse_file(&code).wrap_err_with(|| {
                                format!("invalid Rust documentation example {case}")
                            })?;
                        }
                        examples.rust.push((case, code, mode.to_owned()));
                    }
                    ("python", "compile") => {
                        examples.python.push((product.id.clone(), case, code));
                    }
                    ("sh" | "bash" | "shell", "syntax") => {
                        examples.shell.push((case, code));
                    }
                    ("toml", "syntax") => {
                        toml::from_str::<toml::Value>(&code).wrap_err_with(|| {
                            format!("invalid TOML documentation example {case}")
                        })?;
                        examples.data_blocks += 1;
                    }
                    ("text", "syntax") => {
                        ensure!(!code.trim().is_empty(), "empty text example {case}");
                        examples.data_blocks += 1;
                    }
                    ("rs" | "c" | "form" | "typ", "syntax") => {
                        // The Idenso specification quotes Rust-like Symbolica
                        // patterns plus upstream C and FORM source fragments;
                        // authored Typst pages can likewise show partial import
                        // fragments. They are evidence, not standalone programs.
                        ensure!(!code.trim().is_empty(), "empty source fragment {case}");
                        examples.data_blocks += 1;
                    }
                    _ => ensure!(
                        false,
                        "unsupported documentation example language {language:?} with verification mode {mode:?} in {case}"
                    ),
                }
            }
        }
    }
    Ok(examples)
}

fn component_catalog(
    workspace_root: &Path,
    product: &Product,
    component: &Component,
    language: ApiLanguage,
    stub: Option<&Path>,
) -> Result<DocCatalog> {
    let request = CatalogRequest {
        product_id: product.id.clone(),
        product_title: product.title.clone(),
        component_id: component.id.clone(),
        package: component.package.clone(),
        component_title: component.id.clone(),
        version: manifest_version(&workspace_root.join(&component.version_source))?,
        language,
        module: component.module.clone(),
        features: component.features.clone(),
    };
    export_catalog(&request, workspace_root, stub)
        .wrap_err_with(|| format!("failed to export {} examples", component.id))
}

fn manifest_version(path: &Path) -> Result<String> {
    let source = fs::read_to_string(path)
        .wrap_err_with(|| format!("failed to read version source {}", path.display()))?;
    let manifest: toml::Value = toml::from_str(&source)
        .wrap_err_with(|| format!("failed to parse version source {}", path.display()))?;
    ["package", "project"]
        .iter()
        .find_map(|section| {
            manifest
                .get(section)
                .and_then(|table| table.get("version"))
                .and_then(toml::Value::as_str)
        })
        .map(str::to_owned)
        .wrap_err_with(|| format!("{} has no package/project version", path.display()))
}

fn watch_component_inputs(workspace_root: &Path, component: &Component, stub: Option<&Path>) {
    println!(
        "cargo:rerun-if-changed={}",
        workspace_root.join(&component.version_source).display()
    );
    if let Some(stub) = stub {
        println!("cargo:rerun-if-changed={}", stub.display());
    }
}

fn supported_items(catalog: &DocCatalog) -> Vec<(String, &DocItem)> {
    fn visit<'a>(scope: &'a DocScope, prefix: &str, items: &mut Vec<(String, &'a DocItem)>) {
        for item in scope.items.values().filter(|item| item.supported) {
            items.push((format!("{prefix}::{}", item.id), item));
        }
        for child in scope.scopes.values() {
            visit(child, &format!("{prefix}::{}", child.id), items);
        }
    }

    let mut items = Vec::new();
    visit(&catalog.root, &catalog.root.id, &mut items);
    items
}

fn rust_example_mode(component: &str, item: &str, language: &str) -> Result<RustExampleMode> {
    if language == "ignore" {
        return Ok(RustExampleMode::RustdocIgnored);
    }
    if (component, item, language) == ("gammaloop-api", "gammaloop", "console") {
        return Ok(RustExampleMode::ShellSyntax);
    }
    ensure!(
        language.split(',').next() == Some("rust"),
        "Rust catalog example {component}::{item} uses unsupported language {language}"
    );
    if (component, item) == ("spenso", "Network") {
        return Ok(RustExampleMode::RunIsolated);
    }
    if matches!(
        (component, item),
        (
            "gammalooprs",
            "set_interrupt_handler"
                | "request_interrupt"
                | "is_interrupt_requested"
                | "clear_interrupt_request"
        ) | (
            "gammaloop-api",
            "StateLoadOption::load"
                | "LoadedState::cli_session"
                | "CliSession::execute_command"
                | "RunHistory::load"
                | "State::get_integrand_info"
                | "ImportModel::run"
                | "Generate::run"
                | "Integrate::run"
                | "evaluate_sample"
                | "evaluate_samples"
                | "evaluate_sample_precise"
                | "evaluate_samples_precise"
                | "classify_state_folder"
        ) | ("vakint", "vakint_parse" | "vakint_symbol" | "Vakint")
            | (
                "idenso",
                "bis"
                    | "cof"
                    | "coad"
                    | "gamma"
                    | "gamma0"
                    | "gamma5"
                    | "epsilon"
                    | "color_t"
                    | "color_f"
                    | "representations::initialize"
            )
    ) {
        return Ok(RustExampleMode::CompileOnly);
    }
    ensure!(
        matches!(
            (component, item),
            (
                "gammalooprs",
                "GammaLoopContext"
                    | "HasModel"
                    | "GammaLoopContextContainer"
                    | "HasIntegrand"
                    | "Integrand"
                    | "GlobalSettings"
                    | "RuntimeSettings"
                    | "SingleSampleEvaluationResult"
                    | "BatchSampleEvaluationResult"
                    | "PreciseSingleSampleEvaluationResult"
                    | "PreciseBatchSampleEvaluationResult"
                    | "HistogramSnapshot"
                    | "HistogramSnapshot::merge"
                    | "HistogramSnapshot::rebin"
                    | "HistogramAccumulatorState"
                    | "HistogramAccumulatorState::snapshot"
            ) | (
                "gammaloop-api",
                "StateLoadOption"
                    | "LoadedState"
                    | "CliSession"
                    | "CommandHistory"
                    | "CommandHistory::from_raw_string"
                    | "RunHistory"
                    | "State"
                    | "IntegrandInfo"
                    | "ImportModel"
                    | "Generate"
                    | "Integrate"
                    | "IntegrationOutput"
                    | "EvaluateSamples"
                    | "EvaluateSamplesPrecise"
                    | "StateFolderKind"
                    | "CLISettings"
            ) | (
                "linnet",
                "HedgeGraph"
                    | "HedgeGraphBuilder"
                    | "SuBitGraph"
                    | "SimpleTraversalTree"
                    | "DotGraph"
                    | "HedgePair"
                    | "Flow"
                    | "Orientation"
                    | "SubSetLike"
                    | "SubSetOps"
                    | "Inclusion"
                    | "BaseSubgraph"
                    | "SubGraphLike"
                    | "NodeStorage"
                    | "NodeStorageOps"
            ) | (
                "spenso",
                "TensorStructure"
                    | "DenseTensor"
                    | "SparseTensor"
                    | "Contract"
                    | "Network"
                    | "ContractionStrategy"
                    | "ExecutionStrategy"
                    | "ExecutionMode"
                    | "ExecutionResult"
                    | "TensorNetworkError"
                    | "ParseSettings"
                    | "NetworkParse"
                    | "SymbolicParallelism"
                    | "ParamTensor"
            ) | ("spenso-macros", "SimpleRepresentation")
                | (
                    "spenso-hep-lib",
                    "hep_lib" | "gamma_data_dirac" | "su3_generator_data"
                )
                | (
                    "idenso",
                    "IndexTooling"
                        | "CookSettings"
                        | "Cookable"
                        | "SelectiveExpand"
                        | "bis"
                        | "cof"
                        | "coad"
                        | "gamma"
                        | "gamma0"
                        | "gamma5"
                        | "epsilon"
                        | "color_t"
                        | "color_f"
                        | "representations::initialize"
                )
                | (
                    "vakint",
                    "vakint_parse"
                        | "vakint_symbol"
                        | "Vakint"
                        | "VakintSettings"
                        | "VakintExpression"
                        | "EvaluationOrder"
                        | "VakintExpression::canonicalize"
                        | "VakintExpression::tensor_reduce"
                        | "VakintExpression::evaluate_integral"
                        | "NumericalEvaluationResult"
                )
        ),
        "new Rust catalog example {component}::{item} must be classified as safe-to-run or compile-only"
    );
    Ok(RustExampleMode::Run)
}

fn rust_identifier(case: &str) -> String {
    let mut identifier = String::from("catalog_example_");
    let mut previous_was_separator = false;
    for character in case.chars() {
        if character.is_ascii_alphanumeric() {
            identifier.push(character.to_ascii_lowercase());
            previous_was_separator = false;
        } else if !previous_was_separator {
            identifier.push('_');
            previous_was_separator = true;
        }
    }
    identifier
}

fn indent(source: &str, spaces: usize) -> String {
    let padding = " ".repeat(spaces);
    source
        .lines()
        .map(|line| format!("{padding}{line}"))
        .collect::<Vec<_>>()
        .join("\n")
}
