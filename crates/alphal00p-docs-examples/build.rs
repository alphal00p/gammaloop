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

#[derive(Debug, Deserialize)]
struct Registry {
    product: Vec<Product>,
}

#[derive(Debug, Deserialize)]
struct Product {
    id: String,
    title: String,
    #[serde(default)]
    rust_components: Vec<Component>,
    #[serde(default)]
    python_components: Vec<Component>,
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
    rust: Vec<(String, String)>,
    python: Vec<(String, String)>,
    shell: Vec<(String, String)>,
    data_blocks: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RustExampleMode {
    Run,
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

    let output = PathBuf::from(env::var_os("OUT_DIR").context("missing OUT_DIR")?);
    let mut rust_source = String::new();
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
                        RustExampleMode::Run | RustExampleMode::CompileOnly => {
                            rust_products.insert(product.id.as_str());
                            let function = rust_identifier(&case);
                            writeln!(
                                rust_source,
                                "#[allow(dead_code)]\nfn {function}() -> eyre::Result<()> {{\n{}\n    Ok(())\n}}",
                                indent(&example.code, 4)
                            )?;
                            rust_count += 1;
                            if mode == RustExampleMode::Run {
                                writeln!(
                                    rust_source,
                                    "#[test]\nfn run_{function}() -> eyre::Result<()> {{ {function}() }}"
                                )?;
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
    for (case, code) in &manual.rust {
        let function = rust_identifier(case);
        if syn::parse_file(code).is_ok() {
            writeln!(
                rust_source,
                "#[allow(dead_code, unused_imports)]\nmod {function} {{\n{}\n}}",
                indent(code, 4)
            )?;
        } else {
            writeln!(
                rust_source,
                "#[allow(dead_code)]\nfn {function}() -> eyre::Result<()> {{\n{}\n    Ok(())\n}}",
                indent(code, 4)
            )?;
        }
        rust_count += 1;
    }
    for (case, code) in &manual.python {
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
            for (language, line, code) in fenced_blocks(&source, relative)? {
                let case = format!("manual::{}::{}:{line}", product.id, relative.display());
                match language.as_str() {
                    "rust" => examples.rust.push((case, code)),
                    "python" => examples.python.push((case, code)),
                    "sh" | "bash" | "shell" => examples.shell.push((case, code)),
                    "toml" => {
                        toml::from_str::<toml::Value>(&code).wrap_err_with(|| {
                            format!("invalid TOML documentation example {case}")
                        })?;
                        examples.data_blocks += 1;
                    }
                    "text" => {
                        ensure!(!code.trim().is_empty(), "empty text example {case}");
                        examples.data_blocks += 1;
                    }
                    _ => ensure!(
                        false,
                        "unsupported fenced example language {language:?} in {case}"
                    ),
                }
            }
        }
    }
    Ok(examples)
}

fn fenced_blocks(source: &str, path: &Path) -> Result<Vec<(String, usize, String)>> {
    let lines = source.lines().collect::<Vec<_>>();
    let mut blocks = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        let Some(language) = lines[index].strip_prefix("```") else {
            index += 1;
            continue;
        };
        let language = language.trim().to_owned();
        ensure!(
            !language.is_empty(),
            "untyped fenced block at {}:{}",
            path.display(),
            index + 1
        );
        let start = index + 1;
        index += 1;
        let mut code = Vec::new();
        while index < lines.len() && lines[index] != "```" {
            code.push(lines[index]);
            index += 1;
        }
        ensure!(
            index < lines.len(),
            "unterminated {language} block at {}:{}",
            path.display(),
            start
        );
        blocks.push((language, start, code.join("\n")));
        index += 1;
    }
    Ok(blocks)
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
    if matches!(
        (component, item),
        ("gammalooprs", "set_interrupt_handler")
            | ("gammaloop-api", "StateLoadOption::load")
            | ("gammaloop-api", "LoadedState::cli_session")
            | ("vakint", "vakint_parse" | "vakint_symbol" | "Vakint")
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
            ("gammalooprs", "GammaLoopContext" | "HasModel")
                | ("gammaloop-api", "StateLoadOption" | "CLISettings")
                | (
                    "linnet",
                    "HedgeGraph"
                        | "HedgeGraphBuilder"
                        | "SuBitGraph"
                        | "SimpleTraversalTree"
                        | "DotGraph"
                )
                | (
                    "spenso",
                    "TensorStructure"
                        | "DenseTensor"
                        | "SparseTensor"
                        | "Contract"
                        | "Network"
                        | "ParamTensor"
                )
                | ("spenso-macros", "SimpleRepresentation")
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
