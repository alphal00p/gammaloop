use std::path::PathBuf;

use alphal00p_docs_catalogs::{CatalogRequest, export_catalog, write_catalog};
use alphal00p_docs_schema::ApiLanguage;
use clap::{Parser, ValueEnum};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Language {
    Rust,
    Python,
}

impl From<Language> for ApiLanguage {
    fn from(language: Language) -> Self {
        match language {
            Language::Rust => Self::Rust,
            Language::Python => Self::Python,
        }
    }
}

#[derive(Debug, Parser)]
#[command(about = "Export one isolated alphal00p documentation catalog")]
struct Cli {
    #[arg(long)]
    product: String,
    #[arg(long)]
    product_title: String,
    #[arg(long)]
    component: String,
    #[arg(long)]
    package: String,
    #[arg(long)]
    component_title: String,
    #[arg(long)]
    version: String,
    #[arg(long, value_enum)]
    language: Language,
    #[arg(long)]
    module: Option<String>,
    /// Feature used to compile the documented component; repeat as needed.
    #[arg(long = "feature")]
    features: Vec<String>,
    #[arg(long, default_value = ".")]
    workspace_root: PathBuf,
    #[arg(long)]
    stub: Option<PathBuf>,
    /// Stable repository-relative source path recorded for a staged stub.
    #[arg(long)]
    source_file: Option<String>,
    #[arg(long)]
    output: PathBuf,
}

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();
    let request = CatalogRequest {
        product_id: cli.product,
        product_title: cli.product_title,
        component_id: cli.component,
        package: cli.package,
        component_title: cli.component_title,
        version: cli.version,
        language: cli.language.into(),
        module: cli.module,
        features: cli.features,
    };
    let mut catalog = export_catalog(&request, &cli.workspace_root, cli.stub.as_deref())?;
    if let Some(source_file) = cli.source_file {
        rewrite_source_files(&mut catalog.root, &source_file);
    }
    write_catalog(&cli.output, &catalog)
}

fn rewrite_source_files(scope: &mut alphal00p_docs_schema::DocScope, source_file: &str) {
    if let Some(source) = &mut scope.source {
        source.file = source_file.to_owned();
    }
    for item in scope.items.values_mut() {
        if let Some(source) = &mut item.source {
            source.file = source_file.to_owned();
        }
    }
    for child in scope.scopes.values_mut() {
        rewrite_source_files(child, source_file);
    }
}
