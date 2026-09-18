use std::path::{Path, PathBuf};

use alphal00p_docs_catalogs::{generated::write_or_check_snapshot, vakint_reference};
use clap::Parser;
use eyre::Result;

#[derive(Debug, Parser)]
#[command(about = "Export Vakint's runtime topology and dependency metadata")]
struct Cli {
    #[arg(long, default_value = ".")]
    workspace_root: PathBuf,
    #[arg(long, default_value = "docs/api/generated/vakint-reference.json")]
    output: PathBuf,
    #[arg(long)]
    check: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let output = absolute_from(&cli.workspace_root, &cli.output);
    write_or_check_snapshot(&vakint_reference::export()?, &output, cli.check)
}

fn absolute_from(root: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        root.join(path)
    }
}
