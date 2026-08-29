use std::path::PathBuf;

use alphal00p_docs_builder::{
    BuildChannel, BuildRequest, RustdocCacheMode, SiteBuilder, WatchRequest,
};
use clap::{Parser, Subcommand};

#[derive(Debug, Parser)]
#[command(about = "Build and validate the alphal00p multi-product documentation")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Validate product configuration and generated catalog inputs.
    Check,
    /// Build one product site or all product sites.
    Build {
        /// Product id or "all".
        #[arg(long, default_value = "all")]
        product: String,
        /// Mutable latest output or an immutable tagged snapshot.
        #[arg(long, value_enum, default_value_t = BuildChannel::Latest)]
        channel: BuildChannel,
        /// Root of the ready-to-deploy Pages artifact.
        #[arg(long, default_value = "target/alphal00p-docs")]
        output: PathBuf,
        /// Exact Git tag used when channel is snapshot.
        #[arg(long, env = "GITHUB_REF_NAME")]
        snapshot_tag: Option<String>,
        /// Skip exhaustive cargo rustdoc sidecars.
        #[arg(long)]
        skip_rustdoc: bool,
        /// Skip Typst invocation and emit only generated metadata.
        #[arg(long, hide = true)]
        skip_typst: bool,
        /// Persistent Cargo target used by a live-watch session.
        #[arg(long, hide = true)]
        rustdoc_target_root: Option<PathBuf>,
        /// Directory receiving Typst dependency files from the renderer.
        #[arg(long, hide = true)]
        dependency_output: Option<PathBuf>,
    },
    /// Continuously rebuild with retained Typst state and browser reloads.
    ///
    /// Requires the `persistent-typst` Cargo feature; `just docs-watch` selects
    /// the feature and optimized watcher profile.
    Watch {
        /// Product id or "all".
        #[arg(long, default_value = "all")]
        product: String,
        /// Output updated after successful builds when --no-serve is used.
        #[arg(long, default_value = "target/alphal00p-docs-watch")]
        output: PathBuf,
        /// Local HTTP port used by the live server.
        #[arg(long, default_value_t = 8000)]
        port: u16,
        /// Open the served documentation after the first successful build.
        #[arg(long)]
        open: bool,
        /// Continuously write to --output instead of starting the live server.
        #[arg(long)]
        no_serve: bool,
        /// Skip exhaustive cargo rustdoc sidecars while watching.
        #[arg(long)]
        skip_rustdoc: bool,
    },
}

fn main() -> eyre::Result<()> {
    let cli = Cli::parse();
    let builder = SiteBuilder::discover()?;

    match cli.command {
        Command::Check => builder.check(),
        Command::Build {
            product,
            channel,
            output,
            snapshot_tag,
            skip_rustdoc,
            skip_typst,
            rustdoc_target_root,
            dependency_output,
        } => builder.build(BuildRequest {
            product,
            channel,
            output,
            snapshot_tag,
            include_rustdoc: !skip_rustdoc,
            include_typst: !skip_typst,
            rustdoc_target_root,
            rustdoc_cache: RustdocCacheMode::Disabled,
            dependency_output,
        }),
        Command::Watch {
            product,
            output,
            port,
            open,
            no_serve,
            skip_rustdoc,
        } => builder.watch(WatchRequest {
            product,
            output,
            port,
            open,
            serve: !no_serve,
            include_rustdoc: !skip_rustdoc,
        }),
    }
}
