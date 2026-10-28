use std::path::PathBuf;

#[cfg(feature = "watch")]
use alphal00p_docs_builder::WatchRequest;
use alphal00p_docs_builder::{BuildChannel, BuildRequest, RustdocCacheMode, SiteBuilder};
use clap::{Parser, Subcommand};
use eyre::Result;

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
        /// Exact Git tag required when channel is snapshot.
        #[arg(long)]
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
    /// Requires the `watch` Cargo feature; `just docs-watch` selects
    /// the feature and optimized watcher profile.
    #[cfg(feature = "watch")]
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

fn main() -> Result<()> {
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
        #[cfg(feature = "watch")]
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

#[cfg(test)]
mod tests {
    use std::env;
    use std::process::Command as ProcessCommand;

    use super::*;

    #[test]
    fn snapshot_tags_are_explicit_even_in_ci() {
        const CASE: &str = "ALPHAL00P_DOCS_CLI_TEST_REF";
        if env::var_os(CASE).is_none() {
            // Child environments avoid mutating the parallel test runner's globals.
            for reference in ["main", "feature/docs", "v0.3.4"] {
                let output = ProcessCommand::new(env::current_exe().unwrap())
                    .args(["--exact", "tests::snapshot_tags_are_explicit_even_in_ci"])
                    .env(CASE, reference)
                    .env("GITHUB_REF_NAME", reference)
                    .output()
                    .unwrap();
                assert!(output.status.success(), "reference {reference}: {output:?}");
                let stdout = String::from_utf8(output.stdout).unwrap();
                let summary = stdout.lines().find(|line| line.starts_with("test result:"));
                assert!(
                    summary
                        .is_some_and(|line| line
                            .starts_with("test result: ok. 1 passed; 0 failed; 0 ignored;")),
                    "expected exactly one subprocess test for {reference}: {stdout}"
                );
            }
            return;
        }

        for args in [
            vec!["docs", "build"],
            vec!["docs", "build", "--channel", "latest"],
            vec!["docs", "build", "--channel", "snapshot"],
        ] {
            let cli = Cli::try_parse_from(args).unwrap();
            let Command::Build { snapshot_tag, .. } = cli.command else {
                panic!("expected a build request");
            };
            assert_eq!(snapshot_tag, None);
        }

        let cli = Cli::try_parse_from([
            "docs",
            "build",
            "--channel",
            "snapshot",
            "--snapshot-tag",
            "linnet-v1.2.3",
        ])
        .unwrap();
        let Command::Build {
            channel,
            snapshot_tag,
            ..
        } = cli.command
        else {
            panic!("expected a build request");
        };
        assert_eq!(channel, BuildChannel::Snapshot);
        assert_eq!(snapshot_tag.as_deref(), Some("linnet-v1.2.3"));
    }
}
