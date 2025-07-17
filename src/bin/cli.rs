use _gammaloop::cli::{Cli, Commands};
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
    // --- logging boilerplate (unchanged) ---
    {
        let mut log_builder = env_logger::builder();
        if std::env::var("RUST_LOG").is_err() {
            log_builder.filter_level(log::LevelFilter::Info);
        }
        log_builder
            .format_target(true)
            .format_timestamp(Some(env_logger::TimestampPrecision::Millis))
            .init();
    }

    // Parse once with clap‑derive
    let cli = Cli::parse();
    let settings = cli.get_settings()?;

    match &cli.command {
        Some(Commands::Repl) => {
            // Re‑enter interactive mode. Forward any global flags if you like.
            _gammaloop::cli::repl::start(settings)?;
        }
        _ => {
            // Single‑shot execution path (old behaviour)
            cli.run()?;
        }
    }

    Ok(())
}
