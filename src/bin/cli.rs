use _gammaloop::cli::{Cli, Commands};
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
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
