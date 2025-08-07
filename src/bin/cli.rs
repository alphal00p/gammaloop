use _gammaloop::cli::Cli;
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
    // Parse once with clap‑derive
    let cli = Cli::parse();

    cli.run();

    Ok(())
}
