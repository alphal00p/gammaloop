use _gammaloop::Cli;
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
    // Parse once with clap‑derive
    let cli = Cli::parse();

    if let Err(e) = cli.run() {
        eprintln!("{:?}", e);
    }

    Ok(())
}
