use _gammaloop::Cli;
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
    // Parse once with clap‑derive
    let cli = Cli::parse_env_with_capture()?;

    if let Err(e) = cli.cli.run(cli.input_string) {
        eprintln!("{:?}", e);
    }

    Ok(())
}
