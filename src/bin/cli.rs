use _gammaloop::cli::{state::State, Cli};
use clap::Parser;
use color_eyre::Report;

fn main() -> Result<(), Report> {
    // Parse once with clap‑derive
    let cli = Cli::parse();
    let mut settings = cli.get_settings()?;
    // let mut state = State::load(&cli.state_file,);

    cli.run_with_settings(&mut settings);

    Ok(())
}
