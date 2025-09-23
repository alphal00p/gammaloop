use _gammaloop::{OneShot, Repl};

use color_eyre::Report;

fn main() -> Result<(), Report> {
    // Parse once with clap‑derive
    let cli = OneShot::parse_env_with_capture()?;

    if let Err(e) = cli.cli.run(cli.input_string) {
        eprintln!("{:?}", e);
    }

    Ok(())
}
