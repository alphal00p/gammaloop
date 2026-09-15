use gammaloop_api::OneShot;

use color_eyre::Report;

fn main() -> Result<(), Report> {
    if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
        symbolica::activate_oem_license!("SO-227-gammaloop-2060.01.01-4IXECAD7GCKPJCGBEHWZPZSB3VOYGEWLSHPF54GQ2P2N6XRH4QWG4EHKY2MQQFRSP5SCHRTRNN5S7UJAMCLEKHSIQH3IXKKB3PAK6BI");
    }

    // Parse once with clap‑derive
    let cli = OneShot::parse_env_with_capture()?;

    if let Err(e) = cli.cli.run(cli.input_string) {
        eprintln!("{:?}", e);
    }

    Ok(())
}
