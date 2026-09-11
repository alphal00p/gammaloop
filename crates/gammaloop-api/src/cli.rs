use gammaloop_api::OneShot;

use color_eyre::Report;

fn main() -> Result<(), Report> {
    if option_env!("NO_SYMBOLICA_OEM_LICENSE").is_none() {
        symbolica::activate_oem_license!("SL1.AqlJHAAAAAAAAAAA42dhbW1hbG9vcOIuQQB_MJT0iMEh7ZfmQd1dgxLLkd5e8NDT9N9eJ-QsbhDqxpkIFjJ_ZCPGcWt7L9EgYJZFHkiB9oupQdvArwU");
    }

    // Parse once with clap‑derive
    let cli = OneShot::parse_env_with_capture()?;

    if let Err(e) = cli.cli.run(cli.input_string) {
        eprintln!("{:?}", e);
    }

    Ok(())
}
