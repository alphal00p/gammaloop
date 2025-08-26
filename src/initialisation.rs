use color_eyre::Result;
use symbolica::activate_oem_license;

pub(crate) fn initialise() -> Result<()> {
    println!("{}", env!("CARGO_CRATE_NAME"));
    if option_env!("SYMBOLICA_OEM_LICENSE").is_some() {
        activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
    };
    crate::initialize_reps();

    Ok(())
}
#[cfg(test)]
pub(crate) fn test_initialise() -> Result<()> {
    // env_logger::builder().is_test(true).try_init()?;
    crate::initialize_reps();

    Ok(())
}
