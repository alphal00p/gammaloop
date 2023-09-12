use _gammaloop::cli_functions::cli;
use color_eyre::Report;
use env_logger;
use std::env;

fn main() -> Result<(), Report> {
    let log_builder = &mut env_logger::builder();
    if let Err(_) = env::var("RUST_LOG") {
        log_builder.filter_level(log::LevelFilter::Info);
    }
    log_builder
        .format_target(true)
        .format_timestamp(Some(env_logger::TimestampPrecision::Millis))
        .init();
    cli(&env::args().collect::<Vec<_>>())
}
