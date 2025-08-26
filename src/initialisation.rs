use std::path::Path;

use color_eyre::Result;
use symbolica::activate_oem_license;
use tracing_subscriber::{reload, EnvFilter, Registry};

pub(crate) fn initialise() -> Result<()> {
    println!("{}", env!("CARGO_CRATE_NAME"));
    if option_env!("SYMBOLICA_OEM_LICENSE").is_some() {
        activate_oem_license!("SYMBOLICA_OEM_KEY_23177b25");
    };
    crate::initialize_reps();

    Ok(())
}

// pub(super) fn init_tracing() -> reload::Handle<EnvFilter, Registry> {
//     FILTER_HANDLE
//         .get_or_init(|| {
//             // 2) reloadable filter
//             let spec = std::env::var("RUST_LOG").unwrap_or_else(|_| "debug"..to_string());

//             let (filter_layer, handle) = reload::Layer::new(EnvFilter::new(spec));

//             let _ = std::fs::create_dir_all(dir.as_ref());

//             // e.g. "2025-08-26T21-05-33.123"
//             let ts = Local::now()
//                 .to_rfc3339_opts(SecondsFormat::Millis, true)
//                 .replace(':', "-")
//                 .replace('+', "-"); // keep it filename-friendly

//             let filename = format!("gammalog-{ts}.jsonl");
//             // One file per state (per process), no rotation
//             let file = RollingFileAppender::new(Rotation::NEVER, dir.as_ref(), &filename);
//             let (nb, guard) = NonBlockingBuilder::default()
//                 .buffered_lines_limit(200_000)
//                 .lossy(false)
//                 .finish(file);
//             LOG_GUARD.set(guard).ok();

//             let json = fmt::layer()
//                 .json()
//                 .flatten_event(true)
//                 .with_current_span(true)
//                 .with_span_list(true)
//                 .with_writer(nb);

//             // 4) pretty status to stderr, opt-in via `target="status"`
//             let status_layer = fmt::layer()
//                 .with_target(false)
//                 .event_format(StatusFmt::default())
//                 .with_writer(std::io::stderr)
//                 .with_filter(filter_fn(|m| m.target() == "status"));

//             tracing_subscriber::registry()
//                 .with(filter_layer)
//                 .with(json)
//                 .with(status_layer)
//                 .init();

//             handle
//         })
//         .clone()
// }
#[cfg(test)]
pub(crate) fn test_initialise() -> Result<()> {
    // env_logger::builder().is_test(true).try_init()?;
    crate::initialize_reps();

    Ok(())
}
