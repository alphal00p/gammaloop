//! Maintenance: produce the 5L route preparation add-on of a data checkout.
//!
//! ```text
//! cargo run --release -p vakint --example five_loop_route_preparation -- \
//!     DATA_ROOT OUTPUT_FILE [PREPARATION_WORKERS] [--check]
//! ```
//!
//! Loading the 5L saved application verifies its 8,179 momentum routes
//! exactly (about 140 s on one core). That work depends only on the saved
//! application in `DATA_ROOT/5L/ibps`, so RustRed can record it once into an
//! add-on file (`rustred_app::write_saved_application_route_addon`); the
//! loader then restores the identical routes in well under a second. The
//! artifact itself is not touched or regenerated.
//!
//! This command verifies every route of `DATA_ROOT/5L/ibps` (with
//! `PREPARATION_WORKERS` threads, default 1), writes the add-on to the new
//! file `OUTPUT_FILE`, and with `--check` loads the application once more
//! from the add-on and requires identical routes (owners, roots and complete
//! verified maps). It prints the 5L `manifest.json` entries that publish the
//! file (copy it to `5L/route-preparation.rrbin`):
//!
//! ```json
//! "route_preparation": {"path": "route-preparation.rrbin",
//!                       "format": "rustred.saved-application.route-preparation.v1"},
//! "files": {"route-preparation.rrbin": {"sha256": "...", "bytes": ...}}
//! ```
//!
//! Re-run it whenever the 5L `ibps/` artifact or RustRed's
//! `SAVED_VERIFIED_MAP_SEMANTICS` changes: a missing, corrupt or stale add-on
//! is rejected (one line on stderr) and the load verifies every route as
//! without it.
use std::path::PathBuf;

use rustred_app::{
    SAVED_ROUTE_PREPARATION_SCHEMA, SavedApplicationOptions, check_saved_application_route_addon,
    write_saved_application_route_addon,
};
use serde_json::json;
use sha2::{Digest, Sha256};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result<()> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let check = args.iter().any(|arg| arg == "--check");
    args.retain(|arg| arg != "--check");
    let (root, output, workers) = match args.as_slice() {
        [root, output] => (root, output, 1),
        [root, output, workers] => (root, output, workers.parse()?),
        _ => {
            return Err(
                "usage: five_loop_route_preparation DATA_ROOT OUTPUT_FILE [PREPARATION_WORKERS] [--check]"
                    .into(),
            );
        }
    };
    let application = PathBuf::from(root).join("5L").join("ibps");
    let output = PathBuf::from(output);
    let options = SavedApplicationOptions {
        preparation_workers: workers,
        ..SavedApplicationOptions::default()
    };
    let written = write_saved_application_route_addon::<15>(&application, &output, options)?;
    eprintln!("{}", json!({"event": "written", "report": written}));
    if check {
        let checked = check_saved_application_route_addon::<15>(&application, &output, options)?;
        eprintln!("{}", json!({"event": "checked", "report": checked}));
    }
    let bytes = std::fs::read(&output)?;
    let entries = json!({
        "route_preparation": {
            "path": "route-preparation.rrbin",
            "format": SAVED_ROUTE_PREPARATION_SCHEMA,
        },
        "files": {
            "route-preparation.rrbin": {
                "sha256": format!("{:x}", Sha256::digest(&bytes)),
                "bytes": bytes.len(),
            }
        }
    });
    println!("{}", serde_json::to_string_pretty(&entries)?);
    Ok(())
}
