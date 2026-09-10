//! Neutral, checked-in snapshots produced by feature-isolated runtime exporters.

use std::{fs, path::Path};

use eyre::{Context, Result, bail};
use serde::Serialize;

pub use alphal00p_docs_schema::generated::*;

pub fn write_or_check_snapshot<T: Serialize>(value: &T, output: &Path, check: bool) -> Result<()> {
    let mut rendered = serde_json::to_vec_pretty(value)?;
    rendered.push(b'\n');
    if check {
        let checked_in = fs::read(output)
            .with_context(|| format!("missing generated snapshot {}", output.display()))?;
        if checked_in != rendered {
            bail!(
                "generated snapshot drifted: {}; regenerate it with the corresponding exporter",
                output.display()
            );
        }
        return Ok(());
    }

    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(output, rendered).with_context(|| format!("failed to write {}", output.display()))?;
    Ok(())
}
