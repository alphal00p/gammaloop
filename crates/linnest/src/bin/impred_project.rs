//! Line-delimited host access to Linnest's native constraint owner.
use linnest::impred_projection::ProjectionService;
use serde_json::{json, Value};
use std::io::{self, BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut output = io::stdout().lock();
    let mut service = ProjectionService::default();
    for line in io::stdin().lock().lines() {
        let response = serde_json::from_str::<Value>(&line?)
            .map_err(|error| error.to_string())
            .and_then(|request| service.handle(request))
            .unwrap_or_else(|error| json!({"ok":false,"error":error}));
        writeln!(output, "{response}")?;
        output.flush()?;
    }
    Ok(())
}
