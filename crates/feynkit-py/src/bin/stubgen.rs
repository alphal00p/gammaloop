use std::env;
use std::fs;
use std::path::PathBuf;

fn main() -> pyo3_stub_gen::Result<()> {
    let root = env::var_os("GAMMALOOP_WORKSPACE_ROOT")
        .or_else(|| env::var_os("NEXTEST_WORKSPACE_ROOT"))
        .map(|root| PathBuf::from(root).join("crates/feynkit-py"))
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")));
    let info = feynkit_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.hepkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.hepkit stub module",
            )
        })?;
    fs::write(
        root.join("python/symbolica/community/hepkit/__init__.pyi"),
        normalize_stub_source(&feynkit_py::FeynkitModule::stub_source(module)),
    )?;
    Ok(())
}

fn normalize_stub_source(source: &str) -> String {
    let mut lines = source
        .lines()
        .map(|line| line.trim_end_matches([' ', '\t']))
        .collect::<Vec<_>>();
    while lines.last().is_some_and(|line| line.is_empty()) {
        lines.pop();
    }
    lines.join("\n") + "\n"
}

#[cfg(test)]
mod tests {
    use super::normalize_stub_source;

    #[test]
    fn normalizes_trailing_whitespace_and_end_of_file() {
        assert_eq!(
            normalize_stub_source("class Example:  \n    ...\t\n   \n\n"),
            "class Example:\n    ...\n"
        );
    }
}
