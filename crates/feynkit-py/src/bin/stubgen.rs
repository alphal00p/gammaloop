use std::fs;

const STUB_PATH: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/python/symbolica/community/feynkit/__init__.pyi"
);

fn main() -> pyo3_stub_gen::Result<()> {
    let info = feynkit_py::stub_info()?;
    let module = info
        .modules
        .get("symbolica.community.feynkit")
        .ok_or_else(|| {
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "FeynKit did not contribute a symbolica.community.feynkit stub module",
            )
        })?;
    fs::write(STUB_PATH, normalize_stub_source(&module.to_string()))?;
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
