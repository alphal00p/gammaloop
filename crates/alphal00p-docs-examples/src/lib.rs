//! Compile and runtime acceptance checks for examples in the five product sites.
//!
//! The test module is generated from both neutral catalogs and authored Typst
//! manuals. This keeps checked snippets and rendered snippets at one semantic
//! boundary.

#[cfg(test)]
mod example_markers;

#[cfg(test)]
mod catalog_examples {
    #![allow(dead_code, unused_imports, unused_variables)]

    use std::{env, ffi::OsString, process::Command};

    include!(concat!(env!("OUT_DIR"), "/rust_catalog_examples.rs"));

    #[test]
    fn shell_catalog_examples_parse_without_execution() {
        let output = Command::new("bash")
            .arg("-n")
            .arg("-c")
            .arg(include_str!(concat!(
                env!("OUT_DIR"),
                "/shell_catalog_examples.sh"
            )))
            .output()
            .expect("failed to start bash for catalog syntax validation");
        assert!(
            output.status.success(),
            "shell documentation example validation failed:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    #[test]
    fn python_catalog_examples_compile_without_importing_products() {
        let python = python_interpreter();
        let output = Command::new(&python)
            .arg("-c")
            .arg(include_str!(concat!(
                env!("OUT_DIR"),
                "/python_catalog_examples.py"
            )))
            .output()
            .unwrap_or_else(|error| {
                panic!(
                    "failed to start Python interpreter {}: {error}",
                    python.to_string_lossy()
                )
            });
        assert!(
            output.status.success(),
            "Python documentation example validation failed:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("Python documentation examples"));
    }

    fn python_interpreter() -> OsString {
        env::var_os("PYTHON_BIN_PATH")
            .or_else(|| env::var_os("PYTHON"))
            .unwrap_or_else(|| OsString::from("python3"))
    }
}
