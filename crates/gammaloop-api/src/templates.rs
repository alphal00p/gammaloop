use alphal00p_persistence::OutputPath;
use color_eyre::Result;
use include_dir::Dir;
use std::path::Path;

const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";

/// Embedded drawing assets written beneath a caller-selected output root.
///
/// ```no_run
/// use gammaloop_api::templates::Assets;
/// # use color_eyre::Result;
/// # fn main() -> Result<()> {
/// Assets::extract_templates("drawing-export")?;
/// Assets::extract_justfile("drawing-export")?;
/// assert!(Assets::template_paths().any(|name| name == "figure.typ"));
/// # Ok(()) }
/// ```
pub struct Assets;

include!(concat!(env!("OUT_DIR"), "/embedded_assets.rs"));

impl Assets {
    /// Get all app template paths written under drawings/templates.
    pub fn template_paths() -> impl Iterator<Item = String> {
        GAMMALOOP_TEMPLATE_ASSETS
            .files()
            .map(|file| file.path().to_string_lossy().into_owned())
    }

    /// Extract all drawing templates to drawings/templates relative to target_dir.
    /// See [`Assets`] for an example.
    ///
    /// # Errors
    /// Rejects invalid embedded names and symlinks inside the selected root, and
    /// propagates directory creation and checked file-write failures.
    pub fn extract_templates<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let target = OutputPath::new(target_dir.as_ref())?
            .join("drawings")?
            .join("templates")?;
        target.ensure_directory()?;
        Self::extract_package(&GAMMALOOP_TEMPLATE_ASSETS, &target)?;
        Self::extract_package(
            &LINNEST_PACKAGE_ASSETS,
            &target.join_relative(Path::new(LINNEST_PACKAGE_DIR))?,
        )?;
        Self::extract_package(
            &KURVST_PACKAGE_ASSETS,
            &target.join_relative(Path::new(KURVST_PACKAGE_DIR))?,
        )?;

        Ok(())
    }

    /// Extract the Justfile to the target directory. See [`Assets`] for an example.
    ///
    /// # Errors
    /// Rejects symlinks inside the selected root and propagates creation and write failures.
    pub fn extract_justfile<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let target = OutputPath::new(target_dir.as_ref())?;
        let justfile = target.join("justfile")?;
        target.ensure_directory()?;
        justfile.write(include_bytes!("../../../assets/embedded/drawing/justfile"))?;
        Ok(())
    }

    fn extract_package(assets: &Dir<'_>, target: &OutputPath) -> Result<()> {
        target.ensure_directory()?;
        for file in assets.files() {
            let target_path = target.join_relative(file.path())?;
            if let Some(parent) = file
                .path()
                .parent()
                .filter(|path| !path.as_os_str().is_empty())
            {
                target.join_relative(parent)?.ensure_directory()?;
            }
            target_path.write(file.contents())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gammalooprs::model::ModelGammaLoopExt;
    use std::fs;

    #[test]
    fn extract_templates_uses_canonical_package_layout() -> Result<()> {
        let tempdir = tempfile::tempdir()?;

        Assets::extract_templates(tempdir.path())?;

        let templates = tempdir.path().join("drawings/templates");
        assert!(templates.join("figure.typ").is_file());
        assert!(templates.join("grid.typ").is_file());
        assert!(templates.join("layout.typ").is_file());
        assert!(templates.join("layout-core.typ").is_file());
        assert!(templates.join("physics-edge-style.typ").is_file());
        assert!(templates.join("impl/physics-edge-style.typ").is_file());

        assert!(templates
            .join("crates/linnest/typst/linnest.wasm")
            .is_file());
        assert!(templates
            .join("crates/linnest/typst/src/curve.typ")
            .is_file());
        assert!(templates
            .join("crates/linnest/typst/src/impl/graph.typ")
            .is_file());
        assert!(templates
            .join("crates/linnest/typst/src/render/layout.typ")
            .is_file());
        assert!(!templates
            .join("crates/linnest/typst/src/physics-edge-style.typ")
            .exists());
        assert!(!templates
            .join("crates/linnest/typst/src/impl/physics-edge-style.typ")
            .exists());
        assert!(templates.join("crates/kurvst/typst/kurvst.wasm").is_file());
        assert!(templates.join("crates/kurvst/typst/src/lib.typ").is_file());
        assert!(templates.join("crates/kurvst/typst/src/impl.typ").is_file());

        assert!(!templates.join("linnest.wasm").exists());
        assert!(!templates.join("kurvst.wasm").exists());
        assert!(!templates.join("curve.typ").exists());
        for package in [LINNEST_PACKAGE_DIR, KURVST_PACKAGE_DIR] {
            for excluded in ["docs", "examples", "benches", "README.md", "LICENSE"] {
                assert!(!templates.join(package).join(excluded).exists());
            }
        }

        Assets::extract_justfile(tempdir.path())?;
        assert!(tempdir.path().join("justfile").is_file());

        Ok(())
    }

    #[cfg(unix)]
    #[test]
    #[ignore = "requires Typst, just, and LINNET_TEST_EXECUTABLE pointing to a built CLI"]
    fn exported_model_drawings_render_and_cache_without_python() -> Result<()> {
        use std::collections::BTreeMap;
        use std::os::unix::fs::symlink;
        use std::path::PathBuf;
        use std::process::Command;

        let temp = tempfile::tempdir()?;
        let bin = temp.path().join("bin");
        fs::create_dir(&bin)?;
        let search_path = std::env::var_os("PATH").unwrap_or_default();
        for (name, configured) in [
            (
                "linnet",
                Some(
                    std::env::var_os("LINNET_TEST_EXECUTABLE").expect("set LINNET_TEST_EXECUTABLE"),
                ),
            ),
            ("typst", std::env::var_os("TYPST_TEST_EXECUTABLE")),
            ("just", None),
            ("sh", None),
        ] {
            let executable = configured.map(PathBuf::from).unwrap_or_else(|| {
                std::env::split_paths(&search_path)
                    .map(|directory| directory.join(name))
                    .find(|path| path.is_file())
                    .unwrap_or_else(|| panic!("{name} must be available on PATH"))
            });
            symlink(fs::canonicalize(executable)?, bin.join(name))?;
        }
        let absent = Command::new(bin.join("sh"))
            .env("PATH", &bin)
            .args(["-c", "! command -v python && ! command -v python3"])
            .status()?;
        assert!(absent.success());

        Assets::extract_templates(temp.path())?;
        Assets::extract_justfile(temp.path())?;
        let templates = temp.path().join("drawings/templates");
        let model = gammalooprs::utils::load_generic_model("scalars");
        model.generate_edge_style_template(templates.join("edge-style.typ"))?;
        let data = temp.path().join("processes/amplitudes/example.dot");
        fs::create_dir_all(data.parent().unwrap())?;
        fs::write(
            &data,
            r#"digraph amplitude {
            a [pos="0,0!"]; b [pos="4,0!"];
            incoming [style=invis]; outgoing [style=invis];
            incoming -> a [particle="scalar_0", pos="z:4"];
            a -> b [particle="scalar_0"];
            b -> outgoing [particle="scalar_0", pos="z:-2"];
        }"#,
        )?;

        for scenario in ["initial", "cached", "changed-source", "explicit-inputs"] {
            if scenario == "changed-source" {
                let implementation = templates.join("crates/linnest/typst/src/impl/draw.typ");
                let source = fs::read_to_string(&implementation)?;
                fs::write(
                    implementation,
                    format!("{source}\n// Cache invalidation probe.\n"),
                )?;
            }
            let mut command = Command::new(bin.join("just"));
            command.current_dir(temp.path()).env("PATH", &bin).args([
                "steps=0",
                "columns=1",
                "draw",
            ]);
            if scenario == "explicit-inputs" {
                command.args(["--input", "steps=2", "--input", "columns=2"]);
            }
            let output = command.output()?;
            let expected = if scenario == "cached" {
                "0 built, 1 reused"
            } else {
                "1 built, 0 reused"
            };
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(
                String::from_utf8_lossy(&output.stdout).contains(expected),
                "{}",
                String::from_utf8_lossy(&output.stdout)
            );
            for path in [
                "drawings.pdf",
                "drawings/figs/processes/amplitudes/example.pdf",
            ] {
                assert!(fs::read(temp.path().join(path))?.starts_with(b"%PDF-"));
            }
            let metadata: serde_json::Value = serde_json::from_slice(&fs::read(
                temp.path().join("drawings/.cache/run-metadata.json"),
            )?)?;
            let inputs: Vec<(String, String)> = serde_json::from_value(metadata["input"].clone())?;
            let inputs: BTreeMap<_, _> = inputs.into_iter().collect();
            assert_eq!(
                inputs["steps"],
                if scenario == "explicit-inputs" {
                    "2"
                } else {
                    "0"
                }
            );
            assert_eq!(
                inputs["columns"],
                if scenario == "explicit-inputs" {
                    "2"
                } else {
                    "1"
                }
            );
        }

        let assertions = temp.path().join("model-style.typ");
        fs::write(
            &assertions,
            r#"#import "drawings/templates/edge-style.typ" as physics
#import "drawings/templates/crates/linnest/typst/src/lib.typ": graph
#let g = graph.parse(read("processes/amplitudes/example.dot")).first()
#let particle-label = physics.map.at("scalar_0").label
#let styles = physics.style(label-size: 12pt, label-fill: blue)
#let expected-label = text(size: 12pt, fill: blue)[#particle-label]
#assert(particle-label != none)
#assert(graph.edges(g).len() == 3)
// Drawing callbacks receive flattened fields, while graph.edges stores DOT
// attributes in statements. Use the public mapper's normalized view.
#let _ = graph.map(g, edge: edge => {
  assert((styles.edge-label)(edge.fields) == expected-label)
  none
})
[Model-generated particle labels work directly in Typst.]
"#,
        )?;
        let output = Command::new(bin.join("typst"))
            .env("PATH", &bin)
            .arg("compile")
            .arg("--root")
            .arg(temp.path())
            .arg(&assertions)
            .arg(temp.path().join("model-style.pdf"))
            .output()?;
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        Ok(())
    }
}
