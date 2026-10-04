use color_eyre::Result;
use rust_embed::RustEmbed;
use std::fs;
use std::path::Path;

const LINNEST_PACKAGE_DIR: &str = "crates/linnest/typst";
const KURVST_PACKAGE_DIR: &str = "crates/kurvst/typst";

#[derive(RustEmbed)]
#[folder = "../../assets/embedded"]
pub struct Assets;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../../assets/embedded/drawing/templates"]
#[include = "*.typ"]
#[include = "impl/*.typ"]
struct GammaLoopTemplateAssets;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../linnest/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "linnest.wasm"]
struct LinnestPackageAssets;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../kurvst/typst"]
#[include = "src/*.typ"]
#[include = "src/**/*.typ"]
#[include = "typst.toml"]
#[include = "kurvst.wasm"]
struct KurvstPackageAssets;

impl Assets {
    /// Get all app template paths written under drawings/templates.
    pub fn template_paths() -> impl Iterator<Item = String> {
        GammaLoopTemplateAssets::iter().map(|path| path.to_string())
    }

    /// Extract all drawing templates to drawings/templates relative to target_dir.
    pub fn extract_templates<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let target = target_dir.as_ref().join("drawings/templates");

        fs::create_dir_all(&target)?;

        extract_package::<GammaLoopTemplateAssets>(&target)?;
        extract_package::<LinnestPackageAssets>(&target.join(LINNEST_PACKAGE_DIR))?;
        extract_package::<KurvstPackageAssets>(&target.join(KURVST_PACKAGE_DIR))?;

        Ok(())
    }

    /// Extract the Justfile to the target directory
    pub fn extract_justfile<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let justfile_path = target_dir.as_ref().join("justfile");

        if let Some(file) = Self::get("drawing/justfile") {
            fs::write(&justfile_path, file.data)?;
        }

        Ok(())
    }
}

fn extract_package<E: RustEmbed>(target: &Path) -> Result<()> {
    fs::create_dir_all(target)?;
    for package_path in E::iter() {
        if let Some(file) = E::get(package_path.as_ref()) {
            let target_path = target.join(package_path.as_ref());

            if let Some(parent) = target_path.parent() {
                fs::create_dir_all(parent)?;
            }

            fs::write(&target_path, file.data)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gammalooprs::model::ModelGammaLoopExt;

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
#let styles = physics.style()
#for edge in graph.edges(g) {
  assert(repr((styles.edge-label)(edge)).contains(repr(physics.map.at("scalar_0").label)))
}
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
