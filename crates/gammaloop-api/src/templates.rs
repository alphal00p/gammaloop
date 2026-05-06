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
struct GammaLoopTemplateAssets;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../linnest/typst"]
#[include = "src/*.typ"]
#[include = "typst.toml"]
#[include = "linnest.wasm"]
struct LinnestPackageAssets;

#[derive(RustEmbed)]
#[folder = "$CARGO_MANIFEST_DIR/../kurvst/typst"]
#[include = "src/*.typ"]
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

    #[test]
    fn extract_templates_uses_canonical_package_layout() -> Result<()> {
        let tempdir = tempfile::tempdir()?;

        Assets::extract_templates(tempdir.path())?;

        let templates = tempdir.path().join("drawings/templates");
        assert!(templates.join("figure.typ").is_file());
        assert!(templates.join("grid.typ").is_file());
        assert!(templates.join("layout.typ").is_file());
        assert!(fs::read_to_string(templates.join("grid.typ"))?.contains("page_format"));

        assert!(templates
            .join("crates/linnest/typst/linnest.wasm")
            .is_file());
        assert!(templates
            .join("crates/linnest/typst/src/curve.typ")
            .is_file());
        assert!(templates.join("crates/kurvst/typst/kurvst.wasm").is_file());
        assert!(templates.join("crates/kurvst/typst/src/lib.typ").is_file());

        assert!(!templates.join("linnest.wasm").exists());
        assert!(!templates.join("kurvst.wasm").exists());
        assert!(!templates.join("curve.typ").exists());

        Ok(())
    }
}
