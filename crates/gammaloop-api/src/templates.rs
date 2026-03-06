use color_eyre::Result;
use rust_embed::RustEmbed;
use std::fs;
use std::path::Path;

#[derive(RustEmbed)]
#[folder = "../../assets/embedded"]
pub struct Assets;

impl Assets {
    /// Get all template file paths from drawing/templates
    pub fn template_paths() -> impl Iterator<Item = String> {
        Self::iter()
            .filter(|path| path.starts_with("drawing/templates/"))
            .map(|path| path.to_string())
    }

    /// Extract all template files from drawing/templates to drawings/templates relative to target_dir
    pub fn extract_templates<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let target = target_dir.as_ref().join("drawings/templates");

        // Create the target directory if it doesn't exist
        fs::create_dir_all(&target)?;

        // Extract template files from drawing/templates directory
        for template_path in Self::template_paths() {
            if let Some(file) = Self::get(&template_path) {
                let relative_path = template_path.strip_prefix("drawing/templates/").unwrap();
                let target_path = target.join(relative_path);

                // Create parent directories if they don't exist
                if let Some(parent) = target_path.parent() {
                    fs::create_dir_all(parent)?;
                }

                // Write the file content
                fs::write(&target_path, file.data)?;
            }
        }

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
