use color_eyre::Result;
use rust_embed::RustEmbed;
use std::fs;
use std::path::Path;

#[derive(RustEmbed)]
#[folder = "../assets/linnet-templates"]
pub struct Templates;

impl Templates {
    /// Extract all embedded templates to build/templates relative to target_dir
    pub fn extract_to_build_dir<P: AsRef<Path>>(target_dir: P) -> Result<()> {
        let target = target_dir.as_ref().join("build/templates");

        // Create the target directory if it doesn't exist
        fs::create_dir_all(&target)?;

        // Extract all embedded files
        for file_path in Self::iter() {
            if let Some(file) = Self::get(&file_path) {
                let target_path = target.join(file_path.as_ref());

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
}
