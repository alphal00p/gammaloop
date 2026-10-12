//! Build-time asset selection shared by the Typst renderers.
//!
//! Include Dir has no include/exclude filter. Emit flat entries instead of
//! embedding package documentation, examples, and their generated PDFs.

use std::{fs, io, io::Write, path::Path};

pub fn embed_dir(
    output: &mut impl Write,
    name: &str,
    root: &Path,
    include: impl Fn(&Path) -> bool,
) -> io::Result<()> {
    fn collect(
        directory: &Path,
        root: &Path,
        include: &impl Fn(&Path) -> bool,
        files: &mut Vec<std::path::PathBuf>,
    ) -> io::Result<()> {
        println!("cargo:rerun-if-changed={}", directory.display());
        for entry in fs::read_dir(directory)? {
            let entry = entry?;
            let file_type = entry.file_type()?;
            if file_type.is_dir() {
                collect(&entry.path(), root, include, files)?;
            } else if file_type.is_file()
                && include(entry.path().strip_prefix(root).expect("asset under root"))
            {
                files.push(entry.path());
            }
        }
        Ok(())
    }

    let root = root.canonicalize()?;
    let mut files = Vec::new();
    collect(&root, &root, &include, &mut files)?;
    files.sort();
    writeln!(
        output,
        "static {name}: include_dir::Dir<'static> = include_dir::Dir::new(\"\", &["
    )?;
    for file in files {
        println!("cargo:rerun-if-changed={}", file.display());
        let relative = file
            .strip_prefix(&root)
            .expect("asset under root")
            .to_str()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "non-UTF-8 asset path"))?
            .replace('\\', "/");
        let path = file
            .to_str()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "non-UTF-8 asset path"))?;
        writeln!(
            output,
            "include_dir::DirEntry::File(include_dir::File::new({relative:?}, include_bytes!({path:?}))),"
        )?;
    }
    writeln!(output, "]);\n")
}

pub fn embed_typst_packages(
    output: &mut impl Write,
    crates_dir: &Path,
    additional_files: &[&str],
) -> io::Result<()> {
    for (name, package) in [
        ("LINNEST_PACKAGE_ASSETS", "linnest"),
        ("KURVST_PACKAGE_ASSETS", "kurvst"),
    ] {
        let wasm = format!("{package}.wasm");
        embed_dir(
            output,
            name,
            &crates_dir.join(package).join("typst"),
            |path| {
                (path.starts_with("src")
                    && path.extension().is_some_and(|extension| extension == "typ"))
                    || path == Path::new("typst.toml")
                    || path == Path::new(&wasm)
                    || additional_files.iter().any(|file| path == Path::new(file))
            },
        )?;
    }
    Ok(())
}
