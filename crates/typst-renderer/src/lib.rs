//! Offline Typst compilation with embedded fonts and an explicit package store.
//! The same compiler runs natively and inside Pyodide; it never imports Python
//! packages, invokes a subprocess, or downloads Typst packages.

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    path::Path,
    sync::{LazyLock, Mutex},
};

use typst::{
    Feature, Library, LibraryExt, World, WorldExt,
    diag::{FileError, FileResult, SourceDiagnostic},
    foundations::{Bytes, Datetime, Duration},
    syntax::{FileId, RootedPath, Source, VirtualPath, VirtualRoot},
    text::{Font, FontBook},
    utils::LazyHash,
};
use typst_kit::fonts::FontStore;

static FONTS: LazyLock<FontStore> = LazyLock::new(|| {
    let mut fonts = FontStore::new();
    fonts.extend(typst_kit::fonts::embedded());
    fonts
});

/// One immutable render project. Authored sources override files staged on disk;
/// package imports resolve exclusively within the supplied package directory.
pub struct Document<'a> {
    root: &'a Path,
    packages: &'a Path,
    files: &'a BTreeMap<String, Vec<u8>>,
    sources: Mutex<HashMap<FileId, Source>>,
    library: LazyHash<Library>,
    fonts: Option<FontStore>,
}

impl<'a> Document<'a> {
    /// Compile a virtual document without any package store or plugin assets.
    pub fn compile_sources(
        files: &BTreeMap<String, Vec<u8>>,
        format: &str,
    ) -> Result<Vec<Vec<u8>>, String> {
        let root = tempfile::tempdir().map_err(|e| e.to_string())?;
        Document::new(root.path(), root.path(), files).compile(format)
    }

    /// A self-contained Typst export of an already drawn SVG figure.
    pub fn svg_source(svg: &str) -> String {
        format!("#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes({svg:?}))")
    }

    pub fn new(root: &'a Path, packages: &'a Path, files: &'a BTreeMap<String, Vec<u8>>) -> Self {
        let fonts = std::env::var_os("TYPST_FONT_PATHS").map(|paths| {
            let mut fonts = FontStore::new();
            for path in std::env::split_paths(&paths) {
                fonts.extend(typst_kit::fonts::scan(&path));
            }
            fonts.extend(typst_kit::fonts::embedded());
            fonts
        });
        Self {
            root,
            packages,
            files,
            sources: Mutex::new(HashMap::new()),
            fonts,
            library: LazyHash::new(
                Library::builder()
                    .with_features([Feature::Html].into_iter().collect())
                    .build(),
            ),
        }
    }

    /// Compile HTML/PDF as one document, or SVG/PNG as one document per page.
    pub fn compile(&self, format: &str) -> Result<Vec<Vec<u8>>, String> {
        if !matches!(format, "html" | "svg" | "png" | "pdf") {
            return Err(format!("unsupported Typst output format {format:?}"));
        }
        if format == "html" {
            let document = typst::compile::<typst_html::HtmlDocument>(self)
                .output
                .map_err(|errors| self.diagnostics(&errors))?;
            let html = typst_html::html(&document, &typst_html::HtmlOptions { pretty: true })
                .map_err(|errors| self.diagnostics(&errors))?;
            return Ok(vec![html.into_bytes()]);
        }
        let document = typst::compile::<typst_layout::PagedDocument>(self)
            .output
            .map_err(|errors| self.diagnostics(&errors))?;
        match format {
            "svg" => Ok(document
                .pages()
                .iter()
                .map(|page| typst_svg::svg(page, &typst_svg::SvgOptions::default()).into_bytes())
                .collect()),
            "png" => document
                .pages()
                .iter()
                .map(|page| {
                    typst_render::render(page, &typst_render::RenderOptions::default())
                        .encode_png()
                        .map_err(|error| error.to_string())
                })
                .collect(),
            "pdf" => typst_pdf::pdf(&document, &typst_pdf::PdfOptions::default())
                .map(|pdf| vec![pdf])
                .map_err(|errors| self.diagnostics(&errors)),
            _ => Err(format!("unsupported Typst output format {format:?}")),
        }
    }

    fn diagnostics(&self, errors: &[SourceDiagnostic]) -> String {
        errors
            .iter()
            .map(|error| {
                let location = error.span.id().map_or_else(String::new, |id| {
                    let line = self.source(id).ok().and_then(|source| {
                        let range = self.range(error.span)?;
                        source
                            .lines()
                            .byte_to_line(range.start)
                            .map(|line| line + 1)
                    });
                    format!(
                        "{:?}{}: ",
                        id.vpath(),
                        line.map_or_else(String::new, |line| format!(":{line}"))
                    )
                });
                format!("{location}{}", error.message)
            })
            .collect::<Vec<_>>()
            .join("\n")
    }
}

impl Drop for Document<'_> {
    fn drop(&mut self) {
        // Keep reusable layout/plugin work, but bound the lifetime of old
        // notebook render inputs, including failed compilations.
        comemo::evict(5);
    }
}

impl World for Document<'_> {
    fn library(&self) -> &LazyHash<Library> {
        &self.library
    }
    fn book(&self) -> &LazyHash<FontBook> {
        self.fonts.as_ref().unwrap_or(&FONTS).book()
    }
    fn main(&self) -> FileId {
        RootedPath::new(VirtualRoot::Project, VirtualPath::new("main.typ").unwrap()).intern()
    }
    fn source(&self, id: FileId) -> FileResult<Source> {
        let mut sources = self.sources.lock().unwrap();
        if let Some(source) = sources.get(&id) {
            return Ok(source.clone());
        }
        let bytes = self.file(id)?;
        let text = std::str::from_utf8(&bytes)?;
        let source = Source::new(id, text.into());
        sources.insert(id, source.clone());
        Ok(source)
    }
    fn file(&self, id: FileId) -> FileResult<Bytes> {
        let root = match id.root() {
            VirtualRoot::Project => {
                let key = id.vpath().get_without_slash();
                if let Some(bytes) = self.files.get(key) {
                    return Ok(Bytes::new(bytes.clone()));
                }
                self.root.to_path_buf()
            }
            VirtualRoot::Package(package) => self
                .packages
                .join(package.namespace.as_str())
                .join(package.name.as_str())
                .join(package.version.to_string()),
        };
        let path = id.vpath().realize(&root)?;
        fs::read(&path)
            .map(Bytes::new)
            .map_err(|error| FileError::from_io(error, &path))
    }
    fn font(&self, index: usize) -> Option<Font> {
        self.fonts.as_ref().unwrap_or(&FONTS).font(index)
    }
    fn today(&self, _offset: Option<Duration>) -> Option<Datetime> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_svg_export_compiles_without_packages() {
        let svg = r##"<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"><circle cx="10" cy="10" r="5" fill="#123456"/></svg>"##;
        let files = BTreeMap::from([("main.typ".into(), Document::svg_source(svg).into_bytes())]);
        let result = Document::compile_sources(&files, "svg").unwrap();
        assert!(
            String::from_utf8(result[0].clone())
                .unwrap()
                .contains("<svg")
        );
    }

    #[test]
    fn renders_svg_and_mathml_without_system_fonts_or_packages() {
        let root = tempfile::tempdir().unwrap();
        let files = BTreeMap::from([(
            "main.typ".to_owned(),
            b"#set page(width: auto, height: auto)\n$x^2 + 1$".to_vec(),
        )]);
        let document = Document::new(root.path(), root.path(), &files);
        assert!(
            String::from_utf8(document.compile("svg").unwrap().remove(0))
                .unwrap()
                .contains("<svg")
        );
        assert!(
            String::from_utf8(document.compile("html").unwrap().remove(0))
                .unwrap()
                .contains("<math")
        );
    }
}
