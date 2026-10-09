//! Offline Typst compilation with embedded fonts and an explicit package store.
//! The same compiler runs natively and inside Pyodide; it never imports Python
//! packages, invokes a subprocess, or downloads Typst packages.

use std::{
    collections::{BTreeMap, HashMap},
    fs,
    path::Path,
    sync::{LazyLock, Mutex},
};

use typst::LibraryExt;
use typst::syntax::{FileId, RootedPath, Source, VirtualPath, VirtualRoot};
use typst::utils::LazyHash;
use typst_kit::fonts::FontStore;
use typst_library::{
    Feature, Library, World, WorldExt,
    diag::{FileError, FileResult, SourceDiagnostic},
    foundations::{Bytes, Datetime, Duration, Value},
    text::{Font, FontBook},
};

static FONTS: LazyLock<FontStore> = LazyLock::new(|| {
    let mut fonts = FontStore::new();
    #[cfg(feature = "full")]
    fonts.extend(typst_kit::fonts::embedded());
    // Typst equations request weight 450 (Book), not the 400 Regular face.
    #[cfg(not(feature = "full"))]
    for data in [
        include_bytes!("../fonts/LibertinusSerif-Regular.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-Bold.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-Italic.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-BoldItalic.otf").as_slice(),
        include_bytes!("../fonts/NewCMMath-Book.otf").as_slice(),
        include_bytes!("../fonts/NewCMMath-Bold.otf").as_slice(),
        include_bytes!("../fonts/DejaVuSansMono.ttf").as_slice(),
    ] {
        fonts.extend(Font::iter(Bytes::new(data)).map(|font| {
            let info = font.info().clone();
            (font, info)
        }));
    }
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
    #[cfg(feature = "full")]
    fonts: Option<FontStore>,
}

impl<'a> Document<'a> {
    /// Resolve a named Typst color to its exact SVG hexadecimal paint.
    pub fn named_color_css(name: &str) -> Result<String, String> {
        static LIBRARY: LazyLock<Library> = LazyLock::new(|| Library::builder().build());
        match LIBRARY
            .global
            .scope()
            .get(name)
            .map(|binding| binding.read())
        {
            Some(Value::Color(color)) => Ok(color.to_hex().to_string()),
            _ => Err(format!(
                "unknown Typst color {name:?}; use an RGB or hexadecimal paint"
            )),
        }
    }

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
        #[cfg(feature = "full")]
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
            #[cfg(feature = "full")]
            fonts,
            library: LazyHash::new(
                Library::builder()
                    .with_features([Feature::Html].into_iter().collect())
                    .build(),
            ),
        }
    }

    /// Compile HTML as one document or SVG as one document per page.
    /// The `full` feature additionally enables PDF and PNG output.
    pub fn compile(&self, format: &str) -> Result<Vec<Vec<u8>>, String> {
        let supported = matches!(format, "html" | "svg")
            || cfg!(feature = "full") && matches!(format, "png" | "pdf");
        if !supported {
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
            #[cfg(feature = "full")]
            "png" => document
                .pages()
                .iter()
                .map(|page| {
                    typst_render::render(page, &typst_render::RenderOptions::default())
                        .encode_png()
                        .map_err(|error| error.to_string())
                })
                .collect(),
            #[cfg(feature = "full")]
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
        #[cfg(feature = "full")]
        if let Some(fonts) = &self.fonts {
            return fonts.book();
        }
        FONTS.book()
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
        #[cfg(feature = "full")]
        if let Some(fonts) = &self.fonts {
            return fonts.font(index);
        }
        FONTS.font(index)
    }
    fn today(&self, _offset: Option<Duration>) -> Option<Datetime> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn named_color_svg_paints_use_typst_constants_not_css_names() {
        assert_eq!(Document::named_color_css("red").unwrap(), "#ff4136");
        assert_eq!(Document::named_color_css("blue").unwrap(), "#0074d9");
        assert_eq!(Document::named_color_css("eastern").unwrap(), "#239dad");
        assert!(Document::named_color_css("unknown-mark-color").is_err());
        assert!(Document::named_color_css("rect").is_err());
    }

    #[test]
    fn embedded_fonts_match_default_math_and_text_variants() {
        use typst_library::text::{FontStretch, FontStyle, FontVariant, FontWeight};
        for (family, style, weight) in [
            ("new computer modern math", FontStyle::Normal, 450),
            ("new computer modern math", FontStyle::Normal, 700),
            ("libertinus serif", FontStyle::Normal, 400),
            ("libertinus serif", FontStyle::Normal, 700),
            ("libertinus serif", FontStyle::Italic, 400),
            ("libertinus serif", FontStyle::Italic, 700),
            ("dejavu sans mono", FontStyle::Normal, 400),
        ] {
            let variant =
                FontVariant::new(style, FontWeight::from_number(weight), FontStretch::NORMAL);
            let index = FONTS.book().select(family, variant).unwrap();
            let font = FONTS.font(index).unwrap();
            assert_eq!(font.info().variant, variant, "{family}");
        }
    }

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
    #[test]
    fn renders_math_styles_and_raw_text() {
        let files = BTreeMap::from([(
            "main.typ".into(),
            br#"#set page(width: auto, height: auto, margin: 0pt)
Regular *bold* _italic_ *_bold italic_* `raw <&>`
$bold(x) + italic(y) + cal(L) + bb(R) + sum_(i=1)^n frac(alpha_i, sqrt(beta))$"#
                .to_vec(),
        )]);
        for format in ["svg", "html"] {
            let pages = Document::compile_sources(&files, format).unwrap();
            assert_eq!(pages.len(), 1);
            let rendered = String::from_utf8(pages[0].clone()).unwrap();
            assert!(rendered.contains(if format == "svg" { "<path" } else { "<math" }));
            assert!(!rendered.contains("@preview/"));
        }
    }

    #[cfg(not(feature = "full"))]
    #[test]
    fn lean_compiler_excludes_document_only_features() {
        for (source, expected) in [
            ("#plugin(bytes(()))", "unknown variable: plugin"),
            (
                r#"#bibliography("refs.bib")"#,
                "unknown variable: bibliography",
            ),
            (
                r#"#image(bytes("%PDF-1.4"), format: "pdf")"#,
                "PDF images are not enabled",
            ),
            (
                r#"#raw("x", theme: bytes("theme"))"#,
                "custom syntax highlighting is not enabled",
            ),
        ] {
            let files = BTreeMap::from([("main.typ".into(), source.as_bytes().to_vec())]);
            let error = Document::compile_sources(&files, "svg").unwrap_err();
            assert!(error.contains(expected), "{error}");
        }
        let files = BTreeMap::from([("main.typ".into(), b"$x$".to_vec())]);
        for format in ["pdf", "png"] {
            assert!(
                Document::compile_sources(&files, format)
                    .unwrap_err()
                    .contains("unsupported Typst output format")
            );
        }
    }

    #[cfg(feature = "full")]
    #[test]
    fn full_compiler_keeps_exporters_and_document_features() {
        let files = BTreeMap::from([
            (
                "main.typ".into(),
                br#"#set page(width: 200pt, height: auto)
#assert(type(plugin) == function)
#raw("let x = 1;", lang: "rust")
A citation @math.
#bibliography("refs.bib")"#
                    .to_vec(),
            ),
            (
                "refs.bib".into(),
                br#"@book{math, title={Mathematics}, author={Euler, Leonhard}, year={1748}}"#
                    .to_vec(),
            ),
        ]);
        let pdf = Document::compile_sources(&files, "pdf").unwrap().remove(0);
        assert!(pdf.starts_with(b"%PDF-"));
        let png = Document::compile_sources(&files, "png").unwrap().remove(0);
        assert!(png.starts_with(b"\x89PNG\r\n\x1a\n"));
        // Check the PDF importer and SVG converter as well as the exporters.
        let files = BTreeMap::from([
            (
                "main.typ".into(),
                b"#set page(width: auto, height: auto)\n#image(\"page.pdf\")".to_vec(),
            ),
            ("page.pdf".into(), pdf),
        ]);
        assert!(
            String::from_utf8(Document::compile_sources(&files, "svg").unwrap().remove(0))
                .unwrap()
                .contains("<svg")
        );
    }
}
