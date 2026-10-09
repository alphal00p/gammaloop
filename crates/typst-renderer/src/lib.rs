//! Offline Typst compilation with embedded fonts and explicit package directories.
//! The same compiler runs natively and inside Pyodide; it never imports Python
//! packages, invokes a subprocess, or downloads Typst packages.

use std::{
    collections::{BTreeMap, HashMap},
    env, fs,
    path::{Path, PathBuf},
    sync::{LazyLock, Mutex},
};

use typst::LibraryExt;
use typst::syntax::{FileId, RootedPath, Source, Span, VirtualPath, VirtualRoot};
use typst::utils::LazyHash;
use typst_kit::fonts::FontStore;
use typst_library::{
    Feature, Library, World, WorldExt,
    diag::{FileError, FileResult, PackageError, SourceDiagnostic},
    foundations::{
        Bytes, Datetime, Duration, NativeElement, Recipe, Selector, Styles, Transformation, Value,
    },
    math::EquationElem,
    text::{Font, FontBook, FontFamily, FontList, TextElem},
};

/// The only fonts: Libertinus Serif for text, STIX Two Math for equations and
/// DejaVu Sans Mono for raw text. All are under permissive licenses.
static FONTS: LazyLock<FontStore> = LazyLock::new(|| {
    let mut fonts = FontStore::new();
    for data in [
        include_bytes!("../fonts/LibertinusSerif-Regular.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-Semibold.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-Bold.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-Italic.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-SemiboldItalic.otf").as_slice(),
        include_bytes!("../fonts/LibertinusSerif-BoldItalic.otf").as_slice(),
        include_bytes!("../fonts/STIXTwoMath-Regular.otf").as_slice(),
        include_bytes!("../fonts/DejaVuSansMono.ttf").as_slice(),
        include_bytes!("../fonts/DejaVuSansMono-Bold.ttf").as_slice(),
        include_bytes!("../fonts/DejaVuSansMono-Oblique.ttf").as_slice(),
        include_bytes!("../fonts/DejaVuSansMono-BoldOblique.ttf").as_slice(),
    ] {
        fonts.extend(Font::iter(Bytes::new(data)).map(|font| {
            let info = font.info().clone();
            (font, info)
        }));
    }
    fonts
});

/// One immutable render project. Authored sources override files staged on disk,
/// if the project has a directory; package imports resolve in the first supplied
/// package directory that holds them, each laid out as `namespace/name/version`.
pub struct Document<'a> {
    root: Option<&'a Path>,
    packages: &'a [PathBuf],
    files: &'a BTreeMap<String, Vec<u8>>,
    sources: Mutex<HashMap<FileId, Source>>,
    library: LazyHash<Library>,
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

    /// Compile a virtual document from its sources alone, without a project
    /// directory, package store or plugin assets; nothing touches the disk.
    pub fn compile_sources(
        files: &BTreeMap<String, Vec<u8>>,
        format: &str,
    ) -> Result<Vec<Vec<u8>>, String> {
        Document::with_root(None, &[], files).compile(format)
    }

    /// The package directories of the Typst CLI, in its lookup order: the
    /// package path (`TYPST_PACKAGE_PATH`, else `<data dir>/typst/packages`),
    /// then the package cache (`TYPST_PACKAGE_CACHE_PATH`, else
    /// `<cache dir>/typst/packages`). Unlike the CLI, nothing is downloaded.
    pub fn local_package_dirs() -> Vec<PathBuf> {
        [
            ("TYPST_PACKAGE_PATH", dirs::data_dir()),
            ("TYPST_PACKAGE_CACHE_PATH", dirs::cache_dir()),
        ]
        .into_iter()
        .filter_map(|(variable, default)| {
            env::var_os(variable)
                .map(PathBuf::from)
                .or_else(|| default.map(|directory| directory.join("typst/packages")))
        })
        .collect()
    }

    /// A self-contained Typst export of an already drawn SVG figure.
    pub fn svg_source(svg: &str) -> String {
        format!("#set page(width: auto, height: auto, margin: 0pt)\n#image(bytes({svg:?}))")
    }

    pub fn new(
        root: &'a Path,
        packages: &'a [PathBuf],
        files: &'a BTreeMap<String, Vec<u8>>,
    ) -> Self {
        Self::with_root(Some(root), packages, files)
    }

    fn with_root(
        root: Option<&'a Path>,
        packages: &'a [PathBuf],
        files: &'a BTreeMap<String, Vec<u8>>,
    ) -> Self {
        Self {
            root,
            packages,
            files,
            sources: Mutex::new(HashMap::new()),
            library: LazyHash::new(Self::library()),
        }
    }

    /// The standard library with STIX Two Math as the equation font. Typst's
    /// default, New Computer Modern Math, is not embedded; a font set lacking
    /// it would make equations fall back to a text font and lose their layout.
    fn library() -> Library {
        let mut library = Library::builder()
            .with_features([Feature::Html].into_iter().collect())
            .build();
        let mut math = Styles::new();
        math.set(
            TextElem::font,
            FontList(vec![FontFamily::new("STIX Two Math")]),
        );
        library.styles.push(Recipe::new(
            Some(Selector::Elem(EquationElem::ELEM, None)),
            Transformation::Style(math),
            Span::detached(),
        ));
        library
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
                let Some(root) = self.root else {
                    return Err(FileError::NotFound(key.into()));
                };
                root.to_path_buf()
            }
            VirtualRoot::Package(package) => {
                let relative = Path::new(package.namespace.as_str())
                    .join(package.name.as_str())
                    .join(package.version.to_string());
                self.packages
                    .iter()
                    .map(|directory| directory.join(&relative))
                    .find(|directory| directory.is_dir())
                    .ok_or_else(|| {
                        let searched = self
                            .packages
                            .iter()
                            .map(|directory| directory.display().to_string())
                            .collect::<Vec<_>>()
                            .join(", ");
                        FileError::Package(PackageError::Other(Some(
                            format!(
                                "{package} is not installed (searched: {searched}); \
                                 compile once with the Typst CLI to cache it, or set \
                                 TYPST_PACKAGE_PATH to a directory holding {}",
                                relative.display()
                            )
                            .into(),
                        )))
                    })?
            }
        };
        let path = id.vpath().realize(&root)?;
        fs::read(&path)
            .map(Bytes::new)
            .map_err(|error| FileError::from_io(error, &path))
    }
    fn font(&self, index: usize) -> Option<Font> {
        FONTS.font(index)
    }
    fn today(&self, _offset: Option<Duration>) -> Option<Datetime> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;
    use typst_layout::PagedDocument;
    use typst_library::layout::{Frame, FrameItem};

    /// The families of all glyph runs in a compiled document.
    fn font_families(source: &str) -> BTreeSet<String> {
        fn collect(frame: &Frame, families: &mut BTreeSet<String>) {
            for (_, item) in frame.items() {
                match item {
                    FrameItem::Group(group) => collect(&group.frame, families),
                    FrameItem::Text(text) => {
                        families.insert(text.font.info().family.clone());
                    }
                    _ => {}
                }
            }
        }
        let root = tempfile::tempdir().unwrap();
        let files = BTreeMap::from([("main.typ".into(), source.as_bytes().to_vec())]);
        let document = Document::new(root.path(), &[], &files);
        let paged = typst::compile::<PagedDocument>(&document).output.unwrap();
        let mut families = BTreeSet::new();
        for page in paged.pages() {
            collect(&page.frame, &mut families);
        }
        families
    }

    #[test]
    fn virtual_documents_read_no_project_directory() {
        let files = BTreeMap::from([("main.typ".to_owned(), b"#read(\"Cargo.toml\")".to_vec())]);
        let error = Document::compile_sources(&files, "svg").unwrap_err();
        assert!(error.contains("not found"), "{error}");
    }

    #[test]
    fn packages_resolve_in_the_first_directory_that_holds_them() {
        let stores = tempfile::tempdir().unwrap();
        let directories: Vec<PathBuf> = ["first", "second"]
            .into_iter()
            .map(|name| stores.path().join(name))
            .collect();
        for (directory, packages) in directories
            .iter()
            .zip([&["shared"][..], &["shared", "later"]])
        {
            for package in packages {
                let root = directory.join("preview").join(package).join("1.0.0");
                fs::create_dir_all(&root).unwrap();
                fs::write(
                    root.join("typst.toml"),
                    format!("[package]\nname = \"{package}\"\nversion = \"1.0.0\"\nentrypoint = \"lib.typ\"\n"),
                )
                .unwrap();
                let marker = directory.file_name().unwrap().to_str().unwrap();
                fs::write(
                    root.join("lib.typ"),
                    format!("#let marker = \"{marker}\"\n"),
                )
                .unwrap();
            }
        }
        let files = BTreeMap::from([(
            "main.typ".to_owned(),
            b"#import \"@preview/shared:1.0.0\": marker as shared\n\
              #import \"@preview/later:1.0.0\": marker as later\n\
              #assert.eq(shared, \"first\")\n#assert.eq(later, \"second\")\n[ok]"
                .to_vec(),
        )]);
        let root = tempfile::tempdir().unwrap();
        assert!(
            Document::new(root.path(), &directories, &files)
                .compile("svg")
                .is_ok()
        );

        // A package in no directory names itself, the directories and the remedy.
        let error = Document::new(root.path(), &directories[..1], &files)
            .compile("svg")
            .unwrap_err();
        for expected in [
            "@preview/later:1.0.0 is not installed",
            directories[0].to_str().unwrap(),
            "TYPST_PACKAGE_PATH",
        ] {
            assert!(error.contains(expected), "{error}");
        }
    }

    #[test]
    fn named_color_svg_paints_use_typst_constants_not_css_names() {
        assert_eq!(Document::named_color_css("red").unwrap(), "#ff4136");
        assert_eq!(Document::named_color_css("blue").unwrap(), "#0074d9");
        assert_eq!(Document::named_color_css("eastern").unwrap(), "#239dad");
        assert!(Document::named_color_css("unknown-mark-color").is_err());
        assert!(Document::named_color_css("rect").is_err());
    }

    #[test]
    fn embedded_fonts_cover_text_math_and_raw_variants() {
        use typst_library::text::{FontStretch, FontStyle, FontVariant, FontWeight};
        for (family, style, weight) in [
            ("stix two math", FontStyle::Normal, 400),
            ("libertinus serif", FontStyle::Normal, 400),
            ("libertinus serif", FontStyle::Normal, 600),
            ("libertinus serif", FontStyle::Normal, 700),
            ("libertinus serif", FontStyle::Italic, 400),
            ("libertinus serif", FontStyle::Italic, 600),
            ("libertinus serif", FontStyle::Italic, 700),
            ("dejavu sans mono", FontStyle::Normal, 400),
            ("dejavu sans mono", FontStyle::Normal, 700),
            ("dejavu sans mono", FontStyle::Italic, 400),
            ("dejavu sans mono", FontStyle::Italic, 700),
        ] {
            let variant =
                FontVariant::new(style, FontWeight::from_number(weight), FontStretch::NORMAL);
            let index = FONTS.book().select(family, variant).unwrap();
            let font = FONTS.font(index).unwrap();
            assert_eq!(font.info().variant, variant, "{family}");
        }
    }

    #[test]
    fn typesets_text_equations_and_raw_text_in_their_own_fonts() {
        for (source, family) in [
            ("Text *bold* _italic_", "Libertinus Serif"),
            ("$x + sum_(i=1)^n cal(L)_i$", "STIX Two Math"),
            ("$ bold(x) = integral_0^1 f(t) dif t $", "STIX Two Math"),
            ("`raw`", "DejaVu Sans Mono"),
        ] {
            assert_eq!(
                font_families(source),
                BTreeSet::from([family.to_owned()]),
                "{source}"
            );
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
        let document = Document::new(root.path(), &[], &files);
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
