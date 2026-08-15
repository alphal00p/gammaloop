//! Build and validate the five independently rendered product manuals.

mod server;
mod watch;

pub use watch::WatchRequest;

use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    path::{Component, Path, PathBuf},
    process::Command,
};

use alphal00p_docs_schema::{
    ApiLanguage, DocCatalog, DocFormat, DocItem, DocMember, DocScope, SCHEMA_VERSION,
    generated::{GENERATED_REFERENCE_SCHEMA, GammaLoopReference, VakintReference},
};
use clap::ValueEnum;
use eyre::{Context, ContextCompat, Result, bail, ensure};
use pulldown_cmark::{Event, Options, Parser, html};
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use tempfile::Builder as TempDirBuilder;
use walkdir::WalkDir;

const PRODUCT_IDS: [&str; 5] = ["gammaloop", "linnet", "spenso", "idenso", "vakint"];
const PORTAL_SCHEMA_VERSION: u32 = 2;
const DEVELOPER_SCHEMA_VERSION: u32 = 1;
const PUBLICATION_SCHEMA_VERSION: u32 = 1;
const STRICT_RUSTDOC_FLAGS: &str = "-D rustdoc::broken_intra_doc_links \
    -D rustdoc::invalid_html_tags -D rustdoc::bare_urls";

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum BuildChannel {
    Latest,
    Snapshot,
}

#[derive(Clone, Debug)]
pub struct BuildRequest {
    pub product: String,
    pub channel: BuildChannel,
    pub output: PathBuf,
    pub snapshot_tag: Option<String>,
    pub include_rustdoc: bool,
    pub include_typst: bool,
    /// Persistent Cargo target used only by a long-running watch session.
    pub rustdoc_target_root: Option<PathBuf>,
    /// Directory receiving one NUL-delimited Typst dependency file per product.
    pub dependency_output: Option<PathBuf>,
}

struct ProductBuildOptions<'a> {
    channel: BuildChannel,
    tag: Option<&'a str>,
    output: &'a Path,
    include_typst: bool,
    include_rustdoc: bool,
    rustdoc_target_root: Option<&'a Path>,
    dependency_output: Option<&'a Path>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ProductRegistry {
    schema: u32,
    product: Vec<ProductConfig>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ProductConfig {
    id: String,
    title: String,
    tagline: String,
    logo: Option<String>,
    source: PathBuf,
    citation: CitationConfig,
    #[serde(default)]
    pages: Vec<PageConfig>,
    changelog: Option<PathBuf>,
    #[serde(default)]
    related: Vec<String>,
    #[serde(default)]
    rust_components: Vec<ComponentConfig>,
    #[serde(default)]
    python_components: Vec<ComponentConfig>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct CitationConfig {
    title: String,
    creators: Vec<String>,
    year: u16,
    repository: String,
    doi: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PortalConfig {
    schema: u32,
    eyebrow: String,
    title: String,
    summary: String,
    funding: String,
    funding_url: String,
    #[serde(default)]
    pillar: Vec<PortalPillar>,
    #[serde(default)]
    people: Vec<PortalPerson>,
    #[serde(default)]
    affiliation: Vec<PortalAffiliation>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PortalPillar {
    label: String,
    title: String,
    summary: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PortalPerson {
    id: String,
    name: String,
    initials: String,
    role: String,
    url: String,
    github: String,
    portrait: Option<String>,
    portrait_source: Option<String>,
    inspire_recid: Option<u64>,
    inspire_bai: Option<String>,
    orcid: Option<String>,
    #[serde(default)]
    publications: bool,
    #[serde(default)]
    featured: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PortalAffiliation {
    name: String,
    location: String,
    summary: String,
    url: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperConfig {
    schema: u32,
    title: String,
    summary: String,
    #[serde(default)]
    section: Vec<DeveloperSection>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperSection {
    id: String,
    title: String,
    summary: String,
    #[serde(default)]
    note: Vec<DeveloperNote>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperNote {
    id: String,
    title: String,
    summary: String,
    source: PathBuf,
    kind: String,
    status: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PublicationCache {
    schema: u32,
    source: String,
    updated: String,
    api_url: String,
    authors: Vec<PublicationAuthor>,
    publications: Vec<Publication>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PublicationAuthor {
    id: String,
    name: String,
    inspire_bai: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct Publication {
    id: u64,
    title: String,
    date: String,
    year: u16,
    authors: Vec<String>,
    people: Vec<String>,
    venue: Option<String>,
    doi: Option<String>,
    arxiv: Option<String>,
    citations: u64,
    types: Vec<String>,
    url: String,
    bibtex_url: String,
}

/// One authored chapter in a product's ordered, durable book tree.
///
/// The Typst file owns prose while the registry owns navigation and URLs. This
/// mirrors the separation used by book-oriented Typst sites without coupling
/// the renderer to a particular theme package.
#[derive(Clone, Debug, Deserialize, Serialize)]
struct PageConfig {
    id: String,
    title: String,
    summary: String,
    source: PathBuf,
    symbol: String,
    route: String,
    group: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ComponentConfig {
    id: String,
    package: String,
    kind: String,
    version_source: PathBuf,
    module: Option<String>,
    #[serde(default)]
    features: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
struct ComponentVersion<'a> {
    product: &'a str,
    id: &'a str,
    package: &'a str,
    kind: &'a str,
    module: Option<&'a str>,
    version: String,
    features: &'a [String],
}

#[derive(Clone, Debug, Serialize)]
struct SnapshotMetadata<'a> {
    schema: u32,
    product: &'a str,
    title: &'a str,
    channel: BuildChannel,
    snapshot_tag: Option<&'a str>,
    git_commit: String,
    git_timestamp: u64,
    route: String,
    components: Vec<ComponentVersion<'a>>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct SearchEntry {
    title: String,
    summary: String,
    href: String,
    kind: String,
}

#[derive(Clone, Debug)]
struct SitePage {
    route: String,
    title: String,
    group: String,
}

impl SitePage {
    fn new(route: impl Into<String>, title: impl Into<String>, group: impl Into<String>) -> Self {
        Self {
            route: route.into(),
            title: title.into(),
            group: group.into(),
        }
    }
}

fn supplemental_reference(product: &str) -> Option<(&'static str, &'static str)> {
    match product {
        "gammaloop" => Some(("reference/cli/", "CLI commands and settings")),
        "vakint" => Some(("reference/topologies/", "Topologies and dependencies")),
        _ => None,
    }
}

#[derive(Clone, Debug)]
struct HeadingLink {
    level: u8,
    title: String,
    id: String,
}

pub struct SiteBuilder {
    root: PathBuf,
    api_root: PathBuf,
    registry: ProductRegistry,
    portal: PortalConfig,
    developers: DeveloperConfig,
    publications: PublicationCache,
}

impl SiteBuilder {
    pub fn discover() -> Result<Self> {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .context("documentation builder is not inside a workspace")?
            .to_path_buf();
        Self::load(root)
    }

    fn load(root: PathBuf) -> Result<Self> {
        let registry_path = root.join("docs/products/registry.toml");
        let source = fs::read_to_string(&registry_path)
            .wrap_err_with(|| format!("failed to read {}", registry_path.display()))?;
        let registry = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", registry_path.display()))?;
        let portal_path = root.join("docs/portal.toml");
        let source = fs::read_to_string(&portal_path)
            .wrap_err_with(|| format!("failed to read {}", portal_path.display()))?;
        let portal = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", portal_path.display()))?;
        let developers_path = root.join("docs/developers.toml");
        let source = fs::read_to_string(&developers_path)
            .wrap_err_with(|| format!("failed to read {}", developers_path.display()))?;
        let developers = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", developers_path.display()))?;
        let publications_path = root.join("docs/data/publications.json");
        let publications = serde_json::from_slice(
            &fs::read(&publications_path)
                .wrap_err_with(|| format!("failed to read {}", publications_path.display()))?,
        )
        .wrap_err_with(|| format!("failed to parse {}", publications_path.display()))?;
        let api_root = env::var_os("ALPHAL00P_DOCS_API_ROOT")
            .map(PathBuf::from)
            .unwrap_or_else(|| root.join("docs/api"));
        Ok(Self {
            root,
            api_root,
            registry,
            portal,
            developers,
            publications,
        })
    }

    pub fn check(&self) -> Result<()> {
        ensure!(
            self.registry.schema == SCHEMA_VERSION,
            "registry schema {} does not match catalog schema {}",
            self.registry.schema,
            SCHEMA_VERSION
        );

        let ids = self
            .registry
            .product
            .iter()
            .map(|product| product.id.as_str())
            .collect::<BTreeSet<_>>();
        ensure!(
            ids == PRODUCT_IDS.into_iter().collect(),
            "registry must contain exactly: {}",
            PRODUCT_IDS.join(", ")
        );
        ensure!(
            self.portal.schema == PORTAL_SCHEMA_VERSION,
            "portal schema {} does not match renderer schema {}",
            self.portal.schema,
            PORTAL_SCHEMA_VERSION
        );
        ensure!(
            self.developers.schema == DEVELOPER_SCHEMA_VERSION,
            "developer schema {} does not match renderer schema {}",
            self.developers.schema,
            DEVELOPER_SCHEMA_VERSION
        );
        ensure!(
            !self.developers.title.trim().is_empty()
                && !self.developers.summary.trim().is_empty()
                && !self.developers.section.is_empty(),
            "developer documentation metadata must not be empty"
        );
        let mut developer_section_ids = BTreeSet::new();
        let mut developer_note_ids = BTreeSet::new();
        let mut developer_sources = BTreeSet::new();
        for section in &self.developers.section {
            validate_route_segment(&section.id)?;
            ensure!(
                developer_section_ids.insert(&section.id),
                "duplicate developer section {}",
                section.id
            );
            ensure!(
                !section.title.trim().is_empty()
                    && !section.summary.trim().is_empty()
                    && !section.note.is_empty(),
                "developer section {} is incomplete",
                section.id
            );
            for note in &section.note {
                validate_route_segment(&note.id)?;
                ensure!(
                    developer_note_ids.insert(&note.id),
                    "duplicate developer note {}",
                    note.id
                );
                ensure!(
                    developer_sources.insert(note.source.clone()),
                    "developer source {} is registered more than once",
                    note.source.display()
                );
                ensure!(
                    !note.title.trim().is_empty()
                        && !note.summary.trim().is_empty()
                        && !note.kind.trim().is_empty()
                        && !note.status.trim().is_empty(),
                    "developer note {} is incomplete",
                    note.id
                );
                ensure!(
                    matches!(
                        note.source
                            .extension()
                            .and_then(|extension| extension.to_str()),
                        Some("md" | "html")
                    ),
                    "developer note {} must be Markdown or HTML",
                    note.source.display()
                );
                self.require_file(&note.source)?;
            }
        }
        let architecture_root = self.root.join("docs/architecture");
        let mut architecture_sources = BTreeSet::new();
        for entry in fs::read_dir(&architecture_root)? {
            let path = entry?.path();
            if path.is_file()
                && matches!(
                    path.extension().and_then(|extension| extension.to_str()),
                    Some("md" | "html")
                )
                && path.file_name().and_then(|name| name.to_str()) != Some("architecture.md")
            {
                architecture_sources.insert(path.strip_prefix(&self.root)?.to_path_buf());
            }
        }
        ensure!(
            architecture_sources == developer_sources,
            "docs/architecture notes and docs/developers.toml differ; unclassified: {:?}; missing: {:?}",
            architecture_sources
                .difference(&developer_sources)
                .collect::<Vec<_>>(),
            developer_sources
                .difference(&architecture_sources)
                .collect::<Vec<_>>()
        );
        ensure!(
            !self.portal.title.trim().is_empty()
                && !self.portal.summary.trim().is_empty()
                && !self.portal.eyebrow.trim().is_empty(),
            "portal identity metadata must not be empty"
        );
        ensure!(
            !self.portal.pillar.is_empty()
                && !self.portal.people.is_empty()
                && !self.portal.affiliation.is_empty(),
            "portal must define research pillars, people, and affiliations"
        );
        ensure!(
            self.publications.schema == PUBLICATION_SCHEMA_VERSION,
            "publication schema {} does not match renderer schema {}",
            self.publications.schema,
            PUBLICATION_SCHEMA_VERSION
        );
        ensure!(
            self.publications
                .source
                .starts_with("https://inspirehep.net")
                && self
                    .publications
                    .api_url
                    .starts_with("https://inspirehep.net/api/literature"),
            "publication metadata must come from INSPIRE HEP"
        );
        let mut person_ids = BTreeSet::new();
        let mut scholarly_people = BTreeSet::new();
        for person in &self.portal.people {
            validate_route_segment(&person.id)?;
            ensure!(
                person_ids.insert(&person.id),
                "duplicate person {}",
                person.id
            );
            ensure!(
                !person.name.trim().is_empty()
                    && !person.initials.trim().is_empty()
                    && !person.role.trim().is_empty(),
                "person {} is incomplete",
                person.id
            );
            ensure!(
                person.url.starts_with("https://")
                    && person.github.starts_with("https://github.com/"),
                "person {} must have HTTPS profile and GitHub URLs",
                person.id
            );
            ensure!(
                person.portrait.is_some() == person.portrait_source.is_some(),
                "person {} must pair portrait and portrait_source",
                person.id
            );
            if let Some(portrait) = &person.portrait {
                validate_route_segment(portrait.trim_end_matches(".webp"))?;
                ensure!(
                    portrait.ends_with(".webp"),
                    "person portrait must be WebP: {portrait}"
                );
                self.require_file(&PathBuf::from("docs/assets/people").join(portrait))?;
            }
            if person.publications {
                ensure!(
                    person.inspire_recid.is_some() && person.inspire_bai.is_some(),
                    "scholarly person {} needs stable INSPIRE identifiers",
                    person.id
                );
                scholarly_people.insert(person.id.as_str());
            }
        }
        let cached_authors = self
            .publications
            .authors
            .iter()
            .map(|author| author.id.as_str())
            .collect::<BTreeSet<_>>();
        ensure!(
            cached_authors == scholarly_people,
            "publication cache authors differ from portal scholarly people"
        );
        let mut publication_ids = BTreeSet::new();
        for publication in &self.publications.publications {
            ensure!(
                publication_ids.insert(publication.id),
                "duplicate INSPIRE publication {}",
                publication.id
            );
            ensure!(
                !publication.title.trim().is_empty()
                    && !publication.authors.is_empty()
                    && !publication.people.is_empty()
                    && publication
                        .people
                        .iter()
                        .all(|person| scholarly_people.contains(person.as_str()))
                    && publication
                        .url
                        .starts_with("https://inspirehep.net/literature/")
                    && publication
                        .bibtex_url
                        .starts_with("https://inspirehep.net/api/literature/"),
                "publication {} is incomplete",
                publication.id
            );
        }
        for url in self
            .portal
            .people
            .iter()
            .map(|person| person.url.as_str())
            .chain(
                self.portal
                    .affiliation
                    .iter()
                    .map(|affiliation| affiliation.url.as_str()),
            )
            .chain(std::iter::once(self.portal.funding_url.as_str()))
        {
            ensure!(
                url.starts_with("https://"),
                "portal URL must use HTTPS: {url}"
            );
        }
        self.require_file(Path::new("docs/assets/site.css"))?;
        self.require_file(Path::new("docs/assets/site.js"))?;
        self.require_file(Path::new("docs/assets/local-unitarity-light.svg"))?;
        self.require_file(Path::new("docs/assets/local-unitarity-dark.svg"))?;
        self.require_file(Path::new("docs/assets/spensologo.svg"))?;
        self.require_file(Path::new("assets/gammalooplogo-light.svg"))?;
        self.require_file(Path::new("assets/gammalooplogo-dark.svg"))?;

        let mut component_ids = BTreeSet::new();
        for product in &self.registry.product {
            validate_route_segment(&product.id)?;
            ensure!(
                !product.title.trim().is_empty() && !product.tagline.trim().is_empty(),
                "{} has empty product identity metadata",
                product.id
            );
            ensure!(
                !product.citation.title.trim().is_empty()
                    && !product.citation.creators.is_empty()
                    && product.citation.repository.starts_with("https://"),
                "{} has incomplete citation metadata",
                product.id
            );
            if let Some(doi) = &product.citation.doi {
                ensure!(
                    doi.starts_with("10."),
                    "{} has invalid DOI {doi}",
                    product.id
                );
            }
            ensure!(
                product
                    .logo
                    .as_deref()
                    .is_none_or(|logo| matches!(logo, "gammaloop" | "spenso")),
                "{} has an unknown product logo",
                product.id
            );
            self.require_file(&product.source)?;
            ensure!(
                !product.pages.is_empty(),
                "{} has no authored documentation pages",
                product.id
            );
            let mut page_ids = BTreeSet::new();
            let mut page_routes = BTreeSet::new();
            let mut root_pages = 0;
            let mut tutorial_pages = 0;
            let mut manual_pages = 0;
            for page in &product.pages {
                validate_route_segment(&page.id)?;
                ensure!(
                    page_ids.insert(&page.id),
                    "{} has duplicate page id {}",
                    product.id,
                    page.id
                );
                ensure!(
                    page_routes.insert(&page.route),
                    "{} has duplicate page route {}",
                    product.id,
                    page.route
                );
                ensure!(
                    !page.title.trim().is_empty()
                        && !page.summary.trim().is_empty()
                        && !page.group.trim().is_empty(),
                    "{}:{} has incomplete page metadata",
                    product.id,
                    page.id
                );
                validate_typst_symbol(&page.symbol)?;
                validate_page_route(&page.route)?;
                self.require_file(&page.source)?;
                root_pages += usize::from(page.route.is_empty());
                tutorial_pages += usize::from(page.group == "Tutorial");
                manual_pages += usize::from(page.group == "Manual");
            }
            ensure!(
                root_pages == 1,
                "{} must have exactly one root page",
                product.id
            );
            ensure!(
                tutorial_pages >= 1,
                "{} must have at least one tutorial page",
                product.id
            );
            ensure!(
                manual_pages >= 2,
                "{} must have a manual with subpages",
                product.id
            );
            if let Some(path) = &product.changelog {
                self.require_file(path)?;
            }
            for related in &product.related {
                ensure!(
                    ids.contains(related.as_str()),
                    "{} links unknown product {related}",
                    product.id
                );
                ensure!(
                    related != &product.id,
                    "{} cannot relate to itself",
                    product.id
                );
            }

            ensure!(
                !product.rust_components.is_empty(),
                "{} has no Rust reference component",
                product.id
            );
            ensure!(
                !product.python_components.is_empty(),
                "{} has no Python reference component",
                product.id
            );

            for component in product
                .rust_components
                .iter()
                .chain(&product.python_components)
            {
                validate_route_segment(&component.id)?;
                ensure!(
                    component_ids.insert(format!("{}:{}", product.id, component.id)),
                    "duplicate component {} in {}",
                    component.id,
                    product.id
                );
                let version = self.component_version(component)?;
                ensure!(
                    !version.trim().is_empty(),
                    "{}:{} has an empty version",
                    product.id,
                    component.id
                );
            }
            for component in &product.rust_components {
                for feature in self.component_default_features(component)? {
                    ensure!(
                        component.features.contains(&feature),
                        "{} exact Rust reference matrix must list default feature {feature}",
                        component.package
                    );
                }
                self.component_catalog(product, ApiLanguage::Rust, component)?;
            }
            for component in &product.python_components {
                self.component_catalog(product, ApiLanguage::Python, component)?;
            }
        }
        self.generated_references()?;
        Ok(())
    }

    pub fn build(&self, request: BuildRequest) -> Result<()> {
        self.check()?;
        let tag = match request.channel {
            BuildChannel::Latest => {
                ensure!(
                    request.snapshot_tag.is_none(),
                    "--snapshot-tag is only valid with --channel snapshot"
                );
                None
            }
            BuildChannel::Snapshot => {
                let tag = request
                    .snapshot_tag
                    .as_deref()
                    .context("--snapshot-tag or GITHUB_REF_NAME is required for snapshots")?;
                self.validate_snapshot_tag(tag)?;
                Some(tag)
            }
        };

        let output = absolute_from(&self.root, &request.output);
        ensure_safe_output(&self.root, &output)?;
        fs::create_dir_all(&output)
            .wrap_err_with(|| format!("failed to create {}", output.display()))?;
        clear_stale_staging(&output)?;

        let selected = if request.product == "all" {
            self.registry.product.iter().collect::<Vec<_>>()
        } else {
            vec![
                self.registry
                    .product
                    .iter()
                    .find(|product| product.id == request.product)
                    .wrap_err_with(|| format!("unknown product {}", request.product))?,
            ]
        };

        let options = ProductBuildOptions {
            channel: request.channel,
            tag,
            output: &output,
            include_typst: request.include_typst,
            include_rustdoc: request.include_rustdoc,
            rustdoc_target_root: request.rustdoc_target_root.as_deref(),
            dependency_output: request.dependency_output.as_deref(),
        };
        for product in selected {
            self.build_product(product, &options)?;
        }

        self.write_developer_docs(&output)?;
        self.write_portal(&output, request.channel, tag)?;
        if request.product == "all" {
            self.validate_generated_links(&output, request.include_rustdoc)?;
        }
        Ok(())
    }

    fn build_product(
        &self,
        product: &ProductConfig,
        options: &ProductBuildOptions<'_>,
    ) -> Result<()> {
        let channel_path = match options.channel {
            BuildChannel::Latest => PathBuf::from("latest"),
            BuildChannel::Snapshot => {
                PathBuf::from("snapshots").join(options.tag.expect("validated tag"))
            }
        };
        let destination = options
            .output
            .join("products")
            .join(&product.id)
            .join(&channel_path);
        let staging_root = options.output.join(".staging");
        fs::create_dir_all(&staging_root)?;
        let staging = TempDirBuilder::new()
            .prefix(&format!("{}-", product.id))
            .tempdir_in(&staging_root)?;
        let site = staging.path().join("site");
        fs::create_dir_all(&site)?;
        self.write_site_assets(&site)?;

        let metadata = self.metadata(product, options.channel, options.tag, &channel_path)?;
        fs::write(
            site.join("snapshot.json"),
            serde_json::to_vec_pretty(&metadata)?,
        )?;
        let channel_label = match metadata.channel {
            BuildChannel::Latest => "latest",
            BuildChannel::Snapshot => "snapshot",
        };
        fs::write(
            site.join(".note"),
            format!(
                "alphal00p documentation\nproduct={}\nchannel={}\ntag={}\ncommit={}\nroute={}\n",
                product.id,
                channel_label,
                metadata.snapshot_tag.unwrap_or(""),
                metadata.git_commit,
                metadata.route,
            ),
        )?;
        self.write_component_catalogs(product, &site)?;
        self.write_generated_reference(product, &site)?;
        self.write_changelog_source(product, &site)?;

        if options.include_typst {
            self.render_typst(product, &metadata, &site, options.dependency_output)?;
        } else {
            self.write_fallback_page(product, &metadata, &site)?;
        }
        self.write_python_reference(product, &metadata, &site)?;
        if options.include_rustdoc {
            self.build_rustdoc(product, &site, options.rustdoc_target_root)?;
        } else {
            self.write_rustdoc_placeholder(product, &site)?;
        }
        self.write_rust_reference(product, &metadata, &site)?;
        self.write_reference_hub(product, &site)?;
        self.write_search_index(product, &site)?;
        self.decorate_site_pages(product, &metadata, &site)?;

        if options.channel == BuildChannel::Snapshot && destination.exists() {
            ensure!(
                directories_equal(&site, &destination)?,
                "immutable snapshot differs from {}",
                destination.display()
            );
        } else {
            if destination.exists() {
                fs::remove_dir_all(&destination)
                    .wrap_err_with(|| format!("failed to replace {}", destination.display()))?;
            }
            fs::create_dir_all(
                destination
                    .parent()
                    .context("product destination has no parent")?,
            )?;
            copy_tree(&site, &destination)?;
        }

        if options.channel == BuildChannel::Latest {
            self.write_product_redirect(product, options.output)?;
        }
        Ok(())
    }

    fn metadata<'a>(
        &'a self,
        product: &'a ProductConfig,
        channel: BuildChannel,
        tag: Option<&'a str>,
        channel_path: &Path,
    ) -> Result<SnapshotMetadata<'a>> {
        let mut components = Vec::new();
        for owner in &self.registry.product {
            for component in owner.rust_components.iter().chain(&owner.python_components) {
                components.push(ComponentVersion {
                    product: &owner.id,
                    id: &component.id,
                    package: &component.package,
                    kind: &component.kind,
                    module: component.module.as_deref(),
                    version: self.component_version(component)?,
                    features: &component.features,
                });
            }
        }
        let route = format!(
            "/products/{}/{}/",
            product.id,
            channel_path.to_string_lossy()
        );
        let git_commit = self.git_commit();
        ensure!(
            git_commit != "unknown",
            "cannot determine the documented Git commit; set ALPHAL00P_DOCS_GIT_COMMIT"
        );
        let git_timestamp = self.git_timestamp();
        ensure!(
            git_timestamp > 0,
            "cannot determine the documented commit timestamp; set SOURCE_DATE_EPOCH"
        );
        Ok(SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel,
            snapshot_tag: tag,
            git_commit,
            git_timestamp,
            route,
            components,
        })
    }

    fn write_site_assets(&self, destination: &Path) -> Result<()> {
        let assets = destination.join("assets");
        fs::create_dir_all(&assets)?;
        for (source, name) in [
            ("docs/assets/site.css", "site.css"),
            ("docs/assets/site.js", "site.js"),
            (
                "docs/assets/local-unitarity-light.svg",
                "local-unitarity-light.svg",
            ),
            (
                "docs/assets/local-unitarity-dark.svg",
                "local-unitarity-dark.svg",
            ),
            ("docs/assets/spensologo.svg", "spensologo.svg"),
            ("assets/gammalooplogo-light.svg", "gammalooplogo-light.svg"),
            ("assets/gammalooplogo-dark.svg", "gammalooplogo-dark.svg"),
        ] {
            fs::copy(self.root.join(source), assets.join(name))
                .wrap_err_with(|| format!("failed to copy documentation asset {source}"))?;
        }
        Ok(())
    }

    fn write_portal_assets(&self, destination: &Path) -> Result<()> {
        let people = destination.join("assets/people");
        fs::create_dir_all(&people)?;
        for person in &self.portal.people {
            if let Some(portrait) = &person.portrait {
                fs::copy(
                    self.root.join("docs/assets/people").join(portrait),
                    people.join(portrait),
                )
                .wrap_err_with(|| format!("failed to copy portrait for {}", person.name))?;
            }
        }
        fs::copy(
            self.root.join("docs/data/publications.json"),
            destination.join("assets/publications.json"),
        )?;
        Ok(())
    }

    fn write_component_catalogs(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let destination = site.join("catalogs");
        fs::create_dir_all(&destination)?;
        for (language, component) in product
            .rust_components
            .iter()
            .map(|component| (ApiLanguage::Rust, component))
            .chain(
                product
                    .python_components
                    .iter()
                    .map(|component| (ApiLanguage::Python, component)),
            )
        {
            self.export_component_catalog(
                product,
                language,
                component,
                &destination.join(format!("{}.json", component.id)),
            )?;
        }
        Ok(())
    }

    fn write_changelog_source(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let Some(changelog) = &product.changelog else {
            return Ok(());
        };
        let destination = site.join("reference/content");
        fs::create_dir_all(&destination)?;
        fs::copy(self.root.join(changelog), destination.join("changelog.md"))?;
        Ok(())
    }

    fn changelog_typst(&self, product: &ProductConfig) -> Result<String> {
        let Some(changelog) = &product.changelog else {
            return Ok(String::new());
        };
        let source = fs::read_to_string(self.root.join(changelog))?;
        let mut rendered = format!(
            "= Canonical package changelog <canonical-package-changelog>\n\
             The following release history is imported as typed Markdown from \
             #link(\"reference/content/changelog.md\")[#raw(\"{}\")].\n\n",
            typst_string(&changelog.to_string_lossy())
        );
        rendered.push_str(&markdown_changelog_to_typst(&source));
        Ok(rendered)
    }

    fn generated_references(&self) -> Result<(GammaLoopReference, VakintReference)> {
        let gamma_path = self.api_root.join("generated/gammaloop-reference.json");
        let gamma = serde_json::from_slice::<GammaLoopReference>(&fs::read(&gamma_path)?)
            .wrap_err_with(|| format!("failed to parse {}", gamma_path.display()))?;
        ensure!(
            gamma.schema_version == GENERATED_REFERENCE_SCHEMA,
            "GammaLoop generated reference schema drifted"
        );
        ensure!(
            !gamma.commands.is_empty() && !gamma.settings.is_empty(),
            "GammaLoop generated reference is empty"
        );
        validate_generated_anchors(
            gamma
                .commands
                .iter()
                .map(|command| generated_anchor("command", &command.path)),
        )?;
        validate_generated_anchors(
            gamma
                .settings
                .iter()
                .map(|setting| generated_anchor("setting", &setting.path)),
        )?;
        for source in &gamma.sources {
            self.require_file(Path::new(&source.path))?;
        }

        let vakint_path = self.api_root.join("generated/vakint-reference.json");
        let vakint = serde_json::from_slice::<VakintReference>(&fs::read(&vakint_path)?)
            .wrap_err_with(|| format!("failed to parse {}", vakint_path.display()))?;
        ensure!(
            vakint.schema_version == GENERATED_REFERENCE_SCHEMA,
            "Vakint generated reference schema drifted"
        );
        ensure!(
            !vakint.dependencies.is_empty() && !vakint.topologies.is_empty(),
            "Vakint generated reference is empty"
        );
        validate_generated_anchors(
            vakint
                .topologies
                .iter()
                .map(|topology| generated_anchor("topology", &topology.name)),
        )?;
        for source in &vakint.sources {
            self.require_file(Path::new(&source.path))?;
        }
        Ok((gamma, vakint))
    }

    fn write_generated_reference(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let destination = site.join("reference/generated");
        match product.id.as_str() {
            "gammaloop" => {
                let (reference, _) = self.generated_references()?;
                fs::create_dir_all(&destination)?;
                fs::copy(
                    self.api_root.join("generated/gammaloop-reference.json"),
                    destination.join("gammaloop-reference.json"),
                )?;
                fs::write(
                    destination.join("index.html"),
                    generated_reference_redirect(&product.title, "../cli/", "CLI reference"),
                )?;
                let first_class = site.join("reference/cli");
                fs::create_dir_all(&first_class)?;
                fs::write(
                    first_class.join("index.html"),
                    render_gammaloop_generated_reference(&product.title, &reference),
                )?;
            }
            "vakint" => {
                let (_, reference) = self.generated_references()?;
                fs::create_dir_all(&destination)?;
                fs::copy(
                    self.api_root.join("generated/vakint-reference.json"),
                    destination.join("vakint-reference.json"),
                )?;
                fs::write(
                    destination.join("index.html"),
                    generated_reference_redirect(
                        &product.title,
                        "../topologies/",
                        "topology reference",
                    ),
                )?;
                let first_class = site.join("reference/topologies");
                fs::create_dir_all(&first_class)?;
                fs::write(
                    first_class.join("index.html"),
                    render_vakint_generated_reference(&product.title, &reference),
                )?;
            }
            _ => {}
        }
        Ok(())
    }

    fn component_catalog(
        &self,
        product: &ProductConfig,
        language: ApiLanguage,
        component: &ComponentConfig,
    ) -> Result<DocCatalog> {
        let target = self.root.join("target");
        fs::create_dir_all(&target)?;
        let temporary = TempDirBuilder::new()
            .prefix("alphal00p-catalog-check-")
            .tempdir_in(target)?;
        let output = temporary.path().join(format!("{}.json", component.id));
        self.export_component_catalog(product, language, component, &output)?;
        serde_json::from_slice(&fs::read(&output)?)
            .wrap_err_with(|| format!("failed to merge isolated catalog {}", component.id))
    }

    fn export_component_catalog(
        &self,
        product: &ProductConfig,
        language: ApiLanguage,
        component: &ComponentConfig,
        output: &Path,
    ) -> Result<()> {
        let language = match language {
            ApiLanguage::Rust => "rust",
            ApiLanguage::Python => "python",
            _ => bail!("unsupported component catalog language {language:?}"),
        };
        let mut command = Command::new(env::var_os("CARGO").unwrap_or_else(|| "cargo".into()));
        command
            .current_dir(&self.root)
            .args(["run", "--locked", "--quiet"]);
        if let Some(profile) = env::var_os("ALPHAL00P_DOCS_CARGO_PROFILE") {
            command.arg("--profile").arg(profile);
        }
        command.args([
            "-p",
            "alphal00p-docs-catalogs",
            "--bin",
            "alphal00p-docs-catalogs",
            "--",
            "--product",
            &product.id,
            "--product-title",
            &product.title,
            "--component",
            &component.id,
            "--package",
            &component.package,
            "--component-title",
            &format!("{} {}", product.title, component.id),
            "--version",
            &self.component_version(component)?,
            "--language",
            language,
            "--workspace-root",
        ]);
        command.arg(&self.root);
        if let Some(module) = &component.module {
            command.arg("--module").arg(module);
        }
        for feature in &component.features {
            command.arg("--feature").arg(feature);
        }
        if language == "python" {
            let stub = self.python_stub_source(component).wrap_err_with(|| {
                format!(
                    "missing checked-in Python stub docs/api/python/{}.pyi",
                    component.id
                )
            })?;
            command
                .arg("--stub")
                .arg(stub)
                .arg("--source-file")
                .arg(format!("docs/api/python/{}.pyi", component.id));
        }
        let status = command
            .arg("--output")
            .arg(output)
            .status()
            .wrap_err_with(|| format!("failed to launch catalog exporter {}", component.id))?;
        ensure!(
            status.success(),
            "isolated catalog exporter failed for {}",
            component.id
        );
        Ok(())
    }

    fn catalog_reference_typst(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
    ) -> Result<String> {
        let mut rendered = String::from(
            "= API reference <supported-api-catalog>\n\
             The packages available in this version are listed below. Follow an item link for its complete language reference.\n",
        );
        for component in product
            .rust_components
            .iter()
            .chain(&product.python_components)
        {
            // Everything beyond this boundary consumes the neutral catalog
            // artifact. In particular, no PyO3 inventory is ever linked into
            // the renderer process: each isolated exporter has already emitted
            // its checked stub, and the per-component adapter has serialized it.
            let path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&path)?)
                .wrap_err_with(|| format!("failed to merge {}", path.display()))?;
            rendered.push_str(&format!("\n== #raw(\"{}\")\n", typst_string(&component.id)));
            append_scope_typst(&catalog, &catalog.root, &metadata.git_commit, &mut rendered);
        }
        Ok(rendered)
    }

    fn render_typst(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
        dependency_output: Option<&Path>,
    ) -> Result<()> {
        let target = self.root.join("target");
        fs::create_dir_all(&target)?;
        let work = TempDirBuilder::new()
            .prefix("alphal00p-typst-")
            .tempdir_in(target)?;
        let wrapper = work.path().join("site.typ");
        let source = product.source.to_string_lossy().replace('\\', "/");
        let title = typst_string(&product.title);
        let provenance = provenance_typst(metadata);
        let catalog_reference = self.catalog_reference_typst(product, metadata, site)?;
        let canonical_changelog = self.changelog_typst(product)?;
        let (gamma_reference, vakint_reference) = self.generated_references()?;
        let generated_reference =
            generated_reference_typst(&product.id, &gamma_reference, &vakint_reference);
        let mut page_imports = String::new();
        let mut page_documents = String::new();
        for (index, page) in product.pages.iter().enumerate() {
            let page_source = page.source.to_string_lossy().replace('\\', "/");
            let page_title = typst_string(&format!("{} · {}", product.title, page.title));
            page_imports.push_str(&format!("#import \"/{page_source}\" as chapter{index}\n"));
            page_documents.push_str(&format!(
                "#document(\"{}index.html\", title: [{page_title}])[#chapter{index}.{}]\n",
                page.route, page.symbol
            ));
        }
        let mut document_sections = vec!["#manual", "#catalog-reference"];
        if !generated_reference.is_empty() {
            document_sections.push("#generated-reference");
        }
        if !canonical_changelog.is_empty() {
            document_sections.push("#canonical-changelog");
        }
        document_sections.push("#provenance");
        let document_body = document_sections.join(" #pagebreak() ");
        fs::write(
            &wrapper,
            format!(
                "#import \"/{source}\": manual\n\
                 {page_imports}\
                 #let catalog-reference = [{catalog_reference}]\n\
                 #let generated-reference = [{generated_reference}]\n\
                 #let canonical-changelog = [{canonical_changelog}]\n\
                 #let provenance = [{provenance}]\n\
                 {page_documents}\
                 #document(\"manual.pdf\", title: [{title}])[{document_body}]\n"
            ),
        )?;
        let bundle = work.path().join("bundle");
        let metadata_json = serde_json::to_string(metadata)?;
        let mut command = Command::new("typst");
        command.current_dir(&self.root).args([
            "compile",
            "--root",
            self.root.to_str().context("workspace root is not UTF-8")?,
            "--features",
            "html,bundle",
            "--format",
            "bundle",
            "--creation-timestamp",
            &metadata.git_timestamp.to_string(),
            "--input",
            &format!("metadata={metadata_json}"),
            "--input",
            &format!("git-commit={}", metadata.git_commit),
            "--input",
            match metadata.channel {
                BuildChannel::Latest => "channel=latest",
                BuildChannel::Snapshot => "channel=snapshot",
            },
            "--input",
            &format!("snapshot-tag={}", metadata.snapshot_tag.unwrap_or_default()),
        ]);
        if let Some(directory) = dependency_output {
            fs::create_dir_all(directory)?;
            command
                .args(["--deps-format", "zero", "--deps"])
                .arg(directory.join(format!("{}.deps", product.id)));
        }
        let status = command
            .arg(&wrapper)
            .arg(&bundle)
            .status()
            .wrap_err("failed to launch Typst 0.15; enter the Nix development shell")?;
        ensure!(status.success(), "Typst failed for {}", product.id);
        ensure!(
            bundle.join("index.html").is_file(),
            "Typst emitted no index.html"
        );
        ensure!(
            bundle.join("manual.pdf").is_file(),
            "Typst emitted no manual.pdf"
        );
        for page in &product.pages {
            ensure!(
                bundle.join(&page.route).join("index.html").is_file(),
                "Typst emitted no page {} for {}",
                page.route,
                product.id
            );
        }
        copy_tree(&bundle, site)
    }

    fn build_rustdoc(
        &self,
        product: &ProductConfig,
        site: &Path,
        persistent_target_root: Option<&Path>,
    ) -> Result<()> {
        let target_root = self.root.join("target");
        fs::create_dir_all(&target_root)?;
        let temporary_target;
        let target = if let Some(root) = persistent_target_root {
            // Share Cargo's compiled artifacts across products and watch
            // generations. The Rustdoc output itself is cleared per product so
            // each published sidecar contains only that product's root crates.
            prepare_persistent_rustdoc_target(&self.root, root)?
        } else {
            temporary_target = TempDirBuilder::new()
                .prefix(&format!("alphal00p-rustdoc-{}-", product.id))
                .tempdir_in(target_root)?;
            temporary_target.path().to_path_buf()
        };
        let rustdoc = target.join("doc");
        if persistent_target_root.is_some() {
            clear_persistent_rustdoc_output(&target)?;
        }
        let rustdoc_flags = env::var("RUSTDOCFLAGS").map_or_else(
            |_| STRICT_RUSTDOC_FLAGS.to_owned(),
            |flags| format!("{flags} {STRICT_RUSTDOC_FLAGS}"),
        );
        for component in &product.rust_components {
            let mut command = Command::new("cargo");
            command
                .current_dir(&self.root)
                .env("CARGO_TARGET_DIR", &target)
                .env("RUSTDOCFLAGS", &rustdoc_flags)
                .args([
                    "doc",
                    "--locked",
                    "--no-deps",
                    "--no-default-features",
                    "-p",
                    &component.package,
                ]);
            if !component.features.is_empty() {
                command.args(["--features", &component.features.join(",")]);
            }
            let status = command
                .status()
                .wrap_err_with(|| format!("failed to launch rustdoc for {}", component.package))?;
            ensure!(status.success(), "rustdoc failed for {}", component.package);
        }

        let mut packages = String::new();
        for component in &product.rust_components {
            let crate_name = component.package.replace('-', "_");
            ensure!(
                rustdoc.join(&crate_name).join("index.html").is_file(),
                "rustdoc emitted no crate index for {}",
                component.package
            );
            packages.push_str(&format!(
                "<li><a href=\"reference/rust/{crate_name}/index.html\"><code>{}</code></a></li>",
                escape_html(&component.package)
            ));
        }

        let destination = site.join("reference/rust");
        copy_tree(&rustdoc, &destination)?;
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Rust API",
                &format!(
                    "<p>Choose a crate to browse its complete Rust API.</p><ul>{packages}</ul>"
                ),
            ),
        )?;
        Ok(())
    }

    fn write_python_reference(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
    ) -> Result<()> {
        let destination = site.join("reference/python");
        fs::create_dir_all(&destination)?;
        let mut cards = String::new();
        for component in &product.python_components {
            let stub = self.python_stub_source(component);
            let stub_name = format!("{}.pyi", component.id);
            let stub_text = if let Some(path) = stub {
                let text = fs::read_to_string(&path)
                    .wrap_err_with(|| format!("failed to read {}", path.display()))?;
                fs::write(destination.join(&stub_name), &text)?;
                text
            } else {
                let module = component.module.as_deref().unwrap_or(&component.package);
                let text = format!("# Python interface for {module}\n");
                fs::write(destination.join(&stub_name), &text)?;
                text
            };
            let catalog_path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&catalog_path)?)
                .wrap_err_with(|| format!("failed to parse {}", catalog_path.display()))?;
            catalog.validate()?;
            let export_count = count_catalog_items(&catalog.root);
            let module = component.module.as_deref().unwrap_or(&component.package);
            cards.push_str(&format!(
                "<article class=\"api-component-card\"><h2><a href=\"reference/python/{}/\"><code>{}</code></a></h2><p>{} documented exports · version {}</p><p>{}</p></article>",
                escape_html(&component.id),
                escape_html(module),
                export_count,
                escape_html(&self.component_version(component)?),
                escape_html(&catalog.root.summary.clone().unwrap_or_else(|| {
                    format!("Supported Python surface for {}.", component.id)
                })),
            ));
            let component_page = destination.join(&component.id);
            fs::create_dir_all(&component_page)?;
            fs::write(
                component_page.join("index.html"),
                reference_page(
                    &product.title,
                    &format!("{} Python API", component.id),
                    &render_python_catalog(&catalog, &stub_name, &stub_text, &metadata.git_commit),
                ),
            )?;
        }
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Python API",
                &format!(
                    "<p>Choose a Python module to browse classes, functions, signatures, parameters, examples, and feature requirements. Type-checker <code>.pyi</code> files are also available for download.</p><div class=\"api-component-grid\">{cards}</div>"
                ),
            ),
        )?;
        Ok(())
    }

    fn write_rust_reference(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
    ) -> Result<()> {
        let destination = site.join("reference/rust");
        let supported = destination.join("supported");
        fs::create_dir_all(&supported)?;
        let mut cards = String::new();
        for component in &product.rust_components {
            let catalog_path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&catalog_path)?)
                .wrap_err_with(|| format!("failed to parse {}", catalog_path.display()))?;
            catalog.validate()?;
            let count = count_catalog_items(&catalog.root);
            let crate_name = component.package.replace('-', "_");
            cards.push_str(&format!(
                "<article class=\"api-component-card\"><h2><a href=\"reference/rust/supported/{}/\"><code>{}</code></a></h2><p>{count} curated entries · version {}</p><p>{}</p><a href=\"reference/rust/{crate_name}/\">Full Rustdoc</a></article>",
                escape_html(&component.id),
                escape_html(&component.package),
                escape_html(&self.component_version(component)?),
                escape_html(&catalog.root.summary.clone().unwrap_or_else(|| {
                    format!("Supported Rust surface for {}.", component.id)
                })),
            ));
            let mut body = format!(
                "<p><code>{}</code> · {count} curated entries · version {}</p><p><a href=\"reference/rust/{crate_name}/\">Open the complete generated Rustdoc</a>. The source-backed catalog below highlights the supported public surface with signatures, feature gates, examples, and source links.</p>",
                escape_html(&component.package),
                escape_html(&catalog.component.version),
            );
            let mut path = Vec::new();
            render_python_scope(
                &catalog,
                &catalog.root,
                &mut path,
                &metadata.git_commit,
                2,
                &mut body,
            );
            let component_page = supported.join(&component.id);
            fs::create_dir_all(&component_page)?;
            fs::write(
                component_page.join("index.html"),
                reference_page(
                    &product.title,
                    &format!("{} supported Rust API", component.id),
                    &body,
                ),
            )?;
        }
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Rust API",
                &format!(
                    "<p>Browse curated, source-backed API catalogs first; use the complete Rustdoc sidecars for implementation-level detail.</p><div class=\"api-component-grid\">{cards}</div>"
                ),
            ),
        )?;
        Ok(())
    }

    fn write_reference_hub(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let mut cards = format!(
            "<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Curated + exhaustive</p><h2><a href=\"reference/rust/\">Rust API</a></h2><p>{} source-backed components with complete Rustdoc sidecars.</p></article><article class=\"reference-hub-card\"><p class=\"portal-kicker\">Generated signatures</p><h2><a href=\"reference/python/\">Python API</a></h2><p>{} modules rendered from checked runtime inventories; stubs remain downloadable for type checkers.</p></article>",
            product.rust_components.len(),
            product.python_components.len(),
        );
        if let Some((route, title)) = supplemental_reference(&product.id) {
            cards.push_str(&format!(
                "<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Generated from the implementation</p><h2><a href=\"{route}\">{}</a></h2><p>Version-specific reference data, also available as JSON for tooling.</p></article>",
                escape_html(title),
            ));
        }
        let destination = site.join("reference");
        fs::create_dir_all(&destination)?;
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Reference",
                &format!(
                    "<p>Generated catalogs make the supported interfaces searchable without asking raw source files or Rustdoc to serve as the manual.</p><div class=\"reference-hub-grid\">{cards}</div>"
                ),
            ),
        )?;
        Ok(())
    }

    fn python_stub_source(&self, component: &ComponentConfig) -> Option<PathBuf> {
        let checked_in = self
            .api_root
            .join("python")
            .join(format!("{}.pyi", component.id));
        if checked_in.is_file() {
            return Some(checked_in);
        }
        if component.package == "linnet-py" {
            let path = self.root.join("crates/linnet-py/linnet_py.pyi");
            if path.is_file() {
                return Some(path);
            }
        }
        None
    }

    fn write_rustdoc_placeholder(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let destination = site.join("reference/rust");
        fs::create_dir_all(&destination)?;
        let packages = product
            .rust_components
            .iter()
            .map(|component| format!("<li><code>{}</code></li>", escape_html(&component.package)))
            .collect::<String>();
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Rust API",
                &format!(
                    "<p>Rustdoc generation was skipped for this build.</p><ul>{packages}</ul>"
                ),
            ),
        )?;
        Ok(())
    }

    fn write_fallback_page(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
    ) -> Result<()> {
        let versions = metadata
            .components
            .iter()
            .map(|component| {
                format!(
                    "<li><code>{}</code> {}</li>",
                    escape_html(component.package),
                    escape_html(&component.version)
                )
            })
            .collect::<String>();
        fs::write(
            site.join("index.html"),
            format!(
                "<!doctype html><html><head><meta charset=\"utf-8\"><title>{}</title></head><body><main><h1>{}</h1><p>{}</p><ul>{versions}</ul></main></body></html>",
                escape_html(&product.title),
                escape_html(&product.title),
                escape_html(&product.tagline),
            ),
        )?;
        Ok(())
    }

    fn decorate_site_pages(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
    ) -> Result<()> {
        let mut pages = product
            .pages
            .iter()
            .map(|page| SitePage {
                route: page.route.clone(),
                title: page.title.clone(),
                group: page.group.clone(),
            })
            .collect::<Vec<_>>();
        pages.push(SitePage::new("reference/", "Reference", "Reference"));
        pages.push(SitePage::new(
            "reference/python/",
            "Python API",
            "Reference",
        ));
        for component in &product.python_components {
            pages.push(SitePage::new(
                format!("reference/python/{}/", component.id),
                format!("{} Python API", component.id),
                "Reference",
            ));
        }
        pages.push(SitePage::new("reference/rust/", "Rust API", "Reference"));
        for component in &product.rust_components {
            pages.push(SitePage::new(
                format!("reference/rust/supported/{}/", component.id),
                format!("{} supported Rust API", component.id),
                "Reference",
            ));
        }
        if let Some((route, title)) = supplemental_reference(&product.id) {
            pages.push(SitePage::new(route, title, "Reference"));
        }
        pages.retain(|page| site.join(&page.route).join("index.html").is_file());

        for (index, page) in pages.iter().enumerate() {
            self.decorate_html_page(
                product,
                metadata,
                site,
                page,
                index
                    .checked_sub(1)
                    .and_then(|previous| pages.get(previous)),
                pages.get(index + 1),
            )?;
        }
        Ok(())
    }

    fn decorate_html_page(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
        page: &SitePage,
        previous: Option<&SitePage>,
        next: Option<&SitePage>,
    ) -> Result<()> {
        let path = site.join(&page.route).join("index.html");
        let source = fs::read_to_string(&path)?;
        let mut body = extract_html_body(&source)?.to_owned();
        body = body.replace("<nav><a href=\"../../\">Manual</a></nav>", "");
        if page.route.is_empty() {
            body = format!("{}{body}", product_hero(product));
        } else {
            promote_heading_hierarchy(&mut body);
        }
        let docs_root = page_root_prefix(&page.route);
        body = rewrite_page_links(&body, &docs_root)?;
        let (body, headings) = inject_heading_ids(&body)?;
        let toc = headings
            .iter()
            .filter(|heading| matches!(heading.level, 2 | 3))
            .map(|heading| {
                format!(
                    "<a class=\"toc-link\" data-level=\"{}\" href=\"#{}\">{}</a>",
                    heading.level,
                    escape_html(&heading.id),
                    escape_html(&heading.title)
                )
            })
            .collect::<String>();
        let sidebar = self.site_sidebar(product, metadata, page, &docs_root);
        let product_options = self.product_options(product, metadata, &docs_root);
        let page_navigation = render_page_navigation(previous, next, &docs_root);
        let product_root = if docs_root.is_empty() {
            "./"
        } else {
            &docs_root
        };
        let portal = match metadata.channel {
            BuildChannel::Latest => format!("{docs_root}../../../"),
            BuildChannel::Snapshot => format!("{docs_root}../../../../"),
        };
        let version = match metadata.channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => {
                format!("snapshot {}", metadata.snapshot_tag.unwrap_or("unknown"))
            }
        };
        let toc_markup = if toc.is_empty() {
            String::new()
        } else {
            format!(
                "<aside class=\"docs-toc\" aria-label=\"On this page\"><p class=\"toc-title\">On this page</p>{toc}</aside>"
            )
        };
        let favicon = favicon_links(&format!("{docs_root}assets/"));
        let html = format!(
            "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\">{favicon}<title>{} · {}</title><link rel=\"stylesheet\" href=\"{}assets/site.css\"><script defer src=\"{}assets/site.js\"></script></head><body data-docs-root=\"{}\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"site-header\"><button class=\"header-button menu-button\" type=\"button\" data-menu-toggle aria-label=\"Open navigation\" aria-controls=\"docs-sidebar\" aria-expanded=\"false\">☰</button><a class=\"site-brand\" href=\"{}\"><span class=\"site-brand-mark\">α</span><span class=\"site-brand-name\">αLoop docs</span></a><div class=\"site-header-tools\"><select class=\"product-select\" data-product-select aria-label=\"Select project\">{}</select><button class=\"header-button\" type=\"button\" data-search-open>Search <span class=\"header-button-label\">⌘K</span></button><button class=\"header-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\">◐</button></div></header><div class=\"docs-shell\">{sidebar}<main class=\"docs-main\" id=\"main-content\"><nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"{}\">{}</a> / {} / {}</nav><article class=\"docs-article\">{body}</article>{page_navigation}<footer class=\"page-footer\">{} · <a href=\"{}manual.pdf\">Complete PDF manual</a> · <a href=\"https://github.com/alphal00p/gammaloop/tree/{}\">source {}</a></footer></main>{toc_markup}</div><button class=\"sidebar-backdrop\" type=\"button\" data-sidebar-backdrop aria-label=\"Close navigation\"></button><dialog class=\"search-dialog\" data-search-dialog><form class=\"search-form\" method=\"dialog\"><input class=\"search-input\" type=\"search\" data-search-input placeholder=\"Search this project\" aria-label=\"Search documentation\"><button class=\"header-button\" value=\"close\">Close</button></form><ul class=\"search-results\" data-search-results aria-live=\"polite\"></ul></dialog></body></html>",
            escape_html(&format!("{} — {}", product.tagline, page.title)),
            escape_html(&product.title),
            escape_html(&page.title),
            escape_html(&docs_root),
            escape_html(&docs_root),
            escape_html(&docs_root),
            escape_html(&portal),
            product_options,
            escape_html(product_root),
            escape_html(&product.title),
            escape_html(&page.group),
            escape_html(&page.title),
            escape_html(&version),
            escape_html(&docs_root),
            escape_html(&metadata.git_commit),
            escape_html(&metadata.git_commit.chars().take(12).collect::<String>()),
        );
        fs::write(path, html)?;
        Ok(())
    }

    fn site_sidebar(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        current: &SitePage,
        docs_root: &str,
    ) -> String {
        let mut groups = BTreeMap::<String, Vec<(String, String)>>::new();
        let mut group_order = Vec::new();
        for page in &product.pages {
            if !groups.contains_key(&page.group) {
                group_order.push(page.group.clone());
            }
            groups
                .entry(page.group.clone())
                .or_default()
                .push((page.route.clone(), page.title.clone()));
        }
        group_order.push("Reference".to_owned());
        let reference = groups.entry("Reference".to_owned()).or_default();
        reference.push(("reference/".to_owned(), "Reference overview".to_owned()));
        reference.push(("reference/python/".to_owned(), "Python API".to_owned()));
        reference.push(("reference/rust/".to_owned(), "Rust API".to_owned()));
        if let Some((route, title)) = supplemental_reference(&product.id) {
            reference.push((route.to_owned(), title.to_owned()));
        }

        let mut navigation = String::new();
        for group in group_order {
            let links = groups
                .get(&group)
                .into_iter()
                .flatten()
                .map(|(route, title)| {
                    let href = if route.is_empty() {
                        if docs_root.is_empty() {
                            "./".to_owned()
                        } else {
                            docs_root.to_owned()
                        }
                    } else {
                        format!("{docs_root}{route}")
                    };
                    format!(
                        "<a class=\"sidebar-link\" href=\"{}\"{}>{}</a>",
                        escape_html(&href),
                        if route == &current.route
                            || (matches!(route.as_str(), "reference/python/" | "reference/rust/")
                                && current.route.starts_with(route)
                                && current.route != "reference/")
                        {
                            " aria-current=\"page\""
                        } else {
                            ""
                        },
                        escape_html(title)
                    )
                })
                .collect::<String>();
            navigation.push_str(&format!(
                "<section class=\"sidebar-group\"><p class=\"sidebar-group-title\">{}</p>{links}</section>",
                escape_html(&group)
            ));
        }
        if metadata.channel == BuildChannel::Latest {
            let portal = format!("{docs_root}../../../");
            navigation.push_str(&format!(
                "<section class=\"sidebar-group sidebar-developer-group\"><p class=\"sidebar-group-title\">For developers</p><a class=\"sidebar-link\" href=\"{}developers/\">Architecture &amp; engineering notes</a></section>",
                escape_html(&portal),
            ));
        }
        let version = match metadata.channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => metadata.snapshot_tag.unwrap_or("unknown").to_owned(),
        };
        format!(
            "<aside class=\"docs-sidebar\" id=\"docs-sidebar\" aria-label=\"{} manual\">{navigation}<p class=\"sidebar-meta\"><strong>{}</strong><br><code>{}</code><br><a href=\"{}manual.pdf\">Download PDF</a></p></aside>",
            escape_html(&product.title),
            escape_html(&product.title),
            escape_html(&version),
            escape_html(docs_root)
        )
    }

    fn product_options(
        &self,
        current: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        docs_root: &str,
    ) -> String {
        self.registry
            .product
            .iter()
            .map(|product| {
                let route = match metadata.channel {
                    BuildChannel::Latest => {
                        format!("{docs_root}../../{}/latest/", product.id)
                    }
                    BuildChannel::Snapshot => format!(
                        "{docs_root}../../../{}/snapshots/{}/",
                        product.id,
                        metadata.snapshot_tag.unwrap_or_default()
                    ),
                };
                format!(
                    "<option value=\"{}\"{}>{}</option>",
                    escape_html(&route),
                    if product.id == current.id {
                        " selected"
                    } else {
                        ""
                    },
                    escape_html(&product.title)
                )
            })
            .collect()
    }

    fn write_search_index(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let mut entries = vec![SearchEntry {
            title: product.title.clone(),
            summary: product.tagline.clone(),
            href: "index.html".to_owned(),
            kind: "product".to_owned(),
        }];
        for page in &product.pages {
            let page_href = if page.route.is_empty() {
                "index.html".to_owned()
            } else {
                page.route.clone()
            };
            entries.push(SearchEntry {
                title: page.title.clone(),
                summary: page.summary.clone(),
                href: page_href.clone(),
                kind: page.group.to_lowercase(),
            });
            let source = fs::read_to_string(self.root.join(&page.source))?;
            for line in source.lines() {
                let trimmed = line.trim();
                let level = trimmed
                    .chars()
                    .take_while(|character| *character == '=')
                    .count();
                if level >= 2
                    && let Some(title) = trimmed
                        .get(level..)
                        .and_then(|title| title.strip_prefix(' '))
                {
                    entries.push(SearchEntry {
                        title: title.trim().to_owned(),
                        summary: format!("{} · {}", product.title, page.title),
                        href: format!("{page_href}#{}", slug(title)),
                        kind: page.group.to_lowercase(),
                    });
                }
            }
        }
        for component in product
            .rust_components
            .iter()
            .chain(&product.python_components)
        {
            let path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&path)?)
                .wrap_err_with(|| format!("failed to parse {}", path.display()))?;
            append_catalog_search(&catalog, &catalog.root, &mut entries);
        }
        match product.id.as_str() {
            "gammaloop" => {
                let (reference, _) = self.generated_references()?;
                entries.extend(reference.commands.into_iter().map(|command| SearchEntry {
                    title: command.path.clone(),
                    summary: command.about,
                    href: format!(
                        "reference/cli/#{}",
                        generated_anchor("command", &command.path)
                    ),
                    kind: "command".to_owned(),
                }));
                entries.extend(reference.settings.into_iter().map(|setting| SearchEntry {
                    title: setting.path.clone(),
                    summary: setting.description,
                    href: format!(
                        "reference/cli/#{}",
                        generated_anchor("setting", &setting.path)
                    ),
                    kind: "setting".to_owned(),
                }));
            }
            "vakint" => {
                let (_, reference) = self.generated_references()?;
                entries.extend(
                    reference
                        .topologies
                        .into_iter()
                        .map(|topology| SearchEntry {
                            title: topology.name.clone(),
                            summary: format!(
                                "{}-loop Vakint topology with {} propagator slots",
                                topology.loops, topology.propagator_slots
                            ),
                            href: format!(
                                "reference/topologies/#{}",
                                generated_anchor("topology", &topology.name)
                            ),
                            kind: "topology".to_owned(),
                        }),
                );
                entries.extend(
                    reference
                        .dependencies
                        .into_iter()
                        .map(|dependency| SearchEntry {
                            title: dependency.name,
                            summary: format!(
                                "Minimum supported version {}",
                                dependency.minimum_version
                            ),
                            href: "reference/topologies/#external-dependencies".to_owned(),
                            kind: "dependency".to_owned(),
                        }),
                );
            }
            _ => {}
        }
        fs::write(
            site.join("search-index.json"),
            serde_json::to_vec_pretty(&entries)?,
        )?;
        Ok(())
    }

    fn write_developer_docs(&self, output: &Path) -> Result<()> {
        let staging_root = output.join(".staging");
        fs::create_dir_all(&staging_root)?;
        let staging = TempDirBuilder::new()
            .prefix("developers-")
            .tempdir_in(&staging_root)?;
        let developer_root = staging.path().join("developers");
        fs::create_dir_all(&developer_root)?;
        self.write_site_assets(&developer_root)?;

        let commit = self.git_commit();
        ensure!(
            commit != "unknown",
            "cannot determine the documented Git commit; set ALPHAL00P_DOCS_GIT_COMMIT"
        );
        let notes = self
            .developers
            .section
            .iter()
            .flat_map(|section| section.note.iter())
            .collect::<Vec<_>>();
        let note_routes = notes
            .iter()
            .map(|note| (note.source.clone(), note.id.clone()))
            .collect::<BTreeMap<_, _>>();
        let pages = self
            .developers
            .section
            .iter()
            .flat_map(|section| {
                section.note.iter().map(|note| {
                    SitePage::new(
                        format!("architecture/{}/", note.id),
                        &note.title,
                        &section.title,
                    )
                })
            })
            .collect::<Vec<_>>();
        let mut search = Vec::new();

        for (index, note) in notes.iter().enumerate() {
            let page = &pages[index];
            let docs_root = page_root_prefix(&page.route);
            let destination = developer_root.join(&page.route);
            fs::create_dir_all(&destination)?;
            let source = fs::read_to_string(self.root.join(&note.source))?;
            let source_url = format!(
                "https://github.com/alphal00p/gammaloop/blob/{}/{}",
                commit,
                note.source.to_string_lossy()
            );
            let body = match note
                .source
                .extension()
                .and_then(|extension| extension.to_str())
            {
                Some("md") => {
                    let markdown = strip_first_markdown_title(&source)?;
                    let rendered = render_developer_markdown(markdown);
                    rewrite_developer_source_links(
                        &rendered,
                        &note.source,
                        &commit,
                        &note_routes,
                        &self.root,
                    )?
                }
                Some("html") => {
                    validate_developer_html_note(&source)?;
                    let diagram = rewrite_developer_source_links(
                        &source,
                        &note.source,
                        &commit,
                        &note_routes,
                        &self.root,
                    )?;
                    let diagram = decorate_developer_diagram(&diagram, &note.title, &source_url)?;
                    fs::write(destination.join("diagram.html"), diagram)?;
                    format!(
                        "<div class=\"developer-diagram\"><iframe title=\"{}\" src=\"diagram.html\" loading=\"lazy\" sandbox></iframe><p><a class=\"hero-action\" href=\"diagram.html\">Open the full-page diagram <span aria-hidden=\"true\">↗</span></a></p></div>",
                        escape_html(&note.title)
                    )
                }
                _ => unreachable!("developer note extensions are validated"),
            };
            let status_class = slug(&note.status);
            let article = format!(
                "<header class=\"developer-note-hero\"><p class=\"product-eyebrow\">For developers</p><div class=\"developer-note-title\"><h1>{}</h1><div class=\"developer-note-badges\"><span class=\"developer-status developer-status-{}\">{}</span><span class=\"developer-kind\">{}</span></div></div><p>{}</p><p class=\"developer-source\"><a href=\"{}\">View the source note <span aria-hidden=\"true\">↗</span></a></p></header>{body}",
                escape_html(&note.title),
                escape_html(&status_class),
                escape_html(&note.status),
                escape_html(&note.kind),
                escape_html(&note.summary),
                escape_html(&source_url),
            );
            let (article, headings) = inject_heading_ids(&article)?;
            let toc = developer_toc(&headings);
            let inline_toc = developer_inline_toc(&headings);
            let sidebar = self.developer_sidebar(Some(&note.id), &docs_root);
            let page_navigation = render_page_navigation(
                index
                    .checked_sub(1)
                    .and_then(|previous| pages.get(previous)),
                pages.get(index + 1),
                &docs_root,
            )
            .replace("Manual pagination", "Developer note pagination");
            let portal_root = format!("{docs_root}../");
            let main = format!(
                "<nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"{}\">αLoop</a> / <a href=\"{}\">Developers</a> / {}</nav>{inline_toc}<article class=\"docs-article developer-article\">{article}</article>{page_navigation}<footer class=\"page-footer\">Developer architecture · <a href=\"{}\">source {}</a></footer>",
                escape_html(&portal_root),
                escape_html(&docs_root),
                escape_html(&note.title),
                escape_html(&source_url),
                escape_html(&commit.chars().take(12).collect::<String>()),
            );
            fs::write(
                destination.join("index.html"),
                developer_document(
                    &format!("{} · Developer architecture", note.title),
                    &note.summary,
                    &docs_root,
                    &portal_root,
                    &sidebar,
                    &main,
                    &toc,
                ),
            )?;

            search.push(SearchEntry {
                title: note.title.clone(),
                summary: note.summary.clone(),
                href: page.route.clone(),
                kind: format!("developer {}", note.status.to_lowercase()),
            });
            search.extend(
                headings
                    .into_iter()
                    .filter(|heading| matches!(heading.level, 2 | 3))
                    .map(|heading| SearchEntry {
                        title: heading.title,
                        summary: note.title.clone(),
                        href: format!("{}#{}", page.route, heading.id),
                        kind: "developer heading".to_owned(),
                    }),
            );
        }

        fs::write(
            developer_root.join("search-index.json"),
            serde_json::to_vec_pretty(&search)?,
        )?;
        self.write_developer_hub(&developer_root, &commit)?;
        fs::write(
            developer_root.join(".note"),
            format!(
                "alphal00p developer architecture\ncommit={commit}\nnotes={}\nroute=developers/\n",
                notes.len()
            ),
        )?;
        replace_generated_tree(&developer_root, &output.join("developers"))?;
        Ok(())
    }

    fn write_developer_hub(&self, developer_root: &Path, commit: &str) -> Result<()> {
        let sections = self
            .developers
            .section
            .iter()
            .map(|section| {
                let notes = section
                    .note
                    .iter()
                    .map(|note| {
                        format!(
                            "<article class=\"developer-card\"><div class=\"developer-card-meta\"><span class=\"developer-status developer-status-{}\">{}</span><span>{}</span></div><h3><a href=\"architecture/{}/\">{}</a></h3><p>{}</p><a class=\"portal-text-link\" href=\"architecture/{}/\">Read note <span aria-hidden=\"true\">→</span></a></article>",
                            escape_html(&slug(&note.status)),
                            escape_html(&note.status),
                            escape_html(&note.kind),
                            escape_html(&note.id),
                            escape_html(&note.title),
                            escape_html(&note.summary),
                            escape_html(&note.id),
                        )
                    })
                    .collect::<String>();
                format!(
                    "<section class=\"developer-section\" id=\"{}\"><div class=\"developer-section-heading\"><p class=\"portal-kicker\">Developer notes</p><h2>{}</h2><p>{}</p></div><div class=\"developer-card-grid\">{notes}</div></section>",
                    escape_html(&section.id),
                    escape_html(&section.title),
                    escape_html(&section.summary),
                )
            })
            .collect::<String>();
        let sidebar = self.developer_sidebar(None, "");
        let source_url =
            format!("https://github.com/alphal00p/gammaloop/tree/{commit}/docs/architecture");
        let main = format!(
            "<nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"../\">αLoop</a> / Developers</nav><article class=\"docs-article developer-article\"><header class=\"developer-hero\"><p class=\"product-eyebrow\">For developers</p><h1>{}</h1><p>{}</p><aside class=\"developer-audience\"><strong>Looking for usage documentation?</strong><span>Tutorials and manuals live with each research project. This area documents implementation details, design work, and engineering investigations for contributors.</span><a href=\"../#projects\">Browse project documentation <span aria-hidden=\"true\">→</span></a></aside></header>{sections}</article><footer class=\"page-footer\">{} classified notes · <a href=\"{}\">architecture sources at {}</a></footer>",
            escape_html(&self.developers.title),
            escape_html(&self.developers.summary),
            self.developers
                .section
                .iter()
                .map(|section| section.note.len())
                .sum::<usize>(),
            escape_html(&source_url),
            escape_html(&commit.chars().take(12).collect::<String>()),
        );
        fs::write(
            developer_root.join("index.html"),
            developer_document(
                &self.developers.title,
                &self.developers.summary,
                "",
                "../",
                &sidebar,
                &main,
                "",
            ),
        )?;
        Ok(())
    }

    fn developer_sidebar(&self, current: Option<&str>, docs_root: &str) -> String {
        let overview_current = if current.is_none() {
            " aria-current=\"page\""
        } else {
            ""
        };
        let mut navigation = format!(
            "<section class=\"sidebar-group\"><p class=\"sidebar-group-title\">Developer area</p><a class=\"sidebar-link\" href=\"{}\"{}>Overview</a></section>",
            if docs_root.is_empty() {
                "./"
            } else {
                docs_root
            },
            overview_current,
        );
        for section in &self.developers.section {
            let links = section
                .note
                .iter()
                .map(|note| {
                    format!(
                        "<a class=\"sidebar-link\" href=\"{}architecture/{}/\"{}>{}</a>",
                        escape_html(docs_root),
                        escape_html(&note.id),
                        if current == Some(note.id.as_str()) {
                            " aria-current=\"page\""
                        } else {
                            ""
                        },
                        escape_html(&note.title),
                    )
                })
                .collect::<String>();
            navigation.push_str(&format!(
                "<section class=\"sidebar-group\"><p class=\"sidebar-group-title\">{}</p>{links}</section>",
                escape_html(&section.title),
            ));
        }
        format!(
            "<aside class=\"docs-sidebar developer-sidebar\" id=\"docs-sidebar\" aria-label=\"Developer architecture\">{navigation}<p class=\"sidebar-meta\"><strong>For developers</strong><br>Implementation and engineering notes.<br><a href=\"{}../#projects\">Tutorials &amp; manuals</a></p></aside>",
            escape_html(docs_root),
        )
    }

    fn write_portal(&self, output: &Path, channel: BuildChannel, tag: Option<&str>) -> Result<()> {
        self.write_site_assets(output)?;
        self.write_portal_assets(output)?;
        let channel_route = match channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => format!(
                "snapshots/{}",
                tag.expect("snapshot tag was validated before rendering the portal")
            ),
        };
        let projects = self
            .registry
            .product
            .iter()
            .enumerate()
            .map(|(index, product)| {
                let packages = product
                    .rust_components
                    .iter()
                    .chain(&product.python_components)
                    .map(|component| component.package.as_str())
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .map(|package| {
                        format!(
                            "<span class=\"portal-package\">{}</span>",
                            escape_html(package)
                        )
                    })
                    .collect::<String>();
                let manual_route = product
                    .pages
                    .iter()
                    .find(|page| page.group == "Manual")
                    .map(|page| {
                        page.route
                            .trim_start_matches("manual/")
                            .trim_end_matches('/')
                    })
                    .unwrap_or_default();
                format!(
                    r#"<article class="portal-project-card"><div class="portal-project-meta"><span>{:02}</span><span>Research project</span></div><h3><a href="products/{}/{}/">{}</a></h3><p class="portal-project-summary">{}</p><div class="portal-packages" aria-label="{} crates and modules"><span class="portal-packages-label">Crates &amp; modules</span>{}</div><nav class="portal-card-links" aria-label="{} documentation"><a class="portal-card-primary" href="products/{}/{}/">Overview <span aria-hidden="true">↗</span></a><a href="products/{}/{}/tutorial/">Tutorial</a><a href="products/{}/{}/manual/{}/">Manual</a><a href="products/{}/{}/reference/">Reference</a><a class="portal-card-cite" href="citations/#{}">Cite</a></nav></article>"#,
                    index + 1,
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&product.title),
                    escape_html(&product.tagline),
                    escape_html(&product.title),
                    packages,
                    escape_html(&product.title),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(manual_route),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&product.id),
                )
            })
            .collect::<String>();
        let favicon = favicon_links("assets/");
        fs::write(
            output.join("index.html"),
            format!(
                r##"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="{}"><meta name="theme-color" content="#f9f6f0">{favicon}<title>αLoop · Research software for collider physics</title><link rel="stylesheet" href="assets/site.css"><script defer src="assets/site.js"></script></head><body class="portal-body"><a class="skip-link" href="#main-content">Skip to content</a><header class="portal-header"><a class="portal-brand" href="#overview" aria-label="αLoop home"><span class="portal-brand-logo" aria-hidden="true"></span><span class="portal-brand-copy"><strong>αLoop</strong><small>Local Unitarity research</small></span></a><nav class="portal-nav" aria-label="Primary"><a href="#projects">Projects</a><a href="people/">People</a><a href="publications/">Publications</a><a href="developers/">Developers</a></nav><div class="portal-header-actions"><a class="portal-source-link" href="https://github.com/alphal00p/gammaloop">GitHub <span aria-hidden="true">↗</span></a><button class="portal-theme-button" type="button" data-theme-toggle aria-label="Toggle color theme"><span aria-hidden="true">◐</span></button></div></header><main class="portal-main" id="main-content"><section class="portal-hero portal-section" id="overview"><div class="portal-hero-copy"><p class="portal-kicker">{}</p><h1>{}</h1><p class="portal-lede">{}</p><div class="portal-hero-actions"><a class="portal-button portal-button-primary" href="#projects">Explore the projects <span aria-hidden="true">↓</span></a><a class="portal-button" href="products/gammaloop/{}/tutorial/">Start with GammaLoop <span aria-hidden="true">↗</span></a><a class="portal-button" href="citations/">Cite the software <span aria-hidden="true">↗</span></a></div><dl class="portal-facts"><div><dt>5</dt><dd>connected projects</dd></div><div><dt>2</dt><dd>language ecosystems</dd></div><div><dt>∞</dt><dd>open development</dd></div></dl></div><div class="portal-hero-art" aria-label="αLoop collaboration mark" role="img"><div class="portal-wordmark"></div><p>Local cancellation.<br>Global precision.</p></div></section><section class="portal-section portal-projects" id="projects" aria-labelledby="projects-title"><div class="portal-section-heading"><div><p class="portal-kicker">Research software · 01—05</p><h2 id="projects-title">Projects &amp; crates</h2></div><p>Five connected codebases spanning numerical cross-sections, graph algorithms, tensor networks, symbolic identities, and integral evaluation.</p></div><div class="portal-project-grid">{projects}</div></section></main><footer class="portal-footer"><div><span class="portal-footer-mark" aria-hidden="true"></span><p><strong>αLoop</strong><br>Local Unitarity research software</p></div><nav aria-label="Footer"><a href="#projects">Projects</a><a href="people/">People</a><a href="publications/">Publications</a><a href="citations/">Cite</a><a href="developers/">Developers</a><a href="https://github.com/alphal00p/gammaloop">Source</a></nav><p>Physics, algorithms, and software<br>developed in the open.</p></footer></body></html>"##,
                escape_html(&self.portal.summary),
                escape_html(&self.portal.eyebrow),
                escape_html(&self.portal.title),
                escape_html(&self.portal.summary),
                escape_html(&channel_route),
            ),
        )?;
        self.write_people_page(output)?;
        self.write_publications_page(output)?;
        self.write_citations_page(output)?;
        fs::write(output.join(".nojekyll"), b"")?;
        Ok(())
    }

    fn write_people_page(&self, output: &Path) -> Result<()> {
        let cards = self
            .portal
            .people
            .iter()
            .map(|person| {
                let portrait = person.portrait.as_ref().map_or_else(
                    || {
                        format!(
                            "<span class=\"people-card-initials\" aria-hidden=\"true\">{}</span>",
                            escape_html(&person.initials),
                        )
                    },
                    |portrait| {
                        format!(
                            "<img class=\"people-card-portrait\" src=\"../assets/people/{}\" alt=\"\" width=\"720\" height=\"720\" loading=\"lazy\">",
                            escape_html(portrait),
                        )
                    },
                );
                let mut links = format!(
                    "<a href=\"{}\">Professional profile <span aria-hidden=\"true\">↗</span></a><a href=\"{}\">GitHub <span aria-hidden=\"true\">↗</span></a>",
                    escape_html(&person.url),
                    escape_html(&person.github),
                );
                if let Some(recid) = person.inspire_recid {
                    links.push_str(&format!(
                        "<a href=\"https://inspirehep.net/authors/{recid}\">INSPIRE HEP <span aria-hidden=\"true\">↗</span></a>"
                    ));
                }
                if let Some(orcid) = &person.orcid {
                    links.push_str(&format!(
                        "<a href=\"https://orcid.org/{}\">ORCID <span aria-hidden=\"true\">↗</span></a>",
                        escape_html(orcid),
                    ));
                }
                format!(
                    "<article class=\"people-card\" id=\"{}\">{portrait}<div class=\"people-card-copy\"><p class=\"portal-kicker\">Collaboration</p><h2>{}</h2><p>{}</p><nav aria-label=\"{} profiles\">{links}</nav></div></article>",
                    escape_html(&person.id),
                    escape_html(&person.name),
                    escape_html(&person.role),
                    escape_html(&person.name),
                )
            })
            .collect::<String>();
        let body = format!(
            "<header class=\"portal-page-hero\"><p class=\"portal-kicker\">People</p><h1>People building GammaLoop</h1><p>Researchers and collaborators developing GammaLoop, Local Unitarity methods, and the scientific software that supports them.</p></header><section class=\"people-page-grid\">{cards}</section>"
        );
        let directory = output.join("people");
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            portal_subpage_document(
                "People",
                "Researchers and collaborators building GammaLoop.",
                "people",
                &body,
            ),
        )?;
        Ok(())
    }

    fn write_publications_page(&self, output: &Path) -> Result<()> {
        let author_options = self
            .publications
            .authors
            .iter()
            .map(|author| {
                format!(
                    "<option value=\"{}\">{}</option>",
                    escape_html(&author.id),
                    escape_html(&author.name),
                )
            })
            .collect::<String>();
        let years = self
            .publications
            .publications
            .iter()
            .map(|publication| publication.year)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .rev()
            .map(|year| format!("<option value=\"{year}\">{year}</option>"))
            .collect::<String>();
        let types = self
            .publications
            .publications
            .iter()
            .flat_map(|publication| publication.types.iter())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .map(|kind| {
                format!(
                    "<option value=\"{}\">{}</option>",
                    escape_html(kind),
                    escape_html(&publication_type_label(kind)),
                )
            })
            .collect::<String>();
        let cards = self
            .publications
            .publications
            .iter()
            .map(|publication| {
                let people = publication.people.join(" ");
                let kinds = publication.types.join("|");
                let venue = publication.venue.as_ref().map_or_else(String::new, |venue| {
                    format!("<span>{}</span>", escape_html(venue))
                });
                let doi = publication.doi.as_ref().map_or_else(String::new, |doi| {
                    format!(
                        "<a href=\"https://doi.org/{}\">DOI <span aria-hidden=\"true\">↗</span></a>",
                        escape_html(doi),
                    )
                });
                let arxiv = publication.arxiv.as_ref().map_or_else(String::new, |arxiv| {
                    format!(
                        "<a href=\"https://arxiv.org/abs/{}\">arXiv <span aria-hidden=\"true\">↗</span></a>",
                        escape_html(arxiv),
                    )
                });
                format!(
                    "<article class=\"publication-card\" data-publication data-title=\"{}\" data-people=\"{}\" data-year=\"{}\" data-date=\"{}\" data-types=\"{}\" data-citations=\"{}\"><div class=\"publication-card-meta\"><time datetime=\"{}\">{}</time>{venue}<span>{} citations</span></div><h2><a href=\"{}\">{}</a></h2><p class=\"publication-authors\">{}</p><nav aria-label=\"Citation links\"><a href=\"{}\">INSPIRE</a>{doi}{arxiv}<a href=\"{}\">BibTeX</a></nav></article>",
                    escape_html(&publication.title.to_lowercase()),
                    escape_html(&people),
                    publication.year,
                    escape_html(&publication.date),
                    escape_html(&kinds),
                    publication.citations,
                    escape_html(&publication.date),
                    publication.year,
                    publication.citations,
                    escape_html(&publication.url),
                    escape_html(&publication.title),
                    escape_html(&compact_authors(&publication.authors)),
                    escape_html(&publication.url),
                    escape_html(&publication.bibtex_url),
                )
            })
            .collect::<String>();
        let count = self.publications.publications.len();
        let body = format!(
            "<header class=\"portal-page-hero\"><p class=\"portal-kicker\">Research output</p><h1>Publications</h1><p>Automatically assembled from stable INSPIRE HEP author identifiers, deduplicated across coauthors, and cached for a reproducible documentation build.</p><p class=\"publication-provenance\">Updated {} · <a href=\"{}\">Open the live INSPIRE query</a></p></header><form class=\"publication-filters\" data-publication-filters><label>Search<input type=\"search\" data-publication-search placeholder=\"Title or author\"></label><label>Author<select data-publication-author><option value=\"\">All authors</option>{author_options}</select></label><label>Year<select data-publication-year><option value=\"\">All years</option>{years}</select></label><label>Type<select data-publication-type><option value=\"\">All types</option>{types}</select></label><label>Sort<select data-publication-sort><option value=\"newest\">Newest</option><option value=\"cited\">Most cited</option></select></label><output data-publication-count aria-live=\"polite\">{count} publications</output></form><section class=\"publication-list\" data-publication-list>{cards}</section><noscript><p class=\"portal-page-note\">All records are shown above. Enable JavaScript only if you want to filter or sort them in the browser.</p></noscript>",
            escape_html(&self.publications.updated),
            escape_html(&self.publications.api_url),
        );
        let directory = output.join("publications");
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            portal_subpage_document(
                "Publications",
                "Filterable publications by verified αLoop contributors, sourced from INSPIRE HEP.",
                "publications",
                &body,
            ),
        )?;
        Ok(())
    }

    fn write_citations_page(&self, output: &Path) -> Result<()> {
        let cards = self
            .registry
            .product
            .iter()
            .map(|product| {
                let version = self.component_version(
                    product
                        .rust_components
                        .first()
                        .expect("validated product has a Rust component"),
                )?;
                let citation = software_citation(product, &version);
                let bibtex = software_bibtex(product, &version);
                let persistent = product.citation.doi.as_ref().map_or_else(
                    || {
                        "<p class=\"citation-status\">No registered software DOI is currently configured; this citation uses the versioned source repository.</p>".to_owned()
                    },
                    |doi| {
                        format!(
                            "<p class=\"citation-status\"><a href=\"https://doi.org/{}\">doi:{}</a></p>",
                            escape_html(doi),
                            escape_html(doi),
                        )
                    },
                );
                Ok::<_, eyre::Report>(format!(
                    "<article class=\"citation-card\" id=\"{}\"><p class=\"portal-kicker\">Version {}</p><h2>{}</h2>{persistent}<h3>Suggested citation</h3><pre id=\"citation-{}\"><code>{}</code></pre><button class=\"portal-button\" type=\"button\" data-copy-target=\"citation-{}\">Copy citation</button><details><summary>BibTeX</summary><pre id=\"bibtex-{}\"><code>{}</code></pre><button class=\"portal-button\" type=\"button\" data-copy-target=\"bibtex-{}\">Copy BibTeX</button></details></article>",
                    escape_html(&product.id),
                    escape_html(&version),
                    escape_html(&product.title),
                    escape_html(&product.id),
                    escape_html(&citation),
                    escape_html(&product.id),
                    escape_html(&product.id),
                    escape_html(&bibtex),
                    escape_html(&product.id),
                ))
            })
            .collect::<Result<String>>()?;
        let body = format!(
            "<header class=\"portal-page-hero\"><p class=\"portal-kicker\">Credit the software</p><h1>Cite αLoop projects</h1><p>Use the version-specific records below. Where a Zenodo DOI exists, it is the persistent citation target; otherwise the citation names the versioned source repository without inventing an identifier.</p></header><section class=\"citation-grid\">{cards}</section>"
        );
        let directory = output.join("citations");
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            portal_subpage_document(
                "Citations",
                "Version-specific citations for αLoop research software.",
                "",
                &body,
            ),
        )?;
        Ok(())
    }

    fn write_product_redirect(&self, product: &ProductConfig, output: &Path) -> Result<()> {
        let directory = output.join("products").join(&product.id);
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            format!(
                "<!doctype html><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"0; url=latest/\"><link rel=\"canonical\" href=\"latest/\"><title>{}</title><a href=\"latest/\">Open the latest {} documentation</a>",
                escape_html(&product.title),
                escape_html(&product.title)
            ),
        )?;
        Ok(())
    }

    fn validate_generated_links(&self, output: &Path, include_rustdoc: bool) -> Result<()> {
        let href = Regex::new(r#"(?:href|src)=["']([^"']+)["']"#)?;
        let roots = [output.to_path_buf()];
        let mut failures = vec![];
        for root in roots {
            for entry in WalkDir::new(&root) {
                let entry = entry?;
                if !entry.file_type().is_file()
                    || entry.path().extension().and_then(|value| value.to_str()) != Some("html")
                    || entry
                        .path()
                        .components()
                        .collect::<Vec<_>>()
                        .windows(2)
                        .any(|parts| {
                            parts[0].as_os_str() == "reference" && parts[1].as_os_str() == "rust"
                        })
                {
                    continue;
                }
                let html = fs::read_to_string(entry.path())?;
                for capture in href.captures_iter(&html) {
                    let target = capture.get(1).expect("capture exists").as_str();
                    if !include_rustdoc
                        && target
                            .split(['?', '#'])
                            .next()
                            .is_some_and(|path| path.contains("reference/rust/"))
                    {
                        continue;
                    }
                    validate_local_target(output, entry.path(), target, &mut failures)?;
                }
            }
            for entry in WalkDir::new(&root) {
                let entry = entry?;
                if !entry.file_type().is_file() || entry.file_name() != "search-index.json" {
                    continue;
                }
                let entries = serde_json::from_slice::<Vec<SearchEntry>>(&fs::read(entry.path())?)?;
                let source = entry
                    .path()
                    .parent()
                    .context("search index has no site directory")?
                    .join("index.html");
                for search in entries {
                    if !include_rustdoc && search.kind == "rust-api" {
                        continue;
                    }
                    validate_local_target(output, &source, &search.href, &mut failures)?;
                }
            }
        }
        ensure!(
            failures.is_empty(),
            "generated site has broken links:\n{}",
            failures.join("\n")
        );
        Ok(())
    }

    fn component_version(&self, component: &ComponentConfig) -> Result<String> {
        let path = self.root.join(&component.version_source);
        let source = fs::read_to_string(&path)
            .wrap_err_with(|| format!("failed to read {}", path.display()))?;
        let value = toml::from_str::<toml::Value>(&source)
            .wrap_err_with(|| format!("failed to parse {}", path.display()))?;
        let key = if path.file_name().and_then(|name| name.to_str()) == Some("pyproject.toml") {
            ["project", "version"]
        } else {
            ["package", "version"]
        };
        value
            .get(key[0])
            .and_then(|value| value.get(key[1]))
            .and_then(toml::Value::as_str)
            .map(ToOwned::to_owned)
            .wrap_err_with(|| format!("{} has no {}.{}", path.display(), key[0], key[1]))
    }

    fn component_default_features(&self, component: &ComponentConfig) -> Result<Vec<String>> {
        let path = self.root.join(&component.version_source);
        if path.file_name().and_then(|name| name.to_str()) != Some("Cargo.toml") {
            return Ok(Vec::new());
        }
        let source = fs::read_to_string(&path)
            .wrap_err_with(|| format!("failed to read {}", path.display()))?;
        let value = toml::from_str::<toml::Value>(&source)
            .wrap_err_with(|| format!("failed to parse {}", path.display()))?;
        Ok(value
            .get("features")
            .and_then(|features| features.get("default"))
            .and_then(toml::Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(toml::Value::as_str)
            .map(ToOwned::to_owned)
            .collect())
    }

    fn validate_snapshot_tag(&self, tag: &str) -> Result<()> {
        validate_route_segment(tag)?;
        ensure!(!tag.contains(".."), "snapshot tag cannot contain '..'");
        let (prefix, version) = tag
            .rsplit_once("-v")
            .map(|(prefix, version)| (Some(prefix), version))
            .or_else(|| tag.strip_prefix('v').map(|version| (None, version)))
            .wrap_err_with(|| format!("unrecognized snapshot tag {tag}"))?;
        ensure!(
            looks_like_version(version),
            "invalid release version in {tag}"
        );

        let component = if let Some(prefix) = prefix {
            self.registry
                .product
                .iter()
                .flat_map(|product| {
                    product
                        .rust_components
                        .iter()
                        .chain(&product.python_components)
                })
                .find(|component| component.package == prefix || component.id == prefix)
                .wrap_err_with(|| format!("tag {tag} does not name a registered component"))?
        } else {
            self.registry
                .product
                .iter()
                .find(|product| product.id == "gammaloop")
                .and_then(|product| {
                    product
                        .python_components
                        .iter()
                        .find(|component| component.package == "gammaloop")
                })
                .context("legacy GammaLoop tag has no registered distribution component")?
        };
        {
            let actual = self.component_version(component)?;
            ensure!(
                version == actual,
                "tag {tag} does not match {} version {actual}",
                component.package
            );
        }
        Ok(())
    }

    fn require_file(&self, relative: &Path) -> Result<()> {
        ensure_relative_path(relative)?;
        let path = self.root.join(relative);
        ensure!(
            path.is_file(),
            "required documentation source is missing: {}",
            path.display()
        );
        Ok(())
    }

    fn git_commit(&self) -> String {
        if let Some(commit) = ["ALPHAL00P_DOCS_GIT_COMMIT", "GITHUB_SHA"]
            .iter()
            .find_map(|name| env::var(name).ok())
            .filter(|commit| !commit.trim().is_empty())
        {
            return commit.trim().to_owned();
        }
        Command::new("git")
            .current_dir(&self.root)
            .args(["rev-parse", "HEAD"])
            .output()
            .ok()
            .filter(|output| output.status.success())
            .and_then(|output| String::from_utf8(output.stdout).ok())
            .map(|commit| commit.trim().to_owned())
            .unwrap_or_else(|| "unknown".to_owned())
    }

    fn git_timestamp(&self) -> u64 {
        if let Some(timestamp) = ["ALPHAL00P_DOCS_GIT_TIMESTAMP", "SOURCE_DATE_EPOCH"]
            .iter()
            .find_map(|name| env::var(name).ok())
            .and_then(|timestamp| timestamp.trim().parse().ok())
        {
            return timestamp;
        }
        Command::new("git")
            .current_dir(&self.root)
            .args(["log", "-1", "--format=%ct"])
            .output()
            .ok()
            .filter(|output| output.status.success())
            .and_then(|output| String::from_utf8(output.stdout).ok())
            .and_then(|timestamp| timestamp.trim().parse().ok())
            .unwrap_or(0)
    }
}

fn favicon_links(asset_root: &str) -> String {
    let asset_root = escape_html(asset_root);
    format!(
        "<link rel=\"icon\" type=\"image/svg+xml\" href=\"{asset_root}local-unitarity-light.svg\"><link rel=\"icon\" type=\"image/svg+xml\" href=\"{asset_root}local-unitarity-dark.svg\" media=\"(prefers-color-scheme: dark)\">"
    )
}

fn portal_subpage_document(title: &str, description: &str, active: &str, body: &str) -> String {
    let current = |item: &str| {
        if item == active {
            " aria-current=\"page\""
        } else {
            ""
        }
    };
    let favicon = favicon_links("../assets/");
    format!(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\"><meta name=\"theme-color\" content=\"#f9f6f0\">{favicon}<title>{} · αLoop</title><link rel=\"stylesheet\" href=\"../assets/site.css\"><script defer src=\"../assets/site.js\"></script></head><body class=\"portal-body portal-subpage-body\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"portal-header\"><a class=\"portal-brand\" href=\"../\" aria-label=\"αLoop home\"><span class=\"portal-brand-logo\" aria-hidden=\"true\"></span><span class=\"portal-brand-copy\"><strong>αLoop</strong><small>Local Unitarity research</small></span></a><nav class=\"portal-nav\" aria-label=\"Primary\"><a href=\"../#projects\">Projects</a><a href=\"../people/\"{}>People</a><a href=\"../publications/\"{}>Publications</a><a href=\"../developers/\">Developers</a></nav><div class=\"portal-header-actions\"><a class=\"portal-source-link\" href=\"https://github.com/alphal00p/gammaloop\">GitHub <span aria-hidden=\"true\">↗</span></a><button class=\"portal-theme-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\"><span aria-hidden=\"true\">◐</span></button></div></header><main class=\"portal-main portal-subpage-main\" id=\"main-content\">{body}</main><footer class=\"portal-footer\"><div><span class=\"portal-footer-mark\" aria-hidden=\"true\"></span><p><strong>αLoop</strong><br>Local Unitarity research software</p></div><nav aria-label=\"Footer\"><a href=\"../#projects\">Projects</a><a href=\"../people/\">People</a><a href=\"../publications/\">Publications</a><a href=\"../citations/\">Cite</a><a href=\"../developers/\">Developers</a><a href=\"https://github.com/alphal00p/gammaloop\">Source</a></nav><p>Physics, algorithms, and software<br>developed in the open.</p></footer></body></html>",
        escape_html(description),
        escape_html(title),
        current("people"),
        current("publications"),
    )
}

fn compact_authors(authors: &[String]) -> String {
    if authors.len() <= 8 {
        authors.join(", ")
    } else {
        format!("{}, et al.", authors[..6].join(", "))
    }
}

fn publication_type_label(kind: &str) -> String {
    kind.split(['-', '_'])
        .filter(|word| !word.is_empty())
        .map(|word| {
            let mut characters = word.chars();
            characters.next().map_or_else(String::new, |first| {
                first.to_uppercase().collect::<String>() + characters.as_str()
            })
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn software_citation(product: &ProductConfig, version: &str) -> String {
    let target = product.citation.doi.as_ref().map_or_else(
        || product.citation.repository.clone(),
        |doi| format!("https://doi.org/{doi}"),
    );
    format!(
        "{} ({}). {} (Version {}) [Computer software]. {}",
        product.citation.creators.join(", "),
        product.citation.year,
        product.citation.title,
        version,
        target,
    )
}

fn software_bibtex(product: &ProductConfig, version: &str) -> String {
    let doi = product
        .citation
        .doi
        .as_ref()
        .map_or_else(String::new, |doi| format!("  doi = {{{doi}}},\n"));
    let url = product.citation.doi.as_ref().map_or_else(
        || product.citation.repository.clone(),
        |doi| format!("https://doi.org/{doi}"),
    );
    format!(
        "@software{{{}{},\n  author = {{{}}},\n  title = {{{}}},\n  version = {{{version}}},\n  year = {{{}}},\n{doi}  url = {{{url}}}\n}}",
        product.id,
        product.citation.year,
        product.citation.creators.join(" and "),
        product.citation.title,
        product.citation.year,
    )
}

fn strip_first_markdown_title(source: &str) -> Result<&str> {
    let source = source.strip_prefix('\u{feff}').unwrap_or(source);
    let title_end = source
        .find('\n')
        .map(|index| index + 1)
        .unwrap_or(source.len());
    let title = source[..title_end].trim();
    ensure!(
        title.starts_with("# ") && !title.starts_with("## "),
        "developer Markdown note must begin with one H1 title"
    );
    Ok(&source[title_end..])
}

fn render_developer_markdown(source: &str) -> String {
    let mut options = Options::empty();
    options.insert(Options::ENABLE_TABLES);
    options.insert(Options::ENABLE_FOOTNOTES);
    options.insert(Options::ENABLE_STRIKETHROUGH);
    options.insert(Options::ENABLE_TASKLISTS);
    options.insert(Options::ENABLE_SMART_PUNCTUATION);
    let events = Parser::new_ext(source, options).map(|event| match event {
        Event::Html(value) | Event::InlineHtml(value) => Event::Text(value),
        event => event,
    });
    let mut rendered = String::new();
    html::push_html(&mut rendered, events);
    rendered
}

fn rewrite_developer_source_links(
    html: &str,
    source: &Path,
    commit: &str,
    note_routes: &BTreeMap<PathBuf, String>,
    repository_root: &Path,
) -> Result<String> {
    let attribute = Regex::new(r#"(?P<name>href|src)=\"(?P<target>[^\"]+)\""#)?;
    Ok(attribute
        .replace_all(html, |captures: &regex::Captures<'_>| {
            let name = captures.name("name").expect("attribute name").as_str();
            let target = captures.name("target").expect("attribute target").as_str();
            let external = target.is_empty()
                || target.starts_with('#')
                || target.starts_with('/')
                || target.starts_with("//")
                || target
                    .split(':')
                    .next()
                    .is_some_and(|scheme| matches!(scheme, "http" | "https" | "mailto" | "data"));
            if external {
                return format!("{name}=\"{target}\"");
            }

            let suffix_start = target.find(['?', '#']).unwrap_or(target.len());
            let (path, suffix) = target.split_at(suffix_start);
            let Some(parent) = source.parent() else {
                return captures[0].to_owned();
            };
            let Some(resolved) = normalize_repository_path(&parent.join(path)) else {
                return captures[0].to_owned();
            };
            let rewritten = if let Some(route) = note_routes.get(&resolved) {
                format!("../{route}/{suffix}")
            } else if repository_root.join(&resolved).exists() {
                let kind = if repository_root.join(&resolved).is_dir() {
                    "tree"
                } else {
                    "blob"
                };
                if name == "src" && kind == "blob" {
                    format!(
                        "https://raw.githubusercontent.com/alphal00p/gammaloop/{commit}/{}{suffix}",
                        resolved.to_string_lossy()
                    )
                } else {
                    format!(
                        "https://github.com/alphal00p/gammaloop/{kind}/{commit}/{}{suffix}",
                        resolved.to_string_lossy()
                    )
                }
            } else {
                target.to_owned()
            };
            format!("{name}=\"{}\"", escape_html(&rewritten))
        })
        .into_owned())
}

fn normalize_repository_path(path: &Path) -> Option<PathBuf> {
    let mut normalized = PathBuf::new();
    for component in path.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                if !normalized.pop() {
                    return None;
                }
            }
            Component::Normal(part) => normalized.push(part),
            Component::Prefix(_) | Component::RootDir => return None,
        }
    }
    Some(normalized)
}

fn validate_developer_html_note(source: &str) -> Result<()> {
    let active_element = Regex::new(r"(?is)<\s*(script|iframe|object|embed|form|base|link)\b")?;
    let event_handler = Regex::new(r"(?is)[\s/]on[a-z0-9_-]+\s*=")?;
    let active_meta = Regex::new(r"(?is)<\s*meta\b[^>]*\bhttp-equiv\s*=")?;
    let active_css = Regex::new(r"(?is)(?:@import|url\s*\()")?;
    let url_assignment = Regex::new(r"(?is)\b(?:href|src)\s*=")?;
    let quoted_url = Regex::new(r#"(?is)\b(?:href|src)\s*=\s*\"([^\"]*)\""#)?;
    let character_reference = Regex::new(r"(?is)&(?:#[x]?[0-9a-f]+|[a-z][a-z0-9]+);")?;
    ensure!(
        !active_element.is_match(source)
            && !event_handler.is_match(source)
            && !active_meta.is_match(source)
            && !active_css.is_match(source),
        "developer HTML notes cannot contain scripts, active embeds, forms, external styles, event handlers, active CSS, or HTTP-equivalent metadata"
    );
    ensure!(
        url_assignment.find_iter(source).count() == quoted_url.find_iter(source).count(),
        "developer HTML note href and src attributes must use double-quoted values"
    );
    for capture in quoted_url.captures_iter(source) {
        let target = capture.get(1).expect("URL attribute capture").as_str();
        ensure!(
            !character_reference.is_match(target) && !target.chars().any(char::is_control),
            "developer HTML note URLs cannot contain character references or control characters"
        );
        let browser_normalized_target =
            target.trim_start_matches(|character: char| character <= ' ');
        if let Some((scheme, _)) = browser_normalized_target.split_once(':')
            && !scheme.is_empty()
            && scheme.chars().all(|character| {
                character.is_ascii_alphanumeric() || matches!(character, '+' | '-' | '.')
            })
        {
            ensure!(
                matches!(
                    scheme.to_ascii_lowercase().as_str(),
                    "http" | "https" | "mailto"
                ),
                "developer HTML note uses a forbidden URL scheme: {scheme}"
            );
        }
    }
    Ok(())
}

fn decorate_developer_diagram(source: &str, title: &str, source_url: &str) -> Result<String> {
    ensure!(
        source.contains("</head>"),
        "developer HTML note has no head"
    );
    let banner_style = r#"<meta http-equiv="Content-Security-Policy" content="default-src 'none'; style-src 'unsafe-inline'; img-src data:; font-src data:"><style>.developer-diagram-banner{position:sticky;top:0;z-index:1000;display:flex;flex-wrap:wrap;align-items:center;gap:.65rem 1rem;padding:.7rem 1rem;border-bottom:1px solid #b893c7;background:#3d2645;color:#f9f6f0;font:600 12px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;letter-spacing:.04em}.developer-diagram-banner strong{text-transform:uppercase}.developer-diagram-banner a{color:#eadff0}.developer-diagram-banner a:first-of-type{margin-left:auto}@media(max-width:38rem){.developer-diagram-banner a:first-of-type{margin-left:0}}</style></head>"#;
    let mut rendered = source.replacen("</head>", banner_style, 1);
    let body = rendered
        .find("<body")
        .context("developer HTML note has no body")?;
    let body_start = rendered[body..]
        .find('>')
        .map(|offset| body + offset + 1)
        .context("developer HTML note has an unterminated body")?;
    let banner = format!(
        "<aside class=\"developer-diagram-banner\" aria-label=\"Developer note\"><strong>For developers</strong><span>{}</span><a href=\"../../\">All architecture notes</a><a href=\"{}\">Source ↗</a></aside>",
        escape_html(title),
        escape_html(source_url),
    );
    rendered.insert_str(body_start, &banner);
    Ok(rendered)
}

fn developer_toc(headings: &[HeadingLink]) -> String {
    let links = headings
        .iter()
        .filter(|heading| matches!(heading.level, 2 | 3))
        .map(|heading| {
            format!(
                "<a class=\"toc-link\" data-level=\"{}\" href=\"#{}\">{}</a>",
                heading.level,
                escape_html(&heading.id),
                escape_html(&heading.title),
            )
        })
        .collect::<String>();
    if links.is_empty() {
        String::new()
    } else {
        format!(
            "<aside class=\"docs-toc\" aria-label=\"On this page\"><p class=\"toc-title\">On this page</p>{links}</aside>"
        )
    }
}

fn developer_inline_toc(headings: &[HeadingLink]) -> String {
    let links = headings
        .iter()
        .filter(|heading| matches!(heading.level, 2 | 3))
        .map(|heading| {
            format!(
                "<a data-level=\"{}\" href=\"#{}\">{}</a>",
                heading.level,
                escape_html(&heading.id),
                escape_html(&heading.title),
            )
        })
        .collect::<String>();
    if links.is_empty() {
        String::new()
    } else {
        format!(
            "<details class=\"developer-inline-toc\"><summary>On this page</summary><nav aria-label=\"On this page\">{links}</nav></details>"
        )
    }
}

fn developer_document(
    title: &str,
    description: &str,
    docs_root: &str,
    portal_root: &str,
    sidebar: &str,
    main: &str,
    toc: &str,
) -> String {
    let home = if docs_root.is_empty() {
        "./"
    } else {
        docs_root
    };
    let favicon = favicon_links(&format!("{docs_root}assets/"));
    format!(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\"><meta name=\"theme-color\" content=\"#f9f6f0\">{favicon}<title>{}</title><link rel=\"stylesheet\" href=\"{}assets/site.css\"><script defer src=\"{}assets/site.js\"></script></head><body class=\"developer-body\" data-docs-root=\"{}\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"site-header\"><button class=\"header-button menu-button\" type=\"button\" data-menu-toggle aria-label=\"Open navigation\" aria-controls=\"docs-sidebar\" aria-expanded=\"false\">☰</button><a class=\"site-brand\" href=\"{}\"><span class=\"site-brand-mark\">α</span><span class=\"site-brand-name\">Developer notes</span></a><div class=\"site-header-tools\"><a class=\"header-button header-link\" href=\"{}\" aria-label=\"Research documentation\"><span class=\"header-link-label\">Research docs</span><span aria-hidden=\"true\">←</span></a><button class=\"header-button\" type=\"button\" data-search-open aria-label=\"Search developer documentation\"><span class=\"header-search-label\">Search</span> <span class=\"header-button-label\">⌘K</span><span class=\"header-search-icon\" aria-hidden=\"true\">⌕</span></button><button class=\"header-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\">◐</button></div></header><div class=\"docs-shell\">{sidebar}<main class=\"docs-main developer-main\" id=\"main-content\">{main}</main>{toc}</div><button class=\"sidebar-backdrop\" type=\"button\" data-sidebar-backdrop aria-label=\"Close navigation\"></button><dialog class=\"search-dialog\" data-search-dialog><form class=\"search-form\" method=\"dialog\"><input class=\"search-input\" type=\"search\" data-search-input placeholder=\"Search developer architecture\" aria-label=\"Search developer architecture\"><button class=\"header-button\" value=\"close\">Close</button></form><ul class=\"search-results\" data-search-results aria-live=\"polite\"></ul></dialog></body></html>",
        escape_html(description),
        escape_html(title),
        escape_html(docs_root),
        escape_html(docs_root),
        escape_html(docs_root),
        escape_html(home),
        escape_html(portal_root),
    )
}

fn extract_html_body(html: &str) -> Result<&str> {
    let body = html.find("<body").context("rendered page has no body")?;
    let start = html[body..]
        .find('>')
        .map(|offset| body + offset + 1)
        .context("rendered page has an unterminated body tag")?;
    let end = html
        .rfind("</body>")
        .context("rendered page has no body end")?;
    ensure!(start <= end, "rendered page body is malformed");
    Ok(&html[start..end])
}

fn page_root_prefix(route: &str) -> String {
    "../".repeat(
        route
            .split('/')
            .filter(|segment| !segment.is_empty())
            .count(),
    )
}

fn product_hero(product: &ProductConfig) -> String {
    let first_manual = product
        .pages
        .iter()
        .find(|page| page.group == "Manual")
        .map(|page| page.route.as_str())
        .unwrap_or("manual/");
    let logo = product.logo.as_deref().map_or_else(String::new, |logo| {
        format!(
            "<span class=\"product-hero-logo product-logo-{}\" role=\"img\" aria-label=\"{} logo\"></span>",
            escape_html(logo),
            escape_html(&product.title),
        )
    });
    format!(
        "<header class=\"product-hero\">{logo}<p class=\"product-eyebrow\">Research software documentation</p><h1>{}</h1><p>{}</p><div class=\"hero-actions\"><a class=\"hero-action primary\" href=\"tutorial/\">Start the tutorial</a><a class=\"hero-action\" href=\"{}\">Read the manual</a><a class=\"hero-action\" href=\"reference/\">Browse reference</a></div></header>",
        escape_html(&product.title),
        escape_html(&product.tagline),
        escape_html(first_manual),
    )
}

fn promote_heading_hierarchy(body: &mut String) {
    if body.contains("<h1") {
        return;
    }
    for level in 2..=6 {
        *body = body
            .replace(&format!("<h{level}"), &format!("<h{}", level - 1))
            .replace(&format!("</h{level}>"), &format!("</h{}>", level - 1));
    }
}

fn rewrite_page_links(body: &str, docs_root: &str) -> Result<String> {
    if docs_root.is_empty() {
        return Ok(body.to_owned());
    }
    let attribute = Regex::new(r#"(?P<name>href|src)=\"(?P<target>[^\"]+)\""#)?;
    Ok(attribute
        .replace_all(body, |captures: &regex::Captures<'_>| {
            let name = captures.name("name").expect("attribute name").as_str();
            let target = captures.name("target").expect("attribute target").as_str();
            let external = target.starts_with('#')
                || target.starts_with('/')
                || target.starts_with("//")
                || target
                    .split(':')
                    .next()
                    .is_some_and(|scheme| matches!(scheme, "http" | "https" | "mailto" | "data"));
            if external {
                format!("{name}=\"{target}\"")
            } else {
                format!("{name}=\"{docs_root}{target}\"")
            }
        })
        .into_owned())
}

fn inject_heading_ids(body: &str) -> Result<(String, Vec<HeadingLink>)> {
    let heading = Regex::new(r"(?s)<h(?P<level>[1-6])(?P<attrs>[^>]*)>(?P<body>.*?)</h[1-6]>")?;
    let id_attribute = Regex::new(r#"\sid=\"(?P<id>[^\"]+)\""#)?;
    let tags = Regex::new(r"<[^>]+>")?;
    let mut used = BTreeSet::new();
    let mut headings = Vec::new();
    let mut rendered = String::with_capacity(body.len() + 256);
    let mut cursor = 0;
    for captures in heading.captures_iter(body) {
        let whole = captures.get(0).expect("heading capture exists");
        let level = captures
            .name("level")
            .expect("heading level exists")
            .as_str()
            .parse::<u8>()?;
        let attrs = captures
            .name("attrs")
            .expect("heading attributes exist")
            .as_str();
        let heading_body = captures.name("body").expect("heading body exists").as_str();
        let title = decode_html_text(tags.replace_all(heading_body, "").trim());
        let mut id = id_attribute
            .captures(attrs)
            .and_then(|captures| captures.name("id").map(|id| id.as_str().to_owned()))
            .unwrap_or_else(|| slug(&title));
        if id.is_empty() {
            id = "section".to_owned();
        }
        if !used.insert(id.clone()) {
            let base = id.clone();
            let mut suffix = 2;
            while !used.insert(format!("{base}-{suffix}")) {
                suffix += 1;
            }
            id = format!("{base}-{suffix}");
        }
        rendered.push_str(&body[cursor..whole.start()]);
        if id_attribute.is_match(attrs) {
            rendered.push_str(whole.as_str());
        } else {
            rendered.push_str(&format!(
                "<h{level}{attrs} id=\"{}\">{heading_body}</h{level}>",
                escape_html(&id)
            ));
        }
        headings.push(HeadingLink { level, title, id });
        cursor = whole.end();
    }
    rendered.push_str(&body[cursor..]);
    Ok((rendered, headings))
}

fn render_page_navigation(
    previous: Option<&SitePage>,
    next: Option<&SitePage>,
    docs_root: &str,
) -> String {
    let previous = previous.map_or_else(
        || "<span></span>".to_owned(),
        |page| {
            format!(
                "<a class=\"page-nav-link\" href=\"{}{}\"><span class=\"page-nav-label\">← Previous</span>{}</a>",
                escape_html(docs_root),
                escape_html(&page.route),
                escape_html(&page.title)
            )
        },
    );
    let next = next.map_or_else(
        || "<span></span>".to_owned(),
        |page| {
            format!(
                "<a class=\"page-nav-link next\" href=\"{}{}\"><span class=\"page-nav-label\">Next →</span>{}</a>",
                escape_html(docs_root),
                escape_html(&page.route),
                escape_html(&page.title)
            )
        },
    );
    format!("<nav class=\"page-nav\" aria-label=\"Manual pagination\">{previous}{next}</nav>")
}

fn markdown_changelog_to_typst(markdown: &str) -> String {
    let mut output = String::new();
    let mut code_language = None::<String>;
    let mut code = String::new();
    for line in markdown.lines() {
        if let Some(language) = &code_language {
            if line.trim_start().starts_with("```") {
                output.push_str(&format!(
                    "#raw(\"{}\", block: true, lang: \"{}\")\n\n",
                    typst_string(code.trim_end()),
                    typst_string(language)
                ));
                code_language = None;
                code.clear();
            } else {
                code.push_str(line);
                code.push('\n');
            }
            continue;
        }
        if let Some(language) = line.trim_start().strip_prefix("```") {
            code_language = Some(language.trim().to_owned());
            continue;
        }
        if line.is_empty() {
            output.push('\n');
            continue;
        }
        if let Some((marks, title)) = markdown_heading(line) {
            // The imported changelog is nested below its own canonical heading.
            let level = (marks + 1).min(6);
            output.push_str(&format!(
                "{} {}\n",
                "=".repeat(level),
                markdown_inline_to_typst(title)
            ));
        } else if let Some(item) = line.strip_prefix("- ") {
            output.push_str(&format!("- {}\n", markdown_inline_to_typst(item)));
        } else if let Some(item) = markdown_ordered_item(line) {
            output.push_str(&format!("+ {}\n", markdown_inline_to_typst(item)));
        } else {
            output.push_str(&markdown_inline_to_typst(line));
            output.push('\n');
        }
    }
    if code_language.is_some() {
        output.push_str(&format!(
            "#raw(\"{}\", block: true)\n",
            typst_string(code.trim_end())
        ));
    }
    output
}

fn markdown_heading(line: &str) -> Option<(usize, &str)> {
    let marks = line
        .chars()
        .take_while(|character| *character == '#')
        .count();
    (marks > 0 && marks <= 6)
        .then(|| {
            line.get(marks..)?
                .strip_prefix(' ')
                .map(|title| (marks, title))
        })
        .flatten()
}

fn markdown_ordered_item(line: &str) -> Option<&str> {
    let (number, item) = line.split_once(". ")?;
    (!number.is_empty() && number.chars().all(|character| character.is_ascii_digit()))
        .then_some(item)
}

fn markdown_inline_to_typst(line: &str) -> String {
    let links = Regex::new(r"\[([^\]]+)\]\(([^)]+)\)").expect("static Markdown link regex");
    let mut rendered = String::new();
    let mut cursor = 0;
    for captures in links.captures_iter(line) {
        let whole = captures.get(0).expect("whole Markdown link");
        if whole.start() > cursor {
            rendered.push_str(&format!(
                "#raw(\"{}\")",
                typst_string(&line[cursor..whole.start()])
            ));
        }
        rendered.push_str(&format!(
            "#link(\"{}\")[#raw(\"{}\")]",
            typst_string(&captures[2]),
            typst_string(&captures[1])
        ));
        cursor = whole.end();
    }
    if cursor < line.len() {
        rendered.push_str(&format!("#raw(\"{}\")", typst_string(&line[cursor..])));
    }
    if rendered.is_empty() {
        rendered.push_str("#raw(\"\")");
    }
    rendered
}

fn count_catalog_items(scope: &DocScope) -> usize {
    scope.items.len()
        + scope
            .scopes
            .values()
            .map(count_catalog_items)
            .sum::<usize>()
}

fn render_python_catalog(
    catalog: &DocCatalog,
    stub_name: &str,
    stub_text: &str,
    git_commit: &str,
) -> String {
    let mut body = format!(
        "<p><code>{}</code> · {} documented exports · version {}</p>",
        escape_html(&catalog.component.package),
        count_catalog_items(&catalog.root),
        escape_html(&catalog.component.version)
    );
    if let Some(summary) = &catalog.root.summary {
        body.push_str(&format!("<p>{}</p>", escape_html(summary)));
    }
    if let Some(docs) = &catalog.root.docs {
        body.push_str(&render_doc_text(docs));
    }
    let mut scope_path = Vec::new();
    render_python_scope(
        catalog,
        &catalog.root,
        &mut scope_path,
        git_commit,
        2,
        &mut body,
    );
    body.push_str(&format!(
        "<details class=\"stub-source\"><summary>Type-checker stub source</summary><p><a href=\"reference/python/{}\">Download <code>{}</code></a>. This source is retained for type checkers; the structured reference above is the human-facing documentation.</p><pre><code>{}</code></pre></details>",
        escape_html(stub_name),
        escape_html(stub_name),
        escape_html(stub_text)
    ));
    body
}

fn render_python_scope(
    catalog: &DocCatalog,
    scope: &DocScope,
    path: &mut Vec<String>,
    git_commit: &str,
    level: usize,
    body: &mut String,
) {
    for item in scope.items.values() {
        let anchor = python_item_anchor(path, item);
        let heading = level.clamp(2, 6);
        body.push_str(&format!(
            "<section class=\"api-item\" id=\"{}\"><h{heading}><code>{}</code></h{heading}><span class=\"api-kind\">{}</span>{}",
            escape_html(&anchor),
            escape_html(&item.name),
            escape_html(&format!("{:?}", item.kind)),
            if item.supported {
                ""
            } else {
                "<span class=\"api-kind\">implementation detail</span>"
            }
        ));
        if let Some(docs) = &item.docs {
            body.push_str(&render_doc_text(docs));
        } else if let Some(summary) = &item.summary {
            body.push_str(&format!("<p>{}</p>", escape_html(summary)));
        }
        if let Some(signature) = &item.signature {
            body.push_str(&format!(
                "<pre class=\"api-signature\"><code>{}</code></pre>",
                escape_html(signature)
            ));
        }
        if !item.required_features.is_empty() {
            body.push_str("<p><strong>Requires:</strong> ");
            for feature in &item.required_features {
                body.push_str(&format!(
                    "<span class=\"api-feature\">{}</span>",
                    escape_html(feature)
                ));
            }
            body.push_str("</p>");
        }
        if !item.params.is_empty() {
            body.push_str("<h4>Parameters</h4><table><thead><tr><th>Name</th><th>Type</th><th>Default</th><th>Description</th></tr></thead><tbody>");
            for parameter in &item.params {
                let docs = parameter
                    .docs
                    .as_ref()
                    .map(|docs| escape_html(&docs.body))
                    .unwrap_or_default();
                body.push_str(&format!(
                    "<tr><td><code>{}</code>{}</td><td><code>{}</code></td><td><code>{}</code></td><td>{}</td></tr>",
                    escape_html(&parameter.name),
                    if parameter.required { " <span class=\"api-kind\">required</span>" } else { "" },
                    escape_html(parameter.ty.as_deref().unwrap_or("—")),
                    escape_html(parameter.default.as_deref().unwrap_or("—")),
                    docs,
                ));
            }
            body.push_str("</tbody></table>");
        }
        if let Some(returns) = &item.returns {
            body.push_str("<h4>Returns</h4>");
            body.push_str(&render_doc_text(returns));
        }
        for example in &item.examples {
            body.push_str(&format!(
                "<h4>{}</h4><pre><code data-lang=\"{}\">{}</code></pre>",
                escape_html(&example.title),
                escape_html(&example.language),
                escape_html(&example.code)
            ));
        }
        render_python_members(&item.members, &anchor, 4, body);
        if let Some(source) = &item.source {
            body.push_str(&format!(
                "<p><a href=\"https://github.com/alphal00p/gammaloop/blob/{}/{}#L{}\">Source: <code>{}:{}</code></a></p>",
                escape_html(git_commit),
                escape_html(&source.file),
                source.line,
                escape_html(&source.file),
                source.line,
            ));
        }
        body.push_str("</section>");
    }
    for child in scope.scopes.values() {
        path.push(child.id.clone());
        if child.id != "supported" {
            let heading = level.clamp(2, 6);
            body.push_str(&format!(
                "<h{heading}>{}</h{heading}>",
                escape_html(&child.title)
            ));
            if let Some(summary) = &child.summary {
                body.push_str(&format!("<p>{}</p>", escape_html(summary)));
            }
            if let Some(docs) = &child.docs {
                body.push_str(&render_doc_text(docs));
            }
        }
        render_python_scope(catalog, child, path, git_commit, level + 1, body);
        path.pop();
    }
    let _ = catalog;
}

fn render_python_members(
    members: &[DocMember],
    parent_anchor: &str,
    level: usize,
    body: &mut String,
) {
    let mut anchor_counts = BTreeMap::new();
    for member in members {
        let anchor = next_python_member_anchor(parent_anchor, member, &mut anchor_counts);
        let heading = level.clamp(3, 6);
        body.push_str(&format!(
            "<section class=\"api-member\" id=\"{}\"><h{heading}><code>{}</code></h{heading}><span class=\"api-kind\">{:?}</span>",
            escape_html(&anchor),
            escape_html(&member.name),
            member.kind,
        ));
        if let Some(signature) = &member.signature {
            body.push_str(&format!(
                "<pre class=\"api-member-signature\"><code>{}</code></pre>",
                escape_html(signature)
            ));
        }
        if let Some(docs) = &member.docs {
            body.push_str(&render_doc_text(docs));
        }
        if let Some(default) = &member.default {
            body.push_str(&format!(
                "<p><strong>Default:</strong> <code>{}</code></p>",
                escape_html(default)
            ));
        }
        render_python_members(&member.members, &anchor, level + 1, body);
        body.push_str("</section>");
    }
}

fn render_doc_text(docs: &alphal00p_docs_schema::DocText) -> String {
    let paragraphs = docs
        .body
        .split("\n\n")
        .filter(|paragraph| !paragraph.trim().is_empty())
        .map(|paragraph| {
            escape_html(paragraph.trim())
                .replace("\r\n", "<br>")
                .replace('\n', "<br>")
        })
        .map(|paragraph| format!("<p>{paragraph}</p>"))
        .collect::<String>();
    format!("<div class=\"api-docstring\">{paragraphs}</div>")
}

fn python_item_anchor(path: &[String], item: &DocItem) -> String {
    let mut parts = path.iter().map(String::as_str).collect::<Vec<_>>();
    parts.push(&item.id);
    slug(&parts.join("-"))
}

fn next_python_member_anchor(
    parent_anchor: &str,
    member: &DocMember,
    counts: &mut BTreeMap<String, usize>,
) -> String {
    let base = format!(
        "{parent_anchor}-{}-{}",
        slug(&member.name),
        slug(&format!("{:?}", member.kind))
    );
    let occurrence = counts.entry(base.clone()).or_default();
    *occurrence += 1;
    if *occurrence == 1 {
        base
    } else {
        format!("{base}-{occurrence}")
    }
}

fn append_catalog_search(catalog: &DocCatalog, scope: &DocScope, entries: &mut Vec<SearchEntry>) {
    let mut path = Vec::new();
    append_catalog_search_at(catalog, scope, &mut path, entries);
}

fn append_catalog_search_at(
    catalog: &DocCatalog,
    scope: &DocScope,
    path: &mut Vec<String>,
    entries: &mut Vec<SearchEntry>,
) {
    for item in scope.items.values() {
        let (href, kind, member_anchor) = match catalog.component.language {
            ApiLanguage::Rust => {
                let anchor = python_item_anchor(path, item);
                (
                    format!(
                        "reference/rust/supported/{}/#{anchor}",
                        catalog.component.id
                    ),
                    "rust-api",
                    Some(anchor),
                )
            }
            ApiLanguage::Python => {
                let anchor = python_item_anchor(path, item);
                (
                    format!("reference/python/{}/#{anchor}", catalog.component.id),
                    "python-api",
                    Some(anchor),
                )
            }
            _ => continue,
        };
        entries.push(SearchEntry {
            title: format!("{} · {}", catalog.component.id, item.title),
            summary: item.summary.clone().unwrap_or_default(),
            href: href.clone(),
            kind: kind.to_owned(),
        });
        append_member_search(
            &catalog.component.id,
            &item.title,
            &href,
            kind,
            member_anchor.as_deref(),
            &item.members,
            entries,
        );
    }
    for child in scope.scopes.values() {
        path.push(child.id.clone());
        append_catalog_search_at(catalog, child, path, entries);
        path.pop();
    }
}

fn append_member_search(
    component: &str,
    parent: &str,
    href: &str,
    kind: &str,
    parent_anchor: Option<&str>,
    members: &[DocMember],
    entries: &mut Vec<SearchEntry>,
) {
    let mut anchor_counts = BTreeMap::new();
    for member in members {
        let title = format!("{parent}.{}", member.name);
        let member_anchor = parent_anchor
            .map(|anchor| next_python_member_anchor(anchor, member, &mut anchor_counts));
        let member_href = member_anchor.as_ref().map_or_else(
            || href.to_owned(),
            |anchor| format!("{}#{anchor}", href.split('#').next().unwrap_or(href)),
        );
        let summary = member
            .docs
            .as_ref()
            .map(|docs| docs.body.split_whitespace().collect::<Vec<_>>().join(" "))
            .filter(|docs| !docs.is_empty())
            .or_else(|| member.signature.clone())
            .unwrap_or_default();
        entries.push(SearchEntry {
            title: format!("{component} · {title}"),
            summary,
            href: member_href.clone(),
            kind: kind.to_owned(),
        });
        append_member_search(
            component,
            &title,
            &member_href,
            kind,
            member_anchor.as_deref(),
            &member.members,
            entries,
        );
    }
}

fn append_scope_typst(
    catalog: &DocCatalog,
    scope: &DocScope,
    git_commit: &str,
    output: &mut String,
) {
    if let Some(docs) = &scope.docs {
        match docs.format {
            DocFormat::TypstMarkup => {
                output.push_str(&docs.body);
                output.push_str("\n\n");
            }
            DocFormat::RustMarkdown | DocFormat::PythonDocstring | DocFormat::PlainText => {
                output.push_str(&format!(
                    "#raw(\"{}\", block: true)\n\n",
                    typst_string(&docs.body)
                ));
            }
        }
    }
    for item in scope.items.values().filter(|item| item.supported) {
        output.push_str(&format!("\n=== #raw(\"{}\")\n", typst_string(&item.title)));
        if let Some(summary) = &item.summary {
            output.push_str(&format!("#raw(\"{}\")\n\n", typst_string(summary)));
        }
        if let Some(docs) = &item.docs {
            match docs.format {
                DocFormat::TypstMarkup => {
                    output.push_str(&docs.body);
                    output.push_str("\n\n");
                }
                DocFormat::RustMarkdown | DocFormat::PythonDocstring | DocFormat::PlainText => {
                    output.push_str(&format!(
                        "#raw(\"{}\", block: true)\n\n",
                        typst_string(&docs.body)
                    ))
                }
            }
        }
        if let Some(signature) = &item.signature {
            output.push_str(&format!(
                "#raw(\"{}\", block: true)\n\n",
                typst_string(signature)
            ));
        }
        if !item.required_features.is_empty() {
            output.push_str(&format!(
                "*Requires features:* #raw(\"{}\")\n\n",
                typst_string(&item.required_features.join(", "))
            ));
        } else if !catalog.component.features.is_empty() {
            output.push_str(&format!(
                "*Reference feature matrix:* #raw(\"{}\") (item is available without an additional gate)\n\n",
                typst_string(&catalog.component.features.join(", "))
            ));
        }
        if !item.params.is_empty() {
            let mut params = String::new();
            for param in &item.params {
                params.push_str(&format!(
                    "[#raw(\"{}\")], [#raw(\"{}\")],\n",
                    typst_string(&param.name),
                    typst_string(param.ty.as_deref().unwrap_or("unspecified")),
                ));
            }
            output.push_str(&format!(
                "#table(columns: (1fr, 2fr), table.header([*Parameter*], [*Type*]), {params})\n\n"
            ));
        }
        if let Some(returns) = &item.returns {
            output.push_str(&format!(
                "*Returns:* #raw(\"{}\")\n\n",
                typst_string(&returns.body)
            ));
        }
        for example in &item.examples {
            output.push_str(&format!(
                "*{}*\n#raw(\"{}\", block: true, lang: \"{}\")\n\n",
                typst_string(&example.title),
                typst_string(&example.code),
                typst_string(&example.language),
            ));
        }
        append_members_typst(&item.members, 4, output);
        let reference = match catalog.component.language {
            ApiLanguage::Rust => rustdoc_href(catalog, item),
            ApiLanguage::Python => format!("reference/python/#{}", catalog.component.id),
            _ => String::new(),
        };
        if !reference.is_empty() {
            output.push_str(&format!(
                "#link(\"{}\")[Exhaustive reference]",
                typst_string(&reference)
            ));
        }
        if let Some(source) = &item.source {
            let url = format!(
                "https://github.com/alphal00p/gammaloop/blob/{git_commit}/{}#L{}",
                source.file, source.line
            );
            output.push_str(&format!(" · #link(\"{}\")[Source]\n", typst_string(&url)));
        }
    }
    for child in scope.scopes.values() {
        append_scope_typst(catalog, child, git_commit, output);
    }
}

fn append_members_typst(members: &[DocMember], depth: usize, output: &mut String) {
    for member in members {
        let heading = "=".repeat(depth.min(6));
        output.push_str(&format!(
            "\n{heading} #raw(\"{}\")\n*Kind:* #raw(\"{}\")\n\n",
            typst_string(&member.name),
            typst_string(&format!("{:?}", member.kind)),
        ));
        if let Some(signature) = &member.signature {
            output.push_str(&format!(
                "#raw(\"{}\", block: true)\n\n",
                typst_string(signature)
            ));
        }
        if let Some(default) = &member.default {
            output.push_str(&format!(
                "*Default:* #raw(\"{}\")\n\n",
                typst_string(default)
            ));
        }
        if let Some(docs) = &member.docs {
            match docs.format {
                DocFormat::TypstMarkup => {
                    output.push_str(&docs.body);
                    output.push_str("\n\n");
                }
                DocFormat::RustMarkdown | DocFormat::PythonDocstring | DocFormat::PlainText => {
                    output.push_str(&format!(
                        "#raw(\"{}\", block: true)\n\n",
                        typst_string(&docs.body)
                    ));
                }
            }
        }
        append_members_typst(&member.members, depth + 1, output);
    }
}

fn rustdoc_href(catalog: &DocCatalog, item: &DocItem) -> String {
    let crate_name = catalog.component.package.replace('-', "_");
    let source = item
        .source
        .as_ref()
        .map(|source| source.identifier.as_str())
        .unwrap_or_default();
    let mut path = source.split("::").collect::<Vec<_>>();
    if path.first().is_some_and(|segment| *segment == crate_name) {
        path.remove(0);
    }
    let name = path.pop().unwrap_or(&item.name).trim_end_matches('!');
    let module = if path.len() > 1
        || path.first().is_some_and(|segment| {
            !matches!(item.kind, alphal00p_docs_schema::DocItemKind::Method) && *segment != name
        }) {
        path.join("/")
    } else {
        String::new()
    };
    let base = if module.is_empty() {
        format!("reference/rust/{crate_name}/")
    } else {
        format!("reference/rust/{crate_name}/{module}/")
    };
    use alphal00p_docs_schema::DocItemKind;
    match item.kind {
        DocItemKind::Function => format!("{base}fn.{name}.html"),
        DocItemKind::Trait => format!("{base}trait.{name}.html"),
        DocItemKind::ExportedMacro => format!(
            "{base}{}.{name}.html",
            rustdoc_macro_prefix(item.signature.as_deref())
        ),
        DocItemKind::Method => {
            let owner = path.pop().unwrap_or_default();
            let module = path.join("/");
            let base = if module.is_empty() {
                format!("reference/rust/{crate_name}/")
            } else {
                format!("reference/rust/{crate_name}/{module}/")
            };
            format!("{base}struct.{owner}.html#method.{name}")
        }
        DocItemKind::Type | DocItemKind::Setting => {
            format!("{base}{}.{name}.html", rustdoc_type_prefix(item))
        }
        _ => format!("reference/rust/{crate_name}/index.html"),
    }
}

fn rustdoc_macro_prefix(signature: Option<&str>) -> &'static str {
    if signature.is_some_and(|signature| signature.contains("#[derive(")) {
        "derive"
    } else {
        "macro"
    }
}

fn rustdoc_type_prefix(item: &DocItem) -> &'static str {
    let signature = item.signature.as_deref().unwrap_or_default();
    let tokens = signature
        .split(|character: char| !character.is_ascii_alphanumeric() && character != '_')
        .filter(|token| !token.is_empty())
        .collect::<Vec<_>>();

    // Source-backed signatures retain doc attributes, so looking for a bare
    // keyword misclassifies structs whenever their prose says "type". Match
    // the declaration keyword immediately preceding this item's name instead.
    tokens
        .windows(2)
        .rev()
        .find_map(|pair| {
            (pair[1] == item.name).then(|| match pair[0] {
                "type" => Some("type"),
                "enum" => Some("enum"),
                "union" => Some("union"),
                "struct" => Some("struct"),
                _ => None,
            })?
        })
        .unwrap_or("struct")
}

fn absolute_from(root: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        root.join(path)
    }
}

fn prepare_persistent_rustdoc_target(root: &Path, requested_root: &Path) -> Result<PathBuf> {
    ensure!(
        !requested_root
            .components()
            .any(|component| component == Component::ParentDir),
        "persistent Rustdoc target cannot contain '..'"
    );

    let requested_root = absolute_from(root, requested_root);
    let workspace_target = root.join("target");
    fs::create_dir_all(&workspace_target)
        .wrap_err_with(|| format!("failed to create {}", workspace_target.display()))?;
    let temporary_root = env::temp_dir();

    let allowed_root = if requested_root.starts_with(&workspace_target) {
        &workspace_target
    } else {
        ensure!(
            !requested_root.starts_with(root),
            "persistent Rustdoc target must not be inside the workspace outside {}",
            workspace_target.display()
        );
        ensure!(
            requested_root.starts_with(&temporary_root) && requested_root != temporary_root,
            "persistent Rustdoc target must be a child of {} or {}",
            workspace_target.display(),
            temporary_root.display()
        );
        &temporary_root
    };

    let canonical_root = fs::canonicalize(root)
        .wrap_err_with(|| format!("failed to resolve workspace root {}", root.display()))?;
    let canonical_workspace_target = fs::canonicalize(&workspace_target)
        .wrap_err_with(|| format!("failed to resolve {}", workspace_target.display()))?;
    ensure!(
        canonical_workspace_target.starts_with(&canonical_root)
            && canonical_workspace_target != canonical_root,
        "workspace target escapes workspace root through a symlink"
    );
    let canonical_allowed_root = fs::canonicalize(allowed_root)
        .wrap_err_with(|| format!("failed to resolve {}", allowed_root.display()))?;
    let target = requested_root.join("cargo-target-v1");

    let existing = closest_existing_directory(&target)
        .context("persistent Rustdoc target has no existing directory ancestor")?;
    let canonical_existing = fs::canonicalize(&existing)
        .wrap_err_with(|| format!("failed to resolve {}", existing.display()))?;
    ensure!(
        canonical_existing.starts_with(&canonical_allowed_root),
        "persistent Rustdoc target escapes {} through a symlink",
        allowed_root.display()
    );
    ensure!(
        !canonical_existing.starts_with(&canonical_root)
            || canonical_existing.starts_with(&canonical_workspace_target),
        "persistent Rustdoc target must not resolve inside the workspace outside {}",
        workspace_target.display()
    );

    fs::create_dir_all(&target)
        .wrap_err_with(|| format!("failed to create {}", target.display()))?;
    let canonical_target = fs::canonicalize(&target)
        .wrap_err_with(|| format!("failed to resolve {}", target.display()))?;
    ensure!(
        canonical_target.starts_with(&canonical_allowed_root),
        "persistent Rustdoc target escapes {} through a symlink",
        allowed_root.display()
    );
    ensure!(
        !canonical_target.starts_with(&canonical_root)
            || canonical_target.starts_with(&canonical_workspace_target),
        "persistent Rustdoc target must not resolve inside the workspace outside {}",
        workspace_target.display()
    );
    Ok(target)
}

fn closest_existing_directory(path: &Path) -> Option<PathBuf> {
    path.ancestors()
        .find(|candidate| candidate.is_dir())
        .map(Path::to_path_buf)
}

fn clear_persistent_rustdoc_output(target: &Path) -> Result<()> {
    let rustdoc = target.join("doc");
    let metadata = match fs::symlink_metadata(&rustdoc) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error.into()),
    };

    if metadata.file_type().is_symlink() {
        remove_symlink(&rustdoc)
    } else if metadata.is_dir() {
        let canonical_target = fs::canonicalize(target)
            .wrap_err_with(|| format!("failed to resolve {}", target.display()))?;
        let canonical_rustdoc = fs::canonicalize(&rustdoc)
            .wrap_err_with(|| format!("failed to resolve {}", rustdoc.display()))?;
        ensure!(
            canonical_rustdoc.starts_with(&canonical_target)
                && canonical_rustdoc != canonical_target,
            "cached Rustdoc output escapes {}",
            target.display()
        );
        fs::remove_dir_all(&rustdoc)
    } else {
        fs::remove_file(&rustdoc)
    }
    .wrap_err_with(|| {
        format!(
            "failed to clear cached Rustdoc output {}",
            rustdoc.display()
        )
    })
}

#[cfg(not(windows))]
fn remove_symlink(path: &Path) -> std::io::Result<()> {
    fs::remove_file(path)
}

#[cfg(windows)]
fn remove_symlink(path: &Path) -> std::io::Result<()> {
    if path.is_dir() {
        fs::remove_dir(path)
    } else {
        fs::remove_file(path)
    }
}

fn ensure_safe_output(root: &Path, output: &Path) -> Result<()> {
    ensure!(
        output != Path::new("/"),
        "refusing to write documentation to /"
    );
    ensure!(
        output != root,
        "documentation output cannot be the workspace root"
    );
    ensure!(
        !output
            .components()
            .any(|component| component == Component::ParentDir),
        "documentation output cannot contain '..'"
    );
    Ok(())
}

fn clear_stale_staging(output: &Path) -> Result<()> {
    let staging = output.join(".staging");
    let metadata = match fs::symlink_metadata(&staging) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error.into()),
    };
    if metadata.is_dir() && !metadata.file_type().is_symlink() {
        fs::remove_dir_all(&staging)
    } else {
        fs::remove_file(&staging)
    }
    .wrap_err_with(|| format!("failed to clear stale staging path {}", staging.display()))
}

fn ensure_relative_path(path: &Path) -> Result<()> {
    ensure!(
        !path.is_absolute(),
        "{} must be workspace-relative",
        path.display()
    );
    ensure!(
        !path
            .components()
            .any(|component| component == Component::ParentDir),
        "{} cannot escape the workspace",
        path.display()
    );
    Ok(())
}

fn validate_route_segment(value: &str) -> Result<()> {
    ensure!(!value.is_empty(), "route segment cannot be empty");
    ensure!(
        value
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || "._+-".contains(character)),
        "unsafe route segment {value}"
    );
    Ok(())
}

fn validate_page_route(route: &str) -> Result<()> {
    if route.is_empty() {
        return Ok(());
    }
    ensure!(
        route.ends_with('/'),
        "page route must end with '/': {route}"
    );
    ensure!(
        !route.starts_with('/'),
        "page route must be relative: {route}"
    );
    for segment in route.trim_end_matches('/').split('/') {
        validate_route_segment(segment)?;
    }
    ensure!(
        !route.contains(".."),
        "page route cannot contain '..': {route}"
    );
    Ok(())
}

fn validate_typst_symbol(symbol: &str) -> Result<()> {
    ensure!(!symbol.is_empty(), "Typst export symbol cannot be empty");
    ensure!(
        symbol
            .chars()
            .next()
            .is_some_and(|character| character.is_ascii_alphabetic() || character == '_'),
        "unsafe Typst export symbol {symbol}"
    );
    ensure!(
        symbol.chars().all(|character| {
            character.is_ascii_alphanumeric() || matches!(character, '_' | '-')
        }),
        "unsafe Typst export symbol {symbol}"
    );
    Ok(())
}

fn looks_like_version(value: &str) -> bool {
    let core = value.split(['+', '-']).next().unwrap_or(value);
    let parts = core.split('.').collect::<Vec<_>>();
    parts.len() == 3
        && parts
            .iter()
            .all(|part| !part.is_empty() && part.chars().all(|c| c.is_ascii_digit()))
}

fn typst_string(value: &str) -> String {
    value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
        .replace('\t', "\\t")
}

fn provenance_typst(metadata: &SnapshotMetadata<'_>) -> String {
    let channel = match metadata.channel {
        BuildChannel::Latest => "latest".to_owned(),
        BuildChannel::Snapshot => format!(
            "snapshot {}",
            metadata.snapshot_tag.unwrap_or("missing snapshot tag")
        ),
    };
    let mut components = String::new();
    for component in &metadata.components {
        let features = if component.features.is_empty() {
            "no optional features".to_owned()
        } else {
            format!("features: {}", component.features.join(", "))
        };
        components.push_str(&format!(
            "- #raw(\"{}/{}\") — #raw(\"{}\") (#raw(\"{}\"))\n",
            typst_string(component.product),
            typst_string(component.package),
            typst_string(&component.version),
            typst_string(&features),
        ));
    }
    format!(
        "= Version information <snapshot-version-metadata>\n\
         Save these identifiers with published results and include them in issue reports so that\n\
         collaborators can reproduce the same software environment.\n\n\
         - *Documentation channel:* #raw(\"{}\")\n\
         - *Documentation route:* #raw(\"{}\")\n\
         - *Source revision:* #raw(\"{}\")\n\n\
         == Package versions and enabled features\n\
         {components}",
        typst_string(&channel),
        typst_string(&metadata.route),
        typst_string(&metadata.git_commit),
    )
}

fn generated_reference_typst(
    product: &str,
    gamma: &GammaLoopReference,
    vakint: &VakintReference,
) -> String {
    match product {
        "gammaloop" => {
            let mut commands = String::new();
            for command in &gamma.commands {
                commands.push_str(&format!(
                    "[#raw(\"{}\")], [#raw(\"{}\")],\n",
                    typst_string(&command.path),
                    typst_string(&command.about),
                ));
            }
            let mut settings = String::new();
            for setting in &gamma.settings {
                let default = setting
                    .default
                    .as_ref()
                    .map(Value::to_string)
                    .unwrap_or_else(|| "no serialized default".to_owned());
                settings.push_str(&format!(
                    "[#raw(\"{}\")], [#raw(\"{}\")], [#raw(\"{}\")],\n",
                    typst_string(&setting.path),
                    typst_string(&setting.value_type),
                    typst_string(&default),
                ));
            }
            format!(
                "= CLI commands and settings <generated-cli-settings>\n\
                 Commands and settings available in this version are listed below.\n\n\
                 == Command tree\n\
                 #table(columns: (2fr, 3fr), table.header([*Command*], [*Description*]), {commands})\n\n\
                 == Settings\n\
                 #table(columns: (2fr, 1fr, 2fr), table.header([*Path*], [*Type*], [*Default*]), {settings})"
            )
        }
        "vakint" => {
            let mut dependencies = String::new();
            for dependency in &vakint.dependencies {
                dependencies.push_str(&format!(
                    "[#raw(\"{}\")], [#raw(\"{}\")],\n",
                    typst_string(&dependency.name),
                    typst_string(&dependency.minimum_version),
                ));
            }
            let mut topologies = String::new();
            for topology in &vakint.topologies {
                topologies.push_str(&format!(
                    "[#raw(\"{}\")], [{}], [{}],\n",
                    typst_string(&topology.name),
                    topology.loops,
                    topology.propagator_slots,
                ));
            }
            format!(
                "= Topologies and dependencies <generated-vakint-reference>\n\
                 This version supports the topology patterns and external-tool versions below.\n\n\
                 == External dependencies\n\
                 #table(columns: (1fr, 1fr), table.header([*Dependency*], [*Minimum version*]), {dependencies})\n\n\
                 == Supported topology patterns\n\
                 #table(columns: (2fr, 1fr, 1fr), table.header([*Name*], [*Loops*], [*Propagator slots*]), {topologies})"
            )
        }
        _ => String::new(),
    }
}

fn slug(value: &str) -> String {
    value
        .chars()
        .flat_map(char::to_lowercase)
        .map(|character| {
            if character.is_ascii_alphanumeric() {
                character
            } else {
                '-'
            }
        })
        .collect::<String>()
        .split('-')
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>()
        .join("-")
}

fn generated_reference_redirect(product: &str, target: &str, label: &str) -> String {
    format!(
        "<!doctype html><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"0; url={}\"><link rel=\"canonical\" href=\"{}\"><title>{} · {}</title><a href=\"{}\">Open the {} for {}</a>",
        escape_html(target),
        escape_html(target),
        escape_html(product),
        escape_html(label),
        escape_html(target),
        escape_html(label),
        escape_html(product),
    )
}

fn reference_page(product: &str, title: &str, body: &str) -> String {
    format!(
        "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\"><title>{} · {}</title><style>body{{max-width:75rem;margin:auto;padding:2rem 1rem;color:#172033;font:16px/1.5 system-ui}}a{{color:#5b4bdb}}pre{{overflow:auto;padding:1rem;background:#f1f3f7}}table{{width:100%;border-collapse:collapse;margin:1rem 0 2rem}}th,td{{padding:.45rem .6rem;border:1px solid #d9dce5;text-align:left;vertical-align:top}}th{{background:#f1f3f7}}section{{padding-top:.5rem}}</style></head><body><nav><a href=\"../../\">Manual</a></nav><h1>{} · {}</h1>{body}</body></html>",
        escape_html(product),
        escape_html(title),
        escape_html(product),
        escape_html(title),
    )
}

fn render_gammaloop_generated_reference(product: &str, reference: &GammaLoopReference) -> String {
    let mut body = format!(
        "<p>Commands and settings available in this version are listed below. <a href=\"reference/generated/gammaloop-reference.json\">Download JSON for tooling</a>.</p><p>{} commands · {} settings</p>",
        reference.commands.len(),
        reference.settings.len()
    );
    body.push_str("<h2>Command tree</h2>");
    for command in &reference.commands {
        let aliases = if command.aliases.is_empty() {
            String::new()
        } else {
            format!(
                "<p>Aliases: <code>{}</code></p>",
                escape_html(&command.aliases.join(", "))
            )
        };
        let mut arguments = String::new();
        for argument in &command.arguments {
            let mut names = Vec::new();
            if let Some(short) = argument.short {
                names.push(format!("-{short}"));
            }
            if let Some(long) = &argument.long {
                names.push(format!("--{long}"));
            }
            if names.is_empty() {
                names.push(argument.id.clone());
            }
            let values = if argument.value_names.is_empty() {
                String::new()
            } else {
                format!(" &lt;{}&gt;", escape_html(&argument.value_names.join("|")))
            };
            arguments.push_str(&format!(
                "<tr><td><code>{}{values}</code></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                escape_html(&names.join(", ")),
                escape_html(&argument.help),
                if argument.required { "yes" } else { "no" },
                escape_html(&argument.defaults.join(", ")),
                escape_html(&argument.possible_values.join(", ")),
            ));
        }
        let table = if arguments.is_empty() {
            "<p>No arguments.</p>".to_owned()
        } else {
            format!(
                "<table><thead><tr><th>Argument</th><th>Description</th><th>Required</th><th>Default</th><th>Values</th></tr></thead><tbody>{arguments}</tbody></table>"
            )
        };
        body.push_str(&format!(
            "<section id=\"{}\"><h3><code>{}</code></h3><p>{}</p>{aliases}{table}</section>",
            generated_anchor("command", &command.path),
            escape_html(&command.path),
            escape_html(&command.about),
        ));
    }
    body.push_str("<h2>Settings</h2><table><thead><tr><th>Path</th><th>Type</th><th>Description</th><th>Required</th><th>Default</th><th>Values</th></tr></thead><tbody>");
    for setting in &reference.settings {
        let default = setting
            .default
            .as_ref()
            .map(Value::to_string)
            .unwrap_or_default();
        body.push_str(&format!(
            "<tr id=\"{}\"><td><code>{}</code></td><td>{}</td><td>{}</td><td>{}</td><td><code>{}</code></td><td>{}</td></tr>",
            generated_anchor("setting", &setting.path),
            escape_html(&setting.path),
            escape_html(&setting.value_type),
            escape_html(&setting.description),
            if setting.required { "yes" } else { "no" },
            escape_html(&default),
            escape_html(&setting.possible_values.join(", ")),
        ));
    }
    body.push_str("</tbody></table>");
    reference_page(product, "CLI commands and settings", &body)
}

fn render_vakint_generated_reference(product: &str, reference: &VakintReference) -> String {
    let mut body = "<p>This version supports the external tools and topology patterns below. <a href=\"reference/generated/vakint-reference.json\">Download JSON for tooling</a>.</p><h2 id=\"external-dependencies\">External dependencies</h2><table><thead><tr><th>Dependency</th><th>Minimum version</th></tr></thead><tbody>".to_owned();
    for dependency in &reference.dependencies {
        body.push_str(&format!(
            "<tr><td>{}</td><td><code>{}</code></td></tr>",
            escape_html(&dependency.name),
            escape_html(&dependency.minimum_version),
        ));
    }
    body.push_str("</tbody></table><h2>Supported topology patterns</h2><table><thead><tr><th>Name</th><th>Loops</th><th>Top-level propagator slots</th></tr></thead><tbody>");
    for topology in &reference.topologies {
        body.push_str(&format!(
            "<tr id=\"{}\"><td><code>{}</code></td><td>{}</td><td>{}</td></tr>",
            generated_anchor("topology", &topology.name),
            escape_html(&topology.name),
            topology.loops,
            topology.propagator_slots,
        ));
    }
    body.push_str("</tbody></table>");
    reference_page(product, "Topologies and dependencies", &body)
}

fn generated_anchor(kind: &str, key: &str) -> String {
    let hash = key.bytes().fold(0xcbf29ce484222325_u64, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
    });
    format!("{kind}-{}-{hash:016x}", slug(key))
}

fn validate_generated_anchors(values: impl Iterator<Item = String>) -> Result<()> {
    let mut seen = BTreeSet::new();
    for value in values {
        ensure!(
            seen.insert(value.clone()),
            "duplicate generated anchor {value}"
        );
    }
    Ok(())
}

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn copy_tree(source: &Path, destination: &Path) -> Result<()> {
    fs::create_dir_all(destination)?;
    for entry in WalkDir::new(source) {
        let entry = entry?;
        let relative = entry.path().strip_prefix(source)?;
        if relative.as_os_str().is_empty() {
            continue;
        }
        let target = destination.join(relative);
        if entry.file_type().is_dir() {
            fs::create_dir_all(&target)?;
        } else if entry.file_type().is_file() {
            if let Some(parent) = target.parent() {
                fs::create_dir_all(parent)?;
            }
            fs::copy(entry.path(), &target)?;
        } else {
            bail!(
                "unsupported documentation artifact {}",
                entry.path().display()
            );
        }
    }
    Ok(())
}

fn replace_generated_tree(source: &Path, destination: &Path) -> Result<()> {
    let metadata = match fs::symlink_metadata(destination) {
        Ok(metadata) => Some(metadata),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => return Err(error.into()),
    };
    if let Some(metadata) = metadata {
        if metadata.is_dir() && !metadata.file_type().is_symlink() {
            fs::remove_dir_all(destination)
        } else if metadata.file_type().is_symlink() {
            remove_symlink(destination)
        } else {
            fs::remove_file(destination)
        }
        .wrap_err_with(|| format!("failed to replace {}", destination.display()))?;
    }
    copy_tree(source, destination)
}

fn directory_files(root: &Path) -> Result<BTreeMap<PathBuf, Vec<u8>>> {
    let mut files = BTreeMap::new();
    for entry in WalkDir::new(root) {
        let entry = entry?;
        if entry.file_type().is_file() {
            files.insert(
                entry.path().strip_prefix(root)?.to_path_buf(),
                fs::read(entry.path())?,
            );
        }
    }
    Ok(files)
}

fn directories_equal(left: &Path, right: &Path) -> Result<bool> {
    Ok(directory_files(left)? == directory_files(right)?)
}

fn validate_local_target(
    output: &Path,
    source: &Path,
    href: &str,
    failures: &mut Vec<String>,
) -> Result<()> {
    let Some(resolved) = resolve_local_link(output, source, href) else {
        return Ok(());
    };
    let source = source.strip_prefix(output).unwrap_or(source).display();
    if !resolved.is_file() {
        failures.push(format!("{source} -> {href}"));
        return Ok(());
    }
    let Some((_, fragment)) = href.split_once('#') else {
        return Ok(());
    };
    let fragment = fragment.split('?').next().unwrap_or_default();
    if fragment.is_empty() {
        return Ok(());
    }
    let html = fs::read_to_string(&resolved)?;
    let id_double = format!("id=\"{fragment}\"");
    let id_single = format!("id='{fragment}'");
    let name_double = format!("name=\"{fragment}\"");
    let name_single = format!("name='{fragment}'");
    if ![id_double, id_single, name_double, name_single]
        .iter()
        .any(|anchor| html.contains(anchor))
    {
        failures.push(format!("{source} -> {href} (missing fragment)"));
    }
    Ok(())
}

fn decode_html_text(value: &str) -> String {
    value
        .replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&#39;", "'")
}

fn resolve_local_link(output: &Path, source: &Path, href: &str) -> Option<PathBuf> {
    if href.is_empty()
        || href.starts_with("//")
        || href
            .split(':')
            .next()
            .is_some_and(|scheme| matches!(scheme, "http" | "https" | "mailto" | "data"))
    {
        return None;
    }
    let href = href.split(['?', '#']).next().unwrap_or_default();
    let relative = Path::new(href);
    let joined = if href.is_empty() {
        source.to_path_buf()
    } else if relative.is_absolute() {
        output.join(relative.strip_prefix("/").ok()?)
    } else {
        source.parent()?.join(relative)
    };
    let mut normalized = PathBuf::new();
    for component in joined.components() {
        match component {
            Component::Prefix(prefix) => normalized.push(prefix.as_os_str()),
            Component::RootDir => normalized.push(Path::new("/")),
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            Component::Normal(part) => normalized.push(part),
        }
    }
    if href.ends_with('/') || normalized.is_dir() {
        normalized.push("index.html");
    }
    Some(normalized)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snapshot_tags_are_safe_semver_routes() {
        assert!(looks_like_version("0.17.0"));
        assert!(looks_like_version("1.2.3-rc.1"));
        assert!(!looks_like_version("main"));
        assert!(validate_route_segment("linnet-v0.17.0").is_ok());
        assert!(validate_route_segment("../latest").is_err());
    }

    #[test]
    fn immutable_snapshot_navigation_does_not_link_mutable_developer_notes() {
        let builder = SiteBuilder::discover().unwrap();
        let product = &builder.registry.product[0];
        let current = SitePage::new("", &product.title, "Overview");
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Snapshot,
            snapshot_tag: Some("v0.3.4"),
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/gammaloop/snapshots/v0.3.4/".to_owned(),
            components: Vec::new(),
        };
        let sidebar = builder.site_sidebar(product, &metadata, &current, "");
        assert!(!sidebar.contains("developers/"));
        assert!(!sidebar.contains("For developers"));
    }

    #[test]
    fn every_registered_component_tag_family_matches_its_manifest_version() {
        let builder = SiteBuilder::discover().unwrap();
        for product in &builder.registry.product {
            for component in product
                .rust_components
                .iter()
                .chain(&product.python_components)
            {
                let version = builder.component_version(component).unwrap();
                for prefix in [&component.package, &component.id] {
                    builder
                        .validate_snapshot_tag(&format!("{prefix}-v{version}"))
                        .unwrap();
                }
            }
        }
        let gamma = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "gammaloop")
            .unwrap()
            .python_components
            .iter()
            .find(|component| component.package == "gammaloop")
            .unwrap();
        builder
            .validate_snapshot_tag(&format!("v{}", builder.component_version(gamma).unwrap()))
            .unwrap();
        assert!(builder.validate_snapshot_tag("linnet-v0.0.0").is_err());
    }

    #[test]
    fn slugs_are_stable() {
        assert_eq!(slug("Dirac & color algebra"), "dirac-color-algebra");
    }

    #[test]
    fn canonical_markdown_changelogs_render_headings_links_lists_and_code() {
        let rendered = markdown_changelog_to_typst(
            "# Changelog\n\n## [1.2.3](https://example.test/release)\n\n- fixed it\n\n```rust\nlet value = 1;\n```\n",
        );
        assert!(rendered.contains("== #raw(\"Changelog\")"));
        assert!(rendered.contains("=== #link(\"https://example.test/release\")[#raw(\"1.2.3\")]"));
        assert!(rendered.contains("- #raw(\"fixed it\")"));
        assert!(rendered.contains("block: true, lang: \"rust\""));
    }

    #[test]
    fn developer_markdown_renders_tables_and_escapes_raw_html() {
        let rendered = render_developer_markdown(
            "## Shape\n\n| input | result |\n| --- | --- |\n| tensor | scalar |\n\n<script>alert('no')</script>\n",
        );
        assert!(rendered.contains("<table>"));
        assert!(rendered.contains("<th>input</th>"));
        assert!(rendered.contains("&lt;script&gt;alert('no')&lt;/script&gt;"));
        assert!(!rendered.contains("<script>"));
    }

    #[test]
    fn developer_html_diagrams_reject_active_content() {
        assert!(validate_developer_html_note("<main><a href=\"guide\">Guide</a></main>").is_ok());
        for source in [
            "<script>alert(1)</script>",
            "<main onload=\"alert(1)\"></main>",
            "<svg/onload=alert(1)></svg>",
            "<a href=\"javascript:alert(1)\">Open</a>",
            "<a href=\"   javascript:alert(1)\">Open</a>",
            "<a href=javascript:alert(1)>Open</a>",
            "<a href=\"java&#x73;cript:alert(1)\">Open</a>",
            "<iframe src=\"page.html\"></iframe>",
            "<meta http-equiv=\"refresh\" content=\"0; url=/\">",
            "<style>main { background: url(javascript:alert(1)); }</style>",
        ] {
            assert!(validate_developer_html_note(source).is_err(), "{source}");
        }
    }

    #[test]
    fn rustdoc_uses_the_derive_page_for_derive_macros() {
        assert_eq!(
            rustdoc_macro_prefix(Some(
                "#[derive(SimpleRepresentation)] #[representation(...)]"
            )),
            "derive"
        );
        assert_eq!(rustdoc_macro_prefix(Some("macro_rules! graph")), "macro");
        assert_eq!(rustdoc_macro_prefix(None), "macro");

        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "spenso")
            .unwrap();
        let component = product
            .rust_components
            .iter()
            .find(|component| component.id == "spenso-macros")
            .unwrap();
        let catalog = builder
            .component_catalog(product, ApiLanguage::Rust, component)
            .unwrap();
        let item = catalog
            .root
            .scopes
            .get("supported")
            .unwrap()
            .items
            .values()
            .find(|item| item.id == "SimpleRepresentation")
            .unwrap();
        assert_eq!(
            rustdoc_href(&catalog, item),
            "reference/rust/spenso_macros/derive.SimpleRepresentation.html"
        );
    }

    #[test]
    fn rustdoc_type_prefix_uses_the_declaration_not_doc_prose() {
        let mut structure = DocItem::new(
            "HedgeGraph",
            "HedgeGraph",
            "Half-edge graph",
            alphal00p_docs_schema::DocItemKind::Type,
        );
        structure.signature = Some(
            "#[doc = \"The type used by clients.\"] pub struct HedgeGraph<T> { value: T }".into(),
        );
        assert_eq!(rustdoc_type_prefix(&structure), "struct");

        let mut alias = DocItem::new(
            "Graph",
            "Graph",
            "Graph alias",
            alphal00p_docs_schema::DocItemKind::Type,
        );
        alias.signature =
            Some("/// A structure-shaped type. pub type Graph = HedgeGraph<()>;".into());
        assert_eq!(rustdoc_type_prefix(&alias), "type");

        let mut enumeration = DocItem::new(
            "Flow",
            "Flow",
            "Flow",
            alphal00p_docs_schema::DocItemKind::Type,
        );
        enumeration.signature = Some("/// The type of flow. pub enum Flow { Source, Sink }".into());
        assert_eq!(rustdoc_type_prefix(&enumeration), "enum");
    }

    #[test]
    fn rustdoc_sidecars_reject_invalid_documentation_links_and_markup() {
        for lint in [
            "rustdoc::broken_intra_doc_links",
            "rustdoc::invalid_html_tags",
            "rustdoc::bare_urls",
        ] {
            assert!(STRICT_RUSTDOC_FLAGS.contains(&format!("-D {lint}")));
        }
    }

    #[test]
    fn persistent_rustdoc_targets_accept_dedicated_safe_roots() {
        let temporary = tempfile::tempdir().unwrap();
        let workspace = temporary.path().join("workspace");
        fs::create_dir_all(workspace.join("target")).unwrap();

        let watch_root = workspace.join("target/watch-session/rustdoc");
        assert_eq!(
            prepare_persistent_rustdoc_target(&workspace, &watch_root).unwrap(),
            watch_root.join("cargo-target-v1")
        );

        let nix_root = temporary.path().join("nix-build/alphal00p-docs-rustdoc");
        assert_eq!(
            prepare_persistent_rustdoc_target(&workspace, &nix_root).unwrap(),
            nix_root.join("cargo-target-v1")
        );
    }

    #[test]
    fn persistent_rustdoc_targets_reject_workspace_and_arbitrary_roots() {
        let workspace = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let arbitrary = workspace
            .parent()
            .unwrap()
            .join("alphal00p-docs-unrelated-cache");

        assert!(prepare_persistent_rustdoc_target(workspace, workspace).is_err());
        assert!(prepare_persistent_rustdoc_target(workspace, &arbitrary).is_err());
        assert!(
            prepare_persistent_rustdoc_target(workspace, Path::new("target/watch/../../outside"))
                .is_err()
        );
    }

    #[cfg(unix)]
    #[test]
    fn persistent_rustdoc_target_rejects_symlink_escape() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let workspace = temporary.path().join("workspace");
        fs::create_dir_all(workspace.join("target")).unwrap();
        let repository = Path::new(env!("CARGO_MANIFEST_DIR"))
            .ancestors()
            .nth(2)
            .unwrap();
        let escape = temporary.path().join("escape");
        symlink(repository, &escape).unwrap();

        let error =
            prepare_persistent_rustdoc_target(&workspace, &escape.join("cache")).unwrap_err();
        assert!(format!("{error:#}").contains("escapes"));
    }

    #[cfg(unix)]
    #[test]
    fn cached_rustdoc_symlink_is_unlinked_without_touching_its_target() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let target = temporary.path().join("cargo-target-v1");
        let external = temporary.path().join("external-docs");
        fs::create_dir_all(&target).unwrap();
        fs::create_dir_all(&external).unwrap();
        fs::write(external.join("keep.txt"), "keep").unwrap();
        symlink(&external, target.join("doc")).unwrap();

        clear_persistent_rustdoc_output(&target).unwrap();

        assert!(external.join("keep.txt").is_file());
        assert!(fs::symlink_metadata(target.join("doc")).is_err());
    }

    #[test]
    fn local_links_resolve_lexically() {
        let output = Path::new("/tmp/site");
        let source = output.join("products/gammaloop/latest/index.html");
        assert_eq!(
            resolve_local_link(output, &source, "../../linnet/latest/").unwrap(),
            output.join("products/linnet/latest/index.html")
        );
        assert!(resolve_local_link(output, &source, "https://example.com").is_none());
    }

    #[test]
    fn nested_pages_rebase_product_relative_links() {
        assert_eq!(page_root_prefix(""), "");
        assert_eq!(page_root_prefix("tutorial/"), "../");
        assert_eq!(page_root_prefix("manual/interfaces/"), "../../");
        assert_eq!(page_root_prefix("reference/python/component/"), "../../../");
        let rendered = rewrite_page_links(
            r##"<a href="reference/python/">API</a><a href="#local">Local</a><a href="https://example.test">External</a>"##,
            "../../",
        )
        .unwrap();
        assert!(rendered.contains(r#"href="../../reference/python/""#));
        assert!(rendered.contains(r##"href="#local""##));
        assert!(rendered.contains(r#"href="https://example.test""#));

        let mut headings = "<h2>Tutorial</h2><h3>Step</h3><h4>Detail</h4>".to_owned();
        promote_heading_hierarchy(&mut headings);
        assert_eq!(headings, "<h1>Tutorial</h1><h2>Step</h2><h3>Detail</h3>");
    }

    #[test]
    fn python_catalogs_render_as_structured_reference() {
        use alphal00p_docs_schema::{
            DocComponent, DocMemberKind, DocProduct, DocText, SCHEMA_VERSION,
        };

        let mut root = DocScope::new("module", "Module exports");
        let mut item = DocItem::new(
            "Engine",
            "Engine",
            "Engine",
            alphal00p_docs_schema::DocItemKind::PythonClass,
        );
        item.summary = Some("Runs the documented workflow.".to_owned());
        item.docs = Some(DocText::new(
            DocFormat::PythonDocstring,
            "Runs the documented workflow.\n\nConstruct one engine and reuse it.",
        ));
        item.signature = Some("class Engine:".to_owned());
        let mut method = DocMember::new("run", DocMemberKind::Method);
        method.signature = Some("def run(self, value: int) -> int:".to_owned());
        item.members.push(method);
        let mut overload = DocMember::new("run", DocMemberKind::Method);
        overload.signature = Some("def run(self, value: str) -> str:".to_owned());
        item.members.push(overload);
        item.members
            .push(DocMember::new("global_data", DocMemberKind::Getter));
        item.members
            .push(DocMember::new("global_data", DocMemberKind::Setter));
        root.define_item(item).unwrap();
        let catalog = DocCatalog {
            schema_version: SCHEMA_VERSION,
            product: DocProduct::new("example", "Example"),
            component: DocComponent::new(
                "example-python",
                "example",
                "Example Python",
                "1.0.0",
                ApiLanguage::Python,
            ),
            root,
        };

        let rendered = render_python_catalog(
            &catalog,
            "example-python.pyi",
            "class Engine: ...",
            "0123456789abcdef",
        );
        assert!(rendered.contains("class=\"api-item\""));
        assert!(rendered.contains("class Engine:"));
        assert!(rendered.contains("Construct one engine and reuse it."));
        assert_eq!(rendered.matches("Runs the documented workflow.").count(), 1);
        assert!(rendered.contains("def run(self, value: int)"));
        assert!(rendered.contains("Type-checker stub source"));
        assert_eq!(rendered.matches("id=\"engine-run-method\"").count(), 1);
        assert_eq!(rendered.matches("id=\"engine-run-method-2\"").count(), 1);
        assert_eq!(
            rendered.matches("id=\"engine-global-data-getter\"").count(),
            1
        );
        assert_eq!(
            rendered.matches("id=\"engine-global-data-setter\"").count(),
            1
        );

        let mut entries = Vec::new();
        append_catalog_search(&catalog, &catalog.root, &mut entries);
        assert!(
            entries
                .iter()
                .any(|entry| entry.href == "reference/python/example-python/#engine")
        );
        assert!(
            entries.iter().any(|entry| {
                entry.href == "reference/python/example-python/#engine-run-method"
            })
        );
        assert!(
            entries.iter().any(|entry| {
                entry.href == "reference/python/example-python/#engine-run-method-2"
            })
        );
        assert!(entries.iter().any(|entry| {
            entry.href == "reference/python/example-python/#engine-global-data-getter"
        }));
        assert!(entries.iter().any(|entry| {
            entry.href == "reference/python/example-python/#engine-global-data-setter"
        }));
    }

    #[test]
    fn stale_staging_is_cleared_without_touching_published_files() {
        let temporary = tempfile::tempdir().unwrap();
        let output = temporary.path().join("site");
        let staging = output.join(".staging");
        fs::create_dir_all(staging.join("nested")).unwrap();
        fs::write(staging.join("nested/incomplete.html"), "partial").unwrap();
        fs::write(output.join("published.html"), "last good generation").unwrap();

        clear_stale_staging(&output).unwrap();

        assert!(!staging.exists());
        assert_eq!(
            fs::read_to_string(output.join("published.html")).unwrap(),
            "last good generation"
        );
    }

    #[cfg(unix)]
    #[test]
    fn stale_staging_symlink_is_unlinked_without_removing_its_target() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let output = temporary.path().join("site");
        let external = temporary.path().join("external");
        fs::create_dir_all(&output).unwrap();
        fs::create_dir_all(&external).unwrap();
        fs::write(external.join("keep.txt"), "keep").unwrap();
        symlink(&external, output.join(".staging")).unwrap();

        clear_stale_staging(&output).unwrap();

        assert!(external.join("keep.txt").is_file());
        assert!(fs::symlink_metadata(output.join(".staging")).is_err());
    }

    #[test]
    fn checked_in_registry_is_valid() {
        SiteBuilder::discover().unwrap().check().unwrap();
    }

    #[test]
    fn portal_focuses_on_projects_and_keeps_dedicated_routes() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        builder
            .write_portal(output.path(), BuildChannel::Latest, None)
            .unwrap();

        let html = fs::read_to_string(output.path().join("index.html")).unwrap();
        assert_eq!(html.matches("class=\"portal-project-card\"").count(), 5);
        assert!(html.contains("id=\"projects\""));
        assert!(html.contains("Projects &amp; crates"));
        assert!(html.contains("αLoop"));
        assert!(html.contains("href=\"people/\""));
        assert!(html.contains("href=\"developers/\""));
        assert!(html.contains("href=\"publications/\""));
        assert!(html.contains(
            "rel=\"icon\" type=\"image/svg+xml\" href=\"assets/local-unitarity-light.svg\""
        ));
        assert!(html.contains(
            "href=\"assets/local-unitarity-dark.svg\" media=\"(prefers-color-scheme: dark)\""
        ));
        assert_eq!(html.matches("portal-card-cite").count(), 5);
        for removed in [
            "class=\"portal-pillars\"",
            "id=\"developers\"",
            "id=\"people\"",
            "id=\"affiliations\"",
            "href=\"#affiliations\"",
            "Publicly funded research",
        ] {
            assert!(
                !html.contains(removed),
                "unexpected landing-page copy: {removed}"
            );
        }
        assert!(!html.contains("Product manual"));
        assert!(!html.contains("scientific-computing products"));

        for asset in [
            "local-unitarity-light.svg",
            "local-unitarity-dark.svg",
            "gammalooplogo-light.svg",
            "gammalooplogo-dark.svg",
            "spensologo.svg",
            "publications.json",
        ] {
            assert!(output.path().join("assets").join(asset).is_file());
        }
        assert!(output.path().join("assets/people/valentin.webp").is_file());

        let people = fs::read_to_string(output.path().join("people/index.html")).unwrap();
        assert_eq!(people.matches("class=\"people-card\"").count(), 7);
        assert!(people.contains("https://symbolica.io/about.html"));
        assert!(people.contains("id=\"cedric-sigrist\""));
        assert!(people.contains("People building GammaLoop"));
        assert!(people.contains(
            "Researchers and collaborators developing GammaLoop, Local Unitarity methods"
        ));
        assert!(people.contains("href=\"../assets/local-unitarity-light.svg\""));
        assert!(!people.contains("People building αLoop"));
        for internal_copy in [
            "Contributor scope",
            "default branch",
            "malformed commit identities",
            "Portrait source",
        ] {
            assert!(
                !people.contains(internal_copy),
                "unexpected process copy: {internal_copy}"
            );
        }

        let css = fs::read_to_string(output.path().join("assets/site.css")).unwrap();
        assert!(css.contains(".publication-card:nth-child(2n of :not([hidden]))"));
        assert!(!css.contains(".publication-card:nth-child(2n) {"));
        assert!(css.contains(".people-card-portrait { object-position: left center; }"));
        assert!(css.contains("#ben-ruijl > .people-card-portrait"));
        assert!(css.contains(".product-logo-gammaloop { aspect-ratio:"));
        assert!(css.contains(".product-logo-spenso { aspect-ratio: 637 / 189;"));

        let publications =
            fs::read_to_string(output.path().join("publications/index.html")).unwrap();
        assert_eq!(
            publications.matches("data-publication data-title").count(),
            builder.publications.publications.len()
        );
        assert!(publications.contains("data-publication-author"));
        assert!(publications.contains("Open the live INSPIRE query"));

        let citations = fs::read_to_string(output.path().join("citations/index.html")).unwrap();
        assert_eq!(citations.matches("class=\"citation-card\"").count(), 5);
        assert!(citations.contains("10.5281/zenodo.18429583"));
        assert!(citations.contains("No registered software DOI is currently configured"));
    }

    #[test]
    fn every_product_has_first_class_structured_references() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        for product in &builder.registry.product {
            let site = output.path().join(&product.id);
            fs::create_dir_all(&site).unwrap();
            builder.write_site_assets(&site).unwrap();
            builder.write_component_catalogs(product, &site).unwrap();
            builder.write_generated_reference(product, &site).unwrap();
            let metadata = builder
                .metadata(product, BuildChannel::Latest, None, Path::new("latest"))
                .unwrap();
            builder
                .write_python_reference(product, &metadata, &site)
                .unwrap();
            builder.write_rustdoc_placeholder(product, &site).unwrap();
            builder
                .write_rust_reference(product, &metadata, &site)
                .unwrap();
            builder.write_reference_hub(product, &site).unwrap();

            let hub = fs::read_to_string(site.join("reference/index.html")).unwrap();
            assert!(hub.contains("Rust API"), "{}", product.id);
            assert!(hub.contains("Python API"), "{}", product.id);
            for component in &product.rust_components {
                let structured = site
                    .join("reference/rust/supported")
                    .join(&component.id)
                    .join("index.html");
                assert!(structured.is_file(), "{}:{}", product.id, component.id);
                assert!(
                    fs::read_to_string(structured)
                        .unwrap()
                        .contains("class=\"api-item\""),
                    "{}:{}",
                    product.id,
                    component.id
                );
            }
        }
        assert!(
            output
                .path()
                .join("gammaloop/reference/cli/index.html")
                .is_file()
        );
        assert!(
            output
                .path()
                .join("vakint/reference/topologies/index.html")
                .is_file()
        );
        assert!(product_hero(&builder.registry.product[0]).contains("product-logo-gammaloop"));
        assert!(product_hero(&builder.registry.product[2]).contains("product-logo-spenso"));
    }

    #[test]
    fn developer_architecture_is_classified_searchable_and_source_backed() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let stale = output
            .path()
            .join("developers/architecture/removed-note/index.html");
        fs::create_dir_all(stale.parent().unwrap()).unwrap();
        fs::write(&stale, "stale route").unwrap();
        builder.write_developer_docs(output.path()).unwrap();

        let developer_root = output.path().join("developers");
        assert!(!stale.exists());
        let hub = fs::read_to_string(developer_root.join("index.html")).unwrap();
        assert!(hub.contains("For developers"));
        assert!(hub.contains("href=\"assets/local-unitarity-light.svg\""));
        assert!(hub.contains("Implemented architecture"));
        assert!(hub.contains("Design proposals"));
        assert!(hub.contains("Investigation record"));
        assert!(hub.contains("Performance investigation"));
        assert_eq!(hub.matches("class=\"developer-card\"").count(), 11);

        let current = fs::read_to_string(
            developer_root.join("architecture/gammaloop-architecture/index.html"),
        )
        .unwrap();
        assert!(current.contains("developer-status-implemented"));
        assert!(current.contains("href=\"../uv-renormalization/\""));
        assert!(current.contains("View the source note"));
        assert!(current.contains("class=\"docs-toc\""));
        assert!(current.contains("class=\"developer-inline-toc\""));
        assert!(current.contains("aria-label=\"Developer note pagination\""));

        let investigation = fs::read_to_string(
            developer_root.join("architecture/spenso-improvement-effects/index.html"),
        )
        .unwrap();
        assert!(investigation.contains("<table>"));
        assert!(investigation.contains("developer-status-investigation-record"));

        let diagram = fs::read_to_string(
            developer_root.join("architecture/schoonschip-network/diagram.html"),
        )
        .unwrap();
        assert!(diagram.contains("developer-diagram-banner"));
        assert!(diagram.contains("For developers"));
        assert!(diagram.contains("/blob/"));
        assert!(!diagram.contains("../../crates/idenso/src/"));

        let search = serde_json::from_slice::<Vec<SearchEntry>>(
            &fs::read(developer_root.join("search-index.json")).unwrap(),
        )
        .unwrap();
        assert!(search.iter().any(|entry| {
            entry.href == "architecture/gammaloop-architecture/"
                && entry.kind == "developer implemented"
        }));
        assert!(search.iter().any(|entry| {
            entry.href == "architecture/gammaloop-architecture/#configuration-architecture"
        }));
    }
}
