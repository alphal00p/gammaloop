//! Build and validate the five independently rendered product manuals.

mod server;
mod watch;

pub use watch::WatchRequest;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    env, fs,
    ops::Range,
    path::{Component, Path, PathBuf},
    process::Command,
};

use alphal00p_docs_schema::{
    ApiLanguage, DocCatalog, DocFormat, DocItem, DocMember, DocMemberKind, DocScope,
    SCHEMA_VERSION,
    generated::{
        CliAliasKind, CliArgument, CliArgumentAction, CliCommand, GENERATED_REFERENCE_SCHEMA,
        GammaLoopReference, SettingReference, VakintReference,
    },
};
use chrono::{NaiveDate, Utc};
use clap::ValueEnum;
use eyre::{Context, ContextCompat, Result, bail, ensure};
use percent_encoding::{NON_ALPHANUMERIC, utf8_percent_encode};
use regex::Regex;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use tempfile::{Builder as TempDirBuilder, TempDir};
use walkdir::WalkDir;

const PRODUCT_IDS: [&str; 5] = ["gammaloop", "linnet", "spenso", "idenso", "vakint"];
const PORTAL_TASKS: [(&str, &str, &str); 5] = [
    ("Calculations", "gammaloop", "Calculation application"),
    ("Graphs", "linnet", "Graph-algorithm library"),
    ("Tensors", "spenso", "Tensor-network library"),
    ("Identities", "idenso", "Symbolic-identity library"),
    ("Integrals", "vakint", "Integral-evaluation library"),
];
const PORTAL_GRAPH_IDS: [&str; 11] = [
    "aa-2l-gl00",
    "aa-2l-gl08",
    "aa-3l-gl000",
    "aa-3l-gl150",
    "aa-3l-gl300",
    "gg-hhh-3l",
    "gg-hhh-1l",
    "qq-aaa-pentabox",
    "ad-ad-gluon",
    "epem-bbx",
    "epem-ttbar-cut",
];
const PORTAL_SCHEMA_VERSION: u32 = 2;
const TALKS_SCHEMA_VERSION: u32 = 1;
const DEVELOPER_SCHEMA_VERSION: u32 = 3;
const LEGACY_PROSE_SCHEMA_VERSION: u32 = 1;
const PUBLICATION_SCHEMA_VERSION: u32 = 1;
const PUBLISHED_DOCS_ROOT: &str = "https://alphal00p.github.io/gammaloop";
const STRICT_RUSTDOC_FLAGS: &str = "-D rustdoc::broken_intra_doc_links \
    -D rustdoc::invalid_html_tags -D rustdoc::bare_urls";

fn portal_graph_assets() -> impl Iterator<Item = String> {
    PORTAL_GRAPH_IDS.into_iter().flat_map(|graph| {
        ["light", "dark"]
            .into_iter()
            .map(move |theme| format!("portal-graph-{graph}-{theme}.svg"))
    })
}

fn about_assets() -> impl Iterator<Item = String> {
    ["double-triangle", "local-unitarity-equation"]
        .into_iter()
        .flat_map(|asset| {
            ["light", "dark"]
                .into_iter()
                .map(move |theme| format!("about-{asset}-{theme}.svg"))
        })
}

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
    scope: BuildScope,
    include_typst: bool,
    include_rustdoc: bool,
    rustdoc_target_root: Option<&'a Path>,
    dependency_output: Option<&'a Path>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum BuildScope {
    FullSite,
    ProductPreview,
}

type LinkedPageIndex = HashMap<PathBuf, (HashSet<String>, Option<String>)>;
type LocalPathIndex = HashMap<PathBuf, Option<fs::FileType>>;
type LinkRewriteIndex = HashMap<PathBuf, Vec<LinkRewrite>>;
type LinkValidationState<'a> = (
    &'a mut LinkedPageIndex,
    &'a mut LinkRewriteIndex,
    &'a mut LocalPathIndex,
    Option<&'a str>,
);

#[derive(Clone, Debug)]
struct LinkAttribute {
    attribute_range: Range<usize>,
    target_range: Range<usize>,
    attribute: String,
}

#[derive(Clone, Debug)]
struct LinkRewrite {
    range: Range<usize>,
    expected: String,
    replacement: String,
    target_before: String,
    target_after: Option<String>,
}

#[derive(Clone, Debug)]
struct LinkReference<'a> {
    source: &'a Path,
    href: &'a str,
    attribute: Option<LinkAttribute>,
    rustdoc: bool,
}

struct LinkValidationPatterns {
    script_body: Regex,
    html_comment: Regex,
    tag: Regex,
    href: Regex,
    id: Regex,
    name: Regex,
    named_anchor: Regex,
    redirect: Regex,
}

impl LinkValidationPatterns {
    fn new() -> Result<Self> {
        Ok(Self {
            script_body: Regex::new(r"(?is)(<script\b[^>]*>).*?(</script\s*>)")?,
            html_comment: Regex::new(r"(?is)<!--.*?-->")?,
            tag: Regex::new(r"(?is)<[a-z][a-z0-9-]*\b[^>]*>")?,
            href: Regex::new(
                r#"(?i)(?:^|\s)(?:href|src)\s*=\s*(?:"([^"]+)"|'([^']+)'|([^\s"'=<>`]+))"#,
            )?,
            id: Regex::new(r#"(?i)(?:^|\s)id\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'=<>`]+))"#)?,
            name: Regex::new(r#"(?i)(?:^|\s)name\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s"'=<>`]+))"#)?,
            named_anchor: Regex::new(r"(?i)^<a\b")?,
            redirect: Regex::new(
                r#"(?i)<meta\b[^>]*\bcontent\s*=\s*["'][^"']*\burl=([^"']+)["'][^>]*>"#,
            )?,
        })
    }

    fn rendered_html(&self, html: &str) -> String {
        let mut rendered = html.as_bytes().to_vec();
        for comment in self.html_comment.find_iter(html) {
            rendered[comment.range()].fill(b' ');
        }
        for script in self.script_body.captures_iter(html) {
            let opening = script.get(1).expect("script opening capture");
            let closing = script.get(2).expect("script closing capture");
            rendered[opening.end()..closing.start()].fill(b' ');
        }
        String::from_utf8(rendered).expect("masking HTML with spaces preserves UTF-8")
    }

    fn page_metadata(&self, html: &str) -> (HashSet<String>, Option<String>) {
        let mut anchors = HashSet::new();
        for element in self.tag.find_iter(html) {
            for capture in self.id.captures_iter(element.as_str()) {
                if let Some(value) = (1..=3).find_map(|index| capture.get(index)) {
                    anchors.insert(value.as_str().to_owned());
                }
            }
            if self.named_anchor.is_match(element.as_str()) {
                for capture in self.name.captures_iter(element.as_str()) {
                    if let Some(value) = (1..=3).find_map(|index| capture.get(index)) {
                        anchors.insert(value.as_str().to_owned());
                    }
                }
            }
        }
        let redirect = self
            .redirect
            .captures(html)
            .and_then(|capture| capture.get(1))
            .map(|value| decode_html_text(value.as_str()));
        (anchors, redirect)
    }
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
struct TalksConfig {
    schema: u32,
    #[serde(default)]
    talk: Vec<TalkConfig>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct TalkConfig {
    id: String,
    date: String,
    speaker: String,
    title: String,
    event: String,
    location: String,
    event_url: String,
    slides_url: Option<String>,
    recording_url: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperConfig {
    schema: u32,
    title: String,
    summary: String,
    #[serde(default)]
    owner: Vec<DeveloperOwner>,
    #[serde(default)]
    section: Vec<DeveloperSection>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperOwner {
    id: String,
    name: String,
    contact: String,
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
    lifecycle: String,
    owner: String,
    #[serde(default)]
    reviewed_at: Option<String>,
    #[serde(default)]
    review_ref: Option<String>,
    #[serde(default)]
    products: Vec<String>,
    #[serde(default)]
    topics: Vec<String>,
    #[serde(default)]
    review_triggers: Vec<String>,
    #[serde(default)]
    scope: Vec<DeveloperScope>,
    #[serde(default)]
    implemented_by: Option<String>,
    #[serde(default)]
    superseded_by: Option<String>,
    #[serde(default)]
    evidence_revision: Option<String>,
    #[serde(default)]
    captured_at: Option<String>,
    #[serde(default)]
    evidence_command: Option<String>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct DeveloperScope {
    path: PathBuf,
    symbol: String,
    anchor: String,
    digest: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct LegacyProseConfig {
    schema: u32,
    source: Vec<PathBuf>,
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
    #[serde(default)]
    aliases: Vec<String>,
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
    #[serde(default, skip_serializing_if = "String::is_empty")]
    text: String,
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

fn python_display_module(component: &ComponentConfig) -> &str {
    component
        .module
        .as_deref()
        .filter(|module| !module.contains("._"))
        .unwrap_or(&component.package)
}

fn product_site_pages(product: &ProductConfig) -> Vec<SitePage> {
    let mut groups = BTreeMap::<String, Vec<SitePage>>::new();
    let mut group_order = Vec::new();
    for page in &product.pages {
        if !groups.contains_key(&page.group) {
            group_order.push(page.group.clone());
        }
        groups
            .entry(page.group.clone())
            .or_default()
            .push(SitePage {
                route: page.route.clone(),
                title: page.title.clone(),
                group: page.group.clone(),
            });
    }
    if !groups.contains_key("Reference") {
        group_order.push("Reference".to_owned());
    }
    let reference = groups.entry("Reference".to_owned()).or_default();
    reference.push(SitePage::new("reference/", "Reference", "Reference"));
    reference.push(SitePage::new(
        "reference/python/",
        "Python API",
        "Reference",
    ));
    reference.extend(product.python_components.iter().map(|component| {
        let module = python_display_module(component);
        SitePage::new(
            format!("reference/python/{}/", component.id),
            format!("{module} Python API"),
            "Reference",
        )
    }));
    reference.push(SitePage::new("reference/rust/", "Rust API", "Reference"));
    if let Some((route, title)) = supplemental_reference(&product.id) {
        reference.push(SitePage::new(route, title, "Reference"));
    }
    group_order
        .into_iter()
        .flat_map(|group| groups.remove(&group).unwrap_or_default())
        .collect()
}

#[derive(Clone, Debug)]
struct HeadingLink {
    level: u8,
    title: String,
    id: String,
}

struct LeadingDeveloperTitle {
    title: String,
    range: std::ops::Range<usize>,
}

pub struct SiteBuilder {
    root: PathBuf,
    api_root: PathBuf,
    registry: ProductRegistry,
    portal: PortalConfig,
    talks: TalksConfig,
    developers: DeveloperConfig,
    legacy_prose: LegacyProseConfig,
    publications: PublicationCache,
}

fn interface_guide_title(product_title: &str, surface: &str) -> String {
    match surface {
        "cli" => format!("Using {product_title} from the command line"),
        "python" => format!("Using {product_title} from Python"),
        "rust" => format!("Using {product_title} from Rust"),
        "typst" => format!("Using {product_title} from Typst"),
        _ => format!("Using {product_title}"),
    }
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
        let talks_path = root.join("docs/talks.toml");
        let source = fs::read_to_string(&talks_path)
            .wrap_err_with(|| format!("failed to read {}", talks_path.display()))?;
        let talks = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", talks_path.display()))?;
        let developers_path = root.join("docs/developers.toml");
        let source = fs::read_to_string(&developers_path)
            .wrap_err_with(|| format!("failed to read {}", developers_path.display()))?;
        let developers = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", developers_path.display()))?;
        let legacy_prose_path = root.join("docs/legacy-prose.toml");
        let source = fs::read_to_string(&legacy_prose_path)
            .wrap_err_with(|| format!("failed to read {}", legacy_prose_path.display()))?;
        let legacy_prose = toml::from_str(&source)
            .wrap_err_with(|| format!("failed to parse {}", legacy_prose_path.display()))?;
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
            talks,
            developers,
            legacy_prose,
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
            self.talks.schema == TALKS_SCHEMA_VERSION,
            "talks schema {} does not match renderer schema {}",
            self.talks.schema,
            TALKS_SCHEMA_VERSION
        );
        ensure!(
            self.developers.schema == DEVELOPER_SCHEMA_VERSION,
            "developer schema {} does not match renderer schema {}",
            self.developers.schema,
            DEVELOPER_SCHEMA_VERSION
        );
        self.check_prose_sources()?;
        ensure!(
            !self.developers.title.trim().is_empty()
                && !self.developers.summary.trim().is_empty()
                && !self.developers.owner.is_empty()
                && !self.developers.section.is_empty(),
            "developer documentation metadata must not be empty"
        );
        let mut developer_owner_ids = BTreeSet::new();
        for owner in &self.developers.owner {
            validate_route_segment(&owner.id)?;
            ensure!(
                developer_owner_ids.insert(&owner.id),
                "duplicate developer owner {}",
                owner.id
            );
            ensure!(
                !owner.name.trim().is_empty() && !owner.contact.trim().is_empty(),
                "developer owner {} is incomplete",
                owner.id
            );
        }
        let mut developer_section_ids = BTreeSet::new();
        let mut developer_note_ids = BTreeSet::new();
        let mut developer_sources = BTreeSet::new();
        let mut developer_supersession = BTreeMap::new();
        let today = Utc::now().date_naive();
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
                        && !note.status.trim().is_empty()
                        && !note.lifecycle.trim().is_empty()
                        && !note.owner.trim().is_empty(),
                    "developer note {} is incomplete",
                    note.id
                );
                ensure!(
                    matches!(
                        note.lifecycle.as_str(),
                        "current" | "proposal" | "investigation" | "archived" | "superseded"
                    ),
                    "developer note {} has unsupported lifecycle {}",
                    note.id,
                    note.lifecycle
                );
                ensure!(
                    developer_owner_ids.contains(&note.owner),
                    "developer note {} references unknown owner {}",
                    note.id,
                    note.owner
                );
                if note.owner == "unassigned" {
                    eprintln!(
                        "warning: developer note {} still requires a named owner",
                        note.id
                    );
                }
                ensure!(
                    !note.products.is_empty()
                        && note
                            .products
                            .iter()
                            .all(|product| ids.contains(product.as_str())),
                    "developer note {} has no valid applicable products",
                    note.id
                );
                ensure!(
                    !note.topics.is_empty()
                        && note.topics.iter().all(|topic| !topic.trim().is_empty()),
                    "developer note {} must declare topics",
                    note.id
                );
                let reviewed_at = note
                    .reviewed_at
                    .as_deref()
                    .map(|date| {
                        NaiveDate::parse_from_str(date, "%Y-%m-%d").wrap_err_with(|| {
                            format!("developer note {} has invalid reviewed_at {date}", note.id)
                        })
                    })
                    .transpose()?;
                ensure!(
                    note.reviewed_at.is_some() == note.review_ref.is_some(),
                    "developer note {} must pair reviewed_at with review_ref",
                    note.id
                );
                if matches!(note.lifecycle.as_str(), "current" | "proposal") {
                    if let Some(reviewed_at) = reviewed_at {
                        let age = today.signed_duration_since(reviewed_at).num_days();
                        if note.lifecycle == "current" && age > 60 {
                            eprintln!(
                                "warning: current developer note {} was reviewed {age} days ago (publication blocks after owner-ratified enforcement at 90 days)",
                                note.id
                            );
                        } else if note.lifecycle == "proposal" && age > 180 {
                            eprintln!(
                                "warning: proposal {} has had no disposition for {age} days",
                                note.id
                            );
                        }
                    } else {
                        eprintln!(
                            "warning: {} developer note {} has no review record",
                            note.lifecycle, note.id
                        );
                    }
                    if note.review_triggers.is_empty() {
                        eprintln!(
                            "warning: {} developer note {} has no review triggers",
                            note.lifecycle, note.id
                        );
                    }
                }
                for trigger in &note.review_triggers {
                    ensure!(
                        !trigger.trim().is_empty()
                            && !Path::new(trigger)
                                .components()
                                .any(|component| component == Component::ParentDir),
                        "developer note {} has unsafe review trigger {trigger}",
                        note.id
                    );
                }
                for scope in &note.scope {
                    ensure_relative_path(&scope.path)?;
                    ensure!(
                        !scope.symbol.trim().is_empty()
                            && !scope.anchor.trim().is_empty()
                            && !scope.digest.trim().is_empty(),
                        "developer note {} has an incomplete verified scope",
                        note.id
                    );
                    let source =
                        fs::read_to_string(self.root.join(&scope.path)).wrap_err_with(|| {
                            format!(
                                "developer note {} cannot read verified scope {}",
                                note.id,
                                scope.path.display()
                            )
                        })?;
                    ensure!(
                        source.matches(&scope.anchor).count() == 1,
                        "developer note {} symbol {} anchor is missing or ambiguous in {}",
                        note.id,
                        scope.symbol,
                        scope.path.display()
                    );
                    let digest = Command::new("git")
                        .current_dir(&self.root)
                        .args(["hash-object", "--"])
                        .arg(&scope.path)
                        .output()
                        .wrap_err_with(|| {
                            format!("failed to hash verified scope {}", scope.path.display())
                        })?;
                    ensure!(
                        digest.status.success(),
                        "git hash-object failed for verified scope {}",
                        scope.path.display()
                    );
                    let digest = String::from_utf8(digest.stdout)?;
                    ensure!(
                        digest.trim() == scope.digest,
                        "developer note {} verified scope {} changed (expected {}, found {}); owner review or an explicit no-impact attestation is required",
                        note.id,
                        scope.path.display(),
                        scope.digest,
                        digest.trim()
                    );
                }
                ensure!(
                    note.lifecycle != "current" || !note.scope.is_empty(),
                    "current note must declare a verified code scope"
                );
                if note.lifecycle == "superseded" {
                    let replacement = note.superseded_by.as_deref().filter(|id| !id.is_empty());
                    ensure!(
                        replacement.is_some(),
                        "superseded developer note {} must name superseded_by",
                        note.id
                    );
                    developer_supersession.insert(
                        note.id.clone(),
                        replacement.expect("replacement was checked").to_owned(),
                    );
                }
                if matches!(note.lifecycle.as_str(), "investigation" | "archived")
                    && note.evidence_revision.is_none()
                {
                    eprintln!(
                        "warning: {} developer note {} has no immutable evidence revision",
                        note.lifecycle, note.id
                    );
                }
                ensure!(
                    note.source
                        .extension()
                        .and_then(|extension| extension.to_str())
                        == Some("typ"),
                    "developer note {} must use Typst",
                    note.source.display()
                );
                self.require_file(&note.source)?;
            }
        }
        for (note, replacement) in &developer_supersession {
            ensure!(
                developer_note_ids.contains(replacement),
                "developer note {note} is superseded by missing note {replacement}"
            );
            let mut seen = BTreeSet::new();
            let mut next = note.as_str();
            while let Some(replacement) = developer_supersession.get(next) {
                ensure!(
                    seen.insert(next),
                    "developer supersession cycle contains {next}"
                );
                next = replacement;
            }
        }
        let architecture_root = self.root.join("docs/architecture");
        let mut architecture_sources = BTreeSet::new();
        for entry in fs::read_dir(&architecture_root)? {
            let path = entry?.path();
            if !path.is_file() {
                continue;
            }
            match path.extension().and_then(|extension| extension.to_str()) {
                Some("typ") => {
                    architecture_sources.insert(path.strip_prefix(&self.root)?.to_path_buf());
                }
                Some("md" | "html") => bail!(
                    "developer source {} must be migrated to Typst",
                    path.display()
                ),
                _ => {}
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
        ensure!(
            !self.talks.talk.is_empty(),
            "talk catalogue must not be empty"
        );
        let mut talk_ids = BTreeSet::new();
        for talk in &self.talks.talk {
            validate_route_segment(&talk.id)?;
            ensure!(talk_ids.insert(&talk.id), "duplicate talk {}", talk.id);
            NaiveDate::parse_from_str(&talk.date, "%Y-%m-%d")
                .wrap_err_with(|| format!("talk {} has invalid date {}", talk.id, talk.date))?;
            ensure!(
                person_ids.contains(&talk.speaker),
                "talk {} references unknown speaker {}",
                talk.id,
                talk.speaker
            );
            ensure!(
                !talk.title.trim().is_empty()
                    && !talk.event.trim().is_empty()
                    && !talk.location.trim().is_empty(),
                "talk {} is incomplete",
                talk.id
            );
            for url in std::iter::once(&talk.event_url)
                .chain(talk.slides_url.iter())
                .chain(talk.recording_url.iter())
            {
                ensure!(
                    url.starts_with("https://"),
                    "talk {} URL must use HTTPS: {url}",
                    talk.id
                );
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
        self.require_file(Path::new("docs/assets/rustdoc.css"))?;
        self.require_file(Path::new("docs/assets/STIXTwoMath-Regular.woff2"))?;
        self.require_file(Path::new("docs/assets/STIX-Two-OFL.txt"))?;
        self.require_file(Path::new("scripts/render-docs-svg-assets.sh"))?;
        for source in [
            "docs/assets/typst/theme.typ",
            "docs/assets/typst/marks/local-unitarity.typ",
            "docs/assets/typst/marks/gammaloop.typ",
            "docs/assets/typst/marks/spenso.typ",
            "assets/embedded/drawing/templates/impl/physics-edge-style.typ",
            "assets/embedded/drawing/templates/layout-core.typ",
            "assets/embedded/drawing/templates/physics-edge-style.typ",
            "docs/assets/typst/portal-graphs/figure.typ",
            "docs/assets/typst/portal-graphs/layout.typ",
            "docs/assets/typst/portal-graphs/edge-style.typ",
            "docs/assets/typst/about/double-triangle.typ",
            "docs/assets/typst/about/local-unitarity-equation.typ",
        ] {
            self.require_file(Path::new(source))?;
        }
        self.require_file(Path::new("docs/assets/local-unitarity-light.svg"))?;
        self.require_file(Path::new("docs/assets/local-unitarity-dark.svg"))?;
        for graph in PORTAL_GRAPH_IDS {
            self.require_file(
                &Path::new("docs/assets/typst/portal-graphs/graphs")
                    .join(format!("portal-graph-{graph}.typ")),
            )?;
        }
        for graph in portal_graph_assets() {
            self.require_file(&Path::new("docs/assets/graphs").join(&graph))?;
        }
        for asset in about_assets() {
            self.require_file(&Path::new("docs/assets").join(&asset))?;
        }
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
                product
                    .source
                    .extension()
                    .and_then(|extension| extension.to_str())
                    == Some("typ"),
                "{} product source {} must be canonical Typst",
                product.id,
                product.source.display()
            );
            ensure!(
                !product.pages.is_empty(),
                "{} has no authored documentation pages",
                product.id
            );
            let manual_source = fs::read_to_string(self.root.join(&product.source))?;
            ensure!(
                manual_source.contains("#let manual = product-document(")
                    && manual_source.trim_end().ends_with("#manual"),
                "{} product source must define and emit the complete PDF manual",
                product.id
            );
            let product_source_root = product
                .source
                .parent()
                .context("product source must have a parent directory")?;
            let content_import_count = manual_source
                .lines()
                .filter(|line| line.trim_start().starts_with("#import \"content/"))
                .count();
            ensure!(
                content_import_count == product.pages.len(),
                "{} PDF manual imports {content_import_count} content chapters but the registry has {}",
                product.id,
                product.pages.len()
            );
            let mut previous_manual_page = 0;
            for page in &product.pages {
                let relative_source = page
                    .source
                    .strip_prefix(product_source_root)
                    .wrap_err_with(|| {
                        format!(
                            "{}:{} source {} is outside the product source directory",
                            product.id,
                            page.id,
                            page.source.display()
                        )
                    })?
                    .to_string_lossy()
                    .replace('\\', "/");
                let import = format!("#import \"{relative_source}\": {}", page.symbol);
                ensure!(
                    manual_source.contains(&import),
                    "{} PDF manual does not import registered page {} as {}",
                    product.id,
                    page.id,
                    page.symbol
                );
                let body_entry = Regex::new(&format!(
                    r"#{}(?:-content)?(?:\(|\r?\n)",
                    regex::escape(&page.symbol)
                ))?;
                let body_entries = body_entry.find_iter(&manual_source).collect::<Vec<_>>();
                ensure!(
                    body_entries.len() == 1,
                    "{} PDF manual must include registered page {} exactly once",
                    product.id,
                    page.id
                );
                let body_entry = body_entries[0];
                ensure!(
                    body_entry.start() >= previous_manual_page,
                    "{} PDF manual chapter order differs from the page registry at {}",
                    product.id,
                    page.id
                );
                previous_manual_page = body_entry.end();
            }
            let mut page_ids = BTreeSet::new();
            let mut page_routes = BTreeSet::new();
            let mut page_sources = BTreeSet::new();
            let mut root_pages = 0;
            let mut tutorial_pages = 0;
            let mut guide_pages = 0;
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
                    page_sources.insert(page.source.clone()),
                    "{} registers page source {} more than once",
                    product.id,
                    page.source.display()
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
                for alias in &page.aliases {
                    validate_page_route(alias)?;
                    ensure!(
                        !alias.is_empty() && page_routes.insert(alias),
                        "{}:{} alias {} conflicts with another authored route or alias",
                        product.id,
                        page.id,
                        alias
                    );
                }
                self.require_file(&page.source)?;
                ensure!(
                    page.source
                        .extension()
                        .and_then(|extension| extension.to_str())
                        == Some("typ"),
                    "{}:{} source {} must be canonical Typst",
                    product.id,
                    page.id,
                    page.source.display()
                );
                root_pages += usize::from(page.route.is_empty());
                tutorial_pages += usize::from(page.group == "Tutorial");
                guide_pages += usize::from(page.group == "Guides");
            }
            let content_root = self
                .root
                .join("docs/products")
                .join(&product.id)
                .join("content");
            let mut authored_sources = BTreeSet::new();
            for entry in WalkDir::new(&content_root) {
                let entry = entry?;
                if entry.file_type().is_file()
                    && entry.path().extension().and_then(|value| value.to_str()) == Some("typ")
                {
                    authored_sources.insert(entry.path().strip_prefix(&self.root)?.to_path_buf());
                }
            }
            ensure!(
                authored_sources == page_sources,
                "{} content and page registry differ; unregistered: {:?}; missing: {:?}",
                product.id,
                authored_sources
                    .difference(&page_sources)
                    .collect::<Vec<_>>(),
                page_sources
                    .difference(&authored_sources)
                    .collect::<Vec<_>>()
            );
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
            let quickstart = product
                .pages
                .iter()
                .find(|page| page.id == "quickstart")
                .wrap_err_with(|| format!("{} has no quickstart page", product.id))?;
            ensure!(
                quickstart.title != "Quickstart"
                    && quickstart.route == "quickstart/"
                    && quickstart.group == "Start here"
                    && quickstart.symbol == "quickstart",
                "{} quickstart chooser must use an outcome-specific title and the shared route, navigation group, and Typst symbol",
                product.id
            );
            let mut quickstart_surfaces = Vec::new();
            if product.id == "gammaloop" {
                quickstart_surfaces.push(("cli", "CLI"));
            }
            if !product.python_components.is_empty() {
                quickstart_surfaces.push(("python", "Python"));
            }
            if !product.rust_components.is_empty() {
                quickstart_surfaces.push(("rust", "Rust"));
            }
            if product
                .pages
                .iter()
                .any(|page| page.route == "reference/typst/")
            {
                quickstart_surfaces.push(("typst", "Typst"));
            }
            for (surface, label) in quickstart_surfaces {
                let id = format!("quickstart-{surface}");
                let title = interface_guide_title(&product.title, surface);
                let page = product
                    .pages
                    .iter()
                    .find(|page| page.id == id)
                    .wrap_err_with(|| format!("{} has no {label} quickstart", product.id))?;
                ensure!(
                    page.title == title
                        && page.route == format!("quickstart/{surface}/")
                        && page.group == "Getting started"
                        && page.symbol == id,
                    "{} {label} getting-started guide must use the shared title, route, navigation group, and Typst symbol",
                    product.id
                );
            }
            ensure!(
                guide_pages >= 1,
                "{} must have at least one task-oriented guide",
                product.id
            );
            let version_history = product
                .pages
                .iter()
                .find(|page| page.id == "releases")
                .wrap_err_with(|| format!("{} has no version-history page", product.id))?;
            ensure!(
                version_history.title == "Version history"
                    && version_history.route == "version-history/"
                    && version_history.group == "Version history",
                "{} version history must use the shared title, route, and navigation group",
                product.id
            );
            if let Some(path) = &product.changelog {
                ensure!(
                    path.extension().and_then(|extension| extension.to_str()) == Some("typ"),
                    "{} changelog {} must be canonical Typst",
                    product.id,
                    path.display()
                );
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

    fn check_prose_sources(&self) -> Result<()> {
        ensure!(
            self.legacy_prose.schema == LEGACY_PROSE_SCHEMA_VERSION,
            "legacy prose schema {} does not match checker schema {}",
            self.legacy_prose.schema,
            LEGACY_PROSE_SCHEMA_VERSION
        );
        ensure!(
            self.legacy_prose.source.is_empty(),
            "docs/legacy-prose.toml is retired and cannot accept new sources"
        );
        let declared = self
            .legacy_prose
            .source
            .iter()
            .map(|source| {
                ensure!(
                    !source.is_absolute()
                        && normalize_repository_path(source).as_ref() == Some(source),
                    "legacy prose path {} must be a normalized repository-relative path",
                    source.display()
                );
                let extension = source
                    .extension()
                    .and_then(|extension| extension.to_str())
                    .map(str::to_ascii_lowercase);
                ensure!(
                    matches!(
                        extension.as_deref(),
                        Some(
                            "md" | "markdown" | "mdown" | "mkd" | "mdx" | "html" | "htm" | "xhtml"
                        )
                    ),
                    "legacy prose path {} must be Markdown or HTML",
                    source.display()
                );
                ensure!(
                    !matches!(
                        source.file_name().and_then(|name| name.to_str()),
                        Some("AGENTS.md" | "README.md")
                    ),
                    "{} is a compatibility exception and must not be legacy-listed",
                    source.display()
                );
                self.require_file(source)?;
                Ok(source.clone())
            })
            .collect::<Result<BTreeSet<_>>>()?;
        ensure!(
            declared.len() == self.legacy_prose.source.len(),
            "docs/legacy-prose.toml contains duplicate paths"
        );

        let mut actual = BTreeSet::new();
        for entry in WalkDir::new(&self.root).into_iter().filter_entry(|entry| {
            let top_level = entry
                .path()
                .strip_prefix(&self.root)
                .is_ok_and(|path| path.components().count() == 1);
            !entry.file_type().is_dir()
                || !top_level
                || !matches!(
                    entry.file_name().to_str(),
                    Some(".git" | ".jj" | ".direnv" | ".venv" | "node_modules" | "target")
                )
        }) {
            let entry = entry?;
            if !entry.file_type().is_file() && !entry.file_type().is_symlink() {
                continue;
            }
            let file_name = entry.file_name().to_string_lossy();
            ensure!(
                !file_name.eq_ignore_ascii_case("README.typ"),
                "{} is forbidden; README.md is the compatibility source and must not have a parallel README.typ",
                entry.path().strip_prefix(&self.root)?.display()
            );
            let original_extension = entry
                .path()
                .extension()
                .and_then(|extension| extension.to_str());
            let extension = entry
                .path()
                .extension()
                .and_then(|extension| extension.to_str())
                .map(str::to_ascii_lowercase);
            ensure!(
                extension.as_deref() != Some("typ") || original_extension == Some("typ"),
                "{} uses a non-canonical Typst extension; use lowercase .typ",
                entry.path().strip_prefix(&self.root)?.display()
            );
            if !matches!(
                extension.as_deref(),
                Some(
                    "md" | "markdown"
                        | "mdown"
                        | "mkd"
                        | "mdx"
                        | "html"
                        | "htm"
                        | "xhtml"
                        | "shtml"
                        | "rst"
                        | "rest"
                        | "adoc"
                        | "asciidoc"
                        | "org"
                )
            ) {
                continue;
            }
            let source = entry.path().strip_prefix(&self.root)?.to_path_buf();
            let mut typst = source.clone();
            typst.set_extension("typ");
            ensure!(
                !self.root.join(&typst).is_file(),
                "{} and {} are parallel editable documentation sources",
                source.display(),
                typst.display()
            );
            if !matches!(
                source.file_name().and_then(|name| name.to_str()),
                Some("AGENTS.md" | "README.md")
            ) {
                actual.insert(source);
            }
        }
        ensure!(
            actual == declared,
            "Markdown/HTML sources differ from the shrinking legacy inventory; new: {:?}; migrated or removed: {:?}",
            actual.difference(&declared).collect::<Vec<_>>(),
            declared.difference(&actual).collect::<Vec<_>>()
        );
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

        let (scope, selected) = if request.product == "all" {
            (
                BuildScope::FullSite,
                self.registry.product.iter().collect::<Vec<_>>(),
            )
        } else {
            (
                BuildScope::ProductPreview,
                vec![
                    self.registry
                        .product
                        .iter()
                        .find(|product| product.id == request.product)
                        .wrap_err_with(|| format!("unknown product {}", request.product))?,
                ],
            )
        };

        let requested_output = absolute_from(&self.root, &request.output);
        let resolved_output = ensure_safe_output(&self.root, &requested_output)?;
        let snapshot_build = scope == BuildScope::FullSite && tag.is_some();
        let isolated_work = if scope == BuildScope::ProductPreview || snapshot_build {
            let target = self.root.join("target");
            fs::create_dir_all(&target)?;
            let work = TempDirBuilder::new()
                .prefix(if scope == BuildScope::ProductPreview {
                    "alphal00p-product-preview-"
                } else {
                    "alphal00p-snapshot-build-"
                })
                .tempdir_in(target)?;
            let resolved_work = fs::canonicalize(work.path())?;
            ensure!(
                !resolved_work.starts_with(&resolved_output)
                    && !resolved_output.starts_with(&resolved_work),
                "documentation output overlaps isolated build workspace {}",
                work.path().display()
            );
            Some(work)
        } else {
            None
        };
        let output = isolated_work
            .as_ref()
            .map_or_else(|| requested_output.clone(), |work| work.path().join("site"));
        fs::create_dir_all(&output)
            .wrap_err_with(|| format!("failed to create {}", output.display()))?;
        clear_stale_staging(&output)?;

        let options = ProductBuildOptions {
            channel: request.channel,
            tag,
            output: &output,
            scope,
            include_typst: request.include_typst,
            include_rustdoc: request.include_rustdoc,
            rustdoc_target_root: request.rustdoc_target_root.as_deref(),
            dependency_output: request.dependency_output.as_deref(),
        };
        for product in &selected {
            self.build_product(product, &options)?;
        }

        match scope {
            BuildScope::FullSite => {
                self.write_developer_docs(
                    &output,
                    request.include_typst,
                    request.dependency_output.as_deref(),
                )?;
                self.write_portal(&output, request.channel, tag)?;
                self.write_federated_search_index(&output, request.channel, tag)?;
            }
            BuildScope::ProductPreview => {
                self.write_product_preview(&output, selected[0], request.channel, tag)?;
            }
        }
        clear_stale_staging(&output)?;
        self.validate_generated_links(&output, request.include_rustdoc, tag)?;
        if scope == BuildScope::ProductPreview {
            replace_generated_tree(&output, &requested_output)?;
        } else if snapshot_build {
            let tag = tag.expect("snapshot builds have a validated tag");
            for product in &selected {
                let relative = Path::new("products")
                    .join(&product.id)
                    .join("snapshots")
                    .join(tag);
                let generated = output.join(&relative);
                let existing = requested_output.join(&relative);
                if existing.exists() {
                    ensure!(
                        directories_equal(&generated, &existing)?,
                        "immutable snapshot differs from {}",
                        existing.display()
                    );
                }
            }
            // Validation normalizes generated Rustdoc links. Compare every existing
            // snapshot before carrying the other published product routes into the
            // validated candidate; replacing the tree then prunes stale global files.
            if requested_output.exists() {
                let metadata = fs::symlink_metadata(&requested_output)?;
                ensure!(
                    metadata.is_dir() && !metadata.file_type().is_symlink(),
                    "documentation destination is not a directory: {}",
                    requested_output.display()
                );
                for entry in WalkDir::new(&requested_output) {
                    let entry = entry?;
                    ensure!(
                        entry.file_type().is_dir() || entry.file_type().is_file(),
                        "unsupported documentation artifact {}",
                        entry.path().display()
                    );
                }

                for product in &selected {
                    let existing_product = requested_output.join("products").join(&product.id);
                    let generated_product = output.join("products").join(&product.id);
                    let existing_index = existing_product.join("index.html");
                    let generated_index = generated_product.join("index.html");
                    if existing_index.exists() && !generated_index.exists() {
                        let metadata = fs::symlink_metadata(&existing_index)?;
                        ensure!(
                            metadata.is_file() && !metadata.file_type().is_symlink(),
                            "unsupported documentation artifact {}",
                            existing_index.display()
                        );
                        fs::copy(&existing_index, &generated_index)?;
                    }

                    let existing_latest = existing_product.join("latest");
                    let generated_latest = generated_product.join("latest");
                    if existing_latest.exists() && !generated_latest.exists() {
                        copy_tree(&existing_latest, &generated_latest)?;
                    }

                    let existing_snapshots = existing_product.join("snapshots");
                    if existing_snapshots.exists() {
                        let metadata = fs::symlink_metadata(&existing_snapshots)?;
                        ensure!(
                            metadata.is_dir() && !metadata.file_type().is_symlink(),
                            "unsupported documentation artifact {}",
                            existing_snapshots.display()
                        );
                        for entry in fs::read_dir(existing_snapshots)? {
                            let entry = entry?;
                            if entry.file_name() == Path::new(tag).as_os_str() {
                                continue;
                            }
                            ensure!(
                                entry.file_type()?.is_dir(),
                                "unsupported documentation artifact {}",
                                entry.path().display()
                            );
                            copy_tree(
                                &entry.path(),
                                &generated_product.join("snapshots").join(entry.file_name()),
                            )?;
                        }
                    }
                }
            }
            replace_generated_tree(&output, &requested_output)?;
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
        let mut generated_pages = self.write_generated_reference(product, &site)?;

        if options.include_typst {
            self.render_typst(product, &metadata, &site, options.dependency_output)?;
        } else {
            self.write_fallback_page(product, &metadata, &site)?;
        }
        generated_pages.extend(self.write_python_reference(product, &metadata, &site)?);
        if options.include_rustdoc {
            self.build_rustdoc(product, &site, options.rustdoc_target_root)?;
        } else {
            self.write_rustdoc_placeholder(product, &site)?;
        }
        self.write_rust_reference_with_availability(product, &site, options.include_rustdoc)?;
        self.write_reference_hub(product, &site)?;
        self.decorate_site_pages_with_generated(
            product,
            &metadata,
            &site,
            options.scope,
            &generated_pages,
        )?;
        self.write_search_index(product, &site, options.include_rustdoc)?;

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
                "docs/assets/STIXTwoMath-Regular.woff2",
                "STIXTwoMath-Regular.woff2",
            ),
            ("docs/assets/STIX-Two-OFL.txt", "STIX-Two-OFL.txt"),
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
            let target = assets.join(name);
            match fs::remove_file(&target) {
                Ok(()) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => return Err(error.into()),
            }
            fs::copy(self.root.join(source), target)
                .wrap_err_with(|| format!("failed to copy documentation asset {source}"))?;
        }
        for asset in about_assets() {
            let target = assets.join(&asset);
            match fs::remove_file(&target) {
                Ok(()) => {}
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
                Err(error) => return Err(error.into()),
            }
            fs::copy(self.root.join("docs/assets").join(&asset), target)
                .wrap_err_with(|| format!("failed to copy About asset {asset}"))?;
        }
        Ok(())
    }

    fn write_portal_assets(&self, destination: &Path) -> Result<()> {
        let graphs = destination.join("assets/graphs");
        fs::create_dir_all(&graphs)?;
        for graph in portal_graph_assets() {
            fs::copy(
                self.root.join("docs/assets/graphs").join(&graph),
                graphs.join(&graph),
            )
            .wrap_err_with(|| format!("failed to copy portal graph {graph}"))?;
        }
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

    fn write_generated_reference(
        &self,
        product: &ProductConfig,
        site: &Path,
    ) -> Result<Vec<SitePage>> {
        let destination = site.join("reference/generated");
        let mut pages = Vec::new();
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
                    render_gammaloop_reference_index(&product.title, &reference),
                )?;
                pages.extend(write_gammaloop_reference_pages(
                    &product.title,
                    &reference,
                    site,
                )?);
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
        Ok(pages)
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
            "= Python API reference <supported-api-catalog>\n\
             The supported Python packages available in this version are listed below. The HTML manual provides a separate page for every public class, function, and constant; this printable edition keeps a compact inventory. Rust APIs use canonical Rustdoc in the HTML manual.\n",
        );
        for component in &product.python_components {
            // Everything beyond this boundary consumes the neutral catalog
            // artifact. In particular, no PyO3 inventory is ever linked into
            // the renderer process: each isolated exporter has already emitted
            // its checked stub, and the per-component adapter has serialized it.
            let path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&path)?)
                .wrap_err_with(|| format!("failed to merge {}", path.display()))?;
            rendered.push_str(&format!("\n== #raw(\"{}\")\n", typst_string(&component.id)));
            append_scope_typst(
                &catalog,
                &catalog.root,
                &mut Vec::new(),
                &metadata.git_commit,
                &mut rendered,
            );
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
        document_sections.push("#provenance");
        let document_body = document_sections.join(" #pagebreak() ");
        fs::write(
            &wrapper,
            format!(
                "#import \"/{source}\": manual\n\
                 {page_imports}\
                 #let catalog-reference = [{catalog_reference}]\n\
                 #let generated-reference = [{generated_reference}]\n\
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
        let rustdoc_css = self.root.join("docs/assets/rustdoc.css");
        let rustdoc_before_content = target.join("alphal00p-rustdoc-before-content.html");
        fs::write(
            &rustdoc_before_content,
            format!(
                "<nav class=\"alphal00p-rustdoc-bar\" aria-label=\"αLoop documentation\"><a href=\"{PUBLISHED_DOCS_ROOT}/\">αLoop docs</a><span aria-hidden=\"true\">/</span><span>{} · Rust API</span></nav>",
                escape_html(&product.title),
            ),
        )?;
        let mut rustdoc_flags = env::var("CARGO_ENCODED_RUSTDOCFLAGS")
            .ok()
            .filter(|flags| !flags.is_empty())
            .map(|flags| flags.split('\x1f').map(str::to_owned).collect::<Vec<_>>())
            .or_else(|| {
                env::var("RUSTDOCFLAGS").ok().map(|flags| {
                    flags
                        .split_whitespace()
                        .map(str::to_owned)
                        .collect::<Vec<_>>()
                })
            })
            .unwrap_or_default();
        ensure!(
            !rustdoc_flags.iter().any(|flag| {
                flag == "--extend-css"
                    || flag.starts_with("--extend-css=")
                    || flag.starts_with("-e") && !flag.starts_with("--")
            }),
            "the documentation builder supplies --extend-css; remove it from inherited Rustdoc flags"
        );
        rustdoc_flags.extend(STRICT_RUSTDOC_FLAGS.split_whitespace().map(str::to_owned));
        rustdoc_flags.extend([
            "--extend-css".to_owned(),
            rustdoc_css.to_string_lossy().into_owned(),
            "--html-before-content".to_owned(),
            rustdoc_before_content.to_string_lossy().into_owned(),
        ]);
        let rustdoc_flags = rustdoc_flags.join("\x1f");
        for component in &product.rust_components {
            let mut command = Command::new("cargo");
            command
                .current_dir(&self.root)
                .env("CARGO_TARGET_DIR", &target)
                .env_remove("RUSTDOCFLAGS")
                .env("CARGO_ENCODED_RUSTDOCFLAGS", &rustdoc_flags)
                .args(["doc", "--locked", "--no-deps", "--no-default-features"]);
            if let Some(profile) = env::var_os("ALPHAL00P_DOCS_CARGO_PROFILE") {
                command.arg("--profile").arg(profile);
            }
            command.args(["-p", &component.package]);
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
    ) -> Result<Vec<SitePage>> {
        let destination = site.join("reference/python");
        fs::create_dir_all(&destination)?;
        let mut cards = String::new();
        let mut pages = Vec::new();
        for component in &product.python_components {
            let stub = self.python_stub_source(component);
            let stub_name = format!("{}.pyi", component.id);
            if let Some(path) = stub {
                let text = fs::read_to_string(&path)
                    .wrap_err_with(|| format!("failed to read {}", path.display()))?;
                fs::write(destination.join(&stub_name), &text)?;
            } else {
                let module = component.module.as_deref().unwrap_or(&component.package);
                let text = format!("# Python interface for {module}\n");
                fs::write(destination.join(&stub_name), &text)?;
            }
            let catalog_path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&catalog_path)?)
                .wrap_err_with(|| format!("failed to parse {}", catalog_path.display()))?;
            catalog.validate()?;
            let module = python_display_module(component);
            let routes = python_item_routes(&catalog);
            ensure!(
                routes.values().collect::<BTreeSet<_>>().len() == routes.len(),
                "duplicate Python reference route for {}",
                component.id
            );
            let export_count = routes.len();
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
                    &format!("{} Python API", module),
                    &render_python_catalog_for_module(
                        &catalog,
                        module,
                        &stub_name,
                        &metadata.git_commit,
                    ),
                ),
            )?;
            let supported_items = catalog_surface_items(&catalog.root, true);
            for (index, (path, item)) in supported_items.iter().enumerate() {
                let anchor = api_item_anchor(path, item);
                let route = routes
                    .get(&anchor)
                    .with_context(|| format!("missing Python route for {anchor}"))?;
                let item_page = site.join(route);
                fs::create_dir_all(&item_page)?;
                fs::write(
                    item_page.join("index.html"),
                    reference_page(
                        &product.title,
                        &item.name,
                        &render_python_item_page(
                            &catalog,
                            module,
                            (path, item),
                            &stub_name,
                            &metadata.git_commit,
                            &routes,
                            (
                                index.checked_sub(1).and_then(|previous| {
                                    supported_items.get(previous).map(|(path, item)| {
                                        (
                                            item.name.as_str(),
                                            routes[&api_item_anchor(path, item)].as_str(),
                                        )
                                    })
                                }),
                                supported_items.get(index + 1).map(|(path, item)| {
                                    (
                                        item.name.as_str(),
                                        routes[&api_item_anchor(path, item)].as_str(),
                                    )
                                }),
                            ),
                        ),
                    ),
                )?;
                pages.push(SitePage::new(
                    route.clone(),
                    format!("{}.{}", module, item.name),
                    "Python API",
                ));
            }
        }
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Python API",
                &format!(
                    "<p>Choose a Python module to browse classes, functions, signatures, parameters, examples, and feature requirements. Type-checker <code>.pyi</code> files are also available for download.</p><nav class=\"reference-guide-links\" aria-label=\"Related Python documentation\"><a href=\"reference/interfaces/\">Python interface guide</a></nav><div class=\"api-component-grid\">{cards}</div>"
                ),
            ),
        )?;
        Ok(pages)
    }

    fn write_rust_reference_with_availability(
        &self,
        product: &ProductConfig,
        site: &Path,
        rustdoc_available: bool,
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
            let crate_name = component.package.replace('-', "_");
            let features = if component.features.is_empty() {
                "no optional features".to_owned()
            } else {
                component
                    .features
                    .iter()
                    .map(|feature| format!("<code>{}</code>", escape_html(feature)))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            let rustdoc_link = if rustdoc_available {
                format!(
                    "<a class=\"portal-text-link\" href=\"reference/rust/{crate_name}/\">Open canonical Rustdoc <span aria-hidden=\"true\">→</span></a>"
                )
            } else {
                "<p class=\"reference-missing\">Rustdoc was skipped for this preview build.</p>"
                    .to_owned()
            };
            cards.push_str(&format!(
                "<article class=\"api-component-card rustdoc-card\"><h2><code>{}</code></h2><p>Version {} · {}</p><p>{}</p><p><strong>Documented features:</strong> {features}</p>{rustdoc_link}</article>",
                escape_html(&component.package),
                escape_html(&self.component_version(component)?),
                escape_html(&component.kind),
                escape_html(&catalog.root.summary.clone().unwrap_or_else(|| {
                    format!("Rust API documentation for {}.", component.id)
                })),
            ));
            let component_page = supported.join(&component.id);
            fs::create_dir_all(&component_page)?;
            fs::write(
                component_page.join("index.html"),
                generated_reference_redirect(
                    &product.title,
                    &if rustdoc_available {
                        format!("../../{crate_name}/")
                    } else {
                        "../../".to_owned()
                    },
                    &format!("{} Rustdoc", component.package),
                ),
            )?;
        }
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Rust API",
                &format!(
                    "<p>Rustdoc is the canonical Rust API reference. It already provides the familiar crate hierarchy, exact signatures, trait implementations, feature-gated items, source links, and native symbol search without translating those structures into a second UI.</p><p>Use the authored interface guide for product-level boundaries and choose a crate below when you need item-level detail. Presence in Rustdoc describes the compiled public surface; it is not by itself a compatibility promise. The structured catalogs remain build artifacts for coverage and consistency checks, not a competing public reference.</p><nav class=\"reference-guide-links\" aria-label=\"Related Rust documentation\"><a href=\"reference/interfaces/\">Understand the Rust interfaces</a><a href=\"quickstart/rust/\">{}</a></nav><div class=\"api-component-grid\">{cards}</div>",
                    escape_html(&interface_guide_title(&product.title, "rust")),
                ),
            ),
        )?;
        Ok(())
    }

    fn write_reference_hub(&self, product: &ProductConfig, site: &Path) -> Result<()> {
        let mut cards = product
            .pages
            .iter()
            .find(|page| page.route == "reference/interfaces/")
            .map(|page| {
                format!(
                    "<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Choose an interface</p><h2><a href=\"{}\">{}</a></h2><p>{}</p></article>",
                    escape_html(&page.route),
                    escape_html(&page.title),
                    escape_html(&page.summary),
                )
            })
            .unwrap_or_default();
        cards.push_str(&format!(
            "<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Canonical Rustdoc</p><h2><a href=\"reference/rust/\">Rust API</a></h2><p>{} crates in their native Rust documentation interface, with authored guidance for choosing the right boundary.</p></article><article class=\"reference-hub-card\"><p class=\"portal-kicker\">Python-native reference</p><h2><a href=\"reference/python/\">Python API</a></h2><p>{} modules with one flat page per supported class, function, or constant; stubs remain downloadable for type checkers.</p></article>",
            product.rust_components.len(),
            product.python_components.len(),
        ));
        let has_typst = product
            .pages
            .iter()
            .any(|page| page.route == "reference/typst/");
        if has_typst {
            cards.push_str("<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Typst-native reference</p><h2><a href=\"reference/typst/\">Typst API</a></h2><p>Package exports grouped by graph construction, layout, drawing, domain templates, and subgraph operations.</p></article>");
        }
        if let Some((route, title)) = supplemental_reference(&product.id) {
            cards.push_str(&format!(
                "<article class=\"reference-hub-card\"><p class=\"portal-kicker\">Generated from the implementation</p><h2><a href=\"{route}\">{}</a></h2><p>Version-specific reference data, also available as JSON for tooling.</p></article>",
                escape_html(title),
            ));
        }
        let destination = site.join("reference");
        fs::create_dir_all(&destination)?;
        let mut interface_shapes = vec!["native Rustdoc for crates"];
        if !product.python_components.is_empty() {
            interface_shapes.push("class, function, and constant pages for Python");
        }
        if has_typst {
            interface_shapes.push("focused function and module pages for Typst");
        }
        if product.id == "gammaloop" {
            interface_shapes.push("manpage-style command pages for the CLI");
        }
        fs::write(
            destination.join("index.html"),
            reference_page(
                &product.title,
                "Reference",
                &format!(
                    "<p>Choose the reference shaped for the interface you use: {}. Authored guides connect those exact interfaces to complete workflows.</p><div class=\"reference-hub-grid\">{cards}</div>",
                    interface_shapes.join("; "),
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

    fn decorate_site_pages_with_generated(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
        scope: BuildScope,
        generated_pages: &[SitePage],
    ) -> Result<()> {
        let mut pages = product_site_pages(product);
        pages.retain(|page| site.join(&page.route).join("index.html").is_file());

        for (index, page) in pages.iter().enumerate() {
            self.decorate_html_page(
                product,
                metadata,
                site,
                page,
                (
                    index
                        .checked_sub(1)
                        .and_then(|previous| pages.get(previous)),
                    pages.get(index + 1),
                ),
                (scope, generated_pages),
            )?;
            if let Some(authored) = product
                .pages
                .iter()
                .find(|authored| authored.route == page.route)
            {
                for alias in &authored.aliases {
                    let directory = site.join(alias);
                    fs::create_dir_all(&directory)?;
                    let target = format!("{}{}", page_root_prefix(alias), authored.route);
                    fs::write(
                        directory.join("index.html"),
                        generated_reference_redirect(&product.title, &target, &authored.title),
                    )?;
                }
            }
        }
        for page in generated_pages {
            if site.join(&page.route).join("index.html").is_file() {
                self.decorate_html_page(
                    product,
                    metadata,
                    site,
                    page,
                    (None, None),
                    (scope, generated_pages),
                )?;
            }
        }
        Ok(())
    }

    fn decorate_html_page(
        &self,
        product: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        site: &Path,
        page: &SitePage,
        siblings: (Option<&SitePage>, Option<&SitePage>),
        context: (BuildScope, &[SitePage]),
    ) -> Result<()> {
        let (previous, next) = siblings;
        let (scope, generated_pages) = context;
        let path = site.join(&page.route).join("index.html");
        let source = fs::read_to_string(&path)?;
        let mut body = extract_html_body(&source)?.to_owned();
        body = body.replace("<nav><a href=\"../../\">Documentation home</a></nav>", "");
        if page.route.is_empty() {
            body = format!("{}{body}", product_hero(product));
        } else {
            promote_heading_hierarchy(&mut body);
        }
        let docs_root = page_root_prefix(&page.route);
        body = rewrite_page_links(&body, &docs_root)?;
        if scope == BuildScope::ProductPreview {
            body = self.rewrite_product_preview_links(&body, metadata, page)?;
        }
        for (color, class) in [
            ("#4b69c6", "syntax-call"),
            ("#d73948", "syntax-keyword"),
            ("#b60157", "syntax-number"),
            ("#198810", "syntax-string"),
            ("#6b6b6f", "syntax-raw"),
            ("#8b41b1", "syntax-property"),
            ("#1d6c76", "syntax-label"),
        ] {
            body = body.replace(
                &format!("<span style=\"color: {color}\">"),
                &format!("<span class=\"{class}\">"),
            );
        }
        let invalid_parameter =
            Regex::new(r#"(?s)<p>(?P<heading><h5\b[^>]*>.*?</h5>)(?P<details>.*?)</p>"#)?;
        body = invalid_parameter
            .replace_all(&body, |captures: &regex::Captures<'_>| {
                let heading = captures
                    .name("heading")
                    .expect("Tidy parameter heading")
                    .as_str();
                let details = captures
                    .name("details")
                    .expect("Tidy parameter details")
                    .as_str();
                if details.trim().is_empty() {
                    heading.to_owned()
                } else {
                    format!("{heading}<p>{details}</p>")
                }
            })
            .into_owned();
        body = inject_heading_ids(&body)?.0;
        if let Some((reference_heading, reference_scope)) = match page.route.as_str() {
            "reference/typst/graph/" => Some(("graph", Some("graph"))),
            "reference/typst/layout/" | "reference/typst/drawing/" => Some(("reference", None)),
            "reference/typst/subgraph/" => Some(("subgraph", Some("subgraph"))),
            _ => None,
        } {
            let marker = format!("<h2 id=\"{reference_heading}\">");
            if let Some(reference_start) = body.find(&marker) {
                if let Some(reference_scope) = reference_scope {
                    let symbol_heading = Regex::new(r#"<h3 id=\"[^\"]+\">(?P<title>[^<]+)</h3>"#)?;
                    let qualified = symbol_heading
                        .replace_all(
                            &body[reference_start..],
                            |captures: &regex::Captures<'_>| {
                                let title = captures
                                    .name("title")
                                    .expect("Tidy symbol heading")
                                    .as_str();
                                let symbol = decode_html_text(title);
                                format!(
                                    "<h3 id=\"{reference_scope}.{}\">{title}</h3>",
                                    escape_html(&symbol),
                                )
                            },
                        )
                        .into_owned();
                    body.replace_range(reference_start.., &qualified);
                }
                let wrapper = format!("<div>{marker}");
                body = body.replacen(
                    &wrapper,
                    &format!("<div class=\"typst-api-module\">{marker}"),
                    1,
                );
            }
        }
        if page.route == "reference/typst/layout/" {
            body = body
                .replace("<h3 id=\"layouts-layout\">", "<h3 id=\"layouts.layout\">")
                .replace(
                    "<h3 id=\"layouts-sequence\">",
                    "<h3 id=\"layouts.sequence\">",
                );
        }
        let (body, headings) = inject_heading_ids(&body)?;
        let dense_reference = matches!(page.group.as_str(), "Python API" | "CLI reference");
        let toc_links = headings
            .iter()
            .filter(|heading| heading.level == 2 || heading.level == 3 && !dense_reference)
            .map(|heading| {
                format!(
                    "<a class=\"toc-link\" data-level=\"{}\" href=\"#{}\">{}</a>",
                    heading.level,
                    escape_html(&heading.id),
                    escape_html(&heading.title)
                )
            })
            .collect::<Vec<_>>();
        let toc = if toc_links.len() >= 3 || dense_reference && !toc_links.is_empty() {
            toc_links.concat()
        } else {
            String::new()
        };
        let inline_toc = if toc.is_empty() {
            String::new()
        } else {
            format!(
                "<details class=\"docs-inline-toc\"><summary>On this page</summary><nav aria-label=\"On this page\">{toc}</nav></details>"
            )
        };
        let sidebar =
            self.site_sidebar(product, metadata, page, &docs_root, scope, generated_pages);
        let product_options = self.product_options(product, metadata, &docs_root, scope);
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
        let site_home = if scope == BuildScope::FullSite {
            &portal
        } else {
            product_root
        };
        let version = match metadata.channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => {
                format!("snapshot {}", metadata.snapshot_tag.unwrap_or("unknown"))
            }
        };
        let search_root = if scope == BuildScope::FullSite {
            portal.as_str()
        } else {
            product_root
        };
        let source_path = product
            .pages
            .iter()
            .find(|candidate| candidate.route == page.route)
            .map(|candidate| candidate.source.clone())
            .or_else(|| {
                product
                    .python_components
                    .iter()
                    .find(|component| {
                        page.route
                            .starts_with(&format!("reference/python/{}/", component.id))
                    })
                    .map(|component| PathBuf::from(format!("docs/api/python/{}.pyi", component.id)))
            })
            .unwrap_or_else(|| match (product.id.as_str(), page.route.as_str()) {
                ("gammaloop", route) if route.starts_with("reference/cli/settings/runtime/") => {
                    PathBuf::from("crates/gammalooprs/src/settings/mod.rs")
                }
                ("gammaloop", route) if route.starts_with("reference/cli/settings/") => {
                    PathBuf::from("crates/gammaloop-api/src/lib.rs")
                }
                ("gammaloop", route) if route.starts_with("reference/cli/") => {
                    PathBuf::from("crates/gammaloop-api/src/commands/mod.rs")
                }
                ("vakint", "reference/topologies/") => {
                    PathBuf::from("crates/vakint/src/topologies.rs")
                }
                _ => PathBuf::from("docs/products/registry.toml"),
            });
        let source_url = format!(
            "https://github.com/alphal00p/gammaloop/blob/{}/{}",
            metadata.git_commit,
            source_path.to_string_lossy().replace('\\', "/")
        );
        let documented_route = format!("{}{}", metadata.route, page.route);
        let issue_url = documentation_issue_url(
            &format!("{}: {}", product.title, page.title),
            &documented_route,
            &metadata.git_commit,
        );
        let source_label = if product.python_components.iter().any(|component| {
            page.route
                .starts_with(&format!("reference/python/{}/", component.id))
        }) {
            "View generated type stub at this revision"
        } else if product.id == "gammaloop" && page.route.starts_with("reference/cli/") {
            "View parser or schema source at this revision"
        } else {
            "View page source at this revision"
        };
        let feedback = format!(
            "<nav class=\"page-feedback\" aria-label=\"Documentation feedback\"><a href=\"{}\">{} <span aria-hidden=\"true\">↗</span></a><a href=\"{}\">Report a documentation issue <span aria-hidden=\"true\">↗</span></a></nav>",
            escape_html(&source_url),
            escape_html(source_label),
            escape_html(&issue_url),
        );
        let toc_markup = if toc.is_empty() {
            String::new()
        } else {
            format!(
                "<aside class=\"docs-toc\" aria-label=\"On this page\"><p class=\"toc-title\">On this page</p>{toc}</aside>"
            )
        };
        let reference_language = if page.route.starts_with("reference/python/") {
            "python"
        } else if page.route.starts_with("reference/typst/") {
            "typst"
        } else if page.route.starts_with("reference/rust/") {
            "rust"
        } else if page.route.starts_with("reference/cli/") {
            "cli"
        } else {
            ""
        };
        let breadcrumb_page = if page.group == page.title {
            String::new()
        } else {
            format!(" / {}", escape_html(&page.title))
        };
        let favicon = favicon_links(&format!("{docs_root}assets/"));
        let html = format!(
            "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\">{favicon}<title>{} · {}</title><link rel=\"stylesheet\" href=\"{}assets/site.css\"><script defer src=\"{}assets/site.js\"></script></head><body data-reference-language=\"{}\" data-has-toc=\"{}\" data-docs-root=\"{}\" data-search-index=\"{}search-index.json\" data-search-root=\"{}\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"site-header\"><button class=\"header-button menu-button\" type=\"button\" data-menu-toggle aria-label=\"Open navigation\" aria-controls=\"docs-sidebar\" aria-expanded=\"false\">☰</button><a class=\"site-brand\" href=\"{}\"><span class=\"site-brand-mark\">α</span><span class=\"site-brand-name\">αLoop docs</span></a><div class=\"site-header-tools\"><select class=\"product-select\" data-product-select aria-label=\"Select project\">{}</select><button class=\"header-button\" type=\"button\" data-search-open aria-label=\"Search all documentation\"><span class=\"header-search-label\">Search</span> <span class=\"header-button-label\">⌘K</span><span class=\"header-search-icon\" aria-hidden=\"true\">⌕</span></button><button class=\"header-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\">◐</button></div></header><div class=\"docs-shell\">{sidebar}<main class=\"docs-main\" id=\"main-content\"><nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"{}\">{}</a> / {}{}</nav>{inline_toc}<article class=\"docs-article\">{body}</article>{page_navigation}{feedback}<footer class=\"page-footer\">{} · <a href=\"{}manual.pdf\">Complete PDF manual</a> · documented revision <code>{}</code></footer></main>{toc_markup}</div><button class=\"sidebar-backdrop\" type=\"button\" data-sidebar-backdrop aria-label=\"Close navigation\"></button><dialog class=\"search-dialog\" data-search-dialog><form class=\"search-form\" method=\"dialog\"><input class=\"search-input\" type=\"search\" data-search-input placeholder=\"Search all projects and developer notes\" aria-label=\"Search all documentation\"><button class=\"header-button\" value=\"close\">Close</button></form><ul class=\"search-results\" data-search-results aria-live=\"polite\"></ul></dialog></body></html>",
            escape_html(&format!("{} — {}", product.tagline, page.title)),
            escape_html(&product.title),
            escape_html(&page.title),
            escape_html(&docs_root),
            escape_html(&docs_root),
            escape_html(reference_language),
            if toc.is_empty() { "false" } else { "true" },
            escape_html(&docs_root),
            escape_html(search_root),
            escape_html(search_root),
            escape_html(site_home),
            product_options,
            escape_html(product_root),
            escape_html(&product.title),
            escape_html(&page.group),
            breadcrumb_page,
            escape_html(&version),
            escape_html(&docs_root),
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
        scope: BuildScope,
        generated_pages: &[SitePage],
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
        if !groups.contains_key("Reference") {
            group_order.push("Reference".to_owned());
        }
        let reference = groups.entry("Reference".to_owned()).or_default();
        reference.push(("reference/".to_owned(), "Reference overview".to_owned()));
        reference.push(("reference/python/".to_owned(), "Python API".to_owned()));
        reference.push(("reference/rust/".to_owned(), "Rust API".to_owned()));
        if let Some((route, title)) = supplemental_reference(&product.id) {
            reference.push((route.to_owned(), title.to_owned()));
        }

        let python_navigation = product
            .python_components
            .iter()
            .find_map(|component| {
                let module_route = format!("reference/python/{}/", component.id);
                current
                    .route
                    .starts_with(&module_route)
                    .then_some((component, module_route))
            })
            .map_or_else(String::new, |(component, module_route)| {
                let module_symbols = generated_pages
                    .iter()
                    .filter(|page| {
                        page.group == "Python API" && page.route.starts_with(&module_route)
                    })
                    .map(|page| {
                        let label = page.route[module_route.len()..]
                            .trim_end_matches('/')
                            .replace('/', ".");
                        format!(
                            "<a class=\"sidebar-link sidebar-python-symbol\" href=\"{}{}\"{}><code>{}</code></a>",
                            escape_html(docs_root),
                            escape_html(&page.route),
                            if page.route == current.route {
                                " aria-current=\"page\""
                            } else {
                                ""
                            },
                            escape_html(&label),
                        )
                    })
                    .collect::<String>();
                let module = python_display_module(component);
                format!(
                    "<nav class=\"sidebar-group sidebar-python-group\" aria-label=\"{} Python symbols\"><p class=\"sidebar-group-title\">{} symbols</p><a class=\"sidebar-link\" href=\"{}{}\"{}>Module overview</a>{module_symbols}</nav>",
                    escape_html(module),
                    escape_html(module),
                    escape_html(docs_root),
                    escape_html(&module_route),
                    if current.route == module_route {
                        " aria-current=\"page\""
                    } else {
                        ""
                    },
                )
            });
        let cli_navigation = if !current.route.starts_with("reference/cli/") {
            String::new()
        } else if current.route.starts_with("reference/cli/settings/") {
            let current_namespace = current
                .title
                .strip_suffix(" settings")
                .filter(|namespace| *namespace != "Settings namespaces");
            let namespace_links = generated_pages
                .iter()
                .filter(|page| {
                    if !page.route.starts_with("reference/cli/settings/")
                        || page.route == "reference/cli/settings/"
                    {
                        return false;
                    }
                    let Some(namespace) = page.title.strip_suffix(" settings") else {
                        return false;
                    };
                    let depth = namespace.matches('.').count();
                    depth == 0
                        || current_namespace.is_some_and(|current| {
                            namespace == current
                                || current.starts_with(&format!("{namespace}."))
                                || namespace.rsplit_once('.').is_some_and(|(parent, _)| {
                                    current == parent || current.starts_with(&format!("{parent}."))
                                })
                        })
                })
                .map(|page| {
                    let namespace = page.title.strip_suffix(" settings").unwrap_or(&page.title);
                    let depth = namespace.matches('.').count().min(4);
                    let ancestor = current_namespace
                        .is_some_and(|current| current.starts_with(&format!("{namespace}.")));
                    format!(
                        "<a class=\"sidebar-link sidebar-context-link{}\" data-depth=\"{depth}\" href=\"{}{}\" aria-label=\"{} settings namespace\"{}><code>{}</code></a>",
                        if ancestor { " is-ancestor" } else { "" },
                        escape_html(docs_root),
                        escape_html(&page.route),
                        escape_html(namespace),
                        if page.route == current.route {
                            " aria-current=\"page\""
                        } else {
                            ""
                        },
                        escape_html(namespace.rsplit('.').next().unwrap_or(namespace)),
                    )
                })
                .collect::<String>();
            format!(
                "<nav class=\"sidebar-group sidebar-context-group\" aria-label=\"GammaLoop settings namespaces\"><p class=\"sidebar-group-title\">Settings namespaces</p><a class=\"sidebar-link\" href=\"{}reference/cli/\">All commands</a><a class=\"sidebar-link\" href=\"{}reference/cli/settings/\"{}>All namespaces</a>{namespace_links}</nav>",
                escape_html(docs_root),
                escape_html(docs_root),
                if current.route == "reference/cli/settings/" {
                    " aria-current=\"page\""
                } else {
                    ""
                },
            )
        } else {
            let current_command = current
                .route
                .starts_with("reference/cli/commands/")
                .then_some(current.title.as_str());
            let command_links = generated_pages
                .iter()
                .filter(|page| page.route.starts_with("reference/cli/commands/"))
                .filter(|page| {
                    let depth = page.title.split_whitespace().count().saturating_sub(1);
                    depth <= 1
                        || current_command.is_some_and(|current| {
                            page.title == current
                                || current.starts_with(&format!("{} ", page.title))
                                || page.title.rsplit_once(' ').is_some_and(|(parent, _)| {
                                    current == parent || current.starts_with(&format!("{parent} "))
                                })
                        })
                })
                .map(|page| {
                    let depth = page.title.split_whitespace().count().saturating_sub(1).min(4);
                    let ancestor = current_command
                        .is_some_and(|current| current.starts_with(&format!("{} ", page.title)));
                    format!(
                        "<a class=\"sidebar-link sidebar-context-link{}\" data-depth=\"{depth}\" href=\"{}{}\" aria-label=\"{} command\"{}><code>{}</code></a>",
                        if ancestor { " is-ancestor" } else { "" },
                        escape_html(docs_root),
                        escape_html(&page.route),
                        escape_html(&page.title),
                        if page.route == current.route {
                            " aria-current=\"page\""
                        } else {
                            ""
                        },
                        escape_html(page.title.rsplit(' ').next().unwrap_or(&page.title)),
                    )
                })
                .collect::<String>();
            format!(
                "<nav class=\"sidebar-group sidebar-context-group\" aria-label=\"GammaLoop command tree\"><p class=\"sidebar-group-title\">Command tree</p><a class=\"sidebar-link\" href=\"{}reference/cli/\">All commands</a><a class=\"sidebar-link\" href=\"{}reference/cli/settings/\">Settings namespaces</a>{command_links}</nav>",
                escape_html(docs_root),
                escape_html(docs_root),
            )
        };
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
                        if route == &current.route {
                            " aria-current=\"page\""
                        } else if matches!(
                            route.as_str(),
                            "reference/python/" | "reference/rust/" | "reference/cli/"
                        ) && current.route.starts_with(route)
                        {
                            " aria-current=\"location\""
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
            if group == "Reference" {
                navigation.push_str(&python_navigation);
                navigation.push_str(&cli_navigation);
            }
        }
        if scope == BuildScope::FullSite && metadata.channel == BuildChannel::Latest {
            let portal = format!("{docs_root}../../../");
            navigation.push_str(&format!(
                "<section class=\"sidebar-group sidebar-developer-group\"><p class=\"sidebar-group-title\">For developers</p><a class=\"sidebar-link\" href=\"{}developers/\">Architecture &amp; engineering notes</a></section>",
                escape_html(&portal),
            ));
        }
        let citation_href = match scope {
            BuildScope::FullSite => {
                let portal = match metadata.channel {
                    BuildChannel::Latest => format!("{docs_root}../../../"),
                    BuildChannel::Snapshot => format!("{docs_root}../../../../"),
                };
                format!("{portal}citations/#{}", product.id)
            }
            BuildScope::ProductPreview => {
                format!("{PUBLISHED_DOCS_ROOT}/citations/#{}", product.id)
            }
        };
        let version = match metadata.channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => metadata.snapshot_tag.unwrap_or("unknown").to_owned(),
        };
        format!(
            "<aside class=\"docs-sidebar\" id=\"docs-sidebar\" aria-label=\"{} manual\">{navigation}<p class=\"sidebar-meta\"><strong>{}</strong><br><code>{}</code><br><a href=\"{}manual.pdf\">Download PDF</a><br><a href=\"{}\">Cite {}</a></p></aside>",
            escape_html(&product.title),
            escape_html(&product.title),
            escape_html(&version),
            escape_html(docs_root),
            escape_html(&citation_href),
            escape_html(&product.title),
        )
    }

    fn product_options(
        &self,
        current: &ProductConfig,
        metadata: &SnapshotMetadata<'_>,
        docs_root: &str,
        scope: BuildScope,
    ) -> String {
        self.registry
            .product
            .iter()
            .filter(|product| scope == BuildScope::FullSite || product.id == current.id)
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

    fn rewrite_product_preview_links(
        &self,
        body: &str,
        metadata: &SnapshotMetadata<'_>,
        page: &SitePage,
    ) -> Result<String> {
        let attribute = Regex::new(r#"(?P<name>href|src)=\"(?P<target>[^\"]+)\""#)?;
        let page_root = PathBuf::from(metadata.route.trim_start_matches('/')).join(&page.route);
        Ok(attribute
            .replace_all(body, |captures: &regex::Captures<'_>| {
                let name = captures.name("name").expect("attribute name").as_str();
                let target = captures.name("target").expect("attribute target").as_str();
                let local = !target.is_empty()
                    && !target.starts_with('#')
                    && !target.starts_with('/')
                    && !target.starts_with("//")
                    && !target.split(':').next().is_some_and(|scheme| {
                        matches!(scheme, "http" | "https" | "mailto" | "data")
                    });
                if !local {
                    return captures[0].to_owned();
                }

                let suffix_start = target.find(['?', '#']).unwrap_or(target.len());
                let (path, suffix) = target.split_at(suffix_start);
                let Some(resolved) = normalize_repository_path(&page_root.join(path)) else {
                    return captures[0].to_owned();
                };
                let segments = resolved
                    .iter()
                    .filter_map(|segment| segment.to_str())
                    .collect::<Vec<_>>();
                if segments.len() == 3
                    && segments.first() == Some(&"developers")
                    && segments.get(1) == Some(&"architecture")
                {
                    let developer_id = segments[2];
                    if let Some(note) = self
                        .developers
                        .section
                        .iter()
                        .flat_map(|section| &section.note)
                        .find(|note| note.id == developer_id)
                    {
                        let target = format!(
                            "https://github.com/alphal00p/gammaloop/blob/{}/{}{}",
                            metadata.git_commit,
                            note.source.to_string_lossy().replace('\\', "/"),
                            suffix,
                        );
                        return format!("{name}=\"{}\"", escape_html(&target));
                    }
                    return captures[0].to_owned();
                }
                let Some(product_id) = segments
                    .first()
                    .filter(|segment| **segment == "products")
                    .and_then(|_| segments.get(1))
                else {
                    return captures[0].to_owned();
                };
                if *product_id == metadata.product
                    || !self
                        .registry
                        .product
                        .iter()
                        .any(|product| product.id == *product_id)
                {
                    return captures[0].to_owned();
                }

                let resolved = resolved.to_string_lossy().replace('\\', "/");
                let trailing_slash = if path.ends_with('/') && !resolved.ends_with('/') {
                    "/"
                } else {
                    ""
                };
                format!("{name}=\"{PUBLISHED_DOCS_ROOT}/{resolved}{trailing_slash}{suffix}\"")
            })
            .into_owned())
    }

    fn write_search_index(
        &self,
        product: &ProductConfig,
        site: &Path,
        include_rustdoc: bool,
    ) -> Result<()> {
        let mut entries = vec![SearchEntry {
            title: product.title.clone(),
            summary: product.tagline.clone(),
            href: "index.html".to_owned(),
            kind: "product".to_owned(),
            text: String::new(),
        }];
        for page in &product.pages {
            let page_href = if page.route.is_empty() {
                "index.html".to_owned()
            } else {
                page.route.clone()
            };
            let rendered = site.join(&page.route).join("index.html");
            let html = rendered
                .is_file()
                .then(|| {
                    fs::read_to_string(&rendered)
                        .wrap_err_with(|| format!("failed to read {}", rendered.display()))
                })
                .transpose()?;
            entries.push(SearchEntry {
                title: page.title.clone(),
                summary: page.summary.clone(),
                href: page_href.clone(),
                kind: page.group.to_lowercase(),
                text: html
                    .as_deref()
                    .map(rendered_article_search_text)
                    .transpose()?
                    .unwrap_or_default(),
            });
            if let Some(html) = html {
                append_rendered_page_search(product, page, &page_href, &html, &mut entries)?;
            }
        }
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
            if language == ApiLanguage::Rust && !include_rustdoc {
                continue;
            }
            let display_name = python_display_module(component);
            match language {
                ApiLanguage::Rust => entries.push(SearchEntry {
                    title: format!("{} Rustdoc", component.package),
                    summary: format!(
                        "{} crate · version {}",
                        component.kind,
                        self.component_version(component)?
                    ),
                    href: format!("reference/rust/{}/", component.package.replace('-', "_")),
                    kind: "rust-api".to_owned(),
                    text: component.features.join(" "),
                }),
                ApiLanguage::Python => entries.push(SearchEntry {
                    title: format!("{display_name} Python API"),
                    summary: format!(
                        "Supported Python module · version {}",
                        self.component_version(component)?
                    ),
                    href: format!("reference/python/{}/", component.id),
                    kind: "python-api".to_owned(),
                    text: String::new(),
                }),
                _ => {}
            }
            let path = site.join("catalogs").join(format!("{}.json", component.id));
            let catalog = serde_json::from_slice::<DocCatalog>(&fs::read(&path)?)
                .wrap_err_with(|| format!("failed to parse {}", path.display()))?;
            let catalog_display_name = if language == ApiLanguage::Python {
                display_name
            } else {
                &catalog.component.id
            };
            append_catalog_search_with_component_name(
                &catalog,
                &catalog.root,
                catalog_display_name,
                &mut entries,
            );
        }
        match product.id.as_str() {
            "gammaloop" => {
                let (reference, _) = self.generated_references()?;
                let namespaces = cli_setting_namespaces(&reference);
                entries.extend(
                    reference
                        .commands
                        .into_iter()
                        .filter(|command| !command.hidden && !command.generated_help)
                        .map(|command| {
                            let text = command
                                .arguments
                                .iter()
                                .filter(|argument| !argument.hidden)
                                .map(|argument| {
                                    format!(
                                        "{} {} {} {}",
                                        argument.id,
                                        argument
                                            .long
                                            .as_ref()
                                            .map(|name| format!("--{name}"))
                                            .unwrap_or_default(),
                                        argument
                                            .short
                                            .map(|name| format!("-{name}"))
                                            .unwrap_or_default(),
                                        argument.help,
                                    )
                                })
                                .chain(command.aliases.iter().map(|alias| alias.name.clone()))
                                .chain(std::iter::once(command.usage.clone()))
                                .collect::<Vec<_>>()
                                .join(" ");
                            SearchEntry {
                                title: command.path.clone(),
                                summary: compact_help_summary(&command.about),
                                href: cli_command_route(&command.path),
                                kind: "command".to_owned(),
                                text: format!("{} {text}", command.about),
                            }
                        }),
                );
                entries.extend(reference.settings.into_iter().map(|setting| SearchEntry {
                    title: setting.path.clone(),
                    summary: compact_help_summary(&setting.description),
                    href: format!(
                        "{}#{}",
                        cli_setting_namespace_route(cli_setting_namespace(&setting.path)),
                        generated_anchor("setting", &setting.path)
                    ),
                    kind: "setting".to_owned(),
                    text: format!(
                        "{} {} {} {}",
                        setting.description,
                        setting.value_type,
                        setting.possible_values.join(" "),
                        setting
                            .default
                            .as_ref()
                            .map(Value::to_string)
                            .unwrap_or_default(),
                    ),
                }));
                entries.extend(namespaces.into_iter().map(|namespace| SearchEntry {
                    title: format!("{namespace} settings"),
                    summary: "GammaLoop settings namespace".to_owned(),
                    href: cli_setting_namespace_route(&namespace),
                    kind: "settings namespace".to_owned(),
                    text: namespace,
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
                            text: String::new(),
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
                            text: String::new(),
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

    fn write_federated_search_index(
        &self,
        output: &Path,
        channel: BuildChannel,
        tag: Option<&str>,
    ) -> Result<()> {
        let channel_route = match channel {
            BuildChannel::Latest => "latest".to_owned(),
            BuildChannel::Snapshot => format!(
                "snapshots/{}",
                tag.expect("snapshot tag was validated before federating search")
            ),
        };
        let mut entries = Vec::new();
        for product in &self.registry.product {
            let product_root = output
                .join("products")
                .join(&product.id)
                .join(&channel_route);
            let path = product_root.join("search-index.json");
            let product_entries = serde_json::from_slice::<Vec<SearchEntry>>(&fs::read(&path)?)
                .wrap_err_with(|| format!("failed to parse {}", path.display()))?;
            entries.extend(product_entries.into_iter().map(|mut entry| {
                entry.href = format!("products/{}/{channel_route}/{}", product.id, entry.href);
                entry.kind = format!("{} · {}", product.title, entry.kind);
                entry
            }));
        }

        let developer_index = output.join("developers/search-index.json");
        let developer_entries =
            serde_json::from_slice::<Vec<SearchEntry>>(&fs::read(&developer_index)?)
                .wrap_err_with(|| format!("failed to parse {}", developer_index.display()))?;
        entries.extend(developer_entries.into_iter().map(|mut entry| {
            entry.href = format!("developers/{}", entry.href);
            entry.kind = format!("Developers · {}", entry.kind);
            entry
        }));
        entries.extend(self.registry.product.iter().map(|product| {
            let persistent = product.citation.doi.as_ref().map_or_else(
                || "the versioned source repository".to_owned(),
                |doi| format!("DOI {doi}"),
            );
            SearchEntry {
                title: format!("Cite {}", product.title),
                summary: format!(
                    "Version-specific software citation and BibTeX for {} using {persistent}.",
                    product.title
                ),
                href: format!("citations/#{}", product.id),
                kind: "citation".to_owned(),
                text: format!(
                    "{} {} software citation cite BibTeX",
                    product.citation.title,
                    product
                        .citation
                        .doi
                        .as_deref()
                        .unwrap_or(&product.citation.repository),
                ),
            }
        }));
        entries.push(SearchEntry {
            title: "About the αLoop collaboration".to_owned(),
            summary: self.portal.summary.clone(),
            href: "about/".to_owned(),
            kind: "collaboration".to_owned(),
            text: self
                .portal
                .pillar
                .iter()
                .map(|pillar| format!("{} {}", pillar.title, pillar.summary))
                .collect::<Vec<_>>()
                .join(" "),
        });
        entries.push(SearchEntry {
            title: "Talks".to_owned(),
            summary: "Seminars and conference presentations by αLoop collaborators.".to_owned(),
            href: "talks/".to_owned(),
            kind: "collaboration".to_owned(),
            text: "Local Unitarity GammaLoop numerical perturbation theory".to_owned(),
        });
        entries.extend(self.talks.talk.iter().map(|talk| {
            let speaker = self
                .portal
                .people
                .iter()
                .find(|person| person.id == talk.speaker)
                .expect("talk speaker was validated");
            SearchEntry {
                title: talk.title.clone(),
                summary: format!("{} · {} · {}", speaker.name, talk.event, talk.date),
                href: format!("talks/#{}", talk.id),
                kind: "talk".to_owned(),
                text: format!("{} {} {}", talk.location, speaker.name, talk.event),
            }
        }));

        let mut seen = BTreeSet::new();
        entries.retain(|entry| seen.insert((entry.href.clone(), entry.title.clone())));
        fs::write(
            output.join("search-index.json"),
            serde_json::to_vec_pretty(&entries)?,
        )?;
        Ok(())
    }

    fn render_developer_typst(
        &self,
        notes: &[&DeveloperNote],
        dependency_output: Option<&Path>,
    ) -> Result<Option<TempDir>> {
        if notes.is_empty() {
            return Ok(None);
        }

        let target = self.root.join("target");
        fs::create_dir_all(&target)?;
        let work = TempDirBuilder::new()
            .prefix("alphal00p-developer-typst-")
            .tempdir_in(target)?;
        let wrapper = work.path().join("developers.typ");
        let mut documents = String::new();
        for (index, note) in notes.iter().enumerate() {
            documents.push_str(&format!(
                "#document(\"note-{index}.html\", title: [#raw(\"{}\")])[#include \"/{}\"]\n",
                typst_string(&note.title),
                typst_string(&note.source.to_string_lossy().replace('\\', "/")),
            ));
        }
        fs::write(&wrapper, documents)?;

        let bundle = work.path().join("bundle");
        let timestamp = self.git_timestamp();
        ensure!(
            timestamp > 0,
            "cannot determine the documented commit timestamp; set SOURCE_DATE_EPOCH"
        );
        let commit = self.git_commit();
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
            &timestamp.to_string(),
            "--input",
            &format!("git-commit={commit}"),
        ]);
        if let Some(directory) = dependency_output {
            fs::create_dir_all(directory)?;
            command
                .args(["--deps-format", "zero", "--deps"])
                .arg(directory.join("developers.deps"));
        }
        let status = command
            .arg(&wrapper)
            .arg(&bundle)
            .status()
            .wrap_err("failed to launch Typst 0.15; enter the Nix development shell")?;
        ensure!(status.success(), "Typst failed for developer notes");
        for (index, _) in notes.iter().enumerate() {
            ensure!(
                bundle.join(format!("note-{index}.html")).is_file(),
                "Typst emitted no developer note {index}"
            );
        }
        Ok(Some(work))
    }

    fn write_developer_docs(
        &self,
        output: &Path,
        include_typst: bool,
        dependency_output: Option<&Path>,
    ) -> Result<()> {
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
        let typst_work = if include_typst {
            self.render_developer_typst(&notes, dependency_output)?
        } else {
            None
        };
        let typst_bundle = typst_work.as_ref().map(|work| work.path().join("bundle"));

        for (index, note) in notes.iter().enumerate() {
            let page = &pages[index];
            let docs_root = page_root_prefix(&page.route);
            let destination = developer_root.join(&page.route);
            fs::create_dir_all(&destination)?;
            let source_url = format!(
                "https://github.com/alphal00p/gammaloop/blob/{}/{}",
                commit,
                note.source.to_string_lossy()
            );
            let (body, extra_head) = if let Some(bundle) = &typst_bundle {
                let rendered = fs::read_to_string(bundle.join(format!("note-{index}.html")))?;
                let mut body = extract_html_body(&rendered)?.to_owned();
                validate_developer_typst_body(&body)?;
                let source_title = leading_developer_title(&body)?;
                ensure!(
                    source_title.title == note.title,
                    "developer Typst title {:?} does not match registry title {:?}",
                    source_title.title,
                    note.title
                );
                body.replace_range(source_title.range, "");
                promote_heading_hierarchy(&mut body);
                let body = rewrite_developer_source_links(
                    &body,
                    &note.source,
                    &commit,
                    &note_routes,
                    &self.root,
                )?;
                (body, extract_typst_head_styles(&rendered)?)
            } else {
                (
                    "<p>Typst rendering was disabled for this metadata-only build.</p>".to_owned(),
                    String::new(),
                )
            };
            let status_class = slug(&note.status);
            let owner = self
                .developers
                .owner
                .iter()
                .find(|owner| owner.id == note.owner)
                .expect("developer owners were validated before rendering");
            let review_age = note.reviewed_at.as_deref().and_then(|date| {
                NaiveDate::parse_from_str(date, "%Y-%m-%d")
                    .ok()
                    .map(|date| {
                        Utc::now()
                            .date_naive()
                            .signed_duration_since(date)
                            .num_days()
                    })
            });
            let freshness = match (note.lifecycle.as_str(), review_age) {
                ("current", Some(age)) if age > 90 => "stale",
                ("current", Some(age)) if age > 60 => "review-due",
                ("proposal", Some(age)) if age > 180 => "disposition-due",
                ("current" | "proposal", Some(_)) => "reviewed",
                ("current" | "proposal", None) => "unreviewed",
                _ => "frozen-evidence",
            };
            let review = note.reviewed_at.as_deref().map_or_else(
                || "No review record".to_owned(),
                |date| {
                    format!(
                        "{} · <code>{}</code>",
                        escape_html(date),
                        escape_html(note.review_ref.as_deref().unwrap_or("missing review ref"))
                    )
                },
            );
            let evidence = if matches!(note.lifecycle.as_str(), "investigation" | "archived") {
                format!(
                    "{}{}",
                    note.captured_at
                        .as_deref()
                        .map(|date| format!("captured {}", escape_html(date)))
                        .unwrap_or_else(|| "capture date unavailable".to_owned()),
                    note.evidence_revision
                        .as_deref()
                        .map_or_else(String::new, |revision| {
                            format!(" · <code>{}</code>", escape_html(revision))
                        })
                )
            } else {
                "Not an evidence record".to_owned()
            };
            let scopes = if note.scope.is_empty() {
                "<span class=\"developer-missing\">No verified code scopes recorded</span>"
                    .to_owned()
            } else {
                let items = note
                    .scope
                    .iter()
                    .map(|scope| {
                        format!(
                            "<li><code>{}</code><br><span>{}</span> · digest <code>{}</code></li>",
                            escape_html(&scope.symbol),
                            escape_html(&scope.path.to_string_lossy()),
                            escape_html(&scope.digest.chars().take(12).collect::<String>()),
                        )
                    })
                    .collect::<String>();
                format!("<ul>{items}</ul>")
            };
            let record = format!(
                "<dl class=\"developer-record-meta\" data-freshness=\"{}\"><div><dt>Lifecycle</dt><dd>{}</dd></div><div><dt>Owner</dt><dd><a href=\"{}\">{}</a></dd></div><div><dt>Review</dt><dd>{review}</dd></div><div><dt>Freshness</dt><dd>{}</dd></div><div><dt>Evidence</dt><dd>{evidence}</dd></div><div class=\"developer-record-scopes\"><dt>Verified scopes</dt><dd>{scopes}</dd></div></dl>",
                escape_html(freshness),
                escape_html(&note.lifecycle),
                escape_html(&owner.contact),
                escape_html(&owner.name),
                escape_html(freshness),
            );
            let article = format!(
                "<header class=\"developer-note-hero\"><p class=\"product-eyebrow\">For developers</p><div class=\"developer-note-title\"><h1>{}</h1><div class=\"developer-note-badges\"><span class=\"developer-status developer-status-{}\">{}</span><span class=\"developer-kind\">{}</span><span class=\"developer-lifecycle\">{}</span></div></div><p>{}</p><p class=\"developer-source\"><a href=\"{}\">View the source note <span aria-hidden=\"true\">↗</span></a></p></header>{record}{body}",
                escape_html(&note.title),
                escape_html(&status_class),
                escape_html(&note.status),
                escape_html(&note.kind),
                escape_html(&note.lifecycle),
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
            .replace("Documentation pagination", "Developer note pagination");
            let portal_root = format!("{docs_root}../");
            let issue_url = documentation_issue_url(
                &format!("Developer note: {}", note.title),
                &format!("developers/{}", page.route),
                &commit,
            );
            let feedback = format!(
                "<nav class=\"page-feedback\" aria-label=\"Documentation feedback\"><a href=\"{}\">View page source at this revision <span aria-hidden=\"true\">↗</span></a><a href=\"{}\">Report a documentation issue <span aria-hidden=\"true\">↗</span></a></nav>",
                escape_html(&source_url),
                escape_html(&issue_url),
            );
            let main = format!(
                "<nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"{}\">αLoop</a> / <a href=\"{}\">Developers</a> / {}</nav>{inline_toc}<article class=\"docs-article developer-article\">{article}</article>{page_navigation}{feedback}<footer class=\"page-footer\">Developer architecture · documented revision <code>{}</code></footer>",
                escape_html(&portal_root),
                escape_html(&docs_root),
                escape_html(&note.title),
                escape_html(&commit.chars().take(12).collect::<String>()),
            );
            fs::write(
                destination.join("index.html"),
                developer_document(
                    (
                        &format!("{} · Developer architecture", note.title),
                        &note.summary,
                    ),
                    &docs_root,
                    &portal_root,
                    &extra_head,
                    &sidebar,
                    &main,
                    &toc,
                ),
            )?;

            if note.lifecycle != "superseded" {
                let classification = developer_search_classification(note);
                search.push(SearchEntry {
                    title: note.title.clone(),
                    summary: note.summary.clone(),
                    href: page.route.clone(),
                    kind: format!("developer {classification}"),
                    text: String::new(),
                });
                search.extend(
                    headings
                        .into_iter()
                        .filter(|heading| matches!(heading.level, 2 | 3))
                        .map(|heading| developer_heading_search_entry(note, &page.route, heading)),
                );
            }
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
                let render_card = |note: &DeveloperNote| {
                    format!(
                        "<article class=\"developer-card\"><div class=\"developer-card-meta\"><span class=\"developer-status developer-status-{}\">{}</span><span class=\"developer-lifecycle\">{}</span><span>{}</span></div><h3><a href=\"architecture/{}/\">{}</a></h3><p>{}</p><a class=\"portal-text-link\" href=\"architecture/{}/\">Read note <span aria-hidden=\"true\">→</span></a></article>",
                        escape_html(&slug(&note.status)),
                        escape_html(&note.status),
                        escape_html(&note.lifecycle),
                        escape_html(&note.kind),
                        escape_html(&note.id),
                        escape_html(&note.title),
                        escape_html(&note.summary),
                        escape_html(&note.id),
                    )
                };
                let notes = section
                    .note
                    .iter()
                    .filter(|note| note.lifecycle != "superseded")
                    .map(&render_card)
                    .collect::<String>();
                let superseded = section
                    .note
                    .iter()
                    .filter(|note| note.lifecycle == "superseded")
                    .map(render_card)
                    .collect::<String>();
                let superseded = if superseded.is_empty() {
                    String::new()
                } else {
                    format!(
                        "<details class=\"developer-archive\"><summary>Superseded records</summary><div class=\"developer-card-grid\">{superseded}</div></details>"
                    )
                };
                format!(
                    "<section class=\"developer-section\" id=\"{}\"><div class=\"developer-section-heading\"><p class=\"portal-kicker\">Developer notes</p><h2>{}</h2><p>{}</p></div><div class=\"developer-card-grid\">{notes}{superseded}</div></section>",
                    escape_html(&section.id),
                    escape_html(&section.title),
                    escape_html(&section.summary),
                )
            })
            .collect::<String>();
        let sidebar = self.developer_sidebar(None, "");
        let source_url =
            format!("https://github.com/alphal00p/gammaloop/blob/{commit}/docs/developers.toml");
        let issue_url = documentation_issue_url("Developer architecture", "developers/", commit);
        let feedback = format!(
            "<nav class=\"page-feedback\" aria-label=\"Documentation feedback\"><a href=\"{}\">View page source at this revision <span aria-hidden=\"true\">↗</span></a><a href=\"{}\">Report a documentation issue <span aria-hidden=\"true\">↗</span></a></nav>",
            escape_html(&source_url),
            escape_html(&issue_url),
        );
        let main = format!(
            "<nav class=\"breadcrumbs\" aria-label=\"Breadcrumb\"><a href=\"../\">αLoop</a> / Developers</nav><article class=\"docs-article developer-article\"><header class=\"developer-hero\"><p class=\"product-eyebrow\">For developers</p><h1>{}</h1><p>{}</p><aside class=\"developer-audience\" aria-label=\"Documentation audience\"><strong>Looking for usage documentation?</strong><span>Product guides and reference live with each research project. This area documents implementation details, design work, and engineering investigations for contributors.</span><a href=\"../#projects\">Browse project documentation <span aria-hidden=\"true\">→</span></a></aside></header>{sections}</article>{feedback}<footer class=\"page-footer\">{} classified notes · documented revision <code>{}</code></footer>",
            escape_html(&self.developers.title),
            escape_html(&self.developers.summary),
            self.developers
                .section
                .iter()
                .map(|section| section.note.len())
                .sum::<usize>(),
            escape_html(&commit.chars().take(12).collect::<String>()),
        );
        fs::write(
            developer_root.join("index.html"),
            developer_document(
                (&self.developers.title, &self.developers.summary),
                "",
                "../",
                "",
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
                .filter(|note| note.lifecycle != "superseded")
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
            "<aside class=\"docs-sidebar developer-sidebar\" id=\"docs-sidebar\" aria-label=\"Developer architecture\">{navigation}<p class=\"sidebar-meta\"><strong>For developers</strong><br>Implementation and engineering notes.<br><a href=\"{}../#projects\">Product guides and reference</a></p></aside>",
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
        let tasks = PORTAL_TASKS
            .iter()
            .enumerate()
            .map(|(index, (task, product_id, role))| -> Result<String> {
                let product = self
                    .registry
                    .product
                    .iter()
                    .find(|product| product.id == *product_id)
                    .wrap_err_with(|| format!("portal task {task} references unknown product {product_id}"))?;
                let quickstart = product
                    .pages
                    .iter()
                    .find(|page| page.id == "quickstart")
                    .wrap_err_with(|| format!("{} has no quickstart route", product.id))?;
                Ok(format!(
                    r#"<a class="portal-task-link" data-product="{}" href="products/{}/{}/{}"><span class="portal-task-order">{:02}</span><span class="portal-task-copy"><strong>{}</strong><span>{} · {}</span><small>{}</small></span><span class="portal-task-arrow" aria-hidden="true">→</span></a>"#,
                    escape_html(&product.id),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&quickstart.route),
                    index + 1,
                    escape_html(task),
                    escape_html(&product.title),
                    escape_html(role),
                    escape_html(&product.tagline),
                ))
            })
            .collect::<Result<Vec<_>>>()?
            .join("");
        let projects = self
            .registry
            .product
            .iter()
            .enumerate()
            .map(|(index, product)| -> Result<String> {
                let role = PORTAL_TASKS
                    .iter()
                    .find(|(_, product_id, _)| *product_id == product.id)
                    .map(|(_, _, role)| *role)
                    .wrap_err_with(|| format!("{} has no portal ecosystem role", product.id))?;
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
                let guide_route = product
                    .pages
                    .iter()
                    .find(|page| page.group == "Guides")
                    .map(|page| page.route.as_str())
                    .wrap_err_with(|| format!("{} has no portal guide", product.id))?;
                let quickstart_route = product
                    .pages
                    .iter()
                    .find(|page| page.id == "quickstart")
                    .map(|page| page.route.as_str())
                    .wrap_err_with(|| format!("{} has no portal quickstart", product.id))?;
                Ok(format!(
                    r#"<article class="portal-project-card" data-product="{}"><div class="portal-project-meta"><span>{:02}</span><span>{}</span></div><h3><a href="products/{}/{}/">{}</a></h3><p class="portal-project-summary">{}</p><div class="portal-packages" aria-label="{} crates and modules"><span class="portal-packages-label">Crates &amp; modules</span>{}</div><nav class="portal-card-links" aria-label="{} documentation"><a class="portal-card-primary" href="products/{}/{}/">Overview <span aria-hidden="true">↗</span></a><a href="products/{}/{}/{}">Get started</a><a href="products/{}/{}/{}">Guides</a><a href="products/{}/{}/reference/">Reference</a><a class="portal-card-cite" href="citations/#{}">Cite</a></nav></article>"#,
                    escape_html(&product.id),
                    index + 1,
                    escape_html(role),
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
                    escape_html(quickstart_route),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(guide_route),
                    escape_html(&product.id),
                    escape_html(&channel_route),
                    escape_html(&product.id),
                ))
            })
            .collect::<Result<Vec<_>>>()?
            .join("");
        let process_graphs = PORTAL_GRAPH_IDS
            .iter()
            .map(|graph| {
                format!(
                    r#"<span class="portal-process-graph"><img class="portal-graph-theme portal-graph-theme-light" src="assets/graphs/portal-graph-{graph}-light.svg" alt=""><img class="portal-graph-theme portal-graph-theme-dark" src="assets/graphs/portal-graph-{graph}-dark.svg" alt=""></span>"#,
                )
            })
            .collect::<String>();
        let funding = format!(
            r#"<aside class="portal-funding" aria-labelledby="funding-title"><span class="portal-funding-mark" aria-hidden="true">α</span><div class="portal-funding-copy"><p class="portal-kicker">Funding</p><h2 id="funding-title">Publicly funded research</h2><p>{}</p></div><a class="portal-text-link" href="{}">Funding record <span aria-hidden="true">↗</span></a></aside>"#,
            escape_html(&self.portal.funding),
            escape_html(&self.portal.funding_url),
        );
        let favicon = favicon_links("assets/");
        let html = format!(
            r##"<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="{}"><meta name="theme-color" content="#f9f6f0">{favicon}<title>αLoop · Research software for collider physics</title><link rel="stylesheet" href="assets/site.css"><script defer src="assets/site.js"></script></head><body class="portal-body" data-search-index="search-index.json" data-search-root=""><a class="skip-link" href="#main-content">Skip to content</a><header class="portal-header"><a class="portal-brand" href="#overview" aria-label="αLoop home"><span class="portal-brand-logo" aria-hidden="true"></span><span class="portal-brand-copy"><strong>αLoop</strong><small>Local Unitarity research</small></span></a><nav class="portal-nav" aria-label="Primary"><a href="#tasks">Tasks</a><a href="people/">People</a><a href="publications/">Publications</a><a href="developers/">Developers</a></nav><div class="portal-header-actions"><button class="portal-search-button" type="button" data-search-open>Search <span class="header-button-label">⌘K</span></button><a class="portal-source-link" href="https://github.com/alphal00p/gammaloop">GitHub <span aria-hidden="true">↗</span></a><button class="portal-theme-button" type="button" data-theme-toggle aria-label="Toggle color theme"><span aria-hidden="true">◐</span></button></div></header><main class="portal-main" id="main-content"><section class="portal-hero portal-section" id="overview"><div class="portal-hero-copy"><p class="portal-kicker">{}</p><h1>{}</h1><p class="portal-lede">{}</p><div class="portal-hero-actions"><a class="portal-button portal-button-primary" href="#tasks">Choose a task <span aria-hidden="true">↓</span></a><a class="portal-button" href="products/gammaloop/{}/quickstart/">Start with GammaLoop <span aria-hidden="true">↗</span></a><a class="portal-button" href="citations/">Cite the software <span aria-hidden="true">↗</span></a></div></div><div class="portal-hero-art"><div class="portal-wordmark" aria-label="αLoop collaboration mark" role="img"></div><div class="portal-graph-field" role="img" aria-label="A jumble of Feynman graphs rendered by GammaLoop from real process and test data">{process_graphs}</div><p>Local cancellation.<br>Global precision.</p></div></section><section class="portal-section portal-task-chooser" id="tasks" aria-labelledby="tasks-title"><header class="portal-task-heading"><div><p class="portal-kicker">Choose by task · 01—05</p><h2 id="tasks-title">What do you want to work on?</h2></div><p>Start with the scientific object or operation you have in hand. Each route opens the maintained first workflow for the component that owns it.</p></header><nav class="portal-task-grid" aria-label="Documentation by task">{tasks}</nav></section><section class="portal-section portal-projects" id="projects" aria-labelledby="projects-title"><div class="portal-section-heading"><div><p class="portal-kicker">Research software · 01—05</p><h2 id="projects-title">Projects &amp; crates</h2></div><p>Five connected codebases spanning numerical cross-sections, graph algorithms, tensor networks, symbolic identities, and integral evaluation.</p></div><div class="portal-project-grid">{projects}</div></section>{funding}</main><footer class="portal-footer"><div><span class="portal-footer-mark" aria-hidden="true"></span><p><strong>αLoop</strong><br>Local Unitarity research software</p></div><nav aria-label="Footer"><a href="#tasks">Tasks</a><a href="#projects">Projects</a><a href="people/">People</a><a href="publications/">Publications</a><a href="citations/">Cite</a><a href="developers/">Developers</a><a href="https://github.com/alphal00p/gammaloop">Source</a></nav><p>Physics, algorithms, and software<br>developed in the open.</p></footer><dialog class="search-dialog" data-search-dialog><form class="search-form" method="dialog"><input class="search-input" type="search" data-search-input placeholder="Search all projects and developer notes" aria-label="Search all documentation"><button class="header-button" value="close">Close</button></form><ul class="search-results" data-search-results aria-live="polite"></ul></dialog></body></html>"##,
            escape_html(&self.portal.summary),
            escape_html(&self.portal.eyebrow),
            escape_html(&self.portal.title),
            escape_html(&self.portal.summary),
            escape_html(&channel_route),
        );
        let html = html
            .replace(
                "<nav class=\"portal-nav\" aria-label=\"Primary\"><a href=\"#tasks\">Tasks</a><a href=\"people/\">People</a><a href=\"publications/\">Publications</a><a href=\"developers/\">Developers</a></nav>",
                "<nav class=\"portal-nav\" aria-label=\"Primary\"><a href=\"about/\">About</a><a href=\"#projects\">Projects</a><a href=\"people/\">People</a><a href=\"talks/\">Talks</a><a href=\"publications/\">Publications</a><a href=\"developers/\">Developers</a></nav>",
            )
            .replace(
                "<nav aria-label=\"Footer\"><a href=\"#tasks\">Tasks</a><a href=\"#projects\">Projects</a><a href=\"people/\">People</a><a href=\"publications/\">Publications</a><a href=\"citations/\">Cite</a><a href=\"developers/\">Developers</a><a href=\"https://github.com/alphal00p/gammaloop\">Source</a></nav>",
                "<nav aria-label=\"Footer\"><a href=\"about/\">About</a><a href=\"#projects\">Projects</a><a href=\"people/\">People</a><a href=\"talks/\">Talks</a><a href=\"publications/\">Publications</a><a href=\"citations/\">Cite</a><a href=\"developers/\">Developers</a><a href=\"https://github.com/alphal00p/gammaloop\">Source</a></nav>",
            );
        fs::write(output.join("index.html"), html)?;
        self.write_people_page(output)?;
        self.write_about_page(output)?;
        self.write_talks_page(output)?;
        self.write_publications_page(output)?;
        self.write_citations_page(output)?;
        fs::write(output.join(".nojekyll"), b"")?;
        Ok(())
    }

    fn write_about_page(&self, output: &Path) -> Result<()> {
        let pillars = self
            .portal
            .pillar
            .iter()
            .map(|pillar| {
                format!(
                    "<article class=\"about-pillar\"><p class=\"portal-kicker\">{}</p><h2>{}</h2><p>{}</p></article>",
                    escape_html(&pillar.label),
                    escape_html(&pillar.title),
                    escape_html(&pillar.summary),
                )
            })
            .collect::<String>();
        let affiliations = self
            .portal
            .affiliation
            .iter()
            .map(|affiliation| {
                format!(
                    "<a class=\"about-affiliation\" href=\"{}\"><span>{}</span><strong>{}</strong><small>{}</small><p>{}</p><b aria-hidden=\"true\">↗</b></a>",
                    escape_html(&affiliation.url),
                    escape_html(&affiliation.location),
                    escape_html(&affiliation.name),
                    escape_html(&affiliation.location),
                    escape_html(&affiliation.summary),
                )
            })
            .collect::<String>();
        let body = format!(
            r#"<header class="portal-page-hero about-page-hero"><p class="portal-kicker">About the collaboration</p><h1>Precision through local cancellation.</h1><p>{}</p></header><section class="about-origin"><div class="about-origin-copy"><p class="portal-kicker">Why αLoop</p><h2>Precision is another path to discovery.</h2><p>The lack of obvious sign of new physics phenomenon in collider experiments is an opportunity to take a step back and reflect on the amazing theory we have discovered so far: the Standard Model. In particular, we must now strive to make ever more precise predictions so as to hunt for indirect evidence of new physics in small departure from expectations.</p><p>For this reason, our collaboration is dedicated to theoretical and algorithmic research for the automated computation of cross-sections in Quantum Field Theories at arbitrary perturbative orders. In particular, we develop a new theoretical framework called <em>Local Unitarity</em> (LU) which approaches this problem from an unorthodox way, particularly suited to numerical computations.</p><nav aria-label="Learn about Local Unitarity"><a class="portal-button portal-button-primary" href="https://arxiv.org/abs/2110.15662">Read the introduction <span aria-hidden="true">↗</span></a><a class="portal-button" href="../publications/">Explore publications <span aria-hidden="true">→</span></a></nav></div><aside class="about-equation" aria-label="Schematic Local Unitarity cross-section"><div class="about-equation-illustration"><img src="../assets/about-double-triangle-light.svg" alt="" class="about-equation-graph portal-graph-theme-light"><img src="../assets/about-double-triangle-dark.svg" alt="" class="about-equation-graph portal-graph-theme-dark"></div><div class="about-equation-formula" role="img" aria-label="The differential cross section is a sum over graphs of loop-momentum integrals and a sum over cuts of the Local Unitarity integrand, constrained by the observable."><img src="../assets/about-local-unitarity-equation-light.svg" alt="" class="portal-graph-theme-light"><img src="../assets/about-local-unitarity-equation-dark.svg" alt="" class="portal-graph-theme-dark"></div><small>Real and virtual contributions share one numerical representation.</small></aside></section><section class="about-pillars" aria-labelledby="about-pillars-title"><header><p class="portal-kicker">From method to software</p><h2 id="about-pillars-title">One research programme, connected structures.</h2></header><div>{pillars}</div></section><section class="about-affiliations" aria-labelledby="about-affiliations-title"><header><p class="portal-kicker">Affiliations</p><h2 id="about-affiliations-title">Research across institutions.</h2><p>αLoop connects collider-physics research, mathematical structures, and open scientific-software development.</p></header><div>{affiliations}</div></section><aside class="portal-funding about-funding" aria-labelledby="about-funding-title"><span class="portal-funding-mark" aria-hidden="true">α</span><div class="portal-funding-copy"><p class="portal-kicker">Funding</p><h2 id="about-funding-title">Publicly funded research</h2><p>{}</p></div><a class="portal-text-link" href="{}">Funding record <span aria-hidden="true">↗</span></a></aside><section class="about-next"><p class="portal-kicker">The collaboration</p><h2>Meet the people doing the work.</h2><nav aria-label="Explore the collaboration"><a class="portal-button portal-button-primary" href="../people/">People <span aria-hidden="true">→</span></a><a class="portal-button" href="../talks/">Talks <span aria-hidden="true">→</span></a></nav></section>"#,
            escape_html(&self.portal.summary),
            escape_html(&self.portal.funding),
            escape_html(&self.portal.funding_url),
        );
        let directory = output.join("about");
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            portal_subpage_document(
                "About",
                "The αLoop collaboration develops Local Unitarity methods and open research software for precision collider physics.",
                "about",
                &body,
            ),
        )?;
        Ok(())
    }

    fn write_talks_page(&self, output: &Path) -> Result<()> {
        let people = self
            .portal
            .people
            .iter()
            .map(|person| (person.id.as_str(), person))
            .collect::<BTreeMap<_, _>>();
        let mut talks = self.talks.talk.iter().collect::<Vec<_>>();
        talks.sort_by(|left, right| {
            right
                .date
                .cmp(&left.date)
                .then(left.title.cmp(&right.title))
        });
        let mut years = BTreeMap::<i32, Vec<&TalkConfig>>::new();
        for talk in talks {
            let date = NaiveDate::parse_from_str(&talk.date, "%Y-%m-%d")?;
            years
                .entry(date.format("%Y").to_string().parse()?)
                .or_default()
                .push(talk);
        }
        let timeline = years
            .into_iter()
            .rev()
            .map(|(year, talks)| {
                let cards = talks
                    .into_iter()
                    .map(|talk| {
                        let person = people
                            .get(talk.speaker.as_str())
                            .expect("talk speaker was validated");
                        let date = NaiveDate::parse_from_str(&talk.date, "%Y-%m-%d")
                            .expect("talk date was validated");
                        let mut links = format!(
                            "<a href=\"{}\">Event record <span aria-hidden=\"true\">↗</span></a>",
                            escape_html(&talk.event_url),
                        );
                        if let Some(url) = &talk.slides_url {
                            links.push_str(&format!(
                                "<a href=\"{}\">Slides <span aria-hidden=\"true\">↗</span></a>",
                                escape_html(url),
                            ));
                        }
                        if let Some(url) = &talk.recording_url {
                            links.push_str(&format!(
                                "<a href=\"{}\">Recording <span aria-hidden=\"true\">↗</span></a>",
                                escape_html(url),
                            ));
                        }
                        format!(
                            "<article class=\"talk-card\" id=\"{}\"><div class=\"talk-card-date\"><time datetime=\"{}\">{}</time><span>{}</span></div><div class=\"talk-card-copy\"><p class=\"portal-kicker\"><a href=\"../people/#{}\">{}</a></p><h3>{}</h3><p><strong>{}</strong><br>{}</p><nav aria-label=\"Resources for {}\">{links}</nav></div></article>",
                            escape_html(&talk.id),
                            escape_html(&talk.date),
                            date.format("%d %b"),
                            year,
                            escape_html(&person.id),
                            escape_html(&person.name),
                            escape_html(&talk.title),
                            escape_html(&talk.event),
                            escape_html(&talk.location),
                            escape_html(&talk.title),
                        )
                    })
                    .collect::<String>();
                format!(
                    "<section class=\"talk-year\" aria-labelledby=\"talk-year-{year}\"><h2 id=\"talk-year-{year}\">{year}</h2><div>{cards}</div></section>"
                )
            })
            .collect::<String>();
        let body = format!(
            "<header class=\"portal-page-hero talks-page-hero\"><p class=\"portal-kicker\">Seminars &amp; conferences</p><h1>Talks</h1><p>Selected presentations on Local Unitarity, numerical perturbation theory, GammaLoop, and the scientific software surrounding them.</p><p class=\"talks-provenance\">{} talks · linked to public conference records, slides, and recordings where available.</p></header><div class=\"talk-timeline\">{timeline}</div>",
            self.talks.talk.len(),
        );
        let directory = output.join("talks");
        fs::create_dir_all(&directory)?;
        fs::write(
            directory.join("index.html"),
            portal_subpage_document(
                "Talks",
                "Talks by αLoop collaborators on Local Unitarity, numerical methods, GammaLoop, and related research software.",
                "talks",
                &body,
            ),
        )?;
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

    fn write_product_preview(
        &self,
        output: &Path,
        product: &ProductConfig,
        channel: BuildChannel,
        tag: Option<&str>,
    ) -> Result<()> {
        let route = match channel {
            BuildChannel::Latest => format!("products/{}/latest/", product.id),
            BuildChannel::Snapshot => format!(
                "products/{}/snapshots/{}/",
                product.id,
                tag.expect("snapshot tag was validated before rendering the preview")
            ),
        };
        ensure!(
            output.join(&route).join("index.html").is_file(),
            "product preview has no entry page at {route}"
        );
        fs::write(
            output.join("index.html"),
            generated_reference_redirect(&product.title, &route, "product preview"),
        )?;
        fs::write(
            output.join(".note"),
            format!(
                "alphal00p product preview\nproduct={}\nchannel={}\nroute={route}\n",
                product.id,
                match channel {
                    BuildChannel::Latest => "latest",
                    BuildChannel::Snapshot => "snapshot",
                },
            ),
        )?;
        fs::write(output.join(".nojekyll"), b"")?;
        Ok(())
    }

    fn validate_generated_links(
        &self,
        output: &Path,
        include_rustdoc: bool,
        snapshot_tag: Option<&str>,
    ) -> Result<()> {
        let patterns = LinkValidationPatterns::new()?;
        let documented_revision = self.git_commit();
        let roots = [output.to_path_buf()];
        let mut failures = vec![];
        let mut linked_pages = HashMap::new();
        let mut link_rewrites = LinkRewriteIndex::new();
        let mut local_paths = LocalPathIndex::new();
        for root in roots {
            for entry in WalkDir::new(&root) {
                let entry = entry?;
                let relative = entry.path().strip_prefix(output)?;
                if entry.file_type().is_file()
                    && relative.starts_with("products")
                    && entry.path().extension().and_then(|value| value.to_str()) == Some("typ")
                {
                    let parts = relative.components().collect::<Vec<_>>();
                    ensure!(
                        parts.windows(2).any(|parts| {
                            parts[0].as_os_str() == "assets" && parts[1].as_os_str() == "typst"
                        }),
                        "authored Typst source was published instead of rendered: {}",
                        relative.display()
                    );
                }
                local_paths.insert(entry.path().to_path_buf(), Some(entry.file_type()));
                let rustdoc = entry
                    .path()
                    .components()
                    .collect::<Vec<_>>()
                    .windows(2)
                    .any(|parts| {
                        parts[0].as_os_str() == "reference" && parts[1].as_os_str() == "rust"
                    });
                if !entry.file_type().is_file()
                    || entry.path().extension().and_then(|value| value.to_str()) != Some("html")
                    || (!include_rustdoc && rustdoc)
                {
                    continue;
                }
                let html = fs::read_to_string(entry.path())?;
                let rendered_html = patterns.rendered_html(&html);
                if !linked_pages.contains_key(entry.path()) {
                    linked_pages.insert(
                        entry.path().to_path_buf(),
                        patterns.page_metadata(&rendered_html),
                    );
                }
                for element in patterns.tag.find_iter(&rendered_html) {
                    for capture in patterns.href.captures_iter(element.as_str()) {
                        let target_capture = capture
                            .get(1)
                            .or_else(|| capture.get(2))
                            .or_else(|| capture.get(3))
                            .expect("one URL capture exists");
                        let target = target_capture.as_str();
                        let target_path = target.split(['?', '#']).next().unwrap_or_default();
                        if !include_rustdoc && target_path.contains("reference/rust/") {
                            continue;
                        }
                        if rustdoc
                            && target_path.ends_with(".js")
                            && target_path
                                .split('/')
                                .any(|segment| segment == "trait.impl")
                        {
                            // Rustdoc emits optional async implementor sidecars even when
                            // --no-deps leaves that sidecar absent from this crate-only tree.
                            continue;
                        }
                        let attribute_capture = capture.get(0).expect("URL attribute capture");
                        let offset = element.start();
                        self.validate_generated_target(
                            output,
                            LinkReference {
                                source: entry.path(),
                                href: target,
                                attribute: Some(LinkAttribute {
                                    attribute_range: offset + attribute_capture.start()
                                        ..offset + attribute_capture.end(),
                                    target_range: offset + target_capture.start()
                                        ..offset + target_capture.end(),
                                    attribute: attribute_capture.as_str().to_owned(),
                                }),
                                rustdoc,
                            },
                            &documented_revision,
                            &mut failures,
                            &patterns,
                            (
                                &mut linked_pages,
                                &mut link_rewrites,
                                &mut local_paths,
                                snapshot_tag,
                            ),
                        )?;
                    }
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
                    self.validate_generated_target(
                        output,
                        LinkReference {
                            source: &source,
                            href: &search.href,
                            attribute: None,
                            rustdoc: false,
                        },
                        &documented_revision,
                        &mut failures,
                        &patterns,
                        (
                            &mut linked_pages,
                            &mut link_rewrites,
                            &mut local_paths,
                            snapshot_tag,
                        ),
                    )?;
                }
            }
        }
        for (source, mut replacements) in link_rewrites {
            let mut html = fs::read_to_string(&source)?;
            replacements.sort_by_key(|rewrite| rewrite.range.start);
            ensure!(
                replacements
                    .windows(2)
                    .all(|pair| pair[0].range.end <= pair[1].range.start),
                "overlapping link repairs in {}",
                source.display()
            );
            for rewrite in replacements.iter().rev() {
                ensure!(
                    html.get(rewrite.range.clone()) == Some(rewrite.expected.as_str()),
                    "generated link repair no longer matches {}",
                    source.display()
                );
                html.replace_range(rewrite.range.clone(), &rewrite.replacement);
            }
            let rendered_html = patterns.rendered_html(&html);
            let repaired_targets = patterns
                .tag
                .find_iter(&rendered_html)
                .flat_map(|element| patterns.href.captures_iter(element.as_str()))
                .filter_map(|capture| {
                    (1..=3)
                        .find_map(|index| capture.get(index))
                        .map(|target| target.as_str())
                })
                .collect::<HashSet<_>>();
            for rewrite in &replacements {
                ensure!(
                    !repaired_targets.contains(rewrite.target_before.as_str()),
                    "generated link repair left {} in {}",
                    rewrite.target_before,
                    source.display()
                );
                if let Some(target) = &rewrite.target_after {
                    ensure!(
                        repaired_targets.contains(target.as_str()),
                        "generated link repair did not emit {target} in {}",
                        source.display()
                    );
                }
            }
            fs::write(source, html)?;
        }
        ensure!(
            failures.is_empty(),
            "generated site has broken links:\n{}",
            failures.join("\n")
        );
        Ok(())
    }

    fn validate_generated_target(
        &self,
        output: &Path,
        link: LinkReference<'_>,
        documented_revision: &str,
        failures: &mut Vec<String>,
        patterns: &LinkValidationPatterns,
        state: LinkValidationState<'_>,
    ) -> Result<()> {
        let source_display = link
            .source
            .strip_prefix(output)
            .unwrap_or(link.source)
            .display();
        match repository_source_path(link.href) {
            Ok(Some((revision, _, _))) if revision != documented_revision => {
                failures.push(format!(
                    "{source_display} -> {} (repository source revision {revision} does not match documented revision {documented_revision})",
                    link.href
                ));
                Ok(())
            }
            Ok(Some((_, path, true))) if !self.root.join(&path).is_file() => {
                failures.push(format!(
                    "{source_display} -> {} (missing repository source file {})",
                    link.href,
                    path.display()
                ));
                Ok(())
            }
            Ok(Some((_, path, false))) if !self.root.join(&path).is_dir() => {
                failures.push(format!(
                    "{source_display} -> {} (missing repository source directory {})",
                    link.href,
                    path.display()
                ));
                Ok(())
            }
            Ok(Some(_)) => Ok(()),
            Ok(None) => validate_local_target(output, link, failures, patterns, state),
            Err(error) => {
                failures.push(format!(
                    "{source_display} -> {} (invalid repository source: {error})",
                    link.href
                ));
                Ok(())
            }
        }
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

fn documentation_issue_url(title: &str, route: &str, commit: &str) -> String {
    let title = format!("[docs] {title}");
    let route = route.trim_start_matches('/');
    let body = format!(
        "Documentation page: {PUBLISHED_DOCS_ROOT}/{route}\nDocumented revision: {commit}\n\nWhat should be corrected or clarified?"
    );
    format!(
        "https://github.com/alphal00p/gammaloop/issues/new?labels=documentation&title={}&body={}",
        utf8_percent_encode(&title, NON_ALPHANUMERIC),
        utf8_percent_encode(&body, NON_ALPHANUMERIC),
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
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\"><meta name=\"theme-color\" content=\"#f9f6f0\">{favicon}<title>{} · αLoop</title><link rel=\"stylesheet\" href=\"../assets/site.css\"><script defer src=\"../assets/site.js\"></script></head><body class=\"portal-body portal-subpage-body\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"portal-header\"><a class=\"portal-brand\" href=\"../\" aria-label=\"αLoop home\"><span class=\"portal-brand-logo\" aria-hidden=\"true\"></span><span class=\"portal-brand-copy\"><strong>αLoop</strong><small>Local Unitarity research</small></span></a><nav class=\"portal-nav\" aria-label=\"Primary\"><a href=\"../#projects\">Projects</a><a href=\"../about/\"{}>About</a><a href=\"../people/\"{}>People</a><a href=\"../talks/\"{}>Talks</a><a href=\"../publications/\"{}>Publications</a><a href=\"../developers/\">Developers</a></nav><div class=\"portal-header-actions\"><a class=\"portal-source-link\" href=\"https://github.com/alphal00p/gammaloop\">GitHub <span aria-hidden=\"true\">↗</span></a><button class=\"portal-theme-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\"><span aria-hidden=\"true\">◐</span></button></div></header><main class=\"portal-main portal-subpage-main\" id=\"main-content\">{body}</main><footer class=\"portal-footer\"><div><span class=\"portal-footer-mark\" aria-hidden=\"true\"></span><p><strong>αLoop</strong><br>Local Unitarity research software</p></div><nav aria-label=\"Footer\"><a href=\"../about/\">About</a><a href=\"../#projects\">Projects</a><a href=\"../people/\">People</a><a href=\"../talks/\">Talks</a><a href=\"../publications/\">Publications</a><a href=\"../citations/\">Cite</a><a href=\"../developers/\">Developers</a><a href=\"https://github.com/alphal00p/gammaloop\">Source</a></nav><p>Physics, algorithms, and software<br>developed in the open.</p></footer></body></html>",
        escape_html(description),
        escape_html(title),
        current("about"),
        current("people"),
        current("talks"),
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

fn repository_source_path(href: &str) -> Result<Option<(String, PathBuf, bool)>> {
    let target = href.split(['?', '#']).next().unwrap_or_default();
    let source = [
        ("https://github.com/alphal00p/gammaloop/blob/", true),
        ("https://github.com/alphal00p/gammaloop/tree/", false),
        (
            "https://raw.githubusercontent.com/alphal00p/gammaloop/",
            true,
        ),
    ]
    .into_iter()
    .find_map(|(prefix, file)| target.strip_prefix(prefix).map(|path| (path, file)));
    let Some((remainder, file)) = source else {
        return Ok(None);
    };
    let (revision, path) = remainder.split_once('/').unwrap_or((remainder, ""));
    ensure!(!revision.is_empty(), "missing repository revision");
    let path = normalize_repository_path(Path::new(path))
        .context("repository source path escapes the repository root")?;
    ensure!(
        !file || !path.as_os_str().is_empty(),
        "repository source file path is empty"
    );
    Ok(Some((revision.to_owned(), path, file)))
}

fn validate_developer_typst_body(source: &str) -> Result<()> {
    let active_element =
        Regex::new(r"(?is)<\s*(script|iframe|object|embed|form|base|link|meta|style)\b")?;
    let event_handler = Regex::new(r"(?is)[\s/]on[a-z0-9_-]+\s*=")?;
    let tag = Regex::new(r"(?is)<\s*(?P<tag>[a-z][a-z0-9-]*)\b(?P<attrs>[^>]*)>")?;
    let url_assignment = Regex::new(r"(?is)\b(?:href|src)\s*=")?;
    let quoted_url = Regex::new(r#"(?is)\b(?P<name>href|src)\s*=\s*\"(?P<target>[^\"]*)\""#)?;
    let style = Regex::new(r#"(?is)\bstyle\s*=\s*\"(?P<value>[^\"]*)\""#)?;
    let active_css = Regex::new(r"(?is)(?:@import|url\s*\()")?;
    let character_reference = Regex::new(r"(?is)&(?:#[x]?[0-9a-f]+|[a-z][a-z0-9]+);")?;
    ensure!(
        !active_element.is_match(source) && !event_handler.is_match(source),
        "compiled developer Typst cannot contain active HTML elements or event handlers"
    );
    for tag_capture in tag.captures_iter(source) {
        let element = tag_capture.name("tag").expect("HTML tag capture").as_str();
        let attributes = tag_capture
            .name("attrs")
            .expect("HTML attribute capture")
            .as_str();
        ensure!(
            url_assignment.find_iter(attributes).count()
                == quoted_url.captures_iter(attributes).count(),
            "compiled developer Typst href and src attributes must use double-quoted values"
        );
        for capture in style.captures_iter(attributes) {
            ensure!(
                !active_css.is_match(capture.name("value").expect("style value capture").as_str()),
                "compiled developer Typst cannot contain active inline CSS"
            );
        }
        for capture in quoted_url.captures_iter(attributes) {
            let name = capture.name("name").expect("URL name capture").as_str();
            let target = capture.name("target").expect("URL target capture").as_str();
            ensure!(
                !target.chars().any(char::is_control),
                "compiled developer Typst URLs cannot contain control characters"
            );
            let target = target.trim_start_matches(|character: char| character <= ' ');
            let prefix_end = target
                .find(':')
                .or_else(|| target.find(['/', '?', '#']))
                .unwrap_or(target.len());
            ensure!(
                !character_reference.is_match(&target[..prefix_end]),
                "compiled developer Typst URL schemes cannot contain character references"
            );
            if let Some(data) = target.strip_prefix("data:") {
                ensure!(
                    element.eq_ignore_ascii_case("img")
                        && name.eq_ignore_ascii_case("src")
                        && [
                            "image/png",
                            "image/jpeg",
                            "image/gif",
                            "image/webp",
                            "image/avif"
                        ]
                        .iter()
                        .any(|kind| {
                            data.strip_prefix(kind)
                                .is_some_and(|data| data.starts_with(";base64,"))
                        }),
                    "compiled developer Typst permits only base64 raster data images"
                );
                continue;
            }
            if let Some((scheme, _)) = target.split_once(':')
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
                    "compiled developer Typst uses a forbidden URL scheme: {scheme}"
                );
            }
        }
    }
    Ok(())
}

fn extract_typst_head_styles(source: &str) -> Result<String> {
    let head_start = source
        .find("<head>")
        .map(|start| start + "<head>".len())
        .context("rendered Typst page has no head")?;
    let head_end = source
        .find("</head>")
        .context("rendered Typst page has no head end")?;
    ensure!(head_start <= head_end, "rendered Typst head is malformed");
    let style = Regex::new(r"(?is)<style>(?P<css>.*?)</style>")?;
    let active_css = Regex::new(r"(?is)(?:@import|url\s*\()")?;
    let mut styles = String::new();
    for capture in style.captures_iter(&source[head_start..head_end]) {
        ensure!(
            !active_css.is_match(capture.name("css").expect("Typst style capture").as_str()),
            "compiled developer Typst head contains active CSS"
        );
        styles.push_str(capture.get(0).expect("Typst style element").as_str());
    }
    Ok(styles)
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

fn developer_search_classification(note: &DeveloperNote) -> String {
    let lifecycle = note.lifecycle.to_lowercase();
    let status = note.status.to_lowercase();
    if status == lifecycle
        || status.starts_with(&format!("{lifecycle} "))
        || status.starts_with(&format!("{lifecycle} ·"))
    {
        status
    } else {
        format!("{lifecycle} {status}")
    }
}

fn developer_heading_search_entry(
    note: &DeveloperNote,
    route: &str,
    heading: HeadingLink,
) -> SearchEntry {
    SearchEntry {
        title: heading.title,
        summary: note.title.clone(),
        href: format!("{route}#{}", heading.id),
        kind: format!(
            "developer {} heading",
            developer_search_classification(note)
        ),
        text: String::new(),
    }
}

fn developer_document(
    metadata: (&str, &str),
    docs_root: &str,
    portal_root: &str,
    extra_head: &str,
    sidebar: &str,
    main: &str,
    toc: &str,
) -> String {
    let (title, description) = metadata;
    let home = if docs_root.is_empty() {
        "./"
    } else {
        docs_root
    };
    let favicon = favicon_links(&format!("{docs_root}assets/"));
    format!(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><meta name=\"description\" content=\"{}\"><meta name=\"theme-color\" content=\"#f9f6f0\">{favicon}<title>{}</title><link rel=\"stylesheet\" href=\"{}assets/site.css\">{extra_head}<script defer src=\"{}assets/site.js\"></script></head><body class=\"developer-body\" data-docs-root=\"{}\" data-search-index=\"{}search-index.json\" data-search-root=\"{}\"><a class=\"skip-link\" href=\"#main-content\">Skip to content</a><header class=\"site-header\"><button class=\"header-button menu-button\" type=\"button\" data-menu-toggle aria-label=\"Open navigation\" aria-controls=\"docs-sidebar\" aria-expanded=\"false\">☰</button><a class=\"site-brand\" href=\"{}\"><span class=\"site-brand-mark\">α</span><span class=\"site-brand-name\">Developer notes</span></a><div class=\"site-header-tools\"><a class=\"header-button header-link\" href=\"{}\" aria-label=\"Research documentation\"><span class=\"header-link-label\">Research docs</span><span aria-hidden=\"true\">←</span></a><button class=\"header-button\" type=\"button\" data-search-open aria-label=\"Search all documentation\"><span class=\"header-search-label\">Search</span> <span class=\"header-button-label\">⌘K</span><span class=\"header-search-icon\" aria-hidden=\"true\">⌕</span></button><button class=\"header-button\" type=\"button\" data-theme-toggle aria-label=\"Toggle color theme\">◐</button></div></header><div class=\"docs-shell\">{sidebar}<main class=\"docs-main developer-main\" id=\"main-content\">{main}</main>{toc}</div><button class=\"sidebar-backdrop\" type=\"button\" data-sidebar-backdrop aria-label=\"Close navigation\"></button><dialog class=\"search-dialog\" data-search-dialog><form class=\"search-form\" method=\"dialog\"><input class=\"search-input\" type=\"search\" data-search-input placeholder=\"Search all projects and developer notes\" aria-label=\"Search all documentation\"><button class=\"header-button\" value=\"close\">Close</button></form><ul class=\"search-results\" data-search-results aria-live=\"polite\"></ul></dialog></body></html>",
        escape_html(description),
        escape_html(title),
        escape_html(docs_root),
        escape_html(docs_root),
        escape_html(docs_root),
        escape_html(portal_root),
        escape_html(portal_root),
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

fn leading_developer_title(body: &str) -> Result<LeadingDeveloperTitle> {
    let heading = Regex::new(r"(?s)<h2(?:\s[^>]*)?>(?P<body>.*?)</h2>")?;
    let captures = heading
        .captures(body)
        .context("developer Typst note must begin with one level-one title")?;
    let whole = captures.get(0).expect("developer title capture");
    ensure!(
        body[..whole.start()].trim().is_empty() && !heading.is_match(&body[whole.end()..]),
        "developer Typst note must contain exactly one leading level-one title"
    );
    let tags = Regex::new(r"<[^>]+>")?;
    let title = decode_html_text(
        tags.replace_all(
            captures
                .name("body")
                .expect("developer title body capture")
                .as_str(),
            "",
        )
        .trim(),
    );
    Ok(LeadingDeveloperTitle {
        title,
        range: whole.start()..whole.end(),
    })
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
    let quickstart = product
        .pages
        .iter()
        .find(|page| page.id == "quickstart")
        .expect("product quickstart was validated before rendering");
    let first_guide = product
        .pages
        .iter()
        .find(|page| page.group == "Guides")
        .map(|page| page.route.as_str())
        .unwrap_or("guides/");
    let logo = product.logo.as_deref().map_or_else(String::new, |logo| {
        format!(
            "<span class=\"product-hero-logo product-logo-{}\" role=\"img\" aria-label=\"{} logo\"></span>",
            escape_html(logo),
            escape_html(&product.title),
        )
    });
    format!(
        "<header class=\"product-hero\">{logo}<p class=\"product-eyebrow\">Research software documentation</p><h1>{}</h1><p>{}</p><div class=\"hero-actions\"><a class=\"hero-action primary\" href=\"{}\">Get started</a><a class=\"hero-action\" href=\"{}\">Read the guides</a><a class=\"hero-action\" href=\"reference/\">Browse reference</a></div></header>",
        escape_html(&product.title),
        escape_html(&product.tagline),
        escape_html(&quickstart.route),
        escape_html(first_guide),
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

fn rendered_article_search_text(html: &str) -> Result<String> {
    const ARTICLE: &str = "<article class=\"docs-article\">";
    let start = html
        .find(ARTICLE)
        .map(|start| start + ARTICLE.len())
        .context("rendered page has no documentation article")?;
    let end = html
        .rfind("</article>")
        .context("rendered page has no documentation article end")?;
    ensure!(start <= end, "rendered documentation article is malformed");
    let ignored = Regex::new(
        r"(?is)<nav\b[^>]*>.*?</nav>|<pre\b[^>]*>.*?</pre>|<script\b[^>]*>.*?</script>|<style\b[^>]*>.*?</style>",
    )?;
    let tags = Regex::new(r"(?is)<[^>]+>")?;
    let without_noise = ignored.replace_all(&html[start..end], " ");
    Ok(decode_html_text(&tags.replace_all(&without_noise, " "))
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" "))
}

fn append_rendered_page_search(
    product: &ProductConfig,
    page: &PageConfig,
    page_href: &str,
    html: &str,
    entries: &mut Vec<SearchEntry>,
) -> Result<()> {
    const ARTICLE: &str = "<article class=\"docs-article\">";
    let body = if let Some(start) = html.find(ARTICLE).map(|start| start + ARTICLE.len()) {
        let end = html
            .rfind("</article>")
            .context("rendered page has no documentation article end")?;
        ensure!(start <= end, "rendered documentation article is malformed");
        &html[start..end]
    } else {
        extract_html_body(html)?
    };
    let (_, headings) = inject_heading_ids(body)?;
    let typst_scope = page
        .route
        .strip_prefix("reference/typst/")
        .map(|route| route.trim_matches('/'))
        .filter(|route| !route.is_empty() && !route.contains('/'));
    let module_scope = typst_scope.filter(|scope| matches!(*scope, "graph" | "subgraph"));
    let mut in_generated_reference = false;
    for heading in headings {
        if heading.level == 2 && module_scope == Some(heading.title.as_str()) {
            in_generated_reference = true;
        }
        if heading.level < 2
            || heading.title == page.title
            || page.group == "Typst API" && (heading.level > 3 || heading.title == "Parameters")
        {
            continue;
        }
        let title = if heading.level == 3 && heading.id.contains('.') {
            heading.id.clone()
        } else if let Some(scope) =
            module_scope.filter(|_| in_generated_reference && heading.level == 3)
        {
            format!("{scope}.{}", heading.title)
        } else {
            heading.title
        };
        entries.push(SearchEntry {
            title,
            summary: format!("{} · {}", product.title, page.title),
            href: format!("{page_href}#{}", heading.id),
            kind: page.group.to_lowercase(),
            text: String::new(),
        });
    }
    Ok(())
}

fn render_page_navigation(
    previous: Option<&SitePage>,
    next: Option<&SitePage>,
    docs_root: &str,
) -> String {
    if previous.is_none() && next.is_none() {
        return String::new();
    }
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
    format!(
        "<nav class=\"page-nav\" aria-label=\"Documentation pagination\">{previous}{next}</nav>"
    )
}

fn count_catalog_surface_items(scope: &DocScope, supported: bool) -> usize {
    scope
        .items
        .values()
        .filter(|item| item.supported == supported)
        .count()
        + scope
            .scopes
            .values()
            .map(|scope| count_catalog_surface_items(scope, supported))
            .sum::<usize>()
}

fn catalog_surface_items(scope: &DocScope, supported: bool) -> Vec<(Vec<String>, &DocItem)> {
    fn collect<'a>(
        scope: &'a DocScope,
        path: &mut Vec<String>,
        supported: bool,
        items: &mut Vec<(Vec<String>, &'a DocItem)>,
    ) {
        items.extend(
            scope
                .items
                .values()
                .filter(|item| item.supported == supported)
                .map(|item| (path.clone(), item)),
        );
        for child in scope.scopes.values() {
            path.push(child.id.clone());
            collect(child, path, supported, items);
            path.pop();
        }
    }

    let mut items = Vec::new();
    collect(scope, &mut Vec::new(), supported, &mut items);
    items
}

fn python_item_routes(catalog: &DocCatalog) -> BTreeMap<String, String> {
    let items = catalog_surface_items(&catalog.root, true);
    items
        .into_iter()
        .map(|(path, item)| {
            let anchor = api_item_anchor(&path, item);
            let mut segments = path
                .iter()
                .filter(|segment| segment.as_str() != "exports")
                .cloned()
                .collect::<Vec<_>>();
            segments.push(python_item_disambiguated_segment(item, &item.name));
            (
                anchor,
                format!(
                    "reference/python/{}/{}/",
                    catalog.component.id,
                    segments.join("/")
                ),
            )
        })
        .collect()
}

fn python_item_links(
    catalog: &DocCatalog,
    routes: &BTreeMap<String, String>,
) -> BTreeMap<String, String> {
    let mut candidates = BTreeMap::<String, Vec<String>>::new();
    for (path, item) in catalog_surface_items(&catalog.root, true) {
        candidates
            .entry(item.name.clone())
            .or_default()
            .push(routes[&api_item_anchor(&path, item)].clone());
    }
    candidates
        .into_iter()
        .filter_map(|(name, routes)| {
            (routes.len() == 1).then(|| (name, routes.into_iter().next().unwrap()))
        })
        .collect()
}

fn api_item_kind_label(kind: alphal00p_docs_schema::DocItemKind) -> &'static str {
    match kind {
        alphal00p_docs_schema::DocItemKind::Function
        | alphal00p_docs_schema::DocItemKind::PythonFunction => "Function",
        alphal00p_docs_schema::DocItemKind::Type => "Type",
        alphal00p_docs_schema::DocItemKind::Trait => "Trait",
        alphal00p_docs_schema::DocItemKind::ExportedMacro => "Macro",
        alphal00p_docs_schema::DocItemKind::Method => "Method",
        alphal00p_docs_schema::DocItemKind::Field => "Field",
        alphal00p_docs_schema::DocItemKind::Variant => "Variant",
        alphal00p_docs_schema::DocItemKind::Command => "Command",
        alphal00p_docs_schema::DocItemKind::Setting => "Setting",
        alphal00p_docs_schema::DocItemKind::PythonClass => "Class",
        alphal00p_docs_schema::DocItemKind::PythonConstant => "Constant",
        alphal00p_docs_schema::DocItemKind::ContentPage => "Page",
    }
}

fn link_python_symbols(html: &str, links: &BTreeMap<String, String>) -> String {
    if links.is_empty() {
        return html.to_owned();
    }

    let mut rendered = String::new();
    let mut cursor = 0;
    let mut anchor_depth = 0_u32;
    let mut pre_depth = 0_u32;
    while cursor < html.len() {
        if html[cursor..].starts_with('<') {
            let Some(offset) = html[cursor..].find('>') else {
                rendered.push_str(&html[cursor..]);
                break;
            };
            let end = cursor + offset + 1;
            let tag = &html[cursor..end];
            if tag == "</a>" {
                anchor_depth = anchor_depth.saturating_sub(1);
            } else if tag.starts_with("<a ") || tag == "<a>" {
                anchor_depth += 1;
            } else if tag == "</pre>" {
                pre_depth = pre_depth.saturating_sub(1);
            } else if tag.starts_with("<pre ") || tag == "<pre>" {
                pre_depth += 1;
            }
            rendered.push_str(tag);
            cursor = end;
            continue;
        }
        if html[cursor..].starts_with('&')
            && let Some(offset) = html[cursor..].find(';')
        {
            let end = cursor + offset + 1;
            rendered.push_str(&html[cursor..end]);
            cursor = end;
            continue;
        }
        if anchor_depth == 0 && pre_depth == 0 {
            let token_end = html[cursor..]
                .char_indices()
                .find(|(_, character)| character.is_whitespace() || *character == '<')
                .map_or(html.len(), |(offset, _)| cursor + offset);
            let token = &html[cursor..token_end];
            let uri_scheme = token.split_once(':').is_some_and(|(scheme, _)| {
                scheme
                    .chars()
                    .next()
                    .is_some_and(|character| character.is_ascii_alphabetic())
                    && scheme.chars().all(|character| {
                        character.is_ascii_alphanumeric() || matches!(character, '+' | '-' | '.')
                    })
            });
            if uri_scheme
                || token
                    .get(..4)
                    .is_some_and(|prefix| prefix.eq_ignore_ascii_case("www."))
                || token.contains('/')
            {
                rendered.push_str(token);
                cursor = token_end;
                continue;
            }
        }

        let character = html[cursor..].chars().next().expect("valid HTML cursor");
        if anchor_depth == 0
            && pre_depth == 0
            && (character == '_' || character.is_ascii_alphabetic())
        {
            let end = html[cursor..]
                .char_indices()
                .take_while(|(_, character)| *character == '_' || character.is_ascii_alphanumeric())
                .last()
                .map_or(cursor + character.len_utf8(), |(offset, character)| {
                    cursor + offset + character.len_utf8()
                });
            let identifier = &html[cursor..end];
            if let Some(target) = links.get(identifier) {
                rendered.push_str(&format!(
                    "<a href=\"{}\">{identifier}</a>",
                    escape_html(target),
                ));
            } else {
                rendered.push_str(identifier);
            }
            cursor = end;
            continue;
        }
        rendered.push(character);
        cursor += character.len_utf8();
    }
    rendered
}

fn render_signature_code(signature: &str, links: &BTreeMap<String, String>) -> String {
    link_python_symbols(&escape_html(signature), links)
}

fn python_human_callable_signature(signature: &str) -> String {
    let signature = signature.trim().trim_end_matches(':').trim_end();
    let (prefix, callable) = signature
        .strip_prefix("async def ")
        .map_or(("", signature), |signature| ("async ", signature));
    let callable = callable.strip_prefix("def ").unwrap_or(callable);
    let Some(opening) = callable.find('(') else {
        return format!("{prefix}{callable}");
    };
    let Some(closing) = matching_python_parenthesis(callable, opening) else {
        return format!("{prefix}{callable}");
    };
    let parameters = callable[opening + 1..closing].trim();
    let parameters = ["self", "cls"].into_iter().find_map(|receiver| {
        (parameters == receiver).then_some("").or_else(|| {
            parameters
                .strip_prefix(&format!("{receiver},"))
                .map(str::trim_start)
        })
    });
    format!(
        "{prefix}{}({}){}",
        &callable[..opening],
        parameters.unwrap_or(&callable[opening + 1..closing]),
        &callable[closing + 1..],
    )
}

fn matching_python_parenthesis(signature: &str, opening: usize) -> Option<usize> {
    let mut depth = 0_u32;
    let mut quote = None;
    let mut escaped = false;
    for (offset, character) in signature[opening..].char_indices() {
        if let Some(active_quote) = quote {
            if escaped {
                escaped = false;
            } else if character == '\\' {
                escaped = true;
            } else if character == active_quote {
                quote = None;
            }
            continue;
        }
        if matches!(character, '\'' | '"') {
            quote = Some(character);
        } else if character == '(' {
            depth += 1;
        } else if character == ')' {
            depth -= 1;
            if depth == 0 {
                return Some(opening + offset);
            }
        }
    }
    None
}

fn python_constructor(item: &DocItem) -> Option<&DocMember> {
    (item.kind == alphal00p_docs_schema::DocItemKind::PythonClass)
        .then(|| {
            item.members.iter().find(|member| {
                member.name == "__new__" && member.kind == DocMemberKind::AssociatedFunction
            })
        })
        .flatten()
}

fn python_member_callable_signature(member: &DocMember) -> Option<&str> {
    member.signature.as_deref().or_else(|| {
        member
            .members
            .iter()
            .filter(|member| member.kind == DocMemberKind::Overload)
            .find_map(|member| member.signature.as_deref())
    })
}

fn python_item_display_signature(item: &DocItem) -> Option<String> {
    if item.kind != alphal00p_docs_schema::DocItemKind::PythonClass {
        return item
            .signature
            .as_deref()
            .map(python_human_callable_signature);
    }
    let Some(signature) = python_constructor(item).and_then(python_member_callable_signature)
    else {
        return Some(format!("{}()", item.name));
    };
    let signature = python_human_callable_signature(signature);
    let opening = signature.find('(')?;
    let closing = matching_python_parenthesis(&signature, opening)?;
    Some(format!("{}{}", item.name, &signature[opening..=closing]))
}

fn python_property_display_signature(member: &DocMember) -> String {
    let annotation = match member.kind {
        DocMemberKind::Getter => member.signature.as_deref().and_then(|signature| {
            let signature = python_human_callable_signature(signature);
            let closing = matching_python_parenthesis(&signature, signature.find('(')?)?;
            signature[closing + 1..]
                .trim()
                .strip_prefix("->")
                .map(str::trim)
                .filter(|annotation| !annotation.is_empty())
                .map(str::to_owned)
        }),
        DocMemberKind::Setter => member
            .members
            .iter()
            .find(|member| member.kind == DocMemberKind::Parameter)
            .and_then(|member| member.signature.clone()),
        _ => None,
    };
    annotation.map_or_else(
        || member.name.clone(),
        |annotation| format!("{}: {annotation}", member.name),
    )
}

fn python_member_display_signature(member: &DocMember) -> Option<String> {
    match member.kind {
        DocMemberKind::Getter | DocMemberKind::Setter => {
            Some(python_property_display_signature(member))
        }
        DocMemberKind::Method | DocMemberKind::AssociatedFunction | DocMemberKind::Overload => {
            member
                .signature
                .as_deref()
                .map(python_human_callable_signature)
        }
        _ => member.signature.clone(),
    }
}

fn python_associated_function_label(member: &DocMember) -> &'static str {
    let Some(signature) = python_member_callable_signature(member) else {
        return "Class or static method";
    };
    let Some(opening) = signature.find('(') else {
        return "Class or static method";
    };
    let Some(closing) = matching_python_parenthesis(signature, opening) else {
        return "Class or static method";
    };
    let parameters = signature[opening + 1..closing].trim_start();
    if parameters == "cls" || parameters.starts_with("cls,") {
        "Class method"
    } else {
        "Static method"
    }
}

fn render_python_catalog_for_module(
    catalog: &DocCatalog,
    module: &str,
    stub_name: &str,
    git_commit: &str,
) -> String {
    let supported = catalog_surface_items(&catalog.root, true);
    let implementation_count = count_catalog_surface_items(&catalog.root, false);
    let routes = python_item_routes(catalog);
    let item_links = python_item_links(catalog, &routes);
    let mut body = format!(
        "<p><code>{}</code> · version {}</p><p>Choose a public class, function, or constant below. Each symbol has its own page with a conventional signature-first layout, flat member sections, stable links, and source-backed details.</p><section class=\"reference-coverage\" aria-labelledby=\"python-inventory-title\"><h2 id=\"python-inventory-title\">API inventory</h2><p>The generated inventory distinguishes the supported module API from names that are merely reachable through the extension.</p><div class=\"reference-coverage-grid\"><article><strong>{}</strong><span>supported exports</span><small>one page per supported class, function, or constant</small></article><article><strong>{}</strong><span>implementation-detail exports</span><small>checked for drift but intentionally omitted from the public reference</small></article></div></section>",
        escape_html(module),
        escape_html(&catalog.component.version),
        supported.len(),
        implementation_count,
    );
    if let Some(summary) = &catalog.root.summary {
        body.push_str(&format!("<p>{}</p>", escape_html(summary)));
    }
    if let Some(docs) = &catalog.root.docs {
        body.push_str(&render_doc_text(docs, 2));
    }
    body.push_str(&format!(
        "<nav class=\"reference-guide-links\" aria-label=\"Related Python guides\"><a href=\"quickstart/python/\">{}</a><a href=\"reference/interfaces/\">Python interface guide</a><a href=\"reference/python/\">All Python modules</a></nav>",
        escape_html(&interface_guide_title(&catalog.product.title, "python")),
    ));
    let filter_id = format!("{}-symbol-filter", slug(&catalog.component.id));
    body.push_str(&format!(
        "<div data-reference-filter-root><div class=\"reference-tools\"><label for=\"{}\">Filter public classes, functions, and constants</label><input id=\"{}\" type=\"search\" data-reference-filter placeholder=\"Try a class, function, or constant name\"><output data-reference-filter-count aria-live=\"polite\"></output></div><div class=\"api-symbol-list\" data-reference-filter-scope>",
        escape_html(&filter_id),
        escape_html(&filter_id),
    ));
    for (path, item) in &supported {
        let anchor = api_item_anchor(path, item);
        let route = &routes[&anchor];
        let kind = api_item_kind_label(item.kind);
        let anchored_members = anchored_api_members(&item.members, &anchor);
        let (_, paired_setters) = paired_property_setters(&anchored_members);
        let member_count = anchored_members
            .iter()
            .enumerate()
            .filter(|(index, (member, _))| {
                member.kind != DocMemberKind::Overload
                    && member.name != "__new__"
                    && !paired_setters.contains(index)
            })
            .count();
        let member_meta = if member_count == 0 {
            String::new()
        } else {
            format!(
                "<span class=\"api-symbol-meta\">{member_count} member{}</span>",
                if member_count == 1 { "" } else { "s" },
            )
        };
        let summary = item.summary.clone().or_else(|| {
            item.docs.as_ref().and_then(|docs| {
                docs.body
                    .lines()
                    .map(str::trim)
                    .find(|line| !line.is_empty())
                    .map(str::to_owned)
            })
        });
        body.push_str(&format!(
            "<a class=\"legacy-reference-anchor\" id=\"{}\" href=\"{}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>",
            escape_html(&anchor),
            escape_html(route),
        ));
        append_python_legacy_member_links(&item.members, &anchor, route, &mut body);
        body.push_str(&format!(
            "<article class=\"api-symbol-card\" data-reference-entry data-reference-search=\"{}\"><div><span class=\"api-kind\">{kind}</span><strong class=\"api-symbol-name\"><a href=\"{}\"><code>{}</code></a></strong>{}</div>{member_meta}</article>",
            escape_html(&format!("{} {} {}", item.name, item.title, kind)),
            escape_html(route),
            escape_html(&item.name),
            summary
                .map(|summary| format!("<p>{}</p>", escape_html(&summary)))
                .unwrap_or_default(),
        ));
    }
    body.push_str("</div></div>");
    let source_build_provenance = if catalog.component.features.is_empty() {
        String::new()
    } else {
        format!(
            "<p><strong>Source-build provenance:</strong> this catalog was generated with {}. These are binding build inputs, not Python call arguments.</p>",
            catalog
                .component
                .features
                .iter()
                .map(|feature| format!("<code>{}</code>", escape_html(feature)))
                .collect::<Vec<_>>()
                .join(", "),
        )
    };
    body.push_str(&format!(
        "<aside class=\"stub-source\" aria-label=\"Type-checker stub\"><strong>Type-checker stub source</strong><p><a href=\"reference/python/{}\" download>Download <code>{}</code></a> or <a href=\"https://github.com/alphal00p/gammaloop/blob/{}/docs/api/python/{}\">view its generated source</a>. The structured reference above is the human-facing documentation.</p>{source_build_provenance}</aside>",
        escape_html(stub_name),
        escape_html(stub_name),
        escape_html(git_commit),
        escape_html(stub_name),
    ));
    link_python_symbols(&body, &item_links)
}

fn append_python_legacy_member_links(
    members: &[DocMember],
    parent_anchor: &str,
    route: &str,
    body: &mut String,
) {
    for (member, anchor) in anchored_api_members(members, parent_anchor) {
        body.push_str(&format!(
            "<a class=\"legacy-reference-anchor\" id=\"{}\" href=\"{}#{}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>",
            escape_html(&anchor),
            escape_html(route),
            escape_html(&anchor),
        ));
        append_python_legacy_member_links(&member.members, &anchor, route, body);
    }
}

type ReferenceSibling<'a> = Option<(&'a str, &'a str)>;

fn render_python_member_index(members: &[DocMember], parent_anchor: &str) -> String {
    let anchored = anchored_api_members(members, parent_anchor);
    let (property_setters, paired_setters) = paired_property_setters(&anchored);
    let mut attributes = Vec::new();
    let mut methods = Vec::new();
    let mut class_methods = Vec::new();
    let mut types = Vec::new();
    let mut values = Vec::new();
    for (index, (member, anchor)) in anchored.iter().enumerate() {
        if member.name == "__new__"
            || member.kind == DocMemberKind::Overload
            || paired_setters.contains(&index)
        {
            continue;
        }
        let setter = property_setters
            .get(&index)
            .map(|setter| &anchored[*setter].0);
        let entry = (member, anchor, setter);
        match member.kind {
            DocMemberKind::Getter | DocMemberKind::Setter | DocMemberKind::Field => {
                attributes.push(entry)
            }
            DocMemberKind::Method => methods.push(entry),
            DocMemberKind::AssociatedFunction => class_methods.push(entry),
            DocMemberKind::AssociatedType => types.push(entry),
            DocMemberKind::Variant | DocMemberKind::AssociatedConst => values.push(entry),
            DocMemberKind::Parameter | DocMemberKind::Overload => {}
        }
    }
    let groups = [
        ("Attributes", attributes),
        ("Methods", methods),
        ("Class and static methods", class_methods),
        ("Types", types),
        ("Values", values),
    ];
    if groups.iter().all(|(_, members)| members.is_empty()) {
        return String::new();
    }
    let mut body = "<nav class=\"reference-jump-nav\" aria-labelledby=\"member-overview-title\"><h2 id=\"member-overview-title\">Member overview</h2><div class=\"reference-jump-groups\">".to_owned();
    for (title, members) in groups {
        if members.is_empty() {
            continue;
        }
        body.push_str(&format!(
            "<section class=\"reference-jump-group\"><h3>{title}</h3><ul class=\"reference-jump-list\">"
        ));
        for (member, anchor, setter) in members {
            let label = match member.kind {
                DocMemberKind::Getter if setter.is_some() => "read/write",
                DocMemberKind::Getter => "read-only",
                DocMemberKind::AssociatedFunction => python_associated_function_label(member),
                DocMemberKind::AssociatedType => "type",
                DocMemberKind::AssociatedConst => "constant",
                DocMemberKind::Variant => "value",
                DocMemberKind::Field => "attribute",
                DocMemberKind::Method => "method",
                _ => "member",
            };
            let summary = member
                .docs
                .as_ref()
                .map(|docs| compact_help_summary(&docs.body))
                .filter(|summary| !summary.is_empty())
                .or_else(|| {
                    setter
                        .and_then(|setter| setter.docs.as_ref())
                        .map(|docs| compact_help_summary(&docs.body))
                        .filter(|summary| !summary.is_empty())
                })
                .map(|summary| {
                    format!(
                        "<span class=\"reference-jump-summary\">{}</span>",
                        escape_html(&summary),
                    )
                })
                .unwrap_or_default();
            body.push_str(&format!(
                "<li><a href=\"#{}\"><span><code>{}</code>{summary}</span><small>{}</small></a></li>",
                escape_html(anchor),
                escape_html(&member.name),
                escape_html(label),
            ));
        }
        body.push_str("</ul></section>");
    }
    body.push_str("</div></nav>");
    body
}

fn render_python_item_page(
    catalog: &DocCatalog,
    module: &str,
    symbol: (&[String], &DocItem),
    stub_name: &str,
    git_commit: &str,
    routes: &BTreeMap<String, String>,
    siblings: (ReferenceSibling<'_>, ReferenceSibling<'_>),
) -> String {
    let (path, item) = symbol;
    let (previous, next) = siblings;
    let anchor = api_item_anchor(path, item);
    let anchored_members = anchored_api_members(&item.members, &anchor);
    let constructor = python_constructor(item).and_then(|_| {
        anchored_members.iter().find(|(member, _)| {
            member.name == "__new__" && member.kind == DocMemberKind::AssociatedFunction
        })
    });
    let current_route = &routes[&anchor];
    let mut item_links = python_item_links(catalog, routes);
    item_links.retain(|_, route| route != current_route);
    let kind = api_item_kind_label(item.kind);
    let module_route = format!("reference/python/{}/", catalog.component.id);
    let (workflow_route, workflow_title) = match catalog.product.id.as_str() {
        "gammaloop" => (
            "guides/events-and-observables/",
            "Events and observables workflow",
        ),
        "linnet" => ("guides/algorithms/", "Graph algorithms workflow"),
        "spenso" => ("guides/python/", "Python tensor workflow"),
        "idenso" => ("guides/algebra/", "Algebra and rewrites workflow"),
        "vakint" => ("guides/evaluation/", "Matching and evaluation workflow"),
        _ => ("quickstart/", "Getting started"),
    };
    let mut body = format!(
        "<span id=\"{}\" aria-hidden=\"true\"></span><p class=\"reference-context\"><a href=\"{}\"><code>{}</code></a> <span aria-hidden=\"true\">/</span> {kind}</p><nav class=\"reference-guide-links\" aria-label=\"Python reference navigation\"><a href=\"{}\">All module symbols</a><a href=\"reference/interfaces/\">Python interface guide</a><a href=\"{}\">{}</a></nav>",
        escape_html(&anchor),
        escape_html(&module_route),
        escape_html(module),
        escape_html(&module_route),
        escape_html(workflow_route),
        escape_html(workflow_title),
    );
    if let Some(summary) = item.summary.as_ref().filter(|_| item.docs.is_none()) {
        body.push_str(&format!(
            "<p class=\"reference-lede\">{}</p>",
            escape_html(summary)
        ));
    }
    if let Some(signature) = python_item_display_signature(item) {
        body.push_str(&format!(
            "<div class=\"reference-copy-row\"><pre class=\"api-signature\" id=\"{}-signature\"><code data-lang=\"python\">{}</code></pre><button type=\"button\" data-copy-target=\"{}-signature\">Copy signature</button></div>",
            escape_html(&anchor),
            render_signature_code(&signature, &item_links),
            escape_html(&anchor),
        ));
    }
    let additional_features = item
        .required_features
        .iter()
        .filter(|feature| !catalog.component.features.contains(*feature))
        .collect::<Vec<_>>();
    if !additional_features.is_empty() {
        body.push_str(&format!(
            "<p class=\"api-feature-requirement\"><strong>Additional source-build feature{}:</strong> This binding is generated when {} {} enabled.</p>",
            if additional_features.len() == 1 {
                ""
            } else {
                "s"
            },
            additional_features
                .iter()
                .map(|feature| format!("<code>{}</code>", escape_html(feature)))
                .collect::<Vec<_>>()
                .join(", "),
            if additional_features.len() == 1 {
                "is"
            } else {
                "are"
            },
        ));
    }
    if let Some(docs) = &item.docs {
        body.push_str(&render_doc_text(docs, 2));
    } else if item.summary.is_none() {
        body.push_str(
            "<p class=\"reference-missing\">No description is available for this export.</p>",
        );
    }
    if let Some((constructor, constructor_anchor)) = constructor {
        body.push_str(&format!(
            "<section class=\"api-constructor\" id=\"{}\"><header class=\"reference-detail-heading\"><h2>Constructor</h2><a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to the constructor\">#</a></header>",
            escape_html(constructor_anchor),
            escape_html(constructor_anchor),
        ));
        if let Some(docs) = &constructor.docs {
            body.push_str(&render_doc_text_omitting_parameters(
                docs,
                3,
                constructor
                    .members
                    .iter()
                    .any(|member| member.kind == DocMemberKind::Parameter),
            ));
        }
        render_member_parameters(&constructor.members, 2, &item_links, &mut body);
        body.push_str("</section>");
    }
    if !item.params.is_empty() {
        body.push_str("<h2>Parameters</h2><div class=\"reference-table-wrap\"><table class=\"api-parameter-table\"><thead><tr><th>Name</th><th>Type</th><th>Default</th><th>Description</th></tr></thead><tbody>");
        for parameter in &item.params {
            let default = match parameter.default.as_deref() {
                Some("...") => "runtime default",
                Some(default) => default,
                None => "—",
            };
            let docs = parameter
                .docs
                .as_ref()
                .map(|docs| render_doc_text(docs, 3))
                .unwrap_or_else(|| {
                    "<span aria-label=\"No description available\">—</span>".to_owned()
                });
            body.push_str(&format!(
                "<tr><td><code>{}</code>{}</td><td><code>{}</code></td><td><code>{}</code></td><td>{}</td></tr>",
                escape_html(&parameter.name),
                if parameter.required {
                    " <span class=\"api-kind\">required</span>"
                } else {
                    ""
                },
                render_signature_code(parameter.ty.as_deref().unwrap_or("—"), &item_links),
                escape_html(default),
                docs,
            ));
        }
        body.push_str("</tbody></table></div>");
    }
    if let Some(returns) = &item.returns {
        body.push_str("<h2>Returns</h2>");
        body.push_str(&render_doc_text(returns, 3));
    }
    for example in &item.examples {
        body.push_str(&format!(
            "<h2>{}</h2><pre><code data-lang=\"{}\">{}</code></pre>",
            escape_html(&example.title),
            escape_html(&example.language),
            escape_html(&example.code),
        ));
    }
    body.push_str(&render_python_member_index(&item.members, &anchor));
    if item
        .members
        .iter()
        .any(|member| member.kind != DocMemberKind::Parameter && member.name != "__new__")
    {
        body.push_str("<h2>Member details</h2><div class=\"api-member-list\">");
        render_api_members_flat(&item.members, &anchor, 2, &item_links, &mut body);
        body.push_str("</div>");
    }
    if let Some(source) = &item.source {
        body.push_str(&format!(
            "<p class=\"api-source-link\"><a href=\"https://github.com/alphal00p/gammaloop/blob/{}/{}#L{}\">View generated signature source: <code>{}:{}</code></a></p>",
            escape_html(git_commit),
            escape_html(&source.file),
            source.line,
            escape_html(&source.file),
            source.line,
        ));
    }
    body.push_str(&format!(
        "<aside class=\"stub-source\" aria-label=\"Type-checker stub\"><strong>Type-checker stub</strong><p><a href=\"reference/python/{}\" download>Download <code>{}</code></a>. The module index links every supported public symbol.</p></aside>",
        escape_html(stub_name),
        escape_html(stub_name),
    ));
    if previous.is_some() || next.is_some() {
        body.push_str("<nav class=\"api-symbol-siblings\" aria-label=\"Adjacent Python symbols\">");
        if let Some((name, route)) = previous {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Previous symbol</span><code>{}</code></a>",
                escape_html(route),
                escape_html(name),
            ));
        }
        if let Some((name, route)) = next {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Next symbol</span><code>{}</code></a>",
                escape_html(route),
                escape_html(name),
            ));
        }
        body.push_str("</nav>");
    }
    link_python_symbols(&body, &item_links)
}

fn render_api_members_flat(
    members: &[DocMember],
    parent_anchor: &str,
    level: usize,
    item_links: &BTreeMap<String, String>,
    body: &mut String,
) {
    let member_kind_label = |member: &DocMember| match member.kind {
        alphal00p_docs_schema::DocMemberKind::Parameter => "Parameter",
        alphal00p_docs_schema::DocMemberKind::Overload => "Overload",
        alphal00p_docs_schema::DocMemberKind::Field => "Field",
        alphal00p_docs_schema::DocMemberKind::Variant => "Variant",
        alphal00p_docs_schema::DocMemberKind::AssociatedFunction => {
            python_associated_function_label(member)
        }
        alphal00p_docs_schema::DocMemberKind::AssociatedType => "Associated type",
        alphal00p_docs_schema::DocMemberKind::AssociatedConst => "Associated constant",
        alphal00p_docs_schema::DocMemberKind::Method => "Method",
        alphal00p_docs_schema::DocMemberKind::Getter => "Property · read-only",
        alphal00p_docs_schema::DocMemberKind::Setter => "Property setter",
    };
    render_api_overloads(members, parent_anchor, level, item_links, body);
    let anchored = anchored_api_members(members, parent_anchor);
    let (property_setters, paired_setters) = paired_property_setters(&anchored);
    for (index, (member, anchor)) in anchored.iter().enumerate() {
        if member.kind == DocMemberKind::Overload
            || member.name == "__new__"
            || paired_setters.contains(&index)
        {
            continue;
        }
        let setter = property_setters
            .get(&index)
            .map(|setter| &anchored[*setter]);
        let heading = (level + 1).clamp(3, 6);
        let kind = if setter.is_some() {
            "Property · read/write"
        } else {
            member_kind_label(member)
        };
        if let Some((_, setter_anchor)) = setter {
            body.push_str(&format!(
                "<div class=\"api-member-anchor-alias\" id=\"{}\">",
                escape_html(setter_anchor),
            ));
        }
        body.push_str(&format!(
            "<section class=\"api-member api-member-flat\" id=\"{}\"><header class=\"api-member-heading\"><h{heading}><code>{}</code></h{heading}><span class=\"api-kind\">{kind}</span></header><div class=\"api-member-body\">",
            escape_html(anchor),
            escape_html(&member.name),
        ));
        body.push_str(&format!(
            "<a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to {}\">#</a>{}",
            escape_html(anchor),
            escape_html(&member.name),
            setter.map_or_else(String::new, |(_, setter_anchor)| format!(
                "<a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to the {} setter\"># setter</a>",
                escape_html(setter_anchor),
                escape_html(&member.name),
            )),
        ));
        if let Some(signature) = python_member_display_signature(member) {
            body.push_str(&format!(
                "<pre class=\"api-member-signature\"><code data-lang=\"python\">{}</code></pre>",
                render_signature_code(&signature, item_links),
            ));
        }
        if let Some(docs) = &member.docs {
            body.push_str(&render_doc_text_omitting_parameters(
                docs,
                (heading + 1).clamp(4, 6),
                member
                    .members
                    .iter()
                    .any(|member| member.kind == DocMemberKind::Parameter),
            ));
        } else if !member
            .members
            .iter()
            .any(|member| member.kind == DocMemberKind::Overload)
        {
            body.push_str(
                "<p class=\"reference-missing\">No description is available for this member.</p>",
            );
        }
        if let Some(default) = &member.default {
            body.push_str(&format!(
                "<p><strong>{}:</strong> <code>{}</code></p>",
                if member.kind == DocMemberKind::Variant {
                    "Value"
                } else {
                    "Default"
                },
                escape_html(default)
            ));
        }
        render_member_parameters(&member.members, heading, item_links, body);
        render_api_members_flat(&member.members, anchor, heading, item_links, body);
        if let Some((setter, setter_anchor)) = setter {
            let setter_heading = (heading + 1).clamp(4, 6);
            body.push_str(&format!(
                "<section class=\"api-property-setter\"><h{setter_heading}>Setter</h{setter_heading}>"
            ));
            if let Some(signature) = python_member_display_signature(setter) {
                body.push_str(&format!(
                    "<pre class=\"api-member-signature api-property-setter-signature\"><code data-lang=\"python\">{}</code></pre>",
                    render_signature_code(&signature, item_links),
                ));
            }
            if let Some(docs) = &setter.docs {
                body.push_str(&render_doc_text_omitting_parameters(
                    docs,
                    (setter_heading + 1).clamp(5, 6),
                    setter
                        .members
                        .iter()
                        .any(|member| member.kind == DocMemberKind::Parameter),
                ));
            }
            if let Some(default) = &setter.default {
                body.push_str(&format!(
                    "<p><strong>Default:</strong> <code>{}</code></p>",
                    escape_html(default)
                ));
            }
            render_member_parameters(&setter.members, setter_heading, item_links, body);
            render_api_members_flat(
                &setter.members,
                setter_anchor,
                setter_heading,
                item_links,
                body,
            );
            body.push_str("</section>");
        }
        body.push_str("</div></section>");
        if setter.is_some() {
            body.push_str("</div>");
        }
    }
}

fn render_api_overloads(
    members: &[DocMember],
    parent_anchor: &str,
    level: usize,
    item_links: &BTreeMap<String, String>,
    body: &mut String,
) {
    let overloads = anchored_api_members(members, parent_anchor)
        .into_iter()
        .filter(|(member, _)| member.kind == DocMemberKind::Overload)
        .collect::<Vec<_>>();
    if overloads.is_empty() {
        return;
    }
    let heading = (level + 1).clamp(3, 6);
    let overload_heading = (heading + 1).clamp(4, 6);
    body.push_str(&format!(
        "<section class=\"api-overloads\"><h{heading}>Overloads</h{heading}>"
    ));
    for (index, (overload, anchor)) in overloads.iter().enumerate() {
        body.push_str(&format!(
            "<article class=\"api-overload\" id=\"{}\"><h{overload_heading}>Overload {} <a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to overload {} of {}\">#</a></h{overload_heading}>",
            escape_html(anchor),
            index + 1,
            escape_html(anchor),
            index + 1,
            escape_html(parent_anchor),
        ));
        if let Some(signature) = python_member_display_signature(overload) {
            body.push_str(&format!(
                "<pre class=\"api-member-signature api-overload-signature\"><code data-lang=\"python\">{}</code></pre>",
                render_signature_code(&signature, item_links),
            ));
        }
        if let Some(docs) = &overload.docs {
            body.push_str(&render_doc_text_omitting_parameters(
                docs,
                6,
                overload
                    .members
                    .iter()
                    .any(|member| member.kind == DocMemberKind::Parameter),
            ));
        } else {
            body.push_str("<p class=\"reference-missing\">No additional description is available for this overload.</p>");
        }
        render_member_parameters(&overload.members, overload_heading, item_links, body);
        body.push_str("</article>");
    }
    body.push_str("</section>");
}

fn anchored_api_members<'a>(
    members: &'a [DocMember],
    parent_anchor: &str,
) -> Vec<(&'a DocMember, String)> {
    let mut anchor_counts = BTreeMap::new();
    members
        .iter()
        .filter(|member| member.kind != DocMemberKind::Parameter)
        .map(|member| {
            (
                member,
                next_api_member_anchor(parent_anchor, member, &mut anchor_counts),
            )
        })
        .collect()
}

fn paired_property_setters(
    members: &[(&DocMember, String)],
) -> (BTreeMap<usize, usize>, BTreeSet<usize>) {
    let mut pairs = BTreeMap::new();
    let mut setters = BTreeSet::new();
    for (getter_index, (getter, _)) in members.iter().enumerate() {
        if getter.kind != DocMemberKind::Getter {
            continue;
        }
        if let Some((setter_index, _)) = members.iter().enumerate().find(|(index, (setter, _))| {
            setter.kind == DocMemberKind::Setter
                && setter.name == getter.name
                && !setters.contains(index)
        }) {
            pairs.insert(getter_index, setter_index);
            setters.insert(setter_index);
        }
    }
    (pairs, setters)
}

fn render_member_parameters(
    members: &[DocMember],
    level: usize,
    item_links: &BTreeMap<String, String>,
    body: &mut String,
) {
    let parameters = members
        .iter()
        .filter(|member| member.kind == DocMemberKind::Parameter)
        .collect::<Vec<_>>();
    if parameters.is_empty() {
        return;
    }
    let heading = (level + 1).clamp(3, 6);
    body.push_str(&format!(
        "<h{heading}>Parameters</h{heading}><div class=\"reference-table-wrap\"><table class=\"api-parameter-table\"><thead><tr><th>Name</th><th>Type</th><th>Default</th><th>Description</th></tr></thead><tbody>"
    ));
    for parameter in parameters {
        let default = match parameter.default.as_deref() {
            Some("...") => "runtime default",
            Some(default) => default,
            None => "—",
        };
        let docs = parameter
            .docs
            .as_ref()
            .map(|docs| render_doc_text(docs, (heading + 1).clamp(5, 6)))
            .unwrap_or_else(|| "<span class=\"reference-missing\">—</span>".to_owned());
        body.push_str(&format!(
            "<tr><td><code>{}</code></td><td><code>{}</code></td><td><code>{}</code></td><td>{}</td></tr>",
            escape_html(&parameter.name),
            render_signature_code(parameter.signature.as_deref().unwrap_or("—"), item_links),
            escape_html(default),
            docs,
        ));
    }
    body.push_str("</tbody></table></div>");
}

fn render_doc_inline(text: &str, format: DocFormat) -> String {
    let mut rendered = escape_html(text);
    match format {
        DocFormat::RustMarkdown => {
            let links = Regex::new(r"\[([^\]]+)\]\(([^)]+)\)").expect("static Markdown link regex");
            rendered = links
                .replace_all(&rendered, |captures: &regex::Captures<'_>| {
                    format!("<a href=\"{}\">{}</a>", &captures[2], &captures[1])
                })
                .into_owned();
            let intra_doc = Regex::new(r"\[`([^`]+)`\]").expect("static intra-doc link regex");
            rendered = intra_doc
                .replace_all(&rendered, "<code>$1</code>")
                .into_owned();
            let code = Regex::new(r"`([^`]+)`").expect("static Markdown code regex");
            rendered = code.replace_all(&rendered, "<code>$1</code>").into_owned();
            let strong = Regex::new(r"\*\*([^*]+)\*\*").expect("static Markdown strong regex");
            strong
                .replace_all(&rendered, "<strong>$1</strong>")
                .into_owned()
        }
        DocFormat::PythonDocstring => {
            let code = Regex::new(r"``([^`]+)``").expect("static Python code regex");
            rendered = code.replace_all(&rendered, "<code>$1</code>").into_owned();
            let code = Regex::new(r"`([^`]+)`").expect("static Python code regex");
            rendered = code.replace_all(&rendered, "<code>$1</code>").into_owned();
            let strong = Regex::new(r"\*\*([^*]+)\*\*").expect("static Python strong regex");
            strong
                .replace_all(&rendered, "<strong>$1</strong>")
                .into_owned()
        }
        DocFormat::TypstMarkup => {
            let links = Regex::new(r#"#link\(&quot;(.*?)&quot;\)\[([^\]]+)\]"#)
                .expect("static Typst link regex");
            rendered = links
                .replace_all(&rendered, |captures: &regex::Captures<'_>| {
                    format!("<a href=\"{}\">{}</a>", &captures[1], &captures[2])
                })
                .into_owned();
            let emphasis = Regex::new(r"#emph\[([^\]]+)\]").expect("static Typst emphasis regex");
            rendered = emphasis.replace_all(&rendered, "<em>$1</em>").into_owned();
            let raw = Regex::new(r#"#raw\(&quot;(.*?)&quot;\)"#).expect("static Typst raw regex");
            raw.replace_all(&rendered, "<code>$1</code>").into_owned()
        }
        DocFormat::PlainText => rendered,
    }
}

fn render_prose_blocks(lines: &[&str], format: DocFormat) -> String {
    let render_unfenced = |lines: &[&str]| {
        lines
            .split(|line| line.trim().is_empty())
            .filter(|paragraph| paragraph.iter().any(|line| !line.trim().is_empty()))
            .map(|paragraph| {
                let first_bullet = paragraph
                    .iter()
                    .position(|line| line.trim_start().starts_with("- "));
                let Some(first_bullet) = first_bullet else {
                    return format!(
                        "<p>{}</p>",
                        render_doc_inline(
                            &paragraph
                                .iter()
                                .map(|line| line.trim())
                                .collect::<Vec<_>>()
                                .join(" "),
                            format,
                        )
                    );
                };
                let mut rendered = String::new();
                if first_bullet > 0 {
                    rendered.push_str(&format!(
                        "<p>{}</p>",
                        render_doc_inline(
                            &paragraph[..first_bullet]
                                .iter()
                                .map(|line| line.trim())
                                .collect::<Vec<_>>()
                                .join(" "),
                            format,
                        )
                    ));
                }
                rendered.push_str("<ul>");
                for item in &paragraph[first_bullet..] {
                    rendered.push_str(&format!(
                        "<li>{}</li>",
                        render_doc_inline(item.trim().trim_start_matches("- "), format),
                    ));
                }
                rendered.push_str("</ul>");
                rendered
            })
            .collect::<String>()
    };
    let mut rendered = String::new();
    let mut prose_start = 0;
    let mut index = 0;
    while index < lines.len() {
        let Some(language) = lines[index].trim().strip_prefix("```") else {
            index += 1;
            continue;
        };
        rendered.push_str(&render_unfenced(&lines[prose_start..index]));
        index += 1;
        let code_start = index;
        while index < lines.len() && !lines[index].trim().starts_with("```") {
            index += 1;
        }
        let code_lines = &lines[code_start..index];
        let code = if format == DocFormat::PythonDocstring
            && language.trim() == "python"
            && code_lines
                .iter()
                .find(|line| !line.trim().is_empty())
                .is_some_and(|line| line.trim_start().starts_with(">>>"))
        {
            code_lines
                .iter()
                .filter_map(|line| {
                    let line = line.trim_start();
                    line.strip_prefix(">>>")
                        .or_else(|| line.strip_prefix("..."))
                        .map(|line| line.strip_prefix(' ').unwrap_or(line))
                })
                .collect::<Vec<_>>()
                .join("\n")
        } else {
            code_lines.join("\n")
        };
        rendered.push_str(&format!(
            "<pre><code data-lang=\"{}\">{}</code></pre>",
            escape_html(language.trim()),
            escape_html(&code),
        ));
        index += usize::from(index < lines.len());
        prose_start = index;
    }
    rendered.push_str(&render_unfenced(&lines[prose_start..]));
    rendered
}

fn render_rust_markdown(body: &str, heading_level: usize) -> String {
    let lines = body.lines().collect::<Vec<_>>();
    let mut rendered = String::new();
    let mut prose = Vec::new();
    let mut index = 0;
    while index < lines.len() {
        let line = lines[index].trim();
        if let Some(language) = line.strip_prefix("```") {
            rendered.push_str(&render_prose_blocks(&prose, DocFormat::RustMarkdown));
            prose.clear();
            index += 1;
            let start = index;
            while index < lines.len() && !lines[index].trim().starts_with("```") {
                index += 1;
            }
            rendered.push_str(&format!(
                "<pre><code data-lang=\"{}\">{}</code></pre>",
                escape_html(language.trim()),
                escape_html(&lines[start..index].join("\n")),
            ));
        } else if line.starts_with('#') {
            let hashes = line
                .chars()
                .take_while(|character| *character == '#')
                .count();
            let title = line[hashes..].trim();
            if !title.is_empty() {
                rendered.push_str(&render_prose_blocks(&prose, DocFormat::RustMarkdown));
                prose.clear();
                let level = (heading_level + hashes.saturating_sub(1)).clamp(2, 6);
                rendered.push_str(&format!(
                    "<h{level}>{}</h{level}>",
                    render_doc_inline(title, DocFormat::RustMarkdown),
                ));
            }
        } else {
            prose.push(lines[index]);
        }
        index += 1;
    }
    rendered.push_str(&render_prose_blocks(&prose, DocFormat::RustMarkdown));
    rendered
}

fn render_python_example(content: &[&str]) -> String {
    let start = content
        .iter()
        .position(|line| !line.trim().is_empty())
        .unwrap_or(content.len());
    let end = content
        .iter()
        .rposition(|line| !line.trim().is_empty())
        .map_or(start, |index| index + 1);
    let content = &content[start..end];
    if content
        .iter()
        .any(|line| line.trim_start().starts_with("```"))
    {
        return render_prose_blocks(content, DocFormat::PythonDocstring);
    }

    let render_code = |lines: &[&str]| {
        let indent = lines
            .iter()
            .filter(|line| !line.trim().is_empty())
            .map(|line| line.len() - line.trim_start().len())
            .min()
            .unwrap_or(0);
        let code = lines
            .iter()
            .map(|line| {
                if line.trim().is_empty() {
                    ""
                } else {
                    &line[indent.min(line.len())..]
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        format!(
            "<pre><code data-lang=\"python\">{}</code></pre>",
            escape_html(code.trim_end()),
        )
    };

    if let Some(marker) = content.iter().enumerate().position(|(index, line)| {
        line.trim_end().ends_with("::")
            && content[index + 1..].iter().any(|candidate| {
                !candidate.trim().is_empty() && candidate.len() > candidate.trim_start().len()
            })
    }) {
        let mut rendered = render_prose_blocks(&content[..marker], DocFormat::PythonDocstring);
        let marker_text = content[marker]
            .trim()
            .strip_suffix("::")
            .unwrap_or(content[marker].trim());
        rendered.push_str(&format!(
            "<p>{}:</p>",
            render_doc_inline(marker_text, DocFormat::PythonDocstring),
        ));
        let code_start = (marker + 1..content.len())
            .find(|index| !content[*index].trim().is_empty())
            .unwrap_or(content.len());
        let code_indent = content
            .get(code_start)
            .map_or(0, |line| line.len() - line.trim_start().len());
        let code_end = (code_start + 1..content.len())
            .find(|index| {
                let line = content[*index];
                !line.trim().is_empty() && line.len() - line.trim_start().len() < code_indent
            })
            .unwrap_or(content.len());
        rendered.push_str(&render_code(&content[code_start..code_end]));
        rendered.push_str(&render_prose_blocks(
            &content[code_end..],
            DocFormat::PythonDocstring,
        ));
        return rendered;
    }

    if content
        .first()
        .is_some_and(|line| line.trim_start().starts_with(">>>"))
    {
        let source = content
            .iter()
            .filter_map(|line| {
                let line = line.trim_start();
                line.strip_prefix(">>>")
                    .or_else(|| line.strip_prefix("..."))
                    .map(|line| line.strip_prefix(' ').unwrap_or(line))
            })
            .collect::<Vec<_>>();
        return render_code(&source);
    }
    if let Some(code_start) = content
        .iter()
        .position(|line| !line.trim().is_empty() && line.len() > line.trim_start().len())
    {
        let mut rendered = render_prose_blocks(&content[..code_start], DocFormat::PythonDocstring);
        rendered.push_str(&render_code(&content[code_start..]));
        return rendered;
    }
    render_code(content)
}

fn render_doc_text(docs: &alphal00p_docs_schema::DocText, heading_level: usize) -> String {
    render_doc_text_omitting_parameters(docs, heading_level, false)
}

fn render_doc_text_omitting_parameters(
    docs: &alphal00p_docs_schema::DocText,
    heading_level: usize,
    omit_parameters: bool,
) -> String {
    if docs.format == DocFormat::RustMarkdown {
        return format!(
            "<div class=\"api-docstring\">{}</div>",
            render_rust_markdown(&docs.body, heading_level),
        );
    }
    if docs.format != DocFormat::PythonDocstring {
        return format!(
            "<div class=\"api-docstring\">{}</div>",
            render_prose_blocks(&docs.body.lines().collect::<Vec<_>>(), docs.format),
        );
    }

    let render_paragraphs = |lines: &[&str]| render_prose_blocks(lines, docs.format);

    let lines = docs.body.lines().collect::<Vec<_>>();
    let mut sections = vec![];
    let mut start = 0;
    let mut title = None;
    let mut numpy_section = false;
    let mut index = 0;
    while index < lines.len() {
        let raw_heading = lines[index].trim();
        let heading = raw_heading
            .trim_start_matches('#')
            .trim()
            .trim_end_matches(':');
        let common_heading = matches!(
            heading.to_ascii_lowercase().as_str(),
            "parameters"
                | "arguments"
                | "other parameters"
                | "keyword arguments"
                | "attributes"
                | "variants"
                | "returns"
                | "yields"
                | "raises"
                | "warns"
                | "examples"
                | "notes"
                | "see also"
        );
        let numpy_heading = index + 1 < lines.len()
            && !heading.is_empty()
            && lines[index + 1].trim().len() >= 3
            && lines[index + 1]
                .trim()
                .chars()
                .all(|character| character == '-');
        let markdown_heading = raw_heading.starts_with('#') && common_heading;
        let labeled_heading = raw_heading.ends_with(':') && common_heading;
        if numpy_heading || markdown_heading || labeled_heading {
            if start < index {
                sections.push((title.take(), &lines[start..index], numpy_section));
            }
            title = Some(heading);
            numpy_section = numpy_heading;
            let consumed = if numpy_heading { 2 } else { 1 };
            start = index + consumed;
            index += consumed;
        } else {
            index += 1;
        }
    }
    if start < lines.len() {
        sections.push((title, &lines[start..], numpy_section));
    } else if let Some(title) = title {
        sections.push((Some(title), &lines[lines.len()..], numpy_section));
    }

    let mut rendered = String::from("<div class=\"api-docstring\">");
    for (title, content, numpy_section) in sections {
        let Some(title) = title else {
            rendered.push_str(&render_paragraphs(content));
            continue;
        };
        if omit_parameters
            && matches!(
                title.to_ascii_lowercase().as_str(),
                "parameters" | "arguments" | "other parameters" | "keyword arguments"
            )
        {
            continue;
        }
        let section_class = slug(title);
        rendered.push_str(&format!(
            "<section class=\"api-doc-section api-doc-{}\"><h{heading_level}>{}</h{heading_level}>",
            escape_html(&section_class),
            escape_html(title),
        ));
        match title.to_ascii_lowercase().as_str() {
            "parameters" | "arguments" | "other parameters" | "keyword arguments"
            | "attributes" | "variants" => {
                let variants = title.eq_ignore_ascii_case("variants");
                let mut definitions = vec![];
                let mut term = None;
                let mut description = vec![];
                for line in content
                    .iter()
                    .map(|line| line.trim())
                    .filter(|line| !line.is_empty())
                {
                    if let Some((name, ty)) = line.split_once(" : ") {
                        if let Some((name, ty)) = term.take() {
                            definitions.push((name, ty, std::mem::take(&mut description)));
                        }
                        term = Some((name.trim(), ty.trim()));
                    } else {
                        description.push(line);
                    }
                }
                if let Some((name, ty)) = term {
                    definitions.push((name, ty, description));
                }
                if definitions.is_empty() {
                    rendered.push_str(&render_paragraphs(content));
                } else {
                    rendered.push_str("<dl class=\"api-doc-definitions\">");
                    for (name, ty, description) in definitions {
                        if variants {
                            let mut variant_description = vec![ty];
                            variant_description.extend(description);
                            rendered.push_str(&format!(
                                "<dt><code>{}</code></dt><dd>{}</dd>",
                                escape_html(name),
                                render_paragraphs(&variant_description),
                            ));
                        } else {
                            rendered.push_str(&format!(
                                "<dt><code>{}</code> <span>{}</span></dt><dd>{}</dd>",
                                escape_html(name),
                                escape_html(ty),
                                render_paragraphs(&description),
                            ));
                        }
                    }
                    rendered.push_str("</dl>");
                }
            }
            "returns" | "yields" => {
                if !numpy_section {
                    rendered.push_str(&render_paragraphs(content));
                    rendered.push_str("</section>");
                    continue;
                }
                let content = content
                    .iter()
                    .map(|line| line.trim())
                    .filter(|line| !line.is_empty())
                    .collect::<Vec<_>>();
                if let Some((name, ty)) = content.first().and_then(|line| line.split_once(" : ")) {
                    rendered.push_str(&format!(
                        "<dl class=\"api-doc-definitions\"><dt><code>{}</code> <span>{}</span></dt><dd>{}</dd></dl>",
                        escape_html(name.trim()),
                        escape_html(ty.trim()),
                        render_paragraphs(&content[1..]),
                    ));
                } else if let Some((ty, description)) = content.split_first() {
                    rendered.push_str(&format!(
                        "<p class=\"api-doc-type\"><code>{}</code></p>{}",
                        escape_html(ty),
                        render_paragraphs(description),
                    ));
                }
            }
            "raises" | "warns" => {
                let mut definitions = vec![];
                let mut exception = None;
                let mut description = vec![];
                for line in content
                    .iter()
                    .map(|line| line.trim())
                    .filter(|line| !line.is_empty())
                {
                    let is_exception = line.split(',').all(|candidate| {
                        let candidate = candidate.trim();
                        candidate.ends_with("Error")
                            || candidate.ends_with("Exception")
                            || candidate.ends_with("Warning")
                    });
                    if is_exception {
                        if let Some(exception) = exception.take() {
                            definitions.push((exception, std::mem::take(&mut description)));
                        }
                        exception = Some(line);
                    } else {
                        description.push(line);
                    }
                }
                if let Some(exception) = exception {
                    definitions.push((exception, description));
                }
                if definitions.is_empty() {
                    rendered.push_str(&render_paragraphs(content));
                } else {
                    rendered.push_str("<dl class=\"api-doc-definitions\">");
                    for (exception, description) in definitions {
                        rendered.push_str(&format!(
                            "<dt><code>{}</code></dt><dd>{}</dd>",
                            escape_html(exception),
                            render_paragraphs(&description),
                        ));
                    }
                    rendered.push_str("</dl>");
                }
            }
            "examples" => rendered.push_str(&render_python_example(content)),
            _ => rendered.push_str(&render_paragraphs(content)),
        }
        rendered.push_str("</section>");
    }
    rendered.push_str("</div>");
    rendered
}

fn python_item_disambiguated_segment(item: &DocItem, value: &str) -> String {
    match item.kind {
        alphal00p_docs_schema::DocItemKind::PythonConstant => format!("{value}-constant"),
        alphal00p_docs_schema::DocItemKind::PythonFunction => format!("{value}-function"),
        _ => value.to_owned(),
    }
}

fn api_item_anchor(path: &[String], item: &DocItem) -> String {
    let mut parts = path.to_vec();
    parts.push(python_item_disambiguated_segment(item, &item.id));
    slug(&parts.join("-"))
}

fn next_api_member_anchor(
    parent_anchor: &str,
    member: &DocMember,
    counts: &mut BTreeMap<String, usize>,
) -> String {
    let base = if member.kind == DocMemberKind::Overload {
        format!("{parent_anchor}-overload")
    } else {
        format!(
            "{parent_anchor}-{}-{}",
            slug(&member.name),
            slug(&format!("{:?}", member.kind))
        )
    };
    let occurrence = counts.entry(base.clone()).or_default();
    *occurrence += 1;
    if *occurrence == 1 {
        base
    } else {
        format!("{base}-{occurrence}")
    }
}

fn append_catalog_search_with_component_name(
    catalog: &DocCatalog,
    scope: &DocScope,
    component_name: &str,
    entries: &mut Vec<SearchEntry>,
) {
    let mut path = Vec::new();
    append_catalog_search_at(catalog, scope, component_name, &mut path, entries);
}

fn append_catalog_search_at(
    catalog: &DocCatalog,
    scope: &DocScope,
    component_name: &str,
    path: &mut Vec<String>,
    entries: &mut Vec<SearchEntry>,
) {
    for item in scope.items.values() {
        if !item.supported {
            continue;
        }
        let (href, kind, member_anchor) = match catalog.component.language {
            ApiLanguage::Rust => (rustdoc_href(catalog, item), "rust-api", None),
            ApiLanguage::Python => {
                let anchor = api_item_anchor(path, item);
                (
                    python_item_routes(catalog)[&anchor].clone(),
                    "python-api",
                    Some(anchor),
                )
            }
            _ => continue,
        };
        entries.push(SearchEntry {
            title: if catalog.component.language == ApiLanguage::Python {
                format!("{component_name}.{}", item.name)
            } else {
                format!("{component_name} · {}", item.title)
            },
            summary: compact_help_summary(item.summary.as_deref().unwrap_or_default()),
            href: href.clone(),
            kind: kind.to_owned(),
            text: [
                item.summary.as_deref(),
                item.docs.as_ref().map(|docs| docs.body.as_str()),
                item.signature.as_deref(),
            ]
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .join(" "),
        });
        if catalog.component.language == ApiLanguage::Python {
            append_member_search(
                component_name,
                &item.name,
                &href,
                kind,
                member_anchor.as_deref(),
                &item.members,
                entries,
            );
        }
    }
    for child in scope.scopes.values() {
        path.push(child.id.clone());
        append_catalog_search_at(catalog, child, component_name, path, entries);
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
    let anchored = parent_anchor
        .map(|anchor| anchored_api_members(members, anchor))
        .unwrap_or_else(|| {
            members
                .iter()
                .filter(|member| member.kind != DocMemberKind::Parameter)
                .map(|member| (member, String::new()))
                .collect()
        });
    let (property_setters, paired_setters) = paired_property_setters(&anchored);
    let mut overload_index = 0;
    for (index, (member, member_anchor)) in anchored.iter().enumerate() {
        if paired_setters.contains(&index) {
            continue;
        }
        let setter = property_setters
            .get(&index)
            .map(|setter| &anchored[*setter]);
        let title = if member.name == "__new__" {
            format!("{parent} constructor")
        } else if member.kind == DocMemberKind::Overload {
            overload_index += 1;
            format!("{parent} · overload {overload_index}")
        } else {
            format!("{parent}.{}", member.name)
        };
        let member_href = if member_anchor.is_empty() {
            href.to_owned()
        } else {
            format!("{}#{member_anchor}", href.split('#').next().unwrap_or(href))
        };
        let setter = setter.map(|(setter, _)| *setter);
        let overload_signatures = member
            .members
            .iter()
            .filter(|member| member.kind == DocMemberKind::Overload)
            .filter_map(|member| member.signature.as_deref())
            .collect::<Vec<_>>()
            .join(" | ");
        let summary = [
            member.docs.as_ref().map(|docs| docs.body.as_str()),
            setter.and_then(|setter| setter.docs.as_ref().map(|docs| docs.body.as_str())),
            member.signature.as_deref(),
            setter.and_then(|setter| setter.signature.as_deref()),
            (!overload_signatures.is_empty()).then_some(overload_signatures.as_str()),
        ]
        .into_iter()
        .flatten()
        .map(compact_help_summary)
        .find(|summary| !summary.is_empty())
        .unwrap_or_default();
        let text = [
            member.docs.as_ref().map(|docs| docs.body.as_str()),
            setter.and_then(|setter| setter.docs.as_ref().map(|docs| docs.body.as_str())),
            member.signature.as_deref(),
            setter.and_then(|setter| setter.signature.as_deref()),
            (!overload_signatures.is_empty()).then_some(overload_signatures.as_str()),
        ]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>()
        .join(" ");
        entries.push(SearchEntry {
            title: format!("{component}.{title}"),
            summary,
            href: member_href.clone(),
            kind: kind.to_owned(),
            text,
        });
        append_member_search(
            component,
            &title,
            &member_href,
            kind,
            (!member_anchor.is_empty()).then_some(member_anchor.as_str()),
            &member.members,
            entries,
        );
        if let Some((setter, setter_anchor)) = property_setters
            .get(&index)
            .map(|setter| &anchored[*setter])
        {
            let setter_href = format!("{}#{setter_anchor}", href.split('#').next().unwrap_or(href));
            append_member_search(
                component,
                &title,
                &setter_href,
                kind,
                Some(setter_anchor),
                &setter.members,
                entries,
            );
        }
    }
}

fn append_scope_typst(
    catalog: &DocCatalog,
    scope: &DocScope,
    path: &mut Vec<String>,
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
            ApiLanguage::Python => {
                python_item_routes(catalog)[&api_item_anchor(path, item)].clone()
            }
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
        path.push(child.id.clone());
        append_scope_typst(catalog, child, path, git_commit, output);
        path.pop();
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

    // Match the declaration keyword immediately preceding this item's name.
    // This remains robust for compact signatures and older serialized catalogs
    // whose signature text may still include attributes or prose.
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

fn ensure_safe_output(root: &Path, output: &Path) -> Result<PathBuf> {
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
    let canonical_root = fs::canonicalize(root)
        .wrap_err_with(|| format!("failed to resolve workspace root {}", root.display()))?;
    let existing = closest_existing_directory(output)
        .context("documentation output has no existing directory ancestor")?;
    let canonical_existing = fs::canonicalize(&existing)
        .wrap_err_with(|| format!("failed to resolve {}", existing.display()))?;
    let resolved_output = canonical_existing.join(output.strip_prefix(&existing)?);
    ensure!(
        resolved_output != canonical_root && !canonical_root.starts_with(&resolved_output),
        "documentation output cannot contain the workspace root"
    );
    Ok(resolved_output)
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
            for command in gamma
                .commands
                .iter()
                .filter(|command| !command.hidden && !command.generated_help)
            {
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
        "<!doctype html><meta charset=\"utf-8\"><meta name=\"robots\" content=\"noindex\"><meta http-equiv=\"refresh\" content=\"0; url={}\"><link rel=\"canonical\" href=\"{}\"><title>{} · {}</title><a href=\"{}\">Open the {} for {}</a>",
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
        "<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width\"><title>{} · {}</title><style>body{{max-width:75rem;margin:auto;padding:2rem 1rem;color:#172033;font:16px/1.5 system-ui}}a{{color:#5b4bdb}}pre{{overflow:auto;padding:1rem;background:#f1f3f7}}table{{width:100%;border-collapse:collapse;margin:1rem 0 2rem}}th,td{{padding:.45rem .6rem;border:1px solid #d9dce5;text-align:left;vertical-align:top}}th{{background:#f1f3f7}}section{{padding-top:.5rem}}</style></head><body><nav><a href=\"../../\">Documentation home</a></nav><h1>{}</h1>{body}</body></html>",
        escape_html(product),
        escape_html(title),
        escape_html(title),
    )
}

fn compact_help_summary(help: &str) -> String {
    const MAX_CHARACTERS: usize = 220;
    let summary = help
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .map(|line| {
            line.split(" # ")
                .next()
                .unwrap_or(line)
                .trim_start_matches('#')
                .trim()
                .trim_start_matches("- ")
                .replace('`', "")
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
        })
        .unwrap_or_default();
    let Some(cut) = summary
        .char_indices()
        .nth(MAX_CHARACTERS)
        .map(|(index, _)| index)
    else {
        return summary;
    };
    let prefix = &summary[..cut];
    let boundary = prefix.rfind(char::is_whitespace).unwrap_or(cut);
    format!(
        "{}…",
        prefix[..boundary].trim_end_matches(['.', ',', ':', ';'])
    )
}

fn render_cli_help(help: &str) -> String {
    help.split("\n\n")
        .filter(|block| !block.trim().is_empty())
        .map(|block| {
            if let Some((lead, items)) = block.trim().split_once(": - ") {
                let items = items
                    .split(" - ")
                    .map(|item| {
                        format!(
                            "<li>{}</li>",
                            render_doc_inline(item.trim(), DocFormat::RustMarkdown)
                        )
                    })
                    .collect::<String>();
                format!(
                    "<p>{}:</p><ul>{items}</ul>",
                    render_doc_inline(lead.trim(), DocFormat::RustMarkdown),
                )
            } else {
                render_rust_markdown(block, 6)
            }
        })
        .collect()
}

fn cli_route_segment(value: &str) -> String {
    if value == "!" {
        return "bang".to_owned();
    }
    if let Some(sequence) = value.strip_suffix("[]") {
        return format!("{}-items", slug(sequence));
    }
    slug(value)
}

fn cli_command_route(path: &str) -> String {
    let segments = path
        .split_whitespace()
        .map(cli_route_segment)
        .collect::<Vec<_>>();
    format!("reference/cli/commands/{}/", segments.join("/"))
}

fn cli_setting_namespace(path: &str) -> &str {
    path.rsplit_once('.')
        .map(|(namespace, _)| namespace)
        .unwrap_or(path)
}

fn cli_setting_namespace_route(namespace: &str) -> String {
    format!(
        "reference/cli/settings/{}/",
        namespace
            .split('.')
            .map(cli_route_segment)
            .collect::<Vec<_>>()
            .join("/")
    )
}

fn cli_setting_namespaces(reference: &GammaLoopReference) -> BTreeSet<String> {
    let mut namespaces = BTreeSet::new();
    for setting in &reference.settings {
        let segments = setting.path.split('.').collect::<Vec<_>>();
        for end in 1..segments.len() {
            namespaces.insert(segments[..end].join("."));
        }
    }
    namespaces
}

fn visible_cli_commands(reference: &GammaLoopReference) -> Vec<&CliCommand> {
    reference
        .commands
        .iter()
        .filter(|command| !command.hidden && !command.generated_help)
        .collect()
}

fn render_gammaloop_reference_index(product: &str, reference: &GammaLoopReference) -> String {
    let commands = visible_cli_commands(reference);
    let argument_count = commands
        .iter()
        .flat_map(|command| &command.arguments)
        .filter(|argument| !argument.hidden)
        .count();
    let namespaces = cli_setting_namespaces(reference);
    let mut body = format!(
        "<p>This version-specific reference is generated from the compiled Clap parser and serialized settings schemas. Commands use one manpage-style page each; settings follow their native namespace hierarchy. <a href=\"reference/generated/gammaloop-reference.json\">Download the schema-v{} JSON</a>.</p><p class=\"reference-summary\">{} public commands · {argument_count} public arguments · {} settings · {} settings namespaces</p><nav class=\"reference-guide-links\" aria-label=\"Related task guides\"><a href=\"quickstart/cli/\">{}</a><a href=\"guides/process-generation/\">Generate a process</a><a href=\"guides/diagnostics/\">Diagnose a run</a><a href=\"reference/cli/settings/\">Browse settings namespaces</a></nav><div data-reference-filter-root><div class=\"reference-tools\"><label for=\"cli-reference-filter\">Filter commands, aliases, and options</label><input id=\"cli-reference-filter\" type=\"search\" data-reference-filter placeholder=\"Try generate, integrate, or --help\"><output data-reference-filter-count aria-live=\"polite\"></output></div><div data-reference-filter-scope>",
        reference.schema_version,
        commands.len(),
        reference.settings.len(),
        namespaces.len(),
        escape_html(&interface_guide_title(product, "cli")),
    );
    let mut families = BTreeMap::<String, Vec<&CliCommand>>::new();
    for command in &commands {
        let family = command
            .path
            .split_whitespace()
            .nth(1)
            .unwrap_or("gammaloop")
            .to_owned();
        families.entry(family).or_default().push(command);
    }
    for (family, commands) in families {
        body.push_str(&format!(
            "<section class=\"cli-command-family\" id=\"command-family-{}\" data-reference-group><h2><code>{}</code></h2><div class=\"cli-command-list\">",
            escape_html(&slug(&family)),
            escape_html(&family),
        ));
        for command in commands {
            let route = cli_command_route(&command.path);
            let aliases = command
                .aliases
                .iter()
                .filter(|alias| alias.visible)
                .map(|alias| alias.name.as_str())
                .collect::<Vec<_>>()
                .join(" ");
            let options = command
                .arguments
                .iter()
                .filter(|argument| !argument.hidden)
                .flat_map(|argument| {
                    argument
                        .long
                        .iter()
                        .map(|name| format!("--{name}"))
                        .chain(argument.short.iter().map(|name| format!("-{name}")))
                })
                .collect::<Vec<_>>()
                .join(" ");
            body.push_str(&format!(
                "<a class=\"cli-command-link\" id=\"{}\" href=\"{}\" data-reference-entry data-reference-redirect data-reference-search=\"{}\"><span><code>{}</code><small>{}</small></span><span aria-hidden=\"true\">→</span></a>",
                escape_html(&generated_anchor("command", &command.path)),
                escape_html(&route),
                escape_html(&format!("{} {} {} {}", command.path, command.about, aliases, options)),
                escape_html(&command.path),
                if command.about.trim().is_empty() {
                    "Description missing at the parser boundary".to_owned()
                } else {
                    escape_html(&compact_help_summary(&command.about))
                },
            ));
            for argument in command.arguments.iter().filter(|argument| !argument.hidden) {
                body.push_str(&format!(
                    "<a class=\"legacy-reference-anchor\" id=\"{}\" href=\"{}#{}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>",
                    escape_html(&generated_anchor(
                        "argument",
                        &format!("{}::{}", command.path, argument.id),
                    )),
                    escape_html(&route),
                    escape_html(&generated_anchor(
                        "argument",
                        &format!("{}::{}", command.path, argument.id),
                    )),
                ));
            }
        }
        body.push_str("</div></section>");
    }
    body.push_str("</div></div><section class=\"cli-settings-callout\"><p class=\"portal-kicker\">Configuration reference</p><h2><a href=\"reference/cli/settings/\">Settings namespaces</a></h2><p>Browse parent and child namespaces instead of searching one monolithic settings list.</p></section>");
    let mut legacy_setting_groups = BTreeMap::new();
    for setting in &reference.settings {
        let namespace = cli_setting_namespace(&setting.path);
        let route = cli_setting_namespace_route(namespace);
        body.push_str(&format!(
            "<a class=\"legacy-reference-anchor\" id=\"{}\" href=\"{}#{}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>",
            escape_html(&generated_anchor("setting", &setting.path)),
            escape_html(&route),
            escape_html(&generated_anchor("setting", &setting.path)),
        ));
        let segments = setting.path.split('.').take(2).collect::<Vec<_>>();
        legacy_setting_groups
            .entry(segments.join("."))
            .or_insert_with(|| namespace.to_owned());
    }
    for (group, fallback) in legacy_setting_groups {
        let target = if namespaces.contains(&group) {
            &group
        } else {
            &fallback
        };
        body.push_str(&format!(
            "<a class=\"legacy-reference-anchor\" id=\"settings-group-{}\" href=\"{}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>",
            escape_html(&slug(&group)),
            escape_html(&cli_setting_namespace_route(target)),
        ));
    }
    reference_page(product, "CLI commands and settings", &body)
}

fn render_cli_argument_names(argument: &CliArgument) -> String {
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
    names.extend(argument.aliases.iter().filter(|alias| alias.visible).map(
        |alias| match alias.kind {
            CliAliasKind::Name => alias.name.clone(),
            CliAliasKind::ShortFlag => format!("-{}", alias.name),
            CliAliasKind::LongFlag => format!("--{}", alias.name),
        },
    ));
    names.join(", ")
}

fn render_cli_argument_references(command: &CliCommand, ids: &[String]) -> String {
    ids.iter()
        .map(|id| {
            let Some(argument) = command.arguments.iter().find(|argument| argument.id == *id)
            else {
                return format!("<code>{}</code>", escape_html(id));
            };
            let label = format!(
                "<code>{}</code>",
                escape_html(&render_cli_argument_names(argument))
            );
            if argument.hidden {
                label
            } else {
                let anchor =
                    generated_anchor("argument", &format!("{}::{}", command.path, argument.id));
                format!("<a href=\"#{}\">{label}</a>", escape_html(&anchor))
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}

fn render_cli_argument_index(command: &CliCommand) -> String {
    let groups = [
        (
            "Positional arguments",
            command
                .arguments
                .iter()
                .filter(|argument| !argument.hidden && argument.positional)
                .collect::<Vec<_>>(),
        ),
        (
            "Options",
            command
                .arguments
                .iter()
                .filter(|argument| !argument.hidden && !argument.positional && !argument.inherited)
                .collect::<Vec<_>>(),
        ),
        (
            "Inherited and global options",
            command
                .arguments
                .iter()
                .filter(|argument| !argument.hidden && argument.inherited)
                .collect::<Vec<_>>(),
        ),
    ];
    if groups.iter().all(|(_, arguments)| arguments.is_empty()) {
        return String::new();
    }
    let mut body = "<nav class=\"reference-jump-nav\" aria-labelledby=\"argument-overview-title\"><h2 id=\"argument-overview-title\">Arguments and options</h2><div class=\"reference-jump-groups\">".to_owned();
    for (title, arguments) in groups {
        if arguments.is_empty() {
            continue;
        }
        body.push_str(&format!(
            "<section class=\"reference-jump-group\"><h3>{title}</h3><ul class=\"reference-jump-list\">"
        ));
        for argument in arguments {
            let anchor =
                generated_anchor("argument", &format!("{}::{}", command.path, argument.id));
            let summary = compact_help_summary(&argument.help);
            body.push_str(&format!(
                "<li><a href=\"#{}\"><span><code>{}</code>{}</span><small>{}</small></a></li>",
                escape_html(&anchor),
                escape_html(&render_cli_argument_names(argument)),
                if summary.is_empty() {
                    String::new()
                } else {
                    format!(
                        "<span class=\"reference-jump-summary\">{}</span>",
                        escape_html(&summary),
                    )
                },
                if argument.required {
                    "required"
                } else if argument.inherited {
                    "inherited"
                } else {
                    "optional"
                },
            ));
        }
        body.push_str("</ul></section>");
    }
    body.push_str("</div></nav>");
    body
}

fn render_cli_arguments<'a>(
    command: &CliCommand,
    title: &str,
    arguments: impl Iterator<Item = &'a CliArgument>,
) -> String {
    let arguments = arguments.collect::<Vec<_>>();
    if arguments.is_empty() {
        return String::new();
    }
    let mut body = format!("<section class=\"cli-argument-section\"><h2>{title}</h2>");
    for argument in arguments {
        let anchor = generated_anchor("argument", &format!("{}::{}", command.path, argument.id));
        let action = match argument.action {
            CliArgumentAction::Set => "sets one value",
            CliArgumentAction::Append => "accepts repeated values",
            CliArgumentAction::SetTrue => "sets true",
            CliArgumentAction::SetFalse => "sets false",
            CliArgumentAction::Count => "counts occurrences",
            CliArgumentAction::Help
            | CliArgumentAction::HelpShort
            | CliArgumentAction::HelpLong => "shows help",
            CliArgumentAction::Version => "shows the version",
        };
        let arity = match (argument.arity.min, argument.arity.max) {
            (0, Some(0)) => "no value".to_owned(),
            (min, Some(max)) if min == max => format!("{min} value(s)"),
            (min, Some(max)) => format!("{min}–{max} values"),
            (min, None) => format!("{min} or more values"),
        };
        let names = render_cli_argument_names(argument);
        let mut facts = format!(
            "<div><dt>Behavior</dt><dd>{action}; {arity}{}</dd></div>",
            if argument.required { "; required" } else { "" },
        );
        if argument.global || argument.inherited {
            facts.push_str(&format!(
                "<div><dt>Scope</dt><dd>{}</dd></div>",
                if argument.inherited {
                    "inherited global"
                } else {
                    "global"
                },
            ));
        }
        if argument.exclusive {
            facts.push_str("<div><dt>Exclusive</dt><dd>yes</dd></div>");
        }
        if argument.require_equals {
            facts.push_str("<div><dt>Requires equals sign</dt><dd>yes</dd></div>");
        }
        if !argument.value_names.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Value names</dt><dd>{}</dd></div>",
                argument
                    .value_names
                    .iter()
                    .map(|name| format!("<code>{}</code>", escape_html(name)))
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }
        if let Some(delimiter) = argument.value_delimiter {
            facts.push_str(&format!(
                "<div><dt>Value delimiter</dt><dd><code>{}</code></dd></div>",
                escape_html(&delimiter.to_string()),
            ));
        }
        if let Some(terminator) = &argument.value_terminator {
            facts.push_str(&format!(
                "<div><dt>Value terminator</dt><dd><code>{}</code></dd></div>",
                escape_html(terminator),
            ));
        }
        if !argument.defaults.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Default</dt><dd><code>{}</code></dd></div>",
                escape_html(&argument.defaults.join(", ")),
            ));
        }
        if !argument.default_missing_values.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Default when value omitted</dt><dd>{}</dd></div>",
                argument
                    .default_missing_values
                    .iter()
                    .map(|value| format!("<code>{}</code>", escape_html(value)))
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }
        if !argument.possible_values.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Values</dt><dd>{}</dd></div>",
                argument
                    .possible_values
                    .iter()
                    .map(|value| format!("<code>{}</code>", escape_html(value)))
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }
        if !argument.requires.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Requires</dt><dd>{}</dd></div>",
                render_cli_argument_references(command, &argument.requires),
            ));
        }
        if !argument.conflicts_with.is_empty() {
            facts.push_str(&format!(
                "<div><dt>Conflicts with</dt><dd>{}</dd></div>",
                render_cli_argument_references(command, &argument.conflicts_with),
            ));
        }
        body.push_str(&format!(
            "<article class=\"cli-argument\" id=\"{}\"><header><h3><code>{}</code></h3><a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to {}\">#</a></header>{}<dl class=\"reference-facts\">{facts}</dl></article>",
            escape_html(&anchor),
            escape_html(&names),
            escape_html(&anchor),
            escape_html(&names),
            if argument.help.trim().is_empty() {
                "<p class=\"reference-missing\">No description is available for this argument.</p>".to_owned()
            } else {
                render_cli_help(&argument.help)
            },
        ));
    }
    body.push_str("</section>");
    body
}

fn render_cli_command_page(
    product: &str,
    reference: &GammaLoopReference,
    command: &CliCommand,
) -> String {
    let command_anchor = generated_anchor("command", &command.path);
    let command_line_guide = interface_guide_title(product, "cli");
    let (workflow_route, workflow_title) =
        if command.path == "gammaLoop inspect" || command.path == "gammaLoop approach" {
            (
                "reference/interfaces/#sample-evaluation-contract",
                "Sample evaluation contract",
            )
        } else if command.path.starts_with("gammaLoop generate")
            || command.path == "gammaLoop import graphs"
            || command.path.starts_with("gammaLoop save dot")
            || command.path.starts_with("gammaLoop save uv-forest")
        {
            ("guides/process-generation/", "Process generation workflow")
        } else if command.path.starts_with("gammaLoop integrate")
            || command.path.starts_with("gammaLoop evaluate")
            || command.path.starts_with("gammaLoop inspect")
            || command.path.starts_with("gammaLoop display integrand")
            || command.path.starts_with("gammaLoop display quantities")
            || command.path.starts_with("gammaLoop display observables")
            || command.path.starts_with("gammaLoop display selectors")
        {
            (
                "guides/events-and-observables/",
                "Events and observables workflow",
            )
        } else if command.path.starts_with("gammaLoop set")
            || command.path.starts_with("gammaLoop display settings")
            || command.path.starts_with("gammaLoop approach")
            || command.path.starts_with("gammaLoop bench")
            || command.path.starts_with("gammaLoop profile")
            || command.path.starts_with("gammaLoop batch")
        {
            ("guides/diagnostics/", "Configuration and diagnostics guide")
        } else {
            ("quickstart/cli/", command_line_guide.as_str())
        };
    let related_workflow = match command.path.as_str() {
        "gammaLoop inspect" => {
            "<a href=\"guides/events-and-observables/\">Events and observables workflow</a>"
        }
        "gammaLoop approach" => {
            "<a href=\"guides/diagnostics/\">Configuration and diagnostics guide</a>"
        }
        _ => "",
    };
    let mut body = format!(
        "<span id=\"{}\" aria-hidden=\"true\"></span><p class=\"reference-context\"><a href=\"reference/cli/\">CLI reference</a> <span aria-hidden=\"true\">/</span> Command</p><nav class=\"reference-guide-links\" aria-label=\"CLI reference navigation\"><a href=\"reference/cli/\">All commands</a><a href=\"reference/cli/settings/\">Settings namespaces</a><a href=\"{}\">{}</a>{}</nav>",
        escape_html(&command_anchor),
        escape_html(workflow_route),
        escape_html(workflow_title),
        related_workflow,
    );
    if command.about.trim().is_empty() {
        body.push_str(
            "<p class=\"reference-missing\">No description is available for this command.</p>",
        );
    } else {
        body.push_str(&format!(
            "<div class=\"reference-description\">{}</div>",
            render_cli_help(&command.about),
        ));
    }
    body.push_str(&format!(
        "<h2>Synopsis</h2><div class=\"reference-copy-row\"><pre class=\"api-signature\" id=\"{}-usage\"><code data-lang=\"shell\">{}</code></pre><button type=\"button\" data-copy-target=\"{}-usage\">Copy usage</button></div>",
        escape_html(&command_anchor),
        escape_html(&command.usage),
        escape_html(&command_anchor),
    ));
    let aliases = command
        .aliases
        .iter()
        .filter(|alias| alias.visible)
        .map(|alias| match alias.kind {
            CliAliasKind::Name => alias.name.clone(),
            CliAliasKind::ShortFlag => format!("-{}", alias.name),
            CliAliasKind::LongFlag => format!("--{}", alias.name),
        })
        .collect::<Vec<_>>();
    if !aliases.is_empty() {
        body.push_str(&format!(
            "<p class=\"reference-aliases\"><strong>Aliases:</strong> {}</p>",
            aliases
                .iter()
                .map(|alias| format!("<code>{}</code>", escape_html(alias)))
                .collect::<Vec<_>>()
                .join(", "),
        ));
    }
    let prefix = format!("{} ", command.path);
    let direct_children = visible_cli_commands(reference)
        .into_iter()
        .filter(|candidate| {
            candidate.path.starts_with(&prefix) && !candidate.path[prefix.len()..].contains(' ')
        })
        .collect::<Vec<_>>();
    let parent = command.path.rsplit_once(' ').and_then(|(parent, _)| {
        reference
            .commands
            .iter()
            .find(|candidate| candidate.path == parent && !candidate.hidden)
    });
    if parent.is_some() || !direct_children.is_empty() {
        body.push_str("<section class=\"cli-related\"><h2>Command hierarchy</h2><nav aria-label=\"Related commands\">");
        if let Some(parent) = parent {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Parent</span><code>{}</code></a>",
                escape_html(&cli_command_route(&parent.path)),
                escape_html(&parent.path),
            ));
        }
        for child in direct_children {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Subcommand</span><code>{}</code></a>",
                escape_html(&cli_command_route(&child.path)),
                escape_html(&child.path),
            ));
        }
        body.push_str("</nav></section>");
    }
    body.push_str(&render_cli_argument_index(command));
    body.push_str(&render_cli_arguments(
        command,
        "Positional arguments",
        command
            .arguments
            .iter()
            .filter(|argument| !argument.hidden && argument.positional),
    ));
    body.push_str(&render_cli_arguments(
        command,
        "Options",
        command
            .arguments
            .iter()
            .filter(|argument| !argument.hidden && !argument.positional && !argument.inherited),
    ));
    body.push_str(&render_cli_arguments(
        command,
        "Inherited and global options",
        command
            .arguments
            .iter()
            .filter(|argument| !argument.hidden && argument.inherited),
    ));
    reference_page(product, &command.path, &body)
}

fn render_cli_settings_index(product: &str, reference: &GammaLoopReference) -> String {
    let namespaces = cli_setting_namespaces(reference);
    let mut body = format!(
        "<p>The {} generated settings are organized by their actual schema namespaces. Open a namespace to see its direct settings and child namespaces.</p><nav class=\"reference-guide-links\" aria-label=\"CLI settings navigation\"><a href=\"reference/cli/\">All commands</a><a href=\"guides/diagnostics/\">Precedence and live reload</a></nav><div data-reference-filter-root><div class=\"reference-tools\"><label for=\"settings-namespace-filter\">Filter settings namespaces</label><input id=\"settings-namespace-filter\" type=\"search\" data-reference-filter placeholder=\"Try runtime.integrator or cli.global\"><output data-reference-filter-count aria-live=\"polite\"></output></div><div class=\"settings-namespace-list\" data-reference-filter-scope>",
        reference.settings.len(),
    );
    for namespace in namespaces {
        let direct = reference
            .settings
            .iter()
            .filter(|setting| cli_setting_namespace(&setting.path) == namespace)
            .count();
        let prefix = format!("{namespace}.");
        let descendants = reference
            .settings
            .iter()
            .filter(|setting| setting.path.starts_with(&prefix))
            .count();
        body.push_str(&format!(
            "<a href=\"{}\" data-reference-entry data-reference-search=\"{}\"><code>{}</code><span>{direct} direct · {descendants} total settings</span></a>",
            escape_html(&cli_setting_namespace_route(&namespace)),
            escape_html(&namespace),
            escape_html(&namespace),
        ));
    }
    body.push_str("</div></div>");
    reference_page(product, "Settings namespaces", &body)
}

fn render_cli_setting(setting: &SettingReference, body: &mut String) {
    let default = setting.default.as_ref().map_or_else(
        || "<span>Not set</span>".to_owned(),
        |value| match value {
            Value::Array(_) | Value::Object(_) => format!(
                "<details><summary>Show structured value</summary><pre><code>{}</code></pre></details>",
                escape_html(&serde_json::to_string_pretty(value).unwrap_or_else(|_| value.to_string())),
            ),
            _ => format!("<code>{}</code>", escape_html(&value.to_string())),
        },
    );
    let anchor = generated_anchor("setting", &setting.path);
    let name = setting
        .path
        .rsplit('.')
        .next()
        .unwrap_or(setting.path.as_str());
    body.push_str(&format!(
        "<section class=\"reference-setting reference-setting-flat\" id=\"{}\"><header><h3><code>{}</code></h3><a class=\"reference-permalink\" href=\"#{}\" aria-label=\"Permanent link to {}\">#</a></header><p class=\"reference-setting-path\">{}</p><div class=\"reference-description\">{}</div><dl class=\"reference-facts\"><div><dt>Type</dt><dd><code>{}</code></dd></div><div><dt>Required</dt><dd>{}</dd></div><div class=\"reference-default\"><dt>Default</dt><dd>{}</dd></div><div><dt>Allowed values</dt><dd>{}</dd></div></dl></section>",
        escape_html(&anchor),
        escape_html(name),
        escape_html(&anchor),
        escape_html(&setting.path),
        escape_html(&setting.path),
        if setting.description.trim().is_empty() {
            "<p class=\"reference-missing\">No description is available for this setting.</p>".to_owned()
        } else {
            render_cli_help(&setting.description)
        },
        escape_html(&setting.value_type),
        if setting.required { "yes" } else { "no" },
        default,
        if setting.possible_values.is_empty() {
            "Any value accepted by the declared type".to_owned()
        } else {
            setting
                .possible_values
                .iter()
                .map(|value| format!("<code>{}</code>", escape_html(value)))
                .collect::<Vec<_>>()
                .join(", ")
        },
    ));
}

fn render_cli_settings_namespace_page(
    product: &str,
    reference: &GammaLoopReference,
    namespace: &str,
    namespaces: &BTreeSet<String>,
) -> String {
    let parent = namespace.rsplit_once('.').map(|(parent, _)| parent);
    let prefix = format!("{namespace}.");
    let children = namespaces
        .iter()
        .filter(|candidate| {
            candidate.starts_with(&prefix) && !candidate[prefix.len()..].contains('.')
        })
        .collect::<Vec<_>>();
    let settings = reference
        .settings
        .iter()
        .filter(|setting| cli_setting_namespace(&setting.path) == namespace)
        .collect::<Vec<_>>();
    let mut body = "<p class=\"reference-context\"><a href=\"reference/cli/settings/\">Settings namespaces</a> <span aria-hidden=\"true\">/</span> Namespace</p><nav class=\"reference-guide-links\" aria-label=\"Settings reference navigation\"><a href=\"reference/cli/\">All commands</a><a href=\"reference/cli/settings/\">All namespaces</a><a href=\"guides/diagnostics/\">Configuration and diagnostics guide</a></nav>".to_owned();
    if parent.is_some() || !children.is_empty() {
        body.push_str("<section class=\"settings-hierarchy\"><h2>Namespace hierarchy</h2><nav aria-label=\"Related settings namespaces\">");
        if let Some(parent) = parent {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Parent</span><code>{}</code></a>",
                escape_html(&cli_setting_namespace_route(parent)),
                escape_html(parent),
            ));
        }
        for child in children {
            body.push_str(&format!(
                "<a href=\"{}\"><span>Child</span><code>{}</code></a>",
                escape_html(&cli_setting_namespace_route(child)),
                escape_html(child),
            ));
        }
        body.push_str("</nav></section>");
    }
    if settings.is_empty() {
        body.push_str("<p>This namespace owns child namespaces but no settings directly.</p>");
    } else {
        body.push_str("<nav class=\"reference-jump-nav\" aria-labelledby=\"setting-overview-title\"><h2 id=\"setting-overview-title\">Settings in this namespace</h2><div class=\"reference-jump-groups\"><section class=\"reference-jump-group\"><h3>Setting</h3><ul class=\"reference-jump-list\">");
        for setting in &settings {
            let anchor = generated_anchor("setting", &setting.path);
            let name = setting
                .path
                .rsplit('.')
                .next()
                .unwrap_or(setting.path.as_str());
            let summary = compact_help_summary(&setting.description);
            body.push_str(&format!(
                "<li><a href=\"#{}\"><span><code>{}</code>{}</span><small>{}</small></a></li>",
                escape_html(&anchor),
                escape_html(name),
                if summary.is_empty() {
                    String::new()
                } else {
                    format!(
                        "<span class=\"reference-jump-summary\">{}</span>",
                        escape_html(&summary),
                    )
                },
                escape_html(&setting.value_type),
            ));
        }
        body.push_str(
            "</ul></section></div></nav><h2>Settings</h2><div class=\"reference-setting-list\">",
        );
        for setting in settings {
            render_cli_setting(setting, &mut body);
        }
        body.push_str("</div>");
    }
    reference_page(product, &format!("{namespace} settings"), &body)
}

fn write_gammaloop_reference_pages(
    product: &str,
    reference: &GammaLoopReference,
    site: &Path,
) -> Result<Vec<SitePage>> {
    let mut pages = Vec::new();
    let mut routes = BTreeSet::new();
    for command in visible_cli_commands(reference) {
        let route = cli_command_route(&command.path);
        ensure!(
            routes.insert(route.clone()),
            "duplicate generated CLI route {route}"
        );
        let destination = site.join(&route);
        fs::create_dir_all(&destination)?;
        fs::write(
            destination.join("index.html"),
            render_cli_command_page(product, reference, command),
        )?;
        pages.push(SitePage::new(route, command.path.clone(), "CLI reference"));
    }
    let settings_index = site.join("reference/cli/settings");
    fs::create_dir_all(&settings_index)?;
    fs::write(
        settings_index.join("index.html"),
        render_cli_settings_index(product, reference),
    )?;
    pages.push(SitePage::new(
        "reference/cli/settings/",
        "Settings namespaces",
        "CLI reference",
    ));
    let namespaces = cli_setting_namespaces(reference);
    for namespace in &namespaces {
        let route = cli_setting_namespace_route(namespace);
        ensure!(
            routes.insert(route.clone()),
            "duplicate generated settings route {route}"
        );
        let destination = site.join(&route);
        fs::create_dir_all(&destination)?;
        fs::write(
            destination.join("index.html"),
            render_cli_settings_namespace_page(product, reference, namespace, &namespaces),
        )?;
        pages.push(SitePage::new(
            route,
            format!("{namespace} settings"),
            "CLI reference",
        ));
    }
    Ok(pages)
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
    let source_metadata = fs::symlink_metadata(source)?;
    ensure!(
        source_metadata.is_dir() && !source_metadata.file_type().is_symlink(),
        "documentation source is not a directory: {}",
        source.display()
    );
    match fs::symlink_metadata(destination) {
        Ok(metadata) => {
            ensure!(
                metadata.is_dir() && !metadata.file_type().is_symlink(),
                "documentation destination is not a directory: {}",
                destination.display()
            );
            for entry in WalkDir::new(destination) {
                let entry = entry?;
                ensure!(
                    entry.file_type().is_dir() || entry.file_type().is_file(),
                    "unsupported documentation artifact {}",
                    entry.path().display()
                );
            }
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(error.into()),
    }
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
    let resolved_source = fs::canonicalize(source)
        .wrap_err_with(|| format!("failed to resolve {}", source.display()))?;
    let existing_destination = closest_existing_directory(destination)
        .context("documentation destination has no existing directory ancestor")?;
    let resolved_destination = fs::canonicalize(&existing_destination)
        .wrap_err_with(|| format!("failed to resolve {}", existing_destination.display()))?
        .join(destination.strip_prefix(&existing_destination)?);
    ensure!(
        resolved_source != resolved_destination
            && !resolved_source.starts_with(&resolved_destination)
            && !resolved_destination.starts_with(&resolved_source),
        "cannot replace overlapping documentation trees {} and {}",
        source.display(),
        destination.display()
    );
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
        } else if !entry.file_type().is_dir() {
            bail!(
                "unsupported documentation artifact {}",
                entry.path().display()
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
    link: LinkReference<'_>,
    failures: &mut Vec<String>,
    patterns: &LinkValidationPatterns,
    state: LinkValidationState<'_>,
) -> Result<()> {
    let (linked_pages, link_rewrites, local_paths, snapshot_tag) = state;
    let Some(mut resolved) = resolve_local_link(output, link.source, link.href, local_paths) else {
        return Ok(());
    };
    let source_display = link
        .source
        .strip_prefix(output)
        .unwrap_or(link.source)
        .display();
    if !resolved.starts_with(output) {
        failures.push(format!(
            "{source_display} -> {} (escapes generated output root)",
            link.href
        ));
        return Ok(());
    }
    let mut target_is_file = local_paths
        .entry(resolved.clone())
        .or_insert_with(|| {
            fs::metadata(&resolved)
                .ok()
                .map(|metadata| metadata.file_type())
        })
        .as_ref()
        .is_some_and(|file_type| file_type.is_file());
    if !target_is_file
        && let Some(tag) = snapshot_tag
        && link
            .source
            .strip_prefix(output)
            .is_ok_and(|source| source.starts_with("developers"))
        && let Ok(relative) = resolved.strip_prefix(output)
    {
        let mut components = relative.components();
        if matches!(components.next(), Some(Component::Normal(part)) if part == "products")
            && let Some(Component::Normal(product)) = components.next()
            && matches!(components.next(), Some(Component::Normal(part)) if part == "latest")
        {
            // Snapshot bundles are merged beside the deployed `latest` tree. Validate
            // developer links to that tree against the same route in this snapshot.
            let mut candidate = output
                .join("products")
                .join(product)
                .join("snapshots")
                .join(tag);
            candidate.extend(components.map(|component| component.as_os_str()));
            let candidate_is_directory = local_paths
                .entry(candidate.clone())
                .or_insert_with(|| {
                    fs::metadata(&candidate)
                        .ok()
                        .map(|metadata| metadata.file_type())
                })
                .as_ref()
                .is_some_and(|file_type| file_type.is_dir());
            if candidate_is_directory {
                candidate.push("index.html");
            }
            target_is_file = local_paths
                .entry(candidate.clone())
                .or_insert_with(|| {
                    fs::metadata(&candidate)
                        .ok()
                        .map(|metadata| metadata.file_type())
                })
                .as_ref()
                .is_some_and(|file_type| file_type.is_file());
            if target_is_file {
                resolved = candidate;
            }
        }
    }
    if !target_is_file {
        let target_path = link.href.split(['?', '#']).next().unwrap_or_default();
        let inherited_rust_path = target_path.contains("::")
            || (!target_path.is_empty()
                && !target_path.contains(['/', '.'])
                && !target_path.starts_with('#'));
        if link.rustdoc
            && inherited_rust_path
            && let Some(attribute) = &link.attribute
        {
            // Rustdoc can inline dependency prose whose Rust paths were never
            // resolved to site URLs. Retain the prose, but remove only this
            // confirmed-missing, non-navigable attribute from the copied site.
            link_rewrites
                .entry(link.source.to_path_buf())
                .or_default()
                .push(LinkRewrite {
                    range: attribute.attribute_range.clone(),
                    expected: attribute.attribute.clone(),
                    replacement: String::new(),
                    target_before: link.href.to_owned(),
                    target_after: None,
                });
            return Ok(());
        }
        failures.push(format!("{source_display} -> {}", link.href));
        return Ok(());
    }
    let Some((_, fragment)) = link.href.split_once('#') else {
        return Ok(());
    };
    let fragment = fragment.split('?').next().unwrap_or_default();
    if fragment.is_empty() {
        return Ok(());
    }
    let load_page = |path: &Path| -> Result<(HashSet<String>, Option<String>)> {
        let html = fs::read_to_string(path)?;
        Ok(patterns.page_metadata(&patterns.rendered_html(&html)))
    };
    let mut fragment_exists = |anchors: &HashSet<String>| {
        let rustdoc_source = resolved
            .ancestors()
            .any(|ancestor| ancestor.ends_with(Path::new("reference/rust/src")));
        let rustdoc_target = resolved
            .ancestors()
            .any(|ancestor| ancestor.ends_with(Path::new("reference/rust")));
        let exact_fragment = anchors.contains(fragment)
            || (rustdoc_source
                && fragment.split_once('-').is_some_and(|(start, end)| {
                    start.chars().all(|character| character.is_ascii_digit())
                        && end.chars().all(|character| character.is_ascii_digit())
                        && anchors.contains(start)
                        && anchors.contains(end)
                }));
        let normalized = if exact_fragment {
            Some(fragment.to_owned())
        } else if rustdoc_target
            && let Some(method) = fragment.strip_prefix("tymethod.")
            && anchors.contains(&format!("method.{method}"))
        {
            Some(format!("method.{method}"))
        } else if rustdoc_target
            && let Some((variant, field)) = fragment.split_once(".field.")
            && variant.starts_with("variant.")
            && !field.is_empty()
            && field.chars().all(|character| character.is_ascii_digit())
            && anchors.contains(variant)
        {
            Some(variant.to_owned())
        } else {
            None
        };
        let Some(normalized) = normalized else {
            return false;
        };
        if normalized != fragment {
            let Some(attribute) = &link.attribute else {
                return false;
            };
            let (prefix, raw_fragment) = link
                .href
                .split_once('#')
                .expect("fragment-bearing links were checked above");
            let query = raw_fragment
                .find('?')
                .map_or("", |index| &raw_fragment[index..]);
            link_rewrites
                .entry(link.source.to_path_buf())
                .or_default()
                .push(LinkRewrite {
                    range: attribute.target_range.clone(),
                    expected: link.href.to_owned(),
                    replacement: format!("{prefix}#{normalized}{query}"),
                    target_before: link.href.to_owned(),
                    target_after: Some(format!("{prefix}#{normalized}{query}")),
                });
        }
        true
    };
    if !linked_pages.contains_key(&resolved) {
        linked_pages.insert(resolved.clone(), load_page(&resolved)?);
    }
    let redirect_target = {
        let (anchors, redirect_target) = linked_pages
            .get(&resolved)
            .expect("linked page was inserted before fragment validation");
        if fragment_exists(anchors) {
            return Ok(());
        }
        redirect_target.clone()
    };
    if let Some(redirect_target) = redirect_target
        && let Some(redirected) =
            resolve_local_link(output, &resolved, &redirect_target, local_paths)
        && redirected.starts_with(output)
    {
        let redirected_is_file = local_paths
            .entry(redirected.clone())
            .or_insert_with(|| {
                fs::metadata(&redirected)
                    .ok()
                    .map(|metadata| metadata.file_type())
            })
            .as_ref()
            .is_some_and(|file_type| file_type.is_file());
        if !redirected_is_file {
            failures.push(format!(
                "{source_display} -> {} (missing redirect target {})",
                link.href,
                redirected.display()
            ));
            return Ok(());
        }
        if !linked_pages.contains_key(&redirected) {
            linked_pages.insert(redirected.clone(), load_page(&redirected)?);
        }
        if fragment_exists(
            &linked_pages
                .get(&redirected)
                .expect("redirected page was inserted before fragment validation")
                .0,
        ) {
            return Ok(());
        }
    }
    failures.push(format!(
        "{source_display} -> {} (missing fragment)",
        link.href
    ));
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

fn resolve_local_link(
    output: &Path,
    source: &Path,
    href: &str,
    local_paths: &mut LocalPathIndex,
) -> Option<PathBuf> {
    let external_scheme = href.split_once(':').is_some_and(|(scheme, remainder)| {
        !remainder.starts_with(':')
            && scheme
                .chars()
                .next()
                .is_some_and(|character| character.is_ascii_alphabetic())
            && scheme
                .chars()
                .all(|character| character.is_ascii_alphanumeric() || "+-.".contains(character))
    });
    if href.is_empty() || href.starts_with("//") || external_scheme {
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
    let target_is_directory = local_paths
        .entry(normalized.clone())
        .or_insert_with(|| {
            fs::metadata(&normalized)
                .ok()
                .map(|metadata| metadata.file_type())
        })
        .as_ref()
        .is_some_and(|file_type| file_type.is_dir());
    if href.ends_with('/') || target_is_directory {
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
        let sidebar =
            builder.site_sidebar(product, &metadata, &current, "", BuildScope::FullSite, &[]);
        assert!(!sidebar.contains("developers/"));
        assert!(!sidebar.contains("For developers"));
        assert!(sidebar.contains("href=\"../../../../citations/#gammaloop\">Cite GammaLoop</a>"));
    }

    #[test]
    fn python_sidebar_lists_only_the_current_module_and_marks_one_page() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "gammaloop")
            .unwrap();
        let component = product.python_components.first().unwrap();
        let module_route = format!("reference/python/{}/", component.id);
        let engine_route = format!("{module_route}Engine/");
        let result_route = format!("{module_route}results/EvaluationResult/");
        let generated = vec![
            SitePage::new(&engine_route, "gammaloop.Engine", "Python API"),
            SitePage::new(&result_route, "gammaloop.EvaluationResult", "Python API"),
            SitePage::new(
                "reference/cli/generate/",
                "gammaloop generate",
                "CLI reference",
            ),
        ];
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Latest,
            snapshot_tag: None,
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/gammaloop/latest/".to_owned(),
            components: Vec::new(),
        };
        let current = SitePage::new(&result_route, "EvaluationResult", "Python API");
        let sidebar = builder.site_sidebar(
            product,
            &metadata,
            &current,
            "",
            BuildScope::FullSite,
            &generated,
        );

        assert!(sidebar.contains(&format!(
            "<p class=\"sidebar-group-title\">{} symbols</p>",
            python_display_module(component)
        )));
        assert!(sidebar.contains(&format!(
            "<nav class=\"sidebar-group sidebar-python-group\" aria-label=\"{} Python symbols\">",
            python_display_module(component)
        )));
        assert!(sidebar.contains(&format!("href=\"{module_route}\">Module overview</a>")));
        assert!(sidebar.contains(&format!("href=\"{engine_route}\"><code>Engine</code></a>")));
        assert!(sidebar.contains(
            "href=\"reference/python/gammaloop-python/results/EvaluationResult/\" aria-current=\"page\"><code>results.EvaluationResult</code></a>"
        ));
        assert!(
            sidebar.contains("href=\"reference/python/\" aria-current=\"location\">Python API</a>")
        );
        assert!(!sidebar.contains("gammaloop generate"));
        assert_eq!(sidebar.matches("aria-current=\"page\"").count(), 1);

        let module = SitePage::new(&module_route, "gammaloop Python API", "Reference");
        let module_sidebar = builder.site_sidebar(
            product,
            &metadata,
            &module,
            "",
            BuildScope::FullSite,
            &generated,
        );
        assert!(module_sidebar.contains(&format!(
            "href=\"{module_route}\" aria-current=\"page\">Module overview</a>"
        )));
        assert_eq!(module_sidebar.matches("aria-current=\"page\"").count(), 1);

        let tutorial = SitePage::new("tutorial/", "Tutorial", "Tutorial");
        let tutorial_sidebar = builder.site_sidebar(
            product,
            &metadata,
            &tutorial,
            "",
            BuildScope::FullSite,
            &generated,
        );
        assert!(!tutorial_sidebar.contains("sidebar-python-group"));
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
    fn cli_special_path_segments_have_stable_routes() {
        assert_eq!(
            cli_command_route("gammaLoop !"),
            "reference/cli/commands/gammaloop/bang/"
        );
        assert_eq!(
            cli_command_route("gammaLoop result[]"),
            "reference/cli/commands/gammaloop/result-items/"
        );
    }

    #[test]
    fn product_bundles_render_authored_typst_but_allow_package_assets() {
        let builder = SiteBuilder::discover().unwrap();
        let site = tempfile::tempdir().unwrap();
        let package = site
            .path()
            .join("products/linnet/latest/assets/typst/linnest/src");
        fs::create_dir_all(&package).unwrap();
        fs::write(package.join("lib.typ"), "#let package-asset = true\n").unwrap();
        builder
            .validate_generated_links(site.path(), false, None)
            .unwrap();

        let leaked = site
            .path()
            .join("products/linnet/latest/reference/content/changelog.typ");
        fs::create_dir_all(leaked.parent().unwrap()).unwrap();
        fs::write(&leaked, "= Changelog\n").unwrap();
        let error = builder
            .validate_generated_links(site.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains(
            "authored Typst source was published instead of rendered: products/linnet/latest/reference/content/changelog.typ"
        ));
        assert!(package.join("lib.typ").is_file());
    }

    #[test]
    fn cli_reference_index_and_manpages_preserve_parser_semantics() {
        let builder = SiteBuilder::discover().unwrap();
        let (reference, _) = builder.generated_references().unwrap();
        let commands = visible_cli_commands(&reference);
        let index = render_gammaloop_reference_index("GammaLoop", &reference);
        let (command, argument) = commands
            .iter()
            .find_map(|command| {
                command
                    .arguments
                    .iter()
                    .find(|argument| argument.long.as_deref() == Some("clean-state"))
                    .map(|argument| (*command, argument))
            })
            .unwrap();
        let route = cli_command_route(&command.path);
        let argument_anchor =
            generated_anchor("argument", &format!("{}::{}", command.path, argument.id));
        let command_page = render_cli_command_page("GammaLoop", &reference, command);
        let quit_page = render_cli_command_page(
            "GammaLoop",
            &reference,
            commands
                .iter()
                .find(|command| command.path == "gammaLoop quit")
                .expect("the CLI exposes quit"),
        );
        let all_command_pages = commands
            .iter()
            .map(|command| render_cli_command_page("GammaLoop", &reference, command))
            .collect::<String>();

        assert!(index.contains(&format!("{} public commands", commands.len())));
        assert!(index.contains("class=\"cli-command-link\""));
        assert!(index.contains("data-reference-filter"));
        assert!(
            index.contains("data-reference-entry data-reference-redirect data-reference-search")
        );
        assert!(!index.contains("gammaLoop help"));
        assert!(!index.contains("<details"));
        assert!(
            index.contains("href=\"quickstart/cli/\">Using GammaLoop from the command line</a>")
        );
        assert!(index.contains(&format!(
            "id=\"{}\" href=\"{route}#{argument_anchor}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"",
            argument_anchor,
        )));
        assert!(command_page.contains(&format!(
            "<article class=\"cli-argument\" id=\"{argument_anchor}\""
        )));
        assert!(command_page.contains("--clean-state"));
        assert!(command_page.contains("sets true; no value"));
        assert!(command_page.contains("Copy usage"));
        assert!(!command_page.contains("Copy command"));
        assert!(!command_page.contains("--clean-state &lt;"));
        assert!(!command_page.contains("<details"));
        assert!(!command_page.contains("data-reference-entry"));
        assert!(!command_page.contains("data-reference-search"));
        assert!(
            quit_page
                .contains("href=\"quickstart/cli/\">Using GammaLoop from the command line</a>")
        );
        assert!(all_command_pages.contains(
            "Process reference: <code>#&lt;id&gt;</code>, <code>name:&lt;name&gt;</code>"
        ));
        assert!(all_command_pages.contains("dotted <code>KEY=VALUE</code> assignments"));
    }

    #[test]
    fn multiline_cli_help_uses_compact_summaries_and_block_markup() {
        let help = "Quote-free process specification.\n\nNotes: - `veto` excludes particles. - `only` keeps particles.\n\nA final paragraph\ncontinues on another line.";

        assert_eq!(
            compact_help_summary(help),
            "Quote-free process specification."
        );
        assert_eq!(
            compact_help_summary(
                "Wrap indices with a header. # Arguments - `header`: wrapper # Returns Expression"
            ),
            "Wrap indices with a header."
        );
        let rendered = render_cli_help(help);
        assert!(rendered.contains("<p>Quote-free process specification.</p>"));
        assert!(rendered.contains("<p>Notes:</p><ul>"));
        assert!(rendered.contains("<li><code>veto</code> excludes particles.</li>"));
        assert!(rendered.contains("<li><code>only</code> keeps particles.</li>"));
        assert!(rendered.contains("<p>A final paragraph continues on another line.</p>"));
    }

    #[test]
    fn cli_manpages_preserve_value_and_global_parser_facts() {
        let builder = SiteBuilder::discover().unwrap();
        let (reference, _) = builder.generated_references().unwrap();
        let page = |path: &str| {
            let command = visible_cli_commands(&reference)
                .into_iter()
                .find(|command| command.path == path)
                .unwrap();
            render_cli_command_page("GammaLoop", &reference, command)
        };

        let import_model = page("gammaLoop import model");
        assert!(import_model.contains("<div><dt>Requires equals sign</dt><dd>yes</dd></div>"));
        assert!(
            import_model
                .contains("<div><dt>Value names</dt><dd><code>SIMPLIFY_MODEL</code></dd></div>")
        );

        let inspect = page("gammaLoop inspect");
        let graph_id_anchor = generated_anchor("argument", "gammaLoop inspect::graph_id");
        let discrete_dim_anchor = generated_anchor("argument", "gammaLoop inspect::discrete_dim");
        assert!(inspect.contains("href=\"guides/events-and-observables/\""));
        assert!(inspect.contains("<div><dt>Value names</dt><dd><code>POINT</code></dd></div>"));
        assert!(inspect.contains(&format!(
            "<div><dt>Requires</dt><dd><a href=\"#{graph_id_anchor}\"><code>--graph-id</code></a></dd></div>"
        )));
        assert!(inspect.contains(&format!(
            "<div><dt>Conflicts with</dt><dd><a href=\"#{discrete_dim_anchor}\"><code>-d, --discrete-dim</code></a></dd></div>"
        )));
        assert_eq!(
            inspect
                .matches("<div><dt>Value delimiter</dt><dd><code>,</code></dd></div>")
                .count(),
            2
        );

        let approach = page("gammaLoop approach");
        assert!(approach.contains("href=\"guides/diagnostics/\""));

        let generate = page("gammaLoop generate");
        assert!(generate.contains("href=\"guides/process-generation/\""));
        assert!(generate.contains("<div><dt>Scope</dt><dd>global</dd></div>"));
        let generate_xs = page("gammaLoop generate xs");
        assert!(generate_xs.contains("<div><dt>Scope</dt><dd>inherited global</dd></div>"));

        let save_dot = page("gammaLoop save dot");
        assert!(
            save_dot.contains(
                "<div><dt>Default when value omitted</dt><dd><code>true</code></dd></div>"
            )
        );
    }

    #[test]
    fn cli_settings_use_namespace_pages_and_canonical_fragments() {
        let builder = SiteBuilder::discover().unwrap();
        let (reference, _) = builder.generated_references().unwrap();
        let setting = reference.settings.first().unwrap();
        let namespace = cli_setting_namespace(&setting.path);
        let route = cli_setting_namespace_route(namespace);
        let anchor = generated_anchor("setting", &setting.path);
        let namespaces = cli_setting_namespaces(&reference);
        let index = render_gammaloop_reference_index("GammaLoop", &reference);
        let settings_index = render_cli_settings_index("GammaLoop", &reference);
        let namespace_page =
            render_cli_settings_namespace_page("GammaLoop", &reference, namespace, &namespaces);
        let site = tempfile::tempdir().unwrap();
        let pages = write_gammaloop_reference_pages("GammaLoop", &reference, site.path()).unwrap();

        assert!(settings_index.contains(&format!("href=\"{route}\"")));
        assert!(settings_index.contains("data-reference-entry data-reference-search"));
        assert!(namespace_page.contains("<h2>Settings</h2>"));
        assert!(namespace_page.contains("href=\"guides/diagnostics/\""));
        assert!(namespace_page.contains(&format!(
            "<section class=\"reference-setting reference-setting-flat\" id=\"{anchor}\""
        )));
        assert!(!namespace_page.contains("<details"));
        assert!(!namespace_page.contains("data-reference-entry"));
        assert!(!namespace_page.contains("data-reference-search"));
        assert!(index.contains(&format!(
            "id=\"{anchor}\" href=\"{route}#{anchor}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\""
        )));
        assert!(site.path().join(&route).join("index.html").is_file());
        assert!(pages.iter().any(|page| page.route == route));
    }

    #[test]
    fn python_module_indexes_use_the_public_module_name() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "spenso")
            .unwrap();
        let component = product.python_components.first().unwrap();
        let module = python_display_module(component);
        let catalog = builder
            .component_catalog(product, ApiLanguage::Python, component)
            .unwrap();
        let rendered =
            render_python_catalog_for_module(&catalog, module, "spynso3.pyi", "0123456789abcdef");
        let routes = python_item_routes(&catalog);
        let supported = catalog_surface_items(&catalog.root, true);
        let (path, item) = supported
            .first()
            .map(|(path, item)| (path.as_slice(), *item))
            .unwrap();
        let symbol = render_python_item_page(
            &catalog,
            module,
            (path, item),
            "spynso3.pyi",
            "0123456789abcdef",
            &routes,
            (None, None),
        );

        assert_eq!(module, "symbolica.community.spenso");
        assert_ne!(module, catalog.component.package);
        assert_eq!(item.required_features.as_slice(), ["python_stubgen"]);
        assert!(rendered.contains(&format!("<p><code>{module}</code> · version")));
        assert!(!rendered.contains(&format!(
            "<p><code>{}</code> · version",
            catalog.component.package
        )));
        assert!(rendered.contains(
            "<strong>Source-build provenance:</strong> this catalog was generated with <code>python_stubgen</code>."
        ));
        assert!(!symbol.contains("api-feature-requirement"));
        assert!(symbol.contains("href=\"guides/python/\""));
    }

    #[test]
    fn python_symbol_links_cover_prose_and_inline_code_but_not_existing_links_or_examples() {
        let links = BTreeMap::from([
            (
                "EvaluationResult".to_owned(),
                "reference/python/gammaloop-python/EvaluationResult/".to_owned(),
            ),
            (
                "amp".to_owned(),
                "reference/python/gammaloop-python/amp/".to_owned(),
            ),
        ]);
        let rendered = link_python_symbols(
            "<article><p>EvaluationResult and gammaloop.EvaluationResult, but not EvaluationResults.</p><p><code>EvaluationResult</code> at https://example.test/EvaluationResult?value=amp, HTTPS://example.test/EvaluationResult, ftp://example.test/EvaluationResult, data:text/plain,EvaluationResult, docs.example/EvaluationResult, or reference/EvaluationResult/ &amp; standalone amp.</p><a href=\"manual/\">EvaluationResult</a><pre><code>EvaluationResult()</code></pre></article>",
            &links,
        );

        assert_eq!(
            rendered
                .matches("href=\"reference/python/gammaloop-python/EvaluationResult/\"")
                .count(),
            3
        );
        assert!(rendered.contains("<code><a href=\"reference/python/gammaloop-python/EvaluationResult/\">EvaluationResult</a></code>"));
        assert!(rendered.contains("EvaluationResults"));
        assert!(rendered.contains("https://example.test/EvaluationResult?value=amp"));
        assert!(rendered.contains("HTTPS://example.test/EvaluationResult"));
        assert!(rendered.contains("ftp://example.test/EvaluationResult"));
        assert!(rendered.contains("data:text/plain,EvaluationResult"));
        assert!(rendered.contains("docs.example/EvaluationResult"));
        assert!(rendered.contains("reference/EvaluationResult/"));
        assert!(rendered.contains("<a href=\"manual/\">EvaluationResult</a>"));
        assert!(rendered.contains("<pre><code>EvaluationResult()</code></pre>"));
        assert!(rendered.contains(
            "&amp; standalone <a href=\"reference/python/gammaloop-python/amp/\">amp</a>"
        ));
    }

    #[test]
    fn catalog_doc_text_renders_markup_instead_of_exposing_source_syntax() {
        let rust = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::RustMarkdown,
                "Use [`Engine`] with the [interface guide](reference/interfaces/).\n\n# Safety\nKeep `state` alive.",
            ),
            3,
        );
        assert!(rust.contains("<code>Engine</code>"));
        assert!(rust.contains("href=\"reference/interfaces/\""));
        assert!(rust.contains("<h3>Safety</h3>"));
        assert!(!rust.contains("[`Engine`]"));

        let typst = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::TypstMarkup,
                "Read the #link(\"reference/rust/\")[Rust API] for the #emph[complete] surface.",
            ),
            3,
        );
        assert!(typst.contains("href=\"reference/rust/\""));
        assert!(typst.contains("<em>complete</em>"));
        assert!(!typst.contains("#link"));

        let python = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Construct one engine.\n\n**Lifecycle:** retained.\n\n# Arguments\n- engine: active engine\n\n```python\nengine = Engine()\n```\n\n## Examples\n```python\nengine.run()\n```\n\nReuse the engine.\n\n",
            ),
            3,
        );
        assert!(python.contains("<strong>Lifecycle:</strong> retained."));
        assert!(python.contains("<h3>Arguments</h3>"));
        assert!(python.contains("<code data-lang=\"python\">engine = Engine()</code>"));
        assert!(python.contains("<code data-lang=\"python\">engine.run()</code>"));
        assert!(python.contains("<p>Reuse the engine.</p>"));
        assert!(!python.contains("```"));
    }

    #[test]
    fn markdown_style_python_returns_render_as_prose() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "# Returns\nThe simplified expression is returned without changing the input.",
            ),
            3,
        );

        assert!(rendered.contains(
            "<section class=\"api-doc-section api-doc-returns\"><h3>Returns</h3><p>The simplified expression is returned without changing the input.</p>"
        ));
        assert!(!rendered.contains("class=\"api-doc-type\""));
    }

    #[test]
    fn python_variant_sections_render_as_definitions() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Execution modes.\n\nVariants\n--------\nSingle : Execute one rewrite\nAll : Execute every rewrite",
            ),
            3,
        );

        assert!(rendered.contains("class=\"api-doc-section api-doc-variants\""));
        assert!(rendered.contains("<dt><code>Single</code></dt>"));
        assert!(rendered.contains("<dd><p>Execute one rewrite</p></dd>"));
        assert!(!rendered.contains("<p>Single : Execute one rewrite"));
    }

    #[test]
    fn numpy_examples_render_prose_and_indented_python_separately() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Examples\n--------\nIterate over grouped values::\n\n    for name, values in groups.items():\n        for value in values:\n            print(name, value)",
            ),
            3,
        );

        assert!(rendered.contains("<p>Iterate over grouped values:</p>"));
        assert!(rendered.contains(
            "<code data-lang=\"python\">for name, values in groups.items():\n    for value in values:\n        print(name, value)</code>"
        ));
        assert!(!rendered.contains("<code data-lang=\"python\">Iterate over grouped values"));
    }

    #[test]
    fn numpy_literal_examples_end_before_dedented_prose() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Examples\n--------\nEvaluate one value::\n\n    result = evaluate(1)\n    assert result\n\nThe fixture checks API plumbing.",
            ),
            3,
        );

        assert!(
            rendered
                .contains("<code data-lang=\"python\">result = evaluate(1)\nassert result</code>")
        );
        assert!(rendered.contains("<p>The fixture checks API plumbing.</p>"));
    }

    #[test]
    fn python_doctest_examples_render_as_copyable_source() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Examples\n--------\n>>> values = [\n...     1,\n...     2,\n... ]\n>>> sum(values)\n3",
            ),
            3,
        );

        assert!(rendered.contains(
            "<code data-lang=\"python\">values = [\n    1,\n    2,\n]\nsum(values)</code>"
        ));
        assert!(!rendered.contains("&gt;&gt;&gt;"));
        assert!(!rendered.contains("\n3</code>"));
    }

    #[test]
    fn python_doctest_main_remains_at_module_scope() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Examples\n--------\n>>> def main():\n...     return 1\n>>> main()\n1",
            ),
            3,
        );

        assert!(
            rendered
                .contains("<code data-lang=\"python\">def main():\n    return 1\nmain()</code>")
        );
        assert!(!rendered.contains("&gt;&gt;&gt;"));
    }

    #[test]
    fn fenced_python_doctests_render_source_without_prompts_or_output() {
        let rendered = render_doc_text(
            &alphal00p_docs_schema::DocText::new(
                DocFormat::PythonDocstring,
                "Examples\n--------\n```python\n>>> def main():\n...     return 1\n>>> main()\n1\n```",
            ),
            3,
        );

        assert!(
            rendered
                .contains("<code data-lang=\"python\">def main():\n    return 1\nmain()</code>")
        );
        assert!(!rendered.contains("&gt;&gt;&gt;"));
        assert!(!rendered.contains("\n1</code>"));
    }

    #[test]
    fn compiled_developer_typst_allows_safe_urls_images_and_prose() {
        let source = "<main><a href=\"guide?x=1&amp;y=2\">Guide</a><img src=\"data:image/png;base64,AAAA\"><pre><code>url(example)</code></pre></main>";
        assert!(validate_developer_typst_body(source).is_ok());
    }

    #[test]
    fn compiled_developer_typst_rejects_active_content() {
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
            "<main style=\"background: url(javascript:alert(1))\"></main>",
            "<img src=\"data:image/svg+xml;base64,AAAA\">",
        ] {
            assert!(validate_developer_typst_body(source).is_err(), "{source}");
        }
    }

    #[test]
    fn developer_typst_requires_one_matching_leading_title() {
        let source = "<h2>Architecture &amp; design</h2><h3>Scope</h3>";
        let title = leading_developer_title(source).unwrap();
        assert_eq!(title.title, "Architecture & design");
        assert_eq!(&source[title.range], "<h2>Architecture &amp; design</h2>");
        assert!(leading_developer_title("<p>Before</p><h2>Title</h2>").is_err());
        assert!(leading_developer_title("<h2>One</h2><h2>Two</h2>").is_err());
    }

    #[test]
    fn documentation_issue_routes_have_one_root_separator() {
        let url = documentation_issue_url(
            "GammaLoop tutorial",
            "/products/gammaloop/latest/tutorial/",
            "0123456789abcdef",
        );
        assert!(url.contains("%2Fgammaloop%2Fproducts%2Fgammaloop"));
        assert!(!url.contains("%2Fgammaloop%2F%2Fproducts"));
    }

    #[test]
    fn developer_typst_preserves_safe_compiler_head_styles() {
        let html = "<html><head><title>Note</title><style>math { display: block; }</style></head><body></body></html>";
        assert_eq!(
            extract_typst_head_styles(html).unwrap(),
            "<style>math { display: block; }</style>"
        );
        assert!(
            extract_typst_head_styles(
                "<html><head><style>@import url(https://example.test)</style></head><body></body></html>"
            )
            .is_err()
        );
    }

    #[test]
    fn prose_source_gate_rejects_new_and_parallel_markdown() {
        let temporary = tempfile::tempdir().unwrap();
        fs::write(temporary.path().join("README.md"), "Compatibility readme").unwrap();
        fs::write(temporary.path().join("new-note.md"), "New prose").unwrap();
        let mut builder = SiteBuilder::discover().unwrap();
        builder.root = temporary.path().to_path_buf();
        builder.legacy_prose = LegacyProseConfig {
            schema: LEGACY_PROSE_SCHEMA_VERSION,
            source: Vec::new(),
        };
        assert!(builder.check_prose_sources().is_err());

        builder.legacy_prose.source.push("new-note.md".into());
        assert!(builder.check_prose_sources().is_err());
        builder.legacy_prose.source.clear();
        fs::write(temporary.path().join("new-note.typ"), "= New prose").unwrap();
        assert!(builder.check_prose_sources().is_err());
    }

    #[test]
    fn prose_source_gate_rejects_format_variants_and_readme_typ() {
        let temporary = tempfile::tempdir().unwrap();
        fs::write(temporary.path().join("README.md"), "Compatibility readme").unwrap();
        let mut builder = SiteBuilder::discover().unwrap();
        builder.root = temporary.path().to_path_buf();
        builder.legacy_prose = LegacyProseConfig {
            schema: LEGACY_PROSE_SCHEMA_VERSION,
            source: Vec::new(),
        };

        for source in [
            "new-note.MD",
            "new-note.markdown",
            "new-note.mdx",
            "new-note.HTML",
            "new-note.htm",
            "new-note.shtml",
            "new-note.rst",
            "new-note.adoc",
            "new-note.org",
            "new-note.TYP",
            "readme.md",
            "readme.typ",
        ] {
            let path = temporary.path().join(source);
            fs::write(&path, "New prose").unwrap();
            assert!(builder.check_prose_sources().is_err(), "{source}");
            fs::remove_file(path).unwrap();
        }

        fs::write(temporary.path().join("README.typ"), "= Duplicate readme").unwrap();
        assert!(builder.check_prose_sources().is_err());
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
        let mut entries = Vec::new();
        append_catalog_search_with_component_name(
            &catalog,
            &catalog.root,
            &catalog.component.id,
            &mut entries,
        );
        assert!(entries.iter().any(|entry| {
            entry.href == "reference/rust/spenso_macros/derive.SimpleRepresentation.html"
        }));
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
        let mut local_paths = LocalPathIndex::new();
        assert_eq!(
            resolve_local_link(output, &source, "../../linnet/latest/", &mut local_paths).unwrap(),
            output.join("products/linnet/latest/index.html")
        );
        assert!(
            resolve_local_link(output, &source, "https://example.com", &mut local_paths).is_none()
        );
    }

    #[test]
    fn generated_links_cannot_escape_the_output_root() {
        let builder = SiteBuilder::discover().unwrap();
        let temporary = tempfile::tempdir().unwrap();
        let output = temporary.path().join("site");
        fs::create_dir_all(&output).unwrap();
        fs::write(temporary.path().join("outside.html"), "outside").unwrap();
        fs::write(
            output.join("index.html"),
            r#"<a href="../outside.html">outside</a>"#,
        )
        .unwrap();

        let error = builder
            .validate_generated_links(&output, false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("escapes generated output root"));
    }

    #[test]
    fn snapshot_developer_links_validate_against_the_same_snapshot_route() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let developer = output
            .path()
            .join("developers/architecture/example/index.html");
        let snapshot = output
            .path()
            .join("products/linnet/snapshots/v1.2.3/guides/clinnet/index.html");
        fs::create_dir_all(developer.parent().unwrap()).unwrap();
        fs::create_dir_all(snapshot.parent().unwrap()).unwrap();
        fs::write(&snapshot, r#"<h1 id="python-api">Python API</h1>"#).unwrap();
        fs::write(
            &developer,
            r#"<a href="/products/linnet/latest/guides/clinnet/#python-api">Clinnet</a>"#,
        )
        .unwrap();

        builder
            .validate_generated_links(output.path(), false, Some("v1.2.3"))
            .unwrap();

        fs::write(
            developer,
            r#"<a href="/products/linnet/latest/guides/misspelled/">broken</a>"#,
        )
        .unwrap();
        let error = builder
            .validate_generated_links(output.path(), false, Some("v1.2.3"))
            .unwrap_err();
        assert!(format!("{error:#}").contains("guides/misspelled"));
    }

    #[test]
    fn rustdoc_links_are_validated_only_when_rustdoc_is_included() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let rustdoc = output.path().join("reference/rust/example");
        fs::create_dir_all(&rustdoc).unwrap();
        fs::write(
            rustdoc.join("index.html"),
            r#"<a href="missing.html">missing</a>"#,
        )
        .unwrap();

        builder
            .validate_generated_links(output.path(), false, None)
            .unwrap();
        let error = builder
            .validate_generated_links(output.path(), true, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing.html"));
    }

    #[test]
    fn rustdoc_link_validation_uses_rendered_tags_and_source_line_ranges() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let rustdoc = output.path().join("reference/rust/example");
        let source = output.path().join("reference/rust/src/example");
        fs::create_dir_all(&rustdoc).unwrap();
        fs::create_dir_all(rustdoc.join("traits")).unwrap();
        fs::create_dir_all(&source).unwrap();
        fs::write(
            rustdoc.join("index.html"),
            r##"<script>const template = `<link href="missing/${f}">`;</script>
                <div data-href="missing-data-link.html">not a link</div>
                <a href="../src/example/lib.rs.html#10-12">source</a>
                <a HREF = "enum.Value.html#variant.Grouped.field.0">tuple variant field</a>
                <a HREF=enum.Value.html#variant.Grouped>unquoted link</a>
                <a href="traits/trait.Iterable.html#associatedtype.Data">redirected item</a>
                <a HrEf = trait.Iterable.html#tymethod.iter_flat>trait method</a>
                <a href=guide>extensionless local page</a>
                <a HREF = "crate::Span">inherited type</a>
                <a href="dispatcher#default">inherited module</a>
                <a href="javascript:void(0)">rustdoc control</a>"##,
        )
        .unwrap();
        fs::write(
            rustdoc.join("enum.Value.html"),
            r#"<section id="variant.Grouped">Grouped</section>"#,
        )
        .unwrap();
        fs::write(
            rustdoc.join("traits/trait.Iterable.html"),
            r#"<meta http-equiv="refresh" content="0;URL=../trait.Iterable.html">"#,
        )
        .unwrap();
        fs::write(
            rustdoc.join("trait.Iterable.html"),
            r#"<section id="associatedtype.Data">Data</section><section id="method.iter_flat">iter_flat</section>"#,
        )
        .unwrap();
        fs::write(
            source.join("lib.rs.html"),
            "<a href=#10 id=10 data-nosnippet>10</a><a href=#12 id=12 data-nosnippet>12</a>",
        )
        .unwrap();
        fs::write(rustdoc.join("guide"), "extensionless page").unwrap();

        builder
            .validate_generated_links(output.path(), true, None)
            .unwrap();
        let normalized = fs::read_to_string(rustdoc.join("index.html")).unwrap();
        assert!(normalized.contains("HREF = \"enum.Value.html#variant.Grouped\""));
        assert!(normalized.contains("HrEf = trait.Iterable.html#method.iter_flat"));
        assert!(normalized.contains("href=guide"));
        assert!(normalized.contains("href=\"javascript:void(0)\""));
        assert!(!normalized.contains("crate::Span"));
        assert!(!normalized.contains("href=\"dispatcher#default\""));
    }

    #[test]
    fn generated_link_validation_rejects_phantom_anchors() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        fs::write(
            output.path().join("index.html"),
            r#"<a href="target.html#description">phantom</a>"#,
        )
        .unwrap();
        fs::write(
            output.path().join("target.html"),
            r#"<meta name="description"><div data-id="description"></div><!-- <p id="description">comment</p> -->"#,
        )
        .unwrap();

        let error = builder
            .validate_generated_links(output.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing fragment"));
    }

    #[test]
    fn generated_link_validation_accepts_named_anchors() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        fs::write(
            output.path().join("index.html"),
            r#"<a href="target.html#legacy">legacy anchor</a>"#,
        )
        .unwrap();
        fs::write(
            output.path().join("target.html"),
            r#"<meta name="legacy"><a name="legacy"></a>"#,
        )
        .unwrap();

        builder
            .validate_generated_links(output.path(), false, None)
            .unwrap();
    }

    #[test]
    fn generated_link_validation_rejects_source_ranges_on_ordinary_pages() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        fs::write(
            output.path().join("index.html"),
            r#"<a href="target.html#10-12">not a source range</a>"#,
        )
        .unwrap();
        fs::write(
            output.path().join("target.html"),
            r#"<a id="10"></a><a id="12"></a>"#,
        )
        .unwrap();

        let error = builder
            .validate_generated_links(output.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing fragment"));
    }

    #[test]
    fn generated_search_links_are_validated_without_rewriting_html() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let rustdoc = output.path().join("reference/rust");
        fs::create_dir_all(&rustdoc).unwrap();
        let search_index = rustdoc.join("search-index.json");
        fs::write(rustdoc.join("index.html"), "index").unwrap();
        fs::write(
            rustdoc.join("target.html"),
            r#"<section id="method.run"></section>"#,
        )
        .unwrap();
        let search = serde_json::to_vec(&vec![SearchEntry {
            title: "Run".to_owned(),
            summary: "Run the example".to_owned(),
            href: "target.html#tymethod.run".to_owned(),
            kind: "rust-api".to_owned(),
            text: String::new(),
        }])
        .unwrap();
        fs::write(&search_index, &search).unwrap();

        let error = builder
            .validate_generated_links(output.path(), true, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing fragment"));
        assert_eq!(fs::read(search_index).unwrap(), search);
        assert_eq!(
            fs::read_to_string(rustdoc.join("index.html")).unwrap(),
            "index"
        );
    }

    #[test]
    fn repository_source_links_are_checked_against_the_workspace() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        let revision = builder.git_commit();
        fs::write(
            output.path().join("index.html"),
            format!(
                r#"<a href="https://github.com/alphal00p/gammaloop/blob/{revision}/Cargo.toml#L1">source</a>"#
            ),
        )
        .unwrap();
        builder
            .validate_generated_links(output.path(), false, None)
            .unwrap();

        fs::write(
            output.path().join("index.html"),
            format!(
                r#"<a href="https://github.com/alphal00p/gammaloop/blob/{revision}/does-not-exist.rs#L1">source</a>"#
            ),
        )
        .unwrap();
        let error = builder
            .validate_generated_links(output.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing repository source"));

        for target in [
            format!("https://github.com/alphal00p/gammaloop/blob/{revision}/docs"),
            format!("https://raw.githubusercontent.com/alphal00p/gammaloop/{revision}/docs"),
        ] {
            fs::write(
                output.path().join("index.html"),
                format!(r#"<a href="{target}">source</a>"#),
            )
            .unwrap();
            let error = builder
                .validate_generated_links(output.path(), false, None)
                .unwrap_err();
            assert!(format!("{error:#}").contains("missing repository source file"));
        }

        fs::write(
            output.path().join("index.html"),
            r#"<a href="https://github.com/alphal00p/gammaloop/blob/wrong-revision/Cargo.toml#L1">source</a>"#,
        )
        .unwrap();
        let error = builder
            .validate_generated_links(output.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("does not match documented revision"));
    }

    #[test]
    fn repository_source_paths_are_normalized_and_contained() {
        assert_eq!(
            repository_source_path(
                "https://github.com/alphal00p/gammaloop/blob/revision/crates/example/src/lib.rs#L4"
            )
            .unwrap(),
            Some((
                "revision".to_owned(),
                PathBuf::from("crates/example/src/lib.rs"),
                true
            ))
        );
        assert_eq!(
            repository_source_path(
                "https://raw.githubusercontent.com/alphal00p/gammaloop/revision/docs/README.md"
            )
            .unwrap(),
            Some(("revision".to_owned(), PathBuf::from("docs/README.md"), true))
        );
        assert!(
            repository_source_path("https://github.com/other/repository/blob/revision/file.rs")
                .unwrap()
                .is_none()
        );
        assert!(
            repository_source_path(
                "https://github.com/alphal00p/gammaloop/blob/revision/../../outside.rs"
            )
            .is_err()
        );
        assert!(
            repository_source_path("https://github.com/alphal00p/gammaloop/blob/revision/docs")
                .unwrap()
                .is_some()
        );
        assert!(
            repository_source_path("https://github.com/alphal00p/gammaloop/blob/revision/")
                .is_err()
        );
    }

    #[test]
    fn full_site_developer_links_rebase_from_root_and_nested_pages() {
        assert_eq!(page_root_prefix(""), "");
        assert_eq!(page_root_prefix("tutorial/"), "../");
        assert_eq!(page_root_prefix("reference/interfaces/"), "../../");
        assert_eq!(page_root_prefix("reference/python/component/"), "../../../");

        let developer_link = r#"<a href="../../../developers/architecture/gammaloop-architecture/">Architecture</a>"#;
        let output = Path::new("/tmp/site");
        let mut local_paths = LocalPathIndex::new();
        for (route, expected) in [
            (
                "",
                "../../../developers/architecture/gammaloop-architecture/",
            ),
            (
                "tutorial/",
                "../../../../developers/architecture/gammaloop-architecture/",
            ),
            (
                "reference/interfaces/",
                "../../../../../developers/architecture/gammaloop-architecture/",
            ),
        ] {
            let rendered = rewrite_page_links(developer_link, &page_root_prefix(route)).unwrap();
            assert!(rendered.contains(&format!("href=\"{expected}\"")));
            let source = output
                .join("products/gammaloop/latest")
                .join(route)
                .join("index.html");
            assert_eq!(
                resolve_local_link(output, &source, expected, &mut local_paths).unwrap(),
                output.join("developers/architecture/gammaloop-architecture/index.html"),
            );
        }

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
    fn pagination_keeps_generated_reference_before_version_history() {
        let builder = SiteBuilder::discover().unwrap();
        for product in &builder.registry.product {
            let pages = product_site_pages(product);
            let authored = product
                .pages
                .iter()
                .map(|page| page.route.as_str())
                .collect::<BTreeSet<_>>();
            let last_authored_reference = pages
                .iter()
                .rposition(|page| {
                    page.group == "Reference" && authored.contains(page.route.as_str())
                })
                .unwrap();
            let first_version_history = pages
                .iter()
                .position(|page| page.group == "Version history")
                .unwrap();
            let generated_reference = pages
                .iter()
                .enumerate()
                .filter(|(_, page)| {
                    page.group == "Reference" && !authored.contains(page.route.as_str())
                })
                .map(|(index, _)| index)
                .collect::<Vec<_>>();

            assert!(!generated_reference.is_empty(), "{}", product.id);
            assert!(
                generated_reference.iter().all(|index| {
                    *index > last_authored_reference && *index < first_version_history
                }),
                "{}",
                product.id
            );
        }
    }

    #[test]
    fn rendered_search_text_keeps_prose_without_chrome_or_code_blocks() {
        let html = r#"<!doctype html><html><body><nav>Install chrome noise</nav><article class="docs-article"><h1>Install Linnet</h1><p>Install Linnet for graph work.</p><nav>Draw graph navigation</nav><pre><code>cargo install noisy-example</code></pre><p>Version <code>0.6</code> is current.</p></article><footer>Footer noise</footer></body></html>"#;

        let text = rendered_article_search_text(html).unwrap();

        assert_eq!(
            text,
            "Install Linnet Install Linnet for graph work. Version 0.6 is current."
        );
    }

    #[test]
    fn site_search_uses_body_text_and_de_ranks_non_current_developer_hits() {
        let script = fs::read_to_string(
            SiteBuilder::discover()
                .unwrap()
                .root
                .join("docs/assets/site.js"),
        )
        .unwrap();
        for contract in [
            "const text = normalizeSearch(entry.text);",
            "if (text.includes(term)) score += 3;",
            "if (/command|setting/.test(kind)) score += 120;",
            "if (!/developer current /.test(kind)) score -= 120;",
        ] {
            assert!(
                script.contains(contract),
                "missing search contract: {contract}"
            );
        }
    }

    #[test]
    fn rendered_transitive_typst_headings_are_searchable() {
        let builder = SiteBuilder::discover().unwrap();
        for (product_id, page_id, html, expected_title, expected_href) in [
            (
                "linnet",
                "linnest",
                r#"<!doctype html><html><body><article><h3 id="graph-objects">Graph <code>Objects</code></h3></article></body></html>"#,
                "Graph Objects",
                "guides/linnest/#graph-objects",
            ),
            (
                "gammaloop",
                "kurvst",
                r#"<!doctype html><html><body><article><h2 id="path-wire-format">Path Wire Format</h2></article></body></html>"#,
                "Path Wire Format",
                "guides/kurvst/#path-wire-format",
            ),
            (
                "linnet",
                "clinnet-releases",
                r#"<!doctype html><html><body><article><h2 id="0-1-8-2025-12-04"><a href="https://example.test">0.1.8</a> - 2025-12-04</h2></article></body></html>"#,
                "0.1.8 - 2025-12-04",
                "version-history/clinnet/#0-1-8-2025-12-04",
            ),
        ] {
            let product = builder
                .registry
                .product
                .iter()
                .find(|product| product.id == product_id)
                .unwrap();
            let page = product
                .pages
                .iter()
                .find(|page| page.id == page_id)
                .unwrap();
            let mut entries = Vec::new();

            append_rendered_page_search(product, page, &page.route, html, &mut entries).unwrap();

            assert_eq!(entries.len(), 1);
            assert_eq!(entries[0].title, expected_title);
            assert_eq!(entries[0].href, expected_href);
            assert_eq!(
                entries[0].summary,
                format!("{} · {}", product.title, page.title)
            );
        }
    }

    #[test]
    fn typst_reference_search_uses_qualified_symbols_without_parameter_noise() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "linnet")
            .unwrap();
        let page = product
            .pages
            .iter()
            .find(|page| page.id == "typst-graph")
            .unwrap();
        let html = r#"<!doctype html><html><body><article class="docs-article"><h2 id="graph">graph</h2><h3 id="build">build</h3><h4 id="parameters">Parameters</h4><h5 id="input">input</h5><p>Input graph data.</p></article></body></html>"#;
        let mut entries = Vec::new();

        append_rendered_page_search(product, page, &page.route, html, &mut entries).unwrap();

        assert!(
            entries
                .iter()
                .any(|entry| entry.title == "graph.build" && entry.href.ends_with("#build"))
        );
        assert!(
            entries
                .iter()
                .all(|entry| entry.title != "Parameters" && entry.title != "input")
        );
        assert!(
            rendered_article_search_text(html)
                .unwrap()
                .contains("Input graph data")
        );
    }

    #[test]
    fn typst_module_symbol_anchors_are_qualified_without_qualifying_concepts() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "linnet")
            .unwrap();
        let page = SitePage::new("reference/typst/graph/", "graph module", "Typst API");
        let site = tempfile::tempdir().unwrap();
        let destination = site.path().join(&page.route);
        fs::create_dir_all(&destination).unwrap();
        fs::write(
            destination.join("index.html"),
            r##"<!doctype html><html><body><h1>graph module</h1><h2>Concepts</h2><h3>Graph Specs</h3><div><h2>graph</h2><h3>build</h3><code data-lang="typc"><span style="color: #d73948">none</span></code><p><h5>value</h5><code>int</code></p></div></body></html>"##,
        )
        .unwrap();
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Latest,
            snapshot_tag: None,
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/linnet/latest/".to_owned(),
            components: Vec::new(),
        };

        builder
            .decorate_html_page(
                product,
                &metadata,
                site.path(),
                &page,
                (None, None),
                (BuildScope::FullSite, &[]),
            )
            .unwrap();
        let rendered = fs::read_to_string(destination.join("index.html")).unwrap();

        assert!(rendered.contains("<h3 id=\"graph.build\">build</h3>"));
        assert!(rendered.contains("<h3 id=\"graph-specs\">Graph Specs</h3>"));
        assert!(!rendered.contains("id=\"graph.graph-specs\""));
        assert!(rendered.contains("<div class=\"typst-api-module\"><h2 id=\"graph\">"));
        assert!(rendered.contains("<span class=\"syntax-keyword\">none</span>"));
        assert!(!rendered.contains("style=\"color: #d73948\""));
        assert!(!rendered.contains("<p><h5"));
    }

    #[test]
    fn product_preview_developer_links_use_revision_pinned_sources() {
        let builder = SiteBuilder::discover().unwrap();
        for product in &builder.registry.product {
            let related = builder
                .registry
                .product
                .iter()
                .find(|candidate| candidate.id != product.id)
                .unwrap();
            let developer_note = builder
                .developers
                .section
                .iter()
                .flat_map(|section| &section.note)
                .find(|note| {
                    note.products
                        .iter()
                        .any(|candidate| candidate == &product.id)
                })
                .unwrap();
            let page = SitePage::new("tutorial/", "Tutorial", "Tutorial");
            let metadata = SnapshotMetadata {
                schema: SCHEMA_VERSION,
                product: &product.id,
                title: &product.title,
                channel: BuildChannel::Latest,
                snapshot_tag: None,
                git_commit: "0123456789abcdef".to_owned(),
                git_timestamp: 1,
                route: format!("products/{}/latest/", product.id),
                components: Vec::new(),
            };

            let sidebar = builder.site_sidebar(
                product,
                &metadata,
                &page,
                "../",
                BuildScope::ProductPreview,
                &[],
            );
            assert!(!sidebar.contains("developers/"), "{}", product.id);
            assert!(
                sidebar.contains(&format!(
                    "href=\"{PUBLISHED_DOCS_ROOT}/citations/#{}\">Cite {}</a>",
                    product.id, product.title
                )),
                "{}",
                product.id
            );
            let options =
                builder.product_options(product, &metadata, "../", BuildScope::ProductPreview);
            assert_eq!(options.matches("<option ").count(), 1, "{}", product.id);
            assert!(options.contains(&format!(">{}</option>", product.title)));

            let body = format!(
                "<a href=\"../../../{}/latest/#overview\">Related</a><a href=\"../guides/\">Local</a><a href=\"../../../../developers/architecture/{}/\">Developer</a>",
                related.id, developer_note.id,
            );
            let rendered = builder
                .rewrite_product_preview_links(&body, &metadata, &page)
                .unwrap();
            assert!(
                rendered.contains(&format!(
                    "href=\"{PUBLISHED_DOCS_ROOT}/products/{}/latest/#overview\"",
                    related.id
                )),
                "{}",
                product.id
            );
            assert!(rendered.contains("href=\"../guides/\""), "{}", product.id);
            assert!(
                rendered.contains(&format!(
                    "href=\"https://github.com/alphal00p/gammaloop/blob/0123456789abcdef/{}\"",
                    developer_note.source.to_string_lossy().replace('\\', "/"),
                )),
                "{}",
                product.id,
            );

            let snapshot_link = format!(
                "<a href=\"https://github.com/alphal00p/gammaloop/blob/0123456789abcdef/{}\">Developer at this revision</a>",
                developer_note.source.to_string_lossy().replace('\\', "/"),
            );
            assert_eq!(
                builder
                    .rewrite_product_preview_links(&snapshot_link, &metadata, &page)
                    .unwrap(),
                snapshot_link,
            );
        }
    }

    #[test]
    fn every_product_preview_entrypoint_validates_from_an_empty_output() {
        let builder = SiteBuilder::discover().unwrap();
        for product in &builder.registry.product {
            let output = tempfile::tempdir().unwrap();
            let product_root = output
                .path()
                .join("products")
                .join(&product.id)
                .join("latest");
            fs::create_dir_all(&product_root).unwrap();
            fs::write(
                product_root.join("index.html"),
                "<!doctype html><html><body><h1 id=\"overview\">Overview</h1></body></html>",
            )
            .unwrap();
            fs::write(
                product_root.join("search-index.json"),
                serde_json::to_vec(&vec![SearchEntry {
                    title: product.title.clone(),
                    summary: product.tagline.clone(),
                    href: "index.html#overview".to_owned(),
                    kind: "product".to_owned(),
                    text: String::new(),
                }])
                .unwrap(),
            )
            .unwrap();

            builder
                .write_product_preview(output.path(), product, BuildChannel::Latest, None)
                .unwrap();
            builder
                .validate_generated_links(output.path(), false, None)
                .unwrap();

            let entrypoint = fs::read_to_string(output.path().join("index.html")).unwrap();
            assert!(
                entrypoint.contains(&format!("products/{}/latest/", product.id)),
                "{}",
                product.id
            );
            assert!(
                !entrypoint.contains("portal-project-card"),
                "{}",
                product.id
            );
            assert!(!output.path().join("developers").exists(), "{}", product.id);
            assert!(output.path().join(".nojekyll").is_file());
        }
    }

    #[test]
    fn product_preview_validation_rejects_a_broken_search_fragment() {
        let builder = SiteBuilder::discover().unwrap();
        let product = &builder.registry.product[0];
        let output = tempfile::tempdir().unwrap();
        let product_root = output
            .path()
            .join("products")
            .join(&product.id)
            .join("latest");
        fs::create_dir_all(&product_root).unwrap();
        fs::write(
            product_root.join("index.html"),
            "<!doctype html><html><body><h1 id=\"overview\">Overview</h1></body></html>",
        )
        .unwrap();
        fs::write(
            product_root.join("search-index.json"),
            serde_json::to_vec(&vec![SearchEntry {
                title: product.title.clone(),
                summary: product.tagline.clone(),
                href: "index.html#missing".to_owned(),
                kind: "product".to_owned(),
                text: String::new(),
            }])
            .unwrap(),
        )
        .unwrap();
        builder
            .write_product_preview(output.path(), product, BuildChannel::Latest, None)
            .unwrap();

        let error = builder
            .validate_generated_links(output.path(), false, None)
            .unwrap_err();
        assert!(format!("{error:#}").contains("missing fragment"));
    }

    #[test]
    fn python_callable_signatures_use_call_syntax_without_bound_receivers() {
        assert_eq!(
            python_human_callable_signature(
                "def run(self, values: tuple[int, int] = (1, 2)) -> str:"
            ),
            "run(values: tuple[int, int] = (1, 2)) -> str"
        );
        assert_eq!(
            python_human_callable_signature("async def load(cls, path: str) -> Engine:"),
            "async load(path: str) -> Engine"
        );
        assert_eq!(
            python_human_callable_signature("def version() -> str:"),
            "version() -> str"
        );
    }

    #[test]
    fn python_constants_keep_their_annotated_signature_and_label() {
        let mut item = DocItem::new(
            "AUTO",
            "AUTO",
            "AUTO",
            alphal00p_docs_schema::DocItemKind::PythonConstant,
        );
        item.signature = Some("AUTO: Auto".to_owned());

        assert_eq!(api_item_kind_label(item.kind), "Constant");
        assert_eq!(
            python_item_display_signature(&item).as_deref(),
            Some("AUTO: Auto")
        );
        let class = DocItem::new(
            "Auto",
            "Auto",
            "Auto",
            alphal00p_docs_schema::DocItemKind::PythonClass,
        );
        assert_eq!(api_item_anchor(&[], &item), "auto-constant");
        assert_eq!(api_item_anchor(&[], &class), "auto");
        assert_eq!(
            python_item_disambiguated_segment(&item, &item.name),
            "AUTO-constant"
        );
        let function = DocItem::new(
            "edge",
            "edge",
            "edge",
            alphal00p_docs_schema::DocItemKind::PythonFunction,
        );
        let edge = DocItem::new(
            "Edge",
            "Edge",
            "Edge",
            alphal00p_docs_schema::DocItemKind::PythonClass,
        );
        assert_eq!(api_item_anchor(&[], &function), "edge-function");
        assert_eq!(api_item_anchor(&[], &edge), "edge");
        assert_eq!(
            python_item_disambiguated_segment(&function, &function.name),
            "edge-function"
        );
    }

    #[test]
    fn python_symbol_pages_and_search_use_canonical_routes() {
        use alphal00p_docs_schema::{
            DocComponent, DocExample, DocMemberKind, DocProduct, DocText, SCHEMA_VERSION,
        };

        let mut root = DocScope::new("module", "Module exports");
        let mut internal = DocItem::new(
            "InternalEngine",
            "InternalEngine",
            "InternalEngine",
            alphal00p_docs_schema::DocItemKind::PythonClass,
        );
        internal.supported = false;
        internal.summary = Some("Internal binding machinery.".to_owned());
        internal
            .members
            .push(DocMember::new("debug", DocMemberKind::Method));
        root.define_item(internal).unwrap();
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
        item.required_features = vec!["binding".to_owned(), "accelerated".to_owned()];
        let mut constructor = DocMember::new("__new__", DocMemberKind::AssociatedFunction);
        constructor.signature = Some("def __new__(cls, state: str = None) -> Engine:".to_owned());
        constructor.docs = Some(DocText::new(
            DocFormat::PythonDocstring,
            "Create an engine.\n\nParameters\n----------\nstate : str, optional\n    Initial state.",
        ));
        let mut state = DocMember::new("state", DocMemberKind::Parameter);
        state.signature = Some("str".to_owned());
        state.default = Some("None".to_owned());
        state.docs = Some(DocText::new(DocFormat::PythonDocstring, "Initial state."));
        constructor.members.push(state);
        item.members.push(constructor);
        let mut method = DocMember::new("run", DocMemberKind::Method);
        let mut integer_overload = DocMember::new("run", DocMemberKind::Overload);
        integer_overload.signature = Some("def run(self, value: int) -> LateSupported:".to_owned());
        integer_overload.docs = Some(DocText::new(
            DocFormat::PythonDocstring,
            "Run one value.\n\nParameters\n----------\nvalue : int\n    Input value.\n\nReturns\n-------\nLateSupported\n    Processed value.\n\nRaises\n------\nValueError\n    If the value is invalid.\n\nExamples\n--------\n>>> engine.run(1)\n1\n\nNotes\n-----\nThe engine retains its state.",
        ));
        let mut value = DocMember::new("value", DocMemberKind::Parameter);
        value.signature = Some("int".to_owned());
        value.docs = Some(DocText::new(DocFormat::PythonDocstring, "Input value."));
        integer_overload.members.push(value);
        method.members.push(integer_overload);
        let mut overload = DocMember::new("run", DocMemberKind::Overload);
        overload.signature = Some("def run(self, value: str) -> str:".to_owned());
        method.members.push(overload);
        item.members.push(method);
        let mut property_getter = DocMember::new("global_data", DocMemberKind::Getter);
        property_getter.signature = Some("def global_data(self) -> str:".to_owned());
        property_getter.docs = Some(DocText::new(
            DocFormat::PythonDocstring,
            "Current global data.",
        ));
        item.members.push(property_getter);
        let mut property_setter = DocMember::new("global_data", DocMemberKind::Setter);
        property_setter.signature = Some("def global_data(self, value: str) -> None:".to_owned());
        let mut property_value = DocMember::new("value", DocMemberKind::Parameter);
        property_value.signature = Some("str".to_owned());
        property_setter.members.push(property_value);
        item.members.push(property_setter);
        let mut read_only = DocMember::new("status", DocMemberKind::Getter);
        read_only.signature = Some("def status(self) -> str:".to_owned());
        item.members.push(read_only);
        let mut from_value = DocMember::new("from_value", DocMemberKind::AssociatedFunction);
        from_value.signature = Some("def from_value(cls, value: int) -> Engine:".to_owned());
        item.members.push(from_value);
        let mut version = DocMember::new("version", DocMemberKind::AssociatedFunction);
        version.signature = Some("def version() -> str:".to_owned());
        item.members.push(version);
        root.define_item(item).unwrap();
        let mut curated = DocScope::new("curated", "Curated helpers");
        let mut later_supported = DocItem::new(
            "LateSupported",
            "LateSupported",
            "LateSupported",
            alphal00p_docs_schema::DocItemKind::PythonClass,
        );
        later_supported.summary = Some("A supported export in a later scope.".to_owned());
        later_supported.examples.push(DocExample::new(
            "Example usage",
            "python",
            "LateSupported()",
        ));
        curated.define_item(later_supported).unwrap();
        root.define_scope(curated).unwrap();
        let mut component = DocComponent::new(
            "example-python",
            "example",
            "Example Python",
            "1.0.0",
            ApiLanguage::Python,
        );
        component.features.push("binding".to_owned());
        let catalog = DocCatalog {
            schema_version: SCHEMA_VERSION,
            product: DocProduct::new("example", "Example"),
            component,
            root,
        };

        let module = "public.example";
        let routes = python_item_routes(&catalog);
        let module_page = render_python_catalog_for_module(
            &catalog,
            module,
            "example-python.pyi",
            "0123456789abcdef",
        );
        let supported = catalog_surface_items(&catalog.root, true);
        let (engine_path, engine) = supported
            .iter()
            .find(|(_, item)| item.name == "Engine")
            .map(|(path, item)| (path.as_slice(), *item))
            .unwrap();
        let engine_anchor = api_item_anchor(engine_path, engine);
        let engine_route = &routes[&engine_anchor];
        let engine_page = render_python_item_page(
            &catalog,
            module,
            (engine_path, engine),
            "example-python.pyi",
            "0123456789abcdef",
            &routes,
            (None, None),
        );
        let (late_path, late) = supported
            .iter()
            .find(|(_, item)| item.name == "LateSupported")
            .map(|(path, item)| (path.as_slice(), *item))
            .unwrap();
        let late_anchor = api_item_anchor(late_path, late);
        let late_route = &routes[&late_anchor];
        let late_page = render_python_item_page(
            &catalog,
            module,
            (late_path, late),
            "example-python.pyi",
            "0123456789abcdef",
            &routes,
            (None, None),
        );

        assert_eq!(engine_route, "reference/python/example-python/Engine/");
        assert_eq!(
            late_route,
            "reference/python/example-python/curated/LateSupported/"
        );
        assert!(module_page.contains("data-reference-filter"));
        assert!(module_page.contains("data-reference-entry data-reference-search"));
        assert!(module_page.contains("id=\"python-inventory-title\">API inventory"));
        assert!(!module_page.contains("<progress"));
        assert!(!module_page.contains("full descriptions"));
        assert_eq!(module_page.matches("class=\"api-symbol-card\"").count(), 2);
        assert_eq!(module_page.matches("class=\"api-symbol-meta\"").count(), 1);
        assert!(module_page.contains(">5 members</span>"));
        assert!(!module_page.contains("0 members"));
        assert!(module_page.contains(&format!("href=\"{engine_route}\"")));
        assert!(module_page.contains(&format!("href=\"{late_route}\"")));
        assert!(!module_page.contains("InternalEngine"));
        assert!(!module_page.contains("<details"));
        assert!(module_page.contains("href=\"quickstart/python/\""));
        assert!(module_page.contains("href=\"reference/interfaces/\""));
        assert!(module_page.contains("href=\"reference/python/example-python.pyi\""));
        assert!(module_page.contains(
            "<strong>Source-build provenance:</strong> this catalog was generated with <code>binding</code>."
        ));
        assert!(module_page.contains(&format!(
            "<a class=\"legacy-reference-anchor\" id=\"{engine_anchor}\" href=\"{engine_route}\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\"></a>"
        )));
        assert!(module_page.contains(&format!(
            "id=\"engine-run-method\" href=\"{engine_route}#engine-run-method\" data-reference-redirect tabindex=\"-1\" aria-hidden=\"true\""
        )));

        assert!(engine_page.contains("Engine(state: str = None)"));
        assert!(!engine_page.contains("class Engine:"));
        assert!(engine_page.contains(
            "<strong>Additional source-build feature:</strong> This binding is generated when <code>accelerated</code> is enabled."
        ));
        assert!(!engine_page.contains("<code>binding</code> is enabled"));
        assert!(!engine_page.contains("data-reference-entry"));
        assert!(!engine_page.contains("data-reference-search"));
        assert!(engine_page.contains("Construct one engine and reuse it."));
        assert!(engine_page.contains("href=\"quickstart/\">Getting started</a>"));
        assert_eq!(
            engine_page.matches("Runs the documented workflow.").count(),
            1
        );
        assert!(
            engine_page
                .contains("<section class=\"api-member api-member-flat\" id=\"engine-run-method\"")
        );
        assert!(!engine_page.contains("<details"));
        assert!(engine_page.contains(">Method</span>"));
        assert!(engine_page.contains(">Class method</span>"));
        assert!(engine_page.contains(">Static method</span>"));
        assert!(!engine_page.contains(">Class or static method</span>"));
        assert!(engine_page.contains(">Property · read/write</span>"));
        assert!(engine_page.contains(">Property · read-only</span>"));
        assert!(engine_page.contains("run(value: int)"));
        assert!(engine_page.contains("from_value(value: int)"));
        assert!(engine_page.contains("version() -&gt; str"));
        assert!(!engine_page.contains("def run(self"));
        assert!(!engine_page.contains("from_value(cls"));
        assert!(engine_page.contains("global_data: str"));
        assert!(engine_page.contains("status: str"));
        assert!(!engine_page.contains("def global_data"));
        assert!(!engine_page.contains("def status"));
        assert!(engine_page.contains(&format!(
            "<p class=\"api-doc-type\"><code><a href=\"{late_route}\">LateSupported</a></code></p>"
        )));
        assert!(!engine_page.contains(&format!("href=\"{engine_route}\">Engine</a>")));
        assert!(engine_page.contains("class=\"api-parameter-table\""));
        assert!(engine_page.contains("<td><code>value</code></td><td><code>int</code></td>"));
        assert!(engine_page.contains("class=\"api-doc-section api-doc-returns\""));
        assert!(engine_page.contains("class=\"api-doc-section api-doc-raises\""));
        assert!(engine_page.contains("class=\"api-doc-section api-doc-examples\""));
        assert!(engine_page.contains("class=\"api-doc-section api-doc-notes\""));
        assert_eq!(
            engine_page
                .matches("id=\"engine-run-method-overload\"")
                .count(),
            1
        );
        assert_eq!(
            engine_page
                .matches("id=\"engine-run-method-overload-2\"")
                .count(),
            1
        );
        assert_eq!(engine_page.matches("<h4>Overloads</h4>").count(), 1);
        assert_eq!(engine_page.matches("<h5>Overload ").count(), 2);
        assert_eq!(
            engine_page
                .matches("<h3><code>global_data</code></h3>")
                .count(),
            1
        );
        assert_eq!(
            engine_page
                .matches("id=\"engine-global-data-setter\"")
                .count(),
            1
        );
        assert!(
            engine_page.contains(
                "<div class=\"api-member-anchor-alias\" id=\"engine-global-data-setter\">"
            )
        );
        assert!(engine_page.contains("api-property-setter-signature"));
        assert!(engine_page.contains("<h2>Constructor</h2>"));
        assert!(
            engine_page.find("<h2>Constructor</h2>").unwrap()
                < engine_page.find("Member overview</h2>").unwrap()
        );
        assert!(engine_page.contains("id=\"engine-new-associatedfunction\""));
        assert!(!engine_page.contains("<code>__new__</code>"));
        assert_eq!(engine_page.matches("Initial state.").count(), 1);
        assert!(late_page.contains("<h2>Example usage</h2>"));
        assert!(late_page.contains("LateSupported()"));

        let mut entries = Vec::new();
        append_catalog_search_with_component_name(&catalog, &catalog.root, module, &mut entries);
        assert!(entries.iter().any(|entry| {
            entry.title == "public.example.Engine" && entry.href == *engine_route
        }));
        assert!(entries.iter().any(|entry| {
            entry.title == "public.example.LateSupported" && entry.href == *late_route
        }));
        assert!(
            entries
                .iter()
                .any(|entry| { entry.href == format!("{engine_route}#engine-run-method") })
        );
        assert!(entries.iter().any(|entry| {
            entry.href == format!("{engine_route}#engine-new-associatedfunction")
        }));
        assert!(
            entries.iter().any(|entry| {
                entry.href == format!("{engine_route}#engine-run-method-overload")
            })
        );
        assert!(
            entries.iter().any(|entry| {
                entry.href == format!("{engine_route}#engine-run-method-overload-2")
            })
        );
        assert_eq!(
            entries
                .iter()
                .filter(|entry| entry.title.ends_with("Engine.run"))
                .count(),
            1
        );
        assert_eq!(
            entries
                .iter()
                .filter(|entry| entry.title.contains("Engine.run · overload"))
                .count(),
            2
        );
        assert!(
            entries
                .iter()
                .any(|entry| { entry.href == format!("{engine_route}#engine-global-data-getter") })
        );
        assert!(
            !entries
                .iter()
                .any(|entry| { entry.href == format!("{engine_route}#engine-global-data-setter") })
        );
        assert_eq!(
            entries
                .iter()
                .filter(|entry| entry.title.ends_with("Engine.global_data"))
                .count(),
            1
        );
        assert!(entries.iter().all(|entry| {
            !entry.title.contains("InternalEngine") && !entry.href.contains("internalengine")
        }));
    }

    #[test]
    fn python_member_fragments_have_scroll_and_target_styles() {
        let builder = SiteBuilder::discover().unwrap();
        let css = fs::read_to_string(builder.root.join("docs/assets/site.css")).unwrap();

        assert!(css.contains(
            ".api-member, .api-overload, .api-member-anchor-alias { scroll-margin-top: calc(var(--header-height) + 1rem); }"
        ));
        assert!(css.contains(".api-member-anchor-alias:target > .api-member"));
        assert!(css.contains(
            ".reference-coverage-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(14rem, 1fr));"
        ));
        assert!(css.contains("math[display=\"block\"] { display: block math;"));
        assert!(
            css.contains(".reference-table-wrap > table { display: table; overflow: visible; }")
        );
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
    fn copying_generated_trees_requires_a_directory_source() {
        let temporary = tempfile::tempdir().unwrap();
        let source = temporary.path().join("source.html");
        let destination = temporary.path().join("destination");
        fs::write(&source, "not a generated tree").unwrap();

        let error = copy_tree(&source, &destination).unwrap_err();

        assert!(format!("{error:#}").contains("documentation source is not a directory"));
        assert!(!destination.exists());
    }

    #[test]
    fn replacing_generated_trees_rejects_overlapping_paths() {
        let temporary = tempfile::tempdir().unwrap();
        let destination = temporary.path().join("site");
        let source = destination.join("candidate");
        fs::create_dir_all(&source).unwrap();
        fs::write(destination.join("keep.html"), "keep").unwrap();

        let error = replace_generated_tree(&source, &destination).unwrap_err();

        assert!(format!("{error:#}").contains("overlapping documentation trees"));
        assert_eq!(
            fs::read_to_string(destination.join("keep.html")).unwrap(),
            "keep"
        );
        assert!(source.is_dir());
    }

    #[cfg(unix)]
    #[test]
    fn safe_outputs_reject_a_symlink_alias_to_the_workspace_root() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path().join("workspace");
        let alias = temporary.path().join("alias");
        fs::create_dir(&root).unwrap();
        symlink(&root, &alias).unwrap();

        let error = ensure_safe_output(&root, &alias).unwrap_err();

        assert!(format!("{error:#}").contains("cannot contain the workspace root"));
        assert!(root.is_dir());
    }

    #[cfg(unix)]
    #[test]
    fn replacing_generated_trees_resolves_symlinked_ancestors() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let work = temporary.path().join("work");
        let source = work.join("site");
        let alias = temporary.path().join("alias");
        fs::create_dir_all(&source).unwrap();
        fs::write(source.join("keep.html"), "keep").unwrap();
        symlink(&work, &alias).unwrap();

        let error = replace_generated_tree(&source, &alias.join("site")).unwrap_err();

        assert!(format!("{error:#}").contains("overlapping documentation trees"));
        assert_eq!(
            fs::read_to_string(source.join("keep.html")).unwrap(),
            "keep"
        );
    }

    #[cfg(unix)]
    #[test]
    fn copying_generated_trees_rejects_destination_symlinks() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let source = temporary.path().join("source");
        let destination = temporary.path().join("destination");
        let external = temporary.path().join("external");
        fs::create_dir_all(source.join("developers")).unwrap();
        fs::create_dir_all(&destination).unwrap();
        fs::create_dir_all(&external).unwrap();
        fs::write(source.join("developers/index.html"), "generated").unwrap();
        fs::write(external.join("keep.html"), "keep").unwrap();
        symlink(&external, destination.join("developers")).unwrap();

        let error = copy_tree(&source, &destination).unwrap_err();

        assert!(format!("{error:#}").contains("unsupported documentation artifact"));
        assert_eq!(
            fs::read_to_string(external.join("keep.html")).unwrap(),
            "keep"
        );
        assert!(!external.join("index.html").exists());
    }

    #[cfg(unix)]
    #[test]
    fn immutable_directory_comparison_rejects_symlinks() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let generated = temporary.path().join("generated");
        let existing = temporary.path().join("existing");
        fs::create_dir_all(&generated).unwrap();
        fs::create_dir_all(&existing).unwrap();
        fs::write(generated.join("index.html"), "same").unwrap();
        fs::write(existing.join("index.html"), "same").unwrap();
        symlink(generated.join("index.html"), existing.join("extra.html")).unwrap();

        let error = directories_equal(&generated, &existing).unwrap_err();

        assert!(format!("{error:#}").contains("unsupported documentation artifact"));
    }

    #[test]
    fn checked_in_registry_is_valid() {
        SiteBuilder::discover().unwrap().check().unwrap();
    }

    #[test]
    fn every_product_requires_a_task_chooser_and_each_supported_interface() {
        let mut builder = SiteBuilder::discover().unwrap();
        let quickstart = builder.registry.product[0]
            .pages
            .iter_mut()
            .find(|page| page.id == "quickstart")
            .unwrap();
        quickstart.id = "getting-started".to_owned();

        let error = builder.check().unwrap_err();
        assert!(format!("{error:#}").contains("has no quickstart page"));

        let mut builder = SiteBuilder::discover().unwrap();
        let quickstart = builder.registry.product[0]
            .pages
            .iter_mut()
            .find(|page| page.id == "quickstart")
            .unwrap();
        quickstart.route = "start/".to_owned();

        let error = builder.check().unwrap_err();
        assert!(
            format!("{error:#}").contains("quickstart chooser must use an outcome-specific title")
        );

        let mut builder = SiteBuilder::discover().unwrap();
        builder.registry.product[0]
            .pages
            .iter_mut()
            .find(|page| page.id == "quickstart-cli")
            .unwrap()
            .id = "getting-started-cli".to_owned();

        let error = builder.check().unwrap_err();
        assert!(format!("{error:#}").contains("gammaloop has no CLI quickstart"));
    }

    #[test]
    fn site_assets_replace_read_only_generated_files() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        builder.write_site_assets(output.path()).unwrap();
        let font = output.path().join("assets/STIXTwoMath-Regular.woff2");
        let mut permissions = fs::metadata(&font).unwrap().permissions();
        permissions.set_readonly(true);
        fs::set_permissions(&font, permissions).unwrap();

        builder.write_site_assets(output.path()).unwrap();

        assert_eq!(
            fs::read(font).unwrap(),
            fs::read(builder.root.join("docs/assets/STIXTwoMath-Regular.woff2")).unwrap()
        );
    }

    #[test]
    fn pdf_manual_chapter_order_must_match_the_page_registry() {
        let mut builder = SiteBuilder::discover().unwrap();
        builder.registry.product[0].pages.swap(0, 1);

        let error = builder.check().unwrap_err();
        assert!(format!("{error:#}").contains("PDF manual chapter order differs"));
    }

    #[test]
    fn current_developer_notes_require_a_verified_code_scope() {
        let mut builder = SiteBuilder::discover().unwrap();
        let note = builder
            .developers
            .section
            .iter_mut()
            .flat_map(|section| &mut section.note)
            .find(|note| note.lifecycle == "current")
            .expect("checked-in registry has a current developer note");
        note.scope.clear();

        let error = builder.check().unwrap_err();
        assert!(format!("{error:#}").contains("current note must declare a verified code scope"));
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
        assert!(html.contains("href=\"about/\""));
        assert!(html.contains("href=\"talks/\""));
        assert!(html.contains("href=\"developers/\""));
        assert!(html.contains("href=\"publications/\""));
        assert!(
            html.contains("href=\"products/gammaloop/latest/quickstart/\">Start with GammaLoop")
        );
        assert!(html.contains(
            "class=\"portal-graph-field\" role=\"img\" aria-label=\"A jumble of Feynman graphs rendered by GammaLoop from real process and test data\""
        ));
        assert_eq!(
            html.matches("class=\"portal-process-graph\"").count(),
            PORTAL_GRAPH_IDS.len()
        );
        assert_eq!(
            html.matches("class=\"portal-graph-theme portal-graph-theme-")
                .count(),
            PORTAL_GRAPH_IDS.len() * 2
        );
        assert!(html.contains(
            "class=\"portal-wordmark\" aria-label=\"αLoop collaboration mark\" role=\"img\""
        ));
        assert!(!html.contains(
            "class=\"portal-hero-art\" aria-label=\"αLoop collaboration mark\" role=\"img\""
        ));
        assert!(!html.contains("class=\"portal-graph-card\""));
        assert!(!html.contains("portal-topology-field"));
        assert_eq!(html.matches("class=\"portal-funding\"").count(), 1);
        assert!(html.contains("Publicly funded research"));
        assert!(html.contains(&builder.portal.funding));
        assert!(html.contains(&builder.portal.funding_url));
        assert!(html.contains(
            "rel=\"icon\" type=\"image/svg+xml\" href=\"assets/local-unitarity-light.svg\""
        ));
        assert!(html.contains(
            "href=\"assets/local-unitarity-dark.svg\" media=\"(prefers-color-scheme: dark)\""
        ));
        assert_eq!(html.matches("portal-card-cite").count(), 5);
        for removed in [
            "class=\"portal-pillars\"",
            "class=\"portal-facts\"",
            "id=\"developers\"",
            "id=\"people\"",
            "id=\"affiliations\"",
            "href=\"#affiliations\"",
        ] {
            assert!(
                !html.contains(removed),
                "unexpected landing-page copy: {removed}"
            );
        }
        assert!(!html.contains("Product manual"));
        assert!(!html.contains("scientific-computing products"));

        for asset in [
            "about-double-triangle-light.svg",
            "about-double-triangle-dark.svg",
            "about-local-unitarity-equation-light.svg",
            "about-local-unitarity-equation-dark.svg",
            "local-unitarity-light.svg",
            "local-unitarity-dark.svg",
            "gammalooplogo-light.svg",
            "gammalooplogo-dark.svg",
            "spensologo.svg",
            "STIXTwoMath-Regular.woff2",
            "STIX-Two-OFL.txt",
            "publications.json",
        ] {
            assert!(output.path().join("assets").join(asset).is_file());
        }
        let mut svg_assets = portal_graph_assets()
            .map(|graph| output.path().join("assets/graphs").join(graph))
            .collect::<Vec<_>>();
        svg_assets.extend(
            [
                "about-double-triangle-light.svg",
                "about-double-triangle-dark.svg",
                "about-local-unitarity-equation-light.svg",
                "about-local-unitarity-equation-dark.svg",
                "local-unitarity-light.svg",
                "local-unitarity-dark.svg",
                "gammalooplogo-light.svg",
                "gammalooplogo-dark.svg",
                "spensologo.svg",
            ]
            .map(|asset| output.path().join("assets").join(asset)),
        );
        for asset in svg_assets {
            let svg = fs::read_to_string(&asset).unwrap();
            assert!(
                svg.contains("<svg"),
                "missing SVG root in {}",
                asset.display()
            );
            assert!(
                svg.contains("viewBox="),
                "missing SVG viewBox in {}",
                asset.display()
            );
            for forbidden in ["<script", "<foreignObject", "<image", "href=\"http"] {
                assert!(
                    !svg.contains(forbidden),
                    "unsafe or external SVG content in {}: {forbidden}",
                    asset.display()
                );
            }
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

        let about = fs::read_to_string(output.path().join("about/index.html")).unwrap();
        assert!(about.contains("Precision is another path to discovery"));
        assert!(about.contains("Schematic Local Unitarity cross-section"));
        assert!(about.contains("../assets/about-double-triangle-light.svg"));
        assert!(about.contains("../assets/about-double-triangle-dark.svg"));
        assert!(about.contains("../assets/about-local-unitarity-equation-light.svg"));
        assert!(about.contains("../assets/about-local-unitarity-equation-dark.svg"));
        assert!(about.contains("class=\"about-equation-formula\" role=\"img\""));
        assert!(!about.contains("<span>dσ/d𝒪</span>"));
        assert_eq!(
            about.matches("class=\"about-pillar\"").count(),
            builder.portal.pillar.len()
        );
        assert_eq!(
            about.matches("class=\"about-affiliation\"").count(),
            builder.portal.affiliation.len()
        );
        assert!(about.contains(&builder.portal.funding));

        let talks = fs::read_to_string(output.path().join("talks/index.html")).unwrap();
        assert_eq!(
            talks.matches("class=\"talk-card\"").count(),
            builder.talks.talk.len()
        );
        assert!(talks.contains("The GammaLoop ecosystem"));
        assert!(talks.contains("href=\"../people/#lucien-huber\""));
        assert!(talks.contains("aria-current=\"page\">Talks"));

        let css = fs::read_to_string(output.path().join("assets/site.css")).unwrap();
        assert!(css.contains(".publication-card:nth-child(2n of :not([hidden]))"));
        assert!(!css.contains(".publication-card:nth-child(2n) {"));
        assert!(css.contains(".people-card-portrait { object-position: left center; }"));
        assert!(css.contains("#ben-ruijl > .people-card-portrait"));
        assert!(css.contains(".product-logo-gammaloop { aspect-ratio: 2.016 / 1;"));
        assert!(css.contains("background-size: contain;"));
        assert!(!css.contains("background-size: 137.2% auto;"));
        assert!(css.contains(".product-logo-spenso { aspect-ratio: 637 / 189;"));
        assert!(css.contains(".portal-graph-field"));
        assert!(css.contains("grid-template-columns: clamp(15rem, 20vw, 18rem)"));
        assert!(css.contains(".portal-process-graph:nth-child(11)"));
        assert!(css.contains(":root[data-theme=\"dark\"] .portal-graph-theme-light"));
        assert!(css.contains(":root:not([data-theme]) {"));
        assert!(!css.contains(".portal-process-graph > img { display: block;"));
        assert!(!css.contains(".portal-graph-card"));

        let renderer =
            fs::read_to_string(builder.root.join("scripts/render-docs-svg-assets.sh")).unwrap();
        assert!(renderer.contains("typst compile"));
        for forbidden in ["cargo run", "save dot", "linnet draw", "dot -", "sed -i"] {
            assert!(
                !renderer.contains(forbidden),
                "website SVG renderer uses external generation: {forbidden}"
            );
        }
        for graph in PORTAL_GRAPH_IDS {
            let source = builder
                .root
                .join("docs/assets/typst/portal-graphs/graphs")
                .join(format!("portal-graph-{graph}.typ"));
            let typst = fs::read_to_string(&source).unwrap();
            assert!(
                typst.contains("#render(") && typst.contains("read(") && typst.contains(".dot\""),
                "missing editable Linnest source for {graph}"
            );
        }
        let about_graph = fs::read_to_string(
            builder
                .root
                .join("docs/assets/typst/about/double-triangle.typ"),
        )
        .unwrap();
        assert!(about_graph.contains("tests/resources/graphs/double_triangle.dot"));
        assert!(about_graph.contains("#render(input"));
        assert!(about_graph.contains("cut_curve=blue"));
        assert!(about_graph.contains("cut_curve=red"));
        assert!(about_graph.contains("cut-curves: true"));
        assert!(!about_graph.contains("cut_marker"));
        assert!(about_graph.contains("v3 [label=\\\"\\\", pos="));
        let about_equation = fs::read_to_string(
            builder
                .root
                .join("docs/assets/typst/about/local-unitarity-equation.typ"),
        )
        .unwrap();
        for term in [
            "frac(dif sigma, dif cal(O))",
            "product_(i = 1)^(n_\"loop\") dif^3 bold(k)_i",
            "I_(Gamma,c)^upright(\"LU\")",
            "delta(cal(O)(c, bold(k)_i))",
        ] {
            assert!(about_equation.contains(term), "About equation lost {term}");
        }

        let shared_layout = fs::read_to_string(
            builder
                .root
                .join("assets/embedded/drawing/templates/layout-core.typ"),
        )
        .unwrap();
        let canonical_layout = fs::read_to_string(
            builder
                .root
                .join("assets/embedded/drawing/templates/layout.typ"),
        )
        .unwrap();
        let portal_layout = fs::read_to_string(
            builder
                .root
                .join("docs/assets/typst/portal-graphs/layout.typ"),
        )
        .unwrap();
        let portal_edge_style = fs::read_to_string(
            builder
                .root
                .join("docs/assets/typst/portal-graphs/edge-style.typ"),
        )
        .unwrap();
        for structural_contract in [
            "#let edge-label-style(edge)",
            "#let autogen-external-edge-fields(",
            "#let layout(",
            "#let bind-layout(",
            "graph.parse(input)",
            "graph.style(",
            "for pass in passes {",
            "graph-bytes = apply-layout(",
            "_layout-pass(",
            "draw(",
        ] {
            assert!(
                shared_layout.contains(structural_contract),
                "shared GammaLoop layout lost {structural_contract}"
            );
        }
        assert!(canonical_layout.contains("#import \"layout-core.typ\": bind-layout"));
        assert!(canonical_layout.contains("#import \"physics-edge-style.typ\" as physics"));
        assert!(portal_layout.contains(
            "#import \"../../../../assets/embedded/drawing/templates/layout-core.typ\": ("
        ));
        assert!(portal_layout.contains(
            "#import \"../../../../assets/embedded/drawing/templates/physics-edge-style.typ\" as physics"
        ));
        assert!(portal_edge_style.contains(
            "#import \"../../../../assets/embedded/drawing/templates/physics-edge-style.typ\": ("
        ));
        assert!(portal_layout.contains("cetz.draw.bezier("));
        for (name, adapter) in [
            ("save-dot", canonical_layout.as_str()),
            ("website", portal_layout.as_str()),
        ] {
            assert!(
                adapter.contains("#let layout = bind-layout("),
                "{name} adapter does not bind the shared layout"
            );
            assert!(
                !adapter.contains("graph.parse(input)"),
                "{name} adapter duplicated the shared layout algorithm"
            );
        }
        assert!(portal_layout.contains("diagram-options: ("));
        assert!(portal_layout.contains("node-stroke: palette.ink + 1.45pt"));

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
        for product in &builder.registry.product {
            assert!(citations.contains(&format!(
                "<article class=\"citation-card\" id=\"{}\"",
                product.id
            )));
        }
        assert!(citations.contains("10.5281/zenodo.18429583"));
        assert!(citations.contains("No registered software DOI is currently configured"));
    }

    #[test]
    fn portal_task_chooser_routes_each_goal_to_the_registered_quickstart() {
        let builder = SiteBuilder::discover().unwrap();
        for (channel, tag, channel_route) in [
            (BuildChannel::Latest, None, "latest"),
            (BuildChannel::Snapshot, Some("v1.2.3"), "snapshots/v1.2.3"),
        ] {
            let output = tempfile::tempdir().unwrap();
            builder.write_portal(output.path(), channel, tag).unwrap();
            let html = fs::read_to_string(output.path().join("index.html")).unwrap();

            assert!(html.contains("id=\"tasks\" aria-labelledby=\"tasks-title\""));
            assert_eq!(html.matches("class=\"portal-task-link\"").count(), 5);
            for (task, product_id, role) in PORTAL_TASKS {
                let product = builder
                    .registry
                    .product
                    .iter()
                    .find(|product| product.id == product_id)
                    .unwrap();
                let quickstart = product
                    .pages
                    .iter()
                    .find(|page| page.id == "quickstart")
                    .unwrap();
                assert!(html.contains(&format!(
                    "class=\"portal-task-link\" data-product=\"{product_id}\" href=\"products/{product_id}/{channel_route}/{}\"",
                    quickstart.route
                )));
                assert!(html.contains(&format!("<strong>{task}</strong>")));
                assert!(html.contains(&format!("{} · {role}", product.title)));
            }
        }
    }

    #[test]
    fn portal_project_cards_use_accurate_ecosystem_roles() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        builder
            .write_portal(output.path(), BuildChannel::Latest, None)
            .unwrap();
        let html = fs::read_to_string(output.path().join("index.html")).unwrap();

        assert!(!html.contains("Research project"));
        for (_, product_id, role) in PORTAL_TASKS {
            assert!(html.contains(&format!(
                "<article class=\"portal-project-card\" data-product=\"{product_id}\"><div class=\"portal-project-meta\"><span>"
            )));
            assert!(html.contains(&format!("</span><span>{role}</span></div>")));
            assert!(html.contains(&format!(
                "products/{product_id}/latest/quickstart/\">Get started</a>"
            )));
        }
    }

    #[test]
    fn every_product_orients_rust_to_native_rustdoc_and_keeps_legacy_redirects() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        for product in &builder.registry.product {
            let site = output.path().join(&product.id);
            fs::create_dir_all(&site).unwrap();
            builder.write_site_assets(&site).unwrap();
            builder.write_component_catalogs(product, &site).unwrap();
            let mut generated_pages = builder.write_generated_reference(product, &site).unwrap();
            let metadata = builder
                .metadata(product, BuildChannel::Latest, None, Path::new("latest"))
                .unwrap();
            let python_pages = builder
                .write_python_reference(product, &metadata, &site)
                .unwrap();
            generated_pages.extend(python_pages.iter().cloned());
            builder
                .write_rust_reference_with_availability(product, &site, true)
                .unwrap();
            builder.write_reference_hub(product, &site).unwrap();

            let hub = fs::read_to_string(site.join("reference/index.html")).unwrap();
            assert!(hub.contains("Rust API"), "{}", product.id);
            assert!(hub.contains("Python API"), "{}", product.id);
            if product
                .pages
                .iter()
                .any(|page| page.route == "reference/typst/")
            {
                assert!(hub.contains("href=\"reference/typst/\">Typst API</a>"));
            }
            let rust_index = fs::read_to_string(site.join("reference/rust/index.html")).unwrap();
            assert!(
                rust_index.contains("Rustdoc is the canonical Rust API reference"),
                "{}",
                product.id
            );
            assert!(
                rust_index.contains("href=\"quickstart/rust/\""),
                "{}",
                product.id
            );
            assert!(!rust_index.contains("class=\"api-item\""), "{}", product.id);
            for component in &product.python_components {
                let module_route = format!("reference/python/{}/", component.id);
                let component_pages = python_pages
                    .iter()
                    .filter(|page| page.route.starts_with(&module_route))
                    .collect::<Vec<_>>();
                let module_page = SitePage::new(
                    &module_route,
                    format!("{} Python API", python_display_module(component)),
                    "Reference",
                );
                let sidebar = builder.site_sidebar(
                    product,
                    &metadata,
                    &module_page,
                    "",
                    BuildScope::FullSite,
                    &generated_pages,
                );
                assert_eq!(
                    sidebar.matches("sidebar-python-symbol").count(),
                    component_pages.len(),
                    "{}:{}",
                    product.id,
                    component.id,
                );
                for page in component_pages {
                    assert!(
                        sidebar.contains(&format!("href=\"{}\"", page.route)),
                        "{}:{}:{}",
                        product.id,
                        component.id,
                        page.route,
                    );
                    let symbol =
                        fs::read_to_string(site.join(&page.route).join("index.html")).unwrap();
                    let short_title = page.title.rsplit('.').next().unwrap();
                    assert!(
                        symbol.contains(&format!("<h1>{}</h1>", escape_html(short_title))),
                        "{}:{}:{}",
                        product.id,
                        component.id,
                        page.route,
                    );
                    if symbol.contains("class=\"api-signature\"") {
                        assert!(symbol.contains("<code data-lang=\"python\">"));
                    }
                }
            }
            for component in &product.rust_components {
                let legacy = site
                    .join("reference/rust/supported")
                    .join(&component.id)
                    .join("index.html");
                let crate_name = component.package.replace('-', "_");
                assert!(legacy.is_file(), "{}:{}", product.id, component.id);
                assert!(rust_index.contains(&format!(
                    "href=\"reference/rust/{crate_name}/\">Open canonical Rustdoc"
                )));
                let redirect = fs::read_to_string(legacy).unwrap();
                assert!(redirect.contains("name=\"robots\" content=\"noindex\""));
                assert!(redirect.contains(&format!(
                    "http-equiv=\"refresh\" content=\"0; url=../../{crate_name}/\""
                )));
                assert!(
                    redirect.contains(&format!("rel=\"canonical\" href=\"../../{crate_name}/\""))
                );
                assert!(!redirect.contains("assets/site.css"));
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
        for product in &builder.registry.product {
            let hero = product_hero(product);
            assert!(hero.contains("href=\"quickstart/\""), "{}", product.id);
            assert!(hero.contains("Get started"), "{}", product.id);
        }
    }

    #[test]
    fn generated_search_uses_cli_leaf_routes_and_setting_namespaces() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "gammaloop")
            .unwrap();
        let site = tempfile::tempdir().unwrap();
        builder
            .write_component_catalogs(product, site.path())
            .unwrap();
        builder
            .write_search_index(product, site.path(), false)
            .unwrap();
        let entries = serde_json::from_slice::<Vec<SearchEntry>>(
            &fs::read(site.path().join("search-index.json")).unwrap(),
        )
        .unwrap();
        let (reference, _) = builder.generated_references().unwrap();
        let command = visible_cli_commands(&reference).into_iter().next().unwrap();
        let setting = reference.settings.first().unwrap();
        let namespace = cli_setting_namespace(&setting.path);
        let setting_route = format!(
            "{}#{}",
            cli_setting_namespace_route(namespace),
            generated_anchor("setting", &setting.path)
        );

        assert!(entries.iter().any(|entry| {
            entry.kind == "command"
                && entry.title == command.path
                && entry.href == cli_command_route(&command.path)
        }));
        assert!(entries.iter().any(|entry| {
            entry.kind == "setting" && entry.title == setting.path && entry.href == setting_route
        }));
        assert!(entries.iter().any(|entry| {
            entry.kind == "settings namespace"
                && entry.title == format!("{namespace} settings")
                && entry.href == cli_setting_namespace_route(namespace)
        }));
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
        builder
            .write_developer_docs(output.path(), false, None)
            .unwrap();

        let developer_root = output.path().join("developers");
        assert!(!stale.exists());
        let hub = fs::read_to_string(developer_root.join("index.html")).unwrap();
        assert!(hub.contains("For developers"));
        assert!(hub.contains("href=\"assets/local-unitarity-light.svg\""));
        assert!(hub.contains("Implemented architecture"));
        assert!(hub.contains("Design proposals"));
        assert!(hub.contains("Investigation record"));
        assert!(hub.contains("Performance investigation"));
        let expected_notes = builder
            .developers
            .section
            .iter()
            .map(|section| section.note.len())
            .sum::<usize>();
        assert_eq!(
            hub.matches("class=\"developer-card\"").count(),
            expected_notes
        );
        assert!(hub.contains("GammaLoop drawing style migration"));

        let current = fs::read_to_string(
            developer_root.join("architecture/gammaloop-architecture/index.html"),
        )
        .unwrap();
        assert!(current.contains("developer-status-implemented"));
        assert!(current.contains("View the source note"));
        assert!(current.contains("architecture-current.typ"));
        assert!(current.contains("Typst rendering was disabled"));
        assert!(current.contains("aria-label=\"Developer note pagination\""));

        let investigation = fs::read_to_string(
            developer_root.join("architecture/spenso-improvement-effects/index.html"),
        )
        .unwrap();
        assert!(investigation.contains("developer-status-historical-experiment"));

        let parsing =
            fs::read_to_string(developer_root.join("architecture/schoonschip-network/index.html"))
                .unwrap();
        assert!(parsing.contains("For developers"));
        assert!(parsing.contains("schoonschip-net-parsing.typ"));

        let search = serde_json::from_slice::<Vec<SearchEntry>>(
            &fs::read(developer_root.join("search-index.json")).unwrap(),
        )
        .unwrap();
        assert!(search.iter().any(|entry| {
            entry.href == "architecture/gammaloop-architecture/"
                && entry.kind == "developer current implemented"
        }));
        assert!(search.iter().any(|entry| {
            entry.href == "architecture/documentation-improvement-plan/"
                && entry.kind == "developer proposal in progress"
        }));

        let heading = HeadingLink {
            level: 2,
            title: "Scope".to_owned(),
            id: "scope".to_owned(),
        };
        let current_note = builder
            .developers
            .section
            .iter()
            .flat_map(|section| &section.note)
            .find(|note| note.lifecycle == "current")
            .unwrap();
        let archived_note = builder
            .developers
            .section
            .iter()
            .flat_map(|section| &section.note)
            .find(|note| note.lifecycle == "archived")
            .unwrap();
        let current_heading =
            developer_heading_search_entry(current_note, "architecture/current/", heading.clone());
        let archived_heading =
            developer_heading_search_entry(archived_note, "architecture/archived/", heading);
        assert!(current_heading.kind.starts_with("developer current "));
        assert!(current_heading.kind.ends_with(" heading"));
        assert!(archived_heading.kind.starts_with("developer archived "));
        assert!(archived_heading.kind.ends_with(" heading"));
    }

    #[test]
    fn federated_search_covers_every_product_and_developer_notes() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        for product in &builder.registry.product {
            let product_root = output
                .path()
                .join("products")
                .join(&product.id)
                .join("latest");
            fs::create_dir_all(&product_root).unwrap();
            fs::write(
                product_root.join("search-index.json"),
                serde_json::to_vec(&vec![SearchEntry {
                    title: product.title.clone(),
                    summary: product.tagline.clone(),
                    href: "tutorial/#first-result".to_owned(),
                    kind: "tutorial".to_owned(),
                    text: String::new(),
                }])
                .unwrap(),
            )
            .unwrap();
        }
        let developer_root = output.path().join("developers");
        fs::create_dir_all(&developer_root).unwrap();
        fs::write(
            developer_root.join("search-index.json"),
            serde_json::to_vec(&vec![SearchEntry {
                title: "Current architecture".to_owned(),
                summary: "Implementation boundaries".to_owned(),
                href: "architecture/current/".to_owned(),
                kind: "developer current implemented".to_owned(),
                text: String::new(),
            }])
            .unwrap(),
        )
        .unwrap();

        builder
            .write_federated_search_index(output.path(), BuildChannel::Latest, None)
            .unwrap();
        let search = serde_json::from_slice::<Vec<SearchEntry>>(
            &fs::read(output.path().join("search-index.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(
            search.len(),
            PRODUCT_IDS.len() * 2 + builder.talks.talk.len() + 3
        );
        for product in &builder.registry.product {
            assert!(search.iter().any(|entry| {
                entry.href == format!("products/{}/latest/tutorial/#first-result", product.id)
                    && entry.kind == format!("{} · tutorial", product.title)
            }));
            let citation = search
                .iter()
                .find(|entry| entry.title == format!("Cite {}", product.title))
                .unwrap();
            assert_eq!(citation.href, format!("citations/#{}", product.id));
            assert_eq!(citation.kind, "citation");
            assert!(citation.summary.contains("software citation and BibTeX"));
            assert!(citation.text.contains(&product.citation.title));
        }
        assert!(search.iter().any(|entry| {
            entry.href == "developers/architecture/current/"
                && entry.kind == "Developers · developer current implemented"
        }));
    }

    #[test]
    fn rendered_product_and_reference_pages_offer_revision_pinned_feedback() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "gammaloop")
            .unwrap();
        let site = tempfile::tempdir().unwrap();
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Latest,
            snapshot_tag: None,
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/gammaloop/latest/".to_owned(),
            components: Vec::new(),
        };
        for (page, expected_source, expected_search, expected_citation) in [
            (
                SitePage::new("tutorial/", "Create your first state", "Tutorial"),
                "docs/products/gammaloop/content/tutorial.typ",
                "../../../../search-index.json",
                "../../../../citations/#gammaloop",
            ),
            (
                SitePage::new("reference/cli/", "CLI commands and settings", "Reference"),
                "crates/gammaloop-api/src/commands/mod.rs",
                "../../../../../search-index.json",
                "../../../../../citations/#gammaloop",
            ),
        ] {
            let destination = site.path().join(&page.route);
            fs::create_dir_all(&destination).unwrap();
            fs::write(
                destination.join("index.html"),
                "<!doctype html><html><body><h2>Page title</h2></body></html>",
            )
            .unwrap();
            builder
                .decorate_html_page(
                    product,
                    &metadata,
                    site.path(),
                    &page,
                    (None, None),
                    (BuildScope::FullSite, &[]),
                )
                .unwrap();
            let rendered = fs::read_to_string(destination.join("index.html")).unwrap();
            assert!(rendered.contains(&format!("data-search-index=\"{expected_search}\"")));
            assert!(rendered.contains(&format!("href=\"{expected_citation}\">Cite GammaLoop</a>")));
            assert!(rendered.contains(&format!(
                "https://github.com/alphal00p/gammaloop/blob/0123456789abcdef/{expected_source}"
            )));
            assert!(rendered.contains("issues/new?labels=documentation&amp;title="));
            assert!(rendered.contains("Report a documentation issue"));
        }
    }

    #[test]
    fn dense_reference_tocs_only_list_sections_and_omit_empty_pagination() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "gammaloop")
            .unwrap();
        let component = product.python_components.first().unwrap();
        let page = SitePage::new(
            format!("reference/python/{}/Engine/", component.id),
            "Engine",
            "Python API",
        );
        let site = tempfile::tempdir().unwrap();
        let destination = site.path().join(&page.route);
        fs::create_dir_all(&destination).unwrap();
        fs::write(
            destination.join("index.html"),
            "<!doctype html><html><body><h1>Engine</h1><h2>Members</h2><h3>run</h3></body></html>",
        )
        .unwrap();
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Latest,
            snapshot_tag: None,
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/gammaloop/latest/".to_owned(),
            components: Vec::new(),
        };
        builder
            .decorate_html_page(
                product,
                &metadata,
                site.path(),
                &page,
                (None, None),
                (BuildScope::FullSite, &[]),
            )
            .unwrap();
        let rendered = fs::read_to_string(destination.join("index.html")).unwrap();

        assert!(
            rendered
                .contains("<a class=\"toc-link\" data-level=\"2\" href=\"#members\">Members</a>")
        );
        assert!(!rendered.contains("<a class=\"toc-link\" data-level=\"3\" href=\"#run\">run</a>"));
        assert!(rendered.contains("<h3 id=\"run\">run</h3>"));
        assert!(!rendered.contains("aria-label=\"Documentation pagination\""));
    }

    #[test]
    fn complex_reference_pages_keep_an_inline_outline_below_the_right_rail() {
        let builder = SiteBuilder::discover().unwrap();
        let product = builder
            .registry
            .product
            .iter()
            .find(|product| product.id == "linnet")
            .unwrap();
        let page = SitePage::new("reference/typst/graph/", "graph module", "Typst API");
        let site = tempfile::tempdir().unwrap();
        let destination = site.path().join(&page.route);
        fs::create_dir_all(&destination).unwrap();
        fs::write(
            destination.join("index.html"),
            "<!doctype html><html><body><h1>graph module</h1><h2>Concepts</h2><h2>Reference</h2><h3>build</h3></body></html>",
        )
        .unwrap();
        let metadata = SnapshotMetadata {
            schema: SCHEMA_VERSION,
            product: &product.id,
            title: &product.title,
            channel: BuildChannel::Latest,
            snapshot_tag: None,
            git_commit: "0123456789abcdef".to_owned(),
            git_timestamp: 1,
            route: "products/linnet/latest/".to_owned(),
            components: Vec::new(),
        };

        builder
            .decorate_html_page(
                product,
                &metadata,
                site.path(),
                &page,
                (None, None),
                (BuildScope::FullSite, &[]),
            )
            .unwrap();
        let rendered = fs::read_to_string(destination.join("index.html")).unwrap();

        assert!(rendered.contains("<aside class=\"docs-toc\" aria-label=\"On this page\">"));
        assert!(rendered.contains("<details class=\"docs-inline-toc\">"));
    }

    #[test]
    fn moved_manual_routes_redirect_without_entering_navigation() {
        let builder = SiteBuilder::discover().unwrap();
        let site = tempfile::tempdir().unwrap();
        let mut alias_count = 0;

        for product in &builder.registry.product {
            let product_site = site.path().join(&product.id);
            for page in &product.pages {
                let directory = product_site.join(&page.route);
                fs::create_dir_all(&directory).unwrap();
                fs::write(
                    directory.join("index.html"),
                    "<!doctype html><html><body><h1>Page</h1></body></html>",
                )
                .unwrap();
            }
            let metadata = SnapshotMetadata {
                schema: SCHEMA_VERSION,
                product: &product.id,
                title: &product.title,
                channel: BuildChannel::Latest,
                snapshot_tag: None,
                git_commit: "0123456789abcdef".to_owned(),
                git_timestamp: 1,
                route: format!("products/{}/latest/", product.id),
                components: Vec::new(),
            };
            builder
                .decorate_site_pages_with_generated(
                    product,
                    &metadata,
                    &product_site,
                    BuildScope::FullSite,
                    &[],
                )
                .unwrap();

            let navigation_routes = product_site_pages(product)
                .into_iter()
                .map(|page| page.route)
                .collect::<BTreeSet<_>>();
            for page in &product.pages {
                for alias in &page.aliases {
                    alias_count += 1;
                    assert!(alias.starts_with("manual/"), "{}:{alias}", product.id);
                    assert!(!navigation_routes.contains(alias));
                    let target = format!("{}{}", page_root_prefix(alias), page.route);
                    let redirect =
                        fs::read_to_string(product_site.join(alias).join("index.html")).unwrap();
                    assert!(redirect.contains(&format!(
                        "http-equiv=\"refresh\" content=\"0; url={target}\""
                    )));
                    assert!(redirect.contains(&format!("rel=\"canonical\" href=\"{target}\"")));
                    assert!(!redirect.contains("assets/site.css"));
                }
            }
        }

        assert_eq!(alias_count, 23);
    }

    #[test]
    fn portal_and_developer_pages_use_federated_search_and_feedback() {
        let builder = SiteBuilder::discover().unwrap();
        let output = tempfile::tempdir().unwrap();
        builder
            .write_portal(output.path(), BuildChannel::Latest, None)
            .unwrap();
        let portal = fs::read_to_string(output.path().join("index.html")).unwrap();
        assert!(portal.contains("data-search-index=\"search-index.json\""));
        assert!(portal.contains("class=\"portal-search-button\""));
        assert!(portal.contains("Search all projects and developer notes"));

        builder
            .write_developer_docs(output.path(), false, None)
            .unwrap();
        let hub = fs::read_to_string(output.path().join("developers/index.html")).unwrap();
        assert!(hub.contains("data-search-index=\"../search-index.json\""));
        assert!(hub.contains("Report a documentation issue"));
        assert!(hub.contains("blob/"));
        assert!(hub.contains("/docs/developers.toml"));

        let note = fs::read_to_string(
            output
                .path()
                .join("developers/architecture/gammaloop-architecture/index.html"),
        )
        .unwrap();
        assert!(note.contains("data-search-index=\"../../../search-index.json\""));
        assert!(note.contains("Report a documentation issue"));
        assert!(note.contains("architecture-current.typ"));
        assert!(note.contains("issues/new?labels=documentation&amp;title="));
    }
}
