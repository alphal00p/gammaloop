//! Renderer-independent metadata shared by the alphal00p documentation tools.
//!
//! Product crates do not depend on this crate. Dedicated documentation catalog
//! exporters build [`DocCatalog`] values and serialize them for the site builder.

pub mod generated;

use std::fmt;

use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Version of the serialized catalog schema.
pub const SCHEMA_VERSION: u32 = 3;

/// How a documentation body must be interpreted by a renderer.
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum DocFormat {
    /// Markdown as understood by Rustdoc.
    #[default]
    RustMarkdown,
    /// Trusted Typst markup from an explicitly annotated source.
    TypstMarkup,
    /// A Python docstring, normally using NumPy-style sections.
    PythonDocstring,
    /// Text that must not be interpreted as markup.
    PlainText,
}

/// Public API surface represented by a catalog.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ApiLanguage {
    Rust,
    Python,
    CommandLine,
    Configuration,
    Content,
}

/// Kind of a top-level documented item.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum DocItemKind {
    Function,
    Type,
    Trait,
    ExportedMacro,
    Method,
    Field,
    Variant,
    Command,
    Setting,
    PythonClass,
    PythonFunction,
    ContentPage,
}

/// Stable identity of a separately published documentation site.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocProduct {
    pub id: String,
    pub title: String,
}

impl DocProduct {
    pub fn new(id: impl Into<String>, title: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            title: title.into(),
        }
    }
}

/// Identity and version of one Rust, Python, CLI, or configuration component.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocComponent {
    pub id: String,
    pub package: String,
    pub title: String,
    pub version: String,
    pub language: ApiLanguage,
    /// Feature matrix used to compile this component's exhaustive reference.
    #[serde(default)]
    pub features: Vec<String>,
}

impl DocComponent {
    pub fn new(
        id: impl Into<String>,
        package: impl Into<String>,
        title: impl Into<String>,
        version: impl Into<String>,
        language: ApiLanguage,
    ) -> Self {
        Self {
            id: id.into(),
            package: package.into(),
            title: title.into(),
            version: version.into(),
            language,
            features: vec![],
        }
    }
}

/// Kind of an item member.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum DocMemberKind {
    Parameter,
    /// One signature and its documentation within an overloaded callable.
    Overload,
    Field,
    Variant,
    AssociatedFunction,
    AssociatedType,
    AssociatedConst,
    Method,
    Getter,
    Setter,
}

/// Documentation text together with the format that owns its escaping rules.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocText {
    pub format: DocFormat,
    pub body: String,
}

impl DocText {
    pub fn new(format: DocFormat, body: impl Into<String>) -> Self {
        Self {
            format,
            body: body.into(),
        }
    }
}

/// Source position of a generated descriptor.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SourceLocation {
    /// Stable public identifier, independent of line movement.
    pub identifier: String,
    pub file: String,
    pub line: u32,
    pub column: u32,
}

impl SourceLocation {
    pub fn new(
        identifier: impl Into<String>,
        file: impl Into<String>,
        line: u32,
        column: u32,
    ) -> Self {
        Self {
            identifier: identifier.into(),
            file: file.into(),
            line,
            column,
        }
    }
}

/// A named function or method parameter.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocParam {
    pub name: String,
    pub docs: Option<DocText>,
    pub ty: Option<String>,
    pub default: Option<String>,
    pub required: bool,
}

impl DocParam {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            docs: None,
            ty: None,
            default: None,
            required: true,
        }
    }
}

/// A code example attached to one supported API entry.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocExample {
    pub title: String,
    pub language: String,
    pub code: String,
}

impl DocExample {
    pub fn new(
        title: impl Into<String>,
        language: impl Into<String>,
        code: impl Into<String>,
    ) -> Self {
        Self {
            title: title.into(),
            language: language.into(),
            code: code.into(),
        }
    }
}

/// Parameter, overload, field, variant, or associated item nested under a documented item.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocMember {
    pub name: String,
    pub kind: DocMemberKind,
    pub docs: Option<DocText>,
    pub signature: Option<String>,
    pub default: Option<String>,
    pub members: Vec<DocMember>,
}

impl DocMember {
    pub fn new(name: impl Into<String>, kind: DocMemberKind) -> Self {
        Self {
            name: name.into(),
            kind,
            docs: None,
            signature: None,
            default: None,
            members: vec![],
        }
    }
}

/// A function, type, trait, or macro descriptor.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocItem {
    /// Stable identifier within the containing scope.
    pub id: String,
    /// Source-level or public API name.
    pub name: String,
    /// Human-readable title.
    pub title: String,
    pub kind: DocItemKind,
    /// Whether this item belongs to the product's intentionally supported API.
    ///
    /// Generated exhaustive sidecars may also contain reachable exports whose
    /// compatibility is not promised by the curated support catalog.
    pub supported: bool,
    /// Cargo/Python features that must be enabled for this item to exist.
    #[serde(default)]
    pub required_features: Vec<String>,
    pub docs: Option<DocText>,
    pub summary: Option<String>,
    pub signature: Option<String>,
    pub since: Option<String>,
    pub keywords: Vec<String>,
    pub params: Vec<DocParam>,
    pub returns: Option<DocText>,
    pub examples: Vec<DocExample>,
    pub members: Vec<DocMember>,
    pub source: Option<SourceLocation>,
}

impl DocItem {
    pub fn new(
        id: impl Into<String>,
        name: impl Into<String>,
        title: impl Into<String>,
        kind: DocItemKind,
    ) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            title: title.into(),
            kind,
            supported: true,
            required_features: vec![],
            docs: None,
            summary: None,
            signature: None,
            since: None,
            keywords: vec![],
            params: vec![],
            returns: None,
            examples: vec![],
            members: vec![],
            source: None,
        }
    }
}

/// An explicitly ordered group of items and child scopes.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocScope {
    /// Stable identifier within the parent scope.
    pub id: String,
    pub title: String,
    pub docs: Option<DocText>,
    pub summary: Option<String>,
    pub since: Option<String>,
    pub keywords: Vec<String>,
    pub source: Option<SourceLocation>,
    pub items: IndexMap<String, DocItem>,
    pub scopes: IndexMap<String, DocScope>,
}

impl DocScope {
    pub fn new(id: impl Into<String>, title: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            title: title.into(),
            docs: None,
            summary: None,
            since: None,
            keywords: vec![],
            source: None,
            items: IndexMap::new(),
            scopes: IndexMap::new(),
        }
    }

    /// Insert an item at the end of this scope.
    pub fn define_item(&mut self, item: DocItem) -> Result<&mut DocItem, SchemaError> {
        validate_identifier("item", &item.id)?;
        if self.items.contains_key(&item.id) || self.scopes.contains_key(&item.id) {
            return Err(SchemaError::DuplicateEntry {
                scope: self.id.clone(),
                id: item.id,
            });
        }
        let id = item.id.clone();
        self.items.insert(id.clone(), item);
        Ok(self.items.get_mut(&id).expect("item was just inserted"))
    }

    /// Insert a child scope at the end of this scope.
    pub fn define_scope(&mut self, scope: DocScope) -> Result<&mut DocScope, SchemaError> {
        validate_identifier("scope", &scope.id)?;
        if self.items.contains_key(&scope.id) || self.scopes.contains_key(&scope.id) {
            return Err(SchemaError::DuplicateEntry {
                scope: self.id.clone(),
                id: scope.id,
            });
        }
        let id = scope.id.clone();
        self.scopes.insert(id.clone(), scope);
        Ok(self.scopes.get_mut(&id).expect("scope was just inserted"))
    }

    /// Validate this scope and all descendants after construction or deserialization.
    pub fn validate(&self) -> Result<(), SchemaError> {
        self.validate_at(&self.id)
    }

    fn validate_at(&self, path: &str) -> Result<(), SchemaError> {
        validate_identifier("scope", &self.id)?;
        validate_non_empty("scope title", &self.title)?;

        for (key, item) in &self.items {
            if key != &item.id {
                return Err(SchemaError::MismatchedMapKey {
                    path: path.to_owned(),
                    key: key.clone(),
                    id: item.id.clone(),
                });
            }
            if self.scopes.contains_key(key) {
                return Err(SchemaError::DuplicateEntry {
                    scope: path.to_owned(),
                    id: key.clone(),
                });
            }
            validate_identifier("item", &item.id)?;
            validate_non_empty("item name", &item.name)?;
            validate_non_empty("item title", &item.title)?;
        }

        for (key, scope) in &self.scopes {
            if key != &scope.id {
                return Err(SchemaError::MismatchedMapKey {
                    path: path.to_owned(),
                    key: key.clone(),
                    id: scope.id.clone(),
                });
            }
            let child_path = format!("{path}::{key}");
            scope.validate_at(&child_path)?;
        }

        Ok(())
    }
}

/// Top-level serialization unit emitted by a component-specific exporter.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct DocCatalog {
    pub schema_version: u32,
    pub product: DocProduct,
    pub component: DocComponent,
    pub root: DocScope,
}

impl DocCatalog {
    pub fn new(product: DocProduct, component: DocComponent, root: DocScope) -> Self {
        Self {
            schema_version: SCHEMA_VERSION,
            product,
            component,
            root,
        }
    }

    pub fn validate(&self) -> Result<(), SchemaError> {
        if self.schema_version != SCHEMA_VERSION {
            return Err(SchemaError::UnsupportedSchemaVersion {
                expected: SCHEMA_VERSION,
                found: self.schema_version,
            });
        }
        validate_identifier("product", &self.product.id)?;
        validate_non_empty("product title", &self.product.title)?;
        validate_identifier("component", &self.component.id)?;
        validate_non_empty("component package", &self.component.package)?;
        validate_non_empty("component title", &self.component.title)?;
        validate_non_empty("component version", &self.component.version)?;
        self.root.validate()
    }

    /// Enforce the stronger contract used for the explicitly supported catalog.
    pub fn validate_supported(&self) -> Result<(), SchemaError> {
        self.validate()?;
        let mut item_count = 0;
        validate_supported_scope(&self.root, &mut item_count)?;
        if item_count == 0 {
            return Err(SchemaError::EmptySupportedCatalog {
                component: self.component.id.clone(),
            });
        }
        Ok(())
    }
}

fn validate_supported_scope(scope: &DocScope, item_count: &mut usize) -> Result<(), SchemaError> {
    for item in scope.items.values() {
        if !item.supported {
            continue;
        }
        *item_count += 1;
        let field = |name| SchemaError::IncompleteSupportedItem {
            item: item.id.clone(),
            field: name,
        };
        if item
            .docs
            .as_ref()
            .is_none_or(|docs| docs.body.trim().is_empty())
        {
            return Err(field("docs"));
        }
        if item
            .summary
            .as_ref()
            .is_none_or(|summary| summary.trim().is_empty())
        {
            return Err(field("summary"));
        }
        if item
            .signature
            .as_ref()
            .is_none_or(|signature| signature.trim().is_empty())
        {
            return Err(field("signature"));
        }
        if item.source.as_ref().is_none_or(|source| {
            source.identifier.trim().is_empty() || source.file.trim().is_empty()
        }) {
            return Err(field("source"));
        }
        if item.examples.is_empty()
            || item
                .examples
                .iter()
                .any(|example| example.code.trim().is_empty())
        {
            return Err(field("examples"));
        }
    }
    for child in scope.scopes.values() {
        validate_supported_scope(child, item_count)?;
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, Error, PartialEq)]
pub enum SchemaError {
    #[error("{kind} identifier must not be empty")]
    EmptyIdentifier { kind: &'static str },
    #[error("invalid {kind} identifier `{id}`")]
    InvalidIdentifier { kind: &'static str, id: String },
    #[error("{field} must not be empty")]
    EmptyField { field: &'static str },
    #[error("duplicate documentation entry `{id}` in scope `{scope}`")]
    DuplicateEntry { scope: String, id: String },
    #[error("map key `{key}` does not match entry id `{id}` at `{path}`")]
    MismatchedMapKey {
        path: String,
        key: String,
        id: String,
    },
    #[error("unsupported documentation schema version {found}; expected {expected}")]
    UnsupportedSchemaVersion { expected: u32, found: u32 },
    #[error("supported catalog for `{component}` has no entries")]
    EmptySupportedCatalog { component: String },
    #[error("supported item `{item}` is missing {field}")]
    IncompleteSupportedItem { item: String, field: &'static str },
}

fn validate_identifier(kind: &'static str, id: &str) -> Result<(), SchemaError> {
    if id.is_empty() {
        return Err(SchemaError::EmptyIdentifier { kind });
    }
    if id
        .chars()
        .any(|character| !(character.is_alphanumeric() || "_-.:".contains(character)))
    {
        return Err(SchemaError::InvalidIdentifier {
            kind,
            id: id.to_owned(),
        });
    }
    Ok(())
}

fn validate_non_empty(field: &'static str, value: &str) -> Result<(), SchemaError> {
    if value.trim().is_empty() {
        Err(SchemaError::EmptyField { field })
    } else {
        Ok(())
    }
}

impl fmt::Display for ApiLanguage {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let value = match self {
            Self::Rust => "rust",
            Self::Python => "python",
            Self::CommandLine => "command-line",
            Self::Configuration => "configuration",
            Self::Content => "content",
        };
        formatter.write_str(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn item(id: &str) -> DocItem {
        DocItem::new(id, id, id, DocItemKind::Function)
    }

    #[test]
    fn scopes_preserve_definition_order_through_json() {
        let mut root = DocScope::new("root", "Root");
        root.define_item(item("second")).unwrap();
        root.define_item(item("first")).unwrap();
        root.define_scope(DocScope::new("zeta", "Zeta")).unwrap();
        root.define_scope(DocScope::new("alpha", "Alpha")).unwrap();

        let json = serde_json::to_string(&root).unwrap();
        assert!(json.find("second").unwrap() < json.find("first").unwrap());
        assert!(json.find("zeta").unwrap() < json.find("alpha").unwrap());
        assert_eq!(serde_json::from_str::<DocScope>(&json).unwrap(), root);
    }

    #[test]
    fn an_item_and_scope_cannot_share_a_route_segment() {
        let mut root = DocScope::new("root", "Root");
        root.define_item(item("network")).unwrap();

        assert_eq!(
            root.define_scope(DocScope::new("network", "Network"))
                .unwrap_err(),
            SchemaError::DuplicateEntry {
                scope: "root".to_owned(),
                id: "network".to_owned(),
            }
        );
    }

    #[test]
    fn validation_rejects_tampered_map_keys() {
        let mut root = DocScope::new("root", "Root");
        root.items.insert("alias".to_owned(), item("canonical"));

        assert!(matches!(
            root.validate(),
            Err(SchemaError::MismatchedMapKey { key, id, .. })
                if key == "alias" && id == "canonical"
        ));
    }

    #[test]
    fn catalogs_validate_schema_and_identity() {
        let product = DocProduct::new("spenso", "Spenso");
        let component = DocComponent::new(
            "spenso-py",
            "spynso3",
            "Spenso Python API",
            "0.1.2",
            ApiLanguage::Python,
        );
        let mut catalog = DocCatalog::new(product, component, DocScope::new("spenso", "Spenso"));
        catalog.validate().unwrap();

        catalog.schema_version += 1;
        assert_eq!(
            catalog.validate().unwrap_err(),
            SchemaError::UnsupportedSchemaVersion {
                expected: SCHEMA_VERSION,
                found: SCHEMA_VERSION + 1,
            }
        );
    }
}
