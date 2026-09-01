//! Portable Spenso symbol-metadata attachments.
//!
//! Other Symbolica consumers treat these records as opaque data. Spenso-aware
//! runtimes validate the declarations, reconstruct the corresponding local
//! symbols, and only then import native Atom data. Representation metadata uses
//! `spenso.representation`; structured math displays use
//! `spenso.math-display`.

use std::{collections::BTreeMap, fmt, io::Cursor, str};

use ciborium::value::Value;
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{
        IndexDisplay, IndexPalette, IndexRow, LibraryRep, RepName, RepresentationClass,
        RepresentationMetadata, initialize as initialize_representations,
    },
};
use symbolica::atom::{Atom, AtomView, NamespacedSymbol, Symbol, SymbolBuilder};
use tymbolica_atom_payload::{Attachment, AttachmentKey, AttachmentSet};

/// Attachment schema used for portable Spenso representation declarations.
pub const REPRESENTATION_ATTACHMENT_SCHEMA: &str = "spenso.representation";
/// Current version of [`REPRESENTATION_ATTACHMENT_SCHEMA`].
pub const REPRESENTATION_ATTACHMENT_VERSION: u32 = 1;

/// Attachment schema used for portable structured math displays.
pub const MATH_DISPLAY_ATTACHMENT_SCHEMA: &str = "spenso.math-display";
/// Current version of [`MATH_DISPLAY_ATTACHMENT_SCHEMA`].
pub const MATH_DISPLAY_ATTACHMENT_VERSION: u32 = 1;

const MATH_DISPLAY_SYMBOL_PREFIX: &str = "spenso_index_";
const MATH_DISPLAY_HASH_DOMAIN: &[u8] = b"spenso-display-index-v1\0";

const MAX_DISPLAY_INDEX_DEPTH: usize = 16;
const MAX_DISPLAY_INDEX_NODES: usize = 64;

/// Validate and register all portable Spenso declarations in `attachments`.
///
/// Call this after payload/header inspection and before `Atom::import`. The
/// function validates every declaration and every existing-symbol collision
/// before it registers the first custom representation.
pub fn register_representation_attachments(
    attachments: &AttachmentSet,
) -> std::result::Result<(), String> {
    RepresentationDeclarations::from_attachment_set(attachments)
        .and_then(|declarations| declarations.register_before_atom_import())
        .map_err(|error| error.to_string())
}

/// Validate and register all portable structured math displays in `attachments`.
///
/// Call this before `Atom::import`, for the same reason as
/// [`register_representation_attachments`]: importing first would intern the
/// referenced display symbols without their local Spenso metadata.
pub fn register_math_display_attachments(
    attachments: &AttachmentSet,
) -> std::result::Result<(), String> {
    MathDisplayDeclarations::from_attachment_set(attachments)
        .and_then(|declarations| declarations.register_before_atom_import())
        .map_err(|error| error.to_string())
}

/// Collect all portable Spenso declarations referenced by `atom`.
///
/// The returned set contains both representation metadata and structured math
/// displays. It can be merged with application-specific attachments before the
/// Atom is exported through `tymbolica-atom-payload`.
pub fn attachments_for_atom(atom: &Atom) -> Result<AttachmentSet> {
    let representations = RepresentationDeclarations::referenced_by_atom(atom)?;
    let math_displays = MathDisplayDeclarations::referenced_by_atom(atom)?;
    let mut attachments = AttachmentSet::new();
    representations.append_attachments_to(&mut attachments)?;
    math_displays.append_attachments_to(&mut attachments)?;
    Ok(attachments)
}

/// Failure while inspecting or registering portable Spenso symbol metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepresentationAttachmentError {
    message: String,
}

impl RepresentationAttachmentError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for RepresentationAttachmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for RepresentationAttachmentError {}

/// Error while inspecting or registering a portable structured math display.
///
/// Both registries intentionally share one concrete error type so callers can
/// merge their pre-import validation without erasing error details.
pub type MathDisplayAttachmentError = RepresentationAttachmentError;

type Result<T> = std::result::Result<T, RepresentationAttachmentError>;

/// Stable representation classes carried by the attachment wire format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PortableRepresentationClass {
    SelfDual,
    Dualizable,
    InlineMetric,
}

impl PortableRepresentationClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SelfDual => "self-dual",
            Self::Dualizable => "dualizable",
            Self::InlineMetric => "inline-metric",
        }
    }

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "self-dual" => Ok(Self::SelfDual),
            "dualizable" => Ok(Self::Dualizable),
            "inline-metric" => Ok(Self::InlineMetric),
            other => Err(RepresentationAttachmentError::new(format!(
                "unknown representation class {other:?}"
            ))),
        }
    }
}

impl From<RepresentationClass> for PortableRepresentationClass {
    fn from(value: RepresentationClass) -> Self {
        match value {
            RepresentationClass::SelfDual => Self::SelfDual,
            RepresentationClass::Dualizable => Self::Dualizable,
            RepresentationClass::InlineMetric => Self::InlineMetric,
        }
    }
}

impl From<PortableRepresentationClass> for RepresentationClass {
    fn from(value: PortableRepresentationClass) -> Self {
        match value {
            PortableRepresentationClass::SelfDual => Self::SelfDual,
            PortableRepresentationClass::Dualizable => Self::Dualizable,
            PortableRepresentationClass::InlineMetric => Self::InlineMetric,
        }
    }
}

/// Declarative, runtime-independent description of a Spenso representation.
///
/// The representation's fully-qualified name lives in the attachment key, not
/// in this value.  Its display label is always derived from that name.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RepresentationDeclaration {
    pub class: PortableRepresentationClass,
    pub index_palette: IndexPalette,
    pub index_row: IndexRow,
}

impl RepresentationDeclaration {
    pub fn new(
        class: PortableRepresentationClass,
        index_palette: IndexPalette,
        index_row: IndexRow,
    ) -> Self {
        Self {
            class,
            index_palette,
            index_row,
        }
    }

    /// Recover portable metadata from a locally registered representation.
    pub fn from_symbol(symbol: Symbol) -> Result<Option<(String, Self)>> {
        if !symbol.has_tag(&SPENSO_TAG.representation) || symbol == LibraryRep::Dummy.symbol() {
            return Ok(None);
        }
        let metadata = RepresentationMetadata::from_symbol(symbol).ok_or_else(|| {
            RepresentationAttachmentError::new(format!(
                "representation symbol {} has no portable Spenso metadata",
                symbol.get_name()
            ))
        })?;
        let name = representation_identity(symbol.get_name().as_bytes())?;
        Ok(Some((
            name,
            Self {
                class: metadata.class.into(),
                index_palette: metadata.index_palette,
                index_row: metadata.index_row,
            },
        )))
    }
}

/// Validated declarations, ordered by fully-qualified representation name.
///
/// [`Self::absorb_attachments`] is transactional.  It retains the raw wire
/// records as well as their decoded meaning, so two differently encoded values
/// for one attachment key fail closed instead of being accepted merely because
/// they decode to equal declarations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RepresentationDeclarations {
    raw_attachments: AttachmentSet,
    entries: BTreeMap<String, RepresentationDeclaration>,
}

impl RepresentationDeclarations {
    pub fn new() -> Self {
        Self::default()
    }

    /// Inspect all known representation attachments without touching global
    /// Symbolica or Spenso state. Unknown attachment schemas are ignored.
    pub fn from_attachment_set(attachments: &AttachmentSet) -> Result<Self> {
        let mut declarations = Self::new();
        declarations.absorb_attachments(attachments)?;
        Ok(declarations)
    }

    /// Merge another input's representation attachments transactionally.
    pub fn absorb_attachments(&mut self, attachments: &AttachmentSet) -> Result<()> {
        let mut candidate = self.clone();
        for attachment in attachments.iter() {
            if attachment.schema() != REPRESENTATION_ATTACHMENT_SCHEMA {
                continue;
            }
            if attachment.version() != REPRESENTATION_ATTACHMENT_VERSION {
                return Err(RepresentationAttachmentError::new(format!(
                    "unsupported {REPRESENTATION_ATTACHMENT_SCHEMA} attachment version {}; expected {REPRESENTATION_ATTACHMENT_VERSION}",
                    attachment.version()
                )));
            }

            let name = representation_identity(attachment.identity())?;
            let declaration = decode_representation_declaration(attachment.data())?;
            candidate
                .raw_attachments
                .insert(attachment.to_owned_attachment())
                .map_err(|error| {
                    RepresentationAttachmentError::new(format!(
                        "could not merge representation attachments: {error}"
                    ))
                })?;
            candidate.insert(name, declaration)?;
        }
        *self = candidate;
        Ok(())
    }

    /// Add a programmatically constructed declaration.
    pub fn insert(
        &mut self,
        name: impl Into<String>,
        declaration: RepresentationDeclaration,
    ) -> Result<bool> {
        let name = name.into();
        representation_identity(name.as_bytes())?;
        match self.entries.get(&name) {
            Some(existing) if existing == &declaration => Ok(false),
            Some(_) => Err(RepresentationAttachmentError::new(format!(
                "conflicting {REPRESENTATION_ATTACHMENT_SCHEMA} declarations for {name}"
            ))),
            None => {
                // Constructing the expected metadata also validates the label
                // derived from the short representation name.
                expected_representation_metadata(&name, &declaration)?;
                self.entries.insert(name, declaration);
                Ok(true)
            }
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn get(&self, name: &str) -> Option<&RepresentationDeclaration> {
        self.entries.get(name)
    }

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&str, &RepresentationDeclaration)> + DoubleEndedIterator
    {
        self.entries
            .iter()
            .map(|(name, declaration)| (name.as_str(), declaration))
    }

    /// Validate every collision with process-global Symbolica/Spenso state,
    /// but do not add custom representations yet.
    ///
    /// Builtins are initialized first so inline-metric declarations can refer
    /// to their local metric implementation.
    pub fn preflight_registration(&self) -> Result<()> {
        for (name, declaration) in &self.entries {
            expected_representation_metadata(name, declaration)?;
        }
        initialize_representations();
        for (name, declaration) in &self.entries {
            validate_representation_registration(name, declaration)?;
        }
        Ok(())
    }

    /// Register every declaration before the caller invokes `Atom::import`.
    ///
    /// All predictable collisions are checked before the first custom
    /// registration. Spenso's process-global registry is not transactional, so
    /// an unexpected error from the underlying registration cannot be rolled
    /// back.
    pub fn register_before_atom_import(&self) -> Result<()> {
        self.preflight_registration()?;
        for (name, declaration) in &self.entries {
            register_representation_declaration(name, declaration)?;
        }
        Ok(())
    }

    /// Add canonical current-version records for these declarations to an attachment set.
    pub fn append_attachments_to(&self, attachments: &mut AttachmentSet) -> Result<()> {
        for (name, declaration) in &self.entries {
            attachments
                .insert(representation_attachment(name, declaration)?)
                .map_err(|error| {
                    RepresentationAttachmentError::new(format!(
                        "could not encode representation attachment: {error}"
                    ))
                })?;
        }
        Ok(())
    }

    /// Encode only these representation declarations as an attachment set.
    pub fn to_attachment_set(&self) -> Result<AttachmentSet> {
        let mut attachments = AttachmentSet::new();
        self.append_attachments_to(&mut attachments)?;
        Ok(attachments)
    }

    /// Merge representation declarations referenced anywhere in an Atom.
    pub fn collect_from_atom(&mut self, atom: &Atom) -> Result<()> {
        self.collect_from_view(atom.as_view())
    }

    /// Build declarations for every representation referenced by an Atom.
    pub fn referenced_by_atom(atom: &Atom) -> Result<Self> {
        let mut declarations = Self::new();
        declarations.collect_from_atom(atom)?;
        Ok(declarations)
    }

    fn collect_from_view(&mut self, view: AtomView<'_>) -> Result<()> {
        match view {
            AtomView::Num(_) => Ok(()),
            AtomView::Var(variable) => self.collect_symbol(variable.get_symbol()),
            AtomView::Fun(function) => {
                self.collect_symbol(function.get_symbol())?;
                for argument in function.iter() {
                    self.collect_from_view(argument)?;
                }
                Ok(())
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.collect_from_view(base)?;
                self.collect_from_view(exponent)
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.collect_from_view(factor)?;
                }
                Ok(())
            }
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    self.collect_from_view(term)?;
                }
                Ok(())
            }
        }
    }

    fn collect_symbol(&mut self, symbol: Symbol) -> Result<()> {
        if let Some((name, declaration)) = RepresentationDeclaration::from_symbol(symbol)? {
            self.insert(name, declaration)?;
        }
        Ok(())
    }
}

/// A safe, runtime-independent display tree attached to one generated index
/// symbol.
///
/// The declaration is deliberately restricted to Spenso's [`IndexDisplay`]
/// nodes. It never contains Typst source, callbacks, or arbitrary Symbolica
/// user data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MathDisplayDeclaration {
    pub display: IndexDisplay,
}

impl MathDisplayDeclaration {
    /// Validate and construct a portable display declaration.
    pub fn new(display: IndexDisplay) -> Result<Self> {
        validate_index_display(&display)?;
        Ok(Self { display })
    }

    /// Recover a portable declaration from a locally registered generated
    /// display symbol.
    pub fn from_symbol(symbol: Symbol) -> Result<Option<(String, Self)>> {
        let generated_name = symbol
            .get_name()
            .rsplit("::")
            .next()
            .is_some_and(|name| name.starts_with(MATH_DISPLAY_SYMBOL_PREFIX));

        if !symbol.has_tag(&SPENSO_TAG.index) {
            if generated_name {
                return Err(RepresentationAttachmentError::new(format!(
                    "generated math-display symbol {} has no Spenso index tag",
                    symbol.get_name()
                )));
            }
            return Ok(None);
        }

        let Some(display) = IndexDisplay::from_symbol(symbol) else {
            // Ordinary Spenso index symbols also carry the index tag. Only a
            // generated math-display name promises structured display data.
            if generated_name {
                return Err(RepresentationAttachmentError::new(format!(
                    "generated math-display symbol {} has no structured Spenso display metadata",
                    symbol.get_name()
                )));
            }
            return Ok(None);
        };

        let declaration = Self::new(display)?;
        let name = math_display_identity(symbol.get_name().as_bytes())?;
        validate_math_display_symbol_identity(&name, &declaration)?;
        Ok(Some((name, declaration)))
    }
}

/// Validated structured math displays, ordered by their exact Symbolica symbol
/// identity.
///
/// As with [`RepresentationDeclarations`], attachment absorption is
/// transactional and retains raw records so alternate encodings for one key
/// fail closed.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MathDisplayDeclarations {
    raw_attachments: AttachmentSet,
    entries: BTreeMap<String, MathDisplayDeclaration>,
}

impl MathDisplayDeclarations {
    pub fn new() -> Self {
        Self::default()
    }

    /// Inspect all known math-display attachments without touching global
    /// Symbolica state. Unknown attachment schemas are ignored.
    pub fn from_attachment_set(attachments: &AttachmentSet) -> Result<Self> {
        let mut declarations = Self::new();
        declarations.absorb_attachments(attachments)?;
        Ok(declarations)
    }

    /// Merge another input's math-display attachments transactionally.
    pub fn absorb_attachments(&mut self, attachments: &AttachmentSet) -> Result<()> {
        let mut candidate = self.clone();
        for attachment in attachments.iter() {
            if attachment.schema() != MATH_DISPLAY_ATTACHMENT_SCHEMA {
                continue;
            }
            if attachment.version() != MATH_DISPLAY_ATTACHMENT_VERSION {
                return Err(RepresentationAttachmentError::new(format!(
                    "unsupported {MATH_DISPLAY_ATTACHMENT_SCHEMA} attachment version {}; expected {MATH_DISPLAY_ATTACHMENT_VERSION}",
                    attachment.version()
                )));
            }

            let name = math_display_identity(attachment.identity())?;
            let declaration = decode_math_display_declaration(attachment.data())?;
            candidate
                .raw_attachments
                .insert(attachment.to_owned_attachment())
                .map_err(|error| {
                    RepresentationAttachmentError::new(format!(
                        "could not merge math-display attachments: {error}"
                    ))
                })?;
            candidate.insert(name, declaration)?;
        }
        *self = candidate;
        Ok(())
    }

    /// Add a programmatically constructed declaration.
    pub fn insert(
        &mut self,
        name: impl Into<String>,
        declaration: MathDisplayDeclaration,
    ) -> Result<bool> {
        let name = name.into();
        validate_math_display_symbol_identity(&name, &declaration)?;
        match self.entries.get(&name) {
            Some(existing) if existing == &declaration => Ok(false),
            Some(_) => Err(RepresentationAttachmentError::new(format!(
                "conflicting {MATH_DISPLAY_ATTACHMENT_SCHEMA} declarations for {name}"
            ))),
            None => {
                self.entries.insert(name, declaration);
                Ok(true)
            }
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn get(&self, name: &str) -> Option<&MathDisplayDeclaration> {
        self.entries.get(name)
    }

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&str, &MathDisplayDeclaration)> + DoubleEndedIterator {
        self.entries
            .iter()
            .map(|(name, declaration)| (name.as_str(), declaration))
    }

    /// Validate every collision with process-global Symbolica state without
    /// registering the first new display symbol.
    pub fn preflight_registration(&self) -> Result<()> {
        for (name, declaration) in &self.entries {
            validate_math_display_registration(name, declaration)?;
        }
        Ok(())
    }

    /// Register every structured display symbol before `Atom::import`.
    pub fn register_before_atom_import(&self) -> Result<()> {
        self.preflight_registration()?;
        for (name, declaration) in &self.entries {
            register_math_display_declaration(name, declaration)?;
        }
        Ok(())
    }

    /// Add canonical current-version records for these declarations to an
    /// attachment set.
    pub fn append_attachments_to(&self, attachments: &mut AttachmentSet) -> Result<()> {
        for (name, declaration) in &self.entries {
            attachments
                .insert(math_display_attachment(name, declaration)?)
                .map_err(|error| {
                    RepresentationAttachmentError::new(format!(
                        "could not encode math-display attachment: {error}"
                    ))
                })?;
        }
        Ok(())
    }

    /// Encode only these math-display declarations as an attachment set.
    pub fn to_attachment_set(&self) -> Result<AttachmentSet> {
        let mut attachments = AttachmentSet::new();
        self.append_attachments_to(&mut attachments)?;
        Ok(attachments)
    }

    /// Merge declarations for generated math-display symbols referenced
    /// anywhere in an Atom.
    pub fn collect_from_atom(&mut self, atom: &Atom) -> Result<()> {
        self.collect_from_view(atom.as_view())
    }

    /// Build declarations for every generated math-display symbol referenced
    /// by an Atom.
    pub fn referenced_by_atom(atom: &Atom) -> Result<Self> {
        let mut declarations = Self::new();
        declarations.collect_from_atom(atom)?;
        Ok(declarations)
    }

    fn collect_from_view(&mut self, view: AtomView<'_>) -> Result<()> {
        match view {
            AtomView::Num(_) => Ok(()),
            AtomView::Var(variable) => self.collect_symbol(variable.get_symbol()),
            AtomView::Fun(function) => {
                self.collect_symbol(function.get_symbol())?;
                for argument in function.iter() {
                    self.collect_from_view(argument)?;
                }
                Ok(())
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                self.collect_from_view(base)?;
                self.collect_from_view(exponent)
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    self.collect_from_view(factor)?;
                }
                Ok(())
            }
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    self.collect_from_view(term)?;
                }
                Ok(())
            }
        }
    }

    fn collect_symbol(&mut self, symbol: Symbol) -> Result<()> {
        if let Some((name, declaration)) = MathDisplayDeclaration::from_symbol(symbol)? {
            self.insert(name, declaration)?;
        }
        Ok(())
    }
}

/// Validate and decode a canonical fully-qualified representation identity.
pub fn representation_identity(identity: &[u8]) -> Result<String> {
    qualified_symbol_identity(identity, "representation attachment")
}

/// Validate and decode a canonical fully-qualified math-display symbol
/// identity.
pub fn math_display_identity(identity: &[u8]) -> Result<String> {
    qualified_symbol_identity(identity, "math-display attachment")
}

fn qualified_symbol_identity(identity: &[u8], kind: &str) -> Result<String> {
    let name = str::from_utf8(identity).map_err(|_| {
        RepresentationAttachmentError::new(format!("{kind} identity must be UTF-8"))
    })?;
    let Some((namespace, short_name)) = name.rsplit_once("::") else {
        return Err(RepresentationAttachmentError::new(format!(
            "{kind} identity must be fully qualified"
        )));
    };
    if namespace.is_empty()
        || short_name.is_empty()
        || namespace.split("::").any(|part| part.is_empty())
        || name.trim() != name
        || name.chars().any(char::is_control)
    {
        return Err(RepresentationAttachmentError::new(format!(
            "{kind} identity {name:?} is not canonical"
        )));
    }
    Ok(name.to_owned())
}

/// Qualify a short name with `namespace`, then validate the result.
pub fn canonical_representation_name(name: &str, namespace: &str) -> Result<String> {
    let qualified = if name.contains("::") {
        name.to_owned()
    } else {
        format!("{namespace}::{name}")
    };
    representation_identity(qualified.as_bytes())
}

/// Return the deterministic short Symbolica name for a structured display.
///
/// The hash is stable across Spenso-aware runtimes and independent of their
/// presentation layer.
pub fn math_display_symbol_name(display: &IndexDisplay) -> Result<String> {
    validate_index_display(display)?;
    let mut hasher = blake3::Hasher::new();
    hasher.update(MATH_DISPLAY_HASH_DOMAIN);
    hash_index_display(display, &mut hasher);
    Ok(format!(
        "{MATH_DISPLAY_SYMBOL_PREFIX}{}",
        hasher.finalize().to_hex()
    ))
}

/// Qualify the deterministic display-symbol name with `namespace`.
pub fn canonical_math_display_symbol_name(
    display: &IndexDisplay,
    namespace: &str,
) -> Result<String> {
    let qualified = format!("{namespace}::{}", math_display_symbol_name(display)?);
    math_display_identity(qualified.as_bytes())
}

/// Register one structured display under its deterministic symbol identity.
///
/// The same declaration can later be collected from an Atom and restored in
/// another Spenso-aware runtime through `spenso.math-display` attachments.
pub fn register_math_display_symbol(display: &IndexDisplay, namespace: &str) -> Result<Symbol> {
    let declaration = MathDisplayDeclaration::new(display.clone())?;
    let name = canonical_math_display_symbol_name(display, namespace)?;
    validate_math_display_registration(&name, &declaration)?;
    register_math_display_declaration(&name, &declaration)
}

fn index_row_from_wire(value: &str) -> Result<IndexRow> {
    match value {
        "top" => Ok(IndexRow::Top),
        "bottom" => Ok(IndexRow::Bottom),
        other => Err(RepresentationAttachmentError::new(format!(
            "unknown representation index row {other:?}"
        ))),
    }
}

/// Encode the canonical CBOR DATA field `[class, palette, index-row]`.
pub fn encode_representation_declaration(
    declaration: &RepresentationDeclaration,
) -> Result<Vec<u8>> {
    let value = Value::Array(vec![
        Value::Text(declaration.class.as_str().to_owned()),
        index_palette_to_wire(&declaration.index_palette)?,
        Value::Text(declaration.index_row.as_str().to_owned()),
    ]);
    let mut output = Vec::new();
    ciborium::into_writer(&value, &mut output).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode representation attachment data: {error}"
        ))
    })?;
    Ok(output)
}

/// Decode the canonical DATA field, rejecting trailing bytes.
pub fn decode_representation_declaration(input: &[u8]) -> Result<RepresentationDeclaration> {
    let mut cursor = Cursor::new(input);
    let value = ciborium::from_reader::<Value, _>(&mut cursor).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "representation attachment data must be CBOR-encoded: {error}"
        ))
    })?;
    if cursor.position() != input.len() as u64 {
        return Err(RepresentationAttachmentError::new(
            "representation attachment data has trailing bytes",
        ));
    }
    let Value::Array(fields) = value else {
        return Err(RepresentationAttachmentError::new(
            "representation attachment data must be an array",
        ));
    };
    let [Value::Text(class), palette, Value::Text(index_row)] = fields.as_slice() else {
        return Err(RepresentationAttachmentError::new(
            "representation attachment data must contain exactly class, index palette, and index row",
        ));
    };
    Ok(RepresentationDeclaration {
        class: PortableRepresentationClass::from_str(class)?,
        index_palette: index_palette_from_wire(palette)?,
        index_row: index_row_from_wire(index_row)?,
    })
}

/// Construct one canonical current-version `spenso.representation` attachment.
pub fn representation_attachment(
    name: &str,
    declaration: &RepresentationDeclaration,
) -> Result<Attachment> {
    let name = representation_identity(name.as_bytes())?;
    let key = AttachmentKey::new(
        REPRESENTATION_ATTACHMENT_SCHEMA,
        REPRESENTATION_ATTACHMENT_VERSION,
        name.into_bytes(),
    )
    .map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode representation attachment: {error}"
        ))
    })?;
    Attachment::new(key, encode_representation_declaration(declaration)?).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode representation attachment: {error}"
        ))
    })
}

/// Encode a structured display as canonical CBOR.
pub fn encode_math_display_declaration(declaration: &MathDisplayDeclaration) -> Result<Vec<u8>> {
    let mut nodes = 0;
    let value = index_display_to_wire(&declaration.display, 0, &mut nodes)?;
    let mut output = Vec::new();
    ciborium::into_writer(&value, &mut output).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode math-display attachment data: {error}"
        ))
    })?;
    Ok(output)
}

/// Decode a structured display from exact CBOR, rejecting trailing bytes.
pub fn decode_math_display_declaration(input: &[u8]) -> Result<MathDisplayDeclaration> {
    let mut cursor = Cursor::new(input);
    let value = ciborium::from_reader::<Value, _>(&mut cursor).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "math-display attachment data must be CBOR-encoded: {error}"
        ))
    })?;
    if cursor.position() != input.len() as u64 {
        return Err(RepresentationAttachmentError::new(
            "math-display attachment data has trailing bytes",
        ));
    }
    let mut nodes = 0;
    MathDisplayDeclaration::new(index_display_from_wire(&value, 0, &mut nodes)?)
}

/// Construct one canonical `spenso.math-display` attachment.
pub fn math_display_attachment(
    name: &str,
    declaration: &MathDisplayDeclaration,
) -> Result<Attachment> {
    validate_math_display_symbol_identity(name, declaration)?;
    let key = AttachmentKey::new(
        MATH_DISPLAY_ATTACHMENT_SCHEMA,
        MATH_DISPLAY_ATTACHMENT_VERSION,
        name.as_bytes().to_vec(),
    )
    .map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode math-display attachment: {error}"
        ))
    })?;
    Attachment::new(key, encode_math_display_declaration(declaration)?).map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not encode math-display attachment: {error}"
        ))
    })
}

fn validate_index_display(display: &IndexDisplay) -> Result<()> {
    let mut nodes = 0;
    index_display_to_wire(display, 0, &mut nodes).map(|_| ())
}

fn hash_index_display(display: &IndexDisplay, hasher: &mut blake3::Hasher) {
    match display {
        IndexDisplay::Symbol(name) => {
            hasher.update(&[0]);
            hasher.update(&(name.len() as u64).to_le_bytes());
            hasher.update(name.as_bytes());
        }
        IndexDisplay::Text(value) => {
            // New tags are appended after the original v1 variants so
            // Symbol/Number/Sequence/Attach identities remain stable.
            hasher.update(&[4]);
            hasher.update(&(value.len() as u64).to_le_bytes());
            hasher.update(value.as_bytes());
        }
        IndexDisplay::Number(number) => {
            hasher.update(&[1]);
            hasher.update(&number.to_le_bytes());
        }
        IndexDisplay::Sequence(items) => {
            hasher.update(&[2]);
            hasher.update(&(items.len() as u64).to_le_bytes());
            for item in items {
                hash_index_display(item, hasher);
            }
        }
        IndexDisplay::List(items) => {
            hasher.update(&[5]);
            hasher.update(&(items.len() as u64).to_le_bytes());
            for item in items {
                hash_index_display(item, hasher);
            }
        }
        IndexDisplay::Math { head, arguments } => {
            hasher.update(&[6]);
            hasher.update(&(head.len() as u64).to_le_bytes());
            hasher.update(head.as_bytes());
            hasher.update(&(arguments.len() as u64).to_le_bytes());
            for argument in arguments {
                hash_index_display(argument, hasher);
            }
        }
        IndexDisplay::Attach { base, top, bottom } => {
            hasher.update(&[3]);
            hash_index_display(base, hasher);
            match top {
                Some(top) => {
                    hasher.update(&[1]);
                    hash_index_display(top, hasher);
                }
                None => {
                    hasher.update(&[0]);
                }
            }
            match bottom {
                Some(bottom) => {
                    hasher.update(&[1]);
                    hash_index_display(bottom, hasher);
                }
                None => {
                    hasher.update(&[0]);
                }
            }
        }
    }
}

fn index_display_to_wire(display: &IndexDisplay, depth: usize, nodes: &mut usize) -> Result<Value> {
    if depth > MAX_DISPLAY_INDEX_DEPTH {
        return Err(RepresentationAttachmentError::new(format!(
            "index display exceeds the maximum depth of {MAX_DISPLAY_INDEX_DEPTH}"
        )));
    }
    *nodes += 1;
    if *nodes > MAX_DISPLAY_INDEX_NODES {
        return Err(RepresentationAttachmentError::new(format!(
            "index display exceeds the maximum size of {MAX_DISPLAY_INDEX_NODES} nodes"
        )));
    }
    match display {
        IndexDisplay::Symbol(name) => {
            IndexDisplay::symbol(name.clone())
                .map_err(|error| RepresentationAttachmentError::new(error.to_string()))?;
            Ok(Value::Array(vec![
                Value::Text("symbol".to_owned()),
                Value::Text(name.clone()),
            ]))
        }
        IndexDisplay::Text(value) => {
            IndexDisplay::text(value.clone())
                .map_err(|error| RepresentationAttachmentError::new(error.to_string()))?;
            Ok(Value::Array(vec![
                Value::Text("text".to_owned()),
                Value::Text(value.clone()),
            ]))
        }
        IndexDisplay::Number(number) => Ok(Value::Array(vec![
            Value::Text("number".to_owned()),
            Value::Integer((*number).into()),
        ])),
        IndexDisplay::Sequence(items) => {
            if items.is_empty() || items.len() > MAX_DISPLAY_INDEX_NODES {
                return Err(RepresentationAttachmentError::new(
                    "index display sequence must contain 1 to 64 items",
                ));
            }
            Ok(Value::Array(vec![
                Value::Text("sequence".to_owned()),
                Value::Array(
                    items
                        .iter()
                        .map(|item| index_display_to_wire(item, depth + 1, nodes))
                        .collect::<Result<Vec<_>>>()?,
                ),
            ]))
        }
        IndexDisplay::List(items) => {
            if items.len() > MAX_DISPLAY_INDEX_NODES {
                return Err(RepresentationAttachmentError::new(
                    "index display list cannot contain more than 64 items",
                ));
            }
            Ok(Value::Array(vec![
                Value::Text("list".to_owned()),
                Value::Array(
                    items
                        .iter()
                        .map(|item| index_display_to_wire(item, depth + 1, nodes))
                        .collect::<Result<Vec<_>>>()?,
                ),
            ]))
        }
        IndexDisplay::Math { head, arguments } => {
            IndexDisplay::math(head.clone(), arguments.clone())
                .map_err(|error| RepresentationAttachmentError::new(error.to_string()))?;
            Ok(Value::Array(vec![
                Value::Text("math".to_owned()),
                Value::Text(head.clone()),
                Value::Array(
                    arguments
                        .iter()
                        .map(|argument| index_display_to_wire(argument, depth + 1, nodes))
                        .collect::<Result<Vec<_>>>()?,
                ),
            ]))
        }
        IndexDisplay::Attach { base, top, bottom } => Ok(Value::Array(vec![
            Value::Text("attach".to_owned()),
            index_display_to_wire(base, depth + 1, nodes)?,
            top.as_deref()
                .map(|top| index_display_to_wire(top, depth + 1, nodes))
                .transpose()?
                .unwrap_or(Value::Null),
            bottom
                .as_deref()
                .map(|bottom| index_display_to_wire(bottom, depth + 1, nodes))
                .transpose()?
                .unwrap_or(Value::Null),
        ])),
    }
}

fn index_display_from_wire(value: &Value, depth: usize, nodes: &mut usize) -> Result<IndexDisplay> {
    if depth > MAX_DISPLAY_INDEX_DEPTH {
        return Err(RepresentationAttachmentError::new(format!(
            "index display exceeds the maximum depth of {MAX_DISPLAY_INDEX_DEPTH}"
        )));
    }
    *nodes += 1;
    if *nodes > MAX_DISPLAY_INDEX_NODES {
        return Err(RepresentationAttachmentError::new(format!(
            "index display exceeds the maximum size of {MAX_DISPLAY_INDEX_NODES} nodes"
        )));
    }
    let Value::Array(fields) = value else {
        return Err(RepresentationAttachmentError::new(
            "index display must be an array",
        ));
    };
    let Some(Value::Text(kind)) = fields.first() else {
        return Err(RepresentationAttachmentError::new(
            "index display must start with a text kind",
        ));
    };
    match (kind.as_str(), fields.as_slice()) {
        ("symbol", [_, Value::Text(name)]) => IndexDisplay::symbol(name.clone())
            .map_err(|error| RepresentationAttachmentError::new(error.to_string())),
        ("text", [_, Value::Text(value)]) => IndexDisplay::text(value.clone())
            .map_err(|error| RepresentationAttachmentError::new(error.to_string())),
        ("number", [_, number]) => Ok(IndexDisplay::Number(value_i64(
            number,
            "representation palette number",
        )?)),
        ("sequence", [_, Value::Array(items)]) => {
            if items.is_empty() || items.len() > MAX_DISPLAY_INDEX_NODES {
                return Err(RepresentationAttachmentError::new(
                    "index display sequence must contain 1 to 64 items",
                ));
            }
            Ok(IndexDisplay::Sequence(
                items
                    .iter()
                    .map(|item| index_display_from_wire(item, depth + 1, nodes))
                    .collect::<Result<Vec<_>>>()?,
            ))
        }
        ("list", [_, Value::Array(items)]) => {
            if items.len() > MAX_DISPLAY_INDEX_NODES {
                return Err(RepresentationAttachmentError::new(
                    "index display list cannot contain more than 64 items",
                ));
            }
            Ok(IndexDisplay::List(
                items
                    .iter()
                    .map(|item| index_display_from_wire(item, depth + 1, nodes))
                    .collect::<Result<Vec<_>>>()?,
            ))
        }
        ("math", [_, Value::Text(head), Value::Array(arguments)]) => {
            let arguments = arguments
                .iter()
                .map(|argument| index_display_from_wire(argument, depth + 1, nodes))
                .collect::<Result<Vec<_>>>()?;
            IndexDisplay::math(head.clone(), arguments)
                .map_err(|error| RepresentationAttachmentError::new(error.to_string()))
        }
        ("attach", [_, base, top, bottom]) => Ok(IndexDisplay::Attach {
            base: Box::new(index_display_from_wire(base, depth + 1, nodes)?),
            top: match top {
                Value::Null => None,
                value => Some(Box::new(index_display_from_wire(value, depth + 1, nodes)?)),
            },
            bottom: match bottom {
                Value::Null => None,
                value => Some(Box::new(index_display_from_wire(value, depth + 1, nodes)?)),
            },
        }),
        _ => Err(RepresentationAttachmentError::new(format!(
            "invalid index display node {kind:?}"
        ))),
    }
}

fn index_palette_to_wire(palette: &IndexPalette) -> Result<Value> {
    match palette {
        IndexPalette::Numeric => Ok(Value::Array(vec![Value::Text("numeric".to_owned())])),
        IndexPalette::Cyclic { start, labels } => {
            let start = i64::try_from(*start).map_err(|_| {
                RepresentationAttachmentError::new(
                    "representation palette start must fit in a signed 64-bit integer",
                )
            })?;
            if labels.is_empty() || labels.len() > MAX_DISPLAY_INDEX_NODES {
                return Err(RepresentationAttachmentError::new(
                    "a cyclic index palette must contain 1 to 64 labels",
                ));
            }
            let mut nodes = 0;
            Ok(Value::Array(vec![
                Value::Text("cyclic".to_owned()),
                Value::Integer(start.into()),
                Value::Array(
                    labels
                        .iter()
                        .map(|label| index_display_to_wire(label, 0, &mut nodes))
                        .collect::<Result<Vec<_>>>()?,
                ),
            ]))
        }
    }
}

fn index_palette_from_wire(value: &Value) -> Result<IndexPalette> {
    let Value::Array(fields) = value else {
        return Err(RepresentationAttachmentError::new(
            "representation index palette must be an array",
        ));
    };
    match fields.as_slice() {
        [Value::Text(kind)] if kind == "numeric" => Ok(IndexPalette::Numeric),
        [Value::Text(kind), start, Value::Array(labels)] if kind == "cyclic" => {
            let start = value_i64(start, "representation palette start")?;
            let start = usize::try_from(start).map_err(|_| {
                RepresentationAttachmentError::new(
                    "representation palette start must be non-negative",
                )
            })?;
            let mut nodes = 0;
            IndexPalette::cyclic(
                start,
                labels
                    .iter()
                    .map(|label| index_display_from_wire(label, 0, &mut nodes))
                    .collect::<Result<Vec<_>>>()?,
            )
            .map_err(|error| RepresentationAttachmentError::new(error.to_string()))
        }
        _ => Err(RepresentationAttachmentError::new(
            "invalid representation index palette",
        )),
    }
}

fn value_i64(value: &Value, label: &str) -> Result<i64> {
    match value {
        Value::Integer(value) => (*value)
            .try_into()
            .map_err(|_| RepresentationAttachmentError::new(format!("{label} is out of range"))),
        other => Err(RepresentationAttachmentError::new(format!(
            "{label} must be an integer, got {other:?}"
        ))),
    }
}

fn validate_math_display_symbol_identity(
    name: &str,
    declaration: &MathDisplayDeclaration,
) -> Result<()> {
    math_display_identity(name.as_bytes())?;
    validate_index_display(&declaration.display)?;
    let actual_short_name = name.rsplit("::").next().ok_or_else(|| {
        RepresentationAttachmentError::new(format!("invalid math-display symbol identity {name:?}"))
    })?;
    let expected_short_name = math_display_symbol_name(&declaration.display)?;
    if actual_short_name != expected_short_name {
        return Err(RepresentationAttachmentError::new(format!(
            "math-display symbol {name} does not match its deterministic display identity {expected_short_name}"
        )));
    }
    Ok(())
}

fn validate_math_display_registration(
    name: &str,
    declaration: &MathDisplayDeclaration,
) -> Result<()> {
    validate_math_display_symbol_identity(name, declaration)?;
    if let Some(existing) = Symbol::get_symbol(NamespacedSymbol::parse(name))
        && (!existing.has_tag(&SPENSO_TAG.index)
            || IndexDisplay::from_symbol(existing) != Some(declaration.display.clone()))
    {
        return Err(RepresentationAttachmentError::new(format!(
            "math-display symbol {name} already exists with conflicting metadata"
        )));
    }
    Ok(())
}

fn register_math_display_declaration(
    name: &str,
    declaration: &MathDisplayDeclaration,
) -> Result<Symbol> {
    validate_math_display_registration(name, declaration)?;
    let namespaced = NamespacedSymbol::parse(name);
    if let Some(existing) = Symbol::get_symbol(namespaced.clone()) {
        return Ok(existing);
    }

    SymbolBuilder::new(namespaced)
        .with_tags([SPENSO_TAG.index.clone()])
        .with_user_data(declaration.display.symbol_user_data())
        .build()
        .map_err(|error| {
            RepresentationAttachmentError::new(format!(
                "could not register math-display symbol {name}: {error}"
            ))
        })
}

fn expected_representation_metadata(
    name: &str,
    declaration: &RepresentationDeclaration,
) -> Result<RepresentationMetadata> {
    // Validate programmatically constructed palettes as strictly as decoded
    // wire data. `IndexPalette` and `IndexDisplay` are public enums, so callers
    // are not required to use their checked constructors.
    index_palette_to_wire(&declaration.index_palette)?;
    let label = name.rsplit("::").next().ok_or_else(|| {
        RepresentationAttachmentError::new(format!("invalid representation name {name:?}"))
    })?;
    Ok(RepresentationMetadata {
        class: declaration.class.into(),
        label: IndexDisplay::symbol(label)
            .map_err(|error| RepresentationAttachmentError::new(error.to_string()))?,
        index_palette: declaration.index_palette.clone(),
        index_row: declaration.index_row,
    })
}

fn validate_representation_registration(
    name: &str,
    declaration: &RepresentationDeclaration,
) -> Result<()> {
    let namespaced = NamespacedSymbol::parse(name);
    if let Some(existing) = Symbol::get_symbol(namespaced) {
        let expected = expected_representation_metadata(name, declaration)?;
        let correct_tags = existing.has_tag(&SPENSO_TAG.representation)
            && match declaration.class {
                PortableRepresentationClass::SelfDual
                | PortableRepresentationClass::InlineMetric => {
                    existing.has_tag(&SPENSO_TAG.self_dual)
                        && !existing.has_tag(&SPENSO_TAG.dualizable)
                }
                PortableRepresentationClass::Dualizable => {
                    existing.has_tag(&SPENSO_TAG.dualizable)
                        && !existing.has_tag(&SPENSO_TAG.self_dual)
                }
            };
        if !correct_tags || RepresentationMetadata::from_symbol(existing) != Some(expected) {
            return Err(RepresentationAttachmentError::new(format!(
                "representation symbol {name} already exists with conflicting metadata"
            )));
        }
        if existing.get_print_function().is_none() {
            return Err(RepresentationAttachmentError::new(format!(
                "representation symbol {name} was imported before its local print callback could be registered"
            )));
        }
        if declaration.class == PortableRepresentationClass::InlineMetric {
            LibraryRep::try_from_symbol_coerced(existing).map_err(|error| {
                RepresentationAttachmentError::new(format!(
                    "inline-metric representation {name} is not locally registered: {error}"
                ))
            })?;
        }
    } else if declaration.class == PortableRepresentationClass::InlineMetric {
        return Err(RepresentationAttachmentError::new(format!(
            "inline-metric representation {name} is not a locally registered builtin"
        )));
    }
    Ok(())
}

fn register_representation_declaration(
    name: &str,
    declaration: &RepresentationDeclaration,
) -> Result<LibraryRep> {
    match declaration.class {
        PortableRepresentationClass::SelfDual => {
            LibraryRep::new_self_dual_with_index_palette_and_row(
                name,
                declaration.index_palette.clone(),
                declaration.index_row,
            )
        }
        PortableRepresentationClass::Dualizable => LibraryRep::new_dual_with_index_palette_and_row(
            name,
            declaration.index_palette.clone(),
            declaration.index_row,
        ),
        PortableRepresentationClass::InlineMetric => {
            let symbol = Symbol::get_symbol(NamespacedSymbol::parse(name)).ok_or_else(|| {
                RepresentationAttachmentError::new(format!(
                    "inline-metric representation {name} is not a locally registered builtin"
                ))
            })?;
            return LibraryRep::try_from_symbol_coerced(symbol).map_err(|error| {
                RepresentationAttachmentError::new(format!(
                    "inline-metric representation {name} is not locally registered: {error}"
                ))
            });
        }
    }
    .map_err(|error| {
        RepresentationAttachmentError::new(format!(
            "could not register representation {name}: {error}"
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cyclic_declaration() -> RepresentationDeclaration {
        RepresentationDeclaration::new(
            PortableRepresentationClass::SelfDual,
            IndexPalette::cyclic(
                1,
                [
                    IndexDisplay::symbol("mu").unwrap(),
                    IndexDisplay::symbol("nu")
                        .unwrap()
                        .with_bottom(IndexDisplay::Number(2)),
                ],
            )
            .unwrap(),
            IndexRow::Bottom,
        )
    }

    fn structured_math_display() -> MathDisplayDeclaration {
        let matrix = IndexDisplay::math(
            "mat",
            vec![
                IndexDisplay::List(vec![
                    IndexDisplay::symbol("alpha").unwrap(),
                    IndexDisplay::Number(1),
                ]),
                IndexDisplay::List(vec![
                    IndexDisplay::text("rate").unwrap(),
                    IndexDisplay::math(
                        "frac",
                        vec![
                            IndexDisplay::symbol("x").unwrap(),
                            IndexDisplay::symbol("y").unwrap(),
                        ],
                    )
                    .unwrap(),
                ]),
            ],
        )
        .unwrap();
        MathDisplayDeclaration::new(IndexDisplay::Attach {
            base: Box::new(matrix),
            top: Some(Box::new(IndexDisplay::Sequence(vec![
                IndexDisplay::symbol("mu").unwrap(),
                IndexDisplay::Number(2),
            ]))),
            bottom: None,
        })
        .unwrap()
    }

    #[test]
    fn representation_wire_round_trip_preserves_class_palette_and_row() {
        let declaration = cyclic_declaration();
        let encoded = encode_representation_declaration(&declaration).unwrap();
        assert_eq!(
            decode_representation_declaration(&encoded).unwrap(),
            declaration
        );
    }

    #[test]
    fn representation_wire_rejects_unknown_index_rows() {
        let value = Value::Array(vec![
            Value::Text("self-dual".to_owned()),
            Value::Array(vec![Value::Text("numeric".to_owned())]),
            Value::Text("middle".to_owned()),
        ]);
        let mut encoded = Vec::new();
        ciborium::into_writer(&value, &mut encoded).unwrap();

        assert!(
            decode_representation_declaration(&encoded)
                .unwrap_err()
                .to_string()
                .contains("unknown representation index row")
        );
    }

    #[test]
    fn declaration_data_rejects_trailing_bytes() {
        let mut encoded = encode_representation_declaration(&cyclic_declaration()).unwrap();
        encoded.push(0);
        assert!(
            decode_representation_declaration(&encoded)
                .unwrap_err()
                .to_string()
                .contains("trailing bytes")
        );
    }

    #[test]
    fn attachment_absorption_is_transactional_on_raw_conflict() {
        let name = "example::M";
        let first = cyclic_declaration();
        let second = RepresentationDeclaration::new(
            PortableRepresentationClass::Dualizable,
            IndexPalette::Numeric,
            IndexRow::Top,
        );
        let first_set =
            AttachmentSet::from_attachments([representation_attachment(name, &first).unwrap()])
                .unwrap();
        let second_set =
            AttachmentSet::from_attachments([representation_attachment(name, &second).unwrap()])
                .unwrap();

        let mut declarations = RepresentationDeclarations::from_attachment_set(&first_set).unwrap();
        assert!(declarations.absorb_attachments(&second_set).is_err());
        assert_eq!(declarations.get(name), Some(&first));
    }

    #[test]
    fn unsupported_known_versions_are_rejected_but_unknown_schemas_are_ignored() {
        let declaration = cyclic_declaration();
        let unsupported = Attachment::new(
            AttachmentKey::new(
                REPRESENTATION_ATTACHMENT_SCHEMA,
                REPRESENTATION_ATTACHMENT_VERSION + 1,
                b"example::M".to_vec(),
            )
            .unwrap(),
            encode_representation_declaration(&declaration).unwrap(),
        )
        .unwrap();
        let set = AttachmentSet::from_attachments([unsupported]).unwrap();
        assert!(RepresentationDeclarations::from_attachment_set(&set).is_err());

        let unknown = Attachment::new(
            AttachmentKey::new("example.unknown", 9, b"x".to_vec()).unwrap(),
            vec![1, 2, 3],
        )
        .unwrap();
        let set = AttachmentSet::from_attachments([unknown]).unwrap();
        assert!(
            RepresentationDeclarations::from_attachment_set(&set)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn identities_must_be_fully_qualified_and_canonical() {
        assert!(representation_identity(b"example::M").is_ok());
        assert!(representation_identity(b"M").is_err());
        assert!(representation_identity(b"example::::M").is_err());
        assert!(representation_identity(b" example::M").is_err());
    }

    #[test]
    fn math_display_wire_round_trip_preserves_every_structured_variant() {
        let declaration = structured_math_display();
        let encoded = encode_math_display_declaration(&declaration).unwrap();

        assert_eq!(
            decode_math_display_declaration(&encoded).unwrap(),
            declaration
        );
    }

    #[test]
    fn math_display_wire_rejects_trailing_bytes_and_unapproved_heads() {
        let mut encoded = encode_math_display_declaration(&structured_math_display()).unwrap();
        encoded.push(0);
        assert!(
            decode_math_display_declaration(&encoded)
                .unwrap_err()
                .to_string()
                .contains("trailing bytes")
        );

        let forged = Value::Array(vec![
            Value::Text("math".to_owned()),
            Value::Text("raw-typst".to_owned()),
            Value::Array(vec![]),
        ]);
        let mut encoded = Vec::new();
        ciborium::into_writer(&forged, &mut encoded).unwrap();
        assert!(
            decode_math_display_declaration(&encoded)
                .unwrap_err()
                .to_string()
                .contains("unsupported mathematical display node")
        );
    }

    #[test]
    fn math_display_symbol_identity_is_deterministic_and_structural() {
        let symbol = IndexDisplay::symbol("mu").unwrap();
        let text = IndexDisplay::text("mu").unwrap();
        let first = canonical_math_display_symbol_name(&symbol, "display_identity_test").unwrap();
        let second = canonical_math_display_symbol_name(&symbol, "display_identity_test").unwrap();
        let text_name = canonical_math_display_symbol_name(&text, "display_identity_test").unwrap();

        assert_eq!(first, second);
        assert_ne!(first, text_name);
        assert!(
            first
                .rsplit("::")
                .next()
                .unwrap()
                .starts_with(MATH_DISPLAY_SYMBOL_PREFIX)
        );

        let declaration = MathDisplayDeclaration::new(text).unwrap();
        assert!(math_display_attachment(&first, &declaration).is_err());
    }

    #[test]
    fn math_display_attachments_restore_local_spenso_metadata() {
        let declaration = structured_math_display();
        let name = canonical_math_display_symbol_name(
            &declaration.display,
            "math_display_attachment_restore_test",
        )
        .unwrap();
        let attachments =
            AttachmentSet::from_attachments(
                [math_display_attachment(&name, &declaration).unwrap()],
            )
            .unwrap();
        let declarations = MathDisplayDeclarations::from_attachment_set(&attachments).unwrap();

        assert_eq!(declarations.get(&name), Some(&declaration));
        declarations.register_before_atom_import().unwrap();

        let symbol = Symbol::get_symbol(NamespacedSymbol::parse(&name)).unwrap();
        assert!(symbol.has_tag(&SPENSO_TAG.index));
        assert_eq!(
            IndexDisplay::from_symbol(symbol),
            Some(declaration.display.clone())
        );
    }

    #[test]
    fn math_display_declarations_are_collected_from_nested_atoms() {
        let declaration = structured_math_display();
        let symbol =
            register_math_display_symbol(&declaration.display, "math_display_atom_collection_test")
                .unwrap();
        let atom = symbol.call_args([Atom::num(1), Atom::var(symbol)]);
        let declarations = MathDisplayDeclarations::referenced_by_atom(&atom).unwrap();
        let name = symbol.get_name();

        assert_eq!(declarations.len(), 1);
        assert_eq!(declarations.get(name), Some(&declaration));
        let attachments = declarations.to_attachment_set().unwrap();
        assert_eq!(attachments.len(), 1);
        assert_eq!(
            MathDisplayDeclarations::from_attachment_set(&attachments).unwrap(),
            declarations
        );
    }

    #[test]
    fn math_display_registration_collisions_fail_before_mutation() {
        let expected =
            MathDisplayDeclaration::new(IndexDisplay::symbol("expected").unwrap()).unwrap();
        let conflicting =
            MathDisplayDeclaration::new(IndexDisplay::symbol("conflicting").unwrap()).unwrap();
        let expected_name =
            canonical_math_display_symbol_name(&expected.display, "math_display_collision_test")
                .unwrap();

        SymbolBuilder::new(NamespacedSymbol::parse(&expected_name))
            .with_tags([SPENSO_TAG.index.clone()])
            .with_user_data(conflicting.display.symbol_user_data())
            .build()
            .unwrap();

        let mut declarations = MathDisplayDeclarations::new();
        declarations.insert(expected_name, expected).unwrap();
        assert!(
            declarations
                .register_before_atom_import()
                .unwrap_err()
                .to_string()
                .contains("conflicting metadata")
        );
    }

    #[test]
    fn unsupported_math_display_versions_are_rejected() {
        let declaration = structured_math_display();
        let name =
            canonical_math_display_symbol_name(&declaration.display, "math_display_version_test")
                .unwrap();
        let attachment = Attachment::new(
            AttachmentKey::new(
                MATH_DISPLAY_ATTACHMENT_SCHEMA,
                MATH_DISPLAY_ATTACHMENT_VERSION + 1,
                name.into_bytes(),
            )
            .unwrap(),
            encode_math_display_declaration(&declaration).unwrap(),
        )
        .unwrap();
        let attachments = AttachmentSet::from_attachments([attachment]).unwrap();

        assert!(MathDisplayDeclarations::from_attachment_set(&attachments).is_err());
    }
}
