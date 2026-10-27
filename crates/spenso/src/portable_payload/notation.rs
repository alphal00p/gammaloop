//! Explicit, bounded transport of process-local display notation.

use super::{
    RepresentationAttachmentError, index_display_from_wire, index_display_to_wire,
    qualified_symbol_identity,
};
use ciborium::de::Error as CborDecodeError;
use ciborium::ser::Error as CborEncodeError;
use ciborium::value::Value;
use std::collections::BTreeMap;
use std::io::Cursor;
use std::io::Error as IoError;
use std::sync::Arc;
use symbolica::atom::{Atom, AtomCore};
use symbolica_notation::{
    DisplayText, DisplayTextError, MathMlSubset, MathMlSubsetError, Notation, RegistryError,
    TextOverrides,
};
use symbolica_typst_plugin::payload::{Attachment, AttachmentKey, AttachmentSet, PayloadError};
use thiserror::Error;

/// Opt-in notation schema; see [`NotationDeclarations`] for usage.
pub const NOTATION_ATTACHMENT_SCHEMA: &str = "spenso.notation";
/// Current schema version.
pub const NOTATION_ATTACHMENT_VERSION: u32 = 1;
const MAX_DECLARATIONS: usize = 128;
const MAX_DECLARATION_BYTES: usize = 64 * 1024;

/// Invalid wire data or a conflict with local notation.
#[derive(Debug, Error)]
pub enum NotationAttachmentError {
    /// Invalid common attachment envelope.
    #[error(transparent)]
    Payload(#[from] PayloadError),
    /// Malformed CBOR, retaining its decoding source.
    #[error("invalid notation CBOR: {0}")]
    Decode(#[from] CborDecodeError<IoError>),
    /// Failure to encode a declaration.
    #[error("could not encode notation CBOR: {0}")]
    Encode(#[from] CborEncodeError<IoError>),
    /// Invalid text override.
    #[error(transparent)]
    Text(#[from] DisplayTextError),
    /// Invalid presentation MathML.
    #[error(transparent)]
    MathMl(#[from] MathMlSubsetError),
    /// Malformed or excessive data.
    #[error("invalid notation attachment: {0}")]
    Invalid(String),
    /// A different notation is already registered locally.
    #[error("notation for {name}: {source}")]
    Registry {
        name: String,
        #[source]
        source: RegistryError,
    },
    /// Invalid structured display.
    #[error(transparent)]
    Display(#[from] RepresentationAttachmentError),
}
type Result<T> = std::result::Result<T, NotationAttachmentError>;

/// Validated notation declarations, using the common payload attachment set.
///
/// ```
/// # use spenso::portable_payload::{NotationDeclarations, NotationAttachmentError};
/// # use symbolica::atom::Atom;
/// let declarations = NotationDeclarations::referenced_by_atom(&Atom::num(1))?;
/// assert!(declarations.to_attachment_set()?.is_empty());
/// # Ok::<(), NotationAttachmentError>(())
/// ```
#[derive(Clone, Debug, Default)]
pub struct NotationDeclarations {
    raw: AttachmentSet,
    entries: BTreeMap<String, Arc<Notation>>,
}

impl NotationDeclarations {
    /// Decode and validate every notation record without registering symbols.
    /// See [`Self`] for an export example.
    pub fn from_attachment_set(attachments: &AttachmentSet) -> Result<Self> {
        let mut result = Self::default();
        result.absorb_attachments(attachments)?;
        Ok(result)
    }

    /// Merge records transactionally, rejecting conflicting wire encodings.
    pub fn absorb_attachments(&mut self, attachments: &AttachmentSet) -> Result<()> {
        let mut candidate = self.clone();
        for attachment in attachments
            .iter()
            .filter(|a| a.schema() == NOTATION_ATTACHMENT_SCHEMA)
        {
            if attachment.version() != NOTATION_ATTACHMENT_VERSION {
                return Err(NotationAttachmentError::Invalid(
                    "unsupported spenso.notation version".into(),
                ));
            }
            let name = Self::identity(attachment.identity())?;
            let notation = Self::decode(attachment.data())?;
            candidate.raw.insert(attachment.to_owned_attachment())?;
            candidate.entries.insert(name, Arc::new(notation));
            candidate.check_count()?;
        }
        *self = candidate;
        Ok(())
    }

    /// Collect notation only for symbols actually referenced by the atom.
    pub fn referenced_by_atom(atom: &Atom) -> Result<Self> {
        let mut result = Self::default();
        for symbol in atom.get_all_symbols(true) {
            if let Some(notation) =
                Notation::lookup(symbol).map_err(|source| NotationAttachmentError::Registry {
                    name: symbol.get_name().to_owned(),
                    source,
                })?
            {
                let name = Self::identity(symbol.get_name().as_bytes())?;
                // Encode and decode here too: public IndexDisplay trees can bypass constructors.
                let notation = Self::decode(&Self::encode(&notation)?)?;
                result.entries.insert(name, Arc::new(notation));
                result.check_count()?;
            }
        }
        Ok(result)
    }

    /// Check every existing local registration before installing any declaration.
    pub fn preflight_registration(&self) -> Result<()> {
        for (name, notation) in &self.entries {
            let existing = Notation::lookup_name(name).map_err(|source| {
                NotationAttachmentError::Registry {
                    name: name.to_owned(),
                    source,
                }
            })?;
            if existing
                .as_deref()
                .is_some_and(|existing| !Self::same_display(existing, notation))
            {
                return Err(NotationAttachmentError::Registry {
                    name: name.to_owned(),
                    source: RegistryError::NameConflict {
                        name: name.to_owned(),
                    },
                });
            }
        }
        Ok(())
    }

    /// Install validated notation before Atom import; identical imports are idempotent.
    pub fn register_before_atom_import(&self) -> Result<()> {
        self.preflight_registration()?;
        let declarations = self
            .entries
            .iter()
            .map(|(name, notation)| {
                // Preserve trusted local callback ownership. Foreign declarations
                // never enable callbacks or replace an existing registry handle.
                let existing = Notation::lookup_name(name).map_err(|source| {
                    NotationAttachmentError::Registry {
                        name: name.to_owned(),
                        source,
                    }
                })?;
                if existing
                    .as_deref()
                    .is_some_and(|existing| !Self::same_display(existing, notation))
                {
                    return Err(NotationAttachmentError::Registry {
                        name: name.to_owned(),
                        source: RegistryError::NameConflict {
                            name: name.to_owned(),
                        },
                    });
                }
                Ok((
                    name.to_owned(),
                    existing.unwrap_or_else(|| Arc::clone(notation)),
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        Notation::register_names(declarations).map_err(|source| NotationAttachmentError::Registry {
            name: "portable declarations".into(),
            source,
        })
    }

    fn same_display(left: &Notation, right: &Notation) -> bool {
        left.index_label == right.index_label
            && left.head == right.head
            && left.text.plain.as_ref().map(DisplayText::as_str)
                == right.text.plain.as_ref().map(DisplayText::as_str)
            && left.text.latex.as_ref().map(DisplayText::as_str)
                == right.text.latex.as_ref().map(DisplayText::as_str)
            && left.text.typst.as_ref().map(DisplayText::as_str)
                == right.text.typst.as_ref().map(DisplayText::as_str)
            && left.mathml == right.mathml
    }

    /// Append canonical records to the common attachment set.
    pub fn append_attachments_to(&self, attachments: &mut AttachmentSet) -> Result<()> {
        for (name, notation) in &self.entries {
            let key = AttachmentKey::new(
                NOTATION_ATTACHMENT_SCHEMA,
                NOTATION_ATTACHMENT_VERSION,
                name.as_bytes().to_vec(),
            )?;
            let attachment = Attachment::new(key, Self::encode(notation)?)?;
            attachments.insert(attachment)?;
        }
        Ok(())
    }

    /// Encode this declaration set. See [`Self`] for an example.
    pub fn to_attachment_set(&self) -> Result<AttachmentSet> {
        let mut attachments = AttachmentSet::new();
        self.append_attachments_to(&mut attachments)?;
        Ok(attachments)
    }

    fn check_count(&self) -> Result<()> {
        if self.entries.len() > MAX_DECLARATIONS {
            return Err(NotationAttachmentError::Invalid(
                "more than 128 declarations".into(),
            ));
        }
        Ok(())
    }

    fn identity(bytes: &[u8]) -> Result<String> {
        let name = qualified_symbol_identity(bytes, "notation")?;
        if name.len() > 4096
            || name.split("::").any(|part| {
                let mut chars = part.chars();
                !chars.next().is_some_and(|c| c.is_alphabetic() || c == '_')
                    || !chars.all(|c| c.is_alphanumeric() || c == '_')
            })
        {
            return Err(NotationAttachmentError::Invalid(
                "invalid symbol-name syntax".into(),
            ));
        }
        Ok(name)
    }

    fn encode(notation: &Notation) -> Result<Vec<u8>> {
        let mut nodes = 0;
        let mut fields = Vec::new();
        for (key, tree) in [
            ("index-label", &notation.index_label),
            ("head", &notation.head),
        ] {
            if let Some(tree) = tree {
                fields.push((
                    Value::Text(key.into()),
                    index_display_to_wire(tree, 0, &mut nodes)?,
                ));
            }
        }
        for (key, text) in [
            ("plain", &notation.text.plain),
            ("latex", &notation.text.latex),
            ("typst", &notation.text.typst),
        ] {
            if let Some(text) = text {
                fields.push((Value::Text(key.into()), Value::Text(text.as_str().into())));
            }
        }
        if let Some(mathml) = &notation.mathml {
            let source = mathml.to_mathml();
            // to_mathml supplies a framing mrow; transport its children to preserve equality.
            fields.push((
                Value::Text("mathml".into()),
                Value::Text(source[6..source.len() - 7].into()),
            ));
        }
        // Callback implementations and ownership remain process-local.
        let mut bytes = Vec::new();
        ciborium::into_writer(&Value::Map(fields), &mut bytes)?;
        if bytes.len() > MAX_DECLARATION_BYTES {
            return Err(NotationAttachmentError::Invalid(
                "notation exceeds 64 KiB".into(),
            ));
        }
        Ok(bytes)
    }

    fn decode(bytes: &[u8]) -> Result<Notation> {
        if bytes.len() > MAX_DECLARATION_BYTES {
            return Err(NotationAttachmentError::Invalid(
                "notation exceeds 64 KiB".into(),
            ));
        }
        let mut cursor = Cursor::new(bytes);
        let value: Value = ciborium::from_reader(&mut cursor)?;
        if cursor.position() != bytes.len() as u64 {
            return Err(NotationAttachmentError::Invalid("trailing bytes".into()));
        }
        let Value::Map(fields) = value else {
            return Err(NotationAttachmentError::Invalid(
                "expected notation map".into(),
            ));
        };
        let mut result = Notation {
            text: TextOverrides::default(),
            ..Notation::default()
        };
        let mut seen = Vec::new();
        let mut nodes = 0;
        for (key, value) in fields {
            let Value::Text(key) = key else {
                return Err(NotationAttachmentError::Invalid(
                    "field name must be text".into(),
                ));
            };
            if seen.contains(&key) {
                return Err(NotationAttachmentError::Invalid("duplicate field".into()));
            }
            match key.as_str() {
                "index-label" => {
                    result.index_label = Some(index_display_from_wire(&value, 0, &mut nodes)?)
                }
                "head" => result.head = Some(index_display_from_wire(&value, 0, &mut nodes)?),
                "plain" | "latex" | "typst" | "mathml" => {
                    let Value::Text(text) = value else {
                        return Err(NotationAttachmentError::Invalid(
                            "override must be text".into(),
                        ));
                    };
                    if key == "mathml" {
                        result.mathml = Some(MathMlSubset::parse(&text)?);
                        seen.push(key);
                        continue;
                    }
                    let text = DisplayText::literal(text)?;
                    match key.as_str() {
                        "plain" => result.text.plain = Some(text),
                        "latex" => result.text.latex = Some(text),
                        _ => result.text.typst = Some(text),
                    }
                }
                _ => {
                    return Err(NotationAttachmentError::Invalid(format!(
                        "unknown field {key}"
                    )));
                }
            }
            seen.push(key);
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::portable_payload;
    use symbolica::atom::{NamespacedSymbol, Symbol, SymbolBuilder};
    use symbolica_notation::IndexDisplay;

    fn sample() -> Notation {
        Notation {
            index_label: Some(IndexDisplay::symbol("mu").unwrap()),
            head: Some(IndexDisplay::symbol("eta").unwrap()),
            text: TextOverrides {
                plain: Some(DisplayText::literal("portable label").unwrap()),
                latex: Some(DisplayText::literal("L").unwrap()),
                typst: Some(DisplayText::literal("L").unwrap()),
            },
            mathml: Some(MathMlSubset::parse("<mi mathvariant='bold'>L</mi>").unwrap()),
            custom_print: false,
        }
    }

    fn wire(name: &str, bytes: Vec<u8>) -> AttachmentSet {
        AttachmentSet::from_attachments([Attachment::new(
            AttachmentKey::new(NOTATION_ATTACHMENT_SCHEMA, 1, name.as_bytes().to_vec()).unwrap(),
            bytes,
        )
        .unwrap()])
        .unwrap()
    }

    #[test]
    fn encoding_round_trip_preserves_all_fields() {
        let notation = sample();
        assert_eq!(
            NotationDeclarations::decode(&NotationDeclarations::encode(&notation).unwrap())
                .unwrap(),
            notation
        );
    }

    #[test]
    fn callback_ownership_stays_local_and_cannot_be_imported() {
        let name = "portable_notation_tests::local_callback";
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(name))
            .build()
            .unwrap();
        let mut local = sample();
        local.custom_print = true;
        local.text.typst = Some(DisplayText::new("L").unwrap());
        let original = local.register(symbol).unwrap();
        let declarations = NotationDeclarations::referenced_by_atom(&Atom::var(symbol)).unwrap();
        assert!(!declarations.entries[name].custom_print);
        declarations.register_before_atom_import().unwrap();
        assert!(Arc::ptr_eq(
            &original,
            &Notation::lookup(symbol).unwrap().unwrap()
        ));
        let mut bytes = Vec::new();
        ciborium::into_writer(
            &Value::Map(vec![(
                Value::Text("custom-print".into()),
                Value::Bool(true),
            )]),
            &mut bytes,
        )
        .unwrap();
        assert!(NotationDeclarations::decode(&bytes).is_err());
    }

    #[test]
    fn imported_text_never_becomes_target_source() {
        let mut bytes = Vec::new();
        ciborium::into_writer(
            &Value::Map(vec![
                (
                    Value::Text("typst".into()),
                    Value::Text("#read(\"secret\")".into()),
                ),
                (
                    Value::Text("latex".into()),
                    Value::Text("\\input{secret}".into()),
                ),
            ]),
            &mut bytes,
        )
        .unwrap();
        let notation = NotationDeclarations::decode(&bytes).unwrap();
        assert_eq!(
            notation.text.typst.unwrap().to_typst_source(),
            "upright(\"#read(\\\"secret\\\")\")"
        );
        assert_eq!(
            notation.text.latex.unwrap().to_latex_source(),
            "\\text{\\textbackslash{}input\\{secret\\}}"
        );
    }

    #[test]
    fn plain_atom_export_does_not_serialize_registered_notation() {
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(
            "portable_notation_tests::plain_export",
        ))
        .build()
        .unwrap();
        let atom = Atom::var(symbol);
        let mut before = Vec::new();
        atom.export(&mut before).unwrap();
        sample().register(symbol).unwrap();
        let mut after = Vec::new();
        atom.export(&mut after).unwrap();
        assert_eq!(after, before);
        assert!(
            !portable_payload::attachments_for_atom(&atom)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn registers_unknown_symbol_and_reimports_idempotently() {
        let name = "portable_notation_tests::fresh";
        assert!(Symbol::get_symbol(NamespacedSymbol::parse(name)).is_none());
        let notation = sample();
        let attachments = wire(name, NotationDeclarations::encode(&notation).unwrap());
        let declarations = NotationDeclarations::from_attachment_set(&attachments).unwrap();
        declarations.register_before_atom_import().unwrap();
        declarations.register_before_atom_import().unwrap();
        assert!(Symbol::get_symbol(NamespacedSymbol::parse(name)).is_none());
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(name))
            .build()
            .unwrap();
        assert_eq!(
            Notation::lookup(symbol).unwrap().as_deref(),
            Some(&notation)
        );
        let exported = NotationDeclarations::referenced_by_atom(&Atom::var(symbol)).unwrap();
        assert_eq!(exported.to_attachment_set().unwrap(), attachments);
    }

    #[test]
    fn conflicts_preflight_before_any_registration() {
        let name = "portable_notation_tests::z_conflict";
        let symbol = SymbolBuilder::new(NamespacedSymbol::parse(name))
            .build()
            .unwrap();
        Notation::default().register(symbol).unwrap();
        let mut attachments = wire(
            "portable_notation_tests::a_untouched",
            NotationDeclarations::encode(&sample()).unwrap(),
        );
        attachments
            .merge(&wire(
                name,
                NotationDeclarations::encode(&sample()).unwrap(),
            ))
            .unwrap();
        let declarations = NotationDeclarations::from_attachment_set(&attachments).unwrap();
        assert!(
            matches!(declarations.register_before_atom_import(), Err(NotationAttachmentError::Registry { name: conflict, .. }) if conflict == name || conflict == "portable declarations")
        );
        assert!(
            Symbol::get_symbol(NamespacedSymbol::parse(
                "portable_notation_tests::a_untouched"
            ))
            .is_none()
        );
    }

    #[test]
    fn rejects_untrusted_records_without_registering_anything() {
        for (key, value) in [
            ("mathml", "<mi onclick='boom'>X</mi>".to_owned()),
            ("plain", "\0".to_owned()),
            ("typst", "x".repeat(4097)),
        ] {
            let mut bytes = Vec::new();
            ciborium::into_writer(
                &Value::Map(vec![(Value::Text(key.into()), Value::Text(value))]),
                &mut bytes,
            )
            .unwrap();
            let mut attachments = wire(
                "portable_notation_tests::untouched",
                NotationDeclarations::encode(&sample()).unwrap(),
            );
            attachments
                .merge(&wire("portable_notation_tests::invalid", bytes))
                .unwrap();
            assert!(NotationDeclarations::from_attachment_set(&attachments).is_err());
        }
        assert!(
            Symbol::get_symbol(NamespacedSymbol::parse(
                "portable_notation_tests::untouched"
            ))
            .is_none()
        );
        assert!(NotationDeclarations::decode(&vec![0; MAX_DECLARATION_BYTES + 1]).is_err());
        // A declared array count of u64::MAX must fail on bounded input.
        assert!(
            NotationDeclarations::decode(&[0x9b, 255, 255, 255, 255, 255, 255, 255, 255]).is_err()
        );
        for name in ["bare", "x::::y", "x::bad name", "x::<img>", "x::a\0"] {
            assert!(NotationDeclarations::identity(name.as_bytes()).is_err());
        }
    }

    #[test]
    fn rejects_excessive_declaration_count_and_unknown_versions() {
        let mut attachments = AttachmentSet::new();
        for index in 0..=MAX_DECLARATIONS {
            attachments
                .merge(&wire(
                    &format!("portable_notation_count::s{index}"),
                    vec![0xa0],
                ))
                .unwrap();
        }
        assert!(NotationDeclarations::from_attachment_set(&attachments).is_err());
        assert!(
            Symbol::get_symbol(NamespacedSymbol::parse("portable_notation_count::s0")).is_none()
        );
        let attachments = AttachmentSet::from_attachments([Attachment::new(
            AttachmentKey::new(
                NOTATION_ATTACHMENT_SCHEMA,
                2,
                b"test::unknown_version".to_vec(),
            )
            .unwrap(),
            vec![0xa0],
        )
        .unwrap()])
        .unwrap();
        assert!(NotationDeclarations::from_attachment_set(&attachments).is_err());
    }
}
