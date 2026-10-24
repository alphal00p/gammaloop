use crate::IndexDisplay;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::sync::{Arc, LazyLock, RwLock};
use symbolica::atom::Symbol;
use thiserror::Error;

/// Bounded display text, independent of any target's evaluation machinery.
///
/// This validates size and characters, not LaTeX or Typst syntax. Consumers must
/// treat it as display text and apply their target's escaping or validation.
///
/// ```
/// use symbolica_notation::{DisplayText, DisplayTextError};
/// let text = DisplayText::new("h")?;
/// assert_eq!(text.as_str(), "h");
/// assert!(DisplayText::new("\0").is_err());
/// # Ok::<(), DisplayTextError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DisplayText(String);

impl DisplayText {
    /// Validate at most 4096 UTF-8 bytes without control or XML noncharacters.
    ///
    /// Empty text is valid and explicitly suppresses a label.
    /// See [`DisplayText`] for an example.
    ///
    /// # Errors
    /// Returns [`DisplayTextError`] when the bound or character contract fails.
    pub fn new(value: impl Into<String>) -> Result<Self, DisplayTextError> {
        let value = value.into();
        if value.len() > 4096
            || value
                .chars()
                .any(|c| c.is_control() || matches!(c, '\u{fffe}' | '\u{ffff}'))
        {
            return Err(DisplayTextError);
        }
        Ok(Self(value))
    }

    /// Borrow the validated text. See [`DisplayText`] for an example.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// Invalid size or characters in display text. See [`DisplayText`] for an example.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Error)]
#[error("display text must contain at most 4096 UTF-8 bytes without control or XML noncharacters")]
pub struct DisplayTextError;

/// Optional text overrides for each output target.
///
/// Absence delegates to structured notation or the escaped symbol name.
/// See [`Notation::register`] for an example constructing overrides.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct TextOverrides {
    /// Plain text. See [`Notation::register`].
    pub plain: Option<DisplayText>,
    /// LaTeX display text. See [`Notation::register`].
    pub latex: Option<DisplayText>,
    /// Typst display text. See [`Notation::register`].
    pub typst: Option<DisplayText>,
}

/// Display notation owned by this process, never attached to a symbol's state.
///
/// Index labels and structured heads are separate: `hedge(owner, local)` may
/// use the index label `h` while its ordinary function head has other notation.
/// The public structured trees retain [`IndexDisplay`]'s per-renderer bounds.
/// A future portable attachment may
/// explicitly export notation; symbol strings and state exports never do so
/// implicitly. See [`Self::register`] for a runnable example.
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct Notation {
    /// Optional structured label for an index owner. See [`Self::register`].
    pub index_label: Option<IndexDisplay>,
    /// Optional structured function/variable head. See [`Self::register`].
    pub head: Option<IndexDisplay>,
    /// Optional target-specific text. See [`Self::register`].
    pub text: TextOverrides,
    /// The symbol's callback owns the complete leaf, rather than delegating to
    /// the tensor renderer. This flag never enters serialized symbol state.
    /// See [`Self::register`] for registration and conflict behavior.
    pub custom_print: bool,
}

// The initializer only allocates an empty map. It never calls Symbolica, so it
// can also be forced inside a Symbolica initialize! callback without a probe or
// a symbol-bundle lazy wrapper. No Symbolica operation runs under this lock.
static NOTATION: LazyLock<RwLock<HashMap<Symbol, Arc<Notation>>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

/// Registration conflict or poisoned process-local storage.
///
/// See [`Notation::register`] and [`Notation::lookup`] for successful calls.
///
/// ```
/// use symbolica::symbol;
/// use symbolica_notation::{IndexDisplay, Notation, RegistryError};
/// let symbol = symbol!("notation_doc::conflict");
/// Notation::default().register(symbol)?;
/// let other = Notation { head: Some(IndexDisplay::Number(1)), ..Notation::default() };
/// assert!(matches!(other.register(symbol), Err(RegistryError::Conflict { symbol: s }) if s == symbol));
/// # Ok::<(), RegistryError>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum RegistryError {
    /// A different notation already owns this symbol. See [`RegistryError`].
    #[error("symbol {symbol:?} already has different process-local notation")]
    Conflict {
        /// Symbol whose existing notation was preserved. See [`RegistryError`].
        symbol: Symbol,
    },
    /// A prior panic poisoned storage; no further access is allowed.
    /// See [`Notation::lookup`] for the fallible access contract.
    #[error("process-local notation registry lock is poisoned")]
    Poisoned,
}

impl Notation {
    /// Whether a symbol declares a callback that owns the whole display leaf.
    ///
    /// # Errors
    /// Returns [`RegistryError::Poisoned`] if storage is poisoned.
    /// See [`Self::register`] for an example of the registry lifetime.
    pub fn has_custom_print(symbol: Symbol) -> Result<bool, RegistryError> {
        Ok(Self::lookup(symbol)?.is_some_and(|notation| notation.custom_print))
    }

    /// Register notation or return the existing shared handle when it is equal.
    ///
    /// The symbol must already belong to the current Symbolica state. The
    /// registry does not register symbols, mutate user data, or invoke printers.
    /// The lock covers only comparison and insertion; callers receive an owned
    /// shared handle and can invoke callbacks after this method returns.
    ///
    /// # Errors
    /// Returns [`RegistryError::Conflict`] without changing an existing entry,
    /// or [`RegistryError::Poisoned`] after a panic while holding the lock.
    ///
    /// ```
    /// use std::error::Error;
    /// use std::sync::Arc;
    /// use symbolica::symbol;
    /// use symbolica_notation::{DisplayText, IndexDisplay, Notation, TextOverrides};
    /// let symbol = symbol!("notation_doc::registered");
    /// let notation = Notation {
    ///     index_label: Some(IndexDisplay::symbol("h")?),
    ///     head: Some(IndexDisplay::symbol("eta")?),
    ///     text: TextOverrides {
    ///         plain: Some(DisplayText::new("hedge")?),
    ///         latex: Some(DisplayText::new("h")?),
    ///         typst: Some(DisplayText::new("h")?),
    ///     },
    ///     ..Notation::default()
    /// };
    /// let first = notation.register(symbol)?;
    /// let second = Notation::lookup(symbol)?.expect("registered above");
    /// assert!(Arc::ptr_eq(&first, &second));
    /// # Ok::<(), Box<dyn Error>>(())
    /// ```
    pub fn register(self, symbol: Symbol) -> Result<Arc<Self>, RegistryError> {
        let notation = Arc::new(self);
        let mut entries = NOTATION.write().map_err(|_| RegistryError::Poisoned)?;
        match entries.entry(symbol) {
            Entry::Occupied(entry) if entry.get().as_ref() == notation.as_ref() => {
                Ok(Arc::clone(entry.get()))
            }
            Entry::Occupied(_) => Err(RegistryError::Conflict { symbol }),
            Entry::Vacant(entry) => Ok(Arc::clone(entry.insert(notation))),
        }
    }

    /// Look up a shared immutable handle without invoking Symbolica or a printer.
    ///
    /// `None` means the symbol has no registered notation. Handles remain valid
    /// after releasing the read lock.
    ///
    /// # Errors
    /// Returns [`RegistryError::Poisoned`] if a prior panic poisoned storage.
    ///
    /// ```
    /// use symbolica::symbol;
    /// use symbolica_notation::{Notation, RegistryError};
    /// assert!(Notation::lookup(symbol!("notation_doc::unregistered"))?.is_none());
    /// # Ok::<(), RegistryError>(())
    /// ```
    pub fn lookup(symbol: Symbol) -> Result<Option<Arc<Self>>, RegistryError> {
        let entries = NOTATION.read().map_err(|_| RegistryError::Poisoned)?;
        Ok(entries.get(&symbol).map(Arc::clone))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use symbolica::symbol;

    #[test]
    fn registration_preserves_symbol_data() {
        let symbol = symbol!("notation_tests::registration");
        let before = symbol.get_data();
        let registered = Notation {
            index_label: Some(IndexDisplay::symbol("h").unwrap()),
            ..Notation::default()
        }
        .register(symbol)
        .unwrap();
        assert_eq!(
            registered.index_label.as_ref().unwrap().to_native_string(),
            "h"
        );
        assert!(Arc::ptr_eq(
            &registered,
            &Notation::lookup(symbol).unwrap().unwrap()
        ));
        assert_eq!(symbol.get_data(), before);
    }

    #[test]
    fn equal_registration_is_idempotent() {
        let symbol = symbol!("notation_tests::idempotent");
        let first = Notation::default().register(symbol).unwrap();
        let second = Notation::default().register(symbol).unwrap();
        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn conflicting_registration_retains_original() {
        let symbol = symbol!("notation_tests::conflict");
        let first = Notation::default().register(symbol).unwrap();
        let conflict = Notation {
            head: Some(IndexDisplay::Number(2)),
            ..Notation::default()
        };
        assert_eq!(
            conflict.register(symbol),
            Err(RegistryError::Conflict { symbol })
        );
        assert!(Arc::ptr_eq(
            &first,
            &Notation::lookup(symbol).unwrap().unwrap()
        ));
    }

    #[test]
    fn unregistered_symbol_has_no_notation() {
        assert!(
            Notation::lookup(symbol!("notation_tests::unregistered"))
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn concurrent_readers_share_the_registered_value() {
        let symbol = symbol!("notation_tests::concurrent");
        let notation = Notation::default().register(symbol).unwrap();
        thread::scope(|scope| {
            for _ in 0..8 {
                let expected = &notation;
                scope.spawn(move || {
                    for _ in 0..100 {
                        assert!(Arc::ptr_eq(
                            expected,
                            &Notation::lookup(symbol).unwrap().unwrap()
                        ));
                    }
                });
            }
        });
    }

    #[test]
    fn text_override_bounds_and_characters_are_validated() {
        assert!(DisplayText::new("").is_ok());
        assert!(DisplayText::new("é".repeat(2048)).is_ok());
        assert!(DisplayText::new("é".repeat(2049)).is_err());
        for value in ["\0", "\n", "\u{ffff}", "\u{fffe}"] {
            assert!(DisplayText::new(value).is_err());
        }
    }
}
