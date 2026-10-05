use std::cell::RefCell;
use std::{
    collections::HashMap,
    sync::{LazyLock, RwLock},
};

use symbolica::atom::{AtomView, NamespacedSymbol, Symbol, SymbolBuilder};

use crate::{EntityKind, Model, ModelError};

struct ParameterLabel {
    latex: String,
    typst: String,
}

// Symbols have process-wide identities. Read labels at print time so importing
// a model again updates its presentation without redefining its symbols.
static LABELS: LazyLock<RwLock<HashMap<Symbol, ParameterLabel>>> =
    LazyLock::new(|| RwLock::new(HashMap::new()));

// Each synchronous import stages only its own writes. Nested successful imports
// feed their parent; a failed outer operation drops all of its pending changes.
thread_local! {
    static LABEL_UPDATES: RefCell<Vec<HashMap<Symbol, Option<ParameterLabel>>>> =
        const { RefCell::new(Vec::new()) };
}

struct ParameterDisplayTransaction {
    active: bool,
}

impl ParameterDisplayTransaction {
    fn begin() -> Self {
        LABEL_UPDATES.with(|updates| updates.borrow_mut().push(HashMap::new()));
        Self { active: true }
    }

    fn record(symbol: Symbol, label: Option<ParameterLabel>) {
        let immediate = LABEL_UPDATES.with(|updates| {
            let mut updates = updates.borrow_mut();
            if let Some(current) = updates.last_mut() {
                current.insert(symbol, label);
                None
            } else {
                Some((symbol, label))
            }
        });
        if let Some(update) = immediate {
            Self::publish([update]);
        }
    }

    fn publish(updates: impl IntoIterator<Item = (Symbol, Option<ParameterLabel>)>) {
        let mut labels = LABELS.write().unwrap();
        for (symbol, label) in updates {
            if let Some(label) = label {
                labels.insert(symbol, label);
            } else {
                labels.remove(&symbol);
            }
        }
    }

    fn commit(mut self) {
        let ready = LABEL_UPDATES.with(|updates| {
            let mut updates = updates.borrow_mut();
            let completed = updates
                .pop()
                .expect("parameter display transaction is active");
            self.active = false;
            if let Some(parent) = updates.last_mut() {
                parent.extend(completed);
                None
            } else {
                Some(completed)
            }
        });
        if let Some(ready) = ready {
            Self::publish(ready);
        }
    }
}

impl Drop for ParameterDisplayTransaction {
    fn drop(&mut self) {
        if self.active {
            LABEL_UPDATES.with(|updates| {
                updates
                    .borrow_mut()
                    .pop()
                    .expect("parameter display transaction is active");
            });
        }
    }
}

impl Model {
    /// Publish model-owned parameter labels only if a synchronous operation succeeds.
    ///
    /// Use this around model construction and application-level validation when
    /// the application can reject an otherwise valid model. Existing expressions
    /// keep their committed labels until the operation returns successfully.
    /// Nested successful operations join their parent's transaction. A failure
    /// or panic discards pending updates without reverting other threads' work.
    ///
    /// Parsing still interns immutable symbols and their printing callbacks.
    /// This transaction stages label values on the calling thread; worker-thread
    /// model imports have their own independent registration scope.
    pub fn with_parameter_display_transaction<T, E>(
        operation: impl FnOnce() -> Result<T, E>,
    ) -> Result<T, E> {
        let transaction = ParameterDisplayTransaction::begin();
        let result = operation()?;
        transaction.commit();
        Ok(result)
    }

    /// Declare a parameter's printer before parsing model expressions.
    ///
    /// Import adapters which parse expressions before constructing a `Model`
    /// must call this for every parameter first, including forward references.
    /// Repeated declarations update labels; missing/empty labels restore the
    /// ordinary spelling. Existing user-declared symbol printers are preserved.
    pub fn register_parameter_symbol(
        name: &str,
        texname: Option<&str>,
        typstname: Option<&str>,
    ) -> Result<(), ModelError> {
        let symbol_name = NamespacedSymbol::parse(&format!("UFO::{name}"));
        let symbol = if let Some(symbol) = Symbol::get_symbol(symbol_name.clone()) {
            // Existing user declarations (including their printers) are immutable.
            symbol
        } else {
            SymbolBuilder::new(symbol_name.clone())
                .with_print_function(|view, options, _| {
                    if !options.mode.is_latex() && !options.mode.is_typst() {
                        return None;
                    }
                    let AtomView::Var(variable) = view else {
                        return None;
                    };
                    let labels = LABELS.read().unwrap();
                    let label = labels.get(&variable.get_symbol())?;
                    Some(if options.mode.is_latex() {
                        label.latex.clone()
                    } else {
                        label.typst.clone()
                    })
                })
                .build()
                // Another import can create this same immutable symbol after
                // the first lookup. Reuse that declaration exactly as in the
                // fast path, preserving an existing user printer as well.
                .or_else(|message| Symbol::get_symbol(symbol_name).ok_or(message))
                .map_err(|message| ModelError::SymbolicParse {
                    kind: EntityKind::Parameter,
                    name: name.to_owned(),
                    field: "name",
                    expression: name.to_owned(),
                    message: message.to_string(),
                })?
        };
        let texname = texname.filter(|name| !name.trim().is_empty());
        let typstname = typstname.filter(|name| !name.trim().is_empty());
        if texname.is_some() || typstname.is_some() {
            let latex = texname.unwrap_or(name);
            ParameterDisplayTransaction::record(
                symbol,
                Some(ParameterLabel {
                    latex: latex.to_owned(),
                    typst: typstname
                        .map(str::to_owned)
                        .unwrap_or_else(|| format!("{name:?}")),
                }),
            );
        } else {
            ParameterDisplayTransaction::record(symbol, None);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::panic;
    use std::sync::{Arc, Barrier};
    use std::thread;
    use symbolica::{
        atom::{Atom, AtomCore},
        printer::PrintOptions,
    };

    #[test]
    fn failed_parameter_display_transaction_keeps_committed_labels() {
        let name = "failed_display_transaction";
        Model::register_parameter_symbol(name, Some("before"), Some("before-typst")).unwrap();
        let atom = Atom::var(
            Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap(),
        );
        let outcome = Model::with_parameter_display_transaction(|| {
            Model::register_parameter_symbol(name, Some("pending"), Some("pending-typst")).unwrap();
            assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "before");
            assert_eq!(
                atom.printer(PrintOptions::typst()).to_string(),
                "before-typst"
            );
            Err::<(), _>("model replacement rejected")
        });
        assert!(outcome.is_err());
        assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "before");
        assert_eq!(
            atom.printer(PrintOptions::typst()).to_string(),
            "before-typst"
        );
    }

    #[test]
    fn nested_parameter_display_transactions_publish_only_their_successful_writes() {
        let name = "nested_display_transaction";
        Model::register_parameter_symbol(name, Some("before"), None).unwrap();
        let atom = Atom::var(
            Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap(),
        );
        Model::with_parameter_display_transaction(|| {
            Model::register_parameter_symbol(name, Some("outer"), None)?;
            let rejected: Result<(), ModelError> =
                Model::with_parameter_display_transaction(|| {
                    Model::register_parameter_symbol(name, Some("rejected inner"), None)?;
                    Err(ModelError::SymbolicParse {
                        kind: EntityKind::Parameter,
                        name: name.into(),
                        field: "test",
                        expression: name.into(),
                        message: "reject inner transaction".into(),
                    })
                });
            assert!(rejected.is_err());
            Model::with_parameter_display_transaction(|| {
                Model::register_parameter_symbol(name, Some("accepted inner"), None)
            })?;
            assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "before");
            Ok::<_, ModelError>(())
        })
        .unwrap();
        assert_eq!(
            atom.printer(PrintOptions::latex()).to_string(),
            "accepted inner"
        );
        Model::with_parameter_display_transaction(|| {
            Model::register_parameter_symbol(name, None, None)
        })
        .unwrap();
        assert!(
            !atom
                .printer(PrintOptions::latex())
                .to_string()
                .contains("accepted inner")
        );
    }

    #[test]
    fn successful_inner_display_transaction_is_discarded_with_its_outer_import() {
        let name = "outer_rejected_display_transaction";
        Model::register_parameter_symbol(name, Some("before"), None).unwrap();
        let atom = Atom::var(
            Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap(),
        );
        let outcome = Model::with_parameter_display_transaction(|| {
            Model::with_parameter_display_transaction(|| {
                Model::register_parameter_symbol(name, Some("inner"), None)
            })
            .unwrap();
            Err::<(), _>("outer import rejected")
        });
        assert!(outcome.is_err());
        assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "before");
    }

    #[test]
    fn rejected_parameter_display_transaction_preserves_concurrent_commits() {
        let name = "concurrent_display_transaction";
        Model::register_parameter_symbol(name, Some("before"), None).unwrap();
        let atom = Atom::var(
            Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap(),
        );
        let outcome = Model::with_parameter_display_transaction(|| {
            Model::register_parameter_symbol(name, Some("rejected"), None).unwrap();
            thread::spawn(move || {
                Model::register_parameter_symbol(name, Some("concurrent committed"), None).unwrap();
            })
            .join()
            .unwrap();
            assert_eq!(
                atom.printer(PrintOptions::latex()).to_string(),
                "concurrent committed"
            );
            Err::<(), _>("local validation failed")
        });
        assert!(outcome.is_err());
        assert_eq!(
            atom.printer(PrintOptions::latex()).to_string(),
            "concurrent committed"
        );
    }

    #[test]
    fn panicking_parameter_display_transaction_releases_its_staged_writes() {
        let name = "panicked_display_transaction";
        Model::register_parameter_symbol(name, Some("before"), None).unwrap();
        let atom = Atom::var(
            Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap(),
        );
        let caught = panic::catch_unwind(|| {
            let _ = Model::with_parameter_display_transaction(|| -> Result<(), ModelError> {
                Model::register_parameter_symbol(name, Some("pending"), None)?;
                panic!("reject the model while unwinding");
            });
        });
        assert!(caught.is_err());
        assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "before");
        Model::register_parameter_symbol(name, Some("after"), None).unwrap();
        assert_eq!(atom.printer(PrintOptions::latex()).to_string(), "after");
    }

    #[test]
    fn concurrent_parameter_registration_preserves_process_wide_symbols() {
        let barrier = Arc::new(Barrier::new(16));
        let workers = (0..16)
            .map(|_| {
                let barrier = Arc::clone(&barrier);
                std::thread::spawn(move || {
                    let mut errors = Vec::new();
                    for round in 0..64 {
                        barrier.wait();
                        if let Err(error) = Model::register_parameter_symbol(
                            &format!("concurrent_model_parameter_{round}"),
                            None,
                            None,
                        ) {
                            errors.push(error.to_string());
                        }
                    }
                    errors
                })
            })
            .collect::<Vec<_>>();
        let errors = workers
            .into_iter()
            .flat_map(|thread| thread.join().unwrap())
            .collect::<Vec<_>>();
        assert!(errors.is_empty(), "{errors:?}");
    }

    #[test]
    fn parameter_registration_keeps_existing_user_printer_and_invalid_name_errors() {
        let name = NamespacedSymbol::parse("UFO::preserved_parameter_printer_test");
        let symbol = SymbolBuilder::new(name)
            .with_print_function(|_, _, _| Some("user-printer".into()))
            .build()
            .unwrap();
        Model::register_parameter_symbol(
            "preserved_parameter_printer_test",
            Some("model label"),
            None,
        )
        .unwrap();
        assert_eq!(
            Atom::var(symbol).printer(PrintOptions::latex()).to_string(),
            "user-printer"
        );
        assert!(Model::register_parameter_symbol("invalid parameter name", None, None).is_err());
    }

    #[test]
    fn parameter_reregistration_retains_independent_latex_and_typst_labels() {
        let name = "reregistered_dual_label_parameter";
        Model::register_parameter_symbol(name, Some("old latex"), Some("old typst")).unwrap();
        let symbol = Symbol::get_symbol(NamespacedSymbol::parse(&format!("UFO::{name}"))).unwrap();
        Model::register_parameter_symbol(name, Some("new latex"), Some("new typst")).unwrap();
        let expression = Atom::var(symbol);
        assert_eq!(
            expression.printer(PrintOptions::latex()).to_string(),
            "new latex"
        );
        assert_eq!(
            expression.printer(PrintOptions::typst()).to_string(),
            "new typst"
        );
    }
}
