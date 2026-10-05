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

impl Model {
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
        let mut labels = LABELS.write().unwrap();
        let texname = texname.filter(|name| !name.trim().is_empty());
        let typstname = typstname.filter(|name| !name.trim().is_empty());
        if texname.is_some() || typstname.is_some() {
            let latex = texname.unwrap_or(name);
            labels.insert(
                symbol,
                ParameterLabel {
                    latex: latex.to_owned(),
                    typst: typstname
                        .map(str::to_owned)
                        .unwrap_or_else(|| format!("{name:?}")),
                },
            );
        } else {
            labels.remove(&symbol);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Barrier};
    use symbolica::{
        atom::{Atom, AtomCore},
        printer::PrintOptions,
    };

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
