//! Canonical rank-one momentum heads shared by symbolic kinematics and diagrams.
//!
//! `Q` names a momentum family; graph callers assign its object labels to edges.
//! `K` and `P` name independent loop and external momentum coordinates.
//! The Rust owner does not change their serialized `gammalooprs` namespace.
//!
//! ```
//! use feynkit_kinematics::symbols;
//! use symbolica::atom::AtomCore;
//! use symbolica::printer::PrintOptions;
//! for (head, label) in [
//!     (symbols::momentum(), "q"),
//!     (symbols::loop_momentum(), "k"),
//!     (symbols::external_momentum(), "p"),
//! ] {
//!     assert_eq!(
//!         head.call(27).printer(PrintOptions::typst()).to_string(),
//!         format!("{label}_27"),
//!     );
//! }
//! ```
use spenso::network::tags;
use spenso::network::tags::{COMPONENT_LABELS_SUBSCRIPT_TAG, SPENSO_TAG};
use spenso::structure::representation::{Minkowski, RepName};
use spenso::symbolica_init;
use spenso::symbolica_init::SymbolicaInitLazy;
use std::sync::LazyLock;
use symbolica::atom::{Atom, Symbol};
use symbolica::{initialize, symbol};
use symbolica_notation::{IndexDisplay, Notation};

// Register before parsing can claim bare names. Preserve the initializer-aware
// lazy accessors and initialize the Lorentz representation before tensor heads.
initialize!(|| symbolica_init::in_symbolica_initializer(|| {
    let _ = Minkowski {}.to_symbolic([Atom::Zero]);
    momentum();
    loop_momentum();
    external_momentum();
}));

/// Generic momentum family; graph callers supply edge labels.
/// See the construction and presentation example in [`crate::symbols`].
pub fn momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        let symbol = symbol!(
            "gammalooprs::Q",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned()
            ]
        );
        Notation {
            head: Some(IndexDisplay::symbol("q").expect("the canonical label is valid")),
            ..Notation::default()
        }
        .register(symbol)
        .expect("the canonical graph symbol has one notation");
        symbol
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

/// Independent loop-momentum coordinates.
/// See the construction and presentation example in [`crate::symbols`].
pub fn loop_momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        let symbol = symbol!(
            "gammalooprs::K",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned()
            ]
        );
        Notation {
            head: Some(IndexDisplay::symbol("k").expect("the canonical label is valid")),
            ..Notation::default()
        }
        .register(symbol)
        .expect("the canonical graph symbol has one notation");
        symbol
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

/// Independent external-momentum coordinates.
/// See the construction and presentation example in [`crate::symbols`].
pub fn external_momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        let symbol = symbol!(
            "gammalooprs::P",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned()
            ]
        );
        Notation {
            head: Some(IndexDisplay::symbol("p").expect("the canonical label is valid")),
            ..Notation::default()
        }
        .register(symbol)
        .expect("the canonical graph symbol has one notation");
        symbol
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::atom::AtomCore;
    use symbolica::printer::{PrintMode, PrintOptions};

    #[test]
    fn momentum_identity_survives_parse_and_canonical_roundtrip() {
        // Parse first: the ordinary initializer must claim the tagged heads,
        // even when a caller has not used any of the owning constructors yet.
        for (name, constructor, label) in [
            ("Q", momentum as fn() -> Symbol, "q"),
            ("K", loop_momentum as fn() -> Symbol, "k"),
            ("P", external_momentum as fn() -> Symbol, "p"),
        ] {
            let expression = Atom::parse(
                format!("gammalooprs::{name}(27,spenso::cind(0))"),
                "gammalooprs",
                Default::default(),
            )
            .unwrap();
            let head = constructor();
            assert_eq!(expression.as_view().get_symbol(), Some(head));
            assert_eq!(constructor(), head);
            assert_eq!(head.get_namespace(), "gammalooprs");
            assert_eq!(head.get_stripped_name(), name);
            assert!(head.get_attributes().is_empty());
            assert_eq!(head.get_tags().len(), 3);
            assert!(head.has_tag(&SPENSO_TAG.rank1));
            assert!(head.has_tag(&SPENSO_TAG.tensor));
            assert!(head.has_tag(COMPONENT_LABELS_SUBSCRIPT_TAG));
            let canonical = expression.to_canonical_string();
            let reparsed = Atom::parse(&canonical, "gammalooprs", Default::default()).unwrap();
            assert_eq!(reparsed, expression);
            assert_eq!(reparsed.to_canonical_string(), canonical);
            let bytes = expression.as_view().get_data().to_vec();
            let mathml = expression
                .printer(PrintOptions {
                    mode: PrintMode::MathML,
                    ..PrintOptions::default()
                })
                .to_string();
            assert!(mathml.contains(&format!(">{label}</mi>")), "{mathml}");
            assert_eq!(expression.as_view().get_data(), bytes);
        }
    }

    #[test]
    fn momentum_components_keep_object_indices() {
        for (head, label) in [
            (momentum(), "q"),
            (loop_momentum(), "k"),
            (external_momentum(), "p"),
        ] {
            assert_eq!(
                head.call(27).printer(PrintOptions::typst()).to_string(),
                format!("{label}_27")
            );
            let component = head.call_args([Atom::num(27), AIND_SYMBOLS.cind.call(0)]);
            assert_eq!(
                component.printer(PrintOptions::typst()).to_string(),
                format!("attach(attach({label},b:(27)),t:0)")
            );
            let plain = PrintOptions {
                custom_print_mode: SpensoPrintSettings::typst().into(),
                ..Default::default()
            };
            assert_eq!(
                component.printer(plain).to_string(),
                format!("{label}₂₇^(0)")
            );
            let latex = PrintOptions {
                custom_print_mode: SpensoPrintSettings::typst().into(),
                ..PrintOptions::latex()
            };
            assert_eq!(
                component.printer(latex).to_string(),
                format!("{label}_{{27}}^{{0}}")
            );
            assert!(
                component
                    .pow(2)
                    .printer(PrintOptions::typst())
                    .to_string()
                    .contains("lr((attach(")
            );
        }
    }
}
