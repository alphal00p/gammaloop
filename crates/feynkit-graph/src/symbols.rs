//! Canonical graph symbols shared with GammaLoop.
//!
//! Polarization heads and graph-owned index labels retain their meaning in MathML:
//! ```
//! use feynkit_graph::symbols;
//! use spenso::structure::representation::{Minkowski, RepName};
//! use symbolica::atom::AtomCore;
//! use symbolica::printer::{PrintMode, PrintOptions};
//! let index = symbols::hedge_index().call((3, 1));
//! let slot = Minkowski {}.new_rep(4).to_symbolic([index]);
//! for head in [symbols::epsilon(), symbols::epsilonbar()] {
//!     let expression = head.call((7, slot.clone()));
//!     let html = expression.printer(PrintOptions {
//!         mode: PrintMode::MathML, ..PrintOptions::default()
//!     }).to_string();
//!     assert!(html.contains("<mmultiscripts>"));
//! }
//! ```
use spenso::network::tags;
use spenso::network::tags::{COMPONENT_LABELS_SUBSCRIPT_TAG, SPENSO_TAG, TENSOR_PRINT_HEAD_PREFIX};
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::spenso_print_scripted_indexed;
use spenso::structure::representation::IndexDisplay;
use spenso::symbolica_init::SymbolicaInitLazy;
use spenso::utils;
use std::sync::LazyLock;
use symbolica::atom::{Atom, AtomView, Symbol};
use symbolica::printer::PrintUserData;
use symbolica::symbol;

pub fn momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::Q",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned(),
                "spenso::tensor-label:q".to_owned()
            ]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn edge_index() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::edge",
            print = tags::tensor_print,
            tags = [SPENSO_TAG.index.clone(), "spenso::index-label:e".to_owned()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn vertex_index() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::vertex",
            print = tags::tensor_print,
            tags = [SPENSO_TAG.index.clone(), "spenso::index-label:v".to_owned()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

/// Half-edge-owned abstract indices `hedge(owner, local_index)`, or `hedge(owner)`
/// for local zero.
/// See the polarization example in [`crate::symbols`].
pub fn hedge_index() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::hedge",
            print = |a, opt, state| {
                if opt.mode.is_mathml() {
                    return tags::tensor_print(a, opt, state);
                }
                match opt.custom_print_mode.get("spenso") {
                    Some(PrintUserData::Integer(i)) => {
                        let AtomView::Fun(f) = a else {
                            return None;
                        };
                        let SpensoPrintSettings {
                            index_subscripts, ..
                        } = SpensoPrintSettings::from(*i as usize);

                        let mut out = "".to_string();
                        let mut first = true;
                        for arg in f.iter() {
                            let Ok(i) = isize::try_from(arg) else {
                                return None;
                            };

                            if !first {
                                out.push('.');
                            }
                            first = false;
                            if index_subscripts {
                                out.push_str(&utils::to_superscript(i));
                            } else {
                                out.push_str(&utils::to_subscript(i));
                            }
                        }
                        Some(out)
                    }
                    _ => None,
                }
            },
            tags = [SPENSO_TAG.index.clone(), "spenso::index-label:h".to_owned()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn denominator() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::denom"; Scalar;
            der = |_, arg, out| {
                if arg != 3 {
                    **out = Atom::Zero;
                } else {
                    **out = Atom::num(1);
                }
            }
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn dimension() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| symbol!("gammalooprs::dim"));
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn loop_momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::K",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned(),
                "spenso::tensor-label:k".to_owned()
            ]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn external_momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::P",
            print = tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                COMPONENT_LABELS_SUBSCRIPT_TAG.to_owned(),
                "spenso::tensor-label:p".to_owned()
            ]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn u() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::u",
            print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "u") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn ubar() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::ubar",
            print =
                |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "u̅", "overline(u)") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn v() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::v",
            print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "v") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn vbar() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::vbar",
            print =
                |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "v̅", "overline(v)") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

/// Rank-one polarization tensor with a physical momentum ID and Lorentz slot.
/// See the MathML example in [`crate::symbols`].
pub fn epsilon() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::ϵ",
            print = |a, opt, state| {
                if opt.mode.is_mathml() {
                    return tags::tensor_print(a, opt, state);
                }
                spenso_print_scripted_indexed!(a, opt, "ϵ")
            },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

/// Conjugate polarization tensor; the bar is presentation metadata on its head.
/// Source serialization retains the canonical symbol. See [`crate::symbols`].
pub fn epsilonbar() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::ϵbar",
            print = |a, opt, state| {
                if opt.mode.is_mathml() {
                    return tags::tensor_print(a, opt, state);
                }
                spenso_print_scripted_indexed!(
                    a,
                    opt,
                    if opt.mode.is_latex() {
                        r"\bar{\epsilon}"
                    } else {
                        "ϵ̅"
                    },
                    "macron(epsilon.alt)"
                )
            },
            data = IndexDisplay::math(
                "overline",
                vec![IndexDisplay::symbol("ϵ").expect("epsilon is a valid nonempty symbol label")]
            )
            .expect("overline requires the single epsilon label supplied above")
            .symbol_user_data(),
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
                format!("{TENSOR_PRINT_HEAD_PREFIX}typst:macron(epsilon.alt)"),
                format!(r"{TENSOR_PRINT_HEAD_PREFIX}latex:\bar\epsilon")
            ]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::network::tags::mathml;
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::atom::AtomCore;
    use symbolica::parse;
    use symbolica::printer::{PrintMode, PrintOptions, PrintState};

    #[test]
    fn epsilonbar_canonical_roundtrip_preserves_bar_and_open_axis() {
        // Exercise the ordinary global initializer rather than registering a
        // test-only replacement for the production symbol.
        let original = parse!("gammalooprs::ϵbar(7,spenso::mink(4))");
        assert_eq!(original.as_view().get_symbol(), Some(epsilonbar()));
        let bytes = original.as_view().get_data().to_vec();
        let source = original.to_canonical_string();
        let reparsed = Atom::parse(&source, "gammalooprs", Default::default()).unwrap();
        assert_eq!(reparsed, original);
        assert_eq!(reparsed.to_canonical_string(), source);
        assert!(epsilonbar().has_tag(&SPENSO_TAG.rank1));
        assert!(epsilonbar().has_tag(&SPENSO_TAG.tensor));

        let options = PrintOptions {
            mode: PrintMode::MathML,
            hide_all_namespaces: true,
            ..PrintOptions::new()
        };
        let native = reparsed.printer(options.clone()).to_string();
        let html = mathml::render(reparsed.as_view(), &options, &PrintState::new(), &|_, _| {
            None
        })
        .unwrap();
        let exported =
            mathml::expression(reparsed.as_view(), &options, &PrintState::new(), &|_, _| {
                None
            })
            .unwrap();
        for rendered in [&native, &html, &exported] {
            assert!(
                rendered.contains("<mover accent=\"true\"><mi>ϵ</mi><mo>¯</mo></mover>"),
                "{rendered}"
            );
            assert!(rendered.contains("<mmultiscripts>"), "{rendered}");
            assert!(!rendered.contains("ϵbar"), "{rendered}");
        }
        for rendered in [&html, &exported] {
            assert!(
                rendered.contains("data-spenso-open-axis=\"0\""),
                "{rendered}"
            );
            assert!(rendered.contains("Unresolved axis 0"), "{rendered}");
        }
        assert!(
            tags::prepare_tensor_print(&reparsed)
                .printer(PrintOptions::typst())
                .to_string()
                .contains("macron(epsilon.alt)")
        );
        let latex = PrintOptions {
            custom_print_mode: SpensoPrintSettings::typst().into(),
            ..PrintOptions::latex()
        };
        assert!(
            reparsed
                .printer(latex.clone())
                .to_string()
                .contains(r"\bar{\epsilon}")
        );
        assert!(
            tags::prepare_tensor_print(&reparsed)
                .printer(latex)
                .to_string()
                .contains(r"\bar\epsilon")
        );
        assert_eq!(original.as_view().get_data(), bytes);
        assert_eq!(reparsed.as_view().get_data(), bytes);
        assert_eq!(reparsed.to_canonical_string(), source);
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
