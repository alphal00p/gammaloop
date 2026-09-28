//! Canonical graph symbols shared with GammaLoop.
use spenso::{
    network::tags::SPENSO_TAG,
    shadowing::symbolica_utils::SpensoPrintSettings,
    utils::{to_subscript, to_superscript},
};
use spenso::{spenso_print_scripted_indexed, symbolica_init::SymbolicaInitLazy};
use std::sync::LazyLock;
use symbolica::{
    atom::{Atom, AtomView, Symbol},
    printer::PrintUserData,
    symbol,
};

pub fn momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::Q",
            print = spenso::network::tags::tensor_print,
            tags = [
                SPENSO_TAG.rank1.clone(),
                SPENSO_TAG.tensor.clone(),
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
            print = spenso::network::tags::tensor_print,
            tags = [SPENSO_TAG.index.clone(), "spenso::index-label:e".to_owned()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn vertex_index() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::vertex",
            print = spenso::network::tags::tensor_print,
            tags = [SPENSO_TAG.index.clone(), "spenso::index-label:v".to_owned()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn hedge_index() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::hedge",
            print = |a, opt, _state| {
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
                                out.push_str(&to_superscript(i));
                            } else {
                                out.push_str(&to_subscript(i));
                            }
                        }
                        Some(out)
                    }
                    _ => None,
                }
            },
            tags = [SPENSO_TAG.index.clone()]
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
            print = |a, opt, _state| { spenso::spenso_print_scripted_indexed!(a, opt, "k") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn external_momentum() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::P",
            print = |a, opt, _state| { spenso::spenso_print_scripted_indexed!(a, opt, "p") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
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

pub fn epsilon() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::ϵ",
            print = |a, opt, _state| { spenso_print_scripted_indexed!(a, opt, "ϵ") },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}

pub fn epsilonbar() -> Symbol {
    static SYMBOL: LazyLock<Symbol> = LazyLock::new(|| {
        symbol!(
            "gammalooprs::ϵbar",
            print = |a, opt, _state| {
                spenso_print_scripted_indexed!(a, opt, "ϵ̅", "overline(epsilon.alt)")
            },
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        )
    });
    *SymbolicaInitLazy::new(&SYMBOL)
}
