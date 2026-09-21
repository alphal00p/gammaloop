//! Canonical GammaLoop energy symbols, also used by standalone FeynKit.

use std::sync::LazyLock;

use spenso::{structure::abstract_index::AIND_SYMBOLS, utils::to_subscript};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    printer::{PrintState, PrintUserData},
    symbol,
};

use crate::EdgeId;

macro_rules! spenso_print_simple_indexed {
    ($a:ident, $opt:ident, $symbol:expr) => {
        spenso_print_simple_indexed!($a, $opt, $symbol, $symbol)
    };
    ($a:ident, $opt:ident, $symbol:expr, $typst_symbol:expr) => {{
        match $opt.custom_print_mode.get("spenso") {
            Some(PrintUserData::Integer(_)) => {
                let AtomView::Fun(f) = $a else {
                    return None;
                };

                let mut out = $symbol.to_string();
                let mut args = f.iter();

                let id = args.next().unwrap();
                let Ok(i) = usize::try_from(id) else {
                    return None;
                };

                if $opt.mode.is_typst() {
                    out = $typst_symbol.to_string();
                    out.push('_');
                    out.push_str(&i.to_string());
                } else {
                    out.push_str(&to_subscript(i as isize));
                }
                let mut first = true;
                for arg in args {
                    if first {
                        first = false;
                        out.push('(');
                    } else {
                        out.push(',');
                    }
                    arg.format(&mut out, $opt, PrintState::new()).unwrap();
                }
                if !first {
                    out.push(')');
                }
                Some(out)
            }
            _ => None,
        }
    }};
}

static ON_SHELL: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "gammalooprs::OSE"; Scalar;
        print = |a, opt, _state| {
            spenso_print_simple_indexed!(a, opt, "Eᵒˢ", r#"E^("os")"#)
        },
        der = |_, arg, out| {
            if arg == 1 { **out = Atom::num(1); }
        }
    )
});
static ENERGY: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "gammalooprs::E",
        print = |a, opt, _state| { spenso_print_simple_indexed!(a, opt, "E") }
    )
});

pub fn on_shell() -> Symbol {
    *ON_SHELL
}
pub fn energy() -> Symbol {
    *ENERGY
}
pub fn on_shell_atom(edge: EdgeId) -> Atom {
    function!(on_shell(), edge.index() as i64)
}
pub fn energy_atom(edge: EdgeId) -> Atom {
    function!(energy(), edge.index() as i64)
}
pub fn external_energy_atom(edge: EdgeId) -> Atom {
    function!(
        feynkit_graph::momentum_symbol(),
        edge.index() as i64,
        function!(AIND_SYMBOLS.cind, 0)
    )
}

// Register attributes before an expression parser can create bare energy heads.
symbolica::initialize!(|| {
    on_shell();
    energy();
});
