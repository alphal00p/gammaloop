//! Canonical GammaLoop energy symbols, also used by standalone FeynKit.

use std::sync::LazyLock;

use spenso::network::tags::TENSOR_PRINT_CALLBACK_TAG;
use spenso::shadowing::symbolica_utils::{SpensoPrintBackend, SpensoPrintSettings};
use spenso::structure::abstract_index::AIND_SYMBOLS;
use spenso::symbolica_init::SymbolicaInitLazy;
use spenso::utils;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::printer::{PrintOptions, PrintState};
use symbolica::{function, symbol};

use crate::EdgeId;

// The symbol owns its notation in native text, tensor MathML and graph Typst.
// Canonical/plain serialization deliberately keeps the original function name.
fn print_energy(
    atom: AtomView<'_>,
    options: &PrintOptions,
    state: &PrintState,
    on_shell: bool,
) -> Option<String> {
    let backend = if options.mode.is_latex() {
        SpensoPrintBackend::Latex
    } else {
        SpensoPrintSettings::resolve(options)?.backend
    };
    let AtomView::Fun(function) = atom else {
        return None;
    };
    // CFF uses OSE(edge); GammaLoop's explicit mass-shell marker also
    // carries the radicand in OSE(edge, radicand), before its square root.
    if function.get_nargs() != 1 && (!on_shell || function.get_nargs() != 2) {
        return None;
    }
    let mut arguments = function.iter();
    let index = arguments.next()?;
    let id = index.printer(options.clone()).to_string();
    let radicand = arguments
        .next()
        .map(|argument| argument.printer(options.clone()).to_string());
    Some(match backend {
        SpensoPrintBackend::Typst => {
            let head = if on_shell { "E^(\"os\")" } else { "E" };
            let value = match &radicand {
                Some(radicand) if usize::try_from(index).is_ok() => {
                    format!("{head}_{id}({radicand})")
                }
                Some(radicand) => format!("{head}_({id})({radicand})"),
                None => format!("{head}_({id})"),
            };
            if on_shell && state.in_exp_base {
                format!("lr(({value}))")
            } else {
                value
            }
        }
        SpensoPrintBackend::Latex => {
            let head = if on_shell { r"E^{\mathrm{os}}" } else { "E" };
            let value = match &radicand {
                Some(radicand) => format!(r"{head}_{{{id}}}\left({radicand}\right)"),
                None => format!("{head}_{{{id}}}"),
            };
            if on_shell && state.in_exp_base {
                format!(r"\left({value}\right)")
            } else {
                value
            }
        }
        SpensoPrintBackend::Plain => {
            let head = if on_shell { "Eᵒˢ" } else { "E" };
            let index = usize::try_from(index).ok()?;
            let index = utils::to_subscript(index as isize);
            match &radicand {
                Some(radicand) => format!("{head}{index}({radicand})"),
                None => format!("{head}{index}"),
            }
        }
    })
}

static ON_SHELL: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "gammalooprs::OSE"; Scalar;
        tag = TENSOR_PRINT_CALLBACK_TAG,
        print = |a, opt, state| print_energy(a, opt, state, true),
        der = |_, arg, out| {
            if arg == 1 { **out = Atom::num(1); }
        }
    )
});
static ENERGY: LazyLock<Symbol> = LazyLock::new(|| {
    symbol!(
        "gammalooprs::E",
        tag = TENSOR_PRINT_CALLBACK_TAG,
        print = |a, opt, state| print_energy(a, opt, state, false)
    )
});

pub fn on_shell() -> Symbol {
    *SymbolicaInitLazy::new(&ON_SHELL)
}
pub fn energy() -> Symbol {
    *SymbolicaInitLazy::new(&ENERGY)
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
symbolica::initialize!(|| spenso::symbolica_init::in_symbolica_initializer(|| {
    on_shell();
    energy();
}));

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::parse;

    #[test]
    fn energy_notation_preserves_identity_and_handles_powers() {
        let energy = on_shell_atom(EdgeId(3));
        let canonical = energy.to_canonical_string();
        assert_eq!(energy, parse!("gammalooprs::OSE(3)"));
        assert!(canonical.ends_with("::OSE(3)"));
        assert!(on_shell().has_tag(TENSOR_PRINT_CALLBACK_TAG));
        assert_eq!(
            energy.printer(PrintOptions::typst()).to_string(),
            "E^(\"os\")_(3)"
        );
        assert_eq!(
            energy.printer(PrintOptions::latex()).to_string(),
            r"E^{\mathrm{os}}_{3}"
        );
        let squared = energy
            .clone()
            .pow(2)
            .printer(PrintOptions::typst())
            .to_string();
        assert!(squared.contains("lr((E^(\"os\")_(3)))"), "{squared}");
        assert_eq!(energy.to_canonical_string(), canonical);
        // Unsupported arities keep their ordinary spelling instead of panicking
        // or silently discarding an argument.
        for malformed in [on_shell().call(()), on_shell().call((3, 4, 5))] {
            assert!(
                malformed
                    .printer(PrintOptions::typst())
                    .to_string()
                    .contains("OSE")
            );
        }
    }

    #[test]
    fn mass_shell_radicands_preserve_notation_identity_and_derivative() {
        let energy = on_shell().call((3, 13));
        let canonical = energy.to_canonical_string();
        assert_eq!(energy, parse!("gammalooprs::OSE(3,13)"));
        assert_eq!(parse!(&canonical), energy);
        assert_eq!(
            energy.printer(PrintOptions::typst()).to_string(),
            "E^(\"os\")_3(13)"
        );
        assert_eq!(
            energy.printer(PrintOptions::latex()).to_string(),
            r"E^{\mathrm{os}}_{3}\left(13\right)"
        );
        let mut plain = SpensoPrintSettings::compact().nice_symbolica();
        plain.color_builtin_symbols = false;
        assert_eq!(energy.printer(plain).to_string(), "Eᵒˢ₃(13)");
        assert!(
            energy
                .as_view()
                .pow(2)
                .printer(PrintOptions::typst())
                .to_string()
                .contains("lr((E^(\"os\")_3(13)))")
        );
        assert!(
            energy
                .as_view()
                .pow(2)
                .printer(PrintOptions::latex())
                .to_string()
                .contains(r"\left(E^{\mathrm{os}}_{3}\left(13\right)\right)")
        );
        assert_eq!(energy.to_canonical_string(), canonical);

        let variable = symbol!("feynkit_cff::radicand_variable");
        let radicand = Atom::var(variable).pow(2) + Atom::num(13);
        let marked = on_shell().call((3, radicand.as_view()));
        assert_eq!(
            marked.derivative(variable),
            Atom::num(2) * Atom::var(variable)
        );
        assert_eq!(parse!(marked.to_canonical_string()), marked);
        let radicand_text = radicand.printer(PrintOptions::typst()).to_string();
        assert!(
            marked
                .printer(PrintOptions::typst())
                .to_string()
                .contains(&format!("({radicand_text})"))
        );
    }
}
