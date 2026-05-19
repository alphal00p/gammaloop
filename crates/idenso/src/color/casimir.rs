use spenso::shadowing::TensorCollectExt;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
};

use super::{CS, ColorCasimirSettings};

pub(super) fn color_casimir_basis_impl(
    expression: AtomView<'_>,
    settings: ColorCasimirSettings,
) -> Atom {
    ColorCasimirRewriter { settings }.run(expression)
}

struct ColorCasimirRewriter {
    settings: ColorCasimirSettings,
}

impl ColorCasimirRewriter {
    fn run(&self, expression: AtomView<'_>) -> Atom {
        expression
            .to_owned()
            .replace_map(|arg, _context, out| {
                if let Some(replacement) = self.rewrite_node(arg) {
                    **out = replacement;
                }
            })
            .collect_tensors()
    }

    fn rewrite_node(&self, arg: AtomView<'_>) -> Option<Atom> {
        match arg {
            AtomView::Var(var) => self.rewrite_symbol(var.get_symbol()),
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                self.rewrite_symbol_power(base, exponent)
            }
            _ => None,
        }
    }

    fn rewrite_symbol(&self, sym: Symbol) -> Option<Atom> {
        if self.settings.rewrite_na && Self::is_na_symbol(sym) {
            return Some(Self::adjoint_dimension_in_casimirs());
        }
        if self.settings.rewrite_nc && Self::is_nc_symbol(sym) {
            return Some(Atom::var(CS.ca));
        }
        if self.settings.substitute_tr && Self::is_tr_symbol(sym) {
            return Some(Atom::num(1) / Atom::num(2));
        }
        None
    }

    fn rewrite_symbol_power(&self, base: AtomView<'_>, exponent: AtomView<'_>) -> Option<Atom> {
        let exponent = Self::integer_exponent(exponent)?;
        if self.settings.rewrite_na && Self::is_symbol_var(base, Self::is_na_symbol) {
            return Some(Self::atom_integral_power(
                Self::adjoint_dimension_in_casimirs(),
                exponent,
            ));
        }
        if self.settings.rewrite_nc && Self::is_symbol_var(base, Self::is_nc_symbol) {
            return Some(Self::nc_power_in_casimir_basis(exponent));
        }
        None
    }

    fn nc_power_in_casimir_basis(exponent: i64) -> Atom {
        if exponent == 0 {
            return Atom::num(1);
        }
        if exponent < 0 {
            return Self::atom_integral_power(
                Atom::var(CS.ca) - Atom::num(2) * Atom::var(CS.cf),
                -exponent,
            );
        }

        let nc_squared = Self::adjoint_dimension_in_casimirs() + Atom::num(1);
        let even_part = Self::atom_integral_power(nc_squared, exponent / 2);
        if exponent % 2 == 0 {
            even_part
        } else {
            Atom::var(CS.ca) * even_part
        }
    }

    /// Uses `CA = Nc` and `2 CA CF = Nc^2 - 1`.
    fn adjoint_dimension_in_casimirs() -> Atom {
        Atom::num(2) * Atom::var(CS.ca) * Atom::var(CS.cf)
    }

    fn atom_integral_power(base: Atom, exponent: i64) -> Atom {
        match exponent {
            0 => Atom::num(1),
            1 => base,
            _ => base.pow(Atom::num(exponent)),
        }
    }

    fn integer_exponent(expr: AtomView<'_>) -> Option<i64> {
        let AtomView::Num(number) = expr else {
            return None;
        };
        let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
            return None;
        };
        Some(value)
    }

    fn is_symbol_var(arg: AtomView<'_>, predicate: fn(Symbol) -> bool) -> bool {
        matches!(arg, AtomView::Var(var) if predicate(var.get_symbol()))
    }

    fn is_nc_symbol(sym: Symbol) -> bool {
        sym.get_stripped_name() == "Nc"
    }

    fn is_na_symbol(sym: Symbol) -> bool {
        sym.get_stripped_name() == "NA"
    }

    fn is_tr_symbol(sym: Symbol) -> bool {
        sym.get_stripped_name() == "TR"
    }
}
