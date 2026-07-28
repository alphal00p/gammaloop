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

pub struct CofDimensionInvariantRewriter;

impl CofDimensionInvariantRewriter {
    pub(crate) fn run(&self, expression: AtomView<'_>) -> Atom {
        expression.to_owned().replace_map(|arg, _context, out| {
            if let Some(replacement) = self.rewrite_node(arg) {
                **out = replacement;
            }
        })
    }

    fn rewrite_node(&self, arg: AtomView<'_>) -> Option<Atom> {
        if let AtomView::Var(var) = arg {
            return self.rewrite_legacy_symbol(var.get_symbol());
        }

        if let AtomView::Pow(pow) = arg {
            let (base, exponent) = pow.get_base_exp();
            return self.rewrite_legacy_symbol_power(base, exponent);
        }

        let AtomView::Fun(f) = arg else {
            return None;
        };

        if f.get_symbol() == CS.idx && f.get_nargs() == 2 {
            let args = f.iter().collect::<Vec<_>>();
            return self.rewrite_index(args[0], args[1]);
        }

        if f.get_symbol() == CS.cas && f.get_nargs() == 2 {
            let args = f.iter().collect::<Vec<_>>();
            return self.rewrite_casimir(args[0], args[1]);
        }

        if f.get_symbol() == CS.gram && f.get_nargs() == 3 {
            let args = f.iter().collect::<Vec<_>>();
            return self.rewrite_gram(args[0], args[1], args[2]);
        }

        None
    }

    fn rewrite_legacy_symbol(&self, symbol: Symbol) -> Option<Atom> {
        if symbol == CS.tr {
            return Some(Atom::num(1) / Atom::num(2));
        }
        if symbol == CS.ca {
            return Some(Atom::var(CS.nc));
        }
        if symbol == CS.cf {
            return Some(Self::fundamental_quadratic_casimir(Atom::var(CS.nc)));
        }
        if symbol == CS.na {
            return Some(Atom::var(CS.nc).pow(Atom::num(2)) - Atom::num(1));
        }
        None
    }

    fn rewrite_legacy_symbol_power(
        &self,
        base: AtomView<'_>,
        exponent: AtomView<'_>,
    ) -> Option<Atom> {
        let exponent = Self::integer(exponent)?;
        let replacement = match base {
            AtomView::Var(var) => self.rewrite_legacy_symbol(var.get_symbol())?,
            _ => return None,
        };
        Some(Self::atom_integral_power(replacement, exponent))
    }

    fn rewrite_index(&self, degree: AtomView<'_>, rep: AtomView<'_>) -> Option<Atom> {
        if Self::positive_integer(degree)? != 2 {
            return None;
        }

        Self::fundamental_dimension(rep).map(|_| Atom::num(1) / Atom::num(2))
    }

    fn rewrite_casimir(&self, degree: AtomView<'_>, rep: AtomView<'_>) -> Option<Atom> {
        if Self::positive_integer(degree)? != 2 {
            return None;
        }

        if let Some(dimension) = Self::fundamental_dimension(rep) {
            return Some(Self::fundamental_quadratic_casimir(dimension));
        }

        let dimension = Self::adjoint_dimension(rep)?;
        Self::fundamental_dimension_from_adjoint_dimension(dimension.as_view())
    }

    fn rewrite_gram(
        &self,
        degree: AtomView<'_>,
        left_rep: AtomView<'_>,
        right_rep: AtomView<'_>,
    ) -> Option<Atom> {
        let degree = Self::positive_integer(degree)?;
        let left_dimension = Self::fundamental_dimension(left_rep)?;
        let right_dimension = Self::fundamental_dimension(right_rep)?;
        if left_dimension != right_dimension {
            return None;
        }

        match degree {
            3 => Some(Self::fundamental_gram_three(left_dimension)),
            4 => Some(Self::fundamental_gram_four(left_dimension)),
            _ => None,
        }
    }

    fn fundamental_quadratic_casimir(n: Atom) -> Atom {
        (n.clone().pow(Atom::num(2)) - Atom::num(1)) / (Atom::num(2) * n)
    }

    fn fundamental_gram_three(n: Atom) -> Atom {
        let n_squared = n.clone().pow(Atom::num(2));
        (n_squared.clone() - Atom::num(1)) * (n_squared - Atom::num(4)) / (Atom::num(16) * n)
    }

    fn fundamental_gram_four(n: Atom) -> Atom {
        let n_squared = n.clone().pow(Atom::num(2));
        let n_fourth = n_squared.clone().pow(Atom::num(2));
        (n_squared.clone() - Atom::num(1)) * (n_fourth - Atom::num(6) * n_squared + Atom::num(18))
            / (Atom::num(96) * n.pow(Atom::num(2)))
    }

    fn fundamental_dimension(rep: AtomView<'_>) -> Option<Atom> {
        Self::representation_dimension(rep, CS.fundamental_rep)
    }

    fn adjoint_dimension(rep: AtomView<'_>) -> Option<Atom> {
        Self::representation_dimension(rep, CS.adjoint_rep)
    }

    fn representation_dimension(rep: AtomView<'_>, symbol: Symbol) -> Option<Atom> {
        let AtomView::Fun(f) = rep else {
            return None;
        };
        if f.get_symbol() != symbol || f.get_nargs() != 1 {
            return None;
        }

        f.iter().next().map(|dimension| dimension.to_owned())
    }

    fn fundamental_dimension_from_adjoint_dimension(dimension: AtomView<'_>) -> Option<Atom> {
        if Self::is_symbol(dimension, CS.na) {
            return Some(Atom::var(CS.nc));
        }

        let default_adjoint_dimension = Atom::var(CS.nc).pow(Atom::num(2)) - Atom::num(1);
        if dimension == default_adjoint_dimension.as_view() {
            return Some(Atom::var(CS.nc));
        }

        let dimension_plus_one = dimension.to_owned() + Atom::num(1);
        if let Some(root) = Self::square_root_of_square(dimension_plus_one.as_view()) {
            return Some(root);
        }

        Self::integer_square_root_of_plus_one(dimension)
    }

    fn square_root_of_square(expr: AtomView<'_>) -> Option<Atom> {
        let AtomView::Pow(pow) = expr else {
            return None;
        };
        let (base, exponent) = pow.get_base_exp();
        (Self::positive_integer(exponent)? == 2).then(|| base.to_owned())
    }

    fn integer_square_root_of_plus_one(expr: AtomView<'_>) -> Option<Atom> {
        let natural = Self::natural_number(expr)?;
        let root = integer_sqrt(natural + 1)?;
        (root * root == natural + 1).then(|| Atom::num(root))
    }

    fn positive_integer(expr: AtomView<'_>) -> Option<i64> {
        let value = Self::natural_number(expr)?;
        (value > 0).then_some(value)
    }

    fn integer(expr: AtomView<'_>) -> Option<i64> {
        match expr {
            AtomView::Num(number) => match number.get_coeff_view() {
                CoefficientView::Natural(value, 1, 0, 1) => Some(value),
                _ => None,
            },
            _ => None,
        }
    }

    fn natural_number(expr: AtomView<'_>) -> Option<i64> {
        let AtomView::Num(number) = expr else {
            return None;
        };
        let CoefficientView::Natural(value, 1, 0, 1) = number.get_coeff_view() else {
            return None;
        };
        Some(value)
    }

    fn is_symbol(expr: AtomView<'_>, symbol: Symbol) -> bool {
        matches!(expr, AtomView::Var(var) if var.get_symbol() == symbol)
    }

    fn atom_integral_power(base: Atom, exponent: i64) -> Atom {
        match exponent {
            0 => Atom::num(1),
            1 => base,
            _ => base.pow(Atom::num(exponent)),
        }
    }
}

fn integer_sqrt(value: i64) -> Option<i64> {
    if value < 0 {
        return None;
    }
    let mut root = (value as f64).sqrt() as i64;
    while root
        .saturating_add(1)
        .saturating_mul(root.saturating_add(1))
        <= value
    {
        root += 1;
    }
    while root.saturating_mul(root) > value {
        root -= 1;
    }
    Some(root)
}
