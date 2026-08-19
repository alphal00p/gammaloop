use spenso::shadowing::TensorCollectExt;
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::CoefficientView,
};

use super::{CS, ColorCasimirSettings};

pub(super) fn color_casimir_basis_impl(
    expression: AtomView<'_>,
    fundamental_rep: AtomView<'_>,
    adjoint_rep: AtomView<'_>,
    settings: ColorCasimirSettings,
) -> Atom {
    let Some(rewriter) = ColorCasimirRewriter::new(fundamental_rep, adjoint_rep, settings) else {
        return expression.to_owned();
    };
    rewriter.run(expression)
}

struct ColorCasimirRewriter {
    settings: ColorCasimirSettings,
    fundamental_rep: Atom,
    adjoint_rep: Atom,
    fundamental_dimension: Atom,
    adjoint_dimension: Atom,
}

impl ColorCasimirRewriter {
    fn new(
        fundamental_rep: AtomView<'_>,
        adjoint_rep: AtomView<'_>,
        settings: ColorCasimirSettings,
    ) -> Option<Self> {
        let fundamental_dimension =
            Self::representation_dimension(fundamental_rep, CS.fundamental_rep)?;
        let adjoint_dimension = Self::representation_dimension(adjoint_rep, CS.adjoint_rep)?;
        Some(Self {
            settings,
            fundamental_rep: fundamental_rep.to_owned(),
            adjoint_rep: adjoint_rep.to_owned(),
            fundamental_dimension,
            adjoint_dimension,
        })
    }

    fn run(&self, expression: AtomView<'_>) -> Atom {
        expression
            .to_owned()
            .replace_map(|arg, _context, out| {
                if let Some(replacement) = self.rewrite_node(arg) {
                    **out = replacement;
                    return;
                }
                // Function arguments may carry tensor-interface or occurrence metadata.
                // Treat the function as an atomic coefficient and only rewrite an exact
                // invariant function matched above.
                if matches!(arg, AtomView::Fun(_)) {
                    **out = arg.to_owned();
                }
            })
            .collect_tensors()
    }

    fn rewrite_node(&self, arg: AtomView<'_>) -> Option<Atom> {
        if self.settings.substitute_fundamental_index && arg == self.fundamental_index().as_view() {
            return Some(Atom::num(1) / Atom::num(2));
        }

        if Self::is_symbolic_dimension(&self.adjoint_dimension)
            && arg == self.adjoint_dimension.as_view()
        {
            return Some(self.adjoint_dimension_in_casimirs());
        }
        if self.settings.rewrite_fundamental_dimension
            && Self::is_symbolic_dimension(&self.fundamental_dimension)
            && arg == self.fundamental_dimension.as_view()
        {
            return Some(self.adjoint_casimir());
        }

        let AtomView::Pow(power) = arg else {
            return None;
        };
        let (base, exponent) = power.get_base_exp();
        let exponent = Self::integer_exponent(exponent)?;
        if Self::is_symbolic_dimension(&self.adjoint_dimension)
            && base == self.adjoint_dimension.as_view()
        {
            return Some(Self::atom_integral_power(
                self.adjoint_dimension_in_casimirs(),
                exponent,
            ));
        }
        if self.settings.rewrite_fundamental_dimension
            && Self::is_symbolic_dimension(&self.fundamental_dimension)
            && base == self.fundamental_dimension.as_view()
        {
            return Some(Self::atom_integral_power(self.adjoint_casimir(), exponent));
        }
        None
    }

    /// Uses `d_A T_F = d_F C_F`, with the optional SU(N) conventions
    /// applied while constructing the replacement so the pass is idempotent.
    fn adjoint_dimension_in_casimirs(&self) -> Atom {
        let fundamental_dimension = if self.settings.rewrite_fundamental_dimension
            && Self::is_symbolic_dimension(&self.fundamental_dimension)
        {
            self.adjoint_casimir()
        } else {
            self.fundamental_dimension.clone()
        };
        fundamental_dimension * self.fundamental_casimir() / self.normalized_fundamental_index()
    }

    fn fundamental_casimir(&self) -> Atom {
        CS.cas(Atom::num(2), self.fundamental_rep.clone())
    }

    fn adjoint_casimir(&self) -> Atom {
        CS.cas(Atom::num(2), self.adjoint_rep.clone())
    }

    fn fundamental_index(&self) -> Atom {
        CS.idx(Atom::num(2), self.fundamental_rep.clone())
    }

    fn normalized_fundamental_index(&self) -> Atom {
        if self.settings.substitute_fundamental_index {
            Atom::num(1) / Atom::num(2)
        } else {
            self.fundamental_index()
        }
    }

    fn representation_dimension(rep: AtomView<'_>, symbol: Symbol) -> Option<Atom> {
        let AtomView::Fun(function) = rep else {
            return None;
        };
        if function.get_symbol() != symbol || function.get_nargs() != 1 {
            return None;
        }
        function.iter().next().map(|dimension| dimension.to_owned())
    }

    fn is_symbolic_dimension(dimension: &Atom) -> bool {
        !matches!(dimension.as_view(), AtomView::Num(_))
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
