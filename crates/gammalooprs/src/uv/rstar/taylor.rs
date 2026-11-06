//! Concrete Taylor operations for the graph-first R* preparation.
//!
//! Boundary momenta must already carry the expansion parameter. In particular,
//! a boundary inherited from a parent can be `q + lambda*p`: `q` remains a loop
//! momentum of that parent and must not be scaled here.

use eyre::{Result, WrapErr, ensure};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    id::Replacement,
};

/// A Taylor polynomial of total degree at most `degree` in one common parameter.
///
/// This operation has no pole projection or integration scope. It retains
/// Laurent coefficients and finite pieces that can multiply counterterm poles.
/// Hard-component Taylor expansion uses the same operation at fixed reference
/// mass; only the enclosing compensated expansion calls [`Self::compensate_mass`].
#[derive(Clone, Copy, Debug)]
pub struct TaylorExpansion {
    parameter: Symbol,
    degree: usize,
}

impl TaylorExpansion {
    pub fn new(parameter: Symbol, degree: usize) -> Self {
        Self { parameter, degree }
    }

    /// Replace the marked line's reference mass squared by `(1-lambda^2)*M^2`.
    ///
    /// Apply this to the complete expression, including inherited local
    /// insertions. A contracted marked edge can leave its reference mass in an
    /// insertion, even when no marked propagator remains. The mass parameter is
    /// a symbol representing the squared mass, not the mass itself.
    pub fn compensate_mass(
        &self,
        expression: &Atom,
        reference_mass_squared: Symbol,
        auxiliary_mass_squared: &Atom,
    ) -> Atom {
        let replacement = (Atom::one() - Atom::var(self.parameter).pow(2)) * auxiliary_mass_squared;
        expression.replace(reference_mass_squared).with(replacement)
    }

    /// Scale inherited mass-squared parameters by `lambda^2` in a child task.
    ///
    /// Do this before introducing the child's own auxiliary mass through a
    /// distinct reference parameter. At `lambda=1` the inherited masses are
    /// preserved and the child's auxiliary mass disappears. At zero only that
    /// auxiliary line is massive. Replacements are simultaneous.
    pub fn scale_inherited_masses(&self, expression: &Atom, masses_squared: &[Symbol]) -> Atom {
        let scaling = Atom::var(self.parameter).pow(2);
        let replacements: Vec<_> = masses_squared
            .iter()
            .map(|mass| {
                Replacement::new(Atom::var(*mass).to_pattern(), &scaling * Atom::var(*mass))
            })
            .collect();
        expression.replace_multiple(&replacements)
    }

    /// Expand an explicit expression and retain the common parameter powers.
    ///
    /// Multiply partial polynomials before truncating to the parent's total
    /// degree. Setting the parameter to one in each factor first would retain
    /// spurious cross terms above that degree.
    pub fn polynomial(&self, expression: &Atom) -> Result<Atom> {
        let expansion = expression
            .series(self.parameter, 0, self.degree)
            .wrap_err("could not Taylor-expand the explicit R* expression")?;
        for (power, coefficient) in expansion.terms() {
            ensure!(
                coefficient.is_zero() || (power >= 0 && power.is_integer()),
                "R* Taylor expansion requires an expression regular at zero; found power {power}"
            );
        }
        Ok(expansion.to_atom())
    }

    /// Truncate the complete product or sum, then set its parameter to one.
    pub fn evaluate(&self, expression: &Atom) -> Result<Atom> {
        Ok(self
            .polynomial(expression)?
            .replace(self.parameter)
            .with(Atom::one()))
    }
}

#[cfg(test)]
mod tests {
    use super::TaylorExpansion;
    use symbolica::{
        atom::{Atom, AtomCore},
        parse, symbol,
    };

    fn assert_equal(actual: Atom, expected: Atom) {
        assert!(
            (&actual - &expected).together().is_zero(),
            "actual: {actual}\nexpected: {expected}"
        );
    }

    #[test]
    fn compensated_propagator_includes_mass_compensation_and_actual_power() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 2);
        let denominator = parse!("q2 + 2*lambda*qp + lambda^2*p2 + mu2");
        for power in [1, 2, 3] {
            let marked =
                taylor.compensate_mass(&denominator.pow(-power), symbol!("mu2"), &parse!("M2"));
            let d = parse!("q2 + M2");
            let expected = d.pow(-power) - Atom::num(2 * power) * parse!("qp") * d.pow(-power - 1)
                + Atom::num(power) * parse!("M2 - p2") * d.pow(-power - 1)
                + Atom::num(2 * power * (power + 1)) * parse!("qp^2") * d.pow(-power - 2);
            assert_equal(taylor.evaluate(&marked).unwrap(), expected);
        }
    }

    #[test]
    fn compensated_local_polynomial_loses_only_the_auxiliary_mass() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 4);
        let local = parse!("a*lambda^4*p^4 + b*mu2*lambda^2*p^2 + c*mu2^2");
        let compensated = taylor.compensate_mass(&local, symbol!("mu2"), &parse!("M2"));
        assert_equal(taylor.evaluate(&compensated).unwrap(), parse!("a*p^4"));
    }

    #[test]
    fn child_preserves_inherited_masses_and_removes_its_own_auxiliary_mass() {
        let taylor = TaylorExpansion::new(symbol!("rho"), 2);
        let local = parse!("a*rho^2*p^2 + b*mu2 + c*m2 + d*child2");
        let inherited = taylor.scale_inherited_masses(&local, &[symbol!("mu2"), symbol!("m2")]);
        let compensated = taylor.compensate_mass(&inherited, symbol!("child2"), &parse!("Lambda2"));
        assert_equal(
            compensated.replace(symbol!("rho")).with(Atom::Zero),
            parse!("d*Lambda2"),
        );
        assert_equal(
            taylor.evaluate(&compensated).unwrap(),
            parse!("a*p^2 + b*mu2 + c*m2"),
        );
    }

    #[test]
    fn hard_taylor_keeps_reference_mass_and_remaining_loop_momentum() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 1);
        assert_equal(
            taylor
                .polynomial(&parse!("a*(q+lambda*p)^2 + b*mu2"))
                .unwrap(),
            parse!("a*q^2 + 2*a*lambda*q*p + b*mu2"),
        );
    }

    #[test]
    fn partial_polynomials_are_multiplied_before_total_truncation() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 2);
        let left = taylor.polynomial(&parse!("1+lambda")).unwrap();
        let right = taylor.polynomial(&parse!("1-lambda+lambda^2")).unwrap();
        assert_equal(taylor.evaluate(&(&left * &right)).unwrap(), Atom::one());
        assert_equal(
            taylor.evaluate(&left).unwrap() * taylor.evaluate(&right).unwrap(),
            Atom::num(2),
        );
    }

    #[test]
    fn finite_compensated_tadpole_stays_available_to_multiply_child_poles() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 2);
        // I1 + M^2 I2 = -M^2 + O(epsilon), using the normalization
        // of Eq. (8). Its finite part becomes a pole when multiplied by z/eps.
        let prepared = parse!("z/eps*(I1 + lambda^2*M2*I2)")
            .replace(symbol!("I1"))
            .with(parse!("-M2/eps + M2*(L-1)"))
            .replace(symbol!("I2"))
            .with(parse!("1/eps-L"));
        assert_equal(taylor.evaluate(&prepared).unwrap(), parse!("-z*M2/eps"));
    }

    #[test]
    fn singular_and_fractional_expansions_are_rejected() {
        let taylor = TaylorExpansion::new(symbol!("lambda"), 2);
        assert!(taylor.polynomial(&parse!("1/lambda")).is_err());
        assert!(taylor.polynomial(&parse!("lambda^(1/2)")).is_err());
    }
}
