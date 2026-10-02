use super::*;
use std::{
    borrow::Borrow,
    fmt,
    sync::{Arc, Mutex},
};
use symbolica::{
    domains::{InternalOrdering, RingOps, Set, atom::AtomField},
    printer::{PrintOptions, PrintState},
};

use symbolica::domains::float::{Float, FloatField};

#[test]
fn compound_variables_preserve_rounded_coefficient_multiplication_order() {
    let field = FloatField::from_rep(Float::with_val(53, 0));
    let x = Atom::var(symbol!("emission_boundary::x"));
    let y = Atom::var(symbol!("emission_boundary::y"));
    for [a, b, c] in [[0.1, 0.2, 0.3], [1.0e16, 1.0e-16, 0.3], [1.1, 1.3, 1.7]] {
        let mut variables = [
            Atom::num(Float::with_val(53, a)) * &x,
            Atom::num(Float::with_val(53, b)) * &y,
        ];
        variables.sort();
        let coefficient = Float::with_val(53, c);
        // The general emitter appends its coefficient after variable factors.
        // Moving it first changes rounded arithmetic, even for one monomial.
        let expected = (&variables[0] * &variables[1]) * Atom::num(coefficient.clone());
        let polynomial = MultivariatePolynomial::<_, u8>::from_coefficient_list(
            vec![coefficient],
            vec![1, 1],
            vec![
                PolyVariable::Function(symbol!("emission_boundary::map_a"), variables[0].clone()),
                PolyVariable::Function(symbol!("emission_boundary::map_b"), variables[1].clone()),
            ]
            .into(),
            &field,
        );
        assert_eq!(polynomial.to_expression(), expected, "{a} * {b} * {c}");
    }
}

// The public conversion trait permits expression-valued coefficients. This
// test domain preserves raw callback-bearing atoms until the emitter normalizes.
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Expressions;
impl fmt::Display for Expressions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Expressions")
    }
}
#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct Expression(Atom);
impl InternalOrdering for Expression {
    fn internal_cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}
impl Set for Expressions {
    type Element = Expression;
    fn size(&self) -> Option<Integer> {
        None
    }
}
impl<T: Borrow<Expression>> RingOps<T> for Expressions {
    fn add(&self, a: T, b: T) -> Expression {
        Expression(&a.borrow().0 + &b.borrow().0)
    }
    fn sub(&self, a: T, b: T) -> Expression {
        Expression(&a.borrow().0 - &b.borrow().0)
    }
    fn mul(&self, a: T, b: T) -> Expression {
        Expression(&a.borrow().0 * &b.borrow().0)
    }
    fn neg(&self, a: T) -> Expression {
        Expression(-&a.borrow().0)
    }
    fn add_assign(&self, a: &mut Expression, b: T) {
        a.0 = &a.0 + &b.borrow().0;
    }
    fn sub_assign(&self, a: &mut Expression, b: T) {
        a.0 = &a.0 - &b.borrow().0;
    }
    fn mul_assign(&self, a: &mut Expression, b: T) {
        a.0 = &a.0 * &b.borrow().0;
    }
    fn add_mul_assign(&self, a: &mut Expression, b: T, c: T) {
        a.0 = &a.0 + &b.borrow().0 * &c.borrow().0;
    }
    fn sub_mul_assign(&self, a: &mut Expression, b: T, c: T) {
        a.0 = &a.0 - &b.borrow().0 * &c.borrow().0;
    }
}
impl Ring for Expressions {
    fn zero(&self) -> Expression {
        Expression(Atom::num(0))
    }
    fn one(&self) -> Expression {
        Expression(Atom::num(1))
    }
    fn nth(&self, n: Integer) -> Expression {
        Expression(Atom::num(n))
    }
    fn pow(&self, b: &Expression, e: u64) -> Expression {
        Expression(b.0.pow(e))
    }
    fn is_zero(&self, a: &Expression) -> bool {
        a.0.is_zero()
    }
    fn is_one(&self, a: &Expression) -> bool {
        a.0.is_one()
    }
    fn one_is_gcd_unit() -> bool {
        true
    }
    fn characteristic(&self) -> Integer {
        0.into()
    }
    fn try_inv(&self, a: &Expression) -> Option<Expression> {
        AtomField::default().try_inv(&a.0).map(Expression)
    }
    fn try_div(&self, a: &Expression, b: &Expression) -> Option<Expression> {
        AtomField::default().try_div(&a.0, &b.0).map(Expression)
    }
    fn format<W: fmt::Write>(
        &self,
        a: &Expression,
        o: &PrintOptions,
        s: PrintState,
        w: &mut W,
    ) -> Result<bool, fmt::Error> {
        a.0.as_view().format(w, o, s)
    }
}
impl CoefficientToExpression<Expressions> for Expression {
    fn coefficient_to_expression(&self, _: &Expressions, out: &mut Atom) {
        *out = self.0.clone();
    }
}
#[test]
fn general_emission_preserves_variable_and_coefficient_callback_order() {
    let events = Arc::new(Mutex::new(Vec::<&'static str>::new()));
    let v_events = events.clone();
    let variable = symbol!(
        "conversion_callbacks::variable",
        norm = move |_, out| {
            v_events.lock().unwrap().push("variable");
            **out = Atom::num(2);
        }
    );
    let c_events = events.clone();
    let coefficient = symbol!(
        "conversion_callbacks::coefficient",
        norm = move |_, out| {
            c_events.lock().unwrap().push("coefficient");
            **out = Atom::num(3);
        }
    );
    let mut v = Atom::new();
    v.to_fun(variable);
    let mut c = Atom::new();
    c.to_fun(coefficient);
    let p = MultivariatePolynomial::<_, u8>::from_coefficient_list(
        vec![Expression(c)],
        vec![1],
        vec![PolyVariable::Function(variable, v)].into(),
        &Expressions,
    );
    events.lock().unwrap().clear();
    let out = p.to_expression();
    assert_eq!(out, Atom::num(6));
    assert_eq!(*events.lock().unwrap(), ["variable", "coefficient"]);
}
