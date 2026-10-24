//! The trace recurrence uses the tensor-level exact coefficient collector.

#[cfg(test)]
use symbolica::atom::Atom;
use symbolica::atom::AtomView;
use symbolica::domains::rational::Q;

use super::TraceAlgebra;
use crate::tensor::{CollectedOutput, CollectedPolynomial};

impl TraceAlgebra for &mut CollectedOutput<'_> {
    type Value = CollectedPolynomial;

    fn unit(&mut self, trace_unit: AtomView<'_>) -> CollectedPolynomial {
        self.polynomial(trace_unit, false)
    }

    fn scale(&mut self, coefficient: i32, value: CollectedPolynomial) -> CollectedPolynomial {
        value.mul_coeff(coefficient.into())
    }

    fn metric_product(
        &mut self,
        metric: AtomView<'_>,
        coefficient: i32,
        value: CollectedPolynomial,
    ) -> CollectedPolynomial {
        let metric = self.polynomial(metric, false);
        self.product(metric, value).mul_coeff(coefficient.into())
    }

    fn dimension_product(
        &mut self,
        dimension: AtomView<'_>,
        constant: i32,
        linear: i32,
        value: CollectedPolynomial,
    ) -> CollectedPolynomial {
        let dimension = self.polynomial(dimension, false).mul_coeff(linear.into());
        self.product(
            dimension + CollectedPolynomial::new_one(&Q).mul_coeff(constant.into()),
            value,
        )
    }

    fn add(
        &mut self,
        left: CollectedPolynomial,
        right: CollectedPolynomial,
    ) -> CollectedPolynomial {
        left + right
    }
    fn subtract(
        &mut self,
        left: CollectedPolynomial,
        right: CollectedPolynomial,
    ) -> CollectedPolynomial {
        left - right
    }
    fn sum(&mut self, values: Vec<CollectedPolynomial>) -> CollectedPolynomial {
        values
            .into_iter()
            .fold(CollectedPolynomial::new_zero(&Q), |sum, value| sum + value)
    }

    #[cfg(test)]
    fn finish(self, value: CollectedPolynomial) -> Atom {
        CollectedOutput::finish(self, value).unwrap()
    }
}
