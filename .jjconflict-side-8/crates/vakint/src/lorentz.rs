//! Symbolic Lorentz contractions on factorized expressions.

use std::collections::BTreeMap;

use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use crate::{VakintError, symbols::S};

/// Index multiplicities describe the original expression, before contraction.
/// Keeping consumed indices prevents a power or metric trace from hiding an
/// overcontracted index in its enclosing product.
pub(crate) struct LorentzTensor {
    pub(crate) expression: Atom,
    indices: BTreeMap<Atom, u8>,
}

impl LorentzTensor {
    pub(crate) fn parse(view: AtomView<'_>, dimension: &Atom) -> Result<Self, VakintError> {
        match view {
            AtomView::Add(sum) => {
                let terms = sum
                    .iter()
                    .map(|term| Self::parse(term, dimension))
                    .collect::<Result<Vec<_>, _>>()?;
                Self::sum(terms)
            }
            AtomView::Mul(product) => Self::product(
                product
                    .iter()
                    .map(|factor| Self::parse(factor, dimension))
                    .collect::<Result<Vec<_>, _>>()?,
                dimension,
            ),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let base = Self::parse(base, dimension)?;
                let exponent = Self::parse(exponent, dimension)?;
                if !exponent.is_scalar() {
                    return Err(VakintError::InvalidNumerator(format!(
                        "open Lorentz indices inside a scalar exponent: {view}"
                    )));
                }
                if base.is_scalar() {
                    return Ok(Self {
                        expression: base.expression.pow(exponent.expression),
                        indices: base.indices,
                    });
                }
                if !matches!(i64::try_from(exponent.expression.as_view()), Ok(2)) {
                    return Err(VakintError::InvalidNumerator(format!(
                        "an indexed tensor can only have power two: {view}"
                    )));
                }
                // A tensor square contracts the two copies' free slots. Scalar
                // contractions already completed inside each copy stay local.
                let mut copy = Self::parse(base.expression.as_view(), dimension)?;
                copy.indices.retain(|_, count| *count == 1);
                let other = Self {
                    expression: copy.expression.clone(),
                    indices: copy.indices.clone(),
                };
                Self::product(vec![copy, other], dimension)
            }
            AtomView::Fun(function) if function.get_symbol() == S.tensor => {
                if function.get_nargs() < 2 || function.get(0).get_all_symbols(true).contains(&S.k)
                {
                    return Err(VakintError::InvalidNumerator(
                        "tensor(body, slot, ...) requires a loop-independent body and explicit Lorentz slots".into(),
                    ));
                }
                let mut indices = BTreeMap::new();
                for index in function.iter().skip(1) {
                    let count = indices.entry(index.to_owned()).or_insert(0_u8);
                    *count += 1;
                    if *count > 2 {
                        return Err(VakintError::InvalidNumerator(format!(
                            "Lorentz index {index} is contracted more than once"
                        )));
                    }
                }
                Ok(Self {
                    expression: view.to_owned(),
                    indices,
                })
            }
            AtomView::Fun(function) if function.get_symbol() == S.g => {
                if function.get_nargs() != 2 {
                    return Err(VakintError::InvalidNumerator(format!(
                        "a Lorentz metric requires two indices: {view}"
                    )));
                }
                let left = function.get(0).to_owned();
                let right = function.get(1).to_owned();
                let mut indices = BTreeMap::new();
                *indices.entry(left.clone()).or_insert(0) += 1;
                *indices.entry(right.clone()).or_insert(0) += 1;
                Ok(Self {
                    expression: if left == right {
                        dimension.clone()
                    } else {
                        view.to_owned()
                    },
                    indices,
                })
            }
            AtomView::Fun(function)
                if [S.k, S.p].contains(&function.get_symbol()) && function.get_nargs() == 2 =>
            {
                Ok(Self {
                    expression: view.to_owned(),
                    indices: BTreeMap::from([(function.get(1).to_owned(), 1)]),
                })
            }
            AtomView::Fun(function) => {
                let mut result = FunctionBuilder::new(function.get_symbol());
                for argument in function {
                    let argument = Self::parse(argument, dimension)?;
                    if !argument.is_scalar() {
                        return Err(VakintError::InvalidNumerator(format!(
                            "open Lorentz indices inside a scalar function: {view}"
                        )));
                    }
                    result = result.add_arg(argument.expression);
                }
                Ok(Self::scalar(result.finish()))
            }
            _ => Ok(Self::scalar(view.to_owned())),
        }
    }

    fn scalar(expression: Atom) -> Self {
        Self {
            expression,
            indices: BTreeMap::new(),
        }
    }

    pub(crate) fn is_scalar(&self) -> bool {
        self.indices.values().all(|count| *count == 2)
    }

    fn shared_indices(&self, other: &Self) -> Vec<Atom> {
        self.indices
            .iter()
            .filter(|(_, count)| **count == 1)
            .filter(|(index, _)| other.indices.get(*index) == Some(&1))
            .map(|(index, _)| index.clone())
            .collect()
    }

    fn sum(terms: Vec<Self>) -> Result<Self, VakintError> {
        let mut indices = BTreeMap::new();
        let mut external = None;
        let mut expression = Atom::zero();
        for term in terms {
            let free = term
                .indices
                .iter()
                .filter(|(_, count)| **count == 1)
                .map(|(index, _)| index.clone())
                .collect::<Vec<_>>();
            if external.as_ref().is_some_and(|expected| expected != &free) {
                return Err(VakintError::InvalidNumerator(
                    "summands have different open Lorentz indices".into(),
                ));
            }
            external = Some(free);
            indices.extend(term.indices);
            expression += term.expression;
        }
        Ok(Self {
            expression,
            indices,
        })
    }

    fn product(mut factors: Vec<Self>, dimension: &Atom) -> Result<Self, VakintError> {
        let mut indices = BTreeMap::new();
        for factor in &factors {
            for (index, count) in &factor.indices {
                let total = indices.entry(index.clone()).or_insert(0_u8);
                *total += count;
                if *total > 2 {
                    return Err(VakintError::InvalidNumerator(format!(
                        "Lorentz index {index} is contracted more than once"
                    )));
                }
            }
        }

        // Contract only factors in the same connected index component. Scalar
        // spectators and independent closed contractions never enter the
        // bilinear recursion, so their products of sums remain factorized.
        loop {
            let pair = (0..factors.len())
                .flat_map(|left| (left + 1..factors.len()).map(move |right| (left, right)))
                .map(|(left, right)| {
                    (
                        factors[left].shared_indices(&factors[right]).len(),
                        left,
                        right,
                    )
                })
                .filter(|(shared, _, _)| *shared != 0)
                .max();
            let Some((_, left, right)) = pair else {
                break;
            };
            let right = factors.remove(right);
            let left = factors.remove(left);
            factors.push(left.contract(right, dimension)?);
        }
        Ok(Self {
            expression: factors
                .into_iter()
                .fold(Atom::one(), |product, factor| product * factor.expression),
            indices,
        })
    }

    fn contract(self, other: Self, dimension: &Atom) -> Result<Self, VakintError> {
        if self.expression.is_zero() || other.expression.is_zero() {
            let mut indices = BTreeMap::new();
            for tensor in [&self, &other] {
                for (index, count) in &tensor.indices {
                    if *count == 1 {
                        *indices.entry(index.clone()).or_insert(0) += 1;
                    }
                }
            }
            return Ok(Self {
                expression: Atom::zero(),
                indices,
            });
        }
        if self.expression.get_all_symbols(true).contains(&S.tensor)
            || other.expression.get_all_symbols(true).contains(&S.tensor)
        {
            // A declared tensor's algebra belongs to its caller. Only metrics
            // rename its complete Lorentz slots; other contractions retain the
            // product, including any loop momenta needed by vacuum projection.
            for (metric, tensor) in [(&self, &other), (&other, &self)] {
                if let AtomView::Fun(metric) = metric.expression.as_view()
                    && metric.get_symbol() == S.g
                {
                    let shared = self.shared_indices(&other);
                    let replacement = metric
                        .iter()
                        .find(|index| *index != shared[0].as_view())
                        .unwrap();
                    let expression = tensor
                        .expression
                        .replace(shared[0].to_pattern())
                        .with(replacement.to_owned());
                    return Self::parse(expression.as_view(), dimension);
                }
            }
            let mut indices = self.indices;
            for (index, count) in other.indices {
                *indices.entry(index).or_insert(0) += count;
            }
            return Ok(Self {
                expression: self.expression * other.expression,
                indices,
            });
        }
        for (left, right) in [(&self, &other), (&other, &self)] {
            match left.expression.as_view() {
                AtomView::Mul(product) => {
                    let mut factors = product
                        .iter()
                        .map(|factor| Self::parse(factor, dimension))
                        .collect::<Result<Vec<_>, _>>()?;
                    factors.push(Self::parse(right.expression.as_view(), dimension)?);
                    return Self::product(factors, dimension);
                }
                AtomView::Add(sum) => {
                    return Self::sum(
                        sum.iter()
                            .map(|term| {
                                Self::parse(term, dimension)?.contract(
                                    Self::parse(right.expression.as_view(), dimension)?,
                                    dimension,
                                )
                            })
                            .collect::<Result<Vec<_>, _>>()?,
                    );
                }
                _ => {}
            }
        }
        let shared = self.shared_indices(&other);
        let (AtomView::Fun(left), AtomView::Fun(right)) =
            (self.expression.as_view(), other.expression.as_view())
        else {
            return Err(VakintError::InvalidNumerator(format!(
                "cannot contract {} with {}",
                self.expression, other.expression
            )));
        };
        let expression = if left.get_symbol() == S.g && right.get_symbol() == S.g {
            if shared.len() == 2 {
                dimension.clone()
            } else {
                FunctionBuilder::new(S.g)
                    .add_arg(
                        left.iter()
                            .find(|index| *index != shared[0].as_view())
                            .unwrap(),
                    )
                    .add_arg(
                        right
                            .iter()
                            .find(|index| *index != shared[0].as_view())
                            .unwrap(),
                    )
                    .finish()
            }
        } else if left.get_symbol() == S.g || right.get_symbol() == S.g {
            // metric contraction
            let (metric, vector) = if left.get_symbol() == S.g {
                (left, right)
            } else {
                (right, left)
            };
            FunctionBuilder::new(vector.get_symbol())
                .add_arg(vector.get(0))
                .add_arg(
                    metric
                        .iter()
                        .find(|index| *index != shared[0].as_view())
                        .unwrap(),
                )
                .finish()
        } else {
            // dot products
            S.dot(
                FunctionBuilder::new(left.get_symbol())
                    .add_arg(left.get(0))
                    .finish(),
                FunctionBuilder::new(right.get_symbol())
                    .add_arg(right.get(0))
                    .finish(),
            )
        };
        Self::parse(expression.as_view(), dimension)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn opaque_tensor_slots_validate_factorized_production_sums() {
        let _ = &*S;
        let dimension = Atom::num(4);
        let input = vk_parse!(
            "(a+b)*(c+d)*((tensor(gamma(mu),mu)*k(1,mu)-tensor(gamma(nu),nu)*k(1,nu))*k(1,rho)+tensor(gamma(rho),rho)*dot(k(1),k(1)))"
        ).unwrap();
        let parsed = LorentzTensor::parse(input.as_view(), &dimension).unwrap();
        assert_eq!(parsed.expression, input);
        assert!(!parsed.is_scalar());
        let input = vk_parse!("tensor(gamma(mu),mu)*(k(1,mu)+p(1,mu))").unwrap();
        let parsed = LorentzTensor::parse(input.as_view(), &dimension).unwrap();
        assert_eq!(parsed.expression, input);
        assert!(parsed.is_scalar());
        let invalid = vk_parse!("tensor(gamma(mu),mu)+tensor(gamma(nu),nu)").unwrap();
        assert!(LorentzTensor::parse(invalid.as_view(), &dimension).is_err());
        let hidden_loop = vk_parse!("tensor(gamma(k(1,mu)),mu)").unwrap();
        assert!(LorentzTensor::parse(hidden_loop.as_view(), &dimension).is_err());
    }

    #[test]
    fn metrics_rename_only_declared_complete_lorentz_slots() {
        let _ = &*S;
        let input = vk_parse!(
            "g(mink(D,a),mink(D,c))*tensor(gamma(bis(4,a),bis(4,b),mink(D,a)),mink(D,a))"
        )
        .unwrap();
        let expected = vk_parse!("tensor(gamma(bis(4,a),bis(4,b),mink(D,c)),mink(D,c))").unwrap();
        let parsed = LorentzTensor::parse(input.as_view(), &Atom::num(4)).unwrap();
        assert_eq!(parsed.expression, expected);
        // The bispinor dummy named a is independent of the Minkowski slot a.
        assert_eq!(
            LorentzTensor::parse(parsed.expression.as_view(), &Atom::num(4))
                .unwrap()
                .expression,
            expected,
        );
    }
}
