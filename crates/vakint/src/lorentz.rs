//! Symbolic Lorentz contractions on factorized expressions.

use std::collections::BTreeMap;

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        representation::{BaseRepName, Minkowski},
        slot::SlotMatcher,
    },
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use crate::{VakintError, VakintTerm, symbols::S};

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
            AtomView::Fun(function)
                if function.get_symbol().has_tag(&SPENSO_TAG.tensor)
                    && ![S.k, S.p].contains(&function.get_symbol()) =>
            {
                if view.contains_symbol(S.k) {
                    return Err(VakintError::InvalidNumerator(
                        "loop-dependent native tensor leaves require tensor reduction before scalar evaluation".into(),
                    ));
                }
                let mut indices = BTreeMap::new();
                for argument in function.iter() {
                    if VakintTerm::spenso_minkowski_slot_parts(argument).is_some()
                        || matches!(argument, AtomView::Fun(index) if index.get_symbol() == S.tensor_index)
                    {
                        // Universal projectors replace complete source slots
                        // by their declared marker until the caller restores
                        // native representations. It remains the same port.
                        *indices.entry(argument.to_owned()).or_insert(0_u8) += 1;
                    } else if argument.contains_symbol(Minkowski::selfless_symbol())
                        || argument.contains_symbol(S.tensor_index)
                    {
                        // A compact vector is consumed by a native tensor leaf.
                        // Bindings hidden inside compound chain/trace bodies
                        // require their owning tensor algebra before this scalar
                        // boundary; do not guess that they have no open ports.
                        let compact_vector = matches!(argument, AtomView::Fun(vector)
                            if (vector.get_symbol().has_tag(&SPENSO_TAG.rank1)
                                || vector.get_symbol() == S.p)
                                && vector.iter().any(VakintTerm::is_compact_spenso_minkowski_representation)
                                && !vector.iter().any(|argument| VakintTerm::spenso_minkowski_slot_parts(argument).is_some()));
                        let fixed_basis = SlotMatcher::default()
                            .fixed_basis_vector(argument)
                            .is_some_and(|(_, representation)| {
                                representation.head() == Minkowski::selfless_symbol()
                            });
                        if !compact_vector && !fixed_basis {
                            return Err(VakintError::InvalidNumerator(
                                "compound native tensor interfaces must be materialized before scalar evaluation".into(),
                            ));
                        }
                    }
                }
                if function.get_symbol() == ETS.metric && !indices.is_empty() {
                    return Self::parse(
                        FunctionBuilder::new(S.g)
                            .add_args(function.iter())
                            .finish()
                            .as_view(),
                        dimension,
                    );
                }
                if indices.values().any(|count| *count > 2) {
                    return Err(VakintError::InvalidNumerator(
                        "a native Lorentz index is contracted more than once".into(),
                    ));
                }
                Ok(Self {
                    expression: view.to_owned(),
                    indices,
                })
            }
            AtomView::Fun(function) => {
                let mut result = FunctionBuilder::new(function.get_symbol());
                for argument in function {
                    // A fixed coordinate vector is an operand of an index-free
                    // dot, not a standalone scalar. Its representation must
                    // remain intact until the caller restores native notation.
                    if [S.dot, SPENSO_TAG.dot].contains(&function.get_symbol())
                        && function.get_nargs() == 2
                        && SlotMatcher::default()
                            .fixed_basis_vector(argument)
                            .is_some_and(|(_, representation)| {
                                representation.head() == Minkowski::selfless_symbol()
                            })
                    {
                        if view.contains_symbol(S.k) {
                            return Err(VakintError::InvalidNumerator(
                                "loop-dependent fixed-component dots require tensor reduction before scalar evaluation".into(),
                            ));
                        }
                        result = result.add_arg(argument);
                        continue;
                    }
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
        if [&self.expression, &other.expression]
            .iter()
            .any(|expression| {
                expression.get_all_symbols(true).contains(&S.tensor)
                    || matches!(expression.as_view(), AtomView::Fun(leaf)
                    if leaf.get_symbol().has_tag(&SPENSO_TAG.tensor)
                        && ![S.k, S.p].contains(&leaf.get_symbol()))
            })
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
    fn fixed_coordinate_vectors_are_consumed_only_by_closed_dots() {
        use spenso::structure::{abstract_index::AIND_SYMBOLS, representation::RepName};
        use symbolica::{function, symbol};

        let dimension = Atom::var(symbol!("fixed_coordinate_dimension"));
        let representation = Minkowski {}.to_symbolic([dimension.clone()]);
        let spectator = vk_parse!("(a+b)*(c+d)^11").unwrap();
        for component in 0..4 {
            let basis = function!(
                ETS.delta,
                function!(AIND_SYMBOLS.cind, component),
                &representation
            );
            for dot in [S.dot, SPENSO_TAG.dot] {
                let source = &spectator * function!(dot, &basis, function!(S.p, 3));
                let parsed = LorentzTensor::parse(source.as_view(), &dimension).unwrap();
                assert!(parsed.is_scalar());
                assert_eq!(parsed.expression, source);
                assert_eq!(
                    crate::projection::ScalarTerms::parse(parsed.expression.as_view())
                        .unwrap()
                        .expression(),
                    source
                );
            }
            assert!(LorentzTensor::parse(basis.as_view(), &dimension).is_err());
            let loop_dot = S.dot(&basis, &function!(S.k, 1));
            assert!(LorentzTensor::parse(loop_dot.as_view(), &dimension).is_err());
            let open = S.dot(&basis, &function!(S.p, 3, 2));
            assert!(LorentzTensor::parse(open.as_view(), &dimension).is_err());
        }
        let invalid = function!(ETS.delta, function!(AIND_SYMBOLS.cind, 4), representation);
        assert!(
            LorentzTensor::parse(S.dot(invalid, function!(S.p, 3)).as_view(), &dimension).is_err()
        );
    }

    #[test]
    fn native_projector_markers_remain_lorentz_ports() {
        use spenso::structure::representation::RepName;
        use symbolica::{function, symbol};
        let _ = &*S;
        let dimension = Atom::num(4);
        let mu = Minkowski {}
            .new_rep(4)
            .to_symbolic([Atom::var(symbol!("native_marker_mu"))]);
        let marker = function!(S.tensor_index, symbol!("native_marker_family"), 1, 0);
        let gamma = spenso::tensor_symbol!("native_marker_gamma");
        let source = function!(gamma, &mu) * function!(S.g, &mu, &marker);
        let marked = LorentzTensor::parse(source.as_view(), &dimension).unwrap();
        assert_eq!(marked.expression, function!(gamma, &marker));
        assert!(
            !marked.is_scalar(),
            "the projected gamma retains its declared open Lorentz port"
        );
        let closed = &marked.expression * function!(S.p, 0, &marker) + Atom::one();
        assert!(
            LorentzTensor::parse(closed.as_view(), &dimension)
                .unwrap()
                .is_scalar()
        );
        let hidden = function!(SPENSO_TAG.chain, marked.expression);
        assert!(LorentzTensor::parse(hidden.as_view(), &dimension).is_err());
    }

    #[test]
    fn native_tensor_coefficients_keep_their_connected_lorentz_interface() {
        use spenso::structure::representation::RepName;
        use symbolica::{function, symbol};

        let _ = &*S;
        let dimension_symbol = symbol!("native_coefficient_dimension");
        let dimension = Atom::var(dimension_symbol);
        let mink = Minkowski {}.new_rep(dimension_symbol);
        let mu = mink.to_symbolic([Atom::var(symbol!("native_coefficient_mu"))]);
        let nu = mink.to_symbolic([Atom::var(symbol!("native_loop_nu"))]);
        let gamma = spenso::tensor_symbol!("native_coefficient_gamma");
        let identity = spenso::tensor_symbol!("native_coefficient_identity");
        let spin = function!(gamma, &mu) * function!(S.p, 0, &mu)
            + Atom::var(symbol!("native_coefficient_mass")) * function!(identity);
        let spectator = vk_parse!("(a+b)*(c+d)^11").unwrap();
        let loop_dot = S.dot(function!(S.k, 1), function!(S.p, 1));
        let expected = &spectator * &spin * &loop_dot;
        for scalar_loop in [loop_dot, function!(S.k, 1, &nu) * function!(S.p, 1, &nu)] {
            let input = &spectator * &spin * scalar_loop;
            let parsed = LorentzTensor::parse(input.as_view(), &dimension).unwrap();
            assert!(parsed.is_scalar());
            assert_eq!(parsed.expression, expected);
            let terms = crate::projection::ScalarTerms::parse(parsed.expression.as_view()).unwrap();
            assert_eq!(terms.expression(), expected);
        }
        // A native scalar coefficient in an additive index component must
        // not prevent the attached raw k-p pair from becoming one loop dot.
        let connected = (&spin * function!(S.k, 1, &nu)
            + function!(identity) * function!(S.k, 1, &nu))
            * function!(S.p, 1, &nu);
        let connected = LorentzTensor::parse(connected.as_view(), &dimension).unwrap();
        assert_eq!(
            crate::projection::ScalarTerms::parse(connected.expression.as_view())
                .unwrap()
                .expression(),
            (&spin + function!(identity)) * S.dot(function!(S.k, 1), function!(S.p, 1))
        );
        let compact_slash = function!(gamma, function!(S.p, 0, mink.to_symbolic([])));
        assert!(
            LorentzTensor::parse(compact_slash.as_view(), &dimension)
                .unwrap()
                .is_scalar()
        );
        for invalid in [
            function!(gamma, function!(S.k, 1, mink.to_symbolic([]))),
            function!(SPENSO_TAG.chain, function!(gamma, &mu)),
            function!(gamma, &mu) + function!(gamma, &nu),
        ] {
            assert!(LorentzTensor::parse(invalid.as_view(), &dimension).is_err());
        }
        for unresolved_loop in [
            function!(gamma, &mu) * function!(S.k, 1, &mu),
            function!(S.k, 1, &mu),
            S.dot(function!(S.k, 1), function!(S.k, 1)).pow(-1),
        ] {
            let parsed = LorentzTensor::parse(unresolved_loop.as_view(), &dimension).unwrap();
            assert!(crate::projection::ScalarTerms::parse(parsed.expression.as_view()).is_err());
        }
    }

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
