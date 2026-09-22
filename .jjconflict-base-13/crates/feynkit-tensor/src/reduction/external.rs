//! Covariant reduction by separating the external span from its orthogonal
//! complement. Only the complement is rotationally averaged. Scalar products
//! with basis momenta remain numerator invariants for subsequent IBP reduction.

use super::{
    Outside, TensorMonomial, TensorReducer, TensorReductionError, TensorReductionTerm, dot,
    indexed_vector, outside_pair,
};
use std::collections::HashMap;
use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};

pub(super) struct ExternalProjector {
    basis: Vec<Atom>,
    inverse: Vec<Atom>,
}

impl ExternalProjector {
    fn parallel_pair(&self, left: &Outside, right: &Outside) -> Atom {
        let mut result = Atom::Zero;
        for (i, p) in self.basis.iter().enumerate() {
            for (j, q) in self.basis.iter().enumerate() {
                result += outside_pair(left, &Outside::Vector(p.clone()))
                    * &self.inverse[i * self.basis.len() + j]
                    * outside_pair(&Outside::Vector(q.clone()), right);
            }
        }
        result.together()
    }
}

impl TensorReducer {
    pub(super) fn external_projector(
        &self,
    ) -> Result<Option<ExternalProjector>, TensorReductionError> {
        if self.external_vectors.is_empty() {
            return Ok(None);
        }
        let basis = self.external_vectors.iter().cloned().collect::<Vec<_>>();
        if let Ok(dimension) = usize::try_from(self.dimension.as_view())
            && basis.len() > dimension
        {
            return Err(TensorReductionError::ExternalGram(
                "more basis vectors than Lorentz dimensions".into(),
            ));
        }
        let probe_index = symbol!("FeynKit::ExternalBasisProbeIndex").to_atom();
        for vector in &basis {
            let indexed = indexed_vector(vector, &probe_index);
            let valid = self
                .indexed_vector(indexed.as_view())?
                .is_some_and(|parsed| parsed.compact == *vector);
            if !valid {
                return Err(TensorReductionError::InvalidExternalVector(vector.clone()));
            }
        }
        let variables = (0..basis.len())
            .map(|i| symbol!(&format!("FeynKit::GramCoefficient{i}")).to_atom())
            .collect::<Vec<_>>();
        let equations = basis
            .iter()
            .map(|p| {
                basis
                    .iter()
                    .zip(&variables)
                    .map(|(q, x)| dot(p, q) * x)
                    .sum::<Atom>()
            })
            .collect::<Vec<_>>();
        let (matrix, _) = Atom::system_to_matrix::<u16, _, _>(&equations, &variables)
            .map_err(|error| TensorReductionError::ExternalGram(error.to_string()))?;
        let inverse = matrix
            .inv()
            .map_err(|error| TensorReductionError::ExternalGram(error.to_string()))?
            .into_vec()
            .into_iter()
            .map(|value| value.to_expression())
            .collect();
        Ok(Some(ExternalProjector { basis, inverse }))
    }

    pub(super) fn reduce_external_monomial(
        &self,
        monomial: TensorMonomial,
        projector: &ExternalProjector,
    ) -> Result<Vec<TensorReductionTerm>, TensorReductionError> {
        if monomial.integrated.is_empty() {
            return self.reduce_monomial(monomial);
        }
        let mut transverse_pairs = HashMap::new();
        for vectors in [
            monomial
                .integrated
                .iter()
                .cloned()
                .map(Outside::Vector)
                .collect::<Vec<_>>(),
            monomial.outside.clone(),
        ] {
            for (i, left) in vectors.iter().enumerate() {
                for right in &vectors[i..] {
                    let pair = outside_pair(left, right);
                    let transverse = (&pair - projector.parallel_pair(left, right)).together();
                    transverse_pairs.insert(pair, transverse);
                }
            }
        }
        // Keep longitudinal factors unexpanded. The same configured budget
        // limits intermediate assignments and the final projector terms.
        let transverse_dimension = &self.dimension - Atom::num(projector.basis.len());
        let mut assignments = vec![(Atom::one(), Vec::new(), Vec::new())];
        for (integrated, outside) in monomial.integrated.iter().zip(&monomial.outside) {
            let parallel = projector.parallel_pair(&Outside::Vector(integrated.clone()), outside);
            let transverse_is_zero = transverse_dimension.is_zero()
                || matches!(outside, Outside::Vector(vector) if self.external_vectors.contains(vector));
            let mut next = Vec::new();
            for (factor, integrated_terms, outside_terms) in assignments {
                if !parallel.is_zero() {
                    next.push((
                        &factor * &parallel,
                        integrated_terms.clone(),
                        outside_terms.clone(),
                    ));
                }
                if !transverse_is_zero {
                    let mut integrated_terms = integrated_terms;
                    let mut outside_terms = outside_terms;
                    integrated_terms.push(integrated.clone());
                    outside_terms.push(outside.clone());
                    next.push((factor, integrated_terms, outside_terms));
                }
                if next.len() > self.output_term_limit {
                    return Err(TensorReductionError::OutputLimit {
                        terms: next.len(),
                        limit: self.output_term_limit,
                    });
                }
            }
            assignments = next;
        }
        let mut transverse_reducer = self.clone();
        transverse_reducer.dimension = transverse_dimension;
        let mut result = Vec::new();
        for (parallel, integrated, outside) in assignments {
            if integrated.len() % 2 == 1 {
                continue;
            }
            let terms = transverse_reducer.reduce_monomial(TensorMonomial {
                scalar: Atom::one(),
                integrated,
                outside,
            })?;
            for term in terms {
                let transverse = term.tensor.replace_map(|view, _, output| {
                    if let Some(value) = transverse_pairs.get(&view.to_owned()) {
                        **output = value.clone();
                    }
                });
                let tensor = &monomial.scalar * &parallel * transverse;
                if !tensor.is_zero() {
                    result.push(TensorReductionTerm {
                        coefficient: term.coefficient,
                        tensor,
                        // The invariants now contain longitudinal subtractions,
                        // so vacuum contraction-orbit metadata does not apply.
                        integrated_orbit: None,
                        projector_orbit: None,
                    });
                }
                if result.len() > self.output_term_limit {
                    return Err(TensorReductionError::OutputLimit {
                        terms: result.len(),
                        limit: self.output_term_limit,
                    });
                }
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::network::library::symbolic::ETS;
    use spenso::structure::representation::{Minkowski, RepName};
    use symbolica::function;

    #[test]
    fn one_external_direction_rank_one_and_mixed_rank_two() {
        let d = symbol!("external_test::D");
        let rep = Minkowski {}.new_rep(d);
        let [k, q, p] = ["external_test::k", "external_test::q", "external_test::p"]
            .map(|name| rep.vector(symbol!(name).to_atom().as_view(), []));
        let [mu, nu] = [
            symbol!("external_test::mu").to_atom(),
            symbol!("external_test::nu").to_atom(),
        ];
        let [km, qn, pm, pn] = [
            indexed_vector(&k, &mu),
            indexed_vector(&q, &nu),
            indexed_vector(&p, &mu),
            indexed_vector(&p, &nu),
        ];
        let reducer = TensorReducer::new(d.to_atom())
            .with_integrated_vector(k.clone())
            .with_integrated_vector(q.clone())
            .with_external_vector(p.clone());
        let rank_one = reducer.reduce(km.as_view()).unwrap().into_expression();
        assert!(
            (rank_one - dot(&k, &p) / dot(&p, &p) * &pm)
                .together()
                .is_zero()
        );
        let rank_two = reducer.reduce((&km * &qn).as_view()).unwrap();
        assert!(!rank_two.is_fully_contracted());
        let transverse_metric =
            function!(ETS.metric, rep.pattern(mu), rep.pattern(nu)) - &pm * &pn / dot(&p, &p);
        let expected = (dot(&k, &q) - dot(&k, &p) * dot(&q, &p) / dot(&p, &p)) / (d.to_atom() - 1)
            * transverse_metric
            + dot(&k, &p) * dot(&q, &p) / dot(&p, &p).pow(2) * pm * pn;
        assert!((rank_two.into_expression() - expected).together().is_zero());
    }

    #[test]
    fn rank_three_reuses_transverse_vacuum_moments() {
        let d = symbol!("external_rank3::D");
        let rep = Minkowski {}.new_rep(d);
        let [k, p, r] = [
            "external_rank3::k",
            "external_rank3::p",
            "external_rank3::r",
        ]
        .map(|name| rep.vector(symbol!(name).to_atom().as_view(), []));
        let reducer = TensorReducer::new(d.to_atom())
            .with_integrated_vector(k.clone())
            .with_external_vector(p.clone());
        let mut numerator = Atom::one();
        for name in [
            "external_rank3::mu",
            "external_rank3::nu",
            "external_rank3::rho",
        ] {
            let index = symbol!(name).to_atom();
            numerator *= indexed_vector(&k, &index) * indexed_vector(&r, &index);
        }
        let longitudinal = dot(&k, &p) * dot(&r, &p) / dot(&p, &p);
        let transverse = (dot(&k, &k) - dot(&k, &p).pow(2) / dot(&p, &p))
            * (dot(&r, &r) - dot(&r, &p).pow(2) / dot(&p, &p));
        let expected = longitudinal.clone().pow(3)
            + Atom::num(3) * longitudinal * transverse / (d.to_atom() - 1);
        let result = reducer.reduce(numerator.as_view()).unwrap();
        assert!(result.is_fully_contracted());
        assert!((result.into_expression() - expected).together().is_zero());
    }

    #[test]
    fn two_direction_gram_inverse_and_full_span() {
        let rep = Minkowski {}.new_rep(2);
        let [k, p, r] = ["external_full::k", "external_full::p", "external_full::r"]
            .map(|name| rep.vector(symbol!(name).to_atom().as_view(), []));
        let mu = symbol!("external_full::mu").to_atom();
        let nu = symbol!("external_full::nu").to_atom();
        let reducer = TensorReducer::new(Atom::num(2))
            .with_integrated_vector(k.clone())
            .with_external_vector(p.clone())
            .with_external_vector(r.clone());
        let determinant = dot(&p, &p) * dot(&r, &r) - dot(&p, &r).pow(2);
        let parallel = |index: &Atom| {
            ((dot(&k, &p) * dot(&r, &r) - dot(&k, &r) * dot(&p, &r)) * indexed_vector(&p, index)
                + (dot(&k, &r) * dot(&p, &p) - dot(&k, &p) * dot(&p, &r))
                    * indexed_vector(&r, index))
                / &determinant
        };
        let numerator = indexed_vector(&k, &mu) * indexed_vector(&k, &nu);
        let result = reducer
            .reduce(numerator.as_view())
            .unwrap()
            .into_expression();
        assert!(
            (result - parallel(&mu) * parallel(&nu))
                .together()
                .is_zero()
        );
        let too_many = reducer.with_external_vector(k.clone());
        assert!(matches!(
            too_many.reduce(numerator.as_view()),
            Err(TensorReductionError::ExternalGram(_))
        ));
    }
}
