use std::collections::BTreeMap;

use spenso::structure::dimension::Dimension;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore};
use symbolica::graph::{CanonicalForm, Graph};
use symbolica::tensors::matrix::{Matrix, MatrixError};

use super::{IntegralFamily, IntegralFamilyError, PropagatorMapping};

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum ParametricVertex {
    Parameter,
    Term { second: bool, coefficient: Atom },
}

impl IntegralFamily {
    /// Find a real loop direction with an unconstrained transverse integral.
    ///
    /// The returned nonzero vector `w` is in loop-momentum order and satisfies
    /// `A_a w = 0` for every denominator's quadratic loop matrix `A_a`.
    /// Thus `k_i -> k_i + w_i r_perp` leaves all denominators unchanged for
    /// any vector orthogonal to the external span. Its unrestricted polynomial
    /// integral is scaleless in dimensional regularization, even with masses
    /// and eikonal terms in the remaining coordinates.
    ///
    /// Requires a nonsingular external Gram matrix and a nonempty transverse
    /// space. Symbolic dimension is interpreted generically; a concrete
    /// dimension must exceed the number of external basis vectors. `None`
    /// means no certificate was found, not that the integral is nonzero.
    /// Use [`Self::sector`] first to exclude numerator-only propagators.
    /// This verifies a common real null direction, not merely `det(A)=0`.
    pub fn scaleless_transverse_direction(&self) -> Result<Option<Vec<Atom>>, IntegralFamilyError> {
        if matches!(self.kinematics.dimension(), Dimension::Concrete(d) if d <= self.external_momenta.len())
        {
            return Ok(None);
        }
        if !self.external_momenta.is_empty() {
            let rows = self
                .external_momenta
                .iter()
                .map(|p| {
                    self.external_momenta
                        .iter()
                        .map(|q| Ok(self.kinematics.scalar_product(p, q)? * q))
                        .collect::<Result<Vec<Atom>, IntegralFamilyError>>()
                        .map(|row| row.into_iter().sum::<Atom>())
                })
                .collect::<Result<Vec<_>, _>>()?;
            let (gram, _) = Atom::system_to_matrix::<u16, _, _>(&rows, &self.external_momenta)
                .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?;
            if gram
                .det()
                .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?
                .is_zero()
            {
                return Ok(None);
            }
        }
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        let mut equations = Vec::new();
        for denominator in &self.denominators {
            let coefficients = denominator
                .coefficient_list::<i32>(&self.scalar_products)
                .into_iter()
                .collect::<BTreeMap<_, _>>();
            for p in &self.loop_momenta {
                let row = self
                    .loop_momenta
                    .iter()
                    .map(|q| {
                        coefficients
                            .get(&rep.inner_product(p, q))
                            .cloned()
                            .unwrap_or_default()
                            * q
                            * if p == q { 2 } else { 1 }
                    })
                    .sum::<Atom>();
                equations.push(row);
            }
        }
        // Fix one component to one to exclude the trivial kernel vector.
        // Symbolica supplies the linear solve; no separate row reduction lives here.
        for pivot in &self.loop_momenta {
            equations.push(pivot - 1);
            let (matrix, rhs) = Atom::system_to_matrix::<u16, _, _>(&equations, &self.loop_momenta)
                .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?;
            equations.pop();
            match matrix.solve_any(&rhs) {
                Ok(solution) => {
                    let weights = solution
                        .into_vec()
                        .into_iter()
                        .map(|w| w.to_expression())
                        .collect::<Vec<_>>();
                    if weights.iter().all(|w| w.is_real().is_true()) {
                        return Ok(Some(weights));
                    }
                }
                Err(MatrixError::Inconsistent) => {}
                Err(error) => return Err(IntegralFamilyError::InvalidBasis(error.to_string())),
            }
        }
        Ok(None)
    }

    /// Find a parametric scaling certificate for a scaleless sector.
    ///
    /// For `G = U + F`, solve `sum_i w_i*x_i*dG/dx_i = G` over constant
    /// weights using Symbolica. A solution proves the sector scaleless in
    /// dimensional regularization at generic dimension. `None` means this
    /// criterion did not establish scalelessness; it is not a nonzero-integral
    /// assertion. The certificate uses the supplied parameter order.
    ///
    /// Every family denominator is treated as present. Use [`Self::sector`]
    /// first to exclude zero- and negative-power entries. This criterion
    /// requires a nonsingular quadratic form even though [`Self::symanzik`]
    /// can return the algebraic polynomials for singular forms.
    pub fn scaleless_scaling(
        &self,
        parameters: &[Atom],
    ) -> Result<Option<Vec<Atom>>, IntegralFamilyError> {
        let (u, f) = self.symanzik(parameters)?;
        if u.is_zero() {
            return Err(IntegralFamilyError::InvalidBasis(
                "Parametric scaling requires a nonsingular quadratic loop matrix".into(),
            ));
        }
        let polynomial = (u + f).expand().to_polynomial_in_vars::<u32>(parameters);
        let equations = polynomial
            .into_iter()
            .filter(|term| !term.coefficient.together().is_zero())
            .map(|term| {
                term.exponents
                    .iter()
                    .zip(parameters)
                    .map(|(exponent, weight)| Atom::num(*exponent) * weight)
                    .sum::<Atom>()
                    - 1
            })
            .collect::<Vec<_>>();
        // Parameter names serve only as linear-system unknowns after extracting
        // the exponent vectors; the physical coefficients no longer occur.
        let (matrix, rhs) = Atom::system_to_matrix::<u16, _, _>(&equations, parameters)
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?;
        match matrix.solve_any(&rhs) {
            Ok(solution) => Ok(Some(
                solution
                    .into_vec()
                    .into_iter()
                    .map(|weight| weight.to_expression())
                    .collect(),
            )),
            Err(MatrixError::Inconsistent) => Ok(None),
            Err(error) => Err(IntegralFamilyError::InvalidBasis(error.to_string())),
        }
    }

    /// Find a parameter permutation identifying both Symanzik polynomials.
    ///
    /// Parameters are distinct new symbols or calls, one per source propagator.
    /// Both families must have the same propagator count and external kinematics.
    /// Symbolica canonizes the polynomial incidence graphs, preserving U/F,
    /// exact coefficients and monomial exponents. This detects polynomial
    /// equivalences even when no affine loop-momentum map is available.
    ///
    /// This is an algebraic parameter map, not a contour or prescription check.
    /// There is no implied tensor-numerator or external-momentum substitution.
    /// Singular quadratic forms are rejected because their degenerate U/F
    /// polynomials can discard physical masses and external invariants.
    pub fn parametric_mapping(
        &self,
        target: &Self,
        parameters: &[Atom],
    ) -> Result<Option<PropagatorMapping>, IntegralFamilyError> {
        self.compatible_kinematics(target)?;
        self.validate_labels(parameters)?;
        if self.denominators.len() != target.denominators.len() {
            return Ok(None);
        }
        let source = self.canonical_symanzik(parameters)?;
        let destination = target.canonical_symanzik(parameters)?;
        if source.graph != destination.graph {
            return Ok(None);
        }
        let target_indices = destination.vertex_map[..parameters.len()]
            .iter()
            .enumerate()
            .map(|(i, canonical)| (*canonical, i))
            .collect::<BTreeMap<_, _>>();
        Ok(Some(PropagatorMapping {
            denominator_map: source.vertex_map[..parameters.len()]
                .iter()
                .map(|canonical| target_indices[canonical])
                .collect(),
            target_denominator_count: parameters.len(),
        }))
    }

    fn canonical_symanzik(
        &self,
        parameters: &[Atom],
    ) -> Result<CanonicalForm<ParametricVertex, u32>, IntegralFamilyError> {
        let (u, f) = self.symanzik(parameters)?;
        if u.is_zero() {
            return Err(IntegralFamilyError::InvalidBasis(
                "Parametric mapping requires a nonsingular quadratic loop matrix".into(),
            ));
        }
        let mut graph = Graph::new();
        for _ in parameters {
            graph.add_node(ParametricVertex::Parameter);
        }
        for (second, expression) in [(false, u), (true, f)] {
            let polynomial = expression.to_polynomial_in_vars::<u32>(parameters);
            for term in &polynomial {
                let node = graph.add_node(ParametricVertex::Term {
                    second,
                    coefficient: term.coefficient.together(),
                });
                for (parameter, exponent) in term.exponents.iter().enumerate() {
                    if *exponent != 0 {
                        graph
                            .add_edge(parameter, node, false, *exponent)
                            .expect("both polynomial incidence nodes were inserted above");
                    }
                }
            }
        }
        Ok(graph.canonize())
    }

    /// Compute the Symanzik polynomials `(U, F)` in propagator order.
    ///
    /// For `sum_i x_i D_i = k.M.k + 2 k.Q + J`, the convention is
    /// `U = det(M)` and `F = Q.adj(M).Q - U*J`. Thus a Minkowski
    /// denominator `k^2-m^2` gives `(x, m^2*x^2)`. No integration measure,
    /// propagator prescription or powers are inferred. Parameters must be
    /// distinct symbols or calls absent from the family expressions.
    ///
    /// Symbolica evaluates the determinant and cofactors without requiring an
    /// inverse. Linear eikonal denominators and singular quadratic forms are
    /// accepted. A singular result is algebraic data, not a valid Gaussian
    /// integration formula or a scalelessness certificate.
    pub fn symanzik(&self, parameters: &[Atom]) -> Result<(Atom, Atom), IntegralFamilyError> {
        self.validate_labels(parameters)?;
        let rep = Minkowski {}.new_rep(self.kinematics.dimension());
        let weighted = parameters
            .iter()
            .zip(&self.denominators)
            .map(|(x, d)| x * d)
            .sum::<Atom>()
            .expand();
        let coefficients = weighted
            .coefficient_list::<i32>(&self.scalar_products)
            .into_iter()
            .collect::<BTreeMap<_, _>>();
        let coefficient = |p: &Atom, q: &Atom| {
            coefficients
                .get(&rep.inner_product(p, q))
                .cloned()
                .unwrap_or_else(Atom::new)
        };
        let rows = self
            .loop_momenta
            .iter()
            .map(|p| {
                self.loop_momenta
                    .iter()
                    .map(|q| coefficient(p, q) * q / if p == q { 1 } else { 2 })
                    .sum::<Atom>()
            })
            .collect::<Vec<_>>();
        let (matrix, _) = Atom::system_to_matrix::<u16, _, _>(&rows, &self.loop_momenta)
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?;
        let u = matrix
            .det()
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?
            .to_expression();
        let shifts = self
            .loop_momenta
            .iter()
            .map(|p| {
                self.external_momenta
                    .iter()
                    .map(|q| coefficient(p, q) * q / 2)
                    .sum::<Atom>()
            })
            .collect::<Vec<_>>();
        let constant = coefficients.get(&Atom::one()).cloned().unwrap_or_default();
        let mut f = -&u * constant;
        for (i, p) in shifts.iter().enumerate() {
            for (j, q) in shifts[..=i].iter().enumerate() {
                let product = self.kinematics.scalar_product(p, q)?;
                if product.is_zero() {
                    continue;
                }
                // M is symmetric. The (i,j) adjugate entry is the signed
                // minor deleting row j and column i, including det([])=1.
                let minor = Matrix::from_nested_vec(
                    matrix
                        .iter()
                        .as_slice()
                        .chunks(shifts.len())
                        .enumerate()
                        .filter(|(row, _)| *row != j)
                        .map(|(_, row)| {
                            row.iter()
                                .enumerate()
                                .filter(|(column, _)| *column != i)
                                .map(|(_, value)| value.clone())
                                .collect()
                        })
                        .collect(),
                    matrix.field().clone(),
                )
                .map_err(IntegralFamilyError::InvalidBasis)?;
                let cofactor = minor
                    .det()
                    .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?
                    .to_expression();
                let sign = if (i + j) % 2 == 0 { 1 } else { -1 };
                f += cofactor * product * sign * if i == j { 1 } else { 2 };
            }
        }
        let f = f.together().expand();
        Ok((u.expand(), f))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_kinematics::Kinematics;
    use itertools::Itertools;
    use symbolica::parse;

    #[test]
    fn scaleless_bubble_certificates_obey_euler_identity() {
        let k = parse!("scale_bubble::k");
        let p = parse!("scale_bubble::p");
        let variables = [symbolica::symbol!("x"), symbolica::symbol!("y")];
        let parameters = variables.map(|s| s.to_atom());
        for (invariant, mass_squared, scaleless) in [
            (Atom::new(), Atom::new(), true),
            (parse!("s"), Atom::new(), false),
            (Atom::new(), parse!("m2"), false),
        ] {
            let kin = Kinematics::new()
                .with_momenta([k.clone(), p.clone()])
                .unwrap()
                .with_mass_squared(&p, invariant)
                .unwrap();
            let family = IntegralFamily::new(
                vec![k.clone()],
                vec![p.clone()],
                vec![
                    kin.scalar_product(&k, &k).unwrap() - &mass_squared,
                    kin.scalar_product(&(&k - &p), &(&k - &p)).unwrap() - mass_squared,
                ],
                &kin,
            )
            .unwrap();
            let certificate = family.scaleless_scaling(&parameters).unwrap();
            assert_eq!(certificate.is_some(), scaleless);
            if let Some(weights) = certificate {
                let (u, f) = family.symanzik(&parameters).unwrap();
                let g = u + f;
                let euler = weights
                    .iter()
                    .zip(variables)
                    .zip(&parameters)
                    .map(|((w, v), x)| w * x * g.derivative(v))
                    .sum::<Atom>();
                assert!((euler - g).expand().is_zero());
            }
        }
    }

    #[test]
    fn feyncalc_eikonal_scaleless_examples() {
        // FCLoopPakScalelessQ: (2 k.p)^-1 (k^2-m^2)^-1 is scaleless only for m=0.
        let k = parse!("scale_eikonal::k");
        let p = parse!("scale_eikonal::p");
        let kin = Kinematics::new()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let parameters = [parse!("x"), parse!("y")];
        for (mass_squared, scaleless) in [(Atom::new(), true), (parse!("m2"), false)] {
            let family = IntegralFamily::new(
                vec![k.clone()],
                vec![p.clone()],
                vec![
                    kin.scalar_product(&k, &p).unwrap() * 2,
                    kin.scalar_product(&k, &k).unwrap() - mass_squared,
                ],
                &kin,
            )
            .unwrap();
            let certificate = family.scaleless_scaling(&parameters).unwrap();
            assert_eq!(certificate.is_some(), scaleless);
            if let Some(weights) = certificate {
                assert_eq!(weights, vec![parse!("1/2"), Atom::one()]);
            }
        }
    }

    #[test]
    fn sector_support_excludes_numerators_and_finds_massless_subintegration() {
        let k = parse!("scale_sector::k");
        let l = parse!("scale_sector::l");
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone()])
            .unwrap();
        let kk = kin.scalar_product(&k, &k).unwrap();
        let family = IntegralFamily::new(
            vec![k, l.clone()],
            vec![],
            vec![
                kk.clone(),
                kin.scalar_product(&l, &l).unwrap() - parse!("m2"),
                kk - parse!("M2"),
            ],
            &kin,
        )
        .unwrap();
        assert!(
            family
                .scaleless_scaling(&[parse!("x"), parse!("y"), parse!("z")])
                .unwrap()
                .is_none()
        );
        for powers in [[1, 2, 0], [2, 1, -3]] {
            let sector = family.sector(&powers).unwrap();
            assert_eq!(sector.denominators(), &family.denominators()[..2]);
            assert_eq!(sector.loop_momenta(), family.loop_momenta());
            let weights = sector
                .scaleless_scaling(&[parse!("x"), parse!("y")])
                .unwrap()
                .unwrap();
            assert_eq!(weights, vec![Atom::one(), Atom::new()]);
        }
        assert!(family.sector(&[1]).is_err());
        let empty = family.sector(&[0, -1, 0]).unwrap();
        assert!(empty.denominators().is_empty());
        // A missing quadratic form remains an explicit unsupported case, not
        // a claim that this parametric certificate found no scalelessness.
        assert!(empty.scaleless_scaling(&[]).is_err());
    }

    #[test]
    fn parametric_mapping_beyond_fixed_external_momentum_shifts() {
        // FCLoopFindTopologyMappings: topos4 needs external momentum shifts.
        let [k, l, p, q] = ["ufmap::k", "ufmap::l", "ufmap::p", "ufmap::q"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone(), p.clone(), q.clone()])
            .unwrap()
            .with_mass_squared(&p, Atom::new())
            .unwrap()
            .with_mass_squared(&q, Atom::new())
            .unwrap()
            .with_scalar_product(&p, &q, parse!("s/2"))
            .unwrap();
        let square = |v: Atom| kin.scalar_product(&v, &v).unwrap();
        let mass = parse!("m2");
        let source = IntegralFamily::new(
            vec![k.clone(), l.clone()],
            vec![p.clone(), q.clone()],
            vec![
                square(&k + &p) - &mass,
                square(&k - &l),
                square(&l + &p) - &mass,
                square(&l - &q) - &mass,
                square(l.clone()),
            ],
            &kin,
        )
        .unwrap();
        let target = IntegralFamily::new(
            vec![k.clone(), l.clone()],
            vec![p.clone(), q.clone()],
            vec![
                square(&k - &l) - &mass,
                square(&k - &q),
                square(&l - &q) - &mass,
                square(&l + &p) - &mass,
                square(l.clone()),
            ],
            &kin,
        )
        .unwrap();
        assert!(source.find_mapping(&target, 10_000).unwrap().is_none());
        let parameters = ["x1", "x2", "x3", "x4", "x5"].map(|s| symbolica::symbol!(s).to_atom());
        let mapping = source
            .parametric_mapping(&target, &parameters)
            .unwrap()
            .unwrap();
        let mapped = mapping
            .denominator_map()
            .iter()
            .map(|i| parameters[*i].clone())
            .collect::<Vec<_>>();
        let (source_u, source_f) = source.symanzik(&mapped).unwrap();
        let (target_u, target_f) = target.symanzik(&parameters).unwrap();
        assert!((source_u - target_u).expand().is_zero());
        assert!((source_f - target_f).expand().is_zero());
    }

    #[test]
    fn parametric_permutations_preserve_coefficients_and_powers() {
        let k = parse!("ufperm::k");
        let l = parse!("ufperm::l");
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone()])
            .unwrap();
        let denominators = [&k, &l, &(&k - &l)]
            .into_iter()
            .zip([parse!("m1sq"), parse!("m2sq"), parse!("m3sq")])
            .map(|(p, m)| kin.scalar_product(p, p).unwrap() - m)
            .collect::<Vec<_>>();
        let family = IntegralFamily::new(
            vec![k.clone(), l.clone()],
            vec![],
            denominators.clone(),
            &kin,
        )
        .unwrap();
        let parameters = [parse!("x"), parse!("y"), parse!("z")];
        for permutation in (0..3).permutations(3) {
            let target = IntegralFamily::new(
                vec![k.clone(), l.clone()],
                vec![],
                permutation
                    .iter()
                    .map(|i| denominators[*i].clone())
                    .collect(),
                &kin,
            )
            .unwrap();
            let mapping = family
                .parametric_mapping(&target, &parameters)
                .unwrap()
                .unwrap();
            for (i, j) in mapping.denominator_map().iter().enumerate() {
                assert_eq!(permutation[*j], i);
            }
            let powers = [1, -2, 3];
            assert_eq!(
                mapping.map_powers(&powers).unwrap(),
                permutation.iter().map(|i| powers[*i]).collect::<Vec<_>>()
            );
            assert!(mapping.map_powers(&[1]).is_err());
        }
        let different = IntegralFamily::new(
            vec![k, l],
            vec![],
            denominators.iter().map(|d| d + 1).collect(),
            &kin,
        )
        .unwrap();
        assert!(
            family
                .parametric_mapping(&different, &parameters)
                .unwrap()
                .is_none()
        );
        assert!(
            family
                .parametric_mapping(&family, &[parse!("m1sq"), parse!("y"), parse!("z")])
                .is_err()
        );
    }

    #[test]
    fn vacuum_sunset_polynomials() {
        let k = parse!("ufv::k");
        let l = parse!("ufv::l");
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone()])
            .unwrap();
        let denominators = [&k, &l, &(&k - &l)]
            .into_iter()
            .zip([parse!("m1sq"), parse!("m2sq"), parse!("m3sq")])
            .map(|(p, mass)| kin.scalar_product(p, p).unwrap() - mass)
            .collect();
        let family = IntegralFamily::new(vec![k, l], vec![], denominators, &kin).unwrap();
        let (u, f) = family
            .symanzik(&[parse!("x"), parse!("y"), parse!("z")])
            .unwrap();
        assert_eq!(u, parse!("x*y+x*z+y*z"));
        assert!((f - &u * parse!("x*m1sq+y*m2sq+z*m3sq")).expand().is_zero());
    }

    #[test]
    fn massive_bubble_symanzik_and_shift_invariance() {
        let [k, p, x, y, s, m1, m2] = ["uf::k", "uf::p", "x", "y", "s", "m1sq", "m2sq"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, s.clone())
            .unwrap();
        let mut results = Vec::new();
        for shift in [Atom::new(), p.clone(), &p * 2] {
            let q = &k + shift;
            let family = IntegralFamily::new(
                vec![k.clone()],
                vec![p.clone()],
                vec![
                    kin.scalar_product(&q, &q).unwrap() - &m1,
                    kin.scalar_product(&(&q - &p), &(&q - &p)).unwrap() - &m2,
                ],
                &kin,
            )
            .unwrap();
            let (u, f) = family.symanzik(&[x.clone(), y.clone()]).unwrap();
            assert_eq!(u, &x + &y);
            assert!(
                (f.clone() - ((&x + &y) * (&x * &m1 + &y * &m2) - &s * &x * &y))
                    .expand()
                    .is_zero()
            );
            results.push(f);
            assert!(family.symanzik(&[x.clone(), x.clone()]).is_err());
            assert!(family.symanzik(&[m1.clone(), y.clone()]).is_err());
        }
        assert!(results.windows(2).all(|pair| pair[0] == pair[1]));
    }

    #[test]
    fn feyncalc_massless_two_loop_self_energy_polynomials() {
        // FCFeynmanPrepare manual: parameters follow its displayed propagator table.
        let [k, l, p] = ["uf2::k", "uf2::l", "uf2::p"].map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let denominators = [&k, &l, &(&k - &p), &(&l - &p), &(&k + &l - &p)]
            .into_iter()
            .map(|q| kin.scalar_product(q, q).unwrap())
            .collect();
        let family = IntegralFamily::new(vec![k, l], vec![p], denominators, &kin).unwrap();
        let parameters = ["x1", "x2", "x3", "x4", "x5"].map(|s| symbolica::symbol!(s).to_atom());
        let (u, f) = family.symanzik(&parameters).unwrap();
        assert_eq!(u, parse!("x1*x2+x2*x3+x2*x5+x1*x4+x3*x4+x1*x5+x3*x5+x4*x5"));
        assert!(
            (f - parse!(
                "-s*(x1*x2*x3+x1*x3*x4+x2*x3*x4+x1*x3*x5+x3*x4*x5+x1*x2*x4+x1*x2*x5+x2*x4*x5)"
            ))
            .expand()
            .is_zero()
        );
    }

    #[test]
    fn singular_polynomials_match_a_nonsingular_regulator_limit() {
        let [k, l, p] = ["uf_singular::k", "uf_singular::l", "uf_singular::p"]
            .map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new()
            .with_momenta([k.clone(), l.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let sum = &k + &l;
        let difference = &k - &l;
        let quadratic = kin.scalar_product(&sum, &sum).unwrap() - parse!("m2");
        let linear = kin.scalar_product(&difference, &p).unwrap() + parse!("delta");
        let family = IntegralFamily::new(
            vec![k.clone(), l.clone()],
            vec![p.clone()],
            vec![quadratic.clone(), linear.clone()],
            &kin,
        )
        .unwrap();
        let parameters = [parse!("x"), parse!("y")];
        let (u, f) = family.symanzik(&parameters).unwrap();
        assert!(u.is_zero());
        assert_eq!(f, parse!("s*x*y^2"));
        let eta = parse!("eta");
        let regular = IntegralFamily::new(
            vec![k, l.clone()],
            vec![p],
            vec![
                quadratic,
                linear + &eta * kin.scalar_product(&l, &l).unwrap(),
            ],
            &kin,
        )
        .unwrap();
        let (regular_u, regular_f) = regular.symanzik(&parameters).unwrap();
        assert!(!regular_u.is_zero());
        assert_eq!(regular_u.replace(eta.to_pattern()).with(Atom::new()), u);
        assert_eq!(regular_f.replace(eta.to_pattern()).with(Atom::new()), f);
        assert!(family.scaleless_scaling(&parameters).is_err());
        assert!(family.parametric_mapping(&family, &parameters).is_err());
    }

    #[test]
    fn absent_loop_directions_do_not_create_trivial_scaling_certificates() {
        let [k, l] = ["uf_absent::k", "uf_absent::l"].map(|s| symbolica::symbol!(s).to_atom());
        let kin = Kinematics::new();
        let family = IntegralFamily::new(
            vec![k.clone(), l],
            vec![],
            vec![kin.scalar_product(&k, &k).unwrap() - parse!("m2")],
            &kin,
        )
        .unwrap();
        let parameters = [parse!("x")];
        assert_eq!(
            family.symanzik(&parameters).unwrap(),
            (Atom::new(), Atom::new())
        );
        assert!(family.scaleless_scaling(&parameters).is_err());
        assert!(family.parametric_mapping(&family, &parameters).is_err());
    }

    #[test]
    fn eikonal_and_degenerate_quadratic_forms() {
        let k = parse!("ufe::k");
        let p = parse!("ufe::p");
        let kin = Kinematics::new()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let kk = kin.scalar_product(&k, &k).unwrap();
        let kp = kin.scalar_product(&k, &p).unwrap();
        let linear =
            IntegralFamily::new(vec![k.clone()], vec![p.clone()], vec![kp.clone()], &kin).unwrap();
        assert_eq!(
            linear.symanzik(&[parse!("x")]).unwrap(),
            (Atom::new(), parse!("s*x^2/4")),
        );
        assert!(linear.scaleless_scaling(&[parse!("x")]).is_err());
        assert!(linear.parametric_mapping(&linear, &[parse!("x")]).is_err());
        let mixed =
            IntegralFamily::new(vec![k], vec![p], vec![kk, kp + parse!("delta")], &kin).unwrap();
        let (u, f) = mixed.symanzik(&[parse!("x"), parse!("y")]).unwrap();
        assert_eq!(u, parse!("x"));
        assert_eq!(f, parse!("s*y^2/4-x*y*delta"));
    }

    #[test]
    fn transverse_direction_certifies_mixed_eikonal_sectors() {
        let [k, q, p] = [
            parse!("transverse::k"),
            parse!("transverse::q"),
            parse!("transverse::p"),
        ];
        let kin = Kinematics::new()
            .with_momenta([k.clone(), q.clone(), p.clone()])
            .unwrap()
            .with_mass_squared(&p, parse!("s"))
            .unwrap();
        let denominators = vec![
            kin.scalar_product(&(&k + &q), &(&k + &q)).unwrap() - parse!("m2"),
            kin.scalar_product(&(&k - &q), &p).unwrap() + parse!("delta"),
        ];
        let family = IntegralFamily::new(
            vec![k.clone(), q.clone()],
            vec![p.clone()],
            denominators.clone(),
            &kin,
        )
        .unwrap();
        assert_eq!(
            family.scaleless_transverse_direction().unwrap(),
            Some(vec![Atom::one(), Atom::num(-1)])
        );
        let mut constrained = denominators.clone();
        constrained.push(kin.scalar_product(&(&k - &q), &(&k - &q)).unwrap());
        assert!(
            IntegralFamily::new(
                vec![k.clone(), q.clone()],
                vec![p.clone()],
                constrained,
                &kin
            )
            .unwrap()
            .scaleless_transverse_direction()
            .unwrap()
            .is_none()
        );
        let null_kin = kin.clone().with_mass_squared(&p, Atom::Zero).unwrap();
        assert!(
            IntegralFamily::new(
                vec![k.clone(), q.clone()],
                vec![p.clone()],
                denominators,
                &null_kin
            )
            .unwrap()
            .scaleless_transverse_direction()
            .unwrap()
            .is_none()
        );
        let one_dimensional = Kinematics::in_dimension(&Atom::one())
            .unwrap()
            .with_mass_squared(&p, Atom::one())
            .unwrap();
        let linear = one_dimensional.scalar_product(&k, &p).unwrap();
        assert!(
            IntegralFamily::new(vec![k.clone()], vec![p], vec![linear], &one_dimensional)
                .unwrap()
                .scaleless_transverse_direction()
                .unwrap()
                .is_none()
        );
        let complex = &k + parse!("𝑖") * &q;
        let family = IntegralFamily::new(
            vec![k, q],
            vec![],
            vec![kin.scalar_product(&complex, &complex).unwrap()],
            &kin,
        )
        .unwrap();
        assert!(family.scaleless_transverse_direction().unwrap().is_none());
    }

    #[test]
    fn singular_matrix_pencil_without_common_kernel_is_not_a_certificate() {
        let [k, q, r] = [
            parse!("pencil::k"),
            parse!("pencil::q"),
            parse!("pencil::r"),
        ];
        let kin = Kinematics::new();
        // det([[x,y,z],[y,0,0],[z,0,0]]) = 0, but its null direction
        // depends on the parameters. No loop direction leaves each form fixed.
        let denominators = [&k, &q, &r]
            .map(|p| kin.scalar_product(&k, p).unwrap())
            .to_vec();
        let family = IntegralFamily::new(vec![k, q, r], vec![], denominators, &kin).unwrap();
        let (u, _) = family
            .symanzik(&[parse!("x"), parse!("y"), parse!("z")])
            .unwrap();
        assert!(u.is_zero());
        assert!(family.scaleless_transverse_direction().unwrap().is_none());
    }
}
