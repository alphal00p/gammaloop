//! Exact sufficient certificates for soft/Fermi transversality.
//!
//! On a frozen soft span, the remaining scalar products are constrained
//! linearly by the squared Fermi radii. A nonzero reduced Gram determinant
//! certifies transversality, or absence of real configurations if it is
//! negative. An inconsistent norm equation also certifies absence. Unresolved
//! scalar products are deliberately not used as a positivity certificate.

use bincode_trait_derive::{Decode, Encode};
use symbolica::{domains::atom::AtomField, tensors::matrix::Matrix};

use crate::GammaLoopContext;

use super::*;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct FermiSurfaceDomain {
    pub active_edges: Vec<EdgeIndex>,
    pub soft_bases: Vec<Vec<EdgeIndex>>,
    pub certificate: Atom,
}

impl FermiSurfaceDomain {
    /// A certificate parameter, resolved from the actual runtime mass and
    /// chemical potential before exact rational evaluation.
    pub(crate) fn radius_squared(edge: EdgeIndex) -> Atom {
        function!(symbol!("gammalooprs::fermi_radius_squared"), edge.0)
    }
}

impl FermiSurfaceLocalizer<'_> {
    pub(super) fn soft_domain(
        &self,
        support: &[FermiDistribution],
        basis: &LoopMomentumBasis,
        free: &[LoopIndex],
        soft_bases: Vec<Vec<EdgeIndex>>,
    ) -> Result<Option<FermiSurfaceDomain>> {
        if soft_bases.iter().all(Vec::is_empty) {
            return Ok(None);
        }
        let mut domain = FermiSurfaceDomain {
            active_edges: support.iter().map(|factor| factor.edge).collect(),
            soft_bases,
            certificate: Atom::Zero,
        };
        let field = AtomField {
            statistical_zero_test: false,
            cancel_check_on_division: true,
            custom_normalization: None,
        };
        let routing = support
            .iter()
            .map(|factor| {
                free.iter()
                    .map(|index| {
                        Self::sign_atom(basis.edge_signatures[factor.edge].internal[*index])
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let routing_gram = routing
            .iter()
            .map(|a| {
                routing
                    .iter()
                    .map(|b| Atom::add_many(a.iter().zip(b).map(|(x, y)| x * y)))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let routing_det = Matrix::from_nested_vec(routing_gram.clone(), field.clone())
            .map_err(|error| eyre!(error))?
            .det()
            .map_err(|error| eyre!("{error}"))?;
        if !routing_det.is_zero() {
            // Independent momentum rows imply independent radial normals at
            // every nonzero active radius, including affine external shifts.
            domain.certificate = Atom::one();
            return Ok(Some(domain));
        }
        if !basis.ext_edges.is_empty() {
            // A rank-deficient affine system needs external scalar products;
            // this certificate currently proves homogeneous vacuum geometry.
            return Ok(Some(domain));
        }

        let dot_head = symbol!("gammalooprs::fermi_certificate_scalar_product");
        let scalar_products = (0..free.len())
            .flat_map(|a| (a..free.len()).map(move |b| (a, b)))
            .map(|(a, b)| ((a, b), function!(dot_head, a, b)))
            .collect::<Vec<_>>();
        let unknowns = scalar_products.len();
        let equations = routing
            .iter()
            .zip(support)
            .map(|(row, factor)| {
                scalar_products
                    .iter()
                    .map(|((a, b), _)| Atom::num(if a == b { 1 } else { 2 }) * &row[*a] * &row[*b])
                    .chain([FermiSurfaceDomain::radius_squared(factor.edge)])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut equations =
            Matrix::from_nested_vec(equations, field.clone()).map_err(|error| eyre!(error))?;
        equations.row_reduce(unknowns as u32);
        let mut substitutions = BTreeMap::new();
        let mut consistency = Atom::Zero;
        for row in equations.row_iter() {
            let Some(pivot) = row[..unknowns].iter().position(|entry| !entry.is_zero()) else {
                // A nonzero right-hand side in 0 = rhs rules out this soft
                // stratum. This includes an active momentum forced to zero.
                consistency += row[unknowns].pow(2);
                continue;
            };
            let solution = &row[unknowns]
                - Atom::add_many(
                    scalar_products
                        .iter()
                        .enumerate()
                        .filter(|(index, _)| *index != pivot)
                        .map(|(index, (_, variable))| &row[index] * variable),
                );
            substitutions.insert(scalar_products[pivot].1.clone(), solution);
        }

        let gram = routing
            .iter()
            .enumerate()
            .map(|(a, qa)| {
                routing
                    .iter()
                    .enumerate()
                    .map(|(b, qb)| {
                        &routing_gram[a][b]
                            * Atom::add_many(scalar_products.iter().map(|((i, j), variable)| {
                                let coefficient = if i == j {
                                    &qa[*i] * &qb[*j]
                                } else {
                                    &qa[*i] * &qb[*j] + &qa[*j] * &qb[*i]
                                };
                                coefficient * variable
                            }))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let determinant = Matrix::from_nested_vec(gram, field)
            .map_err(|error| eyre!(error))?
            .det()
            .map_err(|error| eyre!("{error}"))?
            .replace_map(|part, _, out| {
                if let Some(value) = substitutions.get(&part.to_owned()) {
                    **out = value.clone();
                }
            })
            .together()
            .expand();
        if !determinant.contains_symbol(dot_head) {
            consistency += determinant.pow(2);
        }
        domain.certificate = consistency.together().expand();
        Ok(Some(domain))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, utils::load_generic_model};

    #[test]
    fn fermi_surface_soft_domain_certifies_triangle_and_absent_strata() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        let graph = Graph::from_string(
            r#"digraph triangle_soft_domain {
                node [num=1]; edge [num=1];
                A -> B [id=0, particle="d", mass=0];
                A -> B [id=1, particle="u", mass=0];
                A -> B [id=2, particle="s", mass=0];
                B -> A [id=3, particle="g", mass=0];
            }"#,
            &model,
        )?
        .remove(0);
        let support = (0..3)
            .map(|edge| FermiDistribution {
                edge: EdgeIndex(edge),
                orientation: Atom::one(),
                order: 2,
            })
            .collect::<Vec<_>>();
        let soft = EdgeIndex(3);
        let basis = graph.lmb_with_loop_edges(&[soft][..])?;
        let free = basis
            .loop_edges
            .iter_enumerated()
            .filter_map(|(index, edge)| (*edge != soft).then_some(index))
            .collect::<Vec<_>>();
        let localizer = FermiSurfaceLocalizer::new(&graph);
        let domain = localizer
            .soft_domain(&support, &basis, &free, vec![vec![soft]])?
            .unwrap();
        for (radii_squared, certified) in [([1, 1, 1], true), ([1, 1, 4], false), ([1, 1, 9], true)]
        {
            let value = support.iter().zip(radii_squared).fold(
                domain.certificate.clone(),
                |atom, (factor, radius)| {
                    atom.replace(FermiSurfaceDomain::radius_squared(factor.edge))
                        .with(radius)
                },
            );
            assert_eq!(
                !value.is_zero(),
                certified,
                "radii²={radii_squared:?}: {value}"
            );
            if radii_squared == [1, 1, 1] {
                assert_eq!(value, Atom::num((9, 4)));
            }
        }
        let frozen = vec![EdgeIndex(1), EdgeIndex(2), EdgeIndex(3)];
        let frozen_basis = graph.lmb_with_loop_edges(frozen.as_slice())?;
        assert_eq!(
            localizer
                .soft_domain(&support[..1], &frozen_basis, &[], vec![frozen])?
                .unwrap()
                .certificate,
            FermiSurfaceDomain::radius_squared(support[0].edge).pow(2),
            "forcing an active momentum to zero is excluded by its nonzero radius",
        );
        Ok(())
    }

    #[test]
    fn fermi_surface_soft_domain_distinguishes_coincident_and_disjoint_radii() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        let graph = Graph::from_string(
            r#"digraph coincident_soft_domain {
                node [num=1]; edge [num=1];
                A -> B [id=0, particle="d", mass=0];
                A -> B [id=1, particle="u", mass=0];
                B -> A [id=2, particle="g", mass=0];
            }"#,
            &model,
        )?
        .remove(0);
        let support = (0..2)
            .map(|edge| FermiDistribution {
                edge: EdgeIndex(edge),
                orientation: Atom::one(),
                order: 2,
            })
            .collect::<Vec<_>>();
        let soft = EdgeIndex(2);
        let basis = graph.lmb_with_loop_edges(&[soft][..])?;
        let free = basis
            .loop_edges
            .iter_enumerated()
            .filter_map(|(index, edge)| (*edge != soft).then_some(index))
            .collect::<Vec<_>>();
        let domain = FermiSurfaceLocalizer::new(&graph)
            .soft_domain(&support, &basis, &free, vec![vec![soft]])?
            .unwrap();
        for (second_radius_squared, expected) in [(1, 0), (4, 9)] {
            let value = domain
                .certificate
                .replace(FermiSurfaceDomain::radius_squared(EdgeIndex(0)))
                .with(1)
                .replace(FermiSurfaceDomain::radius_squared(EdgeIndex(1)))
                .with(second_radius_squared);
            assert_eq!(value, Atom::num(expected));
        }
        Ok(())
    }
}
