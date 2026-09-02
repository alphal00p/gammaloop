use std::ops::Neg;

use color_eyre::Result;
use eyre::eyre;
use symbolica::atom::Atom;

use crate::{
    cff::{expression::OrientationID, surface::LinearEnergyExpr},
    graph::Graph,
    uv::{
        Integrands,
        approx::{OrientationProjection, local_3d::OrientationIntegrands},
    },
};

/// The energy substitution belonging to one complete generalized residue key.
///
/// Stored production branches recover their map from `selector_host`. Reduced
/// sources retain their explicit map because it is part of their residue key,
/// not temporary projection metadata. Every factor in a branch is mapped with
/// this one authority; the Taylor operator only changes its rescaling series.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum DirectEnergyMap {
    Production,
    Source(Vec<LinearEnergyExpr>),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectResidueKey {
    pub(crate) selector_host: OrientationID,
    pub(crate) energy_map: DirectEnergyMap,
}

impl DirectResidueKey {
    pub(crate) fn production(selector_host: OrientationID) -> Self {
        Self {
            selector_host,
            energy_map: DirectEnergyMap::Production,
        }
    }

    pub(crate) fn source(
        selector_host: OrientationID,
        edge_energy_map: Vec<LinearEnergyExpr>,
    ) -> Self {
        Self {
            selector_host,
            energy_map: DirectEnergyMap::Source(edge_energy_map),
        }
    }

    pub(crate) fn map_numerator(
        &self,
        orientation: OrientationProjection<'_>,
        graph: &Graph,
        numerator: &Atom,
    ) -> Result<Atom> {
        let source_map = match &self.energy_map {
            DirectEnergyMap::Production => None,
            DirectEnergyMap::Source(map) => Some(map.as_slice()),
        };
        orientation.map_numerator(graph, self.selector_host, source_map, numerator)
    }

    pub(crate) fn source_edge_energy_map(&self) -> Option<&[LinearEnergyExpr]> {
        match &self.energy_map {
            DirectEnergyMap::Production => None,
            DirectEnergyMap::Source(map) => Some(map),
        }
    }

    fn selector(&self, materialize_key_selector: bool) -> Atom {
        if materialize_key_selector {
            self.selector_host.atom()
        } else {
            Atom::one()
        }
    }
}

/// Factorized direct-3D bodies grouped by their complete generalized residue
/// key. Distinct source maps hosted by the same production selector remain
/// distinct branches, while exact duplicate keys are coalesced additively.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct DirectResidueBranches(Vec<(DirectResidueKey, Integrands)>);

impl DirectResidueBranches {
    #[cfg(test)]
    pub(crate) fn production(selector_host: OrientationID, integrands: Integrands) -> Result<Self> {
        Self::from_keyed([(DirectResidueKey::production(selector_host), integrands)])
    }

    pub(crate) fn from_transient(source: &OrientationIntegrands) -> Result<Self> {
        Self::from_keyed(source.iter_orientations().map(
            |(selector_host, source_map, integrands)| {
                let key = source_map.map_or_else(
                    || DirectResidueKey::production(selector_host),
                    |map| DirectResidueKey::source(selector_host, map.to_vec()),
                );
                (key, integrands.clone())
            },
        ))
    }

    pub(super) fn from_keyed(
        keyed: impl IntoIterator<Item = (DirectResidueKey, Integrands)>,
    ) -> Result<Self> {
        let mut branches: Vec<(DirectResidueKey, Integrands)> = Vec::new();
        for (key, integrands) in keyed {
            if let Some((_, existing)) = branches
                .iter_mut()
                .find(|(existing_key, _)| *existing_key == key)
            {
                *existing = existing.clone().zip_add(integrands)?;
            } else {
                branches.push((key, integrands));
            }
        }
        let fallback_zero = branches
            .iter()
            .find(|(_, integrands)| integrands.iter().all(|(_, atom)| atom.is_zero()))
            .cloned();
        branches.retain(|(_, integrands)| !integrands.iter().all(|(_, atom)| atom.is_zero()));
        if branches.is_empty() {
            branches.extend(fallback_zero);
        }
        if branches.is_empty() {
            return Err(eyre!("direct local-3D residue branches cannot be empty"));
        }
        Ok(Self(branches))
    }

    pub(crate) fn iter_keys(&self) -> impl Iterator<Item = (&DirectResidueKey, &Integrands)> {
        self.0.iter().map(|(key, integrands)| (key, integrands))
    }

    /// Build the multiplicative identity on this branch family's exact cut
    /// support. A selected raised residue need not contain the uncut root key.
    pub(crate) fn identity_integrands(&self) -> Integrands {
        self.0
            .first()
            .expect("direct residue branches are never empty")
            .1
            .iter()
            .map(|(index, _)| (*index, Atom::one()))
            .collect()
    }

    #[cfg(test)]
    pub(crate) fn factorized_sum(&self) -> Atom {
        self.0
            .iter()
            .flat_map(|(_, integrands)| integrands.iter().map(|(_, atom)| atom))
            .fold(Atom::Zero, |sum, atom| sum + atom)
    }

    #[cfg(test)]
    pub(crate) fn factorized_capacity_envelope(&self) -> Atom {
        self.0
            .iter()
            .flat_map(|(_, integrands)| integrands.iter().map(|(_, atom)| atom))
            .filter(|atom| !atom.is_zero())
            .enumerate()
            .fold(Atom::Zero, |sum, (branch, atom)| {
                let tag = Atom::var(symbolica::symbol!(format!(
                    "__gammaloop_direct_3d_capacity_branch_{branch}"
                )));
                sum + tag * atom
            })
    }

    pub(crate) fn zip_add(&self, other: &Self) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| (key.clone(), integrands.clone()))
                .chain(
                    other
                        .iter_keys()
                        .map(|(key, integrands)| (key.clone(), integrands.clone())),
                ),
        )
    }

    pub(crate) fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| Ok((key.clone(), integrands.zip_mul(other)?)))
                .collect::<Result<Vec<_>>>()?,
        )
    }

    pub(crate) fn map(&self, mut map: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|(key, integrands)| (key.clone(), integrands.map(&mut map)))
                .collect(),
        )
    }

    pub(crate) fn fallible_map(
        &self,
        mut map: impl FnMut(&DirectResidueKey, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        Self::from_keyed(
            self.iter_keys()
                .map(|(key, integrands)| {
                    Ok((key.clone(), integrands.fallible_map(|atom| map(key, atom))?))
                })
                .collect::<Result<Vec<_>>>()?,
        )
    }

    pub(crate) fn multiply_key_mapped(
        &self,
        orientation: OrientationProjection<'_>,
        graph: &Graph,
        factor: &Atom,
    ) -> Result<Self> {
        self.fallible_map(|key, atom| Ok(atom * key.map_numerator(orientation, graph, factor)?))
    }

    /// Materialize residue selectors only at the evaluator boundary.
    pub(crate) fn materialize(&self, materialize_key_selector: bool) -> Result<Integrands> {
        let mut branches = self.iter_keys();
        let (first_key, first) = branches
            .next()
            .ok_or_else(|| eyre!("direct local-3D residue branches cannot be empty"))?;
        let first_selector = first_key.selector(materialize_key_selector);
        let mut result = first.map(|atom| atom * &first_selector);
        for (key, integrands) in branches {
            let selector = key.selector(materialize_key_selector);
            result = result.zip_add(integrands.map(|atom| atom * &selector))?;
        }
        Ok(result)
    }
}

impl Neg for DirectResidueBranches {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|(key, integrands)| (key, -integrands))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
    use symbolica::{atom::AtomView, symbol};
    use typed_index_collections::TiVec;

    use crate::{
        cff::{
            CutCFFIndex,
            expression::{OrientationData, OrientationExpression},
        },
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
        utils::GS,
    };

    use super::*;

    #[test]
    fn one_source_residue_map_is_a_homomorphism_for_all_factorized_factors() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(
            digraph G {
                edge [particle="scalar_1"];
                node [num=1];
                a -> b [id=0];
                a -> b [id=1];
            },
            "scalars"
        )?;
        let production = vec![OrientationExpression {
            data: OrientationData::new(EdgeVec::from_iter([
                Orientation::Default,
                Orientation::Reversed,
            ])),
            loop_energy_map: Vec::new(),
            edge_energy_map: vec![
                LinearEnergyExpr::ose(EdgeIndex(0), 1),
                LinearEnergyExpr::ose(EdgeIndex(1), -1),
            ],
            variants: Vec::new(),
        }]
        .into_iter()
        .collect::<TiVec<OrientationID, _>>();
        let options = graph.denominator_only_cff_3d_expression_options();
        let pattern = OrientationPattern::default();
        let orientation = OrientationProjection::exact(&production, &options, &pattern, false);

        let a = Atom::var(symbol!("direct_3d_test::A"));
        let b = Atom::var(symbol!("direct_3d_test::B"));
        let host = OrientationID(0);
        let index = CutCFFIndex::new_all_none();
        let key = DirectResidueKey::source(
            host,
            vec![LinearEnergyExpr::uniform_scale(2), LinearEnergyExpr::zero()],
        );
        let q0 = GS.emr_mom(EdgeIndex(0), GS.cind(0));

        // Map both factors independently through the same complete residue
        // key. Keeping them factorized must not change that common Q0 -> 2M
        // substitution authority.
        let first = key.map_numerator(orientation, &graph, &(q0.clone() + &a))?;
        let branches = DirectResidueBranches::from_keyed([(
            key,
            [(index, first.clone())].into_iter().collect(),
        )])?;
        let completed = branches.fallible_map(|key, body| {
            let later = key.map_numerator(orientation, &graph, &(q0.clone() + &b))?;
            Ok(body * later)
        })?;
        let sampled_energy = Atom::num(2) * Atom::var(GS.numerator_sampling_scale);
        let later = &sampled_energy + &b;
        assert_eq!(first, &sampled_energy + &a);
        let expected = first * later;

        let explicit = completed.materialize(false)?;
        assert_eq!(explicit.iter().next().unwrap().1, &expected);
        let AtomView::Mul(factors) = explicit.iter().next().unwrap().1.as_view() else {
            panic!("the common source map must retain the factorized product")
        };
        assert_eq!(
            factors
                .iter()
                .filter(|factor| matches!(factor, AtomView::Add(_)))
                .count(),
            2
        );

        let localized = completed.materialize(true)?;
        let localized_body = localized.iter().next().unwrap().1;
        assert_eq!(host.select(localized_body.as_view()), expected);
        assert_eq!(
            OrientationID(1).select(localized_body.as_view()),
            Atom::Zero
        );
        Ok(())
    }

    #[test]
    fn selected_cut_support_survives_root_to_first_sector_identity() -> Result<()> {
        let selected = CutCFFIndex {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: Some(1),
        };
        let body = Atom::var(symbol!("direct_3d_test::selected_cut_body"));
        let root = DirectResidueBranches::production(
            OrientationID(0),
            [(selected, body.clone())].into_iter().collect(),
        )?;

        let identity = root.identity_integrands();
        assert_eq!(
            identity.iter().collect::<Vec<_>>(),
            vec![(&selected, &Atom::one())]
        );
        let first_sector = root.zip_mul_unmapped(&identity)?;
        assert_eq!(first_sector.factorized_sum(), body);
        Ok(())
    }
}
