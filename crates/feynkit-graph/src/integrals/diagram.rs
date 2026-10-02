use std::collections::BTreeMap;

use feynkit_kinematics::Kinematics;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView};

use super::{IntegralFamily, IntegralFamilyError};
use crate::{DiagramError, FeynmanDiagram, symbols};

impl IntegralFamily {
    /// Complete the diagram's internal propagators with auxiliary scalar products.
    ///
    /// Propagators retain their edge order and stored momentum routing, as in
    /// [`FeynmanDiagram::integral_family`]. Preferred independent dot products
    /// are tried in order; redundant entries are skipped and automatic scalar
    /// products fill any remaining directions. Pass an empty slice to choose
    /// the completion automatically. Dependent propagators must be extracted
    /// and partial-fractioned separately before completing their families.
    pub fn from_diagram(
        diagram: &FeynmanDiagram,
        kinematics: &Kinematics,
        independent_dot_products: &[Atom],
    ) -> Result<Self, DiagramError> {
        diagram.integral_family(kinematics, independent_dot_products)
    }
}

impl FeynmanDiagram {
    /// Complete the graph propagators with preferred or automatic scalar products.
    ///
    /// Physical propagators retain their edge order. Candidates are tried in
    /// order, skipping redundant entries; automatic products fill any remaining
    /// directions. Pass an empty slice for automatic completion. Dependent
    /// propagators must be extracted with [`Self::propagator_family`] and
    /// partial-fractioned before completing the resulting families.
    pub fn integral_family(
        &self,
        kinematics: &Kinematics,
        independent_dot_products: &[Atom],
    ) -> Result<IntegralFamily, DiagramError> {
        let family = self.propagator_family(kinematics)?;
        Ok(family.complete(independent_dot_products)?)
    }

    /// Build an integral family from the diagram's internal quadratic propagators.
    ///
    /// Denominators follow ascending internal edge IDs, retaining repeated
    /// propagators and bridges. The stored loop routing supplies momentum names;
    /// dependent external coordinates are eliminated by that same routing.
    /// External carriers and dummy edges are excluded by the shared denominator
    /// builder. Symbolic masses follow its conventions, including UFO ZERO.
    /// Widths, prescriptions and custom UFO denominator formulas are not inferred.
    /// Tree diagrams have no loop-integral family and return an error.
    /// Use [`IntegralFamily::from_diagram`] to append auxiliary scalar products
    /// and obtain a complete independent denominator basis.
    pub fn propagator_family(
        &self,
        kinematics: &Kinematics,
    ) -> Result<IntegralFamily, DiagramError> {
        let basis = self.loop_momentum_basis();
        if basis.loop_edges.is_empty() {
            return Err(IntegralFamilyError::NoLoops.into());
        }
        let loops = (0..basis.loop_edges.len())
            .map(|i| symbols::loop_momentum().call(i))
            .collect::<Vec<_>>();
        let externals = basis
            .external_edges
            .iter()
            .enumerate()
            .filter(|(_, edge)| !basis.dependent_externals.contains(edge))
            .map(|(i, _)| symbols::external_momentum().call(i))
            .collect::<Vec<_>>();
        let kin = kinematics
            .clone()
            .with_momenta(loops.iter().chain(&externals).cloned())
            .map_err(IntegralFamilyError::from)?;
        let annotated = self.denominator_of_in_dimension(
            &self.internal_subgraph(),
            &BTreeMap::new(),
            kin.dimension(),
        )?;
        let factors = match annotated.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        let rep = Minkowski {}.new_rep(kin.dimension());
        let mut denominators = BTreeMap::new();
        for factor in factors {
            let AtomView::Fun(propagator) = factor else {
                return Err(DiagramError::Invariant {
                    operation: "extracting an integral family",
                    message: "expected a tagged propagator from the shared denominator builder"
                        .into(),
                });
            };
            let arguments = propagator.iter().map(|a| a.to_owned()).collect::<Vec<_>>();
            if propagator.get_symbol() != symbols::denominator() || arguments.len() != 4 {
                return Err(DiagramError::Invariant {
                    operation: "extracting an integral family",
                    message: "expected a four-argument propagator annotation".into(),
                });
            }
            let edge = usize::try_from(arguments[0].as_view()).map_err(|error| {
                DiagramError::Invariant {
                    operation: "extracting an integral family",
                    message: format!("invalid propagator edge ID: {error}"),
                }
            })?;
            let momentum = basis.route_expression(&arguments[1]);
            let square = kin
                .scalar_product(&momentum, &momentum)
                .map_err(IntegralFamilyError::from)?;
            // Replace only the quadratic scalar product in the shared formula,
            // so propagator construction and model mass conventions stay owned
            // by denominator_of_in_dimension rather than a second implementation.
            let original_square = rep.inner_product(&arguments[1], &arguments[1]);
            let denominator = arguments[3].replace(original_square).with(square);
            denominators.insert(edge, denominator);
        }
        Ok(IntegralFamily::new(
            loops,
            externals,
            denominators.into_values().collect(),
            &kin,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_model::Model;

    #[test]
    fn dependent_graph_propagators_remain_available_for_partial_fractioning() {
        let diagram = FeynmanDiagram::from_dot(
            Model::phi4(),
            r#"digraph tadpole_insertion {
                ext [style=invis];
                ext -> a [particle="phi"];
                a -> ext [particle="phi"];
                a -> b [particle="phi"];
                a -> b [particle="phi"];
                b -> b [particle="phi"];
            }"#,
        )
        .unwrap();
        let kin = Kinematics::new();
        let raw = diagram.propagator_family(&kin).unwrap();
        assert_eq!(raw.denominators().len(), 3);
        assert!(!raw.is_independent());
        assert!(matches!(
            diagram.integral_family(&kin, &[]),
            Err(DiagramError::IntegralFamily(
                IntegralFamilyError::Dependent { .. }
            ))
        ));
    }

    #[test]
    fn diagram_family_completes_sunrise_with_preferred_products() {
        let diagram = FeynmanDiagram::from_dot(
            Model::phi4(),
            r#"digraph sunrise {
                ext [style=invis];
                ext -> a [particle="phi"];
                a -> b [particle="phi", lmb_id=0];
                a -> b [particle="phi", lmb_id=1];
                a -> b [particle="phi"];
                b -> ext [particle="phi"];
            }"#,
        )
        .unwrap();
        let kin = Kinematics::new();
        let raw = diagram.propagator_family(&kin).unwrap();
        assert_eq!(raw.denominators().len(), 3);
        assert_eq!(raw.scalar_products().len(), 5);
        assert!(!raw.is_complete());
        let automatic = diagram.integral_family(&kin, &[]).unwrap();
        assert!(automatic.is_complete() && automatic.is_independent());
        assert_eq!(&automatic.denominators()[..3], raw.denominators());

        let products = raw
            .loop_momenta()
            .iter()
            .map(|k| {
                raw.kinematics()
                    .scalar_product(k, &raw.external_momenta()[0])
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let preferred = diagram.integral_family(&kin, &products).unwrap();
        let constructed = IntegralFamily::from_diagram(&diagram, &kin, &products).unwrap();
        assert_eq!(preferred.denominators(), constructed.denominators());
        assert_eq!(&preferred.denominators()[..3], raw.denominators());
        assert_eq!(&preferred.denominators()[3..], products);
        assert!(preferred.is_complete() && preferred.is_independent());
        let partial = IntegralFamily::from_diagram(
            &diagram,
            &kin,
            &[raw.denominators()[0].clone(), products[1].clone()],
        )
        .unwrap();
        assert_eq!(partial.denominators()[3], products[1]);
        assert!(partial.is_complete() && partial.is_independent());

        // Validate even candidates after the basis is already complete.
        let mut invalid = products;
        invalid.push(invalid[0].pow(2));
        assert!(IntegralFamily::from_diagram(&diagram, &kin, &invalid).is_err());
    }
}
