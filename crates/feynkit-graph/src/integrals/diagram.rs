use std::collections::BTreeMap;

use feynkit_kinematics::Kinematics;
use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};

use super::{IntegralFamily, IntegralFamilyError};
use crate::{DiagramError, FeynmanDiagram, LoopMomentumBasis, symbols};

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
    /// Edges incident to external vertices are excluded. Symbolic masses use
    /// the model's UFO parameter convention, including literal zero for ZERO.
    /// Widths, prescriptions and custom UFO denominator formulas are not inferred.
    /// Tree diagrams have no loop-integral family and return an error.
    /// Use [`IntegralFamily::from_diagram`] to append auxiliary scalar products
    /// and obtain a complete independent denominator basis.
    pub fn propagator_family(
        &self,
        kinematics: &Kinematics,
    ) -> Result<IntegralFamily, DiagramError> {
        self.validate()?;
        let basis = self.loop_momentum_basis();
        if basis.loop_edges.is_empty() {
            return Err(IntegralFamilyError::NoLoops.into());
        }
        let loops = (0..basis.loop_edges.len())
            .map(|i| symbols::loop_momentum().call(i))
            .collect::<Vec<_>>();
        let all_externals = (0..basis.external_edges.len())
            .map(|i| symbols::external_momentum().call(i))
            .collect::<Vec<_>>();
        let externals = basis
            .external_edges
            .iter()
            .enumerate()
            .filter(|(_, edge)| !basis.dependent_externals.contains(edge))
            .map(|(i, _)| all_externals[i].clone())
            .collect::<Vec<_>>();
        let kin = kinematics
            .clone()
            .with_momenta(loops.iter().chain(&externals).cloned())
            .map_err(IntegralFamilyError::from)?;
        let edges = self
            .edges()
            .map(|(id, _, edge)| (id, edge))
            .collect::<BTreeMap<_, _>>();
        let mut internal = self.topology().internal_edges;
        internal.sort();
        let mut denominators = Vec::with_capacity(internal.len());
        for edge_id in internal {
            // Signatures already eliminate dependent external momenta. Apply
            // the full coordinate arrays so their original indices never shift.
            let signature = basis.edge_signatures.get(&edge_id).ok_or_else(|| {
                DiagramError::InvalidLoopMomentumBasis(format!(
                    "missing signature for edge {}",
                    edge_id.0
                ))
            })?;
            let momentum = signature
                .apply(&loops, &all_externals)
                .map_err(|error| DiagramError::InvalidLoopMomentumBasis(error.to_string()))?
                .unwrap_or(Atom::Zero);
            let square = kin
                .scalar_product(&momentum, &momentum)
                .map_err(IntegralFamilyError::from)?;
            let edge = edges[&edge_id];
            let mass = self
                .model()
                .particle_by_id(edge.particle)?
                .symbolic_mass_parameter(self.model())?;
            denominators.push(square - mass.pow(2));
        }
        Ok(IntegralFamily::new(loops, externals, denominators, &kin)?)
    }
}

impl LoopMomentumBasis {
    /// Substitute stored edge signatures, retaining every trailing tensor argument.
    pub fn momentum_replacements(&self) -> Vec<symbolica::id::Replacement> {
        let arguments = Atom::var(symbol!("feynkit_graph::routing_arguments___"));
        let component = |signature: &feynkit_kinematics::Signature,
                         head: symbolica::atom::Symbol| {
            signature
                .integer_coefficients()
                .into_iter()
                .enumerate()
                .filter(|(_, coefficient)| *coefficient != 0)
                .map(|(index, coefficient)| {
                    Atom::num(coefficient as i64) * head.call((index, &arguments))
                })
                .fold(Atom::Zero, |sum, term| sum + term)
        };
        self.edge_signatures
            .iter()
            .map(|(edge, signature)| {
                let routed = component(&signature.loops, symbols::loop_momentum())
                    + component(&signature.external, symbols::external_momentum());
                symbolica::id::Replacement::new(
                    symbols::momentum().call((edge.0, &arguments)).to_pattern(),
                    routed.to_pattern(),
                )
            })
            .collect()
    }

    /// Route edge momenta without distributing products of numerator factors.
    pub fn route_expression(&self, expression: &Atom) -> Atom {
        expression
            .replace_multiple(self.momentum_replacements())
            .expand_num()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DiagramEdge, DiagramVertex, EdgeId, ExternalState};

    fn diagram(insertion: bool, reverse_tree_edge: bool) -> FeynmanDiagram {
        let model = crate::tests::scalar_model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "family-routing");
        let incoming =
            builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let internal = |name: &str| DiagramVertex {
            name: name.into(),
            interaction: None,
            external: None,
            numerator: Atom::one(),
        };
        let left = builder.add_vertex(internal("left"));
        let right = builder.add_vertex(internal("right"));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("out", 1, ExternalState::Outgoing));
        let edge = || DiagramEdge::new(particle, false);
        builder.add_edge(incoming, left, edge()).unwrap();
        builder.add_edge(left, right, edge()).unwrap();
        builder.add_edge(left, right, edge()).unwrap();
        if insertion {
            builder.add_edge(right, right, edge()).unwrap();
            builder.add_edge(left, outgoing, edge()).unwrap();
        } else {
            if reverse_tree_edge {
                builder.add_edge(right, left, edge()).unwrap();
            } else {
                builder.add_edge(left, right, edge()).unwrap();
            }
            builder.add_edge(right, outgoing, edge()).unwrap();
        }
        let diagram = builder.build().unwrap();
        if insertion {
            diagram
        } else {
            diagram
                .with_loop_momentum_edges(&[EdgeId(1), EdgeId(2)])
                .unwrap()
        }
    }

    #[test]
    fn dependent_graph_propagators_remain_available_for_partial_fractioning() {
        let diagram = diagram(true, false);
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
        let diagram = diagram(false, false);
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

    #[test]
    fn diagram_family_preserves_routing_mass_and_edge_orientation() {
        for reverse in [false, true] {
            let diagram = diagram(false, reverse);
            let family = diagram.propagator_family(&Kinematics::new()).unwrap();
            assert_eq!(
                diagram.loop_momentum_basis().dependent_externals,
                [EdgeId(4)]
            );
            assert_eq!(family.external_momenta().len(), 1);
            let [first, second] = family.loop_momenta() else {
                panic!("sunrise has two loop coordinates");
            };
            let external = &family.external_momenta()[0];
            let expected_momenta = [first.clone(), second.clone(), external - first - second];
            let mass = Atom::var(symbol!("UFO::M"));
            for (denominator, momentum) in family.denominators().iter().zip(expected_momenta) {
                let expected = family
                    .kinematics()
                    .scalar_product(&momentum, &momentum)
                    .unwrap()
                    - mass.pow(2);
                assert!((denominator - expected).expand().is_zero());
            }
            let mu = Atom::var(symbol!("family_test::mu"));
            let nu = Atom::var(symbol!("family_test::nu"));
            let raw = symbols::momentum().call((3, &mu))
                * (Atom::one() + symbols::momentum().call((4, &nu)));
            let routed_tree = symbols::external_momentum().call((0, &mu))
                - symbols::loop_momentum().call((0, &mu))
                - symbols::loop_momentum().call((1, &mu));
            let expected = Atom::num(if reverse { -1 } else { 1 })
                * routed_tree
                * (Atom::one() + symbols::external_momentum().call((0, &nu)));
            let routed = diagram.loop_momentum_basis().route_expression(&raw);
            assert!((&routed - expected).expand().is_zero());
            assert!(
                matches!(routed.as_view(), symbolica::atom::AtomView::Mul(product)
                if product.iter().count() == 2 && product.iter().all(|factor|
                    matches!(factor, symbolica::atom::AtomView::Add(_))))
            );
            let model = diagram.model();
            let mut particle = model.particle("phi").unwrap().clone();
            assert_eq!(particle.symbolic_mass_parameter(model).unwrap(), mass);
            particle.mass = model.parameter_id("ZERO").unwrap();
            assert!(particle.symbolic_mass_parameter(model).unwrap().is_zero());
            assert!(
                particle
                    .symbolic_mass_parameter(&feynkit_model::Model::empty("unrelated"))
                    .is_err()
            );
        }
    }
}
