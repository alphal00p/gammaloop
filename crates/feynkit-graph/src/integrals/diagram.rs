use std::collections::BTreeMap;

use feynkit_kinematics::Kinematics;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{AtomCore, AtomView};

use super::{IntegralFamily, IntegralFamilyError};
use crate::{DiagramError, FeynmanDiagram, symbols};

impl FeynmanDiagram {
    /// Build an integral family from the diagram's internal quadratic propagators.
    ///
    /// Denominators follow ascending internal edge IDs, retaining repeated
    /// propagators and bridges. The stored loop routing supplies momentum names;
    /// dependent external coordinates are eliminated by that same routing.
    /// External carriers and dummy edges are excluded by the shared denominator
    /// builder. Symbolic masses follow its conventions, including UFO ZERO.
    /// Widths, prescriptions and custom UFO denominator formulas are not inferred.
    /// Tree diagrams have no loop-integral family and return an error.
    pub fn integral_family(&self, kinematics: &Kinematics) -> Result<IntegralFamily, DiagramError> {
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
