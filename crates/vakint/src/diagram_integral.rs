//! Checked ingress from FeynKit's stored graph routing, without graph rematching.

use std::collections::HashSet;

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::integer::Integer,
    id::Replacement,
};

use crate::{VakintError, VakintExpression, symbols::S};

/// Explicit scalar substitutions and spectator routing for a graph integral.
///
/// Physical propagators default to power one and auxiliary slots to zero.
/// Substitutions are simultaneous and the family must already contain the
/// substituted physical denominators. No graph symmetry factor or interaction
/// numerator is inserted: the supplied numerator is the complete requested one.
#[derive(Clone, Debug, Default)]
pub struct DiagramIntegralOptions {
    pub powers: Option<Vec<Integer>>,
    pub parameter_substitutions: Vec<(Atom, Atom)>,
    pub external_momenta: Vec<Atom>,
}

impl VakintExpression {
    /// Convert a routed, uncut vacuum graph using its ordered physical denominators.
    ///
    /// Extra family slots are allowed only with nonpositive powers. Integrated
    /// momenta retain family order, mapping to `k(1), k(2), ...`; the explicit
    /// spectator list maps to `p(1), p(2), ...`. Numerators must use native FeynKit
    /// scalar products in the family dimension. Unconverted tensors, undeclared
    /// spectators, on-shell loop assumptions and incompatible routings are refused.
    pub fn from_diagram(
        diagram: &FeynmanDiagram,
        family: &IntegralFamily,
        numerator: &Atom,
        options: &DiagramIntegralOptions,
    ) -> Result<Self, VakintError> {
        let invalid = |message: &str| VakintError::InvalidIntegralFormat(message.into());
        let native_error = |error: &dyn std::fmt::Display| invalid(&error.to_string());
        let basis = diagram.loop_momentum_basis();
        if !basis.external_edges.is_empty()
            || !diagram.cuts().is_empty()
            || !family.external_momenta().is_empty()
        {
            return Err(invalid(
                "only uncut vacuum graphs are supported; numerator spectators are explicit",
            ));
        }
        diagram.validate().map_err(|error| native_error(&error))?;
        let mut edges = diagram.edges().collect::<Vec<_>>();
        edges.sort_by_key(|(id, _, _)| *id);
        if edges.is_empty()
            || edges.iter().any(|(_, endpoints, edge)| {
                edge.is_dummy
                    || edge.external.is_some()
                    || endpoints.source.is_none()
                    || endpoints.target.is_none()
            })
        {
            return Err(invalid(
                "a vacuum topology needs non-dummy, non-dangling propagators",
            ));
        }
        if edges.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(invalid("duplicate graph edge IDs"));
        }

        // Start without assumptions: integrated momenta must remain unrestricted.
        let free = Kinematics::in_dimension(&family.kinematics().dimension().to_symbolic())
            .map_err(|error| native_error(&error))?;
        let physical = diagram
            .propagator_family(&free)
            .map_err(|error| native_error(&error))?;
        let loops = family.loop_momenta();
        if loops != physical.loop_momenta() || loops.len() != diagram.loop_count() {
            return Err(invalid(
                "family loop basis does not match the stored graph routing",
            ));
        }
        let momenta = loops
            .iter()
            .chain(&options.external_momenta)
            .cloned()
            .collect::<Vec<_>>();
        if momenta.iter().collect::<HashSet<_>>().len() != momenta.len() {
            return Err(invalid(
                "integrated and spectator momentum names must be distinct",
            ));
        }
        let heads = momenta
            .iter()
            .map(|momentum| match momentum.as_view() {
                AtomView::Var(variable) => Ok(variable.get_symbol()),
                AtomView::Fun(function) => Ok(function.get_symbol()),
                _ => Err(invalid(
                    "expected a momentum name, not a sum or non-momentum descriptor",
                )),
            })
            .collect::<Result<HashSet<_>, _>>()?;
        let scalar_symbols = |expression: &Atom| {
            expression
                .get_all_symbols(true)
                .iter()
                .all(|symbol| !heads.contains(symbol) && !symbol.get_name().starts_with("spenso::"))
        };
        if options
            .parameter_substitutions
            .iter()
            .any(|(source, target)| {
                !matches!(source.as_view(), AtomView::Var(_))
                    || !scalar_symbols(source)
                    || !scalar_symbols(target)
            })
        {
            return Err(invalid(
                "parameter substitutions must be scalar, not momentum/representation rewrites",
            ));
        }
        if options
            .parameter_substitutions
            .iter()
            .map(|(source, _)| source)
            .collect::<HashSet<_>>()
            .len()
            != options.parameter_substitutions.len()
        {
            return Err(invalid("parameter substitution sources must be distinct"));
        }
        let substitutions = options
            .parameter_substitutions
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let substitute = |expression: &Atom| expression.replace_multiple(&substitutions);
        let denominators = family.denominators();
        if denominators.len() < edges.len()
            || physical.denominators().len() != edges.len()
            || denominators
                .iter()
                .zip(physical.denominators())
                .any(|(actual, wanted)| !(actual - substitute(wanted)).expand().is_zero())
        {
            return Err(invalid(
                "family physical denominator order/sign/masses do not match the graph",
            ));
        }
        let default_powers;
        let powers = match options.powers.as_deref() {
            Some(powers) => powers,
            None => {
                default_powers = (0..denominators.len())
                    .map(|index| Integer::from(i64::from(index < edges.len())))
                    .collect::<Vec<_>>();
                &default_powers
            }
        };
        if powers.len() != denominators.len() {
            return Err(invalid(
                "powers must be one integer (not bool) per family denominator",
            ));
        }
        if powers[..edges.len()]
            .iter()
            .any(|power| power <= &Integer::from(0))
            || powers[edges.len()..]
                .iter()
                .any(|power| power > &Integer::from(0))
        {
            return Err(invalid(
                "physical powers must be positive; auxiliary powers must be nonpositive",
            ));
        }

        // Existing heads install Vakint's tensor tags as well as its dot symmetry.
        let mapped = (0..loops.len())
            .map(|i| S.k.call(i + 1))
            .chain((0..options.external_momenta.len()).map(|i| S.p.call(i + 1)))
            .collect::<Vec<_>>();
        let mut topology = Atom::one();
        for (index, ((id, endpoints, _), power)) in edges.iter().zip(powers).enumerate() {
            let signature = basis
                .edge_signatures
                .get(id)
                .ok_or_else(|| invalid("stored graph routing is missing a propagator"))?;
            if signature.loops.len() != loops.len()
                || signature
                    .external
                    .integer_coefficients()
                    .iter()
                    .any(|value| *value != 0)
            {
                return Err(invalid(
                    "stored graph routing is not the declared vacuum loop basis",
                ));
            }
            let momentum = signature
                .loops
                .apply(&mapped[..loops.len()])
                .map_err(|error| native_error(&error))?
                .unwrap_or_default();
            let original = signature
                .loops
                .apply(loops)
                .map_err(|error| native_error(&error))?
                .unwrap_or_default();
            let square = physical
                .kinematics()
                .scalar_product(&original, &original)
                .map_err(|error| native_error(&error))?;
            let mass_squared = substitute(&(square - &physical.denominators()[index]).expand());
            if !scalar_symbols(&mass_squared) {
                return Err(invalid(
                    "graph denominator does not have a scalar mass squared",
                ));
            }
            let source = endpoints
                .source
                .ok_or_else(|| invalid("propagator has an incomplete incidence"))?;
            let target = endpoints
                .target
                .ok_or_else(|| invalid("propagator has an incomplete incidence"))?;
            topology *= S.prop.call((
                id.0,
                S.edge.call((source.0, target.0)),
                momentum,
                mass_squared,
                Atom::num(power.clone()),
            ));
        }

        let mut scalar = substitute(numerator);
        for (denominator, power) in denominators[edges.len()..]
            .iter()
            .zip(&powers[edges.len()..])
        {
            scalar *= denominator.pow(Atom::num(-power.clone()));
        }
        let free = free
            .with_momenta(momenta.iter().cloned())
            .map_err(|error| native_error(&error))?;
        let mut products = Vec::new();
        for (i, left) in momenta.iter().enumerate() {
            for (j, right) in momenta.iter().enumerate().skip(i) {
                let product = free
                    .scalar_product(left, right)
                    .map_err(|error| native_error(&error))?;
                products.push(Replacement::new(
                    product,
                    S.dot.call((&mapped[i], &mapped[j])),
                ));
            }
        }
        scalar = scalar.expand().replace_multiple(&products);
        if !scalar_symbols(&scalar) {
            return Err(invalid(
                "numerator contains unconverted momentum/tensor notation; use declared scalar products",
            ));
        }
        Self::try_from(scalar * S.topo.call(topology))
    }
}

#[cfg(test)]
mod tests;
