//! Checked ingress from FeynKit's stored graph routing, without graph rematching.

use std::collections::{BTreeMap, HashSet};

use feynkit_graph::{FeynmanDiagram, IntegralFamily};
use feynkit_kinematics::Kinematics;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::integer::Integer,
    id::Replacement,
};

use crate::{VakintError, VakintExpression, symbols::S};

mod constant_factors;

/// Explicit scalar substitutions and spectator routing for a graph integral.
///
/// Physical propagators default to their graph annotations and auxiliary slots
/// to zero. Explicit powers override those annotations.
/// Substitutions are simultaneous and the family must already contain the
/// substituted physical denominators. No graph symmetry factor or interaction
/// numerator is inserted: the supplied numerator is the complete requested one.
#[derive(Clone, Debug, Default)]
pub struct DiagramIntegralOptions {
    pub powers: Option<Vec<Integer>>,
    /// Zero-based family slots in sorted graph-edge order. Defaults to the prefix.
    /// Zero or negative powers pinch the edge; negative powers also enter the numerator.
    pub physical_slots: Option<Vec<usize>>,
    pub parameter_substitutions: Vec<(Atom, Atom)>,
    pub external_momenta: Vec<Atom>,
    /// The name of the external momentum of a propagator-type graph (two
    /// external legs), a bare `p(i)`; `p(1)` by default.
    pub external_momentum: Option<Atom>,
}

impl VakintExpression {
    /// Convert a routed, uncut vacuum graph using its ordered physical denominators.
    ///
    /// Nonphysical family slots are allowed only with nonpositive powers. Physical
    /// zero/negative powers contract their graph edges in the original routing.
    /// Massive zero-flow bridges are first extracted as native scalar factors;
    /// contracting their incidence does not change the integrated loop basis.
    /// Their positive powers must fit the native u32 exponent representation;
    /// this does not promise bounded-size materialization of numeric scalars.
    /// Contractions that remove an independent loop are rejected. Integrated
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
        if basis.external_edges.len() == 2 && diagram.cuts().is_empty() {
            return Self::from_propagator_diagram(diagram, family, numerator, options);
        }
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
        // Vacuum propagators have two internal endpoints. Reject dangling
        // edges before projecting the validated incidence graph.
        if edges.is_empty()
            || edges.iter().any(|(_, endpoints, edge)| {
                edge.external.is_some() || endpoints.source.is_none() || endpoints.target.is_none()
            })
        {
            return Err(invalid("a vacuum topology needs internal propagators"));
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
        let default_slots = (0..edges.len()).collect::<Vec<_>>();
        let physical_slots = options.physical_slots.as_deref().unwrap_or(&default_slots);
        let physical_set = physical_slots.iter().copied().collect::<HashSet<_>>();
        if physical_slots.len() != edges.len()
            || physical_set.len() != physical_slots.len()
            || physical_slots
                .iter()
                .any(|slot| *slot >= denominators.len())
        {
            return Err(invalid(
                "physical slots must be distinct valid family indices, one per graph edge",
            ));
        }
        if denominators.len() < edges.len()
            || physical.denominators().len() != edges.len()
            || physical_slots
                .iter()
                .zip(physical.denominators())
                .any(|(slot, wanted)| {
                    !(&denominators[*slot] - substitute(wanted))
                        .expand()
                        .is_zero()
                })
        {
            return Err(invalid(
                "family physical denominator order/sign/masses do not match the graph",
            ));
        }
        let default_powers;
        let powers = match options.powers.as_deref() {
            Some(powers) => powers,
            None => {
                default_powers = {
                    let mut values = vec![Integer::from(0); denominators.len()];
                    for ((_, _, edge), slot) in edges.iter().zip(physical_slots) {
                        values[*slot] = Integer::from(edge.propagator_power);
                    }
                    values
                };
                &default_powers
            }
        };
        if powers.len() != denominators.len() {
            return Err(invalid(
                "powers must be one integer (not bool) per family denominator",
            ));
        }
        if powers
            .iter()
            .enumerate()
            .any(|(slot, power)| !physical_set.contains(&slot) && power > &Integer::from(0))
        {
            return Err(invalid("auxiliary powers must be nonpositive"));
        }

        let mut constant_slots = Vec::new();
        for ((id, _, _), slot) in edges.iter().zip(physical_slots) {
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
            if signature
                .loops
                .integer_coefficients()
                .iter()
                .all(|value| *value == 0)
            {
                constant_slots.push(*slot);
            }
        }
        let (constant_factor, mut effective_powers) =
            constant_factors::extract(family, powers, &constant_slots)?;

        // Native scalar diagrams may contain repeated inverse propagators
        // (for example the two edges of a bubble with tadpole attachments).
        // Their exact integer powers add before graph contraction; otherwise
        // one integral class would spuriously depend on that graph realization.
        let mut representatives = BTreeMap::<Atom, usize>::new();
        for slot in physical_slots {
            let denominator = denominators[*slot].expand();
            if let Some(representative) = representatives.get(&denominator) {
                let power = effective_powers[*slot].clone();
                effective_powers[*representative] += power;
                effective_powers[*slot] = Integer::from(0);
            } else {
                representatives.insert(denominator, *slot);
            }
        }
        let powers = &effective_powers;

        // Contract incidence only. The original stored routing is also used for
        // every numerator factor, including a cancelled physical denominator.
        let mut vertices = BTreeMap::new();
        for (_, endpoints, _) in &edges {
            for vertex in [
                endpoints.source.expect("validated vacuum endpoint").0,
                endpoints.target.expect("validated vacuum endpoint").0,
            ] {
                vertices.insert(vertex, vertex);
            }
        }
        for ((_, endpoints, _), slot) in edges.iter().zip(physical_slots) {
            if powers[*slot] > 0 {
                continue;
            }
            let source = vertices[&endpoints.source.expect("validated vacuum endpoint").0];
            let target = vertices[&endpoints.target.expect("validated vacuum endpoint").0];
            if source == target {
                return Err(invalid(
                    "a pinch removes an independent vacuum loop; its scaleless integration is not represented by this graph",
                ));
            }
            let (keep, remove) = (source.min(target), source.max(target));
            for vertex in vertices.values_mut() {
                if *vertex == remove {
                    *vertex = keep;
                }
            }
        }

        // Existing heads install Vakint's tensor tags as well as its dot symmetry.
        let mapped = (0..loops.len())
            .map(|i| S.k.call(i + 1))
            .chain((0..options.external_momenta.len()).map(|i| S.p.call(i + 1)))
            .collect::<Vec<_>>();
        let mut topology = Atom::one();
        for (index, ((id, endpoints, _), slot)) in edges.iter().zip(physical_slots).enumerate() {
            let power = &powers[*slot];
            if power <= &Integer::from(0) {
                continue;
            }
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
            let source = endpoints.source.expect("validated vacuum endpoint");
            let target = endpoints.target.expect("validated vacuum endpoint");
            topology *= S.prop.call((
                id.0,
                S.edge.call((vertices[&source.0], vertices[&target.0])),
                momentum,
                mass_squared,
                Atom::num(power.clone()),
            ));
        }

        let mut scalar = substitute(numerator) * constant_factor;
        for (denominator, power) in denominators.iter().zip(powers) {
            if power < &Integer::from(0) {
                scalar *= denominator.pow(Atom::num(-power.clone()));
            }
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
        // Scalar-product replacements also traverse sums, products and powers.
        // Keeping the complete numerator factored avoids materializing a large
        // polynomial merely to rename its already contracted scalar products.
        scalar = scalar.replace_multiple(&products);
        if !scalar_symbols(&scalar) {
            return Err(invalid(
                "numerator contains unconverted momentum/tensor notation; use declared scalar products",
            ));
        }
        Self::try_from(scalar * S.topo.call(topology))
    }

    /// Convert a routed, uncut massless propagator-type graph: exactly two
    /// external legs (docs/p_integrals.md, section 3).
    ///
    /// The external momentum `p(1)` (or `options.external_momentum`) flows into
    /// the graph through the leg with the smallest edge id; the stored routing
    /// supplies its sign, whichever leg FeynKit made dependent. Every graph line
    /// becomes a `prop` with its power, also when it is zero or negative: the
    /// ordinary p-integral normalization then contracts, merges and moves
    /// them to the numerator exactly as for a written Atom. Negative powers of
    /// auxiliary family slots enter the numerator.
    fn from_propagator_diagram(
        diagram: &FeynmanDiagram,
        family: &IntegralFamily,
        numerator: &Atom,
        options: &DiagramIntegralOptions,
    ) -> Result<Self, VakintError> {
        let invalid = |message: &str| VakintError::InvalidIntegralFormat(message.into());
        let native_error = |error: &dyn std::fmt::Display| invalid(&error.to_string());
        if !options.external_momenta.is_empty() {
            return Err(invalid(
                "p-graphs take no spectator momenta: external_momenta must be empty",
            ));
        }
        let head = options
            .external_momentum
            .clone()
            .unwrap_or_else(|| S.p.call(1));
        if !matches!(head.as_view(), AtomView::Fun(function)
            if function.get_symbol() == S.p && function.get_nargs() == 1
                && function.iter().next().and_then(|id| i64::try_from(id).ok()).is_some())
        {
            return Err(invalid("external_momentum must be a bare p(i)"));
        }
        diagram.validate().map_err(|error| native_error(&error))?;
        let basis = diagram.loop_momentum_basis();
        let mut edges = diagram.edges().collect::<Vec<_>>();
        edges.sort_by_key(|(id, _, _)| *id);
        let (legs, mut lines): (Vec<_>, Vec<_>) = edges.into_iter().partition(|(_, ends, edge)| {
            edge.external.is_some() || ends.source.is_none() || ends.target.is_none()
        });
        if legs.len() != 2 || lines.is_empty() {
            return Err(invalid(
                "a p-graph has exactly two external legs and internal propagators",
            ));
        }
        lines.sort_by_key(|(id, _, _)| *id);

        let free = Kinematics::in_dimension(&family.kinematics().dimension().to_symbolic())
            .map_err(|error| native_error(&error))?;
        let physical = diagram
            .propagator_family(&free)
            .map_err(|error| native_error(&error))?;
        let loops = family.loop_momenta();
        if loops != physical.loop_momenta()
            || loops.len() != diagram.loop_count()
            || family.external_momenta() != physical.external_momenta()
            || physical.external_momenta().len() != 1
        {
            return Err(invalid(
                "family loop and external basis does not match the stored graph routing",
            ));
        }
        // The external basis is indexed by leg position; the dependent leg's
        // coordinate never occurs.
        let independent = &physical.external_momenta()[0];
        let externals = basis
            .external_edges
            .iter()
            .map(|edge| {
                if basis.dependent_externals.contains(edge) {
                    Atom::Zero
                } else {
                    independent.clone()
                }
            })
            .collect::<Vec<_>>();
        let signature = |id| {
            basis
                .edge_signatures
                .get(id)
                .ok_or_else(|| invalid("stored graph routing is missing an edge"))
        };
        // p(1) is the momentum flowing into the graph through the first leg.
        let (first, first_ends, _) = &legs[0];
        let along = signature(first)?
            .external
            .apply(&externals)
            .map_err(|error| native_error(&error))?
            .unwrap_or_default();
        let into = if first_ends.source.is_none() {
            along
        } else {
            -along
        };
        // into = c·P with c = ±1, hence P = c·p(1).
        let orientation = if (into.clone() - independent).expand().is_zero() {
            Atom::num(1)
        } else if (into + independent).expand().is_zero() {
            Atom::num(-1)
        } else {
            return Err(invalid(
                "the first external leg does not carry the external momentum",
            ));
        };
        let scaled = &orientation * &head;

        let parameter = |expression: &Atom| {
            expression.get_all_symbols(true).iter().all(|symbol| {
                !symbol.get_name().starts_with("spenso::")
                    && loops
                        .iter()
                        .chain([independent])
                        .all(|momentum| !momentum.get_all_symbols(true).contains(symbol))
            })
        };
        if options
            .parameter_substitutions
            .iter()
            .any(|(source, target)| {
                !matches!(source.as_view(), AtomView::Var(_))
                    || !parameter(source)
                    || !parameter(target)
            })
        {
            return Err(invalid(
                "parameter substitutions must be scalar, not momentum/representation rewrites",
            ));
        }
        let substitutions = options
            .parameter_substitutions
            .iter()
            .map(|(source, target)| Replacement::new(source.clone(), target.clone()))
            .collect::<Vec<_>>();
        let substitute = |expression: &Atom| expression.replace_multiple(&substitutions);
        let denominators = family.denominators();
        let default_slots = (0..lines.len()).collect::<Vec<_>>();
        let physical_slots = options.physical_slots.as_deref().unwrap_or(&default_slots);
        let physical_set = physical_slots.iter().copied().collect::<HashSet<_>>();
        if physical_slots.len() != lines.len()
            || physical_set.len() != physical_slots.len()
            || physical_slots
                .iter()
                .any(|slot| *slot >= denominators.len())
            || physical.denominators().len() != lines.len()
            || physical_slots
                .iter()
                .zip(physical.denominators())
                .any(|(slot, wanted)| {
                    !(&denominators[*slot] - substitute(wanted))
                        .expand()
                        .is_zero()
                })
        {
            return Err(invalid(
                "family physical denominator order/sign/masses do not match the graph",
            ));
        }
        let powers = match options.powers.as_deref() {
            Some(powers) => powers.to_vec(),
            None => {
                let mut values = vec![Integer::from(0); denominators.len()];
                for ((_, _, edge), slot) in lines.iter().zip(physical_slots) {
                    values[*slot] = Integer::from(edge.propagator_power);
                }
                values
            }
        };
        if powers.len() != denominators.len() {
            return Err(invalid(
                "powers must be one integer (not bool) per family denominator",
            ));
        }
        if powers
            .iter()
            .enumerate()
            .any(|(slot, power)| !physical_set.contains(&slot) && power > &Integer::from(0))
        {
            return Err(invalid("auxiliary powers must be nonpositive"));
        }

        let mapped = (0..loops.len())
            .map(|i| S.k.call(i + 1))
            .collect::<Vec<_>>();
        let mut topology = Atom::one();
        for (index, ((id, ends, _), slot)) in lines.iter().zip(physical_slots).enumerate() {
            let line = signature(id)?;
            let momentum = line
                .apply(&mapped, &externals)
                .map_err(|error| native_error(&error))?
                .unwrap_or_default()
                .replace(independent.to_pattern())
                .with(scaled.to_pattern());
            let original = line
                .apply(loops, &externals)
                .map_err(|error| native_error(&error))?
                .unwrap_or_default();
            let square = physical
                .kinematics()
                .scalar_product(&original, &original)
                .map_err(|error| native_error(&error))?;
            if !substitute(&(square - &physical.denominators()[index]).expand())
                .expand()
                .is_zero()
            {
                return Err(invalid("p-graph propagators must be massless"));
            }
            let (Some(source), Some(target)) = (ends.source, ends.target) else {
                unreachable!("internal lines have two endpoints");
            };
            topology *= S.prop.call((
                id.0,
                S.edge.call((source.0, target.0)),
                momentum,
                Atom::Zero,
                Atom::num(powers[*slot].clone()),
            ));
        }

        let mut scalar = substitute(numerator);
        for (slot, (denominator, power)) in denominators.iter().zip(&powers).enumerate() {
            if !physical_set.contains(&slot) && power < &Integer::from(0) {
                scalar *= denominator.pow(Atom::num(-power.clone()));
            }
        }
        let momenta = loops
            .iter()
            .chain([independent])
            .cloned()
            .collect::<Vec<_>>();
        let targets = mapped.iter().chain([&scaled]).cloned().collect::<Vec<_>>();
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
                    S.dot.call((&targets[i], &targets[j])),
                ));
            }
        }
        scalar = scalar.replace_multiple(&products);
        if !parameter(&scalar) {
            return Err(invalid(
                "numerator contains unconverted momentum/tensor notation; use declared scalar products",
            ));
        }
        Self::try_from(scalar * S.topo.call(topology))
    }
}

#[cfg(test)]
mod tests;
