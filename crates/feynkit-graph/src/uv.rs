//! Local massive UV Taylor expansion, before integration or forest subtraction.

use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    symbol,
};

use crate::{DiagramError, EdgeId, FeynmanDiagram, expressions::GraphExpressions, symbols};

impl FeynmanDiagram {
    /// Expand the selected local integrand through its UV degree of divergence.
    ///
    /// Only the selected region's independent loop momenta are hard. For each
    /// loop-dependent propagator, write q = k + p and use
    /// t² / [(k + t p)² - t² m² - (1 - t²) uv_mass²]. The numerator is
    /// evaluated at k/t + p. Keeping powers through t⁰ after multiplication
    /// by t^(-dimension * loops) retains all superficially divergent terms.
    /// The bookkeeping scale is then set to one.
    ///
    /// The result uses edge momenta and tagged `denom` propagators, as does
    /// `denominator_of`. Negate it for an additive local counterterm. A supplied
    /// numerator replaces the selected local numerator and must use edge momenta.
    /// Projectors, overall factors and numerator prefactors are not added.
    /// Empty and tree regions return zero. This is one simultaneous UV limit,
    /// not a recursive forest subtraction or an integrated pole counterterm.
    pub fn uv_expansion_of(
        &self,
        subgraph: &SuBitGraph,
        uv_mass: &Atom,
        dimension: i32,
        numerator: Option<&Atom>,
    ) -> Result<Atom, DiagramError> {
        if subgraph.size() != self.graph.n_hedges() || dimension <= 0 {
            return Err(DiagramError::UvExpansion(
                "expected a selection of this graph and a positive spacetime dimension".into(),
            ));
        }
        let selected = subgraph.intersection(&self.momentum_subgraph());
        if selected.is_empty() {
            return Ok(Atom::Zero);
        }
        let basis = self.momentum_basis_of(&selected)?;
        if basis.loop_edges.is_empty() {
            return Ok(Atom::Zero);
        }
        let scale = symbol!("feynkit_graph::uv_expansion_scale"; Scalar);
        let args = symbol!("feynkit_graph::uv_expansion_args___");
        let loop_pattern = function!(symbols::loop_momentum(), args);
        let rescale = |expression: &Atom| {
            basis
                .route_expression(expression)
                .replace(loop_pattern.to_pattern())
                .with(&loop_pattern / scale)
        };
        let local_numerator;
        let numerator = match numerator {
            Some(numerator) => numerator,
            None => {
                local_numerator =
                    self.numerator_of(&selected, &self.graph.empty_subgraph::<SuBitGraph>());
                &local_numerator
            }
        };
        let metric = Minkowski {}.new_rep(symbols::dimension());
        let scale_squared = Atom::var(scale).pow(2);
        let uv_mass_squared = uv_mass.pow(2);
        let denominator = self.graph.denominator_of(
            &selected.intersection(&self.internal_subgraph()),
            |edge, data| -> Result<Atom, DiagramError> {
                let mass = self
                    .model
                    .particle_by_id(data.particle)?
                    .symbolic_mass(&self.model)
                    .replace(symbol!("UFO::ZERO"))
                    .with(Atom::Zero);
                let edge_momentum = symbols::momentum().call(edge.0);
                let momentum = rescale(&edge_momentum);
                let hard = (&momentum * scale).expand().replace(scale).with(Atom::Zero);
                let quadratic =
                    rescale(&metric.inner_product(&edge_momentum, &edge_momentum)) - mass.pow(2);
                if hard == Atom::Zero {
                    // Bridges carry only soft momentum and remain spectator propagators.
                    return Ok(symbols::denominator().call_args([
                        Atom::num(edge.0),
                        momentum,
                        mass.pow(2),
                        quadratic,
                    ]));
                }
                let quadratic = (quadratic * &scale_squared + &uv_mass_squared * &scale_squared
                    - &uv_mass_squared)
                    .expand();
                Ok(symbols::denominator().call_args([
                    Atom::num(edge.0),
                    hard,
                    uv_mass_squared.clone(),
                    quadratic,
                ]) / &scale_squared)
            },
            |_, _| 1,
        )?;
        let measured = rescale(numerator) / denominator
            * Atom::var(scale).pow(-i64::from(dimension) * basis.loop_edges.len() as i64);
        let expanded = measured
            .series(scale, Atom::Zero, 0)
            .map_err(|error| DiagramError::UvExpansion(error.to_string()))?
            .to_atom()
            .replace(scale)
            .with(Atom::one());
        // Restore the edge coordinates expected by diagram tensor reduction.
        let replacements = [
            (symbols::loop_momentum(), &basis.loop_edges),
            (symbols::external_momentum(), &basis.external_edges),
        ]
        .into_iter()
        .flat_map(|(head, edges)| {
            edges.iter().enumerate().map(move |(index, EdgeId(edge))| {
                Replacement::new(
                    function!(head, index, args).to_pattern(),
                    function!(symbols::momentum(), *edge, args).to_pattern(),
                )
            })
        })
        .collect::<Vec<_>>();
        Ok(expanded.replace_multiple(&replacements))
    }
}
