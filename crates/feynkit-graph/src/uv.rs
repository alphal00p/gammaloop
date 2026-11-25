//! Local massive UV Taylor expansion, before integration or forest subtraction.

use std::collections::BTreeMap;
use std::convert::Infallible;

use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    symbol,
};

use crate::{DOD, DiagramError, EdgeId, FeynmanDiagram, expressions::GraphExpressions, symbols};

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
    /// Signed `edge_powers` override internal propagator powers, defaulting to one;
    /// zero omits a denominator and negative powers move it to the numerator.
    /// Powers do not change the selected region or its independent loop momenta.
    /// Entries outside the selected internal edges are ignored, as in `denominator_of`.
    /// Projectors, overall factors and numerator prefactors are not added.
    /// Empty and tree regions return zero. This is one simultaneous UV limit,
    /// not a recursive forest subtraction or an integrated pole counterterm.
    pub fn uv_expansion_of(
        &self,
        subgraph: &SuBitGraph,
        uv_mass: &Atom,
        dimension: i32,
        numerator: Option<&Atom>,
        edge_powers: &BTreeMap<EdgeId, isize>,
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
                    self.numerator_of(&selected, &self.graph.empty_subgraph::<SuBitGraph>())?;
                local_numerator.expression()
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
            |edge, _| edge_powers.get(&EdgeId(edge.0)).copied().unwrap_or(1),
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

impl FeynmanDiagram {
    /// Classifying UV degree of divergence of a selected region in `dimension`
    /// spacetime dimensions, by rescaling only its hard loop momenta.
    ///
    /// The region's local numerator (as [`Self::numerator_of`], with dummy edges
    /// excluded) is routed in the region's own momentum basis, whose external
    /// coordinates are the momenta crossing its boundary. Every independent loop
    /// momentum `k` becomes `k/t`; external momenta stay fixed. Each internal
    /// propagator of the region contributes a massless `1/q²`, and the measure
    /// contributes `t^(-dimension * loops)`. The result is the leading power of
    /// `1/t`: zero is logarithmic, positive is power divergent and negative is
    /// convergent. Unlike [`Self::superficial_degree_of_divergence_of`], numerator
    /// factors carrying only external momenta do not count. Masses do not change
    /// the degree.
    ///
    /// All lines are routed to common loop coordinates before rescaling, so hard
    /// parts of distinct lines can cancel in a local factor: a vertex factor
    /// `(q₂⁰ + q₃⁰)²` on two parallel bubble lines oriented alike counts as zero.
    /// GammaLoop's runtime counting tags each line's hard momentum, does not
    /// cancel them, and gives two there.
    ///
    /// A local numerator that vanishes has no leading power; it is assigned degree
    /// zero, so such a region counts as logarithmically divergent with a vanishing
    /// integrand. Vanishing is detected when the routed numerator is zero or one
    /// routed edge or vertex factor expands to zero. Indices are never contracted:
    /// a numerator that vanishes only by contraction across factors, such as
    /// `g_μν (q₀^μ p^ν − q₀^ν p^μ)` split over two vertices, is not detected and can
    /// come out at any degree, including a negative one, so such a region may be
    /// classified as convergent although its integrand is zero. Contracting is not
    /// done because it is too costly to run on every candidate region.
    ///
    /// # Errors
    ///
    /// Fails when the region cannot be routed or its leading power is not an
    /// integer.
    pub fn uv_degree_of_divergence_of(
        &self,
        subgraph: &SuBitGraph,
        dimension: i32,
    ) -> Result<i32, DiagramError> {
        let selected = subgraph.intersection(&self.momentum_subgraph());
        let basis = self.momentum_basis_of(&selected)?;
        let scale = symbol!("feynkit_graph::uv_degree_scale"; Scalar);
        let args = symbol!("feynkit_graph::uv_degree_args___");
        let loop_pattern = function!(symbols::loop_momentum(), args);
        let rescale = |expression: &Atom| {
            basis
                .route_expression(expression)
                .replace(loop_pattern.to_pattern())
                .with(&loop_pattern / scale)
        };
        let numerator = self.numerator_of(&selected, &self.graph.empty_subgraph::<SuBitGraph>())?;
        let metric = Minkowski {}.new_rep(symbols::dimension());
        let routed = rescale(numerator.expression());
        // The local numerator is a product of edge and vertex factors, so it
        // vanishes when one routed factor expands to zero. Expanding the factors
        // separately exposes exact cancellations, such as an identity of
        // squares, without expanding their product.
        let vanishing_factor = self
            .graph
            .try_fold_numerator(
                &selected,
                &self.graph.empty_subgraph::<SuBitGraph>(),
                false,
                |_, vertex| {
                    basis
                        .route_expression(vertex.numerator.expression())
                        .expand()
                        .is_zero()
                },
                |edge| {
                    basis
                        .route_expression(edge.numerator.expression())
                        .expand()
                        .is_zero()
                },
                |edge| edge.is_dummy,
                |vanishes, factor| Ok::<_, Infallible>(vanishes || factor),
            )
            .unwrap_or_else(|never| match never {});
        if routed.is_zero() || vanishing_factor {
            return Ok(0);
        }
        let mut integrand =
            routed * Atom::var(scale).pow(-i64::from(dimension) * basis.loop_edges.len() as i64);
        for (pair, edge, _) in self
            .graph
            .iter_edges_of(&selected.intersection(&self.internal_subgraph()))
        {
            if pair.is_paired() {
                let momentum = symbols::momentum().call(edge.0);
                integrand /= rescale(&metric.inner_product(&momentum, &momentum));
            }
        }
        integrand.trailing_exponent(scale)
    }
}
