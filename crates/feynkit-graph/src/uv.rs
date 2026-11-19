//! Local massive UV Taylor expansion, before integration or forest subtraction.

use std::collections::{BTreeMap, HashMap};

use idenso::tensor::SymbolicTensor;
use idenso::tensor::inference::TensorInferenceError;
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike, SubSetOps};
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::id::Replacement;
use symbolica::{function, symbol};

use crate::expressions::GraphExpressions;
use crate::{DiagramError, EdgeId, FeynmanDiagram, LoopMomentumBasis, symbols};

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
        let Some(numerator) = numerator else {
            return Ok(self
                .uv_expansion_tensor_of(subgraph, uv_mass, dimension, None, edge_powers)?
                .into_expression());
        };
        let Some(expansion) = self.uv_expansion(subgraph, uv_mass, dimension, edge_powers)? else {
            return Ok(Atom::Zero);
        };
        let routed = expansion.basis.route_expression(numerator);
        let scaled = UvExpansion::rescale_loops(&routed, expansion.scale);
        Ok(expansion
            .expand(&scaled)?
            .replace_multiple(&expansion.restore))
    }

    /// Expand an admitted numerator while retaining its tensor constructors.
    /// Routing, loop rescaling and restoration transport the original raw ports.
    /// Tensor constructors independent of the Taylor scale remain coefficients
    /// of the scalar series. Scale dependence inside a tensor constructor would
    /// require unsupported tensor derivatives and returns an explicit error.
    ///
    /// # Errors
    /// Returns invalid selections, tensor transport errors, unsupported tensor
    /// derivatives, or errors from the scalar series expansion.
    pub fn uv_expansion_tensor_of(
        &self,
        subgraph: &SuBitGraph,
        uv_mass: &Atom,
        dimension: i32,
        numerator: Option<&SymbolicTensor<PartialStructure>>,
        edge_powers: &BTreeMap<EdgeId, isize>,
    ) -> Result<SymbolicTensor<PartialStructure>, DiagramError> {
        let Some(expansion) = self.uv_expansion(subgraph, uv_mass, dimension, edge_powers)? else {
            return numerator
                .map_or_else(
                    || {
                        SymbolicTensor::from_produced_parts(
                            Atom::Zero,
                            PartialStructure::from_logical_slots([]),
                        )
                    },
                    |source| source.with_algebra_result(Atom::Zero),
                )
                .map_err(DiagramError::from);
        };
        let local_numerator;
        let numerator = match numerator {
            Some(numerator) => numerator,
            None => {
                local_numerator = self.numerator_of(
                    &expansion.selected,
                    &self.graph.empty_subgraph::<SuBitGraph>(),
                )?;
                &local_numerator
            }
        };
        let numerator = expansion.basis.route_tensor(numerator)?;
        let numerator = numerator.map_produced_expression(&HashMap::new(), |node| {
            Ok::<_, TensorInferenceError>(
                matches!(node, AtomView::Fun(function) if function.get_symbol() == symbols::loop_momentum())
                    .then(|| node / expansion.scale),
            )
        })?;
        let mut dependent_constructor = None;
        numerator.expression().visitor(&mut |node| {
            let AtomView::Fun(function) = node else {
                return true;
            };
            if function.get_symbol().is_scalar() {
                return false;
            }
            if function.get_symbol().has_tag(&SPENSO_TAG.tensor)
                && node.get_all_symbols(true).contains(&expansion.scale)
            {
                dependent_constructor = Some(function.get_symbol());
                return false;
            }
            dependent_constructor.is_none()
        });
        if let Some(head) = dependent_constructor {
            return Err(DiagramError::UvExpansion(format!(
                "Taylor scale occurs inside tensor constructor `{head}`; derivatives of tensor metadata are unsupported",
            )));
        }
        let expanded = numerator.with_materialized_expression(
            "UV Taylor expansion",
            expansion.expand(numerator.expression())?,
        )?;
        Ok(expanded.map_produced_expression(&HashMap::new(), |node| {
            if !matches!(node, AtomView::Fun(function)
                if function.get_symbol() == symbols::loop_momentum()
                    || function.get_symbol() == symbols::external_momentum())
            {
                return Ok(None);
            }
            let restored = node.replace_multiple(&expansion.restore);
            Ok::<_, TensorInferenceError>((restored.as_view() != node).then_some(restored))
        })?)
    }

    fn uv_expansion(
        &self,
        subgraph: &SuBitGraph,
        uv_mass: &Atom,
        dimension: i32,
        edge_powers: &BTreeMap<EdgeId, isize>,
    ) -> Result<Option<UvExpansion>, DiagramError> {
        if subgraph.size() != self.graph.n_hedges() || dimension <= 0 {
            return Err(DiagramError::UvExpansion(
                "expected a selection of this graph and a positive spacetime dimension".into(),
            ));
        }
        let selected = subgraph.intersection(&self.momentum_subgraph());
        if selected.is_empty() {
            return Ok(None);
        }
        let basis = self.momentum_basis_of(&selected)?;
        if basis.loop_edges.is_empty() {
            return Ok(None);
        }
        let scale = symbol!("feynkit_graph::uv_expansion_scale"; Scalar);
        let args = symbol!("feynkit_graph::uv_expansion_args___");
        let rescale = |expression: &Atom| {
            UvExpansion::rescale_loops(&basis.route_expression(expression), scale)
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
        // Restore the edge coordinates expected by diagram tensor reduction.
        let restore = [
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
        let measure = Atom::var(scale).pow(-i64::from(dimension) * basis.loop_edges.len() as i64)
            / denominator;
        Ok(Some(UvExpansion {
            selected,
            basis,
            scale,
            measure,
            restore,
        }))
    }
}

struct UvExpansion {
    selected: SuBitGraph,
    basis: LoopMomentumBasis,
    scale: Symbol,
    measure: Atom,
    restore: Vec<Replacement>,
}

impl UvExpansion {
    fn rescale_loops(routed: &Atom, scale: Symbol) -> Atom {
        let args = symbol!("feynkit_graph::uv_expansion_args___");
        let pattern = function!(symbols::loop_momentum(), args);
        routed.replace(pattern.to_pattern()).with(&pattern / scale)
    }

    fn expand(&self, scaled_numerator: &Atom) -> Result<Atom, DiagramError> {
        Ok((scaled_numerator * &self.measure)
            .series(self.scale, Atom::Zero, 0)
            .map_err(|error| DiagramError::UvExpansion(error.to_string()))?
            .to_atom()
            .replace(self.scale)
            .with(Atom::one()))
    }
}
