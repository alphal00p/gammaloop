use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Deref,
};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::eyre::{Context, Result, eyre};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use spenso::{
    algebra::complex::Complex,
    network::{ExecutionResult, SequentialRef, SmallestDegree},
    structure::abstract_index::AIND_SYMBOLS,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    domains::float::SingleFloat,
    evaluate::FunctionMap,
    function,
    prelude::*,
    symbol,
};

use crate::{
    GammaLoopContext,
    cff::esurface::EsurfaceID,
    graph::{Graph, LoopMomentumBasis},
    integrands::evaluation::EvaluationMetaData,
    momentum::{
        ThreeMomentum,
        sample::{ExternalFourMomenta, ExternalThreeMomenta, LoopMomenta, MomentumSample},
    },
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::{
        EvaluatorSettings, ThresholdCountertermEvaluatorRegistration, ThresholdCountertermVariantId,
    },
    utils::{F, FUN_LIB, FloatLike, GS, TENSORLIB},
};

use super::evaluators::{EvaluatorBackendPolicy, GenericEvaluator, evaluate_evaluator_single};

/// The global kinematic view at which a multiplier input is evaluated.
///
/// `Star` is assembled for the expanded algebraic component being evaluated. Participating
/// counterterm sides use their selected roots; an absent side has no independently solved root,
/// so its spectator parent-LMB coordinates undergo the identity transformation. Its original
/// integrand factor is still present. The complete sample is rebuilt from the participating roots,
/// so cross-dependent edge momenta or E-surfaces need not equal their base-sample values. An
/// iterated component carries both selected roots without exposing a user-level Cartesian product.
#[derive(Clone, Copy, Debug, Encode, Decode, Eq, Ord, PartialEq, PartialOrd)]
pub enum ThresholdMultiplierPoint {
    Effective,
    Star,
}

impl ThresholdMultiplierPoint {
    fn symbol(self) -> Symbol {
        match self {
            Self::Effective => symbol!("effective"),
            Self::Star => symbol!("star"),
        }
    }

    fn atom(self) -> Atom {
        Atom::var(self.symbol())
    }

    fn parse(atom: AtomView<'_>) -> Result<Self> {
        let AtomView::Var(variable) = atom else {
            return Err(eyre!(
                "expected multiplier frame `effective` or `star`, got `{atom}`"
            ));
        };
        if variable.get_symbol() == Self::Effective.symbol() {
            Ok(Self::Effective)
        } else if variable.get_symbol() == Self::Star.symbol() {
            Ok(Self::Star)
        } else {
            Err(eyre!(
                "unknown multiplier frame `{atom}`; expected `effective` or `star`"
            ))
        }
    }
}

/// One entry in the fixed cut-owned multiplier input layout.
#[derive(Clone, Debug, Encode, Decode, Eq, PartialEq)]
pub enum ThresholdMultiplierInput {
    ModelParameter {
        index: usize,
    },
    AdditionalParameter {
        index: usize,
    },
    ExternalMomentum {
        position: usize,
        component: usize,
    },
    EdgeMomentum {
        point: ThresholdMultiplierPoint,
        edge: usize,
        component: usize,
    },
    Esurface {
        point: ThresholdMultiplierPoint,
        esurface: usize,
    },
}

/// The identity and runtime binding data for one E-surface equation.
#[derive(Clone, Debug, Encode, Decode, Eq, Ord, PartialEq, PartialOrd)]
pub struct ThresholdMultiplierEsurface {
    pub edges: Vec<usize>,
    pub external_shift: Vec<(usize, i64)>,
}

/// A deterministic parameter layout shared by all multiplier evaluators under one cut. Model
/// parameters come first, followed by the graph's DOT `params` in declaration order.
#[derive(Clone, Debug, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ThresholdMultiplierLayout {
    model_parameters: Vec<Atom>,
    additional_parameters: Vec<Atom>,
    external_count: usize,
    edges: Vec<usize>,
    esurfaces: Vec<ThresholdMultiplierEsurface>,
    inputs: Vec<ThresholdMultiplierInput>,
    parameters: Vec<Atom>,
}

impl ThresholdMultiplierLayout {
    pub fn new(
        model_parameters: Vec<Atom>,
        additional_parameters: Vec<Atom>,
        external_count: usize,
        mut edges: Vec<usize>,
        esurfaces: Vec<ThresholdMultiplierEsurface>,
    ) -> Result<Self> {
        for (index, parameter) in model_parameters.iter().enumerate() {
            if model_parameters[..index].contains(parameter) {
                return Err(eyre!(
                    "duplicate model parameter `{parameter}` in threshold-multiplier layout"
                ));
            }
        }
        for (index, parameter) in additional_parameters.iter().enumerate() {
            if additional_parameters[..index].contains(parameter) {
                return Err(eyre!(
                    "duplicate additional parameter `{parameter}` in threshold-multiplier layout"
                ));
            }
            if model_parameters.contains(parameter) {
                return Err(eyre!(
                    "threshold-multiplier parameter `{parameter}` is declared as both a model and additional parameter"
                ));
            }
        }

        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!("duplicate graph edge in threshold-multiplier layout"));
        }
        let edge_set = edges.iter().copied().collect::<BTreeSet<_>>();

        let mut canonical_esurfaces = Vec::with_capacity(esurfaces.len());
        for mut esurface in esurfaces {
            if esurface.edges.is_empty() {
                return Err(eyre!(
                    "an E-surface multiplier input must contain at least one edge"
                ));
            }
            esurface.edges.sort_unstable();
            if esurface.edges.windows(2).any(|pair| pair[0] == pair[1]) {
                return Err(eyre!(
                    "E-surface multiplier input `{:?}` contains a duplicate energy edge",
                    esurface.edges,
                ));
            }
            if let Some(edge) = esurface
                .edges
                .iter()
                .chain(esurface.external_shift.iter().map(|(edge, _)| edge))
                .find(|edge| !edge_set.contains(edge))
            {
                return Err(eyre!(
                    "E-surface multiplier input `{:?}` refers to unknown graph edge {edge}",
                    esurface.edges,
                ));
            }

            let mut shift = BTreeMap::<usize, i64>::new();
            for (edge, coefficient) in esurface.external_shift {
                *shift.entry(edge).or_default() += coefficient;
            }
            esurface.external_shift = shift
                .into_iter()
                .filter(|(_, coefficient)| *coefficient != 0)
                .collect();
            canonical_esurfaces.push(esurface);
        }
        canonical_esurfaces.sort_unstable();
        let mut unique_esurfaces =
            Vec::<ThresholdMultiplierEsurface>::with_capacity(canonical_esurfaces.len());
        for esurface in canonical_esurfaces {
            if unique_esurfaces.last() == Some(&esurface) {
                continue;
            }
            unique_esurfaces.push(esurface);
        }

        let mut layout = Self {
            model_parameters,
            additional_parameters,
            external_count,
            edges,
            esurfaces: unique_esurfaces,
            inputs: Vec::new(),
            parameters: Vec::new(),
        };
        layout.populate_inputs();
        Ok(layout)
    }

    /// Build a layout from the graph's generation-time parameter ordering and a cut-owned
    /// E-surface catalog. E-surface IDs are deduplicated, while equations with equal energy edges
    /// but different external shifts remain distinct so an `eta(...)` use can report ambiguity.
    pub fn from_graph_esurfaces(
        graph: &Graph,
        esurface_ids: impl IntoIterator<Item = EsurfaceID>,
    ) -> Result<Self> {
        let mut esurface_ids = esurface_ids.into_iter().collect::<Vec<_>>();
        esurface_ids.sort_unstable();
        esurface_ids.dedup();
        let esurfaces = esurface_ids
            .into_iter()
            .map(|esurface_id| {
                let esurface = graph
                    .surface_cache
                    .esurface_cache
                    .raw
                    .get(esurface_id.0)
                    .ok_or_else(|| {
                        eyre!(
                            "threshold-multiplier catalog refers to missing E-surface {} in graph '{}'",
                            esurface_id.0,
                            graph.name,
                        )
                    })?;
                Ok(ThresholdMultiplierEsurface {
                    edges: esurface.energies.iter().map(|edge| edge.0).collect(),
                    external_shift: esurface
                        .external_shift
                        .iter()
                        .map(|(edge, coefficient)| (edge.0, *coefficient))
                        .collect(),
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Self::new(
            graph.param_builder.pairs.model_parameters.params.clone(),
            graph.param_builder.pairs.additional_params.params.clone(),
            graph.loop_momentum_basis.ext_edges.len(),
            (0..graph.n_edges()).collect(),
            esurfaces,
        )
    }

    fn populate_inputs(&mut self) {
        for (index, parameter) in self.model_parameters.iter().enumerate() {
            self.inputs
                .push(ThresholdMultiplierInput::ModelParameter { index });
            self.parameters.push(parameter.clone());
        }

        for (index, parameter) in self.additional_parameters.iter().enumerate() {
            self.inputs
                .push(ThresholdMultiplierInput::AdditionalParameter { index });
            self.parameters.push(parameter.clone());
        }

        for position in 0..self.external_count {
            for component in 0..4 {
                self.inputs
                    .push(ThresholdMultiplierInput::ExternalMomentum {
                        position,
                        component,
                    });
                self.parameters.push(function!(
                    GS.external_mom,
                    position as i64,
                    AIND_SYMBOLS.cind.call_args([component as i64])
                ));
            }
        }

        for point in [
            ThresholdMultiplierPoint::Effective,
            ThresholdMultiplierPoint::Star,
        ] {
            for &edge in &self.edges {
                // Q3 has a Lorentz index, but its temporal component is identically zero and
                // therefore does not consume an evaluator input slot.
                for component in 1..4 {
                    self.inputs.push(ThresholdMultiplierInput::EdgeMomentum {
                        point,
                        edge,
                        component,
                    });
                    self.parameters.push(function!(
                        GS.emr_vec,
                        point.atom(),
                        edge as i64,
                        AIND_SYMBOLS.cind.call_args([component as i64])
                    ));
                }
            }
        }

        for point in [
            ThresholdMultiplierPoint::Effective,
            ThresholdMultiplierPoint::Star,
        ] {
            for esurface in 0..self.esurfaces.len() {
                self.inputs
                    .push(ThresholdMultiplierInput::Esurface { point, esurface });
                self.parameters.push(eta_parameter_atom(point, esurface));
            }
        }
    }

    pub fn inputs(&self) -> &[ThresholdMultiplierInput] {
        &self.inputs
    }

    pub fn parameters(&self) -> &[Atom] {
        &self.parameters
    }

    pub fn model_parameters(&self) -> &[Atom] {
        &self.model_parameters
    }

    pub fn additional_parameters(&self) -> &[Atom] {
        &self.additional_parameters
    }

    pub fn external_count(&self) -> usize {
        self.external_count
    }

    pub fn edges(&self) -> &[usize] {
        &self.edges
    }

    pub fn esurfaces(&self) -> &[ThresholdMultiplierEsurface] {
        &self.esurfaces
    }

    pub fn len(&self) -> usize {
        self.inputs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inputs.is_empty()
    }

    pub fn parse_expression(&self, source: &str) -> Result<ThresholdMultiplierExpression> {
        // The public Q3 symbol has a global normalizer that turns concrete spatial Q3
        // components into Q components. Parse multiplier Q3 calls through a private tensor
        // symbol so that the multiplier ABI can reject genuine Q while preserving Q3 aliases.
        let private_q3 = private_q3_symbol();
        let rewritten = rewrite_q3_function_name(source, private_q3.get_name());
        let parsed = try_parse!(rewritten.as_str()).map_err(|error| {
            eyre!("failed to parse threshold-multiplier expression `{source}`: {error}")
        })?;
        let normalized = self.normalize_before_contraction(&parsed, private_q3)?;
        let scalar = scalarize(normalized.as_view()).with_context(|| {
            format!("threshold-multiplier expression `{source}` is not a scalar tensor expression")
        })?;
        let scalar = self.normalize_after_contraction(&scalar, private_q3)?;
        Ok(ThresholdMultiplierExpression { scalar })
    }

    fn normalize_before_contraction(&self, atom: &Atom, private_q3: Symbol) -> Result<Atom> {
        let ascii_eta = symbol!("eta");
        let mut error = None;
        let normalized = atom.replace_map(|view, _, output| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = view else {
                return;
            };
            let function_symbol = function.get_symbol();
            if function_symbol == GS.emr_mom || function_symbol == GS.loop_mom {
                error = Some(eyre!(
                    "threshold multipliers do not expose `{}`; use `Q3` for graph-edge spatial momenta",
                    function_symbol.get_name()
                ));
                return;
            }

            if function_symbol == private_q3 || function_symbol == GS.emr_vec {
                match self.normalize_q3_call(function.iter().collect(), private_q3) {
                    Ok(replacement) => **output = replacement,
                    Err(call_error) => error = Some(call_error),
                }
            } else if function_symbol == ascii_eta || function_symbol == GS.eta {
                match self.normalize_eta_call(function.iter().collect()) {
                    Ok(replacement) => **output = replacement,
                    Err(call_error) => error = Some(call_error),
                }
            } else if function_symbol == GS.external_mom {
                let arguments = function.iter().collect::<Vec<_>>();
                if arguments.len() != 2 {
                    error = Some(eyre!(
                        "P expects `(external_position, lorentz_index)`, got {} arguments",
                        arguments.len()
                    ));
                } else {
                    match parse_index(arguments[0], "external position") {
                        Ok(position) if position < self.external_count => {}
                        Ok(position) => {
                            error = Some(eyre!(
                                "P external position {position} is out of range for {} external momenta",
                                self.external_count
                            ));
                        }
                        Err(call_error) => error = Some(call_error),
                    }
                }
            }
        });
        error.map_or(Ok(normalized), Err)
    }

    fn normalize_q3_call(&self, arguments: Vec<AtomView<'_>>, symbol: Symbol) -> Result<Atom> {
        let (point, edge, index) = match arguments.as_slice() {
            [edge, index] => (ThresholdMultiplierPoint::Effective, *edge, *index),
            [point, edge, index] => (ThresholdMultiplierPoint::parse(*point)?, *edge, *index),
            _ => {
                return Err(eyre!(
                    "Q3 expects `(edge, lorentz_index)` or `(effective|star, edge, lorentz_index)`, got {} arguments",
                    arguments.len()
                ));
            }
        };
        let edge = parse_index(edge, "Q3 graph edge")?;
        if self.edges.binary_search(&edge).is_err() {
            return Err(eyre!("Q3 refers to unknown graph edge {edge}"));
        }
        Ok(symbol.call_args([
            point.atom().as_view(),
            Atom::num(edge as i64).as_view(),
            index,
        ]))
    }

    fn normalize_eta_call(&self, arguments: Vec<AtomView<'_>>) -> Result<Atom> {
        if arguments.is_empty() {
            return Err(eyre!("eta expects at least one graph edge"));
        }

        let (point, edge_arguments) = match arguments[0] {
            AtomView::Var(_) => (
                ThresholdMultiplierPoint::parse(arguments[0])?,
                &arguments[1..],
            ),
            _ => (ThresholdMultiplierPoint::Effective, arguments.as_slice()),
        };
        if edge_arguments.is_empty() {
            return Err(eyre!("eta expects at least one graph edge"));
        }
        let mut edges = edge_arguments
            .iter()
            .map(|edge| parse_index(*edge, "eta graph edge"))
            .collect::<Result<Vec<_>>>()?;
        edges.sort_unstable();
        if edges.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(eyre!(
                "eta edge set `{edges:?}` contains a duplicate graph edge"
            ));
        }
        let matching = self
            .esurfaces
            .iter()
            .enumerate()
            .filter(|(_, esurface)| esurface.edges == edges)
            .collect::<Vec<_>>();
        if matching.is_empty() {
            return Err(eyre!(
                "eta edge set `{edges:?}` is not an E-surface in this cut"
            ));
        }
        if matching.len() > 1 {
            return Err(eyre!(
                "eta edge set `{edges:?}` is ambiguous in this cut; matching external shifts are {:?}",
                matching
                    .iter()
                    .map(|(_, esurface)| &esurface.external_shift)
                    .collect::<Vec<_>>(),
            ));
        }
        Ok(eta_parameter_atom(point, matching[0].0))
    }

    fn normalize_after_contraction(&self, atom: &Atom, private_q3: Symbol) -> Result<Atom> {
        let mut error = None;
        let normalized = atom.replace_map(|view, _, output| {
            if error.is_some() {
                return;
            }
            let AtomView::Fun(function) = view else {
                return;
            };
            let function_symbol = function.get_symbol();
            if function_symbol == GS.emr_mom || function_symbol == GS.loop_mom {
                error = Some(eyre!(
                    "threshold multipliers do not expose `{}`; use `Q3` for graph-edge spatial momenta",
                    function_symbol.get_name()
                ));
                return;
            }

            if function_symbol == private_q3 || function_symbol == GS.emr_vec {
                let arguments = function.iter().collect::<Vec<_>>();
                if arguments.len() != 3 {
                    error = Some(eyre!(
                        "internal Q3 contraction produced {} arguments instead of 3",
                        arguments.len()
                    ));
                    return;
                }
                let result = (|| {
                    let point = ThresholdMultiplierPoint::parse(arguments[0])?;
                    let edge = parse_index(arguments[1], "Q3 graph edge")?;
                    let component = parse_component(arguments[2], "Q3")?;
                    if component == 0 {
                        Ok(Atom::Zero)
                    } else {
                        Ok(function!(
                            GS.emr_vec,
                            point.atom(),
                            edge as i64,
                            AIND_SYMBOLS.cind.call_args([component as i64])
                        ))
                    }
                })();
                match result {
                    Ok(replacement) => **output = replacement,
                    Err(call_error) => error = Some(call_error),
                }
            } else if function_symbol == GS.external_mom {
                let arguments = function.iter().collect::<Vec<_>>();
                if arguments.len() != 2 {
                    error = Some(eyre!(
                        "internal P contraction produced {} arguments instead of 2",
                        arguments.len()
                    ));
                    return;
                }
                let result = (|| {
                    let position = parse_index(arguments[0], "external position")?;
                    let component = parse_component(arguments[1], "P")?;
                    Ok(function!(
                        GS.external_mom,
                        position as i64,
                        AIND_SYMBOLS.cind.call_args([component as i64])
                    ))
                })();
                match result {
                    Ok(replacement) => **output = replacement,
                    Err(call_error) => error = Some(call_error),
                }
            }
        });
        error.map_or(Ok(normalized), Err)
    }

    pub fn build_evaluator(
        &self,
        expression: &ThresholdMultiplierExpression,
        settings: &EvaluatorSettings,
    ) -> Result<ThresholdMultiplierEvaluator> {
        let evaluator = GenericEvaluator::new_from_raw_params(
            [expression.scalar.clone()],
            &self.parameters,
            &FunctionMap::new(),
            Vec::new(),
            settings.optimization_settings(),
            None,
            settings,
        )?
        .into_eager_only();
        Ok(ThresholdMultiplierEvaluator {
            evaluator,
            input_count: self.len(),
            expression: expression.scalar.clone(),
        })
    }
}

/// A parsed, validated and scalarized multiplier expression.
#[derive(Clone, Debug)]
pub struct ThresholdMultiplierExpression {
    scalar: Atom,
}

impl ThresholdMultiplierExpression {
    pub fn scalar(&self) -> &Atom {
        &self.scalar
    }
}

/// The serializable eager-only evaluator for one scalar multiplier.
#[derive(Clone, Debug, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ThresholdMultiplierEvaluator {
    evaluator: GenericEvaluator,
    input_count: usize,
    expression: Atom,
}

impl ThresholdMultiplierEvaluator {
    pub fn evaluate<T: FloatLike>(
        &mut self,
        values: &ThresholdMultiplierInputValues<T>,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
    ) -> Result<F<T>> {
        if values.len() != self.input_count {
            return Err(eyre!(
                "threshold-multiplier evaluator expected {} inputs, got {}",
                self.input_count,
                values.len()
            ));
        }
        let value = evaluate_evaluator_single(
            &mut self.evaluator,
            values.as_slice(),
            evaluation_metadata,
            record_primary_timing,
        );
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err(eyre!(
                "threshold multiplier evaluated to non-finite value {value}"
            ));
        }
        if !value.im.is_zero() {
            return Err(eyre!(
                "threshold multiplier evaluated to non-real value {value}"
            ));
        }
        Ok(value.re)
    }

    pub fn generic_evaluator(&self) -> &GenericEvaluator {
        &self.evaluator
    }

    pub fn generic_evaluator_mut(&mut self) -> &mut GenericEvaluator {
        &mut self.evaluator
    }

    pub fn expression(&self) -> &Atom {
        &self.expression
    }

    pub fn input_count(&self) -> usize {
        self.input_count
    }
}

#[derive(Clone, Copy, Debug, Encode, Decode, Eq, PartialEq)]
pub struct ThresholdMultiplierEvaluatorId(pub usize);

/// The multiplier evaluator, if any, associated with one resolved threshold variant.
#[derive(Clone, Copy, Debug, Encode, Decode, Eq, PartialEq)]
pub struct ThresholdMultiplierVariantReference {
    pub variant_id: ThresholdCountertermVariantId,
    pub evaluator_id: Option<ThresholdMultiplierEvaluatorId>,
}

/// One cut-local multiplier registry. It is stored behind `Option` by the LU counterterm so the
/// legacy identity path retains no layout, evaluator, or workspace allocation.
#[derive(Clone, Debug, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ThresholdMultiplierEvaluatorCollection {
    layout: ThresholdMultiplierLayout,
    evaluators: Vec<ThresholdMultiplierEvaluator>,
    left_variants: Vec<ThresholdMultiplierVariantReference>,
    right_variants: Vec<ThresholdMultiplierVariantReference>,
}

impl ThresholdMultiplierEvaluatorCollection {
    /// Builds and deduplicates scalar evaluators. `None` is returned when every variant has an
    /// identity multiplier, allowing the caller to preserve the no-directive fast path.
    pub fn build(
        layout: ThresholdMultiplierLayout,
        left: Vec<(
            ThresholdCountertermVariantId,
            Option<ThresholdMultiplierExpression>,
        )>,
        right: Vec<(
            ThresholdCountertermVariantId,
            Option<ThresholdMultiplierExpression>,
        )>,
        settings: &EvaluatorSettings,
    ) -> Result<Option<Self>> {
        let mut canonical_expressions = Vec::<Atom>::new();
        let mut evaluators = Vec::new();
        let mut intern = |expression: Option<ThresholdMultiplierExpression>| -> Result<_> {
            let Some(expression) = expression else {
                return Ok(None);
            };
            let evaluator_id = if let Some(index) = canonical_expressions
                .iter()
                .position(|canonical| canonical == expression.scalar())
            {
                ThresholdMultiplierEvaluatorId(index)
            } else {
                let evaluator_id = ThresholdMultiplierEvaluatorId(evaluators.len());
                evaluators.push(layout.build_evaluator(&expression, settings)?);
                canonical_expressions.push(expression.scalar);
                evaluator_id
            };
            Ok(Some(evaluator_id))
        };

        let mut build_references = |variants: Vec<(
            ThresholdCountertermVariantId,
            Option<ThresholdMultiplierExpression>,
        )>| {
            variants
                .into_iter()
                .map(|(variant_id, expression)| {
                    Ok(ThresholdMultiplierVariantReference {
                        variant_id,
                        evaluator_id: intern(expression)?,
                    })
                })
                .collect::<Result<Vec<_>>>()
        };
        let left_variants = build_references(left)?;
        let right_variants = build_references(right)?;
        if evaluators.is_empty() {
            return Ok(None);
        }

        let collection = Self {
            layout,
            evaluators,
            left_variants,
            right_variants,
        };
        collection.validate(
            collection.left_variants.len(),
            collection.right_variants.len(),
        )?;
        Ok(Some(collection))
    }

    pub fn validate(&self, expected_left: usize, expected_right: usize) -> Result<()> {
        if self.left_variants.len() != expected_left || self.right_variants.len() != expected_right
        {
            return Err(eyre!(
                "threshold-multiplier registry dimensions are {}x{}, expected {expected_left}x{expected_right}",
                self.left_variants.len(),
                self.right_variants.len(),
            ));
        }
        if self.evaluators.is_empty() {
            return Err(eyre!(
                "an empty threshold-multiplier registry must be represented by `None`"
            ));
        }

        let mut variant_ids = BTreeSet::new();
        let mut referenced_evaluators = vec![false; self.evaluators.len()];
        for reference in self.left_variants.iter().chain(&self.right_variants) {
            if !variant_ids.insert(reference.variant_id) {
                return Err(eyre!(
                    "threshold-counterterm variant {} occurs more than once in one multiplier registry",
                    reference.variant_id.0,
                ));
            }
            if let Some(evaluator_id) = reference.evaluator_id {
                let referenced = referenced_evaluators.get_mut(evaluator_id.0).ok_or_else(|| {
                    eyre!(
                        "threshold-counterterm variant {} references missing multiplier evaluator {} (registry has {})",
                        reference.variant_id.0,
                        evaluator_id.0,
                        self.evaluators.len(),
                    )
                })?;
                *referenced = true;
            }
        }
        if let Some(orphan) = referenced_evaluators
            .iter()
            .position(|referenced| !referenced)
        {
            return Err(eyre!(
                "threshold-multiplier evaluator {orphan} is not referenced by any variant"
            ));
        }

        for (index, evaluator) in self.evaluators.iter().enumerate() {
            if evaluator.input_count() != self.layout.len() {
                return Err(eyre!(
                    "threshold-multiplier evaluator {index} expects {} inputs but its layout has {}",
                    evaluator.input_count(),
                    self.layout.len(),
                ));
            }
            if evaluator.generic_evaluator().exprs_len != 1
                || evaluator.generic_evaluator().dual_shape.is_some()
                || evaluator.generic_evaluator().backend_policy != EvaluatorBackendPolicy::EagerOnly
            {
                return Err(eyre!(
                    "threshold-multiplier evaluator {index} must be scalar, non-dual, and eager-only"
                ));
            }
        }
        Ok(())
    }

    pub fn layout(&self) -> &ThresholdMultiplierLayout {
        &self.layout
    }

    pub fn evaluators(&self) -> &[ThresholdMultiplierEvaluator] {
        &self.evaluators
    }

    pub fn evaluators_mut(&mut self) -> &mut [ThresholdMultiplierEvaluator] {
        &mut self.evaluators
    }

    pub fn left_variants(&self) -> &[ThresholdMultiplierVariantReference] {
        &self.left_variants
    }

    pub fn right_variants(&self) -> &[ThresholdMultiplierVariantReference] {
        &self.right_variants
    }

    pub fn evaluator_id_for_variant(
        &self,
        variant_id: ThresholdCountertermVariantId,
    ) -> Result<Option<ThresholdMultiplierEvaluatorId>> {
        self.left_variants
            .iter()
            .chain(&self.right_variants)
            .find(|reference| reference.variant_id == variant_id)
            .map(|reference| reference.evaluator_id)
            .ok_or_else(|| {
                eyre!(
                    "threshold-counterterm variant {} is absent from the multiplier registry",
                    variant_id.0,
                )
            })
    }

    /// Project the cut-local evaluator collection into neutral static-metadata registrations.
    /// The graph registry assigns final evaluator IDs after concatenating these in cut order.
    pub(crate) fn metadata_registrations(
        &self,
        cut_group_id: Option<usize>,
    ) -> Vec<ThresholdCountertermEvaluatorRegistration> {
        self.evaluators
            .iter()
            .enumerate()
            .map(|(collection_evaluator_id, evaluator)| {
                let variant_ids = self
                    .left_variants
                    .iter()
                    .chain(&self.right_variants)
                    .filter_map(|reference| {
                        (reference.evaluator_id
                            == Some(ThresholdMultiplierEvaluatorId(collection_evaluator_id)))
                        .then_some(reference.variant_id.0)
                    })
                    .collect();
                ThresholdCountertermEvaluatorRegistration {
                    cut_group_id,
                    collection_evaluator_id,
                    expression: evaluator.expression().to_string(),
                    variant_ids,
                }
            })
            .collect()
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for evaluator in &mut self.evaluators {
            f(evaluator.generic_evaluator_mut())?;
        }
        Ok(())
    }
}

/// Reusable storage for the values matching a [`ThresholdMultiplierLayout`].
#[derive(Clone, Debug)]
pub struct ThresholdMultiplierInputValues<T: FloatLike> {
    values: Vec<Complex<F<T>>>,
}

impl<T: FloatLike> ThresholdMultiplierInputValues<T> {
    pub fn new(layout: &ThresholdMultiplierLayout, zero: F<T>) -> Self {
        Self {
            values: vec![Complex::new_re(zero); layout.len()],
        }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn as_slice(&self) -> &[Complex<F<T>>] {
        &self.values
    }

    pub fn as_mut_slice(&mut self) -> &mut [Complex<F<T>>] {
        &mut self.values
    }

    pub fn set(&mut self, index: usize, value: Complex<F<T>>) -> Result<()> {
        let input_count = self.values.len();
        let slot = self.values.get_mut(index).ok_or_else(|| {
            eyre!(
                "threshold-multiplier input index {index} is out of range for {} inputs",
                input_count
            )
        })?;
        *slot = value;
        Ok(())
    }

    pub fn set_real(&mut self, index: usize, value: F<T>) -> Result<()> {
        self.set(index, Complex::new_re(value))
    }
}

struct ThresholdMultiplierKinematicPoint<T: FloatLike> {
    loop_momenta: LoopMomenta<F<T>>,
    external_momenta: ExternalFourMomenta<F<T>>,
    edge_momenta: Vec<ThreeMomentum<F<T>>>,
    esurface_values: Vec<F<T>>,
}

impl<T: FloatLike> ThresholdMultiplierKinematicPoint<T> {
    fn matches(&self, sample: &MomentumSample<T>) -> bool {
        &self.loop_momenta == sample.loop_moms() && &self.external_momenta == sample.external_moms()
    }
}

/// Reusable cut-local storage for binding the fixed multiplier ABI to expanded components.
///
/// Both samples must already be expressed in `generation_lmb`. For a single-sided component the
/// caller constructs `star` by applying that side's root; the absent side has no separate root, so
/// its spectator parent-LMB coordinates use the identity transformation while its original
/// integrand factor remains present. For an iterated component, `star` contains the selected pair's
/// merged roots. Each distinct complete sample is cached as one indivisible kinematic view: edge
/// momenta and E-surfaces are always recomputed after all participating roots have been merged and
/// are never assembled by splicing per-edge values from separate roots.
pub struct ThresholdMultiplierInputWorkspace<T: FloatLike> {
    layout_inputs: Vec<ThresholdMultiplierInput>,
    layout_esurfaces: Vec<ThresholdMultiplierEsurface>,
    values: ThresholdMultiplierInputValues<T>,
    model_values: Option<Vec<Complex<F<f64>>>>,
    additional_values: Option<Vec<F<T>>>,
    external_momenta: Option<ExternalFourMomenta<F<T>>>,
    generation_lmb: Option<LoopMomentumBasis>,
    edge_masses: Option<EdgeVec<F<T>>>,
    kinematic_points: Vec<ThresholdMultiplierKinematicPoint<T>>,
    external_spatial_scratch: ExternalThreeMomenta<F<T>>,
}

impl<T: FloatLike> ThresholdMultiplierInputWorkspace<T> {
    pub fn new(layout: &ThresholdMultiplierLayout, zero: F<T>) -> Self {
        let mut external_spatial_scratch = ExternalThreeMomenta::new();
        external_spatial_scratch.raw.reserve(layout.external_count);
        Self {
            layout_inputs: layout.inputs.clone(),
            layout_esurfaces: layout.esurfaces.clone(),
            values: ThresholdMultiplierInputValues::new(layout, zero),
            model_values: None,
            additional_values: None,
            external_momenta: None,
            generation_lmb: None,
            edge_masses: None,
            // A one-sided context needs base/root and an iterated context needs
            // base/left/right/merged. Additional variants grow this cut-local cache as needed.
            kinematic_points: Vec::with_capacity(4),
            external_spatial_scratch,
        }
    }

    fn kinematic_point(
        &mut self,
        layout: &ThresholdMultiplierLayout,
        generation_lmb: &LoopMomentumBasis,
        edge_masses: &EdgeVec<F<T>>,
        sample: &MomentumSample<T>,
    ) -> Result<usize> {
        if let Some(index) = self
            .kinematic_points
            .iter()
            .position(|point| point.matches(sample))
        {
            return Ok(index);
        }

        let mut edge_momenta = Vec::with_capacity(layout.edges.len());
        let mut esurface_values = Vec::with_capacity(layout.esurfaces.len());
        fill_kinematic_point(
            layout,
            generation_lmb,
            edge_masses,
            sample,
            &mut self.external_spatial_scratch,
            &mut edge_momenta,
            &mut esurface_values,
        )?;
        self.kinematic_points
            .push(ThresholdMultiplierKinematicPoint {
                loop_momenta: sample.loop_moms().clone(),
                external_momenta: sample.external_moms().clone(),
                edge_momenta,
                esurface_values,
            });
        Ok(self.kinematic_points.len() - 1)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn bind<'a>(
        &'a mut self,
        layout: &ThresholdMultiplierLayout,
        generation_lmb: &LoopMomentumBasis,
        edge_masses: &EdgeVec<F<T>>,
        model_values: &[Complex<F<f64>>],
        additional_values: &[F<T>],
        effective: &MomentumSample<T>,
        star: &MomentumSample<T>,
    ) -> Result<&'a ThresholdMultiplierInputValues<T>> {
        if self.layout_inputs != layout.inputs || self.layout_esurfaces != layout.esurfaces {
            return Err(eyre!(
                "threshold-multiplier workspace was constructed for a different input layout"
            ));
        }
        if model_values.len() != layout.model_parameters.len() {
            return Err(eyre!(
                "threshold-multiplier layout expects {} model values, got {}",
                layout.model_parameters.len(),
                model_values.len(),
            ));
        }
        if additional_values.len() != layout.additional_parameters.len() {
            return Err(eyre!(
                "threshold-multiplier layout expects {} additional parameter values, got {}",
                layout.additional_parameters.len(),
                additional_values.len(),
            ));
        }

        if self.generation_lmb.as_ref() != Some(generation_lmb)
            || self.edge_masses.as_ref() != Some(edge_masses)
        {
            self.generation_lmb = Some(generation_lmb.clone());
            self.edge_masses = Some(edge_masses.clone());
            self.kinematic_points.clear();
        }

        if self.model_values.as_deref() != Some(model_values) {
            for (input_index, input) in layout.inputs.iter().enumerate() {
                if let ThresholdMultiplierInput::ModelParameter { index } = *input {
                    self.values.values[input_index] = Complex {
                        re: F::from_ff64(model_values[index].re),
                        im: F::from_ff64(model_values[index].im),
                    };
                }
            }
            self.model_values = Some(model_values.to_vec());
        }
        if self.additional_values.as_deref() != Some(additional_values) {
            for (input_index, input) in layout.inputs.iter().enumerate() {
                if let ThresholdMultiplierInput::AdditionalParameter { index } = *input {
                    self.values.values[input_index] =
                        Complex::new_re(additional_values[index].clone());
                }
            }
            self.additional_values = Some(additional_values.to_vec());
        }
        if self.external_momenta.as_ref() != Some(effective.external_moms()) {
            for (input_index, input) in layout.inputs.iter().enumerate() {
                if let ThresholdMultiplierInput::ExternalMomentum {
                    position,
                    component,
                } = *input
                {
                    self.values.values[input_index] = Complex::new_re(four_momentum_component(
                        &effective.external_moms()
                            [crate::momentum::sample::ExternalIndex(position)],
                        component,
                    ));
                }
            }
            self.external_momenta = Some(effective.external_moms().clone());
        }

        let effective_point =
            self.kinematic_point(layout, generation_lmb, edge_masses, effective)?;
        let star_point = if std::ptr::eq(effective, star) {
            effective_point
        } else {
            self.kinematic_point(layout, generation_lmb, edge_masses, star)?
        };

        let points = &self.kinematic_points;
        for (input_index, input) in layout.inputs.iter().enumerate() {
            let value = match *input {
                ThresholdMultiplierInput::ModelParameter { .. }
                | ThresholdMultiplierInput::AdditionalParameter { .. }
                | ThresholdMultiplierInput::ExternalMomentum { .. } => continue,
                ThresholdMultiplierInput::EdgeMomentum {
                    point,
                    edge,
                    component,
                } => {
                    let edge_position = layout.edges.binary_search(&edge).map_err(|_| {
                        eyre!("threshold-multiplier input refers to missing graph edge {edge}")
                    })?;
                    let point = match point {
                        ThresholdMultiplierPoint::Effective => effective_point,
                        ThresholdMultiplierPoint::Star => star_point,
                    };
                    let momentum = &points[point].edge_momenta[edge_position];
                    Complex::new_re(three_momentum_component(momentum, component)?)
                }
                ThresholdMultiplierInput::Esurface { point, esurface } => {
                    let point = match point {
                        ThresholdMultiplierPoint::Effective => effective_point,
                        ThresholdMultiplierPoint::Star => star_point,
                    };
                    Complex::new_re(points[point].esurface_values[esurface].clone())
                }
            };
            self.values.values[input_index] = value;
        }
        Ok(&self.values)
    }

    pub fn values(&self) -> &ThresholdMultiplierInputValues<T> {
        &self.values
    }
}

fn fill_kinematic_point<T: FloatLike>(
    layout: &ThresholdMultiplierLayout,
    generation_lmb: &LoopMomentumBasis,
    edge_masses: &EdgeVec<F<T>>,
    sample: &MomentumSample<T>,
    external_spatial: &mut ExternalThreeMomenta<F<T>>,
    edge_momenta: &mut Vec<ThreeMomentum<F<T>>>,
    esurface_values: &mut Vec<F<T>>,
) -> Result<()> {
    if sample.loop_moms().0.len() != generation_lmb.loop_edges.len() {
        return Err(eyre!(
            "threshold-multiplier sample has {} loop momenta but the generation LMB has {}",
            sample.loop_moms().0.len(),
            generation_lmb.loop_edges.len(),
        ));
    }
    if sample.external_moms().len() != layout.external_count {
        return Err(eyre!(
            "threshold-multiplier sample has {} external momenta but the layout expects {}",
            sample.external_moms().len(),
            layout.external_count,
        ));
    }

    external_spatial.clear();
    external_spatial.extend(
        sample
            .external_moms()
            .iter()
            .map(|momentum| momentum.spatial.clone()),
    );
    edge_momenta.clear();
    for &edge in &layout.edges {
        let edge = EdgeIndex(edge);
        let signature = generation_lmb.edge_signatures.get(edge).ok_or_else(|| {
            eyre!(
                "threshold-multiplier layout graph edge {} is absent from the generation LMB",
                edge.0,
            )
        })?;
        edge_momenta.push(signature.compute_momentum(sample.loop_moms(), external_spatial));
    }

    esurface_values.clear();
    for equation in &layout.esurfaces {
        let mut value = sample.sample.zero();
        for &edge in &equation.edges {
            let edge_position = layout.edges.binary_search(&edge).map_err(|_| {
                eyre!("threshold-multiplier E-surface refers to missing graph edge {edge}")
            })?;
            let mass = edge_masses.get(EdgeIndex(edge)).ok_or_else(|| {
                eyre!("threshold-multiplier masses are missing graph edge {edge}")
            })?;
            let momentum = &edge_momenta[edge_position];
            value += (momentum.norm_squared() + mass * mass).sqrt();
        }
        for &(edge, coefficient) in &equation.external_shift {
            let signature = generation_lmb
                .edge_signatures
                .get(EdgeIndex(edge))
                .ok_or_else(|| {
                    eyre!(
                        "threshold-multiplier E-surface shift refers to graph edge {edge} absent from the generation LMB"
                    )
                })?;
            let temporal = signature
                .external
                .try_apply(&sample.external_moms().raw)
                .map(|momentum| momentum.temporal.value)
                .unwrap_or_else(|| value.zero());
            let coefficient = value.from_i64(coefficient);
            value += coefficient * temporal;
        }
        esurface_values.push(value);
    }
    Ok(())
}

fn four_momentum_component<T: FloatLike>(
    momentum: &crate::momentum::FourMomentum<F<T>>,
    component: usize,
) -> F<T> {
    match component {
        0 => momentum.temporal.value.clone(),
        1 => momentum.spatial.px.clone(),
        2 => momentum.spatial.py.clone(),
        3 => momentum.spatial.pz.clone(),
        _ => unreachable!("threshold-multiplier layout validates Lorentz components"),
    }
}

fn three_momentum_component<T: FloatLike>(
    momentum: &ThreeMomentum<F<T>>,
    component: usize,
) -> Result<F<T>> {
    match component {
        1 => Ok(momentum.px.clone()),
        2 => Ok(momentum.py.clone()),
        3 => Ok(momentum.pz.clone()),
        _ => Err(eyre!(
            "threshold-multiplier Q3 component {component} is not spatial"
        )),
    }
}

fn eta_parameter_atom(point: ThresholdMultiplierPoint, esurface: usize) -> Atom {
    private_eta_symbol().call_args([point.atom(), Atom::num(esurface as i64)])
}

fn private_eta_symbol() -> Symbol {
    symbol!("gammalooprs::threshold_multiplier_eta")
}

fn private_q3_symbol() -> Symbol {
    use spenso::network::tags::SPENSO_TAG;

    symbol!(
        "gammalooprs::threshold_multiplier_Q3",
        tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
    )
}

fn rewrite_q3_function_name(source: &str, private_name: &str) -> String {
    let mut rewritten = String::with_capacity(source.len());
    let mut rest = source;
    while let Some(position) = rest.find("Q3") {
        let (prefix, candidate) = rest.split_at(position);
        rewritten.push_str(prefix);
        let after_name = &candidate[2..];
        let previous_is_identifier = prefix.chars().next_back().is_some_and(|character| {
            character.is_alphanumeric() || character == '_' || character == ':'
        });
        let is_call = after_name.trim_start().starts_with('(');
        if previous_is_identifier || !is_call {
            rewritten.push_str("Q3");
        } else {
            rewritten.push_str(private_name);
        }
        rest = after_name;
    }
    rewritten.push_str(rest);
    rewritten
}

fn parse_index(atom: AtomView<'_>, description: &str) -> Result<usize> {
    let value = i64::try_from(atom)
        .map_err(|_| eyre!("expected an integer {description}, got `{atom}`"))?;
    usize::try_from(value).map_err(|_| eyre!("expected a non-negative {description}, got {value}"))
}

fn parse_component(atom: AtomView<'_>, tensor: &str) -> Result<usize> {
    let AtomView::Fun(component) = atom else {
        return Err(eyre!(
            "{tensor} contraction did not produce a concrete Lorentz component: `{atom}`"
        ));
    };
    // Textual expressions may resolve `cind` in the default namespace, whereas Spenso emits
    // its registered `cind` symbol. Accept the public spelling and canonicalize it afterwards.
    if component.get_symbol().get_name().rsplit("::").next() != Some("cind")
        || component.get_nargs() != 1
    {
        return Err(eyre!(
            "{tensor} contraction did not produce a concrete Lorentz component: `{atom}`"
        ));
    }
    let value = parse_index(component.iter().next().unwrap(), "Lorentz component")?;
    if value >= 4 {
        return Err(eyre!(
            "{tensor} Lorentz component {value} is out of range for four dimensions"
        ));
    }
    Ok(value)
}

fn scalarize(atom: AtomView<'_>) -> Result<Atom> {
    let mut network = atom.parse_into_net().map_err(|error| eyre!(error))?;
    network
        .execute::<SequentialRef, SmallestDegree, _, _, _>(
            TENSORLIB.read().unwrap().deref(),
            FUN_LIB.deref(),
        )
        .map_err(|error| eyre!(error))?;
    network
        .execute::<SequentialRef, SmallestDegree, _, _, _>(
            TENSORLIB.read().unwrap().deref(),
            FUN_LIB.deref(),
        )
        .map_err(|error| eyre!(error))?;
    network
        .result_scalar()
        .map(|result| match result {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(value) => value.into_owned(),
        })
        .map_err(|error| eyre!(error))
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::initialisation::test_initialise;
    use crate::{
        GammaLoopContextContainer, dot,
        graph::parse::from_dot::IntoGraph,
        integrands::process::evaluators::{ActiveF64Backend, EvaluatorBackendPolicy},
        momentum::{
            FourMomentum,
            sample::{BareMomentumSample, ExternalFourMomenta, LoopMomenta},
        },
        settings::{
            RuntimeSettings,
            global::{CompilationOptimizationLevel, FrozenCompilationMode},
        },
        utils::{ArbPrec, f128, load_generic_model},
    };
    use symbolica::state::State;

    fn initialized_layout(
        model_parameters: Vec<Atom>,
        external_count: usize,
        edges: Vec<usize>,
        esurfaces: Vec<ThresholdMultiplierEsurface>,
    ) -> ThresholdMultiplierLayout {
        test_initialise().unwrap();
        ThresholdMultiplierLayout::new(
            model_parameters,
            Vec::new(),
            external_count,
            edges,
            esurfaces,
        )
        .unwrap()
    }

    fn esurface(
        edges: Vec<usize>,
        external_shift: Vec<(usize, i64)>,
    ) -> ThresholdMultiplierEsurface {
        ThresholdMultiplierEsurface {
            edges,
            external_shift,
        }
    }

    fn assert_context_cache_reuse<T: FloatLike>(
        layout: &ThresholdMultiplierLayout,
        generation_lmb: &LoopMomentumBasis,
        edge_count: usize,
        base: &MomentumSample<T>,
        root: &MomentumSample<T>,
    ) {
        let zero = base.zero();
        let masses = (0..edge_count)
            .map(|_| zero.clone())
            .collect::<EdgeVec<_>>();
        let mut workspace = ThresholdMultiplierInputWorkspace::new(layout, zero);
        workspace
            .bind(layout, generation_lmb, &masses, &[], &[], base, root)
            .unwrap();
        assert_eq!(workspace.kinematic_points.len(), 2);
        workspace
            .bind(layout, generation_lmb, &masses, &[], &[], root, root)
            .unwrap();
        assert_eq!(workspace.kinematic_points.len(), 2);
    }

    #[test]
    fn layout_is_deterministic_and_detects_ambiguous_esurfaces() {
        let layout = initialized_layout(
            vec![parse!("UFO::MT")],
            1,
            vec![3, 2, 0, 1],
            vec![
                esurface(vec![2, 1], vec![(0, 1), (3, -1), (0, 2)]),
                esurface(vec![0], Vec::new()),
                esurface(vec![1, 2], vec![(3, -1), (0, 3)]),
            ],
        );
        assert_eq!(layout.edges, vec![0, 1, 2, 3]);
        assert_eq!(
            layout.esurfaces,
            vec![
                esurface(vec![0], Vec::new()),
                esurface(vec![1, 2], vec![(0, 3), (3, -1)]),
            ]
        );
        assert_eq!(layout.len(), 1 + 4 + 2 * 4 * 3 + 2 * 2);

        let ambiguous_layout = ThresholdMultiplierLayout::new(
            Vec::new(),
            Vec::new(),
            0,
            vec![0, 1],
            vec![
                esurface(vec![0, 1], vec![(0, 1)]),
                esurface(vec![1, 0], vec![(0, 2)]),
            ],
        )
        .unwrap();
        let error = ambiguous_layout.parse_expression("eta(0, 1)").unwrap_err();
        assert!(error.to_string().contains("ambiguous in this cut"));
    }

    #[test]
    fn q3_alias_preserves_temporal_zero_and_minkowski_sign() {
        let layout = initialized_layout(Vec::new(), 0, vec![0, 1], Vec::new());
        let temporal = layout
            .parse_expression("Q3(0, cind(0)) + Q3(0, cind(1))")
            .unwrap();
        assert!(!temporal.scalar().to_string().contains("cind(0)"));

        let expression = layout
            .parse_expression("Q3(0, spenso::mink(4, 7)) * Q3(1, spenso::mink(4, 7))")
            .unwrap();
        let mut evaluator = layout
            .build_evaluator(&expression, &EvaluatorSettings::default())
            .unwrap();
        assert_eq!(
            evaluator.generic_evaluator().backend_policy,
            EvaluatorBackendPolicy::EagerOnly
        );

        let mut values = ThresholdMultiplierInputValues::new(&layout, F(0.0));
        for (index, input) in layout.inputs().iter().enumerate() {
            let ThresholdMultiplierInput::EdgeMomentum { component, .. } = input else {
                continue;
            };
            values
                .set_real(index, F([0.0, 2.0, 3.0, 6.0][*component]))
                .unwrap();
        }
        let result = evaluator
            .evaluate(&values, &mut EvaluationMetaData::new_empty(), false)
            .unwrap();
        assert_eq!(result, F(-49.0));
    }

    #[test]
    fn eta_views_share_the_fixed_layout() {
        let layout = initialized_layout(
            Vec::new(),
            0,
            vec![0, 1],
            vec![esurface(vec![1, 0], Vec::new())],
        );
        let expression = layout
            .parse_expression("eta(1, 0) + eta(star, 0, 1)")
            .unwrap();
        let mut evaluator = layout
            .build_evaluator(&expression, &EvaluatorSettings::default())
            .unwrap();
        let mut values = ThresholdMultiplierInputValues::new(&layout, F(0.0));
        for (index, input) in layout.inputs().iter().enumerate() {
            if let ThresholdMultiplierInput::Esurface { point, .. } = input {
                values
                    .set_real(
                        index,
                        F(match point {
                            ThresholdMultiplierPoint::Effective => 3.0,
                            ThresholdMultiplierPoint::Star => 4.0,
                        }),
                    )
                    .unwrap();
            }
        }
        let result = evaluator
            .evaluate(&values, &mut EvaluationMetaData::new_empty(), false)
            .unwrap();
        assert_eq!(result, F(7.0));
    }

    #[test]
    fn model_parameters_and_external_momenta_are_bound_in_existing_order() {
        let layout = initialized_layout(vec![parse!("UFO::MT")], 1, Vec::new(), Vec::new());
        let expression = layout.parse_expression("UFO::MT + P(0, cind(0))").unwrap();
        let mut evaluator = layout
            .build_evaluator(&expression, &EvaluatorSettings::default())
            .unwrap();
        let mut values = ThresholdMultiplierInputValues::new(&layout, F(0.0));
        for (index, input) in layout.inputs().iter().enumerate() {
            let value = match input {
                ThresholdMultiplierInput::ModelParameter { index: 0 } => 5.0,
                ThresholdMultiplierInput::ExternalMomentum {
                    position: 0,
                    component: 0,
                } => 7.0,
                _ => 0.0,
            };
            values.set_real(index, F(value)).unwrap();
        }
        let result = evaluator
            .evaluate(&values, &mut EvaluationMetaData::new_empty(), false)
            .unwrap();
        assert_eq!(result, F(12.0));
    }

    #[test]
    fn additional_parameters_follow_model_parameters_and_work_in_all_precisions() {
        test_initialise().unwrap();
        let layout = ThresholdMultiplierLayout::new(
            vec![parse!("UFO::MT")],
            vec![parse!("alpha"), parse!("beta")],
            0,
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        assert!(matches!(
            layout.inputs()[0],
            ThresholdMultiplierInput::ModelParameter { index: 0 }
        ));
        assert!(matches!(
            layout.inputs()[1],
            ThresholdMultiplierInput::AdditionalParameter { index: 0 }
        ));
        assert!(matches!(
            layout.inputs()[2],
            ThresholdMultiplierInput::AdditionalParameter { index: 1 }
        ));
        assert_eq!(
            layout.additional_parameters(),
            &[parse!("alpha"), parse!("beta")]
        );

        let expression = layout
            .parse_expression("UFO::MT + 2 * alpha + 3 * beta")
            .unwrap();
        let mut evaluator = layout
            .build_evaluator(&expression, &EvaluatorSettings::default())
            .unwrap();
        let mut f64_values = ThresholdMultiplierInputValues::new(&layout, F(0.0));
        let mut f128_values = ThresholdMultiplierInputValues::new(&layout, F(f128::from_f64(0.0)));
        let mut arb_values =
            ThresholdMultiplierInputValues::new(&layout, F(ArbPrec::from_f64(0.0)));
        f64_values.set_real(0, F(2.0)).unwrap();
        f128_values.set_real(0, F(f128::from_f64(2.0))).unwrap();
        arb_values.set_real(0, F(ArbPrec::from_f64(2.0))).unwrap();
        let mut runtime = RuntimeSettings::default();
        runtime.general.additional_param_values = vec![3.0, 5.0];
        for (index, value) in runtime.additional_params::<f64>().into_iter().enumerate() {
            f64_values.set_real(index + 1, value).unwrap();
        }
        for (index, value) in runtime.additional_params::<f128>().into_iter().enumerate() {
            f128_values.set_real(index + 1, value).unwrap();
        }
        for (index, value) in runtime
            .additional_params::<ArbPrec>()
            .into_iter()
            .enumerate()
        {
            arb_values.set_real(index + 1, value).unwrap();
        }
        let mut metadata = EvaluationMetaData::new_empty();
        assert_eq!(
            evaluator
                .evaluate(&f64_values, &mut metadata, false)
                .unwrap(),
            F(23.0)
        );
        assert_eq!(
            evaluator
                .evaluate(&f128_values, &mut metadata, false)
                .unwrap()
                .0
                .into_f64(),
            23.0
        );
        assert_eq!(
            evaluator
                .evaluate(&arb_values, &mut metadata, false)
                .unwrap()
                .0
                .into_f64(),
            23.0
        );

        let duplicate = ThresholdMultiplierLayout::new(
            Vec::new(),
            vec![parse!("alpha"), parse!("alpha")],
            0,
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();
        assert!(
            duplicate
                .to_string()
                .contains("duplicate additional parameter")
        );
        let collision = ThresholdMultiplierLayout::new(
            vec![parse!("alpha")],
            vec![parse!("alpha")],
            0,
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();
        assert!(
            collision
                .to_string()
                .contains("both a model and additional")
        );
    }

    #[test]
    fn graph_declared_additional_parameters_enter_the_cut_layout_in_dot_order() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph multiplier_additional_parameters {
            graph [params="alpha;beta";]
            ext [style=invis]
            v1 [num=1]
            v2 [num=1]
            edge [num=1 mass=0]
            ext -> v1 [id=0]
            v1 -> v2 [id=1]
            v1 -> v2 [id=2]
            v2 -> ext [id=3]
        })
        .unwrap();
        let layout = ThresholdMultiplierLayout::from_graph_esurfaces(&graph, Vec::new()).unwrap();
        assert_eq!(
            layout.additional_parameters(),
            &[parse!("alpha"), parse!("beta")]
        );
        assert!(layout.parse_expression("alpha + 2 * beta").is_ok());
    }

    #[test]
    fn forbidden_momenta_and_invalid_frames_are_rejected() {
        let layout = initialized_layout(Vec::new(), 0, vec![0], Vec::new());
        for expression in ["Q(0, cind(1))", "K(0, cind(1))"] {
            let error = layout.parse_expression(expression).unwrap_err();
            assert!(error.to_string().contains("use `Q3`"));
        }
        let error = layout.parse_expression("Q3(now, 0, cind(1))").unwrap_err();
        assert!(error.to_string().contains("unknown multiplier frame"));
    }

    #[test]
    fn collection_deduplicates_expressions_and_keeps_identity_references() {
        let layout = initialized_layout(Vec::new(), 0, Vec::new(), Vec::new());
        let expression = layout.parse_expression("7").unwrap();
        let collection = ThresholdMultiplierEvaluatorCollection::build(
            layout,
            vec![
                (ThresholdCountertermVariantId(0), Some(expression.clone())),
                (ThresholdCountertermVariantId(1), None),
            ],
            vec![(ThresholdCountertermVariantId(2), Some(expression))],
            &EvaluatorSettings::default(),
        )
        .unwrap()
        .unwrap();

        assert_eq!(collection.evaluators().len(), 1);
        assert_eq!(
            collection.left_variants(),
            [
                ThresholdMultiplierVariantReference {
                    variant_id: ThresholdCountertermVariantId(0),
                    evaluator_id: Some(ThresholdMultiplierEvaluatorId(0)),
                },
                ThresholdMultiplierVariantReference {
                    variant_id: ThresholdCountertermVariantId(1),
                    evaluator_id: None,
                },
            ]
        );
        assert_eq!(
            collection.right_variants(),
            [ThresholdMultiplierVariantReference {
                variant_id: ThresholdCountertermVariantId(2),
                evaluator_id: Some(ThresholdMultiplierEvaluatorId(0)),
            }]
        );
        let registrations = collection.metadata_registrations(Some(4));
        assert_eq!(registrations.len(), 1);
        assert_eq!(registrations[0].cut_group_id, Some(4));
        assert_eq!(registrations[0].collection_evaluator_id, 0);
        assert_eq!(registrations[0].variant_ids, vec![0, 2]);
        assert_eq!(registrations[0].expression, "7");
    }

    #[test]
    fn eager_only_collection_roundtrips_and_evaluates_in_all_precisions() {
        let layout = initialized_layout(Vec::new(), 0, Vec::new(), Vec::new());
        let expression = layout.parse_expression("7").unwrap();
        let collection = ThresholdMultiplierEvaluatorCollection::build(
            layout,
            vec![(ThresholdCountertermVariantId(0), Some(expression))],
            Vec::new(),
            &EvaluatorSettings::default(),
        )
        .unwrap()
        .unwrap();

        let encoded = bincode::encode_to_vec(&collection, bincode::config::standard()).unwrap();
        let model = load_generic_model("sm");
        let mut state_bytes = Vec::new();
        State::export(&mut state_bytes).unwrap();
        let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
        let context = GammaLoopContextContainer {
            model: &model,
            state_map: &state_map,
        };
        let (mut decoded, _): (ThresholdMultiplierEvaluatorCollection, _) =
            bincode::decode_from_slice_with_context(&encoded, bincode::config::standard(), context)
                .unwrap();
        decoded.validate(1, 0).unwrap();

        let f64_values = ThresholdMultiplierInputValues::new(decoded.layout(), F(0.0));
        let f128_values =
            ThresholdMultiplierInputValues::new(decoded.layout(), F(f128::from_f64(0.0)));
        let arb_values =
            ThresholdMultiplierInputValues::new(decoded.layout(), F(ArbPrec::from_f64(0.0)));
        {
            let evaluator = &mut decoded.evaluators_mut()[0];
            let mut metadata = EvaluationMetaData::new_empty();
            assert_eq!(
                evaluator
                    .evaluate(&f64_values, &mut metadata, false)
                    .unwrap(),
                F(7.0)
            );
            assert_eq!(
                evaluator
                    .evaluate(&f128_values, &mut metadata, false)
                    .unwrap()
                    .0
                    .into_f64(),
                7.0
            );
            assert_eq!(
                evaluator
                    .evaluate(&arb_values, &mut metadata, false)
                    .unwrap()
                    .0
                    .into_f64(),
                7.0
            );

            let generic = evaluator.generic_evaluator_mut();
            generic
                .activate_symjit(CompilationOptimizationLevel::O2)
                .unwrap();
            assert_eq!(generic.active_f64_backend(), ActiveF64Backend::Eager);
            generic
                .activate_external_from_artifact(ActiveF64Backend::Cpp)
                .unwrap();
            assert_eq!(generic.active_f64_backend(), ActiveF64Backend::Eager);
            generic
                .activate_external_from_artifact(ActiveF64Backend::Assembly)
                .unwrap();
            assert_eq!(generic.active_f64_backend(), ActiveF64Backend::Eager);
            generic
                .compile_external(
                    "unused-threshold-multiplier.cpp",
                    "unused_threshold_multiplier",
                    "unused-threshold-multiplier.so",
                    &FrozenCompilationMode::Eager,
                )
                .unwrap();
            assert_eq!(generic.active_f64_backend(), ActiveF64Backend::Eager);
        }

        let expected_visits = decoded.evaluators().len();
        let mut visits = 0;
        decoded
            .for_each_generic_evaluator_mut(|generic| {
                visits += 1;
                assert_eq!(generic.backend_policy, EvaluatorBackendPolicy::EagerOnly);
                assert_eq!(generic.active_f64_backend(), ActiveF64Backend::Eager);
                Ok(())
            })
            .unwrap();
        assert_eq!(visits, expected_visits);
    }

    #[test]
    fn workspace_binds_global_effective_and_star_views_in_generation_lmb() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph threshold_multiplier_workspace {
            graph [params="graph_weight"]
            ext [style=invis]
            v1 [num=1]
            v2 [num=1]
            edge [num=1 mass=0]
            ext -> v1 [id=0]
            v1 -> v2 [id=1]
            v1 -> v2 [id=2]
            v2 -> ext [id=3]
        })
        .unwrap();
        let energy_edge = graph.loop_momentum_basis.loop_edges.raw[0];
        let shift_edge = graph.loop_momentum_basis.ext_edges.raw[0];
        let layout = ThresholdMultiplierLayout::new(
            vec![parse!("UFO::MT")],
            graph.param_builder.pairs.additional_params.params.clone(),
            graph.loop_momentum_basis.ext_edges.len(),
            (0..graph.n_edges()).collect(),
            vec![esurface(vec![energy_edge.0], vec![(shift_edge.0, 1)])],
        )
        .unwrap();
        let expression = layout
            .parse_expression(&format!(
                "UFO::MT + graph_weight + P(0, cind(0)) + Q3({}, cind(1)) + Q3(star, {}, cind(1)) + eta({}) + eta(star, {})",
                energy_edge.0, energy_edge.0, energy_edge.0, energy_edge.0,
            ))
            .unwrap();
        let mut evaluator = layout
            .build_evaluator(&expression, &EvaluatorSettings::default())
            .unwrap();

        let make_sample = |loop_scale: f64| MomentumSample {
            sample: BareMomentumSample {
                loop_moms: (0..graph.loop_momentum_basis.loop_edges.len())
                    .map(|index| ThreeMomentum::new(F(loop_scale + index as f64), F(0.25), F(-0.5)))
                    .collect::<LoopMomenta<_>>(),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: (0..graph.loop_momentum_basis.ext_edges.len())
                    .map(|position| {
                        FourMomentum::from_args(
                            F(5.0 + position as f64),
                            F(0.1 * position as f64),
                            F(0.0),
                            F(0.0),
                        )
                    })
                    .collect::<ExternalFourMomenta<_>>(),
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let effective = make_sample(1.0);
        let star = make_sample(2.0);
        let masses = (0..graph.n_edges())
            .map(|edge| F(if edge == energy_edge.0 { 3.0 } else { 0.0 }))
            .collect::<EdgeVec<_>>();
        let model_values = [Complex::new_re(F(2.0))];
        let mut runtime = RuntimeSettings::default();
        runtime.general.additional_param_values = vec![11.0];
        let additional_values = runtime.additional_params();
        let mut workspace = ThresholdMultiplierInputWorkspace::new(&layout, F(0.0));
        let error = workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &model_values,
                &[],
                &effective,
                &star,
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("expects 1 additional parameter values")
        );
        let values = workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &model_values,
                &additional_values,
                &effective,
                &star,
            )
            .unwrap();

        let expected = layout
            .inputs()
            .iter()
            .enumerate()
            .filter(|(_, input)| match input {
                ThresholdMultiplierInput::ModelParameter { index: 0 }
                | ThresholdMultiplierInput::AdditionalParameter { index: 0 }
                | ThresholdMultiplierInput::ExternalMomentum {
                    position: 0,
                    component: 0,
                }
                | ThresholdMultiplierInput::Esurface { esurface: 0, .. } => true,
                ThresholdMultiplierInput::EdgeMomentum {
                    edge, component: 1, ..
                } => *edge == energy_edge.0,
                _ => false,
            })
            .fold(F(0.0), |sum, (index, _)| sum + values.as_slice()[index].re);
        let result = evaluator
            .evaluate(values, &mut EvaluationMetaData::new_empty(), false)
            .unwrap();
        assert!((result.0 - expected.0).abs() < 1.0e-12);

        let eta_values = layout
            .inputs()
            .iter()
            .enumerate()
            .filter_map(|(index, input)| {
                matches!(input, ThresholdMultiplierInput::Esurface { .. })
                    .then_some(values.as_slice()[index].re.0)
            })
            .collect::<Vec<_>>();
        assert_eq!(eta_values.len(), 2);
        assert_ne!(eta_values[0], eta_values[1]);

        let changed_model_values = [Complex::new_re(F(13.0))];
        let changed_additional_values = [F(17.0)];
        let changed_values = workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &changed_model_values,
                &changed_additional_values,
                &effective,
                &star,
            )
            .unwrap();
        for (index, input) in layout.inputs().iter().enumerate() {
            match input {
                ThresholdMultiplierInput::ModelParameter { index: 0 } => {
                    assert_eq!(changed_values.as_slice()[index].re, F(13.0));
                }
                ThresholdMultiplierInput::AdditionalParameter { index: 0 } => {
                    assert_eq!(changed_values.as_slice()[index].re, F(17.0));
                }
                _ => {}
            }
        }
        assert_eq!(workspace.kinematic_points.len(), 2);
    }

    #[test]
    fn workspace_reuses_exact_kinematic_contexts_despite_cache_id_aliases() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph threshold_multiplier_context_cache {
            ext [style=invis]
            v1 [num=1]
            v2 [num=1]
            edge [num=1 mass=0]
            ext -> v1 [id=0]
            v1 -> v2 [id=1]
            v1 -> v2 [id=2]
            v2 -> ext [id=3]
        })
        .unwrap();
        let layout = ThresholdMultiplierLayout::new(
            Vec::new(),
            Vec::new(),
            graph.loop_momentum_basis.ext_edges.len(),
            (0..graph.n_edges()).collect(),
            Vec::new(),
        )
        .unwrap();
        let make_sample = |x: f64, cache_id: usize| MomentumSample {
            sample: BareMomentumSample {
                loop_moms: vec![ThreeMomentum::new(F(x), F(2.0), F(3.0))]
                    .into_iter()
                    .collect(),
                dual_loop_moms: None,
                // Root transformations may either inherit or replace cache IDs. The multiplier
                // cache is therefore keyed by the exact non-dual kinematics it consumes.
                loop_mom_cache_id: cache_id,
                loop_mom_base_cache_id: 0,
                external_moms: (0..graph.loop_momentum_basis.ext_edges.len())
                    .map(|_| FourMomentum::from_args(F(10.0), F(0.0), F(0.0), F(0.0)))
                    .collect(),
                external_mom_cache_id: cache_id,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let base = make_sample(1.0, 0);
        let root = make_sample(4.0, 0);
        let root_with_new_ids = make_sample(4.0, 17);
        let other_root = make_sample(7.0, 0);
        let masses = (0..graph.n_edges()).map(|_| F(0.0)).collect();
        let mut workspace = ThresholdMultiplierInputWorkspace::new(&layout, F(0.0));

        let first = workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &[],
                &[],
                &base,
                &root,
            )
            .unwrap()
            .as_slice()
            .to_vec();
        assert_eq!(workspace.kinematic_points.len(), 2);
        let repeated = workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &[],
                &[],
                &base,
                &root,
            )
            .unwrap()
            .as_slice()
            .to_vec();
        assert_eq!(workspace.kinematic_points.len(), 2);
        assert_eq!(first, repeated);

        workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &[],
                &[],
                &base,
                &root_with_new_ids,
            )
            .unwrap();
        assert_eq!(workspace.kinematic_points.len(), 2);

        workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &[],
                &[],
                &root,
                &root,
            )
            .unwrap();
        assert_eq!(workspace.kinematic_points.len(), 2);
        workspace
            .bind(
                &layout,
                &graph.loop_momentum_basis,
                &masses,
                &[],
                &[],
                &base,
                &other_root,
            )
            .unwrap();
        assert_eq!(workspace.kinematic_points.len(), 3);

        let base_f128 = base.higher_precision();
        let root_f128 = root.higher_precision();
        assert_context_cache_reuse(
            &layout,
            &graph.loop_momentum_basis,
            graph.n_edges(),
            &base_f128,
            &root_f128,
        );
        assert_context_cache_reuse(
            &layout,
            &graph.loop_momentum_basis,
            graph.n_edges(),
            &base_f128.higher_precision(),
            &root_f128.higher_precision(),
        );
    }

    #[test]
    fn workspace_recomputes_complete_merged_views_before_caching() {
        test_initialise().unwrap();
        let graph: Graph = dot!(digraph threshold_multiplier_merged_cache {
            ext [style=invis]
            v1 [num=1]
            v2 [num=1]
            edge [num=1 mass=0]
            ext -> v1 [id=0]
            v1 -> v2 [id=1]
            v1 -> v2 [id=2]
            v1 -> v2 [id=3]
            v2 -> ext [id=4]
        })
        .unwrap();
        assert_eq!(graph.loop_momentum_basis.loop_edges.len(), 2);
        let energy_edges = graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>();
        let layout = ThresholdMultiplierLayout::new(
            Vec::new(),
            Vec::new(),
            graph.loop_momentum_basis.ext_edges.len(),
            (0..graph.n_edges()).collect(),
            vec![esurface(energy_edges, Vec::new())],
        )
        .unwrap();
        let external_momenta = (0..graph.loop_momentum_basis.ext_edges.len())
            .map(|_| FourMomentum::from_args(F(20.0), F(0.5), F(-0.25), F(0.75)))
            .collect::<ExternalFourMomenta<_>>();
        let make_sample = |loop_momenta: [[f64; 3]; 2]| MomentumSample {
            sample: BareMomentumSample {
                loop_moms: loop_momenta
                    .map(|[x, y, z]| ThreeMomentum::new(F(x), F(y), F(z)))
                    .into_iter()
                    .collect(),
                dual_loop_moms: None,
                loop_mom_cache_id: 0,
                loop_mom_base_cache_id: 0,
                external_moms: external_momenta.clone(),
                external_mom_cache_id: 0,
                external_mom_base_cache_id: 0,
                jacobian: F(1.0),
                orientation: None,
                parameterization_branch: None,
            },
        };
        let base = make_sample([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let left = make_sample([[11.0, 12.0, 13.0], [4.0, 5.0, 6.0]]);
        let right = make_sample([[1.0, 2.0, 3.0], [24.0, 25.0, 26.0]]);
        let merged = make_sample([[11.0, 12.0, 13.0], [24.0, 25.0, 26.0]]);
        let masses = (0..graph.n_edges()).map(|_| F(0.0)).collect();
        let mut workspace = ThresholdMultiplierInputWorkspace::new(&layout, F(0.0));
        let star_values = |workspace: &mut ThresholdMultiplierInputWorkspace<f64>,
                           star: &MomentumSample<f64>| {
            workspace
                .bind(
                    &layout,
                    &graph.loop_momentum_basis,
                    &masses,
                    &[],
                    &[],
                    &base,
                    star,
                )
                .unwrap()
                .as_slice()
                .to_vec()
        };
        let left_values = star_values(&mut workspace, &left);
        let right_values = star_values(&mut workspace, &right);
        let merged_values = star_values(&mut workspace, &merged);
        assert_eq!(workspace.kinematic_points.len(), 4);

        let external_spatial = merged
            .external_moms()
            .iter()
            .map(|momentum| momentum.spatial)
            .collect::<ExternalThreeMomenta<_>>();
        let expected_edges = layout
            .edges()
            .iter()
            .map(|&edge| {
                graph.loop_momentum_basis.edge_signatures[EdgeIndex(edge)]
                    .compute_momentum(merged.loop_moms(), &external_spatial)
            })
            .collect::<Vec<_>>();
        for (input_index, input) in layout.inputs().iter().enumerate() {
            match *input {
                ThresholdMultiplierInput::EdgeMomentum {
                    point: ThresholdMultiplierPoint::Star,
                    edge,
                    component,
                } => {
                    let edge_position = layout.edges().binary_search(&edge).unwrap();
                    let expected =
                        three_momentum_component(&expected_edges[edge_position], component)
                            .unwrap();
                    assert_eq!(merged_values[input_index].re, expected);
                }
                ThresholdMultiplierInput::Esurface {
                    point: ThresholdMultiplierPoint::Star,
                    esurface: 0,
                } => {
                    let expected = layout.esurfaces()[0]
                        .edges
                        .iter()
                        .map(|edge| {
                            let edge_position = layout.edges().binary_search(edge).unwrap();
                            expected_edges[edge_position].norm_squared().sqrt()
                        })
                        .fold(F(0.0), |sum, energy| sum + energy);
                    assert_eq!(merged_values[input_index].re, expected);
                }
                _ => {}
            }
        }
        assert!(layout.inputs().iter().enumerate().any(|(index, input)| {
            matches!(
                input,
                ThresholdMultiplierInput::EdgeMomentum {
                    point: ThresholdMultiplierPoint::Star,
                    ..
                }
            ) && merged_values[index] != left_values[index]
                && merged_values[index] != right_values[index]
        }));

        let repeated = star_values(&mut workspace, &merged);
        assert_eq!(workspace.kinematic_points.len(), 4);
        assert_eq!(repeated, merged_values);
    }
}
