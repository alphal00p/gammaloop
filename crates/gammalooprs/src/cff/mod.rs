use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
    hash::{Hash, Hasher},
    sync::Arc,
    time::{Duration, Instant},
};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::{EdgeIndex, Hedge},
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    symbol,
};

use crate::{
    cff::{
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, normalize_cut_edge_support_with_raised_edge_groups,
            normalize_three_d_expression_cut_support_with_raised_edge_groups,
            select_indexed_cff_residues,
        },
        generation::{ExactCffPreparationKey, PreparationKey, PreparationValue},
        orientations::GraphOrientation,
        surface::{GammaLoopLinearEnergyExpr, GammaLoopSurfaceCache, LinearEnergyExpr},
    },
    graph::{
        ExactUvSubLmbFrame, FeynmanGraph, FourDDenominator, Graph, GraphThreeDSource,
        LoopMomentumBasis,
        cuts::CutSet,
        get_cff_inverse_energy_product_impl,
        three_d_source::{ExactSourceEnergyMapper, ExactSourceMappingContext},
    },
    numerator::energy_degree::{
        EnergyPowerAnalyzer, EnergyPowerAssignmentPlan, PlannedEnergyExpression,
    },
    settings::global::OrientationPattern,
    utils::GS,
    uv::{
        Integrands,
        approx::{local_4d::CanonicalUvDenominatorClass, projected_4d::Local4dProjectionContext},
    },
};
use color_eyre::Result;
use three_dimensional_reps::{
    CffEnergyFactorOwnership, CffGlobalPrefactorSign, Generate3DExpressionOptions,
    GeneratedThreeDExpression,
};

pub mod orientations;
//pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;
mod vertex_set;
pub(crate) use vertex_set::VertexSet;

pub(crate) struct CFFOrientationTerm {
    pub(crate) expression: Atom,
    pub(crate) orientation: OrientationExpression,
    pub(crate) production_orientation_id: Option<OrientationID>,
}

pub struct CFFTerm {
    // Ordinary CFF maps retain production identity for direct-3D localization;
    // exact 4D maps instead remain source-local and share their parent mapper.
    pub(crate) orientations: Vec<CFFOrientationTerm>,
    exact_source_numerator: Option<Arc<PlannedExactSourceNumerator>>,
}

#[derive(PartialEq, Eq, Hash)]
struct ExactNumeratorTemplateBinding {
    context: ExactSourceMappingContext,
    assignment: Arc<EnergyPowerAssignmentPlan>,
}

/// Hash the complete binding once. Equality still checks the exact structure;
/// neither a hash collision nor a cache eviction changes the chosen assignment.
#[derive(Clone)]
pub(crate) struct ExactNumeratorTemplateKey {
    binding: Arc<ExactNumeratorTemplateBinding>,
    hash: u64,
}

impl PartialEq for ExactNumeratorTemplateKey {
    fn eq(&self, other: &Self) -> bool {
        self.hash == other.hash && self.binding == other.binding
    }
}

impl Eq for ExactNumeratorTemplateKey {}

impl Hash for ExactNumeratorTemplateKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash.hash(state);
    }
}

/// A small allocation binds memo rows to one immutable prepared template.
/// Keeping it alive prevents address reuse without retaining an evicted AST.
#[derive(Clone)]
struct ExactNumeratorTemplateIdentity(Arc<()>);

impl PartialEq for ExactNumeratorTemplateIdentity {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl Eq for ExactNumeratorTemplateIdentity {}

impl Hash for ExactNumeratorTemplateIdentity {
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state);
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub(crate) struct ExactNumeratorRowKey {
    template: ExactNumeratorTemplateIdentity,
    node: usize,
    samples: Vec<Atom>,
}

pub(crate) struct PlannedExactSourceNumerator {
    mapper: Arc<ExactSourceEnergyMapper>,
    binding: ExactNumeratorTemplateKey,
    identity: ExactNumeratorTemplateIdentity,
    template: PlannedEnergyExpression,
    parameters: Vec<Atom>,
    dependencies: Vec<Vec<usize>>,
    node_ends: Vec<usize>,
    constants: Vec<Option<Atom>>,
}

impl PlannedExactSourceNumerator {
    /// Localization can merge simple LTD factors into a raised pole. Bound the
    /// immutable assigned template only then, including fixed affine carriers
    /// which depend directly on source-loop energies rather than occurrences.
    fn residue_energy_degree_bounds(
        &self,
        loop_count: usize,
        edge_count: usize,
    ) -> Result<Vec<(usize, usize)>> {
        if self.parameters.len() < loop_count || self.parameters.len() - loop_count > edge_count {
            return Err(eyre::eyre!(
                "exact numerator template samples do not fit the residue energy maps"
            ));
        }
        let analyzer =
            EnergyPowerAnalyzer::for_scalar_parameters(self.parameters.iter().enumerate().map(
                |(index, parameter)| {
                    // Laurent selection appends loop slots after the unchanged
                    // occurrence namespace. Template preparation uses loops first.
                    let slot = if index < loop_count {
                        edge_count + index
                    } else {
                        index - loop_count
                    };
                    (parameter.clone(), EdgeIndex(slot))
                },
            ));
        let numerator = self
            .template
            .map(&mut |_, factor, _| Ok::<_, std::convert::Infallible>(factor.clone()))
            .unwrap();
        Ok(analyzer.analyze_atom(&numerator)?.into_generation_bounds())
    }

    /// Keep the certified ordinary numerator once and describe every residue
    /// through scalar coefficients of its full affine energy maps. In
    /// particular, sampling nodes multiply M rather than replacing it, and
    /// repeated denominator occurrences retain independent coefficient slots.
    fn prepare_parametric(
        &self,
        orientations: &[CFFOrientationTerm],
        scope: &Atom,
    ) -> Result<(Atom, Vec<Atom>, Vec<Vec<Atom>>)> {
        let coefficient = symbol!("gammalooprs::uv::numerator_coefficient"; Scalar);
        let mut captured = false;
        let _ = self.template.map(&mut |_, factor, _| {
            let _ = factor.replace_map(|view, _, _| {
                if let symbolica::atom::AtomView::Fun(function) = view
                    && function.get_symbol() == coefficient
                    && function.get_nargs() == 2
                    && function.get(0) == scope.as_view()
                {
                    captured = true;
                }
            });
            Ok::<_, std::convert::Infallible>(Atom::one())
        });
        if captured {
            return Err(eyre::eyre!(
                "exact numerator coefficient scope {scope} is already in use"
            ));
        }
        let required = &self.dependencies[0];
        let Some(first) = orientations.first() else {
            return Err(eyre::eyre!("exact numerator family has no residue rows"));
        };
        let loop_count = first.orientation.loop_energy_map.len();
        let edge_count = first.orientation.edge_energy_map.len();
        let mut support = BTreeMap::new();
        for (row, orientation) in orientations.iter().enumerate() {
            let orientation = &orientation.orientation;
            if orientation.loop_energy_map.len() != loop_count
                || orientation.edge_energy_map.len() != edge_count
            {
                return Err(eyre::eyre!(
                    "exact numerator family has inconsistent affine map lengths"
                ));
            }
            // Validate the source's exact coordinate and occurrence bindings
            // before preparing even an identically zero sample row.
            self.mapper.sample_energies(
                &orientation.loop_energy_map,
                &orientation.edge_energy_map,
                required,
            )?;
            for &slot in required {
                let energy = if slot < loop_count {
                    &orientation.loop_energy_map[slot]
                } else {
                    &orientation.edge_energy_map[slot - loop_count]
                }
                .clone()
                .canonical();
                let coefficients = energy
                    .internal_terms
                    .iter()
                    .map(|(edge, value)| (0, usize::from(*edge), value))
                    .chain(
                        energy
                            .external_terms
                            .iter()
                            .map(|(edge, value)| (1, usize::from(*edge), value)),
                    )
                    .chain([
                        (2, 0, &energy.uniform_scale_coeff),
                        (3, 0, &energy.constant),
                    ]);
                for (kind, edge, value) in coefficients {
                    let value = Atom::num(value.clone());
                    if !value.is_zero() {
                        support
                            .entry((slot, kind, edge))
                            .or_insert_with(|| vec![Atom::Zero; orientations.len()])[row] = value;
                    }
                }
            }
        }
        let mut parameters = Vec::with_capacity(support.len());
        let mut rows = vec![Vec::with_capacity(support.len()); orientations.len()];
        let mut samples = required
            .iter()
            .copied()
            .map(|slot| (slot, Atom::Zero))
            .collect::<BTreeMap<_, _>>();
        for ((slot, kind, edge), values) in support {
            let mut unit = LinearEnergyExpr::zero();
            match kind {
                0 => unit.internal_terms.push((EdgeIndex(edge), 1.into())),
                1 => unit.external_terms.push((EdgeIndex(edge), 1.into())),
                2 => unit.uniform_scale_coeff = 1.into(),
                3 => unit.constant = 1.into(),
                _ => unreachable!(),
            }
            // The existing source mapper remains the only authority for exact
            // OSE aliases and affine external shifts, including UV class IDs.
            let basis = unit
                .to_atom_gs(&[])
                .replace_multiple(self.mapper.exact_ose_replacements());
            // The scope is part of the formal itself, so expression-local
            // persistence retains it even after this family is recomposed.
            let parameter = coefficient.call_args([scope.clone(), Atom::num(parameters.len())]);
            *samples.get_mut(&slot).unwrap() += &parameter * basis;
            parameters.push(parameter);
            for (arguments, value) in rows.iter_mut().zip(values) {
                arguments.push(value);
            }
        }
        let mapped = self.map_template(&self.template, &samples, &mut 0, &mut None);
        Ok((
            self.mapper.set_inactive_loop_energies_to_zero(mapped),
            parameters,
            rows,
        ))
    }

    /// Return the certified template and its separate build/cache costs so
    /// candidate selection can charge discarded preparation without duplication.
    fn prepare(
        mapper: Arc<ExactSourceEnergyMapper>,
        assignment: Arc<EnergyPowerAssignmentPlan>,
        context: Option<&mut Local4dProjectionContext>,
    ) -> Result<(Arc<Self>, Duration, Duration)> {
        let started = Instant::now();
        let binding = ExactNumeratorTemplateBinding {
            context: mapper.mapping_context()?,
            assignment,
        };
        let mut hash = std::hash::DefaultHasher::new();
        binding.hash(&mut hash);
        let key = ExactNumeratorTemplateKey {
            binding: Arc::new(binding),
            hash: hash.finish(),
        };
        let cache_key = PreparationKey::Numerator(key.clone());
        let mut context = context;
        if let Some(context) = context.as_deref_mut()
            && let Some(PreparationValue::Numerator(template)) =
                context.preparations.get(&cache_key)
        {
            let template = Arc::clone(template);
            drop(cache_key);
            drop(key);
            let cache_time = started.elapsed();
            context.preparations.cache_time += cache_time;
            crate::debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "numerator_template",
                cache_hit = true,
                elapsed_ms = cache_time.as_secs_f64() * 1000.0,
                cache_work_ms = cache_time.as_secs_f64() * 1000.0,
                template_build_ms = 0.0,
                "Reused immutable exact-source numerator template"
            );
            return Ok((template, Duration::ZERO, cache_time));
        }
        let build_started = Instant::now();
        let (template, parameters) = mapper.prepare_numerator(&key.binding.assignment)?;
        let mut prepared = Self {
            mapper,
            binding: key.clone(),
            identity: ExactNumeratorTemplateIdentity(Arc::new(())),
            template,
            parameters,
            dependencies: Vec::new(),
            node_ends: Vec::new(),
            constants: Vec::new(),
        };
        Self::parameter_dependencies(
            &prepared.template,
            &prepared.parameters,
            &mut prepared.dependencies,
            &mut prepared.node_ends,
            &mut prepared.constants,
        );
        let build_time = build_started.elapsed();
        let bytes = prepared.accounted_bytes()?;
        let prepared = Arc::new(prepared);
        if let Some(context) = context.as_deref_mut() {
            context.numerator_template_builds += 1;
            context.preparations.insert(
                cache_key,
                PreparationValue::Numerator(Arc::clone(&prepared)),
                bytes,
            );
        }
        let elapsed = started.elapsed();
        let cache_time = elapsed.saturating_sub(build_time);
        if let Some(context) = context {
            context.preparations.cache_time += cache_time;
        }
        crate::debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "numerator_template",
            cache_hit = false,
            template_nodes = prepared.dependencies.len(),
            parameters = prepared.parameters.len(),
            accounted_bytes = bytes,
            elapsed_ms = elapsed.as_secs_f64() * 1000.0,
            cache_work_ms = cache_time.as_secs_f64() * 1000.0,
            template_build_ms = build_time.as_secs_f64() * 1000.0,
            "Prepared immutable exact-source numerator template"
        );
        Ok((prepared, build_time, cache_time))
    }

    fn parameter_dependencies(
        expression: &PlannedEnergyExpression,
        parameters: &[Atom],
        nodes: &mut Vec<Vec<usize>>,
        node_ends: &mut Vec<usize>,
        constants: &mut Vec<Option<Atom>>,
    ) -> usize {
        let node = nodes.len();
        nodes.push(Vec::new());
        node_ends.push(0);
        constants.push(None);
        let (dependencies, constant) = match expression {
            PlannedEnergyExpression::Factor { expression, .. } => {
                let mut used = BTreeSet::new();
                let _ = expression.replace_map(|view, _, _| {
                    if let Some(index) = parameters
                        .iter()
                        .position(|parameter| parameter.as_view() == view)
                    {
                        used.insert(index);
                    }
                });
                let constant = used.is_empty().then(|| expression.clone());
                (used, constant)
            }
            PlannedEnergyExpression::Add(children)
            | PlannedEnergyExpression::Mul(children)
            | PlannedEnergyExpression::MultilinearFunction {
                arguments: children,
                ..
            } => {
                let children = children
                    .iter()
                    .map(|child| {
                        Self::parameter_dependencies(child, parameters, nodes, node_ends, constants)
                    })
                    .collect::<Vec<_>>();
                let dependencies = children
                    .iter()
                    .flat_map(|child| nodes[*child].iter().copied())
                    .collect::<BTreeSet<_>>();
                let constant = dependencies.is_empty().then(|| match expression {
                    PlannedEnergyExpression::Add(_) => {
                        children.iter().fold(Atom::Zero, |sum, child| {
                            sum + constants[*child].as_ref().unwrap()
                        })
                    }
                    PlannedEnergyExpression::Mul(_) => {
                        children.iter().fold(Atom::one(), |product, child| {
                            product * constants[*child].as_ref().unwrap()
                        })
                    }
                    PlannedEnergyExpression::MultilinearFunction { symbol, .. } => {
                        let mut builder = FunctionBuilder::new(*symbol);
                        for child in &children {
                            builder =
                                builder.add_arg(constants[*child].as_ref().unwrap().as_view());
                        }
                        builder.finish()
                    }
                    _ => unreachable!(),
                });
                (dependencies, constant)
            }
            PlannedEnergyExpression::Repeat { base, exponent } => {
                let child =
                    Self::parameter_dependencies(base, parameters, nodes, node_ends, constants);
                (
                    nodes[child].iter().copied().collect(),
                    constants[child]
                        .as_ref()
                        .map(|constant| constant.pow(*exponent as u64)),
                )
            }
        };
        nodes[node] = dependencies.iter().copied().collect();
        node_ends[node] = nodes.len();
        constants[node] = constant;
        node
    }

    fn accounted_bytes(&self) -> Result<usize> {
        // The binding and the owned mapper both retain source coordinates and
        // replacements. Arc key copies share the binding payload.
        Ok(std::mem::size_of::<Self>()
            + 128
            + 2 * std::mem::size_of::<PreparationKey>()
            + std::mem::size_of::<PreparationValue>()
            + self.binding.binding.context.accounted_bytes()
            + self.mapper.accounted_bytes()?
            + self.binding.binding.assignment.accounted_bytes()
            + self.template.accounted_bytes()
            + self.parameters.capacity() * std::mem::size_of::<Atom>()
            + self
                .parameters
                .iter()
                .map(|atom| atom.as_view().get_byte_size())
                .sum::<usize>()
            + self.dependencies.capacity() * std::mem::size_of::<Vec<usize>>()
            + self
                .dependencies
                .iter()
                .map(|indices| indices.capacity() * std::mem::size_of::<usize>())
                .sum::<usize>()
            + self.node_ends.capacity() * std::mem::size_of::<usize>()
            + self.constants.capacity() * std::mem::size_of::<Option<Atom>>()
            + self
                .constants
                .iter()
                .flatten()
                .map(|atom| atom.as_view().get_byte_size())
                .sum::<usize>())
    }

    #[cfg(test)]
    fn sample(
        &self,
        orientation: &OrientationExpression,
        context: Option<&mut Local4dProjectionContext>,
    ) -> Result<Atom> {
        let started = Instant::now();
        let required = &self.dependencies[0];
        let values = self.mapper.sample_energies(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            required,
        )?;
        let samples = required
            .iter()
            .copied()
            .zip(values)
            .collect::<BTreeMap<_, _>>();
        let mut cache = context.map(|context| &mut context.numerator_rows);
        let mapped = self.map_template(&self.template, &samples, &mut 0, &mut cache);
        let mapped = self.mapper.set_inactive_loop_energies_to_zero(mapped);
        crate::debug_tags!(#generation, #uv, #local, #four_d, #profile;
            stage = "numerator_row_mapping",
            required_energies = required.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Mapped one exact-source residue sample"
        );
        Ok(mapped)
    }

    fn map_template(
        &self,
        expression: &PlannedEnergyExpression,
        samples: &BTreeMap<usize, Atom>,
        next_node: &mut usize,
        cache: &mut Option<&mut generation::GenerationCache<ExactNumeratorRowKey, Atom>>,
    ) -> Atom {
        let node = *next_node;
        *next_node += 1;
        if let Some(constant) = &self.constants[node] {
            *next_node = self.node_ends[node];
            return constant.clone();
        }
        let lookup_started = cache.as_ref().map(|_| Instant::now());
        let key = cache.as_ref().map(|_| ExactNumeratorRowKey {
            template: self.identity.clone(),
            node,
            samples: self.dependencies[node]
                .iter()
                .map(|index| samples[index].clone())
                .collect(),
        });
        if let Some(cache) = cache.as_deref_mut() {
            let value = cache.get(key.as_ref().unwrap()).cloned();
            if let Some(value) = value {
                *next_node = self.node_ends[node];
                drop(key);
                cache.cache_time += lookup_started.unwrap().elapsed();
                return value;
            }
            cache.cache_time += lookup_started.unwrap().elapsed();
        }
        let mapped = match expression {
            PlannedEnergyExpression::Factor { expression, .. } => {
                let dependencies = &self.dependencies[node];
                if dependencies.is_empty() {
                    return expression.clone();
                }
                expression.replace_map(|view, _, output| {
                    for index in dependencies {
                        if self.parameters[*index].as_view() == view {
                            **output = samples[index].clone();
                            return;
                        }
                    }
                })
            }
            PlannedEnergyExpression::Add(children) => {
                children.iter().fold(Atom::Zero, |sum, child| {
                    sum + self.map_template(child, samples, next_node, cache)
                })
            }
            PlannedEnergyExpression::Mul(children) => {
                children.iter().fold(Atom::one(), |product, child| {
                    product * self.map_template(child, samples, next_node, cache)
                })
            }
            PlannedEnergyExpression::MultilinearFunction { symbol, arguments } => {
                let mut builder = FunctionBuilder::new(*symbol);
                for child in arguments {
                    builder = builder.add_arg(self.map_template(child, samples, next_node, cache));
                }
                builder.finish()
            }
            PlannedEnergyExpression::Repeat { base, exponent } => self
                .map_template(base, samples, next_node, cache)
                .pow(*exponent as u64),
        };
        if let Some(cache) = cache.as_deref_mut() {
            let insertion_started = Instant::now();
            let key = key.unwrap();
            let bytes = 2
                * (std::mem::size_of::<ExactNumeratorRowKey>()
                    + key.samples.capacity() * std::mem::size_of::<Atom>()
                    + key
                        .samples
                        .iter()
                        .map(|atom| atom.as_view().get_byte_size())
                        .sum::<usize>())
                + std::mem::size_of::<Atom>()
                + mapped.as_view().get_byte_size()
                + 96;
            if cache.can_retain(bytes) {
                cache.insert(key, mapped.clone(), bytes);
            } else {
                drop(key);
            }
            cache.cache_time += insertion_started.elapsed();
        }
        mapped
    }
}

impl CFFTerm {
    pub(crate) fn prepare_exact_source_numerator(
        &self,
        scope: &Atom,
    ) -> Result<(Atom, Vec<Atom>, Vec<Vec<Atom>>)> {
        self.exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator plan"))?
            .prepare_parametric(&self.orientations, scope)
    }

    #[cfg(test)]
    pub(crate) fn map_exact_source_atom(
        &self,
        orientation: &OrientationExpression,
        atom: &Atom,
    ) -> Result<Atom> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator map"))?;
        planned.mapper.map_numerator(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            atom,
        )
    }

    #[cfg(test)]
    pub(crate) fn map_exact_source_physical_loop_lift_energies(
        &self,
        orientation: &OrientationExpression,
        physical_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Result<Vec<(EdgeIndex, Atom)>> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator map"))?;
        planned.mapper.map_physical_owner_loop_lift_energies(
            &orientation.loop_energy_map,
            &orientation.edge_energy_map,
            physical_edges,
        )
    }

    #[cfg(test)]
    pub(crate) fn map_exact_source_numerator(
        &self,
        orientation: &OrientationExpression,
        context: Option<&mut Local4dProjectionContext>,
    ) -> Result<Atom> {
        let planned = self
            .exact_source_numerator
            .as_ref()
            .ok_or_else(|| eyre::eyre!("ordinary CFF term has no exact-source numerator plan"))?;
        planned.sample(orientation, context)
    }

    pub fn expression_with_selectors(&self) -> Atom {
        self.orientations
            .iter()
            .map(|term| {
                let selector = term.production_orientation_id.map_or_else(
                    || term.orientation.data.orientation.orientation_thetas(),
                    OrientationID::atom,
                );
                term.expression.clone() * selector
            })
            .reduce(|left, right| left + right)
            .unwrap_or(Atom::Zero)
    }
}

#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Encode, Decode, Serialize, Deserialize,
)]
// This describes the combinations of residues that are selected.
pub struct CutCFFIndex {
    pub left_threshold_order: Option<usize>,
    pub right_threshold_order: Option<usize>,
    pub lu_cut_order: Option<usize>,
}

impl CutCFFIndex {
    pub fn new_all_none() -> Self {
        Self {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: None,
        }
    }
}

impl Display for CutCFFIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut parts = vec![];
        if let Some(order) = self.lu_cut_order {
            parts.push(format!("lu_cut_{}", order));
        }

        if let Some(order) = self.left_threshold_order {
            parts.push(format!("left_th_{}", order));
        }

        if let Some(order) = self.right_threshold_order {
            parts.push(format!("right_th_{}", order));
        }

        if parts.is_empty() {
            write!(f, "")
        } else {
            write!(f, "{}", parts.join("_"))
        }
    }
}

/// Namespace used by the energy-degree bounds supplied to a CFF source.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum CffEnergyBoundSourceKind {
    /// Source IDs are physical graph-edge IDs.
    PhysicalGraph,
    /// Source IDs are exact four-dimensional denominator-occurrence IDs.
    ExactFourD,
}

/// Transient diagnostic record for one CFF generation source.
///
/// Exact four-dimensional Taylor terms can contain several denominator
/// occurrences owned by one physical edge. The CFF generator must receive the
/// selected certified assignment in that occurrence namespace, while the physical-parent
/// bounds remain useful for checking the numerator analysis. Keeping both
/// explicitly prevents occurrence IDs from being mistaken for physical edges.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct CffEnergyDegreeBoundReport {
    pub source_kind: CffEnergyBoundSourceKind,
    pub physical_parent_bounds: Vec<(usize, usize)>,
    pub assigned_cff_source_bounds: Vec<(usize, usize)>,
}

pub struct CutCFF {
    pub terms: BTreeMap<CutCFFIndex, CFFTerm>,
    pub(crate) energy_degree_bound_report: CffEnergyDegreeBoundReport,
    // Terms retain the shared-core contour convention for exact-CFF users and
    // normalization oracles. GammaLoop production consumes this typed bridge
    // only when localizing either the direct-3D or exact-4D route.
    production_prefactor_bridge: CffGlobalPrefactorSign,
}

impl CutCFF {
    /// Convert a generated CFF into GammaLoop's scalar-denominator convention.
    ///
    /// Together with the energy factors below this yields the signed
    /// `dq0/(2*pi*i)` contour. Physical amplitudes restore `dq0/(2*pi)`
    /// through a separate `i` per integrated loop energy.
    ///
    /// `three-dimensional-reps` writes every source with positive local
    /// `1/(2E_i)` factors. For an ordinary CFF, surface conversion removes
    /// those factors and `Graph::cff` restores the historical global
    /// `1/prod(-2E_i)` product. Its `(-1)^N` is therefore already explicit in
    /// the converted expression. A generalized CFF must retain its factors on
    /// each variant because contact terms need not have the same half-edge
    /// support; in that case the same `(-1)^N` conversion is supplied here.
    ///
    /// An ordinary component keeps its established core conversion. A
    /// generalized component instead uses the numerator-bound-independent
    /// scalar-denominator frame. Its generalized core sign is already encoded
    /// in the raw Laurent functional and must not be multiplied again: doing
    /// so would let an unused numerator-rank allowance change the value of an
    /// otherwise identical scalar integrand.
    fn gamma_loop_prefactor_conversion<E, H>(
        generated: &GeneratedThreeDExpression<E, H>,
    ) -> CffGlobalPrefactorSign {
        // LTD already carries the signed dq0/(2*pi*i) contour, including its
        // basis-local positive energy factors. Only CFF needs the source-frame
        // bridge below; the physical i-per-loop normalization is shared.
        if generated.representation == three_dimensional_reps::RepresentationMode::Ltd {
            return CffGlobalPrefactorSign::default();
        }
        let retained_positive_energy_factors =
            generated.energy_factor_ownership == CffEnergyFactorOwnership::VariantLocal;
        generated.energy_factor_components.iter().fold(
            CffGlobalPrefactorSign::default(),
            |conversion, component| {
                let source_frame = match component.ownership {
                    CffEnergyFactorOwnership::GlobalSourceProduct => {
                        component.core_global_prefactor_sign
                    }
                    CffEnergyFactorOwnership::VariantLocal => {
                        component.denominator_only_global_prefactor_sign
                    }
                };
                conversion
                    .product(CffGlobalPrefactorSign::from_exponent(
                        component.internal_edge_ids.len()
                            * usize::from(retained_positive_energy_factors),
                    ))
                    .product(source_frame)
            },
        )
    }

    pub(crate) const fn production_prefactor_factor(&self) -> i64 {
        self.production_prefactor_bridge.factor()
    }

    pub fn expression_with_selectors(&self) -> Integrands {
        let production_prefactor = Atom::num(self.production_prefactor_factor());
        self.terms
            .iter()
            .map(|(index, term)| {
                (
                    *index,
                    term.expression_with_selectors() * &production_prefactor,
                )
            })
            .collect()
    }
}

impl Graph {
    pub(crate) fn cff_from_production_expression(
        &self,
        production: &GeneratedThreeDExpression<esurface::Esurface, hsurface::Hsurface>,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CutCFF> {
        let production_prefactor_bridge = CutCFF::gamma_loop_prefactor_conversion(production);
        let contract_subgraph = self.tree_edges.subtract(&self.initial_state_cut);
        let contract_edges = self
            .iter_edges_of(&contract_subgraph)
            .filter_map(|(pair, edge_id, _)| pair.is_paired().then_some(edge_id))
            .collect::<Vec<_>>();
        let mut cff = production.expression.clone();
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut cff,
            &self.get_raised_edge_groups(),
        );
        let residues = select_indexed_cff_residues(
            cff,
            cutset,
            production.representation,
            &self.get_raised_edge_groups(),
            || {
                if !production.source_energy_degree_bounds.is_empty() {
                    return Ok(production.source_energy_degree_bounds.clone());
                }
                // Simple LTD generation needs no degree information. Localizing
                // several physical surfaces can merge its remaining factors into
                // a raised pole, whose numerator jet needs bounds only now.
                let numerator = self.production_numerator_atom_for_full_3d_expression();
                Ok(self
                    .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                        [&numerator],
                        self.iter_edges_of(&self.initial_state_cut)
                            .chain(self.iter_edges_of(&self.tree_edges))
                            .map(|(_, edge, _)| edge),
                        1,
                    )?)
            },
        )?;
        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);
        let cff_loop_number = self
            .get_loop_number()
            .saturating_sub(self.cyclotomatic_number(&contract_subgraph));
        let cff_phase = Atom::i().pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let cff_energy_factor = match production.energy_factor_ownership {
            CffEnergyFactorOwnership::GlobalSourceProduct => {
                get_cff_inverse_energy_product_impl(self, &graph_without_is_cut, &contract_edges)
            }
            CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
        };

        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: Vec::new(),
                exact_source_numerator: None,
            };
            for (orientation_index, orientation) in expr.orientations.into_iter().enumerate() {
                if !orientation_pattern.filter_orientation(&orientation.data.orientation) {
                    continue;
                }
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    * &cff_energy_factor
                    * &cff_normalization;
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation,
                    production_orientation_id: (production.representation
                        == three_dimensional_reps::RepresentationMode::Cff)
                        .then_some(OrientationID(orientation_index)),
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        Ok(CutCFF {
            terms,
            energy_degree_bound_report: CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                physical_parent_bounds: production.source_energy_degree_bounds.clone(),
                assigned_cff_source_bounds: production.source_energy_degree_bounds.clone(),
            },
            production_prefactor_bridge,
        })
    }

    #[cfg(test)]
    pub(crate) fn cff_from_4d_denominators(
        &mut self,
        denominators: &[FourDDenominator],
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_edges(
            denominators,
            [],
            cutset,
            options,
            analysis_numerator,
            None,
        )
    }

    #[cfg(test)]
    pub(crate) fn cff_from_4d_denominators_in_uv_edges(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        context: Option<&mut Local4dProjectionContext>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_edges_and_boundaries(
            denominators,
            uv_edges,
            [],
            cutset,
            options,
            analysis_numerator,
            context,
        )
    }

    /// Generate an exact CFF while retaining the crown of a non-vacuum UV
    /// source. Vacuum Taylor terms use the boundary-free wrapper above.
    #[cfg(test)]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cff_from_4d_denominators_in_uv_edges_and_boundaries(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        context: Option<&mut Local4dProjectionContext>,
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_coordinates(
            denominators,
            uv_edges,
            uv_boundary_hedges,
            None,
            cutset,
            options,
            analysis_numerator,
            context,
            &[],
        )
    }

    /// Generate a proper-subgraph exact CFF in its sub-LMB coordinates. Crown
    /// momenta are external in this source and remain available to the later
    /// outer CFF instead of being treated as inactive child loop energies.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn cff_from_4d_denominators_in_uv_sub_lmb(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        sub_lmb: &LoopMomentumBasis,
        frame: ExactUvSubLmbFrame,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        context: Option<&mut Local4dProjectionContext>,
        classes: &[CanonicalUvDenominatorClass],
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        self.cff_from_4d_denominators_in_uv_coordinates(
            denominators,
            uv_edges,
            uv_boundary_hedges,
            Some((sub_lmb, frame)),
            cutset,
            options,
            analysis_numerator,
            context,
            classes,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn cff_from_4d_denominators_in_uv_coordinates(
        &mut self,
        denominators: &[FourDDenominator],
        uv_edges: impl IntoIterator<Item = EdgeIndex>,
        uv_boundary_hedges: impl IntoIterator<Item = Hedge>,
        coordinates: Option<(&LoopMomentumBasis, ExactUvSubLmbFrame)>,
        cutset: &CutSet,
        options: &Generate3DExpressionOptions,
        analysis_numerator: &Atom,
        context: Option<&mut Local4dProjectionContext>,
        classes: &[CanonicalUvDenominatorClass],
    ) -> Result<(CutCFF, linnet::half_edge::subgraph::SuBitGraph)> {
        let (
            generated,
            physical_surfaces,
            physical_energy_edges,
            physical_ose_coordinates,
            exact_source_numerator,
            inverse_energy_product,
            cff_loop_number,
            contract_subgraph,
            energy_degree_bound_report,
            physical_cut_support_edges,
            production_prefactor_bridge,
        ) = {
            let mut context = context;
            let source_started = Instant::now();
            let key = Arc::new(ExactCffPreparationKey {
                denominators: denominators.to_vec(),
                uv_edges: uv_edges
                    .into_iter()
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
                boundary_hedges: uv_boundary_hedges
                    .into_iter()
                    .collect::<BTreeSet<_>>()
                    .into_iter()
                    .collect(),
                coordinates: coordinates.map(|(lmb, frame)| (lmb.clone(), frame)),
                classes: classes.to_vec(),
                numerator: analysis_numerator.clone(),
                options: options.clone(),
            });
            let cache_key = PreparationKey::Source(Arc::clone(&key));
            let cached = context
                .as_deref_mut()
                .and_then(|context| context.preparations.get(&cache_key))
                .map(|entry| match entry {
                    PreparationValue::Source(preparation) => Arc::clone(preparation),
                    _ => unreachable!("preparation key and value kinds agree"),
                });
            let cache_hit = cached.is_some();
            let mut source_build_time = Duration::ZERO;
            let mut analysis_preparation_time = Duration::ZERO;
            let preparation = if let Some(preparation) = cached {
                drop(cache_key);
                drop(key);
                preparation
            } else {
                let reconstruction_started = Instant::now();
                let source = if let Some((sub_lmb, frame)) = &key.coordinates {
                    GraphThreeDSource::from_exact_denominators_in_uv_sub_lmb(
                        self,
                        &key.denominators,
                        key.uv_edges.iter().copied(),
                        key.boundary_hedges.iter().copied(),
                        sub_lmb,
                        *frame,
                    )?
                } else {
                    GraphThreeDSource::from_exact_denominators_in_uv_edges_and_boundaries(
                        self,
                        &key.denominators,
                        key.uv_edges.iter().copied(),
                        key.boundary_hedges.iter().copied(),
                    )?
                };
                source_build_time = reconstruction_started.elapsed();
                if let Some(context) = context.as_deref_mut() {
                    context.shared_preparation_time += source_build_time;
                }
                // Capacity analysis and certified assignment planning below
                // depend on the requested representation and keep its timing.
                let preparation_started = Instant::now();
                let preparation = Arc::new(self.prepare_3d_expression_for_4d_term(
                    &source,
                    options,
                    analysis_numerator,
                    classes,
                )?);
                analysis_preparation_time = preparation_started.elapsed();
                if let Some(context) = context.as_deref_mut() {
                    let bytes = key.accounted_bytes()
                        + preparation.accounted_bytes()?
                        + 2 * std::mem::size_of::<PreparationKey>()
                        + std::mem::size_of::<PreparationValue>()
                        + 128;
                    context.source_preparation_builds += 1;
                    context.preparations.insert(
                        cache_key,
                        PreparationValue::Source(Arc::clone(&preparation)),
                        bytes,
                    );
                }
                drop(source);
                drop(key);
                preparation
            };
            let request_time = source_started.elapsed();
            let cache_time =
                request_time.saturating_sub(source_build_time + analysis_preparation_time);
            if let Some(context) = context.as_deref_mut() {
                context.preparations.cache_time += cache_time;
            }
            crate::debug_tags!(#generation, #uv, #local, #four_d, #profile;
                stage = "source_reconstruction",
                graph = %self.name,
                denominator_occurrences = denominators.len(),
                classes = classes.len(),
                cache_hit,
                elapsed_ms = source_build_time.as_secs_f64() * 1000.0,
                analysis_preparation_ms = analysis_preparation_time.as_secs_f64() * 1000.0,
                request_elapsed_ms = request_time.as_secs_f64() * 1000.0,
                cache_work_ms = cache_time.as_secs_f64() * 1000.0,
                "Prepared certified exact denominator source incidence and numerator analysis"
            );
            let (generated, exact_source_numerator, _, energy_degree_bound_report) =
                self.generate_3d_expression_for_4d_term(&preparation, context)?;
            let physical_surfaces = generated
                .expression
                .surfaces
                .linear_surface_cache
                .iter()
                .map(|surface| {
                    surface
                        .expression
                        .internal_terms
                        .iter()
                        .all(|(edge, _)| {
                            preparation
                                .physical_surface_edges
                                .contains(&usize::from(*edge))
                        })
                        .then(|| {
                            let mut surface = surface.clone();
                            surface.expression = surface.expression.remap_energy_edges(
                                &preparation.physical_energy_edges.internal,
                                &BTreeMap::new(),
                            );
                            surface
                        })
                })
                .collect::<Vec<_>>();
            let physical_energy_edges = preparation.physical_energy_edges.clone();
            let physical_ose_coordinates = physical_energy_edges
                .internal
                .iter()
                .filter(|(occurrence, _)| preparation.physical_surface_edges.contains(occurrence))
                .map(|(&occurrence, &physical)| (occurrence, physical))
                .collect::<BTreeMap<_, _>>();
            let production_prefactor_bridge = CutCFF::gamma_loop_prefactor_conversion(&generated);
            (
                generated,
                physical_surfaces,
                physical_energy_edges,
                physical_ose_coordinates,
                exact_source_numerator,
                preparation.inverse_energy_product.clone(),
                preparation.active_loop_count,
                preparation.contract_subgraph.clone(),
                energy_degree_bound_report,
                preparation.physical_cut_support_edges.clone(),
                production_prefactor_bridge,
            )
        };
        let (generated, surface_ownership) =
            self.convert_4d_expression_surfaces(generated, &physical_surfaces)?;
        let energy_factor_ownership = generated.energy_factor_ownership;
        // Component metadata records the precise typed convention consumed by
        // the conversion above; no incidence or momentum-sign reconstruction
        // is performed after generation.
        crate::debug_tags!(#generation, #cff, #inspect;
            denominator_count = denominators.len(),
            production_prefactor_bridge = ?production_prefactor_bridge,
            aggregate_ownership = ?energy_factor_ownership,
            components = ?generated.energy_factor_components,
            physical_energy_edges = ?physical_energy_edges,
            surface_ownership = ?surface_ownership,
            "Exact CFF energy-factor component metadata: context={:?}, bridge={}, ownership={:?}, components={:?}, physical_energy_edges={:?}",
            options.cff_generation_context,
            production_prefactor_bridge.factor(),
            energy_factor_ownership,
            generated.energy_factor_components,
            physical_energy_edges,
        );
        let mut cff = generated.expression;
        if options.representation == three_dimensional_reps::RepresentationMode::Ltd {
            // The Laurent chart uses physical OSEs. Carry precisely the
            // certified physical occurrence coordinates into every value it
            // differentiates, retaining independent numerator argument slots
            // and nonphysical UV energies in their exact-source namespace.
            for orientation in &mut cff.orientations {
                orientation.remap_on_shell_energy_values(&physical_ose_coordinates);
            }
            for surface in &mut cff.surfaces.linear_surface_cache {
                surface.expression =
                    std::mem::replace(&mut surface.expression, LinearEnergyExpr::zero())
                        .remap_internal_edges(&physical_ose_coordinates);
            }
        }
        // Residue support belongs to physical Cutkosky alternatives. Its
        // provenance projection is separate from energy-coordinate transport:
        // numerator argument slots remain occurrence-local in both modes.
        let physical_edges = |edge: linnet::half_edge::involution::EdgeIndex| {
            let edge_id = usize::from(edge);
            if edge_id < physical_energy_edges.orientation_edge_count {
                return Ok(vec![edge]);
            }
            physical_cut_support_edges
                .get(&edge_id)
                .cloned()
                .ok_or_else(|| {
                    eyre::eyre!("exact CFF cut support contains unmapped occurrence edge {edge_id}")
                })
        };
        let raised_edge_groups = self.get_raised_edge_groups();
        let retain_physical_support_with_raised_representatives = |support: &mut Vec<
            linnet::half_edge::involution::EdgeIndex,
        >| {
            let representatives =
                normalize_cut_edge_support_with_raised_edge_groups(support, &raised_edge_groups);
            support.extend(representatives);
            support.sort_unstable();
            support.dedup();
        };
        for orientation in cff.orientations.iter_mut() {
            for variant in &mut orientation.variants {
                let mut denominator_edges = variant
                    .denominator_edges
                    .iter()
                    .copied()
                    .map(&physical_edges)
                    .collect::<Result<Vec<_>>>()?
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>();
                retain_physical_support_with_raised_representatives(&mut denominator_edges);
                variant.denominator_edges = denominator_edges;
                variant.denominator_edge_support_signs =
                    std::mem::take(&mut variant.denominator_edge_support_signs)
                        .into_iter()
                        .try_fold(BTreeMap::new(), |mut mapped, (support, sign)| {
                            let mut support = support
                                .into_iter()
                                .map(&physical_edges)
                                .collect::<Result<Vec<_>>>()?
                                .into_iter()
                                .flatten()
                                .collect::<Vec<_>>();
                            retain_physical_support_with_raised_representatives(&mut support);
                            *mapped.entry(support).or_insert(1) *= sign;
                            Ok::<_, color_eyre::Report>(mapped)
                        })?;
            }
        }
        let (loop_count, edge_count) = cff
            .orientations
            .first()
            .map(|orientation| {
                (
                    orientation.loop_energy_map.len(),
                    orientation.edge_energy_map.len(),
                )
            })
            .unwrap_or_default();
        let residues = select_indexed_cff_residues(
            cff,
            cutset,
            options.representation,
            &self.get_raised_edge_groups(),
            || exact_source_numerator.residue_energy_degree_bounds(loop_count, edge_count),
        )?;
        let cff_phase = Atom::i().pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let mut terms = BTreeMap::new();
        for (cut_cff_index, expr) in residues {
            let cff_energy_factor = match energy_factor_ownership {
                CffEnergyFactorOwnership::GlobalSourceProduct => inverse_energy_product.clone(),
                CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
            };
            crate::debug_tags!(#generation, #cff, #inspect;
                cut_index = ?cut_cff_index,
                log.cff_energy_factor = cff_energy_factor,
                "Exact CFF component-local energy-factor bridge"
            );
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: Vec::new(),
                exact_source_numerator: Some(exact_source_numerator.clone()),
            };
            for orientation in expr.orientations {
                let expression = orientation
                    .to_atom_gs()
                    .replace_multiple(&replacement_rules)
                    .replace_multiple(exact_source_numerator.mapper.exact_ose_replacements())
                    * &cff_energy_factor
                    * &cff_normalization;
                // Distinct exact factors remain in `expression`. The
                // orientation is deliberately source-local: its affine map
                // evaluates the complete parent numerator directly, rather
                // than being remapped to a production OrientationID.
                cff_term.orientations.push(CFFOrientationTerm {
                    expression,
                    orientation,
                    production_orientation_id: None,
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }
        Ok((
            CutCFF {
                terms,
                energy_degree_bound_report,
                production_prefactor_bridge,
            },
            contract_subgraph,
        ))
    }

    pub fn cff<S: SubGraphLike + SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<CutCFF> {
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }
        contract_edges.sort_unstable();
        contract_edges.dedup();

        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);

        // Reduced UV graphs use the production representation while retaining
        // their own affine numerator maps and independent residue sums.
        let generated = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            options,
            analysis_numerator,
        )?;
        self.cff_from_generated_expression(
            generated,
            contract_subgraph,
            cutset,
            orientation_pattern,
            options,
            analysis_numerator,
        )
    }

    pub(crate) fn cff_from_generated_expression<S: SubGraphLike + SubSetLike>(
        &self,
        generated: GeneratedThreeDExpression<esurface::Esurface, hsurface::Hsurface>,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<CutCFF> {
        let mut contract_edges = self
            .iter_edges_of(contract_subgraph)
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .collect::<Vec<_>>();
        contract_edges.sort_unstable();
        contract_edges.dedup();
        let energy_factor_ownership = generated.energy_factor_ownership;
        let source_energy_degree_bounds = generated.source_energy_degree_bounds.clone();
        let production_prefactor_bridge = CutCFF::gamma_loop_prefactor_conversion(&generated);
        let mut cff = generated.expression;
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut cff,
            &self.get_raised_edge_groups(),
        );

        let residues = select_indexed_cff_residues(
            cff,
            cutset,
            options.representation,
            &self.get_raised_edge_groups(),
            || {
                if options.energy_degree_bounds.is_some() || !source_energy_degree_bounds.is_empty()
                {
                    return Ok(source_energy_degree_bounds.clone());
                }
                let numerator = analysis_numerator.ok_or(
                three_dimensional_reps::generation::GenerationError::LtdResidueRequiresEnergyBounds,
            )?;
                Ok(self
                    .automatic_numerator_energy_degree_bounds_in_atoms_excluding_with_min_degree(
                        [numerator],
                        self.iter_edges_of(&self.initial_state_cut)
                            .chain(self.iter_edges_of(&self.tree_edges))
                            .map(|(_, edge, _)| edge),
                        1,
                    )?)
            },
        )?;

        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);
        // The CFF carries the measure normalization for the loop variables that remain after
        // contracting a UV subgraph. Fully contracted integrated CTs therefore get no extra CFF
        // measure factor, while ordinary root terms get the full graph-loop factor.
        let cff_loop_number = self
            .get_loop_number()
            .saturating_sub(self.cyclotomatic_number(contract_subgraph));
        let cff_phase = Atom::i().pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        let cff_energy_factor = match energy_factor_ownership {
            CffEnergyFactorOwnership::GlobalSourceProduct => {
                get_cff_inverse_energy_product_impl(self, &graph_without_is_cut, &contract_edges)
            }
            CffEnergyFactorOwnership::VariantLocal => Atom::num(1),
        };
        crate::debug_tags!(#cff, #trace;
            stage = "graph_cff_normalization",
            graph = %self.name,
            cff_loop_number = cff_loop_number,
            production_prefactor_bridge = ?production_prefactor_bridge,
            log.cff_normalization = cff_normalization,
            "Graph CFF normalization: graph={}, context={:?}, loops={}, bridge={}",
            self.name,
            options.cff_generation_context,
            cff_loop_number,
            production_prefactor_bridge.factor(),
        );

        let mut terms = BTreeMap::new();

        for (cut_cff_index, expr) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: vec![],
                exact_source_numerator: None,
            };
            for orientation in expr.orientations.iter().filter(|orientation| {
                orientation_pattern.filter_orientation(&orientation.data.orientation)
            }) {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                ose_expr *= &cff_energy_factor;
                ose_expr *= cff_normalization.clone();

                crate::debug_tags!(#cff, #trace;
                    stage = "graph_cff_term_expr",
                    graph = %self.name,
                    cut_index = ?cut_cff_index,
                    log.expr = ose_expr,
                    "Graph CFF term expression"
                );
                // println!("ose expr :{}", ose_expr);
                cff_term.orientations.push(CFFOrientationTerm {
                    expression: ose_expr,
                    orientation: orientation.clone(),
                    production_orientation_id: None,
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF {
            terms,
            energy_degree_bound_report: CffEnergyDegreeBoundReport {
                source_kind: CffEnergyBoundSourceKind::PhysicalGraph,
                physical_parent_bounds: source_energy_degree_bounds.clone(),
                assigned_cff_source_bounds: source_energy_degree_bounds,
            },
            production_prefactor_bridge,
        };
        Ok(cut_cff)
    }
}

#[cfg(test)]
mod tests;
