use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::{SuBitGraph, SubGraphLike, SubSetOps},
};
use spenso::structure::concrete_index::ExpandedIndex;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::Replacement,
};
use three_dimensional_reps::{
    EnergyEdgeIndexMap, Generate3DExpressionOptions, GraphOrientation, MomentumSignature,
    NumeratorSamplingScaleMode, OrientationSelector, ParsedGraph, RepresentationMode,
    graph_io::{
        ParsedGraphExternalEdge, ParsedGraphInitialStateCutEdge, ParsedGraphInternalEdge,
        initial_state_cut_external_alias,
    },
    repeated_groups,
};

use crate::{
    cff::{
        expression::{
            GammaLoopCFFVariant, GammaLoopGraphOrientation, GammaLoopOrientationExpression,
            OrientationID,
        },
        generation::generated_cff_expression_uses_variant_half_edges,
        surface::GammaLoopSurfaceCache,
    },
    graph::{FeynmanGraph, Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::{Sign, SignOrZero, sample::LoopIndex},
    numerator::energy_degree::EnergyPowerAnalyzer,
    settings::global::{GenerationSettings, ThreeDRepresentation},
    utils::{GS, symbolica_ext::LogPrint},
    uv::{UltravioletGraph, approx::ForestNodeLike},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Expanded4DApprox {
    pub sign: Sign,
    pub terms: Vec<Atom>,
}

#[derive(Clone, Debug)]
struct ExpandedDenominator {
    source_edge: EdgeIndex,
    mass_squared: Atom,
    full_expr: Atom,
}

#[derive(Clone, Debug)]
struct ExpandedTerm {
    numerator: Atom,
    denominators: Vec<ExpandedDenominator>,
}

#[derive(Clone, Debug)]
struct ParsedDenominatorEdge {
    local_edge_id: usize,
    source_edge: EdgeIndex,
    mass_squared: Atom,
    full_expr: Atom,
    is_preserved_4d: bool,
}

struct Expanded4DParsedSource {
    parsed: ParsedGraph,
    edge_map: EnergyEdgeIndexMap,
    denominator_edges: Vec<ParsedDenominatorEdge>,
}

impl Expanded4DApprox {
    pub fn expr(&self) -> (&[Atom], Sign) {
        (&self.terms, self.sign)
    }

    pub fn root(graph: &Graph) -> Self {
        let denominator_graph = graph.full_filter().subtract(&graph.initial_state_cut);
        Self {
            sign: Sign::Positive,
            terms: vec![Atom::num(1) / graph.denominator(&denominator_graph, |_| 1)],
        }
    }
}

pub fn expanded_4d_uv_start<S: ForestNodeLike>(
    graph: &Graph,
    current: &S,
    given: &S,
    integrand: &Atom,
) -> Result<Atom> {
    let reduced = current.reduced_subgraph(given);
    let numerator = graph
        .numerator(&reduced, given.subgraph())
        .to_d_dim(4)
        .get_single_atom()
        .map_err(|error| eyre!("graph numerator is not a single symbolic atom: {error}"))?;
    Ok(numerator * integrand)
}

pub fn expanded_4d_uv_kernel<S: ForestNodeLike>(
    graph: &Graph,
    current: &S,
    given: &S,
    integrand: &Atom,
) -> Result<Atom> {
    let reduced = current.reduced_subgraph(given);
    let n_loops = graph.n_loops(current.subgraph()) - graph.n_loops(given.subgraph());
    expanded_4d_uv_rescaled(graph, &reduced, n_loops, current.lmb(), integrand)
}

fn expanded_4d_uv_rescaled(
    graph: &Graph,
    replacement_subgraph: &SuBitGraph,
    n_loops: usize,
    lmb: &LoopMomentumBasis,
    atom: &Atom,
) -> Result<Atom> {
    let atomarg = graph.uv_rescaled(replacement_subgraph, n_loops, lmb, atom);
    let series = atomarg
        .series(GS.rescale, Atom::Zero, 0.into(), true)
        .map_err(|error| eyre!("expanded 4D local UV series expansion failed: {error}"))?;
    Ok(series.to_atom().replace(GS.rescale).with(Atom::num(1)))
}

#[allow(clippy::too_many_arguments)]
pub fn expanded_4d_terms_to_3d_parametric_integrands(
    graph: &mut Graph,
    atom: &Atom,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
    _valid_orientations: &[EdgeVec<Orientation>],
    contract_subgraph: &impl SubGraphLike<Base = linnet::half_edge::subgraph::SuBitGraph>,
) -> Result<Vec<Atom>> {
    let terms = extract_expanded_terms(atom)?;
    if terms.is_empty() {
        return Ok(vec![Atom::Zero]);
    }

    let mut out: Option<Vec<Atom>> = None;
    for term in terms {
        if let Some(generated) = project_expanded_4d_term_to_3d_parametric_integrands(
            graph,
            term,
            cutset,
            representation,
            settings,
            contract_subgraph,
        )? {
            accumulate_expanded_4d_term(&mut out, generated)?;
        }
    }

    Ok(out.unwrap_or_else(|| vec![Atom::Zero]))
}

#[allow(clippy::too_many_arguments)]
fn project_expanded_4d_term_to_3d_parametric_integrands(
    graph: &mut Graph,
    term: ExpandedTerm,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
    contract_subgraph: &impl SubGraphLike<Base = linnet::half_edge::subgraph::SuBitGraph>,
) -> Result<Option<Vec<Atom>>> {
    if term.denominators.is_empty() {
        if cutset_has_residue_selector(cutset) {
            return Ok(None);
        }
        return Ok(Some(vec![term.numerator]));
    }

    let source = build_expanded_4d_parsed_source(graph, &term.denominators, contract_subgraph)?;
    let expression = generate_expanded_4d_source_expression(
        graph,
        &source,
        &term.numerator,
        representation,
        settings,
    )?;
    let residues = select_cut_residues(expression, cutset)?;

    residues
        .into_iter()
        .map(|residue| {
            expanded_expression_parametric_atom(
                graph,
                &residue,
                &term.numerator,
                &source,
                cutset,
                representation,
                settings,
            )
        })
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

fn generate_expanded_4d_source_expression(
    graph: &mut Graph,
    source: &Expanded4DParsedSource,
    numerator: &Atom,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
) -> Result<crate::cff::expression::ThreeDExpression<OrientationID>> {
    let options =
        expression_options_for_expanded_term(source, numerator, representation, settings, graph)?;
    let raw_expression =
        three_dimensional_reps::generate_3d_expression_from_parsed(&source.parsed, &options)
            .map_err(|error| {
                eyre!(
                    "generalized 3D expression generation failed for expanded 4D local UV term: {error}\n{}",
                    expanded_source_summary(&source.parsed)
                )
            })?
            .remap_energy_edge_indices(&source.edge_map)
            .fuse_compatible_variants();

    let initial_state_cut_edges = graph
        .iter_edges_of(&graph.initial_state_cut)
        .map(|(_, edge_id, _)| edge_id)
        .collect_vec();
    graph.convert_generated_expression_surfaces(
        raw_expression,
        representation_mode(representation),
        generated_cff_expression_uses_variant_half_edges(&options),
        &graph.get_esurface_canonization(&graph.loop_momentum_basis),
        &initial_state_cut_edges,
    )
}

fn accumulate_expanded_4d_term(out: &mut Option<Vec<Atom>>, generated: Vec<Atom>) -> Result<()> {
    let target = out.get_or_insert_with(|| vec![Atom::Zero; generated.len()]);
    if target.len() != generated.len() {
        return Err(eyre!(
            "expanded 4D local UV term generated {} residue integrands, but previous terms generated {}",
            generated.len(),
            target.len()
        ));
    }
    for (sum, term) in target.iter_mut().zip(generated) {
        *sum += term;
    }
    Ok(())
}

fn cutset_has_residue_selector(cutset: &CutSet) -> bool {
    cutset.residue_selector.right_th_cut.is_some()
        || cutset.residue_selector.left_th_cut.is_some()
        || cutset.residue_selector.lu_cut.is_some()
}

fn select_cut_residues(
    expression: crate::cff::expression::ThreeDExpression<OrientationID>,
    cutset: &CutSet,
) -> Result<Vec<crate::cff::expression::ThreeDExpression<OrientationID>>> {
    let residue = if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        let residues = expression.select_esurface_residue(right_threshold);
        expect_single_threshold_residue(residues, "right")?
    } else {
        expression
    };

    let residue = if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        let residues = residue.select_esurface_residue(left_threshold);
        expect_single_threshold_residue(residues, "left")?
    } else {
        residue
    };

    Ok(
        if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
            residue.select_esurface_residue(lu_cut)
        } else {
            vec![residue]
        },
    )
}

fn expect_single_threshold_residue(
    mut residues: Vec<crate::cff::expression::ThreeDExpression<OrientationID>>,
    side: &str,
) -> Result<crate::cff::expression::ThreeDExpression<OrientationID>> {
    if residues.len() != 1 {
        return Err(eyre!(
            "{side} threshold residue produced {} expressions; expanded 4D local UV generation supports exactly one expression at this stage",
            residues.len()
        ));
    }
    residues
        .pop()
        .ok_or_else(|| eyre!("{side} threshold residue did not produce an expression"))
}

fn expanded_expression_parametric_atom(
    graph: &Graph,
    expression: &crate::cff::expression::ThreeDExpression<OrientationID>,
    numerator: &Atom,
    source: &Expanded4DParsedSource,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
) -> Result<Atom> {
    let replacement_rules = expression.surfaces.get_all_replacements_gs(&[]);
    let ose_replacements = expanded_ose_replacements(source);
    let residual_denominator_factor = residual_denominator_factor_for_expanded_source(source);
    let inverse_energy_product = if representation == ThreeDRepresentation::Cff
        && !crate::cff::cff_expression_uses_local_half_edges(expression)
    {
        cff_inverse_energy_product_for_expanded_source(source)
    } else {
        Atom::num(1)
    };
    let ltd_parity_correction =
        ltd_cross_section_residue_parity_correction(graph, source, cutset, representation);
    let all_orientations = crate::settings::global::OrientationPattern::default();
    let pattern = if settings.explicit_orientation_sum_only {
        &all_orientations
    } else {
        &settings.orientation_pattern
    };

    let mut sum = Atom::Zero;
    for orientation in expression.orientations.iter() {
        if !pattern.filter_orientation(orientation.orientation()) {
            continue;
        }
        let mut atom = numerator.replace_multiple(&orientation.energy_replacements_gs(graph));
        atom *= orientation
            .variants
            .iter()
            .map(GammaLoopCFFVariant::to_atom_gs)
            .reduce(|acc, term| acc + term)
            .unwrap_or_else(Atom::new);
        atom = atom.replace_multiple(&replacement_rules);
        if !settings.explicit_orientation_sum_only {
            atom *= orientation.orientation_thetas_gs();
        }
        sum += atom;
    }

    Ok(
        (sum * inverse_energy_product * residual_denominator_factor * ltd_parity_correction)
            .replace_multiple(&ose_replacements)
            .replace(GS.dim)
            .with(4)
            .simplify_color()
            .expand_dots()?
            .collect_factors(),
    )
}

fn ltd_cross_section_residue_parity_correction(
    graph: &Graph,
    source: &Expanded4DParsedSource,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
) -> Atom {
    if representation != ThreeDRepresentation::Ltd || cutset.residue_selector.lu_cut.is_none() {
        return Atom::num(1);
    }

    // Cross-section LU residues are converted to GammaLoop's CFF production
    // convention after the forest is assembled. Expanded-4D local UV terms are
    // generated on reduced parsed sources, whose generalized-CFF global sign
    // can differ from the full graph, especially for repeated denominators.
    // Correct each reduced source to the full-graph convention before the
    // outer cross-section parity factor is applied.
    let full_exponent = graph.cff_global_sign_exponent_for_3d_expression();
    let source_exponent = cff_global_sign_exponent_for_expanded_source(source);
    if (full_exponent + source_exponent).is_multiple_of(2) {
        Atom::num(1)
    } else {
        Atom::num(-1)
    }
}

fn cff_inverse_energy_product_for_expanded_source(source: &Expanded4DParsedSource) -> Atom {
    source
        .denominator_edges
        .iter()
        .filter(|edge| !edge.is_preserved_4d)
        .map(|edge| {
            Atom::num(1) / (-Atom::num(2) * crate::utils::ose_atom_from_index(edge.source_edge))
        })
        .reduce(|acc, factor| acc * factor)
        .unwrap_or_else(|| Atom::num(1))
}

fn cff_global_sign_exponent_for_expanded_source(source: &Expanded4DParsedSource) -> usize {
    let preserved_edges = source
        .denominator_edges
        .iter()
        .filter_map(|edge| edge.is_preserved_4d.then_some(edge.local_edge_id))
        .collect::<BTreeSet<_>>();
    let duplicate_signature_excess = repeated_groups(&source.parsed)
        .into_iter()
        .map(|group| {
            group
                .edge_ids
                .into_iter()
                .filter(|edge_id| !preserved_edges.contains(edge_id))
                .count()
                .saturating_sub(1)
        })
        .sum::<usize>();
    source.parsed.loop_names.len().saturating_sub(1) + duplicate_signature_excess
}

fn residual_denominator_factor_for_expanded_source(source: &Expanded4DParsedSource) -> Atom {
    source
        .denominator_edges
        .iter()
        .filter(|edge| edge.is_preserved_4d)
        .map(|edge| Atom::num(1) / &edge.full_expr)
        .reduce(|acc, factor| acc * factor)
        .unwrap_or_else(|| Atom::num(1))
}

fn expanded_ose_replacements(source: &Expanded4DParsedSource) -> Vec<Replacement> {
    source
        .denominator_edges
        .iter()
        .map(|edge| {
            let ose = crate::utils::ose_atom_from_index(edge.source_edge);
            Replacement::new(ose.to_pattern(), expanded_ose_atom(edge).to_pattern())
        })
        .collect()
}

fn expanded_ose_atom(edge: &ParsedDenominatorEdge) -> Atom {
    let mut q3q3 = Atom::Zero;
    for spatial_index in 1..=3 {
        let component = GS.emr_mom(
            edge.source_edge,
            Atom::from(ExpandedIndex::from_iter([spatial_index])),
        );
        q3q3 += component.clone() * component;
    }
    (edge.mass_squared.clone() + q3q3).sqrt()
}

fn expression_options_for_expanded_term(
    source: &Expanded4DParsedSource,
    numerator: &Atom,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
    graph: &Graph,
) -> Result<Generate3DExpressionOptions> {
    let local_to_source = source
        .denominator_edges
        .iter()
        .map(|edge| (edge.local_edge_id, edge.source_edge))
        .collect::<BTreeMap<_, _>>();
    let active_edges = local_to_source.values().copied().collect::<BTreeSet<_>>();
    let degree_by_source = EnergyPowerAnalyzer::with_internal_edges(
        graph.loop_momentum_basis.loop_edges.iter().copied(),
        active_edges,
    )
    .analyze_atom_upper_bound(numerator)?;

    let mut energy_degree_bounds = Vec::new();
    for edge in &source.denominator_edges {
        if edge.is_preserved_4d {
            continue;
        }
        let degree = degree_by_source
            .iter()
            .find_map(|(source_edge, degree)| (source_edge == edge.source_edge).then_some(degree))
            .unwrap_or(0);
        energy_degree_bounds.push((edge.local_edge_id, degree));
    }

    Ok(Generate3DExpressionOptions {
        representation: representation_mode(representation),
        energy_degree_bounds,
        numerator_sampling_scale: numerator_sampling_scale_mode(
            settings.uniform_numerator_sampling_scale,
        ),
        include_cff_duplicate_signature_excess_sign: true,
        preserve_internal_edges_as_four_d_denominators: source
            .denominator_edges
            .iter()
            .filter_map(|edge| edge.is_preserved_4d.then_some(edge.local_edge_id))
            .collect(),
    })
}

fn build_expanded_4d_parsed_source(
    graph: &Graph,
    denominators: &[ExpandedDenominator],
    contract_subgraph: &impl SubGraphLike<Base = linnet::half_edge::subgraph::SuBitGraph>,
) -> Result<Expanded4DParsedSource> {
    let mut parent = (0..graph.n_nodes()).collect_vec();
    for (pair, _edge_id, edge_data) in graph.underlying.iter_edges() {
        if edge_data.data.is_dummy || !contract_subgraph.includes(&pair) {
            continue;
        }
        if let HedgePair::Paired { source, sink } = pair {
            union_parent(
                &mut parent,
                usize::from(graph.node_id(source)),
                usize::from(graph.node_id(sink)),
            );
        }
    }

    let node_ids = graph
        .underlying
        .iter_nodes()
        .map(|(node_id, _, _)| node_id)
        .sorted()
        .collect::<Vec<_>>();
    let mut root_to_internal = BTreeMap::<usize, usize>::new();
    for node_id in &node_ids {
        let root = find_parent(&mut parent, usize::from(*node_id));
        if !root_to_internal.contains_key(&root) {
            root_to_internal.insert(root, root_to_internal.len());
        }
    }
    let mut node_to_internal = BTreeMap::new();
    for node_id in &node_ids {
        let root = find_parent(&mut parent, usize::from(*node_id));
        node_to_internal.insert(*node_id, root_to_internal[&root]);
    }

    let active_loop_columns = active_loop_columns_for_denominators(graph, denominators);
    let loop_names = active_loop_columns
        .iter()
        .map(|loop_index| {
            graph.underlying[graph.loop_momentum_basis.loop_edges[LoopIndex::from(*loop_index)]]
                .name
                .value
                .clone()
        })
        .collect::<Vec<_>>();
    let external_names = graph
        .loop_momentum_basis
        .ext_edges
        .iter()
        .map(|edge_id| graph.underlying[*edge_id].name.value.clone())
        .collect::<Vec<_>>();

    let mut internal_edges = Vec::new();
    let mut denominator_edges = Vec::new();
    let mut initial_state_cut_edges = Vec::new();
    let mut edge_map_internal = BTreeMap::new();
    let denominator_source_edges = denominators
        .iter()
        .map(|denominator| denominator.source_edge)
        .collect::<BTreeSet<_>>();

    for (_, edge_index, _) in graph
        .underlying
        .iter_edges_of(&graph.initial_state_cut)
        .sorted_by_key(|(_, edge_index, _)| *edge_index)
    {
        let edge_data = &graph.underlying[edge_index];
        if edge_data.is_dummy || denominator_source_edges.contains(&edge_index) {
            continue;
        }
        let (_, pair) = graph[&edge_index];
        let HedgePair::Paired { source, sink } = pair else {
            return Err(eyre!(
                "expanded 4D initial-state cut edge {edge_index} is not a paired internal edge",
            ));
        };
        let local_edge_id = internal_edges.len();
        let source_node = *node_to_internal
            .get(&graph.node_id(source))
            .ok_or_else(|| eyre!("missing source-node mapping for edge {edge_index}"))?;
        let sink_node = *node_to_internal
            .get(&graph.node_id(sink))
            .ok_or_else(|| eyre!("missing sink-node mapping for edge {edge_index}"))?;
        let signature = &graph.loop_momentum_basis.edge_signatures[edge_index];
        let momentum_signature = MomentumSignature {
            loop_signature: active_loop_columns
                .iter()
                .map(|loop_index| sign_to_i32(signature.internal[LoopIndex::from(*loop_index)]))
                .collect(),
            external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
        };
        let (external_id, external_sign) =
            initial_state_cut_external_alias(usize::from(edge_index), &momentum_signature)?;
        internal_edges.push(ParsedGraphInternalEdge {
            edge_id: local_edge_id,
            tail: source_node,
            head: sink_node,
            label: graph.underlying[edge_index].name.value.clone(),
            mass_key: Some(edge_data.particle.mass_atom().to_canonical_string()),
            signature: momentum_signature,
            had_pow: false,
        });
        initial_state_cut_edges.push(ParsedGraphInitialStateCutEdge {
            edge_id: local_edge_id,
            external_id,
            external_sign,
        });
        edge_map_internal.insert(local_edge_id, usize::from(edge_index));
    }

    for denominator in denominators {
        let (_, pair) = graph[&denominator.source_edge];
        let HedgePair::Paired { source, sink } = pair else {
            return Err(eyre!(
                "expanded 4D denominator edge {} is not a paired internal edge",
                denominator.source_edge
            ));
        };
        let local_edge_id = internal_edges.len();
        let source_node = *node_to_internal
            .get(&graph.node_id(source))
            .ok_or_else(|| {
                eyre!(
                    "missing source-node mapping for edge {}",
                    denominator.source_edge
                )
            })?;
        let sink_node = *node_to_internal.get(&graph.node_id(sink)).ok_or_else(|| {
            eyre!(
                "missing sink-node mapping for edge {}",
                denominator.source_edge
            )
        })?;
        let signature = &graph.loop_momentum_basis.edge_signatures[denominator.source_edge];
        let momentum_signature = MomentumSignature {
            loop_signature: active_loop_columns
                .iter()
                .map(|loop_index| sign_to_i32(signature.internal[LoopIndex::from(*loop_index)]))
                .collect(),
            external_signature: (&signature.external).into_iter().map(sign_to_i32).collect(),
        };
        let is_preserved_4d = momentum_signature
            .loop_signature
            .iter()
            .all(|sign| *sign == 0);
        internal_edges.push(ParsedGraphInternalEdge {
            edge_id: local_edge_id,
            tail: source_node,
            head: sink_node,
            label: graph.underlying[denominator.source_edge].name.value.clone(),
            mass_key: Some(denominator.mass_squared.to_canonical_string()),
            signature: momentum_signature,
            had_pow: false,
        });
        denominator_edges.push(ParsedDenominatorEdge {
            local_edge_id,
            source_edge: denominator.source_edge,
            mass_squared: denominator.mass_squared.clone(),
            full_expr: denominator.full_expr.clone(),
            is_preserved_4d,
        });
        edge_map_internal.insert(local_edge_id, usize::from(denominator.source_edge));
    }

    let external_edges = graph
        .underlying
        .iter_edges()
        .filter_map(|(pair, edge_id, _)| {
            let HedgePair::Unpaired { hedge, flow } = pair else {
                return None;
            };
            let node = *node_to_internal.get(&graph.node_id(hedge))?;
            let signature = &graph.loop_momentum_basis.edge_signatures[edge_id];
            let (source, destination) = match flow {
                linnet::half_edge::involution::Flow::Source => (Some(node), None),
                linnet::half_edge::involution::Flow::Sink => (None, Some(node)),
            };
            Some(ParsedGraphExternalEdge {
                edge_id: 10_000_000usize + usize::from(edge_id),
                source,
                destination,
                label: graph.underlying[edge_id].name.value.clone(),
                external_coefficients: (&signature.external).into_iter().map(sign_to_i32).collect(),
            })
        })
        .collect::<Vec<_>>();

    let parsed = ParsedGraph {
        internal_edges,
        external_edges,
        initial_state_cut_edges,
        loop_names,
        external_names,
        node_name_to_internal: root_to_internal
            .into_iter()
            .map(|(root, node)| (format!("n{root}"), node))
            .collect(),
    };

    let edge_map = EnergyEdgeIndexMap {
        internal: edge_map_internal,
        external: graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .enumerate()
            .map(|(external_id, edge_id)| (external_id, usize::from(*edge_id)))
            .collect(),
        orientation_edge_count: graph.underlying.n_edges(),
    };

    Ok(Expanded4DParsedSource {
        parsed,
        edge_map,
        denominator_edges,
    })
}

fn active_loop_columns_for_denominators(
    graph: &Graph,
    denominators: &[ExpandedDenominator],
) -> Vec<usize> {
    let rows = denominators
        .iter()
        .map(|denominator| {
            graph.loop_momentum_basis.edge_signatures[denominator.source_edge]
                .internal
                .iter()
                .map(|sign| sign_to_i32(*sign))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    independent_loop_columns(&rows)
}

fn extract_expanded_terms(atom: &Atom) -> Result<Vec<ExpandedTerm>> {
    let expanded = atom.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect::<Vec<_>>(),
        _ => vec![expanded],
    };

    terms
        .into_iter()
        .map(|term| {
            let mut denominators = Vec::new();
            let numerator = collect_denominators_from_view(term.as_view(), &mut denominators)?;
            Ok(ExpandedTerm {
                numerator: numerator.collect_factors(),
                denominators,
            })
        })
        .collect()
}

fn collect_denominators_from_view(
    view: AtomView<'_>,
    denominators: &mut Vec<ExpandedDenominator>,
) -> Result<Atom> {
    match view {
        AtomView::Mul(mul) => {
            let mut numerator = Atom::num(1);
            for factor in mul.iter() {
                numerator *= collect_denominators_from_view(factor, denominators)?;
            }
            Ok(numerator)
        }
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            if let Some(denominator) = denominator_from_view(base)? {
                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(eyre!(
                        "expanded 4D local UV denominator has non-integer power `{}`",
                        exponent.to_owned().log_print(None)
                    ));
                };
                if exponent >= 0 {
                    return Ok(view.to_owned());
                }
                for _ in 0..exponent.unsigned_abs() {
                    denominators.push(denominator.clone());
                }
                Ok(Atom::num(1))
            } else {
                Ok(view.to_owned())
            }
        }
        _ => Ok(view.to_owned()),
    }
}

fn denominator_from_view(view: AtomView<'_>) -> Result<Option<ExpandedDenominator>> {
    let AtomView::Fun(function) = view else {
        return Ok(None);
    };
    if function.get_symbol() != GS.den {
        return Ok(None);
    }
    if function.get_nargs() != 4 {
        return Err(eyre!(
            "expected expanded 4D denominator wrapper with four arguments, found {}",
            function.get_nargs()
        ));
    }
    let source_edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
        eyre!(
            "expanded 4D denominator wrapper has non-integer edge id `{}`",
            function.get(0).to_owned().log_print(None)
        )
    })?);
    Ok(Some(ExpandedDenominator {
        source_edge,
        mass_squared: function.get(2).to_owned(),
        full_expr: function.get(3).to_owned(),
    }))
}

fn numerator_sampling_scale_mode(
    setting: crate::settings::global::UniformNumeratorSamplingScale,
) -> NumeratorSamplingScaleMode {
    match setting {
        crate::settings::global::UniformNumeratorSamplingScale::None => {
            NumeratorSamplingScaleMode::None
        }
        crate::settings::global::UniformNumeratorSamplingScale::BeyondQuadratic => {
            NumeratorSamplingScaleMode::BeyondQuadratic
        }
        crate::settings::global::UniformNumeratorSamplingScale::All => {
            NumeratorSamplingScaleMode::All
        }
    }
}

fn representation_mode(representation: ThreeDRepresentation) -> RepresentationMode {
    match representation {
        ThreeDRepresentation::Cff => RepresentationMode::Cff,
        ThreeDRepresentation::Ltd => RepresentationMode::Ltd,
    }
}

fn expanded_source_summary(parsed: &ParsedGraph) -> String {
    let internal_edges = parsed
        .internal_edges
        .iter()
        .map(|edge| {
            format!(
                "edge {}: {} -> {}, loops={:?}, externals={:?}, mass={:?}",
                edge.edge_id,
                edge.tail,
                edge.head,
                edge.signature.loop_signature,
                edge.signature.external_signature,
                edge.mass_key
            )
        })
        .join("; ");
    format!(
        "expanded 4D UV 3D source has {} loop names {:?}, {} internal edges [{}]",
        parsed.loop_names.len(),
        parsed.loop_names,
        parsed.internal_edges.len(),
        internal_edges
    )
}

fn find_parent(parent: &mut [usize], node: usize) -> usize {
    let parent_node = parent[node];
    if parent_node == node {
        node
    } else {
        let root = find_parent(parent, parent_node);
        parent[node] = root;
        root
    }
}

fn union_parent(parent: &mut [usize], left: usize, right: usize) {
    let left_root = find_parent(parent, left);
    let right_root = find_parent(parent, right);
    if left_root != right_root {
        parent[right_root] = left_root;
    }
}

fn independent_loop_columns(rows: &[Vec<i32>]) -> Vec<usize> {
    let Some(column_count) = rows.iter().map(Vec::len).max() else {
        return Vec::new();
    };
    let mut selected = Vec::new();
    let mut rank = 0;
    for column in 0..column_count {
        let mut candidate = selected.clone();
        candidate.push(column);
        let candidate_rank = integer_matrix_column_rank(rows, &candidate);
        if candidate_rank > rank {
            selected.push(column);
            rank = candidate_rank;
        }
    }
    selected
}

fn integer_matrix_column_rank(rows: &[Vec<i32>], columns: &[usize]) -> usize {
    let mut matrix = rows
        .iter()
        .map(|row| {
            columns
                .iter()
                .map(|column| {
                    symbolica::domains::rational::Rational::from(
                        row.get(*column).copied().unwrap_or_default(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .filter(|row| row.iter().any(|entry| *entry != 0))
        .collect::<Vec<_>>();

    let mut rank = 0;
    let mut pivot_column = 0;
    while rank < matrix.len() && pivot_column < columns.len() {
        let Some(pivot_row) = (rank..matrix.len()).find(|row| matrix[*row][pivot_column] != 0)
        else {
            pivot_column += 1;
            continue;
        };
        matrix.swap(rank, pivot_row);
        let pivot = matrix[rank][pivot_column].clone();
        for entry in matrix[rank].iter_mut().skip(pivot_column) {
            *entry /= pivot.clone();
        }
        let pivot_row = matrix[rank].clone();
        for (row, row_values) in matrix.iter_mut().enumerate() {
            if row == rank {
                continue;
            }
            let factor = row_values[pivot_column].clone();
            if factor == 0 {
                continue;
            }
            for (entry, pivot_entry) in row_values
                .iter_mut()
                .zip(pivot_row.iter())
                .skip(pivot_column)
            {
                *entry -= factor.clone() * pivot_entry.clone();
            }
        }
        rank += 1;
        pivot_column += 1;
    }
    rank
}

fn sign_to_i32(sign: SignOrZero) -> i32 {
    match sign {
        SignOrZero::Minus => -1,
        SignOrZero::Zero => 0,
        SignOrZero::Plus => 1,
    }
}
