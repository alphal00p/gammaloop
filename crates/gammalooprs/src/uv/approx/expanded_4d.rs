use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair, Orientation},
    subgraph::{Inclusion, SuBitGraph, SubSetLike, SubSetOps},
};
use spenso::structure::concrete_index::ExpandedIndex;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
};
use three_dimensional_reps::{
    EnergyEdgeIndexMap, Generate3DExpressionOptions, MomentumSignature, NumeratorSamplingScaleMode,
    ParsedGraph, RepresentationMode,
    graph_io::{ParsedGraphExternalEdge, ParsedGraphInternalEdge},
    repeated_groups,
};

use crate::{
    cff::{
        esurface::{Esurface, EsurfaceID, RaisedEsurfaceGroup},
        expression::{
            GammaLoopCFFVariant, GammaLoopOrientationExpression, OrientationID,
            RaisedEsurfaceGroupView, ThreeDExpression, localize_three_d_expression_on_esurface,
            normalize_cut_edge_support_with_raised_edge_groups,
            normalize_three_d_expression_cut_support_with_raised_edge_groups,
            remove_ltd_global_contact_completions_from_local_residue,
            select_lu_cut_residue_for_representation,
        },
        generation::generated_cff_expression_uses_variant_half_edges,
        surface::GammaLoopSurfaceCache,
    },
    graph::{FeynmanGraph, Graph, LoopMomentumBasis, cuts::CutSet},
    momentum::{Sign, SignOrZero, sample::LoopIndex},
    numerator::energy_degree::EnergyPowerAnalyzer,
    settings::global::{GenerationSettings, ThreeDRepresentation},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::{UltravioletGraph, approx::ForestNodeLike},
};
use typed_index_collections::TiVec;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Expanded4DApprox {
    pub sign: Sign,
    pub terms: Vec<Atom>,
}

#[derive(Clone, Debug)]
struct ExpandedDenominator {
    source_edge: EdgeIndex,
    momentum: Atom,
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
    momentum: Atom,
    mass_squared: Atom,
    full_expr: Atom,
    is_preserved_4d: bool,
}

struct Expanded4DParsedSource {
    parsed: ParsedGraph,
    edge_map: EnergyEdgeIndexMap,
    denominator_edges: Vec<ParsedDenominatorEdge>,
}

#[derive(Clone, Copy)]
pub struct Expanded4DSourceContext<'a> {
    cograph_contract_subgraph: &'a SuBitGraph,
    uv_leading_subgraph: Option<&'a SuBitGraph>,
}

impl<'a> Expanded4DSourceContext<'a> {
    pub fn cograph_only(cograph_contract_subgraph: &'a SuBitGraph) -> Self {
        Self {
            cograph_contract_subgraph,
            uv_leading_subgraph: None,
        }
    }

    pub fn local_uv(
        cograph_contract_subgraph: &'a SuBitGraph,
        uv_leading_subgraph: &'a SuBitGraph,
    ) -> Self {
        Self {
            cograph_contract_subgraph,
            uv_leading_subgraph: Some(uv_leading_subgraph),
        }
    }

    fn label(self) -> String {
        match self.uv_leading_subgraph {
            Some(uv_subgraph) => format!(
                "cograph_{}_uv_{}",
                self.cograph_contract_subgraph.string_label(),
                uv_subgraph.string_label()
            ),
            None => self.cograph_contract_subgraph.string_label(),
        }
    }
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
    if let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_EXPANDED_4D_RESCALING") {
        std::fs::create_dir_all(&dump_dir).ok();
        let path = std::path::Path::new(&dump_dir)
            .join(format!("{}_expanded_4d_rescaling.txt", graph.name));
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .ok();
        if let Some(file) = file.as_mut() {
            use std::io::Write as _;
            writeln!(
                file,
                "replacement_subgraph={} n_loops={} input={} rescaled={}",
                replacement_subgraph.string_label(),
                n_loops,
                atom.to_canonical_string(),
                atomarg.to_canonical_string()
            )
            .ok();
        }
    }
    let series = atomarg
        .series(GS.rescale, Atom::Zero, 0)
        .map_err(|error| eyre!("expanded 4D local UV series expansion failed: {error}"))?;
    let result = series.to_atom().replace(GS.rescale).with(Atom::num(1));
    if let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_EXPANDED_4D_RESCALING") {
        let path = std::path::Path::new(&dump_dir)
            .join(format!("{}_expanded_4d_rescaling.txt", graph.name));
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .ok();
        if let Some(file) = file.as_mut() {
            use std::io::Write as _;
            writeln!(file, "series_result={}", result.to_canonical_string()).ok();
        }
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub fn expanded_4d_terms_to_3d_parametric_integrands(
    graph: &mut Graph,
    atom: &Atom,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
    valid_orientations: &[EdgeVec<Orientation>],
    source_context: Expanded4DSourceContext<'_>,
    root_expression: Option<&ThreeDExpression<OrientationID>>,
) -> Result<Vec<Atom>> {
    let terms = extract_expanded_terms(atom)?;
    if let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_EXPANDED_4D_PROJECTION") {
        std::fs::create_dir_all(&dump_dir).ok();
        let path = std::path::Path::new(&dump_dir).join(format!(
            "{}_{:?}_expanded_4d_projection.txt",
            graph.name, representation
        ));
        let mut dump = String::new();
        use std::fmt::Write;
        writeln!(
            dump,
            "terms={} has_selector={} contract_subgraph={}",
            terms.len(),
            cutset_has_residue_selector(cutset),
            source_context.label()
        )
        .ok();
        for (term_id, term) in terms.iter().enumerate() {
            writeln!(
                dump,
                "term {term_id}: numerator_zero={} denominators={}",
                term.numerator.is_zero(),
                term.denominators
                    .iter()
                    .map(|denominator| denominator.source_edge.0.to_string())
                    .join(",")
            )
            .ok();
        }
        let mut file = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .ok();
        if let Some(file) = file.as_mut() {
            use std::io::Write as _;
            writeln!(file, "{dump}").ok();
        }
    }
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
            valid_orientations,
            source_context,
            root_expression,
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
    valid_orientations: &[EdgeVec<Orientation>],
    source_context: Expanded4DSourceContext<'_>,
    root_expression: Option<&ThreeDExpression<OrientationID>>,
) -> Result<Option<Vec<Atom>>> {
    if term.denominators.is_empty() {
        if cutset_has_residue_selector(cutset) {
            return Ok(None);
        }
        return Ok(Some(vec![term.numerator]));
    }

    let source = build_expanded_4d_parsed_source(graph, &term.denominators, source_context)?;
    dump_expanded_4d_parsed_source_if_requested(
        graph,
        &source,
        &term.numerator,
        cutset,
        representation,
        settings,
        source_context,
    )?;
    let mut expression = generate_expanded_4d_source_expression(
        graph,
        &source,
        &term.numerator,
        representation,
        settings,
    )?;
    let raised_edge_groups = graph.get_raised_edge_groups();
    normalize_three_d_expression_cut_support_with_raised_edge_groups(
        &mut expression,
        &raised_edge_groups,
    );

    let source_cutset = remap_cutset_to_expression_surfaces(
        graph,
        cutset,
        root_expression,
        &expression.surfaces.esurface_cache,
        &raised_edge_groups,
        source_context,
        representation,
    )?;

    let residues = select_cut_residues(expression, &source_cutset, representation)?;
    dump_expanded_4d_selected_residues_if_requested(
        graph,
        &source,
        &residues,
        representation,
        source_context,
    );

    residues
        .into_iter()
        .map(|residue| {
            expanded_expression_parametric_atom(
                graph,
                &residue,
                &term.numerator,
                &source,
                &source_cutset,
                representation,
                settings,
                valid_orientations,
            )
        })
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

fn dump_expanded_4d_selected_residues_if_requested(
    graph: &Graph,
    source: &Expanded4DParsedSource,
    residues: &[crate::cff::expression::ThreeDExpression<OrientationID>],
    representation: ThreeDRepresentation,
    source_context: Expanded4DSourceContext<'_>,
) {
    let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_EXPANDED_4D_SELECTED_RESIDUES") else {
        return;
    };
    std::fs::create_dir_all(&dump_dir).ok();
    let path = std::path::Path::new(&dump_dir).join(format!(
        "{}_{:?}_{}_selected_residues.txt",
        graph.name,
        representation,
        source_context.label()
    ));
    let mut dump = String::new();
    use std::fmt::Write;
    writeln!(dump, "graph {}", graph.name).ok();
    writeln!(dump, "representation {:?}", representation).ok();
    writeln!(dump, "contract_subgraph {}", source_context.label()).ok();
    writeln!(dump, "source {}", expanded_source_summary(&source.parsed)).ok();
    writeln!(dump, "residue_count {}", residues.len()).ok();
    for (residue_id, residue) in residues.iter().enumerate() {
        writeln!(
            dump,
            "residue {residue_id}: orientations={} esurfaces={} hsurfaces={} linear_surfaces={}",
            residue.orientations.len(),
            residue.surfaces.esurface_cache.len(),
            residue.surfaces.hsurface_cache.len(),
            residue.surfaces.linear_surface_cache.len()
        )
        .ok();
        writeln!(dump, "  esurfaces {:?}", residue.surfaces.esurface_cache).ok();
        writeln!(dump, "  hsurfaces {:?}", residue.surfaces.hsurface_cache).ok();
        for (orientation_id, orientation) in residue.orientations.iter_enumerated() {
            writeln!(
                dump,
                "  orientation {} label {:?} variants={} unfolded={}",
                usize::from(orientation_id),
                orientation.data.label,
                orientation.variants.len(),
                orientation.num_unfolded_terms()
            )
            .ok();
            for (variant_id, variant) in orientation.variants.iter().enumerate() {
                writeln!(
                    dump,
                    "    variant {variant_id} origin {:?} prefactor {} half_edges {:?} denominator_edges {:?} denominator_edge_support_signs {:?} numerator_surfaces {:?} denominator_nodes {}",
                    variant.origin,
                    variant.prefactor.to_canonical_string(),
                    variant.half_edges,
                    variant.denominator_edges,
                    variant.denominator_edge_support_signs,
                    variant.numerator_surfaces,
                    variant.denominator.get_num_nodes()
                )
                .ok();
                writeln!(dump, "      denominator {:?}", variant.denominator).ok();
            }
        }
    }
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .ok();
    if let Some(file) = file.as_mut() {
        use std::io::Write as _;
        writeln!(file, "{dump}").ok();
    }
}

fn dump_expanded_4d_parsed_source_if_requested(
    graph: &Graph,
    source: &Expanded4DParsedSource,
    numerator: &Atom,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
    settings: &GenerationSettings,
    source_context: Expanded4DSourceContext<'_>,
) -> Result<()> {
    let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_EXPANDED_4D_SOURCE") else {
        return Ok(());
    };
    let options =
        expression_options_for_expanded_term(source, numerator, representation, settings, graph)?;
    std::fs::create_dir_all(&dump_dir).ok();
    let path = std::path::Path::new(&dump_dir).join(format!(
        "{}_{:?}_{}_source.txt",
        graph.name,
        representation,
        source_context.label()
    ));
    let mut dump = String::new();
    use std::fmt::Write;
    writeln!(dump, "graph {}", graph.name).ok();
    writeln!(dump, "representation {:?}", representation).ok();
    writeln!(dump, "contract_subgraph {}", source_context.label()).ok();
    writeln!(dump, "has_selector {}", cutset_has_residue_selector(cutset)).ok();
    writeln!(dump, "lu_cut {:?}", cutset.residue_selector.lu_cut).ok();
    if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
        write_raised_group_surfaces(&mut dump, graph, "lu_cut", lu_cut);
    }
    writeln!(
        dump,
        "lu_cut_edge_sets {:?}",
        cutset.residue_selector.lu_cut_edge_sets
    )
    .ok();
    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        write_raised_group_surfaces(&mut dump, graph, "left_th_cut", left_threshold);
    }
    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        write_raised_group_surfaces(&mut dump, graph, "right_th_cut", right_threshold);
    }
    writeln!(
        dump,
        "energy_degree_bounds {:?}",
        options.energy_degree_bounds
    )
    .ok();
    writeln!(
        dump,
        "preserve_internal_edges_as_four_d_denominators {:?}",
        options.preserve_internal_edges_as_four_d_denominators
    )
    .ok();
    writeln!(dump, "numerator {}", numerator.to_canonical_string()).ok();
    writeln!(dump, "source {}", expanded_source_summary(&source.parsed)).ok();
    writeln!(dump, "denominator_edges").ok();
    for edge in &source.denominator_edges {
        writeln!(
            dump,
            "  local={} source={} preserved={} mass={} momentum={} full={}",
            edge.local_edge_id,
            edge.source_edge.0,
            edge.is_preserved_4d,
            edge.mass_squared.to_canonical_string(),
            edge.momentum.to_canonical_string(),
            edge.full_expr.to_canonical_string()
        )
        .ok();
    }
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .ok();
    if let Some(file) = file.as_mut() {
        use std::io::Write as _;
        writeln!(file, "{dump}").ok();
    }
    Ok(())
}

fn write_raised_group_surfaces(
    dump: &mut String,
    graph: &Graph,
    label: &str,
    group: &impl RaisedEsurfaceGroupView,
) {
    use std::fmt::Write;
    for esurface_id in group.esurface_ids() {
        if let Some(esurface) = graph.surface_cache.esurface_cache.get(*esurface_id) {
            writeln!(
                dump,
                "{label} surface {:?}: energies={:?} external_shift={:?}",
                esurface_id, esurface.energies, esurface.external_shift
            )
            .ok();
        }
    }
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
    let use_generated_cff_half_edges =
        generated_cff_expression_uses_variant_half_edges(&options, &raw_expression);
    graph.convert_generated_expression_surfaces(
        raw_expression,
        representation_mode(representation),
        use_generated_cff_half_edges,
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

fn remap_cutset_to_expression_surfaces(
    graph: &Graph,
    cutset: &CutSet,
    root_expression: Option<&ThreeDExpression<OrientationID>>,
    generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
    raised_edge_groups: &[Vec<EdgeIndex>],
    source_context: Expanded4DSourceContext<'_>,
    representation: ThreeDRepresentation,
) -> Result<CutSet> {
    if !cutset_has_residue_selector(cutset) {
        return Ok(cutset.clone());
    }

    let root_expression = root_expression.ok_or_else(|| {
        eyre!(
            "expanded 4D local UV residue selection for graph {} requires the production root 3D expression surface cache",
            graph.name
        )
    })?;
    let normalized_generated_esurfaces = generated_esurfaces
        .iter()
        .map(|esurface| {
            Graph::normalize_esurface_with_raised_edge_groups(esurface, &raised_edge_groups)
        })
        .collect::<TiVec<EsurfaceID, _>>();
    let reference_esurfaces = &root_expression.surfaces.esurface_cache;

    let remap_group = |group: &RaisedEsurfaceGroup, label: &str| -> Result<RaisedEsurfaceGroup> {
        let esurface_ids = group
                .esurface_ids
                .iter()
                .map(|reference_id| {
                    let reference = reference_esurfaces.get(*reference_id).ok_or_else(|| {
                        eyre!(
                            "{label} residue selector in graph {} references E-surface {:?}, but the production root expression contains only {} E-surfaces",
                            graph.name,
                            reference_id,
                            reference_esurfaces.len()
                        )
                    })?;
                    find_matching_generated_esurface(
                        graph,
                        reference,
                        generated_esurfaces,
                        &normalized_generated_esurfaces,
                        &raised_edge_groups,
                        format!(
                            "{label} residue selector for {:?} expanded-4D source {}",
                            representation,
                            source_context.label()
                        ),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
        Ok(RaisedEsurfaceGroup {
            esurface_ids,
            max_occurence: group.max_occurence,
        })
    };

    let mut remapped = cutset.clone();
    let remapped_lu_cut = cutset
        .residue_selector
        .lu_cut
        .as_ref()
        .map(|group| remap_group(group, "Cutkosky"))
        .transpose()?;
    remapped.residue_selector.ltd_lu_cut_esurface_signs = if let (
        Some(original_lu_cut),
        Some(remapped_lu_cut),
    ) = (
        cutset.residue_selector.lu_cut.as_ref(),
        remapped_lu_cut.as_ref(),
    ) {
        cutset
            .residue_selector
            .ltd_lu_cut_esurface_signs
            .iter()
            .map(|(esurface_id, sign)| {
                let position = original_lu_cut
                    .esurface_ids
                    .iter()
                    .position(|original_id| original_id == esurface_id)
                    .ok_or_else(|| {
                        eyre!(
                            "Cutkosky LTD orientation sign references E-surface {esurface_id:?}, but the original selector group is {:?}",
                            original_lu_cut.esurface_ids
                        )
                    })?;
                Ok((remapped_lu_cut.esurface_ids[position], *sign))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };
    remapped.residue_selector.lu_cut = remapped_lu_cut;
    remapped.residue_selector.lu_cut_edge_sets = cutset
        .residue_selector
        .lu_cut_edge_sets
        .iter()
        .map(|cut_edges| {
            normalize_cut_edge_support_with_raised_edge_groups(cut_edges, raised_edge_groups)
        })
        .collect();
    remapped.residue_selector.left_th_cut = cutset
        .residue_selector
        .left_th_cut
        .as_ref()
        .map(|group| remap_group(group, "left threshold"))
        .transpose()?;
    remapped.residue_selector.right_th_cut = cutset
        .residue_selector
        .right_th_cut
        .as_ref()
        .map(|group| remap_group(group, "right threshold"))
        .transpose()?;
    Ok(remapped)
}

fn find_matching_generated_esurface(
    graph: &Graph,
    requested: &Esurface,
    generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
    normalized_generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
    raised_edge_groups: &[Vec<EdgeIndex>],
    context: impl std::fmt::Display,
) -> Result<EsurfaceID> {
    let normalized_requested =
        Graph::normalize_esurface_with_raised_edge_groups(requested, raised_edge_groups);
    if let Some(position) = generated_esurfaces
        .iter()
        .position(|generated| generated == requested)
    {
        return Ok(position.into());
    }
    generated_esurfaces
        .iter()
        .zip(normalized_generated_esurfaces.iter())
        .position(|(_, normalized_generated)| normalized_generated == &normalized_requested)
        .ok_or_else(|| {
            eyre!(
                "{context} E-surface {requested:?} for graph {} was not present in the expanded 4D generated source expression.\n\
                 Normalized requested E-surface: {:?}\n\
                 Generated E-surfaces: {:?}",
                graph.name,
                normalized_requested,
                generated_esurfaces,
            )
        })
        .map(Into::into)
}

fn select_cut_residues(
    expression: crate::cff::expression::ThreeDExpression<OrientationID>,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
) -> Result<Vec<crate::cff::expression::ThreeDExpression<OrientationID>>> {
    let mut residues = vec![expression];

    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        residues = residues
            .into_iter()
            .map(|residue| {
                let residues = residue.select_esurface_residue(right_threshold);
                expect_single_threshold_residue(residues, "right")
            })
            .collect::<Result<Vec<_>>>()?;
    }

    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        residues = residues
            .into_iter()
            .map(|residue| {
                let residues = residue.select_esurface_residue(left_threshold);
                expect_single_threshold_residue(residues, "left")
            })
            .collect::<Result<Vec<_>>>()?;
    };

    if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
        residues = residues
            .into_iter()
            .flat_map(|expression| {
                select_lu_cut_residue_for_representation(
                    expression,
                    lu_cut,
                    &cutset.residue_selector.lu_cut_edge_sets,
                    &cutset.residue_selector.ltd_lu_cut_esurface_signs,
                    representation,
                )
            })
            .collect();
    }

    if representation == ThreeDRepresentation::Ltd {
        for residue in &mut residues {
            remove_ltd_global_contact_completions_from_local_residue(residue);
            localize_ltd_threshold_residue_if_needed(residue, cutset)?;
        }
    }

    Ok(residues)
}

fn localize_ltd_threshold_residue_if_needed(
    residue: &mut crate::cff::expression::ThreeDExpression<OrientationID>,
    cutset: &CutSet,
) -> Result<()> {
    if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
        for esurface_id in right_threshold.esurface_ids() {
            localize_three_d_expression_on_esurface(residue, *esurface_id)?;
        }
    }
    if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
        for esurface_id in left_threshold.esurface_ids() {
            localize_three_d_expression_on_esurface(residue, *esurface_id)?;
        }
    }
    Ok(())
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
    valid_orientations: &[EdgeVec<Orientation>],
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
    let source_global_sign_factor =
        cross_section_residue_source_global_sign_factor(graph, source, cutset, representation);
    let mut sum = Atom::Zero;
    for orientation in expression.orientations.iter() {
        let mut atom = numerator.replace_multiple(orientation.energy_replacements_gs(graph));
        atom *= orientation
            .variants
            .iter()
            .map(GammaLoopCFFVariant::to_atom_gs)
            .reduce(|acc, term| acc + term)
            .unwrap_or_else(Atom::new);
        atom = atom.replace_multiple(&replacement_rules);
        sum += atom;
    }
    if !settings.explicit_orientation_sum_only {
        if valid_orientations.is_empty() {
            return Err(eyre!(
                "expanded 4D local UV term for graph {} cannot be embedded in an empty production orientation space",
                graph.name
            ));
        }
        // The orientations generated for an expanded-4D UV source are the
        // internal energy-residue decomposition of that reduced source. They
        // are not the production graph's orientation variables. Sum the source
        // orientations explicitly, then distribute the result uniformly over
        // the outer orientation projector used by ordinary CFF evaluation.
        sum /= Atom::num(valid_orientations.len() as i64);
    }

    Ok(
        (sum * inverse_energy_product * residual_denominator_factor * source_global_sign_factor)
            .replace_multiple(&ose_replacements)
            .replace(GS.dim)
            .with(4)
            .simplify_color()
            .expand_dots()?
            .collect_factors(),
    )
}

fn cross_section_residue_source_global_sign_factor(
    graph: &Graph,
    source: &Expanded4DParsedSource,
    cutset: &CutSet,
    representation: ThreeDRepresentation,
) -> Atom {
    if cutset.residue_selector.lu_cut.is_none() {
        return Atom::num(1);
    }
    if representation == ThreeDRepresentation::Ltd {
        return Atom::num(cutset.residue_selector.ltd_residue_prefactor_sign());
    }

    // Cross-section LU residues are assembled in GammaLoop's full-graph 3D
    // production convention. Expanded-4D local UV terms are projected through
    // reduced parsed sources, whose generalized 3D-rep global sign convention
    // is set by their loop count and repeated-signature excess. Convert each
    // reduced source to the full-graph convention before combining forest
    // terms; this is an algebraic convention change, not a representation-
    // dependent tuning of individual cuts.
    let full_exponent = graph.three_d_global_sign_exponent();
    let source_exponent = three_d_global_sign_exponent_for_expanded_source(source);
    let convention_bridge = if (full_exponent + source_exponent).is_multiple_of(2) {
        1
    } else {
        -1
    };
    Atom::num(convention_bridge)
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

fn three_d_global_sign_exponent_for_expanded_source(source: &Expanded4DParsedSource) -> usize {
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
        let component = momentum_component_atom(&edge.momentum, spatial_index);
        q3q3 += component.clone() * component;
    }
    (edge.mass_squared.clone() + q3q3).sqrt()
}

fn momentum_component_atom(momentum: &Atom, spatial_index: usize) -> Atom {
    let component_index = Atom::from(ExpandedIndex::from_iter([spatial_index]));
    momentum
        .replace(function!(GS.emr_mom, W_.i_))
        .with(function!(GS.emr_mom, W_.i_, component_index))
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
    source_context: Expanded4DSourceContext<'_>,
) -> Result<Expanded4DParsedSource> {
    let mut parent = (0..graph.n_nodes()).collect_vec();
    for (pair, _edge_id, edge_data) in graph.underlying.iter_edges() {
        if edge_data.data.is_dummy || !source_context.cograph_contract_subgraph.includes(&pair) {
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
    let mut next_node_id = root_to_internal.len();
    let uv_node_to_internal = uv_leading_node_map(
        graph,
        denominators,
        source_context.uv_leading_subgraph,
        &mut next_node_id,
    )?;

    let active_loop_columns = active_loop_columns_for_denominators(graph, denominators)?;
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
    let mut edge_map_internal = BTreeMap::new();

    for denominator in denominators {
        let (_, pair) = graph[&denominator.source_edge];
        let HedgePair::Paired { source, sink } = pair else {
            return Err(eyre!(
                "expanded 4D denominator edge {} is not a paired internal edge",
                denominator.source_edge
            ));
        };
        let local_edge_id = internal_edges.len();
        let (source_node, sink_node) = if denominator_uses_uv_loop_basis(denominator) {
            let Some(uv_subgraph) = source_context.uv_leading_subgraph else {
                return Err(eyre!(
                    "expanded 4D local UV denominator {} uses the UV loop basis, but no UV source subgraph was provided",
                    denominator.source_edge
                ));
            };
            if !uv_subgraph.includes(&pair) {
                return Err(eyre!(
                    "expanded 4D local UV denominator {} uses the UV loop basis, but its source edge is not part of the UV-leading subgraph {}",
                    denominator.source_edge,
                    uv_subgraph.string_label()
                ));
            }
            let source_node =
                *uv_node_to_internal
                    .get(&graph.node_id(source))
                    .ok_or_else(|| {
                        eyre!(
                            "missing UV source-node mapping for edge {} in UV-leading subgraph {}",
                            denominator.source_edge,
                            uv_subgraph.string_label()
                        )
                    })?;
            let sink_node = *uv_node_to_internal
                .get(&graph.node_id(sink))
                .ok_or_else(|| {
                    eyre!(
                        "missing UV sink-node mapping for edge {} in UV-leading subgraph {}",
                        denominator.source_edge,
                        uv_subgraph.string_label()
                    )
                })?;
            (source_node, sink_node)
        } else {
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
            (source_node, sink_node)
        };
        let momentum_signature =
            denominator_momentum_signature(graph, denominator, &active_loop_columns)?;
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
            momentum: denominator.momentum.clone(),
            mass_squared: denominator.mass_squared.clone(),
            full_expr: denominator.full_expr.clone(),
            is_preserved_4d,
        });
        edge_map_internal.insert(local_edge_id, usize::from(denominator.source_edge));
    }

    let mut external_edges = graph
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
    complete_external_balance_edges(&mut external_edges, &internal_edges, &external_names);

    let parsed = ParsedGraph {
        internal_edges,
        external_edges,
        initial_state_cut_edges: Vec::new(),
        loop_names,
        external_names,
        node_name_to_internal: root_to_internal
            .into_iter()
            .map(|(root, node)| (format!("n{root}"), node))
            .chain(
                uv_node_to_internal
                    .into_iter()
                    .map(|(node_id, node)| (format!("uv{}", usize::from(node_id)), node)),
            )
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

fn uv_leading_node_map(
    graph: &Graph,
    denominators: &[ExpandedDenominator],
    uv_subgraph: Option<&SuBitGraph>,
    next_node_id: &mut usize,
) -> Result<BTreeMap<linnet::half_edge::NodeIndex, usize>> {
    let Some(uv_subgraph) = uv_subgraph else {
        return Ok(BTreeMap::new());
    };

    let mut uv_nodes = BTreeSet::new();
    for denominator in denominators
        .iter()
        .filter(|denominator| denominator_uses_uv_loop_basis(denominator))
    {
        let (_, pair) = graph[&denominator.source_edge];
        let HedgePair::Paired { source, sink } = pair else {
            return Err(eyre!(
                "expanded 4D UV-leading denominator edge {} is not a paired internal edge",
                denominator.source_edge
            ));
        };
        if !uv_subgraph.includes(&pair) {
            return Err(eyre!(
                "expanded 4D UV-leading denominator edge {} is outside UV-leading subgraph {}",
                denominator.source_edge,
                uv_subgraph.string_label()
            ));
        }
        uv_nodes.insert(graph.node_id(source));
        uv_nodes.insert(graph.node_id(sink));
    }

    let mut uv_node_to_internal = BTreeMap::new();
    for node in uv_nodes {
        uv_node_to_internal.insert(node, *next_node_id);
        *next_node_id += 1;
    }
    Ok(uv_node_to_internal)
}

fn complete_external_balance_edges(
    external_edges: &mut Vec<ParsedGraphExternalEdge>,
    internal_edges: &[ParsedGraphInternalEdge],
    external_names: &[String],
) {
    let external_count = external_names.len();
    let mut balances = BTreeMap::<usize, Vec<i32>>::new();
    for edge in internal_edges {
        balances
            .entry(edge.tail)
            .or_insert_with(|| vec![0; external_count]);
        balances
            .entry(edge.head)
            .or_insert_with(|| vec![0; external_count]);
        for (external_id, coeff) in edge.signature.external_signature.iter().enumerate() {
            balances.get_mut(&edge.tail).unwrap()[external_id] -= coeff;
            balances.get_mut(&edge.head).unwrap()[external_id] += coeff;
        }
    }
    for edge in external_edges.iter() {
        if let Some(source) = edge.source {
            balances
                .entry(source)
                .or_insert_with(|| vec![0; external_count]);
            for (external_id, coeff) in edge.external_coefficients.iter().enumerate() {
                balances.get_mut(&source).unwrap()[external_id] += coeff;
            }
        }
        if let Some(destination) = edge.destination {
            balances
                .entry(destination)
                .or_insert_with(|| vec![0; external_count]);
            for (external_id, coeff) in edge.external_coefficients.iter().enumerate() {
                balances.get_mut(&destination).unwrap()[external_id] += coeff;
            }
        }
    }

    let mut next_external_id = external_edges
        .iter()
        .map(|edge| edge.edge_id)
        .max()
        .map(|edge_id| edge_id + 1)
        .unwrap_or(0);
    for (node, balance) in balances {
        for (external_id, coeff) in balance.into_iter().enumerate() {
            let name = external_names
                .get(external_id)
                .cloned()
                .unwrap_or_else(|| format!("p{external_id}"));
            for _ in 0..coeff.max(0) {
                let mut external_coefficients = vec![0; external_count];
                external_coefficients[external_id] = -1;
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_id,
                    source: Some(node),
                    destination: None,
                    label: format!("-{name}"),
                    external_coefficients,
                });
                next_external_id += 1;
            }
            for _ in 0..(-coeff).max(0) {
                let mut external_coefficients = vec![0; external_count];
                external_coefficients[external_id] = 1;
                external_edges.push(ParsedGraphExternalEdge {
                    edge_id: next_external_id,
                    source: None,
                    destination: Some(node),
                    label: name.clone(),
                    external_coefficients,
                });
                next_external_id += 1;
            }
        }
    }
}

fn active_loop_columns_for_denominators(
    graph: &Graph,
    denominators: &[ExpandedDenominator],
) -> Result<Vec<usize>> {
    let rows = denominators
        .iter()
        .map(|denominator| {
            Ok(full_denominator_momentum_signature(graph, denominator)?
                .loop_signature
                .into_iter()
                .collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(independent_loop_columns(&rows))
}

fn denominator_momentum_signature(
    graph: &Graph,
    denominator: &ExpandedDenominator,
    active_loop_columns: &[usize],
) -> Result<MomentumSignature> {
    let full_signature = full_denominator_momentum_signature(graph, denominator)?;
    let signature = MomentumSignature {
        loop_signature: active_loop_columns
            .iter()
            .map(|loop_index| {
                full_signature
                    .loop_signature
                    .get(*loop_index)
                    .copied()
                    .unwrap_or_default()
            })
            .collect(),
        external_signature: full_signature.external_signature,
    };
    Ok(signature)
}

fn full_denominator_momentum_signature(
    graph: &Graph,
    denominator: &ExpandedDenominator,
) -> Result<MomentumSignature> {
    let mut loop_signature = vec![0; graph.loop_momentum_basis.loop_edges.len()];
    let mut external_signature = vec![0; graph.loop_momentum_basis.ext_edges.len()];
    accumulate_momentum_signature_from_view(
        graph,
        denominator.momentum.as_view(),
        1,
        denominator_uses_uv_loop_basis(denominator),
        &mut loop_signature,
        &mut external_signature,
    )
    .map_err(|error| {
        eyre!(
            "could not extract expanded 4D denominator momentum signature for edge {} from `{}`: {error}",
            denominator.source_edge,
            denominator.momentum.log_print(None)
        )
    })?;
    Ok(MomentumSignature {
        loop_signature,
        external_signature,
    })
}

fn accumulate_momentum_signature_from_view(
    graph: &Graph,
    view: AtomView<'_>,
    coefficient: i32,
    use_uv_loop_basis: bool,
    loop_signature: &mut [i32],
    external_signature: &mut [i32],
) -> Result<()> {
    match view {
        AtomView::Add(add) => {
            for term in add.iter() {
                accumulate_momentum_signature_from_view(
                    graph,
                    term,
                    coefficient,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            }
        }
        AtomView::Mul(mul) => {
            let mut scalar = coefficient;
            let mut momentum_factor = None;
            for factor in mul.iter() {
                if let Ok(integer) = i64::try_from(factor) {
                    scalar = scalar.checked_mul(i32::try_from(integer)?).ok_or_else(|| {
                        eyre!("integer coefficient overflow in expanded denominator momentum")
                    })?;
                } else if momentum_factor.replace(factor).is_some() {
                    return Err(eyre!(
                        "expected a linear momentum expression, found product `{}`",
                        view.to_owned().log_print(None)
                    ));
                }
            }
            if let Some(momentum_factor) = momentum_factor {
                accumulate_momentum_signature_from_view(
                    graph,
                    momentum_factor,
                    scalar,
                    use_uv_loop_basis,
                    loop_signature,
                    external_signature,
                )?;
            } else if scalar != 0 {
                return Err(eyre!(
                    "expected a momentum factor, found scalar `{}`",
                    view.to_owned().log_print(None)
                ));
            }
        }
        AtomView::Fun(function) => {
            if function.get_symbol() != GS.emr_mom || function.get_nargs() != 1 {
                return Err(eyre!(
                    "expected an EMR momentum, found `{}`",
                    view.to_owned().log_print(None)
                ));
            }
            let edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                eyre!(
                    "EMR momentum has non-integer edge id `{}`",
                    function.get(0).to_owned().log_print(None)
                )
            })?);
            if use_uv_loop_basis {
                let loop_index = graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .position(|loop_edge| *loop_edge == edge)
                    .ok_or_else(|| {
                        eyre!(
                            "UV-leading denominator momentum references edge {edge}, which is not a loop-basis edge"
                        )
                    })?;
                loop_signature[loop_index] += coefficient;
            } else {
                let signature = graph
                    .loop_momentum_basis
                    .edge_signatures
                    .get(edge)
                    .ok_or_else(|| {
                        eyre!("edge {edge} is not present in the loop-momentum basis")
                    })?;
                for (loop_index, sign) in signature.internal.iter_enumerated() {
                    loop_signature[usize::from(loop_index)] += coefficient * sign_to_i32(*sign);
                }
                for (external_index, sign) in signature.external.iter_enumerated() {
                    external_signature[usize::from(external_index)] +=
                        coefficient * sign_to_i32(*sign);
                }
            }
        }
        AtomView::Num(_) => {
            if !view.to_owned().is_zero() {
                return Err(eyre!(
                    "expected a momentum expression, found scalar `{}`",
                    view.to_owned().log_print(None)
                ));
            }
        }
        _ => {
            return Err(eyre!(
                "expected a linear momentum expression, found `{}`",
                view.to_owned().log_print(None)
            ));
        }
    }
    Ok(())
}

fn denominator_uses_uv_loop_basis(denominator: &ExpandedDenominator) -> bool {
    denominator.mass_squared == Atom::var(GS.m_uv).pow(2)
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
        momentum: function.get(1).to_owned(),
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
