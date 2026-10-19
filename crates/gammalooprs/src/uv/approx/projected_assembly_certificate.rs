//! Test-only observations of the actual final projected assembly.
//!
//! Keys are retained before argument pruning/materialization loses information.
//! The captured chunks already contain their marker, tree and frozen factors;
//! replay must only apply the ordinary final normalizer.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use color_eyre::{
    Result,
    eyre::{bail, ensure, eyre},
};
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::subgraph::{SuBitGraph, SubSetOps};
use spenso::structure::IndexedStructure;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    id::Replacement,
    symbol,
};

use super::FinalIntegrandBuilder;
use crate::integrands::process::param_builder::NumeratorDefinition;
use crate::uv::approx::direct_3d::DirectResidueKey;
use crate::uv::{
    UVgenerationSettings, UltravioletGraph,
    approx::{ForestNodeLike, local_3d::Localizer, projected_4d::Projected4dCts},
    marker::UvMarker,
};
use crate::{
    cff::CutCFFIndex, graph::Graph, numerator::symbolica_ext::NumeratorAtomExt, utils::GS,
    uv::Integrands,
};

#[derive(Clone, Debug)]
pub(crate) enum ProjectedAssemblySource {
    Local(Box<ProjectedLocalAssembly>),
    Integrated { frozen_integrands: Integrands },
}

#[derive(Clone, Debug)]
pub(crate) struct ProjectedLocalAssembly {
    pub(crate) sector: usize,
    pub(crate) child: Arc<NumeratorDefinition>,
    pub(crate) prepared: Arc<NumeratorDefinition>,
    pub(crate) routed_numerator: SymbolicTensor<IndexedStructure>,
    pub(crate) routed_carrier: Atom,
    pub(crate) call_weights: BTreeMap<Atom, Atom>,
    /// Actual weight-to-call replacements, before scalar cancellation.
    pub(crate) prepared_calls: BTreeMap<Atom, Atom>,
    pub(crate) outer_arguments: Vec<Atom>,
    pub(crate) frozen_factor: Atom,
    pub(crate) outer: Integrands,
    pub(crate) mapped_carrier: Atom,
}

#[derive(Clone, Debug)]
pub(crate) struct ProjectedAssemblyEmission {
    pub(crate) source: ProjectedAssemblySource,
    pub(crate) key: DirectResidueKey,
    pub(crate) tree_factor: Atom,
    pub(crate) mapped: Integrands,
}

#[derive(Clone, Debug)]
pub(crate) enum ProjectedAssemblySkip {
    ZeroSector { sector: usize },
    ZeroCarrier { sector: usize, child: Atom },
    EmptyOuter { sector: usize, child: Atom },
    EmptyIntegrated,
}

#[derive(Clone, Default, Debug)]
pub(crate) struct ProjectedAssemblyCertificate {
    pub(crate) sector_count: usize,
    pub(crate) visited_families: Vec<(usize, Arc<NumeratorDefinition>)>,
    pub(crate) local_key_domains: Vec<(usize, Atom, Vec<DirectResidueKey>)>,
    pub(crate) integrated_key_domain: Vec<DirectResidueKey>,
    pub(crate) integrated_is_zero: Option<bool>,
    pub(crate) emissions: Vec<ProjectedAssemblyEmission>,
    pub(crate) skips: Vec<ProjectedAssemblySkip>,
}

impl ProjectedAssemblyCertificate {
    /// Certify the supported nonzero, unpruned assembly without resolving a
    /// numerator body. Inner child-contour provenance is a separate certificate.
    pub(crate) fn certify<S: ForestNodeLike>(
        &self,
        graph: &Graph,
        current: &S,
        projected: &Projected4dCts,
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
        actual: &Integrands,
    ) -> Result<()> {
        let mut empty_integrated = 0usize;
        for skip in &self.skips {
            match skip {
                ProjectedAssemblySkip::ZeroSector { sector } => {
                    bail!("assembly certificate does not support zero sector {sector}")
                }
                ProjectedAssemblySkip::ZeroCarrier { sector, child } => bail!(
                    "assembly certificate does not support zero carrier in sector {sector}: {child}"
                ),
                ProjectedAssemblySkip::EmptyOuter { sector, child } => bail!(
                    "assembly certificate does not support empty outer domain in sector {sector}: {child}"
                ),
                ProjectedAssemblySkip::EmptyIntegrated => empty_integrated += 1,
            }
        }
        ensure!(
            self.sector_count == projected.sectors().len() && self.sector_count > 0,
            "assembly sector inventory differs from projected input"
        );
        let expected_families = projected
            .sectors()
            .iter()
            .enumerate()
            .flat_map(|(sector, input)| {
                input
                    .coefficient
                    .numerators()
                    .iter()
                    .map(move |child| (sector, child))
            })
            .collect::<Vec<_>>();
        ensure!(
            !expected_families.is_empty() && self.visited_families.len() == expected_families.len(),
            "assembly family inventory is incomplete"
        );
        ensure!(
            self.local_key_domains.len() == expected_families.len(),
            "assembly local key-domain inventory is incomplete"
        );
        for (((sector, child), (visited_sector, visited)), (domain_sector, domain_child, keys)) in
            expected_families
                .iter()
                .zip(&self.visited_families)
                .zip(&self.local_key_domains)
        {
            ensure!(
                sector == visited_sector && Arc::ptr_eq(child, visited),
                "assembly visited a different original child family"
            );
            ensure!(
                sector == domain_sector && child.lhs == *domain_child && !keys.is_empty(),
                "assembly key domain has a different sector/family owner"
            );
            ensure!(
                keys.iter()
                    .enumerate()
                    .all(|(i, key)| !keys[..i].contains(key)),
                "assembly local key inventory contains a duplicate"
            );
            let emitted = self
                .emissions
                .iter()
                .filter_map(|emission| match &emission.source {
                    ProjectedAssemblySource::Local(local)
                        if local.sector == *sector && Arc::ptr_eq(&local.child, child) =>
                    {
                        Some(&emission.key)
                    }
                    _ => None,
                })
                .collect::<Vec<_>>();
            ensure!(
                emitted == keys.iter().collect::<Vec<_>>(),
                "assembly did not emit every actual local key exactly once in order for sector {sector}"
            );
        }
        let integrated_keys = self
            .emissions
            .iter()
            .filter_map(|emission| {
                matches!(&emission.source, ProjectedAssemblySource::Integrated { .. })
                    .then_some(&emission.key)
            })
            .collect::<Vec<_>>();
        ensure!(
            integrated_keys == self.integrated_key_domain.iter().collect::<Vec<_>>()
                && empty_integrated == usize::from(self.integrated_key_domain.is_empty()),
            "assembly integrated domain/emissions and exact empty-branch decision disagree"
        );
        ensure!(
            self.integrated_key_domain
                .iter()
                .enumerate()
                .all(|(i, key)| !self.integrated_key_domain[..i].contains(key)),
            "assembly integrated key inventory contains a duplicate"
        );

        let allowed_keys = localizer.cutset.residue_selector.generate_allowed_keys();
        let integrated_is_zero = self
            .integrated_is_zero
            .ok_or_else(|| eyre!("assembly did not record its actual integrated input"))?;
        let allowed = allowed_keys.iter().copied().collect::<BTreeSet<_>>();
        ensure!(
            !allowed.is_empty(),
            "assembly certificate requires a nonempty cut domain"
        );
        let check_cuts = |integrands: &Integrands| -> Result<()> {
            ensure!(
                integrands
                    .iter()
                    .map(|(key, _)| *key)
                    .collect::<BTreeSet<_>>()
                    == allowed,
                "assembly chunk has incomplete cut-key support"
            );
            Ok(())
        };
        let orientations = localizer.orientation.exact_orientations()?;
        let hosts = localizer.orientation.orientation_ids()?;
        let tree = GS.wrap_tree_denoms(graph.denominator(
            &graph.tree_edges.subtract(&graph.initial_state_cut),
            &graph.model,
            |_| -1,
        ));
        let full_graph = graph.full_filter();
        let marker = UvMarker::new(settings);
        let family = symbol!("gammalooprs::uv::numerator_family");
        for emission in &self.emissions {
            ensure!(
                hosts.contains(&emission.key.selector_host),
                "assembly key has an unselected production host"
            );
            let production = orientations
                .get(emission.key.selector_host)
                .ok_or_else(|| eyre!("assembly key has no production host"))?;
            let map = emission
                .key
                .source_edge_energy_map()
                .unwrap_or(&production.edge_energy_map);
            ensure!(
                map.len() == graph.underlying.n_edges(),
                "assembly key has {} edge rows; actual graph has {}",
                map.len(),
                graph.underlying.n_edges()
            );
            ensure!(
                emission.tree_factor == tree,
                "assembly tree factor differs from the full actual production tree"
            );
            check_cuts(&emission.mapped)?;
            match &emission.source {
                ProjectedAssemblySource::Integrated { frozen_integrands } => {
                    check_cuts(frozen_integrands)?;
                    let kernel = if integrated_is_zero {
                        Atom::one()
                    } else {
                        GS.localizing_integrand(current.lmb())
                    };
                    ensure!(
                        frozen_integrands.iter().all(|(_, atom)| *atom == kernel),
                        "integrated kernel differs from its actual zero/LMB policy"
                    );
                    if integrated_is_zero {
                        ensure!(
                            emission.key.source_edge_energy_map().is_none()
                                && emission.mapped.iter().all(|(_, atom)| atom.is_zero()),
                            "the zero integrated input must retain production-key zeros with a unit kernel"
                        );
                    }
                }
                ProjectedAssemblySource::Local(local) => {
                    let input = projected.sectors().get(local.sector).ok_or_else(|| {
                        eyre!("assembly local sector is outside the input domain")
                    })?;
                    ensure!(
                        input
                            .coefficient
                            .numerators()
                            .iter()
                            .any(|child| Arc::ptr_eq(child, &local.child)),
                        "emission child is absent from its actual sector"
                    );
                    ensure!(
                        local.frozen_factor == input.frozen_factor,
                        "assembly changed a sector's frozen kernel"
                    );
                    ensure!(
                        local.child.tags.len() == 1 && local.prepared.tags.len() == 1,
                        "assembly certificate currently requires one fixed family tag"
                    );
                    let child_parameters = local
                        .child
                        .args
                        .iter()
                        .cloned()
                        .map(Atom::from)
                        .collect::<Vec<_>>();
                    let parameters = local
                        .prepared
                        .args
                        .iter()
                        .cloned()
                        .map(Atom::from)
                        .collect::<Vec<_>>();
                    ensure!(
                        parameters.starts_with(&child_parameters)
                            && parameters.len()
                                == child_parameters.len() + local.outer_arguments.len(),
                        "prepared family lost ordered child/outer parameters"
                    );
                    ensure!(
                        local.prepared.inlining == local.child.inlining
                            && local.prepared.is_alias == local.child.is_alias,
                        "prepared family changed definition policy"
                    );
                    let mut calls = BTreeSet::new();
                    for (_, root) in input.coefficient.iter() {
                        root.visitor(&mut |view| {
                            if let AtomView::Fun(call) = view
                                && call.get_symbol() == family
                            {
                                if call.get_nargs() > 0
                                    && call.get(0) == local.child.tags[0].as_view()
                                {
                                    calls.insert(view.to_owned());
                                }
                                return false;
                            }
                            true
                        });
                    }
                    ensure!(
                        !calls.is_empty() && calls == local.call_weights.keys().cloned().collect(),
                        "recorded child-call occurrence inventory differs from the actual sector"
                    );
                    ensure!(
                        local.call_weights.values().collect::<BTreeSet<_>>().len() == calls.len(),
                        "child calls share a diagnostic weight"
                    );
                    ensure!(
                        local.prepared_calls.keys().collect::<BTreeSet<_>>()
                            == local.call_weights.values().collect::<BTreeSet<_>>(),
                        "actual prepared-call inventory differs from the complete weight inventory"
                    );
                    let mut replacements = Vec::with_capacity(calls.len());
                    for (call, weight) in &local.call_weights {
                        let AtomView::Fun(call_view) = call.as_view() else {
                            bail!("recorded child call is not a function");
                        };
                        let child_arguments = call_view
                            .iter()
                            .skip(local.child.tags.len())
                            .map(|argument| argument.to_owned())
                            .collect::<Vec<_>>();
                        ensure!(
                            child_arguments.len() == child_parameters.len(),
                            "recorded child call has an incorrect ordered tuple length"
                        );
                        let mapped_arguments = child_arguments
                            .iter()
                            .map(|argument| {
                                emission
                                    .key
                                    .map_numerator(localizer.orientation, graph, argument)
                            })
                            .collect::<Result<Vec<_>>>()?;
                        let arguments = mapped_arguments
                            .into_iter()
                            .chain(local.outer_arguments.iter().cloned())
                            .collect::<Vec<_>>();
                        let expected_call = local.prepared.lhs.replace_multiple(
                            parameters
                                .iter()
                                .zip(&arguments)
                                .map(|(parameter, argument)| {
                                    Replacement::new(parameter.to_pattern(), argument.clone())
                                }),
                        );
                        let actual_call = local.prepared_calls.get(weight).ok_or_else(|| {
                            eyre!("actual prepared call is missing for a retained child occurrence")
                        })?;
                        ensure!(
                            actual_call == &expected_call,
                            "prepared child-call tuple differs before scalar cancellation in sector {}, host {}: actual {}; expected {}",
                            local.sector,
                            emission.key.selector_host.0,
                            actual_call
                                .to_plain_string()
                                .chars()
                                .take(160)
                                .collect::<String>(),
                            expected_call
                                .to_plain_string()
                                .chars()
                                .take(160)
                                .collect::<String>()
                        );
                        replacements.push(Replacement::new(weight.to_pattern(), expected_call));
                    }
                    let expected_carrier = emission
                        .key
                        .map_numerator(localizer.orientation, graph, &local.routed_carrier)?
                        .replace_multiple(replacements);
                    ensure!(
                        expected_carrier == local.mapped_carrier,
                        "mapped child-call tuple/carrier differs for sector {}, host {} (actual {} bytes, expected {} bytes)",
                        local.sector,
                        emission.key.selector_host.0,
                        local.mapped_carrier.as_view().get_byte_size(),
                        expected_carrier.as_view().get_byte_size()
                    );
                    let specialized = local.prepared.rhs.with_checked_expression(
                        local.prepared.rhs.expression().replace_multiple(
                            parameters[child_parameters.len()..]
                                .iter()
                                .zip(&local.outer_arguments)
                                .map(|(parameter, argument)| {
                                    Replacement::new(parameter.to_pattern(), argument.clone())
                                }),
                        ),
                    )?;
                    let expected = local.routed_numerator.with_checked_expression(
                        emission.key.map_numerator(
                            localizer.orientation,
                            graph,
                            local.routed_numerator.expression(),
                        )?,
                    )?;
                    ensure!(
                        specialized == expected,
                        "prepared outer argument tuple differs from its actual full key in sector {}, host {}",
                        local.sector,
                        emission.key.selector_host.0
                    );
                    check_cuts(&local.outer)?;
                    let reconstructed = local
                        .outer
                        .map(|atom| {
                            marker.prefix(
                                &full_graph,
                                current.subgraph(),
                                &(atom * &expected_carrier * &local.frozen_factor),
                            )
                        })
                        .with_numerators(
                            local
                                .outer
                                .numerators()
                                .iter()
                                .cloned()
                                .chain([Arc::clone(&local.prepared)]),
                        )?;
                    ensure!(
                        reconstructed == emission.mapped,
                        "recorded local chunk differs from its factorized inputs in sector {}",
                        local.sector
                    );
                }
            }
        }
        let reduced = full_graph
            .subtract(current.subgraph())
            .subtract(&graph.initial_state_cut);
        check_cuts(actual)?;
        ensure!(
            self.finalized_view(graph, &reduced, &allowed_keys)? == *actual,
            "captured assembly chunks do not replay to the actual final factorized payload"
        );
        crate::debug_tags!(#uv, #local, #four_d, #cff, #summary;
            stage = "projected_assembly_certificate", sectors = self.sector_count,
            local_families = self.local_key_domains.len(),
            local_keys = self.local_key_domains.iter().map(|(_, _, keys)| keys.len()).sum::<usize>(),
            integrated_keys = self.integrated_key_domain.len(),
            integrated_input_zero = integrated_is_zero,
            exact_empty_integrated = empty_integrated == 1,
            prepared_call_tuples_before_cancellation = true,
            file.local_key_domains = ?self.local_key_domains,
            file.integrated_key_domain = ?self.integrated_key_domain,
            final_factorized_replay = true, original_child_arguments_mapped = true,
            "Actual projected assembly preserves its complete outer key and ordered-call inventory"
        );
        Ok(())
    }

    /// Replay captured chunks through only the production finalizer suffix.
    /// All residue maps, markers, tree factors and kernels are already owned.
    pub(crate) fn finalized_view(
        &self,
        graph: &Graph,
        reduced: &SuBitGraph,
        allowed_keys: &[CutCFFIndex],
    ) -> Result<Integrands> {
        Ok(Self::finalize_chunks(
            graph,
            reduced,
            allowed_keys,
            self.emissions
                .iter()
                .map(|emission| emission.mapped.clone()),
        )?
        .map(|atom| atom.collect_compact_factors()))
    }

    /// Partition only the already captured emission, retaining every zero and
    /// full key by ordinal. Exact recomposition is required before linear Phi.
    pub(crate) fn finalized_views(
        &self,
        graph: &Graph,
        reduced: &SuBitGraph,
        allowed_keys: &[CutCFFIndex],
        actual: &Integrands,
    ) -> Result<Vec<Integrands>> {
        let views = self
            .emissions
            .iter()
            .map(|emission| {
                Self::finalize_chunks(graph, reduced, allowed_keys, [emission.mapped.clone()])
            })
            .collect::<Result<Vec<_>>>()?;
        let allowed_zero: Integrands = allowed_keys.iter().map(|key| (*key, Atom::Zero)).collect();
        let recomposed = allowed_zero
            .zip_add(views.iter().cloned())?
            .map(|atom| atom.collect_compact_factors());
        if recomposed != *actual {
            crate::debug_tags!(#uv, #local, #four_d, #cff, #inspect;
                stage = "projected_assembly_partition_refused",
                same_numerator_catalog = recomposed.numerators() == actual.numerators(),
                same_scalar_catalog = recomposed.scalar_definitions() == actual.scalar_definitions(),
                file.recomposed_roots = ?recomposed.iter().collect::<Vec<_>>(),
                file.actual_roots = ?actual.iter().collect::<Vec<_>>(),
                "Actual emission partition failed exact finalizer recomposition"
            );
        }
        ensure!(
            recomposed == *actual,
            "separately finalized actual emissions do not exactly recompose the returned factorized payload"
        );
        Ok(views)
    }

    /// One owner for the unchanged ordinary finalizer suffix. No key mapping,
    /// numerator materialization or denominator clearing happens here.
    fn finalize_chunks(
        graph: &Graph,
        reduced: &SuBitGraph,
        allowed_keys: &[CutCFFIndex],
        chunks: impl IntoIterator<Item = Integrands>,
    ) -> Result<Integrands> {
        let allowed_zero: Integrands = allowed_keys.iter().map(|key| (*key, Atom::Zero)).collect();
        let normalize_scalar = FinalIntegrandBuilder::scalar_normalizer(graph, reduced);
        let result = allowed_zero
            .zip_add(chunks)?
            .resolved_scalars()?
            .fallible_map(|atom| {
                Ok(normalize_scalar(atom)
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum)
                    .replace(GS.dim_epsilon)
                    .with(0))
            })?
            .map_numerators(FinalIntegrandBuilder::final_normalizer(graph, reduced))?;
        // The production collector sees the complete sum once. Applying it to
        // each chunk first can choose a different, equally factorized grouping.
        // Keep the additive chunks before that final collection boundary.
        Ok(result)
    }
}
