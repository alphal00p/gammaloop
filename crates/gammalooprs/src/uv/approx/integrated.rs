use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::debug_instrument;
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::{AGS, GammaSimplifySettings},
    representations::Bispinor,
    tensor::{
        ContractSettings, SymbolicNetExt, SymbolicNetParse, SymbolicTensor,
        inference::InterfaceInference,
    },
};

use linnet::half_edge::{
    HedgeGraph,
    builder::HedgeGraphBuilder,
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
};
use spenso::{
    network::{
        library::symbolic::ETS,
        parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
        tags::SPENSO_TAG,
    },
    structure::{
        IndexedStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::{DummyAind, ParseableAind, SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::{atom::AtomField, integer::Z, rational::Q},
    function,
    id::Replacement,
    parse, parse_lit,
    poly::{PolyVariable, series::Series},
};
use symbolica_utils::ReplaceBuilderExt;
use vakint::{Vakint, VakintExpression, vakint_symbol};

use crate::{
    debug_tags,
    graph::{Graph, LMBext, LoopMomentumBasis},
    numerator::{aind::Aind, ufo::UFO},
    utils::{GS, W_},
    uv::{
        ApproximationType, UltravioletGraph,
        approx::{
            ForestNodeLike, Rooted, UVCtx,
            local_4d::{FourDSector, Local4dCts},
        },
        marker::{UvMarker, UvOperation},
        settings::VakintSettings,
        uv_graph::UVE,
    },
};

/// Laurent projections of an integrated counterterm.
///
/// The native series engine computes the Laurent projections transiently. Both
/// retained projections keep the completed tensor's logical interface, including zero.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) struct IntegratedCts {
    pole: SymbolicTensor<IndexedStructure>,
    finite_counterterm: SymbolicTensor<IndexedStructure>,
    scale_power: i64,
}

impl IntegratedCts {
    fn from_series(
        source: &SymbolicTensor<IndexedStructure>,
        expansion: &Series<AtomField>,
        scale_power: i64,
    ) -> Result<Self> {
        Ok(Self {
            pole: source.with_checked_expression(truncate(expansion, false))?,
            finite_counterterm: source.with_checked_expression(-truncate(expansion, true))?,
            scale_power,
        })
    }

    pub(crate) fn factorized_product<'a>(
        factors: impl IntoIterator<Item = &'a Self>,
        depth: usize,
    ) -> Result<Self> {
        let mut factors = factors.into_iter();
        let first = factors
            .next()
            .ok_or_else(|| eyre!("a factorized integrated counterterm cannot be empty"))?;
        let mut pole = first.pole.clone();
        let mut finite_counterterm = first.finite_counterterm.clone();
        let mut scale_power = first.scale_power;

        for factor in factors {
            pole = pole.multiply(factor.pole.clone())?;
            finite_counterterm = finite_counterterm.multiply(factor.finite_counterterm.clone())?;
            scale_power += factor.scale_power;
        }

        // Native Laurent truncation owns precision and powers. Its retained
        // projections inherit the complete product interface, even if both vanish.
        let combined = pole.try_sub(&finite_counterterm)?;
        let expansion = series(combined.expression(), depth)?;
        Self::from_series(&combined, &expansion, scale_power)
    }

    fn projected(&self, finite: bool) -> Result<SymbolicTensor<IndexedStructure>> {
        // The enclosing forest Taylor operation rescales the vacuum mass in
        // this coefficient along with its active momenta. The bookkeeping
        // factor restores only consumed loop measures; physical projections
        // and the next integration set it to one.
        let projection = if finite {
            &self.finite_counterterm
        } else {
            &self.pole
        };
        Ok(projection.with_checked_expression(
            projection.expression() * Atom::var(GS.integrated_loop_scale).pow(self.scale_power),
        )?)
    }

    pub(crate) fn pole(&self) -> Result<SymbolicTensor<IndexedStructure>> {
        self.projected(false)
    }

    pub(crate) fn finite_counterterm(&self) -> Result<SymbolicTensor<IndexedStructure>> {
        self.projected(true)
    }

    pub(crate) fn physical_pole(&self) -> Result<SymbolicTensor<IndexedStructure>> {
        Ok(self.pole.with_checked_expression(
            self.pole
                .expression()
                .replace(GS.integrated_loop_scale)
                .with(Atom::one()),
        )?)
    }

    pub(crate) fn physical_finite_counterterm(&self) -> Result<SymbolicTensor<IndexedStructure>> {
        Ok(self.finite_counterterm.with_checked_expression(
            self.finite_counterterm
                .expression()
                .replace(GS.integrated_loop_scale)
                .with(Atom::one()),
        )?)
    }
}

fn series(expr: &Atom, depth: usize) -> Result<Series<AtomField>> {
    Ok(expr.series(GS.dim_epsilon, 0, depth)?)
}

fn truncate(series: &Series<AtomField>, finite: bool) -> Atom {
    let mut truncated = Atom::Zero;

    for (power, p) in series.terms() {
        if (power >= 0) == finite {
            truncated += p * Atom::var(GS.dim_epsilon).pow(power);
        }
    }

    truncated
}

impl Rooted for IntegratedCts {
    fn root() -> Self {
        let zero =
            SymbolicTensor::infer_indexed(Atom::Zero).expect("zero has a scalar tensor interface");
        Self {
            pole: zero.clone(),
            finite_counterterm: zero,
            scale_power: 0,
        }
    }
}

impl Graph {
    /// Certify the minimal-QCD two-loop insertion for which the parent soft
    /// jet integrates to zero: a quark self-energy in an equal-mass quark
    /// bubble with two external gluons. Scheme and completed-child ownership
    /// are checked by the caller; this certifies topology and actual rules.
    pub(crate) fn has_qcd_self_energy_insertion_ward_identity(
        &self,
        model: &crate::model::Model,
        parent: &SuBitGraph,
        child: &SuBitGraph,
    ) -> bool {
        if parent == child
            || !parent.includes(child)
            || self.n_loops(parent) != 2
            || self.n_loops(child) != 1
            || self.underlying.connected_components(parent).len() != 1
            || self.underlying.connected_components(child).len() != 1
        {
            return false;
        }
        let child_crown = self.dummy_less_full_crown(child);
        let boundary = child_crown.included_iter().collect::<Vec<_>>();
        if boundary.len() != 2 || self.flow(boundary[0]) == self.flow(boundary[1]) {
            return false;
        }
        let reference = &self[self[&boundary[0]]];
        let Some(quark) = reference
            .particle
            .particle()
            .and_then(|id| model.particle_by_id(id).ok())
        else {
            return false;
        };
        let mass = reference.mass_atom(model);
        let flavour = quark.pdg_code.abs();
        if mass.is_zero()
            || !quark.is_spinor()
            || quark.color.abs() != 3
            || !(1..=6).contains(&flavour)
        {
            return false;
        }
        let is_quark = |edge: &crate::graph::Edge| {
            edge.mass_atom(model) == mass
                && edge
                    .particle
                    .particle()
                    .and_then(|id| model.particle_by_id(id).ok())
                    .is_some_and(|particle| {
                        particle.is_spinor()
                            && particle.color.abs() == 3
                            && particle.pdg_code.abs() == flavour
                    })
        };
        let is_gluon = |edge: &crate::graph::Edge| {
            edge.mass_atom(model).is_zero()
                && edge
                    .particle
                    .particle()
                    .and_then(|id| model.particle_by_id(id).ok())
                    .is_some_and(|particle| {
                        particle.pdg_code == 21 && particle.is_vector() && particle.color == 8
                    })
        };
        if !boundary.iter().all(|hedge| is_quark(&self[self[hedge]])) {
            return false;
        }
        let parent_crown = self.dummy_less_full_crown(parent);
        if parent_crown.included_iter().count() != 2
            || !parent_crown
                .included_iter()
                .all(|hedge| is_gluon(&self[self[&hedge]]))
        {
            return false;
        }
        let child_edges = self.iter_edges_of(child).collect::<Vec<_>>();
        if child_edges.len() != 2
            || child_edges
                .iter()
                .filter(|(_, _, edge)| is_quark(edge.data))
                .count()
                != 1
            || child_edges
                .iter()
                .filter(|(_, _, edge)| is_gluon(edge.data))
                .count()
                != 1
        {
            return false;
        }
        let cograph = parent.subtract(child);
        let cograph_edges = self.iter_edges_of(&cograph).collect::<Vec<_>>();
        if cograph_edges.len() != 3 || !cograph_edges.iter().all(|(_, _, edge)| is_quark(edge.data))
        {
            return false;
        }
        // Compare actual tensors in the source frame. Native and legacy
        // propagators can differ only by compact slashes or dummy labels;
        // the shared typed algebra certifies that normal-form equivalence.
        let same_tensor = |expected: Atom, actual: &SymbolicTensor<IndexedStructure>| {
            let certificate = (|| -> Result<bool> {
                let expected =
                    crate::graph::parse::NumeratorInterface::admit(expected.clone(), None)?
                        .try_into_indexed()?;
                let difference = expected.try_sub(actual)?;
                Ok(simplify(&difference)?
                    .coefficients_are_zero(spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
                    .is_true())
            })();
            if !certificate.as_ref().is_ok_and(|value| *value) {
                debug_tags!(#uv, #integrated, #ward, #inspect;
                    stage = "qcd_insertion_ward_tensor_certificate",
                    graph = %self.name,
                    result = ?certificate,
                    file.expected = %expected.to_canonical_string(),
                    file.actual = %actual.expression().to_canonical_string(),
                    "Actual source tensor did not certify the minimal QCD Ward rule"
                );
            }
            certificate.unwrap_or(false)
        };
        if !self.iter_nodes_of(parent).all(|(node, hedges, vertex)| {
            let Some(rule) = vertex
                .vertex_rule
                .and_then(|id| model.vertex_rule_by_id(id).ok())
            else {
                return false;
            };
            if rule.particles.len() != 3 {
                return false;
            }
            let position = |pdg| {
                rule.particles
                    .iter()
                    .position(|id| {
                        model
                            .particle_by_id(*id)
                            .is_ok_and(|particle| particle.pdg_code == pdg)
                    })
                    .map(|index| index + 1)
            };
            let (Some(anti), Some(fermion), Some(gluon)) =
                (position(-flavour), position(flavour), position(21))
            else {
                return false;
            };
            let minimal_spin = function!(
                UFO.gamma,
                UFO.idx(gluon, 1),
                UFO.idx(anti, 1),
                UFO.idx(fermion, 1)
            );
            let minimal_color = function!(UFO.t, gluon, fermion, anti);
            if rule.lorentz_structures.len() != 1
                || !model
                    .lorentz_structure_by_id(rule.lorentz_structures[0])
                    .is_ok_and(|lorentz| lorentz.structure == minimal_spin)
                || rule.color_structures != [minimal_color.clone()]
                || rule.couplings.len() != 1
                || rule.couplings[0].len() != 1
            {
                return false;
            }
            let Some(coupling) = rule.couplings[0][0].and_then(|id| model.coupling_by_id(id).ok())
            else {
                return false;
            };
            let Ok(node_label) = u16::try_from(node.0) else {
                return false;
            };
            let mut hedges = hedges.collect::<Vec<_>>();
            hedges.sort_by_key(|hedge| self.underlying[*hedge].ufo_order.value);
            if hedges.len() != 3
                || hedges.iter().enumerate().any(|(slot, hedge)| {
                    usize::from(self.underlying[*hedge].ufo_order.value) != slot
                })
            {
                return false;
            }
            let spin_slots = hedges
                .iter()
                .map(|hedge| {
                    &self.underlying[*hedge]
                        .num_indices
                        .spin_indices
                        .vertex_indices
                })
                .collect::<Vec<_>>();
            let color_slots = hedges
                .iter()
                .map(|hedge| {
                    &self.underlying[*hedge]
                        .num_indices
                        .color_indices
                        .vertex_indices
                })
                .collect::<Vec<_>>();
            let momenta = hedges
                .iter()
                .map(|hedge| (self.flow(*hedge), self[hedge]))
                .collect::<Vec<_>>();
            let Ok(spin) = UFO.reindex_spin(&spin_slots, &momenta, minimal_spin, |index| {
                Aind::Vertex(node_label, index as u16)
            }) else {
                return false;
            };
            let Ok(color) = UFO.reindex_color(&color_slots, minimal_color, |index| {
                Aind::Vertex(node_label, index as u16)
            }) else {
                return false;
            };
            same_tensor(
                Atom::var(symbolica::symbol!(&format!("UFO::{}", coupling.name))) * color * spin,
                &vertex.num.value,
            )
        }) {
            return false;
        }
        self.iter_edges_of(parent).all(|(pair, edge_id, edge)| {
            let HedgePair::Paired { source, sink } = pair else {
                return false;
            };
            if edge
                .data
                .extra_data
                .vakint_edge_power
                .is_some_and(|power| power != 1)
            {
                return false;
            }
            let first = UFO.idx(1, 1);
            let second = UFO.idx(2, 1);
            let propagator = if is_quark(edge.data) {
                let (left, right, sign) = if model
                    .particle_by_id(edge.data.particle.particle().unwrap())
                    .unwrap()
                    .pdg_code
                    > 0
                {
                    (second, first, 1)
                } else {
                    (first, second, -1)
                };
                Atom::i()
                    * (Atom::num(sign) * function!(UFO.pslash, &left, &right)
                        + &mass * function!(UFO.identity, left, right))
            } else if is_gluon(edge.data) {
                -Atom::i() * function!(UFO.metric, first, second)
            } else {
                return false;
            };
            let spin_slots = [
                &self.underlying[source]
                    .num_indices
                    .spin_indices
                    .edge_indices,
                &self.underlying[sink].num_indices.spin_indices.edge_indices,
            ];
            UFO.reindex_spin(
                &spin_slots,
                &[(Flow::Source, edge_id), (Flow::Sink, edge_id)],
                propagator,
                |index| Aind::Edge(edge_id.0 as u16, index as u16),
            )
            .is_ok_and(|spin| {
                same_tensor(
                    spin * self.underlying[source].color_kronekers(&self.underlying[sink]),
                    &edge.data.num.value,
                )
            })
        })
    }
}

pub(super) fn simplify(
    integrand: &SymbolicTensor<IndexedStructure>,
) -> Result<SymbolicTensor<IndexedStructure>> {
    let started = std::time::Instant::now();
    debug_tags!(#uv, #integrated, #diagnostic_stage;
        stage = 1u8, finished = false,
        input_bytes = integrand.expression().as_view().get_byte_size(),
        "Start retained analytic tensor simplification"
    );
    if integrand.expression().contains_symbol(AGS.sigma) {
        return Err(eyre!(
            "Sigma tensors are not supported by d-dimensional analytic UV numerator algebra"
        ));
    }
    // Only the analytically integrated UV subgraph reaches this boundary. The
    // shared tensor owner keeps compact scalar products and color spectators
    // factorized while closing spin chains, so powered dots need no temporary
    // expansion/protection registry here. The Vakint bridge exports the resolved
    // expression only after this carrier has completed native tensor algebra.
    let simplified = integrand
        .simplify_algebra(&idenso::tensor::AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        })?
        .contract(idenso::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })?;
    let dotted = simplified.to_dots()?;
    debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #gamma;
        log.expr = dotted.expression(),
        "Shared Dirac simplification of the retained analytic UV tensor"
    );

    if dotted
        .expression()
        .replace(function!(
            SPENSO_TAG.trace,
            Bispinor {}.to_symbolic([W_.d_]),
            W_.x___
        ))
        .matches()
    {
        return Err(eyre!(
            "an unresolved Dirac trace remains after d-dimensional UV numerator algebra"
        ));
    }
    if dotted
        .expression()
        .replace(function!(
            SPENSO_TAG.trace,
            Minkowski {}.to_symbolic([W_.d_]),
            W_.x___
        ))
        .matches()
    {
        return Err(eyre!(
            "an unresolved Lorentz trace remains after d-dimensional UV numerator algebra"
        ));
    }
    debug_tags!(#uv, #integrated, #diagnostic_stage;
        stage = 1u8, finished = true,
        elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
        output_bytes = dotted.expression().as_view().get_byte_size(),
        "Completed retained analytic tensor simplification"
    );
    Ok(dotted)
}

pub(crate) struct Integrated<'a> {
    pub vakint: &'a Vakint,
    pub vakint_settings: &'a vakint::VakintSettings,
}

impl Integrated<'_> {
    /// Keep numeric public ports distinct from Vakint's numeric dummy ABI.
    /// The returned map restores the original logical positions after the
    /// backend result has been checked against the protected carrier.
    pub(crate) fn protect_numeric_ports(
        source: &SymbolicTensor<IndexedStructure>,
    ) -> Result<(
        SymbolicTensor<IndexedStructure>,
        std::collections::HashMap<usize, AbstractIndex>,
    )> {
        let originals = source
            .structure()
            .logical_slots()
            .into_iter()
            .enumerate()
            .filter_map(|(position, slot)| {
                matches!(slot.aind.to_atom().as_view(), AtomView::Num(_))
                    .then_some((position, slot.aind))
            })
            .collect::<std::collections::HashMap<_, _>>();
        if originals.is_empty() {
            return Ok((source.clone(), originals));
        }
        let reserved = SymbolicTensor::<IndexedStructure>::reserved_dummies([source]);
        let protected = originals
            .keys()
            .map(|&position| (position, reserved.fresh_index()))
            .collect();
        Ok((source.reindex_interface_ports(&protected)?, originals))
    }

    pub(crate) fn new<'a>(
        vakint: &'a Vakint,
        vakint_settings: &'a vakint::VakintSettings,
    ) -> Integrated<'a> {
        Integrated {
            vakint,
            vakint_settings,
        }
    }

    pub(crate) fn run<S: super::ForestNodeLike, M: super::ForestNodeLike>(
        &self,
        integrand: &Local4dCts,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
        marker_current: &M,
        marker_given: &M,
    ) -> Result<IntegratedCts> {
        let graph = ctx.graph;

        let n_loops = graph.n_loops(current.subgraph());

        let scheme = current.renormalization_scheme();
        match scheme {
            ApproximationType::MUV
            | ApproximationType::PolePart
            | ApproximationType::IR
            | ApproximationType::OS => {
                // The local operation owns its integration policy. In the OS
                // scheme this is the signed vacuum remainder -(U - UO); the
                // physical O branch already supplies the on-shell mass shift
                // and must never enter the finite MS-bar addback. An enclosing
                // soft operation can supply a vacuum-only payload after its
                // integrated soft jet has been certified to vanish. Other
                // components retain their complete signed expression.
                let mut integrated = integrand.atom().with_checked_expression(Atom::Zero)?;
                for sector in integrand.active_sectors() {
                    let physical = sector
                        .integration_atom()
                        .expression()
                        .replace(GS.integrated_loop_scale)
                        .with(Atom::one());
                    let physical = GS.resolve_os_momenta(&physical);
                    let physical = sector.integration_atom().with_algebra_result(physical)?;
                    let simplified = simplify(&physical)?;
                    if !simplified.is_zero() {
                        let coefficient = self
                            .integrate(&simplified, sector, ctx, current, given)
                            .wrap_err_with(|| {
                                format!(
                                    "while integrating {scheme} counterterm {} with retained component frames {:?}",
                                    current.subgraph().string_label(),
                                    sector.active_components,
                                )
                            })?;
                        integrated = integrated.try_add(&coefficient)?;
                    }
                }
                let marker = UvMarker::new(ctx.settings);
                let marked = marker.apply(
                    UvOperation::Integrate,
                    marker_current.subgraph(),
                    marker_given.subgraph(),
                    integrated.expression(),
                );
                let integrated = integrated.with_checked_expression(marked)?;
                let expansion_depth =
                    usize::try_from(self.vakint_settings.number_of_terms_in_epsilon_expansion)
                        .wrap_err("Vakint epsilon expansion depth must be nonnegative")?;
                let expanded = series(integrated.expression(), expansion_depth.max(n_loops + 1))?
                    .map_coeff(|coefficient| {
                        marker.apply(
                            UvOperation::Series,
                            marker_current.subgraph(),
                            marker_given.subgraph(),
                            coefficient,
                        )
                    });
                let expansion = expanded.map_coeff(|coefficient| {
                    marker.apply(
                        UvOperation::Truncate,
                        marker_current.subgraph(),
                        marker_given.subgraph(),
                        coefficient,
                    )
                });

                // Retain the consumed loop measures for subsequent UV rescalings. Keep this
                // marker independent of mUV so enclosing limits still rescale the vacuum mass.
                IntegratedCts::from_series(&integrated, &expansion, 4 * n_loops as i64)
            }
            ApproximationType::Unsubtracted => {
                panic!("should have been kept out of the wood");
            }
        }
    }

    #[debug_instrument(
        current = %current.log_display(),
        given = %given.log_display(),
        reduced,
    )]
    fn integrate<S: ForestNodeLike>(
        &self,
        integrand: &SymbolicTensor<IndexedStructure>,
        sector: &FourDSector,
        ctx: &UVCtx<'_>,
        current: &S,
        given: &S,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        let (protected_integrand, source_ports) = Self::protect_numeric_ports(integrand)?;
        let integrand = &protected_integrand;
        let graph = ctx.graph;
        let reduced = current.reduced_subgraph(given);
        let settings = ctx.settings;
        let reduced_label = reduced.string_label();
        tracing::Span::current().record("reduced", reduced_label.as_str());
        debug_tags!(#uv, #integrated, #vakint, #trace, #input;
            log.integrand = integrand.expression(),
            "Integrating and truncating"
        );
        // Match the vacuum result to the unintegrated CFF loop measure through
        // Vakint's configured per-loop normalization hook. Forest-subtraction
        // signs are folded into the expression at the integrated CT composition
        // sites.
        let integrand_vakint = to_vakint_integrand(
            integrand,
            graph,
            current.subgraph(),
            given.subgraph(),
            &sector.active_components,
            &settings.vakint,
            false,
        )?;

        let propagator = function!(
            vakint::symbols::S.prop,
            W_.a_,
            W_.b_,
            W_.mom_,
            W_.mass_,
            W_.e_
        )
        .to_pattern();
        let retained_loop_count = sector
            .active_components
            .iter()
            .flat_map(|(_, _, lmb)| lmb.loop_edges.iter())
            .collect::<std::collections::HashSet<_>>()
            .len();
        let mut res = Atom::Zero;
        for source_term in integrand_vakint.0 {
            // A product of independent one-loop denominators can carry a
            // coupled numerator. Integrate one factor at a time, treating all
            // spectator loops as fixed vectors for its angular projection.
            let mut loop_factors = std::collections::BTreeMap::<i64, Atom>::new();
            let mut masses = std::collections::HashSet::new();
            let mut independent = true;
            for matched in source_term.integral.pattern_match(&propagator, None, None) {
                let momentum = &matched[&W_.mom_];
                let loop_id = [1, -1].into_iter().find_map(|sign| {
                    let momentum = momentum * Atom::num(sign);
                    let AtomView::Fun(vector) = momentum.as_view() else {
                        return None;
                    };
                    (vector.get_symbol() == vakint::symbols::S.k && vector.get_nargs() == 1)
                        .then(|| i64::try_from(vector.get(0)).ok())
                        .flatten()
                });
                let Some(loop_id) = loop_id else {
                    independent = false;
                    break;
                };
                masses.insert(matched[&W_.mass_].clone());
                *loop_factors.entry(loop_id).or_insert_with(Atom::one) *= propagator
                    .replace_wildcards(&matched)
                    .map_err(|error| eyre!(error))?;
            }
            let factorized = independent
                && retained_loop_count > 1
                && loop_factors.len() == retained_loop_count
                && masses.len() > 1;
            let loop_ids = loop_factors.keys().copied().collect::<Vec<_>>();
            let stages = if factorized {
                loop_factors
                    .into_iter()
                    .map(|(loop_id, propagators)| {
                        (
                            function!(vakint::symbols::S.topo, propagators),
                            loop_ids
                                .iter()
                                .copied()
                                .filter(|id| *id != loop_id)
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>()
            } else {
                vec![(source_term.integral.clone(), Vec::new())]
            };
            let mut coefficient = source_term.numerator.clone();
            if factorized {
                // Preserve even custom normalization hooks that are not a
                // product of identical per-loop factors. The converter's
                // additional normalization already belongs to the full measure.
                let normalization = self
                    .vakint_settings
                    .get_integral_normalization_factor_atom()?;
                coefficient *= normalization
                    .replace(vakint::symbols::S.n_loops.to_pattern())
                    .with(Atom::num(retained_loop_count))
                    / normalization
                        .replace(vakint::symbols::S.n_loops.to_pattern())
                        .with(Atom::one())
                        .pow(retained_loop_count);
            }
            // Keep the same conservative epsilon target at every stage. Later
            // vacuum poles can consume positive orders, while Vakint reserves
            // coefficient and normalization poles independently.
            for (integral, spectators) in stages {
                let stage_loop_count = if factorized { 1 } else { retained_loop_count };
                let external = function!(vakint::symbols::S.p, W_.i_, W_.x___).to_pattern();
                let mut used_ids = coefficient
                    .pattern_match(&external, None, None)
                    .map(|matched| i64::try_from(matched[&W_.i_].as_view()))
                    .collect::<Result<std::collections::HashSet<_>, _>>()
                    .map_err(|error| eyre!(error))?;
                let mut next_id = 0i64;
                let mut restore_spectators = Vec::new();
                for spectator in spectators {
                    while used_ids.contains(&next_id) {
                        next_id = next_id.checked_add(1).ok_or_else(|| {
                            eyre!("no free external momentum ID for vacuum spectator")
                        })?;
                    }
                    used_ids.insert(next_id);
                    let source = function!(vakint::symbols::S.k, spectator, W_.x___);
                    let target = function!(vakint::symbols::S.p, next_id, W_.x___);
                    coefficient = coefficient
                        .replace(source.to_pattern())
                        .with(target.to_pattern());
                    restore_spectators.push(Replacement::new(target, source));
                }
                let mut stage_term = source_term.clone();
                stage_term.integral = integral;
                stage_term.numerator = coefficient;
                let mut integrand_vakint = VakintExpression(vec![stage_term]);
                let mut same_momentum = std::collections::BTreeMap::<Atom, Vec<_>>::new();
                for matched in integrand_vakint.0[0]
                    .integral
                    .pattern_match(&propagator, None, None)
                {
                    let momentum = &matched[&W_.mom_];
                    same_momentum
                        .entry(momentum.clone().min(-momentum))
                        .or_default()
                        .push((
                            usize::try_from(matched[&W_.a_].as_view())
                                .map_err(|error| eyre!(error))?,
                            matched[&W_.mass_].clone(),
                            matched[&W_.e_].clone(),
                        ));
                }
                for (momentum, propagators) in same_momentum {
                    let (_, first_mass, _) = &propagators[0];
                    if propagators.iter().all(|(_, mass, _)| mass == first_mass)
                        || momentum.contains_symbol(vakint::symbols::S.p)
                        || propagators.iter().any(|(_, mass, power)| {
                            mass.contains_symbol(W_.x_)
                                || !i64::try_from(power.as_view()).is_ok_and(|power| power > 0)
                        })
                    {
                        continue;
                    }

                    // A mixed-mass one-loop vacuum is a sum of single-mass
                    // tadpoles: 1/(D_a D_b) = (1/D_a - 1/D_b)/(m_a²-m_b²).
                    // The same identity applies to any equal-momentum group
                    // with other vacuum propagators held fixed. Partial-fraction
                    // only the denominator, including raised powers; the graph
                    // numerator remains an unchanged factor.
                    let loop_square = Atom::var(W_.x_);
                    let denominator = propagators
                        .iter()
                        .fold(Atom::one(), |product, (_, mass, power)| {
                            product * (&loop_square - mass).pow(-power)
                        });
                    let fractions = denominator
                        .try_to_rational_polynomial::<_, _, u16>(&Q, &Z, [W_.x_])?
                        .apart_factored_denominators(0);
                    for term in std::mem::take(&mut integrand_vakint.0) {
                        let powers = term
                            .integral
                            .pattern_match(&propagator, None, None)
                            .map(|matched| {
                                Ok((
                                    usize::try_from(matched[&W_.a_].as_view())
                                        .map_err(|error| eyre!(error))?,
                                    matched[&W_.e_].clone(),
                                ))
                            })
                            .collect::<Result<std::collections::BTreeMap<_, _>>>()?;
                        let graph = vakint::graph::Graph::new_from_atom(
                            term.integral.as_view(),
                            *powers.keys().next_back().unwrap(),
                        )?;
                        for (coefficient, denominator, power) in &fractions {
                            let denominator = denominator.to_expression();
                            let constant = denominator.replace(W_.x_).with(Atom::Zero);
                            let leading = denominator.replace(W_.x_).with(Atom::one()) - &constant;
                            let (survivor, _, _) = propagators.iter().find(|(_, mass, _)| {
                                (&denominator - &leading * (&loop_square - mass)).expand().is_zero()
                            }).ok_or_else(|| eyre!("vacuum partial fraction has an unexpected denominator {denominator}"))?;
                            let coefficient = coefficient.to_expression() / leading.pow(*power);
                            if coefficient.contains_symbol(W_.x_) {
                                return Err(eyre!(
                                    "vacuum partial fraction retained loop dependence in its coefficient"
                                ));
                            }
                            let mut contracted = graph.clone();
                            for (id, _, _) in &propagators {
                                if id != survivor {
                                    contracted.contract_edges(&[*id].into_iter().collect());
                                }
                            }
                            // Removing a serial line contracts its endpoints.
                            // Dropping it without that contraction would silently
                            // lower a sunset's graph rank although both loop
                            // generators still occur in its denominators.
                            if contracted.to_symbolica_graph(false).num_loops() != stage_loop_count
                            {
                                return Err(eyre!(
                                    "vacuum partial fraction changed the retained {stage_loop_count}-loop integration domain for component {reduced_label}"
                                ));
                            }
                            let mut fraction = term.clone();
                            fraction.numerator *= coefficient;
                            fraction.integral = function!(
                                vakint::symbols::S.topo,
                                contracted
                                    .edges
                                    .values()
                                    .fold(Atom::one(), |product, edge| {
                                        product
                                            * function!(
                                                vakint::symbols::S.prop,
                                                edge.id,
                                                function!(
                                                    vakint::symbols::S.edge,
                                                    edge.left_node_id,
                                                    edge.right_node_id
                                                ),
                                                &edge.momentum,
                                                &edge.mass,
                                                if edge.id == *survivor {
                                                    Atom::num(*power)
                                                } else {
                                                    powers[&edge.id].clone()
                                                }
                                            )
                                    })
                            );
                            integrand_vakint.0.push(fraction);
                        }
                    }
                }
                // Only a vacuum topology whose actual masses all vanish is scaleless.
                // Certify polynomial loop dependence as well: an external scale hidden
                // inside a nonpolynomial numerator must not be erased with the soft jet.
                integrand_vakint.0.retain(|term| {
                    let propagators = term
                        .integral
                        .pattern_match(&propagator, None, None)
                        .map(|matched| {
                            matched[&W_.mass_].is_zero()
                                && matched[&W_.mom_].contains_symbol(vakint::symbols::S.k)
                        })
                        .collect::<Vec<_>>();
                    if propagators.is_empty()
                        || propagators.iter().any(|massless_loop| !massless_loop)
                        || term.integral.contains_symbol(vakint::symbols::S.p)
                    {
                        return true;
                    }
                    let mut pending = vec![term.numerator.as_view()];
                    let mut scaleless = true;
                    while let Some(view) = pending.pop() {
                        if !view.contains_symbol(vakint::symbols::S.k) {
                            continue;
                        }
                        match view {
                            AtomView::Add(sum) => pending.extend(sum.iter()),
                            AtomView::Mul(product) => pending.extend(product.iter()),
                            AtomView::Pow(power) => {
                                let (base, exponent) = power.get_base_exp();
                                if i64::try_from(exponent).is_ok_and(|power| power >= 0) {
                                    pending.push(base);
                                } else {
                                    scaleless = false;
                                    break;
                                }
                            }
                            AtomView::Fun(function)
                                if function.get_symbol() == vakint::symbols::S.k
                                    && function.get_nargs() > 0
                                    && usize::try_from(function.get(0)).is_ok() => {}
                            AtomView::Fun(function)
                                if function.get_symbol() == vakint::symbols::S.dot
                                    && function.get_nargs() == 2
                                    && function.iter().all(|argument| {
                                        matches!(argument, AtomView::Fun(vector)
                                            if [vakint::symbols::S.k, vakint::symbols::S.p]
                                                .contains(&vector.get_symbol())
                                                && vector.get_nargs() == 1
                                                && usize::try_from(vector.get(0)).is_ok())
                                    }) => {}
                            _ => {
                                scaleless = false;
                                break;
                            }
                        }
                    }
                    if scaleless {
                        debug_tags!(#uv, #integrated, #vakint, #trace;
                            stage = "certified_scaleless_vacuum",
                            log.integral = term.integral,
                            log.numerator = term.numerator,
                            "Dropping a massless vacuum term with polynomial loop numerator"
                        );
                    }
                    !scaleless
                });
                if integrand_vakint.0.is_empty() {
                    coefficient = Atom::Zero;
                    break;
                }

                for (term_index, t) in integrand_vakint.0.iter().enumerate() {
                    debug_tags!(#uv,#integrated,#vakint,#trace,#to_vakint;
                        term_index = %term_index,
                        log.integral = t.integral,
                        log.numerator = t.numerator,
                        "Vakint term as input"
                    );
                }
                debug_tags!(#uv,#integrated,#vakint;settings = ?&self.vakint_settings,"Vakint args");

                // let mut res = vakint
                //     .0
                //     .evaluate(&vakint.1, integrand_vakint.as_view())
                //     .unwrap();

                // The analytic backends identify a vacuum mass through a
                // symbol or its square. Alias other exact nonzero mass values
                // only in propagator mass slots, then restore the values below.
                let mut mass_aliases = std::collections::BTreeMap::new();
                let mut used_symbols = integrand_vakint
                    .0
                    .iter()
                    .flat_map(|term| [&term.integral, &term.numerator])
                    .flat_map(|atom| atom.get_all_symbols(true))
                    .collect::<std::collections::HashSet<_>>();
                let mut mass_serial = 0;
                for term in &mut integrand_vakint.0 {
                    let mut replacements = Vec::new();
                    for matched in term.integral.pattern_match(&propagator, None, None) {
                        let mass = &matched[&W_.mass_];
                        if mass.contains_symbol(vakint::symbols::S.k)
                            || mass.contains_symbol(vakint::symbols::S.p)
                        {
                            return Err(eyre!(
                                "a vacuum propagator mass must be independent of momentum"
                            ));
                        }
                        if mass.is_zero()
                            || matches!(mass.as_view(), AtomView::Var(_))
                            || matches!(mass.as_view(), AtomView::Pow(power)
                                if matches!(power.get_base_exp().0, AtomView::Var(_))
                                    && power.get_base_exp().1 == Atom::num(2).as_view())
                        {
                            continue;
                        }
                        let alias = mass_aliases.entry(mass.clone()).or_insert_with(|| {
                            loop {
                                let symbol = vakint_symbol!(format!(
                                    "integrated_mass_squared_{mass_serial}"
                                ));
                                mass_serial += 1;
                                if used_symbols.insert(symbol) {
                                    break Atom::var(symbol);
                                }
                            }
                        });
                        replacements.push(Replacement::new(
                            propagator
                                .replace_wildcards(&matched)
                                .map_err(|error| eyre!(error))?,
                            function!(
                                vakint::symbols::S.prop,
                                &matched[&W_.a_],
                                &matched[&W_.b_],
                                &matched[&W_.mom_],
                                &*alias,
                                &matched[&W_.e_]
                            ),
                        ));
                    }
                    term.integral = term.integral.replace_multiple(&replacements);
                }
                integrand_vakint.canonicalize(
                    self.vakint_settings,
                    &self.vakint.topologies,
                    false,
                )?;
                for (term_index, t) in integrand_vakint.0.iter().enumerate() {
                    debug_tags!(#uv,#integrated,#vakint,#trace,#canonicalize;
                        term_index = %term_index,
                        log.integral = t.integral,
                        log.numerator = t.numerator,
                        "Vakint term after canonicalization"
                    );
                }
                integrand_vakint.tensor_reduce(self.vakint, self.vakint_settings)?;
                for (term_index, t) in integrand_vakint.0.iter().enumerate() {
                    debug_tags!(#uv,#integrated,#vakint,#trace,#tensor_reduce;
                        term_index = %term_index,
                        log.integral = t.integral,
                        log.numerator = t.numerator,
                        "Vakint term after tensor reduction"
                    );
                }
                for (term_index, term) in integrand_vakint.0.iter_mut().enumerate() {
                    term.numerator =
                        self.simplify_projected_numerator(&term.numerator, current.topo_order())?;
                    debug_tags!(#uv, #integrated, #vakint, #trace;
                        stage = "projected_numerator_after_d_dimensional_algebra",
                        term_index = %term_index,
                        log.numerator = term.numerator,
                        "Completed projected numerator algebra before Laurent expansion"
                    );
                }
                let started = std::time::Instant::now();
                debug_tags!(#uv, #integrated, #diagnostic_stage;
                    stage = 2u8, finished = false,
                    term_count = integrand_vakint.0.len(),
                    "Start scalar vacuum integral evaluation"
                );
                integrand_vakint.evaluate_integral(self.vakint, self.vakint_settings)?;
                debug_tags!(#uv, #integrated, #diagnostic_stage;
                    stage = 2u8, finished = true,
                    elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                    output_bytes = integrand_vakint.0.iter()
                        .map(|term| term.numerator.as_view().get_byte_size()).sum::<usize>(),
                    "Completed scalar vacuum integral evaluation"
                );
                for (term_index, t) in integrand_vakint.0.iter().enumerate() {
                    debug_tags!(#uv,#integrated,#vakint,#trace,#evaluate;
                        term_index = %term_index,
                        log.integral = t.integral,
                        log.numerator = t.numerator,
                        "Vakint term after evaluation"
                    );
                }

                coefficient = Atom::from(integrand_vakint).replace_multiple(
                    mass_aliases
                        .into_iter()
                        .map(|(mass, alias)| Replacement::new(alias, mass))
                        .collect::<Vec<_>>(),
                );
                if coefficient.contains_symbol(vakint::symbols::S.k) {
                    return Err(eyre!(
                        "integrated vacuum factor retained an active loop momentum"
                    ));
                }
                coefficient = coefficient.replace_multiple(&restore_spectators);
            }
            res += coefficient;
        }

        debug_tags!(#uv,#integrated,#vakint,#trace,#raw;
            log.res = res,
            "Raw post vakint "
        );

        let started = std::time::Instant::now();
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 3u8, finished = false,
            input_bytes = res.as_view().get_byte_size(),
            "Start integrated coefficient cleanup"
        );
        let res = Self::restore_numerator(res, current.topo_order());

        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // Close internal momentum pairs before independently integrated factors
        // are multiplied. Their local UV dummy labels can coincide; completed
        // scalar products must not expose those labels to the next factor.
        let dotted = integrand
            .with_checked_expression(cooking.try_cook(res.as_view())?)?
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })?
            .to_dots()?;
        let regulated =
            dotted.with_algebra_result(Self::dimensionally_regularized(dotted.expression()))?;

        debug_tags!(#uv, #integrated, #vakint, #inspect, #trace, #replace;
            log.res = regulated.expression(),
            "Replaced post vakint "
        );

        // This strips as many dummies as possible after undoing chains and traces,
        // so that terms can merge later on.
        let collected = regulated.contract(ContractSettings {
            representations: Some(&[Bispinor {}.into()]),
            metrics: false,
            rank_one: false,
            collect_chains: true,
            collect_traces: true,
            ..Default::default()
        })?;
        // Collection already returns isolated one-factor words in indexed form.
        debug_tags!(#uv, #integrated, #vakint, #profile, #trace, #collect;
            log.expr = collected.expression(),
            "Integrated UV cleanup through shared chain collection"
        );

        // println!("\nIntegrated CT:\n{}\n", res);
        // Normalize the fully integrated UV subgraph before Laurent projection
        // and reinsertion: identical scalar masters can arrive as 2*(A+B) or
        // 2*A+2*B. No cograph numerator is present at this analytic boundary.
        let result = collected.expanded(None, false)?;
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 3u8, finished = true,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            output_bytes = result.expression().as_view().get_byte_size(),
            "Completed integrated coefficient cleanup"
        );
        if source_ports.is_empty() {
            Ok(result)
        } else {
            Ok(result.reindex_interface_ports(&source_ports)?)
        }
    }

    fn simplify_projected_numerator(
        &self,
        numerator: &Atom,
        topology_order: usize,
    ) -> Result<Atom> {
        let started = std::time::Instant::now();
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 4u8, finished = false,
            topology_order,
            input_bytes = numerator.as_view().get_byte_size(),
            "Start projected numerator cleanup"
        );
        // Tensor projection can identify Lorentz slots belonging to different
        // gamma matrices. Complete their d-dimensional algebra before Vakint
        // expands scalar coefficients in epsilon; otherwise an evanescent
        // contraction can multiply a pole only after its finite term was lost.
        // The projected expression contains native gamma/fixed-vector leaves
        // alongside Vakint momenta and metrics. Restore their shared typed
        // interface before contraction; the scalar-only Vakint dot parser
        // cannot account for ports carried by those native tensor leaves.
        let numerator = Self::restore_numerator(numerator.clone(), topology_order);
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let numerator = SymbolicTensor::infer_indexed(cooking.try_cook(numerator.as_view())?)?;
        let numerator = simplify(&numerator)?;
        Self::ensure_resolved_lorentz_contractions(numerator.expression())?;
        let numerator = numerator
            .with_algebra_result(Self::dimensionally_regularized(numerator.expression()))?;
        let materialized = numerator
            .expression()
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: false,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: true,
                    chain: true,
                },
                ..Default::default()
            })?
            .simple_execute::<()>()?;
        let numerator = numerator.with_checked_expression(materialized)?;
        let result = Self::to_vakint_numerator(&numerator)?;
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 4u8, finished = true,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            output_bytes = result.as_view().get_byte_size(),
            "Completed projected numerator cleanup"
        );
        Ok(result)
    }

    /// Inspect a borrowed native expression without changing its retained owner.
    fn ensure_resolved_lorentz_contractions(numerator: &Atom) -> Result<()> {
        let minkowski = Minkowski {}.to_symbolic([GS.dim]).get_symbol().unwrap();
        if !numerator.contains_symbol(minkowski) {
            return Ok(());
        }
        let started = std::time::Instant::now();
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 5u8, finished = false,
            input_bytes = numerator.as_view().get_byte_size(),
            "Start structural residual-contraction check"
        );
        // Expose compact inter-chain contractions, retaining every analytic
        // branch and scalar coefficient. The shared occurrence observer checks
        // branch-local pairs without distributing products of sums.
        // Native index allocation reserves the existing encoded identities.
        // Neither source indices nor newly materialized dummy names cross an
        // encoding boundary, so distinct labels stay distinct during the check.
        let explicit = numerator
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode::full(),
                    trace: false,
                    chain: true,
                },
                ..Default::default()
            })
            .map_err(|error| eyre!("invalid analytic UV Lorentz tensor notation: {error}"))?
            .simple_execute::<()>()?;
        let mut slots = SlotMatcher::default();
        let mut unresolved_inner_product = false;
        explicit.visitor(&mut |view| {
            let AtomView::Fun(function) = view else {
                return true;
            };
            if [SPENSO_TAG.dot, ETS.metric].contains(&function.get_symbol())
                && view.contains_symbol(minkowski)
                && !function
                    .iter()
                    .all(|argument| slots.parse::<Minkowski, AbstractIndex>(argument).is_ok())
            {
                // Full shorthand parsing must expose compact Lorentz products
                // before excluded momentum leaves can hide their ports.
                unresolved_inner_product = true;
            }
            // Function payloads remain opaque, as in the tensor leaf observer.
            false
        });
        if unresolved_inner_product {
            return Err(eyre!(
                "an unresolved compact Lorentz inner product remains after d-dimensional UV numerator algebra"
            ));
        }
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "residual_lorentz_contraction_materialization",
            file.source = %numerator.to_plain_string(),
            file.explicit = %explicit.to_plain_string(),
            "Check residual contractions in the native index scope"
        );
        let momenta = [
            GS.emr_mom,
            GS.loop_mom,
            vakint::symbols::S.p,
            vakint::symbols::S.k,
        ];
        let mut selectors = SlotMatcher::default();
        if InterfaceInference::has_explicit_index_pairs(
            explicit.as_view(),
            Minkowski {}.into(),
            |value| {
                let AtomView::Fun(leaf) = value else {
                    return false;
                };
                if momenta.contains(&leaf.get_symbol()) {
                    return true;
                }
                if leaf.get_symbol() != ETS.delta || leaf.get_nargs() != 2 {
                    return false;
                }
                let SlotMatch::Explicit(slot) = selectors.classify(leaf.get(1)) else {
                    return false;
                };
                let representation = slot.representation();
                // A validated fixed physical basis projection is not a free
                // D-dimensional tensor contraction. Keep D and scalar epsilon
                // dependence untouched; unchecked delta shapes remain visible.
                representation.head() == minkowski
                    && representation.is_self_dual()
                    && selectors
                        .fixed_basis_component(leaf.get(0), representation)
                        .is_some_and(|component| component < 4)
            },
        )
        .map_err(|error| eyre!("invalid residual d-dimensional tensor structure: {error}"))?
        {
            return Err(eyre!(
                "an unresolved internal Lorentz tensor contraction remains after d-dimensional UV numerator algebra; a d-dimensional tensor-operator reduction is required before analytic integration"
            ));
        }
        debug_tags!(#uv, #integrated, #diagnostic_stage;
            stage = 5u8, finished = true,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            output_bytes = explicit.as_view().get_byte_size(),
            "Completed structural residual-contraction check"
        );
        Ok(())
    }

    fn dimensionally_regularized(numerator: &Atom) -> Atom {
        let minkowski = Minkowski {}.to_symbolic([GS.dim]).get_symbol().unwrap();
        let dimension = Atom::var(GS.dim);
        let regulated = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        numerator.replace_map(|view, _, out| {
            if view.get_symbol() == Some(minkowski) {
                // A retained open Lorentz basis keeps its representation. All
                // scalar occurrences, including inside function coefficients,
                // must instead participate in the Laurent expansion.
                **out = view.to_owned();
            } else if view == dimension.as_view() {
                **out = regulated.clone();
            }
        })
    }

    /// Restore the raw Vakint tensor ABI. Public numeric source ports must be
    /// transported to distinct native labels before entering the backend;
    /// unlabelled numeric indices here belong to its topology-local dummies.
    pub(crate) fn restore_numerator(mut res: Atom, topology_order: usize) -> Atom {
        // Projector indices may close two opaque tensor leaves. Restore their
        // complete Minkowski representation before exposing the tensor bodies.
        let mink = Minkowski {}.new_rep(GS.dim);
        let minkowski = mink.to_symbolic([]).get_symbol().unwrap();
        let regulated_dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        // The projector works in 4-2*epsilon while opaque source slots retain
        // d. These are the same Lorentz representation; unify only metadata,
        // leaving every scalar evanescent coefficient for the later expansion.
        res = res
            .replace(function!(minkowski, &regulated_dimension, W_.x___))
            .with(function!(minkowski, GS.dim, W_.x___))
            .replace(function!(
                minkowski,
                W_.d_,
                vakint::symbols::S.dot_dummy_ind(W_.j_)
            ))
            .when(W_.j_.filter(|r| r.is_integer().is_true()))
            .with(function!(
                minkowski,
                W_.d_,
                function!(GS.uvaind, topology_order, W_.j_)
            ));

        let mut projected_indices = std::collections::BTreeMap::new();
        res = res.replace_map(|view, _, out| {
            if matches!(view, AtomView::Fun(index) if index.get_symbol() == vakint::symbols::S.tensor_index) {
                let replacement = projected_indices.entry(view.to_owned()).or_insert_with(|| loop {
                    let candidate = mink.to_symbolic([Aind::new_dummy().to_atom()]);
                    if res.pattern_match(&candidate.to_pattern(), None, None).next().is_none() {
                        break candidate;
                    }
                });
                **out = replacement.clone();
            }
        });

        res = res
            .replace(parse_lit!(vakint::cl2))
            .with(parse_lit!(cl2))
            .replace(parse_lit!(vakint::sqrt3))
            .with(parse_lit!(sqrt(3)));

        let vk_metric = vakint_symbol!("g");

        // Vakint momenta omit their compact ports, whereas fixed native
        // vectors retain them. Restore both operands of the same dot without
        // removing or reconstructing the fixed vector's actual tensor leaf.
        res = res.replace_map(|view, _, out| {
            let AtomView::Fun(dot) = view else { return };
            if dot.get_symbol() != vakint::symbols::S.dot || dot.get_nargs() != 2 {
                return;
            }
            let restored = [dot.get(0), dot.get(1)].map(|vector| {
                if let AtomView::Fun(momentum) = vector
                    && [vakint::symbols::S.k, vakint::symbols::S.p].contains(&momentum.get_symbol())
                    && momentum.get_nargs() == 1
                {
                    let head = if momentum.get_symbol() == vakint::symbols::S.k {
                        GS.loop_mom
                    } else {
                        GS.emr_mom
                    };
                    function!(head, momentum.get(0), mink.to_symbolic([]))
                } else {
                    vector.to_owned()
                }
            });
            **out = function!(SPENSO_TAG.dot, &restored[0], &restored[1]);
        });

        // Raw indexed momenta also occur when dot notation is disabled. Bind
        // loop and external indices in the same topology-local chart before
        // the native tensor owner contracts them with projected metrics.
        for (source, target) in [
            (vakint::symbols::S.k, GS.loop_mom),
            (vakint::symbols::S.p, GS.emr_mom),
        ] {
            res = res
                .replace(function!(source, W_.i_, W_.j_))
                .when(W_.j_.filter(|r| r.is_integer().is_true()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([function!(GS.uvaind, topology_order, W_.j_)])
                ))
                .replace(function!(
                    source,
                    W_.i_,
                    vakint::symbols::S.dot_dummy_ind(W_.j_)
                ))
                .when(W_.j_.filter(|r| r.is_integer().is_true()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([function!(GS.uvaind, topology_order, W_.j_)])
                ))
                // Native projection can return the original abstract label
                // without its representation. Preserve that source identity;
                // only raw numeric backend dummies use the topology-local chart.
                .replace(function!(source, W_.i_, W_.j_))
                .when(W_.j_.filter(|index| Aind::try_from(index).is_ok()))
                .with(function!(
                    source,
                    W_.i_,
                    mink.to_symbolic([Atom::var(W_.j_)])
                ))
                .replace(function!(source, W_.x___))
                .with(function!(target, W_.x___));
        }
        res = res.replace_map(|view, _, out| {
            let AtomView::Fun(metric) = view else { return };
            if metric.get_symbol() != vk_metric || metric.get_nargs() != 2 {
                return;
            }
            let slots = [metric.get(0), metric.get(1)].map(|slot| {
                let raw_index = if let AtomView::Fun(index) = slot
                    && index.get_symbol() == vakint::symbols::S.dot_dummy_ind
                    && index.get_nargs() == 1
                {
                    index.get(0)
                } else {
                    slot
                };
                if raw_index.is_integer().is_true() {
                    mink.to_symbolic([function!(GS.uvaind, topology_order, raw_index)])
                } else if Aind::try_from(raw_index).is_ok() {
                    mink.to_symbolic([raw_index.to_owned()])
                } else {
                    slot.to_owned()
                }
            });
            **out = function!(ETS.metric, &slots[0], &slots[1]);
        });

        res = res.replace(vakint::symbols::S.cmplx_i).with(Atom::i());

        res
    }

    fn to_vakint_numerator(numerator: &SymbolicTensor<IndexedStructure>) -> Result<Atom> {
        // A concrete component is a fixed-vector projection, not an Einstein
        // index. Keep repeated temporal/spatial components distinct from dummy
        // multiplicity before Vakint supplies its formal-dimensional slots.
        let representation = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let slots = SlotMatcher::default();
        let mut component_error = None;
        let expression = numerator.expression().replace_map(|view, _, output| {
            let AtomView::Fun(momentum) = view else {
                return;
            };
            if ![GS.loop_mom, GS.emr_mom].contains(&momentum.get_symbol())
                || momentum.get_nargs() != 2
            {
                return;
            }
            let Some(component) = slots.concrete_component(momentum.get(1)) else {
                return;
            };
            let component = match component {
                Ok(component) if component < 4 => component,
                other => {
                    component_error = Some(eyre!(
                        "invalid physical Minkowski momentum component: {other:?}"
                    ));
                    return;
                }
            };
            let basis = function!(GS.delta_vec, GS.cind(component), &representation);
            let vector = function!(momentum.get_symbol(), momentum.get(0), &representation);
            let sign = if (Minkowski {}).is_neg(component) {
                -1
            } else {
                1
            };
            **output = Atom::num(sign) * function!(SPENSO_TAG.dot, basis, vector);
        });
        if let Some(error) = component_error {
            return Err(error);
        }
        // Normalize scalar coefficients inside compact dots before translating
        // their vector leaves. The shared notation owner preserves powers and
        // indexed contractions without opening the surrounding numerator.
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let dotted = numerator.with_checked_expression(expression)?.to_dots()?;
        let numerator = cooking.uncook(dotted.expression().as_view());
        // Vakint's existing bridge carries opaque tensor slots through its
        // Lorentz projectors. Keep the actual leaves: wrapping a duplicate body
        // and slot list lets contraction relabel only the latter, disconnecting
        // the restored gamma matrices from their projected metric.
        Ok(numerator
            .replace(function!(GS.loop_mom, W_.x___))
            .with(function!(vakint::symbols::S.k, W_.x___))
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___))
            .replace_map(|view, _, out| {
                let AtomView::Fun(dot) = view else { return };
                if dot.get_symbol() != SPENSO_TAG.dot || dot.get_nargs() != 2 {
                    return;
                }
                let vectors = [dot.get(0), dot.get(1)].map(|vector| {
                    if let AtomView::Fun(momentum) = vector
                        && [vakint::symbols::S.k, vakint::symbols::S.p]
                            .contains(&momentum.get_symbol())
                        && momentum.get_nargs() == 2
                        && momentum.get(1)
                            == (Minkowski {}).new_rep(GS.dim).to_symbolic([]).as_view()
                    {
                        function!(momentum.get_symbol(), momentum.get(0))
                    } else {
                        // Only Vakint momenta have an implicit representation.
                        // A fixed vector must retain its native compact port so
                        // angular projection restores the actual tensor leaf.
                        vector.to_owned()
                    }
                });
                **out = vakint::symbols::S.dot(&vectors[0], &vectors[1]);
            })
            .replace(function!(
                ETS.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            ))
            .with(function!(
                vakint::symbols::S.metric,
                Minkowski {}.to_symbolic([W_.a__]),
                Minkowski {}.to_symbolic([W_.b__])
            )))
    }
}

#[debug_instrument]
pub(crate) fn to_vakint_integrand<
    E: UVE,
    V,
    H,
    S: SubGraphLike + SubSetLike<Base = SuBitGraph>,
    SS: SubGraphLike,
>(
    integrand: &SymbolicTensor<IndexedStructure>,
    graph: &HedgeGraph<E, V, H>,
    reduced: &S,
    dependent_subgraph: &SS,
    active_components: &[(SuBitGraph, SuBitGraph, LoopMomentumBasis)],
    settings: &VakintSettings,
    substitute_masses_to_m_uv: bool,
) -> Result<VakintExpression> {
    let additional_normalization =
        SymbolicTensor::infer_indexed(parse!(&settings.additional_normalization))?;
    if !additional_normalization.is_scalar() {
        return Err(eyre!(
            "additional_normalization must have a scalar tensor interface"
        ));
    }
    let reduced_label = reduced.string_label();
    let dependent_subgraph_label = dependent_subgraph.string_label();
    let source_lmbs = active_components
        .iter()
        .map(|(_, _, lmb)| lmb)
        .collect::<Vec<_>>();
    let mut retained_loop_edges = source_lmbs
        .iter()
        .flat_map(|lmb| lmb.loop_edges.iter().copied())
        .collect::<Vec<_>>();
    retained_loop_edges.sort();
    retained_loop_edges.dedup();
    // Resolve color contractions before opening shorthand traces. Repeated
    // color spectators then share their reduced form instead of acquiring
    // distinct dummy indices during Lorentz preparation. Keep the kinematic
    // coefficients factorized throughout this independent color pass.
    let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
    let input = GS.erase_uv_momentum_provenance(integrand.expression());
    let color = integrand.with_checked_expression(input)?.simplify_algebra(
        &idenso::tensor::AlgebraSettings {
            color: Some(ColorSimplifySettings {
                ..Default::default()
            }),
            ..Default::default()
        },
    )?;
    let lorentz_input = cooking.uncook(color.expression().as_view());
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_before_lorentz_materialization",
        reduced = %reduced_label,
        file.expr = %lorentz_input.to_plain_string(),
        "Admitted color result before Lorentz network materialization"
    );
    let mut integrand_vakint = lorentz_input
        .parse_to_symbolic_net::<Aind>(&ParseSettings {
            shorthand_parsing: ShorthandParsing::Expand {
                schoonschip: SchoonschipExpansionMode {
                    inner_products: false,
                    expand_inside_chains: true,
                    expand_schoonship: true,
                },
                trace: true,
                chain: true,
            },
            ..Default::default()
        })?
        .simple_execute::<()>()?;
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_undo_shorthands",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        substitute_masses_to_m_uv = substitute_masses_to_m_uv,
        log.integrand = integrand_vakint,
        "Vakint trace after undo shorthands"
    );
    //Atom::Zero

    // The denominator-to-propagator replacements below strip the momentum
    // wrapper without distributing the numerator. The former standalone
    // denominator rewrite was:
    // .replace(function!(
    //     GS.den,
    //     W_.prop_,
    //     function!(GS.emr_mom, W_.prop_, W_.mom_),
    //     W_.x__
    // ))
    // .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__))

    // Nested counterterms can expose a boundary metric next to the propagator
    // metric of the reduced graph. Contract those metric-only structures before
    // the expression is split into Vakint terms.
    let metric_contracted = color
        .with_checked_expression(cooking.try_cook(integrand_vakint.as_view())?)?
        .contract(ContractSettings {
            rank_one: false,
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })?;
    integrand_vakint = cooking.uncook(metric_contracted.expression().as_view());
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_simplify_metrics",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after metric simplification"
    );

    // Separate denominator monomials before contracting their graph incidence.
    // In A*B*(c*B+d*B^2), fusing only the outside A*B would leave the
    // inner B attached to nodes that have already been contracted. Collect
    // only complete denominator atoms; numerator sums stay factored coefficients.
    let denominator_pattern = function!(GS.den, W_.x___).to_pattern();
    let denominators = integrand_vakint
        .pattern_match(&denominator_pattern, None, None)
        .map(|matched| denominator_pattern.replace_wildcards(&matched).unwrap())
        .collect::<std::collections::HashSet<_>>();
    integrand_vakint = integrand_vakint
        .coefficient_list::<i32>(&denominators.into_iter().collect::<Vec<_>>())
        .into_iter()
        .map(|(denominator, numerator)| denominator * numerator)
        .sum();
    let mut propagator_id = 1;
    let mut propagator_replacements = Vec::new();
    let mut incidence_replacements = Vec::new();
    let mut assigned_edges = std::collections::HashSet::new();

    let vk_prop = vakint::symbols::S.prop;
    let vk_edge = vakint_symbol!("edge");
    let vk_topo = vakint_symbol!("topo");

    // let contracted_nodes: BTreeSet<NodeIndex> = graph
    //     .iter_nodes_of(dependent_subgraph)
    //     .map(|(nid, c, v)| nid)
    //     .collect();

    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(reduced),
        "Den to prop for"
    );
    // let first = contracted_nodes.first();

    for (owners, source_scope, _) in active_components {
        // Shrink only this component's omitted prefix. Active child propagators
        // retain their own incidence; its quotient attaches at the contracted
        // vertex. Contract disjoint prefixes separately, preserving their loop
        // domains instead of identifying unrelated component vertices.
        let mut contracted_nodes = std::collections::HashMap::new();
        for prefix in graph.connected_components(&source_scope.subtract(owners)) {
            let mut nodes = graph.iter_nodes_of(&prefix).map(|(node, _, _)| node);
            if let Some(first) = nodes.next() {
                for node in nodes {
                    contracted_nodes.insert(node, first);
                }
            }
        }
        for (pair, index, _data) in graph.iter_edges_of(owners) {
            if let HedgePair::Paired { source, sink } = pair {
                if !assigned_edges.insert(index) {
                    return Err(eyre!(
                        "Vakint denominator owner {index} belongs to multiple active components"
                    ));
                }
                // The former global contraction selected `first` from all
                // dependent nodes. Resolve each endpoint with this owner's
                // prefix instead, leaving active child endpoints unchanged.
                let source = graph.node_id(source);
                let sink = graph.node_id(sink);
                let contracted_source = contracted_nodes.get(&source).copied().unwrap_or(source);
                let contracted_sink = contracted_nodes.get(&sink).copied().unwrap_or(sink);
                propagator_replacements.push(Replacement::new(
                    function!(
                        GS.den,
                        usize::from(index) as i64,
                        W_.mom_,
                        W_.mass_,
                        W_.x___
                    ),
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(vk_edge, usize::from(source), usize::from(sink)),
                        W_.mom_,
                        if substitute_masses_to_m_uv {
                            Atom::var(GS.m_uv_vacuum).pow(2)
                        } else {
                            Atom::var(W_.mass_)
                        },
                        1
                    ),
                ));
                incidence_replacements.push(Replacement::new(
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(vk_edge, usize::from(source), usize::from(sink)),
                        W_.x___
                    ),
                    function!(
                        vk_prop,
                        propagator_id,
                        function!(
                            vk_edge,
                            usize::from(contracted_source),
                            usize::from(contracted_sink)
                        ),
                        W_.x___
                    ),
                ));
                propagator_id += 1;
            }
        }
    }
    // Edge IDs make these denominator replacements disjoint. Convert all
    // propagators together, then encode their powers once for the whole atom.
    integrand_vakint = integrand_vakint
        .replace_multiple(&propagator_replacements)
        .replace(function!(vk_prop, W_.x___, 1).pow(Atom::var(W_.e_)))
        .with(function!(vk_prop, W_.x___, -Atom::var(W_.e_)));
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_den_to_prop",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after denominator-to-propagator conversion"
    );

    debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
        reduced = %graph.dot(dependent_subgraph),
        "Shrinking each component's omitted prefix for vakint"
    );
    // Shrink vertices only on propagators belonging to that quotient.
    integrand_vakint = integrand_vakint.replace_multiple(&incidence_replacements);
    if integrand_vakint.contains_symbol(GS.den) {
        return Err(eyre!(
            "Vakint component {reduced_label} contains a denominator without an active owner"
        ));
    }
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_shrink_subgraph",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after shrinking subgraph"
    );

    // flip edges to positive momentum
    // FIXME: how will this work for sums of momenta?
    integrand_vakint = integrand_vakint
        .replace(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.a_, W_.b_),
            -Atom::var(W_.y_),
            W_.e___
        ))
        .repeat()
        .with(function!(
            vk_prop,
            W_.x_,
            function!(vk_edge, W_.b_, W_.a_),
            W_.y_,
            W_.e___
        ));

    // fuse raised edges
    integrand_vakint = integrand_vakint
        .replace(
            function!(
                vk_prop,
                W_.x_,
                function!(vk_edge, W_.a_, W_.b_),
                W_.x___,
                W_.e_
            ) * function!(
                vk_prop,
                W_.y_,
                function!(vk_edge, W_.b_, W_.c_),
                W_.x___,
                W_.f_
            ),
        )
        .repeat()
        .with(function!(
            vk_prop,
            W_.y_,
            function!(vk_edge, W_.a_, W_.c_),
            W_.x___,
            W_.e_ + W_.f_
        ));
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_after_flip_fuse",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = integrand_vakint,
        "Vakint trace after edge flip and fuse"
    );

    // println!(
    //     "Integrand pre vakint: {:}",
    //     VakintExpression::try_from(
    //         integrand_vakint
    //             .replace(function!(vk_prop, W_.x__))
    //             .with(function!(vk_topo, function!(vk_prop, W_.x__)))
    //             .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
    //             .repeat()
    //             .with(function!(vk_topo, W_.x_ * W_.y_))
    //     )
    //     .unwrap()
    // );

    let vakint_input_atom = integrand_vakint
        .replace(function!(vk_prop, W_.x__))
        .with(function!(vk_topo, function!(vk_prop, W_.x__)))
        .replace(function!(vk_topo, W_.x_) * function!(vk_topo, W_.y_))
        .repeat()
        .with(function!(vk_topo, W_.x_ * W_.y_));
    debug_tags!(#uv, #integrated, #vakint, #trace;
        stage = "to_vakint_integrand_before_split_terms",
        reduced = %reduced_label,
        dependent_subgraph = %dependent_subgraph_label,
        log.integrand = vakint_input_atom,
        "Vakint trace before split terms"
    );

    let mut a = VakintExpression::try_from(vakint_input_atom)
        .wrap_err("could not split integrand into Vakint terms")?;

    for (term_index, t) in a.0.iter_mut().enumerate() {
        debug_tags!(#uv, #integrated, #vakint, #inspect, #trace;
            stage = "to_vakint_integrand_term_initial",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Starting integral"
        );

        let mut graph = HedgeGraphBuilder::new();
        //prop(<id>, edge(<node_source>,<node_sink>), <mom>, <mass_squared>, <power>)
        let pat = function!(
            vk_prop,
            W_.a_,
            function!(vakint::symbols::S.edge, W_.i_, W_.j_),
            W_.c_,
            W_.d_,
            W_.e_
        )
        .to_pattern();
        let mut nodemap = std::collections::HashMap::new();

        struct ContractibleEdge {
            mom: Atom,
            mass_squared: Atom,
            power: i32,
        }

        for m in t.integral.pattern_match(&pat, None, None) {
            let i: usize = m[&W_.i_].as_view().try_into().unwrap();
            let j: usize = m[&W_.j_].as_view().try_into().unwrap();
            nodemap.entry(i).or_insert_with(|| graph.add_node(()));
            nodemap.entry(j).or_insert_with(|| graph.add_node(()));

            graph.add_edge(
                nodemap[&i],
                nodemap[&j],
                ContractibleEdge {
                    mom: m[&W_.c_].clone(),
                    mass_squared: m[&W_.d_].clone(),
                    power: m[&W_.e_].as_view().try_into().unwrap(),
                },
                false,
            );
        }

        let mut system = vec![];
        let mut momentum_variables = vec![];

        let mut graph: HedgeGraph<ContractibleEdge, ()> = graph.build();
        let uncontracted_propagator_count =
            graph.iter_edges().filter(|(p, _, _)| p.is_paired()).count();
        let uncontracted_propagator_power_sum = graph
            .iter_edges()
            .filter_map(|(p, _, e)| p.is_paired().then_some(e.data.power))
            .sum::<i32>();

        while let Some(same_mass_two_bond) = graph.a_bond(&|c| {
            let mut count = 0;
            let mut mass_squared = None;
            for (_, _, d) in graph.iter_edges_of(c) {
                count += 1;

                if let Some(m) = &mass_squared
                    && m != &d.data.mass_squared
                {
                    return false;
                } else {
                    mass_squared = Some(d.data.mass_squared.clone());
                }
                if count > 2 {
                    return false;
                }
            }
            count == 2
        }) {
            let mut iter = same_mass_two_bond.included_iter();
            let first = graph[&iter.next().unwrap()];
            let second = graph[&iter.next().unwrap()];
            graph[first].power += graph[second].power;
            let mut to_contract: SuBitGraph = graph.empty_subgraph();
            to_contract.add(graph[&second].1);
            graph.contract_subgraph(&to_contract, ());
        }
        let mut nodes_to_merge = vec![];

        for c in graph.connected_components(&graph.full_filter()) {
            let Some((nid, _, _)) = graph.iter_nodes_of(&c).next() else {
                continue;
            };
            nodes_to_merge.push(nid);
        }

        if !nodes_to_merge.is_empty() {
            graph.identify_nodes(&nodes_to_merge, ());
        }

        graph.forget_identification_history();
        debug_tags!(#uv, #integrated, #vakint, #graph, #dump;
            log.graph = %graph.base_dot(),
            "Graph"
        );

        let mut new_integral: Atom = 1.into();
        for (p, eid, e) in graph.iter_edges() {
            let HedgePair::Paired { source, sink } = p else {
                continue;
            };
            new_integral *= function!(
                vk_prop,
                eid.0 + 1, //vakint propagator ids are 1-indexed
                function!(
                    vakint::symbols::S.edge,
                    graph.node_id(source).0,
                    graph.node_id(sink).0
                ),
                &e.data.mom,
                &e.data
                    .mass_squared
                    .replace(GS.m_uv_expansion)
                    .with(GS.m_uv_vacuum),
                e.data.power
            )
        }

        // println!("{}->{}", t.integral, new_integral);
        t.integral = function!(vakint::symbols::S.topo, new_integral);
        let nloops = graph.cyclotomatic_number(&graph.full_filter());
        // A cancelled denominator can remove an integration direction without
        // leaving that coordinate in the numerator. Never silently replace the
        // retained active measure by a lower-loop vacuum integral. Redundant
        // EMR variables remain allowed when they span the full retained domain.
        if nloops != retained_loop_edges.len() {
            return Err(eyre!(
                "Vakint term {term_index} for component {reduced_label} after prefix {dependent_subgraph_label} has {nloops} loop generators, expected {} retained active generators {:?}; a changed integration domain requires an explicit scaleless or factorization certificate",
                retained_loop_edges.len(),
                retained_loop_edges,
            ));
        }
        let contracted_propagator_count =
            graph.iter_edges().filter(|(p, _, _)| p.is_paired()).count();
        let contracted_propagator_power_sum = graph
            .iter_edges()
            .filter_map(|(p, _, e)| p.is_paired().then_some(e.data.power))
            .sum::<i32>();
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_graph_rebuild",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            nloops = nloops,
            uncontracted_propagator_count = uncontracted_propagator_count,
            uncontracted_propagator_power_sum = uncontracted_propagator_power_sum,
            contracted_propagator_count = contracted_propagator_count,
            contracted_propagator_power_sum = contracted_propagator_power_sum,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );

        let lmb = graph.lmb();
        let mom_pat = function!(GS.emr_mom, W_.a_).to_pattern();
        for (p, e, ed) in graph.iter_edges() {
            if p.is_paired() {
                // println!("{e}");
                let loop_expr = lmb.loop_atom::<Atom>(e, GS.loop_mom, &[], false);

                ed.data
                    .mom
                    .pattern_match(&mom_pat, None, None)
                    .for_each(|m| {
                        // An affine vacuum routing can contain fixed source
                        // external momenta. Solve only for hard coordinates;
                        // otherwise H=Q-P can solve for P and leave Q in the
                        // integrated numerator as a spurious external variable.
                        // A child crown carrier can become an active loop of
                        // its parent. Retained hard coordinates take precedence
                        // over an external classification in another frame.
                        if usize::try_from(m[&W_.a_].as_view()).is_ok_and(|edge| {
                            let edge = EdgeIndex(edge);
                            source_lmbs.iter().any(|lmb| lmb.ext_from(edge).is_some())
                                && !source_lmbs.iter().any(|lmb| lmb.loop_edges.contains(&edge))
                        }) {
                            return;
                        }
                        let var = mom_pat.replace_wildcards(&m).unwrap();
                        if !momentum_variables.iter().any(|existing| existing == &var) {
                            momentum_variables.push(var);
                        }
                    });

                // println!("{loop_expr}");

                // println!("{external_expr}");
                let is_zero = &ed.data.mom - loop_expr;
                // println!("Momentum check for edge {}: {}", e, is_zero);
                system.push(is_zero);
            }
        }

        let add_additional_args = [
            Replacement::new(
                function!(GS.emr_mom, W_.i_).to_pattern(),
                function!(GS.emr_mom, W_.i_, W_.a___),
            )
            .allow_new_wildcards_on_rhs(true),
            Replacement::new(
                function!(GS.loop_mom, W_.i_).to_pattern(),
                function!(GS.loop_mom, W_.i_, W_.a___),
            )
            .allow_new_wildcards_on_rhs(true),
        ];
        let momentum_solution =
            VakintMomentumSolution::solve(&system, &momentum_variables, &add_additional_args)
                .wrap_err("could not solve momentum system for Vakint integrand")?;
        let numerator = SymbolicTensor::infer_indexed(cooking.try_cook(t.numerator.as_view())?)?;
        let mut numerator = momentum_solution.rewrite_numerator(&numerator)?;
        t.integral = momentum_solution.rewrite_integral(&t.integral, &add_additional_args);
        momentum_solution.ensure_free_variables_eliminated(
            term_index,
            "integral",
            &t.integral,
            &add_additional_args,
        )?;
        momentum_solution.ensure_free_variables_eliminated(
            term_index,
            "numerator",
            numerator.expression(),
            &add_additional_args,
        )?;
        let retained_hard_coordinates = retained_loop_edges
            .iter()
            .map(|edge| function!(GS.emr_mom, usize::from(*edge)))
            .collect::<Vec<_>>();
        momentum_solution.ensure_variables_eliminated(
            term_index,
            "numerator",
            numerator.expression(),
            &retained_hard_coordinates,
            &add_additional_args,
        )?;
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_momentum_solve",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = numerator.expression(),
            "Vakint trace"
        );

        // debug!(
        //     "Graph from vakint expression:\n{}\n{}",
        //     graph.dot_lmb_of(&graph.full_filter(), &graph.lmb()),
        //     term
        // );

        let normalization =
            additional_normalization.pow(&SymbolicTensor::infer_indexed(Atom::num(nloops))?)?;
        numerator = numerator.multiply(normalization)?;
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_loop_normalization",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            nloops = nloops,
            log.additional_normalization = additional_normalization.expression(),
            log.numerator = numerator.expression(),
            "Vakint trace after loop normalization"
        );

        // Vakint needs explicit tensor indices; only translate metric shorthands
        // to dot notation here, without reintroducing Schoonschip rank-1 factors.
        let numerator = numerator.to_dots()?;
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_metric_shorthand_to_dot",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = numerator.expression(),
            "Vakint trace"
        );
        t.integral = t
            .integral
            .replace(function!(GS.loop_mom, W_.x___))
            .with(function!(vakint::symbols::S.k, W_.x___))
            .replace(function!(GS.emr_mom, W_.x___))
            .with(function!(vakint::symbols::S.p, W_.x___));
        t.numerator = Integrated::to_vakint_numerator(&numerator)?;
        debug_tags!(#uv, #integrated, #vakint, #trace;
            stage = "to_vakint_integrand_term_after_vakint_symbols",
            term_index = %term_index,
            reduced = %reduced_label,
            dependent_subgraph = %dependent_subgraph_label,
            log.integral = t.integral,
            log.numerator = t.numerator,
            "Vakint trace"
        );
    }

    Ok(a)
}

struct VakintMomentumSolution {
    replacements: Vec<Replacement>,
    free_variables: Vec<Atom>,
}

impl VakintMomentumSolution {
    fn solve(
        system: &[Atom],
        variables: &[Atom],
        add_additional_args: &[Replacement],
    ) -> Result<Self> {
        if variables.is_empty() {
            return Ok(Self {
                replacements: vec![],
                free_variables: vec![],
            });
        }
        if system.is_empty() {
            return Ok(Self {
                replacements: vec![],
                free_variables: variables.to_vec(),
            });
        }

        let solutions = Atom::solve(system)
            .wrt_with_exponent::<u8, _>(variables)
            .map_err(|source| eyre!("{source}"))?;
        let [solution] = solutions.iter().as_slice() else {
            return Err(eyre!(
                "expected one Vakint momentum solution, got {} branches",
                solutions.len()
            ));
        };
        let solution = variables
            .iter()
            .map(|variable| {
                let polynomial_variable =
                    PolyVariable::try_from(variable.clone()).map_err(|source| eyre!("{source}"))?;
                if solution.free_variables().contains(&polynomial_variable) {
                    Ok(variable.clone())
                } else {
                    solution.get(&polynomial_variable).cloned().ok_or_else(|| {
                        eyre!("Vakint momentum solution has no value for {polynomial_variable}")
                    })
                }
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self::from_solution(
            &solution,
            variables,
            add_additional_args,
        ))
    }

    fn from_solution(
        solution: &[Atom],
        variables: &[Atom],
        add_additional_args: &[Replacement],
    ) -> Self {
        debug_assert_eq!(solution.len(), variables.len());
        let mut replacements = vec![];
        let mut free_variables = vec![];

        for index in (0..variables.len()).rev() {
            let replacement = &solution[index];
            let variable = &variables[index];
            let replacement = replacement.replace_multiple(&replacements);
            if replacement == variable {
                free_variables.push(variable.clone());
            } else {
                replacements.push(Replacement::new(
                    variable.replace_multiple(add_additional_args).to_pattern(),
                    replacement
                        .replace_multiple(add_additional_args)
                        .to_pattern(),
                ));
            }
        }

        Self {
            replacements,
            free_variables,
        }
    }

    fn ensure_free_variables_eliminated(
        &self,
        term_index: usize,
        expression_kind: &str,
        expression: &Atom,
        add_additional_args: &[Replacement],
    ) -> Result<()> {
        self.ensure_variables_eliminated(
            term_index,
            expression_kind,
            expression,
            self.free_variables.iter(),
            add_additional_args,
        )
    }

    fn ensure_variables_eliminated<'a>(
        &self,
        term_index: usize,
        expression_kind: &str,
        expression: &Atom,
        variables: impl IntoIterator<Item = &'a Atom>,
        add_additional_args: &[Replacement],
    ) -> Result<()> {
        for variable in variables {
            let bare_pattern = variable.to_pattern();
            let indexed_pattern = variable.replace_multiple(add_additional_args).to_pattern();
            if expression
                .pattern_match(&bare_pattern, None, None)
                .next()
                .is_some()
                || expression
                    .pattern_match(&indexed_pattern, None, None)
                    .next()
                    .is_some()
            {
                return Err(eyre!(
                    "Underdetermined Vakint momentum solve left free variable {} in {} of term {}",
                    variable,
                    expression_kind,
                    term_index
                ));
            }
        }
        Ok(())
    }

    fn rewrite_numerator(
        &self,
        numerator: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        // Production numerators have already passed metric-only preparation.
        // Affine momentum replacement preserves port labels and introduces no
        // metric edges. Keep vectors explicit while normalizing the new products.
        let expression = numerator.expression().replace_multiple(&self.replacements);
        let normalized =
            numerator
                .with_checked_expression(expression)?
                .contract(ContractSettings {
                    rank_one: false,
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })?;
        Ok(normalized)
    }

    fn rewrite_integral(&self, integral: &Atom, add_additional_args: &[Replacement]) -> Atom {
        let free_variable_zero_replacements =
            self.free_variable_zero_replacements(add_additional_args);
        Self::normalize_topology_momenta(
            integral
                .replace_multiple(&self.replacements)
                .replace_multiple(&free_variable_zero_replacements),
        )
    }

    fn free_variable_zero_replacements(
        &self,
        add_additional_args: &[Replacement],
    ) -> Vec<Replacement> {
        let mut replacements = Vec::with_capacity(2 * self.free_variables.len());
        for variable in &self.free_variables {
            replacements.push(Replacement::new(variable.to_pattern(), Atom::Zero));
            let indexed_variable = variable.replace_multiple(add_additional_args);
            if indexed_variable != *variable {
                replacements.push(Replacement::new(indexed_variable.to_pattern(), Atom::Zero));
            }
        }
        replacements
    }

    fn normalize_topology_momenta(expression: Atom) -> Atom {
        if !expression
            .replace(function!(vakint::symbols::S.topo, W_.x_))
            .matches()
        {
            return expression;
        }

        expression
            .replace(function!(
                vakint::symbols::S.prop,
                W_.edgeid_,
                W_.x_,
                W_.mom_,
                W_.mass_,
                W_.prop_
            ))
            .with_map(|matches| {
                function!(
                    vakint::symbols::S.prop,
                    matches.get(W_.edgeid_).unwrap().to_atom(),
                    matches.get(W_.x_).unwrap().to_atom(),
                    matches.get(W_.mom_).unwrap().to_atom().expand(),
                    matches.get(W_.mass_).unwrap().to_atom(),
                    matches.get(W_.prop_).unwrap().to_atom()
                )
            })
    }
}

#[cfg(test)]
mod tests {
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    use symbolica::symbol;

    use crate::{initialisation::test_initialise, utils::symbolica_ext::Q_I};

    #[test]
    fn integrated_projections_preserve_open_ports_when_a_laurent_part_vanishes() -> Result<()> {
        test_initialise()?;
        let vector = crate::numerator::AppliedFeynmanRule::new(parse!(
            "gammalooprs::Q(0,spenso::mink(4,1))"
        ))?
        .expr;
        let epsilon = Atom::var(GS.dim_epsilon);
        let source = vector.with_checked_expression(vector.expression() * epsilon.pow(-1))?;
        let integrated = IntegratedCts::from_series(&source, &series(source.expression(), 2)?, 4)?;
        assert_eq!(integrated.pole()?.structure(), vector.structure());
        assert_eq!(
            integrated.pole()?.expression(),
            &(vector.expression() * epsilon.pow(-1) * Atom::var(GS.integrated_loop_scale).pow(4)),
        );
        assert_eq!(
            integrated.physical_pole()?.expression(),
            source.expression()
        );
        let finite = integrated.physical_finite_counterterm()?;
        assert!(finite.is_zero());
        assert_eq!(finite.rank(), 1);
        assert_eq!(finite.structure(), vector.structure());

        let zero = vector.with_checked_expression(Atom::Zero)?;
        let integrated_zero = IntegratedCts::from_series(&zero, &series(zero.expression(), 2)?, 4)?;
        for projection in [
            integrated_zero.pole()?,
            integrated_zero.finite_counterterm()?,
            integrated_zero.physical_pole()?,
            integrated_zero.physical_finite_counterterm()?,
        ] {
            assert!(projection.is_zero());
            assert_eq!(projection.structure(), vector.structure());
        }
        Ok(())
    }

    #[test]
    fn integrated_products_keep_contracted_and_zero_tensor_interfaces() -> Result<()> {
        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let left = crate::numerator::AppliedFeynmanRule::new(parse!(
            "gammalooprs::Q(0,spenso::mink(4,1))"
        ))?
        .expr;
        let left_source =
            left.with_checked_expression(left.expression() * (epsilon.pow(-1) + 2))?;
        let left_ct =
            IntegratedCts::from_series(&left_source, &series(left_source.expression(), 2)?, 4)?;
        for (index, expected_rank) in [(1, 0), (2, 2)] {
            let right = crate::numerator::AppliedFeynmanRule::new(parse!(&format!(
                "gammalooprs::Q(1,spenso::mink(4,{index}))"
            )))?
            .expr;
            let right_source =
                right.with_checked_expression(right.expression() * (epsilon.pow(-1) + 3))?;
            let right_ct = IntegratedCts::from_series(
                &right_source,
                &series(right_source.expression(), 2)?,
                4,
            )?;
            let tensor_product = left.clone().multiply(right)?;
            assert_eq!(tensor_product.rank(), expected_rank);
            let product = IntegratedCts::factorized_product([&left_ct, &right_ct], 2)?;
            let pole = product.physical_pole()?;
            let finite = product.physical_finite_counterterm()?;
            assert_eq!(pole.structure(), tensor_product.structure());
            assert_eq!(finite.structure(), tensor_product.structure());
            // Native Laurent collection can restore an indexed closed pair.
            // Compare the exact coefficients in the shared contraction notation.
            let tensor_product = tensor_product.to_dots()?;
            let pole = pole
                .contract(ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })?
                .to_dots()?;
            let finite = finite
                .contract(ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })?
                .to_dots()?;
            assert_eq!(
                pole.expression(),
                &(tensor_product.expression() * epsilon.pow(-2))
            );
            assert_eq!(
                finite.expression(),
                &(tensor_product.expression() * Atom::num(6))
            );

            let zero = left.with_checked_expression(Atom::Zero)?;
            let zero_ct = IntegratedCts::from_series(&zero, &series(zero.expression(), 2)?, 4)?;
            let zero_product = IntegratedCts::factorized_product([&zero_ct, &right_ct], 2)?;
            for projection in [
                zero_product.physical_pole()?,
                zero_product.physical_finite_counterterm()?,
            ] {
                assert!(projection.is_zero());
                assert_eq!(projection.structure(), tensor_product.structure());
            }
        }
        Ok(())
    }

    // Fixtures enter at the source-expression boundary; production sectors
    // already own an admitted tensor and never pass through these adapters.
    fn simplify(expression: &Atom) -> Result<Atom> {
        use crate::numerator::{AppliedFeynmanRule, GetSingleAtom};
        let input = AppliedFeynmanRule::new(expression.clone())?;
        let output = super::simplify(&input.expr)?;
        Ok(AppliedFeynmanRule::from_tensor(output).get_single_atom()?)
    }

    fn to_vakint_numerator(expression: &Atom) -> Result<Atom> {
        let input = crate::numerator::AppliedFeynmanRule::new(expression.clone())?;
        Integrated::to_vakint_numerator(&input.expr)
    }

    fn to_vakint_integrand<
        E: UVE,
        V,
        H,
        S: SubGraphLike + SubSetLike<Base = SuBitGraph>,
        SS: SubGraphLike,
    >(
        expression: &Atom,
        graph: &HedgeGraph<E, V, H>,
        reduced: &S,
        dependent_subgraph: &SS,
        active_components: &[(SuBitGraph, SuBitGraph, LoopMomentumBasis)],
        settings: &VakintSettings,
        substitute_masses_to_m_uv: bool,
    ) -> Result<VakintExpression> {
        let input = crate::numerator::AppliedFeynmanRule::new(expression.clone())?;
        super::to_vakint_integrand(
            &input.expr,
            graph,
            reduced,
            dependent_subgraph,
            active_components,
            settings,
            substitute_masses_to_m_uv,
        )
    }

    use super::*;

    fn exact_one_loop_constants(atom: &Atom) -> Atom {
        // MATAD's symbolic mode retains Gamma/PolyGamma atoms, whereas AlphaLoop
        // uses z3. Its numerical constant substitution would approximate zeta(3),
        // so compare the two backends using only these exact identities. They
        // cover the Gamma(1+epsilon) expansion needed through epsilon^2 of I2.
        atom.replace(function!(vakint_symbol!("Gamma"), 1))
            .with(Atom::one())
            .replace(function!(vakint_symbol!("PolyGamma"), 0, 1))
            .with(-Atom::var(vakint_symbol!("EulerGamma")))
            .replace(function!(vakint_symbol!("PolyGamma"), 1, 1))
            .with(Atom::var(symbolica::atom::Symbol::PI).pow(2) / Atom::num(6))
            .replace(function!(vakint_symbol!("PolyGamma"), 2, 1))
            .with(-Atom::num(2) * Atom::var(vakint_symbol!("z3")))
    }

    #[test]
    fn fixed_momentum_components_have_the_formal_dimensional_wick_projection() {
        test_initialise().unwrap();
        let vakint = crate::utils::vakint().unwrap();
        let settings = VakintSettings::default().true_settings();
        let spectator = parse!("(component_spectator_a+component_spectator_b)^11");
        let dimension = Atom::var(vakint_symbol!(&settings.epsilon_symbol)) * -2 + 4;
        let square = vakint::symbols::S.dot(
            &function!(vakint::symbols::S.k, 1),
            &function!(vakint::symbols::S.k, 1),
        );
        for (axes, coefficient) in [
            (vec![0, 0], Atom::one() / &dimension),
            (vec![1, 1], -Atom::one() / &dimension),
            (vec![0, 1], Atom::Zero),
            (
                vec![0, 0, 0, 0],
                Atom::num(3) / (&dimension * (&dimension + 2)),
            ),
            (
                vec![0, 0, 1, 1],
                -Atom::one() / (&dimension * (&dimension + 2)),
            ),
        ] {
            let source = &spectator
                * Atom::mul_many(
                    axes.iter()
                        .map(|axis| function!(GS.loop_mom, 1, GS.cind(*axis))),
                );
            let mut term = vakint::VakintTerm {
                integral: vakint::vakint_parse!("topo(I1L(component_mass_squared,3))").unwrap(),
                numerator: to_vakint_numerator(&source).unwrap(),
                vectors: vec![("k".into(), 1)],
            };
            term.tensor_reduce(vakint, &settings).unwrap();
            let resolved = Integrated::new(vakint, &settings)
                .simplify_projected_numerator(&term.numerator, 0)
                .unwrap();
            // Cancel the unchanged spectator and invariant powers literally;
            // only the rational coefficient in D is brought to a common denominator.
            let actual = (resolved / (&spectator * square.pow(axes.len() / 2))).collect_factors();
            assert!(
                (&actual - &coefficient).together().is_zero(),
                "axes={axes:?}: {actual} != {coefficient}"
            );
            assert_eq!(
                source,
                &spectator
                    * Atom::mul_many(axes.iter().map(|axis| function!(
                        GS.loop_mom,
                        1,
                        GS.cind(*axis)
                    )))
            );
        }
    }

    #[test]
    fn physical_four_components_keep_joint_dimensional_angular_moments() -> Result<()> {
        use crate::numerator::symbolica_ext::NumeratorAtomExt;

        test_initialise()?;
        let vakint = crate::utils::vakint()?;
        let settings = VakintSettings::default().true_settings();
        assert_eq!(
            settings.tensor_reduction_method,
            vakint::TensorReductionMethod::FeynKit
        );
        let integrated = Integrated::new(vakint, &settings);
        let epsilon_symbol = vakint_symbol!(&settings.epsilon_symbol);
        let epsilon = Atom::var(epsilon_symbol);
        assert_eq!(epsilon_symbol, GS.dim_epsilon);
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let spectator = parse!("(angular_lift_spectator_a+angular_lift_spectator_b)^11");
        let component = |loop_id, axis| function!(GS.loop_mom, loop_id, GS.cind(axis));
        let physical_pair = |left, right| {
            Atom::add_many((0..4).map(|axis| {
                let sign = if (Minkowski {}).is_neg(axis) { -1 } else { 1 };
                Atom::num(sign) * component(left, axis) * component(right, axis)
            }))
        };
        let k_four = physical_pair(1, 1);
        let l_four = physical_pair(2, 2);
        let mixed_four = physical_pair(1, 2);
        let k = function!(vakint::symbols::S.k, 1);
        let l = function!(vakint::symbols::S.k, 2);
        let kk = vakint::symbols::S.dot(&k, &k);
        let ll = vakint::symbols::S.dot(&l, &l);
        let kl = vakint::symbols::S.dot(&k, &l);
        let invariant_a = &kk * &ll;
        let invariant_b = kl.pow(2);
        let keys = [kk.clone(), ll.clone(), kl.clone()];
        let assert_coefficients_zero = |label: &str, expression: &Atom| -> Result<()> {
            let rows = expression
                .coefficient_list_exact(&keys)
                .ok_or_else(|| eyre!("{label}: bounded invariant collection refused"))?;
            for (_, coefficient) in rows {
                assert!(coefficient.together().is_zero(), "{label}: {coefficient}");
            }
            Ok(())
        };

        // Joint O(D) isotropy gives T_mnrs=a*g_mn*g_rs
        // + b*(g_mr*g_ns+g_ms*g_nr). The two full-D contractions
        // determine a,b from kk*ll and (k.l)^2, independently of a projector.
        let denominator = &dimension * (&dimension - Atom::one()) * (&dimension + Atom::num(2));
        let a = ((&dimension + Atom::one()) * &invariant_a - Atom::num(2) * &invariant_b)
            / &denominator;
        let b = (&dimension * &invariant_b - &invariant_a) / &denominator;
        assert_coefficients_zero(
            "isotropy full trace",
            &(dimension.pow(2) * &a + Atom::num(2) * &dimension * &b - &invariant_a),
        )?;
        assert_coefficients_zero(
            "isotropy crossed trace",
            &(&dimension * &a + &dimension * (&dimension + Atom::one()) * &b - &invariant_b),
        )?;
        let rank = Atom::num(4);
        let expected_kl = rank.pow(2) * &a + Atom::num(2) * &rank * &b;
        let expected_mixed_square = &rank * &a + &rank * (&rank + Atom::one()) * &b;
        assert_coefficients_zero(
            "physical-four double trace",
            &(&expected_kl
                - (Atom::num(8) * (Atom::num(2) * &dimension + Atom::one()) * &invariant_a
                    + Atom::num(8) * (&dimension - Atom::num(4)) * &invariant_b)
                    / &denominator),
        )?;
        assert_coefficients_zero(
            "physical-four crossed trace",
            &(&expected_mixed_square
                - (Atom::num(4) * (&dimension - Atom::num(4)) * &invariant_a
                    + Atom::num(4) * (Atom::num(5) * &dimension - Atom::num(2)) * &invariant_b)
                    / &denominator),
        )?;

        // For z=n.P.n, these are <z^2> and <z^3> for a rank-r
        // projector. Conditional isotropy, Tr(P)=r and P^2=P give
        // B=(r E2-E3)/(D-1), C=(E2-E3)/(D-1). Keep the same joint
        // angular dimension for both loop vectors at rank six.
        let e2 = &rank * (&rank + Atom::num(2)) / (&dimension * (&dimension + Atom::num(2)));
        let e3 = &e2 * (&rank + Atom::num(4)) / (&dimension + Atom::num(4));
        let b6 = (&rank * &e2 - &e3) / (&dimension - Atom::one());
        let c6 = (&e2 - &e3) / (&dimension - Atom::one());
        let invariant_six_a = &kk * ll.pow(2);
        let invariant_six_b = &ll * kl.pow(2);
        let expected_kl_squared = &b6 * &invariant_six_a + (&e3 - &b6) * &invariant_six_b;
        let expected_mixed_square_l = &c6 * &invariant_six_a + (&e3 - &c6) * &invariant_six_b;
        let denominator_six = &denominator * (&dimension + Atom::num(4));
        assert_coefficients_zero(
            "physical-four rank-six K L squared",
            &(&expected_kl_squared
                - (Atom::num(96) * (&dimension + Atom::num(2)) * &invariant_six_a
                    + Atom::num(96) * (&dimension - Atom::num(4)) * &invariant_six_b)
                    / &denominator_six),
        )?;
        assert_coefficients_zero(
            "physical-four rank-six M squared L",
            &(&expected_mixed_square_l
                - (Atom::num(24) * (&dimension - Atom::num(4)) * &invariant_six_a
                    + Atom::num(24)
                        * (Atom::num(7) * &dimension - Atom::num(4))
                        * &invariant_six_b)
                    / &denominator_six),
        )?;

        for (label, physical, expected) in [
            ("K", k_four.clone(), Atom::num(4) * &kk / &dimension),
            ("L", l_four.clone(), Atom::num(4) * &ll / &dimension),
            ("K L", &k_four * &l_four, expected_kl.clone()),
            (
                "M squared",
                mixed_four.pow(2),
                expected_mixed_square.clone(),
            ),
            (
                "K L squared",
                &k_four * l_four.pow(2),
                expected_kl_squared.clone(),
            ),
            (
                "M squared L",
                mixed_four.pow(2) * &l_four,
                expected_mixed_square_l.clone(),
            ),
        ] {
            let source = &spectator * physical;
            let retained_source = source.clone();
            let mut term = vakint::VakintTerm {
                // This connected vacuum topology requires one JOINT angular
                // average; its mixed moments do not factor into powers of 4/D.
                integral: vakint::vakint_parse!("topo(I2L(angular_lift_mass_squared,1,1,1))")
                    .unwrap(),
                numerator: to_vakint_numerator(&source)?,
                vectors: vec![("k".into(), 1), ("k".into(), 2)],
            };
            term.tensor_reduce(vakint, &settings)?;
            let resolved = integrated.simplify_projected_numerator(&term.numerator, 0)?;
            assert!(
                resolved
                    .pattern_match(&spectator.to_pattern(), None, None)
                    .next()
                    .is_some(),
                "{label}: powered scalar spectator must remain a literal factor"
            );
            // Strip only the unchanged spectator. Collect complete radial
            // invariants; together receives their rational coefficients only.
            let actual = resolved.replace(spectator.to_pattern()).with(1);
            assert_coefficients_zero(label, &(actual - expected))?;
            assert_eq!(source, retained_source, "{label}: source mutation");
        }

        // Freezing the angular dimension at four before Laurent expansion
        // loses these independently nonzero finite rank-four and rank-six coefficients.
        for (label, expected, wrong, finite_control, double_pole_control) in [
            (
                "K L",
                expected_kl,
                invariant_a.clone(),
                Atom::num((19, 18)) * &invariant_a - Atom::num((2, 9)) * &invariant_b,
                Atom::num((31, 36)) * &invariant_a - Atom::num((1, 3)) * &invariant_b,
            ),
            (
                "M squared",
                expected_mixed_square,
                invariant_b.clone(),
                -Atom::num((1, 9)) * &invariant_a + Atom::num((17, 18)) * &invariant_b,
                -Atom::num((1, 6)) * &invariant_a + Atom::num((25, 36)) * &invariant_b,
            ),
            (
                "K L squared",
                expected_kl_squared,
                invariant_six_a.clone(),
                Atom::num((17, 12)) * &invariant_six_a - Atom::num((1, 3)) * &invariant_six_b,
                Atom::num((199, 144)) * &invariant_six_a - Atom::num((7, 12)) * &invariant_six_b,
            ),
            (
                "M squared L",
                expected_mixed_square_l,
                invariant_six_b.clone(),
                -Atom::num((1, 12)) * &invariant_six_a + Atom::num((7, 6)) * &invariant_six_b,
                -Atom::num((7, 48)) * &invariant_six_a + Atom::num((17, 18)) * &invariant_six_b,
            ),
        ] {
            assert_coefficients_zero(
                "early-four angular value",
                &(expected.replace(epsilon_symbol).with(0) - &wrong),
            )?;
            let difference = expected - wrong;
            let lost = truncate(&series(&(&difference / &epsilon), 0)?, true);
            assert_coefficients_zero(label, &(lost - &finite_control))?;
            // A double enclosing pole consumes the independently retained
            // epsilon-squared angular coefficient as a finite contribution.
            let double_pole = series(&(&difference / epsilon.pow(2)), 0)?;
            assert_coefficients_zero(
                "double-pole angular Laurent series",
                &(double_pole.to_atom() - &finite_control / &epsilon - &double_pole_control),
            )?;
            assert_coefficients_zero(
                "double-pole finite angular loss",
                &(truncate(&double_pole, true) - &double_pole_control),
            )?;
            for control in [&finite_control, &double_pole_control] {
                let rows = control
                    .coefficient_list_exact(&keys)
                    .ok_or_else(|| eyre!("{label}: finite negative-control collection refused"))?;
                assert_eq!(rows.len(), 2);
                assert!(
                    rows.iter()
                        .all(|(_, coefficient)| !coefficient.together().is_zero())
                );
            }
        }
        // These admit angular mechanisms, without asserting that an actual W
        // input needs rank six. They are not an actual D09 R2 numerator, UV/IR
        // separation or finite renormalization-condition proof.
        Ok(())
    }

    #[test]
    fn projected_loop_dummies_do_not_capture_retained_numeric_ports() {
        test_initialise().unwrap();
        let vakint = crate::utils::vakint().unwrap();
        let settings = VakintSettings::default().true_settings();
        let dimension = Atom::var(vakint_symbol!(&settings.epsilon_symbol)) * -2 + 4;
        let spectator = parse!("(retained_port_a+retained_port_b)^11");
        let loop_momentum = function!(vakint::symbols::S.k, 1);
        let source = &spectator
            * 2
            * vakint::symbols::S.dot(&function!(vakint::symbols::S.p, 3), &loop_momentum)
            * function!(
                vakint::symbols::S.k,
                1,
                Minkowski {}.to_symbolic([Atom::var(GS.dim), Atom::num(2)])
            );
        let mut term = vakint::VakintTerm {
            integral: vakint::vakint_parse!("topo(I1L(retained_port_mass_squared,3))").unwrap(),
            numerator: source.clone(),
            vectors: vec![("k".into(), 1)],
        };
        term.tensor_reduce(vakint, &settings).unwrap();
        let resolved = Integrated::new(vakint, &settings)
            .simplify_projected_numerator(&term.numerator, 0)
            .unwrap();
        let retained = function!(
            vakint::symbols::S.p,
            3,
            Minkowski {}.to_symbolic([Atom::var(GS.dim), Atom::num(2)])
        );
        let square = vakint::symbols::S.dot(&loop_momentum, &loop_momentum);
        let coefficient = (resolved / (&spectator * retained * square)).collect_factors();
        assert!(
            (coefficient - Atom::num(2) / dimension)
                .together()
                .is_zero()
        );
        assert!(source.contains(&spectator));
    }

    #[test]
    fn qcd_insertion_ward_certificate_checks_rules_and_propagators() -> Result<()> {
        test_initialise()?;
        let model = crate::utils::load_generic_model("sm");
        let source = include_str!(
            "../../../../../tests/resources/graphs/paper_appendix_b1_nested_gluon_self_energy.dot"
        );
        let mut graphs = Graph::import_dot(source, &model)?;
        assert_eq!(graphs.len(), 1);
        let graph = graphs.remove(0);
        let external: SuBitGraph = graph.external_filter();
        let parent = graph.full_filter().subtract(&external);
        let mut child = graph.empty_subgraph::<SuBitGraph>();
        for edge in [EdgeIndex(3), EdgeIndex(4)] {
            child.add(graph[&edge].1);
        }
        assert!(graph.has_qcd_self_energy_insertion_ward_identity(&model, &parent, &child));

        let node = graph.iter_nodes_of(&parent).next().unwrap().0;
        let mut explicit = graph.clone();
        explicit.underlying[node].num.autogenerated = false;
        let original = &explicit.underlying[node].num.value;
        explicit.underlying[node].num.value = original
            .with_checked_expression(original.expression() * parse!("ward_noncanonical_vertex"))?;
        assert!(!explicit.has_qcd_self_energy_insertion_ward_identity(&model, &parent, &child));
        // The metadata bit cannot certify the actual rule: a custom factor
        // marked generated must fail the same native tensor comparison.
        explicit.underlying[node].num.autogenerated = true;
        assert!(!explicit.has_qcd_self_energy_insertion_ward_identity(&model, &parent, &child));

        let mut altered_propagator = graph.clone();
        let original = &altered_propagator.underlying[EdgeIndex(2)].num.value;
        altered_propagator.underlying[EdgeIndex(2)].num.value = original.with_checked_expression(
            original.expression() * parse!("ward_noncanonical_propagator"),
        )?;
        assert!(
            !altered_propagator
                .has_qcd_self_energy_insertion_ward_identity(&model, &parent, &child)
        );

        // Canonical model records are immutable. Rebuild through their
        // validated serialization boundary, then instantiate the changed rule
        // in a newly imported graph instead of mutating an Arc or a runtime ID.
        let rule = model.vertex_rule_by_id(graph.underlying[node].vertex_rule.unwrap())?;
        let lorentz_name = &model
            .lorentz_structure_by_id(rule.lorentz_structures[0])?
            .name;
        let mut definition = serde_json::to_value(&model)?;
        let lorentz = definition["lorentz_structures"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|entry| entry["name"] == *lorentz_name)
            .unwrap();
        lorentz["structure"] = serde_json::to_value(
            (function!(UFO.identity, UFO.idx(1, 1), UFO.idx(2, 1))
                * function!(UFO.momentum, UFO.idx(3, 1), 1))
            .to_canonical_string(),
        )?;
        let nonminimal_model: crate::model::Model = serde_json::from_value(definition)?;
        let mut graphs = Graph::import_dot(source, &nonminimal_model)?;
        assert_eq!(graphs.len(), 1);
        let mut nonminimal = graphs.remove(0);
        assert_ne!(
            nonminimal.underlying[node].num.value,
            graph.underlying[node].num.value
        );
        nonminimal.underlying[node].num.autogenerated = true;
        assert!(nonminimal.underlying[node].num.autogenerated);
        assert!(!nonminimal.has_qcd_self_energy_insertion_ward_identity(
            &nonminimal_model,
            &parent,
            &child
        ));

        let mut graphs = Graph::import_dot(
            include_str!(
                "../../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
            ),
            &model,
        )?;
        assert_eq!(graphs.len(), 1);
        let quark_parent = graphs.remove(0);
        let external: SuBitGraph = quark_parent.external_filter();
        assert!(!quark_parent.has_qcd_self_energy_insertion_ward_identity(
            &model,
            &quark_parent.full_filter().subtract(&external),
            &child
        ));
        Ok(())
    }

    #[test]
    fn equal_mass_gluon_bubble_and_mass_insertion_obey_dimensional_ward_identity() -> Result<()> {
        test_initialise()?;
        let settings = vakint::VakintSettings {
            number_of_terms_in_epsilon_expansion: 5,
            ..VakintSettings::default().true_settings()
        };
        let vakint = crate::utils::vakint()?;
        let integrator = Integrated::new(vakint, &settings);
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let momentum = function!(GS.loop_mom, 1, mink.to_symbolic([]));
        let radial = function!(SPENSO_TAG.dot, &momentum, &momentum);
        let mass = parse!("ward_quark_mass");
        let gamma =
            |index: &Atom| function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, index);
        let trace_two = function!(
            SPENSO_TAG.trace,
            Bispinor {}.new_rep(4).to_symbolic([]),
            gamma(&mu),
            gamma(&nu)
        );
        let trace_four = function!(
            SPENSO_TAG.trace,
            Bispinor {}.new_rep(4).to_symbolic([]),
            gamma(&mu),
            gamma(&momentum),
            gamma(&nu),
            gamma(&momentum)
        );
        let bubble = &trace_four + mass.pow(2) * &trace_two;
        // S(r)^2 = dS(r)/dm. This is one insertion on the equal-mass
        // fermion loop, equivalently half the mass derivative at q=0.
        let mass_insertion =
            &mass * (&radial + mass.pow(2)) * &trace_two + Atom::num(2) * &mass * &trace_four;

        for (label, numerator, power) in
            [("bubble", bubble, 2), ("mass insertion", mass_insertion, 3)]
        {
            let mut term = vakint::VakintTerm {
                integral: function!(
                    vakint::symbols::S.topo,
                    function!(vakint_symbol!("I1L"), mass.pow(2), power)
                ),
                numerator: to_vakint_numerator(&simplify(&numerator)?)?,
                vectors: vec![("k".into(), 1)],
            };
            term.tensor_reduce(vakint, &settings)?;
            debug_tags!(#uv, #integrated, #vakint, #trace;
                stage = "ward_reference_after_tensor_projection",
                label,
                file.numerator = %term.numerator.to_canonical_string(),
                file.restored = %Integrated::restore_numerator(term.numerator.clone(), 0).to_canonical_string(),
                "Dimensional Ward reference at the native tensor boundary"
            );
            term.numerator = integrator
                .simplify_projected_numerator(&term.numerator, 0)
                .wrap_err_with(|| {
                    format!("{label}: dimensional algebra after tensor projection")
                })?;
            term.evaluate_integral(vakint, &settings)?;
            assert!(
                series(&term.numerator, 2)?.to_atom().expand().is_zero(),
                "the {label} constant soft coefficient must vanish before multiplication by any Laurent-series child: {}",
                term.numerator
            );
        }

        // Independently certify the all-d identity before choosing an epsilon
        // depth. The two IBP recurrences are exact meromorphic identities.
        let d = Atom::var(GS.dim);
        let i1 = parse!("ward_I1");
        let i2 = (&d - Atom::num(2)) * &i1 / (Atom::num(2) * mass.pow(2));
        let i3 = (&d - Atom::num(4)) * &i2 / (Atom::num(4) * mass.pow(2));
        for identity in [
            (Atom::num(2) / &d - Atom::one()) * &i1 + Atom::num(2) * mass.pow(2) / &d * &i2,
            (Atom::num(4) / &d - Atom::one()) * &i2 + Atom::num(4) * mass.pow(2) / &d * &i3,
        ] {
            assert!(
                identity
                    .try_to_rational_polynomial::<_, _, u32>(&*Q_I, &*Q_I, None)?
                    .numerator
                    .is_zero()
            );
        }
        // The linear soft jet has an odd number of loop vectors and vanishes
        // under the same d-dimensional angular projector, not by setting d=4.
        let mut odd = vakint::VakintTerm {
            integral: function!(
                vakint::symbols::S.topo,
                function!(vakint_symbol!("I1L"), mass.pow(2), 3)
            ),
            numerator: to_vakint_numerator(&(radial * function!(GS.loop_mom, 1, &mu)))?,
            vectors: vec![("k".into(), 1)],
        };
        odd.tensor_reduce(vakint, &settings)?;
        assert!(odd.numerator.is_zero());
        Ok(())
    }

    #[test]
    fn physical_on_shell_finite_prefix_keeps_epsilon_squared_for_two_enclosing_poles() -> Result<()>
    {
        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let zeta_three = Atom::var(vakint_symbol!("z3"));
        let normalization = Atom::i() / (Atom::num(16) * pi.pow(2));
        for (backend, evaluation_order) in [
            (
                "AlphaLoop",
                vakint::EvaluationOrder(vec![vakint::EvaluationMethod::AlphaLoop(
                    vakint::AlphaLoopOptions {
                        susbstitute_masters: false,
                    },
                )]),
            ),
            (
                "MATAD",
                vakint::EvaluationOrder::matad_only(Some(vakint::MATADOptions {
                    susbstitute_masters: false,
                    direct_numerical_substition: false,
                    ..Default::default()
                })),
            ),
        ] {
            let settings = vakint::VakintSettings {
                // A one-loop child times two enclosing poles needs precisely
                // epsilon^-1, epsilon^0, epsilon^1 and epsilon^2.
                number_of_terms_in_epsilon_expansion: 4,
                evaluation_order,
                ..VakintSettings::default().true_settings()
            };
            let mut vacuum = vakint::VakintTerm {
                integral: function!(
                    vakint::symbols::S.topo,
                    function!(vakint_symbol!("I1L"), Atom::var(GS.m_uv_vacuum).pow(2), 2)
                ),
                // Before the forest minus V=-(1-epsilon) I2 (p-slash-m).
                // The signed local counterterm therefore integrates -V.
                numerator: Atom::one() - &epsilon,
                vectors: vec![("k".into(), 1)],
            };
            vacuum.evaluate_integral(crate::utils::vakint()?, &settings)?;
            let child_expansion = series(
                &exact_one_loop_constants(&vacuum.numerator)
                    .replace(GS.mu_r_sq)
                    .with(Atom::one())
                    .replace(GS.m_uv_vacuum)
                    .with(Atom::one()),
                2,
            )?;
            let child_source = SymbolicTensor::infer_indexed(child_expansion.to_atom())?;
            let child = IntegratedCts::from_series(&child_source, &child_expansion, 4)?;
            let finite_child = child.physical_finite_counterterm()?;
            // log(exp(gamma_E epsilon) Gamma(1+epsilon)) has coefficients
            // zeta(2)/2 and -zeta(3)/3. At mu=M, the finite addback is
            // N[1-zeta(2) epsilon/2+(zeta(2)/2+zeta(3)/3) epsilon^2].
            for (enclosing_poles, child_coefficient) in [
                (1, -pi.pow(2) / Atom::num(12)),
                (2, pi.pow(2) / Atom::num(12) + &zeta_three / Atom::num(3)),
            ] {
                let parent_poles = (-&normalization * epsilon.pow(-1)).pow(enclosing_poles);
                let complete_source = finite_child
                    .with_checked_expression(finite_child.expression() * &parent_poles)?;
                let complete = IntegratedCts::from_series(
                    &complete_source,
                    &series(complete_source.expression(), 0)?,
                    4 * (1 + enclosing_poles),
                )?;
                let insufficient_child =
                    series(finite_child.expression(), (enclosing_poles - 1) as usize)?.to_atom();
                let insufficient_source =
                    finite_child.with_checked_expression(insufficient_child * &parent_poles)?;
                let insufficient = IntegratedCts::from_series(
                    &insufficient_source,
                    &series(insufficient_source.expression(), 0)?,
                    complete.scale_power,
                )?;
                let finite = complete.physical_finite_counterterm()?;
                let expected =
                    -&normalization * (-&normalization).pow(enclosing_poles) * child_coefficient;
                assert!(
                    (finite.expression() - &expected).expand().is_zero(),
                    "{backend}: a one-loop OS prefix lost the finite term from {enclosing_poles} enclosing poles"
                );
                assert!(!expected.is_zero());
                assert!(
                    insufficient.physical_finite_counterterm()?.is_zero(),
                    "{backend}: truncating below epsilon^{enclosing_poles} must lose this nonzero finite contribution"
                );
            }
        }
        Ok(())
    }

    #[test]
    fn nested_os_child_gluon_parent_keeps_local_mass_and_vacuum_scheme_input() -> Result<()> {
        use crate::uv::{
            UVgenerationSettings,
            approx::local_4d::{
                FourDTerm, Full4dCts, on_shell_forest_tests::nested_gluon_fixture, uv_limit,
            },
        };

        test_initialise()?;
        let (graph, parent, child, root) = nested_gluon_fixture()?;
        assert!(graph.has_qcd_self_energy_insertion_ward_identity(
            &graph.model,
            parent.subgraph(),
            child.subgraph()
        ));
        let settings = UVgenerationSettings {
            generate_integrated: true,
            add_marker: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
        let seed = Full4dCts::from_factorized_local(&Local4dCts::root());
        let child_local = uv_limit(&seed, &ctx, &child, &root, &child, &root)?;
        let vakint_settings = vakint::VakintSettings {
            number_of_terms_in_epsilon_expansion: 4,
            ..settings.vakint.true_settings()
        };
        let child_integrated = Integrated::new(crate::utils::vakint()?, &vakint_settings).run(
            &child_local,
            &ctx,
            &child,
            &root,
            &child,
            &root,
        )?;
        assert!(
            !child_integrated
                .physical_finite_counterterm()?
                .expression()
                .derivative(GS.dim_epsilon)
                .replace(GS.dim_epsilon)
                .with(Atom::Zero)
                .is_zero()
        );
        let completed = Full4dCts::recursion_input(
            &child_local,
            &child_integrated,
            ApproximationType::OS,
            false,
            child.lmb(),
        )?;
        let parent_local = uv_limit(&completed, &ctx, &parent, &child, &parent, &child)?;
        assert!(
            parent_local.atom().expression().contains_symbol(GS.os_mom),
            "the local parent H must retain the physical child on-shell contribution"
        );
        let mass = graph[EdgeIndex(3)].mass_atom(&graph.model);
        let mut frozen_prefixes = 0;
        for sector in parent_local.active_sectors() {
            assert_ne!(
                sector.integration_atom(),
                &sector.atom,
                "the completed-child soft branch stays local under the certified Ward policy"
            );
            let analytic = GS.resolve_os_momenta(sector.integration_atom().expression());
            assert!(!analytic.contains_symbol(GS.os_mom));
            let terms = FourDTerm::from_view(analytic.as_view())?;
            assert!(!terms.is_empty());
            for denominator in terms.iter().flat_map(|term| &term.denominators) {
                assert!(
                    denominator.mass_squared == Atom::var(GS.m_uv_expansion).pow(2)
                        || denominator.mass_squared == Atom::var(GS.m_uv_vacuum).pow(2),
                    "the analytic policy retained a non-vacuum mass: {}",
                    denominator.mass_squared,
                );
                assert!(
                    denominator
                        .full_expr
                        .pattern_match(&mass.to_pattern(), None, None)
                        .next()
                        .is_none()
                );
                assert!(!denominator.momentum.contains_symbol(GS.os_mom));
            }
            if !sector.frozen_lmbs().is_empty() {
                frozen_prefixes += 1;
                assert_eq!(sector.frozen_lmbs(), std::slice::from_ref(child.lmb()));
                assert!(
                    sector.atom.expression().contains_symbol(GS.dim_epsilon),
                    "the enclosing Taylor operation must retain the finite child's positive epsilon orders"
                );
                assert_eq!(sector.active_components.len(), 1);
                assert!(
                    sector.active_components[0]
                        .2
                        .loop_edges
                        .iter()
                        .all(|edge| { !child.lmb().loop_edges.contains(edge) })
                );
            }
        }
        assert_eq!(
            frozen_prefixes, 1,
            "the integrated child coefficient has one localization owner"
        );

        // The same physical graph without a completed child is not covered
        // by this Ward proof: its bare two-loop soft branch remains analytic.
        let bare_parent = uv_limit(&seed, &ctx, &parent, &root, &parent, &root)?;
        for sector in bare_parent.active_sectors() {
            assert_eq!(sector.integration_atom(), &sector.atom);
        }
        Ok(())
    }

    #[test]
    fn analytic_os_temporal_dot_keeps_scalar_and_numerator_factors() {
        test_initialise().unwrap();
        let vakint = crate::utils::vakint().unwrap();
        let settings = VakintSettings::default().true_settings();
        let vector_slot = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let temporal = GS.energy_delta(&vector_slot);
        let mass = parse!("os_mass");
        let spectator = parse!("(os_a+os_b)^3*(os_c+os_d)");
        let reserved = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy_at(1_000_000).to_atom()]);
        let reserved_vector = SPENSO_TAG.rank_one_tensor_symbol("os_reserved_vector");
        let source_contraction =
            function!(GS.loop_mom, 7, &reserved) * function!(reserved_vector, &reserved);
        let second_dot = function!(
            SPENSO_TAG.dot,
            &temporal,
            function!(GS.emr_mom, 8, &vector_slot)
        );
        for (source_head, converted_head) in [
            (GS.loop_mom, vakint::symbols::S.k),
            (GS.emr_mom, vakint::symbols::S.p),
        ] {
            let vector = function!(source_head, 7, &vector_slot);
            for coefficient in [Atom::one(), mass.clone(), parse!("os_mass*(os_x+os_y)")] {
                for power in [1, 2] {
                    for (left, right) in [
                        (&coefficient * &temporal, vector.clone()),
                        (vector.clone(), &coefficient * &temporal),
                    ] {
                        let source = &spectator
                            * &source_contraction
                            * &second_dot
                            * function!(SPENSO_TAG.dot, left, right).pow(power);
                        let converted = to_vakint_numerator(&source).unwrap();
                        let expected = &spectator
                            * function!(vakint::symbols::S.k, 7, &reserved)
                            * function!(reserved_vector, &reserved)
                            * vakint::symbols::S
                                .dot(&temporal, &function!(vakint::symbols::S.p, 8))
                            * coefficient.pow(power)
                            * vakint::symbols::S
                                .dot(&temporal, &function!(converted_head, 7))
                                .pow(power);
                        assert_eq!(
                            converted, expected,
                            "vector={vector}, coefficient={coefficient}, power={power}"
                        );
                        assert!(!converted.contains_symbol(SPENSO_TAG.dot));
                        assert_eq!(
                            to_vakint_numerator(&Integrated::restore_numerator(
                                converted.clone(),
                                0
                            ))
                            .unwrap(),
                            converted,
                            "native fixed-vector and momentum ports roundtrip exactly"
                        );
                        assert!(
                            converted
                                .pattern_match(&spectator.to_pattern(), None, None)
                                .next()
                                .is_some()
                        );
                        let mut term = vakint::VakintTerm {
                            integral: vakint::vakint_parse!("topo(I1L(os_mass_squared,3))")
                                .unwrap(),
                            numerator: converted,
                            vectors: vec![("k".into(), 7)],
                        };
                        // Vakint opens the compact powers with distinct dummy
                        // scopes while preserving the source's reserved index.
                        term.tensor_reduce(vakint, &settings).unwrap();
                    }
                }
            }
        }
    }

    #[test]
    fn projected_temporal_and_external_slashes_share_a_native_interface() -> Result<()> {
        use idenso::shorthands::UndoShorthands;

        test_initialise()?;
        let vakint = Vakint::new()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let c = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let compact = mink.to_symbolic([]);
        let temporal = GS.energy_delta(&compact);
        let k = function!(GS.loop_mom, 1, &compact);
        let p = function!(GS.emr_mom, 0, &compact);
        let mass = Atom::var(symbol!("projected_temporal::mass"; Scalar));
        let spectator = (Atom::var(symbol!("projected_temporal::x"; Scalar))
            + Atom::var(symbol!("projected_temporal::y"; Scalar)))
        .pow(11);
        let epsilon = Atom::var(GS.dim_epsilon);
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let source = &spectator
            * (Atom::var(GS.dim) - Atom::num(4))
            * (function!(SPENSO_TAG.dot, &k, Atom::num(3) * &mass * &temporal)
                * function!(GS.loop_mom, 1, &mu)
                * function!(AGS.gamma, &a, &b, &mu)
                * function!(AGS.gamma0, &b, &c)
                + function!(SPENSO_TAG.dot, &k, &p)
                    * function!(GS.loop_mom, 1, &nu)
                    * function!(AGS.gamma, &a, &c, &nu));
        // Independently, <k_mu k_nu> = k^2 g_mu_nu/d and
        // slash(e0) gamma0 = 1. Both terms therefore share the open (a,c)
        // spinor interface; the evanescent factor is kept before pole expansion.
        let reference = &spectator * Atom::num(-2) * &epsilon / &dimension
            * function!(SPENSO_TAG.dot, &k, &k)
            * (Atom::num(3) * &mass * function!(ETS.metric, &a, &c)
                + function!(AGS.gamma, &a, &c, &p));
        let expected = to_vakint_numerator(&reference)?;
        // FeynKit consumes the native compact fixed-vector form directly.
        // The independent cross-backend comparison also materializes exactly
        // those vector arguments, as required by AlphaLoop's existing input
        // notation, without distributing surrounding numerator factors.
        for (tensor_reduction_method, compact_fixed_vectors) in [
            (vakint::TensorReductionMethod::FeynKit, true),
            (vakint::TensorReductionMethod::FeynKit, false),
            (vakint::TensorReductionMethod::AlphaLoop, false),
        ] {
            for project_onto_tensor_integrals in [false, true] {
                for use_dot_product_notation in [false, true] {
                    let settings = vakint::VakintSettings {
                        number_of_terms_in_epsilon_expansion: 3,
                        tensor_reduction_method,
                        project_onto_tensor_integrals,
                        use_dot_product_notation,
                        ..VakintSettings::default().true_settings()
                    };
                    let numerator = if compact_fixed_vectors {
                        to_vakint_numerator(&source)?
                    } else {
                        source
                            .undo_dots::<Aind>()?
                            .replace(function!(GS.loop_mom, W_.x___))
                            .with(function!(vakint::symbols::S.k, W_.x___))
                            .replace(function!(GS.emr_mom, W_.x___))
                            .with(function!(vakint::symbols::S.p, W_.x___))
                    };
                    let mut term = vakint::VakintTerm {
                        integral: vakint::vakint_parse!("topo(I1L(projected_mass_squared,3))")
                            .unwrap(),
                        numerator,
                        vectors: vec![("k".into(), 1)],
                    };
                    debug_tags!(#uv, #integrated, #vakint, #inspect;
                        stage = "mixed_native_tensor_projection_input",
                        backend = ?tensor_reduction_method,
                        compact = compact_fixed_vectors,
                        projected = project_onto_tensor_integrals,
                        dots = use_dot_product_notation,
                        file.numerator = %term.numerator.to_canonical_string(),
                        "Exact mixed fixed-vector input before tensor projection"
                    );
                    term.tensor_reduce(&vakint, &settings).wrap_err_with(|| format!(
                        "mixed native input: backend={tensor_reduction_method:?}, compact={compact_fixed_vectors}, projected={project_onto_tensor_integrals}, dots={use_dot_product_notation}"
                    ))?;
                    let resolved = Integrated::new(&vakint, &settings)
                        .simplify_projected_numerator(&term.numerator, 0)?;
                    assert!(
                        resolved
                            .pattern_match(&spectator.to_pattern(), None, None)
                            .next()
                            .is_some()
                    );
                    // Remove the identical opaque spectator before exact scalar
                    // coefficient cancellation; never distribute its powered sum.
                    let difference = resolved.replace(spectator.to_pattern()).with(1)
                        - expected.replace(spectator.to_pattern()).with(1);
                    let difference = Integrated::restore_numerator(difference, 0);
                    let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
                    assert!(
                        SymbolicTensor::infer_indexed(
                            cooking.try_cook(simplify(&difference)?.as_view())?
                        )?
                        .coefficients_are_zero(
                            spenso::shadowing::TensorCollectFilter::<0>::Tensors
                        )?
                        .is_true(),
                        "backend={tensor_reduction_method:?}, compact={compact_fixed_vectors}, projected={project_onto_tensor_integrals}, dots={use_dot_product_notation}"
                    );
                    let mut scalar = vakint::VakintTerm {
                        integral: term.integral.clone(),
                        numerator: vakint::symbols::S.dot(
                            function!(vakint::symbols::S.k, 1),
                            function!(vakint::symbols::S.k, 1),
                        ),
                        vectors: term.vectors.clone(),
                    };
                    term.numerator = resolved;
                    term.evaluate_integral(&vakint, &settings)?;
                    scalar.evaluate_integral(&vakint, &settings)?;
                    let spin = Atom::num(3) * &mass * function!(ETS.metric, &a, &c)
                        + function!(AGS.gamma, &a, &c, &p);
                    // Expand only the scalar Laurent coefficient. The spin
                    // matrix and powered spectator remain literal factors in
                    // both the result and this independent reference.
                    let scalar_reference = series(
                        &(Atom::num(-2) * &epsilon / &dimension * &scalar.numerator),
                        1,
                    )?
                    .to_atom();
                    let expected_integral = &spectator * spin * scalar_reference;
                    assert!(
                        term.numerator
                            .pattern_match(&spectator.to_pattern(), None, None)
                            .next()
                            .is_some()
                    );
                    let difference = Integrated::restore_numerator(
                        term.numerator.replace(spectator.to_pattern()).with(1),
                        0,
                    ) - expected_integral.replace(spectator.to_pattern()).with(1);
                    let difference = simplify(&difference)?;
                    assert!(
                        SymbolicTensor::infer_indexed(cooking.try_cook(difference.as_view())?)?
                            .coefficients_are_zero(
                                spenso::shadowing::TensorCollectFilter::<0>::Tensors
                            )?
                            .is_true(),
                        "backend={tensor_reduction_method:?}, compact={compact_fixed_vectors}, projected={project_onto_tensor_integrals}, dots={use_dot_product_notation}: {difference}"
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn projected_fixed_velocity_clifford_algebra_retains_evanescent_factors() -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let settings = VakintSettings::default().true_settings();
        let integrated = Integrated::new(&vakint, &settings);
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let c = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let velocity = [
            Atom::num((5, 4)),
            Atom::num((1, 4)),
            Atom::num((1, 2)),
            Atom::num((1, 2)),
        ];
        let slash = |left: &Atom, right: &Atom, slot: &Atom| {
            &velocity[0] * function!(AGS.gamma0, left, right)
                + Atom::add_many(
                    (1..4)
                        .map(|axis| &velocity[axis] * function!(GS.delta_vec, GS.cind(axis), slot)),
                ) * function!(AGS.gamma, left, right, slot)
        };
        // This is the fixed-vector product exposed by angular projection:
        // spatial basis vectors sit outside the ordered spinor factors. Its
        // reduction needs the same chain/vector boundary as the OS remainder.
        let spin = slash(&a, &b, &mu) * slash(&b, &c, &nu);
        let identity = function!(ETS.metric, &a, &c);
        let spectator = (Atom::var(symbol!("projected_fixed_velocity::x"; Scalar))
            + Atom::var(symbol!("projected_fixed_velocity::y"; Scalar)))
        .pow(11);
        let epsilon = Atom::var(GS.dim_epsilon);
        let d_minus_four = Atom::var(GS.dim) - Atom::num(4);
        for power in 0..=2 {
            let numerator = &spectator * d_minus_four.pow(power) * &spin;
            let projected = to_vakint_numerator(&numerator)?;
            let resolved = integrated.simplify_projected_numerator(&projected, 0)?;
            let expected = &spectator * (Atom::num(-2) * &epsilon).pow(power) * &identity;
            // A literal equality certifies both slash(v)^2 = 1 and retained
            // scalar factorization; no expanded diagnostic copy is needed.
            assert_eq!(resolved, expected, "evanescent power={power}");
            assert!(
                resolved
                    .pattern_match(&spectator.to_pattern(), None, None)
                    .next()
                    .is_some(),
                "the powered spectator must remain a single factor"
            );
            for pole_order in 1..=2 {
                let evaluated = series(&(&resolved * epsilon.pow(-pole_order)), 1)?.to_atom();
                let reference = &spectator
                    * Atom::num(-2).pow(power)
                    * epsilon.pow(power - pole_order)
                    * &identity;
                assert_eq!(evaluated, series(&reference, 1)?.to_atom());
                if power == pole_order {
                    assert_eq!(truncate(&series(&reference, 1)?, true), reference);
                    assert!(!reference.is_zero(), "pole times evanescent terms survive");
                    let prematurely_four_dimensional = numerator.replace(GS.dim).with(4);
                    assert!(prematurely_four_dimensional.is_zero());
                }
            }
        }
        Ok(())
    }

    #[test]
    fn production_on_shell_scheme_integral_matches_bare_fermion_reference() -> Result<()> {
        use crate::uv::{
            UVgenerationSettings,
            approx::local_4d::{Full4dCts, on_shell_tests::fixture, uv_limit},
        };

        test_initialise()?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let (graph, current, root) = fixture(false)?;
        for velocity in [
            ["1", "0", "0", "0"],
            ["5/4", "3/4", "0", "0"],
            ["5/4", "1/4", "1/2", "1/2"],
        ] {
            let settings = UVgenerationSettings {
                on_shell_subtraction_velocity: velocity.map(str::to_owned),
                generate_integrated: true,
                add_marker: false,
                ..Default::default()
            };
            let ctx = UVCtx::new(&graph, &graph.model, &settings);
            let local = uv_limit(
                &Full4dCts::from_factorized_local(&Local4dCts::root()),
                &ctx,
                &current,
                &root,
                &current,
                &root,
            )?;
            let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);

            // Determine only the graph's Feynman-rule phase, coupling, color, and
            // open-spin normalization from its bare numerator N=C[(2-d)k/+dm].
            // No Taylor result enters the oracle. The integrated field matrix is
            // C(p/-m), independently of which internal line carries the loop k.
            let numerator = graph
                .numerator(current.subgraph(), root.subgraph())?
                .to_d_dim(GS.dim)?
                .color_simplify()?
                .state
                .expr;
            let numerator = cooking.uncook(numerator.expression().as_view());
            let massless = numerator.replace(mass.to_pattern()).with(Atom::Zero);
            let external_slash = massless
                .replace(function!(GS.emr_mom, 2, W_.x___))
                .with(function!(GS.emr_mom, 0, W_.x___).to_pattern());
            let field_matrix = simplify(
                &(external_slash / (Atom::num(2) - Atom::var(GS.dim))
                    - (&numerator - &massless) / Atom::var(GS.dim)),
            )?;
            let field_matrix = field_matrix.replace(GS.dim).with(4);
            assert!(!field_matrix.is_zero());
            assert!(field_matrix.contains_symbol(GS.emr_mom));

            let pi = Atom::var(symbolica::atom::Symbol::PI);
            let normalization = Atom::i() / (Atom::num(16) * pi.pow(2));
            let epsilon = Atom::var(GS.dim_epsilon);
            let zeta_three = Atom::var(vakint_symbol!("z3"));
            for (backend, evaluation_order) in [
                (
                    "AlphaLoop",
                    vakint::EvaluationOrder(vec![vakint::EvaluationMethod::AlphaLoop(
                        vakint::AlphaLoopOptions {
                            susbstitute_masters: false,
                        },
                    )]),
                ),
                (
                    "MATAD",
                    vakint::EvaluationOrder::matad_only(Some(vakint::MATADOptions {
                        susbstitute_masters: false,
                        direct_numerical_substition: false,
                        ..Default::default()
                    })),
                ),
            ] {
                for project_onto_tensor_integrals in [true, false] {
                    for enclosing_poles in [1, 2] {
                        let vakint_settings = vakint::VakintSettings {
                            // Whole-forest depth: three terms for one enclosing pole,
                            // four for two. The one-loop run must preserve that request.
                            number_of_terms_in_epsilon_expansion: enclosing_poles as i64 + 2,
                            project_onto_tensor_integrals,
                            evaluation_order: evaluation_order.clone(),
                            ..settings.vakint.true_settings()
                        };
                        let integrated = Integrated::new(crate::utils::vakint()?, &vakint_settings)
                            .run(&local, &ctx, &current, &root, &current, &root)?;
                        let pole = integrated.physical_pole()?;
                        let pole = cooking.uncook(pole.expression().as_view());
                        assert!(
                            !pole.contains_symbol(GS.m_uv_vacuum)
                                && !pole.contains_symbol(GS.m_uv_expansion),
                            "the integrated field pole must be auxiliary-scale independent before choosing a scale"
                        );
                        let pole = pole.replace(GS.dim).with(4);
                        let expected_pole = &normalization * &field_matrix * epsilon.pow(-1);
                        let mut comparisons =
                            vec![("signed field pole".to_owned(), &pole - expected_pole, 0)];
                        for vacuum_mass in [1, 2] {
                            // For real positive scales, L = log(mu²)-log(M²) is also
                            // -log(M²/mu²). Use the ratio form without assuming a generic
                            // complex-log product identity inside the symbolic engine.
                            let logarithm = -function!(
                                symbolica::atom::Symbol::LOG,
                                Atom::num(vacuum_mass).pow(2) / Atom::var(GS.mu_r_sq)
                            );
                            let finite = integrated.physical_finite_counterterm()?;
                            let finite = cooking
                                .uncook(finite.expression().as_view())
                                .replace(GS.m_uv_vacuum)
                                .with(Atom::num(vacuum_mass))
                                .replace(GS.dim)
                                .with(4);
                            let expected_finite = &normalization
                                * &field_matrix
                                * (Atom::one() - &logarithm
                                    + &epsilon
                                        * (&logarithm
                                            - logarithm.pow(2) / Atom::num(2)
                                            - pi.pow(2) / Atom::num(12))
                                    + epsilon.pow(2)
                                        * (-logarithm.pow(3) / Atom::num(6)
                                            + logarithm.pow(2) / Atom::num(2)
                                            + (Atom::one() - &logarithm) * pi.pow(2)
                                                / Atom::num(12)
                                            + &zeta_three / Atom::num(3)));
                            // Spin coefficients four, but the same D-dimensional
                            // angular/radial measure: -V4 = I2 * C(p-slash-m).
                            // The independent vacuum test below derives this from
                            // I2 and I3; this Gamma series is not read from the CT.
                            let vacuum_master = epsilon.pow(-1)
                                + &logarithm
                                + &epsilon
                                    * (logarithm.pow(2) / Atom::num(2) + pi.pow(2) / Atom::num(12))
                                + epsilon.pow(2)
                                    * (logarithm.pow(3) / Atom::num(6)
                                        + &logarithm * pi.pow(2) / Atom::num(12)
                                        - &zeta_three / Atom::num(3));
                            let signed_four_vacuum =
                                &normalization * &field_matrix * &vacuum_master;
                            // Actual -VD = pole - finite; hence VD - V4 below.
                            comparisons.push((
                                format!("dimensional vacuum bridge at M={vacuum_mass}"),
                                &finite - &pole + signed_four_vacuum
                                    - &normalization * &field_matrix * &epsilon * vacuum_master,
                                enclosing_poles,
                            ));
                            comparisons.push((
                                format!("finite scheme correction at M={vacuum_mass}"),
                                finite - expected_finite,
                                enclosing_poles,
                            ));
                        }
                        for (label, difference, depth) in comparisons {
                            let difference =
                                series(&simplify(&exact_one_loop_constants(&difference))?, depth)?
                                    .to_atom();
                            assert!(
                                difference
                                    .try_to_rational_polynomial::<_, _, u32>(&*Q_I, &*Q_I, None)?
                                    .numerator
                                    .is_zero(),
                                "{backend}, velocity={velocity:?}, project={project_onto_tensor_integrals}, enclosing poles={enclosing_poles}: {label} differs from the independent one-loop fermion reference: {difference}"
                            );
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn boosted_on_shell_production_branches_obey_the_integrated_pole_condition() -> Result<()> {
        use crate::{
            numerator::symbolica_ext::NumeratorAtomExt,
            uv::{
                UVgenerationSettings,
                approx::local_4d::{
                    Full4dCts,
                    on_shell_tests::{
                        SpinReduction, assert_rational_zero, fixture, gamma_velocity_at_port,
                    },
                    uv_limit,
                },
                overlap_numerators::NumeratorJetCache,
            },
        };
        use idenso::{color::ColorSimplifier, shorthands::UndoShorthands};
        use linnet::half_edge::involution::Hedge;
        use spenso::structure::slot::IsAbstractSlot;

        test_initialise()?;
        let (graph, current, root) = fixture(false)?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let numerator_source = graph
            .numerator(current.subgraph(), root.subgraph())?
            .to_d_dim(GS.dim)?
            .color_simplify()?
            .state
            .expr;
        let bare = numerator_source.with_checked_expression(
            numerator_source.expression()
                * graph.denominator(current.subgraph(), &graph.model, |_| -1),
        )?;
        let numerator = cooking.uncook(numerator_source.expression().as_view());
        let mass = graph[EdgeIndex(2)].mass_atom(&graph.model);
        let epsilon = Atom::var(GS.dim_epsilon);
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let normalization =
            Atom::i() / (Atom::num(16) * Atom::var(symbolica::atom::Symbol::PI).pow(2));
        // This common scalar master is exp(gamma_E*epsilon) Gamma(epsilon)
        // (mu^2/m^2)^epsilon. Keep it independent: the pole condition below is
        // an exact meromorphic identity, not a cancellation at one scale.
        let physical_master = Atom::var(symbol!("boosted_os_reference::physical_master"));
        // With the massive quark routed as k, Feynman combination and a shift
        // give Delta=x^2 m^2 and <k>=(1-x)p. These two independent moments are
        // int_0^1 dx x^(-2epsilon) and int_0^1 dx (1-x)x^(-2epsilon).
        let j0 = Atom::one() / (Atom::one() - Atom::num(2) * &epsilon);
        let j1 = &j0 / (Atom::num(2) - Atom::num(2) * &epsilon);
        let evaluate_vector = |atom: &Atom, edge: usize, value: &Atom| {
            GS.resolve_os_momenta(&GS.erase_uv_momentum_provenance(atom).replace_map(
                |view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() > 0
                        && usize::try_from(momentum.get(0)).ok() == Some(edge)
                    {
                        let indices = momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>();
                        **output = GS.indexed_momentum(value, &indices);
                    }
                },
            ))
        };
        let regulate = |atom: &Atom| -> Result<Atom> {
            Ok(Integrated::dimensionally_regularized(&simplify(atom)?).map_mink_dim(4))
        };
        let assert_integrated_zero = |label: &str, atom: &Atom| -> Result<()> {
            let reduced = simplify(atom)?;
            let certificate =
                NumeratorJetCache::default().canonicalize(&reduced.undo_dots::<Aind>()?)?;
            assert!(
                certificate.is_zero,
                "{label}: nonzero finite spinor tensor certificate {}",
                certificate.key,
            );
            Ok(())
        };
        let numerator = regulate(&numerator)?;
        let at_zero = evaluate_vector(&numerator, 2, &Atom::Zero);
        let mass_matrix = &at_zero / &dimension;
        let expected_mass = &mass_matrix * (Atom::num(3) - Atom::num(2) * &epsilon) * &j0;
        // Freeze only scalar spin coefficients after native contraction and
        // dimensional regularization. The same D-dimensional j0/j1 moments
        // below continue the four-dimensional numerator; do not freeze them.
        let four_numerator = numerator.replace(GS.dim_epsilon).with(0);
        let four_at_zero = at_zero.replace(GS.dim_epsilon).with(0);
        let loop_vector = function!(GS.emr_mom, 2);

        // A nonzero trace of one exact physical specialization proves that the
        // complete matrix is not identically zero. At mu=m, the physical master
        // is 1/epsilon + O(epsilon): its zero-order coefficient vanishes.
        // J0=1+2*epsilon+... and (2-d)*J1=-J0 therefore give
        // Fin[R/normalization] = -2*C*(p-slash-m). Its spin trace, with the
        // fixture's normalized color closure, is twice the trace of C*m.
        let ports = [0, 1].map(|hedge| {
            graph.underlying[Hedge(hedge)]
                .edge_spin_slots()
                .next()
                .expect("the primitive fermion has one spinor port per boundary edge")
                .to_atom()
        });
        let trace_projection = function!(ETS.metric, &ports[0], &ports[1])
            * cooking.uncook(graph.global_prefactor.num.expression().as_view());
        let strong_coupling = Atom::from(crate::model::UFOSymbol::from("G"));
        let vertex_coupling = Atom::from(crate::model::UFOSymbol::from("GC_11"));
        let coupling_at_unit_g = graph
            .model
            .coupling("GC_11")?
            .expression
            .replace(strong_coupling.to_pattern())
            .with(1);
        assert_eq!(
            coupling_at_unit_g,
            Atom::i(),
            "the model's quark-gluon coupling is i*g_s at the exact point g_s=1"
        );
        let specialize_witness = |atom: &Atom| {
            atom.to_cof_dimension_invariants()
                .replace(mass.to_pattern())
                .with(1)
                .replace(vertex_coupling.to_pattern())
                .with(coupling_at_unit_g.clone())
        };
        let mass_trace = specialize_witness(&(&trace_projection * &mass_matrix))
            .replace(GS.dim_epsilon)
            .with(0);
        let mut witness_tensors = NumeratorJetCache::default();
        let mass_trace = witness_tensors.canonicalize(&mass_trace.undo_dots::<Aind>()?)?;
        assert_eq!(
            mass_trace.key,
            Atom::one(),
            "the reference trace is an exact scalar"
        );
        assert!(
            matches!(mass_trace.scalar_prefactor.as_view(), AtomView::Num(_))
                && !mass_trace.scalar_prefactor.is_zero(),
            "the independent bare mass matrix must give a nonzero exact number at m=g_s=1"
        );
        let expected_finite_trace = Atom::num(2) * &mass_trace.scalar_prefactor;

        for velocity in [
            ["1", "0", "0", "0"],
            ["5/4", "3/4", "0", "0"],
            ["5/4", "1/4", "1/2", "1/2"],
        ] {
            let settings = UVgenerationSettings {
                on_shell_subtraction_velocity: velocity.map(str::to_owned),
                generate_integrated: true,
                add_marker: false,
                ..Default::default()
            };
            let velocity = settings.on_shell_velocity()?;
            let branches = graph.on_shell_branches(
                &graph.model,
                &bare,
                current.subgraph(),
                current.lmb(),
                &velocity,
            )?;
            let mut on_shell_moment = Atom::Zero;
            let mut four_on_shell_moment = Atom::Zero;
            for (sigma, branch) in [1, -1].into_iter().zip(branches) {
                // Certify the factorized source and full denominator data
                // before replacing the two physical propagators by moments.
                let sector = FourDSector::new(
                    branch.clone(),
                    vec![(
                        current.subgraph().clone(),
                        current.subgraph().clone(),
                        current.lmb().clone(),
                    )],
                    Vec::new(),
                );
                let _certificate = sector.canonical_projection(&graph)?;
                let point =
                    GS.on_shell_momentum(EdgeIndex(0), &(Atom::num(sigma) * &mass), &velocity);
                assert!(
                    !branch
                        .expression()
                        .contains_symbol(GS.integrated_loop_scale)
                );
                for term in sector.physical_terms()? {
                    assert_eq!(term.denominators.len(), 2);
                    assert_eq!(
                        term.denominators
                            .iter()
                            .map(|denominator| denominator.source_edge)
                            .collect::<std::collections::BTreeSet<_>>(),
                        std::collections::BTreeSet::from([EdgeIndex(2), EdgeIndex(3)]),
                    );
                    for denominator in &term.denominators {
                        let (momentum, mass_squared) = if denominator.source_edge == EdgeIndex(2) {
                            (loop_vector.clone(), mass.pow(2))
                        } else {
                            (&loop_vector - &point, Atom::Zero)
                        };
                        assert_eq!(denominator.mass_squared, mass_squared);
                        let expected = (0..4).fold(-mass_squared, |sum, mu| {
                            sum + Atom::num(if mu == 0 { 1 } else { -1 })
                                * GS.indexed_momentum(&momentum, &[GS.cind(mu)]).pow(2)
                        });
                        // This is only a quadratic denominator polynomial;
                        // expand its fixed components, never the numerator.
                        let expected = GS.resolve_os_momenta(&expected).expand();
                        assert_rational_zero(
                            "Feynman-parameter source denominator",
                            &graph,
                            &current,
                            &(GS.erase_uv_momentum_provenance(&denominator.full_expr) - expected),
                            SpinReduction::SymbolicFirst,
                        )?;
                    }
                    let numerator =
                        regulate(&cooking.uncook(term.numerator.expression().as_view()))?;
                    let zero = evaluate_vector(&numerator, 2, &Atom::Zero);
                    let endpoint = evaluate_vector(&numerator, 2, &point);
                    on_shell_moment += &zero * &j0 + (endpoint - zero) * &j1;
                    // The companion uses this same certified production branch
                    // and both complete source denominators, not a second O oracle.
                    let four = numerator.replace(GS.dim_epsilon).with(0);
                    let four_zero = evaluate_vector(&four, 2, &Atom::Zero);
                    let four_endpoint = evaluate_vector(&four, 2, &point);
                    four_on_shell_moment += &four_zero * &j0 + (four_endpoint - four_zero) * &j1;
                }
            }
            assert_integrated_zero(
                "integrated production O is the boost-independent mass counterterm",
                &(&on_shell_moment - &expected_mass),
            )?;

            let os_difference =
                series(&((&on_shell_moment - &four_on_shell_moment) / &epsilon), 0)?.to_atom();
            assert_integrated_zero(
                "the standalone dimensional O difference is minus the mass matrix",
                &(&os_difference + &mass_matrix),
            )?;
            let os_difference_trace = specialize_witness(&(&trace_projection * &os_difference));
            let os_difference_trace =
                witness_tensors.canonicalize(&os_difference_trace.undo_dots::<Aind>()?)?;
            assert_eq!(os_difference_trace.key, Atom::one());
            assert_eq!(
                os_difference_trace.scalar_prefactor,
                -&mass_trace.scalar_prefactor
            );
            assert!(!os_difference_trace.scalar_prefactor.is_zero());

            let ctx = UVCtx::new(&graph, &graph.model, &settings);
            let local = uv_limit(
                &Full4dCts::from_factorized_local(&Local4dCts::root()),
                &ctx,
                &current,
                &root,
                &current,
                &root,
            )?;
            let vakint_settings = settings.vakint.true_settings();
            let integrated = Integrated::new(crate::utils::vakint()?, &vakint_settings)
                .run(&local, &ctx, &current, &root, &current, &root)?;
            let signed_field_pole = integrated.physical_pole()?;
            let signed_field_pole = cooking.uncook(signed_field_pole.expression().as_view());
            let actual_finite_addback = integrated.physical_finite_counterterm()?;
            let actual_finite_addback =
                cooking.uncook(actual_finite_addback.expression().as_view());
            for external_velocity in [
                ["1", "0", "0", "0"],
                ["5/4", "3/4", "0", "0"],
                ["5/4", "1/4", "1/2", "1/2"],
            ] {
                let external_velocity = UVgenerationSettings {
                    on_shell_subtraction_velocity: external_velocity.map(str::to_owned),
                    ..Default::default()
                }
                .on_shell_velocity()?;
                let point = GS.on_shell_momentum(EdgeIndex(0), &mass, &external_velocity);
                let bare_moment =
                    &at_zero * &j0 + (evaluate_vector(&numerator, 2, &point) - &at_zero) * &j1;
                // The local forest has -O-V, and the analytic addback is
                // +Fin(V). Its remaining field contribution is -Pole(V),
                // exactly the signed pole returned by the production policy.
                let pole =
                    evaluate_vector(&evaluate_vector(&signed_field_pole, 0, &point), 1, &point);
                let renormalized =
                    &normalization * &physical_master * (&bare_moment - &on_shell_moment)
                        + regulate(&pole)?;
                let right_projected = (&renormalized
                    + gamma_velocity_at_port(&graph, &renormalized, 0, &external_velocity))
                    / Atom::num(2);
                assert_integrated_zero(
                    "integrated boosted on-shell spinor condition",
                    &right_projected,
                )?;
                // Only the required epsilon^0 coefficient is extracted.
                // O(epsilon) terms of the physical master cannot contribute to
                // this finite coefficient, since the moment matrix is regular.
                let at_reference_scale = renormalized
                    .replace(physical_master.to_pattern())
                    .with(epsilon.pow(-1))
                    / &normalization;
                let finite_matrix = truncate(&series(&at_reference_scale, 0)?, true);
                let finite_trace = specialize_witness(&(&trace_projection * finite_matrix));
                let finite_trace =
                    witness_tensors.canonicalize(&finite_trace.undo_dots::<Aind>()?)?;
                assert_eq!(
                    finite_trace.key,
                    Atom::one(),
                    "the physical finite trace must reduce to an exact scalar"
                );
                assert_eq!(
                    finite_trace.scalar_prefactor, expected_finite_trace,
                    "a nonzero finite trace proves the complete MS-bar-field matrix survives"
                );

                let four_bare_moment = &four_at_zero * &j0
                    + (evaluate_vector(&four_numerator, 2, &point) - &four_at_zero) * &j1;
                let external_slash = (evaluate_vector(&numerator, 2, &point) - &at_zero)
                    / (Atom::num(2) - &dimension);
                let field_matrix = simplify(&(&external_slash - &mass_matrix))?;
                let bare_difference =
                    series(&((&bare_moment - &four_bare_moment) / &epsilon), 0)?.to_atom();
                assert_integrated_zero(
                    "the dimensional bare difference is C(p-slash-2m)",
                    &(bare_difference - &external_slash + Atom::num(2) * &mass_matrix),
                )?;
                let physical_log = -function!(
                    symbolica::atom::Symbol::LOG,
                    mass.pow(2) / Atom::var(GS.mu_r_sq)
                );
                // Only the finite open matrix is compared here. These common
                // masters need the pole and constant, not a D=4 angular average.
                let bare_master = epsilon.pow(-1) + &physical_log;
                let formal_finite = truncate(
                    &series(
                        &(renormalized
                            .replace(physical_master.to_pattern())
                            .with(bare_master.clone())
                            / &normalization),
                        0,
                    )?,
                    true,
                );
                assert_integrated_zero(
                    "the formal finite open matrix has MS-bar field normalization",
                    &(&formal_finite + (&physical_log + Atom::num(2)) * &field_matrix),
                )?;
                for vacuum_mass in [1, 2] {
                    let vacuum_log = -function!(
                        symbolica::atom::Symbol::LOG,
                        Atom::num(vacuum_mass).pow(2) / Atom::var(GS.mu_r_sq)
                    );
                    // The existing fixed_boosted_os_matches_compact_operator_and_local_shell_projection
                    // certifies the actual physical V4 source. Together with that
                    // prerequisite, the radial test transfers -V4 = I2*C(p/-m).
                    // Its contribution here is local; the addback is the actual
                    // same-production Integrated value, with its one loop scale consumed.
                    let local_four = &bare_master * (&four_bare_moment - &four_on_shell_moment)
                        + (epsilon.pow(-1) + &vacuum_log) * &field_matrix;
                    let addback = actual_finite_addback
                        .replace(GS.m_uv_vacuum)
                        .with(Atom::num(vacuum_mass));
                    let addback = evaluate_vector(&evaluate_vector(&addback, 0, &point), 1, &point);
                    let emitted_finite = truncate(
                        &series(&(&local_four + regulate(&addback)? / &normalization), 0)?,
                        true,
                    );
                    assert_integrated_zero(
                        "emitted spin-four local terms plus actual dimensional addback",
                        &(&emitted_finite - &formal_finite),
                    )?;
                    let four_addback = -&vacuum_log * &field_matrix;
                    let missing_rational =
                        truncate(&series(&(local_four + four_addback), 0)?, true);
                    assert_integrated_zero(
                        "prematurely spin-four addback misses exactly the field matrix",
                        &(&emitted_finite - missing_rational - &field_matrix),
                    )?;
                    // It is not enough to project this negative onto an on-shell
                    // spinor, where it vanishes. Keep its nonzero open-matrix trace.
                    let field_trace = specialize_witness(&(&trace_projection * &field_matrix));
                    let field_trace =
                        witness_tensors.canonicalize(&field_trace.undo_dots::<Aind>()?)?;
                    assert_eq!(field_trace.key, Atom::one());
                    assert_eq!(field_trace.scalar_prefactor, -&mass_trace.scalar_prefactor);
                    assert!(!field_trace.scalar_prefactor.is_zero());
                }
                // A finite bridge is not an all-epsilon child identity. At
                // mu=M=m, the exact local difference starts at 3*epsilon*C(p/-m).
                // This scalar +1/epsilon residue is a mechanism control only.
                let local_difference = (&bare_moment - &on_shell_moment - &four_bare_moment
                    + &four_on_shell_moment)
                    / &epsilon
                    - &field_matrix;
                assert_integrated_zero(
                    "the local dimensional difference has zero finite coefficient",
                    &series(&local_difference, 0)?.to_atom(),
                )?;
                let enclosing_pole_finite =
                    truncate(&series(&(local_difference / &epsilon), 0)?, true);
                assert_integrated_zero(
                    "a positive-epsilon difference survives a scalar enclosing pole",
                    &(enclosing_pole_finite - Atom::num(3) * &field_matrix),
                )?;
            }
        }
        Ok(())
    }

    #[test]
    fn on_shell_vacuum_remainder_matches_dimensional_one_loop_reference() -> Result<()> {
        test_initialise()?;
        let vakint = crate::utils::vakint()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let dimension = Atom::num(4) - Atom::num(2) * &epsilon;
        let mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let zeta_three = Atom::var(vakint_symbol!("z3"));
        let normalization = Atom::i() / (Atom::num(16) * pi.pow(2));

        // Independently contract the Feynman-gauge numerator before comparing
        // its angular average. Internal gamma algebra is d-dimensional; the
        // on-shell projector and external spinor slots are four-dimensional.
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let rho = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let [a, b, c, d] = std::array::from_fn(|_| bis.to_symbolic([Aind::new_dummy().to_atom()]));
        let contracted = simplify(
            &(function!(AGS.gamma, &a, &b, &mu)
                * function!(AGS.gamma, &b, &c, &rho)
                * function!(AGS.gamma, &c, &d, &mu)),
        )?;
        let expected =
            simplify(&((Atom::num(2) - Atom::var(GS.dim)) * function!(AGS.gamma, &a, &d, &rho)))?;
        assert!((contracted - expected).expand().is_zero());

        // Feynman parameters give int_0^1 dt t^(-2 epsilon)[2+(2-2 epsilon)t].
        // Derive its coefficients from exact logarithmic moments
        // int_0^1 dt t^a log(t)^n = (-1)^n n!/(a+1)^(n+1).
        // The scalar master has one pole, so epsilon^3 of this weight is
        // necessary and sufficient for the physical mass term at epsilon^2.
        let mut on_shell_weight = Atom::Zero;
        for order in 0..=3_i64 {
            let factorial = (1..=order).product::<i64>();
            let exponential_coefficient =
                Atom::num((-2_i64).pow(order as u32)) / Atom::num(factorial);
            for radial_power in [0_i64, 1] {
                let moment = Atom::num((-1_i64).pow(order as u32) * factorial)
                    / Atom::num(radial_power + 1).pow(order + 1);
                let contribution =
                    Atom::num(2) * &exponential_coefficient * moment * epsilon.pow(order);
                on_shell_weight += &contribution;
                if radial_power == 1 && order < 3 {
                    on_shell_weight -= contribution * &epsilon;
                }
            }
        }
        let closed_weight =
            (Atom::num(3) - Atom::num(2) * &epsilon) / (Atom::one() - Atom::num(2) * &epsilon);
        assert!(
            series(&(&on_shell_weight - closed_weight), 3)?
                .to_atom()
                .is_zero()
        );

        for (backend, evaluation_order) in [
            (
                "AlphaLoop",
                vakint::EvaluationOrder(vec![vakint::EvaluationMethod::AlphaLoop(
                    vakint::AlphaLoopOptions {
                        susbstitute_masters: false,
                    },
                )]),
            ),
            (
                "MATAD",
                vakint::EvaluationOrder::matad_only(Some(vakint::MATADOptions {
                    susbstitute_masters: false,
                    direct_numerical_substition: false,
                    ..Default::default()
                })),
            ),
        ] {
            let settings = vakint::VakintSettings {
                number_of_terms_in_epsilon_expansion: 4,
                evaluation_order,
                ..VakintSettings::default().true_settings()
            };
            let evaluate = |numerator: Atom, power| -> Result<Atom> {
                let mut term = vakint::VakintTerm {
                    integral: function!(
                        vakint::symbols::S.topo,
                        function!(vakint_symbol!("I1L"), &mass_squared, power)
                    ),
                    numerator,
                    vectors: vec![("k".into(), 1)],
                };
                term.evaluate_integral(vakint, &settings)?;
                Ok(exact_one_loop_constants(&term.numerator))
            };

            // V=-(1-epsilon) I2 (p-slash-m) precedes the forest minus.
            // The signed local counterterm -V has radial coefficient
            // (d-2)[I2-2/d (I2+M^2 I3)]. Evaluate the two denominator
            // topologies independently instead of assuming the IBP identity.
            // Evaluate each scalar topology once and reuse it for both spin
            // dimensions. Angular dimension D remains unchanged in either lift.
            let i2 = evaluate(Atom::one(), 2)?;
            let i3 = evaluate(Atom::one(), 3)?;
            let tensor_reduced = series(
                &((&dimension - Atom::num(2)) * (Atom::one() - Atom::num(2) / &dimension) * &i2
                    - Atom::num(2) * (&dimension - Atom::num(2)) * &mass_squared / &dimension
                        * &i3),
                2,
            )?
            .to_atom();
            // The quark-loop chart gives -V4 = 4/D * <k^2/Dvac^3>.
            // Its radial numerator is I2+M^2 I3; this also tests the tensor
            // integration identity without substituting a counterterm result.
            let signed_four_vacuum = series(
                &(Atom::num(4) / &dimension * (&i2 + &mass_squared * &i3)),
                2,
            )?
            .to_atom();
            assert!(
                series(&(&signed_four_vacuum - &i2), 2)?
                    .to_atom()
                    .expand()
                    .is_zero()
            );
            assert!(
                series(&(&signed_four_vacuum - &tensor_reduced - &epsilon * &i2), 2)?
                    .to_atom()
                    .expand()
                    .is_zero()
            );
            // Deliberately alter only the lifted-four angular average. The
            // formal-D tensor_reduced reference must remain unchanged.
            let wrong_angular_four = &i2 + &mass_squared * &i3;
            let wrong_bridge = truncate(
                &series(
                    &((wrong_angular_four - &tensor_reduced) / &normalization),
                    0,
                )?,
                true,
            );
            assert!(
                (wrong_bridge - Atom::num((1, 2))).expand().is_zero(),
                "freezing the D angular average early must lose half the rational bridge"
            );
            let scalar_reference = evaluate(Atom::one() - &epsilon, 2)?;
            assert!(
                series(&(&tensor_reduced - &scalar_reference), 2)?
                    .to_atom()
                    .expand()
                    .is_zero(),
                "{backend}: the signed OS vacuum remainder must equal (1-epsilon) I2 through epsilon^2"
            );

            // exp(gamma_E epsilon) Gamma(epsilon) at mu=M follows from
            // log Gamma(1+epsilon)=-gamma_E epsilon+zeta(2)epsilon^2/2
            // -zeta(3)epsilon^3/3+... . Multiplication by (1-epsilon)
            // also contributes -zeta(2)/2 to the epsilon^2 coefficient.
            let at_reference_scale = tensor_reduced
                .replace(GS.mu_r_sq)
                .with(mass_squared.clone());
            let expected = &normalization
                * (epsilon.pow(-1) - Atom::one() + pi.pow(2) * &epsilon / Atom::num(12)
                    - (pi.pow(2) / Atom::num(12) + &zeta_three / Atom::num(3)) * epsilon.pow(2));
            assert!(
                series(&(at_reference_scale - expected), 2)?
                    .to_atom()
                    .expand()
                    .is_zero(),
                "{backend}: the vacuum Gamma-function reference differs through epsilon^2"
            );

            let physical_mass = Atom::var(symbol!("os_reference_physical_mass"));
            let logarithm = -function!(
                symbolica::atom::Symbol::LOG,
                physical_mass.pow(2) / Atom::var(GS.mu_r_sq)
            );
            let mass_reference = evaluate(on_shell_weight.clone(), 2)?
                .replace(GS.m_uv_vacuum)
                .with(physical_mass);
            // O = exp(gamma_E epsilon) Gamma(epsilon) (mu^2/m^2)^epsilon
            // times the independently integrated Feynman-parameter weight.
            let expected_mass = &normalization
                * (Atom::num(3) * epsilon.pow(-1)
                    + Atom::num(3) * &logarithm
                    + Atom::num(4)
                    + &epsilon
                        * (Atom::num(3) * logarithm.pow(2) / Atom::num(2)
                            + Atom::num(4) * &logarithm
                            + Atom::num(8)
                            + pi.pow(2) / Atom::num(4))
                    + epsilon.pow(2)
                        * (logarithm.pow(3) / Atom::num(2)
                            + Atom::num(2) * logarithm.pow(2)
                            + Atom::num(8) * &logarithm
                            + Atom::num(16)
                            + &logarithm * pi.pow(2) / Atom::num(4)
                            + pi.pow(2) / Atom::num(3)
                            - &zeta_three));
            let difference = series(&(mass_reference - expected_mass), 2)?.to_atom();
            assert!(
                difference
                    .try_to_rational_polynomial::<_, _, u32>(&*Q_I, &*Q_I, None)?
                    .numerator
                    .is_zero(),
                "{backend}: the physical OS mass reference differs through epsilon^2: {difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn exact_os_coefficients_use_fundamental_color_dimension_identities() -> Result<()> {
        test_initialise()?;
        let spectator = parse_lit!((color_a + color_b) * (color_c - color_d) ^ 11);
        for dimension in [3, 5] {
            let casimir = idenso::color_cas!(2, idenso::cof!(dimension as usize));
            let index = idenso::color_idx!(2, idenso::cof!(dimension as usize));
            let dimension_ratio = Atom::num(dimension * dimension - 1) / Atom::num(dimension);
            let settings = idenso::tensor::AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            };
            let difference = casimir.pow(2) - &dimension_ratio * &casimir * index;
            assert!(
                SymbolicTensor::infer_indexed(&spectator * difference)?
                    .simplify_algebra(&settings)?
                    .into_expression()
                    .is_zero()
            );
            assert_eq!(
                SymbolicTensor::infer_indexed(&spectator * casimir)?
                    .simplify_algebra(&settings)?
                    .into_expression(),
                &spectator * dimension_ratio / Atom::num(2),
            );
        }
        Ok(())
    }

    #[test]
    fn analytic_spin_simplification_preserves_shared_color_factor() {
        use idenso::{bis, coad, cof, color_f, color_idx, color_t, gamma};
        use spenso::{mink, trace, trace_sym};

        test_initialise().unwrap();
        let a = coad!(8, color_a);
        let b = coad!(8, color_b);
        let c = coad!(8, color_c);
        let color = trace_sym!(cof!(3), color_t!(&a), color_t!(&b), color_t!(&c))
            + Atom::num(2) * color_f!(&a, &b, &c) * color_idx!(2, cof!(3));
        let spin = trace!(
            bis!(4),
            gamma!(mink!(4, mu)),
            gamma!(mink!(4, nu)),
            gamma!(mink!(4, rho)),
            gamma!(mink!(4, sigma))
        );
        let reduced_spin = simplify(&spin).unwrap();
        assert!(!reduced_spin.is_zero());
        assert!(!reduced_spin.contains_symbol(AGS.gamma));
        // Compare the literal factorization after the complete analytic spin
        // pipeline, without expanding the numerator for the assertion.
        assert_eq!(simplify(&(&color * spin)).unwrap(), color * reduced_spin);
    }

    #[test]
    fn analytic_spin_expansion_preserves_powered_scalar_products() {
        test_initialise().unwrap();
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let p = function!(GS.loop_mom, 1, mink.to_symbolic([]));
        let k = function!(GS.loop_mom, 2, mink.to_symbolic([]));
        let spin = function!(AGS.gamma, &a, &b, &mu) * function!(AGS.gamma, &b, &a, &nu);
        for product in [spenso::dot!(&p, &k), spenso::g!(&p, &k)] {
            for power in [2, 3] {
                // A closed Dirac loop makes the analytic spin-expansion path
                // active. Its scalar spectator must retain its own contractions.
                let expected =
                    Atom::num(4) * spenso::g!(&mu, &nu) * spenso::dot!(&p, &k).pow(power);
                assert_eq!(simplify(&(product.pow(power) * &spin)).unwrap(), expected);
                // These compact metrics are tensor contractions, so they must
                // remain visible even beside the protected scalar spectator.
                assert_eq!(
                    simplify(&(product.pow(power) * &spin * spenso::g!(&mu, &p))).unwrap(),
                    Atom::num(4) * function!(GS.loop_mom, 1, &nu) * spenso::dot!(&p, &k).pow(power)
                );
                assert_eq!(
                    simplify(
                        &(product.pow(power) * &spin * spenso::g!(&mu, &p) * spenso::g!(&nu, &k)),
                    )
                    .unwrap(),
                    Atom::num(4) * spenso::dot!(&p, &k).pow(power + 1)
                );
            }
        }
    }

    #[test]
    fn vakint_projectors_preserve_dirac_ports_in_both_backends() {
        test_initialise().unwrap();
        let vakint = Vakint::new().unwrap();
        let mink = Minkowski {}.new_rep(GS.dim);
        let bis = Bispinor {}.new_rep(4);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let c = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let radial = vakint::symbols::S.dot(
            function!(vakint::symbols::S.k, 1),
            function!(vakint::symbols::S.k, 2),
        );
        for tensor_reduction_method in [
            vakint::TensorReductionMethod::FeynKit,
            vakint::TensorReductionMethod::AlphaLoop,
        ] {
            for project_onto_tensor_integrals in [true, false] {
                for use_dot_product_notation in [true, false] {
                    let settings = vakint::VakintSettings {
                        tensor_reduction_method,
                        project_onto_tensor_integrals,
                        use_dot_product_notation,
                        ..VakintSettings::default().true_settings()
                    };
                    let integrated = Integrated::new(&vakint, &settings);
                    for closed in [true, false] {
                        let numerator = function!(GS.loop_mom, 1, &mu)
                            * function!(GS.loop_mom, 2, &nu)
                            * function!(AGS.gamma, &a, &b, &mu)
                            * function!(AGS.gamma, &b, if closed { &a } else { &c }, &nu);
                        let mut term = vakint::VakintTerm {
                            integral: vakint::vakint_parse!("topo(I2L(mUVsq,1,1,1))").unwrap(),
                            numerator: to_vakint_numerator(&numerator).unwrap(),
                            vectors: vec![("k".into(), 1), ("k".into(), 2)],
                        };
                        term.tensor_reduce(&vakint, &settings).unwrap();
                        let result = integrated
                            .simplify_projected_numerator(&term.numerator, 0)
                            .unwrap();
                        let spin = if closed {
                            Atom::num(4)
                        } else {
                            function!(ETS.metric, &a, &c)
                        };
                        // Projectors may retain D/D as a scalar coefficient.
                        // Compare the complete standalone Dirac identity exactly.
                        assert_eq!(
                            (&result - &spin * &radial).cancel(),
                            Atom::Zero,
                            "{tensor_reduction_method:?}, projected={project_onto_tensor_integrals}, dots={use_dot_product_notation}, closed={closed}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn projected_dirac_algebra_precedes_single_and_double_pole_expansion() {
        test_initialise().unwrap();
        let vakint = Vakint::new().unwrap();
        let epsilon = Atom::var(GS.dim_epsilon);
        let radial = vakint::symbols::S.dot(
            function!(vakint::symbols::S.k, 1),
            function!(vakint::symbols::S.k, 2),
        );
        let mink = Minkowski {}.new_rep(GS.dim);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let c = bis.to_symbolic([Aind::new_dummy().to_atom()]);

        for project_onto_tensor_integrals in [true, false] {
            for use_dot_product_notation in [false, true] {
                let settings = vakint::VakintSettings {
                    use_dot_product_notation,
                    project_onto_tensor_integrals,
                    ..VakintSettings::default().true_settings()
                };
                let integrated = Integrated::new(&vakint, &settings);
                for closed in [false, true] {
                    // Distinct loop slashes have no repeated gamma index before
                    // angular projection. Its projector creates the contraction.
                    let numerator = function!(GS.loop_mom, 1, &mu)
                        * function!(GS.loop_mom, 2, &nu)
                        * function!(AGS.gamma, &a, &b, &mu)
                        * function!(AGS.gamma, &b, if closed { &a } else { &c }, &nu);
                    let mut term = vakint::VakintTerm {
                        integral: vakint::vakint_parse!("topo(I2L(mUVsq,1,1,1))").unwrap(),
                        numerator: to_vakint_numerator(&numerator).unwrap(),
                        vectors: vec![("k".into(), 1), ("k".into(), 2)],
                    };
                    term.tensor_reduce(&vakint, &settings).unwrap();
                    let spin = if closed {
                        Atom::num(4)
                    } else {
                        function!(ETS.metric, &a, &c)
                    };
                    let late_four_dimensional = simplify(
                        &Integrated::restore_numerator(
                            Vakint::convert_to_dot_notation(&settings, term.numerator.as_view())
                                .unwrap(),
                            0,
                        )
                        .replace(GS.dim)
                        .with(4),
                    )
                    .unwrap();
                    assert!(
                        !truncate(
                            &series(
                                &(late_four_dimensional
                                    .replace(radial.to_pattern())
                                    .with(epsilon.pow(-1))
                                    - &spin * epsilon.pow(-1)),
                                1,
                            )
                            .unwrap(),
                            true,
                        )
                        .is_zero(),
                        "the old late-4D contraction must miss a nonzero finite term"
                    );
                    for evanescent_power in 0..=2 {
                        let projected = &term.numerator
                            * (Atom::var(GS.dim) - Atom::num(4)).pow(evanescent_power);
                        let resolved = integrated
                            .simplify_projected_numerator(&projected, 0)
                            .unwrap();
                        assert!(!resolved.contains_symbol(AGS.gamma));
                        assert!(!resolved.contains_symbol(SPENSO_TAG.trace));
                        assert!(!resolved.contains_symbol(GS.dim));
                        let coefficient = &spin * (Atom::num(-2) * &epsilon).pow(evanescent_power);
                        for pole in [
                            epsilon.pow(-1) + Atom::num(5) + Atom::num(7) * &epsilon,
                            Atom::num(2) * epsilon.pow(-2)
                                + Atom::num(3) * epsilon.pow(-1)
                                + Atom::num(5)
                                + Atom::num(7) * &epsilon,
                        ] {
                            // Mock only the scalar master value: the angular
                            // projection and d-dimensional Dirac algebra above use
                            // the production pipeline. This isolates the pole ×
                            // evanescent finite terms from master normalization.
                            let evaluated = resolved
                                .replace(radial.to_pattern())
                                .with(pole.to_pattern());
                            assert!(
                                series(&(evaluated - &coefficient * pole), 2)
                                    .unwrap()
                                    .to_atom()
                                    .is_zero(),
                                "project={project_onto_tensor_integrals}, closed={closed}, dot={use_dot_product_notation}, evanescent power={evanescent_power}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn projected_dirac_numerator_matches_scalar_master_before_backend_truncation() {
        test_initialise().unwrap();
        let vakint = Vakint::new().unwrap();
        let mut results = Vec::new();
        for project_onto_tensor_integrals in [true, false] {
            let settings = vakint::VakintSettings {
                // A two-loop vacuum master is retained only through its finite
                // term, so omitting the numerator's evanescent factors fails here.
                number_of_terms_in_epsilon_expansion: 3,
                project_onto_tensor_integrals,
                ..VakintSettings::default().true_settings()
            };
            let integrated = Integrated::new(&vakint, &settings);
            let mink = Minkowski {}.new_rep(GS.dim);
            let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
            let nu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
            let bis = Bispinor {}.new_rep(4);
            let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
            let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
            let numerator = function!(GS.loop_mom, 1, &mu)
                * function!(GS.loop_mom, 2, &nu)
                * function!(AGS.gamma, &a, &b, &mu)
                * function!(AGS.gamma, &b, &a, &nu);
            let mut projected = vakint::VakintTerm {
                integral: vakint::vakint_parse!("topo(I2L(mUVsq,1,1,1))").unwrap(),
                numerator: to_vakint_numerator(&numerator).unwrap(),
                vectors: vec![("k".into(), 1), ("k".into(), 2)],
            };
            projected.tensor_reduce(&vakint, &settings).unwrap();
            let d_minus_four = Atom::var(GS.dim) - Atom::num(4);
            projected.numerator = integrated
                .simplify_projected_numerator(
                    &(&projected.numerator * (&d_minus_four + d_minus_four.pow(2))),
                    0,
                )
                .unwrap();
            let epsilon = Atom::var(GS.dim_epsilon);
            let mut scalar = vakint::VakintTerm {
                numerator: Atom::num(4)
                    * (Atom::num(-2) * &epsilon + Atom::num(4) * epsilon.pow(2))
                    * vakint::symbols::S.dot(
                        function!(vakint::symbols::S.k, 1),
                        function!(vakint::symbols::S.k, 2),
                    ),
                ..projected.clone()
            };
            let resolved_input = projected.numerator.clone();
            // The projector contains 1/(4-2ε) before symbolic-d gamma
            // contraction. Cancel the analytic oracle's rational denominators
            // exactly, using the unfactorized polynomial representation.
            assert!(
                (&resolved_input - &scalar.numerator)
                    .try_to_rational_polynomial::<_, _, u32>(&*Q_I, &*Q_I, None)
                    .expect("the analytic scalar input difference is rational over Q(i)")
                    .numerator
                    .is_zero(),
                "project={project_onto_tensor_integrals}: resolved input {resolved_input} differs from scalar oracle {}",
                scalar.numerator
            );
            projected.evaluate_integral(&vakint, &settings).unwrap();
            scalar.evaluate_integral(&vakint, &settings).unwrap();
            assert!(
                !truncate(&series(&scalar.numerator, 1).unwrap(), true).is_zero(),
                "the independently integrated evanescent scalar master has a nonzero finite term"
            );
            let difference = series(&(&projected.numerator - &scalar.numerator), 1)
                .unwrap()
                .to_atom();
            // These are analytically integrated subgraph coefficients. Expand
            // their difference to compare exact Laurent coefficients, including
            // algebraically equal backend output with different factorization.
            assert!(
                difference.expand().is_zero(),
                "project={project_onto_tensor_integrals}: resolved input {resolved_input}; projected result {}; scalar result {}; difference {difference}; expanded difference {}",
                projected.numerator,
                scalar.numerator,
                difference.expand()
            );
            results.push(projected.numerator);
        }
        assert!(
            series(&(&results[0] - &results[1]), 1)
                .unwrap()
                .to_atom()
                .expand()
                .is_zero()
        );
    }

    #[test]
    fn restored_numeric_source_ports_remain_distinct_from_backend_dummies() -> Result<()> {
        test_initialise()?;
        let source = crate::numerator::AppliedFeynmanRule::new(parse!(
            "spenso::g(spenso::mink(gammalooprs::dim,0),spenso::mink(gammalooprs::dim,1))"
        ))?
        .expr;
        let (protected, originals) = Integrated::protect_numeric_ports(&source)?;
        let labels = protected
            .structure()
            .logical_slots()
            .into_iter()
            .map(|slot| slot.aind.to_atom())
            .collect::<Vec<_>>();
        let dummy = vakint::symbols::S.dot_dummy_ind(0);
        let backend = function!(vakint::symbols::S.g, &labels[0], &dummy)
            * function!(vakint::symbols::S.g, &dummy, &labels[1]);
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        for backend in [backend, Atom::Zero] {
            let is_zero = backend.is_zero();
            let restored = Integrated::restore_numerator(backend, 19);
            let result = protected
                .with_checked_expression(cooking.try_cook(restored.as_view())?)?
                .contract(ContractSettings::default())?
                .reindex_interface_ports(&originals)?;
            assert_eq!(
                result.structure().logical_slots(),
                source.structure().logical_slots()
            );
            assert_eq!(result.rank(), 2);
            assert_eq!(result.is_zero(), is_zero);
            if !is_zero {
                assert_eq!(result.expression(), source.expression());
            }
        }
        Ok(())
    }

    #[test]
    fn integrated_backend_preserves_numeric_open_ports_in_both_projection_modes() -> Result<()> {
        use crate::{
            finalized_runtime_dot,
            graph::Graph,
            numerator::AppliedFeynmanRule,
            uv::{Spinney, UVgenerationSettings, hedge_poset::OwnedForestNode},
        };
        use linnet::half_edge::subgraph::InternalSubGraph;

        test_initialise()?;
        let graph: Graph = finalized_runtime_dot!(digraph numeric_open_vacuum {
            graph [projector=1]
            edge [particle=scalar_1 num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"]
        }, "scalars")?;
        let full = graph.full_filter();
        let current = OwnedForestNode {
            spinney: Spinney::with_scheme(
                InternalSubGraph::cleaned_filter_optimist(full.clone(), graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
                ApproximationType::MUV,
                0,
            )?
            .expect("the tadpole has its full compatible loop basis"),
            topo_order: 23,
        };
        let given = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let denominator =
            AppliedFeynmanRule::new(graph.denominator(&full, &graph.model, |_| -3))?.expr;
        let numerator = AppliedFeynmanRule::new(parse!(
            "gammalooprs::Q(0,spenso::mink(gammalooprs::dim,0))*gammalooprs::Q(0,spenso::mink(gammalooprs::dim,1))"
        ))?.expr;
        let source = numerator.multiply(denominator)?;
        let metric = AppliedFeynmanRule::new(parse!(
            "spenso::g(spenso::mink(gammalooprs::dim,0),spenso::mink(gammalooprs::dim,1))"
        ))?
        .expr;
        let scalar_source = source
            .clone()
            .multiply(metric.clone())?
            .contract(ContractSettings::default())?;
        assert!(scalar_source.is_scalar());
        let settings = UVgenerationSettings {
            add_marker: false,
            ..Default::default()
        };
        let ctx = UVCtx::new(&graph, &graph.model, &settings);
        let sector = FourDSector::new(
            source.clone(),
            vec![(full.clone(), full, graph.loop_momentum_basis.clone())],
            Vec::new(),
        );
        for project_onto_tensor_integrals in [false, true] {
            let backend_settings = vakint::VakintSettings {
                project_onto_tensor_integrals,
                // This logarithmically divergent integral starts at epsilon^-1.
                // Four coefficients supply the epsilon^2 order checked below,
                // including the evanescent term from the closing D-dimensional metric.
                number_of_terms_in_epsilon_expansion: 4,
                ..settings.vakint.true_settings()
            };
            let integrated = Integrated::new(crate::utils::vakint()?, &backend_settings);
            let result = integrated.integrate(&source, &sector, &ctx, &current, &given)?;
            assert_eq!(result.rank(), 2);
            assert_eq!(
                result.structure().logical_slots(),
                source.structure().logical_slots()
            );
            assert!(!result.is_zero());
            let scalar = integrated.integrate(&scalar_source, &sector, &ctx, &current, &given)?;
            assert!(scalar.is_scalar());
            assert!(!scalar.is_zero());
            let closed = result
                .multiply(metric.clone())?
                .contract(ContractSettings::default())?;
            assert!(closed.is_scalar());
            // Tensor metadata stays in symbolic D. Only this completed scalar
            // contraction uses D=4-2*epsilon, including every evanescent term.
            let residual =
                Integrated::dimensionally_regularized(closed.expression()) - scalar.expression();
            assert!(
                series(&residual, 2)?.to_atom().expand().is_zero(),
                "project={project_onto_tensor_integrals}: {residual}"
            );
            let zero = source.with_checked_expression(Atom::Zero)?;
            let zero = integrated.integrate(&zero, &sector, &ctx, &current, &given)?;
            assert!(zero.is_zero());
            assert_eq!(
                zero.structure().logical_slots(),
                source.structure().logical_slots()
            );
        }
        Ok(())
    }

    #[test]
    fn restored_symbolic_vakint_slots_preserve_the_native_source_interface() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spectator = parse_lit!((restore_a + restore_b) ^ 11);
        for label in [
            Aind::new_dummy(),
            Aind::Vertex(2, 3),
            Aind::UVTerm(5, 7),
            Aind::Symbol(symbol!("restored_source_slot")),
        ] {
            let index = label.to_atom();
            let other = Aind::new_dummy().to_atom();
            let left = mink.to_symbolic([index.clone()]);
            let right = mink.to_symbolic([other.clone()]);
            for represented_left in [false, true] {
                let source_left = if represented_left { &left } else { &index };
                let backend = &spectator
                    * (function!(vakint::symbols::S.g, source_left, &other)
                        + function!(vakint::symbols::S.k, 1, source_left)
                            * function!(vakint::symbols::S.p, 0, &other));
                let expected = &spectator
                    * (function!(ETS.metric, &left, &right)
                        + function!(GS.loop_mom, 1, &left) * function!(GS.emr_mom, 0, &right));
                let restored = Integrated::restore_numerator(backend, 9);
                assert_eq!(
                    restored, expected,
                    "label={label:?}, represented={represented_left}"
                );
                let canonical_backend = to_vakint_numerator(&expected)?;
                assert_eq!(to_vakint_numerator(&restored)?, canonical_backend);
                assert_eq!(
                    Integrated::restore_numerator(canonical_backend, 9),
                    expected
                );
            }
        }
        Ok(())
    }

    #[test]
    fn restored_projector_indices_are_fresh_after_dimension_unification() {
        test_initialise().unwrap();
        let Aind::Dummy(previous) = Aind::new_dummy() else {
            unreachable!()
        };
        let reserved_index = Aind::new_dummy_at(previous + 1).to_atom();
        let native = Minkowski {}.new_rep(GS.dim).to_symbolic([reserved_index]);
        let regulated = native
            .replace(GS.dim)
            .with(Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon));
        let first = spenso::tensor_symbol!("reserved_projector_port");
        let second = spenso::tensor_symbol!("fresh_projector_port");
        let input = function!(first, regulated)
            * function!(second, function!(vakint::symbols::S.tensor_index, 1, 0));
        let restored = Integrated::restore_numerator(input, 0);
        assert!(!restored.contains_symbol(vakint::symbols::S.tensor_index));
        assert!(
            restored
                .pattern_match(&function!(first, &native).to_pattern(), None, None)
                .next()
                .is_some()
        );
        assert!(
            restored
                .pattern_match(&function!(second, &native).to_pattern(), None, None)
                .next()
                .is_none()
        );
    }

    #[test]
    fn scalar_dimensions_are_substituted_inside_functions_but_not_lorentz_slots() {
        test_initialise().unwrap();
        let slot = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let logarithm = symbolica::symbol!("log");
        let tensor = symbolica::symbol!("dimension_test_tensor");
        let input = function!(logarithm, GS.dim) * function!(tensor, &slot);
        let expected = function!(
            logarithm,
            Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon)
        ) * function!(tensor, &slot);
        assert_eq!(Integrated::dimensionally_regularized(&input), expected);
        let dimension = Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon);
        let projected_slot = slot.replace(GS.dim).with(dimension);
        let evanescent = Atom::var(GS.dim) - Atom::num(4);
        assert_eq!(
            Integrated::restore_numerator(&evanescent * function!(tensor, projected_slot), 0),
            evanescent * function!(tensor, &slot),
            "restoration unifies Lorentz metadata without changing scalar d dependence"
        );
    }

    #[test]
    fn analytic_uv_spin_algebra_rejects_unsupported_residual_contractions() {
        test_initialise().unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let b = bis.to_symbolic([Aind::new_dummy().to_atom()]);
        let mu = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let nu = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let unknown = function!(
            spenso::tensor_symbol!("unsupported_dirac_factor"),
            SPENSO_TAG.chain_in,
            SPENSO_TAG.chain_out
        );
        let gamma = function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, &mu);
        for (input, message) in [
            (function!(AGS.sigma, &a, &b, &mu, &nu), "Sigma tensors"),
            (
                function!(SPENSO_TAG.trace, bis.to_symbolic([]), &unknown),
                "unresolved Dirac trace",
            ),
            (
                function!(
                    SPENSO_TAG.trace,
                    Minkowski {}.new_rep(GS.dim).to_symbolic([]),
                    &unknown
                ),
                "unresolved Lorentz trace",
            ),
        ] {
            assert!(simplify(&input).unwrap_err().to_string().contains(message));
        }
        let residual = function!(SPENSO_TAG.chain, &a, &b, &gamma, &unknown, &gamma);
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking
                    .try_cook(simplify(&residual).unwrap().as_view())
                    .unwrap()
            )
            .is_err()
        );
        // An open gamma is a valid d-dimensional tensor basis element.
        let open = simplify(&function!(AGS.gamma, &a, &b, &mu)).unwrap();
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(open.as_view()).unwrap()
            )
            .is_ok()
        );
    }

    #[test]
    fn restored_vakint_dots_have_ports_and_roundtrip_without_distribution() {
        test_initialise().unwrap();
        let spectator = (Atom::var(symbol!("restore_dot::x")) + 1).pow(5);
        for left in [vakint::symbols::S.k, vakint::symbols::S.p] {
            for right in [vakint::symbols::S.k, vakint::symbols::S.p] {
                let dot = function!(vakint::symbols::S.dot, left.call(1), right.call(2));
                let source = &spectator * dot.pow(2);
                let restored = Integrated::restore_numerator(source.clone(), 0);
                let admitted = SymbolicTensor::infer_indexed(restored.clone()).unwrap();
                assert_eq!(admitted.expression(), &restored);
                assert_eq!(to_vakint_numerator(&restored).unwrap(), source);
                assert_eq!(
                    to_vakint_numerator(&simplify(&restored).unwrap()).unwrap(),
                    source
                );
            }
        }
    }

    #[test]
    fn residual_lorentz_check_preserves_native_color_dummy_indices() -> Result<()> {
        test_initialise()?;
        let adjoint = idenso::representations::ColorAdjoint {}.new_rep(8);
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        for native_dummy_names in [false, true] {
            let [a, b, c, x, y] = [0, 1, 2, 1_000_000, 1_000_001].map(|index| {
                adjoint.to_symbolic([if native_dummy_names {
                    AbstractIndex::new_dummy_at(index).to_atom()
                } else {
                    Aind::Dummy(index).to_atom()
                }])
            });
            let word = spenso::chain!(
                &a,
                &b,
                function!(
                    idenso::color::CS.f,
                    SPENSO_TAG.chain_in,
                    SPENSO_TAG.chain_out,
                    &x
                ),
                function!(
                    idenso::color::CS.f,
                    SPENSO_TAG.chain_in,
                    SPENSO_TAG.chain_out,
                    &y
                )
            );
            let momentum = GS.emr_mom(
                EdgeIndex(0),
                Minkowski {}
                    .new_rep(GS.dim)
                    .to_symbolic([Aind::Hedge(9, 1).to_atom()]),
            );
            let source = word * idenso::color_f!(&x, &y, &c) * momentum;
            let encoded = cooking.try_cook(source.as_view())?;
            let native = SymbolicTensor::infer_indexed(encoded.clone())?;
            Integrated::ensure_resolved_lorentz_contractions(&encoded)?;
            Integrated::ensure_resolved_lorentz_contractions(native.expression())?;
            assert_eq!(cooking.uncook(native.expression().as_view()), source);

            // An actual third occurrence still fails after reserving native and
            // encoded source dummies. Only invented collisions are prevented.
            let extra = idenso::color_f!(
                &x,
                adjoint.to_symbolic([Aind::Hedge(10, 1).to_atom()]),
                adjoint.to_symbolic([Aind::Hedge(11, 1).to_atom()])
            );
            let malformed = native.expression() * cooking.try_cook(extra.as_view())?;
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(&malformed)
                    .unwrap_err()
                    .to_string()
                    .contains("more than two compatible tensor ports")
            );
        }
        Ok(())
    }

    #[test]
    fn residual_tensor_contractions_are_checked_per_analytic_branch() {
        test_initialise().unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let [a, b, c, d] =
            std::array::from_fn::<_, 4, _>(|_| bis.to_symbolic([Aind::new_dummy().to_atom()]));
        let mink = Minkowski {}.new_rep(GS.dim);
        let [mu, nu, rho] =
            std::array::from_fn::<_, 3, _>(|_| mink.to_symbolic([Aind::new_dummy().to_atom()]));
        let left = function!(AGS.gamma, &a, &b, &mu);
        let right = function!(AGS.gamma, &c, &d, &mu);
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook((&left * &right).as_view()).unwrap()
            )
            .is_err()
        );
        // Compact metrics contract Lorentz axes while retaining spinor ports;
        // dot is reserved for rank-one operands in the typed tensor API.
        let compact = function!(
            ETS.metric,
            function!(AGS.gamma, &a, &b, mink.to_symbolic([])),
            function!(AGS.gamma, &c, &d, mink.to_symbolic([]))
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(compact.as_view()).unwrap()
            )
            .is_err()
        );
        let same_chain = function!(
            ETS.metric,
            function!(AGS.gamma, &a, &b, mink.to_symbolic([])),
            function!(AGS.gamma, &b, &c, mink.to_symbolic([]))
        );
        let same_chain = simplify(&same_chain).unwrap();
        assert!(
            (&same_chain - Atom::var(GS.dim) * function!(ETS.metric, &a, &c))
                .expand()
                .is_zero(),
            "a compact contraction within one gamma chain must close to d times its spinor identity: {same_chain}"
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(same_chain.as_view()).unwrap()
            )
            .is_ok()
        );
        let free =
            &left * function!(AGS.gamma, &c, &d, &nu) + function!(AGS.gamma, &a, &b, &nu) * &right;
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(free.as_view()).unwrap()
            )
            .is_ok()
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking
                    .try_cook(function!(ETS.metric, &mu, &nu).as_view())
                    .unwrap()
            )
            .is_ok()
        );
        let ambiguous = function!(
            ETS.metric,
            function!(
                AGS.sigma,
                &a,
                &b,
                mink.to_symbolic([]),
                mink.to_symbolic([])
            ),
            function!(
                AGS.sigma,
                &c,
                &d,
                mink.to_symbolic([]),
                mink.to_symbolic([])
            )
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(ambiguous.as_view()).unwrap()
            )
            .is_err()
        );
        let slash = function!(
            AGS.gamma,
            &a,
            &b,
            function!(GS.emr_mom, 0, mink.to_symbolic([]))
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(slash.as_view()).unwrap()
            )
            .is_ok()
        );
        let slash_product = &slash
            * function!(
                AGS.gamma,
                &c,
                &d,
                function!(GS.emr_mom, 1, mink.to_symbolic([]))
            );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(slash_product.as_view()).unwrap()
            )
            .is_ok()
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking
                    .try_cook((&left / (Atom::var(GS.dim) - Atom::one())).as_view())
                    .unwrap()
            )
            .is_ok()
        );

        // Diagnostic tensor operator, not a physical graph/locality oracle:
        // crossed spinor closure of two contracted triple-gamma chains gives
        // 4d(-d²+6d-4). Its evanescent term changes the finite part of a pole.
        let factors = [&mu, &nu, &rho]
            .map(|index| function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, index));
        let operator = function!(
            SPENSO_TAG.chain,
            &a,
            &b,
            &factors[0],
            &factors[1],
            &factors[2]
        ) * function!(
            SPENSO_TAG.chain,
            &c,
            &d,
            &factors[0],
            &factors[1],
            &factors[2]
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(operator.as_view()).unwrap()
            )
            .is_err()
        );
        let closed = function!(
            SPENSO_TAG.trace,
            bis.to_symbolic([]),
            &factors[0],
            &factors[1],
            &factors[2],
            &factors[0],
            &factors[1],
            &factors[2]
        );
        let closed = Integrated::dimensionally_regularized(&simplify(&closed).unwrap());
        let epsilon = Atom::var(GS.dim_epsilon);
        assert_eq!(
            truncate(&series(&(closed / epsilon), 0).unwrap(), true),
            Atom::num(32)
        );
    }

    #[test]
    fn residual_lorentz_diagnostic_keeps_factored_weighted_branches() {
        test_initialise().unwrap();
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let [a, b, c, d] =
            std::array::from_fn::<_, 4, _>(|_| bis.to_symbolic([Aind::new_dummy().to_atom()]));
        let mink = Minkowski {}.new_rep(GS.dim);
        let [mu, nu] =
            std::array::from_fn::<_, 2, _>(|_| mink.to_symbolic([Aind::new_dummy().to_atom()]));
        let left = function!(AGS.gamma, &a, &b, &mu);
        let right = function!(AGS.gamma, &c, &d, &mu);
        let metadata = function!(
            symbolica::symbol!("residual_lorentz_metadata"; Scalar),
            function!(
                SPENSO_TAG.dot,
                function!(AGS.gamma, &a, &b, mink.to_symbolic([])),
                function!(AGS.gamma, &c, &d, mink.to_symbolic([]))
            )
        );
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(metadata.as_view()).unwrap()
            )
            .is_ok()
        );
        let coefficient = parse!("(residual_x-residual_y)*(residual_x+residual_y)^30");
        let alternative = function!(
            SPENSO_TAG.tensor_symbol("residual_lorentz_alternative"),
            &a,
            &b,
            &mu
        );
        let factored = &coefficient * (&left + alternative) * &right;
        let original = factored.clone();
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(factored.as_view()).unwrap()
            )
            .unwrap_err()
            .to_string()
            .contains("unresolved internal Lorentz tensor contraction")
        );
        assert_eq!(factored, original);
        let free = &coefficient
            * (&left * function!(AGS.gamma, &c, &d, &nu)
                + function!(AGS.gamma, &a, &b, &nu) * &right);
        let original = free.clone();
        assert!(
            Integrated::ensure_resolved_lorentz_contractions(
                &cooking.try_cook(free.as_view()).unwrap()
            )
            .is_ok()
        );
        assert_eq!(free, original);
    }

    #[test]
    fn integrated_counterterm_projects_one_laurent_expansion() {
        test_initialise().unwrap();

        let epsilon = Atom::var(GS.dim_epsilon);
        let expansion = series(
            &(Atom::num(2) * epsilon.pow(-2)
                + Atom::num(3) * epsilon.pow(-1)
                + Atom::num(5)
                + Atom::num(7) * &epsilon),
            2,
        )
        .unwrap();
        let source = SymbolicTensor::infer_indexed(expansion.to_atom()).unwrap();
        let integrated = IntegratedCts::from_series(&source, &expansion, 4).unwrap();
        let scale = Atom::var(GS.integrated_loop_scale).pow(4);

        assert_eq!(
            integrated.pole().unwrap().expression().expand(),
            ((Atom::num(2) * epsilon.pow(-2) + Atom::num(3) * epsilon.pow(-1)) * &scale).expand()
        );
        assert_eq!(
            integrated
                .finite_counterterm()
                .unwrap()
                .expression()
                .expand(),
            (-(Atom::num(5) + Atom::num(7) * epsilon) * scale).expand()
        );
        assert_eq!(
            integrated
                .physical_finite_counterterm()
                .unwrap()
                .expression(),
            &(-(Atom::num(5) + Atom::num(7) * Atom::var(GS.dim_epsilon)))
        );
    }

    #[test]
    fn finite_projection_retains_child_epsilon_terms_until_parent_poles_multiply() {
        test_initialise().unwrap();
        let epsilon = Atom::var(GS.dim_epsilon);
        let child_source = SymbolicTensor::infer_indexed(
            Atom::num(2) + Atom::num(3) * &epsilon + Atom::num(5) * epsilon.pow(2),
        )
        .unwrap();
        let child = IntegratedCts::from_series(
            &child_source,
            &series(child_source.expression(), 3).unwrap(),
            4,
        )
        .unwrap();
        let parent_poles =
            SymbolicTensor::infer_indexed(epsilon.pow(-2) + Atom::num(7) * epsilon.pow(-1))
                .unwrap();
        let nested_source = child
            .physical_finite_counterterm()
            .unwrap()
            .multiply(parent_poles)
            .unwrap();
        let nested = IntegratedCts::from_series(
            &nested_source,
            &series(nested_source.expression(), 1).unwrap(),
            8,
        )
        .unwrap();
        // The finite projection retains nonnegative powers for further nesting.
        let finite = nested.physical_finite_counterterm().unwrap();
        assert_eq!(
            finite.expression().replace(GS.dim_epsilon).with(Atom::Zero),
            Atom::num(26)
        );
        assert!(
            (finite.expression() - Atom::num(26) - Atom::num(35) * &epsilon)
                .expand()
                .is_zero()
        );
        assert_eq!(
            nested.physical_pole().unwrap().expression(),
            &(-Atom::num(2) * epsilon.pow(-2) - Atom::num(17) * epsilon.pow(-1)),
        );
    }

    #[test]
    fn factorized_product_projects_each_component() {
        test_initialise().unwrap();

        let epsilon = Atom::var(GS.dim_epsilon);
        let factors = [2, 3, 5].map(|finite| {
            let source =
                SymbolicTensor::infer_indexed(epsilon.pow(-1) + Atom::num(finite)).unwrap();
            IntegratedCts::from_series(&source, &series(source.expression(), 1).unwrap(), 0)
                .unwrap()
        });
        let product = IntegratedCts::factorized_product(&factors[..2], 1).unwrap();

        assert_eq!(
            product.physical_pole().unwrap().expression(),
            &epsilon.pow(-2)
        );
        assert_eq!(
            product.physical_finite_counterterm().unwrap().expression(),
            &Atom::num(6)
        );

        let product = IntegratedCts::factorized_product(&factors, 1).unwrap();
        assert_eq!(
            product.physical_pole().unwrap().expression(),
            &epsilon.pow(-3)
        );
        assert_eq!(
            product.physical_finite_counterterm().unwrap().expression(),
            &Atom::num(-30)
        );
    }

    // #[test]
    // fn integrated_triangle_norm_is_euclidean() {
    //     test_initialise().unwrap();

    //     let edge = EdgeIndex(7);
    //     let euclidean_norm = GS.emr_mom(edge, GS.cind(1)).pow(2)
    //         + GS.emr_mom(edge, GS.cind(2)).pow(2)
    //         + GS.emr_mom(edge, GS.cind(3)).pow(2);
    //     let minkowski = Minkowski {}.new_rep(4).to_symbolic([]);
    //     let minkowski_norm = function!(
    //         SPENSO_TAG.dot,
    //         GS.emr_vec(edge, minkowski.as_view()),
    //         GS.emr_vec(edge, minkowski.as_view())
    //     );

    //     assert_eq!(
    //         euclidean_norm,
    //         GS.emr_mom(edge, GS.cind(1)).pow(2)
    //             + GS.emr_mom(edge, GS.cind(2)).pow(2)
    //             + GS.emr_mom(edge, GS.cind(3)).pow(2)
    //     );
    //     assert_ne!(euclidean_norm, minkowski_norm);
    // }

    #[test]
    fn vakint_dot_conversion_keeps_loop_momentum_tagged_until_to_dots() {
        test_initialise().unwrap();

        let mink = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let numerator = function!(
            ETS.metric,
            function!(GS.emr_mom, 0, mink.clone()),
            function!(GS.loop_mom, 1, mink.clone())
        );

        let converted = to_vakint_numerator(&numerator).unwrap();

        assert_eq!(
            converted,
            vakint::symbols::S.dot(
                function!(vakint::symbols::S.p, 0),
                function!(vakint::symbols::S.k, 1)
            )
        );
    }

    #[test]
    fn affine_vakint_numerator_preserves_explicit_pairs_and_factorization() {
        test_initialise().unwrap();
        let q = function!(GS.emr_mom, 0);
        let p = function!(GS.emr_mom, 9);
        let k = function!(GS.loop_mom, 1);
        let add_arguments = [
            Replacement::new(
                function!(GS.emr_mom, W_.i_).to_pattern(),
                function!(GS.emr_mom, W_.i_, W_.a___),
            )
            .allow_new_wildcards_on_rhs(true),
            Replacement::new(
                function!(GS.loop_mom, W_.i_).to_pattern(),
                function!(GS.loop_mom, W_.i_, W_.a___),
            )
            .allow_new_wildcards_on_rhs(true),
        ];
        let solution =
            VakintMomentumSolution::solve(&[&q - &k - &p], &[q], &add_arguments).unwrap();
        let mu = Minkowski {}
            .new_rep(GS.dim)
            .to_symbolic([Aind::new_dummy().to_atom()]);
        let q = function!(GS.emr_mom, 0, &mu);
        let p = function!(GS.emr_mom, 9, &mu);
        let k = function!(GS.loop_mom, 1, &mu);
        let spectator = parse_lit!((affine_metric_x + affine_metric_y) ^ 30);
        for power in [1, 2, -1] {
            let numerator = &spectator * (&q * &p + Atom::one()).pow(power);
            let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
            let prepared =
                SymbolicTensor::infer_indexed(cooking.try_cook(numerator.as_view()).unwrap())
                    .unwrap()
                    .contract(ContractSettings {
                        rank_one: false,
                        collect_chains: false,
                        collect_traces: false,
                        ..Default::default()
                    })
                    .unwrap();
            let rewritten = solution.rewrite_numerator(&prepared).unwrap();
            assert_eq!(
                cooking.uncook(rewritten.expression().as_view()),
                &spectator * ((&k + &p) * &p + Atom::one()).pow(power)
            );
        }
    }

    #[test]
    fn underdetermined_vakint_momentum_solve_tracks_free_variables() {
        test_initialise().unwrap();

        let q0 = function!(GS.emr_mom, 0);
        let q1 = function!(GS.emr_mom, 1);
        let q2 = function!(GS.emr_mom, 2);
        let k0 = function!(GS.loop_mom, 0);
        let k1 = function!(GS.loop_mom, 1);
        let edge_3 = -&q0 - &q1 - &q2;
        let system = vec![&q0 - &k0, &edge_3 - &k1];
        let variables = vec![q0.clone(), q1.clone(), q2.clone()];

        let solution = VakintMomentumSolution::solve(&system, &variables, &[]).unwrap();

        assert_eq!(solution.free_variables, vec![q2]);
        assert!(
            solution
                .ensure_free_variables_eliminated(0, "numerator", &solution.free_variables[0], &[])
                .is_err()
        );
    }

    #[test]
    fn underdetermined_vakint_momentum_solve_projects_topology() {
        test_initialise().unwrap();

        let q1 = function!(GS.emr_mom, 1);
        let q2 = function!(GS.emr_mom, 2);
        let k0 = function!(GS.loop_mom, 0);
        let system = vec![-&q1 - &q2 - &k0];
        let variables = vec![q1.clone(), q2.clone()];
        let topology = function!(
            vakint::symbols::S.topo,
            function!(
                vakint::symbols::S.prop,
                1,
                function!(vakint::symbols::S.edge, 0, 0),
                -&q1 - &q2,
                GS.m_uv_vacuum,
                1
            )
        );

        let solution = VakintMomentumSolution::solve(&system, &variables, &[]).unwrap();

        solution
            .ensure_free_variables_eliminated(0, "numerator", &Atom::from(1), &[])
            .unwrap();
        assert_eq!(
            solution.rewrite_integral(&topology, &[]),
            function!(
                vakint::symbols::S.topo,
                function!(
                    vakint::symbols::S.prop,
                    1,
                    function!(vakint::symbols::S.edge, 0, 0),
                    k0,
                    GS.m_uv_vacuum,
                    1
                )
            )
        );
    }

    #[test]
    fn empty_vakint_momentum_solve_treats_variables_as_free() {
        test_initialise().unwrap();

        let q0 = function!(GS.emr_mom, 0);
        let no_variables = VakintMomentumSolution::solve(&[], &[], &[]).unwrap();
        let free_variable =
            VakintMomentumSolution::solve(&[], std::slice::from_ref(&q0), &[]).unwrap();

        assert!(no_variables.free_variables.is_empty());
        assert_eq!(free_variable.free_variables, vec![q0]);
    }

    #[test]
    fn vakint_preparation_reduces_color_before_opening_traces() {
        use crate::{finalized_runtime_dot, graph::Graph};
        use idenso::{coad, cof, color_cas, color_f, color_idx, color_t};
        use spenso::{g, trace_sym};

        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph color_vakint_tadpole {
            projector=1
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> a [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"]
            a -> outgoing [id=2 source="{ufo_order:3}"]
        })
        .unwrap();
        let [a, b, c, d] = [
            "color_vakint_a",
            "color_vakint_b",
            "color_vakint_c",
            "color_vakint_d",
        ]
        .map(|name| coad!(8, Atom::var(symbol!(name))));
        let trace = trace_sym!(cof!(3), color_t!(&a), color_t!(&b), color_t!(&c));
        let index = color_idx!(2, cof!(3));
        let color = color_f!(&d, &b, &c) * (trace + Atom::i() / 2 * color_f!(&a, &b, &c) * &index);
        let coefficient = parse!("(color_x+color_y)^3*(color_u+color_v)");
        let denominator = function!(
            GS.den,
            1,
            function!(GS.emr_mom, 1),
            Atom::var(GS.m_uv_vacuum).pow(2)
        );
        let actual = to_vakint_integrand(
            &(color * &coefficient / denominator.pow(2)),
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                graph.loop_momentum_basis.clone(),
            )],
            &VakintSettings {
                additional_normalization: "1".to_string(),
                ..Default::default()
            },
            true,
        )
        .unwrap();

        // The symmetric trace annihilates against f; the remaining f*f
        // contraction is C_A times the adjoint metric, independent of kinematics.
        let expected = Atom::i() / 2 * color_cas!(2, coad!(8)) * index * g!(&a, &d) * coefficient;
        assert_eq!(actual.0.len(), 1);
        assert_eq!(actual.0[0].numerator, expected);
    }

    #[test]
    fn vakint_additional_normalization_requires_a_scalar_interface() {
        use crate::{finalized_runtime_dot, graph::Graph};

        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph normalization_tadpole {
            projector=1
            edge [num=1 mass=2]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> a [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"]
            a -> outgoing [id=2 source="{ufo_order:3}"]
        })
        .unwrap();
        let momentum = function!(GS.emr_mom, 1);
        let vector = function!(GS.emr_mom, 1, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let denominator = function!(
            GS.den,
            1,
            momentum,
            4,
            function!(SPENSO_TAG.dot, &vector, &vector) - 4
        );
        let run = |normalization: &str| {
            to_vakint_integrand(
                &denominator.pow(-3),
                &graph,
                &graph.full_filter(),
                &graph.empty_subgraph::<SuBitGraph>(),
                &[(
                    graph.full_filter(),
                    graph.full_filter(),
                    graph.loop_momentum_basis.clone(),
                )],
                &VakintSettings {
                    additional_normalization: normalization.to_owned(),
                    ..Default::default()
                },
                false,
            )
        };
        let baseline = run("1").unwrap();
        let normalized = run("3").unwrap();
        assert_eq!(baseline.0.len(), 1);
        assert_eq!(normalized.0.len(), 1);
        assert_eq!(normalized.0[0].integral, baseline.0[0].integral);
        assert_eq!(normalized.0[0].numerator, 3 * &baseline.0[0].numerator);
        let error = run("spenso::g(spenso::mink(4,75101),spenso::mink(4,75103))")
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("additional_normalization must have a scalar tensor interface"),
            "{error}"
        );
    }

    #[test]
    fn vakint_conversion_preserves_squared_masses() {
        use crate::{finalized_runtime_dot, graph::Graph};

        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph vakint_mass_tadpole { graph [projector=1]
            edge [num=1 mass=2]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> a [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"]
            a -> outgoing [id=2 source="{ufo_order:3}"]
        })
        .unwrap();
        let momentum = function!(GS.emr_mom, 1);
        let indexed_momentum =
            function!(GS.emr_mom, 1, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
        let settings = VakintSettings {
            additional_normalization: "1".to_string(),
            ..Default::default()
        };
        // Both wrappers carry squared masses. A nonunit physical mass detects
        // an accidental second square, while the UV substitution remains MUV^2.
        for mass in [Atom::num(2), parse!("m_phys")] {
            let mass_squared = mass.pow(2);
            let denominator = function!(
                GS.den,
                1,
                &momentum,
                &mass_squared,
                function!(SPENSO_TAG.dot, &indexed_momentum, &indexed_momentum) - &mass_squared
            );
            for substitute_masses_to_m_uv in [false, true] {
                let actual = to_vakint_integrand(
                    &denominator.pow(-3),
                    &graph,
                    &graph.full_filter(),
                    &graph.empty_subgraph::<SuBitGraph>(),
                    &[(
                        graph.full_filter(),
                        graph.full_filter(),
                        graph.loop_momentum_basis.clone(),
                    )],
                    &settings,
                    substitute_masses_to_m_uv,
                )
                .unwrap();
                assert_eq!(actual.0.len(), 1);
                assert_eq!(actual.0[0].numerator, Atom::one());
                assert_eq!(
                    actual.0[0].integral,
                    function!(
                        vakint::symbols::S.topo,
                        function!(
                            vakint::symbols::S.prop,
                            1,
                            function!(vakint::symbols::S.edge, 0, 0),
                            function!(vakint::symbols::S.k, 0),
                            if substitute_masses_to_m_uv {
                                Atom::var(GS.m_uv_vacuum).pow(2)
                            } else {
                                mass_squared.clone()
                            },
                            3
                        )
                    )
                );
            }
        }
    }

    #[test]
    fn vakint_conversion_rejects_a_lost_active_loop() {
        use crate::{finalized_runtime_dot, graph::Graph};

        test_initialise().unwrap();
        let graph: Graph =
            finalized_runtime_dot!(digraph cancelled_tadpole_loop { graph [projector=1]
                edge [num=1 mass=2]
                node [num=1]
                a -> a [id=0 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"]
                a -> a [id=1 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:3}"]
            })
            .unwrap();
        let remaining = graph.get_edge_subgraph(EdgeIndex(1));
        // This is the denominator domain after the other tadpole propagator
        // cancels against its edge-local numerator. The original two-loop
        // measure cannot be replaced by the remaining massive one-loop value.
        let input = graph.denominator(&remaining, &graph.model, |_| -3);
        let error = to_vakint_integrand(
            &input,
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                graph.loop_momentum_basis.clone(),
            )],
            &VakintSettings::default(),
            false,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("has 1 loop generators, expected 2")
        );
        let remaining_lmb = graph.lmb_of(&remaining);
        let converted = to_vakint_integrand(
            &input,
            &graph,
            &graph.full_filter(),
            &graph.empty_subgraph::<SuBitGraph>(),
            &[(
                graph.full_filter(),
                graph.full_filter(),
                remaining_lmb.clone(),
            )],
            &VakintSettings::default(),
            false,
        )
        .unwrap();
        assert_eq!(converted.0.len(), 1);
    }

    #[test]
    fn affine_vakint_routing_keeps_external_momenta_fixed() {
        use crate::{finalized_runtime_dot, graph::Graph};

        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph affine_vakint_tadpole {
            projector=1
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            a -> a [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:2}"]
            a -> outgoing [id=2 source="{ufo_order:3}"]
        })
        .unwrap();
        let hard = function!(GS.emr_mom, 1) - function!(GS.emr_mom, 0);
        let minkowski = Minkowski {}.new_rep(GS.dim).to_symbolic([]);
        let indexed_hard =
            function!(GS.emr_mom, 1, &minkowski) - function!(GS.emr_mom, 0, &minkowski);
        let numerator = function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, 0, &minkowski),
            &indexed_hard
        ) * function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, 2, &minkowski),
            &indexed_hard
        );
        let settings = VakintSettings {
            additional_normalization: "1".to_string(),
            ..Default::default()
        };
        // A retained parent loop overrides an earlier component's external
        // classification. This adapter diagnostic changes coordinate metadata
        // only; both signs must still preserve the fixed physical externals.
        let mut child_frame = graph.empty_lmb();
        child_frame.ext_edges.push(EdgeIndex(1));
        for source_lmbs in [
            vec![&graph.loop_momentum_basis],
            vec![&child_frame, &graph.loop_momentum_basis],
        ] {
            for sign in [-1, 1] {
                let denominator = function!(
                    GS.den,
                    1,
                    Atom::num(sign) * &hard,
                    Atom::var(GS.m_uv_vacuum).pow(2),
                    function!(SPENSO_TAG.dot, &indexed_hard, &indexed_hard)
                        - Atom::var(GS.m_uv_vacuum).pow(2)
                );
                let actual = to_vakint_integrand(
                    &(&numerator / denominator.pow(3)),
                    &graph,
                    &graph.full_filter(),
                    &graph.empty_subgraph::<SuBitGraph>(),
                    &source_lmbs
                        .iter()
                        .map(|lmb| {
                            let owners = if lmb.loop_edges.is_empty() {
                                graph.empty_subgraph::<SuBitGraph>()
                            } else {
                                graph.full_filter()
                            };
                            (owners.clone(), owners, (*lmb).clone())
                        })
                        .collect::<Vec<_>>(),
                    &settings,
                    true,
                )
                .unwrap();
                // Compare the complete rational integrand through Vakint's public
                // expression conversion, independent of topology IDs and term layout.
                let actual = Atom::from(actual)
                    .replace(function!(vakint::symbols::S.topo, W_.x_))
                    .with(W_.x_)
                    .replace(function!(
                        vakint::symbols::S.prop,
                        W_.a_,
                        W_.b_,
                        W_.c_,
                        W_.d_,
                        W_.e_
                    ))
                    .with(
                        (function!(vakint::symbols::S.dot, W_.c_, W_.c_) - Atom::var(W_.d_))
                            .pow(-Atom::var(W_.e_)),
                    );
                let expected_numerator = function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.p, 0),
                    function!(vakint::symbols::S.k, 0)
                ) * function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.p, 2),
                    function!(vakint::symbols::S.k, 0)
                );
                let radial = function!(
                    vakint::symbols::S.dot,
                    function!(vakint::symbols::S.k, 0),
                    function!(vakint::symbols::S.k, 0)
                ) - Atom::var(GS.m_uv_vacuum).pow(2);
                assert!(
                    (actual.collect_factors()
                        - (expected_numerator / radial.pow(3)).collect_factors())
                    .is_zero(),
                    "the complete integrand must retain the fixed external momenta for either D(H) spelling"
                );
            }
        }
    }

    #[test]
    fn nested_active_sunset_keeps_its_two_loop_incidence() {
        use crate::{finalized_runtime_dot, graph::Graph, momentum::sample::LoopIndex};

        test_initialise().unwrap();
        let graph: Graph =
            finalized_runtime_dot!(digraph nested_active_sunset { graph [projector=1]
                edge [num=1 mass=1]
                node [num=1]
                ext [style=invis]
                ext -> a:0 [id=0 sink="{ufo_order:0}"]
                d:1 -> ext [id=1 source="{ufo_order:0}"]
                a -> b [id=2 source="{ufo_order:1}" sink="{ufo_order:0}"]
                b -> c [id=3 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:0}"]
                b -> c [id=4 source="{ufo_order:2}" sink="{ufo_order:1}"]
                c -> d [id=5 source="{ufo_order:2}" sink="{ufo_order:1}"]
                a -> d [id=6 lmb_id=2 source="{ufo_order:2}" sink="{ufo_order:2}"]
                b -> c [id=7 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:3}"]
            })
            .unwrap();
        let full = graph
            .full_filter()
            .subtract(&graph.external_filter::<SuBitGraph>());
        let mut child = graph.empty_subgraph::<SuBitGraph>();
        for edge in [3, 4, 7] {
            child.add(graph[&EdgeIndex(edge)].1);
        }
        let quotient = full.subtract(&child);
        let child_lmb = graph
            .try_compatible_sub_lmb(
                &child,
                graph.dummy_less_full_crown(&child),
                &graph.loop_momentum_basis,
            )
            .unwrap();
        let mut quotient_lmb = graph.loop_momentum_basis.clone();
        for index in (0..quotient_lmb.loop_edges.len()).rev() {
            let index = LoopIndex(index);
            if child_lmb
                .loop_edges
                .contains(&quotient_lmb.loop_edges[index])
            {
                quotient_lmb.put_loop_to_ext(index);
            }
        }
        let [k, l, p] = [3, 7, 6].map(|edge| function!(GS.emr_mom, edge));
        let mass_squared = Atom::var(GS.m_uv_vacuum).pow(2);
        let child_denominators = function!(GS.den, 3, &k, &mass_squared)
            * function!(GS.den, 7, &l, &mass_squared)
            * function!(GS.den, 4, -&k - &l, &mass_squared).pow(3);
        let quotient_denominators = function!(GS.den, 2, -&p, &mass_squared)
            * function!(GS.den, 5, -&p, &mass_squared)
            * function!(GS.den, 6, &p, &mass_squared);
        let numerator = parse!("(a+b)*(c+d)");
        for child_active in [true, false] {
            let mut components = Vec::new();
            let input = if child_active {
                components.push((child.clone(), child.clone(), child_lmb.clone()));
                &numerator / (&child_denominators * &quotient_denominators)
            } else {
                &numerator / &quotient_denominators
            };
            components.push((quotient.clone(), full.clone(), quotient_lmb.clone()));
            let converted = to_vakint_integrand(
                &input,
                &graph,
                &full,
                &child,
                &components,
                &VakintSettings::default(),
                false,
            )
            .unwrap();
            assert_eq!(converted.0.len(), 1);
            let term = &converted.0[0];
            assert_eq!(term.numerator, numerator);
            let propagator = function!(
                vakint::symbols::S.prop,
                W_.a_,
                W_.b_,
                W_.mom_,
                W_.mass_,
                W_.e_
            )
            .to_pattern();
            let mut powers = term
                .integral
                .pattern_match(&propagator, None, None)
                .map(|matched| {
                    assert_eq!(matched[&W_.mass_], mass_squared);
                    assert!(!matched[&W_.mom_].contains_symbol(vakint::symbols::S.p));
                    i32::try_from(matched[&W_.e_].as_view()).unwrap()
                })
                .collect::<Vec<_>>();
            powers.sort();
            assert_eq!(
                powers,
                if child_active {
                    vec![1, 1, 3, 3]
                } else {
                    vec![3]
                }
            );
        }
    }

    #[test]
    fn nested_vacuum_denominators_preserve_incidence_and_factorized_numerators() {
        use crate::{finalized_runtime_dot, graph::Graph};

        test_initialise().unwrap();
        let graph: Graph = finalized_runtime_dot!(digraph nested_vacuum_bubble {
            projector=1
            edge [num=1 mass=1]
            node [num=1]
            incoming [style=invis]
            outgoing [style=invis]
            incoming -> a [id=0 sink="{ufo_order:0}"]
            b -> outgoing [id=1 source="{ufo_order:0}"]
            a -> b [id=2 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"]
            a -> b [id=3 source="{ufo_order:2}" sink="{ufo_order:2}"]
        })
        .unwrap();
        let denominators = [(2, 1), (3, -1)].map(|(edge, sign)| {
            function!(
                GS.den,
                edge,
                Atom::num(sign) * function!(GS.emr_mom, 2),
                Atom::var(GS.m_uv_vacuum).pow(2)
            )
        });
        let numerator = parse!("(a+b)*(c+d)");
        let [first, second] = [parse!("s"), parse!("t")];
        // Include powers beyond signed eight-bit range, matching the i32
        // propagator powers used by the incidence conversion below.
        for base_power in [1, 127] {
            let input = &numerator / (&denominators[0] * denominators[1].pow(base_power))
                * (&first / &denominators[1] + &second / denominators[1].pow(2));
            let actual = to_vakint_integrand(
                &input,
                &graph,
                &graph.full_filter(),
                &graph.empty_subgraph::<SuBitGraph>(),
                &[(
                    graph.full_filter(),
                    graph.full_filter(),
                    graph.loop_momentum_basis.clone(),
                )],
                &VakintSettings {
                    // This conversion oracle uses unit loop normalization in every prefix.
                    additional_normalization: "1".to_string(),
                    ..VakintSettings::default()
                },
                true,
            )
            .unwrap();
            let actual = Atom::from(actual)
                .replace(function!(vakint::symbols::S.topo, W_.x_))
                .with(W_.x_)
                .replace(function!(
                    vakint::symbols::S.prop,
                    W_.a_,
                    W_.b_,
                    W_.c_,
                    W_.d_,
                    W_.e_
                ))
                .with(
                    (function!(vakint::symbols::S.dot, W_.c_, W_.c_) - Atom::var(W_.d_))
                        .pow(-Atom::var(W_.e_)),
                );
            let radial = function!(
                vakint::symbols::S.dot,
                function!(vakint::symbols::S.k, 0),
                function!(vakint::symbols::S.k, 0)
            ) - Atom::var(GS.m_uv_vacuum).pow(2);
            let expected = &numerator
                * (&first / radial.pow(base_power + 2) + &second / radial.pow(base_power + 3));
            assert_eq!(actual.collect_factors(), expected.collect_factors());
        }
    }

    #[test]
    fn analytic_fixed_basis_projection_preserves_captured_word_and_factored_term() -> Result<()> {
        test_initialise()?;
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        for fixture in [
            include_str!("../../../tests/resources/uv_fixed_basis/d02_guard_fixed_word.sym"),
            include_str!("../../../tests/resources/uv_fixed_basis/d02_guard_factored_term.sym"),
        ] {
            let source = Atom::parse(fixture, "gammalooprs", Default::default())
                .map_err(|error| eyre!(error))?;
            let native = SymbolicTensor::infer_indexed(cooking.try_cook(source.as_view())?)?;
            let original = native.clone();
            Integrated::ensure_resolved_lorentz_contractions(native.expression())?;
            assert_eq!(native, original);
            assert_eq!(cooking.uncook(native.expression().as_view()), source);
        }
        Ok(())
    }

    #[test]
    fn analytic_fixed_basis_bounds_preserve_formal_dimension_and_laurent_depth() -> Result<()> {
        test_initialise()?;
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::Hedge(881, 1).to_atom()]);
        let b = bis.to_symbolic([Aind::Hedge(882, 1).to_atom()]);
        let power = (Atom::var(symbol!("analytic_fixed_basis::x"))
            + Atom::var(symbol!("analytic_fixed_basis::y")))
        .pow(9);
        let inverse = (Atom::var(symbol!("analytic_fixed_basis::u"))
            + Atom::var(symbol!("analytic_fixed_basis::v")))
        .pow(-3);
        let spectator = function!(symbol!("analytic_fixed_basis::opaque"), &power) * &inverse;
        let epsilon = Atom::var(GS.dim_epsilon);
        for dimension in [Atom::var(GS.dim), Atom::num(4)] {
            let rep = Minkowski {}.to_symbolic([&dimension]);
            for i in 0..4 {
                for j in 0..4 {
                    let fixed = |component| {
                        function!(ETS.delta, function!(AIND_SYMBOLS.cind, component), &rep)
                    };
                    let word = spenso::chain!(
                        &a,
                        &b,
                        function!(
                            AGS.gamma,
                            SPENSO_TAG.chain_in,
                            SPENSO_TAG.chain_out,
                            fixed(i)
                        ),
                        function!(
                            AGS.gamma,
                            SPENSO_TAG.chain_in,
                            SPENSO_TAG.chain_out,
                            fixed(j)
                        )
                    );
                    for pole in [1, 2] {
                        let coefficient = (Atom::var(GS.dim) - Atom::num(4)) * Atom::var(GS.dim)
                            / epsilon.pow(pole);
                        let source = &spectator * &coefficient * &word;
                        let native =
                            SymbolicTensor::infer_indexed(cooking.try_cook(source.as_view())?)?;
                        let original = native.clone();
                        Integrated::ensure_resolved_lorentz_contractions(native.expression())?;
                        assert_eq!(native, original);
                        assert!(
                            cooking
                                .uncook(native.expression().as_view())
                                .contains(power.as_view())
                        );
                        assert!(
                            cooking
                                .uncook(native.expression().as_view())
                                .contains(inverse.as_view())
                        );
                        // Scalar D remains independent through admission. Only
                        // this reference explicitly performs dimensional substitution.
                        let scalar = coefficient
                            .replace(GS.dim)
                            .with(Atom::num(4) - Atom::num(2) * &epsilon)
                            .series(GS.dim_epsilon, 0, 0)?
                            .to_atom();
                        assert_eq!(
                            scalar,
                            if pole == 1 {
                                Atom::num(-8)
                            } else {
                                Atom::num(-8) / &epsilon + Atom::num(4)
                            }
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn analytic_fixed_basis_projection_keeps_unknown_selectors_and_real_sums_visible() -> Result<()>
    {
        test_initialise()?;
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let [a, b, c, d, e, f] =
            [891, 892, 893, 894, 895, 896].map(|i| bis.to_symbolic([Aind::Hedge(i, 1).to_atom()]));
        for dimension in [Atom::var(GS.dim), Atom::num(4)] {
            let rep = Minkowski {}.to_symbolic([&dimension]);
            for coordinate in [
                function!(AIND_SYMBOLS.cind, -1),
                function!(AIND_SYMBOLS.cind, 4),
                function!(AIND_SYMBOLS.cind, symbol!("analytic_fixed_basis::unknown")),
                function!(AIND_SYMBOLS.cind, 1, 2),
            ] {
                let source = function!(AGS.gamma, &a, &b, function!(ETS.delta, &coordinate, &rep));
                let source = cooking.try_cook(source.as_view())?;
                assert!(Integrated::ensure_resolved_lorentz_contractions(&source).is_err());
            }
            let mu = Minkowski {}
                .to_symbolic([&dimension, &Atom::var(symbol!("analytic_fixed_basis::mu"))]);
            let left = function!(AGS.gamma, &a, &b, &mu);
            let right = function!(AGS.gamma, &c, &d, &mu);
            let malformed = &left
                * function!(ETS.delta, function!(AIND_SYMBOLS.cind, 1), &mu)
                * function!(ETS.delta, function!(AIND_SYMBOLS.cind, 2), &mu);
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(malformed.as_view())?
                )
                .is_err()
            );

            let fixed = function!(
                AGS.gamma,
                &a,
                &b,
                function!(ETS.delta, function!(AIND_SYMBOLS.cind, 1), &rep)
            );
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook((&left * &right).as_view())?
                )
                .is_err()
            );
            // A valid fixed selector in another additive branch cannot hide
            // an actual free D-dimensional contraction.
            let separate_fixed = function!(
                AGS.gamma,
                &e,
                &f,
                function!(ETS.delta, function!(AIND_SYMBOLS.cind, 2), &rep)
            );
            let independent_sum = &separate_fixed * &left * &right;
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(independent_sum.as_view())?
                )
                .is_err()
            );
            let nu = Minkowski {}
                .to_symbolic([&dimension, &Atom::var(symbol!("analytic_fixed_basis::nu"))]);
            let metric_pair = &separate_fixed
                * &left
                * function!(ETS.metric, &mu, &nu)
                * function!(AGS.gamma, &c, &d, &nu);
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(metric_pair.as_view())?
                )
                .is_err()
            );
            let ordinary_vector =
                SPENSO_TAG.rank_one_tensor_symbol("analytic_fixed_basis::unclassified_vector");
            let unclassified = function!(
                AGS.gamma,
                &a,
                &b,
                function!(ordinary_vector, function!(AIND_SYMBOLS.cind, 1), &rep)
            );
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(unclassified.as_view())?
                )
                .is_err()
            );
            let source = &left * &right + &fixed * function!(ETS.metric, &c, &d);
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(source.as_view())?
                )
                .is_err()
            );
            let compact = function!(
                ETS.metric,
                function!(AGS.gamma, &a, &b, &rep),
                function!(AGS.gamma, &c, &d, &rep)
            );
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(compact.as_view())?
                )
                .is_err()
            );
            let dot = function!(
                SPENSO_TAG.dot,
                function!(AGS.gamma, &a, &b, &rep),
                function!(AGS.gamma, &c, &d, &rep)
            );
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(&cooking.try_cook(dot.as_view())?)
                    .is_err()
            );
        }
        // Generic fixed-basis bounds and the physical analytic subspace are
        // distinct: (2,2) is out of bounds; (6,4) is bounded but nonphysical.
        for (dimension, component) in [(2, 2), (6, 4)] {
            let rep = Minkowski {}.to_symbolic([Atom::num(dimension)]);
            let source = function!(
                AGS.gamma,
                &a,
                &b,
                function!(ETS.delta, function!(AIND_SYMBOLS.cind, component), &rep)
            );
            assert!(
                Integrated::ensure_resolved_lorentz_contractions(
                    &cooking.try_cook(source.as_view())?
                )
                .is_err()
            );
        }
        Ok(())
    }

    #[test]
    fn analytic_fixed_gamma_one_two_has_independent_nonzero_weyl_components() -> Result<()> {
        use crate::{
            numerator::symbolica_ext::NumeratorAtomExt,
            utils::{FUN_LIB, TENSORLIB},
        };
        use idenso::shorthands::UndoShorthands;
        use spenso::network::{ExecutionResult, Sequential, SmallestDegree};
        test_initialise()?;
        let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
        let bis = Bispinor {}.new_rep(4);
        let a = bis.to_symbolic([Aind::Hedge(901, 1).to_atom()]);
        let b = bis.to_symbolic([Aind::Hedge(902, 1).to_atom()]);
        let rep = Minkowski {}.new_rep(4).to_symbolic([]);
        let fixed = |i| function!(ETS.delta, function!(AIND_SYMBOLS.cind, i), &rep);
        let word = spenso::chain!(
            &a,
            &b,
            function!(
                AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                fixed(1)
            ),
            function!(
                AGS.gamma,
                SPENSO_TAG.chain_in,
                SPENSO_TAG.chain_out,
                fixed(2)
            )
        );
        // Weyl gamma1*gamma2 = -i diag(sigma3,sigma3). Both spatial
        // covector signs cancel; this is a nonzero matrix, not a vanishing word.
        for row in 0..4 {
            for column in 0..4 {
                let component = &word
                    * function!(ETS.delta, function!(AIND_SYMBOLS.cind, row), &a)
                    * function!(ETS.delta, function!(AIND_SYMBOLS.cind, column), &b);
                let explicit = component.undo_dots::<Aind>()?;
                let mut network = explicit.parse_into_net()?;
                network.execute::<Sequential, SmallestDegree, _, _, _>(
                    &*TENSORLIB.read().unwrap(),
                    &*FUN_LIB,
                )?;
                let actual = match network.result_scalar()? {
                    ExecutionResult::One => Atom::one(),
                    ExecutionResult::Zero => Atom::Zero,
                    ExecutionResult::Val(value) => value.into_owned(),
                };
                let expected = if row != column {
                    Atom::Zero
                } else if row % 2 == 0 {
                    -Atom::i()
                } else {
                    Atom::i()
                };
                assert_eq!(actual, expected, "row={row}, column={column}");
            }
        }
        // Run the independent nonzero reference even on the expected-RED base.
        Integrated::ensure_resolved_lorentz_contractions(&cooking.try_cook(word.as_view())?)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "integrated_native_mass_tests.rs"]
pub(in crate::uv) mod native_mass_tests;
