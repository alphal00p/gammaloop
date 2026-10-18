//! Actual quark-loop gluon prefixes in the native two-loop top self-energy.

use super::*;
use crate::{
    graph::{Graph, LMBext, lmb::LMBwithEdges},
    numerator::AppliedFeynmanRule,
    utils::{GS, W_},
    uv::{
        ApproximationType, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike, Rooted, UVCtx,
            local_4d::{FourDTerm, Full4dCts, Local4dCts, uv_limit},
        },
        hedge_poset::OwnedForestNode,
        settings::RenormalizationPrescriptionSettings,
    },
};
use linnet::half_edge::subgraph::{SubSetLike, SubSetOps};
use spenso::structure::slot::IsAbstractSlot;

/// Actual IR child and quotient data; no transverse expectation is stored here.
struct CertifiedGluonPrefix {
    pub(super) diagram: FeynmanDiagram,
    pub(super) graph: Graph,
    pub(super) child: OwnedForestNode,
    pub(super) species: isize,
    pub(super) mass: Atom,
    pub(super) weight: Atom,
    pub(super) spin: Vec<Atom>,
    pub(super) color: Vec<Atom>,
    pub(super) crown_edges: Vec<linnet::half_edge::involution::EdgeIndex>,
    pub(super) q: Atom,
    pub(super) local: Local4dCts,
    pub(super) integrated: super::super::IntegratedCts,
    pub(super) denominators: Vec<crate::graph::FourDDenominator>,
    pub(super) actual_child: SymbolicTensor<IndexedStructure>,
    pub(super) actual_numerator: Atom,
}

impl QcdMassReference {
    fn certified_gluon_prefixes(
        &self,
        internal_species: &BTreeSet<isize>,
    ) -> Result<Vec<CertifiedGluonPrefix>> {
        let d = Atom::var(GS.dim);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let g = Atom::var(symbol!("UFO::G"));
        let measure = Atom::i() / (Atom::num(16) * pi.pow(2));
        let a = g.pow(2) / (Atom::num(16) * pi.pow(2));
        let settings = UVgenerationSettings {
            generate_integrated: true,
            add_marker: false,
            renormalization_prescription: RenormalizationPrescriptionSettings {
                log_divergent: ApproximationType::MUV,
                massive_power_divergent: ApproximationType::OS,
                massless_power_divergent: ApproximationType::IR,
                ..Default::default()
            },
            ..Default::default()
        };
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let quadratic = |momentum: &Atom, mass_squared: &Atom| {
            let vector =
                GS.indexed_momentum(momentum, &[Minkowski {}.new_rep(GS.dim).to_symbolic([])]);
            function!(spenso::network::tags::SPENSO_TAG.dot, &vector, &vector) - mass_squared
        };
        let mut records = Vec::new();
        for diagram in self.sources(2)? {
            let graph = Graph::from_feynkit(&diagram, None, true)?;
            let model = &graph.model;
            let empty = graph.empty_subgraph();
            let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
            let root = OwnedForestNode {
                spinney: Spinney::empty(&graph),
                topo_order: 0,
            };
            let classified =
                graph.classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?;
            for spinney in &classified {
                if spinney.n_components() != 1
                    || spinney.dod != 2
                    || graph.n_loops(spinney.filter()) != 1
                    || spinney.renormalization_scheme != ApproximationType::IR
                {
                    continue;
                }
                let identifier = graph.ct_identifier(spinney.filter(), model);
                let Some(flavours) = &identifier.internal_pdg_set else {
                    continue;
                };
                if flavours.len() != 1 || !flavours.iter().all(|pdg| internal_species.contains(pdg))
                {
                    continue;
                }
                ensure!(
                    identifier
                        .external_pdg_set
                        .iter()
                        .all(|pdg| pdg.abs() == 21)
                );
                let flavour = *flavours.first().unwrap();
                let child = OwnedForestNode {
                    spinney: spinney.clone(),
                    topo_order: 1,
                };
                let edges = graph.paired_edges(child.subgraph());
                ensure!(
                    (edges.len() == 2 || (edges.len() == 1 && flavour == 21))
                        && child.lmb().loop_edges.len() == 1
                );
                let mass = graph[edges[0]].particle.mass_atom(model);
                ensure!(
                    edges
                        .iter()
                        .all(|edge| graph[*edge].particle.mass_atom(model) == mass)
                );
                let crown = graph.dummy_less_full_crown(child.subgraph());
                let boundary = crown.included_iter().collect::<Vec<_>>();
                ensure!(boundary.len() == 2);
                let crown_edges = boundary
                    .iter()
                    .map(|hedge| graph[hedge])
                    .collect::<Vec<_>>();
                ensure!(
                    crown_edges
                        .iter()
                        .all(|edge| graph[*edge].particle.mass_atom(model).is_zero())
                );
                let spin = boundary
                    .iter()
                    .map(|hedge| {
                        let slots = graph.underlying[*hedge]
                            .edge_spin_slots()
                            .collect::<Vec<_>>();
                        ensure!(slots.len() == 1);
                        Ok(slots[0].to_atom().map_mink_dim(d.clone()))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let color = boundary
                    .iter()
                    .map(|hedge| {
                        let slots = graph.underlying[*hedge]
                            .edge_color_slots()
                            .collect::<Vec<_>>();
                        ensure!(slots.len() == 1);
                        Ok(slots[0].to_atom())
                    })
                    .collect::<Result<Vec<_>>>()?;
                let q = function!(GS.emr_mom, crown_edges[0].0);
                let ctx = UVCtx::new(&graph, model, &settings);
                let local = uv_limit(
                    &Full4dCts::from_factorized_local(&Local4dCts::root()),
                    &ctx,
                    &child,
                    &root,
                    &child,
                    &root,
                )?;
                let [sector] = local.active_sectors() else {
                    return Err(eyre!("one actual IR child sector required"));
                };
                ensure!(sector.active_components.len() == 1 && sector.frozen_lmbs().is_empty());
                ensure!(&sector.active_components[0].0 == child.subgraph());
                ensure!(&sector.active_components[0].1 == child.subgraph());
                ensure!(sector.active_components[0].2.loop_edges == child.lmb().loop_edges);
                ensure!(
                    sector.integration_atom().same_value(&sector.atom),
                    "bare IR child integrates its complete actual U+S-US expression"
                );
                let _source_certificate = sector.canonical_projection(&graph)?;
                let backend = vakint::VakintSettings {
                    number_of_terms_in_epsilon_expansion: 3,
                    evaluation_order: vakint::EvaluationOrder(vec![
                        vakint::EvaluationMethod::AlphaLoop(vakint::AlphaLoopOptions {
                            susbstitute_masters: false,
                        }),
                    ]),
                    ..settings.vakint.true_settings()
                };
                let integrated = super::super::Integrated::new(crate::utils::vakint()?, &backend)
                    .run(&local, &ctx, &child, &root, &child, &root)?;
                let completed = Full4dCts::recursion_input(
                    &local,
                    &integrated,
                    ApproximationType::IR,
                    false,
                    child.lmb(),
                )?;
                let frozen = completed
                    .sectors()
                    .filter(|sector| !sector.frozen_lmbs().is_empty())
                    .collect::<Vec<_>>();
                let physical_finite = integrated.physical_finite_counterterm()?;
                if edges.len() == 1 {
                    // The actual massless seagull is integrated before its zero
                    // completion is pruned. Retain its nonzero local source.
                    ensure!(flavour == 21 && mass.is_zero());
                    ensure!(integrated.physical_pole()?.is_zero());
                    ensure!(physical_finite.is_zero());
                    ensure!(frozen.is_empty());
                } else {
                    ensure!(frozen.len() == 1 && frozen[0].active_components.is_empty());
                    ensure!(frozen[0].frozen_lmbs() == std::slice::from_ref(child.lmb()));
                }
                let finite = frozen.iter().try_fold(
                    physical_finite.with_checked_expression(Atom::Zero)?,
                    |sum, sector| sum.try_add(&sector.atom),
                )?;
                let reduced = super::super::simplify(&finite.try_sub(
                    &physical_finite.with_checked_expression(
                        physical_finite.expression() * Atom::var(GS.integrated_loop_scale).pow(4),
                    )?,
                )?)?;
                Self::certify_finite_tagged_rows_zero(
                    "the frozen finite insertion owns one consumed loop measure",
                    reduced.with_checked_expression(
                        reduced.expression().to_cof_dimension_invariants(),
                    )?,
                )?;

                let projection = self.project(&diagram)?;
                ensure!(
                    diagram.numerator_prefactor().is_scalar()
                        && diagram.numerator_prefactor().expression().is_one()
                );
                let actual_child = integrated.physical_pole()?;
                let actual_child = actual_child.with_checked_expression(
                    model
                        .expand_couplings(&(actual_child.expression() * &projection.kernel_weight))
                        .to_cof_dimension_invariants(),
                )?;
                ensure!(
                    !actual_child.expression().contains_symbol(GS.m_uv_vacuum)
                        && !actual_child.expression().contains_symbol(GS.m_uv_expansion)
                        && !actual_child.expression().contains_symbol(GS.mu_r_sq)
                );
                let full_integral = integrated.physical_pole()?.try_sub(&physical_finite)?;
                let finite = finite.with_checked_expression(
                    finite
                        .expression()
                        .replace(GS.integrated_loop_scale)
                        .with(1),
                )?;
                let difference = full_integral
                    .try_add(&finite)?
                    .try_sub(&integrated.physical_pole()?)?;
                let reduced = super::super::simplify(&difference.with_checked_expression(
                    super::super::series(difference.expression(), 1)?.to_atom(),
                )?)?;
                Self::certify_finite_tagged_rows_zero(
                    "local integral plus actual completion through epsilon one",
                    reduced.with_checked_expression(
                        reduced.expression().to_cof_dimension_invariants(),
                    )?,
                )?;

                let parent = classified
                    .iter()
                    .find(|spinney| {
                        spinney.n_components() == 1 && graph.n_loops(spinney.filter()) == 2
                    })
                    .ok_or_else(|| eyre!("native two-loop parent required"))?;
                let parent = OwnedForestNode {
                    spinney: parent.clone(),
                    topo_order: 2,
                };
                let reduced = parent.reduced_subgraph(&child);
                let outside = graph.paired_edges(&reduced);
                ensure!(outside.len() == 3);
                let quarks = outside
                    .iter()
                    .copied()
                    .filter(|edge| !graph[*edge].particle.mass_atom(model).is_zero())
                    .collect::<Vec<_>>();
                ensure!(quarks.len() == 1);
                let quark = quarks[0];
                ensure!(
                    outside
                        .iter()
                        .filter(|edge| **edge != quark)
                        .copied()
                        .collect::<BTreeSet<_>>()
                        == crown_edges.iter().copied().collect()
                );
                for edge in &outside {
                    ensure!(
                        child.lmb().edge_signatures[*edge]
                            .internal
                            .iter()
                            .all(|sign| sign.is_zero()),
                        "the cograph must be independent of the consumed hard coordinate"
                    );
                }
                let child_raw = graph
                    .numerator(child.subgraph(), &empty)?
                    .to_d_dim(GS.dim)?
                    .state
                    .expr;
                let cograph = graph
                    .numerator(&reduced, child.subgraph())?
                    .to_d_dim(GS.dim)?
                    .state
                    .expr;
                let full = graph
                    .numerator(parent.subgraph(), &empty)?
                    .to_d_dim(GS.dim)?
                    .state
                    .expr;
                ensure!(
                    child_raw
                        .multiply(cograph.clone())?
                        .try_sub(&full)?
                        .is_zero(),
                    "native source/cograph factor ownership"
                );
                let quotient = graph.lmb_with_loop_edges(
                    &[child.lmb().loop_edges.first().copied().unwrap(), quark][..],
                )?;
                let chart =
                    graph.uv_wrapped_replacement(&graph.full_filter(), &quotient, &[W_.x___]);
                let k = function!(GS.emr_mom, quark.0);
                let p = function!(GS.emr_mom, 0);
                let scalar_chart =
                    graph.uv_wrapped_replacement(&graph.full_filter(), &quotient, &[] as &[Atom]);
                ensure!(k.replace_multiple(&scalar_chart) == k);
                for edge in &crown_edges {
                    let routed = function!(GS.emr_mom, edge.0).replace_multiple(&scalar_chart);
                    ensure!(
                        routed == &k - &p || routed == &p - &k,
                        "signed cograph routing must retain the external shift"
                    );
                }
                let mass_top = graph[quark].particle.mass_atom(model);
                let denominator = graph.denominator(&reduced, model, |_| -1);
                let terms = FourDTerm::from_view(denominator.as_view())?;
                ensure!(terms.len() == 1 && terms[0].denominators.len() == 3);
                for den in &terms[0].denominators {
                    let mass_squared = if den.source_edge == quark {
                        mass_top.pow(2)
                    } else {
                        Atom::Zero
                    };
                    ensure!(den.mass_squared == mass_squared);
                    let reduced = super::super::simplify(
                        &AppliedFeynmanRule::new(
                            GS.erase_uv_momentum_provenance(&den.full_expr)
                                - quadratic(
                                    &function!(GS.emr_mom, den.source_edge.0),
                                    &mass_squared,
                                ),
                        )?
                        .expr,
                    )?;
                    Self::certify_finite_tagged_rows_zero(
                        "cograph full polynomial retains its power and prescription",
                        reduced.with_checked_expression(
                            reduced.expression().to_cof_dimension_invariants(),
                        )?,
                    )?;
                }
                ensure!(
                    terms[0]
                        .denominators
                        .iter()
                        .map(|den| den.source_edge)
                        .collect::<BTreeSet<_>>()
                        == outside.iter().copied().collect()
                );

                let denominators = terms
                    .into_iter()
                    .next()
                    .expect("the single complete cograph term was certified")
                    .denominators;

                let (projector, projected_external) = self.mass_projector(&diagram, &d)?;
                ensure!(
                    projected_external == p,
                    "the mass projector owns the same external momentum as the quotient chart"
                );
                let projector = projector.map_representation_dimensions(|_, dimension| {
                    dimension.replace(CS.nc).with(3)
                })?;
                let projector = projector
                    .with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
                let projector = projector.try_into_indexed()?;
                let closed = actual_child
                    .clone()
                    .multiply(cograph)?
                    .multiply(projector)?;
                let closed_expression = model
                    .expand_couplings(closed.expression())
                    .replace(mass_top.to_pattern())
                    .with(1);
                let tensor = closed.with_checked_expression(
                    GS.erase_uv_momentum_provenance(&closed_expression)
                        .replace_multiple(&chart),
                )?;
                ensure!(tensor.structure().canonical().is_scalar());
                let tensor = tensor
                    .simplify_algebra(&AlgebraSettings {
                        gamma: Some(GammaSimplifySettings::canonical()),
                        color: Some(
                            ColorSimplifySettings::default().with_cof_dimension_invariants(),
                        ),
                        ..Default::default()
                    })?
                    .contract(ContractSettings::default())?;
                ensure!(tensor.contraction_complete());
                let traced = cooking
                    .uncook(tensor.to_dots()?.expression().as_view())
                    .to_cof_dimension_invariants();
                let pk = Atom::var(symbol!("native_gluon_prefix::pk"));
                let kk = Atom::var(symbol!("native_gluon_prefix::kk"));
                let kin = Kinematics::in_dimension(&d)?
                    .with_mass_squared(&p, Atom::one())?
                    .with_scalar_product(&p, &k, pk.clone())?
                    .with_scalar_product(&k, &k, kk.clone())?;
                let scalar = kin.apply(&traced);
                ensure!(!scalar.contains_symbol(GS.emr_mom));
                let actual_numerator = super::super::Integrated::dimensionally_regularized(
                    &(Atom::i() * &measure * scalar / a.pow(2)),
                );
                records.push(CertifiedGluonPrefix {
                    diagram: diagram.clone(),
                    graph: graph.clone(),
                    child,
                    species: flavour,
                    mass,
                    weight: projection.kernel_weight,
                    spin,
                    color,
                    crown_edges,
                    q,
                    local,
                    integrated,
                    denominators,
                    actual_child,
                    actual_numerator,
                });
            }
        }
        Ok(records)
    }
}

#[test]
fn native_fermion_gluon_prefixes_have_the_ms_field_pole_and_cograph_density() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certify_fermion_gluon_prefixes()?;
    Ok(())
}

impl QcdMassReference {
    pub(super) fn certify_fermion_gluon_prefixes(&self) -> Result<CertifiedMassPrefixes> {
        let mut one_loop_sources = self.sources(1)?;
        self.project(&one_loop_sources.remove(0))?
            .certify_one_loop()?;
        longitudinal_gluon_mass_trace_has_only_an_odd_vacuum_remainder()?;
        let d = Atom::var(GS.dim);
        let epsilon = Atom::var(GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let g = Atom::var(symbol!("UFO::G"));
        let a = g.pow(2) / (Atom::num(16) * pi.pow(2));
        let mut seen = BTreeSet::new();
        let mut light_density = None;
        let mut transferred = CertifiedMassPrefixes::new();
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        for record in self.certified_gluon_prefixes(&(1..=6).collect())? {
            let graph = &record.graph;
            let child = &record.child;
            let flavour = record.species;
            ensure!(
                seen.insert(flavour),
                "each actual native flavor occurs once"
            );
            ensure!(graph.paired_edges(child.subgraph()).len() == 2);
            ensure!(record.mass.is_zero() == (flavour != 6));
            ensure!(
                record.weight == Atom::num(-1),
                "retain the closed fermion-loop sign"
            );
            let vector =
                GS.indexed_momentum(&record.q, &[Minkowski {}.new_rep(GS.dim).to_symbolic([])]);
            let transverse = (function!(spenso::network::tags::SPENSO_TAG.dot, &vector, &vector)
                * spenso::g!(&record.spin[0], &record.spin[1])
                - GS.indexed_momentum(&record.q, &[record.spin[0].clone()])
                    * GS.indexed_momentum(&record.q, &[record.spin[1].clone()]))
                * spenso::g!(&record.color[0], &record.color[1]);
            let expected_child = AppliedFeynmanRule::new(
                Atom::i() * &a * Atom::num((2, 3)) / &epsilon * transverse,
            )?
            .expr;
            let routing =
                graph.uv_wrapped_replacement(&graph.full_filter(), child.lmb(), &[W_.x___]);
            let reduced = super::super::simplify(
                &record
                    .actual_child
                    .with_checked_expression(
                        GS.erase_uv_momentum_provenance(record.actual_child.expression())
                            .replace_multiple(&routing),
                    )?
                    .try_sub(&expected_child.with_checked_expression(
                        expected_child.expression().replace_multiple(&routing),
                    )?)?,
            )?;
            Self::certify_finite_tagged_rows_zero(
                "actual open gluon pole with native crown and closed-loop sign",
                reduced
                    .with_checked_expression(reduced.expression().to_cof_dimension_invariants())?,
            )?;
            let pk = Atom::var(symbol!("native_gluon_prefix::pk"));
            let kk = Atom::var(symbol!("native_gluon_prefix::kk"));
            let dk = &kk - Atom::one();
            let dg = &kk - Atom::num(2) * &pk + Atom::one();
            let n = (Atom::num(2) - &d) * &pk + &d;
            let expected_numerator = super::super::Integrated::dimensionally_regularized(
                &(Atom::num((8, 9)) / &epsilon * (&dg * &n - (&pk - Atom::one()) * &dk)),
            );
            ensure!(
                super::super::series(&(&record.actual_numerator - expected_numerator), 1)?
                    .to_atom()
                    .together()
                    .is_zero(),
                "actual gluon proper prefix differs from transverse plus longitudinal Dk Dg^2 density"
            );
            if flavour != 6 {
                if let Some(previous) = &light_density {
                    ensure!(
                        record.actual_numerator == *previous,
                        "five actual light completions agree in the certified chart"
                    );
                } else {
                    light_density = Some(record.actual_numerator.clone());
                }
            }
            let f = Atom::num((4, 3))
                * (Atom::num(3) / &epsilon
                    + Atom::num(4)
                    + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(4)));
            let integral = super::super::series(&(Atom::num((2, 3)) / &epsilon * f), 0)?.to_atom();
            let key = (
                record.diagram.name().to_owned(),
                graph
                    .paired_edges(child.subgraph())
                    .iter()
                    .map(|edge| edge.0)
                    .collect(),
            );
            let unique_source = transferred.insert(key, integral.clone()).is_none();
            ensure!(unique_source);
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_fermion_gluon_prefix_certificate", source = %record.diagram.name(), flavour,
                file.child_edges = ?graph.paired_edges(child.subgraph()), file.crown_edges = ?record.crown_edges,
                file.actual_pole = %cooking.uncook(record.actual_child.expression().as_view()).to_canonical_string(),
                file.actual_numerator = %record.actual_numerator.to_canonical_string(),
                file.transferred_integral = %integral.to_canonical_string(),
                epsilon_input_depth = 1, integrated_longitudinal_identity_checked_here = true,
                parent_overall_operator_computed = false,
                "Actual fermion-loop prefix, independent cograph density and analytic Ward transfer"
            );
        }
        ensure!(seen == (1..=6).collect::<BTreeSet<_>>());
        Ok(transferred)
    }
}

#[test]
fn native_pure_gauge_gluon_prefixes_have_the_grouped_ms_pole_and_mass_density() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certify_pure_gauge_gluon_prefixes()?;
    Ok(())
}

impl QcdMassReference {
    pub(super) fn certify_pure_gauge_gluon_prefixes(&self) -> Result<CertifiedMassPrefixes> {
        use spenso::{
            network::tags::SPENSO_TAG,
            structure::{abstract_index::AbstractIndex, representation::LibraryRep, slot::Slot},
        };
        use std::collections::HashMap;
        use symbolica::atom::FunctionBuilder;

        let mut one_loop_sources = self.sources(1)?;
        self.project(&one_loop_sources.remove(0))?
            .certify_one_loop()?;
        longitudinal_gluon_mass_trace_has_only_an_odd_vacuum_remainder()?;
        let records = self.certified_gluon_prefixes(&BTreeSet::from([21, 9000005]))?;
        ensure!(records.len() == 3);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let d = Atom::var(GS.dim);
        let epsilon = Atom::var(GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let coupling = Atom::var(symbol!("UFO::G"));
        let a = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
        let q_head = symbol!(
            "native_gauge_prefix::q",
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        );
        let q = Atom::var(q_head);
        let q2 = Atom::var(symbol!("native_gauge_prefix::q_squared"));
        let kin = Kinematics::in_dimension(&d)?.with_mass_squared(&q, q2.clone())?;
        let indices = [
            symbol!("native_gauge_prefix::mu"),
            symbol!("native_gauge_prefix::nu"),
            symbol!("native_gauge_prefix::a"),
            symbol!("native_gauge_prefix::b"),
        ];
        let original_ports = records[0]
            .spin
            .iter()
            .chain(&records[0].color)
            .collect::<Vec<_>>();
        let common_ports = original_ports
            .iter()
            .zip(indices)
            .map(|(port, index)| {
                let cooked = cooking.try_cook(port.as_view())?;
                let slot = Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view())?;
                Ok(slot
                    .rep()
                    .slot::<AbstractIndex, _>(AbstractIndex::from(index))
                    .to_atom())
            })
            .collect::<Result<Vec<_>>>()?;
        let common_ports = common_ports
            .iter()
            .map(|port| cooking.uncook(port.as_view()))
            .collect::<Vec<_>>();
        for (i, port) in common_ports.iter().enumerate() {
            let AtomView::Fun(port) = port.as_view() else {
                return Err(eyre!("explicit native crown slot required"));
            };
            ensure!(port.get_nargs() == 2);
            if i < 2 {
                ensure!(port.get_symbol() == symbol!("spenso::mink") && port.get(0) == d.as_view());
            } else {
                ensure!(
                    port.get_symbol() == CS.adjoint_rep && port.get(0) == Atom::num(8).as_view()
                );
            }
        }
        let metric = spenso::g!(&common_ports[0], &common_ports[1]);
        let color = spenso::g!(&common_ports[2], &common_ports[3]);
        let qq = function!(q_head, &common_ports[0]) * function!(q_head, &common_ports[1]);
        // Only the standalone rank probe needs an explicit free axis. Use
        // the actual named crown port; compact q remains inside closed dots.
        let q_probe = SymbolicTensor::infer_indexed(
            cooking.try_cook(function!(q_head, &common_ports[0]).as_view())?,
        )?;
        let q_port = cooking.try_cook(common_ports[0].as_view())?;
        let q_port = Slot::<LibraryRep, AbstractIndex>::try_from(q_port.as_view())?;
        ensure!(q_probe.rank() == 1);
        ensure!(q_probe.structure().logical_slots() == vec![q_port]);
        let q_vector = Minkowski {}.new_rep(GS.dim).vector(q.as_view(), []);
        let q_squared = function!(SPENSO_TAG.dot, &q_vector, &q_vector);
        let mut aligned_poles = Vec::new();
        let mut sum_a = Atom::Zero;
        let mut sum_b = Atom::Zero;
        let mut density_sum = Atom::Zero;
        let mut transferred = BTreeMap::new();
        let mut inventory = BTreeSet::new();
        for record in records {
            let graph = &record.graph;
            let child = &record.child;
            let edges = graph.paired_edges(child.subgraph());
            let (kind, name, weight) = match (record.species, edges.len()) {
                (9000005, 2) => ("ghost", "D00", Atom::num(-1)),
                (21, 2) => ("gluon", "D01", Atom::num((1, 2))),
                (21, 1) => ("tadpole", "D03", Atom::num((1, 2))),
                _ => return Err(eyre!("unexpected native pure-gauge child")),
            };
            ensure!(record.diagram.name() == name && record.weight == weight);
            ensure!(inventory.insert(kind) && record.mass.is_zero());
            ensure!(
                !record.local.atom().is_zero(),
                "D03 and both bubble sources remain local"
            );
            ensure!(record.denominators.len() == 3);
            if kind == "tadpole" {
                ensure!(record.integrated.physical_pole()?.is_zero());
                ensure!(record.integrated.physical_finite_counterterm()?.is_zero());
                ensure!(record.actual_numerator.is_zero());
            }

            // Certify the signed native crown momentum before using one neutral
            // coordinate. This is a chart for comparison, never new EMR ownership.
            let chart =
                graph.uv_wrapped_replacement(&graph.full_filter(), child.lmb(), &[] as &[Atom]);
            let q_native = record.q.replace_multiple(&chart);
            let names = child
                .lmb()
                .loop_edges
                .iter()
                .chain(child.lmb().ext_edges.iter())
                .map(|edge| function!(GS.emr_mom, edge.0))
                .collect::<BTreeSet<_>>();
            let linear = Kinematics::in_dimension(&d)?.with_momenta(names)?;
            let terms = linear.linear_terms(&q_native)?;
            let [(carrier, sign)] = terms.as_slice() else {
                return Err(eyre!("one crown momentum required"));
            };
            ensure!(*sign == Atom::one() || *sign == Atom::num(-1));
            for edge in &record.crown_edges {
                let routed = function!(GS.emr_mom, edge.0).replace_multiple(&chart);
                ensure!(
                    (&routed - &q_native).expand().is_zero()
                        || (&routed + &q_native).expand().is_zero()
                );
            }
            let actual_source = cooking.uncook(record.actual_child.expression().as_view());
            let routed = GS
                .erase_uv_momentum_provenance(&actual_source)
                .replace_multiple(graph.uv_wrapped_replacement(
                    &graph.full_filter(),
                    child.lmb(),
                    &[W_.x___],
                ));
            let neutral = routed.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    assert_eq!(
                        function!(GS.emr_mom, momentum.get(0)),
                        *carrier,
                        "completed child must not retain a consumed hard momentum"
                    );
                    let indices = momentum
                        .iter()
                        .skip(1)
                        .map(|arg| arg.to_owned())
                        .collect::<Vec<_>>();
                    **output = sign * FunctionBuilder::new(q_head).add_args(indices).finish();
                }
            });
            ensure!(!neutral.contains_symbol(GS.emr_mom));
            let reversed_momentum = neutral.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == q_head
                {
                    let indices = momentum
                        .iter()
                        .map(|arg| arg.to_owned())
                        .collect::<Vec<_>>();
                    **output = sign * GS.indexed_momentum(carrier, &indices);
                }
            });
            ensure!(
                reversed_momentum == routed,
                "signed crown-coordinate map reconstructs the actual integrated source"
            );
            ensure!(indices.iter().all(|name| !neutral.contains_symbol(*name)));
            let source_interface = graph
                .numerator(child.subgraph(), &graph.empty_subgraph())?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            ensure!(
                record
                    .actual_child
                    .structure()
                    .canonical()
                    .external_structure()
                    == source_interface
                        .structure()
                        .canonical()
                        .external_structure(),
                "the actual integrated child retains its named source interface, including D03 zero"
            );
            let value = cooking.try_cook(neutral.as_view())?;
            // Coordinate comparison changes only momentum leaves. The actual
            // carrier retains its public interface and proof state, even at zero.
            let tensor = record.actual_child.with_checked_expression(value)?;
            ensure!(
                tensor.structure().logical_slots()
                    == record.actual_child.structure().logical_slots()
            );
            let ports = record.spin.iter().chain(&record.color).collect::<Vec<_>>();
            let logical = tensor.structure().logical_slots();
            let mut map = HashMap::new();
            let mut inverse = HashMap::new();
            for (source, target) in ports.iter().zip(&common_ports) {
                let source = cooking.try_cook(source.as_view())?;
                let target = cooking.try_cook(target.as_view())?;
                let target = Slot::<LibraryRep, AbstractIndex>::try_from(target.as_view())?;
                let position = logical
                    .iter()
                    .position(|slot| slot.to_atom() == source)
                    .ok_or_else(|| eyre!("native crown port missing from the integrated child"))?;
                ensure!(logical[position].rep() == target.rep());
                let old = logical[position].aind();
                ensure!(map.insert(position, target.aind()).is_none());
                inverse.insert(position, old);
            }
            ensure!(map.len() == 4 && tensor.rank() == 4);
            let aligned = tensor.reindex_interface_ports(&map)?;
            ensure!(
                aligned
                    .reindex_interface_ports(&inverse)?
                    .same_value(&tensor),
                "public crown reindexing must reconstruct the exact source"
            );
            let aligned = aligned
                .with_checked_expression(aligned.expression().to_cof_dimension_invariants())?;
            let aligned_source = cooking.uncook(aligned.expression().as_view());

            // Two scalar contractions determine A and B. No tensor division or
            // denominator clearing of a graph numerator enters this projection.
            let scalar_trace = |projector: &Atom| -> Result<Atom> {
                let projector =
                    SymbolicTensor::infer_indexed(cooking.try_cook(projector.as_view())?)?;
                let closed = aligned.clone().multiply(projector)?;
                ensure!(closed.structure().canonical().is_scalar());
                let contracted = closed
                    .simplify_algebra(&AlgebraSettings {
                        color: Some(
                            ColorSimplifySettings::default().with_cof_dimension_invariants(),
                        ),
                        ..Default::default()
                    })?
                    .contract(ContractSettings::default())?;
                ensure!(contracted.contraction_complete());
                let scalar = cooking
                    .uncook(contracted.to_dots()?.expression().as_view())
                    .to_cof_dimension_invariants();
                let scalar = kin.apply(&scalar);
                ensure!(!scalar.contains_symbol(q_head));
                Ok(scalar.collect_compact_factors().together())
            };
            let trace = scalar_trace(&(&metric * &color))? / (Atom::num(8) * &q2);
            let longitudinal = scalar_trace(&(&qq * &color))? / (Atom::num(8) * q2.pow(2));
            let coefficient_a = ((&trace - &longitudinal) / (&d - Atom::one())).together();
            let coefficient_b = ((&d * &longitudinal - &trace) / (&d - Atom::one())).together();
            ensure!(!coefficient_a.contains(q2.as_view()) && !coefficient_b.contains(q2.as_view()));
            let reconstructed =
                (&coefficient_a * &q_squared * &metric + &coefficient_b * &qq) * &color;
            // The independently reconstructed zero still has the certified
            // four-port boundary; do not infer a scalar interface for D03.
            let expected =
                aligned.with_checked_expression(cooking.try_cook(reconstructed.as_view())?)?;
            let reduced = super::super::simplify(&aligned.try_sub(&expected)?)?;
            Self::certify_finite_tagged_rows_zero(
                "actual individual gauge pole reconstructs from A and B",
                reduced
                    .with_checked_expression(reduced.expression().to_cof_dimension_invariants())?,
            )?;
            let pk = Atom::var(symbol!("native_gluon_prefix::pk"));
            let kk = Atom::var(symbol!("native_gluon_prefix::kk"));
            let dk = &kk - Atom::one();
            let dg = &kk - Atom::num(2) * &pk + Atom::one();
            let n = (Atom::num(2) - &d) * &pk + &d;
            let field_a = &coefficient_a / (Atom::i() * &a);
            let field_b = &coefficient_b / (Atom::i() * &a);
            // Independent bare-rule/Feynman-parameter pole residues, including
            // native loop signs and AutG once (see bare-rule-derivation.typ).
            let (expected_a, expected_b) = match kind {
                "ghost" => (Atom::num((-1, 4)) / &epsilon, Atom::num((-1, 2)) / &epsilon),
                "gluon" => (
                    Atom::num((-19, 4)) / &epsilon,
                    Atom::num((11, 2)) / &epsilon,
                ),
                "tadpole" => (Atom::Zero, Atom::Zero),
                _ => unreachable!(),
            };
            ensure!((&field_a - &expected_a).together().is_zero());
            ensure!((&field_b - &expected_b).together().is_zero());
            let expected_density = super::super::Integrated::dimensionally_regularized(
                &(Atom::num((4, 3)) * (&field_a * &dg * &n + &field_b * (&pk - Atom::one()) * &dk)),
            );
            ensure!(
                super::super::series(&(&record.actual_numerator - expected_density), 1)?
                    .to_atom()
                    .together()
                    .is_zero(),
                "actual gauge cograph does not match independent A/B mass density"
            );
            let f = Atom::num((4, 3))
                * (Atom::num(3) / &epsilon
                    + Atom::num(4)
                    + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(4)));
            let integral = super::super::series(
                &super::super::Integrated::dimensionally_regularized(&(&field_a * f)),
                0,
            )?
            .to_atom();
            let key = (
                record.diagram.name().to_owned(),
                edges.iter().map(|edge| edge.0).collect(),
            );
            let previous = transferred.insert(key, integral.clone());
            ensure!(previous.is_none());
            aligned_poles.push(aligned);
            sum_a += &coefficient_a;
            sum_b += &coefficient_b;
            density_sum += &record.actual_numerator;
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_pure_gauge_gluon_prefix", source = %record.diagram.name(), kind,
                file.weight = %record.weight.to_canonical_string(), file.child_edges = ?edges,
                file.native_pole = %actual_source.to_canonical_string(), file.aligned_pole = %aligned_source.to_canonical_string(),
                file.coefficient_a = %coefficient_a.to_canonical_string(), file.coefficient_b = %coefficient_b.to_canonical_string(),
                file.density = %record.actual_numerator.to_canonical_string(), file.transferred_integral = %integral.to_canonical_string(),
                file.local = %cooking.uncook(record.local.atom().expression().as_view()).to_canonical_string(),
                individual_transversality_required = false,
                "Actual weighted gauge child; tensor reconstruction and cograph density precede grouped transversality"
            );
        }
        ensure!(inventory == BTreeSet::from(["ghost", "gluon", "tadpole"]));
        let expected_a = -Atom::i() * &a * Atom::num(5) / &epsilon;
        ensure!((&sum_a - &expected_a).together().is_zero());
        ensure!((&sum_b + &expected_a).together().is_zero());
        ensure!(aligned_poles.len() == 3);
        let sum = aligned_poles.iter().try_fold(
            aligned_poles[0].with_checked_expression(Atom::Zero)?,
            |sum, value| sum.try_add(value),
        )?;
        let expected = sum.with_checked_expression(
            cooking.try_cook((&expected_a * (&q_squared * &metric - &qq) * &color).as_view())?,
        )?;
        let reduced = super::super::simplify(&sum.try_sub(&expected)?)?;
        Self::certify_finite_tagged_rows_zero(
            "grouped native pure-gauge field pole is transverse",
            reduced.with_checked_expression(reduced.expression().to_cof_dimension_invariants())?,
        )?;
        let f = Atom::num((4, 3))
            * (Atom::num(3) / &epsilon
                + Atom::num(4)
                + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(4)));
        let expected_integral = super::super::series(&(-Atom::num(5) / &epsilon * f), 0)?.to_atom();
        ensure!(
            (transferred.values().sum::<Atom>() - &expected_integral)
                .together()
                .is_zero()
        );
        crate::debug_tags!(#uv, #integrated, #summary;
            stage = "native_pure_gauge_gluon_prefix_group", sources = 3,
            file.actual_density_sum = %density_sum.to_canonical_string(),
            file.transferred_integral = %expected_integral.to_canonical_string(),
            epsilon_input_depth = 1, parent_overall_operator_computed = false,
            "D00+D01+D03 retain original weights; grouped pole and analytic Ward transfer give -5F/epsilon at SU3"
        );
        Ok(transferred)
    }
}

#[test]
fn longitudinal_gluon_mass_trace_has_only_an_odd_vacuum_remainder() -> Result<()> {
    use idenso::{dirac::AGS, representations::Bispinor};
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    use vakint::vakint_symbol;

    crate::initialisation::test_initialise()?;
    let d = Atom::var(GS.dim);
    let mink = Minkowski {}.new_rep(GS.dim);
    let spin = Bispinor {}.new_rep(4);
    let spin_index = spin.to_symbolic([Atom::var(symbol!("gluon_ward::spin_index"))]);
    assert_eq!(
        super::super::simplify(
            &AppliedFeynmanRule::new(function!(ETS.metric, &spin_index, &spin_index))?.expr
        )?
        .expression(),
        &Atom::num(4),
        "Tr(1)=4 does not specialize the Lorentz dimension"
    );
    let p = function!(GS.emr_mom, 0);
    let k = function!(GS.emr_mom, 1);
    let q = function!(GS.loop_mom, 1);
    let p_vec = GS.indexed_momentum(&p, &[mink.to_symbolic([])]);
    let k_vec = GS.indexed_momentum(&k, &[mink.to_symbolic([])]);
    let q_vec = mink.vector(q.as_view(), []);
    let gamma =
        |vector: &Atom| function!(AGS.gamma, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out, vector);
    // All four finite traces in Tr[(p-slash+1) q-slash (k-slash+1) q-slash]
    // are retained. In particular the mass term is not discarded with the
    // odd words. The two powers of q use the same signed Minkowski vector.
    let four = function!(
        SPENSO_TAG.trace,
        spin.to_symbolic([]),
        gamma(&p_vec),
        gamma(&q_vec),
        gamma(&k_vec),
        gamma(&q_vec)
    );
    let odd_p = function!(
        SPENSO_TAG.trace,
        spin.to_symbolic([]),
        gamma(&p_vec),
        gamma(&q_vec),
        gamma(&q_vec)
    );
    let odd_k = function!(
        SPENSO_TAG.trace,
        spin.to_symbolic([]),
        gamma(&q_vec),
        gamma(&k_vec),
        gamma(&q_vec)
    );
    let two = function!(
        SPENSO_TAG.trace,
        spin.to_symbolic([]),
        gamma(&q_vec),
        gamma(&q_vec)
    );
    assert!(super::super::simplify(&AppliedFeynmanRule::new(odd_p.clone())?.expr)?.is_zero());
    assert!(super::super::simplify(&AppliedFeynmanRule::new(odd_k.clone())?.expr)?.is_zero());
    let spectator = function!(
        symbol!("gluon_ward::spectator"; Scalar),
        Atom::var(symbol!("gluon_ward::spectator_argument"; Scalar))
    )
    .pow(7);
    let source = &spectator * (&four + &odd_p + &odd_k + &two) / Atom::num(4);
    let traced = super::super::simplify(&AppliedFeynmanRule::new(source)?.expr)?;
    assert!(
        traced
            .expression()
            .pattern_match(&spectator.to_pattern(), None, None)
            .next()
            .is_some()
    );
    // Remove only the identical opaque spectator before the finite scalar
    // coefficient certificate. No graph numerator or denominator is expanded.
    let traced = traced
        .with_checked_expression(traced.expression().replace(spectator.to_pattern()).with(1))?;
    assert!(traced.is_scalar());
    let cooking = CookSettings::indices().with_mode(CookMode::ReversibleEncoding);
    let traced = cooking.uncook(traced.expression().as_view());
    let pk = Atom::var(symbol!("gluon_ward::pk"; Scalar));
    let kk = Atom::var(symbol!("gluon_ward::kk"; Scalar));
    let pq = Atom::var(symbol!("gluon_ward::pq"; Scalar));
    let kq = Atom::var(symbol!("gluon_ward::kq"; Scalar));
    let qq = Atom::var(symbol!("gluon_ward::qq"; Scalar));
    let independent = Kinematics::in_dimension(&d)?
        .with_scalar_product(&p, &k, pk.clone())?
        .with_scalar_product(&p, &q, pq.clone())?
        .with_scalar_product(&k, &q, kq.clone())?
        .with_scalar_product(&q, &q, qq.clone())?;
    let scalar = independent.apply(&traced);
    assert!(!scalar.contains_symbol(GS.emr_mom) && !scalar.contains_symbol(GS.loop_mom));
    let trace_formula = Atom::num(2) * &pq * &kq + (Atom::one() - &pk) * &qq;
    assert!((&scalar - trace_formula).together().is_zero());

    // Exact Gram relations for q=k-p, p^2=1. No internal d-dimensional
    // completeness relation or d=4 specialization enters this trace.
    let on_shell = Kinematics::in_dimension(&d)?
        .with_mass_squared(&p, Atom::one())?
        .with_scalar_product(&p, &k, pk.clone())?
        .with_scalar_product(&k, &k, kk.clone())?
        .with_scalar_product(&p, &q, &pk - Atom::one())?
        .with_scalar_product(&k, &q, &kk - &pk)?
        .with_scalar_product(&q, &q, &kk - Atom::num(2) * &pk + Atom::one())?;
    let routed = on_shell.apply(&traced);
    let dk = &kk - Atom::one();
    assert!((&routed - (&pk - Atom::one()) * &dk).together().is_zero());
    let pq_dot = function!(SPENSO_TAG.dot, &p_vec, &q_vec);
    assert_eq!(on_shell.apply(&pq_dot), &pk - Atom::one());
    // This is a polynomial numerator factorization. Cancelling Dk in the
    // dimensionally/analytically continued integral leaves p.q/(q^2)^2.
    // It is not the identity Dphys/(Dphys+i0)=1 at formal nonzero i0:
    // the exact finite-prescription quotient retains -i0/(Dphys+i0).
    let settings = vakint::VakintSettings {
        tensor_reduction_method: vakint::TensorReductionMethod::FeynKit,
        use_dot_product_notation: true,
        ..crate::uv::settings::VakintSettings::default().true_settings()
    };
    let vakint = crate::utils::vakint()?;
    let integrator = super::super::Integrated::new(vakint, &settings);
    for mass_squared in [
        Atom::Zero,
        Atom::var(symbol!("gluon_ward::aux_mass_squared"; Scalar)),
    ] {
        let integral = function!(
            vakint::symbols::S.topo,
            function!(vakint_symbol!("I1L"), &mass_squared, 2)
        );
        let mut odd = vakint::VakintTerm {
            integral: integral.clone(),
            numerator: super::super::Integrated::to_vakint_numerator(
                &AppliedFeynmanRule::new(&spectator * &pq_dot)?.expr,
            )?,
            vectors: vec![("k".into(), 1)],
        };
        odd.tensor_reduce(vakint, &settings)?;
        assert!(
            odd.numerator.is_zero(),
            "native odd angular moment, mass squared {mass_squared}"
        );
        assert_eq!(
            odd.integral, integral,
            "no scalar master integral was evaluated"
        );

        // A nonzero even control distinguishes odd angular projection from
        // discarding the whole scaleless integral. The projector retains its
        // exact d=4-2epsilon dependence and the same opaque powered spectator.
        let mut even = vakint::VakintTerm {
            integral: integral.clone(),
            numerator: super::super::Integrated::to_vakint_numerator(
                &AppliedFeynmanRule::new(&spectator * pq_dot.pow(2))?.expr,
            )?,
            vectors: vec![("k".into(), 1)],
        };
        even.tensor_reduce(vakint, &settings)?;
        let reduced = integrator.simplify_projected_numerator(&even.numerator, 0)?;
        assert!(!reduced.is_zero());
        assert!(
            reduced
                .pattern_match(&spectator.to_pattern(), None, None)
                .next()
                .is_some()
        );
        let expected = &spectator
            * function!(SPENSO_TAG.dot, &p_vec, &p_vec)
            * function!(SPENSO_TAG.dot, &q_vec, &q_vec)
            / &d;
        let expected = super::super::Integrated::to_vakint_numerator(
            &AppliedFeynmanRule::new(super::super::Integrated::dimensionally_regularized(
                &expected,
            ))?
            .expr,
        )?;
        assert!(
            (reduced.replace(spectator.to_pattern()).with(1)
                - expected.replace(spectator.to_pattern()).with(1))
            .together()
            .is_zero()
        );
        assert_eq!(even.integral, integral);
    }
    Ok(())
}
