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
                ensure!(edges.len() == 2 && child.lmb().loop_edges.len() == 1);
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
                ensure!(frozen.len() == 1 && frozen[0].active_components.is_empty());
                ensure!(frozen[0].frozen_lmbs() == std::slice::from_ref(child.lmb()));
                let finite = &frozen[0].atom;
                let physical_finite = integrated.physical_finite_counterterm()?;
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
    fn certify_fermion_gluon_prefixes(&self) -> Result<()> {
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
        Ok(())
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
