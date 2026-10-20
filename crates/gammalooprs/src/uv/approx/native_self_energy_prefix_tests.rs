//! Actual native D09 proper prefix, before any overall parent subtraction.
//! The causal propagators are kept separate from numerator pinch identities.

use super::*;
use crate::{
    graph::Graph,
    numerator::{AppliedFeynmanRule, aind::Aind},
    utils::{GS, W_},
    uv::{
        ApproximationType, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            ForestNodeLike, Rooted, UVCtx,
            local_4d::{FourDSector, FourDTerm, Full4dCts, Local4dCts, uv_limit},
        },
        hedge_poset::OwnedForestNode,
        settings::{CTIdentifier, RenormalizationPrescriptionSettings},
    },
};
use idenso::tensor::{SymbolicNetExt, SymbolicNetParse};
use linnet::half_edge::{
    involution::{EdgeIndex, Flow},
    subgraph::{SubSetLike, SubSetOps},
};
use spenso::network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing};

#[test]
fn native_os_quark_proper_prefix_matches_mass_and_field_insertions() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certify_os_quark_proper_prefix()?;
    Ok(())
}

// This source-only gate does not interchange an OS limit with CFF residues.
// The actual D09 prefix owner below supplies its native graph and child ports.
#[test]
fn native_d09_complete_child_has_pointwise_matched_shell_endpoints() -> Result<()> {
    use crate::{graph::parse::NumeratorInterface, uv::overlap_numerators::NumeratorJetCache};
    use idenso::shorthands::UndoShorthands;

    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let mut records = reference.certify_os_quark_proper_prefix()?;
    ensure!(records.len() == 1);
    let record = records
        .remove(&("D09".to_owned(), vec![2, 3]))
        .ok_or_else(|| eyre!("the actual D09 OS child is missing"))?;
    let graph = &record.graph;
    let child = &record.given;
    let empty = graph.empty_subgraph();
    let root = OwnedForestNode {
        spinney: Spinney::empty(graph),
        topo_order: 0,
    };
    ensure!(graph.paired_edges(child.subgraph()) == [EdgeIndex(2), EdgeIndex(3)]);
    ensure!(child.lmb().loop_edges.iter().copied().eq([EdgeIndex(2)]));
    let crown = graph.dummy_less_full_crown(child.subgraph());
    ensure!(
        crown
            .included_iter()
            .map(|hedge| graph[&hedge])
            .collect::<BTreeSet<_>>()
            == BTreeSet::from([EdgeIndex(4), EdgeIndex(5)])
    );
    // This is the same source factor and prescribed denominator product as the
    // prefix owner; the external amplitude and its ordering factor never enter.
    let source = graph
        .numerator(child.subgraph(), &empty)?
        .to_d_dim(GS.dim)?
        .color_simplify()?
        .state
        .expr;
    let bare = source.with_checked_expression(
        source.expression() * graph.denominator(child.subgraph(), &graph.model, |_| -1),
    )?;
    let frozen_kernel = function!(symbol!("native_os_prefix::pointwise_frozen_kernel"));
    let mut endpoint_count = 0;
    let mut nonzero_control_ran = false;
    for components in [["1", "0", "0", "0"], ["5/4", "1/4", "1/2", "1/2"]] {
        let mut settings = record.settings.clone();
        settings.on_shell_subtraction_velocity = components.map(str::to_owned);
        let velocity = settings.on_shell_velocity()?;
        let ctx = UVCtx::new(graph, &graph.model, &settings);
        let completed = if components == ["1", "0", "0", "0"] {
            record.completed.clone()
        } else {
            // Rebuild the actual boosted operator and its own analytic policy;
            // the rest-frame child coefficient is never substituted for it.
            let local = uv_limit(
                &Full4dCts::from_factorized_local(&Local4dCts::root()),
                &ctx,
                child,
                &root,
                child,
                &root,
            )?;
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
                .run(&local, &ctx, child, &root, child, &root)?;
            Full4dCts::recursion_input(
                &local,
                &integrated,
                ApproximationType::OS,
                false,
                child.lmb(),
            )?
        };
        let mut local_remainder = bare.clone();
        let mut frozen_remainder = bare.with_checked_expression(Atom::Zero)?;
        let mut active_count = 0;
        let mut operation_lmb = None;
        let mut frozen_count = 0;
        for sector in completed.sectors() {
            if sector.frozen_lmbs().is_empty() {
                active_count += 1;
                ensure!(sector.active_components.len() == 1);
                ensure!(&sector.active_components[0].0 == child.subgraph());
                ensure!(&sector.active_components[0].1 == child.subgraph());
                ensure!(sector.active_components[0].2.loop_edges == child.lmb().loop_edges);
                operation_lmb = Some(&sector.active_components[0].2);
                local_remainder = local_remainder.try_add(&sector.atom)?;
            } else {
                frozen_count += 1;
                ensure!(sector.active_components.is_empty());
                ensure!(sector.frozen_lmbs() == std::slice::from_ref(child.lmb()));
                // Keep the localization kernel independent and frozen. The
                // consumed-measure marker is not a physical loop coordinate.
                frozen_remainder = frozen_remainder.try_add(
                    &sector
                        .atom
                        .with_checked_expression(&frozen_kernel * sector.atom.expression())?,
                )?;
            }
        }
        ensure!((active_count, frozen_count) == (1, 1));
        let operation_lmb =
            operation_lmb.ok_or_else(|| eyre!("missing actual child operation LMB"))?;
        let complete = local_remainder.try_add(&frozen_remainder)?;
        for (label, payload) in [
            (
                "actual native bare plus signed local child",
                &local_remainder,
            ),
            ("actual native frozen analytic child", &frozen_remainder),
            ("actual complete native child", &complete),
        ] {
            // This component operation changes only the child's own boundary
            // via its retained sub-LMB. The surrounding graph is untouched.
            let endpoints = graph.on_shell_branches(
                &graph.model,
                payload,
                child.subgraph(),
                operation_lmb,
                &velocity,
            )?;
            for endpoint in endpoints {
                crate::uv::approx::local_4d::on_shell_tests::assert_rational_zero(
                    label,
                    graph,
                    child,
                    &NumeratorInterface::export(&endpoint).0,
                    crate::uv::approx::local_4d::on_shell_tests::SpinReduction::FiniteComponents,
                )?;
                endpoint_count += 1;
            }
        }
        // The same exact denominator-class oracle must distinguish the raw O
        // source. This control neither clears a forest denominator nor relies
        // on a numerical nonzero witness.
        if !nonzero_control_ran {
            let raw_endpoint = graph
                .on_shell_branches(
                    &graph.model,
                    &bare,
                    child.subgraph(),
                    operation_lmb,
                    &velocity,
                )?
                .into_iter()
                .next()
                .ok_or_else(|| eyre!("missing raw O branch"))?;
            let coefficients = crate::uv::approx::local_4d::on_shell_tests::rational_coefficients(
                graph,
                child,
                &NumeratorInterface::export(&raw_endpoint).0,
            )?;
            let mut cache = NumeratorJetCache::default();
            let mut found_nonzero = false;
            for numerator in coefficients.into_values() {
                let exact = cache.canonicalize(&numerator.undo_dots::<Aind>()?)?;
                found_nonzero |= !exact.is_zero;
            }
            ensure!(
                found_nonzero,
                "the actual raw child must fail the zero-endpoint certificate"
            );
            nonzero_control_ran = true;
        }
    }
    ensure!(endpoint_count == 12 && nonzero_control_ran);
    crate::debug_tags!(#uv, #local, #on_shell, #summary;
        stage = "native_d09_pointwise_child_endpoints",
        endpoint_count, nonzero_control_ran,
        physical_source_owners = ?[EdgeIndex(2), EdgeIndex(3)],
        boundary_owners = ?[EdgeIndex(4), EdgeIndex(5)],
        frozen_kernel_retained = true,
        cff_limit_exchange_proved = false,
        "Actual native D09 child has matched rest and boosted pointwise shell endpoints"
    );
    Ok(())
}

impl QcdMassReference {
    pub(super) fn certify_os_quark_proper_prefix(&self) -> Result<CertifiedMassPrefixes> {
        let reference = self;
        let mut one_loop_sources = reference.sources(1)?;
        let one_loop = reference.project(&one_loop_sources.remove(0))?;
        one_loop.certify_one_loop()?;
        let d = Atom::var(GS.dim);
        let epsilon = Atom::var(GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let g = Atom::var(symbol!("UFO::G"));
        let cf = Atom::num((4, 3));
        let measure = Atom::i() / (Atom::num(16) * pi.pow(2));
        let a = g.pow(2) / (Atom::num(16) * pi.pow(2));
        let f = Atom::num(3) / &epsilon
            + Atom::num(4)
            + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(4));
        let h = Atom::num(6) / &epsilon
            + Atom::num(2)
            + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(2));
        let f_color = &cf * &f;
        let h_color = &cf * &h;
        let z1 = -&f_color;
        let w1 = -&cf / &epsilon;
        // exp(gamma_E epsilon) Gamma(epsilon), at mu=m=1. Its independent
        // native one-loop normalization and epsilon-one coefficient are already
        // certified by the D0 export/reference; the numerator moments below are
        // applied to the actual generated O branches rather than this target.
        let physical_master = epsilon.pow(-1) + pi.pow(2) * &epsilon / Atom::num(12);
        let j0 = Atom::one() / (Atom::one() - Atom::num(2) * &epsilon);
        let j1 = &j0 / (Atom::num(2) - Atom::num(2) * &epsilon);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // This acts only on a certified analytic projection copy. Source ownership
        // and denominator routing are checked before the provenance is removed.
        let evaluate_vector = |tensor: &SymbolicTensor<IndexedStructure>,
                               edge: EdgeIndex,
                               value: &Atom|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let expression = GS.resolve_os_momenta(
                &GS.erase_uv_momentum_provenance(tensor.expression())
                    .replace_map(|view, _, output| {
                        if let AtomView::Fun(momentum) = view
                            && momentum.get_symbol() == GS.emr_mom
                            && momentum.get_nargs() > 0
                            && usize::try_from(momentum.get(0)).ok() == Some(edge.0)
                        {
                            let indices = momentum
                                .iter()
                                .skip(1)
                                .map(|index| index.to_owned())
                                .collect::<Vec<_>>();
                            **output = GS.indexed_momentum(value, &indices);
                        }
                    }),
            );
            Ok(tensor.with_checked_expression(expression)?)
        };
        // Open only compact vector arguments through the existing native
        // parser. Products and sums of the graph numerator stay factorized.
        let expose_vectors = |tensor: &SymbolicTensor<IndexedStructure>|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let tensor = super::super::simplify(&tensor.with_checked_expression(
                GS.erase_uv_momentum_provenance(tensor.expression()),
            )?)?;
            let expanded = cooking.uncook(tensor.expression().as_view())
                .parse_to_symbolic_net::<Aind>(&ParseSettings {
                    shorthand_parsing: ShorthandParsing::Expand {
                        schoonschip: SchoonschipExpansionMode {
                            inner_products: false,
                            expand_inside_chains: true,
                            expand_schoonship: true,
                        },
                        trace: false,
                        chain: false,
                    },
                    ..Default::default()
                })?
                .simple_execute::<()>()?;
            Ok(tensor.with_checked_expression(cooking.try_cook(expanded.as_view())?)?)
        };
        let quadratic = |momentum: &Atom, mass_squared: &Atom| {
            let represented =
                GS.indexed_momentum(momentum, &[Minkowski {}.new_rep(GS.dim).to_symbolic([])]);
            function!(
                spenso::network::tags::SPENSO_TAG.dot,
                &represented,
                &represented
            ) - mass_squared
        };
        let mut settings = UVgenerationSettings {
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
        settings.vakint.normalization = "MSbar".to_owned();
        settings.vakint.additional_normalization = "1".to_owned();
        let velocity = settings.on_shell_velocity()?;
        let mut prefixes = CertifiedMassPrefixes::new();
        for diagram in reference.sources(2)? {
            let graph = Graph::from_feynkit(&diagram, None, true)?;
            let model = &graph.model;
            let empty = graph.empty_subgraph();
            let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
            let root = OwnedForestNode {
                spinney: Spinney::empty(&graph),
                topo_order: 0,
            };
            for spinney in
                graph.classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?
            {
                if spinney.n_components() != 1
                    || spinney.dod != 1
                    || graph.n_loops(spinney.filter()) != 1
                    || spinney.renormalization_scheme != ApproximationType::OS
                {
                    continue;
                }
                ensure!(
                    graph.ct_identifier(spinney.filter(), model)
                        == CTIdentifier::new(
                            BTreeSet::from([-6, 6]),
                            Some(BTreeSet::from([6, 21]))
                        )
                );
                ensure!(
                    diagram.name() == "D09",
                    "native inventory moved its unique quark child"
                );
                let child = OwnedForestNode {
                    spinney,
                    topo_order: 1,
                };
                let child_edges = graph.paired_edges(child.subgraph());
                ensure!(
                    child_edges.iter().copied().collect::<BTreeSet<_>>()
                        == BTreeSet::from([EdgeIndex(2), EdgeIndex(3)])
                );
                ensure!(
                    child.lmb().loop_edges.iter().copied().collect::<Vec<_>>() == [EdgeIndex(2)]
                );
                let rainbow = reference.project(&diagram)?;
                rainbow.certify_rainbow(NativeQuarkTrace::Mass)?;
                ensure!(
                    rainbow.kernel_weight.is_one()
                        && (diagram.numerator_prefactor().is_scalar()
                            && diagram.numerator_prefactor().expression().is_one())
                );

                let crown = graph.dummy_less_full_crown(child.subgraph());
                let boundary = crown.included_iter().collect::<Vec<_>>();
                ensure!(boundary.len() == 2);
                let incoming = *boundary
                    .iter()
                    .find(|hedge| graph.flow(**hedge) == Flow::Sink)
                    .unwrap();
                let outgoing = *boundary
                    .iter()
                    .find(|hedge| graph.flow(**hedge) == Flow::Source)
                    .unwrap();
                let owner = graph[&incoming];
                ensure!(owner == EdgeIndex(4) && graph[&outgoing] == EdgeIndex(5));
                let mass = graph[owner].particle.mass_atom(model);
                ensure!(mass == graph[EdgeIndex(2)].particle.mass_atom(model));
                ensure!(mass == graph[EdgeIndex(5)].particle.mass_atom(model));
                let loop_vector = function!(GS.emr_mom, 2);
                let parent_vector = function!(GS.emr_mom, owner.0);
                let source_numerator = graph
                    .numerator(child.subgraph(), &empty)?
                    .to_d_dim(GS.dim)?
                    .color_simplify()?
                    .state
                    .expr;
                let bare = source_numerator.with_checked_expression(
                    source_numerator.expression()
                        * graph.denominator(child.subgraph(), model, |_| -1),
                )?;
                let ctx = UVCtx::new(&graph, model, &settings);
                let local = uv_limit(
                    &Full4dCts::from_factorized_local(&Local4dCts::root()),
                    &ctx,
                    &child,
                    &root,
                    &child,
                    &root,
                )?;
                let [active] = local.active_sectors() else {
                    return Err(eyre!("one native OS child sector required"));
                };
                ensure!(active.active_components.len() == 1 && active.frozen_lmbs().is_empty());
                ensure!(&active.active_components[0].0 == child.subgraph());
                ensure!(&active.active_components[0].1 == child.subgraph());
                ensure!(active.active_components[0].2.loop_edges == child.lmb().loop_edges);

                let branches = graph.on_shell_branches(
                    model,
                    &bare,
                    child.subgraph(),
                    child.lmb(),
                    &velocity,
                )?;
                Self::certify_finite_tagged_rows_zero(
                    "the actual local minus scheme is exactly its signed O part",
                    super::super::simplify(
                        &local
                            .atom()
                            .try_sub(active.integration_atom())?
                            .try_add(&branches[0])?
                            .try_add(&branches[1])?,
                    )?,
                )?;
                let zero = evaluate_vector(&source_numerator, EdgeIndex(2), &Atom::Zero)?;
                let mass_matrix =
                    super::super::simplify(&zero.with_checked_expression(zero.expression() / &d)?)?;
                let field_matrix = super::super::simplify(
                    &(evaluate_vector(
                        &source_numerator.with_checked_expression(
                            (source_numerator.expression() - zero.expression())
                                / (Atom::num(2) - &d),
                        )?,
                        EdgeIndex(2),
                        &parent_vector,
                    )?
                    .try_sub(&mass_matrix)?),
                )?;
                // Admit the independent one-loop J0/J1 parameter reference on
                // this actual native source. M and K retain its original phase;
                // the existing field_matrix is K-M, never an OS field counterterm.
                let kinetic_matrix = field_matrix.try_add(&mass_matrix)?;
                let complete_reference = mass_matrix.with_checked_expression(
                    &d * mass_matrix.expression()
                        + (Atom::num(2) - &d) * kinetic_matrix.expression(),
                )?;
                Self::certify_finite_tagged_rows_zero(
                    "native child complete affine numerator binds its J0/J1 matrix carriers",
                    super::super::simplify(
                        &evaluate_vector(&source_numerator, EdgeIndex(2), &parent_vector)?
                            .try_sub(&complete_reference)?,
                    )?,
                )?;
                ensure!(!mass_matrix.is_zero() && !kinetic_matrix.is_zero());
                let x = Atom::var(symbol!("native_os_prefix::feynman_x"; Scalar));
                let one_minus_x = Atom::one() - &x;
                let exposed_source = expose_vectors(&source_numerator)?;
                let exposed_zero = evaluate_vector(&exposed_source, EdgeIndex(2), &Atom::Zero)?;
                let exposed_parent =
                    evaluate_vector(&exposed_source, EdgeIndex(2), &parent_vector)?;
                let shifted = evaluate_vector(
                    &exposed_source,
                    EdgeIndex(2),
                    &(&loop_vector + &one_minus_x * &parent_vector),
                )?;
                let shifted_reference = exposed_source.try_add(
                    &exposed_parent
                        .try_sub(&exposed_zero)?
                        .with_checked_expression(
                            &one_minus_x
                                * (exposed_parent.expression() - exposed_zero.expression()),
                        )?,
                )?;
                Self::certify_finite_tagged_rows_zero(
                    "native numerator loop shift leaves only the linear odd child moment",
                    super::super::simplify(&shifted.try_sub(&shifted_reference)?)?,
                )?;

                // This positive scalar regulator is local to the original
                // child's gluon. It cannot change its Feynman-gauge numerator.
                let lambda = Atom::from(crate::model::UFOSymbol::from("MH"));
                let _parameter = model.parameter("MH")?;
                let mut regulated_graph = graph.clone();
                let gluon = &mut regulated_graph.underlying[EdgeIndex(3)];
                let original_particle = gluon.particle.particle();
                let original_numerator = gluon.num.clone();
                gluon.particle = gluon.particle.clone().override_mass(Some(lambda.clone()));
                gluon.mass = crate::graph::edge::EdgeMass::from_atom(
                    lambda.clone(),
                    model,
                    &graph.param_builder,
                )?;
                ensure!(gluon.particle.particle() == original_particle);
                ensure!(gluon.num == original_numerator);
                let regulated_source = regulated_graph
                    .numerator(child.subgraph(), &empty)?
                    .to_d_dim(GS.dim)?
                    .color_simplify()?
                    .state
                    .expr;
                ensure!(regulated_source.same_value(&source_numerator));
                let regulated_bare = regulated_source.with_checked_expression(
                    regulated_source.expression()
                        * regulated_graph.denominator(child.subgraph(), model, |_| -1),
                )?;
                let parameter_sector = FourDSector::new(
                    regulated_bare,
                    vec![(
                        child.subgraph().clone(),
                        child.subgraph().clone(),
                        child.lmb().clone(),
                    )],
                    Vec::new(),
                );
                let parameter_terms = parameter_sector.physical_terms()?;
                ensure!(parameter_terms.len() == 1);
                let poles = &parameter_terms[0].denominators;
                ensure!(poles.len() == 2);
                let quark = poles
                    .iter()
                    .find(|pole| pole.source_edge == EdgeIndex(2))
                    .ok_or_else(|| eyre!("native parameter reference lost source quark two"))?;
                let gluon = poles
                    .iter()
                    .find(|pole| pole.source_edge == EdgeIndex(3))
                    .ok_or_else(|| eyre!("native parameter reference lost source gluon three"))?;
                ensure!(quark.mass_squared == mass.pow(2));
                ensure!(gluon.mass_squared == lambda.pow(2));
                let quark_polynomial = quark.polynomial_in_lmb(child.lmb())?;
                let gluon_polynomial = gluon.polynomial_in_lmb(child.lmb())?;
                let mut reference_pole = (*quark).clone();
                reference_pole.full_expr = quadratic(&loop_vector, &mass.pow(2));
                let quark_reference = reference_pole.polynomial_in_lmb(child.lmb())?;
                let quark_identity = (&quark_polynomial - &quark_reference).expand().is_zero();
                ensure!(quark_identity);
                reference_pole.full_expr =
                    quadratic(&(&loop_vector - &parent_vector), &lambda.pow(2));
                let gluon_reference = reference_pole.polynomial_in_lmb(child.lmb())?;
                let gluon_identity = (&gluon_polynomial - &gluon_reference).expand().is_zero();
                ensure!(gluon_identity);
                let parent_squared = quadratic(&parent_vector, &Atom::Zero);
                let delta = &x * mass.pow(2) + &one_minus_x * lambda.pow(2)
                    - &x * &one_minus_x * &parent_squared;
                reference_pole.full_expr =
                    quadratic(&(&loop_vector - &one_minus_x * &parent_vector), &delta);
                // Only these bounded quadratic propagator polynomials expand.
                let completed_square = reference_pole.polynomial_in_lmb(child.lmb())?;
                let parameter_identity =
                    (&x * &quark_polynomial + &one_minus_x * &gluon_polynomial - &completed_square)
                        .expand()
                        .is_zero();
                ensure!(parameter_identity);
                let below_shell = Atom::var(symbol!("native_os_prefix::below_shell_delta"; Scalar));
                let delta_below = &x * mass.pow(2) + &one_minus_x * lambda.pow(2)
                    - &x * &one_minus_x * (mass.pow(2) - &below_shell);
                let delta_os = x.pow(2) * mass.pow(2) + &one_minus_x * lambda.pow(2);
                let below_shell_identity =
                    (&delta_below - &delta_os - &x * &one_minus_x * &below_shell)
                        .expand()
                        .is_zero();
                ensure!(below_shell_identity);
                ensure!(!(&one_minus_x * lambda.pow(2)).is_zero());

                for reference_velocity in [
                    [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
                    [
                        Atom::num((5, 4)),
                        Atom::num((1, 4)),
                        Atom::num((1, 2)),
                        Atom::num((1, 2)),
                    ],
                ] {
                    for sigma in [1, -1] {
                        let point = GS.on_shell_momentum(
                            owner,
                            &(Atom::num(sigma) * &mass),
                            &reference_velocity,
                        );
                        let endpoint = evaluate_vector(&field_matrix, owner, &point)?;
                        let raw = cooking.uncook(endpoint.expression().as_view());
                        // The source retains formal D. Embed only the helper's
                        // new Minkowski-4 projector block in that same frame.
                        let retained_frame = raw.map_mink_dim(GS.dim);
                        ensure!(retained_frame == raw);
                        // Source flow chooses row for the left product and
                        // column for the right product. No matrices commute.
                        for hedge in [outgoing, incoming] {
                            let slash =
                                crate::uv::approx::local_4d::on_shell_tests::gamma_velocity_at_port(
                                    &graph,
                                    &raw,
                                    hedge.0,
                                    &reference_velocity,
                                ).map_mink_dim(GS.dim);
                            let projected = endpoint.with_checked_expression(cooking.try_cook(
                                ((&raw + Atom::num(sigma) * slash) / Atom::num(2)).as_view(),
                            )?)?;
                            // Expose compact vectors on the complete ordered
                            // product before admitting the fixed physical block.
                            let projected = expose_vectors(&projected)?;
                            // This isolated endpoint has no active momentum or
                            // unconstrained Lorentz direction. Admit its fixed
                            // physical subspace, keeping all scalar D coefficients.
                            let raw = cooking.uncook(projected.expression().as_view());
                            ensure!(!raw.contains_symbol(GS.emr_mom));
                            ensure!(!raw.contains_symbol(GS.loop_mom));
                            ensure!(!raw.contains_symbol(GS.os_mom));
                            let mut matcher = spenso::structure::slot::SlotMatcher::default();
                            let mut lorentz = BTreeSet::new();
                            let mut fixed = BTreeSet::new();
                            let mut invalid = false;
                            raw.visitor(&mut |view| {
                                if matcher.parse::<Minkowski, Aind>(view).is_ok() {
                                    lorentz.insert(view.to_owned());
                                    return false;
                                }
                                if let AtomView::Fun(selector) = view
                                    && selector.get_symbol() == GS.delta_vec
                                    && selector.get_nargs() == 2
                                    && let Some((representation, _)) =
                                        matcher.port_representation(selector.get(1))
                                    && representation
                                        == spenso::structure::representation::LibraryRep::from(
                                            Minkowski {},
                                        )
                                {
                                    invalid |= !matches!(
                                        matcher.concrete_component(selector.get(0)),
                                        Some(Ok(0..=3))
                                    );
                                    fixed.insert(selector.get(1).to_owned());
                                }
                                true
                            });
                            ensure!(
                                !invalid && lorentz == fixed,
                                "analytic-child endpoint must contain only fixed physical Lorentz directions"
                            );
                            let physical = projected.with_checked_expression(
                                cooking.try_cook(raw.map_mink_dim(4).as_view())?,
                            )?;
                            let exact =
                                Self::materialize_native_four_component_numerator(&physical)?;
                            Self::certify_finite_tagged_rows_zero(
                                "native complete-child field carrier has both ordered endpoint zeros",
                                exact,
                            )?;
                            if sigma == 1
                                && hedge == outgoing
                                && reference_velocity[1..].iter().all(|value| value.is_zero())
                            {
                                let mass_source =
                                    cooking.uncook(mass_matrix.expression().as_view());
                                ensure!(
                                    !mass_source.contains_symbol(
                                        spenso::structure::representation::LibraryRep::from(
                                            Minkowski {}
                                        )
                                        .symbol()
                                    )
                                );
                                let shifted = physical.try_add(&mass_matrix)?;
                                let shifted =
                                    Self::materialize_native_four_component_numerator(&shifted)?;
                                let rejected = Self::certify_finite_tagged_rows_zero(
                                    "same actual endpoint rejects an added nonzero native mass carrier",
                                    shifted,
                                ).is_err();
                                ensure!(rejected);
                            }
                        }
                    }
                }
                let beta_mass = super::super::series(
                    &(&physical_master * (Atom::num(3) - Atom::num(2) * &epsilon) * &j0),
                    1,
                )?
                .to_atom();
                Self::certify_finite_tagged_rows_zero(
                    "J0/J1 endpoint beta integrals retain actual native OS mass normalization through epsilon one",
                    mass_matrix.with_checked_expression(
                        &measure * mass_matrix.expression() * (&beta_mass - &f),
                    )?,
                )?;
                let endpoint_field = &physical_master * (Atom::num(2) - &d) * &j1 + epsilon.pow(-1);
                let endpoint_field =
                    super::super::Integrated::dimensionally_regularized(&endpoint_field);
                let endpoint_field = super::super::series(&endpoint_field, 1)?.to_atom();
                let expected_field =
                    -Atom::num(2) - &epsilon * (Atom::num(4) + pi.pow(2) / Atom::num(12));
                Self::certify_finite_tagged_rows_zero(
                    "native MS-field child endpoint has the independent finite and epsilon-one parameter coefficients",
                    field_matrix.with_checked_expression(
                        &measure * field_matrix.expression() * (endpoint_field - expected_field),
                    )?,
                )?;
                crate::debug_tags!(#uv, #integrated, #summary;
                    stage = "native_os_child_parametric_source_reference",
                    original_quark_source = 2, original_gluon_source = 3,
                    feynman_parameter_denominator_exact = true,
                    numerator_shift_exact = true, source_phase_retained = true,
                    ordered_endpoint_products = 8, independent_beta_depth = 1,
                    incoming_hedge = incoming.0, outgoing_hedge = outgoing.0,
                    finite_native_endpoint_components = true, mass_shift_control = true,
                    lambda_source_local = true, unregulated_fubini_proved = false,
                    "Actual native child binds the independent J0/J1 source and endpoint reference"
                );
                let mut moment = source_numerator.with_checked_expression(Atom::Zero)?;
                for (sigma, branch) in [1, -1].into_iter().zip(branches) {
                    let sector = FourDSector::new(
                        branch.clone(),
                        vec![(
                            child.subgraph().clone(),
                            child.subgraph().clone(),
                            child.lmb().clone(),
                        )],
                        Vec::new(),
                    );
                    let _certificate = sector.canonical_projection(&graph)?;
                    let point = GS.on_shell_momentum(owner, &(Atom::num(sigma) * &mass), &velocity);
                    ensure!(
                        !branch
                            .expression()
                            .contains_symbol(GS.integrated_loop_scale)
                    );
                    for term in sector.physical_terms()? {
                        ensure!(term.denominators.len() == 2);
                        ensure!(
                            term.denominators
                                .iter()
                                .map(|den| den.source_edge)
                                .collect::<BTreeSet<_>>()
                                == BTreeSet::from([EdgeIndex(2), EdgeIndex(3)])
                        );
                        for den in &term.denominators {
                            let (momentum, mass_squared) = if den.source_edge == EdgeIndex(2) {
                                (loop_vector.clone(), mass.pow(2))
                            } else {
                                (&loop_vector - &point, Atom::Zero)
                            };
                            ensure!(den.mass_squared == mass_squared);
                            let difference = GS.resolve_os_momenta(
                                &GS.erase_uv_momentum_provenance(&den.full_expr),
                            ) - GS
                                .resolve_os_momenta(&quadratic(&momentum, &mass_squared));
                            // Keep the loop scalar product formal-d. The physical
                            // fixed point has four-component support in that frame.
                            Self::certify_finite_tagged_rows_zero(
                                "actual O denominator polynomial and prescription",
                                super::super::simplify(&AppliedFeynmanRule::new(difference)?.expr)?,
                            )?;
                        }
                        // Expose only compact vector arguments before the scalar
                        // affine/moment check, as in the independent one-loop
                        // oracle. A compact gamma(lambda*q) is not a scalar
                        // coefficient until its linear argument is exposed.
                        let numerator = expose_vectors(&term.numerator)?;
                        let at_zero = evaluate_vector(&numerator, EdgeIndex(2), &Atom::Zero)?;
                        let lambda = Atom::var(symbol!("native_os_prefix::linear_loop_scale"));
                        Self::certify_finite_tagged_rows_zero(
                            "actual O numerator is affine in its child loop",
                            super::super::simplify(
                                &evaluate_vector(
                                    &numerator,
                                    EdgeIndex(2),
                                    &(&lambda * &loop_vector),
                                )?
                                .try_sub(&at_zero)?
                                .try_sub(
                                    &numerator.with_checked_expression(
                                        &lambda * (numerator.expression() - at_zero.expression()),
                                    )?,
                                )?,
                            )?,
                        )?;
                        let at_endpoint = evaluate_vector(&numerator, EdgeIndex(2), &point)?;
                        moment = moment.try_add(&at_zero.with_checked_expression(
                            at_zero.expression() * &j0
                                + (at_endpoint.expression() - at_zero.expression()) * &j1,
                        )?)?;
                    }
                }
                let regulate = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
                    let reduced = super::super::simplify(tensor)?;
                    Ok(reduced.with_checked_expression(super::super::Integrated::dimensionally_regularized(reduced.expression()))?)
                };
                Self::certify_finite_tagged_rows_zero(
                    "actual O moment is the independent scalar mass anchor",
                    super::super::simplify(&regulate(&moment.try_sub(
                        &mass_matrix.with_checked_expression(
                            mass_matrix.expression()
                                * (Atom::num(3) - Atom::num(2) * &epsilon)
                                * &j0,
                        )?,
                    )?)?)?,
                )?;
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
                let prefix = Full4dCts::recursion_input(
                    &local,
                    &integrated,
                    ApproximationType::OS,
                    false,
                    child.lmb(),
                )?;
                let sectors = prefix.sectors().collect::<Vec<_>>();
                ensure!(sectors.len() == 2);
                let frozen = sectors
                    .iter()
                    .filter(|sector| !sector.frozen_lmbs().is_empty())
                    .collect::<Vec<_>>();
                ensure!(frozen.len() == 1 && frozen[0].active_components.is_empty());
                ensure!(frozen[0].frozen_lmbs() == std::slice::from_ref(child.lmb()));
                Self::certify_finite_tagged_rows_zero(
                    "one completed child measure and one frozen localization domain",
                    super::super::simplify(
                        &frozen[0].atom.try_sub(
                            &integrated
                                .physical_finite_counterterm()?
                                .with_checked_expression(
                                    integrated.physical_finite_counterterm()?.expression()
                                        * Atom::var(GS.integrated_loop_scale).pow(4),
                                )?,
                        )?,
                    )?,
                )?;
                let specialize = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
                    let expression = model
                        .expand_couplings(tensor.expression())
                        .replace(mass.to_pattern())
                        .with(1)
                        .replace(GS.m_uv_vacuum)
                        .with(1)
                        .replace(GS.m_uv_expansion)
                        .with(1)
                        .replace(GS.mu_r_sq)
                        .with(1)
                        .replace(GS.integrated_loop_scale)
                        .with(1);
                    Ok(tensor.with_checked_expression(expression)?)
                };
                let vacuum_integral = integrated
                    .physical_pole()?
                    .try_sub(&integrated.physical_finite_counterterm()?)?;
                let finite = specialize(&frozen[0].atom)?;
                let actual_input = moment.with_checked_expression(
                    -&measure * &physical_master * moment.expression()
                        + vacuum_integral.expression()
                        + frozen[0].atom.expression(),
                )?;
                let actual_input = specialize(&regulate(&actual_input)?)?;
                let actual = actual_input.with_checked_expression(
                    super::super::series(actual_input.expression(), 1)?.to_atom(),
                )?;
                let expected = specialize(&regulate(&mass_matrix.with_checked_expression(
                    &measure
                        * (-mass_matrix.expression() * &f + field_matrix.expression() / &epsilon),
                )?)?)?;
                let difference = actual.try_sub(&expected)?;
                Self::certify_finite_tagged_rows_zero(
                    "completed actual OS child, including mass and MS field through epsilon one",
                    super::super::simplify(&difference.with_checked_expression(
                        super::super::series(difference.expression(), 1)?.to_atom(),
                    )?)?,
                )?;
                let expected_finite =
                    specialize(&regulate(&field_matrix.with_checked_expression(
                        &measure
                            * field_matrix.expression()
                            * (Atom::one() - pi.pow(2) * &epsilon / Atom::num(12)),
                    )?)?)?;
                let difference = finite.try_sub(&expected_finite)?;
                Self::certify_finite_tagged_rows_zero(
                    "actual finite completion is independently nontrivial",
                    super::super::simplify(&difference.with_checked_expression(
                        super::super::series(difference.expression(), 1)?.to_atom(),
                    )?)?,
                )?;
                ensure!(!expected_finite.is_zero());
                let difference = actual
                    .try_sub(&finite)?
                    .try_sub(&expected)?
                    .try_add(&expected_finite)?;
                Self::certify_finite_tagged_rows_zero(
                    "omitting the actual completion leaves the known finite field remainder",
                    super::super::simplify(&difference.with_checked_expression(
                        super::super::series(difference.expression(), 1)?.to_atom(),
                    )?)?,
                )?;

                let parent = graph
                    .classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?
                    .into_iter()
                    .find(|spinney| {
                        spinney.n_components() == 1 && graph.n_loops(spinney.filter()) == 2
                    })
                    .ok_or_else(|| eyre!("native D09 parent missing"))?;
                let parent = OwnedForestNode {
                    spinney: parent,
                    topo_order: 2,
                };
                let reduced = parent.reduced_subgraph(&child);
                let outside = graph.paired_edges(&reduced);
                ensure!(
                    outside.iter().copied().collect::<BTreeSet<_>>()
                        == BTreeSet::from([EdgeIndex(4), EdgeIndex(5), EdgeIndex(6)])
                );
                for edge in &outside {
                    ensure!(
                        active.active_components[0].2.edge_signatures[*edge]
                            .internal
                            .iter()
                            .all(|sign| sign.is_zero()),
                        "untouched cograph depends on consumed child hard coordinate"
                    );
                }
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
                let uncontracted_child = graph
                    .numerator(child.subgraph(), &empty)?
                    .to_d_dim(GS.dim)?
                    .state
                    .expr;
                ensure!(
                    uncontracted_child
                        .multiply(cograph.clone())?
                        .try_sub(&full)?
                        .is_zero(),
                    "native source/cograph factor ownership"
                );
                let cograph_denominator = graph.denominator(&reduced, model, |_| -1);
                let terms = FourDTerm::from_view(cograph_denominator.as_view())?;
                ensure!(terms.len() == 1 && terms[0].denominators.len() == 3);
                ensure!(
                    terms[0]
                        .denominators
                        .iter()
                        .map(|den| den.source_edge)
                        .collect::<BTreeSet<_>>()
                        == outside.iter().copied().collect()
                );
                for den in &terms[0].denominators {
                    let native_momentum = function!(GS.emr_mom, den.source_edge.0);
                    ensure!(GS.erase_uv_momentum_provenance(&den.momentum) == native_momentum);
                    let mass_squared = if den.source_edge == EdgeIndex(6) {
                        Atom::Zero
                    } else {
                        mass.pow(2)
                    };
                    ensure!(den.mass_squared == mass_squared);
                    Self::certify_finite_tagged_rows_zero(
                        "cograph full polynomial and finite prescription",
                        super::super::simplify(
                            &AppliedFeynmanRule::new(
                                GS.erase_uv_momentum_provenance(&den.full_expr)
                                    - quadratic(&native_momentum, &mass_squared),
                            )?
                            .expr,
                        )?,
                    )?;
                }
                // Native numerator momenta, not denominator evenness, choose the
                // quark orientation. The original D09 family already certifies its
                // doubled quark and both physical child propagators independently.
                let basis = diagram.loop_momentum_basis();
                let k_native = &rainbow.family.loop_momenta()[1];
                let p_native = &rainbow.family.external_momenta()[0];
                ensure!(basis.route_expression(&function!(GS.emr_mom, 4)) == *k_native);
                ensure!(basis.route_expression(&function!(GS.emr_mom, 5)) == *k_native);
                ensure!(basis.route_expression(&function!(GS.emr_mom, 0)) == *p_native);
                let gluon = basis.route_expression(&function!(GS.emr_mom, 6));
                ensure!(gluon == k_native - p_native || gluon == p_native - k_native);
                let (projector, reference_p) = reference.mass_projector(&diagram, &d)?;
                let projector = projector.map_representation_dimensions(|_, dimension| {
                    dimension.replace(CS.nc).with(3)
                })?;
                let projector = projector
                    .with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
                let projector = projector.try_into_indexed()?;
                let pk = Atom::var(symbol!("native_os_prefix::pk"));
                let kk = Atom::var(symbol!("native_os_prefix::kk"));
                let external = function!(GS.emr_mom, 0);
                ensure!(reference_p == external);
                // The independent source projector supplies I and pslash, not a
                // produced child coefficient. Reverse each typed endpoint pair
                // to obtain the OPEN native interface, retaining an inverse.
                let zero_projector = evaluate_vector(&projector, EdgeIndex(0), &Atom::Zero)?;
                let closing_identity = zero_projector
                    .with_checked_expression(Atom::num(12) * zero_projector.expression())?;
                let closing_kinetic = projector.try_sub(&zero_projector)?;
                let closing_kinetic = closing_kinetic
                    .with_checked_expression(Atom::num(12) * closing_kinetic.expression())?;
                let open_identity = Self::open_native_quark_dual(&closing_identity, &full)?;
                let open_kinetic = Self::open_native_quark_dual(&closing_kinetic, &full)?;
                let open_loop = evaluate_vector(&open_kinetic, EdgeIndex(0), &parent_vector)?;
                // The normalized mass trace calibrates both independently
                // constructed matrices, but is not used to infer their equality.
                for matrix in [&open_identity, &open_kinetic] {
                    let trace =
                        super::super::simplify(&projector.clone().multiply(matrix.clone())?)?;
                    ensure!(trace.is_scalar());
                    let trace = Kinematics::in_dimension(&d)?
                        .with_mass_squared(&external, Atom::one())?
                        .apply(&cooking.uncook(trace.expression().as_view()));
                    ensure!((&trace - Atom::one()).together().is_zero());
                }
                let four = Atom::num(4);
                let finite_without_epsilon = finite
                    .with_checked_expression(finite.expression().replace(GS.dim_epsilon).with(0))?;
                // A complete open word keeps its formal representation. Only
                // its finite scalar coefficient participates in epsilon series.
                let truncate_open = |input: &SymbolicTensor<IndexedStructure>,
                                     spin_dimension: &Atom|
                 -> Result<SymbolicTensor<IndexedStructure>> {
                    use crate::numerator::symbolica_ext::ScalarRationalSize;
                    use spenso::structure::partial::PartialStructure;
                    let kin = Kinematics::in_dimension(spin_dimension)?
                        .with_mass_squared(&external, Atom::one())?
                        .with_scalar_product(&external, &parent_vector, pk.clone())?
                        .with_scalar_product(&parent_vector, &parent_vector, kk.clone())?;
                    let partial: SymbolicTensor<PartialStructure> = input.clone().into();
                    let rows = partial.coefficient_list(
                        spenso::shadowing::TensorCollectFilter::<0>::TaggedTensors,
                    )?;
                    let recomposed =
                        input.with_checked_expression(Atom::add_many(rows.iter().map(
                            |(word, coefficient)| word.expression() * coefficient.expression(),
                        )))?;
                    Self::certify_finite_tagged_rows_zero(
                        "D09 open coefficient decomposition retains the exact native carrier",
                        input.try_sub(&recomposed)?,
                    )?;
                    let mut contracted = input.with_checked_expression(Atom::Zero)?;
                    for (word, coefficient) in rows {
                        if word.expression().is_zero() || coefficient.expression().is_zero() {
                            continue;
                        }
                        ensure!(coefficient.is_scalar());
                        ensure!(
                            word.expression().as_view().get_byte_size() <= 8192,
                            "D09 complete open word exceeds the bounded native contraction scope"
                        );
                        let word = word
                            .simplify_algebra(&AlgebraSettings {
                                gamma: Some(GammaSimplifySettings::canonical()),
                                color: Some(
                                    ColorSimplifySettings::default()
                                        .with_cof_dimension_invariants(),
                                ),
                                ..Default::default()
                            })?
                            .contract(ContractSettings::default())?;
                        ensure!(
                            word.contraction_complete(),
                            "D09 bounded complete open word did not finish native contraction"
                        );
                        let word = word.to_dots()?.try_into_indexed()?;
                        ensure!(
                            word.structure().canonical().external_structure()
                                == input.structure().canonical().external_structure(),
                            "D09 complete word changed the native open interface"
                        );
                        let expression = kin.apply(
                            &cooking
                                .uncook(word.expression().as_view())
                                .to_cof_dimension_invariants(),
                        );
                        let word = input
                            .with_checked_expression(cooking.try_cook(expression.as_view())?)?;
                        contracted = contracted.try_add(&word.with_checked_expression(
                            word.expression() * coefficient.expression(),
                        )?)?;
                    }
                    let partial: SymbolicTensor<PartialStructure> = contracted.into();
                    let mut result = input.with_checked_expression(Atom::Zero)?;
                    for (word, coefficient) in partial.coefficient_list(
                        spenso::shadowing::TensorCollectFilter::<0>::TaggedTensors,
                    )? {
                        if word.expression().is_zero() || coefficient.expression().is_zero() {
                            continue;
                        }
                        ensure!(coefficient.is_scalar());
                        let scalar = coefficient.expression();
                        let rows = scalar.coefficient_list_exact(&[
                            d.clone(), kk.clone(), pk.clone(),
                            Atom::var(GS.dim_epsilon), Atom::var(symbolica::atom::Symbol::PI),
                        ]).ok_or_else(|| eyre!("D09 finite open coefficient exceeds the unchanged Laurent scope"))?;
                        let scalar = Atom::add_many(
                            rows.into_iter()
                                .map(|(monomial, coefficient)| monomial * coefficient),
                        );
                        ensure!(
                            !scalar.contains_symbol(GS.den)
                                && ScalarRationalSize::estimate(scalar.as_view()).is_some(),
                            "D09 open epsilon truncation requires a bounded finite scalar coefficient"
                        );
                        let scalar = super::super::series(
                            &super::super::Integrated::dimensionally_regularized(&scalar),
                            1,
                        )?
                        .to_atom();
                        // The complete word, including mink(D), is never fed to
                        // dimensional replacement or the scalar Laurent engine.
                        result = result
                            .try_add(&input.with_checked_expression(word.expression() * scalar)?)?;
                    }
                    ensure!(
                        result.structure().logical_slots() == input.structure().logical_slots()
                    );
                    Ok(result)
                };
                let mut projected_inputs = std::array::from_fn::<_, 6, _>(|_| Atom::Zero);
                let mut actual_open_inputs = Vec::new();
                for (input_index, ((child_value, spin_dimension), projected_numerator)) in [
                    (&actual, &d),
                    (&actual, &four),
                    (&finite, &d),
                    (&finite, &four),
                    (&finite_without_epsilon, &d),
                    (&finite_without_epsilon, &four),
                ]
                .into_iter()
                .zip(&mut projected_inputs)
                .enumerate()
                {
                    // Only the already integrated child's scalar epsilon series is
                    // changed in the omission controls. Spin-four changes the outer
                    // contraction metadata after child integration; the remaining
                    // scalar momentum polynomials still use the same D-dimensional
                    // Feynman-parameter integration below.
                    // This defines a formal spin-coefficient companion in the
                    // fixed kk/pk family. It does not certify the dimensional
                    // angular lift of physical-four loop momenta: that lift
                    // must separately retain projected 4/D and mixed moments.
                    let attached = specialize(&child_value.clone().multiply(cograph.clone())?)?;
                    if input_index < 2 {
                        let open = if spin_dimension == &four {
                            attached.map_representation_dimensions(|rep, dimension| {
                                if rep == (Minkowski {}).into() {
                                    Atom::num(4)
                                } else {
                                    dimension.to_owned()
                                }
                            })?
                        } else {
                            attached.clone()
                        };
                        // The certified doubled quark has Q5=Q4. Identify its
                        // vector before Gamma reduction can create scalar dots.
                        let open = open.with_checked_expression(
                            open.expression()
                                .replace(function!(GS.emr_mom, 5, W_.x___))
                                .with(function!(GS.emr_mom, 4, W_.x___)),
                        )?;
                        let open = open
                            .simplify_algebra(&AlgebraSettings {
                                gamma: Some(GammaSimplifySettings::canonical()),
                                color: Some(
                                    ColorSimplifySettings::default()
                                        .with_cof_dimension_invariants(),
                                ),
                                ..Default::default()
                            })?
                            .contract(ContractSettings::default())?;
                        ensure!(
                            open.rank() == 4
                                && open.structure().canonical().external_structure()
                                    == full.structure().canonical().external_structure(),
                            "actual D09 attachment must retain its four native spin/color ports: rank={}, status={:?}, ports={:?}",
                            open.rank(),
                            open.reduction_status(),
                            open.structure().logical_slots()
                        );
                        // Deferred retains an exact tensor. The independent
                        // complete density equality below must still prove zero;
                        // a cached planner status is not that equality oracle.
                        crate::debug_tags!(#uv, #integrated, #inspect;
                            stage = "native_d09_open_density_before_reference", input_index,
                            contraction_complete = open.contraction_complete(),
                            reduction_status = ?open.reduction_status(),
                            file.open_expression = %cooking.uncook(open.expression().as_view()).to_canonical_string(),
                            "Retained exact native open tensor before its independent density comparison"
                        );
                        let dotted = open.to_dots()?;
                        ensure!(
                            dotted.structure().logical_slots() == open.structure().logical_slots()
                        );
                        let open = dotted;
                        let expression = cooking
                            .uncook(open.expression().as_view())
                            .to_cof_dimension_invariants();
                        let kin = Kinematics::in_dimension(spin_dimension)?
                            .with_mass_squared(&external, Atom::one())?
                            .with_scalar_product(&external, &parent_vector, pk.clone())?
                            .with_scalar_product(&parent_vector, &parent_vector, kk.clone())?;
                        let expression = Atom::i()
                            * &measure
                            * &rainbow.kernel_weight
                            * diagram.numerator_prefactor().expression()
                            * kin.apply(&expression)
                            / a.pow(2);
                        let normalized =
                            open.with_checked_expression(cooking.try_cook(expression.as_view())?)?;
                        actual_open_inputs.push(truncate_open(&normalized, spin_dimension)?);
                    }
                    let closed = attached.multiply(projector.clone())?;
                    let tensor = if spin_dimension == &four {
                        closed.map_representation_dimensions(|rep, dimension| {
                            if rep == (Minkowski {}).into() {
                                Atom::num(4)
                            } else {
                                dimension.to_owned()
                            }
                        })?
                    } else {
                        closed
                    };
                    ensure!(
                        tensor.structure().canonical().is_scalar(),
                        "actual attached prefix mass trace must close all native ports"
                    );
                    let tensor = tensor
                        .simplify_algebra(&AlgebraSettings {
                            gamma: Some(GammaSimplifySettings::canonical()),
                            color: Some(
                                ColorSimplifySettings::default().with_cof_dimension_invariants(),
                            ),
                            ..Default::default()
                        })?
                        .contract(ContractSettings::default())?;
                    ensure!(
                        tensor.contraction_complete(),
                        "actual D09 finite contraction must complete"
                    );
                    let traced = cooking
                        .uncook(tensor.to_dots()?.expression().as_view())
                        .to_cof_dimension_invariants();
                    let traced = traced
                        .replace(function!(GS.emr_mom, 5, W_.x___))
                        .with(function!(GS.emr_mom, 4, W_.x___));
                    let kin = Kinematics::in_dimension(spin_dimension)?
                        .with_mass_squared(&external, Atom::one())?
                        .with_scalar_product(&external, &parent_vector, pk.clone())?
                        .with_scalar_product(&parent_vector, &parent_vector, kk.clone())?;
                    let scalar = kin.apply(&traced);
                    ensure!(
                        !scalar.contains_symbol(GS.emr_mom) && !scalar.contains_symbol(GS.os_mom)
                    );
                    // One loop is already integrated in actual. Only the remaining
                    // Minkowski measure and the single Gamma2=-i Sigma conversion enter.
                    *projected_numerator = super::super::Integrated::dimensionally_regularized(
                        &(Atom::i()
                            * &measure
                            * &rainbow.kernel_weight
                            * diagram.numerator_prefactor().expression()
                            * scalar
                            / a.pow(2)),
                    );
                }
                let [
                    actual_numerator,
                    actual_outer_four_numerator,
                    finite_d_numerator,
                    finite_four_numerator,
                    truncated_d_numerator,
                    truncated_four_numerator,
                ] = projected_inputs;
                let dk = &kk - Atom::one();
                let n = (Atom::num(2) - &d) * &pk + &d;
                let n_four = -Atom::num(2) * &pk + Atom::num(4);
                let expected_numerator = super::super::Integrated::dimensionally_regularized(
                    &(&cf * (&z1 * (&d * &dk + Atom::num(2) * &n) - &w1 * &n * &dk)),
                );
                // Compare only the finite closed numerator. Its separately certified
                // full denominator is Dk² Dg; no causal denominator is cancelled.
                ensure!(
                    super::super::series(&(&actual_numerator - &expected_numerator), 1)?
                        .to_atom()
                        .together()
                        .is_zero(),
                    "actual D09 attachment differs from fixed-external mass/field insertion density"
                );
                // Keep the already integrated D-dimensional child scalar series.
                // Only the outer native gamma contraction is four-dimensional;
                // its independent moments below still use D-dimensional measure.
                let expected_outer_four_numerator =
                    super::super::Integrated::dimensionally_regularized(
                        &(&cf
                            * (&z1 * (Atom::num(4) * &dk + Atom::num(2) * &n_four)
                                - &w1 * &n_four * &dk)),
                    );
                ensure!(
                    super::super::series(
                        &(&actual_outer_four_numerator - &expected_outer_four_numerator),
                        1,
                    )?
                    .to_atom()
                    .together()
                    .is_zero(),
                    "actual completed D09 child with outer spin four differs from the independent insertion density"
                );
                // These complete native OPEN equalities independently retain
                // the scalar and kslash coefficients. A single mass trace would
                // collapse them and miss field-like residuals.
                ensure!(actual_open_inputs.len() == 2);
                for (actual_open, spin_dimension) in actual_open_inputs.iter().zip([&d, &four]) {
                    let identity_coefficient = &cf
                        * (spin_dimension * &z1 * (&kk + Atom::one()) - spin_dimension * &w1 * &dk);
                    let kinetic_coefficient =
                        &cf * ((Atom::num(2) - spin_dimension) * (Atom::num(2) * &z1 - &w1 * &dk));
                    let reference = open_identity.with_checked_expression(
                        identity_coefficient * open_identity.expression()
                            + kinetic_coefficient * open_loop.expression(),
                    )?;
                    let reference = if spin_dimension == &four {
                        reference.map_representation_dimensions(|rep, dimension| {
                            if rep == (Minkowski {}).into() {
                                Atom::num(4)
                            } else {
                                dimension.to_owned()
                            }
                        })?
                    } else {
                        reference
                    };
                    let reference = truncate_open(&reference, spin_dimension)?;
                    Self::certify_finite_tagged_rows_zero(
                        "actual complete D09 open density has independent mass and kinetic coefficients",
                        super::super::simplify(&actual_open.try_sub(&reference)?)?,
                    )?;
                }
                // The same native D09 child/cograph/source and complete pole
                // certificates above apply to every entry. In particular this
                // is the real Dk^2 Dg cograph, not a synthetic enclosing pole.
                let unit = Atom::one();
                let child_factor = &unit - pi.pow(2) * &epsilon / Atom::num(12);
                for (label, observed, child_weight, spin_numerator) in [
                    (
                        "finite child, D spin",
                        &finite_d_numerator,
                        &child_factor,
                        &n,
                    ),
                    (
                        "finite child, four spin",
                        &finite_four_numerator,
                        &child_factor,
                        &n_four,
                    ),
                    (
                        "child epsilon omitted, D spin",
                        &truncated_d_numerator,
                        &unit,
                        &n,
                    ),
                    (
                        "child epsilon omitted, four spin",
                        &truncated_four_numerator,
                        &unit,
                        &n_four,
                    ),
                ] {
                    let expected = super::super::Integrated::dimensionally_regularized(
                        &(cf.pow(2) * child_weight * spin_numerator * &dk),
                    );
                    ensure!(
                        super::super::series(&(observed - expected), 1)?
                            .to_atom()
                            .together()
                            .is_zero(),
                        "{label}: actual frozen child differs from the native field-insertion density"
                    );
                }
                // The independently certified native one-loop family supplies
                // this integral. The numerator pinch is used only at the same
                // analytic Feynman boundary as the existing F/H transfer; no
                // finite i0 is cancelled in the retained source polynomials.
                let f_d_moment = super::super::series(
                    &super::super::Integrated::dimensionally_regularized(
                        &(&physical_master * ((Atom::num(2) - &d) * &j1 + &d * &j0)),
                    ),
                    1,
                )?
                .to_atom();
                ensure!((&f_d_moment - &f).together().is_zero());
                let f_four = super::super::series(
                    &(&physical_master * (-Atom::num(2) * &j1 + Atom::num(4) * &j0)),
                    1,
                )?
                .to_atom();
                let expected_f_four = Atom::num(3) / &epsilon
                    + Atom::num(5)
                    + &epsilon * (Atom::num(9) + pi.pow(2) / Atom::num(4));
                ensure!(
                    super::super::series(&(&f_four - expected_f_four), 1)?
                        .to_atom()
                        .together()
                        .is_zero(),
                    "spin-four uses the same D-dimensional scalar/vector moments"
                );
                let finite_d_transfer =
                    super::super::series(&(cf.pow(2) * &child_factor * &f_d_moment), 0)?.to_atom();
                let finite_four_transfer =
                    super::super::series(&(cf.pow(2) * &child_factor * &f_four), 0)?.to_atom();
                let truncated_d_transfer =
                    super::super::series(&(cf.pow(2) * &f_d_moment), 0)?.to_atom();
                let truncated_four_transfer =
                    super::super::series(&(cf.pow(2) * &f_four), 0)?.to_atom();
                for (observed, finite_part) in [
                    (&finite_d_transfer, Atom::num(4) - pi.pow(2) / Atom::num(4)),
                    (
                        &finite_four_transfer,
                        Atom::num(5) - pi.pow(2) / Atom::num(4),
                    ),
                    (&truncated_d_transfer, Atom::num(4)),
                    (&truncated_four_transfer, Atom::num(5)),
                ] {
                    let expected = cf.pow(2) * (Atom::num(3) / &epsilon + finite_part);
                    ensure!((observed - &expected).together().is_zero());
                }
                let lost_child_epsilon = -cf.pow(2) * pi.pow(2) / Atom::num(4);
                ensure!(!lost_child_epsilon.is_zero());
                ensure!(
                    (&finite_d_transfer - &truncated_d_transfer - &lost_child_epsilon)
                        .together()
                        .is_zero()
                );
                ensure!(
                    (&finite_d_transfer - &finite_four_transfer + cf.pow(2))
                        .together()
                        .is_zero()
                );
                // Differentiating only the certified outer spin coefficient
                // in H_s gives Dg/(Dk^2 Dg). At the analytic Feynman boundary
                // its tadpole moment is the SAME master, not a D=4 average.
                // Thus H_D-H_4=-2*epsilon*M. Do not cancel a retained i0.
                let h_four =
                    super::super::series(&(&h + Atom::num(2) * &epsilon * &physical_master), 1)?
                        .to_atom();
                let expected_h_four = Atom::num(6) / &epsilon
                    + Atom::num(4)
                    + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(2));
                ensure!((&h_four - &expected_h_four).together().is_zero());
                let proper_outer_four =
                    super::super::series(&(&cf * (&z1 * &h_four - &w1 * &f_four)), 0)?.to_atom();
                let expected_proper_outer_four = -cf.pow(2)
                    * (Atom::num(15) / epsilon.pow(2)
                        + Atom::num(31) / &epsilon
                        + Atom::num(79)
                        + Atom::num(11) * pi.pow(2) / Atom::num(4));
                ensure!(
                    (&proper_outer_four - &expected_proper_outer_four)
                        .together()
                        .is_zero()
                );
                let proper_d = super::super::series(&(&cf * (&z1 * &h - &w1 * &f)), 0)?.to_atom();
                let proper_difference = cf.pow(2) * (Atom::num(5) / &epsilon + Atom::num(7));
                ensure!(
                    (&proper_d - &proper_outer_four - &proper_difference)
                        .together()
                        .is_zero()
                );
                // Analytic moments act on the now-certified entire open
                // density. They retain D-dimensional integration for either
                // outer spin choice and the actual D-dimensional child series.
                let mut open_transfers = Vec::new();
                for spin_dimension in [&d, &four] {
                    let f_open = [
                        &physical_master * (Atom::num(2) - spin_dimension) * &j1,
                        &physical_master * spin_dimension * &j0,
                    ];
                    let h_open = [
                        &physical_master * (Atom::num(2) - spin_dimension) * &j0,
                        &physical_master * spin_dimension * (&j0 + Atom::one()),
                    ];
                    let mut coefficients = Vec::new();
                    for (f_coefficient, h_coefficient) in f_open.into_iter().zip(h_open) {
                        coefficients.push(
                            super::super::series(
                                &super::super::Integrated::dimensionally_regularized(
                                    &(&cf * (&z1 * h_coefficient - &w1 * f_coefficient)),
                                ),
                                0,
                            )?
                            .to_atom(),
                        );
                    }
                    open_transfers.push(coefficients);
                }
                let expected_open_d = [
                    cf.pow(2)
                        * (Atom::num(5) / epsilon.pow(2)
                            + Atom::num(12) / &epsilon
                            + Atom::num(32)
                            + Atom::num(11) * pi.pow(2) / Atom::num(12)),
                    -cf.pow(2)
                        * (Atom::num(20) / epsilon.pow(2)
                            + Atom::num(38) / &epsilon
                            + Atom::num(104)
                            + Atom::num(11) * pi.pow(2) / Atom::num(3)),
                ];
                let expected_open_four = [
                    cf.pow(2)
                        * (Atom::num(5) / epsilon.pow(2)
                            + Atom::num(17) / &epsilon
                            + Atom::num(49)
                            + Atom::num(11) * pi.pow(2) / Atom::num(12)),
                    -cf.pow(2)
                        * (Atom::num(20) / epsilon.pow(2)
                            + Atom::num(48) / &epsilon
                            + Atom::num(128)
                            + Atom::num(11) * pi.pow(2) / Atom::num(3)),
                ];
                let expected_open_difference = [
                    -cf.pow(2) * (Atom::num(5) / &epsilon + Atom::num(17)),
                    cf.pow(2) * (Atom::num(10) / &epsilon + Atom::num(24)),
                ];
                for component in 0..2 {
                    ensure!(
                        (&open_transfers[0][component] - &expected_open_d[component])
                            .together()
                            .is_zero()
                    );
                    ensure!(
                        (&open_transfers[1][component] - &expected_open_four[component])
                            .together()
                            .is_zero()
                    );
                    ensure!(
                        (&open_transfers[0][component]
                            - &open_transfers[1][component]
                            - &expected_open_difference[component])
                            .together()
                            .is_zero()
                    );
                }
                ensure!(
                    (Atom::add_many(expected_open_difference.iter()) - &proper_difference)
                        .together()
                        .is_zero()
                );
                let invisible_control =
                    super::super::simplify(&open_kinetic.try_sub(&open_identity)?)?;
                ensure!(!invisible_control.is_zero());
                ensure!(Self::certify_finite_tagged_rows_zero(
                    "open kinetic-minus-mass control must be rejected independently of its zero mass trace",
                    invisible_control.clone(),
                ).is_err());
                // A rejected zero certificate alone is not a nonzero proof.
                // This independently constructed closing dual sees the field-like
                // direction even though the normalized mass projector cannot.
                let control_dual = closing_identity.try_sub(&closing_kinetic)?;
                let dual_trace =
                    super::super::simplify(&control_dual.multiply(invisible_control.clone())?)?;
                ensure!(dual_trace.is_scalar());
                let dual_trace = Kinematics::in_dimension(&d)?
                    .with_mass_squared(&external, Atom::one())?
                    .apply(&cooking.uncook(dual_trace.expression().as_view()));
                ensure!(
                    (&dual_trace + Atom::num(24)).together().is_zero(),
                    "the independent open-control dual trace must equal -24"
                );
                ensure!(
                    (&dual_trace / Atom::num(12) + Atom::num(2))
                        .together()
                        .is_zero(),
                    "the normalized independent open-control dual trace must equal -2"
                );
                let invisible_trace = super::super::simplify(
                    &projector.clone().multiply(invisible_control.clone())?,
                )?;
                let invisible_trace = Kinematics::in_dimension(&d)?
                    .with_mass_squared(&external, Atom::one())?
                    .apply(&cooking.uncook(invisible_trace.expression().as_view()));
                ensure!(invisible_trace.together().is_zero());
                crate::debug_tags!(#uv, #integrated, #summary;
                    stage = "native_d09_actual_open_dimensional_bridge", source = %diagram.name(),
                    file.actual_open_D = %cooking.uncook(actual_open_inputs[0].expression().as_view()).to_canonical_string(),
                    file.actual_open_outer_four = %cooking.uncook(actual_open_inputs[1].expression().as_view()).to_canonical_string(),
                    file.open_D_transfers = ?open_transfers[0],
                    file.open_outer_four_transfers = ?open_transfers[1],
                    file.open_D_minus_outer_four = ?expected_open_difference,
                    both_open_coefficients_certified = true,
                    field_like_mass_invisible_control_rejected = true,
                    file.field_like_control_dual_trace = %dual_trace.to_canonical_string(),
                    field_like_control_nonzero_dual_witness = true,
                    full_d09_emitted_R2_condition_proved = false,
                    "Actual native D09 open attachment retains the independent dimensional mass and kinetic transfers"
                );
                crate::debug_tags!(#uv, #integrated, #summary;
                    stage = "native_d09_finite_child_dimensional_bridge",
                    source = %diagram.name(),
                    file.child_edges = ?child_edges, file.cograph_edges = ?outside,
                    file.finite_child = %cooking.uncook(finite.expression().as_view()).to_canonical_string(),
                    file.actual_outer_four_numerator = %actual_outer_four_numerator.to_canonical_string(),
                    file.proper_outer_four_transfer = %proper_outer_four.to_canonical_string(),
                    file.proper_D_minus_outer_four = %proper_difference.to_canonical_string(),
                    actual_completed_prefix_outer_spin_four = true,
                    file.finite_d_numerator = %finite_d_numerator.to_canonical_string(),
                    file.finite_four_numerator = %finite_four_numerator.to_canonical_string(),
                    file.finite_d_transfer = %finite_d_transfer.to_canonical_string(),
                    file.finite_four_transfer = %finite_four_transfer.to_canonical_string(),
                    file.truncated_d_transfer = %truncated_d_transfer.to_canonical_string(),
                    file.truncated_four_transfer = %truncated_four_transfer.to_canonical_string(),
                    file.lost_child_epsilon = %lost_child_epsilon.to_canonical_string(),
                    child_epsilon_depth = 1, frozen_child_unit_normalizations = 1,
                    actual_cograph_pole = true, full_d09_emitted_condition_proved = false,
                    "Actual D09 finite-child sector retains the dimensional parent and child-depth contributions"
                );

                let transferred = &z1 * &h_color - &w1 * &f_color;
                let independent = -cf.pow(2)
                    * (Atom::num(15) / epsilon.pow(2)
                        + Atom::num(26) / &epsilon
                        + Atom::num(72)
                        + Atom::num(11) * pi.pow(2) / Atom::num(4));
                ensure!(
                    super::super::series(&(&transferred - &independent), 0)?
                        .to_atom()
                        .together()
                        .is_zero()
                );
                crate::debug_tags!(#uv, #integrated, #summary;
                    stage = "native_os_quark_proper_prefix_certificate", source = %diagram.name(),
                    file.child_edges = ?child_edges, file.cograph_edges = ?outside,
                    file.completed_child = %cooking.uncook(actual.expression().as_view()).to_canonical_string(),
                    file.actual_numerator = %actual_numerator.to_canonical_string(),
                    file.reference_numerator = %expected_numerator.to_canonical_string(),
                    epsilon_input_depth = 1, parent_overall_operator_computed = false,
                    emitted_four_d_condition_proved = false,
                    "Actual D09 proper prefix transferred to independently certified F and H families"
                );
                let key = (
                    diagram.name().to_owned(),
                    child_edges.iter().map(|edge| edge.0).collect(),
                );
                ensure!(
                    prefixes
                        .insert(
                            key,
                            CertifiedMassPrefix {
                                coefficient: transferred,
                                diagram: diagram.clone(),
                                graph: graph.clone(),
                                given: OwnedForestNode {
                                    spinney: child.spinney.clone(),
                                    topo_order: child.topo_order,
                                },
                                completed: prefix.clone(),
                                settings: settings.clone(),
                            }
                        )
                        .is_none(),
                    "the actual quark prefix occurs once"
                );
            }
        }
        ensure!(prefixes.len() == 1);
        Ok(prefixes)
    }
}

#[test]
fn native_gluon_tadpole_ir_child_cancels_locally_before_its_scaleless_integral() -> Result<()> {
    use crate::uv::approx::local_4d::Local4dBranch;

    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
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
    let exact_zero = |label: &str, difference: &SymbolicTensor<IndexedStructure>| -> Result<()> {
        let tensor = super::super::simplify(difference)?;
        ensure!(
            tensor
                .coefficients_are_zero(spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
                .is_true(),
            "{label}: exact native tensor identity is unresolved"
        );
        Ok(())
    };
    let mut count = 0;
    for diagram in reference.sources(2)? {
        let graph = Graph::from_feynkit(&diagram, None, true)?;
        let model = &graph.model;
        let empty = graph.empty_subgraph();
        let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        for spinney in
            graph.classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?
        {
            if spinney.n_components() != 1
                || spinney.dod != 2
                || spinney.renormalization_scheme != ApproximationType::IR
                || graph.n_loops(spinney.filter()) != 1
                || graph.paired_edges(spinney.filter()).len() != 1
            {
                continue;
            }
            ensure!(
                diagram.name() == "D03",
                "the native inventory's isolated gluon tadpole moved"
            );
            let child = OwnedForestNode {
                spinney,
                topo_order: 1,
            };
            ensure!(
                graph.ct_identifier(child.subgraph(), model)
                    == CTIdentifier::new(BTreeSet::from([21]), Some(BTreeSet::from([21])))
            );
            let edges = graph.paired_edges(child.subgraph());
            ensure!(edges == [EdgeIndex(5)]);
            ensure!(child.lmb().loop_edges.iter().copied().collect::<Vec<_>>() == [EdgeIndex(5)]);
            let crown = graph.dummy_less_full_crown(child.subgraph());
            ensure!(
                crown
                    .included_iter()
                    .map(|hedge| graph[&hedge])
                    .collect::<BTreeSet<_>>()
                    == BTreeSet::from([EdgeIndex(3), EdgeIndex(4)])
            );
            ensure!(graph[EdgeIndex(5)].particle.mass_atom(model).is_zero());

            let raw_numerator = graph
                .numerator(child.subgraph(), &empty)?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            let numerator = graph
                .numerator(child.subgraph(), &empty)?
                .to_d_dim(GS.dim)?
                .color_simplify()?
                .state
                .expr;
            ensure!(!numerator.is_zero());
            ensure!(
                !numerator.expression().contains_symbol(GS.emr_mom)
                    && !numerator.expression().contains_symbol(GS.loop_mom)
                    && !numerator.expression().contains_symbol(GS.os_mom),
                "the actual four-gluon source has no external or loop momentum numerator"
            );
            let denominator = graph.denominator(child.subgraph(), model, |_| -1);
            let terms = FourDTerm::from_view(denominator.as_view())?;
            ensure!(
                terms.len() == 1
                    && terms[0].denominators.len() == 1
                    && (terms[0].numerator.is_scalar() && terms[0].numerator.expression().is_one())
            );
            let pole = &terms[0].denominators[0];
            let loop_momentum = function!(GS.emr_mom, 5);
            ensure!(pole.source_edge == EdgeIndex(5) && pole.mass_squared.is_zero());
            ensure!(GS.erase_uv_momentum_provenance(&pole.momentum) == loop_momentum);
            exact_zero(
                "actual massless tadpole polynomial and prescription",
                &AppliedFeynmanRule::new(
                    GS.erase_uv_momentum_provenance(&pole.full_expr)
                        - Minkowski {}
                            .new_rep(GS.dim)
                            .inner_product(&loop_momentum, &loop_momentum),
                )?
                .expr,
            )?;
            let bare = numerator.with_checked_expression(numerator.expression() * &denominator)?;

            // The independently certified source is boundary-independent:
            // its quartic numerator has no momentum and its only pole is the
            // unshifted child coordinate. Thus S1 acts as the identity. Execute
            // the full native U2+S1-U2S1 below rather than replacing that operator
            // with its expected answer. The private soft-rescaling owner stays
            // private; this test supplies the independent source-level proof.
            ensure!(
                child.lmb().edge_signatures[EdgeIndex(5)]
                    .external
                    .iter()
                    .all(|sign| sign.is_zero())
            );
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
                return Err(eyre!("retain one actual, nonzero tadpole sector"));
            };
            ensure!(!sector.atom.is_zero() && sector.frozen_lmbs().is_empty());
            ensure!(
                sector.active_components.len() == 1
                    && &sector.active_components[0].0 == child.subgraph()
                    && &sector.active_components[0].1 == child.subgraph()
            );
            let branches = local
                .provenance_branches()
                .iter()
                .map(|branch| branch.branch)
                .collect::<Vec<_>>();
            ensure!(branches.len() == 4);
            for branch in [
                Local4dBranch::U,
                Local4dBranch::S,
                Local4dBranch::US,
                Local4dBranch::Combined,
            ] {
                ensure!(branches.contains(&branch));
            }
            let _certificate = sector.canonical_projection(&graph)?;
            let local_projection = local.atom().with_checked_expression(
                GS.erase_uv_momentum_provenance(local.atom().expression()),
            )?;
            exact_zero(
                "actual U2+S1-U2S1 with one forest minus cancels the bare tadpole",
                &local_projection.try_add(&bare)?,
            )?;
            // Bare IR policy integrates this massless S term. Do not replace its
            // actual nonzero input by zero merely because its integral vanishes.
            ensure!(!sector.integration_atom().is_zero());
            exact_zero(
                "stored analytic input retains the signed physical tadpole",
                &sector.integration_atom().try_sub(local.atom())?,
            )?;

            let parent = graph
                .classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?
                .into_iter()
                .find(|spinney| spinney.n_components() == 1 && graph.n_loops(spinney.filter()) == 2)
                .ok_or_else(|| eyre!("native D03 parent missing"))?;
            let parent = OwnedForestNode {
                spinney: parent,
                topo_order: 2,
            };
            let reduced = parent.reduced_subgraph(&child);
            let outside = graph.paired_edges(&reduced);
            ensure!(
                outside.iter().copied().collect::<BTreeSet<_>>()
                    == BTreeSet::from([EdgeIndex(2), EdgeIndex(3), EdgeIndex(4)])
            );
            for edge in &outside {
                ensure!(
                    child.lmb().edge_signatures[*edge]
                        .internal
                        .iter()
                        .all(|sign| sign.is_zero())
                );
            }
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
                raw_numerator
                    .multiply(cograph.clone())?
                    .try_sub(&full)?
                    .is_zero(),
                "retain actual source/cograph factors"
            );
            let cograph_denominator = graph.denominator(&reduced, model, |_| -1);
            let complete_cancellation = local_projection.try_add(&bare)?.multiply(cograph)?;
            exact_zero(
                "D03 source plus actual proper prefix cancels before integration",
                &complete_cancellation.with_checked_expression(
                    complete_cancellation.expression() * &cograph_denominator,
                )?,
            )?;

            // Only after all local identities: run the existing independent
            // family scaling certificate and the actual analytic child policy.
            let projected = reference.project(&diagram)?;
            ensure!(
                projected.scaleless && !projected.scalar.is_zero() && !projected.raw_phi.is_zero()
            );
            ensure!(
                projected.kernel_weight == Atom::num((1, 2))
                    && (diagram.numerator_prefactor().is_scalar()
                        && diagram.numerator_prefactor().expression().is_one())
            );
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
            ensure!(integrated.physical_pole()?.is_zero());
            ensure!(integrated.physical_finite_counterterm()?.is_zero());
            let prefix = Full4dCts::recursion_input(
                &local,
                &integrated,
                ApproximationType::IR,
                false,
                child.lmb(),
            )?;
            ensure!(prefix.atom().same_value(local.atom()));
            for completion in prefix
                .sectors()
                .filter(|sector| !sector.frozen_lmbs().is_empty())
            {
                ensure!(completion.atom.is_zero() && completion.active_components.is_empty());
                ensure!(completion.frozen_lmbs() == std::slice::from_ref(child.lmb()));
            }
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_gluon_tadpole_proper_prefix_certificate", source = %diagram.name(),
                file.source_overall_factor = %diagram.overall_factor().to_canonical_string(),
                file.bare = %cooking.uncook(bare.expression().as_view()).to_canonical_string(), file.local = %cooking.uncook(local.atom().expression().as_view()).to_canonical_string(),
                file.integration_input = %cooking.uncook(sector.integration_atom().expression().as_view()).to_canonical_string(),
                local_cancellation_before_integration = true, integrated_prefix_zero = true,
                parent_overall_operator_computed = false,
                "Nonzero native D03 retained through its actual local and analytic operators"
            );
            count += 1;
        }
    }
    ensure!(count == 1);
    Ok(())
}
