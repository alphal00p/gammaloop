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
                rainbow.certify_rainbow()?;
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
                        let numerator =
                            super::super::simplify(&term.numerator.with_checked_expression(
                                GS.erase_uv_momentum_provenance(term.numerator.expression()),
                            )?)?;
                        let expanded = cooking
                            .uncook(numerator.expression().as_view())
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
                        let numerator = numerator
                            .with_checked_expression(cooking.try_cook(expanded.as_view())?)?;
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
                let (projector, _) = reference.mass_projector(&diagram, &d)?;
                let projector = projector.map_representation_dimensions(|_, dimension| {
                    dimension.replace(CS.nc).with(3)
                })?;
                let projector = projector
                    .with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
                let projector = projector.try_into_indexed()?;
                let pk = Atom::var(symbol!("native_os_prefix::pk"));
                let kk = Atom::var(symbol!("native_os_prefix::kk"));
                let external = function!(GS.emr_mom, 0);
                let four = Atom::num(4);
                let finite_without_epsilon = finite
                    .with_checked_expression(finite.expression().replace(GS.dim_epsilon).with(0))?;
                let mut projected_inputs = std::array::from_fn::<_, 5, _>(|_| Atom::Zero);
                for ((child_value, spin_dimension), projected_numerator) in [
                    (&actual, &d),
                    (&finite, &d),
                    (&finite, &four),
                    (&finite_without_epsilon, &d),
                    (&finite_without_epsilon, &four),
                ]
                .into_iter()
                .zip(&mut projected_inputs)
                {
                    // Only the already integrated child's scalar epsilon series is
                    // changed in the omission controls. Spin-four changes the outer
                    // contraction metadata after child integration; the remaining
                    // scalar momentum polynomials still use the same D-dimensional
                    // Feynman-parameter integration below.
                    let closed = specialize(
                        &child_value
                            .clone()
                            .multiply(cograph.clone())?
                            .multiply(projector.clone())?,
                    )?;
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
                    finite_d_numerator,
                    finite_four_numerator,
                    truncated_d_numerator,
                    truncated_four_numerator,
                ] = projected_inputs;
                let dk = &kk - Atom::one();
                let n = (Atom::num(2) - &d) * &pk + &d;
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
                // The same native D09 child/cograph/source and complete pole
                // certificates above apply to every entry. In particular this
                // is the real Dk^2 Dg cograph, not a synthetic enclosing pole.
                let unit = Atom::one();
                let child_factor = &unit - pi.pow(2) * &epsilon / Atom::num(12);
                let n_four = -Atom::num(2) * &pk + Atom::num(4);
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
                crate::debug_tags!(#uv, #integrated, #summary;
                    stage = "native_d09_finite_child_dimensional_bridge",
                    source = %diagram.name(),
                    file.child_edges = ?child_edges, file.cograph_edges = ?outside,
                    file.finite_child = %cooking.uncook(finite.expression().as_view()).to_canonical_string(),
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
