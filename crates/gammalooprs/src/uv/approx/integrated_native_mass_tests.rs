//! Exact native QCD mass projections and independent integrated reference checks.
use crate::{
    model::{InputParamCard, ModelGammaLoopExt},
    numerator::symbolica_ext::NumeratorAtomExt,
};
use color_eyre::{
    Result,
    eyre::{ensure, eyre},
};
use feynkit_generator::{
    FeynmanDiagramRulesExt, GenerationFilter, GenerationOptions, NumeratorGrouping, Process,
};
use feynkit_graph::expressions::evaluate_overall_factor;
use feynkit_graph::{
    DiagramEndpoint, DiagramHalfEdge, ExternalState, FeynmanDiagram, IntegralFamily, symbols,
};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::color::{CS, ColorCasimirSettings, ColorSimplifier, ColorSimplifySettings};
use idenso::dirac::GammaSimplifySettings;
use idenso::representations::{ColorAdjoint, ColorFundamental};
use idenso::tensor::{AlgebraSettings, ContractSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings};
use spenso::structure::{
    IndexedStructure, TensorStructure,
    partial::PartialStructure,
    representation::{Minkowski, RepName},
};
use std::collections::{BTreeMap, BTreeSet};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Replacement,
    symbol,
};

/// One explicit vacuum model and normalization convention for the native reference.
pub(in crate::uv) struct QcdMassReference {
    model: Model,
}

type CertifiedMassPole = (usize, Atom, Atom, Atom);

/// Completed finite projection; raw_phi retains internal MT for the derivative oracle.
struct MassProjection {
    raw_phi: Atom,
    scalar: Atom,
    family: IntegralFamily,
    colors: Vec<(Atom, Atom)>,
    // Edge identity, routed momentum, squared mass, full unit-mass polynomial.
    poles: Vec<CertifiedMassPole>,
    kernel_weight: Atom,
    numerator_prefactor: Atom,
    coupling_definitions: Vec<(String, Atom)>,
    scaleless: bool,
}

impl QcdMassReference {
    pub(in crate::uv) fn new() -> Result<Self> {
        let mut restriction: serde_json::Value = serde_json::from_str(include_str!(
            "../../../../../assets/models/json/sm/restrict_light_quarks_massless.json"
        ))?;
        for name in [
            "MU", "MD", "MS", "MC", "MB", "ymu", "ymd", "yms", "ymc", "ymb", "WT",
        ] {
            restriction[name] = serde_json::json!([0.0, 0.0]);
        }
        restriction["MT"] = serde_json::json!([1.0, 0.0]);
        for name in ["WT", "muB", "muLe", "muLmu", "muLtau", "muQ"] {
            ensure!(restriction[name] == serde_json::json!([0.0, 0.0]));
        }
        let mut model = crate::utils::load_generic_model("sm");
        model.simplify(&mut InputParamCard::from_str(
            restriction.to_string(),
            "json",
        )?)?;
        for name in ["u", "d", "s", "c", "b"] {
            ensure!(model.particle(name)?.is_massless(&model));
        }
        ensure!(!model.particle("t")?.is_massless(&model));
        Ok(Self { model })
    }

    pub(in crate::uv) fn sources(&self, loops: usize) -> Result<Vec<FeynmanDiagram>> {
        let options = GenerationOptions::default()
            .threads(1)
            .max_vertices(2 * loops)
            .with_loop_count(loops, loops)?
            .allow_self_loops(true)
            .numerator_grouping(NumeratorGrouping::None)
            .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
            .with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
                ("QCD".to_owned(), (2 * loops, Some(2 * loops))),
                ("QED".to_owned(), (0, Some(0))),
            ])));
        let result = Process::new(["t"], ["t"]).generate_diagrams(self.model.clone(), &options)?;
        result.validate_groups()?;
        ensure!(result.report.completed);
        ensure!(result.diagrams.len() == if loops == 1 { 1 } else { 12 });
        for (index, diagram) in result.diagrams.iter().enumerate() {
            diagram.validate()?;
            ensure!(diagram.loop_momentum_basis().loop_edges.len() == loops);
            let groups = result
                .groups
                .iter()
                .filter(|group| group.master == index)
                .collect::<Vec<_>>();
            ensure!(groups.len() == 1 && groups[0].members.len() == 1);
            let member = &groups[0].members[0];
            ensure!(member.source_id == diagram.id() && member.ratio.is_one());
            ensure!(member.overall_factor == *diagram.overall_factor());
        }
        Ok(result.diagrams)
    }

    fn mass_projector(
        &self,
        diagram: &FeynmanDiagram,
        d: &Atom,
    ) -> Result<(SymbolicTensor<PartialStructure>, Atom)> {
        let model = &self.model;
        let nc = Atom::var(CS.nc);
        let dimension = Kinematics::in_dimension(d)?.dimension();
        ensure!(
            diagram
                .external_state_projector()?
                .same_value(diagram.projector()),
            "native external-state owner changed"
        );
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);

        // Recover actual finalized spinor ports from the native wavefunctions.
        // The mass trace replaces the wavefunctions; it never multiplies them in.
        let mut ports = BTreeMap::new();
        // Source notation is inspected only to recover finalized ownership.
        // The actual projector and source retain their native carriers.
        let source_projector = cooking.uncook(diagram.projector().expression().as_view());
        source_projector.visitor(&mut |view| {
            if let AtomView::Fun(wave) = view
                && [symbols::u(), symbols::ubar()].contains(&wave.get_symbol())
                && wave.get_nargs() == 2
                && let Ok(edge) = usize::try_from(wave.get(0))
            {
                ports.insert(edge, (wave.get_symbol(), wave.get(1).to_owned()));
            }
            true
        });
        ensure!(
            ports.len() == 2,
            "the native source must have two spinor ports"
        );
        let mut incoming = None;
        let mut outgoing = None;
        for (edge, endpoints, data) in diagram.edges() {
            let Some(leg) = &data.external else { continue };
            ensure!(model.particle_by_id(data.particle)?.pdg_code == 6);
            let endpoint = match (endpoints.source, endpoints.target) {
                (Some(_), None) => DiagramEndpoint::Source,
                (None, Some(_)) => DiagramEndpoint::Target,
                _ => {
                    return Err(eyre!(
                        "mass reference requires unpaired amplitude external edges"
                    ));
                }
            };
            let hedge = diagram
                .half_edge_id(DiagramHalfEdge { edge, endpoint })
                .unwrap();
            let (wave, spin) = ports
                .get(&edge.0)
                .ok_or_else(|| eyre!("missing native wavefunction port"))?;
            let AtomView::Fun(rep) = spin.as_view() else {
                return Err(eyre!("spinor slot is not explicit"));
            };
            ensure!(rep.get_symbol() == symbol!("spenso::bis") && rep.get_nargs() == 2);
            ensure!(
                rep.get(0) == Atom::num(4).as_view(),
                "Tr(1)=4 remains independent of d"
            );
            let native_index = symbols::hedge_index().call((hedge.0, 1_i64));
            ensure!(
                rep.get(1) == native_index.as_view(),
                "wavefunction/finalized half-edge ownership differs"
            );
            match leg.state {
                ExternalState::Incoming => {
                    ensure!(*wave == symbols::u() && incoming.is_none());
                    incoming = Some((edge, spin.clone(), native_index));
                }
                ExternalState::Outgoing => {
                    ensure!(*wave == symbols::ubar() && outgoing.is_none());
                    outgoing = Some((edge, spin.clone(), native_index));
                }
            }
        }
        let (incoming_edge, incoming_spin, incoming_index) =
            incoming.ok_or_else(|| eyre!("incoming t missing"))?;
        let (_, outgoing_spin, outgoing_index) =
            outgoing.ok_or_else(|| eyre!("outgoing t missing"))?;
        let p_native = symbols::momentum().call(incoming_edge.0);
        let indices =
            spenso::network::parsing::ParseState::<crate::numerator::aind::Aind>::default();
        let source_numerator = cooking.uncook(diagram.numerator().expression().as_view());
        indices.reserve_indices(source_numerator.as_view());
        let lorentz = Minkowski {}
            .new_rep(dimension)
            .to_symbolic([Atom::from(indices.fresh_index())]);
        let p_vector = function!(symbols::momentum(), incoming_edge.0, &lorentz);
        let color_closure = spenso::g!(
            function!(CS.fundamental_rep, &nc, outgoing_index),
            spenso::dind!(function!(CS.fundamental_rep, &nc, incoming_index))
        );
        // External pole mass is exactly one, while internal MT stays symbolic
        // until the one-loop fixed-external mass derivative is certified.
        let mass_projector = (&p_vector
            * idenso::gamma!(incoming_spin.as_view(), outgoing_spin.as_view(), &lorentz)
            + spenso::g!(incoming_spin.as_view(), outgoing_spin.as_view()))
            * color_closure
            / (Atom::num(4) * &nc);
        Ok((
            SymbolicTensor::infer(cooking.try_cook(mass_projector.as_view())?)?,
            p_native,
        ))
    }

    fn certify_denominators(
        &self,
        diagram: &FeynmanDiagram,
        family: &IntegralFamily,
        scalar: &Atom,
        coordinates: &[(Atom, Atom)],
    ) -> Result<(Vec<CertifiedMassPole>, bool)> {
        let model = &self.model;
        let basis = diagram.loop_momentum_basis();
        let mass_symbol = symbol!("UFO::MT");
        let internal_edges = diagram
            .edges()
            .filter(|(_, _, data)| data.external.is_none() && !data.is_dummy)
            .collect::<Vec<_>>();
        ensure!(internal_edges.len() == family.denominators().len());
        let mut poles = Vec::new();
        for ((edge, _, data), denominator) in internal_edges.iter().zip(family.denominators()) {
            let particle = model.particle_by_id(data.particle)?;
            ensure!(
                particle.propagator.is_none(),
                "custom UFO denominators are outside this vacuum reference"
            );
            let width = model.parameter_by_id(particle.width)?;
            ensure!(width.name == "ZERO" || width.expression.as_ref().is_some_and(Atom::is_zero));
            let mass_squared = particle.symbolic_mass(model).pow(2);
            let momentum = basis.route_expression(&symbols::momentum().call(edge.0));
            let expected =
                family.kinematics().scalar_product(&momentum, &momentum)? - &mass_squared;
            ensure!(
                (denominator - expected).together().is_zero(),
                "native edge/pole routing mismatch"
            );
            poles.push((
                edge.0,
                momentum,
                mass_squared.replace(mass_symbol).with(1),
                denominator.replace(mass_symbol).with(1),
            ));
        }
        let mut scaleless = false;
        for (edge, endpoints, data) in &internal_edges {
            if endpoints.source.is_none()
                || endpoints.source != endpoints.target
                || !model.particle_by_id(data.particle)?.is_massless(model)
            {
                continue;
            }
            let position = basis
                .loop_edges
                .iter()
                .position(|candidate| candidate == edge)
                .ok_or_else(|| {
                    eyre!("the massless self-loop needs its native independent coordinate")
                })?;
            let loop_momentum = &family.loop_momenta()[position];
            let routed = basis.route_expression(&symbols::momentum().call(edge.0));
            ensure!(
                (&routed - loop_momentum).is_zero() || (&routed + loop_momentum).is_zero(),
                "the self-loop must be an unshifted native loop coordinate"
            );
            for (other, _, _) in &internal_edges {
                if other != edge {
                    let momentum = basis.route_expression(&symbols::momentum().call(other.0));
                    ensure!(
                        !momentum.contains(loop_momentum.as_view()),
                        "the proposed scaleless self-loop is coupled to another propagator"
                    );
                }
            }
            for (product, variable) in coordinates {
                if product.contains(loop_momentum.as_view()) {
                    ensure!(
                        !scalar.contains(variable.as_view()),
                        "isolated massless loop must not occur in the scalar numerator"
                    );
                }
            }
            let parameter_symbols = (0..family.denominators().len())
                .map(|index| symbol!(&format!("native_mass_reference::feynman_x{index}")))
                .collect::<Vec<_>>();
            let parameters = parameter_symbols
                .iter()
                .copied()
                .map(Atom::var)
                .collect::<Vec<_>>();
            let weights = family.scaleless_scaling(&parameters)?.ok_or_else(|| {
                eyre!("native parametric scaling did not certify the massless self-loop")
            })?;
            let (u, f) = family.symanzik(&parameters)?;
            let lee_pomeransky = &u + &f;
            let euler = weights
                .iter()
                .zip(&parameters)
                .zip(&parameter_symbols)
                .map(|((weight, parameter), symbol)| {
                    weight * parameter * lee_pomeransky.derivative(*symbol)
                })
                .sum::<Atom>();
            ensure!(
                (euler - &lee_pomeransky).together().is_zero(),
                "native scaleless weights must satisfy the exact parametric Euler identity"
            );
            ensure!(
                !scaleless,
                "this reference has a single isolated massless tadpole"
            );
            scaleless = true;
        }
        // This is an integral-zero certificate only: the source and its
        // nonzero local numerator remain in the generated inventory.
        Ok((poles, scaleless))
    }

    fn project(&self, diagram: &FeynmanDiagram) -> Result<MassProjection> {
        let model = &self.model;
        let loops = diagram.loop_momentum_basis().loop_edges.len();
        let name = diagram.name();
        let d = Atom::var(symbol!("native_mass_reference::d"));
        let nc = Atom::var(CS.nc);
        let cf = Atom::var(symbol!("native_mass_reference::CF"));
        let ca = Atom::var(symbol!("native_mass_reference::CA"));
        let tf = Atom::var(symbol!("native_mass_reference::TF"));
        let mass_symbol = symbol!("UFO::MT");
        let coupling = Atom::var(symbol!("UFO::G"));
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let fundamental = cooking.try_cook(function!(CS.fundamental_rep, &nc).as_view())?;
        let adjoint =
            cooking.try_cook(function!(CS.adjoint_rep, nc.pow(2) - Atom::one()).as_view())?;
        let color_rules = [
            Replacement::new(
                function!(CS.cas, 2, &fundamental).to_pattern(),
                cf.to_pattern(),
            ),
            Replacement::new(function!(CS.cas, 2, &adjoint).to_pattern(), ca.to_pattern()),
            Replacement::new(
                function!(CS.idx, 2, &fundamental).to_pattern(),
                tf.to_pattern(),
            ),
        ];
        let external_ordering = symbol!("feynkit_generator_factor::ExternalFermionOrderingSign");
        let (mass_projector, p_native) = self.mass_projector(diagram, &d)?;
        let dimensional = diagram
            .numerator()
            .map_representation_dimensions(|rep, current| {
                if rep == (Minkowski {}).into() {
                    d.clone()
                } else if rep.base() == (ColorFundamental {}).into() {
                    nc.clone()
                } else if rep == (ColorAdjoint {}).into() {
                    nc.pow(2) - Atom::one()
                } else {
                    current.to_owned()
                }
            })?;
        let mut lorentz_slots = 0usize;
        let mut wrong_dimension = false;
        let dimensional_source = cooking.uncook(dimensional.expression().as_view());
        dimensional_source.visitor(&mut |view| {
            if let AtomView::Fun(rep) = view
                && rep.get_symbol() == symbol!("spenso::mink")
            {
                lorentz_slots += 1;
                wrong_dimension |= rep.get(0) != d.as_view();
            }
            true
        });
        ensure!(
            !wrong_dimension && (lorentz_slots > 0 || dimensional.is_zero()),
            "d must reach numerator indices before contraction"
        );

        let admitted = dimensional.multiply(&mass_projector)?;
        let projected_source = cooking.uncook(admitted.expression().as_view());
        ensure!(
            cooking.uncook(cooking.try_cook(projected_source.as_view())?.as_view())
                == projected_source,
            "reversible native index/dimension admission"
        );
        ensure!(
            admitted.structure().canonical().is_scalar(),
            "mass trace must close the native spin/color interface"
        );
        let identities = admitted.simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default()),
            ..Default::default()
        })?;
        let contracted = identities.contract(ContractSettings::default())?;
        ensure!(
            contracted.contraction_complete(),
            "native mass projection must finish without re-admission"
        );

        let traced = contracted.to_dots()?;
        ensure!(traced.structure().canonical().is_scalar());
        let casimir_trace = traced
            .expression()
            .to_color_casimir_with(
                fundamental.as_view(),
                adjoint.as_view(),
                ColorCasimirSettings::default().without_fundamental_dimension_rewrite(),
            )
            .replace_multiple(&color_rules);
        let trace = cooking.uncook(casimir_trace.as_view());

        ensure!(
            !trace.contains_symbol(CS.cas) && !trace.contains_symbol(CS.idx),
            "unresolved color invariant"
        );

        let basis = diagram.loop_momentum_basis();
        let p_routed = basis.route_expression(&p_native);
        let kin = Kinematics::in_dimension(&d)?.with_mass_squared(&p_routed, Atom::one())?;
        let family = diagram.propagator_family(&kin)?;
        ensure!(family.external_momenta().len() == 1);
        let mut edge_kin = family.kinematics().clone();
        let routed_edges = diagram
            .edges()
            .map(|(edge, _, _)| {
                let native = symbols::momentum().call(edge.0);
                (edge.0, native.clone(), basis.route_expression(&native))
            })
            .collect::<Vec<_>>();
        for (_, left, left_routed) in &routed_edges {
            for (_, right, right_routed) in &routed_edges {
                edge_kin = edge_kin.with_scalar_product(
                    left,
                    right,
                    family
                        .kinematics()
                        .scalar_product(left_routed, right_routed)?,
                )?;
            }
        }
        let phi = edge_kin.apply(&trace).collect_compact_factors();
        ensure!(
            !phi.contains_symbol(symbols::momentum()),
            "edge scalar products were not completely routed"
        );
        ensure!(
            !phi.contains_symbol(CS.nc),
            "external color trace normalization did not cancel Nc"
        );

        let expanded_phi = model.expand_couplings(&phi);
        let coupling_definitions = model
            .couplings()
            .iter()
            .filter(|coupling| {
                diagram
                    .numerator()
                    .expression()
                    .contains_symbol(symbol!(&format!("UFO::{}", coupling.name)))
            })
            .map(|coupling| (coupling.name.clone(), coupling.expression.clone()))
            .collect::<Vec<_>>();
        // Convert the generated external field order to the amputated
        // Grassmann kernel exactly as the maintained FeynCalc regressions.
        let raw_factor = diagram.overall_factor();
        let mut removed = 0usize;
        let kernel_factor_source = raw_factor.replace_map(|view, _, out| {
            if let AtomView::Fun(factor) = view
                && factor.get_symbol() == external_ordering
            {
                removed += 1;
                **out = Atom::one();
            }
        });
        ensure!(removed == 1);
        let overall = evaluate_overall_factor(raw_factor.as_view());
        let kernel_factor = evaluate_overall_factor(kernel_factor_source.as_view());
        let ordering_factor = (&overall / &kernel_factor).together();
        ensure!(
            ordering_factor == Atom::num(-1),
            "native t-to-t external pairing has odd field order"
        );
        ensure!(
            evaluate_overall_factor((&kernel_factor_source * &ordering_factor).as_view())
                == overall
        );

        let coordinates = family
            .scalar_products()
            .iter()
            .enumerate()
            .map(|(index, product)| {
                (
                    product.clone(),
                    Atom::var(symbol!(&format!("native_mass_reference::sp{index}"))),
                )
            })
            .collect::<Vec<_>>();
        let coordinate_rules = coordinates
            .iter()
            .map(|(product, variable)| {
                Replacement::new(product.to_pattern(), variable.to_pattern())
            })
            .collect::<Vec<_>>();
        // i converts Γ2=-iΣ to Σ; the other i^L are the physical loop
        // measures relative to the PySecDec normalization, not graph signs.
        ensure!(
            diagram.numerator_prefactor().is_scalar(),
            "the independent mass reference requires a scalar global numerator prefactor"
        );
        let numerator_prefactor = diagram.numerator_prefactor().expression();
        let scalar = (Atom::i().pow((loops + 1) as i64)
            * &kernel_factor
            * numerator_prefactor
            * &expanded_phi
            / coupling.pow((2 * loops) as i64))
        .replace(mass_symbol)
        .with(1)
        .replace_multiple(&coordinate_rules)
        .collect_compact_factors();
        for head in [
            symbols::momentum(),
            symbols::loop_momentum(),
            symbols::external_momentum(),
            spenso::network::tags::SPENSO_TAG.dot,
            spenso::network::library::symbolic::ETS.metric,
            idenso::dirac::AGS.gamma,
            spenso::network::tags::SPENSO_TAG.trace,
            CS.t,
            CS.f,
            CS.cas,
            CS.idx,
        ] {
            ensure!(
                !scalar.contains_symbol(head),
                "unresolved tensor/momentum head in mass projection {name}: {head}"
            );
        }
        ensure!(
            !scalar.contains(coupling.as_view()),
            "coupling degree does not match loop order"
        );
        let allowed_variables = [d.clone(), cf.clone(), ca.clone(), tf.clone()]
            .into_iter()
            .chain(coordinates.iter().map(|(_, variable)| variable.clone()))
            .flat_map(|variable| variable.get_all_symbols(false))
            .collect::<BTreeSet<_>>();
        ensure!(
            scalar
                .get_all_symbols(false)
                .into_iter()
                .all(|symbol| allowed_variables.contains(&symbol)),
            "unexpected model or momentum parameter survived scalar normalization"
        );
        let color_rows = scalar
            .coefficient_list_exact(&[cf.clone(), ca.clone(), tf.clone()])
            .ok_or_else(|| eyre!("bounded exact color coefficient extraction declined"))?;
        let allowed_colors = if loops == 1 {
            vec![cf.clone()]
        } else {
            vec![cf.pow(2), &cf * &ca, &cf * &tf]
        };
        for (color, coefficient) in &color_rows {
            ensure!(
                coefficient.is_zero() || allowed_colors.contains(color),
                "unexpected native color channel"
            );
        }
        ensure!(
            (color_rows
                .iter()
                .map(|(key, coefficient)| key * coefficient)
                .sum::<Atom>()
                - &scalar)
                .together()
                .is_zero(),
            "exact color recomposition"
        );

        let (poles, scaleless) =
            self.certify_denominators(diagram, &family, &scalar, &coordinates)?;
        Ok(MassProjection {
            raw_phi: expanded_phi,
            scalar,
            family,
            colors: color_rows,
            poles,
            kernel_weight: kernel_factor,
            numerator_prefactor: numerator_prefactor.clone(),
            coupling_definitions,
            scaleless,
        })
    }
}
impl MassProjection {
    fn certify_one_loop(&self) -> Result<()> {
        let d = Atom::var(symbol!("native_mass_reference::d"));
        let cf = Atom::var(symbol!("native_mass_reference::CF"));
        let mass_symbol = symbol!("UFO::MT");
        let internal_mass = Atom::var(mass_symbol);
        let coupling = Atom::var(symbol!("UFO::G"));
        let family = &self.family;
        let p = &family.external_momenta()[0];
        let expanded_phi = self.raw_phi.clone();
        let kernel_factor = self.kernel_weight.clone();

        ensure!(family.denominators().len() == 2);
        let k = &family.loop_momenta()[0];
        let pk = family.kinematics().scalar_product(p, k)?;
        let expected = (Atom::num(2) - &d) * &pk + &d * &internal_mass;
        ensure!(
            (&expanded_phi + coupling.pow(2) * &cf * &expected)
                .together()
                .is_zero(),
            "independent one-loop d-dimensional mass numerator"
        );
        ensure!(kernel_factor == Atom::one() && self.numerator_prefactor.is_one());
        let dk = &family.denominators()[0];
        let dg = &family.denominators()[1];
        ensure!(
            (dk - family.kinematics().scalar_product(k, k)? + internal_mass.pow(2))
                .together()
                .is_zero()
        );
        ensure!(
            (dg - family.kinematics().scalar_product(&(k - p), &(k - p))?)
                .together()
                .is_zero()
        );
        let native_stripped = -&expanded_phi / (coupling.pow(2) * &cf);
        let mass_insertion = (&internal_mass * (native_stripped / dk / dg).derivative(mass_symbol))
            .replace(mass_symbol)
            .with(1);
        let dk_unit = dk.replace(mass_symbol).with(1);
        let expected_insertion = &d / (&dk_unit * dg)
            + Atom::num(2) * ((Atom::num(2) - &d) * &pk + &d) / (dk_unit.pow(2) * dg);
        ensure!(
            (&mass_insertion - &expected_insertion).together().is_zero(),
            "fixed-external internal mass derivative"
        );

        Ok(())
    }
    fn certify_rainbow(&self) -> Result<()> {
        let d = Atom::var(symbol!("native_mass_reference::d"));
        let cf = Atom::var(symbol!("native_mass_reference::CF"));
        let mass_symbol = symbol!("UFO::MT");
        let coupling = Atom::var(symbol!("UFO::G"));
        let family = &self.family;
        let p = &family.external_momenta()[0];
        let expanded_phi = self.raw_phi.clone();
        let kernel_factor = self.kernel_weight.clone();

        // The independent denominator multiset below must establish this
        // native chart: inner l=K(0), outer k=K(1).
        let l = &family.loop_momenta()[0];
        let k = &family.loop_momenta()[1];
        let kin = family.kinematics();
        let kk = kin.scalar_product(k, k)?;
        let kl = kin.scalar_product(k, l)?;
        let pk = kin.scalar_product(p, k)?;
        let pl = kin.scalar_product(p, l)?;
        let expected = (Atom::num(2) - &d).pow(2)
            * (Atom::num(2) * &kl * &pk + (Atom::one() - &kk) * &pl)
            + Atom::num(2) * &d * (Atom::num(2) - &d) * (&pk + &kl)
            + d.pow(2) * (&kk + Atom::one());
        let actual = expanded_phi.replace(mass_symbol).with(1);
        ensure!(
            (actual - Atom::i() * coupling.pow(4) * cf.pow(2) * &expected)
                .together()
                .is_zero(),
            "independent rainbow Dirac/color numerator"
        );
        let mut expected_denominators = vec![
            kin.scalar_product(l, l)? - Atom::one(),
            kin.scalar_product(&(k - l), &(k - l))?,
            kk.clone() - Atom::one(),
            kk - Atom::one(),
            kin.scalar_product(&(k - p), &(k - p))?,
        ];
        let mut actual_denominators = family
            .denominators()
            .iter()
            .map(|denominator| denominator.replace(mass_symbol).with(1))
            .collect::<Vec<_>>();
        expected_denominators.sort();
        actual_denominators.sort();
        ensure!(
            actual_denominators == expected_denominators,
            "rainbow full denominator chart and doubled quark"
        );
        ensure!(kernel_factor == Atom::one());

        Ok(())
    }
}

impl QcdMassReference {
    fn certified_projections(&self) -> Result<Vec<(FeynmanDiagram, MassProjection)>> {
        let reference = self;
        let mut projections = Vec::new();
        let cf = Atom::var(symbol!("native_mass_reference::CF"));
        let ca = Atom::var(symbol!("native_mass_reference::CA"));
        let tf = Atom::var(symbol!("native_mass_reference::TF"));
        let mut incidence = BTreeMap::<Atom, usize>::new();
        let mut light = Vec::new();
        let (mut sources, mut rainbows, mut crossed, mut zeros) = (0, 0, 0, 0);
        for loops in [1, 2] {
            for diagram in reference.sources(loops)? {
                let projected = reference.project(&diagram)?;
                ensure!(
                    !projected.scalar.is_zero(),
                    "even the scaleless tadpole has a nonzero local numerator"
                );
                sources += 1;
                zeros += usize::from(projected.scaleless);
                for (color, coefficient) in &projected.colors {
                    if !coefficient.is_zero() {
                        *incidence.entry(color.clone()).or_default() += 1;
                    }
                }
                if loops == 1 {
                    projected.certify_one_loop()?;
                }
                let abelian = projected.colors.iter().find(|(key, _)| *key == cf.pow(2));
                let nonabelian = projected.colors.iter().find(|(key, _)| *key == &cf * &ca);
                match (abelian, nonabelian) {
                    (Some(_), None) => {
                        projected.certify_rainbow()?;
                        rainbows += 1;
                    }
                    (Some((_, a)), Some((_, n))) => {
                        ensure!(
                            (a + Atom::num(2) * n).together().is_zero(),
                            "crossed color factor must be CF*(CF-CA/2)"
                        );
                        crossed += 1;
                    }
                    _ => {}
                }
                let flavors = diagram
                    .edges()
                    .filter(|(_, _, edge)| edge.external.is_none() && !edge.is_dummy)
                    .map(|(_, _, edge)| {
                        reference
                            .model
                            .particle_by_id(edge.particle)
                            .map(|p| p.pdg_code.abs())
                    })
                    .collect::<Result<BTreeSet<_>, _>>()?
                    .into_iter()
                    .filter(|pdg| (1..=5).contains(pdg))
                    .collect::<Vec<_>>();
                if !flavors.is_empty() {
                    ensure!(flavors.len() == 1);
                    light.push((flavors[0], projections.len()));
                }
                projections.push((diagram, projected));
            }
        }
        ensure!((sources, rainbows, crossed, zeros) == (13, 1, 1, 1));
        ensure!(
            incidence
                == BTreeMap::from([
                    (cf.clone(), 1),
                    (cf.pow(2), 2),
                    (&cf * &ca, 5),
                    (&cf * &tf, 6)
                ])
        );
        light.sort_by_key(|(pdg, _)| *pdg);
        ensure!(light.len() == 5 && light.iter().map(|(pdg, _)| *pdg).eq(1..=5));
        // These are five distinct native source incidences, not an extra flavor
        // multiplier. Exact equality permits later integration reuse only once.
        let (_, first) = &projections[light[0].1];
        for (_, index) in &light[1..] {
            let (_, projected) = &projections[*index];
            ensure!(projected.scalar == first.scalar);
            ensure!(projected.colors == first.colors && projected.poles == first.poles);
            ensure!(projected.family.loop_momenta() == first.family.loop_momenta());
            ensure!(projected.family.external_momenta() == first.family.external_momenta());
            ensure!(projected.family.scalar_products() == first.family.scalar_products());
            ensure!(projected.kernel_weight == first.kernel_weight);
            ensure!(projected.numerator_prefactor == first.numerator_prefactor);
            ensure!(projected.coupling_definitions == first.coupling_definitions);
        }
        Ok(projections)
    }
}

#[test]
fn native_qcd_mass_projections_certify_all_one_and_two_loop_sources() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certified_projections()?;
    Ok(())
}

impl MassProjection {
    /// Serialize only the certified finite scalar projection and its native poles.
    fn numeric_input(&self, diagram: &FeynmanDiagram, model: &Model) -> Result<serde_json::Value> {
        use serde_json::json;
        use symbolica_utils::AtomPrintExt;
        let family = &self.family;
        let momentum_names = family
            .loop_momenta()
            .iter()
            .enumerate()
            .map(|(index, momentum)| {
                (
                    momentum.clone(),
                    Atom::var(symbol!(&format!("native_mass_reference::k{index}"))),
                )
            })
            .chain(family.external_momenta().iter().map(|momentum| {
                (
                    momentum.clone(),
                    Atom::var(symbol!("native_mass_reference::p")),
                )
            }))
            .collect::<Vec<_>>();
        let momentum_rules = momentum_names
            .iter()
            .map(|(native, named)| Replacement::new(native.to_pattern(), named.to_pattern()))
            .collect::<Vec<_>>();
        let mut coordinates = Vec::new();
        let mut coordinate_rules = Vec::new();
        for (index, product) in family.scalar_products().iter().enumerate() {
            let variable = Atom::var(symbol!(&format!("native_mass_reference::sp{index}")));
            let mut pair = None;
            for (i, (left, left_name)) in momentum_names.iter().enumerate() {
                for (right, right_name) in momentum_names.iter().skip(i) {
                    if family.kinematics().scalar_product(left, right)? == *product {
                        ensure!(
                            pair.is_none(),
                            "scalar coordinate must have a unique momentum pair"
                        );
                        pair = Some((
                            left_name.to_bare_ordered_string(),
                            right_name.to_bare_ordered_string(),
                        ));
                    }
                }
            }
            let (left, right) =
                pair.ok_or_else(|| eyre!("scalar coordinate has no exact momentum pair"))?;
            coordinates.push(json!({"variable": variable.to_bare_ordered_string(), "left": left, "right": right}));
            coordinate_rules.push(Replacement::new(
                product.to_pattern(),
                variable.to_pattern(),
            ));
        }
        let cf = Atom::var(symbol!("native_mass_reference::CF"));
        let ca = Atom::var(symbol!("native_mass_reference::CA"));
        let tf = Atom::var(symbol!("native_mass_reference::TF"));
        let color_rows = self
            .colors
            .iter()
            .map(|(color, coefficient)| {
                let channel = if *color == cf {
                    "CF"
                } else if *color == cf.pow(2) {
                    "CF2"
                } else if *color == &cf * &ca {
                    "CFCA"
                } else if *color == &cf * &tf {
                    "CFTF"
                } else {
                    return Err(eyre!("uncertified color channel"));
                };
                Ok(json!({"channel": channel, "coefficient": coefficient.to_bare_ordered_string()}))
            })
            .collect::<Result<Vec<_>>>()?;
        let propagators = self.poles.iter().map(|(edge, momentum, mass_squared, polynomial)| {
            ensure!(mass_squared.is_zero() || mass_squared.is_one(), "reference masses must be zero or one");
            Ok(json!({
                "edge": edge,
                "momentum": momentum.replace_multiple(&momentum_rules).to_bare_ordered_string(),
                "mass_squared": mass_squared.to_bare_ordered_string(),
                "power": 1,
                "polynomial": polynomial.replace_multiple(&coordinate_rules).to_bare_ordered_string(),
            }))
        }).collect::<Result<Vec<_>>>()?;
        let light_flavors = diagram
            .edges()
            .filter(|(_, _, edge)| edge.external.is_none() && !edge.is_dummy)
            .map(|(_, _, edge)| {
                model
                    .particle_by_id(edge.particle)
                    .map(|p| p.pdg_code.abs())
            })
            .collect::<Result<BTreeSet<_>, _>>()?
            .into_iter()
            .filter(|pdg| (1..=5).contains(pdg))
            .collect::<Vec<_>>();
        ensure!(light_flavors.len() <= 1);
        let loops = family.loop_momenta().len();
        Ok(json!({
            "name": diagram.name(),
            "loops": loops,
            "light_flavor": light_flavors.first(),
            "integral_zero_certified": self.scaleless,
            "color_rows": color_rows,
            "loop_momenta": momentum_names[..loops].iter().map(|(_, n)| n.to_bare_ordered_string()).collect::<Vec<_>>(),
            "external_momenta": momentum_names[loops..].iter().map(|(_, n)| n.to_bare_ordered_string()).collect::<Vec<_>>(),
            "scalar_products": coordinates,
            "propagators": propagators,
        }))
    }
}

#[test]
fn native_two_loop_qcd_bare_mass_matches_primary_reference() -> Result<()> {
    use color_eyre::eyre::WrapErr;
    use std::{
        fs,
        process::Command,
        time::{SystemTime, UNIX_EPOCH},
    };
    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let sources = reference
        .certified_projections()?
        .iter()
        .map(|(diagram, projected)| projected.numeric_input(diagram, &reference.model))
        .collect::<Result<Vec<_>>>()?;
    let directory = std::env::temp_dir().join(format!(
        "gammaloop-os-mass-reference-{}-{}",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    fs::create_dir(&directory)?;
    let input = directory.join("native-sources.json");
    let script = directory.join("os_mass_reference.py");
    let output = directory.join("results");
    fs::write(
        &input,
        serde_json::to_vec_pretty(&serde_json::json!({"sources": sources}))?,
    )?;
    fs::write(
        &script,
        include_str!("../../../tests/resources/os_mass_reference.py"),
    )?;
    let python = std::env::var("PYTHON_BIN_PATH")
        .unwrap_or_else(|_| "gammaloop-os-mass-reference-python".to_owned());
    let log = fs::File::create(directory.join("backend.log"))?;
    let status = Command::new(python)
        .arg(&script)
        .arg(&input)
        .arg(&output)
        .stdout(log.try_clone()?)
        .stderr(log)
        .status()
        .wrap_err_with(|| {
            format!(
                "failed to launch pinned reference backend; inputs retained at {}",
                directory.display()
            )
        })?;
    ensure!(
        status.success(),
        "independent two-loop bare-mass reference failed; retained inputs and backend log at {}",
        directory.display()
    );
    let summary: serde_json::Value =
        serde_json::from_slice(&fs::read(output.join("summary.json"))?)?;
    ensure!(summary["one_loop_anchors_passed"] == true);
    ensure!(summary["bare_four_color_reference_passed"] == true);
    ensure!(summary["independent_insertions_reconstruct_Zm_reference"] == true);
    ensure!(summary["SU3_reference_passed"] == true);
    ensure!(summary["production_OS_forest_integrated_condition_evaluated"] == false);
    crate::debug_tags!(#uv, #integrated, #summary;
        stage = "native_two_loop_qcd_bare_mass_reference",
        file.reference = %summary,
        "Independent bare reference and insertions agree; production forest remains a separate condition"
    );
    fs::remove_dir_all(directory)?;
    Ok(())
}

impl QcdMassReference {
    // This checks valid standalone MUV children, never the guarded OS parent.
    fn certify_integrated_qg_children(
        &self,
        internal_pdgs: &[i64],
        vertex_series: &Atom,
    ) -> Result<()> {
        use crate::{
            graph::Graph,
            numerator::{aind::Aind, ufo::UFO},
            utils::GS,
            uv::{
                ApproximationType, Spinney, UVgenerationSettings, UltravioletGraph,
                approx::{
                    ForestNodeLike, Rooted, UVCtx,
                    local_4d::{Full4dCts, Local4dCts, uv_limit},
                },
                hedge_poset::OwnedForestNode,
                settings::{CTIdentifier, RenormalizationPrescriptionSettings},
            },
        };
        use linnet::half_edge::{
            involution::{Flow, Hedge},
            subgraph::{SubSetLike, SubSetOps},
        };
        use std::collections::BTreeSet;
        ensure!(internal_pdgs == [6, 6, 21] || internal_pdgs == [6, 21, 21]);
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
        // CTIdentifier signs even self-conjugate outgoing gluons. Retain
        // both native orientations, then attach through model conjugation.
        let selectors =
            [-21, 21].map(|gluon| CTIdentifier::new(BTreeSet::from([-6, 6, gluon]), None));
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let coupling = Atom::var(symbol!("UFO::G"));
        let loop_factor = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut count = 0usize;
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
            // The actual immutable model rule owns both the source vertex and the
            // tree transplanted onto the child crown. UFO owns slot/flow reindexing.
            let rule_tensor = |rule: &feynkit_model::VertexRule,
                               hedges: &[Hedge],
                               owner: u16|
             -> Result<SymbolicTensor<IndexedStructure>> {
                eyre::ensure!(rule.particles.len() == 3 && hedges.len() == 3);
                eyre::ensure!(
                    rule.color_structures.len() == 1 && rule.lorentz_structures.len() == 1
                );
                eyre::ensure!(rule.couplings.len() == 1 && rule.couplings[0].len() == 1);
                let coupling = model.coupling_by_id(
                    rule.couplings[0][0]
                        .ok_or_else(|| eyre!("canonical cubic rule needs its coupling"))?,
                )?;
                let spin_slots = hedges
                    .iter()
                    .map(|hedge| {
                        &graph.underlying[*hedge]
                            .num_indices
                            .spin_indices
                            .vertex_indices
                    })
                    .collect::<Vec<_>>();
                let color_slots = hedges
                    .iter()
                    .map(|hedge| {
                        &graph.underlying[*hedge]
                            .num_indices
                            .color_indices
                            .vertex_indices
                    })
                    .collect::<Vec<_>>();
                let momenta = hedges
                    .iter()
                    .map(|hedge| (graph.flow(*hedge), graph[hedge]))
                    .collect::<Vec<_>>();
                let spin = UFO.reindex_spin(
                    &spin_slots,
                    &momenta,
                    model
                        .lorentz_structure_by_id(rule.lorentz_structures[0])?
                        .structure
                        .clone(),
                    |index| Aind::Vertex(owner, index as u16),
                )?;
                let color =
                    UFO.reindex_color(&color_slots, rule.color_structures[0].clone(), |index| {
                        Aind::Vertex(owner, index as u16)
                    })?;
                let source = Atom::var(symbol!(&format!("UFO::{}", coupling.name))) * color * spin;
                Ok(SymbolicTensor::infer_indexed(
                    cooking.try_cook(source.as_view())?,
                )?)
            };
            let same_tensor = |difference: &SymbolicTensor<IndexedStructure>| -> Result<bool> {
                let value = difference.simplify_algebra(&idenso::tensor::AlgebraSettings {
                    gamma: Some(GammaSimplifySettings::canonical()),
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                })?;
                Ok(value
                    .coefficients_are_zero(spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
                    .is_true())
            };
            for spinney in
                graph.classified_spinneys(&uncut, model, &settings, &graph.loop_momentum_basis)?
            {
                if spinney.n_components() != 1
                    || spinney.dod != 0
                    || graph.n_loops(spinney.filter()) != 1
                    || !selectors.iter().any(|selector| {
                        selector.matches(&graph.ct_identifier(spinney.filter(), model))
                    })
                {
                    continue;
                }
                eyre::ensure!(spinney.renormalization_scheme == ApproximationType::MUV);
                let current = OwnedForestNode {
                    spinney,
                    topo_order: 1,
                };
                let crown = graph.dummy_less_full_crown(current.subgraph());
                let boundary = crown.included_iter().collect::<Vec<_>>();
                eyre::ensure!(boundary.len() == 3);
                let boundary_pdg = |hedge: Hedge| -> Result<i64> {
                    let particle = graph[graph[&hedge]]
                        .particle
                        .particle()
                        .ok_or_else(|| eyre!("physical crown particle required"))?;
                    let particle = model.particle_by_id(particle)?;
                    Ok(if graph.flow(hedge) == Flow::Source {
                        model.antiparticle(particle)?.pdg_code
                    } else {
                        particle.pdg_code
                    })
                };
                let mut actual_internal_pdgs = Vec::new();
                for (_, edge, data) in graph.iter_edges_of(current.subgraph()) {
                    eyre::ensure!(
                        data.data
                            .extra_data
                            .vakint_edge_power
                            .is_none_or(|power| power == 1)
                    );
                    let particle = model.particle_by_id(
                        data.data
                            .particle
                            .particle()
                            .ok_or_else(|| eyre!("physical propagator required"))?,
                    )?;
                    eyre::ensure!(particle.propagator.is_none());
                    actual_internal_pdgs.push(particle.pdg_code.abs());
                    if particle.pdg_code == 21 {
                        ensure!(graph[edge].particle.mass_atom(model).is_zero());
                    }
                }
                actual_internal_pdgs.sort_unstable();
                ensure!(actual_internal_pdgs == [6, 6, 21] || actual_internal_pdgs == [6, 21, 21]);
                if actual_internal_pdgs != internal_pdgs {
                    continue;
                }
                let abelian = internal_pdgs == [6, 6, 21];
                let mut qg_rules = Vec::new();
                let mut three_gluon_vertices = 0;
                let mut vertex_owners = Vec::new();
                for (node, hedges, vertex) in graph.iter_nodes_of(current.subgraph()) {
                    let rule_id = vertex
                        .vertex_rule
                        .ok_or_else(|| eyre!("actual native vertex rule required"))?;
                    let rule = model.vertex_rule_by_id(rule_id)?;
                    let pdgs = rule
                        .particles
                        .iter()
                        .map(|id| model.particle_by_id(*id).map(|particle| particle.pdg_code))
                        .collect::<Result<Vec<_>, _>>()?;
                    let mut hedges = hedges.collect::<Vec<_>>();
                    hedges.sort_by_key(|hedge| graph.underlying[*hedge].ufo_order.value);
                    eyre::ensure!(
                        hedges
                            .iter()
                            .enumerate()
                            .all(|(position, hedge)| usize::from(
                                graph.underlying[*hedge].ufo_order.value
                            ) == position)
                    );
                    eyre::ensure!(
                        same_tensor(
                            &rule_tensor(rule, &hedges, u16::try_from(node.0)?)?
                                .try_sub(&vertex.num.value)?
                        )?,
                        "actual vertex tensor differs from its immutable model rule"
                    );
                    vertex_owners.push(serde_json::json!({
                    "node":node.0, "rule_id":rule_id.index(), "rule_name":rule.name,
                    "model_order_pdgs":pdgs,
                    "hedges":hedges.iter().map(|hedge|hedge.0).collect::<Vec<_>>(),
                    "flows":hedges.iter().map(|hedge|format!("{:?}",graph.flow(*hedge))).collect::<Vec<_>>(),
                    "stored_tensor":cooking.uncook(vertex.num.value.expression().as_view()).to_canonical_string(),
                }));
                    if pdgs.iter().copied().collect::<BTreeSet<_>>() == BTreeSet::from([-6, 6, 21])
                    {
                        let anti = pdgs.iter().position(|pdg| *pdg == -6).unwrap() + 1;
                        let quark = pdgs.iter().position(|pdg| *pdg == 6).unwrap() + 1;
                        let gluon = pdgs.iter().position(|pdg| *pdg == 21).unwrap() + 1;
                        eyre::ensure!(
                            rule.color_structures == [function!(UFO.t, gluon, quark, anti)]
                        );
                        eyre::ensure!(
                            model
                                .lorentz_structure_by_id(rule.lorentz_structures[0])?
                                .structure
                                == function!(
                                    UFO.gamma,
                                    UFO.idx(gluon, 1),
                                    UFO.idx(anti, 1),
                                    UFO.idx(quark, 1)
                                )
                        );
                        qg_rules.push(rule);
                    } else {
                        eyre::ensure!(pdgs == [21, 21, 21]);
                        three_gluon_vertices += 1;
                    }
                }
                eyre::ensure!(
                    (qg_rules.len(), three_gluon_vertices) == if abelian { (3, 0) } else { (2, 1) }
                );
                let tree_rule = qg_rules[0];
                eyre::ensure!(qg_rules.iter().all(|rule| rule.name == tree_rule.name));
                let ordered_crown = tree_rule
                    .particles
                    .iter()
                    .map(|particle| {
                        let pdg = model.particle_by_id(*particle)?.pdg_code;
                        let matching = boundary
                            .iter()
                            .copied()
                            .filter_map(|hedge| {
                                boundary_pdg(hedge)
                                    .ok()
                                    .filter(|actual| *actual == pdg)
                                    .map(|_| hedge)
                            })
                            .collect::<Vec<_>>();
                        eyre::ensure!(
                            matching.len() == 1,
                            "unique native field-order crown matching"
                        );
                        Ok(matching[0])
                    })
                    .collect::<Result<Vec<_>>>()?;
                let tree = rule_tensor(tree_rule, &ordered_crown, 0)?
                    .map_representation_dimensions(|rep, current| {
                        if rep == (Minkowski {}).into() {
                            Atom::var(GS.dim)
                        } else {
                            current.to_owned()
                        }
                    })?;
                let tree_coupling = model.coupling_by_id(tree_rule.couplings[0][0].unwrap())?;
                eyre::ensure!(
                    model.expand_couplings(&Atom::var(symbol!(&format!(
                        "UFO::{}",
                        tree_coupling.name
                    )))) == Atom::i() * &coupling,
                    "the canonical amputated tree phase must be iG, not fitted from the loop"
                );
                // The shared graph traversal factors the physical source exactly.
                // This checks crown attachment without dividing by an open tensor.
                let full = graph.numerator(&graph.full_filter(), &empty)?.state.expr;
                let child = graph.numerator(current.subgraph(), &empty)?.state.expr;
                let cograph = graph
                    .numerator(&graph.full_filter(), current.subgraph())?
                    .state
                    .expr;
                eyre::ensure!(
                    child.multiply(cograph)?.try_sub(&full)?.is_zero(),
                    "exact native source/cograph factor ownership"
                );
                let ctx = UVCtx::new(&graph, model, &settings);
                let input = Full4dCts::from_factorized_local(&Local4dCts::root());
                let local = uv_limit(&input, &ctx, &current, &root, &current, &root)?;
                let backend = vakint::VakintSettings {
                    number_of_terms_in_epsilon_expansion: 3,
                    evaluation_order: vakint::EvaluationOrder(vec![
                        vakint::EvaluationMethod::AlphaLoop(vakint::AlphaLoopOptions {
                            susbstitute_masters: false,
                        }),
                    ]),
                    ..settings.vakint.true_settings()
                };
                let integrated = super::Integrated::new(crate::utils::vakint()?, &backend)
                    .run(&local, &ctx, &current, &root, &current, &root)?;
                // uv_limit stores -U0. Recombine the signed integral, then compare
                // every coefficient through epsilon^1 at mu=M=1.
                // Symbolic AlphaLoop's one-loop master already has exact pi^2;
                // no MATAD Gamma/PolyGamma identity or floating substitution is used.
                let completed = integrated
                    .physical_pole()?
                    .try_sub(&integrated.physical_finite_counterterm()?)?;
                let actual = completed.with_checked_expression(
                    model
                        .expand_couplings(completed.expression())
                        .replace(GS.m_uv_vacuum)
                        .with(1)
                        .replace(GS.m_uv_expansion)
                        .with(1)
                        .replace(GS.mu_r_sq)
                        .with(1),
                )?;
                let expected = tree.with_checked_expression(
                    -&loop_factor * model.expand_couplings(tree.expression()) * vertex_series,
                )?;
                let difference = actual.try_sub(&expected)?;
                let difference = difference.with_checked_expression(
                    super::series(difference.expression(), 1)?.to_atom(),
                )?;
                let equivalent = same_tensor(&difference)?;
                count += 1;
                records.push(serde_json::json!({
                "source":diagram.name(), "id":diagram.id().to_string(),
                "edges":graph.paired_edges(current.subgraph()).iter().map(|edge|edge.0).collect::<Vec<_>>(),
                "kind":if abelian { "Abelian" } else { "non-Abelian" },
                "internal_pdgs":actual_internal_pdgs, "vertices":vertex_owners,
                "crown":ordered_crown.iter().map(|hedge|hedge.0).collect::<Vec<_>>(),
                "lmb":current.lmb().loop_edges.iter().map(|edge|edge.0).collect::<Vec<_>>(),
                "tree":cooking.uncook(tree.expression().as_view()).to_canonical_string(),
                "actual":cooking.uncook(actual.expression().as_view()).to_canonical_string(),
                "expected":cooking.uncook(expected.expression().as_view()).to_canonical_string(), "parent_computed":false,
                "source_overall_factor_not_inserted":diagram.overall_factor().to_canonical_string(),
                "source_numerator_prefactor_not_inserted":cooking.uncook(diagram.numerator_prefactor().expression().as_view()).to_canonical_string(),
                "epsilon_terms_requested":3, "highest_compared_epsilon_power":1,
                "equivalent":equivalent,
            }));
                crate::debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_qg_vertex_child_integrated_reference",
                    file.record = %serde_json::to_string(records.last().unwrap())?,
                    "Native isolated U0 child compared through epsilon^1"
                );
                ensure!(
                    equivalent,
                    "isolated native logarithmic child differs from independent open-tree Laurent tensor: {} {:?}",
                    diagram.name(),
                    actual_internal_pdgs
                );
            }
        }
        // Native inventory must certify the two distinct attachments for each
        // triangle type; selection above never uses diagram IDs or array positions.
        assert_eq!(count, 2);
        crate::debug_tags!(#uv, #integrated, #inspect;
            stage = "native_qg_vertex_children_integrated_reference",
            file.records = %serde_json::to_string(&records)?,
            "Isolated one-loop MUV children; full OS-parent admission remains unchanged"
        );
        Ok(())
    }
}

#[test]
fn native_abelian_qg_children_match_integrated_open_tree_through_epsilon() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
    let pi = Atom::var(symbolica::atom::Symbol::PI);
    // (CF - CA/2) S_M (1-epsilon)^2 / epsilon at mu=M=1 and SU(3).
    let coefficient = Atom::num((-1, 6))
        * (epsilon.pow(-1) - Atom::num(2) + &epsilon * (Atom::one() + pi.pow(2) / Atom::num(12)));
    QcdMassReference::new()?.certify_integrated_qg_children(&[6, 6, 21], &coefficient)
}

#[test]
fn native_nonabelian_qg_children_match_integrated_open_tree_through_epsilon() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
    let pi = Atom::var(symbolica::atom::Symbol::PI);
    // CA S_M (3-2epsilon)/(2epsilon), same independent scale convention.
    let coefficient = Atom::num(3)
        * (Atom::num((3, 2)) * epsilon.pow(-1) - Atom::one() + &epsilon * pi.pow(2) / Atom::num(8));
    QcdMassReference::new()?.certify_integrated_qg_children(&[6, 21, 21], &coefficient)
}
