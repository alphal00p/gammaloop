//! Numerical integration of the production-emitted, complete native OS forest.
//! The temporal closing word is attached by the ordinary final cograph owner.

use super::*;
use crate::{
    graph::{Graph, edge::EdgeMass, feynman_graph::FeynmanGraph},
    integrands::process::{
        ParamBuilder,
        evaluators::{EvaluatorStack, GenericEvaluator, GenericEvaluatorFloat},
    },
    model::UFOSymbol,
    numerator::AppliedFeynmanRule,
    settings::global::GenerationSettings,
    utils::{ArbPrec, F, GS},
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings,
        UltravioletGraph,
        approx::{
            CutStructure, OrientationProjection, local_4d::FourDSector,
            projected_4d::Local4dProjectionContext,
        },
        hedge_poset::Wood,
        orchestrator::validate_local_counterterm_schemes,
        settings::FinalIntegrandDimension,
        uv_graph::UVE,
    },
};
use idenso::{dirac::AGS, representations::Bispinor};
use linnet::half_edge::subgraph::SubSetOps;
use spenso::structure::representation::LibraryRep;
use spenso::{
    algebra::complex::Complex,
    network::{ExecutionResult, Sequential, SmallestDegree, library::symbolic::ETS},
};
use symbolica::numerical_integration::{
    ContinuousGrid, MonteCarloRng, Sample, StatisticsAccumulator,
};
use symbolica_utils::AtomPrintExt;

struct EmittedVectorIntegral {
    evaluator: GenericEvaluator,
    values: Vec<Complex<F<f64>>>,
    coordinates: Vec<[usize; 3]>,
    normalization: f64,
    source: serde_json::Value,
    evaluations: usize,
    precise_evaluations: usize,
    precise_seconds: f64,
    uv_mass_slots: [usize; 2],
    localization_slot: usize,
    physical_slots: [usize; 3],
    gluon_mass: f64,
}

impl EmittedVectorIntegral {
    fn new(reference: &QcdMassReference, loops: usize) -> Result<Self> {
        let started = std::time::Instant::now();
        let diagram = reference
            .sources(loops)?
            .into_iter()
            .find(|source| loops == 1 || source.name() == "D09")
            .ok_or_else(|| eyre!("native rainbow source missing"))?;
        let projection = reference.project_quark_trace(&diagram, NativeQuarkTrace::Vector)?;
        ensure!(projection.kernel_weight.is_one() && projection.numerator_prefactor.is_one());
        ensure!(projection.external_ordering == Atom::num(-1));
        let mut graph = Graph::from_feynkit(&diagram, None, true)?;
        if loops == 2 {
            ensure!(
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .eq([2, 4]),
                "the D09 soft charts require the certified native inner l=e2 and outer k=e4 basis"
            );
            projection.certify_rainbow(NativeQuarkTrace::Vector)?;
        }
        let full = graph.full_filter();
        let empty = graph.empty_subgraph();
        let original_numerator = graph.numerator(&full, &empty)?.state.expr;
        let original_bare = original_numerator.with_checked_expression(
            original_numerator.expression() * graph.denominator(&full, &graph.model, |_| -1),
        )?;
        let original_terms =
            FourDSector::new(original_bare, Vec::new(), Vec::new()).physical_terms()?;
        ensure!(original_terms.len() == 1);
        let original_weight = graph.global_prefactor.num.clone();
        let original_factor = graph.overall_factor.clone();

        // Keep lambda as a physical typed mass until every enclosing Taylor
        // operation has run. Only its numerical parameter slot is set to one.
        let lambda = Atom::from(UFOSymbol::from("MH"));
        let mut gluons = Vec::new();
        for edge in graph.paired_edges(&full) {
            if graph[edge].particle_pdg_code(&graph.model) == Some(21) {
                ensure!(graph[edge].mass_atom(&graph.model).is_zero());
                let numerator = graph[edge].num.clone();
                let target = &mut graph.underlying[edge];
                target.particle = target.particle.clone().override_mass(Some(lambda.clone()));
                target.mass =
                    EdgeMass::from_atom(lambda.clone(), &graph.model, &graph.param_builder)?;
                ensure!(target.num == numerator);
                gluons.push(edge);
            }
        }
        ensure!(gluons.len() == loops);
        ensure!(original_numerator.same_value(&graph.numerator(&full, &empty)?.state.expr));
        let regulated_numerator = graph.numerator(&full, &empty)?.state.expr;
        let regulated_bare = regulated_numerator.with_checked_expression(
            regulated_numerator.expression() * graph.denominator(&full, &graph.model, |_| -1),
        )?;
        let regulated_terms =
            FourDSector::new(regulated_bare, Vec::new(), Vec::new()).physical_terms()?;
        ensure!(regulated_terms.len() == 1);
        ensure!(
            original_terms[0]
                .numerator
                .same_value(&regulated_terms[0].numerator)
        );
        ensure!(original_terms[0].denominators.len() == regulated_terms[0].denominators.len());
        let mut typed_denominators = Vec::new();
        for (old, regulated) in original_terms[0]
            .denominators
            .iter()
            .zip(&regulated_terms[0].denominators)
        {
            let shifted = gluons.contains(&old.source_edge);
            let shift = if shifted { lambda.pow(2) } else { Atom::Zero };
            ensure!(old.source_edge == regulated.source_edge && old.momentum == regulated.momentum);
            ensure!(regulated.mass_squared == &old.mass_squared + &shift);
            ensure!(
                (&regulated.full_expr - &old.full_expr + &shift)
                    .expand()
                    .is_zero()
            );
            ensure!(
                regulated.full_expr.replace(lambda.to_pattern()).with(0) == old.full_expr,
                "typed gluon-mass zero specialization changed the original prescribed polynomial"
            );
            ensure!(
                regulated.mass_squared.replace(lambda.to_pattern()).with(0) == old.mass_squared
            );
            typed_denominators.push(serde_json::json!({"edge":old.source_edge.0,
                "momentum":old.momentum.to_canonical_string(),"mass_squared":regulated.mass_squared.to_canonical_string(),
                "polynomial":regulated.full_expr.to_canonical_string(),"shifted":shifted,"power":1}));
        }
        let native_vector_source = projection.numeric_input(&diagram, &reference.model)?;
        let mut regulated_vector_source = native_vector_source.clone();
        let poles = regulated_vector_source["propagators"]
            .as_array_mut()
            .ok_or_else(|| eyre!("native vector propagator export is not an array"))?;
        ensure!(poles.len() == typed_denominators.len());
        for pole in poles {
            let edge = usize::try_from(
                pole["edge"]
                    .as_u64()
                    .ok_or_else(|| eyre!("native vector edge is not an integer"))?,
            )?;
            ensure!(pole["power"] == 1);
            if gluons.iter().any(|gluon| gluon.0 == edge) {
                ensure!(pole["mass_squared"] == "0");
                let polynomial = pole["polynomial"]
                    .as_str()
                    .ok_or_else(|| eyre!("native polynomial is not a string"))?;
                pole["polynomial"] = serde_json::json!(format!("({polynomial})-1"));
                pole["mass_squared"] = serde_json::json!("1");
            }
        }

        // Build the finite native closing word from the source's own external
        // spin/color metric. It is independent of all counterterm coefficients.
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let (mass_dual, external_p) = reference.mass_projector(&diagram, &Atom::var(GS.dim))?;
        let mass_dual = mass_dual
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let mass_dual = mass_dual
            .with_checked_expression(mass_dual.expression().replace(CS.nc).with(3))?
            .try_into_indexed()?;
        let identity = mass_dual.with_checked_expression(
            Atom::num(12)
                * mass_dual.expression().replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                    {
                        **output = Atom::Zero;
                    }
                }),
        )?;
        let spin_rep = LibraryRep::from(Bispinor {}).symbol();
        let mut replaced = 0;
        let temporal =
            cooking
                .uncook(identity.expression().as_view())
                .replace_map(|view, _, output| {
                    if let AtomView::Fun(metric) = view
                        && metric.get_symbol() == ETS.metric
                        && metric.get_nargs() == 2
                        && (0..2).all(|index| {
                            matches!(metric.get(index), AtomView::Fun(rep)
                        if rep.get_symbol() == spin_rep && rep.get_nargs() == 2
                            && rep.get(0) == Atom::num(4).as_view())
                        })
                    {
                        replaced += 1;
                        **output = function!(AGS.gamma0, metric.get(0), metric.get(1));
                    }
                });
        ensure!(replaced == 1);
        let temporal = AppliedFeynmanRule::new(temporal)?.expr;
        let temporal = identity.with_checked_expression(temporal.expression().clone())?;
        let duals = [identity.clone(), temporal.clone()]
            .into_iter()
            .map(|word| Ok(word.with_checked_expression(word.expression() / Atom::num(12))?))
            .collect::<Result<Vec<_>>>()?;
        let open_basis = [identity, temporal]
            .iter()
            .map(|word| QcdMassReference::open_native_quark_dual(word, &original_numerator))
            .collect::<Result<Vec<_>>>()?;
        for (row, dual) in duals.iter().enumerate() {
            for (column, basis) in open_basis.iter().enumerate() {
                let trace = dual.clone().multiply(basis.clone())?;
                ensure!(trace.is_scalar());
                // Execute only this constant 4x4/3x3 closing word. No physical
                // numerator, Taylor coefficient or denominator enters the Gram.
                let source = cooking.uncook(trace.expression().as_view());
                let mut network = source.parse_into_net()?;
                network.execute::<Sequential, SmallestDegree, _, _, _>(
                    &*crate::utils::TENSORLIB.read().unwrap(),
                    &*crate::utils::FUN_LIB,
                )?;
                let scalar = match network.result_scalar()? {
                    ExecutionResult::One => Atom::one(),
                    ExecutionResult::Zero => Atom::Zero,
                    ExecutionResult::Val(value) => value.into_owned(),
                };
                ensure!(
                    (scalar - Atom::num(i32::from(row == column)))
                        .collect_compact_factors()
                        .is_zero(),
                    "emitted native temporal/identity Gram matrix entry ({row},{column})"
                );
            }
        }
        // This global factor is attached to each COMPLETE numerator product by
        // FinalIntegrandBuilder, after its child and parent operators have run.
        graph.global_prefactor.projector = duals[1].clone();
        ensure!(
            graph.global_prefactor.num == original_weight
                && graph.overall_factor == original_factor
        );
        ensure!(
            graph
                .paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut))
                .is_empty()
        );

        let model = graph.model.clone();
        let settings = GenerationSettings {
            explicit_orientation_sum_only: true,
            uv: UVgenerationSettings {
                generate_integrated: true,
                add_marker: false,
                local_uv_cts_from_expanded_4d_integrands: true,
                final_integrand: FinalIntegrandDimension::ThreeD,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    log_divergent: ApproximationType::MUV,
                    massive_power_divergent: ApproximationType::OS,
                    massless_power_divergent: ApproximationType::IR,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let numerator = graph.numerator(&full, &empty)?.state.expr;
        let options = graph.denominator_only_cff_3d_expression_options();
        let surfaces = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let production = graph.generate_3d_expression_for_integrand(
            &[],
            &surfaces,
            &options,
            Some(numerator.expression()),
        )?;
        let orientation = OrientationProjection::exact_expression(
            &production,
            &options,
            &settings.orientation_pattern,
            true,
        );
        let cuts = CutStructure::empty(&graph);
        validate_local_counterterm_schemes(&graph, &cuts, &settings.uv)?;
        let mut forests = Wood::new(cuts, &graph, &model, &settings)?.unfold();
        ensure!(forests.graph.n_nodes() == if loops == 1 { 2 } else { 4 });
        let mut histories = BTreeSet::new();
        for (node, _, _) in forests.graph.iter_nodes() {
            let mut history = Vec::new();
            let mut cursor = node;
            while !forests.graph[cursor].key.is_empty() {
                let covered = forests.graph[cursor]
                    .covers()
                    .ok_or_else(|| eyre!("nonempty native forest history lost its cover"))?;
                history.push(
                    graph
                        .paired_edges(&covered)
                        .into_iter()
                        .map(|edge| edge.0)
                        .collect::<Vec<_>>(),
                );
                cursor = forests
                    .graph
                    .unique_parent(cursor)
                    .ok_or_else(|| eyre!("connected native rainbow history lost its parent"))?
                    .0;
            }
            history.reverse();
            ensure!(
                histories.insert(history),
                "duplicate complete native rainbow history"
            );
        }
        let parent_edges = graph
            .paired_edges(&full)
            .into_iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>();
        let expected = if loops == 1 {
            BTreeSet::from([vec![], vec![parent_edges]])
        } else {
            ensure!(parent_edges == [2, 3, 4, 5, 6]);
            BTreeSet::from([
                vec![],
                vec![vec![2, 3]],
                vec![parent_edges.clone()],
                vec![vec![2, 3], parent_edges],
            ])
        };
        ensure!(
            histories == expected,
            "actual native wood omitted or repeated a required complete history"
        );
        let mut projection_context = Local4dProjectionContext {
            certify_denominator_reconstruction: true,
            ..Default::default()
        };
        forests.compute(
            &mut graph,
            &model,
            crate::utils::vakint()?,
            orientation,
            &settings.uv,
            &mut projection_context,
        )?;
        ensure!(projection_context.denominator_reconstruction_certificates > 0);
        ensure!(
            projection_context.denominator_reconstruction_certificates
                == projection_context.source_preparation_builds
        );
        let mut emitted = forests.orientation_parametric_exprs(&graph, &settings.uv)?;
        ensure!(emitted.len() == 1);
        let emitted = emitted.remove(0).integrands.resolved_scalars()?;
        let physical = emitted
            .fallible_map(|value| Ok(value.replace(GS.integrated_loop_scale).with(1)))?
            .map_numerators(|tensor| {
                Ok(tensor.with_checked_expression(
                    tensor
                        .expression()
                        .replace(GS.integrated_loop_scale)
                        .with(1),
                )?)
            })?;
        let roots = [physical.iter().map(|(_, value)| value).sum::<Atom>()];
        let builder = ParamBuilder::<f64>::new(
            &graph,
            &model,
            &graph.loop_momentum_basis,
            Vec::<Atom>::new(),
        )?;
        let evaluator_settings = crate::processes::EvaluatorSettings {
            inline_numerator_functions: false,
            do_algebra: false,
            compile: false,
            store_atom: true,
            ..Default::default()
        };
        let (stack, timings) = EvaluatorStack::new_explicit_sum_with_timings(
            &roots,
            &builder,
            physical.numerators(),
            None,
            &evaluator_settings,
        )?;
        let params = (&builder.pairs)
            .into_iter()
            .flat_map(|pair| pair.params.iter())
            .collect::<Vec<_>>();
        let mut values = builder.values[0].clone();
        let mut coordinates = Vec::new();
        for edge in &graph.loop_momentum_basis.loop_edges {
            let mut indices = [0; 3];
            for component in 1..=3 {
                let coordinate = GS.emr_mom(*edge, GS.cind(component));
                indices[component - 1] = params
                    .iter()
                    .position(|atom| **atom == coordinate)
                    .ok_or_else(|| {
                        eyre!("native integration loop parameter missing: {coordinate}")
                    })?;
            }
            coordinates.push(indices);
        }
        let routed_p = diagram.loop_momentum_basis().route_expression(&external_p);
        ensure!(
            params
                .iter()
                .filter(|parameter| ***parameter == lambda)
                .count()
                == 1,
            "the typed gluon mass needs exactly one numerical parameter slot"
        );
        for (index, parameter) in params.iter().enumerate() {
            if [
                GS.m_uv_vacuum,
                GS.m_uv_expansion,
                GS.renormalization_localization_scale,
                GS.mu_r_sq,
                GS.numerator_sampling_scale,
            ]
            .iter()
            .any(|head| **parameter == Atom::var(*head))
                || **parameter == lambda
            {
                values[index] = Complex::new_re(F(1.0));
            }
            if builder.pairs.external_energies.params.contains(parameter) {
                let AtomView::Fun(momentum) = parameter.as_view() else {
                    return Err(eyre!("native external energy is not a momentum component"));
                };
                let routed = diagram
                    .loop_momentum_basis()
                    .route_expression(&function!(GS.emr_mom, momentum.get(0)));
                let sign = if (&routed - &routed_p).expand().is_zero() {
                    1.0
                } else if (&routed + &routed_p).expand().is_zero() {
                    -1.0
                } else {
                    return Err(eyre!(
                        "native external energy has unmatched physical momentum"
                    ));
                };
                values[index] = Complex::new_re(F(sign));
            }
        }
        let mt_position = params
            .iter()
            .position(|parameter| **parameter == Atom::from(UFOSymbol::from("MT")))
            .ok_or_else(|| eyre!("native physical MT parameter missing"))?;
        ensure!(
            values[mt_position] == Complex::new_re(F(1.0)),
            "the independent unit-mass source was changed"
        );
        let parameter_index = |parameter: &Atom| {
            params
                .iter()
                .position(|candidate| **candidate == *parameter)
                .ok_or_else(|| {
                    eyre!("required numerical integration parameter is missing: {parameter}")
                })
        };
        let uv_mass_slots = [
            parameter_index(&Atom::var(GS.m_uv_vacuum))?,
            parameter_index(&Atom::var(GS.m_uv_expansion))?,
        ];
        let localization_slot = parameter_index(&Atom::var(GS.renormalization_localization_scale))?;
        let physical_slots = [
            mt_position,
            parameter_index(&lambda)?,
            parameter_index(&Atom::var(GS.mu_r_sq))?,
        ];
        // a=G²/(16 pi²), SU(3) CF=4/3. Preserve all source coupling phases.
        let g_position = params
            .iter()
            .position(|parameter| **parameter == Atom::from(UFOSymbol::from("G")))
            .ok_or_else(|| eyre!("native QCD coupling G parameter missing"))?;
        let g = values[g_position].re.0;
        ensure!(g > 0.0 && values[g_position].im.0 == 0.0);
        let a = g * g / (16.0 * std::f64::consts::PI.powi(2));
        let normalization = -1.0 / (a * (4.0 / 3.0)).powi(i32::try_from(loops)?);
        let source = serde_json::json!({
            "source":diagram.name(), "loops":loops, "typed_gluon_edges":gluons.iter().map(|edge| edge.0).collect::<Vec<_>>(),
            "mass":1, "gluon_mass":1, "uv_mass":1, "localization_scale":1, "mu_r_squared":1,
            "forest_owner":"production HedgePoset", "global_source_ordering":-1,
            "temporal_gram":[[1,0],[0,1]], "tree_factor":1,
            "root_bytes":roots[0].as_view().get_byte_size(), "numerator_definitions":physical.numerators().len(),
            "numerator_bytes":physical.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
            "construction_seconds":started.elapsed().as_secs_f64(), "evaluator_timings":format!("{timings:?}"),
            "coupling_g":g, "normalization_i_times":normalization,
            "native_vector_source":native_vector_source.clone(),
            "unregulated_source":{
                "projection":"pslash/(4*Nc*p_squared)",
                "source_kernel_weight":projection.kernel_weight.to_bare_ordered_string(),
                "external_ordering":projection.external_ordering.to_bare_ordered_string(),
                "numerator_prefactor":projection.numerator_prefactor.to_bare_ordered_string(),
                "propagator_prescription":"Feynman+i0",
                "source":native_vector_source.clone(),
            },
            "regulated_source":{
                "projection":"pslash/(4*Nc*p_squared)",
                "source_kernel_weight":projection.kernel_weight.to_bare_ordered_string(),
                "external_ordering":projection.external_ordering.to_bare_ordered_string(),
                "numerator_prefactor":projection.numerator_prefactor.to_bare_ordered_string(),
                "propagator_prescription":"Feynman+i0",
                "source":regulated_vector_source,
            },
            "typed_denominator_inventory":typed_denominators,
            "typed_denominator_shift_identity":true,"native_full_D_numerator_unchanged":true,
            "typed_zero_mass_denominator_specialization":true,
            "massless_specialization_scope":"same native quark OS histories; O and degree-one Taylor before numerical MH=0",
            "vacuum_invariance_reference":"native_d09_positive_gluon_mass_preserves_all_three_vacuum_inputs",
            "typed_gluon_mass_parameter_slots":1,"physical_MT_parameter_value":1,
            "complete_forest_histories":histories,
            "denominator_reconstruction_certificates":projection_context.denominator_reconstruction_certificates,
            "source_preparation_builds":projection_context.source_preparation_builds,
            "extra_R2_term_added":false,
        });
        Ok(Self {
            evaluator: stack.single_parametric,
            values,
            coordinates,
            normalization,
            source,
            evaluations: 0,
            precise_evaluations: 0,
            precise_seconds: 0.0,
            uv_mass_slots,
            localization_slot,
            physical_slots,
            gluon_mass: 1.0,
        })
    }

    fn set_auxiliary_scales(&mut self, uv_mass: f64, localization_scale: f64) -> Result<()> {
        ensure!(uv_mass > 0.0 && localization_scale > 0.0);
        for position in self.uv_mass_slots {
            self.values[position] = Complex::new_re(F(uv_mass));
        }
        self.values[self.localization_slot] = Complex::new_re(F(localization_scale));
        ensure!(
            self.physical_slots
                .iter()
                .zip([1.0, self.gluon_mass, 1.0])
                .all(|(position, value)| self.values[*position] == Complex::new_re(F(value))),
            "an auxiliary-scale variation changed the physical mass, gluon mass or renormalization scale"
        );
        Ok(())
    }

    fn set_gluon_mass(&mut self, mass: f64) -> Result<()> {
        ensure!(mass == 0.0 || mass == 1.0);
        ensure!(
            self.values[self.physical_slots[0]] == Complex::new_re(F(1.0))
                && self.values[self.physical_slots[2]] == Complex::new_re(F(1.0))
        );
        self.gluon_mass = mass;
        self.values[self.physical_slots[1]] = Complex::new_re(F(mass));
        Ok(())
    }

    fn evaluate(&mut self, momenta: &[[f64; 3]]) -> Result<Complex<F<f64>>> {
        ensure!(momenta.len() == self.coordinates.len());
        for (indices, momentum) in self.coordinates.iter().zip(momenta) {
            for (index, value) in indices.iter().zip(momentum) {
                self.values[*index] = Complex::new_re(F(*value));
            }
        }
        self.evaluations += 1;
        // Hard-region cancellation is performed on the COMPLETE forest in the
        // same evaluator. Keep the original finite f64 sample as the input
        // boundary, and increase arithmetic precision rather than cutting tails.
        let hard = momenta
            .iter()
            .flatten()
            .any(|coordinate| coordinate.abs() > 1.0e3);
        let norm = |momentum: &[f64; 3]| {
            momentum
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt()
        };
        let soft = self.gluon_mass == 0.0 && {
            let mut minimum = momenta.iter().map(norm).fold(f64::INFINITY, f64::min);
            if momenta.len() == 2 {
                minimum = minimum.min(norm(&std::array::from_fn(|component| {
                    momenta[0][component] - momenta[1][component]
                })));
            }
            minimum < 1.0e-3
        };
        let precise = hard || soft;
        let (real, imaginary) = if precise {
            let started = std::time::Instant::now();
            let values = self
                .values
                .iter()
                .map(|value| {
                    Complex::new(
                        F::<ArbPrec>::from_ff64(value.re),
                        F::<ArbPrec>::from_ff64(value.im),
                    )
                })
                .collect::<Vec<_>>();
            let value = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(
                &mut self.evaluator,
            )(&values);
            self.precise_evaluations += 1;
            self.precise_seconds += started.elapsed().as_secs_f64();
            (value.re.into_f64(), value.im.into_f64())
        } else {
            let value = <f64 as GenericEvaluatorFloat>::get_evaluator_single(&mut self.evaluator)(
                &self.values,
            );
            (value.re.0, value.im.0)
        };
        ensure!(
            real.is_finite() && imaginary.is_finite(),
            "nonfinite complete emitted forest at {momenta:?}"
        );
        Ok(Complex::new(
            F(-self.normalization * imaginary),
            F(self.normalization * real),
        ))
    }

    /// Two linear physical charts, each split into its two Hepp sectors.
    /// The channel partition is evaluated on the native l,k point, and the
    /// complete emitted forest is evaluated once for each weighted image.
    fn soft_chart_samples(xs: &[f64]) -> Result<Vec<([[f64; 3]; 2], f64)>> {
        ensure!(xs.len() == 3 && xs.iter().all(|x| *x > 0.0 && *x < 1.0));
        let radius = xs[0] / (1.0 - xs[0]);
        let ratio = xs[1];
        let z = 2.0 * xs[2] - 1.0;
        let sine = (1.0 - z * z).sqrt();
        let mut samples = Vec::with_capacity(4);
        for channel in 0..2 {
            for sector in 0..2 {
                let (h, k) = if sector == 0 {
                    (radius, radius * ratio)
                } else {
                    (radius * ratio, radius)
                };
                let hard = [0.0, 0.0, k];
                let radial = [h * sine, 0.0, h * z];
                let l = if channel == 0 {
                    radial
                } else {
                    std::array::from_fn(|component| radial[component] + hard[component])
                };
                let q: [f64; 3] = std::array::from_fn(|component| l[component] - hard[component]);
                let l2 = l.iter().map(|value| value * value).sum::<f64>();
                let q2 = q.iter().map(|value| value * value).sum::<f64>();
                let l_scale = l2.sqrt().max(k);
                let q_scale = q2.sqrt().max(k);
                let d_l = l2 * l_scale * (1.0 + l_scale).powi(2);
                let d_q = q2 * q_scale * (1.0 + q_scale).powi(2);
                let partition = (if channel == 0 { d_q } else { d_l }) / (d_l + d_q);
                let jacobian = 16.0 * std::f64::consts::PI.powi(2) * h * h * k * k * radius
                    / (1.0 - xs[0]).powi(2);
                ensure!(partition.is_finite() && jacobian.is_finite());
                samples.push(([l, hard], partition * jacobian));
            }
        }
        Ok(samples)
    }

    fn integrate(
        &mut self,
        seed: u64,
        samples_per_iteration: usize,
        iterations: usize,
    ) -> Result<serde_json::Value> {
        let started = std::time::Instant::now();
        let initial_evaluations = self.evaluations;
        let initial_precise = self.precise_evaluations;
        let initial_precise_seconds = self.precise_seconds;
        let dimensions = if self.coordinates.len() == 1 { 1 } else { 3 };
        let mut grid = ContinuousGrid::<f64>::new(dimensions, 64, 100, None, false)
            .map_err(|error| eyre!(error))?;
        let mut rng = MonteCarloRng::new(seed, 0);
        let mut real = StatisticsAccumulator::<f64>::new();
        let mut imaginary = StatisticsAccumulator::<f64>::new();
        let mut history = Vec::new();
        for iteration in 0..iterations {
            for _ in 0..samples_per_iteration {
                let mut sample = Sample::new();
                grid.sample(&mut rng, &mut sample);
                let Sample::Continuous(weight, xs) = &sample else {
                    return Err(eyre!("continuous grid emitted noncontinuous sample"));
                };
                let (re, im) = if dimensions == 3 && self.gluon_mass == 0.0 {
                    let mut real = 0.0;
                    let mut imaginary = 0.0;
                    // These are four explicit integration images, not a
                    // randomly selected channel. There is no factor of 1/4.
                    for (momenta, weight) in Self::soft_chart_samples(xs)? {
                        let value = self.evaluate(&momenta)?;
                        real += weight * value.re.0;
                        imaginary += weight * value.im.0;
                    }
                    (real, imaginary)
                } else {
                    let r = xs[0] / (1.0 - xs[0]);
                    let jacobian_r = r * r / (1.0 - xs[0]).powi(2);
                    let (momenta, jacobian) = if dimensions == 1 {
                        (vec![[0.0, 0.0, r]], 4.0 * std::f64::consts::PI * jacobian_r)
                    } else {
                        let s = xs[1] / (1.0 - xs[1]);
                        let z = 2.0 * xs[2] - 1.0;
                        (
                            vec![[0.0, 0.0, r], [s * (1.0 - z * z).sqrt(), 0.0, s * z]],
                            16.0 * std::f64::consts::PI.powi(2) * jacobian_r * s * s
                                / (1.0 - xs[1]).powi(2),
                        )
                    };
                    let value = self.evaluate(&momenta)?;
                    (value.re.0 * jacobian, value.im.0 * jacobian)
                };
                grid.add_training_sample(&sample, re.abs())
                    .map_err(|error| eyre!(error))?;
                real.add_sample(re * weight, Some(&sample));
                imaginary.add_sample(im * weight, Some(&sample));
            }
            real.update_iter(false);
            imaginary.update_iter(false);
            grid.update(0.5);
            history.push(
                serde_json::json!({"iteration":iteration,"value":real.avg,"error":real.err,
                "imaginary":imaginary.avg,"imaginary_error":imaginary.err}),
            );
            crate::debug_tags!(#uv,#integrated,#summary;
                stage="native_emitted_vector_integration_iteration", seed, iteration,
                value=real.avg, error=real.err, imaginary=imaginary.avg,
                "Integrated the complete production-emitted native temporal coefficient"
            );
        }
        Ok(
            serde_json::json!({"seed":seed,"stream":0,"samples_per_iteration":samples_per_iteration,
            "iterations":iterations,"samples":real.processed_samples,"value":real.avg,"error":real.err,
            "imaginary":imaginary.avg,"imaginary_error":imaginary.err,"chi_squared":real.chi_sq,
            "seconds":started.elapsed().as_secs_f64(),"history":history,"dimension":dimensions,
            "evaluations":self.evaluations-initial_evaluations,
            "precise_evaluations":self.precise_evaluations-initial_precise,
            "precise_seconds":self.precise_seconds-initial_precise_seconds,
            "uv_mass":self.values[self.uv_mass_slots[0]].re.0,
            "localization_scale":self.values[self.localization_slot].re.0,
            "gluon_mass":self.gluon_mass,
            "integration_map":if self.gluon_mass == 0.0 && dimensions == 3 {"two soft charts times two Hepp sectors, explicit sum"} else {"product radial"},
            "precision_policy":if self.gluon_mass == 0.0 {"ArbPrec for max absolute loop component > 1000 or minimum native soft norm < 0.001; f64 otherwise"}
                else {"ArbPrec for max absolute loop component > 1000; f64 otherwise"}}),
        )
    }

    fn one_loop_calibration(
        reference: &QcdMassReference,
        gluon_mass: f64,
    ) -> Result<serde_json::Value> {
        let mut calibration = Self::new(reference, 1)?;
        calibration.set_gluon_mass(gluon_mass)?;
        let seed = if gluon_mass == 0.0 {
            2026102110
        } else {
            2026102101
        };
        let measured = calibration.integrate(seed, 20_000, 6)?;
        let expected = if gluon_mass == 0.0 {
            -2.0
        } else {
            -1.0 + std::f64::consts::PI / 3.0_f64.sqrt()
        };
        let value = measured["value"].as_f64().unwrap();
        let error = measured["error"].as_f64().unwrap();
        ensure!(
            (value - expected).abs() < 5.0 * error + 1.0e-4,
            "one-loop emitted vector normalization: {value} +/- {error}, expected {expected}"
        );
        Ok(
            serde_json::json!({"source":calibration.source,"measured":measured,"reference":expected}),
        )
    }
}

#[test]
fn native_emitted_one_loop_vector_matches_massless_and_positive_mass_references() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    for gluon_mass in [1.0, 0.0] {
        let measured = EmittedVectorIntegral::one_loop_calibration(&reference, gluon_mass)?;
        crate::debug_tags!(#uv,#integrated,#summary;
            stage="native_emitted_one_loop_vector_calibration", file.measured=%serde_json::to_string(&measured)?,
            "Production-emitted one-loop field coefficient matches its independent reference"
        );
    }
    Ok(())
}

pub(super) fn two_loop(reference: &QcdMassReference) -> Result<serde_json::Value> {
    let one_loop = EmittedVectorIntegral::one_loop_calibration(reference, 1.0)?;
    let mut integral = EmittedVectorIntegral::new(reference, 2)?;
    let momenta = [[0.0, 0.0, 0.7], [0.8, 0.0, -0.4]];
    let baseline = integral.evaluate(&momenta)?;
    let rotated = integral.evaluate(&[[0.7, 0.0, 0.0], [-0.4, 0.8, 0.0]])?;
    ensure!(
        (baseline.re.0 - rotated.re.0).abs() <= 1.0e-9 * (1.0 + baseline.re.0.abs()),
        "rest-frame rotation changed the complete scalar"
    );
    let precise_values = integral
        .values
        .iter()
        .map(|value| {
            Complex::new(
                F::<ArbPrec>::from_ff64(value.re),
                F::<ArbPrec>::from_ff64(value.im),
            )
        })
        .collect::<Vec<_>>();
    let precise = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(&mut integral.evaluator)(
        &precise_values,
    );
    let precise_re = -integral.normalization * precise.im.into_f64();
    ensure!((rotated.re.0 - precise_re).abs() < 1.0e-9 * (1.0 + precise_re.abs()));
    let mut tail_checks = Vec::new();
    for radius in [1.0e2, 1.0e3, 1.0e5] {
        for hierarchy in [false, true] {
            let r = radius;
            let s = if hierarchy { 0.7 } else { radius * 0.8 };
            let point = [[0.0, 0.0, r], [s * 0.8, 0.0, s * 0.6]];
            let actual = integral.evaluate(&point)?;
            let binary = <f64 as GenericEvaluatorFloat>::get_evaluator_single(
                &mut integral.evaluator,
            )(&integral.values);
            let values = integral
                .values
                .iter()
                .map(|value| {
                    Complex::new(
                        F::<ArbPrec>::from_ff64(value.re),
                        F::<ArbPrec>::from_ff64(value.im),
                    )
                })
                .collect::<Vec<_>>();
            let precise = <ArbPrec as GenericEvaluatorFloat>::get_evaluator_single(
                &mut integral.evaluator,
            )(&values);
            let reference = -integral.normalization * precise.im.into_f64();
            let binary = -integral.normalization * binary.im.0;
            let relative = (actual.re.0 - reference).abs() / reference.abs().max(1.0e-300);
            if radius <= 1.0e3 {
                ensure!(
                    relative < 1.0e-5,
                    "binary complete-forest tail precision below switch: r={r}, s={s}, relative={relative}"
                );
            } else {
                ensure!(
                    actual.re.0 == reference,
                    "hard-tail policy did not use the same retained arbitrary-precision evaluator"
                );
            }
            tail_checks.push(serde_json::json!({"r":r,"s":s,"hierarchical":hierarchy,
                "policy_value":actual.re.0,"f64_value":binary,"precise_value":reference,"policy_relative_error":relative}));
        }
    }
    let first = integral.integrate(2026102102, 100_000, 10)?;
    let second = integral.integrate(2026102103, 100_000, 10)?;
    integral.set_auxiliary_scales(2.0, 1.0)?;
    let changed_uv_mass = integral.integrate(2026102104, 100_000, 10)?;
    integral.set_auxiliary_scales(1.0, 2.0)?;
    let changed_localization = integral.integrate(2026102105, 100_000, 10)?;
    integral.set_auxiliary_scales(1.0, 1.0)?;
    let massless_one_loop = EmittedVectorIntegral::one_loop_calibration(reference, 0.0)?;
    integral.set_gluon_mass(0.0)?;
    let d_l = Atom::var(symbol!("native_emitted_soft_partition::D_l"; Scalar));
    let d_q = Atom::var(symbol!("native_emitted_soft_partition::D_q"; Scalar));
    ensure!(
        (&d_q / (&d_l + &d_q) + &d_l / (&d_l + &d_q) - Atom::one())
            .together()
            .is_zero(),
        "the two native soft charts do not form an exact partition of unity"
    );
    let mut chart_checks = Vec::new();
    for xs in [[0.37, 0.41, 0.63], [0.012, 0.019, 0.92], [0.83, 0.94, 0.04]] {
        let samples = EmittedVectorIntegral::soft_chart_samples(&xs)?;
        ensure!(samples.len() == 4);
        for (image, (momenta, weighted_jacobian)) in samples.into_iter().enumerate() {
            let radius = xs[0] / (1.0 - xs[0]);
            let ratio = xs[1];
            let norm =
                |vector: &[f64; 3]| vector.iter().map(|value| value * value).sum::<f64>().sqrt();
            let [l, k] = momenta;
            let q: [f64; 3] = std::array::from_fn(|component| l[component] - k[component]);
            let h = if image < 2 { l } else { q };
            let expected_h = if image % 2 == 0 {
                radius
            } else {
                radius * ratio
            };
            let expected_k = if image % 2 == 0 {
                radius * ratio
            } else {
                radius
            };
            ensure!((norm(&h) - expected_h).abs() < 1.0e-14 * (1.0 + radius));
            ensure!((norm(&k) - expected_k).abs() < 1.0e-14 * (1.0 + radius));
            let l_scale = norm(&l).max(norm(&k));
            let q_scale = norm(&q).max(norm(&k));
            let d_l = norm(&l).powi(2) * l_scale * (1.0 + l_scale).powi(2);
            let d_q = norm(&q).powi(2) * q_scale * (1.0 + q_scale).powi(2);
            ensure!((d_q / (d_l + d_q) + d_l / (d_l + d_q) - 1.0).abs() < 1.0e-15);
            let value = integral.evaluate(&momenta)?;
            let rotated = momenta.map(|vector| [vector[2], vector[0], vector[1]]);
            let orbit = integral.evaluate(&rotated)?;
            ensure!((value.re.0 - orbit.re.0).abs() <= 1.0e-9 * (1.0 + value.re.0.abs()));
            chart_checks.push(
                serde_json::json!({"xs":xs,"image":image,"native_l":l,"native_k":k,
                "native_q":q,"weighted_jacobian":weighted_jacobian,"value":value.re.0,
                "orbit_value":orbit.re.0,"partition_sum":d_q/(d_l+d_q)+d_l/(d_l+d_q)}),
            );
        }
    }
    let unregulated_first = integral.integrate(2026102111, 100_000, 20)?;
    let unregulated_second = integral.integrate(2026102112, 100_000, 20)?;
    integral.set_auxiliary_scales(2.0, 1.0)?;
    let unregulated_uv_mass = integral.integrate(2026102113, 100_000, 20)?;
    integral.set_auxiliary_scales(1.0, 2.0)?;
    let unregulated_localization = integral.integrate(2026102114, 100_000, 20)?;
    integral.set_auxiliary_scales(1.0, 1.0)?;
    Ok(
        serde_json::json!({"source":integral.source,"one_loop_calibration":one_loop,
        "pointwise_rotation_and_precision_checks":true,"tail_precision_checks":tail_checks,
        "integrations":[first,second,changed_uv_mass,changed_localization],
        "unregulated_integrations":[unregulated_first,unregulated_second,unregulated_uv_mass,unregulated_localization],
        "unregulated_one_loop_calibration":massless_one_loop,
        "unregulated_chart_checks":chart_checks,"exact_soft_chart_partition":true,
        "measure":"plain d3k d3l, rest-frame 8pi2 radial reduction",
        "threshold_regulator":"both typed physical gluon masses=lambda=1",
        "absolute_unregulated_integrability_claimed":false}),
    )
}
