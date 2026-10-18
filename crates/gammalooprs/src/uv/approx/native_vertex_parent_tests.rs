//! Private proof boundary for actual logarithmic children inside OS parents.
//! Public mixed-prescription admission is deliberately unchanged.

use super::*;
use crate::{
    initialisation::test_initialise,
    numerator::AppliedFeynmanRule,
    uv::{
        approx::{integrated, integrated::native_mass_tests::QcdMassReference},
        hedge_poset::OwnedForestNode,
        orchestrator::validate_local_counterterm_schemes,
    },
};
use idenso::{
    color::{CS, ColorSimplifySettings},
    dirac::GammaSimplifySettings,
    representations::Bispinor,
};
use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use spenso::structure::{
    TensorStructure,
    abstract_index::{AIND_SYMBOLS, AbstractIndex},
    partial::PartialStructure,
    representation::{LibraryRep, Minkowski, RepName},
    slot::{IsAbstractSlot, Slot, SlotMatcher},
};

#[test]
fn native_logarithmic_vertex_prefixes_preserve_parent_mass_projection() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_parent_mass_projection(["1", "0", "0", "0"])
}

#[test]
fn native_logarithmic_vertex_prefixes_preserve_moving_parent_mass_projection() -> Result<()> {
    test_initialise()?;
    QcdMassReference::new()?.certify_parent_mass_projection(["5/4", "1/4", "1/2", "1/2"])
}

#[test]
fn native_parent_scalar_certificate_preserves_factored_cancellation() -> Result<()> {
    test_initialise()?;
    let dimension = Atom::var(GS.dim);
    let imaginary = Atom::i();
    // D01 coefficient row 57, with its five complete denominator functions
    // replaced by independent scalar aliases. Keep their inverses and the
    // original four complex factors; do not clear a common denominator.
    let inverses = (0..5).map(|index| {
        Atom::var(symbol!(
            &format!("native_parent_scalar_regression::s{index}"); Scalar
        ))
        .pow(-1)
    });
    let common = Atom::var(symbol!("UFO::GC_10")).pow(2)
        * Atom::var(symbol!("UFO::GC_11")).pow(2)
        * Atom::mul_many(inverses);
    let scalar = Atom::add_many([
        (Atom::num(-2) * &imaginary * &dimension + Atom::num(4) * &imaginary)
            * Atom::num(4)
            * &common,
        (-&imaginary * &dimension + Atom::num(2) * &imaginary) * Atom::num(8) * &common,
        (Atom::num(2) * &imaginary * &dimension - Atom::num(4) * &imaginary)
            * Atom::num(4)
            * &common,
        (&imaginary * &dimension - Atom::num(2) * &imaginary) * Atom::num(8) * &common,
    ]);
    let metric = function!(
        ETS.metric,
        idenso::bis!(
            4,
            Atom::var(symbol!("native_parent_scalar_regression::row"))
        ),
        idenso::bis!(
            4,
            Atom::var(symbol!("native_parent_scalar_regression::column"))
        )
    );
    QcdMassReference::certify_parent_tensor_zero(
        "captured scalar coefficient on an open metric",
        &AppliedFeynmanRule::new(&scalar * &metric)?.expr,
    )?;
    assert!(
        QcdMassReference::certify_parent_tensor_zero(
            "unit perturbation of the same scalar coefficient",
            &AppliedFeynmanRule::new((&scalar + Atom::one()) * &metric)?.expr,
        )
        .is_err(),
        "the private certificate must reject a nonzero metric contribution"
    );
    Ok(())
}

#[test]
fn native_parent_scalar_certificate_accepts_ranked_zero() -> Result<()> {
    test_initialise()?;
    let adjoint = idenso::representations::ColorAdjoint {}.new_rep(8);
    let [a, b, c, d] = ["a", "b", "c", "d"].map(|name| {
        adjoint.to_symbolic([Atom::var(symbol!(&format!(
            "native_parent_ranked_zero::{name}"
        )))])
    });
    let metric = function!(ETS.metric, &c, &d);
    let product = function!(CS.f, &a, &b, &c) * function!(CS.f, &a, &b, &d);
    // The final SU(3) color reduction gives a zero with two retained ports.
    // Its coefficient need not have a scalar interface to certify this zero.
    QcdMassReference::certify_parent_tensor_zero(
        "open adjoint Casimir identity",
        &AppliedFeynmanRule::new(&product - Atom::num(3) * &metric)?.expr,
    )?;
    assert!(
        QcdMassReference::certify_parent_tensor_zero(
            "nonzero adjoint metric remainder",
            &AppliedFeynmanRule::new(product - Atom::num(2) * metric)?.expr,
        )
        .is_err()
    );
    Ok(())
}

impl QcdMassReference {
    fn certify_parent_tensor_zero(
        label: &str,
        difference: &SymbolicTensor<IndexedStructure>,
    ) -> Result<()> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_vertex_parent_exact_input", label,
            file.difference = %cooking.uncook(difference.expression().as_view()).to_canonical_string(),
            "Retained factorized input to the parent proof certificate"
        );
        let reduced = integrated::simplify(difference)
            .wrap_err_with(|| format!("{label}: parent proof algebra"))?;
        let value = reduced.simplify_algebra(&idenso::tensor::AlgebraSettings {
            // The source has native SU(3) ports and the independent
            // vertex reference uses CF=4/3, CA=3, as in the child test.
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            ..Default::default()
        })?;
        // The coefficient engine allows partial words inside the complete tensor.
        // Move the reduced carrier into that boundary, retaining layout and proofs.
        let value: SymbolicTensor<PartialStructure> = value.into();
        // Native collection keeps each complete tensor word opaque. Normalize
        // only its scalar coefficient, as in the independent bare-source proof;
        // complete denominator functions and powers retain their owners.
        for (row, (word, coefficient)) in value
            .coefficient_list(spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
            .into_iter()
            .enumerate()
        {
            // A certified zero can retain logical tensor ports. Preserve the
            // native certificate's zero-factor rule before testing scalar rank.
            if word.expression().is_zero() || coefficient.expression().is_zero() {
                continue;
            }
            eyre::ensure!(
                coefficient.is_scalar(),
                "{label}: nonzero coefficient {row} retains tensor ports"
            );
            let normalized = coefficient.expression().collect_compact_factors();
            if !normalized.is_zero() {
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_parent_exact_coefficient", label, row,
                    file.coefficient = %coefficient.expression().to_canonical_string(),
                    file.normalized = %normalized.to_canonical_string(),
                    "Unresolved scalar coefficient retains its factorized source"
                );
            }
            eyre::ensure!(
                normalized.is_zero(),
                "{label}: exact scalar coefficient {row} is not zero"
            );
        }
        Ok(())
    }

    /// Reversible external-frame certificate only: no internal-d completeness.
    pub(in crate::uv) fn reduce_parent_external_frame(
        label: &str,
        difference: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let source_tensor = difference;
        let source_expression = cooking.uncook(source_tensor.expression().as_view());
        let difference = &source_expression;
        // Only the certified component OS operators introduce this fixed frame.
        // Embed its four Gram relations into formal d, without assuming that
        // these four vectors span the internal Lorentz space.
        let d = Atom::var(GS.dim);
        let rep = Minkowski {}.to_symbolic([&d]);
        let frame = symbol!(
            "native_vertex_parent::frame",
            tags = [SPENSO_TAG.rank1.clone(), SPENSO_TAG.tensor.clone()]
        );
        let names = (0..4).map(|i| function!(frame, i)).collect::<Vec<_>>();
        let vectors = (0..4)
            .map(|i| function!(frame, i, &rep))
            .collect::<Vec<_>>();
        let selectors = (0..4)
            .map(|i| function!(ETS.delta, function!(AIND_SYMBOLS.cind, i), &rep))
            .collect::<Vec<_>>();
        let mut kinematics = feynkit_kinematics::Kinematics::in_dimension(&d)?;
        for i in 0..4 {
            for j in 0..4 {
                let gram = Atom::num(if i != j {
                    0
                } else if i == 0 {
                    1
                } else {
                    -1
                });
                assert_eq!(
                    function!(SPENSO_TAG.dot, &selectors[i], &selectors[j]),
                    gram
                );
                kinematics = kinematics.with_scalar_product(&names[i], &names[j], gram)?;
            }
        }
        // Embed the raw class before any gamma-zero-sensitive reduction.
        // A literal inverse below certifies the whole factorized source;
        // no completeness relation for the internal d-space is assumed.
        let momentum_leaves = |value: &Atom| {
            let mut leaves = BTreeSet::new();
            value.visitor(&mut |view| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    assert!(momentum.get_nargs() >= 2);
                    let AtomView::Fun(port) = momentum.get(momentum.get_nargs() - 1) else {
                        panic!("the frame certificate does not convert fixed internal momentum components")
                    };
                    assert_eq!(port.get_symbol(), LibraryRep::from(Minkowski {}).symbol());
                    assert!(matches!(port.get_nargs(), 1 | 2));
                    assert_eq!(port.get(0), d.as_view());
                    leaves.insert(view.to_owned());
                }
                true
            });
            leaves
        };
        assert!(!difference.contains_symbol(frame));
        let momenta = momentum_leaves(difference);
        let embedded = difference.replace_map(|view, _, output| {
            if let AtomView::Fun(value) = view {
                if value.get_symbol() == AGS.gamma0 {
                    assert_eq!(value.get_nargs(), 2);
                    **output = function!(AGS.gamma, value.get(0), value.get(1), &vectors[0]);
                } else if value.get_symbol() == ETS.delta {
                    assert_eq!(value.get_nargs(), 2);
                    let AtomView::Fun(coordinate) = value.get(0) else {
                        panic!("the external frame requires a native fixed coordinate")
                    };
                    assert_eq!(coordinate.get_symbol(), AIND_SYMBOLS.cind);
                    assert_eq!(coordinate.get_nargs(), 1);
                    let i =
                        usize::try_from(coordinate.get(0)).expect("literal external coordinate");
                    assert!(i < 4);
                    let AtomView::Fun(port) = value.get(1) else {
                        panic!("the external frame requires an explicit Minkowski port")
                    };
                    assert_eq!(port.get_symbol(), LibraryRep::from(Minkowski {}).symbol());
                    assert!(matches!(port.get_nargs(), 1 | 2));
                    assert_eq!(port.get(0), d.as_view());
                    **output = function!(frame, i, value.get(1));
                }
            }
        });
        assert_eq!(momentum_leaves(&embedded), momenta);
        let reverse = |value: &Atom| {
            value.replace_map(|view, _, output| {
                if let AtomView::Fun(value) = view {
                    if value.get_symbol() == AGS.gamma
                        && value.get_nargs() == 3
                        && value.get(2) == vectors[0].as_view()
                    {
                        **output = function!(AGS.gamma0, value.get(0), value.get(1));
                    } else if value.get_symbol() == frame {
                        assert_eq!(value.get_nargs(), 2);
                        **output = function!(
                            ETS.delta,
                            function!(AIND_SYMBOLS.cind, value.get(0)),
                            value.get(1)
                        );
                    }
                }
            })
        };
        assert_eq!(
            reverse(&embedded),
            *difference,
            "literal source reconstruction precedes all formal-d Clifford algebra"
        );
        let embedded_tensor =
            source_tensor.with_checked_expression(cooking.try_cook(embedded.as_view())?)?;
        assert_eq!(
            source_tensor.structure().logical_slots(),
            embedded_tensor.structure().logical_slots()
        );
        // All public ports are explicit. Commutative rebuilding may change
        // their discovery layout, but never their named tensor interface.
        assert_eq!(
            source_tensor.structure().canonical().external_structure(),
            embedded_tensor.structure().canonical().external_structure()
        );
        let algebra = idenso::tensor::AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            epsilon: true,
            ..Default::default()
        };
        let ordinary_tensor = embedded_tensor.simplify_algebra(&algebra)?.to_dots()?;
        let ordinary = cooking.uncook(ordinary_tensor.expression().as_view());
        // Chain traversal may transpose the embedded temporal matrix. Only
        // E0 inherits the native gamma0 endpoint symmetry; no generic gamma
        // or internal momentum is made symmetric by this certificate.
        let mut temporal_transposes = 0;
        let ordered = ordinary.replace_map(|view, _, output| {
            if let AtomView::Fun(value) = view
                && value.get_symbol() == AGS.gamma
                && value.get_nargs() == 3
                && value.get(2) == vectors[0].as_view()
            {
                let native = function!(AGS.gamma0, value.get(0), value.get(1));
                let AtomView::Fun(temporal) = native.as_view() else {
                    unreachable!()
                };
                let normalized =
                    function!(AGS.gamma, temporal.get(0), temporal.get(1), &vectors[0]);
                if normalized.as_view() != view {
                    temporal_transposes += 1;
                    **output = normalized;
                }
            }
        });
        assert_eq!(reverse(&ordered), reverse(&ordinary));
        assert_eq!(momentum_leaves(&ordered), momentum_leaves(&ordinary));
        let completed_tensor = ordinary_tensor
            .with_checked_expression(cooking.try_cook(ordered.as_view())?)?
            .simplify_algebra(&algebra)?
            .to_dots()?;
        let completed = cooking.uncook(completed_tensor.expression().as_view());
        let constrained = kinematics.apply(&completed);
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_vertex_parent_external_frame_certificate", label,
            file.unembedded = %difference.to_canonical_string(),
            file.embedded = %embedded.to_canonical_string(),
            file.constrained = %constrained.to_canonical_string(),
            unchanged_momentum_leaves = momenta.len(), temporal_transposes,
            "Reversible external frame, native temporal symmetry and Gram relations retain internal d"
        );
        Ok(completed_tensor.with_checked_expression(cooking.try_cook(constrained.as_view())?)?)
    }

    /// Physical rational-finalizer boundary only; no energy residues or evaluator.
    pub(in crate::uv) fn certify_d09_rational_finalizer_mass_projection(
        &self,
        diagram: &feynkit_graph::FeynmanDiagram,
        ctx: &UVCtx<'_>,
        given: &OwnedForestNode,
        prefix: &Full4dCts,
        actual: &(OwnedForestNode, Local4dCts, IntegratedCts),
    ) -> Result<()> {
        use crate::uv::approx::final_integrand::FinalIntegrandBuilder;

        let (parent, local, integrated) = actual;
        let graph = ctx.graph;
        eyre::ensure!(diagram.name() == "D09" && graph.name == "D09");
        eyre::ensure!(parent.renormalization_scheme() == ApproximationType::OS);
        eyre::ensure!(
            ctx.settings.on_shell_velocity()? == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]
        );
        eyre::ensure!(local.recursive_completion().is_empty());
        let reduced = graph
            .full_filter()
            .subtract(parent.subgraph())
            .subtract(&graph.initial_state_cut);
        let normalize = FinalIntegrandBuilder::final_normalizer(graph, &reduced);
        let normalize_scalar = FinalIntegrandBuilder::scalar_normalizer(graph, &reduced);
        let scalar_tail = |atom: &Atom| {
            normalize_scalar(atom)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim_epsilon)
                .with(0)
        };
        let boundary = graph
            .dummy_less_full_crown(parent.subgraph())
            .included_iter()
            .collect::<Vec<_>>();
        eyre::ensure!(boundary.len() == 2);
        let incoming = boundary
            .iter()
            .find(|hedge| graph.flow(**hedge) == Flow::Sink)
            .ok_or_else(|| eyre!("D09 incoming parent port missing"))?;
        let owner = graph[incoming];
        let mass = graph[owner].particle.mass_atom(&graph.model);
        let (projector, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
        let basis = diagram.loop_momentum_basis();
        let routed_p = basis.route_expression(&reference_p);
        eyre::ensure!(
            (basis.route_expression(&function!(GS.emr_mom, owner.0)) - &routed_p)
                .expand()
                .is_zero()
        );
        let mut signs = BTreeMap::new();
        for hedge in &boundary {
            let edge = graph[hedge];
            let routed = basis.route_expression(&function!(GS.emr_mom, edge.0));
            let sign = if (&routed - &routed_p).expand().is_zero() {
                1
            } else if (&routed + &routed_p).expand().is_zero() {
                -1
            } else {
                return Err(eyre!(
                    "D09 boundary is not the certified signed external momentum"
                ));
            };
            signs.insert(edge, sign);
        }
        let velocity = ctx.settings.on_shell_velocity()?;
        let at_point = |tensor: &SymbolicTensor<IndexedStructure>, lmb: &LoopMomentumBasis| {
            let routed = GS
                .erase_uv_momentum_provenance(tensor.expression())
                .replace_multiple(graph.uv_wrapped_replacement(parent.subgraph(), lmb, &[W_.x___]));
            let mapped = routed.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && momentum.get_nargs() > 0
                    && let Ok(edge) = usize::try_from(momentum.get(0))
                    && let Some(sign) = signs.get(&EdgeIndex(edge))
                {
                    let point = GS.on_shell_momentum(owner, &(Atom::num(*sign) * &mass), &velocity);
                    let indices = momentum
                        .iter()
                        .skip(1)
                        .map(|index| index.to_owned())
                        .collect::<Vec<_>>();
                    **output = GS.indexed_momentum(&point, &indices);
                }
            });
            tensor.with_checked_expression(
                GS.resolve_os_momenta(&mapped)
                    .replace(mass.to_pattern())
                    .with(1),
            )
        };
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector =
            projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
        let projector = at_point(
            &normalize(&projector.try_into_indexed()?)?,
            &graph.loop_momentum_basis,
        )?;
        let phi = |tensor: &SymbolicTensor<IndexedStructure>,
                   lmb: &LoopMomentumBasis|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let projected = projector.clone().multiply(at_point(tensor, lmb)?)?;
            eyre::ensure!(
                projected.structure().canonical().is_scalar(),
                "D09 finalized mass probe leaves ports: projector={:?}; operand={:?}; product={:?}",
                projector.structure().logical_slots(),
                tensor.structure().logical_slots(),
                projected.structure().logical_slots(),
            );
            Ok(projected)
        };
        // Use the existing physical-shell owner to consume loop bookkeeping.
        // Final tensor bodies and scalar denominators retain their separate owners.
        let emit = |sector: &FourDSector| -> Result<SymbolicTensor<IndexedStructure>> {
            let emit_term = |term: FourDTerm| -> Result<SymbolicTensor<IndexedStructure>> {
                let numerator = normalize(&term.numerator)?;
                let denominators =
                    Atom::mul_many(term.denominators.into_iter().map(|denominator| {
                        GS.den(
                            usize::from(denominator.source_edge),
                            denominator.momentum,
                            denominator.mass_squared,
                            denominator.full_expr,
                        )
                    }));
                Ok(numerator
                    .with_checked_expression(numerator.expression() / scalar_tail(&denominators))?)
            };
            let mut terms = sector.physical_terms()?.into_iter();
            let first = terms
                .next()
                .ok_or_else(|| eyre!("D09 physical operand has no rational terms"))?;
            let first = emit_term(first)?;
            terms.try_fold(first, |sum, term| -> Result<_> {
                Ok(sum.try_add(&emit_term(term)?)?)
            })
        };
        // The same finalizer and point map must preserve an independently
        // known mass insertion. Use vertex color slots, dual to the removed
        // external propagators, exactly as in the native amputation controls.
        let spin = boundary
            .iter()
            .map(|hedge| {
                graph.underlying[*hedge]
                    .edge_spin_slots()
                    .next()
                    .map(|slot| slot.to_atom())
                    .ok_or_else(|| eyre!("D09 boundary spin slot missing"))
            })
            .collect::<Result<Vec<_>>>()?;
        let color = boundary
            .iter()
            .map(|hedge| {
                graph.underlying[*hedge]
                    .num_indices
                    .color_indices
                    .vertex_indices
                    .external_structure_iter()
                    .next()
                    .map(|slot| slot.to_atom())
                    .ok_or_else(|| eyre!("D09 boundary vertex color slot missing"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mass_control = AppliedFeynmanRule::new(
            &mass
                * function!(ETS.metric, &spin[0], &spin[1])
                * function!(ETS.metric, &color[0], &color[1]),
        )?
        .expr;
        eyre::ensure!(
            mass_control.structure().canonical().external_structure()
                == local.atom().structure().canonical().external_structure(),
            "the authored mass probe must use the actual parent's indexed interface"
        );
        let projected_mass_control = phi(&normalize(&mass_control)?, &graph.loop_momentum_basis)?;
        let calibration = crate::uv::overlap_numerators::NumeratorJetCache::default()
            .canonicalize(&on_shell_tests::scalar_component(
                projected_mass_control.expression(),
            )?)?;
        eyre::ensure!(
            !calibration.is_zero
                && calibration.key.is_one()
                && calibration.scalar_prefactor.is_one(),
            "D09 finalized mass projector must give exactly one on the native mass identity"
        );
        // A single explicit coordinate witness proves that these reconstructed
        // scalar kernels are not identically zero. No D-dimensional integral
        // or actual Projected4dSector transport is inferred from this check.
        let kernel_at_origin = |kernel: &Atom, loops: usize| -> Result<()> {
            eyre::ensure!(!kernel.is_zero());
            let origin = kernel
                .replace(function!(GS.localizing_integrand, W_.a_))
                .with(W_.a_)
                .replace(function!(GS.emr_mom, W_.a_, W_.b_))
                .with(0)
                .replace(GS.renormalization_localization_scale)
                .with(1);
            let inverse_origin =
                Atom::var(symbolica::atom::Symbol::PI).pow(i64::try_from(2 * loops)?);
            eyre::ensure!(
                (&origin * &inverse_origin)
                    .collect_compact_factors()
                    .is_one()
            );
            Ok(())
        };
        let co_graph = grow_factor(ctx, parent, given)?;
        let mut domains = std::collections::HashSet::new();
        let mut matched = BTreeSet::new();
        let mut frozen_domains = 0;
        for input in prefix.sectors() {
            eyre::ensure!(
                !input.atom.is_zero(),
                "D09 cannot match a pruned input by position"
            );
            let single = Full4dCts(FourDSectors::new(vec![input.clone()], Vec::new())?);
            let coordinate_lmb = compatible_lmb(ctx, parent, given, &single)?;
            let retained = input
                .active_components
                .iter()
                .flat_map(|(_, _, lmb)| lmb.loop_edges.iter())
                .chain(
                    input
                        .frozen_lmbs
                        .iter()
                        .flat_map(|lmb| lmb.loop_edges.iter()),
                )
                .copied()
                .collect::<BTreeSet<_>>();
            let mut quotient = coordinate_lmb.clone();
            if !given.subgraph().is_empty() {
                for index in (0..quotient.loop_edges.len()).rev() {
                    let index = LoopIndex(index);
                    if retained.contains(&quotient.loop_edges[index]) {
                        quotient.put_loop_to_ext(index);
                    }
                }
                graph.canonicalize_lmb_external_order(&mut quotient);
            }
            eyre::ensure!(
                quotient.loop_edges.len()
                    == graph.n_loops(parent.subgraph()) - graph.n_loops(given.subgraph())
            );
            let mut components = input.active_components.clone();
            components.push((
                parent.reduced_subgraph(given),
                parent.subgraph().clone(),
                quotient,
            ));
            let key = (components, input.frozen_lmbs.clone());
            eyre::ensure!(
                domains.insert(key.clone()),
                "D09 input domains must be distinct"
            );
            let outputs = local
                .active_sectors()
                .iter()
                .enumerate()
                .filter(|(_, output)| {
                    output.active_components == key.0 && output.frozen_lmbs == key.1
                })
                .collect::<Vec<_>>();
            eyre::ensure!(
                outputs.len() == 1,
                "D09 source lacks one exact full-domain output"
            );
            let (index, output) = outputs[0];
            eyre::ensure!(matched.insert(index), "D09 output domain was reused");
            let source = FourDSector {
                atom: finalize(input.atom.clone().multiply(co_graph.clone())?)?,
                ..input.clone()
            };
            let source = emit(&source)?;
            let local = emit(output)?;
            let frozen = scalar_tail(&Atom::mul_many(
                input
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| GS.localizing_integrand(lmb)),
            ));
            let frozen_tensor = SymbolicTensor::infer_indexed(frozen.clone())?;
            eyre::ensure!(at_point(&frozen_tensor, &coordinate_lmb)?.expression() == &frozen);
            kernel_at_origin(
                &frozen,
                input
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| lmb.loop_edges.len())
                    .sum(),
            )?;
            if !input.frozen_lmbs.is_empty() {
                frozen_domains += 1;
                eyre::ensure!(input.active_components.is_empty());
                eyre::ensure!(input.frozen_lmbs.as_slice() == std::slice::from_ref(given.lmb()));
                eyre::ensure!(
                    !normalize(&input.atom)?.is_zero() && !source.is_zero(),
                    "actual finite-child carrier must survive finalization"
                );
            }
            let source = phi(&source, &coordinate_lmb)?;
            let local = phi(&local, &coordinate_lmb)?;
            let difference = source.try_add(&local)?;
            on_shell_tests::assert_rational_zero(
                "D09 separately finalized source and signed parent local output",
                graph,
                parent,
                &(difference.expression() * &frozen),
                on_shell_tests::SpinReduction::FiniteComponents,
            )?;
            // An omitted/misnormalized mass subtraction cannot hide behind
            // the matched zero: the independently finalized wrong local
            // operand must leave the unit times the same nonzero kernel.
            let perturbed = FourDSector {
                atom: output.atom.try_add(&mass_control)?,
                ..output.clone()
            };
            let perturbed = phi(&emit(&perturbed)?, &coordinate_lmb)?;
            let wrong = source.try_add(&perturbed)?;
            let wrong = wrong.with_checked_expression(wrong.expression() * &frozen)?;
            let expected = projected_mass_control
                .with_checked_expression(projected_mass_control.expression() * &frozen)?;
            let controlled_residual = wrong.try_sub(&expected)?;
            on_shell_tests::assert_rational_zero(
                "D09 finalized mass perturbation leaves the certified nonzero frozen kernel",
                graph,
                parent,
                controlled_residual.expression(),
                on_shell_tests::SpinReduction::FiniteComponents,
            )?;
        }
        eyre::ensure!(matched.len() == local.active_sectors().len());
        eyre::ensure!(domains.len() == if given.subgraph().is_empty() { 1 } else { 2 });
        eyre::ensure!(frozen_domains == usize::from(!given.subgraph().is_empty()));
        let physical_finite = integrated.physical_finite_counterterm()?;
        let finite = normalize(&physical_finite)?;
        if given.subgraph().is_empty() {
            eyre::ensure!(
                !finite.is_zero(),
                "actual nonzero finite field correction must survive finalization"
            );
        }
        let parent_kernel = scalar_tail(&GS.localizing_integrand(parent.lmb()));
        kernel_at_origin(&parent_kernel, parent.lmb().loop_edges.len())?;
        let kernel_tensor = SymbolicTensor::infer_indexed(parent_kernel.clone())?;
        eyre::ensure!(
            at_point(&kernel_tensor, &graph.loop_momentum_basis)?.expression() == &parent_kernel
        );
        // Exercise the actual recursive finite-carrier owner relative to the
        // actual Integrated result, including consumed loop bookkeeping.
        // A zero result has no surviving recursive payload; it cannot be zipped.
        let completion = Full4dCts::recursion_input(
            local,
            integrated,
            ApproximationType::OS,
            false,
            parent.lmb(),
        )?;
        let completions = &completion.0.recursive_completion;
        if physical_finite.is_zero() {
            eyre::ensure!(completions.is_empty());
        } else {
            eyre::ensure!(completions.len() == 1);
            let completion = &completions[0];
            eyre::ensure!(completion.active_components.is_empty());
            eyre::ensure!(completion.frozen_lmbs.as_slice() == std::slice::from_ref(parent.lmb()));
            let transported = emit(completion)?;
            Self::certify_parent_tensor_zero(
                "D09 signed recursive finite carrier matches the actual physical Integrated output",
                &transported.try_sub(&finite)?,
            )?;
            let kernel = scalar_tail(&Atom::mul_many(
                completion
                    .frozen_lmbs
                    .iter()
                    .map(|lmb| GS.localizing_integrand(lmb)),
            ));
            eyre::ensure!(kernel == parent_kernel);
        }
        let finite_projection = phi(&finite, &graph.loop_momentum_basis)?;
        on_shell_tests::assert_rational_zero(
            "D09 actual finalized finite parent correction",
            graph,
            parent,
            &(finite_projection.expression() * &parent_kernel),
            on_shell_tests::SpinReduction::FiniteComponents,
        )?;
        let perturbed_finite = normalize(&physical_finite.try_add(&mass_control)?)?;
        let wrong = phi(&perturbed_finite, &graph.loop_momentum_basis)?;
        let wrong = wrong.with_checked_expression(wrong.expression() * &parent_kernel)?;
        let expected = projected_mass_control
            .with_checked_expression(projected_mass_control.expression() * &parent_kernel)?;
        on_shell_tests::assert_rational_zero(
            "D09 finalized finite mass perturbation leaves the certified nonzero parent kernel",
            graph,
            parent,
            wrong.try_sub(&expected)?.expression(),
            on_shell_tests::SpinReduction::FiniteComponents,
        )?;
        // First emitted prefix: construct the actual parent CFF and certify
        // its retained boundaries. No source/local mass cancellation is asserted.
        {
            use crate::{
                graph::feynman_graph::FeynmanGraph,
                settings::global::OrientationPattern,
                uv::{
                    Integrands,
                    approx::{
                        CutStructure, OrientationProjection,
                        local_3d::Localizer,
                        projected_4d::{Local4dProjectionContext, Projected4dApproximation},
                    },
                    settings::FinalIntegrandDimension,
                },
            };
            use symbolica::id::Replacement;

            let retained_local = local.clone();
            let retained_finite = physical_finite.clone();
            let mut projection_graph = graph.clone();
            let mut settings = ctx.settings.clone();
            settings.local_uv_cts_from_expanded_4d_integrands = true;
            settings.final_integrand = FinalIntegrandDimension::ThreeD;
            eyre::ensure!(settings.generate_integrated && !settings.add_marker);
            let pattern = OrientationPattern::default();
            let numerator = projection_graph
                .numerator(&graph.full_filter(), &graph.empty_subgraph())?
                .state
                .expr;
            let options = projection_graph.denominator_only_cff_3d_expression_options();
            let surfaces = projection_graph.get_esurface_canonization(&graph.loop_momentum_basis);
            let contracted_edges =
                graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
            let production = projection_graph.generate_3d_expression_for_integrand(
                &contracted_edges,
                &surfaces,
                &options,
                Some(numerator.expression()),
            )?;
            let orientation =
                OrientationProjection::exact_expression(&production, &options, &pattern, true);
            eyre::ensure!(!orientation.orientation_ids()?.is_empty());
            let cuts = CutStructure::empty(&projection_graph);
            let localizer = Localizer::new(&cuts.cuts[0], orientation);
            let allowed_keys = cuts.cuts[0].residue_selector.generate_allowed_keys();
            eyre::ensure!(!allowed_keys.is_empty());

            let family = symbol!("gammalooprs::uv::numerator_family");
            let calls = |integrands: &Integrands| {
                let mut calls = BTreeSet::new();
                for (_, root) in integrands.iter() {
                    root.visitor(&mut |view| {
                        if let AtomView::Fun(call) = view
                            && call.get_symbol() == family
                        {
                            calls.insert(view.to_owned());
                        }
                        true
                    });
                }
                calls
            };
            let check_store = |integrands: &Integrands| -> Result<()> {
                // Reuse the store's flatness, binding and collision admission.
                let checked = integrands
                    .clone()
                    .with_numerators(integrands.numerators().iter().cloned())?;
                eyre::ensure!(&checked == integrands);
                for call in calls(integrands) {
                    let AtomView::Fun(call) = call.as_view() else {
                        unreachable!()
                    };
                    let matches = integrands
                        .numerators()
                        .iter()
                        .filter(|entry| {
                            call.get_nargs() == entry.tags.len() + entry.args.len()
                                && entry
                                    .tags
                                    .iter()
                                    .enumerate()
                                    .all(|(index, tag)| call.get(index) == tag.as_view())
                        })
                        .count();
                    eyre::ensure!(
                        matches == 1,
                        "D09 CFF call lacks one complete retained binding"
                    );
                }
                for (_, root) in integrands.iter() {
                    eyre::ensure!(
                        ![AGS.gamma, AGS.gamma0, GS.u, GS.ubar]
                            .into_iter()
                            .any(|head| root.contains_symbol(head)),
                        "D09 CFF dispatch root contains an explicit spin tensor"
                    );
                }
                Ok(())
            };

            // The producer emits one sector per already grouped full domain,
            // in this stable order. Never zip the ungrouped Taylor sectors.
            let groups = local.projection_sectors()?;
            eyre::ensure!(!groups.is_empty());
            let mut full_domains = std::collections::HashSet::new();
            for group in &groups {
                eyre::ensure!(
                    full_domains
                        .insert((group.active_components.clone(), group.frozen_lmbs.clone())),
                    "D09 projection retained duplicate full domains"
                );
            }
            let mut projection_context = Local4dProjectionContext::default();
            eyre::ensure!(projection_context.source_preparation_builds == 0);
            eyre::ensure!(projection_context.denominator_reconstruction_certificates == 0);
            projection_context.certify_denominator_reconstruction = true;
            let projected =
                Projected4dApproximation::new(localizer, &mut projection_graph, &settings)
                    .project_local_4d(local, &mut projection_context)?;
            eyre::ensure!(projected.sectors().len() == groups.len());
            let certificates = projection_context.denominator_reconstruction_certificates;
            eyre::ensure!(certificates > 0);
            eyre::ensure!(certificates == projection_context.source_preparation_builds);
            eyre::ensure!(projection_context.projected_component_requests >= certificates);
            let mut projected_families = 0;
            for (index, group) in groups.iter().enumerate() {
                let sector = &projected.sectors()[index];
                let expected_kernel = group.frozen_lmbs.iter().fold(Atom::one(), |product, lmb| {
                    product * GS.localizing_integrand(lmb)
                });
                eyre::ensure!(sector.frozen_factor == expected_kernel);
                // The preceding rational gate supplies the nonzero origin
                // witness for these exact same frozen LMB kernels.
                eyre::ensure!(!sector.frozen_factor.is_zero());
                kernel_at_origin(
                    &scalar_tail(&sector.frozen_factor),
                    group
                        .frozen_lmbs
                        .iter()
                        .map(|lmb| lmb.loop_edges.len())
                        .sum(),
                )?;
                check_store(&sector.coefficient)?;
                projected_families += sector.coefficient.numerators().len();
            }
            eyre::ensure!(projected_families > 0);

            // Certify the actual finite localizer's ordered carrier transport
            // relative to this very Integrated output, before finalization/Phi.
            eyre::ensure!(
                !physical_finite
                    .expression()
                    .contains_symbol(GS.integrated_loop_scale)
                    && !physical_finite
                        .expression()
                        .contains_symbol(GS.localizing_integrand)
            );
            let localized = localizer.localize(&physical_finite, &mut projection_graph, parent)?;
            let mut finite_calls = 0;
            if !physical_finite.is_zero() {
                let kernel = GS.localizing_integrand(parent.lmb());
                eyre::ensure!(!kernel.is_zero());
                eyre::ensure!(localized.frozen_integrands.iter().count() == allowed_keys.len());
                for (key, actual_kernel) in localized.frozen_integrands.iter() {
                    eyre::ensure!(allowed_keys.contains(key) && actual_kernel == &kernel);
                    kernel_at_origin(&scalar_tail(actual_kernel), parent.lmb().loop_edges.len())?;
                }
                for (host, source_map, integrands) in localized.active.iter_orientations() {
                    check_store(integrands)?;
                    eyre::ensure!(integrands.numerators().len() == 1);
                    let entry = &integrands.numerators()[0];
                    let expected = orientation.map_numerator(
                        &projection_graph,
                        host,
                        source_map,
                        physical_finite.expression(),
                    )?;
                    let expected = physical_finite.with_checked_expression(expected)?;
                    for call in calls(integrands) {
                        let AtomView::Fun(call) = call.as_view() else {
                            unreachable!()
                        };
                        let replacements = entry.args.iter().enumerate().map(|(index, formal)| {
                            Replacement::new(
                                Atom::from(formal.clone()).to_pattern(),
                                call.get(entry.tags.len() + index).to_owned(),
                            )
                        });
                        // Specialize one retained typed body as a proof view;
                        // never substitute the body into its dispatch roots.
                        let transported = entry.rhs.with_checked_expression(
                            entry.rhs.expression().replace_multiple(replacements),
                        )?;
                        Self::certify_parent_tensor_zero(
                            "D09 actual CFF finite carrier and ordered branch arguments",
                            &transported.try_sub(&expected)?,
                        )?;
                        finite_calls += 1;
                    }
                }
                eyre::ensure!(finite_calls > 0);
            }
            // This is the actual final builder with both retained operands.
            // into_integrands keeps every numerator body separate from roots.
            let finalized = FinalIntegrandBuilder::new(localizer, &settings)
                .build_projected(&mut projection_graph, parent, &projected, integrated)?
                .into_integrands();
            check_store(&finalized)?;
            eyre::ensure!(finalized.scalar_definitions().is_empty());
            eyre::ensure!(!finalized.numerators().is_empty());
            eyre::ensure!(finalized.iter().any(|(_, root)| !root.is_zero()));
            eyre::ensure!(finalized.iter().count() == allowed_keys.len());
            eyre::ensure!(finalized.iter().all(|(key, _)| allowed_keys.contains(key)));
            eyre::ensure!(local == &retained_local);
            eyre::ensure!(integrated.physical_finite_counterterm()? == retained_finite);
            debug_tags!(#uv, #local, #four_d, #cff, #integrated, #summary;
                stage = "native_d09_parent_cff_reconstruction",
                source = %graph.name, given = %given.subgraph().string_label(),
                full_domains = groups.len(), certificates,
                file.grouped_domains = ?groups.iter().map(|group| (&group.active_components, &group.frozen_lmbs)).collect::<Vec<_>>(),
                component_requests = projection_context.projected_component_requests,
                projected_families, finite_calls,
                final_families = finalized.numerators().len(),
                final_root_bytes = finalized.iter().map(|(_, atom)| atom.as_view().get_byte_size()).sum::<usize>(),
                final_body_bytes = finalized.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
                actual_parent_cff_constructed = true, finite_carrier_transport_checked = true,
                emitted_mass_cancellation_proved = false, dimensional_R2_condition_proved = false,
                "Actual D09 parent CFF retains its certified factorized operands"
            );
            // Construct the actual source lanes next. Keep the external matrix
            // open in this separate projection view; no mass projector or
            // matched-point specialization enters source/assignment selection.
            {
                use crate::numerator::symbolica_ext::NumeratorAtomExt;
                use crate::uv::approx::direct_3d::{Direct3dCts, DirectResidueBranches};
                use std::{collections::HashMap, sync::Arc};

                let retained_prefix = prefix.clone();
                let retained_global_projector = projection_graph.global_prefactor.projector.clone();
                let mut source_graph = projection_graph.clone();
                source_graph.global_prefactor.projector =
                    SymbolicTensor::infer_indexed(Atom::one())?;
                eyre::ensure!(
                    source_graph.global_prefactor.num == projection_graph.global_prefactor.num
                );
                eyre::ensure!(source_graph.overall_factor == projection_graph.overall_factor);
                eyre::ensure!(
                    source_graph.loop_momentum_basis == projection_graph.loop_momentum_basis
                );
                // The production expression above uses this full graph's
                // numerator without its external global projector. Reusing it
                // retains the actual stored energy bounds/maps for this view.
                eyre::ensure!(
                    !retained_global_projector
                        .expression()
                        .contains_symbol(GS.loop_mom)
                );
                let external_owners = boundary
                    .iter()
                    .map(|hedge| graph[hedge].0)
                    .collect::<BTreeSet<_>>();
                let mut unexpected_external = Vec::new();
                retained_global_projector.expression().visitor(&mut |view| {
                    if let AtomView::Fun(factor) = view
                        && [GS.u, GS.ubar, GS.emr_mom].contains(&factor.get_symbol())
                        && (factor.get_nargs() == 0
                            || usize::try_from(factor.get(0))
                                .ok()
                                .is_none_or(|edge| !external_owners.contains(&edge)))
                    {
                        unexpected_external.push(view.to_owned());
                    }
                    true
                });
                eyre::ensure!(
                    unexpected_external.is_empty(),
                    "the external projection view cannot remove a hard-loop-dependent factor: {unexpected_external:?}"
                );
                // Reuse the already certified parent coefficients in the SAME
                // open external graph view as the source. This is an actual
                // final-builder call with the original Integrated result.
                let open_parent = FinalIntegrandBuilder::new(localizer, &settings)
                    .build_projected(&mut source_graph, parent, &projected, integrated)?
                    .into_integrands();
                check_store(&open_parent)?;
                eyre::ensure!(open_parent.scalar_definitions().is_empty());
                eyre::ensure!(!open_parent.numerators().is_empty());
                eyre::ensure!(open_parent.iter().count() == allowed_keys.len());
                eyre::ensure!(
                    open_parent
                        .iter()
                        .all(|(key, _)| allowed_keys.contains(key))
                );
                let source_result = if given.subgraph().is_empty() {
                    let direct = Direct3dCts::root(&source_graph, localizer)?;
                    FinalIntegrandBuilder::new(localizer, &settings)
                        .build_direct(&mut source_graph, given, &direct, &IntegratedCts::root())?
                        .into_integrands()
                } else {
                    let sectors = prefix.sectors().collect::<Vec<_>>();
                    let active = sectors
                        .iter()
                        .filter(|sector| sector.frozen_lmbs().is_empty())
                        .map(|sector| (*sector).clone())
                        .collect::<Vec<_>>();
                    let finite = sectors
                        .iter()
                        .filter(|sector| !sector.frozen_lmbs().is_empty())
                        .copied()
                        .collect::<Vec<_>>();
                    eyre::ensure!(active.len() == 1 && finite.len() == 1);
                    eyre::ensure!(finite[0].active_components.is_empty());
                    eyre::ensure!(finite[0].frozen_lmbs() == std::slice::from_ref(given.lmb()));
                    eyre::ensure!(graph.n_loops(given.subgraph()) == 1);
                    let active_source = Local4dCts(FourDSectors::new(active, Vec::new())?);
                    Self::certify_parent_tensor_zero(
                        "D09 actual source prefix active/finite recomposition",
                        &prefix
                            .atom()
                            .try_sub(&active_source.atom().try_add(&finite[0].atom)?)?,
                    )?;
                    let mut physical_child_finite =
                        finite[0].atom.with_checked_expression(Atom::Zero)?;
                    let terms = finite[0].physical_terms()?;
                    eyre::ensure!(!terms.is_empty());
                    for term in terms {
                        eyre::ensure!(
                            term.denominators.is_empty(),
                            "a retained finite child cannot own physical propagators"
                        );
                        physical_child_finite = physical_child_finite.try_add(&term.numerator)?;
                    }
                    eyre::ensure!(!physical_child_finite.is_zero());
                    eyre::ensure!(
                        !physical_child_finite
                            .expression()
                            .contains_symbol(GS.localizing_integrand)
                    );
                    let bookkeeping = physical_child_finite.with_checked_expression(
                        physical_child_finite.expression()
                            * Atom::var(GS.integrated_loop_scale).pow(4),
                    )?;
                    Self::certify_parent_tensor_zero(
                        "D09 actual finite child retains exactly its consumed one-loop measure",
                        &finite[0].atom.try_sub(&bookkeeping)?,
                    )?;
                    let active_groups = active_source.projection_sectors()?;
                    let mut child_context = Local4dProjectionContext {
                        certify_denominator_reconstruction: true,
                        ..Default::default()
                    };
                    eyre::ensure!(child_context.source_preparation_builds == 0);
                    eyre::ensure!(child_context.denominator_reconstruction_certificates == 0);
                    let child_projected =
                        Projected4dApproximation::new(localizer, &mut source_graph, &settings)
                            .project_local_4d(&active_source, &mut child_context)?;
                    eyre::ensure!(child_projected.sectors().len() == active_groups.len());
                    eyre::ensure!(child_context.source_preparation_builds > 0);
                    eyre::ensure!(
                        child_context.denominator_reconstruction_certificates
                            == child_context.source_preparation_builds
                    );
                    for (ordinal, group) in active_groups.iter().enumerate() {
                        let expected = group
                            .frozen_lmbs()
                            .iter()
                            .fold(Atom::one(), |product, lmb| {
                                product * GS.localizing_integrand(lmb)
                            });
                        eyre::ensure!(child_projected.sectors()[ordinal].frozen_factor == expected);
                    }
                    // This zero means no ADDITIONAL finite child term here.
                    // The actual retained finite completion is localized below.
                    let active_result = FinalIntegrandBuilder::new(localizer, &settings)
                        .build_projected(
                            &mut source_graph,
                            given,
                            &child_projected,
                            &IntegratedCts::root(),
                        )?
                        .into_integrands();
                    check_store(&active_result)?;
                    let finite_localized =
                        localizer.localize(&physical_child_finite, &mut source_graph, given)?;
                    let actual_kernel = GS.localizing_integrand(given.lmb());
                    eyre::ensure!(!actual_kernel.is_zero());
                    eyre::ensure!(
                        finite_localized.frozen_integrands.iter().count() == allowed_keys.len()
                    );
                    for (key, kernel) in finite_localized.frozen_integrands.iter() {
                        eyre::ensure!(allowed_keys.contains(key) && kernel == &actual_kernel);
                    }
                    let finite_branches =
                        DirectResidueBranches::from_transient(&finite_localized.combine()?)?;
                    let mut finite_entries = HashMap::new();
                    for (key, integrands) in finite_branches.iter_keys() {
                        check_store(integrands)?;
                        eyre::ensure!(integrands.numerators().len() == 1);
                        eyre::ensure!(
                            finite_entries
                                .insert(key.clone(), Arc::clone(&integrands.numerators()[0]))
                                .is_none()
                        );
                    }
                    let reduced_source = graph
                        .full_filter()
                        .subtract(given.subgraph())
                        .subtract(&graph.initial_state_cut);
                    let outside = source_graph
                        .numerator(&reduced_source, given.subgraph())?
                        .state
                        .expr
                        .simplify_algebra(&idenso::tensor::AlgebraSettings {
                            color: Some(ColorSimplifySettings::default()),
                            ..Default::default()
                        })?
                        .multiply(source_graph.global_numerator()?)?;
                    eyre::ensure!(
                        !outside.is_scalar(),
                        "the actual D09 finite child has a tensor cograph"
                    );
                    let complete_source =
                        physical_child_finite.clone().multiply(outside.clone())?;
                    eyre::ensure!(
                        complete_source.structure().canonical().external_structure()
                            == numerator.structure().canonical().external_structure()
                    );
                    let completed_branches = finite_branches.multiply_key_mapped(
                        orientation,
                        &source_graph,
                        &outside,
                        DirectResidueBranches::numerator_scope(),
                    )?;
                    let mut witnessed = 0usize;
                    let mut source_keys = 0usize;
                    for (key, integrands) in completed_branches.iter_keys() {
                        source_keys += 1;
                        check_store(integrands)?;
                        let finite_entry = finite_entries
                            .get(key)
                            .ok_or_else(|| eyre!("finite child lost its exact residue key"))?;
                        eyre::ensure!(
                            integrands
                                .numerators()
                                .iter()
                                .any(|entry| Arc::ptr_eq(entry, finite_entry))
                        );
                        let outside_entries = integrands
                            .numerators()
                            .iter()
                            .filter(|entry| !Arc::ptr_eq(entry, finite_entry))
                            .collect::<Vec<_>>();
                        eyre::ensure!(
                            outside_entries.len() == 1 && integrands.numerators().len() == 2
                        );
                        let outside_entry = outside_entries[0];
                        let actual_calls = calls(integrands).into_iter().collect::<Vec<_>>();
                        for (cut, root) in integrands.iter() {
                            let Some(rows) = root.coefficient_list_exact(&actual_calls) else {
                                debug_tags!(#uv, #local, #cff, #summary;
                                    stage = "native_d09_source_monomial_refused", file.residue_key = ?key,
                                    file.cut = ?cut, file.root = %root.to_canonical_string(),
                                    calls = actual_calls.len(), root_bytes = root.as_view().get_byte_size(),
                                    "Existing bounded complete-family-call polynomial owner refused; no fallback"
                                );
                                return Err(eyre!(
                                    "D09 source complete-call polynomial extraction refused"
                                ));
                            };
                            let reconstructed = rows
                                .iter()
                                .map(|(monomial, coefficient)| monomial * coefficient)
                                .sum::<Atom>();
                            let residual = (root - &reconstructed)
                                .coefficient_list_exact(&actual_calls)
                                .ok_or_else(|| {
                                    eyre!("D09 source complete-call reconstruction refused")
                                })?;
                            eyre::ensure!(
                                residual
                                    .iter()
                                    .all(|(_, coefficient)| coefficient.is_zero())
                            );
                            for (monomial, coefficient) in rows {
                                eyre::ensure!(
                                    !coefficient.contains_symbol(family),
                                    "D09 complete-call extraction left an unbound family in its coefficient"
                                );
                                if coefficient.is_zero() {
                                    continue;
                                }
                                let AtomView::Mul(product) = monomial.as_view() else {
                                    return Err(eyre!(
                                        "D09 actual finite source did not retain its complete two-factor call monomial"
                                    ));
                                };
                                let factors = product
                                    .iter()
                                    .map(|factor| factor.to_owned())
                                    .collect::<Vec<_>>();
                                eyre::ensure!(factors.len() == 2);
                                let matches = |call: &Atom, entry: &crate::integrands::process::param_builder::NumeratorDefinition| {
                                    matches!(call.as_view(), AtomView::Fun(call) if call.get_symbol() == family
                                        && call.get_nargs() == entry.tags.len() + entry.args.len()
                                        && entry.tags.iter().enumerate().all(|(index, tag)| call.get(index) == tag.as_view()))
                                };
                                let child_call = factors
                                    .iter()
                                    .filter(|call| matches(call, finite_entry))
                                    .collect::<Vec<_>>();
                                let outside_call = factors
                                    .iter()
                                    .filter(|call| matches(call, outside_entry))
                                    .collect::<Vec<_>>();
                                eyre::ensure!(
                                    child_call.len() == 1
                                        && outside_call.len() == 1
                                        && child_call[0] != outside_call[0]
                                );
                                // Every actual nonzero-coefficient monomial,
                                // under every actual key/cut, gets the typed
                                // transport proof before materialization.
                                {
                                    let bind = |call: &Atom, entry: &crate::integrands::process::param_builder::NumeratorDefinition| -> Result<_> {
                                        let AtomView::Fun(function) = call.as_view() else { unreachable!() };
                                        let replacements = entry.args.iter().enumerate().map(|(index, formal)| Replacement::new(
                                            Atom::from(formal.clone()).to_pattern(), function.get(entry.tags.len() + index).to_owned(),
                                        )).collect::<Vec<_>>();
                                        eyre::ensure!(entry.lhs.replace_multiple(&replacements) == *call);
                                        Ok(entry.rhs.with_checked_expression(entry.rhs.expression().replace_multiple(&replacements))?)
                                    };
                                    let child_body = bind(child_call[0], finite_entry)?;
                                    let outside_body = bind(outside_call[0], outside_entry)?;
                                    for (label, body, original) in [
                                        (
                                            "actual finite child call",
                                            &child_body,
                                            &physical_child_finite,
                                        ),
                                        (
                                            "actual finite child cograph call",
                                            &outside_body,
                                            &outside,
                                        ),
                                    ] {
                                        let expected = original.with_checked_expression(
                                            key.map_numerator(
                                                orientation,
                                                &source_graph,
                                                original.expression(),
                                            )?,
                                        )?;
                                        Self::certify_parent_tensor_zero(
                                            label,
                                            &body.try_sub(&expected)?,
                                        )?;
                                    }
                                    let complete_body = child_body.multiply(outside_body)?;
                                    eyre::ensure!(
                                        complete_body.structure().canonical().external_structure()
                                            == complete_source
                                                .structure()
                                                .canonical()
                                                .external_structure()
                                    );
                                    let expected = complete_source.with_checked_expression(
                                        key.map_numerator(
                                            orientation,
                                            &source_graph,
                                            complete_source.expression(),
                                        )?,
                                    )?;
                                    Self::certify_parent_tensor_zero(
                                        "actual same-key finite child/cograph complete monomial",
                                        &complete_body.try_sub(&expected)?,
                                    )?;
                                    debug_tags!(#uv, #local, #four_d, #cff, #summary;
                                        stage = "native_d09_source_complete_monomial", file.residue_key = ?key,
                                        file.cut = ?cut, file.monomial = %monomial.to_canonical_string(),
                                        file.child_call = %child_call[0].to_canonical_string(),
                                        file.cograph_call = %outside_call[0].to_canonical_string(),
                                        file.frozen_domain = ?finite[0].frozen_lmbs(),
                                        file.active_domains = ?active_groups.iter().map(|group| (&group.active_components, group.frozen_lmbs())).collect::<Vec<_>>(),
                                        file.kernel = %actual_kernel.to_canonical_string(),
                                        whole_parent_rank = complete_body.rank(), actual_call_pair_bound = true,
                                        mass_projector_applied = false,
                                        "Two actual partial families form one native open parent matrix under the same exact source key"
                                    );
                                    witnessed += 1;
                                }
                            }
                        }
                    }
                    eyre::ensure!(source_keys == finite_entries.len() && witnessed > 0);
                    let normalize_source =
                        FinalIntegrandBuilder::final_normalizer(&source_graph, &reduced_source);
                    let normalize_source_scalar =
                        FinalIntegrandBuilder::scalar_normalizer(&source_graph, &reduced_source);
                    let finite_result = completed_branches
                        .resolved_scalars()?
                        .fallible_map(|_, atom| {
                            Ok(normalize_source_scalar(atom)
                                .replace(GS.m_uv_expansion)
                                .with(GS.m_uv_vacuum)
                                .replace(GS.dim_epsilon)
                                .with(0))
                        })?
                        .map_numerators(|entry| normalize_source(&entry.rhs))?
                        .materialize(false)?;
                    check_store(&finite_result)?;
                    active_result.zip_add([finite_result])?
                };
                check_store(&source_result)?;
                eyre::ensure!(source_result.scalar_definitions().is_empty());
                eyre::ensure!(!source_result.numerators().is_empty());
                eyre::ensure!(source_result.iter().count() == allowed_keys.len());
                eyre::ensure!(
                    source_result
                        .iter()
                        .all(|(key, _)| allowed_keys.contains(key))
                );
                eyre::ensure!(prefix == &retained_prefix);
                eyre::ensure!(local == &retained_local);
                eyre::ensure!(integrated.physical_finite_counterterm()? == retained_finite);
                eyre::ensure!(
                    projection_graph.global_prefactor.projector == retained_global_projector
                );
                debug_tags!(#uv, #local, #four_d, #cff, #summary;
                    stage = "native_d09_actual_source_cff", given = %given.subgraph().string_label(),
                    families = source_result.numerators().len(),
                    open_parent_families = open_parent.numerators().len(),
                    root_bytes = source_result.iter().map(|(_, root)| root.as_view().get_byte_size()).sum::<usize>(),
                    body_bytes = source_result.numerators().iter().map(|entry| entry.rhs.expression().as_view().get_byte_size()).sum::<usize>(),
                    actual_source_constructed = true, projection_only_open_external_matrix = true,
                    analytic_backend_repeated = false, emitted_mass_cancellation_proved = false,
                    "Actual D09 source CFF retains local and finite child contributions; no Phi or full mass claim"
                );
            }
        }
        debug_tags!(#uv, #integrated, #summary;
            stage = "native_d09_rational_finalizer_mass_projection",
            source = %diagram.name(), given = %given.subgraph().string_label(),
            domains = domains.len(), frozen_domains,
            per_full_domain_source_local_zero = true, actual_finalized_finite_mass_zero = true,
            mass_projector_unit_calibration = true, mass_perturbation_nonzero_control = true,
            signed_parent_carrier_transport = true, reconstructed_kernels_nonzero = true,
            child_epsilon_depth_checked_separately = true,
            emitted_cff_condition_proved = false, dimensional_R2_condition_proved = false,
            "Actual D09 operands pass the matched physical rational-finalizer boundary"
        );
        Ok(())
    }

    pub(in crate::uv) fn certify_parent_mass_projection_input(
        &self,
        diagram: &feynkit_graph::FeynmanDiagram,
        ctx: &UVCtx<'_>,
        parent: &OwnedForestNode,
        given: &OwnedForestNode,
        prefix: &Full4dCts,
    ) -> Result<Local4dCts> {
        let reference = self;
        let graph = ctx.graph;
        let child = given;
        eyre::ensure!(
            parent.renormalization_scheme() == ApproximationType::OS && parent.dod() == 1
        );
        eyre::ensure!(graph.n_loops(parent.subgraph()) == 2);
        eyre::ensure!(
            graph.paired_edges(parent.subgraph()) == graph.paired_edges(&graph.full_filter())
        );
        eyre::ensure!(ctx.settings.generate_integrated && !ctx.settings.add_marker);
        eyre::ensure!(child.subgraph().is_empty() || graph.n_loops(child.subgraph()) == 1);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let exact_zero = Self::certify_parent_tensor_zero;
        let frame_zero =
            |label: &str, difference: &SymbolicTensor<IndexedStructure>| -> Result<()> {
                exact_zero(
                    label,
                    &Self::reduce_parent_external_frame(label, difference)?,
                )
            };
        assert_eq!(GS.emr_mom, feynkit_graph::symbols::momentum());
        assert_eq!(GS.hedgeaind, feynkit_graph::symbols::hedge_index());
        let native_projector = diagram.projector().clone().try_into_indexed()?;
        assert!(
            graph
                .global_prefactor
                .projector
                .same_value(&native_projector)
        );
        let mut native_ports = Vec::new();
        let projector_source = cooking.uncook(diagram.projector().expression().as_view());
        projector_source.visitor(&mut |view| {
            if let AtomView::Fun(wavefunction) = view
                && [GS.u, GS.ubar].contains(&wavefunction.get_symbol())
                && wavefunction.get_nargs() == 2
                && let Ok(edge) = usize::try_from(wavefunction.get(0))
            {
                native_ports.push((EdgeIndex(edge), wavefunction.get(1).to_owned()));
            }
            true
        });
        assert_eq!(native_ports.len(), 2);
        for (edge, native_port) in native_ports {
            let runtime_ports = graph.underlying[graph[&edge].1.any_hedge()]
                .edge_spin_slots()
                .map(|slot| slot.to_atom())
                .collect::<Vec<_>>();
            assert_eq!(runtime_ports, [native_port]);
        }
        let crown = graph.dummy_less_full_crown(parent.subgraph());
        let boundary = crown.included_iter().collect::<Vec<_>>();
        assert_eq!(boundary.len(), 2);
        let incoming = *boundary
            .iter()
            .find(|h| graph.flow(**h) == Flow::Sink)
            .unwrap();
        let owner = graph[&incoming];
        let mass = graph[owner].particle.mass_atom(&graph.model);
        let external = boundary.iter().map(|h| graph[h]).collect::<Vec<_>>();
        let (projector, _) = reference.mass_projector(diagram, &Atom::var(GS.dim))?;
        // The generic reference promotes its source color ports to Nc.
        // This runtime source instead retains its native SU(3) representations.
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector =
            projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
        let projector = projector.try_into_indexed()?;
        let velocity = ctx.settings.on_shell_velocity()?;
        let lmb = compatible_lmb(ctx, parent, child, prefix)?;
        let source = grow(prefix, ctx, parent, child)?;
        let bare_model_source = graph
            .numerator(parent.subgraph(), &graph.empty_subgraph())?
            .to_d_dim(GS.dim)?
            .state
            .expr;
        eyre::ensure!(
            !bare_model_source.expression().contains_symbol(AGS.gamma0)
                && !bare_model_source.expression().contains_symbol(ETS.delta),
            "the canonical model itself has no fixed subtraction frame"
        );
        if child.renormalization_scheme() == ApproximationType::OS {
            eyre::ensure!(graph.n_loops(child.subgraph()) == 1 && child.dod() == 1);
            eyre::ensure!(prefix.0.has_on_shell_ancestry);
            eyre::ensure!(
                velocity == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero],
                "continuation of the certified primitive OS child is initially rest-only"
            );
        } else {
            eyre::ensure!(
                !source.expression().contains_symbol(AGS.gamma0)
                    && !source.expression().contains_symbol(ETS.delta),
                "a non-OS prefix must not introduce an external fixed frame"
            );
        }
        let at_point = |tensor: &SymbolicTensor<IndexedStructure>,
                        sign: i32|
         -> Result<SymbolicTensor<IndexedStructure>> {
            let atom = tensor.expression();
            let point = GS.on_shell_momentum(owner, &(Atom::num(sign) * &mass), &velocity);
            let routed = GS.erase_uv_momentum_provenance(atom).replace_multiple(
                graph.uv_wrapped_replacement(parent.subgraph(), &lmb, &[W_.x___]),
            );
            let expression = GS
                .resolve_os_momenta(&routed.replace_map(|view, _, output| {
                    if let AtomView::Fun(momentum) = view
                        && momentum.get_symbol() == GS.emr_mom
                        && momentum.get_nargs() > 0
                        && usize::try_from(momentum.get(0))
                            .is_ok_and(|edge| external.contains(&EdgeIndex(edge)))
                    {
                        let indices = momentum
                            .iter()
                            .skip(1)
                            .map(|index| index.to_owned())
                            .collect::<Vec<_>>();
                        **output = GS.indexed_momentum(&point, &indices);
                    }
                }))
                .replace(mass.to_pattern())
                .with(1);
            Ok(tensor.with_checked_expression(expression)?)
        };
        let projector = at_point(&projector, 1)?;
        let phi = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
            Ok(projector.clone().multiply(at_point(tensor, 1)?)?)
        };
        let [positive, negative] =
            graph.on_shell_branches(&graph.model, &source, parent.subgraph(), &lmb, &velocity)?;
        for (label, value) in [
            ("source", phi(&source)?),
            ("positive", phi(&positive)?),
            ("negative", phi(&negative)?),
        ] {
            let admitted = &value;
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_mass_channel", label,
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.velocity = ?velocity, rank = admitted.rank(),
                file.value = %cooking.uncook(value.expression().as_view()).to_canonical_string(),
                "Admitting each actual mass-projected branch before equality"
            );
            eyre::ensure!(
                admitted.structure().canonical().is_scalar(),
                "{label}: mass projector left unmatched native spin/color ports"
            );
        }
        let outgoing = *boundary
            .iter()
            .find(|h| graph.flow(**h) == Flow::Source)
            .unwrap();
        let row = graph.underlying[outgoing]
            .edge_spin_slots()
            .next()
            .unwrap()
            .to_atom();
        let spinor_ports = |atom: &Atom| {
            let mut matcher = SlotMatcher::default();
            let mut ports = BTreeSet::new();
            atom.visitor(&mut |view| {
                if matcher.parse::<Bispinor, Aind>(view).is_ok() {
                    ports.insert(view.to_owned());
                    false
                } else {
                    true
                }
            });
            ports
        };
        // Only this isolated external block is represented in its finite
        // physical basis. Every indexed Lorentz port must be fixed by a
        // delta coordinate 0..3; no internal source or loop vector enters.
        let physical_external_block = |block: &Atom| -> Result<Atom> {
            eyre::ensure!(
                !block.contains_symbol(GS.emr_mom)
                    && !block.contains_symbol(GS.den)
                    && !block.contains_symbol(SPENSO_TAG.trace),
                "the finite block cannot contain internal momentum, poles or traces"
            );
            let mut matcher = SlotMatcher::default();
            let mut lorentz = BTreeSet::new();
            let mut fixed = BTreeSet::new();
            let mut invalid = false;
            block.visitor(&mut |view| {
                if matcher.parse::<Minkowski, Aind>(view).is_ok() {
                    lorentz.insert(view.to_owned());
                    return false;
                }
                if let AtomView::Fun(function) = view {
                    let head = function.get_symbol();
                    invalid |= head.has_tag(&SPENSO_TAG.tensor)
                        && ![ETS.metric, AGS.gamma, AGS.gamma0, GS.delta_vec].contains(&head);
                    if head == GS.delta_vec && function.get_nargs() == 2 {
                        invalid |=
                            !matches!(matcher.concrete_component(function.get(0)), Some(Ok(0..=3)));
                        fixed.insert(function.get(1).to_owned());
                    }
                    if head == CS.fundamental_rep {
                        invalid |= function.get(0) != Atom::num(3).as_view();
                    }
                }
                true
            });
            eyre::ensure!(
                !invalid && lorentz == fixed,
                "external projector must use only its explicit physical 0..3 basis"
            );
            let without_lorentz_ports = block.replace_map(|view, _, output| {
                if let AtomView::Fun(port) = view
                    && port.get_symbol() == LibraryRep::from(Minkowski {}).symbol()
                {
                    **output = Atom::one();
                }
            });
            eyre::ensure!(
                !without_lorentz_ports.contains_symbol(GS.dim),
                "the physical block cannot carry a contracted formal-d Lorentz factor"
            );
            Ok(block.map_mink_dim(4))
        };
        let finite_projector =
            physical_external_block(&cooking.uncook(projector.expression().as_view()))?;
        let block_index = idenso::bis!(4, Atom::from(Aind::new_dummy()));
        eyre::ensure!(
            !cooking
                .uncook(source.expression().as_view())
                .contains(block_index.as_view())
        );
        let identity = function!(
            spenso::network::library::symbolic::ETS.metric,
            &row,
            &block_index
        );
        let gamma = on_shell_tests::gamma_velocity_at_port(graph, &identity, outgoing.0, &velocity);
        let gamma = physical_external_block(&gamma)?;
        for sign in [1, -1] {
            let p = (&identity + Atom::num(sign) * &gamma) / Atom::num(2);
            let expected = if sign == 1 {
                finite_projector
                    .replace(row.to_pattern())
                    .with(block_index.clone())
            } else {
                Atom::Zero
            };
            // This is only the already witnessed finite external block.
            // Native exact library components recognize transposed chain
            // spellings without assigning four dimensions to the source X.
            use spenso::{
                iterators::IteratableTensor,
                network::{ExecutionResult, Sequential, SmallestDegree},
            };
            let difference = physical_external_block(&(&finite_projector * p - expected))?;
            let mut network = difference.as_view().parse_into_net()?;
            let library = crate::utils::TENSORLIB.read().unwrap();
            network.execute::<Sequential, SmallestDegree, _, _, _>(
                &*library,
                &*crate::utils::FUN_LIB,
            )?;
            let stored_entries = match network.result_tensor(&*library)? {
                ExecutionResult::Zero => 0,
                ExecutionResult::One => {
                    return Err(eyre!("finite external projector difference is one"));
                }
                ExecutionResult::Val(tensor) => {
                    assert_eq!(tensor.size()?, 4 * 4 * 3 * 3);
                    let mut entries = 0;
                    for (_, value) in tensor.iter_flat() {
                        eyre::ensure!(
                            value.is_zero(),
                            "finite ordered external Phi P_{sign} has a nonzero exact component: {value}"
                        );
                        entries += 1;
                    }
                    entries
                }
            };
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_finite_external_components", sign,
                graph = %graph.name, child = %child.subgraph().string_label(),
                component_space = 144, stored_entries,
                file.difference = %difference.to_canonical_string(),
                "All native exact spin and SU(3) components of the external projector difference vanish"
            );
        }

        // Extract the identity coefficient of the ACTUAL projector action.
        // Its new spinor port did not exist in X. Setting just the matrices
        // from the graph output to that port to zero recovers X_sigma/2.
        // Reconstruction below certifies this decomposition before the
        // independent source/routing comparison. No child tensor is replaced.
        let certify_branches = |label: &str,
                                input_native: &SymbolicTensor<IndexedStructure>,
                                branches: [&SymbolicTensor<IndexedStructure>; 2]|
         -> Result<()> {
            let untouched_native = input_native.clone();
            let input_source = cooking.uncook(input_native.expression().as_view());
            let input = &input_source;
            let untouched = input.clone();
            let original_ports = spinor_ports(input);
            for (sign, branch_native) in [1, -1].into_iter().zip(branches) {
                let branch_source = cooking.uncook(branch_native.expression().as_view());
                let branch = &branch_source;
                let new_ports = spinor_ports(branch)
                    .difference(&original_ports)
                    .cloned()
                    .collect::<Vec<_>>();
                let [fresh] = new_ports.as_slice() else {
                    return Err(eyre!(
                        "{label}: actual projector must introduce one fresh spinor interface"
                    ));
                };
                let mut selected = 0;
                let identity_part = branch.replace_map(|view, _, output| {
                    if let AtomView::Fun(matrix) = view {
                        let ordered = matrix.get_nargs() >= 2
                            && matrix.get(0) == row.as_view()
                            && matrix.get(1) == fresh.as_view();
                        // Native gamma0 is symmetric; ordinary gamma endpoint
                        // order remains strict, including the spatial block.
                        let temporal = matrix.get_symbol() == AGS.gamma0
                            && matrix.get_nargs() == 2
                            && (ordered
                                || (matrix.get(1) == row.as_view()
                                    && matrix.get(0) == fresh.as_view()));
                        if temporal
                            || (matrix.get_symbol() == AGS.gamma
                                && matrix.get_nargs() == 3
                                && ordered)
                        {
                            selected += 1;
                            **output = Atom::Zero;
                        }
                    }
                });
                eyre::ensure!(
                    selected > 0,
                    "{label}: no actual ordered external projector block found"
                );
                let evaluated = Atom::num(2) * identity_part;
                eyre::ensure!(
                    !evaluated.contains(fresh.as_view()),
                    "{label}: projector coefficient extraction retained the new spinor port"
                );
                let evaluated_native =
                    input_native.with_checked_expression(cooking.try_cook(evaluated.as_view())?)?;
                let row_native = cooking.try_cook(row.as_view())?;
                let row_native = Slot::<LibraryRep, AbstractIndex>::try_from(row_native.as_view())?;
                let fresh_native = cooking.try_cook(fresh.as_view())?;
                let fresh_native =
                    Slot::<LibraryRep, AbstractIndex>::try_from(fresh_native.as_view())?;
                assert_eq!(row_native.rep(), fresh_native.rep());
                let row_position = evaluated_native
                    .structure()
                    .logical_slots()
                    .iter()
                    .position(|slot| slot.to_atom() == row_native.to_atom())
                    .ok_or_else(|| {
                        eyre!("{label}: certified free row is not an actual logical port")
                    })?;
                let renamed_native = evaluated_native.reindex_interface_ports(
                    &std::collections::HashMap::from([(row_position, fresh_native.aind())]),
                )?;
                let renamed = cooking.uncook(renamed_native.expression().as_view());
                assert_eq!(
                    renamed.replace(fresh.to_pattern()).with(row.clone()),
                    evaluated,
                    "only the certified free endpoint may be renamed"
                );
                let identity =
                    function!(spenso::network::library::symbolic::ETS.metric, &row, fresh);
                let block4 =
                    on_shell_tests::gamma_velocity_at_port(graph, &identity, outgoing.0, &velocity);
                let block4 = physical_external_block(&block4)?;
                // Embed ONLY the witnessed physical external block into the
                // unchanged formal-d interface of the actual source.
                let block_d = block4.map_mink_dim(GS.dim);
                let spinor = [LibraryRep::from(Bispinor {})];
                let attachment = ContractSettings {
                    representations: Some(&spinor),
                    rank_one: false,
                    collect_chains: false,
                    collect_traces: false,
                    expand: false,
                    ..Default::default()
                };
                // Remove the helper's auxiliary identity metric before sewing
                // this small block to the unchanged factorized source. Whole-
                // expression contraction can choose a different source layout.
                let cooked_block = cooking.try_cook(block_d.as_view())?;
                assert_eq!(cooking.uncook(cooked_block.as_view()), block_d);
                let block = SymbolicTensor::infer_indexed(cooked_block)?;
                let normalized_block = block.contract(attachment)?;
                assert_eq!(
                    block.structure().canonical().external_structure(),
                    normalized_block
                        .structure()
                        .canonical()
                        .external_structure()
                );
                let normalized_block_source =
                    cooking.uncook(normalized_block.expression().as_view());
                let projected_native = normalized_block.multiply(renamed_native)?;
                let projected_native = projected_native
                    .with_checked_expression(Atom::num(sign) * projected_native.expression())?;
                let template_native = evaluated_native.try_add(&projected_native)?;
                let template_native = template_native
                    .with_checked_expression(template_native.expression() / Atom::num(2))?;
                assert_eq!(
                    branch_native.structure().logical_slots(),
                    template_native.structure().logical_slots()
                );
                let template = cooking.uncook(template_native.expression().as_view());
                let literal_attachment = branch_native.same_value(&template_native);
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_actual_projector_template", label, sign, literal_attachment,
                    graph = %graph.name, child = %child.subgraph().string_label(),
                    file.actual_branch = %branch.to_canonical_string(),
                    file.actual_source = %evaluated.to_canonical_string(),
                    file.external_block_source = %block_d.to_canonical_string(),
                    file.external_block = %normalized_block_source.to_canonical_string(),
                    file.template = %template.to_canonical_string(),
                    "Certifying actual ordered projector factor and evaluated source separately"
                );
                // Identical operands already give an exact attachment proof.
                // Still admit the tensor and its interface, without repeating
                // expensive normalization of the complete source. Otherwise
                // contract only spinor metrics and canonicalize internal indices.
                let mut attached = Vec::new();
                for value in [branch_native, &template_native]
                    .into_iter()
                    .take(if literal_attachment { 1 } else { 2 })
                {
                    let source = cooking.uncook(value.expression().as_view());
                    let cooked = cooking.try_cook(source.as_view())?;
                    assert_eq!(cooking.uncook(cooked.as_view()), source);
                    let tensor = value.with_checked_expression(value.expression().clone())?;
                    let logical = tensor.structure().logical_slots();
                    let slots = tensor.structure().canonical().external_structure();
                    let normalized = if literal_attachment {
                        tensor
                    } else {
                        tensor.contract(attachment)?.canonize()?
                    };
                    assert_eq!(
                        slots,
                        normalized.structure().canonical().external_structure()
                    );
                    assert_eq!(logical, normalized.structure().logical_slots());
                    attached.push(normalized);
                }
                if !literal_attachment {
                    assert_eq!(
                        attached[0].structure().canonical().external_structure(),
                        attached[1].structure().canonical().external_structure()
                    );
                    eyre::ensure!(
                        attached[0].expression() == attached[1].expression(),
                        "actual branch differs from the native ordered external-block template"
                    );
                }
                let actual_source = at_point(&evaluated_native, sign)?;
                let expected_source = at_point(input_native, sign)?;
                exact_zero(
                    "actual evaluated X_sigma equals independently routed source",
                    &actual_source.try_sub(&expected_source)?,
                )?;
            }
            assert_eq!(
                *input, untouched,
                "formal-d graph source stays byte-exact during external-block proof"
            );
            assert!(input_native.same_value(&untouched_native));
            Ok(())
        };
        certify_branches("complete parent input", &source, [&positive, &negative])?;

        // Compare physical rational terms in one certified chart. Owner
        // metadata is retained below as a separate witness, not mistaken
        // for a scalar denominator. Only the existing denominator owner
        // expands its complete small polynomial; numerators retain d and
        // their factorized source circuits.
        let common_chart_zero = |label: &str,
                                 left: &SymbolicTensor<IndexedStructure>,
                                 right: &SymbolicTensor<IndexedStructure>|
         -> Result<()> {
            let mut coefficients = BTreeMap::<Vec<Atom>, SymbolicTensor<IndexedStructure>>::new();
            let mut witnesses = Vec::new();
            for (side, sign, expression) in [("left", 1, left), ("right", -1, right)] {
                let raw_terms = FourDTerm::from_tensor(expression)?;
                for (term, value) in raw_terms.iter().enumerate() {
                    for denominator in &value.denominators {
                        eyre::ensure!(
                            denominator.source_edge.0 < graph.underlying.n_edges()
                                && parent
                                    .subgraph()
                                    .includes(&graph[&denominator.source_edge].1),
                            "{label}: denominator lost its native parent source owner"
                        );
                        witnesses.push(serde_json::json!({
                            "side":side,"term":term,"source_edge":denominator.source_edge.0,
                            "momentum":denominator.momentum.to_canonical_string(),
                            "mass_squared":denominator.mass_squared.to_canonical_string(),
                            "full_polynomial":denominator.full_expr.to_canonical_string(),
                        }));
                    }
                }
                // Admit compact source contractions before routing their
                // vector arguments. Then delta_i(Q) is an explicit dot,
                // so Q -> Q4-Q2 stays a linear vector contraction rather
                // than putting that sum into rank-one scalar metadata.
                let admitted = expression;
                let dotted = admitted.to_dots()?;
                assert_eq!(admitted.structure(), dotted.structure());
                let physical_expression = GS.resolve_os_momenta(
                    &GS.erase_uv_momentum_provenance(dotted.expression())
                        .replace_multiple(graph.uv_wrapped_replacement(
                            parent.subgraph(),
                            &lmb,
                            &[W_.x___],
                        )),
                );
                let physical = dotted.with_checked_expression(physical_expression)?;
                for term in FourDTerm::from_tensor(&physical)? {
                    let mut polynomials = term
                        .denominators
                        .iter()
                        .map(|denominator| denominator.polynomial_in_lmb(&lmb))
                        .collect::<Result<Vec<_>>>()?;
                    polynomials.sort(); // repeated entries retain propagator powers
                    let numerator = term
                        .numerator
                        .with_checked_expression(Atom::num(sign) * term.numerator.expression())?;
                    if let Some(sum) = coefficients.get_mut(&polynomials) {
                        *sum = sum.try_add(&numerator)?;
                    } else {
                        coefficients.insert(polynomials, numerator);
                    }
                }
            }
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_parent_common_chart", label,
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.left = %cooking.uncook(left.expression().as_view()).to_canonical_string(), file.right = %cooking.uncook(right.expression().as_view()).to_canonical_string(),
                file.owners = %serde_json::to_string(&witnesses)?,
                file.lmb = ?lmb, denominator_classes = coefficients.len(),
                "Preserved native owners before formal-d denominator-class comparison"
            );
            for (class, (denominators, numerator)) in coefficients.into_iter().enumerate() {
                debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_vertex_parent_denominator_class", label, class,
                    file.denominators = ?denominators,
                    file.numerator = %cooking.uncook(numerator.expression().as_view()).to_canonical_string(),
                    "Only identical full denominator polynomials share a numerator certificate"
                );
                frame_zero(&format!("{label}, denominator class {class}"), &numerator)?;
            }
            Ok(())
        };
        let ordinary = t_raw(&source, ctx, parent, child, &lmb)?;
        if child.renormalization_scheme() != ApproximationType::OS {
            eyre::ensure!(
                !ordinary.expression().contains_symbol(AGS.gamma0)
                    && !ordinary.expression().contains_symbol(ETS.delta),
                "ordinary vacuum Taylor grading cannot introduce a fixed frame"
            );
        }
        let on_shell = positive.try_add(&negative)?;
        let on_shell_resolved =
            on_shell.with_checked_expression(GS.resolve_os_momenta(on_shell.expression()))?;
        let overlap = t_raw(&on_shell_resolved, ctx, parent, child, &lmb)?;
        let [ordinary_positive, ordinary_negative] =
            graph.on_shell_branches(&graph.model, &ordinary, parent.subgraph(), &lmb, &velocity)?;
        certify_branches(
            "ordinary parent coefficient",
            &ordinary,
            [&ordinary_positive, &ordinary_negative],
        )?;
        common_chart_zero(
            "actual open-tensor U1 O equals O U1",
            &overlap,
            &ordinary_positive.try_add(&ordinary_negative)?,
        )?;
        let actual = uv_limit(prefix, ctx, parent, child, parent, child)?;
        common_chart_zero(
            "actual parent U1/O/U1O composition",
            actual.atom(),
            &overlap.try_sub(&ordinary)?.try_sub(&on_shell)?,
        )?;
        let mut vacuum_sectors = actual.active_sectors().iter().map(|s| s.integration_atom());
        let mut vacuum = vacuum_sectors
            .next()
            .ok_or_else(|| eyre!("actual parent has no analytic sector"))?
            .clone();
        for sector in vacuum_sectors {
            vacuum = vacuum.try_add(sector)?;
        }
        common_chart_zero(
            "actual analytic input is only the signed vacuum remainder",
            &vacuum,
            &overlap.try_sub(&ordinary)?,
        )?;
        // The certified finite Phi P block, actual branch templates,
        // source gates, and U1/O commutation imply Phi[-(U1+O-U1O)X]
        // = -Phi[X], without evaluating mixed-temporal formal-d traces.
        assert!(actual.recursive_completion().is_empty());
        assert_eq!(actual.active_sectors().len(), prefix.sectors().count());
        debug_tags!(#uv, #integrated, #inspect;
            stage = "native_mass_parent_input_certificate", source = %diagram.name(),
            given = %child.subgraph().string_label(), given_scheme = %child.renormalization_scheme(),
            file.velocity = ?velocity, file.parent_lmb = ?lmb,
            actual_parent_compared = true, actual_analytic_input_compared = true,
            integrated_parent_finite_computed = false, public_dispatch_admitted = false,
            "Actual parent/projector/source identities; actual finite integration is a separate next gate"
        );
        Ok(actual)
    }

    fn certify_parent_mass_projection(&self, subtraction_velocity: [&str; 4]) -> Result<()> {
        let reference = self;
        let mut certified = reference.certify_integrated_qg_children(&[6, 6, 21])?;
        certified.extend(reference.certify_integrated_qg_children(&[6, 21, 21])?);
        assert_eq!(certified.len(), 4);
        let exact_zero = Self::certify_parent_tensor_zero;
        for record in certified {
            let graph = &record.graph;
            let child = &record.child;
            let mut settings = record.settings.clone();
            let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
            let mut parents = graph
                .classified_spinneys(&uncut, &graph.model, &settings, &graph.loop_momentum_basis)?
                .into_iter()
                .filter(|spinney| {
                    spinney.n_components() == 1
                        && spinney.dod == 1
                        && spinney.renormalization_scheme == ApproximationType::OS
                        && graph.n_loops(spinney.filter()) == 2
                        && spinney.filter().includes(child.subgraph())
                })
                .collect::<Vec<_>>();
            assert_eq!(parents.len(), 1);
            let parent = OwnedForestNode {
                spinney: parents.remove(0),
                topo_order: 2,
            };
            let guard_error = validate_local_counterterm_schemes(
                graph,
                &super::super::CutStructure::empty(graph),
                &settings,
            )
            .expect_err("this private proof must not admit the public OS/MUV forest");
            let guard_error = format!("{guard_error:#}");
            assert!(
                guard_error.contains("on-shell/OS component")
                    && guard_error.contains("unsupported nested or disconnected scheme policy")
            );
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_prefix_parent_proof_start",
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.child_edges = ?graph.paired_edges(child.subgraph()),
                file.child_lmb = ?child.lmb(), public_guard_error = %guard_error,
                "Inspecting actual certified child before private parent composition"
            );

            let prefix = Full4dCts::recursion_input(
                &record.local,
                &record.integrated,
                ApproximationType::MUV,
                false,
                child.lmb(),
            )?;
            let sectors = prefix.sectors().collect::<Vec<_>>();
            assert_eq!(
                sectors.len(),
                2,
                "one local child plus one actual finite completion"
            );
            let active = sectors
                .iter()
                .filter(|s| s.frozen_lmbs().is_empty())
                .collect::<Vec<_>>();
            let frozen = sectors
                .iter()
                .filter(|s| !s.frozen_lmbs().is_empty())
                .collect::<Vec<_>>();
            assert_eq!(active.len(), 1);
            assert_eq!(frozen.len(), 1);
            assert_eq!(active[0].active_components.len(), 1);
            assert_eq!(&active[0].active_components[0].0, child.subgraph());
            assert_eq!(&active[0].active_components[0].1, child.subgraph());
            assert_eq!(active[0].active_components[0].2.loop_edges.len(), 1);
            assert!(frozen[0].active_components.is_empty());
            assert_eq!(frozen[0].frozen_lmbs(), std::slice::from_ref(child.lmb()));
            let physical_finite = record.integrated.physical_finite_counterterm()?;
            let bookkeeping_finite = physical_finite.with_checked_expression(
                physical_finite.expression() * Atom::var(GS.integrated_loop_scale).pow(4),
            )?;
            exact_zero(
                "consumed one-loop measure marker",
                &frozen[0].atom.try_sub(&bookkeeping_finite)?,
            )?;
            let child_integral = record
                .integrated
                .physical_pole()?
                .try_sub(&physical_finite)?;
            let frozen_physical = frozen[0].atom.with_checked_expression(
                frozen[0]
                    .atom
                    .expression()
                    .replace(GS.integrated_loop_scale)
                    .with(1),
            )?;
            let completed_integral = child_integral.try_add(&frozen_physical)?;
            let completed_integral = completed_integral.with_checked_expression(
                graph
                    .model
                    .expand_couplings(completed_integral.expression())
                    .replace(GS.m_uv_vacuum)
                    .with(1)
                    .replace(GS.m_uv_expansion)
                    .with(1)
                    .replace(GS.mu_r_sq)
                    .with(1),
            )?;
            let pole_difference = completed_integral.try_sub(&record.pole_reference)?;
            let pole_difference = pole_difference.with_checked_expression(
                pole_difference
                    .expression()
                    .series(GS.dim_epsilon, 0, 1)?
                    .to_atom(),
            )?;
            exact_zero(
                "completed native vertex equals independent pole tensor",
                &pole_difference,
            )?;

            // Use actual component/cograph ownership; no synthetic replacement of
            // the completed local child enters the parent operation.
            let reduced = parent.reduced_subgraph(child);
            let cograph_edges = graph.paired_edges(&reduced);
            assert_eq!(cograph_edges.len(), 2);
            let mut cograph_pdgs = cograph_edges
                .iter()
                .map(|edge| {
                    graph
                        .model
                        .particle_by_id(graph[*edge].particle.particle().unwrap())
                        .map(|particle| particle.pdg_code.abs())
                })
                .collect::<Result<Vec<_>, _>>()?;
            cograph_pdgs.sort_unstable();
            assert_eq!(cograph_pdgs, [6, 21]);
            let child_frame = &active[0].active_components[0].2;
            for edge in &cograph_edges {
                assert!(
                    child_frame.edge_signatures[*edge]
                        .internal
                        .iter()
                        .all(|sign| sign.is_zero()),
                    "cograph edge must be independent of the actual child hard integration coordinate"
                );
            }
            let cograph_numerator = graph
                .numerator(&reduced, child.subgraph())?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            let cograph_denominator = graph.denominator(&reduced, &graph.model, |_| 1);
            let denominator_terms = FourDTerm::from_view(cograph_denominator.pow(-1).as_view())?;
            assert_eq!(denominator_terms.len(), 1);
            assert_eq!(denominator_terms[0].denominators.len(), 2);
            assert_eq!(
                denominator_terms[0]
                    .denominators
                    .iter()
                    .map(|d| d.source_edge)
                    .collect::<BTreeSet<_>>(),
                cograph_edges.iter().copied().collect()
            );

            // The child boundary gluon becomes an INTERNAL summed parent slot.
            // Check its native formal representation before the parent finalizer.
            let finite_tensor = &frozen[0].atom;
            let mut slots = SlotMatcher::default();
            let mut child_lorentz_ports = 0;
            for slot in finite_tensor.structure().logical_slots() {
                let atom = slot.to_atom();
                if slots
                    .parse::<Minkowski, AbstractIndex>(atom.as_view())
                    .is_ok()
                {
                    let AtomView::Fun(port) = atom.as_view() else {
                        unreachable!()
                    };
                    assert_eq!(port.get(0), Atom::var(GS.dim).as_view());
                    child_lorentz_ports += 1;
                }
            }
            assert_eq!(child_lorentz_ports, 1);
            let attached = prefix.atom().clone().multiply(cograph_numerator)?;
            let attached =
                attached.with_checked_expression(attached.expression() / &cograph_denominator)?;
            for slot in attached.structure().logical_slots() {
                let atom = slot.to_atom();
                assert!(
                    slots
                        .parse::<Minkowski, AbstractIndex>(atom.as_view())
                        .is_err(),
                    "the actual child gluon port must close against its cograph partner"
                );
            }

            settings.on_shell_subtraction_velocity = subtraction_velocity.map(str::to_owned);
            let ctx = UVCtx::new(graph, &graph.model, &settings);
            let lmb = compatible_lmb(&ctx, &parent, child, &prefix)?;
            // Audit the actual hard grading before extracting the parent's U1.
            // No frozen localization kernel is part of this coefficient.
            let t = Atom::var(GS.rescale);
            let scaled = frozen[0].atom.with_checked_expression(graph.uv_rescaled(
                parent.subgraph(),
                2,
                &lmb,
                &lmb,
                frozen[0].atom.expression(),
            ))?;
            let expected_scaled = frozen[0].atom.with_checked_expression(
                frozen[0]
                    .atom
                    .expression()
                    .replace(GS.m_uv_vacuum)
                    .with(Atom::var(GS.m_uv_vacuum) / &t)
                    .replace(GS.integrated_loop_scale)
                    .with(Atom::var(GS.integrated_loop_scale) * &t)
                    * t.pow(-8),
            )?;
            exact_zero(
                "completed-child mass and consumed-measure hard grading",
                &scaled.try_sub(&expected_scaled)?,
            )?;
            let actual = self.certify_parent_mass_projection_input(
                &record.diagram,
                &ctx,
                &parent,
                child,
                &prefix,
            )?;
            let velocity = settings.on_shell_velocity()?;
            assert_eq!(actual.active_sectors().len(), 2);
            assert!(actual.recursive_completion().is_empty());
            assert_eq!(
                actual
                    .active_sectors()
                    .iter()
                    .filter(|s| !s.frozen_lmbs().is_empty())
                    .count(),
                1
            );
            let output_frozen = actual
                .active_sectors()
                .iter()
                .find(|sector| !sector.frozen_lmbs().is_empty())
                .expect("the actual parent retains its one frozen completion");
            assert_eq!(
                output_frozen.frozen_lmbs(),
                std::slice::from_ref(child.lmb()),
                "the actual parent output retains the complete child localization basis"
            );
            debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_prefix_parent_mass_projection",
                graph = %graph.name, child = %child.subgraph().string_label(),
                file.velocity = ?velocity, file.parent_lmb = ?lmb,
                active_sectors = actual.active_sectors().len(), frozen_sectors = 1,
                public_dispatch_admitted = false,
                "Actual native child prefix satisfies matched parent mass-projection identities"
            );
        }
        Ok(())
    }
}
