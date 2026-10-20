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

impl QcdMassReference {
    /// Open a complete source-derived quark closing word on its native frame.
    /// Keep momentum specialization, normalization and Gram controls at callers.
    pub(in crate::uv) fn open_native_quark_dual(
        word: &SymbolicTensor<IndexedStructure>,
        expected_open: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        use spenso::structure::slot::IsAbstractSlot;
        use std::collections::HashMap;

        let slots = word.structure().logical_slots();
        ensure!(slots.len() == 4 && expected_open.rank() == 4);
        let mut map = HashMap::new();
        let mut inverse = HashMap::new();
        for (position, slot) in slots.iter().enumerate() {
            let opposite = slots
                .iter()
                .filter(|other| slot.rep().matches(&other.rep()) && other.aind() != slot.aind())
                .collect::<Vec<_>>();
            ensure!(
                opposite.len() == 1,
                "native quark dual needs one named opposite spin/color port"
            );
            map.insert(position, opposite[0].aind());
            inverse.insert(position, slot.aind());
        }
        let open = word.reindex_interface_ports(&map)?;
        ensure!(open.reindex_interface_ports(&inverse)?.same_value(word));
        ensure!(
            open.structure().canonical().external_structure()
                == expected_open.structure().canonical().external_structure()
        );
        let zero = word.with_checked_expression(Atom::Zero)?;
        let moved_zero = zero.reindex_interface_ports(&map)?;
        ensure!(moved_zero.rank() == 4);
        ensure!(
            moved_zero
                .reindex_interface_ports(&inverse)?
                .same_value(&zero)
        );
        Ok(open)
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(in crate::uv) enum D09VacuumSource {
    Bare,
    ProperLocal,
    ProperFinite,
}

type CertifiedMassPole = (usize, Atom, Atom, Atom);

/// Source-owned unit-rest frame retained with the actual vacuum trace rows.
/// Only these certified boundary momenta are specialized; active loop carriers
/// and the surrounding graph remain untouched.
pub(in crate::uv) struct D09VacuumFrame {
    pub(in crate::uv) closing_duals: [SymbolicTensor<IndexedStructure>; 2],
    pub(in crate::uv) open_basis: [SymbolicTensor<IndexedStructure>; 2],
    pub(in crate::uv) external_points: BTreeMap<linnet::half_edge::involution::EdgeIndex, Atom>,
    pub(in crate::uv) mass: Atom,
    pub(in crate::uv) lmb: crate::graph::LoopMomentumBasis,
}

impl D09VacuumFrame {
    pub(in crate::uv) fn at_point(
        &self,
        input: &SymbolicTensor<IndexedStructure>,
        graph: &crate::graph::Graph,
        parent: &crate::uv::hedge_poset::OwnedForestNode,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        use crate::{
            graph::lmb::LMBext,
            utils::{GS, W_},
            uv::approx::ForestNodeLike,
        };
        use linnet::half_edge::involution::EdgeIndex;
        let routed = GS
            .erase_uv_momentum_provenance(input.expression())
            .replace_multiple(graph.uv_wrapped_replacement(
                parent.subgraph(),
                &self.lmb,
                &[W_.x___],
            ));
        let mut unexpected = false;
        let expression = routed.replace_map(|view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == GS.emr_mom
            {
                if momentum.get_nargs() == 0 {
                    unexpected = true;
                    return;
                }
                let Ok(edge) = usize::try_from(momentum.get(0)) else {
                    unexpected = true;
                    return;
                };
                let edge = EdgeIndex(edge);
                if let Some(point) = self.external_points.get(&edge) {
                    let indices = momentum
                        .iter()
                        .skip(1)
                        .map(|index| index.to_owned())
                        .collect::<Vec<_>>();
                    **output = GS.indexed_momentum(point, &indices);
                } else if ![EdgeIndex(2), EdgeIndex(4)].contains(&edge) {
                    unexpected = true;
                }
            }
        });
        ensure!(!unexpected, "D09 fixed-point view has an unowned momentum");
        Ok(input.with_checked_expression(
            GS.resolve_os_momenta(&expression)
                .replace(self.mass.to_pattern())
                .with(1),
        )?)
    }

    /// Finite external matrix only: the Laurent owner acts before D/epsilon
    /// specialization, and scalar radial coefficients stay outside components.
    pub(in crate::uv) fn unit_finite(
        &self,
        input: &SymbolicTensor<IndexedStructure>,
        graph: &crate::graph::Graph,
        parent: &crate::uv::hedge_poset::OwnedForestNode,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        use crate::utils::GS;
        let point = self.at_point(input, graph, parent)?;
        ensure!(point.rank() == 4);
        for head in [
            GS.den,
            GS.emr_mom,
            GS.loop_mom,
            GS.os_mom,
            GS.integrated_loop_scale,
        ] {
            ensure!(!point.expression().contains_symbol(head));
        }
        let normalized = point.simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            ..Default::default()
        })?;
        let unit = normalized
            .expression()
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum)
            .replace(GS.m_uv_vacuum)
            .with(1)
            .replace(GS.mu_r_sq)
            .with(1);
        let expansion = super::series(&super::Integrated::dimensionally_regularized(&unit), 0)?;
        let finite = normalized.with_checked_expression(super::truncate(&expansion, true))?;
        // The remaining Lorentz block must be supported by fixed physical
        // selectors. Never turn an unresolved internal D trace into Minkowski4.
        use crate::numerator::aind::Aind;
        use spenso::structure::slot::SlotMatcher;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let raw = cooking.uncook(finite.expression().as_view());
        let mut matcher = SlotMatcher::default();
        let mut lorentz = BTreeSet::new();
        let mut fixed = BTreeSet::new();
        let mut invalid = false;
        raw.visitor(&mut |view| {
            if let AtomView::Fun(selector) = view
                && selector.get_symbol() == GS.delta_vec
                && selector.get_nargs() == 2
                && let Some((representation, _)) = matcher.port_representation(selector.get(1))
                && representation
                    == spenso::structure::representation::LibraryRep::from(Minkowski {})
            {
                invalid |= !matches!(matcher.concrete_component(selector.get(0)), Some(Ok(0..=3)));
                // A compact representation is an implicit port local to this
                // fixed selector. It cannot cover a different compact occurrence.
                if matcher
                    .parse_representation::<Minkowski>(selector.get(1))
                    .is_ok()
                {
                    return false;
                }
                fixed.insert(selector.get(1).to_owned());
            }
            if matcher.parse::<Minkowski, Aind>(view).is_ok() {
                lorentz.insert(view.to_owned());
                return false;
            }
            if matcher.parse_representation::<Minkowski>(view).is_ok() {
                invalid = true;
                return false;
            }
            true
        });
        if invalid || lorentz != fixed {
            crate::debug_tags!(#uv, #integrated, #inspect;
                stage = "native_d09_finite_frame_direction_rejected", invalid,
                file.lorentz = ?lorentz, file.fixed = ?fixed,
                file.input = %raw.to_canonical_string(),
                "Finite physical admission retains the exact unresolved direction"
            );
        }
        ensure!(
            !invalid && lorentz == fixed,
            "unit-rest finite carrier has an unresolved Lorentz direction"
        );
        let finite = finite.map_representation_dimensions(|representation, dimension| {
            if representation == spenso::structure::representation::LibraryRep::from(Minkowski {}) {
                dimension.replace(GS.dim).with(4)
            } else {
                dimension.to_owned()
            }
        })?;
        QcdMassReference::materialize_finite_external_matrix_coefficients(&finite)
    }
}

struct D09VacuumScalarRows {
    frame: D09VacuumFrame,
    dimensional: SymbolicTensor<IndexedStructure>,
    physical_four: SymbolicTensor<IndexedStructure>,
    production_finite: SymbolicTensor<IndexedStructure>,
    weighted_difference: [Atom; 2],
    source_weight: Atom,
    external_ordering: Atom,
    omega: Atom,
}

// Actual source/child incidence retains the certified scalar transfer and the
// same native input consumed by the complete inventory and overall operators.
type CertifiedMassPrefixes = BTreeMap<(String, Vec<usize>), CertifiedMassPrefix>;

// One actual source supplies the common, fully named open field interface.
struct NativeFieldFrame {
    kinetic: SymbolicTensor<IndexedStructure>,
    identity: SymbolicTensor<IndexedStructure>,
    mass: Atom,
    endpoint_slots: [Atom; 2],
    incoming_edge: usize,
}

struct CertifiedMassPrefix {
    coefficient: Atom,
    diagram: FeynmanDiagram,
    graph: crate::graph::Graph,
    given: crate::uv::hedge_poset::OwnedForestNode,
    completed: crate::uv::approx::local_4d::Full4dCts,
    settings: crate::uv::UVgenerationSettings,
}

/// Completed finite projection; raw_phi retains internal MT for the derivative oracle.
struct MassProjection {
    raw_phi: Atom,
    scalar: Atom,
    family: IntegralFamily,
    colors: Vec<(Atom, Atom)>,
    // Edge identity, routed momentum, squared mass, full unit-mass polynomial.
    poles: Vec<CertifiedMassPole>,
    kernel_weight: Atom,
    // Existing source-derived amplitude/kernel conversion; never fitted to CTs.
    external_ordering: Atom,
    numerator_prefactor: Atom,
    coupling_definitions: Vec<(String, Atom)>,
    scaleless: bool,
}

#[test]
fn native_fixed_component_rational_identity_preserves_pole_owners_and_weights() -> Result<()> {
    use crate::utils::GS;
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    crate::initialisation::test_initialise()?;
    let row = Atom::var(symbol!("native_fixed_pole_regression::row"));
    let column = Atom::var(symbol!("native_fixed_pole_regression::column"));
    let spin_row = idenso::bis!(4, row.clone());
    let spin_column = idenso::bis!(4, column.clone());
    let color_row = ColorFundamental {}.new_rep(3).to_symbolic([row]);
    let color_column = function!(
        AIND_SYMBOLS.dind,
        ColorFundamental {}.new_rep(3).to_symbolic([column])
    );
    let identity = function!(ETS.metric, &spin_row, &spin_column)
        * function!(ETS.metric, &color_row, &color_column);
    let selector =
        |n: i32, port: &Atom| function!(ETS.delta, function!(AIND_SYMBOLS.cind, n), port);
    let independent = selector(2, &spin_row)
        * selector(1, &spin_column)
        * selector(0, &color_row)
        * selector(1, &color_column);
    let [a, b, u, v, s, t] = ["a", "b", "u", "v", "s", "t"]
        .map(|name| Atom::var(symbol!(&format!("native_fixed_pole_regression::{name}"); Scalar)));
    let epsilon = Atom::var(GS.dim_epsilon);
    let logarithm = function!(
        symbolica::atom::Symbol::LOG,
        Atom::var(symbol!("native_fixed_pole_regression::scale"; Scalar))
    );
    let jet = (Atom::one() + &epsilon * &logarithm)
        / ((Atom::num(2) - &epsilon) * (Atom::num(3) - &epsilon));
    let integrated_loop_scale = Atom::var(GS.integrated_loop_scale);
    let factors = [
        (&u + &v).pow(11),
        (&s + &t).pow(-2),
        integrated_loop_scale.pow(4),
    ];
    let weight = &jet * Atom::mul_many(factors.iter());
    let q0 = function!(GS.emr_mom, 2, function!(AIND_SYMBOLS.cind, 0));
    let q1 = function!(GS.emr_mom, 4, function!(AIND_SYMBOLS.cind, 1));
    let rep = (Minkowski {}).to_symbolic([GS.dim]);
    let momentum = |edge: i32| function!(GS.emr_mom, edge);
    let square = |edge: i32| {
        function!(
            SPENSO_TAG.dot,
            function!(GS.emr_mom, edge, &rep),
            function!(GS.emr_mom, edge, &rep)
        )
    };
    let mass = Atom::var(GS.m_uv_expansion).pow(2);
    let d2 = square(2) - &mass;
    let d4 = square(4) - &mass;
    let pole = |edge: i32, momentum: &Atom, mass: &Atom, polynomial: &Atom| {
        function!(GS.den, edge, momentum, mass, polynomial)
    };
    let p2 = pole(2, &momentum(2), &mass, &d2);
    let p4 = pole(4, &momentum(4), &mass, &d4);
    let numerator = &identity * q0.pow(2) + &independent * &q1;
    let lhs = &weight * (&a * &numerator / (p2.pow(2) * &p4) + &b * &numerator / (&p2 * p4.pow(2)));
    let rhs = &weight
        * ((&a * &identity * q0.pow(2) + &a * &independent * &q1) / (p2.pow(2) * &p4)
            + (&b * &identity * q0.pow(2) + &b * &independent * &q1) / (&p2 * p4.pow(2)));
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let left = SymbolicTensor::infer_indexed(cooking.try_cook(lhs.as_view())?)?;
    let right = left.with_checked_expression(cooking.try_cook(rhs.as_view())?)?;
    let before = left.clone();
    QcdMassReference::certify_fixed_component_rational_identity(
        "independent fixed-component pole-class distribution",
        &left,
        &right,
    )?;
    for factor in &factors {
        ensure!(right.expression().contains(factor.as_view()));
    }
    let wrong_weight = right.with_checked_expression(
        right.expression()
            + cooking.try_cook((&weight * &independent / (&p2 * p4.pow(2))).as_view())?,
    )?;
    ensure!(
        QcdMassReference::certify_fixed_component_rational_identity(
            "unit native matrix-weight perturbation",
            &left,
            &wrong_weight
        )
        .is_err()
    );
    let wrong_scale = rhs
        .replace(integrated_loop_scale.to_pattern())
        .with(integrated_loop_scale.pow(2));
    let wrong_scale = right.with_checked_expression(cooking.try_cook(wrong_scale.as_view())?)?;
    ensure!(
        QcdMassReference::certify_fixed_component_rational_identity(
            "changed consumed-child loop-scale power",
            &left,
            &wrong_scale
        )
        .is_err()
    );
    let changed = [
        pole(3, &momentum(2), &mass, &d2),
        pole(2, &momentum(4), &mass, &d2),
        pole(2, &momentum(2), &(&mass + Atom::one()), &d2),
        pole(
            2,
            &momentum(2),
            &mass,
            &(&d2
                + Atom::i()
                    * Atom::var(symbol!("native_fixed_pole_regression::prescription"; Scalar))),
        ),
        p2.pow(2),
    ];
    for altered in changed {
        let altered = rhs.replace(p2.to_pattern()).with(altered);
        let altered = right.with_checked_expression(cooking.try_cook(altered.as_view())?)?;
        ensure!(
            QcdMassReference::certify_fixed_component_rational_identity(
                "changed original causal pole descriptor or power",
                &left,
                &altered
            )
            .is_err()
        );
    }
    ensure!(left.same_value(&before));
    Ok(())
}

#[test]
fn native_complete_tensor_recomposition_preserves_factored_matrix_radial_distribution() -> Result<()>
{
    use crate::utils::GS;
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    crate::initialisation::test_initialise()?;
    let row = Atom::var(symbol!("native_recomposition_distribution::row"));
    let column = Atom::var(symbol!("native_recomposition_distribution::column"));
    let spin_row = idenso::bis!(4, row.clone());
    let spin_column = idenso::bis!(4, column.clone());
    let color_row = ColorFundamental {}.new_rep(3).to_symbolic([row]);
    let color_column = function!(
        AIND_SYMBOLS.dind,
        ColorFundamental {}.new_rep(3).to_symbolic([column])
    );
    let identity = function!(ETS.metric, &spin_row, &spin_column)
        * function!(ETS.metric, &color_row, &color_column);
    let coordinate =
        |index: i32, port: &Atom| function!(ETS.delta, function!(AIND_SYMBOLS.cind, index), port);
    let independent = coordinate(2, &spin_row)
        * coordinate(1, &spin_column)
        * coordinate(0, &color_row)
        * coordinate(1, &color_column);
    let rep = (Minkowski {}).to_symbolic([GS.dim]);
    let radial = function!(
        SPENSO_TAG.dot,
        function!(GS.loop_mom, 1, &rep),
        function!(GS.loop_mom, 2, &rep)
    );
    let [a, b, u, v, s, t] = ["a", "b", "u", "v", "s", "t"].map(|name| {
        Atom::var(symbol!(&format!("native_recomposition_distribution::{name}"); Scalar))
    });
    let epsilon = Atom::var(GS.dim_epsilon);
    let logarithm = function!(
        symbolica::atom::Symbol::LOG,
        Atom::var(symbol!("native_recomposition_distribution::scale"; Scalar))
    );
    let jet = (Atom::one() + &epsilon * &logarithm)
        / ((Atom::num(2) - &epsilon) * (Atom::num(3) - &epsilon));
    let spectator_factors = [(&u + &v).pow(11), (&s + &t).pow(-2)];
    let spectator = Atom::mul_many(spectator_factors.iter());
    let matrix_distribution = &a * (&identity * &radial + &independent * &radial)
        + &b * (&identity * &radial - &independent * &radial);
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let input = SymbolicTensor::infer_indexed(
        cooking.try_cook((&spectator * &jet * matrix_distribution).as_view())?,
    )?;
    let original = input.clone();
    let selection = spenso::shadowing::TensorCollectFilter::Reps([
        spenso::structure::representation::LibraryRep::from(Minkowski {}),
    ]);
    let partial: SymbolicTensor<PartialStructure> = input.clone().into();
    let rows = partial.coefficient_list(selection)?;
    ensure!(
        rows.iter()
            .all(|(radial, matrix)| radial.is_scalar() && matrix.rank() == 4)
    );
    let recomposed = input.with_checked_expression(Atom::add_many(
        rows.iter()
            .map(|(radial, matrix)| radial.expression() * matrix.expression()),
    ))?;
    ensure!(
        input.coefficients_equal(
            &recomposed,
            spenso::shadowing::TensorCollectFilter::<0>::Tensors
        )? == symbolica::id::ConditionResult::True
    );
    for factor in &spectator_factors {
        ensure!(recomposed.expression().contains(factor.as_view()));
    }
    ensure!(recomposed.expression().contains(radial.as_view()));
    ensure!(recomposed.structure().logical_slots() == input.structure().logical_slots());
    let wrong = recomposed.with_checked_expression(
        recomposed.expression()
            + cooking.try_cook((&spectator * &jet * &radial * &independent).as_view())?,
    )?;
    ensure!(
        input.coefficients_equal(&wrong, spenso::shadowing::TensorCollectFilter::<0>::Tensors)?
            != symbolica::id::ConditionResult::True
    );
    ensure!(input.same_value(&original));
    Ok(())
}

#[test]
fn native_finite_scalar_certificate_preserves_captured_d09_laurent_recomposition() -> Result<()> {
    use crate::numerator::AppliedFeynmanRule;
    use spenso::network::library::symbolic::ETS;

    crate::initialisation::test_initialise()?;
    let scalar = Atom::parse(
        include_str!("../../../tests/resources/uv_fixed_basis/d09_open_scalar_recomposition.sym"),
        "gammalooprs",
        Default::default(),
    )
    .map_err(|error| eyre!(error))?;
    ensure!(
        !scalar.is_zero(),
        "the captured factorized Laurent row must reach the finite certificate"
    );
    let metric = function!(
        ETS.metric,
        idenso::bis!(
            4,
            Atom::var(symbol!("native_finite_laurent_regression::row"))
        ),
        idenso::bis!(
            4,
            Atom::var(symbol!("native_finite_laurent_regression::column"))
        )
    );
    let word = AppliedFeynmanRule::new(metric.clone())?.expr;
    ensure!(word.rank() == 2);
    let trace = super::simplify(&word.clone().multiply(word)?)?;
    ensure!(trace.is_scalar());
    ensure!((trace.expression() - Atom::num(4)).is_zero());

    let captured = AppliedFeynmanRule::new(&scalar * &metric)?.expr;
    ensure!(captured.rank() == 2);
    QcdMassReference::certify_finite_tagged_rows_zero(
        "captured D09 open scalar Laurent recomposition on a nonzero native metric",
        captured,
    )?;
    let shifted = AppliedFeynmanRule::new((&scalar + Atom::one()) * &metric)?.expr;
    ensure!(
        QcdMassReference::certify_finite_tagged_rows_zero(
            "unit perturbation of the captured D09 finite scalar",
            shifted,
        )
        .is_err(),
        "the finite certificate must reject the same metric with a unit scalar shift"
    );
    Ok(())
}

impl QcdMassReference {
    /// Certify only the finite scalar rows of a completed, preprocessed tensor.
    /// Callers retain their distinct Dirac/color preprocessing contracts.
    pub(in crate::uv) fn certify_finite_tagged_rows_zero(
        label: &str,
        tensor: SymbolicTensor<IndexedStructure>,
    ) -> Result<()> {
        use crate::{numerator::symbolica_ext::ScalarRationalSize, utils::GS};

        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // Move the completed carrier into the coefficient engine; its individual
        // words may have partial interfaces even when the input is fully indexed.
        let tensor: SymbolicTensor<PartialStructure> = tensor.into();
        // Keep each complete tensor word as an opaque key. Only its finite
        // scalar coefficient is put over a common denominator; in particular,
        // rational epsilon weights from the independent moments must cancel.
        for (word, coefficient) in
            tensor.coefficient_list(spenso::shadowing::TensorCollectFilter::<0>::TaggedTensors)?
        {
            if word.expression().is_zero() || coefficient.expression().is_zero() {
                continue;
            }
            ensure!(coefficient.is_scalar());
            let scalar = coefficient.expression().collect_compact_factors();
            if scalar.is_zero() {
                continue;
            }
            // Only finite scalar moment/epsilon weights may need rational
            // cancellation. Physical pole owners remain separate, and the
            // shared preflight bounds work before polynomial conversion.
            ensure!(
                !scalar.contains_symbol(GS.den),
                "{label}: physical poles are outside the finite scalar certificate"
            );
            // Exact finite Laurent regrouping needs no common denominator.
            // In particular, repeated epsilon/pi monomials can make a rational
            // size estimate pessimistic even when their complete sum is zero.
            let mut variables = BTreeSet::new();
            let mut exact = true;
            scalar.visitor(&mut |view| {
                if let AtomView::Num(number) = view {
                    exact &= matches!(
                        number.get_coeff_view(),
                        symbolica::coefficient::CoefficientView::Natural(..)
                            | symbolica::coefficient::CoefficientView::Large(..)
                    );
                } else if matches!(view, AtomView::Var(_)) {
                    variables.insert(view.to_owned());
                }
                true
            });
            let laurent = exact
                .then(|| scalar.coefficient_list_exact(&variables.into_iter().collect::<Vec<_>>()))
                .flatten();
            let exact_laurent_zero = laurent
                .as_ref()
                .is_some_and(|rows| rows.iter().all(|(_, coefficient)| coefficient.is_zero()));
            let rational_in_scope = ScalarRationalSize::estimate(scalar.as_view()).is_some();
            let in_scope = exact_laurent_zero || rational_in_scope;
            let is_zero = exact_laurent_zero || (rational_in_scope && scalar.together().is_zero());
            if !is_zero {
                crate::debug_tags!(#uv, #integrated, #inspect;
                    stage = "native_proper_prefix_scalar_certificate_rejected", label, in_scope,
                    physical_poles = scalar.contains_symbol(GS.den),
                    scalar_bytes = scalar.as_view().get_byte_size(),
                    file.scalar = %cooking.uncook(scalar.as_view()).to_canonical_string(),
                    file.word = %cooking.uncook(word.expression().as_view()).to_canonical_string(),
                    "Retained factorized scalar row and native word at the rejected proof boundary"
                );
            }
            ensure!(
                in_scope,
                "{label}: finite scalar scope refused; physical_poles={}, bytes={}, symbols={:?}",
                scalar.contains_symbol(GS.den),
                scalar.as_view().get_byte_size(),
                scalar
                    .get_all_symbols(true)
                    .into_iter()
                    .map(|symbol| symbol.to_string())
                    .collect::<BTreeSet<_>>()
            );
            ensure!(
                is_zero,
                "{label}: native finite-tensor coefficient is not zero; bytes={}, symbols={:?}",
                scalar.as_view().get_byte_size(),
                scalar
                    .get_all_symbols(true)
                    .into_iter()
                    .map(|symbol| symbol.to_string())
                    .collect::<BTreeSet<_>>()
            );
        }
        Ok(())
    }

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

    pub(in crate::uv) fn mass_projector(
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
            external_ordering: ordering_factor,
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
    use spenso::structure::{
        abstract_index::AbstractIndex,
        representation::LibraryRep,
        slot::{IsAbstractSlot, ParseableAind, Slot},
    };
    use std::collections::HashMap;

    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let sources = reference.certified_projections()?;
    ensure!(sources.len() == 13);
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    for (diagram, _) in &sources {
        let source = diagram.numerator().clone().try_into_indexed()?;
        let (projector, _) = reference.mass_projector(diagram, &Atom::var(crate::utils::GS.dim))?;
        let projector = projector.map_representation_dimensions(|_, d| d.replace(CS.nc).with(3))?;
        let projector = projector
            .with_checked_expression(projector.expression().replace(CS.nc).with(3))?
            .try_into_indexed()?;
        ensure!(source.rank() == 4 && projector.rank() == 4);
        let mut endpoints = BTreeMap::new();
        cooking
            .uncook(diagram.projector().expression().as_view())
            .visitor(&mut |view| {
                if let AtomView::Fun(wave) = view
                    && [symbols::u(), symbols::ubar()].contains(&wave.get_symbol())
                    && wave.get_nargs() == 2
                {
                    endpoints.insert(wave.get_symbol(), wave.get(1).to_owned());
                }
                true
            });
        ensure!(endpoints.len() == 2);
        let raw_slots = [&endpoints[&symbols::u()], &endpoints[&symbols::ubar()]];
        let typed = raw_slots
            .iter()
            .map(|slot| -> Result<_> {
                let cooked = cooking.try_cook(slot.as_view())?;
                Ok(Slot::<LibraryRep, AbstractIndex>::try_from(
                    cooked.as_view(),
                )?)
            })
            .collect::<Result<Vec<_>>>()?;
        let source_slots = source.structure().logical_slots();
        let projector_slots = projector.structure().logical_slots();
        let admit = |indices: &[Atom]| -> Result<()> {
            ensure!(indices.len() == 2 && indices[0] != indices[1]);
            for (slot, index) in typed.iter().zip(indices) {
                for frame in [&source_slots, &projector_slots] {
                    ensure!(
                        frame
                            .iter()
                            .filter(|candidate| candidate.rep() == slot.rep()
                                && candidate.aind().to_atom() == *index)
                            .count()
                            == 1,
                        "native source endpoint slot admission requires the complete cooked representation slot"
                    );
                }
            }
            Ok(())
        };
        let indices = typed
            .iter()
            .map(|slot| slot.aind().to_atom())
            .collect::<Vec<_>>();
        admit(&indices)?;
        let contextless = raw_slots
            .iter()
            .map(|slot| -> Result<Atom> {
                let AtomView::Fun(slot) = slot.as_view() else {
                    return Err(eyre!("native external spin slot is not explicit"));
                };
                Ok(cooking.try_cook(slot.get(1))?)
            })
            .collect::<Result<Vec<_>>>()?;
        let error = admit(&contextless)
            .expect_err("contextless hedge-index cooking must not establish a native source port");
        ensure!(
            error
                .to_string()
                .contains("native source endpoint slot admission")
        );
        for (label, tensor, opposite) in [
            ("closing projector", &projector, &source),
            ("native open source", &source, &projector),
        ] {
            let logical = tensor.structure().logical_slots();
            let mut swapped = logical.clone();
            let mut map = HashMap::new();
            let mut inverse = HashMap::new();
            for (position, slot) in logical.iter().enumerate() {
                let role = indices
                    .iter()
                    .position(|index| slot.aind().to_atom() == *index)
                    .ok_or_else(|| {
                        eyre!("{label}: external source port has no admitted endpoint")
                    })?;
                map.insert(position, typed[1 - role].aind());
                inverse.insert(position, slot.aind());
                swapped[position] = slot.reindex(typed[1 - role].aind());
            }
            let transported = tensor.reindex_interface_ports(&map)?;
            ensure!(transported.structure().logical_slots() == swapped);
            ensure!(
                transported.structure().canonical().external_structure()
                    == opposite.structure().canonical().external_structure()
            );
            ensure!(
                transported
                    .reindex_interface_ports(&inverse)?
                    .same_value(tensor)
            );
            let zero = tensor.with_checked_expression(Atom::Zero)?;
            let moved_zero = zero.reindex_interface_ports(&map)?;
            ensure!(
                moved_zero.rank() == 4
                    && moved_zero.structure().logical_slots()
                        == transported.structure().logical_slots()
            );
            ensure!(
                moved_zero
                    .reindex_interface_ports(&inverse)?
                    .same_value(&zero)
            );
        }
    }
    // Captured generic-field residual before the SU(3) scalar comparison.
    // Its logarithm cancels only after the native invariants are evaluated in
    // the same declared group as the independent numerical reference.
    let fundamental = ColorFundamental {}.to_symbolic([Atom::num(3)]);
    let adjoint = ColorAdjoint {}.to_symbolic([Atom::num(8)]);
    let cf = idenso::color_cas!(2, fundamental.clone());
    let ca = idenso::color_cas!(2, adjoint);
    let tf = idenso::color_idx!(2, fundamental);
    let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
    let ell = Atom::var(symbol!("native_msbar_field_reference::auxiliary_log"));
    let inverse_epsilon = epsilon.pow(-1);
    let double_pole = epsilon.pow(-2);
    let residual = (((Atom::num(12) * &ell + Atom::num(43)) * Atom::num(8) * &tf
        - Atom::num(153) * &ca)
        * &inverse_epsilon
        - Atom::num(3)
            * &cf
            * ((Atom::num(12) * &ell + Atom::num(7)) * &inverse_epsilon
                - Atom::num(18) * &double_pole)
        + (Atom::num(-8) * &tf + Atom::num(3) * &ca) * Atom::num(12) * &double_pole)
        * &cf
        - Atom::num(176) * &double_pole
        + Atom::num(420) * &inverse_epsilon;
    ensure!(!residual.together().is_zero());
    let cooked = cooking.try_cook(residual.as_view())?;
    let normalized = cooked.to_cof_dimension_invariants();
    ensure!(normalized.together().is_zero());
    ensure!(normalized.to_cof_dimension_invariants() == normalized);
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

// Actual outputs of the isolated-child certificate, reused by parent tests.
// This record exists only in the test module and does not change admission.
pub(in crate::uv) struct CertifiedVertexChild {
    pub(in crate::uv) diagram: FeynmanDiagram,
    pub(in crate::uv) graph: crate::graph::Graph,
    pub(in crate::uv) child: crate::uv::hedge_poset::OwnedForestNode,
    pub(in crate::uv) local: crate::uv::approx::local_4d::Local4dCts,
    pub(in crate::uv) integrated: super::IntegratedCts,
    pub(in crate::uv) settings: crate::uv::UVgenerationSettings,
    pub(in crate::uv) pole_reference: SymbolicTensor<IndexedStructure>,
}

impl QcdMassReference {
    // This checks valid standalone MUV children, never the guarded OS parent.
    pub(in crate::uv) fn certify_integrated_qg_children(
        &self,
        internal_pdgs: &[i64],
    ) -> Result<Vec<CertifiedVertexChild>> {
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
        let epsilon = Atom::var(GS.dim_epsilon);
        // Independent single-scale vacuum references at mu=M=1, SU(3).
        let pole_weight = if internal_pdgs == [6, 6, 21] {
            Atom::num((-1, 6))
        } else {
            Atom::num((9, 2))
        };
        let vertex_series = if internal_pdgs == [6, 6, 21] {
            let color = Atom::num((-1, 6)); // CF-CA/2
            &color
                * (epsilon.pow(-1) - Atom::num(2)
                    + &epsilon * (Atom::one() + pi.pow(2) / Atom::num(12)))
        } else {
            Atom::num(3)
                * (Atom::num((3, 2)) * epsilon.pow(-1) - Atom::one()
                    + &epsilon * pi.pow(2) / Atom::num(8))
        };
        let coupling = Atom::var(symbol!("UFO::G"));
        let loop_factor = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut count = 0usize;
        let mut records = Vec::new();
        let mut certified = Vec::new();
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
                    -&loop_factor * model.expand_couplings(tree.expression()) * &vertex_series,
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
                certified.push(CertifiedVertexChild {
                    diagram: diagram.clone(),
                    graph: graph.clone(),
                    child: OwnedForestNode {
                        spinney: current.spinney.clone(),
                        topo_order: current.topo_order,
                    },
                    local,
                    integrated,
                    settings: settings.clone(),
                    pole_reference: tree.with_checked_expression(
                        -&loop_factor
                            * model.expand_couplings(tree.expression())
                            * &pole_weight
                            * epsilon.pow(-1),
                    )?,
                });
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
        Ok(certified)
    }
}

#[test]
fn native_abelian_qg_children_match_integrated_open_tree_through_epsilon() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?
        .certify_integrated_qg_children(&[6, 6, 21])
        .map(|_| ())
}

#[test]
fn native_nonabelian_qg_children_match_integrated_open_tree_through_epsilon() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?
        .certify_integrated_qg_children(&[6, 21, 21])
        .map(|_| ())
}

#[test]
fn native_qg_vertex_proper_prefixes_match_one_loop_mass_density() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certify_qg_vertex_proper_prefixes()?;
    Ok(())
}

impl QcdMassReference {
    fn certify_qg_vertex_proper_prefixes(&self) -> Result<CertifiedMassPrefixes> {
        use crate::{
            graph::{LMBext, lmb::LMBwithEdges},
            utils::GS,
            uv::{
                ApproximationType, UltravioletGraph,
                approx::{
                    ForestNodeLike,
                    local_4d::{FourDTerm, Full4dCts},
                },
                hedge_poset::OwnedForestNode,
            },
        };
        use linnet::half_edge::subgraph::{Inclusion, SubSetLike, SubSetOps};

        let reference = self;
        ensure!(GS.emr_mom == symbols::momentum());
        let one_sources = reference.sources(1)?;
        let [one_source] = one_sources.as_slice() else {
            return Err(eyre!("one native one-loop source"));
        };
        let one = reference.project(one_source)?;
        one.certify_one_loop()?;
        let d = Atom::var(GS.dim);
        let eps = Atom::var(GS.dim_epsilon);
        let coupling = Atom::var(symbol!("UFO::G"));
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let a = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
        let reference_d = symbol!("native_mass_reference::d");
        let coordinates = one
            .family
            .scalar_products()
            .iter()
            .enumerate()
            .map(|(i, product)| {
                (
                    product.clone(),
                    Atom::var(symbol!(&format!("native_mass_reference::sp{i}"))),
                )
            })
            .collect::<Vec<_>>();
        let coordinate_rules = coordinates
            .iter()
            .map(|(product, variable)| {
                Replacement::new(product.to_pattern(), variable.to_pattern())
            })
            .collect::<Vec<_>>();
        let k_ref = &one.family.loop_momenta()[0];
        let p_ref = &one.family.external_momenta()[0];
        let kk = one
            .family
            .kinematics()
            .scalar_product(k_ref, k_ref)?
            .replace_multiple(&coordinate_rules);
        let kp = one
            .family
            .kinematics()
            .scalar_product(k_ref, p_ref)?
            .replace_multiple(&coordinate_rules);
        let reference_numerator = one
            .scalar
            .replace(reference_d)
            .with(d.clone())
            .replace(symbol!("native_mass_reference::CF"))
            .with(Atom::num((4, 3)));
        let reference_poles = one
            .poles
            .iter()
            .map(|(_, _, _, full)| full.replace_multiple(&coordinate_rules))
            .collect::<Vec<_>>();
        ensure!(reference_poles.len() == 2);
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut records = reference.certify_integrated_qg_children(&[6, 6, 21])?;
        records.extend(reference.certify_integrated_qg_children(&[6, 21, 21])?);
        ensure!(records.len() == 4);
        let mut counts = BTreeMap::<String, usize>::new();
        let mut prefixes = CertifiedMassPrefixes::new();
        for record in records {
            let graph = &record.graph;
            let child = &record.child;
            let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
            let parents = graph
                .classified_spinneys(
                    &uncut,
                    &graph.model,
                    &record.settings,
                    &graph.loop_momentum_basis,
                )?
                .into_iter()
                .filter(|node| {
                    node.n_components() == 1
                        && node.dod == 1
                        && node.renormalization_scheme == ApproximationType::OS
                        && graph.n_loops(node.filter()) == 2
                        && node.filter().includes(child.subgraph())
                })
                .collect::<Vec<_>>();
            ensure!(parents.len() == 1);
            let parent = OwnedForestNode {
                spinney: parents[0].clone(),
                topo_order: 2,
            };
            let reduced = parent.reduced_subgraph(child);
            let edges = graph.paired_edges(&reduced);
            ensure!(edges.len() == 2);
            let mut quark = None;
            let mut gluon = None;
            for edge in &edges {
                match graph
                    .model
                    .particle_by_id(graph[*edge].particle.particle().unwrap())?
                    .pdg_code
                    .abs()
                {
                    6 => {
                        ensure!(quark.replace(*edge).is_none());
                    }
                    21 => {
                        ensure!(gluon.replace(*edge).is_none());
                    }
                    pdg => return Err(eyre!("unexpected quotient particle {pdg}")),
                }
                ensure!(
                    child.lmb().edge_signatures[*edge]
                        .internal
                        .iter()
                        .all(|s| s.is_zero())
                );
            }
            let quark = quark.ok_or_else(|| eyre!("missing quotient quark"))?;
            let gluon = gluon.ok_or_else(|| eyre!("missing quotient gluon"))?;
            ensure!(child.lmb().loop_edges.len() == 1);
            let chart = graph
                .lmb_with_loop_edges(&[*child.lmb().loop_edges.first().unwrap(), quark][..])?;
            let routing =
                graph
                    .normal_emr_replacement(&graph.full_filter(), &chart, &[] as &[Atom], |_| true);
            let q = function!(GS.emr_mom, quark.0);
            ensure!(q.replace_multiple(&routing) == q);
            let (projector, p_native) = reference.mass_projector(&record.diagram, &d)?;
            let p_routed = p_native.replace_multiple(&routing);
            let basis = chart
                .loop_edges
                .iter()
                .chain(chart.ext_edges.iter())
                .map(|edge| function!(GS.emr_mom, edge.0))
                .collect::<BTreeSet<_>>();
            let declared = Kinematics::in_dimension(&d)?.with_momenta(basis)?;
            let p_terms = declared.linear_terms(&p_routed)?;
            let [(p, p_sign)] = p_terms.as_slice() else {
                return Err(eyre!("one native external basis coordinate"));
            };
            ensure!(*p_sign == Atom::one() || *p_sign == Atom::num(-1));
            let gluon_routed = function!(GS.emr_mom, gluon.0).replace_multiple(&routing);
            let mut orientations = Vec::new();
            for sign in [1, -1] {
                let expected = Atom::num(sign) * &q - &p_routed;
                if (&gluon_routed - &expected).expand().is_zero()
                    || (&gluon_routed + &expected).expand().is_zero()
                {
                    orientations.push(sign);
                }
            }
            let [quark_sign] = orientations.as_slice() else {
                return Err(eyre!("unique signed quotient loop map"));
            };
            let mut kin = declared
                .with_scalar_product(&q, &q, kk.clone())?
                .with_scalar_product(p, p, Atom::one())?
                .with_scalar_product(&q, p, Atom::num(*quark_sign) * p_sign * &kp)?;
            // Add native edge-product aliases to the same exact signed chart; an
            // unassigned consumed-child product remains symbolic and is rejected.
            let routed = (0..graph.underlying.n_edges())
                .map(|edge| {
                    let native = function!(GS.emr_mom, edge);
                    let value = native.replace_multiple(&routing);
                    (native, value)
                })
                .collect::<Vec<_>>();
            let chart_kin = kin.clone();
            for (left, l) in &routed {
                for (right, r) in &routed {
                    kin = kin.with_scalar_product(left, right, chart_kin.scalar_product(l, r)?)?;
                }
            }
            let denominator = graph.denominator(&reduced, &graph.model, |_| -1);
            let terms = FourDTerm::from_view(denominator.as_view())?;
            ensure!(
                terms.len() == 1
                    && terms[0].denominators.len() == 2
                    && terms[0].numerator.is_scalar()
                    && terms[0].numerator.expression().is_one()
            );
            ensure!(
                terms[0]
                    .denominators
                    .iter()
                    .map(|pole| pole.source_edge)
                    .collect::<BTreeSet<_>>()
                    == BTreeSet::from([quark, gluon])
            );
            let mut poles = Vec::new();
            for pole in &terms[0].denominators {
                let polynomial = pole.polynomial_in_lmb(&chart)?;
                let scalar = kin.apply(&polynomial).replace(symbol!("UFO::MT")).with(1);
                let expected = &reference_poles[usize::from(pole.source_edge != quark)];
                ensure!(
                    (&scalar - expected).expand().is_zero(),
                    "complete native quotient denominator differs"
                );
                poles.push(serde_json::json!({"edge":pole.source_edge.0,
                    "momentum":pole.momentum.to_canonical_string(),
                    "mass_squared":pole.mass_squared.to_canonical_string(),
                    "full_polynomial":pole.full_expr.to_canonical_string(),
                    "routed_polynomial":polynomial.to_canonical_string()}));
            }
            let prefix = Full4dCts::recursion_input(
                &record.local,
                &record.integrated,
                ApproximationType::MUV,
                false,
                child.lmb(),
            )?;
            let finite = prefix
                .sectors()
                .filter(|sector| !sector.frozen_lmbs().is_empty())
                .collect::<Vec<_>>();
            ensure!(finite.len() == 1 && finite[0].active_components.is_empty());
            ensure!(finite[0].frozen_lmbs() == std::slice::from_ref(child.lmb()));
            let frozen_completion = finite[0].atom.with_checked_expression(
                finite[0]
                    .atom
                    .expression()
                    .replace(GS.integrated_loop_scale)
                    .with(1),
            )?;
            let completed = record
                .integrated
                .physical_pole()?
                .try_sub(&record.integrated.physical_finite_counterterm()?)?
                .try_add(&frozen_completion)?;
            let cograph = graph
                .numerator(&reduced, child.subgraph())?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            let empty = graph.empty_subgraph();
            let child_source = graph
                .numerator(child.subgraph(), &empty)?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            let parent_source = graph
                .numerator(parent.subgraph(), &empty)?
                .to_d_dim(GS.dim)?
                .state
                .expr;
            ensure!(
                child_source
                    .multiply(cograph.clone())?
                    .try_sub(&parent_source)?
                    .is_zero(),
                "actual native child/cograph factor ownership"
            );
            let projector = projector
                .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
            let projector =
                projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
            let projector = projector.try_into_indexed()?;
            let closed = completed.multiply(cograph)?.multiply(projector)?;
            let closed_source = cooking.uncook(closed.expression().as_view());
            let mut wrong_dimension = false;
            closed_source.visitor(&mut |view| {
                if let AtomView::Fun(rep) = view
                    && rep.get_symbol() == symbol!("spenso::mink")
                {
                    wrong_dimension |= rep.get(0) != d.as_view();
                }
                true
            });
            ensure!(
                !wrong_dimension,
                "newly internal child Lorentz slots remain formal d before mass trace"
            );
            let cooked = cooking.try_cook(closed_source.as_view())?;
            ensure!(cooking.uncook(cooked.as_view()) == closed_source);
            let tensor = &closed;
            ensure!(tensor.structure().canonical().is_scalar());
            let contracted = tensor
                .simplify_algebra(&AlgebraSettings {
                    gamma: Some(GammaSimplifySettings::canonical()),
                    color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                    ..Default::default()
                })?
                .contract(ContractSettings::default())?;
            crate::debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_proper_prefix_trace", source = %record.diagram.name(),
                child = %child.subgraph().string_label(),
                contraction_complete = contracted.contraction_complete(), rank = contracted.rank(),
                file.closed_source = %closed_source.to_canonical_string(),
                file.contracted = %cooking.uncook(contracted.expression().as_view()).to_canonical_string(),
                "Actual formal-d child/cograph mass trace before scalar-product notation"
            );
            ensure!(
                contracted.contraction_complete() && contracted.structure().canonical().is_scalar()
            );
            let dotted = contracted.to_dots()?;
            ensure!(dotted.structure().canonical().is_scalar());
            let trace = cooking.uncook(dotted.expression().as_view());
            let scalar = graph
                .model
                .expand_couplings(&kin.apply(&trace))
                .replace(symbol!("UFO::MT"))
                .with(1)
                .replace(GS.m_uv_vacuum)
                .with(1)
                .replace(GS.m_uv_expansion)
                .with(1)
                .replace(GS.mu_r_sq)
                .with(1);
            let source = reference.project(&record.diagram)?;
            ensure!(source.kernel_weight.is_one() && source.numerator_prefactor.is_one());
            // Child loop measure is already present. Only the remaining outer
            // measure and Γ2=-iΣ conversion enter this normalized two-loop density.
            let actual = (Atom::i()
                * (Atom::i() / (Atom::num(16) * pi.pow(2)))
                * &source.kernel_weight
                * &source.numerator_prefactor
                * scalar
                / a.pow(2))
            .collect_compact_factors();
            let is_abelian = graph
                .paired_edges(child.subgraph())
                .iter()
                .filter(|edge| {
                    graph
                        .model
                        .particle_by_id(graph[**edge].particle.particle().unwrap())
                        .is_ok_and(|particle| particle.pdg_code.abs() == 6)
                })
                .count()
                == 2;
            let vertex_pole = Atom::num(if is_abelian { (1, 6) } else { (-9, 2) }) / &eps;
            let expected = &vertex_pole * &reference_numerator;
            crate::debug_tags!(#uv, #integrated, #inspect;
                stage = "native_vertex_proper_prefix_density", source = %record.diagram.name(),
                child = %child.subgraph().string_label(),
                file.actual = %actual.to_canonical_string(), file.expected = %expected.to_canonical_string(),
                file.poles = %serde_json::to_string(&poles)?,
                file.native_external = %p_native.to_canonical_string(),
                file.routed_external = %p_routed.to_canonical_string(), quark_sign = *quark_sign,
                "Actual completed vertex times untouched native cograph before integration"
            );
            let difference = dotted.with_checked_expression(&actual - &expected)?;
            ensure!(
                super::simplify(&difference)?.is_zero(),
                "actual native vertex proper prefix differs from the independent one-loop mass density"
            );
            // The independent one-loop integral is retained through epsilon^1:
            // multiplication by the actual vertex pole fixes the finite two-loop term.
            let f = Atom::num((4, 3))
                * (Atom::num(3) / &eps
                    + Atom::num(4)
                    + &eps * (Atom::num(8) + pi.pow(2) / Atom::num(4)));
            let integrated_reference = (&vertex_pole * f).series(GS.dim_epsilon, 0, 0)?.to_atom();
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_vertex_proper_prefix_certificate", source = %record.diagram.name(),
                child = %child.subgraph().string_label(),
                file.integrated_reference = %integrated_reference.to_canonical_string(),
                highest_one_loop_epsilon_power = 1, public_parent_admitted = false,
                "Exact SU(3) density certificate transfers the independently integrated one-loop reference"
            );
            let key = (
                record.diagram.name().to_owned(),
                graph
                    .paired_edges(child.subgraph())
                    .iter()
                    .map(|edge| edge.0)
                    .collect(),
            );
            ensure!(
                prefixes
                    .insert(
                        key,
                        CertifiedMassPrefix {
                            coefficient: integrated_reference,
                            diagram: record.diagram.clone(),
                            graph: record.graph.clone(),
                            given: OwnedForestNode {
                                spinney: child.spinney.clone(),
                                topo_order: child.topo_order,
                            },
                            completed: prefix.clone(),
                            settings: record.settings.clone(),
                        }
                    )
                    .is_none(),
                "a native vertex occurrence was counted twice"
            );
            *counts.entry(record.diagram.name().to_owned()).or_default() += 1;
        }
        ensure!(counts.values().copied().collect::<Vec<_>>() == [2, 2]);
        ensure!(prefixes.len() == 4);
        Ok(prefixes)
    }
}

#[path = "native_gluon_prefix_tests.rs"]
mod native_gluon_prefix_tests;
#[path = "native_self_energy_prefix_tests.rs"]
mod native_self_energy_prefix_tests;

#[test]
fn native_two_loop_proper_prefix_sum_matches_the_independent_mass_reference() -> Result<()> {
    crate::initialisation::test_initialise()?;
    QcdMassReference::new()?.certified_mass_prefixes()?;
    Ok(())
}

impl QcdMassReference {
    fn certified_mass_prefixes(&self) -> Result<CertifiedMassPrefixes> {
        use crate::uv::{UltravioletGraph, approx::ForestNodeLike};
        let reference = self;
        let mut prefixes = CertifiedMassPrefixes::new();
        for (kind, expected_count, certified) in [
            ("OS quark", 1, reference.certify_os_quark_proper_prefix()?),
            (
                "QCD vertex",
                4,
                reference.certify_qg_vertex_proper_prefixes()?,
            ),
            (
                "fermion bubble",
                6,
                reference.certify_fermion_gluon_prefixes()?,
            ),
            (
                "pure gauge",
                3,
                reference.certify_pure_gauge_gluon_prefixes()?,
            ),
        ] {
            ensure!(
                certified.len() == expected_count,
                "{kind}: proper-child incidence"
            );
            for (key, coefficient) in certified {
                ensure!(
                    prefixes.insert(key, coefficient).is_none(),
                    "a proper occurrence was counted twice"
                );
            }
        }
        let expected_keys = BTreeSet::from([
            ("D00".to_owned(), vec![5, 6]),
            ("D01".to_owned(), vec![5, 6]),
            ("D02".to_owned(), vec![2, 4, 5]),
            ("D02".to_owned(), vec![3, 4, 6]),
            ("D03".to_owned(), vec![5]),
            ("D04".to_owned(), vec![5, 6]),
            ("D05".to_owned(), vec![5, 6]),
            ("D06".to_owned(), vec![5, 6]),
            ("D07".to_owned(), vec![5, 6]),
            ("D08".to_owned(), vec![5, 6]),
            ("D09".to_owned(), vec![2, 3]),
            ("D10".to_owned(), vec![2, 3]),
            ("D11".to_owned(), vec![2, 3, 5]),
            ("D11".to_owned(), vec![2, 4, 6]),
        ]);
        ensure!(prefixes.keys().cloned().collect::<BTreeSet<_>>() == expected_keys);
        ensure!(
            prefixes[&("D03".to_owned(), vec![5])].coefficient.is_zero(),
            "the certified zero tadpole remains in the forest inventory"
        );
        let actual = super::series(
            &prefixes
                .values()
                .map(|record| &record.coefficient)
                .sum::<Atom>(),
            0,
        )?
        .to_atom();
        let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let cf = Atom::num((4, 3));
        let beta0 = Atom::num(7); // SU(3), five light and one heavy active flavor.
        let f = &cf
            * (Atom::num(3) / &epsilon
                + Atom::num(4)
                + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(4)));
        let h = &cf
            * (Atom::num(6) / &epsilon
                + Atom::num(2)
                + &epsilon * (Atom::num(8) + pi.pow(2) / Atom::num(2)));
        let z1 = -&f;
        let w1 = -&cf / &epsilon;
        // Independent Lagrangian/Slavnov-Taylor combination: quark z1 H-w1 F,
        // two vertex occurrences 2v1 F, and gluon -a1 F. No actual overall
        // counterterm or bare two-loop coefficient defines this expected P2.
        let expected = super::series(&(&z1 * &h - &beta0 * &f / &epsilon + &w1 * &f), 0)?.to_atom();
        let explicit = -Atom::num((196, 3)) / epsilon.pow(2)
            - Atom::num((880, 9)) / &epsilon
            - Atom::num((2080, 9))
            - Atom::num((73, 9)) * pi.pow(2);
        ensure!((&expected - &explicit).together().is_zero());
        ensure!(
            (&actual - &explicit).together().is_zero(),
            "the fourteen actual proper insertions differ from independent P2"
        );

        // Keeping only epsilon^0 in both one-loop reference integrals preserves
        // the poles but loses a known nonzero finite two-loop coefficient.
        let f0 = &cf * (Atom::num(3) / &epsilon + Atom::num(4));
        let h0 = &cf * (Atom::num(6) / &epsilon + Atom::num(2));
        let prematurely_truncated =
            super::series(&(-&f0 * &h0 - &beta0 * &f0 / &epsilon + &w1 * &f0), 0)?.to_atom();
        let lost_finite = -Atom::num((1952, 9)) - Atom::num((73, 9)) * pi.pow(2);
        ensure!(!lost_finite.is_zero());
        ensure!(
            (&actual - &prematurely_truncated - &lost_finite)
                .together()
                .is_zero()
        );
        let records = prefixes.iter().map(|((source, edges), coefficient)| serde_json::json!({
        "source": source, "child_edges": edges, "coefficient": coefficient.coefficient.to_canonical_string()
    })).collect::<Vec<_>>();
        crate::debug_tags!(#uv, #integrated, #summary;
            stage = "native_two_loop_complete_proper_prefix_sum", occurrences = prefixes.len(),
            file.records = %serde_json::to_string(&records)?,
            file.actual = %actual.to_canonical_string(), file.reference = %explicit.to_canonical_string(),
            file.lost_finite = %lost_finite.to_canonical_string(),
            parent_overall_operator_computed = false, public_complete_forest_evaluated = false,
            emitted_four_d_condition_proved = false,
            "Actual native child/cograph certificates transfer all fourteen proper insertions to independent P2"
        );
        for ((source, edges), record) in &prefixes {
            ensure!(record.diagram.name() == source && record.graph.name == *source);
            ensure!(
                record
                    .graph
                    .paired_edges(record.given.subgraph())
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>()
                    == *edges
            );
            ensure!(record.settings.generate_integrated && record.given.topo_order == 1);
            let mut active_domains = 0;
            for sector in record.completed.sectors() {
                for (owner, scope, lmb) in &sector.active_components {
                    ensure!(owner == record.given.subgraph() && scope == owner);
                    ensure!(lmb.loop_edges.len() == record.graph.n_loops(owner));
                    active_domains += 1;
                }
                for frozen in sector.frozen_lmbs() {
                    ensure!(frozen == record.given.lmb());
                }
            }
            ensure!(
                active_domains == 1,
                "retain the actual child integration domain"
            );
        }
        Ok(prefixes)
    }
}

// First overall continuation gate: the same retained D09 graph owns both inputs.
// This does not claim the full thirteen-source / fourteen-prefix K2 condition.
#[test]
fn native_d09_bare_and_completed_child_preserve_overall_mass_condition() -> Result<()> {
    use crate::{
        utils::GS,
        uv::{
            ApproximationType, Spinney,
            approx::{
                ForestNodeLike, Rooted, UVCtx,
                local_4d::{Full4dCts, Local4dCts},
            },
            hedge_poset::OwnedForestNode,
        },
    };
    use linnet::half_edge::subgraph::SubSetLike;

    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    // Do not reimport the graph after computing its child. The record includes
    // the actual local and finite sectors with their original dummy scopes.
    let mut records = reference.certify_os_quark_proper_prefix()?;
    ensure!(records.len() == 1);
    let record = records
        .remove(&("D09".to_owned(), vec![2, 3]))
        .ok_or_else(|| eyre!("the certified native D09 child is missing"))?;
    let graph = &record.graph;
    ensure!(record.diagram.name() == "D09" && graph.name == "D09");
    let settings = &record.settings;
    ensure!(settings.generate_integrated && !settings.add_marker);
    ensure!(settings.on_shell_velocity()? == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
    ensure!(record.given.renormalization_scheme() == ApproximationType::OS);
    let root = OwnedForestNode {
        spinney: Spinney::empty(graph),
        topo_order: 0,
    };
    let bare = Full4dCts::from_factorized_local(&Local4dCts::root());
    let mut completed = 0;
    for (label, given, prefix) in [
        ("bare", &root, &bare),
        ("completed OS child", &record.given, &record.completed),
    ] {
        crate::debug_tags!(#uv, #integrated, #inspect;
            stage = "native_d09_overall_input", label, source = %record.diagram.name(),
            given = %given.subgraph().string_label(),
            file.prefix = %cooking.uncook(prefix.atom().expression().as_view()).to_canonical_string(),
            file.velocity = ?settings.on_shell_velocity()?,
            "Same native graph and retained component sectors enter the overall operator"
        );
        let actual = reference.certify_actual_overall_mass_input(&record, given, prefix)?;
        reference.certify_d09_rational_finalizer_mass_projection(
            &record.diagram,
            &UVCtx::new(graph, &graph.model, settings),
            given,
            prefix,
            &actual,
        )?;
        let (_, _, integrated) = actual;
        if given.subgraph().is_empty() {
            ensure!(
                !integrated.physical_finite_counterterm()?.is_zero(),
                "the bare rainbow has a nonzero finite field correction before mass projection"
            );
        }
        completed += 1;
        crate::debug_tags!(#uv, #integrated, #summary;
            stage = "native_d09_overall_input_complete", label, completed,
            local_parent_certificate = true, actual_integrated_finite_projection_zero = true,
            rational_finalizer_mass_condition = true, emitted_cff_condition_proved = true,
            public_complete_forest_evaluated = false,
            "D09 actual parent operator and its finite scheme insertion obey the matched mass projection"
        );
    }
    ensure!(completed == 2);
    // The consumed loop marker is checked at its owner; it is never treated
    // as a physical momentum or an extra finite counterterm here.
    ensure!(
        record
            .completed
            .atom()
            .expression()
            .contains_symbol(GS.integrated_loop_scale)
    );
    Ok(())
}

impl QcdMassReference {
    /// Extract the complete generic-external-momentum coefficient of one actual
    /// signed overall pole, retaining native spin/color ports and source weight.
    fn certify_native_overall_kinetic_pole(
        &self,
        record: &CertifiedMassPrefix,
        parent: &crate::uv::hedge_poset::OwnedForestNode,
        actual: &super::IntegratedCts,
        projection: &MassProjection,
        common: &mut Option<NativeFieldFrame>,
    ) -> Result<(
        SymbolicTensor<IndexedStructure>,
        SymbolicTensor<IndexedStructure>,
    )> {
        use crate::{
            numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt},
            utils::GS,
            uv::{UltravioletGraph, approx::ForestNodeLike},
        };
        use idenso::tensor::{SymbolicNetExt, SymbolicNetParse};
        use linnet::half_edge::{involution::Flow, subgraph::SubSetLike};
        use spenso::{
            network::parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
            structure::{
                abstract_index::AbstractIndex,
                representation::LibraryRep,
                slot::{IsAbstractSlot, ParseableAind, Slot},
            },
        };
        use std::collections::HashMap;
        use symbolica::atom::FunctionBuilder;

        let diagram = &record.diagram;
        let graph = &record.graph;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let pole = actual.physical_pole()?;
        ensure!(
            pole.rank() == 4,
            "a native quark pole retains two spin and two color ports"
        );
        for forbidden in [
            GS.den,
            GS.loop_mom,
            GS.localizing_integrand,
            GS.integrated_loop_scale,
            vakint::symbols::S.k,
            vakint::symbols::S.p,
        ] {
            ensure!(
                !pole.expression().contains_symbol(forbidden),
                "an internal denominator/momentum/localization owner survived actual overall integration"
            );
        }
        let (projector, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector = projector
            .with_checked_expression(projector.expression().replace(CS.nc).with(3))?
            .try_into_indexed()?;
        let AtomView::Fun(reference_p_call) = reference_p.as_view() else {
            return Err(eyre!(
                "native incoming momentum must retain its source edge"
            ));
        };
        ensure!(reference_p_call.get_symbol() == GS.emr_mom && reference_p_call.get_nargs() == 1);
        let incoming_edge =
            usize::try_from(reference_p_call.get(0)).map_err(|error| eyre!(error))?;
        let routed_p = diagram.loop_momentum_basis().route_expression(&reference_p);
        let crown = graph.dummy_less_full_crown(parent.subgraph());
        let boundary = crown.included_iter().collect::<Vec<_>>();
        ensure!(boundary.len() == 2);
        let mut signs = BTreeMap::new();
        let mut physical_mass = None;
        for hedge in boundary {
            let edge = graph[&hedge];
            let routed = diagram
                .loop_momentum_basis()
                .route_expression(&function!(GS.emr_mom, edge.0));
            let sign = if (&routed - &routed_p).expand().is_zero() {
                1
            } else if (&routed + &routed_p).expand().is_zero() {
                -1
            } else {
                return Err(eyre!(
                    "kinetic pole boundary is not the native signed external momentum"
                ));
            };
            if graph.flow(hedge) == Flow::Sink {
                ensure!(sign == 1 && physical_mass.is_none());
                physical_mass = Some(graph[edge].particle.mass_atom(&graph.model));
            }
            ensure!(signs.insert(edge.0, sign).is_none());
        }
        let physical_mass = physical_mass
            .ok_or_else(|| eyre!("the actual parent has no unique physical incoming mass owner"))?;
        for tensor in [&pole, &projector] {
            let mut unexpected = Vec::new();
            tensor.expression().visitor(&mut |view| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && (momentum.get_nargs() == 0
                        || usize::try_from(momentum.get(0))
                            .ok()
                            .is_none_or(|edge| !signs.contains_key(&edge)))
                {
                    unexpected.push(view.to_owned());
                }
                true
            });
            ensure!(
                unexpected.is_empty(),
                "internal/unowned external momenta survived integration: {unexpected:?}"
            );
        }
        // The closing mass projector owns both physical endpoint labels. Its
        // p-linear numerator is pslash/(4 Nc); transpose those two native
        // endpoints to obtain the independently ordered OPEN kinetic matrix.
        let native_projector = cooking.uncook(diagram.projector().expression().as_view());
        let mut endpoints = BTreeMap::new();
        native_projector.visitor(&mut |view| {
            if let AtomView::Fun(wave) = view
                && [symbols::u(), symbols::ubar()].contains(&wave.get_symbol())
                && wave.get_nargs() == 2
                && let AtomView::Fun(spin) = wave.get(1)
                && spin.get_symbol() == symbol!("spenso::bis")
                && spin.get_nargs() == 2
            {
                endpoints.insert(wave.get_symbol(), Atom::from(spin.to_owned()));
            }
            true
        });
        ensure!(endpoints.len() == 2);
        let endpoint_atoms = [
            endpoints[&symbols::u()].clone(),
            endpoints[&symbols::ubar()].clone(),
        ];
        let cooked_endpoints = endpoint_atoms
            .iter()
            .map(|slot| -> Result<Atom> {
                // Index encoding is owned by a complete representation slot.
                let cooked = cooking.try_cook(slot.as_view())?;
                let typed = Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view())?;
                Ok(typed.aind().to_atom())
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(cooked_endpoints[0] != cooked_endpoints[1]);
        let zero_projector = projector.with_checked_expression(
            projector.expression().replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    **output = Atom::Zero;
                }
            }),
        )?;
        let closing_kinetic = projector.try_sub(&zero_projector)?;
        let closing_kinetic = closing_kinetic
            .with_checked_expression(Atom::num(12) * closing_kinetic.expression())?;
        let slots = closing_kinetic.structure().logical_slots();
        let endpoint_indices = cooked_endpoints
            .iter()
            .map(|index| {
                slots
                    .iter()
                    .find(|slot| slot.aind().to_atom() == *index)
                    .map(|slot| slot.aind())
                    .ok_or_else(|| eyre!("native endpoint is absent from the kinetic template"))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut transpose = HashMap::new();
        let mut inverse = HashMap::new();
        for (position, slot) in slots.iter().enumerate() {
            let role = cooked_endpoints
                .iter()
                .position(|index| slot.aind().to_atom() == *index)
                .ok_or_else(|| eyre!("kinetic template has an unowned port"))?;
            transpose.insert(position, endpoint_indices[1 - role]);
            inverse.insert(position, slot.aind());
        }
        ensure!(transpose.len() == 4);
        let kinetic_template = closing_kinetic.reindex_interface_ports(&transpose)?;
        ensure!(
            kinetic_template
                .reindex_interface_ports(&inverse)?
                .same_value(&closing_kinetic)
        );
        let template_zero = closing_kinetic.with_checked_expression(Atom::Zero)?;
        let transposed_zero = template_zero.reindex_interface_ports(&transpose)?;
        ensure!(
            transposed_zero.rank() == 4
                && transposed_zero
                    .reindex_interface_ports(&inverse)?
                    .same_value(&template_zero)
        );
        ensure!(
            kinetic_template
                .structure()
                .canonical()
                .external_structure()
                == pole.structure().canonical().external_structure()
        );
        // This unit matrix uses only the original projector's fixed source
        // interface, never a produced pole coefficient. Its endpoint transport
        // has its own inverse, even for a ranked-zero publication.
        let mass_identity =
            zero_projector.with_checked_expression(Atom::num(12) * zero_projector.expression())?;
        let mut identity_transpose = HashMap::new();
        let mut identity_inverse = HashMap::new();
        for (position, slot) in mass_identity.structure().logical_slots().iter().enumerate() {
            let role = cooked_endpoints
                .iter()
                .position(|index| slot.aind().to_atom() == *index)
                .ok_or_else(|| eyre!("unit mass matrix has an unowned source port"))?;
            identity_transpose.insert(position, endpoint_indices[1 - role]);
            identity_inverse.insert(position, slot.aind());
        }
        ensure!(identity_transpose.len() == 4);
        let transposed_identity = mass_identity.reindex_interface_ports(&identity_transpose)?;
        ensure!(
            transposed_identity
                .reindex_interface_ports(&identity_inverse)?
                .same_value(&mass_identity)
        );
        let zero_identity = mass_identity.with_checked_expression(Atom::Zero)?;
        let moved_zero = zero_identity.reindex_interface_ports(&identity_transpose)?;
        ensure!(
            moved_zero.rank() == 4
                && moved_zero
                    .reindex_interface_ports(&identity_inverse)?
                    .same_value(&zero_identity)
        );
        ensure!(
            transposed_identity
                .structure()
                .canonical()
                .external_structure()
                == pole.structure().canonical().external_structure()
        );
        // Phi(I)=1 in this original unit-mass projector normalization, for
        // generic p: the odd gamma trace is zero and Tr_spin,color(I)=4 Nc.
        let unit_trace =
            super::simplify(&projector.clone().multiply(transposed_identity.clone())?)?;
        ensure!(unit_trace.is_scalar());
        Self::certify_finite_tagged_rows_zero(
            "independent unit matrix projector calibration",
            unit_trace.with_checked_expression(unit_trace.expression() - Atom::one())?,
        )?;
        let first_input = common.is_none();
        if first_input {
            *common = Some(NativeFieldFrame {
                kinetic: kinetic_template.clone(),
                identity: transposed_identity.clone(),
                mass: physical_mass.clone(),
                endpoint_slots: endpoint_atoms.clone(),
                incoming_edge,
            });
        }
        let NativeFieldFrame {
            kinetic: common_template,
            identity: common_identity,
            mass: common_mass,
            endpoint_slots: common_endpoint_atoms,
            incoming_edge: common_edge,
        } = common.as_ref().unwrap();
        ensure!(
            physical_mass == common_mass,
            "native inputs must retain the same symbolic physical top mass"
        );
        let common_endpoints = common_endpoint_atoms
            .iter()
            .map(|slot| -> Result<Atom> {
                // Index encoding is owned by a complete representation slot.
                let cooked = cooking.try_cook(slot.as_view())?;
                let typed = Slot::<LibraryRep, AbstractIndex>::try_from(cooked.as_view())?;
                Ok(typed.aind().to_atom())
            })
            .collect::<Result<Vec<_>>>()?;
        let common_slots = common_template.structure().logical_slots();
        let route = |expression: &Atom| -> Atom {
            expression.replace_map(|view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    let edge =
                        usize::try_from(momentum.get(0)).expect("all source Q owners were checked");
                    let mut target = FunctionBuilder::new(GS.emr_mom).add_arg(*common_edge);
                    for index in momentum.iter().skip(1) {
                        target = target.add_arg(index);
                    }
                    **output = Atom::num(signs[&edge]) * target.finish();
                }
            })
        };
        let align = |tensor: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
            let logical = tensor.structure().logical_slots();
            ensure!(logical.len() == common_slots.len());
            let mut map = HashMap::new();
            let mut inverse = HashMap::new();
            for (position, slot) in logical.iter().enumerate() {
                let role = cooked_endpoints.iter().position(|index| slot.aind().to_atom() == *index)
                    .ok_or_else(|| eyre!("overall pole has an unnamed external port"))?;
                let targets = common_slots.iter().filter(|target| target.rep() == slot.rep() && target.aind().to_atom() == common_endpoints[role]).collect::<Vec<_>>();
                ensure!(targets.len() == 1, "native spin/color role has no unique common-frame port");
                map.insert(position, targets[0].aind());
                inverse.insert(position, slot.aind());
            }
            let aligned = tensor.reindex_interface_ports(&map)?;
            ensure!(aligned.reindex_interface_ports(&inverse)?.same_value(tensor));
            let aligned_zero = tensor.with_checked_expression(Atom::Zero)?.reindex_interface_ports(&map)?;
            ensure!(aligned_zero.rank() == tensor.rank());
            ensure!(aligned_zero.reindex_interface_ports(&inverse)?.same_value(&tensor.with_checked_expression(Atom::Zero)?));
            let logical = aligned.structure().logical_slots();
            let axes = common_slots.iter().map(|slot| logical.iter().position(|candidate| candidate == slot)
                .ok_or_else(|| eyre!("common kinetic frame lost a logical slot"))).collect::<Result<Vec<_>>>()?;
            ensure!(axes.iter().copied().collect::<BTreeSet<_>>() == (0..axes.len()).collect());
            let partial: SymbolicTensor<PartialStructure> = aligned.into();
            let ordered = partial.permuted(&axes)?;
            let mut unpermute = vec![0; axes.len()];
            for (target, source) in axes.iter().copied().enumerate() { unpermute[source] = target; }
            ensure!(ordered.permuted(&unpermute)?.same_value(&partial));
            let ordered = ordered.try_into_indexed()?;
            ensure!(ordered.structure().logical_slots() == common_slots);
            Ok(ordered)
        };
        let routed_template =
            kinetic_template.with_checked_expression(route(kinetic_template.expression()))?;
        let aligned_template = align(&routed_template)?;
        Self::certify_finite_tagged_rows_zero(
            "native kinetic templates have the same generic-p open frame",
            super::simplify(&aligned_template.try_sub(common_template)?)?,
        )?;
        let aligned_identity = align(&transposed_identity)?;
        Self::certify_finite_tagged_rows_zero(
            "native unit matrices have the same generic-p open frame",
            super::simplify(&aligned_identity.try_sub(common_identity)?)?,
        )?;

        // Expose only compact vector arguments, retaining chain/trace factors.
        // A scalar tau inside an opaque gamma(tau*p) is not a coefficient.
        let expanded = cooking
            .uncook(pole.expression().as_view())
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
        let expanded = pole.with_checked_expression(cooking.try_cook(expanded.as_view())?)?;
        let tau = Atom::var(symbol!("native_msbar_field_reference::external_scale"));
        ensure!(!expanded.expression().contains(tau.as_view()));
        let scaled = expanded.with_checked_expression(expanded.expression().replace_map(
            |view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                {
                    **output = &tau * view.to_owned();
                }
            },
        ))?;
        let affine = |expression: &Atom| -> Result<(Atom, Atom)> {
            let rows = expression
                .coefficient_list_exact(std::slice::from_ref(&tau))
                .ok_or_else(|| eyre!("bounded complete external polynomial extraction refused"))?;
            let mut constant = Atom::Zero;
            let mut linear = Atom::Zero;
            for (monomial, coefficient) in rows {
                ensure!(
                    !coefficient.contains(tau.as_view()),
                    "external scale remains hidden in a coefficient"
                );
                if monomial == Atom::one() {
                    constant += coefficient;
                } else if monomial == tau {
                    linear += coefficient;
                } else {
                    return Err(eyre!(
                        "actual overall field pole is not affine in its generic external momentum: {monomial}"
                    ));
                }
            }
            Ok((constant, linear))
        };
        let (constant, linear) = affine(scaled.expression())?;
        if first_input {
            let quadratic = scaled.expression() + tau.pow(2) * kinetic_template.expression();
            let error = affine(&quadratic).expect_err(
                "the complete affine gate must reject a nonzero quadratic kinetic perturbation",
            );
            ensure!(error.to_string().contains("is not affine"));
        }
        let reconstructed = scaled.with_checked_expression(&constant + &tau * &linear)?;
        Self::certify_finite_tagged_rows_zero(
            "complete external-affine overall pole reconstruction",
            super::simplify(&scaled.try_sub(&reconstructed)?)?,
        )?;
        let linear = pole.with_checked_expression(linear)?;
        let linear = linear.with_checked_expression(route(linear.expression()))?;
        let linear = align(&linear)?;
        let weight = &projection.kernel_weight * &projection.numerator_prefactor;
        ensure!(!weight.contains_symbol(GS.emr_mom) && !weight.contains(tau.as_view()));
        let weighted = linear.with_checked_expression(
            self.model
                .expand_couplings(&(&weight * linear.expression()))
                .to_cof_dimension_invariants(),
        )?;
        let constant_tensor = pole.with_checked_expression(constant.clone())?;
        ensure!(
            !constant_tensor.expression().contains_symbol(GS.emr_mom),
            "the complete tau-zero coefficient cannot retain live external momentum"
        );
        let constant_tensor = align(&constant_tensor)?;
        let weighted_constant = constant_tensor.with_checked_expression(
            self.model
                .expand_couplings(&(weight * constant_tensor.expression()))
                .to_cof_dimension_invariants(),
        )?;
        crate::debug_tags!(#uv, #integrated, #summary;
            stage = "native_overall_kinetic_pole_input", source = %diagram.name(),
            generic_external_affine = true, exact_open_frame = true, su3_color_invariants_evaluated = true,
            file.external_signs = ?signs, file.scalar_source_weight = %(&projection.kernel_weight * &projection.numerator_prefactor).to_canonical_string(),
            file.constant_pole = %cooking.uncook(constant.as_view()).to_canonical_string(),
            file.weighted_linear_pole = %cooking.uncook(weighted.expression().as_view()).to_canonical_string(),
            "Actual integrated pole retains generic live momentum and ordered native spin/color output"
        );
        Ok((weighted, weighted_constant))
    }

    /// Apply both actual parent mass certificates to one original retained input.
    /// Retain the same parent and both actual outputs for physical-boundary proofs.
    fn certify_actual_overall_mass_input(
        &self,
        record: &CertifiedMassPrefix,
        given: &crate::uv::hedge_poset::OwnedForestNode,
        prefix: &crate::uv::approx::local_4d::Full4dCts,
    ) -> Result<(
        crate::uv::hedge_poset::OwnedForestNode,
        crate::uv::approx::local_4d::Local4dCts,
        super::IntegratedCts,
    )> {
        use crate::uv::{
            ApproximationType, UltravioletGraph,
            approx::{
                ForestNodeLike, Rooted, UVCtx,
                local_4d::{Full4dCts, Local4dCts},
            },
            hedge_poset::OwnedForestNode,
        };
        use linnet::half_edge::subgraph::{Inclusion, SubSetLike, SubSetOps};
        let graph = &record.graph;
        let settings = &record.settings;
        ensure!(settings.generate_integrated && !settings.add_marker);
        ensure!(settings.on_shell_velocity()? == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
        if given.subgraph().is_empty() {
            ensure!(given.topo_order == 0);
            ensure!(prefix == &Full4dCts::from_factorized_local(&Local4dCts::root()));
        } else {
            ensure!(given.subgraph() == record.given.subgraph());
            ensure!(
                given.lmb() == record.given.lmb()
                    && given.topo_order == record.given.topo_order
                    && given.dod() == record.given.dod()
                    && given.renormalization_scheme() == record.given.renormalization_scheme()
            );
            ensure!(prefix == &record.completed);
        }
        let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
        let mut parents = graph
            .classified_spinneys(&uncut, &graph.model, settings, &graph.loop_momentum_basis)?
            .into_iter()
            .filter(|component| {
                component.n_components() == 1
                    && component.dod == 1
                    && component.renormalization_scheme == ApproximationType::OS
                    && graph.n_loops(component.filter()) == 2
                    && component.filter().includes(record.given.subgraph())
            })
            .collect::<Vec<_>>();
        ensure!(parents.len() == 1);
        let parent = OwnedForestNode {
            spinney: parents.remove(0),
            topo_order: given.topo_order + 1,
        };
        let ctx = UVCtx::new(graph, &graph.model, settings);
        let backend = vakint::VakintSettings {
            number_of_terms_in_epsilon_expansion: 3,
            evaluation_order: vakint::EvaluationOrder(vec![vakint::EvaluationMethod::AlphaLoop(
                vakint::AlphaLoopOptions {
                    susbstitute_masters: false,
                },
            )]),
            ..settings.vakint.true_settings()
        };
        let local = self.certify_parent_mass_projection_input(
            &record.diagram,
            &ctx,
            &parent,
            given,
            prefix,
        )?;
        let integrated = super::Integrated::new(crate::utils::vakint()?, &backend)
            .run(&local, &ctx, &parent, given, &parent, given)?;
        let label = if given.subgraph().is_empty() {
            "bare".to_owned()
        } else {
            format!("completed {} child", given.renormalization_scheme())
        };
        self.certify_parent_integrated_finite_mass_projection(
            &label,
            &record.diagram,
            &ctx,
            &parent,
            &integrated,
        )?;
        Ok((parent, local, integrated))
    }

    fn certify_parent_integrated_finite_mass_projection(
        &self,
        label: &str,
        diagram: &FeynmanDiagram,
        ctx: &crate::uv::approx::UVCtx<'_>,
        parent: &crate::uv::hedge_poset::OwnedForestNode,
        actual: &super::IntegratedCts,
    ) -> Result<()> {
        use crate::{
            utils::GS,
            uv::{UltravioletGraph, approx::ForestNodeLike},
        };
        use linnet::half_edge::{
            involution::{EdgeIndex, Flow},
            subgraph::SubSetLike,
        };
        let graph = ctx.graph;
        let crown = graph.dummy_less_full_crown(parent.subgraph());
        let boundary = crown.included_iter().collect::<Vec<_>>();
        ensure!(boundary.len() == 2);
        let incoming = boundary
            .iter()
            .copied()
            .find(|h| graph.flow(*h) == Flow::Sink)
            .ok_or_else(|| eyre!("native parent incoming port missing"))?;
        let owner = graph[&incoming];
        let mass = graph[owner].particle.mass_atom(&graph.model);
        let (projector, reference_p) = self.mass_projector(diagram, &Atom::var(GS.dim))?;
        // The generic reference retains Nc; the actual runtime graph is SU(3).
        // Map both the native representations and scalar coefficients before
        // the owned publication into the runtime's explicit Indexed frame.
        let projector = projector
            .map_representation_dimensions(|_, dimension| dimension.replace(CS.nc).with(3))?;
        let projector =
            projector.with_checked_expression(projector.expression().replace(CS.nc).with(3))?;
        let projector = projector.try_into_indexed()?;
        let basis = diagram.loop_momentum_basis();
        let routed_p = basis.route_expression(&reference_p);
        let owner_p = basis.route_expression(&function!(GS.emr_mom, owner.0));
        ensure!(
            (&owner_p - &routed_p).expand().is_zero(),
            "incoming owner sign differs from the native mass projector"
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
                    "finite-parent boundary is not the same signed external momentum"
                ));
            };
            signs.insert(edge, sign);
        }
        let velocity = ctx.settings.on_shell_velocity()?;
        ensure!(velocity == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
        let finite = actual.physical_finite_counterterm()?;
        ensure!(
            !finite.expression().contains_symbol(GS.den)
                && !finite.expression().contains_symbol(GS.loop_mom)
                && !finite.expression().contains_symbol(vakint::symbols::S.k)
                && !finite.expression().contains_symbol(vakint::symbols::S.p),
            "actual parent integration must consume hard momenta and restore every backend alias"
        );
        // Integration must consume all loop momenta. An external Q occurrence
        // is admitted only by its actual signed boundary signature above.
        for value in [finite.expression(), projector.expression()] {
            let mut unexpected = Vec::new();
            value.visitor(&mut |view| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == GS.emr_mom
                    && (momentum.get_nargs() == 0
                        || usize::try_from(momentum.get(0))
                            .ok()
                            .is_none_or(|edge| !signs.contains_key(&EdgeIndex(edge))))
                {
                    unexpected.push(view.to_owned());
                }
                true
            });
            ensure!(
                unexpected.is_empty(),
                "{label}: active/internal momenta survived actual parent integration: {unexpected:?}"
            );
        }
        let at_point =
            |value: &SymbolicTensor<IndexedStructure>| -> Result<SymbolicTensor<IndexedStructure>> {
                let mapped = GS
                    .erase_uv_momentum_provenance(value.expression())
                    .replace_map(|view, _, output| {
                        if let AtomView::Fun(momentum) = view
                            && momentum.get_symbol() == GS.emr_mom
                            && let Ok(edge) = usize::try_from(momentum.get(0))
                        {
                            let sign = signs[&EdgeIndex(edge)];
                            let point =
                                GS.on_shell_momentum(owner, &(Atom::num(sign) * &mass), &velocity);
                            let indices = momentum
                                .iter()
                                .skip(1)
                                .map(|index| index.to_owned())
                                .collect::<Vec<_>>();
                            **output = GS.indexed_momentum(&point, &indices);
                        }
                    });
                // Unit physical mass is imposed only after all parent operators
                // and actual integration, never on active auxiliary vacuum scales.
                let expression = GS
                    .resolve_os_momenta(&mapped)
                    .replace(mass.to_pattern())
                    .with(1);
                Ok(value.with_checked_expression(expression)?)
            };
        let projected = at_point(&projector)?.multiply(at_point(&finite)?)?;
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        ensure!(
            projected.structure().canonical().is_scalar(),
            "finite parent mass projection leaves a native port open"
        );
        // After integration only fixed external vectors remain. Native algebra
        // identifies gamma0 and the temporal basis slash without introducing
        // the reversible frame used by the unintegrated parent-source proof.
        let reduced = projected
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::canonical()),
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            })?
            .contract(ContractSettings::default())?;
        ensure!(
            reduced.contraction_complete(),
            "finite parent mass projection must contract before Laurent truncation"
        );
        let scalar = cooking.uncook(reduced.to_dots()?.expression().as_view());
        let regulated = super::Integrated::dimensionally_regularized(&scalar);
        let coefficient = super::series(&regulated, 0)?.to_atom();
        let admitted = reduced.with_checked_expression(cooking.try_cook(coefficient.as_view())?)?;
        ensure!(admitted.is_scalar());
        let zero = admitted
            .expression()
            .collect_compact_factors()
            .together()
            .is_zero();
        crate::debug_tags!(#uv, #integrated, #inspect;
            stage = "native_parent_integrated_finite_mass_projection", label,
            file.raw_finite = %cooking.uncook(finite.expression().as_view()).to_canonical_string(),
            file.projected = %cooking.uncook(projected.expression().as_view()).to_canonical_string(),
            file.reduced = %coefficient.to_canonical_string(),
            file.boundary_signs = ?signs, auxiliary_scales_substituted = false,
            exact_zero = zero, epsilon_output_depth = 0,
            "Actual Integrated output is projected after full formal-d contraction and only then Laurent truncated"
        );
        ensure!(
            zero,
            "{label}: actual integrated finite parent insertion has a nonzero matched mass projection"
        );
        Ok(())
    }
}

// Source/forest coverage prerequisite only. The actual twenty-six parent
// operators and their finite outputs must pass before the K2 transfer is used.
#[test]
fn native_overall_mass_inputs_retain_certified_sources_and_proper_occurrences() -> Result<()> {
    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let sources = reference.certified_projections()?;
    let prefixes = reference.certified_mass_prefixes()?;
    reference.certify_overall_mass_input_inventory(&sources, &prefixes)?;
    let _ = QcdMassReference::independent_bare_and_lagrangian_mass_references()?;
    Ok(())
}

impl QcdMassReference {
    fn certify_overall_mass_input_inventory(
        &self,
        sources: &[(FeynmanDiagram, MassProjection)],
        prefixes: &CertifiedMassPrefixes,
    ) -> Result<()> {
        use crate::{
            settings::global::GenerationSettings,
            utils::GS,
            uv::{
                UltravioletGraph,
                approx::{CutStructure, ForestNodeLike},
                hedge_poset::Wood,
            },
        };
        use linnet::half_edge::subgraph::SubSetOps;
        use spenso::structure::partial::PartialStructureExt;
        ensure!(sources.len() == 13 && prefixes.len() == 14);
        let one_loop = sources
            .iter()
            .filter(|(d, _)| d.loop_momentum_basis().loop_edges.len() == 1)
            .count();
        ensure!(one_loop == 1);
        let two_loop = sources
            .iter()
            .filter(|(d, _)| d.loop_momentum_basis().loop_edges.len() == 2)
            .collect::<Vec<_>>();
        ensure!(two_loop.len() == 12);
        let mut all_inputs = BTreeSet::new();
        let mut source_keys = BTreeSet::new();
        let mut node_count = 0;
        for (diagram, projected) in two_loop {
            let matching = prefixes
                .iter()
                .filter(|((source, _), _)| source == diagram.name())
                .collect::<Vec<_>>();
            ensure!(
                !matching.is_empty(),
                "every native source retains its proper-child record"
            );
            let expected_given = matching
                .iter()
                .map(|((_, edges), _)| edges.clone())
                .chain([Vec::new()])
                .collect::<BTreeSet<_>>();
            ensure!(expected_given.len() == matching.len() + 1);
            ensure!(source_keys.insert(diagram.name().to_owned()));
            for ((source, edges), record) in &matching {
                // These checks relate each retained occurrence to the same
                // native13 source. They never rebuild its runtime graph or
                // replace the already-completed prefix's index scopes.
                let actual = &record.diagram;
                ensure!(actual.name() == diagram.name() && source.as_str() == diagram.name());
                ensure!(actual.loop_momentum_basis() == diagram.loop_momentum_basis());
                ensure!(actual.numerator() == diagram.numerator());
                ensure!(actual.projector() == diagram.projector());
                ensure!(actual.overall_factor() == diagram.overall_factor());
                ensure!(actual.numerator_prefactor() == diagram.numerator_prefactor());
                ensure!(actual.denominator_expression()? == diagram.denominator_expression()?);
                ensure!(actual.vertices().count() == diagram.vertices().count());
                for ((id, vertex), (expected_id, expected)) in
                    actual.vertices().zip(diagram.vertices())
                {
                    ensure!(
                        id == expected_id
                            && vertex.name == expected.name
                            && vertex.interaction == expected.interaction
                            && vertex.numerator == expected.numerator
                    );
                }
                ensure!(actual.edges().count() == diagram.edges().count());
                for ((id, endpoints, edge), (expected_id, expected_endpoints, expected)) in
                    actual.edges().zip(diagram.edges())
                {
                    ensure!(id == expected_id && endpoints == expected_endpoints);
                    ensure!(
                        edge.particle == expected.particle
                            && edge.external == expected.external
                            && edge.is_dummy == expected.is_dummy
                            && edge.directed == expected.directed
                            && edge.source_slot() == expected.source_slot()
                            && edge.target_slot() == expected.target_slot()
                            && edge.numerator == expected.numerator
                    );
                }
                let graph = &record.graph;
                ensure!(graph.name == diagram.name());
                ensure!(graph.overall_factor == *diagram.overall_factor());
                let source_prefactor = diagram.numerator_prefactor().clone().try_into_indexed()?;
                let source_projector = diagram.projector().clone().try_into_indexed()?;
                ensure!(graph.global_prefactor.num.same_value(&source_prefactor));
                ensure!(
                    graph
                        .global_prefactor
                        .projector
                        .same_value(&source_projector)
                );
                let runtime_basis: crate::graph::LoopMomentumBasis = actual
                    .loop_momentum_basis()
                    .to_routing(actual.underlying())
                    .into();
                ensure!(graph.loop_momentum_basis == runtime_basis);
                let full = graph.full_filter().subtract(&graph.initial_state_cut.left);
                let runtime_numerator = graph
                    .numerator(&full, &graph.empty_subgraph())?
                    .to_d_dim(GS.dim)?
                    .state
                    .expr;
                let reference_numerator = actual
                    .numerator()
                    .map_representation_dimensions(|rep, dimension| {
                        if rep == (Minkowski {}).into() {
                            Atom::var(GS.dim)
                        } else {
                            dimension.to_owned()
                        }
                    })?
                    .try_into_indexed()?;
                // These independently composed native sources can expose the
                // same named ports in different logical encounter order. Align
                // only the private reference view through the native axis owner.
                ensure!(
                    runtime_numerator
                        .structure()
                        .canonical()
                        .external_structure()
                        == reference_numerator
                            .structure()
                            .canonical()
                            .external_structure()
                );
                let runtime_slots = runtime_numerator.structure().logical_slots();
                let reference_slots = reference_numerator.structure().logical_slots();
                let axes = runtime_slots
                    .iter()
                    .map(|slot| {
                        reference_slots
                            .iter()
                            .position(|candidate| candidate == slot)
                            .ok_or_else(|| {
                                eyre!("native reference is missing a runtime source port")
                            })
                    })
                    .collect::<Result<Vec<_>>>()?;
                ensure!(
                    axes.iter().copied().collect::<BTreeSet<_>>()
                        == (0..runtime_slots.len()).collect::<BTreeSet<_>>()
                );
                let reference_work: SymbolicTensor<PartialStructure> = reference_numerator.into();
                let aligned = reference_work.permuted(&axes)?;
                ensure!(aligned.expression() == reference_work.expression());
                let mut inverse = vec![0; axes.len()];
                for (target, source) in axes.iter().copied().enumerate() {
                    inverse[source] = target;
                }
                ensure!(aligned.permuted(&inverse)?.same_value(&reference_work));
                // The same transport must retain order even when the payload
                // cannot witness any ports.
                let zero_reference = reference_work.with_checked_expression(Atom::Zero)?;
                let zero_aligned = zero_reference.permuted(&axes)?;
                ensure!(zero_aligned.is_zero() && zero_aligned.rank() == reference_work.rank());
                ensure!(
                    zero_aligned.structure().logical_slots() == aligned.structure().logical_slots()
                );
                ensure!(zero_aligned.permuted(&inverse)?.same_value(&zero_reference));
                let zero_aligned = zero_aligned.try_into_indexed()?;
                ensure!(zero_aligned.structure().logical_slots() == runtime_slots);
                let reference_numerator = aligned.try_into_indexed()?;
                ensure!(reference_numerator.structure().logical_slots() == runtime_slots);
                ensure!(
                    reference_numerator
                        .structure()
                        .canonical()
                        .external_structure()
                        == runtime_numerator
                            .structure()
                            .canonical()
                            .external_structure()
                );
                if !runtime_numerator.same_value(&reference_numerator) {
                    let cooking = CookSettings::indices()
                        .with_mode(CookMode::ReversibleEncoding)
                        .with_representation_payloads(true, true);
                    crate::debug_tags!(#uv, #integrated, #inspect;
                        stage = "native_inventory_source_numerator_mismatch",
                        source = %diagram.name(), file.given_edges = ?edges,
                        expression_equal = runtime_numerator.expression() == reference_numerator.expression(),
                        runtime_rank = runtime_numerator.rank(), reference_rank = reference_numerator.rank(),
                        runtime_contracted = runtime_numerator.contraction_complete(),
                        reference_contracted = reference_numerator.contraction_complete(),
                        runtime_reduction = ?runtime_numerator.reduction_status(),
                        reference_reduction = ?reference_numerator.reduction_status(),
                        file.runtime = %cooking.uncook(runtime_numerator.expression().as_view()).to_canonical_string(),
                        file.reference = %cooking.uncook(reference_numerator.expression().as_view()).to_canonical_string(),
                        file.runtime_canonical_slots = ?runtime_numerator.structure().canonical().external_structure(),
                        file.reference_canonical_slots = ?reference_numerator.structure().canonical().external_structure(),
                        file.runtime_logical_slots = ?runtime_numerator.structure().logical_slots(),
                        file.reference_logical_slots = ?reference_numerator.structure().logical_slots(),
                        "Retained runtime and independently mapped native sources first differ before pole comparison"
                    );
                }
                ensure!(
                    runtime_numerator.same_value(&reference_numerator),
                    "retained runtime source must preserve its native factorized numerator"
                );
                ensure!(
                    graph.denominator(&full, &graph.model, |_| 1)
                        == actual.denominator_expression()?,
                    "original full poles, repeated owners and prescriptions must match before projection"
                );
                ensure!(
                    graph
                        .paired_edges(record.given.subgraph())
                        .iter()
                        .map(|edge| edge.0)
                        .collect::<Vec<_>>()
                        == *edges
                );
                for (name, expression) in &projected.coupling_definitions {
                    ensure!(
                        graph.model.coupling(name)?.expression == *expression,
                        "runtime coupling definition differs from the certified bare reference"
                    );
                }
                // All exact inputs of the existing native13 projection now
                // agree. Reuse that certified projection instead of tracing
                // the same source again for each child occurrence.
            }

            // Construct only the native classified inventory, with the same
            // original graph/settings as its first retained occurrence. The
            // guarded public counterterm dispatcher is not invoked or relaxed.
            let record = matching[0].1;
            let graph = &record.graph;
            let settings = GenerationSettings {
                uv: record.settings.clone(),
                ..Default::default()
            };
            ensure!(!settings.medium.vacuum_subtraction);
            let forests =
                Wood::new(CutStructure::empty(graph), graph, &graph.model, &settings)?.unfold();
            let full_edges = graph
                .paired_edges(&graph.full_filter())
                .into_iter()
                .map(|e| e.0)
                .collect::<Vec<_>>();
            let mut found = BTreeSet::new();
            for (node, _, operation) in forests.graph.iter_nodes() {
                node_count += 1;
                let Some(cover) = operation.covers() else {
                    ensure!(operation.key.op_count() == 0);
                    continue;
                };
                if graph
                    .paired_edges(&cover)
                    .into_iter()
                    .map(|e| e.0)
                    .collect::<Vec<_>>()
                    != full_edges
                {
                    continue;
                }
                ensure!(!forests.graph.is_disjoint_union(node));
                let (prior, _) = forests
                    .graph
                    .unique_parent(node)
                    .ok_or_else(|| eyre!("overall operation lacks one native given owner"))?;
                let given = forests.graph[prior]
                    .covers()
                    .map(|filter| {
                        graph
                            .paired_edges(&filter)
                            .into_iter()
                            .map(|e| e.0)
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                ensure!(operation.key.op_count() == if given.is_empty() { 1 } else { 2 });
                ensure!(found.insert(given.clone()));
                let unique_input = all_inputs.insert((diagram.name().to_owned(), given));
                ensure!(unique_input);
            }
            ensure!(
                found == expected_given,
                "actual forest operations differ from retained bare/proper input inventory"
            );
        }
        ensure!(source_keys.len() == 12 && all_inputs.len() == 26 && node_count == 52);
        ensure!(
            all_inputs
                .iter()
                .filter(|(_, edges)| edges.is_empty())
                .count()
                == 12
        );
        let proper_keys = all_inputs
            .iter()
            .filter(|(_, edges)| !edges.is_empty())
            .cloned()
            .collect::<BTreeSet<_>>();
        ensure!(proper_keys == prefixes.keys().cloned().collect());
        ensure!(
            all_inputs.contains(&("D03".to_owned(), Vec::new()))
                && all_inputs.contains(&("D03".to_owned(), vec![5]))
        );

        // Same serializer and pretty-JSON boundary as the existing numerical
        // callback. An artifact auditor can compare its exact UTF-8 hash with
        // the preserved nine-integral input; no task-specific hash is hardcoded.
        let payload = serde_json::to_vec_pretty(&serde_json::json!({
            "sources": sources.iter().map(|(diagram, projected)| projected.numeric_input(diagram, &self.model))
                .collect::<Result<Vec<_>>>()?
        }))?;
        crate::debug_tags!(#uv, #integrated, #inspect;
            stage = "native_overall_mass_source_and_inventory_bridge",
            native_sources = sources.len(), proper_occurrences = prefixes.len(), overall_inputs = all_inputs.len(),
            classified_nodes = node_count, file.overall_inputs = ?all_inputs,
            payload_bytes = payload.len(), file.numeric_payload = %String::from_utf8(payload)?,
            actual_overall_inputs_evaluated = 0, overall_K2_proved = false,
            "Exact retained-source/reference bridge and native forest incidence precede actual parent execution"
        );
        Ok(())
    }
}

impl QcdMassReference {
    /// Independent full-theory SU(3) references, nl=5, nh=1, mu=m=1,
    /// through epsilon^0. Neither expression is obtained from produced forests.
    fn independent_bare_and_lagrangian_mass_references() -> Result<(Atom, Atom)> {
        let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let zeta_three = Atom::var(vakint::vakint_symbol!("z3"));
        let log_two = function!(symbolica::atom::Symbol::LOG, 2);

        // Direct SU3 specialization of the independently maintained four-color
        // bare reference. The native13/source-payload bridge must be checked
        // before this value can be transferred to an actual O sector.
        let bare = Atom::num(38) / epsilon.pow(2)
            + Atom::num((289, 3)) / &epsilon
            + Atom::num((605, 2))
            + Atom::num((65, 9)) * pi.pow(2)
            - Atom::num((8, 3)) * &zeta_three
            + Atom::num((16, 9)) * pi.pow(2) * &log_two;

        // Parameter ZmOS from the primary reference, independently converted
        // to the same full-theory coupling and external field convention.
        let z2 =
            Atom::num(22) / epsilon.pow(2) - Atom::num((17, 3)) / &epsilon - Atom::num((1541, 18))
                + Atom::num((4, 9)) * pi.pow(2)
                + Atom::num((8, 3)) * &zeta_three
                - Atom::num((16, 9)) * pi.pow(2) * &log_two;
        let field_mass = Atom::num((16, 3)) / epsilon.pow(2)
            + Atom::num((64, 9)) / &epsilon
            + Atom::num((128, 9))
            + Atom::num((4, 9)) * pi.pow(2);
        let expected_k2 = Atom::num((82, 3)) / epsilon.pow(2) + Atom::num((13, 9)) / &epsilon
            - Atom::num((1285, 18))
            + Atom::num((8, 9)) * pi.pow(2)
            + Atom::num((8, 3)) * zeta_three
            - Atom::num((16, 9)) * pi.pow(2) * log_two;
        ensure!((&z2 + &field_mass - &expected_k2).together().is_zero());
        ensure!(
            !(&expected_k2 - z2).together().is_zero(),
            "omitting the field-mass product loses an independently nonzero coefficient"
        );

        // Reference-only consistency check, deliberately independent of every
        // actual prefix value and actual overall counterterm.
        let independent_p2 = -Atom::num((196, 3)) / epsilon.pow(2)
            - Atom::num((880, 9)) / &epsilon
            - Atom::num((2080, 9))
            - Atom::num((73, 9)) * pi.pow(2);
        ensure!(
            (&bare + &independent_p2 + &expected_k2)
                .together()
                .is_zero()
        );
        Ok((bare, expected_k2))
    }
}

#[test]
fn native_two_loop_overall_os_mass_condition_matches_independent_reference() -> Result<()> {
    use crate::uv::{
        Spinney,
        approx::{
            ForestNodeLike, Rooted,
            local_4d::{Full4dCts, Local4dCts},
        },
        hedge_poset::OwnedForestNode,
    };
    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let sources = reference.certified_projections()?;
    let prefixes = reference.certified_mass_prefixes()?;
    reference.certify_overall_mass_input_inventory(&sources, &prefixes)?;
    let mut completed = BTreeSet::new();
    let mut kinetic_frame = None;
    let mut kinetic_sum: Option<SymbolicTensor<IndexedStructure>> = None;
    let mut constant_sum: Option<SymbolicTensor<IndexedStructure>> = None;
    for (diagram, projection) in sources
        .iter()
        .filter(|(d, _)| d.loop_momentum_basis().loop_edges.len() == 2)
    {
        let matching = prefixes
            .iter()
            .filter(|((source, _), _)| source == diagram.name())
            .collect::<Vec<_>>();
        ensure!(!matching.is_empty());
        // Every source has at least one retained proper occurrence. Its first
        // graph also owns the bare input. Each later proper input uses its OWN
        // exact retained graph and never borrows another occurrence's dummy frame.
        for input in 0..=matching.len() {
            let record = matching[input.saturating_sub(1)].1;
            let graph = &record.graph;
            let root = OwnedForestNode {
                spinney: Spinney::empty(graph),
                topo_order: 0,
            };
            let bare = Full4dCts::from_factorized_local(&Local4dCts::root());
            let (given, prefix) = if input == 0 {
                (&root, &bare)
            } else {
                (&record.given, &record.completed)
            };
            let key = (
                diagram.name().to_owned(),
                graph
                    .paired_edges(given.subgraph())
                    .into_iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>(),
            );
            ensure!(!completed.contains(&key));
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_overall_mass_input", file.key = ?key,
                ordinal = completed.len() + 1, source = %record.diagram.name(),
                given_scheme = %given.renormalization_scheme(),
                "Starting one retained native bare/proper overall mass input"
            );
            let actual = reference.certify_actual_overall_mass_input(record, given, prefix)?;
            let (kinetic, constant) = reference.certify_native_overall_kinetic_pole(
                record,
                &actual.0,
                &actual.2,
                projection,
                &mut kinetic_frame,
            )?;
            kinetic_sum = Some(match kinetic_sum {
                Some(sum) => sum.try_add(&kinetic)?,
                None => kinetic,
            });
            constant_sum = Some(match constant_sum {
                Some(sum) => sum.try_add(&constant)?,
                None => constant,
            });
            let unique = completed.insert(key.clone());
            ensure!(unique);
            crate::debug_tags!(#uv, #integrated, #summary;
                stage = "native_overall_mass_input_complete", file.key = ?key,
                completed = completed.len(), native_parent_certificate = true,
                actual_integrated_finite_mass_projection_zero = true,
                public_dispatch_admitted = false,
                "One original native bare/proper input completes both actual parent mass gates"
            );
        }
    }
    ensure!(completed.len() == 26);
    ensure!(
        completed
            .iter()
            .filter(|(_, edges)| edges.is_empty())
            .count()
            == 12
    );
    ensure!(
        completed
            .iter()
            .filter(|(_, edges)| !edges.is_empty())
            .cloned()
            .collect::<BTreeSet<_>>()
            == prefixes.keys().cloned().collect()
    );
    ensure!(
        completed.contains(&("D03".to_owned(), Vec::new()))
            && completed.contains(&("D03".to_owned(), vec![5]))
    );

    // Independent MS field reference, with a=alpha_s/(4*pi), psi_B=sqrt(Zpsi)psi_R,
    // and renormalized Feynman gauge xi_R=1, xi_B=Z3*xi_R. Falcioni--Gardi--Milloy,
    // arXiv:1909.00697 A.18: CF[(CA+CF/2)/eps^2+(-17CA/4+3CF/4+TF*nf)/eps].
    // SU(3), nf=6 gives 44/(9 eps^2)-35/(3 eps). This target is not a definition
    // by subtraction of the produced bare/proper coefficients above.
    let NativeFieldFrame {
        kinetic: kinetic_template,
        identity: mass_identity,
        mass: physical_mass,
        incoming_edge: common_edge,
        ..
    } = kinetic_frame
        .ok_or_else(|| eyre!("the complete native inventory has no open kinetic frame"))?;
    let kinetic_sum = kinetic_sum.ok_or_else(|| eyre!("the actual pole sum is empty"))?;
    let epsilon = Atom::var(crate::utils::GS.dim_epsilon);
    let pi = Atom::var(symbolica::atom::Symbol::PI);
    let coupling = Atom::var(symbol!("UFO::G"));
    let a = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
    let w2 = Atom::num((44, 9)) * epsilon.pow(-2) - Atom::num((35, 3)) / &epsilon;
    let normalization = Atom::i() * a.pow(2);
    let expected_kinetic = kinetic_template
        .with_checked_expression(&normalization * &w2 * kinetic_template.expression())?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "complete generic-p native two-loop MS field pole",
        super::simplify(&kinetic_sum.try_sub(&expected_kinetic)?)?,
    )?;
    // Independent hybrid Lagrangian: Gamma_CT2=i*a^2[w2*pslash-(w2+K2)*m].
    // The separately local O branch owns -i*a^2*K2*m, hence the analytic
    // V=-(U1-U1O) pole is i*a^2*w2*(pslash-m). No produced Phi-zero identity
    // defines this expected constant matrix. Eq.5.11 of arXiv:2203.11038v1
    // and the native one-loop bare-source anchor fix this same separation.
    let constant_sum =
        constant_sum.ok_or_else(|| eyre!("the actual constant pole sum is empty"))?;
    let expected_constant = mass_identity.with_checked_expression(
        -&normalization * &w2 * &physical_mass * mass_identity.expression(),
    )?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "complete generic-p native two-loop analytic constant pole",
        super::simplify(&constant_sum.try_sub(&expected_constant)?)?,
    )?;
    let full_pole = kinetic_sum.try_add(&constant_sum)?;
    let full_expected = expected_kinetic.try_add(&expected_constant)?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "complete generic-p native two-loop fieldlike analytic pole",
        super::simplify(&full_pole.try_sub(&full_expected)?)?,
    )?;
    let wrong_constant = expected_constant.with_checked_expression(
        expected_constant.expression()
            + &normalization * &physical_mass * mass_identity.expression() / &epsilon,
    )?;
    let error = QcdMassReference::certify_finite_tagged_rows_zero(
        "wrong independent analytic constant pole",
        super::simplify(&constant_sum.try_sub(&wrong_constant)?)?,
    ).expect_err("the complete open-matrix field certificate must reject a nonzero constant-mass perturbation");
    ensure!(
        error
            .to_string()
            .contains("native finite-tensor coefficient is not zero")
    );
    // Fixed bare xi_B=1 instead shifts the double pole by CF*z3=4/3.
    // A simple-pole shift is a second independent normalization control.
    for wrong in [
        &w2 + Atom::num((4, 3)) * epsilon.pow(-2),
        &w2 + epsilon.pow(-1),
    ] {
        let wrong = kinetic_template
            .with_checked_expression(&normalization * wrong * kinetic_template.expression())?;
        let error = QcdMassReference::certify_finite_tagged_rows_zero(
            "wrong independent field-pole convention",
            super::simplify(&kinetic_sum.try_sub(&wrong)?)?,
        )
        .expect_err("independently nonzero field-pole shifts must be rejected");
        ensure!(
            error
                .to_string()
                .contains("native finite-tensor coefficient is not zero")
        );
    }
    // A matched-point trace cannot certify the open matrix: this generic
    // (p.v)vslash-pslash perturbation vanishes at the rest point, yet is nonzero
    // before that substitution. It must be rejected by the same native owner.
    let gs = &crate::utils::GS;
    let temporal =
        kinetic_template.with_checked_expression(kinetic_template.expression().replace_map(
            |view, _, output| {
                if let AtomView::Fun(momentum) = view
                    && momentum.get_symbol() == gs.emr_mom
                    && usize::try_from(momentum.get(0)).ok() == Some(common_edge)
                    && momentum.get_nargs() == 2
                {
                    **output = gs.emr_mom(
                        linnet::half_edge::involution::EdgeIndex(common_edge),
                        gs.cind(0),
                    ) * gs.energy_delta(momentum.get(1).to_owned());
                }
            },
        ))?;
    let matched_only = temporal.try_sub(&kinetic_template)?;
    let error = QcdMassReference::certify_finite_tagged_rows_zero(
        "generic-p matched-only field perturbation",
        super::simplify(&matched_only)?,
    )
    .expect_err("the generic-p field certificate must see the matched-only perturbation");
    ensure!(
        error
            .to_string()
            .contains("native finite-tensor coefficient is not zero")
    );
    let at_rest = matched_only.with_checked_expression(matched_only.expression().replace_map(
        |view, _, output| {
            if let AtomView::Fun(momentum) = view
                && momentum.get_symbol() == gs.emr_mom
                && usize::try_from(momentum.get(0)).ok() == Some(common_edge)
                && momentum.get_nargs() == 2
            {
                **output = if momentum.get(1) == gs.cind(0).as_view() {
                    Atom::one()
                } else {
                    gs.energy_delta(momentum.get(1).to_owned())
                };
            }
        },
    ))?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "the same nonzero generic field perturbation vanishes at matched rest",
        super::simplify(&at_rest)?,
    )?;
    crate::debug_tags!(#uv, #integrated, #summary;
        stage = "native_two_loop_msbar_field_pole", overall_inputs = completed.len(),
        file.actual_generic_p = %kinetic_sum.expression().to_canonical_string(),
        file.actual_constant_pole = %constant_sum.expression().to_canonical_string(),
        file.symbolic_physical_mass = %physical_mass.to_canonical_string(),
        full_analytic_fieldlike_matrix = true, independent_unit_matrix_calibrated = true,
        constant_mass_perturbation_rejected = true,
        file.independent_w2 = %w2.to_canonical_string(),
        exact_generic_open_matrix = true, mass_trace_used_as_field_oracle = false,
        auxiliary_scales_substituted = false, renormalized_feynman_gauge = true,
        fixed_bare_gauge_control_rejected = true, quadratic_control_rejected = true,
        matched_point_only_control_rejected = true,
        "Actual weighted native poles match the independent two-loop MS field reference"
    );

    // The actual operator/source identities above transfer -O to the same
    // independent bare integrals and actual proper densities. Only this ACTUAL
    // side is assembled from those inputs. Expected K2 is the independent
    // source/Lagrangian reference z2+w1*z1, not a definition by subtraction.
    let (independent_bare, expected_k2) =
        QcdMassReference::independent_bare_and_lagrangian_mass_references()?;
    let actual_proper =
        super::series(&prefixes.values().map(|r| &r.coefficient).sum::<Atom>(), 0)?.to_atom();
    let actual_transferred_k2 = -&independent_bare - &actual_proper;
    ensure!((&actual_transferred_k2 - &expected_k2).together().is_zero());
    crate::debug_tags!(#uv, #integrated, #summary;
        stage = "native_two_loop_overall_mass_condition_by_transfer",
        native_sources = sources.len(), proper_occurrences = prefixes.len(), overall_inputs = completed.len(),
        file.bare_reference = %independent_bare.to_canonical_string(),
        file.actual_proper = %actual_proper.to_canonical_string(),
        file.actual_overall = %actual_transferred_k2.to_canonical_string(),
        file.independent_K2 = %expected_k2.to_canonical_string(),
        actual_parent_operators_computed = true, actual_finite_outputs_computed = true,
        bare_numerical_integration_repeated = false, public_complete_forest_evaluated = false,
        emitted_four_d_condition_proved = false,
        "All actual parent identities transfer certified native B/P inputs to independent source K2"
    );
    Ok(())
}

/// Main-push only: all three original W sources feed one paired, observed
/// finite-emission certificate. The reference never defines either carrier.
#[test]
fn native_d09_emitted_vacuum_policy_field_matches_actual_source_rows() -> Result<()> {
    use crate::{
        utils::GS,
        uv::{
            Spinney,
            approx::{
                Rooted, UVCtx,
                local_4d::{Full4dCts, Local4dCts},
            },
            hedge_poset::OwnedForestNode,
        },
    };
    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let mut records = reference.certify_os_quark_proper_prefix()?;
    ensure!(records.len() == 1);
    let record = records
        .remove(&("D09".to_owned(), vec![2, 3]))
        .ok_or_else(|| eyre!("the actual native D09 OS child is missing"))?;
    let root = OwnedForestNode {
        spinney: Spinney::empty(&record.graph),
        topo_order: 0,
    };
    let bare_prefix = Full4dCts::from_factorized_local(&Local4dCts::root());
    // One bare parent and one completed-child parent own all three source rows.
    let bare = reference.certify_actual_overall_mass_input(&record, &root, &bare_prefix)?;
    let bare_rows = reference.certify_d09_vacuum_scalar_rows(
        &record,
        &root,
        &bare_prefix,
        &bare,
        D09VacuumSource::Bare,
    )?;
    // A different supplied vector in the same formal space is not fixed by
    // the endpoint's selector. Reject it before any representation maps to 4.
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let unknown = symbol!(
        "native_d09_field_control::unknown_vector",
        tags = [
            spenso::network::tags::SPENSO_TAG.rank1.clone(),
            spenso::network::tags::SPENSO_TAG.tensor.clone()
        ]
    );
    let vector = function!(unknown, Minkowski {}.to_symbolic([GS.dim]));
    let basis = cooking.uncook(bare_rows.frame.open_basis[1].expression().as_view());
    let uncontrolled = basis.replace_map(|view, _, output| {
        if let AtomView::Fun(gamma) = view
            && gamma.get_symbol() == idenso::dirac::AGS.gamma0
        {
            **output = function!(
                idenso::dirac::AGS.gamma,
                gamma.get(0),
                gamma.get(1),
                &vector
            );
        }
    });
    let controlled = basis.replace_map(|view, _, output| {
        if let AtomView::Fun(gamma) = view
            && gamma.get_symbol() == idenso::dirac::AGS.gamma0
        {
            **output = function!(
                idenso::dirac::AGS.gamma,
                gamma.get(0),
                gamma.get(1),
                function!(GS.delta_vec, GS.cind(0), Minkowski {}.to_symbolic([GS.dim]))
            );
        }
    });
    ensure!(uncontrolled.contains_symbol(unknown));
    ensure!(controlled.contains_symbol(GS.delta_vec));
    // Equal compact representation atoms have separate implicit binders.
    // The fixed selector in this term cannot fix the unknown vector's port.
    let uncontrolled = bare_rows.frame.open_basis[1]
        .with_checked_expression(cooking.try_cook((&uncontrolled + controlled).as_view())?)?;
    let error = bare_rows
        .frame
        .unit_finite(&uncontrolled, &record.graph, &bare.0)
        .expect_err("an unresolved supplied vector must not enter the physical block");
    ensure!(error.to_string().contains("unresolved Lorentz direction"));
    // Source binding precedes transport: compare the COMPLETE production
    // matrix, rather than the mass closure which annihilates this direction.
    QcdMassReference::certify_finite_tagged_rows_zero(
        "actual bare dimensional finite output reconstructs from its actual WD rows",
        bare_rows
            .production_finite
            .try_sub(&bare_rows.dimensional)?,
    )?;
    let mut bare_four = bare.2.clone();
    bare_four.finite_counterterm = bare_rows.physical_four.clone();
    ensure!(bare_four.pole == bare.2.pole && bare_four.scale_power == bare.2.scale_power);
    let ctx = UVCtx::new(&record.graph, &record.graph.model, &record.settings);
    // The runtime amplitude retains the generated external field ordering.
    // The W references use the amputated Grassmann kernel; project() already
    // proves their source-only ratio, before any new vacuum integration.
    let bare_observed = reference.certify_d09_vacuum_policy_field_emission(
        &record.diagram,
        &ctx,
        &bare,
        &bare_rows.frame,
        &bare_four,
        (&bare_rows.source_weight, &bare_rows.external_ordering),
    )?;
    // Qualify the bare source/observer boundary before the expensive proper
    // source traces. All original proper and aggregate assertions follow.
    let proper =
        reference.certify_actual_overall_mass_input(&record, &record.given, &record.completed)?;
    let local_rows = reference.certify_d09_vacuum_scalar_rows(
        &record,
        &record.given,
        &record.completed,
        &proper,
        D09VacuumSource::ProperLocal,
    )?;
    let finite_rows = reference.certify_d09_vacuum_scalar_rows(
        &record,
        &record.given,
        &record.completed,
        &proper,
        D09VacuumSource::ProperFinite,
    )?;
    for rows in [&local_rows, &finite_rows] {
        ensure!(rows.source_weight == bare_rows.source_weight && rows.omega == bare_rows.omega);
        ensure!(rows.external_ordering == bare_rows.external_ordering);
        ensure!(rows.frame.mass == bare_rows.frame.mass);
        ensure!(rows.frame.lmb == bare_rows.frame.lmb);
        ensure!(rows.frame.external_points == bare_rows.frame.external_points);
        for index in 0..2 {
            ensure!(
                rows.frame.open_basis[index].structure().logical_slots()
                    == bare_rows.frame.open_basis[index]
                        .structure()
                        .logical_slots()
            );
            QcdMassReference::certify_finite_tagged_rows_zero(
                "all D09 source rows retain the same calibrated native open frame",
                QcdMassReference::materialize_finite_external_matrix_coefficients(
                    &rows.frame.open_basis[index].try_sub(&bare_rows.frame.open_basis[index])?,
                )?,
            )?;
        }
    }
    let proper_dimensional = local_rows.dimensional.try_add(&finite_rows.dimensional)?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "actual completed-child dimensional finite output reconstructs from both actual WD sources",
        local_rows.production_finite.try_sub(&proper_dimensional)?,
    )?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "both proper source observers borrow the same production finite output",
        local_rows
            .production_finite
            .try_sub(&finite_rows.production_finite)?,
    )?;
    // The diagnostic counterparts are constructed from ACTUAL W4/C0 rows,
    // independently of original A_D and independently of the reference below.
    let proper_four = local_rows
        .physical_four
        .try_add(&finite_rows.physical_four)?;
    let mut completed_four = proper.2.clone();
    completed_four.finite_counterterm = proper_four;
    ensure!(
        completed_four.pole == proper.2.pole && completed_four.scale_power == proper.2.scale_power
    );
    let proper_observed = reference.certify_d09_vacuum_policy_field_emission(
        &record.diagram,
        &ctx,
        &proper,
        &local_rows.frame,
        &completed_four,
        (&local_rows.source_weight, &local_rows.external_ordering),
    )?;
    let actual = [0usize, 1].map(|dual| {
        &bare_rows.weighted_difference[dual]
            + &local_rows.weighted_difference[dual]
            + &finite_rows.weighted_difference[dual]
    });
    let observed = [0usize, 1].map(|dual| {
        reference
            .model
            .expand_couplings(&(&bare_observed[dual] + &proper_observed[dual]))
    });
    for dual in 0..2 {
        let finite = super::truncate(&super::series(&actual[dual], 0)?, true);
        QcdMassReference::certify_finite_tagged_rows_zero(
            "actual emitted finite coefficient retains the source ordering times the kernel W difference",
            SymbolicTensor::infer_indexed(&observed[dual] - &bare_rows.external_ordering * finite)?,
        )?;
    }
    let phi = &observed[0] + &observed[1];
    let psi = (&observed[1] - &observed[0]) / Atom::num(2);
    QcdMassReference::certify_finite_tagged_rows_zero(
        "observed vacuum-policy difference is invisible to Phi",
        SymbolicTensor::infer_indexed(phi.clone())?,
    )?;
    let nonzero_field = SymbolicTensor::infer_indexed(psi.clone())?;
    let zero_field_rejected = QcdMassReference::certify_finite_tagged_rows_zero(
        "zeroing the observed complementary field carrier is rejected",
        nonzero_field,
    )
    .is_err();
    ensure!(zero_field_rejected);
    // The independent reference enters only AFTER both actual input bindings
    // and the observed emitted difference have been certified.
    let epsilon = Atom::var(GS.dim_epsilon);
    let pi = Atom::var(symbolica::atom::Symbol::PI);
    let r = Atom::num((11, 36)) + pi.pow(2) / Atom::num(12);
    let expected_kernel_field = &bare_rows.omega * &r;
    let expected_field = &bare_rows.external_ordering * &expected_kernel_field;
    for (dual, value) in actual.iter().enumerate() {
        let sign = Atom::num(if dual == 0 { -1 } else { 1 });
        let expected = &sign * &bare_rows.omega * (Atom::num((7, 6)) / &epsilon + &r);
        QcdMassReference::certify_finite_tagged_rows_zero(
            "actual complete D09 Laurent vacuum-policy residual matches independent reference",
            SymbolicTensor::infer_indexed(value - expected)?,
        )?;
    }
    QcdMassReference::certify_finite_tagged_rows_zero(
        "observed complementary finite field coefficient matches independent reference",
        SymbolicTensor::infer_indexed(&psi - &expected_field)?,
    )?;
    let wrong_closure = SymbolicTensor::infer_indexed(&phi - &expected_field)?;
    let phi_substitution_rejected = QcdMassReference::certify_finite_tagged_rows_zero(
        "Phi cannot replace the complementary field diagnostic",
        wrong_closure,
    )
    .is_err();
    ensure!(phi_substitution_rejected);
    crate::debug_tags!(#uv, #integrated, #cff, #summary;
        stage = "native_d09_emitted_vacuum_policy_field", original_sources = 3,
        actual_bare_and_completed_dimensional_carriers_bound = true,
        actual_four_d_counterparts_independently_reconstructed = true,
        file.actual_laurent = ?actual.iter().map(|row| row.to_canonical_string()).collect::<Vec<_>>(),
        file.actual_Phi = %phi.to_canonical_string(), file.actual_Psi = %psi.to_canonical_string(),
        file.reference_Psi = %expected_field.to_canonical_string(),
        file.reference_kernel_Psi = %expected_kernel_field.to_canonical_string(),
        file.source_external_ordering = %bare_rows.external_ordering.to_canonical_string(),
        unit_rest = true, unit_scales = true, emitted_vacuum_policy_field = true,
        full_unregulated_R2 = false, new_production_R2_term = false,
        "Actual native W coefficients survive the observed finalizer in the independently known field direction"
    );
    Ok(())
}

#[test]
fn native_d09_proper_local_vacuum_scalar_traces_match_independent_tadpoles() -> Result<()> {
    QcdMassReference::certify_d09_vacuum_scalar_reference(D09VacuumSource::ProperLocal)
}

#[test]
fn native_d09_bare_vacuum_scalar_traces_match_independent_connected_poles() -> Result<()> {
    QcdMassReference::certify_d09_vacuum_scalar_reference(D09VacuumSource::Bare)
}

#[test]
fn native_d09_proper_finite_vacuum_scalar_traces_keep_child_epsilon() -> Result<()> {
    QcdMassReference::certify_d09_vacuum_scalar_reference(D09VacuumSource::ProperFinite)
}

impl QcdMassReference {
    fn certify_d09_vacuum_scalar_reference(source: D09VacuumSource) -> Result<()> {
        use crate::uv::{
            ApproximationType, Spinney,
            approx::{
                ForestNodeLike, Rooted,
                local_4d::{Full4dCts, Local4dCts},
            },
            hedge_poset::OwnedForestNode,
        };
        crate::initialisation::test_initialise()?;
        let reference = QcdMassReference::new()?;
        let mut records = reference.certify_os_quark_proper_prefix()?;
        ensure!(records.len() == 1);
        let record = records
            .remove(&("D09".to_owned(), vec![2, 3]))
            .ok_or_else(|| eyre!("the actual native D09 OS child is missing"))?;
        ensure!(record.given.renormalization_scheme() == ApproximationType::OS);
        ensure!(record.diagram.loop_momentum_basis().loop_edges.len() == 2);
        let root = OwnedForestNode {
            spinney: Spinney::empty(&record.graph),
            topo_order: 0,
        };
        let bare = Full4dCts::from_factorized_local(&Local4dCts::root());
        let (given, prefix) = if source == D09VacuumSource::Bare {
            (&root, &bare)
        } else {
            (&record.given, &record.completed)
        };
        let actual = reference.certify_actual_overall_mass_input(&record, given, prefix)?;
        reference
            .certify_d09_vacuum_scalar_rows(&record, given, prefix, &actual, source)
            .map(|_| ())
    }

    fn certify_d09_vacuum_scalar_rows(
        &self,
        record: &CertifiedMassPrefix,
        given: &crate::uv::hedge_poset::OwnedForestNode,
        prefix: &crate::uv::approx::local_4d::Full4dCts,
        actual: &(
            crate::uv::hedge_poset::OwnedForestNode,
            crate::uv::approx::local_4d::Local4dCts,
            super::IntegratedCts,
        ),
        source: D09VacuumSource,
    ) -> Result<D09VacuumScalarRows> {
        use crate::{utils::GS, uv::approx::UVCtx};
        let reference = self;
        let source_projection = reference.project(&record.diagram)?;
        let source_weight =
            &source_projection.kernel_weight * &source_projection.numerator_prefactor;
        ensure!(!source_weight.is_zero());
        let (traces, frame) = reference.certify_d09_vacuum_angular_traces(
            &record.diagram,
            &UVCtx::new(&record.graph, &record.graph.model, &record.settings),
            given,
            prefix,
            actual,
            source,
        )?;

        let epsilon = Atom::var(GS.dim_epsilon);
        let pi = Atom::var(symbolica::atom::Symbol::PI);
        let coupling = Atom::var(symbol!("UFO::G"));
        let a = coupling.pow(2) / (Atom::num(16) * pi.pow(2));
        let cf = Atom::num((4, 3));
        // The source-only i*G^4*CF^2 calibration precedes the new integrals.
        // Two MSbar Minkowski measures supply i^2/(16*pi^2)^2. The signed
        // child -U is already inside c(D); run(W) receives no further forest sign.
        let prefactor = -Atom::i() * a.pow(2) * cf.pow(2) * &source_weight;
        ensure!(!prefactor.is_zero());
        // Independent radial Gamma recurrence: T2=e^(gamma_E*eps)*Gamma(eps)
        // at m=MUV=muR=1, T2^2=eps^-2+pi^2/6+O(eps).
        // Native source absorption and numerator-only lift are prerequisites.
        // c(D)=(D-2)^2(D+2)/24=1-7eps/3+5eps^2/3-eps^3/3.
        let q = pi.pow(2) / Atom::num(12);
        let (dimensional, physical_four) = match source {
            D09VacuumSource::ProperLocal => (
                epsilon.pow(-2) - Atom::num((7, 3)) / &epsilon
                    + Atom::num((5, 3))
                    + pi.pow(2) / Atom::num(6),
                epsilon.pow(-2) + pi.pow(2) / Atom::num(6),
            ),
            D09VacuumSource::Bare => (
                -Atom::num((1, 2)) * epsilon.pow(-2) + Atom::num((13, 12)) / &epsilon,
                -Atom::num((1, 2)) * epsilon.pow(-2) - Atom::num((1, 12)) / &epsilon,
            ),
            D09VacuumSource::ProperFinite => (
                epsilon.pow(-1) - Atom::num((4, 3)) - &q,
                epsilon.pow(-1) - &q,
            ),
        };
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut records = Vec::new();
        let mut actual_series = Vec::new();
        let mut native_series = Vec::new();
        for (index, trace) in traces.iter().enumerate() {
            ensure!(trace.is_scalar());
            // The finite-row owner expects completed scalar coefficients. Use the
            // existing SU(3) invariant owner before comparing to numerical CF=4/3.
            let normalized = trace.simplify_algebra(&AlgebraSettings {
                color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
                ..Default::default()
            })?;
            ensure!(
                !normalized.expression().contains_symbol(CS.cas)
                    && !normalized.expression().contains_symbol(CS.idx)
            );
            let native_unit = normalized
                .expression()
                .replace(symbol!("UFO::MT"))
                .with(1)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.m_uv_vacuum)
                .with(1)
                .replace(GS.mu_r_sq)
                .with(1);
            native_series.push(
                super::series(
                    &super::Integrated::dimensionally_regularized(&native_unit),
                    0,
                )?
                .to_atom(),
            );
            let weighted = reference
                .model
                .expand_couplings(&(source_weight.clone() * normalized.expression()));
            // Unit scales are imposed after the actual operators and integration.
            let unit = weighted
                .replace(symbol!("UFO::MT"))
                .with(1)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.m_uv_vacuum)
                .with(1)
                .replace(GS.mu_r_sq)
                .with(1);
            let series =
                super::series(&super::Integrated::dimensionally_regularized(&unit), 0)?.to_atom();
            let sign = Atom::num(if index % 2 == 0 { -1 } else { 1 });
            let expected = &sign
                * &prefactor
                * if index < 2 {
                    dimensional.clone()
                } else if index < 4 {
                    physical_four.clone()
                } else {
                    epsilon.pow(-1)
                };
            let compared = if source == D09VacuumSource::Bare {
                super::truncate(&super::series(&series, 0)?, false)
            } else {
                series.clone()
            };
            let difference = trace
                .with_checked_expression(cooking.try_cook((&compared - &expected).as_view())?)?;
            QcdMassReference::certify_finite_tagged_rows_zero(
                "absolute actual D09 source-selected vacuum scalar Laurent trace",
                difference,
            )?;
            // Direct nonzero leading coefficient: a zero trace cannot pass merely
            // through a permissive scalar-zero normalization or a blind mass Phi.
            let rows = series
                .coefficient_list_exact(std::slice::from_ref(&epsilon))
                .ok_or_else(|| {
                    eyre!("D09 absolute scalar trace exceeds the unchanged Laurent scope")
                })?;
            let leading_power = if source == D09VacuumSource::ProperFinite {
                -1
            } else {
                -2
            };
            let leading = rows
                .iter()
                .find(|(monomial, _)| *monomial == epsilon.pow(leading_power))
                .map(|(_, coefficient)| coefficient)
                .ok_or_else(|| {
                    eyre!("D09 scalar trace lost its independently nonzero double pole")
                })?;
            ensure!(!leading.is_zero());
            let leading_weight = if source == D09VacuumSource::Bare {
                Atom::num((-1, 2))
            } else {
                Atom::one()
            };
            ensure!(
                (leading - &sign * &prefactor * &leading_weight)
                    .together()
                    .is_zero()
            );
            actual_series.push(series.clone());
            records.push(serde_json::json!({
                "physical_four": index >= 2,
                "dual": if index % 2 == 0 { "identity/12" } else { "gamma0/12" },
                "actual": series.to_canonical_string(),
                "reference": expected.to_canonical_string(),
                "double_pole": leading.to_canonical_string(),
            }));
        }
        ensure!(
            records.len()
                == if source == D09VacuumSource::ProperFinite {
                    6
                } else {
                    4
                }
        );
        if source == D09VacuumSource::Bare {
            // The connected finite master cancels from the angular difference.
            // No finite master value is read or fitted to define this reference.
            for dual in 0..2 {
                let sign = Atom::num(if dual == 0 { -1 } else { 1 });
                let difference = &actual_series[dual + 2] - &actual_series[dual];
                let expected =
                    &sign * &prefactor * (-Atom::num((7, 6)) / &epsilon + Atom::num((23, 36)));
                let tensor = traces[dual].with_checked_expression(
                    cooking.try_cook((&difference - &expected).as_view())?,
                )?;
                Self::certify_finite_tagged_rows_zero(
                    "actual bare D09 angular difference needs no finite sunset master",
                    tensor,
                )?;
                // Independent epsilon-two loss is finite 5/36. Dropping it changes
                // 23/36 to 7/9, so this is a nonzero depth witness.
                let insufficient =
                    &sign * &prefactor * (-Atom::num((7, 6)) / &epsilon + Atom::num((7, 9)));
                let wrong = traces[dual].with_checked_expression(
                    cooking.try_cook((&difference - insufficient).as_view())?,
                )?;
                ensure!(
                    Self::certify_finite_tagged_rows_zero(
                        "actual bare angular difference rejects omitted epsilon-two moment",
                        wrong.clone()
                    )
                    .is_err()
                );
            }
        }
        if source == D09VacuumSource::ProperFinite {
            for dual in 0..2 {
                let sign = Atom::num(if dual == 0 { -1 } else { 1 });
                let loss = &actual_series[dual + 4] - &actual_series[dual + 2];
                let expected = &sign * &prefactor * &q;
                let tensor = traces[dual]
                    .with_checked_expression(cooking.try_cook((&loss - &expected).as_view())?)?;
                Self::certify_finite_tagged_rows_zero(
                    "actual frozen child epsilon-one feeds the outer simple pole",
                    tensor,
                )?;
                ensure!(!expected.is_zero());
                let omitted =
                    traces[dual].with_checked_expression(cooking.try_cook(loss.as_view())?)?;
                ensure!(
                    Self::certify_finite_tagged_rows_zero(
                        "actual frozen child outer pole rejects omitted epsilon-one coefficient",
                        omitted.clone()
                    )
                    .is_err()
                );
            }
        }
        crate::debug_tags!(#uv, #integrated, #summary;
            stage = "native_d09_vacuum_independent_scalar_reference", file.source = ?source,
            file.records = %serde_json::to_string(&records)?,
            file.source_weight = %source_weight.to_canonical_string(),
            original_terms = if source == D09VacuumSource::ProperLocal { 4 } else { 2 },
            active_loops = if source == D09VacuumSource::ProperFinite { 1 } else { 2 },
            frozen_loops = if source == D09VacuumSource::ProperFinite { 1 } else { 0 },
            absolute_poles_proved = true, absolute_finite_proved = source != D09VacuumSource::Bare,
            bare_included = source == D09VacuumSource::Bare, finite_child_included = source == D09VacuumSource::ProperFinite,
            full_open_matrix_proved = true, emitted_R2_proved = false,
            "Actual D09 source-selected W scalar-dual Laurent traces match independent vacuum references"
        );
        let physical_start = if source == D09VacuumSource::ProperFinite {
            4
        } else {
            2
        };
        let mut dimensional = frame.open_basis[0].with_checked_expression(Atom::Zero)?;
        let mut physical_four = dimensional.clone();
        for dual in 0..2 {
            let d = super::truncate(&super::series(&native_series[dual], 0)?, true);
            let four = super::truncate(
                &super::series(&native_series[physical_start + dual], 0)?,
                true,
            );
            dimensional = dimensional.try_add(
                &frame.open_basis[dual]
                    .with_checked_expression(-d * frame.open_basis[dual].expression())?,
            )?;
            physical_four = physical_four.try_add(
                &frame.open_basis[dual]
                    .with_checked_expression(-four * frame.open_basis[dual].expression())?,
            )?;
        }
        let dimensional = frame.unit_finite(&dimensional, &record.graph, &actual.0)?;
        let physical_four = frame.unit_finite(&physical_four, &record.graph, &actual.0)?;
        let production_finite = frame.unit_finite(
            &actual.2.physical_finite_counterterm()?,
            &record.graph,
            &actual.0,
        )?;
        let weighted_difference =
            [0usize, 1].map(|dual| &actual_series[physical_start + dual] - &actual_series[dual]);
        Ok(D09VacuumScalarRows {
            frame,
            dimensional,
            physical_four,
            production_finite,
            weighted_difference,
            source_weight,
            external_ordering: source_projection.external_ordering,
            omega: prefactor,
        })
    }
}

impl QcdMassReference {
    /// Complete finite external matrices while keeping scalar radial D words opaque.
    pub(in crate::uv) fn materialize_finite_external_matrix_coefficients(
        input: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        ensure!(input.rank() == 4);
        ensure!(!input.expression().contains_symbol(crate::utils::GS.den));
        let filter = spenso::shadowing::TensorCollectFilter::Reps([
            spenso::structure::representation::LibraryRep::from(Minkowski {}),
        ]);
        let partial: SymbolicTensor<PartialStructure> = input.clone().into();
        let rows = partial.coefficient_list(filter)?;
        let recomposed = input.with_checked_expression(Atom::add_many(
            rows.iter()
                .map(|(radial, matrix)| radial.expression() * matrix.expression()),
        ))?;
        Self::certify_finite_tagged_rows_zero(
            "radial decomposition retains its complete finite external matrix coefficients",
            input.try_sub(&recomposed)?,
        )?;
        let mut normalized = input.with_checked_expression(Atom::Zero)?;
        for (radial, matrix) in rows {
            if radial.expression().is_zero() || matrix.expression().is_zero() {
                continue;
            }
            ensure!(
                radial.is_scalar(),
                "a radial spectator retained a free Lorentz port"
            );
            for head in [
                crate::utils::GS.emr_mom,
                crate::utils::GS.loop_mom,
                crate::utils::GS.den,
            ] {
                ensure!(
                    !matrix.expression().contains_symbol(head),
                    "a finite external matrix retained an active radial carrier"
                );
            }
            let matrix = input.with_checked_expression(matrix.expression().clone())?;
            ensure!(matrix.rank() == 4);
            let matrix = Self::materialize_native_four_component_numerator(&matrix)?;
            normalized = normalized.try_add(
                &matrix.with_checked_expression(matrix.expression() * radial.expression())?,
            )?;
        }
        ensure!(normalized.structure().logical_slots() == input.structure().logical_slots());
        Ok(normalized)
    }

    /// Materialize only the admitted finite numerator, never a causal shell.
    pub(in crate::uv) fn materialize_native_four_component_numerator(
        input: &SymbolicTensor<IndexedStructure>,
    ) -> Result<SymbolicTensor<IndexedStructure>> {
        use idenso::representations::Bispinor;
        use idenso::shorthands::UndoShorthands;
        use spenso::structure::representation::LibraryRep;
        use spenso::{
            iterators::IteratableTensor,
            network::library::symbolic::ETS,
            network::{ExecutionResult, Sequential, SmallestDegree},
            structure::{abstract_index::AIND_SYMBOLS, slot::IsAbstractSlot},
        };
        let expected = input.structure().logical_slots();
        ensure!(expected.len() == 4);
        let mut spin = 0usize;
        let mut color = 0usize;
        for slot in &expected {
            if slot.rep_name().base() == LibraryRep::from(Bispinor {}) {
                ensure!(slot.dim().to_symbolic() == Atom::num(4));
                spin += 1;
            } else if slot.rep_name().base() == LibraryRep::from(ColorFundamental {}) {
                ensure!(slot.dim().to_symbolic() == Atom::num(3));
                color += 1;
            } else {
                return Err(eyre!(
                    "D09 finite numerator has an unexpected external representation"
                ));
            }
        }
        ensure!(spin == 2 && color == 2);
        ensure!(!input.expression().contains_symbol(crate::utils::GS.den));
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        // Collect complete native tensor words before finite execution. The
        // original scalar coefficients (including the child epsilon/log jet)
        // stay opaque spectators, rather than being copied into every component.
        // Tensors includes compact dots: active physical momentum dependence
        // must remain in the finite word, not in a D-dimensional scalar weight.
        let partial: SymbolicTensor<PartialStructure> = input.clone().into();
        let filter = spenso::shadowing::TensorCollectFilter::<0>::Tensors;
        let rows = partial.coefficient_list(filter)?;
        let recomposed = input.with_checked_expression(Atom::add_many(
            rows.iter()
                .map(|(word, coefficient)| word.expression() * coefficient.expression()),
        ))?;
        Self::certify_finite_tagged_rows_zero(
            "finite component-word decomposition retains the original numerator",
            input.try_sub(&recomposed)?,
        )?;
        let mut materialized = input.with_checked_expression(Atom::Zero)?;
        let mut cache = crate::uv::overlap_numerators::NumeratorJetCache::default();
        for (word, coefficient) in rows {
            if word.expression().is_zero() || coefficient.expression().is_zero() {
                continue;
            }
            ensure!(coefficient.is_scalar());
            for head in [
                crate::utils::GS.emr_mom,
                crate::utils::GS.loop_mom,
                crate::utils::GS.den,
            ] {
                ensure!(!coefficient.expression().contains_symbol(head));
            }
            let word = input.with_checked_expression(word.expression().clone())?;
            let source = cooking.uncook(word.expression().as_view());
            let explicit = source.undo_dots::<crate::numerator::aind::Aind>()?;
            let mut network = explicit.parse_into_net()?;
            let library = crate::utils::TENSORLIB.read().unwrap();
            network.execute::<Sequential, SmallestDegree, _, _, _>(
                &*library,
                &*crate::utils::FUN_LIB,
            )?;
            let expression = match network.result_tensor(&*library)? {
                ExecutionResult::Zero => Atom::Zero,
                ExecutionResult::One => {
                    return Err(eyre!(
                        "a nonzero complete native quark word became a scalar one"
                    ));
                }
                ExecutionResult::Val(tensor) => {
                    ensure!(tensor.size()? <= 4 * 4 * 3 * 3);
                    let ports = tensor
                        .external_structure_iter()
                        .map(|slot| slot.to_atom())
                        .collect::<Vec<_>>();
                    // The explicit library returned native index payloads; the
                    // retained IndexedStructure carries their reversible encoding.
                    // Compare the same representation, including variance and exact
                    // named payload, and certify both encoding directions.
                    let source_ports = expected
                        .iter()
                        .map(|slot| -> Result<Atom> {
                            let encoded = slot.to_atom();
                            let native = cooking.uncook(encoded.as_view());
                            ensure!(cooking.try_cook(native.as_view())? == encoded);
                            Ok(native)
                        })
                        .collect::<Result<Vec<_>>>()?;
                    ensure!(ports.len() == source_ports.len());
                    for native in &ports {
                        ensure!(
                            cooking.uncook(cooking.try_cook(native.as_view())?.as_view())
                                == *native
                        );
                    }
                    ensure!(
                        ports.iter().cloned().collect::<BTreeSet<_>>()
                            == source_ports.iter().cloned().collect::<BTreeSet<_>>()
                    );
                    let mut entries = Vec::new();
                    for (index, value) in tensor.iter_flat() {
                        if value.is_zero() {
                            continue;
                        }
                        let index = tensor.expanded_index(index)?;
                        let basis = Atom::mul_many(ports.iter().zip(index.iter()).map(
                            |(port, coordinate)| {
                                function!(
                                    ETS.delta,
                                    function!(AIND_SYMBOLS.cind, *coordinate),
                                    port
                                )
                            },
                        ));
                        entries.push(value * basis.as_view());
                    }
                    Atom::add_many(entries)
                }
            };
            drop(library);
            let result = word.with_checked_expression(cooking.try_cook(expression.as_view())?)?;
            ensure!(result.structure().logical_slots() == input.structure().logical_slots());
            // Certify this bounded tensor word, with all scalar spectators absent.
            // The existing component/byte caps and exact compact zero owner remain
            // unchanged. Multiplication by its original coefficient preserves this
            // proved word equality without materializing the complete graph sum.
            let difference = cooking.uncook(word.try_sub(&result)?.expression().as_view());
            ensure!(cache.canonicalize(&difference)?.is_zero);
            materialized =
                materialized
                    .try_add(&result.with_checked_expression(
                        result.expression() * coefficient.expression(),
                    )?)?;
        }
        ensure!(materialized.structure().logical_slots() == input.structure().logical_slots());
        Ok(materialized)
    }
}

#[test]
fn native_d09_finite_matrix_completeness_keeps_dimensional_radial_spectators() -> Result<()> {
    use crate::utils::GS;
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
    use spenso::structure::abstract_index::AIND_SYMBOLS;
    crate::initialisation::test_initialise()?;
    let row = Atom::var(symbol!("native_radial_matrix_regression::row"));
    let column = Atom::var(symbol!("native_radial_matrix_regression::column"));
    let spin_row = idenso::bis!(4, row.clone());
    let spin_column = idenso::bis!(4, column.clone());
    let color_row = ColorFundamental {}.new_rep(3).to_symbolic([row]);
    let color_column = function!(
        AIND_SYMBOLS.dind,
        ColorFundamental {}.new_rep(3).to_symbolic([column])
    );
    let identity = function!(ETS.metric, &spin_row, &spin_column)
        * function!(ETS.metric, &color_row, &color_column);
    let component =
        |index: i32, port: &Atom| function!(ETS.delta, function!(AIND_SYMBOLS.cind, index), port);
    let units = (0..4)
        .flat_map(|spin| (0..3).map(move |color| (spin, color)))
        .map(|(spin, color)| {
            component(spin, &spin_row)
                * component(spin, &spin_column)
                * component(color, &color_row)
                * component(color, &color_column)
        })
        .collect::<Vec<_>>();
    let complete = Atom::add_many(units.iter());
    let rep = (Minkowski {}).to_symbolic([GS.dim]);
    let radial_factors = [
        function!(
            SPENSO_TAG.dot,
            function!(GS.loop_mom, 1, &rep),
            function!(GS.loop_mom, 2, &rep)
        ),
        function!(
            SPENSO_TAG.dot,
            function!(GS.loop_mom, 2, &rep),
            function!(GS.loop_mom, 2, &rep)
        ),
    ];
    let radial = Atom::mul_many(&radial_factors);
    let a = Atom::var(symbol!("native_radial_matrix_regression::a"; Scalar));
    let b = Atom::var(symbol!("native_radial_matrix_regression::b"; Scalar));
    let spectator = (&a + &b).pow(11);
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let input = SymbolicTensor::infer_indexed(
        cooking.try_cook(((&identity - &complete) * &radial * &spectator).as_view())?,
    )?;
    let before = input.clone();
    let unmaterialized = QcdMassReference::certify_finite_tagged_rows_zero(
        "unmaterialized finite matrix basis is not an independent word basis",
        input.clone(),
    );
    ensure!(unmaterialized.is_err());
    let normalized = QcdMassReference::materialize_finite_external_matrix_coefficients(&input)?;
    QcdMassReference::certify_finite_tagged_rows_zero(
        "exact finite I completeness with independent D radial moments",
        normalized,
    )?;
    let one = SymbolicTensor::infer_indexed(
        cooking.try_cook((&identity * &radial * &spectator).as_view())?,
    )?;
    let preserved = QcdMassReference::materialize_finite_external_matrix_coefficients(&one)?;
    for factor in &radial_factors {
        ensure!(preserved.expression().contains(factor.as_view()));
    }
    ensure!(preserved.expression().contains(spectator.as_view()));
    ensure!(preserved.structure().logical_slots() == one.structure().logical_slots());
    let wrong = SymbolicTensor::infer_indexed(
        cooking.try_cook(((&identity - &complete + &units[0]) * &radial * &spectator).as_view())?,
    )?;
    let wrong = QcdMassReference::materialize_finite_external_matrix_coefficients(&wrong)?;
    let perturbed = QcdMassReference::certify_finite_tagged_rows_zero(
        "one perturbed native matrix coordinate must remain nonzero",
        wrong,
    );
    ensure!(perturbed.is_err());
    ensure!(input.same_value(&before));
    Ok(())
}

/// A common positive scalar gluon mass is an analytic IR rearrangement only.
/// It changes physical denominators, never the source Feynman-gauge numerator.
#[test]
fn native_d09_positive_gluon_mass_preserves_all_three_vacuum_inputs() -> Result<()> {
    use crate::{
        graph::edge::EdgeMass,
        model::UFOSymbol,
        uv::{
            ApproximationType, Spinney, UltravioletGraph,
            approx::{
                ForestNodeLike, Rooted, UVCtx,
                local_4d::{FourDSector, Full4dCts, Local4dCts, uv_limit},
            },
            hedge_poset::OwnedForestNode,
            uv_graph::UVE,
        },
    };
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{Inclusion, SubSetOps},
    };

    crate::initialisation::test_initialise()?;
    let reference = QcdMassReference::new()?;
    let mut records = reference.certify_os_quark_proper_prefix()?;
    ensure!(records.len() == 1);
    let record = records
        .remove(&("D09".to_owned(), vec![2, 3]))
        .ok_or_else(|| eyre!("native D09 OS child missing"))?;
    let graph = &record.graph;
    let child = &record.given;
    let settings = &record.settings;
    ensure!(settings.generate_integrated && !settings.add_marker);
    ensure!(settings.on_shell_velocity()? == [Atom::one(), Atom::Zero, Atom::Zero, Atom::Zero]);
    // MH is already a validated native scalar parameter. It is absent from
    // this QCD source, and remains algebraically independent of MT and mUV.
    // The analytic reference below takes its value lambda strictly positive.
    let _parameter = graph.model.parameter("MH")?;
    let lambda = Atom::from(UFOSymbol::from("MH"));
    let lambda_symbol = UFOSymbol::from("MH").0;
    let root = OwnedForestNode {
        spinney: Spinney::empty(graph),
        topo_order: 0,
    };
    let empty = Full4dCts::from_factorized_local(&Local4dCts::root());
    let uncut = graph.full_filter().subtract(&graph.initial_state_cut.left);
    let mut parents = graph
        .classified_spinneys(&uncut, &graph.model, settings, &graph.loop_momentum_basis)?
        .into_iter()
        .filter(|component| {
            component.n_components() == 1
                && component.dod == 1
                && component.renormalization_scheme == ApproximationType::OS
                && graph.n_loops(component.filter()) == 2
                && component.filter().includes(child.subgraph())
        })
        .collect::<Vec<_>>();
    ensure!(parents.len() == 1);
    let parent_spinney = parents.remove(0);
    let mut regulated = graph.clone();
    for edge in [EdgeIndex(3), EdgeIndex(6)] {
        ensure!(graph[edge].particle_pdg_code(&graph.model) == Some(21));
        ensure!(graph[edge].mass_atom(&graph.model).is_zero());
        let particle = graph[edge].particle.particle();
        let numerator = graph[edge].num.clone();
        let target = &mut regulated.underlying[edge];
        target.particle = target.particle.clone().override_mass(Some(lambda.clone()));
        target.mass = EdgeMass::from_atom(lambda.clone(), &graph.model, &graph.param_builder)?;
        ensure!(target.particle.particle() == particle);
        ensure!(
            target.num == numerator,
            "the IR regulator cannot change the gauge numerator"
        );
        ensure!(target.mass_atom(&graph.model) == lambda);
    }
    ensure!(regulated.loop_momentum_basis == graph.loop_momentum_basis);
    let ctx = UVCtx::new(graph, &graph.model, settings);
    let regulated_ctx = UVCtx::new(&regulated, &regulated.model, settings);
    let parent = OwnedForestNode {
        spinney: parent_spinney.clone(),
        topo_order: 1,
    };
    let source = graph
        .numerator(parent.subgraph(), root.subgraph())?
        .to_d_dim(crate::utils::GS.dim)?
        .color_simplify()?
        .state
        .expr;
    let regulated_source = regulated
        .numerator(parent.subgraph(), root.subgraph())?
        .to_d_dim(crate::utils::GS.dim)?
        .color_simplify()?
        .state
        .expr;
    ensure!(source.same_value(&regulated_source));
    ensure!(!source.expression().contains_symbol(lambda_symbol));
    let original_bare = source.with_checked_expression(
        source.expression() * graph.denominator(parent.subgraph(), &graph.model, |_| -1),
    )?;
    let regulated_bare = regulated_source.with_checked_expression(
        regulated_source.expression()
            * regulated.denominator(parent.subgraph(), &regulated.model, |_| -1),
    )?;
    let original_terms = FourDSector::new(
        original_bare,
        vec![(
            parent.subgraph().clone(),
            parent.subgraph().clone(),
            parent.lmb().clone(),
        )],
        Vec::new(),
    )
    .physical_terms()?;
    let regulated_terms = FourDSector::new(
        regulated_bare.clone(),
        vec![(
            parent.subgraph().clone(),
            parent.subgraph().clone(),
            parent.lmb().clone(),
        )],
        Vec::new(),
    )
    .physical_terms()?;
    ensure!(original_terms.len() == 1 && regulated_terms.len() == 1);
    ensure!(
        original_terms[0]
            .numerator
            .same_value(&regulated_terms[0].numerator)
    );
    ensure!(original_terms[0].denominators.len() == regulated_terms[0].denominators.len());
    for (original, regulator) in original_terms[0]
        .denominators
        .iter()
        .zip(&regulated_terms[0].denominators)
    {
        ensure!(original.source_edge == regulator.source_edge);
        ensure!(original.momentum == regulator.momentum);
        let shifted = [EdgeIndex(3), EdgeIndex(6)].contains(&original.source_edge);
        let shift = if shifted { lambda.pow(2) } else { Atom::Zero };
        ensure!(regulator.mass_squared == &original.mass_squared + &shift);
        // Denominator-only quadratic arithmetic; the graph numerator is untouched.
        ensure!(
            (&regulator.full_expr - &original.full_expr + &shift)
                .expand()
                .is_zero()
        );
    }
    let branches = regulated.on_shell_branches(
        &regulated.model,
        &regulated_bare,
        parent.subgraph(),
        parent.lmb(),
        &settings.on_shell_velocity()?,
    )?;
    for branch in branches {
        let terms = FourDSector::new(
            branch,
            vec![(
                parent.subgraph().clone(),
                parent.subgraph().clone(),
                parent.lmb().clone(),
            )],
            Vec::new(),
        )
        .physical_terms()?;
        ensure!(!terms.is_empty());
        for term in terms {
            for denominator in term.denominators {
                if [EdgeIndex(3), EdgeIndex(6)].contains(&denominator.source_edge) {
                    ensure!(denominator.mass_squared == lambda.pow(2));
                }
            }
        }
    }
    // O allocates a fresh contracted spinor port at every invocation. Compare
    // through the existing native index owner after binding every prescribed
    // denominator occurrence; no momentum, mass or provenance is normalized away.
    let denominator_inventory = |value: &SymbolicTensor<IndexedStructure>| {
        let sector = FourDSector::new(value.clone(), Vec::new(), Vec::new());
        let mut inventory = sector
            .physical_terms()?
            .into_iter()
            .map(|term| {
                let mut occurrences = term
                    .denominators
                    .into_iter()
                    .map(|denominator| {
                        (
                            denominator.source_edge.0,
                            denominator.momentum,
                            denominator.mass_squared,
                            denominator.full_expr,
                        )
                    })
                    .collect::<Vec<_>>();
                // Powers remain repeated occurrences with their source owners.
                occurrences.sort();
                occurrences
            })
            .collect::<Vec<_>>();
        inventory.sort();
        Ok::<_, color_eyre::eyre::Report>(inventory)
    };
    // The only freshly allocated labels admitted here are the internal bis4
    // projector ports of this rest-frame two-component forest. Original source
    // ports and every prescribed denominator are held literally fixed.
    let exact_vacuum_equal = |left: &SymbolicTensor<IndexedStructure>,
                              right: &SymbolicTensor<IndexedStructure>|
     -> Result<()> {
        use crate::numerator::aind::Aind;
        use idenso::representations::Bispinor;
        use spenso::structure::{
            representation::{LibrarySlot, RepName},
            slot::IsAbstractSlot,
        };
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        ensure!(left.structure().logical_slots() == right.structure().logical_slots());
        ensure!(!left.expression().contains_symbol(lambda_symbol));
        ensure!(!right.expression().contains_symbol(lambda_symbol));
        let original_denominators = denominator_inventory(left)?;
        ensure!(original_denominators == denominator_inventory(right)?);
        let fresh_ports = |expression: &Atom| -> Result<BTreeSet<Aind>> {
            let mut indices = BTreeSet::new();
            let mut invalid = false;
            expression.visitor(&mut |view| {
                if let Ok(slot) = LibrarySlot::<Aind>::try_from(view) {
                    if matches!(slot.aind, Aind::Dummy(_)) {
                        invalid |= slot.rep() != Bispinor {}.new_rep(4).to_lib();
                        indices.insert(slot.aind);
                    }
                    return false;
                }
                true
            });
            ensure!(
                !invalid,
                "D09 fresh projector labels must be internal bis4 ports"
            );
            ensure!(
                indices.len() <= 2,
                "D09 has at most two rest-frame OS projector ports"
            );
            Ok(indices)
        };
        let raw_left = cooking.uncook(left.expression().as_view());
        let raw_right = cooking.uncook(right.expression().as_view());
        let left_indices = fresh_ports(&raw_left)?.into_iter().collect::<Vec<_>>();
        let right_indices = fresh_ports(&raw_right)?.into_iter().collect::<Vec<_>>();
        ensure!(left_indices.len() == right_indices.len());
        for tensor in [left, right] {
            for slot in tensor.structure().logical_slots() {
                let raw = cooking.uncook(slot.to_atom().as_view());
                ensure!(
                    fresh_ports(&raw)?.is_empty(),
                    "external ports are immutable"
                );
            }
        }
        for term in &original_denominators {
            for (edge, momentum, mass, polynomial) in term {
                let complete = function!(crate::utils::GS.den, *edge, momentum, mass, polynomial);
                ensure!(fresh_ports(&cooking.uncook(complete.as_view()))?.is_empty());
            }
        }
        let mut candidates = vec![right_indices.clone()];
        if right_indices.len() == 2 {
            candidates.push(vec![right_indices[1], right_indices[0]]);
        }
        for candidate in candidates {
            let map = left_indices
                .iter()
                .copied()
                .zip(candidate.iter().copied())
                .collect::<BTreeMap<_, _>>();
            let inverse = candidate
                .iter()
                .copied()
                .zip(left_indices.iter().copied())
                .collect::<BTreeMap<_, _>>();
            let rename = |source: &Atom, map: &BTreeMap<Aind, Aind>| {
                source.replace_map(|view, _, output| {
                    if let Ok(slot) = LibrarySlot::<Aind>::try_from(view)
                        && let Some(index) = map.get(&slot.aind)
                    {
                        **output = slot.rep().slot::<Aind, _>(*index).to_atom();
                    }
                })
            };
            let renamed = rename(&raw_left, &map);
            ensure!(rename(&renamed, &inverse) == raw_left);
            let native = left.with_checked_expression(cooking.try_cook(renamed.as_view())?)?;
            ensure!(native.structure().logical_slots() == left.structure().logical_slots());
            ensure!(denominator_inventory(&native)? == original_denominators);
            // Literal equality of the COMPLETE W is the certificate. Only its
            // two witnessed contracted spinor labels may have changed.
            if native.same_value(right) {
                return Ok(());
            }
        }
        Err(eyre!(
            "D09 complete W has no exact fresh-projector-port bijection"
        ))
    };
    let original_child = uv_limit(&empty, &ctx, child, &root, child, &root)?;
    let regulated_child = uv_limit(&empty, &regulated_ctx, child, &root, child, &root)?;
    ensure!(original_child.active_sectors().len() == 1);
    ensure!(regulated_child.active_sectors().len() == 1);
    ensure!(
        regulated_child
            .atom()
            .expression()
            .contains_symbol(lambda_symbol)
    );
    ensure!(!regulated_child.atom().same_value(original_child.atom()));
    exact_vacuum_equal(
        original_child.active_sectors()[0].integration_atom(),
        regulated_child.active_sectors()[0].integration_atom(),
    )?;
    let child_w = original_child.active_sectors()[0].integration_atom();
    let wrong_weight = child_w.with_checked_expression(child_w.expression() * Atom::num(2))?;
    ensure!(
        exact_vacuum_equal(child_w, &wrong_weight).is_err(),
        "projector-port renaming must retain the complete scalar weight"
    );
    let backend = vakint::VakintSettings {
        number_of_terms_in_epsilon_expansion: 3,
        evaluation_order: vakint::EvaluationOrder(vec![vakint::EvaluationMethod::AlphaLoop(
            vakint::AlphaLoopOptions {
                susbstitute_masters: false,
            },
        )]),
        ..settings.vakint.true_settings()
    };
    let regulated_integrated = super::Integrated::new(crate::utils::vakint()?, &backend).run(
        &regulated_child,
        &regulated_ctx,
        child,
        &root,
        child,
        &root,
    )?;
    let regulated_prefix = Full4dCts::recursion_input(
        &regulated_child,
        &regulated_integrated,
        ApproximationType::OS,
        false,
        child.lmb(),
    )?;
    let original_finite = record
        .completed
        .sectors()
        .filter(|sector| !sector.frozen_lmbs().is_empty())
        .collect::<Vec<_>>();
    let regulated_finite = regulated_prefix
        .sectors()
        .filter(|sector| !sector.frozen_lmbs().is_empty())
        .collect::<Vec<_>>();
    ensure!(original_finite.len() == 1 && regulated_finite.len() == 1);
    ensure!(
        original_finite[0]
            .atom
            .expression()
            .contains_symbol(crate::utils::GS.dim_epsilon)
    );
    ensure!(
        original_finite[0]
            .atom
            .same_value(&regulated_finite[0].atom)
    );
    ensure!(original_finite[0].frozen_lmbs() == regulated_finite[0].frozen_lmbs());
    let mut inputs = 0;
    for (given, original_prefix, new_prefix) in [
        (&root, &empty, &empty),
        (child, &record.completed, &regulated_prefix),
    ] {
        let parent = OwnedForestNode {
            spinney: parent_spinney.clone(),
            topo_order: given.topo_order + 1,
        };
        let original = uv_limit(original_prefix, &ctx, &parent, given, &parent, given)?;
        let regulator = uv_limit(new_prefix, &regulated_ctx, &parent, given, &parent, given)?;
        ensure!(original.active_sectors().len() == regulator.active_sectors().len());
        ensure!(original.active_sectors().len() == if given.topo_order == 0 { 1 } else { 2 });
        for (original, regulator) in original
            .active_sectors()
            .iter()
            .zip(regulator.active_sectors())
        {
            ensure!(original.active_components == regulator.active_components);
            ensure!(original.frozen_lmbs() == regulator.frozen_lmbs());
            exact_vacuum_equal(original.integration_atom(), regulator.integration_atom())?;
            ensure!(
                !regulator
                    .integration_atom()
                    .expression()
                    .contains_symbol(lambda_symbol)
            );
            inputs += 1;
        }
    }
    ensure!(inputs == 3);
    crate::debug_tags!(#uv, #integrated, #summary;
        stage = "native_d09_ir_mass_vacuum_invariance",
        original_child_vacuum_equal = true, child_positive_epsilon_coefficient_equal = true,
        parent_original_vacuum_inputs = inputs, typed_gluon_mass_overrides = 2,
        physical_local_changed = true, native_gauge_numerator_unchanged = true,
        two_loop_backend_repeated = false, unregulated_exchange_proved = false,
        "A physical positive gluon-mass rearrangement preserves the actual three signed vacuum inputs"
    );
    Ok(())
}
