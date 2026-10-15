use ahash::HashMap;
use idenso::dirac::AGS;
use idenso::representations::Bispinor;
use linnet::half_edge::involution::EdgeIndex;
use spenso::{
    algebra::complex::Complex,
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::panicing::ErroringLibrary,
        library::symbolic::{ETS, ExplicitKey, TensorLibrary},
        parsing::ParseSettings,
    },
    structure::{
        Canonicalized, TensorDataLayout,
        representation::{LibraryRep, Minkowski, RepName},
    },
    tensors::{
        data::{DataTensor, DenseTensor, SparseOrDense, StorageTensor},
        parametric::{ParamTensor, atomcore::TensorAtomMaps},
    },
    vector, vector_symbol,
};
use spenso_hep_lib::hep_lib_atom;
use std::sync::OnceLock;
use symbolica::id::ReplaceWith;
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};
use tracing::debug;

use crate::{
    momentum::{FourMomentum, Polarization, Sign, SignOrZero},
    numerator::{ParsingNet, aind::Aind},
    settings::global::VectorPolarizationSumGauge,
    utils::{ApproxEq, F, GS, W_},
};

use super::{Model, Particle, ParticleGammaLoopExt, UFOSymbol};

type ComplexF64 = Complex<F<f64>>;
type SymComplexF64 = symbolica::domains::float::Complex<F<f64>>;
type SymComplexPlainF64 = symbolica::domains::float::Complex<f64>;
type TestTensorLibrary = TensorLibrary<ParamTensor<ExplicitKey<Aind>>, Aind>;

static TEST_INITIALIZED: OnceLock<()> = OnceLock::new();

fn ensure_test_initialized() {
    TEST_INITIALIZED.get_or_init(|| {
        crate::initialisation::test_initialise().unwrap();
    });
}

fn model() -> &'static Model {
    static MODEL: OnceLock<Model> = OnceLock::new();
    MODEL.get_or_init(Model::standard_model)
}

fn insert_explicit_complex_tensor(
    lib: &mut TestTensorLibrary,
    key: Canonicalized<ExplicitKey<Aind>>,
    data: Vec<ComplexF64>,
) {
    let layout = TensorDataLayout::from_canonicalized(&key)
        .expect("test tensor dimensions must be concrete");
    let storage_data = layout
        .reorder_to_storage(data)
        .expect("test tensor data must match its logical shape");
    let tensor: Canonicalized<ParamTensor<_>> = key.map_canonical(|structure| {
        DenseTensor::from_storage_data(storage_data, structure)
            .unwrap()
            .map_data(|c| Atom::num(c.re) + Atom::num(c.im) * Atom::i())
            .into()
    });
    lib.insert_explicit(tensor);
}

fn insert_rank_one_tensor(
    lib: &mut TestTensorLibrary,
    name: symbolica::atom::Symbol,
    rep: LibraryRep,
    data: Vec<ComplexF64>,
) {
    insert_explicit_complex_tensor(
        lib,
        ExplicitKey::from_iter([rep.new_rep(4)], name, None),
        data,
    );
}

fn insert_identity_tensor(
    lib: &mut TestTensorLibrary,
    name: symbolica::atom::Symbol,
    rep: LibraryRep,
) {
    let key = ExplicitKey::from_iter([rep.new_rep(4), rep.new_rep(4)], name, None);
    lib.insert_explicit(key.map_canonical(TestTensorLibrary::identity));
}

fn insert_spinor_tensor(
    lib: &mut TestTensorLibrary,
    name: symbolica::atom::Symbol,
    polarization: Polarization<ComplexF64>,
) {
    insert_rank_one_tensor(
        lib,
        name,
        LibraryRep::from(Bispinor {}),
        polarization.tensor.data,
    );
}

fn insert_vector_tensor(
    lib: &mut TestTensorLibrary,
    name: symbolica::atom::Symbol,
    polarization: Polarization<ComplexF64>,
) {
    insert_rank_one_tensor(
        lib,
        name,
        LibraryRep::from(Minkowski {}),
        polarization.tensor.data,
    );
}

fn fermion_sum_test_library(
    momentum: &FourMomentum<F<f64>>,
    eid: EdgeIndex,
    antiparticle: bool,
) -> TestTensorLibrary {
    let mut lib: TestTensorLibrary = hep_lib_atom();
    let positive = if antiparticle {
        momentum.v(Sign::Positive)
    } else {
        momentum.u(Sign::Positive)
    };
    let negative = if antiparticle {
        momentum.v(Sign::Negative)
    } else {
        momentum.u(Sign::Negative)
    };
    // `GS.emr_mom(eid, mu)` contracts with the HEP-library gamma tensor, so we must register the
    // contravariant momentum components p^mu. Lowering them here would apply the Minkowski signs
    // twice once the Lorentz contraction is executed by Spenso.
    insert_explicit_complex_tensor(
        &mut lib,
        ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.emr_mom,
            Some(vec![Atom::num(eid.0)]),
        ),
        vec![
            Complex::new_re(momentum.temporal.value),
            Complex::new_re(momentum.spatial.px),
            Complex::new_re(momentum.spatial.py),
            Complex::new_re(momentum.spatial.pz),
        ],
    );
    insert_spinor_tensor(&mut lib, vector_symbol!(spin_plus), positive.clone());
    insert_spinor_tensor(&mut lib, vector_symbol!(spin_plus_bar), positive.bar());
    insert_spinor_tensor(&mut lib, vector_symbol!(spin_minus), negative.clone());
    insert_spinor_tensor(&mut lib, vector_symbol!(spin_minus_bar), negative.bar());
    insert_identity_tensor(&mut lib, ETS.metric, LibraryRep::from(Bispinor {}));

    lib
}

fn vector_sum_test_library(
    momentum: &FourMomentum<F<f64>>,
    eid: EdgeIndex,
    include_longitudinal: bool,
) -> TestTensorLibrary {
    let mut lib: TestTensorLibrary = hep_lib_atom();

    let plus = momentum.eps_pol(SignOrZero::Plus);
    let minus = momentum.eps_pol(SignOrZero::Minus);

    insert_vector_tensor(&mut lib, vector_symbol!(eps_plus), plus.clone());
    insert_vector_tensor(&mut lib, vector_symbol!(eps_plus_bar), plus.bar());
    insert_vector_tensor(&mut lib, vector_symbol!(eps_minus), minus.clone());
    insert_vector_tensor(&mut lib, vector_symbol!(eps_minus_bar), minus.bar());

    if include_longitudinal {
        let zero = momentum.eps_pol(SignOrZero::Zero);
        insert_vector_tensor(&mut lib, vector_symbol!(eps_zero), zero.clone());
        insert_vector_tensor(&mut lib, vector_symbol!(eps_zero_bar), zero.bar());
    }

    // `GS.emr_mom(eid, mu)` is the full contravariant external momentum. We register p^mu directly
    // so the HEP-library metric handles the lowering when the tensor network contracts indices.
    insert_explicit_complex_tensor(
        &mut lib,
        ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.emr_mom,
            Some(vec![Atom::num(eid.0)]),
        ),
        vec![
            Complex::new_re(momentum.temporal.value),
            Complex::new_re(momentum.spatial.px),
            Complex::new_re(momentum.spatial.py),
            Complex::new_re(momentum.spatial.pz),
        ],
    );
    // `GS.emr_vec(eid, mink(4))` is the spatial part of the momentum embedded as a Lorentz
    // vector with vanishing time component. The axial-gauge reference vector and denominator use
    // this same Minkowski representation; its self-contraction is therefore `-|q⃗|²`.
    insert_explicit_complex_tensor(
        &mut lib,
        ExplicitKey::from_iter(
            [Minkowski {}.new_rep(4)],
            GS.emr_vec,
            Some(vec![Atom::num(eid.0)]),
        ),
        vec![
            Complex::new_re(F(0.0)),
            Complex::new_re(momentum.spatial.px),
            Complex::new_re(momentum.spatial.py),
            Complex::new_re(momentum.spatial.pz),
        ],
    );
    lib
}

fn bispinor_slot(index: i64) -> Atom {
    Bispinor {}.new_rep(4).to_symbolic([Atom::num(index)])
}

fn minkowski_slot(index: i64) -> Atom {
    Minkowski {}.new_rep(4).to_symbolic([Atom::num(index)])
}

fn fermion_sum_lhs_expr() -> Atom {
    let row = bispinor_slot(1);
    let col = bispinor_slot(2);

    vector!(spin_plus, row.clone()) * vector!(spin_plus_bar, col.clone())
        + vector!(spin_minus, row) * vector!(spin_minus_bar, col)
}

fn vector_sum_lhs_expr(include_longitudinal: bool) -> Atom {
    let row = minkowski_slot(1);
    let col = minkowski_slot(2);
    let transverse = vector!(eps_plus, row.clone()) * vector!(eps_plus_bar, col.clone())
        + vector!(eps_minus, row.clone()) * vector!(eps_minus_bar, col.clone());

    if include_longitudinal {
        transverse + vector!(eps_zero, row) * vector!(eps_zero_bar, col)
    } else {
        transverse
    }
}

fn vector_sum_rhs_expr(
    particle: &Particle,
    eid: EdgeIndex,
    gauge: VectorPolarizationSumGauge,
    average: bool,
) -> Atom {
    let row = minkowski_slot(1);
    let col = minkowski_slot(2);
    (function!(GS.epsilon, eid.0, row) * function!(GS.epsilonbar, eid.0, col)).replace_multiple(&[
        particle
            .polarization_sum(model(), eid, average, gauge)
            .unwrap()
            .unwrap(),
    ])
}

fn evaluate_tensor_network(
    expression: Atom,
    lib: &TestTensorLibrary,
) -> DenseTensor<ComplexF64, spenso::network::parsing::ShadowedStructure<Aind>> {
    let empty_constants: HashMap<Atom, SymComplexF64> = HashMap::default();
    evaluate_tensor_network_with_constants(expression, lib, &empty_constants)
}

fn evaluate_tensor_network_with_constants(
    expression: Atom,
    lib: &TestTensorLibrary,
    constants: &HashMap<Atom, SymComplexF64>,
) -> DenseTensor<ComplexF64, spenso::network::parsing::ShadowedStructure<Aind>> {
    let mut network =
        ParsingNet::try_from_view(expression.as_view(), lib, &ParseSettings::default()).unwrap();

    network
        .execute::<Sequential, SmallestDegree, _, _, _>(lib, &ErroringLibrary::new())
        .unwrap();

    let ExecutionResult::Val(result) = network.result_tensor(lib).unwrap() else {
        panic!("polarization-sum test network did not evaluate to a tensor");
    };

    let constants_f64: HashMap<Atom, SymComplexPlainF64> = constants
        .iter()
        .map(|(atom, value)| {
            (
                atom.clone(),
                SymComplexPlainF64::new(value.re.0, value.im.0),
            )
        })
        .collect();
    let dense = result
        .into_owned()
        .evaluate(&constants_f64)
        .unwrap()
        .map_data(|c| Complex::new(F(c.re), F(c.im)))
        .to_dense();

    let DataTensor::Dense(dense) = dense else {
        panic!("polarization-sum test tensor did not densify into a dense tensor");
    };

    dense
}

fn assert_tensor_close<S: PartialEq + std::fmt::Debug>(
    actual: &DenseTensor<ComplexF64, S>,
    expected: &DenseTensor<ComplexF64, S>,
) {
    assert_eq!(actual.data.len(), expected.data.len());
    let tolerance = F(1.0e-12);
    for (entry, (actual, expected)) in actual.data.iter().zip(expected.data.iter()).enumerate() {
        actual.approx_eq_res(expected, &tolerance).unwrap_or_else(|err| {
            panic!(
                "tensor entry {entry} mismatch: actual={actual}, expected={expected}, error={err}"
            )
        });
    }
}

fn fermion_sum_rhs_expr(particle: &Particle, eid: EdgeIndex) -> Atom {
    let row_slot = bispinor_slot(1);
    let col_slot = bispinor_slot(2);
    if particle.is_antiparticle() {
        function!(GS.v, eid.0, row_slot.clone()) * function!(GS.vbar, eid.0, col_slot.clone())
    } else {
        function!(GS.u, eid.0, row_slot.clone()) * function!(GS.ubar, eid.0, col_slot.clone())
    }
    .replace_multiple(&[particle
        .polarization_sum(
            model(),
            eid,
            false,
            VectorPolarizationSumGauge::LightLikeAxial,
        )
        .unwrap()
        .unwrap()])
}

fn assert_fermion_sum_matches_rule(momentum: FourMomentum<F<f64>>, mass: i64) {
    ensure_test_initialized();
    let mut constants = HashMap::default();
    constants.insert(
        Atom::from(UFOSymbol::from("MT")),
        SymComplexF64::new(F(mass as f64), F(0.0)),
    );

    for antiparticle in [false, true] {
        let particle = if antiparticle {
            dummy_spinor_particle(-11)
        } else {
            dummy_spinor_particle(11)
        };
        let lib = fermion_sum_test_library(&momentum, EdgeIndex(7), antiparticle);
        let lhs = evaluate_tensor_network(fermion_sum_lhs_expr(), &lib);
        let rhs = evaluate_tensor_network_with_constants(
            fermion_sum_rhs_expr(&particle, EdgeIndex(7)),
            &lib,
            &constants,
        );
        debug!(
            antiparticle,
            mass,
            momentum = %momentum,
            "fermion polarization sum tensors\nlhs:\n{lhs}\n\nrhs:\n{rhs}"
        );
        assert_tensor_close(&lhs, &rhs);
    }
}

fn replacement_rhs_atom(replacement: symbolica::id::Replacement) -> Atom {
    let ReplaceWith::Pattern(rhs) = replacement.rhs else {
        panic!("polarization sum rhs should be a pattern");
    };
    rhs.to_atom().unwrap()
}

fn simplify_fermion_rhs_for_symbolic_checks(rhs: Atom, eid: EdgeIndex) -> Atom {
    rhs.replace(
        function!(GS.emr_mom, eid.0 as i64, W_.m_) * function!(AGS.gamma, W_.a_, W_.b_, W_.m_),
    )
    .with(symbol!("PSLASH"))
    .replace(Atom::from(UFOSymbol::from("MT")) * function!(ETS.metric, W_.a_, W_.b_))
    .with(symbol!("MASS_TERM"))
}

fn assert_vector_sum_matches_rule(
    momentum: FourMomentum<F<f64>>,
    mass_parameter: UFOSymbol,
    gauge: VectorPolarizationSumGauge,
    include_longitudinal: bool,
) {
    ensure_test_initialized();
    let particle = dummy_vector_particle(22, mass_parameter);
    let mut constants = HashMap::default();
    constants.insert(
        Atom::from(UFOSymbol::from("MT")),
        SymComplexF64::new(momentum.norm(), F(0.0)),
    );

    let eid = EdgeIndex(7);
    // The axial reference vector uses the scalar Q(eid, cind(0)) as well as
    // the full momentum tensor. Give both representations the same physical
    // time component; registering the tensor does not define this scalar.
    constants.insert(
        GS.emr_mom(eid, GS.cind(0)),
        SymComplexF64::new(momentum.temporal.value, F(0.0)),
    );
    let lib = vector_sum_test_library(&momentum, eid, include_longitudinal);
    let lhs = evaluate_tensor_network(vector_sum_lhs_expr(include_longitudinal), &lib);
    let rhs = evaluate_tensor_network_with_constants(
        vector_sum_rhs_expr(&particle, eid, gauge, false),
        &lib,
        &constants,
    );
    debug!(
        ?gauge,
        include_longitudinal,
        momentum = %momentum,
        "vector polarization sum tensors\nlhs:\n{lhs}\n\nrhs:\n{rhs}"
    );
    assert_tensor_close(&lhs, &rhs);
}

fn dummy_spinor_particle(pdg_code: i64) -> Particle {
    test_particle(pdg_code, 2, false, "psi")
}

fn dummy_vector_particle(pdg_code: i64, mass: UFOSymbol) -> Particle {
    test_particle(pdg_code, 3, mass.is_zero(), "vec")
}

fn test_particle(pdg_code: i64, spin: i64, massless: bool, name: &str) -> Particle {
    let mut particle = model()
        .particle_by_pdg(if pdg_code < 0 { -11 } else { 11 })
        .unwrap()
        .clone();
    particle.pdg_code = pdg_code;
    particle.name = name.into();
    particle.spin = spin;
    particle.mass = model()
        .parameter_id(if massless { "ZERO" } else { "MT" })
        .unwrap();
    particle
}

#[test]
fn fermion_polarization_sum_symbolic_mass_term_sign_matches_completeness_relation() {
    ensure_test_initialized();
    let eid = EdgeIndex(7);
    let particle = dummy_spinor_particle(11);
    let antiparticle = dummy_spinor_particle(-11);

    let particle_rhs = replacement_rhs_atom(
        particle
            .polarization_sum(
                model(),
                eid,
                false,
                VectorPolarizationSumGauge::LightLikeAxial,
            )
            .unwrap()
            .unwrap(),
    );
    let antiparticle_rhs = replacement_rhs_atom(
        antiparticle
            .polarization_sum(
                model(),
                eid,
                false,
                VectorPolarizationSumGauge::LightLikeAxial,
            )
            .unwrap()
            .unwrap(),
    );

    let particle_simplified = particle_rhs
        .replace(function!(GS.emr_mom, eid.0 as i64, W_.m_))
        .with(Atom::num(0))
        .replace(Atom::from(UFOSymbol::from("MT")))
        .with(symbol!("MASS"))
        .replace(function!(ETS.metric, W_.a_, W_.b_))
        .with(symbol!("ID"));
    let antiparticle_simplified = antiparticle_rhs
        .replace(function!(GS.emr_mom, eid.0 as i64, W_.m_))
        .with(Atom::num(0))
        .replace(Atom::from(UFOSymbol::from("MT")))
        .with(symbol!("MASS"))
        .replace(function!(ETS.metric, W_.a_, W_.b_))
        .with(symbol!("ID"));

    let expected_particle = Atom::var(symbol!("MASS")) * Atom::var(symbol!("ID"));
    let expected_antiparticle = -Atom::var(symbol!("MASS")) * Atom::var(symbol!("ID"));

    assert_eq!(
        particle_simplified.to_canonical_string(),
        expected_particle.to_canonical_string()
    );
    assert_eq!(
        antiparticle_simplified.to_canonical_string(),
        expected_antiparticle.to_canonical_string()
    );
}

#[test]
fn fermion_polarization_sum_average_includes_half_factor() {
    ensure_test_initialized();
    let eid = EdgeIndex(7);
    let particle = dummy_spinor_particle(11);

    let summed = simplify_fermion_rhs_for_symbolic_checks(
        replacement_rhs_atom(
            particle
                .polarization_sum(
                    model(),
                    eid,
                    false,
                    VectorPolarizationSumGauge::LightLikeAxial,
                )
                .unwrap()
                .unwrap(),
        ),
        eid,
    );
    let averaged = simplify_fermion_rhs_for_symbolic_checks(
        replacement_rhs_atom(
            particle
                .polarization_sum(
                    model(),
                    eid,
                    true,
                    VectorPolarizationSumGauge::LightLikeAxial,
                )
                .unwrap()
                .unwrap(),
        ),
        eid,
    );

    assert_eq!(
        averaged.to_canonical_string(),
        (summed / Atom::num(2)).to_canonical_string()
    );
}

#[test]
fn vector_polarization_sum_average_uses_two_for_massless_and_three_for_massive() {
    ensure_test_initialized();
    let eid = EdgeIndex(7);
    let massless = dummy_vector_particle(22, UFOSymbol::zero());
    let massive = dummy_vector_particle(23, UFOSymbol::from("MT"));

    let massless_rhs = replacement_rhs_atom(
        massless
            .polarization_sum(model(), eid, true, VectorPolarizationSumGauge::Feynman)
            .unwrap()
            .unwrap(),
    );
    let massive_rhs = replacement_rhs_atom(
        massive
            .polarization_sum(model(), eid, true, VectorPolarizationSumGauge::Feynman)
            .unwrap()
            .unwrap(),
    );

    assert_eq!(
        massless_rhs.to_canonical_string(),
        (-function!(ETS.metric, W_.a_, W_.b_) / Atom::num(2)).to_canonical_string()
    );
    assert_eq!(
        massive_rhs.to_canonical_string(),
        (-function!(ETS.metric, W_.a_, W_.b_) / Atom::num(3)).to_canonical_string()
    );
}

#[test]
fn vector_polarization_sum_feynman_gauge_is_always_minus_metric() {
    ensure_test_initialized();
    let eid = EdgeIndex(7);
    let massless = dummy_vector_particle(22, UFOSymbol::zero());
    let massive = dummy_vector_particle(23, UFOSymbol::from("MT"));

    let massless_rhs = replacement_rhs_atom(
        massless
            .polarization_sum(model(), eid, false, VectorPolarizationSumGauge::Feynman)
            .unwrap()
            .unwrap(),
    );
    let massive_rhs = replacement_rhs_atom(
        massive
            .polarization_sum(model(), eid, false, VectorPolarizationSumGauge::Feynman)
            .unwrap()
            .unwrap(),
    );

    assert_eq!(
        massless_rhs.to_canonical_string(),
        (-function!(ETS.metric, W_.a_, W_.b_)).to_canonical_string()
    );
    assert_eq!(
        massive_rhs.to_canonical_string(),
        (-function!(ETS.metric, W_.a_, W_.b_)).to_canonical_string()
    );
}

#[test]
fn higher_spin_polarization_sums_report_structured_unsupported_errors() {
    ensure_test_initialized();
    for spin in [4, 5] {
        for massless in [false, true] {
            let particle = test_particle(9_000_001, spin, massless, "unsupported_spin");
            let error = particle
                .polarization_sum(
                    model(),
                    EdgeIndex(7),
                    true,
                    VectorPolarizationSumGauge::LightLikeAxial,
                )
                .unwrap_err();
            assert!(matches!(
                error.downcast_ref::<feynkit_amplitude::SpinSumError>(),
                Some(feynkit_amplitude::SpinSumError::UnsupportedSpin(actual)) if *actual == spin
            ));
        }
    }
}

#[test]
fn fermion_polarization_sum_massless_z_axis() {
    assert_fermion_sum_matches_rule(FourMomentum::from_args(F(4.0), F(0.0), F(0.0), F(4.0)), 0);
}

#[test]
fn fermion_polarization_sum_massless_antiz_axis() {
    assert_fermion_sum_matches_rule(FourMomentum::from_args(F(4.0), F(0.0), F(0.0), F(-4.0)), 0);
}

#[test]
fn fermion_polarization_sum_massive_z_axis() {
    assert_fermion_sum_matches_rule(FourMomentum::from_args(F(5.0), F(0.0), F(0.0), F(4.0)), 3);
}

#[test]
fn fermion_polarization_sum_massless_transverse() {
    assert_fermion_sum_matches_rule(FourMomentum::from_args(F(4.0), F(0.0), F(4.0), F(0.0)), 0);
}

#[test]
fn fermion_polarization_sum_massive_transverse() {
    assert_fermion_sum_matches_rule(FourMomentum::from_args(F(5.0), F(3.0), F(0.0), F(0.0)), 4);
}

#[test]
fn vector_polarization_sum_massive_z_axis() {
    assert_vector_sum_matches_rule(
        FourMomentum::from_args(F(5.0), F(0.0), F(0.0), F(4.0)),
        UFOSymbol::from("MT"),
        VectorPolarizationSumGauge::LightLikeAxial,
        true,
    );
}

#[test]
fn vector_polarization_sum_massive_transverse() {
    assert_vector_sum_matches_rule(
        FourMomentum::from_args(F(5.0), F(3.0), F(0.0), F(0.0)),
        UFOSymbol::from("MT"),
        VectorPolarizationSumGauge::LightLikeAxial,
        true,
    );
}

#[test]
fn vector_polarization_sum_massless_z_axis() {
    assert_vector_sum_matches_rule(
        FourMomentum::from_args(F(4.0), F(0.0), F(0.0), F(4.0)),
        UFOSymbol::zero(),
        VectorPolarizationSumGauge::LightLikeAxial,
        false,
    );
}

#[test]
fn vector_polarization_sum_massless_antiz_axis() {
    assert_vector_sum_matches_rule(
        FourMomentum::from_args(F(4.0), F(0.0), F(0.0), F(-4.0)),
        UFOSymbol::zero(),
        VectorPolarizationSumGauge::LightLikeAxial,
        false,
    );
}

#[test]
fn vector_polarization_sum_massless_transverse() {
    assert_vector_sum_matches_rule(
        FourMomentum::from_args(F(4.0), F(0.0), F(4.0), F(0.0)),
        UFOSymbol::zero(),
        VectorPolarizationSumGauge::LightLikeAxial,
        false,
    );
}
