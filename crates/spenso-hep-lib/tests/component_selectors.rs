//! Exact coordinate extraction for the native rank-one component selector.

use idenso::{
    representations::{Bispinor, ColorFundamental, initialize},
    tensor::SymbolicTensor,
};
use spenso::{
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey},
        parsing::{ParseSettings, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        concrete_index::DualConciousIndex,
        representation::{Minkowski, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    tensors::{
        data::DenseTensor,
        parametric::{MixedTensor, ParamTensor},
    },
};
use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

#[test]
fn delta_selects_native_shadow_coordinates_in_both_variances() {
    initialize();
    for selected in 0..3 {
        for actual in 0..3 {
            for coordinate in [
                DualConciousIndex::SelfDual(actual),
                DualConciousIndex::Up(actual),
                DualConciousIndex::Down(actual),
            ] {
                let scalar = function!(
                    ETS.delta,
                    function!(AIND_SYMBOLS.cind, selected),
                    function!(AIND_SYMBOLS.cind, Atom::from(coordinate))
                );
                assert_eq!(scalar, Atom::num(i32::from(selected == actual)));
            }
        }
    }
}

#[test]
fn delta_matches_independent_unit_vectors_and_preserves_metric_pairing() {
    initialize();
    let fundamental = ColorFundamental {}.new_rep(3).to_lib();
    let representations = [
        Minkowski {}.new_rep(4).to_lib(),
        Bispinor {}.new_rep(4).to_lib(),
        fundamental,
        fundamental.dual(),
    ];
    let settings =
        ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
    let vector = SPENSO_TAG.rank_one_tensor_symbol("component_selector_reference::unit");
    let spectator = (Atom::var(symbol!("component_selector_reference::x"; Scalar))
        + Atom::var(symbol!("component_selector_reference::y"; Scalar)))
    .pow(9);
    for representation in representations {
        let dimension = usize::try_from(representation.dim).unwrap();
        let slot = representation.slot::<AbstractIndex, _>(AbstractIndex::Normal(821743));
        for populated in 0..dimension {
            // These independent library data do not use the delta callback.
            let mut library = hep_lib_atom();
            let key = ExplicitKey::<AbstractIndex>::from_iter([representation], vector, None);
            let entries = (0..dimension)
                .map(|component| Atom::num(i32::from(component == populated)))
                .collect();
            library.insert_explicit(key.map_canonical(|structure| {
                MixedTensor::Param(ParamTensor::param(
                    DenseTensor::from_storage_data(entries, structure)
                        .unwrap()
                        .into(),
                ))
            }));
            for selected in 0..dimension {
                let selector = function!(
                    ETS.delta,
                    function!(AIND_SYMBOLS.cind, selected),
                    slot.dual().to_atom()
                );
                assert_eq!(SymbolicTensor::infer(selector.clone()).unwrap().rank(), 1);
                let expression = &spectator * selector * function!(vector, slot.to_atom());
                assert!(
                    SymbolicTensor::infer(expression.clone())
                        .unwrap()
                        .is_scalar()
                );
                let mut network = HepNet::<AbstractIndex>::try_from_view(
                    expression.as_view(),
                    &library,
                    &settings,
                )
                .unwrap();
                network
                    .execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)
                    .unwrap();
                let actual = match network.result_scalar().unwrap() {
                    ExecutionResult::One => Atom::one(),
                    ExecutionResult::Zero => Atom::Zero,
                    ExecutionResult::Val(value) => value.into_owned(),
                };
                let sign = if representation.is_neg(selected) {
                    -1
                } else {
                    1
                };
                let raw_component = if selected == populated {
                    spectator.clone()
                } else {
                    Atom::Zero
                };
                assert_eq!(
                    actual,
                    Atom::num(sign) * &raw_component,
                    "wrong pairing for {representation}, selected={selected}, populated={populated}"
                );
                assert_eq!(
                    Atom::num(sign) * actual,
                    raw_component,
                    "metric-corrected selector must recover the raw coordinate"
                );
            }
        }
    }
}
