use ahash::{HashMap, HashMapExt};
use idenso::{Cookable, representations::Bispinor};
use spenso::{
    algebra::upgrading_arithmetic::FallibleSub,
    iterators::IteratableTensor,
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        library::symbolic::{ETS, ExplicitKey},
        parsing::{ParseSettings, ShadowedStructure, StrictTensorFilter},
    },
    shadowing::Concretize,
    structure::{
        TensorStructure,
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
    tensors::{
        data::{DenseTensor, SparseOrDense},
        parametric::{MixedTensor, ParamTensor},
    },
};
use spenso_hep_lib::{FUN_LIB, HepNet, hep_lib_atom};

use symbolica::{
    atom::{Atom, AtomCore},
    function, parse_lit, symbol,
};

use crate::common::{gamma, gamma0, gammaadj, gammaconj, p, q, test_initialize, u, ub};

mod common;
#[test]
fn validate() {
    let _ = (
        spenso::vector_symbol!("spenso::p"),
        spenso::vector_symbol!("spenso::q"),
    );
    test_initialize();
    let mut const_map = HashMap::new();
    let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::p"),
        None,
    )
    .into_canonical()
    .to_shell()
    .concretize()
    .unwrap();

    for (i, a) in pt.iter_flat() {
        const_map.insert(a.clone(), Atom::num(usize::from(i)));
    }

    let qt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::q"),
        None,
    )
    .into_canonical()
    .to_shell()
    .concretize()
    .unwrap();

    for (i, a) in qt.iter_flat() {
        const_map.insert(a.clone(), Atom::num(usize::from(i) + 1));
    }

    // gamma.reindex_storage(&[1, 2, 3]).unwrap().apply()

    let expr =
        p(1) * (p(3) + q(3)) * gamma(1, 2, 1) * gamma(2, 3, 2) * gamma(3, 4, 3) * gamma(4, 1, 4);

    validate_gamma(expr, const_map.clone());
    let _expr = p(1) * p(1);

    let bis = Bispinor {}.new_rep(2);
    let mink = Minkowski {}.new_rep(2);
    let _expr = function!(
        symbol!("A"),
        bis.slot::<AbstractIndex, _>(1).to_atom(),
        bis.slot::<AbstractIndex, _>(2).to_atom(),
        mink.slot::<AbstractIndex, _>(1).to_atom()
    ) * function!(
        symbol!("A"),
        bis.slot::<AbstractIndex, _>(2).to_atom(),
        bis.slot::<AbstractIndex, _>(1).to_atom(),
        mink.slot::<AbstractIndex, _>(2).to_atom()
    );

    let expr = gamma(1, 2, 1) * gamma(2, 1, 1);
    validate_gamma(expr, const_map.clone());
    let expr = gamma(1, 2, 2) * gamma(2, 1, 1) + gamma(1, 2, 1) * gamma(2, 1, 2);
    validate_gamma(expr, const_map.clone());

    let expr = gamma0(1, 2) * gamma(2, 3, 1) * gamma0(3, 4) - gammaconj(4, 1, 1);
    // let expr2 = gamma0(1, 2) * gammaconj(3, 2, 1) * gamma0(3, 4);

    validate_gamma(expr, const_map.clone());

    let expr = gammaadj(1, 2, 1) - gamma0(1, 3) * gamma(3, 4, 1) * gamma0(4, 2);
    // let expr2 = gamma0(1, 2) * gammaconj(3, 2, 1) * gamma0(3, 4);

    validate_gamma(expr, const_map.clone());

    let expr = gammaadj(1, 2, 1) - gammaconj(2, 1, 1);
    // let expr2 = gamma0(1, 2) * gammaconj(3, 2, 1) * gamma0(3, 4);

    validate_gamma(expr, const_map.clone());
    let _a = u(1, 1) * gamma(1, 2, 1) * ub(2, 2);
    // validate_gamma(
    //     a.conj()
    //         .replace(function!(Symbol::CONJ, symbol!("a__")))
    //         .with(function!(
    //             symbol!("conju", tag = SPENSO_TAG.tag),
    //             symbol!("a__")
    //         ))
    //         * a,
    //     const_map.clone(),
    // );

    // let a = Atom::num(1) / u(1, 1) * gamma(1, 2, 1) * ub(2, 2);
    // validate_gamma(
    //     a.conj()
    //         .replace(function!(Symbol::CONJ, symbol!("a__")))
    //         .with(function!(
    //             symbol!("conju", tag = SPENSO_TAG.tag),
    //             symbol!("a__")
    //         ))
    //         * a,
    //     const_map.clone(),
    // );

    // validate_gamma(expr2, const_map.clone());
    // let expr = gamma(1, 2, 2) * gamma(2, 1, 1) + gamma(1, 2, 1) * gamma(2, 1, 2);
    // // + gamma(1, 2, 1) * gamma(2, 1, 1);

    // // let expr = A(1, 2, 0) * B(2, 1, 3);
    // validate_gamma(expr, const_map.clone());
    // let expr = gamma(1, 2, 1);

    // validate_gamma(expr, const_map.clone());
    // let expr = gamma(2, 1, 1);

    // validate_gamma(expr, const_map.clone());
    // assert_eq!(pt, qt);

    // let a: DataTensor<_, _> = DenseTensor::fill(
    //     OrderedStructure::<Bispinor>::from_iter([Bispinor {}.new_slot(2, 2)]).into_canonical(),
    //     Complex::new(-1., 0.),
    // )
    // .into();
    // let b: DataTensor<_, _> = DenseTensor::fill(
    //     OrderedStructure::<Bispinor>::from_iter([Bispinor {}.new_slot(2, 2)]).into_canonical(),
    //     Complex::new(1., 0.),
    // )
    // .into();
    // assert_eq!(
    //     ParamOrConcrete::<_, OrderedStructure<Bispinor>>::Concrete(a),
    //     ParamOrConcrete::<_, OrderedStructure<Bispinor>>::Concrete(b)
    // );
}

fn validate_gamma(expr: Atom, const_map: HashMap<Atom, Atom>) {
    let mut library = hep_lib_atom::<AbstractIndex, MixedTensor<f64, ExplicitKey<AbstractIndex>>>();
    // Keep metric leaves exact too: the mixed library's generic factory uses f64.
    for (representation, lorentzian) in [
        (Minkowski {}.new_rep(4).to_lib(), true),
        (Bispinor {}.new_rep(4).to_lib(), false),
    ] {
        let key = ExplicitKey::from_iter([representation; 2], ETS.metric, None);
        let components = (0..16)
            .map(|flat| {
                Atom::num(if flat / 4 != flat % 4 {
                    0
                } else if lorentzian && flat > 0 {
                    -1
                } else {
                    1
                })
            })
            .collect();
        library.insert_explicit(key.map_canonical(|structure| {
            MixedTensor::Param(ParamTensor::param(
                DenseTensor::from_storage_data(components, structure)
                    .unwrap()
                    .into(),
            ))
        }));
    }
    let settings =
        ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
    let mut net =
        HepNet::<AbstractIndex>::try_from_view(expr.as_view(), &library, &settings).unwrap();

    println!("Expression: {}", expr);
    let simplified = idenso::tensor::SymbolicTensor::infer(expr)
        .unwrap()
        .simplify_gamma(idenso::dirac::GammaSimplifySettings::default())
        .unwrap();

    println!("Simplified to {}", simplified.root().expression());
    // Keep the declared ports when simplification produces a typed zero.
    let (mut net_simplified, scalar_definitions): (spenso_hep_lib::HepNet<AbstractIndex>, _) =
        simplified
            .to_network(
                &library,
                &*FUN_LIB,
                &settings,
                None,
                |net| net.execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB),
                Ok,
                |_| false,
            )
            .unwrap();
    assert!(scalar_definitions.get_aliases().is_empty());

    println!("{}", net_simplified.dot_pretty());
    net.execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)
        .unwrap();
    net_simplified
        .execute::<Sequential, SmallestDegree, _, _, _>(&library, &*FUN_LIB)
        .unwrap();

    if let ExecutionResult::Val(v) = net.result_tensor(&library).unwrap() {
        if let ExecutionResult::Val(v2) = net_simplified.result_tensor(&library).unwrap() {
            let mut res = v.into_owned();
            println!("{res}");
            let mut res_simplified = v2.into_owned();
            // Apply the same integer assignments exactly before comparing components.
            // Rounding each route separately can turn exact cancellations into nonzero floats.
            for (component, value) in &const_map {
                res = res.replace(component.to_pattern()).with(value.to_pattern());
                res_simplified = res_simplified
                    .replace(component.to_pattern())
                    .with(value.to_pattern());
            }
            res = res.to_dense();
            res_simplified = res_simplified.to_dense();

            let mut sub = res.sub_fallible(&res_simplified).unwrap();
            sub.to_param();
            let sub = sub.try_into_parametric().unwrap();
            for (index, component) in sub.iter_flat() {
                assert!(
                    component.is_zero(),
                    "component {index}: {component}\noriginal: {res}\nsimplified: {res_simplified}"
                );
            }
        } else {
            panic!("Expected tensor result");
        }
    } else {
        panic!("Expected tensor result");
    }
}

mod failing {
    use super::*;

    #[test]
    fn gl_03() {
        let _ = (
            spenso::vector_symbol!("spenso::P"),
            spenso::vector_symbol!("spenso::K"),
        );
        test_initialize();
        let mut const_map = HashMap::new();
        let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
            [Minkowski {}.new_slot(4, 1)],
            symbol!("spenso::P"),
            Some(vec![Atom::num(0)]),
        )
        .into_canonical()
        .to_shell()
        .concretize()
        .unwrap();

        for (i, a) in pt.iter_flat() {
            const_map.insert(a.clone(), Atom::num(usize::from(i)));
        }

        let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
            [Minkowski {}.new_slot(4, 1)],
            symbol!("spenso::K"),
            Some(vec![Atom::num(1)]),
        )
        .into_canonical()
        .to_shell()
        .concretize()
        .unwrap();

        for (i, a) in pt.iter_flat() {
            const_map.insert(a.clone(), Atom::num(usize::from(i)));
        }

        let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
            [Minkowski {}.new_slot(4, 1)],
            symbol!("spenso::K"),
            Some(vec![Atom::num(0)]),
        )
        .into_canonical()
        .to_shell()
        .concretize()
        .unwrap();

        for (i, a) in pt.iter_flat() {
            const_map.insert(a.clone(), Atom::num(usize::from(i)));
        }

        const_map.insert(parse_lit!(spenso::MC), Atom::num(11232));

        const_map.insert(parse_lit!(spenso::MW), Atom::num(1231));

        let expr = parse_lit!(
            1 / 6
                ^ 4
                ^ -2 * (MC * g(bis(4, hedge(1)), bis(4, hedge(2)))
                    - K(0, mink(4, edge(1, 1)))
                        * spenso::gamma(bis(4, hedge(1)), bis(4, hedge(2)), mink(4, edge(1, 1))))
                    * (-K(0, mink(4, edge(3, 1))) - K(1, mink(4, edge(3, 1))))
                    * (-g(mink(4, hedge(7)), mink(4, hedge(8))) + MW
                        ^ -2 * (-P(0, mink(4, hedge(7))) - K(1, mink(4, hedge(7))))
                            * (-P(0, mink(4, hedge(8))) - K(1, mink(4, hedge(8)))))
                    * (P(0, mink(4, edge(5, 1)))
                        + K(0, mink(4, edge(5, 1)))
                        + K(1, mink(4, edge(5, 1))))
                    * g(mink(4, hedge(0)), mink(4, hedge(8)))
                    * spenso::gamma(bis(4, hedge(10)), bis(4, hedge(6)), mink(4, hedge(11)))
                    * spenso::gamma(bis(4, hedge(2)), bis(4, vertex(1, 1)), mink(4, hedge(7)))
                    * spenso::gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge(3, 1)))
                    * spenso::gamma(bis(4, hedge(9)), bis(4, hedge(10)), mink(4, edge(5, 1)))
                    * projm(bis(4, hedge(5)), bis(4, hedge(1)))
                    * projm(bis(4, vertex(1, 1)), bis(4, hedge(9)))
                    * (1 / 2)
                ^ (1 / 2),
            default_namespace = "spenso"
        );

        validate_gamma(expr.cook_indices(), const_map.clone());
    }
}
