use ahash::{HashMap, HashMapExt};
use idenso::{IndexTooling, gamma::GammaSimplifier, representations::Bispinor};
use spenso::{
    algebra::upgrading_arithmetic::FallibleSub,
    iterators::IteratableTensor,
    network::{
        ExecutionResult, Sequential, SmallestDegree,
        parsing::{ParseSettings, ShadowedStructure},
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
        parametric::atomcore::TensorAtomMaps,
    },
};
use spenso_hep_lib::{FUN_LIB, HEP_LIB};

use symbolica::{atom::Atom, function, id::ConditionResult, parse_lit, symbol};

use crate::common::{HepAtomExt, gamma, gamma0, gammaadj, gammaconj, p, q, test_initialize, u, ub};

mod common;
#[test]
fn validate() {
    test_initialize();
    let mut const_map = HashMap::new();
    let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::p"),
        None,
    )
    .structure
    .to_shell()
    .concretize(None);

    for (i, a) in pt.iter_flat() {
        const_map.insert(
            a.clone(),
            symbolica::domains::float::Complex::new(usize::from(i) as f64 * 1., 0.),
        );
    }

    let qt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::q"),
        None,
    )
    .structure
    .to_shell()
    .concretize(None);

    for (i, a) in qt.iter_flat() {
        const_map.insert(
            a.clone(),
            symbolica::domains::float::Complex::new((usize::from(i) + 1) as f64 * 1., 0.),
        );
    }

    // gamma.reindex([1,2,3]).unwrap().map_structure(|a|)

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
    //     OrderedStructure::<Bispinor>::from_iter([Bispinor {}.new_slot(2, 2)]).structure,
    //     Complex::new(-1., 0.),
    // )
    // .into();
    // let b: DataTensor<_, _> = DenseTensor::fill(
    //     OrderedStructure::<Bispinor>::from_iter([Bispinor {}.new_slot(2, 2)]).structure,
    //     Complex::new(1., 0.),
    // )
    // .into();
    // assert_eq!(
    //     ParamOrConcrete::<_, OrderedStructure<Bispinor>>::Concrete(a),
    //     ParamOrConcrete::<_, OrderedStructure<Bispinor>>::Concrete(b)
    // );
}

#[test]
fn gl_03() {
    test_initialize();
    let mut const_map = HashMap::new();
    let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::P"),
        Some(vec![Atom::num(0)]),
    )
    .structure
    .to_shell()
    .concretize(None);

    for (i, a) in pt.iter_flat() {
        const_map.insert(
            a.clone(),
            symbolica::domains::float::Complex::new(usize::from(i) as f64 * 1., 0.),
        );
    }

    let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::K"),
        Some(vec![Atom::num(1)]),
    )
    .structure
    .to_shell()
    .concretize(None);

    for (i, a) in pt.iter_flat() {
        const_map.insert(
            a.clone(),
            symbolica::domains::float::Complex::new(usize::from(i) as f64 * 1., 0.),
        );
    }

    let pt: DenseTensor<Atom, _> = ShadowedStructure::<AbstractIndex>::from_iter(
        [Minkowski {}.new_slot(4, 1)],
        symbol!("spenso::K"),
        Some(vec![Atom::num(0)]),
    )
    .structure
    .to_shell()
    .concretize(None);

    for (i, a) in pt.iter_flat() {
        const_map.insert(
            a.clone(),
            symbolica::domains::float::Complex::new(usize::from(i) as f64 * 1., 0.),
        );
    }

    const_map.insert(
        parse_lit!(spenso::MC),
        symbolica::domains::float::Complex::new(11232., 0.),
    );

    const_map.insert(
        parse_lit!(spenso::MW),
        symbolica::domains::float::Complex::new(1231., 0.),
    );

    let expr = parse_lit!(
        1 / 6
            ^ 4
            ^ -2 * (MC * g(bis(4, hedge(1)), bis(4, hedge(2)))
                - K(0, mink(4, edge(1, 1)))
                    * gamma(bis(4, hedge(1)), bis(4, hedge(2)), mink(4, edge(1, 1))))
                * (-K(0, mink(4, edge(3, 1))) - K(1, mink(4, edge(3, 1))))
                * (-g(mink(4, hedge(7)), mink(4, hedge(8))) + MW
                    ^ -2 * (-P(0, mink(4, hedge(7))) - K(1, mink(4, hedge(7))))
                        * (-P(0, mink(4, hedge(8))) - K(1, mink(4, hedge(8)))))
                * (P(0, mink(4, edge(5, 1)))
                    + K(0, mink(4, edge(5, 1)))
                    + K(1, mink(4, edge(5, 1))))
                * g(mink(4, hedge(0)), mink(4, hedge(8)))
                * gamma(bis(4, hedge(10)), bis(4, hedge(6)), mink(4, hedge(11)))
                * gamma(bis(4, hedge(2)), bis(4, vertex(1, 1)), mink(4, hedge(7)))
                * gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge(3, 1)))
                * gamma(bis(4, hedge(9)), bis(4, hedge(10)), mink(4, edge(5, 1)))
                * projm(bis(4, hedge(5)), bis(4, hedge(1)))
                * projm(bis(4, vertex(1, 1)), bis(4, hedge(9)))
                * (1 / 2)
            ^ (1 / 2),
        default_namespace = "spenso"
    );

    validate_gamma(expr.cook_indices(), const_map.clone());
}
fn validate_gamma(expr: Atom, const_map: HashMap<Atom, symbolica::domains::float::Complex<f64>>) {
    let mut net = expr.parse_to_hep_net(&ParseSettings::default()).unwrap();

    println!("Expression: {}", expr);
    let simplified = expr.simplify_gamma();

    println!("Simplified to {}", simplified);
    let mut net_simplified = simplified
        .parse_to_hep_net(&ParseSettings::default())
        .unwrap();

    println!("{}", net_simplified.dot_pretty());
    net.execute::<Sequential, SmallestDegree, _, _, _>(&*HEP_LIB, &*FUN_LIB)
        .unwrap();
    net_simplified
        .execute::<Sequential, SmallestDegree, _, _, _>(&*HEP_LIB, &*FUN_LIB)
        .unwrap();

    let function_map = HashMap::new();

    if let ExecutionResult::Val(v) = net.result_tensor(&*HEP_LIB).unwrap() {
        if let ExecutionResult::Val(v2) = net_simplified.result_tensor(&*HEP_LIB).unwrap() {
            let mut res = v.into_owned();
            println!("{res}");
            let mut res_simplified = v2.into_owned();
            res.evaluate_complex(|c| c.into(), &const_map, &function_map);
            res_simplified.evaluate_complex(|c| c.into(), &const_map, &function_map);
            res = res.to_dense();
            res_simplified = res_simplified.to_dense();

            let mut sub = res.sub_fallible(&res_simplified).unwrap();
            sub.to_param();
            let sub = sub.try_into_parametric().unwrap();
            let zero = sub
                .zero_test(10, 0.01)
                .iter_flat()
                .fold(ConditionResult::True, |a, (_, b)| a & *b);

            match zero {
                ConditionResult::False => panic!(
                    "Should be zero but \n{}\n minus simplified\n{}\n is \n{}",
                    res, res_simplified, sub
                ),
                ConditionResult::Inconclusive => panic!("Inconclusive"),
                _ => {
                    println!("Works:res\n{}res_simplified\n{}", res, res_simplified)
                }
            }
        } else {
            panic!("Expected tensor result");
        }
    } else {
        panic!("Expected tensor result");
    }
}
