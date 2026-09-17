use idenso::{dirac::GammaSimplifier, shorthands::schoonschip::Schoonschip};
use symbolica::{
    atom::{Atom, AtomView},
    coefficient::CoefficientView,
    function, parse_lit, symbol,
};

use crate::{
    initialisation::test_initialise,
    numerator::aind::Aind,
    utils::{FUN_LIB, TENSORLIB},
};
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::{
    iterators::IteratableTensor,
    network::{
        library::{
            FunctionLibrary,
            function_lib::INBUILTS,
            symbolic::{ETS, ExplicitKey},
        },
        parsing::ShadowedStructure,
    },
    structure::{
        HasStructure, ScalarTensor,
        representation::{LibraryRep, RepName},
    },
    tensors::parametric::ParamTensor,
};

#[test]
fn symbolic_tensor_functions_preserve_unregistered_functions_and_conjugation() {
    test_initialise().unwrap();
    let value = Atom::num(2) + Atom::i();
    let tensor = ParamTensor::<ShadowedStructure<Aind>>::new_scalar(value.clone());
    let function = symbol!("tensor_library_test::component_function");
    let wrapped = FUN_LIB.apply(&function, tensor.clone()).unwrap();
    assert_eq!(wrapped.scalar(), Some(function!(function, value)));

    let conjugated = FUN_LIB.apply(&INBUILTS.conj, tensor).unwrap();
    assert_eq!(conjugated.scalar(), Some(Atom::num(2) - Atom::i()));
}

#[test]
fn symbolic_hep_library_keeps_metric_and_identity_coefficients_exact() {
    test_initialise().unwrap();
    let library = TENSORLIB.read().unwrap();

    // Inspect the source coefficients before evaluator preprocessing can rationalize them.
    for rep in LibraryRep::all_representations() {
        let key = ExplicitKey::<Aind>::from_iter(
            [rep.new_rep(4), rep.dual().new_rep(4)],
            ETS.metric,
            None,
        );
        let tensor = library.get(&key.structure).unwrap();
        for (_, coefficient) in tensor.iter_flat() {
            let AtomView::Num(number) = coefficient else {
                panic!("non-numeric metric coefficient for {rep:?}: {coefficient}");
            };
            assert!(
                matches!(number.get_coeff_view(), CoefficientView::Natural(..)),
                "inexact metric coefficient for {rep:?}: {coefficient}"
            );
        }
    }
}

#[test]
fn algebra() {
    test_initialise().unwrap();

    let expr = parse_lit!(
        UFO::GC_11
            ^ 2 * Q(7, spenso::mink(dim, edge(7, 1)))
                * spenso::g(spenso::mink(dim, hedge(14)), spenso::mink(dim, hedge(15)))
                * spenso::g(spenso::coad(8, hedge(14)), spenso::coad(8, hedge(15)))
                * spenso::g(
                    spenso::dind(spenso::cof(3, hedge(13))),
                    spenso::cof(3, hedge(12))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(4)),
                    spenso::bis(4, hedge(13)),
                    spenso::mink(dim, hedge(15))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(12)),
                    spenso::bis(4, hedge(9)),
                    spenso::mink(dim, hedge(14))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(13)),
                    spenso::bis(4, hedge(12)),
                    spenso::mink(dim, edge(7, 1))
                )
                * spenso::t(
                    spenso::coad(8, hedge(14)),
                    spenso::cof(3, hedge(9)),
                    spenso::dind(spenso::cof(3, hedge(12)))
                )
                * spenso::t(
                    spenso::coad(8, hedge(15)),
                    spenso::cof(3, hedge(13)),
                    spenso::dind(spenso::cof(3, hedge(4)))
                )
    );

    println!("{}", expr.log_print(Some(120)));

    println!("{}", expr.schoonschip().log_print(Some(120)));
    println!(
        "{}",
        expr.schoonschip().simplify_gamma().log_print(Some(120))
    );
    println!(
        "{}",
        expr.schoonschip_net::<Aind>()
            .simplify_gamma()
            .log_print(Some(120))
    );
}
