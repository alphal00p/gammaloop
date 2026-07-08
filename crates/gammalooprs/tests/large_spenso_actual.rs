use gammalooprs::{
    initialisation::test_initialise,
    numerator::{ParsingNet, aind::Aind},
    utils::{F, FUN_LIB, GS, TENSORLIB},
};
use spenso::{
    network::{
        ExecutionResult, MinResultRank,
        parsing::{ParseSettings, ShadowedStructure, ShorthandParsing, StructureInferenceMode},
    },
    structure::{HasStructure, TensorStructure},
    tensors::parametric::MixedTensor,
};
use symbolica::{
    atom::Atom, function, parser::ParseSettings as SymbolicaParseSettings, wrap_input,
};

type ActualTensor = MixedTensor<F<f64>, ShadowedStructure<Aind>>;

fn parse_inline_expression(input: &str) -> Atom {
    Atom::parse_with_default_namespace(wrap_input!(input), SymbolicaParseSettings::default())
        .expect("inline diagnostic expression should parse")
}

fn parse_actual_net(label: &str, expr: &Atom) -> ParsingNet {
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque {
            inference: StructureInferenceMode::Fast,
        },
        ..Default::default()
    };
    let lib = TENSORLIB.read().unwrap();
    ParsingNet::try_from_view(expr.as_view(), &*lib, &settings)
        .unwrap_or_else(|error| panic!("{label} actual parse failed: {error}"))
}

fn execute_actual_net_min_result_rank_parallel_tensor(
    label: &str,
    mut net: ParsingNet,
) -> ActualTensor {
    let lib = TENSORLIB.read().unwrap();
    net.execute_parallel::<MinResultRank, _, _, _>(&*lib, &*FUN_LIB)
        .unwrap_or_else(|error| {
            panic!("{label} parallel min-result-rank tensor execution failed: {error}")
        });

    let result = net.result_tensor(&*lib).unwrap_or_else(|error| {
        panic!("{label} parallel min-result-rank tensor result failed: {error}")
    });
    match result {
        ExecutionResult::Val(value) => value.into_owned(),
        ExecutionResult::One => panic!("{label} produced tensor result One"),
        ExecutionResult::Zero => panic!("{label} produced tensor result Zero"),
    }
}

#[test]
fn min_result_rank_disconnected_tensor_product_mwe() {
    test_initialise().expect("GammaLoop initialization should succeed");

    let mu = parse_inline_expression("spenso::mink(4,mu)");
    let rho = parse_inline_expression("spenso::mink(4,rho)");
    let expr = function!(GS.emr_mom, 1i64, mu) * function!(GS.emr_mom, 2i64, rho);
    let tensor = execute_actual_net_min_result_rank_parallel_tensor(
        "min_result_rank_disconnected_tensor_product_mwe",
        parse_actual_net("min_result_rank_disconnected_tensor_product_mwe", &expr),
    );

    assert_eq!(tensor.structure().order(), 2);
    assert_eq!(tensor.structure().size().unwrap(), 16);
}
