use std::time::Instant;

use idenso::{
    metric::MetricSimplifier,
    representations::initialize,
    tensor::{SymbolicNetParse, SymbolicTensor, contract::MetricSimplify},
};
use spenso::{
    network::{
        SingleSmallestDegree, SmallestDegree, SmallestDegreeIter, Steps, StructureLessDisplay,
        contract::SingleLargestDegree,
        library::{DummyLibrary, function_lib::Wrap},
        parsing::{NetworkParse, ParseSettings},
        store::TensorScalarStore,
        tags::SPENSO_TAG,
    },
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::{Match, MatchSettings, MatchStack, Pattern, Replacement},
    parse, symbol,
    transformer::Transformer,
};

use idenso::tensor::SymbolicNetExt;
fn run_network_informed() {
    initialize();
    let d = symbol!("d"; Symmetric, Linear); // metric

    let mink = Minkowski {}.new_rep(4);

    let k0 = symbol!("k0");
    let muw1_ = symbol!(
        "muw1_",
        tags = [&SPENSO_TAG.representation, &SPENSO_TAG.self_dual]
    );
    let k1_ = symbol!("k1_");
    let x___ = symbol!("x___");

    let mut r = parse!(
        "vx(1,-k(0), k(0)-k(1), k(1), k(10), spenso::mink(4,1), spenso::mink(4,8))
    vx(2,-k(1), k(2), k(1)-k(2), spenso::mink(4,1), spenso::mink(4,2), spenso::mink(4,9))
    vx(3,-k(2), k(3), k(2)-k(3), spenso::mink(4,2), spenso::mink(4,3), spenso::mink(4,10))
    vx(4,-k(3), k(4), k(3)-k(4), spenso::mink(4,3), spenso::mink(4,4), spenso::mink(4,11))
    vx(5,-k(4), k(0), k(4)-k(0), spenso::mink(4,4), k(20), spenso::mink(4,5))
    vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), spenso::mink(4,5), spenso::mink(4,11), spenso::mink(4,6))
    vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), spenso::mink(4,6), spenso::mink(4,10), spenso::mink(4,7))
    vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), spenso::mink(4,7), spenso::mink(4,9), spenso::mink(4,8))"
    );

    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let _k1_ = symbol!("k1_");
    let _muw1_ = symbol!("muw1_", tags = ["spenso::index"]);
    let _x___ = symbol!("x___");

    let _r = parse!(
        "vx(1,-k(0), k(0)-k(1), k(1), k(10), mu1, mu8)
     vx(2,-k(1), k(2), k(1)-k(2), mu1, mu2, mu9)
     vx(3,-k(2), k(3), k(2)-k(3), mu2, mu3, mu10)
     vx(4,-k(3), k(4), k(3)-k(4), mu3, mu4, mu11)
     vx(5,-k(4), k(0), k(4)-k(0), mu4, k(20), mu5)
     vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), mu5, mu11, mu6)
     vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), mu6, mu10, mu7)
     vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), mu7, mu9, mu8)"
    );

    let mut r = parse!(
        "vx(1,-k0, k0-k1, k1, k10, mu1, mu8)
    vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)
    vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)
    vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)
    vx(5,-k4, k0, k4-k0, mu4, k20, mu5)
    vx(6,-k4+k0, -k3+k4, k3-k0, mu5, mu11, mu6)
    vx(7,-k3+k0, -k2+k3, k2-k0, mu6, mu10, mu7)
    vx(8,-k2+k0, -k1+k2, k1-k0, mu7, mu9, mu8)"
    );

    // let mut r = parse!(
    //     "vx(1,-k(0), k(0)-k(1), k(1), k(10), mu1, mu8)
    // vx(2,-k(1), k(2), k(1)-k(2), mu1, mu2, mu9)
    // vx(3,-k(2), k(3), k(2)-k(3), mu2, mu3, mu10)
    // vx(4,-k(3), k(4), k(3)-k(4), mu3, mu4, mu11)
    // vx(5,-k(4), k(0), k(4)-k(0), mu4, k(20), mu5)
    // vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), mu5, mu11, mu6)
    // vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), mu6, mu10, mu7)
    // vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), mu7, mu9, mu8)"
    // );

    let gluon_rule = parse!(
        "(- spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k1_, mu2_)
                + spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k1_, mu3_)
                + spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k2_, mu1_)
                - spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k2_, mu3_)
                - spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k3_, mu1_)
                + spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k3_, mu2_)
                )"
    )
    .to_pattern();

    let rhs_subs = Pattern::Transformer(Box::new((
        Some(gluon_rule),
        vec![Transformer::Map(Box::new(move |input, _state, out| {
            *out = input.expand().simplify_metrics(); //.to_dots();
            Ok(())
        }))],
    )));

    for i in 0..8 {
        // r = r.to_dots();
        r = r
            .replace(parse!(format!(
                "vx({}, k1_, k2_, k3_, mu1_, mu2_, mu3_)",
                i + 1
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with(&rhs_subs);
        // println!("{:>}", r);
    }

    let mut a = r
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
            depth_limit: Some(1),
            take_first_term_from_sum: false,
            ..Default::default()
        })
        .unwrap();

    println!(
        "{}",
        a.graph.dot_impl(
            |i| {
                let ss = &a.store.get_scalar(i);
                format!("{}:{}", i, ss)
            },
            |k| k.display(),
            |t| {
                let tt = &a.store.get_tensor(t);
                tt.expression.to_string()
            },
            |fk| fk.to_string(),
        )
    );

    // let lib = DummyLibrary::<SymbolicTensor<AbstractIndex>>::new();

    // a.execute::<Steps<1>, SingleSmallestDegree<false, MetricSimplify>, _, _, _>(&lib, &Wrap {})
    //     .unwrap();
    // a.execute::<Steps<1>, SingleSmallestDegree<false, MetricSimplify>, _, _, _>(&lib, &Wrap {})
    //     .unwrap();

    // println!(
    //     "{}",
    //     a.graph.dot_impl(
    //         |i| {
    //             let ss = &a.store.get_scalar(i);
    //             format!("{}:{}", i, ss)
    //         },
    //         |k| k.display(),
    //         |t| {
    //             let tt = &a.store.get_tensor(t);
    //             tt.expression.to_string()
    //         },
    //         |fk| fk.to_string(),
    //     )
    // );

    // println!("",

    let res = a.simple_execute::<MetricSimplify>();
}

#[test]
fn network_informed() {
    run_network_informed();
}
