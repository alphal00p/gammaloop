use spenso::network::parsing::{Parse, ParseSettings, ShadowedStructure};
use spenso::network::{
    ContractScalars, ExecutionResult, Network, Sequential, SingleSmallestDegree, SmallestDegree,
    Steps, TensorOrScalarOrKey, tags::SPENSO_TAG,
};
use symbolica::atom::{AtomCore, Symbol};

use symbolica::symbol;

use spenso::network::library::DummyLibrary;
use spenso::network::library::panicing::ErroringLibrary;
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use spenso::structure::abstract_index::AbstractIndex;

// use log::trace;

use symbolica::atom::Atom;

// Generate TestSymbols with all alphabet characters and some multi-character symbols
symbol_set!(TestSymbols, TS;
    a b c d e f h i j k l m n o p q r s t u v w x y z
    A B C D E F G H I J K L M N O P Q R S T U V W X Y Z
    vbar
    ebar
    edge_1_1 edge_2_1 edge_3_1 edge_4_1 edge_5_1 edge_6_1 edge_7_1 edge_8_1 edge_9_1 edge_10_1 edge_11_1 edge_12_1 edge_13_1 edge_14_1 edge_15_1 edge_16_1 edge_17_1 edge_18_1 edge_19_1 edge_20_1
    hedge0 hedge1 hedge2 hedge3 hedge4 hedge5 hedge6 hedge7 hedge8 hedge9 hedge10 hedge11 hedge12 hedge13 hedge14 hedge15 hedge16 hedge17 hedge18 hedge19 hedge20
    hedge_0 hedge_1 hedge_2 hedge_3 hedge_4 hedge_5 hedge_6 hedge_7 hedge_8 hedge_9 hedge_10 hedge_11 hedge_12 hedge_13 hedge_14 hedge_15 hedge_16 hedge_17 hedge_18 hedge_19 hedge_20 mul
);

pub fn test_initialize() {
    initialize();
    // let dim  =
    let _ = TS.p;
}

use core::panic;

use spenso::{
    shadowing::symbolica_utils::TypstSettings,
    structure::{
        ToSymbolic,
        representation::{Euclidean, Lorentz, Minkowski, RepName},
        slot::IsAbstractSlot,
    },
};

use insta::assert_snapshot;
use symbolica::{parse, parse_lit};

use crate::representations::initialize;
use crate::symbol_set;
use crate::tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor};

#[test]
fn parse_scalar() {
    let expr = parse!("1");

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn parse_pow() {
    test_initialize();
    let _sqrt = symbol!("sqrt_scalar", tag = SPENSO_TAG.tag);
    let expr = parse!("ee^3*sqrt_scalar((m+d(spenso::mink(4,1))^2)^(-5))");

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    println!("{}", net.snapshot_dot());
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn parse_ratio() {
    test_initialize();
    let expr = parse_lit!(
        (P(1, spenso::mink(4, 1)) * P(2, spenso::mink(4, 1)))
            / f(mul(P(3, spenso::mink(4, 1)) * P(4, spenso::mink(4, 1))))
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    assert_snapshot!(net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "∏"];
      1	 [label = "S:(f(mul(P(3,mink(4,1))*P(4,mink(4,1)))))^(-1)"];
      2	 [label = "T:P(1,mink(4,1))"];
      3	 [label = "T:P(2,mink(4,1))"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      3:7:s	-> 0:1:s	 [id=1  color="red:blue;0.5"];
      2:5:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
      1:4:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
      2:6:s	-> 3:8:s	 [id=4 dir=none  color="red:blue;0.5" label="mink4|1"];
    }
    "#);

    assert_eq!(net.simple_execute::<()>(), expr);
}
#[test]
fn parse_div() {
    let expr = parse!(
        "c*a/bracket(d(spenso::mink(4,1))* b(spenso::mink(4,1)))(bracket(d(spenso::mink(4,1))* b(spenso::mink(4,1)))^-3)(bracket(d(spenso::mink(4,1))* b(spenso::mink(4,1)))^-2)"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_snapshot!(net.simple_execute::<()>().to_bare_ordered_string(), @"(bracket(b(mink(4,1))*d(mink(4,1))))^(-6)*a*c");

    let expr = parse_lit!(st(Q(1, mink(4, 1)) * Q(2, mink(4, 1))) ^ -1);
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    println!("{}", net.snapshot_dot());
    assert_snapshot!(net.simple_execute::<()>().to_bare_ordered_string(), @"(st(Q(1,mink(4,1))*Q(2,mink(4,1))))^(-1)");
}

#[test]
fn parse_scalar_tensor() {
    test_initialize();
    let expr = parse!("(
            (
                  -1*gammalooprs::{}::mUV^2+gammalooprs::{}::Q(6,spenso::mink(4,gammalooprs::{}::uv_mink_1337))
                    *gammalooprs::{}::Q(7,spenso::mink(4,gammalooprs::{}::uv_mink_1337))
                 )
            )*2");
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    println!("{}", net.dot_pretty());
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn parse_val() {
    test_initialize();
    let expr = parse_lit!(
        ((Q(5, cind(0)))
            ^ 2 + (Q(5, cind(1)))
            ^ 2 * -1 + (Q(5, cind(2)))
            ^ 2 * -1 + (Q(5, cind(3)))
            ^ 2 * -1)
            ^ (-1)
                * ((Q(6, cind(0)))
                    ^ 2 + (Q(6, cind(1)))
                    ^ 2 * -1 + (Q(6, cind(2)))
                    ^ 2 * -1 + (Q(6, cind(3)))
                    ^ 2 * -1)
            ^ (-1) * 1𝑖 / 3
            ^ 3 * Q(5, mink(4, edge_5_1))
                * Q(6, mink(4, edge_6_1))
                * u(1, bis(4, hedge_1))
                * vbar(2, bis(4, hedge_2))
                * ebar(0, mink(4, hedge_0))
                * ebar(3, mink(4, hedge_3))
                * ebar(4, mink(4, hedge_4))
                * g(cof(3, hedge_1), dind(cof(3, hedge_2)))
                * g(cof(3, hedge_2), dind(cof(3, hedge_1)))
                * gamma(bis(4, hedge_2), bis(4, hedge_8), mink(4, hedge_0))
                * gamma(bis(4, hedge_5), bis(4, hedge_1), mink(4, hedge_3))
                * gamma(bis(4, hedge_6), bis(4, hedge_5), mink(4, edge_5_1))
                * gamma(bis(4, hedge_7), bis(4, hedge_6), mink(4, hedge_4))
                * gamma(bis(4, hedge_8), bis(4, hedge_7), mink(4, edge_6_1)),
        default_namespace = "spenso"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_snapshot!(net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "∏"];
      1	 [label = "S:((Q(5,cind(0)))^2+(Q(5,cind(1)))^2*-1+(Q(5,cind(2)))^2*-1+(Q(5,cind(3)))^2*-1)^(-1)*((Q(6,cind(0)))^2+(Q(6,cind(1)))^2*-1+(Q(6,cind(2)))^2*-1+(Q(6,cind(3)))^2*-1)^(-1)*1𝑖/27"];
      2	 [label = "T:g(cof(3,hedge_1),dind(cof(3,hedge_2)))"];
      3	 [label = "T:g(cof(3,hedge_2),dind(cof(3,hedge_1)))"];
      4	 [label = "T:gamma(bis(4,hedge_2),bis(4,hedge_8),mink(4,hedge_0))"];
      5	 [label = "T:gamma(bis(4,hedge_6),bis(4,hedge_5),mink(4,edge_5_1))"];
      6	 [label = "T:gamma(bis(4,hedge_7),bis(4,hedge_6),mink(4,hedge_4))"];
      7	 [label = "T:gamma(bis(4,hedge_5),bis(4,hedge_1),mink(4,hedge_3))"];
      8	 [label = "T:gamma(bis(4,hedge_8),bis(4,hedge_7),mink(4,edge_6_1))"];
      9	 [label = "T:Q(5,mink(4,edge_5_1))"];
      10	 [label = "T:Q(6,mink(4,edge_6_1))"];
      11	 [label = "T:u(1,bis(4,hedge_1))"];
      12	 [label = "T:vbar(2,bis(4,hedge_2))"];
      13	 [label = "T:ebar(0,mink(4,hedge_0))"];
      14	 [label = "T:ebar(3,mink(4,hedge_3))"];
      15	 [label = "T:ebar(4,mink(4,hedge_4))"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      15:55:s	-> 0:1:s	 [id=1  color="red:blue;0.5"];
      14:53:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
      13:51:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
      12:49:s	-> 0:4:s	 [id=4  color="red:blue;0.5"];
      11:47:s	-> 0:5:s	 [id=5  color="red:blue;0.5"];
      10:45:s	-> 0:6:s	 [id=6  color="red:blue;0.5"];
      9:43:s	-> 0:7:s	 [id=7  color="red:blue;0.5"];
      8:39:s	-> 0:8:s	 [id=8  color="red:blue;0.5"];
      7:35:s	-> 0:9:s	 [id=9  color="red:blue;0.5"];
      6:31:s	-> 0:10:s	 [id=10  color="red:blue;0.5"];
      5:27:s	-> 0:11:s	 [id=11  color="red:blue;0.5"];
      4:23:s	-> 0:12:s	 [id=12  color="red:blue;0.5"];
      3:20:s	-> 0:13:s	 [id=13  color="red:blue;0.5"];
      2:17:s	-> 0:14:s	 [id=14  color="red:blue;0.5"];
      1:16:s	-> 0:15:s	 [id=15  color="red:blue;0.5"];
      2:19:s	-> 3:21:s	 [id=16 dir=back  color="red:blue;0.5" label="cof🠓3|_hedge_2"];
      2:18:s	-> 3:22:s	 [id=17  color="red:blue;0.5" label="cof🠑3|^hedge_1"];
      4:26:s	-> 13:52:s	 [id=18 dir=none  color="red:blue;0.5" label="mink4|hedge_0"];
      4:24:s	-> 12:50:s	 [id=19 dir=none  color="red:blue;0.5" label="bis4|hedge_2"];
      4:25:s	-> 8:41:s	 [id=20 dir=none  color="red:blue;0.5" label="bis4|hedge_8"];
      5:30:s	-> 9:44:s	 [id=21 dir=none  color="red:blue;0.5" label="mink4|edge_5_1"];
      5:28:s	-> 6:32:s	 [id=22 dir=none  color="red:blue;0.5" label="bis4|hedge_6"];
      5:29:s	-> 7:37:s	 [id=23 dir=none  color="red:blue;0.5" label="bis4|hedge_5"];
      6:33:s	-> 8:40:s	 [id=24 dir=none  color="red:blue;0.5" label="bis4|hedge_7"];
      6:34:s	-> 15:56:s	 [id=25 dir=none  color="red:blue;0.5" label="mink4|hedge_4"];
      7:38:s	-> 14:54:s	 [id=26 dir=none  color="red:blue;0.5" label="mink4|hedge_3"];
      7:36:s	-> 11:48:s	 [id=27 dir=none  color="red:blue;0.5" label="bis4|hedge_1"];
      8:42:s	-> 10:46:s	 [id=28 dir=none  color="red:blue;0.5" label="mink4|edge_6_1"];
    }
    "#);
    assert_eq!(net.simple_execute::<()>(), expr);
    let mut out = String::new();
    expr.typst_fmt(&mut out, &TypstSettings::lowering())
        .unwrap();
    println!("{}", out)
}

#[test]
fn parse_scalar_tensors_step_by() {
    initialize();
    let expr = parse!(
        "c*a*b(spenso::mink(4,1))*d(spenso::mink(4,2))*d(spenso::mink(4,1))*d(spenso::mink(4,2))"
    );

    let mut net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    let lib = DummyLibrary::<_>::new();
    let fnlib = ErroringLibrary::<Symbol>::new();

    let mut netc = net.clone();
    assert_snapshot!(
        net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "∏"];
      1	 [label = "S:a*c"];
      2	 [label = "T:b(mink(4,1))"];
      3	 [label = "T:d(mink(4,1))"];
      4	 [label = "^( 2 )"];
      5	 [label = "T:d(mink(4,2))"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      4:10:s	-> 0:1:s	 [id=1  color="red:blue;0.5"];
      3:8:s	-> 0:2:s	 [id=2  color="red:blue;0.5"];
      2:6:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
      1:5:s	-> 0:4:s	 [id=4  color="red:blue;0.5"];
      2:7:s	-> 3:9:s	 [id=5 dir=none  color="red:blue;0.5" label="mink4|1"];
      5:14:s	-> 4:12:s	 [id=6 dir=none  color="red:blue;0.5" label="mink4|2"];
      5:13:s	-> 4:11:s	 [id=7  color="red:blue;0.5"];
    }
    "#
    );
    net.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
        .unwrap();
    assert_snapshot!(
        net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "∏"];
      1	 [label = "S:a*c"];
      2	 [label = "T:d(mink(4,1))"];
      3	 [label = "T:b(mink(4,1))"];
      4	 [label = "S:(d(mink(4,2)))^2"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      4:10:s	-> 0:2:s	 [id=1  color="red:blue;0.5"];
      1:5:s	-> 0:4:s	 [id=2  color="red:blue;0.5"];
      2:6:s	-> 0:3:s	 [id=3  color="red:blue;0.5"];
      2:7:s	-> 3:8:s	 [id=4 dir=none  color="red:blue;0.5" label="mink4|1"];
      3:9:s	-> 0:1:s	 [id=5  color="red:blue;0.5"];
    }
    "#
    );
    net.execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
        .unwrap();
    assert_snapshot!(
        net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      3	 [label = "∏"];
      0	 [label = "S:(d(mink(4,2)))^2"];
      1	 [label = "S:a*c"];
      2	 [label = "T:b(mink(4,1))*d(mink(4,1))"];
      ext0	 [style=invis];
      3:0:s	-> ext0	 [id=0 color="red"];
      1:5:s	-> 3:4:s	 [id=1  color="red:blue;0.5"];
      2:6:s	-> 3:3:s	 [id=2  color="red:blue;0.5"];
      0:1:s	-> 3:2:s	 [id=3  color="red:blue;0.5"];
    }
    "#
    );
    net.execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
        .unwrap();
    assert_snapshot!(
        net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      3	 [label = "∏"];
      0	 [label = "S:(d(mink(4,2)))^2"];
      1	 [label = "S:a*c"];
      2	 [label = "T:b(mink(4,1))*d(mink(4,1))"];
      ext0	 [style=invis];
      3:0:s	-> ext0	 [id=0 color="red"];
      0:1:s	-> 3:2:s	 [id=1  color="red:blue;0.5"];
      1:5:s	-> 3:4:s	 [id=2  color="red:blue;0.5"];
      2:6:s	-> 3:3:s	 [id=3  color="red:blue;0.5"];
    }
    "#
    );
    net.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
        .unwrap();
    assert_snapshot!(
        net.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "T:(d(mink(4,2)))^2*a*b(mink(4,1))*c*d(mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
        }
        "#
    );
    netc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
        .unwrap();
    assert_snapshot!(
        netc.snapshot_dot(),@r#"
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "T:(d(mink(4,2)))^2*a*b(mink(4,1))*c*d(mink(4,1))"];
          ext0	 [style=invis];
          0:0:s	-> ext0	 [id=0 color="red"];
        }
        "#
    );
    if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) = net.result().unwrap()
    {
        if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor {
            tensor: tensor2, ..
        }) = netc.result().unwrap()
        {
            assert_eq!(tensor2.expression, tensor.expression);
        } else {
            panic!("Not scalar")
        }
    } else {
        panic!("Not scalar")
    }
}

#[test]
fn parse_scalar_expr() {
    let expr = parse!("(y+x(spenso::mink(4,1))*y(spenso::mink(4,1))) *(1+1+2*x*(3*sin(r))/t)");
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn parse_tensor_expr() {
    let tensor1 = ShadowedStructure::<AbstractIndex>::from_iter(
        [
            Minkowski {}.new_slot(4, 1).to_lib(),
            Euclidean {}.new_slot(4, 2).to_lib(),
            Minkowski {}.new_slot(symbol!("d"), symbol!("mu")).to_lib(),
        ],
        symbol!("T"),
        None,
    )
    .structure
    .to_symbolic(None)
    .unwrap();

    let tensor2 = ShadowedStructure::<AbstractIndex>::from_iter(
        [
            Minkowski {}.new_slot(4, 1).to_lib(),
            Lorentz {}.new_slot(7, 1).to_lib(),
            Lorentz {}.new_slot(3, 2).to_lib(),
        ],
        symbol!("TT"),
        None,
    )
    .structure
    .to_symbolic(None)
    .unwrap();

    let tensor3 = ShadowedStructure::<AbstractIndex>::from_iter(
        [
            Lorentz {}.dual().new_slot(7, 1).to_lib(),
            Euclidean {}.new_slot(4, 2).to_lib(),
            Minkowski {}.new_slot(symbol!("d"), symbol!("mu")).to_lib(),
        ],
        symbol!("TTT"),
        None,
    )
    .structure
    .to_symbolic(None)
    .unwrap();

    let tensor4 = ShadowedStructure::<AbstractIndex>::from_iter(
        [Lorentz {}.new_slot(3, 2).to_lib()],
        symbol!("L"),
        None,
    )
    .structure
    .to_symbolic(None)
    .unwrap();

    let tensor5 = ShadowedStructure::<AbstractIndex>::from_iter(
        [Lorentz {}.dual().new_slot(3, 2).to_lib()],
        symbol!("P"),
        None,
    )
    .structure
    .to_symbolic(None)
    .unwrap();

    let expr = (parse!("a*sin(x/2)") * tensor1 * tensor2 * tensor3 + tensor4) * tensor5;

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

// -G^2*(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7))+g(mink(4,5),mink(4,7))*Q(2,mink(4,6))+g(mink(4,5),mink(4,7))*Q(4,mink(4,6))-g(mink(4,6),mink(4,7))*Q(3,mink(4,5))-g(mink(4,6),mink(4,7))*Q(4,mink(4,5)))*id(mink(4,2),mink(4,5))*id(mink(4,3),mink(4,6))*id(euc(4,0),euc(4,5))*id(euc(4,1),euc(4,4))*g(mink(4,4),mink(4,7))*vbar(1,euc(4,1))*u(0,euc(4,0))*ebar(2,mink(4,2))*ebar(3,mink(4,3))*gamma(mink(4,4),euc(4,5),euc(4,4))

#[test]
fn parse_big_tensors() {
    initialize();
    let expr = parse!(
        "-G^2*(-g(mink(4,5),mink(4,6))*Q(2,mink(4,7))+g(mink(4,5),mink(4,6))*Q(3,mink(4,7))+g(mink(4,5),mink(4,7))*Q(2,mink(4,6))+g(mink(4,5),mink(4,7))*Q(4,mink(4,6))-g(mink(4,6),mink(4,7))*Q(3,mink(4,5))-g(mink(4,6),mink(4,7))*Q(4,mink(4,5)))*g(mink(4,2),mink(4,5))*g(mink(4,3),mink(4,6))*g(euc(4,0),euc(4,5))*g(euc(4,1),euc(4,4))*g(mink(4,4),mink(4,7))*vbar(1,euc(4,1))*u(0,euc(4,0))*ebar(2,mink(4,2))*ebar(3,mink(4,3))*gamma(euc(4,5),euc(4,4),mink(4,4))",
        default_namespace = "spenso"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn equal_duals() {
    initialize();
    let expr = parse_lit!(
        ((Q(5, cind(0)))
            ^ 2 + (Q(5, cind(1)))
            ^ 2 * -1 + (Q(5, cind(2)))
            ^ 2 * -1 + (Q(5, cind(3)))
            ^ 2 * -1)
            ^ (-1)
                * ((Q(6, cind(0)))
                    ^ 2 + (Q(6, cind(1)))
                    ^ 2 * -1 + (Q(6, cind(2)))
                    ^ 2 * -1 + (Q(6, cind(3)))
                    ^ 2 * -1)
            ^ (-1)
                * ((Q(7, cind(0)))
                    ^ 2 + (Q(7, cind(1)))
                    ^ 2 * -1 + (Q(7, cind(2)))
                    ^ 2 * -1 + (Q(7, cind(3)))
                    ^ 2 * -1)
            ^ (-1)
                * ((Q(8, cind(0)))
                    ^ 2 + (Q(8, cind(1)))
                    ^ 2 * -1 + (Q(8, cind(2)))
                    ^ 2 * -1 + (Q(8, cind(3)))
                    ^ 2 * -1)
            ^ (-1)
                * ((Q(9, cind(0)))
                    ^ 2 + (Q(9, cind(1)))
                    ^ 2 * -1 + (Q(9, cind(2)))
                    ^ 2 * -1 + (Q(9, cind(3)))
                    ^ 2 * -1)
            ^ (-1) * -1𝑖 / 3 * UFO::GC_11
            ^ 2 * UFO::GC_1
            ^ 3 * UFO::Gamma(mink(4, hedge_13), bis(4, hedge_11), 2)
                * Q(5, mink(4, edge_5_1))
                * Q(6, mink(4, edge_6_1))
                * Q(7, mink(4, edge_7_1))
                * Q(8, mink(4, edge_8_1))
                * u(1)
                * vbar(2, bis(4, hedge_2))
                * ebar(0, mink(4, hedge_0))
                * ebar(3, mink(4, hedge_3))
                * ebar(4, mink(4, hedge_4))
                * g(dind(cof(3, hedge_1)), cof(3, hedge_2))
                * g(mink(4, hedge_13), mink(4, hedge_14))
                * gamma(bis(4, hedge_10), bis(4, hedge_9), mink(4, edge_8_1))
                * gamma(bis(4, hedge_12), bis(4, hedge_11), mink(4, edge_7_1))
                * gamma(bis(4, hedge_2), bis(4, hedge_10), mink(4, hedge_14))
                * gamma(bis(4, hedge_5), bis(4, hedge_12), mink(4, hedge_3))
                * gamma(bis(4, hedge_6), bis(4, hedge_5), mink(4, edge_5_1))
                * gamma(bis(4, hedge_7), bis(4, hedge_6), mink(4, hedge_4))
                * gamma(bis(4, hedge_8), bis(4, hedge_7), mink(4, edge_6_1))
                * gamma(bis(4, hedge_9), bis(4, hedge_8), mink(4, hedge_0))
                * t(coad(8, hedge_13), cof(2), dind(cof(3, hedge_11)))
                * t(coad(8, hedge_13), cof(3, hedge_11), dind(cof(3, hedge_2))),
        default_namespace = "spenso"
    );

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn gammaloop_six_photon() {
    initialize();
    let expr = parse_lit!(
        -64 / 243 * ee
            ^ 6 * (MT * g(euc(4, hedge3), euc(4, hedge4))
                + gamma(euc(4, hedge3), euc(4, hedge4), mink(4, edge_3_1))
                    * Q(3, mink(4, edge_3_1)))
                * (MT * g(euc(4, hedge6), euc(4, hedge7))
                    + gamma(euc(4, hedge6), euc(4, hedge7), mink(4, edge_5_1))
                        * Q(5, mink(4, edge_5_1)))
                * (MT * g(euc(4, hedge9), euc(4, hedge10))
                    + gamma(euc(4, hedge9), euc(4, hedge10), mink(4, edge_7_1))
                        * Q(7, mink(4, edge_7_1)))
                * (MT * g(euc(4, hedge11), euc(4, hedge12))
                    + gamma(euc(4, hedge11), euc(4, hedge12), mink(4, edge_8_1))
                        * Q(8, mink(4, edge_8_1)))
                * (MT * g(euc(4, hedge13), euc(4, hedge14))
                    + gamma(euc(4, hedge13), euc(4, hedge14), mink(4, edge_9_1))
                        * Q(9, mink(4, edge_9_1)))
                * (MT * g(euc(4, hedge16), euc(4, hedge17))
                    + gamma(euc(4, hedge16), euc(4, hedge17), mink(4, edge_11_1))
                        * Q(11, mink(4, edge_11_1)))
                * gamma(euc(4, hedge4), euc(4, hedge6), mink(4, hedge5))
                * gamma(euc(4, hedge7), euc(4, hedge9), mink(4, hedge8))
                * gamma(euc(4, hedge10), euc(4, hedge11), mink(4, hedge0))
                * gamma(euc(4, hedge12), euc(4, hedge13), mink(4, hedge1))
                * gamma(euc(4, hedge14), euc(4, hedge16), mink(4, hedge15))
                * gamma(euc(4, hedge17), euc(4, hedge3), mink(4, hedge2))
                * eps(0, mink(4, hedge0))
                * eps(1, mink(4, hedge1))
                * epsbar(2, mink(4, hedge2))
                * epsbar(4, mink(4, hedge5))
                * epsbar(6, mink(4, hedge8))
                * epsbar(10, mink(4, hedge15)),
        default_namespace = "spenso"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn parse_neg_tensors() {
    initialize();
    let expr = parse!(
        "-d(mink(4,6),mink(4,5))*Q(2,mink(4,7))+d(mink(4,6),mink(4,5))*Q(3,mink(4,7))",
        default_namespace = "spenso"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn many_sums() {
    initialize();
    let expr = parse_lit!(
        (P(4, mink(4, r_2)) + N(4, mink(4, r_2)))
            * (P(5, mink(4, r_3)) + N(5, mink(4, r_3)))
            * (A(mink(4, r_2), mink(4, r_3)) + B(mink(4, r_3), mink(4, r_2))),
        default_namespace = "spenso"
    );
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn contract_problem() {
    initialize();

    let expr = parse_lit!(
        (-1 * Q(EMRID(0, 4), mink(4, l_20))
            * gamma(euc(4, l_3), euc(4, l_6), mink(4, l_0))
            * gamma(euc(4, l_5), euc(4, l_2), mink(4, l_20))
            * gamma(euc(4, l_6), euc(4, l_5), mink(4, l_1))
            + 2 * Q(EMRID(0, 4), mink(4, l_1)) * gamma(euc(4, l_3), euc(4, l_2), mink(4, l_0)))
            * 1𝑖
            * G
            ^ 2,
        default_namespace = "spenso"
    );

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
// #[should_panic]
fn parse_problem() {
    initialize();
    let expr = parse_lit!(
        1 / 18 * ee
            ^ 2 * G
            ^ 4 * ((mUV
                ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                    + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (1 / 2)
                    + (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (1 / 2))
            ^ -1 * ((mUV
                ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                    + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (1 / 2)
                    + (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (1 / 2))
            ^ -1 * (-Q3(3, mink(4, edge_4_1))
                + (mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                ^ (1 / 2) * sigma(4) * delta(cind(0), mink(4, edge_4_1)))
                * (-Q3(3, mink(4, edge_5_1))
                    + (mUV
                        ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                            + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
                    ^ (1 / 2) * sigma(5) * delta(cind(0), mink(4, edge_5_1)))
                * (mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(3, Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
            ^ (-1 / 2)
                * (mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(4, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
            ^ (-1 / 2)
                * (mUV
                    ^ 2 - (mUV ^ 2 - dot(Q3(3), Q3(3))) - dot(Q3(3), Q3(3))
                        + OSE(5, -Q3(3), mUV ^ 2, mUV ^ 2 - dot(Q3(3), Q3(3))))
            ^ (-1 / 2)
                * g(mink(4, hedge_3), mink(4, hedge_4))
                * gamma(bis(4, hedge(1)), bis(4, hedge(8)), mink(4, hedge_4))
            ^ 2 * gamma(bis(4, hedge(5)), bis(4, hedge(0)), mink(4, hedge_3))
            ^ 2 * gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge_4_1))
                * gamma(bis(4, hedge(7)), bis(4, hedge(6)), mink(4, hedge_2))
            ^ 2 * gamma(bis(4, hedge(8)), bis(4, hedge(7)), mink(4, edge_5_1)),
        default_namespace = "spenso"
    );

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}
#[test]
// #[should_panic]
fn infinite_execution() {
    test_initialize();
    let _s = symbol!("spenso::dim");

    let _ = parse_lit!(
        (N(0, mink(dim, l(0))) * P(0, mink(dim, r(0)))
            + N(0, mink(dim, r(0))) * P(0, mink(dim, l(0))))
            * N(0, mink(dim, dummy_ss(0, 1)))
            * P(0, mink(dim, dummy_ss(0, 1)))
            + -1 * N(0, mink(dim, dummy_ss(0, 3)))
                * N(0, mink(dim, dummy_ss(0, 4)))
                * P(0, mink(dim, dummy_ss(0, 3)))
                * P(0, mink(dim, dummy_ss(0, 4)))
                * g(mink(dim, l(0)), mink(dim, r(0))),
        default_namespace = "spenso"
    )
    .replace(parse_lit!(dim))
    .with(Atom::num(4));

    let expr = parse_lit!(
        -1𝑖 * G
            ^ 3 * (g(mink(4, l_6), mink(4, l_8)) * g(mink(4, l_7), mink(4, l_9))
                - g(mink(4, l_6), mink(4, l_9)) * g(mink(4, l_8), mink(4, l_7)))
                * g(mink(4, l_0), mink(4, l_6))
                * g(mink(4, l_1), mink(4, l_7))
                * g(mink(4, l_4), mink(4, l_8))
                * g(mink(4, l_5), mink(4, l_9))
                * g(bis(4, l_2), bis(4, l_5))
                * g(bis(4, l_3), bis(4, l_6))
                * gamma(bis(4, l_6), bis(4, l_5), mink(4, l_5)),
        default_namespace = "spenso"
    );

    let lib = DummyLibrary::<_>::new();
    let fnlib = ErroringLibrary::<Symbol>::new();
    let mut net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();

    let mut net_iter = net.clone();

    loop {
        let old = net_iter.clone();
        net_iter
            .execute::<Steps<1>, SingleSmallestDegree<false>, _, _, _>(&lib, &fnlib)
            .unwrap();
        net_iter
            .execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
            .unwrap();
        if net_iter == old {
            break;
        }
    }
    assert_snapshot!(
        net_iter.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "T:(-1*g(mink(4,l_6),mink(4,l_9))*g(mink(4,l_7),mink(4,l_8))+g(mink(4,l_6),mink(4,l_8))*g(mink(4,l_7),mink(4,l_9)))*-1𝑖*G^3*g(bis(4,l_2),bis(4,l_5))*g(bis(4,l_3),bis(4,l_6))*g(mink(4,l_0),mink(4,l_6))*g(mink(4,l_1),mink(4,l_7))*g(mink(4,l_4),mink(4,l_8))*g(mink(4,l_5),mink(4,l_9))*gamma(bis(4,l_6),bis(4,l_5),mink(4,l_5))"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      ext1	 [style=invis];
      0:1:s	-> ext1	 [id=1 dir=none color="red" label="mink4|l_4"];
      ext2	 [style=invis];
      0:2:s	-> ext2	 [id=2 dir=none color="red" label="bis4|l_3"];
      ext3	 [style=invis];
      0:3:s	-> ext3	 [id=3 dir=none color="red" label="mink4|l_0"];
      ext4	 [style=invis];
      0:4:s	-> ext4	 [id=4 dir=none color="red" label="mink4|l_1"];
      ext5	 [style=invis];
      0:5:s	-> ext5	 [id=5 dir=none color="red" label="bis4|l_2"];
    }
    "#
    );
    net.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
        .unwrap();
}

#[test]
fn gammaloop_input() {
    initialize();
    let expr = parse_lit!(
        16 / 81 * ee
            ^ 4 * (MT * g(euc(4, hedge_3), euc(4, hedge_4))
                    + gamma(euc(4, hedge_3), euc(4, hedge_4), mink(4, edge_3_1))
                        * Q(3, mink(4, edge_3_1)))
                    * (MT * g(euc(4, hedge_6), euc(4, hedge_7))
                    + gamma(euc(4, hedge_6), euc(4, hedge_7), mink(4, edge_5_1))
                            * Q(5, mink(4, edge_5_1)))
                    * (MT * g(euc(4, hedge_10), euc(4, hedge_11))
                    + gamma(euc(4, hedge_10), euc(4, hedge_11), mink(4, edge_8_1))
                            * Q(8, mink(4, edge_8_1)))
                    // * (-1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, 1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, -1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, 1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (-Q(1, cind(0)) + Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, -1, 1, 1, 1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, 1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (-Q(1, cind(0)) + Q(2, cind(0)) + Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, 1, 1, -1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (-Q(1, cind(0))
                    //         + Q(2, cind(0))
                    //         + Q(4, cind(0))
                    //         + Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, 1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1 - 1 / 8 * (Q(1, cind(0)) - Q(7, cind(0)) + OSE(8))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(7, cind(0)) + OSE(3))
                    //     ^ -1 * (Q(1, cind(0)) - Q(2, cind(0)) - Q(4, cind(0)) - Q(7, cind(0))
                    //         + OSE(5))
                    //     ^ -1 * delta_sigma(1, 1, 1, -1, 1, -1, 1, 1, -1) * OSE(3)
                    //     ^ -1 * OSE(5)
                    //     ^ -1 * OSE(8)
                    //     ^ -1)
                    * g(dind(cof(3, hedge_8)), cof(3, hedge_1))
                    * gamma(euc(4, hedge_1), euc(4, hedge_10), mink(4, hedge_9))
                    * gamma(euc(4, hedge_4), euc(4, hedge_6), mink(4, hedge_5))
                    * gamma(euc(4, hedge_7), euc(4, hedge_8), mink(4, hedge_0))
                    * gamma(euc(4, hedge_11), euc(4, hedge_3), mink(4, hedge_2))
                    * ubar(6, euc(4, hedge_8))
                    * u(1, euc(4, hedge_1))
                    * eps(0, mink(4, hedge_0))
                    * eps(2, mink(4, hedge_2))
                    * eps(4, mink(4, hedge_5))
                    * eps(7, mink(4, hedge_9)),
        default_namespace = "spenso"
    );

    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(net.simple_execute::<()>(), expr);
}
#[test]
fn wrapping() {
    initialize();
    let expr = parse_lit!(
        A * g(euc(4, hedge_3), euc(4, hedge_5)),
        default_namespace = "spenso"
    );

    let expr2 = parse_lit!(
        B * gg(euc(4, hedge_3), euc(4, hedge_4)),
        default_namespace = "spenso"
    );
    let expr3 = parse_lit!(
        C * ggg(euc(4, hedge_5), euc(4, hedge_4)),
        default_namespace = "spenso"
    );

    let lib = DummyLibrary::<SymbolicTensor>::new();
    let fnlib = ErroringLibrary::<Symbol>::new();
    let settings = &ParseSettings::default();
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(settings)
        .unwrap();
    let net2 = expr2
        .parse_to_symbolic_net::<AbstractIndex>(settings)
        .unwrap();
    let net3 = expr3
        .parse_to_symbolic_net::<AbstractIndex>(settings)
        .unwrap();

    let mut acc = Network::one();
    println!("{}", expr);

    acc *= net;
    acc *= net2;
    acc *= net3;
    println!("{}", acc.snapshot_dot());

    acc.merge_ops();

    println!("{}", acc.snapshot_dot());

    // acc.execute::<Sequential, SmallestDegree, _, _,_>(&lib,&fnlib)
    //     .unwrap();

    acc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
        .unwrap();
    println!("{}", acc.snapshot_dot());
    let obt: Atom = acc.result_scalar().unwrap().into();
    assert_eq!(obt, expr * expr2 * expr3)
}
#[test]
fn scalar_mult() {
    initialize();
    let expr = parse_lit!(A * B * C * g(euc(4, hedge_3)), default_namespace = "spenso");

    let expr2 = parse_lit!(B * gg(euc(4, hedge_3)), default_namespace = "spenso");
    let expr3 = parse_lit!(
        C * ggg(euc(4, hedge_5)) * g(euc(4, hedge_4)),
        default_namespace = "spenso"
    );
    let expr4 = parse_lit!(
        A * B * ggg(euc(4, hedge_5)) * g(euc(4, hedge_4)),
        default_namespace = "spenso"
    );
    let fnlib = ErroringLibrary::<Symbol>::new();
    let lib = DummyLibrary::<SymbolicTensor>::new();

    for ex in [expr, expr2, expr3, expr4] {
        println!("Expr:{ex}");
        let mut acc = ex
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        // acc *= net3;
        println!("{}", acc.snapshot_dot());

        acc.execute::<Steps<1>, ContractScalars, _, _, _>(&lib, &fnlib)
            .unwrap();

        println!("{}", acc.snapshot_dot());
        acc.execute::<Sequential, SmallestDegree, _, _, _>(&lib, &fnlib)
            .unwrap();

        println!("{}", acc.snapshot_dot());

        if let ExecutionResult::Val(TensorOrScalarOrKey::Tensor { tensor, .. }) =
            acc.result().unwrap()
        {
            // println!("YaY:{}", (&expr - &tensor.expression).expand());
            assert_eq!(ex, tensor.expression);
        } else {
            panic!("Not tensor")
        }
    }
}

#[test]
fn dot_derivative() {
    initialize();

    let l = symbol!("lambda";Scalar);
    let expr = parse_lit!(
        // (2 * dot(Q3(0), Q3(2), mink(4)) + dot(Q3(0), Q3(0), mink(4)))
        dot(mink(4), Q(1), Q(3)) * dot(mink(4), Q(3), Q(1) + Q(2)),
        default_namespace = "spenso"
    );

    // 2*dot(Q3(0),Q3(2),mink(4)) is Scalar vs dot(Q3(0),Q3(0),mink(4)

    let e = expr
        .replace(parse_lit!(Q(1), default_namespace = "spenso"))
        .with(parse_lit!(Q(1), default_namespace = "spenso") * l)
        .expand();
    println!("{}", e);

    let net = e
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
            parse_inner_products: true,
            ..Default::default()
        })
        .unwrap();

    assert_snapshot!(net.snapshot_dot(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "S:(dot(Q(1),Q(3),mink(4)))^2*lambda^2+dot(Q(1),Q(3),mink(4))*dot(Q(2),Q(3),mink(4))*lambda"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
    }
    "#);

    assert_snapshot!(net.simple_execute::<()>().to_bare_ordered_string(),@"dot(mink(4),Q(1)*lambda,Q(3))*dot(mink(4),Q(3),Q(1)*lambda+Q(2))");
}

#[test]
fn symbolic_structure_parsing() {
    initialize();

    let a = parse_lit!(
        (g(mink(4, mu1), mink(4, mu8)) * g(mink(4), k(0), k(10))
            - 2 * g(mink(4, mu1), mink(4, mu8)) * g(mink(4), k(1), k(10))
            + k(0, mink(4, mu1)) * k(10, mink(4, mu8))
            - 2 * k(0, mink(4, mu8)) * k(10, mink(4, mu1))
            + k(1, mink(4, mu1)) * k(10, mink(4, mu8))
            + k(1, mink(4, mu8)) * k(10, mink(4, mu1)))
            * (-g(mink(4, mu4), mink(4, mu5)) * g(mink(4), k(0), k(20))
                + 2 * g(mink(4, mu4), mink(4, mu5)) * g(mink(4), k(4), k(20))
                + 2 * k(0, mink(4, mu4)) * k(20, mink(4, mu5))
                - k(0, mink(4, mu5)) * k(20, mink(4, mu4))
                - k(4, mink(4, mu4)) * k(20, mink(4, mu5))
                - k(4, mink(4, mu5)) * k(20, mink(4, mu4))),
        default_namespace = "spenso"
    );

    let _a = SymbolicTensor::<AbstractIndex>::parse(a.as_view()).unwrap();
}
