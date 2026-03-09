use insta::assert_snapshot;
use spenso::network::parsing::ParseSettings;
use symbolica::parse_lit;

use crate::common::{HepAtomExt, NetExt, test_initialize};

mod common;
#[test]
fn exterior_prod_simple() {
    test_initialize();

    let expr = parse_lit!(
        spenso::g(
            spenso::cof(3, _gammaloop::hedge_3),
            spenso::dind(spenso::cof(3, _gammaloop::hedge_4))
        ) * spenso::g(
            spenso::cof(3, gammalooprs::hedge_1),
            spenso::dind(spenso::cof(3, gammalooprs::hedge_2))
        )
    );

    println!("{expr}");

    let net = expr.parse_to_hep_net(&ParseSettings::default()).unwrap();

    assert_snapshot!(net.dot_pretty(),@r#"
    digraph {
      node	 [shape=circle,height=0.1,label=""];
      overlap = "scale";
      layout = "neato";

      0	 [label = "∏"];
      1	 [label= "L:g()"];
      2	 [label= "L:g()"];
      ext0	 [style=invis];
      0:0:s	-> ext0	 [id=0 color="red"];
      ext1	 [style=invis];
      2:8:s	-> ext1	 [id=1 dir=back color="red" label="cof🠓3|_hedge_2"];
      ext2	 [style=invis];
      1:5:s	-> ext2	 [id=2 dir=back color="red" label="cof🠓3|_hedge_4"];
      ext3	 [style=invis];
      1:4:s	-> ext3	 [id=3 color="red" label="cof🠑3|^hedge_3"];
      1:3:s	-> 0:2:s	 [id=4  color="red:blue;0.5"];
      ext5	 [style=invis];
      2:7:s	-> ext5	 [id=5 color="red" label="cof🠑3|^hedge_1"];
      2:6:s	-> 0:1:s	 [id=6  color="red:blue;0.5"];
    }
    "#);

    let r = net.execute_and_res().unwrap();

    println!("{}", r)
}
