use std::sync::LazyLock;

use crate::utils::vakint_macros::vk_parse;
use ahash::HashMap;
use symbolica::{atom::Atom, id::Condition};

#[rustfmt::skip]
pub static DIRECT_SUBSTITUTIONS: LazyLock<HashMap<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>> = LazyLock::new(|| {
    HashMap::<Atom, (Atom, symbolica::id::Condition<symbolica::id::PatternRestriction>)>::from_iter([
    // Zeta[3] value
    ( vk_parse!("z3").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0000.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // Sqrt[3] value
    ( vk_parse!("sqrt3").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0001.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // cl2 is obtained in Mathematica with ResourceFunction["ClausenCl"][2, N[\[Pi]/3, 20000]]/2
    ( vk_parse!("cl2").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0002.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // // Used to be rationalized as 28375860251062315/42167665933747557. Need expression for more precision!
    // ( vk_parse!("MASTER2m22EP").unwrap(), ( vk_parse!("0.67292935529430366752247093678895`32").unwrap(), Condition::default() ) ),
    // Used to be rationalized as 28375860251062315/42167665933747557.
    // Extracted from matad using this input
    // let integral = symbolica::try_parse!(
    //     "(
    //          vakint::a_var * vakint::a_fun(vakint::x) * vakint::dot(vakint::p(1),vakint::k(0))*vakint::k(0,mink(4,1))*vakint::p(3,mink(4,1))
    //          *vakint::k(0,mink(4,77))*vakint::k(0,mink(4,88))
    //          +vakint::p(1,mink(4,77))*vakint::p(4,mink(4,88))
    //     )*vakint::topo(\
    //         vakint::prop(1,vakint::edge(0,0),vakint::k(0),MUVsq,1)\
    // )"
    // )
    // .unwrap();
    // With this custom normalization
    // integral_normalization_factor: vakint::LoopNormalizationFactor::Custom(
    //     "(\
    //         𝑖*(𝜋^((4-2*eps)/2))\
    //      * (exp(-EulerGamma))^(eps)\
    //      * (exp(-log_mu_sq))^(eps)\
    //      )^(-n_loops)"
    //         .to_string(),
    // ),
    //  and reading the O(eps) coefficient
    ( vk_parse!("MASTER2m22EP").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0003.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // Used to be rationalized as -23459561604811217099/24730222849782180
    // In Mathematica:
    // N[-360 HR[1, -1, 1, 0, om1] + 900 HR[1, -1, 1, 1, om1] +
    //    740 HR[1, 0, 1, 0, om1] - 720 HR[1, 1, -1, 0, om1] +
    //    1800 HR[1, 1, -1, 1, om1] + 2700 HR[1, 1, 1, -1, om1] -
    //    33587 HR[1, 1, 1, 1, om1] /. nhlpre, 20000]
    ( vk_parse!("MASTER3m111111FINPIECE").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0004.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // Used to be rationalized as 693238079765918839/11778746848918154
    // Obtained with:
    // N[-1*(-(229/3) + 6 Sqrt[3] Ls3 - (17 \[Pi]^2)/12 + \[Pi]^3/Sqrt[3] -
    //      6 Sqrt[3] Ls2 (-5 + Log[3]) +
    //      7 Zeta[3] /. {Ls2 -> N[Ls[2, 0, 2 Pi/3], 30000],
    //      Ls3 -> N[Ls[3, 0, 2 Pi/3], 30000]}), 20000]
    // Some convergence issue, maybe not all 20k digits are stable below
    ( vk_parse!("MASTER3m011111FIN").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0005.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // Used to be rationalized as 531237557421164097/8442111747388969
    // Obtained with:
    // N[1/240 (-5670 + 1050 \[Pi]^2 + 19 \[Pi]^4 + 13440 HR[1, 1, 0, om1] -
    //      1840 Zeta[3]) /. nhlpre, 20000]
    ( vk_parse!("MASTER3m011101EP").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0006.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
    // Used to be rationalized as 650828837758949887/14793500421610002
    // Obtained with:
    // N[1/2 Sqrt[
    //    3] (26 Sqrt[3] + Sqrt[3] \[Pi]^2 + 16 HI[1, 0, om1] -
    //      8 HI[-1, 1, 0, om1] - 8 HI[1, -1, 0, om1] + 26 HI[1, 1, 1, om1] -
    //       Sqrt[3] Zeta[3]) /. nhlpim, 2000]
    ( vk_parse!("MASTER3m001111FIN").unwrap(), ( vk_parse!(crate::data_assets::ALPHALOOP_NUMERICS_0007.text().unwrap().as_ref()).unwrap(), Condition::default() ) ),
        ])
});
