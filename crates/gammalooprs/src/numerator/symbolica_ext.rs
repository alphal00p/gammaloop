use std::ops::Deref;

use idenso::{
    color::{CS, ColorSimplifier},
    representations::{ColorAdjoint, ColorFundamental},
};
use spenso::{
    network::parsing::ParseSettings,
    structure::representation::{Minkowski, RepName},
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::{float::Complex as SymComplex, rational::Rational},
    function,
};

use crate::utils::{GS, TENSORLIB, W_};

use super::ParsingNet;
pub type ParsingNetError = spenso::network::TensorNetworkError<
    spenso::structure::IndexlessNamedStructure<
        Symbol,
        Vec<Atom>,
        spenso::structure::representation::LibraryRep,
        super::aind::Aind,
    >,
    Symbol,
>;

pub trait AtomCoreExt {
    fn to_param_color(&self) -> Atom;
    // fn wrap_color(&self, symbol: Symbol) -> Atom;
    fn kill_color(&self) -> Atom;

    fn floatify(&self, prec: u32) -> Atom;

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom;

    fn unwrap_function(&self, symbol: Symbol) -> Atom;

    #[allow(clippy::result_large_err)]
    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError>;
}

impl AtomCoreExt for Atom {
    fn to_param_color(&self) -> Atom {
        self.as_view().to_param_color()
    }
    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
    }
    fn floatify(&self, prec: u32) -> Atom {
        self.as_view().floatify(prec)
    }

    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.as_view().wrap_color(symbol)
    // }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.as_view().unwrap_function(symbol)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        self.as_view().parse_into_net()
    }
}

impl AtomCoreExt for AtomView<'_> {
    fn kill_color(&self) -> Atom {
        self.wrap_color(GS.killing_func)
    }

    fn to_param_color(&self) -> Atom {
        let adj = ColorAdjoint {};
        let fund = ColorFundamental {};
        self.replace(adj.to_symbolic([W_.d_, W_.a_]))
            .with(adj.to_symbolic([CS.nc * CS.nc - 1, Atom::var(W_.a_)]))
            .replace(fund.to_symbolic([W_.d_, W_.a_]))
            .with(fund.to_symbolic([CS.nc, W_.a_]))
    }
    fn floatify(&self, prec: u32) -> Atom {
        self.map_coefficient(|c| match c {
            CoefficientView::Natural(r, d, ri, di) => {
                if (ri == 0 || di == 1) && (r == 0 || d == 1) {
                    Coefficient::Complex(SymComplex::new(
                        Rational::new(r, d),
                        Rational::new(ri, di),
                    ))
                } else {
                    Coefficient::Float(SymComplex::new(
                        Rational::from((r, d)).to_multi_prec_float(prec),
                        Rational::from((ri, di)).to_multi_prec_float(prec),
                    ))
                }
            }
            CoefficientView::Large(r, ri) => Coefficient::Float(SymComplex::new(
                r.to_rat().to_multi_prec_float(prec),
                ri.to_rat().to_multi_prec_float(prec),
            )),
            _ => c.to_owned(),
        })
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.replace(Minkowski {}.to_symbolic([W_.d_, W_.a___]))
            .with(Minkowski {}.to_symbolic([dim.into().into_owned(), Atom::var(W_.a___)]))
    }
    // fn wrap_color(&self, symbol: Symbol) -> Atom {
    //     self.expand_color()
    //         .into_iter()
    //         .fold(Atom::Zero, |a, (c, s)| a + function!(symbol, c) * s)
    // }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.replace(function!(symbol, W_.a___)).with(W_.a___)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        ParsingNet::try_from_view(
            *self,
            TENSORLIB.read().unwrap().deref(),
            &ParseSettings {
                parse_inner_products: false,
                ..Default::default()
            },
        )
    }

    // fn parse_into_only_lib_net<T: TensorLibraryData + Clone + Default>(
    //     &self,
    //     one: T,
    //     zero: T,
    // ) -> Result<ParsingNet, ParsingNetError> {
    //     let mut lib = hep_lib(one, zero);

    //     ParsingNet::try_from_view(*self, &lib);
    // }
}

#[cfg(test)]
mod tests {
    use idenso::IndexTooling;
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::parsing::SPENSO_TAG,
        structure::{
            representation::{Minkowski, RepName},
            slot::{DummyAind, IsAbstractSlot, Slot},
        },
    };
    use symbolica::{
        atom::{Atom, AtomCore},
        function, parse_lit, symbol,
    };

    use crate::{
        dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        numerator::aind::Aind,
        utils::GS,
        uv::UltravioletGraph,
    };

    use super::AtomCoreExt;

    #[test]
    fn dummy_parsing() {
        test_initialise().unwrap();

        let e_mass = parse_lit!(M_e);

        let m2 = &e_mass * &e_mass;

        let mink: Slot<Minkowski, Aind> = Minkowski {}.new_rep(4).slot(Aind::new_dummy());

        let e = EdgeIndex(0);

        let sqrt = symbol!("sqrt_scalar", tag = SPENSO_TAG.tag);

        let a = function!(
            sqrt,
            (GS.emr_vec_index(e, mink.to_atom()) * GS.emr_vec_index(e, mink.to_atom()) + m2)
                .pow(Atom::num(2))
        );

        let net = a.parse_into_net().unwrap();

        println!("{}", net.dot_pretty())
    }

    #[test]
    fn canonize_color() {
        test_initialise().unwrap();
        let gls: Vec<Graph> = dot!(
            digraph{
            num = "1";

            ext0 [style=invis];
            2:0-> ext0 [id=0 dir=none is_cut=0  particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 dir=none is_cut=0  particle="a"];
            0:2-> 1:3 [id=2   particle="d"];
            0:4-> 1:5 [id=3 dir=none   particle="g"];
            3:6-> 0:7 [id=4   particle="d"];
            1:8-> 2:9 [id=5   particle="d"];
            2:10-> 3:11 [id=6   particle="d"];
        }

        digraph GL8{
            num = 1;
        0[int_id=V_74];
        1[int_id=V_74];
        2[int_id=V_71];
        3[int_id=V_71];
        ext0 [style=invis];
        2:0-> ext0 [id=0 dir=none is_cut=0  particle=a];
        ext1 [style=invis];
        ext1-> 3:1 [id=1 dir=none is_cut=0  particle=a];
        0:2-> 1:3 [id=2   particle=d];
        0:4-> 1:5 [id=3 dir=none   particle=g];
        0:6-> 3:7 [id=4 dir=back   particle="d~"];
        1:8-> 2:9 [id=5   particle=d];
        2:10-> 3:11 [id=6   particle=d];
        }

        )
        .unwrap();

        for g in gls {
            let mut numerator = g.numerator(&g.no_dummy(), &g.empty_subgraph());

            // TODO Check if we include overall factor in main
            numerator.state.expr *= &g.global_prefactor.num * &g.global_prefactor.projector; // * &gl5.overall_factor;
            // numerator.state.expr = numerator.state.expr.replace_multiple(&cpl_reps);

            let numerator_color_simplified = numerator
                .clone()
                .color_simplify()
                .get_single_atom()
                .unwrap()
                .canonize(Aind::Dummy);

            println!("numerator_color_simplified:{numerator_color_simplified}");
            println!("numerator:{}", numerator.state.expr);
        }
    }

    #[test]
    fn canonizations() {
        test_initialise().unwrap();

        let a = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(2))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(2))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(8))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
        );
        println!("a:{}", a);
        let b = parse_lit!(
            ((-2 * spenso::projp(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            ) + spenso::projm(
                spenso::bis(4, gammalooprs::edge(0)),
                spenso::bis(4, gammalooprs::hedge(6))
            )) * -1𝑖
                / 6
                * UFO::sw
                ^ 2 + -1𝑖 / 2 * UFO::cw
                ^ 2 * spenso::projm(
                    spenso::bis(4, gammalooprs::edge(0)),
                    spenso::bis(4, gammalooprs::hedge(6))
                ))
                * ((-2
                    * spenso::projp(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    )
                    + spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                    * -1𝑖
                    / 6
                    * UFO::sw
                    ^ 2 + -1𝑖 / 2 * UFO::cw
                    ^ 2 * spenso::projm(
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::bis(4, gammalooprs::hedge(3))
                    ))
                * (-1 * UFO::MZ
                    ^ 2 * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    ) + gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(4)))
                        * gammalooprs::K(1, spenso::mink(4, gammalooprs::hedge(5))))
                * (-1 * gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                * 1
                / 3
                * UFO::MZ
                ^ (-2) * UFO::cw
                ^ (-2) * UFO::ee
                ^ 4 * UFO::sw
                ^ (-2)
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::e(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ebar(0, spenso::mink(4, gammalooprs::hedge(0)))
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::hedge(0))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::edge(0)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::mink(4, gammalooprs::edge(2, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::mink(4, gammalooprs::edge(4, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(0)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
        );

        println!("b:{}", b);

        println!("ratio:{}", &a / &b);

        let ac = a.canonize(Aind::Dummy);
        let bc = b.canonize(Aind::Dummy);
        println!("ac:{}", ac);
        println!("bc:{}", bc);
        println!("ratio canonized:{}", ac / bc);
    }

    // #[test]
    // fn test_can() {
    //     let a = parse_lit!(T(a, b, c) * T(c, d, e) * T(d, b, f)(K(e) + P(e)) * (K(f) + P(f)));
    //     let b = parse_lit!(T(a, d) * T(d, c));

    //     let indices = vec![
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(a), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(b), 1),
    //         (parse_lit!(c), 1),
    //         (parse_lit!(d), 1),
    //         (parse_lit!(e), 1),
    //         (parse_lit!(f), 1),
    //     ];

    //     let ac = a.canonize_tensors(&indices);
    //     println!("{}", ac.unwrap());
    //     let bc = b.canonize_tensors(&indices);
    //     println!("{}", bc.unwrap());
    // }
}
