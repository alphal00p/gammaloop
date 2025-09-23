use std::ops::Deref;

use idenso::color::SelectiveExpand;
use itertools::Itertools;
use spenso::{
    network::{
        library::{DummyKey, DummyLibrary},
        store::{NetworkStore, TensorScalarStore},
        Network, StructureLessDisplay,
    },
    structure::{
        representation::{Minkowski, RepName},
        slot::{DualSlotTo, IsAbstractSlot},
        HasName, TensorStructure,
    },
    tensors::symbolic::SymbolicTensor,
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::float::Complex as SymComplex,
    domains::rational::Rational,
    function,
    id::Replacement,
};
use tracing::{debug, event_enabled, Level};

use crate::{
    numerator::aind::Aind,
    utils::{TENSORLIB, W_},
};

use super::ParsingNet;
pub type ParsingNetError = spenso::network::TensorNetworkError<
    spenso::structure::IndexlessNamedStructure<
        Symbol,
        Vec<Atom>,
        spenso::structure::representation::LibraryRep,
        super::aind::Aind,
    >,
>;

pub trait AtomCoreExt {
    fn wrap_color(&self, symbol: Symbol) -> Atom;

    fn floatify(&self, prec: u32) -> Atom;

    fn canonize_spenso(&self) -> Atom;

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom;

    fn unwrap_function(&self, symbol: Symbol) -> Atom;

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError>;
}

impl AtomCoreExt for Atom {
    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
    }
    fn floatify(&self, prec: u32) -> Atom {
        self.as_view().floatify(prec)
    }

    fn wrap_color(&self, symbol: Symbol) -> Atom {
        self.as_view().wrap_color(symbol)
    }

    fn canonize_spenso(&self) -> Atom {
        self.as_view().canonize_spenso()
    }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.as_view().unwrap_function(symbol)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        self.as_view().parse_into_net()
    }
}

impl AtomCoreExt for AtomView<'_> {
    fn floatify(&self, prec: u32) -> Atom {
        self.map_coefficient(|c| match c {
            CoefficientView::Natural(r, d, ri, di) => Coefficient::Float(SymComplex::new(
                Rational::from((r, d)).to_multi_prec_float(prec),
                Rational::from((ri, di)).to_multi_prec_float(prec),
            )),
            CoefficientView::Large(r, ri) => Coefficient::Float(SymComplex::new(
                r.to_rat().to_multi_prec_float(prec),
                ri.to_rat().to_multi_prec_float(prec),
            )),
            _ => c.to_owned(),
        })
    }

    fn canonize_spenso(&self) -> Atom {
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();
        let mut net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, DummyKey, Aind>::try_from_view(
                *self, &lib,
            )
            .unwrap();

        debug!(net=%net.graph.dot_impl(
            |i| {
                let ss = &net.store.get_scalar(i);
                format!("{:#}", ss)
            },
            |_| "AAAA".to_string(),
            |t| {
                let tt = &net.store.get_tensor(t);
                format!("{}({})",
                    tt.name().map(|a|a.to_string()).unwrap_or("NA".to_string()),
                    tt.args().map(|a|a.iter().map(|a|format!("{}",a)).join(",")).unwrap_or(String::new()))
            },
        ),"Network for canonization");

        let mut redual_reps = vec![];

        let mut indices = vec![];

        for t in net.store.tensors.iter_mut() {
            let mut reps = vec![];
            debug!(name=%t.name().unwrap());
            let mut pat = FunctionBuilder::new(t.name().unwrap());
            let mut rhs = pat.clone();
            for s in t.structure.external_structure_iter() {
                debug!(slot=%s.to_string(),rep=%s.rep().to_string());
                if !s.rep_name().is_self_dual() && s.rep_name().is_dual() {
                    pat = pat.add_arg(s.rep().dual().to_symbolic([Atom::var(W_.a_)]));
                    indices.push((s.dual().to_atom(), s.dual().rep()));
                    reps.push(Replacement::new(
                        s.rep().to_symbolic([Atom::var(W_.a_)]).to_pattern(),
                        s.rep().dual().to_symbolic([Atom::var(W_.a_)]),
                    ));
                } else {
                    pat = pat.add_arg(s.rep().to_symbolic([Atom::var(W_.a_)]));
                    indices.push((s.to_atom(), s.rep()));
                }
                rhs = rhs.add_arg(s.rep().to_symbolic([Atom::var(W_.a_)]));
            }
            if !reps.is_empty() {
                redual_reps.push(Replacement::new(pat.finish().to_pattern(), rhs.finish()));
                t.expression = t.expression.replace_multiple(&reps);
            }
        }

        debug!(len=%indices.len());

        if event_enabled!(Level::DEBUG) {
            let indices = indices
                .iter()
                .map(|(i, r)| format!("{}:{}", i.to_canonical_string(), r.to_string()))
                .join("\n");
            debug!(indices=%indices);
        }

        self.canonize_tensors(&indices)
            .unwrap()
            .replace_multiple(&redual_reps)
    }

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.replace(Minkowski {}.to_symbolic([W_.d_, W_.a_]))
            .with(Minkowski {}.to_symbolic([dim.into().into_owned(), Atom::var(W_.a_)]))
    }
    fn wrap_color(&self, symbol: Symbol) -> Atom {
        self.expand_color()
            .into_iter()
            .fold(Atom::Zero, |a, (c, s)| a + function!(symbol, c) * s)
    }

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.replace(function!(symbol, W_.a___)).with(W_.a___)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        ParsingNet::try_from_view(*self, TENSORLIB.read().unwrap().deref())
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
    use libc::YESEXPR;
    use symbolica::{atom::AtomCore, parse_lit};

    use crate::{
        dot,
        graph::{
            parse::{self, IntoGraph},
            FeynmanGraph, Graph,
        },
        initialisation::test_initialise,
        uv::UltravioletGraph,
    };

    use super::AtomCoreExt;

    #[test]
    fn canonize_color() {
        test_initialise().unwrap();
        let gls: Vec<Graph> = dot!(
            digraph{
            num = "1";

            ext0 [style=invis];
            2:0-> ext0 [id=0 dir=none is_cut=0 is_dummy=false particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 dir=none is_cut=0 is_dummy=false particle="a"];
            0:2-> 1:3 [id=2  is_dummy=false particle="d"];
            0:4-> 1:5 [id=3 dir=none  is_dummy=false particle="g"];
            3:6-> 0:7 [id=4  is_dummy=false particle="d"];
            1:8-> 2:9 [id=5  is_dummy=false particle="d"];
            2:10-> 3:11 [id=6  is_dummy=false particle="d"];
        }

        digraph GL8{
            num = 1;
        0[int_id=V_74];
        1[int_id=V_74];
        2[int_id=V_71];
        3[int_id=V_71];
        ext0 [style=invis];
        2:0-> ext0 [id=0 dir=none is_cut=0 is_dummy=false particle=a];
        ext1 [style=invis];
        ext1-> 3:1 [id=1 dir=none is_cut=0 is_dummy=false particle=a];
        0:2-> 1:3 [id=2  is_dummy=false particle=d];
        0:4-> 1:5 [id=3 dir=none  is_dummy=false particle=g];
        0:6-> 3:7 [id=4 dir=back  is_dummy=false particle="d~"];
        1:8-> 2:9 [id=5  is_dummy=false particle=d];
        2:10-> 3:11 [id=6  is_dummy=false particle=d];
        }

        )
        .unwrap();

        for g in gls {
            let mut numerator = g.numerator(&g.no_dummy());

            // TODO Check if we include overall factor in main
            numerator.state.expr *= &g.global_prefactor.num * &g.global_prefactor.projector; // * &gl5.overall_factor;
                                                                                             // numerator.state.expr = numerator.state.expr.replace_multiple(&cpl_reps);

            let numerator_color_simplified = numerator
                .clone()
                .color_simplify()
                .get_single_atom()
                .unwrap()
                .canonize_spenso();

            println!("numerator_color_simplified:{numerator_color_simplified}");
            println!("numerator:{}", numerator.state.expr);
        }
    }

    #[test]
    fn canonizations() {
        test_initialise().unwrap();
        let a = parse_lit!(
            1 / 27 * UFO::ee
                ^ 4 * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                    * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                        + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                    * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                        + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1)))
                        - gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1))))
                    * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::hedge(0))
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
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::ϵ(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ϵbar(0, spenso::mink(4, gammalooprs::hedge(0)))
        );
        println!("a:{}", a);
        let b = parse_lit!(
            1 / 27 * UFO::ee
                ^ 4 * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(1, 1)))
                    + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(1, 1))))
                    * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(4, 1)))
                        + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(4, 1))))
                    * (gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(5, 1)))
                        + gammalooprs::K(1, spenso::mink(4, gammalooprs::edge(5, 1)))
                        - gammalooprs::P(0, spenso::mink(4, gammalooprs::edge(5, 1))))
                    * spenso::g(
                        spenso::mink(4, gammalooprs::hedge(4)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(2)),
                        spenso::bis(4, gammalooprs::hedge(6)),
                        spenso::mink(4, gammalooprs::hedge(4))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(7)),
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::mink(4, gammalooprs::hedge(1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::bis(4, gammalooprs::hedge(3)),
                        spenso::mink(4, gammalooprs::hedge(5))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::mink(4, gammalooprs::hedge(0))
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
                        spenso::bis(4, gammalooprs::hedge(9)),
                        spenso::bis(4, gammalooprs::hedge(8)),
                        spenso::mink(4, gammalooprs::edge(1, 1))
                    )
                    * spenso::gamma(
                        spenso::bis(4, gammalooprs::hedge(11)),
                        spenso::bis(4, gammalooprs::hedge(10)),
                        spenso::mink(4, gammalooprs::edge(5, 1))
                    )
                    * gammalooprs::K(0, spenso::mink(4, gammalooprs::edge(2, 1)))
                    * gammalooprs::ϵ(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ϵbar(0, spenso::mink(4, gammalooprs::hedge(0)))
        );

        println!("b:{}", b);

        println!("ratio:{}", &a / &b);

        let ac = a.canonize_spenso();
        let bc = b.canonize_spenso();
        println!("ac:{}", ac);
        println!("bc:{}", bc);
        println!("ratio canonized:{}", ac / bc);
    }

    #[test]
    fn test_can() {
        let a = parse_lit!(T(a, b, c) * T(c, d, e) * T(d, b, f)(K(e) + P(e)) * (K(f) + P(f)));
        let b = parse_lit!(T(a, d) * T(d, c));

        let indices = vec![
            (parse_lit!(a), 1),
            (parse_lit!(a), 1),
            (parse_lit!(a), 1),
            (parse_lit!(b), 1),
            (parse_lit!(b), 1),
            (parse_lit!(c), 1),
            (parse_lit!(d), 1),
            (parse_lit!(e), 1),
            (parse_lit!(f), 1),
        ];

        let ac = a.canonize_tensors(&indices);
        println!("{}", ac.unwrap());
        let bc = b.canonize_tensors(&indices);
        println!("{}", bc.unwrap());

        let a = parse_lit!(
            1 / 27 * UFO::ee
                ^ 4 * (K(0, m4e_1_1) + K(1, m4e_1_1))
                    * (K(0, m4e_4_1) + K(1, m4e_4_1))
                    * (K(0, m4e_5_1) + K(1, m4e_5_1) - P(0, m4e_5_1))
                    * g(m4he_4, m4he_5)
                    * gamma(bis_4(3), bis_4(8), m4he_5)
                    * gamma(bis_4(6), bis_4(2), m4he_4)
                    * gamma(bis_4(9), bis_4(10), m4he_0)
                    * gamma(bis_4(11), bis_4(7), m4he_1)
                    * gamma(bis_4(2), bis_4(3), m4e_2_1)
                    * gamma(bis_4(7), bis_4(6), m4e_4_1)
                    * gamma(bis_4(8), bis_4(9), m4e_1_1)
                    * gamma(bis_4(10), bis_4(11), m4e_5_1)
                    * K(0, m4e_2_1)
                    * ϵ(0, m4he_1)
                    * ϵbar(0, m4he_0)
        );

        let indices = vec![
            (parse_lit!(a), 1),
            (parse_lit!(b), 1),
            (parse_lit!(c), 1),
            (parse_lit!(d), 1),
            (parse_lit!(e), 1),
            (parse_lit!(f), 1),
        ];
    }
}
