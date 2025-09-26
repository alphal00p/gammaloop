use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Deref,
};

use bitvec::vec::BitVec;
use idenso::color::SelectiveExpand;
use itertools::Itertools;
use linnet::half_edge::subgraph::ModifySubgraph;
use spenso::{
    network::{
        graph::{NetworkEdge, NetworkLeaf, NetworkNode, NetworkOp},
        library::{DummyKey, DummyLibrary},
        store::{NetworkStore, TensorScalarStore},
        Network,
    },
    structure::{
        representation::{LibraryRep, Minkowski, RepName},
        slot::IsAbstractSlot,
        HasName, TensorStructure,
    },
    tensors::symbolic::SymbolicTensor,
};

use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    coefficient::{Coefficient, CoefficientView},
    domains::{float::Complex as SymComplex, rational::Rational},
    function,
    id::Replacement,
};
use tracing::debug;

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

    fn canonize_tensors_nice<'a>(
        &'a self,
        is_index: &[Replacement],
        new_dummy: impl FnMut(usize, &Atom) -> Atom,
    ) -> Atom;

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

    fn canonize_tensors_nice<'a>(
        &'a self,
        is_index: &[Replacement],
        new_dummy: impl FnMut(usize, &Atom) -> Atom,
    ) -> Atom {
        self.as_view().canonize_tensors_nice(is_index, new_dummy)
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
    fn canonize_tensors_nice<'a>(
        &'a self,
        is_index: &[Replacement],
        mut new_dummy: impl FnMut(usize, &Atom) -> Atom,
    ) -> Atom {
        let mut indices_sorted = BTreeSet::new();

        for i in 0..(is_index.len()) {
            let rep = &is_index[i..(i + 1)];
            let r = &is_index[i];
            for m in self.pattern_match(&r.pat, None, None) {
                let a = r.pat.replace_wildcards(&m);
                let group = a.replace_multiple(&rep);
                indices_sorted.insert((group, a));
            }
        }

        let mut indices = vec![];

        let mut slot_location_replacements = vec![];

        let mut last_group = Atom::Zero;

        let mut groups = vec![];

        for (slot_loc, (g, i)) in indices_sorted.into_iter().enumerate() {
            if last_group != g {
                groups.push((slot_loc, g.clone()));
                last_group = g.clone();
            };
            let last_shift = groups.last().unwrap().0;
            slot_location_replacements.push(Replacement::new(
                i.to_pattern(),
                new_dummy(slot_loc - last_shift, &g),
            ));
            indices.push((i, g));
        }

        let can = self.canonize_tensors(&indices).unwrap();

        can //.replace_multiple(&slot_location_replacements)
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

    fn canonize_spenso(&self) -> Atom {
        fn new_dummy(i: usize, b: &Atom) -> Atom {
            let mut index = Atom::Zero;
            if let AtomView::Fun(f) = b.as_view() {
                if f.get_nargs() == 1 {
                    if let AtomView::Num(n) = f.iter().next().unwrap() {
                        let dim = i64::try_from(f.iter().next().unwrap()).unwrap() as usize;
                        index = LibraryRep::try_from_symbol_coerced(f.get_symbol())
                            .unwrap()
                            .new_slot::<Aind, _, Aind>(dim, Aind::Dummy(i))
                            .to_atom();
                    }
                }
            }

            index
        }

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

        for t in net.store.tensors.iter_mut() {
            let mut reps = vec![];
            debug!(name=%t.name().unwrap());
            let mut pat = FunctionBuilder::new(t.name().unwrap());
            let mut rhs = pat.clone();
            for s in t.structure.external_structure_iter() {
                debug!(slot=%s.to_string(),rep=%s.rep().to_string());
                if !s.rep_name().is_self_dual() && s.rep_name().is_dual() {
                    pat = pat.add_arg(s.rep().dual().to_symbolic([Atom::var(W_.a_)]));

                    reps.push(Replacement::new(
                        s.rep().to_symbolic([Atom::var(W_.a_)]).to_pattern(),
                        s.rep().dual().to_symbolic([Atom::var(W_.a_)]),
                    ));
                } else {
                    pat = pat.add_arg(s.rep().to_symbolic([Atom::var(W_.a_)]));
                }
                rhs = rhs.add_arg(s.rep().to_symbolic([Atom::var(W_.a_)]));
            }
            if !reps.is_empty() {
                redual_reps.push(Replacement::new(pat.finish().to_pattern(), rhs.finish()));
                t.expression = t.expression.replace_multiple(&reps);
            }
        }

        let mut dummies = BTreeMap::new();

        let mut only_products_and_associated: BitVec = net.graph.graph.empty_subgraph();

        for (_, crown, d) in net.graph.graph.iter_nodes() {
            if let NetworkNode::Op(NetworkOp::Product) = d {
                for h in crown {
                    let n = net.graph.graph.node_id(h);
                    if let Some(cc) = net.graph.graph.involved_node_crown(h) {
                        for h in cc {
                            only_products_and_associated.add(h);
                        }
                    }
                    for h in net.graph.graph.iter_crown(n) {
                        only_products_and_associated.add(h);
                    }
                    only_products_and_associated.add(h);
                }
            }
        }

        debug!(net=%net.graph.graph.dot_impl(
                &only_products_and_associated,
                "",
                &|_| None,
                &|e| {
                    if let NetworkEdge::Slot(s) = e {
                        Some(format!("label=\"{s}\""))
                    } else {
                        None
                    }
                },
                &|n| match n {
                    NetworkNode::Leaf(l) => match l {
                        NetworkLeaf::LibraryKey(l) => {
                            Some(format!("label= \"L:\""))
                        }
                        NetworkLeaf::LocalTensor(l) => {
                            let tt = &net.store.get_tensor(*l);
                            let a = format!("{}({})",
                                tt.name().map(|a|a.to_string()).unwrap_or("NA".to_string()),
                                tt.args().map(|a|a.iter().map(|a|format!("{}",a)).join(",")).unwrap_or(String::new()));
                            Some(format!("label = \"T:{a}\""))
                        }
                        NetworkLeaf::Scalar(s) => { let ss = &net.store.get_scalar(*s);
                            Some(format!("label = \"S:{ss}\""))},
                    },
                    NetworkNode::Op(o) => Some(format!("label = \"{o}\"")),
                },
            ),"only_products_and_associated"
        );

        for (p, e, d) in net.graph.graph.iter_edges_of(&only_products_and_associated) {
            if p.is_paired() {
                if let NetworkEdge::Slot(s) = d.data {
                    dummies.insert(s.to_atom(), s.rep().to_symbolic([]));
                }
            }
        }

        let index_ident_pat: Vec<Replacement> = dummies
            .iter()
            .map(|(k, v)| {
                debug!(k=%k,v=%v.to_string(),"Index identified for canonization");
                Replacement::new(k.to_pattern(), v.clone())
            })
            .collect();

        // mink(W_.dim,W_.i)=>mink(W_.dim)

        self.canonize_tensors_nice(&index_ident_pat, new_dummy)
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
    use symbolica::{atom::AtomCore, parse_lit};

    use crate::{
        dot,
        graph::{parse::IntoGraph, FeynmanGraph, Graph},
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
                    * gammalooprs::ϵ(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ϵbar(0, spenso::mink(4, gammalooprs::hedge(0)))
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
                    * gammalooprs::ϵ(0, spenso::mink(4, gammalooprs::hedge(1)))
                    * gammalooprs::ϵbar(0, spenso::mink(4, gammalooprs::hedge(0)))
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
