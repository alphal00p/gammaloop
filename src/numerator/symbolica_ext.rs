use std::ops::Deref;

use idenso::color::SelectiveExpand;
use spenso::{
    network::{library::DummyLibrary, store::NetworkStore, Network},
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

use crate::utils::{TENSORLIB, W_};

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
        let lib = DummyLibrary::<SymbolicTensor>::new();
        let mut net =
            Network::<NetworkStore<SymbolicTensor, Atom>, _>::try_from_view(*self, &lib).unwrap();

        let mut redual_reps = vec![];

        let mut indices = vec![];

        for t in net.store.tensors.iter_mut() {
            let mut reps = vec![];
            let mut pat = FunctionBuilder::new(t.name().unwrap());
            let mut rhs = pat.clone();
            for s in t.structure.external_structure_iter() {
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
}
