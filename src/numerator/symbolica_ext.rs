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

    fn canonize_spenso(&self) -> Atom;

    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom;

    fn unwrap_function(&self, symbol: Symbol) -> Atom;

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError>;
}

impl AtomCoreExt for Atom {
    fn map_mink_dim<'a>(&self, dim: impl Into<AtomOrView<'a>>) -> Atom {
        self.as_view().map_mink_dim(dim)
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
