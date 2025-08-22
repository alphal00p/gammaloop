use std::ops::Deref;

use idenso::color::SelectiveExpand;
use spenso::structure::representation::{Minkowski, RepName};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, Symbol},
    function,
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

    fn unwrap_function(&self, symbol: Symbol) -> Atom {
        self.as_view().unwrap_function(symbol)
    }

    fn parse_into_net(&self) -> Result<ParsingNet, ParsingNetError> {
        self.as_view().parse_into_net()
    }
}

impl AtomCoreExt for AtomView<'_> {
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
}
