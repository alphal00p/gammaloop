use spenso::{
    network::tags::SPENSO_TAG as T,
    structure::representation::{LibraryRep, RepName, Representation},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Match,
};

use crate::W_;

pub trait Chain {
    fn collect_chains(&self, representation: LibraryRep) -> Atom;

    fn normalize_chains(&self) -> Atom;

    fn chainify(&self, representation: LibraryRep) -> Atom;
}

impl Chain for Atom {
    fn collect_chains(&self, representation: LibraryRep) -> Atom {
        self.as_view().collect_chains(representation)
    }

    fn chainify(&self, representation: LibraryRep) -> Atom {
        self.as_view().chainify(representation)
    }

    fn normalize_chains(&self) -> Atom {
        self.as_view().normalize_chains()
    }
}
impl<'a> Chain for AtomView<'a> {
    fn collect_chains(&self, representation: LibraryRep) -> Atom {
        let in_index = representation.to_symbolic([W_.d_, W_.i_]);
        let dummy_out = representation.dual().to_symbolic([W_.d_, W_.j_]);
        let dummy_in = representation.to_symbolic([W_.d_, W_.j_]);
        let out_index = representation.to_symbolic([W_.d_, W_.k_]);

        let product = function!(
            T.chain,
            in_index,
            dummy_out,
            W_.d___,
            function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___)
        ) * function!(
            T.chain,
            dummy_in,
            out_index,
            function!(W_.b_, W_.e___, T.chain_in, W_.f___, T.chain_out, W_.g___),
            W_.h___
        );

        // println!("{}", product);
        let collected = function!(
            T.chain,
            in_index,
            out_index,
            W_.d___,
            function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___),
            function!(W_.b_, W_.e___, T.chain_in, W_.f___, T.chain_out, W_.g___),
            W_.h___
        );
        // println!("{}", collected);

        self.replace(product)
            .repeat()
            .with(collected)
            .normalize_chains()
    }

    fn chainify(&self, representation: LibraryRep) -> Atom {
        let in_index = representation.to_symbolic([W_.d_, W_.i_]);

        let out_index = representation.dual().to_symbolic([W_.d_, W_.j_]);

        self.replace(function!(
            W_.a_, W_.a___, in_index, W_.b___, out_index, W_.c___
        ))
        .when(W_.a_.filter_match(|a| {
            if let Match::FunctionName(a) = a
                && *a != T.chain
            {
                true
            } else {
                false
            }
        }))
        .with(function!(
            T.chain,
            in_index,
            out_index,
            function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___)
        ))
    }

    fn normalize_chains(&self) -> Atom {
        let rep = T.self_dual_::<1, _>([W_.d_, W_.i_]);
        let stripped_rep = T.self_dual_::<1, _>([W_.d_]);

        let traced_chain = function!(T.chain, rep, rep, W_.a___);
        self.replace(traced_chain)
            .with(function!(T.trace, stripped_rep, W_.a___))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::parse_lit;

    use crate::representations::{Bispinor, initialize};

    use super::*;

    #[test]
    fn collect_gamma_chains() {
        initialize();
        let gammas = parse_lit!(
            gamma(bis(4, 3), bis(4, 4), p(2, mink(4)))
                * gamma(bis(4, 4), bis(4, 5), mink(4, mu))
                * gamma(bis(4, 5), bis(4, 3), p(3, mink(4))),
            default_namespace = "spenso"
        );
        let rep = Bispinor {}.into();
        let normalized = gammas.chainify(rep).chainify(rep);
        println!("{}", normalized);
        let collected = normalized.collect_chains(rep);

        println!("{}", collected);
    }
}
