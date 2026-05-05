use spenso::{
    network::tags::SPENSO_TAG as T,
    structure::representation::{LibraryRep, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::Match,
};

use crate::W_;

pub trait Chain {
    /// combine adjacent chain expressions into a single chain expression.
    ///
    /// `chain(rep(d,i), rep(d,j), ...,A(..,in,out...))*chain(rep(d,j), rep(d,k), B(..,in,out...),...)`
    ///
    ///  becomes
    ///
    /// `chain(rep(d,i), rep(d,k), ...,A(..,in,out...),B(..,in,out...),...)`
    fn collect_chains(&self, representation: LibraryRep) -> Atom;

    /// Turns traced out chains into traces e.g.
    ///
    /// `chain(rep(d,i), rep(d,i), ...)` becomes `trace(rep(d),...)`
    fn normalize_chains(&self) -> Atom;

    /// turns tensors with two indices of the representation into chain expressions, for collecting using [`Atom::collect_chains`]
    ///
    /// `A(..,rep(d,i),rep(d,j),...)` becomes `chain(rep(d,i),rep(d,j),A(..,in,out...))`
    ///
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
    use insta::assert_snapshot;
    use spenso::{chain, shadowing::symbolica_utils::AtomCoreExt, slot};
    use symbolica::parse_lit;

    use crate::gamma;
    use crate::representations::{Bispinor, initialize};
    use crate::test_support::TestReps;

    use super::*;

    #[test]
    fn collect_gamma_chains_and_close_trace() {
        initialize();
        let gammas = parse_lit!(
            gamma(bis(4, 3), bis(4, 4), p(2, mink(4)))
                * gamma(bis(4, 4), bis(4, 5), mink(4, mu))
                * gamma(bis(4, 5), bis(4, 3), p(3, mink(4))),
            default_namespace = "spenso"
        );
        let rep = Bispinor {}.into();
        let normalized = gammas.chainify(rep).chainify(rep);
        let collected = normalized.collect_chains(rep);

        assert_snapshot!(collected.to_bare_ordered_string(), @"trace(bis(4),gamma(in,out,p(2,mink(4))),gamma(in,out,mink(4,mu)),gamma(in,out,p(3,mink(4))))");
    }

    #[test]
    fn collect_two_open_chains() {
        let r = TestReps::new();
        let chains = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
            gamma!(slot!(r.mink4, nu)),
        ) * chain!(
            slot!(r.bis4, b),
            slot!(r.bis4, c),
            gamma!(parse_lit!(p(1, mink(4)), default_namespace = "spenso")),
        );
        let rep = Bispinor {}.into();

        assert_snapshot!(chains.collect_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,c),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)),gamma(in,out,p(1,mink(4))))");
    }

    #[test]
    #[ignore = "pending reverse-orientation chain collection"]
    fn collect_oppositely_oriented_chain() {
        let r = TestReps::new();
        let chains = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        ) * chain!(
            slot!(r.bis4, c),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, nu)),
            gamma!(parse_lit!(p(1, mink(4)), default_namespace = "spenso")),
        );
        let rep = Bispinor {}.into();

        assert_snapshot!(chains.collect_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,c),gamma(in,out,mink(4,mu)),gamma(out,in,p(1,mink(4))),gamma(out,in,mink(4,nu)))");
    }
}
