use std::sync::LazyLock;

use spenso::{
    chain, dualizable_, dualizable_dual_,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    self_dual_,
    shadowing::Collectable,
    structure::representation::{LibraryRep, RepName},
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
    id::{Match, Replacement},
};
use symbolica_utils::{PatternReplacement, ReplaceBuilderExt};

use crate::W_;
static SINGLE_LENGTH_NORM: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let sd_rep = self_dual_!(1; W_.d_, W_.i_);

    let sd_rep2 = self_dual_!(2; W_.d_, W_.j_);
    let d_rep = dualizable_!(1; W_.d_, W_.i_);
    let dd_rep = dualizable_dual_!(2; W_.d_, W_.j_);

    [
        Replacement::new(
            chain!(
                &sd_rep,
                &sd_rep2,
                function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___)
            )
            .to_pattern(),
            function!(W_.a_, W_.a___, &sd_rep, W_.b___, &sd_rep2, W_.c___),
        ),
        Replacement::new(
            chain!(
                &d_rep,
                &dd_rep,
                function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___)
            )
            .to_pattern(),
            function!(W_.a_, W_.a___, &d_rep, W_.b___, &dd_rep, W_.c___),
        ),
    ]
});

static CHAIN_NORMALIZATIONS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let sd_rep = self_dual_!(1; W_.d_, W_.i_);

    let sd_stripped_rep = self_dual_!(1; W_.d_);
    let d_rep = dualizable_!(1; W_.d_, W_.i_);
    let dd_rep = dualizable_dual_!(1; W_.d_, W_.i_);
    let d_stripped_rep = dualizable_!(1; W_.d_);

    [
        Replacement::new(
            // A chain whose endpoints are the same slot is a closed line:
            // chain(rep(d,i), rep(d,i), factors...) -> trace(rep(d), factors...).
            chain!(&sd_rep, &sd_rep, W_.a___).to_pattern(),
            trace!(sd_stripped_rep, W_.a___),
        ),
        Replacement::new(
            chain!(&d_rep, &dd_rep, W_.a___).to_pattern(),
            trace!(d_stripped_rep, W_.a___),
        ),
    ]
});

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

    /// and single length chains back into the corresponding tensor
    fn undo_single_length(&self) -> Atom;

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

    fn undo_single_length(&self) -> Atom {
        self.as_view().undo_single_length()
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
        // Every rule below requires a chain. Avoid polynomial collection of
        // unrelated scalar factors when there is no chain to compose.
        if !self.contains_symbol(T.chain) {
            return self.to_owned();
        }

        let in_index = representation.to_symbolic([W_.d_, W_.i_]);
        let dummy_out = representation.dual().to_symbolic([W_.d_, W_.j_]);
        let dummy_in = representation.to_symbolic([W_.d_, W_.j_]);
        let out_index = representation.dual().to_symbolic([W_.d_, W_.k_]);
        let chain = function!(T.chain, &in_index, &out_index, W_.x___).to_pattern();

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

        // Use the shared collector's opaque coefficients here as well: chain
        // composition must not statistically simplify the momentum numerator.
        // Select the ports used by composition, not incidental representations
        // inside the payload or chains belonging to another representation.
        self.collect_with_map(|atom| {
            atom.get_symbol() == Some(T.chain) && atom.replace(&chain).max_level(0).matches()
        })
        .unwrap_collect()
        .replace(product)
        .repeat()
        .with(collected)
        .normalize_chains()
    }

    fn chainify(&self, representation: LibraryRep) -> Atom {
        let in_index = representation.to_symbolic([W_.d_, W_.i_]);

        let out_index = representation.dual().to_symbolic([W_.d_, W_.j_]);

        // Existing ports belong to an enclosing chain. A second representation
        // must retain its explicit slots rather than reuse those untyped ports.
        let no_ports = |matched: &Match<'_>| {
            let args = match matched {
                Match::Multiple(_, args) => args.as_slice(),
                Match::Single(arg) => std::slice::from_ref(arg),
                Match::FunctionName(_) => return false,
            };
            args.iter()
                .all(|arg| !arg.contains_symbol(T.chain_in) && !arg.contains_symbol(T.chain_out))
        };
        let replacement = Replacement::new(
            function!(W_.a_, W_.a___, in_index, W_.b___, out_index, W_.c___).to_pattern(),
            function!(
                T.chain,
                in_index,
                out_index,
                function!(W_.a_, W_.a___, T.chain_in, W_.b___, T.chain_out, W_.c___)
            ),
        )
        .when(
            W_.a_.filter_match(
                |a| matches!(a, Match::FunctionName(a) if *a != T.chain && *a != ETS.metric),
            ) & W_.a___.filter_match(no_ports)
                & W_.b___.filter_match(no_ports)
                & W_.c___.filter_match(no_ports),
        )
        .max_level(0);
        self.replace_map(|atom, _, out| {
            let AtomView::Fun(fun) = atom else { return };
            if [T.chain, T.trace].contains(&fun.get_symbol()) {
                // Chain and trace share one untyped in/out binder. In particular,
                // an explicit spin tensor inside a color word must not acquire
                // a second nested binder; its slots remain explicit for parsing.
                out.set_from_view(&atom);
            } else {
                let rewritten = atom.replace_multiple([&replacement]);
                if rewritten.as_view() != atom {
                    **out = rewritten;
                }
            }
        })
    }

    fn normalize_chains(&self) -> Atom {
        self.to_owned()
            .replace_multiple_repeat(CHAIN_NORMALIZATIONS.as_ref())
    }

    fn undo_single_length(&self) -> Atom {
        self.to_owned()
            .replace_multiple_repeat(SINGLE_LENGTH_NORM.as_ref())
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use spenso::g;
    use spenso::{chain, slot};
    use symbolica::{parse, parse_lit};
    use symbolica_utils::AtomPrintExt;

    use crate::representations::{Bispinor, ColorFundamental};
    use crate::test_support::{TestReps, test_initialize};
    use crate::{bis, gamma};

    use super::*;

    #[test]
    fn mixed_representation_chains_preserve_slots_and_round_trip() {
        use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
        use spenso::{network::parsing::ParseSettings, structure::abstract_index::AbstractIndex};

        test_initialize();
        let input = parse_lit!(
            (x + y) ^ 8 * mixed_chain_tensor(bis(4, i), bis(4, j), cof(3, a), dind(cof(3, b))),
            default_namespace = "spenso"
        );
        let settings = ParseSettings::default();
        let mut expected_slots = input
            .parse_to_symbolic_net::<AbstractIndex>(&settings)
            .unwrap()
            .graph
            .dangling_indices();
        expected_slots.sort();
        assert_eq!(expected_slots.len(), 4);
        let spin = Bispinor {}.into();
        let color = ColorFundamental {}.into();
        for (first, second) in [(spin, color), (color, spin)] {
            let chained = input.chainify(first).chainify(second);
            let net = chained
                .parse_to_symbolic_net::<AbstractIndex>(&settings)
                .unwrap();
            let mut slots = net.graph.dangling_indices();
            slots.sort();
            assert_eq!(slots, expected_slots);
            // Parsing restores both tensor interfaces and the untouched scalar
            // factor, without distributing any numerator products or powers.
            assert_eq!(net.simple_execute::<()>().unwrap(), input);
        }
    }

    #[test]
    fn chainify_keeps_existing_binder_scopes_opaque() {
        use crate::{color::CS, dirac::gamma_tensor};
        use spenso::network::parsing::AtomStructureExt;

        test_initialize();
        let left = parse_lit!(spenso::bis(4, 3));
        let right = parse_lit!(spenso::bis(4, 4));
        let mu = parse_lit!(spenso::mink(4, 7));
        let gamma = gamma_tensor(left.clone(), right.clone(), mu.clone());
        let generator = CS.chain_t(parse_lit!(spenso::coad(8, 9)));
        let color = trace!(
            parse_lit!(spenso::cof(3)),
            generator.clone(),
            generator * gamma,
        );
        let spectator = parse_lit!((scope_a + scope_b) ^ 8);
        let outer = gamma_tensor(right.clone(), left.clone(), mu.clone());
        let expected_outer = chain!(right, left, gamma!(mu));
        let spin = Bispinor {}.into();

        assert_eq!(color.chainify(spin), color);
        let result = (&spectator * &color * outer).chainify(spin);
        assert_eq!(result, spectator * color * expected_outer);
        assert!(result.validate_chain_like_nesting().is_ok());
    }

    #[test]
    fn collect_chains_preserves_factorized_scalars_without_chains() {
        test_initialize();
        let scalar = parse_lit!((a + b) ^ 8 * (c + d) ^ 8 / (1 + a * c + b * d));
        for representation in [Bispinor {}.into(), ColorFundamental {}.into()] {
            // Structural equality checks the original factorization, not only
            // equality after expanding or evaluating the scalar expression.
            assert_eq!(scalar.collect_chains(representation), scalar);
        }
    }

    #[test]
    fn collect_gamma_chains_and_close_trace() {
        test_initialize();
        let gammas = parse_lit!(
            gamma(bis(4, 3), bis(4, 4), p(2, mink(4)))
                * gamma(bis(4, 4), bis(4, 5), mink(4, mu))
                * gamma(bis(4, 5), bis(4, 3), p(3, mink(4))),
            default_namespace = "spenso"
        );
        let rep = Bispinor {}.into();
        let normalized = gammas.chainify(rep).chainify(rep);
        let collected = normalized.collect_chains(rep);

        assert_snapshot!(collected.to_bare_ordered_string(), @"trace(bis(4),cyclic(gamma(in,out,mink(4,mu)),gamma(in,out,p(3,mink(4))),gamma(in,out,p(2,mink(4)))))");
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
    fn collect_chains_preserves_opaque_factorized_spectators() {
        let r = TestReps::new();
        let first = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        );
        let second = chain!(
            slot!(r.bis4, b),
            slot!(r.bis4, c),
            gamma!(slot!(r.mink4, nu)),
        );
        let composed = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, c),
            gamma!(slot!(r.mink4, mu)),
            gamma!(slot!(r.mink4, nu)),
        );
        let spectator = parse_lit!((x + y) ^ 8 * (u + v) ^ n / (1 + x * u + y * v));
        let input = &spectator * first * second;

        // Exact Atom equality retains the coefficient's sums and powers while
        // checking the ordered chain payload independently of scalar algebra.
        assert_eq!(
            input.collect_chains(Bispinor {}.into()),
            &spectator * composed
        );
        assert_eq!(input.collect_chains(ColorFundamental {}.into()), input);
    }

    #[test]
    fn collect_chains_keeps_other_representations_factored() {
        test_initialize();
        let spin = parse!(
            "chain(bis(4, i), bis(4, j), gamma(in, out, mink(4, mu)))
                * chain(bis(4, j), bis(4, k), gamma(in, out, mink(4, nu)))",
            default_namespace = "spenso"
        );
        let spin_composed = parse!(
            "chain(
                bis(4, i), bis(4, k),
                gamma(in, out, mink(4, mu)), gamma(in, out, mink(4, nu))
            )",
            default_namespace = "spenso"
        );
        let color = parse!(
            "chain(cof(3, a), dind(cof(3, b)), t(coad(8, alpha), in, out))
                * (chain(cof(3, b), dind(cof(3, c)), t(coad(8, beta), in, out))
                    + chain(cof(3, b), dind(cof(3, c)), t(coad(8, delta), in, out)))",
            default_namespace = "spenso"
        );
        let color_composed = parse!(
            "chain(
                cof(3, a), dind(cof(3, c)),
                t(coad(8, alpha), in, out), t(coad(8, beta), in, out)
            ) + chain(
                cof(3, a), dind(cof(3, c)),
                t(coad(8, alpha), in, out), t(coad(8, delta), in, out)
            )",
            default_namespace = "spenso"
        );
        let input = &spin * &color;

        assert_eq!(color.collect_chains(Bispinor {}.into()), color);
        assert_eq!(
            input.collect_chains(Bispinor {}.into()),
            &spin_composed * &color
        );
        assert_eq!(
            input.collect_chains(ColorFundamental {}.into()),
            &spin * &color_composed
        );
    }

    #[test]
    fn collect_chains_selects_ports_instead_of_payload_representations() {
        test_initialize();
        let input = parse!(
            "chain(cof(3, a), dind(cof(3, b)), mixed(in, out, bis(4, i), bis(4, j)))
                * (chain(cof(3, b), dind(cof(3, c)), mixed(in, out, bis(4, j), bis(4, k)))
                    + chain(cof(3, b), dind(cof(3, c)), other(in, out, bis(4, j), bis(4, k))))",
            default_namespace = "spenso"
        );

        // The explicit spin slots belong to the mixed tensor payload; they do
        // not turn its fundamental-color chain ports into spin-chain ports.
        assert_eq!(input.collect_chains(Bispinor {}.into()), input);
    }

    #[test]
    fn undo_single_gamma_chain() {
        let r = TestReps::new();
        let chain = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        );

        assert_snapshot!(chain.undo_single_length().to_bare_ordered_string(), @"gamma(bis(4,a),bis(4,b),mink(4,mu))");
    }

    #[test]
    fn collect_chains_leaves_single_gamma_chain_for_explicit_undo() {
        let r = TestReps::new();
        let chain = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        );
        let rep = Bispinor {}.into();

        assert_snapshot!(chain.collect_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))");
    }

    #[test]
    fn collect_chains_composes_common_prefix_terms_before_factoring() {
        let r = TestReps::new();
        let shared_prefix = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        );
        let first_tail = chain!(
            slot!(r.bis4, b),
            slot!(r.bis4, c),
            gamma!(slot!(r.mink4, nu)),
        );
        let second_tail = chain!(
            slot!(r.bis4, b),
            slot!(r.bis4, d),
            gamma!(slot!(r.mink4, rho)),
        );
        let chains = shared_prefix.clone() * first_tail + shared_prefix * second_tail;
        let rep = Bispinor {}.into();

        assert_snapshot!(chains.collect_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,c),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)))+chain(bis(4,a),bis(4,d),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,rho)))");
    }

    #[test]
    fn chainify_does_not_hit_metrics() {
        let a = g!(bis!(4, 1), bis!(4, 2));

        assert_eq!(a.chainify(Bispinor {}.into()), a);
    }

    #[test]
    fn normalize_works_on_gammaloop_input() {
        let expr = parse!(
            "spenso::chain(spenso::bis(4,gammalooprs::hedge(17)),spenso::bis(4,gammalooprs::hedge(17)),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(2))),spenso::gamma(spenso::in,spenso::out,gammalooprs::Q(9,spenso::mink(gammalooprs::dim))),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(1))),spenso::gamma(spenso::in,spenso::out,gammalooprs::Q(9,spenso::mink(gammalooprs::dim))),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(7))),spenso::gamma(spenso::in,spenso::out,gammalooprs::Q(4,spenso::mink(gammalooprs::dim))),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(7))),spenso::gamma(spenso::in,spenso::out,gammalooprs::Q(9,spenso::mink(gammalooprs::dim))),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(0))),spenso::gamma(spenso::in,spenso::out,gammalooprs::Q(9,spenso::mink(gammalooprs::dim))),spenso::gamma(spenso::in,spenso::out,spenso::mink(gammalooprs::dim,gammalooprs::hedge(3))))"
        );

        println!("{}", expr.normalize_chains())
    }
}
