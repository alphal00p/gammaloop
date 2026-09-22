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
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
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
    ///
    /// In self-dual spaces, common starts or ends can also be joined by
    /// transposing one chain: reverse its factors and exchange `in`/`out`.
    /// This is an ordinary transpose and does not conjugate scalar entries.
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

        let mut joins = vec![(
            chain!(&in_index, &dummy_out, W_.a___) * chain!(&dummy_in, &out_index, W_.b___),
            false,
            false,
        )];
        if representation.is_self_dual() {
            // Common-end contractions need a transpose, not a conjugate:
            // A(i,j) B(k,j) = (A B^T)(i,k), and similarly for common starts.
            joins.extend([
                (
                    chain!(&in_index, &dummy_out, W_.a___) * chain!(&out_index, &dummy_in, W_.b___),
                    false,
                    true,
                ),
                (
                    chain!(&dummy_out, &in_index, W_.a___) * chain!(&dummy_in, &out_index, W_.b___),
                    true,
                    false,
                ),
            ]);
        }
        // Use the shared collector's opaque coefficients here as well: chain
        // composition must not statistically simplify the momentum numerator.
        // Select the ports used by composition, not incidental representations
        // inside the payload or chains belonging to another representation.
        let mut result = self
            .collect_with_map(|atom| {
                atom.get_symbol() == Some(T.chain) && atom.replace(&chain).max_level(0).matches()
            })
            .unwrap_collect();
        loop {
            let previous = result.clone();
            for (product, transpose_left, transpose_right) in &joins {
                // println!("{}", product);
                let start = in_index.to_pattern();
                let end = out_index.to_pattern();
                let transpose_left = *transpose_left;
                let transpose_right = *transpose_right;
                result = result
                    .replace(product.to_pattern())
                    .repeat()
                    .with_map(move |matches| {
                        let mut collected = FunctionBuilder::new(T.chain)
                            .add_arg(start.replace_wildcards_with_matches(matches))
                            .add_arg(end.replace_wildcards_with_matches(matches));
                        for (sequence, transpose) in
                            [(W_.a___, transpose_left), (W_.b___, transpose_right)]
                        {
                            let mut factors = match matches.get(sequence).unwrap() {
                                Match::Single(factor) => vec![*factor],
                                Match::Multiple(_, factors) => factors.to_vec(),
                                Match::FunctionName(_) => {
                                    unreachable!("a factor sequence is not a function name")
                                }
                            };
                            if transpose {
                                factors.reverse();
                            }
                            for factor in factors {
                                let factor = if transpose {
                                    factor.replace_map(|atom, _, output| {
                                        if let AtomView::Var(symbol) = atom {
                                            if symbol.get_symbol() == T.chain_in {
                                                **output = Atom::var(T.chain_out);
                                            } else if symbol.get_symbol() == T.chain_out {
                                                **output = Atom::var(T.chain_in);
                                            }
                                        }
                                    })
                                } else {
                                    factor.to_owned()
                                };
                                collected = collected.add_arg(factor);
                            }
                        }
                        // println!("{}", collected);
                        collected.finish()
                    });
            }
            result = result.normalize_chains();
            // Every join removes one chain; normalization may close it as a trace.
            if result == previous {
                return result;
            }
        }
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

    #[test]
    fn self_dual_common_end_chains_preserve_complex_matrix_components() {
        use crate::{dirac::spinor_matrix_structure, tensor::SymbolicTensor};
        use spenso::{
            network::{
                ExecutionResult, Sequential, SmallestDegree,
                library::{
                    function_lib::PanicMissingConcrete,
                    symbolic::{ExplicitKey, TensorLibrary},
                },
            },
            structure::{TensorStructure, abstract_index::AbstractIndex},
            tensors::{
                data::{DenseTensor, GetTensorData},
                parametric::{MixedTensor, ParamOrConcrete},
            },
        };
        test_initialize();
        let names = [
            spenso::tensor_symbol!("chain_numeric_A"),
            spenso::tensor_symbol!("chain_numeric_B"),
            spenso::tensor_symbol!("chain_numeric_C"),
        ];
        let data = [
            vec![Atom::num(1), Atom::i(), Atom::num(2), Atom::num(3)],
            vec![Atom::num(2), Atom::num(1), Atom::num(-1), Atom::num(4)],
            vec![Atom::num(3), Atom::num(-2), Atom::num(1), Atom::num(2)],
        ];
        let mut library =
            TensorLibrary::<MixedTensor<f64, ExplicitKey<AbstractIndex>>, AbstractIndex>::new();
        for (name, components) in names.into_iter().zip(data) {
            let key = spinor_matrix_structure::<AbstractIndex>(name, 2);
            library.insert_explicit(key.map_canonical(|structure| {
                ParamOrConcrete::param(
                    DenseTensor::from_storage_data(components, structure)
                        .unwrap()
                        .into(),
                )
            }));
        }
        let functions = PanicMissingConcrete::new_lib::<
            spenso::tensors::complex::RealOrComplexTensor<
                f64,
                spenso::network::parsing::ShadowedStructure<AbstractIndex>,
            >,
        >();
        let rep: LibraryRep = Bispinor {}.into();
        let left = bis!(2, chain_numeric_i);
        let right = bis!(2, chain_numeric_j);
        let dummy = bis!(2, chain_numeric_k);
        let factors = names.map(|name| function!(name, T.chain_in, T.chain_out));
        let [a, b, c] = &factors;
        // Independently multiplied matrices: (A B) C^T and (A B)^T C.
        // A contains i, so an accidental Hermitian conjugate is also detected.
        let cases = [
            (
                chain!(&left, &dummy, a, b) * chain!(&right, &dummy, c),
                [
                    Atom::num(4) - Atom::num(11) * Atom::i(),
                    Atom::num(4) + Atom::num(7) * Atom::i(),
                    Atom::num(-25),
                    Atom::num(29),
                ],
            ),
            (
                chain!(&dummy, &left, a, b) * chain!(&dummy, &right, c),
                [
                    Atom::num(7) - Atom::num(3) * Atom::i(),
                    Atom::num(-2) + Atom::num(2) * Atom::i(),
                    Atom::num(17) + Atom::num(12) * Atom::i(),
                    Atom::num(26) - Atom::num(8) * Atom::i(),
                ],
            ),
        ];
        for (expression, expected) in cases {
            let collected = expression.collect_chains(rep);
            assert_ne!(collected, expression);
            assert_eq!(collected.collect_chains(rep), collected);
            let mut reference_indices = None;
            for candidate in [expression, collected] {
                let mut network = SymbolicTensor::<AbstractIndex>::empty(candidate)
                    .to_network(&library)
                    .unwrap();
                network
                    .execute::<Sequential, SmallestDegree, _, _, _>(&library, &functions)
                    .unwrap();
                let ExecutionResult::Val(tensor) = network.result_tensor(&library).unwrap() else {
                    panic!("the matrix product must have two open indices");
                };
                let tensor = tensor.into_owned().try_into_parametric().unwrap();
                let indices = tensor.external_indices();
                if let Some(reference) = &reference_indices {
                    assert_eq!(&indices, reference);
                } else {
                    reference_indices = Some(indices);
                }
                for row in 0..2 {
                    for column in 0..2 {
                        assert_eq!(
                            tensor.get_owned([row, column]).unwrap(),
                            expected[2 * row + column]
                        );
                    }
                }
            }
        }
    }
}
