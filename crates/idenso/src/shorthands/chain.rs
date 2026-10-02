use std::sync::LazyLock;

use spenso::{
    chain, dualizable_, dualizable_dual_,
    network::{library::symbolic::ETS, parsing::ParseState, tags::SPENSO_TAG as T},
    self_dual_,
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibraryRep, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function,
    id::{Condition, Match, PatternRestriction, Replacement},
};
use symbolica_utils::PatternReplacement;

use crate::W_;
#[cfg(test)]
thread_local! {
    pub(crate) static CHAINIFY_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}
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

pub(crate) trait Chain {
    /// Join adjacent factors in one selected chain monomial.
    /// Distribution and foreign coefficients belong to the shared typed collector.
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
    fn join_chains(&self, representation: LibraryRep) -> Atom;

    /// Collect only the selected representation, optionally retaining closed words
    /// as exact chains. A disabled collection only converts existing closed chains.
    /// `retain_singletons` keeps isolated matrices in chain notation for algebra
    /// kernels; structural contraction otherwise preserves their indexed form.
    fn collect_chains(
        &self,
        representation: LibraryRep,
        collect: bool,
        traces: bool,
        retain_singletons: bool,
    ) -> Atom;

    /// Collapse identity words to metrics and closed chains to traces, e.g.
    ///
    /// `chain(rep(d,i), rep(d,i), ...)` becomes `trace(rep(d),...)`
    fn normalize_chains(&self) -> Atom;

    /// Normalize empty or identity-only words, optionally in one representation.
    fn normalize_chain_identities(&self, representation: Option<LibraryRep>) -> Atom;

    /// and single length chains back into the corresponding tensor
    fn undo_single_length(&self) -> Atom;

    /// Turn two compatible argument positions into a chain for [`Atom::join_chains`].
    /// A position may be an explicit slot, an unresolved representation, or a
    /// supplied compact rank-one tensor; supplied endpoints remain bound.
    ///
    /// `A(..,rep(d,i),rep(d,j),...)` becomes `chain(rep(d,i),rep(d,j),A(..,in,out...))`
    ///
    fn chainify(&self, representation: LibraryRep) -> Atom;
}

impl Chain for Atom {
    fn join_chains(&self, representation: LibraryRep) -> Atom {
        self.as_view().join_chains(representation)
    }

    fn collect_chains(
        &self,
        representation: LibraryRep,
        collect: bool,
        traces: bool,
        retain_singletons: bool,
    ) -> Atom {
        self.as_view()
            .collect_chains(representation, collect, traces, retain_singletons)
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

    fn normalize_chain_identities(&self, representation: Option<LibraryRep>) -> Atom {
        self.as_view().normalize_chain_identities(representation)
    }
}
impl<'a> Chain for AtomView<'a> {
    fn join_chains(&self, representation: LibraryRep) -> Atom {
        self.collect_chains(representation, true, true, true)
            .normalize_chains()
    }

    fn collect_chains(
        &self,
        representation: LibraryRep,
        collect: bool,
        traces: bool,
        retain_singletons: bool,
    ) -> Atom {
        // Match full slot compatibility, including dimensions and duality.
        // Names alone never authorize a closure.
        let closed = |node: AtomView<'_>, slots: &mut SlotMatcher| {
            use spenso::structure::slot::IsAbstractSlot;
            let AtomView::Fun(chain) = node else {
                return None;
            };
            if chain.get_symbol() != T.chain || chain.get_nargs() < 2 {
                return None;
            }
            let mut arguments = chain.iter();
            let start = slots
                .parse::<LibraryRep, AbstractIndex>(arguments.next().unwrap())
                .ok()?;
            let end = slots
                .parse::<LibraryRep, AbstractIndex>(arguments.next().unwrap())
                .ok()?;
            if start.rep().rep != representation
                || !start.rep().matches(&end.rep())
                || start.aind != end.aind
                || !(representation.is_self_dual()
                    || representation.is_base() && end.rep().rep.is_dual())
            {
                return None;
            }
            Some(spenso::shadowing::trace(
                start.rep().to_symbolic([]),
                arguments,
            ))
        };
        let normalize = |value: Atom| {
            if !traces {
                return value;
            }
            let mut slots = SlotMatcher::default();
            let mut closes = false;
            value.visitor(&mut |node| {
                closes |= closed(node, &mut slots).is_some();
                !closes
            });
            if !closes {
                return value;
            }
            value.replace_map(|node, _, output| {
                if let Some(trace) = closed(node, &mut slots) {
                    **output = trace;
                }
            })
        };
        if !collect {
            return normalize(self.to_owned());
        }
        let source = self.chainify(representation);
        // Every local join requires a chain. The typed collector has already
        // selected this monomial and retained unrelated factors outside it.
        if !source.contains_symbol(T.chain) {
            return source;
        }

        let in_index = Atom::var(W_.i_);
        let dummy_out = representation.dual().to_symbolic([W_.d_, W_.j_]);
        let dummy_in = representation.to_symbolic([W_.d_, W_.j_]);
        let out_index = Atom::var(W_.k_);
        let endpoints = W_.i_.filter_cmp(W_.d_, move |endpoint, dimension| {
            let (Match::Single(endpoint), Match::Single(dimension)) = (endpoint, dimension) else {
                return false;
            };
            SlotMatcher::default()
                .port_representation(*endpoint)
                .is_some_and(|(rep, dim)| rep == representation && dim == *dimension)
        }) & W_.k_.filter_cmp(W_.d_, move |endpoint, dimension| {
            let (Match::Single(endpoint), Match::Single(dimension)) = (endpoint, dimension) else {
                return false;
            };
            SlotMatcher::default()
                .port_representation(*endpoint)
                .is_some_and(|(rep, dim)| rep == representation.dual() && dim == *dimension)
        });

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
        // Joins match two chains of one product. Each product with two or
        // more chains reaches its fixed point among its own chain factors,
        // children first, instead of matching every product of the whole
        // expression once per pass.
        let mut joined = false;
        let mut fixed_point = |mut result: Atom| loop {
            let previous = result.clone();
            for (product, transpose_left, transpose_right) in &joins {
                let start = in_index.to_pattern();
                let end = out_index.to_pattern();
                let transpose_left = *transpose_left;
                let transpose_right = *transpose_right;
                result = result
                    .replace(product.to_pattern())
                    .when(endpoints.clone())
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
                                    let mut slots = SlotMatcher::default();
                                    let mut opaque_level = None;
                                    factor.replace_map(|atom, context, output| {
                                        if opaque_level
                                            .is_some_and(|level| context.function_level > level)
                                        {
                                            return;
                                        }
                                        opaque_level = None;
                                        if !matches!(slots.classify(atom), SlotMatch::Other) {
                                            opaque_level = Some(context.function_level);
                                            return;
                                        }
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
                        collected.finish()
                    });
            }
            joined |= result != previous;
            result = normalize(result);
            // Every join removes one chain; normalization may close it as a trace.
            if result == previous {
                return result;
            }
        };
        let joined_products = joining_products(source.as_view(), &mut fixed_point);
        // Close single chains outside joined products too.
        let result = normalize(joined_products.unwrap_or_else(|| source.clone()));
        if retain_singletons {
            // The existing pass has already collected every matrix.
            // Algebra kernels keep those words, so neither unfold them
            // here nor run chainify again to reconstruct them afterward.
            return result;
        }
        if !joined && (!traces || result == source) {
            // No connection was joined or closed. Keep the original
            // notation, including existing one-factor chains: unfolding
            // them here competes with the algebra kernels' collection.
            return self.to_owned();
        }
        // An isolated matrix has no product to collect. Retain its exact
        // indexed form; a closed word must survive when traces are disabled.
        result.replace_map(|node, _, output| {
            let AtomView::Fun(chain) = node else {
                return;
            };
            if chain.get_symbol() != T.chain || chain.get_nargs() != 3 {
                return;
            }
            let mut arguments = chain.iter();
            let start = arguments.next().unwrap();
            let end = arguments.next().unwrap();
            use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
            let mut slots = SlotMatcher::default();
            let Ok(start) = slots.parse::<LibraryRep, AbstractIndex>(start) else {
                return;
            };
            let Ok(end) = slots.parse::<LibraryRep, AbstractIndex>(end) else {
                return;
            };
            if start.rep().rep == representation && start.aind != end.aind {
                **output = node.undo_single_length();
            }
        })
    }

    fn chainify(&self, representation: LibraryRep) -> Atom {
        // Both endpoints must use this representation. Terminal metric/epsilon
        // expressions often have no such slots; avoid building a replacement
        // and invoking its matcher at every unrelated function in that case.
        #[cfg(test)]
        CHAINIFY_CALLS.with(|count| count.set(count.get() + 1));
        if !self.contains_symbol(representation.symbol()) {
            return self.to_owned();
        }

        let in_index = Atom::var(W_.i_);
        let out_index = Atom::var(W_.j_);
        let endpoints = chainify_endpoints(representation);
        let replacement = chainify_rule(representation);
        let opaque_heads = [T.chain, T.trace];
        let existing_chain = chain!(&in_index, &out_index, W_.a___).to_pattern();
        // An unchanged pruning assignment in replace_map rebuilds ancestors
        // and invokes their normalizers. Keep borrowed children until this
        // existing chain rule actually changes a node.
        fn rewrite<'a>(
            atom: AtomView<'a>,
            slots: &mut SlotMatcher,
            opaque_heads: &[Symbol],
            apply: &impl Fn(AtomView<'_>) -> Atom,
            open: &impl Fn(symbolica::atom::representation::FunView<'_>) -> Option<Atom>,
        ) -> AtomOrView<'a> {
            if !matches!(slots.classify(atom), SlotMatch::Other) {
                return atom.into();
            }
            let children = match atom {
                AtomView::Fun(function) => {
                    if let Some(opened) = open(function) {
                        return rewrite(opened.as_view(), slots, opaque_heads, apply, open)
                            .into_owned()
                            .into();
                    }
                    if opaque_heads.contains(&function.get_symbol()) {
                        return atom.into();
                    }
                    let changed = apply(atom);
                    if changed.as_view() != atom {
                        return changed.into();
                    }
                    function.iter().collect::<Vec<_>>()
                }
                AtomView::Add(sum) => sum.iter().collect(),
                AtomView::Mul(product) => product.iter().collect(),
                AtomView::Pow(power) => {
                    let (base, exponent) = power.get_base_exp();
                    vec![base, exponent]
                }
                _ => return atom.into(),
            };
            let children = children
                .into_iter()
                .map(|child| rewrite(child, slots, opaque_heads, apply, open))
                .collect::<Vec<_>>();
            if children
                .iter()
                .all(|child| matches!(child, AtomOrView::View(_)))
            {
                return atom.into();
            }
            match atom {
                AtomView::Fun(function) => {
                    let mut output = FunctionBuilder::new(function.get_symbol());
                    for child in children {
                        output = output.add_arg(child);
                    }
                    output.finish().into()
                }
                AtomView::Add(_) => Atom::add_many(children).into(),
                AtomView::Mul(_) => Atom::mul_many(children).into(),
                AtomView::Pow(_) => children[0].as_view().pow(children[1].as_view()).into(),
                _ => unreachable!("only compound expressions have children"),
            }
        }
        let state = std::cell::OnceCell::new();
        rewrite(
            *self,
            &mut SlotMatcher::default(),
            &opaque_heads,
            &|atom| match chainified_slot_function(atom, representation) {
                Some(chained) => chained.unwrap_or_else(|| atom.to_owned()),
                None => atom.replace_multiple([&replacement]),
            },
            &|function| {
                if function.get_symbol() != T.dot && function.get_symbol() != ETS.metric {
                    return None;
                }
                let (left, right, _, _) = spenso::structure::slot::SlotMatcher::default()
                    .compact_inner_product_parts::<AbstractIndex>(function)?;
                // Select with the existing endpoint rule, including its opaque-head
                // and placeholder restrictions. Existing chains already own their
                // endpoints.
                // Either operand may carry the selected endpoints, as in gamma·p.
                let selected = [left, right].into_iter().any(|operand| {
                    operand
                        .pattern_match(&replacement.pat, replacement.conditions.as_ref(), None)
                        .next()
                        .is_some()
                        || operand
                            .pattern_match(&existing_chain, Some(&endpoints), None)
                            .next()
                            .is_some()
                });
                if !selected {
                    return None;
                }
                let state = state.get_or_init(|| {
                    let state = ParseState::<AbstractIndex>::default();
                    state.reserve_indices(*self);
                    state
                });
                state.materialize_inner_product(function)
            },
        )
        .into_owned()
    }

    fn normalize_chains(&self) -> Atom {
        self.normalize_chain_identities(None)
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
    fn join_chains_preserves_factorized_scalars_without_chains() {
        test_initialize();
        let scalar = parse_lit!((a + b) ^ 8 * (c + d) ^ 8 / (1 + a * c + b * d));
        for representation in [Bispinor {}.into(), ColorFundamental {}.into()] {
            // Structural equality checks the original factorization, not only
            // equality after expanding or evaluating the scalar expression.
            assert_eq!(scalar.join_chains(representation), scalar);
        }
    }

    #[test]
    fn chainify_preserves_terminal_tensors_and_scalar_factorization() {
        test_initialize();
        let terminal = parse_lit!(
            ((x + y) ^ 8)
                * (g(mink(4, mu), p(mink(4)))
                    + epsilon(mink(4, mu), mink(4, nu), mink(4, rho), mink(4, sigma)))
                / (1 + x * y),
            default_namespace = "spenso"
        );
        for representation in [Bispinor {}.into(), ColorFundamental {}.into()] {
            assert_eq!(terminal.chainify(representation), terminal);
        }
    }

    #[test]
    fn chainify_finds_generic_tensors_inside_sums_and_powers() {
        test_initialize();
        let scalar = parse_lit!((x + y) ^ 8);
        let head = symbol!("chainify_generic_tensor");
        let metadata = symbol!("chainify_metadata");
        for representation in [
            LibraryRep::from(Bispinor {}),
            LibraryRep::from(ColorFundamental {}),
            LibraryRep::from(ColorFundamental {}).dual(),
        ] {
            let dimension = parse_lit!(n ^ 2 - 1);
            let start = representation.to_symbolic([dimension.clone(), parse_lit!(label(a))]);
            let end = representation
                .dual()
                .to_symbolic([dimension, parse_lit!(label(b))]);
            let tensor = function!(head, &scalar, &start, metadata, &end);
            let chained = chain!(
                start,
                end,
                function!(head, &scalar, T.chain_in, metadata, T.chain_out)
            );
            let input = (&scalar + tensor).pow(2);
            let expected = (&scalar + chained).pow(2);
            assert_eq!(input.chainify(representation), expected);
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
        let collected = normalized.join_chains(rep);

        assert_snapshot!(collected.to_bare_ordered_string(), @"trace(bis(4),cyclic(gamma(in,out,mink(4,mu)),gamma(in,out,p(3,mink(4))),gamma(in,out,p(2,mink(4)))))");
    }

    /// The self-dual fast paths of chainify and of the single-length undo
    /// give what their rules give: structure constants with slots in every
    /// order, mixed dimensions, a slot of another representation, and the
    /// undo of each result.
    #[test]
    fn self_dual_fast_paths_agree_with_their_rules() {
        test_initialize();
        let adjoint = crate::representations::ColorAdjoint {};
        let slot = |dimension: i64, name: &str| {
            adjoint.to_symbolic([
                Atom::num(dimension),
                Atom::var(symbol!(&format!("idenso::chain_fast::{name}"))),
            ])
        };
        let fundamental = ColorFundamental {}
            .to_symbolic([Atom::num(3), Atom::var(symbol!("idenso::chain_fast::i"))]);
        let head = symbol!("idenso::chain_fast::T");
        let [a, b, c, d] = ["a", "b", "c", "d"].map(|name| slot(8, name));
        let small = slot(3, "s");
        let mut inputs = vec![
            crate::color_f!(&a, &b, &c),
            crate::color_f!(&c, &a, &b),
            crate::color_f!(&b, &c, &a),
        ];
        for arguments in [
            vec![&a, &b, &c, &d],
            vec![&a, &small, &b],
            vec![&small, &a, &b],
            vec![&a, &b, &small],
            vec![&a, &fundamental, &b],
            vec![&a, &small],
            vec![&a],
        ] {
            inputs.push(
                arguments
                    .into_iter()
                    .fold(FunctionBuilder::new(head), |builder, argument| {
                        builder.add_arg(argument)
                    })
                    .finish(),
            );
        }
        let representation = LibraryRep::from(adjoint);
        let rule = chainify_rule(representation);
        for input in inputs {
            let fast = chainified_slot_function(input.as_view(), representation)
                .expect("explicit slots are decided")
                .unwrap_or_else(|| input.clone());
            assert_eq!(fast, input.replace_multiple([&rule]), "{input}");
            if let Some(undone) = undone_self_dual_single_length(fast.as_view()) {
                assert_eq!(
                    undone,
                    fast.replace_multiple_repeat(SINGLE_LENGTH_NORM.as_ref()),
                    "{fast}"
                );
                assert_eq!(undone, input, "{fast}");
            }
        }
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

        assert_snapshot!(chains.join_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,c),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)),gamma(in,out,p(1,mink(4))))");
    }

    #[test]
    fn join_chains_preserves_opaque_factorized_spectators() {
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
        assert_eq!(input.join_chains(Bispinor {}.into()), &spectator * composed);
        assert_eq!(input.join_chains(ColorFundamental {}.into()), input);
    }

    #[test]
    fn joining_indexed_words_keeps_singletons_without_repeating_chainify() {
        test_initialize();
        let rep = Bispinor {}.into();
        let [a, b, c, d, e] =
            [98971, 98972, 98973, 98974, 98975].map(|index| bis!(4, Atom::num(index)));
        let [mu, nu, rho] = [98976, 98977, 98978].map(|index| spenso::mink!(4, Atom::num(index)));
        let first = gamma!(&a, &b, &mu);
        let second = gamma!(&b, &c, &nu);
        let isolated = gamma!(&d, &e, &rho);
        let spectator = parse_lit!((chain_collection_x + chain_collection_y) ^ 9);
        let expected_word = chain!(&a, &c, gamma!(&mu), gamma!(&nu));
        let expected_single = chain!(&d, &e, gamma!(&rho));
        let input = &spectator * &first * &second * &isolated;
        CHAINIFY_CALLS.with(|calls| calls.set(0));
        let result = input.join_chains(rep);
        assert_eq!(result, &spectator * &expected_word * &expected_single);
        CHAINIFY_CALLS.with(|calls| assert_eq!(calls.get(), 1));
        assert_eq!(result.join_chains(rep), result);
        assert_eq!(
            input.collect_chains(rep, true, true, false),
            &spectator * &expected_word * &isolated,
        );
        CHAINIFY_CALLS.with(|calls| calls.set(0));
        assert_eq!(isolated.join_chains(rep), expected_single);
        CHAINIFY_CALLS.with(|calls| assert_eq!(calls.get(), 1));

        // The same shared pass closes a cycle, preserving its ordered factors.
        let closed = &first * gamma!(&b, &a, &nu);
        CHAINIFY_CALLS.with(|calls| calls.set(0));
        assert_eq!(
            closed.join_chains(rep),
            trace!(bis!(4), gamma!(&mu), gamma!(&nu))
        );
        CHAINIFY_CALLS.with(|calls| assert_eq!(calls.get(), 1));
    }

    #[test]
    fn structural_collection_preserves_existing_isolated_chain_notation() {
        test_initialize();
        let source = parse!(
            "chain(bis(4,i),bis(4,j),spenso::gamma(in,out,mink(4,mu)))",
            default_namespace = "spenso"
        );
        let representation = Bispinor {}.into();
        assert_eq!(
            source.collect_chains(representation, true, true, false),
            source
        );
        let explicit = source.undo_single_length();
        assert_eq!(
            explicit.collect_chains(representation, true, true, false),
            explicit
        );
    }

    #[test]
    fn join_chains_keeps_other_representations_factored() {
        // This test selects family prerequisites only, preserving its explicit output scope.

        test_initialize();
        let spin = parse!(
            "chain(bis(4, i), bis(4, j), spenso::gamma(in, out, mink(4, mu)))
                * chain(bis(4, j), bis(4, k), spenso::gamma(in, out, mink(4, nu)))",
            default_namespace = "spenso"
        );
        let spin_composed = parse!(
            "chain(
                bis(4, i), bis(4, k),
                spenso::gamma(in, out, mink(4, mu)), spenso::gamma(in, out, mink(4, nu))
            )",
            default_namespace = "spenso"
        );
        let color = parse!(
            "chain(cof(3, a), dind(cof(3, b)), t(coad(8, alpha), in, out))
                * (chain(cof(3, b), dind(cof(3, c)), t(coad(8, beta), in, out), t(coad(8, delta), in, out))
                    + chain(cof(3, b), dind(cof(3, c)), t(coad(8, delta), in, out), t(coad(8, beta), in, out)))",
            default_namespace = "spenso"
        );
        let color_composed = parse!(
            "chain(
                cof(3, a), dind(cof(3, c)),
                t(coad(8, alpha), in, out), t(coad(8, beta), in, out), t(coad(8, delta), in, out)
            ) + chain(
                cof(3, a), dind(cof(3, c)),
                t(coad(8, alpha), in, out), t(coad(8, delta), in, out), t(coad(8, beta), in, out)
            )",
            default_namespace = "spenso"
        );
        // Different free adjoint labels do not form one typed tensor sum.
        let incompatible = parse!(
            "chain(cof(3, a), dind(cof(3, b)), t(coad(8, alpha), in, out))
                * (chain(cof(3, b), dind(cof(3, c)), t(coad(8, beta), in, out))
                    + chain(cof(3, b), dind(cof(3, c)), t(coad(8, delta), in, out)))",
            default_namespace = "spenso"
        );
        assert!(crate::tensor::SymbolicTensor::infer(&spin * incompatible).is_err());
        let input = &spin * &color;

        assert_eq!(color.join_chains(Bispinor {}.into()), color);
        assert_eq!(
            input.join_chains(Bispinor {}.into()),
            &spin_composed * &color
        );
        // Local chain joining does not distribute sums. The shared typed
        // collector owns selection and retains the foreign spin factor intact.
        assert_eq!(input.join_chains(ColorFundamental {}.into()), input);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer(input)
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    contract: crate::tensor::AlgebraContraction::None,
                    color: Some(crate::color::ColorSimplifySettings {
                        ..Default::default()
                    }),
                    ..Default::default()
                })
                .unwrap()
                .expression,
            &spin * &color_composed
        );
    }

    #[test]
    fn join_chains_selects_ports_instead_of_payload_representations() {
        test_initialize();
        let input = parse!(
            "chain(cof(3, a), dind(cof(3, b)), mixed(in, out, bis(4, i), bis(4, j)))
                * (chain(cof(3, b), dind(cof(3, c)), mixed(in, out, bis(4, j), bis(4, k)))
                    + chain(cof(3, b), dind(cof(3, c)), other(in, out, bis(4, j), bis(4, k))))",
            default_namespace = "spenso"
        );

        // The explicit spin slots belong to the mixed tensor payload; they do
        // not turn its fundamental-color chain ports into spin-chain ports.
        assert_eq!(input.join_chains(Bispinor {}.into()), input);
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
    fn join_chains_leaves_single_gamma_chain_for_explicit_undo() {
        let r = TestReps::new();
        let chain = chain!(
            slot!(r.bis4, a),
            slot!(r.bis4, b),
            gamma!(slot!(r.mink4, mu)),
        );
        let rep = Bispinor {}.into();

        assert_snapshot!(chain.join_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))");
    }

    #[test]
    fn join_chains_composes_common_prefix_terms_before_factoring() {
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

        assert_snapshot!(chains.join_chains(rep).to_bare_ordered_string(), @"chain(bis(4,a),bis(4,c),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,nu)))+chain(bis(4,a),bis(4,d),gamma(in,out,mink(4,mu)),gamma(in,out,mink(4,rho)))");
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
            let collected = expression.join_chains(rep);
            assert_ne!(collected, expression);
            assert_eq!(collected.join_chains(rep), collected);
            let mut reference_indices = None;
            for candidate in [expression, collected] {
                let mut network =
                    SymbolicTensor::<OrderedStructure<LibraryRep, AbstractIndex>>::empty(candidate)
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
