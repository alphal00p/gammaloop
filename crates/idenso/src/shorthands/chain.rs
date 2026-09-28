use std::sync::LazyLock;

use spenso::{
    chain, dualizable_, dualizable_dual_,
    network::{library::symbolic::ETS, parsing::ParseState, tags::SPENSO_TAG as T},
    self_dual_,
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibraryRep, RepName},
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder, Symbol},
    function,
    id::{Match, Replacement},
};
use symbolica_utils::PatternReplacement;

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

    /// Collapse identity words to metrics and closed chains to traces, e.g.
    ///
    /// `chain(rep(d,i), rep(d,i), ...)` becomes `trace(rep(d),...)`
    fn normalize_chains(&self) -> Atom;

    /// and single length chains back into the corresponding tensor
    fn undo_single_length(&self) -> Atom;

    /// turns tensors with two indices of the representation into chain expressions, for collecting using [`Atom::join_chains`]
    ///
    /// `A(..,rep(d,i),rep(d,j),...)` becomes `chain(rep(d,i),rep(d,j),A(..,in,out...))`
    ///
    fn chainify(&self, representation: LibraryRep) -> Atom;
}

impl Chain for Atom {
    fn join_chains(&self, representation: LibraryRep) -> Atom {
        self.as_view().join_chains(representation)
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
    fn join_chains(&self, representation: LibraryRep) -> Atom {
        // Every local join requires a chain. The typed collector has already
        // selected this monomial and retained unrelated factors outside it.
        if !self.contains_symbol(T.chain) {
            return self.to_owned();
        }

        let in_index = representation.to_symbolic([W_.d_, W_.i_]);
        let dummy_out = representation.dual().to_symbolic([W_.d_, W_.j_]);
        let dummy_in = representation.to_symbolic([W_.d_, W_.j_]);
        let out_index = representation.dual().to_symbolic([W_.d_, W_.k_]);

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
        let mut result = self.to_owned();
        loop {
            let previous = result.clone();
            for (product, transpose_left, transpose_right) in &joins {
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
        // Both endpoints must use this representation. Terminal metric/epsilon
        // expressions often have no such slots; avoid building a replacement
        // and invoking its matcher at every unrelated function in that case.
        if !self.contains_symbol(representation.symbol()) {
            return self.to_owned();
        }

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
        let tensor_alias = spenso::tensor_symbol!("idenso::tensor_alias");
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
            W_.a_.filter_match(move |a| {
                matches!(a, Match::FunctionName(a)
                    if *a != T.chain && *a != T.trace && *a != ETS.metric && *a != tensor_alias)
            }) & W_.a___.filter_match(no_ports)
                & W_.b___.filter_match(no_ports)
                & W_.c___.filter_match(no_ports),
        )
        .max_level(0);
        let opaque_heads = [T.chain, T.trace, tensor_alias];
        let existing_chain = chain!(&in_index, &out_index, W_.a___).to_pattern();
        // An unchanged pruning assignment in replace_map rebuilds ancestors
        // and invokes their normalizers. Keep borrowed children until this
        // existing chain rule actually changes a node.
        fn rewrite<'a>(
            atom: AtomView<'a>,
            opaque_heads: &[Symbol],
            apply: &impl Fn(AtomView<'_>) -> Atom,
            open: &impl Fn(symbolica::atom::representation::FunView<'_>) -> Option<Atom>,
        ) -> AtomOrView<'a> {
            let children = match atom {
                AtomView::Fun(function) => {
                    if let Some(opened) = open(function) {
                        return rewrite(opened.as_view(), opaque_heads, apply, open)
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
                .map(|child| rewrite(child, opaque_heads, apply, open))
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
            &opaque_heads,
            &|atom| atom.replace_multiple([&replacement]),
            &|function| {
                if function.get_symbol() != T.dot && function.get_symbol() != ETS.metric {
                    return None;
                }
                let (left, right, _, _) = spenso::structure::slot::SlotMatcher::default()
                    .compact_inner_product_parts::<AbstractIndex>(function)?;
                // This raw operation cannot register a new literal alias spelling.
                // Either operand may be relabelled by compact materialization, even
                // when only the other one supplies the selected chain endpoints.
                if left.contains_symbol(tensor_alias) || right.contains_symbol(tensor_alias) {
                    return None;
                }
                // Select with the existing endpoint rule, including its opaque-head
                // and placeholder restrictions. Existing chains already own their
                // endpoints; tensor alias ports must stay literal registry keys.
                // Either operand may carry the selected endpoints, as in gamma·p.
                let selected = [left, right].into_iter().any(|operand| {
                    operand
                        .pattern_match(&replacement.pat, replacement.conditions.as_ref(), None)
                        .next()
                        .is_some()
                        || operand
                            .pattern_match(&existing_chain, None, None)
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
        self.replace_map(|atom, _, out| {
            let AtomView::Fun(chain) = atom else { return };
            if chain.get_symbol() != T.chain || chain.get_nargs() < 2 {
                return;
            }
            let mut factors = chain.iter();
            let start = factors.next().unwrap();
            let end = factors.next().unwrap();
            // Empty words and words consisting only of identity lines have
            // the physical endpoint metric, in every representation.
            if factors.all(|factor| {
                let AtomView::Fun(metric) = factor else { return false };
                if metric.get_symbol() != ETS.metric || metric.get_nargs() != 2 {
                    return false;
                }
                let mut endpoints = metric.iter();
                matches!(endpoints.next(), Some(AtomView::Var(v)) if v.get_symbol() == T.chain_in)
                    && matches!(endpoints.next(), Some(AtomView::Var(v)) if v.get_symbol() == T.chain_out)
            }) {
                **out = function!(ETS.metric, start, end);
            }
        })
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
    use spenso::structure::OrderedStructure;
    use spenso::{chain, slot};
    use symbolica::{parse, parse_lit, symbol};
    use symbolica_utils::AtomPrintExt;

    use crate::representations::{Bispinor, ColorFundamental};
    use crate::test_support::{TestReps, test_initialize};
    use crate::{bis, gamma};

    use super::*;

    #[test]
    fn compact_dot_chainifies_mixed_gamma_vector_operands() {
        use crate::{dirac::GammaSimplifySettings, tensor::SymbolicTensor};
        use spenso::structure::partial::PartialStructure;

        let reps = test_initialize();
        let left = bis!(4, 93701);
        let inner = bis!(4, 93702);
        let right = bis!(4, 93703);
        for compact in [reps.mink4.to_symbolic([]), reps.mink_d.to_symbolic([])] {
            let momentum = spenso::p!(&compact);
            let expected = function!(T.dot, &momentum, &momentum) * g!(&left, &right);
            for chained in [false, true] {
                let gamma = |start: &Atom, end: &Atom| {
                    if chained {
                        chain!(start, end, gamma!(&compact))
                    } else {
                        gamma!(start.as_view(), end.as_view(), &compact)
                    }
                };
                let expression = function!(T.dot, gamma(&left, &inner), &momentum)
                    * function!(T.dot, &momentum, gamma(&inner, &right));
                let source = SymbolicTensor::<PartialStructure>::infer(expression.clone()).unwrap();
                assert_eq!(source.rank(), 2);
                assert_eq!(
                    source.structure,
                    SymbolicTensor::<PartialStructure>::infer(g!(&left, &right))
                        .unwrap()
                        .structure
                );
                let chainified = expression.chainify(Bispinor {}.into());
                assert!(!chainified.contains_symbol(T.dot));
                let result = source
                    .simplify_gamma(GammaSimplifySettings::default())
                    .unwrap();
                assert_eq!(result.root().structure, source.structure);
                assert_eq!(
                    result.resolved().unwrap().to_dots().unwrap().expression,
                    expected,
                    "compact={compact}, chained={chained}\ninput={expression}\nchainify={chainified}\njoin={}\nmetric-first={:?}\nchains={:?}\nexplicit-gamma={:?}",
                    chainified.join_chains(Bispinor {}.into()),
                    source
                        .contract(Default::default())
                        .and_then(|value| value.resolved())
                        .map(|value| value.expression),
                    source
                        .simplify_gamma(GammaSimplifySettings {
                            output: crate::dirac::GammaOutput::Chains,
                            ..Default::default()
                        })
                        .and_then(|value| value.resolved())
                        .map(|value| value.expression),
                    source
                        .undo_dots()
                        .and_then(|value| value.simplify_gamma(GammaSimplifySettings::default()))
                        .and_then(|value| value.resolved())
                        .map(|value| value.expression),
                );
                let aliased = source.alias_handle().unwrap();
                let dag = std::sync::Arc::new(
                    aliased
                        .clone()
                        .with_aliases([(aliased, source.clone())])
                        .unwrap(),
                );
                assert_eq!(
                    dag.simplify_gamma(GammaSimplifySettings::default())
                        .unwrap()
                        .resolved()
                        .unwrap()
                        .to_dots()
                        .unwrap()
                        .expression,
                    expected,
                );
            }
        }
    }

    #[test]
    fn compact_dot_chainification_keeps_scalar_and_foreign_products_opaque() {
        use std::sync::{Arc, Mutex};

        let reps = test_initialize();
        let compact = reps.mink4.to_symbolic([]);
        let scalar_dot = function!(T.dot, spenso::p!(&compact), spenso::q!(&compact));
        let foreign = function!(
            symbol!("chainify_foreign_dot"),
            spenso::p!(&compact),
            spenso::q!(&compact)
        );
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let callback = symbolica::symbol!("chainify_compact_scalar_callback"; Scalar;
            norm = move |node, _| recorded.lock().unwrap().push(node.to_owned())
        );
        let spectator = function!(callback, &scalar_dot) * &foreign;
        let raw = gamma!(bis!(4, 93711), bis!(4, 93712), spenso::mink!(4, 93713));
        let input = &spectator * &raw;
        let expected = &spectator * raw.chainify(Bispinor {}.into());
        calls.lock().unwrap().clear();
        assert_eq!(input.chainify(Bispinor {}.into()), expected);
        assert!(calls.lock().unwrap().is_empty());
    }

    #[test]
    fn compact_dot_chainification_preserves_weighted_scalar_dot_coefficients() {
        use std::sync::{Arc, Mutex};

        let reps = test_initialize();
        let compact = reps.mink4.to_symbolic([]);
        let scalar_dot = function!(T.dot, spenso::p!(&compact), spenso::q!(&compact));
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let callback = symbolica::symbol!("chainify_weighted_compact_callback"; Scalar;
            norm = move |node, _| recorded.lock().unwrap().push(node.to_owned())
        );
        let coefficient = function!(callback, &scalar_dot)
            * (Atom::var(symbol!("compact_weight_x")) + Atom::var(symbol!("compact_weight_y")))
                .pow(3);
        let input = function!(
            T.dot,
            &coefficient * gamma!(bis!(4, 93731), bis!(4, 93732), &compact),
            spenso::p!(&compact)
        );
        calls.lock().unwrap().clear();
        let actual = input.chainify(Bispinor {}.into());
        assert_ne!(actual, input);
        assert!(calls.lock().unwrap().is_empty());
        assert!(
            actual
                .replace(coefficient.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );
        let mut dots = Vec::new();
        actual.visitor(&mut |node| {
            if matches!(node, AtomView::Fun(fun) if fun.get_symbol() == T.dot) {
                dots.push(node.to_owned());
            }
            true
        });
        assert_eq!(dots, vec![scalar_dot]);
    }

    #[test]
    fn compact_dot_chain_words_retain_explicit_spectators() {
        use crate::{dirac::GammaSimplifySettings, tensor::SymbolicTensor};
        use spenso::structure::partial::PartialStructure;

        let reps = test_initialize();
        let left = bis!(4, 93741);
        let right = bis!(4, 93742);
        for rep in [reps.mink4, reps.mink_d] {
            let compact = rep.to_symbolic([]);
            let spectator = rep.to_symbolic([Atom::num(93743)]);
            let momentum = spenso::p!(&compact);
            let forward = chain!(&left, &right, gamma!(&compact), gamma!(&spectator));
            let reverse = chain!(&left, &right, gamma!(&spectator), gamma!(&compact));
            let expression =
                function!(T.dot, forward, &momentum) + function!(T.dot, &momentum, reverse);
            let source = SymbolicTensor::<PartialStructure>::infer(expression.clone()).unwrap();
            assert_eq!(source.rank(), 3);
            let chained = expression.chainify(Bispinor {}.into());
            assert!(!chained.contains_symbol(T.dot));
            assert_eq!(
                SymbolicTensor::<PartialStructure>::infer(chained)
                    .unwrap()
                    .structure,
                source.structure
            );
            let actual = source
                .simplify_gamma(GammaSimplifySettings::canonical())
                .unwrap();
            assert_eq!(actual.root().structure, source.structure);
            assert_eq!(
                actual.resolved().unwrap().expression,
                Atom::num(2) * spenso::p!(&spectator) * g!(&left, &right)
            );
        }
    }

    #[test]
    fn compact_dot_chainification_preserves_literal_alias_ports() {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::PartialStructure;

        let reps = test_initialize();
        let compact = reps.mink4.to_symbolic([]);
        let body = SymbolicTensor::<PartialStructure>::infer(gamma!(
            bis!(4, 93721),
            bis!(4, 93722),
            &compact
        ))
        .unwrap();
        let handle = body.alias_handle().unwrap();
        let input = function!(T.dot, &handle.expression, spenso::p!(&compact));
        assert_eq!(input.chainify(Bispinor {}.into()), input);

        let gamma = gamma!(bis!(4, 93722), bis!(4, 93723), &compact);
        for (input, expected) in [
            (
                function!(T.dot, &handle.expression, &gamma),
                function!(
                    T.dot,
                    &handle.expression,
                    gamma.chainify(Bispinor {}.into())
                ),
            ),
            (
                function!(T.dot, &gamma, &handle.expression),
                function!(
                    T.dot,
                    gamma.chainify(Bispinor {}.into()),
                    &handle.expression
                ),
            ),
        ] {
            let actual = input.chainify(Bispinor {}.into());
            assert_eq!(actual, expected);
            assert_eq!(actual.chainify(Bispinor {}.into()), actual);
        }

        // A scalar alias is an opaque coefficient, never a compact port owner.
        let scalar =
            SymbolicTensor::<PartialStructure>::infer(Atom::var(symbol!("compact_alias_weight")))
                .unwrap()
                .alias_handle()
                .unwrap();
        let weighted = function!(T.dot, &scalar.expression * &gamma, spenso::p!(&compact));
        let actual = weighted.chainify(Bispinor {}.into());
        assert!(!actual.contains_symbol(T.dot));
        assert!(
            actual
                .replace(scalar.expression.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );
    }

    #[test]
    fn chain_collection_keeps_literal_alias_ports_opaque() {
        use crate::tensor::SymbolicTensor;
        use spenso::structure::partial::PartialStructure;

        test_initialize();
        let body = SymbolicTensor::<PartialStructure>::infer(gamma!(
            bis!(4, 93601),
            bis!(4, 93603),
            spenso::mink!(4, 93605)
        ))
        .unwrap();
        let handle = body.alias_handle().unwrap();
        assert_eq!(
            handle.expression.chainify(Bispinor {}.into()),
            handle.expression
        );
        let original = handle.expression.clone();
        let source = std::sync::Arc::new(
            handle
                .clone()
                .with_aliases([(handle, body.clone())])
                .unwrap(),
        );
        let collected = source
            .simplify_gamma(crate::dirac::GammaSimplifySettings {
                output: crate::dirac::GammaOutput::Chains,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(collected.root().expression, original);
        assert_eq!(
            collected
                .resolved()
                .unwrap()
                .expression
                .undo_single_length(),
            body.expression
        );
    }

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
    fn chainify_pruning_keeps_scalar_normalizers_untouched() {
        use crate::shorthands::UndoShorthands;
        use spenso::structure::abstract_index::AbstractIndex;
        use std::sync::{Arc, Mutex};
        let r = test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let callback = symbolica::symbol!("chainify_scalar_callback"; Scalar;
            norm = move |node, _| recorded.lock().unwrap().push(node.to_owned())
        );
        let word = trace!(
            r.bis4.to_symbolic([]),
            gamma!(slot!(r.mink4, 99411)),
            gamma!(slot!(r.mink4, 99412))
        );
        let input = function!(callback, function!(callback, &word));
        calls.lock().unwrap().clear();
        assert_eq!(input.chainify(Bispinor {}.into()), input);
        assert!(calls.lock().unwrap().is_empty());
        let raw = gamma!(
            slot!(r.bis4, 99413),
            slot!(r.bis4, 99414),
            slot!(r.mink4, 99415)
        );
        let input = function!(callback, &raw);
        calls.lock().unwrap().clear();
        let result = input.chainify(Bispinor {}.into());
        // The one changed argument invokes its containing callback once.
        assert_eq!(calls.lock().unwrap().len(), 1);
        // Undo preserves Scalar-tagged wrappers as opaque. Check the changed
        // argument at its own tensor boundary instead of changing that policy.
        let AtomView::Fun(wrapped) = result.as_view() else {
            panic!("the scalar callback must retain its wrapper");
        };
        let body = wrapped.iter().next().unwrap();
        assert_eq!(body.undo_chain::<AbstractIndex>().unwrap(), raw);
        assert_eq!(*calls.lock().unwrap(), vec![result.clone()]);
    }

    #[test]
    fn chainify_changed_rounded_ancestors_match_the_original_map_primitive() {
        use spenso::shadowing::IntoAtom;
        use std::sync::{Arc, Mutex};
        let r = test_initialize();
        let calls = Arc::new(Mutex::new(Vec::new()));
        let recorded = Arc::clone(&calls);
        let callback = symbolica::symbol!("chainify_changed_callback"; Scalar;
            norm = move |node, _| recorded.lock().unwrap().push(node.to_owned())
        );
        let start = slot!(r.bis4, 99421).into_atom();
        let end = slot!(r.bis4, 99422).into_atom();
        let first = slot!(r.mink4, 99423).into_atom();
        let second = slot!(r.mink4, 99424).into_atom();
        let raw = [gamma!(&start, &end, &first), gamma!(&start, &end, &second)];
        let expected_words = [
            chain!(&start, &end, gamma!(&first)),
            chain!(&start, &end, gamma!(&second)),
        ];
        let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
        let body = (&rounded * &raw[0] + rounded.pow(2) * &raw[1]) * (&raw[0] + &rounded * &raw[1]);
        let input = function!(callback, body);
        calls.lock().unwrap().clear();
        let expected = input.replace_map(|node, _, output| {
            if let Some(position) = raw.iter().position(|value| value.as_view() == node) {
                **output = expected_words[position].clone();
            }
        });
        let expected_calls = std::mem::take(&mut *calls.lock().unwrap());
        assert_eq!(expected_calls.len(), 1);
        let actual = input.chainify(Bispinor {}.into());
        assert_eq!(actual, expected);
        assert_eq!(*calls.lock().unwrap(), expected_calls);
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
    fn join_chains_keeps_other_representations_factored() {
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
                .simplify_color(crate::color::ColorSimplifySettings {
                    simplify_non_color: false,
                    ..Default::default()
                })
                .unwrap()
                .resolved()
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
