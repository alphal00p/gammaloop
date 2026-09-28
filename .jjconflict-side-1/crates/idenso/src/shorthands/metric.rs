use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, Mutex},
};
#[cfg(test)]
use symbolica::atom::FunctionBuilder;

use crate::{
    IndexToolingError,
    shorthands::schoonschip::DotNormalizer,
    tensor::{
        SymbolicNet, SymbolicTensor,
        inference::{InterfaceInference, TensorInferenceError},
    },
};

#[cfg(test)]
use super::schoonschip::Schoonschip;
use super::schoonschip::{SchoonschipSettings, SchoonschipWithSettings};
use spenso::{
    network::{
        library::{DummyLibrary, symbolic::ETS},
        parsing::{ParamNet, ParseSettings, StrictTensorFilter},
        tags::SPENSO_TAG,
    },
    structure::{
        HasName, OrderedStructure, Reindexed, TensorStructure, ToSymbolic,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        partial::PartialStructure,
        representation::{ExtendibleReps, LibraryRep, LibrarySlot, RepName},
        slot::{AbsInd, DummyAind, IsAbstractSlot, ParseableAind, SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomType, AtomView, Symbol},
    function,
    id::{Condition, Match, MatchSettings, MatchStack, PatternRestriction, WildcardRestriction},
};
use symbolica_utils::{IntoArgs, IntoSymbol};

use eyre::Result;

use crate::rep_symbols::RS;

pub fn wrap_indices_impl(view: AtomView, header: Symbol) -> Atom {
    AbstractIndex::wrap_expression(view, header)
}

fn dangling_indices<Aind: ParseableAind + AbsInd + DummyAind>(
    view: AtomView,
) -> Result<Vec<Atom>, String> {
    let library = DummyLibrary::<SymbolicTensor<OrderedStructure<LibraryRep, Aind>>>::new();
    let slots = SymbolicNet::<Aind>::try_external_slots::<OrderedStructure<LibraryRep, Aind>, _>(
        view,
        &library,
        &ParseSettings {
            take_first_term_from_sum: true,
            strict_tensor_filter: StrictTensorFilter::ContainsReps,
            ..Default::default()
        },
    )
    .map_err(|error| error.to_string())?;
    Ok(slots.into_iter().map(|slot| slot.to_atom()).collect())
}

pub fn list_dangling_impl<Aind: ParseableAind + AbsInd + DummyAind>(
    view: AtomView,
) -> Result<Vec<Atom>, IndexToolingError> {
    dangling_indices::<Aind>(view).map_err(|reason| IndexToolingError::ListDangling { reason })
}

pub fn wrap_dummies_impl<Aind: ParseableAind + AbsInd + DummyAind>(
    view: AtomView,
    header: Symbol,
) -> Result<Atom, IndexToolingError> {
    let externals: HashSet<_> = dangling_indices::<Aind>(view)
        .map_err(|reason| IndexToolingError::WrapDummies { reason })?
        .into_iter()
        .collect();

    let mut expr = view.to_owned();
    let settings = MatchSettings::new().min_level(0).max_level(Some(0));

    for i in LibraryRep::all_self_duals().chain(LibraryRep::all_inline_metrics()) {
        let ipat = i.to_symbolic([RS.d_, RS.a_]).to_pattern();
        expr = expr.replace_map(|term, ctx, out| {
            if ctx.function_level < 2
                && ctx.function_level > 0
                && let Some(c) = term.pattern_match(&ipat, None, &settings).next()
            {
                let atom = ipat.replace_wildcards(&c).unwrap();
                if !externals.contains(&atom) {
                    **out =
                        i.to_symbolic([c[&RS.d_].clone(), function!(header, c[&RS.a_].clone())]);
                }
            }
        });
    }
    for i in LibraryRep::all_dualizables() {
        let ipat = i.to_symbolic([RS.d_, RS.a_]).to_pattern();
        let ipat_dual = i.dual().to_symbolic([RS.d_, RS.a_]).to_pattern();

        expr = expr.replace_map(|term, ctx, out| {
            if ctx.function_level < 2 && ctx.function_level > 0 {
                if let Some(c) = term.pattern_match(&ipat, None, &settings).next() {
                    let atom = ipat.replace_wildcards(&c).unwrap();
                    if !externals.contains(&atom) {
                        **out = i
                            .to_symbolic([c[&RS.d_].clone(), function!(header, c[&RS.a_].clone())]);
                    }
                } else if let Some(c) = term.pattern_match(&ipat_dual, None, &settings).next() {
                    let atom = ipat_dual.replace_wildcards(&c).unwrap();
                    if !externals.contains(&atom) {
                        **out = i
                            .dual()
                            .to_symbolic([c[&RS.d_].clone(), function!(header, c[&RS.a_].clone())]);
                    }
                }
            }
        });
    }

    Ok(expr)
}

pub(crate) fn not_slot(sym: Symbol) -> Condition<PatternRestriction> {
    sym.restrict(WildcardRestriction::IsAtomType(AtomType::Var))
        | sym.restrict(WildcardRestriction::IsAtomType(AtomType::Num))
        | sym.restrict(WildcardRestriction::filter(|a| match a {
            Match::FunctionName(f) => !f.has_tag(&SPENSO_TAG.representation),
            Match::Multiple(_, views) => {
                !views.iter().any(|a| a.has_attributes_of(SPENSO_TAG.rep_))
            }
            Match::Single(s) => !s.has_attributes_of(SPENSO_TAG.rep_),
        }))
}

// pub fn not_aind(sym: Symbol) -> Condition<PatternRestriction> {
//     sym.restrict(WildcardRestriction::IsAtomType(AtomType::Var))
//         | sym.restrict(WildcardRestriction::IsAtomType(AtomType::Num))
//         | sym.restrict(WildcardRestriction::filter(|a| match a {
//             Match::FunctionName(f) => {
//                 println!("FunctionName{f}");
//                 LibraryRep::all_representations().all(|r| r.symbol() != *f)
//             }
//             Match::Multiple(_, views) => {
//                 println!("Multiple:");
//                 for v in views {
//                     print!("{v}");
//                 }
//                 views
//                     .iter()
//                     .all(|a| LibrarySlot::<Parsind>::try_from(*a).is_err())
//             }
//             Match::Single(s) => {
//                 println!("Single{s}");
//                 LibrarySlot::<Parsind>::try_from(*s).is_err()
//             }
//         }))
// }

#[cfg(test)]
pub(crate) fn to_dots_impl(expr: AtomView) -> Atom {
    fn append_rep(atom: Atom, rep: &Atom) -> Atom {
        match atom.as_view() {
            AtomView::Fun(fun) => {
                let mut rebuilt = FunctionBuilder::new(fun.get_symbol());
                for arg in fun.iter() {
                    rebuilt = rebuilt.add_arg(arg);
                }
                rebuilt.add_arg(rep).finish()
            }
            AtomView::Var(var) => FunctionBuilder::new(var.get_symbol()).add_arg(rep).finish(),
            _ => atom,
        }
    }

    fn func_with_rep(m: &MatchStack<'_>, fun_wild: Symbol, arg_wild: Symbol, rep: &Atom) -> Atom {
        match m.get(arg_wild).unwrap() {
            Match::FunctionName(_) => {
                panic!("Can't be a function")
            }
            Match::Single(a) => {
                let Match::FunctionName(f) = m.get(fun_wild).unwrap() else {
                    panic!("Not a function");
                };
                FunctionBuilder::new(*f).add_arg(a).add_arg(rep).finish()
            }

            Match::Multiple(_, args) => {
                let atom = if args.is_empty() {
                    m.get(fun_wild).unwrap().to_atom()
                } else if let Match::FunctionName(f) = m.get(fun_wild).unwrap() {
                    FunctionBuilder::new(*f).add_args(args).finish()
                } else {
                    panic!("Not a function");
                };

                append_rep(atom, rep)
            }
        }
    }

    // A metric carrying a vector as one argument is a contraction, not a
    // rank-one function with that vector as a scalar label.
    let [vector_f, vector_g] = [RS.f_, RS.g_].map(|head| {
        head.restrict(WildcardRestriction::filter(
            |matched| !matches!(matched, Match::FunctionName(symbol) if *symbol == ETS.metric),
        ))
    });

    expr.replace(
        function!(
            RS.f_,
            RS.a___,
            SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_])
        ) * function!(
            RS.g_,
            RS.b___,
            SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_])
        ),
    )
    .min_level(0)
    .max_level(Some(0))
    .when(not_slot(RS.a___) & not_slot(RS.b___) & vector_f.clone() & vector_g.clone())
    .repeat()
    .with_map(move |m| {
        let rep = SPENSO_TAG
            .self_dual_::<0, _>([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);
        let f = func_with_rep(m, RS.f_, RS.a___, &rep);
        let g = func_with_rep(m, RS.g_, RS.b___, &rep);

        function!(SPENSO_TAG.dot, f, g)
    })
    .replace(
        function!(
            RS.f_,
            RS.a___,
            SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_])
        )
        .pow(2),
    )
    .min_level(0)
    .max_level(Some(0))
    .when(not_slot(RS.a___) & vector_f.clone())
    .repeat()
    .with_map(move |m| {
        let rep = SPENSO_TAG
            .self_dual_::<0, _>([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);
        let f = func_with_rep(m, RS.f_, RS.a___, &rep);

        function!(SPENSO_TAG.dot, &f, &f)
    })
    .replace(
        function!(
            RS.f_,
            RS.a___,
            SPENSO_TAG.dualizable_::<0, _>([RS.d_, RS.i_])
        ) * function!(
            RS.g_,
            RS.b___,
            SPENSO_TAG.dualizable_dual_::<0, _>([RS.d_, RS.i_])
        ),
    )
    .min_level(0)
    .max_level(Some(0))
    .when(not_slot(RS.a___) & not_slot(RS.b___) & vector_f.clone() & vector_g.clone())
    .repeat()
    .with_map(move |m| {
        let rep = SPENSO_TAG
            .dualizable_::<0, _>([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);
        let dual_rep = SPENSO_TAG
            .dualizable_dual_::<0, _>([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);
        let f = func_with_rep(m, RS.f_, RS.a___, &rep);
        let g = func_with_rep(m, RS.g_, RS.b___, &dual_rep);

        function!(SPENSO_TAG.dot, f, g)
    })
    .replace(function!(ETS.metric, RS.f_, RS.g_))
    .min_level(0)
    .max_level(Some(0))
    .when(not_slot(RS.f_) & not_slot(RS.g_))
    .repeat()
    .with(function!(SPENSO_TAG.dot, RS.f_, RS.g_))
}

fn matched_i64(m: &MatchStack<'_>, symbol: Symbol) -> Option<i64> {
    match m.get(symbol)? {
        Match::Single(value) => i64::try_from(*value).ok(),
        _ => None,
    }
}

fn simplify_generated_metric_components(expr: Atom) -> Atom {
    let pat = function!(ETS.metric, function!(AIND_SYMBOLS.cind, RS.i_, RS.j_)).to_pattern();
    let result = expr.replace(pat.clone()).with_map(move |m| {
        let Some(i) = matched_i64(m, RS.i_) else {
            return pat.replace_wildcards_with_matches(m);
        };
        let Some(j) = matched_i64(m, RS.j_) else {
            return pat.replace_wildcards_with_matches(m);
        };

        if i != j {
            Atom::Zero
        } else if i == 0 {
            Atom::num(1)
        } else {
            Atom::num(-1)
        }
    });
    SchoonschipWithSettings {
        settings: &SchoonschipSettings::default().without_rank1_tensors(),
    }
    .run(result.as_view(), &mut Vec::new())
}

impl SymbolicTensor<PartialStructure> {
    /// Evaluate finite components of compact dots, including dots in opaque
    /// scalar payloads. Symbolic dimensions stay unchanged. This is distinct
    /// from `undo_dots`, which only opens symbolic indices, and from `contract`.
    /// Neither the surrounding expression nor its numerator is expanded.
    /// The result is a raw component expression, not a symbolic tensor proof:
    /// concretized heads may carry `cind` arguments instead of symbolic ports.
    pub fn expand_dots(&self) -> Result<Atom, TensorInferenceError> {
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let compact = DotNormalizer::metric_shorthand_to_dot(self.expression().as_view());
        let pat = function!(SPENSO_TAG.dot, RS.f_, RS.g_).to_pattern();
        let metric_pat = function!(ETS.metric, RS.f_, RS.g_).to_pattern();
        // Only intrinsically normalized, callback-free component operations may
        // run ahead of their replacement position. Custom metric sign functions
        // are callbacks too, even when all tensor heads have no normalizer.
        let mut dots = Vec::new();
        let mut metrics = Vec::new();
        let mut seen = HashSet::new();
        let mut matches = compact.pattern_match(&pat, None, None);
        while let Some(matched) = matches.next_detailed() {
            let dot = pat.replace_wildcards_with_matches(matched.match_stack);
            if Self::component_batch_is_intrinsic(dot.as_view()) && seen.insert(dot.clone()) {
                dots.push(dot);
                metrics.push(metric_pat.replace_wildcards_with_matches(matched.match_stack));
            }
        }
        drop(matches);
        let inputs = metrics.iter().map(Atom::as_view).collect::<Vec<_>>();
        let results = ParamNet::<AbstractIndex>::evaluate_scalar_batch(&inputs, &settings)
            .map_err(|error| TensorInferenceError::Invalid(error.to_string()))?;
        let prepared: HashMap<_, _> = dots
            .into_iter()
            .zip(results)
            .map(|(dot, result)| {
                let value =
                    result.map_or_else(|_| dot.clone(), simplify_generated_metric_components);
                (dot, value)
            })
            .collect();
        let failure = Arc::new(Mutex::new(None));
        let observed_failure = Arc::clone(&failure);
        let expression = compact.as_view().replace(pat.clone()).with_map(move |a| {
            let filled = pat.replace_wildcards_with_matches(a);
            if let Some(result) = prepared.get(&filled) {
                return result.clone();
            }
            if observed_failure.lock().unwrap().is_some() {
                return filled;
            }
            // The same component owner handles sensitive inputs one at a time
            // here, preserving parse -> execute -> next match and enclosing
            // normalizer order. No speculative callback or result cache applies.
            let metric_filled = metric_pat.replace_wildcards_with_matches(a);
            match ParamNet::<AbstractIndex>::evaluate_scalar_batch(
                &[metric_filled.as_view()],
                &settings,
            ) {
                Ok(mut values) => values
                    .pop()
                    .expect("one input has one result")
                    .map_or(filled, simplify_generated_metric_components),
                Err(error) => {
                    *observed_failure.lock().unwrap() = Some(error.to_string());
                    filled
                }
            }
        });
        if let Some(error) = failure.lock().unwrap().take() {
            return Err(TensorInferenceError::Invalid(error));
        }
        Ok(expression)
    }

    fn component_batch_is_intrinsic(value: AtomView<'_>) -> bool {
        if !InterfaceInference::normalization_is_intrinsic(value) {
            return false;
        }
        let mut intrinsic = true;
        let mut slots = SlotMatcher::default();
        value.visitor(&mut |node| {
            if let AtomView::Fun(function) = node {
                intrinsic &= !function.get_symbol().has_tag(&SPENSO_TAG.broadcast);
                let representation = match slots.classify(node) {
                    SlotMatch::Explicit(slot) => slots.representation(slot).ok(),
                    _ => slots
                        .parse_representation::<LibraryRep>(node)
                        .ok()
                        .map(|rep| rep.rep),
                };
                if let Some(LibraryRep::InlineMetric(_)) = representation {
                    intrinsic &= representation == Some(ExtendibleReps::MINKOWSKI);
                }
            }
            intrinsic
        });
        intrinsic
    }
}

pub trait PermuteWithMetric {
    fn permute_with_metric(self) -> Atom;
}

impl<N, Aind: AbsInd + DummyAind + ParseableAind> PermuteWithMetric for Reindexed<N>
where
    N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
    N::Name: IntoSymbol + Clone,
    N::Args: IntoArgs,
{
    fn permute_with_metric(self) -> Atom {
        let value = self
            .map_target(|a| SymbolicTensor::from_named(&a).unwrap())
            .apply();
        SchoonschipWithSettings {
            settings: &SchoonschipSettings::default().without_rank1_tensors(),
        }
        .run(value.expression.as_view(), &mut Vec::new())
    }
}

#[cfg(test)]
mod test {

    use crate::{Cookable, representations::Bispinor, test_support::test_initialize};
    use symbolica_utils::AtomPrintExt;

    use super::*;

    use spenso::{
        network::parsing::{ShadowedStructure, StructureFromAtom},
        structure::{
            IndexlessNamedStructure,
            abstract_index::AbstractIndex,
            representation::{Euclidean, Lorentz},
        },
    };
    use symbolica::{parse_lit, symbol};

    #[test]
    fn dangling_observation_preserves_shorthands_branches_and_errors() {
        use crate::tensor::SymbolicNetParse;
        test_initialize();
        let p = spenso::vector_symbol!("dangling_observation::p");
        let q = spenso::vector_symbol!("dangling_observation::q");
        let rep = parse_lit!(spenso::mink(4));
        let p1 = function!(p, parse_lit!(spenso::mink(4, 1)));
        let q1 = function!(q, parse_lit!(spenso::mink(4, 1)));
        let p2 = function!(p, parse_lit!(spenso::mink(4, 2)));
        let compact_p = function!(p, &rep);
        let compact_q = function!(q, &rep);
        let settings = ParseSettings {
            take_first_term_from_sum: true,
            ..Default::default()
        };
        let chain = symbolica::parse!(
            "spenso::chain(spenso::bis(4,3),spenso::bis(4,7),spenso::gamma(spenso::in,spenso::out,spenso::mink(4,5)))"
        );
        let trace = symbolica::parse!(
            "spenso::trace(spenso::bis(4),spenso::gamma(spenso::in,spenso::out,spenso::mink(4,5)),spenso::gamma(spenso::in,spenso::out,spenso::mink(4,9)))"
        );
        for input in [
            &p1 * &q1 * &p2,
            (&p1 + &q1) * &p2,
            &p1 + &p2,
            p1.clone().pow(2),
            p1.clone().pow(3),
            spenso::bracket!(&compact_q, &compact_p),
            function!(SPENSO_TAG.dot, &compact_p, &compact_q),
            chain,
            trace,
        ] {
            let expected = input
                .parse_to_symbolic_net::<AbstractIndex>(&settings)
                .unwrap()
                .graph
                .dangling_indices();
            assert_eq!(
                dangling_indices::<AbstractIndex>(input.as_view()).unwrap(),
                expected
                    .iter()
                    .map(IsAbstractSlot::to_atom)
                    .collect::<Vec<_>>(),
                "{input}"
            );
        }
        let library = DummyLibrary::<SymbolicTensor<OrderedStructure>>::new();
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        assert!(matches!(
            SymbolicNet::<AbstractIndex>::try_external_slots::<OrderedStructure, _>(
                (&p1 + &p2).as_view(),
                &library,
                &settings
            ),
            Err(spenso::network::TensorNetworkError::IncompatibleSummand(_))
        ));
        let invalid = p1.pow(-1);
        assert!(matches!(
            invalid.parse_to_symbolic_net::<AbstractIndex>(&settings),
            Err(spenso::network::TensorNetworkError::NegativeExponentNonScalar(_))
        ));
        assert!(matches!(
            SymbolicNet::<AbstractIndex>::try_external_slots::<OrderedStructure, _>(
                invalid.as_view(),
                &library,
                &settings
            ),
            Err(spenso::network::TensorNetworkError::NegativeExponentNonScalar(_))
        ));
    }

    #[test]
    fn dangling_observation_preserves_materialization_callback_order() {
        use crate::tensor::SymbolicNetParse;
        test_initialize();
        let events = Arc::new(Mutex::new(Vec::new()));
        let observed = Arc::clone(&events);
        let p = spenso::vector_symbol!(
            "dangling_observation_callback::p",
            norm = move |value, _output| {
                observed.lock().unwrap().push(value.to_owned());
            }
        );
        let q = spenso::vector_symbol!("dangling_observation_callback::q");
        let rep = parse_lit!(spenso::mink(4));
        let input = function!(SPENSO_TAG.dot, function!(p, &rep), function!(q, &rep));
        events.lock().unwrap().clear();
        let expected = input
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                take_first_term_from_sum: true,
                ..Default::default()
            })
            .unwrap()
            .graph
            .dangling_indices();
        let expected_events = std::mem::take(&mut *events.lock().unwrap());
        let actual = dangling_indices::<AbstractIndex>(input.as_view()).unwrap();
        assert_eq!(
            actual,
            expected
                .iter()
                .map(IsAbstractSlot::to_atom)
                .collect::<Vec<_>>()
        );
        assert_eq!(*events.lock().unwrap(), expected_events);
        assert!(!expected_events.is_empty());
    }

    #[test]
    fn dangling_observation_preserves_custom_index_grammar() {
        use crate::tensor::SymbolicNetParse;
        use spenso::structure::slot::SlotError;
        #[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
        struct PayloadIndex(usize);
        impl std::fmt::Display for PayloadIndex {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "payload({})", self.0)
            }
        }
        impl AbsInd for PayloadIndex {}
        impl ParseableAind for PayloadIndex {
            type Error = SlotError;
            fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
                if let AtomView::Fun(fun) = view
                    && fun.get_symbol() == symbol!("dangling_payload")
                    && fun.get_nargs() == 1
                {
                    return usize::try_from(fun.iter().next().unwrap())
                        .map(Self)
                        .map_err(|_| SlotError::NotNatural);
                }
                Err(SlotError::Composite)
            }
            fn to_atom(&self) -> Atom {
                function!(symbol!("dangling_payload"), self.0)
            }
        }
        impl DummyAind for PayloadIndex {
            fn new_dummy() -> Self {
                Self(usize::MAX)
            }
            fn new_dummy_at(i: usize) -> Self {
                Self(i)
            }
            fn is_dummy(&self) -> bool {
                self.0 == usize::MAX
            }
        }
        test_initialize();
        let p = spenso::vector_symbol!("dangling_custom::p");
        let q = spenso::vector_symbol!("dangling_custom::q");
        let slot = |i| {
            spenso::structure::representation::Minkowski {}
                .new_rep(4)
                .to_symbolic([PayloadIndex(i).to_atom()])
        };
        let input = (function!(p, slot(7)) + function!(q, slot(7))) * function!(p, slot(9));
        let expected = input
            .parse_to_symbolic_net::<PayloadIndex>(&ParseSettings {
                take_first_term_from_sum: true,
                ..Default::default()
            })
            .unwrap()
            .graph
            .dangling_indices();
        assert_eq!(
            dangling_indices::<PayloadIndex>(input.as_view()).unwrap(),
            expected
                .iter()
                .map(IsAbstractSlot::to_atom)
                .collect::<Vec<_>>()
        );
        assert!(dangling_indices::<AbstractIndex>(input.as_view()).is_err());
    }

    // The old finite-component route is retained only as an independent
    // scheduling oracle. It uses one ordinary network per original dot.
    fn independent_components(input: &Atom) -> Atom {
        use spenso::network::parsing::NetworkParse;
        let settings =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let compact = input.to_dots();
        let pattern = function!(SPENSO_TAG.dot, RS.f_, RS.g_).to_pattern();
        let metric = function!(ETS.metric, RS.f_, RS.g_).to_pattern();
        compact.replace(pattern.clone()).with_map(move |matched| {
            let original = pattern.replace_wildcards_with_matches(matched);
            let input = metric.replace_wildcards_with_matches(matched);
            match input.parse_to_atom_net::<AbstractIndex>(&settings) {
                Err(_) => original,
                Ok(mut network) => {
                    network.simple_execute();
                    network.result_scalar().map_or(original, |scalar| {
                        simplify_generated_metric_components(Atom::from(scalar))
                    })
                }
            }
        })
    }

    #[test]
    fn finite_dots_keep_unsupported_dimensions_and_opaque_scalar_payloads() {
        test_initialize();
        let p = spenso::vector_symbol!("finite_dot_scope::p");
        let q = spenso::vector_symbol!("finite_dot_scope::q");
        let dot = |rep: Atom| function!(SPENSO_TAG.dot, function!(p, &rep), function!(q, &rep));
        let finite = dot(parse_lit!(spenso::mink(2)));
        let symbolic = dot(parse_lit!(spenso::mink(finite_dot_scope::D)));
        let wrapper = symbol!("finite_dot_scope::opaque");
        let source = function!(wrapper, &finite, &symbolic)
            * (parse_lit!(finite_dot_scope::x) + parse_lit!(finite_dot_scope::y));
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        let result = value.expand_dots().unwrap();
        assert_eq!(result, independent_components(&source));
        let finite_result = independent_components(&finite);
        assert_eq!(
            result,
            (function!(wrapper, &finite_result, &symbolic)
                * (parse_lit!(finite_dot_scope::x) + parse_lit!(finite_dot_scope::y)))
        );
        assert_eq!(value.expand_dots().unwrap(), result);
        assert_eq!(value.expression(), &source);
    }

    #[test]
    fn finite_dots_preserve_mixed_callback_and_outer_normalizer_order() {
        test_initialize();
        let events = Arc::new(Mutex::new(Vec::new()));
        let matrix = spenso::tensor_symbol!("finite_dot_schedule::M");
        let vector = spenso::vector_symbol!("finite_dot_schedule::V");
        let other = spenso::vector_symbol!("finite_dot_schedule::B");
        let plain = spenso::vector_symbol!("finite_dot_schedule::plain");
        let observed = Arc::clone(&events);
        let callback = spenso::vector_symbol!(
            "finite_dot_schedule::P",
            norm = move |value, output| {
                let AtomView::Fun(function) = value else {
                    return;
                };
                let argument = function.iter().next().unwrap();
                observed
                    .lock()
                    .unwrap()
                    .push(("vector", argument.to_owned()));
                if SlotMatcher::default()
                    .parse::<LibraryRep, AbstractIndex>(argument)
                    .is_ok()
                {
                    let internal = parse_lit!(spenso::mink(2, 77));
                    **output =
                        function!(matrix, argument, &internal) * function!(vector, &internal);
                }
            }
        );
        let observed = Arc::clone(&events);
        let outer = symbol!(
            "finite_dot_schedule::outer",
            norm = move |value, _output| {
                observed.lock().unwrap().push(("outer", value.to_owned()));
            }
        );
        let rep = parse_lit!(spenso::mink(2));
        let dot = |head, label| {
            function!(
                SPENSO_TAG.dot,
                function!(head, &rep),
                function!(other, label, &rep)
            )
        };
        let source = function!(outer, dot(callback, 1), dot(plain, 2), dot(callback, 3));
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        events.lock().unwrap().clear();
        let expected = independent_components(&source);
        let expected_events = std::mem::take(&mut *events.lock().unwrap());
        let result = value.expand_dots().unwrap();
        assert_eq!(result, expected);
        assert_eq!(*events.lock().unwrap(), expected_events);
        assert!(expected_events.iter().any(|(kind, _)| *kind == "outer"));
        let concrete_calls = expected_events
            .iter()
            .filter(|(kind, argument)| {
                *kind == "vector"
                    && SlotMatcher::default()
                        .parse::<LibraryRep, AbstractIndex>(argument.as_view())
                        .is_ok()
            })
            .map(|(_, argument)| argument)
            .collect::<Vec<_>>();
        assert_eq!(concrete_calls.len(), 2);
        assert_eq!(concrete_calls[0], concrete_calls[1]);
        assert!(!result.to_string().contains("mink"));
        assert!(!SymbolicTensor::component_batch_is_intrinsic(
            dot(callback, 1).as_view()
        ));
        assert!(SymbolicTensor::component_batch_is_intrinsic(
            dot(plain, 2).as_view()
        ));
    }

    #[test]
    fn finite_dots_keep_ordered_open_interfaces_and_typed_zero() {
        test_initialize();
        let p = spenso::vector_symbol!("finite_dot_ports::p");
        let q = spenso::vector_symbol!("finite_dot_ports::q");
        let a = spenso::vector_symbol!("finite_dot_ports::a");
        let b = spenso::vector_symbol!("finite_dot_ports::b");
        let rep = parse_lit!(spenso::mink(2));
        let scalar = function!(SPENSO_TAG.dot, function!(p, &rep), function!(q, &rep));
        let ordered = spenso::bracket!(function!(b, &rep), function!(a, &rep));
        let value = SymbolicTensor::infer(&scalar * &ordered).unwrap();
        let result = value.expand_dots().unwrap();
        assert_eq!(result, independent_components(&scalar) * &ordered);
        assert_eq!(value.expression(), &(&scalar * &ordered));
        let zero = SymbolicTensor::from_normalized_parts(Atom::Zero, value.structure().clone());
        let result = zero.expand_dots().unwrap();
        assert_eq!(result, Atom::Zero);
        assert_eq!(zero.structure(), value.structure());
    }

    #[test]
    fn finite_dots_keep_raw_callback_output_at_materialization_boundary() {
        test_initialize();
        let p = spenso::vector_symbol!("finite_dot_rank::p");
        let q = spenso::vector_symbol!("finite_dot_rank::q");
        let head = spenso::tensor_symbol!(
            "finite_dot_rank::T",
            norm = |value, output| {
                let AtomView::Fun(function) = value else {
                    return;
                };
                let first = function.iter().next().unwrap();
                if !first.contains_symbol(SPENSO_TAG.dot) {
                    **output = Atom::num(1);
                }
            }
        );
        let rep = parse_lit!(spenso::mink(2));
        let dot = function!(SPENSO_TAG.dot, function!(p, &rep), function!(q, &rep));
        let source = function!(head, &dot, parse_lit!(spenso::mink(2, finite_dot_rank::i)));
        let value = SymbolicTensor::infer(source.clone()).unwrap();
        assert_eq!(
            value.expand_dots().unwrap(),
            independent_components(&source)
        );
        assert_eq!(value.expand_dots().unwrap(), Atom::num(1));
        assert_eq!(value.expression(), &source);
        assert_eq!(value.structure().canonical().order(), 1);
    }

    #[test]
    fn cook() {
        test_initialize();
        let expr = parse_lit!(
            spenso::g(spenso::mink(4, f(0)), spenso::dind(spenso::cof(4, f(1))))
                * p(spenso::mink(4, 1))
        )
        .cook_indices();

        println!("{}", expr);
    }

    #[test]
    fn metric_contract() {
        test_initialize();
        let expr =
            parse_lit!(spenso::g(spenso::mink(4, 0), spenso::mink(4, 1)) * p(spenso::mink(4, 1)))
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());

        assert_eq!(expr, parse_lit!(p(spenso::mink(4, 0))), "got {:#}", expr);
    }

    #[test]
    fn permute() {
        test_initialize();
        let f = IndexlessNamedStructure::<Symbol, ()>::from_iter(
            [
                Lorentz {}.new_rep(8).to_lib(),
                Lorentz {}.new_rep(2).cast(),
                Lorentz {}.new_rep(2).cast(),
                Euclidean {}.new_rep(2).cast(),
                Lorentz {}.new_rep(4).cast(),
                Lorentz {}.new_rep(2).cast(),
                Lorentz {}.new_rep(7).cast(),
            ],
            symbol!("test"),
            None,
        );
        let logical_indices: [AbstractIndex; 7] = [
            6.into(),
            4.into(),
            5.into(),
            2.into(),
            3.into(),
            1.into(),
            0.into(),
        ];
        let storage_indices = f.layout().logical_to_canonical(&logical_indices);
        let layout = f.layout().clone();
        let order = f.canonical().order();
        let f = f
            .into_canonical()
            .reindex_storage(&storage_indices)
            .unwrap()
            .map_target(|a| SymbolicTensor::from_named(&a).unwrap());

        let f_p = f.apply();

        let simplified = f_p
            .expression
            .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());
        let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

        assert_eq!(order, f_parsed.canonical().order());
        let logical_positions = (0..order).collect::<Vec<_>>();
        let canonical_positions = layout.logical_to_canonical(&logical_positions);
        assert_eq!(
            layout.canonical_to_logical(&canonical_positions),
            logical_positions
        );
    }

    #[test]
    fn id_trace() {
        test_initialize();
        let bis = Bispinor {}.new_rep(symbol!("dim"));

        let expr = bis
            .g(9, 9)
            .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors());

        assert_eq!(expr, Atom::var(symbol!("dim")), "got {:#}", expr);
    }
    #[test]
    fn dotsx() {
        test_initialize();

        let a = parse_lit!(P(label, spenso::mink(4, 2)) * P(spenso::mink(4, 2))).to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(P(label,mink(4)),P(mink(4)))");
        insta::assert_snapshot!(SymbolicTensor::infer(a).unwrap().expand_dots().unwrap().to_bare_ordered_string(),@"-1*P(cind(1))*P(label,cind(1))+-1*P(cind(2))*P(label,cind(2))+-1*P(cind(3))*P(label,cind(3))+P(cind(0))*P(label,cind(0))");

        let a = parse_lit!(Q(spenso::mink(4, mu1)) ^ 2).to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(Q(mink(4)),Q(mink(4)))");

        let a = parse_lit!(spenso::g(
            K(1, spenso::mink(4)) + K(2, spenso::mink(4)),
            P(3, spenso::mink(4))
        ))
        .to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(K(1,mink(4)),P(3,mink(4)))+dot(K(2,mink(4)),P(3,mink(4)))");
    }

    #[test]
    fn true_cooking() {
        test_initialize();
        let expr = parse_lit!(spenso::g(spenso::mink(4, true_cooking_index(0))));

        assert_eq!(expr.cook().uncook(), expr);
    }
}
