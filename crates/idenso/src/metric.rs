use std::{collections::HashSet, sync::LazyLock};

use crate::{
    schoonschip::Schoonschip,
    tensor::{SymbolicNetParse, SymbolicTensor},
};
use spenso::{
    network::{
        Network,
        library::{DummyLibrary, symbolic::ETS},
        parsing::{NetworkParse, ParseSettings, StrictTensorFilter},
        store::NetworkStore,
        tags::SPENSO_TAG,
    },
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    structure::{
        HasName, PermutedStructure, TensorStructure, ToSymbolic,
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        permuted::Perm,
        representation::{LibraryRep, LibrarySlot, RepName},
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
    symbolica_atom::TensorCollectExt,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomType, AtomView, FunctionBuilder, Symbol},
    function,
    id::{
        Condition, FilterFn, Match, MatchSettings, MatchStack, PatternRestriction, Replacement,
        WildcardRestriction,
    },
    symbol,
};

use eyre::Result;

use super::rep_symbols::RS;

pub struct MetricSymbols {
    pub dim: Symbol,
    pub dummy: Symbol,
}

pub static MS: LazyLock<MetricSymbols> = LazyLock::new(|| MetricSymbols {
    dim: symbol!("spenso::dim"),
    dummy: symbol!("spenso::dummy"),
});

pub fn canonize_impl(view: AtomView) -> Atom {
    let lib = DummyLibrary::<SymbolicTensor>::new();
    let settings =
        ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
    let mut net = Network::<NetworkStore<SymbolicTensor, Atom>, _, Symbol>::try_from_view::<
        SymbolicTensor,
        _,
    >(view, &lib, &settings)
    .unwrap();

    let mut redual_reps = vec![];

    let mut indices = vec![];

    for t in net.store.tensors.iter_mut() {
        let mut reps = vec![];
        let mut pat = FunctionBuilder::new(t.name().unwrap());
        let mut rhs = pat.clone();
        for s in t.structure.external_structure_iter() {
            if !s.rep_name().is_self_dual() && s.rep_name().is_dual() {
                pat = pat.add_arg(s.rep().dual().to_symbolic([Atom::var(RS.a_)]));
                indices.push((s.dual().to_atom(), s.dual().rep()));
                reps.push(Replacement::new(
                    s.rep().to_symbolic([Atom::var(RS.a_)]).to_pattern(),
                    s.rep().dual().to_symbolic([Atom::var(RS.a_)]),
                ));
            } else {
                pat = pat.add_arg(s.rep().to_symbolic([Atom::var(RS.a_)]));
                indices.push((s.to_atom(), s.rep()));
            }
            rhs = rhs.add_arg(s.rep().to_symbolic([Atom::var(RS.a_)]));
        }
        if !reps.is_empty() {
            redual_reps.push(Replacement::new(pat.finish().to_pattern(), rhs.finish()));
            t.expression = t.expression.replace_multiple(&reps);
        }
    }

    view.canonize_tensors(indices)
        .unwrap()
        .canonical_form
        .replace_multiple(&redual_reps)
}

pub fn not_wraped_aind(header: Symbol) -> impl FilterFn + 'static {
    move |a| match a {
        Match::FunctionName(f) => *f != header,
        Match::Single(a) => {
            if let Some(n) = a.get_symbol() {
                n != header
            } else {
                true
            }
        }
        _ => false,
    }
}

pub fn wrap_indices_impl(view: AtomView, header: Symbol) -> Atom {
    let mut expr = view.collect_tensors();
    let dim = RS.d_;
    let dima = Atom::var(dim);
    let settings = MatchSettings {
        level_range: (0, Some(1)),
        ..Default::default()
    };

    let mut reps = vec![];
    for i in LibraryRep::all_self_duals().chain(LibraryRep::all_inline_metrics()) {
        reps.push(
            Replacement::new(
                i.to_symbolic([dim, RS.a_]).to_pattern(),
                i.to_symbolic([dima.clone(), function!(header, Atom::var(RS.a_))]),
            )
            .with_conditions(RS.a_.filter_match(not_wraped_aind(header)))
            .with_settings(settings.clone()),
        );
    }

    for i in LibraryRep::all_dualizables() {
        let di = i.dual();
        reps.push(
            Replacement::new(
                i.to_symbolic([dim, RS.a_]).to_pattern(),
                i.to_symbolic([dima.clone(), function!(header, Atom::var(RS.a_))]),
            )
            .with_conditions(RS.a_.filter_match(not_wraped_aind(header)))
            .with_settings(settings.clone()),
        );
        reps.push(
            Replacement::new(
                di.to_symbolic([dim, RS.a_]).to_pattern(),
                di.to_symbolic([dima.clone(), function!(header, Atom::var(RS.a_))]),
            )
            .with_conditions(RS.a_.filter_match(not_wraped_aind(header)))
            .with_settings(settings.clone()),
        );
    }
    let mut atom = Atom::new();
    while expr.replace_multiple_into(&reps, &mut atom) {
        std::mem::swap(&mut expr, &mut atom);
    }
    expr
}

pub fn list_dangling_impl<Aind: ParseableAind + AbsInd + DummyAind>(view: AtomView) -> Vec<Atom> {
    let a = view
        .parse_to_symbolic_net::<Aind>(&ParseSettings {
            take_first_term_from_sum: true,
            ..Default::default()
        })
        .unwrap();

    a.graph
        .dangling_indices()
        .into_iter()
        .map(|a| a.to_atom())
        .collect()
}

pub fn wrap_dummies_impl<Aind: ParseableAind + AbsInd + DummyAind>(
    view: AtomView,
    header: Symbol,
) -> Atom {
    let externals: HashSet<_> = list_dangling_impl::<Aind>(view).into_iter().collect();

    let mut expr = view.to_owned();
    let settings = MatchSettings {
        level_range: (0, Some(0)),
        ..Default::default()
    };

    for i in LibraryRep::all_self_duals().chain(LibraryRep::all_inline_metrics()) {
        let ipat = i.to_symbolic([RS.d_, RS.a_]).to_pattern();
        expr = expr.replace_map(|term, ctx, out| {
            if ctx.function_level < 2
                && ctx.function_level > 0
                && let Some(c) = term.pattern_match(&ipat, None, &settings).next()
            {
                let atom = ipat.replace_wildcards(&c);
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
                    let atom = ipat.replace_wildcards(&c);
                    if !externals.contains(&atom) {
                        **out = i
                            .to_symbolic([c[&RS.d_].clone(), function!(header, c[&RS.a_].clone())]);
                    }
                } else if let Some(c) = term.pattern_match(&ipat_dual, None, &settings).next() {
                    let atom = ipat_dual.replace_wildcards(&c);
                    if !externals.contains(&atom) {
                        **out = i
                            .dual()
                            .to_symbolic([c[&RS.d_].clone(), function!(header, c[&RS.a_].clone())]);
                    }
                }
            }
        });
    }

    expr
}

pub fn not_slot(sym: Symbol) -> Condition<PatternRestriction> {
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

pub fn to_dots_impl(expr: AtomView) -> Atom {
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
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___) & not_slot(RS.b___))
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
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___))
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
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___) & not_slot(RS.b___))
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
    .level_range((0, Some(0)))
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
    expr.replace(pat.clone()).with_map(move |m| {
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
    })
}

/// Trait for simplifying expressions involving metric tensors and converting
/// index contractions to compact Schoonschip scalar products.
///
/// Provides methods for contracting indices with metric tensors (`g(mu, nu)`) or
/// identity tensors (`id(mu, nu)`), and for replacing contracted index patterns
/// (like `p(mu)*q(mu)`) with `dot(p(rep), q(rep))`.
pub trait MetricSimplifier {
    /// Simplifies contractions involving metric tensors (`g` or `metric`) and identity tensors (`id` or `𝟙`).
    ///
    /// Applies rules like `g(mu, nu) * p(nu) -> p(mu)`, `g(mu, mu) -> D`, etc.
    ///
    /// # Returns
    /// An [`Atom`] representing the expression after metric simplification.
    fn simplify_metrics(&self) -> Atom;

    fn expand_dots(&self) -> Result<Atom>;

    /// Converts contracted index patterns into compact scalar-product notation.
    ///
    /// Replaces expressions like `p(mu) * q(mu)` or `p(mu) * M(mu, nu) * q(nu)` (implicitly via metric rules)
    /// with `dot(p(rep), q(rep))`. Assumes standard representations for vectors and tensors involved
    /// in the contractions.
    ///
    /// # Returns
    /// An [`Atom`] where contractions have been replaced by compact scalar products where possible.
    fn to_dots(&self) -> Atom;
}

impl MetricSimplifier for Atom {
    fn to_dots(&self) -> Atom {
        self.as_view().to_dots()
    }

    fn expand_dots(&self) -> Result<Atom> {
        self.as_view().expand_dots()
    }
    fn simplify_metrics(&self) -> Atom {
        self.schoonschip()
    }
}

impl MetricSimplifier for AtomView<'_> {
    fn to_dots(&self) -> Atom {
        to_dots_impl(*self)
    }

    fn expand_dots(&self) -> Result<Atom> {
        let set =
            ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps);
        let compact = self.to_dots();
        let pat = function!(SPENSO_TAG.dot, RS.f_, RS.g_).to_pattern();
        let metric_pat = function!(ETS.metric, RS.f_, RS.g_).to_pattern();
        Ok(compact.as_view().replace(pat.clone()).with_map(move |a| {
            let filled = pat.replace_wildcards_with_matches(a);
            let metric_filled = metric_pat.replace_wildcards_with_matches(a);

            match metric_filled.parse_to_atom_net::<AbstractIndex>(&set) {
                Err(_) => filled,
                Ok(mut net) => {
                    net.simple_execute();
                    match net.result_scalar() {
                        Ok(scalar) => simplify_generated_metric_components(Atom::from(scalar))
                            .simplify_metrics(),
                        Err(_) => filled,
                    }
                }
            }
        }))
    }
    fn simplify_metrics(&self) -> Atom {
        self.schoonschip()
    }
}

pub trait PermuteWithMetric {
    fn permute_with_metric(self) -> Atom;
}

impl<N, Aind: AbsInd + DummyAind + ParseableAind> PermuteWithMetric for PermutedStructure<N>
where
    N: ToSymbolic + HasName + TensorStructure<Slot = LibrarySlot<Aind>>,
    N::Name: IntoSymbol + Clone,
    N::Args: IntoArgs,
{
    fn permute_with_metric(self) -> Atom {
        self.map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
            .permute_inds()
            .expression
            .simplify_metrics()
    }
}

#[cfg(test)]
mod test {

    use crate::{Cookable, representations::Bispinor, test_support::test_initialize};

    use super::*;

    use spenso::{
        network::parsing::{ShadowedStructure, StructureFromAtom},
        shadowing::symbolica_utils::AtomCoreExt,
        structure::{
            IndexlessNamedStructure,
            abstract_index::AbstractIndex,
            permuted::Perm,
            representation::{Euclidean, Lorentz},
        },
    };
    use symbolica::parse_lit;

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
                .simplify_metrics();

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
        )
        .clone()
        .reindex([6, 4, 5, 2, 3, 1, 0])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap());

        let f_p = f.clone().permute_inds();

        let simplified = f_p.expression.simplify_metrics();
        let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

        assert_eq!(
            f.structure.structure.order(),
            f_parsed.structure.structure.order()
        );
        assert!(f_parsed.rep_permutation.is_identity());

        let f_p = f.clone().permute_reps_wrapped().permute_inds();

        let simplified = f_p.expression.simplify_metrics();
        let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

        assert_eq!(
            f.structure.structure.order(),
            f_parsed.structure.structure.order()
        );
    }

    #[test]
    fn id_trace() {
        test_initialize();
        let bis = Bispinor {}.new_rep(symbol!("dim"));

        let expr = bis.g(9, 9).simplify_metrics();

        assert_eq!(expr, Atom::var(symbol!("dim")), "got {:#}", expr);
    }
    #[test]
    fn dots() {
        test_initialize();

        let a =
            parse_lit!(P(spenso::mink(4, -1 * g(2)), spenso::mink(4, 2)) * P(spenso::mink(4, 2)))
                .to_dots();
        assert_eq!(
            a,
            parse_lit!(P(spenso::mink(4, 2)) * P(spenso::mink(4, -g(2)), spenso::mink(4, 2)))
        );
        let a = parse_lit!(P(label, spenso::mink(4, 2)) * P(spenso::mink(4, 2))).to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(P(label,mink(4)),P(mink(4)))");
        insta::assert_snapshot!(a.expand_dots().unwrap().to_bare_ordered_string(),@"-1*P(cind(1))*P(label,cind(1))+-1*P(cind(2))*P(label,cind(2))+-1*P(cind(3))*P(label,cind(3))+P(cind(0))*P(label,cind(0))");

        let a = parse_lit!(k1(spenso::mink(4, mu1)) ^ 2).to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(k1(mink(4)),k1(mink(4)))");

        let a = parse_lit!(spenso::g(
            k(1, spenso::mink(4)) + k(2, spenso::mink(4)),
            p(3, spenso::mink(4))
        ))
        .to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"dot(k(1,mink(4)),p(3,mink(4)))+dot(k(2,mink(4)),p(3,mink(4)))");
    }

    #[test]
    fn true_cooking() {
        test_initialize();
        let expr = parse_lit!(spenso::g(spenso::mink(4, true_cooking_index(0))));

        assert_eq!(expr.cook().uncook(), expr);
    }
}
