use std::{
    collections::HashSet,
    sync::{Arc, LazyLock, atomic::AtomicBool},
};

use spenso::{
    network::{
        Network,
        library::{DummyLibrary, symbolic::ETS},
        parsing::{NetworkParse, ParseSettings, SPENSO_TAG},
        store::NetworkStore,
    },
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    structure::{
        HasName, PermutedStructure, TensorStructure, ToSymbolic,
        abstract_index::AbstractIndex,
        permuted::Perm,
        representation::{LibraryRep, LibrarySlot, RepName},
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomType, AtomView, FunctionBuilder, Symbol, representation::FunView},
    coefficient::CoefficientView,
    function,
    id::{
        Condition, FilterFn, Match, MatchSettings, MatchStack, PatternRestriction, Replacement,
        WildcardRestriction,
    },
    symbol,
};
use tracing::warn;

use crate::tensor::{SymbolicNetParse, SymbolicTensor};

use crate::parsing_ind::Parsind;
use eyre::{Result, eyre};

use super::rep_symbols::RS;

pub struct MetricSymbols {
    pub dim: Symbol,
    // pub dot: Symbol,
    pub dummy: Symbol,
}

pub static MS: LazyLock<MetricSymbols> = LazyLock::new(|| MetricSymbols {
    dim: symbol!("spenso::dim"),
    // dot: symbol!("spenso::dot";Symmetric, Linear),
    dummy: symbol!("spenso::dummy"),
});

#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum CookingError {
    Add,
    Mul,
    Pow,
    RatCoeff,
    FiniteField,
    Float,
}

pub fn canonize_impl(view: AtomView) -> Atom {
    let lib = DummyLibrary::<SymbolicTensor>::new();
    let mut net = Network::<NetworkStore<SymbolicTensor, Atom>, _, Symbol>::try_from_view::<
        SymbolicTensor,
        _,
    >(view, &lib, &ParseSettings::default())
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

pub fn cook_function_view(view: AtomView) -> Result<Atom, CookingError> {
    match view {
        AtomView::Var(_) | AtomView::Num(_) => Ok(view.to_owned()),
        AtomView::Mul(_) => Err(CookingError::Mul),
        AtomView::Add(_) => Err(CookingError::Add),
        AtomView::Pow(_) => Err(CookingError::Pow),
        AtomView::Fun(f) => {
            let s = cook_function_impl(f)?;
            Ok(Atom::var(s))
        }
    }
}

pub fn cook_function_impl(fun: FunView) -> Result<Symbol, CookingError> {
    let mut name = fun.get_symbol().get_name().to_string();

    for arg in fun.iter() {
        name.push('_');
        match arg {
            AtomView::Fun(f) => {
                let arg_sym = cook_function_impl(f)?;
                name.push_str(arg_sym.get_stripped_name());
            }
            AtomView::Num(n) => match n.get_coeff_view() {
                CoefficientView::Indeterminate => name.push_str("ind"),
                CoefficientView::Infinity(_) => name.push_str("inf"),
                CoefficientView::FiniteField(_, _) => {
                    return Err(CookingError::FiniteField);
                }
                CoefficientView::Natural(n, d, imnum, imden) => {
                    name.push_str(&n.to_string());
                    if d != 1 {
                        name.push(':');
                        name.push_str(&d.to_string());
                    }
                    if imnum != 0 {
                        name.push('i');
                        name.push_str(&imnum.to_string());
                        if d != 1 {
                            name.push(':');
                            name.push_str(&imden.to_string());
                        }
                    }
                }
                CoefficientView::Float(_, _) => {
                    return Err(CookingError::Float);
                }
                CoefficientView::Large(r, imr) => {
                    let rat = r.to_rat();
                    name.push_str(&rat.numerator().to_string());
                    if !rat.is_integer() {
                        name.push(':');
                        name.push_str(&rat.denominator().to_string());
                    }

                    if !imr.is_zero() {
                        let rat = imr.to_rat();
                        name.push('i');
                        name.push_str(&rat.numerator().to_string());
                        if !rat.is_integer() {
                            name.push(':');
                            name.push_str(&rat.denominator().to_string());
                        }
                    }
                }
                CoefficientView::RationalPolynomial(_) => {
                    return Err(CookingError::RatCoeff);
                }
            },
            AtomView::Var(s) => {
                name.push_str(s.get_symbol().get_stripped_name());
            }
            AtomView::Pow(_) => {
                return Err(CookingError::Pow);
            }
            AtomView::Add(_) => {
                return Err(CookingError::Add);
            }
            AtomView::Mul(_) => {
                return Err(CookingError::Mul);
            }
        }
    }

    Ok(symbol!(&name))
}

pub fn cook_indices_impl(view: AtomView) -> Atom {
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
                && let Ok(aind) = cook_function_view(c[&RS.a_].as_view())
            {
                **out = i.to_symbolic([c[&RS.d_].clone(), aind]);
            }
        });
    }

    for i in LibraryRep::all_dualizables() {
        let ipat = i.to_symbolic([RS.d_, RS.a_]).to_pattern();
        // println!("{ipat}");
        let ipat_dual = i.dual().to_symbolic([RS.d_, RS.a_]).to_pattern();
        expr = expr.replace_map(|term, ctx, out| {
            if ctx.function_level == 1
                && let Some(c) = term.pattern_match(&ipat_dual, None, &settings).next()
                && let Ok(aind) = cook_function_view(c[&RS.a_].as_view())
            {
                // println!("{aind}");
                **out = i.dual().to_symbolic([c[&RS.d_].clone(), aind]);
            }
        });
        expr = expr.replace_map(|term, ctx, out| {
            if ctx.function_level == 1
                && let Some(c) = term.pattern_match(&ipat, None, &settings).next()
                && let Ok(aind) = cook_function_view(c[&RS.a_].as_view())
            {
                **out = i.to_symbolic([c[&RS.d_].clone(), aind]);
            }
        });
    }

    expr
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
    let mut expr = view.expand();
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

pub fn list_dangling_impl<Aind: ParseableAind + AbsInd>(view: AtomView) -> Vec<Atom> {
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

pub fn wrap_dummies_impl<Aind: ParseableAind + AbsInd>(view: AtomView, header: Symbol) -> Atom {
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

pub fn simplify_metrics_impl(view: AtomView) -> Atom {
    view.replace(
        function!(ETS.metric, RS.a_, SPENSO_TAG.self_dual_)
            * function!(RS.f_, RS.a___, SPENSO_TAG.self_dual_, RS.b___),
    )
    .repeat()
    .level_range((0, Some(0)))
    .with(function!(RS.f_, RS.a___, RS.a_, RS.b___))
    // // g(a_,dual_)*f_(a___,dind(dual_),b___) = f_(a___,a_,b___)
    // .replace(
    //     function!(ETS.metric, RS.a_, SPENSO_TAG.dualizable_([RS.a__]))
    //         * function!(
    //             RS.f_,
    //             RS.a___,
    //             SPENSO_TAG.dualizable_dual_([RS.a__]),
    //             RS.b___
    //         ),
    // )
    // .repeat()
    // .level_range((0, Some(0)))
    // .with(function!(RS.f_, RS.a___, RS.a_, RS.b___))
    // // g(a_,dind(dual_))*f_(a___,dual_,b___) = f_(a___,dual_,b___)
    // .replace(
    //     function!(ETS.metric, RS.a_, SPENSO_TAG.dualizable_dual_([RS.a__]))
    //         * function!(RS.f_, RS.a___, SPENSO_TAG.dualizable_([RS.a__]), RS.b___),
    // )
    // .repeat()
    // .level_range((0, Some(0)))
    // .with(function!(RS.f_, RS.a___, RS.a_, RS.b___))
    .replace(function!(
        ETS.metric,
        SPENSO_TAG.self_dual_([RS.d_, RS.i_]),
        SPENSO_TAG.self_dual_([RS.d_, RS.i_])
    ))
    .repeat()
    .level_range((0, Some(0)))
    .with(Atom::var(RS.d_))
    .replace(
        function!(
            ETS.metric,
            SPENSO_TAG.self_dual_([RS.d_, RS.i_]),
            SPENSO_TAG.self_dual_([RS.d_, RS.j_])
        )
        .npow(2),
    )
    .repeat()
    .level_range((0, Some(0)))
    .with(Atom::var(RS.d_))
    // .replace(function!(
    //     ETS.metric,
    //     SPENSO_TAG.dualizable_([RS.d_, RS.i_]),
    //     SPENSO_TAG.dualizable_dual_([RS.d_, RS.i_])
    // ))
    // .repeat()
    // .level_range((0, Some(0)))
    // .with(Atom::var(RS.d_))
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
    fn func_without_index(m: &MatchStack<'_>, fun_wild: Symbol, arg_wild: Symbol) -> Atom {
        match m.get(arg_wild).unwrap() {
            Match::FunctionName(a) => {
                panic!("Can't be a function")
            }
            Match::Single(a) => {
                let Match::FunctionName(f) = m.get(fun_wild).unwrap() else {
                    panic!("Not a function");
                };
                FunctionBuilder::new(*f).add_arg(a).finish()
            }

            Match::Multiple(a, args) => {
                if args.is_empty() {
                    m.get(fun_wild).unwrap().to_atom()
                } else if let Match::FunctionName(f) = m.get(fun_wild).unwrap() {
                    FunctionBuilder::new(*f).add_args(args).finish()
                } else {
                    panic!("Not a function");
                }
            }
        }
    }

    expr.replace(
        function!(RS.f_, RS.a___, SPENSO_TAG.self_dual_([RS.d_, RS.i_]))
            * function!(RS.g_, RS.b___, SPENSO_TAG.self_dual_([RS.d_, RS.i_])),
    )
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___) & not_slot(RS.b___))
    .repeat()
    .with_map(move |m| {
        let f = func_without_index(m, RS.f_, RS.a___);
        let g = func_without_index(m, RS.g_, RS.b___);
        // println!("{}", f);

        let rep = SPENSO_TAG
            .self_dual_([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);

        function!(ETS.metric, rep, f, g)
    })
    .replace(function!(RS.f_, RS.a___, SPENSO_TAG.self_dual_([RS.d_, RS.i_])).npow(2))
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___))
    .repeat()
    .with_map(move |m| {
        let f = func_without_index(m, RS.f_, RS.a___);

        let rep = SPENSO_TAG
            .self_dual_([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);

        function!(ETS.metric, rep, &f, &f)
    })
    .replace(
        function!(RS.f_, RS.a___, SPENSO_TAG.dualizable_([RS.d_, RS.i_]))
            * function!(RS.g_, RS.b___, SPENSO_TAG.dualizable_dual_([RS.d_, RS.i_])),
    )
    .level_range((0, Some(0)))
    .when(not_slot(RS.a___) & not_slot(RS.b___))
    .repeat()
    .with_map(move |m| {
        let f = func_without_index(m, RS.f_, RS.a___);
        let g = func_without_index(m, RS.g_, RS.b___);

        let rep = SPENSO_TAG
            .dualizable_([RS.d_])
            .to_pattern()
            .replace_wildcards_with_matches(m);

        function!(ETS.metric, rep, f, g)
    })
}

/// Trait for simplifying expressions involving metric tensors and converting
/// index contractions to dot product notation.
///
/// Provides methods for contracting indices with metric tensors (`g(mu, nu)`) or
/// identity tensors (`id(mu, nu)`), and for replacing contracted index patterns
/// (like `p(mu)*q(mu)`) with dot products (`dot(p, q)`).
pub trait MetricSimplifier {
    /// Simplifies contractions involving metric tensors (`g` or `metric`) and identity tensors (`id` or `𝟙`).
    ///
    /// Applies rules like `g(mu, nu) * p(nu) -> p(mu)`, `g(mu, mu) -> D`, etc.
    ///
    /// # Returns
    /// An [`Atom`] representing the expression after metric simplification.
    fn simplify_metrics(&self) -> Atom;

    fn expand_dots(&self) -> Result<Atom>;

    /// Converts contracted index patterns into dot product notation `dot(...)`.
    ///
    /// Replaces expressions like `p(mu) * q(mu)` or `p(mu) * M(mu, nu) * q(nu)` (implicitly via metric rules)
    /// with `dot(p, q)`. Assumes standard representations for vectors and tensors involved
    /// in the contractions.
    ///
    /// # Returns
    /// An [`Atom`] where contractions have been replaced by `dot` functions where possible.
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
        simplify_metrics_impl(self.as_view())
    }
}

impl MetricSimplifier for AtomView<'_> {
    fn to_dots(&self) -> Atom {
        to_dots_impl(*self)
    }

    fn expand_dots(&self) -> Result<Atom> {
        let set = ParseSettings::default();
        let pat = function!(SPENSO_TAG.dot, RS.f_, RS.g_, RS.h_).to_pattern();
        let has_errored = Arc::new(AtomicBool::new(false));
        let has_errored_for_closure = has_errored.clone();
        let out = self.replace(pat.clone()).with_map(move |a| {
            let filled = pat.replace_wildcards_with_matches(a);

            match filled.parse_to_atom_net::<AbstractIndex>(&set) {
                Err(e) => {
                    warn!("Failed to parse network from {}:{}", filled, e);
                    has_errored_for_closure.store(true, std::sync::atomic::Ordering::Relaxed);
                    filled
                }
                Ok(mut net) => {
                    net.simple_execute();
                    net.result_scalar().unwrap().into()
                }
            }
        });

        if has_errored.load(std::sync::atomic::Ordering::Relaxed) {
            Err(eyre!("Failed to parse network"))
        } else {
            Ok(out)
        }
    }
    fn simplify_metrics(&self) -> Atom {
        simplify_metrics_impl(*self)
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

    use crate::{IndexTooling, representations::Bispinor, test::test_initialize};

    use super::*;

    use spenso::{
        network::parsing::ShadowedStructure,
        shadowing::symbolica_utils::AtomCoreExt,
        structure::{
            IndexlessNamedStructure, PermutedStructure,
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

        println!("{}\n", f.structure.expression);

        let f_p = f.clone().permute_inds();

        println!("{}\n", f_p.expression);
        println!("{}\n", f_p.expression.simplify_metrics());
        let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
            &f_p.expression.simplify_metrics(),
        )
        .unwrap();

        assert_eq!(f.index_permutation, f_parsed.index_permutation);
        assert!(f_parsed.rep_permutation.is_identity());

        let f_p = f.clone().permute_reps_wrapped().permute_inds();

        let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
            &f_p.expression.simplify_metrics(),
        )
        .unwrap();

        println!("{}\n", f_p.expression);
        println!("{}\n", f_p.expression.simplify_metrics());
        assert_eq!(f.index_permutation, f_parsed.index_permutation);
        assert_eq!(f.rep_permutation, f_parsed.rep_permutation);
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
        let a =
            parse_lit!(P(wrong::mink(4, -1 * g(2)), spenso::mink(4, 2)) * P(spenso::mink(4, 2)))
                .to_dots();
        assert_eq!(
            a,
            parse_lit!(spenso::g(
                spenso::mink(4),
                idenso::P,
                idenso::P(wrong::mink(4, -idenso::g(2)))
            ))
        );
        insta::assert_snapshot!(a.expand_dots().unwrap().to_bare_ordered_string(),@"-1*P(cind(1))*P(mink(4,-1*g(2)),cind(1))+-1*P(cind(2))*P(mink(4,-1*g(2)),cind(2))+-1*P(cind(3))*P(mink(4,-1*g(2)),cind(3))+P(cind(0))*P(mink(4,-1*g(2)),cind(0))");

        let a = parse_lit!(k1(spenso::mink(4, mu1)) ^ 2).to_dots();
        insta::assert_snapshot!(a.to_bare_ordered_string(),@"g(k1,k1,mink(4))");
    }

    #[test]
    fn true_cooking() {
        test_initialize();
        let expr = parse_lit!(spenso::g(spenso::mink(4, f(0))));

        println!("{}", expr.cook());
        println!("{}", expr.cook().uncook());
        // println!("{}", symbol!(&a, data = UserData::Atom(expr)));
    }
}
