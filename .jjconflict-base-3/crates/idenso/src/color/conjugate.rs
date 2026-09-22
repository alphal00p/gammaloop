use std::sync::LazyLock;

use spenso::{
    network::{library::function_lib::INBUILTS, tags::SPENSO_TAG as T},
    shadowing,
    structure::representation::RepName,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    id::Replacement,
};

use crate::{IndexTooling, color_t, rep_symbols::RS, representations::ColorFundamental};

use super::{
    CS,
    simplify::{color_generator_adjoint_view, fundamental_chain_dimension_view},
};

static COLOR_CONJ_GENERATOR_TRANSPOSITIONS: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let cof = ColorFundamental {};
    let coaf = ColorFundamental {}.dual();

    [Replacement::new(
        color_t!(
            RS.i__,
            cof.to_symbolic([RS.d_, RS.a_]),
            coaf.to_symbolic([RS.d_, RS.b_]),
        )
        .to_pattern(),
        color_t!(
            RS.i__,
            coaf.to_symbolic([RS.d_, RS.b_]),
            cof.to_symbolic([RS.d_, RS.a_]),
        ),
    )]
});

static COLOR_CONJ_REPRESENTATION_SWAPS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let cof = ColorFundamental {};
    let coaf = ColorFundamental {}.dual();

    [
        Replacement::new(
            coaf.to_symbolic([RS.d_, RS.a_]).to_pattern(),
            cof.to_symbolic([RS.d_, RS.a_]),
        ),
        Replacement::new(
            cof.to_symbolic([RS.d_, RS.a_]).to_pattern(),
            coaf.to_symbolic([RS.d_, RS.a_]),
        ),
    ]
});

struct ColorConjugator;

impl ColorConjugator {
    // Adjoint a compact generator word without expanding its projectors.
    fn run(expression: AtomView<'_>) -> Option<Atom> {
        let word = match expression {
            AtomView::Num(_) | AtomView::Var(_) => return Some(expression.spenso_conj()),
            AtomView::Mul(product) => {
                return product.iter().try_fold(Atom::one(), |result, factor| {
                    Some(result * Self::run(factor)?)
                });
            }
            AtomView::Add(sum) => {
                return sum
                    .iter()
                    .try_fold(Atom::Zero, |result, term| Some(result + Self::run(term)?));
            }
            AtomView::Fun(word) => word,
            _ => return None,
        };
        if color_generator_adjoint_view(expression).is_some() {
            return Some(expression.to_owned());
        }
        if word.get_symbol() == INBUILTS.conj
            && word.get_nargs() == 1
            && matches!(word.iter().next(), Some(AtomView::Var(_)))
        {
            return Some(expression.spenso_conj());
        }
        let args = word.iter().collect::<Vec<_>>();
        let (prefix, factors) = if word.get_symbol() == T.chain {
            let [start, end, factors @ ..] = args.as_slice() else {
                return None;
            };
            let _ = fundamental_chain_dimension_view(*start, *end)?;
            // Hermitian generators reverse the word and transpose its
            // endpoints. The shared representation swap below dualizes them.
            (vec![*end, *start], factors.to_vec())
        } else if let Some((rep, factors)) = shadowing::trace_parts(word) {
            if !matches!(rep, AtomView::Fun(f) if f.get_symbol() == CS.fundamental_rep && f.get_nargs() == 1)
            {
                return None;
            }
            (vec![rep], factors)
        } else if [*shadowing::SYM, *shadowing::ANTISYM, *shadowing::CYCLIC]
            .contains(&word.get_symbol())
        {
            (vec![], args)
        } else {
            return None;
        };
        let factors = factors
            .into_iter()
            .rev()
            .map(Self::run)
            .collect::<Option<Vec<_>>>()?;
        Some(if word.get_symbol() == T.trace {
            shadowing::trace(prefix[0], factors)
        } else {
            prefix
                .into_iter()
                .map(|arg| arg.to_owned())
                .chain(factors)
                .fold(
                    FunctionBuilder::new(word.get_symbol()),
                    |builder, factor| builder.add_arg(factor),
                )
                .finish()
        })
    }
}

pub fn color_conj_impl(expression: AtomView<'_>) -> Atom {
    expression
        .replace_map(|arg, _context, out| {
            let AtomView::Fun(conjugate) = arg else {
                return;
            };
            if conjugate.get_symbol() == INBUILTS.conj
                && conjugate.get_nargs() == 1
                && let Some(AtomView::Fun(word)) = conjugate.iter().next()
                && [T.chain, T.trace].contains(&word.get_symbol())
                && let Some(result) = ColorConjugator::run(AtomView::Fun(word))
            {
                **out = result;
            }
        })
        .replace_multiple(&*COLOR_CONJ_GENERATOR_TRANSPOSITIONS)
        .replace_multiple(&*COLOR_CONJ_REPRESENTATION_SWAPS)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        IndexTooling,
        color::ColorSimplifier,
        shorthands::{UndoShorthands, metric::MetricSimplifier},
        test_support::test_initialize,
    };
    use spenso::{shadowing::ProjectorExpander, structure::abstract_index::AbstractIndex};
    use symbolica::parse;

    #[test]
    fn adjoint_color_generator_has_positive_su3_norm() {
        test_initialize();
        let generator = parse!(
            "t(coad(8,a),cof(3,i),dind(cof(3,j)))",
            default_namespace = "spenso"
        );
        let adjoint = generator.dirac_adjoint::<AbstractIndex>().unwrap();
        assert_eq!(
            adjoint,
            parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))",
                default_namespace = "spenso"
            )
        );
        let norm = (generator * adjoint)
            .simplify_color()
            .to_cof_dimension_invariants()
            .simplify_metrics();
        assert_eq!(norm, Atom::num(4));
    }

    #[test]
    fn color_conjugation_dualizes_open_network_and_is_involutive() {
        test_initialize();
        let network = Atom::i()
            * parse!(
                "t(coad(8,a),cof(3,i),dind(cof(3,j)))*g(cof(3,j),dind(cof(3,k)))",
                default_namespace = "spenso"
            );
        let conjugate = network.spenso_conj();
        let expected = -Atom::i()
            * parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))*g(cof(3,k),dind(cof(3,j)))",
                default_namespace = "spenso"
            );
        assert_eq!(conjugate, expected);
        assert_eq!(conjugate.spenso_conj(), network);
    }

    #[test]
    fn compact_color_chain_conjugation_matches_explicit_network() {
        test_initialize();
        let chain = parse!(
            "chain(cof(3,i),dind(cof(3,j)),t(coad(8,a),in,out),t(coad(8,b),in,out))",
            default_namespace = "spenso"
        );
        let conjugate = chain.spenso_conj();
        assert_eq!(
            conjugate,
            parse!(
                "chain(cof(3,j),dind(cof(3,i)),t(coad(8,b),in,out),t(coad(8,a),in,out))",
                default_namespace = "spenso"
            )
        );
        assert_eq!(conjugate.spenso_conj(), chain);
        let explicit = chain.undo_chain::<AbstractIndex>().unwrap();
        assert_eq!(
            conjugate.simplify_color(),
            explicit.spenso_conj().simplify_color()
        );
        let norm = (chain * conjugate)
            .simplify_color()
            .to_cof_dimension_invariants()
            .simplify_metrics();
        assert_eq!(norm, parse!("16/3"));
    }

    #[test]
    fn compact_color_trace_conjugation_reverses_order() {
        test_initialize();
        let trace = parse!(
            "trace(cof(3),cyclic(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
            default_namespace = "spenso"
        );
        let conjugate = trace.spenso_conj();
        let reversed = parse!(
            "trace(cof(3),cyclic(t(coad(8,c),in,out),t(coad(8,b),in,out),t(coad(8,a),in,out)))",
            default_namespace = "spenso"
        );
        assert_ne!(
            trace, reversed,
            "a three-generator trace is not generally real"
        );
        assert_eq!(conjugate, reversed);
        assert_eq!(conjugate.spenso_conj(), trace);
        assert_eq!((Atom::i() * &trace).spenso_conj(), -Atom::i() * &conjugate);
        let explicit = trace.undo_trace::<AbstractIndex>().unwrap();
        assert_eq!(
            conjugate.simplify_color(),
            explicit.spenso_conj().simplify_color()
        );
    }

    #[test]
    fn projected_color_words_match_expanded_conjugation() {
        test_initialize();
        for source in [
            "chain(cof(3,i),dind(cof(3,j)),sym(t(coad(8,a),in,out),t(coad(8,b),in,out)),t(coad(8,c),in,out))",
            "chain(cof(3,i),dind(cof(3,j)),antisym(t(coad(8,a),in,out),t(coad(8,b),in,out)),t(coad(8,c),in,out))",
            "trace(cof(3),sym(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
            "trace(cof(3),cyclic(antisym(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out))))",
            "trace(cof(3),sym(t(coad(8,a),in,out),cyclic(t(coad(8,b),in,out),t(coad(8,c),in,out),t(coad(8,d),in,out))))",
            "trace(cof(3),sym((1+1𝑖)*t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
            "trace(cof(3),sym(z*t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
        ] {
            let expression = parse!(source, default_namespace = "spenso");
            let conjugate = expression.spenso_conj();
            assert_eq!(
                conjugate.expand_projectors().expand(),
                expression.expand_projectors().spenso_conj().expand(),
                "compact projector conjugation: {source}"
            );
            assert_eq!(
                conjugate.spenso_conj().expand_projectors().expand(),
                expression.expand_projectors().expand(),
                "projector conjugation must be involutive: {source}"
            );
        }
    }

    #[test]
    fn color_projector_trace_reality_follows_reversal_parity() {
        test_initialize();
        let symmetric = parse!(
            "trace(cof(3),sym(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
            default_namespace = "spenso"
        );
        let antisymmetric = parse!(
            "trace(cof(3),cyclic(antisym(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out))))",
            default_namespace = "spenso"
        );
        assert_eq!(symmetric.spenso_conj(), symmetric);
        assert_eq!(
            antisymmetric.spenso_conj().expand_projectors().expand(),
            (-antisymmetric.expand_projectors()).expand()
        );
        // Unknown matrix factors must not be assumed Hermitian.
        let unknown = parse!(
            "trace(cof(3),cyclic(A(in,out)))",
            default_namespace = "spenso"
        );
        assert_eq!(unknown.spenso_conj(), INBUILTS.conj(unknown));
    }

    #[test]
    fn color_conjugation_preserves_scalar_representation_invariants() {
        test_initialize();
        let invariant = parse!("cas(2,cof(3))", default_namespace = "spenso");
        assert_eq!(invariant.spenso_conj(), invariant);
    }
}
