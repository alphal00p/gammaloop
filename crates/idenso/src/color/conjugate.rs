use std::sync::Arc;
use std::sync::LazyLock;

use crate::tensor::SymbolicTensor;
use crate::tensor::inference::TensorInferenceError;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::slot::DualSlotTo;

use spenso::network::library::function_lib::INBUILTS;
use spenso::network::parsing::{AtomStructureExt, StrictTensorFilter};
use spenso::network::tags::SPENSO_TAG as T;
use spenso::shadowing;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::representation::{LibraryRep, RepName};
use spenso::structure::slot::{IsAbstractSlot, Slot};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::id::Replacement;

use crate::epsilon::EPSILON_SYMBOL;
use crate::rep_symbols::RS;
use crate::representations::ColorFundamental;
use crate::{IndexTooling, color_t};

use super::CS;
use super::simplify::{color_generator_adjoint_view, fundamental_chain_dimension_view};

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
    fn run(
        expression: AtomView<'_>,
        observer: &mut impl FnMut(AtomView<'_>, AtomView<'_>),
    ) -> Option<Atom> {
        let mut scalar = true;
        expression.visitor(&mut |value| {
            if matches!(value, AtomView::Fun(fun) if fun.get_symbol().is_scalar()) {
                return false;
            }
            scalar &= !matches!(value, AtomView::Var(var) if [T.chain_in, T.chain_out].contains(&var.get_symbol()))
                && !matches!(value, AtomView::Fun(fun) if [*shadowing::SYM, *shadowing::ANTISYM, *shadowing::CYCLIC].contains(&fun.get_symbol()))
                && !value.is_tensorial(StrictTensorFilter::Tagged)
                && !value.is_tensorial(StrictTensorFilter::ContainsReps);
            scalar
        });
        if scalar {
            // Symbolica owns scalar conjugation, including powers and functions.
            // Matrix placeholders and projectors retain their tensor semantics.
            return Some(expression.spenso_conj());
        }
        let word = match expression {
            AtomView::Mul(product) => {
                return product.iter().try_fold(Atom::one(), |result, factor| {
                    Some(result * Self::run(factor, observer)?)
                });
            }
            AtomView::Add(sum) => {
                return sum.iter().try_fold(Atom::Zero, |result, term| {
                    Some(result + Self::run(term, observer)?)
                });
            }
            AtomView::Fun(word) => word,
            _ => return None,
        };
        if color_generator_adjoint_view(expression).is_some() {
            return Some(expression.to_owned());
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
            .map(|factor| Self::run(factor, observer))
            .collect::<Option<Vec<_>>>()?;
        let result = if word.get_symbol() == T.trace {
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
        };
        observer(expression, result.as_view());
        Some(result)
    }

    fn words(
        expression: AtomView<'_>,
        observer: &mut impl FnMut(AtomView<'_>, AtomView<'_>),
    ) -> Atom {
        expression
            .replace_map(|arg, _context, out| {
                let AtomView::Fun(conjugate) = arg else {
                    return;
                };
                if conjugate.get_symbol() != INBUILTS.conj || conjugate.get_nargs() != 1 {
                    return;
                }
                let Some(AtomView::Fun(word)) = conjugate.iter().next() else {
                    return;
                };
                if word.get_symbol() == *EPSILON_SYMBOL && word.get_nargs() == 3 {
                    let triplet = ColorFundamental {}.new_rep(3).cast::<LibraryRep>();
                    let reps = word
                        .iter()
                        .map(|arg| {
                            Slot::<LibraryRep, AbstractIndex>::try_from(arg).map(|slot| slot.rep())
                        })
                        .collect::<Result<Vec<_>, _>>();
                    if let Ok(reps) = reps
                        && (reps.iter().all(|rep| *rep == triplet)
                            || reps.iter().all(|rep| *rep == triplet.dual()))
                    {
                        **out = word.as_view().to_owned();
                    }
                } else if [T.chain, T.trace].contains(&word.get_symbol())
                    && let Some(result) = Self::run(word.as_view(), observer)
                {
                    **out = result;
                }
            })
            .replace_map(|value, _, out| {
                let AtomView::Fun(function) = value else {
                    return;
                };
                if function.get_symbol() != CS.t || function.get_nargs() != 3 {
                    return;
                }
                let result = value.replace_multiple(&*COLOR_CONJ_GENERATOR_TRANSPOSITIONS);
                if result.as_view() != value {
                    observer(value, result.as_view());
                    **out = result;
                }
            })
    }
}

pub fn color_conj_impl(expression: AtomView<'_>) -> Atom {
    ColorConjugator::words(expression, &mut |_, _| {})
        .replace_multiple(&*COLOR_CONJ_REPRESENTATION_SWAPS)
}

impl SymbolicTensor<PartialStructure> {
    /// Conjugate the actual producer tree once while carrying each primitive's
    /// admitted interface. Colour-word reversal and representation duality are
    /// separate typed transformations, so newly emitted templates stay known.
    pub(crate) fn spenso_conjugated(&self) -> Result<Self, TensorInferenceError> {
        let (expression, mapped) = self.map_admitted_expression(
            |value| {
                Ok(matches!(
                    value,
                    AtomView::Num(_) | AtomView::Var(_) | AtomView::Fun(_)
                )
                .then(|| {
                    value
                        .conj()
                        .replace(Atom::var(RS.a__).conj())
                        .with(INBUILTS.conj(RS.a__))
                }))
            },
            |interface| Ok(interface.clone()),
        )?;
        let native = Self::from_validated_parts(
            expression,
            self.structure.clone(),
            Arc::clone(
                self.proofs
                    .occurrences
                    .get()
                    .expect("conjugation source carries index history"),
            ),
        );
        let native = self.finish_mapped_result("native conjugation", native, mapped)?;
        let mut failure = None;
        let expression =
            ColorConjugator::words(native.expression().as_view(), &mut |source, target| {
                if let Some(interface) = native.admitted_interface(source) {
                    native.proofs.retain_interface(target, interface);
                }
                if let AtomView::Fun(function) = target
                    && function.get_symbol() == CS.t
                    && function.get_nargs() == 3
                    && let Err(error) = native.proofs.retain_leaf_ports(
                        target,
                        function.iter().map(|arg| arg.to_owned()).collect(),
                    )
                {
                    failure = Some(error);
                }
            });
        if let Some(error) = failure {
            return Err(error);
        }
        // Transposition alone temporarily reverses base and dual positions
        // in the generator signature. Finish duality before publishing it.
        let words = expression;
        let fundamental = LibraryRep::from(ColorFundamental {});
        let dual_slot = |slot: Slot<LibraryRep, AbstractIndex>| {
            if slot.rep_name().base() == fundamental.base() {
                slot.dual()
            } else {
                slot
            }
        };
        let dual_interface = |interface: &PartialStructure| {
            Ok(PartialStructure::from_logical_slots(
                interface.logical_slots().into_iter().map(|slot| {
                    if slot.rep_name().base() == fundamental.base() {
                        slot.dual()
                    } else {
                        slot
                    }
                }),
            ))
        };
        let structure = dual_interface(native.structure())?;
        let (expression, mapped) = native.map_producer_draft(
            words.as_view(),
            |value| {
                if value.get_symbol().is_some_and(|head| head.is_scalar()) {
                    return Ok(Some(value.to_owned()));
                }
                Ok(Slot::<LibraryRep, AbstractIndex>::try_from(value)
                    .ok()
                    .filter(|slot| slot.rep_name().base() == fundamental.base())
                    .map(|slot| slot.dual().to_atom()))
            },
            dual_interface,
        )?;
        let occurrences = native
            .explicit_index_occurrences()?
            .mapped(|slot| Ok::<_, TensorInferenceError>(dual_slot(slot)))?;
        let result = Self::from_validated_parts(expression, structure, Arc::new(occurrences));
        native.finish_mapped_result("colour conjugation", result, mapped)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::IndexTooling;
    use crate::color::ColorSimplifier;
    use crate::epsilon::EpsilonSimplifier;
    use crate::shorthands::UndoShorthands;
    use crate::tensor::inference::tests::INFERENCE_CALLS;
    use crate::test_support::test_initialize;
    use spenso::shadowing::ProjectorExpander;
    use spenso::structure::abstract_index::AbstractIndex;
    use symbolica::parse;

    #[test]
    fn typed_colour_conjugation_retains_subregions_without_inference() {
        test_initialize();
        for expression in [
            parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))",
                default_namespace = "spenso"
            ),
            parse!(
                "chain(cof(3,j),dind(cof(3,i)),t(coad(8,a),in,out),t(coad(8,b),in,out))",
                default_namespace = "spenso"
            ),
        ] {
            let expected = expression.spenso_conj();
            let source = SymbolicTensor::infer(expression).unwrap();
            INFERENCE_CALLS.with(|count| count.set(0));
            let conjugated = source.spenso_conjugated().unwrap();
            assert_eq!(conjugated.expression(), &expected);
            let roundtrip = conjugated.spenso_conjugated().unwrap();
            assert_eq!(roundtrip.expression(), source.expression());
            assert_eq!(roundtrip.structure(), source.structure());
            assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        }
    }

    #[test]
    fn color_epsilon_conjugation_dualizes_without_reversing_indices() {
        test_initialize();
        let _ = *EPSILON_SYMBOL;
        let epsilon = parse!(
            "epsilon(cof(3,i),cof(3,j),cof(3,k))",
            default_namespace = "spenso"
        );
        let barred = parse!(
            "epsilon(dind(cof(3,i)),dind(cof(3,j)),dind(cof(3,k)))",
            default_namespace = "spenso"
        );
        assert_eq!(epsilon.spenso_conj(), barred);
        assert_eq!(barred.spenso_conj(), epsilon);
        assert_eq!((Atom::i() * &epsilon).spenso_conj(), -Atom::i() * &barred);
        assert_eq!((&epsilon * &barred).simplify_epsilon(), Atom::num(6));
    }

    #[test]
    fn color_epsilon_conjugation_preserves_other_conventions() {
        test_initialize();
        let _ = *EPSILON_SYMBOL;
        for source in [
            "epsilon(mink(4,i),mink(4,j),mink(4,k),mink(4,l))",
            "epsilon(cof(4,i),cof(4,j),cof(4,k))",
            "epsilon(cof(3,i),cof(3,j))",
            "epsilon(cof(3,i),cof(3,j),dind(cof(3,k)))",
        ] {
            let expression = parse!(source, default_namespace = "spenso");
            let conjugate = expression.spenso_conj();
            assert!(conjugate.contains_symbol(INBUILTS.conj));
            assert_eq!(conjugate.spenso_conj(), expression);
        }
    }

    #[test]
    fn adjoint_color_generator_has_positive_su3_norm() {
        test_initialize();
        let generator = parse!(
            "t(coad(8,a),cof(3,i),dind(cof(3,j)))",
            default_namespace = "spenso"
        );
        let adjoint = generator.dirac_adjoint::<AbstractIndex>(false).unwrap();
        assert_eq!(
            adjoint,
            parse!(
                "t(coad(8,a),cof(3,j),dind(cof(3,i)))",
                default_namespace = "spenso"
            )
        );
        let norm = crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((generator * adjoint).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(crate::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .into_expression()
                .to_cof_dimension_invariants())
            .as_atom_view()
            .to_owned(),
        )
        .unwrap()
        .contract(crate::tensor::ContractSettings::default().without_rank_one_tensors())
        .unwrap()
        .into_expression();
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
            crate::tensor::SymbolicTensor::infer((conjugate).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(crate::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .into_expression(),
            crate::tensor::SymbolicTensor::infer(
                (explicit.spenso_conj()).as_atom_view().to_owned()
            )
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression()
        );
        let norm = crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((chain * conjugate).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(crate::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .into_expression()
                .to_cof_dimension_invariants())
            .as_atom_view()
            .to_owned(),
        )
        .unwrap()
        .contract(crate::tensor::ContractSettings::default().without_rank_one_tensors())
        .unwrap()
        .into_expression();
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
            crate::tensor::SymbolicTensor::infer((conjugate).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(crate::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .into_expression(),
            crate::tensor::SymbolicTensor::infer(
                (explicit.spenso_conj()).as_atom_view().to_owned()
            )
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression()
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
    fn compact_color_words_conjugate_scalar_powers_and_functions() {
        test_initialize();
        let generator = parse!("t(coad(8,a),in,out)", default_namespace = "spenso");
        for source in [
            "chain(cof(3,i),dind(cof(3,j)),t(coad(8,a),in,out),t(coad(8,b),in,out))",
            "trace(cof(3),cyclic(t(coad(8,a),in,out),t(coad(8,b),in,out),t(coad(8,c),in,out)))",
            "trace(cof(3),sym(t(coad(8,a),in,out),antisym(t(coad(8,b),in,out),t(coad(8,c),in,out))))",
            "trace(cof(3),cyclic(sym(z,w),t(coad(8,a),in,out),t(coad(8,b),in,out)))",
            "chain(cof(3,i),dind(cof(3,j)),cyclic(z,w),t(coad(8,a),in,out),t(coad(8,b),in,out))",
            "trace(cof(3),cyclic(z*sym(w,x),t(coad(8,a),in,out),t(coad(8,b),in,out)))",
        ] {
            let word = parse!(source, default_namespace = "spenso");
            for scalar in [
                "z^2",
                "z^-1",
                "(z+w)^2",
                "z^(1/2)",
                "z^w",
                "h(z)",
                "h(z)^-1",
                "spenso::conj(h(z))",
                "cas(2,cof(3))",
            ] {
                let coefficient = parse!(scalar, default_namespace = "spenso");
                let weighted = word
                    .replace(generator.to_pattern())
                    .with(&coefficient * &generator);
                let expected = word
                    .spenso_conj()
                    .replace(generator.to_pattern())
                    .with(coefficient.spenso_conj() * &generator);
                assert_eq!(weighted.spenso_conj(), expected, "{source}, {scalar}");
                assert_eq!(
                    expected.spenso_conj(),
                    weighted,
                    "involution: {source}, {scalar}"
                );
                assert_eq!(
                    expected.expand_projectors().expand(),
                    weighted.expand_projectors().spenso_conj().expand(),
                    "projector expansion: {source}, {scalar}"
                );
            }
        }
    }

    #[test]
    fn compact_color_conjugation_keeps_unknown_matrix_factors_explicit() {
        test_initialize();
        let word = parse!(
            "trace(cof(3),cyclic(A,t(coad(8,a),in,out),t(coad(8,b),in,out)))",
            default_namespace = "spenso"
        );
        for matrix in [
            "A(in,out)",
            "A(coad(8,c))",
            "wrapper(A(in,out))",
            "A(in,out)^2",
            "t(coad(8,c),in,out)^2",
            "z^A(in,out)",
        ] {
            let matrix = parse!(matrix, default_namespace = "spenso");
            let expression = word
                .replace(parse!("A", default_namespace = "spenso"))
                .with(matrix);
            assert_eq!(expression.spenso_conj(), INBUILTS.conj(&expression));
        }
        let tagged = symbolica::function!(spenso::tensor_symbol!("unknown_color_matrix"));
        let expression = word
            .replace(parse!("A", default_namespace = "spenso"))
            .with(tagged);
        assert_eq!(expression.spenso_conj(), INBUILTS.conj(&expression));
    }

    #[test]
    fn color_conjugation_preserves_scalar_representation_invariants() {
        test_initialize();
        let invariant = parse!("cas(2,cof(3))", default_namespace = "spenso");
        assert_eq!(invariant.spenso_conj(), invariant);
    }
}
