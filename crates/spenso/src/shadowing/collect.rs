mod tape;
pub use tape::{TermLeaf, TermTape};

use std::collections::BTreeMap;

use crate::{
    network::{
        library::symbolic::ETS,
        parsing::{StrictTensorFilter, structure_inference::TensorialSyntax},
        tags::SPENSO_TAG,
    },
    shadowing::static_symbols::W_,
    structure::{representation::LibraryRep, slot::SlotMatcher},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    coefficient::CoefficientView,
    id::{AliasedAtom, Context},
    symbol,
    utils::Settable,
};

crate::symbolica_init_lazy_static! {
    pub static COLLECT, COLLECT_INNER: Symbol = || symbol!("spenso::collect");
}

/// Selects which tensorial subexpressions are protected during collection.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TensorCollectFilter<const N: usize> {
    /// Collect every tensor head that is tagged with the tensor tag.
    TaggedTensors,
    /// Collect every syntactically tensorial subexpression.
    Tensors,
    /// Collect tensorial subexpressions that contain this representation.
    Reps([LibraryRep; N]),
    /// Collect only metric tensors.
    Metrics,
    /// Collect only chains and traces of tensors.
    ChainsAndTraces,
}

pub trait Collectable {
    /// Collect selected leaves while retaining opaque coefficients as alias definitions.
    fn collect_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> AliasedAtom;
    fn expand_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom;
    fn collect_collects(self) -> Atom;
    fn map_collects<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(self, map: F) -> Atom;
    fn wrap_in_collect(self) -> Atom;
    fn unwrap_collect(self) -> Atom;
}
impl Collectable for Atom {
    fn collect_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> AliasedAtom {
        self.as_view().collect_with_map(map)
    }

    fn expand_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom {
        self.as_view().expand_with_map(map)
    }

    fn map_collects<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(self, map: F) -> Atom {
        self.as_view().map_collects(map)
    }

    fn collect_collects(self) -> Atom {
        self.as_view().collect_collects()
    }

    fn unwrap_collect(self) -> Atom {
        self.as_view().unwrap_collect()
    }

    fn wrap_in_collect(self) -> Atom {
        self.as_view().wrap_in_collect()
    }
}
impl Collectable for AtomView<'_> {
    fn collect_collects(self) -> Atom {
        // A repeated tensor belongs to the complete collected monomial. Keep
        // its power compressed while exposing all copies to mapping callbacks.
        self.replace_map(|arg, _context, out| {
            let AtomView::Pow(power) = arg else { return };
            let (AtomView::Fun(base), AtomView::Num(exponent)) = power.get_base_exp() else {
                return;
            };
            let coefficient = exponent.get_coeff_view();
            if base.get_symbol() != *COLLECT
                || base.get_nargs() != 1
                || !matches!(
                    coefficient,
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                )
                || !coefficient.is_integer()
                || coefficient.to_owned().is_negative()
                || coefficient.to_owned().is_zero()
            {
                return;
            }
            **out = base
                .iter()
                .next()
                .unwrap()
                .pow(AtomView::Num(exponent))
                .wrap_in_collect();
        })
        .replace(COLLECT.call(W_.a_) * COLLECT.call(W_.b_))
        .repeat()
        .with(COLLECT.call(W_.a_ * W_.b_))
    }

    fn map_collects<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        self,
        mut map: F,
    ) -> Atom {
        self.replace_map(|arg, context, out| {
            let AtomView::Fun(a) = arg else { return };
            if a.get_symbol() == *COLLECT {
                map(arg, context, out)
            }
        })
    }
    fn expand_with_map(self, mut matches: impl FnMut(AtomView<'_>) -> bool) -> Atom {
        let mut hit = false;
        let wrapped = self.replace_map(|arg, _context, out| {
            if matches(arg) {
                hit = true;
                **out = arg.wrap_in_collect()
            }
        });

        if !hit {
            return self.to_owned();
        }
        wrapped.expand_in(*COLLECT)
    }
    fn collect_with_map(self, mut matches: impl FnMut(AtomView<'_>) -> bool) -> AliasedAtom {
        let mut hit = false;
        let wrapped = self.replace_map(|arg, _context, out| {
            if matches(arg) {
                hit = true;
                **out = arg.wrap_in_collect()
            }
        });

        if !hit {
            return self.to_owned().into();
        }
        // Keep complete unselected coefficients opaque. Symbolica's polynomial
        // collector otherwise statistically tests their growing sums for zero,
        // repeatedly traversing the whole momentum numerator. Collection needs
        // only the selected tensor factors; it is not scalar simplification.
        let mut used = wrapped.get_all_symbols(true);
        let mut aliases = BTreeMap::<Atom, Atom>::new();
        let mut serial = 0usize;
        let protected = wrapped.replace_map(|arg, _context, out| {
            if matches!(arg, AtomView::Fun(fun) if fun.get_symbol() == *COLLECT) {
                // Setting the unchanged output also stops traversal into the
                // selected tensor, preserving its complete slots and payload.
                **out = arg.to_owned();
            } else if !matches!(arg, AtomView::Num(_) | AtomView::Var(_))
                && !arg.contains_symbol(*COLLECT)
            {
                let alias = aliases.entry(arg.to_owned()).or_insert_with(|| {
                    loop {
                        let symbol = symbol!(&format!("spenso::collect_coefficient_{serial}"));
                        serial += 1;
                        if used.insert(symbol) {
                            break Atom::var(symbol);
                        }
                    }
                });
                **out = alias.clone();
            }
        });
        let mut protected = AliasedAtom::from(protected);
        for (original, alias) in aliases {
            protected.register_alias(alias, original);
        }
        // Keep coefficient definitions on the existing alias owner. Callers
        // resolve explicitly before callbacks that require complete terms.
        protected.map_root(|root| {
            let AtomView::Mul(product) = root.as_view() else {
                return root.collect_symbol::<i16>(*COLLECT).collect_collects();
            };
            // Keep a shared coefficient outside the tensor sum. Collect the
            // entire remaining product so callbacks still see every factor
            // of a contraction, including tensors multiplying a sum.
            let (tensors, coefficients): (Vec<_>, Vec<_>) = product
                .iter()
                .partition(|factor| factor.contains_symbol(*COLLECT));
            tensors
                .into_iter()
                .product::<Atom>()
                .collect_symbol::<i16>(*COLLECT)
                .collect_collects()
                * coefficients.into_iter().product::<Atom>()
        })
    }

    fn unwrap_collect(self) -> Atom {
        fn collect_inner(arg: AtomView<'_>) -> Option<AtomView<'_>> {
            let AtomView::Fun(fun) = arg else {
                return None;
            };
            if fun.get_symbol() != *COLLECT || fun.get_nargs() != 1 {
                return None;
            }
            fun.iter().next()
        }

        self.replace_map(|arg, _context, out| {
            if let Some(inner) = collect_inner(arg) {
                **out = inner.to_owned();
                return;
            }

            if let AtomView::Pow(pow) = arg {
                let (base, exponent) = pow.get_base_exp();
                if let Some(inner) = collect_inner(base) {
                    **out = inner.to_owned().pow(exponent.to_owned());
                }
            }
        })
    }

    fn wrap_in_collect(self) -> Atom {
        FunctionBuilder::new(*COLLECT).add_arg(self).finish()
    }
}

impl<const N: usize> TensorCollectFilter<N> {
    fn collect(self, expression: AtomView<'_>) -> Atom {
        expression
            .collect_with_map(|a| self.matches(a))
            .into_inner()
            .unwrap_collect()
    }

    fn collect_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        self,
        expression: AtomView<'_>,
        map: F,
    ) -> Atom {
        expression
            .collect_with_map(|a| self.matches(a))
            .into_inner()
            .map_collects(map)
            .unwrap_collect()
    }

    fn expand_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        self,
        expression: AtomView<'_>,
        map: F,
    ) -> Atom {
        expression
            .expand_with_map(|a| self.matches(a))
            .map_collects(map)
            .unwrap_collect()
    }

    /// Match a complete tensor leaf or an opaque scalar power scope.
    pub fn matches(self, arg: AtomView<'_>) -> bool {
        if let AtomView::Pow(power) = arg {
            let (base, exponent) = power.get_base_exp();
            // Positive integer powers belong to the collector's copy frontier.
            // Other powers can be opaque scalar leaves in a parsed network;
            // select the whole scope rather than hiding its domain in a weight.
            if let AtomView::Num(number) = exponent {
                let coefficient = number.get_coeff_view();
                if matches!(
                    coefficient,
                    CoefficientView::Natural(..) | CoefficientView::Large(..)
                ) && coefficient.is_integer()
                    && !coefficient.to_owned().is_negative()
                {
                    return false;
                }
            }
            let mut selected = false;
            base.visitor(&mut |node| {
                if selected {
                    return false;
                }
                if let AtomView::Fun(_) = node {
                    selected = self.matches(node);
                    return false;
                }
                matches!(node, AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_))
            });
            return selected;
        }
        let AtomView::Fun(fun) = arg else {
            return false;
        };

        match self {
            Self::Metrics => fun.get_symbol() == ETS.metric,
            Self::ChainsAndTraces => {
                fun.get_symbol() == SPENSO_TAG.chain || fun.get_symbol() == SPENSO_TAG.trace
            }
            Self::Tensors => TensorialSyntax::function_is_tensorial(
                fun,
                StrictTensorFilter::ContainsReps,
                &SlotMatcher::default(),
            ),
            Self::TaggedTensors => TensorialSyntax::function_is_tensorial(
                fun,
                StrictTensorFilter::Tagged,
                &SlotMatcher::default(),
            ),
            Self::Reps(rep) => Self::function_contains_rep(fun, &rep),
        }
    }

    pub(crate) fn function_contains_rep(fun: FunView<'_>, reps: &[LibraryRep]) -> bool {
        let symbol = fun.get_symbol();

        if symbol == SPENSO_TAG.pure_scalar {
            return false;
        }

        if symbol == SPENSO_TAG.bracket {
            return false;
        }

        if symbol.has_tag(&SPENSO_TAG.broadcast) {
            let args = fun.iter().collect::<Vec<_>>();
            return matches!(args.as_slice(), [AtomView::Fun(arg)] if Self::function_contains_rep(*arg, reps));
        }

        // Match rep(a__) against the whole argument, as matches() does through
        // partial(false). Only the recognized shorthand branches below recurse.
        if fun.iter().any(|arg| {
            matches!(arg, AtomView::Fun(inner)
                if inner.get_nargs() > 0
                    && reps.iter().any(|rep| inner.get_symbol() == rep.symbol()))
        }) {
            return true;
        }

        // Compact products retain their representation on the two vector
        // operands, rather than directly on the metric/dot arguments. Keep
        // this descent restricted to the registered product owners.
        if symbol == SPENSO_TAG.dot || symbol == ETS.metric {
            return fun.iter().any(|argument| {
                matches!(argument, AtomView::Fun(vector)
                    if Self::function_contains_rep(vector, reps))
            });
        }

        // A scalar function owns its interface, but a requested domain can
        // occur in its metadata. Select the whole function so the local rewrite
        // and its normalizer keep their original order and checked boundary.
        if symbol.is_scalar() {
            return fun.iter().any(|argument| {
                let mut selected = false;
                argument.visitor(&mut |node| {
                    if selected {
                        return false;
                    }
                    let AtomView::Fun(nested) = node else {
                        return true;
                    };
                    selected = Self::function_contains_rep(nested, reps);
                    !selected
                        && !nested.get_symbol().is_scalar()
                        && nested.get_symbol() != SPENSO_TAG.pure_scalar
                        && nested.get_symbol() != SPENSO_TAG.bracket
                });
                selected
            });
        }

        if symbol == SPENSO_TAG.chain {
            return fun.iter().skip(2).any(
                |arg| matches!(arg, AtomView::Fun(arg) if Self::function_contains_rep(arg, reps)),
            );
        }
        if symbol == SPENSO_TAG.trace {
            return fun.iter().skip(1).any(
                |arg| matches!(arg, AtomView::Fun(arg) if Self::function_contains_rep(arg, reps)),
            );
        }

        // A registered tensor can carry an explicitly scalar metadata field.
        // Its domain still needs local rewriting after an enclosing chain or
        // trace disappears. Keep the entire tensor as the checked callback
        // boundary, and leave unregistered function scopes opaque.
        symbol.has_tag(&SPENSO_TAG.tensor)
            && fun.iter().any(|argument| {
                matches!(argument, AtomView::Fun(metadata)
                    if metadata.get_symbol().is_scalar()
                        && Self::function_contains_rep(metadata, reps))
            })
    }
}

/// Collect tensor factors without full expression expansion.
pub trait TensorCollectExt {
    /// Collect common tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn collect_tensors(&self) -> Atom;

    /// Collect common tensor leaves that contain `rep` as one of their slot representations.
    fn collect_rep(&self, rep: LibraryRep) -> Atom;

    /// Collect common tensor leaves that contain `rep` as one of their slot representations, using a custom map function.
    fn collect_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom;

    /// Collect common tensor leaves that contain any of the `reps` as one of their slot representations.
    fn collect_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom;

    /// Collect common metric tensors.
    fn collect_metrics(&self) -> Atom;

    /// Expand common tensor leaves that contain `rep` as one of their slot representations, using a custom map function.
    fn expand_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom;
}

impl TensorCollectExt for Atom {
    fn collect_tensors(&self) -> Atom {
        self.as_view().collect_tensors()
    }

    fn collect_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom {
        self.as_view().collect_rep_with_map(rep, map)
    }

    fn collect_rep(&self, rep: LibraryRep) -> Atom {
        self.as_view().collect_rep(rep)
    }

    fn collect_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom {
        self.as_view().collect_reps(reps)
    }

    fn collect_metrics(&self) -> Atom {
        self.as_view().collect_metrics()
    }

    fn expand_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom {
        self.as_view().expand_rep_with_map(rep, map)
    }
}

impl TensorCollectExt for AtomView<'_> {
    fn collect_tensors(&self) -> Atom {
        TensorCollectFilter::<0>::Tensors.collect(*self)
    }

    fn collect_rep(&self, rep: LibraryRep) -> Atom {
        TensorCollectFilter::Reps([rep]).collect(*self)
    }

    fn collect_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom {
        TensorCollectFilter::Reps([rep]).collect_with_map(*self, map)
    }

    fn collect_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom {
        TensorCollectFilter::Reps(reps).collect(*self)
    }

    fn collect_metrics(&self) -> Atom {
        TensorCollectFilter::<0>::Metrics.collect(*self)
    }

    fn expand_rep_with_map<F: FnMut(AtomView, &Context, &mut Settable<'_, Atom>)>(
        &self,
        rep: LibraryRep,
        map: F,
    ) -> Atom {
        TensorCollectFilter::Reps([rep]).expand_with_map(*self, map)
    }
}

#[cfg(test)]
mod scalar_domain_tests {
    use super::*;
    use crate::structure::representation::{Euclidean, Minkowski};

    #[test]
    fn representation_filter_selects_scalar_metadata_without_crossing_explicit_boundaries() {
        crate::structure::representation::initialize();
        let rep = LibraryRep::from(Minkowski {});
        let other = LibraryRep::from(Euclidean {});
        let scalar = symbol!("collect_domain_scalar"; Scalar);
        let nested = symbol!("collect_domain_nested");
        let leaf = FunctionBuilder::new(symbol!("collect_domain_leaf"))
            .add_arg(FunctionBuilder::new(rep.symbol()).add_arg(4).finish())
            .finish();
        let scalar_leaf = FunctionBuilder::new(scalar).add_arg(&leaf).finish();
        let nested_sum = FunctionBuilder::new(scalar)
            .add_arg(FunctionBuilder::new(nested).add_arg(&leaf).finish() + Atom::one())
            .finish();
        for input in [&scalar_leaf, &nested_sum] {
            assert!(TensorCollectFilter::Reps([rep]).matches(input.as_view()));
            assert!(!TensorCollectFilter::Reps([other]).matches(input.as_view()));
        }
        for boundary in [SPENSO_TAG.pure_scalar, SPENSO_TAG.bracket] {
            let wrapped = FunctionBuilder::new(boundary).add_arg(&leaf).finish();
            assert!(!TensorCollectFilter::Reps([rep]).matches(wrapped.as_view()));
            let input = FunctionBuilder::new(scalar).add_arg(wrapped).finish();
            assert!(!TensorCollectFilter::Reps([rep]).matches(input.as_view()));
        }
        let ordinary = FunctionBuilder::new(nested).add_arg(scalar_leaf).finish();
        assert!(!TensorCollectFilter::Reps([rep]).matches(ordinary.as_view()));
    }

    #[test]
    fn representation_filter_selects_scalar_metadata_on_registered_tensor_leaves() {
        crate::structure::representation::initialize();
        let mink = LibraryRep::from(Minkowski {});
        let euc = LibraryRep::from(Euclidean {});
        let scalar = symbol!("collect_tensor_metadata_scalar"; Scalar);
        let inner = FunctionBuilder::new(
            SPENSO_TAG.rank_one_tensor_symbol("collect_tensor_metadata_inner"),
        )
        .add_arg(crate::mink!(4))
        .finish();
        let outer = SPENSO_TAG.rank_one_tensor_symbol("collect_tensor_metadata_outer");
        let representation = FunctionBuilder::new(euc.symbol()).add_arg(4).finish();
        let make_outer = |metadata: Atom| {
            FunctionBuilder::new(outer)
                .add_arg(FunctionBuilder::new(scalar).add_arg(metadata).finish())
                .add_arg(&representation)
                .finish()
        };
        let tensor = make_outer(inner.clone());
        let filter = TensorCollectFilter::Reps([mink]);
        assert!(filter.matches(tensor.as_view()));
        let compact = FunctionBuilder::new(outer)
            .add_arg(&representation)
            .finish();
        let product = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(&tensor)
            .add_arg(&compact)
            .finish();
        assert!(filter.matches(product.as_view()));
        assert!(!filter.matches(compact.as_view()));
        for boundary in [SPENSO_TAG.pure_scalar, SPENSO_TAG.bracket] {
            let hidden = FunctionBuilder::new(boundary).add_arg(&inner).finish();
            assert!(!filter.matches(make_outer(hidden).as_view()));
        }
    }

    #[test]
    fn representation_filter_selects_compact_products_without_opening_foreign_functions() {
        crate::structure::representation::initialize();
        let mink = LibraryRep::from(Minkowski {});
        let euclidean = LibraryRep::from(Euclidean {});
        let p = SPENSO_TAG.rank_one_tensor_symbol("collect_compact_product_p");
        let q = SPENSO_TAG.rank_one_tensor_symbol("collect_compact_product_q");
        let opaque = symbol!("collect_compact_product_opaque");
        for rep in [mink, euclidean] {
            let representation = FunctionBuilder::new(rep.symbol()).add_arg(4).finish();
            let p = FunctionBuilder::new(p).add_arg(&representation).finish();
            let q = FunctionBuilder::new(q).add_arg(&representation).finish();
            for owner in [SPENSO_TAG.dot, ETS.metric] {
                let product = FunctionBuilder::new(owner).add_arg(&p).add_arg(&q).finish();
                assert!(TensorCollectFilter::Reps([rep]).matches(product.as_view()));
                let other = if rep == mink { euclidean } else { mink };
                assert!(!TensorCollectFilter::Reps([other]).matches(product.as_view()));
                for boundary in [opaque, SPENSO_TAG.pure_scalar, SPENSO_TAG.bracket] {
                    let wrapped = FunctionBuilder::new(boundary).add_arg(&product).finish();
                    assert!(!TensorCollectFilter::Reps([rep]).matches(wrapped.as_view()));
                }
            }
        }
    }

    #[test]
    fn filtered_collection_retains_aliases_until_explicit_resolution() {
        crate::structure::representation::initialize();
        let tensor =
            FunctionBuilder::new(SPENSO_TAG.rank_one_tensor_symbol("collect_alias_tensor"))
                .add_arg(crate::mink!(4, 99101))
                .finish();
        let coefficient = (Atom::one() + Atom::var(symbol!("collect_alias_x"))).pow(9);
        let source = &coefficient * &tensor;
        let collected = source
            .as_view()
            .collect_with_map(|part| part == tensor.as_view());
        assert!(!collected.get_aliases().is_empty());
        assert_ne!(collected.get_root(), &source);
        let mut observed = Vec::new();
        let resolved = collected
            .into_inner()
            .map_collects(|wrapped, _, _| {
                observed.push(
                    wrapped
                        .as_fun_view()
                        .unwrap()
                        .iter()
                        .next()
                        .unwrap()
                        .to_owned(),
                );
            })
            .unwrap_collect();
        assert_eq!(observed, vec![tensor]);
        assert_eq!(resolved, source);
        let unselected = source.as_view().collect_with_map(|_| false);
        assert!(unselected.get_aliases().is_empty());
        assert_eq!(unselected.into_inner(), source);
    }
    #[test]
    fn representation_filter_keeps_nonpositive_and_fractional_power_scopes() {
        crate::structure::representation::initialize();
        let rep = LibraryRep::from(Minkowski {});
        let vector = FunctionBuilder::new(SPENSO_TAG.rank_one_tensor_symbol("collect_power_p"))
            .add_arg(crate::mink!(4))
            .finish();
        let dot = FunctionBuilder::new(SPENSO_TAG.dot)
            .add_arg(&vector)
            .add_arg(&vector)
            .finish();
        let filter = TensorCollectFilter::Reps([rep]);
        for exponent in [Atom::num(-2), Atom::num((1, 2)), Atom::num((-1, 2))] {
            assert!(filter.matches((Atom::one() + &dot).pow(&exponent).as_view()));
            for boundary in [
                symbol!("collect_power_foreign"),
                SPENSO_TAG.pure_scalar,
                SPENSO_TAG.bracket,
            ] {
                let hidden = FunctionBuilder::new(boundary).add_arg(&dot).finish();
                assert!(!filter.matches((Atom::one() + hidden).pow(&exponent).as_view()));
            }
        }
        assert!(!filter.matches((Atom::one() + dot).pow(2).as_view()));
    }
}
