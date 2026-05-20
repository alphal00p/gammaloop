use std::sync::LazyLock;

use crate::{
    network::{
        library::symbolic::ETS,
        parsing::{StrictTensorFilter, structure_inference::TensorialSyntax},
        tags::SPENSO_TAG,
    },
    shadowing::{static_symbols::W_, symbolica_utils::ReplaceBuilderExt},
    structure::representation::LibraryRep,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, representation::FunView},
    function, symbol,
};

pub static COLLECT: LazyLock<Symbol> = LazyLock::new(|| symbol!("spenso::collect"));

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
    fn collect_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom;
    fn expand_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom;
    fn wrap_in_collect(self) -> Atom;
    fn unwrap_collect(self) -> Atom;
}
impl Collectable for Atom {
    fn collect_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom {
        self.as_view().collect_with_map(map)
    }

    fn expand_with_map(self, map: impl FnMut(AtomView<'_>) -> bool) -> Atom {
        self.as_view().expand_with_map(map)
    }

    fn unwrap_collect(self) -> Atom {
        self.as_view().unwrap_collect()
    }

    fn wrap_in_collect(self) -> Atom {
        self.as_view().wrap_in_collect()
    }
}
impl Collectable for AtomView<'_> {
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
        let collected = wrapped.expand_in(*COLLECT);
        collected.unwrap_collect()
    }
    fn collect_with_map(self, mut matches: impl FnMut(AtomView<'_>) -> bool) -> Atom {
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
        let collected = wrapped.collect_symbol::<i16>(*COLLECT, None, None);
        collected.unwrap_collect()
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
            }
        })
    }

    fn wrap_in_collect(self) -> Atom {
        FunctionBuilder::new(*COLLECT).add_arg(self).finish()
    }
}

impl<const N: usize> TensorCollectFilter<N> {
    fn collect(self, expression: AtomView<'_>) -> Atom {
        expression.collect_with_map(|a| self.matches(a))
    }

    fn expand(self, expression: AtomView<'_>) -> Atom {
        expression.expand_with_map(|a| self.matches(a))
    }

    fn matches(self, arg: AtomView<'_>) -> bool {
        let AtomView::Fun(fun) = arg else {
            return false;
        };

        match self {
            Self::Metrics => fun.get_symbol() == ETS.metric,
            Self::ChainsAndTraces => {
                fun.get_symbol() == SPENSO_TAG.chain || fun.get_symbol() == SPENSO_TAG.trace
            }
            Self::Tensors => {
                TensorialSyntax::function_is_tensorial(fun, StrictTensorFilter::ContainsReps)
            }
            Self::TaggedTensors => {
                TensorialSyntax::function_is_tensorial(fun, StrictTensorFilter::Tagged)
            }
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

        for a in fun.iter() {
            for r in reps {
                if a.replace(function!(r.symbol(), W_.a__)).matches() {
                    return true;
                }
            }
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

        false
    }
}

/// Collect tensor factors without full expression expansion.
pub trait TensorCollectExt {
    /// Collect common tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn collect_tensors(&self) -> Atom;

    /// Collect common tagged tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn collect_tagged_tensors(&self) -> Atom;

    /// Collect common tensor leaves that contain `rep` as one of their slot representations.
    fn collect_rep(&self, rep: LibraryRep) -> Atom;

    /// Collect common tensor leaves that contain any of the `reps` as one of their slot representations.
    fn collect_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom;

    /// Collect common metric tensors.
    fn collect_metrics(&self) -> Atom;

    /// Collect common chain and trace tensors.
    fn collect_chains_and_traces(&self) -> Atom;

    /// Expand common tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn expand_tensors(&self) -> Atom;

    /// Expand common tagged tensor leaves by temporarily wrapping them in `spenso::collect(...)`.
    fn expand_tagged_tensors(&self) -> Atom;

    /// Expand common tensor leaves that contain `rep` as one of their slot representations.
    fn expand_rep(&self, rep: LibraryRep) -> Atom;

    /// Expand common tensor leaves that contain any of the `reps` as one of their slot representations.
    fn expand_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom;

    /// Expand common metric tensors.
    fn expand_metrics(&self) -> Atom;

    /// Expand common chain and trace tensors.
    fn expand_chains_and_traces(&self) -> Atom;
}

impl TensorCollectExt for Atom {
    fn collect_tensors(&self) -> Atom {
        self.as_view().collect_tensors()
    }

    fn collect_tagged_tensors(&self) -> Atom {
        self.as_view().collect_tagged_tensors()
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

    fn collect_chains_and_traces(&self) -> Atom {
        self.as_view().collect_chains_and_traces()
    }

    fn expand_tensors(&self) -> Atom {
        self.as_view().expand_tensors()
    }

    fn expand_tagged_tensors(&self) -> Atom {
        self.as_view().expand_tagged_tensors()
    }

    fn expand_rep(&self, rep: LibraryRep) -> Atom {
        self.as_view().expand_rep(rep)
    }

    fn expand_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom {
        self.as_view().expand_reps(reps)
    }

    fn expand_metrics(&self) -> Atom {
        self.as_view().expand_metrics()
    }

    fn expand_chains_and_traces(&self) -> Atom {
        self.as_view().expand_chains_and_traces()
    }
}

impl TensorCollectExt for AtomView<'_> {
    fn collect_chains_and_traces(&self) -> Atom {
        TensorCollectFilter::<0>::ChainsAndTraces.collect(*self)
    }
    fn collect_tensors(&self) -> Atom {
        TensorCollectFilter::<0>::Tensors.collect(*self)
    }

    fn collect_tagged_tensors(&self) -> Atom {
        TensorCollectFilter::<0>::TaggedTensors.collect(*self)
    }

    fn collect_rep(&self, rep: LibraryRep) -> Atom {
        TensorCollectFilter::Reps([rep]).collect(*self)
    }

    fn collect_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom {
        TensorCollectFilter::Reps(reps).collect(*self)
    }

    fn collect_metrics(&self) -> Atom {
        TensorCollectFilter::<0>::Metrics.collect(*self)
    }

    fn expand_chains_and_traces(&self) -> Atom {
        TensorCollectFilter::<0>::ChainsAndTraces.expand(*self)
    }
    fn expand_tensors(&self) -> Atom {
        TensorCollectFilter::<0>::Tensors.expand(*self)
    }

    fn expand_tagged_tensors(&self) -> Atom {
        TensorCollectFilter::<0>::TaggedTensors.expand(*self)
    }

    fn expand_rep(&self, rep: LibraryRep) -> Atom {
        TensorCollectFilter::Reps([rep]).expand(*self)
    }

    fn expand_reps<const N: usize>(&self, reps: [LibraryRep; N]) -> Atom {
        TensorCollectFilter::Reps(reps).expand(*self)
    }

    fn expand_metrics(&self) -> Atom {
        TensorCollectFilter::<0>::Metrics.expand(*self)
    }
}
