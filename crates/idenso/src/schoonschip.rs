use linnet::half_edge::subgraph::SubSetLike;
use spenso::{
    algebra::ScalarMul,
    contraction::{Contract, ContractionError},
    network::{
        ContractScalars, ContractionStrategy, ExecutionResult, Sequential, SingleSmallestDegree,
        TensorNetworkError, TensorOrScalarOrKey,
        graph::NetworkGraph,
        library::{DummyKey, DummyLibrary, function_lib::Wrap, symbolic::ETS},
        parsing::ParseSettings,
        store::NetworkStore,
        tags::SPENSO_TAG as T,
    },
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        StructureContract, TensorStructure,
        permuted::PermuteTensor,
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
};

use crate::{
    W_,
    tensor::{SymbolicNetParse, SymbolicTensor},
};

const TRACE_SCHOONSCHIP: bool = false;

pub struct Schoonschipify<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool>;

fn is_sum(expr: &Atom) -> bool {
    matches!(expr.as_view(), AtomView::Add(_))
}

struct SchoonschipSmallestDegree<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
>;

impl<const EXPANDSUMS: bool, const RECURSE: bool, const DEPTH_FIRST: bool>
    SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
{
    fn simplify_scalar_tensors<Aind: AbsInd + DummyAind + ParseableAind + 'static>(
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
    ) {
        if !RECURSE {
            return;
        }

        let settings = if DEPTH_FIRST {
            SchoonschipSettings::depth_first(Some(1)).without_parse_inner_products()
        } else {
            SchoonschipSettings::breadth_first(Some(1)).without_parse_inner_products()
        };

        for tensor in &mut executor.tensors {
            if tensor.structure.is_scalar() && tensor.is_composite {
                tensor.expression = tensor
                    .expression
                    .schoonschip_with_net::<EXPANDSUMS, true, Aind>(&settings);
                tensor.is_composite = false;
                tensor.is_metric = false;
            }
        }
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
>
    ContractionStrategy<
        NetworkStore<SymbolicTensor<Aind>, Atom>,
        DummyLibrary<SymbolicTensor<Aind>>,
        DummyKey,
        symbolica::atom::Symbol,
        Aind,
    > for SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>
where
    SymbolicTensor<Aind>: ScalarMul<Atom, Output = SymbolicTensor<Aind>>
        + PermuteTensor<Permuted = SymbolicTensor<Aind>>,
{
    fn contract(
        executor: &mut NetworkStore<SymbolicTensor<Aind>, Atom>,
        graph: NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>,
        lib: &DummyLibrary<SymbolicTensor<Aind>>,
    ) -> Result<
        (NetworkGraph<DummyKey, symbolica::atom::Symbol, Aind>, bool),
        TensorNetworkError<DummyKey, symbolica::atom::Symbol>,
    > {
        Self::simplify_scalar_tensors(executor);
        let (mut graph, mut didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        while {
            let (newgraph, smth) = SingleSmallestDegree::<
                false,
                Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
            >::contract(executor, graph, lib)?;
            graph = newgraph;
            smth
        } {
            didsmth = true;
        }

        Self::simplify_scalar_tensors(executor);
        let (graph, scalar_didsmth) = ContractScalars::<
            Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
        >::contract(executor, graph, lib)?;

        Ok((graph, didsmth || scalar_didsmth))
    }
}

impl<
    const EXPANDSUMS: bool,
    const RECURSE: bool,
    const DEPTH_FIRST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
> Contract<SymbolicTensor<Aind>, Schoonschipify<EXPANDSUMS, RECURSE, DEPTH_FIRST>>
    for SymbolicTensor<Aind>
{
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
        if TRACE_SCHOONSCHIP {
            println!(
                "Contracting  {} {}rank {} with rank {} {} {}: \n{}\nwith\n{}\n gives:",
                if self.is_composite { "composite " } else { "" },
                if self.is_metric { "metric " } else { "" },
                self.structure.order(),
                if other.is_composite { "composite " } else { "" },
                if other.is_metric { "metric " } else { "" },
                other.structure.order(),
                self.expression,
                other.expression
            );
        }

        let recursive_settings = || {
            if DEPTH_FIRST {
                SchoonschipSettings::depth_first(Some(1)).without_parse_inner_products()
            } else {
                SchoonschipSettings::breadth_first(Some(1)).without_parse_inner_products()
            }
        };

        let recursive_schoonschip = |expr: &Atom| {
            expr.schoonschip_with_net::<EXPANDSUMS, true, Aind>(&recursive_settings())
        };

        let (sexpr, oexpr) = if RECURSE && DEPTH_FIRST {
            (
                recursive_schoonschip(&self.expression),
                recursive_schoonschip(&other.expression),
            )
        } else {
            (self.expression.clone(), other.expression.clone())
        };

        let finish = |mut result: SymbolicTensor<Aind>, recurse_result: bool| {
            if recurse_result {
                result.expression = recursive_schoonschip(&result.expression);
            }
            Ok(result)
        };

        if self.structure.is_scalar() || other.structure.is_scalar() {
            let (sexpr, oexpr) = if RECURSE && !DEPTH_FIRST {
                (
                    recursive_schoonschip(&self.expression),
                    recursive_schoonschip(&other.expression),
                )
            } else {
                (sexpr, oexpr)
            };
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return finish(
                SymbolicTensor {
                    structure,
                    is_composite: true,
                    is_metric: false,
                    expression: &sexpr * &oexpr,
                },
                false,
            );
        }

        if !self.is_composite && self.structure.order() == 1 {
            let slot = self.structure.get_slot(0).unwrap();
            let stripped = slot.rep().base().to_symbolic([]);

            let expr = sexpr.replace(slot.to_atom()).with(stripped);
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return finish(
                SymbolicTensor {
                    structure,
                    is_composite: true,
                    is_metric: other.is_metric,
                    expression: oexpr
                        .replace(slot.dual().to_atom())
                        .with(expr)
                        .normalize_dots(),
                },
                RECURSE && !DEPTH_FIRST,
            );
        } else if !other.is_composite && other.structure.order() == 1 {
            let slot = other.structure.get_slot(0).unwrap();
            let stripped = slot.rep().base().to_symbolic([]);

            let expr = oexpr.replace(slot.to_atom()).with(stripped);
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return finish(
                SymbolicTensor {
                    structure,
                    is_composite: true,
                    is_metric: self.is_metric,
                    expression: sexpr
                        .replace(slot.dual().to_atom())
                        .with(expr)
                        .normalize_dots(),
                },
                RECURSE && !DEPTH_FIRST,
            );
        } else {
            let expression = &oexpr * &sexpr;
            let (structure, pos_self, _, _) = self.structure.merge(&other.structure)?;
            let expression =
                if EXPANDSUMS && pos_self.n_included() > 0 && (is_sum(&sexpr) || is_sum(&oexpr)) {
                    expression
                        .expand()
                        .schoonschip_with_net::<false, true, Aind>(&recursive_settings())
                } else {
                    expression
                };

            finish(
                Self {
                    structure,
                    is_composite: true,
                    is_metric: false,
                    expression,
                },
                RECURSE && !DEPTH_FIRST,
            )
        }
    }
}

pub struct SchoonschipSettings {
    depth_limit: Option<usize>,
    mode: SchoonschipMode,
    parse_inner_products: bool,
    expand_contracted_sums: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SchoonschipMode {
    SinglePass,
    Recursive(SchoonschipTraversal),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchoonschipTraversal {
    DepthFirst,
    BreadthFirst,
}

impl Default for SchoonschipSettings {
    fn default() -> Self {
        Self::partial()
    }
}

impl SchoonschipSettings {
    pub fn new(depth_limit: Option<usize>) -> Self {
        Self::depth_first(depth_limit)
    }

    pub fn depth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
        }
    }

    pub fn breadth_first(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst),
            parse_inner_products: true,
            expand_contracted_sums: false,
        }
    }

    pub fn single_pass(depth_limit: Option<usize>) -> Self {
        Self {
            depth_limit,
            mode: SchoonschipMode::SinglePass,
            parse_inner_products: true,
            expand_contracted_sums: false,
        }
    }

    pub fn partial() -> Self {
        Self::new(Some(1))
    }

    pub fn with_depth(depth_limit: usize) -> Self {
        Self::new(Some(depth_limit))
    }

    pub fn breadth_first_with_depth(depth_limit: usize) -> Self {
        Self::breadth_first(Some(depth_limit))
    }

    pub fn full() -> Self {
        Self::single_pass(None)
    }

    pub fn with_expanded_contracted_sums(mut self) -> Self {
        self.expand_contracted_sums = true;
        self
    }

    fn without_parse_inner_products(mut self) -> Self {
        self.parse_inner_products = false;
        self
    }
}

pub trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn normalize_dots(&self) -> Atom;
    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom;
}

impl Schoonschip for Atom {
    fn schoonschip(&self) -> Atom {
        self.as_view().schoonschip()
    }

    fn normalize_dots(&self) -> Atom {
        self.as_view().normalize_dots()
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_once_with_net::<EXPANDSUMS, RECURSE, DEPTH_FIRST, Aind>(settings)
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        self.as_view()
            .schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(settings)
    }
}

impl Schoonschip for AtomView<'_> {
    fn normalize_dots(&self) -> Atom {
        let stripped = T.rep_::<0, _>([W_.d_]);
        self.replace(T.rank1_::<0, _>([
            Atom::var(W_.c___),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ]))
        .with(ETS.metric(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &stripped]),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &stripped]),
        ))
    }
    fn schoonschip(&self) -> Atom {
        let index_cond = T.index_fiter(W_.i_);
        let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);

        let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
        let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
        let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
        let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

        // first replace vector(rep(d,i)) * vector(rep(d,i)) with g(vector(rep(d)),vector(rep(d)))
        self.replace(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual]),
        )
        .when(&index_cond)
        .with(ETS.metric(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual_stripped]),
        ))
        .replace(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_dual]),
        )
        .when(&index_cond)
        .with(ETS.metric(
            T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
            T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_stripped]),
        ))
        // Now replace vector(rep(d,i)) * T(..,rep(d,i),..) with T(..,vector(rep(d)),..)
        .replace(
            function!(W_.a_, W_.a___, &self_dual, W_.b___)
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual]),
        )
        .when(&index_cond)
        .repeat()
        .with(function!(
            W_.a_,
            W_.a___,
            T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
            W_.b___
        ))
        .replace(
            function!(W_.a_, W_.a___, &dualizable, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]),
        )
        .when(&index_cond)
        .repeat()
        .with(function!(
            W_.a_,
            W_.a___,
            T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
            W_.b___
        ))
        .replace(
            function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]),
        )
        .repeat()
        .with(function!(
            W_.a_,
            W_.a___,
            T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
            W_.b___
        ))
    }

    fn schoonschip_once_with_net<
        const EXPANDSUMS: bool,
        const RECURSE: bool,
        const DEPTH_FIRST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        let mut net = (*self)
            .parse_to_symbolic_net::<Aind>(&ParseSettings {
                depth_limit: settings.depth_limit,
                take_first_term_from_sum: false,
                parse_inner_products: settings.parse_inner_products,
                parse_composite_scalars_as_tensors: RECURSE,
                ..Default::default()
            })
            .unwrap();
        let lib = DummyLibrary::<SymbolicTensor<Aind>>::new();

        net.execute::<
            Sequential,
            SchoonschipSmallestDegree<EXPANDSUMS, RECURSE, DEPTH_FIRST>,
            _,
            _,
            _,
        >(&lib, &Wrap {})
        .unwrap();

        match net.result().unwrap() {
            ExecutionResult::One => Atom::num(1),
            ExecutionResult::Zero => Atom::Zero,
            ExecutionResult::Val(a) => match a {
                TensorOrScalarOrKey::Key { .. } => panic!("unexpected library key result"),
                TensorOrScalarOrKey::Scalar(s) => s.clone(),
                TensorOrScalarOrKey::Tensor { tensor, .. } => tensor.expression.clone(),
            },
        }
    }

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        let new = if settings.expand_contracted_sums {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<true, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<true, true, false, Aind>(settings)
                }
            }
        } else {
            match (settings.mode, DEEPEST) {
                (SchoonschipMode::SinglePass, _) | (_, false) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, false, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::DepthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, true, Aind>(settings)
                }
                (SchoonschipMode::Recursive(SchoonschipTraversal::BreadthFirst), true) => {
                    self.schoonschip_once_with_net::<EXPANDSUMS, true, false, Aind>(settings)
                }
            }
        };

        if TRACE_SCHOONSCHIP {
            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
        }

        new
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;
    use spenso::{
        shadowing::symbolica_utils::AtomCoreExt,
        structure::{
            abstract_index::AbstractIndex,
            representation::{Minkowski, RepName, Representation},
            slot::IsAbstractSlot,
        },
    };
    use symbolica::symbol;

    use super::*;

    #[test]
    fn simple_dot() {
        let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
        let p = T.rank_one_tensor_symbol("P");
        let q = T.rank_one_tensor_symbol("Q");

        let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p1_2 = function!(p, 1, mink.slot::<AbstractIndex, _>(2).to_atom());
        let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());

        let q2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(1).to_atom()
        );
        let q2_2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(2).to_atom()
        );

        let q3 = function!(q, 3, mink.slot::<AbstractIndex, _>(1).to_atom());
        let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

        let result = (&p1 * &q2).schoonschip();
        assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = (&p1 * (&q2 + &p2 * &q3_2 * &q2_2))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial());
        assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),P(2,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = (&p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2)))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
        assert_snapshot!(result.to_bare_ordered_string(),@"(g(P(2,mink(D)),Q(2,bla,mink(D)))+g(Q(2,bla,mink(D)),Q(3,mink(D))))*g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let expr = (p1 + q3 * p1_2 * q2_2) * (q2 + p2);

        let result =
            expr.schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
        assert_snapshot!(result.to_bare_ordered_string(),@"(P(1,mink(D,1))+Q(3,mink(D,1))*g(P(1,mink(D)),Q(2,bla,mink(D))))*(P(2,mink(D,1))+Q(2,bla,mink(D,1)))");

        let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
            &SchoonschipSettings::full().with_expanded_contracted_sums(),
        );
        let result = result.to_bare_ordered_string();
        assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
        assert!(!result.contains("mink(D,1)"));

        let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
            &SchoonschipSettings::partial().with_expanded_contracted_sums(),
        );
        let result = result.to_bare_ordered_string();
        assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
        assert!(!result.contains("mink(D,1)"));
    }

    #[test]
    fn benchmark_modes_output() {
        let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
        let p = T.rank_one_tensor_symbol("P");
        let q = T.rank_one_tensor_symbol("Q");

        let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
        let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());

        let q2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(1).to_atom()
        );
        let q2_2 = function!(
            q,
            2,
            symbol!("bla"),
            mink.slot::<AbstractIndex, _>(2).to_atom()
        );
        let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

        let expr = &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2));

        let single_pass_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::single_pass(
                Some(1),
            ))
            .to_bare_ordered_string();
        let depth_first_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial())
            .to_bare_ordered_string();
        let breadth_first_depth_one = expr
            .clone()
            .schoonschip_with_net::<false, true, AbstractIndex>(
                &SchoonschipSettings::breadth_first(Some(1)),
            )
            .to_bare_ordered_string();
        let full_top = expr
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full())
            .to_bare_ordered_string();

        assert_ne!(single_pass_depth_one, full_top);
        assert_eq!(depth_first_depth_one, full_top);
        assert_eq!(breadth_first_depth_one, full_top);
    }
}
