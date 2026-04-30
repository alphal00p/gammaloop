use spenso::{
    contraction::{Contract, ContractionError},
    network::{library::symbolic::ETS, parsing::ParseSettings, tags::SPENSO_TAG as T},
    shadowing::symbolica_utils::SpensoPrintSettings,
    structure::{
        StructureContract, TensorStructure,
        slot::{AbsInd, DualSlotTo, DummyAind, IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
};

use crate::{
    W_,
    tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor},
};

pub struct Schoonschipify<const EXPANDSUMS: bool, const DEEPEST: bool>;

impl<
    const EXPANDSUMS: bool,
    const DEEPEST: bool,
    Aind: AbsInd + DummyAind + ParseableAind + 'static,
> Contract<SymbolicTensor<Aind>, Schoonschipify<EXPANDSUMS, DEEPEST>> for SymbolicTensor<Aind>
{
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
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

        let (sexpr, oexpr) = if DEEPEST {
            (
                self.expression
                    .schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(&SchoonschipSettings {
                        repeat: true,
                    }),
                other
                    .expression
                    .schoonschip_with_net::<EXPANDSUMS, DEEPEST, Aind>(&SchoonschipSettings {
                        repeat: true,
                    }),
            )
        } else {
            (self.expression.clone(), other.expression.clone())
        };

        if self.structure.is_scalar() || other.structure.is_scalar() {
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return Ok(SymbolicTensor {
                structure,
                is_composite: true,
                is_metric: false,
                expression: &sexpr * &oexpr,
            });
        }

        if !self.is_composite && self.structure.order() == 1 {
            let slot = self.structure.get_slot(0).unwrap();
            let stripped = slot.rep().base().to_symbolic([]);

            let expr = sexpr.replace(slot.to_atom()).with(stripped);
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return Ok(SymbolicTensor {
                structure,
                is_composite: true,
                is_metric: other.is_metric,
                expression: oexpr
                    .replace(slot.dual().to_atom())
                    .with(expr)
                    .normalize_dots(),
            });
        } else if !other.is_composite && other.structure.order() == 1 {
            let slot = other.structure.get_slot(0).unwrap();
            let stripped = slot.rep().base().to_symbolic([]);

            let expr = oexpr.replace(slot.to_atom()).with(stripped);
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            return Ok(SymbolicTensor {
                structure,
                is_composite: true,
                is_metric: self.is_metric,
                expression: sexpr
                    .replace(slot.dual().to_atom())
                    .with(expr)
                    .normalize_dots(),
            });
        } else {
            let expression = &oexpr * &sexpr;
            let (structure, _, _, _) = self.structure.merge(&other.structure)?;

            Ok(Self {
                structure,
                is_composite: true,
                is_metric: false,
                expression,
            })
        }
    }
}

pub struct SchoonschipSettings {
    repeat: bool,
}

impl Default for SchoonschipSettings {
    fn default() -> Self {
        Self { repeat: false }
    }
}
pub trait Schoonschip {
    fn schoonschip(&self) -> Atom;

    fn normalize_dots(&self) -> Atom;
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

    fn schoonschip_with_net<
        const EXPANDSUMS: bool,
        const DEEPEST: bool,
        Aind: AbsInd + DummyAind + ParseableAind + 'static,
    >(
        &self,
        settings: &SchoonschipSettings,
    ) -> Atom {
        let mut res: Option<Atom> = None;
        while {
            let (mut cont, new) = if let Some(old) = res {
                let net = old
                    .as_view()
                    .parse_to_symbolic_net::<Aind>(&ParseSettings {
                        depth_limit: Some(1),
                        take_first_term_from_sum: false,
                        ..Default::default()
                    })
                    .unwrap();
                let new = net.simple_execute::<Schoonschipify<EXPANDSUMS, DEEPEST>>();
                (old.as_view() != new.as_view(), new)
            } else {
                let net = (*self)
                    .parse_to_symbolic_net::<Aind>(&ParseSettings {
                        depth_limit: Some(1),
                        take_first_term_from_sum: false,
                        ..Default::default()
                    })
                    .unwrap();
                let new = net.simple_execute::<Schoonschipify<EXPANDSUMS, DEEPEST>>();
                (*self != new.as_view(), new)
            };

            println!(
                "New: {}",
                new.printer(SpensoPrintSettings::compact().nice_symbolica())
            );
            res = Some(new);
            cont &= settings.repeat;
            cont
        } {}
        res.unwrap()
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
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings {
                repeat: true,
            });
        assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),P(2,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = (&p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2)))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings {
                repeat: true,
            });
        assert_snapshot!(result.to_bare_ordered_string(),@"(P(2,mink(D,2))*Q(2,bla,mink(D,2))+Q(2,bla,mink(D,2))*Q(3,mink(D,2)))*g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

        let result = ((p1 + q3 * p1_2 * q2_2) * (q2 + p2))
            .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings {
                repeat: true,
            });
        assert_snapshot!(result.to_bare_ordered_string(),@"(P(1,mink(D,1))+Q(3,mink(D,1))*g(P(1,mink(D)),Q(2,bla,mink(D))))*(P(2,mink(D,1))+Q(2,bla,mink(D,1)))");
    }
}
