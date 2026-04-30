use crate::rep_symbols::RS;

use super::*;
use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use symbolica::atom::Atom;

pub struct MetricSimplify;

impl<Aind: AbsInd + ParseableAind> Contract<SymbolicTensor<Aind>, MetricSimplify>
    for SymbolicTensor<Aind>
{
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
        println!("HI");
        // println!(
        //     "Contracting\n{}\nwith\n{}",
        //     self.expression, other.expression
        // );
        // let self_net = self
        //     .expression
        //     .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
        //         depth_limit: Some(1),
        //         take_first_term_from_sum: false,
        //         ..Default::default()
        //     })
        //     .unwrap();
        // let other_net = other
        //     .expression
        //     .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
        //         depth_limit: Some(1),
        //         take_first_term_from_sum: false,
        //         ..Default::default()
        //     })
        //     .unwrap();

        let mut self_res = self.expression.clone(); //.simple_execute::<MetricSimplify>().expand();
        let mut other_res = other.expression.clone(); //other_net.simple_execute::<MetricSimplify>().expand();

        let (structure, self_pos, _, _) = self.structure.merge(&other.structure)?;

        let mut to_rep = vec![];
        for i in self_pos.included_iter() {
            let rep = self.structure.get_rep(i).unwrap();
            let slot = self.structure.get_slot(i).unwrap();
            let slot_atom = slot.to_atom();

            let lhs = function!(symbol!("f_"), &slot_atom).to_pattern();
            let rhs = function!(
                ETS.metric,
                slot_atom,
                function!(symbol!("f_"), rep.to_symbolic([]))
            )
            .to_pattern();
            // println!("Replacing {} with {}", lhs, rhs);
            self_res = self_res.replace(&lhs).with(&rhs);
            other_res = other_res.replace(&lhs).with(&rhs);

            // a.pattern(symbol!("a_"));
            to_rep.push((
                function!(ETS.metric, slot_atom, symbol!("a_")).to_pattern(),
                symbol!("a_"),
                self.structure.get_slot(i).unwrap(),
            ));
        }

        // let prod = (self_res * other_res).expand();

        // for (pat, ind, matched_slot) in &to_rep {
        //     prod.replace(
        //         (function!(ETS.metric, matched_slot.to_atom(), symbol!("a_")) * symbol!("x___"))
        //             .to_pattern(),
        //     )
        //     .level_range((0, Some(0)))
        //     .repeat()
        //     .with_map(move |m| {
        //         // let a1 = m.get(muw1_).unwrap().to_atom();
        //         // let a2 = m.get(k1_).unwrap().to_atom();
        //         // let dest = m.get(x___).unwrap().to_atom(); // PREVENT!

        //         // dest.replace(a1)
        //         //     .level_range((1, Some(1)))
        //         //     .rhs_cache_size(0)
        //         //     .with(a2)

        //         // if let Some(Match::Single(a1)) = m.get(muw1_) {
        //         //     if let Some(Match::Single(a2)) = m.get(k1_) {
        //         //         l

        //         //         return dest.replace_map(|input, _, out| {
        //         //             if input == *a1 {
        //         //                 out.set_from_view(a2);
        //         //             }
        //         //         });
        //         //     }
        //         // }

        //         // Atom::Zero
        //     });
        // }
        let mut expression = Atom::Zero;

        for t in self_res.expand().terms() {
            let mut t = t.to_owned();
            let mut other = other_res.clone();
            // println!("Processing term: {}", t);
            // let mut t = t.to_owned();
            for (pat, ind, matched_slot) in &to_rep {
                // println!("For slot {}", matched_slot);
                let slot_atom = matched_slot.to_atom();

                let Some(m) = t.replace(pat).once().match_iter().next() else {
                    continue;
                };

                let other_pat = slot_atom.to_pattern();
                let other_rep = m[ind].to_pattern();

                println!("{}->{}", other_pat, other_rep);

                other = other.replace(other_pat).with(other_rep);

                let self_pat = pat.replace_wildcards(&m).to_pattern();
                println!("{}->1", self_pat);
                t = t.replace(self_pat).with(Atom::one());
            }

            // for rep in &final_rep {
            //     println!("Final rep: {}", rep);
            // }

            // for rep in &initial_rep {
            //     println!("Initial rep: {}", rep);
            // }

            let other_net = other
                .expand()
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                    depth_limit: Some(2),
                    take_first_term_from_sum: false,
                    ..Default::default()
                })
                .unwrap();
            let mut metric_rep = t * other_net.simple_execute::<MetricSimplify>();

            for (pat, ind, matched_slot) in &to_rep {
                println!("For slot {}", matched_slot);
                let Some(m) = metric_rep
                    .replace(matched_slot.to_atom())
                    .once()
                    .match_iter()
                    .next()
                else {
                    continue;
                };

                panic!("Found here:{:>}", metric_rep.expand());
            }
            // .replace_multiple(&final_rep);
            //
            // println!("Metric rep: {}", metric_rep);

            expression += metric_rep;

            // t.replace_map(||)
        }

        Ok(SymbolicTensor {
            expression: expression
                .replace(function!(
                    ETS.metric,
                    SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_]),
                    SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_])
                ))
                .repeat()
                .level_range((0, Some(0)))
                .with(Atom::var(RS.d_))
                .replace(
                    function!(
                        ETS.metric,
                        SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.i_]),
                        SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.j_])
                    )
                    .npow(2),
                )
                .repeat()
                .level_range((0, Some(0)))
                .with(Atom::var(RS.d_))
                .replace(
                    function!(
                        ETS.metric,
                        RS.a_,
                        SPENSO_TAG.self_dual_::<0, _>([RS.d_, RS.j_])
                    )
                    .npow(2),
                )
                .repeat()
                .level_range((0, Some(0)))
                .with(function!(ETS.metric, RS.a_, RS.a_)),
            is_composite: true,
            is_metric: false,
            structure,
        })
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use spenso::network::library::symbolic::ETS;
    use symbolica::parse;

    #[test]
    fn parse() {
        let _ = ETS.metric;
        let expr = parse!("g(mink(4,6),mink(4,7))");

        let structure = SymbolicTensor::from_permuted(
            &PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(expr).unwrap(),
        )
        .unwrap();

        Contract::<SymbolicTensor, ()>::contract(&structure, &structure).unwrap();
    }

    #[test]
    fn simplify_metrics() {
        let _ = ETS.metric;
        let expr = parse!("k(spenso::mink(d,7))*p(spenso::mink(d,7))");

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                depth_limit: Some(1),
                take_first_term_from_sum: false,
                ..Default::default()
            })
            .unwrap();

        println!("{}", net.dot_pretty());

        println!("{}", net.simple_execute::<MetricSimplify>())
    }
}
