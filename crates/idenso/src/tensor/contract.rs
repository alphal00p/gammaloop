use super::*;
use crate::metric::MetricSimplifier;
use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};
use spenso::structure::SlotIndex;

pub struct MetricSimplify;

impl<Aind: AbsInd + ParseableAind> SymbolicTensor<Aind> {
    fn single_contracted_pos(positions: &SubSet<SlotIndex>) -> Option<SlotIndex> {
        (positions.n_included() == 1)
            .then(|| positions.included_iter().next())
            .flatten()
    }

    fn metric_free_slot(&self, contracted_pos: SlotIndex) -> Option<LibrarySlot<Aind>> {
        if self.structure.order() != 2 {
            return None;
        }

        (0..self.structure.order())
            .map(SlotIndex::from)
            .find(|&pos| pos != contracted_pos)
            .and_then(|pos| self.structure.get_slot(pos))
    }

    fn contract_metric_into(
        metric: &Self,
        tensor: &Self,
        metric_positions: &SubSet<SlotIndex>,
        tensor_positions: &SubSet<SlotIndex>,
        structure: OrderedStructure<LibraryRep, Aind>,
    ) -> Option<Self> {
        if metric.is_composite || !metric.is_metric {
            return None;
        }

        let metric_pos = Self::single_contracted_pos(metric_positions)?;
        let tensor_pos = Self::single_contracted_pos(tensor_positions)?;
        let free_metric_slot = metric.metric_free_slot(metric_pos)?;
        let contracted_tensor_slot = tensor.structure.get_slot(tensor_pos)?;

        Some(Self {
            structure,
            is_composite: tensor.is_composite,
            is_metric: tensor.is_metric,
            expression: tensor
                .expression
                .replace(contracted_tensor_slot.to_atom())
                .with(free_metric_slot.to_atom())
                .simplify_metrics(),
        })
    }
}

impl<Aind: AbsInd + ParseableAind> Contract<SymbolicTensor<Aind>, MetricSimplify>
    for SymbolicTensor<Aind>
{
    type LCM = SymbolicTensor<Aind>;
    fn contract(&self, other: &SymbolicTensor<Aind>) -> Result<Self::LCM, ContractionError> {
        let (structure, self_pos, other_pos, _) = self.structure.merge(&other.structure)?;

        if let Some(result) =
            Self::contract_metric_into(self, other, &self_pos, &other_pos, structure.clone())
        {
            return Ok(result);
        }

        if let Some(result) =
            Self::contract_metric_into(other, self, &other_pos, &self_pos, structure.clone())
        {
            return Ok(result);
        }

        Ok(Self {
            expression: (&other.expression * &self.expression).simplify_metrics(),
            is_composite: true,
            is_metric: false,
            structure,
        })
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    use crate::representations::initialize;
    use spenso::network::library::symbolic::ETS;
    use spenso::network::parsing::StructureFromAtom;
    use symbolica::{parse, symbol};

    #[test]
    fn parse() {
        let _ = ETS.metric;
        let expr = parse!("g(mink(4,6),mink(4,7))");

        let parsed = ShadowedStructure::<AbstractIndex>::parse(expr.as_view()).unwrap();
        let structure = SymbolicTensor::from_permuted(&parsed).unwrap();

        Contract::<SymbolicTensor, ()>::contract(&structure, &structure).unwrap();
    }

    #[test]
    fn simplify_metrics() {
        initialize();
        symbol!("metric_simplify_mu", "metric_simplify_nu"; tags = ["spenso::index"]);
        let expr = parse!(
            "spenso::g(spenso::mink(4,metric_simplify_mu),spenso::mink(4,metric_simplify_nu))*p(spenso::mink(4,metric_simplify_nu))"
        );

        let net = expr
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                depth_limit: Some(1),
                take_first_term_from_sum: false,
                ..Default::default()
            })
            .unwrap();

        assert_eq!(
            net.simple_execute::<MetricSimplify>(),
            parse!("p(spenso::mink(4,metric_simplify_mu))")
        );
    }
}
