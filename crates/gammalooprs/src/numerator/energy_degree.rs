use std::{
    collections::{BTreeMap, BTreeSet},
    ops::MulAssign,
};

use linnet::half_edge::involution::EdgeIndex;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomView, Symbol};
use thiserror::Error;

use crate::{
    graph::Graph,
    utils::{GS, symbolica_ext::LogPrint},
    uv::UltravioletGraph,
};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct EnergyPowerCapMap {
    degrees: BTreeMap<EdgeIndex, usize>,
}

impl EnergyPowerCapMap {
    pub fn is_empty(&self) -> bool {
        self.degrees.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (EdgeIndex, usize)> + '_ {
        self.degrees.iter().map(|(edge, degree)| (*edge, *degree))
    }

    pub fn into_generation_bounds(self) -> Vec<(usize, usize)> {
        self.degrees
            .into_iter()
            .map(|(edge, degree)| (usize::from(edge), degree))
            .collect()
    }

    fn unit(edge: EdgeIndex) -> Self {
        let mut degrees = BTreeMap::new();
        degrees.insert(edge, 1);
        Self { degrees }
    }

    fn add_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            *self.degrees.entry(edge).or_insert(0) += degree;
        }
    }

    fn max_assign(&mut self, other: Self) {
        for (edge, degree) in other.degrees {
            let slot = self.degrees.entry(edge).or_insert(0);
            *slot = (*slot).max(degree);
        }
    }

    fn scale(&mut self, exponent: usize) {
        for degree in self.degrees.values_mut() {
            degree.mul_assign(exponent);
        }
    }
}

#[derive(Debug, Error)]
pub enum EnergyPowerAnalysisError {
    #[error("loop momentum K({loop_index}, ...) has no carrier edge in this loop momentum basis")]
    MissingLoopCarrierEdge { loop_index: usize },
    #[error(
        "non-polynomial energy dependence encountered in numerator subexpression `{expression}` with exponent `{exponent}`"
    )]
    NonPolynomialEnergyPower {
        expression: String,
        exponent: String,
    },
}

pub struct EnergyPowerAnalyzer {
    loop_edges: Vec<EdgeIndex>,
    internal_edges: Option<BTreeSet<EdgeIndex>>,
    minkowski_symbol: Symbol,
}

impl EnergyPowerAnalyzer {
    pub fn new(loop_edges: impl IntoIterator<Item = EdgeIndex>) -> Self {
        Self {
            loop_edges: loop_edges.into_iter().collect(),
            internal_edges: None,
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn with_internal_edges(
        loop_edges: impl IntoIterator<Item = EdgeIndex>,
        internal_edges: impl IntoIterator<Item = EdgeIndex>,
    ) -> Self {
        Self {
            loop_edges: loop_edges.into_iter().collect(),
            internal_edges: Some(internal_edges.into_iter().collect()),
            minkowski_symbol: LibraryRep::from(Minkowski {}).symbol(),
        }
    }

    pub fn analyze_atom(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(expression.as_view())
    }

    fn analyze_view(
        &self,
        expression: AtomView<'_>,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Num(_) | AtomView::Var(_) => Ok(EnergyPowerCapMap::default()),
            AtomView::Add(add) => {
                let mut max_degree = EnergyPowerCapMap::default();
                for term in add.iter() {
                    max_degree.max_assign(self.analyze_view(term)?);
                }
                Ok(max_degree)
            }
            AtomView::Mul(mul) => {
                let mut product_degree = EnergyPowerCapMap::default();
                for factor in mul.iter() {
                    product_degree.add_assign(self.analyze_view(factor)?);
                }
                Ok(product_degree)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut base_degree = self.analyze_view(base)?;
                if base_degree.is_empty() {
                    return Ok(base_degree);
                }

                let Ok(exponent) = i64::try_from(exponent) else {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                };
                if exponent < 0 {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_string(),
                    });
                }
                base_degree.scale(exponent as usize);
                Ok(base_degree)
            }
            AtomView::Fun(function) => {
                if function.get_symbol() == GS.emr_vec {
                    return Ok(EnergyPowerCapMap::default());
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.emr_mom {
                    let edge = EdgeIndex(usize::try_from(function.get(0)).unwrap_or(usize::MAX));
                    return if self.is_internal_edge(edge) && self.is_energy_index(function.get(1)) {
                        Ok(EnergyPowerCapMap::unit(edge))
                    } else {
                        Ok(EnergyPowerCapMap::default())
                    };
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.loop_mom {
                    let loop_index = usize::try_from(function.get(0)).unwrap_or(usize::MAX);
                    return if self.is_energy_index(function.get(1)) {
                        let edge = self.loop_edges.get(loop_index).copied().ok_or(
                            EnergyPowerAnalysisError::MissingLoopCarrierEdge { loop_index },
                        )?;
                        Ok(EnergyPowerCapMap::unit(edge))
                    } else {
                        Ok(EnergyPowerCapMap::default())
                    };
                }

                let mut degree = EnergyPowerCapMap::default();
                for argument in function.iter() {
                    degree.add_assign(self.analyze_view(argument)?);
                }
                Ok(degree)
            }
        }
    }

    fn is_energy_index(&self, index: AtomView<'_>) -> bool {
        let AtomView::Fun(function) = index else {
            return false;
        };

        if function.get_symbol() == AIND_SYMBOLS.cind {
            return function.get_nargs() == 1
                && usize::try_from(function.get(0)).is_ok_and(|index| index == 0);
        }

        function.get_symbol() == self.minkowski_symbol && function.get_nargs() == 2
    }

    fn is_internal_edge(&self, edge: EdgeIndex) -> bool {
        self.internal_edges
            .as_ref()
            .is_none_or(|internal_edges| internal_edges.contains(&edge))
    }
}

impl Graph {
    pub fn full_numerator_atom(&self) -> Atom {
        self.numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .expect("Graph numerator should be available")
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor
    }

    pub fn numerator_energy_power_caps(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        let internal_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge));
        EnergyPowerAnalyzer::with_internal_edges(
            self.loop_momentum_basis.loop_edges.iter().copied(),
            internal_edges,
        )
        .analyze_atom(&self.full_numerator_atom())
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::structure::{
        abstract_index::AIND_SYMBOLS,
        representation::{LibraryRep, Minkowski},
    };
    use symbolica::{
        atom::{Atom, AtomCore, FunctionBuilder},
        function, symbol,
    };

    use crate::{
        numerator::energy_degree::EnergyPowerAnalyzer,
        utils::{GS, symbolica_ext::CallSymbol},
    };

    fn mink_index(label: &str) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::var(symbol!(label)).as_view())
            .finish()
    }

    fn bounds(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    fn bounds_with_internal_edges(expression: Atom) -> Vec<(usize, usize)> {
        EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(3)])
            .analyze_atom(&expression)
            .unwrap()
            .into_generation_bounds()
    }

    #[test]
    fn emr_minkowski_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_external_edge_minkowski_index_does_not_count_when_internal_edges_are_known() {
        let expression = function!(GS.emr_mom, 0, mink_index("mu"));
        assert!(bounds_with_internal_edges(expression).is_empty());
    }

    #[test]
    fn emr_temporal_concrete_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.f([0]));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_spatial_concrete_index_does_not_count() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.f([2]));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn spatial_emr_vector_does_not_count() {
        let expression = function!(GS.emr_vec, 3, mink_index("mu"));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn lmb_energy_maps_to_carrier_edge() {
        let expression = function!(GS.loop_mom, 1, mink_index("mu"));
        assert_eq!(bounds(expression), vec![(9, 1)]);
    }

    #[test]
    fn sums_take_edgewise_max_and_products_add() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = q3.clone() * k1.clone() + q3.pow(Atom::num(3));
        assert_eq!(bounds(expression), vec![(3, 3), (9, 1)]);
    }
}
