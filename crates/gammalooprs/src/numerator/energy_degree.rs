use std::{
    collections::{BTreeMap, BTreeSet},
    ops::MulAssign,
};

use linnet::half_edge::involution::EdgeIndex;
use spenso::shadowing::symbolica_utils::LogPrint;
use spenso::structure::{
    abstract_index::AIND_SYMBOLS,
    representation::{LibraryRep, Minkowski},
};
use symbolica::atom::{Atom, AtomView, Symbol};
use thiserror::Error;

use crate::{graph::Graph, utils::GS, uv::UltravioletGraph};

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
        "energy-power analysis expected an integer edge id in Q(edge, index), found `{argument}`"
    )]
    InvalidEmrEdgeArgument { argument: String },
    #[error(
        "energy-power analysis expected an integer loop id in K(loop, index), found `{argument}`"
    )]
    InvalidLoopMomentumArgument { argument: String },
    #[error(
        "non-polynomial energy dependence encountered in numerator subexpression `{expression}` with exponent `{exponent}`"
    )]
    NonPolynomialEnergyPower {
        expression: String,
        exponent: String,
    },
    #[error(
        "energy dependence encountered in function without declared multilinearity: `{expression}`"
    )]
    OpaqueEnergyFunction { expression: String },
}

pub struct EnergyPowerAnalyzer {
    loop_edges: Vec<EdgeIndex>,
    internal_edges: Option<BTreeSet<EdgeIndex>>,
    minkowski_symbol: Symbol,
}

#[derive(Clone, Copy)]
enum EnergyPowerAnalysisMode {
    StrictPolynomial,
    ConservativeUpperBound,
}

#[derive(Clone, Copy)]
enum LorentzIndexKind {
    Temporal,
    Spatial,
    Abstract,
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
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::StrictPolynomial,
        )
    }

    pub fn analyze_atom_upper_bound(
        &self,
        expression: &Atom,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_view(
            expression.as_view(),
            EnergyPowerAnalysisMode::ConservativeUpperBound,
        )
    }

    fn analyze_view(
        &self,
        expression: AtomView<'_>,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        match expression {
            AtomView::Num(_) | AtomView::Var(_) => Ok(EnergyPowerCapMap::default()),
            AtomView::Add(add) => {
                let mut max_degree = EnergyPowerCapMap::default();
                for term in add.iter() {
                    max_degree.max_assign(self.analyze_view(term, mode)?);
                }
                Ok(max_degree)
            }
            AtomView::Mul(mul) => {
                let mut product_degree = EnergyPowerCapMap::default();
                for factor in mul.iter() {
                    product_degree.add_assign(self.analyze_view(factor, mode)?);
                }
                Ok(product_degree)
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                let mut base_degree = self.analyze_view(base, mode)?;
                if !self.analyze_view(exponent, mode)?.is_empty() {
                    return Err(EnergyPowerAnalysisError::NonPolynomialEnergyPower {
                        expression: base.to_owned().log_print(None),
                        exponent: exponent.to_owned().log_print(None),
                    });
                }
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
                    if matches!(mode, EnergyPowerAnalysisMode::ConservativeUpperBound) {
                        return Ok(EnergyPowerCapMap::default());
                    }
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
                    let edge = EdgeIndex(usize::try_from(function.get(0)).map_err(|_| {
                        EnergyPowerAnalysisError::InvalidEmrEdgeArgument {
                            argument: function.get(0).to_owned().log_print(None),
                        }
                    })?);
                    return Ok(if self.is_internal_edge(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                if function.get_nargs() == 2 && function.get_symbol() == GS.loop_mom {
                    let loop_index = usize::try_from(function.get(0)).map_err(|_| {
                        EnergyPowerAnalysisError::InvalidLoopMomentumArgument {
                            argument: function.get(0).to_owned().log_print(None),
                        }
                    })?;
                    let edge =
                        self.loop_edges.get(loop_index).copied().ok_or(
                            EnergyPowerAnalysisError::MissingLoopCarrierEdge { loop_index },
                        )?;
                    return Ok(if self.is_internal_edge(edge) {
                        self.component_degree(edge, function.get(1))
                    } else {
                        EnergyPowerCapMap::default()
                    });
                }

                let argument_degrees = function
                    .iter()
                    .map(|argument| self.analyze_view(argument, mode))
                    .collect::<Result<Vec<_>, _>>()?;
                if argument_degrees.iter().any(|degree| !degree.is_empty())
                    && !function.get_symbol().is_linear()
                {
                    return Err(EnergyPowerAnalysisError::OpaqueEnergyFunction {
                        expression: Atom::from(function.to_owned()).log_print(None),
                    });
                }
                if function.get_symbol().is_linear() {
                    let mut degree = EnergyPowerCapMap::default();
                    for argument_degree in argument_degrees {
                        degree.add_assign(argument_degree);
                    }
                    Ok(degree)
                } else {
                    Ok(EnergyPowerCapMap::default())
                }
            }
        }
    }

    fn component_degree(&self, edge: EdgeIndex, index: AtomView<'_>) -> EnergyPowerCapMap {
        match self.lorentz_index_kind(index) {
            LorentzIndexKind::Temporal | LorentzIndexKind::Abstract => {
                EnergyPowerCapMap::unit(edge)
            }
            LorentzIndexKind::Spatial => EnergyPowerCapMap::default(),
        }
    }

    fn lorentz_index_kind(&self, index: AtomView<'_>) -> LorentzIndexKind {
        if let Ok(index) = usize::try_from(index) {
            return if index == 0 {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        let AtomView::Fun(function) = index else {
            return LorentzIndexKind::Spatial;
        };

        if function.get_symbol() == AIND_SYMBOLS.cind {
            return if function.get_nargs() == 1
                && usize::try_from(function.get(0)).is_ok_and(|index| index == 0)
            {
                LorentzIndexKind::Temporal
            } else {
                LorentzIndexKind::Spatial
            };
        }

        if function.get_symbol() == self.minkowski_symbol && function.get_nargs() == 2 {
            // `mink(4, n)` is an abstract Lorentz-index label in the numerator
            // algebra, not a concrete component selection. Concrete components
            // use `cind(n)` and are handled above. Any abstract Minkowski index
            // can contract onto the temporal component, so it contributes to the
            // conservative energy-power bound.
            return LorentzIndexKind::Abstract;
        }

        LorentzIndexKind::Spatial
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
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::StrictPolynomial)
    }

    pub(crate) fn numerator_energy_power_upper_bounds(
        &self,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        self.analyze_numerator_energy_powers(EnergyPowerAnalysisMode::ConservativeUpperBound)
    }

    fn analyze_numerator_energy_powers(
        &self,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        // Keep edge-energy variables attached to their source denominator.
        // Rewriting them to a loop-momentum basis loses cancellation
        // information in non-coordinate UV directions, e.g. when q_i-q_j is
        // kept fixed while q_i and q_j are both large.
        let numerator = self.full_numerator_atom();
        self.analyze_numerator_atom_energy_powers(&numerator, mode)
    }

    fn analyze_numerator_atom_energy_powers(
        &self,
        numerator: &Atom,
        mode: EnergyPowerAnalysisMode,
    ) -> Result<EnergyPowerCapMap, EnergyPowerAnalysisError> {
        let internal_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy).then_some(edge)
            });

        EnergyPowerAnalyzer::with_internal_edges(
            self.loop_momentum_basis.loop_edges.iter().copied(),
            internal_edges,
        )
        .analyze_view(numerator.as_view(), mode)
    }

    pub fn automatic_numerator_energy_degree_bounds(
        &self,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        Ok(self
            .numerator_energy_power_upper_bounds()?
            .into_generation_bounds())
    }

    pub(crate) fn automatic_numerator_energy_degree_bounds_in_atom_excluding_with_min_degree(
        &self,
        numerator: &Atom,
        excluded_edges: impl IntoIterator<Item = EdgeIndex>,
        min_degree: usize,
    ) -> Result<Vec<(usize, usize)>, EnergyPowerAnalysisError> {
        let excluded_edges = excluded_edges.into_iter().collect::<BTreeSet<_>>();
        let active_edges = self
            .underlying
            .iter_edges()
            .filter_map(|(pair, edge, edge_data)| {
                (pair.is_paired() && !edge_data.data.is_dummy && !excluded_edges.contains(&edge))
                    .then_some(edge)
            });
        Ok(EnergyPowerAnalyzer::with_internal_edges(
            self.loop_momentum_basis.loop_edges.iter().copied(),
            active_edges,
        )
        .analyze_atom(numerator)?
        .iter()
        .filter_map(|(edge, degree)| (degree >= min_degree).then_some((edge.into(), degree)))
        .collect())
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
        numerator::energy_degree::{EnergyPowerAnalysisError, EnergyPowerAnalyzer},
        utils::GS,
    };

    fn mink_index(label: &str) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::var(symbol!(label)).as_view())
            .finish()
    }

    fn mink_component(index: i64) -> Atom {
        FunctionBuilder::new(LibraryRep::from(Minkowski {}).symbol())
            .add_arg(Atom::num(4).as_view())
            .add_arg(Atom::num(index).as_view())
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
    fn emr_numeric_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(1));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_numeric_zero_minkowski_label_counts_as_abstract_energy_power() {
        let expression = function!(GS.emr_mom, 3, mink_component(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_external_edge_minkowski_index_does_not_count_when_internal_edges_are_known() {
        let expression = function!(GS.emr_mom, 0, mink_index("mu"));
        assert!(bounds_with_internal_edges(expression).is_empty());
    }

    #[test]
    fn emr_temporal_concrete_index_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(0));
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_spatial_concrete_index_does_not_count() {
        let expression = function!(GS.emr_mom, 3, AIND_SYMBOLS.cind.call(2));
        assert!(bounds(expression).is_empty());
    }

    #[test]
    fn emr_direct_temporal_component_counts_as_energy_power() {
        let expression = function!(GS.emr_mom, 3, 0);
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn emr_direct_spatial_component_does_not_count() {
        let expression = function!(GS.emr_mom, 3, 1);
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

    #[test]
    fn powers_of_affine_atoms_produce_edgewise_bounds() {
        let k0 = function!(GS.loop_mom, 0, mink_index("mu"));
        let k1 = function!(GS.loop_mom, 1, mink_index("nu"));
        let expression = (k0 + k1).pow(Atom::num(2));
        assert_eq!(bounds(expression), vec![(7, 2), (9, 2)]);
    }

    #[test]
    fn factorized_products_compose_energy_bounds_without_expansion() {
        let q3 = function!(GS.emr_mom, 3, mink_index("mu"));
        let expression = (q3.clone() + Atom::num(1)) * (q3 + Atom::num(2));
        assert_eq!(bounds(expression), vec![(3, 2)]);
    }

    #[test]
    fn strict_analysis_rejects_negative_energy_powers() {
        let expression = function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1));
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_ignores_energy_powers_outside_the_active_edge_set() {
        let analyzer =
            EnergyPowerAnalyzer::with_internal_edges([EdgeIndex(7), EdgeIndex(9)], [EdgeIndex(5)]);
        for expression in [
            function!(GS.emr_mom, 3, mink_index("mu")).pow(Atom::num(-1)),
            function!(GS.loop_mom, 0, mink_index("mu")).pow(Atom::num(-1)),
        ] {
            assert!(analyzer.analyze_atom(&expression).unwrap().is_empty());
        }
    }

    #[test]
    fn energy_dependent_exponents_are_rejected_even_for_constant_bases() {
        let exponent = function!(GS.emr_mom, 3, 0);
        let expression = Atom::num(2).pow(exponent);
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::NonPolynomialEnergyPower { .. }
        ));
    }

    #[test]
    fn strict_analysis_rejects_opaque_energy_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function"))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn zero_arguments_do_not_hide_energy_in_opaque_functions() {
        let expression = FunctionBuilder::new(symbol!("opaque_energy_function_with_zero"))
            .add_arg(Atom::zero().as_view())
            .add_arg(function!(GS.emr_mom, 3, 0).as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom_upper_bound(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::OpaqueEnergyFunction { .. }
        ));
    }

    #[test]
    fn declared_multilinear_functions_preserve_energy_bounds() {
        let expression = FunctionBuilder::new(symbol!("multilinear_energy_function"; Linear))
            .add_arg(function!(GS.emr_mom, 3, mink_index("mu")).as_view())
            .finish();
        assert_eq!(bounds(expression), vec![(3, 1)]);
    }

    #[test]
    fn malformed_emr_edge_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.emr_mom)
            .add_arg(Atom::var(symbol!("edge")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidEmrEdgeArgument { .. }
        ));
    }

    #[test]
    fn malformed_loop_argument_is_reported() {
        let expression = FunctionBuilder::new(GS.loop_mom)
            .add_arg(Atom::var(symbol!("loop_id")).as_view())
            .add_arg(mink_index("mu").as_view())
            .finish();
        let error = EnergyPowerAnalyzer::new([EdgeIndex(7), EdgeIndex(9)])
            .analyze_atom(&expression)
            .unwrap_err();
        assert!(matches!(
            error,
            EnergyPowerAnalysisError::InvalidLoopMomentumArgument { .. }
        ));
    }
}
