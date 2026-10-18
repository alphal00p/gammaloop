use std::collections::BTreeMap;

use color_eyre::eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec, Orientation};
use symbolica::atom::{Atom, AtomCore, AtomView, Indeterminate};
use three_dimensional_reps::{
    ThermalDistributionCall, ThermalDistributionFactor, thermal::OCCUPATION_FORM,
};

use crate::{
    cff::{expression::OrientationID, orientations::GraphOrientation, surface::HybridSurfaceRef},
    graph::Graph,
    momentum::SignOrZero,
    numerator::symbolica_ext::NumeratorAtomExt,
    processes::{build_derivative_structure_atom, params_for_derivative_order},
    utils::GS,
    uv::uv_graph::UVE,
};

use super::FermiSurfaceProduct;

/// A coefficient and its independent, oriented zero-temperature delta support.
/// Numerator functions remain opaque; their existing definitions accompany the
/// coefficient into evaluator construction.
pub(crate) struct FermiSurfaceSector {
    pub(crate) product: FermiSurfaceProduct,
    pub(crate) orientations: Vec<i32>,
    pub(crate) coefficient: Atom,
}

impl FermiSurfaceSector {
    fn has_fermi_surface(graph: &Graph, edge_id: EdgeIndex) -> Result<bool> {
        let (_, _, edge) = graph
            .iter_edges()
            .find(|(_, edge, _)| *edge == edge_id)
            .ok_or_else(|| eyre!("Unknown thermal distribution edge {edge_id}"))?;
        Ok(edge.data.is_fermion() && edge.data.chemical_potential_atom().is_some())
    }

    /// The ordinary dual evaluator differentiates cold steps only away from
    /// their boundary. Retain the missing distributional part of a double LU
    /// pole as a simple-cut coefficient before building that evaluator.
    pub(crate) fn lu_boundary(
        graph: &Graph,
        expression: &Atom,
        order: usize,
        surface: HybridSurfaceRef<'_>,
        catalogue: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
    ) -> Result<Atom> {
        if order <= 1 {
            return Ok(Atom::Zero);
        }
        let mut branches = Vec::new();
        for (_, branch, selector) in Self::oriented_branches(expression, catalogue)? {
            // Thermal wrappers contain only the weight polynomial. Normalize
            // contour offsets and cancel its polynomial terms here, before a
            // derivative can attach them to the factorized smooth coefficient.
            // Neither the surrounding coefficient nor graph numerators expand.
            let branch = branch
                .replace_multiple(&*OCCUPATION_FORM)
                .replace_map(|atom, _, out| {
                    if let AtomView::Fun(function) = atom
                        && function.get_symbol() == GS.thermal_weight_wrapper
                        && function.get_nargs() == 1
                    {
                        let weight = function.iter().next().unwrap().expand();
                        **out = symbolica::function!(GS.thermal_weight_wrapper, weight);
                    }
                });
            let mut active = Vec::new();
            for (key, call) in Self::calls(&branch)? {
                if !call.temperature_flag.is_zero() || !Self::has_fermi_surface(graph, call.edge)? {
                    continue;
                }
                ensure!(
                    call.derivative_order == 0,
                    "Raised cold LU coefficients with existing Fermi distributions require higher joint intersections"
                );
                active.push((key, call));
            }
            if !active.is_empty() {
                ensure!(
                    order == 2,
                    "Cold occupation derivatives of LU poles above order two require higher joint intersections"
                );
                branches.push((branch, selector, active));
            }
        }
        if branches.is_empty() {
            return Ok(Atom::Zero);
        }
        let (energies, _) = FermiSurfaceProduct::cut_equation(surface)?;
        let eta_derivative = energies.into_iter().fold(Atom::Zero, |sum, (edge, sign)| {
            sum + Atom::num(sign)
                * graph
                    .loop_momentum_basis
                    .radial_energy_derivative_atom(edge)
        });
        let parameters = params_for_derivative_order(2);
        // Extract the f' coefficient from the shared Laurent owner itself.
        // Multiplying by eta' converts its radial residue weight into a
        // coefficient of the simple cut delta used by the joint localizer.
        let laurent_multiplier = (build_derivative_structure_atom(2, -1)
            .replace(parameters[1].to_pattern())
            .with(Atom::one())
            .replace(parameters[0].to_pattern())
            .with(Atom::Zero)
            * &parameters[2])
            .replace(parameters[2].to_pattern())
            .with(eta_derivative);
        let mut boundary = Atom::Zero;
        for (branch, selector, active) in branches {
            let smooth = branch.unwrap_function(GS.thermal_weight_wrapper);
            for (key, call) in active {
                let variable = Indeterminate::try_from(key)
                    .map_err(|error| eyre!("Invalid thermal differentiation variable: {error}"))?;
                let coefficient = smooth.derivative(&variable);
                if coefficient.is_zero() {
                    continue;
                }
                boundary += coefficient
                    * ThermalDistributionCall::atom(
                        call.edge.0 as i64,
                        1,
                        0,
                        call.thermal_sign,
                        call.orientation_sign,
                    )
                    * graph
                        .loop_momentum_basis
                        .radial_energy_derivative_atom(call.edge)
                    * &laurent_multiplier
                    * &selector;
            }
        }
        Ok(boundary)
    }

    fn calls(expression: &Atom) -> Result<Vec<(Atom, ThermalDistributionCall)>> {
        let mut calls = Vec::new();
        let mut error = None;
        expression.visitor(&mut |atom| {
            if error.is_some() {
                return false;
            }
            if let AtomView::Fun(function) = atom
                && function.get_symbol() == GS.thermal_weight_wrapper
                && function.get_nargs() != 1
            {
                error = Some(eyre!(
                    "Thermal weight wrapper must have one argument, got {atom}"
                ));
                return false;
            }
            if let AtomView::Fun(function) = atom
                && function.get_symbol() == GS.thermal_distribution
            {
                let call = ThermalDistributionCall::try_from(atom)
                    .map_err(Into::into)
                    .and_then(|call| {
                        ensure!(
                            call.temperature_flag.is_zero() || call.temperature_flag.is_one(),
                            "Thermal distribution temperature flag must be zero or one, got {}",
                            call.temperature_flag
                        );
                        Ok(call)
                    });
                match call {
                    Ok(call) => {
                        if !calls
                            .iter()
                            .any(|(key, _): &(Atom, _)| key.as_view() == atom)
                        {
                            calls.push((atom.to_owned(), call));
                        }
                    }
                    Err(parse_error) => error = Some(parse_error),
                }
            }
            true
        });
        if let Some(error) = error {
            return Err(error);
        }
        Ok(calls)
    }

    /// Detect distribution derivatives before resolving physical orientations.
    pub(crate) fn is_present(expression: &Atom) -> Result<bool> {
        Ok(Self::calls(expression)?
            .iter()
            .any(|(_, call)| call.temperature_flag.is_zero() && call.derivative_order > 0))
    }

    fn oriented_branches(
        expression: &Atom,
        catalogue: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
    ) -> Result<Vec<(Option<OrientationID>, Atom, Atom)>> {
        let branches = if let Some((orientations, production_ids)) = catalogue {
            ensure!(
                orientations.len() == production_ids.len(),
                "Fermi extraction requires matching physical and production orientation catalogues"
            );
            orientations
                .iter()
                .zip(production_ids)
                .enumerate()
                .map(|(index, (orientation, production_id))| {
                    // Remove inactive residue maps before resolving physical
                    // directions, just as the ordinary evaluator does.
                    let selected = production_id.select(expression);
                    (
                        Some(OrientationID(index)),
                        GS.collect_orientation_if(orientation.select(selected)),
                        production_id.atom(),
                    )
                })
                .collect()
        } else {
            vec![(None, expression.clone(), Atom::one())]
        };
        Ok(branches)
    }

    /// Resolve the same physical and residue-map channels as the smooth
    /// evaluator before extracting distribution support. Each coefficient
    /// retains its production selector, including repeated physical directions.
    // The tuple keeps the regular expression and orientation-tagged supports
    // together without adding another persistent sector representation.
    #[allow(clippy::type_complexity)]
    pub(crate) fn extract_orientations(
        graph: &Graph,
        expression: &Atom,
        catalogue: Option<(&[EdgeVec<Orientation>], &[OrientationID])>,
    ) -> Result<(Atom, Vec<(Option<OrientationID>, Self)>)> {
        if !Self::is_present(expression)? {
            return Ok((expression.clone(), Vec::new()));
        }
        let branches = Self::oriented_branches(expression, catalogue)?;
        let mut bulk = Atom::Zero;
        let mut localized = Vec::new();
        for (orientation_id, branch, selector) in branches {
            let (regular, sectors) = Self::extract(graph, &branch)?;
            bulk += regular * &selector;
            for mut sector in sectors {
                sector.coefficient *= &selector;
                localized.push((orientation_id, sector));
            }
        }
        Ok((bulk, localized))
    }

    /// Separate regular terms and Fermi supports after residue and physical
    /// orientation selection. The fifth N argument must be a definite sign;
    /// distinct signs in an explicit sum remain distinct supports.
    pub(crate) fn extract(graph: &Graph, expression: &Atom) -> Result<(Atom, Vec<Self>)> {
        // Only fermions with a chemical potential can have a zero-temperature
        // Fermi surface. A missing chemical potential is zero, and warm-up
        // restricts bosons to |mu| < m or m = mu = 0. Either occupation is then
        // constant on the shell E > 0, so its positive-order derivatives vanish.
        let calls = Self::calls(expression)?;
        let mut expression = expression.unwrap_function(GS.thermal_weight_wrapper);
        let mut active = Vec::new();
        for (key, call) in calls {
            if !call.temperature_flag.is_zero() || call.derivative_order == 0 {
                continue;
            }
            if !Self::has_fermi_surface(graph, call.edge)? {
                expression = expression.replace(key.to_pattern()).with(Atom::Zero);
                continue;
            }
            let sign = i64::try_from(call.thermal_sign.as_view())
                .map_err(|_| eyre!("Fermi-surface thermal sign must be resolved in {key}"))?;
            let orientation = i64::try_from(call.orientation_sign.as_view()).map_err(|_| {
                eyre!("Fermi-surface chemical orientation must be resolved in {key}")
            })?;
            ensure!(
                matches!(sign, -1..=1) && matches!(orientation, -1 | 1),
                "Fermi-surface thermal sign must be -1, 0 or +1 and chemical orientation must be +1 or -1 in {key}"
            );
            active.push((
                key,
                ThermalDistributionFactor {
                    edge_id: call.edge,
                    sign: sign as i32,
                    derivative_order: call.derivative_order,
                },
                orientation as i32,
            ));
        }
        if active.is_empty() || expression.is_zero() {
            return Ok((expression, Vec::new()));
        }
        let contains_active = |expression: AtomView<'_>| {
            let mut found = false;
            expression.visitor(&mut |part| {
                found |= active.iter().any(|(key, _, _)| key.as_view() == part);
                !found
            });
            found
        };
        let mut invalid = None;
        expression.visitor(&mut |part| {
            if invalid.is_some() || active.iter().any(|(key, _, _)| key.as_view() == part) {
                return false;
            }
            if matches!(part, AtomView::Fun(_) | AtomView::Pow(_)) && contains_active(part) {
                invalid = Some(part.to_owned());
                return false;
            }
            true
        });
        ensure!(
            invalid.is_none(),
            "Fermi distributions must occur polynomially outside functions, inverses and powers: {}",
            invalid.unwrap_or(Atom::Zero)
        );

        let keys = active
            .iter()
            .map(|(key, _, _)| key.clone())
            .collect::<Vec<_>>();
        let mut bulk = Atom::Zero;
        let mut grouped = BTreeMap::<Vec<(ThermalDistributionFactor, i32)>, Atom>::new();
        // Collect only tagged distributions: graph numerator bodies and their
        // factorized symbolic references never enter polynomial expansion.
        for (key, coefficient) in expression.coefficient_list::<u32>(&keys) {
            ensure!(
                !contains_active(coefficient.as_view()),
                "Fermi distribution remains hidden in coefficient {coefficient}"
            );
            if key.is_one() {
                bulk += coefficient;
                continue;
            }
            // Expansion can multiply equal factors; a squared delta is undefined.
            let power = |factor: AtomView<'_>| matches!(factor, AtomView::Pow(_));
            ensure!(
                !power(key.as_view())
                    && !matches!(key.as_view(), AtomView::Mul(factors) if factors.iter().any(power)),
                "Repeated Fermi distribution {key} has no defined product"
            );
            let mut support = active
                .iter()
                .filter(|(candidate, _, _)| {
                    let mut found = false;
                    key.visitor(&mut |part| {
                        found |= candidate.as_view() == part;
                        !found
                    });
                    found
                })
                .map(|(_, factor, orientation)| (*factor, *orientation))
                .collect::<Vec<_>>();
            support.sort();
            *grouped.entry(support).or_insert(Atom::Zero) += coefficient;
        }
        let mut sectors = Vec::with_capacity(grouped.len());
        for (support, mut coefficient) in grouped {
            if coefficient.is_zero() {
                continue;
            }
            let (factors, orientations): (Vec<_>, Vec<_>) = support.into_iter().unzip();
            let product = FermiSurfaceProduct::new(graph, &factors)?;
            ensure!(
                !coefficient.contains_symbol(GS.heaviside),
                "An explicit step in a Fermi coefficient has no certified independent routing"
            );
            for (key, step) in Self::calls(&coefficient)? {
                if !step.temperature_flag.is_zero()
                    || step.derivative_order != 0
                    || !Self::has_fermi_surface(graph, step.edge)?
                {
                    continue;
                }
                // Active shells have nu > m >= 0. An exactly opposite
                // oriented chemical potential has no occupied positive-energy
                // states anywhere on this support. Retain N's contour offset.
                let step_mu = &step.orientation_sign
                    * graph[step.edge]
                        .chemical_potential_atom()
                        .expect("validated cold step has a chemical potential");
                if factors
                    .iter()
                    .zip(&orientations)
                    .any(|(factor, orientation)| {
                        let shell_mu = Atom::num(*orientation)
                            * graph[factor.edge_id]
                                .chemical_potential_atom()
                                .expect("validated shell has a chemical potential");
                        (&step_mu + shell_mu).expand().is_zero()
                    })
                {
                    coefficient = coefficient
                        .replace(key.to_pattern())
                        .with(ThermalDistributionFactor::vacuum_atom(0, step.thermal_sign));
                    if coefficient.is_zero() {
                        break;
                    }
                    continue;
                }
                let routing = product
                    .lmb()
                    .edge_signatures
                    .get(step.edge)
                    .ok_or_else(|| {
                        eyre!("Missing routing for thermal step on edge {}", step.edge)
                    })?;
                ensure!(
                    routing
                        .internal
                        .iter()
                        .take(factors.len())
                        .all(|sign| *sign == SignOrZero::Zero),
                    "Thermal step on edge {} depends on an active Fermi momentum; its boundary intersection is unsupported in production localization",
                    step.edge
                );
            }
            if coefficient.is_zero() {
                continue;
            }
            sectors.push(Self {
                product,
                orientations,
                coefficient,
            });
        }
        Ok((bulk, sectors))
    }
}

#[cfg(test)]
mod tests {
    use symbolica::{function, symbol};

    use super::*;
    use crate::{
        dot, graph::parse::IntoGraph, initialisation::test_initialise,
        integrands::process::fermi_surface::routing,
    };

    #[test]
    fn fermi_sector_extraction_preserves_factorized_coefficients_and_regular_terms() -> Result<()> {
        let graph = routing::test_graph()?;
        let x = Atom::var(symbol!("fermi_sector_x"));
        let y = Atom::var(symbol!("fermi_sector_y"));
        let numerator = function!(symbol!("fermi_sector_numerator"), &x, &y);
        let coefficient = numerator * (&x + Atom::num(1)) * (&y + Atom::num(2));
        // Keep a nontrivial occupation on the independent momentum. The
        // opposite chemical orientation is identically empty on this shell.
        let independent_step = ThermalDistributionCall::atom(3, 0, 0, -1, 1);
        let bulk = Atom::num(7) * ThermalDistributionCall::atom(1, 3, 1, 1, 1);
        for derivative_order in [2, 3] {
            for sign in [-1, 0, 1] {
                let delta = ThermalDistributionCall::atom(1, derivative_order as i64, 0, sign, 1);
                let expression = &bulk
                    + function!(GS.thermal_weight_wrapper, &delta)
                        * &coefficient
                        * &independent_step;
                let (actual_bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
                assert_eq!(actual_bulk, bulk);
                assert_eq!(sectors.len(), 1);
                assert_eq!(sectors[0].coefficient, &coefficient * &independent_step);
                assert_eq!(sectors[0].orientations, [1]);
                assert_eq!(
                    sectors[0].product.factors(),
                    &[ThermalDistributionFactor {
                        edge_id: EdgeIndex(1),
                        sign,
                        derivative_order,
                    }]
                );
                assert_eq!(sectors[0].product.derivative_orders, [derivative_order - 1]);
                assert_eq!(sectors[0].product.shape.len(), derivative_order);
                assert_eq!(
                    actual_bulk + &sectors[0].coefficient * delta,
                    expression.unwrap_function(GS.thermal_weight_wrapper)
                );
            }
        }
        Ok(())
    }

    #[test]
    fn fermi_sectors_keep_opposite_chemical_supports_and_independent_products() -> Result<()> {
        let graph = routing::test_graph()?;
        let plus = ThermalDistributionCall::atom(1, 1, 0, 1, 1);
        let minus = ThermalDistributionCall::atom(1, 1, 0, 1, -1);
        let second = ThermalDistributionCall::atom(3, 2, 0, -1, -1);
        let expression =
            Atom::num(2) * &plus + Atom::num(3) * &minus + Atom::num(5) * &plus * &second;
        let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
        assert!(bulk.is_zero());
        assert_eq!(sectors.len(), 3);
        let reconstructed = sectors.iter().fold(bulk, |sum, sector| {
            sum + sector
                .product
                .factors()
                .iter()
                .zip(&sector.orientations)
                .fold(
                    sector.coefficient.clone(),
                    |coefficient, (factor, orientation)| {
                        coefficient
                            * ThermalDistributionCall::atom(
                                factor.edge_id.0 as i64,
                                factor.derivative_order as i64,
                                0,
                                factor.sign,
                                *orientation,
                            )
                    },
                )
        });
        assert_eq!(reconstructed, expression);
        assert!(sectors.iter().any(|sector| sector.orientations == [1, -1]));
        Ok(())
    }

    #[test]
    fn fermi_sector_detector_accepts_unresolved_orientation_and_ignores_regular_weights()
    -> Result<()> {
        let graph = routing::test_graph()?;
        let unresolved = ThermalDistributionCall::atom(1, 1, 0, 1, GS.sign(EdgeIndex(1)));
        assert!(FermiSurfaceSector::is_present(&unresolved)?);
        assert!(FermiSurfaceSector::extract(&graph, &unresolved).is_err());
        for expression in [
            ThermalDistributionCall::atom(1, 0, 0, 1, 1),
            ThermalDistributionCall::atom(1, 2, 1, 1, 1),
        ] {
            assert!(!FermiSurfaceSector::is_present(&expression)?);
            let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
            assert_eq!(bulk, expression);
            assert!(sectors.is_empty());
        }
        for derivative_order in [1, 3] {
            let (bulk, sectors) = FermiSurfaceSector::extract(
                &graph,
                &ThermalDistributionCall::atom(4, derivative_order, 0, 1, 1),
            )?;
            assert!(bulk.is_zero());
            assert!(sectors.is_empty());
        }
        Ok(())
    }

    #[test]
    fn fermi_sector_extraction_drops_boson_derivatives() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph boson_cycle {
            node [num=1]
            edge [num=1 particle="W+"]
            A -> B [id=0 lmb_id=0]
            B -> A [id=1]
        })?;
        assert!(graph[EdgeIndex(0)].chemical_potential_atom().is_some());
        let step = ThermalDistributionCall::atom(1, 0, 0, 1, 1);
        let expression = ThermalDistributionCall::atom(0, 2, 0, 1, 1) * &step + &step;
        let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
        assert_eq!(bulk, step);
        assert!(sectors.is_empty());
        Ok(())
    }

    #[test]
    fn fermi_sectors_treat_literal_zero_chemical_potential_as_absent() -> Result<()> {
        let mut graph = routing::test_graph()?;
        let edge = EdgeIndex(1);
        let mut particle = graph[edge].particle().unwrap();
        for chemical_potential in [
            None,
            Some(crate::model::ParameterName(crate::model::UFOSymbol::zero())),
        ] {
            std::sync::Arc::make_mut(&mut particle.0).chemical_potential = chemical_potential;
            graph.underlying[edge].particle = particle.clone().into();
            for derivative_order in [1, 3] {
                let expression =
                    ThermalDistributionCall::atom(1, derivative_order, 0, 1, 1) + Atom::num(7);
                let (bulk, sectors) = FermiSurfaceSector::extract(&graph, &expression)?;
                assert_eq!(bulk, Atom::num(7));
                assert!(sectors.is_empty());
                let factor = ThermalDistributionFactor {
                    edge_id: edge,
                    sign: 1,
                    derivative_order: derivative_order as usize,
                };
                let error = FermiSurfaceProduct::new(&graph, &[factor]).err().unwrap();
                assert!(error.to_string().contains("has no chemical potential"));
            }
        }
        Ok(())
    }

    #[test]
    fn fermi_sector_extraction_rejects_undefined_products_and_uncertified_steps() -> Result<()> {
        let graph = routing::test_graph()?;
        let delta = ThermalDistributionCall::atom(1, 1, 0, 1, 1);
        // Only polynomial expansion exposes this square.
        let expanded = (&delta + Atom::var(symbol!("fermi_sector_x"))) * (&delta + Atom::num(2));
        let Err(error) = FermiSurfaceSector::extract(&graph, &expanded) else {
            panic!("accepted {expanded}");
        };
        assert!(
            error.to_string().contains("Repeated Fermi distribution"),
            "{error}"
        );
        for malformed in [
            delta.clone().pow(2),
            delta.clone().pow(-1),
            function!(symbol!("fermi_sector_hidden"), &delta),
            &delta * ThermalDistributionCall::atom(0, 1, 0, 1, 1),
            &delta * ThermalDistributionCall::atom(1, 2, 0, 1, 1),
            &delta * ThermalDistributionCall::atom(1, 0, 0, 1, 1),
            &delta * ThermalDistributionCall::atom(0, 0, 0, 1, 1),
            &delta * GS.heaviside(Atom::var(symbol!("fermi_sector_hidden_step"))),
            ThermalDistributionCall::atom(99, 1, 0, 1, 1),
            ThermalDistributionCall::atom(1, 1, 0, 2, 1),
            ThermalDistributionCall::atom(1, 1, 0, 1, 0),
            ThermalDistributionCall::atom(1, -1, 0, 1, 1),
            ThermalDistributionCall::atom(1, 1, 2, 1, 1),
            function!(GS.thermal_distribution, 1, 1, 0, 1),
        ] {
            assert!(
                FermiSurfaceSector::extract(&graph, &malformed).is_err(),
                "accepted {malformed}"
            );
        }
        Ok(())
    }

    #[test]
    fn opposite_chemical_potential_steps_keep_the_contour_offset() -> Result<()> {
        let graph = routing::test_graph()?;
        for orientation in [-1, 1] {
            for thermal_sign in [-1, 0, 1] {
                let delta = ThermalDistributionCall::atom(1, 1, 0, 1, orientation);
                let opposite = ThermalDistributionCall::atom(0, 0, 0, thermal_sign, -orientation);
                let (bulk, sectors) =
                    FermiSurfaceSector::extract(&graph, &(Atom::num(7) * &delta * opposite))?;
                assert!(bulk.is_zero());
                if thermal_sign == -1 {
                    assert!(sectors.is_empty());
                } else {
                    assert_eq!(sectors.len(), 1);
                    assert_eq!(sectors[0].orientations, vec![orientation]);
                    assert_eq!(
                        sectors[0].coefficient,
                        Atom::num(7 * (1 + thermal_sign)) / Atom::num(2),
                    );
                }
                let same = ThermalDistributionCall::atom(0, 0, 0, thermal_sign, orientation);
                assert!(FermiSurfaceSector::extract(&graph, &(delta * same)).is_err());
            }
        }
        Ok(())
    }

    #[test]
    fn cold_lu_boundary_preserves_numerator_and_original_affine_dilation() -> Result<()> {
        use crate::{
            cff::{VertexSet, esurface::Esurface, hsurface::Hsurface},
            integrands::process::param_builder::ParamBuilderGraph,
        };
        test_initialise()?;
        let graph: Graph = dot!(digraph cold_lu_bubble {
            node [num=1]
            edge [num=1]
            ext [style=invis is_cut=0]
            ext -> A [id=0 particle="H"]
            A -> B [id=1 lmb_id=0 particle="b"]
            B -> A [id=2 particle="b"]
            B -> ext [id=3 particle="H"]
        })?;
        let surface = Esurface {
            energies: vec![EdgeIndex(1), EdgeIndex(2)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        // An actual edge-local numerator sum must remain a factor; an opaque
        // function alone would not detect accidental whole-coefficient expansion.
        let numerator =
            GS.emr_mom(EdgeIndex(1), GS.cind(1)).pow(2) + graph[EdgeIndex(1)].mass_atom().pow(2);
        let prefactor =
            &numerator * Atom::var(GS.rescale_star).pow(3) * Atom::var(GS.hfunction_lu_cut);
        let d1 = graph
            .loop_momentum_basis
            .radial_energy_derivative_atom(EdgeIndex(1));
        let d2 = graph
            .loop_momentum_basis
            .radial_energy_derivative_atom(EdgeIndex(2));
        for orientation in [-1, 1] {
            let occupation = ThermalDistributionCall::atom(1, 0, 0, -1, orientation);
            let coefficient = &prefactor * function!(GS.thermal_weight_wrapper, occupation);
            let lifted = FermiSurfaceSector::lu_boundary(
                &graph,
                &coefficient,
                2,
                HybridSurfaceRef::Esurface(&surface),
                None,
            )?;
            let expected =
                &prefactor * ThermalDistributionCall::atom(1, 1, 0, -1, orientation) * &d1
                    / (&d1 + &d2);
            assert_eq!(lifted, expected);
            assert!(
                FermiSurfaceSector::lu_boundary(
                    &graph,
                    &coefficient,
                    3,
                    HybridSurfaceRef::Esurface(&surface),
                    None,
                )
                .is_err()
            );
        }
        let crossed = Hsurface {
            positive_energies: vec![EdgeIndex(1)],
            negative_energies: vec![EdgeIndex(2)],
            external_shift: vec![(EdgeIndex(0), -1)],
            vertex_set: VertexSet::dummy(),
        };
        for orientation in [-1, 1] {
            let [first_plus, first_minus, second_plus, second_minus] =
                [(1, 1), (1, -1), (2, 1), (2, -1)].map(|(edge, sign)| {
                    ThermalDistributionCall::atom(edge, 0, 0, sign, orientation)
                });
            let weight = first_plus * &second_minus - &first_minus * second_plus;
            let canonical = (&first_minus + Atom::one()) * &second_minus
                - &first_minus * (&second_minus + Atom::one());
            let expected = &prefactor
                * (ThermalDistributionCall::atom(2, 1, 0, -1, orientation) * &d2
                    - ThermalDistributionCall::atom(1, 1, 0, -1, orientation) * &d1)
                / (&d1 - &d2);
            let mut reference = None;
            for weight in [weight, canonical, second_minus - first_minus] {
                let source = &prefactor * function!(GS.thermal_weight_wrapper, weight);
                let lifted = FermiSurfaceSector::lu_boundary(
                    &graph,
                    &source,
                    2,
                    HybridSurfaceRef::Hsurface(&crossed),
                    None,
                )?;
                assert!(
                    (lifted.clone() - &expected)
                        .replace(numerator.to_pattern())
                        .with(Atom::one())
                        .expand()
                        .is_zero()
                );
                assert!(
                    FermiSurfaceSector::calls(&lifted)?
                        .iter()
                        .all(|(_, call)| call.derivative_order == 1)
                );
                let (_, sectors) = FermiSurfaceSector::extract(&graph, &lifted)?;
                assert_eq!(sectors.len(), 2);
                for sector in sectors {
                    assert!(FermiSurfaceSector::calls(&sector.coefficient)?.is_empty());
                    assert!(sector.coefficient.contains(&numerator));
                }
                if let Some(reference) = &reference {
                    assert_eq!(
                        &lifted, reference,
                        "contour canonicalization must be idempotent"
                    );
                } else {
                    reference = Some(lifted);
                }
            }
        }
        let warm = &prefactor * ThermalDistributionCall::atom(1, 0, 1, -1, 1);
        assert!(
            FermiSurfaceSector::lu_boundary(
                &graph,
                &warm,
                2,
                HybridSurfaceRef::Esurface(&surface),
                None,
            )?
            .is_zero()
        );
        let existing_delta = &prefactor * ThermalDistributionCall::atom(1, 1, 0, -1, 1);
        assert!(
            FermiSurfaceSector::lu_boundary(
                &graph,
                &existing_delta,
                2,
                HybridSurfaceRef::Esurface(&surface),
                None,
            )
            .is_err()
        );

        // Resolve physical signs only after selecting the production map.
        // An unsupported derivative on an inactive map must disappear first.
        let positive: EdgeVec<_> =
            vec![Orientation::Default; graph.loop_momentum_basis.edge_signatures.iter().count()]
                .into();
        let mut negative = positive.clone();
        negative[EdgeIndex(1)] = Orientation::Reversed;
        let catalogue = [positive, negative];
        let ids = [OrientationID(3), OrientationID(9)];
        let unresolved = ThermalDistributionCall::atom(1, 0, 0, -1, GS.sign(EdgeIndex(1)));
        let coefficient = &prefactor * (ids[0].atom() + Atom::num(2) * ids[1].atom()) * unresolved
            + OrientationID(15).atom() * existing_delta;
        let lifted = FermiSurfaceSector::lu_boundary(
            &graph,
            &coefficient,
            2,
            HybridSurfaceRef::Esurface(&surface),
            Some((&catalogue, &ids)),
        )?;
        for (id, orientation, multiple) in [(ids[0], 1, 1), (ids[1], -1, 2)] {
            let selected = id.select(&lifted);
            let expected = Atom::num(multiple)
                * &prefactor
                * ThermalDistributionCall::atom(1, 1, 0, -1, orientation)
                * &d1
                / (&d1 + &d2);
            assert_eq!(selected, expected);
        }

        // Independently differentiate the shared explicit on-shell expression.
        // Only this energy polynomial is expanded; the graph numerator stays opaque.
        let dilation = symbol!("cold_lu_dilation");
        let mut scaled_energy = graph.explicit_ose_atom(EdgeIndex(2));
        let mut actual = d2
            .replace(GS.ose(EdgeIndex(2)))
            .with(graph.explicit_ose_atom(EdgeIndex(2)));
        for component in 1..=3 {
            let index = GS.cind(component);
            let internal = graph.loop_momentum_basis.loop_atom(
                EdgeIndex(2),
                GS.emr_mom,
                std::slice::from_ref(&index),
                true,
            );
            let external = graph.loop_momentum_basis.ext_atom(
                EdgeIndex(2),
                GS.emr_mom,
                std::slice::from_ref(&index),
                true,
            );
            let edge_momentum = GS.emr_mom(EdgeIndex(2), index);
            scaled_energy = scaled_energy
                .replace(edge_momentum.to_pattern())
                .with(Atom::var(dilation) * &internal + &external);
            actual = actual.replace(edge_momentum).with(internal + external);
        }
        let expected = scaled_energy
            .derivative(dilation)
            .replace(dilation)
            .with(Atom::one());
        assert!((actual - expected).expand().is_zero());
        Ok(())
    }
}
