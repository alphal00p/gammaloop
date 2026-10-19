//! A two-energy cut localized in the polar angle of one Fermi sphere.

use color_eyre::eyre::{Result, ensure, eyre};
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use rug::Float;
use symbolica::domains::float::SingleFloat;

use crate::{
    cff::{
        esurface::ExternalShift,
        radial::{DirectedBound, ENCLOSURE_PRECISION, EnergySurfaceRay},
        surface::HybridSurfaceRef,
    },
    graph::Graph,
    integrands::process::sampling::maps::SamplingEvaluationError,
    momentum::{
        SignOrZero, ThreeMomentum,
        sample::{ExternalFourMomenta, ExternalThreeMomenta, LoopIndex, LoopMomenta},
    },
    utils::{F, FloatLike, compute_t_part_of_shift_part, hyperdual_utils::new_constant},
    uv::uv_graph::UVE,
};

use super::{FermiCutLocalization, FermiSurfaceProduct, PreparedFermiSurface};

impl FermiSurfaceProduct {
    pub(super) fn cut_equation(
        surface: HybridSurfaceRef<'_>,
    ) -> Result<(Vec<(EdgeIndex, i64)>, &ExternalShift)> {
        match surface {
            HybridSurfaceRef::Esurface(surface) => Ok((
                surface.energies.iter().map(|edge| (*edge, 1)).collect(),
                &surface.external_shift,
            )),
            HybridSurfaceRef::Hsurface(surface) => Ok((
                surface
                    .positive_energies
                    .iter()
                    .map(|edge| (*edge, 1))
                    .chain(surface.negative_energies.iter().map(|edge| (*edge, -1)))
                    .collect(),
                &surface.external_shift,
            )),
            _ => Err(eyre!("Fermi cut localization requires an E- or H-surface")),
        }
    }

    fn angular_cut_edges(
        &self,
        graph: &Graph,
        energies: &[(EdgeIndex, i64)],
    ) -> Result<[(EdgeIndex, i64); 2]> {
        ensure!(
            self.factors.len() == 1 && self.factors[0].derivative_order == 1,
            "Angular cut/Fermi localization supports one delta surface; higher distribution derivatives require observable-boundary intersections"
        );
        ensure!(
            energies.len() == 2,
            "Angular cut/Fermi localization requires exactly two signed energies"
        );
        let shell_mass = graph[self.factors[0].edge_id].mass_atom();
        let first = energies.iter().position(|(edge, _)| {
            let route = &self.lmb.edge_signatures[*edge];
            route.internal[LoopIndex(0)].is_sign()
                && route.internal.iter().skip(1).all(|sign| !sign.is_sign())
                && route.external.iter().all(|sign| !sign.is_sign())
                && graph[*edge].mass_atom() == shell_mass
        });
        let first = first.ok_or_else(|| {
            eyre!("The cut has no energy matching the Fermi shell's exact routing and mass")
        })?;
        let other = 1 - first;
        ensure!(
            self.lmb.edge_signatures[energies[other].0].internal[LoopIndex(0)].is_sign(),
            "The second cut energy must depend on the active Fermi momentum"
        );
        Ok([energies[first], energies[other]])
    }

    pub(super) fn classify_cut(
        &self,
        graph: &Graph,
        surface: HybridSurfaceRef<'_>,
    ) -> Result<FermiCutLocalization> {
        let (energies, _) = Self::cut_equation(surface)?;
        if energies.iter().all(|(edge, _)| {
            self.lmb.edge_signatures[*edge]
                .internal
                .iter()
                .take(self.factors.len())
                .all(|sign| *sign == SignOrZero::Zero)
        }) {
            return Ok(FermiCutLocalization::Independent);
        }
        self.angular_cut_edges(graph, &energies)?;
        Ok(FermiCutLocalization::CutAngle)
    }

    /// Preserve the sampled azimuth and integrate the cut delta in cos(theta).
    /// The unused original cosine has normalized density 1/2 on [-1,1].
    /// The shell's radial profile is already present in `prepared`; no second
    /// LU profile, radial Jacobian or compensating small-slope division occurs.
    pub(super) fn localize_cut_angle<T: FloatLike>(
        &self,
        prepared: &mut PreparedFermiSurface<T>,
        graph: &Graph,
        surface: HybridSurfaceRef<'_>,
        external: &ExternalFourMomenta<F<T>>,
        masses: &EdgeVec<F<T>>,
        residual_scale: &F<T>,
    ) -> Result<bool> {
        let (energies, shift) = Self::cut_equation(surface)?;
        let [(_, shell_sign), (other_edge, other_sign)] =
            self.angular_cut_edges(graph, &energies)?;
        let nu = &prepared.shell_energies[0];
        let zero = nu.zero();
        let one = nu.one();
        // Keep the same native shell mass used by the radial measure, including
        // masses defined by composite model expressions.
        let mass = prepared.shell_masses[0].clone();
        let other_mass = masses[other_edge].abs();
        let radius = (nu - &mass).sqrt() * (nu + &mass).sqrt();
        let mut inactive = LoopMomenta::from_iter(
            prepared
                .momenta
                .iter()
                .map(|momentum| momentum.map_ref(&|value| value.values[0].clone())),
        );
        let original = inactive[LoopIndex(0)].clone();
        inactive[LoopIndex(0)] = ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone());
        let spatial: ExternalThreeMomenta<F<T>> =
            external.iter().map(|p| p.spatial.clone()).collect();
        let route = &self.lmb.edge_signatures[other_edge];
        // Route a directly, excluding the active coordinate. Subtracting
        // epsilon*q from a completed edge momentum can lose a small shift.
        let offset = route.compute_momentum(&inactive, &spatial);
        let epsilon = route.internal[LoopIndex(0)] * 1_i64;
        let shift_value = shift.iter().fold(zero.clone(), |sum, (edge, sign)| {
            sum + nu.from_i64(*sign)
                * compute_t_part_of_shift_part(&self.lmb.edge_signatures[*edge].external, external)
        });
        let other_energy = (-&shift_value - nu.from_i64(shell_sign) * nu) * nu.from_i64(other_sign);

        if [
            &mass,
            &other_mass,
            &radius,
            &offset.px,
            &offset.py,
            &offset.pz,
            &shift_value,
            &other_energy,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
        {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "angular Fermi cut",
                detail: "nonfinite routed energy or spatial offset".into(),
            }
            .into());
        }
        let native = Float::enclose_native;
        let square = Float::enclosed_square;
        let multiply = Float::enclosed_mul;
        let radial_squared = Float::enclosed_sum(
            [(1, square(&native(nu))), (-1, square(&native(&mass)))].into_iter(),
        );
        let offset_squared = Float::enclosed_sum(
            [&offset.px, &offset.py, &offset.pz]
                .into_iter()
                .map(|value| (1, square(&native(value)))),
        );
        let energy = Float::enclosed_sum(
            [
                (-other_sign, native(&shift_value)),
                (-shell_sign * other_sign, native(nu)),
            ]
            .into_iter(),
        );
        if energy[1] < native(&other_mass)[0] {
            return Ok(false);
        }
        let uncertain = |detail: &str| SamplingEvaluationError::UncertainGeometry {
            detail: format!("angular Fermi cut: {detail}"),
        };
        if energy[0] <= 0 || radial_squared[0] <= 0 {
            return Err(uncertain("unresolved positive energy or regular shell radius").into());
        }
        let numerator = Float::enclosed_sum(
            [
                (1, square(&energy)),
                (-1, square(&native(&other_mass))),
                (-1, radial_squared.clone()),
                (-1, offset_squared.clone()),
            ]
            .into_iter(),
        );
        let denominator_squared = multiply(
            &native(&nu.from_i64(4)),
            &multiply(&radial_squared, &offset_squared),
        );
        let domain =
            Float::enclosed_sum([(1, square(&numerator)), (-1, denominator_squared)].into_iter());
        if domain[0] > 0 {
            return Ok(false);
        }
        if domain[1] >= 0 || offset_squared[0] <= 0 {
            return Err(uncertain("tangent or unresolved angular intersection").into());
        }
        let offset_norm = offset.norm_squared().sqrt();
        let cosine =
            (other_energy.square() - other_mass.square() - radius.square() - offset.norm_squared())
                / (nu.from_i64(2 * epsilon) * &offset_norm * &radius);
        if !cosine.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "angular Fermi cut",
                detail: "nonfinite angular projection".into(),
            }
            .into());
        }
        if cosine.abs() >= one {
            return Err(uncertain(
                "certified interior angle is not representable at this precision",
            )
            .into());
        }
        let axis = offset.map_ref(&|value| value / &offset_norm);
        // Subtracting the parallel projection leaves an O(epsilon*|q|)
        // longitudinal error, amplified when a translated sampling point is
        // nearly axial. The second cross product keeps the azimuthal vector
        // perpendicular even when the first cross product is small.
        let transverse = axis.cross_product(&original).cross_product(&axis);
        let transverse_norm = transverse.norm_squared().sqrt();
        if !transverse_norm.is_finite() || transverse_norm <= zero {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "angular Fermi cut",
                detail: "sampling point lies on the azimuthal coordinate pole".into(),
            }
            .into());
        }
        let sine = (one.clone() - cosine.square()).sqrt();
        let physical = axis.map_ref(&|value| value * &cosine * &radius)
            + transverse.map_ref(&|value| value / &transverse_norm * &sine * &radius);
        let template = &prepared.momenta.0[0].px;
        prepared.momenta[LoopIndex(0)] = physical.map_ref(&|value| new_constant(template, value));

        let mapped = LoopMomenta::from_iter(
            prepared
                .momenta
                .iter()
                .map(|momentum| momentum.map_ref(&|value| value.values[0].clone())),
        );
        let zeros = LoopMomenta::from_iter(
            mapped
                .iter()
                .map(|_| ThreeMomentum::new(zero.clone(), zero.clone(), zero.clone())),
        );
        let source = EnergySurfaceRay::routed(
            energies.iter().copied(),
            shift,
            &zeros,
            &mapped,
            external,
            masses,
            &self.lmb,
        );
        let unit = Float::with_val(ENCLOSURE_PRECISION, 1);
        let unit_range = [unit.clone(), unit];
        let (residual, _) = source.evaluate_interval(&unit_range, false)?;
        ensure!(
            residual_scale.is_finite() && residual_scale > &zero,
            "Angular Fermi cut requires a finite positive residual scale"
        );
        let tolerance = native(&(nu.epsilon() * nu.from_i64(64) * residual_scale));
        if residual[0].clone().abs().max(&residual[1].clone().abs()) > tolerance[0] {
            return Err(uncertain("mapped point fails the original cut equation").into());
        }
        let shell_residual = (physical.norm_squared() + mass.square()).sqrt() - nu;
        if shell_residual.abs() > nu.epsilon() * nu.from_i64(64) * nu.abs() {
            return Err(uncertain("mapped point fails its Fermi-shell equation").into());
        }
        // Laurent differentiation of a cold occupation supplies E'_F/eta'
        // in the ORIGINAL graph dilation. Its symbolic coefficient uses that
        // source frame, so certify the same derivative here. This explicitly
        // excludes radial tangencies even when the angular chart is regular.
        let parent = mapped.lmb_transform(&self.lmb, &graph.loop_momentum_basis, &spatial);
        let dilation = EnergySurfaceRay::routed(
            energies.into_iter(),
            shift,
            &parent,
            &zeros,
            external,
            masses,
            &graph.loop_momentum_basis,
        );
        let (_, slope) = dilation.evaluate_interval(&unit_range, true)?;
        if slope[0] <= 0 && slope[1] >= 0 {
            return Err(uncertain(
                "original LU derivative is zero or unresolved on the joint support",
            )
            .into());
        }
        let template = &prepared.momenta.0[0].px;
        prepared.factors.extend([
            (new_constant(template, &other_energy), false),
            (new_constant(template, &(one / nu.from_i64(2))), false),
            (new_constant(template, &offset_norm), true),
            (new_constant(template, &radius), true),
        ]);
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use spenso::algebra::complex::Complex;
    use symbolica::{domains::dual::HyperDual, prelude::Real};
    use three_dimensional_reps::ThermalDistributionFactor;
    use typed_index_collections::ti_vec;

    use super::*;
    use crate::{
        cff::{esurface::Esurface, hsurface::Hsurface},
        dot,
        graph::parse::IntoGraph,
        initialisation::test_initialise,
        momentum::FourMomentum,
        settings::runtime::{HFunction, HFunctionSettings},
        utils::{h, h_dual},
    };

    #[test]
    fn angular_fermi_cut_reproduces_two_body_shell_phase_space() -> Result<()> {
        test_initialise()?;
        let graph: Graph = dot!(digraph angular_fermi_bubble {
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
            vertex_set: crate::cff::VertexSet::dummy(),
        };
        let product = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: 1,
            }],
        )?;
        assert!(
            product.classify_cut(&graph, HybridSurfaceRef::Esurface(&surface))?
                == FermiCutLocalization::CutAngle
        );
        let masses: EdgeVec<_> = vec![F(5.0), F(1.5), F(1.5)].into();
        let external = ti_vec![FourMomentum::from_args(F(6.25), F(0.0), F(0.0), F(3.75))];
        let profile = HFunctionSettings {
            function: HFunction::Exponential,
            ..Default::default()
        };
        let mut point = product
            .prepare(
                &LoopMomenta(vec![ThreeMomentum::new(F(0.6), F(0.8), F(0.0))]),
                &[F(1.5)],
                &[F(2.5)],
                |scale| Ok(h_dual(scale, None, None, &profile)),
            )?
            .unwrap();
        assert!(product.localize_cut_angle(
            &mut point,
            &graph,
            HybridSurfaceRef::Esurface(&surface),
            &external,
            &masses,
            &F(5.0),
        )?);
        let q = point.momenta[LoopIndex(0)].map_ref(&|value| value.values[0]);
        assert!((q.norm_squared().0 - 4.0).abs() < 1e-14);
        assert!((q.pz.0 / 2.0 - 5.0 / 12.0).abs() < 1e-14);
        assert!((q.py.0 / q.px.0 - 4.0 / 3.0).abs() < 1e-14);
        let value = product.localize_prepared(
            point,
            HyperDual::from_values(product.shape.clone(), vec![Complex::new_re(F(1.0))]),
        )?;
        assert!((value.re.0 - 2.5 * h(&F(2.0), None, None, &profile).0).abs() < 1e-14);

        // Direct coarea identity: integral dPhi_2 delta(E_b-mu)=1/(8*pi*|P|).
        // This includes the unused angular coordinate's normalized density.
        let mut integral = 0.0;
        let step = 8.0 / 512.0;
        for index in 0..512 {
            let scale = (index as f64 + 0.5) * step;
            let radius = 2.0 / scale;
            let mut point = product
                .prepare(
                    &LoopMomenta(vec![ThreeMomentum::new(F(radius), F(0.0), F(0.0))]),
                    &[F(1.5)],
                    &[F(2.5)],
                    |scale| Ok(h_dual(scale, None, None, &profile)),
                )?
                .unwrap();
            assert!(product.localize_cut_angle(
                &mut point,
                &graph,
                HybridSurfaceRef::Esurface(&surface),
                &external,
                &masses,
                &F(5.0),
            )?);
            let coefficient = 1.0 / (16.0 * std::f64::consts::PI.powi(2) * 2.5 * 3.75);
            let value = product.localize_prepared(
                point,
                HyperDual::from_values(
                    product.shape.clone(),
                    vec![Complex::new_re(F(coefficient))],
                ),
            )?;
            integral += value.re.0 * 4.0 * std::f64::consts::PI * 8.0 / scale.powi(4) * step;
        }
        let expected = 1.0 / (8.0 * std::f64::consts::PI * 3.75);
        assert!((integral / expected - 1.0).abs() < 1e-12);

        for (mu, empty) in [
            (1.6, true),
            (4.7, true),
            (6.25, true),
            (1.625, false),
            (4.625, false),
        ] {
            let mut point = product
                .prepare(
                    &LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0))]),
                    &[F(1.5)],
                    &[F(mu)],
                    |scale| Ok(h_dual(scale, None, None, &profile)),
                )?
                .unwrap();
            let result = product.localize_cut_angle(
                &mut point,
                &graph,
                HybridSurfaceRef::Esurface(&surface),
                &external,
                &masses,
                &F(5.0),
            );
            if empty {
                assert!(!result?, "mu={mu}");
            } else {
                assert!(
                    result.is_err(),
                    "tangent mu={mu} must not become a physical zero"
                );
            }
        }
        let spacelike = ti_vec![FourMomentum::from_args(F(0.5), F(0.0), F(0.0), F(2.0))];
        for reversed in [false, true] {
            let signed = Hsurface {
                positive_energies: vec![EdgeIndex(if reversed { 2 } else { 1 })],
                negative_energies: vec![EdgeIndex(if reversed { 1 } else { 2 })],
                external_shift: vec![(EdgeIndex(0), if reversed { 1 } else { -1 })],
                vertex_set: crate::cff::VertexSet::dummy(),
            };
            let mut point = product
                .prepare(
                    &LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0))]),
                    &[F(1.5)],
                    &[F(2.5)],
                    |scale| Ok(h_dual(scale, None, None, &profile)),
                )?
                .unwrap();
            assert!(product.localize_cut_angle(
                &mut point,
                &graph,
                HybridSurfaceRef::Hsurface(&signed),
                &spacelike,
                &masses,
                &F(5.0),
            )?);
            assert!((point.momenta[LoopIndex(0)].pz.values[0].0 / 2.0 - 25.0 / 32.0).abs() < 1e-14);
        }
        let raised = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: 2,
            }],
        )?;
        assert!(
            raised
                .classify_cut(&graph, HybridSurfaceRef::Esurface(&surface))
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn angular_fermi_cut_preserves_rotated_shells_near_the_azimuthal_axis() -> Result<()> {
        fn check<T: FloatLike>(graph: &Graph, surface: &Esurface) -> Result<()> {
            let one = F::<T>::default().one();
            let zero = one.zero();
            let axis_norm = one.from_i64(14).sqrt();
            let axis = ThreeMomentum::new(
                &one / &axis_norm,
                one.from_i64(2) / &axis_norm,
                one.from_i64(3) / &axis_norm,
            );
            let probe = axis.map_ref(&|value| value * one.from_i64(15) / one.from_i64(4));
            let external = ti_vec![FourMomentum::from_args(
                one.from_i64(25) / one.from_i64(4),
                probe.px.clone(),
                probe.py.clone(),
                probe.pz.clone(),
            )];
            let mass = one.from_i64(3) / one.from_i64(2);
            let mu = one.from_i64(5) / one.from_i64(2);
            let masses: EdgeVec<_> = vec![one.from_i64(5), mass.clone(), mass.clone()].into();
            // Both physical daughters can be the occupied fermion. The
            // shifted daughter P-k becomes nearly axial at small input k.
            for shell_edge in [EdgeIndex(1), EdgeIndex(2)] {
                let product = FermiSurfaceProduct::new(
                    graph,
                    &[ThermalDistributionFactor {
                        edge_id: shell_edge,
                        sign: 1,
                        derivative_order: 1,
                    }],
                )?;
                for displacement in [0_i64, 1, 100, 10000] {
                    let input = if displacement == 0 {
                        // First fine-grid radial point and angular node from
                        // the translated negative-density Yukawa benchmark.
                        let radius = one.from_i64(2) / one.from_i64(255);
                        let boundary =
                            one.from_i64(35) / (one.from_i64(9) * one.from_i64(21).sqrt());
                        let cosine = -&one
                            + (one.clone() - boundary)
                                * (one.clone() + F::from_f64(-0.9894009349916499))
                                / one.from_i64(2);
                        let sine = (one.clone() - cosine.square()).sqrt();
                        let transverse_norm = one.from_i64(5).sqrt();
                        let transverse = ThreeMomentum::new(
                            one.from_i64(2) / &transverse_norm,
                            -&one / transverse_norm,
                            zero.clone(),
                        );
                        axis.map_ref(&|value| value * &radius * &cosine)
                            + transverse.map_ref(&|value| value * &radius * &sine)
                    } else {
                        let small = one.from_i64(displacement) / one.from_i64(1_000_000_000);
                        ThreeMomentum::new(small, zero.clone(), zero.clone())
                    };
                    let parent = LoopMomenta(vec![input]);
                    let spatial = external.iter().map(|p| p.spatial.clone()).collect();
                    let adapted =
                        parent.lmb_transform(&graph.loop_momentum_basis, product.lmb(), &spatial);
                    let mut point = product
                        .prepare(
                            &adapted,
                            std::slice::from_ref(&mass),
                            std::slice::from_ref(&mu),
                            |scale| Ok((-scale.clone()).exp()),
                        )?
                        .unwrap();
                    assert!(product.localize_cut_angle(
                        &mut point,
                        graph,
                        HybridSurfaceRef::Esurface(surface),
                        &external,
                        &masses,
                        &(one.from_i64(25) / one.from_i64(4)),
                    )?);
                    let occupied =
                        point.momenta[LoopIndex(0)].map_ref(&|value| value.values[0].clone());
                    let energy = (occupied.norm_squared() + mass.square()).sqrt();
                    assert!((energy - &mu).abs() <= one.epsilon() * one.from_i64(64) * &mu);
                    let aligned =
                        &occupied.px * &axis.px + &occupied.py * &axis.py + &occupied.pz * &axis.pz;
                    let expected = one.from_i64(if shell_edge == EdgeIndex(1) { 5 } else { -5 })
                        / one.from_i64(6);
                    assert!((aligned - expected).abs() <= one.epsilon() * one.from_i64(64));
                }
            }
            Ok(())
        }
        test_initialise()?;
        let graph: Graph = dot!(digraph rotated_fermi_bubble {
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
            vertex_set: crate::cff::VertexSet::dummy(),
        };
        check::<f64>(&graph, &surface)?;
        check::<crate::utils::f128>(&graph, &surface)?;
        check::<crate::utils::ArbPrec>(&graph, &surface)?;
        Ok(())
    }
}
