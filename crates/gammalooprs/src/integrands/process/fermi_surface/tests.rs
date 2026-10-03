use super::*;
use crate::{
    integrands::process::param_builder::ParamBuilderGraph,
    momentum::ThreeMomentum,
    settings::runtime::{HFunction, HFunctionSettings},
    utils::{
        ArbPrec, h, h_dual,
        symbols::{GS, ThermalDistributionLimit},
    },
    uv::uv_graph::UVE,
};
use symbolica::atom::{Atom, AtomCore};

#[test]
fn fermi_surface_derivatives_include_profile_and_spatial_jacobian() -> Result<()> {
    let graph = routing::test_graph()?;
    let profile = HFunctionSettings {
        function: HFunction::Exponential,
        ..Default::default()
    };
    let momenta = LoopMomenta(vec![
        ThreeMomentum::new(F(2.0), F(0.0), F(0.0)),
        ThreeMomentum::new(F(7.0), F(-2.0), F(1.0)),
    ]);
    // Here p_F=4 and t*=2. Independently differentiating
    // h(sqrt(E^2-9)/2) E(E^2-9)/16 gives these ordinary derivatives.
    for (order, multiplier) in [(1, 5.0), (2, 67.0 / 8.0), (3, 10.0)] {
        for sign in [-1, 1] {
            let product = FermiSurfaceProduct::new(
                &graph,
                &[ThermalDistributionFactor {
                    edge_id: EdgeIndex(1),
                    sign,
                    derivative_order: order,
                }],
            )?;
            let actual = product.localize(
                &momenta,
                &[F(3.0)],
                &[F(5.0)],
                |_, t| Ok(h_dual(t, None, None, &profile)),
                |mapped| {
                    assert!((mapped.0[0].px.values[0].0 - 4.0).abs() < 2e-15);
                    assert_eq!(mapped.0[1].px.values[0], F(7.0));
                    assert!(
                        mapped.0[1].px.values[1..]
                            .iter()
                            .all(|value| *value == F(0.0))
                    );
                    let jet = &mapped.0[0].px;
                    let value = new_constant(jet, &jet.values[0].one());
                    Ok(HyperDual::from_values(
                        product.shape.clone(),
                        value.values.into_iter().map(Complex::new_re).collect(),
                    ))
                },
            )?;
            let expected = multiplier * h(&F(2.0), None, None, &profile).0;
            assert!(
                (actual.re.0 - expected).abs() < 2e-14,
                "order={order}, sign={sign}: {actual:?} != {expected}"
            );
            assert_eq!(actual.im, F(0.0));
        }
    }
    Ok(())
}

#[test]
fn fermi_surface_integrated_derivatives_are_profile_independent() -> Result<()> {
    let graph = routing::test_graph()?;
    let (mass, nu) = (3.0_f64, 5.0_f64);
    let radius = (nu * nu - mass * mass).sqrt();
    let sphere_oracles = [
        nu * radius,
        -(2.0 * nu * nu - mass * mass) / radius,
        nu * (2.0 * nu * nu - 3.0 * mass * mass) / radius.powi(3),
    ];
    for (index, sphere_oracle) in sphere_oracles.into_iter().enumerate() {
        let product = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: index + 1,
            }],
        )?;
        for function in [HFunction::Exponential, HFunction::PolyExponential] {
            for sigma in [0.5, 2.0] {
                let profile = HFunctionSettings {
                    function: function.clone(),
                    sigma,
                    ..Default::default()
                };
                // Integrate the unchanged d^3k measure using t=p_F/|k|.
                // Midpoints avoid the coordinate endpoints; Gaussian tails
                // beyond eight profile widths are negligible.
                let step = 8.0 * sigma / 1024.0;
                let mut integral = 0.0;
                for i in 0..1024 {
                    let t = (i as f64 + 0.5) * step;
                    let sampling_radius = radius / t;
                    let momenta = LoopMomenta(vec![
                        ThreeMomentum::new(F(sampling_radius), F(0.0), F(0.0)),
                        ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
                    ]);
                    let value = product.localize(
                        &momenta,
                        &[F(mass)],
                        &[F(nu)],
                        |_, scale| Ok(h_dual(scale, None, None, &profile)),
                        |mapped| {
                            let jet = &mapped.0[0].px;
                            let value = new_constant(jet, &jet.values[0].one());
                            Ok(HyperDual::from_values(
                                product.shape.clone(),
                                value.values.into_iter().map(Complex::new_re).collect(),
                            ))
                        },
                    )?;
                    integral += value.re.0 * radius.powi(3) / t.powi(4) * step;
                }
                assert!(
                    (integral - sphere_oracle).abs() < 2e-10,
                    "order={}, {function:?}, sigma={sigma}: {integral} != {sphere_oracle}",
                    index + 1
                );
            }
        }
    }
    Ok(())
}

#[test]
fn fermi_surface_matches_integrated_finite_temperature_derivatives() -> Result<()> {
    let graph = routing::test_graph()?;
    let nu = 1.0_f64;
    let profile = HFunctionSettings {
        function: HFunction::Exponential,
        ..Default::default()
    };
    for derivative_order in [2, 3] {
        let product = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order,
            }],
        )?;
        let expression = graph
            .explicit_thermal_distribution_atom(
                EdgeIndex(1),
                derivative_order,
                Atom::num(1),
                Atom::num(1),
                ThermalDistributionLimit::Default,
            )
            .unwrap();
        let chemical_potential = graph[EdgeIndex(1)].chemical_potential_atom().unwrap();
        let mut thermal = expression
            .evaluator(&[
                GS.ose(EdgeIndex(1)),
                Atom::var(GS.inverse_temperature),
                chemical_potential,
            ])
            .function_map(graph.param_builder.fn_map.clone())
            .build()?
            .map_coeff(&|coefficient| coefficient.re.to_f64());
        let expected_zero = if derivative_order == 2 {
            -2.0 * nu
        } else {
            2.0
        };
        let step = 8.0 / 1024.0;
        let mut zero_temperature = 0.0;
        for i in 0..1024 {
            let t = (i as f64 + 0.5) * step;
            let momenta = LoopMomenta(vec![
                ThreeMomentum::new(F(nu / t), F(0.0), F(0.0)),
                ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ]);
            let value = product.localize(
                &momenta,
                &[F(0.0)],
                &[F(nu)],
                |_, scale| Ok(h_dual(scale, None, None, &profile)),
                |mapped| {
                    let jet = &mapped.0[0].px;
                    let value = new_constant(jet, &jet.values[0].one());
                    Ok(HyperDual::from_values(
                        product.shape.clone(),
                        value.values.into_iter().map(Complex::new_re).collect(),
                    ))
                },
            )?;
            zero_temperature += value.re.0 * nu.powi(3) / t.powi(4) * step;
        }
        assert!((zero_temperature - expected_zero).abs() < 2e-12);
        let mut previous_error = f64::INFINITY;
        for temperature in [0.2_f64, 0.1, 0.05] {
            let lower = -nu / temperature;
            let step = (40.0 - lower) / 8192.0;
            let finite_temperature = (0..8192)
                .map(|i| {
                    let x = lower + (i as f64 + 0.5) * step;
                    let energy = nu + temperature * x;
                    energy
                        * energy
                        * thermal.evaluate_single(&[energy, temperature.recip(), nu])
                        * temperature
                        * step
                })
                .sum::<f64>();
            // These are exact integrals of E^2 W_T'' and E^2 W_T'''
            // on E>=0; the common angular factor 4*pi is suppressed.
            let expected_finite = if derivative_order == 2 {
                -2.0 * (nu + temperature * (-nu / temperature).exp().ln_1p())
            } else {
                2.0 / (1.0 + (-nu / temperature).exp())
            };
            assert!(
                (finite_temperature - expected_finite).abs() < 2e-9,
                "order={derivative_order}, T={temperature}: {finite_temperature} != {expected_finite}"
            );
            let error = (finite_temperature - zero_temperature).abs();
            assert!(
                error < previous_error,
                "order={derivative_order}: finite-T integrated actions must converge"
            );
            previous_error = error;
        }
        assert!(previous_error < 1e-8);
    }
    Ok(())
}

#[test]
fn fermi_surface_independent_mixed_derivatives_commute() -> Result<()> {
    let graph = routing::test_graph()?;
    let profile = HFunctionSettings {
        function: HFunction::Exponential,
        ..Default::default()
    };
    let momenta = LoopMomenta(vec![
        ThreeMomentum::new(F(1.0), F(0.0), F(0.0)),
        ThreeMomentum::new(F(1.0), F(0.0), F(0.0)),
    ]);
    let expected = -104.0 * h(&F(1.0), None, None, &profile).0 * h(&F(2.0), None, None, &profile).0;
    for reversed in [false, true] {
        let mut factors = [
            ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: 2,
            },
            ThermalDistributionFactor {
                edge_id: EdgeIndex(3),
                sign: -1,
                derivative_order: 3,
            },
        ];
        let mut chemical_potentials = [F(1.0), F(2.0)];
        if reversed {
            factors.reverse();
            chemical_potentials.reverse();
        }
        let product = FermiSurfaceProduct::new(&graph, &factors)?;
        let actual = product.localize(
            &momenta,
            &[F(0.0), F(0.0)],
            &chemical_potentials,
            |_, scale| Ok(h_dual(scale, None, None, &profile)),
            |mapped| {
                let (first, second) = if reversed {
                    (&mapped.0[1].px, &mapped.0[0].px)
                } else {
                    (&mapped.0[0].px, &mapped.0[1].px)
                };
                // A coupled complex coefficient tests cross derivatives and
                // both observable components, not just a product of constants.
                let value = first.clone() + second.clone() * new_constant(second, &F(3.0));
                Ok(HyperDual::from_values(
                    product.shape.clone(),
                    value
                        .values
                        .into_iter()
                        .map(|value| Complex::new(value, value * F(2.0)))
                        .collect(),
                ))
            },
        )?;
        assert!(
            (actual.re.0 - expected).abs() < 2e-13,
            "reversed={reversed}: {actual:?} != {expected}"
        );
        assert!((actual.im.0 - 2.0 * expected).abs() < 4e-13);
    }
    Ok(())
}

#[test]
fn fermi_surface_differentiates_a_conditional_cut_root() -> Result<()> {
    let graph = routing::test_graph()?;
    let product = FermiSurfaceProduct::new(
        &graph,
        &[ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 2,
        }],
    )?;
    let profile = HFunctionSettings {
        function: HFunction::Exponential,
        ..Default::default()
    };
    let actual = product.localize(
        &LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0)); 2]),
        &[F(0.0)],
        &[F(1.0)],
        |_, scale| Ok(h_dual(scale, None, None, &profile)),
        |mapped| {
            // Toy massless cut: r_cut + E_F - 4 = 0. Its localized
            // radial weight is r_cut^2 / |dg_cut/dr_cut| = (4-E_F)^2.
            let energy = &mapped.0[0].px;
            let cut_root = new_constant(energy, &F(4.0)) - energy;
            let value = cut_root.clone() * &cut_root;
            Ok(HyperDual::from_values(
                product.shape.clone(),
                value.values.into_iter().map(Complex::new_re).collect(),
            ))
        },
    )?;
    let expected = -3.0 * h(&F(1.0), None, None, &profile).0;
    assert!(
        (actual.re.0 - expected).abs() < 2e-14,
        "{actual:?} != {expected}"
    );
    Ok(())
}

#[test]
fn fermi_surface_preserves_arbitrary_precision() -> Result<()> {
    let graph = routing::test_graph()?;
    let product = FermiSurfaceProduct::new(
        &graph,
        &[ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 3,
        }],
    )?;
    let zero = F::<ArbPrec>::default();
    let one = zero.one();
    let two = zero.from_usize(2);
    let actual = product.localize(
        &LoopMomenta(vec![
            ThreeMomentum::new(
                two.clone(),
                zero.clone(),
                zero.clone()
            );
            2
        ]),
        &[zero.from_usize(3)],
        &[zero.from_usize(5)],
        // h(t)=exp(-t) is normalized without binary64 constants.
        |_, scale| Ok((-scale.clone()).exp()),
        |mapped| {
            let jet = &mapped.0[0].px;
            let value = new_constant(jet, &jet.values[0].one());
            Ok(HyperDual::from_values(
                product.shape.clone(),
                value.values.into_iter().map(Complex::new_re).collect(),
            ))
        },
    )?;
    // d_E^2 [exp(-sqrt(E^2-9)/2) E(E^2-9)/16] at E=5.
    let expected = -(one.from_usize(125) / one.from_usize(128)) * (-two).exp();
    assert!(
        (actual.re - &expected).abs() < one.epsilon() * one.from_usize(1000),
        "arbitrary-precision derivative differs from its rational/exponential oracle"
    );
    assert_eq!(actual.im, zero);
    Ok(())
}

#[test]
fn fermi_surface_distinguishes_empty_support_from_invalid_shells() -> Result<()> {
    let graph = routing::test_graph()?;
    let product = FermiSurfaceProduct::new(
        &graph,
        &[ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 2,
        }],
    )?;
    let momenta = LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0)); 2]);
    for chemical_potential in [-5.0, 2.0] {
        let actual = product.localize(
            &momenta,
            &[F(3.0)],
            &[F(chemical_potential)],
            |_, _| panic!("empty support must not evaluate the profile"),
            |_| panic!("empty support must not evaluate the coefficient"),
        )?;
        assert_eq!(actual, Complex::new_re(F(0.0)));
    }
    for (mass, nu, expected) in [
        (3.0, 3.0, "Degenerate Fermi surface"),
        (-1.0, 3.0, "finite nonnegative mass"),
        (0.0, f64::INFINITY, "finite oriented chemical potential"),
        (f64::NAN, 3.0, "finite nonnegative mass"),
    ] {
        let error = product
            .localize(
                &momenta,
                &[F(mass)],
                &[F(nu)],
                |_, _| panic!("invalid shell must not evaluate the profile"),
                |_| panic!("invalid shell must not evaluate the coefficient"),
            )
            .unwrap_err();
        assert!(error.to_string().contains(expected), "{error}");
    }
    for radial_momentum in [0.0, f64::NAN, f64::INFINITY] {
        let mut invalid = momenta.clone();
        invalid.0[0].px = F(radial_momentum);
        assert!(
            product
                .localize(
                    &invalid,
                    &[F(0.0)],
                    &[F(1.0)],
                    |_, _| panic!("invalid momentum must not evaluate the profile"),
                    |_| panic!("invalid momentum must not evaluate the coefficient"),
                )
                .is_err()
        );
    }
    assert!(
        product
            .localize(
                &momenta,
                &[],
                &[F(1.0)],
                |_, _| panic!("invalid dimensions must not evaluate the profile"),
                |_| panic!("invalid dimensions must not evaluate the coefficient"),
            )
            .is_err()
    );
    Ok(())
}

#[test]
fn fermi_surface_preserves_representable_fully_weighted_scales() -> Result<()> {
    let graph = routing::test_graph()?;
    let product = FermiSurfaceProduct::new(
        &graph,
        &[ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 1,
        }],
    )?;
    for (sampling_radius, nu, coefficient, expected) in
        [(1e80, 1.0, 1e100, 1e-220), (1.0, 1e-200, 1e300, 1e-300)]
    {
        let actual = product.localize(
            &LoopMomenta(vec![
                ThreeMomentum::new(F(sampling_radius), F(0.0), F(0.0)),
                ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
            ]),
            &[F(0.0)],
            &[F(nu)],
            |_, scale| Ok((-scale.clone()).exp()),
            |mapped| {
                let projected = mapped.0[0].px.values[0].0;
                assert!(projected > 0.0);
                assert!((projected / nu - 1.0).abs() < 2e-15);
                let value = new_constant(&mapped.0[0].px, &F(coefficient));
                Ok(HyperDual::from_values(
                    product.shape.clone(),
                    value.values.into_iter().map(Complex::new_re).collect(),
                ))
            },
        )?;
        // For a massless shell the complete value is
        // coefficient * exp(-nu/|k|) * nu^3/|k|^4. Its factors cannot
        // safely be rounded separately at either of these scales.
        assert!(
            (actual.re.0 / expected - 1.0).abs() < 2e-13,
            "k={sampling_radius}, nu={nu}, coefficient={coefficient}: {} != {expected}",
            actual.re.0
        );
        assert_eq!(actual.im, F(0.0));
    }
    Ok(())
}

#[test]
fn fermi_surface_balances_profile_and_geometry_for_each_complex_component() -> Result<()> {
    let graph = routing::test_graph()?;
    let product = FermiSurfaceProduct::new(
        &graph,
        &[ThermalDistributionFactor {
            edge_id: EdgeIndex(1),
            sign: 1,
            derivative_order: 1,
        }],
    )?;
    let momenta = LoopMomenta(vec![
        ThreeMomentum::new(F(1e-80), F(0.0), F(0.0)),
        ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
    ]);
    for (coefficient, expected) in [
        ((1e-200, 0.0), (1e-40, 0.0)),
        ((1e-200, 1e100), (1e-40, 1e260)),
        ((1e100, 1e-200), (1e260, 1e-40)),
    ] {
        let actual = product.localize(
            &momenta,
            &[F(0.0)],
            &[F(1.0)],
            |_, scale| {
                // This positive profile has integral one on (0, infinity).
                let one = new_constant(scale, &scale.values[0].one());
                let denominator = scale.clone() + &one;
                Ok(one / (denominator.clone() * &denominator))
            },
            |mapped| {
                let value = new_constant(&mapped.0[0].px, &F(1.0));
                Ok(HyperDual::from_values(
                    product.shape.clone(),
                    value
                        .values
                        .into_iter()
                        .map(|entry| {
                            Complex::new(entry * F(coefficient.0), entry * F(coefficient.1))
                        })
                        .collect(),
                ))
            },
        )?;
        // t=1e80: h(t) is approximately 1e-160 and t^2/|k|^2=1e320.
        // Multiplying the small coefficient by h first loses its finite final
        // contribution. A large other component must not determine that order.
        for (value, target) in [(actual.re.0, expected.0), (actual.im.0, expected.1)] {
            if target == 0.0 {
                assert_eq!(value, 0.0);
            } else {
                assert!(
                    (value / target - 1.0).abs() < 3e-13,
                    "coefficient={coefficient:?}: {actual:?} != {expected:?}"
                );
            }
        }
    }
    Ok(())
}

#[test]
fn fermi_surface_products_preserve_inactive_derivative_axes() -> Result<()> {
    let graph = routing::test_graph()?;
    let momenta = LoopMomenta(vec![ThreeMomentum::new(F(1.0), F(0.0), F(0.0)); 2]);
    for (second_order, multiplier) in [(1, 8.0), (2, -4.0)] {
        let product = FermiSurfaceProduct::new(
            &graph,
            &[
                ThermalDistributionFactor {
                    edge_id: EdgeIndex(1),
                    sign: 1,
                    derivative_order: 1,
                },
                ThermalDistributionFactor {
                    edge_id: EdgeIndex(3),
                    sign: 1,
                    derivative_order: second_order,
                },
            ],
        )?;
        let actual = product.localize(
            &momenta,
            &[F(0.0), F(0.0)],
            &[F(1.0), F(2.0)],
            |_, scale| Ok((-scale.clone()).exp()),
            |mapped| {
                let jet = &mapped.0[0].px;
                let value = new_constant(jet, &jet.values[0].one());
                Ok(HyperDual::from_values(
                    product.shape.clone(),
                    value.values.into_iter().map(Complex::new_re).collect(),
                ))
            },
        )?;
        // The two factors are E^3 exp(-E) at E=1 and E=2. The
        // first axis has no derivative even when the second one does.
        let expected = multiplier * (-3.0_f64).exp();
        assert!(
            (actual.re.0 - expected).abs() < 2e-15,
            "orders=[1,{second_order}]: {actual:?} != {expected}"
        );
        assert_eq!(actual.im, F(0.0));
    }
    Ok(())
}
