//! Additional high-precision radial controls of complete local forest sums.
//! Exact operator, hard/soft Taylor and scalar residue certificates live in
//! local_4d::on_shell_forest_tests; these sampled slopes do not replace them.

use super::on_shell_tests::DirectForestFixture;
use crate::utils::{ArbPrec, F};
use color_eyre::Result;
use symbolica::prelude::*;

#[derive(Clone, Copy, Debug)]
enum SoftStratum {
    PhysicalChildGluon,
    ChildOnShellGluon,
    OuterGluon,
    SimultaneousGluons,
    RainbowOuterHierarchy(i32),
    RainbowPhysicalChildHierarchy(i32),
    RainbowFixedChildHierarchy(i32),
    DisjointLeftGluon,
    DisjointRightGluon,
    DisjointGluons,
    ThreeLoopOuterOnShellGluon,
    ThreeLoopMiddleOnShellGluon,
    ThreeLoopInnerOnShellGluon,
    ThreeLoopFlat {
        label: &'static str,
        centers: [[i32; 3]; 3],
        normal_slots: [Option<usize>; 3],
    },
}

impl SoftStratum {
    fn momenta(self, radius: &Atom, ray: usize) -> Vec<[Atom; 3]> {
        let direction = if ray == 0 { [1, 0, 0] } else { [0, 1, 0] };
        let transverse = if ray == 0 { [0, 1, 0] } else { [1, 0, 0] };
        let soft = std::array::from_fn::<_, 3, _>(|index| Atom::num(direction[index]) * radius);
        let hard_parent = [Atom::num(5), Atom::Zero, Atom::Zero];
        let hard_child = [Atom::Zero, Atom::num(9), Atom::Zero];
        match self {
            Self::PhysicalChildGluon => vec![
                hard_parent.clone(),
                std::array::from_fn(|index| &hard_parent[index] + &soft[index]),
            ],
            Self::ChildOnShellGluon | Self::DisjointRightGluon => vec![hard_parent, soft],
            Self::OuterGluon | Self::DisjointLeftGluon => vec![soft, hard_child],
            Self::SimultaneousGluons => vec![
                soft.clone(),
                std::array::from_fn(|index| &soft[index] + Atom::num(transverse[index]) * radius),
            ],
            Self::RainbowOuterHierarchy(power)
            | Self::RainbowPhysicalChildHierarchy(power)
            | Self::RainbowFixedChildHierarchy(power) => {
                // Exact rational unit directions with two transverse angles.
                let first = [Atom::one(), Atom::Zero, Atom::Zero];
                let second = if ray == 0 {
                    [Atom::num((3, 5)), Atom::num((4, 5)), Atom::Zero]
                } else {
                    [Atom::num((-3, 5)), Atom::Zero, Atom::num((4, 5))]
                };
                let large = std::array::from_fn::<_, 3, _>(|axis| &first[axis] * radius);
                let small =
                    std::array::from_fn::<_, 3, _>(|axis| &second[axis] * radius.pow(power));
                match self {
                    Self::RainbowOuterHierarchy(_) => vec![
                        std::array::from_fn(|axis| &first[axis] * radius.pow(power)),
                        std::array::from_fn(|axis| &second[axis] * radius),
                    ],
                    Self::RainbowPhysicalChildHierarchy(_) => vec![
                        large.clone(),
                        std::array::from_fn(|axis| &large[axis] + &small[axis]),
                    ],
                    Self::RainbowFixedChildHierarchy(_) => vec![large, small],
                    _ => unreachable!(),
                }
            }
            // Each component owns its own subtraction point. The second
            // soft gluon is measured from zero, without a shift by the first
            // component's loop momentum.
            Self::DisjointGluons => vec![
                soft,
                std::array::from_fn(|index| Atom::num(transverse[index]) * radius),
            ],
            Self::ThreeLoopOuterOnShellGluon => {
                vec![soft, hard_child, [Atom::Zero, Atom::Zero, Atom::num(16)]]
            }
            Self::ThreeLoopMiddleOnShellGluon => {
                vec![hard_parent, soft, [Atom::Zero, Atom::Zero, Atom::num(16)]]
            }
            Self::ThreeLoopInnerOnShellGluon => vec![hard_parent, hard_child, soft],
            Self::ThreeLoopFlat {
                centers,
                normal_slots,
                ..
            } => {
                // Two independent transverse direction sets. Every independent
                // three-vector normal scales uniformly with the same radius.
                let directions = if ray == 0 {
                    [[2, 1, -1], [-1, 3, 2], [1, -2, 3]]
                } else {
                    [[1, 3, 2], [-2, 1, 4], [3, -1, 2]]
                };
                centers
                    .into_iter()
                    .zip(normal_slots)
                    .map(|(center, slot)| {
                        std::array::from_fn(|axis| {
                            let fixed = Atom::num(center[axis]);
                            match slot {
                                Some(slot) => fixed + Atom::num(directions[slot][axis]) * radius,
                                None => fixed,
                            }
                        })
                    })
                    .collect()
            }
        }
    }

    fn radial_measure_power(self) -> i32 {
        match self {
            // Both three-dimensional loop vectors scale together. The radial
            // measure of this six-dimensional stratum is rho^5 d rho.
            Self::SimultaneousGluons | Self::DisjointGluons => 5,
            Self::RainbowOuterHierarchy(power)
            | Self::RainbowPhysicalChildHierarchy(power)
            | Self::RainbowFixedChildHierarchy(power) => 3 * power + 2,
            Self::ThreeLoopFlat { normal_slots, .. } => {
                // Rank r means 3r independent real normal coordinates.
                3 * normal_slots.iter().filter(|slot| slot.is_some()).count() as i32 - 1
            }
            _ => 2,
        }
    }
}

fn check_soft_rays(
    fixture: DirectForestFixture,
    batches: &[([i32; 4], &[SoftStratum])],
) -> Result<()> {
    crate::initialisation::test_initialise()?;
    // Reuse the complete-family exact source certificates and numeric engine.
    // The matched production measurements select explicit direct for the
    // soft-gluon insertion and projected 4D for the other nested fixtures.
    // Small pointwise equivalence tests retain all three representations.
    let route = if matches!(fixture, DirectForestFixture::SoftGluonInQuark) {
        1
    } else {
        2
    };
    let mut evaluator = fixture.full_forest_evaluator(route, false, ["1", "0", "0", "0"])?;
    // Numerical external data do not alter the prepared native forest.
    // Reuse its retained families and compiled scalar output probe across
    // the coincident-center and generic-external batches.
    for (external, strata) in batches {
        let external = (*external).map(Atom::num);
        let unit = F(ArbPrec::default()).one();
        let zero = unit.zero();
        for &stratum in *strata {
            for ray in 0..2 {
                let mut weighted = Vec::new();
                for exponent in 3..=7 {
                    let denominator = 1usize << exponent;
                    let radius = Atom::one() / Atom::num(denominator);
                    let values = evaluator.evaluate(&external, &stratum.momenta(&radius, ray))?;
                    let mut norm = zero.clone();
                    for value in values.values() {
                        let component = value.clone().norm().re;
                        assert!(
                            component.is_finite(),
                            "{fixture:?}, {stratum:?}, ray {ray}: finite tensor value"
                        );
                        norm = norm.max(component);
                    }
                    assert!(
                        norm > zero,
                        "{fixture:?}, {stratum:?}, ray {ray}: nonzero full-forest control"
                    );
                    let numeric_radius = &unit / unit.from_usize(denominator);
                    weighted.push(norm * numeric_radius.powi(stratum.radial_measure_power()));
                }
                let flat = match stratum {
                    SoftStratum::ThreeLoopFlat { label, .. } => label,
                    _ => "existing radial stratum",
                };
                crate::debug_tags!(#uv, #local, #soft, #summary;
                    stage = "complete_local_forest_soft_ray",
                    fixture = ?fixture, stratum = ?stratum, flat, route, ray,
                    radial_measure_power = stratum.radial_measure_power(),
                    file.external = ?external, file.weighted = ?weighted,
                    "Complete local CFF radial history at rho=2^-e, e=3..7"
                );
                // On the last two halvings, an uncanceled additional integer pole
                // doubles the expected ratio. The disjoint product is rho^-4,
                // so its rho^5-weighted norm should halve; the other strata only
                // require bounded radial weights. The allowance for subleading
                // terms excludes an additional integer pole in either case.
                // This is a finite-ray control, not a general symbolic IR proof.
                let ratio_denominator = if matches!(stratum, SoftStratum::DisjointGluons) {
                    4
                } else {
                    2
                };
                let maximum_growth = unit.from_usize(3) / unit.from_usize(ratio_denominator);
                for tail in weighted[2..].windows(2) {
                    assert!(
                        tail[1] <= &maximum_growth * &tail[0],
                        "{fixture:?}, {stratum:?}, ray {ray}: radial weight grows toward the soft point: {weighted:?}",
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn complete_os_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInQuark,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::PhysicalChildGluon,
                SoftStratum::ChildOnShellGluon,
                SoftStratum::OuterGluon,
                SoftStratum::SimultaneousGluons,
            ],
        )],
    )
}

#[test]
fn complete_os_in_ir_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::OsQuarkInGluon,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::PhysicalChildGluon,
                SoftStratum::ChildOnShellGluon,
            ],
        )],
    )
}

#[test]
fn complete_ir_in_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    check_soft_rays(
        DirectForestFixture::SoftGluonInQuark,
        &[([5, 0, 0, 0], &[SoftStratum::OuterGluon])],
    )
}

#[test]
fn complete_disjoint_os_forest_has_bounded_radial_soft_controls() -> Result<()> {
    // At the fixed rest-frame OS point each one-loop residue is at most
    // rho^-2. The live bridge has fixed q=(5,0,0,0), m=12, hence no soft
    // pole; the product is at most rho^-4 when both loops scale. Its six-
    // dimensional radial measure is rho^5 d rho. These sampled rays test
    // that bound for the complete (bare + CT)_left * bridge * (bare + CT)_right
    // sum; exact operator and scalar-residue identities are separate tests.
    check_soft_rays(
        DirectForestFixture::DisjointOsQuarks,
        &[(
            [5, 0, 0, 0],
            &[
                SoftStratum::DisjointLeftGluon,
                SoftStratum::DisjointRightGluon,
                SoftStratum::DisjointGluons,
            ],
        )],
    )
}

#[test]
fn complete_three_loop_local_os_forest_has_bounded_uniform_soft_flats() -> Result<()> {
    // Retain the original q=(5,0,0,0) degeneracy control first: the physical
    // outer center coincides with the rest point. The other two carriers
    // stay generic while each of p,k,l separately approaches zero.
    let rest_centers = [
        SoftStratum::ThreeLoopInnerOnShellGluon,
        SoftStratum::ThreeLoopMiddleOnShellGluon,
        SoftStratum::ThreeLoopOuterOnShellGluon,
    ];
    // The factory certifies the physical gluon planes A:p=q, B:k=p, C:l=k
    // and the component-owned rest-point planes D:p=0, E:k=0, F:l=0.
    // A nonzero spatial q separates physical from artificial OS centers.
    // These are the 20 distinct consistent nonempty affine intersections:
    // rank 1: 6, rank 2: 10, rank 3: 4. DE also contains B; EF contains C.
    // This is a uniform finite-ray check of all eight LOCAL forest families.
    // It does not yet include finite integrated prefixes, hierarchical
    // rho/rho^2 approaches, or an all-angle/global integrability claim.
    let q = [1, 2, 3];
    let p = [5, 0, 0];
    let k = [0, 9, 0];
    let l = [0, 0, 16];
    let zero = [0, 0, 0];
    let flats = [
        ("A", [q, k, l], [Some(0), None, None]),
        ("B", [p, p, l], [None, Some(0), None]),
        ("C", [p, k, k], [None, None, Some(0)]),
        ("D", [zero, k, l], [Some(0), None, None]),
        ("E", [p, zero, l], [None, Some(0), None]),
        ("F", [p, k, zero], [None, None, Some(0)]),
        ("AB", [q, q, l], [Some(0), Some(1), None]),
        ("AC", [q, k, k], [Some(0), None, Some(1)]),
        ("BC", [p, p, p], [None, Some(0), Some(1)]),
        ("DE", [zero, zero, l], [Some(0), Some(1), None]),
        ("DF", [zero, k, zero], [Some(0), None, Some(1)]),
        ("EF", [p, zero, zero], [None, Some(0), Some(1)]),
        ("AE", [q, zero, l], [Some(0), Some(1), None]),
        ("AF", [q, k, zero], [Some(0), None, Some(1)]),
        ("BF", [p, p, zero], [None, Some(0), Some(1)]),
        ("CD", [zero, k, k], [Some(0), None, Some(1)]),
        ("ABC", [q, q, q], [Some(0), Some(1), Some(2)]),
        ("DEF", [zero, zero, zero], [Some(0), Some(1), Some(2)]),
        ("ABF", [q, q, zero], [Some(0), Some(1), Some(2)]),
        ("AEF", [q, zero, zero], [Some(0), Some(1), Some(2)]),
    ];
    let strata = flats.map(
        |(label, centers, normal_slots)| SoftStratum::ThreeLoopFlat {
            label,
            centers,
            normal_slots,
        },
    );
    check_soft_rays(
        DirectForestFixture::OsQuarkInQuarkWithSpectator,
        &[
            ([5, 0, 0, 0], rest_centers.as_slice()),
            ([5, 1, 2, 3], strata.as_slice()),
        ],
    )
}

#[test]
fn rainbow_hierarchy_coordinates_resolve_all_three_soft_normals_exactly() -> Result<()> {
    use crate::{
        utils::GS,
        uv::{
            UVgenerationSettings, UltravioletGraph,
            approx::{ForestNodeLike, local_4d::on_shell_forest_tests::nested_quark_fixture},
            uv_graph::UVE,
        },
    };
    use linnet::half_edge::{
        involution::{EdgeIndex, Flow},
        subgraph::SubSetLike,
    };

    crate::initialisation::test_initialise()?;
    // Geometry only: the three normals are r=e2-p, q=e3-e2 and u=e3-p.
    // In particular u<<r is q approximately -r, not q approximately zero.
    for power in [2, 3] {
        for stratum in [
            SoftStratum::RainbowOuterHierarchy(power),
            SoftStratum::RainbowPhysicalChildHierarchy(power),
            SoftStratum::RainbowFixedChildHierarchy(power),
        ] {
            assert_eq!(stratum.radial_measure_power(), 3 * power + 2);
            for ray in 0..2 {
                for exponent in 3..=7 {
                    let radius = Atom::one() / Atom::num(1usize << exponent);
                    let loops = stratum.momenta(&radius, ray);
                    let r = &loops[0];
                    let u = &loops[1];
                    let q: [Atom; 3] = std::array::from_fn(|axis| &u[axis] - &r[axis]);
                    let norm2 = |v: &[Atom; 3]| v.iter().map(|x| x.pow(2)).sum::<Atom>();
                    assert!(!norm2(r).is_zero());
                    assert!(!norm2(&q).is_zero());
                    assert!(!norm2(u).is_zero());
                    assert!(
                        (norm2(r)
                            - match stratum {
                                SoftStratum::RainbowOuterHierarchy(_) => radius.pow(2 * power),
                                _ => radius.pow(2),
                            })
                        .is_zero()
                    );
                    let small = match stratum {
                        SoftStratum::RainbowOuterHierarchy(_) => r,
                        SoftStratum::RainbowPhysicalChildHierarchy(_) => &q,
                        SoftStratum::RainbowFixedChildHierarchy(_) => u,
                        _ => unreachable!(),
                    };
                    assert!((norm2(small) - radius.pow(2 * power)).is_zero());
                    for axis in 0..3 {
                        assert!((&q[axis] + &r[axis] - &u[axis]).is_zero());
                    }
                }
            }
        }
    }
    let (graph, parent, child, _) = nested_quark_fixture(false)?;
    assert_eq!(
        graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        [EdgeIndex(2), EdgeIndex(3)]
    );
    // Certify the physical inner-gluon plane from the stored source routing.
    let physical = graph
        .loop_momentum_basis
        .loop_atom::<Atom>(EdgeIndex(4), GS.emr_mom, &[], true)
        + graph
            .loop_momentum_basis
            .ext_atom::<Atom>(EdgeIndex(4), GS.emr_mom, &[], true);
    let expected = function!(GS.emr_mom, 3) - function!(GS.emr_mom, 2);
    assert!(
        (&physical - &expected).expand().is_zero() || (&physical + &expected).expand().is_zero()
    );
    let settings = UVgenerationSettings {
        on_shell_subtraction_velocity: ["5/4", "1/4", "1/2", "1/2"].map(str::to_owned),
        ..Default::default()
    };
    let velocity = settings.on_shell_velocity()?;
    let future = [15, 3, 6, 6].map(Atom::num);
    let mass = graph[EdgeIndex(0)].mass_atom(&graph.model);
    assert!(!mass.is_zero());
    assert_eq!(mass, graph[EdgeIndex(2)].mass_atom(&graph.model));
    for (owner, node) in [(EdgeIndex(0), &parent), (EdgeIndex(2), &child)] {
        let crown = graph.dummy_less_full_crown(node.subgraph());
        let incoming = crown
            .included_iter()
            .filter(|hedge| graph.flow(*hedge) == Flow::Sink)
            .collect::<Vec<_>>();
        assert_eq!(incoming.len(), 1);
        assert_eq!(graph[&incoming[0]], owner);
        for sign in [1, -1] {
            // The two actual component owners retain the same fixed velocity;
            // only the branch's signed physical mass chooses its center.
            let point = GS
                .on_shell_momentum(owner, &(Atom::num(sign) * &mass), &velocity)
                .replace(mass.to_pattern())
                .with(12);
            for (axis, component) in future.iter().enumerate() {
                assert_eq!(
                    GS.resolve_os_momenta(&GS.indexed_momentum(&point, &[GS.cind(axis)])),
                    Atom::num(sign) * component
                );
            }
        }
    }
    let variable = symbol!("os_fixed_center_geometry::rho");
    let radius = Atom::var(variable);
    let mut charts = 0;
    for live_sign in [1, -1] {
        let external = future
            .clone()
            .map(|component| Atom::num(live_sign) * component);
        for relative_signs in [[-1, -1], [1, -1], [-1, 1]] {
            let center_signs = relative_signs.map(|sign| live_sign * sign);
            let centers = center_signs.map(|sign| {
                std::array::from_fn::<_, 3, _>(|axis| Atom::num(sign) * &future[axis + 1])
            });
            let offset: [Atom; 3] =
                std::array::from_fn(|axis| &centers[1][axis] - &centers[0][axis]);
            let mixed = center_signs[0] != center_signs[1];
            let norm2 = |v: &[Atom; 3]| v.iter().map(|x| x.pow(2)).sum::<Atom>();
            assert_eq!(norm2(&offset), Atom::num(if mixed { 324 } else { 0 }));
            for power in [2, 3] {
                for stratum in [
                    SoftStratum::RainbowOuterHierarchy(power),
                    SoftStratum::RainbowPhysicalChildHierarchy(power),
                    SoftStratum::RainbowFixedChildHierarchy(power),
                ] {
                    if mixed && matches!(stratum, SoftStratum::RainbowPhysicalChildHierarchy(_)) {
                        continue;
                    }
                    for ray in 0..2 {
                        let normals = stratum.momenta(&radius, ray);
                        let loops = std::array::from_fn::<_, 2, _>(|index| {
                            std::array::from_fn::<_, 3, _>(|axis| {
                                &centers[index][axis] + &normals[index][axis]
                            })
                        });
                        let r: [Atom; 3] =
                            std::array::from_fn(|axis| &loops[0][axis] - &centers[0][axis]);
                        let w: [Atom; 3] =
                            std::array::from_fn(|axis| &loops[1][axis] - &centers[1][axis]);
                        let q: [Atom; 3] =
                            std::array::from_fn(|axis| &loops[1][axis] - &loops[0][axis]);
                        assert_eq!(r, normals[0]);
                        assert_eq!(w, normals[1]);
                        for index in 0..2 {
                            for axis in 0..3 {
                                assert_eq!(
                                    loops[index][axis].replace(variable).with(0),
                                    centers[index][axis]
                                );
                            }
                        }
                        for axis in 0..3 {
                            // Normalize only this finite scalar rho polynomial;
                            // no graph numerator enters the geometry certificate.
                            assert!(
                                (&q[axis] - &offset[axis] - (&w[axis] - &r[axis]))
                                    .expand()
                                    .is_zero()
                            );
                            assert_eq!(q[axis].replace(variable).with(0), offset[axis]);
                            assert_eq!(
                                external[axis + 1],
                                Atom::num(live_sign) * &future[axis + 1]
                            );
                        }
                        let small = match stratum {
                            SoftStratum::RainbowOuterHierarchy(_) => r,
                            SoftStratum::RainbowPhysicalChildHierarchy(_) => q,
                            SoftStratum::RainbowFixedChildHierarchy(_) => w,
                            _ => unreachable!(),
                        };
                        assert!((norm2(&small) - radius.pow(2 * power)).expand().is_zero());
                        charts += 1;
                    }
                }
            }
        }
    }
    assert_eq!(charts, 56);
    Ok(())
}

mod slow {
    use super::*;

    #[test]
    fn complete_os_rainbow_has_bounded_hierarchical_soft_controls() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let fixture = DirectForestFixture::OsQuarkInQuark;
        // The integrated public rainbow agrees across all three routes. Its
        // projected route had the fastest evaluator in the retained 2026-10-20
        // correctness run; retain the smaller equivalence cases for all routes.
        const ROUTE: usize = 2;
        let unit = F(ArbPrec::default()).one();
        let zero = unit.zero();
        for (future, velocity) in [
            ([12, 0, 0, 0], ["1", "0", "0", "0"]),
            ([15, 3, 6, 6], ["5/4", "1/4", "1/2", "1/2"]),
        ] {
            // One actual finite-child forest and fixed component frame, reused for
            // both signed external endpoints and every hierarchy/angular point.
            let mut evaluator = fixture.full_forest_evaluator(ROUTE, true, velocity)?;
            for sign in [1, -1] {
                let external = future.map(|value| Atom::num(sign * value));
                for power in [2, 3] {
                    for stratum in [
                        SoftStratum::RainbowOuterHierarchy(power),
                        SoftStratum::RainbowPhysicalChildHierarchy(power),
                        SoftStratum::RainbowFixedChildHierarchy(power),
                    ] {
                        for ray in 0..2 {
                            let mut weighted = Vec::new();
                            for exponent in 3..=7 {
                                let denominator = 1usize << exponent;
                                let radius = Atom::one() / Atom::num(denominator);
                                let loops = stratum
                                    .momenta(&radius, ray)
                                    .into_iter()
                                    .map(|normal| {
                                        std::array::from_fn(|axis| {
                                            &external[axis + 1] + &normal[axis]
                                        })
                                    })
                                    .collect::<Vec<_>>();
                                let values = evaluator.evaluate(&external, &loops)?;
                                let norm =
                                    values.values().try_fold(zero.clone(), |norm, value| {
                                        let component = value.clone().norm().re;
                                        eyre::ensure!(
                                            component.is_finite(),
                                            "finite complete hierarchical rainbow tensor"
                                        );
                                        Ok::<_, eyre::Report>(norm.max(component))
                                    })?;
                                assert!(
                                    norm > zero,
                                    "nonzero complete hierarchical rainbow control"
                                );
                                let numeric_radius = &unit / unit.from_usize(denominator);
                                weighted.push(
                                    norm * numeric_radius.powi(stratum.radial_measure_power()),
                                );
                            }
                            crate::debug_tags!(#uv, #local, #soft, #summary;
                                stage = "complete_rainbow_hierarchical_soft_control",
                                route = ROUTE, stratum = ?stratum, power, ray, sign,
                                radial_measure_power = stratum.radial_measure_power(),
                                expected_margin = power - 1,
                                file.velocity = ?velocity, file.external = ?external,
                                file.weighted = ?weighted,
                                "Finite-grid causal-power-count control; not a uniform IR/epsilon proof"
                            );
                            // Independent chart power count: with r and u as the
                            // two vector coordinates, rho/rho^a has Jacobian
                            // rho^(3a+2). The candidate causal bounds r^-2 u^-3,
                            // q^-2 r^-3 and u^-2 r^-3 yield rho^(a-1) after that
                            // measure. The first uses the matched outer endpoint
                            // zero and one inner derivative. This grid TESTS that
                            // necessary power-count hypothesis on native CFF data;
                            // it does not prove its uniform angular remainder.
                            // No exponent is fitted to the produced values.
                            let maximum_growth =
                                unit.from_usize(3) / unit.from_usize(1usize << power);
                            for tail in weighted[2..].windows(2) {
                                assert!(
                                    tail[1] <= &maximum_growth * &tail[0],
                                    "{stratum:?}, ray {ray}, sign {sign}: hierarchical weight violates the derived candidate rho^{} margin: {weighted:?}",
                                    power - 1
                                );
                            }
                        }
                    }
                }
            }
            if future[1..].iter().all(|component| *component == 0) {
                continue;
            }
            // Keep the original 48 matched-center groups above unchanged. These
            // additional 56 groups fix the live external point while choosing
            // the two component-owned OS centers independently. At rest all
            // spatial centers coincide, so only the nonaxial frame is repeated.
            for sign in [1, -1] {
                let external = future.map(|component| Atom::num(sign * component));
                for (center_class, relative_signs) in [
                    ("opposite", [-1, -1]),
                    ("mixed_outer_matched", [1, -1]),
                    ("mixed_child_matched", [-1, 1]),
                ] {
                    let center_signs = relative_signs.map(|relative| sign * relative);
                    let centers = center_signs.map(|center_sign| {
                        std::array::from_fn::<_, 3, _>(|axis| {
                            Atom::num(center_sign * future[axis + 1])
                        })
                    });
                    let offset: [Atom; 3] =
                        std::array::from_fn(|axis| &centers[1][axis] - &centers[0][axis]);
                    let mixed = center_signs[0] != center_signs[1];
                    assert_eq!(
                        offset
                            .iter()
                            .map(|component| component.pow(2))
                            .sum::<Atom>(),
                        Atom::num(if mixed { 324 } else { 0 })
                    );
                    for power in [2, 3] {
                        for stratum in [
                            SoftStratum::RainbowOuterHierarchy(power),
                            SoftStratum::RainbowPhysicalChildHierarchy(power),
                            SoftStratum::RainbowFixedChildHierarchy(power),
                        ] {
                            if mixed
                                && matches!(stratum, SoftStratum::RainbowPhysicalChildHierarchy(_))
                            {
                                continue;
                            }
                            let expected_margin = if mixed { power } else { power - 1 };
                            for ray in 0..2 {
                                let mut weighted = Vec::new();
                                for exponent in 3..=7 {
                                    let denominator = 1usize << exponent;
                                    let radius = Atom::one() / Atom::num(denominator);
                                    let loops = stratum
                                        .momenta(&radius, ray)
                                        .into_iter()
                                        .enumerate()
                                        .map(|(index, normal)| {
                                            std::array::from_fn(|axis| {
                                                &centers[index][axis] + &normal[axis]
                                            })
                                        })
                                        .collect::<Vec<_>>();
                                    let values = evaluator.evaluate(&external, &loops)?;
                                    let norm =
                                        values.values().try_fold(zero.clone(), |norm, value| {
                                            let component = value.clone().norm().re;
                                            eyre::ensure!(
                                                component.is_finite(),
                                                "finite complete fixed-center rainbow tensor"
                                            );
                                            Ok::<_, eyre::Report>(norm.max(component))
                                        })?;
                                    assert!(
                                        norm > zero,
                                        "nonzero complete fixed-center rainbow control"
                                    );
                                    let numeric_radius = &unit / unit.from_usize(denominator);
                                    weighted.push(
                                        norm * numeric_radius.powi(stratum.radial_measure_power()),
                                    );
                                }
                                crate::debug_tags!(#uv, #local, #soft, #summary;
                                    stage="complete_rainbow_fixed_center_soft_control", route=ROUTE,
                                    center_class, stratum=?stratum, power, ray, sign,
                                    center_sigma=center_signs[0], center_tau=center_signs[1],
                                    relative_center_sigma=relative_signs[0], relative_center_tau=relative_signs[1],
                                    radial_measure_power=stratum.radial_measure_power(), expected_margin,
                                    file.velocity=?velocity, file.external=?external, file.centers=?centers,
                                    file.physical_child_center_offset=?offset, file.weighted=?weighted,
                                    source_assumptions=if mixed { "opposite ordered P_sigma P_tau endpoint; inner fixed-O soft pole and doubled-quark CFF leading coefficient; candidate r^-2 w^-2" } else { "complete child at the opposite matched endpoint; same-center three-normal candidate bounds" },
                                    "Independent fixed-center finite-grid control; not a uniform CFF/epsilon theorem"
                                );
                                // The mixed physical inner gluon remains hard:
                                // q=(tau-sigma)*p_vec+w-r. Ordered P_sigma P_tau=0
                                // removes the leading outer doubled-quark pole.
                                // The candidate r^-2 w^-2 leaves margin a after
                                // rho^(3a+2); opposite same-center cases retain
                                // the original a-1 hypothesis. Neither is a
                                // proof of a uniform angular CFF remainder.
                                let maximum_growth = unit.from_usize(3)
                                    / unit.from_usize(1usize << (expected_margin + 1));
                                for tail in weighted[2..].windows(2) {
                                    assert!(
                                        tail[1] <= &maximum_growth * &tail[0],
                                        "{center_class}, {stratum:?}, ray {ray}, live sign {sign}: fixed-center weight violates candidate rho^{expected_margin}: {weighted:?}"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
